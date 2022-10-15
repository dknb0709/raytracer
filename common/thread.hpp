#include <condition_variable>
#include <functional>
#include <future>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>

template <class F, class... Args>
using ResultType = typename std::result_of<F(Args...)>::type;

template <class F, class... Args>
auto makeTask(F&& f, Args&&... args)
    -> std::shared_ptr<std::packaged_task<ResultType<F, Args...>()>> {
  return std::make_shared<std::packaged_task<ResultType<F, Args...>()>>(
      std::bind(f, args...));
}

class ThreadPool;

class Worker {
 public:
  explicit Worker(ThreadPool& pool) : pool(pool) {}
  ~Worker() {}
  void operator()();

 private:
  ThreadPool& pool;
};

class ThreadPool {
 public:
  ThreadPool(std::size_t nThreads) : stop{false} {
    for (int i = 0; i < nThreads; ++i) {
      workers.emplace_back(Worker(*this));
    }
  }

  ~ThreadPool() {
    {
      std::lock_guard<std::mutex> lock(queue_mutex);
      stop = true;
    }
    condition.notify_all();
    for (auto& worker : workers) {
      worker.join();
    }
  }

  template <typename F, typename... Args>
  auto enqueue(F&& f, Args&&... args) -> std::future<ResultType<F, Args...>> {
    auto task = makeTask(f, args...);
    {
      std::unique_lock<std::mutex> lock(queue_mutex);
      tasks.emplace_back([task]() { (*task)(); });
    }
    condition.notify_one();
    return task->get_future();
  }

 private:
  friend class Worker;

  std::vector<std::thread> workers;
  std::deque<std::function<void()>> tasks;

  std::mutex queue_mutex;
  std::condition_variable condition;
  bool stop;
};

void Worker::operator()() {
  while (true) {
    std::function<void()> task;
    {
      std::unique_lock<std::mutex> lock(pool.queue_mutex);
      pool.condition.wait(
          lock, [this]() { return pool.stop || !pool.tasks.empty(); });
      if (pool.stop) return;
      task = std::move(pool.tasks.front());
      pool.tasks.pop_front();
    }
    task();
  }
}
