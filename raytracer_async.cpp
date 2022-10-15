#include <algorithm>
#include <cfloat>
#include <ctime>
#include <future>
#include <random>
#include <thread>

#include "common/Vec3.hpp"
#include "common/bitmap.hpp"

const int threadCount = 12;
const int samples = 120;
const int W = 600, H = 400;

inline float randf() {
  static bool initialized = false;
  if (!initialized) {
    srand(time(nullptr));
    initialized = true;
  }
  return static_cast<float>(rand()) / static_cast<float>(RAND_MAX + 1.0);
}

inline Vec3<float> randVec() { return {randf(), randf(), randf()}; }

inline Vec3<float> randUnitSphere() {
  Vec3<float> v;
  do {
    v = 2.0f * randVec() - Vec3<float>(1.0f);
  } while (v.length() >= 1.0f);
  return v;
}

struct Ray {
  Vec3<float> origin;
  Vec3<float> direction;
  Ray() = default;
  Ray(const Vec3<float>& o, const Vec3<float>& d) : origin{o}, direction{d} {}
  Vec3<float> at(float t) const { return origin + t * direction; }
};

struct Material;

struct HitRec {
  float t;
  Vec3<float> p;
  Vec3<float> n;
  Material* material;

  HitRec() : t{0.0f}, p{}, n{} {}
};

struct ScatterRec {
  Ray ray;
  Vec3<float> albedo;
};

struct Material {
  virtual bool scatter(const Ray& ray, const HitRec& hrec,
                       ScatterRec& srec) const = 0;
  virtual Vec3<float> emit(const Ray& ray, const HitRec& hrec) const {
    return Vec3<float>{0.0f};
  };
};

struct Lambertian : Material {
  Vec3<float> albedo;
  Lambertian(const Vec3<float>& albedo) : albedo{albedo} {}
  virtual bool scatter(const Ray& ray, const HitRec& hrec,
                       ScatterRec& srec) const {
    Vec3<float> target = hrec.p + hrec.n + randUnitSphere();
    Vec3<float> direction = target - hrec.p;
    srec.ray = Ray{hrec.p, direction};
    srec.albedo = albedo;
    return true;
  }
};

struct Metal : Material {
  Vec3<float> albedo;
  Metal(const Vec3<float>& albedo) : albedo{albedo} {}
  virtual bool scatter(const Ray& ray, const HitRec& hrec,
                       ScatterRec& srec) const {
    auto d = ray.direction.normalize();
    auto n = hrec.n;
    Vec3<float> reflected = d - 2.0f * d.dot(n) * n;
    srec.ray = Ray{hrec.p, reflected};
    srec.albedo = albedo;
    return reflected.dot(hrec.n) > 0;
  }
};

struct DiffuseLight : Material {
  Vec3<float> color;
  DiffuseLight(const Vec3<float>& color) : color{color} {}
  virtual bool scatter(const Ray& ray, const HitRec& hrec,
                       ScatterRec& srec) const {
    return false;
  }
  virtual Vec3<float> emit(const Ray& ray, const HitRec& hrec) const {
    return color;
  };
};

struct Dielectric : Material {
  float ri;
  Dielectric(float ri) : ri{ri} {}
  bool refract(const Vec3<float>& v, const Vec3<float>& n, float ratio,
               Vec3<float>& refracted) const {
    Vec3<float> uv = v.normalize();
    float dt = uv.dot(n);
    float D = 1.0f - ratio * ratio * (1.0f - dt * dt);
    if (D > 0.0f) {
      refracted = -ratio * (uv - dt * n) - std::sqrt(D) * n;
      return true;
    }
    return false;
  }

  float schlick(float cosine, float ri) const {
    float r0 = ((1.0f - ri) / (1.0f + ri)) * ((1.0f - ri) / (1.0f + ri));
    return r0 + (1.0f - r0) * (1.0f - cosine) * (1.0f - cosine) *
                    (1.0f - cosine) * (1.0f - cosine) * (1.0f - cosine);
  }

  virtual bool scatter(const Ray& ray, const HitRec& hrec,
                       ScatterRec& srec) const {
    Vec3<float> n;
    Vec3<float> reflected =
        ray.direction - 2.0f * ray.direction.dot(hrec.n) * hrec.n;
    float cosine;
    float ratio;
    if (ray.direction.dot(hrec.n) > 0.0f) {
      n = -hrec.n;
      ratio = ri;
      cosine = ri * ray.direction.dot(hrec.n) / ray.direction.length();
    } else {
      n = hrec.n;
      ratio = 1.0f / ri;
      cosine = -ray.direction.dot(hrec.n) / ray.direction.length();
    }

    srec.albedo = Vec3<float>(1.0f);

    float reflect_prob;
    Vec3<float> refracted;
    if (refract(-ray.direction, n, ratio, refracted)) {
      reflect_prob = schlick(cosine, ri);
    } else {
      reflect_prob = 1.0f;
    }

    if (randf() < reflect_prob) {
      srec.ray = Ray{hrec.p, reflected};
    } else {
      srec.ray = Ray{hrec.p, refracted};
    }

    return true;
  }
};

struct Shape {
  virtual bool hit(const Ray& ray, float t0, float t1, HitRec& rec) = 0;
};

struct Sphere : Shape {
  Vec3<float> center;
  float radius;
  Material* material;
  Sphere(const Vec3<float>& center, float radius, Material* material)
      : center{center}, radius{radius}, material{material} {}
  virtual bool hit(const Ray& ray, float t0, float t1, HitRec& rec) override {
    Vec3<float> co = ray.origin - center;
    float a = ray.direction.dot(ray.direction);
    float b = 2.0f * ray.direction.dot(co);
    float c = co.dot(co) - radius * radius;
    float d = b * b - 4.0f * a * c;
    if (d > 0) {
      float t = (-b - std::sqrt(d)) / (2.0f * a);
      if (t0 < t && t < t1) {
        rec.t = t;
        rec.p = ray.at(t);
        rec.n = (rec.p - center) / radius;
        rec.material = material;
        return true;
      }

      t = (-b + std::sqrt(d)) / (2.0f * a);
      if (t0 < t && t < t1) {
        rec.t = t;
        rec.p = ray.at(t);
        rec.n = (rec.p - center) / radius;
        rec.material = material;
        return true;
      }
    }
    return false;
  }
};

struct Triangle : Shape {
  Vec3<float> p1, p2, p3;
  Triangle(const Vec3<float>& p1, const Vec3<float>& p2, const Vec3<float>& p3)
      : p1{p1}, p2{p2}, p3{p3} {}

  float denom(const Vec3<float>& a, const Vec3<float>& b,
              const Vec3<float>& c) {
    return a.x * b.y * c.z + a.y * b.z * c.x + a.z * b.x * c.y -
           a.x * b.z * c.y - a.y * b.x * c.z - a.z * b.y * c.x;
  }

  virtual bool hit(const Ray& ray, float t0, float t1, HitRec& rec) override {
    Vec3<float> edge1 = p2 - p1;
    Vec3<float> edge2 = p3 - p1;
    float u = 0.0f, v = 0.0f;
    Vec3<float> v1 = p1 + u * edge1 + v * edge2;
    Vec3<float> v2 = ray.origin + t0 * ray.direction;
    float den = denom(edge1, edge2, -ray.direction);
    return false;
  }
};

struct Camera {
  Vec3<float> origin;
  Vec3<float> view[3];

  Camera() = default;
  Camera(const Vec3<float>& lookFrom, const Vec3<float>& lookAt,
         const Vec3<float>& vUp) {
    float vFov = 50.0f, aspect = float(W) / float(H);  // 3 / 2
    origin = lookFrom;
    Vec3<float> w = (lookFrom - lookAt).normalize();
    Vec3<float> u = vUp.cross(w).normalize();
    Vec3<float> v = w.cross(u);
    float halfH = tanf(3.141592f * vFov / 360.0f);
    float halfW = aspect * halfH;
    view[2] = origin - halfW * u - halfH * v - w;
    view[0] = 2.0f * halfW * u;
    view[1] = 2.0f * halfH * v;
  }

  Ray getRay(float u, float v) {
    return {origin, u * view[0] + v * view[1] + view[2] - origin};
  }
};

struct ShapeList : Shape {
  std::vector<Shape*> shapes;
  ShapeList() {}
  void add(Shape* shape) { shapes.push_back(shape); }
  virtual bool hit(const Ray& ray, float t0, float t1, HitRec& rec) override {
    HitRec rec2;
    bool hitAnything = false;
    float far = t1;
    for (auto& shape : shapes) {
      if (shape->hit(ray, t0, far, rec2)) {
        hitAnything = true;
        far = rec2.t;
        rec = rec2;
      }
    }
    return hitAnything;
  }
};

ShapeList shapes;

void init() {
  float scale = 1000.0f;
  auto red = new Lambertian({0.8f, 0.0f, 0.0f});
  auto blue = new Lambertian({0.1f, 0.2f, 0.5f});
  auto green = new Lambertian({0.1f, 0.5f, 0.2f});
  auto white = new Lambertian({0.9f, 0.9f, 0.9f});
  auto metal = new Metal(Vec3<float>{0.8f});
  auto metal1 = new Metal(Vec3<float>{0.7f, 0.5f, 0.5f});
  auto metal2 = new Metal(Vec3<float>{0.5f, 0.7f, 0.5f});
  auto metal3 = new Metal(Vec3<float>{0.5f, 0.5f, 0.7f});
  auto ground = new Lambertian({0.3f, 0.8f, 0.0f});
  auto light = new DiffuseLight({0.8f, 0.8f, 0.8f});
  auto skylight = new DiffuseLight({0.5f, 0.7f, 1.0f});
  auto dielectric = new Dielectric(1.5f);

  // ground
  shapes.add(new Sphere(Vec3<float>{0.0f, -scale - 0.3f, 0.0f}, scale, ground));
  // ceiling
  shapes.add(new Sphere(Vec3<float>{0.0f, scale + 3.0f, 0.0f}, scale, light));
  // right wall
  shapes.add(new Sphere(Vec3<float>{-scale - 3.0f, 0.0f, 0.0f}, scale, red));
  // left wall
  shapes.add(new Sphere(Vec3<float>{scale + 3.0f, 0.0f, 0.0f}, scale, blue));
  // front wall
  shapes.add(new Sphere(Vec3<float>{0.0f, 0.0f, scale + 7.0f}, scale, white));
  // back wall
  shapes.add(new Sphere(Vec3<float>{0.0f, 0.0f, -scale - 7.0f}, scale, white));

  std::vector<Material*> materials = {red,    blue,   green,  white,     metal,
                                      metal1, metal2, metal3, dielectric};
  for (int i = 0; i < 12; ++i) {
    shapes.add(new Sphere(Vec3<float>{5.0f * randf() - 2.5f, 2.0f * randf(),
                                      5.0f * randf() - 2.5f},
                          0.3f, materials[i % materials.size()]));
  }
  shapes.add(new Sphere(Vec3<float>{-1.5f, 0.1f, -1.5f}, 0.32f, dielectric));
}

Vec3<float> color(const Ray& ray, int depth = 0) {
  if (depth > 50) return Vec3<float>{0.000001f};
  HitRec rec;
  if (shapes.hit(ray, 0.001f, FLT_MAX, rec)) {
    Vec3<float> emitted = rec.material->emit(ray, rec);
    ScatterRec srec;
    if (rec.material->scatter(ray, rec, srec)) {
      Vec3<float> c = color(srec.ray, depth + 1);
      return {
          emitted.x + srec.albedo.x * c.x,
          emitted.y + srec.albedo.y * c.y,
          emitted.z + srec.albedo.z * c.z,
      };
    } else {
      return emitted;
    }
  }
  auto v = ray.direction.normalize();
  float t = 0.5f * (v.y + 1.0f);
  return Vec3<float>::lerp({0.9f}, {0.5f, 0.7f, 1.0f}, t);
}

uint8_t clip(float v) {
  v *= 255;
  if (v > 255) return 255;
  if (v < 0) return 0;
  return v;
}

int main() {
  init();

  Bitmap bmp(W, H);
  Camera camera({-2.0f, 1.5f, -5.0f}, {0.2f, 0.0f, -0.2f}, {0, -1, 0});

  // render
  std::vector<std::future<void>> futures;
  auto render = [&](int x1, int x2) {
    for (int y = 0; y < H; ++y) {
      for (int x = x1; x < x2; ++x) {
        Vec3<float> c{0.0f};
        for (int i = 0; i < samples; ++i) {
          float u = float(x + randf()) / float(W);
          float v = float(y + randf()) / float(H);
          Ray ray = camera.getRay(u, v);
          c = c + color(ray);
        }
        c /= samples;
        bmp.setColor(x, y, Color{clip(c.x), clip(c.y), clip(c.z)});
      }
    }
  };

  int w = W / threadCount;
  for (int threadID = 0; threadID < threadCount; ++threadID) {
    futures.push_back(
        std::async(std::launch::async, render, threadID * w, threadID * w + w));
  }

  for (int threadID = 0; threadID < threadCount; ++threadID)
    futures[threadID].get();

  bmp.save("result_async.bmp");
}
