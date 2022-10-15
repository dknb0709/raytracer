#pragma once

#include <fstream>

static const unsigned int BITMAP_FILE_HEADER_SIZE = 14;
static const unsigned int BITMAP_INFO_HEADER_SIZE = 40;
static const unsigned int BITMAP_HEADER_SIZE =
    BITMAP_FILE_HEADER_SIZE + BITMAP_INFO_HEADER_SIZE;

#pragma pack(2)
struct BitmapFileHeader {
  unsigned short bfType;
  unsigned int bfSize;
  unsigned short bfReserved1;
  unsigned short bfReserved2;
  unsigned int bfOffBits;
};

#pragma pack()
struct BitmapInfoHeader {
  unsigned int biSize;
  unsigned int biWidth;
  unsigned int biHeight;
  unsigned short biPlanes;
  unsigned short biBitCount;
  unsigned int biCompression;
  unsigned int biImageSize;
  unsigned int biXPixPerMeter;
  unsigned int biYPixPerMeter;
  unsigned int biColorUsed;
  unsigned int biColorImportant;
};

struct Color {
  uint8_t r, g, b;
};

struct Image {
  Image(int w, int h) : w{w}, h{h} { data = new Color[w * h]; }
  virtual ~Image() { delete data; }

  Color getColor(int x, int y) { return data[x + w * y]; }
  void setColor(int x, int y, Color c) { data[x + w * y] = c; }

  Color* data;
  int w, h;
};

class Bitmap {
 public:
  Bitmap(const std::string& path) { load(path); }
  Bitmap(int w, int h) { m_image = new Image(w, h); }
  virtual ~Bitmap() { delete m_image; }

  void load(const std::string& path) {
    std::ifstream fin(path, std::ios::binary);
    if (!fin.is_open()) return;
    // load header
    char headerBuffer[BITMAP_HEADER_SIZE];
    fin.read(headerBuffer, sizeof(char) * BITMAP_HEADER_SIZE);
    BitmapFileHeader* fileHeader = (BitmapFileHeader*)headerBuffer;
    BitmapInfoHeader* infoHeader =
        (BitmapInfoHeader*)(headerBuffer + BITMAP_FILE_HEADER_SIZE);

    // load image data
    int w = infoHeader->biWidth;
    int h = infoHeader->biHeight;
    int bpp = infoHeader->biBitCount;
    int compression = infoHeader->biCompression;
    if (w * h > 33177600 || bpp != 24 || compression != 0) {
      return;
    }
    m_image = new Image(w, h);
    int stride = (w * 3 + 3) / 4 * 4;
    char row[stride];
    for (int y = 0; y < h; ++y) {
      fin.read(row, sizeof(char) * stride);
      for (int x = 0; x < w; ++x) {
        Color& pixel = m_image->data[x + w * (h - 1 - y)];
        pixel.b = row[x * 3 + 0];
        pixel.g = row[x * 3 + 1];
        pixel.r = row[x * 3 + 2];
      }
    }
  }

  void save(const std::string& path) {
    std::ofstream fout(path, std::ios::binary | std::ios::trunc);
    if (!fout.is_open()) return;
    // generate header buffer
    char headerBuffer[BITMAP_HEADER_SIZE];
    BitmapFileHeader* fileHeader = (BitmapFileHeader*)headerBuffer;
    BitmapInfoHeader* infoHeader =
        (BitmapInfoHeader*)(headerBuffer + BITMAP_FILE_HEADER_SIZE);

    int stride = (3 * m_image->w + 3) / 4 * 4;
    int imageSize = stride * m_image->h;
    int fileSize = imageSize + BITMAP_HEADER_SIZE;

    fileHeader->bfType = 0x4D42;
    fileHeader->bfSize = fileSize;
    fileHeader->bfReserved1 = 0x00;
    fileHeader->bfReserved2 = 0x00;
    fileHeader->bfOffBits = BITMAP_HEADER_SIZE;
    infoHeader->biSize = BITMAP_INFO_HEADER_SIZE;
    infoHeader->biWidth = m_image->w;
    infoHeader->biHeight = m_image->h;
    infoHeader->biPlanes = 0;
    infoHeader->biBitCount = 24;
    infoHeader->biCompression = 0;
    infoHeader->biImageSize = imageSize;
    infoHeader->biXPixPerMeter = 0;
    infoHeader->biYPixPerMeter = 0;
    infoHeader->biColorUsed = 0;
    infoHeader->biColorImportant = 0;

    // write header
    fout.write(headerBuffer, sizeof(char) * BITMAP_HEADER_SIZE);

    // write image data
    char row[stride];
    for (int y = 0; y < m_image->h; ++y) {
      for (int x = 0; x < m_image->w; ++x) {
        Color c = m_image->data[x + m_image->w * (m_image->h - 1 - y)];
        row[3 * x + 0] = c.b;
        row[3 * x + 1] = c.g;
        row[3 * x + 2] = c.r;
      }
      fout.write(row, stride);
    }
  }

  int getWidth() { return m_image->w; }

  int getHeight() { return m_image->h; }

  Color getColor(int x, int y) {
    int w = m_image->w;
    int h = m_image->h;
    if (x < 0 || x > w || y < 0 || y > h) return Color{0, 0, 0};
    return m_image->getColor(x, y);
  }

  void setColor(int x, int y, Color c) {
    int w = m_image->w;
    int h = m_image->h;
    if (x < 0 || x > w || y < 0 || y > h) return;
    m_image->setColor(x, y, c);
  }

 private:
  Image* m_image;
};
