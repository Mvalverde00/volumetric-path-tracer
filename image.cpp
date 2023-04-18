#include "image.h"

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

// Taken from PBRT
inline float gamma_correct(float value) {
  if (value <= 0.0031308f) return 12.92f * value;
  return 1.055f * std::pow(value, (float)(1.f / 2.4f)) - 0.055f;
}

Image::Image(int w, int h) : width(w), height(h), pixels(NULL) {}

Image::~Image() {
  if (pixels != NULL) {
    delete[] pixels;
  }
}

void Image::init() {
  if (width * height > 0 && pixels == NULL) {
    pixels = new Color[width * height];
  }
}

void Image::set_pixel(int x, int y, Color c) {
  if (x >= width || y >= height || x < 0 || y < 0) {
    std::cerr << "Tried to write pixel out of bounds";
    return;
  }

  pixels[x + y * width] = c;
}

Color Image::get_pixel(int x, int y) {
  if (x >= width || y >= height || x < 0 || y < 0) {
    std::cerr << "Tried to read pixel out of bounds";
    return Color(0.0, 0.0, 0.0);
  }

  return pixels[x + y * width];
}

void Image::write_to_file(std::string filename) {
  int channels = 3; // r, g, b
  uint8_t* rgb_pixels = new uint8_t[width * height * channels];
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      Color c = get_pixel(x, y);
      c = glm::clamp(c, 0.0f, 1.0f);
      //std::cout << x << ", " << y << ";  color " << c[0] << "/" << c[1] << "/" << c[2] << "\n";
      /*
      uint8_t r = (uint8_t)(c.r * 255);
      uint8_t g = (uint8_t)(c.g * 255);
      uint8_t b = (uint8_t)(c.b * 255);
      */
      /*
      float gamma_c = 1.0 / 2.5f;
      uint8_t r = (uint8_t)(std::powf(c.r, gamma_c) * 255);
      uint8_t g = (uint8_t)(std::powf(c.g, gamma_c) * 255);
      uint8_t b = (uint8_t)(std::powf(c.b, gamma_c) * 255);
      */
      uint8_t r = (uint8_t)(gamma_correct(c.r) * 255.0f + 0.5f);
      uint8_t g = (uint8_t)(gamma_correct(c.g) * 255.0f + 0.5f);
      uint8_t b = (uint8_t)(gamma_correct(c.b) * 255.0f + 0.5f);


      // Invert y axis when drawing.
      int draw_y = height - y - 1;
      rgb_pixels[3 * (x + draw_y * width)] = r;
      rgb_pixels[3 * (x + draw_y * width) + 1] = g;
      rgb_pixels[3 * (x + draw_y * width) + 2] = b;
    }
  }

  stbi_write_png(filename.c_str(), width, height, channels, rgb_pixels, width * channels);

  delete[] rgb_pixels;
}