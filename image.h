#pragma once

#include <string>
#include <iostream>
#include <glm/glm.hpp>

using Color = glm::vec3;

class Image {

  int width;
  int height;

  /* 1D Array of pixels of length width * height.  Can be thought of as 2D
    array with clever indexing. */
  Color* pixels;

public:
  Image() : width(0), height(0), pixels(NULL) {};
  Image(int w, int h);

  ~Image(); // Needs to free the pixels

  void init();
 
  void set_pixel(int x, int y, Color c);
  Color get_pixel(int x, int y);
  void write_to_file(std::string filename);

  int get_width() { return width; };
  int get_height() { return height; };
};
