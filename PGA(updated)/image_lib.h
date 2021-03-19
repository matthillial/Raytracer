//CSCI 5607 Header-only Image Library

#ifndef IMAGE_LIB_H
#define IMAGE_LIB_H

#include "stb_image.h"
#include "stb_image_write.h"
#include <stdint.h> //for uint8_t

struct Color{
  float r,g,b;

  Color(float r, float g, float b) : r(r), g(g), b(b) {}
  Color() : r(0), g(0), b(0) {}
};

struct Image{
  int width, height;
  Color* pixels;

  Image(int w, int h) : width(w), height(h) {
    pixels = new Color[width*height];
  }

  //Copy constructor - Called on: Image img1 = img2; //Making new image
  Image(const Image& cpy){
    width = cpy.width;
    height = cpy.height;
    memcpy(pixels, cpy.pixels, width*height*sizeof(Color));
  }

  //Assignment operator - Called on:  img1 = img2; //Overwriting existing data
  Image&  operator=(const Image& rhs){
    width = rhs.width;
    height = rhs.height;
    memcpy(pixels, rhs.pixels, width*height*sizeof(Color));
    return *this;
  }

  //Creat an image from a file
  Image(const char* fname) {
    int numComponents; //(e.g., Y, YA, RGB, or RGBA)
    unsigned char *data = stbi_load(fname, &width, &height, &numComponents, 4);
    
    if (data == NULL){
      printf("Error loading image: '%s'\n", fname);
      exit(-1);
    }

    pixels = new Color[width*height];

    for (int i = 0; i < width; i++){
      for (int j = 0; j < height; j++){
        pixels[(i+j*width)] = Color(data[4*(i+j*width)+0]/255.0,data[4*(i+j*width)+1]/255.0,data[4*(i+j*width)+2]/255.0);
      }
    }

    stbi_image_free(data);
  }

  void setPixel(int i, int j, Color c){
    pixels[i+j*width] = c;
  }

  Color& getPixel(int i, int j){
    return pixels[i+j*width];
  }

  uint8_t* toBytes(){
    uint8_t* rawPixels = new uint8_t[width*height*4];
    for (int i = 0; i < width; i++){
      for (int j = 0; j < height; j++){
        Color col = getPixel(i,j);
        rawPixels[4*(i+j*width)+0] = uint8_t(fmin(col.r,1)*255);
        rawPixels[4*(i+j*width)+1] = uint8_t(fmin(col.g,1)*255);
        rawPixels[4*(i+j*width)+2] = uint8_t(fmin(col.b,1)*255);
        rawPixels[4*(i+j*width)+3] = 255; //alpha
      }
    }
    return rawPixels;
  }

  void write(const char* fname){

    uint8_t* rawBytes = toBytes();
    
    int lastc = strlen(fname);

    switch (fname[lastc-1]){
      case 'g': //jpeg (or jpg) or png
        if (fname[lastc-2] == 'p' || fname[lastc-2] == 'e') //jpeg or jpg
            stbi_write_jpg(fname, width, height, 4, rawBytes, 95);  //95% jpeg quality
        else //png
            stbi_write_png(fname, width, height, 4, rawBytes, width*4);
        break;
      case 'a': //tga (targa)
        stbi_write_tga(fname, width, height, 4, rawBytes);
        break;
      case 'p': //bmp
      default:
        stbi_write_bmp(fname, width, height, 4, rawBytes);
    }

    delete[] rawBytes;
  }

  ~Image(){delete[] pixels;}
};

#endif