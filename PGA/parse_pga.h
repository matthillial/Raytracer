
//Set the global scene parameter variables
//TODO: Set the scene parameters based on the values in the scene file

#ifndef PARSE_PGA_H
#define PARSE_PGA_H

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <vector>

#define MAXLINE 256

//using namespace std;
using std::vector;

//Camera & Scene Parmaters (Global Variables)
//Here we set default values, override them in parseSceneFile()

//Image Parmaters
int img_width = 800, img_height = 600;
std::string imgName = "raytraced.png";



//Camera Parmaters
Point3D eye = Point3D(0,0,0); 
Dir3D forward = Dir3D(0,0,-1).normalized();
Dir3D up = Dir3D(0,1,0).normalized();
Dir3D right = Dir3D(-1,0,0).normalized();
float halfAngleVFOV = 35; 

//Scene (Sphere) Parmaters
class Sphere {
  public:
  Point3D pos;
  float r;
};

vector<Sphere> spheres;

//Light Stuff
class Light {
  public:
  Point3D pos;
  Color color;
};

vector<Light> directionalLights;
Color ambient = Color(0, 0, 0);

void parseSceneFile(std::string fileName){
  //TODO: Override the default values with new data from the file "fileName"
  std::ifstream inFile(fileName, std::ios::in);
  char line[MAXLINE];
  char subLine[20];

  while (inFile.getline(line, MAXLINE)) {
    if (line[0] == '#') {
      continue;
    } 
    
    std::string str(line);
    std::istringstream ss(str);
    std::string word;
    ss >> word;

    if (word == "sphere:") {
      float x, y, z, r;
      ss >> x;
      ss >> y;
      ss >> z;
      ss >> r;
      spheres.push_back(Sphere{Point3D(x, y, z), r});
    }
    else if (word == "ambient_light:") {
      float r, g, b;
      ss >> r;
      ss >> g;
      ss >> b;
      ambient = Color(r, g, b);
    }
    else if (word == "directional_light:") {
      float x, y, z, r, g, b;
      ss >> r;
      ss >> g;
      ss >> b;
      ss >> x;
      ss >> y;
      ss >> z;
      pointLights.push_back(Light{Point3D(x, y, z), Color(r, g, b)});
    }

    else if (word == "image_resolution:") {
      int w, h;
      ss >> w;
      ss >> h;
      img_width = w;
      img_height = h;
    }
    else if (word == "output_image:") {
      ss >> imgName;
    }
    else if (word == "camera_pos:") {
      float x, y, z;
      ss >> x;
      ss >> y;
      ss >> z;
      eye = Point3D(x, y, z);
    }
    else if (word == "camera_fwd:") {
      float fx, fy, fz;
      ss >> fx;
      ss >> fy;
      ss >> fz;
      forward = Dir3D(fx, fy, fz).normalized();
    }
    else if (word == "camera_up:") {
      float ux, uy, uz;
      ss >> ux;
      ss >> uy;
      ss >> uz;
      up = Dir3D(ux, uy, uz).normalized();
    }
    else if (word == "camera_fov_ha:") {
      float ha;
      ss >> ha;
      halfAngleVFOV = ha;
    }
  }

  if (abs(dot(up, forward)) == 1) {
    printf("ERROR: Up and forward vectors are parallel, no orthogonal basis can be found");
    Dir3D forward = Dir3D(0,0,-1).normalized();
    Dir3D up = Dir3D(0,1,0).normalized();
  }

  right = cross(up, forward).normalized();
  up = cross(forward, right).normalized();

  //TODO: Create an orthagonal camera basis, based on the provided up and right vectors
  printf("Orthagonal Camera Basis:\n");
  forward.print("forward");
  right.print("right");
  up.print("up");
}

#endif