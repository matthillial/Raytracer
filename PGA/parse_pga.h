
//Set the global scene parameter variables
//TODO: Set the scene parameters based on the values in the scene file

#ifndef PARSE_PGA_H
#define PARSE_PGA_H

#include <cstdio>
#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>

#include "image_lib.h" 

#define MAXLINE 256

//using namespace std;

//Camera & Scene Parmaters (Global Variables)
//Here we set default values, override them in parseSceneFile()

//Image Parmaters
int img_width = 640, img_height = 480;
std::string imgName = "raytraced.bmp";
Color background = Color(0,0,0);



//Camera Parmaters
Point3D eye = Point3D(0,0,0); 
Dir3D forward = Dir3D(0,0,1).normalized();
Dir3D up = Dir3D(0,1,0).normalized();
Dir3D right = Dir3D(-1,0,0).normalized();
float halfAngleVFOV = 45; 

//Scene (Sphere) Parmaters
Point3D spherePos = Point3D(0,0,2);
float sphereRadius = 1; 

//Material Parameters
float ar = 0;
float ag = 0;
float ab = 0;
float dr = 1;
float dg = 1;
float db = 1;
float sr = 0;
float sg = 0;
float sb = 0;
float ns = 5;
float tr = 0;
float tg = 0;
float tb = 0;
float ior = 1;



//Lighting Parameters
float pointlight = 0;  		//point light
float plr = 0; 
float plg = 0;
float plb = 0;
float plx = 0;
float ply = 0;
float plz = 0;

float ambr = 0;				//ambient light (default 0)
float ambg = 0;
float ambb = 0;





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
      spherePos = Point3D(x, y, z);
      sphereRadius = r;
    }
    else if (word == "film_resolution:") {
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
	else if (word == "background:") {
	  float r, g, b;
	  ss >> r;
	  ss >> g;
	  ss >> b;
	  background = Color(r,g,b);
	}
	else if (word == "material:") {
	  ss >> ar;
	  ss >> ag;
	  ss >> ab;
	  ss >> dr;
	  ss >> dg;
	  ss >> db;
	  ss >> sr;
	  ss >> sg;
	  ss >> sb;
	  ss >> ns;
	  ss >> tr;
	  ss >> tg;
	  ss >> tb;
	  ss >> ior;
		
	}
	else if (word == "point_light:") {
	  pointlight = 1;
	  ss >> plr;
	  ss >> plg;
	  ss >> plb;
	  ss >> plx;
	  ss >> ply;
	  ss >> plz;
	}
	else if (word == "ambient_light:") {
	  ss >> ambr;
	  ss >> ambg;
	  ss >> ambb;
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
