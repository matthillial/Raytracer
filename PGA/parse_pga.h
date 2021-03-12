
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

#include "image_lib.h" 

#define MAXLINE 256

//using namespace std;
using std::vector;

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

//Materials
class Material {
  public:
  Color ambient;
  Color diffuse;
  Color specular;
  float ns;
  Color transmissive;
  float ior;
};

vector<Material> materials = {Material{Color(1,1,1),Color(1,1,1),Color(0,0,0),5, Color(0,0,0), 1}};
int currentMaterial = 0;
// materials.push_back(Material{0,0,0,1,1,1,0,0,0,0,0,0});

//Scene (Sphere) Parmaters
class Sphere {
  public:
  Point3D pos;
  float r;
  int matIndex;
};

vector<Sphere> spheres;

//Light Stuff
class DirLight {
  public:
  Point3D dir;
  Color color;
};

class PointLight {
  public:
  Point3D pos;
  Color color;
};

vector<DirLight> directionalLights;
vector<PointLight> pointLights;
Color ambient = Color(0, 0, 0);

//Material Parameters


float ar = 0;	//ambient
float ag = 0;
float ab = 0;

float dr = 1;	//diffuse
float dg = 1;
float db = 1;

float sr = 0;	//specular
float sg = 0;
float sb = 0;

float ns = 5;	//phong cosine power for specular

float tr = 0;	//transmissive
float tg = 0;
float tb = 0;	

float ior = 1;	//index of refraction



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


float dirr = 0;
float dirg = 0;
float dirb = 0;
Dir3D dirdir;

int max_depth = 2;



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
      spheres.push_back(Sphere{Point3D(x, y, z), r, currentMaterial});
      printf("Material: %d\n", spheres[spheres.size()-1].matIndex);
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
    printf("Material set");
    currentMaterial += 1;
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
		materials.push_back(Material{Color(ar, ag, ab), Color(dr, dg, db), Color(sr, sg, sb), ns, Color(tr, tg, tb), ior});
	}
	else if (word == "point_light:") {
	  pointlight = 1;
	  ss >> plr;
	  ss >> plg;
	  ss >> plb;
	  ss >> plx;
	  ss >> ply;
	  ss >> plz;
    pointLights.push_back(PointLight{Point3D(plx, ply, plz), Color(plr, plg, plb)});
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
      directionalLights.push_back(DirLight{Point3D(x, y, z), Color(r, g, b)});
    }
  else if (word == "max_depth:") {
    int d;
    ss >> d;
    printf("%d", d);
    max_depth = d;
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
