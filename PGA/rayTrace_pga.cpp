//CSCI 5607 HW3 - Rays & Files
//This HW has three steps:
// 1. Compile and run the program (the program takes a single command line argument)
// 2. Understand the code in this file (rayTrace_pga.cpp), in particular be sure to understand:
//     -How ray-sphere intersection works
//     -How the rays are being generated
//     -The pipeline from rays, to intersection, to pixel color
//    After you finish this step, and understand the math, take the HW quiz on canvas
// 3. Update the file parse_pga.h so that the function parseSceneFile() reads the passed in file
//     and sets the relevant global variables for the rest of the code to product to correct image

//To Compile: g++ -fsanitize=address -std=c++11 rayTrace_pga.cpp

//For Visual Studios
#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

//Images Lib includes:
#define STB_IMAGE_IMPLEMENTATION //only place once in one .cpp file
#define STB_IMAGE_WRITE_IMPLEMENTATION //only place once in one .cpp files
#include "image_lib.h" //Defines an image class and a color class

//#3D PGA
#include "PGA_3D.h"

//High resolution timer
#include <chrono>

//Scene file parser
#include "parse_pga.h"

#include <limits>
#define INF std::numeric_limits<float>::infinity()

bool raySphereIntersect_fast(Point3D rayStart, Line3D rayLine, Point3D sphereCenter, float sphereRadius){
  Dir3D dir = rayLine.dir();
  float a = dot(dir,dir);
  Dir3D toStart = (rayStart - sphereCenter);
  float b = 2 * dot(dir,toStart);
  float c = dot(toStart,toStart) - sphereRadius*sphereRadius;
  float discr = b*b - 4*a*c;
  if (discr < 0) return false;
  else{
    float t0 = (-b + sqrt(discr))/(2*a);
    float t1 = (-b - sqrt(discr))/(2*a);
    if (t0 > 0 || t1 > 0) return true;
  }
  return false;
}

float raySphereIntersect(Point3D rayStart, Line3D rayLine, Point3D sphereCenter, float sphereRadius){
  Point3D projPoint = dot(rayLine,sphereCenter)*rayLine;      //Project to find closest point between circle center and line [proj(sphereCenter,rayLine);]
  float distSqr = projPoint.distToSqr(sphereCenter);          //Point-line distance (squared)
  float d = distSqr/(sphereRadius*sphereRadius);             //If distance is larger than radius, then...
  if (d > 1) return INF;                                     //... the ray missed the sphere
  float w = sphereRadius*sqrt(1-d);                          //Pythagorean theorem to determine dist between proj point and intersection points
  Point3D p1 = projPoint - rayLine.dir()*w;                   //Add/subtract above distance to find hit points
  Point3D p2 = projPoint + rayLine.dir()*w; 
  float d1 = (p1 - rayStart).magnitude();
  float d2 = (p2 - rayStart).magnitude();

  if (dot((p1-rayStart),rayLine.dir()) >= 0 ) return d1;     //Is the first point in same direction as the ray line?
  if (dot((p2-rayStart),rayLine.dir()) >= 0 ) return d2;     //Is the second point in same direction as the ray line?
  return INF;
}

int main(int argc, char** argv){
  
  //Read command line paramaters to get scene file
  if (argc != 2){
     std::cout << "Usage: ./a.out scenefile\n";
     return(0);
  }
  std::string secenFileName = argv[1];

  //Parse Scene File
  parseSceneFile(secenFileName);

  float imgW = img_width, imgH = img_height;
  float halfW = imgW/2, halfH = imgH/2;
  float d = halfH / tanf(halfAngleVFOV * (M_PI / 180.0f));

  Image outputImg = Image(img_width,img_height);
  auto t_start = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < img_width; i++){
    for (int j = 0; j < img_height; j++){
      //TODO: In what way does this assumes the basis is orthonormal?
      float u = (halfW - (imgW)*((i+0.5)/imgW));
      float v = (halfH - (imgH)*((j+0.5)/imgH));
      Point3D p = eye - d*forward + u*right + v*up;
      Dir3D rayDir = (p - eye); 
      Line3D rayLine = vee(eye,rayDir).normalized();  //Normalizing here is optional

      // Calculate closest intersection
      float closeD = INF; // Distance to nearest sphere
      int closeIndex = -1; // Index of nearest sphere
      for(int x = 0; x < spheres.size(); x++) {
        float d = raySphereIntersect(eye,rayLine,spheres[x].pos,spheres[x].r);
        if (d < closeD) {
          closeD = d;
          closeIndex = x;
        }
      }

      Color color;
      if (closeIndex >= 0) {
        Point3D closePoint = eye + (closeD - 0.01) * rayLine.dir();
        Dir3D normal = closePoint - spheres[closeIndex].pos;
        Material material = materials[spheres[closeIndex].matIndex];
        Color mAmbient = material.ambient;
        Color mDiffuse = material.diffuse;
        Color mSpecular = material.specular;
        Color mTransmissive = material.transmissive;

        // Set composite color to ambient light to start
        color = Color(mAmbient.r * ambient.r, mAmbient.g * ambient.g, mAmbient.b * ambient.b);

        // Check directional lights
        for(int lightI = 0; lightI < directionalLights.size(); lightI++) {
          DirectionalLight currentLight = directionalLights[lightI];
          Dir3D lightDir = -1 * currentLight.dir;
          Line3D lightLine = vee(closePoint, lightDir).normalized();
          bool blocked = false;
          for(int sphereI = 0; sphereI < spheres.size(); sphereI ++) {
            Sphere currentSphere = spheres[sphereI];
            if (raySphereIntersect(closePoint, lightLine, currentSphere.pos, currentSphere.r) < INF) {
              blocked = true;
              break;
            }
          }

          // Add light into composite color if it isn't blocked
          if (!blocked) {
            float difdot = (0 < dot(lightDir, normal)) ? dot(normal, lightDir) : 0;
            float r = color.r + (currentLight.color.r * mDiffuse.r * difdot);
            float g = color.g + (currentLight.color.g * mDiffuse.g * difdot);
            float b = color.b + (currentLight.color.b * mDiffuse.b * difdot);
            color = Color(r, g, b);
          }
        }
      }
      else {
        color = background;
      }
      
      outputImg.setPixel(i,j, color);
      //outputImg.setPixel(i,j, Color(fabs(i/imgW),fabs(j/imgH),fabs(0))); //TODO: Try this, what is it visualizing?
    }
  }
  auto t_end = std::chrono::high_resolution_clock::now();
  printf("Rendering took %.2f ms\n",std::chrono::duration<double, std::milli>(t_end-t_start).count());

  outputImg.write(imgName.c_str());
  return 0;
}
