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


float rayPlaneIntersect(Point3D rayStart, Line3D rayLine, Point3D v1, Point3D v2, Point3D v3){
  Plane3D triPlane = vee(v1, v2, v3);	//make plane
  Point3D hitPoint = Point3D(wedge(rayLine,triPlane));
  if ((dot(vee(v1,v2,hitPoint),triPlane) >= 0) && (dot(vee(v2,v3,hitPoint),triPlane) >= 0) && (dot(vee(v3,v1,hitPoint), triPlane) >= 0)) {
    float d = (hitPoint - rayStart).magnitude();
    return d;
  }
  return INF;
}

Color rayCast(Point3D start, Dir3D rayDir, Line3D rayLine, int recursionDepth, bool debug) {
	  //variables for differentiating between sphere and triangle and triangle(normals specified)
	  float dif = -1;
		
	
      // Calculate closest intersection for spheres
      float closeD = INF; // Distance to nearest sphere
      int closeIndex = -1; // Index of nearest sphere
      for(int x = 0; x < spheres.size(); x++) {
        float d = raySphereIntersect(start,rayLine,spheres[x].pos,spheres[x].r);
        if (d < closeD) {
          closeD = d;
          closeIndex = x;
        }
      }
      
      
      // Calculate closest intersection for triangles(non specified normals)
      float closeDTriangle = INF;	// Distance to nearest triangle
      int closeTriangleIndex = -1;  // Index of nearest triangle
      for(int x = 0; x < triangles.size(); x++) {
		 float d = rayPlaneIntersect(start, rayLine, vertex[triangles[x].v1], vertex[triangles[x].v2], vertex[triangles[x].v3]);
		 if (d < closeDTriangle) {
		   closeDTriangle = d;
		   closeTriangleIndex = x;
		 }
	  }
	  
	  // Calculate closest interesection for triangles(specified normals)
	  float closeDNormalTriangle = INF;	// Distance to nearest triangle
      int closeNormalTriangleIndex = -1;  // Index of nearest triangle
      for(int x = 0; x < normaltriangles.size(); x++) {
		 float d = rayPlaneIntersect(start, rayLine, vertex[normaltriangles[x].v1], vertex[normaltriangles[x].v2], vertex[normaldtriangles[x].v3]);
		 if (d < closeDTriangle) {
		   closeDNormalTriangle = d;
		   closeNormalTriangleIndex = x;
		 }
	  }
	  
	  
	  // Calculate which intersection is the closest between all the primitives || dif = 0 for spheres || dif = 1 for triangle (normals not specified) || dif = 2 for traingles (normals specified)
	  if (closeD < closeDTriangle && closeD < closeDNormalTriangle) {
		dif = 0;
	  }
	  else if (closeDTriangle < closeD && closeDTriangle < closeDNormalTriangle) {
		dif = 1;
	  }
	  else if (closeDNormalTriangle < closeD && closeDNormalTriangle < closeDTriangle) {
		dif = 2;
	  }
	  
	  
      
      
      

      // Calculate color of point for sphere
      Color color;
      if (closeIndex >= 0 && dif == 0) {
        Point3D closePoint = start + (closeD - 0.001) * rayLine.dir();
        Point3D insidePoint = start + (closeD + 0.0001) * rayLine.dir();
        Dir3D normal = (closePoint - spheres[closeIndex].pos).normalized();

        Material material = materials[spheres[closeIndex].matIndex];
        Color mAmbient = material.ambient;
        Color mDiffuse = material.diffuse;
        Color mSpecular = material.specular;
        Color mTransmissive = material.transmissive;

        // Set composite color to ambient light to start
        color = Color(mAmbient.r * ambient.r, mAmbient.g * ambient.g, mAmbient.b * ambient.b);

        // Check directional lights
        for(int lightI = 0; lightI < directionalLights.size(); lightI++) {
          DirLight currentLight = directionalLights[lightI];
          Dir3D lightDir = -1 * currentLight.dir;
          Line3D lightLine = vee(closePoint, lightDir).normalized();

          // Check all spheres to see if light is blocked
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
            Dir3D viewDir = (start - closePoint).normalized();
            Dir3D halfway = (lightDir + viewDir).normalized();
            float difdot = (0 < dot(lightDir, normal)) ? dot(normal, lightDir) : 0;
            float specdot = (0 < pow(dot(halfway,normal),ns)) ? pow(dot(halfway,normal),ns) : 0;

            float r = color.r + (currentLight.color.r * mDiffuse.r * difdot);
            float g = color.g + (currentLight.color.g * mDiffuse.g * difdot);
            float b = color.b + (currentLight.color.b * mDiffuse.b * difdot);
            color = Color(r, g, b);
          }
        }

        // Check point lights
        for(int lightI = 0; lightI < pointLights.size(); lightI++) {
          PointLight currentLight = pointLights[lightI];
          Color lightColor = currentLight.color;
          Dir3D lightDir = (currentLight.pos - closePoint).normalized();
          Line3D lightLine = vee(closePoint, lightDir).normalized();

          // Check all spheres to see if light is blocked
          bool blocked = false;
          for(int sphereI = 0; sphereI < spheres.size(); sphereI ++) {
            Sphere currentSphere = spheres[sphereI];
            if (raySphereIntersect(closePoint, lightLine, currentSphere.pos, currentSphere.r) < INF) {
              blocked = true;
              break;
            }
          }
          
          // Add light to composite color if light isn't blocked
          if (!blocked) {
            Dir3D viewDir = (start - closePoint).normalized();
            Dir3D halfway = (lightDir + viewDir).normalized();
            float falloff = pow(closePoint.distTo(currentLight.pos),2);
      
            //cosine falloffs
            float difdot = (0 < dot(lightDir,normal)) ? dot(normal,lightDir) : 0;
            float specdot = (0 < pow(dot(halfway,normal),material.ns)) ? pow(dot(halfway,normal),material.ns) : 0;
      
            //point light r g b values
            float r = color.r + (mDiffuse.r * (lightColor.r / falloff) * difdot) + (mSpecular.r * (lightColor.r / falloff) * specdot);
            float g = color.g + (mDiffuse.g * (lightColor.g / falloff) * difdot) + (mSpecular.g * (lightColor.g / falloff) * specdot);
            float b = color.b + (mDiffuse.b * (lightColor.b / falloff) * difdot) + (mSpecular.b * (lightColor.b / falloff) * specdot);
            color = Color(r, g, b);
          }
        }

        float incidentAngle = acos(dot(vee(closePoint,(start-closePoint)).normalized(), vee(closePoint, normal).normalized()));

        // Calculate transparency
        if ((mTransmissive.r > 0 || mTransmissive.g > 0 || mTransmissive.b > 0) && recursionDepth < max_depth) {
          
          float refractionAngle = asin((1/material.ior) * sin(incidentAngle));
          Dir3D tDir = ((1/material.ior) * cos(incidentAngle) - cos(refractionAngle)) * normal - (1/material.ior) * (start - closePoint).normalized();
          Color tColor = rayCast(insidePoint, tDir, vee(insidePoint, tDir).normalized(), recursionDepth+1, false);
          //Color tColor = rayCast(insidePoint, tDir, vee(closePoint, tDir).normalized(), recursionDepth+1, false);
          //Color tColor = rayCast(insidePoint, rayDir, rayLine, recursionDepth+1, false);
          float r = color.r + mTransmissive.r * tColor.r;
          float g = color.g + mTransmissive.g * tColor.g;
          float b = color.b + mTransmissive.b * tColor.b;
          color = Color(r, g, b);
        }

        // Calculate reflection
        if ((mSpecular.r > 0 || mSpecular.g > 0 || mSpecular.b > 0) && recursionDepth < max_depth) {
          Dir3D reflectionDir = rayDir + dot(rayDir, normal)*normal*-2;
          Color rColor = rayCast(closePoint, reflectionDir, vee(closePoint, reflectionDir).normalized(), recursionDepth+1, false);
        }

      }
      
      //calculate color for a Triangle(non specified normal)
	  else if (closeTriangleIndex >= 0 && dif == 1) {
	    // apply coloring 
	   
	  }
	  
	  
	  //calculate color for a Triangle (specified normal)
	  else if (closeNormalTriangleIndex >= 0 && dif == 2) {
	    // apply coloring 
	   
	  }
      
      else {
        color = background;
      }

      return color;
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

  // Loop through each pixel
  for (int i = 0; i < img_width; i++){
    for (int j = 0; j < img_height; j++){
      //TODO: In what way does this assumes the basis is orthonormal?
      float u = (halfW - (imgW)*((i+0.5)/imgW));
      float v = (halfH - (imgH)*((j+0.5)/imgH));
      Point3D p = eye - d*forward + u*right + v*up;
      Dir3D rayDir = (p - eye); 
      Line3D rayLine = vee(eye,rayDir).normalized();  //Normalizing here is optional

      Color color = rayCast(eye, rayDir, rayLine, 0, false);
      
      outputImg.setPixel(i,j, color);
      //outputImg.setPixel(i,j, Color(fabs(i/imgW),fabs(j/imgH),fabs(0))); //TODO: Try this, what is it visualizing?
    }
  }
  auto t_end = std::chrono::high_resolution_clock::now();
  printf("Rendering took %.2f ms\n",std::chrono::duration<double, std::milli>(t_end-t_start).count());

  outputImg.write(imgName.c_str());
  return 0;
}
