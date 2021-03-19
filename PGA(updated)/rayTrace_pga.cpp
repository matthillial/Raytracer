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

bool sameSide(Point3D p1, Point3D p2, Point3D a, Point3D b) {
  Dir3D cp1 = cross(b-a, p1-a);
  Dir3D cp2 = cross(b-a, p2-a);
  return dot(cp1, cp2) >= 0;
}

float rayPlaneIntersect(Point3D rayStart, Line3D rayLine, Point3D v1, Point3D v2, Point3D v3){
  Plane3D triPlane = vee(v1, v2, v3);	//make plane
  Point3D hitPoint = Point3D(wedge(rayLine,triPlane));
  // if ((dot(vee(v1,v2,hitPoint),triPlane) >= 0) && (dot(vee(v2,v3,hitPoint),triPlane) >= 0) && (dot(vee(v3,v1,hitPoint), triPlane) >= 0)) {
  //   float d = (hitPoint - rayStart).magnitude();
  //   return d;
  // }
  if (sameSide(hitPoint, v1, v2, v3) && sameSide(hitPoint, v2, v1, v3) && sameSide(hitPoint, v3, v1, v2) && dot(rayLine.dir(), hitPoint-rayStart) > 0) {
    float d = (hitPoint - rayStart).magnitude();
    return d;
  }
  return INF;
}

float areaTriangle(Point3D p1, Point3D p2, Point3D p3) {
  return 0.5 * vee(p1, p2, p3).magnitude();
}

vector<float> rayPlaneIntersectBary(Point3D rayStart, Line3D rayLine, Point3D v1, Point3D v2, Point3D v3) {
  v1 = v1.normalized();
  v2 = v2.normalized();
  v3 = v3.normalized();
  Plane3D triPlane = vee(v1, v2, v3);
  float area = 0.5 * triPlane.magnitude();
  Point3D hitPoint = Point3D(wedge(rayLine, triPlane));
  float b1, b2, b3;
  b1 = areaTriangle(hitPoint, v2, v3) / area;
  b2 = areaTriangle(hitPoint, v3, v1) / area;
  b3 = areaTriangle(hitPoint, v1, v2) / area;
  float d = (hitPoint - rayStart).magnitude();
  //float d0 = rayPlaneIntersect(rayStart, rayLine, v1, v2, v3);
  //if (Point3D(b1*v1 + b2*v2 + b3*v3) != hitPoint)
  // if (d0 != INF) {
  //   printf("MISMATCHED INTERSECTION: bary: %f, %f\t regular: %f\n", d, b1+b2+b3, d0);
  //   printf("\tv1: %f, %f, %f; v2: %f, %f, %f; v3: %f, %f, %f; hit: %f, %f, %f\n", v1.x, v1.y, v1.z, v2.x, v2.y, v2.z, v3.x, v3.y, v3.z, hitPoint.x, hitPoint.y, hitPoint.z);
  // }
  if ((b1 >= 0 && b1 <= 1) && (b2 >= 0 && b2 <= 1) && (b2 >= 0 && b2 <= 1) && (b1 + b2 + b3 - 1 < 0.00001) && (b1 + b2 + b3 - 1 > -0.00001) && dot(rayLine.dir(), hitPoint-rayStart) > 0) { //point in triangle
    
    return {d, b1, b2, b3};
    //vector<float> ret = {d, b1, b2, b3};
    //return ret;
    //return std::pair<float, vector<float>>
  } 
  return {INF, 0, 0, 0};
}

// check to see if any object is blocking
bool blocked(Point3D closePoint, Line3D lightLine, float lightDist) {
  // spheres
  for (int sphereI = 0; sphereI < spheres.size(); sphereI++) {
    Sphere currentSphere = spheres[sphereI];
    if (raySphereIntersect(closePoint, lightLine, currentSphere.pos, currentSphere.r) < lightDist) {
      return true;
    }
  }
  // triangles
  for (int triangleI = 0; triangleI < triangles.size(); triangleI++) {
    Triangle currentTriangle = triangles[triangleI];
    if (rayPlaneIntersect(closePoint, lightLine, currentTriangle.v1, currentTriangle.v2, currentTriangle.v3) < lightDist) {
      return true;
    }
  }
  // normal triangles
  for (int nTriangleI = 0; nTriangleI < normaltriangles.size(); nTriangleI++) {
    NormalTriangle currentNTriangle = normaltriangles[nTriangleI];
    if (rayPlaneIntersect(closePoint, lightLine, currentNTriangle.v1, currentNTriangle.v2, currentNTriangle.v3) < lightDist) {
      return true;
    }
  }
  return false;
}

Dir3D* refract(Dir3D in, Dir3D normal, float n, bool debug) {
  in = in.normalized();
  normal = normal.normalized();
  if (debug) printf("\tn: %f\n", n);
  float ndotu = dot(normal, in);
  float root = 1 - (n*n) + (ndotu*ndotu)*(n*n);
  if (root < 0) return NULL;
  return new Dir3D((( ( ((ndotu > 0) - (ndotu < 0)) * sqrt(root) ) - (ndotu*n) ) * normal) + (n*in));
}

Color light(Material material, Point3D closePoint, Dir3D normal, Point3D start, bool debug) {
  Color mAmbient = material.ambient;
  Color mDiffuse = material.diffuse;
  Color mSpecular = material.specular;
  Color mTransmissive = material.transmissive;
  // Set composite color to ambient light to start
  Color color = Color(mAmbient.r * ambient.r, mAmbient.g * ambient.g, mAmbient.b * ambient.b);

  // Check directional lights
  for (int lightI = 0; lightI < directionalLights.size(); lightI++)
  {
    DirLight currentLight = directionalLights[lightI];
    Dir3D lightDir = -1 * currentLight.dir;
    lightDir = lightDir.normalized();
    if (debug) printf("lightDir: %f %f %f\n", lightDir.x, lightDir.y, lightDir.z);
    Line3D lightLine = vee(closePoint, lightDir).normalized();

    // Add light into composite color if it isn't blocked
    if (!blocked(closePoint, lightLine, INF))
    {
      Dir3D viewDir = (start - closePoint).normalized();
      Dir3D halfway = (lightDir + viewDir).normalized();
      float difdot = (0 < dot(lightDir, normal)) ? dot(normal, lightDir) : 0;
      if (debug) {printf("difdot: %f, normal: %f %f %f\n",difdot, normal.x, normal.y, normal.z);}
      float specdot = (0 < pow(dot(halfway, normal), ns)) ? pow(dot(halfway, normal), ns) : 0;

      float r = color.r + (currentLight.color.r * mDiffuse.r * difdot);
      float g = color.g + (currentLight.color.g * mDiffuse.g * difdot);
      float b = color.b + (currentLight.color.b * mDiffuse.b * difdot);
      color = Color(r, g, b);
    }
  }

  // Check point lights
  for (int lightI = 0; lightI < pointLights.size(); lightI++)
  {
    PointLight currentLight = pointLights[lightI];
    Color lightColor = currentLight.color;
    Dir3D lightDir = (currentLight.pos - closePoint).normalized();
    float lightDist = (currentLight.pos - closePoint).magnitude();
    Line3D lightLine = vee(closePoint, lightDir).normalized();

    // Add light to composite color if light isn't blocked
    if (!blocked(closePoint, lightLine, lightDist))
    {
      Dir3D viewDir = (start - closePoint).normalized();
      Dir3D halfway = (lightDir + viewDir).normalized();
      float falloff = pow(closePoint.distTo(currentLight.pos), 2);

      //cosine falloffs
      float difdot = (0 < dot(lightDir, normal)) ? dot(normal, lightDir) : 0;
      float specdot = (0 < pow(dot(halfway, normal), material.ns)) ? pow(dot(halfway, normal), material.ns) : 0;

      //point light r g b values
      float r = color.r + (mDiffuse.r * (lightColor.r / falloff) * difdot) + (mSpecular.r * (lightColor.r / falloff) * specdot);
      float g = color.g + (mDiffuse.g * (lightColor.g / falloff) * difdot) + (mSpecular.g * (lightColor.g / falloff) * specdot);
      float b = color.b + (mDiffuse.b * (lightColor.b / falloff) * difdot) + (mSpecular.b * (lightColor.b / falloff) * specdot);
      color = Color(r, g, b);
    }
  }

  // Check spot lights
  for (int lightI = 0; lightI < spotLights.size(); lightI++) {
    SpotLight currentLight = spotLights[lightI];
    Color lightColor = currentLight.color;
    Dir3D lightDir = (currentLight.pos - closePoint).normalized();
    float lightDist = (currentLight.pos - closePoint).magnitude();
    Line3D lightLine = vee(closePoint, lightDir).normalized();

    //printf("lightDirx: %f   lightDiry: %f   lightDirz: %f \n", currentLight.dir.x, currentLight.dir.y, currentLight.dir.z);
    //printf("lightDirToHitx: %f   lightDirToHity: %f   lightDirToHitz: %f \n", lightDir.x, lightDir.y, lightDir.z);

    // Add light to composite color if light isn't blocked
    if (!blocked(closePoint, lightLine, lightDist)) {
      Dir3D viewDir = (start - closePoint).normalized();
      Dir3D halfway = (lightDir + viewDir).normalized();
      float falloff = pow(closePoint.distTo(currentLight.pos), 2);

      //determine falloff for light intensity using angles given
      float angle = acos(dot(currentLight.dir, (lightDir * -1)) / (currentLight.dir.magnitude() * lightDir.magnitude())) * (180 / M_PI);
      //printf("angle: %f   angle1: %f    angle2: %f \n\n", angle, currentLight.angle1, currentLight.angle2);
      if (angle < currentLight.angle1) { //angles less than angle1 should behave like point light
        falloff = 1 / pow(closePoint.distTo(currentLight.pos), 2);
        //printf("yeet");
      }
      else if (angle > currentLight.angle2) { //angles greating than angle2 should contribute nothing;
        falloff = 1 / INF;
        //printf("yeet");
      }
      else {
        falloff = 1 / pow(closePoint.distTo(currentLight.pos), 2) * (1 - ((angle - currentLight.angle1) / (currentLight.angle2 - currentLight.angle1)));
      }

      //cosine falloffs
      float difdot = (0 < dot(lightDir, normal)) ? dot(normal, lightDir) : 0;
      float specdot = (0 < pow(dot(halfway, normal), material.ns)) ? pow(dot(halfway, normal), material.ns) : 0;

      //point light r g b values
      float r = color.r + (mDiffuse.r * (lightColor.r * falloff) * difdot) + (mSpecular.r * (lightColor.r * falloff) * specdot);
      float g = color.g + (mDiffuse.g * (lightColor.g * falloff) * difdot) + (mSpecular.g * (lightColor.g * falloff) * specdot);
      float b = color.b + (mDiffuse.b * (lightColor.b * falloff) * difdot) + (mSpecular.b * (lightColor.b * falloff) * specdot);
      color = Color(r, g, b);
    }
  }
  return color;
}

Color rayCast(Point3D start, Dir3D rayDir, Line3D rayLine, int recursionDepth, float startIor, bool debug) {
  //variables for differentiating between sphere and triangle and triangle(normals specified)
  float dif = -1;

  // Calculate closest intersection for spheres
  float closeD = INF;  // Distance to nearest sphere
  int closeIndex = -1; // Index of nearest sphere
  for (int x = 0; x < spheres.size(); x++) {
    float d = raySphereIntersect(start, rayLine, spheres[x].pos, spheres[x].r);
    if (d < closeD) {
      closeD = d;
      closeIndex = x;
    }
  }

  // Calculate closest intersection for triangles(non specified normals)
  float closeDTriangle = INF;  // Distance to nearest triangle
  int closeTriangleIndex = -1; // Index of nearest triangle
  for (int x = 0; x < triangles.size(); x++) {
    float d = rayPlaneIntersect(start, rayLine, vertex[triangles[x].v1], vertex[triangles[x].v2], vertex[triangles[x].v3]);
    // vector<float> intersect = rayPlaneIntersectBary(start, rayLine, vertex[triangles[x].v1], vertex[triangles[x].v2], vertex[triangles[x].v3]);
    // float d = intersect[0];
    if (d < closeDTriangle) {
      closeDTriangle = d;
      closeTriangleIndex = x;
    }
  }

  // Calculate closest interesection for triangles(specified normals)
  float closeDNormalTriangle = INF;  // Distance to nearest triangle
  int closeNormalTriangleIndex = -1; // Index of nearest triangle
  float b1, b2, b3;
  for (int x = 0; x < normaltriangles.size(); x++) {
    float d0 = rayPlaneIntersect(start, rayLine, vertex[normaltriangles[x].v1], vertex[normaltriangles[x].v2], vertex[normaltriangles[x].v3]);
    vector<float> intersect = rayPlaneIntersectBary(start, rayLine, vertex[normaltriangles[x].v1], vertex[normaltriangles[x].v2], vertex[normaltriangles[x].v3]);
    float d = intersect[0];
    // if (d != d0) {
    //   printf("MISMATCHED INTERSECTION: bary: %f, %f\t regular: %f\n", d, intersect[0]+intersect[1]+intersect[2], d0);
    // }
    if (d < closeDTriangle) {
      closeDNormalTriangle = d;
      closeNormalTriangleIndex = x;
      b1 = intersect[1];
      b2 = intersect[2];
      b3 = intersect[3];
    }
  }

  // Calculate which intersection is the closest between all the primitives || dif = 0 for spheres || dif = 1 for triangle (normals not specified) || dif = 2 for traingles (normals specified)
  if (closeD < closeDTriangle && closeD < closeDNormalTriangle) dif = 0;
  else if (closeDTriangle < closeD && closeDTriangle < closeDNormalTriangle) dif = 1;
  else if (closeDNormalTriangle < closeD && closeDNormalTriangle < closeDTriangle) dif = 2;
  if (debug) printf("dif: %f\n", dif);

  

  // Calculate color of point for sphere
  Color color, mAmbient, mDiffuse, mSpecular, mTransmissive;
  Material material;
  Dir3D normal;
  Point3D closePoint, insidePoint;
  if (closeIndex >= 0 && dif == 0) { //sphere
    closePoint = start + (closeD - 0.001) * rayLine.dir();
    insidePoint = start + (closeD + 0.0001) * rayLine.dir();
    normal = (closePoint - spheres[closeIndex].pos).normalized();
    material = materials[spheres[closeIndex].matIndex];
  }
  //calculate color for a Triangle(non specified normal
  else if (closeTriangleIndex >= 0 && dif == 1) { 
    closePoint = start + (closeDTriangle - 0.001) * rayLine.dir();
    //printf("point: %f, %f, %f\n", closePoint.x, closePoint.y, closePoint.z);
    insidePoint = start + (closeDTriangle + 0.0001) * rayLine.dir();
    Triangle tri = triangles[closeTriangleIndex];
	  normal = cross(vertex[tri.v1]-vertex[tri.v2], vertex[tri.v3]-vertex[tri.v2]).normalized();
    if (dot(normal, rayDir) > 0) normal = -1*normal;
    //printf("normal: %f, %f, %f\n", normal.x, normal.y, normal.z);
    material = materials[tri.matIndex];
    //if (debug) printf("Material: %d\n", triangles[closeTriangleIndex].matIndex);
	}
	//calculate color for a Triangle (specified normal)
	else if (closeNormalTriangleIndex >= 0 && dif == 2) {
    closePoint = start + (closeDNormalTriangle - 0.001) * rayLine.dir();
    insidePoint = start + (closeDNormalTriangle + 0.0001) * rayLine.dir();
    NormalTriangle tri = normaltriangles[closeNormalTriangleIndex];
    normal = b1 * normals[tri.n1] + b2 * normals[tri.n2] + b3 * normals[tri.n3];
	  material = materials[tri.matIndex];
	}
  else return background; // return background color if no hit

  mAmbient = material.ambient;
  mDiffuse = material.diffuse;
  mSpecular = material.specular;
  mTransmissive = material.transmissive;

  color = light(material, closePoint, normal, start, debug);
  if (debug) printf("Color: %f, %f, %f\n", color.r, color.g, color.b);

    // Calculate transparency
    if ((mTransmissive.r > 0 || mTransmissive.g > 0 || mTransmissive.b > 0) && recursionDepth < max_depth) {
      rayDir = rayDir.normalized();

      // reflection
      Dir3D reflectionDir = rayDir + dot(rayDir, normal) * normal * -2;
      Color rColor = rayCast(closePoint, reflectionDir, vee(closePoint, reflectionDir).normalized(), recursionDepth + 1, -1, false);

      // refraction
      float n = 1 / material.ior;
      if (debug) printf("n: %f\n", n);
      Dir3D *tDir;
      float c;
      Color beerFalloff = Color(1, 1, 1);
      float insideDist = (closePoint - start).magnitude();
      bool tir = false; //total internal reflection

      if (dot(rayDir, normal) < 0) { // entering material
        tDir = refract(rayDir, normal, n, debug);
        c = dot(-1 * rayDir, normal);
        //beerFalloff = Color(1, 1, 1);
      }
      else { // leaving material
        //beerFalloff = something
        tDir = refract(rayDir, -1 * normal, 1 / n, debug);
        if (tDir == NULL) tir = true; // total internal reflection
        else c = dot(rayDir, normal);
      }
      float r0 = (n - 1) * (n - 1) / (n + 1) * (n + 1);
      float schlickR = (tir) ? 1 : r0 + (1 - r0) * pow((1 - c), 5);

      Color tColor = Color(0, 0, 0);
      if (!tir) tColor = rayCast(insidePoint, *tDir, vee(insidePoint, *tDir).normalized(), recursionDepth + 1, material.ior, debug);
      float r = color.r + beerFalloff.r * (schlickR * mTransmissive.r * rColor.r + (1 - schlickR) * mTransmissive.r * tColor.r);
      float g = color.g + beerFalloff.g * (schlickR * mTransmissive.g * rColor.g + (1 - schlickR) * mTransmissive.g * tColor.g);
      float b = color.b + beerFalloff.b * (schlickR * mTransmissive.b * rColor.b + (1 - schlickR) * mTransmissive.b * tColor.b);
      color = Color(r, g, b);
    }

    // Calculate reflection
    if ((mSpecular.r > 0 || mSpecular.g > 0 || mSpecular.b > 0) && recursionDepth < max_depth) {
      Dir3D reflectionDir = rayDir + dot(rayDir, normal) * normal * -2;
      Color rColor = rayCast(closePoint, reflectionDir, vee(closePoint, reflectionDir).normalized(), recursionDepth + 1, startIor, false);
      float r = color.r + mSpecular.r * rColor.r;
      float g = color.g + mSpecular.g * rColor.g;
      float b = color.b + mSpecular.b * rColor.b;
      color = Color(r, g, b);
    }
  if (debug) color = Color(0, 0, 0);
  return color;
}

int main(int argc, char** argv){
  //printf("\nArea test: %f\n\n", areaTriangle(Point3D(0, 0, 0), Point3D(2, 0, 0), Point3D(0, 1, 0)));
  //printf("Intersection test: %f\n", )
  
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
  #pragma omp parallel for num_threads(5)
  for (int i = 0; i < img_width; i++){
    for (int j = 0; j < img_height; j++){
      //TODO: In what way does this assumes the basis is orthonormal?
      float u = (halfW - (imgW)*((i+0.5)/imgW));
      float v = (halfH - (imgH)*((j+0.5)/imgH));
      Point3D p = eye - d*forward + u*right + v*up;
      Dir3D rayDir = (p - eye); 
      Line3D rayLine = vee(eye,rayDir).normalized();  //Normalizing here is optional

      Color color = rayCast(eye, rayDir, rayLine, 0, 1, (i==200 && j==440) ? true : false);
      
      outputImg.setPixel(i,j, color);
      //outputImg.setPixel(i,j, Color(fabs(i/imgW),fabs(j/imgH),fabs(0))); //TODO: Try this, what is it visualizing?
    }
  }
  auto t_end = std::chrono::high_resolution_clock::now();
  printf("Rendering took %.2f ms\n",std::chrono::duration<double, std::milli>(t_end-t_start).count());

  outputImg.write(imgName.c_str());
  return 0;
}
