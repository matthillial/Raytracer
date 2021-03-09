//3D PGA Library
// Stephen J. Guy, sjguy@umn.edu

#ifndef PGA_3D_H
#define PGA_3D_H

#include "multivector.h"

// This header-only library defines the follow primatives:

struct Point3D; //A 3D point (x,y,z) [homogeneous coordinate 1]
struct Dir3D; //A 3D displacement (x,y,z) [homogeneous coordinate 0]
struct HomogeneousPoint3D; //A homogeneous 3D point (x,y,z,w)
struct Line3D; //A line in 3D with Plucker-style coordinates P-Q = (x,y,z) and PxQ = (mx,my,mz)
struct IdealLine3D; //An "ideal" line at infinity representing only directionality (x,y,z)
struct Plane3D; //A plane with the equation: xX + yY + zZ + w = 0 [plane normal is (x,y,z), w indicates distance from orig.]
struct Motor3D; //An element which applies rotaton or tranlations, can be composed together [isomorphic to dual quaternions]

// ---------------------------------
//  Define Primatives
// ---------------------------------


struct Plane3D{
  float x,y,z,w;

  Plane3D(float x=0, float y=0, float z=0, float w=0): x(x), y(y), z(z), w(w)  {}

  Plane3D(MultiVector mv) : x(mv.x), y(mv.y), z(mv.z), w(mv.w) {}
  operator MultiVector(){ return MultiVector(0,w,x,y,z,0,0,0,0,0,0,0,0,0,0,0); }

  float magnitude(){
    return sqrt(x*x+y*y+z*z);
  }

  float magnitudeSqr(){
    return x*x+y*y+z*z;
  }

  Plane3D normalized(){
    float mag = magnitude();
    return Plane3D(x/mag, y/mag, z/mag, w/mag);
  }

  operator std::string() const {
    char result[150];
    sprintf(result, "Plane %.2fx+%.2fy+%.2fz+%f=0",x,y,z,w);
    return result;
  }

  void print(const char* title=""){
    printf("%s - %s\n",title, std::string(*this).c_str());
  }

};

//Lines
//Represents a line in plucker coordinates for a line runing through points P and Q
//  (x,y,z) is the lines "direction": P-Q
//  (mx,my,mz) defines the plane the line lies on as PxQ

struct IdealLine3D{
  float x,y,z;

  IdealLine3D(float x=0, float y=0, float z=0): x(x), y(y), z(z)  {}

  IdealLine3D(MultiVector mv) : x(mv.wx), y(mv.wy), z(mv.wz) {}
  operator MultiVector(){ return MultiVector(0,0,0,0,0,x,y,z,0,0,0,0,0,0,0,0); }

  float magnitude(){
    return sqrt(x*x+y*y+z*z);
  }

  float magnitudeSqr(){
    return x*x+y*y+z*z;
  }

  IdealLine3D normalized(){
    float mag = magnitude();
    return IdealLine3D(x/mag, y/mag, z/mag);
  }

  operator std::string() const {
    char result[100];
    sprintf(result, "Direction (%.2fx,%.2fy,%.2fz)",x,y,z);
    return result;
  }

  void print(const char* title=""){
    printf("%s - %s\n",title, std::string(*this).c_str());
  }

};

struct Line3D{
  float x,y,z,mx,my,mz;

  Line3D(float x=0, float y=0, float z=0, float mx=0, float my=0, float mz=0): x(x), y(y), z(z), mx(mx), my(my), mz(mz)  {}

  Line3D(MultiVector mv) : mx(mv.wx), my(mv.wy), mz(mv.wz), x(mv.yz), y(mv.zx), z(mv.xy) {}

  operator MultiVector(){ return MultiVector(0,0,0,0,0,mx,my,mz,z,y,x,0,0,0,0,0); }

  Dir3D dir();

  float magnitude(){
    return sqrt(x*x+y*y+z*z);
  }

  float magnitudeSqr(){
    return x*x+y*y+z*z;
  }

  Line3D normalized(){
    float mag = magnitude();
    return Line3D(x/mag, y/mag, z/mag, mx/mag, my/mag, mz/mag);
  }

  operator std::string() const {
    char result[200];
    sprintf(result, "Direction (%.2fx,%.2fy,%.2fz),  Moment: (%.2fx,%.2fy,%.2fz)", x, y, z, mx, my, mz);
    return result;
  }

  void print(const char* title=""){
    printf("%s - %s\n",title, std::string(*this).c_str());
  }

};

//Points and Displacements
//-------------------
struct Point3D{
  float x,y,z;

  Point3D(float x=0, float y=0, float z=0): x(x), y(y), z(z) {}

  Point3D scale(float f){
    return Point3D(x*f,y*f,z*f);
  }

  Point3D(MultiVector mv) : x(mv.wzy/mv.xyz), y(mv.wxz/mv.xyz), z(mv.wyx/mv.xyz) {};
  operator MultiVector(){return MultiVector(0,0,0,0,0,0,0,0,0,0,0,z,y,x,1,0);};

  Dir3D operator-(Point3D rhs);
  Point3D operator+(Dir3D rhs);
  Point3D operator-(Dir3D rhs);
  Point3D operator+(Point3D rhs);

  Point3D normalized(){
    return Point3D(x,y,z); //Points are already normalized
  }

  float distTo(Point3D p){
    return sqrt((x-p.x)*(x-p.x) + (y-p.y)*(y-p.y) + (z-p.z)*(z-p.z));
  }

  float distToSqr(Point3D p){
    return (x-p.x)*(x-p.x) + (y-p.y)*(y-p.y) + (z-p.z)*(z-p.z);
  }

  operator std::string() const {
    char result[100];
    sprintf(result, "Point: (%.2f, %.2f, %.2f)",x,y,z);
    return result;
  }

  void print(const char* title=""){
    printf("%s - %s\n",title, std::string(*this).c_str());
  }
  
};

struct Dir3D{
  float x,y,z;

  Dir3D(float x=0, float y=0, float z=0): x(x), y(y), z(z) {}

  Dir3D(MultiVector mv) : x(mv.wzy), y(mv.wxz), z(mv.wyx) {}
  operator MultiVector(){ return MultiVector(0,0,0,0,0,0,0,0,0,0,0,z,y,x,0,0); }

  Point3D operator+(Point3D rhs);
  Dir3D operator+(Dir3D rhs);
  Dir3D operator-(Dir3D rhs);

  float magnitude(){
    return sqrt(x*x+y*y+z*z);
  }

  float magnitudeSqr(){
    return x*x+y*y+z*z;
  }

  Dir3D normalized(){
    float mag = magnitude();
    return Dir3D(x/mag, y/mag, z/mag);
  }

  operator std::string() const {
    char result[100];
    sprintf(result, "Dir: (%.2f, %.2f, %.2f)",x,y,z);
    return result;
  }

  void print(const char* title=""){
    printf("%s - %s\n",title, std::string(*this).c_str());
  }
};

struct HomogeneousPoint3D{
  float x,y,z,w;

  HomogeneousPoint3D(float x=0, float y=0, float z=0, float w=1): x(x), y(y), z(z), w(w) {}

  HomogeneousPoint3D(MultiVector mv) : x(mv.wzy), y(mv.wxz), z(mv.wyx), w(mv.xyz) {}
  HomogeneousPoint3D(Point3D p) : x(p.x), y(p.y), z(p.z), w(1) {}
  HomogeneousPoint3D(Dir3D d) : x(d.x), y(d.y), z(d.z), w(0) {}
  operator MultiVector(){ return MultiVector(0,0,0,0,0,0,0,0,0,0,0,z,y,x,w,0); }

  HomogeneousPoint3D scale(float f){
    return HomogeneousPoint3D(x*f,y*f,z*f);
  }

  HomogeneousPoint3D operator-(HomogeneousPoint3D rhs);
  HomogeneousPoint3D operator+(HomogeneousPoint3D rhs);

  float magnitude(){ 
    return w;
  }

  float magnitudeSqr(){
    return w*w;
  }

  Point3D normalized(){
    assert(w > 0);
    return Point3D(x/w, y/w, z/w);
  }

  operator std::string() const {
    char result[100];
    sprintf(result, "H. Point: (%.2f, %.2f, %.2f, w=%.2f)",x,y,z,w);
    return result;
  }

  void print(const char* title=""){
    printf("%s - %s\n",title, std::string(*this).c_str());
  }
};

//Motors
//----------
//Equivalent to Dual Quaternions
struct Motor3D{
  float s; //Scalar
  float wx, wy, wz, yz, zx, xy; //BiVector Components
  float wxyz; //Pseudoscalar

  Motor3D(float s=1, float wx=0, float wy=0, float wz=0, float yz=0, float zx=0, float xy=0, float wxyz = 0): s(s), wx(wx), wy(wy), wz(wz), yz(yz), zx(zx), xy(xy), wxyz(wxyz) {}
  Motor3D(float s, Line3D l): s(s), wx(l.mx), wy(l.my), wz(l.mz), yz(l.x), zx(l.y), xy(l.z) {}
  Motor3D(MultiVector mv): s(mv.s), wx(mv.wx), wy(mv.wy), wz(mv.wz), yz(mv.yz), zx(mv.zx), xy(mv.xy), wxyz(mv.wxyz) {};
  operator MultiVector(){ return MultiVector(s,0,0,0,0,wx,wy,wz,xy,zx,yz,0,0,0,0,0); }

  Point3D operator()(Point3D p){ return MultiVector(*this).transform(p); }
  Dir3D operator()(Dir3D p){ return MultiVector(*this).transform(p); }
  HomogeneousPoint3D operator()(HomogeneousPoint3D p){ return MultiVector(*this).transform(p); }
  IdealLine3D operator()(IdealLine3D p){ return MultiVector(*this).transform(p); }
  Line3D operator()(Line3D p){ return MultiVector(*this).transform(p); }
  Plane3D operator()(Plane3D p){ return MultiVector(*this).transform(p); }


  operator std::string() const {
    char result[100];
    sprintf(result, "Motor: Scalar: %.2f, Bivector: %.2f %.2f %.2f  %.2f %.2f %.2f  PseudoScalar: %f \n", s, wx, wy, wz, yz, zx, xy, wxyz);
    return result;
  }

  void print(const char* title=""){
    printf("%s %s\n",title, std::string(*this).c_str());
  }
};

//Motor that translates perp. to perp_direction by dist
inline Motor3D Translator3D(float dist, IdealLine3D perp_direction) {
  return Motor3D(1,dist/2*MultiVector(perp_direction));
}

//Motor that rotates by angle around line
inline Motor3D Rotator3D(float angle, Line3D line) {
  float c = std::cos(angle/2.0);
  float s = std::sin(angle/2.0);
  return Motor3D(c,s*line);
}

//Printing Operators:
//--------------
std::ostream& operator<<(std::ostream& os, const Point3D& mv){
  return os << std::string(mv);
}

std::ostream& operator<<(std::ostream& os, const Dir3D& mv){
  return os << std::string(mv);
}

std::ostream& operator<<(std::ostream& os, const HomogeneousPoint3D& mv){
  return os << std::string(mv);
}

std::ostream& operator<<(std::ostream& os, const IdealLine3D& mv){
  return os << std::string(mv);
}

std::ostream& operator<<(std::ostream& os, const Line3D& mv){
  return os << std::string(mv);
}

std::ostream& operator<<(std::ostream& os, const Plane3D& mv){
  return os << std::string(mv);
}

std::ostream& operator<<(std::ostream& os, const Motor3D& mv){
  return os << std::string(mv);
}

//--------
// Operations that act on various primitives
//--------


//Operations on Points:
inline Point3D Point3D::operator+(Point3D rhs){
  return Point3D(x+rhs.x,y+rhs.y,z+rhs.z);
}

inline Dir3D Point3D::operator-(Point3D rhs){
  return Dir3D(x-rhs.x,y-rhs.y,z-rhs.z);
}

inline Point3D Point3D::operator+(Dir3D rhs){
  return Point3D(x+rhs.x,y+rhs.y,z+rhs.z);
}

inline Point3D Point3D::operator-(Dir3D rhs){
  return Point3D(x-rhs.x,y-rhs.y,z-rhs.z);
}



//Operations on Directions:
inline Point3D Dir3D::operator+(Point3D rhs){
  return Point3D(x+rhs.x,y+rhs.y,z+rhs.z);
}

inline Dir3D Dir3D::operator+(Dir3D rhs){
  return Dir3D(x+rhs.x,y+rhs.y,z+rhs.z);
}

inline Dir3D Dir3D::operator-(Dir3D rhs){
  return Dir3D(x-rhs.x,y-rhs.y,z-rhs.z);
}

inline Dir3D operator*(Dir3D d, float f){
  return Dir3D(d.x*f,d.y*f,d.z*f);
}

inline Dir3D operator*(float f,Dir3D d){
  return Dir3D(d.x*f,d.y*f,d.z*f);
}

//Operations on Homogeneous Points
inline HomogeneousPoint3D HomogeneousPoint3D::operator-(HomogeneousPoint3D rhs){
  return MultiVector(*this) - rhs;
}

inline HomogeneousPoint3D HomogeneousPoint3D::operator+(HomogeneousPoint3D rhs){
  return MultiVector(*this) + rhs;
}

//Operations on Lines
Dir3D Line3D::dir(){
  return Dir3D(x,y,z);
}

//Cross Product (Only for Dir3D)
inline Dir3D cross(Dir3D a, Dir3D b){
  return Dir3D(a.y*b.z - b.y*a.z, a.z*b.x - a.x*b.z, a.x*b.y - b.x*a.y);
}

//Wedge/Meet Operations
inline Line3D wedge(Point3D lhs, Point3D rhs){
  return MultiVector(lhs).wedge(rhs);
}

inline Plane3D wedge(Point3D a, Point3D b, Point3D c){
  return MultiVector(a).vee(b).vee(c);
}

inline HomogeneousPoint3D wedge(Line3D l, Plane3D p){
  //return MultiVector(l).wedge(p);
  float wyx = - l.mx*p.y + l.my*p.x - l.z*p.w;
  float wxz = + l.mx*p.z - l.mz*p.x - l.y*p.w;
  float wzy = - l.my*p.z + l.mz*p.y - l.x*p.w;
  float xyz = + l.z*p.z + l.y*p.y + l.x*p.x;
  return HomogeneousPoint3D(wzy,wxz,wyx,xyz);
}

inline HomogeneousPoint3D wedge(Plane3D lhs, Line3D rhs){
  return MultiVector(lhs).wedge(rhs);
}


// Times

inline MultiVector operator*(Plane3D p, Line3D l){
    float w = - p.x*l.mx - p.y*l.my - p.z*l.mz;
    float x = - p.y*l.z + p.z*l.y;
    float y = p.x*l.z - p.z*l.x;
    float z = - p.x*l.y + p.y*l.x;
    float wyx = -p.w*l.z + p.x*l.my - p.y*l.mx;
    float wxz = -p.w*l.y - p.x*l.mz + p.z*l.mx;
    float wzy = -p.w*l.x + p.y*l.mz - p.z*l.my;
    float xyz = p.x*l.x + p.y*l.y + p.z*l.z;
    return MultiVector(0, w, x, y, z, 0, 0, 0, 0, 0, 0, wyx, wxz, wzy, xyz, 0);
}


//Compute -plp
inline MultiVector sandwhich(Plane3D p, Line3D l){
  //return -1*Multivector(p)*MultiVVector(l)*MultiVVector(p)
  
  //Compute p*l
  float w = - p.x*l.mx - p.y*l.my - p.z*l.mz;
  float x = - p.y*l.z + p.z*l.y;
  float y = p.x*l.z - p.z*l.x;
  float z = - p.x*l.y + p.y*l.x;
  float wyx = -p.w*l.z + p.x*l.my - p.y*l.mx;
  float wxz = -p.w*l.y - p.x*l.mz + p.z*l.mx;
  float wzy = -p.w*l.x + p.y*l.mz - p.z*l.my;
  float xyz = p.x*l.x + p.y*l.y + p.z*l.z;
  
  //Compute (p*l)*p
  float wx = w*p.x - x*p.w - wyx*p.y + wxz*p.z;
  float wy = w*p.y - y*p.w + wyx*p.x - wzy*p.z;
  float wz = w*p.z - z*p.w - wxz*p.x + wzy*p.y;
  float xy = x*p.y - y*p.x + xyz*p.z;
  float zx = z*p.x - x*p.z + xyz*p.y;
  float yz = y*p.z - z*p.y + xyz*p.x;
  return MultiVector(0, 0, 0,0,0, -wx, -wy, -wz, -xy, -zx, -yz, 0, 0, 0 ,0, 0);
}

// Vee/Join Operations

inline Line3D vee(Point3D lhs, Point3D rhs){
  //return MultiVector(lhs).vee(rhs);
  float wx = rhs.x - lhs.x;
  float wy = rhs.y - lhs.y;
  float wz = rhs.z - lhs.z;
  float xy = lhs.x*rhs.y - lhs.y*rhs.x;
  float zx = lhs.z*rhs.x - lhs.x*rhs.z;
  float yz = lhs.y*rhs.z - lhs.z*rhs.y;
  return Line3D(wx, wy, wz, yz, zx, xy); //swapped for dual
}

inline Line3D vee(Point3D p, Dir3D d){
  //return MultiVector(p).vee(d);
  float wx = d.x;
  float wy = d.y;
  float wz = d.z;
  float xy = p.x*d.y - p.y*d.x;
  float zx = p.z*d.x - p.x*d.z;
  float yz = p.y*d.z - p.z*d.y;
  return Line3D(wx,wy,wz,yz,zx,xy); //swapped b/c vee is dual
}

inline Plane3D vee(Point3D a, Point3D b, Point3D c){
  //return MultiVector(a).vee(b).vee(c);
  //Vee to points a and b to make a line
  float wx = b.x - a.x;
  float wy = b.y - a.y;
  float wz = b.z - a.z;
  float xy = a.x*b.y - a.y*b.x;
  float zx = a.z*b.x - a.x*b.z;
  float yz = a.y*b.z - a.z*b.y;

  //Vee line and point c to get a plane
  float wyx = - wx*c.y + wy*c.x - xy;
  float wxz = + wx*c.z - wz*c.x - zx;
  float wzy = - wy*c.z + wz*c.y - yz;
  float xyz = + xy*c.z + zx*c.y + yz;

  return Plane3D(wzy,wxz,wyx,xyz); //Return dual as it's a vee operation
}

inline Plane3D vee(Point3D a, Point3D b, HomogeneousPoint3D c){
  //return MultiVector(a).vee(b).vee(c);
  float wx = b.x - a.x;
  float wy = b.y - a.y;
  float wz = b.z - a.z;
  float xy = a.x*b.y - a.y*b.x;
  float zx = a.z*b.x - a.x*b.z;
  float yz = a.y*b.z - a.z*b.y;

  //Vee line and point c to get a plane
  float wyx = - wx*c.y + wy*c.x - xy*c.w;
  float wxz = + wx*c.z - wz*c.x - zx*c.w;
  float wzy = - wy*c.z + wz*c.y - yz*c.w;
  float xyz = + xy*c.z + zx*c.y + yz*c.x;

  return Plane3D(wzy,wxz,wyx,xyz); //Return dual as it's a vee operation
}

inline Plane3D vee(Line3D l, Point3D p){
  return MultiVector(l).vee(p);
}


//Dot product
inline float dot(Plane3D a, Plane3D b){
  //return (MultiVector(a).dot(b)).s;
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

inline float dot(Line3D a, Line3D b){
  //return (MultiVector(lhs).dot(rhs)).s;
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

inline float dot(Dir3D a, Dir3D b){
  //return MultiVector(a).dot(b);
  return a.x*b.x + a.y*b.y + a.z*b.z;
}

inline Plane3D dot(Point3D p, Line3D l){
  //return MultiVector(p).dot(l);
  float w = l.z*p.z + l.y*p.y + l.x*p.x;
  float x = -l.x;
  float y = -l.y;
  float z = -l.z;
  return Plane3D(x,y,z,w);
}

inline Plane3D dot(Line3D l, Point3D p){
  //return MultiVector(lhs).dot(rhs);
  float w = l.z*p.z + l.y*p.y + l.x*p.x;
  float x = -l.x;
  float y = -l.y;
  float z = -l.z;
  return Plane3D(x,y,z,w);
}

inline Line3D dot(Point3D lhs, Plane3D rhs){
  return MultiVector(lhs).dot(rhs);
}

inline Line3D dot(Plane3D lhs, Point3D rhs){
  return MultiVector(lhs).dot(rhs);
}

inline Plane3D dot(Line3D lhs, Plane3D rhs){
  return MultiVector(lhs).dot(rhs);
}

inline Plane3D dot(Plane3D lhs, Line3D rhs){
  return MultiVector(lhs).dot(rhs);
}


//--------------------------

inline Point3D proj(Point3D p, Line3D l){
  //return dot(l,p)*l;
  
  //Dot point and line to form a plane
  float pw = l.z*p.z + l.y*p.y + l.x*p.x;
  float px = -l.x;
  float py = -l.y;
  float pz = -l.z;

  //Times plane and line to get point (no vector component b/c plane and line are orthognol)
  float wyx = -pw*l.z + px*l.my - py*l.mx;
  float wxz = -pw*l.y - px*l.mz + pz*l.mx;
  float wzy = -pw*l.x + py*l.mz - pz*l.my;
  float xyz = px*l.x + py*l.y + pz*l.z;
  return Point3D(wzy/xyz,wxz/xyz,wyx/xyz);
}

#endif