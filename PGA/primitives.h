//CSCI 5607 PGA Library
//This header defines the following 2D PGA elements:
// -points, directions, lines, motors, multivectors

#ifndef PRIMITIVES_H
#define PRIMITIVES_H

struct MultiVector;
struct Dir2D;

//Lines
struct Line2D{
  float x,y,w;

  Line2D(float x=0, float y=0, float w=0): x(x), y(y), w(w)  {}

  Line2D(MultiVector mv);

  float magnitude(){
    return sqrt(x*x+y*y);
  }

  float magnitudeSqr(){
    return x*x+y*y;
  }

  Line2D normalized(){
    float mag = magnitude();
    return Line2D(x/mag, y/mag, w/mag);
  }

  operator std::string() const {
    char result[100];
    sprintf(result, "%.2fx+%.2fy+%.2f=0",x,y,w);
    return result;
  }

  void print(const char* title=""){
    printf("Line:%s %s\n",title, std::string(*this).c_str());
  }

};

//Points and Directions
//-------------------
struct Point2D{
  float x,y;

  Point2D(float x=0, float y=0): x(x), y(y)  {}

  Point2D scale(float f){
    return Point2D(x*f,y*f);
  }

  Point2D(MultiVector mv);

  Dir2D operator-(Point2D rhs);
  Point2D operator+(Dir2D rhs);
  Point2D operator+(Point2D rhs);

  Point2D normalized(){
    return Point2D(x,y); //Points are already normalized
  }

  operator std::string() const {
    char result[100];
    sprintf(result, "(%.2f,%.2f)",x,y);
    return result;
  }

  void print(const char* title=""){
    printf("Point:%s %s\n",title, std::string(*this).c_str());
  }
  
};

struct Dir2D{
  float x,y;

  Dir2D(float x=0, float y=0): x(x), y(y)  {}

  Dir2D(MultiVector mv);

  Point2D operator+(Point2D rhs);
  Dir2D operator+(Dir2D rhs);
  Dir2D operator-(Dir2D rhs);

  float magnitude(){
    return sqrt(x*x+y*y);
  }

  float magnitudeSqr(){
    return x*x+y*y;
  }

  Dir2D normalized(){
    float mag = magnitude();
    return Dir2D(x/mag, y/mag);
  }

  Dir2D perp(){
    return Dir2D(y, -x);
  }

  operator std::string() const {
    char result[100];
    sprintf(result, "(%.2f,%.2f)",x,y);
    return result;
  }

  void print(const char* title=""){
    printf("Dir: %s %s\n",title, std::string(*this).c_str());
  }
};
//-------------------


struct Motor2D{
  float s; //Scalar
  float yw, wx, xy; //BiVector Components

  Motor2D(float s=1, float yw=0, float wx=0, float xy=0): s(s), yw(yw), wx(wx), xy(xy) {}
  Motor2D(MultiVector mv);

  void print(const char* title=""){
    printf("Motor: %s, Scalar: %.2f, Bivector: %.2f %.2f %.2f\n",title, s, yw, wx, xy);
  }
};

Motor2D Translator2D(float dist, Dir2D perp_direction) {
  return Motor2D(1,dist/2.0*perp_direction.x,dist/2.0*perp_direction.y,0);
}

Motor2D Translator2D(Dir2D displacement) {
  float dist = displacement.magnitude();
  displacement = displacement.normalized();
  if (dist > 0){
    return Translator2D(dist,displacement.perp()); //TODO: Add perp to slides somewhere
  }
  return Motor2D(); //Identity transformation
}

Motor2D Rotator2D(float angle, Point2D point) {
  float c = std::cos(angle/2.0);
  float s = std::sin(angle/2.0);
  return Motor2D(c,s*point.x,s*point.y,s);
}

#endif
