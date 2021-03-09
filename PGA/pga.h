//CSCI 5611 PGA Library
//This header defines the following 2D PGA operations:
// -wedge product, dot product, vee product, dual, reverse
//over the following elements:
// -points, directions, lines, motors, multivectors

#ifndef PGA_H
#define PGA_H

#include "multivector.h"
#include "primitives.h"

//--------
// Operations that act on various primitives
//--------

//Operations on Directions
inline Point2D Dir2D::operator+(Point2D rhs){
  return MultiVector(*this) + rhs;
}

inline Dir2D Dir2D::operator+(Dir2D rhs){
  return MultiVector(*this) + rhs;
}

inline Dir2D Dir2D::operator-(Dir2D rhs){
  return MultiVector(*this) - rhs;
}

//Operations on Points
inline Dir2D Point2D::operator-(Point2D rhs){
  return MultiVector(*this) - rhs;
}

inline Point2D Point2D::operator+(Dir2D rhs){
  return MultiVector(*this) + rhs;
}

inline Point2D Point2D::operator+(Point2D rhs){
  return MultiVector(*this) + rhs;
}


//Wedge Product 
//TODO:  WRONG: Could return a point or a direction!
inline Point2D wedge(Line2D lhs, Line2D rhs){
  return MultiVector(lhs).wedge(rhs);
}


//Vee/joint Operators 
inline Line2D vee(Point2D lhs, Point2D rhs){
  return MultiVector(lhs).vee(rhs);
}

//Quiz: Why no vee() that takes in 4 points?
inline float vee(Point2D a, Point2D b, Point2D c){
  return (MultiVector(a).vee(b).vee(c)).s;
}

inline Line2D vee(Point2D lhs, Dir2D rhs){
  return MultiVector(lhs).vee(rhs);
}

inline float vee(Line2D lhs, Point2D rhs){
  return (MultiVector(lhs).vee(rhs)).s;
}

inline float vee(Point2D lhs, Line2D rhs){
  return (MultiVector(lhs).vee(rhs)).s;
}

inline float vee(Line2D lhs, Dir2D rhs){
  return (MultiVector(lhs).vee(rhs)).s;
}

inline float vee(Dir2D lhs, Line2D rhs){
  return (MultiVector(lhs).vee(rhs)).s;
}

//Dot product
inline float dot(Line2D lhs, Line2D rhs){
  return (MultiVector(lhs).dot(rhs)).s;
}

inline Line2D dot(Point2D lhs, Line2D rhs){
  return MultiVector(lhs).dot(rhs);
}

inline Line2D dot(Line2D lhs, Point2D rhs){
  return MultiVector(lhs).dot(rhs);
}

//--------
// Constructors to turn MultiVectors in points, directions, lines, and motors
//--------

inline Point2D::Point2D(MultiVector mv){
  x = mv.yw/mv.xy;
  y = mv.wx/mv.xy;
}

inline Dir2D::Dir2D(MultiVector mv){
  x = mv.yw;
  y = mv.wx;
}

inline Motor2D::Motor2D(MultiVector mv){
  s = mv.s;
  yw = mv.yw;
  wx = mv.wx;
  xy = mv.xy;
}

inline Line2D::Line2D(MultiVector mv){
  x = mv.x;
  y = mv.y;
  w = mv.w;
}


//--------
// Helper Functions
//--------
float clamp(float f, float min, float max) {
  const float t = f < min ? min : f;
  return f > max ? max : f;
}

int sign(float f){
  if (f >= 0) return 1;
  return -1;
}

#endif