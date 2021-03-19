//3D PGA Library
// Stephen J. Guy, sjguy@umn.edu

#ifndef PRINT_H
#define PRINT_H

#include "multivector.h"
#include "PGA_3D.h"

//Printing Operators:
//--------------
inline std::ostream& operator<<(std::ostream& os, const Point3D& mv){
  return os << std::string(mv);
}

inline std::ostream& operator<<(std::ostream& os, const Dir3D& mv){
  return os << std::string(mv);
}

inline std::ostream& operator<<(std::ostream& os, const HomogeneousPoint3D& mv){
  return os << std::string(mv);
}

inline std::ostream& operator<<(std::ostream& os, const IdealLine3D& mv){
  return os << std::string(mv);
}

inline std::ostream& operator<<(std::ostream& os, const Line3D& mv){
  return os << std::string(mv);
}

inline std::ostream& operator<<(std::ostream& os, const Plane3D& mv){
  return os << std::string(mv);
}

inline std::ostream& operator<<(std::ostream& os, const Motor3D& mv){
  return os << std::string(mv);
}

inline std::ostream& operator<<(std::ostream& os, const MultiVector& mv){
  return os << std::string(mv);
}

#endif