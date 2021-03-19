//CSCI 5607 3D PGA Library
//This header defines the following 3D PGA operations:
// -wedge product, dot product, vee product, dual, reverse

#ifndef MULTIVECTOR_H
#define MULTIVECTOR_H


#include <cmath>

#include <string>


inline std::string roundToString(float num, int p=3){
  std::string num_text = std::to_string(num+0.5*powf(10,-(p+1)));
  return num_text.substr(0, num_text.find(".")+(p+1));
}
  
//Multivector Struct
struct MultiVector{
  //Constructor
  MultiVector(float s=0, float w=0, float x=0, float y=0, float z=0,
              float wx=0, float wy=0, float wz=0, float xy=0, float zx=0, float yz=0, 
              float wyx=0, float wxz=0, float wzy=0, float xyz=0, float wxyz=0) : 
  s(s), w(w), x(x), y(y), z(z), wx(wx), wy(wy), wz(wz), xy(xy), zx(zx), yz(yz),
  wyx(wyx), wxz(wxz), wzy(wzy), xyz(xyz), wxyz(wxyz) {}

  float s; //Scalar
  float w, x,y,z; //Vector Components
  float wx, wy, wz, xy, zx, yz; //BiVector Components
  float wyx, wxz, wzy, xyz; //TriVector Components
  float wxyz; //Pseudoscalar

  MultiVector times(MultiVector rhs){
    float _s = s*rhs.s + x*rhs.x + y*rhs.y + z*rhs.z - xy*rhs.xy - zx*rhs.zx - yz*rhs.yz - xyz*rhs.xyz;
    float _w = s*rhs.w + w*rhs.s - x*rhs.wx - y*rhs.wy - z*rhs.wz + wx*rhs.x + wy*rhs.y + wz*rhs.z + xy*rhs.wyx +
               zx*rhs.wxz + yz*rhs.wzy + wyx*rhs.xy + wxz*rhs.zx + wzy*rhs.yz + xyz*rhs.wxyz - wxyz*rhs.wzy;
    float _x = s*rhs.x + x*rhs.s - y*rhs.xy + z*rhs.zx + xy*rhs.y - zx*rhs.z - yz*rhs.xyz - xyz*rhs.yz;
    float _y = s*rhs.y + x*rhs.xy + y*rhs.s - z*rhs.yz - xy*rhs.x - zx*rhs.xyz + yz*rhs.z - xyz*rhs.zx;
    float _z = s*rhs.z - x*rhs.zx + y*rhs.yz + z*rhs.s - xy*rhs.xyz + zx*rhs.x - yz*rhs.y - xyz*rhs.xy;
    float _wx = s*rhs.wx + w*rhs.x - x*rhs.w + wx*rhs.s - y*rhs.wyx + z*rhs.wxz - yz*rhs.wxyz - wyx*rhs.y + wxz*rhs.z - wxyz*rhs.yz;
    float _wy = s*rhs.wy + w*rhs.y - y*rhs.w + wy*rhs.s + x*rhs.wyx - z*rhs.wzy - zx*rhs.wxyz + wyx*rhs.x - wzy*rhs.z - wxyz*rhs.zx;
    float _wz = s*rhs.wz + w*rhs.z - z*rhs.w + wz*rhs.s - x*rhs.wxz + y*rhs.wzy - xy*rhs.wxyz - wxz*rhs.x + wzy*rhs.y - wxyz*rhs.xy;
    float _xy = s*rhs.xy + x*rhs.y - y*rhs.x + xy*rhs.s + z*rhs.xyz + zx*rhs.yz - yz*rhs.zx + xyz*rhs.z;
    float _zx = s*rhs.zx + z*rhs.x - x*rhs.z + zx*rhs.s + y*rhs.xyz - xy*rhs.yz + yz*rhs.xy + xyz*rhs.y;
    float _yz = s*rhs.yz + y*rhs.z - z*rhs.y + yz*rhs.s + x*rhs.xyz - zx*rhs.xy + xy*rhs.zx + + xyz*rhs.x;
    float _wyx = s*rhs.wyx - w*rhs.xy + x*rhs.wy - y*rhs.wx + z*rhs.wxyz - wx*rhs.y + wy*rhs.x - xy*rhs.w + zx*rhs.wzy - yz*rhs.wxz + wyx*rhs.s - wzy*rhs.zx + wxz*rhs.yz - wxyz*rhs.z;
    float _wxz = s*rhs.wxz - w*rhs.zx - x*rhs.wz + y*rhs.wxyz + z*rhs.wx + wx*rhs.z - wz*rhs.x - xy*rhs.wzy - zx*rhs.w + yz*rhs.wyx + wxz*rhs.s + wzy*rhs.xy - wyx*rhs.yz - wxyz*rhs.y;
    float _wzy = s*rhs.wzy - w*rhs.yz + x*rhs.wxyz + y*rhs.wz - z*rhs.wy - wy*rhs.z + wz*rhs.y + xy*rhs.wxz - zx*rhs.wyx - yz*rhs.w + wzy*rhs.s - wxz*rhs.xy + wyx*rhs.zx - wxyz*rhs.x;
    float _xyz = s*rhs.xyz + x*rhs.yz + y*rhs.zx + z*rhs.xy + xy*rhs.z + zx*rhs.y + yz*rhs.x + xyz*rhs.s;
    float _wxyz = s*rhs.wxyz + w*rhs.xyz + x*rhs.wzy + y*rhs.wxz + z*rhs.wyx + wx*rhs.yz + wy*rhs.zx + wz*rhs.xy +
                  xy*rhs.wz + zx*rhs.wy + yz*rhs.wx - wyx*rhs.z - wxz*rhs.y - wzy*rhs.x - xyz*rhs.w + wxyz*rhs.s;
    return MultiVector(_s, _w, _x, _y, _z, _wx, _wy, _wz, _xy, _zx, _yz, _wyx, _wxz, _wzy, _xyz, _wxyz);
  }

  MultiVector add(MultiVector rhs){
    float _s = s+rhs.s;
    float _w = w+rhs.w;
    float _x = x+rhs.x;
    float _y = y+rhs.y;
    float _z = z+rhs.z;
    float _wx = wx+rhs.wx;
    float _wy = wy+rhs.wy;
    float _wz = wz+rhs.wz;
    float _xy = xy+rhs.xy;
    float _zx = zx+rhs.zx;
    float _yz = yz+rhs.yz;
    float _wyx = wyx+rhs.wyx;
    float _wxz = wxz+rhs.wxz;
    float _wzy = wzy+rhs.wzy;
    float _xyz = xyz+rhs.xyz;
    float _wxyz = wxyz+rhs.wxyz;
    return MultiVector(_s, _w, _x, _y, _z, _wx, _wy, _wz, _xy, _zx, _yz, _wyx, _wxz, _wzy, _xyz, _wxyz);
  }

  MultiVector sub(MultiVector rhs){
    float _s = s-rhs.s;
    float _w = w-rhs.w;
    float _x = x-rhs.x;
    float _y = y-rhs.y;
    float _z = z-rhs.z;
    float _wx = wx-rhs.wx;
    float _wy = wy-rhs.wy;
    float _wz = wz-rhs.wz;
    float _xy = xy-rhs.xy;
    float _zx = zx-rhs.zx;
    float _yz = yz-rhs.yz;
    float _wyx = wyx-rhs.wyx;
    float _wxz = wxz-rhs.wxz;
    float _wzy = wzy-rhs.wzy;
    float _xyz = xyz-rhs.xyz;
    float _wxyz = wxyz-rhs.wxyz;
    return MultiVector(_s, _w, _x, _y, _z, _wx, _wy, _wz, _xy, _zx, _yz, _wyx, _wxz, _wzy, _xyz, _wxyz);
  }

  MultiVector mul(float f){
    return MultiVector(s*f, w*f, x*f, y*f, z*f, wx*f, wy*f, wz*f, xy*f, zx*f, yz*f, wyx*f, wxz*f, wzy*f, xyz*f, wxyz*f);
  }

  MultiVector div(float f){
    return MultiVector(s/f, w/f, x/f, y/f, z/f, wx/f, wy/f, wz/f, xy/f, zx/f, yz/f, wyx/f, wxz/f, wzy/f, xyz/f, wxyz/f);
  }

  MultiVector wedge(MultiVector rhs){
    float _s = s*rhs.s;
    float _w = s*rhs.w + w*rhs.s;
    float _x = s*rhs.x + x*rhs.s;
    float _y = s*rhs.y + y*rhs.s;
    float _z = s*rhs.z + z*rhs.s;
    float _wx = s*rhs.wx + w*rhs.x - x*rhs.w + wx*rhs.s;
    float _wy = s*rhs.wy + w*rhs.y - y*rhs.w + wy*rhs.s;
    float _wz = s*rhs.wz + w*rhs.z - z*rhs.w + wz*rhs.s;
    float _xy = s*rhs.xy + x*rhs.y - y*rhs.x + xy*rhs.s;
    float _zx = s*rhs.zx + z*rhs.x - x*rhs.z + zx*rhs.s;
    float _yz = s*rhs.yz + y*rhs.z - z*rhs.y + yz*rhs.s;
    float _wyx = s*rhs.wyx - w*rhs.xy + x*rhs.wy - y*rhs.wx - wx*rhs.y + wy*rhs.x - xy*rhs.w + wyx*rhs.s;
    float _wxz = s*rhs.wxz - w*rhs.zx - x*rhs.wz + z*rhs.wx + wx*rhs.z - wz*rhs.x - zx*rhs.w + wxz*rhs.s;
    float _wzy = s*rhs.wzy - w*rhs.yz + y*rhs.wz - z*rhs.wy - wy*rhs.z + wz*rhs.y - yz*rhs.w + wzy*rhs.s;
    float _xyz = s*rhs.xyz + x*rhs.yz + y*rhs.zx + z*rhs.xy + xy*rhs.z + zx*rhs.y + yz*rhs.x + xyz*rhs.s;
    float _wxyz = s*rhs.wxyz + w*rhs.xyz + x*rhs.wzy + y*rhs.wxz + z*rhs.wyx + wx*rhs.yz + wy*rhs.zx + wz*rhs.xy +
                  xy*rhs.wz + zx*rhs.wy + yz*rhs.wx - wyx*rhs.z - wxz*rhs.y - wzy*rhs.x - xyz*rhs.w + wxyz*rhs.s;
    return MultiVector(_s, _w, _x, _y, _z, _wx, _wy, _wz, _xy, _zx, _yz, _wyx, _wxz, _wzy, _xyz, _wxyz);
  }

  MultiVector dot(MultiVector rhs){ // dot(ei,ej) = 1 if i==j, 0 elsewise
    float _s = s*rhs.s + x*rhs.x + y*rhs.y + z*rhs.z + xy*rhs.xy + zx*rhs.zx + yz*rhs.yz + xyz*rhs.xyz;
    float _w = s*rhs.w + w*rhs.s - x*rhs.wx - y*rhs.wy - z*rhs.wz + wx*rhs.x + wy*rhs.y + wz*rhs.z + xy*rhs.wyx +
               zx*rhs.wxz + yz*rhs.wzy + wyx*rhs.xy + wxz*rhs.zx + wzy*rhs.yz + xyz*rhs.wxyz - wxyz*rhs.wzy;
    float _x = s*rhs.x + x*rhs.s - y*rhs.xy + z*rhs.zx + xy*rhs.y - zx*rhs.z - yz*rhs.xyz - xyz*rhs.yz;
    float _y = s*rhs.y + x*rhs.xy + y*rhs.s - z*rhs.yz - xy*rhs.x - zx*rhs.xyz + yz*rhs.z - xyz*rhs.zx;
    float _z = s*rhs.z - x*rhs.zx + y*rhs.yz + z*rhs.s - xy*rhs.xyz + zx*rhs.x - yz*rhs.y - xyz*rhs.xy;
    float _wx = s*rhs.wx - y*rhs.wyx + z*rhs.wxz + wx*rhs.s - yz*rhs.wxyz - wyx*rhs.y + wxz*rhs.z - wxyz*rhs.yz;
    float _wy = s*rhs.wy + x*rhs.wyx - z*rhs.wzy + wy*rhs.s - zx*rhs.wxyz + wyx*rhs.x - wzy*rhs.z - wxyz*rhs.zx;
    float _wz = s*rhs.wz - x*rhs.wxz + y*rhs.wzy + wz*rhs.s - xy*rhs.wxyz - wxz*rhs.x + wzy*rhs.y - wxyz*rhs.xy;
    float _xy = s*rhs.xy + z*rhs.xyz + xy*rhs.s + xyz*rhs.z;
    float _zx = s*rhs.zx + y*rhs.xyz + zx*rhs.s + xyz*rhs.y;
    float _yz = s*rhs.yz + x*rhs.xyz + yz*rhs.s + xyz*rhs.x;
    float _wyx = s*rhs.wyx + z*rhs.wxyz + wyx*rhs.s;
    float _wxz = s*rhs.wxz + y*rhs.wxyz + wxz*rhs.s;
    float _wzy = s*rhs.wzy + x*rhs.wxyz + wzy*rhs.s;
    float _xyz = s*rhs.xyz + xyz*rhs.s;
    float _wxyz = s*rhs.wxyz + rhs.wxyz*rhs.s; 
    return MultiVector(_s, _w, _x, _y, _z, _wx, _wy, _wz, _xy, _zx, _yz, _wyx, _wxz, _wzy, _xyz, _wxyz);
  }

  //dual(v): v*dual(v) = wxyz
  MultiVector dual(){
    return MultiVector(wxyz, xyz, wzy, wxz, wyx, yz, zx, xy, wz, wy, wx, z, y, x, w, s);
  }

  MultiVector vee(MultiVector rhs){
    return dual().wedge(rhs.dual()).dual();
  }

  //reverse(v): v*reverse(v) = |v|^2
  MultiVector reverse(){
    return MultiVector(s, w, x, y, z, -wx, -wy, -wz, -xy, -zx, -yz, -wyx, -wxz, -wzy, -xyz, wxyz);
  }

  float magnitude(){
    return std::sqrt(this->times(reverse()).s);
  }

  float magnitudeSqr(){
    return this->times(reverse()).s;
  }

  MultiVector normalized(){
    return this->times(1/magnitude());
  }

  MultiVector transform(MultiVector m){
    MultiVector self = *this;
    return m.times(self).times(m.reverse());
  }

  operator std::string() const {
    std::string retString = "";
    float p = 2; //precession
    if (fabs(s) >= 0.001) retString += roundToString(s,p) + " + ";
    if (fabs(x) >= 0.001) retString += roundToString(x,p) + "x +";
    if (fabs(y) >= 0.001) retString += roundToString(y,p) + "y +";
    if (fabs(z) >= 0.001) retString += roundToString(z,p) + "z +";
    if (fabs(w) >= 0.001) retString += roundToString(w,p) + "w +";
    if (fabs(wx) >= 0.001) retString += roundToString(wx,p) + "wx +";
    if (fabs(wy) >= 0.001) retString += roundToString(wy,p) + "wy +";
    if (fabs(wz) >= 0.001) retString += roundToString(wz,p) + "wz +";
    if (fabs(xy) >= 0.001) retString += roundToString(xy,p) + "xy +";
    if (fabs(zx) >= 0.001) retString += roundToString(zx,p) + "zx +";
    if (fabs(yz) >= 0.001) retString += roundToString(yz,p) + "yz +";
    if (fabs(wyx) >= 0.001) retString += roundToString(wyx,p) + "wyx +";
    if (fabs(wxz) >= 0.001) retString += roundToString(wxz,p) + "wxz +";
    if (fabs(wzy) >= 0.001) retString += roundToString(wzy,p) + "wzy +";
    if (fabs(xyz) >= 0.001) retString += roundToString(xyz,p) + "xyz +";
    if (fabs(wxyz) >= 0.001) retString += roundToString(wxyz,p) + "wxyz +";
    if (retString == "") retString = "0 +";
    retString = retString.substr(0,retString.length()-2);
    return retString;
  }

  void print(const char* title=""){
    printf("MutiVector:%s %s\n",title, std::string(*this).c_str());
  }

};


//MultiVector products
inline MultiVector wedge(MultiVector lhs, MultiVector rhs){
  return lhs.wedge(rhs);
}

inline MultiVector dot(MultiVector lhs, MultiVector rhs){
  return lhs.dot(rhs);
}

inline MultiVector vee(MultiVector lhs, MultiVector rhs){
  return lhs.vee(rhs);
}

//Various Operator Overloads
inline MultiVector operator*(MultiVector lhs, MultiVector rhs){
  return lhs.times(rhs);
}

inline MultiVector operator+(MultiVector lhs, MultiVector rhs){
  return lhs.add(rhs);
}

inline MultiVector operator-(MultiVector lhs, MultiVector rhs){
  return lhs.sub(rhs);
}

inline MultiVector operator/(MultiVector lhs, float f){
  return lhs.div(f);
}

inline MultiVector operator*(MultiVector lhs, float f){
  return lhs.mul(f);
}

inline MultiVector operator*(float f,MultiVector rhs){
  return rhs.mul(f);
}

//TODO: You may find these useful
// inline bool operator>(MultiVector lhs, float f){
//   return lhs.s > f;
// }

// inline bool operator>(float f,MultiVector rhs){
//   return f > rhs;
// }

// inline bool operator<(MultiVector lhs, float f){
//   return lhs.s < f;
// }

// inline bool operator<(float f,MultiVector rhs){
//   return f < rhs;
// }

//Motor
inline MultiVector transform(MultiVector lhs, MultiVector m){
  return m*lhs*m.reverse();
}



#endif