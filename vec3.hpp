// vec3.hpp
// A simple class to represent a vector in R^3.  No .cpp
// file because everything is inlined.
// 
// Author: Bob Sumner (sumnerb@inf.ethz.ch) 4 May 2006
// 
// Adapted from the "vl" vector library by Andrew Willmont.  See
// http://www.cs.cmu.edu/afs/cs/user/ajw/www/software/index.html#VL

#ifndef VEC3_HPP
#define VEC3_HPP

#include <math.h>

// ---------- Vec3 ------------------------------------------------------

class Vec3 {

public:

  // Constructors
  inline Vec3();                                 // Null vector
  inline Vec3(double x, double y, double z);     // [x y z]
  inline Vec3(const Vec3 &v);                    // copy constructor

  // Accessor functions
  inline int size() const { return(3); }
  inline double       &operator [] (int i);               
  inline const double &operator [] (int i) const;

  // Return pointer to data
  inline double *Ref() const;

  // Assignment operators
  inline Vec3   &operator =  (const Vec3 &a);                   
  inline Vec3   &operator += (const Vec3 &a);
  inline Vec3   &operator -= (const Vec3 &a);
  inline Vec3   &operator *= (const Vec3 &a);
  inline Vec3   &operator *= (double s);
  inline Vec3   &operator /= (const Vec3 &a);
  inline Vec3   &operator /= (double s);

  // Arithmetic operators
  inline Vec3   operator + (const Vec3 &a) const;       // v + a
  inline Vec3   operator - (const Vec3 &a) const;       // v - a
  inline Vec3   operator - () const;                    // -v
  inline Vec3   operator * (double s) const;            // v * s
  inline Vec3   operator / (double s) const;            // v / s

  // Initializers
  inline Vec3   &MakeZero();                            // Zero vector

  inline Vec3   &Normalize();                           // normalize vector

protected:
  
  double _data[3];
};

// ---------- Vec3 operators --------------------------------------------

inline double	dot(const Vec3 &a, const Vec3 &b);    // v . a
inline double	len(const Vec3 &v);                   // || v ||
inline double	sqrlen(const Vec3 &v);                // v . v
inline void	normalize(Vec3 &v);                     // v = norm(v)
inline Vec3	cross(const Vec3 &a, const Vec3 &b);    // a x b


// ---------- Vec3 inlines ----------------------------------------------

inline double &Vec3::operator [] (int i) {
  return(_data[i]);
}

inline const double &Vec3::operator [] (int i) const {
  return(_data[i]);
}

inline Vec3::Vec3() {
  _data[0] = 0.0;
  _data[1] = 0.0;
  _data[2] = 0.0;
}

inline Vec3::Vec3(double x, double y, double z) {
  _data[0] = x;
  _data[1] = y;
  _data[2] = z;
}

inline Vec3::Vec3(const Vec3 &v)  {
  _data[0] = v[0];
  _data[1] = v[1];
  _data[2] = v[2];
}

inline double *Vec3::Ref() const {
  return ((double *) _data);
}

inline Vec3 &Vec3::operator = (const Vec3 &v) {
  _data[0] = v[0];
  _data[1] = v[1];
  _data[2] = v[2];
  
  return(*this);
}

inline Vec3 &Vec3::operator += (const Vec3 &v) {
  _data[0] += v[0];
  _data[1] += v[1];
  _data[2] += v[2];
  
  return(*this);
}

inline Vec3 &Vec3::operator -= (const Vec3 &v) {
  _data[0] -= v[0];
  _data[1] -= v[1];
  _data[2] -= v[2];
  
  return(*this);
}

inline Vec3 &Vec3::operator *= (const Vec3 &a) {
  _data[0] *= a[0];
  _data[1] *= a[1];
  _data[2] *= a[2];
  
  return(*this);
}

inline Vec3 &Vec3::operator *= (double s) {
  _data[0] *= s;
  _data[1] *= s;
  _data[2] *= s;
	
  return(*this);
}

inline Vec3 &Vec3::operator /= (const Vec3 &a) {
  _data[0] /= a[0];
  _data[1] /= a[1];
  _data[2] /= a[2];
	
  return(*this);
}

inline Vec3 &Vec3::operator /= (double s) {
  _data[0] /= s;
  _data[1] /= s;
  _data[2] /= s;
	
  return(*this);
}

inline Vec3 Vec3::operator + (const Vec3 &a) const {
  Vec3 result;
  
  result[0] = _data[0] + a[0];
  result[1] = _data[1] + a[1];
  result[2] = _data[2] + a[2];
	
  return(result);
}

inline Vec3 Vec3::operator - (const Vec3 &a) const {
  Vec3 result;
	
  result[0] = _data[0] - a[0];
  result[1] = _data[1] - a[1];
  result[2] = _data[2] - a[2];
  
  return(result);
}

inline Vec3 Vec3::operator - () const {
  Vec3 result;
	
  result[0] = -_data[0];
  result[1] = -_data[1];
  result[2] = -_data[2];
	
  return(result);
}

inline Vec3 Vec3::operator * (double s) const {
  Vec3 result;
	
  result[0] = _data[0] * s;
  result[1] = _data[1] * s;
  result[2] = _data[2] * s;
	
  return(result);
}

inline Vec3 Vec3::operator / (double s) const {
  Vec3 result;
	
  result[0] = _data[0] / s;
  result[1] = _data[1] / s;
  result[2] = _data[2] / s;
	
  return(result);
}

inline Vec3 operator * (double s, const Vec3 &v) {
  return(v * s);
}

inline Vec3 &Vec3::MakeZero() {
  _data[0] = 0.0; _data[1] = 0.0; _data[2] = 0.0;
  return(*this);
}

inline Vec3 &Vec3::Normalize() {
  *this /= len(*this);
  return(*this);
}

inline double dot(const Vec3 &a, const Vec3 &b) {
  return(a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
}

inline double len(const Vec3 &v) {
  return(sqrt(dot(v, v)));
}

inline double sqrlen(const Vec3 &v) {
  return(dot(v, v));
}

inline void normalize(Vec3 &v)	{
  v /= len(v);
}

inline Vec3 cross(const Vec3 &a, const Vec3 &b) {
  Vec3 result;

  result[0] = a[1] * b[2] - a[2] * b[1];	
  result[1] = a[2] * b[0] - a[0] * b[2];	
  result[2] = a[0] * b[1] - a[1] * b[0];	

  return(result);
}

#endif
