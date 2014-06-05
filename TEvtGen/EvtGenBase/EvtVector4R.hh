//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: EvtGen/EvtVector4R.hh
//
// Description: Class to describe real 4 vectors
//
// Modification history:
//
//    DJL/RYD     September 25, 1996         Module created
//
//------------------------------------------------------------------------

#ifndef EVTVECTOR4R_HH
#define EVTVECTOR4R_HH

#include <iostream>
#include <math.h>

class EvtVector3R;

class EvtVector4R {

  

  inline friend EvtVector4R operator*(double d,const EvtVector4R& v2); 
  inline friend EvtVector4R operator*(const EvtVector4R& v2,double d); 
  inline friend EvtVector4R operator/(const EvtVector4R& v2,double d); 
  inline friend double operator*(const EvtVector4R& v1,const EvtVector4R& v2); 
  inline friend EvtVector4R operator+(const EvtVector4R& v1,const EvtVector4R& v2); 
  inline friend EvtVector4R operator-(const EvtVector4R& v1,const EvtVector4R& v2); 
  
public:
  EvtVector4R();
  EvtVector4R(double e,double px,double py ,double pz);
  inline void set(int i,double d);
  inline void set(double e,double px,double py ,double pz);
  inline EvtVector4R& operator*=(double c);
  inline EvtVector4R& operator/=(double c);
  inline EvtVector4R& operator=(const EvtVector4R& v2);
  inline EvtVector4R& operator+=(const EvtVector4R& v2);
  inline EvtVector4R& operator-=(const EvtVector4R& v2);
  inline double get(int i) const;
  inline double cont(const EvtVector4R& v4) const;
  friend std::ostream& operator<<(std::ostream& s, const EvtVector4R& v);  
  double mass2() const;     
  double mass() const;
  void applyRotateEuler(double alpha,double beta,double gamma);
  void applyBoostTo(const EvtVector4R& p4, bool inverse = false);
  void applyBoostTo(const EvtVector3R& boost, bool inverse = false);
  EvtVector4R cross(const EvtVector4R& v2);
  double dot(const EvtVector4R& v2) const;
  double d3mag() const;

  // Added by AJB - calculate scalars in the rest frame of the current object
  double scalartripler3( const EvtVector4R& p1, const EvtVector4R& p2,
          const EvtVector4R& p3 ) const;
  double dotr3( const EvtVector4R& p1, const EvtVector4R& p2 ) const;
  double mag2r3( const EvtVector4R& p1 ) const;
  double magr3( const EvtVector4R& p1 ) const;


private:

  double v[4];

  inline double Square( double x ) const { return x*x; }

};

EvtVector4R rotateEuler(const EvtVector4R& rs,
				 double alpha,double beta,double gamma);
EvtVector4R boostTo(const EvtVector4R& rs,
		     const EvtVector4R& p4, bool inverse = false);
EvtVector4R boostTo(const EvtVector4R& rs,
		     const EvtVector3R& boost, bool inverse = false);

inline EvtVector4R& EvtVector4R::operator=(const EvtVector4R& v2){

  v[0]=v2.v[0];
  v[1]=v2.v[1];
  v[2]=v2.v[2];
  v[3]=v2.v[3];
  
  return *this; 
}

inline EvtVector4R& EvtVector4R::operator+=(const EvtVector4R& v2){

  v[0]+=v2.v[0];
  v[1]+=v2.v[1];
  v[2]+=v2.v[2];
  v[3]+=v2.v[3];
  
  return *this; 
}

inline EvtVector4R& EvtVector4R::operator-=(const EvtVector4R& v2){

  v[0]-=v2.v[0];
  v[1]-=v2.v[1];
  v[2]-=v2.v[2];
  v[3]-=v2.v[3];
  
  return *this; 
}

inline double EvtVector4R::mass2() const{

  return v[0]*v[0]-v[1]*v[1]-v[2]*v[2]-v[3]*v[3];
}

inline EvtVector4R operator*(double c,const EvtVector4R& v2){
  
  return EvtVector4R(v2)*=c;
}

inline EvtVector4R operator*(const EvtVector4R& v2,double c){
  
  return EvtVector4R(v2)*=c;
}

inline EvtVector4R operator/(const EvtVector4R& v2,double c){
  
  return EvtVector4R(v2)/=c;
}

inline EvtVector4R& EvtVector4R::operator*=(double c){

  v[0]*=c;  
  v[1]*=c;  
  v[2]*=c;  
  v[3]*=c;  

  return *this;
}

inline EvtVector4R& EvtVector4R::operator/=(double c){

  double cinv=1.0/c;  
  v[0]*=cinv;  
  v[1]*=cinv;  
  v[2]*=cinv;  
  v[3]*=cinv;  

  return *this;
}

inline double operator*(const EvtVector4R& v1,const EvtVector4R& v2){

  return v1.v[0]*v2.v[0]-v1.v[1]*v2.v[1]-
         v1.v[2]*v2.v[2]-v1.v[3]*v2.v[3];
}

inline double EvtVector4R::cont(const EvtVector4R& v4) const {
  
  return v[0]*v4.v[0]-v[1]*v4.v[1]-
         v[2]*v4.v[2]-v[3]*v4.v[3];
}

inline EvtVector4R operator-(const EvtVector4R& v1,const EvtVector4R& v2){
  
  return EvtVector4R(v1)-=v2;
}

inline EvtVector4R operator+(const EvtVector4R& v1,const EvtVector4R& v2){
  
  return EvtVector4R(v1)+=v2;
}

inline double EvtVector4R::get(int i) const {
  return v[i];
}

inline void EvtVector4R::set(int i,double d){
  
  v[i]=d;
}

inline void EvtVector4R::set(double e,double p1,double p2, double p3){

  v[0]=e;
  v[1]=p1;
  v[2]=p2;
  v[3]=p3;
}

#endif

