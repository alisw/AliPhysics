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
// Module: EvtGen/EvtVector3R.hh
//
// Description: Class to describe real 3 vectors
//
// Modification history:
//
//    RYD    Sept. 5, 1997      Module created
//
//------------------------------------------------------------------------

#ifndef EVTVECTOR3R_HH
#define EVTVECTOR3R_HH

#include <iosfwd>

class EvtVector3R {

  friend EvtVector3R rotateEuler(const EvtVector3R& v,
				 double phi,double theta,double ksi);

  inline friend EvtVector3R operator*(double c,const EvtVector3R& v2); 
  inline friend double operator*(const EvtVector3R& v1,const EvtVector3R& v2); 
  inline friend EvtVector3R operator+(const EvtVector3R& v1,const EvtVector3R& v2);
  inline friend EvtVector3R operator-(const EvtVector3R& v1,const EvtVector3R& v2);
  inline friend EvtVector3R operator*(const EvtVector3R& v1,double c);
  inline friend EvtVector3R operator/(const EvtVector3R& v1,double c);
  friend EvtVector3R cross(const EvtVector3R& v1,const EvtVector3R& v2);
  
public:
  EvtVector3R();
  EvtVector3R(double x,double y ,double z);
  virtual ~EvtVector3R(); 
  inline EvtVector3R& operator*=(const double c);
  inline EvtVector3R& operator/=(const double c);
  inline EvtVector3R& operator+=(const EvtVector3R& v2);
  inline EvtVector3R& operator-=(const EvtVector3R& v2);
  inline void set(int i,double d);
  inline void set(double x,double y ,double z);
  void applyRotateEuler(double phi,double theta,double ksi);
  inline double get(int i) const;
  friend std::ostream& operator<<(std::ostream& s,const EvtVector3R& v);
  double dot(const EvtVector3R& v2);
  double d3mag() const;

private:

  double v[3];

};

inline EvtVector3R& EvtVector3R::operator*=(const double c){

  v[0]*=c;
  v[1]*=c;
  v[2]*=c;
  return *this;
}

inline EvtVector3R& EvtVector3R::operator/=(const double c){

  v[0]/=c;
  v[1]/=c;
  v[2]/=c;
  return *this;
}

inline EvtVector3R& EvtVector3R::operator+=(const EvtVector3R& v2){

  v[0]+=v2.v[0];
  v[1]+=v2.v[1];
  v[2]+=v2.v[2];
  return *this;
}

inline EvtVector3R& EvtVector3R::operator-=(const EvtVector3R& v2){

  v[0]-=v2.v[0];
  v[1]-=v2.v[1];
  v[2]-=v2.v[2];
  return *this;
}

inline EvtVector3R operator*(double c,const EvtVector3R& v2){
  
  return EvtVector3R(v2)*=c;
}

inline EvtVector3R operator*(const EvtVector3R& v1,double c){
  
  return EvtVector3R(v1)*=c;
}

inline EvtVector3R operator/(const EvtVector3R& v1,double c){

  return EvtVector3R(v1)/=c; 
}

inline double operator*(const EvtVector3R& v1,const EvtVector3R& v2){

  return v1.v[0]*v2.v[0]+v1.v[1]*v2.v[1]+v1.v[2]*v2.v[2];
}

inline EvtVector3R operator+(const EvtVector3R& v1,const EvtVector3R& v2) {
  
  return EvtVector3R(v1)+=v2; 
}

inline EvtVector3R operator-(const EvtVector3R& v1,const EvtVector3R& v2) {
  
  return EvtVector3R(v1)-=v2; 

}

inline double EvtVector3R::get(int i) const {
  return v[i];
}

inline void EvtVector3R::set(int i,double d){
  
  v[i]=d;
}

inline void EvtVector3R::set(double x,double y, double z){

  v[0]=x;
  v[1]=y;
  v[2]=z;
}

#endif

