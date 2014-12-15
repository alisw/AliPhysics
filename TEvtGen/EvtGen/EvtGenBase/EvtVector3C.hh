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
// Module: EvtGen/EvtVector3C.hh
//
// Description: Class for complex 3 vectors
//
// Modification history:
//
//    RYD     September 6, 1997         Module created
//
//------------------------------------------------------------------------

#ifndef EVTVECTOR3C_N
#define EVTVECTOR3C_N

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtVector3R.hh"
#include <iosfwd>

class EvtVector3C {

  friend EvtVector3C rotateEuler(const EvtVector3C& v,
				 double phi,double theta,double ksi);
  
  inline friend EvtVector3C operator*(const EvtComplex& c,const EvtVector3C& v2); 
  inline friend EvtVector3C operator*(const EvtComplex& c,const EvtVector3R& v2); 
  inline friend EvtComplex operator*(const EvtVector3R& v1,const EvtVector3C& v2); 
  inline friend EvtComplex operator*(const EvtVector3C& v1,const EvtVector3R& v2); 
  inline friend EvtComplex operator*(const EvtVector3C& v1,const EvtVector3C& v2); 
  inline friend EvtVector3C operator+(const EvtVector3C& v1,const EvtVector3C& v2);
  inline friend EvtVector3C operator-(const EvtVector3C& v1,const EvtVector3C& v2);
  inline friend EvtVector3C operator*(const EvtVector3C& v1,const EvtComplex& c);
  

public:

  EvtVector3C();
  EvtVector3C(const EvtComplex&,const EvtComplex&,const EvtComplex&);
  virtual ~EvtVector3C();
  inline void set(const int,const EvtComplex&);
  inline void set(const EvtComplex&,const EvtComplex&,const EvtComplex&);
  inline void set(double,double,double);
  inline EvtVector3C& operator*=(const EvtComplex& c);
  inline EvtVector3C& operator/=(const EvtComplex& c);
  inline EvtVector3C& operator+=(const EvtVector3C& v2);
  inline EvtVector3C& operator-=(const EvtVector3C& v2);
  inline EvtVector3C(const EvtVector3R& v1); 
  void applyRotateEuler(double phi,double theta,double ksi);
  inline const EvtComplex& get(int) const;
  inline EvtVector3C conj() const;
  EvtVector3C cross(const EvtVector3C& v2);
  friend std::ostream& operator<<(std::ostream& c,const EvtVector3C& v); 
  double dot( const EvtVector3C& p2 );  
private:

  EvtComplex v[3];
};

inline EvtVector3C::EvtVector3C(const EvtVector3R& v1){

  v[0]=EvtComplex(v1.get(0),0.0);
  v[1]=EvtComplex(v1.get(1),0.0);
  v[2]=EvtComplex(v1.get(2),0.0);

} 

inline void EvtVector3C::set(const int i,const EvtComplex& c){

  v[i]=c;
  
}

inline void EvtVector3C::set(const EvtComplex& x,const EvtComplex& y, 
			     const EvtComplex& z){

   v[0]=x; v[1]=y; v[2]=z; 
}

inline void EvtVector3C::set(double x,
			  double y,double z){

     v[0]=EvtComplex(x); v[1]=EvtComplex(y); v[2]=EvtComplex(z);
}

inline const EvtComplex& EvtVector3C::get(int i) const {

   return v[i];
}

inline EvtVector3C& EvtVector3C::operator*=(const EvtComplex& c){

  v[0]*=c;
  v[1]*=c;
  v[2]*=c;
  return *this;
}

inline EvtVector3C& EvtVector3C::operator/=(const EvtComplex& c){

  v[0]/=c;
  v[1]/=c;
  v[2]/=c;
  return *this;
}

inline EvtVector3C& EvtVector3C::operator+=(const EvtVector3C& v2){

  v[0]+=v2.v[0];
  v[1]+=v2.v[1];
  v[2]+=v2.v[2];
  return *this;
}

inline EvtVector3C& EvtVector3C::operator-=(const EvtVector3C& v2){

  v[0]-=v2.v[0];
  v[1]-=v2.v[1];
  v[2]-=v2.v[2];
  return *this;
}

inline EvtVector3C operator+(const EvtVector3C& v1,const EvtVector3C& v2) {

  return EvtVector3C(v1)+=v2;
}

inline EvtVector3C operator-(const EvtVector3C& v1,const EvtVector3C& v2) {

  return EvtVector3C(v1)-=v2;
}

inline EvtVector3C operator*(const EvtVector3C& v1,const EvtComplex& c) {

  return EvtVector3C(v1)*=c;
}

inline EvtVector3C operator*(const EvtComplex& c,const EvtVector3C& v2){

  return EvtVector3C(v2)*=c;
}

inline EvtVector3C operator*(const EvtComplex& c,const EvtVector3R& v2){

  return EvtVector3C(v2)*=c;
}

inline EvtComplex operator*(const EvtVector3R& v1,const EvtVector3C& v2){
 
  return v1.get(0)*v2.v[0]+v1.get(1)*v2.v[1]+v1.get(2)*v2.v[2];
}

inline EvtComplex operator*(const EvtVector3C& v1,const EvtVector3R& v2){
 
  return v1.v[0]*v2.get(0)+v1.v[1]*v2.get(1)+v1.v[2]*v2.get(2);
}

inline EvtComplex operator*(const EvtVector3C& v1,const EvtVector3C& v2){
 
  return v1.v[0]*v2.v[0]+v1.v[1]*v2.v[1]+v1.v[2]*v2.v[2];
}

inline EvtVector3C EvtVector3C::conj() const { 

  return EvtVector3C(::conj(v[0]),::conj(v[1]),
		  ::conj(v[2]));
}

#endif

