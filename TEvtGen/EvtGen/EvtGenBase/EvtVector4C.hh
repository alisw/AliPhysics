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
// Module: EvtGen/EvtVector4C.hh
//
// Description: Class for complex 4 vectors
//
// Modification history:
//
//    DJL/RYD     September 25, 1996         Module created
//
//------------------------------------------------------------------------

#ifndef EVTVECTOR4C_HH
#define EVTVECTOR4C_HH

#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtVector3C.hh"
#include "EvtGenBase/EvtVector4R.hh"

#include <iosfwd>

class EvtVector4C {

  friend EvtVector4C rotateEuler(const EvtVector4C& e,
				 double alpha,double beta,double gamma);
  friend EvtVector4C boostTo(const EvtVector4C& e,
			     const EvtVector4R p4);
  friend EvtVector4C boostTo(const EvtVector4C& e,
			     const EvtVector3R boost);  
  inline friend EvtVector4C operator*(double d,const EvtVector4C& v2); 
  inline friend EvtVector4C operator*(const EvtComplex& c,const EvtVector4C& v2); 
  inline friend EvtVector4C operator*(const EvtVector4C& v2,const EvtComplex& c); 
  inline friend EvtVector4C operator*(const EvtComplex& c,const EvtVector4R& v2); 
  inline friend EvtComplex operator*(const EvtVector4R& v1,const EvtVector4C& v2); 
  inline friend EvtComplex operator*(const EvtVector4C& v1,const EvtVector4R& v2); 
  inline friend EvtComplex operator*(const EvtVector4C& v1,const EvtVector4C& v2); 
  friend EvtVector4C operator+(const EvtVector4C& v1,const EvtVector4C& v2);
  friend EvtVector4C operator-(const EvtVector4C& v1,const EvtVector4C& v2);
  
public:

  EvtVector4C();
  EvtVector4C(const EvtComplex&,const EvtComplex&,
	      const EvtComplex&,const EvtComplex&);
  virtual ~EvtVector4C();
  inline void set(int,const EvtComplex&);
  inline void set(const EvtComplex&,const EvtComplex&,
		  const EvtComplex&,const EvtComplex&);
  inline void set(double,double,double,double);
  inline EvtVector4C(const EvtVector4R& v1); 
  inline const EvtComplex& get(int) const;
  inline EvtComplex cont(const EvtVector4C& v4) const;
  inline EvtVector4C conj() const;
  EvtVector3C vec() const;
  inline EvtVector4C& operator=(const EvtVector4C& v2);
  inline EvtVector4C& operator-=(const EvtVector4C& v2);
  inline EvtVector4C& operator+=(const EvtVector4C& v2);
  inline EvtVector4C& operator*=(const EvtComplex& c);
  void applyRotateEuler(double alpha,double beta,double gamma);
  void applyBoostTo(const EvtVector4R& p4);
  void applyBoostTo(const EvtVector3R& boost);
  friend std::ostream& operator<<(std::ostream& s, const EvtVector4C& v);
  double dot( const EvtVector4C& p2 );  
private:

  EvtComplex v[4];

};

inline EvtVector4C& EvtVector4C::operator=(const EvtVector4C& v2){

  v[0]=v2.v[0];
  v[1]=v2.v[1];
  v[2]=v2.v[2];
  v[3]=v2.v[3];

  return *this;
}

inline EvtVector4C& EvtVector4C::operator+=(const EvtVector4C& v2){

  v[0]+=v2.v[0];
  v[1]+=v2.v[1];
  v[2]+=v2.v[2];
  v[3]+=v2.v[3];

  return *this;
}

inline EvtVector4C& EvtVector4C::operator-=(const EvtVector4C& v2){

  v[0]-=v2.v[0];
  v[1]-=v2.v[1];
  v[2]-=v2.v[2];
  v[3]-=v2.v[3];

  return *this;
}

inline void EvtVector4C::set(int i,const EvtComplex& c){

  v[i]=c;
}

inline EvtVector3C EvtVector4C::vec() const {

  return EvtVector3C(v[1],v[2],v[3]);
}

inline void EvtVector4C::set(const EvtComplex& e,const EvtComplex& p1,
			     const EvtComplex& p2,const EvtComplex& p3){

   v[0]=e; v[1]=p1; v[2]=p2; v[3]=p3;
}

inline void EvtVector4C::set(double e,double p1,
			  double p2,double p3){

   v[0]=EvtComplex(e); v[1]=EvtComplex(p1); v[2]=EvtComplex(p2); v[3]=EvtComplex(p3);
}

inline const EvtComplex& EvtVector4C::get(int i) const {

   return v[i];
}

inline EvtVector4C operator+(const EvtVector4C& v1,const EvtVector4C& v2) {

  return EvtVector4C(v1)+=v2;
}

inline EvtVector4C operator-(const EvtVector4C& v1,const EvtVector4C& v2) {

  return EvtVector4C(v1)-=v2;
}

inline EvtComplex EvtVector4C::cont(const EvtVector4C& v4) const {

  return v[0]*v4.v[0]-v[1]*v4.v[1]-
         v[2]*v4.v[2]-v[3]*v4.v[3];
}

inline EvtVector4C& EvtVector4C::operator*=(const EvtComplex& c) {

  v[0]*=c;
  v[1]*=c;
  v[2]*=c;
  v[3]*=c;

  return *this;
}

inline EvtVector4C operator*(double d,const EvtVector4C& v2){

  return EvtVector4C(v2.v[0]*d,v2.v[1]*d,v2.v[2]*d,v2.v[3]*d);
}

inline EvtVector4C operator*(const EvtComplex& c,const EvtVector4C& v2){

  return EvtVector4C(v2)*=c;
}

inline EvtVector4C operator*(const EvtVector4C& v2,const EvtComplex& c){

  return EvtVector4C(v2)*=c;
}

inline EvtVector4C operator*(const EvtComplex& c,const EvtVector4R& v2){

  return EvtVector4C(c*v2.get(0),c*v2.get(1),c*v2.get(2),c*v2.get(3));
}

inline EvtVector4C::EvtVector4C(const EvtVector4R& v1){
 
  v[0]=EvtComplex(v1.get(0)); v[1]=EvtComplex(v1.get(1));
  v[2]=EvtComplex(v1.get(2)); v[3]=EvtComplex(v1.get(3));
}

inline EvtComplex operator*(const EvtVector4R& v1,const EvtVector4C& v2){
 
  return v1.get(0)*v2.v[0]-v1.get(1)*v2.v[1]-
         v1.get(2)*v2.v[2]-v1.get(3)*v2.v[3];
}

inline EvtComplex operator*(const EvtVector4C& v1,const EvtVector4R& v2){
 
  return v1.v[0]*v2.get(0)-v1.v[1]*v2.get(1)-
         v1.v[2]*v2.get(2)-v1.v[3]*v2.get(3);
}

inline EvtComplex operator*(const EvtVector4C& v1,const EvtVector4C& v2){
 
  return v1.v[0]*v2.v[0]-v1.v[1]*v2.v[1]-
         v1.v[2]*v2.v[2]-v1.v[3]*v2.v[3];
}

inline EvtVector4C EvtVector4C::conj() const { 

  return EvtVector4C(::conj(v[0]),::conj(v[1]),
		  ::conj(v[2]),::conj(v[3]));
}

#endif

