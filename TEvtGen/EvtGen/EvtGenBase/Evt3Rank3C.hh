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
// Module: EvtGen/Evt3Rank3C.hh
//
// Description:Class to handle complex 3rd rank 3D tensors
//
// Modification history:
//
//    RYD     September 14, 1997         Module created
//
//------------------------------------------------------------------------

#ifndef EVT3RANK3C_HH
#define EVT3RANK3C_HH

#include <iostream>
#include "EvtGenBase/EvtComplex.hh"

class EvtTensor3C;
class EvtVector3C;
class EvtVector3R;


class  Evt3Rank3C ;
inline Evt3Rank3C operator*(const EvtComplex& c,const Evt3Rank3C& t2);
inline Evt3Rank3C operator*(const double d,const Evt3Rank3C& t2);
inline Evt3Rank3C operator*(const Evt3Rank3C& t2,const EvtComplex& c);
inline Evt3Rank3C operator*(const Evt3Rank3C& t2,const double d);
inline Evt3Rank3C operator+(const Evt3Rank3C& t1,const Evt3Rank3C& t2);
inline Evt3Rank3C operator-(const Evt3Rank3C& t1,const Evt3Rank3C& t2);
Evt3Rank3C directProd(const EvtVector3C& c1,const EvtVector3C& c2,
		      const EvtVector3C& c3); 
Evt3Rank3C conj(const Evt3Rank3C& t2);


class Evt3Rank3C {

  friend Evt3Rank3C operator*(const EvtComplex& c,const Evt3Rank3C& t2);
  friend Evt3Rank3C operator*(const double d,const Evt3Rank3C& t2);
  friend Evt3Rank3C operator*(const Evt3Rank3C& t2,const EvtComplex& c);
  friend Evt3Rank3C operator*(const Evt3Rank3C& t2,const double d);
  friend Evt3Rank3C operator+(const Evt3Rank3C& t1,const Evt3Rank3C& t2);
  friend Evt3Rank3C operator-(const Evt3Rank3C& t1,const Evt3Rank3C& t2);
  friend Evt3Rank3C directProd(const EvtVector3C& c1,const EvtVector3C& c2,
			      const EvtVector3C& c3); 
  friend Evt3Rank3C conj(const Evt3Rank3C& t2);

  friend std::ostream& operator<<(std::ostream& s, const Evt3Rank3C& t2);
  
public:
  Evt3Rank3C();
  Evt3Rank3C(const Evt3Rank3C& t1 );
  virtual ~Evt3Rank3C();
  Evt3Rank3C& operator=(const Evt3Rank3C& t1);
  inline void set(int i,int j,int k,const EvtComplex& c);
  inline const EvtComplex& get(int i, int j, int k) const;
  void zero();
  

  Evt3Rank3C& operator+=(const Evt3Rank3C& t2);
  Evt3Rank3C& operator-=(const Evt3Rank3C& t2);
  Evt3Rank3C& operator*=(const double d);
  Evt3Rank3C& operator*=(const EvtComplex& c);
  Evt3Rank3C conj() const;
  EvtTensor3C cont1(const EvtVector3C& v) const; 
  EvtTensor3C cont2(const EvtVector3C& v) const; 
  EvtTensor3C cont3(const EvtVector3C& v) const; 
  EvtTensor3C cont1(const EvtVector3R& v) const; 
  EvtTensor3C cont2(const EvtVector3R& v) const; 
  EvtTensor3C cont3(const EvtVector3R& v) const; 
  
  
private:

  EvtComplex t[3][3][3];

};


inline Evt3Rank3C operator*(const EvtComplex& c,const Evt3Rank3C& t2){
   return Evt3Rank3C(t2)*=c;
}

inline Evt3Rank3C operator*(const double d,const Evt3Rank3C& t2){
   return Evt3Rank3C(t2)*=d;
}

inline Evt3Rank3C operator*(const Evt3Rank3C& t2,const EvtComplex& c){
   return Evt3Rank3C(t2)*=c;
}

inline Evt3Rank3C operator*(const Evt3Rank3C& t2,const double d){
   return Evt3Rank3C(t2)*=d;
}

inline Evt3Rank3C operator+(const Evt3Rank3C& t1,const Evt3Rank3C& t2){
   return Evt3Rank3C(t1)+=t2;
}

inline Evt3Rank3C operator-(const Evt3Rank3C& t1,const Evt3Rank3C& t2){
   return Evt3Rank3C(t1)-=t2;
}

inline void Evt3Rank3C::set(int i,int j,int k,const EvtComplex& c){
   t[i][j][k]=c;
}

inline const EvtComplex& Evt3Rank3C::get(int i,int j,int k) const{
   return t[i][j][k];
}


#endif


