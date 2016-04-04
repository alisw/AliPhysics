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
// Module: Evt3Rank3C.cc
//
// Description: Implementation of complex 3Rank 3D tensors.
//
// Modification history:
//
//    RYD     September 14, 1996         Module created
//
//------------------------------------------------------------------------
#include "EvtGenBase/EvtPatches.hh"

#include <iostream>
#include <math.h>
#include "EvtGenBase/Evt3Rank3C.hh"
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtVector3C.hh"
#include "EvtGenBase/EvtTensor3C.hh"
#include "EvtGenBase/EvtReport.hh"



Evt3Rank3C::Evt3Rank3C( const Evt3Rank3C& t1 ) {

  int i,j,k;
  
  for(i=0;i<3;i++) {
    for(j=0;j<3;j++) {
      for(k=0;k<3;j++) {
	t[i][j][k] = t1.t[i][j][k];
      }
    }
  }
}

Evt3Rank3C::~Evt3Rank3C() { }


Evt3Rank3C& Evt3Rank3C::operator=(const Evt3Rank3C& t1) {
  int i,j,k;
  
  for(i=0;i<3;i++) {
    for(j=0;j<3;j++) {
      for(k=0;k<3;k++) {
	t[i][j][k] = t1.t[i][j][k];
      }
    }
  }
  return *this;
}


Evt3Rank3C Evt3Rank3C::conj() const {
  Evt3Rank3C temp;
  
  int i,j,k;
  
  for(i=0;i<3;i++) {
    for(j=0;j<3;j++) {
      for(k=0;k<3;k++) {
	temp.set(i,j,k,::conj(t[i][j][k]));
      }
    }
  }
  return temp;
}

void Evt3Rank3C::zero(){
  int i,j,k;
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      for(k=0;k<3;k++){
	t[i][j][k]=EvtComplex(0.0,0.0);
      }
    }
  }
}


Evt3Rank3C::Evt3Rank3C(){
  
  int i,j,k;
  
  for(i=0;i<3;i++){
    for(j=0;j<3;j++){
      for(k=0;k<3;k++){
	t[i][j][k]=EvtComplex(0.0,0.0);
      }
    }
  }

}

std::ostream& operator<<(std::ostream& s, const Evt3Rank3C& t2){
  int i,j,k;
  for(k=0;k<3;k++){
    for(i=0;i<3;i++){
      for(j=0;j<3;j++){
	report(Severity::Info,"EvtGen") <<t2.t[k][i][j];
      }
      report(Severity::Info,"EvtGen") << "\n";
    }
  }
  return s;
}

Evt3Rank3C& Evt3Rank3C::operator+=(const Evt3Rank3C& t2){
  
  int i,j,k;
  
  for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
      for (k=0;k<3;k++) {
	t[i][j][k]+=t2.t[i][j][k];
      }
    }
  }
  return *this;
}

Evt3Rank3C& Evt3Rank3C::operator-=(const Evt3Rank3C& t2) {

  int i,j,k;
  
  for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
      for (k=0;k<3;k++) {
	t[i][j][k]-=t2.t[i][j][k];
      }
    }
  }

  return *this;

}


Evt3Rank3C& Evt3Rank3C::operator*=(const EvtComplex& c){

  int i,j,k;
  
  for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
      for (k=0;k<3;k++) {
	t[i][j][k]*=c;
      }
    }
  }
  return *this;
}

Evt3Rank3C& Evt3Rank3C::operator*=(const double c){
  int i,j,k;
  
  for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
      for (k=0;k<3;k++) {
	t[i][j][k]*=c;
      }
    }
  }

  return *this;

}

Evt3Rank3C conj(const Evt3Rank3C& t2) { 
  Evt3Rank3C temp;
  
  int i,j,k;

  for(i=0;i<3;i++){ 
    for(j=0;j<3;j++){ 
      for(k=0;k<3;k++){ 
	temp.t[i][j][k]=::conj(t2.t[i][j][k]);
      }
    }
  }
  return temp;
}

EvtTensor3C Evt3Rank3C::cont1(const EvtVector3C& v) const {
  EvtTensor3C temp;
  
  int i,k;
  
  for(i=0;i<3;i++){
    for(k=0;k<3;k++){
      temp.set(i,k,t[0][i][k]*v.get(0)+t[1][i][k]*v.get(1)
	     +t[2][i][k]*v.get(2));
    }
  }
  return temp;
} 


EvtTensor3C Evt3Rank3C::cont2(const EvtVector3C& v) const {
  EvtTensor3C temp;
  
  int i,k;
  
  for(i=0;i<3;i++){
    for(k=0;k<3;k++){
      temp.set(i,k,t[i][0][k]*v.get(0)+t[i][1][k]*v.get(1)
	     +t[i][2][k]*v.get(2));
    }
  }
  return temp;
} 


EvtTensor3C Evt3Rank3C::cont3(const EvtVector3C& v) const {
  EvtTensor3C temp;
  
  int i,k;
  
  for(i=0;i<3;i++){
    for(k=0;k<3;k++){
      temp.set(i,k,t[i][k][0]*v.get(0)+t[i][k][1]*v.get(1)
	     +t[i][k][2]*v.get(2));
    }
  }
  return temp;
} 

EvtTensor3C Evt3Rank3C::cont1(const EvtVector3R& v) const {
  EvtTensor3C temp;
  
  int i,k;
  
  for(i=0;i<3;i++){
    for(k=0;k<3;k++){
      temp.set(i,k,t[0][i][k]*v.get(0)+t[1][i][k]*v.get(1)
	     +t[2][i][k]*v.get(2));
    }
  }
  return temp;
} 


EvtTensor3C Evt3Rank3C::cont2(const EvtVector3R& v) const {
  EvtTensor3C temp;
  
  int i,k;
  
  for(i=0;i<3;i++){
    for(k=0;k<3;k++){
      temp.set(i,k,t[i][0][k]*v.get(0)+t[i][1][k]*v.get(1)
	     +t[i][2][k]*v.get(2));
    }
  }
  return temp;
} 


EvtTensor3C Evt3Rank3C::cont3(const EvtVector3R& v) const {
  EvtTensor3C temp;
  
  int i,k;
  
  for(i=0;i<3;i++){
    for(k=0;k<3;k++){
      temp.set(i,k,t[i][k][0]*v.get(0)+t[i][k][1]*v.get(1)
	     +t[i][k][2]*v.get(2));
    }
  }
  return temp;
} 

