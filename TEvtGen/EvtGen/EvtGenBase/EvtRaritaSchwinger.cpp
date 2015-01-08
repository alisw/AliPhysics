//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2000      Caltech, UCSB
//
// Module: EvtGen/EvtRaritaSchwing.hh
//
// Description:Class to handle spin 3/2
//
// Modification history:
//
//    RYD     April 23, 2000         Module created
//
//------------------------------------------------------------------------
#include "EvtGenBase/EvtPatches.hh"


#include "EvtGenBase/EvtRaritaSchwinger.hh"
using std::endl;
using std::ostream;

EvtRaritaSchwinger::~EvtRaritaSchwinger(){}


EvtRaritaSchwinger rotateEuler(const EvtRaritaSchwinger& rs,
			       double alpha,double beta,double gamma){

  EvtRaritaSchwinger tmp(rs);
  tmp.applyRotateEuler(alpha,beta,gamma);
  return tmp;

}

EvtRaritaSchwinger boostTo(const EvtRaritaSchwinger& rs,
			   const EvtVector4R p4){

  EvtRaritaSchwinger tmp(rs);
  tmp.applyBoostTo(p4);
  return tmp;

}

EvtRaritaSchwinger boostTo(const EvtRaritaSchwinger& rs,
			   const EvtVector3R boost){

  EvtRaritaSchwinger tmp(rs);
  tmp.applyBoostTo(boost);
  return tmp;

}


void EvtRaritaSchwinger::set(int i,int j,const EvtComplex& sp){_rs[i][j]=sp;}

EvtComplex EvtRaritaSchwinger::get(int i,int j) const {return _rs[i][j];} 

void EvtRaritaSchwinger::applyRotateEuler(double alpha,double beta,
					  double gamma){

  //inefficient but simple to code...
  EvtVector4C v0=getVector(0);
  EvtVector4C v1=getVector(1);
  EvtVector4C v2=getVector(2);
  EvtVector4C v3=getVector(3);
  v0.applyRotateEuler(alpha,beta,gamma);
  v1.applyRotateEuler(alpha,beta,gamma);
  v2.applyRotateEuler(alpha,beta,gamma);
  v3.applyRotateEuler(alpha,beta,gamma);
  setVector(0,v0);
  setVector(1,v1);
  setVector(2,v2);
  setVector(3,v3);
  EvtDiracSpinor sp0=getSpinor(0);
  EvtDiracSpinor sp1=getSpinor(1);
  EvtDiracSpinor sp2=getSpinor(2);
  EvtDiracSpinor sp3=getSpinor(3);
  sp0.applyRotateEuler(alpha,beta,gamma);
  sp1.applyRotateEuler(alpha,beta,gamma);
  sp2.applyRotateEuler(alpha,beta,gamma);
  sp3.applyRotateEuler(alpha,beta,gamma);
  setSpinor(0,sp0);
  setSpinor(1,sp1);
  setSpinor(2,sp2);
  setSpinor(3,sp3);

}
  

void EvtRaritaSchwinger::applyBoostTo(const EvtVector4R p4){

  double e=p4.get(0);

  EvtVector3R boost(p4.get(1)/e,p4.get(2)/e,p4.get(3)/e);

  applyBoostTo(boost);
  
  return;

}
  

void EvtRaritaSchwinger::applyBoostTo(const EvtVector3R boost){

  //inefficient but simple to code...
  EvtVector4C v0=getVector(0);
  EvtVector4C v1=getVector(1);
  EvtVector4C v2=getVector(2);
  EvtVector4C v3=getVector(3);
  v0.applyBoostTo(boost);
  v1.applyBoostTo(boost);
  v2.applyBoostTo(boost);
  v3.applyBoostTo(boost);
  setVector(0,v0);
  setVector(1,v1);
  setVector(2,v2);
  setVector(3,v3);
  EvtDiracSpinor sp0=getSpinor(0);
  EvtDiracSpinor sp1=getSpinor(1);
  EvtDiracSpinor sp2=getSpinor(2);
  EvtDiracSpinor sp3=getSpinor(3);
  sp0.applyBoostTo(boost);
  sp1.applyBoostTo(boost);
  sp2.applyBoostTo(boost);
  sp3.applyBoostTo(boost);
  setSpinor(0,sp0);
  setSpinor(1,sp1);
  setSpinor(2,sp2);
  setSpinor(3,sp3);


}


ostream& operator<<(ostream& s, const EvtRaritaSchwinger& rs){

  int i,j;
  s<< endl;
  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      s << rs._rs[i][j];
    }
    s << endl;
  }
  return s;
  
}



EvtVector4C EvtRaritaSchwinger::getVector(int i) const{

  EvtVector4C tmp(_rs[i][0],_rs[i][1],_rs[i][2],_rs[i][3]);
  return tmp;

}
 
EvtDiracSpinor EvtRaritaSchwinger::getSpinor(int i) const{

  EvtDiracSpinor tmp;
  tmp.set(_rs[0][i],_rs[1][i],_rs[2][i],_rs[3][i]);
  return tmp;

}

void EvtRaritaSchwinger::setVector(int i,const EvtVector4C& v){
  
  _rs[i][0]=v.get(0);
  _rs[i][1]=v.get(1);
  _rs[i][2]=v.get(2);
  _rs[i][3]=v.get(3);

}

void EvtRaritaSchwinger::setSpinor(int i,const EvtDiracSpinor& sp){

  _rs[0][i]=sp.get_spinor(0);
  _rs[1][i]=sp.get_spinor(1);
  _rs[2][i]=sp.get_spinor(2);
  _rs[3][i]=sp.get_spinor(3);

}


EvtRaritaSchwinger dirProd(EvtVector4R v,EvtDiracSpinor u){

  int i,j;

  EvtRaritaSchwinger tmp;

  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      tmp._rs[i][j]=u.get_spinor(i)*v.get(j);
    }
  }

  return tmp;

}


EvtRaritaSchwinger dirProd(EvtVector4C v,EvtDiracSpinor u){

  int i,j;

  EvtRaritaSchwinger tmp;

  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      tmp._rs[i][j]=u.get_spinor(i)*v.get(j);
    }
  }

  return tmp;

}


EvtComplex operator*(const EvtRaritaSchwinger& u1,
		     const EvtRaritaSchwinger& u2){

  int i,j;
  EvtComplex tmp=0.0;

  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      tmp+=conj(u1._rs[i][j])*u2._rs[i][j];
    }
  }

  return tmp;

}



EvtRaritaSchwinger& EvtRaritaSchwinger::operator+=(const EvtRaritaSchwinger& u2){

  int i,j;

  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      _rs[i][j]+=u2._rs[i][j];
    }
  }
  
  return *this; 
}

EvtRaritaSchwinger operator+(const EvtRaritaSchwinger& u1,
				const EvtRaritaSchwinger& u2){
  
  return EvtRaritaSchwinger(u1)+=u2;

}

EvtRaritaSchwinger& EvtRaritaSchwinger::operator-=(const EvtRaritaSchwinger& u2){

  int i,j;

  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      _rs[i][j]+=u2._rs[i][j];
    }
  }
  
  return *this; 
}

EvtRaritaSchwinger operator-(const EvtRaritaSchwinger& u1,
				const EvtRaritaSchwinger& u2){
  
  return EvtRaritaSchwinger(u1)-=u2;

}


