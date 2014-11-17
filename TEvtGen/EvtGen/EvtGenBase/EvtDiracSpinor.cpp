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
// Module: EvtDiracSpinor.cc
//
// Description:  Class to describe (EvtDiracParticle) spinors.
//
// Modification history:
//
//    DJL/RYD     September 25, 1996         Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <math.h>
#include <assert.h>
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtGammaMatrix.hh"
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtTensor4C.hh"
using std::ostream;


EvtDiracSpinor::~EvtDiracSpinor(){}

EvtDiracSpinor::EvtDiracSpinor(const EvtComplex& sp0,const EvtComplex& sp1,
				    const EvtComplex& sp2,const EvtComplex& sp3){
  set(sp0,sp1,sp2,sp3);
}

void EvtDiracSpinor::set(const EvtComplex& sp0,const EvtComplex& sp1,
			 const EvtComplex& sp2,const EvtComplex& sp3){

  spinor[0]=sp0;spinor[1]=sp1;spinor[2]=sp2;spinor[3]=sp3;
}

void EvtDiracSpinor::set_spinor(int i,const EvtComplex& sp){

  spinor[i]=sp;
}

ostream& operator<<(ostream& s, const EvtDiracSpinor& sp){

  s <<"["<<sp.spinor[0]<<","<<sp.spinor[1]<<","
    <<sp.spinor[2]<<","<<sp.spinor[3]<<"]";
  return s;

}


const EvtComplex& EvtDiracSpinor::get_spinor(int i) const { 
   
  return spinor[i];

}

EvtDiracSpinor rotateEuler(const EvtDiracSpinor& sp,
			   double alpha,double beta,double gamma){

  EvtDiracSpinor tmp(sp);
  tmp.applyRotateEuler(alpha,beta,gamma);
  return tmp;

}

EvtDiracSpinor boostTo(const EvtDiracSpinor& sp,
		       const EvtVector4R p4){

  EvtDiracSpinor tmp(sp);
  tmp.applyBoostTo(p4);
  return tmp;

}

EvtDiracSpinor boostTo(const EvtDiracSpinor& sp,
		       const EvtVector3R boost){

  EvtDiracSpinor tmp(sp);
  tmp.applyBoostTo(boost);
  return tmp;

}

void EvtDiracSpinor::applyBoostTo(const EvtVector4R& p4){

  double e=p4.get(0);

  EvtVector3R boost(p4.get(1)/e,p4.get(2)/e,p4.get(3)/e);

  applyBoostTo(boost);

  return;

}



void EvtDiracSpinor::applyBoostTo(const EvtVector3R& boost) {

  double bx,by,bz,gamma,b2,f1,f2;
  EvtComplex spinorp[4];

  bx=boost.get(0);
  by=boost.get(1);
  bz=boost.get(2);
  b2=bx*bx+by*by+bz*bz;

  if (b2==0.0){
    return;
  }

  //assert(b2<1.0);

  gamma=1.0;
  if (b2 < 1.0) {gamma = 1.0/sqrt(1.0-b2);}
  
  f1=sqrt((gamma+1.0)/2.0);
  f2=f1*gamma/(gamma+1.0);

  spinorp[0]=f1*spinor[0]+f2*bz*spinor[2]+
    f2*EvtComplex(bx,-by)*spinor[3];
  spinorp[1]=f1*spinor[1]+f2*EvtComplex(bx,by)*spinor[2]-
    f2*bz*spinor[3];
  spinorp[2]=f2*bz*spinor[0]+f2*EvtComplex(bx,-by)*spinor[1]+
    f1*spinor[2];
  spinorp[3]=f2*EvtComplex(bx,by)*spinor[0]-
    f2*bz*spinor[1]+f1*spinor[3];
  
  spinor[0]=spinorp[0];
  spinor[1]=spinorp[1];
  spinor[2]=spinorp[2];
  spinor[3]=spinorp[3];

  return;
}

void EvtDiracSpinor::applyRotateEuler(double alpha,double beta,
				      double gamma) {

  EvtComplex retVal[4];
  
  double cb2=cos(0.5*beta);
  double sb2=sin(0.5*beta);
  double capg2=cos(0.5*(alpha+gamma));
  double camg2=cos(0.5*(alpha-gamma));
  double sapg2=sin(0.5*(alpha+gamma));
  double samg2=sin(0.5*(alpha-gamma));

  EvtComplex m11(cb2*capg2,-cb2*sapg2);
  EvtComplex m12(-sb2*camg2,sb2*samg2);
  EvtComplex m21(sb2*camg2,sb2*samg2);
  EvtComplex m22(cb2*capg2,cb2*sapg2);

  retVal[0]=m11*spinor[0]+m12*spinor[1];
  retVal[1]=m21*spinor[0]+m22*spinor[1];
  retVal[2]=m11*spinor[2]+m12*spinor[3];
  retVal[3]=m21*spinor[2]+m22*spinor[3];

  spinor[0]=retVal[0];
  spinor[1]=retVal[1];
  spinor[2]=retVal[2];
  spinor[3]=retVal[3];

  return;
}



EvtDiracSpinor EvtDiracSpinor::conj() const {

  EvtDiracSpinor sp;

  for ( int i=0; i<4; i++)
    sp.set_spinor(i,::conj(spinor[i]));
  
  return sp;
}

EvtVector4C EvtLeptonVACurrent(const EvtDiracSpinor& d,const EvtDiracSpinor& dp){

  //Old code; below is a new specialized code that does it more efficiently.
  //EvtGammaMatrix mat;
  //EvtVector4C temp;
  //mat.va0();
  //temp.set(0,d*(mat*dp));
  //mat.va1();
  //temp.set(1,d*(mat*dp));
  //mat.va2();
  //temp.set(2,d*(mat*dp));
  //mat.va3();
  //temp.set(3,d*(mat*dp));
  //return temp;
 

  EvtComplex u02=::conj(d.spinor[0]-d.spinor[2]);  
  EvtComplex u13=::conj(d.spinor[1]-d.spinor[3]);  

  EvtComplex v02=dp.spinor[0]-dp.spinor[2];
  EvtComplex v13=dp.spinor[1]-dp.spinor[3];

  EvtComplex a=u02*v02;
  EvtComplex b=u13*v13;

  EvtComplex c=u02*v13;
  EvtComplex e=u13*v02;

  return EvtVector4C(a+b,-(c+e),EvtComplex(0,1)*(c-e),b-a);

  
}

EvtVector4C EvtLeptonVCurrent(const EvtDiracSpinor& d,const EvtDiracSpinor& dp){

  EvtVector4C temp;

  // no conjugate here; done in the multiplication
  // yes this is stupid and fooled me to for a long time (ryd)

  temp.set(0,d*(EvtGammaMatrix::v0()*dp));
  temp.set(1,d*(EvtGammaMatrix::v1()*dp));
  temp.set(2,d*(EvtGammaMatrix::v2()*dp));
  temp.set(3,d*(EvtGammaMatrix::v3()*dp));
  
  return temp;
}


EvtVector4C EvtLeptonACurrent(const EvtDiracSpinor& d,const EvtDiracSpinor& dp){

  EvtVector4C temp;

  EvtGammaMatrix mat;

  // no conjugate here; done in the multiplication
  // yes this is stupid and fooled me to for a long time (ryd)

  mat = EvtGammaMatrix::v0()-EvtGammaMatrix::va0();
  temp.set(0,d*(mat*dp));

  mat = EvtGammaMatrix::v1()-EvtGammaMatrix::va1();
  temp.set(1,d*(mat*dp));

  mat = EvtGammaMatrix::v2()-EvtGammaMatrix::va2();
  temp.set(2,d*(mat*dp));

  mat = EvtGammaMatrix::v3()-EvtGammaMatrix::va3();
  temp.set(3,d*(mat*dp));
  
  return temp;
}

EvtComplex EvtLeptonSCurrent(const EvtDiracSpinor& d,const EvtDiracSpinor& dp){

  EvtComplex temp;

  // no conjugate here; done in the multiplication
  // yes this is stupid and fooled me to for a long time (ryd)

  temp=d*(EvtGammaMatrix::g0()*dp);
  
  return temp;
}

EvtComplex EvtLeptonPCurrent(const EvtDiracSpinor& d,const EvtDiracSpinor& dp){

  EvtComplex temp;

  // no conjugate here; done in the multiplication
  // yes this is stupid and fooled me to for a long time (ryd)
  static EvtGammaMatrix m=EvtGammaMatrix::g0()*EvtGammaMatrix::g5();
  temp=d*(m*dp);
  
  return temp;
}

EvtTensor4C EvtLeptonTCurrent(const EvtDiracSpinor& d,const EvtDiracSpinor& dp){

  EvtTensor4C temp;
  temp.zero();
  EvtComplex i2(0,0.5);

  static EvtGammaMatrix mat01=EvtGammaMatrix::g0()*
    (EvtGammaMatrix::g0()*EvtGammaMatrix::g1()-
     EvtGammaMatrix::g1()*EvtGammaMatrix::g0());
  static EvtGammaMatrix mat02=EvtGammaMatrix::g0()*
    (EvtGammaMatrix::g0()*EvtGammaMatrix::g2()-
     EvtGammaMatrix::g2()*EvtGammaMatrix::g0());
  static EvtGammaMatrix mat03=EvtGammaMatrix::g0()*
    (EvtGammaMatrix::g0()*EvtGammaMatrix::g3()-
    EvtGammaMatrix::g3()*EvtGammaMatrix::g0());
  static EvtGammaMatrix mat12=EvtGammaMatrix::g0()*
    (EvtGammaMatrix::g1()*EvtGammaMatrix::g2()-
    EvtGammaMatrix::g2()*EvtGammaMatrix::g1());
  static EvtGammaMatrix mat13=EvtGammaMatrix::g0()*
    (EvtGammaMatrix::g1()*EvtGammaMatrix::g3()-
     EvtGammaMatrix::g3()*EvtGammaMatrix::g1());
  static EvtGammaMatrix mat23=EvtGammaMatrix::g0()*
    (EvtGammaMatrix::g2()*EvtGammaMatrix::g3()-
     EvtGammaMatrix::g3()*EvtGammaMatrix::g2());

 
  temp.set(0,1,i2*(d*(mat01*dp)));
  temp.set(1,0,-temp.get(0,1));

  temp.set(0,2,i2*(d*(mat02*dp)));
  temp.set(2,0,-temp.get(0,2));

  temp.set(0,3,i2*(d*(mat03*dp)));
  temp.set(3,0,-temp.get(0,3));

  temp.set(1,2,i2*(d*(mat12*dp)));
  temp.set(2,1,-temp.get(1,2));

  temp.set(1,3,i2*(d*(mat13*dp)));
  temp.set(3,1,-temp.get(1,3));

  temp.set(2,3,i2*(d*(mat23*dp)));
  temp.set(3,2,-temp.get(2,3));
  
  return temp;
}


EvtDiracSpinor operator*(const EvtComplex& c, const EvtDiracSpinor& d) {
     EvtDiracSpinor result;
     result.spinor[0] = c*d.spinor[0];
     result.spinor[1] = c*d.spinor[1];
     result.spinor[2] = c*d.spinor[2];
     result.spinor[3] = c*d.spinor[3];
     return result;
 }

EvtDiracSpinor EvtDiracSpinor::adjoint() const
{
    EvtDiracSpinor d = this->conj(); // first conjugate, then multiply with gamma0
    EvtGammaMatrix g0 = EvtGammaMatrix::g0();
    EvtDiracSpinor result; // automatically initialized to 0

    for (int i=0; i<4; ++i)
        for (int j=0; j<4; ++j)
            result.spinor[i] += d.spinor[j] * g0._gamma[i][j];

    return result;
}

EvtComplex operator*(const EvtDiracSpinor& d,const EvtDiracSpinor& dp){

  int i;
  EvtComplex temp;
  
  temp=EvtComplex(0.0,0.0);
  
  for(i=0;i<4;i++){
    temp += conj( d.get_spinor(i) ) * dp.get_spinor( i ) ;
  }
  return temp;
}
