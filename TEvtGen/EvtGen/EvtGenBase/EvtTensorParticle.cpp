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
// Module: EvtTensorParticle.cc
//
// Description: Class to describe spin 2 particles.
//
// Modification history:
//
//    DJL/RYD   September 25,1996           Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtTensorParticle.hh"
#include "EvtGenBase/EvtReport.hh"

EvtTensorParticle::~EvtTensorParticle(){}

void EvtTensorParticle::init(EvtId part_n,const EvtVector4R& p4){
 
  init(part_n,p4.get(0),p4.get(1)
       ,p4.get(2),p4.get(3));

  setLifetime();

  
}

void EvtTensorParticle::init(EvtId part_n,double e,double px,double py,double pz){

  _validP4=true;
  setp(e,px,py,pz);
  setpart_num(part_n);
  
  eps[0].setdiag(0.0,-1.0/sqrt(6.0),-1.0/sqrt(6.0),
	       2.0/sqrt(6.0));
  eps[1].setdiag(0.0,1.0/sqrt(2.0),-1.0/sqrt(2.0),0.0);
  eps[2].setdiag(0.0,0.0,0.0,0.0);
  eps[3].setdiag(0.0,0.0,0.0,0.0);
  eps[4].setdiag(0.0,0.0,0.0,0.0);

  eps[2].set(1,2,EvtComplex(1.0/sqrt(2.0),0.0));
  eps[2].set(2,1,EvtComplex(1.0/sqrt(2.0),0.0));
  eps[3].set(1,3,EvtComplex(1.0/sqrt(2.0),0.0));
  eps[3].set(3,1,EvtComplex(1.0/sqrt(2.0),0.0));
  eps[4].set(2,3,EvtComplex(1.0/sqrt(2.0),0.0));
  eps[4].set(3,2,EvtComplex(1.0/sqrt(2.0),0.0));

  setLifetime();
  
}


void EvtTensorParticle::init(EvtId part_n,const EvtVector4R& p4,
			     const EvtTensor4C& epsin1, 
			     const EvtTensor4C& epsin2,
			     const EvtTensor4C& epsin3, 
			     const EvtTensor4C& epsin4,
			     const EvtTensor4C& epsin5){
 
  _validP4=true;
  setp(p4);
  setpart_num(part_n);

  eps[0]=epsin1;
  eps[1]=epsin2;
  eps[2]=epsin3;
  eps[3]=epsin4;
  eps[4]=epsin5;

  setLifetime();
  
}



EvtTensor4C EvtTensorParticle::epsTensorParent(int i) const {

  EvtTensor4C temp=eps[i];

  temp.applyBoostTo(this->getP4());
  return temp;
  
} //epsParent


EvtTensor4C EvtTensorParticle::epsTensor(int i) const {
   
  return eps[i];

} //eps



EvtSpinDensity EvtTensorParticle::rotateToHelicityBasis() const{


  static EvtVector4C eplus(0.0,-1.0/sqrt(2.0),EvtComplex(0.0,-1.0/sqrt(2.0)),0.0);
  static EvtVector4C ezero(0.0,0.0,0.0,1.0);
  static EvtVector4C eminus(0.0,1.0/sqrt(2.0),EvtComplex(0.0,-1.0/sqrt(2.0)),0.0);
  
  static EvtTensor4C dPpp(EvtGenFunctions::directProd(eplus,eplus));
  static EvtTensor4C dPp0(EvtGenFunctions::directProd(eplus,ezero));
  static EvtTensor4C dP0p(EvtGenFunctions::directProd(ezero,eplus));
  static EvtTensor4C dPpm(EvtGenFunctions::directProd(eplus,eminus));
  static EvtTensor4C dP00(EvtGenFunctions::directProd(ezero,ezero));
  static EvtTensor4C dPmp(EvtGenFunctions::directProd(eminus,eplus));
  static EvtTensor4C dPmm(EvtGenFunctions::directProd(eminus,eminus));
  static EvtTensor4C dPm0(EvtGenFunctions::directProd(eminus,ezero));
  static EvtTensor4C dP0m(EvtGenFunctions::directProd(ezero,eminus));

  static EvtTensor4C es0(conj(dPpp));
  static EvtTensor4C es1(conj((1/sqrt(2.0))*dPp0 +(1/sqrt(2.0))*dP0p));
  static EvtTensor4C es2(conj((1/sqrt(6.0))*dPpm +(2/sqrt(6.0))*dP00 +(1/sqrt(6.0))*dPmp));
  static EvtTensor4C es3(conj((1/sqrt(2.0))*dPm0  +(1/sqrt(2.0))*dP0m));
  static EvtTensor4C es4(conj(dPmm));


  EvtSpinDensity R;
  R.setDim(5);

  for (int j=0; j<5; j++) {
    R.set(0,j,cont(es0,eps[j]));
    R.set(1,j,cont(es1,eps[j]));
    R.set(2,j,cont(es2,eps[j]));
    R.set(3,j,cont(es3,eps[j]));
    R.set(4,j,cont(es4,eps[j]));
  }
  return R;

}


EvtSpinDensity EvtTensorParticle::rotateToHelicityBasis(double alpha,
							double beta,
							double gamma) const{

  EvtTensor4C es[5];

  static EvtVector4C eplus(0.0,-1.0/sqrt(2.0),EvtComplex(0.0,-1.0/sqrt(2.0)),0.0);
  static EvtVector4C ezero(0.0,0.0,0.0,1.0);
  static EvtVector4C eminus(0.0,1.0/sqrt(2.0),EvtComplex(0.0,-1.0/sqrt(2.0)),0.0);

  eplus.applyRotateEuler(alpha,beta,gamma);
  ezero.applyRotateEuler(alpha,beta,gamma);
  eminus.applyRotateEuler(alpha,beta,gamma);

  for (int i=0; i<5; i++) es[i].zero();    
  
  es[0]=EvtGenFunctions::directProd(eplus,eplus);
  es[1] =(1/sqrt(2.0))*EvtGenFunctions::directProd(eplus,ezero)
    +(1/sqrt(2.0))*EvtGenFunctions::directProd(ezero,eplus);
  es[2] =(1/sqrt(6.0))*EvtGenFunctions::directProd(eplus,eminus)
    +(2/sqrt(6.0))*EvtGenFunctions::directProd(ezero,ezero)
    +(1/sqrt(6.0))*EvtGenFunctions::directProd(eminus,eplus);
  es[3] =(1/sqrt(2.0))*EvtGenFunctions::directProd(eminus,ezero)
    +(1/sqrt(2.0))*EvtGenFunctions::directProd(ezero,eminus);
  es[4]=EvtGenFunctions::directProd(eminus,eminus);

  for (int i=0; i<5; i++) es[i]=conj(es[i]);    

  EvtSpinDensity R;
  R.setDim(5);

  for (int i=0; i<5; i++)
    for (int j=0; j<5; j++)
	R.set(i,j,cont(es[i],eps[j]));

  return R;

}

  






