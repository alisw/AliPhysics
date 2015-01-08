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
// Module: EvtHighSpinParticle.cc
//
// Description: Class to describe particles with spin>2.
//
// Modification history:
//
//    RYD   August 8, 2000           Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtPatches.hh"
#include <iostream>
#include <math.h>
#include <assert.h>
#include "EvtGenBase/EvtHighSpinParticle.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtSpinDensity.hh"
#include "EvtGenBase/EvtdFunction.hh"


EvtHighSpinParticle::~EvtHighSpinParticle() {}


void EvtHighSpinParticle::init(EvtId id,const EvtVector4R& p4){

  _validP4=true;
  setp(p4);
  setpart_num(id);

  setLifetime();

}

EvtSpinDensity EvtHighSpinParticle::rotateToHelicityBasis() const{

  int n=EvtSpinType::getSpinStates(EvtPDL::getSpinType(getId()));

  EvtSpinDensity R;
  R.setDiag(n);

  return R;

}



EvtSpinDensity EvtHighSpinParticle::rotateToHelicityBasis(double alpha,
							  double beta,
							  double gamma) const{

  int i,j;
  
  int n=EvtSpinType::getSpinStates(EvtPDL::getSpinType(getId()));

  EvtSpinDensity R;
  
  R.setDim(n);

  int J2=EvtSpinType::getSpin2(EvtPDL::getSpinType(getId()));

  assert(n==J2+1);

  int *lambda2;

  lambda2=new int[J2+1];

  for(i=0;i<J2+1;i++){
    lambda2[i]=J2-i*2;
  }


  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      R.set(i,j,EvtdFunction::d(J2,lambda2[j],lambda2[i],beta)*
	exp(EvtComplex(0.0,0.5*(alpha*lambda2[i]-gamma*lambda2[j]))));
    }
  }

  delete [] lambda2;

  return R;


}


