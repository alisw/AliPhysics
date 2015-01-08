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
// Module: EvtVectorParticle.cc
//
// Description: Class to describe spin 1 particles
//
// Modification history:
//
//    DJL/RYD   September 25, 1996           Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtVectorParticle.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"

EvtVectorParticle::~EvtVectorParticle(){}


void EvtVectorParticle::init(EvtId part_n,double e,double px,double py,double pz){

  _validP4=true;
  setp(e,px,py,pz);
  setpart_num(part_n);
  
  _eps[0].set(0.0,1.0,0.0,0.0);
  _eps[1].set(0.0,0.0,1.0,0.0);
  _eps[2].set(0.0,0.0,0.0,1.0);

  setLifetime();
}

void EvtVectorParticle::init(EvtId part_n,const EvtVector4R& p4){

  _validP4=true;
  setp(p4);
  setpart_num(part_n);
  
  _eps[0].set(0.0,1.0,0.0,0.0);
  _eps[1].set(0.0,0.0,1.0,0.0);
  _eps[2].set(0.0,0.0,0.0,1.0);
  setLifetime();
}

void EvtVectorParticle::init(EvtId part_n,const EvtVector4R& p4,
			     const EvtVector4C & epsin1,
			     const EvtVector4C & epsin2,
			     const EvtVector4C & epsin3){

  _validP4=true;
  setp(p4);
  setpart_num(part_n);
  
  _eps[0]=epsin1;
  _eps[1]=epsin2;
  _eps[2]=epsin3;

  setLifetime();
}



EvtSpinDensity EvtVectorParticle::rotateToHelicityBasis() const{

  static EvtVector4C eplus(0.0,-1.0/sqrt(2.0),EvtComplex(0.0,-1.0/sqrt(2.0)),0.0);
  static EvtVector4C ezero(0.0,0.0,0.0,1.0);
  static EvtVector4C eminus(0.0,1.0/sqrt(2.0),EvtComplex(0.0,-1.0/sqrt(2.0)),0.0);
  
  static EvtVector4C eplusC(eplus.conj());
  static EvtVector4C ezeroC(ezero.conj());
  static EvtVector4C eminusC(eminus.conj());

  EvtSpinDensity R;
  R.setDim(3);

  for ( int i=0; i<3; i++ ) {
    R.set(0,i,(eplusC)*_eps[i]);
    R.set(1,i,(ezeroC)*_eps[i]);
    R.set(2,i,(eminusC)*_eps[i]);
  }

  return R;

}


EvtSpinDensity EvtVectorParticle::rotateToHelicityBasis(double alpha,
							double beta,
							double gamma) const{

  EvtVector4C eplus(0.0,-1.0/sqrt(2.0),EvtComplex(0.0,-1.0/sqrt(2.0)),0.0);
  EvtVector4C ezero(0.0,0.0,0.0,1.0);
  EvtVector4C eminus(0.0,1.0/sqrt(2.0),EvtComplex(0.0,-1.0/sqrt(2.0)),0.0);
  
  eplus.applyRotateEuler(alpha,beta,gamma);
  ezero.applyRotateEuler(alpha,beta,gamma);
  eminus.applyRotateEuler(alpha,beta,gamma);

  EvtSpinDensity R;
  R.setDim(3);


  for ( int i=0; i<3; i++ ) {
    R.set(0,i,(eplus.conj())*_eps[i]);
    R.set(1,i,(ezero.conj())*_eps[i]);
    R.set(2,i,(eminus.conj())*_eps[i]);
  }

  return R;

}


