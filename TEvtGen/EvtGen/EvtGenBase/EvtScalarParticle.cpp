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
// Module: EvtScalarParticle.cc
//
// Description: Class to describe scalar particles
//
// Modification history:
//
//    DJL/RYD   September 25, 1996           Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <iostream>
#include <math.h>
#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtScalarParticle.hh"
#include "EvtGenBase/EvtVector4R.hh"


void EvtScalarParticle::init(EvtId part_n,double e,double px,double py,double pz){

  _validP4=true;
  setp(e,px,py,pz);
  setpart_num(part_n);

  setLifetime();

}

EvtScalarParticle::~EvtScalarParticle() {}


void EvtScalarParticle::init(EvtId part_n,const EvtVector4R& p4){

  _validP4=true;
  setp(p4);
  setpart_num(part_n);

  setLifetime();

}

EvtSpinDensity EvtScalarParticle::rotateToHelicityBasis() const{

  EvtSpinDensity R;
  R.setDim(1);
      
  R.set(0,0,1.0);

  return R;

}


EvtSpinDensity EvtScalarParticle::rotateToHelicityBasis(double,
						       double,
						       double) const{

  EvtSpinDensity R;
  R.setDim(1);
      
  R.set(0,0,1.0);

  return R;

}

