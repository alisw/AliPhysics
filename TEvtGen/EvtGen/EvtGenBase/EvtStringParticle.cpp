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
// Module: EvtStringParticle.cc
//
// Description: Class to describe the partons that are produced in JetSet.
//
// Modification history:
//
//    RYD       Febuary 27,1998       Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include "EvtGenBase/EvtStringParticle.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtReport.hh"


EvtStringParticle::~EvtStringParticle(){

  if (_npartons!=0){
    
    delete [] _p4partons;
    delete [] _idpartons;
    
  }

}


EvtStringParticle::EvtStringParticle(){

  _p4partons=0;
  _idpartons=0;
  _npartons=0;

  return;
  
}

void EvtStringParticle::init(EvtId id, const EvtVector4R& p4){

  _validP4=true;
  setp(p4);
  setpart_num(id);

}


void EvtStringParticle::initPartons(int npartons,
			      EvtVector4R* p4partons,EvtId* idpartons){

  _p4partons = new EvtVector4R[npartons];
  _idpartons = new EvtId[npartons];

  int i;

  _npartons=npartons;

  for(i=0;i<npartons;i++){

    _p4partons[i]=p4partons[i];
    _idpartons[i]=idpartons[i];

  }

}

int EvtStringParticle::getNPartons(){

  return _npartons;

}

EvtId EvtStringParticle::getIdParton(int i){

  return _idpartons[i];

}

EvtVector4R EvtStringParticle::getP4Parton(int i){

  return _p4partons[i];

}
  
EvtSpinDensity EvtStringParticle::rotateToHelicityBasis() const{

  report(Severity::Error,"EvtGen") << "rotateToHelicityBasis not implemented for strin particle."; 
  report(Severity::Error,"EvtGen") << "Will terminate execution."; 

  ::abort();

  EvtSpinDensity rho;
  return rho;

}

EvtSpinDensity EvtStringParticle::rotateToHelicityBasis(double,
							double,
							double) const{

  report(Severity::Error,"EvtGen") << "rotateToHelicityBasis(alpha,beta,gamma) not implemented for string particle."; 
  report(Severity::Error,"EvtGen") << "Will terminate execution."; 

  ::abort();

  EvtSpinDensity rho;
  return rho;

}


