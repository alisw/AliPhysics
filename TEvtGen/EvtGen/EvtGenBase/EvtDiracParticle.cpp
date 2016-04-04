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
// Module: EvtDiracParticle.cc
//
// Description: Class to describe spin 1/2 particles.
//
// Modification history:
//
//    DJL/RYD     September 25, 1996         Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtDiracParticle.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtSpinDensity.hh"
#include "EvtGenBase/EvtGammaMatrix.hh"
using std::endl;

EvtDiracParticle::~EvtDiracParticle(){}


EvtDiracParticle::EvtDiracParticle(){

  return;
  
}

void EvtDiracParticle::init(EvtId part_n,const EvtVector4R& p4){

  _validP4=true; 
  setp(p4);
  setpart_num(part_n);

  if (EvtPDL::getStdHep(part_n)==0){
    report(Severity::Error,"EvtGen") << "Error in EvtDiracParticle::init, part_n="
                           << part_n.getId()<<endl;
    ::abort();
  }

  if (EvtPDL::getStdHep(part_n)>0){  

    _spinorRest[0].set(EvtComplex(sqrt(2.0*mass()),0.0),EvtComplex(0.0,0.0),
		       EvtComplex(0.0,0.0),EvtComplex(0.0,0.0));
    _spinorRest[1].set(EvtComplex(0.0,0.0),EvtComplex(sqrt(2.0*mass()),0.0),
		       EvtComplex(0.0,0.0),EvtComplex(0.0,0.0));
  
    _spinorParent[0]=boostTo(_spinorRest[0],p4);
    _spinorParent[1]=boostTo(_spinorRest[1],p4);


  }
  else{

    _spinorRest[0].set(EvtComplex(0.0,0.0),EvtComplex(0.0,0.0),
		     EvtComplex(sqrt(2.0*mass()),0.0),EvtComplex(0.0,0.0));
    _spinorRest[1].set(EvtComplex(0.0,0.0),EvtComplex(0.0,0.0),
		     EvtComplex(0.0,0.0),EvtComplex(sqrt(2.0*mass()),0.0));
  
    _spinorParent[0]=boostTo(_spinorRest[0],p4);
    _spinorParent[1]=boostTo(_spinorRest[1],p4);



  }

  setLifetime();
}


void EvtDiracParticle::init(EvtId part_n,const EvtVector4R& p4,
			    const EvtDiracSpinor & prod1, 
			    const EvtDiracSpinor & prod2,
			    const EvtDiracSpinor & rest1, 
			    const EvtDiracSpinor & rest2){
  
  _validP4=true; 
  setp(p4);
  setpart_num(part_n);
  
  if (EvtPDL::getStdHep(part_n)==0){
    report(Severity::Error,"EvtGen") << "Error in EvtDiracParticle::init, part_n="
                           << part_n.getId()<<std::endl;
    ::abort();
  }
  _spinorRest[0]=rest1;
  _spinorRest[1]=rest2;
  _spinorParent[0]=prod1;
  _spinorParent[1]=prod2;
  
  setLifetime();
}
 


EvtSpinDensity EvtDiracParticle::rotateToHelicityBasis() const{

  EvtDiracSpinor spplus;
  EvtDiracSpinor spminus;
      
  double sqmt2=sqrt(2.*(getP4().mass()));
      
  if (EvtPDL::getStdHep(getId())>0){  
    spplus.set(1.0,0.0,0.0,0.0);
    spminus.set(0.0,1.0,0.0,0.0);
  } else {
    spplus.set(0.0,0.0,0.0,1.0);
    spminus.set(0.0,0.0,1.0,0.0);
  }

      
  EvtSpinDensity R;
  R.setDim(2);
      
  for (int i=0; i<2; i++) {
    if (EvtPDL::getStdHep(getId())>0){  
      R.set(0,i,(spplus*_spinorRest[i])/sqmt2);
      R.set(1,i,(spminus*_spinorRest[i])/sqmt2);
    } else {
      R.set(0,i,(_spinorRest[i]*spplus)/sqmt2);
      R.set(1,i,(_spinorRest[i]*spminus)/sqmt2);
    } 
  }

  return R;

}


EvtSpinDensity EvtDiracParticle::rotateToHelicityBasis(double alpha,
						       double beta,
						       double gamma) const{


  EvtDiracSpinor spplus;
  EvtDiracSpinor spminus;
      
  double sqmt2=sqrt(2.*(getP4().mass()));
      
  if (EvtPDL::getStdHep(getId())>0){  
    spplus.set(1.0,0.0,0.0,0.0);
    spminus.set(0.0,1.0,0.0,0.0);
  } else {
    spplus.set(0.0,0.0,0.0,1.0);
    spminus.set(0.0,0.0,1.0,0.0);
  }
      
  spplus.applyRotateEuler(alpha,beta,gamma);
  spminus.applyRotateEuler(alpha,beta,gamma);

  EvtSpinDensity R;
  R.setDim(2);
      
  for (int i=0; i<2; i++) {
    if (EvtPDL::getStdHep(getId())>0){  
      R.set(0,i,(spplus*_spinorRest[i])/sqmt2);
      R.set(1,i,(spminus*_spinorRest[i])/sqmt2);
    } else {
      R.set(0,i,(_spinorRest[i]*spplus)/sqmt2);
      R.set(1,i,(_spinorRest[i]*spminus)/sqmt2);
    } 
  }

  return R;

}



