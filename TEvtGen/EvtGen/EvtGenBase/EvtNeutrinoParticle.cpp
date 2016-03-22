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
// Module: EvtNeutrinoParticle.cc
//
// Description: Class to describe neutrinos
//
// Modification history:
//
//    DJL/RYD     September 25, 1996         Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenBase/EvtNeutrinoParticle.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
using std::endl;


EvtNeutrinoParticle::~EvtNeutrinoParticle(){}

EvtNeutrinoParticle::EvtNeutrinoParticle(){

  return;
}

void EvtNeutrinoParticle::init(EvtId part_n,const EvtVector4R& p4){
  
  _validP4=true;
  setp(p4);
  setpart_num(part_n);
   
  double e,px,py,pz;
  e=p4.get(0);
  px=p4.get(1);
  py=p4.get(2);
  pz=p4.get(3);

  if (EvtPDL::getStdHep(part_n)==0){
    report(Severity::Error,"EvtGen") << "Error in EvtNeutrinoParticle::init, part_n="
			   << part_n.getId()<<endl;
  }

  if (EvtPDL::getStdHep(part_n)>0){  
  
    double beta,alpha,p2,norm;
  
    // See Sakurai p. 167-169
    // and Renton p. 126
  
    p2=px*px+py*py+pz*pz;
  
    beta=acos(pz/sqrt(p2));
    alpha=atan2(py,px);
  
    norm=sqrt(2*e);
  
    double cosb,sinb,cosa,sina;
  
    cosb=cos(0.5*beta);
    sinb=sin(0.5*beta);
  
    cosa=cos(0.5*alpha);
    sina=sin(0.5*alpha);
  
    spinor_parent.set(-norm*sinb*EvtComplex(cosa,-sina),
	  	    norm*cosb*EvtComplex(cosa,sina),
		    norm*sinb*EvtComplex(cosa,-sina),
		    -norm*cosb*EvtComplex(cosa,sina));

  }
  else{

    px=-p4.get(1);
    py=-p4.get(2);
    pz=-p4.get(3);
   
    double pn,sqrpn;

    pn=e;
    sqrpn=sqrt(pn-pz);
   
    spinor_parent.set((1.0/sqrpn)*EvtComplex(px,-py),
                      EvtComplex(sqrpn,0.0),
	 	      (-1.0/sqrpn)*EvtComplex(px,-py),
		      -EvtComplex(sqrpn,0.0)); 


  }

  setLifetime();

} 


EvtDiracSpinor EvtNeutrinoParticle::spParentNeutrino() const {
  
  return spinor_parent;
}

EvtDiracSpinor EvtNeutrinoParticle::spNeutrino() const {

  report(Severity::Error,"EvtGen") << "Tried to get neutrino spinor in restframe"; 
  report(Severity::Error,"EvtGen") << "Will terminate execution."; 

  ::abort();

  return spinor_rest;
}


EvtSpinDensity EvtNeutrinoParticle::rotateToHelicityBasis() const{

  report(Severity::Error,"EvtGen") << "rotateToHelicityBasis not implemented for neutrino."; 
  report(Severity::Error,"EvtGen") << "Will terminate execution."; 

  ::abort();

  EvtSpinDensity rho;
  return rho;
  
}

EvtSpinDensity EvtNeutrinoParticle::rotateToHelicityBasis(double,
							  double,
							  double) const{

  report(Severity::Error,"EvtGen") << "rotateToHelicityBasis(alpha,beta,gama) not implemented for neutrino."; 
  report(Severity::Error,"EvtGen") << "Will terminate execution."; 

  ::abort();

  EvtSpinDensity R;
  R.setDiag(1);
      
  return R;

}



