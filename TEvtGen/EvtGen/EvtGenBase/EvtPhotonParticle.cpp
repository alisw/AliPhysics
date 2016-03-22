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
// Module: EvtPhotonParticle.cc
//
// Description: Class to describe massless vectors
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
#include "EvtGenBase/EvtPhotonParticle.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtReport.hh"
using std::endl;

EvtPhotonParticle::EvtPhotonParticle() {
}

EvtPhotonParticle::~EvtPhotonParticle(){}

void EvtPhotonParticle::init(EvtId part_n,const EvtVector4R& p4){

  init(part_n,p4.get(0),p4.get(1),
       p4.get(2),p4.get(3));
}

void EvtPhotonParticle::init(EvtId part_n,double e,double px,double py,double pz){

  _validP4=true;
  setp(e,px,py,pz);
  setpart_num(part_n);

  setLifetime();
  
  //defere calculation of basis vectors untill they are needed!
  _evalBasis=0;

}


EvtVector4C EvtPhotonParticle::epsParentPhoton(int i){

  if (!_evalBasis){

    _evalBasis=1;
    eps1.set(EvtComplex(0.0,0.0), EvtComplex(-1.0/sqrt(2.0),0.0),
	     EvtComplex(0.0,-1.0/sqrt(2.0)), EvtComplex(0.0,0.0));
    eps2.set(EvtComplex(0.0,0.0), EvtComplex(1.0/sqrt(2.0),0.0),
	     EvtComplex(0.0,-1.0/sqrt(2.0)), EvtComplex(0.0,0.0));
    
    // These are for photon along z axis.  Rotate to get 
    // correct direction...
    
    double phi,theta;

    EvtVector4R p=this->getP4();

    double px=p.get(1);
    double py=p.get(2);
    double pz=p.get(3);

    phi = atan2(py,px);
    theta = acos(pz/sqrt(px*px+py*py+pz*pz));
    eps1.applyRotateEuler(phi,theta,-phi);
    eps2.applyRotateEuler(phi,theta,-phi);
    
  }


  EvtVector4C temp;
  
  switch(i) {

  case 0:
    temp=eps1;
    break;
  case 1:
    temp=eps2;
    break;
  default:
    report(Severity::Error,"EvtGen") << "EvtPhotonParticle.cc: Asked "
			   << "for state:"<<i<<endl;
    ::abort();
    break;
  }

  return temp;
}

EvtVector4C EvtPhotonParticle::epsPhoton(int ){

  report(Severity::Error,"EvtGen") << "EvtPhotonParticle.cc: Can not get "
			 << "state in photons restframe."<<endl;;
  ::abort();
  return EvtVector4C();

}


EvtSpinDensity EvtPhotonParticle::rotateToHelicityBasis() const {

  EvtVector4C eplus(0.0,-1.0/sqrt(2.0),
		    EvtComplex(0.0,-1.0/sqrt(2.0)),0.0);
  EvtVector4C eminus(0.0,1.0/sqrt(2.0),
		     EvtComplex(0.0,-1.0/sqrt(2.0)),0.0);

  //Really uggly have to cast away constness because the
  //function epsParentPhoton caches the state vectors...
  EvtVector4C e1=((EvtParticle*)this)->epsParentPhoton(0);
  EvtVector4C e2=((EvtParticle*)this)->epsParentPhoton(1);


  EvtSpinDensity R;
  R.setDim(2);

  R.set(0,0,(eplus.conj())*e1);
  R.set(0,1,(eplus.conj())*e2);
  
  R.set(1,0,(eminus.conj())*e1);
  R.set(1,1,(eminus.conj())*e2);
  
  return R;
	
}


EvtSpinDensity EvtPhotonParticle::rotateToHelicityBasis(double alpha,
							double beta,
							double gamma) const{

  EvtVector4C eplus(0.0,-1.0/sqrt(2.0),
		    EvtComplex(0.0,-1.0/sqrt(2.0)),0.0);
  EvtVector4C eminus(0.0,1.0/sqrt(2.0),
		     EvtComplex(0.0,-1.0/sqrt(2.0)),0.0);
  
  eplus.applyRotateEuler(alpha,beta,gamma);
  eminus.applyRotateEuler(alpha,beta,gamma);
  
  
  //Really uggly have to cast away constness because the
  //function epsParentPhoton caches the state vectors...
  EvtVector4C e1=((EvtParticle*)this)->epsParentPhoton(0);
  EvtVector4C e2=((EvtParticle*)this)->epsParentPhoton(1);
  
  EvtSpinDensity R;
  R.setDim(2);
  
  R.set(0,0,(eplus.conj())*e1);
  R.set(0,1,(eplus.conj())*e2);
  
  R.set(1,0,(eminus.conj())*e1);
  R.set(1,1,(eminus.conj())*e2);
  
  return R;

}


