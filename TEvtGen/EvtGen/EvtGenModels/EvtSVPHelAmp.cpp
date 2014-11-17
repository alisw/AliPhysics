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
// Module: EvtSVPHelAmp.cc
//
// Description: Routine to decay scalar -> vectors+photon
//              by specifying the helicity amplitudes
//
// Modification history:
//
//    RYD       July 26, 1997       Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtVector4R.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtComplex.hh"
#include "EvtGenModels/EvtSVPHelAmp.hh"
#include <string>

EvtSVPHelAmp::~EvtSVPHelAmp() {}

std::string EvtSVPHelAmp::getName(){

  return "SVP_HELAMP";     

}


EvtDecayBase* EvtSVPHelAmp::clone(){

  return new EvtSVPHelAmp;

}

void EvtSVPHelAmp::initProbMax(){

  setProbMax(getArg(0)*getArg(0)+getArg(2)*getArg(2));

}


void EvtSVPHelAmp::init(){

  // check that there are 4 arguments
  checkNArg(4);
  checkNDaug(2);

  checkSpinParent(EvtSpinType::SCALAR);

  checkSpinDaughter(0,EvtSpinType::VECTOR);
  checkSpinDaughter(1,EvtSpinType::PHOTON);

}

void EvtSVPHelAmp::decay( EvtParticle *p ){

  EvtComplex hp(getArg(0)*cos(getArg(1)),getArg(0)*sin(getArg(1)));
  EvtComplex hm(getArg(2)*cos(getArg(3)),getArg(2)*sin(getArg(3)));


  //  Routine to decay a vector into a vector and scalar.  Started
  //  by ryd on Oct 17, 1996.

      
  // This routine is adopted from EvtSVVHel and since there is
  // a photon that can not have helicity 0 this is put in by 
  // setting the h0 amplitude to 0.
  EvtComplex h0=EvtComplex(0.0,0.0);

  EvtParticle *v1,*ph;

  p->initializePhaseSpace(getNDaug(),getDaugs());
  v1 = p->getDaug(0);
  ph = p->getDaug(1);
  EvtVector4R momv1 = v1->getP4();
  EvtVector4R momph = ph->getP4();

  EvtTensor4C d,g;

  g.setdiag(1.0,-1.0,-1.0,-1.0);

  EvtVector4R v,vp;

  v=momv1/momv1.d3mag();
  vp=(momv1+momph)/(momv1+momph).mass();   

  d=((1.0/sqrt(3.0))*(h0-(hp+hm))*(-1.0/sqrt(3.0)))*g+
    ((1.0/sqrt(2.0))*(hp-hm)*EvtComplex(0.0,1.0)*(sqrt(1.0/2.0)))*dual(EvtGenFunctions::directProd(v,vp))+
    (sqrt(2.0/3.0)*(h0+0.5*(hp+hm))*sqrt(3.0/2.0))*(EvtGenFunctions::directProd(v,v)+(1.0/3.0)*g);

  EvtVector4C ep0,ep1,ep2;  
  
  ep0=d.cont1(v1->eps(0).conj());
  ep1=d.cont1(v1->eps(1).conj());
  ep2=d.cont1(v1->eps(2).conj());

  EvtVector4C ep20,ep21,ep22;

  ep20=ph->epsParentPhoton(0).conj();  
  ep21=ph->epsParentPhoton(1).conj();  

  vertex(0,0,ep0*ep20);
  vertex(0,1,ep0*ep21);
  
  vertex(1,0,ep1*ep20);
  vertex(1,1,ep1*ep21);
   
  vertex(2,0,ep2*ep20);
  vertex(2,1,ep2*ep21);

				 
  return ;

}

