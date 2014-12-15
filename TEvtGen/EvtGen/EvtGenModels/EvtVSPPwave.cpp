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
// Module: EvtVSPPwave.cc
//
// Description: Routine to decay vector-> scalar photon in P-wave
//
// Modification history:
//
//    RYD       November 24, 1996       Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtTensor4C.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenModels/EvtVSPPwave.hh"
#include "EvtGenBase/EvtReport.hh"
#include <string>
#include "EvtGenBase/EvtVector4C.hh"

EvtVSPPwave::~EvtVSPPwave() {}

std::string EvtVSPPwave::getName(){

  return "VSP_PWAVE";     

}


EvtDecayBase* EvtVSPPwave::clone(){

  return new EvtVSPPwave;

}

void EvtVSPPwave::init(){

  // check that there are 0 arguments
  checkNArg(0);
  checkNDaug(2);

  checkSpinParent(EvtSpinType::VECTOR);
  
  checkSpinDaughter(0,EvtSpinType::SCALAR);
  checkSpinDaughter(1,EvtSpinType::PHOTON);

}

void EvtVSPPwave::initProbMax(){

  setProbMax(1);
}

void EvtVSPPwave::decay( EvtParticle *p){

  p->initializePhaseSpace(getNDaug(),getDaugs());

  EvtParticle *gamma;
  gamma = p->getDaug(1);

  double m_p=p->mass();
  EvtVector4R momgamma = gamma->getP4();

  //work in the parent ,p,  rest frame.
  EvtVector4R p4_p;
  p4_p.set(m_p, 0.0, 0.0, 0.0 );

  //  Put phase space results into the daughters.
 
  EvtTensor4C tds;

  double norm=1/(m_p*momgamma.d3mag());

  tds = dual(EvtGenFunctions::directProd(norm*p4_p,momgamma));

  vertex(0,0,(tds.cont1( p->eps(0))).cont(
	 gamma->epsParentPhoton(0).conj() ) );
  vertex(0,1,(tds.cont1( p->eps(0))).cont(
              gamma->epsParentPhoton(1).conj() ) );

  vertex(1,0,(tds.cont1( p->eps(1))).cont(
              gamma->epsParentPhoton(0).conj() ) );
  vertex(1,1,(tds.cont1( p->eps(1))).cont(
              gamma->epsParentPhoton(1).conj() ) );

  vertex(2,0,(tds.cont1( p->eps(2))).cont(
              gamma->epsParentPhoton(0).conj() ) );
  vertex(2,1,(tds.cont1( p->eps(2))).cont(
              gamma->epsParentPhoton(1).conj() ) );

  return ;
}

