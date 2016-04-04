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
// Module: EvtPi0Dalitz.cc
//
// Description: pi0 -> e+ e- gamma 
//
// Modification history:
//
//    DJL/RYD    June 30, 1998          Module created
//
//------------------------------------------------------------------------
//
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include <fstream>
#include <stdio.h>
#include <string>
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtEta2MuMuGamma.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtDiracSpinor.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtTensor4C.hh"
using std::fstream;

EvtEta2MuMuGamma::~EvtEta2MuMuGamma() {}

std::string EvtEta2MuMuGamma::getName(){

    return "ETA2MUMUGAMMA";

}

EvtDecayBase* EvtEta2MuMuGamma::clone(){

  return new EvtEta2MuMuGamma;

}


void EvtEta2MuMuGamma::initProbMax(){

  setProbMax(3.5);

}


void EvtEta2MuMuGamma::init(){

  // check that there are 0 arguments
  checkNArg(0);
  checkNDaug(3);

    
  checkSpinParent(EvtSpinType::SCALAR);

  checkSpinDaughter(0,EvtSpinType::DIRAC);
  checkSpinDaughter(1,EvtSpinType::DIRAC);
  checkSpinDaughter(2,EvtSpinType::PHOTON);

}


void EvtEta2MuMuGamma::decay( EvtParticle *p){

  EvtParticle *ep, *em, *gamma;
  setWeight(p->initializePhaseSpace(getNDaug(),getDaugs(),false,
				    0.00000002,0,1));

  ep=p->getDaug(0);
  em=p->getDaug(1);
  gamma=p->getDaug(2);

  // the next four lines generates events with a weight such that
  // the efficiency for selecting them is good. The parameter below of
  // 0.1 is the size of the peak at low q^2 (in arbitrary units).
  // The value of 0.1 is appropriate for muons. 
  // when you use this remember to remove the cut on q^2!
   

  //ep em invariant mass^2
  double m2=(ep->getP4()+em->getP4()).mass2();
  EvtVector4R q=ep->getP4()+em->getP4();
  //Just use the prob summed over spins...

  EvtTensor4C w,v;

  v=2.0*(gamma->getP4()*q)*EvtGenFunctions::directProd(q,gamma->getP4()) 
    - (gamma->getP4()*q)*(gamma->getP4()*q)*EvtTensor4C::g()
    -m2*EvtGenFunctions::directProd(gamma->getP4(),gamma->getP4());
 
  w=4.0*( EvtGenFunctions::directProd(ep->getP4(),em->getP4()) + 
          EvtGenFunctions::directProd(em->getP4(),ep->getP4())
	   -EvtTensor4C::g()*(ep->getP4()*em->getP4()-ep->getP4().mass2()));

  double prob=(real(cont(v,w)))/(m2*m2);
  prob *=(1.0/( (0.768*0.768-m2)*(0.768*0.768-m2)
           +0.768*0.768*0.151*0.151));

  //  report(Severity::Info,"EvtGen") << "prob is "<<prob<<endl;
  setProb(prob);

  return;
}


