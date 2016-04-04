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
// Module: EvtGen/EvtDecayIncoherent.cc
//
// Description:
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------
#include "EvtGenBase/EvtPatches.hh"

#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtDecayIncoherent.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtRadCorr.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtPDL.hh"


void EvtDecayIncoherent::makeDecay(EvtParticle* p, bool recursive){

  //initialize this the hard way..
  //Lange June 26, 2000
  for (size_t i=0; i<static_cast<unsigned int>(MAX_DAUG); i++ ) { 
    spinDensitySet[i]=0;
  }

  _daugsDecayedByParentModel=false;

  decay(p);
  p->setDecayProb(1.0);

  EvtSpinDensity rho;

  rho.setDiag(p->getSpinStates());

  p->setSpinDensityBackward(rho);

  if (getPHOTOS() || EvtRadCorr::alwaysRadCorr()) {
    EvtRadCorr::doRadCorr(p);
  }

  if(!recursive) return;

  //Now decay the daughters.

  if ( !daugsDecayedByParentModel()) {
    
    for(size_t i=0;i<p->getNDaug();i++){
      //Need to set the spin density of the daughters to be
      //diagonal.
      rho.setDiag(p->getDaug(i)->getSpinStates());
      //if (p->getDaug(i)->getNDaug()==0){
      //only do this if the user has not already set the 
      //spin density matrix herself.
      //Lange June 26, 2000
      if ( isDaughterSpinDensitySet(i)==0 ) { 
	p->getDaug(i)->setSpinDensityForward(rho);
      }
      else{
	//report(Severity::Info,"EvtGen") << "spinDensitymatrix already set!!!\n";
	EvtSpinDensity temp=p->getDaug(i)->getSpinDensityForward();
	//	report(Severity::Info,"EvtGen") <<temp<<endl;
      }
      //Now decay the daughter.  Really!
      p->getDaug(i)->decay();
    } 
  }
			    
}







