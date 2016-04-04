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
// Module: EvtGen/EvtDecayProb.cc
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
#include "EvtGenBase/EvtDecayProb.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtRadCorr.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
using std::endl;

void EvtDecayProb::makeDecay(EvtParticle* p, bool recursive){

  int ntimes=10000;

  double dummy;

  do{
    _weight=1.0;
    _daugsDecayedByParentModel=false;

    decay(p);

    ntimes--;
    
    _prob = _prob/_weight;
    
    dummy=getProbMax(_prob)*EvtRandom::Flat();
    p->setDecayProb(_prob/getProbMax(_prob));

  }while(ntimes&&(_prob<dummy));

  if (ntimes==0){
    report(Severity::Debug,"EvtGen") << "Tried accept/reject:10000"
			   <<" times, and rejected all the times!"<<endl;
    report(Severity::Debug,"EvtGen") << "Is therefore accepting the last event!"<<endl;
    report(Severity::Debug,"EvtGen") << "Decay of particle:"<<
      EvtPDL::name(p->getId()).c_str()<<"(channel:"<<
      p->getChannel()<<") with mass "<<p->mass()<<endl;
    
    for(size_t ii=0;ii<p->getNDaug();ii++){
      report(Severity::Debug,"EvtGen") <<"Daughter "<<ii<<":"<<
	EvtPDL::name(p->getDaug(ii)->getId()).c_str()<<" with mass "<<
	p->getDaug(ii)->mass()<<endl;
    }				   
  }


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
      p->getDaug(i)->setSpinDensityForward(rho);
      
      //Now decay the daughter.  Really!
      p->getDaug(i)->decay();
    } 
  }
			    
}




