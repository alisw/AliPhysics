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
// Module: EvtGen/EvtDecayAmp.cc
//
// Description: Baseclass for models that calculates amplitudes
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------
#include "EvtGenBase/EvtPatches.hh"



#include "EvtGenBase/EvtDecayBase.hh"
#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtRadCorr.hh"
#include "EvtGenBase/EvtAmp.hh"
#include "EvtGenBase/EvtReport.hh"
using std::endl;


void EvtDecayAmp::makeDecay(EvtParticle* p, bool recursive){

  //original default value
  int ntimes=10000;
  
  int more;

  EvtSpinDensity rho;
  double prob,prob_max;

  _amp2.init(p->getId(),getNDaug(),getDaugs());

  do{

    _daugsDecayedByParentModel=false;
    _weight = 1.0;
    decay(p);

    rho=_amp2.getSpinDensity();

    prob=p->getSpinDensityForward().normalizedProb(rho);

    if (prob<0.0) {
      report(Severity::Error,"EvtGen")<<"Negative prob:"<<p->getId().getId()
			    <<" "<<p->getChannel()<<endl;

      report(Severity::Error,"EvtGen") << "rho_forward:"<<endl;
      report(Severity::Error,"EvtGen") << p->getSpinDensityForward();
      report(Severity::Error,"EvtGen") << "rho decay:"<<endl;
      report(Severity::Error,"EvtGen") << rho <<endl;
    }

    if (prob!=prob) {

      report(Severity::Debug,"EvtGen") << "Forward density matrix:"<<endl;
      report(Severity::Debug,"EvtGen") << p->getSpinDensityForward();

      report(Severity::Debug,"EvtGen") << "Decay density matrix:"<<endl;
      report(Severity::Debug,"EvtGen") << rho;

      report(Severity::Debug,"EvtGen") << "prob:"<<prob<<endl;
      
      report(Severity::Debug,"EvtGen") << "Particle:"
			     <<EvtPDL::name(p->getId()).c_str()<<endl;
      report(Severity::Debug,"EvtGen") << "channel        :"<<p->getChannel()<<endl;
      report(Severity::Debug,"EvtGen") << "Momentum:" << p->getP4() << " " << p->mass() << endl;
      if( p->getParent()!=0){
	report(Severity::Debug,"EvtGen") << "parent:"
			       <<EvtPDL::name(
				p->getParent()->getId()).c_str()<<endl;
	report(Severity::Debug,"EvtGen") << "parent channel        :"
			       <<p->getParent()->getChannel()<<endl;

        size_t i;
	report(Severity::Debug,"EvtGen") << "parent daughters  :";
        for (i=0;i<p->getParent()->getNDaug();i++){
	  report(Severity::Debug,"") << EvtPDL::name(
			    p->getParent()->getDaug(i)->getId()).c_str()
				 << " ";
        }
	report(Severity::Debug,"") << endl;

	report(Severity::Debug,"EvtGen") << "daughters  :";
        for (size_t i=0;i<p->getNDaug();i++){
	  report(Severity::Debug,"") << EvtPDL::name(
			    p->getDaug(i)->getId()).c_str()
				 << " ";
        }
	report(Severity::Debug,"") << endl;

	report(Severity::Debug,"EvtGen") << "daughter momenta  :" << endl;;
        for (size_t i=0;i<p->getNDaug();i++){
	  report(Severity::Debug,"") << p->getDaug(i)->getP4() << " " << p->getDaug(i)->mass();
	  report(Severity::Debug,"") << endl;
        }

      }
    }


    prob/=_weight;

    prob_max = getProbMax(prob);
    p->setDecayProb(prob/prob_max);

    more=prob<EvtRandom::Flat(prob_max);
    
    ntimes--;

  }while(ntimes&&more);

  if (ntimes==0){
    report(Severity::Debug,"EvtGen") << "Tried accept/reject: 10000" 
			   <<" times, and rejected all the times!"<<endl;
   
    report(Severity::Debug,"EvtGen")<<p->getSpinDensityForward()<<endl;
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

  EvtSpinDensity rho_list[10];

  rho_list[0]=p->getSpinDensityForward();

  EvtAmp ampcont;

  if (_amp2._pstates!=1){
    ampcont=_amp2.contract(0,p->getSpinDensityForward());
  }
  else{
    ampcont=_amp2;
  }


  // it may be that the parent decay model has already
  // done the decay - this should be rare and the
  // model better know what it is doing..
  
  if ( !daugsDecayedByParentModel() ){

    if(recursive) {   
    
      for(size_t i=0;i<p->getNDaug();i++){

	rho.setDim(_amp2.dstates[i]);
      
	if (_amp2.dstates[i]==1) {
	  rho.set(0,0,EvtComplex(1.0,0.0));
	}
	else{
	  rho=ampcont.contract(_amp2._dnontrivial[i],_amp2);
	}
	
	if (!rho.check()) {
	  
	  report(Severity::Error,"EvtGen") << "-------start error-------"<<endl;
	  report(Severity::Error,"EvtGen")<<"forward rho failed Check:"<<
	    EvtPDL::name(p->getId()).c_str()<<" "<<p->getChannel()<<" "<<i<<endl;
	  
	  p->printTree();

	  for (size_t idaug = 0; idaug < p->getNDaug(); idaug++) {
	    EvtParticle* daughter = p->getDaug(idaug);
	    if (daughter != 0) {daughter->printTree();}
	  }

	  EvtParticle* pParent = p->getParent();
	  if (pParent != 0) {
	    report(Severity::Error,"EvtGen")<<"Parent:"<<EvtPDL::name(pParent->getId()).c_str()<<endl;

	    EvtParticle* grandParent = pParent->getParent();

	    if (grandParent != 0) {
	      report(Severity::Error,"EvtGen")<<"GrandParent:"<<EvtPDL::name(grandParent->getId()).c_str()<<endl;
	    }
	  }

	  report(Severity::Error,"EvtGen") << " EvtSpinDensity rho: " << rho;
	  
	  _amp2.dump();

	  for(size_t ii=0;ii<i+1;ii++){
	    report(Severity::Error,"EvtGen") << "rho_list[" << ii << "] = " << rho_list[ii];
	  }

	  report(Severity::Error,"EvtGen") << "-------Done with error-------"<<endl;  

	}
      
	p->getDaug(i)->setSpinDensityForward(rho);
	p->getDaug(i)->decay();
	
	rho_list[i+1]=p->getDaug(i)->getSpinDensityBackward();
	
	if (_amp2.dstates[i]!=1){
	ampcont=ampcont.contract(_amp2._dnontrivial[i],rho_list[i+1]);
	}
      
	
      }
      
      p->setSpinDensityBackward(_amp2.getBackwardSpinDensity(rho_list));
      
      
      if (!p->getSpinDensityBackward().check()) {
	
	report(Severity::Error,"EvtGen")<<"rho_backward failed Check"<<
	  p->getId().getId()<<" "<<p->getChannel()<<endl;
      
	report(Severity::Error,"EvtGen") << p->getSpinDensityBackward();
      
      }
    }    
  }


  if (getPHOTOS() || EvtRadCorr::alwaysRadCorr()) {
    int n_daug_orig=p->getNDaug();
    EvtRadCorr::doRadCorr(p);
    int n_daug_new=p->getNDaug();
    for (int i=n_daug_orig;i<n_daug_new;i++){
      p->getDaug(i)->decay();
    }
  }

}










