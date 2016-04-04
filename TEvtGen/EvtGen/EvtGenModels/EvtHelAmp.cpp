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
// Module: EvtHelAmp.cc
//
// Description: Decay model for implementation of generic 2 body
//              decay specified by the helicity amplitudes
//
//
// Modification history:
//
//    RYD       March 14, 1999       Module created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenModels/EvtHelAmp.hh"
#include "EvtGenBase/EvtId.hh"
#include <string>
#include "EvtGenBase/EvtConst.hh"
#include "EvtGenBase/EvtEvalHelAmp.hh"
using std::endl;


EvtHelAmp::~EvtHelAmp() {

  delete _evalHelAmp;

}

std::string EvtHelAmp::getName(){

  return "HELAMP";     

}


EvtDecayBase* EvtHelAmp::clone(){

  return new EvtHelAmp;

}

void EvtHelAmp::init(){

  checkNDaug(2);


  //find out how many states each particle have
  int _nA=EvtSpinType::getSpinStates(EvtPDL::getSpinType(getParentId()));
  int _nB=EvtSpinType::getSpinStates(EvtPDL::getSpinType(getDaug(0)));
  int _nC=EvtSpinType::getSpinStates(EvtPDL::getSpinType(getDaug(1)));

  if (verbose()){
    report(Severity::Info,"EvtGen")<<"_nA,_nB,_nC:"
			 <<_nA<<","<<_nB<<","<<_nC<<endl;
  }

  //find out what 2 times the spin is
  int _JA2=EvtSpinType::getSpin2(EvtPDL::getSpinType(getParentId()));
  int _JB2=EvtSpinType::getSpin2(EvtPDL::getSpinType(getDaug(0)));
  int _JC2=EvtSpinType::getSpin2(EvtPDL::getSpinType(getDaug(1)));

  if (verbose()){
    report(Severity::Info,"EvtGen")<<"_JA2,_JB2,_JC2:"
			 <<_JA2<<","<<_JB2<<","<<_JC2<<endl;
  }

  //allocate memory
  int* _lambdaA2=new int[_nA];
  int* _lambdaB2=new int[_nB];
  int* _lambdaC2=new int[_nC];

  EvtComplexPtr* _HBC=new EvtComplexPtr[_nB];
  for(int ib=0;ib<_nB;ib++){
    _HBC[ib]=new EvtComplex[_nC];
  }

  int i;
  //find the allowed helicities (actually 2*times the helicity!)

  fillHelicity(_lambdaA2,_nA,_JA2,getParentId());
  fillHelicity(_lambdaB2,_nB,_JB2,getDaug(0));
  fillHelicity(_lambdaC2,_nC,_JC2,getDaug(1));

  if (verbose()){
    report(Severity::Info,"EvtGen")<<"Helicity states of particle A:"<<endl;
    for(i=0;i<_nA;i++){
      report(Severity::Info,"EvtGen")<<_lambdaA2[i]<<endl;
    }

    report(Severity::Info,"EvtGen")<<"Helicity states of particle B:"<<endl;
    for(i=0;i<_nB;i++){
      report(Severity::Info,"EvtGen")<<_lambdaB2[i]<<endl;
    }

    report(Severity::Info,"EvtGen")<<"Helicity states of particle C:"<<endl;
    for(i=0;i<_nC;i++){
      report(Severity::Info,"EvtGen")<<_lambdaC2[i]<<endl;
    }
  }

  //now read in the helicity amplitudes

  int argcounter=0;

  for(int ib=0;ib<_nB;ib++){
    for(int ic=0;ic<_nC;ic++){
      _HBC[ib][ic]=0.0;
      if (abs(_lambdaB2[ib]-_lambdaC2[ic])<=_JA2) argcounter+=2;
    }
  }

  checkNArg(argcounter);

  argcounter=0;

  for(int ib=0;ib<_nB;ib++){
    for(int ic=0;ic<_nC;ic++){
      if (abs(_lambdaB2[ib]-_lambdaC2[ic])<=_JA2) {
	_HBC[ib][ic]=getArg(argcounter)*exp(EvtComplex(0.0,getArg(argcounter+1)));;
	argcounter+=2;
	if (verbose()){
	  report(Severity::Info,"EvtGen")<<"_HBC["<<ib<<"]["<<ic<<"]="
			       <<_HBC[ib][ic]<<endl;
	}
      }
    }
  }

  _evalHelAmp=new EvtEvalHelAmp(getParentId(),
				getDaug(0),
				getDaug(1),
				_HBC);

  // Note: these are not class data members but local variables.
  delete [] _lambdaA2;
  delete [] _lambdaB2;
  delete [] _lambdaC2;
  for(int ib=0;ib<_nB;ib++){    
    delete [] _HBC[ib];
  }
  delete [] _HBC;  // _HBC is copied in ctor of EvtEvalHelAmp above.

}


void EvtHelAmp::initProbMax(){

  double maxprob=_evalHelAmp->probMax();

  if (verbose()){
    report(Severity::Info,"EvtGen")<<"Calculated probmax"<<maxprob<<endl;
  }

  setProbMax(maxprob);

}


void EvtHelAmp::decay( EvtParticle *p){

  //first generate simple phase space
  p->initializePhaseSpace(getNDaug(),getDaugs());

  _evalHelAmp->evalAmp(p,_amp2);
    
  return ;

}


void EvtHelAmp::fillHelicity(int* lambda2,int n,int J2, EvtId id){
  
  int i;
  
  //photon is special case!
  if (n==2&&J2==2) {
    lambda2[0]=2;
    lambda2[1]=-2;
    return;
  }

  //and so is the neutrino!
  if (n==1&&J2==1) {
    if (EvtPDL::getStdHep(id)>0){
	//particle i.e. lefthanded
        lambda2[0]=-1;
    }else{
	//anti particle i.e. righthanded
        lambda2[0]=1;
    }
    return;
  }

  assert(n==J2+1);

  for(i=0;i<n;i++){
    lambda2[i]=n-i*2-1;
  }

  return;

}







