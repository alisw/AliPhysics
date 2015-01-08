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
// Module: EvtSVP.cc
//
// Description: Routine to implement radiative decay chi_c0 -> psi gamma
//
//
// Modification history:
//	AVL	Jul 6, 2012	modle created
//
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtTensorParticle.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenBase/EvtTensor4C.hh"


#include "EvtGenModels/EvtSVP.hh"


#include <string>
#include <iostream>

using namespace std;



EvtSVP::~EvtSVP() {
  //    cout<<"(* AVL EvtSVP::destructor getProbMax(-1) = "<<getProbMax(-1)<<" *)"<<endl;
  //    cout<<"(* AVL EvtSVP::destructor "<<ncall<<" calls *)"<<endl;

}

std::string EvtSVP::getName(){
  return "SVP";     
}


EvtDecayBase* EvtSVP::clone(){
  //  cout<<" (* AVL: === EvtSVP::clone() ============ *)"<<endl;
  return new EvtSVP;

}

void EvtSVP::decay( EvtParticle *root ){
  //  cout<<"(* AVL EvtSVP::decay getProbMax(-1) = "<<getProbMax(-1)<<" *)"<<endl;
  ncall++;
  //  cout<<" (* AVL  EvtSVP::decay() ============ *)"<<endl;
  root ->initializePhaseSpace(getNDaug(),getDaugs());
  
  EvtVector4R p = root->getDaug(1)->getP4(), // J/psi momentum
    k = root->getDaug(0)->getP4();           // Photon momentum
  for(int iPsi = 0; iPsi < 4; iPsi++) {
    for(int iGamma = 0; iGamma < 1; iGamma++) {
      EvtVector4C epsPsi = root->getDaug(1)->epsParent(iPsi).conj();
      EvtVector4C epsGamma = root->getDaug(0)->epsParentPhoton(iGamma).conj();

      //      EvtComplex amp = epsPsi*epsGamma - (epsPsi*k)/(epsGamma*p)/(k*p);
      EvtComplex amp = (epsPsi*epsGamma) - (epsPsi*k)*(epsGamma*p)/(k*p);
      
      //      cout<<"EvtSVP::decay():  (k*p) = "<<(k*p)<<endl;
      //cout<<"EvtSVP::decay():  amp = "<<amp<<endl;
      
      vertex(iGamma, iPsi, amp);
      };
    };
  
}


void EvtSVP::init(){
  //  cout<<" (* AVL: ==== EvtSVP::init() ============ *)"<<endl;

    ncall = 0;
  
  checkNArg(0);
  checkNDaug(2);


  checkSpinParent(EvtSpinType::SCALAR);

  checkSpinDaughter(0,EvtSpinType::PHOTON);
  checkSpinDaughter(1,EvtSpinType::VECTOR);

}

void EvtSVP::initProbMax() {
    setProbMax(1.2);
};

