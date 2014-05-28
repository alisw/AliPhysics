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
// Module: EvtTVP.cc
//
// Description: Routine to implement radiative decay chi_c2 -> psi gamma
//			matrix element from [S.P Baranov et al, PRD 85, 014034 (2012)]
//
// Modification history:
//	AVL	6 July, 2012	Module created
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


#include "EvtGenModels/EvtTVP.hh"


#include <string>
#include <iostream>

using namespace std;



EvtTVP::~EvtTVP() {
//   cout<<"(* AVL EvtTVP::destructor getProbMax(-1) = "<<getProbMax(-1)<<" *)"<<endl;
//   cout<<"(* AVL EvtTVP::destructor "<<ncall<<" calls *)"<<endl;
}

std::string EvtTVP::getName(){
  return "TVP";     
}


EvtDecayBase* EvtTVP::clone(){
//   cout<<" (* AVL: === EvtTVP::clone() ============ *)"<<endl;
  return new EvtTVP;

}

void EvtTVP::decay( EvtParticle *root ){
  ncall++;
//   cout<<" (* AVL  EvtTVP::decay() ============ *)"<<endl;
  double amp2=0;
  root ->initializePhaseSpace(getNDaug(),getDaugs());
  
  EvtVector4R p = root->getDaug(1)->getP4(), // J/psi momentum
    k = root->getDaug(0)->getP4();           // Photon momentum
/*  
    cout<<"(* AVL *) p="<<p<<endl;
    cout<<"(* AVL *) k="<<k<<endl;*/
    
  for(int iPsi = 0; iPsi < 4; iPsi++) {
    for(int iGamma = 0; iGamma < 1; iGamma++) {
      for(int iChi = 0; iChi<4; iChi++) {
	  EvtTensor4C epsChi = root->epsTensor(iChi);
	  EvtVector4C epsPsi = root->getDaug(1)->epsParent(iPsi).conj();
	  EvtVector4C epsGamma = root->getDaug(0)->epsParentPhoton(iGamma).conj();

	  // [Baranov, (11)
	  // matr = p^mu epsPsi^a epsChi_{a b} ( k_mu epsGamma_b  - k_b epsGamma_mu


	  EvtVector4C eee = epsChi.cont1(epsPsi);
	  EvtVector4C vvv = (p*k)*eee - (k*eee)*p;
// 	  cout <<" (* AVL: ginv "<<(vvv*k)<<"  *) "<<endl;
	  EvtComplex amp = vvv*epsGamma;

// 	  cout << "(* AVL *) amp="<<amp<<endl;
	  vertex(iChi, iGamma, iPsi, amp);
	  amp2 = amp2 + abs2(amp);
      };
    };
  };
//   cout <<"(* AVL: amp2 = "<<amp2<<"*)"<<endl;
  
}


void EvtTVP::init(){
//   cout<<" (* AVL: ==== EvtTVP::init() ============ *)"<<endl;

    ncall = 0;
  
  checkNArg(0);
  checkNDaug(2);


  checkSpinParent(EvtSpinType::TENSOR);

  checkSpinDaughter(0,EvtSpinType::PHOTON);
  checkSpinDaughter(1,EvtSpinType::VECTOR);

}

void EvtTVP::initProbMax() {
  setProbMax(1.);
};

