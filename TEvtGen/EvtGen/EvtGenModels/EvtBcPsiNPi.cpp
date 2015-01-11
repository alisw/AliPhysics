//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package. If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//
// Module: EvtGenModels/EvtBcBsNPi.hh
//
// Description: Decay model for Bc -> J/psi + npi
//
// Modification history:
//
//    A.Berezhnoy, A.Likhoded, A.Luchinsky  April 2011   Module created
//
//------------------------------------------------------------------------

#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenModels/EvtBcPsiNPi.hh"
#include "EvtGenBase/EvtSpinType.hh"

EvtBcPsiNPi::EvtBcPsiNPi() {

  _beta=-0.108; _mRho=0.775; _gammaRho=0.149;
  _mRhopr=1.364; _gammaRhopr=0.400; _mA1=1.23; _gammaA1=0.4;

  FA0_N=5.9; FA0_c1= 0.049; FA0_c2= 0.0015;
  FAm_N=0.0; FAm_c1=0.0; FAm_c2=0.0;
  FAp_N=-0.074; FAp_c1= 0.049; FAp_c2= 0.0015;
  FV_N=0.11; FV_c1= 0.049; FV_c2= 0.0015; 

}

EvtBcPsiNPi::~EvtBcPsiNPi() {

}

std::string EvtBcPsiNPi::getName() {

  return "BC_PSI_NPI";

}

EvtDecayBase* EvtBcPsiNPi::clone() {

  return new EvtBcPsiNPi;

}

void EvtBcPsiNPi::init() {

  checkNArg(0);

  // check spins
  checkSpinParent(EvtSpinType::SCALAR);
  checkSpinDaughter(0,EvtSpinType::VECTOR);
  // the others are scalar
  for (int i=1; i<=(getNDaug()-1);i++) {
    checkSpinDaughter(i,EvtSpinType::SCALAR);
  }

}

void EvtBcPsiNPi::initProbMax() {

  setProbMax(100.);
  if( getNDaug() == 2 ) {
    setProbMax(330.);
  } else if( getNDaug() == 3 ) {
    setProbMax(11000.); // checked with 30k events
  } else if( getNDaug() == 4 ) {
    setProbMax(36000.);
  }

}
