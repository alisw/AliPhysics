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
// Description: Decay model for Bc -> Bs + npi
//
// Modification history:
//
//    A.Berezhnoy, A.Likhoded, A.Luchinsky  July 2011   Module created
//
//------------------------------------------------------------------------

#include "EvtGenBase/EvtPatches.hh"
#include "EvtGenBase/EvtSpinType.hh"
#include "EvtGenModels/EvtBcBsNPi.hh"

EvtBcBsNPi::EvtBcBsNPi() {

  _beta=-0.108; _mRho=0.775; _gammaRho=0.149;
  _mRhopr=1.364; _gammaRhopr=0.400; _mA1=1.23; _gammaA1=0.4;
  //		Fp_N=1.3; Fp_c1=0.30; Fp_c2=0.069;
  Fp_N=3*1.3; Fp_c1=0.30; Fp_c2=0.069;
  Fm_N=0.0; Fm_c1=0.0; Fm_c2=0.0;

}

EvtBcBsNPi::~EvtBcBsNPi() {
}

std::string EvtBcBsNPi::getName() {

  return "BC_BS_NPI";

}

EvtDecayBase* EvtBcBsNPi::clone() {

  return new EvtBcBsNPi;

}

void EvtBcBsNPi::init() {

  checkNArg(0);

  // check spins
  checkSpinParent(EvtSpinType::SCALAR);
  checkSpinDaughter(0,EvtSpinType::SCALAR);
  // the others are scalar
  for (int i=1; i<=(getNDaug()-1);i++) {
    checkSpinDaughter(i,EvtSpinType::SCALAR);
  }

}

void EvtBcBsNPi::initProbMax() {

  if ( getNDaug() == 2 ) {
    setProbMax(250.);
  } else if ( getNDaug() == 3 ) {
    setProbMax(25000.);// checked at 30k events
  } else if( getNDaug() == 4 ) {
    setProbMax(45000.); // checked at 30k events
  }

}
