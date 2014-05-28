//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2001      Royal Holloway, University of London
//
// Module: EvtSVVNONCPEIGEN.cc
//
// Description: Routine to decay scalar -> vector vector 
//              and has CP violation.
//
//              This model does all the ckm-suppressed decays and mixing for you. It randomly 'overwrites' 
//              any reco or tagging state as set in the Y(4S) decay model (VSS_(B)MIX) with its own generated states.
//
//              As such, the corresponding dec file requires only one decay-mode description, for example:
//              Decay MyB0
//              1.000    rho+ MyD*-       SVV_NONCPEIGEN dm beta gamma 0.322 0.31 0.941 0 0.107 1.42 0.02 0 0.02 0 0.02 0 ;
//              EndDecay
//              and furthermore Y(4S) only needs to decay to B0's (or B0bar's).
//              The decay above should be a CKM-favored mode (eg. B0->D*-rho+ or B0bar->D*+rho-).
//              All ckm-suppressed decays and the mixing are derived from this line in the ::Decay function.
//
//              There are 15 or 27 arguments. The first three are dm, phase1
//              and phase2. dm is the B0-B0bar mass difference. Phases 1
//              and 2 are the CKM weak phases relevant for the particular mode, 
//              eg for B-->DstRho phase1 is beta and phase2 is gamma.
//
//              The next arguments are the 2 amplitudes (= 12 input parameters) 
//              in the order: A_f, Abar_f. In the example above, the 'A_f' amplitude now 
//              stands for the ckm-favored decay 'B0->D*-rho+', and 'Abar_f' stands for 'B0bar->D*-rho+'
//
//              Each amplitude has its 3 helicity states in the order +, 0, -, which are each 
//              specified by a magnitude and a strong phase.
//
//              The last 2 arguments A_fbar and Abar_fbar (=12 input parameters) are not necessary, 
//              but can included if one wants to set them differently from A_f, Abar_f.
//
//              Mind you that Hbar_+- = H_-+ (ignoring the weak phase, which flips sign).
//              It is custumary to select one set of helicity states (eg H_+-) and to adopt these for
//              the CP-conjugate decays as well (ie. depict Hbar_-+ with H_+-), which is the interpretation
//              we use for the input-parameters above. 
//              However, the angular decay in EvtGen is just a formula in which helicity amplitudes are 'plugged' in,
//              making no difference between B0 or B0bar decays. In the model below we (thus) account for the +- 
//              flipping between B0 and B0bar.
//              
//
// Modification history:
//    Ajit Kurup 9 March 2001        Module created (from EvtSVSNONCPEIGEN)
//    Max Baak 01/16/2004            Fix of Helicity amplitude ordering.
//                                   Decay also works for B0bar decays.
//------------------------------------------------------------------------
// 
#include "EvtGenBase/EvtPatches.hh"
#include <stdlib.h>
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtRandom.hh"
#include "EvtGenBase/EvtGenKine.hh"
#include "EvtGenBase/EvtCPUtil.hh"
#include "EvtGenBase/EvtPDL.hh"
#include "EvtGenBase/EvtReport.hh"
#include "EvtGenBase/EvtVector4C.hh"
#include "EvtGenModels/EvtSVVNONCPEIGEN.hh"
#include <string>
#include "EvtGenModels/EvtSVVHelAmp.hh"
#include "EvtGenBase/EvtConst.hh"

EvtSVVNONCPEIGEN::~EvtSVVNONCPEIGEN() {}

std::string EvtSVVNONCPEIGEN::getName(){

  return "SVV_NONCPEIGEN";     

}


EvtDecayBase* EvtSVVNONCPEIGEN::clone(){

  return new EvtSVVNONCPEIGEN;

}

void EvtSVVNONCPEIGEN::init(){

  // check that there are 27 arguments
  checkNArg(27,15);
  checkNDaug(2);

  checkSpinDaughter(0,EvtSpinType::VECTOR);
  checkSpinDaughter(1,EvtSpinType::VECTOR);

  //  The ordering of A_f is :
  //  A_f[0-2] = A_f
  //  A_f[3-5] = Abar_f
  //  A_f[6-8] = A_fbar 
  //  A_f[9-11] = Abar_fbar
  //  
  //  Each of the 4 amplitudes include the 3 different helicity states in 
  //  the order +, 0, -. See more about helicity amplitude ordering in ::decay

  int i=0;
  int j=(getNArg()-3)/2;

  for(i=0; i<j; ++i){
    _A_f[i] = getArg((2*i)+3) * EvtComplex( cos(getArg((2*i)+4)),sin(getArg((2*i)+4)) );
  }

  //  If only 6 amplitudes are specified, calculate the last 6 from the first 6:
  if(6 == j){
    for(i = 0; i < 3; ++i){
      _A_f[6+i] = _A_f[3+i];
      _A_f[9+i] = _A_f[i];
    }
  }
}

void EvtSVVNONCPEIGEN::initProbMax() {
  double probMax = 0;
  for (int i = 0; i < 12; ++i){
    double amp = abs(_A_f[i]);
    probMax += amp * amp;
  }

  setProbMax(probMax); 
}

void EvtSVVNONCPEIGEN::decay( EvtParticle *p){

  //added by Lange Jan4,2000
  static EvtId B0=EvtPDL::getId("B0");
  static EvtId B0B=EvtPDL::getId("anti-B0");

  double t;
  EvtId other_b;
  EvtId daugs[2];


  // MB: flip selects the final of the decay
  int flip = ((p->getId() == B0) ? 0 : 1);
  daugs[0]=getDaug(0);
  daugs[1]=getDaug(1);
  p->initializePhaseSpace(2,daugs);

  EvtCPUtil::getInstance()->OtherB(p,t,other_b,0.5);

  EvtComplex amp[3];

  double dmt2 = getArg(0) * t / (2 * EvtConst::c);
  double phiCKM = (2.0 * getArg(1) + getArg(2));   // 2b+g
  EvtComplex ePlusIPhi(cos(phiCKM), sin(phiCKM));
  EvtComplex eMinusIPhi(cos(-phiCKM), sin(-phiCKM));

  // flip == 0 : D*-rho+
  // flip == 1 : D*+rho-

  if (!flip) {
    if (other_b==B0B){
      // At t=0 we have a B0
      for (int i=0; i<3; ++i) {
	amp[i] = _A_f[i]*cos(dmt2) + eMinusIPhi*EvtComplex(0.0,sin(dmt2))*_A_f[i+3];
      }
    }
    if (other_b==B0){
      // At t=0 we have a B0bar
      for(int i=0; i<3; ++i) {
	amp[i] = _A_f[i]*ePlusIPhi*EvtComplex(0.0,sin(dmt2)) + _A_f[i+3]*cos(dmt2);
      }
    }
  } else{
    if (other_b==B0B){
      // At t=0 we have a B0

      // M.Baak 01/16/2004
      // Note: \bar{H}+- = H-+ 
      // If one wants to use the correct helicities for B0 and B0bar decays but the same formula-notation (as done in EvtSVV_HelAmp), 
      // count the B0bar helicities backwards. (Equivalently, one could flip the chi angle.)

      for(int i=0; i<3; ++i) { 
	amp[i] = _A_f[8-i]*cos(dmt2) + eMinusIPhi*EvtComplex(0.0,sin(dmt2))*_A_f[11-i];
      }
    }
    if (other_b==B0){
      // At t=0 we have a B0bar
      for(int i=0; i<3; ++i) {
	amp[i] = _A_f[8-i] * ePlusIPhi * EvtComplex(0.0,sin(dmt2)) + _A_f[11-i]*cos(dmt2);
      }
    }
  }
  
  EvtSVVHelAmp::SVVHel(p,_amp2,daugs[0],daugs[1],amp[0],amp[1],amp[2]);

  return ;
}

std::string EvtSVVNONCPEIGEN::getParamName(int i) {
  switch(i) {
  case 0:
    return "deltaM";
  case 1:
    return "weakPhase1";
  case 2:
    return "weakPhase2";
  case 3:
    return "AfPlusHelAmp";
  case 4:
    return "AfPlusHelAmpPhase";
  case 5:
    return "AfZeroHelAmp";
  case 6:
    return "AfZeroHelAmpPhase";
  case 7:
    return "AfMinusHelAmp";
  case 8:
    return "AfMinusHelAmpPhase";
  case 9:
    return "AbarfPlusHelAmp";
  case 10:
    return "AbarfPlusHelAmpPhase";
  case 11:
    return "AbarfZeroHelAmp";
  case 12:
    return "AbarfZeroHelAmpPhase";
  case 13:
    return "AbarfMinusHelAmp";
  case 14:
    return "AbarfMinusHelAmpPhase";
  case 15:
    return "AfbarPlusHelAmp";
  case 16:
    return "AfbarPlusHelAmpPhase";
  case 17:
    return "AfbarZeroHelAmp";
  case 18:
    return "AfbarZeroHelAmpPhase";
  case 19:
    return "AfbarMinusHelAmp";
  case 20:
    return "AfbarMinusHelAmpPhase";
  case 21:
    return "AbarfbarPlusHelAmp";
  case 22:
    return "AbarfbarPlusHelAmpPhase";
  case 23:
    return "AbarfbarZeroHelAmp";
  case 24:
    return "AbarfbarZeroHelAmpPhase";
  case 25:
    return "AbarfbarMinusHelAmp";
  case 26:
    return "AbarfbarMinusHelAmpPhase";
  default:
    return "";
  }
}

std::string EvtSVVNONCPEIGEN::getParamDefault(int i) {
  switch(i) {
  case 3:
    return "1.0";
  case 4:
    return "0.0";
  case 5:
    return "1.0";
  case 6:
    return "0.0";
  case 7:
    return "1.0";
  case 8:
    return "0.0";
  case 9:
    return "1.0";
  case 10:
    return "0.0";
  case 11:
    return "1.0";
  case 12:
    return "0.0";
  case 13:
    return "1.0";
  case 14:
    return "0.0";
  default:
    return "";
  }
}
