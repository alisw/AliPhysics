//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2002    INFN-Pisa
//
// Module: EvtGen/EvtVSSBMixCPT.hh
//
// Description:
//    Routine to decay vector-> scalar scalar with coherent BB-like mixing
//                              including CPT effects
//    Based on VSSBMIX
//
// Modification history:
//
//    F. Sandrelli, Fernando M-V March 03, 2002 
//
//------------------------------------------------------------------------

#ifndef EVTVSSBMIXCPT_HH
#define EVTVSSBMIXCPT_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtParticle.hh"
#include "EvtGenBase/EvtComplex.hh"

class EvtVSSBMixCPT : public EvtDecayAmp  {
public:
  EvtVSSBMixCPT() {}
  virtual ~EvtVSSBMixCPT();

  std::string getName();
  EvtDecayBase* clone();

  void decay(EvtParticle *p); 
  void init();
  void initProbMax();

  int nRealDaughters() {return 2;}
private:
  double _freq;   // mixing frequency in hbar/mm
  double _dGamma;
  EvtComplex _qoverp;
  EvtComplex _poverq;
  EvtComplex _z; 
  double _chib0_b0bar;
  double _chib0bar_b0;

  EvtComplex _A_f;
  EvtComplex _Abar_f;
  
  EvtComplex _A_fbar;
  EvtComplex _Abar_fbar;

  std::string getParamName(int i);
  std::string getParamDefault(int i);
};

#endif
