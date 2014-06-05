//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 12001      Caltech
//
// Module: EvtGen/EvtSSDCP.hh
//
// Description: This module is part of the unification of simulation of CP violation in 
//              B decays. This model handles decays of the type B->SD where D is either
//              a spin 0, 1, or 2 particle. See long writeup for more information.
//
// Modification history:
//
//    DJL/RYD     August 12, 2001         Module created
//
//------------------------------------------------------------------------

#ifndef EVTSSDCP_HH
#define EVTSSDCP_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

class EvtSSDCP:public  EvtDecayAmp  {

public:

  EvtSSDCP() {}
  virtual ~EvtSSDCP();
  
  std::string getName();
  EvtDecayBase* clone();

  void initProbMax();
  void init();
  void decay(EvtParticle *p); 

  std::string getParamName(int i);
  std::string getParamDefault(int i);

private:

  //Arguments

  double _dm;

  double _dgog;

  EvtComplex _qoverp;
  EvtComplex _poverq;
  EvtComplex _z;  //FS CPTV parameter

  // FS commented next line becuse not used
  //  int _cp; 

  EvtComplex _A_f;
  EvtComplex _Abar_f;
  
  EvtComplex _A_fbar;
  EvtComplex _Abar_fbar;

  //Derived quantities

  double _gamma;
  double _dgamma;

  bool _eigenstate;

};

#endif
