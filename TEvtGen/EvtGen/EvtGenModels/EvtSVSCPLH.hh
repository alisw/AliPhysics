//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 1999      Caltech, UCSB
//
// Module: EvtGen/EvtSVSCPLH.hh
//
// Description:
//
// Modification history:
//
//    Ryd       March 29, 2001         Module created
//
//------------------------------------------------------------------------

#ifndef EVTSVSCPLH_HH
#define EVTSVSCPLH_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtComplex.hh"

class EvtParticle;

class EvtSVSCPLH:public  EvtDecayAmp  {

public:

  EvtSVSCPLH() {}
  virtual ~EvtSVSCPLH();

  std::string getName();
  EvtDecayBase* clone();

  void initProbMax();
  void init();

  void decay(EvtParticle *p); 


private:

  EvtComplex _Af,_Abarf;
  EvtComplex _qop,_poq;
  
  double _dm;
  double _dgamma;

  


};

#endif




