//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2001        Royal Holloway, University of London
//
// Module: EvtGen/EvtSVVNONCPEIGEN.hh
//
// Description:  See EvtSVVNONCPEIGEN.cc
//
// Modification history:
//
//    Ajit Kurup   9 March 2001      Module created (from EvtSVSNONCPEIGEN.hh)
//
//------------------------------------------------------------------------

#ifndef EVTSVVNONCPEIGEN_HH
#define EVTSVVNONCPEIGEN_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;

class EvtSVVNONCPEIGEN:public  EvtDecayAmp  {

public:

  EvtSVVNONCPEIGEN() {}
  virtual ~EvtSVVNONCPEIGEN();

  std::string getName();
  EvtDecayBase* clone();

  void initProbMax();
  void init();

  void decay(EvtParticle *p); 

  std::string getParamName(int i);
  std::string getParamDefault(int i);

private:

  EvtComplex _A_f[12];
};

#endif
