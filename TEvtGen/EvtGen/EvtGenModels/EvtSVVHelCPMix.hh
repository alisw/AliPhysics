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
// Module: EvtGen/EvtSVVHelCPMix.hh
//
// Description:
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module (SVV_HELAMP) created
//    CATMORE	  March 2004	  	  Amendments made t
//------------------------------------------------------------------------

#ifndef EVTSVVHELCPMIX_HH
#define EVTSVVHELCPMIX_HH

#include "EvtGenBase/EvtDecayAmp.hh"

//Class to handle decays of the form SCALAR -> VECTOR VECTOR
//according the the helicity amplitudes specified by the
//user.  There are 6 arguements, orders as amplitude then
//phase for H+, H0, and H-, in that order.

class EvtAmp;
class EvtParticle;
class EvtId;

class EvtSVVHelCPMix:public  EvtDecayAmp  {

public:

  EvtSVVHelCPMix() {}
  virtual ~EvtSVVHelCPMix();

  std::string getName();
  EvtDecayBase* clone();

  void init();

  EvtComplex hp;
  EvtComplex h0;
  EvtComplex hm;
  double averageM;
  double deltaM;
  double gamma;
  double deltagamma;
  EvtComplex strongphase1;
  EvtComplex strongphase2;
  EvtComplex weakmixingphase;
  EvtComplex weakdirectphase;

  void initProbMax();

  void decay(EvtParticle *p); 

  std::string getParamName(int i);
  std::string getParamDefault(int i);
};

#endif
