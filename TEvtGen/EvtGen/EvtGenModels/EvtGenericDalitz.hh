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
// Module: EvtGen/EvtGenericDalitz.hh
//
// Description: Model to describe a generic dalitz decay
//
// Modification history:
//
//    DCC     16 December, 2011         Module created
//
//------------------------------------------------------------------------

#ifndef EVTGENERICDALITZ_HH
#define EVTGENERICDALITZ_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtFlatte.hh"
#include "EvtGenBase/EvtDalitzReso.hh"
#include <string>
#include <vector>

class EvtParticle;

class EvtGenericDalitz:public  EvtDecayAmp  {

public:

  EvtGenericDalitz() {}
  virtual ~EvtGenericDalitz() {}

  std::string getName();
  EvtDecayBase* clone();

  void init();
  void initProbMax() {};//prob max will be set in init

  void decay(EvtParticle *p); 

  std::string getParamName(int i);

private:

  int _d1,_d2,_d3;
  std::vector<std::pair<EvtComplex,EvtDalitzReso> > _resonances;
};

#endif
