//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See EvtGen/COPYRIGHT
//      Copyright (C) 2003      Caltech
//
// Module: EvtGen/EvtbTosllAli.hh
//
// Description:Implementation of the b->sll decays according to Ali '01 et al.
//
// Modification history:
//
//    Ryd     March 30, 2003         Module created
//
//------------------------------------------------------------------------

#ifndef EVTBTOSLLALI_HH
#define EVTBTOSLLALI_HH


#include "EvtGenBase/EvtDecayAmp.hh"

class EvtbTosllFF;
class EvtbTosllAmp;
class EvtParticle;

class EvtbTosllAli:public  EvtDecayAmp  {

public:

  EvtbTosllAli(): _aliffmodel(0), _calcamp(0) {}
  virtual ~EvtbTosllAli();

  std::string getName();
  EvtDecayBase* clone();

  void decay(EvtParticle *p);
  void init();
  void initProbMax();

private:
  EvtbTosllFF *_aliffmodel;
  EvtbTosllAmp *_calcamp;
  double _poleSize;
};

#endif

