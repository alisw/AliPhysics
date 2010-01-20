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
// Module: EvtGen/EvtFlatQ2.hh
//
// Description:   B->Xu l nu with flat q2 distribution
//
// Modification history:
//
//    David Cote, U. de Montreal, 11/02/2003    Module created
//
//------------------------------------------------------------------------

#ifndef EVTFLATQ2_HH
#define EVTFLATQ2_HH

#include "EvtGenBase/EvtDecayProb.hh"

class EvtParticle;

class EvtFlatQ2: public  EvtDecayProb  {

public:

  EvtFlatQ2() {}
  virtual ~EvtFlatQ2();

  std::string getName();
  EvtDecayBase* clone();

  void init();
  void initProbMax();

  void decay(EvtParticle *p);
};

#endif

