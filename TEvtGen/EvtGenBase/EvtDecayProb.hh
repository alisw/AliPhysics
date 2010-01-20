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
// Module: EvtGen/EvtDecayProb.hh
//
// Description:
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EvtDecayProb_HH
#define EvtDecayProb_HH

#include "EvtGenBase/EvtDecayBase.hh"

class EvtParticle;

class EvtDecayProb : public EvtDecayBase{

public:

  void makeDecay(EvtParticle* p, bool recursive=true);

  void setProb(double prob) { _prob=prob;}
  double getProb() { return _prob;}
  inline void setWeight(double weight) {_weight=weight;}

  virtual ~EvtDecayProb() {}


private:

  double _prob;
  double _weight;
  
};



#endif




