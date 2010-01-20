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
// Module: EvtGen/EvtGoityRoberts.hh
//
// Description:
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTGOITYROBERTS_HH
#define EVTGOITYROBERTS_HH

#include "EvtGenBase/EvtDecayAmp.hh"

class EvtParticle;
class EvtId;

class EvtGoityRoberts:public  EvtDecayAmp  {

public:

  EvtGoityRoberts() {}
  virtual ~EvtGoityRoberts();

  std::string getName();
  EvtDecayBase* clone();
  
  void init();
  void decay(EvtParticle *p); 
  void initProbMax();
  
 private:

  void DecayBDstarpilnuGR(EvtParticle *pb,EvtId ndstar,
			  EvtId nlep, EvtId nnu);

  void DecayBDpilnuGR(EvtParticle *pb,EvtId nd,
		      EvtId nlep, EvtId nnu);



};

#endif
