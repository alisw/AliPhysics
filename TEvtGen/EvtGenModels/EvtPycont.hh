//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information: See BelEvtGen/COPYRIGHT
//      Copyright (C) 1998      Caltech, UCSB
//
// Module: BelEvtGen/EvtJscont.hh
//
// Description:
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//    RS          October 28, 2002        copied from JETSET module
//
//------------------------------------------------------------------------

#ifndef EVTPYCONT_HH
#define EVTPYCONT_HH

#include "EvtGenBase/EvtDecayIncoherent.hh"
#include "EvtGenBase/EvtParticle.hh"

class EvtPycont:public  EvtDecayIncoherent  {

public:
  
  EvtPycont() {}
  virtual ~EvtPycont();
  
  std::string getName();
  EvtDecayBase* clone();
  
  void init();
  void initProbMax();
  
  void decay(EvtParticle *p); 
 
};

#endif
