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
// Module: EvtGen/EvtDDalitz.hh
//
// Description:
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTDDALITZ_HH
#define EVTDDALITZ_HH

#include "EvtGenBase/EvtDecayAmp.hh"
#include "EvtGenBase/EvtFlatte.hh"
#include <vector>

class EvtParticle;

class EvtDDalitz:public  EvtDecayAmp  {

public:

  EvtDDalitz() {}
  virtual ~EvtDDalitz();

  std::string getName();
  EvtDecayBase* clone();

  void init();
  void initProbMax();

  void decay(EvtParticle *p); 

private:

  int _d1,_d2,_d3,_flag;

  EvtComplex amplDtoK0PiPi(EvtVector4R p4_p,  EvtVector4R moms1, 
  			   EvtVector4R moms2, EvtVector4R moms3);
  EvtComplex amplDtoK0KK(EvtVector4R p4_p,  EvtVector4R moms1, 
  			 EvtVector4R moms2, EvtVector4R moms3);
  
  vector<EvtFlatteParam> _kkpi_params;

};

#endif
