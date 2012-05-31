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
// Module: EvtGen/EvtSemiLeptonicScalarAmp.hh
//
// Description:
//
// Modification history:
//
//    DJL/RYD     August 11, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTSEMILEPTONICSCALARAMP_HH
#define EVTSEMILEPTONICSCALARAMP_HH

#include "EvtGenBase/EvtSemiLeptonicAmp.hh"

class EvtAmp;
class EvtParticle;
class EvtSemiLeptonicFF;

class EvtSemiLeptonicScalarAmp : public EvtSemiLeptonicAmp {

  //Daughters are initialized and have been added to the parent.
  //No need to carry around the daughters seperately!
  void CalcAmp( EvtParticle *parent, EvtAmp& amp,
		EvtSemiLeptonicFF *FormFactors );

};

#endif

