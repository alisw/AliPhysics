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
// Module: EvtGen/EvtbTosllScalarAmp.hh
//
// Description:
//
// Modification history:
//
//    RYD     January 5, 2000         Module created
//
//------------------------------------------------------------------------

#ifndef EVTBTOSLLSCALARAMP_HH
#define EVTBTOSLLSCALARAMP_HH

#include "EvtGenModels/EvtbTosllAmp.hh"
class EvtParticle;
class EvtAmp;
class bTosllFF;

class EvtbTosllScalarAmp : public EvtbTosllAmp {

public:

 EvtbTosllScalarAmp(double c7, double c9, double c10):_c7(c7),_c9(c9),_c10(c10){}

  //Daughters are initialized and have been added to the parent.
  //No need to carry around the daughters seperately!
  void CalcAmp( EvtParticle *parent, EvtAmp& amp,
		EvtbTosllFF *formFactors );

private:
  double _c7,_c9,_c10;

};

#endif

