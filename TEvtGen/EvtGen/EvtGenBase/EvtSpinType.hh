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
// Module: EvtGen/EvtSpinType.hh
//
// Description: Class for enumarating the different types of
//              particles and the number of states they have.
//
// Modification history:
//
//    RYD     August 12, 1998         Module created
//
//------------------------------------------------------------------------

#ifndef EVTSPINTYPE_HH
#define EVTSPINTYPE_HH

#include "EvtGenBase/EvtReport.hh"

class EvtSpinType{

public:

  enum spintype { SCALAR,VECTOR,TENSOR,DIRAC,PHOTON,NEUTRINO,STRING,
                  RARITASCHWINGER,SPIN3,SPIN4,SPIN5HALF,SPIN7HALF};

  static int getSpin2(spintype stype);

  static int getSpinStates(spintype stype);

private:

}; 

#endif









