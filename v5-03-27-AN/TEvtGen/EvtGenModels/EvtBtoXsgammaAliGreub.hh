//--------------------------------------------------------------------------
//
//
// Copyright Information: See EvtGen/COPYRIGHT
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CLEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
//
// Module: EvtGen/EvtBtoXsgammaAliGreubFcn.hh
//
// Description:
// Class to generate non-resonant two-body b->s,gamma decays.
//
// Modification history:
//
//    Mark Ian Williams       July 20, 2000       Module created
//    Jane Tinslay            March 21, 2000      Separated from EvtBtoXsgamma 
//                                                class to allow choice of input models.
//
//------------------------------------------------------------------------

#ifndef EVTBTOXSGAMMAALIGREUB_HH
#define EVTBTOXSGAMMAALIGREUB_HH
#include "EvtGenModels/EvtBtoXsgammaAbsModel.hh"

class EvtBtoXsgammaAliGreub : public EvtBtoXsgammaAbsModel {

public:
  
  EvtBtoXsgammaAliGreub() {}

  virtual ~EvtBtoXsgammaAliGreub();

  void init(int, double*);

  double GetMass(int code);

private:

};

#endif

