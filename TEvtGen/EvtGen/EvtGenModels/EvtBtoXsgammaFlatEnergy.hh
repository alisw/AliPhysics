//--------------------------------------------------------------------------
//
// Environment:
//      This software is part of the EvtGen package developed jointly
//      for the BaBar and CEO collaborations.  If you use all or part
//      of it, please give an appropriate acknowledgement.
//
// Copyright Information:
//      Copyright (C) 2001      Brunel University, University of Wisconsin
//
// Module: EvtGen/EvtBtoXsgammaFlatEnergy.hh
//
// Description:
//       Implimentation of a fixed hadronic mass to measure spectrum
//
// Modification history:
//
//       Jim Libby     October 11, 2002  Module created
//
//------------------------------------------------------------------------

#ifndef EVTBTOXSGAMMAFLATENERGY_HH
#define EVTBTOXSGAMMAFLATENERGY_HH

#include "EvtGenModels/EvtBtoXsgammaAbsModel.hh"

class EvtBtoXsgammaFlatEnergy : public EvtBtoXsgammaAbsModel {

public:
  
  EvtBtoXsgammaFlatEnergy() {}

  virtual ~EvtBtoXsgammaFlatEnergy();

  void init(int, double*);

  double GetMass(int code);

private:
  //Input parameters
  double _eMin;
  double _eMax;
  //
  double _mB0;
  double _eRange;
};

#endif



