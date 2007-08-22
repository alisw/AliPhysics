#ifndef ALIGENBEAMGASNEW_H
#define ALIGENBEAMGASNEW_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliGenCocktail.h"

//
// Class for the simulation of beam gas events with correct timing.
// By default HIJING is used as a generator for pO collisions.
//
// Author: Jochen Klein
//

class AliGenBeamGasNew : public AliGenCocktail
{
 public:
  AliGenBeamGasNew();
  virtual ~AliGenBeamGasNew();

  virtual void Generate();
  void VertexInternal();
  virtual void Init();

  void SetTimeWindow(Float_t twindow);
  bool SetRate(Float_t rate);
  void SetZWindow(Float_t zwindow);

 protected:
  Float_t fItime;   // time of bg-interaction
  Float_t fTwindow; // time-window in which tpc-gate is open
  Float_t fZwindow; // extension of simulation in z-direction in cm
  Float_t fRate;    // rate for bg-interaction in Hz/m

 private:
  AliGenBeamGasNew& operator=(const AliGenBeamGasNew &rhs);
  AliGenBeamGasNew(const AliGenBeamGasNew& rhs);
  ClassDef(AliGenBeamGasNew,1);

};
#endif
