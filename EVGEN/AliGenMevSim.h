#ifndef ALIGENMEVSIM_H
#define ALIGENMEVSIM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Implementation of the interface for TMevsim
// Author: Sylwester Radomski <radomski@if.pw.edu.pl>


#include "TString.h" 

#include "AliGenerator.h"

#include "AliMevSimConfig.h"
#include "AliMevSimParticle.h"


class AliGenMevSim : public AliGenerator { 

  AliMevSimConfig *fConfig;
 
public:

  AliGenMevSim();
  AliGenMevSim(AliMevSimConfig *config);

  virtual ~AliGenMevSim();

  // 
  virtual void SetConfig(AliMevSimConfig *config);
  virtual void AddParticleType(AliMevSimParticle *type);

  virtual void Init();
  virtual void Generate();

  ClassDef(AliGenMevSim,1) // Interface class for AliMevsim
    
};
#endif
