#ifndef ALIGENMEVSIM_H
#define ALIGENMEVSIM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Implementation of the interface for TMevsim
// Author: Sylwester Radomski <radomski@if.pw.edu.pl>


#include "AliGenerator.h"

class AliMevSimConfig;
class AliMevSimParticle;

class TMevSim;

class AliGenMevSim : public AliGenerator { 

public:
  AliGenMevSim();
  AliGenMevSim(AliMevSimConfig *config);

  virtual ~AliGenMevSim();

  // 
  virtual void SetConfig(AliMevSimConfig *config);
  virtual void AddParticleType(AliMevSimParticle *type);

  virtual void Init();
  virtual void Generate();

protected:
  TMevSim * fMevSim;           // Pointer to the MevSim generator
  AliMevSimConfig *fConfig;    // Pointer to the Configuration class

 private:
  AliGenMevSim(const AliGenMevSim & gen);
  AliGenMevSim & operator=(const AliGenMevSim & gen);
      
  ClassDef(AliGenMevSim,1) // Interface class for AliMevsim
    
};
#endif
