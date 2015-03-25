#ifndef ALIGENDIME_H
#define ALIGENDIME_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

// Interface to AliRoot of the DIME generator.
// Author: Mikael.Mieskolainen@cern.ch

//- Root Includes
#include <TString.h>
#include <TParticle.h>
#include "TDime.h"

//- AliRoot Includes
#include "AliGenMC.h"

class AliGenDime : public AliGenMC {
public:
  AliGenDime();
  AliGenDime(Int_t npart);

  virtual ~AliGenDime();

  virtual void Init();
  virtual void Generate();

  TDime *GetTDime() {
    return (TDime*) fDMgenerator;
  }
  Bool_t NoDaughters(const TParticle *p) const {
    return (p->GetFirstDaughter() < 0);
  }
  TDime* GetDimeGenerator() const {
    return fDMgenerator;
  }

private:
  AliGenDime(const AliGenDime &p);
  AliGenDime& operator=(const AliGenDime &p);

  TDime* fDMgenerator;  //! Pointer to Dime Generator

  ClassDef(AliGenDime, 1);    // DIME parameterisation generator
};

#endif

