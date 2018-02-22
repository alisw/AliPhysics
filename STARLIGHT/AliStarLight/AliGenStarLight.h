#ifndef ALIGENSTARLIGHT_H
#define ALIGENSTARLIGHT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: $ */

// Interface to AliRoot of the STARlight generator.
// Author: Christoph.Mayer@cern.ch, Bjorn.Nilsen@cern.ch

//- Root Includes
#include <TString.h>
#include <TParticle.h>
#include "TStarLight.h"

//- AliRoot Includes
#include "AliGenMC.h"

class AliSLEventHeader;

class AliGenStarLight : public AliGenMC {
public:
  AliGenStarLight();
  AliGenStarLight(Int_t npart);

  virtual ~AliGenStarLight();

  const char* GetSVNVersion() const { return TStarLight::GetSVNVersion(); }

  void ImportConfigurationFromFile(const char* filename);
  void SetParameter(const char* line);
  virtual void Init();
  virtual void Generate();

  void SetRapidityMotherRange(Double_t yMin, Double_t yMax) {
    fRapidityMotherMin = yMin;
    fRapidityMotherMax = yMax;
  }
  void SetEtaChildRange(Double_t etaMin, Double_t etaMax) {
    fEtaChildMin = etaMin;
    fEtaChildMax = etaMax;
  }

  TStarLight *GetTStarLight() {
    return (TStarLight*)fSLgenerator;
  }
  Bool_t NoDaughters(const TParticle *p) const {
    return (p->GetFirstDaughter()<0);
  }
  TStarLight* GetStarLightGenerator() const {
    return fSLgenerator;
  }

 private:
  AliGenStarLight(const AliGenStarLight &p);
  AliGenStarLight& operator=(const AliGenStarLight &p);

  Double_t    fRapidityMotherMin; // Max < Min: no cut
  Double_t    fRapidityMotherMax;
  Double_t    fEtaChildMin;       // Max < Min: no cut
  Double_t    fEtaChildMax;
  TStarLight *fSLgenerator;       //! Pointer to StarLight Generator.

  AliSLEventHeader *fHeader;      //!

  ClassDef(AliGenStarLight,7);
} ;

#endif
