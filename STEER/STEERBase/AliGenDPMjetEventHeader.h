#ifndef ALIGENDPMJETEVENTHEADER_H
#define ALIGENDPMJETEVENTHEADER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TLorentzVector.h>

#include "AliGenEventHeader.h"
#include "AliCollisionGeometry.h"

class AliGenDPMjetEventHeader : public AliGenEventHeader, public AliCollisionGeometry
{
 public:
    AliGenDPMjetEventHeader(const char* name);
    AliGenDPMjetEventHeader();
    virtual ~AliGenDPMjetEventHeader() {}

  // Getters
  Float_t TotalEnergy()  const {return fTotalEnergy;} 
  Int_t   Trials()       const {return fTrials;}
  Int_t   ProcessType()  const {return fProcessType;}
  void    SetNDiffractive(Int_t sd1, Int_t sd2, Int_t sdd) {fNSD1 = sd1; fNSD2 = sd2; fNDD = sdd;}
  void    GetNDiffractive(Int_t& sd1, Int_t& sd2, Int_t& sdd) {sd1 = fNSD1; sd2 = fNSD2; sdd = fNDD;}
  // Setters
  void SetTotalEnergy(Float_t energy)  {fTotalEnergy = energy;}
  void SetTrials(Int_t trials)         {fTrials      = trials;}
  void SetProcessType(Int_t type)      {fProcessType = type;}
	  
protected:
  Float_t fTotalEnergy;              // Total energy of produced particles
  Int_t   fTrials;                   // Number of trials to fulfill trigger condition
  Int_t   fProcessType;              // Process Type 
  Int_t   fNSD1;                     // number of SD1 in pA, AA 
  Int_t   fNSD2;                     // number of SD2 in pA, AA
  Int_t   fNDD;                      // number of DD  in pA, AA
  ClassDef(AliGenDPMjetEventHeader,1) // Event header for dpmjet event
};

#endif
