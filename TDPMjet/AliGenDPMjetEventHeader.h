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

  AliGenDPMjetEventHeader(const char* name){;}
  AliGenDPMjetEventHeader(){;}
  virtual ~AliGenDPMjetEventHeader() {}

  // Getters
  Float_t TotalEnergy()  {return fTotalEnergy;} 
  Int_t   Trials() {return fTrials;}

  // Setters
  void SetTotalEnergy(Float_t energy)  {fTotalEnergy=energy;}
  void SetTrials(Int_t trials) {fTrials = trials;}
	  
protected:
  Float_t fTotalEnergy;              // Total energy of produced particles
  Int_t   fTrials;                   // Number of trials to fulfill trigger condition
  
  ClassDef(AliGenDPMjetEventHeader,1) // Event header for dpmjet event
};

#endif
