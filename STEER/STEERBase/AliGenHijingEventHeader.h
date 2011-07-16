#ifndef ALIGENHIJINGEVENTHEADER_H
#define ALIGENHIJINGEVENTHEADER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TLorentzVector.h>

#include "AliGenEventHeader.h"
#include "AliCollisionGeometry.h"

class AliGenHijingEventHeader : public AliGenEventHeader, public AliCollisionGeometry
{
 public:
    AliGenHijingEventHeader(const char* name);
  AliGenHijingEventHeader();
  virtual ~AliGenHijingEventHeader() {}
  // Getters
  Float_t TotalEnergy() const {return fTotalEnergy;} 
  Int_t   Trials() const {return fTrials;}
  
	  
  // Setters
  void SetTotalEnergy(Float_t energy)  {fTotalEnergy=energy;}
  void SetJets(const TLorentzVector* jet1, const TLorentzVector* jet2,
	       const TLorentzVector* jet3, const TLorentzVector* jet4)
      {fJet1 = *jet1; fJet2 = *jet2; fJetFsr1 = *jet3; fJetFsr2 = *jet4;}
  void GetJets(TLorentzVector& jet1, TLorentzVector& jet2,
	       TLorentzVector& jet3, TLorentzVector& jet4) const  
      {jet1 = fJet1; jet2 = fJet2; jet3 = fJetFsr1; jet4 = fJetFsr2;}
  void SetTrials(Int_t trials) {fTrials = trials;}
	  
protected:
  Float_t fTotalEnergy;              // Total energy of produced particles
  Int_t   fTrials;                   // Number of trials to fulfill trigger condition
  
  TLorentzVector  fJet1;             // 4-Momentum-Vector of first   triggered jet  
  TLorentzVector  fJet2;             // 4-Momentum-Vector of second  triggered jet     
  TLorentzVector  fJetFsr1;          // 4-Momentum-Vector of first   triggered jet  
  TLorentzVector  fJetFsr2;          // 4-Momentum-Vector of second  triggered jet     
  
  ClassDef(AliGenHijingEventHeader,5) // Event header for hijing event
};

#endif
