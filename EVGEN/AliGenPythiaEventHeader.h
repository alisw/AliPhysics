#ifndef ALIGENPYTHIAEVENTHEADER_H
#define ALIGENPYTHIAEVENTHEADER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliGenEventHeader.h"


class AliGenPythiaEventHeader : public AliGenEventHeader
{
 public:
  AliGenPythiaEventHeader(const char* name){fNJets = 0;}
  AliGenPythiaEventHeader(){fNJets = 0;}
  virtual ~AliGenPythiaEventHeader() {}
  // Getters
  Int_t ProcessType()  {return fProcessType;}
  // Setters
  void  SetProcessType(Int_t type)  {fProcessType = type;}
  Int_t Trials() {return fTrials;}
  void  SetTrials(Int_t trials) {fTrials = trials;}
  void  AddJet(Float_t px, Float_t py, Float_t pz, Float_t e);
  Int_t NTriggerJets() {return fNJets;}
  void  TriggerJet(Int_t i, Float_t p[4]);
  
	  
protected:
  Int_t   fProcessType;               // PYTHIA process id for this event 
  Int_t   fTrials;                    // Number of trials to fulfill trigger condition
  Int_t   fNJets;                     // Number of triggered jets
  Float_t fJets[4][10];               // Trigger jets   
  ClassDef(AliGenPythiaEventHeader,2) // Event header for Pythia event
};

#endif
