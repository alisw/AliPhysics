#ifndef ALIGENPYTHIAEVENTHEADER_H
#define ALIGENPYTHIAEVENTHEADER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliGenEventHeader.h"


class AliGenPythiaEventHeader : public AliGenEventHeader
{
 public:
  AliGenPythiaEventHeader(const char* name){;}
  AliGenPythiaEventHeader(){;}
  virtual ~AliGenPythiaEventHeader() {}
  // Getters
  Int_t ProcessType()  {return fProcessType;}
  // Setters
  void  SetProcessType(Int_t type)  {fProcessType = type;}
  Int_t Trials() {return fTrials;}
  void  SetTrials(Int_t trials) {fTrials = trials;}
protected:
  Int_t   fProcessType;              // PYTHIA process id for this event 
  Int_t   fTrials;                   // Number of trials to fulfill trigger condition
  ClassDef(AliGenPythiaEventHeader,1) // Event header for Pythia event
};

#endif
