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
  Int_t   ProcessType()  {return fProcessType;}
  // Setters
  void SetProcessType(Int_t type)  {fProcessType = type;}
protected:
  Int_t fProcessType;
  ClassDef(AliGenPythiaEventHeader,1) // Event header for Pythia event
};

#endif
