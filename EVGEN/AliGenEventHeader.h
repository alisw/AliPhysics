#ifndef ALIGENEVENTHEADER_H
#define ALIGENEVENTHEADER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TNamed.h>
#include <TDatime.h>

// Event header base class for generator. 
// Stores as a minimum the date, run number, event number, number of particles produced  
// and the impact parameter
// Author: andreas.morsch@cern.ch

class AliGenEventHeader : public TNamed
{
 public:

  AliGenEventHeader(const char* name="Undefined");
  virtual ~AliGenEventHeader() {}
  // Getters
  virtual TDatime Date() {return fDate;}
  virtual Int_t   RunNumber()       {return fRunNumber;}
  virtual Int_t   EventNumber()     {return fEventNumber;}    
  virtual Int_t   NProduced()       {return fNProduced;}
  virtual Float_t ImpactParameter() {return fImpactParameter;}  
  // Setters
  virtual void   SetDate(TDatime date)             {fDate=date;}
  virtual void   SetRunNumber(Int_t run)           {fRunNumber=run;}
  virtual void   SetEventNumber(Int_t event)       {fEventNumber=event;}    
  virtual void   SetNProduced(Int_t nprod)         {fNProduced=nprod;}
  virtual void   SetImpactParameter(Float_t b)     {fImpactParameter=b;}  
  
protected:
  TDatime   fDate;                      // Date
  Int_t     fRunNumber;                 // Run   Number
  Int_t     fEventNumber;               // Event Number
  Int_t     fNProduced;                 // Number stable or undecayed particles
  Float_t   fImpactParameter;           // Impact Parameter

  ClassDef(AliGenEventHeader,1) // Event header for primary event
};

#endif
