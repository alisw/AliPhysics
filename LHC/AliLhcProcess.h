#ifndef ALILHCPROCESS_H
#define ALILHCPROCESS_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
#include "AliLhcMonitor.h"
#include <TNamed.h>

class AliLHC;

class AliLhcProcess : public TNamed, public AliLhcMonitor
{
 public:
  AliLhcProcess(AliLHC* lhc, const char* name, const char* title);
  virtual void  SetAccelerator(AliLHC* acc) {fAccelerator = acc;}
  AliLhcProcess(const AliLhcProcess &process);
  virtual ~AliLhcProcess();
  virtual void Init(){;}
  virtual void Evolve(Float_t dt);
  virtual void  SetMonitor(Int_t /*n*/) {;}
  virtual void  Record(){;}
  virtual void  DrawPlots(){;}
  AliLhcProcess & operator=(const AliLhcProcess & rhs);
  
 protected:
  AliLHC* fAccelerator;         // Accelerator
  //
  ClassDef(AliLhcProcess,1) // LHC Process Base Class
};

#endif
