#ifndef ALIEMCALALIGNDATA_H
#define ALIEMCALALIGNDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

///////////////////////////////////////////////////////////////////////
//  Class for EMCAL alignment parameters - go to standard tools      //
//  Apply allignment to super modules only                           //
///////////////////////////////////////////////////////////////////////

#include "TNamed.h"

class AliAlignObjMatrix;

class AliEMCALAlignData: public TNamed {

 public:
  AliEMCALAlignData();
  AliEMCALAlignData(const char* name);
  AliEMCALAlignData(const AliEMCALAlignData &alignda);
  AliEMCALAlignData& operator= (const AliEMCALAlignData &alignda);
  virtual ~AliEMCALAlignData();

  void Reset();
  virtual void Print(Option_t *option = "") const; // *MENU*

  // Getters
  Int_t   GetNSuperModules() const {return fNSuperModules;}
  AliAlignObjMatrix *GetSuperModuleMatrix(Int_t module) const
  {
    if(module>=0&&module<fNSuperModules) return fSuperModuleMatrix[module];
    else                                 return 0;
  }

  // Setters
  void SetNSuperModules(Int_t nSuperModules) {fNSuperModules = nSuperModules;}
  void SetSuperModuleMatrix(Int_t module, AliAlignObjMatrix *matrix) 
  {
    if(module>=0&&module<fNSuperModules) fSuperModuleMatrix[module] = matrix;
  }

 protected:
  Int_t   fNSuperModules;                    // number of EMCAL supermodules (max=12)
  AliAlignObjMatrix *fSuperModuleMatrix[12]; //matrix info for supermodules

  ClassDef(AliEMCALAlignData,1)    // EMCAL Alignment data
};

#endif
