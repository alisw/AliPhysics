#ifndef ALIEMCALALIGNDATA_H
#define ALIEMCALALIGNDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  class for EMCAL alignment parameters       //
////////////////////////////////////////////////

#include "TNamed.h"
#include "AliEMCAL.h"

class AliEMCALAlignData: public TNamed {

 public:
  AliEMCALAlignData();
  AliEMCALAlignData(const char* name);
  AliEMCALAlignData(const AliEMCALAlignData &alignda);
  AliEMCALAlignData& operator= (const AliEMCALAlignData &alignda);
  virtual ~AliEMCALAlignData();
  void Reset();
  virtual void Print(Option_t *option = "") const; 

  // Getters
  Int_t   GetNSuperModules() const {return fNSuperModules;}
  Float_t GetSuperModuleCenter(Int_t module, Int_t axis) const {
    return fSuperModuleCenter[module][axis];}
  Float_t GetSuperModuleAngle(Int_t module, Int_t axis, Int_t angle) const {
    return fSuperModuleAngle[module][axis][angle];}

  // Setters
  void SetNSuperModules(Int_t nSuperModules) {fNSuperModules = nSuperModules;}
  void SetSuperModuleCenter(Int_t module, Int_t axis, Float_t coord) {
    fSuperModuleCenter[module][axis] = coord;}
  void SetSuperModuleAngle(Int_t module, Int_t axis, Int_t angle, Float_t value) {
    fSuperModuleAngle[module][axis][angle] = value;}

 protected:
  Int_t   fNSuperModules;             // number of EMCAL supermodules (max=12)
  Float_t fSuperModuleCenter[12][3];  // xyz-position of the supermodule center
  Float_t fSuperModuleAngle[12][3][2];// polar and azymuth angles for 3 axes 
                                      // of supermodules

  ClassDef(AliEMCALAlignData,1)    // EMCAL Alignment data
};

#endif
