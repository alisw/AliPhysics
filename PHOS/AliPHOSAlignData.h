#ifndef ALIPHOSALIGNDATA_H
#define ALIPHOSALIGNDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  class for PHOS alignment parameters       //
////////////////////////////////////////////////

#include "TNamed.h"

class AliPHOSAlignData: public TNamed {

 public:
  AliPHOSAlignData();
  AliPHOSAlignData(const char* name);
  AliPHOSAlignData(const AliPHOSAlignData &alignda);
  AliPHOSAlignData& operator= (const AliPHOSAlignData &alignda);
  virtual ~AliPHOSAlignData();
  void Reset();
  virtual void Print(Option_t *option = "") const; 

  // Getters
  Int_t   GetNModules() const {return fNModules;}
  Float_t GetModuleCenter(Int_t module, Int_t axis) const {
    return fModuleCenter[module][axis];}
  Float_t GetModuleAngle(Int_t module, Int_t axis, Int_t angle) const {
    return fModuleAngle[module][axis][angle];}

  // Setters
  void SetNModules(Int_t nModules) {fNModules = nModules;}
  void SetModuleCenter(Int_t module, Int_t axis, Float_t coord) {
    fModuleCenter[module][axis] = coord;}
  void SetModuleAngle(Int_t module, Int_t axis, Int_t angle, Float_t value) {
    fModuleAngle[module][axis][angle] = value;}

 protected:
  Int_t   fNModules;             // number of PHOS modules (max=5)
  Float_t fModuleCenter[5][3];   // xyz-position of the module center
  Float_t fModuleAngle[5][3][2]; // polar and azymuth angles for 3 axes of modules

  ClassDef(AliPHOSAlignData,1)    // PHOS Alignment data
};

#endif
