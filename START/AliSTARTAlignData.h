#ifndef ALISTARTALIGNDATA_H
#define ALISTARTALIGNDATA_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  class for START algnment                  //
////////////////////////////////////////////////

#include "TNamed.h"
#include "AliSTART.h"

class AliSTARTAlignData: public TNamed {

 public:
  AliSTARTAlignData();
  AliSTARTAlignData(const char* name);
  AliSTARTAlignData(const AliSTARTAlignData &alignda);
  AliSTARTAlignData& operator= (const AliSTARTAlignData &alignda);
  virtual ~AliSTARTAlignData();
  void Reset();
  virtual void Print(Option_t *option = "") const; 
  //
  Float_t GetZposition(Int_t i) const {return fSTARTzPosition[i];}
  //
  void SetZposition( Float_t valueC, Float_t valueA) {
    fSTARTzPosition[0]=valueC, fSTARTzPosition[1]=valueA;}

 protected:
  Float_t  fSTARTzPosition[2] ;  // z-position of the two STARTs
  //
  ClassDef(AliSTARTAlignData,1)    // START Alignment data
};

#endif
