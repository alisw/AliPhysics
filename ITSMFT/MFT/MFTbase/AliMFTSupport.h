#ifndef AliMFTSupport_H
#define AliMFTSupport_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//=============================================================================================
//
//      Class describing geometry of one MFT half-disk support
//
//      Contact author: raphael.tieulent@cern.ch
//
//=============================================================================================

#include "TNamed.h"
#include "AliLog.h"
#include "TGeoVolume.h"

//=============================================================================================


class AliMFTSupport : public TNamed {
  
public:
  
  AliMFTSupport();
  AliMFTSupport(Double_t zIn, Double_t zOut, Double_t rMin, Double_t rMax, Bool_t isBottom); 
  
  virtual ~AliMFTSupport();
  
  TGeoVolume * CreateVolume();
  
  
protected:
  
  Double_t fZin;
  Double_t fZout;
  Double_t fRmin;
  Double_t fRmax;
  Bool_t   fIsBottom;


private:
  
  ClassDef(AliMFTSupport,1)
  
};

//=============================================================================================

#endif

