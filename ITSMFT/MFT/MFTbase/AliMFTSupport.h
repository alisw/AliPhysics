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
  
  virtual ~AliMFTSupport();
  
  TGeoVolumeAssembly* CreateVolume(Int_t half, Int_t disk);
  TGeoVolumeAssembly* CreatePCBs(Int_t half, Int_t disk);
  TGeoVolume* CreateSupport(Int_t half, Int_t disk);

  
protected:
  
  TGeoVolumeAssembly * fSupportVolume;
  Double_t fSupportThickness;
  Double_t fPCBThickness;

private:
  
  ClassDef(AliMFTSupport,1)
  
};

//=============================================================================================

#endif

