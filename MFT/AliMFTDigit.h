#ifndef AliMFTDigit_H
#define AliMFTDigit_H

/* Copyright(c) 2004-2006, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//====================================================================================================================================================
//
//      Digit description for the ALICE Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "AliDigit.h"

//====================================================================================================================================================

class AliMFTDigit: public AliDigit {

public:

  AliMFTDigit();
//   AliMFTDigit(const AliMFTDigit &d);

  virtual ~AliMFTDigit() {}

  void AddOffset2TrackID(Int_t offset) { for(Int_t iTr=0; iTr<fNMCTracks; iTr++) if (fMCLabel[iTr]>0) fMCLabel[iTr]+=offset; }  // needed for merging (?)
    
  void SetPlane(Int_t plane) { fPlane = plane; }
  void SetDetElemID(Int_t detElemID) { fDetElemID = detElemID; }
  void SetPixID(Int_t pixelX, Int_t pixelY, Int_t pixelZ) { fPixelX = pixelX;  fPixelY = pixelY; fPixelZ = pixelZ; }
  void SetPixCenter(Double_t pixelCenterX, Double_t pixelCenterY, Double_t pixelCenterZ) { 
    fPixelCenterX = pixelCenterX;  
    fPixelCenterY = pixelCenterY; 
    fPixelCenterZ = pixelCenterZ; 
  }
  void SetPixWidth(Double_t pixelWidthX, Double_t pixelWidthY, Double_t pixelWidthZ) { 
    fPixelWidthX = pixelWidthX;  
    fPixelWidthY = pixelWidthY;
    fPixelWidthZ = pixelWidthZ; 
  }
  void SetEloss(Double_t sig) { fEloss = sig; fNElectrons = fEloss/fElossPerElectron; }

  void  AddMCLabel(Int_t label) { if (fNMCTracks>=fNMaxMCTracks) return; else fMCLabel[fNMCTracks++]=label; }
  Int_t GetNMCTracks() const { return fNMCTracks; }
  Int_t GetMCLabel(Int_t track) const { if (track<fNMCTracks && track>=0) return fMCLabel[track]; else return -1; }

  Double_t  GetEloss()         const { return fEloss; }
  Double_t  GetNElectrons()    const { return fNElectrons; }
  Int_t     GetPlane()         const { return fPlane; }
  Int_t     GetDetElemID()     const { return fDetElemID; }
  Int_t     GetPixelX()        const { return fPixelX; }
  Int_t     GetPixelY()        const { return fPixelY; }
  Int_t     GetPixelZ()        const { return fPixelZ; }
  Double_t  GetPixelCenterX()  const { return fPixelCenterX; }
  Double_t  GetPixelCenterY()  const { return fPixelCenterY; }
  Double_t  GetPixelCenterZ()  const { return fPixelCenterZ; }
  Double_t  GetPixelWidthX()   const { return fPixelWidthX; }
  Double_t  GetPixelWidthY()   const { return fPixelWidthY; }
  Double_t  GetPixelWidthZ()   const { return fPixelWidthZ; }
    
protected:
    
  static const Double_t fElossPerElectron = 3.62e-09;
  static const Int_t fNMaxMCTracks = 10;

  Int_t fNMCTracks;

  Int_t fPixelX;
  Int_t fPixelY;  
  Int_t fPixelZ;
  Double_t fPixelCenterX;
  Double_t fPixelCenterY;
  Double_t fPixelCenterZ;
  Double_t fPixelWidthX;
  Double_t fPixelWidthY;
  Double_t fPixelWidthZ;
  Int_t fPlane;     
  Int_t fDetElemID;     
  Double_t fEloss;       // total signal as Eloss in the medium
  Double_t fNElectrons; 
  
  Int_t fMCLabel[fNMaxMCTracks];

  ClassDef(AliMFTDigit,3)

};

//====================================================================================================================================================

#endif



