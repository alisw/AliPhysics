#ifndef AliMFTSegmentation_H
#define AliMFTSegmentation_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//====================================================================================================================================================
//
//      Segmentation class for the planes of the ALICE Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "TFile.h"
#include "TNtuple.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "AliMFTPlane.h"
#include "TMath.h"

//====================================================================================================================================================

class AliMFTSegmentation : public TObject {

public:
  
  AliMFTSegmentation();
  AliMFTSegmentation(const Char_t *nameGeomFile);

  virtual ~AliMFTSegmentation() {}

  THnSparseC* GetDetElem(Int_t detElemID);

  Int_t GetDetElemID(Int_t plane, Int_t detElem) { return fNMaxDetElemPerPlane*plane + detElem; }
    
  Bool_t Hit2PixelID(Double_t xHit, Double_t yHit, Int_t detElemID, Int_t &xPixel, Int_t &yPixel);  

  Double_t GetPixelSizeX(Int_t detElemID) { THnSparseC *detElem = GetDetElem(detElemID); return detElem->GetAxis(0)->GetBinWidth(1); }
  Double_t GetPixelSizeY(Int_t detElemID) { THnSparseC *detElem = GetDetElem(detElemID); return detElem->GetAxis(1)->GetBinWidth(1); }
  Double_t GetPixelSizeZ(Int_t detElemID) { THnSparseC *detElem = GetDetElem(detElemID); return detElem->GetAxis(2)->GetBinWidth(1); }

  Double_t GetPixelCenterX(Int_t detElemID, Int_t iPixel) { THnSparseC *detElem = GetDetElem(detElemID); return detElem->GetAxis(0)->GetBinCenter(iPixel+1); }
  Double_t GetPixelCenterY(Int_t detElemID, Int_t iPixel) { THnSparseC *detElem = GetDetElem(detElemID); return detElem->GetAxis(1)->GetBinCenter(iPixel+1); }
  Double_t GetPixelCenterZ(Int_t detElemID, Int_t iPixel) { THnSparseC *detElem = GetDetElem(detElemID); return -1.*(detElem->GetAxis(2)->GetBinCenter(iPixel+1)); }

  Int_t GetNPlanes() { return fMFTPlanes->GetEntries(); }

  AliMFTPlane* GetPlane(Int_t iPlane) { if (iPlane>=0 && iPlane<fMFTPlanes->GetEntries()) return (AliMFTPlane*) fMFTPlanes->At(iPlane); else return NULL; }
 
protected:

  static const Int_t fNMaxPlanes = 20;                // max number of MFT planes
  static const Double_t fRadiusMin = 2.225;           // minimum radial distance of the MFT sensors. To be carefully coordinated with fDetElemSuperposition
  static const Double_t fDetElemSuperposition = 0.05; // superposition between bands tasselling the MFT planes, for having a full acceptance coverage
                                                      // even in case of 10 degrees inclined tracks
  static const Double_t fHeightDetElem = 0.5;         // height of the active volume bands composing the planes
  static const Double_t fSupportExtMargin = 0.3;      // minimum border size between the end of the support plane and the sensors

  static const Int_t fNMaxDetElemPerPlane = 1000;

  TClonesArray *fMFTPlanes;

private:

  AliMFTSegmentation(const AliMFTSegmentation &source);
  AliMFTSegmentation& operator=(const AliMFTSegmentation &source);

  ClassDef(AliMFTSegmentation,1)
    
};

//====================================================================================================================================================

#endif

