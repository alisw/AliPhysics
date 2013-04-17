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
#include "AliMFTConstants.h"

//====================================================================================================================================================

class AliMFTSegmentation : public TObject {

public:
  
  AliMFTSegmentation();
  AliMFTSegmentation(const Char_t *nameGeomFile);
  
  virtual ~AliMFTSegmentation();
  virtual void Clear(const Option_t* /*opt*/);

  THnSparseC* GetDetElem(Int_t detElemID) const;

  Int_t GetDetElemID(Int_t plane, Int_t detElem) const { return fNMaxDetElemPerPlane*plane + detElem; }
    
  Bool_t Hit2PixelID(Double_t xHit, Double_t yHit, Int_t detElemID, Int_t &xPixel, Int_t &yPixel);  

  Double_t GetPixelSizeX(Int_t detElemID) const { THnSparseC *detElem = GetDetElem(detElemID); return detElem->GetAxis(0)->GetBinWidth(1); }
  Double_t GetPixelSizeY(Int_t detElemID) const { THnSparseC *detElem = GetDetElem(detElemID); return detElem->GetAxis(1)->GetBinWidth(1); }
  Double_t GetPixelSizeZ(Int_t detElemID) const { THnSparseC *detElem = GetDetElem(detElemID); return detElem->GetAxis(2)->GetBinWidth(1); }

  Double_t GetPixelCenterX(Int_t detElemID, Int_t iPixel) const { THnSparseC *detElem = GetDetElem(detElemID); return detElem->GetAxis(0)->GetBinCenter(iPixel+1); }
  Double_t GetPixelCenterY(Int_t detElemID, Int_t iPixel) const { THnSparseC *detElem = GetDetElem(detElemID); return detElem->GetAxis(1)->GetBinCenter(iPixel+1); }
  Double_t GetPixelCenterZ(Int_t detElemID, Int_t iPixel) const { THnSparseC *detElem = GetDetElem(detElemID); return -1.*(detElem->GetAxis(2)->GetBinCenter(iPixel+1)); }

  Int_t GetNPlanes() const { return fMFTPlanes->GetEntries(); }

  AliMFTPlane* GetPlane(Int_t iPlane) const { if (iPlane>=0 && iPlane<fMFTPlanes->GetEntries()) return (AliMFTPlane*) fMFTPlanes->At(iPlane); else return NULL; }

  Bool_t DoesPixelExist(Int_t detElemID, Int_t xPixel, Int_t yPixel);
 
protected:

  static const Int_t fNMaxPlanes          = AliMFTConstants::fNMaxPlanes;                // max number of MFT planes
  static const Int_t fNMaxDetElemPerPlane = AliMFTConstants::fNMaxDetElemPerPlane;

  TClonesArray *fMFTPlanes;

private:

  AliMFTSegmentation(const AliMFTSegmentation &source);
  AliMFTSegmentation& operator=(const AliMFTSegmentation &source);

  ClassDef(AliMFTSegmentation,1)
    
};

//====================================================================================================================================================

#endif

