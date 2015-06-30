#ifndef AliMFTHalfSegmentation_H
#define AliMFTHalfSegmentation_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//====================================================================================================================================================
//
//      Segmentation class for each half of the ALICE Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "TFile.h"
#include "THnSparse.h"
#include "TXMLEngine.h"
#include "TAxis.h"
#include "TVector3.h"
#include "TRotation.h"
#include "AliMFTSegmentation.h"
#include "AliMFTConstants.h"
#include "AliMFTVSegmentation.h"

//====================================================================================================================================================

class AliMFTHalfDiskSegmentation;

class AliMFTHalfSegmentation : public AliMFTVSegmentation {

public:
  
  AliMFTHalfSegmentation();
  AliMFTHalfSegmentation(const Char_t *initFile, const Short_t id);
  AliMFTHalfSegmentation(const AliMFTHalfSegmentation &source);
  void CreateHalfDisks(TXMLEngine* xml, XMLNodePointer_t node);
  XMLNodePointer_t SetHalfPosRot(TXMLEngine* xml, XMLNodePointer_t node);
  void FindHalf(TXMLEngine* xml, XMLNodePointer_t node, XMLNodePointer_t &retnode);

  
  virtual ~AliMFTHalfSegmentation();
  virtual void Clear(const Option_t* /*opt*/);
  
  Bool_t GetID() const {return (GetUniqueID()>>12);};

  THnSparseC* GetDetElem(Int_t detElemID) const;

  Int_t GetDetElemGlobalID(Int_t plane, Int_t detElem) const { return fNMaxDetElemPerPlane*plane + detElem; }
  Int_t GetDetElemLocalID(Int_t detElem) const { return detElem%fNMaxDetElemPerPlane; }
    
  Bool_t Hit2PixelID(Double_t xHit, Double_t yHit, Int_t detElemID, Int_t &xPixel, Int_t &yPixel);  

  Double_t GetPixelSizeX(Int_t detElemID) const { THnSparseC *detElem = GetDetElem(detElemID); return detElem->GetAxis(0)->GetBinWidth(1); }
  Double_t GetPixelSizeY(Int_t detElemID) const { THnSparseC *detElem = GetDetElem(detElemID); return detElem->GetAxis(1)->GetBinWidth(1); }
  Double_t GetPixelSizeZ(Int_t detElemID) const { THnSparseC *detElem = GetDetElem(detElemID); return detElem->GetAxis(2)->GetBinWidth(1); }

  Double_t GetPixelCenterX(Int_t detElemID, Int_t iPixel) const { THnSparseC *detElem = GetDetElem(detElemID); return detElem->GetAxis(0)->GetBinCenter(iPixel+1); }
  Double_t GetPixelCenterY(Int_t detElemID, Int_t iPixel) const { THnSparseC *detElem = GetDetElem(detElemID); return detElem->GetAxis(1)->GetBinCenter(iPixel+1); }
  Double_t GetPixelCenterZ(Int_t detElemID, Int_t iPixel) const { THnSparseC *detElem = GetDetElem(detElemID); return -1.*(detElem->GetAxis(2)->GetBinCenter(iPixel+1)); }

  Int_t GetNHalfDisks() const { return fMFTHalfDisks->GetEntries(); }

  AliMFTHalfDiskSegmentation* GetHalfDisk(Int_t iDisk) const { if (iDisk>=0 && iDisk<fMFTHalfDisks->GetEntries()) return (AliMFTHalfDiskSegmentation*) fMFTHalfDisks->At(iDisk); else return NULL; }

  Bool_t DoesPixelExist(Int_t detElemID, Int_t xPixel, Int_t yPixel);
 
protected:

  static const Int_t fNMaxDisks           = AliMFTConstants::fNMaxPlanes;                // max number of MFT planes
  static const Int_t fNMaxDetElemPerPlane = AliMFTConstants::fNMaxDetElemPerPlane;

  TClonesArray *fMFTHalfDisks;

private:

  ClassDef(AliMFTHalfSegmentation,1)
    
};

//====================================================================================================================================================

#endif

