#ifndef AliMFTSegmentation_H
#define AliMFTSegmentation_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//====================================================================================================================================================
//
//      Class for the virtual segmentation of the ALICE Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "TNamed.h"
#include "TClonesArray.h"
#include "AliLog.h"
#include "AliMFTVSegmentation.h"
class AliMFTHalfSegmentation;
class AliMFTPlane;


//====================================================================================================================================================

static const Char_t *fNameHalves[2] = { "MFT_Top", "MFT_Bottom" };


class AliMFTSegmentation : public TNamed {

public:
  
  enum THalfMFTType { kBottom, kTop };

  AliMFTSegmentation();
  AliMFTSegmentation(const Char_t *nameGeomFile);
  
  virtual ~AliMFTSegmentation();
  virtual void Clear(const Option_t* /*opt*/);

  void SetHalf(Int_t iHalf) {};
  AliMFTHalfSegmentation* GetHalf(Int_t iHalf) const { if (iHalf==kTop || iHalf==kBottom) return (AliMFTHalfSegmentation*) fMFTHalves->At(iHalf); else return NULL; }

  Int_t GetDetElemLocalID(Int_t half, Int_t disk, Int_t ladder, Int_t sensor) const;
    
  Bool_t Hit2PixelID(Double_t xHit, Double_t yHit, Double_t zHit, Int_t half, Int_t disk, Int_t ladder, Int_t sensor, Int_t &xPixel, Int_t &yPixel);


private:

  TClonesArray *fMFTHalves;

  ClassDef(AliMFTSegmentation,1)
    
};

//====================================================================================================================================================

#endif

