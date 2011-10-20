#ifndef AliMFTReconstructor_H
#define AliMFTReconstructor_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//====================================================================================================================================================
//
//      Reconstructor class for the ALICE Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "TObjArray.h"
#include "TTree.h"
#include "AliMFTSegmentation.h"
#include "AliReconstructor.h"
#include "AliMFTClusterFinder.h"

//====================================================================================================================================================

class AliRawReader;

class AliMFTReconstructor: public AliReconstructor { 

public:

  AliMFTReconstructor();
  virtual ~AliMFTReconstructor();
  virtual void Init();

  virtual void ResetDigits(); 
  virtual void ResetDigits(Int_t plane);

  virtual void  Reconstruct(TTree *digitsTree, TTree *clustersTree) const; 
  virtual void  Reconstruct(AliRawReader* /*rawdata*/, TTree* /*clustersTree*/) const { AliInfo("Not implemented"); } 

  //  static const AliMFTRecoParam* GetRecoParam() { return dynamic_cast<const AliMFTRecoParam*>(AliReconstructor::GetRecoParam(0)); }

private:
 
  AliMFTReconstructor(const AliMFTReconstructor&);              // Not implemented
  AliMFTReconstructor &operator=(const AliMFTReconstructor&);   // Not implemented

  static const Int_t fNMaxDigitPerPlane = 10000;

  TObjArray  *fDigits;     
  Int_t      fNPlanes;

  ClassDef(AliMFTReconstructor, 1)        // class for the MFT reconstruction

};

//====================================================================================================================================================

#endif
