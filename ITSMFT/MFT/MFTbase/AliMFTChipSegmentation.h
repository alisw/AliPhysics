#ifndef AliMFTChipSegmentation_H
#define AliMFTChipSegmentation_H 

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//====================================================================================================================================================
//
//      Class for the description of the structure for the planes of the ALICE Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "AliMFTVSegmentation.h"
#include "AliMFTConstants.h"

//====================================================================================================================================================

class AliMFTChipSegmentation : public AliMFTVSegmentation {

public:

  AliMFTChipSegmentation();
  AliMFTChipSegmentation(UInt_t uniqueID);
  
  virtual ~AliMFTChipSegmentation(){};
  virtual void Clear(const Option_t* /*opt*/) {;}
  virtual void Print(Option_t* /*option*/);

  Bool_t Hit2PixelID(Double_t xHit, Double_t yHit, Int_t &xPixel, Int_t &yPixel);
  
private:
  
  ClassDef(AliMFTChipSegmentation, 1)

};

//====================================================================================================================================================
	
#endif

