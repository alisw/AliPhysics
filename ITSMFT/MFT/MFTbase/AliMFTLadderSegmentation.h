#ifndef AliMFTLadderSegmentation_H
#define AliMFTLadderSegmentation_H 

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//====================================================================================================================================================
//
//      Class for the description of the virtual segmentation of the ladders of the ALICE Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "AliMFTVSegmentation.h"
#include "TClonesArray.h"
#include "AliMFTChipSegmentation.h"

//====================================================================================================================================================

class AliMFTLadderSegmentation : public AliMFTVSegmentation {

public:

  AliMFTLadderSegmentation();
  AliMFTLadderSegmentation(UInt_t uniqueID);
  AliMFTLadderSegmentation(const AliMFTLadderSegmentation& ladder);
  AliMFTLadderSegmentation& operator=(const AliMFTLadderSegmentation& ladder);

  virtual ~AliMFTLadderSegmentation() { if(fChips){fChips->Delete(); delete fChips; fChips=NULL;} }
  virtual void Print(Option_t* opt="");
  virtual void Clear(const Option_t* /*opt*/) { if(fChips){fChips->Clear();} }
  
  AliMFTChipSegmentation* GetSensor(Int_t sensor);


  void CreateSensors();
  Int_t GetNSensors() const { return fNSensors; };
  Int_t GetNumberOfChips() const { return fNSensors; };
  void SetNSensors(Int_t val) {fNSensors = val;};
  
  Bool_t GetChipActiveOrigin(Int_t chipNumber, Double_t *origin);
  Bool_t GetChipActiveLength(Int_t chipNumber, Double_t *length);
  Bool_t GetChipReadoutOrigin(Int_t chipNumber, Double_t *origin);
  Bool_t GetChipReadoutLength(Int_t chipNumber, Double_t *length);

  AliMFTChipSegmentation* GetChip(Int_t chipNumber);

private:
  
  // measures in cm

  Int_t fNSensors;
  TClonesArray *fChips;

  ClassDef(AliMFTLadderSegmentation, 1)

};

//====================================================================================================================================================
	
#endif

