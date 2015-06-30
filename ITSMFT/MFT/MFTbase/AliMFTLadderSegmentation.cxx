/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

//====================================================================================================================================================
//
//      Class for the description of the virtual segmentation of the ladders of the ALICE Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "AliLog.h"
#include "AliMFTConstants.h"
#include "AliMFTChipSegmentation.h"
#include "AliMFTLadderSegmentation.h"
#include "AliMFTGeometry.h"

ClassImp(AliMFTLadderSegmentation)

//====================================================================================================================================================

AliMFTLadderSegmentation::AliMFTLadderSegmentation():
  AliMFTVSegmentation(),
  fChips(NULL)
{

  // default constructor

}

//====================================================================================================================================================

AliMFTLadderSegmentation::AliMFTLadderSegmentation(UInt_t uniqueID):
  AliMFTVSegmentation(),
  fChips(NULL)
{
  SetUniqueID(uniqueID);

  AliMFTGeometry * mftGeom = AliMFTGeometry::Instance();
  
  SetName(Form("MFT_L_%d_%d_%d",
               mftGeom->GetHalfMFTID(GetUniqueID()),
               mftGeom->GetHalfDiskID(GetUniqueID()),
               mftGeom->GetLadderID(GetUniqueID()) ));

  // constructor
  
}

//====================================================================================================================================================

AliMFTLadderSegmentation::AliMFTLadderSegmentation(const AliMFTLadderSegmentation& ladder):
  AliMFTVSegmentation(ladder),
  fNSensors(ladder.fNSensors)
{
  // copy constructor
  
  if (ladder.fChips) fChips = new TClonesArray(*(ladder.fChips));
  else       fChips = new TClonesArray("AliMFTChipSegmentation",fNSensors);

  fChips -> SetOwner(kTRUE);

	
}

//====================================================================================================================================================

AliMFTLadderSegmentation& AliMFTLadderSegmentation::operator=(const AliMFTLadderSegmentation& ladder) {
  
  // Assignment operator
  
  // check assignement to self
  if (this != &ladder) {
    
    // base class assignement
    TNamed::operator=(ladder);
    
    // clear memory
    Clear("");
    
    if (ladder.fChips) fChips = new TClonesArray(*(ladder.fChips));
    else fChips = new TClonesArray("AliMFTChipSegmentation",fNSensors);
    fChips -> SetOwner(kTRUE);

  }
  
  return *this;
  
}

//====================================================================================================================================================

void AliMFTLadderSegmentation::CreateSensors() {
  
  if (!fChips) {
    fChips = new TClonesArray("AliMFTChipSegmentation",fNSensors);
    fChips -> SetOwner(kTRUE);
  }

  AliMFTGeometry * mftGeom = AliMFTGeometry::Instance();

  for (Int_t iSensor=0; iSensor<fNSensors; iSensor++) {
    UInt_t sensorUniqueID = mftGeom->GetObjectID(AliMFTGeometry::kSensorType,
                                                 mftGeom->GetHalfMFTID(GetUniqueID()),
                                                 mftGeom->GetHalfDiskID(GetUniqueID()),
                                                 mftGeom->GetLadderID(GetUniqueID()),
                                                 iSensor);
    
    AliMFTChipSegmentation *chip = new AliMFTChipSegmentation(sensorUniqueID);

    new ((*fChips)[iSensor]) AliMFTChipSegmentation(*chip);
    delete chip;
  }

}


//====================================================================================================================================================

AliMFTChipSegmentation* AliMFTLadderSegmentation::GetSensor(Int_t sensorID) {
  
  if (sensorID<0 || sensorID>=fNSensors) return NULL;
  
  AliMFTChipSegmentation *chip = (AliMFTChipSegmentation*) fChips->At(sensorID);
  
  return chip;
  
}

//==================================================================================================================
void AliMFTLadderSegmentation::Print(Option_t* opt){
  
  AliInfo(Form("Ladder %s (Unique ID = %d)",GetName(),GetUniqueID()));
  GetTransformation()->Print();
  AliInfo(Form("N Sensors = %d",GetNSensors()));
  if(opt && (strstr(opt,"sensor")||strstr(opt,"s"))){
    for (int i=0; i<GetNSensors(); i++)  GetSensor(i)->Print("");

  }
  
}

