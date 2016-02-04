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
// $Id$

//-----------------------------------------------------------------------------
/// \class AliMFTLadderSegmentation
///
/// Description of the virtual segmentation of a ladder
///
// author Raphael Tieulent <raphael.tieulent@cern.ch>
//-----------------------------------------------------------------------------

#include "AliLog.h"
#include "AliMFTConstants.h"
#include "AliMFTLadderSegmentation.h"
#include "AliMFTGeometry.h"

/// \cond CLASSIMP
ClassImp(AliMFTLadderSegmentation);
/// \endcond

//====================================================================================================================================================
/// Default constructor

AliMFTLadderSegmentation::AliMFTLadderSegmentation():
  AliMFTVSegmentation(),
  fChips(NULL)
{


}

//====================================================================================================================================================
/// Constructor
/// \param [in] uniqueID UInt_t: Unique ID of the Ladder to build
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
/// Copy Constructor
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
/// Creates the Sensors Segmentation array on the Ladder
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
/// Returns pointer to a sensor segmentation
/// \param [in] sensorID Int_t: ID of the sensor on the ladder
AliMFTChipSegmentation* AliMFTLadderSegmentation::GetSensor(Int_t sensorID) const {
  
  if (sensorID<0 || sensorID>=fNSensors) return NULL;
  
  AliMFTChipSegmentation *chip = (AliMFTChipSegmentation*) fChips->At(sensorID);
  
  return chip;
  
}


//==================================================================================================================
/// Print out Ladder information (position, orientation, # of sensors)
/// \param [in] opt "s" or "sensor" -> The individual sensor information will be printed out as well
void AliMFTLadderSegmentation::Print(Option_t* opt){
  
  AliInfo(Form("Ladder %s (Unique ID = %d)",GetName(),GetUniqueID()));
  GetTransformation()->Print();
  AliInfo(Form("N Sensors = %d",GetNSensors()));
  if(opt && (strstr(opt,"sensor")||strstr(opt,"s"))){
    for (int i=0; i<GetNSensors(); i++)  GetSensor(i)->Print("");

  }
  
}

