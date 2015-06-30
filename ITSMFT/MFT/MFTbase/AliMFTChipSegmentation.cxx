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
//      Class for the description of the virtual segmentation of the chips of the ALICE Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "AliMFTConstants.h"
#include "AliMFTChipSegmentation.h"
#include "AliMFTGeometry.h"
#include "AliLog.h"

ClassImp(AliMFTChipSegmentation)

//====================================================================================================================================================

AliMFTChipSegmentation::AliMFTChipSegmentation():
  AliMFTVSegmentation()
{

  // Default constructor

}

//====================================================================================================================================================

AliMFTChipSegmentation::AliMFTChipSegmentation(UInt_t uniqueID):
  AliMFTVSegmentation()
{
  // constructor
  AliMFTGeometry * mftGeom = AliMFTGeometry::Instance();

  SetUniqueID(uniqueID);

  SetName(Form("MFT_S_%d_%d_%d_%d",
               mftGeom->GetHalfMFTID(GetUniqueID()),
               mftGeom->GetHalfDiskID(GetUniqueID()),
               mftGeom->GetLadderID(GetUniqueID()),
               mftGeom->GetSensorID(GetUniqueID()) ));

  Double_t pos[3];
  pos[0] = mftGeom->GetSensorID(GetUniqueID())*(AliMFTConstants::kSensorLength + AliMFTConstants::kSensorInterspace)
                 + AliMFTConstants::kSensorSideOffset;
  pos[1] = AliMFTConstants::kSensorTopOffset ;
  pos[2] = AliMFTConstants::kFlexThickness ;
  SetPosition(pos);
  
  AliDebug(2,Form("Creating %s, UniqueID = %d, Position = (%.2f, %.2f, %.2f)",GetName(), GetUniqueID(), pos[0], pos[1], pos[2]));

  
}


//====================================================================================================================================================

Bool_t AliMFTChipSegmentation::Hit2PixelID(Double_t xHit, Double_t yHit, Int_t &xPixel, Int_t &yPixel) {
  // TODO Need to work on the Misalignment
  
  Double_t xHitLocal = xHit-AliMFTConstants::kSensorMargin;
  Double_t yHitLocal = yHit-(AliMFTConstants::kSensorMargin + AliMFTConstants::kSensorHeight - AliMFTConstants::kSensorActiveHeight);

//  Double_t xHitLocal = (xHit-(fActiveOrigin[0]+fMisalignmentShift[0]))*fSignLength[0];
//  Double_t yHitLocal = (yHit-(fActiveOrigin[1]+fMisalignmentShift[1]))*fSignLength[1];
  AliDebug(2,Form("Hit %f %f --> Pixel Pitch %f  %f ",xHitLocal,yHitLocal,AliMFTConstants::kXPixelPitch,AliMFTConstants::kYPixelPitch));

  if (xHitLocal<0. || xHitLocal>AliMFTConstants::kSensorActiveWidth) return kFALSE;
  if (yHitLocal<0. || yHitLocal>AliMFTConstants::kSensorActiveHeight) return kFALSE;

  xPixel = Int_t( xHitLocal / AliMFTConstants::kXPixelPitch );
  yPixel = Int_t( yHitLocal / AliMFTConstants::kYPixelPitch );
  AliDebug(1,Form("--> Pixel %d ; %d ",xPixel,yPixel));

  return kTRUE;

}

//==================================================================================================================
void AliMFTChipSegmentation::Print(Option_t* /*option*/){
  
  AliInfo(Form("Sensor %s (Unique ID = %d)",GetName(),GetUniqueID()));
  GetTransformation()->Print();
}
