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
/// \class AliMFTChipSegmentation
///
/// Description of the virtual segmentation of the chips
///
// author Raphael Tieulent <raphael.tieulent@cern.ch>
//-----------------------------------------------------------------------------

#include "AliLog.h"

#include "AliMFTConstants.h"
#include "AliMFTChipSegmentation.h"
#include "AliMFTGeometry.h"

/// \cond CLASSIMP
ClassImp(AliMFTChipSegmentation);
/// \endcond

//====================================================================================================================================================
/// Default constructor

AliMFTChipSegmentation::AliMFTChipSegmentation():
  AliMFTVSegmentation()
{
}

//====================================================================================================================================================
/// Constructor
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
  pos[0] = mftGeom->GetSensorID(GetUniqueID())*(AliMFTGeometry::kSensorLength + AliMFTGeometry::kSensorInterspace)
                  + AliMFTGeometry::kSensorSideOffset;
  pos[1] = AliMFTGeometry::kSensorTopOffset ;
  pos[2] = AliMFTGeometry::kFlexThickness ;
  SetPosition(pos);
  
  AliDebug(2,Form("Creating %s, UniqueID = %d, Position = (%.2f, %.2f, %.2f)",GetName(), GetUniqueID(), pos[0], pos[1], pos[2]));
  
}


//====================================================================================================================================================
/// Returns the pixel ID corresponding to a hit at (x,y) in the Sensor  frame
///
/// \param [in] xHit Double_t : x Position of the Hit
/// \param [in] yHit Double_t : y Position of the Hit
///
/// \param [out] xPixel Int_t : x position of the pixel hit on the sensor matrix
/// \param [out] yPixel Int_t : y position of the pixel hit on the sensor matrix
/// \retval <kTRUE> if hit into the active part of the sensor
/// \retval <kFALSE> if hit outside the active part
//

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
  AliDebug(1,Form("--> Hit in Pixel %d ; %d ",xPixel,yPixel));

  return kTRUE;

}

//==================================================================================================================
/// \brief Print out Sensor information (Name, ID, position, orientation)
void AliMFTChipSegmentation::Print(Option_t* /*option*/){
  
  AliInfo(Form("Sensor %s (Unique ID = %d)",GetName(),GetUniqueID()));
  GetTransformation()->Print();
}
