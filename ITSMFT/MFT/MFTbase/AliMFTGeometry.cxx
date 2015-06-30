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
/* $Id$ */

//-----------------------------------------------------------------------------
/// \class AliMFTGeometry
///
/// Geometry mananger for the MFT
///
// author Raphael Tieulent <raphael.tieulent@cern.ch>
//-----------------------------------------------------------------------------

#include "TSystem.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"

#include "AliMFTGeometry.h"
#include "AliMFTConstants.h"
#include "AliMFTGeometryBuilder.h"
#include "AliMFTSegmentation.h"
#include "AliMFTHalfSegmentation.h"
#include "AliMFTHalfDiskSegmentation.h"
#include "AliMFTLadderSegmentation.h"
#include "AliMFTChipSegmentation.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliMFTGeometry);
/// \endcond


//____________________________________________________________________
AliMFTGeometry* AliMFTGeometry::fgInstance = 0;

//____________________________________________________________________
/// \brief Singleton access
AliMFTGeometry*
AliMFTGeometry::Instance()
{
  if (!fgInstance) fgInstance = new AliMFTGeometry("MFT");
  return fgInstance;
}

//____________________________________________________________________
/// \brief Constructor
AliMFTGeometry::AliMFTGeometry()
: AliGeometry("MFT", "Muon Forward Tracker"),
fBuilder(NULL),
fSegmentation(NULL),
fSensorVolumeId(0)
{
  
}

//____________________________________________________________________
/// \brief Default Constructor
AliMFTGeometry::AliMFTGeometry(const char* )
: AliGeometry("MFT", "Muon Forward Tracker"),
fBuilder(NULL),
fSegmentation(NULL),
fSensorVolumeId(0)
{

}
//____________________________________________________________________

AliMFTGeometry::~AliMFTGeometry() {
  delete fBuilder;
  delete fSegmentation;
}

//____________________________________________________________________
/// \brief Builf both the Virtual segmentation and the real volumes
/// Real part is delegates to AliMFTGeometryBuilder class
void AliMFTGeometry::Build()
{
  LoadSegmentation();

  if (!fBuilder) fBuilder = new AliMFTGeometryBuilder();
  fBuilder->BuildGeometry();
  delete fBuilder;
  
}
//____________________________________________________________________
/// \brief Creates the virtual Segmentation from the XML file
/// $(ALICE_ROOT)/ITSMFT/MFT/data/AliMFTGeometry.xml
void AliMFTGeometry::LoadSegmentation() {

  if(!fSegmentation) fSegmentation = new AliMFTSegmentation(gSystem->ExpandPathName("$(ALICE_ROOT)/ITSMFT/MFT/data/AliMFTGeometry.xml" ));
  
}
//____________________________________________________________________
/// \brief Returns the object Unique ID
/// \param [in] type: Type of the object (see AliMFTGeometry::ObjectTypes)
/// \param [in] half: Half-MFT ID
/// \param [in] disk: Half-Disk ID
/// \param [in] ladder: Ladder ID
/// \param [in] chip: Sensor ID

UInt_t AliMFTGeometry::GetObjectID(ObjectTypes type, Int_t half, Int_t disk, Int_t ladder, Int_t chip) const{

  UInt_t uniqueID = (type<<13) +  (half<<12) + (disk<<9) + (ladder<<4) + chip;

  return uniqueID;
}
//____________________________________________________________________
/// \brief Returns the pixel ID corresponding to a hit at (x,y,z) in the ALICE global frame
///
/// \param [in] xHit Double_t : x Position of the Hit
/// \param [in] yHit Double_t : y Position of the Hit
/// \param [in] zHit Double_t : z Position of the Hit
/// \param [in] detElemID Int_t : Sensor Unique ID in which the hit occured
///
/// \param [out] xPixel Int_t : x position of the pixel hit on the sensor matrix
/// \param [out] yPixel Int_t : y position of the pixel hit on the sensor matrix
/// \retval <kTRUE> if hit into the active part of the sensor
/// \retval <kFALSE> if hit outside the active part

Bool_t AliMFTGeometry::Hit2PixelID(Double_t xHit, Double_t yHit, Double_t zHit, Int_t detElemID, Int_t &xPixel, Int_t &yPixel) const{

  return (fSegmentation->Hit2PixelID(xHit, yHit, zHit,
                                     GetHalfMFTID(detElemID), GetHalfDiskID(detElemID), GetLadderID(detElemID), GetSensorID(detElemID),
                                     xPixel,  yPixel));
}

//____________________________________________________________________
/// \brief Returns the center of the pixel position in the ALICE global frame
///
/// \param [in] xPixel Int_t : x position of the pixel hit on the sensor matrix
/// \param [in] yPixel Int_t : y position of the pixel hit on the sensor matrix
/// \param [in] detElemID Int_t : Sensor Unique ID in which the hit occured
/// \param [out] xCenter,yCenter,zCenter Double_t : (x,y,z) Position of the Hit in ALICE global frame

void AliMFTGeometry::GetPixelCenter(Int_t xPixel, Int_t yPixel, Int_t detElemID, Double_t &xCenter, Double_t &yCenter, Double_t &zCenter ) const{

  Double_t local[3];
  local[0] = (0.5+xPixel) * AliMFTConstants::kXPixelPitch + AliMFTConstants::kSensorMargin;
  local[1] = (0.5+yPixel) * AliMFTConstants::kYPixelPitch + (AliMFTConstants::kSensorHeight-AliMFTConstants::kSensorActiveHeight+ AliMFTConstants::kSensorMargin);
  local[2] = AliMFTConstants::kSensorThickness/2.;

  Double_t master[3];
  
  AliMFTHalfSegmentation * halfSeg = fSegmentation->GetHalf(GetHalfMFTID(detElemID));
  AliMFTHalfDiskSegmentation * diskSeg = halfSeg->GetHalfDisk(GetHalfDiskID(detElemID));
  AliMFTLadderSegmentation * ladderSeg = diskSeg->GetLadder(GetLadderID(detElemID));
  AliMFTChipSegmentation * chipSeg = ladderSeg->GetSensor(GetSensorID(detElemID));

  chipSeg->GetTransformation()->LocalToMaster(local, master);
  for (int i=0; i<3; i++) local[i] = master[i];
  ladderSeg->GetTransformation()->LocalToMaster(local, master);
  for (int i=0; i<3; i++) local[i] = master[i];
  diskSeg->GetTransformation()->LocalToMaster(local, master);
  for (int i=0; i<3; i++) local[i] = master[i];
  halfSeg->GetTransformation()->LocalToMaster(local, master);

  xCenter = master[0];
  yCenter = master[1];
  zCenter = master[2];

}
//____________________________________________________________________
/// \brief Returns the number of sensors on the entire disk (top+bottom)
/// \param [in] diskId Int_t: Disk ID = [0,4]
Int_t AliMFTGeometry::GetDiskNSensors(Int_t diskId) const{

  Int_t nSensors = 0;
  for (int iHalf=0; iHalf<2; iHalf++) {
    AliMFTHalfDiskSegmentation * diskSeg = fSegmentation->GetHalf(iHalf)->GetHalfDisk(diskId);
    if(diskSeg) nSensors += diskSeg->GetNChips();

  }
  return nSensors;
}

//____________________________________________________________________
/// \brief Returns the local ID of the sensor on the disk
/// \param [in] detElemID Int_t: Sensor Unique ID

Int_t AliMFTGeometry::GetDetElemLocalID(Int_t detElemID) const{
  
  return  fSegmentation->GetDetElemLocalID(GetHalfMFTID(detElemID), GetHalfDiskID(detElemID), GetLadderID(detElemID), GetSensorID(detElemID));
  
}

