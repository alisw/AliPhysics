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

//____________________________________________________________________
//
#include "AliMFTGeometry.h"
#include "AliMFTConstants.h"
#include "AliMFTGeometryBuilder.h"
#include "AliMFTSegmentation.h"
#include "AliMFTHalfSegmentation.h"
#include "AliMFTHalfDiskSegmentation.h"
#include "AliMFTLadderSegmentation.h"
#include "AliMFTChipSegmentation.h"
#include <TGeoManager.h>
#include <TGeoVolume.h>
#include "TSystem.h"
#include "AliLog.h"

/// \cond CLASSIMP
ClassImp(AliMFTGeometry);
/// \endcond


//____________________________________________________________________
AliMFTGeometry* AliMFTGeometry::fgInstance = 0;

//____________________________________________________________________
AliMFTGeometry*
AliMFTGeometry::Instance()
{
  //
  // singleton access
  //
  // Return:
  //    Singleton
  //
  if (!fgInstance) fgInstance = new AliMFTGeometry("MFT");
  return fgInstance;
}

//____________________________________________________________________
AliMFTGeometry::AliMFTGeometry()
: AliGeometry("MFT", "Muon Forward Tracker"),
fBuilder(NULL),
fSegmentation(NULL),
fSensorVolumeId(0)
{
  // PROTECTED
  //
  // CTOR
  //
  
}

//____________________________________________________________________
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
void AliMFTGeometry::Build()
{
  //
  // Make the geometry.  This delegates to AliMFTGeometryBuilder
  //
  LoadSegmentation();

  if (!fBuilder) fBuilder = new AliMFTGeometryBuilder();
  fBuilder->BuildGeometry();
  delete fBuilder;
  
}
//____________________________________________________________________

void AliMFTGeometry::LoadSegmentation() {
  ///
  /// Loads the MFT Segmentation
  ///
  
  if(!fSegmentation) fSegmentation = new AliMFTSegmentation(gSystem->ExpandPathName("$(ALICE_ROOT)/ITSMFT/MFT/data/AliMFTGeometry.xml" ));
  
  
}
//____________________________________________________________________

UInt_t AliMFTGeometry::GetObjectID(ObjectTypes type, Int_t half, Int_t disk, Int_t ladder, Int_t chip) const{

  UInt_t uniqueID = (type<<13) +  (half<<12) + (disk<<9) + (ladder<<4) + chip;

  return uniqueID;
}
//____________________________________________________________________

Bool_t AliMFTGeometry::Hit2PixelID(Double_t xHit, Double_t yHit, Double_t zHit, Int_t detElemID, Int_t &xPixel, Int_t &yPixel) const{

  return (fSegmentation->Hit2PixelID(xHit, yHit, zHit,
                                     GetHalfMFTID(detElemID), GetHalfDiskID(detElemID), GetLadderID(detElemID), GetSensorID(detElemID),
                                     xPixel,  yPixel));
}

//____________________________________________________________________

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

Int_t AliMFTGeometry::GetDiskNSensors(Int_t diskId) const{

  Int_t nSensors = 0;
  for (int iHalf=0; iHalf<2; iHalf++) {
    AliMFTHalfDiskSegmentation * diskSeg = fSegmentation->GetHalf(iHalf)->GetHalfDisk(diskId);
    if(diskSeg) nSensors += diskSeg->GetNChips();

  }
  return nSensors;
}

//____________________________________________________________________

Int_t AliMFTGeometry::GetDetElemLocalID(Int_t detElemID) const{

  Int_t half, disk, ladder, sensor;
  half   = GetHalfMFTID(detElemID);
  disk   = GetHalfDiskID(detElemID);
  ladder = GetLadderID(detElemID);
  sensor = GetSensorID(detElemID);
  
  return  fSegmentation->GetDetElemLocalID(half, disk, ladder, sensor);
  
}

