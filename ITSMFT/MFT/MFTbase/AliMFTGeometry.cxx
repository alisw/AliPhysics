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

const Double_t AliMFTGeometry::kSensorLength=3.; //[cm]
const Double_t AliMFTGeometry::kSensorHeight=1.5; //[cm]
const Double_t AliMFTGeometry::kXPixelPitch=29.250e-4; // 29.15 micron // TODO : Check that
const Double_t AliMFTGeometry::kYPixelPitch=26.880e-4; // 26.88 micron // TODO : Check that
const Double_t AliMFTGeometry::kSensorMargin=29.120e-4; // 29.12 micron // TODO : Check that

const Double_t AliMFTGeometry::kSensorActiveWidth  = kNPixelX * kXPixelPitch; //[cm]
const Double_t AliMFTGeometry::kSensorActiveHeight = kNPixelY * kYPixelPitch; //[cm]

const Double_t AliMFTGeometry::kSensorInterspace = 0.01; //[cm]  Offset between two adjacent chip on a ladder
const Double_t AliMFTGeometry::kSensorSideOffset = 0.04; // [cm] Side Offset between the ladder edge and the chip edge
const Double_t AliMFTGeometry::kSensorTopOffset = 0.04; // [cm] Top Offset between the ladder edge and the chip edge
const Double_t AliMFTGeometry::kLadderOffsetToEnd = 4.7; // [cm] Offset between the last Chip and the end of the ladder toward the DAQ connector
const Double_t AliMFTGeometry::kSensorThickness = 50.e-4; // 50 microns

const Double_t AliMFTGeometry::fHeightActive = 1.3;
const Double_t AliMFTGeometry::fHeightReadout = 0.2;

// Allmost everything you wanted to know about the FPC
const Double_t AliMFTGeometry::kLineWidth= 100.e-4;         // line width, 100 microns
const Double_t AliMFTGeometry::kVarnishThickness= 20.e-4;   // 20 micron
const Double_t AliMFTGeometry::kAluThickness = 25.e-4;      // 25 microns
const Double_t AliMFTGeometry::kKaptonThickness = 75.e-4;   // 75 microns
const Double_t AliMFTGeometry::kFlexThickness = kKaptonThickness + 2*kAluThickness + 2*kVarnishThickness; // total thickness of a FPC
const Double_t AliMFTGeometry::kFlexHeight = 1.68;
const Double_t AliMFTGeometry::kClearance=300.e-4;      // 300 microns clearance without any conducting metal all around the FPC
const Double_t AliMFTGeometry::kRadiusHole1=0.125;      // diameter of the FPC crew, closest to the FPC electric connector
const Double_t AliMFTGeometry::kRadiusHole2=0.1;        // diameter of the FPC pin locator, after the previous hole crew
const Double_t AliMFTGeometry::kHoleShift1=2.8;        // shift of the FPC crew
const Double_t AliMFTGeometry::kHoleShift2=3.6;        // shift of the FPC pin locator
const Double_t AliMFTGeometry::kConnectorOffset=0.4;    // distance between the connector and the start of the FPC
const Double_t AliMFTGeometry::kCapacitorDx=0.05;
const Double_t AliMFTGeometry::kCapacitorDy=0.02;
const Double_t AliMFTGeometry::kCapacitorDz=0.1;
const Double_t AliMFTGeometry::kConnectorLength=0.07; 
const Double_t AliMFTGeometry::kConnectorWidth=0.025;
const Double_t AliMFTGeometry::kConnectorThickness=0.020;
const Double_t AliMFTGeometry::kShiftDDGNDline=0.4; // positionning of the line to separate AVDD/DVDD et AGND/DGND on the FPC
const Double_t AliMFTGeometry::kShiftline=0.025; // positionning of the line along the FPC side
const Double_t AliMFTGeometry::kEpsilon=0.0001; // to see the removed volumes produced by TGeoSubtraction

const Double_t AliMFTGeometry::kGlueThickness=100.e-4; // 100 microns
const Double_t AliMFTGeometry::kGlueEdge=300.e-4; // in case the glue is not spreaded on the whole surface of the sensor


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

