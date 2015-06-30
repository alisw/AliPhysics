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

//=============================================================================================
//
//      Class describing geometry of one MFT Ladder
//
//      Contact author: raphael.tieulent@cern.ch
//
//=============================================================================================
#include "TGeoMatrix.h"
#include "TGeoManager.h"
#include "TGeoBBox.h"

#include "AliMFTLadderSegmentation.h"
#include "AliMFTChipSegmentation.h"
#include "AliMFTFlex.h"
#include "AliMFTChip.h"
#include "AliMFTLadder.h"
#include "AliMFTConstants.h"
#include "AliMFTGeometry.h"

ClassImp(AliMFTLadder)

// Units are cm
const Double_t AliMFTLadder::kLadderDeltaY = AliMFTConstants::kSensorHeight + 2.*AliMFTConstants::kSensorTopOffset;
const Double_t AliMFTLadder::kLadderDeltaZ = AliMFTConstants::kFlexThickness + AliMFTConstants::kSensorThickness; // TODO: Adjust that value when adding glue layer

//=============================================================================================

AliMFTLadder::AliMFTLadder():
TNamed(), fMFTFlex(NULL){
  
  // default constructor
  
}
//=============================================================================================

AliMFTLadder::AliMFTLadder(AliMFTLadderSegmentation *segmentation):TNamed(segmentation->GetName(),segmentation->GetName()),
fSegmentation(segmentation), fMFTFlex(NULL)
{
  AliInfo(Form("Creating : %s", GetName()));
  Int_t nChips = fSegmentation->GetNumberOfChips();
  TGeoMedium *medVac  = gGeoManager->GetMedium("MFT_Air$");
  Double_t ladderDeltaX = AliMFTConstants::kLadderOffsetToEnd + AliMFTConstants::kSensorSideOffset + nChips * (AliMFTConstants::kSensorLength + AliMFTConstants::kSensorInterspace) ;
  fLadderVolume = gGeoManager->MakeBox(GetName(), medVac, ladderDeltaX/2., kLadderDeltaY/2., kLadderDeltaZ/2.);
  
  //fLadderVolume = new TGeoVolumeAssembly(GetName());

  AliWarning("TODO ### Flex To be added ###");
}


//=============================================================================================

AliMFTLadder::~AliMFTLadder() {
  delete fMFTFlex;
  
}

//=============================================================================================

TGeoVolume * AliMFTLadder::CreateVolume() {
//  AliWarning("TODO ### Positinning To be worked on ###");

  // Create the flex
  
//  fMFTFlex = new AliMFTFlex(fSegmentation);
//  TGeoVolumeAssembly * flexVol = fMFTFlex->MakeFlex(fSegmentation->GetNumberOfChips());
//  fLadderVolume->AddNode(flexVol,1);
  
  
  // Create the CMOS Sensors
  
  CreateSensors();
  

  return fLadderVolume;
  
}

//=============================================================================================

void AliMFTLadder::CreateSensors() {
  // Create Shapes
  
  // The sensor part
  TGeoBBox *sensor = new TGeoBBox(AliMFTConstants::kSensorLength/2.,AliMFTConstants::kSensorActiveHeight/2.,  AliMFTConstants::kSensorThickness/2.);
  
  // The readout part
  TGeoBBox *readout = new TGeoBBox(AliMFTConstants::kSensorLength/2.,(AliMFTConstants::kSensorHeight-AliMFTConstants::kSensorActiveHeight)/2.,  AliMFTConstants::kSensorThickness/2.);
  
  // Get Mediums
  TGeoMedium *medSensorSi  = gGeoManager->GetMedium("MFT_Si$");
  TGeoMedium *medReadoutSi = gGeoManager->GetMedium("MFT_Readout$");
  TGeoMedium *medAir  = gGeoManager->GetMedium("MFT_Air$");
  
  
  AliMFTGeometry * mftGeom = AliMFTGeometry::Instance();
  
  AliInfo(Form("Object Type = %d",mftGeom->GetObjectType(fSegmentation->GetUniqueID())));
  AliInfo(Form("Half-MFT    = %d",mftGeom->GetHalfMFTID(fSegmentation->GetUniqueID())));
  AliInfo(Form("Half-Disk   = %d",mftGeom->GetHalfDiskID(fSegmentation->GetUniqueID())));
  AliInfo(Form("Ladder      = %d",mftGeom->GetLadderID(fSegmentation->GetUniqueID())));
  
  
  
  TString namePrefix = Form("MFT_S_%d_%d_%d",
                                    mftGeom->GetHalfMFTID(fSegmentation->GetUniqueID()),
                                    mftGeom->GetHalfDiskID(fSegmentation->GetUniqueID()),
                                    mftGeom->GetLadderID(fSegmentation->GetUniqueID()) );
  
  TGeoVolume * chipVol = gGeoManager->MakeBox(namePrefix.Data(), medAir,AliMFTConstants::kSensorLength/2.,AliMFTConstants::kSensorHeight/2.,  AliMFTConstants::kSensorThickness/2. );

  // Create Volumes
  // Chip Volume
  chipVol->SetVisibility(kTRUE);
  
  // The sensor Volume
  TGeoVolume *sensorVol = new TGeoVolume("MFTSensor", sensor, medSensorSi);
  sensorVol->SetVisibility(kTRUE);
  sensorVol->SetLineColor(kGreen+1);
  sensorVol->SetLineWidth(1);
  sensorVol->SetFillColor(sensorVol->GetLineColor());
  sensorVol->SetFillStyle(4000); // 0% transparent
  if(!mftGeom->GetSensorVolumeID()){
    mftGeom->SetSensorVolumeID(sensorVol->GetNumber());
  } else if (mftGeom->GetSensorVolumeID() != sensorVol->GetNumber()){
    AliFatal(Form("Difference Sensor VOLUME ID in TGeo !!!!"));
  }
  
  // The Readout Volume
  TGeoVolume *readoutVol = new TGeoVolume("Readout", readout, medReadoutSi);
  readoutVol->SetVisibility(kTRUE);
  readoutVol->SetLineColor(kRed-6);
  readoutVol->SetLineWidth(1);
  readoutVol->SetFillColor(readoutVol->GetLineColor());
  readoutVol->SetFillStyle(4000); // 0% transparent

  // Building up the chip
  chipVol->AddNode(readoutVol, 1, new TGeoTranslation(0.,-AliMFTConstants::kSensorHeight/2.+readout->GetDY(),  0.));
  chipVol->AddNode(sensorVol, 1, new TGeoTranslation( 0., AliMFTConstants::kSensorHeight/2.-sensor->GetDY(),0.));

  for (int ichip =0; ichip<fSegmentation->GetNumberOfChips(); ichip++) {
    AliMFTChipSegmentation * chipSeg = fSegmentation->GetSensor(ichip);
    TGeoCombiTrans * chipPos = chipSeg->GetTransformation();
    // Position of the center on the chip in the chip coordinate system
    Double_t pos[3] ={AliMFTConstants::kSensorLength/2., AliMFTConstants::kSensorHeight/2., AliMFTConstants::kSensorThickness/2.};
    Double_t master[3];
    chipPos->LocalToMaster(pos, master);
    
    TGeoBBox* shape = (TGeoBBox*)fLadderVolume->GetShape();
    master[0] -= shape->GetDX();
    master[1] -= shape->GetDY();
    master[2] -= shape->GetDZ();
    AliInfo(Form("Adding Chip %s_%d ",namePrefix.Data(),ichip));
    fLadderVolume->AddNode(chipVol, ichip, new TGeoTranslation(master[0],master[1],master[2]));

  }

}