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
/// \class AliMFTLadder
///
/// Ladder Builder
///
// author Raphael Tieulent <raphael.tieulent@cern.ch>
//-----------------------------------------------------------------------------

#include "TGeoMatrix.h"
#include "TGeoManager.h"
#include "TGeoBBox.h"
#include "TGeoTube.h"

#include "AliLog.h"

#include "AliMFTLadderSegmentation.h"
#include "AliMFTChipSegmentation.h"
#include "AliMFTFlex.h"
#include "AliMFTChip.h"
#include "AliMFTLadder.h"
#include "AliMFTConstants.h"
#include "AliMFTGeometry.h"
#include "TGeoCompositeShape.h"
#include "TGeoBoolNode.h"

/// \cond CLASSIMP
ClassImp(AliMFTLadder);
/// \endcond

// Units are cm
const Double_t AliMFTLadder::kLadderDeltaY = AliMFTGeometry::kSensorHeight + 2.*AliMFTGeometry::kSensorTopOffset;
const Double_t AliMFTLadder::kLadderDeltaZ = AliMFTGeometry::kFlexThickness + AliMFTGeometry::kSensorThickness; // TODO: Adjust that value when adding glue layer

//=============================================================================================
/// \brief Default constructor

AliMFTLadder::AliMFTLadder():
TNamed(), fMFTFlex(NULL){
  
  
}
//=============================================================================================
/// \brief Constructor
AliMFTLadder::AliMFTLadder(AliMFTLadderSegmentation *segmentation):TNamed(segmentation->GetName(),segmentation->GetName()),fSegmentation(segmentation), fMFTFlex(NULL)
{
  AliDebug(1, Form("Creating : %s", GetName()));
  fLadderVolume = new TGeoVolumeAssembly(GetName());
  

}


//=============================================================================================

AliMFTLadder::~AliMFTLadder() {
  delete fMFTFlex;
  
}

//=============================================================================================
/// \brief Build the ladder
TGeoVolume * AliMFTLadder::CreateVolume() {

  Int_t nChips = fSegmentation->GetNSensors();

  // Create the flex
  fMFTFlex = new AliMFTFlex(fSegmentation);         
  Double_t kFlexLength = nChips*(AliMFTGeometry::kSensorLength+AliMFTGeometry::kSensorInterspace)+AliMFTGeometry::kLadderOffsetToEnd + AliMFTGeometry::kSensorSideOffset;
  Double_t kShiftY = 2*AliMFTGeometry::kSensorTopOffset+AliMFTGeometry::kSensorHeight-AliMFTGeometry::kFlexHeight/2; // strange
  TGeoVolumeAssembly * flexVol = fMFTFlex->MakeFlex(fSegmentation->GetNSensors(), kFlexLength);                               
  fLadderVolume->AddNode(flexVol, 1, new TGeoTranslation(kFlexLength/2+AliMFTGeometry::kSensorSideOffset/2, kShiftY, AliMFTGeometry::kFlexThickness/2-AliMFTGeometry::kRohacell));     

  // Create the CMOS Sensors
  CreateSensors();

  return fLadderVolume;
  
}



//=============================================================================================
/// \brief Build the sensors
void AliMFTLadder::CreateSensors() {
  // Create Shapes
  
  // The sensor part
  TGeoBBox *sensor = new TGeoBBox(AliMFTGeometry::kSensorLength/2., AliMFTGeometry::kSensorActiveHeight/2., AliMFTGeometry::kSensorThickness/2.);
  
  // The readout part
  TGeoBBox *readout = new TGeoBBox(AliMFTGeometry::kSensorLength/2.,(AliMFTGeometry::kSensorHeight-AliMFTGeometry::kSensorActiveHeight)/2.,  AliMFTGeometry::kSensorThickness/2.);
  
  // Get Mediums
  TGeoMedium *medSensorSi  = gGeoManager->GetMedium("MFT_Si$");
  TGeoMedium *medReadoutSi = gGeoManager->GetMedium("MFT_Readout$");
  TGeoMedium *medAir  = gGeoManager->GetMedium("MFT_Air$");
  //TGeoMedium *kMedGlue = gGeoManager->GetMedium("MFT_Epoxy$"); 
  TGeoMedium *kMedGlue = gGeoManager->GetMedium("MFT_SE4445$"); 

  
  AliMFTGeometry * mftGeom = AliMFTGeometry::Instance();
  
  TString namePrefix = Form("MFT_S_%d_%d_%d",
                                    mftGeom->GetHalfMFTID(fSegmentation->GetUniqueID()),
                                    mftGeom->GetHalfDiskID(fSegmentation->GetUniqueID()),
                                    mftGeom->GetLadderID(fSegmentation->GetUniqueID()) );
  
  TGeoVolume * chipVol = gGeoManager->MakeBox(namePrefix.Data(), medAir,AliMFTGeometry::kSensorLength/2.,AliMFTGeometry::kSensorHeight/2.,  AliMFTGeometry::kSensorThickness/2. );
  TGeoVolume * glue = gGeoManager->MakeBox(namePrefix.Data(), kMedGlue, (AliMFTGeometry::kSensorLength-AliMFTGeometry::kGlueEdge)/2., 
					   (AliMFTGeometry::kSensorHeight-AliMFTGeometry::kGlueEdge)/2., AliMFTGeometry::kGlueThickness/2.);
  glue->SetVisibility(kTRUE);
  glue->SetLineColor(kRed-10);
  glue->SetLineWidth(1);
  glue->SetFillColor(glue->GetLineColor());
  glue->SetFillStyle(4000); // 0% transparent

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
    AliFatal(Form("Different Sensor VOLUME ID in TGeo !!!!"));
  }
  
  // The Readout Volume
  TGeoVolume *readoutVol = new TGeoVolume("Readout", readout, medReadoutSi);
  readoutVol->SetVisibility(kTRUE);
  readoutVol->SetLineColor(kRed-6);
  readoutVol->SetLineWidth(1);
  readoutVol->SetFillColor(readoutVol->GetLineColor());
  readoutVol->SetFillStyle(4000); // 0% transparent

  // Building up the chip
  chipVol->AddNode(readoutVol, 1, new TGeoTranslation(0.,-AliMFTGeometry::kSensorHeight/2.+readout->GetDY(),  0.));
  chipVol->AddNode(sensorVol, 1, new TGeoTranslation( 0., AliMFTGeometry::kSensorHeight/2.-sensor->GetDY(),0.));

  for (int ichip =0; ichip<fSegmentation->GetNSensors(); ichip++) {
    AliMFTChipSegmentation * chipSeg = fSegmentation->GetSensor(ichip);
    TGeoCombiTrans * chipPos = chipSeg->GetTransformation();
    TGeoCombiTrans * chipPosGlue = chipSeg->GetTransformation();
    // Position of the center on the chip in the chip coordinate system
    Double_t pos[3] ={AliMFTGeometry::kSensorLength/2., AliMFTGeometry::kSensorHeight/2., AliMFTGeometry::kSensorThickness/2. - AliMFTGeometry::kGlueThickness - AliMFTGeometry::kRohacell};
    Double_t posglue[3] ={AliMFTGeometry::kSensorLength/2., AliMFTGeometry::kSensorHeight/2., AliMFTGeometry::kGlueThickness/2 - AliMFTGeometry::kRohacell};
    Double_t master[3];
    Double_t masterglue[3];
    chipPos->LocalToMaster(pos, master);
    chipPosGlue->LocalToMaster(posglue, masterglue);
    
    TGeoBBox* shape = (TGeoBBox*)fLadderVolume->GetShape();
    master[0] -= shape->GetDX();
    master[1] -= shape->GetDY();
    master[2] -= shape->GetDZ();

    masterglue[0] -= shape->GetDX();
    masterglue[1] -= shape->GetDY();
    masterglue[2] -= shape->GetDZ();

    AliDebug(1,Form("Adding Chip %s_%d ",namePrefix.Data(),ichip));
    fLadderVolume->AddNode(chipVol, ichip, new TGeoTranslation(master[0],master[1],master[2]));
    fLadderVolume->AddNode(glue, ichip, new TGeoTranslation(masterglue[0],masterglue[1],masterglue[2]));

  }

}
