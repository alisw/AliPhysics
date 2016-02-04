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
//      Class describing geometry of MFT CMOS MAP Chip
//
//      Contact author: raphael.tieulent@cern.ch
//
//=============================================================================================

#include "TGeoManager.h"
#include "TGeoBBox.h"

#include "AliMFTLadderSegmentation.h"
#include "AliMFTChipSegmentation.h"
#include "AliMFTChip.h"

ClassImp(AliMFTChip)

//=============================================================================================

AliMFTChip::AliMFTChip():
TNamed()
{
  
  // default constructor
  
}
//=============================================================================================

AliMFTChip::AliMFTChip(AliMFTChipSegmentation *segmentation, const char * ladderName):TNamed(ladderName,ladderName)
{
 
  

}


//=============================================================================================

AliMFTChip::~AliMFTChip() {
  
}


//=============================================================================================
void AliMFTChip::GetPosition(AliMFTLadderSegmentation * ladderSeg, Int_t iChip, Double_t *pos){
//  Double_t * fFlexDimensions = new Double_t[3];
//  ladderSeg->GetFlexLength(fFlexDimensions);
//  
//  AliInfo(Form("fFlexDimensions %f %f %f",fFlexDimensions[0],fFlexDimensions[1], fFlexDimensions[2]));
//  
//  pos[0] = AliMFTConstants::fChipSideOffset + AliMFTConstants::fChipWidth/2. + iChip*(AliMFTConstants::fChipWidth+AliMFTConstants::fChipInterspace);
//  pos[1] = -(AliMFTConstants::fChipTopOffset + fChipHeight/2.) ;
//  pos[2] =  fFlexDimensions[2] + AliMFTConstants::fChipThickness/2.;
//  AliWarning ("---- Z position of Chip to be worked out --- ");
//  if (!ladderSeg->IsLeftType()) pos[0]  *= -1.;
  
  
}


//=============================================================================================
TGeoVolume * AliMFTChip::CreateVolume(){
  
//  // Create Shapes
//  
//  // The sensor part
//  TGeoBBox *sensor = new TGeoBBox(AliMFTConstants::kSensorLength/2.,AliMFTConstants::kSensorHeight/2.,  AliMFTConstants::fChipThickness/2.);
//  
//  // The readout part
//  TGeoBBox *readout = new TGeoBBox(AliMFTConstants::fChipWidth/2.,(fChipHeight-fSensorHeight)/2.,  AliMFTConstants::fChipThickness/2.);
//  
//  // Get Mediums
//  TGeoMedium *medSensorSi  = gGeoManager->GetMedium("MFT_Si");
//  TGeoMedium *medReadoutSi = gGeoManager->GetMedium("MFT_Readout");
//
//  // Create Volumes
//  // Chip Volume
//  TGeoVolumeAssembly *chipVol = new TGeoVolumeAssembly("Chip");
//  chipVol->SetVisibility(kTRUE);
//
//  // The sensor Volume
//  TGeoVolume *sensorVol = new TGeoVolume("Sensor", sensor, medSensorSi);
//  sensorVol->SetVisibility(kTRUE);
//  sensorVol->SetLineColor(kGreen+1);
//
//  // The Readout Volume
//  TGeoVolume *readoutVol = new TGeoVolume("Readout", readout, medReadoutSi);
//  readoutVol->SetVisibility(kTRUE);
//  readoutVol->SetLineColor(kRed+2);
//
//  // Building up the chip
//  chipVol->AddNode(readoutVol, 1, new TGeoTranslation(0.,-fChipHeight/2.+readout->GetDY(),  0.));
//  chipVol->AddNode(sensorVol, 1, new TGeoTranslation( 0., fChipHeight/2.-sensor->GetDY(),0.));
//
//  
//  
//  return chipVol;
  
}
