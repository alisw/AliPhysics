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
//      Segmentation class for the planes of the ALICE Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "TFile.h"
#include "TNtuple.h"
#include "TClonesArray.h"
#include "TMath.h"
#include "AliMFTPlane.h"
#include "AliMFTSegmentation.h"

ClassImp(AliMFTSegmentation)

//====================================================================================================================================================

AliMFTSegmentation::AliMFTSegmentation(): 
  TObject(),
  fMFTPlanes(0)
{ 

  // TO BE CHECKED
  
  // default constructor

  fMFTPlanes = new TClonesArray("AliMFTPlane", fNMaxPlanes);

}

//====================================================================================================================================================

AliMFTSegmentation::AliMFTSegmentation(const Char_t *nameGeomFile): 
  TObject(),
  fMFTPlanes(0)
{ 

  fMFTPlanes = new TClonesArray("AliMFTPlane", fNMaxPlanes);

  Float_t zCenter, rMin, rMax, pixelSizeX, pixelSizeY, thicknessActive, thicknessSupport, thicknessReadout;
  Float_t equivalentSilicon, equivalentSiliconBeforeFront, equivalentSiliconBeforeBack;

  TFile *geomFile = new TFile(nameGeomFile);
  TNtuple *geomNtuple = (TNtuple*) geomFile->Get("AliMFTGeometry");

  geomNtuple -> SetBranchAddress("zCenter", &zCenter);
  geomNtuple -> SetBranchAddress("rMin",    &rMin);
  geomNtuple -> SetBranchAddress("rMax",    &rMax);
  geomNtuple -> SetBranchAddress("pixelSizeX", &pixelSizeX);
  geomNtuple -> SetBranchAddress("pixelSizeY", &pixelSizeY);
  geomNtuple -> SetBranchAddress("thicknessActive",  &thicknessActive);
  geomNtuple -> SetBranchAddress("thicknessSupport", &thicknessSupport);
  geomNtuple -> SetBranchAddress("thicknessReadout", &thicknessReadout);
  geomNtuple -> SetBranchAddress("equivalentSilicon",            &equivalentSilicon);
  geomNtuple -> SetBranchAddress("equivalentSiliconBeforeFront", &equivalentSiliconBeforeFront);
  geomNtuple -> SetBranchAddress("equivalentSiliconBeforeBack",  &equivalentSiliconBeforeBack);
  
  Int_t nPlanes = geomNtuple->GetEntries();

  for (Int_t iPlane=0; iPlane<nPlanes; iPlane++) {

    // Create new plane

    printf("Setting segmentation for MFT plane #%02d\n", iPlane);

    geomNtuple -> GetEntry(iPlane);
    zCenter = TMath::Abs(zCenter);

    AliMFTPlane *plane = new AliMFTPlane(Form("MFTPlane_%02d", iPlane), Form("MFTPlane_%02d", iPlane));
    plane -> Init(iPlane, zCenter, rMin, rMax, pixelSizeX, pixelSizeY, thicknessActive, thicknessSupport, thicknessReadout);
    plane -> SetEquivalentSilicon(equivalentSilicon);
    plane -> SetEquivalentSiliconBeforeFront(equivalentSiliconBeforeFront);
    plane -> SetEquivalentSiliconBeforeBack(equivalentSiliconBeforeBack);
    plane -> CreateStructure();
    
    new ((*fMFTPlanes)[fMFTPlanes->GetEntries()]) AliMFTPlane(*plane);

  }
  
  delete geomFile;

  printf("MFT segmentation set!\n");

}

//====================================================================================================================================================

THnSparseC* AliMFTSegmentation::GetDetElem(Int_t detElemID) {
      
  // Find det elem

  Int_t planeNb = detElemID/fNMaxDetElemPerPlane;
  Int_t detElemNb = detElemID - planeNb*fNMaxDetElemPerPlane;
  
  THnSparseC *detElem = GetPlane(planeNb)->GetActiveElement(detElemNb);

  return detElem;

}

//====================================================================================================================================================

Bool_t AliMFTSegmentation::Hit2PixelID(Double_t xHit, Double_t yHit, Int_t detElemID, Int_t &xPixel, Int_t &yPixel) {

  THnSparseC *detElem = GetDetElem(detElemID);

  if ( xHit<detElem->GetAxis(0)->GetXmin() ||
       xHit>detElem->GetAxis(0)->GetXmax() ||
       yHit<detElem->GetAxis(1)->GetXmin() ||
       yHit>detElem->GetAxis(1)->GetXmax() ) return kFALSE;

  xPixel = detElem->GetAxis(0)->FindBin(xHit) - 1;
  yPixel = detElem->GetAxis(1)->FindBin(yHit) - 1;

  return kTRUE;

}

//====================================================================================================================================================

