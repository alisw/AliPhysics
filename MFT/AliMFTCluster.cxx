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
//      Class for the description of the clusters of the ALICE Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "TObject.h"
#include "AliMUONRawCluster.h"
#include "AliMUONVCluster.h"
#include "AliMFTDigit.h"
#include "TMath.h"
#include "AliMFTCluster.h"

ClassImp(AliMFTCluster)

//====================================================================================================================================================

AliMFTCluster::AliMFTCluster():
  TObject(),
  fX(0), 
  fY(0), 
  fZ(0),
  fErrX(0), 
  fErrY(0), 
  fErrZ(0),
  fNElectrons(0),
  fNMCTracks(0),
  fPlane(-1),
  fDetElemID(-1),
  fSize(0),
  fTrackChi2(0),
  fLocalChi2(0),
  fDigitsInCluster(0),
  fIsClusterEditable(kTRUE),
  fIsClusterFront(kTRUE)
{

  // default constructor

  for (Int_t iTrack=0; iTrack<fNMaxMCTracks; iTrack++) fMCLabel[iTrack] = -1;

  fDigitsInCluster = new TClonesArray("AliMFTDigit", fNMaxDigitsPerCluster);
  fDigitsInCluster -> SetOwner(kTRUE);
}

//====================================================================================================================================================

AliMFTCluster::AliMFTCluster(const AliMFTCluster& cluster): 
  TObject(cluster),
  fX(cluster.fX), 
  fY(cluster.fY), 
  fZ(cluster.fZ),
  fErrX(cluster.fErrX), 
  fErrY(cluster.fErrY), 
  fErrZ(cluster.fErrZ),
  fNElectrons(cluster.fNElectrons),
  fNMCTracks(cluster.fNMCTracks),
  fPlane(cluster.fPlane),
  fDetElemID(cluster.fDetElemID),
  fSize(cluster.fSize),
  fTrackChi2(cluster.fTrackChi2),
  fLocalChi2(cluster.fLocalChi2),
  fDigitsInCluster(NULL),
  fIsClusterEditable(cluster.fIsClusterEditable),
  fIsClusterFront(cluster.fIsClusterFront)
{

  // copy constructor
  for (Int_t iTrack=0; iTrack<fNMaxMCTracks; iTrack++) fMCLabel[iTrack] = (cluster.fMCLabel)[iTrack];
  if (cluster.fDigitsInCluster) {
    fDigitsInCluster = new TClonesArray(*(cluster.fDigitsInCluster));
    fDigitsInCluster -> SetOwner(kTRUE);
  }
  else {
    fDigitsInCluster = new TClonesArray("AliMFTDigit", fNMaxDigitsPerCluster);
    fDigitsInCluster -> SetOwner(kTRUE);
  }    
 
}

//====================================================================================================================================================

AliMFTCluster& AliMFTCluster::operator=(const AliMFTCluster& cluster) {

  // Asignment operator

  // check assignement to self
  if (this == &cluster) return *this;

  // base class assignement
  TObject::operator=(cluster);
  
  // clear memory
  Clear("");
  
  fX                 = cluster.fX; 
  fY                 = cluster.fY; 
  fZ                 = cluster.fZ;
  fErrX              = cluster.fErrX; 
  fErrY              = cluster.fErrY; 
  fErrZ              = cluster.fErrZ;
  fNElectrons        = cluster.fNElectrons;
  fNMCTracks         = cluster.fNMCTracks;
  fPlane             = cluster.fPlane;
  fDetElemID         = cluster.fDetElemID;
  fSize              = cluster.fSize;
  fTrackChi2         = cluster.fTrackChi2;
  fLocalChi2         = cluster.fLocalChi2;
  fIsClusterEditable = cluster.fIsClusterEditable;
  fIsClusterFront    = cluster.fIsClusterFront;

  for (Int_t iTrack=0; iTrack<fNMaxMCTracks; iTrack++) fMCLabel[iTrack] = (cluster.fMCLabel)[iTrack];
  fDigitsInCluster      = new TClonesArray(*(cluster.fDigitsInCluster));
  fDigitsInCluster->SetOwner(kTRUE);

  return *this;

}

//====================================================================================================================================================

Double_t AliMFTCluster::GetDistanceFromPixel(AliMFTDigit *pixel) {

  // the distance is expressed in units of pixels!!!
  // useful to decide if the pixel is compatible with the current digits array

  Double_t distance = -1;

  if (!fSize) return distance;

  if (pixel->GetDetElemID()!=fDetElemID || pixel->GetPlane()!=fPlane) return 9999.;

  for (Int_t iDigit=0; iDigit<fDigitsInCluster->GetEntries(); iDigit++) {
    AliMFTDigit *tmpDig = (AliMFTDigit*) fDigitsInCluster->At(iDigit);
    Int_t distX = TMath::Abs(tmpDig->GetPixelX() - pixel->GetPixelX());
    Int_t distY = TMath::Abs(tmpDig->GetPixelY() - pixel->GetPixelY());
    if (distX<=1 &&  distY<=1) return 0;
    if (!iDigit) distance = TMath::Sqrt(distX*distX + distY*distY);
    else         distance = TMath::Min(distance, TMath::Sqrt(distX*distX + distY*distY));
  }

  return distance;

}

//====================================================================================================================================================

Bool_t AliMFTCluster::AddPixel(AliMFTDigit *pixel) {

  if (!fIsClusterEditable || fSize>=fNMaxDigitsPerCluster) return kFALSE;
  if (fSize && (pixel->GetPlane()!=fPlane || pixel->GetDetElemID()!=fDetElemID)) return kFALSE;

  new ((*fDigitsInCluster)[fDigitsInCluster->GetEntries()]) AliMFTDigit(*pixel);

  if (!fSize) {
    SetPlane(pixel->GetPlane());
    SetDetElemID(pixel->GetDetElemID());
  }

  fSize++;

  return kTRUE;

}

//====================================================================================================================================================

void AliMFTCluster::TerminateCluster() {

  Double_t xCenters[fNMaxDigitsPerCluster] = {0};
  Double_t yCenters[fNMaxDigitsPerCluster] = {0};
  Double_t nElectrons = 0.;

  for (Int_t iDigit=0; iDigit<fDigitsInCluster->GetEntries(); iDigit++) {
    AliMFTDigit *tmpDig = (AliMFTDigit*) fDigitsInCluster->At(iDigit);
    xCenters[iDigit] = tmpDig->GetPixelCenterX();
    yCenters[iDigit] = tmpDig->GetPixelCenterY();
    nElectrons      += tmpDig->GetNElectrons();
    for (Int_t iTrack=0; iTrack<tmpDig->GetNMCTracks(); iTrack++) AddMCLabel(tmpDig->GetMCLabel(iTrack));
  }

  SetX(TMath::Mean(fDigitsInCluster->GetEntries(), xCenters));
  SetY(TMath::Mean(fDigitsInCluster->GetEntries(), yCenters));
  SetZ(((AliMFTDigit*) fDigitsInCluster->At(0))->GetPixelCenterZ());

  Double_t minErrX = ((AliMFTDigit*) fDigitsInCluster->At(0))->GetPixelWidthX() / TMath::Sqrt(12.);
  Double_t minErrY = ((AliMFTDigit*) fDigitsInCluster->At(0))->GetPixelWidthY() / TMath::Sqrt(12.);
  Double_t minErrZ = ((AliMFTDigit*) fDigitsInCluster->At(0))->GetPixelWidthZ() / TMath::Sqrt(12.);
  SetErrX( TMath::Max(TMath::RMS(fDigitsInCluster->GetEntries(), xCenters), minErrX) );
  SetErrY( TMath::Max(TMath::RMS(fDigitsInCluster->GetEntries(), yCenters), minErrY) );
  SetErrZ( minErrZ );
    
  SetNElectrons(nElectrons);

  fIsClusterEditable = kFALSE;

}
  
//====================================================================================================================================================

void AliMFTCluster::AddMCLabel(Int_t label) { 

  if (fNMCTracks>=fNMaxMCTracks) return; 

  for (Int_t iTrack=0; iTrack<fNMCTracks; iTrack++) if (label==fMCLabel[iTrack]) return;

  fMCLabel[fNMCTracks++]=label; 

}

//====================================================================================================================================================

AliMUONRawCluster* AliMFTCluster::CreateMUONCluster() {

  AliMUONRawCluster *cluster = new AliMUONRawCluster();
  
  cluster->SetXYZ(GetX(), GetY(), GetZ());
  cluster->SetErrXY(GetErrX(),GetErrY());
  cluster->SetDetElemId(100);   // to get the cluster compatible with the AliMUONTrack::AddTrackParamAtCluster(...) method

  return cluster;

}

//====================================================================================================================================================
