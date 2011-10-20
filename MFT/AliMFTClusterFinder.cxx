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
//      Class for finding and building the clusters of the ALICE Muon Forward Tracker
//
//      Contact author: antonio.uras@cern.ch
//
//====================================================================================================================================================

#include "AliLog.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "AliMFTDigit.h"
#include "AliMFTCluster.h"
#include "AliMFTSegmentation.h"
#include "TTree.h"
#include "TMath.h"
#include "AliMFTClusterFinder.h"

ClassImp(AliMFTClusterFinder)

//====================================================================================================================================================

AliMFTClusterFinder::AliMFTClusterFinder() : 
  TObject(),
  fDigitsInCluster(0),
  fCurrentDig(0),
  fSegmentation(0),
  fNPlanes(0)
{

  // Default constructor
  
  for (Int_t iPlane=0; iPlane<fNMaxPlanes; iPlane++) fClustersPerPlane[iPlane] = NULL;
  fDigitsInCluster = new TClonesArray("AliMFTDigit", fNMaxDigitsPerCluster);

}

//====================================================================================================================================================

AliMFTClusterFinder::~AliMFTClusterFinder() {

  AliDebug(1, "Deleting AliMFTClusterFinder...");

  for (Int_t iPlane=0; iPlane<fNMaxPlanes; iPlane++) delete fClustersPerPlane[iPlane];

  AliDebug(1, "... done!");

}

//====================================================================================================================================================

void AliMFTClusterFinder::Init(Char_t *nameGeomFile) {

  fSegmentation = new AliMFTSegmentation(nameGeomFile);
  fNPlanes = fSegmentation -> GetNPlanes();
  
}

//====================================================================================================================================================

void AliMFTClusterFinder::StartEvent() {

  // Cleaning up and preparation for the clustering procedure

  AliDebug(1, "Starting Event...");
  
  for (Int_t iPlane=0; iPlane<fNPlanes; iPlane++) {
    fClustersPerPlane[iPlane]->Clear();
  }

  AliDebug(1, "... done!");

}

//====================================================================================================================================================

void AliMFTClusterFinder::DigitsToClusters(const TObjArray *pDigitList) {

  // where the clusterization is performed

  AliInfo("Starting Clusterization for MFT");
  AliDebug(1, Form("nPlanes = %d", fNPlanes));

  StartEvent(); 
  
  for (Int_t iPlane=0; iPlane<fNPlanes; iPlane++) {

    AliDebug(1, Form("Plane %02d", iPlane));

    TClonesArray *myDigitList = (TClonesArray*) pDigitList->At(iPlane);
    //    myDigitList->Sort();

    AliDebug(1, Form("myDigitList->GetEntries() = %d", myDigitList->GetEntries()));

    while (myDigitList->GetEntries()) {

      fDigitsInCluster->Clear();
      
      Bool_t clusterUpdated=kTRUE;

      //---------------------------------------------------------------------------------------------------------

      while (clusterUpdated) {    // repeat the loop on the digits until no new compatible digit is found

	clusterUpdated = kFALSE;

	for (Int_t iDig=0; iDig<myDigitList->GetEntries(); iDig++) {
	  fCurrentDig = (AliMFTDigit*) myDigitList->At(iDig);
	  if (fDigitsInCluster->GetEntries()<fNMaxDigitsPerCluster) {
	    if (fDigitsInCluster->GetEntries()==0) {
	      new ((*fDigitsInCluster)[fDigitsInCluster->GetEntries()]) AliMFTDigit(*fCurrentDig);
	      myDigitList->Remove(fCurrentDig);
	      clusterUpdated = kTRUE;
	    }
	    else if (IsCurrentDigitCompatible()) {
	      new ((*fDigitsInCluster)[fDigitsInCluster->GetEntries()]) AliMFTDigit(*fCurrentDig);
	      myDigitList->Remove(fCurrentDig);
	      clusterUpdated = kTRUE;
	    }
	  }
	}
	
	if (clusterUpdated) myDigitList->Compress();

      }

      //---------------------------------------------------------------------------------------------------------

      AliDebug(1, "Building new cluster");

      BuildNewCluster(iPlane);  // here the new cluster is built

    }

  }

}

//====================================================================================================================================================

Bool_t AliMFTClusterFinder::IsCurrentDigitCompatible() {

  // where it is decided if the current digit (fCurrentDig) is compatible with the current digits array (fDigitsInCluster)

  for (Int_t iDigit=0; iDigit<fDigitsInCluster->GetEntries(); iDigit++) {
    AliMFTDigit *tmpDig = (AliMFTDigit*) fDigitsInCluster->At(iDigit);
    Int_t distX = TMath::Abs(tmpDig->GetPixelX() - fCurrentDig->GetPixelX());
    Int_t distY = TMath::Abs(tmpDig->GetPixelY() - fCurrentDig->GetPixelY());
    if (distX<=1 &&  distY<=1) return kTRUE;
  }

  return kFALSE;
    
}

//====================================================================================================================================================

void AliMFTClusterFinder::BuildNewCluster(Int_t plane) {

  // where a new cluster is built, starting from the array of compatible digits (fDigitsInCluster)

  AliDebug(1, Form("Starting cluster building from %d digits", fDigitsInCluster->GetEntries()));

  AliMFTCluster *newCluster = new AliMFTCluster();

  Double_t xCenters[fNMaxDigitsPerCluster] = {0};
  Double_t yCenters[fNMaxDigitsPerCluster] = {0};
  Double_t nElectrons = 0.;

  for (Int_t iDigit=0; iDigit<fDigitsInCluster->GetEntries(); iDigit++) {
    AliMFTDigit *tmpDig = (AliMFTDigit*) fDigitsInCluster->At(iDigit);
    xCenters[iDigit] = tmpDig->GetPixelCenterX();
    yCenters[iDigit] = tmpDig->GetPixelCenterY();
    nElectrons      += tmpDig->GetNElectrons();
    for (Int_t iTrack=0; iTrack<tmpDig->GetNMCTracks(); iTrack++) newCluster->AddMCLabel(tmpDig->GetMCLabel(iTrack));
  }

  newCluster -> SetX(TMath::Mean(fDigitsInCluster->GetEntries(), xCenters));
  newCluster -> SetY(TMath::Mean(fDigitsInCluster->GetEntries(), yCenters));
  newCluster -> SetZ(((AliMFTDigit*) fDigitsInCluster->At(0))->GetPixelCenterZ());

  Double_t minErrX = ((AliMFTDigit*) fDigitsInCluster->At(0))->GetPixelWidthX() / TMath::Sqrt(12.);
  Double_t minErrY = ((AliMFTDigit*) fDigitsInCluster->At(0))->GetPixelWidthY() / TMath::Sqrt(12.);
  Double_t minErrZ = ((AliMFTDigit*) fDigitsInCluster->At(0))->GetPixelWidthZ() / TMath::Sqrt(12.);
  newCluster -> SetErrX( TMath::Max(TMath::RMS(fDigitsInCluster->GetEntries(), xCenters), minErrX) );
  newCluster -> SetErrY( TMath::Max(TMath::RMS(fDigitsInCluster->GetEntries(), yCenters), minErrY) );
  newCluster -> SetErrZ( minErrZ );
    
  newCluster -> SetNElectrons(nElectrons);
  newCluster -> SetPlane(plane);

  newCluster -> SetSize(fDigitsInCluster->GetEntries());

  AliDebug(1, Form("Adding cluster #%02d to tree (%f, %f, %f)", 
		   fClustersPerPlane[plane]->GetEntries(), newCluster->GetX(), newCluster->GetY(), newCluster->GetZ()));

  new ((*fClustersPerPlane[plane])[fClustersPerPlane[plane]->GetEntries()]) AliMFTCluster(*newCluster);

}

//====================================================================================================================================================

void AliMFTClusterFinder::MakeClusterBranch(TTree *treeCluster) {

  // Creating the cluster branches, one for each plane (see AliMFTReconstructor::Reconstruct)

  AliDebug(1, "Making Cluster Branch");

  CreateClusters();

  if (treeCluster) {
    for(Int_t iPlane=0; iPlane<fNPlanes; iPlane++) {
      AliDebug(1, Form("Setting Branch Plane_%02d for treeCluster",iPlane));
      TBranch *branch = treeCluster->GetBranch(Form("Plane_%02d",iPlane));
      if (branch) continue;
      AliDebug(1, Form("Branch Plane_%02d does not exist, creating!",iPlane));
      branch = treeCluster->Branch(Form("Plane_%02d",iPlane), &(fClustersPerPlane[iPlane])); 
    }
  }

}

//====================================================================================================================================================

void AliMFTClusterFinder::SetClusterTreeAddress(TTree *treeCluster) {

  // Addressing the cluster branches, one for each plane (see AliMFTReconstructor::Reconstruct)

  if (treeCluster && treeCluster->GetBranch("Plane_00")) {
    CreateClusters();
    for(Int_t iPlane=0; iPlane<fNPlanes; iPlane++) {
      if (treeCluster->GetBranch(Form("Plane_%02d",iPlane))) {
	treeCluster->SetBranchAddress(Form("Plane_%02d",iPlane), &(fClustersPerPlane[iPlane]));
      }
      else AliError(Form("No branch available with name Plane_%02d", iPlane));
    }
  }

}

//====================================================================================================================================================

void AliMFTClusterFinder::CreateClusters() {

  // create cluster list

  AliDebug(1, Form("Creating clusters list: nPlanes = %d",fNPlanes));

  if (fClustersPerPlane[0]) return;
  
  for(Int_t iPlane=0; iPlane<fNPlanes; iPlane++) {
    AliDebug(1, Form("plane %02d", iPlane));
    fClustersPerPlane[iPlane] = new TClonesArray("AliMFTCluster");
  }

}

//====================================================================================================================================================
