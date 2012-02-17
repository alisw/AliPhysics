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
#include "AliMFTConstants.h"
#include "AliMFTClusterFinder.h"

const Double_t AliMFTClusterFinder::fCutForAvailableDigits = AliMFTConstants::fCutForAvailableDigits;
const Double_t AliMFTClusterFinder::fCutForAttachingDigits = AliMFTConstants::fCutForAttachingDigits;


ClassImp(AliMFTClusterFinder)

//====================================================================================================================================================

AliMFTClusterFinder::AliMFTClusterFinder() : 
  TObject(),
  fDigitsInCluster(0),
  fCurrentDigit(0),
  fCurrentCluster(0),
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

void AliMFTClusterFinder::Init(const Char_t *nameGeomFile) {

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
  Bool_t isDigAvailableForNewCluster = kTRUE;
  
  for (Int_t iPlane=0; iPlane<fNPlanes; iPlane++) {

    AliDebug(1, Form("Plane %02d", iPlane));

    TClonesArray *myDigitList = (TClonesArray*) pDigitList->At(iPlane);

    AliDebug(1, Form("myDigitList->GetEntries() = %d", myDigitList->GetEntries()));

    Int_t cycleOverDigits = 0;
    Double_t myCutForAvailableDigits = fCutForAvailableDigits;
    
    while (myDigitList->GetEntries()) {

      for (Int_t iDig=0; iDig<myDigitList->GetEntries(); iDig++) {

	AliDebug(1, Form("Check %d: Digit %5d of %5d", cycleOverDigits, iDig, myDigitList->GetEntries()));

	fCurrentDigit = (AliMFTDigit*) myDigitList->At(iDig);
	isDigAvailableForNewCluster = kTRUE;

	for (Int_t iCluster=0; iCluster<fClustersPerPlane[iPlane]->GetEntries(); iCluster++) {
	  fCurrentCluster = (AliMFTCluster*) fClustersPerPlane[iPlane]->At(iCluster);
	  AliDebug(2, Form("Distance between cluster and digit = %f",fCurrentCluster->GetDistanceFromPixel(fCurrentDigit))); 
	  if (fCurrentCluster->GetDistanceFromPixel(fCurrentDigit) < fCutForAttachingDigits) { 
	    fCurrentCluster->AddPixel(fCurrentDigit); 
	    myDigitList->Remove(fCurrentDigit);
	    myDigitList->Compress();
	    iDig--;
	    isDigAvailableForNewCluster = kFALSE;
	    break; 
	  }
	  if (fCurrentCluster->GetDistanceFromPixel(fCurrentDigit) < myCutForAvailableDigits) isDigAvailableForNewCluster=kFALSE;
	}

	if (isDigAvailableForNewCluster) {
	  AliMFTCluster *newCluster = new AliMFTCluster();
	  newCluster->AddPixel(fCurrentDigit);
	  myDigitList->Remove(fCurrentDigit);
	  myDigitList->Compress();
	  iDig--;
	  new ((*fClustersPerPlane[iPlane])[fClustersPerPlane[iPlane]->GetEntries()]) AliMFTCluster(*newCluster);
	}

      }   // end of cycle over the digits

      if (cycleOverDigits) myCutForAvailableDigits -= 0.5;
      cycleOverDigits++;

    }   // no more digits to check in current plane!

    for (Int_t iCluster=0; iCluster<fClustersPerPlane[iPlane]->GetEntries(); iCluster++) {
      fCurrentCluster = (AliMFTCluster*) fClustersPerPlane[iPlane]->At(iCluster);
      fCurrentCluster -> TerminateCluster();
    }

    AliDebug(1, Form("Found %d clusters in plane %02d", fClustersPerPlane[iPlane]->GetEntries(), iPlane));

  }  // end of cycle over the planes

}

//====================================================================================================================================================

void AliMFTClusterFinder::MakeClusterBranch(TTree *treeCluster) {

  // Creating the cluster branches, one for each plane (see AliMFTReconstructor::Reconstruct)

  AliDebug(1, "Making Cluster Branch");

  CreateClusters();

  if (treeCluster) {
    for(Int_t iPlane=0; iPlane<fNPlanes; iPlane++) {
      AliDebug(1, Form("Setting Branch Plane_%02d for treeCluster",iPlane));
      if (treeCluster->GetBranch(Form("Plane_%02d",iPlane))) continue;
      AliDebug(1, Form("Branch Plane_%02d does not exist, creating!",iPlane));
      treeCluster->Branch(Form("Plane_%02d",iPlane), &(fClustersPerPlane[iPlane])); 
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
