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
#include "AliMFTPlane.h"
#include "AliMFTSegmentation.h"
#include "TTree.h"
#include "TMath.h"
#include "AliMFTConstants.h"
#include "AliMFTClusterFinder.h"
#include "TRandom.h"
#include "AliMFTGeometry.h"
#include "AliCodeTimer.h"

const Double_t AliMFTClusterFinder::fCutForAvailableDigits = AliMFTConstants::fCutForAvailableDigits;
const Double_t AliMFTClusterFinder::fCutForAttachingDigits = AliMFTConstants::fCutForAttachingDigits;
const Double_t AliMFTClusterFinder::fMisalignmentMagnitude = AliMFTConstants::fMisalignmentMagnitude;

ClassImp(AliMFTClusterFinder)

//====================================================================================================================================================

AliMFTClusterFinder::AliMFTClusterFinder() : 
  TObject(),
  fDigitsInCluster(0),
  fCurrentDigit(0),
  fCurrentCluster(0),
  fSegmentation(0),
  fNPlanes(0),
  fApplyMisalignment(kFALSE),
  fStopWatch(0)
{

  // Default constructor
  AliInfo("start");
  
  for (Int_t iPlane=0; iPlane<fNMaxPlanes; iPlane++) fClustersPerPlane[iPlane] = NULL;
  fDigitsInCluster = new TClonesArray("AliMFTDigit", fNMaxDigitsPerCluster);
  fDigitsInCluster -> SetOwner(kTRUE);
  fStopWatch = new TStopwatch();
  fStopWatch -> Reset();
  AliDebug(1, "... done!");

}

//====================================================================================================================================================

AliMFTClusterFinder::~AliMFTClusterFinder() {
  
  AliDebug(1, "Deleting AliMFTClusterFinder...");
  
  for (Int_t iPlane=0; iPlane<fNMaxPlanes; iPlane++) {
    if (fClustersPerPlane[iPlane]) fClustersPerPlane[iPlane]->Delete(); delete fClustersPerPlane[iPlane]; fClustersPerPlane[iPlane] = 0x0;
  }
  if (fDigitsInCluster) fDigitsInCluster->Delete(); delete fDigitsInCluster; fDigitsInCluster = NULL;

  delete fStopWatch;
  
  AliDebug(1, "... done!");
  
}

//====================================================================================================================================================

void AliMFTClusterFinder::Clear(const Option_t* /*opt*/) {
  
  AliDebug(1, "Clearing AliMFTClusterFinder...");
  
  for (Int_t iPlane=0; iPlane<fNMaxPlanes; iPlane++) {
    if(fClustersPerPlane[iPlane]) fClustersPerPlane[iPlane]->Clear("C"); 
  }
  if(fDigitsInCluster) fDigitsInCluster->Clear("C");	
  
  AliDebug(1, "... done!");
  
}

//====================================================================================================================================================

void AliMFTClusterFinder::Init(const Char_t *nameGeomFile) {

  AliInfo("start");
  AliMFTGeometry *mftGeo = AliMFTGeometry::Instance();
  fSegmentation = mftGeo->GetSegmentation();
  fNPlanes = AliMFTConstants::kNDisks;
  AliInfo("Done");
}

//====================================================================================================================================================

void AliMFTClusterFinder::StartEvent() {

  // Cleaning up and preparation for the clustering procedure

  AliDebug(1, "Starting Event...");
  
  for (Int_t iPlane=0; iPlane<fNPlanes; iPlane++) {
    fClustersPerPlane[iPlane]->Delete();
  }

  AliDebug(1, "... done!");

}

//====================================================================================================================================================

void AliMFTClusterFinder::DigitsToClusters(const TObjArray *pDigitList) {
  AliCodeTimerAuto("",0);

  // where the clusterization is performed
  AliMFTGeometry *mftGeo = AliMFTGeometry::Instance();

  AliInfo("Starting Clusterization for MFT");
  AliDebug(1, Form("nPlanes = %d", fNPlanes));

  if (!fStopWatch) fStopWatch = new TStopwatch();
  fStopWatch -> Reset();

  StartEvent(); 
  Bool_t isDigAvailableForNewCluster = kTRUE;
  
  TClonesArray *myDigitList = 0;

  for (Int_t iPlane=0; iPlane<AliMFTConstants::kNDisks; iPlane++) {


    const Int_t nDetElem = mftGeo->GetDiskNSensors(iPlane);
    AliDebug(1, Form("Plane %02d : # Det Elem = %d ", iPlane,nDetElem));
    
    TClonesArray *clustersPerDetElem[nDetElem];
    for (Int_t iDetElem=0; iDetElem<nDetElem; iDetElem++) clustersPerDetElem[iDetElem] = new TClonesArray("AliMFTCluster");

    myDigitList = (TClonesArray*) pDigitList->At(iPlane);

    AliInfo( Form("myDigitList->GetEntries() = %d", myDigitList->GetEntries()));

    Int_t cycleOverDigits = 0;
    Double_t myCutForAvailableDigits = 0;
    
    Int_t currentDetElem = -1;
    Int_t currentDetElemLocal = -1;
    Bool_t areThereSkippedDigits = kFALSE;

    fStopWatch -> Start();

    while (myDigitList->GetEntries()) {

      for (Int_t iDig=0; iDig<myDigitList->GetEntries(); iDig++) {

	fCurrentDigit = (AliMFTDigit*) myDigitList->At(iDig);

	if (!iDig) {
	  if (fCurrentDigit->GetDetElemID() != currentDetElem) {      
	    // first iteration over the digits of a specific detection element
	    currentDetElem = fCurrentDigit->GetDetElemID();
      currentDetElemLocal = mftGeo->GetDetElemLocalID(currentDetElem);
	    cycleOverDigits = 0;
	    myCutForAvailableDigits = fCutForAvailableDigits;
	  }
	  else if (fCurrentDigit->GetDetElemID()==currentDetElem && areThereSkippedDigits) {
	    // second (or further) iteration over the digits of a specific detection element
	    cycleOverDigits++;
	    myCutForAvailableDigits -= 0.5;
	  }
	  areThereSkippedDigits = kFALSE;
	}
	else {
	  areThereSkippedDigits = kTRUE;
	  if (fCurrentDigit->GetDetElemID() != currentDetElem) break;
	}

	isDigAvailableForNewCluster = kTRUE;

	for (Int_t iCluster=0; iCluster<clustersPerDetElem[currentDetElemLocal]->GetEntries(); iCluster++) {
	  fCurrentCluster = (AliMFTCluster*) clustersPerDetElem[currentDetElemLocal]->At(iCluster);
	  if (fCurrentCluster->GetDistanceFromPixel(fCurrentDigit) < fCutForAttachingDigits) { 
	    fCurrentCluster->AddPixel(fCurrentDigit); 
	    myDigitList->Remove(fCurrentDigit);
	    myDigitList->Compress();
	    iDig--;
	    isDigAvailableForNewCluster = kFALSE;
	    break; 
	  }
	  if (fCurrentCluster->GetDistanceFromPixel(fCurrentDigit) < myCutForAvailableDigits) {
	    areThereSkippedDigits = kTRUE;
	    isDigAvailableForNewCluster=kFALSE;
	  }
	}

	if (isDigAvailableForNewCluster) {
	  AliMFTCluster *newCluster = new AliMFTCluster();
	  newCluster->AddPixel(fCurrentDigit);
	  myDigitList->Remove(fCurrentDigit);
	  myDigitList->Compress();
	  iDig--;
	  new ((*clustersPerDetElem[currentDetElemLocal])[clustersPerDetElem[currentDetElemLocal]->GetEntries()]) AliMFTCluster(*newCluster);
	  delete newCluster;
	}
	
      }   // end of cycle over the digits

    }   // no more digits to check in current plane!

    fStopWatch -> Print("m");

    AliInfo(Form("Plane %d: clusters found in %f seconds\n",iPlane,fStopWatch->CpuTime()));
    fStopWatch->Start();

    // Now we merge the cluster lists coming from each detection element, to build the cluster list of the plane

    Double_t misalignmentPhi = 2.*TMath::Pi()*gRandom->Rndm();
    Double_t misalignmentX   = fMisalignmentMagnitude*TMath::Cos(misalignmentPhi);
    Double_t misalignmentY   = fMisalignmentMagnitude*TMath::Sin(misalignmentPhi);

    AliMFTCluster *newCluster = NULL;
    for (Int_t iDetElem=0; iDetElem<nDetElem; iDetElem++) {
      for (Int_t iCluster=0; iCluster<clustersPerDetElem[iDetElem]->GetEntries(); iCluster++) {
	newCluster = (AliMFTCluster*) (clustersPerDetElem[iDetElem]->At(iCluster));
	newCluster -> TerminateCluster();
	if (fApplyMisalignment) {
	  newCluster -> SetClusterEditable(kTRUE);
	  newCluster -> SetX(newCluster->GetX()+misalignmentX);
	  newCluster -> SetY(newCluster->GetY()+misalignmentY);
	  newCluster -> SetClusterEditable(kFALSE);
	}
	new ((*fClustersPerPlane[iPlane])[fClustersPerPlane[iPlane]->GetEntries()]) AliMFTCluster(*newCluster);
      }
    }

    printf("%d Clusters in plane %02d merged in %f seconds\n", fClustersPerPlane[iPlane]->GetEntries(), iPlane, fStopWatch->CpuTime());

    for (Int_t iDetElem=0; iDetElem<nDetElem; iDetElem++) {
      clustersPerDetElem[iDetElem] -> Delete();
      delete clustersPerDetElem[iDetElem];
    }

    myDigitList -> Delete();

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
    fClustersPerPlane[iPlane] -> SetOwner(kTRUE);

  }

}

//====================================================================================================================================================
