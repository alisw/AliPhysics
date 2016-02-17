#include "AliMFTCluster.h"
#include "TFile.h"
#include "AliMFTCluster.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "AliMFTConstants.h"
#include "TBranch.h"
#include "TRandom.h"

TFile *fileIn=0, *fileOut=0;

const Int_t labelMCOffset = 1000000;
Double_t misalignmentX[AliMFTConstants::fNMaxPlanes] = {0};
Double_t misalignmentY[AliMFTConstants::fNMaxPlanes] = {0};

//====================================================================================================================================================

void AddMisalignmentToClusters(Char_t *nameDir=".",
			       Int_t seed = 12345,
			       Double_t misalignment = 0.0015) {

  gRandom -> SetSeed(seed);

  TClonesArray *fRecPointsPerPlaneIn[AliMFTConstants::fNMaxPlanes] = {0};
  TClonesArray *fRecPointsPerPlaneOut[AliMFTConstants::fNMaxPlanes] = {0};

  fileIn = new TFile(Form("%s/MFT.RecPoints.root",nameDir));
  fileOut = new TFile(Form("%s/MFT.RecPoints.Misaligned.root",nameDir), "recreate");
  
  Int_t iEv=0;
  while (fileIn->cd(Form("Event%d",iEv))) {

    for (Int_t iPlane=0; iPlane<AliMFTConstants::fNMaxPlanes; iPlane++) {
      fRecPointsPerPlaneIn[iPlane] = new TClonesArray("AliMFTCluster");
      fRecPointsPerPlaneOut[iPlane] = new TClonesArray("AliMFTCluster");
    }

    printf("Event %d\n",iEv);
    TTree *treeIn = (TTree*) gDirectory->Get("TreeR");
    treeIn -> SetName("TreeR_In");
    fileOut-> cd();
    TTree *treeOut = new TTree("TreeR", "Reconstructed Points Container");

    Int_t iPlane=0;
    while (treeIn->GetBranch(Form("Plane_%02d",iPlane))) {
      //      printf("Plane %02d\n",iPlane);
      Double_t misalignmentPhi = 2.*TMath::Pi()*gRandom->Rndm();
      misalignmentX[iPlane] = misalignment*TMath::Cos(misalignmentPhi);
      misalignmentY[iPlane] = misalignment*TMath::Sin(misalignmentPhi);
      treeIn ->SetBranchAddress(Form("Plane_%02d",iPlane), &(fRecPointsPerPlaneIn[iPlane]));
      treeOut->Branch(Form("Plane_%02d",iPlane), &(fRecPointsPerPlaneOut[iPlane])); 
      //      treeOut->SetBranchAddress(Form("Plane_%02d",iPlane), &(fRecPointsPerPlaneOut[iPlane]));
      iPlane++;
    }

    iPlane=0;
    treeIn -> GetEntry(0);
    while (treeIn->GetBranch(Form("Plane_%02d",iPlane))) {
      Int_t nClusters = fRecPointsPerPlaneIn[iPlane]->GetEntries();
      for (Int_t iCluster=0; iCluster<nClusters; iCluster++) {
	//	printf("Cluster %4d\n",iCluster);
	AliMFTCluster *newCluster = (AliMFTCluster*) fRecPointsPerPlaneIn[iPlane]->At(iCluster);
	newCluster->SetClusterEditable(kTRUE);
	newCluster->SetX(newCluster->GetX()+misalignmentX[iPlane]);
	newCluster->SetY(newCluster->GetY()+misalignmentY[iPlane]);
	new ((*fRecPointsPerPlaneOut[iPlane])[fRecPointsPerPlaneOut[iPlane]->GetEntries()]) AliMFTCluster(*newCluster);
      }
      iPlane++;
    }

    treeOut -> Fill();
    fileOut -> mkdir(Form("Event%d",iEv));
    fileOut -> cd(Form("Event%d",iEv));
    treeOut -> Write();
    for (Int_t jPlane=0; jPlane<AliMFTConstants::fNMaxPlanes; jPlane++) {
      delete fRecPointsPerPlaneIn[jPlane];
      delete fRecPointsPerPlaneOut[jPlane];
    }

    iEv++;

  }

  fileOut -> Close();

}

//====================================================================================================================================================
