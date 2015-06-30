#include "AliMFTCluster.h"
#include "TFile.h"
#include "AliMFTCluster.h"
#include "TObjArray.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "AliMFTConstants.h"
#include "TBranch.h"
#include "TRandom.h"
#include "TMath.h"

TFile *fileInPileUp=0, *fileInSignal, *fileInUnderlying, *fileOutPileUp=0, *fileOutTotal;

enum {kMergePileUpClusters_pp, kMergePileUpClusters_PbPb, kMergeAllClusters_pp, kMergeAllClusters_PbPb};
enum {k_pp, k_PbPb};

const Char_t *fPileUpDir[2] = {"$UnderlyingDirPP/spreadVtx",      // local directory with pile-up event p-p
			       "$UnderlyingDirPbPb/spreadVtx"};   // local directory with pile-up event Pb-Pb

const Int_t fNRunsPileUp_pp   = 100;
const Int_t fNRunsPileUp_PbPb =   3;

const Int_t fNRunsPileUpAvailable[2] = {500, 500};    // pp, PbPb

TClonesArray *fRecPointsPerPlaneIn[AliMFTConstants::fNMaxPlanes] = {0};
TClonesArray *fRecPointsPerPlaneOut[AliMFTConstants::fNMaxPlanes] = {0};

void MergePileUpClusters(Int_t *runs, Int_t nRunsPileUp, const Char_t *pileUpDir, Int_t nEvents);
void MergeAllClusters();
			       
//====================================================================================================================================================

void MergeClustersMFT(Int_t seed = 12345,
		      Int_t option = 9999,
		      Int_t nEvents = 9999,
		      Int_t maxPileUp = 9999) {
  
  gRandom->SetSeed(seed);

  if (option==kMergePileUpClusters_pp) {

    Int_t runs_pp[fNRunsPileUp_pp] = {0};
    Int_t nFiles_pp = 0;
    while (nFiles_pp<TMath::Min(fNRunsPileUp_pp,maxPileUp)) {
      Int_t run = gRandom->Integer(fNRunsPileUpAvailable[k_pp]);
      fileInPileUp = new TFile(Form("%s/run_%d/MFT.RecPoints.MCShifted.root",fPileUpDir[k_pp],run));
      if (!fileInPileUp->IsOpen() || !fileInPileUp->cd("Event0")) {
	delete fileInPileUp;
	continue;
      }
      Bool_t runAlreadyLoaded=kFALSE;
      for (Int_t iAddedRun=0; iAddedRun<nFiles_pp; iAddedRun++) {
	if (run==runs_pp[iAddedRun]) {
	  runAlreadyLoaded=kTRUE;
	  break;
	}
      }
      if (!runAlreadyLoaded) runs_pp[nFiles_pp++] = run;
      delete fileInPileUp;
    }
    MergePileUpClusters(runs_pp, TMath::Min(fNRunsPileUp_pp,maxPileUp), fPileUpDir[k_pp], nEvents);
  }
    
  else if (option==kMergePileUpClusters_PbPb) {
    Int_t runs_PbPb[fNRunsPileUp_PbPb] = {0};
    Int_t nFiles_PbPb = 0;
    while (nFiles_PbPb<TMath::Min(fNRunsPileUp_PbPb,maxPileUp)) {
      Int_t run = gRandom->Integer(fNRunsPileUpAvailable[k_PbPb]);
      fileInPileUp = new TFile(Form("%s/run_%d/MFT.RecPoints.MCShifted.root",fPileUpDir[k_PbPb],run));
      if (!fileInPileUp->IsOpen() || !fileInPileUp->cd("Event0")) {
	delete fileInPileUp;
	continue;
      }
      Bool_t runAlreadyLoaded=kFALSE;
      for (Int_t iAddedRun=0; iAddedRun<nFiles_PbPb; iAddedRun++) {
	if (run==runs_PbPb[iAddedRun]) {
	  runAlreadyLoaded=kTRUE;
	  break;
	}
      }
      if (!runAlreadyLoaded) runs_PbPb[nFiles_PbPb++] = run;
      delete fileInPileUp;
    }
    MergePileUpClusters(runs_PbPb, TMath::Min(fNRunsPileUp_PbPb,maxPileUp), fPileUpDir[k_PbPb], nEvents);
  }

  if (option>=kMergeAllClusters_pp) MergeAllClusters();

}

//====================================================================================================================================================

void MergePileUpClusters(Int_t *runs, Int_t nRunsPileUp, const Char_t *pileUpDir, Int_t nEvents) {

  fileOutPileUp = new TFile("./MFT.RecPoints.PileUp.root", "recreate");
    
  Bool_t eventExist = kTRUE;
  Int_t iEv=0;

  while (eventExist && iEv<nEvents) {

    for (Int_t iPlane=0; iPlane<AliMFTConstants::fNMaxPlanes; iPlane++) {
      fRecPointsPerPlaneOut[iPlane] = new TClonesArray("AliMFTCluster");
    }

    if (iEv) fileOutPileUp = new TFile("./MFT.RecPoints.PileUp.root", "update");

    printf("Merging Event %d\n",iEv);
    
    fileOutPileUp-> cd();
    TTree *treeOutPileUp = new TTree("TreeR_OutPileUp", "Reconstructed Points Container");

    Int_t nFilesAdded = 0;
    while (nFilesAdded<nRunsPileUp) {
      fileInPileUp = new TFile(Form("%s/run_%d/MFT.RecPoints.MCShifted.root",pileUpDir,runs[nFilesAdded++]));
      printf("Merging file %s\n", fileInPileUp->GetName());
      if (!fileInPileUp->cd(Form("Event%d",iEv))) {
	eventExist = kFALSE;
	break;
      }
      for (Int_t iPlane=0; iPlane<AliMFTConstants::fNMaxPlanes; iPlane++) {
	fRecPointsPerPlaneIn[iPlane] = new TClonesArray("AliMFTCluster");
      }
      TTree *treeInPileUp = (TTree*) gDirectory->Get("TreeR");
      treeInPileUp -> SetName("TreeR_InPileUp");
      Int_t iPlane=0;
      while (treeInPileUp->GetBranch(Form("Plane_%02d",iPlane))) {
	//	printf("Plane %02d\n",iPlane);
	treeInPileUp ->SetBranchAddress(Form("Plane_%02d",iPlane), &(fRecPointsPerPlaneIn[iPlane]));
	TBranch *branch = treeOutPileUp->GetBranch(Form("Plane_%02d",iPlane));
	if (!branch) treeOutPileUp->Branch(Form("Plane_%02d",iPlane), &(fRecPointsPerPlaneOut[iPlane])); 
	iPlane++;
      }
      iPlane=0;
      treeInPileUp -> GetEntry(0);
      while (treeInPileUp->GetBranch(Form("Plane_%02d",iPlane))) {
	Int_t nClusters = fRecPointsPerPlaneIn[iPlane]->GetEntries();
	for (Int_t iCluster=0; iCluster<nClusters; iCluster++) {
	  //	  printf("Cluster %4d\n",iCluster);
	  AliMFTCluster *newCluster = (AliMFTCluster*) fRecPointsPerPlaneIn[iPlane]->At(iCluster);
	  new ((*fRecPointsPerPlaneOut[iPlane])[fRecPointsPerPlaneOut[iPlane]->GetEntries()]) AliMFTCluster(*newCluster);
	}
	iPlane++;
      }
      for (Int_t jPlane=0; jPlane<AliMFTConstants::fNMaxPlanes; jPlane++) {
	delete fRecPointsPerPlaneIn[jPlane];
      }
      delete treeInPileUp;
      fileInPileUp -> Close();
      delete fileInPileUp;
    }
    if (eventExist) {
      treeOutPileUp -> Fill();
      fileOutPileUp -> mkdir(Form("Event%d",iEv));
      fileOutPileUp -> cd(Form("Event%d",iEv));
      treeOutPileUp -> SetName("TreeR");
      treeOutPileUp -> Write();
      delete treeOutPileUp;
      fileOutPileUp -> Close();
      delete fileOutPileUp;
      for (Int_t iPlane=0; iPlane<AliMFTConstants::fNMaxPlanes; iPlane++) {
	delete fRecPointsPerPlaneOut[iPlane];
      }
      iEv++;
    }
  }  

}

//====================================================================================================================================================

void MergeAllClusters() {

  fileOutTotal = new TFile("./MFT.RecPoints.root", "recreate");
    
  Bool_t eventExist = kTRUE;
  Int_t iEv=0;

  while (eventExist) {

    for (Int_t iPlane=0; iPlane<AliMFTConstants::fNMaxPlanes; iPlane++) {
      fRecPointsPerPlaneOut[iPlane] = new TClonesArray("AliMFTCluster");
    }

    if (iEv) fileOutTotal = new TFile("./MFT.RecPoints.root", "update");

    printf("Merging Event %d\n",iEv);
    
    fileOutTotal-> cd();
    TTree *treeOutTotal = new TTree("TreeR_Out", "Reconstructed Points Container");

    // Adding Pile-Up Clusters

    fileInPileUp = new TFile("./MFT.RecPoints.PileUp.root");
    printf("Merging file %s\n", fileInPileUp->GetName());
    if (!fileInPileUp->cd(Form("Event%d",iEv))) {
      eventExist = kFALSE;
      break;
    }
    for (Int_t jPlane=0; jPlane<AliMFTConstants::fNMaxPlanes; jPlane++) {
      fRecPointsPerPlaneIn[jPlane] = new TClonesArray("AliMFTCluster");
    }
    TTree *treeInPileUp = (TTree*) gDirectory->Get("TreeR");
    treeInPileUp -> SetName("TreeR_InPileUp");
    Int_t iPlane=0;
    while (treeInPileUp->GetBranch(Form("Plane_%02d",iPlane))) {
      //	printf("Plane %02d\n",iPlane);
      treeInPileUp ->SetBranchAddress(Form("Plane_%02d",iPlane), &(fRecPointsPerPlaneIn[iPlane]));
      TBranch *branch = treeOutTotal->GetBranch(Form("Plane_%02d",iPlane));
      if (!branch) treeOutTotal->Branch(Form("Plane_%02d",iPlane), &(fRecPointsPerPlaneOut[iPlane])); 
      iPlane++;
    }
    iPlane=0;
    treeInPileUp -> GetEntry(0);
    while (treeInPileUp->GetBranch(Form("Plane_%02d",iPlane))) {
      Int_t nClusters = fRecPointsPerPlaneIn[iPlane]->GetEntries();
      for (Int_t iCluster=0; iCluster<nClusters; iCluster++) {
	//	  printf("Cluster %4d\n",iCluster);
	AliMFTCluster *newCluster = (AliMFTCluster*) fRecPointsPerPlaneIn[iPlane]->At(iCluster);
	new ((*fRecPointsPerPlaneOut[iPlane])[fRecPointsPerPlaneOut[iPlane]->GetEntries()]) AliMFTCluster(*newCluster);
      }
      iPlane++;
    }
    for (Int_t jPlane=0; jPlane<AliMFTConstants::fNMaxPlanes; jPlane++) {
      delete fRecPointsPerPlaneIn[jPlane];
    }
    delete treeInPileUp;
    fileInPileUp -> Close();
    delete fileInPileUp;

    // Adding Underlying Event Clusters

    fileInUnderlying = new TFile("./MFT.RecPoints.Underlying.root");
    printf("Merging file %s\n", fileInUnderlying->GetName());
    if (!fileInUnderlying->cd(Form("Event%d",iEv))) {
      eventExist = kFALSE;
      break;
    }
    for (Int_t jPlane=0; jPlane<AliMFTConstants::fNMaxPlanes; jPlane++) {
      fRecPointsPerPlaneIn[jPlane] = new TClonesArray("AliMFTCluster");
    }
    TTree *treeInUnderlying = (TTree*) gDirectory->Get("TreeR");
    treeInUnderlying -> SetName("TreeR_InUnderlying");
    iPlane=0;
    while (treeInUnderlying->GetBranch(Form("Plane_%02d",iPlane))) {
      //	printf("Plane %02d\n",iPlane);
      treeInUnderlying ->SetBranchAddress(Form("Plane_%02d",iPlane), &(fRecPointsPerPlaneIn[iPlane]));
      TBranch *branch = treeOutTotal->GetBranch(Form("Plane_%02d",iPlane));
      if (!branch) treeOutTotal->Branch(Form("Plane_%02d",iPlane), &(fRecPointsPerPlaneOut[iPlane])); 
      iPlane++;
    }
    iPlane=0;
    treeInUnderlying -> GetEntry(0);
    while (treeInUnderlying->GetBranch(Form("Plane_%02d",iPlane))) {
      Int_t nClusters = fRecPointsPerPlaneIn[iPlane]->GetEntries();
      for (Int_t iCluster=0; iCluster<nClusters; iCluster++) {
	//	  printf("Cluster %4d\n",iCluster);
	AliMFTCluster *newCluster = (AliMFTCluster*) fRecPointsPerPlaneIn[iPlane]->At(iCluster);
	new ((*fRecPointsPerPlaneOut[iPlane])[fRecPointsPerPlaneOut[iPlane]->GetEntries()]) AliMFTCluster(*newCluster);
      }
      iPlane++;
    }
    for (Int_t jPlane=0; jPlane<AliMFTConstants::fNMaxPlanes; jPlane++) {
      delete fRecPointsPerPlaneIn[jPlane];
    }
    delete treeInUnderlying;
    fileInUnderlying -> Close();
    delete fileInUnderlying;

    // Adding Signal Clusters

    fileInSignal = new TFile("./MFT.RecPoints.Signal.root");
    printf("Merging file %s\n", fileInSignal->GetName());
    if (!fileInSignal->cd(Form("Event%d",iEv))) {
      eventExist = kFALSE;
      break;
    }
    for (Int_t jPlane=0; jPlane<AliMFTConstants::fNMaxPlanes; jPlane++) {
      fRecPointsPerPlaneIn[jPlane] = new TClonesArray("AliMFTCluster");
    }
    TTree *treeInSignal = (TTree*) gDirectory->Get("TreeR");
    treeInSignal -> SetName("TreeR_InSignal");
    iPlane=0;
    while (treeInSignal->GetBranch(Form("Plane_%02d",iPlane))) {
      //	printf("Plane %02d\n",iPlane);
      treeInSignal ->SetBranchAddress(Form("Plane_%02d",iPlane), &(fRecPointsPerPlaneIn[iPlane]));
      TBranch *branch = treeOutTotal->GetBranch(Form("Plane_%02d",iPlane));
      if (!branch) treeOutTotal->Branch(Form("Plane_%02d",iPlane), &(fRecPointsPerPlaneOut[iPlane])); 
      iPlane++;
    }
    iPlane=0;
    treeInSignal -> GetEntry(0);
    while (treeInSignal->GetBranch(Form("Plane_%02d",iPlane))) {
      Int_t nClusters = fRecPointsPerPlaneIn[iPlane]->GetEntries();
      for (Int_t iCluster=0; iCluster<nClusters; iCluster++) {
	//	  printf("Cluster %4d\n",iCluster);
	AliMFTCluster *newCluster = (AliMFTCluster*) fRecPointsPerPlaneIn[iPlane]->At(iCluster);
	new ((*fRecPointsPerPlaneOut[iPlane])[fRecPointsPerPlaneOut[iPlane]->GetEntries()]) AliMFTCluster(*newCluster);
      }
      iPlane++;
    }
    for (Int_t jPlane=0; jPlane<AliMFTConstants::fNMaxPlanes; jPlane++) {
      delete fRecPointsPerPlaneIn[jPlane];
    }
    delete treeInSignal;
    fileInSignal -> Close();
    delete fileInSignal;

    if (eventExist) {
      treeOutTotal -> Fill();
      fileOutTotal -> mkdir(Form("Event%d",iEv));
      fileOutTotal -> cd(Form("Event%d",iEv));
      treeOutTotal -> SetName("TreeR");
      treeOutTotal -> Write();
      delete treeOutTotal;
      fileOutTotal -> Close();
      delete fileOutTotal;
      for (Int_t jPlane=0; jPlane<AliMFTConstants::fNMaxPlanes; jPlane++) {
	delete fRecPointsPerPlaneOut[jPlane];
      }
      iEv++;
    }

  }   

}

//====================================================================================================================================================

//  LocalWords:  nFiles
