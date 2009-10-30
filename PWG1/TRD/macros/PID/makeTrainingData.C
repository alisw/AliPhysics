// This macro splits the file created by the AliTRDpidRefMaker to smaller training 
// samples according to the momentum of the track in order to increase the training 
// speed for the TMultiLayerPerceptrons.
// The procedure is the following:
// 1. Create a directory where the training should take place
// 2. Copy the TRD.TaskPidRefMakerNN.root and the TRD.TaskPidRefMakerLQ.root in this directory.
// 3. Run makeTrainingDataNN(). This creates new directories: 0.6GeV, 0.8GeV, ..., 10.0GeV and 
//    create a subset of the training data according to the momentum.
// 4. Run makeDataLQ(). Does the same as make TraiingDataNN for the LQ data without the creation
//    of the directories.
// 5. Run CreateDummy() to create a TRD.TaskPidRefMaker.root file. This is necessary for the 
//    monitoring of the training progress.
// 6. Go to the subdirectories and run the training.

#ifndef __CINT__
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TObjArray.h"
#include "TGraphErrors.h"

#include "AliPID.h"
#include "AliTRDpidUtil.h"
#include "Cal/AliTRDCalPID.h"
#include "qaRec/AliTRDpidRefMasker.h"
#endif

Int_t makeTrainingDataNN(){

  Int_t fLayer;
  Float_t fMom, fv0pid[AliPID::kSPECIES], fdEdx[AliTRDpidUtil::kNNslices];

  AliTRDpidUtil Util;
  TFile *fIn = new TFile("../TRD.TaskPidRefMakerNN.root","READ");
  TTree *tIn = (TTree*) fIn -> Get("NN");

  tIn -> SetBranchAddress("fLayer", &fLayer);
  tIn -> SetBranchAddress("fMom", &fMom);
  tIn -> SetBranchAddress("fv0pid", fv0pid);
  tIn -> SetBranchAddress("fdEdx", fdEdx);

  for(Int_t iMomBin = 0; iMomBin < AliTRDCalPID::kNMom; iMomBin++){

    gSystem->Exec(Form("mkdir -v ./%3.1fGeV",AliTRDCalPID::GetMomentum(iMomBin)));

    printf("Extracting training set for momentum bin %3.1f GeV/c\n", AliTRDCalPID::GetMomentum(iMomBin));
    TFile *fOut = new TFile(Form("./%3.1fGeV/TRD.TaskPidRefMakerNN.root",AliTRDCalPID::GetMomentum(iMomBin)),"RECREATE");
    TTree *tOut = new TTree("NN", "Reference data for NN"); 
    tOut -> Branch("fLayer", &fLayer, "fLayer/I");
    tOut -> Branch("fMom", &fMom, "fMom/F");
    tOut -> Branch("fv0pid", fv0pid, Form("fv0pid[%d]/F",AliPID::kSPECIES));
    tOut -> Branch("fdEdx", fdEdx, Form("fdEdx[%d]/F", AliTRDpidUtil::kNNslices));
    
    for(Int_t iEv = 0; iEv < (tIn -> GetEntries()); iEv++){
      fLayer = 0;
      fMom = 0.0;
      for(Int_t i = 0; i < AliPID::kSPECIES; i++) fv0pid[i] = 0.0;
      for(Int_t i = 0; i < AliTRDpidUtil::kNNslices; i++) fdEdx[i] = 0.0;
      tIn -> GetEntry(iEv);
      if(Util.GetMomentumBin(fMom) == iMomBin){
	tOut -> Fill();
      }
    }
    
    tOut -> Write();
    tOut -> Delete();
    fOut -> Close();

  }
  
  printf("Extraction completed!");
  return 1;

}


Int_t makeDataLQ(){

  Int_t fLayer;
  Float_t fMom, fv0pid[AliPID::kSPECIES], fdEdx[AliTRDpidUtil::kLQslices];

  AliTRDpidUtil Util;
  TFile *fIn = new TFile("../TRD.TaskPidRefMakerLQ.root","READ");
  TTree *tIn = (TTree*) fIn -> Get("LQ");

  tIn -> SetBranchAddress("fLayer", &fLayer);
  tIn -> SetBranchAddress("fMom", &fMom);
  tIn -> SetBranchAddress("fv0pid", fv0pid);
  tIn -> SetBranchAddress("fdEdx", fdEdx);

  for(Int_t iMomBin = 0; iMomBin < AliTRDCalPID::kNMom; iMomBin++){
    printf("Extracting training set for momentum bin %3.1f GeV/c\n", AliTRDCalPID::GetMomentum(iMomBin));
    TFile *fOut = new TFile(Form("./%3.1fGeV/TRD.TaskPidRefMakerLQ.root",AliTRDCalPID::GetMomentum(iMomBin)),"RECREATE");
    TTree *tOut = new TTree("LQ", "Reference data for LQ"); 
    tOut -> Branch("fLayer", &fLayer, "fLayer/I");
    tOut -> Branch("fMom", &fMom, "fMom/F");
    tOut -> Branch("fv0pid", fv0pid, Form("fv0pid[%d]/F",AliPID::kSPECIES));
    tOut -> Branch("fdEdx", fdEdx, Form("fdEdx[%d]/F", AliTRDpidUtil::kLQslices));
    
    for(Int_t iEv = 0; iEv < (tIn -> GetEntries()); iEv++){
      fLayer = 0;
      fMom = 0.0;
      for(Int_t i = 0; i < AliPID::kSPECIES; i++) fv0pid[i] = 0.0;
      for(Int_t i = 0; i < AliTRDpidUtil::kLQslices; i++) fdEdx[i] = 0.0;
      tIn -> GetEntry(iEv);
      if(Util.GetMomentumBin(fMom) == iMomBin){
	tOut -> Fill();
      }
    }
    
    tOut -> Write();
    tOut -> Delete();
    fOut -> Close();

  }
  
  printf("Extraction completed!");
  return 1;

}


Int_t CreateDummy(){

  for(Int_t iMomBin = 0; iMomBin < AliTRDCalPID::kNMom; iMomBin++){
    printf("Creating dummy for momentum bin %3.1f GeV/c\n", AliTRDCalPID::GetMomentum(iMomBin));
    TFile *fOut = new TFile(Form("./%3.1fGeV/TRD.TaskPidRefMaker.root",AliTRDCalPID::GetMomentum(iMomBin)),"RECREATE");
    TObjArray *fContainer = new TObjArray();

    TGraphErrors *gEffisTrain = new TGraphErrors(50);
//     TGraphErrors *gEffisTrain = new TGraphErrors(AliTRDpidRefMaker::kMoniTrain);
    gEffisTrain -> SetLineColor(4);
    gEffisTrain -> SetMarkerColor(4);
    gEffisTrain -> SetMarkerStyle(29);
    gEffisTrain -> SetMarkerSize(1);

    TGraphErrors *gEffisTest = new TGraphErrors(50);
//     TGraphErrors *gEffisTest = new TGraphErrors(AliTRDpidRefMaker::kMoniTrain);
    gEffisTest -> SetLineColor(2);
    gEffisTest -> SetMarkerColor(2);
    gEffisTest -> SetMarkerStyle(29);
    gEffisTest -> SetMarkerSize(1);
    
    fContainer -> AddAt(gEffisTrain,0);
    fContainer -> AddAt(gEffisTest,1);
//     fContainer -> AddAt(gEffisTrain,AliTRDpidRefMaker::kGraphTrain);
//     fContainer -> AddAt(gEffisTest,AliTRDpidRefMaker::kGraphTest);
  
  }

}
