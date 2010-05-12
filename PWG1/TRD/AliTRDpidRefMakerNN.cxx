/*************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercialf purposes is hereby granted  *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: AliTRDpidRefMakerNN.cxx 27496 2008-07-22 08:35:45Z cblume $ */

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  Builds the reference tree for the training of neural networks         //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include "TSystem.h"
#include "TDatime.h"
#include "TPDGCode.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TTree.h"
#include "TEventList.h"
#include "TMultiLayerPerceptron.h"

#include "AliPID.h"
#include "AliESDtrack.h"
#include "AliTrackReference.h"

#include "AliTRDtrackV1.h"
#include "AliTRDReconstructor.h"
#include "AliTRDpidUtil.h"
#include "AliTRDpidRefMakerNN.h"
#include "AliTRDpidUtil.h"

#include "Cal/AliTRDCalPID.h"
#include "Cal/AliTRDCalPIDNN.h"
#include "info/AliTRDtrackInfo.h"
#include "info/AliTRDv0Info.h"
#include "info/AliTRDpidInfo.h"

ClassImp(AliTRDpidRefMakerNN)

//________________________________________________________________________
  AliTRDpidRefMakerNN::AliTRDpidRefMakerNN() 
  :AliTRDpidRefMaker()
  ,fNet(NULL)
  ,fTrainMomBin(kAll)
  ,fEpochs(1000)
  ,fMinTrain(100)
  ,fDate(0)
  ,fDoTraining(0)
  ,fContinueTraining(0)
  ,fTrainPath(NULL)
  ,fScale(0)
  ,fLy(0)
  ,fNtrkl(0)
  ,fRef(NULL)
{
  //
  // Default constructor
  //
  SetNameTitle("refMakerNN", "PID(NN) Reference Maker");
}

//________________________________________________________________________
  AliTRDpidRefMakerNN::AliTRDpidRefMakerNN(const char *name) 
  :AliTRDpidRefMaker(name, "PID(NN) Reference Maker")
  ,fNet(NULL)
  ,fTrainMomBin(kAll)
  ,fEpochs(1000)
  ,fMinTrain(100)
  ,fDate(0)
  ,fDoTraining(0)
  ,fContinueTraining(0)
  ,fTrainPath(NULL)
  ,fScale(0)
  ,fLy(0)
  ,fNtrkl(0)
  ,fRef(NULL)
{
  //
  // Default constructor
  //

  memset(fTrain, 0, AliTRDCalPID::kNMom*sizeof(TEventList*));
  memset(fTest, 0, AliTRDCalPID::kNMom*sizeof(TEventList*));
  memset(fTrainData, 0, AliTRDCalPID::kNMom*sizeof(TTree*));

  SetAbundance(.67);
  SetScaledEdx(Float_t(AliTRDCalPIDNN::kMLPscale));
  TDatime datime;
  fDate = datime.GetDate();
}


//________________________________________________________________________
AliTRDpidRefMakerNN::~AliTRDpidRefMakerNN() 
{
}


//________________________________________________________________________
void AliTRDpidRefMakerNN::MakeTrainTestTrees()
{
  // Create output file and tree
  // Called once

  fRef = new TFile("TRD.CalibPIDrefMakerNN.root", "RECREATE");
  for(Int_t ip = 0; ip < AliTRDCalPID::kNMom; ip++){
    fTrainData[ip] = new TTree(Form("fTrainData_%d", ip), Form("NN Reference Data for MomBin %d", ip));
    fTrainData[ip] -> Branch("fdEdx", fdEdx, Form("fdEdx[%d]/F", AliTRDpidUtil::kNNslices));
    fTrainData[ip] -> Branch("fPID", fPID, Form("fPID[%d]/F", AliPID::kSPECIES));
    fTrainData[ip] -> Branch("fLy", &fLy, "fLy/I");
    fTrainData[ip] -> Branch("fNtrkl", &fNtrkl, "fNtrkl/I");
    
    fTrain[ip] = new TEventList(Form("fTrainMom%d", ip), Form("Training list for momentum intervall %d", ip));
    fTest[ip] = new TEventList(Form("fTestMom%d", ip), Form("Test list for momentum intervall %d", ip));
  }
}



//________________________________________________________________________
Bool_t AliTRDpidRefMakerNN::PostProcess()
{
  // Draw result to the screen
  // Called once at the end of the query

  TFile *fCalib = TFile::Open(Form("TRD.CalibPIDrefMaker.root"));
  if (!fCalib) {
    AliError("Calibration file not available");
    return kFALSE;
  }
  fData = (TTree*)fCalib->Get("refMakerNN");
  if (!fData) {
    AliError("Calibration data not available");
    return kFALSE;
  }
  TObjArray *o = NULL;
  if(!(o = (TObjArray*)fCalib->Get("MonitorNN"))) {
    AliWarning("Missing monitoring container.");
    return kFALSE;
  }
  fContainer = (TObjArray*)o->Clone("monitor");

  
  if (!fRef) {
    AliDebug(2, "Loading file TRD.CalibPIDrefMakerNN.root");
    LoadFile("TRD.CalibPIDrefMakerNN.root");
  }
  else AliDebug(2, "file available");

  if (!fRef) {
    MakeTrainTestTrees();
  
    // Convert the CaliPIDrefMaker tree to 11 (different momentum bin) trees for NN training

    LinkPIDdata();
    for(Int_t ip=0; ip < AliTRDCalPID::kNMom; ip++){ 
      for(Int_t is=0; is < AliPID::kSPECIES; is++) {
	memset(fPID, 0, AliPID::kSPECIES*sizeof(Float_t));
	fPID[is] = 1;
	Int_t n(0); // index of data
	for(Int_t itrk = 0; itrk<fData->GetEntries() && n<kMaxStat; itrk++){
	  if(!(fData->GetEntry(itrk))) continue;
	  if(fPIDdataArray->GetPID()!=is) continue;
    fNtrkl = fPIDdataArray->GetNtracklets();
	  for(Int_t ily=fPIDdataArray->GetNtracklets(); ily--;){
	    fLy = ily;
	    if(fPIDdataArray->GetData(ily)->Momentum()!= ip) continue;
	    memset(fdEdx, 0, AliTRDpidUtil::kNNslices*sizeof(Float_t));
	    for(Int_t islice=AliTRDCalPID::kNSlicesNN; islice--;){
	      fdEdx[islice]+=fPIDdataArray->GetData(ily)->fdEdx[islice];
	      fdEdx[islice]/=fScale;
	    }
	    fTrainData[ip] -> Fill();
	    n++;
	  }
	}
	AliDebug(2, Form("%d %d %d", ip, is, n));
      }
    }

  
    fRef -> cd();
    for(Int_t ip = 0; ip < AliTRDCalPID::kNMom; ip++){
      fTrainData[ip] -> Write();
    }

  }
  else AliDebug(2, "file available");


  // build the training and the test list for the neural networks
  for(Int_t ip = 0; ip < AliTRDCalPID::kNMom; ip++){
    MakeTrainingLists(ip);        
  }
  if(!fDoTraining) return kTRUE;



  // train the neural networks
  gSystem->Exec(Form("mkdir ./Networks_%d/",fDate));
  AliDebug(2, Form("TrainMomBin [%d] [%d]", fTrainMomBin, kAll));

  // train single network for a single momentum (recommended)
  if(!(fTrainMomBin == kAll)){
    if(fTrain[fTrainMomBin] -> GetN() < fMinTrain){
      AliError("Warning in AliTRDpidRefMakerNN::PostProcess : Not enough events for training available! Please check Data sample!");
      return kFALSE;
    }
    MakeRefs(fTrainMomBin);
    MonitorTraining(fTrainMomBin);
  }
  // train all momenta
  else{
    for(Int_t iMomBin = 0; iMomBin < AliTRDCalPID::kNMom; iMomBin++){
      if(fTrain[iMomBin] -> GetN() < fMinTrain){
	AliError(Form("Warning in AliTRDpidRefMakerNN::PostProcess : Not enough events for training available for momentum bin [%d]! Please check Data sample!", iMomBin));
	continue;
      }
      MakeRefs(fTrainMomBin);
      MonitorTraining(iMomBin);
    }
  }

  return kTRUE; // testing protection
}


//________________________________________________________________________
void AliTRDpidRefMakerNN::MakeTrainingLists(Int_t mombin) 
{
  //
  // build the training lists for the neural networks
  //

  if (!fRef) {
    LoadFile(Form("TRD.Calib%s.root", GetName()));
  }

  if (!fRef) {
    AliError("ERROR file for building training list not available");
    return;
  }

  AliDebug(2, "\n Making training lists! \n");

  Int_t nPart[AliPID::kSPECIES];
  memset(nPart, 0, AliPID::kSPECIES*sizeof(Int_t));

  // set needed branches
  fTrainData[mombin] -> SetBranchAddress("fdEdx", fdEdx);
  fTrainData[mombin] -> SetBranchAddress("fPID", fPID);
  fTrainData[mombin] -> SetBranchAddress("fLy", &fLy);
  fTrainData[mombin] -> SetBranchAddress("fNtrkl", &fNtrkl);
  
  // start first loop to check total number of each particle type
  for(Int_t iEv=0; iEv < fTrainData[mombin] -> GetEntries(); iEv++){
    fTrainData[mombin] -> GetEntry(iEv);

    // use only events with goes through 6 layers TRD
    if(fNtrkl != AliTRDgeometry::kNlayer) continue;

    if(fPID[AliPID::kElectron] == 1)
      nPart[AliPID::kElectron]++;
    else if(fPID[AliPID::kMuon] == 1)
      nPart[AliPID::kMuon]++;
    else if(fPID[AliPID::kPion] == 1)
      nPart[AliPID::kPion]++;
    else if(fPID[AliPID::kKaon] == 1)
      nPart[AliPID::kKaon]++;
    else if(fPID[AliPID::kProton] == 1)
      nPart[AliPID::kProton]++;
  }

  AliDebug(2, "Particle multiplicities:");
  AliDebug(2, Form("Momentum[%d]  Elecs[%d] Muons[%d] Pions[%d] Kaons[%d] Protons[%d]", mombin, nPart[AliPID::kElectron], nPart[AliPID::kMuon], nPart[AliPID::kPion], nPart[AliPID::kKaon], nPart[AliPID::kProton]));

  

  //   // implement counter of training and test sample size
  Int_t iTrain = 0, iTest = 0;

  // set training sample size per momentum interval to 2/3 
  // of smallest particle counter and test sample to 1/3
  iTrain = nPart[0];
  for(Int_t iPart = 1; iPart < AliPID::kSPECIES; iPart++){
    // exclude muons and kaons if not availyable
    // this is neeeded since we do not have v0 candiates
    if((iPart == AliPID::kMuon || iPart == AliPID::kKaon) && (nPart[AliPID::kMuon] == 0 || nPart[AliPID::kKaon] == 0)) continue;
    if(iTrain > nPart[iPart])
      iTrain = nPart[iPart];
  } 
  iTest = Int_t( iTrain * (1-fFreq));
  iTrain = Int_t(iTrain * fFreq);
  AliDebug(2, Form("Momentum[%d]  Train[%d] Test[%d]", mombin, iTrain, iTest));



  // reset couters
  memset(nPart, 0, AliPID::kSPECIES*sizeof(Int_t));

  // start second loop to set the event lists
  for(Int_t iEv = 0; iEv < fTrainData[mombin] -> GetEntries(); iEv++){
    fTrainData[mombin] -> GetEntry(iEv);

    // set event list
    for(Int_t is = 0; is < AliPID::kSPECIES; is++){
      if(nPart[is] < iTrain && fPID[is] == 1){
       fTrain[mombin] -> Enter(iEv);
        nPart[is]++;
      } else if(nPart[is] < iTest+iTrain && fPID[is] == 1){
        fTest[mombin] -> Enter(iEv);
        nPart[is]++;
      } else continue;
    }
  }
  
  AliDebug(2, "Particle multiplicities in both lists:");
  AliDebug(2, Form("Momentum[%d]  Elecs[%d] Muons[%d] Pions[%d] Kaons[%d] Protons[%d]", mombin, nPart[AliPID::kElectron], nPart[AliPID::kMuon], nPart[AliPID::kPion], nPart[AliPID::kKaon], nPart[AliPID::kProton]));

  return;
}


//________________________________________________________________________
void AliTRDpidRefMakerNN::MakeRefs(Int_t mombin) 
{
  //
  // train the neural networks
  //
  
  
  if (!fTrainData[mombin]) LoadFile(Form("TRD.Calib%s.root", GetName()));

  if (!fTrainData[mombin]) {
    AliError("Tree for training list not available");
    return;
  }

  TDatime datime;
  fDate = datime.GetDate();

  AliDebug(2, Form("Training momentum bin %d", mombin));

  // set variable to monitor the training and to save the development of the networks
  Int_t nEpochs = fEpochs/kMoniTrain;       
  AliDebug(2, Form("Training %d times %d epochs", kMoniTrain, nEpochs));

  // make directories to save the networks 
  gSystem->Exec(Form("rm -r ./Networks_%d/MomBin_%d",fDate, mombin));
  gSystem->Exec(Form("mkdir ./Networks_%d/MomBin_%d",fDate, mombin));

  // variable to check if network can load weights from previous training
  Bool_t bFirstLoop = kTRUE;

  // train networks over several loops and save them after each loop
  for(Int_t iLoop = 0; iLoop < kMoniTrain; iLoop++){
    fTrainData[mombin] -> SetEventList(fTrain[mombin]);
    fTrainData[mombin] -> SetEventList(fTest[mombin]);
      
    AliDebug(2, Form("Momentum[%d] Trainingloop[%d]", mombin, iLoop));
      
    // check if network is already implemented
    if(bFirstLoop == kTRUE){
      fNet = new TMultiLayerPerceptron("fdEdx[0],fdEdx[1],fdEdx[2],fdEdx[3],fdEdx[4],fdEdx[5],fdEdx[6],fdEdx[7]:15:7:fPID[0],fPID[1],fPID[2],fPID[3],fPID[4]!",fTrainData[mombin],fTrain[mombin],fTest[mombin]);
      fNet -> SetLearningMethod(TMultiLayerPerceptron::kStochastic);       // set learning method
      fNet -> TMultiLayerPerceptron::SetEta(0.001);                        // set learning speed
      if(!fContinueTraining){
	if(AliLog::GetDebugLevel("","AliTRDpidRefMakerNN")>=2) fNet -> Train(nEpochs,"text update=10, graph");
	else fNet -> Train(nEpochs,"");
      }
      else{
	fNet -> LoadWeights(Form("./Networks_%d/MomBin_%d/Net_%d",fTrainPath, mombin, kMoniTrain - 1));
	if(AliLog::GetDebugLevel("","AliTRDpidRefMakerNN")>=2) fNet -> Train(nEpochs,"text update=10, graph+");      
	else fNet -> Train(nEpochs,"+");                   
      }
      bFirstLoop = kFALSE;
    }
    else{    
      if(AliLog::GetDebugLevel("","AliTRDpidRefMakerNN")>=2) fNet -> Train(nEpochs,"text update=10, graph+");      
      else fNet -> Train(nEpochs,"+");                   
    }
      
    // save weights for monitoring of the training
    fNet -> DumpWeights(Form("./Networks_%d/MomBin_%d/Net_%d",fDate, mombin, iLoop));
  }   // end training loop
}



//________________________________________________________________________
void AliTRDpidRefMakerNN::MonitorTraining(Int_t mombin) 
{
  //
  // train the neural networks
  //
  
  if (!fTrainData[mombin]) LoadFile(Form("TRD.Calib%s.root", GetName()));
  if (!fTrainData[mombin]) {
    AliError("Tree for training list not available");
    return;
  }

  // init networks and set event list
  for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
  fNet = new TMultiLayerPerceptron("fdEdx[0],fdEdx[1],fdEdx[2],fdEdx[3],fdEdx[4],fdEdx[5],fdEdx[6],fdEdx[7]:15:7:fPID[0],fPID[1],fPID[2],fPID[3],fPID[4]!",fTrainData[mombin],fTrain[mombin],fTest[mombin]);   
  fTrainData[mombin] -> SetEventList(fTrain[mombin]);
  fTrainData[mombin] -> SetEventList(fTest[mombin]);
  }

  // implement variables for likelihoods
  Float_t like[AliPID::kSPECIES][AliTRDgeometry::kNlayer];
  memset(like, 0, AliPID::kSPECIES*AliTRDgeometry::kNlayer*sizeof(Float_t));
  Float_t likeAll[AliPID::kSPECIES], totProb;

  Double_t pionEffiTrain[kMoniTrain], pionEffiErrTrain[kMoniTrain];
  Double_t pionEffiTest[kMoniTrain], pionEffiErrTest[kMoniTrain];
  memset(pionEffiTrain, 0, kMoniTrain*sizeof(Double_t));
  memset(pionEffiErrTrain, 0, kMoniTrain*sizeof(Double_t));
  memset(pionEffiTest, 0, kMoniTrain*sizeof(Double_t));
  memset(pionEffiErrTest, 0, kMoniTrain*sizeof(Double_t));

  // init histos
  const Float_t epsilon = 1/(2*(AliTRDpidUtil::kBins-1));     // get nice histos with bin center at 0 and 1
  TH1F *hElecs, *hPions;
  hElecs = new TH1F("hElecs","Likelihood for electrons", AliTRDpidUtil::kBins, 0.-epsilon, 1.+epsilon);
  hPions = new TH1F("hPions","Likelihood for pions", AliTRDpidUtil::kBins, 0.-epsilon, 1.+epsilon);

  TGraphErrors *gEffisTrain = new TGraphErrors(kMoniTrain);
  gEffisTrain -> SetLineColor(4);
  gEffisTrain -> SetMarkerColor(4);
  gEffisTrain -> SetMarkerStyle(29);
  gEffisTrain -> SetMarkerSize(1);

  TGraphErrors *gEffisTest = new TGraphErrors(kMoniTrain);
  gEffisTest -> SetLineColor(2);
  gEffisTest -> SetMarkerColor(2);
  gEffisTest -> SetMarkerStyle(29);
  gEffisTest -> SetMarkerSize(1);

  AliTRDpidUtil *util = new AliTRDpidUtil();
  
  // monitor training progress
  for(Int_t iLoop = 0; iLoop < kMoniTrain; iLoop++){

    // load weights
    fNet -> LoadWeights(Form("./Networks_%d/MomBin_%d/Net_%d",fDate, mombin, iLoop));
    AliDebug(2, Form("./Networks_%d/MomBin_%d/Net_%d",fDate, mombin, iLoop));

    // event loop training list

    for(Int_t is = 0; is < AliPID::kSPECIES; is++){

      if(!((is == AliPID::kElectron) || (is == AliPID::kPion))) continue;

      Int_t iChamb = 0;
      // reset particle probabilities
      for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
	likeAll[iPart] = 1./AliPID::kSPECIES;
      }
      totProb = 0.;

      AliDebug(2, Form("%d",fTrain[mombin] -> GetN()));
      for(Int_t iEvent = 0; iEvent < fTrain[mombin] -> GetN(); iEvent++ ){
	fTrainData[mombin] -> GetEntry(fTrain[mombin] -> GetEntry(iEvent));
	// use event only if it is electron or pion
	if(!(fPID[is] == 1.0)) continue;
	// get the probabilities for each particle type in each chamber
	for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
	  like[iPart][iChamb] = fNet -> Result(fTrain[mombin] -> GetEntry(iEvent), iPart);
	  likeAll[iPart] *=  like[iPart][iChamb];
	}
	//end chamber loop
	iChamb++;

	// get total probability and normalize it
	for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
	  totProb += likeAll[iPart];
	}
	for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
	  likeAll[iPart] /= totProb;
	}

	if(iChamb == 5){
	  // fill likelihood distributions
	  if(fPID[AliPID::kElectron] == 1)      
	    hElecs -> Fill(likeAll[AliPID::kElectron]);
	  if(fPID[AliPID::kPion] == 1)      
	    hPions -> Fill(likeAll[AliPID::kElectron]);
	  iChamb = 0;
	  // reset particle probabilities
	  for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
	    likeAll[iPart] = 1./AliPID::kSPECIES;
	  }
	  totProb = 0.;
	}
      } // end event loop
    } // end species loop
    

    // calculate the pion efficiency and fill the graph
    util -> CalculatePionEffi(hElecs, hPions);
    pionEffiTrain[iLoop] = util -> GetPionEfficiency();
    pionEffiErrTrain[iLoop] = util -> GetError();

    gEffisTrain -> SetPoint(iLoop, iLoop+1, pionEffiTrain[iLoop]);
    gEffisTrain -> SetPointError(iLoop, 0, pionEffiErrTrain[iLoop]);
    hElecs -> Reset();
    hPions -> Reset();
    AliDebug(2, Form("TrainingLoop[%d] PionEfficiency[%f +/- %f]", iLoop, pionEffiTrain[iLoop], pionEffiErrTrain[iLoop]));
    // end training loop
    


    // monitor validation progress
    for(Int_t is = 0; is < AliPID::kSPECIES; is++){

      if(!((is == AliPID::kElectron) || (is == AliPID::kPion))) continue;

      Int_t iChamb = 0;
      // reset particle probabilities
      for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
	likeAll[iPart] = 1./AliPID::kSPECIES;
      }
      totProb = 0.;

      for(Int_t iEvent = 0; iEvent < fTest[mombin] -> GetN(); iEvent++ ){
	fTrainData[mombin] -> GetEntry(fTest[mombin] -> GetEntry(iEvent));
	// use event only if it is electron or pion
	if(!(fPID[is] == 1.0)) continue;

	// get the probabilities for each particle type in each chamber
	for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
	  like[iPart][iChamb] = fNet -> Result(fTest[mombin] -> GetEntry(iEvent), iPart);
	  likeAll[iPart] *=  like[iPart][iChamb];
	}
	iChamb++;

	// get total probability and normalize it
	for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
	  totProb += likeAll[iPart];
	}
	for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
	  likeAll[iPart] /= totProb;
	}

	if(iChamb == 5){
	  // fill likelihood distributions
	  if(fPID[AliPID::kElectron] == 1)      
	    hElecs -> Fill(likeAll[AliPID::kElectron]);
	  if(fPID[AliPID::kPion] == 1)      
	    hPions -> Fill(likeAll[AliPID::kElectron]);
	  iChamb = 0;
	  // reset particle probabilities
	  for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
	    likeAll[iPart] = 1./AliPID::kSPECIES;
	  }
	  totProb = 0.;
	}
      } // end event loop
    } // end species loop
    
    // calculate the pion efficiency and fill the graph
    util -> CalculatePionEffi(hElecs, hPions);
    pionEffiTest[iLoop] = util -> GetPionEfficiency();
    pionEffiErrTest[iLoop] = util -> GetError();

    gEffisTest -> SetPoint(iLoop, iLoop+1, pionEffiTest[iLoop]);
    gEffisTest -> SetPointError(iLoop, 0, pionEffiErrTest[iLoop]);
    hElecs -> Reset();
    hPions -> Reset();
    AliDebug(2, Form("TestLoop[%d] PionEfficiency[%f +/- %f] \n", iLoop, pionEffiTest[iLoop], pionEffiErrTest[iLoop]));
    
  } //   end validation loop

  util -> Delete();

  gEffisTest -> Draw("PAL");
  gEffisTrain -> Draw("PL");

}


//________________________________________________________________________
Bool_t AliTRDpidRefMakerNN::LoadFile(const Char_t *InFileNN) 
{
  //
  // Loads the files and sets the event list
  // for neural network training.
  // Useable for training outside of the makeResults.C macro
  //

  fRef = TFile::Open(Form("%s", InFileNN));
  if(!fRef) return 0;
  for(Int_t ip = 0; ip < AliTRDCalPID::kNMom; ip++){
    fTrainData[ip] = (TTree*)fRef -> Get(Form("fTrainData_%d", ip));
  }

  for(Int_t iMom = 0; iMom < AliTRDCalPID::kNMom; iMom++){
      fTrain[iMom] = new TEventList(Form("fTrainMom%d", iMom), Form("Training list for momentum intervall %d", iMom));
      fTest[iMom] = new TEventList(Form("fTestMom%d", iMom), Form("Test list for momentum intervall %d", iMom));
  }
  return 1;
}


