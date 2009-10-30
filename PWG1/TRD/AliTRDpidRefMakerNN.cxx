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

ClassImp(AliTRDpidRefMakerNN)

//________________________________________________________________________
AliTRDpidRefMakerNN::AliTRDpidRefMakerNN() 
  :AliTRDpidRefMaker("PidRefMakerNN", "PID(NN) Reference Maker")
//   :AliTRDrecoTask("PidRefMakerNN", "PID(NN) Reference Maker")
  ,fTrainMomBin(kAll)
  ,fEpochs(1000)
  ,fMinTrain(100)
  ,fDate(0)
  ,fDoTraining(0)
  ,fContinueTraining(0)
  ,fTrainPath(0x0)
  ,fScale(0)
{
  //
  // Default constructor
  //

  memset(fTrain, 0, AliTRDCalPID::kNMom*AliTRDgeometry::kNlayer*sizeof(TEventList*));
  memset(fTest, 0, AliTRDCalPID::kNMom*AliTRDgeometry::kNlayer*sizeof(TEventList*));
  memset(fNet, 0, AliTRDgeometry::kNlayer*sizeof(TMultiLayerPerceptron*));

  SetAbundance(.67);
  SetScaledEdx(Float_t(AliTRDCalPIDNN::kMLPscale));
  TDatime datime;
  fDate = datime.GetDate();

  DefineInput(1, TObjArray::Class());
  DefineOutput(1, TTree::Class());
}


//________________________________________________________________________
AliTRDpidRefMakerNN::~AliTRDpidRefMakerNN() 
{
}


//________________________________________________________________________
void AliTRDpidRefMakerNN::CreateOutputObjects()
{
  // Create histograms
  // Called once

  AliTRDpidRefMaker::CreateOutputObjects();
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

  fContainer -> AddAt(gEffisTrain,kGraphTrain);
  fContainer -> AddAt(gEffisTest,kGraphTest);
}



//________________________________________________________________________
Bool_t AliTRDpidRefMakerNN::PostProcess()
{
  // Draw result to the screen
  // Called once at the end of the query

  // build the training andthe test list for the neural networks
  MakeTrainingLists();        
  if(!fDoTraining) return kTRUE;

  // train the neural networks and build the refrence histos for 2-dim LQ
  gSystem->Exec(Form("mkdir ./Networks_%d/",fDate));
  AliDebug(2, Form("TrainMomBin [%d] [%d]", fTrainMomBin, kAll));

  // train single network for a single momentum (recommended)
  if(!(fTrainMomBin == kAll)){
    if(fTrain[fTrainMomBin][0] -> GetN() < fMinTrain){
      AliDebug(2, Form("Warning in AliTRDpidRefMakerNN::PostProcess : Not enough events for training available! Please check Data sample!"));
      return kFALSE;
    }
    MakeRefs(fTrainMomBin);
//     TrainNetworks(fTrainMomBin);
    MonitorTraining(fTrainMomBin);
  }
  // train all momenta
  else{
    for(Int_t iMomBin = 0; iMomBin < AliTRDCalPID::kNMom; iMomBin++){
      if(fTrain[iMomBin][0] -> GetN() < fMinTrain){
        AliDebug(2, Form("Warning in AliTRDpidRefMakerNN::PostProcess : Not enough events for training available for momentum bin [%d]! Please check Data sample!", iMomBin));
  continue;
      }
      MakeRefs(fTrainMomBin);
//       TrainNetworks(iMomBin);
      MonitorTraining(iMomBin);
    }
  }

  return kTRUE; // testing protection
}


//________________________________________________________________________
void AliTRDpidRefMakerNN::MakeTrainingLists() 
{
  //
  // build the training lists for the neural networks
  //

  if (!fData) {
    LoadFile("TRD.CalibPidRefMakerNN.root");
  }

  if (!fData) {
    Printf("ERROR tree for training list not available");
    return;
  }

  AliDebug(2, "  Making training lists! \n");

  Int_t nPart[AliPID::kSPECIES][AliTRDCalPID::kNMom];
  memset(nPart, 0, AliPID::kSPECIES*AliTRDCalPID::kNMom*sizeof(Int_t));

  // set needed branches
  LinkPIDdata();

  // start first loop to check total number of each particle type
  for(Int_t iEv=0; iEv < fData -> GetEntries(); iEv++){
    fData -> GetEntry(iEv);

    // use only events with goes through 6 layers TRD
    if(fPIDdataArray->fNtracklets != AliTRDgeometry::kNlayer) continue;

    for(Int_t ily=AliTRDgeometry::kNlayer; ily--;) nPart[fPIDbin][fPIDdataArray->fData[ily].fPLbin & 0xf]++;
  }

  AliDebug(2, "Particle multiplicities:");
  for(Int_t iMomBin = 0; iMomBin <AliTRDCalPID::kNMom; iMomBin++)
    AliDebug(2, Form("Momentum[%d]  Elecs[%d] Muons[%d] Pions[%d] Kaons[%d] Protons[%d]", iMomBin, nPart[AliPID::kElectron][iMomBin], nPart[AliPID::kMuon][iMomBin], nPart[AliPID::kPion][iMomBin], nPart[AliPID::kKaon][iMomBin], nPart[AliPID::kProton][iMomBin]));

  // implement counter of training and test sample size
  Int_t iTrain[AliTRDCalPID::kNMom], iTest[AliTRDCalPID::kNMom];
  memset(iTrain, 0, AliTRDCalPID::kNMom*sizeof(Int_t));
  memset(iTest, 0, AliTRDCalPID::kNMom*sizeof(Int_t));

  // set training sample size per momentum interval to 2/3 
  // of smallest particle counter and test sample to 1/3
  for(Int_t iMomBin = 0; iMomBin < AliTRDCalPID::kNMom; iMomBin++){
    iTrain[iMomBin] = nPart[0][iMomBin];
    for(Int_t iPart = 1; iPart < AliPID::kSPECIES; iPart++){
      if(iTrain[iMomBin] > nPart[iPart][iMomBin])
  iTrain[iMomBin] = nPart[iPart][iMomBin];
    } 
    iTrain[iMomBin] = Int_t(iTrain[iMomBin] * .66);
    iTest[iMomBin] = Int_t( iTrain[iMomBin] * .5);
    AliDebug(2, Form("Momentum[%d]  Train[%d] Test[%d]", iMomBin, iTrain[iMomBin], iTest[iMomBin]));
  }


  // reset couters
  memset(nPart, 0, AliPID::kSPECIES*AliTRDCalPID::kNMom*sizeof(Int_t));

  // start second loop to set the event lists
  for(Int_t iEv = 0; iEv < fData -> GetEntries(); iEv++){
    fData -> GetEntry(iEv);

    // use only events with goes through 6 layers TRD
    if(fPIDdataArray->fNtracklets != AliTRDgeometry::kNlayer) continue;

    for(Int_t ily=AliTRDgeometry::kNlayer; ily--;){ 
      Int_t iMomBin = fPIDdataArray->fData[ily].fPLbin & 0xf;
    
      // set event list
      if(nPart[fPIDbin][iMomBin] < iTrain[iMomBin]){
        fTrain[iMomBin][ily] -> Enter(iEv + ily);
        nPart[fPIDbin][iMomBin]++;
      } else if(nPart[fPIDbin][iMomBin] < iTest[iMomBin]+iTrain[iMomBin]){
        fTest[iMomBin][ily] -> Enter(iEv + ily);
        nPart[fPIDbin][iMomBin]++;
      } else continue;
    }
  }
  
  AliDebug(2, "Particle multiplicities in both lists:");
  for(Int_t iMomBin = 0; iMomBin <AliTRDCalPID::kNMom; iMomBin++)
    AliDebug(2, Form("Momentum[%d]  Elecs[%d] Muons[%d] Pions[%d] Kaons[%d] Protons[%d]", iMomBin, nPart[AliPID::kElectron][iMomBin], nPart[AliPID::kMuon][iMomBin], nPart[AliPID::kPion][iMomBin], nPart[AliPID::kKaon][iMomBin], nPart[AliPID::kProton][iMomBin]));
}


//________________________________________________________________________
void AliTRDpidRefMakerNN::MakeRefs(Int_t mombin) 
{
  //
  // train the neural networks
  //
  
  
  if (!fData) LoadFile(Form("TRD.Calib%s.root", GetName()));

  if (!fData) {
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
  Bool_t bFirstLoop[AliTRDgeometry::kNlayer];
  memset(bFirstLoop, kTRUE, AliTRDgeometry::kNlayer*sizeof(Bool_t));

  // train networks over several loops and save them after each loop
  for(Int_t iLoop = 0; iLoop < kMoniTrain; iLoop++){
    // loop over chambers
    for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
      // set the event lists
      fData -> SetEventList(fTrain[mombin][iChamb]);
      fData -> SetEventList(fTest[mombin][iChamb]);
      
      AliDebug(2, Form("Trainingloop[%d] Chamber[%d]", iLoop, iChamb));
      
      // check if network is already implemented
      if(bFirstLoop[iChamb] == kTRUE){
  fNet[iChamb] = new TMultiLayerPerceptron("fdEdx[0],fdEdx[1],fdEdx[2],fdEdx[3],fdEdx[4],fdEdx[5],fdEdx[6],fdEdx[7]:15:7:fPID[0],fPID[1],fPID[2],fPID[3],fPID[4]!",fData,fTrain[mombin][iChamb],fTest[mombin][iChamb]);
  fNet[iChamb] -> SetLearningMethod(TMultiLayerPerceptron::kStochastic);       // set learning method
  fNet[iChamb] -> TMultiLayerPerceptron::SetEta(0.001);                        // set learning speed
  if(!fContinueTraining){
    if(DebugLevel()>=2) fNet[iChamb] -> Train(nEpochs,"text update=10, graph");
    else fNet[iChamb] -> Train(nEpochs,"");
  }
  else{
    fNet[iChamb] -> LoadWeights(Form("./Networks_%d/MomBin_%d/Net%d_%d",fTrainPath, mombin, iChamb, kMoniTrain - 1));
    if(DebugLevel()>=2) fNet[iChamb] -> Train(nEpochs,"text update=10, graph+");      
    else fNet[iChamb] -> Train(nEpochs,"+");                   
  }
  bFirstLoop[iChamb] = kFALSE;
      }
      else{    
  if(DebugLevel()>=2) fNet[iChamb] -> Train(nEpochs,"text update=10, graph+");      
  else fNet[iChamb] -> Train(nEpochs,"+");                   
      }
      
      // save weights for monitoring of the training
      fNet[iChamb] -> DumpWeights(Form("./Networks_%d/MomBin_%d/Net%d_%d",fDate, mombin, iChamb, iLoop));
    } // end chamber loop
  }   // end training loop
}



//________________________________________________________________________
void AliTRDpidRefMakerNN::MonitorTraining(Int_t mombin) 
{
  //
  // train the neural networks
  //
  
  if (!fData) LoadFile(Form("TRD.Calib%s.root", GetName()));
  if (!fData) {
    AliError("Tree for training list not available");
    return;
  }

  // init networks and set event list
  for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
  fNet[iChamb] = new TMultiLayerPerceptron("fdEdx[0],fdEdx[1],fdEdx[2],fdEdx[3],fdEdx[4],fdEdx[5],fdEdx[6],fdEdx[7]:15:7:fPID[0],fPID[1],fPID[2],fPID[3],fPID[4]!",fData,fTrain[mombin][iChamb],fTest[mombin][iChamb]);   
  fData -> SetEventList(fTrain[mombin][iChamb]);
  fData -> SetEventList(fTest[mombin][iChamb]);
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

  TGraphErrors *gEffisTrain=0x0, *gEffisTest=0x0;
  gEffisTrain = (TGraphErrors*)fContainer->At(kGraphTrain);
  gEffisTest = (TGraphErrors*)fContainer->At(kGraphTest);

  AliTRDpidUtil *util = new AliTRDpidUtil();
  
  // monitor training progress
  for(Int_t iLoop = 0; iLoop < kMoniTrain; iLoop++){

    // load weights
    for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
      fNet[iChamb] -> LoadWeights(Form("./Networks_%d/MomBin_%d/Net%d_%d",fDate, mombin, iChamb, iLoop));
    }

    // event loop training list
    for(Int_t iEvent = 0; iEvent < fTrain[mombin][0] -> GetN(); iEvent++ ){

      // reset particle probabilities
      for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
  likeAll[iPart] = 1./AliPID::kSPECIES;
      }
      totProb = 0.;

      fData -> GetEntry(fTrain[mombin][0] -> GetEntry(iEvent));
      // use event only if it is electron or pion
      if(!((fPID[AliPID::kElectron] == 1.0) || (fPID[AliPID::kPion] == 1.0))) continue;

      // get the probabilities for each particle type in each chamber
      for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
  for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
    like[iPart][iChamb] = fNet[iChamb] -> Result(fTrain[mombin][iChamb] -> GetEntry(iEvent), iPart);
    likeAll[iPart] *=  like[iPart][iChamb];
  }
      }

      // get total probability and normalize it
      for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
  totProb += likeAll[iPart];
      }
      for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
  likeAll[iPart] /= totProb;
      }

      // fill likelihood distributions
      if(fPID[AliPID::kElectron] == 1)      
  hElecs -> Fill(likeAll[AliPID::kElectron]);
      if(fPID[AliPID::kPion] == 1)      
  hPions -> Fill(likeAll[AliPID::kElectron]);
    } // end event loop


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
    


    // event loop test list
    for(Int_t iEvent = 0; iEvent < fTest[mombin][0] -> GetN(); iEvent++ ){

      // reset particle probabilities
      for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
  likeAll[iPart] = 1./AliTRDgeometry::kNlayer;
      }
      totProb = 0.;

      fData -> GetEntry(fTest[mombin][0] -> GetEntry(iEvent));
      // use event only if it is electron or pion
      if(!((fPID[AliPID::kElectron] == 1.0) || (fPID[AliPID::kPion] == 1.0))) continue;

      // get the probabilities for each particle type in each chamber
      for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
  for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
    like[iPart][iChamb] = fNet[iChamb] -> Result(fTest[mombin][iChamb] -> GetEntry(iEvent), iPart);
    likeAll[iPart] *=  like[iPart][iChamb];
  }
      }

      // get total probability and normalize it
      for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
  totProb += likeAll[iPart];
      }
      for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
  likeAll[iPart] /= totProb;
      }

      // fill likelihood distributions
      if(fPID[AliPID::kElectron] == 1)      
  hElecs -> Fill(likeAll[AliPID::kElectron]);
      if(fPID[AliPID::kPion] == 1)      
  hPions -> Fill(likeAll[AliPID::kElectron]);
    } // end event loop

    // calculate the pion efficiency and fill the graph
    util -> CalculatePionEffi(hElecs, hPions);
    pionEffiTest[iLoop] = util -> GetPionEfficiency();
    pionEffiErrTest[iLoop] = util -> GetError();

    gEffisTest -> SetPoint(iLoop, iLoop+1, pionEffiTest[iLoop]);
    gEffisTest -> SetPointError(iLoop, 0, pionEffiErrTest[iLoop]);
    hElecs -> Reset();
    hPions -> Reset();
    AliDebug(2, Form("TestLoop[%d] PionEfficiency[%f +/- %f] \n", iLoop, pionEffiTest[iLoop], pionEffiErrTest[iLoop]));
    
  } //   end training loop

  util -> Delete();

  gEffisTest -> Draw("PAL");
  gEffisTrain -> Draw("PL");

}


//________________________________________________________________________
void AliTRDpidRefMakerNN::LoadFile(const Char_t *InFileNN) 
{
  //
  // Loads the files and sets the event list
  // for neural network training and 
  // building of the 2-dim reference histograms.
  // Useable for training outside of the makeResults.C macro
  //

  TFile *fInFileNN;
  fInFileNN = new TFile(InFileNN, "READ");
  fData = (TTree*)fInFileNN -> Get("NN");

  for(Int_t iMom = 0; iMom < AliTRDCalPID::kNMom; iMom++){
    for(Int_t ily = 0; ily < AliTRDgeometry::kNlayer; ily++){
      fTrain[iMom][ily] = new TEventList(Form("fTrainMom%d_%d", iMom, ily), Form("Training list for momentum intervall %d and plane %d", iMom, ily));
      fTest[iMom][ily] = new TEventList(Form("fTestMom%d_%d", iMom, ily), Form("Test list for momentum intervall %d and plane %d", iMom, ily));
    }
  }
}


// //________________________________________________________________________
// void AliTRDpidRefMakerNN::LoadContainer(const Char_t *InFileCont) 
// {

//   //
//   // Loads the container if no container is there.
//   // Useable for training outside of the makeResults.C macro
//   //

//   TFile *fInFileCont;
//   fInFileCont = new TFile(InFileCont, "READ");
//   fContainer = (TObjArray*)fInFileCont -> Get("PidRefMaker");

// }


// //________________________________________________________________________
// void AliTRDpidRefMakerNN::CreateGraphs()
// {
//   // Create histograms
//   // Called once

//   OpenFile(0, "RECREATE");
//   fContainer = new TObjArray();
//   fContainer->AddAt(new TH1F("hPDG","hPDG",AliPID::kSPECIES,-0.5,5.5),0);

//   TGraphErrors *gEffisTrain = new TGraphErrors(kMoniTrain);
//   gEffisTrain -> SetLineColor(4);
//   gEffisTrain -> SetMarkerColor(4);
//   gEffisTrain -> SetMarkerStyle(29);
//   gEffisTrain -> SetMarkerSize(2);

//   TGraphErrors *gEffisTest = new TGraphErrors(kMoniTrain);
//   gEffisTest -> SetLineColor(2);
//   gEffisTest -> SetMarkerColor(2);
//   gEffisTest -> SetMarkerSize(2);

//   fContainer -> AddAt(gEffisTrain,kGraphTrain);
//   fContainer -> AddAt(gEffisTest,kGraphTest);
// }

