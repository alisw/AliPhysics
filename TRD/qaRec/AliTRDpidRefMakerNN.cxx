#include "TSystem.h"
#include "TDatime.h"
#include "TEventList.h"
#include "TGraphErrors.h"
#include "TTreeStream.h"
#include "TH1F.h"
#include "TMultiLayerPerceptron.h"

#include "AliLog.h"
#include "AliPID.h"
#include "AliESDEvent.h"

#include "../Cal/AliTRDCalPIDNN.h"

#include "AliTRDpidUtil.h"

#include "AliTRDpidRefMakerNN.h"
#include "info/AliTRDtrackInfo.h"

// builds the reference tree for the training of neural networks


ClassImp(AliTRDpidRefMakerNN)

//________________________________________________________________________
AliTRDpidRefMakerNN::AliTRDpidRefMakerNN() 
  :AliTRDpidRefMaker(UChar_t(AliTRDpidUtil::kNNslices), "PidRefMakerNN", "PID(NN) Reference Maker")
  ,fTrainMomBin(kAll)
  ,fEpochs(1000)
  ,fMinTrain(100)
  ,fDate(0)
  ,fDoTraining(0)
  ,fContinueTraining(0)
  ,fTrainPath(0x0)
{
  //
  // Default constructor
  //
  SetAbundance(.67, .33);
  SetScaledEdx(Float_t(AliTRDCalPIDNN::kMLPscale));
  TDatime datime;
  fDate = datime.GetDate();
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
  MakeTrainingList();        
  if(!fDoTraining) return kTRUE;

  // train the neural networks and build the refrence histos for 2-dim LQ
  gSystem->Exec(Form("mkdir ./Networks_%d/",fDate));
  AliDebug(3, Form("TrainMomBin [%d] [%d]", fTrainMomBin, kAll));

  // train single network for a single momentum (recommended)
  if(!(fTrainMomBin == kAll)){
    if(fTrain[fTrainMomBin][0] -> GetN() < fMinTrain){
      AliWarning("Not enough events for training available! Please check Data sample!");
      return kFALSE;
    }
    MakeRefs(fTrainMomBin);
    MonitorTraining(fTrainMomBin);
  }
  // train all momenta
  else{
    for(Int_t iMomBin = 0; iMomBin < AliTRDCalPID::kNMom; iMomBin++){
      if(fTrain[iMomBin][0] -> GetN() < fMinTrain){
        AliWarning(Form("Not enough events for training available for momentum bin [%d]! Please check Data sample!", iMomBin));
        continue;
      }
      MakeRefs(fTrainMomBin);
      MonitorTraining(iMomBin);
    }
  }

  return kTRUE; // testing protection
}



//________________________________________________________________________
void AliTRDpidRefMakerNN::MakeRefs(Int_t mombin) 
{
  //
  // train the neural networks
  //
  
  
  if (!fData) LoadFile(Form("TRD.%s.root", GetName()));

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
	  if(fDebugLevel>=2) fNet[iChamb] -> Train(nEpochs,"text update=10, graph");
	  else fNet[iChamb] -> Train(nEpochs,"");
	}
	else{
	  fNet[iChamb] -> LoadWeights(Form("./Networks_%d/MomBin_%d/Net%d_%d",fTrainPath, mombin, iChamb, kMoniTrain - 1));
	  if(fDebugLevel>=2) fNet[iChamb] -> Train(nEpochs,"text update=10, graph+");      
	  else fNet[iChamb] -> Train(nEpochs,"+");                   
	}
	bFirstLoop[iChamb] = kFALSE;
      }
      else{    
	if(fDebugLevel>=2) fNet[iChamb] -> Train(nEpochs,"text update=10, graph+");      
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
  
  if(!fContainer) LoadContainer(Form("TRD.Task%s.root", GetName()));

  if(!fContainer){
    AliError("ERROR container not available");
    return;
  }

  if (!fData) LoadFile(Form("TRD.Task%s.root", GetName()));
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
  Float_t Like[AliPID::kSPECIES][AliTRDgeometry::kNlayer];
  memset(Like, 0, AliPID::kSPECIES*AliTRDgeometry::kNlayer*sizeof(Float_t));
  Float_t LikeAll[AliPID::kSPECIES], TotProb;

  Double_t PionEffiTrain[kMoniTrain], PionEffiErrTrain[kMoniTrain];
  Double_t PionEffiTest[kMoniTrain], PionEffiErrTest[kMoniTrain];
  memset(PionEffiTrain, 0, kMoniTrain*sizeof(Double_t));
  memset(PionEffiErrTrain, 0, kMoniTrain*sizeof(Double_t));
  memset(PionEffiTest, 0, kMoniTrain*sizeof(Double_t));
  memset(PionEffiErrTest, 0, kMoniTrain*sizeof(Double_t));

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
	LikeAll[iPart] = 1./AliPID::kSPECIES;
      }
      TotProb = 0.;

      fData -> GetEntry(fTrain[mombin][0] -> GetEntry(iEvent));
      // use event only if it is electron or pion
      if(!((fPID[AliPID::kElectron] == 1.0) || (fPID[AliPID::kPion] == 1.0))) continue;

      // get the probabilities for each particle type in each chamber
      for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
	for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
	  Like[iPart][iChamb] = fNet[iChamb] -> Result(fTrain[mombin][iChamb] -> GetEntry(iEvent), iPart);
	  LikeAll[iPart] *=  Like[iPart][iChamb];
	}
      }

      // get total probability and normalize it
      for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
	TotProb += LikeAll[iPart];
      }
      for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
	LikeAll[iPart] /= TotProb;
      }

      // fill likelihood distributions
      if(fPID[AliPID::kElectron] == 1)      
	hElecs -> Fill(LikeAll[AliPID::kElectron]);
      if(fPID[AliPID::kPion] == 1)      
	hPions -> Fill(LikeAll[AliPID::kElectron]);
    } // end event loop


    // calculate the pion efficiency and fill the graph
    util -> CalculatePionEffi(hElecs, hPions);
    PionEffiTrain[iLoop] = util -> GetPionEfficiency();
    PionEffiErrTrain[iLoop] = util -> GetError();

    gEffisTrain -> SetPoint(iLoop, iLoop+1, PionEffiTrain[iLoop]);
    gEffisTrain -> SetPointError(iLoop, 0, PionEffiErrTrain[iLoop]);
    hElecs -> Reset();
    hPions -> Reset();
    if(fDebugLevel>=2) Printf("TrainingLoop[%d] PionEfficiency[%f +/- %f]", iLoop, PionEffiTrain[iLoop], PionEffiErrTrain[iLoop]);
    // end training loop
    


    // event loop test list
    for(Int_t iEvent = 0; iEvent < fTest[mombin][0] -> GetN(); iEvent++ ){

      // reset particle probabilities
      for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
	LikeAll[iPart] = 1./AliTRDgeometry::kNlayer;
      }
      TotProb = 0.;

      fData -> GetEntry(fTest[mombin][0] -> GetEntry(iEvent));
      // use event only if it is electron or pion
      if(!((fPID[AliPID::kElectron] == 1.0) || (fPID[AliPID::kPion] == 1.0))) continue;

      // get the probabilities for each particle type in each chamber
      for(Int_t iChamb = 0; iChamb < AliTRDgeometry::kNlayer; iChamb++){
	for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
	  Like[iPart][iChamb] = fNet[iChamb] -> Result(fTest[mombin][iChamb] -> GetEntry(iEvent), iPart);
	  LikeAll[iPart] *=  Like[iPart][iChamb];
	}
      }

      // get total probability and normalize it
      for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
	TotProb += LikeAll[iPart];
      }
      for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
	LikeAll[iPart] /= TotProb;
      }

      // fill likelihood distributions
      if(fPID[AliPID::kElectron] == 1)      
	hElecs -> Fill(LikeAll[AliPID::kElectron]);
      if(fPID[AliPID::kPion] == 1)      
	hPions -> Fill(LikeAll[AliPID::kElectron]);
    } // end event loop

    // calculate the pion efficiency and fill the graph
    util -> CalculatePionEffi(hElecs, hPions);
    PionEffiTest[iLoop] = util -> GetPionEfficiency();
    PionEffiErrTest[iLoop] = util -> GetError();

    gEffisTest -> SetPoint(iLoop, iLoop+1, PionEffiTest[iLoop]);
    gEffisTest -> SetPointError(iLoop, 0, PionEffiErrTest[iLoop]);
    hElecs -> Reset();
    hPions -> Reset();
    if(fDebugLevel>=2) Printf("TestLoop[%d] PionEfficiency[%f +/- %f] \n", iLoop, PionEffiTest[iLoop], PionEffiErrTest[iLoop]);
    
  } //   end training loop
 
  util -> Delete();

  gEffisTest -> Draw("PAL");
  gEffisTrain -> Draw("PL");

}


