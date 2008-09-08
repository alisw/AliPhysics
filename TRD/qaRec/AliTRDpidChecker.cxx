#include "TPDGCode.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include <TClonesArray.h>
#include <TObjArray.h>
#include <TList.h>

#include "AliPID.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliTrackReference.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliTRDtrackV1.h"
#include "AliTRDReconstructor.h"
#include "AliCDBManager.h"
#include "../Cal/AliTRDCalPID.h"


#include "AliTRDpidChecker.h"
#include "AliTRDtrackInfo/AliTRDtrackInfo.h"

// calculate pion efficiency at 90% electron efficiency for 11 momentum bins
// this task should be used with simulated data only

ClassImp(AliTRDpidChecker)

//________________________________________________________________________
AliTRDpidChecker::AliTRDpidChecker(const char *name) 
  :AliAnalysisTask(name, "")
  ,fObjectContainer(0x0)
  ,fTracks(0x0)
  ,fReconstructor(0x0)
//   ,fDebugLevel(1)
//   ,fDebugStream(0x0)
{
  //
  // Default constructor
  //

  fReconstructor = new AliTRDReconstructor();
  fReconstructor->SetRecoParam(AliTRDrecoParam::GetLowFluxParam());

  DefineInput(0, TObjArray::Class());
  DefineOutput(0, TObjArray::Class());
}


//________________________________________________________________________
AliTRDpidChecker::~AliTRDpidChecker() 
{
  if(fObjectContainer){ 
    //fObjectContainer->Delete();
    //delete fObjectContainer;
  }
}



//________________________________________________________________________
void AliTRDpidChecker::ConnectInputData(Option_t *) 
{
  //
  // Connect input data
  //

  fTracks = dynamic_cast<TObjArray*>(GetInputData(0));
}


//________________________________________________________________________
void AliTRDpidChecker::CreateOutputObjects()
{
  // Create histograms
  // Called once

  OpenFile(0, "RECREATE");
  fObjectContainer = new TObjArray();
  fObjectContainer->Add(new TH1F("hMom", "momentum distribution", AliTRDCalPID::kNMom, 0.5, 11.5));


  // histos of the electron probability of all 5 particle species and 11 momenta for the 2-dim LQ method 
  const Int_t kBins = 12000;         // binning of the histos and eficiency calculation
  for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
    for(Int_t iMom = 0; iMom < AliTRDCalPID::kNMom; iMom++){
      fObjectContainer->Add(new TH1F(Form("PID%d_%d_LQ",iPart,iMom), Form("PID distribution for %d_%d", iPart, iMom), kBins, -0.1, 1.1));
    }
  }

  // histos of the electron probability of all 5 particle species and 11 momenta for the neural network method 
  for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
    for(Int_t iMom = 0; iMom < AliTRDCalPID::kNMom; iMom++){
      fObjectContainer->Add(new TH1F(Form("PID%d_%d_NN",iPart,iMom), Form("PID distribution for %d_%d", iPart, iMom), kBins, -0.1, 1.1));
    }
  }

  // frame and TGraph of the pion efficiencies
  fObjectContainer -> Add(new TH2F("hFrame", "", 10, 0.4, 12., 10, 0.0005, 0.1));
  fObjectContainer -> Add(new TGraph(AliTRDCalPID::kNMom));
  fObjectContainer -> Add(new TGraphErrors(AliTRDCalPID::kNMom));
  fObjectContainer -> Add(new TGraph(AliTRDCalPID::kNMom));
  fObjectContainer -> Add(new TGraphErrors(AliTRDCalPID::kNMom));
}

//________________________________________________________________________
void AliTRDpidChecker::Exec(Option_t *) 
{
  // Main loop
  // Called for each event


//   if(!AliTracker::GetFieldMap()){
//     AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., AliMagFMaps::k5kG);
//     AliTracker::SetFieldMap(field, kTRUE);
//   }

  TH1F *hMom = (TH1F*)fObjectContainer->UncheckedAt(0);	
  TH1F *hPIDLQ[AliPID::kSPECIES][AliTRDCalPID::kNMom];
  TH1F *hPIDNN[AliPID::kSPECIES][AliTRDCalPID::kNMom];

  for(Int_t iPart = 0; iPart < AliPID::kSPECIES; iPart++){
    for(Int_t iMom = 0; iMom < AliTRDCalPID::kNMom; iMom++){
      hPIDLQ[iPart][iMom] = (TH1F*)fObjectContainer->At(iPart*AliTRDCalPID::kNMom+iMom+1);
      hPIDNN[iPart][iMom] = (TH1F*)fObjectContainer->At(iPart*AliTRDCalPID::kNMom+iMom+1+AliPID::kSPECIES*AliTRDCalPID::kNMom);
    }
  }
	
  Int_t labelsacc[10000]; 
  memset(labelsacc, 0, sizeof(Int_t) * 10000);
  
  Float_t mom;
  ULong_t status;
  Int_t nTRD = 0;


  AliTRDtrackInfo     *track = 0x0;
  AliTRDtrackV1 *TRDtrack = 0x0;
  AliTrackReference     *ref = 0x0;
  AliExternalTrackParam *esd = 0x0;
  for(Int_t itrk=0; itrk<fTracks->GetEntriesFast(); itrk++){
    track = (AliTRDtrackInfo*)fTracks->UncheckedAt(itrk);
    if(!track->HasESDtrack()) continue;
    status = track->GetStatus();
    if(!(status&AliESDtrack::kTPCout)) continue;

    if(!(TRDtrack = track->GetTRDtrack())) continue; 
    //&&(track->GetNumberOfClustersRefit()
      // use only tracks taht hit 6 chambers
    if(!(TRDtrack->GetNumberOfTracklets() == AliTRDCalPID::kNPlane)) continue;
     
    ref = track->GetTrackRef(0);
    esd = track->GetOuterParam();
    mom = ref ? ref->P(): esd->P();

    labelsacc[nTRD] = track->GetLabel();
    nTRD++;
      
    // fill momentum histo to have the momentum distribution
    hMom -> Fill(mom);


    // set the 11 momentum bins
    Int_t iMomBin = -1;
    if(mom < .7) iMomBin = 0;
    else if(mom < .9) iMomBin = 1;
    else if(mom < 1.25) iMomBin = 2;
    else if(mom < 1.75) iMomBin = 3;
    else if(mom < 2.5) iMomBin = 4;
    else if(mom < 3.5) iMomBin = 5;
    else if(mom < 4.5) iMomBin = 6;
    else if(mom < 5.5) iMomBin = 7;
    else if(mom < 7.0) iMomBin = 8;
    else if(mom < 9.0) iMomBin = 9;
    else  iMomBin = 10;

    // set the reconstructor
    TRDtrack -> SetReconstructor(fReconstructor);

    
    // calculate the probabilities and fill histograms for electrons using 2-dim LQ
    fReconstructor -> SetOption("!nn");
    TRDtrack -> CookPID();

    switch(track->GetPDG()){
    case kElectron:
    case kPositron:
      hPIDLQ[AliPID::kElectron][iMomBin] -> Fill(TRDtrack -> GetPID(0));
      break;
    case kMuonPlus:
    case kMuonMinus:
      hPIDLQ[AliPID::kMuon][iMomBin] -> Fill(TRDtrack -> GetPID(0));
      break;
    case kPiPlus:
    case kPiMinus:
      hPIDLQ[AliPID::kPion][iMomBin] -> Fill(TRDtrack -> GetPID(0));
      break;
    case kKPlus:
    case kKMinus:
      hPIDLQ[AliPID::kKaon][iMomBin] -> Fill(TRDtrack -> GetPID(0));
      break;
    case kProton:
    case kProtonBar:
      hPIDLQ[AliPID::kMuon][iMomBin] -> Fill(TRDtrack -> GetPID(0));
      break;
    }


    // calculate the probabilities and fill histograms for electrons using NN
    fReconstructor -> SetOption("nn");
    TRDtrack->CookPID();
    switch(track->GetPDG()){
    case kElectron:
    case kPositron:
      hPIDNN[AliPID::kElectron][iMomBin] -> Fill(TRDtrack -> GetPID(0));
      break;
    case kMuonPlus:
    case kMuonMinus:
      hPIDNN[AliPID::kMuon][iMomBin] -> Fill(TRDtrack -> GetPID(0));
      break;
    case kPiPlus:
    case kPiMinus:
      hPIDNN[AliPID::kPion][iMomBin] -> Fill(TRDtrack -> GetPID(0));
      break;
    case kKPlus:
    case kKMinus:
      hPIDNN[AliPID::kKaon][iMomBin] -> Fill(TRDtrack -> GetPID(0));
      break;
    case kProton:
    case kProtonBar:
      hPIDNN[AliPID::kMuon][iMomBin] -> Fill(TRDtrack -> GetPID(0));
      break;
    }
  }

  PostData(0, fObjectContainer);
}

//________________________________________________________________________
void AliTRDpidChecker::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  fObjectContainer = dynamic_cast<TObjArray*>(GetOutputData(0));
  if (!fObjectContainer) {
    Printf("ERROR: list not available");
    return;
  }

  
  // container for the pion efficiencies and the errors
  Double_t  PionEffiLQ[AliTRDCalPID::kNMom], 
            PionEffiNN[AliTRDCalPID::kNMom],
            PionEffiErrorLQ[AliTRDCalPID::kNMom], 
            PionEffiErrorNN[AliTRDCalPID::kNMom];


  // calculate the pion efficiencies and the errors for 90% electron efficiency (2-dim LQ)
  for(Int_t iMom = 0; iMom < AliTRDCalPID::kNMom; iMom++){
    PionEffiLQ[iMom] = GetPionEfficiency(iMom+1,iMom+23);
    PionEffiErrorLQ[iMom] = GetError(iMom+1,iMom+23);
    Printf("Pion Efficiency for 2-dim LQ is : %f +/- %f\n\n", PionEffiLQ[iMom], PionEffiErrorLQ[iMom]);
  }
  
  // calculate the pion efficiencies and the errors for 90% electron efficiency (NN)
  for(Int_t iMom = 0; iMom < AliTRDCalPID::kNMom; iMom++){
    PionEffiNN[iMom] = GetPionEfficiency(iMom+56,iMom+78);
    PionEffiErrorNN[iMom] = GetError(iMom+56,iMom+78);
    Printf("Pion Efficiency for NN is : %f +/- %f\n\n", PionEffiNN[iMom], PionEffiErrorNN[iMom]);
  }
  

  // create TGraph to plot the pion efficiencies
  TGraph *gEffisLQ=0x0, *gEffisNN=0x0;
  TGraphErrors *gEffisErrLQ=0x0, *gEffisErrNN=0x0;

  gEffisLQ = (TGraph*)fObjectContainer->At(112);
  gEffisErrLQ = (TGraphErrors*)fObjectContainer->At(113);
  gEffisNN = (TGraph*)fObjectContainer->At(114);
  gEffisErrNN = (TGraphErrors*)fObjectContainer->At(115);

  for(Int_t iBin = 0; iBin < AliTRDCalPID::kNMom; iBin++){
    Float_t momentum = AliTRDCalPID::GetMomentum(iBin);
    gEffisLQ->SetPoint(iBin, momentum, PionEffiLQ[iBin]);
    gEffisErrLQ->SetPoint(iBin, momentum, PionEffiLQ[iBin]);
    gEffisErrLQ->SetPointError(iBin, 0., PionEffiErrorLQ[iBin]);

    gEffisNN->SetPoint(iBin, momentum, PionEffiNN[iBin]);
    gEffisErrNN->SetPoint(iBin, momentum, PionEffiNN[iBin]);
    gEffisErrNN->SetPointError(iBin, 0., PionEffiErrorNN[iBin]);
  }

  gEffisLQ -> SetNameTitle("gEffisLQ", "Efficiencies of the 2-dim LQ method");
  gEffisErrLQ -> SetNameTitle("gEffisLQErr", "Efficiencies and Errors of the 2-dim LQ method");
  gEffisNN -> SetNameTitle("gEffisNN", "Efficiencies of the NN method");
  gEffisErrNN -> SetNameTitle("gEffisNNErr", "Efficiencies and Errors of the NN method");
}


//________________________________________________________________________
Double_t  AliTRDpidChecker::GetPionEfficiency(Int_t Index1, Int_t Index2){

  Float_t Multi = 0.9;           // electron efficiency
  Int_t abin, bbin;              
  Double_t SumElecs, SumPions;   // integrated sum of elecs and pions in the histos
  Double_t aBinSum, bBinSum;     // content of a single bin in the histos
  
  TH1F *Histo1 = (TH1F*)fObjectContainer->At(Index1);  // electron histo
  TH1F *Histo2 = (TH1F*)fObjectContainer->At(Index2);  // pion histo


  SumElecs = 0.;
  if(!(Histo1 -> GetEntries() > 0 && Histo2 -> GetEntries() > 0)){
    Printf("Warning: Histo momentum intervall %d has no Entries!", Index1-1);
    return -1.;
  }
  Histo1 -> Scale(1./Histo1->GetEntries());
  Histo2 -> Scale(1./Histo2->GetEntries());


  // calculate threshold for pion efficiency
  for (abin = (Histo1 -> GetNbinsX()); abin >= 0; abin--){  
    aBinSum = 0;
    aBinSum = Histo1 -> GetBinContent(abin);
    if(!(aBinSum == 0)){
      SumElecs = SumElecs + aBinSum;
    }
    if (SumElecs >= Multi){
      break;
    }
  }

    
  // calculate pion efficiency
  SumPions = 0.;
  for (bbin = (Histo2 -> GetNbinsX()); bbin >= abin; bbin--){	
    bBinSum = 0;
    bBinSum = Histo2 -> GetBinContent(bbin);
    if(!(bBinSum == 0)){
      SumPions = SumPions + bBinSum;
    }
  }
  
  
  // print the electron efficiency and its cuts
  Printf("Cut for momentum intervall %d and electron efficiency of %f for: 0.%d", Index1-1, SumElecs, abin-1000);
  Printf("(%.0f electrons and %.0f pions)",Histo1 -> GetEntries(), Histo2 -> GetEntries());


  // return the pion efficiency
  return SumPions;

} 
  

//________________________________________________________________________
Double_t  AliTRDpidChecker::GetError(Int_t Index1, Int_t Index2){

  
  const Int_t kBins = 12000;         // binning of the histos and eficiency calculation
  const Float_t Multi = 0.9;                          // electron efficiency
  Int_t abinE, bbinE, cbinE;                    
  Double_t SumElecsE[kBins], SumPionsE[kBins];  // array of the integrated sum in each bin
  Double_t aBinSumE, bBinSumE;                  // content of a single bin
  Double_t EleEffi, PioEffi;                    // electron and pion efficiency
  Bool_t bFirst = 1;                            // checks if threshold is crossed for the first time
  Double_t fError = 0.;                         // error of the pion efficiency


  TH1F *Histo1 = (TH1F*)fObjectContainer->At(Index1);  // electron histo
  TH1F *Histo2 = (TH1F*)fObjectContainer->At(Index2);  // pion histo

  if(!(Histo1 -> GetEntries() > 0 && Histo2 -> GetEntries() > 0)){
    Printf("Warning: Histo momentum intervall %d has no Entries!", Index1-1);
    return -1.;
  }

  for(Int_t iBin = 0; iBin < kBins; iBin++){
    SumElecsE[iBin] = 0.;
    SumPionsE[iBin] = 0.;
  }

  EleEffi = 0.;
  PioEffi = 0.;
  cbinE = -1;


  // calculate electron efficiency of each bin
  for (abinE = (Histo1 -> GetNbinsX())-2; abinE >= 0; abinE--){  
    aBinSumE = 0;
    aBinSumE = Histo1 -> GetBinContent(abinE);
    
    SumElecsE[abinE] = SumElecsE[abinE+1] + aBinSumE;
    if((SumElecsE[abinE] >= Multi) && (bFirst == 1)){
      bFirst = 0;
      cbinE = abinE;
      EleEffi = (SumElecsE[cbinE]); 
    }
  }
  

  // calculate pion efficiency of each bin
  for (bbinE = (Histo2 -> GetNbinsX())-2; bbinE >= abinE; bbinE--){	
    bBinSumE = 0;
    bBinSumE = Histo2 -> GetBinContent(bbinE);

    SumPionsE[bbinE] = SumPionsE[bbinE+1] + bBinSumE;
    if(bbinE == cbinE){
      PioEffi = (SumPionsE[cbinE]);
    }
  }
  

  // pion efficiency vs electron efficiency
  TGraph *gEffis = new TGraph(kBins, SumElecsE, SumPionsE);

  // use fit function to get derivate of the TGraph for error calculation
  TF1 *f1 = new TF1("f1","[0]*x*x+[1]*x+[2]", Multi-.05, Multi+.05);
  gEffis -> Fit("f1","Q","",Multi-.05, Multi+.05);
  Printf("Derivative at %4.2f : %f", Multi, f1 -> Derivative(Multi));

  // return the error of the pion efficiency
  fError = sqrt(((1/Histo2 -> GetEntries())*PioEffi*(1-PioEffi))+((f1 -> Derivative(Multi))*(f1 -> Derivative(Multi))*(1/Histo1 -> GetEntries())*EleEffi*(1-EleEffi)));
  return fError;
}
