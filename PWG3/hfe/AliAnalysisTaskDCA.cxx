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
//
// The analysis task:
// impact parameter resolution and pull study
// for tracks which survivied the particle cuts
// 
// 
// Authors:
//  Hongyan Yang <hongyan@physi.uni-heidelberg.de>
//  Carlo Bombonati <carlo.bombonati@cern.ch>
//
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TList.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TParticle.h>
#include <TString.h>

#include <TCanvas.h>

#include "AliCFManager.h"
#include "AliMCEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"

#include "AliHFEcuts.h"
#include "AliHFEdca.h"

#include "AliAnalysisTaskDCA.h"


//____________________________________________________________
AliAnalysisTaskDCA::AliAnalysisTaskDCA():
  AliAnalysisTask("Impact Parameter Resolution and Pull Analysis", "")
  , fPlugins(0)
  , fESD(0x0)
  , fMC(0x0)
  , fCuts(0x0)
  , fCFM(0x0)
  , fDCA(0x0)
  , fPixelStatus(0x0)
  , fNclustersITS(0x0)
  , fNEvents(0x0)
  , fResidualList(0x0)
  , fPullList(0x0)
  , fOutput(0x0)
{
  //
  // Dummy constructor
  //
  DefineInput(0, TChain::Class());
  DefineOutput(0, TH1I::Class());   // event
  DefineOutput(1, TList::Class());  // output

  SetHasMCData();

}

//____________________________________________________________
AliAnalysisTaskDCA::AliAnalysisTaskDCA(const char * name):
  AliAnalysisTask(name, "")
  , fPlugins(0)
  , fESD(0x0)
  , fMC(0x0)
  , fCuts(0x0)
  , fCFM(0x0)
  , fDCA(0x0)
  , fPixelStatus(0x0)
  , fNclustersITS(0x0)
  , fNEvents(0x0)
  , fResidualList(0x0)
  , fPullList(0x0)
  , fOutput(0x0)
{
  //
  // Default constructor
  //
  DefineInput(0, TChain::Class());
  DefineOutput(0, TH1I::Class());
  DefineOutput(1, TList::Class());

  SetHasMCData();
}

//____________________________________________________________
AliAnalysisTaskDCA::AliAnalysisTaskDCA(const AliAnalysisTaskDCA &ref):
  AliAnalysisTask(ref)
  , fPlugins(ref.fPlugins)
  , fESD(ref.fESD)
  , fMC(ref.fMC)
  , fCuts(ref.fCuts)
  , fCFM(ref.fCFM)
  , fDCA(ref.fDCA)
  , fPixelStatus(ref.fPixelStatus)
  , fNclustersITS(ref.fNclustersITS)
  , fNEvents(ref.fNEvents)
  , fResidualList(ref.fResidualList)
  , fPullList(ref.fPullList)
  , fOutput(ref.fOutput)
{
  //
  // Copy Constructor
  //
}

//____________________________________________________________
AliAnalysisTaskDCA &AliAnalysisTaskDCA::operator=(const AliAnalysisTaskDCA &ref){
  //
  // Assignment operator
  //
  if(this == &ref) return *this;
  AliAnalysisTask::operator=(ref);
  fPlugins = ref.fPlugins;
  fESD = ref.fESD;
  fMC = ref.fMC;
  fCuts = ref.fCuts;
  fCFM = ref.fCFM;
  fDCA = ref.fDCA;
  fPixelStatus = ref.fPixelStatus;
  fNclustersITS = ref.fNclustersITS;
  fNEvents = ref.fNEvents;
  fOutput = ref.fOutput;
  fResidualList = ref.fResidualList;
  fPullList = ref.fPullList;

  return *this;
}

//____________________________________________________________
AliAnalysisTaskDCA::~AliAnalysisTaskDCA(){
  //
  // Destructor
  //

  if(fESD) delete fESD;
  if(fMC) delete fMC;

  if(fCFM) delete fCFM;

  if(fDCA) delete fDCA;  

  if(fOutput){ 
    fOutput->Clear();
    delete fOutput;
  }
  
  if(fResidualList){ 
    fResidualList->Clear();
    delete fResidualList;
  }

  if(fPullList){ 
    fPullList->Clear();
    delete fPullList;
  }
  
  if(fNEvents) delete fNEvents;
}

//____________________________________________________________
void AliAnalysisTaskDCA::ConnectInputData(Option_t *){
  //
  // Connecting the input
  
  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!esdH){      
    AliError("No ESD input handler");
    return;
  } else {
    fESD = esdH->GetEvent();
  }
  AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if(!mcH){       
    AliError("No MC truth handler");
    return;
  } else {
    fMC = mcH->MCEvent();
  }
}

//____________________________________________________________
void AliAnalysisTaskDCA::CreateOutputObjects(){
  // create output objects
  // fNEvents
  // residual and pull

  fNEvents = new TH1I("nEvents", "Number of Events in the Analysis", 2, 0, 2); // Number of Events neccessary for the analysis and not a QA histogram
  if(!fOutput) fOutput = new TList;
  // Initialize correction Framework and Cuts
  fCFM = new AliCFManager;
  MakeParticleContainer();
  // Temporary fix: Initialize particle cuts with 0x0
  for(Int_t istep = 0; istep < fCFM->GetParticleContainer()->GetNStep(); istep++)
    fCFM->SetParticleCutsList(istep, 0x0);
  if(!fCuts){
    AliWarning("Cuts not available. Default cuts will be used");
    fCuts = new AliHFEcuts;
    fCuts->CreateStandardCuts();
  }
  
  fCuts->SetCutITSpixel(fPixelStatus);
  fCuts->Initialize(fCFM);
  

  // dca study----------------------------------
  if(!fDCA) fDCA = new AliHFEdca;
  if(!fResidualList) fResidualList = new TList();
  if(!fPullList) fPullList = new TList();

  fDCA->CreateHistogramsResidual(fResidualList);
  fDCA->CreateHistogramsPull(fPullList);
  
  // add output objects to the List
  fOutput->AddAt(fResidualList,0);
  fOutput->AddAt(fPullList,1);

}

//____________________________________________________________
void AliAnalysisTaskDCA::Exec(Option_t *){
  //
  // Run the analysis
  // 

  AliDebug(3, "Processing ESD events");

  if(!fESD){
    AliError("No ESD Event");
    return;
  }
  if(!fMC){
    AliError("No MC Event");
    return;
  }
  if(!fCuts){
    AliError("HFE cuts not available");
    return;
  }

  //
  // Loop ESD
  //
  AliESDtrack *track = 0x0;  
  fCFM->SetRecEventInfo(fESD);
  // event cut level
  if(!fCFM->CheckEventCuts(AliHFEcuts::kEventStepReconstructed, fESD)) return;

  for(Int_t itrack = 0; itrack < fESD->GetNumberOfTracks(); itrack++){

    track = fESD->GetTrack(itrack);

    // RecPrim: primary cuts
    if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepRecPrim, track)) continue;
    // RecKine: ITSTPC cuts  
    if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepRecKineITSTPC, track)) continue;
    // HFEcuts: ITS layers cuts
    if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepHFEcutsITS, track)) continue;
    
    if(track->GetITSclusters(0)<=fNclustersITS) continue;  // require 6 hits on all pixel layers

    if(GetPlugin(kImpactPar)) {
      //      Printf("analysis on impact parameter is ON");
      fDCA->FillHistograms(fESD, track, fMC);
    }   
    
  }  
  fNEvents->Fill(1);
  
  PostData(0, fNEvents);
  PostData(1, fOutput);
}

//____________________________________________________________
void AliAnalysisTaskDCA::Terminate(Option_t *){
  //
  // Terminate not implemented at the moment
  //  
  
  if(GetPlugin(kPostProcess)){
    fOutput = dynamic_cast<TList *>(GetOutputData(1));
    if(!fOutput){
      AliError("Results not available");
      return;
    }
    PostProcess();
  }
  
}


//____________________________________________________________
void AliAnalysisTaskDCA::Load(TString filename){
  // no need for postprocessing for the moment
  TFile *input = TFile::Open(filename.Data());
  if(!input || input->IsZombie()){
    AliError("Cannot read file");
    return;
  }

  /* 
  TH1 *htmp = dynamic_cast<TH1I *>(input->Get("nEvents"));
  if(htmp)
    fNEvents = dynamic_cast<TH1I *>(htmp->Clone());
  else
    AliError("Event Counter histogram not found"); 
  */
  input->Close();
  delete input;
  
  
}

//____________________________________________________________
void AliAnalysisTaskDCA::PostProcess(){
  // do post processing
  // should do fitting here for dca resolution
  // moved to an external macro to do the job
  
  Load("impactPar.root");
  TCanvas *c1 = new TCanvas("c1", "number of analyzed events", 300, 400);
  fNEvents->Draw();
  c1->SaveAs("temp.png");
 
}




//____________________________________________________________
void AliAnalysisTaskDCA::PrintStatus() const {
  
  //
  // Print Analysis status
  //
  printf("\n\tAnalysis Settings\n\t========================================\n");
  printf("\timpact parameter analysis is %s\n", GetPlugin(kImpactPar)?"ON":"OFF");
  printf("\tcuts: %s\n", (fCuts != NULL) ? "YES" : "NO");
  printf("\t ");
  printf("\n");
}

//__________________________________________                                                  
void AliAnalysisTaskDCA::SwitchOnPlugin(Int_t plug){
  //                                            
  // Switch on Plugin          
  // Available:                                  
  //  - analyze impact parameter
  //  - Post Processing                                                                      
  
  switch(plug){
  case kImpactPar: 
    SETBIT(fPlugins, plug); 
    break;
  case kPostProcess: 
    SETBIT(fPlugins, plug); 
    break;
  default: 
    AliError("Unknown Plugin");
  };
}


//____________________________________________________________
void AliAnalysisTaskDCA::MakeParticleContainer(){
  //
  // Create the particle container (borrowed from AliAnalysisTaskHFE)
  //
  const Int_t kNvar   = 3 ; //number of variables on the grid:pt,eta, phi
  const Double_t kPtmin = 0.1, kPtmax = 10.;
  const Double_t kEtamin = -0.9, kEtamax = 0.9;
  const Double_t kPhimin = 0., kPhimax = 2. * TMath::Pi();

  //arrays for the number of bins in each dimension
  Int_t iBin[kNvar];
  iBin[0] = 40; //bins in pt
  iBin[1] =  8; //bins in eta 
  iBin[2] = 18; // bins in phi

  //arrays for lower bounds :
  Double_t* binEdges[kNvar];
  for(Int_t ivar = 0; ivar < kNvar; ivar++)
    binEdges[ivar] = new Double_t[iBin[ivar] + 1];

  //values for bin lower bounds
  for(Int_t i=0; i<=iBin[0]; i++) binEdges[0][i]=(Double_t)TMath::Power(10,TMath::Log10(kPtmin) + (TMath::Log10(kPtmax)-TMath::Log10(kPtmin))/iBin[0]*(Double_t)i);  
  for(Int_t i=0; i<=iBin[1]; i++) binEdges[1][i]=(Double_t)kEtamin  + (kEtamax-kEtamin)/iBin[1]*(Double_t)i;
  for(Int_t i=0; i<=iBin[2]; i++) binEdges[2][i]=(Double_t)kPhimin  + (kPhimax-kPhimin)/iBin[2]*(Double_t)i;

  //one "container" for MC
  AliCFContainer* container = new AliCFContainer("container","container for tracks", (AliHFEcuts::kNcutStepsTrack + 1 + 2*(AliHFEcuts::kNcutStepsESDtrack + 1)), kNvar, iBin);

  //setting the bin limits
  for(Int_t ivar = 0; ivar < kNvar; ivar++)
    container -> SetBinLimits(ivar, binEdges[ivar]);
  fCFM->SetParticleContainer(container);

  //create correlation matrix for unfolding
  Int_t thnDim[2*kNvar];
  for (int k=0; k<kNvar; k++) {
    //first half  : reconstructed 
    //second half : MC
    thnDim[k]      = iBin[k];
    thnDim[k+kNvar] = iBin[k];
  }


  // add more containers for correction purpose


}
