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
/*
 * The analysis task:
 * Filling an AliCFContainer with the quantities pt, eta and phi
 * for tracks which survivied the particle cuts (MC resp. ESD tracks)
 * Track selection is done using the AliHFE package
 * 
 * Author:
 *  Raphaelle Bailhache <R.Bailhache@gsi.de>
 *  Markus Fasel <M.Fasel@gsi.de>
 *  MinJung Kweon <minjung@physi.uni-heidelberg.de>
 */
#include <TAxis.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TH1F.h>
#include <TH1I.h>
#include <TH2F.h>
#include <TIterator.h>
#include <TList.h>
#include <TLegend.h>
#include <TMath.h>
#include <TObjArray.h>
#include <TParticle.h>
#include <TProfile.h>
#include <TString.h>
#include <TTree.h>

#include "AliCFContainer.h"
#include "AliCFManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"
#include "AliPID.h"
#include "AliStack.h"

#include "AliHFEpid.h"
#include "AliHFEcuts.h"
#include "AliHFEmcQA.h"
#include "AliHFEsecVtx.h"
#include "AliAnalysisTaskHFE.h"

//____________________________________________________________
AliAnalysisTaskHFE::AliAnalysisTaskHFE():
  AliAnalysisTask("PID efficiency Analysis", "")
  , fQAlevel(0)
  , fPIDdetectors("")
  , fESD(0x0)
  , fMC(0x0)
  , fCFM(0x0)
  , fCorrelation(0x0)
  , fFakeElectrons(0x0)
  , fPID(0x0)
  , fCuts(0x0)
  , fSecVtx(0x0)
  , fMCQA(0x0)
  , fNEvents(0x0)
  , fNElectronTracksEvent(0x0)
  , fQA(0x0)
  , fOutput(0x0)
  , fHistMCQA(0x0)
  , fHistSECVTX(0x0)
{
  //
  // Default constructor
  //
  DefineInput(0, TChain::Class());
  DefineOutput(0, TH1I::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());

  // Initialize cuts
  fCuts = new AliHFEcuts;
  fPID = new AliHFEpid;
}

//____________________________________________________________
AliAnalysisTaskHFE::AliAnalysisTaskHFE(const AliAnalysisTaskHFE &ref):
  AliAnalysisTask(ref)
  , fQAlevel(ref.fQAlevel)
  , fPIDdetectors(ref.fPIDdetectors)
  , fESD(ref.fESD)
  , fMC(ref.fMC)
  , fCFM(ref.fCFM)
  , fCorrelation(ref.fCorrelation)
  , fFakeElectrons(ref.fFakeElectrons)
  , fPID(ref.fPID)
  , fCuts(ref.fCuts)
  , fSecVtx(ref.fSecVtx)
  , fMCQA(ref.fMCQA)
  , fNEvents(ref.fNEvents)
  , fNElectronTracksEvent(ref.fNElectronTracksEvent)
  , fQA(ref.fQA)
  , fOutput(ref.fOutput)
  , fHistMCQA(ref.fHistMCQA)
  , fHistSECVTX(ref.fHistSECVTX)
{
  //
  // Copy Constructor
  //
}

//____________________________________________________________
AliAnalysisTaskHFE &AliAnalysisTaskHFE::operator=(const AliAnalysisTaskHFE &ref){
  //
  // Assignment operator
  //
  if(this == &ref) return *this;
  AliAnalysisTask::operator=(ref);
  fQAlevel = ref.fQAlevel;
  fPIDdetectors = ref.fPIDdetectors;
  fESD = ref.fESD;
  fMC = ref.fMC;
  fCFM = ref.fCFM;
  fCorrelation = ref.fCorrelation;
  fFakeElectrons = ref.fFakeElectrons;
  fPID = ref.fPID;
  fCuts = ref.fCuts;
  fSecVtx = ref.fSecVtx;
  fMCQA = ref.fMCQA;
  fNEvents = ref.fNEvents;
  fNElectronTracksEvent = ref.fNElectronTracksEvent;
  fQA = ref.fQA;
  fOutput = ref.fOutput;
  fHistMCQA = ref.fHistMCQA;
  fHistSECVTX = ref.fHistSECVTX;
  return *this;
}

//____________________________________________________________
AliAnalysisTaskHFE::~AliAnalysisTaskHFE(){
  //
  // Destructor
  //
  if(fESD) delete fESD;
  if(fMC) delete fMC;
  if(fPID) delete fPID;
  if(fQA){
    fQA->Clear();
    delete fQA;
  }
  if(fOutput){ 
    fOutput->Clear();
    delete fOutput;
  }
  if(fHistMCQA){
    fHistMCQA->Clear();
    delete fHistMCQA;
  }
  if(fHistSECVTX){
    fHistSECVTX->Clear();
    delete fHistSECVTX;
  }
  if(fCuts) delete fCuts;
  if(fSecVtx) delete fSecVtx;
  if(fMCQA) delete fMCQA;
  if(fNEvents) delete fNEvents;
  if(fCorrelation) delete fCorrelation;
  if(fFakeElectrons) delete fFakeElectrons;
}

//____________________________________________________________
void AliAnalysisTaskHFE::ConnectInputData(Option_t *){
/*  TTree *esdchain = dynamic_cast<TChain *>(GetInputData(0));
  if(!esdchain){
    AliError("ESD chain empty");
    return;
  } else {
    esdchain->SetBranchStatus("Tracks", 1);
  }
*/
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
void AliAnalysisTaskHFE::CreateOutputObjects(){
  fNEvents = new TH1I("nEvents", "Number of Events in the Analysis", 2, 0, 2); // Number of Events neccessary for the analysis and not a QA histogram
  fNElectronTracksEvent = new TH1I("nElectronTracksEvent", "Number of Electron Candidates", 100, 0, 100);
  // First Step: TRD alone
  if(!fQA) fQA = new TList;
  fQA->AddAt(new TProfile("conr", "Electron PID contamination", 20, 0, 20), 0);
  fQA->AddAt(new TH1F("alpha_rec", "Alpha from reconstructed tracks with TRD hits", 36, -TMath::Pi(), TMath::Pi()), 1);
  fQA->AddAt(new TH1F("alpha_sim", "Alpha from simulated electron tracks", 36, -TMath::Pi(), TMath::Pi()), 2);
  fQA->AddAt(new TH1F("nElectron", "Number of electrons", 100, 0, 100), 3);
  fQA->AddAt(new TProfile("pidquality", "TRD PID quality as function of momentum", 20, 0, 20), 4);
  fQA->AddAt(new TProfile("ntrdclusters", "Number of TRD clusters as function of momentum", 20, 0, 20), 5);
  fQA->AddAt(new TH1F("chi2TRD","#chi2 per TRD cluster", 20, 0, 20), 6);

  if(!fOutput) fOutput = new TList;
  // Initialize correction Framework and Cuts
  fCFM = new AliCFManager;
  MakeParticleContainer();
  // Temporary fix: Initialize particle cuts with 0x0
  for(Int_t istep = 0; istep < fCFM->GetParticleContainer()->GetNStep(); istep++)
    fCFM->SetParticleCutsList(istep, 0x0);
  if(IsQAOn(kCUTqa)){
    AliInfo("QA on for Cuts");
    fCuts->SetDebugMode();
    fQA->Add(fCuts->GetQAhistograms());
  }
  fCuts->CreateStandardCuts();
  fCuts->SetMinNTrackletsTRD(0);  // Minimal requirement to get a minimum biased electron sample (only TPC and ITS requirements allowed)
  fCuts->Initialize(fCFM);
  // add output objects to the List
  fOutput->AddAt(fCFM->GetParticleContainer(), 0);
  fOutput->AddAt(fCorrelation, 1);
  fOutput->AddAt(fFakeElectrons, 2);
  fOutput->AddAt(fNElectronTracksEvent, 3);

  // Initialize PID
  if(IsQAOn(kPIDqa)){
    AliInfo("PID QA switched on");
    fPID->SetQAOn();
    fQA->Add(fPID->GetQAhistograms());
  }
  fPID->SetHasMCData(kTRUE);
  if(!fPIDdetectors.Length()) AddPIDdetector("TPC");
  fPID->InitializePID(fPIDdetectors.Data());     // Only restrictions to TPC allowed 

  // mcQA----------------------------------
  if (IsQAOn(kMCqa)) {
    AliInfo("MC QA on");
    if(!fMCQA) fMCQA = new AliHFEmcQA;
    if(!fHistMCQA) fHistMCQA = new TList();
    fHistMCQA->SetName("MCqa");
    fMCQA->CreateHistograms(AliHFEmcQA::kCharm,0,"mcqa_");               // create histograms for charm
    fMCQA->CreateHistograms(AliHFEmcQA::kBeauty,0,"mcqa_");              // create histograms for beauty
    fMCQA->CreateHistograms(AliHFEmcQA::kCharm,1,"mcqa_barrel_");        // create histograms for charm 
    fMCQA->CreateHistograms(AliHFEmcQA::kBeauty,1,"mcqa_barrel_");       // create histograms for beauty
    TIter next(gDirectory->GetList());
    TObject *obj;
    int counter = 0;
    TString objname;
    while ((obj = next.Next())) {
      objname = obj->GetName();
      TObjArray *toks = objname.Tokenize("_");
      if (toks->GetEntriesFast()){
        TObjString *fpart = (TObjString *)(toks->UncheckedAt(0));
        if ((fpart->String()).CompareTo("mcqa") == 0) fHistMCQA->AddAt(obj, counter++);
      }
    }
    fQA->Add(fHistMCQA);
  } 

  // secvtx----------------------------------
  if (IsSecVtxOn()) {
    AliInfo("Secondary Vertex Analysis on");
    fSecVtx = new AliHFEsecVtx;

    if(!fHistSECVTX) fHistSECVTX = new TList();
    fHistSECVTX->SetName("SecVtx");
    fSecVtx->CreateHistograms("secvtx_");
    TIter next(gDirectory->GetList());
    TObject *obj;
    int counter = 0;
    TString objname;
    while ((obj = next.Next())) {
      objname = obj->GetName();
      TObjArray *toks = objname.Tokenize("_");
      if (toks->GetEntriesFast()){
        TObjString *fpart = (TObjString *)(toks->UncheckedAt(0));
        if ((fpart->String()).CompareTo("secvtx") == 0) fHistSECVTX->AddAt(obj, counter++);
      }
    }
    fOutput->Add(fHistSECVTX);
  } 
}

//____________________________________________________________
void AliAnalysisTaskHFE::Exec(Option_t *){
  //
  // Run the analysis
  // 
  if(!fESD){
    AliError("No ESD Event");
    return;
  }
  if(!fMC){
    AliError("No MC Event");
    return;
  }
  fCFM->SetEventInfo(fMC);
  fPID->SetMCEvent(fMC);

  //fCFM->CheckEventCuts(AliCFManager::kEvtGenCuts, fMC);
  Double_t container[6];

  // Loop over the Monte Carlo tracks to see whether we have overlooked any track
  AliMCParticle *mctrack = 0x0;
  Int_t nElectrons = 0;

  if (IsSecVtxOn()) { 
    fSecVtx->SetEvent(fESD);
    fSecVtx->SetStack(fMC->Stack());
  }

  // run mc QA ------------------------------------------------
  if (IsQAOn(kMCqa)) {
    AliDebug(2, "Running MC QA");

    fMCQA->SetStack(fMC->Stack());
    fMCQA->Init();

    Int_t nPrims = fMC->Stack()->GetNprimary();
    Int_t nMCTracks = fMC->Stack()->GetNtrack();

    // loop over primary particles for quark and heavy hadrons
    for (Int_t igen = 0; igen < nPrims; igen++){
      fMCQA->GetQuarkKine(igen, AliHFEmcQA::kCharm);
      fMCQA->GetQuarkKine(igen, AliHFEmcQA::kBeauty);
      fMCQA->GetHadronKine(igen, AliHFEmcQA::kCharm);
      fMCQA->GetHadronKine(igen, AliHFEmcQA::kBeauty);
    }
    fMCQA->EndOfEventAna(AliHFEmcQA::kCharm);
    fMCQA->EndOfEventAna(AliHFEmcQA::kBeauty);

    // loop over all tracks for decayed electrons
    for (Int_t igen = 0; igen < nMCTracks; igen++){
      fMCQA->GetDecayedKine(igen, AliHFEmcQA::kCharm,  AliHFEmcQA::kElectronPDG, 0);
      fMCQA->GetDecayedKine(igen, AliHFEmcQA::kBeauty, AliHFEmcQA::kElectronPDG, 0);
      fMCQA->GetDecayedKine(igen, AliHFEmcQA::kCharm,  AliHFEmcQA::kElectronPDG, 1, kTRUE); // barrel
      fMCQA->GetDecayedKine(igen, AliHFEmcQA::kBeauty, AliHFEmcQA::kElectronPDG, 1, kTRUE); // barrel
    }

  } // end of MC QA loop
  // -----------------------------------------------------------------

  //
  // Loop MC
  //
  for(Int_t imc = fMC->GetNumberOfTracks(); imc--;){
    mctrack = fMC->GetTrack(imc);

    container[0] = mctrack->Pt();
    container[1] = mctrack->Eta();
    container[2] = mctrack->Phi();

    if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepMCGenerated, mctrack)) continue;
    fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepMCGenerated);
    (dynamic_cast<TH1F *>(fQA->At(2)))->Fill(mctrack->Phi() - TMath::Pi());
    if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepMCInAcceptance, mctrack)) continue;
    // find the label in the vector
    fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepMCInAcceptance);
    nElectrons++;
  }
  (dynamic_cast<TH1F *>(fQA->At(3)))->Fill(nElectrons);

  // fCFM->CheckEventCuts(AliCFManager::kEvtRecCuts, fESD);
        

  //
  // Loop ESD
  //

  Int_t nbtrackrec = fESD->GetNumberOfTracks();
  Int_t nElectronCandidates = 0;
  AliESDtrack *track = 0x0, *htrack = 0x0;
  Int_t pid = 0;
  // For double counted tracks
  Int_t *indexRecKineTPC = new Int_t[nbtrackrec];
  Int_t *indexRecKineITS = new Int_t[nbtrackrec];
  Int_t *indexRecPrim = new Int_t[nbtrackrec];
  Int_t *indexHFEcutsITS = new Int_t[nbtrackrec];
  Int_t *indexHFEcutsTPC = new Int_t[nbtrackrec];
  Int_t *indexHFEcutsTRD = new Int_t[nbtrackrec];
  Int_t *indexPid = new Int_t[nbtrackrec];

  for(Int_t nb = 0; nb < nbtrackrec; nb++){
    indexRecKineTPC[nb] = -1;
    indexRecKineITS[nb] = -1;
    indexRecPrim[nb]    = -1;
    indexHFEcutsITS[nb] = -1;
    indexHFEcutsTPC[nb] = -1;
    indexHFEcutsTRD[nb] = -1;
    indexPid[nb]        = -1;
  } 

  Bool_t alreadyseen = kFALSE;

  for(Int_t itrack = 0; itrack < fESD->GetNumberOfTracks(); itrack++){
    
    track = fESD->GetTrack(itrack);
          
    container[0] = track->Pt();
    container[1] = track->Eta();
    container[2] = track->Phi();

    Bool_t signal = kTRUE;
          
    // RecKine: TPC cuts  
    if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepRecKineTPC, track)) continue;
    fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepRecKineTPC);
    

    // Check if it is signal electrons
    if(!(mctrack = fMC->GetTrack(TMath::Abs(track->GetLabel())))) continue;
    
    container[3] = mctrack->Pt();
    container[4] = mctrack->Eta();
    container[5] = mctrack->Phi();
    
    if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepMCGenerated, mctrack)) signal = kFALSE;
    
    if(signal) {
      alreadyseen = kFALSE;
      for(Int_t nb = 0; nb < itrack; nb++){
        if(indexRecKineTPC[nb] == TMath::Abs(track->GetLabel())) alreadyseen = kTRUE;
      }
      indexRecKineTPC[itrack] = TMath::Abs(track->GetLabel());
      
      fCFM->GetParticleContainer()->Fill(container, (AliHFEcuts::kStepRecKineTPC + AliHFEcuts::kNcutESDSteps + 1));
      fCFM->GetParticleContainer()->Fill(&container[3], (AliHFEcuts::kStepRecKineTPC + 2*(AliHFEcuts::kNcutESDSteps + 1)));
      if(alreadyseen) {
        fCFM->GetParticleContainer()->Fill(&container[3], (AliHFEcuts::kStepRecKineTPC + 4*(AliHFEcuts::kNcutESDSteps + 1)));
      }
      else {
        fCFM->GetParticleContainer()->Fill(&container[3], (AliHFEcuts::kStepRecKineTPC + 3*(AliHFEcuts::kNcutESDSteps + 1)));
      }
    }

    
    // RecKine: ITS cuts
    if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepRecKineITS, track)) continue;
    fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepRecKineITS);
    if(signal) {
      alreadyseen = kFALSE;
      for(Int_t nb = 0; nb < itrack; nb++){
        if(indexRecKineITS[nb] == TMath::Abs(track->GetLabel())) alreadyseen = kTRUE;
      }
      indexRecKineITS[itrack] = TMath::Abs(track->GetLabel());
      fCFM->GetParticleContainer()->Fill(container, (AliHFEcuts::kStepRecKineITS + AliHFEcuts::kNcutESDSteps + 1));
      fCFM->GetParticleContainer()->Fill(&container[3], (AliHFEcuts::kStepRecKineITS + 2*(AliHFEcuts::kNcutESDSteps + 1)));
      if(alreadyseen) {
        fCFM->GetParticleContainer()->Fill(&container[3], (AliHFEcuts::kStepRecKineITS + 4*(AliHFEcuts::kNcutESDSteps + 1)));
      }
      else {
        fCFM->GetParticleContainer()->Fill(&container[3], (AliHFEcuts::kStepRecKineITS + 3*(AliHFEcuts::kNcutESDSteps + 1)));
      }
    }


    // Check TRD criterions (outside the correction framework)
    if(track->GetTRDncls()){
      (dynamic_cast<TH1F *>(fQA->At(6)))->Fill(track->GetTRDchi2()/track->GetTRDncls());
      (dynamic_cast<TH1F *>(fQA->At(1)))->Fill(track->GetAlpha());    // Check the acceptance without tight cuts
      (dynamic_cast<TProfile *>(fQA->At(4)))->Fill(container[0], track->GetTRDpidQuality());
      (dynamic_cast<TProfile *>(fQA->At(5)))->Fill(container[0], track->GetTRDncls());
    }

    
    // RecPrim
    if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepRecPrim, track)) continue;
    fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepRecPrim);
    if(signal) {
      fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepRecPrim + AliHFEcuts::kNcutESDSteps + 1);
      fCFM->GetParticleContainer()->Fill(&container[3], AliHFEcuts::kStepRecPrim + 2*(AliHFEcuts::kNcutESDSteps + 1));
      alreadyseen = kFALSE;
      for(Int_t nb = 0; nb < itrack; nb++){
        if(indexRecPrim[nb] == TMath::Abs(track->GetLabel())) alreadyseen = kTRUE;
      }
      indexRecPrim[itrack] = TMath::Abs(track->GetLabel());
      if(alreadyseen) {
        fCFM->GetParticleContainer()->Fill(&container[3], AliHFEcuts::kStepRecPrim + 4*(AliHFEcuts::kNcutESDSteps + 1));
      }
      else {
        fCFM->GetParticleContainer()->Fill(&container[3], AliHFEcuts::kStepRecPrim + 3*(AliHFEcuts::kNcutESDSteps + 1));
      }
    }


    // HFEcuts: ITS layers cuts
    if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepHFEcutsITS, track)) continue;
    fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepHFEcutsITS);
    if(signal) {
      fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepHFEcutsITS + AliHFEcuts::kNcutESDSteps + 1);
      fCFM->GetParticleContainer()->Fill(&container[3], AliHFEcuts::kStepHFEcutsITS + 2*(AliHFEcuts::kNcutESDSteps + 1));
      alreadyseen = kFALSE;
      for(Int_t nb = 0; nb < itrack; nb++){
        if(indexHFEcutsITS[nb] == TMath::Abs(track->GetLabel())) alreadyseen = kTRUE;
      }
      indexHFEcutsITS[itrack] = TMath::Abs(track->GetLabel());
      if(alreadyseen) {
        fCFM->GetParticleContainer()->Fill(&container[3], (AliHFEcuts::kStepHFEcutsITS + 4*(AliHFEcuts::kNcutESDSteps + 1)));
      }
      else {
        fCFM->GetParticleContainer()->Fill(&container[3], (AliHFEcuts::kStepHFEcutsITS + 3*(AliHFEcuts::kNcutESDSteps + 1)));
      }
    }

    // HFEcuts: TPC ratio clusters
    if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepHFEcutsTPC, track)) continue;
    fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepHFEcutsTPC);
    if(signal) {
      fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepHFEcutsTPC + AliHFEcuts::kNcutESDSteps + 1);
      fCFM->GetParticleContainer()->Fill(&container[3], AliHFEcuts::kStepHFEcutsTPC + 2*(AliHFEcuts::kNcutESDSteps + 1));
      alreadyseen = kFALSE;
      for(Int_t nb = 0; nb < itrack; nb++){
        if(indexHFEcutsTPC[nb] == TMath::Abs(track->GetLabel())) alreadyseen = kTRUE;
      }
      indexHFEcutsTPC[itrack] = TMath::Abs(track->GetLabel());    
      if(alreadyseen) {
        fCFM->GetParticleContainer()->Fill(&container[3], (AliHFEcuts::kStepHFEcutsTPC + 4*(AliHFEcuts::kNcutESDSteps + 1)));
      }
      else {
        fCFM->GetParticleContainer()->Fill(&container[3], (AliHFEcuts::kStepHFEcutsTPC + 3*(AliHFEcuts::kNcutESDSteps + 1)));
      }
    }
    
    // HFEcuts: Nb of tracklets TRD
    if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepHFEcutsTRD, track)) continue;
    fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepHFEcutsTRD);
    if(signal) {
      fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepHFEcutsTRD + AliHFEcuts::kNcutESDSteps + 1);
      fCFM->GetParticleContainer()->Fill(&container[3], AliHFEcuts::kStepHFEcutsTRD + 2*(AliHFEcuts::kNcutESDSteps + 1));
      alreadyseen = kFALSE;
      for(Int_t nb = 0; nb < itrack; nb++){
        if(indexHFEcutsTRD[nb] == TMath::Abs(track->GetLabel())) alreadyseen = kTRUE;
      }
      indexHFEcutsTRD[itrack] = TMath::Abs(track->GetLabel());
      if(alreadyseen) {
        fCFM->GetParticleContainer()->Fill(&container[3], (AliHFEcuts::kStepHFEcutsTRD + 4*(AliHFEcuts::kNcutESDSteps + 1)));
      }
      else {
        fCFM->GetParticleContainer()->Fill(&container[3], (AliHFEcuts::kStepHFEcutsTRD + 3*(AliHFEcuts::kNcutESDSteps + 1)));
      }
    }


    // track accepted, do PID
    if(!fPID->IsSelected(track)) continue;
    nElectronCandidates++;


    // Fill Containers
    fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepHFEcutsTRD + 1);
    if(signal) {
      fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepHFEcutsTRD + AliHFEcuts::kNcutESDSteps + 2);
      fCFM->GetParticleContainer()->Fill(&container[3], AliHFEcuts::kStepHFEcutsTRD + 2*(AliHFEcuts::kNcutESDSteps + 1)+1);
      alreadyseen = kFALSE;
      for(Int_t nb = 0; nb < itrack; nb++){
        if(indexPid[nb] == TMath::Abs(track->GetLabel())) alreadyseen = kTRUE;
      }
      indexPid[itrack] = TMath::Abs(track->GetLabel());
      if(alreadyseen) {
        fCFM->GetParticleContainer()->Fill(&container[3], (AliHFEcuts::kStepHFEcutsTRD + 1 + 4*(AliHFEcuts::kNcutESDSteps + 1)));
      }
      else {
        fCFM->GetParticleContainer()->Fill(&container[3], (AliHFEcuts::kStepHFEcutsTRD + 1 + 3*(AliHFEcuts::kNcutESDSteps + 1)));
      }
      // dimensions 3&4&5 : pt,eta,phi (MC)
      fCorrelation->Fill(container);
    }


    // Track selected: distinguish between true and fake
    if((pid = TMath::Abs(mctrack->Particle()->GetPdgCode())) == 11){
      // pair analysis [mj]
      if (IsSecVtxOn()) {
        AliDebug(2, "Running Secondary Vertex Analysis");
        fSecVtx->InitAnaPair();
        for(Int_t jtrack = 0; jtrack < fESD->GetNumberOfTracks(); jtrack++){
          htrack = fESD->GetTrack(jtrack);
          if ( itrack == jtrack ) continue;
          //if( fPID->IsSelected(htrack) && (itrack < jtrack)) continue;
          if( abs(fSecVtx->GetMCPID(track)) == 11 && (itrack < jtrack)) continue;
          fSecVtx->AnaPair(track, htrack, jtrack);
        }
        // based on the partner of e info, you run secandary vertexing function
        fSecVtx->RunSECVTX(track);
      } // end of pair analysis
    } else {
      // Fill THnSparse with the information for Fake Electrons
      switch(pid){
      case 13:    container[3] = AliPID::kMuon;
      case 211:   container[3] = AliPID::kPion;
      case 321:   container[3] = AliPID::kKaon;
      case 2212:  container[3] = AliPID::kProton;
      }
      fFakeElectrons->Fill(container);
    }
  }
  fNEvents->Fill(1);
  fNElectronTracksEvent->Fill(nElectronCandidates);

  // release memory
  delete[] indexRecKineTPC;
  delete[] indexRecKineITS;
  delete[] indexRecPrim;
  delete[] indexHFEcutsITS;
  delete[] indexHFEcutsTPC;
  delete[] indexHFEcutsTRD;
  delete[] indexPid;
  
  // Done!!!
  PostData(0, fNEvents);
  PostData(1, fOutput);
  PostData(2, fQA);
}

//____________________________________________________________
void AliAnalysisTaskHFE::Terminate(Option_t *){
  //
  // Terminate not implemented at the moment
  //
}

//____________________________________________________________
void AliAnalysisTaskHFE::MakeParticleContainer(){
  //
  // Create the particle container for the correction framework manager and 
  // link it
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
  AliCFContainer* container = new AliCFContainer("container","container for tracks", (AliHFEcuts::kNcutSteps + 1 + 4*(AliHFEcuts::kNcutESDSteps + 1)), kNvar, iBin);

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

  fCorrelation = new THnSparseF("correlation","THnSparse with correlations",2*kNvar,thnDim);
  for (int k=0; k<kNvar; k++) {
    fCorrelation->SetBinEdges(k,binEdges[k]);
    fCorrelation->SetBinEdges(k+kNvar,binEdges[k]);
  }
  fCorrelation->Sumw2();

  // Add a histogram for Fake electrons
  thnDim[3] = AliPID::kSPECIES;
  fFakeElectrons = new THnSparseF("fakeEkectrons", "Output for Fake Electrons", kNvar + 1, thnDim);
  for(Int_t idim = 0; idim < kNvar; idim++)
    fFakeElectrons->SetBinEdges(idim, binEdges[idim]);
}

void AliAnalysisTaskHFE::AddPIDdetector(Char_t *detector){
  if(!fPIDdetectors.Length()) 
    fPIDdetectors = detector;
  else
    fPIDdetectors += Form(":%s", detector);
}

//____________________________________________________________
void AliAnalysisTaskHFE::PrintStatus(){
  printf("\n\tAnalysis Settings\n\t========================================\n\n");
  printf("\tSecondary Vertex finding: %s\n", IsSecVtxOn() ? "YES" : "NO");
  printf("\tPrimary Vertex resolution: %s\n", IsPriVtxOn() ? "YES" : "NO");
  printf("\n");
  printf("\tParticle Identification Detectors:\n");
  TObjArray *detectors = fPIDdetectors.Tokenize(":");
  for(Int_t idet = 0; idet < detectors->GetEntries(); idet++)
    printf("\t\t%s\n", (dynamic_cast<TObjString *>(detectors->At(idet)))->String().Data());
  printf("\n");
  printf("\tQA: \n");
  printf("\t\tPID: %s\n", IsQAOn(kPIDqa) ? "YES" :  "NO");
  printf("\t\tCUTS: %s\n", IsQAOn(kCUTqa) ? "YES" : "NO");
  printf("\t\tMC: %s\n", IsQAOn(kMCqa) ? "YES" : "NO");
  printf("\n");
}
