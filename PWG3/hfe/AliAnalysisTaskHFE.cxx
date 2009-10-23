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
 *  Matus Kalisky <matus.kalisky@cern.ch>
 *  MinJung Kweon <minjung@physi.uni-heidelberg.de>
 */
#include <TAxis.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TDirectory.h>
#include <TFile.h>
#include <TH1D.h>
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
#include "AliCFEffGrid.h"
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
  , fPIDstrategy(0)
  , fESD(0x0)
  , fMC(0x0)
  , fCFM(0x0)
  , fCorrelation(0x0)
  , fPIDperformance(0x0)
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
  fPID = new AliHFEpid;

  SetHasMCData();
}

//____________________________________________________________
AliAnalysisTaskHFE::AliAnalysisTaskHFE(const char * name):
  AliAnalysisTask(name, "")
  , fQAlevel(0)
  , fPIDdetectors("")
  , fPIDstrategy(0)
  , fESD(0x0)
  , fMC(0x0)
  , fCFM(0x0)
  , fCorrelation(0x0)
  , fPIDperformance(0x0)
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
  fPID = new AliHFEpid;

  SetHasMCData();
}

//____________________________________________________________
AliAnalysisTaskHFE::AliAnalysisTaskHFE(const AliAnalysisTaskHFE &ref):
  AliAnalysisTask(ref)
  , fQAlevel(ref.fQAlevel)
  , fPIDdetectors(ref.fPIDdetectors)
  , fPIDstrategy(ref.fPIDstrategy)
  , fESD(ref.fESD)
  , fMC(ref.fMC)
  , fCFM(ref.fCFM)
  , fCorrelation(ref.fCorrelation)
  , fPIDperformance(ref.fPIDperformance)
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
  fPIDstrategy = ref.fPIDstrategy;
  fESD = ref.fESD;
  fMC = ref.fMC;
  fCFM = ref.fCFM;
  fCorrelation = ref.fCorrelation;
  fPIDperformance = ref.fPIDperformance;
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
  if(fSecVtx) delete fSecVtx;
  if(fMCQA) delete fMCQA;
  if(fNEvents) delete fNEvents;
  if(fCorrelation){
    fCorrelation->Clear();
    delete fCorrelation;
  }
  if(fPIDperformance) delete fPIDperformance;
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
  if(!fCuts){
    AliWarning("Cuts not available. Default cuts will be used");
    fCuts = new AliHFEcuts;
    fCuts->CreateStandardCuts();
  }
  fCuts->Initialize(fCFM);
  if(fCuts->IsInDebugMode()) fQA->Add(fCuts->GetQAhistograms());
 
  // add output objects to the List
  fOutput->AddAt(fCFM->GetParticleContainer(), 0);
  fOutput->AddAt(fCorrelation, 1);
  fOutput->AddAt(fPIDperformance, 2);
  fOutput->AddAt(fNElectronTracksEvent, 3);

  // Initialize PID
  if(IsQAOn(kPIDqa)){
    AliInfo("PID QA switched on");
    //fPID->SetDebugLevel(2);
    fPID->SetQAOn();
    fQA->Add(fPID->GetQAhistograms());
  }
  fPID->SetHasMCData(HasMCData());
  if(!fPIDdetectors.Length() && ! fPIDstrategy) AddPIDdetector("TPC");
  if(fPIDstrategy)
    fPID->InitializePID(Form("Strategy%d", fPIDstrategy));
  else
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
    fMCQA->CreateHistograms(AliHFEmcQA::kCharm,2,"mcqa_unitY_");         // create histograms for charm 
    fMCQA->CreateHistograms(AliHFEmcQA::kBeauty,2,"mcqa_unitY_");        // create histograms for beauty
    fMCQA->CreateHistograms(AliHFEmcQA::kCharm,3,"mcqa_reccut_");        // create histograms for charm 
    fMCQA->CreateHistograms(AliHFEmcQA::kBeauty,3,"mcqa_reccut_");       // create histograms for beauty
    fMCQA->CreateHistograms(AliHFEmcQA::kCharm,4,"mcqa_recpidcut_");     // create histograms for charm 
    fMCQA->CreateHistograms(AliHFEmcQA::kBeauty,4,"mcqa_recpidcut_");    // create histograms for beauty
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
  if(!fCuts){
    AliError("HFE cuts not available");
    return;
  }

  //fCFM->CheckEventCuts(AliCFManager::kEvtGenCuts, fMC);
  Double_t container[6];
  // container for the output THnSparse
  Double_t dataE[5]; // [pT, eta, Phi, type, 'C' or 'B']

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

    Int_t nMCTracks = fMC->Stack()->GetNtrack();

    // loop over all tracks for decayed electrons
    for (Int_t igen = 0; igen < nMCTracks; igen++){
      TParticle* mcpart = fMC->Stack()->Particle(igen);
      fMCQA->GetQuarkKine(mcpart, igen, AliHFEmcQA::kCharm);
      fMCQA->GetQuarkKine(mcpart, igen, AliHFEmcQA::kBeauty);
      fMCQA->GetHadronKine(mcpart, AliHFEmcQA::kCharm);
      fMCQA->GetHadronKine(mcpart, AliHFEmcQA::kBeauty);
      fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kCharm,  AliHFEmcQA::kElectronPDG, 0); // no accept cut
      fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kBeauty, AliHFEmcQA::kElectronPDG, 0); // no accept cut
      if (TMath::Abs(mcpart->Eta()) < 0.9) {
        fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kCharm,  AliHFEmcQA::kElectronPDG, 1); // accept |eta|<0.9
        fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kBeauty, AliHFEmcQA::kElectronPDG, 1); // accept |eta|<0.9
      }
      if (TMath::Abs(GetRapidity(mcpart)) < 0.5) {
        fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kCharm,  AliHFEmcQA::kElectronPDG, 2); // accept |y|<0.5
        fMCQA->GetDecayedKine(mcpart, AliHFEmcQA::kBeauty, AliHFEmcQA::kElectronPDG, 2); // accept |y|<0.5
      }
    }
    fMCQA->EndOfEventAna(AliHFEmcQA::kCharm);
    fMCQA->EndOfEventAna(AliHFEmcQA::kBeauty);

  } // end of MC QA loop
  // -----------------------------------------------------------------

  //
  // Loop MC
  //
  fCuts->SetEventInfo(fMC);
  for(Int_t imc = fMC->GetNumberOfTracks(); imc--;){
    mctrack = dynamic_cast<AliMCParticle *>(fMC->GetTrack(imc));

    container[0] = mctrack->Pt();
    container[1] = mctrack->Eta();
    container[2] = mctrack->Phi();

    if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepMCGenerated, mctrack)) continue;
    fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepMCGenerated);
    if(IsSignalElectron(mctrack)) fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepMCsignal);
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
  Int_t nElectronCandidates = 0;
  AliESDtrack *track = 0x0, *htrack = 0x0;
  Int_t pid = 0;
  // For double counted tracks
  LabelContainer cont(fESD->GetNumberOfTracks());
  Bool_t alreadyseen = kFALSE;

  Bool_t signal = kTRUE;

  fCuts->SetEventInfo(fESD);
  for(Int_t itrack = 0; itrack < fESD->GetNumberOfTracks(); itrack++){
    
    track = fESD->GetTrack(itrack);
          
    container[0] = track->Pt();
    container[1] = track->Eta();
    container[2] = track->Phi();

    dataE[0] = track->Pt();
    dataE[1] = track->Eta();
    dataE[2] = track->Phi();
    dataE[3] = -1;
    dataE[4] = -1;

    signal = kTRUE;
          
    // RecKine: ITSTPC cuts  
    if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepRecKineITSTPC, track)) continue;
    
    // Check if it is electrons near the vertex
    if(!(mctrack = dynamic_cast<AliMCParticle *>(fMC->GetTrack(TMath::Abs(track->GetLabel()))))) continue;
    TParticle* mctrack4QA = fMC->Stack()->Particle(TMath::Abs(track->GetLabel()));

    container[3] = mctrack->Pt();
    container[4] = mctrack->Eta();
    container[5] = mctrack->Phi();
    
    if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepMCGenerated, mctrack)) signal = kFALSE;
    
    if(signal) {
      alreadyseen = cont.Find(TMath::Abs(track->GetLabel()));
      cont.Append(TMath::Abs(track->GetLabel()));
      
      fCFM->GetParticleContainer()->Fill(&container[3], AliHFEcuts::kStepRecKineITSTPC);
      fCFM->GetParticleContainer()->Fill(&container[0], (AliHFEcuts::kStepRecKineITSTPC + 2*(AliHFEcuts::kNcutESDSteps + 1)));
      if(alreadyseen) {
        fCFM->GetParticleContainer()->Fill(&container[3], (AliHFEcuts::kStepRecKineITSTPC + (AliHFEcuts::kNcutESDSteps + 1)));
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
    if(signal) {
      fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepRecPrim + 2*(AliHFEcuts::kNcutESDSteps + 1));
      fCFM->GetParticleContainer()->Fill(&container[3], AliHFEcuts::kStepRecPrim);
      if(alreadyseen) {
        fCFM->GetParticleContainer()->Fill(&container[3], AliHFEcuts::kStepRecPrim + AliHFEcuts::kNcutESDSteps + 1);
      }
    }


    // HFEcuts: ITS layers cuts
    if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepHFEcutsITS, track)) continue;
    if(signal) {
      fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepHFEcutsITS + 2*(AliHFEcuts::kNcutESDSteps + 1));
      fCFM->GetParticleContainer()->Fill(&container[3], AliHFEcuts::kStepHFEcutsITS);
      if(alreadyseen) {
        fCFM->GetParticleContainer()->Fill(&container[3], AliHFEcuts::kStepHFEcutsITS + AliHFEcuts::kNcutESDSteps + 1);
      }
    }

    // HFEcuts: Nb of tracklets TRD0
    if(!fCFM->CheckParticleCuts(AliHFEcuts::kStepHFEcutsTRD, track)) continue;
    if(signal) {
      fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepHFEcutsTRD + 2*(AliHFEcuts::kNcutESDSteps + 1));
      fCFM->GetParticleContainer()->Fill(&container[3], AliHFEcuts::kStepHFEcutsTRD);
      if(alreadyseen) {
        fCFM->GetParticleContainer()->Fill(&container[3], (AliHFEcuts::kStepHFEcutsTRD + (AliHFEcuts::kNcutESDSteps + 1)));
      }
      // dimensions 3&4&5 : pt,eta,phi (MC)
      ((THnSparseF *)fCorrelation->At(0))->Fill(container);
    }

   if (IsQAOn(kMCqa)) {
     // mc qa for after the reconstruction cuts  
     AliDebug(2, "Running MC QA");
     fMCQA->GetDecayedKine(mctrack4QA, AliHFEmcQA::kCharm,  AliHFEmcQA::kElectronPDG, 3);  // charm
     fMCQA->GetDecayedKine(mctrack4QA, AliHFEmcQA::kBeauty,  AliHFEmcQA::kElectronPDG, 3); // beauty 
    } 

    // track accepted, do PID
    AliHFEpidObject hfetrack;
    hfetrack.fAnalysisType = AliHFEpidObject::kESDanalysis;
    hfetrack.fRecTrack = track;
    if(HasMCData()) hfetrack.fMCtrack = mctrack;
    if(!fPID->IsSelected(&hfetrack)) continue;
    nElectronCandidates++;

   if (IsQAOn(kMCqa)) {
     // mc qa for after the reconstruction and pid cuts  
     AliDebug(2, "Running MC QA");
     fMCQA->GetDecayedKine(mctrack4QA, AliHFEmcQA::kCharm,  AliHFEmcQA::kElectronPDG, 4);  // charm
     fMCQA->GetDecayedKine(mctrack4QA, AliHFEmcQA::kBeauty,  AliHFEmcQA::kElectronPDG, 4); // beauty 
    } 

    // Fill Containers
    if(signal) {
      fCFM->GetParticleContainer()->Fill(container, AliHFEcuts::kStepHFEcutsTRD +1 + 2*(AliHFEcuts::kNcutESDSteps + 1));
      fCFM->GetParticleContainer()->Fill(&container[3], AliHFEcuts::kStepHFEcutsTRD + 1);
      if(alreadyseen) {
        fCFM->GetParticleContainer()->Fill(&container[3], (AliHFEcuts::kStepHFEcutsTRD + 1 + (AliHFEcuts::kNcutESDSteps + 1)));
      }
      // dimensions 3&4&5 : pt,eta,phi (MC)
      ((THnSparseF *)fCorrelation->At(1))->Fill(container);
    }

    // Track selected: distinguish between true and fake
    AliDebug(1, Form("Candidate Selected, filling THnSparse, PID: %d\n", mctrack->Particle()->GetPdgCode()));
    if((pid = TMath::Abs(mctrack->Particle()->GetPdgCode())) == 11){
      Int_t type = IsSignalElectron(track);
      AliDebug(1, Form("Type: %d\n", type));
      if(type){
	      dataE[4] = type; // beauty[1] or charm[2]
	      dataE[3] = 2;  // signal electron
      }
      else{
	      dataE[3] = 1; // not a signal electron
	      dataE[4] = 0;
      }
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
    } 
    else {
      // Fill THnSparse with the information for Fake Electrons
      dataE[3] = 0;
      dataE[4] = 0;
    }
    // fill the performance THnSparse, if the mc origin could be defined
    if(dataE[3] > -1){
      AliDebug(1, Form("Entries: [%.3f|%.3f|%.3f|%f|%f]\n", dataE[0],dataE[1],dataE[2],dataE[3],dataE[4]));
      fPIDperformance->Fill(dataE);
    }
  }
  
 
  fNEvents->Fill(1);
  fNElectronTracksEvent->Fill(nElectronCandidates);
  
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
  if(IsRunningPostProcess()){
    fOutput = dynamic_cast<TList *>(GetOutputData(1));
    if(!fOutput){
      AliError("Results not available");
      return;
    }
    PostProcess();
  }
}

//____________________________________________________________
void AliAnalysisTaskHFE::Load(TString filename){
  //
  // Load Results into the task
  //
  fQA = NULL; fOutput = NULL; fNEvents = NULL;
  TFile *input = TFile::Open(filename.Data());
  if(!input || input->IsZombie()){
    AliError("Cannot read file");
    return;
  }
  TH1 *htmp = dynamic_cast<TH1I *>(input->Get("nEvents"));
  if(htmp)
    fNEvents = dynamic_cast<TH1I *>(htmp->Clone());
  else
    AliError("Event Counter histogram not found"); 
  TList *ltmp = dynamic_cast<TList *>(input->Get("Results"));
  if(ltmp)
    fOutput = dynamic_cast<TList *>(ltmp->Clone());
  else
    AliError("Output Histograms not found");
  ltmp = dynamic_cast<TList *>(input->Get("QA"));
  if(ltmp)
    fQA = dynamic_cast<TList *>(ltmp->Clone());
  else
    AliError("QA histograms not found");
  input->Close();
  delete input;
  Int_t nObjects = 0;
  if(fNEvents) nObjects++;  
  if(fOutput) nObjects++;
  if(fQA) nObjects++;
  AliInfo(Form("Loaded %d Objects into task", nObjects));
}

//____________________________________________________________
void AliAnalysisTaskHFE::PostProcess(){
  //
  // Plotting Ratio histograms
  // + All electrons / all candidates (Purity for Electrons)
  // + All signal electrons / all electrons (Purity for signals)
  // For this the following pt-histograms have to be projected from the THnSparse
  // + All Electron candidates
  // + All Real electrons
  // + All Signal Electrons
  // + All misidentified electrons
  // Additionally we draw Efficiency histograms:
  // + Reconstructible Electrons
  // + Signal Electrons
  // + Selected Electrons
  // + Selected Electrons in Acceptance
  //

  fPIDperformance = dynamic_cast<THnSparseF *>(fOutput->FindObject("PIDperformance"));
  if(!fPIDperformance){
    AliError("PID Performance histogram not found in the output container");
    return;
  }
  // Make projection
  // always project to pt dimension
  // get the histograms under our control
  TH1 *allCandidates = NULL, *allElectrons = NULL, *allSignals = NULL, *allFakes = NULL;
  allCandidates = fPIDperformance->Projection(0);
  allCandidates->SetName("hAllCandidates");
  allCandidates->SetTitle("All Candidates");
  allCandidates->Sumw2();
  Int_t firstDim3 = fPIDperformance->GetAxis(3)->GetFirst(), lastDim3 = fPIDperformance->GetAxis(3)->GetLast();
  fPIDperformance->GetAxis(3)->SetRange(firstDim3 + 1, lastDim3);
  allElectrons = fPIDperformance->Projection(0);
  allElectrons->Sumw2();
  allElectrons->SetName("hAllElectrons");
  allElectrons->SetTitle("All Electrons");
  Int_t firstDim4 = fPIDperformance->GetAxis(4)->GetFirst(), lastDim4 = fPIDperformance->GetAxis(4)->GetLast();
  fPIDperformance->GetAxis(4)->SetRange(firstDim4 + 1, lastDim4);
  allSignals = fPIDperformance->Projection(0);
  allSignals->Sumw2();
  allSignals->SetName("hAllSignals");
  allSignals->SetTitle("All Signal Electrons");
  fPIDperformance->GetAxis(4)->SetRange(firstDim4, lastDim4);       // Reset 4th axis
  fPIDperformance->GetAxis(3)->SetRange(firstDim3, firstDim3);  // Select fakes
  allFakes = fPIDperformance->Projection(0);
  allFakes->Sumw2();
  allFakes->SetName("hAllFakes");
  allFakes->SetTitle("All Fakes");
  fPIDperformance->GetAxis(3)->SetRange(firstDim3, lastDim3);       // Reset also 3rd axis

  // Make Ratios
  TH1D *electronPurity = dynamic_cast<TH1D *>(allElectrons->Clone());
  electronPurity->Divide(allCandidates);
  electronPurity->SetName("electronPurity");
  electronPurity->SetTitle("Electron Purity");
  electronPurity->GetXaxis()->SetTitle("p_{T} / GeV/c");
  electronPurity->GetYaxis()->SetTitle("Purity / %");
  TH1D *signalPurity = dynamic_cast<TH1D *>(allSignals->Clone());
  signalPurity->Divide(allElectrons);
  signalPurity->SetName("signalPurity");
  signalPurity->SetTitle("Purity of Electrons coming from Heavy flavours");
  signalPurity->GetXaxis()->SetTitle("p_{T} / GeV/c");
  signalPurity->GetYaxis()->SetTitle("Purity / %");
  TH1D *fakeContamination = dynamic_cast<TH1D *>(allFakes->Clone());
  fakeContamination->Divide(allCandidates);
  fakeContamination->SetName("fakeContamination");
  fakeContamination->SetTitle("Contamination of misidentified hadrons");
  fakeContamination->GetXaxis()->SetTitle("p_{T} / GeV/c");
  fakeContamination->GetYaxis()->SetTitle("Purity / %");

  // Graphics settings
  const Int_t nHistos = 7;
  TH1 *hptr[7] = {allCandidates, allElectrons, allSignals, allFakes, electronPurity, signalPurity, fakeContamination}, *histo;
  for(Int_t ihist = 0; ihist < nHistos; ihist++){
    (histo = hptr[ihist])->SetStats(kFALSE);
    histo->SetLineColor(kBlack);
    histo->SetMarkerColor(kBlue);
    histo->SetMarkerStyle(22);
  }

  // Make percent output
  electronPurity->Scale(100);
  signalPurity->Scale(100);
  fakeContamination->Scale(100);
  
  // Draw output
  TCanvas *cSpectra = new TCanvas("cSpectra","pt Spectra", 800, 600);
  cSpectra->Divide(2,2);
  TCanvas *cRatios = new TCanvas("cRatios", "Ratio Plots", 800, 600);
  cRatios->Divide(2,2);
  cSpectra->cd(1);
  gPad->SetLogy();
  allCandidates->GetXaxis()->SetTitle("p_{T} / GeV/c");
  allCandidates->GetYaxis()->SetTitle("dN/dp_{T} 1/GeV/c");
  allCandidates->Draw("e");
  cSpectra->cd(2);
  gPad->SetLogy();
  allElectrons->GetXaxis()->SetTitle("p_{T} / GeV/c");
  allElectrons->GetYaxis()->SetTitle("dN/dp_{T} 1/GeV/c");
  allElectrons->Draw("e");
  cSpectra->cd(3);
  gPad->SetLogy();
  allSignals->GetXaxis()->SetTitle("p_{T} / GeV/c");
  allSignals->GetYaxis()->SetTitle("dN/dp_{T} 1/GeV/c");
  allSignals->Draw("e");
  cSpectra->cd(4);
  gPad->SetLogy();
  allFakes->GetXaxis()->SetTitle("p_{T} / GeV/c");
  allFakes->GetYaxis()->SetTitle("dN/dp_{T} 1/GeV/c");
  allFakes->Draw("e");
  cRatios->cd(1);
  electronPurity->Draw("e");
  cRatios->cd(2);
  signalPurity->Draw("e");
  cRatios->cd(3);
  fakeContamination->Draw("e");

  // Now draw the Efficiency
  // We show: 
  // + InAcceptance / Generated
  // + Signal / Generated
  // + Selected / Generated
  // + Selected / InAcceptance (Reconstructible)
  TCanvas *cEff = new TCanvas("cEff", "Efficiency", 800, 600);
  cEff->Divide(2,2);
  AliCFContainer *effc = dynamic_cast<AliCFContainer *>(fOutput->At(0));
  if(!effc){
    AliError("Efficiency Container not found in Output Container");
    return;
  }
  AliCFEffGrid *effCalc = new AliCFEffGrid("effCalc", "Efficiency Calculation Grid", *effc);
  effCalc->CalculateEfficiency(AliHFEcuts::kStepMCInAcceptance, AliHFEcuts::kStepMCGenerated);
  TH1 *effReconstructibleP = effCalc->Project(0);
  effReconstructibleP->SetName("effReconstructibleP");
  effReconstructibleP->SetTitle("Efficiency of reconstructible tracks");
  effReconstructibleP->GetXaxis()->SetTitle("p_{T} / GeV/c");
  effReconstructibleP->GetYaxis()->SetTitle("Efficiency");
  effReconstructibleP->GetYaxis()->SetRangeUser(0.,1.);
  effReconstructibleP->SetMarkerStyle(22);
  effReconstructibleP->SetMarkerColor(kBlue);
  effReconstructibleP->SetLineColor(kBlack);
  effReconstructibleP->SetStats(kFALSE);
  cEff->cd(1);
  effReconstructibleP->Draw("e");
  effCalc->CalculateEfficiency(AliHFEcuts::kStepMCsignal, AliHFEcuts::kStepMCGenerated);
  TH1 *effSignal = effCalc->Project(0);
  effSignal->SetName("effSignal");
  effSignal->SetTitle("Efficiency of Signal Electrons");
  effSignal->GetXaxis()->SetTitle("p_{T} / GeV/c");
  effSignal->GetYaxis()->SetTitle("Efficiency");
  effSignal->GetYaxis()->SetRangeUser(0., 1.);
  effSignal->SetMarkerStyle(22);
  effSignal->SetMarkerColor(kBlue);
  effSignal->SetLineColor(kBlack);
  effSignal->SetStats(kFALSE);
  cEff->cd(2);
  effSignal->Draw("e");
  effCalc->CalculateEfficiency(AliHFEcuts::kStepHFEcutsTRD + 1, AliHFEcuts::kStepMCGenerated);
  TH1 *effPIDP = effCalc->Project(0);
  effPIDP->SetName("effPIDP");
  effPIDP->SetTitle("Efficiency of selected tracks");
  effPIDP->GetXaxis()->SetTitle("p_{T} / GeV/c");
  effPIDP->GetYaxis()->SetTitle("Efficiency");
  effPIDP->GetYaxis()->SetRangeUser(0.,1.);
  effPIDP->SetMarkerStyle(22);
  effPIDP->SetMarkerColor(kBlue);
  effPIDP->SetLineColor(kBlack);
  effPIDP->SetStats(kFALSE);
  cEff->cd(3);
  effPIDP->Draw("e");
  effCalc->CalculateEfficiency(AliHFEcuts::kStepHFEcutsTRD + 1, AliHFEcuts::kStepMCInAcceptance);
  TH1 *effPIDAcc = effCalc->Project(0);
  effPIDAcc->SetName("effPIDAcc");
  effPIDAcc->SetTitle("Efficiency of selected tracks in acceptance");
  effPIDAcc->GetXaxis()->SetTitle("p_{T} / GeV/c");
  effPIDAcc->GetYaxis()->SetTitle("Efficiency");
  effPIDAcc->GetYaxis()->SetRangeUser(0.,1.);
  effPIDAcc->SetMarkerStyle(22);
  effPIDAcc->SetMarkerColor(kBlue);
  effPIDAcc->SetLineColor(kBlack);
  effPIDAcc->SetStats(kFALSE);
  cEff->cd(4);
  effPIDAcc->Draw("e");
  delete effCalc;

  // cleanup
  //delete allCandidates; delete allElectrons; delete allSignals; delete allFakes;
  //delete electronPurity; delete signalPurity; delete fakeContamination;
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
  AliCFContainer* container = new AliCFContainer("container","container for tracks", (AliHFEcuts::kNcutSteps + 1 + 2*(AliHFEcuts::kNcutESDSteps + 1)), kNvar, iBin);

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

  if(!fCorrelation) fCorrelation = new TList();
  fCorrelation->SetName("correlation");

  THnSparseF *correlation0 = new THnSparseF("correlationstepbeforePID","THnSparse with correlations",2*kNvar,thnDim);
  THnSparseF *correlation1 = new THnSparseF("correlationstepafterPID","THnSparse with correlations",2*kNvar,thnDim);
  for (int k=0; k<kNvar; k++) {
    correlation0->SetBinEdges(k,binEdges[k]);
    correlation0->SetBinEdges(k+kNvar,binEdges[k]);
    correlation1->SetBinEdges(k,binEdges[k]);
    correlation1->SetBinEdges(k+kNvar,binEdges[k]);
  }
  correlation0->Sumw2();
  correlation1->Sumw2();
  
  fCorrelation->AddAt(correlation0,0);
  fCorrelation->AddAt(correlation1,1);
  
  // Add a histogram for Fake electrons
  const Int_t nDim=5;
  Int_t nBin[nDim] = {40, 8, 18, 3, 3};
  Double_t* binEdges2[nDim];
  for(Int_t ivar = 0; ivar < nDim; ivar++)
    binEdges2[ivar] = new Double_t[nBin[ivar] + 1];

  //values for bin lower bounds
  for(Int_t i=0; i<=nBin[0]; i++) binEdges2[0][i]=(Double_t)TMath::Power(10,TMath::Log10(kPtmin) + (TMath::Log10(kPtmax)-TMath::Log10(kPtmin))/nBin[0]*(Double_t)i);  
  for(Int_t i=0; i<=nBin[1]; i++) binEdges2[1][i]=(Double_t)kEtamin  + (kEtamax-kEtamin)/nBin[1]*(Double_t)i;
  for(Int_t i=0; i<=nBin[2]; i++) binEdges2[2][i]=(Double_t)kPhimin  + (kPhimax-kPhimin)/nBin[2]*(Double_t)i;
  for(Int_t i=0; i<=nBin[3]; i++) binEdges2[3][i] = i;
  for(Int_t i=0; i<=nBin[4]; i++) binEdges2[4][i] = i;

  fPIDperformance = new THnSparseF("PIDperformance", "PID performance; pT [GeV/c]; theta [rad]; phi [rad] type (0 - not el, 1 - other el, 2 - HF el; flavor (0 - no, 1 - charm, 2 - bottom)", nDim, nBin);
  for(Int_t idim = 0; idim < nDim; idim++)
    fPIDperformance->SetBinEdges(idim, binEdges2[idim]);
}

void AliAnalysisTaskHFE::AddPIDdetector(Char_t *detector){
  if(!fPIDdetectors.Length()) 
    fPIDdetectors = detector;
  else
    fPIDdetectors += Form(":%s", detector);
}

//____________________________________________________________
void AliAnalysisTaskHFE::PrintStatus(){
  //
  // Print Analysis status
  //
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
  printf("\t\tCUTS: %s\n", (fCuts != NULL && fCuts->IsInDebugMode()) ? "YES" : "NO");
  printf("\t\tMC: %s\n", IsQAOn(kMCqa) ? "YES" : "NO");
  printf("\n");
}

//____________________________________________________________
AliAnalysisTaskHFE::LabelContainer::LabelContainer(Int_t capacity):
  fContainer(NULL),
  fBegin(NULL),
  fEnd(NULL),
  fLast(NULL),
  fCurrent(NULL)
{
  //
  // Default constructor
  //
  fContainer = new Int_t[capacity];
  fBegin = &fContainer[0];
  fEnd = &fContainer[capacity - 1];
  fLast = fCurrent = fBegin;
}

//____________________________________________________________
Bool_t AliAnalysisTaskHFE::LabelContainer::Append(Int_t label){
  //
  // Add Label to the container
  //
  if(fLast > fEnd) return kFALSE;
  *fLast++ = label;
  return kTRUE;
}

//____________________________________________________________
Bool_t AliAnalysisTaskHFE::LabelContainer::Find(Int_t label){
  //
  // Find track in the list of labels
  //
  for(Int_t *entry = fBegin; entry <= fLast; entry++) 
    if(*entry == label) return kTRUE;
  return kFALSE;
}

//____________________________________________________________
Int_t AliAnalysisTaskHFE::LabelContainer::Next(){
  //
  // Mimic iterator
  //
  if(fCurrent > fLast) return -1; 
  return *fCurrent++;
}

//____________________________________________________________
Int_t AliAnalysisTaskHFE::IsSignalElectron(AliVParticle *fTrack) const{
  //
  // Checks whether the identified electron track is coming from heavy flavour
  // returns 0 in case of no signal, 1 in case of charm and 2 in case of Bottom
  //
  enum{
    kNoSignal = 0,
    kCharm = 1,
    kBeauty = 2
  };
  AliMCParticle *mctrack = NULL;
  TString objname = fTrack->IsA()->GetName();
  if(!objname.CompareTo("AliESDtrack")){
    AliDebug(2, "Checking signal for ESD track");
    AliESDtrack *esdtrack = dynamic_cast<AliESDtrack *>(fTrack);
    mctrack = dynamic_cast<AliMCParticle *>(fMC->GetTrack(TMath::Abs(esdtrack->GetLabel())));
  }
  else if(!objname.CompareTo("AliAODtrack")){
    AliError("AOD Analysis not implemented yet");
    return kNoSignal;
  }
  else if(!objname.CompareTo("AliMCParticle")){
    AliDebug(2, "Checking signal for MC track");
    mctrack = dynamic_cast<AliMCParticle *>(fTrack);
  }
  else{
    AliError("Input object not supported");
    return kNoSignal;
  }
  if(!mctrack) return kNoSignal;
  TParticle *ecand = mctrack->Particle(); 
  if(TMath::Abs(ecand->GetPdgCode()) != 11) return kNoSignal; // electron candidate not true electron
  Int_t motherLabel = TMath::Abs(ecand->GetFirstMother());
  AliDebug(3, Form("mother label: %d\n", motherLabel));
  if(!motherLabel) return kNoSignal; // mother track unknown
  AliMCParticle *motherTrack = dynamic_cast<AliMCParticle *>(fMC->GetTrack(motherLabel));
  if(!motherTrack) return kNoSignal;
  TParticle *mparticle = motherTrack->Particle();
  Int_t pid = TMath::Abs(mparticle->GetPdgCode());
  AliDebug(3, Form("PDG code: %d\n", pid));

  // identify signal according to Pdg Code 
  if((pid % 1000) / 100 == 4) return kCharm;    // charmed meson, 3rd position in pdg code == 4
  if(pid / 1000 == 4) return kCharm;            // charmed baryon, 4th position in pdg code == 4
  if((pid % 1000) / 100 == 5) return kBeauty;   // beauty meson, 3rd position in pdg code == 5
  if(pid / 1000 == 5) return kBeauty;           // beauty baryon, 4th position in pdg code == 5   
  return kNoSignal;
}

//__________________________________________
Float_t AliAnalysisTaskHFE::GetRapidity(TParticle *part)
{
      // return rapidity

       Float_t rapidity;
       if(!((part->Energy() - part->Pz())*(part->Energy() + part->Pz())>0)) rapidity=-999;
       else rapidity = 0.5*(TMath::Log((part->Energy()+part->Pz()) / (part->Energy()-part->Pz())));
       return rapidity;
}

