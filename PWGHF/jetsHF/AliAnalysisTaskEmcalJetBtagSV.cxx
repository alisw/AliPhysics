/*************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

/* Class AliAnalysisTaskEmcalJetBtagSV:                                    *
 * AliAnalysisTaskSE for the extraction of the B-jet spectrum        *
 * using the properties of secondary verices as tagging observables. */

/* Mailto: andrea.rossi@cern.ch, elena.bruna@to.infn.it, svallero@to.infn.it, s.lapointe@cern.ch */

#include <TH1F.h>
#include <TH2F.h>
#include <TAxis.h>
#include <TArrayI.h>
#include <TArrayD.h>
#include <THnSparse.h>
#include <TDatabasePDG.h>
#include <TMath.h>
#include <TROOT.h>
#include <TList.h>
#include <TChain.h>
#include <TFile.h>
#include <TKey.h>
#include <TProfile.h>

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliGenEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenerator.h"
#include "AliAODEvent.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoDecay.h"
#include "AliAnalysisDataSlot.h"
#include "AliAnalysisDataContainer.h"
#include "AliAODTrack.h"
#include "AliAODHandler.h"
#include "AliESDtrack.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "AliVertexerTracks.h"
#include "AliAODMCParticle.h"
#include "AliAODPid.h"
#include "AliTPCPIDResponse.h"
#include "AliAODMCHeader.h"
#include "AliAnalysisTaskEmcalJetBtagSV.h"
#include "AliRDHFCutsD0toKpi.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliNormalizationCounter.h"
#include "AliLog.h"
#include "AliHFJetsContainerVertex.h"
#include "AliAnalysisHelperJetTasks.h"

#include "AliAnalysisUtils.h"

#include "AliVCluster.h"
#include "AliAODCaloCluster.h"
#include "AliESDCaloCluster.h"
#include "AliVTrack.h"
#include "AliEmcalJet.h"
#include "AliRhoParameter.h"
#include "AliLog.h"
#include "AliJetContainer.h"
#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliPicoTrack.h"

class TCanvas;
class TTree;
class TChain;
class AliCFContainer;

ClassImp(AliAnalysisTaskEmcalJetBtagSV)

AliAnalysisTaskEmcalJetBtagSV::AliAnalysisTaskEmcalJetBtagSV() 
: AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetBtagSV", kTRUE),
fhJetVtx(0),
fhJetVtxData(0),
fhQaVtx(0),
fEvent(0x0),
fCorrMode(kTRUE),
fRecoJetsBranch(),
fMcJetsBranch(),
fSelectPtHard(kFALSE),
fUseTriggerData(kFALSE),
fPtHardMin(),
fFlavor(),
fTaggingRadius(),
fDoQAVtx(kFALSE),
fFiredClass(),
fNentries(0),
fNtriggers(0),
fTagger(0),
fCutsHFjets(0),
fAnalysisUtils(0x0),
fbJetArray(0),
fArrayMC(0),
fJetArray(0x0),
fMcJetArray(0x0),
fJetContName(""),
fTrackContName(""),
fRhoTaskName(""),
fMcJetContName(""),
fMcTrackContName(""),
fMcRhoTaskName(""),
fGTIp(0),fGTIn(0), fTrackBuffSize(19000),
fHistTrials(0),
fHistXsection(0),
fHistEvents(0),
fOutputList(0x0),fMCTracksCont(0)
{

  // default constructor

}

//________________________________________________________________________
AliAnalysisTaskEmcalJetBtagSV::AliAnalysisTaskEmcalJetBtagSV(const char *name) 
: AliAnalysisTaskEmcalJet(name,kTRUE),
fhJetVtx(0),
fhJetVtxData(0),
fhQaVtx(0),
fEvent(0x0),
fCorrMode(kTRUE),
fRecoJetsBranch(),
fMcJetsBranch(),
fSelectPtHard(kFALSE),
fUseTriggerData(kFALSE),
fPtHardMin(),
fFlavor(),
fTaggingRadius(),
fDoQAVtx(kFALSE),
fFiredClass(),
fNentries(0),
fNtriggers(0),
fTagger(0),
fCutsHFjets(0),
fAnalysisUtils(0x0),
fbJetArray(0),
fArrayMC(0),
fJetArray(0x0),
fMcJetArray(0x0),
fJetContName(""),
fTrackContName(""),
fRhoTaskName(""),
fMcJetContName(""),
fMcTrackContName(""),
fMcRhoTaskName(""),
fGTIp(0),fGTIn(0),fTrackBuffSize(19000),
fHistTrials(0),
fHistXsection(0),
fHistEvents(0),
fOutputList(0x0),fMCTracksCont(0)
{
  // standard constructor
  AliInfo(MSGINFO("+++ Executing Constructor +++"));

  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskEmcalJetBtagSV::~AliAnalysisTaskEmcalJetBtagSV(){

  // destructor
  AliInfo(MSGINFO("+++ Executing Destructor +++"));
  
  // Do not delete outputs in proof mode or merging will fail
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    if (fOutputList) delete fOutputList;
    // if (fhJets) delete fhJets;
    if (fhQaVtx) delete fhQaVtx;
    if (fhJetVtx) delete fhJetVtx;
    if (fhJetVtxData) delete fhJetVtxData;
  }

  if (fTagger) delete fTagger;
  if (fCutsHFjets) delete fCutsHFjets;

  // Array, note the [] with the delete
  if (fGTIp)
    delete[] fGTIp;
  fGTIp = 0;
  if (fGTIn)
    delete[] fGTIn;
  fGTIn = 0;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::Init()
{
  // Initialization
  AliInfo(MSGINFO("+++ Executing Init +++"));

  AliLog::SetGlobalDebugLevel(AliLog::kError);

}
//________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::UserCreateOutputObjects()
{
  // AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  AliInfo(MSGINFO("+++ Executing UserCreateOutputObjects +++"));

  // Initialize output list of containers
  if (fOutputList != NULL) {
    delete fOutputList;
    fOutputList = NULL;
  }
  
  if (!fOutputList) {
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
  }
  
  // Control histogram
  fNentries = new TH1F("nentriesChFr", "Analyzed sample properties", 9, -.5, 8.5);
  fNentries->GetXaxis()->SetBinLabel(1, "nEventsAnal");
  fNentries->GetXaxis()->SetBinLabel(2, "nEvSel");
  fNentries->GetXaxis()->SetBinLabel(3, "nEvGoodVtx");
  fNentries->GetXaxis()->SetBinLabel(4, "nEvPile-up spd");
  fNentries->GetXaxis()->SetBinLabel(5, "nEvPile-up mv");
  fNentries->GetXaxis()->SetBinLabel(6, "nTracksEv");
  fNentries->GetXaxis()->SetBinLabel(7, "nJetsCand");
  fNentries->GetXaxis()->SetBinLabel(8, "nJetsTagged");
  fNentries->GetXaxis()->SetBinLabel(9, "nUnexpError");
  fOutputList->Add(fNentries);

  // Control histogram
  fNtriggers = new TH1F("ntriggers", "Analyzed sample properties", 5, -.5, 5.5);
  fNtriggers->GetXaxis()->SetBinLabel(1, "nJ1");
  fNtriggers->GetXaxis()->SetBinLabel(3, "nJ2");
  fNtriggers->GetXaxis()->SetBinLabel(5, "notJ1orJ2");
  fOutputList->Add(fNtriggers);

  fHistTrials = new TH1F("fHistTrials", "fHistTrials", 11, 0, 11);
  fHistTrials->GetXaxis()->SetTitle("p_{T} hard bin");
  fHistTrials->GetYaxis()->SetTitle("trials");
  fOutputList->Add(fHistTrials);

  fHistEvents = new TH1F("fHistEvents", "fHistEvents", 11, 0, 11);
  fHistEvents->GetXaxis()->SetTitle("p_{T} hard bin");
  fHistEvents->GetYaxis()->SetTitle("total events");
  fOutputList->Add(fHistEvents);

  fHistXsection = new TProfile("fHistXsection", "fHistXsection", 11, 0, 11);
  fHistXsection->GetXaxis()->SetTitle("p_{T} hard bin");
  fHistXsection->GetYaxis()->SetTitle("xsection");
  fOutputList->Add(fHistXsection);

  
  // const Int_t ptHardLo[11] = { 0, 5,11,21,36,57, 84,117,152,191,234};
  // const Int_t ptHardHi[11] = { 5,11,21,36,57,84,117,152,191,234,1000000};
    
  // for (Int_t i = 1; i < 12; i++) {
  //   fHistTrials->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));
  //   fHistXsection->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));
  //   fHistEvents->GetXaxis()->SetBinLabel(i, Form("%d-%d",ptHardLo[i-1],ptHardHi[i-1]));
  // }

  // Create the containers for jet and vertex properties
  // all jets
  // fhJets= new AliHFJetsContainerVertex("kJets",AliHFJetsContainerVertex::kJets);
  // fOutputList->Add(fhJets);
  // vertices within the jet
  fhJetVtx = new AliHFJetsContainerVertex("kJetVtx", AliHFJetsContainerVertex::kJetVtx);
  fOutputList->Add(fhJetVtx);
  // vertices within the jet - data
  fhJetVtxData = new AliHFJetsContainerVertex("kJetVtxData", AliHFJetsContainerVertex::kJetVtxData);
  fOutputList->Add(fhJetVtxData);
  // vertices QA
  fhQaVtx = new AliHFJetsContainerVertex("kQaVtx", AliHFJetsContainerVertex::kQaVtx);
  fOutputList->Add(fhQaVtx);

  fGTIp = new AliAODTrack *[fTrackBuffSize]; // Array of pointers
  fGTIn = new AliAODTrack *[fTrackBuffSize]; // Array of pointers
  
  PostData(1, fOutputList);

}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::UserExec(Option_t */*option*/)
{

  AliInfo(MSGINFO("+++ Executing UserExec +++"));

  // Execute analysis for current event
  if (fCorrMode) AnalyseCorrectionsMode(); // must be MC, all steps are filled for container kBJets (only)
  else AnalyseDataMode(); // can also be MC, only step kCFStepReco is filled also for kBJets
  
  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::AnalyseDataMode(){

  if (!fbJetArray) fbJetArray = new TClonesArray("AliAODVertex", 0);
  Double_t arrDispersion[10000];

  // AOD input event
  AliAODEvent *aod = dynamic_cast<AliAODEvent *> (InputEvent());
  if (!aod) {
    AliInfo(MSGINFO("Input AOD not available, trying with output handler..."));
    if ( AODEvent() && IsStandardAOD() ) {
      // In case there is an AOD filter writing a standard AOD, use the AOD
      // event in memory rather than the input event.
      aod = dynamic_cast<AliAODEvent *> (AODEvent());
    }
    else AliError(MSGERROR("No AOD handler found or no standard AOD!"));
  }
  
  fEvent = aod;
  
  if ( !GetArrays() ) return;
  
  // ALL EVENTS
  fNentries->Fill(0); // EventsAnal
  
  AliAODVertex *vtx1;
  vtx1 = (AliAODVertex *) aod->GetPrimaryVertex();

  if (!fCutsHFjets->IsEventSelected(aod)) {
    AliDebug(AliLog::kDebug, MSGDEBUG("Event did not pass event selection from AliRDHFJetsCuts!"));
    return;
  }

  fNentries->Fill(1); // event selected, pileup, trigger, etc...

  Int_t mult = 0;

  // Init steps
  AliHFJetsContainer::CFSteps step;
 
  // Convert to AliESDVertex // mettere in metodo separato nel task, mi servira' anche dopo TODO 
  Double_t primvtx[3], primcov[6];
  vtx1->GetXYZ(primvtx);
  vtx1->GetCovarianceMatrix(primcov);
  Int_t nPrimContr=vtx1->GetNContributors();
  Double_t chi2=vtx1->GetChi2();
  AliESDVertex *v1 = new AliESDVertex(primvtx, primcov, chi2, nPrimContr);
  Double_t magzkG = (Double_t) aod->GetMagneticField();


  // Reset the reference array to the global tracks..
  ResetTrackReference();
  // ..and set it
  for (Int_t iTrack=0; iTrack<aod->GetNumberOfTracks(); iTrack++) {
    AliAODTrack *track = static_cast<AliAODTrack *> (aod->GetTrack(iTrack));
    if (!track) continue;
    
    // Store the reference of the global tracks
    StoreTrackReference(track);
  }

  Int_t nJets = fJetArray->GetEntries();

  AliEmcalJet *jet;
  Int_t nvtx=0;
 
  for(Int_t jetcand=0;jetcand<nJets;jetcand++){
    nvtx=0;
    jet = (AliEmcalJet *) fJetArray->UncheckedAt(jetcand);
    TClonesArray *fTrackArray = dynamic_cast<TClonesArray *> (fEvent->FindListObject("PicoTracks"));
    if (!fTrackArray) return;
    //    Int_t numtracks = fTrackArray->GetEntries();
    if ( !fCutsHFjets->IsJetSelected(jet) ) {
      AliDebug(AliLog::kDebug,Form(MSGDEBUG("Jet not selected: pT=%f, eta=%f!"), jet->Pt(), jet->Eta()));
      continue;
    }
    
    Double_t ptJet = jet->Pt();
    if (fDoBkgRej) {
      Double_t rho, areaJet;
      rho  = GetExternalRho();
      //loop over jets:
      
      areaJet = jet->Area();
      ptJet -= areaJet * rho;
    }
    step = AliHFJetsContainer::kCFStepReco;
    // // Run b-tagger
    nvtx = fTagger->FindVertices(jet, fGTIp, fGTIn, fTrackArray, aod, v1, magzkG, fbJetArray, arrDispersion);
 
    fhJetVtxData->FillStepJetVtxData(step, mult, jet, fbJetArray, nvtx, vtx1, arrDispersion, ptJet);
    fbJetArray->Clear();
  }
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::AnalyseCorrectionsMode()
{
  if ( !fbJetArray ) fbJetArray = new TClonesArray("AliAODVertex", 0);
  Double_t arrDispersion[10000];

  // AOD input event
  AliAODEvent *aod = dynamic_cast<AliAODEvent *> ( InputEvent() );
  if ( !aod ) {
    AliInfo(MSGINFO("Input AOD not available, trying with output handler..."));
    if ( AODEvent() && IsStandardAOD() ) {
      // In case there is an AOD filter writing a standard AOD, use the AOD
      // event in memory rather than the input event.
      aod = dynamic_cast<AliAODEvent *> ( AODEvent() );
    }
    else AliError(MSGERROR("No AOD handler found or no standard AOD!"));
  }

  fEvent = aod;

  if ( !GetArrays() ) return;
 
  // ALL EVENTS
  fNentries->Fill(0); // EventsAnal 

  AliAODVertex *vtx1;
  vtx1 = (AliAODVertex *) aod->GetPrimaryVertex();
  
  if ( !fCutsHFjets->IsEventSelected(aod) ) {
    AliDebug(AliLog::kDebug, MSGDEBUG("Event did not pass event selection from AliRDHFJetsCuts!"));
    return;
  }

  fNentries->Fill(1); // event selected, pileup, trigger, etc...

  // fill ntrials and xsec, at the momemt for 1 pt hard bin
  AliAODMCHeader *aodmcHeader = 0x0;

  // load MC header
  aodmcHeader = (AliAODMCHeader *)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  if ( !aodmcHeader ) AliError(MSGERROR("MC header branch not found!"));
  
  if ( fSelectPtHard ) {
    TString title = aodmcHeader->GetGeneratorName();
    if (!title.Contains(fPtHardMin) || !title.Contains(fFlavor)) return;
  }
  
  fNentries->Fill(2);    //is the pt-hard event
  
  // Multiplicity MC (for vertex reco correction)
  // Get array of MC particles
  fArrayMC = (TClonesArray *)aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if ( !fArrayMC ) AliError(MSGERROR("MC particles branch not found!"));
  
  // Count number of primary MC partcles
  Int_t nMC = fArrayMC->GetEntries();
  Int_t multMC = 0;
  for (Int_t i = 0; i < nMC; ++i) {
    AliAODMCParticle *part = (AliAODMCParticle *) fArrayMC->At(i);
    if ( part->IsPhysicalPrimary() ) multMC++;
  }
  
  // Init steps
  AliHFJetsContainer::CFSteps step;
  // Loop on MC jets
  Int_t nMCJets = fMcJetArray->GetEntries();

  Double_t rhoMC=0, rho=0;
  if (fDoBkgRej) {
    rhoMC  = GetMcExternalRho();
    rho  = GetExternalRho();
  }
      
  AliEmcalJet *jetMC;
  for (Int_t jetcand = 0; jetcand < nMCJets; jetcand++) {
    jetMC = (AliEmcalJet *)fMcJetArray->UncheckedAt(jetcand);
    if ( !jetMC ) continue;

    // restrict jet eta and pT ranges
    if ( !fCutsHFjets->IsJetSelected(jetMC) ) {
      AliDebug(AliLog::kDebug,Form(MSGDEBUG("JetMC not selected: pT=%f, eta=%f!"), jetMC->Pt(),jetMC->Eta()));
      continue;
    }
    
    TClonesArray *fTrackArray = dynamic_cast<TClonesArray *> (fEvent->FindListObject("MCParticlesSelected"));
    
    // Get jet flavour from 2 methods
    Double_t partonnatMC[2] = {-1, -1};
    Double_t ptpartMC[2] = {-1, -1};
  
    GetFlavour2Methods(jetMC, partonnatMC, ptpartMC, fTaggingRadius);
    
    // Fill container tagger
    step = AliHFJetsContainer::kCFStepEventSelected;
    
    // At this point we do not need to fill the secondary vertex QA container 
    
    Double_t areaJetMC;
    areaJetMC = jetMC->Area();
    Double_t ptJetMC = jetMC->Pt() - areaJetMC * rhoMC;
  
    fhJetVtx->FillStepJetVtx(step, multMC, jetMC, 0, 0, 0, 0, partonnatMC, ptpartMC, arrDispersion, ptJetMC);

  } // end loop on jets
  
  // Convert to AliESDVertex 
  Double_t primvtx[3],primcov[6];
  vtx1->GetXYZ(primvtx);
  vtx1->GetCovarianceMatrix(primcov);
  Int_t nPrimContr = vtx1->GetNContributors();
  Double_t chi2 = vtx1->GetChi2();
  AliESDVertex *v1 = new AliESDVertex(primvtx, primcov, chi2, nPrimContr);
  Double_t magzkG = (Double_t) aod->GetMagneticField();

  // Reset the reference array to the global tracks..
  ResetTrackReference();
  // ..and set it
    
  for (Int_t iTrack = 0; iTrack<aod->GetNumberOfTracks(); iTrack++) {
    AliAODTrack *track = static_cast<AliAODTrack *> (aod->GetTrack(iTrack));
    if (!track) continue;

    // Store the reference of the global tracks
    StoreTrackReference(track);
  }

  TClonesArray *arrayMC = 0x0;
  Double_t vtxTrue[3];

  // load MC particles
  arrayMC = (TClonesArray *) aod->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if (!arrayMC) AliError(MSGERROR("MC particles branch not found!"));

  // MC primary vertex
  aodmcHeader->GetVertex(vtxTrue);

  // Loop on jets (clusterized on RECO particles)
   
  Int_t nJets = fJetArray->GetEntries();
  AliEmcalJet *jet;

  Int_t nvtx = 0;
 
  for (Int_t jetcand = 0; jetcand<nJets; jetcand++) {
    nvtx = 0;
    jet = (AliEmcalJet *) fJetArray->UncheckedAt(jetcand);

    if (!fCutsHFjets->IsJetSelected(jet)) {
      AliDebug(AliLog::kDebug, Form(MSGDEBUG("Jet not selected: pT=%f, eta=%f!"), jet->Pt(),jet->Eta()));
      continue;
    }
    
    TClonesArray *fTrackArrayMc = dynamic_cast<TClonesArray *> (fEvent->FindListObject("MCParticlesSelected"));
    TClonesArray *fTrackArrayRec = dynamic_cast<TClonesArray *> (fEvent->FindListObject("PicoTracks"));
        
    // Get jet flavour from 3 methods
    Double_t partonnat[2] = {-1, -1};
    Double_t ptpart[2] = {-1, -1};
    GetFlavour2Methods(jet, partonnat, ptpart, fTaggingRadius);
    
    Double_t ptJet = jet->Pt();
    Double_t areaJet;
    areaJet = jet->Area();
    ptJet -= areaJet * rho;
    
    
    step = AliHFJetsContainer::kCFStepReco;
    // // Run vertex tagging
    nvtx = fTagger->FindVertices(jet, fGTIp, fGTIn, fTrackArrayRec, aod, v1, magzkG, fbJetArray, arrDispersion);
   
    if ( fDoQAVtx ) {
      fhQaVtx->FillStepQaVtx(step, multMC, jet, fbJetArray, arrDispersion, nvtx, vtx1, arrayMC, partonnat, ptJet);
    }
    // Fill jet-with-vertex container
    fhJetVtx->FillStepJetVtx(step, multMC, jet, fbJetArray, nvtx, vtx1, arrayMC, partonnat, ptpart, arrDispersion, ptJet);
    
    AliEmcalJet *matchedjet = NULL;
    
    matchedjet = jet->ClosestJet();
    if(!matchedjet) continue;
    
    GetFlavour2Methods(matchedjet, partonnat, ptpart, fTaggingRadius);

    Double_t ptJetMC = matchedjet->Pt();
    Double_t areaJetMC;
    areaJetMC = matchedjet->Area();
    ptJetMC -= areaJetMC * rhoMC;
      
    // step kCFStepMatchedAny
    step = AliHFJetsContainer::kCFStepMatchedAny;
    
    if (fDoQAVtx) {
      fhQaVtx->FillStepQaVtx(step, multMC, jet, fbJetArray, arrDispersion, nvtx, vtx1, arrayMC, partonnat, ptJetMC);
    }
    fhJetVtx->FillStepJetVtx(step, multMC, matchedjet, fbJetArray, nvtx, vtx1, arrayMC, partonnat, ptpart, arrDispersion, ptJetMC);
    
    fbJetArray->Clear();
  }
  delete v1;
  
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::Terminate(const Option_t *)
{
  //TERMINATE METHOD: NOTHING TO DO
  AliInfo(MSGINFO("+++ Executing Terminate +++"));

}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::GetFlavour2Methods(AliEmcalJet *jet, Double_t (&partonnat)[2], Double_t (&ptpart)[2], Double_t radius){

  // 2 methods to associate jet to mother parton

  /* Nature of the parton (methods 1)      *
   * 1 = gluon (pdg 21)                    *
   * 2 = light quark (pdg < 4)             *
   * 3 = c (pdg 4)                         *
   * 4 = b (pdg 5)                         *
   * Nature of the meson/barion (method 2) *
   * 2 = light                             *
   * 3 = with c                            *
   * 4 = with b                            */

  // Initialize output values
  for (Int_t i=0; i<2; i++){
    partonnat[i]=0;
    ptpart[i] = -1.;
  }
  AliAODMCParticle *parton[2];
  
  parton[0] = (AliAODMCParticle *) fTagger->IsMCJetParton(fArrayMC, jet, radius); // method 2
  
  parton[1] = (AliAODMCParticle *) fTagger->IsMCJetMeson(fArrayMC,jet,radius); // method 3

  Printf(" %p  %p ", parton[0], parton[1]);
  
  if( parton[0] ) {
    Int_t pdg = TMath::Abs(parton[0]->PdgCode());
    //if(pdg==4 || pdg==5)
    AliInfo(Form(MSGINFO("parton method -> pdg parton: %d"), pdg));
    if (pdg == 21)     partonnat[0] = 1;
    else if (pdg < 4)  partonnat[0] = 2;
    else if (pdg == 4) partonnat[0] = 3;
    else if (pdg == 5) partonnat[0] = 4;
    ptpart[0] = parton[0]->Pt();
  }

  if (parton[1]) {
    Int_t pdg = TMath::Abs(parton[1]->PdgCode());
    AliInfo(Form(MSGINFO("meson method -> pdg parton: %d"), pdg));
    if ( (pdg >= 400 && pdg <= 500) || (pdg >= 4000 && pdg <= 5000)) partonnat[1]=3;
    else {
      if ( (pdg >= 500 && pdg <= 600) || (pdg >= 5000 && pdg <= 6000)) partonnat[1]=4;
    }
    ptpart[1]=parton[1]->Pt();
  }
  else partonnat[1]=2;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetBtagSV::GetArrays()
{
  // Get jet collection
  if (!fJetContName.IsNull()) {
    AliInfo(Form(MSGINFO("Retrieve jets %s!"), fJetContName.Data()));
    fJetArray = dynamic_cast<TClonesArray *> (fEvent->FindListObject(fJetContName));
    if (!fJetArray) {
      AliError(Form(MSGERROR("%s: Could not retrieve jets %s!"), GetName(), fJetContName.Data()));
      return kFALSE;
    }
  }
  // Get Mc jet collection
  
  if (!fMcJetContName.IsNull() && fCorrMode) {
    AliInfo(Form(MSGINFO("Retrieve jets %s!"), fMcJetContName.Data()));
    fMcJetArray = dynamic_cast<TClonesArray *> (fEvent->FindListObject(fMcJetContName));
    if (!fMcJetArray) {
      AliError(Form(MSGERROR("%s: Could not retrieve MC jets %s!"), GetName(), fMcJetContName.Data()));
      return kFALSE;
    }
  }
  
  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::StoreTrackReference(AliAODTrack *track)
{
  // Stores the pointer to the global (positive) track
  
  // Check that the id is positive
  if (track->GetID()> -1) {
    //    printf("Warning: track has negative ID: %d\n",track->GetID());
    // return;
    // Check id is not too big for buffer
    if (track->GetID() >= fTrackBuffSize) {
      Printf(MSGWARNING("Warning: track ID too big for buffer: ID: %d, buffer %d\n")
             ,track->GetID(), fTrackBuffSize);
      return;
    }

    // Warn if we overwrite a track
    if (fGTIp[track->GetID()]) {
      // Seems like there are FilterMap 0 tracks
      // that have zero TPCNcls, don't store these!
      if( !track->GetFilterMap() && !track->GetTPCNcls())
        return;
      
      // Imagine the other way around,
      // the zero map zero clusters track
      // is stored and the good one wants 
      // to be added. We ommit the warning
      // and just overwrite the 'bad' track
      if ( fGTIp[track->GetID()]->GetFilterMap() || fGTIp[track->GetID()]->GetTPCNcls()) {
        // If we come here, there's a problem
        AliWarning(MSGWARNING("Warning! global track info already there!"));
        Printf(MSGINFO("         TPCNcls track1 %u track2 %u"),
               (fGTIp[track->GetID()])->GetTPCNcls(),
               track->GetTPCNcls()
               );
        Printf(MSGINFO("         FilterMap track1 %u track2 %u\n"),
               (fGTIp[track->GetID()])->GetFilterMap(),
               track->GetFilterMap()
               );
      }
    } // Two tracks same id
    
    // Assign the pointer
    (fGTIp[track->GetID()]) = track;
  }
  

  // Check that the id is negative
  if (track->GetID() < 0) {
    Int_t idx = TMath::Abs(track->GetID());
    // printf("Warning: track has negative ID: %d\n",track->GetID());
    // return;
    // Check id is not too big for buffer
    if ( idx >= fTrackBuffSize) {
      Printf(MSGWARNING("Warning: track ID too big for buffer: ID: %d, buffer %d\n"),
             track->GetID(),
             fTrackBuffSize
             );
      return;
    }
    
    // Warn if we overwrite a track
    if (fGTIn[idx]) {
      // Seems like there are FilterMap 0 tracks
      // that have zero TPCNcls, don't store these!
      if( !track->GetFilterMap() && !track->GetTPCNcls() )
        return;

      // Imagine the other way around,
      // the zero map zero clusters track
      // is stored and the good one wants 
      // to be added. We omit the warning
      // and just overwrite the 'bad' track
      if ( fGTIn[idx]->GetFilterMap() || fGTIn[idx]->GetTPCNcls()) {
        // If we come here, there's a problem
        Printf(MSGWARNING("Warning! global track info already there!"));
        Printf(MSGINFO("         TPCNcls track1 %u track2 %u"),
               fGTIn[idx]->GetTPCNcls(),
               track->GetTPCNcls()
               );
        Printf(MSGINFO("         FilterMap track1 %u track2 %u\n"),
               fGTIn[idx]->GetFilterMap(),
               track->GetFilterMap()
               );
      }
    } // Two tracks same id
    
    // Assign the pointer
    fGTIn[idx] = track;
  }
}

//________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::ResetTrackReference()
{
  // Sets all the pointers to zero. To be called at
  // the beginning or end of an event
  for (UShort_t i = 0; i<fTrackBuffSize; i++) {
    fGTIp[i]=0;
    fGTIn[i]=0;
  }
  
}

//------------------------------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskEmcalJetBtagSV::PythiaInfoFromFile(const char* currFile, Float_t &fXsec, Float_t &fTrials)
{
  //
  // Get the cross section and the trails either from pyxsec.root or from pysec_hists.root
  // Get the pt hard bin from the file path
  // This is to called in Notify and should provide the path to the AOD/ESD file
  // (Partially copied from AliAnalysisHelperJetTasks)

  TString file(currFile);  
  fXsec = 0;
  fTrials = 1;

  if (file.Contains(".zip#")) {
    Ssiz_t pos1 = file.Index("root_archive",12,0,TString::kExact);
    Ssiz_t pos = file.Index("#",1,pos1,TString::kExact);
    Ssiz_t pos2 = file.Index(".root",5,TString::kExact);
    file.Replace(pos+1,pos2-pos1,"");
  } else {
    // not an archive take the basename....
    file.ReplaceAll(gSystem->BaseName(file.Data()),"");
  }
  AliDebug(1, Form(MSGDEBUG("File name: %s"), file.Data()));

  // problem that we cannot really test the existance of a file in a archive so we have to live with open error message from root
  TFile *fxsec = TFile::Open(Form("%s%s",file.Data(), "pyxsec.root"));
  
  if (!fxsec) {
    // next trial fetch the histgram file
    fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec_hists.root"));
    if (!fxsec) {
      // not a severe condition but inciate that we have no information
      return kFALSE;
    }
    else {
      // find the tlist we want to be independtent of the name so use the Tkey
      TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0);
      if (!key) {
        fxsec->Close();
        return kFALSE;
      }
      
      TList *list = dynamic_cast<TList *> (key->ReadObj());
      if (!list) {
        fxsec->Close();
        return kFALSE;
      }
      fXsec = ((TProfile *)list->FindObject("h1Xsec"))->GetBinContent(1);
      fTrials  = ((TH1F*)list->FindObject("h1Trials"))->GetBinContent(1);
      fxsec->Close();
    }
  }
  else { // no tree pyxsec.root
    TTree *xtree = (TTree*)fxsec->Get("Xsection");
    if (!xtree) {
      fxsec->Close();
      return kFALSE;
    }
    UInt_t   ntrials  = 0;
    Double_t  xsection  = 0;
    xtree->SetBranchAddress("xsection",&xsection);
    xtree->SetBranchAddress("ntrials",&ntrials);
    xtree->GetEntry(0);
    fTrials = ntrials;
    fXsec = xsection;
    fxsec->Close();
  }
  
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetBtagSV::UserNotify()
{
  // Called when file changes.

  if (!fCorrMode)
    return kTRUE;

  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  if (!tree) {
    AliError(Form(MSGERROR("%s - UserNotify: No current tree!"), GetName()));
    return kFALSE;
  }
  AliAODEvent *aod = dynamic_cast<AliAODEvent*> (InputEvent());

  AliAODMCHeader *aodmcHeader=0x0;
  // load MC header
  aodmcHeader =
    (AliAODMCHeader*)aod->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  if (!aodmcHeader) AliError(MSGERROR("MC header branch not found!"));

  TString title = aodmcHeader->GetGeneratorName();
  if ( fSelectPtHard ) {
     if (!title.Contains(fPtHardMin) || !title.Contains(fFlavor)) return kFALSE; 
  }
   
  Float_t xsection    = 0;
  Float_t trials      = 0;
  Int_t   pthardbin   = 0;

  // get pthard from the title, where bin 1 corresponds to pt-hard min of 10, bin 2 to pt-hard min 11(lf), 18(hf-b), 20(hf-c) and so on
  pthardbin = fPtHardMin.Atoi()/10.;

  TFile *curfile = tree->GetCurrentFile();
  if (!curfile) {
    AliError(Form(MSGERROR("%s - UserNotify: No current file!"), GetName()));
    return kFALSE;
  }

  TChain *chain = dynamic_cast<TChain*>(tree);
  if (chain) tree = chain->GetTree();

  Int_t nevents = tree->GetEntriesFast();

  PythiaInfoFromFile(curfile->GetName(), xsection, trials);

  fHistTrials->Fill(pthardbin, trials);
  fHistXsection->Fill(pthardbin, xsection);
  fHistEvents->Fill(pthardbin, nevents);

  return kTRUE;
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetBtagSV::GetExternalRho()
{
  // Get rho from event using CMS approach
  AliRhoParameter *rho = NULL;
  if (!fRhoTaskName.IsNull()) {
    rho = dynamic_cast<AliRhoParameter *> (InputEvent()->FindListObject(fRhoTaskName.Data()));
    if (!rho) {
      AliWarning(Form(MSGWARNING("%s: Could not retrieve rho with name %s!"), GetName(), fRhoTaskName.Data()));
      return 0.;
    }
  }
  else {
    AliWarning(MSGWARNING("No Rho task name provided"));
    return 0.;
  }
  
  return rho->GetVal();
}

//________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetBtagSV::GetMcExternalRho()
{
  // Get rho from event using CMS approach
  AliRhoParameter *rho = NULL;
  if (!fMcRhoTaskName.IsNull()) {
    rho = dynamic_cast<AliRhoParameter *> (InputEvent()->FindListObject(fMcRhoTaskName.Data()));
    if (!rho) {
      AliWarning(Form(MSGWARNING("%s: Could not retrieve rho with name %s!"), GetName(), fMcRhoTaskName.Data()));
      return 0.;
    }
  }
  else {
    AliWarning(MSGWARNING("No Rho task name provided"));
    return 0.;
  }
  
  return rho->GetVal();
}
