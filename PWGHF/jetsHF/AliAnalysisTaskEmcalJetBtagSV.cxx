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

/* Class AliAnalysisTaskEmcalJetBtagSV:                               *
 * AliAnalysisTaskSE for the extraction of the B-jet spectrum        *
 * using the properties of secondary verices as tagging observables. */

/*
 Mailto: andrea.rossi@cern.ch, elena.bruna@to.infn.it,
 svallero@to.infn.it, s.lapointe@cern.ch
 ycorrale@cern.ch
 */


//--Root--
#include <TFile.h>
#include <TH1F.h>
#include <TKey.h>
#include <TList.h>
#include <TProfile.h>
#include <TSystem.h>

//--AliRoot--
#include "AliAnalysisManager.h"
#include "AliAnalysisUtils.h"
#include "AliAODEvent.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "AliLog.h"
#include "AliRhoParameter.h"

//--AliHFJetsClass--
#include "AliHFJetsTaggingVertex.h"
#include "AliHFJetsContainerVertex.h"

#include "AliAnalysisTaskEmcalJetBtagSV.h"

//_____________________________________________________________________________________

ClassImp(AliAnalysisTaskEmcalJetBtagSV)

//_____________________________________________________________________________________
AliAnalysisTaskEmcalJetBtagSV::AliAnalysisTaskEmcalJetBtagSV() :
AliAnalysisTaskEmcalJet("AliAnalysisTaskEmcalJetBtagSV", kTRUE),
fCorrMode(kFALSE),
fDoBkgRej(kTRUE),
fDoQAVtx(kFALSE),
fDoFillV0Trks(kFALSE),
fUseTriggerData(kFALSE),
fRecJetsBranch(),
fGenJetsBranch(),
fGenNamePattern(""),
fJetContName(""),
fTrkContName(""),
fRhoTaskName(""),
fMcJetContName(""),
fMcTrkContName(""),
fMcRhoTaskName(""),
fFiredClass(""),
fTaggingRadius(0.4),
fMCWeight(1.),
fMCXsec(0.),
fMCavgTrials(0.),
fCurrFileName(""),
fCheckMCCrossSection(kFALSE),
fIsMCInfoFilled(kFALSE),
fDebug(AliLog::kInfo),
fOutputList(NULL),
fhJetVtxSim(NULL),
fhJetVtxData(NULL),
fhQaVtx(NULL),
fhEntries(NULL),
fhXsec(NULL),
fhTrials(NULL),
fEvent(NULL),
fTagger(NULL),
fCutsHFjets(NULL),
fAnalysisUtils(NULL),
fMCTracksCont(NULL),
fbJetArray(NULL),
fArrayMC(NULL),
fJetArray(NULL),
fMcJetArray(NULL)
{
  // default constructor
}

//_____________________________________________________________________________________
AliAnalysisTaskEmcalJetBtagSV::AliAnalysisTaskEmcalJetBtagSV(const char *name):
AliAnalysisTaskEmcalJet(name, kTRUE),
fCorrMode(kFALSE),
fDoBkgRej(kTRUE),
fDoQAVtx(kFALSE),
fDoFillV0Trks(kFALSE),
fUseTriggerData(kFALSE),
fRecJetsBranch(),
fGenJetsBranch(),
fGenNamePattern(""),
fJetContName(""),
fTrkContName(""),
fRhoTaskName(""),
fMcJetContName(""),
fMcTrkContName(""),
fMcRhoTaskName(""),
fFiredClass(""),
fTaggingRadius(0.4),
fMCWeight(1.),
fMCXsec(0.),
fMCavgTrials(0.),
fCurrFileName(""),
fCheckMCCrossSection(kFALSE),
fIsMCInfoFilled(kFALSE),
fDebug(AliLog::kInfo),
fOutputList(NULL),
fhJetVtxSim(NULL),
fhJetVtxData(NULL),
fhQaVtx(NULL),
fhEntries(NULL),
fhXsec(NULL),
fhTrials(NULL),
fEvent(NULL),
fTagger(NULL),
fCutsHFjets(NULL),
fAnalysisUtils(NULL),
fMCTracksCont(NULL),
fbJetArray(NULL),
fArrayMC(NULL),
fJetArray(NULL),
fMcJetArray(NULL)
{
  // standard constructor
  AliInfo(MSGINFO("+++ Executing Constructor +++"));

  DefineOutput(1, TList::Class());
}

//_____________________________________________________________________________________
AliAnalysisTaskEmcalJetBtagSV::~AliAnalysisTaskEmcalJetBtagSV() {

  // destructor
  AliInfo(MSGINFO("+++ Executing Destructor +++"));
  
  // Do not delete outputs in proof mode or merging will fail
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {

    if (fOutputList)  delete fOutputList;
    if (fhJetVtxSim)  delete fhJetVtxSim;
    if (fhJetVtxData) delete fhJetVtxData;
    if (fhQaVtx)      delete fhQaVtx;
  }
  
  if (fTagger)     delete fTagger;
  if (fCutsHFjets) delete fCutsHFjets;
}

//_____________________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::Init() {
  
  // Initialization
  AliInfo(MSGINFO("+++ Executing Init +++"));
}

//_____________________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::UserCreateOutputObjects() {
  
  // AliAnalysisTaskEmcalJet::UserCreateOutputObjects();

  AliInfo(MSGINFO("+++ Executing UserCreateOutputObjects +++"));

  // Initialize output list of containers
  if (fOutputList != NULL) {
    delete fOutputList; fOutputList = NULL;
  }
  
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  // vertices within the jet - Correction Mode (MC)
  if (fCorrMode) {
    fhJetVtxSim = new AliHFJetsContainerVertex("kJetVtxSim", AliHFJetsContainerVertex::kJetVtxSim);
    fOutputList->Add(fhJetVtxSim);
  }
  else {
    // vertices within the jet - Data
    fhJetVtxData = new AliHFJetsContainerVertex("kJetVtxData", AliHFJetsContainerVertex::kJetVtxData);
    fOutputList->Add(fhJetVtxData);
  }
  // vertices QA
  if (fDoQAVtx) {
    // vertices within the jet - QA
    fhQaVtx = new AliHFJetsContainerVertex("kQaVtx", AliHFJetsContainerVertex::kQaVtx);
    fOutputList->Add(fhQaVtx);
  }

  // Control histogram
  fhEntries = new TH1F("hEntries", "Analyzed sample properties", 9, -.5, 8.5);
  fhEntries->GetXaxis()->SetBinLabel(1, "nEventsAnal");
  fhEntries->GetXaxis()->SetBinLabel(2, "nEvSel");
  fhEntries->GetXaxis()->SetBinLabel(3, "nEvGoodVtx");
  fhEntries->GetXaxis()->SetBinLabel(4, "nEvPile-up spd");
  fhEntries->GetXaxis()->SetBinLabel(5, "nEvPile-up mv");
  fhEntries->GetXaxis()->SetBinLabel(6, "nTracksEv");
  fhEntries->GetXaxis()->SetBinLabel(7, "nJetsCand");
  fhEntries->GetXaxis()->SetBinLabel(8, "nJetsTagged");
  fhEntries->GetXaxis()->SetBinLabel(9, "nUnexpError");
  fOutputList->Add(fhEntries);

  fhXsec = new TH1F("hXsec", "xsec from pyxsec.root", 2, -0.5, 1.5);
  fhXsec->GetXaxis()->SetBinLabel(1, "AllFiles");
  fhXsec->GetXaxis()->SetBinLabel(2, Form("SelEvent_%s", fGenNamePattern.Data()));
  fhXsec->GetXaxis()->SetTitle("p_{T} hard bin");
  fhXsec->GetYaxis()->SetTitle("#<sigma>");
  fOutputList->Add(fhXsec);

  fhTrials = new TH1F("hTrials", "trials root file", 2, -0.5, 1.5);
  fhTrials->GetXaxis()->SetBinLabel(1, "AllFiles");
  fhTrials->GetXaxis()->SetBinLabel(2, Form("SelEvent_%s", fGenNamePattern.Data()));
  fhTrials->GetXaxis()->SetTitle("p_{T} hard bin");
  fhTrials->GetYaxis()->SetTitle("#sum{ntrials}");
  fOutputList->Add(fhTrials);
  
  PostData(1, fOutputList);
}

//_____________________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::UserExec(Option_t */*option*/) {
  
  AliLog::SetGlobalLogLevel(fDebug);
  
  AliInfo(MSGINFO("+++ Executing UserExec +++"));

  // Execute analysis for current event
  if (fCorrMode)
    AnalyseCorrectionsMode(); // must be MC, all steps are filled for container kBJets (only)
  else
    AnalyseDataMode();        // can also be MC, only step kCFStepReco is filled also for kBJets
  
  PostData(1, fOutputList);
}

//_____________________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::AnalyseDataMode()
{
  // AOD input event
  fEvent = dynamic_cast<AliAODEvent *> (InputEvent());
  if (!fEvent) {
    
    AliInfo(MSGINFO("Input AOD not available, trying with output handler..."));
    if (AODEvent() && IsStandardAOD()) {
      // In case there is an AOD filter writing a standard AOD, use the AOD
      // event in memory rather than the input event.
      fEvent = dynamic_cast<AliAODEvent *>(AODEvent());
    }
    else
      AliError(MSGERROR("No AOD handler found or no standard AOD!"));
  }
  
  // ALL EVENTS
  fhEntries->Fill(0); // EventsAnal

  if (!fCutsHFjets->IsEventSelected((AliAODEvent *)fEvent)) {
    
    AliDebug(AliLog::kDebug, MSGDEBUG("Event did not pass event selection from AliRDHFJetsCuts!"));
    return;
  }

  fhEntries->Fill(1); // event selected, pileup, trigger, etc...

  if (!GetArrays())
    return;
  
  if (!fbJetArray)
    fbJetArray = new TClonesArray("AliAODVertex", 0);
  
  //Fill V0 tracks AODTrack map array
  
  std::vector<Int_t> TrkIDs;
  if (fDoFillV0Trks) FillV0trks(TrkIDs);

  map_AliAODTrk               fAODgTrkMap;       //  Map of AOD trks (std::map<Int_t, AliAODTrack * >)
  fAODgTrkMap.clear();
  
  for (Int_t iTrk = 0; iTrk < fEvent->GetNumberOfTracks(); ++iTrk) {
    
    AliAODTrack *track = static_cast<AliAODTrack *>(fEvent->GetTrack(iTrk));
    
    Int_t trkID = track->GetID();
    if (trkID < 0) continue;
    
    Bool_t trkBelongToV0 = (fDoFillV0Trks) ? IsV0(trkID, TrkIDs) : kFALSE;
    fAODgTrkMap[trkID] = std::make_pair(track, trkBelongToV0);
  }
  
  //End AODTrack Map
  
  // Convert to AliESDVertex // mettere in metodo separato nel task, mi servira' anche dopo TODO
  AliAODVertex *pVtx = (AliAODVertex *)fEvent->GetPrimaryVertex();
  
  Double_t pvXYZ[3], pvCov[6];
  
  pVtx->GetXYZ(pvXYZ);
  pVtx->GetCovarianceMatrix(pvCov);
  
  AliESDVertex *esdVtx = new AliESDVertex(pvXYZ, pvCov, pVtx->GetChi2(), pVtx->GetNContributors());
  
  Double_t magzkG = (Double_t)fEvent->GetMagneticField();
  
  vctr_pair_dbl_int aVtxDisp;
  aVtxDisp.reserve(10);
  
  Int_t nJets = fJetArray->GetEntries();
  
  Double_t rho = (fDoBkgRej) ? GetExternalRho() : 0.;
  
  AliEmcalJet *jet;
  for (Int_t jetcand = 0; jetcand < nJets; ++jetcand) {
    
    jet = (AliEmcalJet *) fJetArray->UncheckedAt(jetcand);
    
    TClonesArray *fTrackArray = dynamic_cast<TClonesArray *>(fEvent->FindListObject("PicoTracks"));
    if (!fTrackArray)
      return;
    
    if (!fCutsHFjets->IsJetSelected(jet)) {
      
      AliDebugF(AliLog::kDebug, MSGDEBUG("Jet not selected: pT=%f, eta=%f!"), jet->Pt(), jet->Eta());
      continue;
    }
    
    Double_t ptJet_wBkgRej = jet->Pt() - (jet->Area() * rho);
    
    aVtxDisp.clear();
    // Run b-tagger
    Int_t nVtx = fTagger->FindVertices(jet, fAODgTrkMap, fTrackArray, (AliAODEvent *)fEvent,
                                       esdVtx, magzkG, fbJetArray, aVtxDisp);
    
    fhJetVtxData->FillStepJetVtxData(AliHFJetsContainer::kCFStepReco, nVtx, 0, ptJet_wBkgRej, aVtxDisp, fbJetArray, pVtx, jet, fMCWeight);
    
    fbJetArray->Clear();
  }
  
  delete esdVtx;
}

//_____________________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::AnalyseCorrectionsMode() {
  
  // AOD input event
  fEvent = dynamic_cast<AliAODEvent *> (InputEvent());
  if (!fEvent) {
    
    AliInfo(MSGINFO("Input AOD not available, trying with output handler..."));
    if (AODEvent() && IsStandardAOD()) {
      // In case there is an AOD filter writing a standard AOD, use the AOD
      // event in memory rather than the input event.
      fEvent = dynamic_cast<AliAODEvent *>(AODEvent());
    }
    else
      AliError(MSGERROR("No AOD handler found or no standard AOD!"));
  }
  
  // load MC header
  AliAODMCHeader *aodmcHeader = (AliAODMCHeader *)fEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  if (!aodmcHeader) {
    
    AliError(MSGERROR("MC header branch not found!"));
    return;
  }
  
  if (fCheckMCCrossSection && (!fGenNamePattern.IsNull() || !fGenNamePattern.IsWhitespace())) {
    
    TString title = aodmcHeader->GetGeneratorName();
    if (!title.Contains(fGenNamePattern))
      return;
    
    if (!fIsMCInfoFilled) {
      
      fhXsec->Fill(  1., fMCXsec);
      fhTrials->Fill(1., fMCavgTrials);
      fIsMCInfoFilled = kTRUE;
    }
  }
  
  // ALL EVENTS
  fhEntries->Fill(0); // EventsAnal
  
  if (!fCutsHFjets->IsEventSelected((AliAODEvent *)fEvent)) {
    
    AliDebug(AliLog::kDebug, MSGDEBUG("Event did not pass event selection from AliRDHFJetsCuts!"));
    return;
  }
  
  fhEntries->Fill(1); // event selected, pileup, trigger, etc...
  
  if (!GetArrays())
    return;
  
  if (!fbJetArray)
    fbJetArray = new TClonesArray("AliAODVertex", 0);

  // Get array of MC particles
  fArrayMC = (TClonesArray *)fEvent->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if (!fArrayMC)
    AliError( MSGERROR("MC particles branch not found!") );
  
  // Multiplicity MC (for vertex reco correction)
  // Count number of primary MC partcles
  //  Int_t multMC = 0;
  //  if (0) {
  //    Int_t nMC = fArrayMC->GetEntries();
  //    for (Int_t i = 0; i < nMC; ++i) {
  //
  //      AliAODMCParticle *part = (AliAODMCParticle *)fArrayMC->At(i);
  //      if (part->IsPhysicalPrimary())
  //        multMC++;
  //    }
  //  }
  vctr_pair_dbl_int aVtxDisp;
  aVtxDisp.reserve(10);
  
  // Loop on MC jets
  Int_t nMCJets = fMcJetArray->GetEntries();
  
  Double_t rhoMC = (fDoBkgRej) ? GetExternalRho(kTRUE)  : 0.;
  Double_t rho   = (fDoBkgRej) ? GetExternalRho(kFALSE) : 0.;
  
  AliEmcalJet *jetMC;
  for (Int_t jetcand = 0; jetcand < nMCJets; ++jetcand) {
    
    jetMC = (AliEmcalJet *)fMcJetArray->UncheckedAt(jetcand);
    if (!jetMC)
      continue;

    // restrict jet eta and pT ranges
    if (!fCutsHFjets->IsJetSelected(jetMC)) {
      
      AliDebugF(AliLog::kDebug, MSGDEBUG("JetMC not selected: pT=%f, eta=%f!"), jetMC->Pt(), jetMC->Eta());
      continue;
    }
    
    // Get jet flavour from 2 methods
    Double_t partonnatMC[2] = {-1, -1};
    Double_t ptpartMC[2]    = {-1, -1};

    GetFlavour2Methods(jetMC, partonnatMC, ptpartMC, fTaggingRadius);
    
    Double_t ptJetGen_wBkgRej = jetMC->Pt() - (jetMC->Area() * rhoMC);

    // Fill container tagger
    // At this point we do not need to fill the secondary vertex QA container
    fhJetVtxSim->FillStepJetVtxSim(AliHFJetsContainer::kCFStepEventSelected, 0, 0, ptJetGen_wBkgRej, aVtxDisp, NULL, NULL, jetMC, NULL, partonnatMC, ptpartMC, fMCWeight);
  } // end loop on jets

  //Fill V0 tracks AODTrack map array
  
  vector <Int_t> TrkIDs;
  if (fDoFillV0Trks) FillV0trks(TrkIDs);
  
  map_AliAODTrk               fAODgTrkMap;       //  Map of AOD trks (std::map<Int_t, AliAODTrack * >)
  fAODgTrkMap.clear();
  
  for (Int_t iTrk = 0; iTrk < fEvent->GetNumberOfTracks(); ++iTrk) {
    
    AliAODTrack *track = static_cast<AliAODTrack *>(fEvent->GetTrack(iTrk));
    Int_t trkID = track->GetID();
    if (trkID < 0) continue;
    
    Bool_t trkBelongToV0 = (fDoFillV0Trks) ? IsV0(trkID, TrkIDs) : kFALSE;
    fAODgTrkMap[trkID] = make_pair(track, trkBelongToV0);
  }

  //End AODTrack Map
  
  // Convert to AliESDVertex // mettere in metodo separato nel task, mi servira' anche dopo TODO
  AliAODVertex *pVtx = (AliAODVertex *)fEvent->GetPrimaryVertex();
  
  Double_t pvXYZ[3], pvCov[6];
  
  pVtx->GetXYZ(pvXYZ);
  pVtx->GetCovarianceMatrix(pvCov);
  
  AliESDVertex *esdVtx = new AliESDVertex(pvXYZ, pvCov, pVtx->GetChi2(), pVtx->GetNContributors());
  
  Double_t magzkG = (Double_t)fEvent->GetMagneticField();
  
  TClonesArray *arrayMC = 0x0;
  Double_t vtxTrue[3];

  // load MC particles
  arrayMC = (TClonesArray *)fEvent->GetList()->FindObject(AliAODMCParticle::StdBranchName());
  if (!arrayMC)
    AliError(MSGERROR("MC particles branch not found!"));

  // MC primary vertex
  aodmcHeader->GetVertex(vtxTrue);

  // Loop on jets (clusterized on RECO particles)
   
  Int_t nJets = fJetArray->GetEntries();
  
  AliEmcalJet *jet;
  for (Int_t jetcand = 0; jetcand < nJets; ++jetcand) {
    
    jet = (AliEmcalJet *)fJetArray->UncheckedAt(jetcand);

    if (!fCutsHFjets->IsJetSelected(jet)) {
      
      AliDebugF(AliLog::kDebug, MSGDEBUG("Jet not selected: pT=%f, eta=%f!"), jet->Pt(), jet->Eta());
      continue;
    }
    
    TClonesArray *fTrackArrayRec = dynamic_cast<TClonesArray *>(fEvent->FindListObject("PicoTracks"));
    
    // Get jet flavour from 3 methods
    Double_t partonnat[2] = {-1, -1};
    Double_t ptpart[2]    = {-1, -1};
    
    GetFlavour2Methods(jet, partonnat, ptpart, fTaggingRadius);
    
    Double_t ptJet_wBkgRej = jet->Pt() - (jet->Area() * rho);
    
    // // Run vertex tagging
    Int_t nVtx = fTagger->FindVertices(jet, fAODgTrkMap, fTrackArrayRec, (AliAODEvent *)fEvent,
                                       esdVtx, magzkG, fbJetArray, aVtxDisp);

    if (fDoQAVtx) {
      fhQaVtx->FillStepQaVtx(AliHFJetsContainer::kCFStepReco, nVtx, 0, pVtx, jet, fbJetArray, arrayMC, aVtxDisp, partonnat,fMCWeight);
    }
    // Fill jet-with-vertex container
    fhJetVtxSim->FillStepJetVtxSim(AliHFJetsContainer::kCFStepReco, nVtx, 0, ptJet_wBkgRej, aVtxDisp, fbJetArray, pVtx, jet, arrayMC, partonnat, ptpart, fMCWeight);
    
    AliEmcalJet *matchedjet = NULL;
    
    matchedjet = jet->ClosestJet();
    if (!matchedjet)
      continue;

    GetFlavour2Methods(matchedjet, partonnat, ptpart, fTaggingRadius);

    Double_t ptJetMC_wBkgRej = matchedjet->Pt() - (matchedjet->Area() * rhoMC);
      
    // step kCFStepMatchedAny
    if (fDoQAVtx) {
      fhQaVtx->FillStepQaVtx(AliHFJetsContainer::kCFStepMatchedAny, nVtx, 0, pVtx, jet, fbJetArray, arrayMC, aVtxDisp, partonnat, fMCWeight);
    }
    fhJetVtxSim->FillStepJetVtxSim(AliHFJetsContainer::kCFStepMatchedAny, nVtx, 0, ptJetMC_wBkgRej, aVtxDisp, fbJetArray, pVtx, matchedjet, arrayMC, partonnat, ptpart, fMCWeight);
    
    fbJetArray->Clear();
  }
  
  delete esdVtx;
}

//_____________________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::Terminate(const Option_t *) {
  
  //TERMINATE METHOD: NOTHING TO DO
  AliInfo(MSGINFO("+++ Executing Terminate +++"));
}

//_____________________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::GetPythiaCrossSection() {
  
  // Fetch the aod also from the input in,
  // have todo it in notify
  
  Float_t xsection  = 0;
  Float_t trials    = 1;
  Float_t avgTrials = 0;
  
  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  if (!tree) return;
  
  TFile *curfile = tree->GetCurrentFile();
  if (!curfile) return;
  
  // Check if file not accessed previously, if so
  // return the previously calculated weight
  if (fCurrFileName == curfile->GetName()) return;
  
  fCurrFileName = TString(curfile->GetName());
  
  if (!fhXsec || !fhTrials) {
    AliWarningF(MSGWARNING("%s%d No Histogram fhXsec"),(char*)__FILE__,__LINE__);
    return;
  }
  
  Bool_t ok = GetPythiaInfoFromFile(fCurrFileName,xsection,trials);
  
  if (!ok) {
    
    AliWarning(MSGWARNING("Parameters from file not recovered properly"));
    return;
  }
  
  fMCXsec = xsection;
  fhXsec->Fill(0., xsection);
  
  // average number of trials
  Float_t nEntries = (Float_t)tree->GetTree()->GetEntries();
  
  if (trials >= nEntries && nEntries > 0.) avgTrials = trials/nEntries;
  
  fMCavgTrials = avgTrials;
  fhTrials->Fill(0., avgTrials);
  
  Printf(MSGINFO("xs %e, trial %e, avg trials %2.2f, events per file %e"),
           xsection, trials, avgTrials, nEntries);
  
  AliDebugF(1, MSGDEBUG("Reading File %s"), curfile->GetName());
  
  if (avgTrials > 0.) {
    
    fMCWeight =  xsection / avgTrials ;
    
    Printf(MSGINFO("MC Weight: %e"), fMCWeight);
  }
  else {
    
    AliWarningF(MSGWARNING("Average number of trials is NULL!! Set weight to 1: xs : %e, trials %e, entries %e"),
                xsection,trials,nEntries);
    
    fMCWeight = 1.;
  }
  
  return;
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetBtagSV::GetPythiaInfoFromFile(TString file,
                                                           Float_t &xsec,
                                                           Float_t &trials)
{
  xsec   = 0;
  trials = 1;
  
  if (file.Contains("root_archive.zip#")) {
    
    Ssiz_t pos1 = file.Index("root_archive",12,0,TString::kExact);
    Ssiz_t pos  = file.Index("#",1,pos1,TString::kExact);
    Ssiz_t pos2 = file.Index(".root",5,TString::kExact);
    file.Replace(pos+1,pos2-pos1,"");
  }
  else {
    // not an archive take the basename....
    file.ReplaceAll(gSystem->BaseName(file.Data()), "");
  }
  
  //Printf("%s",file.Data());
  
  TFile *fxsec = TFile::Open(Form("%s%s",file.Data(), "pyxsec.root"));
  if (!fxsec || fxsec->IsZombie()) {
    // next trial fetch the histgram file
    fxsec = TFile::Open(Form("%s%s",file.Data(),"pyxsec_hists.root"));
    if (!fxsec || fxsec->IsZombie()) {
      // not a severe condition but inciate that we have no information
      return kFALSE;
    }
    else {
      // find the tlist we want to be independtent of the name so use the Tkey
      TKey *key = (TKey *)fxsec->GetListOfKeys()->At(0);
      if(!key) {
        fxsec->Close();
        return kFALSE;
      }
      
      TList *list = dynamic_cast<TList *>(key->ReadObj());
      if (!list) {
        fxsec->Close();
        return kFALSE;
      }
      
      xsec    = ((TProfile *)list->FindObject("h1Xsec"))  ->GetBinContent(1);
      trials  = ((TH1F *)    list->FindObject("h1Trials"))->GetBinContent(1);
      fxsec->Close();
    }
  } // no tree pyxsec.root
  else {
    TTree *xtree = (TTree *)fxsec->Get("Xsection");
    if (!xtree) {
      fxsec->Close();
      return kFALSE;
    }
    
    UInt_t   ntrials  = 0;
    Double_t xsection = 0;
    xtree->SetBranchAddress("xsection", &xsection);
    xtree->SetBranchAddress("ntrials",  &ntrials);
    xtree->GetEntry(0);
    trials = ntrials;
    xsec = xsection;
    fxsec->Close();
  }
  
  return kTRUE;
}

//_____________________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::GetFlavour2Methods(AliEmcalJet *jet,
                                                      Double_t (&partonnat)[2],
                                                      Double_t (&ptpart)[2],
                                                      Double_t radius)
{
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
  for (Int_t i = 0; i < 2; ++i) {
    partonnat[i] =   0;
    ptpart[i]    = -1.;
  }
  
  AliAODMCParticle *parton[2];
  
  parton[0] = (AliAODMCParticle *) fTagger->IsMCJetParton(fArrayMC, jet, radius); // method 2
  parton[1] = (AliAODMCParticle *) fTagger->IsMCJetMeson(fArrayMC, jet, radius);  // method 3
  
  if (parton[0]) {
    
    Int_t pdg = TMath::Abs(parton[0]->PdgCode());
    //if(pdg==4 || pdg==5)
    AliInfoF(MSGINFO("parton method -> pdg parton: %d"), pdg);
    
    if      (pdg == 21) partonnat[0] = 1;
    else if (pdg  < 4 ) partonnat[0] = 2;
    else if (pdg == 4 ) partonnat[0] = 3;
    else if (pdg == 5 ) partonnat[0] = 4;
    
    ptpart[0] = parton[0]->Pt();
  }
  else {

    AliWarning(MSGWARNING("No parton method output"));
  }

  if (!parton[1])
    partonnat[1] = 2;
  else {
    Int_t pdg = TMath::Abs(parton[1]->PdgCode());
    
    AliInfoF(MSGINFO("meson method -> pdg parton: %d"), pdg);
    
    if      ((pdg >= 400 && pdg <= 500) || (pdg >= 4000 && pdg <= 5000)) partonnat[1] = 3;
    else if ((pdg >= 500 && pdg <= 600) || (pdg >= 5000 && pdg <= 6000)) partonnat[1] = 4;
    
    ptpart[1]=parton[1]->Pt();
  }
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetBtagSV::GetArrays() {
  
  // Get jet collection
  if (!fJetContName.IsNull()) {
    
    AliInfoF(MSGINFO("Retrieve jets %s!"), fJetContName.Data());
    
    fJetArray = dynamic_cast<TClonesArray *>(fEvent->FindListObject(fJetContName));
    if (!fJetArray) {
      
      AliErrorF(MSGERROR("%s: Could not retrieve jets %s!"), GetName(), fJetContName.Data());
      return kFALSE;
    }
  }
  
  // Get Mc jet collection
  if (!fMcJetContName.IsNull() && fCorrMode) {
    
    AliInfoF(MSGINFO("Retrieve jets %s!"), fMcJetContName.Data());
    
    fMcJetArray = dynamic_cast<TClonesArray *> (fEvent->FindListObject(fMcJetContName));
    if (!fMcJetArray) {
      
      AliErrorF(MSGERROR("%s: Could not retrieve MC jets %s!"), GetName(), fMcJetContName.Data());
      return kFALSE;
    }
  }
  
  return kTRUE;
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetBtagSV::FillV0trks(std::vector<Int_t> &TrkIDs) {
  
  //Fill V0 tracks AODTrack map array
  if (!fEvent)
    return kFALSE;
  
  AliAODEvent *aod = (AliAODEvent *)fEvent;
  
  Int_t nV0s = aod->GetNumberOfV0s();
  
  TrkIDs.clear();
  TrkIDs.reserve(nV0s);
  
  for (Int_t iV0 = 0; iV0 < nV0s; ++iV0) {
    
    AliAODv0 *aodV0 = (AliAODv0 *)aod->GetV0(iV0);
    if (!aodV0) continue;
    
    AliAODTrack *pTrack = (AliAODTrack *)aodV0->GetDaughter(0);
    AliAODTrack *nTrack = (AliAODTrack *)aodV0->GetDaughter(1);
    
    TrkIDs.push_back(pTrack->GetID());
    TrkIDs.push_back(nTrack->GetID());
  }

  return kTRUE;
}
//_____________________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetBtagSV::IsV0(Int_t trkId, std::vector<Int_t> &vTrks) {
  
  vector<Int_t>::iterator it = std::find(vTrks.begin(), vTrks.end(), trkId);
  if (it != vTrks.end())
    return kTRUE;
  
  return kFALSE;
}

//_____________________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetBtagSV::GetExternalRho(Bool_t isMC) {

  // Get rho from event using CMS approach
  AliRhoParameter *rho = NULL;
  
  TString rhoname = (!isMC) ? fRhoTaskName : fMcRhoTaskName;
  if (!rhoname.IsNull()) {
    
    rho = dynamic_cast<AliRhoParameter *>(InputEvent()->FindListObject(rhoname.Data()));

    if (!rho) {
      
      AliWarningF(MSGWARNING("%s: Could not retrieve rho with name %s!"), GetName(), rhoname.Data());
      return 0.;
    }
  }
  else {
    
    AliWarning(MSGWARNING("No Rho task name provided"));
    return 0.;
  }
  
  return rho->GetVal();
}
