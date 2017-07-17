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
#include <THn.h>
#include <THnSparse.h>
#include <TList.h>
#include <TProfile.h>
#include <TSystem.h>
#include <TRandom3.h>

//--AliRoot--
#include "AliAnalysisManager.h"
#include "AliAnalysisUtils.h"
#include "AliAODEvent.h"
#include "AliAODMCHeader.h"
#include "AliMultSelection.h"
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
  AliAnalysisTaskEmcalJet("AnalysisTaskEmcalJetBtagSV", kTRUE),
  fCorrMode(kFALSE),
  fDoBkgRej(kTRUE),
  fDoQAVtx(kFALSE),
  fDoFillV0Trks(kFALSE),
  fDoDetRespMtx(kFALSE),
  fDoOnlyMtxAna(kFALSE),
  fUseTriggerData(kFALSE),
  fRecJetsBranch(),
  fGenJetsBranch(),
  fPtHardName("pthard"),
  fGenNamePattern(""),
  fJetContName(""),
  fTrkContName(""),
  fRhoTaskName(""),
  fMCJetContName(""),
  fMCTrkContName(""),
  fMCRhoTaskName(""),
  fTaggingRadius(0.4),
  fMCWeight(1.),
  fInitialized(kFALSE),
  fMCXsec(0.),
  fMCAvgTrials(0.),
  fZNApercentile(0.),
  fCurrFileName(""),
  fCheckMCCrossSection(kFALSE),
  fSkipWeightInfo(kFALSE),
  fUseWeight(kFALSE),
  fOutputList(NULL),
  fhJetVtxSim(NULL),
  fhJetVtxData(NULL),
  fhQaVtx(NULL),
  fhEntries(NULL),
  fhEvtRej(NULL),
  fhZNApercentQa(NULL),
  fhHFjetQa(NULL),
  fhRhoQa(NULL),
  fhMCRhoQa(NULL),
  fhDeltaPt(NULL),
  fhnDetRespMtx(NULL),
  fhnGenerated(NULL),
  fhXsec(NULL),
  fhTrials(NULL),
  fEvent(NULL),
  fMCHeader(NULL),
  fMultSelection(NULL),
  fTagger(NULL),
  fCutsHFjets(NULL),
  fAnalysisUtils(NULL),
  fMCTracksCont(NULL),
  fRecJetArray(NULL),
  fRecTrkArray(NULL),
  fMCJetArray(NULL),
  fMCPartArray(NULL),
  fHFvertexing(NULL),
  fV0gTrkMap(NULL),
  fRandom(new TRandom3(0)),
  fGlLogLevel(AliLog::kError),
  fLcDebLevel(1),
  fStartBin(0)
{
  // default constructor
}

//_____________________________________________________________________________________
AliAnalysisTaskEmcalJetBtagSV::AliAnalysisTaskEmcalJetBtagSV(const char* name):
  AliAnalysisTaskEmcalJet(name, kTRUE),
  fCorrMode(kFALSE),
  fDoBkgRej(kTRUE),
  fDoQAVtx(kFALSE),
  fDoFillV0Trks(kFALSE),
  fDoDetRespMtx(kFALSE),
  fDoOnlyMtxAna(kFALSE),
  fUseTriggerData(kFALSE),
  fRecJetsBranch(),
  fGenJetsBranch(),
  fPtHardName("pthard"),
  fGenNamePattern(""),
  fJetContName(""),
  fTrkContName(""),
  fRhoTaskName(""),
  fMCJetContName(""),
  fMCTrkContName(""),
  fMCRhoTaskName(""),
  fTaggingRadius(0.4),
  fMCWeight(1.),
  fInitialized(kFALSE),
  fMCXsec(0.),
  fMCAvgTrials(0.),
  fZNApercentile(0.),
  fCurrFileName(""),
  fCheckMCCrossSection(kFALSE),
  fSkipWeightInfo(kFALSE),
  fUseWeight(kFALSE),
  fOutputList(NULL),
  fhJetVtxSim(NULL),
  fhJetVtxData(NULL),
  fhQaVtx(NULL),
  fhnDetRespMtx(NULL),
  fhnGenerated(NULL),
  fhEntries(NULL),
  fhEvtRej(NULL),
  fhZNApercentQa(NULL),
  fhHFjetQa(NULL),
  fhRhoQa(NULL),
  fhMCRhoQa(NULL),  
  fhDeltaPt(NULL),
  fhXsec(NULL),
  fhTrials(NULL),
  fEvent(NULL),
  fMCHeader(NULL),
  fMultSelection(NULL),
  fTagger(NULL),
  fCutsHFjets(NULL),
  fAnalysisUtils(NULL),
  fMCTracksCont(NULL),
  fRecJetArray(NULL),
  fRecTrkArray(NULL),
  fMCJetArray(NULL),
  fMCPartArray(NULL),
  fHFvertexing(NULL),
  fV0gTrkMap(NULL),
  fRandom(new TRandom3(0)),
  fGlLogLevel(AliLog::kError),
  fLcDebLevel(1),
  fStartBin(0)
{
  // standard constructor
  AliInfo(MSGINFO("+++ Executing Constructor +++"));

  DefineOutput(1, TList::Class());
}

//_____________________________________________________________________________________
AliAnalysisTaskEmcalJetBtagSV::~AliAnalysisTaskEmcalJetBtagSV()
{
  // destructor
  AliInfo(MSGINFO("+++ Executing Destructor +++"));

  // Do not delete outputs in proof mode or merging will fail
  if (!AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
    if (fOutputList)  delete fOutputList;
    if (fHFvertexing) delete fHFvertexing;
    if (fV0gTrkMap)   delete fV0gTrkMap;
  }
  
  if (fTagger)     delete fTagger;
  if (fCutsHFjets) delete fCutsHFjets;

  if (fRandom) delete fRandom;
}

//_____________________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::UserCreateOutputObjects()
{
  AliInfo(MSGINFO("+++ Executing UserCreateOutputObjects +++"));

  // Initialize output list of containers
  if (fOutputList != NULL) {
    delete fOutputList; fOutputList = NULL;
  }

  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  // vertices within the jet - Correction Mode (MC)
  if (fCorrMode) {
    if (!fDoOnlyMtxAna) {
      fhJetVtxSim = new AliHFJetsContainerVertex("kJetVtxSim", AliHFJetsContainerVertex::kJetVtxSim);
      fOutputList->Add(fhJetVtxSim);
    }

    if (fDoDetRespMtx) {
      // detector response matrix for unfolding (from Gyulnara)
      // dimensions: pt_reco, pt_gen, eta_reco, eta_gen, flavor{g=1, L=2, C=3, B=4} BH and BP
	  const int kNbins = 6;
      Int_t bins[kNbins]    = {200, 200,  20, 20,    5,   5};
      Double_t xmin[kNbins] = {  0,   0, -1., -1., -.5, -.5};
      Double_t xmax[kNbins] = {200, 200,  1.,  1., 4.5, 4.5};
      fhnDetRespMtx = new THnSparseF("fhnDetRespMtx", "Detector response matrix", kNbins, bins, xmin, xmax);
      fOutputList->Add(fhnDetRespMtx);
      
      // MC generated histogram is needed to calculate efficiency during unfolding
      // dimensions: pt_gen, eta_gen, flavor{g=1, L=2, C=3, B=4} BH and BP
	  const Int_t kNhbins = 4;
      Int_t binsh[kNhbins]  =   {200,  20,    5,   5};
      Double_t xminh[kNhbins] = {  0,  -1., -.5, -.5};
      Double_t xmaxh[kNhbins] = {200,   1., 4.5, 4.5};
      fhnGenerated = new THnF("fhnGenerated", "MC Generated histogram", kNhbins, binsh, xminh, xmaxh);
      fOutputList->Add(fhnGenerated);
    }
  } else {
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
  fhEntries->GetXaxis()->SetBinLabel(2, "nEvPhySel");
  fhEntries->GetXaxis()->SetBinLabel(3, "nEvGoodJetArray");
  fhEntries->GetXaxis()->SetBinLabel(4, "nEvPile-up spd");
  fhEntries->GetXaxis()->SetBinLabel(5, "nEvPile-up mv");
  fhEntries->GetXaxis()->SetBinLabel(6, "nTracksEv");
  fhEntries->GetXaxis()->SetBinLabel(7, "nJetsCand");
  fhEntries->GetXaxis()->SetBinLabel(8, "nJetsTagged");
  fhEntries->GetXaxis()->SetBinLabel(9, "nUnexpError");
  fOutputList->Add(fhEntries);

  fhEvtRej  = new TH1F("fhEvtRej", "Event rejection criteria", 10, -.5, 10.5);
  fOutputList->AddLast(fhEvtRej);

  fhZNApercentQa = new TH1F("fhZNApercentQa", "ZNA multiplicity percentile;percent;dN/d(percent)", 100, 0., 100.);
  fOutputList->Add(fhZNApercentQa);

  fhHFjetQa = new TH1F("fhHFjetQa", "Check some QA for HF Jets", 12, 0.5, 12.5);
  fhHFjetQa->GetXaxis()->SetBinLabel(1,  "nWrgProngValue");      //Must be 0;
  fhHFjetQa->GetXaxis()->SetBinLabel(2,  "nJetiRejInFindVtx");   //Must be 0;
  fhHFjetQa->GetXaxis()->SetBinLabel(3,  "nJet_wTrks<2");
  fhHFjetQa->GetXaxis()->SetBinLabel(4,  "nJet_wGoodTrk<2"); fStartBin = 5;
  fhHFjetQa->GetXaxis()->SetBinLabel(5,  "nTrk_wWrongIndex");
  fhHFjetQa->GetXaxis()->SetBinLabel(6,  "nTrk_All");
  fhHFjetQa->GetXaxis()->SetBinLabel(7,  "nTrk_wFB_4");
  fhHFjetQa->GetXaxis()->SetBinLabel(8,  "nTrk_wFB_9");
  fhHFjetQa->GetXaxis()->SetBinLabel(9,  "nTrk_wNegID");
  fhHFjetQa->GetXaxis()->SetBinLabel(10, "nTrk_wFB_4&9");
  fhHFjetQa->GetXaxis()->SetBinLabel(11, "nTrk_wFB_!4!9");
  fhHFjetQa->GetXaxis()->SetBinLabel(12, "nTrkRejected");
  fOutputList->AddLast(fhHFjetQa);

  fhRhoQa = new TH1F("fhRhoQa", "Rho distribution;Probability Density;pt", 600, -10, 50);
  fOutputList->Add(fhRhoQa);
  if (fCorrMode) {
    fhMCRhoQa = new TH1F("fMCRhoQa", "Rho distribution;Probability Density;pt", 600, -10, 50);
    fOutputList->Add(fhMCRhoQa);
  }

  if (fDoRndmCone) {
	fhDeltaPt = new TH1F("fDeltaPt", "DeltaPt distribution", 500, -125, +125);
	fOutputList->Add(fhDeltaPt);
  }

  
  fhXsec = new TProfile("hXsec", "xsec from pyxsec.root", 1, 0.5, 1.5);
  fhXsec->GetXaxis()->SetBinLabel(1, Form("SelEvent_%s", fPtHardName.Data()));
  fhXsec->GetXaxis()->SetTitle("p_{T} hard bin");
  fhXsec->GetYaxis()->SetTitle("#<sigma>");
  fOutputList->Add(fhXsec);

  fhTrials = new TH1F("hTrials", "trials root file", 1, 0.5, 1.5);
  fhTrials->GetXaxis()->SetBinLabel(1, Form("SelEvent_%s", fPtHardName.Data()));
  fhTrials->GetXaxis()->SetTitle("p_{T} hard bin");
  fhTrials->GetYaxis()->SetTitle("#sum{ntrials}");
  fOutputList->Add(fhTrials);


  PostData(1, fOutputList);
}

//_____________________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::ExecOnce()
{
  AliInfo(MSGINFO("+++ Executing ExecOnce +++"));
  if (fGlLogLevel) {
    AliLog::SetGlobalLogLevel((AliLog::EType_t)fGlLogLevel);
  }

  if ((fGlLogLevel > AliLog::kInfo) && (fLcDebLevel)) {
    AliLog::SetGlobalDebugLevel(0);
    AliLog::SetClassDebugLevel("AliHFJetsTaggingVertex",        fLcDebLevel);
    AliLog::SetClassDebugLevel("AliAnalysisTaskEmcalJetBtagSV", fLcDebLevel);
  }

  fInitialized = kTRUE;
  return;
}

//_____________________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::UserExec(Option_t* /*option*/)
{
  if (!fInitialized)
    ExecOnce();

  AliDebug(1, MSGINFO("+++ Executing UserExec +++"));

  // AOD input event
  fEvent = dynamic_cast<AliAODEvent*> (InputEvent());
  if (!fEvent) {
    AliWarning(MSGWARNING("Input AOD not available, trying with output handler..."));
    if (AODEvent() && IsStandardAOD()) {
      // In case there is an AOD filter writing a standard AOD, use the AOD
      // event in memory rather than the input event.
      fEvent = dynamic_cast<AliAODEvent*>(AODEvent());
    } else {
      AliError(MSGERROR("No AOD handler found or no standard AOD!"));
      return;
    }
  }

  // load MC header
  fMCHeader = (AliAODMCHeader*)fEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
  if (!fMCHeader && fCorrMode) {
    AliError(MSGERROR("MC header branch not found!"));
    return;
  }
  


  if (!fGenNamePattern.IsNull() || !fGenNamePattern.IsWhitespace()) {
    TString title = (fMCHeader) ? fMCHeader->GetGeneratorName() : "";
    if (!title.IsNull() && !title.Contains(fGenNamePattern)) {
      AliDebugF(3, MSGWARNING("Pattern not found in MC header title %s"), title.Data());
      return;
    }
  }


  // get multiplicity and centrality percentile
  if (fEvent) 
    fMultSelection = (AliMultSelection*) fEvent->FindListObject("MultSelection");
    
  if(fMultSelection) { 
	fZNApercentile = fMultSelection->GetMultiplicityPercentile("ZNA");
    fhZNApercentQa->Fill(fZNApercentile);
  }
  else {
	fZNApercentile = -1;
	AliWarning("AliMultSelection object not found!");
  }

  // ALL EVENTS
  fhEntries->Fill(0); // EventsAnal
  if (fCheckMCCrossSection && !fSkipWeightInfo) {
    fhXsec->Fill(  1., fMCXsec);
    fhTrials->Fill(1., fMCAvgTrials);
    fSkipWeightInfo = kTRUE;
  }

  if (!fCutsHFjets->IsEventSelected((AliAODEvent*)fEvent)) {
    AliDebug(5, MSGDEBUG("Event did not pass event selection from AliRDHFJetsCuts!"));
    fhEvtRej->Fill(fCutsHFjets->GetWhyRejection(), 1);
    return;
  }

  fhEntries->Fill(1); // event selected, pileup, trigger, etc...

  if (!GetArrays()) {
    PostData(1, fOutputList);
    return;
  }
  fhEntries->Fill(2); // events with Jet arrays

  fHFvertexing = (!fHFvertexing) ? new TClonesArray("AliAODVertex", 0) : fHFvertexing;
  fHFvertexing->SetOwner(kTRUE);

  if (fDoFillV0Trks && !FillMapOfV0gTrkIDs()) {
    AliError(MSGERROR("Error filling V0 tracks info"));
    PostData(1, fOutputList);
    return;
  }

  // Execute analysis for current event
  if (fCorrMode)
    AnalyseCorrectionsMode(); // must be MC, all steps are filled for container kBJets (only)
  else
    AnalyseDataMode();        // can also be MC, only step kCFStepReco is filled also for kBJets

  PostData(1, fOutputList);
  return;
}

//_____________________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::AnalyseDataMode()
{
  // Convert to AliESDVertex // mettere in metodo separato nel task, mi servira' anche dopo TODO
  AliAODVertex* pVtx = (AliAODVertex*)fEvent->GetPrimaryVertex();

  Double_t pvXYZ[3], pvCov[6];

  pVtx->GetXYZ(pvXYZ);
  pVtx->GetCovarianceMatrix(pvCov);

  AliESDVertex* esdVtx = new AliESDVertex(pvXYZ, pvCov, pVtx->GetChi2(), pVtx->GetNContributors());

  Double_t magzkG = (Double_t)fEvent->GetMagneticField();

  Int_t nJets = fRecJetArray->GetEntries();

  Double_t rho = (fDoBkgRej) ? GetExternalRho(kFALSE) : 0.;
  fhRhoQa->Fill(rho);

  if (fDoRndmCone) {
	Double_t deltapt = GetDeltaPtRandomCone(fTaggingRadius, rho);
	fhDeltaPt->Fill(deltapt);
  }

  vctr_pair_dbl_int aVtxDisp;
  aVtxDisp.reserve(5);   // reserve space for 5 vertex sigma position

  AliEmcalJet* jet;
  for (Int_t jetcand = 0; jetcand < nJets; ++jetcand) {
    jet = (AliEmcalJet*) fRecJetArray->UncheckedAt(jetcand);

    if (!fCutsHFjets->IsJetSelected(jet)) {
      AliDebugF(5, MSGDEBUG("--> Jet with pt=%3.2f and eta=%3.2f not selected in FindVertices"),
                jet->Pt(), jet->Eta());
      continue;
    }
    Double_t ptJet_wBkgRej = jet->Pt() - (jet->Area() * rho);
    // Run b-tagger
    Int_t  nDauRejCount = 0;
    Int_t nVtx = fTagger->FindVertices(jet,
                                       fRecTrkArray,
                                       (AliAODEvent*)fEvent,
                                       esdVtx,
                                       magzkG,
                                       fHFvertexing,
                                       fV0gTrkMap,
                                       aVtxDisp,
                                       nDauRejCount);
    fhHFjetQa->Fill(12, nDauRejCount);
    if (nVtx < 0) {
      fhHFjetQa->Fill(-1 * nVtx);
      continue;
    }
    fhJetVtxData->FillStepJetVtxData(AliHFJetsContainer::kCFStepReco,
                                     nVtx,
                                     fZNApercentile,
                                     ptJet_wBkgRej,
                                     aVtxDisp,
                                     fHFvertexing,
                                     pVtx,
                                     jet,
                                     fMCWeight);

    fHFvertexing->Clear();
    aVtxDisp.clear();
  }

  delete esdVtx;
}

//_____________________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::AnalyseCorrectionsMode()
{
  // Convert to AliESDVertex // mettere in metodo separato nel task, mi servira' anche dopo TODO
  AliAODVertex* pVtx = (AliAODVertex*)fEvent->GetPrimaryVertex();

  Double_t pvXYZ[3], pvCov[6];

  pVtx->GetXYZ(pvXYZ);
  pVtx->GetCovarianceMatrix(pvCov);

  AliESDVertex* esdVtx = new AliESDVertex(pvXYZ, pvCov, pVtx->GetChi2(), pVtx->GetNContributors());

  Double_t magzkG = (Double_t)fEvent->GetMagneticField();

  // MC primary vertex
  Double_t vtxTrue[3];
  fMCHeader->GetVertex(vtxTrue);

  Int_t nMCJets = fMCJetArray->GetEntries();
  Int_t nJets   = fRecJetArray->GetEntries();

  Double_t rhoMC = (fDoBkgRej) ? GetExternalRho(kTRUE)  : 0.;
  Double_t rho   = (fDoBkgRej) ? GetExternalRho(kFALSE) : 0.;
  fhRhoQa->Fill(rho, fMCWeight);
  fhMCRhoQa->Fill(rhoMC, fMCWeight);

  if (fDoRndmCone) {
	Double_t deltapt = GetDeltaPtRandomCone(fTaggingRadius, rho);
	fhDeltaPt->Fill(deltapt, fMCWeight);
  }

  vctr_pair_dbl_int aVtxDisp;
  aVtxDisp.reserve(5);

  // Loop on MC jets
  AliEmcalJet* jetMC;
  for (Int_t jetcand = 0; jetcand < nMCJets; ++jetcand) {

    jetMC = (AliEmcalJet*)fMCJetArray->UncheckedAt(jetcand);
    if (!jetMC) continue;

    // restrict jet eta and pT ranges
    if (!fCutsHFjets->IsJetSelected(jetMC)) {
      AliDebugF(5, MSGDEBUG("JetMC not selected: pT=%f, eta=%f!"), jetMC->Pt(), jetMC->Eta());
      continue;
    }

    // Get jet flavour from 2 methods
    Double_t partonnatMC[2] = { -1, -1};
    Double_t ptpartMC[2]    = { -1, -1};

    GetFlavour2Methods(jetMC, partonnatMC, ptpartMC, fTaggingRadius);

    Double_t ptJetGen_wBkgRej = jetMC->Pt() - (jetMC->Area() * rhoMC);

    if (!fDoOnlyMtxAna) {
      // Fill container tagger
      // At this point we do not need to fill the secondary vertex QA container
       fhJetVtxSim->FillStepJetVtxSim(AliHFJetsContainer::kCFStepEventSelected, 0, 0, ptJetGen_wBkgRej, aVtxDisp,
                                     NULL, NULL, jetMC, NULL, partonnatMC, ptpartMC, fMCWeight);

       aVtxDisp.clear(); //just on case
    } 
    
    if (fDoDetRespMtx) {
      Double_t vector[4] = {ptJetGen_wBkgRej, jetMC->Eta(), partonnatMC[0], partonnatMC[1]};
      fhnGenerated->Fill(vector, fMCWeight);
	}
    
  } // end loop on jets
  // Loop on jets (clusterized on RECO particles)
  AliEmcalJet* jet;
  for (Int_t jetcand = 0; jetcand < nJets; ++jetcand) {

    jet = (AliEmcalJet*)fRecJetArray->UncheckedAt(jetcand);
    if (!fCutsHFjets->IsJetSelected(jet)) {
      AliDebugF(5, MSGDEBUG("Jet not selected: pT=%f, eta=%f!"), jet->Pt(), jet->Eta());
      continue;
    }

    // Get jet flavour from 3 methods
    Double_t partonnat[2] = { -1, -1};
    Double_t ptpart[2]    = { -1, -1};

    GetFlavour2Methods(jet, partonnat, ptpart, fTaggingRadius);

    Double_t ptJet_wBkgRej = jet->Pt() - (jet->Area() * rho);
    Int_t nVtx = 0;
    if (!fDoOnlyMtxAna) {
      // // Run vertex tagging
      Int_t nDauRejCount = 0;
      nVtx = fTagger->FindVertices(jet,
                                   fRecTrkArray,
                                   (AliAODEvent*)fEvent,
                                   esdVtx,
                                   magzkG,
                                   fHFvertexing,
                                   fV0gTrkMap,
                                   aVtxDisp,
                                   nDauRejCount);
      fhHFjetQa->Fill(12, nDauRejCount);
      if (nVtx < 0) {
        fhHFjetQa->Fill(-1 * nVtx);
        continue;
      }

      // Fill jet-with-vertex container
      fhJetVtxSim->FillStepJetVtxSim(AliHFJetsContainer::kCFStepReco,
                                     nVtx,
                                     fZNApercentile,
                                     ptJet_wBkgRej,
                                     aVtxDisp,
                                     fHFvertexing,
                                     pVtx,
                                     jet,
                                     fMCPartArray,
                                     partonnat,
                                     ptpart,
                                     fMCWeight);

      if (fDoQAVtx) {
        fhQaVtx->FillStepQaVtx(AliHFJetsContainer::kCFStepReco, nVtx, 0, pVtx, jet,
                               fHFvertexing, fMCPartArray, aVtxDisp, partonnat, fMCWeight);
      }
    }



    AliEmcalJet* matchedjet = jet->ClosestJet();
    if (!matchedjet) continue;

    GetFlavour2Methods(matchedjet, partonnat, ptpart, fTaggingRadius);

    Double_t ptJetMC_wBkgRej = matchedjet->Pt() - (matchedjet->Area() * rhoMC);

    if (!fDoOnlyMtxAna) {
      // step kCFStepMatchedAny
      fhJetVtxSim->FillStepJetVtxSim(AliHFJetsContainer::kCFStepMatchedAny,
                                     nVtx,
                                     fZNApercentile,
                                     ptJetMC_wBkgRej,
                                     aVtxDisp,
                                     fHFvertexing,
                                     pVtx,
                                     matchedjet,
                                     fMCPartArray,
                                     partonnat,
                                     ptpart,
                                     fMCWeight);

      if (fDoQAVtx) {
        fhQaVtx->FillStepQaVtx(AliHFJetsContainer::kCFStepMatchedAny, nVtx, 0, pVtx, jet,
                               fHFvertexing, fMCPartArray, aVtxDisp, partonnat, fMCWeight);
      }
    }
    fHFvertexing->Clear();
    aVtxDisp.clear();

    // Fill detector response multidimensional histo
    if (fDoDetRespMtx) {
      Double_t vector[6] = {ptJet_wBkgRej, ptJetMC_wBkgRej, jet->Eta(), matchedjet->Eta(), partonnat[0], partonnat[1]};
      fhnDetRespMtx->Fill(vector, fMCWeight);
    }

  }

  delete esdVtx;
}

//_____________________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::Terminate(const Option_t*)
{
  //TERMINATE METHOD: NOTHING TO DO
  AliInfo(MSGINFO("+++ Executing Terminate +++"));
}

//_____________________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::GetPythiaCrossSection()
{
  // Fetch the aod also from the input in,
  // have todo it in notify

  Float_t xsection  = 0;
  Float_t trials    = 1;

  TTree* tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  if (!tree) {
    AliWarning(MSGWARNING("Analysis Manager could not found the tree...."));
    return;
  }

  TFile* curfile = tree->GetCurrentFile();
  if (!curfile) {
    AliWarning(MSGWARNING("Current file not found"));
    return;
  }

  // Check if file not accessed previously, if so
  // return the previously calculated weight
  if (fCurrFileName == curfile->GetName()) {
    AliDebug(3, MSGDEBUG("File already read. Previous values used"));
    return;
  }
  fCurrFileName = TString(curfile->GetName());

  Bool_t ok = GetPythiaInfoFromFile(fCurrFileName, xsection, trials);
  if (!ok) {
    AliWarning(MSGWARNING("Parameters from file not recovered properly"));
    return;
  }

  fMCXsec = xsection;
  fMCAvgTrials = trials;

  // average number of trials
  Float_t nEntries = (Float_t)tree->GetTree()->GetEntries();

  //  if (trials >= nEntries && nEntries > 0.) fMCavgTrials /= nEntries;
  AliDebugF(3, MSGINFO("xs %e, trial %e, avg trials %2.2f, events per file %e"),
            xsection, trials, fMCAvgTrials, nEntries);

  AliDebugF(3, MSGDEBUG("Reading File %s"), curfile->GetName());
  if (fMCAvgTrials > 0.) {
    fMCWeight =  (fUseWeight) ? fMCXsec / fMCAvgTrials : 1.;
  } else {
    AliWarningF(MSGWARNING("Average number of trials is NULL!! Set weight to 1: xs : %e, trials %e, entries %e"),
                xsection, trials, nEntries);
    fMCWeight = 1.;
  }

  AliDebugF(3, MSGINFO("MC Weight: %f"), fMCWeight);
  return;
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetBtagSV::GetPythiaInfoFromFile(TString file,
                                                            Float_t& xsec,
                                                            Float_t& trials)
{
  xsec   = 0;
  trials = 1;

  if (file.Contains("root_archive.zip#")) {
    Ssiz_t pos1 = file.Index("root_archive", 12, 0, TString::kExact);
    Ssiz_t pos  = file.Index("#", 1, pos1, TString::kExact);
    Ssiz_t pos2 = file.Index(".root", 5, TString::kExact);
    file.Replace(pos + 1, pos2 - pos1, "");
  } else {
    // not an archive take the basename....
    file.ReplaceAll(gSystem->BaseName(file.Data()), "");
  }

  //Printf("%s",file.Data());

  TFile* fxsec = TFile::Open(Form("%s%s", file.Data(), "pyxsec.root"));
  if (!fxsec || fxsec->IsZombie()) {
    // next trial fetch the histgram file
    fxsec = TFile::Open(Form("%s%s", file.Data(), "pyxsec_hists.root"));
    if (!fxsec || fxsec->IsZombie()) {
      // not a severe condition but inciate that we have no information
      return kFALSE;
    } else {
      // find the tlist we want to be independtent of the name so use the Tkey
      TKey* key = (TKey*)fxsec->GetListOfKeys()->At(0);
      if (!key) {
        fxsec->Close();
        return kFALSE;
      }

      TList* list = dynamic_cast<TList*>(key->ReadObj());
      if (!list) {
        fxsec->Close();
        return kFALSE;
      }

      xsec    = ((TProfile*)list->FindObject("h1Xsec"))  ->GetBinContent(1);
      trials  = ((TH1F*)    list->FindObject("h1Trials"))->GetBinContent(1);
      fxsec->Close();
    }
  } // no tree pyxsec.root
  else {
    TTree* xtree = (TTree*)fxsec->Get("Xsection");
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
void AliAnalysisTaskEmcalJetBtagSV::GetFlavour2Methods(AliEmcalJet* jet,
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

  AliAODMCParticle* parton[2];

  parton[0] = (AliAODMCParticle*) fTagger->IsMCJetParton(fMCPartArray, jet, radius);  // method 2
  parton[1] = (AliAODMCParticle*) fTagger->IsMCJetMeson(fMCPartArray, jet, radius);   // method 3

  if (parton[0]) {
    Int_t pdg = TMath::Abs(parton[0]->PdgCode());
    //if(pdg==4 || pdg==5)
    AliDebugF(6, MSGDEBUG("parton method -> pdg parton: %d"), pdg);

    if      (pdg == 21) partonnat[0] = 1;
    else if (pdg  < 4 ) partonnat[0] = 2;
    else if (pdg == 4 ) partonnat[0] = 3;
    else if (pdg == 5 ) partonnat[0] = 4;

    ptpart[0] = parton[0]->Pt();
  } else {

    AliWarning(MSGWARNING("No parton method output"));
  }

  if (!parton[1])
    partonnat[1] = 2;
  else {
    Int_t pdg = TMath::Abs(parton[1]->PdgCode());

    AliDebugF(6, MSGDEBUG("meson method -> pdg parton: %d"), pdg);

    if      ((pdg >= 400 && pdg <= 500) || (pdg >= 4000 && pdg <= 5000)) partonnat[1] = 3;
    else if ((pdg >= 500 && pdg <= 600) || (pdg >= 5000 && pdg <= 6000)) partonnat[1] = 4;

    ptpart[1] = parton[1]->Pt();
  }
}

//_____________________________________________________________________________________
void AliAnalysisTaskEmcalJetBtagSV::CheckTrackQAinJets()
{
  if (fRecJetArray && fRecTrkArray) {
    AliEmcalJet* jet;
    for (Int_t jetcand = 0; jetcand < fRecJetArray->GetEntries(); ++jetcand) {
      jet = (AliEmcalJet*) fRecJetArray->UncheckedAt(jetcand);

      Int_t nTrksInJet = jet->GetNumberOfTracks();
      for (Int_t j = 0; j < nTrksInJet; ++j) {
        //utilize dynamic cast and then check pointer
        AliAODTrack* jTrk   = ((AliAODTrack*)jet->TrackAt(j, fRecTrkArray));
        if (!jTrk) {
          AliWarningF(MSGWARNING("Track in Jet with index %d/%d not found. Total number of AODtracks %d"),
                      j, nTrksInJet, fRecTrkArray->GetEntries());
          fhHFjetQa->Fill(fStartBin);  //fStartBin
          continue;
        }
        fhHFjetQa->Fill(fStartBin + 1);  //fStartBin + 1
        Int_t jTrkID = jTrk->GetID();
        UInt_t trkFilterMap = jTrk->GetFilterMap();
        bool hasFilterBit4 = TESTBIT(trkFilterMap, 4);
        bool hasFilterBit9 = TESTBIT(trkFilterMap, 9);

        if (hasFilterBit4)                    fhHFjetQa->Fill(fStartBin + 2); //bin fStartBin + 2
        if (hasFilterBit9)                    fhHFjetQa->Fill(fStartBin + 3); //bin fStartBin + 3
        if (jTrkID < 0)                       fhHFjetQa->Fill(fStartBin + 4); //bin fStartBin + 4
        if (hasFilterBit4 && hasFilterBit9)   fhHFjetQa->Fill(fStartBin + 5); //bin fStartBin + 5
        if (!hasFilterBit4 && !hasFilterBit9) fhHFjetQa->Fill(fStartBin + 6); //bin fStartBin + 6
      } //loop over tracks
    }   //loop over jets
  } else {
    AliError(MSGERROR("Either Jet or Tracks array found."));
  }//if arrays

  return;
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetBtagSV::GetArrays()
{
  // Get jet and track collections
  if (!fJetContName.IsNull()) {
    AliDebugF(4, MSGDEBUG("Retrieve jets %s!"), fJetContName.Data());

    fRecJetArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fJetContName));
    if (!fRecJetArray) {
      AliErrorF(MSGERROR("%s: Could not retrieve jets %s!"), GetName(), fJetContName.Data());
      return kFALSE;
    }
  }

  if (!fTrkContName.IsNull()) {
    AliDebugF(4, MSGDEBUG("Retrieve tracks %s!"), fTrkContName.Data());

    fRecTrkArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fTrkContName));
    if (!fRecTrkArray) {
      AliErrorF(MSGERROR("%s: Could not retrieve tracks %s!"), GetName(), fTrkContName.Data());
      return kFALSE;
    }
  }

  //Get MC jet and particles collections
  if (fCorrMode) {

    if (!fMCJetContName.IsNull()) {
      AliDebugF(4, MSGDEBUG("Retrieve MC jets %s!"), fMCJetContName.Data());

      fMCJetArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fMCJetContName));
      if (!fMCJetArray) {
        AliErrorF(MSGERROR("%s: Could not retrieve MC jets %s!"), GetName(), fMCJetContName.Data());
        return kFALSE;
      }
    }

    TString fMCPartContName(AliAODMCParticle::StdBranchName());
    AliDebugF(4, MSGDEBUG("Retrieve MC particles %s!"), fMCPartContName.Data());

    fMCPartArray = dynamic_cast<TClonesArray*>(fEvent->FindListObject(fMCPartContName));
    if (!fMCPartArray) {
      AliError( MSGERROR("MC particles branch not found!"));
      return kFALSE;
    }
  }

  CheckTrackQAinJets();

  return kTRUE;
}

//-------------------------------------------------------------------------------------
Bool_t AliAnalysisTaskEmcalJetBtagSV::FillMapOfV0gTrkIDs()
{

  fV0gTrkMap = (!fV0gTrkMap) ? new map_int_bool() : fV0gTrkMap;  //  Map of V0 trks (std::map< Int_t, Bool_t >              map_int_bool;)
  fV0gTrkMap->clear();

  //Fill V0 tracks AODTrack map array

  std::vector<Int_t> TrkIDs;
  FillVecOfV0gTrkIDs(TrkIDs);

  for (Int_t iTrk = 0; iTrk < fEvent->GetNumberOfTracks(); ++iTrk) {

    AliAODTrack* track = static_cast<AliAODTrack*>(fEvent->GetTrack(iTrk));

    Int_t trkID = track->GetID();

    //IMPORTANT:
    //Use same good track selection as it's implemented
    //in AliHFJetsTaggingVertex::FindVertex
    if (trkID < 0) continue;
    if (!track->GetFilterMap() && !track->GetTPCNcls()) continue;

    Bool_t trkBelongToV0 = IsAODtrkBelongToV0(TrkIDs, trkID);
    fV0gTrkMap->insert(std::pair<Int_t, Bool_t> (trkID, trkBelongToV0));
  }

  return kTRUE;
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetBtagSV::FillVecOfV0gTrkIDs(std::vector<Int_t>& vctrTrkIDs)
{

  //Fill V0 tracks AODTrack map array
  if (!fEvent)
    return kFALSE;

  AliAODEvent* aodEvent = (AliAODEvent*)fEvent;

  Int_t nV0s = aodEvent->GetNumberOfV0s();

  vctrTrkIDs.clear();
  vctrTrkIDs.reserve(nV0s);

  for (Int_t iV0 = 0; iV0 < nV0s; ++iV0) {

    AliAODv0* aodV0 = (AliAODv0*)aodEvent->GetV0(iV0);
    if (!aodV0) continue;

    AliAODTrack* pTrack = (AliAODTrack*)aodV0->GetDaughter(0);
    AliAODTrack* nTrack = (AliAODTrack*)aodV0->GetDaughter(1);

    vctrTrkIDs.push_back(pTrack->GetID());
    vctrTrkIDs.push_back(nTrack->GetID());
  }

  return kTRUE;
}

//_____________________________________________________________________________________
Bool_t AliAnalysisTaskEmcalJetBtagSV::IsAODtrkBelongToV0(std::vector<Int_t>& vctrTrkIDs, Int_t trkID)
{
  vector<Int_t>::iterator it = std::find(vctrTrkIDs.begin(), vctrTrkIDs.end(), trkID);
  if (it != vctrTrkIDs.end())
    return kTRUE;

  return kFALSE;
}

//_____________________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetBtagSV::GetExternalRho(Bool_t isMC)
{
  // Get rho from event using CMS approach
  AliRhoParameter* rho = NULL;
  TString rhoname = (!isMC) ? fRhoTaskName : fMCRhoTaskName;
  if (!rhoname.IsNull()) {
    rho = dynamic_cast<AliRhoParameter*>(InputEvent()->FindListObject(rhoname.Data()));
    if (!rho) {
      AliWarningF(MSGWARNING("%s: Could not retrieve rho with name %s!"), GetName(), rhoname.Data());
      return 0.;
    }
  } else {
    AliWarningF(MSGWARNING("No %s Rho task name provided"), (!isMC ? "DATA" : "MC"));
    return 0.;
  }

  return rho->GetVal();
}

//_____________________________________________________________________________________
Double_t AliAnalysisTaskEmcalJetBtagSV::GetDeltaPtRandomCone(Double_t jetradius, Double_t rhovalue)
{
	Double_t ptcut = 0.15;
	Double_t etacut = 0.9;
	Double_t minConeEta = jetradius - etacut;
	Double_t maxConeEta = etacut - jetradius;

	// throw random cone
	Double_t coneEta = minConeEta + fRandom->Rndm()*(maxConeEta - minConeEta);
	Double_t conePhi = fRandom->Rndm()*TMath::TwoPi();

	// collect track pt within cone
	Double_t conePt = 0.;
	for (Int_t i = 0; i < fEvent->GetNumberOfTracks(); i++) {
		AliAODTrack* trk = static_cast<AliAODTrack*>(fEvent->GetTrack(i));

		// track filter hardwired...
		UInt_t trkFilterMap = trk->GetFilterMap();  
        if (!TESTBIT(trkFilterMap, 4) && !TESTBIT(trkFilterMap, 9)) continue;

		if ( (fabs(trk->Eta()) < etacut) && (trk->Pt() > ptcut) ) {
			Double_t dphi = TVector2::Phi_mpi_pi((trk->Phi() - conePhi));
			Double_t deta = trk->Eta() - coneEta;
			Double_t dist = sqrt(deta*deta + dphi*dphi);
			if (dist < jetradius) conePt += trk->Pt();
		}
	} // track loop
		
	if (conePt > ptcut) // sanity check: at least one track found
	  return conePt - jetradius*jetradius*TMath::Pi() * rhovalue;
	else 
	  return -999;
}
