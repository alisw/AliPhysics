/*
 * AliOtonOmegaAnalysis.cxx
 *
 *  Created on: 24 Nov 2017
 *      Author: bernhardhohlweger. Adapted for Proton-Omega by Oton Vazquez Doce.
 */
#include <vector>
#include "AliLog.h"
#include "AliOtonOmegaAnalysis.h"
#include "TClonesArray.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include <iostream>
ClassImp(AliOtonOmegaAnalysis)
AliOtonOmegaAnalysis::AliOtonOmegaAnalysis()
    : fMVPileUp(false),
      fEvtCutQA(false),
      fQA(),
      fFemtoTrack(nullptr),
      fFemtov0(nullptr),
      fFemtoCasc(nullptr),
      fEvent(nullptr),
      fEvtCuts(nullptr),
      fTrackCuts(nullptr),
      fAntiTrackCuts(nullptr),
      fv0Cuts(nullptr),
      fAntiv0Cuts(nullptr),
      fCascCuts(nullptr),
      fAntiCascCuts(nullptr),
      fCascOmegaCuts(nullptr),
      fAntiCascOmegaCuts(nullptr),
      fPairCleaner(nullptr),
      fControlSample(nullptr),
      fTrackBufferSize(0),
      fGTI(nullptr),
      fConfig(nullptr),
      fPartColl(nullptr),
      fPIDResponse(NULL),
      fIsMC(false) {

}

AliOtonOmegaAnalysis::AliOtonOmegaAnalysis(
    const AliOtonOmegaAnalysis& analysis)
    : fMVPileUp(analysis.fMVPileUp),
      fEvtCutQA(analysis.fEvtCutQA),
      fQA(nullptr),
      fFemtoTrack(nullptr),
      fFemtov0(nullptr),
      fFemtoCasc(nullptr),
      fEvent(nullptr),
      fEvtCuts(analysis.fEvtCuts),
      fTrackCuts(analysis.fTrackCuts),
      fAntiTrackCuts(analysis.fAntiTrackCuts),
      fv0Cuts(analysis.fv0Cuts),
      fAntiv0Cuts(analysis.fAntiv0Cuts),
      fCascCuts(analysis.fCascCuts),
      fAntiCascCuts(analysis.fAntiCascCuts),
      fCascOmegaCuts(analysis.fCascOmegaCuts),
      fAntiCascOmegaCuts(analysis.fAntiCascOmegaCuts),
      fPairCleaner(nullptr),
      fControlSample(nullptr),
      fTrackBufferSize(analysis.fTrackBufferSize),
      fGTI(nullptr),
      fConfig(analysis.fConfig),
      fPartColl(nullptr),
      fPIDResponse(NULL),
      fIsMC(analysis.fIsMC) {
  Init(fIsMC, AliVEvent::kINT7);
}

AliOtonOmegaAnalysis& AliOtonOmegaAnalysis::operator=(
    const AliOtonOmegaAnalysis& analysis) {
  if (this != &analysis) {
    this->fMVPileUp = analysis.fMVPileUp;
    this->fEvtCutQA = analysis.fEvtCutQA;
    this->fQA = nullptr;
    this->fFemtoTrack = nullptr;
    this->fFemtov0 = nullptr;
    this->fFemtoCasc = nullptr;
    this->fEvent = nullptr;
    this->fEvtCuts = analysis.fEvtCuts;
    this->fTrackCuts = analysis.fTrackCuts;
    this->fAntiTrackCuts = analysis.fAntiTrackCuts;
    this->fv0Cuts = analysis.fv0Cuts;
    this->fAntiv0Cuts = analysis.fAntiv0Cuts;
    this->fCascCuts = analysis.fCascCuts;
    this->fAntiCascCuts = analysis.fAntiCascCuts;
    this->fCascOmegaCuts = analysis.fCascOmegaCuts;
    this->fAntiCascOmegaCuts = analysis.fAntiCascOmegaCuts;
    this->fPairCleaner = nullptr;
    this->fControlSample = nullptr;
    this->fTrackBufferSize = analysis.fTrackBufferSize;
    this->fGTI = nullptr;
    this->fConfig = analysis.fConfig;
    this->fPartColl = nullptr;
    this->fIsMC = analysis.fIsMC;
    this->Init(this->fIsMC, AliVEvent::kINT7);
  }
  return *this;
}

AliOtonOmegaAnalysis::~AliOtonOmegaAnalysis() {
  if (fEvent) {
    delete fEvent;
  }
  if (fFemtoTrack) {
    delete fFemtoTrack;
  }
  if (fFemtov0) {
    delete fFemtov0;
  }
  if (fFemtoCasc) {
    delete fFemtoCasc;
  }
  if (fPairCleaner) {
    delete fPairCleaner;
  }
  if (fPartColl) {
    delete fPartColl;
  }
  if(fPIDResponse) delete fPIDResponse;
  if (fControlSample) {
    delete fControlSample;
  }
}

void AliOtonOmegaAnalysis::Init(bool isMonteCarlo, UInt_t trigger) {
  fIsMC = isMonteCarlo;
  fFemtoTrack = new AliFemtoDreamTrack();
  fFemtoTrack->SetUseMCInfo(isMonteCarlo);

  fFemtov0 = new AliFemtoDreamv0();
//do not set v0s
/*
  fFemtov0->SetPDGCode(fv0Cuts->GetPDGv0());
  fFemtov0->SetUseMCInfo(isMonteCarlo);
  fFemtov0->SetPDGDaughterPos(fv0Cuts->GetPDGPosDaug());  //order +sign doesnt play a role
  fFemtov0->GetPosDaughter()->SetUseMCInfo(isMonteCarlo);
  fFemtov0->SetPDGDaughterNeg(fv0Cuts->GetPDGNegDaug());  //only used for MC Matching
  fFemtov0->GetNegDaughter()->SetUseMCInfo(isMonteCarlo);
*/


//For MC this will have to be probably changed, in order to have Xi and Omega (ask B. Hohlweger):
  fFemtoCasc = new AliOtonOmegaCascade();
  fFemtoCasc->SetUseMCInfo(isMonteCarlo);
  //PDG Codes should be set assuming Xi- to also work for Xi+
  fFemtoCasc->SetPDGCode(fCascCuts->GetPDGCodeCasc());
  fFemtoCasc->SetPDGDaugPos(fCascCuts->GetPDGCodePosDaug());
  fFemtoCasc->GetPosDaug()->SetUseMCInfo(isMonteCarlo);
  fFemtoCasc->SetPDGDaugNeg(fCascCuts->GetPDGCodeNegDaug());
  fFemtoCasc->GetNegDaug()->SetUseMCInfo(isMonteCarlo);
  fFemtoCasc->SetPDGDaugBach(fCascCuts->GetPDGCodeBach());
  fFemtoCasc->GetBach()->SetUseMCInfo(isMonteCarlo);
  fFemtoCasc->Setv0PDGCode(fCascCuts->GetPDGv0());


  fEvtCuts->InitQA();
  fTrackCuts->Init();
  fAntiTrackCuts->Init();
//  fv0Cuts->Init();
//  fAntiv0Cuts->Init();
  fCascCuts->Init();
  fAntiCascCuts->Init();
 if(!fCascOmegaCuts) { AliFatal("omegaCuts not set");}
  fCascOmegaCuts->Init();
  fAntiCascOmegaCuts->Init();
  fGTI = new AliAODTrack*[fTrackBufferSize];
  fEvent = new AliFemtoDreamEvent(fMVPileUp, fEvtCutQA, trigger);
  fEvent->SetMultiplicityEstimator(fConfig->GetMultiplicityEstimator());
  bool MinBooking = !((!fConfig->GetMinimalBookingME())
      || (!fConfig->GetMinimalBookingSample()));
  fPairCleaner = new AliFemtoDreamPairCleaner(4, 4, MinBooking);

  if (!MinBooking) {
    fQA = new TList();
    fQA->SetOwner();
    fQA->SetName("QA");
    fQA->Add(fPairCleaner->GetHistList());
    if (fEvtCutQA) {
      fQA->Add(fEvent->GetEvtCutList());
    }
  }
  if (fConfig->GetUseEventMixing()) {
    fPartColl = new AliFemtoDreamPartCollection(fConfig,
                                                fConfig->GetMinimalBookingME());
  }
  if (fConfig->GetUsePhiSpinning()) {
    fControlSample = new AliFemtoDreamControlSample(
        fConfig);
  }

  InitializeTreeBooking();

  return;
}


void AliOtonOmegaAnalysis::InitializeTreeBooking() {

 AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
 AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
 fPIDResponse = inputHandler->GetPIDResponse();

 //book tree:
 foTTree = new TTree("oTTree","a simple TTree");
 foTTree->SetDirectory(0); // This is to force a memory-resident Tree, and avoid errors. // ????? is this necessary? does it create memory problems?
 foTTree->Branch("RunNumber",&fTRunNumber,"fTRunNumber/I");
 foTTree->Branch("Vz",&fTVz,"fTVz/F");
 foTTree->Branch("nCascade",&fTnCascade,"fTnCascade/I");
 foTTree->Branch("CascadePt",&fTCascadePt,"fTCascadePt[fTnCascade]/F");
 foTTree->Branch("CascadeCharge",&fTCascadeCharge,"fTCascadeCharge[fTnCascade]/S");
 foTTree->Branch("CascadeDCA",&fTCascadeDCA,"fTCascadeDCA[fTnCascade]/F");
 foTTree->Branch("CascadeDaughtersDCA",&fTCascadeDaughtersDCA,"fTCascadeDaughtersDCA[fTnCascade]/F");
 foTTree->Branch("CascadeXiMass",&fTCascadeXiMass,"fTCascadeXiMass[fTnCascade]/F");
 foTTree->Branch("CascadeOmegaMass",&fTCascadeOmegaMass,"fTCascadeOmegaMass[fTnCascade]/F");
 foTTree->Branch("CascadeVr",&fTCascadeVr,"fTCascadeVr[fTnCascade]/F");
 foTTree->Branch("CascadePA",&fTCascadePA,"fTCascadePA[fTnCascade]/F");
 foTTree->Branch("LambdaPt",&fTLambdaPt,"fTLambdaPt[fTnCascade]/F");
 foTTree->Branch("LambdaDCA",&fTLambdaDCA,"fTLambdaDCA[fTnCascade]/F");
 foTTree->Branch("LambdaDaughtersDCA",&fTLambdaDaughtersDCA,"fTLambdaDaughtersDCA[fTnCascade]/F");
 foTTree->Branch("LambdaMass",&fTLambdaMass,"fTLambdaMass[fTnCascade]/F");
 foTTree->Branch("LambdaK0Mass",&fTLambdaK0Mass,"fTLambdaK0Mass[fTnCascade]/F");
 foTTree->Branch("LambdaVr",&fTLambdaVr,"fTLambdaVr[fTnCascade]/F");
 foTTree->Branch("LambdaPA",&fTLambdaPA,"fTLambdaPA[fTnCascade]/F");
 foTTree->Branch("TrackP",&fTTrackP,"fTTrackP[fTnCascade][3]/F");
 foTTree->Branch("TrackTPCmom",&fTTrackTPCmom,"fTTrackTPCmom[fTnCascade][3]/F");
 //foTTree->Branch("TrackEta",&fTTrackEta,"fTTrackEta[fTnCascade][3]/F");
 foTTree->Branch("TrackCharge",&fTTrackCharge,"fTTrackCharge[fTnCascade][3]/S");
 foTTree->Branch("TrackDCA",&fTTrackDCA,"fTTrackDCA[fTnCascade][3]/F");
 foTTree->Branch("TrackITSspi",&fTTrackITSspi,"fTTrackITSspi[fTnCascade][3]/F");
 foTTree->Branch("TrackITSsk",&fTTrackITSsk,"fTTrackITSsk[fTnCascade][3]/F");
 foTTree->Branch("TrackITSsp",&fTTrackITSsp,"fTTrackITSsp[fTnCascade][3]/F");
 foTTree->Branch("TrackTPCspi",&fTTrackTPCspi,"fTTrackTPCspi[fTnCascade][3]/F");
 foTTree->Branch("TrackTPCsk",&fTTrackTPCsk,"fTTrackTPCsk[fTnCascade][3]/F");
 foTTree->Branch("TrackTPCsp",&fTTrackTPCsp,"fTTrackTPCsp[fTnCascade][3]/F");
 foTTree->Branch("TrackTOFspi",&fTTrackTOFspi,"fTTrackTOFspi[fTnCascade][3]/F");
 foTTree->Branch("TrackTOFsk",&fTTrackTOFsk,"fTTrackTOFsk[fTnCascade][3]/F");
 foTTree->Branch("TrackTOFsp",&fTTrackTOFsp,"fTTrackTOFsp[fTnCascade][3]/F");
 foTTree->Branch("TrackITStime",&fTTrackITStime,"fTTrackITStime[fTnCascade][3]/O");
 foTTree->Branch("TrackTOFtime",&fTTrackTOFtime,"fTTrackTOFtime[fTnCascade][3]/O");
 foTTree->Branch("TrackTPConly",&fTTrackTPConly,"fTTrackTPConly[fTnCascade][3]/O");
 foTTree->Branch("TrackITScomplementary",&fTTrackITScomplementary,"fTTrackITScomplementary[fTnCascade][3]/O");
 foTTree->Branch("TrackITSpure",&fTTrackITSpure,"fTTrackITSpure[fTnCascade][3]/O");
 foTTree->Branch("TrackGLOBAL",&fTTrackGLOBAL,"fTTrackGLOBAL[fTnCascade][3]/O");

 fomegaTTree = new TTree("omegaTTree","a simple TTree of omega leafs and life truths");
 fomegaTTree->SetDirectory(0); // This is to force a memory-resident Tree, and avoid errors. // ????? is this necessary? does it create memory problems?
 fomegaTTree->Branch("RunNumber",&fTRunNumber,"fTRunNumber/I");
 //fomegaTTree->Branch("Vx",&fTVx,"fTVx/F");
 //fomegaTTree->Branch("Vy",&fTVy,"fTVy/F");
 fomegaTTree->Branch("Vz",&fTVz,"fTVz/F");
 fomegaTTree->Branch("Mult",&fTMult,"fTMult/I");
 //protons:
 fomegaTTree->Branch("nProton",&fTnProton,"fTnProton/I");
 //fomegaTTree->Branch("ProtonP",&fTProtonP,"fTProtonP[fTnProton]/F");
 //fomegaTTree->Branch("ProtonPt",&fTProtonPt,"fTProtonPt[fTnProton]/F");
 fomegaTTree->Branch("ProtonPx",&fTProtonPx,"fTProtonPx[fTnProton]/F");
 fomegaTTree->Branch("ProtonPy",&fTProtonPy,"fTProtonPy[fTnProton]/F");
 fomegaTTree->Branch("ProtonPz",&fTProtonPz,"fTProtonPz[fTnProton]/F");
 fomegaTTree->Branch("ProtonVPx",&fTProtonVPx,"fTProtonVPx[fTnProton]/F");
 fomegaTTree->Branch("ProtonVPy",&fTProtonVPy,"fTProtonVPy[fTnProton]/F");
 fomegaTTree->Branch("ProtonVPz",&fTProtonVPz,"fTProtonVPz[fTnProton]/F");
 fomegaTTree->Branch("ProtonEta",&fTProtonEta,"fTProtonEta[fTnProton]/F");
 //fomegaTTree->Branch("ProtonmT",&fTProtonmT,"fTProtonmT[fTnProton]/F");
 //fomegaTTree->Branch("ProtonTPCmom",&fTProtonTPCmom,"fTProtonTPCmom[fTnProton]/F");
 fomegaTTree->Branch("ProtonCharge",&fTProtonCharge,"fTProtonCharge[fTnProton]/S");
 fomegaTTree->Branch("ProtonTPCsp",&fTProtonTPCsp,"fTProtonTPCsp[fTnProton]/F");
 fomegaTTree->Branch("ProtonTOFsp",&fTProtonTOFsp,"fTProtonTOFsp[fTnProton]/F");
 fomegaTTree->Branch("ProtonDCA",&fTProtonDCA,"fTProtonDCA[fTnProton]/F");
 fomegaTTree->Branch("ProtonNcl",&fTProtonNcl,"fTProtonNcl[fTnProton]/I");
 //fomegaTTree->Branch("ProtonCrF",&fTProtonCrF,"fTProtonCrF[fTnProton]/F");
 //fomegaTTree->Branch("ProtonShared",&fTProtonShared,"fTProtonShared[fTnProton]/I");
 //fomegaTTree->Branch("ProtonTPCchi2",&fTProtonTPCchi2,"fTProtonTPCchi2[fTnProton]/F");
 fomegaTTree->Branch("ProtonITStime",&fTProtonITStime,"fTProtonITStime[fTnProton]/O");
 fomegaTTree->Branch("ProtonTOFtime",&fTProtonTOFtime,"fTProtonTOFtime[fTnProton]/O");
 fomegaTTree->Branch("ProtonTPConly",&fTProtonTPConly,"fTProtonTPConly[fTnProton]/O");
 //fomegaTTree->Branch("ProtonITScomplementary",&fTProtonITScomplementary,"fTProtonITScomplementary[fTnProton]/O");
 //fomegaTTree->Branch("ProtonITSpure",&fTProtonITSpure,"fTProtonITSpure[fTnProton]/O");
 fomegaTTree->Branch("ProtonGLOBAL",&fTProtonGLOBAL,"fTProtonGLOBAL[fTnProton]/O");
 //fomegaTTree->Branch("ProtonFilterBit",&fTProtonFilterBit,"fTProtonFilterBit[fTnProton]/i");
 //omegas:
 fomegaTTree->Branch("nCascade",&fTnCascade,"fTnCascade/I");
 //fomegaTTree->Branch("CascadePx",&fTCascadePx,"fTCascadePx[fTnCascade]/F");
 //fomegaTTree->Branch("CascadePy",&fTCascadePy,"fTCascadePy[fTnCascade]/F");
 //fomegaTTree->Branch("CascadePz",&fTCascadePz,"fTCascadePz[fTnCascade]/F");
 //fomegaTTree->Branch("CascadePt",&fTCascadePt,"fTCascadePt[fTnCascade]/F");
 //fomegaTTree->Branch("CascademT",&fTCascademT,"fTCascademT[fTnCascade]/F");
 fomegaTTree->Branch("CascadeCharge",&fTCascadeCharge,"fTCascadeCharge[fTnCascade]/S");
 fomegaTTree->Branch("CascadeDCA",&fTCascadeDCA,"fTCascadeDCA[fTnCascade]/F");
 fomegaTTree->Branch("CascadeDaughtersDCA",&fTCascadeDaughtersDCA,"fTCascadeDaughtersDCA[fTnCascade]/F");
 fomegaTTree->Branch("CascadeXiMass",&fTCascadeXiMass,"fTCascadeXiMass[fTnCascade]/F");
 fomegaTTree->Branch("CascadeOmegaMass",&fTCascadeOmegaMass,"fTCascadeOmegaMass[fTnCascade]/F");
 fomegaTTree->Branch("CascadeVr",&fTCascadeVr,"fTCascadeVr[fTnCascade]/F");
 //fomegaTTree->Branch("CascadeVx",&fTCascadeVx,"fTCascadeVx[fTnCascade]/F");
 //fomegaTTree->Branch("CascadeVy",&fTCascadeVy,"fTCascadeVy[fTnCascade]/F");
 //fomegaTTree->Branch("CascadeVz",&fTCascadeVz,"fTCascadeVz[fTnCascade]/F");
 fomegaTTree->Branch("CascadePA",&fTCascadePA,"fTCascadePA[fTnCascade]/F");
 //fomegaTTree->Branch("LambdaPt",&fTLambdaPt,"fTLambdaPt[fTnCascade]/F");
 fomegaTTree->Branch("LambdaDCA",&fTLambdaDCA,"fTLambdaDCA[fTnCascade]/F");
 fomegaTTree->Branch("LambdaDaughtersDCA",&fTLambdaDaughtersDCA,"fTLambdaDaughtersDCA[fTnCascade]/F");
 fomegaTTree->Branch("LambdaMass",&fTLambdaMass,"fTLambdaMass[fTnCascade]/F");
 //fomegaTTree->Branch("LambdaK0Mass",&fTLambdaK0Mass,"fTLambdaK0Mass[fTnCascade]/F");
 //fomegaTTree->Branch("LambdaVx",&fTLambdaVx,"fTLambdaVx[fTnCascade]/F");
 //fomegaTTree->Branch("LambdaVy",&fTLambdaVy,"fTLambdaVy[fTnCascade]/F");
 //fomegaTTree->Branch("LambdaVz",&fTLambdaVz,"fTLambdaVz[fTnCascade]/F");
 fomegaTTree->Branch("LambdaVr",&fTLambdaVr,"fTLambdaVr[fTnCascade]/F");
 fomegaTTree->Branch("LambdaPA",&fTLambdaPA,"fTLambdaPA[fTnCascade]/F");
 //fomegaTTree->Branch("TrackP",&fTTrackP,"fTTrackP[fTnCascade][3]/F");
 fomegaTTree->Branch("TrackPx",&fTTrackPx,"fTTrackPx[fTnCascade][3]/F");
 fomegaTTree->Branch("TrackPy",&fTTrackPy,"fTTrackPy[fTnCascade][3]/F");
 fomegaTTree->Branch("TrackPz",&fTTrackPz,"fTTrackPz[fTnCascade][3]/F");
 //fomegaTTree->Branch("TrackTPCmom",&fTTrackTPCmom,"fTTrackTPCmom[fTnCascade][3]/F");
 fomegaTTree->Branch("TrackEta",&fTTrackEta,"fTTrackEta[fTnCascade][3]/F");
 fomegaTTree->Branch("TrackCharge",&fTTrackCharge,"fTTrackCharge[fTnCascade][3]/S");
 fomegaTTree->Branch("TrackDCA",&fTTrackDCA,"fTTrackDCA[fTnCascade][3]/F");
 //fomegaTTree->Branch("TrackITSspi",&fTTrackITSspi,"fTTrackITSspi[fTnCascade][3]/F");
 //fomegaTTree->Branch("TrackITSsk",&fTTrackITSsk,"fTTrackITSsk[fTnCascade][3]/F");
 //fomegaTTree->Branch("TrackITSsp",&fTTrackITSsp,"fTTrackITSsp[fTnCascade][3]/F");
 fomegaTTree->Branch("TrackTPCspi",&fTTrackTPCspi,"fTTrackTPCspi[fTnCascade][3]/F");
 fomegaTTree->Branch("TrackTPCsk",&fTTrackTPCsk,"fTTrackTPCsk[fTnCascade][3]/F");
 fomegaTTree->Branch("TrackTPCsp",&fTTrackTPCsp,"fTTrackTPCsp[fTnCascade][3]/F");
 //fomegaTTree->Branch("TrackTOFspi",&fTTrackTOFspi,"fTTrackTOFspi[fTnCascade][3]/F");
 //fomegaTTree->Branch("TrackTOFsk",&fTTrackTOFsk,"fTTrackTOFsk[fTnCascade][3]/F");
 //fomegaTTree->Branch("TrackTOFsp",&fTTrackTOFsp,"fTTrackTOFsp[fTnCascade][3]/F");
 fomegaTTree->Branch("TrackNcl",&fTTrackNcl,"fTTrackNcl[fTnCascade][3]/I");
 fomegaTTree->Branch("TrackCrF",&fTTrackCrF,"fTTrackCrF[fTnCascade][3]/F");
 fomegaTTree->Branch("TrackShared",&fTTrackShared,"fTTrackShared[fTnCascade][3]/I");
 //fomegaTTree->Branch("TrackTPCchi2",&fTTrackTPCchi2,"fTTrackTPCchi2[fTnCascade][3]/F");
 fomegaTTree->Branch("TrackITStime",&fTTrackITStime,"fTTrackITStime[fTnCascade][3]/O");
 fomegaTTree->Branch("TrackTOFtime",&fTTrackTOFtime,"fTTrackTOFtime[fTnCascade][3]/O");
 fomegaTTree->Branch("TrackTPConly",&fTTrackTPConly,"fTTrackTPConly[fTnCascade][3]/O");
 //fomegaTTree->Branch("TrackITScomplementary",&fTTrackITScomplementary,"fTTrackITScomplementary[fTnCascade][3]/O");
 //fomegaTTree->Branch("TrackITSpure",&fTTrackITSpure,"fTTrackITSpure[fTnCascade][3]/O");
 fomegaTTree->Branch("TrackGLOBAL",&fTTrackGLOBAL,"fTTrackGLOBAL[fTnCascade][3]/O");
 //fomegaTTree->Branch("TrackFilterBit",&fTTrackFilterBit,"fTTrackFilterBit[fTnCascade][3]/i");

 //intialize the random:
 frndm = new TRandom3();

}


void AliOtonOmegaAnalysis::ResetGlobalTrackReference() {
  //This method was inherited form H. Beck analysis

  // Sets all the pointers to zero. To be called at
  // the beginning or end of an event
  for (UShort_t i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}
void AliOtonOmegaAnalysis::StoreGlobalTrackReference(AliAODTrack *track) {
  //This method was inherited form H. Beck analysis

  //bhohlweg@cern.ch: We ask for the Unique Track ID that points back to the
  //ESD. Seems like global tracks have a positive ID, Tracks with Filterbit
  //128 only have negative ID, this is used to match the Tracks later to their
  //global counterparts

  // Stores the pointer to the global track

  // This was AOD073
  // // Don't use the filter bits 2 (ITS standalone) and 128 TPC only
  // // Remove this return statement and you'll see they don't have
  // // any TPC signal
  // if(track->TestFilterBit(128) || track->TestFilterBit(2))
  //   return;
  // This is AOD086
  // Another set of tracks was introduced: Global constrained.
  // We only want filter bit 1 <-- NO! we also want no
  // filter bit at all, which are the v0 tracks
  //  if(!track->TestFilterBit(1))
  //    return;

  // There are also tracks without any filter bit, i.e. filter map 0,
  // at the beginning of the event: they have ~id 1 to 5, 1 to 12
  // This are tracks that didn't survive the primary track filter but
  // got written cause they are V0 daughters

  // Check whether the track has some info
  // I don't know: there are tracks with filter bit 0
  // and no TPC signal. ITS standalone V0 daughters?
  // if(!track->GetTPCsignal()){
  //   printf("Warning: track has no TPC signal, "
  //     //    "not adding it's info! "
  //     "ID: %d FilterMap: %d\n"
  //     ,track->GetID(),track->GetFilterMap());
  //   //    return;
  // }

  // Check that the id is positive
  const int trackID = track->GetID();
  if (trackID < 0) {
    return;
  }

  // Check id is not too big for buffer
  if (trackID >= fTrackBufferSize) {
    printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n", trackID,
           fTrackBufferSize);
    return;
  }

  // Warn if we overwrite a track
  if (fGTI[trackID]) {
    // Seems like there are FilterMap 0 tracks
    // that have zero TPCNcls, don't store these!
    if ((!track->GetFilterMap()) && (!track->GetTPCNcls())) {
      return;
    }
    // Imagine the other way around, the zero map zero clusters track
    // is stored and the good one wants to be added. We ommit the warning
    // and just overwrite the 'bad' track
    if (fGTI[trackID]->GetFilterMap() || fGTI[trackID]->GetTPCNcls()) {
      // If we come here, there's a problem
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
             (fGTI[trackID])->GetTPCNcls(), track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
             (fGTI[trackID])->GetFilterMap(), track->GetFilterMap());
    }
  }  // Two tracks same id

  // // There are tracks with filter bit 0,
  // // do they have TPCNcls stored?
  // if(!track->GetFilterMap()){
  //   printf("Filter map is zero, TPCNcls: %u\n"
  //     ,track->GetTPCNcls());
  // }

  // Assign the pointer
  (fGTI[trackID]) = track;
}



void AliOtonOmegaAnalysis::Make(AliAODEvent *evt, bool OmegaTreeFlag) {
  if (!evt) {
    AliFatal("No Input Event");
  }
  fEvent->SetEvent(evt);
  if (!fEvtCuts->isSelected(fEvent)) {
    return;
  }

  //Initialize Tree
  InitializeTreeValues();


  //start filling Tree, event properties and other things:
  fTRunNumber = evt->GetRunNumber();
  Double_t PrimVtx[3];
  evt->GetPrimaryVertex()->GetXYZ(PrimVtx);
  //fTV[0]=PrimVtx[0];
  //fTV[1]=PrimVtx[1];
  fTVz=PrimVtx[2];
  fTVx=PrimVtx[0];
  fTVy=PrimVtx[1];

 
 
  //prepare aod tracks:
  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < evt->GetNumberOfTracks(); ++iTrack) {
    AliAODTrack *track = static_cast<AliAODTrack*>(evt->GetTrack(iTrack));
    if (!track) {
      AliFatal("No Standard AOD");
      return;
    }
    StoreGlobalTrackReference(track);
  }


  //cascade loop
  std::vector<AliFemtoDreamBasePart> XiDecays;
  std::vector<AliFemtoDreamBasePart> AntiXiDecays;
  std::vector<AliFemtoDreamBasePart> XiOmegaDecays;
  std::vector<AliFemtoDreamBasePart> AntiXiOmegaDecays;
  int numcascades = evt->GetNumberOfCascades();
  for (int iXi = 0; iXi < numcascades; ++iXi) {
    AliAODcascade *xi = evt->GetCascade(iXi);
    if (!xi)
      continue;

    //By default Omegas are not written in the tree
    FillOmegaTree_omega = kFALSE;
    FillOmegaTree_aomega = kFALSE;

    fFemtoCasc->SetCascade(evt, xi);
    if (fCascCuts->isSelected(fFemtoCasc)) {
      XiDecays.push_back(*fFemtoCasc);
    }
    if (fAntiCascCuts->isSelected(fFemtoCasc)) {
      AntiXiDecays.push_back(*fFemtoCasc);
    }
    if (fCascOmegaCuts->isSelected(fFemtoCasc)) {
      XiOmegaDecays.push_back(*fFemtoCasc);
      FillOmegaTree_omega = kTRUE;
    }
    if (fAntiCascOmegaCuts->isSelected(fFemtoCasc)) {
      AntiXiOmegaDecays.push_back(*fFemtoCasc);
      FillOmegaTree_aomega = kTRUE;
    }

    //if an omega passed histos cuts, fill the tree:
    if(FillOmegaTree_omega||FillOmegaTree_aomega){
     Bool_t CascadeFilled = FillTreeCascadeAOD(evt,xi);
     fTnCascade++;
    }
  }//cascades loop


  //tracks loop
  //For proton-Omega exclusive analysis, go on with the proton loop only if we have selected at least one omega candidate:
  std::vector<AliFemtoDreamBasePart> Particles;
  std::vector<AliFemtoDreamBasePart> AntiParticles;
  if(fTnCascade>0){
   fFemtoTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
   for (int iTrack = 0; iTrack < evt->GetNumberOfTracks(); ++iTrack) {
     AliAODTrack *track = static_cast<AliAODTrack*>(evt->GetTrack(iTrack));
     if (!track) {
       AliFatal("No Standard AOD");
       return;
     }
     fFemtoTrack->SetTrack(track, fEvent->GetMultiplicity());
     if (fTrackCuts->isSelected(fFemtoTrack)) {
       Particles.push_back(*fFemtoTrack);
     }
     if (fAntiTrackCuts->isSelected(fFemtoTrack)) {
       AntiParticles.push_back(*fFemtoTrack);
     }
   }
  }


  //skip v0's
  /*
  std::vector<AliFemtoDreamBasePart> Decays;
  std::vector<AliFemtoDreamBasePart> AntiDecays;
  //  Look for the lambda, store it in an event
  //  Get a V0 from the event:
  TClonesArray *v01 = static_cast<TClonesArray*>(evt->GetV0s());
  //number of V0s:

  fFemtov0->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  int entriesV0 = v01->GetEntriesFast();
  for (int iv0 = 0; iv0 < entriesV0; iv0++) {
    AliAODv0 *v0 = evt->GetV0(iv0);
    fFemtov0->Setv0(evt, v0, fEvent->GetMultiplicity());
    if (fv0Cuts->isSelected(fFemtov0)) {
      Decays.push_back(*fFemtov0);
    }
    if (fAntiv0Cuts->isSelected(fFemtov0)) {
      AntiDecays.push_back(*fFemtov0);
    }
  }
  //  std::cout << "=====================================\n" ;
  //  std::cout << "=====================================\n" ;
  //  std::cout << "==========New event==================\n" ;
  //  std::cout << "=====================================\n" ;
  //  std::cout << "=====================================\n" ;
  */



  //forget for the moment about MC:
  /*
  //loop once over the MC stack to calculate Efficiency/Purity
  if (fIsMC) {
    AliAODInputHandler *eventHandler =
        dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()
            ->GetInputEventHandler());
    AliMCEvent* fMC = eventHandler->MCEvent();

    for (int iPart = 0; iPart < (fMC->GetNumberOfTracks()); iPart++) {
      AliAODMCParticle *mcPart = (AliAODMCParticle*) fMC->GetTrack(iPart);
      if (TMath::Abs(mcPart->Eta()) < 0.8 && mcPart->IsPhysicalPrimary()) {
        if (mcPart->GetPdgCode() == fTrackCuts->GetPDGCode()) {
          fTrackCuts->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fAntiTrackCuts->GetPDGCode()) {
          fAntiTrackCuts->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fv0Cuts->GetPDGv0()) {
          fv0Cuts->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fAntiv0Cuts->GetPDGv0()) {
          fAntiv0Cuts->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fCascCuts->GetPDGCodeCasc()) {
          fCascCuts->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fAntiCascCuts->GetPDGCodeCasc()) {
          fAntiCascCuts->FillGenerated(mcPart->Pt());
        }
      }
    }
  }
  */

  //if we have selected at least one cascade and tree flag is on, fill the omega tree:
  if(OmegaTreeFlag&&fTnCascade>0) fomegaTTree->Fill();


  //pairs stuff
  fPairCleaner->ResetArray();
  //fPairCleaner->CleanTrackAndDecay(&Particles, &Decays, 0);
  fPairCleaner->CleanTrackAndDecay(&Particles, &XiDecays, 2);
  fPairCleaner->CleanTrackAndDecay(&Particles, &XiOmegaDecays, 0);
  //fPairCleaner->CleanTrackAndDecay(&AntiParticles, &AntiDecays, 1);
  fPairCleaner->CleanTrackAndDecay(&AntiParticles, &AntiXiDecays, 3);
  fPairCleaner->CleanTrackAndDecay(&AntiParticles, &AntiXiOmegaDecays, 1);

  //fPairCleaner->CleanDecay(&Decays, 0);
  //fPairCleaner->CleanDecay(&AntiDecays, 1);
  fPairCleaner->CleanDecay(&XiDecays, 2);
  fPairCleaner->CleanDecay(&AntiXiDecays, 3);
  fPairCleaner->CleanDecay(&XiOmegaDecays, 0);
  fPairCleaner->CleanDecay(&AntiXiOmegaDecays, 1);

  fPairCleaner->StoreParticle(Particles);
  fPairCleaner->StoreParticle(AntiParticles);
  //fPairCleaner->StoreParticle(Decays);
  //fPairCleaner->StoreParticle(AntiDecays);
  fPairCleaner->StoreParticle(XiDecays);
  fPairCleaner->StoreParticle(AntiXiDecays);
  fPairCleaner->StoreParticle(XiOmegaDecays); // comes in the same order as in the addtask
  fPairCleaner->StoreParticle(AntiXiOmegaDecays);// comes in the same order as in the addtask


  if (fConfig->GetUseEventMixing()) {
    fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(),
                        fEvent->GetMultiplicity(), fEvent->GetV0MCentrality());
  }
  if (fConfig->GetUsePhiSpinning()) {
    fControlSample->SetEvent(fPairCleaner->GetCleanParticles(),fEvent);
  }
}

void AliOtonOmegaAnalysis::Make(AliESDEvent *evt, AliMCEvent *mcEvent, bool CascadeTreeFlag, bool OmegaTreeFlag, Int_t Cut) {
  if (!evt) {
    AliFatal("No Input Event");
  }
  fEvent->SetEvent(evt);
  if (!fEvtCuts->isSelected(fEvent)) {
    return;
  }

  //Initialize Tree
  InitializeTreeValues();


  //start filling Tree, event properties and other things:
  fTRunNumber = evt->GetRunNumber();
  Double_t PrimVtx[3];
  evt->GetPrimaryVertex()->GetXYZ(PrimVtx);
  //fTV[0]=PrimVtx[0];
  //fTV[1]=PrimVtx[1];
  fTVz=PrimVtx[2];
  fTVx=PrimVtx[0];
  fTVy=PrimVtx[1];
  fTMult = fEvent->GetMultiplicity();

  //skip v0's
/*
  std::vector<AliFemtoDreamBasePart> Decays;
  std::vector<AliFemtoDreamBasePart> AntiDecays;

  for (int iv0 = 0; iv0 < evt->GetNumberOfV0s(); ++iv0) {
    AliESDv0 *v0 = evt->GetV0(iv0);
    fFemtov0->Setv0(evt, mcEvent, v0, fEvent->GetMultiplicity());
    if (fv0Cuts->isSelected(fFemtov0)) {
      Decays.push_back(*fFemtov0);
    }
    if (fAntiv0Cuts->isSelected(fFemtov0)) {
      AntiDecays.push_back(*fFemtov0);
    }
  }
*/

  //proton loop
  std::vector<AliFemtoDreamBasePart> Particles;
  std::vector<AliFemtoDreamBasePart> AntiParticles;
  for (int iTrack = 0; iTrack < evt->GetNumberOfTracks(); ++iTrack) {
    AliESDtrack *track = static_cast<AliESDtrack *>(evt->GetTrack(iTrack));
    fFemtoTrack->SetTrack(track, mcEvent, fEvent->GetMultiplicity());
    if (fTrackCuts->isSelected(fFemtoTrack)) {
      Particles.push_back(*fFemtoTrack); 
      Bool_t protonfilled = FillProtonTrack(evt, iTrack);
      fTnProton++;
    }
    if (fAntiTrackCuts->isSelected(fFemtoTrack)) {
      AntiParticles.push_back(*fFemtoTrack);
      Bool_t protonfilled = FillProtonTrack(evt, iTrack);
      fTnProton++;
    }
  }


  //Cascades loop
  std::vector<AliFemtoDreamBasePart> XiOmegaDecays;
  std::vector<AliFemtoDreamBasePart> AntiXiOmegaDecays;
  std::vector<AliFemtoDreamBasePart> XiDecays;
  std::vector<AliFemtoDreamBasePart> AntiXiDecays;
  for (Int_t nCascade = 0; nCascade < evt->GetNumberOfCascades(); ++nCascade) {
    AliESDcascade *esdCascade = evt->GetCascade(nCascade);

    //By default Omegas are not written in the tree
    FillOmegaTree_omega = kFALSE;
    FillOmegaTree_aomega = kFALSE;

    //Xi block goes now to bkg block:
    fFemtoCasc->SetCascade(evt, mcEvent, esdCascade);
    if (fCascCuts->isSelected(fFemtoCasc)) {
      XiDecays.push_back(*fFemtoCasc);
      //FillOmegaTree_omega = kTRUE;
    }
    if (fAntiCascCuts->isSelected(fFemtoCasc)) {
      AntiXiDecays.push_back(*fFemtoCasc);
      //FillOmegaTree_aomega = kTRUE;
    }

    //signal omegas:
    if (fCascOmegaCuts->isSelected(fFemtoCasc)) {
      XiOmegaDecays.push_back(*fFemtoCasc);
      FillOmegaTree_omega = kTRUE;
    }
    if (fAntiCascOmegaCuts->isSelected(fFemtoCasc)) {
      AntiXiOmegaDecays.push_back(*fFemtoCasc);
      FillOmegaTree_aomega = kTRUE;
    }

    //if an omega passed histos cuts, fill the tree:
    if(FillOmegaTree_omega||FillOmegaTree_aomega){
     Bool_t CascadeFilled = FillTreeCascade(evt, esdCascade);
     fTnCascade++;
    }
  }//cascades loop


  //if we have selected at least one cascade and tree flag is on, fill the omega tree:
  //if(OmegaTreeFlag&&fTnCascade>0) fomegaTTree->Fill();
  // -> Now fill also 1% of events with protons
  Float_t r3=frndm->Rndm();
  if(OmegaTreeFlag&&( fTnCascade>0 || (fTnProton&&r3<.01) )) fomegaTTree->Fill();

  //pairing stuff:
  fPairCleaner->ResetArray();
  //fPairCleaner->CleanTrackAndDecay(&Particles, &Decays, 0);
  fPairCleaner->CleanTrackAndDecay(&Particles, &XiDecays, 2);
  fPairCleaner->CleanTrackAndDecay(&Particles, &XiOmegaDecays, 0);
  //fPairCleaner->CleanTrackAndDecay(&AntiParticles, &AntiDecays, 1);
  fPairCleaner->CleanTrackAndDecay(&AntiParticles, &AntiXiDecays, 3);
  fPairCleaner->CleanTrackAndDecay(&AntiParticles, &AntiXiOmegaDecays, 1);

  //fPairCleaner->CleanDecay(&Decays, 0);
  //fPairCleaner->CleanDecay(&AntiDecays, 1);
  fPairCleaner->CleanDecay(&XiDecays, 2);
  fPairCleaner->CleanDecay(&AntiXiDecays, 3);
  fPairCleaner->CleanDecay(&XiOmegaDecays, 0);
  fPairCleaner->CleanDecay(&AntiXiOmegaDecays, 1);

  fPairCleaner->StoreParticle(Particles);
  fPairCleaner->StoreParticle(AntiParticles);
  //fPairCleaner->StoreParticle(Decays);
  //fPairCleaner->StoreParticle(AntiDecays);
  fPairCleaner->StoreParticle(XiDecays);
  fPairCleaner->StoreParticle(AntiXiDecays);
  fPairCleaner->StoreParticle(XiOmegaDecays); // comes in the same order as in the addtask
  fPairCleaner->StoreParticle(AntiXiOmegaDecays); // comes in the same order as in the addtask

  if (fConfig->GetUseEventMixing()) {
    fPartColl->SetEvent(fPairCleaner->GetCleanParticles(), fEvent->GetZVertex(),
                        fEvent->GetMultiplicity(), fEvent->GetV0MCentrality());
  }
  if (fConfig->GetUsePhiSpinning()) {
    fControlSample->SetEvent(fPairCleaner->GetCleanParticles(),fEvent);
  }
}






//_______________________________________________________________________________________________________
Bool_t AliOtonOmegaAnalysis::FillTreeCascade(AliESDEvent *evt, AliESDcascade *casc) {

 Bool_t Filled = kFALSE;
 Double_t PrimVtx[3];
 evt->GetPrimaryVertex()->GetXYZ(PrimVtx);

 //First of all check if the V0 of the Cascade is good:
 int idxPosFromV0Dghter = casc->GetPindex();
 int idxNegFromV0Dghter = casc->GetNindex();
 AliESDv0 * currentV0 = 0x0;
 int idxV0FromCascade = -1;
 for (int iV0 = 0; iV0 < evt->GetNumberOfV0s(); ++iV0) {
  currentV0 = evt->GetV0(iV0);
  // Make sure not to use online v0 here
  // (the cascade is reconstructed with offline v0 only)
  if (currentV0->GetOnFlyStatus() == kTRUE) continue;
  int posCurrentV0 = currentV0->GetPindex();
  int negCurrentV0 = currentV0->GetNindex();
  if (posCurrentV0 == idxPosFromV0Dghter&& negCurrentV0 == idxNegFromV0Dghter) {
   idxV0FromCascade = iV0;
   break;
  }
 }

 if (idxV0FromCascade >= 0) { //we got a correct Lambda
 //----------------------------------------------------


  //fill the tracks
  Bool_t tr0 = FillTreeTrack(0, casc->GetPindex() , idxV0FromCascade, evt , casc);
  Bool_t tr1 = FillTreeTrack(1, casc->GetNindex() , idxV0FromCascade, evt , casc);
  Bool_t tr2 = FillTreeTrack(2, casc->GetBindex() , idxV0FromCascade, evt , casc);

  //tracks are filled
  //-----------------
  if(tr0&&tr1&&tr2) {
cout<<"tracks filled "<<endl;

   //get the Lambda:
   AliESDv0* v0 = evt->GetV0(idxV0FromCascade);

   //Check the Xi/Omega mass mass:
   //Lambda E:
   Double_t Ev0 = sqrt( pow(TDatabasePDG::Instance()->GetParticle(3122)->Mass(),2) + (v0->Px()*v0->Px()+v0->Py()*v0->Py()+v0->Pz()*v0->Pz()));
   //Bach E:
   Double_t EBach = sqrt( pow(TDatabasePDG::Instance()->GetParticle(211)->Mass(),2) + fTTrackP[fTnCascade][2]);
   //Cascade momentum:
   Double_t px,py,pz;
   px = fTTrackPx[fTnCascade][0] + fTTrackPx[fTnCascade][1] + fTTrackPx[fTnCascade][2];
   py = fTTrackPy[fTnCascade][0] + fTTrackPy[fTnCascade][1] + fTTrackPy[fTnCascade][2];
   pz = fTTrackPz[fTnCascade][0] + fTTrackPz[fTnCascade][1] + fTTrackPz[fTnCascade][2];
   fTCascadePx[fTnCascade] = px;
   fTCascadePy[fTnCascade] = py;
   fTCascadePz[fTnCascade] = pz;
   fTCascadeP[fTnCascade] = sqrt(px*px+py*py+pz*pz);
   fTCascadePt[fTnCascade] = sqrt(px*px+py*py);
   fTCascademT[fTnCascade] = sqrt( pow(1.67245,2) +   px*px+py*py );
   Double_t mXi = sqrt( pow(Ev0 + EBach, 2) - pow(fTCascadeP[fTnCascade],2) );
   //for bachelor kaon:
   EBach = sqrt( pow(TDatabasePDG::Instance()->GetParticle(321)->Mass(),2) + pow(fTTrackP[fTnCascade][2],2));
   Double_t mOmega = sqrt( pow(Ev0 + EBach, 2) - pow(fTCascadeP[fTnCascade],2));

   //good Xi/Omega mass:
   //-------------------
   if((mXi>1.20&&mXi<1.45)||(mOmega>1.55&&mOmega<1.80)){

    //fill lambda
    fTLambdaPx[fTnCascade] = v0->Px();
    fTLambdaPy[fTnCascade] = v0->Py();
    fTLambdaPz[fTnCascade] = v0->Pz();
    fTLambdaPt[fTnCascade] = sqrt(pow(v0->Px(),2)+pow(v0->Py(),2));
    fTLambdaDCA[fTnCascade] = v0->GetD(PrimVtx[0], PrimVtx[1], PrimVtx[2]);
    fTLambdaDaughtersDCA[fTnCascade] = v0->GetDcaV0Daughters();
    v0->ChangeMassHypothesis(kK0Short);
    fTLambdaK0Mass[fTnCascade]=v0->GetEffMass();
    if(casc->Charge()>0) {
     v0->ChangeMassHypothesis(kLambda0Bar);
    } else {
     v0->ChangeMassHypothesis(kLambda0);
    }
    fTLambdaMass[fTnCascade]=v0->GetEffMass();
    fTLambdaVr[fTnCascade] = sqrt(pow(v0->Xv(),2)+pow(v0->Yv(),2));
    fTLambdaVx[fTnCascade] = v0->Xv();
    fTLambdaVy[fTnCascade] = v0->Yv();
    fTLambdaVz[fTnCascade] = v0->Zv();
    fTLambdaPA[fTnCascade] = v0->GetV0CosineOfPointingAngle(PrimVtx[0], PrimVtx[1], PrimVtx[2]);

    //fill cascade
    fTCascadeCharge[fTnCascade] = casc->Charge();
    fTCascadeDaughtersDCA[fTnCascade] = casc->GetDcaXiDaughters();
    fTCascadeXiMass[fTnCascade] = mXi;
    fTCascadeOmegaMass[fTnCascade] = mOmega;
    double decayPosXi[3] = { 0. };
    casc->GetXYZcascade(decayPosXi[0], decayPosXi[1], decayPosXi[2]);
    fTCascadeVr[fTnCascade] = sqrt(pow(decayPosXi[0],2)+pow(decayPosXi[1],2));
    fTCascadeVx[fTnCascade] = decayPosXi[0];
    fTCascadeVy[fTnCascade] = decayPosXi[1];
    fTCascadeVz[fTnCascade] = decayPosXi[2];
    //calc cascade dca, following recipe of AOD
    fTCascadeDCA[fTnCascade] = casc->GetDcascade(PrimVtx[0], PrimVtx[1], PrimVtx[2]);
    fTCascadePA[fTnCascade] = casc->GetCascadeCosineOfPointingAngle(PrimVtx[0], PrimVtx[1], PrimVtx[2]);

    Filled = kTRUE;

   }//good Xi/Omega mass
  } //tracks filled
 } //correct lambda
 return Filled;
}






//_______________________________________________________________________________________________________________________________
Bool_t AliOtonOmegaAnalysis::FillTreeTrack(Int_t jj, Int_t idtrack, Int_t V0index, AliESDEvent *evt, AliESDcascade *casc) {
//here jj=0 refers to pos, 1 neg, 2 bach

 Bool_t Filled = kFALSE;

 //primary vertex
 Double_t PrimVtx[3];
 evt->GetPrimaryVertex()->GetXYZ(PrimVtx);

 //get ESD track
 AliESDtrack *track = evt->GetTrack(idtrack);

 //fill tpc inner momentum
 fTTrackTPCmom[fTnCascade][jj]= track->GetTPCmomentum();

 //REQUIRE EVERYTHING TO HAVE TPCMOM > 50MeV
 if(track->GetTPCmomentum()>.050){

 //fill dca
 fTTrackDCA[fTnCascade][jj]= track->GetD(PrimVtx[0], PrimVtx[1], evt->GetMagneticField());

 //fill momentum in the correct vertex position:
 double Mom[3] = { 0. };
 //the momenta of the daughters have to be taken at the v0 vertex,
 //the bachelor momenta at the cascade vertex
 if(jj<2){
  AliESDv0 * V0 = 0x0;
  V0 = evt->GetV0(V0index);
  if(jj==0)V0->GetPPxPyPz(Mom[0],Mom[1],Mom[2]);
  if(jj==1)V0->GetNPxPyPz(Mom[0],Mom[1],Mom[2]);
 }else{
  casc->GetBPxPyPz(Mom[0], Mom[1], Mom[2]);
 }
 fTTrackPx[fTnCascade][jj]=Mom[0];
 fTTrackPy[fTnCascade][jj]=Mom[1];
 fTTrackPz[fTnCascade][jj]=Mom[2];
 fTTrackP[fTnCascade][jj]=sqrt(pow(Mom[0],2)+pow(Mom[1],2)+pow(Mom[2],2));

 if(sqrt(pow(Mom[0],2)+pow(Mom[1],2))>0.
  && fabs(fTTrackP[fTnCascade][jj])<1000.
  ){

   //fill Eta:
   fTTrackEta[fTnCascade][jj]= track->Eta();

   //fill charge:
   fTTrackCharge[fTnCascade][jj]=track->Charge();

   //fill #TPCclusters
   fTTrackNcl[fTnCascade][jj]  = track->GetTPCNcls();

   //fill TPC CrF : Crossed rows / findable
   if(track->GetTPCNclsF()>0) fTTrackCrF[fTnCascade][jj]  = track->GetTPCCrossedRows()/track->GetTPCNclsF();

   //fill chi2 TPC
   fTTrackTPCchi2[fTnCascade][jj]  = track->GetTPCchi2();

   //fill #shared TPC clusters
   fTTrackShared[fTnCascade][jj]  = track->GetTPCnclsS();

   //does it have a hit in SPD, SSD or TOF?
   if(track->HasPointOnITSLayer(0)||track->HasPointOnITSLayer(1)||track->HasPointOnITSLayer(4)||track->HasPointOnITSLayer(5)) fTTrackITStime[fTnCascade][jj] = kTRUE;
   if(track->GetTOFBunchCrossing() == 0) fTTrackTOFtime[fTnCascade][jj] = kTRUE;

   //-----------------------------------------------------------------------------------------
   //-----------------------------------------------------------------------------------------
   //-----------------------------------------------------------------------------------------
   //The question here is, is that enough that they have ITS or TOF time (in particular ITS),
   //or do I have to do SOMETHING with this time, i.e. check that is in
   //time with the vertex or with the other tracks of the cascade etc.
   //-----------------------------------------------------------------------------------------
   //-----------------------------------------------------------------------------------------
   //-----------------------------------------------------------------------------------------

   //TPCtrack only?
   if(track->GetStatus()&AliESDtrack::kTPCrefit && track->GetStatus()&AliESDtrack::kTPCin && track->GetStatus()&AliESDtrack::kITSrefit && track->GetStatus()&AliESDtrack::kITSin) fTTrackGLOBAL[fTnCascade][jj]=kTRUE;
   if(track->GetStatus()&AliESDtrack::kTPCrefit && track->GetStatus()&AliESDtrack::kTPCin && !(track->GetStatus()&AliESDtrack::kITSrefit) && !(track->GetStatus()&AliESDtrack::kITSin)) fTTrackTPConly[fTnCascade][jj]=kTRUE;
   if(!(track->GetStatus()&AliESDtrack::kTPCin) && track->GetStatus()&AliESDtrack::kITSrefit && track->GetStatus()&AliESDtrack::kITSin && !(track->GetStatus()&AliESDtrack::kITSpureSA)) fTTrackITScomplementary[fTnCascade][jj]=kTRUE;
   if(track->GetStatus()&AliESDtrack::kITSrefit && track->GetStatus()&AliESDtrack::kITSpureSA) fTTrackITSpure[fTnCascade][jj]=kTRUE;


  //Info not included so far, to be considered: - number of its clusters,  - number of tpc shared clusters,  - its shared clusters,  - track index
  //STILL TO BE CHECKED: - is the dca, momentum, etc evaluated for each track in the correct way?
  //                     - For tpc only tracks, using or not the vertex, etc...
  //                     - The info here stored coincides with the info evaluated by the vertexer?

   //pid
   AliPIDResponse::EDetPidStatus statusPID = fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS, track);
   if (statusPID == AliPIDResponse::kDetPidOk)  {
     fTTrackITSspi[fTnCascade][jj] = fPIDResponse->NumberOfSigmas(AliPIDResponse::kITS, track, AliPID::kPion);
     fTTrackITSsk[fTnCascade][jj] = fPIDResponse->NumberOfSigmas(AliPIDResponse::kITS, track, AliPID::kKaon);
     fTTrackITSsp[fTnCascade][jj] = fPIDResponse->NumberOfSigmas(AliPIDResponse::kITS, track, AliPID::kProton);
    }
   statusPID = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track);
   if (statusPID == AliPIDResponse::kDetPidOk) {
     fTTrackTPCspi[fTnCascade][jj] = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, AliPID::kPion);
     fTTrackTPCsk[fTnCascade][jj] = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, AliPID::kKaon);
     fTTrackTPCsp[fTnCascade][jj] = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, AliPID::kProton);
    }
   statusPID = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track);
   if (statusPID == AliPIDResponse::kDetPidOk)  {
     fTTrackTOFspi[fTnCascade][jj] = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF, track, AliPID::kPion);
     fTTrackTOFsk[fTnCascade][jj] = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF, track, AliPID::kKaon);
     fTTrackTOFsp[fTnCascade][jj] = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF, track, AliPID::kProton);
    }


if(fTTrackTPCsk[fTnCascade][jj]<-999.&&fTTrackITSsk[fTnCascade][jj]>-999.){
 cout<<" the track "<<jj<<" has no TPC (pid="<<fTTrackTPCsk[fTnCascade][jj]<<" and it has ITS, pid="<<fTTrackITSsk[fTnCascade][jj]<<endl;
}

cout<<"track "<<jj<<" set to filled "<<endl;
  Filled = kTRUE;

 }//good momentum
}//tpcmom>50 for all 

 return Filled;
}






//_______________________________________________________________________________________________________________________________
Bool_t AliOtonOmegaAnalysis::FillProtonTrack(AliESDEvent *evt, Int_t idtrack) {
//here jj=0 refers to pos, 1 neg, 2 bach

 Bool_t Filled = kFALSE;

 //get ESD track
 AliESDtrack *track = evt->GetTrack(idtrack);

 //fill tpc inner momentum
 fTProtonTPCmom[fTnProton]= track->GetTPCmomentum();

 fTProtonPx[fTnProton]=track->Px();
 fTProtonPy[fTnProton]=track->Py();
 fTProtonPz[fTnProton]=track->Pz();
 fTProtonP[fTnProton]=sqrt(pow(track->Px(),2)+pow(track->Py(),2)+pow(track->Pz(),2));
 fTProtonPt[fTnProton]=sqrt(pow(track->Px(),2)+pow(track->Py(),2));
 fTProtonmT[fTnProton]= sqrt( pow(0.9382721,2) + pow(track->Px(),2) + pow(track->Py(),2) );


 //fill momentum with vtx constraint:
 AliESDtrack *fESDTPCOnlyTrack;
 fESDTPCOnlyTrack = AliESDtrackCuts::GetTPCOnlyTrack(evt, track->GetID());
 double p[3] = { 0. };
 const AliESDVertex *vtxSPD = evt->GetPrimaryVertexSPD();
 AliExternalTrackParam exParam;
 Bool_t relate = false;
 relate = fESDTPCOnlyTrack->RelateToVertexTPC(vtxSPD, evt->GetMagneticField(), 1e30, &exParam);
 // set the constrained parameters to the track
 fESDTPCOnlyTrack->Set(exParam.GetX(), exParam.GetAlpha(),exParam.GetParameter(), exParam.GetCovariance());
 fESDTPCOnlyTrack->GetPxPyPz(p);
 fTProtonVPx[fTnProton]=p[0];
 fTProtonVPy[fTnProton]=p[1];
 fTProtonVPz[fTnProton]=p[2];



 //fill Eta:
 fTProtonEta[fTnProton]= track->Eta();


 //fill charge GIVEN:
 fTProtonCharge[fTnProton]=track->Charge();

 //DCA
   Double_t PrimVtx[3];
   evt->GetPrimaryVertex()->GetXYZ(PrimVtx);
   fTProtonDCA[fTnProton]= track->GetD(PrimVtx[0], PrimVtx[1], evt->GetMagneticField());


 //fill #TPCclusters
 fTProtonNcl[fTnProton]  = track->GetTPCNcls();

 //fill TPC CrF : Crossed rows / findable
 if(track->GetTPCNclsF()>0) fTProtonCrF[fTnProton]  = track->GetTPCCrossedRows()/track->GetTPCNclsF();

 //fill chi2 TPC
 fTProtonTPCchi2[fTnProton]  = track->GetTPCchi2();

 //fill #shared TPC clusters
 fTProtonShared[fTnProton]  = track->GetTPCnclsS();

 //does it have a hit in SPD, SSD or TOF?
 if(track->HasPointOnITSLayer(0)||track->HasPointOnITSLayer(1)||track->HasPointOnITSLayer(4)||track->HasPointOnITSLayer(5)) fTProtonITStime[fTnProton] = kTRUE;
 if(track->GetTOFBunchCrossing() == 0) fTProtonTOFtime[fTnProton] = kTRUE;

 //TPCProton only?
 if(track->GetStatus()&AliESDtrack::kTPCrefit && track->GetStatus()&AliESDtrack::kTPCin && track->GetStatus()&AliESDtrack::kITSrefit && track->GetStatus()&AliESDtrack::kITSin) fTProtonGLOBAL[fTnProton]=kTRUE;
 if(track->GetStatus()&AliESDtrack::kTPCrefit && track->GetStatus()&AliESDtrack::kTPCin && !(track->GetStatus()&AliESDtrack::kITSrefit) && !(track->GetStatus()&AliESDtrack::kITSin)) fTProtonTPConly[fTnProton]=kTRUE;
 if(!(track->GetStatus()&AliESDtrack::kTPCin) && track->GetStatus()&AliESDtrack::kITSrefit && track->GetStatus()&AliESDtrack::kITSin && !(track->GetStatus()&AliESDtrack::kITSpureSA)) fTProtonITScomplementary[fTnProton]=kTRUE;
 if(track->GetStatus()&AliESDtrack::kITSrefit && track->GetStatus()&AliESDtrack::kITSpureSA) fTProtonITSpure[fTnProton]=kTRUE;

 //pid
 AliPIDResponse::EDetPidStatus statusPID = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track);
 if (statusPID == AliPIDResponse::kDetPidOk) {
  fTProtonTPCsp[fTnProton] = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, AliPID::kProton);
 }
 statusPID = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track);
 if (statusPID == AliPIDResponse::kDetPidOk) {
  fTProtonTOFsp[fTnProton] = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF, track, AliPID::kProton);
 }

 Filled = kTRUE;

 return Filled;
}


//_______________________________________________________________________________________________________________________________
void AliOtonOmegaAnalysis::InitializeTreeValues(){
  for(int ii=0;ii<MAXPROTONS;ii++){
   fTProtonP[ii]=-100000.;
   fTProtonPt[ii]=-100000.;
   fTProtonmT[ii]=-100000.;
   fTProtonEta[ii]=-100000.;
   fTProtonPx[ii]=-100000.;
   fTProtonPy[ii]=-100000.;
   fTProtonPz[ii]=-100000.;
   fTProtonVPx[ii]=-100000.;
   fTProtonVPy[ii]=-100000.;
   fTProtonVPz[ii]=-100000.;
   fTProtonTPCmom[ii]=-100000.;
   fTProtonCharge[ii]=-10;
   fTProtonDCA[ii]=-100000.;
   fTProtonTPCsp[ii]=-100000.;
   fTProtonTOFsp[ii]=-100000.;
   fTProtonNcl[ii]=-100000;
   fTProtonCrF[ii]=-100000.;
   fTProtonShared[ii]=-100000;
   fTProtonTPCchi2[ii]=-100000.;
   fTProtonITStime[ii]=kFALSE;
   fTProtonTOFtime[ii]=kFALSE;
   fTProtonTPConly[ii]=kFALSE;
   fTProtonITScomplementary[ii]=kFALSE;
   fTProtonITSpure[ii]=kFALSE;
   fTProtonGLOBAL[ii]=kFALSE;
   fTProtonFilterBit[ii]=0;
  }
  for(int ii=0;ii<MAXCASCADES;ii++){
   fTCascadePx[ii]=-100000.;
   fTCascadePy[ii]=-100000.;
   fTCascadePz[ii]=-100000.;
   fTCascadeP[ii]=-100000.;
   fTCascadePt[ii]=-100000.;
   fTCascademT[ii]=-100000.;
   fTCascadeCharge[ii]=-10;
   fTCascadeDCA[ii]=-100000.;
   fTCascadeDaughtersDCA[ii]=-100000.;
   fTCascadeXiMass[ii]=-100000.;
   fTCascadeOmegaMass[ii]=-100000.;
   fTCascadeVr[ii]=-100000.;
   fTCascadeVx[ii]=-100000.;
   fTCascadeVy[ii]=-100000.;
   fTCascadeVz[ii]=-100000.;
   fTCascadePA[ii]=-100000.;
   fTLambdaPt[ii]=-100000.;
   fTLambdaPx[ii]=-100000.;
   fTLambdaPy[ii]=-100000.;
   fTLambdaPz[ii]=-100000.;
   fTLambdaDCA[ii]=-100000.;
   fTLambdaDaughtersDCA[ii]=-100000.;
   fTLambdaMass[ii]=-100000.;
   fTLambdaK0Mass[ii]=-100000.;
   fTLambdaVr[ii]=-100000.;
   fTLambdaVx[ii]=-100000.;
   fTLambdaVy[ii]=-100000.;
   fTLambdaVz[ii]=-100000.;
   fTLambdaPA[ii]=-100000.;

   for(int jj=0;jj<3;jj++){
    fTTrackP[ii][jj]=-100000.;
    fTTrackPx[ii][jj]=-100000.;
    fTTrackPy[ii][jj]=-100000.;
    fTTrackPz[ii][jj]=-100000.;
    fTTrackTPCmom[ii][jj]=-100000.;
    fTTrackEta[ii][jj]=-100000.;
    fTTrackCharge[ii][jj]=-10;
    fTTrackDCA[ii][jj]=-100000.;
    fTTrackITSspi[ii][jj]=-100000.;
    fTTrackITSsk[ii][jj]=-100000.;
    fTTrackITSsp[ii][jj]=-100000.;
    fTTrackTPCspi[ii][jj]=-100000.;
    fTTrackTPCsk[ii][jj]=-100000.;
    fTTrackTPCsp[ii][jj]=-100000.;
    fTTrackTOFspi[ii][jj]=-100000.;
    fTTrackTOFsk[ii][jj]=-100000.;
    fTTrackTOFsp[ii][jj]=-100000.;
    fTTrackNcl[ii][jj]=-100000;
    fTTrackCrF[ii][jj]=-100000.;
    fTTrackShared[ii][jj]=-100000;
    fTTrackTPCchi2[ii][jj]=-100000.;
    fTTrackITStime[ii][jj]=kFALSE;
    fTTrackTOFtime[ii][jj]=kFALSE;
    fTTrackTPConly[ii][jj]=kFALSE;
    fTTrackITScomplementary[ii][jj]=kFALSE;
    fTTrackITSpure[ii][jj]=kFALSE;
    fTTrackGLOBAL[ii][jj]=kFALSE;
    fTTrackFilterBit[ii][jj]=0;
   }
  }
  fTnCascade=0;
  fTnProton=0;
}







//_______________________________________________________________________________________________________
Bool_t AliOtonOmegaAnalysis::FillTreeCascadeAOD(AliAODEvent *evt, AliAODcascade *casc) {

 Bool_t Filled = kFALSE;

 Double_t PrimVtx[3];
 evt->GetPrimaryVertex()->GetXYZ(PrimVtx);

 //fill the tracks
 Bool_t tr0 = FillTreeTrackAOD(0, casc);
 Bool_t tr1 = FillTreeTrackAOD(1, casc);
 Bool_t tr2 = FillTreeTrackAOD(2, casc);

 //tracks are filled
 //-----------------
 if(tr0&&tr1&&tr2) {

   fTCascadePx[fTnCascade] = casc->MomXiX();
   fTCascadePy[fTnCascade] = casc->MomXiY();
   fTCascadePz[fTnCascade] = casc->MomXiZ();
   fTCascadePt[fTnCascade] = sqrt(casc->MomXiX()*casc->MomXiX()+casc->MomXiY()*casc->MomXiY());
   fTCascadeP[fTnCascade] = sqrt(casc->MomXiX()*casc->MomXiX() + casc->MomXiY()*casc->MomXiY() + casc->MomXiZ()*casc->MomXiZ() );
   fTCascadeXiMass[fTnCascade] = casc->MassXi();
   fTCascadeOmegaMass[fTnCascade] = casc->MassOmega();
   Double_t mXi = casc->MassXi();
   Double_t mOmega = casc->MassOmega();

   //good Xi/Omega mass:
   //-------------------
   if((mXi>1.20&&mXi<1.45)||(mOmega>1.55&&mOmega<1.80)){

    //fill lambda
    fTLambdaPx[fTnCascade] = casc->MomV0X();
    fTLambdaPy[fTnCascade] = casc->MomV0Y();
    fTLambdaPz[fTnCascade] = casc->MomV0Z();
    fTLambdaPt[fTnCascade] = sqrt(pow(casc->MomV0X(),2)+pow(casc->MomV0Y(),2));
    fTLambdaDCA[fTnCascade] =   casc->DcaV0ToPrimVertex();
    fTLambdaDaughtersDCA[fTnCascade] = casc->DcaV0Daughters();;
    fTLambdaK0Mass[fTnCascade]= casc->MassK0Short();
    if(casc->Charge()>0) {
     fTLambdaMass[fTnCascade]=casc->MassAntiLambda();
    } else {
     fTLambdaMass[fTnCascade]=casc->MassLambda();
    }
    fTLambdaVr[fTnCascade] = sqrt(pow(casc->DecayVertexV0X(),2)+pow(casc->DecayVertexV0Y(),2));
    fTLambdaVx[fTnCascade] = casc->DecayVertexV0X();
    fTLambdaVy[fTnCascade] = casc->DecayVertexV0Y();
    fTLambdaVz[fTnCascade] = casc->DecayVertexV0Z();
    fTLambdaPA[fTnCascade] = casc->CosPointingAngle(evt->GetPrimaryVertex());

    //fill cascade
    fTCascadeCharge[fTnCascade] = casc->ChargeXi();
    fTCascadeDaughtersDCA[fTnCascade] = casc->DcaXiDaughters();
    fTCascadeVr[fTnCascade] = sqrt(pow(casc->DecayVertexXiX(),2)+pow(casc->DecayVertexXiY(),2));
    fTCascadeVx[fTnCascade] = casc->DecayVertexXiX();
    fTCascadeVy[fTnCascade] = casc->DecayVertexXiY();
    fTCascadeVz[fTnCascade] = casc->DecayVertexXiZ();
    //calc cascade dca, following recipe of AOD
    fTCascadeDCA[fTnCascade] = casc->DcaXiToPrimVertex(PrimVtx[0], PrimVtx[1], PrimVtx[2]);
    fTCascadePA[fTnCascade] = casc->CosPointingAngleXi(PrimVtx[0], PrimVtx[1], PrimVtx[2]);

    Filled = kTRUE;

   }//good Xi/Omega mass
  } //tracks filled
 return Filled;
}



//_______________________________________________________________________________________________________________________________
Bool_t AliOtonOmegaAnalysis::FillTreeTrackAOD(Int_t jj, AliAODcascade *casc) {
 //here jj=0 refers to pos, 1 neg, 2 bach
 Bool_t Filled = kFALSE;
 //get AOD track
 AliAODTrack *track;
 if(jj<2) {
  track = dynamic_cast<AliAODTrack*>(casc->GetDaughter(jj));
 }else{
  track = dynamic_cast<AliAODTrack*>(casc->GetDecayVertexXi()->GetDaughter(0));
 }

 //fill tpc inner momentum
 fTTrackTPCmom[fTnCascade][jj]= track->GetTPCmomentum();

 //REQUIRE EVERYTHING TO HAVE TPCMOM > 50MeV
 if(track->GetTPCmomentum()>.050){

  //fill dca
  fTTrackDCA[fTnCascade][jj]= track->DCA();

  //fill momentum
  fTTrackPx[fTnCascade][jj]=track->Px();
  fTTrackPy[fTnCascade][jj]=track->Py();
  fTTrackPz[fTnCascade][jj]=track->Pz();
  fTTrackP[fTnCascade][jj]=sqrt(pow(track->Px(),2)+pow(track->Py(),2)+pow(track->Pz(),2));

  if(
   fabs(fTTrackPx[fTnCascade][jj])<1000.
   && fabs(fTTrackPy[fTnCascade][jj])<1000.
   && fabs(fTTrackPz[fTnCascade][jj])<1000.
   ){

    //fill charge:
    fTTrackCharge[fTnCascade][jj]=track->Charge();

    //fill #TPCclusters
    fTTrackNcl[fTnCascade][jj]  = track->GetTPCNcls();

    //fill TPC CrF : Crossed rows / findable
    if(track->GetTPCNclsF()>0) fTTrackCrF[fTnCascade][jj]  = track->GetTPCClusterInfo(2,1)/track->GetTPCNclsF();


   //does it have a hit in SPD, SSD or TOF?
   if(track->HasPointOnITSLayer(0)||track->HasPointOnITSLayer(1)||track->HasPointOnITSLayer(4)||track->HasPointOnITSLayer(5)) fTTrackITStime[fTnCascade][jj] = kTRUE;
   if(track->GetTOFBunchCrossing() == 0) fTTrackTOFtime[fTnCascade][jj] = kTRUE;

   //-----------------------------------------------------------------------------------------
   //-----------------------------------------------------------------------------------------
   //-----------------------------------------------------------------------------------------
   //The question here is, is that enough that they have ITS or TOF time (in particular ITS),
   //or do I have to do SOMETHING with this time, i.e. check that is in
   //time with the vertex or with the other tracks of the cascade etc.
   //-----------------------------------------------------------------------------------------
   //-----------------------------------------------------------------------------------------
   //-----------------------------------------------------------------------------------------

   //TPCtrack only?
   if(track->GetStatus()&AliAODTrack::kTPCrefit && track->GetStatus()&AliAODTrack::kTPCin && track->GetStatus()&AliAODTrack::kITSrefit && track->GetStatus()&AliAODTrack::kITSin) fTTrackGLOBAL[fTnCascade][jj]=kTRUE;
   if(track->GetStatus()&AliAODTrack::kTPCrefit && track->GetStatus()&AliAODTrack::kTPCin && !(track->GetStatus()&AliAODTrack::kITSrefit) && !(track->GetStatus()&AliAODTrack::kITSin)) fTTrackTPConly[fTnCascade][jj]=kTRUE;
   if(!(track->GetStatus()&AliAODTrack::kTPCin) && track->GetStatus()&AliAODTrack::kITSrefit && track->GetStatus()&AliAODTrack::kITSin && !(track->GetStatus()&AliAODTrack::kITSpureSA)) fTTrackITScomplementary[fTnCascade][jj]=kTRUE;
   if(track->GetStatus()&AliAODTrack::kITSrefit && track->GetStatus()&AliAODTrack::kITSpureSA) fTTrackITSpure[fTnCascade][jj]=kTRUE;


  //Info not included so far, to be considered: - number of its clusters,  - number of tpc shared clusters,  - its shared clusters,  - track index
  //STILL TO BE CHECKED: - is the dca, momentum, etc evaluated for each track in the correct way?
  //                     - For tpc only tracks, using or not the vertex, etc...
  //                     - The info here stored coincides with the info evaluated by the vertexer?

   //pid
   AliPIDResponse::EDetPidStatus statusPID = fPIDResponse->CheckPIDStatus(AliPIDResponse::kITS, track);
   if (statusPID == AliPIDResponse::kDetPidOk)  {
     fTTrackITSspi[fTnCascade][jj] = fPIDResponse->NumberOfSigmas(AliPIDResponse::kITS, track, AliPID::kPion);
     fTTrackITSsk[fTnCascade][jj] = fPIDResponse->NumberOfSigmas(AliPIDResponse::kITS, track, AliPID::kKaon);
     fTTrackITSsp[fTnCascade][jj] = fPIDResponse->NumberOfSigmas(AliPIDResponse::kITS, track, AliPID::kProton);
    }
   statusPID = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTPC, track);
   if (statusPID == AliPIDResponse::kDetPidOk) {
     fTTrackTPCspi[fTnCascade][jj] = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, AliPID::kPion);
     fTTrackTPCsk[fTnCascade][jj] = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, AliPID::kKaon);
     fTTrackTPCsp[fTnCascade][jj] = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTPC, track, AliPID::kProton);
    }
   statusPID = fPIDResponse->CheckPIDStatus(AliPIDResponse::kTOF, track);
   if (statusPID == AliPIDResponse::kDetPidOk)  {
     fTTrackTOFspi[fTnCascade][jj] = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF, track, AliPID::kPion);
     fTTrackTOFsk[fTnCascade][jj] = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF, track, AliPID::kKaon);
     fTTrackTOFsp[fTnCascade][jj] = fPIDResponse->NumberOfSigmas(AliPIDResponse::kTOF, track, AliPID::kProton);
    }


  Filled = kTRUE;

 }//good momentum
}//tpcmom>50 for all 

 return Filled;
}




//_______________________________________________________________________________________________________________________________
Bool_t AliOtonOmegaAnalysis::FillProtonTrackAOD(AliAODTrack *Proton) {
 Bool_t Filled = kFALSE;

 //fill tpc inner momentum
 fTProtonTPCmom[fTnCascade]= Proton->GetTPCmomentum();

 //REQUIRE EVERYTHING TO HAVE TPCMOM > 50MeV
 if(Proton->GetTPCmomentum()>.050){

  //fill dca
  fTProtonDCA[fTnCascade]= Proton->DCA();

  //fill momentum
  fTProtonPx[fTnCascade]=Proton->Px();
  fTProtonPy[fTnCascade]=Proton->Py();
  fTProtonPz[fTnCascade]=Proton->Pz();
  fTProtonP[fTnCascade]=sqrt(pow(Proton->Px(),2)+pow(Proton->Py(),2)+pow(Proton->Pz(),2));
  fTProtonPt[fTnCascade]=sqrt(pow(Proton->Px(),2)+pow(Proton->Py(),2));

  if(
   fabs(fTProtonPx[fTnCascade])<1000.
   && fabs(fTProtonPy[fTnCascade])<1000.
   && fabs(fTProtonPz[fTnCascade])<1000.
   ){

    //fill charge:
    fTProtonCharge[fTnCascade]=Proton->Charge();

    //fill #TPCclusters
    fTProtonNcl[fTnCascade]  = Proton->GetTPCNcls();

    //fill TPC CrF : Crossed rows / findable
    if(Proton->GetTPCNclsF()>0) fTProtonCrF[fTnCascade]  = Proton->GetTPCClusterInfo(2,1)/Proton->GetTPCNclsF();


   //does it have a hit in SPD, SSD or TOF?
   if(Proton->HasPointOnITSLayer(0)||Proton->HasPointOnITSLayer(1)||Proton->HasPointOnITSLayer(4)||Proton->HasPointOnITSLayer(5)) fTProtonITStime[fTnCascade] = kTRUE;
   if(Proton->GetTOFBunchCrossing() == 0) fTProtonTOFtime[fTnCascade] = kTRUE;

   //TPCProton only?
   if(Proton->GetStatus()&AliAODTrack::kTPCrefit && Proton->GetStatus()&AliAODTrack::kTPCin && Proton->GetStatus()&AliAODTrack::kITSrefit && Proton->GetStatus()&AliAODTrack::kITSin) fTProtonGLOBAL[fTnCascade]=kTRUE;
   if(Proton->GetStatus()&AliAODTrack::kTPCrefit && Proton->GetStatus()&AliAODTrack::kTPCin && !(Proton->GetStatus()&AliAODTrack::kITSrefit) && !(Proton->GetStatus()&AliAODTrack::kITSin)) fTProtonTPConly[fTnCascade]=kTRUE;
   if(!(Proton->GetStatus()&AliAODTrack::kTPCin) && Proton->GetStatus()&AliAODTrack::kITSrefit && Proton->GetStatus()&AliAODTrack::kITSin && !(Proton->GetStatus()&AliAODTrack::kITSpureSA)) fTProtonITScomplementary[fTnCascade]=kTRUE;
   if(Proton->GetStatus()&AliAODTrack::kITSrefit && Proton->GetStatus()&AliAODTrack::kITSpureSA) fTProtonITSpure[fTnCascade]=kTRUE;

  Filled = kTRUE;

 }//good momentum
}//tpcmom>50 for all 

 return Filled;
}

