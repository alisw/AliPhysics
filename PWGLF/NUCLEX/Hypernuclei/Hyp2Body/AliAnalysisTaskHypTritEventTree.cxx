#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "THnSparse.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TRandom2.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"
#include "AliInputEventHandler.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliESDtrack.h"
#include "AliESDpid.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliESDv0.h"
#include "AliESDVertex.h"
#include "AliVertexerTracks.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCVertex.h"
#include "AliStack.h"
#include "TPDGCode.h"
#include "AliEventCuts.h"
#include "AliAnalysisTaskHypTritEventTree.h"
#include "AliReducedHypTritEvent.h"

using namespace std;
/// \cond CLASSIMP
ClassImp(AliAnalysisTaskHypTritEventTree)
/// \endcond

// Default Constructor
AliAnalysisTaskHypTritEventTree::AliAnalysisTaskHypTritEventTree()
  :AliAnalysisTaskSE("AliAnalysisTaskHypTritEventTree"),
  fInputHandler(0),
  fPID(0),
  fESDevent(0),
  fReducedEvent(0),
  fReducedEventMCGen(0),
  fStack(),
  fV0(),
  fV0Array(),
  fHistdEdx(0),
  fHistdEdxV0(0),
  fHistNumEvents(0),
  fHistTrigger(0),
  fHistV0(0),
  fHistMcGen(0),
  fHistMcRec(0),
  fTree(0),
  fTreeMCGen(0),
  fMCGenRecArray(),
  fHistogramList(NULL),
  fMCGenRec(),
  fMomPos(),
  fMomNeg(),
  fPrimaryVertex(),
  fMagneticField(),
  fNV0Cand(),
  fMcGenRecCounter(),
  fPidQa(0),
  fUseAnalysisTrkSel(kTRUE),
  fPIDCheckOnly(kFALSE),
  fMCtrue(0),
  fEventCuts(),
  fPeriod(0),
  fTriggerMask(0),
  fBetheSplines(kFALSE),
  fBetheParamsHe(),
  fBetheParamsT() {

  }

// Constructor
AliAnalysisTaskHypTritEventTree::AliAnalysisTaskHypTritEventTree(const char *name)
  :AliAnalysisTaskSE(name),
  fInputHandler(0),
  fPID(0),
  fESDevent(0),
  fReducedEvent(0),
  fReducedEventMCGen(0),
  fStack(),
  fV0(),
  fV0Array(),
  fHistdEdx(0),
  fHistdEdxV0(0),
  fHistNumEvents(0),
  fHistTrigger(0),
  fHistV0(0),
  fHistMcGen(0),
  fHistMcRec(0),
  fTree(0),
  fTreeMCGen(0),
  fMCGenRecArray(),
  fHistogramList(NULL),
  fMCGenRec(),
  fMomPos(),
  fMomNeg(),
  fPrimaryVertex(),
  fMagneticField(),
  fNV0Cand(),
  fMcGenRecCounter(),
  fPidQa(0),
  fUseAnalysisTrkSel(kTRUE),
  fPIDCheckOnly(kFALSE),
  fMCtrue(0),
  fEventCuts(),
  fPeriod(00),
  fTriggerMask(),
  fBetheSplines(kFALSE),
  fBetheParamsHe(),
  fBetheParamsT()
  {
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TTree::Class());
    DefineOutput(3, TTree::Class());
  }

// Destructor
AliAnalysisTaskHypTritEventTree::~AliAnalysisTaskHypTritEventTree() {
  if (fMCGenRecArray) fMCGenRecArray->Delete();
}

void AliAnalysisTaskHypTritEventTree::UserCreateOutputObjects() {
  fInputHandler = dynamic_cast<AliESDInputHandler*>
    (AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  if(!fInputHandler) {
    AliError("Could not get ESD InputHandler.\n");
    return;
  }
  fPID = fInputHandler->GetESDpid();
  if (!fPID) {
    AliError("Could not get PID response.\n");
    return;
  }
  fHistdEdx = new TH2F("fHistdEdX","dE/dx;#frac{#it{p}}{z} (GeV/#it{c});TPC Signal (a.u.)",1000,-5.0,5.0,1000,0.0,1500);
  fHistdEdxV0 = new TH2F("fHistdEdXV0","dE/dx;#frac{#it{p}}{z} (GeV/#it{c});TPC Signal (a.u.)",1000,-5.0,5.0,1000,0.0,1500);

  fHistNumEvents = new TH1F("fHistNumEvents","Number of Events",2,0,2);
  fHistNumEvents->GetXaxis()->SetBinLabel(1,"before PhysSel");
  fHistNumEvents->GetXaxis()->SetBinLabel(2,"after PhysSel");

  fHistTrigger = new TH1F("fHistTrigger","Trigger",7,0,7);
  fHistTrigger->GetXaxis()->SetBinLabel(1,"other");
  fHistTrigger->GetXaxis()->SetBinLabel(2,"kINT7");
  fHistTrigger->GetXaxis()->SetBinLabel(3,"kHighMultV0");
  fHistTrigger->GetXaxis()->SetBinLabel(4,"kHighMultSPD");
  fHistTrigger->GetXaxis()->SetBinLabel(5,"HNU");
  fHistTrigger->GetXaxis()->SetBinLabel(6,"HQU");
  fHistTrigger->GetXaxis()->SetBinLabel(7,"HJT");
  fHistV0 = new TH1F("fHistV0","Trigger V0s",7,0,7);
  fHistV0->GetXaxis()->SetBinLabel(1,"other");
  fHistV0->GetXaxis()->SetBinLabel(2,"kINT7");
  fHistV0->GetXaxis()->SetBinLabel(3,"kHighMultV0");
  fHistV0->GetXaxis()->SetBinLabel(4,"kHighMultSPD");
  fHistV0->GetXaxis()->SetBinLabel(5,"HNU");
  fHistV0->GetXaxis()->SetBinLabel(6,"HQU");
  fHistV0->GetXaxis()->SetBinLabel(7,"HJT");

  fHistMcGen = new TH1F("fHistMcGen","mc generated; ct (cm);counts",40,0,40);
  fHistMcRec = new TH1F("fHistMcRec","mc reconstructed; ct (cm);counts",40,0,40);

//  TF1 *tBethe = new TF1("tBethe","[7]*[5]*[5]*AliExternalTrackParam::BetheBlochAleph([5]*x/([6]),[0],[1],[2],[3],[4])",0,25);
//  tBethe->SetParameters(fBetheParamsT[0], fBetheParamsT[1],
//    fBetheParamsT[2],fBetheParamsT[3],
//    fBetheParamsT[4], 1, AliPID::ParticleMass(AliPID::kTriton), 1);
//  TF1 *tBetheU = new TF1("tBetheU","[7]*[5]*[5]*AliExternalTrackParam::BetheBlochAleph([5]*x/([6]),[0],[1],[2],[3],[4])",0,25);
//  tBetheU->SetParameters(fBetheParamsT[0], fBetheParamsT[1],
//    fBetheParamsT[2],fBetheParamsT[3],
//    fBetheParamsT[4], 1, AliPID::ParticleMass(AliPID::kTriton), 1+3*fBetheParamsT[5]);
//  TF1 *tBetheL = new TF1("tBetheL","[7]*[5]*[5]*AliExternalTrackParam::BetheBlochAleph([5]*x/([6]),[0],[1],[2],[3],[4])",0,25);
//  tBetheL->SetParameters(fBetheParamsT[0], fBetheParamsT[1],
//    fBetheParamsT[2],fBetheParamsT[3],
//    fBetheParamsT[4], 1, AliPID::ParticleMass(AliPID::kTriton), 1-3*fBetheParamsT[5]);
//
//  TF1 *bethe = new TF1("heBethe","[7]*[5]*[5]*AliExternalTrackParam::BetheBlochAleph([5]*x/([6]),[0],[1],[2],[3],[4])",0,25);
//  bethe->SetParameters(fBetheParamsHe[0], fBetheParamsHe[1],
//    fBetheParamsHe[2],fBetheParamsHe[3],
//    fBetheParamsHe[4], 2, AliPID::ParticleMass(AliPID::kHe3), 1);
//  TF1 *betheU = new TF1("heBetheU","[7]*[5]*[5]*AliExternalTrackParam::BetheBlochAleph([5]*x/([6]),[0],[1],[2],[3],[4])",0,25);
//  betheU->SetParameters(fBetheParamsHe[0], fBetheParamsHe[1],
//    fBetheParamsHe[2],fBetheParamsHe[3],
//    fBetheParamsHe[4], 2, AliPID::ParticleMass(AliPID::kHe3), 1+3*fBetheParamsHe[5]);
//  TF1 *betheL = new TF1("heBetheL","[7]*[5]*[5]*AliExternalTrackParam::BetheBlochAleph([5]*x/([6]),[0],[1],[2],[3],[4])",0,25);
//  betheL->SetParameters(fBetheParamsHe[0], fBetheParamsHe[1],
//    fBetheParamsHe[2],fBetheParamsHe[3],
//    fBetheParamsHe[4], 2, AliPID::ParticleMass(AliPID::kHe3), 1-3*fBetheParamsHe[5]);

  fHistogramList = new TList();
  fHistogramList->SetOwner(kTRUE);
  fHistogramList->SetName(GetName());
  fHistogramList->Add(fHistdEdx);
  fHistogramList->Add(fHistdEdxV0);
  fHistogramList->Add(fHistNumEvents);
  fHistogramList->Add(fHistTrigger);
  fHistogramList->Add(fHistV0);
  fHistogramList->Add(fHistMcGen);
  fHistogramList->Add(fHistMcRec);
  //fHistogramList->Add(tBethe);
  //fHistogramList->Add(tBetheU);
  //fHistogramList->Add(tBetheL);
  //fHistogramList->Add(bethe);
  //fHistogramList->Add(betheU);
  //fHistogramList->Add(betheL);

  fEventCuts.AddQAplotsToList(fHistogramList);

  fTree = new TTree("tree","fTree");
  fReducedEvent = new AliReducedHypTritEvent();
  fTree->Branch("event","AliReducedHypTritEvent",&fReducedEvent,32000,99);
  fTreeMCGen = new TTree("tree_mc", "fTreeMCGen");
  fReducedEventMCGen = new AliReducedHypTritEvent();
  fTreeMCGen->Branch("event","AliReducedHypTritEvent",&fReducedEventMCGen,32000,99);
  PostData(1, fHistogramList);
  PostData(2, fTree);
  PostData(3, fTreeMCGen);
}

void AliAnalysisTaskHypTritEventTree::UserExec(Option_t *) {
  // MC
  fMCtrue = kTRUE;
  AliMCEventHandler* mcEventHandler = dynamic_cast<AliMCEventHandler*>
      (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!mcEventHandler) {
    fMCtrue = kFALSE;
  }
  AliMCEvent* mcEvent = 0x0;
  if (mcEventHandler) mcEvent = mcEventHandler->MCEvent();
  if (!mcEvent) {
    if (fMCtrue) return;
  }
    if (fMCtrue) {
    fStack = mcEvent->Stack();
    if (!fStack) return;
  }
  // Data

  fESDevent = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESDevent) {
    AliError("Could not get ESD Event.\n");
    return;
  }
  if (!fPID) {
    AliError("Could not get PID response.\n");
    return;
  }
  for (int i = 0; i < 40; i++) {
    fMCGenRec[i] = -1;
  }
  fHistNumEvents->Fill(0);
  Float_t centrality = -1;
  const AliESDVertex *vertex = fESDevent->GetPrimaryVertexTracks();
  fEventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kINT7 | AliVEvent::kTRD | AliVEvent::kHighMultV0 | AliVEvent::kHighMultSPD);

  if(!fEventCuts.AcceptEvent(fESDevent)) {
    PostData(1,fHistogramList);
    return;
  }

  Int_t runNumber = fESDevent->GetRunNumber();
  SetBetheBlochParams(runNumber);
  TriggerSelection();
  SetMultiplicity();
  fHistNumEvents->Fill(1);
  fReducedEvent->fCentrality = centrality;
  fMagneticField  = fESDevent->GetMagneticField();
  fPrimaryVertex.SetXYZ(vertex->GetX(),vertex->GetY(),vertex->GetZ());
  fReducedEvent->fVertexPosition = fPrimaryVertex;
  fReducedEvent->fRunNumber = runNumber;
  fMcGenRecCounter = 0;
  fV0Array = fReducedEvent->fV0s;
  fNV0Cand = 0;
  fMCGenRecArray = new TObjArray();
  AliESDtrackCuts trackCutsV0("AlitrackCutsV0", "AlitrackCutsV0");

  if(fUseAnalysisTrkSel){
    trackCutsV0.SetEtaRange(-0.9,0.9);
    trackCutsV0.SetAcceptKinkDaughters(kFALSE);
    trackCutsV0.SetRequireTPCRefit(kTRUE);
    trackCutsV0.SetMaxChi2PerClusterTPC(5);
    trackCutsV0.SetMinNClustersTPC(60);
  } else {
      trackCutsV0.SetAcceptKinkDaughters(kFALSE);
      trackCutsV0.SetMinNClustersTPC(80);
      trackCutsV0.SetMaxChi2PerClusterITS(10);// TO BE INVESTIGATED !!!!!!!!!!!!!!
      trackCutsV0.SetMaxChi2PerClusterTPC(5);
      trackCutsV0.SetRequireTPCRefit(kTRUE);
      trackCutsV0.SetRequireITSRefit(kTRUE);
      trackCutsV0.SetMinNClustersITS(2);
      trackCutsV0.SetMaxDCAToVertexXY(0.1);
      trackCutsV0.SetMaxDCAToVertexZ(0.5);
      trackCutsV0.SetEtaRange(-0.8,0.8);
  }
  // Pidqa loop
  if (fPidQa) {
    AliESDtrackCuts* trackCutsPid = new AliESDtrackCuts("trackCutsPid", "trackCutsPid");
    trackCutsPid = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts ();
    trackCutsPid->SetEtaRange(-0.9,0.9);
    for (Int_t itrack = 0; itrack < fESDevent->GetNumberOfTracks(); itrack++) {
      AliESDtrack* track = fESDevent->GetTrack(itrack);
      if (!trackCutsPid->AcceptTrack(track)) continue;
      Double_t momentum = track->GetInnerParam()->GetP();
      fHistdEdx->Fill(momentum * track->GetSign(), track->GetTPCsignal());
    }
    delete trackCutsPid;
  }

  // V0 loop
  for (Int_t ivertex = 0; ivertex < fESDevent->GetNumberOfV0s(); ivertex++) {
	fHistV0->Fill(fReducedEvent->fTrigger);
    fV0 = fESDevent->GetV0(ivertex);
    Bool_t v0ChargeCorrect = kTRUE;
    AliESDtrack* trackN = fESDevent->GetTrack(fV0->GetIndex(0));
    AliESDtrack* trackP = fESDevent->GetTrack(fV0->GetIndex(1));
    // Checks charge because of bug in V0 interface.
    if (trackN->GetSign() > 0 ) {
      trackN = fESDevent->GetTrack(fV0->GetIndex(1));
      trackP = fESDevent->GetTrack(fV0->GetIndex(0));
      v0ChargeCorrect = kFALSE;
    }
    if (!trackCutsV0.AcceptTrack(trackN)) continue;
    if (!trackCutsV0.AcceptTrack(trackP)) continue;
    fHistdEdxV0->Fill(trackP->GetInnerParam()->GetP() * trackP->GetSign(), trackP->GetTPCsignal());
    fHistdEdxV0->Fill(trackN->GetInnerParam()->GetP() * trackN->GetSign(), trackN->GetTPCsignal());
    if(fPIDCheckOnly) continue;
    if (trackN->GetTPCsignal() > 1500 || trackP->GetTPCsignal() > 1500) continue;
    if (trackN->GetInnerParam()->GetP() > 5 || trackP->GetInnerParam()->GetP() > 5) continue;
    Bool_t pionPositive     = kFALSE;
    Bool_t pionNegative     = kFALSE;
    // Bool_t protonPositive   = kFALSE;
    // Bool_t protonNegative   = kFALSE;
    // Bool_t deuteronPositive = kFALSE;
    // Bool_t deuteronNegative = kFALSE;
    // Bool_t tritonPositive   = kFALSE;
    // Bool_t tritonNegative   = kFALSE;
    Bool_t helium3Positive  = kFALSE;
    Bool_t helium3Negative  = kFALSE;
    if (TMath::Abs(fPID->NumberOfSigmasTPC(trackP, AliPID::kPion)) < 3) {
      pionPositive = kTRUE;
    }
    if (TMath::Abs(fPID->NumberOfSigmasTPC(trackN, AliPID::kPion)) < 3) {
      pionNegative = kTRUE;
    }
    if (fBetheSplines) {
//      if (TMath::Abs(fPID->NumberOfSigmasTPC(trackP, AliPID::kTriton)) < 3) {
//        tritonPositive = kTRUE;
//      } else if (TMath::Abs(fPID->NumberOfSigmasTPC(trackN, AliPID::kTriton)) < 3) {
//        tritonNegative = kTRUE;
//      }
      if (TMath::Abs(fPID->NumberOfSigmasTPC(trackP, AliPID::kHe3)) < 4) {
        helium3Positive = kTRUE;
      } else if (TMath::Abs(fPID->NumberOfSigmasTPC(trackN, AliPID::kHe3)) < 4) {
         helium3Negative = kTRUE;
      }
    } else {
      if (TMath::Abs(Bethe(*trackP, AliPID::ParticleMass(AliPID::kHe3),  2, fBetheParamsHe)) < 4) {
        helium3Positive = kTRUE;
      } else if (TMath::Abs(Bethe(*trackN, AliPID::ParticleMass(AliPID::kHe3), 2, fBetheParamsHe)) < 4) {
        helium3Negative = kTRUE;
      }
//      if (TMath::Abs(Bethe(*trackP, AliPID::ParticleMass(AliPID::kTriton), 1, fBetheParamsT)) < 3) {
//        tritonPositive = kTRUE;
//      } else if (TMath::Abs(Bethe(*trackN, AliPID::ParticleMass(AliPID::kTriton), 1, fBetheParamsT)) < 3) {
//        tritonNegative = kTRUE;
//      }
    }
    if (helium3Positive && pionNegative) {
      SetMomentum( 2, v0ChargeCorrect);
      CalculateV0(*trackN, *trackP,  AliPID::kPion, AliPID::kHe3);
    }
    if (helium3Negative && pionPositive) {
      SetMomentum(-2, v0ChargeCorrect);
      CalculateV0(*trackN, *trackP, AliPID::kHe3, AliPID::kPion);
    }
  }
  fReducedEvent->fNumberV0s = (fNV0Cand);
  if (fNV0Cand) fTree->Fill();
  fReducedEvent->ClearEvent();
  if (fMCtrue) {
    const AliMCVertex* mcVertex = (const AliMCVertex*) mcEvent->GetPrimaryVertex();
    fPrimaryVertex.SetXYZ(mcVertex->GetX(), mcVertex->GetY(), mcVertex->GetZ());
    MCStackLoop(fStack);
  }
  PostData(1, fHistogramList);
  PostData(2, fTree);
  PostData(3, fTreeMCGen);
}
//_____________________________________________________________________________
void AliAnalysisTaskHypTritEventTree::Terminate(const Option_t*) {
  if (!GetOutputData(0)) return;
}
//_____________________________________________________________________________
/// Sets momentum of a pair of identified tracks with a given charge
/// \param charge charge of nuclei
/// \param v0Charge corrects for mistake in AliESDv0. It is true when the charge was correctly set
void AliAnalysisTaskHypTritEventTree::SetMomentum(Int_t charge, Bool_t v0ChargeCorrect) {
  TVector3 momentumVector(0,0,0);
  if (charge > 0) {
    fV0->GetPPxPyPz(momentumVector(0), momentumVector(1), momentumVector(2));
    fMomPos.SetVect(charge * momentumVector);
    fV0->GetNPxPyPz(momentumVector(0), momentumVector(1), momentumVector(2));
    fMomNeg.SetVect(momentumVector);
    if (!v0ChargeCorrect) {
      fV0->GetNPxPyPz(momentumVector(0), momentumVector(1), momentumVector(2));
      fMomPos.SetVect(charge * momentumVector);
      fV0->GetPPxPyPz(momentumVector(0), momentumVector(1), momentumVector(2));
      fMomNeg.SetVect(momentumVector);
    }
  } else {
    fV0->GetPPxPyPz(momentumVector(0), momentumVector(1), momentumVector(2));
    fMomPos.SetVect(momentumVector);
    fV0->GetNPxPyPz(momentumVector(0), momentumVector(1), momentumVector(2));
    fMomNeg.SetVect(-charge * momentumVector);
    if (!v0ChargeCorrect) {
      fV0->GetNPxPyPz(momentumVector(0), momentumVector(1), momentumVector(2));
      fMomPos.SetVect(momentumVector);
      fV0->GetPPxPyPz(momentumVector(0), momentumVector(1), momentumVector(2));
      fMomNeg.SetVect(-charge * momentumVector);
    }
  }
}
//_____________________________________________________________________________
/// Calculates V0 parameters and sets them in the reduced event
/// \param trackN  track with negative charge
/// \param trackP  track with positive charge
/// \param typeNeg Particle type hypothesis of negative particle used in calculation
/// \param typePos Particle type hypothesis of positive particle used in calculation
void AliAnalysisTaskHypTritEventTree::CalculateV0(const AliESDtrack& trackN, const AliESDtrack& trackP, AliPID::EParticleType typeNeg, AliPID::EParticleType typePos) {
  fMomPos.SetE(TMath::Sqrt(AliPID::ParticleMass(typePos) * AliPID::ParticleMass(typePos) + fMomPos.Vect().Mag2()));
  fMomNeg.SetE(TMath::Sqrt(AliPID::ParticleMass(typeNeg) * AliPID::ParticleMass(typeNeg) + fMomNeg.Vect().Mag2()));
  Double_t v0M  = (fMomNeg + fMomPos).M();
  Double_t v0Pt = (fMomNeg + fMomPos).Pt();
  Double_t v0P  = (fMomNeg + fMomPos).P();
  Double_t v0Y  = (fMomNeg + fMomPos).Rapidity();
  Int_t charge = -AliPID::ParticleCharge(typeNeg) + AliPID::ParticleCharge(typePos);
  Float_t dcav0 = fV0->GetDcaV0Daughters();
  Float_t cosineOfPointingAngle = fV0->GetV0CosineOfPointingAngle();
  TVector3 secondaryVertex(fV0->Xv(), fV0->Yv(), fV0->Zv());
  secondaryVertex = secondaryVertex - fPrimaryVertex;
  AliReducedHypTritV0* reducedV0 = dynamic_cast<AliReducedHypTritV0*>(fV0Array->ConstructedAt(fNV0Cand));
  fNV0Cand = fNV0Cand + 1;
  AliReducedHypTritTrack* reducedPi = reducedV0->Pi();
  AliReducedHypTritTrack* reducedHe = reducedV0->He();
  reducedV0->fPosition = secondaryVertex;
  TVector3 momentumVector(0,0,0);
  fV0->GetPxPyPz(momentumVector(0), momentumVector(1), momentumVector(2));
  reducedV0->fPvect = momentumVector;
  reducedV0->fCharge = charge;
  reducedV0->fM = v0M;
  reducedV0->fP= v0P;
  reducedV0->fPt = v0Pt;
  reducedV0->fDcaV0 = dcav0;
  reducedV0->fCosPointingAngle = cosineOfPointingAngle;
  reducedV0->fDecayLength = secondaryVertex.Mag() * v0M / v0P;
  reducedV0->fMcTruth = 0;
  reducedV0->fRapidity = v0Y;
  reducedV0->fParticleSpecies = typeNeg * 100 + typePos;
  reducedV0->fOnFlyStatus = fV0->GetOnFlyStatus();
  if (charge < 0) {
    reducedHe->fP = fMomNeg;
    reducedPi->fP = fMomPos;
    reducedHe->fDedx = trackN.GetTPCsignal();
    reducedPi->fDedx = trackP.GetTPCsignal();
    if (fBetheSplines) {
      reducedHe->fDedxSigma = fPID->NumberOfSigmasTPC(&trackN, typeNeg);
      reducedHe->fDedxSigmaTriton = fPID->NumberOfSigmasTPC(&trackN, AliPID::kTriton);
    } else {
      reducedHe->fDedxSigma = Bethe(trackN, AliPID::ParticleMass(typeNeg), 2, fBetheParamsHe);
      reducedHe->fDedxSigmaTriton = Bethe(trackN, AliPID::ParticleMass(AliPID::kTriton), 1, fBetheParamsT);
    }
    reducedPi->fDedxSigma = fPID->NumberOfSigmasTPC(&trackP, typePos);
    reducedHe->fDca = TMath::Abs(trackN.GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
    reducedPi->fDca = TMath::Abs(trackP.GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
    reducedPi->fEta = trackP.Eta();
    reducedHe->fEta = trackN.Eta();
    reducedPi->fPhi = trackP.Phi();
    reducedHe->fPhi = trackN.Phi();
    reducedPi->fGeoLength = GeoLength(trackP);
    reducedHe->fGeoLength = GeoLength(trackN);
    reducedHe->fPtrack = trackN.GetInnerParam()->GetP();
    reducedPi->fPtrack = trackP.GetInnerParam()->GetP();
    reducedHe->fTpcNClusters = trackN.GetTPCNcls();
    reducedPi->fTpcNClusters = trackP.GetTPCNcls();
  }
  if (charge > 0) {
    reducedHe->fP = fMomPos;
    reducedPi->fP = fMomNeg;
    reducedHe->fDedx = trackP.GetTPCsignal();
    reducedPi->fDedx = trackN.GetTPCsignal();
    if (fBetheSplines) {
      reducedHe->fDedxSigma = fPID->NumberOfSigmasTPC(&trackP, typePos);
      reducedHe->fDedxSigmaTriton = fPID->NumberOfSigmasTPC(&trackP, AliPID::kTriton);
    } else {
      reducedHe->fDedxSigma = Bethe(trackP, AliPID::ParticleMass(typePos), 2, fBetheParamsHe);
      reducedHe->fDedxSigmaTriton = Bethe(trackP, AliPID::ParticleMass(AliPID::kTriton), 1, fBetheParamsT);
    }
    reducedPi->fDedxSigma = fPID->NumberOfSigmasTPC(&trackP, typePos);
    reducedHe->fDca = TMath::Abs(trackN.GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
    reducedPi->fDca = TMath::Abs(trackP.GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
    reducedPi->fEta = trackP.Eta();
    reducedPi->fDedxSigma = fPID->NumberOfSigmasTPC(&trackN, typeNeg);
    reducedHe->fDca = TMath::Abs(trackP.GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
    reducedPi->fDca = TMath::Abs(trackN.GetD(fPrimaryVertex.X(), fPrimaryVertex.Y(), fMagneticField));
    reducedPi->fEta = trackN.Eta();
    reducedHe->fEta = trackP.Eta();
    reducedPi->fPhi = trackN.Phi();
    reducedHe->fPhi = trackP.Phi();
    reducedHe->fGeoLength = GeoLength(trackP);
    reducedPi->fGeoLength = GeoLength(trackN);
    reducedHe->fPtrack = trackP.GetInnerParam()->GetP();
    reducedPi->fPtrack = trackN.GetInnerParam()->GetP();
    reducedHe->fTpcNClusters = trackP.GetTPCNcls();
    reducedPi->fTpcNClusters = trackN.GetTPCNcls();
  }
  if (fMCtrue && ((typePos == AliPID::kHe3 && typeNeg == AliPID::kPion) || (typePos == AliPID::kPion && typeNeg == AliPID::kHe3))) {
    Int_t labelP = trackP.GetLabel();
    Int_t labelN = trackN.GetLabel();
    TParticle *daughterParticleP = fStack->Particle(TMath::Abs(labelP));
    TParticle *daughterParticleN = fStack->Particle(TMath::Abs(labelN));
    TParticle *particleMotherP = fStack->Particle(TMath::Abs(daughterParticleP->GetFirstMother()));
    TParticle *particleMotherN = fStack->Particle(TMath::Abs(daughterParticleN->GetFirstMother()));
    Int_t labelMotherP = daughterParticleP->GetFirstMother();
    Int_t labelMotherN = daughterParticleN->GetFirstMother();
    if (((particleMotherN->GetPdgCode() == 1010010030 && particleMotherP->GetPdgCode() == 1010010030)   ||
          (particleMotherN->GetPdgCode() == -1010010030 && particleMotherP->GetPdgCode() == -1010010030)) &&
        (labelMotherN == labelMotherP)) {
      reducedV0->fMcTruth = 1;
      fMCGenRecArray->AddAtAndExpand(reducedV0, fMcGenRecCounter);
      fMCGenRec[fMcGenRecCounter] = labelMotherP;
      fMcGenRecCounter++;
      //if (McCuts(*reducedV0, *reducedHe, *reducedPi)) 
		fHistMcRec->Fill(reducedV0->fDecayLength);
    }
  }
}
//_____________________________________________________________________________
/// Loops over MC stack and matches generated particles with reconstructed particles
/// \param stack MC stack
void AliAnalysisTaskHypTritEventTree::MCStackLoop(AliStack *stack) {
  TClonesArray *v0Array = (TClonesArray*) fReducedEventMCGen->fV0s;
  Int_t nV0Gen = 0;
  for (Int_t istack = 0; istack < stack->GetNtrack(); istack++) {
    const TParticle *tparticleMother = stack->Particle(istack);
    Long_t pdgCodeMother = tparticleMother->GetPdgCode();
    if (pdgCodeMother == 1010010030 || pdgCodeMother == -1010010030) {
      TLorentzVector momentumDaughter1;
      TLorentzVector momentumDaughter2;
      Int_t labelSecondDaughter = tparticleMother->GetDaughter(1);
      Int_t labelFirstDaughter  = labelSecondDaughter - 1;
      TParticle *daughterparticle1 = stack->Particle(TMath::Abs(labelFirstDaughter));
      TParticle *daughterparticle2 = stack->Particle(TMath::Abs(labelSecondDaughter));
      if ((daughterparticle1->GetPdgCode() == 1000020030  /*Helium3*/     &&
            daughterparticle2->GetPdgCode() == -211)        /*PionMinus*/   ||
          (daughterparticle1->GetPdgCode() == -1000020030 /*AntiHelium3*/ &&
           daughterparticle2->GetPdgCode() == 211))        /*PionPlus*/    {
        momentumDaughter1.SetPxPyPzE(daughterparticle1->Px(), daughterparticle1->Py(),
            daughterparticle1->Pz(), daughterparticle1->Energy());
        momentumDaughter2.SetPxPyPzE(daughterparticle2->Px(), daughterparticle2->Py(),
            daughterparticle2->Pz(), daughterparticle2->Energy());
        AliReducedHypTritV0 *reducedV0 = (AliReducedHypTritV0*)v0Array->ConstructedAt(nV0Gen);
        AliReducedHypTritTrack *reducedPi = reducedV0->Pi();
        AliReducedHypTritTrack *reducedHe = reducedV0->He();
        Double_t posx = daughterparticle1->Vx();
        Double_t posy = daughterparticle1->Vy();
        Double_t posz = daughterparticle1->Vz();
        Double_t disx = posx - fPrimaryVertex.X();
        Double_t disy = posy - fPrimaryVertex.Y();
        Double_t disz = posz - fPrimaryVertex.Z();
        Double_t distance = TMath::Sqrt(disx*disx + disy*disy + disz*disz );
        reducedV0->fM = tparticleMother->GetCalcMass();
        reducedV0->fP = tparticleMother->P();
        reducedV0->fPt = tparticleMother->Pt();
		reducedV0->fRapidity = (momentumDaughter1 + momentumDaughter2).Rapidity();
        reducedV0->fDecayLength = distance * tparticleMother->GetCalcMass() / tparticleMother->P();
        reducedV0->fMcTruth = 0;
		reducedV0->fCharge = -1;
		if (pdgCodeMother == 1010010030) 
			reducedV0->fCharge = 1;
        fHistMcGen->Fill(reducedV0->fDecayLength);
        nV0Gen = nV0Gen +1;
      }
		}
	}
	fReducedEventMCGen->fNumberV0s = nV0Gen;
	fTreeMCGen->Fill();
	fReducedEventMCGen->ClearEvent();
	fMCGenRecArray->Clear();
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskHypTritEventTree::TriggerSelection() {
  fReducedEvent->fTrigger = 0;
  if (!fMCtrue){
  TString classes = fESDevent->GetFiredTriggerClasses();
  fReducedEvent->fTriggerClasses = classes;
  if ((fInputHandler->IsEventSelected() & AliVEvent::kINT7)) fReducedEvent->fTrigger = 1;
  if ((fInputHandler->IsEventSelected() & AliVEvent::kHighMultV0)) fReducedEvent->fTrigger = 2;
  if ((fInputHandler->IsEventSelected() & AliVEvent::kHighMultSPD)) fReducedEvent->fTrigger = 3;
  if (classes.Contains("HNU")) fReducedEvent->fTrigger = 4;
  if (classes.Contains("HQU")) fReducedEvent->fTrigger = 5;
  if (classes.Contains("HJT")) fReducedEvent->fTrigger = 6;
  fHistTrigger->Fill(fReducedEvent->fTrigger);
  } else {
  	// MC: simulate TRD trigger
		Int_t nTrdTracks = fESDevent->GetNumberOfTrdTracks();
		if (nTrdTracks > 0) {
			for (Int_t iTrack = 0; iTrack < nTrdTracks; ++iTrack) {
				AliESDTrdTrack* trdTrack = fESDevent->GetTrdTrack(iTrack);
				if (!trdTrack) continue;
				// simulate HNU
				if((trdTrack->GetPID() >= 255 && trdTrack->GetNTracklets() == 4) || 
					(trdTrack->GetPID() >= 235 && trdTrack->GetNTracklets() > 4)) {
					fReducedEvent->fTrigger = 4;
				}
				// simulate HQU
				if (TMath::Abs(trdTrack->GetPt()) >= 256 &&
					trdTrack->GetPID() >= 130 && 
					trdTrack->GetNTracklets() >= 5 && 
					(trdTrack->GetLayerMask() & 1) ){	
						Float_t sag = GetInvPtDevFromBC(trdTrack->GetB(), trdTrack->GetC());
						if (sag < 0.2 && sag > -0.2) {
							fReducedEvent->fTrigger = 5;
						}
				}
			}
		}
	fHistTrigger->Fill(fReducedEvent->fTrigger);
	}	
	// additional information for high multiplicity trigger 
	AliESDVZERO *vzero = fESDevent->GetVZEROData();
	fReducedEvent->fV0Multiplicity = 0;
	for (Int_t ii = 0; ii < 64; ii++){
		fReducedEvent->fV0Multiplicity += vzero->GetMultiplicity(ii);
	}	
	AliMultiplicity *multSPD = fESDevent->GetMultiplicity();
	fReducedEvent->fSPDCluster	= multSPD->GetNumberOfSPDClusters();
	fReducedEvent->fSPDTracklets = multSPD->GetNumberOfTracklets();
	fReducedEvent->fSPDFiredChips0 = multSPD->GetNumberOfFiredChips(0);
	fReducedEvent->fSPDFiredChips1 = multSPD->GetNumberOfFiredChips(1);
  
  Bool_t isTriggered = kTRUE;
  if (fReducedEvent->fTrigger == 0) isTriggered = kFALSE;
  return isTriggered;
}
//_____________________________________________________________________________
Float_t AliAnalysisTaskHypTritEventTree::GetInvPtDevFromBC(Int_t b, Int_t c) {
	//returns d(1/Pt) in c/GeV 
	//in case of no gtu simulation -> return maximum 0.5
	if(b==0 && c==0) return 0.5;
	Int_t tmp = (((b & 0xfff) << 12) ^ 0x800000) - 0x800000;
	tmp += (c & 0xfff);
	Float_t invPtDev = tmp * 0.000001;
	return invPtDev;
}
//_____________________________________________________________________________
void AliAnalysisTaskHypTritEventTree::SetMultiplicity() {
	AliMultSelection *MultSelection = (AliMultSelection*) fESDevent->FindListObject("MultSelection");
	if (MultSelection) {
		fReducedEvent->fMultV0M = MultSelection->GetMultiplicityPercentile("V0M");
		fReducedEvent->fMultOfV0M = MultSelection->GetMultiplicityPercentile("OnlineV0M");
		fReducedEvent->fMultSPDTracklet = MultSelection->GetMultiplicityPercentile("SPDClusters");
		fReducedEvent->fMultSPDCluster = MultSelection->GetMultiplicityPercentile("SPDTracklets");
		fReducedEvent->fMultRef05 = MultSelection->GetMultiplicityPercentile("RefMult05");
		fReducedEvent->fMultRef08 = MultSelection->GetMultiplicityPercentile("RefMult08");
	}
}
//_____________________________________________________________________________
/// Calculates number of sigma deviation from expected dE/dx in TPC
/// \param track particle track
/// \param mass mass hypothesis of particle
/// \param charge particle charge hypothesis
/// \param params Parameters of Aleph parametrization of Bethe Energy-loss
Double_t AliAnalysisTaskHypTritEventTree::Bethe(const AliESDtrack& track, Double_t mass, Int_t charge, Double_t* params){
  Double_t expected = charge*charge*AliExternalTrackParam::BetheBlochAleph(charge*track.GetInnerParam()->GetP()/mass,params[0],params[1],params[2],params[3],params[4]);
  Double_t sigma = expected*params[5];
  if (TMath::IsNaN(expected)) return -999;
  return (track.GetTPCsignal() - expected) / sigma;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskHypTritEventTree::McCuts(const AliReducedHypTritV0& v0, const AliReducedHypTritTrack& he, const AliReducedHypTritTrack& pi) {
  Bool_t cut = v0.fDcaV0 < 0.03; // &&
  return cut;
//    he.fDca < 4 &&
//    he.fDedxSigmaTriton > 5 &&
//    he.fDedxSigma > -3 &&
//    v0.fP > 2 && v0.fP < 10 &&
//    v0.fCosPointingAngle > 0.996 &&
//    v0.fCosPointingAngle < 1.0 &&
//    ((he.fP.P() > 1.3 && v0.fCharge > 0) || v0.fCharge < 0 );
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskHypTritEventTree::GeoLength(const AliESDtrack& track) {
  Double_t deadZoneWidth = 3;
  Double_t lengthInActiveZone = track.GetLengthInActiveZone(1, deadZoneWidth, 220, track.GetESDEvent()->GetMagneticField());
  return lengthInActiveZone;
}
//_____________________________________________________________________________
void AliAnalysisTaskHypTritEventTree::SetBetheBlochParams(Int_t runNumber) {
	// set Bethe-Bloch parameter
	if (runNumber >= 252235 && runNumber <= 264347 ) { // 2016 pp
		if(!fMCtrue) { // Data
			// LHC16 + LHC18
			// He3
			fBetheParamsT[0] = 0.427978;
			fBetheParamsT[1] = 105.46;
			fBetheParamsT[2] =-7.08642e-07;
			fBetheParamsT[3] = 2.23332;
			fBetheParamsT[4] = 18.8231;
			fBetheParamsT[5] = 0.06;
			// Triton
			fBetheParamsHe[0] = 1.81085;
			fBetheParamsHe[1] = 29.4656;
			fBetheParamsHe[2] = 0.0458225;
			fBetheParamsHe[3] = 2.08689;
			fBetheParamsHe[4] = 2.28772;
			fBetheParamsHe[5] = 0.06;
		} else { // MC
			if (runNumber >= 262424 || runNumber <= 256418 ) {
				//LHC18a2b (->LHC16)
				// He3
				fBetheParamsHe[0] = 3.05245;
				fBetheParamsHe[1] = 15.7252;
				fBetheParamsHe[2] = -0.00453331;
				fBetheParamsHe[3] = 2.17241;
				fBetheParamsHe[4] = 2.88422;
				fBetheParamsHe[5] = 0.0834274;
				// Triton
				fBetheParamsT[0] = 2.74259;
				fBetheParamsT[1] = 18.3295;
				fBetheParamsT[2] = 5.91594;
				fBetheParamsT[3] = 1.93471;
				fBetheParamsT[4] = 0.292147;
				fBetheParamsT[5] = 0.0728241;
			}
			if (runNumber >= 256941 && runNumber <= 258537 ) { 
				// LHC18a2b2 (LHC16k)
				// He3
				fBetheParamsHe[0] = 2.80527;
				fBetheParamsHe[1] = 14.2379;
				fBetheParamsHe[2] = 0.0232811;
				fBetheParamsHe[3] = 2.11464;
				fBetheParamsHe[4] = 1.615;
				fBetheParamsHe[5] = 0.0815227;
				// Triton
				fBetheParamsT[0] = 1.31603;
				fBetheParamsT[1] = 36.1798;
				fBetheParamsT[2] = 493.036;
				fBetheParamsT[3] = 2.10841;
				fBetheParamsT[4] = 7.43391;
				fBetheParamsT[5] = 0.0769041;
			}
			if (runNumber >= 258962 && runNumber <= 259888 ) {
				//LHC18a2b3 (->LHC16l)
				// He3
				fBetheParamsHe[0] = 2.80121;
				fBetheParamsHe[1] = 14.2397;
				fBetheParamsHe[2] = 0.0100894;
				fBetheParamsHe[3] = 2.10396;
				fBetheParamsHe[4] = 1.41608;
				fBetheParamsHe[5] = 0.0817429;
				// Triton
				fBetheParamsT[0] = 4.80597;
				fBetheParamsT[1] = 13.8813;
				fBetheParamsT[2] = 189.651;
				fBetheParamsT[3] = 2.05969;
				fBetheParamsT[4] = 4.38013;
				fBetheParamsT[5] = 0.077593;
			} 
		}
	}
	if (runNumber >= 270581 && runNumber <= 282704) { // 2017 pp
		if(!fMCtrue) {
			//LHC17
			// He3
			fBetheParamsHe[0] = 3.20025;
			fBetheParamsHe[1] = 16.4971;
			fBetheParamsHe[2] = -0.0116571;
			fBetheParamsHe[3] = 2.3152;
			fBetheParamsHe[4] = 3.11135;
			fBetheParamsHe[5] = 0.06;
			// Triton
			fBetheParamsT[0] = 0.420434;
			fBetheParamsT[1] = 106.102;
			fBetheParamsT[2] = -3.15587e-07;
			fBetheParamsT[3] = 2.32499;
			fBetheParamsT[4] = 21.3439;
			fBetheParamsT[5] = 0.06;	
		} else {
			// LHC18a2a (->LHC17)
			// He3
			fBetheParamsHe[0] = 3.12796;
			fBetheParamsHe[1] = 16.1359;
			fBetheParamsHe[2] = -0.00682978;
			fBetheParamsHe[3] = 2.26624;
			fBetheParamsHe[4] = 2.58652;
			fBetheParamsHe[5] = 0.0847009;
			// Triton
			fBetheParamsT[0] = 2.8303;
			fBetheParamsT[1] = 15.4337;
			fBetheParamsT[2] = 3.18352;
			fBetheParamsT[3] = 2.20975;
			fBetheParamsT[4] = 0.218244;
			fBetheParamsT[5] = 0.0780191;	
		}
	}
	if (runNumber >= 285009 && runNumber <= 294925) { // 2018 pp
		if(!fMCtrue) {
			// LHC16 + LHC18
			// He3
			fBetheParamsT[0] = 0.427978;
			fBetheParamsT[1] = 105.46;
			fBetheParamsT[2] =-7.08642e-07;
			fBetheParamsT[3] = 2.23332;
			fBetheParamsT[4] = 18.8231;
			fBetheParamsT[5] = 0.06;
			// Triton
			fBetheParamsHe[0] = 1.81085;
			fBetheParamsHe[1] = 29.4656;
			fBetheParamsHe[2] = 0.0458225;
			fBetheParamsHe[3] = 2.08689;
			fBetheParamsHe[4] = 2.28772;
			fBetheParamsHe[5] = 0.06;
		} else {
			//LHC18a2d (->LHC18)
			// He3
			fBetheParamsHe[0] = 3.07104;
			fBetheParamsHe[1] = 15.8085;
			fBetheParamsHe[2] = 0.0150992;
			fBetheParamsHe[3] = 2.13909;
			fBetheParamsHe[4] = 2.59495;
			fBetheParamsHe[5] = 0.0865179;
			// Triton
			fBetheParamsT[0] = 2.54486;
			fBetheParamsT[1] = 17.1203;
			fBetheParamsT[2] = -0.0452007;
			fBetheParamsT[3] = 2.00988;
			fBetheParamsT[4] = 0.849292;
			fBetheParamsT[5] = 0.0768715;		
		}
	}
}
