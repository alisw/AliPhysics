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
  fMCtrue(0),
  fEventCuts(),
  fPeriod(),
  fBetheSplines(),
  fBetheParamsHe(),
  fBetheParamsT(),
  fPidQa(0),
  fUseAnalysisTrkSel(kTRUE),
  fPIDCheckOnly(kFALSE) {

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
  fMCtrue(0),
  fEventCuts(),
  fPeriod(),
  fBetheSplines(),
  fBetheParamsHe(),
  fBetheParamsT(),
  fPidQa(0),
  fUseAnalysisTrkSel(kTRUE),
  fPIDCheckOnly(kFALSE) {
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
  fHistdEdx = new TH2F("fHistdEdX","dE/dx;#frac{#it{p}}{z} (GeV/#it{c});TPC Signal (a.u.)",1000,-8.0,8.0,1000,0.0,1200);
  fHistdEdxV0 = new TH2F("fHistdEdXV0","dE/dx;#frac{#it{p}}{z} (GeV/#it{c});TPC Signal (a.u.)",1000,-8.0,8.0,1000,0.0,1200);
  fHistNumEvents = new TH1F("fHistNumEvents","Number of Events",2,0,2);
  fHistNumEvents->GetXaxis()->SetBinLabel(1,"before PhysSel");
  fHistNumEvents->GetXaxis()->SetBinLabel(2,"after PhysSel");

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
  if (fPeriod == 2010 || fPeriod == 2011) {
      if (vertex->GetNContributors() < 1) {
        vertex = fESDevent->GetPrimaryVertexSPD();
        if (vertex->GetNContributors() < 1) {
          PostData(1,fHistogramList);
          return;
        }
      }
    if (TMath::Abs(vertex->GetZ()) > 10) {
      PostData(1, fHistogramList);
      return;
    }
    centrality = fESDevent->GetCentrality()->GetCentralityPercentile("V0M");
    if(!fMCtrue){
      if (centrality < 0.0 || centrality > 100.0 ) {
        return;
      }
    }
  }
  if (fPeriod == 2015) {
    if(!fEventCuts.AcceptEvent(fESDevent)) {
      PostData(1,fHistogramList);
      return;
    }
    // 0 = V0M
    centrality = fEventCuts.GetCentrality(0);
    if(!fMCtrue){
      if (centrality < 0.0 || centrality > 100.0 ) {
        return;
      }
    }
  }
  if (fPeriod == 2016) {
    fEventCuts.SetManualMode();
    fEventCuts.SetupRun2pp();

    if(!fEventCuts.AcceptEvent(fESDevent)) {
      PostData(1,fHistogramList);
      return;
    }
    // 0 = V0M
    centrality = fEventCuts.GetCentrality(0);
  }

  Int_t runNumber = fESDevent->GetRunNumber();
  TriggerSelection();
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
    trackCutsPid = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kTRUE,0);
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
    if (trackN->GetInnerParam()->GetP() > 10 || trackP->GetInnerParam()->GetP() > 10) continue;
    Bool_t pionPositive     = kFALSE;
    Bool_t pionNegative     = kFALSE;
    Bool_t protonPositive   = kFALSE;
    Bool_t protonNegative   = kFALSE;
    Bool_t deuteronPositive = kFALSE;
    Bool_t deuteronNegative = kFALSE;
    Bool_t tritonPositive   = kFALSE;
    Bool_t tritonNegative   = kFALSE;
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
  fTree->Fill();
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

void AliAnalysisTaskHypTritEventTree::Terminate(const Option_t*) {
  if (!GetOutputData(0)) return;
}

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
      if (McCuts(*reducedV0, *reducedHe, *reducedPi)) fHistMcRec->Fill(reducedV0->fDecayLength);
    }
  }
}

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
        reducedV0->fDecayLength = distance * tparticleMother->GetCalcMass() / tparticleMother->P();
        reducedV0->fMcTruth = 0;
        fHistMcGen->Fill(reducedV0->fDecayLength);
        for(Int_t ii = 0; ii < 40; ii++) {
          if (fMCGenRec[ii] == istack) {
            AliReducedHypTritV0 *mcReducedV0 = (AliReducedHypTritV0*) fMCGenRecArray->At(ii);
            AliReducedHypTritTrack *mcReducedPi = mcReducedV0->Pi();
            AliReducedHypTritTrack *mcReducedHe = mcReducedV0->He();
            reducedV0->fDecayLength = mcReducedV0->fDecayLength;
            reducedV0->fDcaV0 = mcReducedV0->fDcaV0;
            reducedV0->fCosPointingAngle = mcReducedV0->fCosPointingAngle;
            reducedPi->fP = mcReducedPi->fP;
            reducedHe->fP = mcReducedHe->fP;
            reducedPi->fEta = mcReducedPi->fEta;
            reducedHe->fEta = mcReducedHe->fEta;
            reducedPi->fTpcNClusters = mcReducedPi->fTpcNClusters;
            reducedHe->fTpcNClusters = mcReducedHe->fTpcNClusters;
            reducedV0->fRapidity = mcReducedV0->fRapidity;
            reducedV0->fM = mcReducedV0->fM;
            reducedV0->fMcTruth = mcReducedV0->fMcTruth;
            reducedV0->fParticleSpecies = mcReducedV0->fParticleSpecies;
            reducedV0->fCharge = mcReducedV0->fCharge;
            reducedHe->fDedxSigmaTriton = mcReducedHe->fDedxSigmaTriton;
            reducedHe->fDedxSigma = mcReducedHe->fDedxSigma;
            reducedHe->fPtrack = mcReducedHe->fPtrack;
            reducedPi->fPtrack = mcReducedHe->fPtrack;
            reducedV0->fOnFlyStatus = mcReducedV0->fOnFlyStatus;
          }
        }
        nV0Gen = nV0Gen +1;
      }
		}
	}
	fReducedEventMCGen->fNumberV0s = nV0Gen;
	fTreeMCGen->Fill();
	fReducedEventMCGen->ClearEvent();
	fMCGenRecArray->Clear();
}

/// Set trigger information in reduced event
/// \return returns kTRUE is successful.
Bool_t AliAnalysisTaskHypTritEventTree::TriggerSelection() {
  fReducedEvent->fTrigger = 0;
  TString classes = fESDevent->GetFiredTriggerClasses();

  if ((fInputHandler->IsEventSelected() & AliVEvent::kINT7)) fReducedEvent->fTrigger = 1;
  if ((fInputHandler->IsEventSelected() & AliVEvent::kCentral)) fReducedEvent->fTrigger = 2;
  if ((fInputHandler->IsEventSelected() & AliVEvent::kSemiCentral)) fReducedEvent->fTrigger =  3;
  if(classes.Contains("HNU") || classes.Contains("HQU") || classes.Contains("HSE")) fReducedEvent->fTrigger = 4;
  Bool_t isTriggered = kTRUE;
  if (fReducedEvent->fTrigger == 0) isTriggered = kFALSE;
  return isTriggered;
}
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

Double_t AliAnalysisTaskHypTritEventTree::GeoLength(const AliESDtrack& track) {
  Double_t deadZoneWidth = 3;
  Double_t lengthInActiveZone = track.GetLengthInActiveZone(1, deadZoneWidth, 220, track.GetESDEvent()->GetMagneticField());
  return lengthInActiveZone;
}
