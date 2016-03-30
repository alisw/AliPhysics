#include "AliAnalysisTaskAbsorptionStudies.h"
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TParticle.h>
#include <TArrayD.h>

#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"

#include <climits>

enum {
  kPi = 211,
  kKa = 321,
  kPr = 2112,
  kDe = 1000010020
};

ClassImp(AliAnalysisTaskAbsorptionStudies)

AliAnalysisTaskAbsorptionStudies::AliAnalysisTaskAbsorptionStudies(const char* name)
:AliAnalysisTaskSE(name),
  mTOFtimeAnalysis(false),
  mMCtrue(false),
  mMagField(0),
  mCuts(),
  mCentrality(5,0,50),
  mMassAxis(800,0,4),
  mPtAxis(20,1,3),
  mPhiAxis(64,0,TMath::TwoPi()),
  mTimeAxis(3000,12500,42500),
  mDCAxyAxis(200,-1,1),
  mOutput(),
  mPIDresponse(),
  mEventCounter(0x0),
  mPhiPtMass(0x0),
  mPhiPtTime(0x0),
  mPhiPtTimeParticles(0x0),
  mPhiPtDCAxy(0x0),
  mPhiPtDCAxyWeak(0x0),
  mPhiPtDCAxyMaterial(0x0),
  mPhiPtGen(0x0),
  mPhiPtRec(0x0)
{
  mCuts.SetMinNClustersTPC(60);
  mCuts.SetMaxChi2PerClusterTPC(6);
  mCuts.SetAcceptKinkDaughters(false);
  mCuts.SetRequireTPCRefit(true);
  mCuts.SetRequireITSRefit(true);
  mCuts.SetMaxDCAToVertexZ(1);
  mCuts.SetMaxDCAToVertexXY(0.5);
  mCuts.SetMaxChi2PerClusterITS(36);
  mCuts.SetEtaRange(-0.8,0.8);
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());

}


AliAnalysisTaskAbsorptionStudies::~AliAnalysisTaskAbsorptionStudies() {
  if (mPIDresponse) {
    delete mPIDresponse;
    mPIDresponse = 0x0;
  }

  if (mOutput) {
    delete mOutput;
    mOutput = 0x0;
  }
}

void AliAnalysisTaskAbsorptionStudies::UserCreateOutputObjects() {
  // Creation of the histograms, this is called once
  //
  mOutput = new TList();
  mOutput->SetOwner();

  /// Positive/Negative, species, centrality, phi, pt, dcaxy
  int nBinsDCAxy[6]      = {  2  ,  4  , mCentrality.GetNbins(), mPhiAxis.GetNbins(), mPtAxis.GetNbins(), mDCAxyAxis.GetNbins()};
  double minBinsDCAxy[6] = { -0.5, -0.5, mCentrality.GetXmin(), mPhiAxis.GetXmin(), mPtAxis.GetXmin(), mDCAxyAxis.GetXmin()};
  double maxBinsDCAxy[6] = {  1.5,  3.5, mCentrality.GetXmax(), mPhiAxis.GetXmax(), mPtAxis.GetXmax(), mDCAxyAxis.GetXmax()};
  mPhiPtDCAxy = new THnSparseF("mPhiPtDCAxy","",6,nBinsDCAxy,minBinsDCAxy,maxBinsDCAxy);
  mOutput->Add(mPhiPtDCAxy);

  mEventCounter = new TH1F("mEventCounter",";Centrality;Events",mCentrality.GetNbins(),mCentrality.GetXmin(),mCentrality.GetXmax());
  mOutput->Add(mEventCounter);
  if (mMCtrue) {
    mPhiPtDCAxyWeak = new THnSparseF("mPhiPtDCAxyWeak","",6,nBinsDCAxy,minBinsDCAxy,maxBinsDCAxy);
    mPhiPtDCAxyMaterial = new THnSparseF("mPhiPtDCAxyMaterial","",6,nBinsDCAxy,minBinsDCAxy,maxBinsDCAxy);
    mOutput->Add(mPhiPtDCAxy);
    mOutput->Add(mPhiPtDCAxyWeak);
    mOutput->Add(mPhiPtDCAxyMaterial);

    /// Positive/Negative, species, centrality, phi, pt
    mPhiPtGen = new THnSparseF("mPhiPtGen","",5,nBinsDCAxy,minBinsDCAxy,maxBinsDCAxy);
    mPhiPtRec = new THnSparseF("mPhiPtRec","",5,nBinsDCAxy,minBinsDCAxy,maxBinsDCAxy);
    mOutput->Add(mPhiPtGen);
    mOutput->Add(mPhiPtRec);
  } else {
    /// Positive/Negative, centrality, phi, pt, mass
    int nBinsMass[6]      = {  2  ,  2  , mCentrality.GetNbins(), mPhiAxis.GetNbins(), mPtAxis.GetNbins(), mMassAxis.GetNbins()};
    double minBinsMass[6] = { -0.5, -0.5, mCentrality.GetXmin(), mPhiAxis.GetXmin(), mPtAxis.GetXmin(), mMassAxis.GetXmin()};
    double maxBinsMass[6] = {  1.5,  1.5, mCentrality.GetXmax(), mPhiAxis.GetXmax(), mPtAxis.GetXmax(), mMassAxis.GetXmax()};
    mPhiPtMass = new THnSparseF("mPhiPtMass","",6,nBinsMass,minBinsMass,maxBinsMass);
    mOutput->Add(mPhiPtMass);

    if (mTOFtimeAnalysis) {
      /// Positive/Negative, proper time / mismatch ,centrality, phi, pt, time
      int nBinsTime[6]      = {  2  ,  2  , mCentrality.GetNbins(), mPhiAxis.GetNbins(), mPtAxis.GetNbins(), mTimeAxis.GetNbins()};
      double minBinsTime[6] = { -0.5, -0.5, mCentrality.GetXmin(), mPhiAxis.GetXmin(), mPtAxis.GetXmin(), mTimeAxis.GetXmin()};
      double maxBinsTime[6] = {  1.5,  1.5, mCentrality.GetXmax(), mPhiAxis.GetXmax(), mPtAxis.GetXmax(), mTimeAxis.GetXmax()};
      mPhiPtTime = new THnSparseF("mPhiPtTime","",6,nBinsTime,minBinsTime,maxBinsTime);
      mOutput->Add(mPhiPtTime);

      /// Positive/Negative, species, centrality, phi, pt, time
      int nBinsTimeSpecies[6]      = {  2  ,  4  , mCentrality.GetNbins(), mPhiAxis.GetNbins(), mPtAxis.GetNbins(), mTimeAxis.GetNbins()};
      double minBinsTimeSpecies[6] = { -0.5, -0.5, mCentrality.GetXmin(), mPhiAxis.GetXmin(), mPtAxis.GetXmin(), mTimeAxis.GetXmin()};
      double maxBinsTimeSpecies[6] = {  1.5,  3.5, mCentrality.GetXmax(), mPhiAxis.GetXmax(), mPtAxis.GetXmax(), mTimeAxis.GetXmax()};
      mPhiPtTimeParticles = new THnSparseF("mPhiPtTimeParticles","",6,nBinsTimeSpecies,minBinsTimeSpecies,maxBinsTimeSpecies);
      mOutput->Add(mPhiPtTimeParticles);
    }
  }
  PostData(1,mOutput);
  //
  // Get PID response object
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  if(!man)
    AliFatal("Could not find manager");
  AliInputEventHandler* handler = dynamic_cast<AliInputEventHandler*> (man->GetInputEventHandler());
  if(!handler)
    AliFatal("No input event handler");
  mPIDresponse = dynamic_cast<AliPIDResponse*>(handler->GetPIDResponse());
  if (!mPIDresponse)
    AliFatal("PIDResponse object was not created"); // Escalated to fatal. Task is unusable without PID response.
}

void AliAnalysisTaskAbsorptionStudies::UserExec(Option_t *){
  AliESDEvent *event = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!event) return;

  const AliESDVertex *vertex = event->GetPrimaryVertex();
  if (!vertex) return;
  if (fabs(vertex->GetZ()) > 10.f) return;
  if (vertex->GetNContributors() < 1) return;

  AliCentrality *f_centrality = event->GetCentrality();
  if (!f_centrality)  AliFatal("Missing centrality object.");
  float centrality = char(f_centrality->GetCentralityPercentile("V0M"));
  if (centrality < mCentrality.GetXmin() || centrality > mCentrality.GetXmax()) return;
  mEventCounter->Fill(centrality);
  if (event->GetMagneticField() * mMagField < 0) return;

  //
  // Check Monte Carlo information and other access first:
  //
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!eventHandler && mMCtrue) {
    AliFatal("Missing MC event handler");
  }
  //
  AliMCEvent* mcEvent = 0x0;
  AliStack* stack = 0x0;
  if (eventHandler) mcEvent = eventHandler->MCEvent();
  if (!mcEvent && mMCtrue)
    AliFatal("Missing MC event");
  if (mMCtrue) {
    stack = mcEvent->Stack();
    if (!stack)
      AliFatal("Missing MC stack");
  }

  for (int iEv = 0;iEv < event->GetNumberOfTracks(); ++iEv) {
    AliESDtrack *track = event->GetTrack(iEv);

    /// Tracking quality cuts
    if (!track->GetInnerParam()) continue;
    if (!mCuts.AcceptTrack(track)) continue;

    /// Charge
    int iCharge = 0;
    if (track->Charge() < 0) iCharge = 1;

    /// Cut tracks without TOF matching at moderate high pT
    float len = track->GetIntegratedLength();
    if (!(track->GetStatus() & AliESDtrack::kTOFout) ||
        !(track->GetStatus() & AliESDtrack::kTIME) ||
        len < 350.f)
      continue;

    const float time0 = mPIDresponse->GetTOFResponse().GetStartTime(track->P());
    const float time = track->GetTOFsignal() - time0;
    if (time < 1.e-3)
      continue;
    float beta = len / (time * 2.99792457999999984e-02);
    if (beta < 1.e-10)
      continue;

    float mass = track->GetInnerParam()->GetP() * track->GetInnerParam()->GetP() * (1.f / (beta * beta) - 1.f);

    float dcaxy = 10.f,dcaz = 10.f;
    track->GetImpactParameters(dcaxy,dcaz);
    double fill[6] = {double(iCharge), 0., centrality, track->Phi(), track->Pt(), 0.};
    if (mMCtrue) {
      int label = track->GetLabel();
      TParticle* part = stack->Particle(abs(label));
      if (!part) continue;
      fill[1] = GetParticleId(part->GetPdgCode());
      mPhiPtRec->Fill(fill);
      fill[5] = dcaxy;
      if (stack->IsPhysicalPrimary(abs(label))) mPhiPtDCAxy->Fill(fill);
      if (stack->IsSecondaryFromMaterial(abs(label))) mPhiPtDCAxyMaterial->Fill(fill);
      if (stack->IsSecondaryFromWeakDecay(abs(label))) mPhiPtDCAxyWeak->Fill(fill);
      /*int TOFlabel[3];
        track->GetTOFLabel(TOFlabel);
        bool TOFmismatch = true;
        for (int i = 0; i < 3; ++i)
        if (TOFlabel[i] > 0 && TOFlabel[i] == abs(label))
        TOFmismatch = false;
        */
    } else {
      float nSigmaDe = mPIDresponse->NumberOfSigmasTPC(track,AliPID::kDeuteron);
      if (fabs(mPIDresponse->NumberOfSigmasTPC(track,AliPID::kHe3)) < 3) continue;

      fill[1] = 0.;
      fill[5] = mass;
      if (fabs(nSigmaDe) < 2.5) {
        fill[1] = 1.;
      }
      mPhiPtMass->Fill(fill);

      fill[5] = dcaxy;
      mPhiPtDCAxy->Fill(fill);

      if (mTOFtimeAnalysis) {
        AliTOFPIDResponse* tofPID = &(mPIDresponse->GetTOFResponse());
        fill[1] = 0.;
        fill[5] = time;
        mPhiPtTime->Fill(fill);
        fill[1] = 1.;
        fill[5] = tofPID->GetMismatchRandomValue(track->Eta());
        mPhiPtTime->Fill(fill);

        fill[1] = 0.;
        fill[5] = tofPID->GetExpectedSignal(track,AliPID::kPion);
        mPhiPtTimeParticles->Fill(fill);
        fill[1] = 1.;
        fill[5] = tofPID->GetExpectedSignal(track,AliPID::kKaon);
        mPhiPtTimeParticles->Fill(fill);
        fill[1] = 2.;
        fill[5] = tofPID->GetExpectedSignal(track,AliPID::kProton);
        mPhiPtTimeParticles->Fill(fill);
        fill[1] = 3.;
        fill[5] = tofPID->GetExpectedSignal(track,AliPID::kDeuteron);
        mPhiPtTimeParticles->Fill(fill);
      }
    }

  }

  if (mMCtrue) {
    for (int iP = 0; iP < stack->GetNtrack(); ++iP) {
      TParticle* track = stack->Particle(iP);
      if (!track)
        continue;

      if (fabs(track->Y()) > 0.5)
        continue;

      if (!stack->IsPhysicalPrimary(iP))
        continue;

      int iCharge = track->GetPdgCode() > 0 ? 0 : 1;
      double fill[5] = {double(iCharge), (double)GetParticleId(track->GetPdgCode()), centrality, track->Phi(), track->Pt()};
      mPhiPtGen->Fill(fill);
    }
  }
  PostData(1,mOutput);
}

int AliAnalysisTaskAbsorptionStudies::GetParticleId(int pdg) {
  switch (abs(pdg)) {
    case kPi:
      return 0;
    case kKa:
      return 1;
    case kPr:
      return 2;
    case kDe:
      return 3;
    default:
      return -1;
  }
}

