#include <TChain.h>
#include <TParticle.h>

#include "AliAnalysisCODEXtask.h"
#include "AliAnalysisCODEX.h"

#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliMultSelection.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDTOFCluster.h"
#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliTOFPIDResponse.h"
#include "TH2I.h"

#include <climits>

ClassImp(AliAnalysisCODEXtask);

using namespace AliAnalysisCODEX;

AliAnalysisCODEXtask::AliAnalysisCODEXtask(const char* name)
  :AliAnalysisTaskSE(name)
  ,mMCtrue(false)
  ,mCentralityMode(0)
  ,Cuts()
  ,mOutput(0x0)
  ,mTree(0x0)
  ,mPIDresponse(0x0)
  ,mHeader()
  ,mTracks()
  ,mTimeChan(0x0)
  ,mEventCuts(false)
{
  Cuts.SetMinNClustersTPC(60);
  Cuts.SetMaxChi2PerClusterTPC(6);
  Cuts.SetAcceptKinkDaughters(false);
  Cuts.SetRequireTPCRefit(true);
  Cuts.SetRequireITSRefit(false);
  Cuts.SetMaxDCAToVertexZ(3);
  Cuts.SetMaxDCAToVertexXY(3);
  Cuts.SetMaxChi2PerClusterITS(100000000.);
  Cuts.SetMaxChi2TPCConstrainedGlobal(100000000.);
  Cuts.SetEtaRange(-0.8,0.8);
  DefineInput(0, TChain::Class());
  DefineOutput(1, TTree::Class());
  DefineOutput(2, TList::Class());
}


AliAnalysisCODEXtask::~AliAnalysisCODEXtask() {
  if (mPIDresponse) {
    delete mPIDresponse;
    mPIDresponse = 0x0;
  }

  if (mTree) {
    delete mTree;
    mTree = 0x0;
  }

  if (mOutput){
    delete mOutput;
    mOutput = 0x0;
  }

}

void AliAnalysisCODEXtask::UserCreateOutputObjects() {
  // Creation of the histograms, this is called once
  //
  OpenFile(1);
  mTree = new TTree("AliCODEX","Alice COmpressed Dataset for EXotica");
  mTree->Branch("Header",&mHeader);
  mTree->Branch("Tracks",&mTracks);
  //
  mTree->SetAutoSave(100000000);
  PostData(1,mTree);

  OpenFile(2);
  if (!mOutput) mOutput = new TList();
  mOutput->SetOwner();

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

  //
  //Define output histograms
  //
  mTimeChan = new TH2I("mTimeChan", "Histogram with Channel/Time correlation;Channel;TOF (ps)", 3276, 0. -.5, 157248. -.5, 5000, 10000, 80000);
  mOutput->AddLast(mTimeChan);

  mEventCuts.AddQAplotsToList(mOutput);
  PostData(2,mOutput);

}

void AliAnalysisCODEXtask::UserExec(Option_t *){
  mHeader.mEventMask = 0;

  AliESDEvent *event = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!event) return;
  if (!mEventCuts.AcceptEvent(event)) {
    PostData(2, mOutput);
    return;
  }

  AliTOFPIDResponse& tofPID = mPIDresponse->GetTOFResponse();
  for (int iB = 0; iB < 10; ++iB) {
    mHeader.mT0event[iB] = tofPID.GetT0bin(iB);
    mHeader.mT0resolution[iB] = tofPID.GetT0binRes(iB);
    mHeader.mT0mask[iB] = tofPID.GetT0binMask(iB);
  }

  const AliESDVertex *vertex = static_cast<const AliESDVertex*>(mEventCuts.GetPrimaryVertex());
  mHeader.SetVertexZposition(vertex->GetZ());

  const unsigned short standard_mask = (Cuts.GetRequireTPCRefit()) ? kTPCrefit : 0;

  float cent = mEventCuts.GetCentrality(0);
  mHeader.SetCentrality(cent);

  if (event->GetMagneticField() < 0) mHeader.mEventMask |= kNegativeB;

  /// Check Monte Carlo information and other access first:
  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!eventHandler) {
    if (mMCtrue)
      AliFatal("You asked for MC analysis, but I don't find any MCEventHandler... did you forget to add it to your analysis manager?");
  }
  //
  AliMCEvent* mcEvent = nullptr;
  if (eventHandler) mcEvent = eventHandler->MCEvent();
  if (!mcEvent && mMCtrue)
    AliFatal("Missing MC event");
  if (mMCtrue) {
    mHeader.mEventMask |= kMCevent;
  }

  bool first = true;
  mTracks.clear();
  Track t;
  for (int iEv = 0;iEv < event->GetNumberOfTracks(); ++iEv) {
    AliESDtrack *track = event->GetTrack(iEv);

    /// Tracking quality cuts
    if (!track->GetInnerParam()) continue;
    if (!Cuts.AcceptTrack(track)) continue;

    /// PID cuts
    float sigEl = mPIDresponse->NumberOfSigmas(AliPIDResponse::kTPC,track,AliPID::kElectron);
    float sigPi = mPIDresponse->NumberOfSigmas(AliPIDResponse::kTPC,track,AliPID::kPion);
    float sigKa = mPIDresponse->NumberOfSigmas(AliPIDResponse::kTPC,track,AliPID::kKaon);
    float sigPr = mPIDresponse->NumberOfSigmas(AliPIDResponse::kTPC,track,AliPID::kProton);
    float sigDe = mPIDresponse->NumberOfSigmas(AliPIDResponse::kTPC,track,AliPID::kDeuteron);
    float sigH3  = mPIDresponse->NumberOfSigmas(AliPIDResponse::kTPC,track,AliPID::kTriton);
    float sigHe3 = mPIDresponse->NumberOfSigmas(AliPIDResponse::kTPC,track,AliPID::kHe3);
    float sigHe4 = mPIDresponse->NumberOfSigmas(AliPIDResponse::kTPC,track,AliPID::kAlpha);

    /// Cut tracks without TOF matching at moderate high pT
    float time = 0.;
    int   channel = -1;
    if (!(track->GetStatus() & AliESDtrack::kTOFout) || !(track->GetStatus() & AliESDtrack::kTIME)) {
      time = -1.f;
    } else {
      float len = track->GetIntegratedLength();
      if (len < 350.f) {
        time = -1.f;
      } else {
        float time0 = mPIDresponse->GetTOFResponse().GetStartTime(track->P());
        time = track->GetTOFsignal() - time0;

        // TOF number of clusters
        t.SetTOFNcls(track->GetNTOFclusters());

        // Momentum at the TOF
        const double t_d = tofPID.GetExpectedSignal(track, AliPID::kDeuteron);
        const double beta_d = len / (t_d * kCtof);
        t.SetTOFmomentum((fabs(beta_d - 1.f) < 1.e-24) ? track->GetTPCmomentum() : kDeuteronMass * beta_d / sqrt(1. - (beta_d * beta_d)));

        channel = track->GetTOFCalChannel();
        // Get mismatch signal
        if (t.GetTOFNcls() > 1) {
          int idx = track->GetTOFclusterArray()[t.TOFnClusters - 1];
          AliESDTOFCluster* clust = (AliESDTOFCluster *)event->GetESDTOFClusters()->At(idx);
          int clsindex = -1;
          for(Int_t iMTr = 0; iMTr < clust->GetNMatchableTracks(); iMTr++){
            if(clust->GetTrackIndex(iMTr) == track->GetID()){
              clsindex = iMTr;
              break;
            }
          }
          if (clsindex > -1) {
            float TOFmismatch = clust->GetTime(clsindex) * clust->GetLength(clsindex) / len;
            mTimeChan->Fill(channel, TOFmismatch);
          }
        }
      }
    }
    t.SetTOFchannel(channel);
    t.SetTOFsignal(time);
    t.SetIntegratedLength(track->GetIntegratedLength());

    if (track->GetNumberOfTRDClusters() > 0)
      t.mask |= AliAnalysisCODEX::kTRDout;

    /// Put the right information in the right place
    t.eta = track->Eta();
    t.phi = track->Phi();
    t.pTPC = track->GetInnerParam()->GetP();
    t.pT = track->Pt();
    if (track->Charge() < 0) t.pT = -t.pT;

    /// PID information
    //for (int i = 0; i < 8; ++i) t.TPCsigmas[i] = SCHAR_MAX;
    t.SetTPCsigma(kEl,sigEl);
    t.SetTPCsigma(kPi,sigPi);
    t.SetTPCsigma(kKa,sigKa);
    t.SetTPCsigma(kPr,sigPr);
    t.SetTPCsigma(kDe,sigDe);
    t.SetTPCsigma(kH3,sigH3);
    t.SetTPCsigma(kHe3,sigHe3);
    t.SetTPCsigma(kHe4,sigHe4);
    t.TPCsignal = track->GetTPCsignal();

    /// Number of clusters
    t.TPCnClusters = track->GetTPCNcls();
    /// Number of findable clusters
    t.TPCnClustersF = track->GetTPCNclsF();
    /// Number of crossed rows
    t.TPCnXedRows = track->GetTPCCrossedRows();
    t.TPCnSignal = track->GetTPCsignalN();
    t.ITSmap = track->GetITSClusterMap();
    if ((track->GetStatus() & AliESDtrack::kITSrefit) == AliESDtrack::kITSrefit) t.ITSmap |= AliAnalysisCODEX::kITSrefit;
    /// Track length in active region
    t.ActiveLength = track->GetLengthInActiveZone(0, 3, 220., event->GetMagneticField());

    /// Mask
    t.mask = standard_mask | kIsReconstructed;

    /// Binned information
    float cov[3],dca[2];
    track->GetImpactParameters(dca,cov);
    t.SetDCAxy(dca[0]);
    t.SetDCAz(dca[1]);
    t.SetTPCChi2NDF(track->GetTPCchi2() / track->GetTPCNcls());
    t.SetITSChi2NDF((track->GetITSNcls()) ? track->GetITSchi2() / track->GetITSNcls() : 1.e9);
    t.SetGoldenChi2NDF(track->GetChi2TPCConstrainedVsGlobal(vertex));

    if (mMCtrue) {
      int label = track->GetLabel();
      TParticle* part = mcEvent->Particle(abs(label));
      if (!part) continue;
      int particle_mask = GetParticleMask(part);
      if (!particle_mask) continue;
      t.mask |= particle_mask;

      if (label < 0) t.mask |= kIsFake;

      if (time > 0) {
        int TOFlabel[3];
        track->GetTOFLabel(TOFlabel);
        bool TOFmismatch = true;
        for (int i = 0; i < 3; ++i)
          if (TOFlabel[i] > 0 && TOFlabel[i] == abs(label))
            TOFmismatch = false;
        if (TOFmismatch) t.mask |= kTOFmismatch;
      }

      if (mcEvent->IsPhysicalPrimary(abs(label))) t.mask |= kIsPrimary;
      else if (mcEvent->IsSecondaryFromMaterial(abs(label))) t.mask |= kIsSecondaryFromMaterial;
    }

    /// If everything went OK pushing the track
    mTracks.push_back(t);
  }

  if (mMCtrue) {
    t.TPCnClusters = 0;
    t.TPCnSignal = 0;
    for (int i = 0; i < 8; ++i) t.TPCsigmas[i] = SCHAR_MAX;
    t.TOFsignal = -1;
    t.ITSmap = 0;
    t.DCAxy = 0;
    t.DCAz = 0;
    t.TPCchi2NDF = 0;
    t.ITSchi2NDF = 0;
    for (int iP = 0; iP < mcEvent->GetNumberOfTracks(); ++iP) {
      TParticle* particle = mcEvent->Particle(iP);
      if (!particle)
        continue;

      if (particle->Energy() <= fabs(particle->Pz())) continue; // fix improper TParticle->Y() behaviour
      if (fabs(particle->Y()) > 1.)
        continue;

      int particle_mask = GetParticleMask(particle);
      if (!particle_mask)
        continue;
      t.mask = particle_mask;
      if (mcEvent->IsPhysicalPrimary(iP)) t.mask |= kIsPrimary;
      else if (mcEvent->IsSecondaryFromMaterial(iP)) t.mask |= kIsSecondaryFromMaterial;

      t.eta = particle->Eta();
      t.phi = particle->Phi();
      t.pT = (particle->GetPdgCode() > 0) ? particle->Pt() : -particle->Pt();
      t.pTPC = particle->P();

      /// Pushing the MC particle
      mTracks.push_back(t);
    }
  }
  mTree->Fill();
  PostData(1,mTree);

  PostData(2,mOutput);
}

long AliAnalysisCODEXtask::GetParticleMask(TParticle *part) {
  int pdg = part->GetPdgCode();
  switch (abs(pdg)) {
    case 11: return kEl;
    case 211: return kPi;
    case 321: return kKa;
    case 2212: return kPr;
    case 1000010020: return kDe;
    case 1000010030: return kH3;
    case 1000020030: return kHe3;
    case 1000020040: return kHe4;
  }
  return 0;
}
