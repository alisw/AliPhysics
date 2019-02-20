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
#include "AliMultSelectionTask.h"
#include "TH2I.h"

#include <climits>

ClassImp(AliAnalysisCODEXtask);

using namespace AliAnalysisCODEX;

AliAnalysisCODEXtask::AliAnalysisCODEXtask(const char* name) :
AliAnalysisTaskSE(name),
mMCtrue{false},
mCentralityMode{0},
Cuts{},
mEventCuts{},
mPtCut{0.1},
mITSsaPtCut{0.4},
mPOI{255},
mEventPOI{},
mNsigmaTPCselectionPOI{5.},
mNsigmaTOFselectionPOI{10.},
mStartingPtTOFselection{1.e4},
mSkipEmptyEvents{false},
mITSstandalone{true},
mOutput{nullptr},
mTree{nullptr},
mPIDresponse{nullptr},
mHeader{},
mTracks{},
mTimeChan{nullptr},
mToDiscard{}
{
  // Cuts.SetMinNClustersTPC(60);
  // Cuts.SetMaxChi2PerClusterTPC(6);
  Cuts.SetAcceptKinkDaughters(false);
  Cuts.SetRequireITSRefit(true);
  if (!mITSstandalone) Cuts.SetRequireTPCRefit(true);
  Cuts.SetMaxDCAToVertexZ(2);
  Cuts.SetMaxDCAToVertexXY(1.5);
  Cuts.SetMaxChi2PerClusterITS(36);
  Cuts.SetMaxChi2TPCConstrainedGlobal(100000000.);
  Cuts.SetEtaRange(-0.8,0.8);
  DefineInput(0, TChain::Class());
  DefineOutput(1, TTree::Class());
  DefineOutput(2, TList::Class());
}


AliAnalysisCODEXtask::~AliAnalysisCODEXtask() {
  delete mPIDresponse;
  delete mTree;
  delete mOutput;
}

void AliAnalysisCODEXtask::UserCreateOutputObjects() {
  // Creation of the histograms, this is called once
  //
  OpenFile(1);
  mTree = new TTree("AliCODEX","Alice COmpressed Dataset for EXotica");
  mTree->Branch("Header",&mHeader);
  mTree->Branch("Tracks",&mTracks);
  Discard();
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
  const AliPID::EParticleType particle_species[8] = {AliPID::kElectron,AliPID::kPion,AliPID::kKaon,AliPID::kProton,AliPID::kDeuteron,AliPID::kTriton,AliPID::kHe3,AliPID::kAlpha};

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

  if (AliMultSelectionTask::IsINELgtZERO(event)) {
    mHeader.mEventMask |= kInelGt0;
  }

  if (mEventCuts.PassedCut(AliEventCuts::kTriggerClasses))
    mHeader.mEventMask |= kTriggerClasses;

  bool EventWithPOI = !bool(mEventPOI);

  mTracks.clear();
  Track t;
  for (int iEv = 0;iEv < event->GetNumberOfTracks(); ++iEv) {
    AliESDtrack *track = event->GetTrack(iEv);

    t.ITSmap = track->GetITSClusterMap();
    if ((track->GetStatus() & AliESDtrack::kITSrefit) == AliESDtrack::kITSrefit) t.ITSmap |= AliAnalysisCODEX::kITSrefit;
    const unsigned short standard_mask = ((track->GetStatus() & AliESDtrack::kTPCrefit) == AliESDtrack::kTPCrefit) ? kTPCrefit : 0;

    /// Tracking quality cuts
    if (!Cuts.AcceptTrack(track)) continue;
    if (track->Pt() < mPtCut) continue;

    /// PID cuts
    float sig[8] = {999.f,999.f,999.f,999.f,999.f,999.f,999.f,999.f};
    bool reject = !mMCtrue; /// In the MC the cut on the TPC pid is replaced by a cut on the true MC particle
    bool tpcRefit = (track->GetStatus() & AliESDtrack::kTPCrefit) == AliESDtrack::kTPCrefit;
    for (int iS = 0; iS < 8; ++iS) {
      if (tpcRefit && track->GetInnerParam() && track->GetTPCNcls() > 60) {
        sig[iS] = mPIDresponse->NumberOfSigmas(AliPIDResponse::kTPC,track,particle_species[iS]);
      } else if (mITSstandalone && !tpcRefit && track->Pt() < mITSsaPtCut) {
        sig[iS] = mPIDresponse->NumberOfSigmas(AliPIDResponse::kITS,track,particle_species[iS]);
      }
      if (std::abs(sig[iS]) < mNsigmaTPCselectionPOI) {
        reject = !bool(mPOI & BIT(iS));
        EventWithPOI = bool(mEventPOI & BIT(iS)) || EventWithPOI;
      }
    }
    if (reject) continue;

    reject = track->Pt() >= mStartingPtTOFselection;
    for (int iS = 0; iS < 8; ++iS) {
      double sigTOF = mPIDresponse->NumberOfSigmas(AliPIDResponse::kTOF,track,particle_species[iS]);
      if (std::abs(sigTOF) < mNsigmaTOFselectionPOI) {
        reject = !bool(mPOI & BIT(iS));
        EventWithPOI = bool(mEventPOI & BIT(iS)) || EventWithPOI;
      }
    }
    if (reject) continue;

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
        t.SetTOFmomentum((fabs(beta_d - 1.f) < 1.e-12) ? track->GetTPCmomentum() : kDeuteronMass * beta_d / sqrt(1. - (beta_d * beta_d)));

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

    /// Put the right information in the right place
    t.eta = track->Eta();
    t.phi = track->Phi();
    t.pTPC = track->GetInnerParam()->GetP();
    t.pT = track->Pt();
    if (track->Charge() < 0) t.pT = -t.pT;

    /// PID information
    for (int iS = 0; iS < 8; ++iS) {
      t.SetTPCsigma(static_cast<AliAnalysisCODEX::BitMask>(BIT(iS)),sig[iS]);
    }
    t.TPCsignal = track->GetTPCsignal();

    if ((track->GetStatus() & AliESDtrack::kTPCrefit) == AliESDtrack::kTPCrefit) {
      /// Number of clusters
      t.TPCnClusters = track->GetTPCNcls();
      /// Number of findable clusters
      t.TPCnClustersF = track->GetTPCNclsF();
      /// Number of crossed rows
      t.TPCnXedRows = track->GetTPCCrossedRows();
      t.TPCnSignal = track->GetTPCsignalN();
    }

    /// Track length in active region
    t.ActiveLength = track->GetLengthInActiveZone(0, 3, 220., event->GetMagneticField());

    /// Mask
    t.mask = standard_mask | kIsReconstructed;

    if (track->GetStatus() & AliVTrack::kTRDrefit) t.mask |= AliAnalysisCODEX::kTRDrefit;
    
    /// Binned information
    float cov[3],dca[2];
    double ITSsamp[4];
    track->GetImpactParameters(dca,cov);
    track->GetITSdEdxSamples(ITSsamp);
    for(int iL = 0; iL < 4; iL++) t.ITSSignal[iL] = ITSsamp[iL];
    t.SetDCAxy(dca[0]);
    t.SetDCAz(dca[1]);
    t.SetTPCChi2NDF(track->GetTPCchi2() / track->GetTPCNcls());
    t.SetITSChi2NDF((track->GetITSNcls()) ? track->GetITSchi2() / track->GetITSNcls() : 1.e9);
    t.SetGoldenChi2NDF(track->GetChi2TPCConstrainedVsGlobal(vertex));
    t.SetTRDnTracklets(track->GetTRDntracklets());

    if (mMCtrue) {
      int label = track->GetLabel();
      TParticle* part = mcEvent->Particle(abs(label));
      if (!part) continue;
      if (part->Pt() < mPtCut) continue;
      int particle_mask = GetParticleMask(part);
      if (!particle_mask || !(particle_mask & mPOI)) continue; /// Reject all the particles absent in the mPOI mask
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
      if (particle->Pt() < mPtCut) continue;

      if (particle->Energy() <= fabs(particle->Pz())) continue; // fix improper TParticle->Y() behaviour
      if (fabs(particle->Y()) > 1.)
      continue;

      int particle_mask = GetParticleMask(particle);
      if (!particle_mask || !(particle_mask & mPOI))
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
  if (!(mSkipEmptyEvents && mTracks.empty()) && EventWithPOI) {
    mTree->Fill();
  }
  PostData(1,mTree);

  PostData(2,mOutput);
}

long AliAnalysisCODEXtask::GetParticleMask(TParticle *part) {
  int pdg = part->GetPdgCode();
  switch (std::abs(pdg)) {
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

void AliAnalysisCODEXtask::Discard(){
  if(mToDiscard.IsNull() || mToDiscard.IsWhitespace()) return;
  TObjArray *arr = mToDiscard.Tokenize(" ");
  for(Int_t i = 0; i < arr->GetEntries(); i++) mTree->SetBranchStatus(static_cast<TObjString*>(arr->At(i))->GetName(), 0);
  mToDiscard = "";
}
