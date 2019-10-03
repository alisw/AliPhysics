#include "AliAnalysisTaskLFefficiencies.h"

#include <array>
#include <string>
#include <algorithm>

// ROOT includes
#include <TAxis.h>
#include <TChain.h>
#include <TH3D.h>
#include <TF1.h>
#include <TList.h>
#include <TMath.h>
#include <TParticle.h>
#include <TClonesArray.h>
#include <TTree.h>
#include <TRandom3.h>

// ALIROOT includes
#include "AliAnalysisManager.h"
#include "AliPIDResponse.h"
#include "AliPDG.h"
#include "AliMultSelection.h"
#include "AliTPCPIDResponse.h"
#include "AliTOFPIDResponse.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "AliVEvent.h"
#include "AliVParticle.h"
#include "AliMCEvent.h"
#include "AliInputEventHandler.h"
#include "AliVEventHandler.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "AliAODVertex.h"

using TMath::TwoPi;

const std::string AliAnalysisTaskLFefficiencies::fPosNeg[2] = {"neg","pos"};
const int AliAnalysisTaskLFefficiencies::fNcuts = 8;
const std::string AliAnalysisTaskLFefficiencies::fCutNames[8] = {"FB4","FB5","FB5+PID TPC", "FB5 + TOF matching", "FB5 + PID TOF", "FB5 + TOF matching - TOF mismatch", "FB5 + TOF matching - TOF mismatch + TOF pid", "FB5 + hasTOF + TOF mismatch"};

///\cond CLASSIMP
ClassImp(AliAnalysisTaskLFefficiencies);
///\endcond

AliAnalysisTaskLFefficiencies::AliAnalysisTaskLFefficiencies(TString taskname) :
  AliAnalysisTaskSE(taskname.Data()),
  fEventCut{false},
  fUseMCtruthParams{false}
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

/// Standard destructor
///
AliAnalysisTaskLFefficiencies::~AliAnalysisTaskLFefficiencies(){
  if (AliAnalysisManager::GetAnalysisManager()->IsProofMode()) return;
  if (fOutputList) delete fOutputList;
}

/// This function creates all the histograms and all the objects in general used during the analysis
/// \return void
///
void AliAnalysisTaskLFefficiencies::UserCreateOutputObjects() {

  fOutputList = new TList();
  fOutputList->SetOwner(true);

  fNumberOfRecoPrimaryTracks = new TH1D("fNumberOfRecoPrimaryTracks",";Number of tracks;Number",4001,-0.5,4000.5);
  fOutputList->Add(fNumberOfRecoPrimaryTracks);
  for (int iSpecies = 0; iSpecies < AliPID::kSPECIESC; iSpecies++) {
    for (int iCharge = 0; iCharge < 2; ++iCharge) {
      fGeneratedYPhiPt[iSpecies][iCharge] = new TH3D(Form("Gen_%s_%s",AliPID::ParticleShortName(iSpecies),fPosNeg[iCharge].data()),
        ";y;#varphi;#it{p}_{T} (GeV/#it{c})",9,-0.9,0.9,16,0.,TwoPi(),60,0.,6.);
      fOutputList->Add(fGeneratedYPhiPt[iSpecies][iCharge]);
      fGeneratedEtaPhiPt[iSpecies][iCharge] = new TH3D(Form("GenEta_%s_%s",AliPID::ParticleShortName(iSpecies),fPosNeg[iCharge].data()),
        ";#eta;#varphi;#it{p}_{T} (GeV/#it{c})",10,-1.,1.,16,0.,TwoPi(),60,0.,6.);
      fOutputList->Add(fGeneratedEtaPhiPt[iSpecies][iCharge]);
      for (int iCut = 0; iCut < fNcuts; ++iCut) {
        fReconstructedYPhiPt[iSpecies][iCharge][iCut] = new TH3D(Form("Rec_%s_%s_%i",AliPID::ParticleShortName(iSpecies),fPosNeg[iCharge].data(),iCut),
          Form("%s;y;#varphi;#it{p}_{T} (GeV/#it{c})",fCutNames[iCut].data()),9,-0.9,0.9,16,0.,TwoPi(),60,0.,6.);
        fOutputList->Add(fReconstructedYPhiPt[iSpecies][iCharge][iCut]);
        fReconstructedEtaPhiPt[iSpecies][iCharge][iCut] = new TH3D(Form("RecEta_%s_%s_%i",AliPID::ParticleShortName(iSpecies),fPosNeg[iCharge].data(),iCut),
          Form("%s;#eta;#varphi;#it{p}_{T} (GeV/#it{c})",fCutNames[iCut].data()),10,-1.,1.,16,0.,TwoPi(),60,0.,6.);
        fOutputList->Add(fReconstructedEtaPhiPt[iSpecies][iCharge][iCut]);
      }
      fNsigmaTOFvsPt[iSpecies][iCharge] = new TH2D(Form("nSigmaTOF_%s_%s",AliPID::ParticleShortName(iSpecies),fPosNeg[iCharge].data()),";#it{p}_{T} (GeV/#it{c}); n#sigma_{TOF}",60,0.,6.,1001,-100.1,100.1);
      fOutputList->Add(fNsigmaTOFvsPt[iSpecies][iCharge]);
    }
  }
  fEventCut.AddQAplotsToList(fOutputList);

  PostData(1,fOutputList);
}

/// This is the function that is evaluated for each event. The analysis code stays here.
///
/// \param options Deprecated parameter
/// \return void
///
void AliAnalysisTaskLFefficiencies::UserExec(Option_t *){
  AliVEvent *ev = InputEvent();
  bool EventAccepted = fEventCut.AcceptEvent(ev);

  if (!EventAccepted) {
    PostData(1, fOutputList);
    return;
  }

  /// To perform the majority of the analysis - and also this one - the standard PID handler is
  /// required.
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* handl = (AliInputEventHandler*)mgr->GetInputEventHandler();
  AliPIDResponse* pid = handl->GetPIDResponse();
  if (!pid) {
    ::Fatal("AliAnalysisTaskLFefficiencies::UserExec","Missing PID response. Did you attach the AliPIDresponseTask to your analysis?");
  }

  AliMCEvent* mcEv = MCEvent();
  if (!mcEv)
    ::Fatal("AliAnalysisTaskLFefficiencies::UserExec","MC analysis requested on a sample without the MC particle array.");

  for (int iMC = 0; iMC < mcEv->GetNumberOfTracks(); ++iMC) {
    AliVParticle* part = mcEv->GetTrack(iMC);
    if (!part->IsPhysicalPrimary()) continue;
    const int pdg = std::abs(part->PdgCode());
    const int iCharge = part->Charge() > 0 ? 1 : 0;
    for (int iSpecies = 0; iSpecies < AliPID::kSPECIESC; ++iSpecies) {
      if (pdg == AliPID::ParticleCode(iSpecies)) {
        fGeneratedYPhiPt[iSpecies][iCharge]->Fill(part->Y(), part->Phi(), part->Pt());
        fGeneratedEtaPhiPt[iSpecies][iCharge]->Fill(part->Eta(), part->Phi(), part->Pt());
        break;
      }
    }
  }


  /// Checking how many deuterons in acceptance are reconstructed well
  TLorentzVector v;
  int nPrimaries{0};
  for (int iT = 0; iT < (int)ev->GetNumberOfTracks(); ++iT) {
    /// Get the track and do the minimal cuts
    AliAODTrack *track = dynamic_cast<AliAODTrack*>(ev->GetTrack(iT));
    if (track->GetID() < 0) continue;
    if (!track->TestFilterBit(BIT(4))) continue;

    AliVParticle *part = (AliAODMCParticle*)mcEv->GetTrack(TMath::Abs(track->GetLabel()));
    if (!part) continue;
    if (!part->IsPhysicalPrimary()) continue;
    nPrimaries++;
    const int iCharge = part->Charge() > 0 ? 1 : 0;
    int iSpecies = -1;
    for (int iS = 0; iS < AliPID::kSPECIESC; ++iS) {
      if (std::abs(part->PdgCode()) == AliPID::ParticleCode(iS)) {
        iSpecies = iS;
        break;
      }
    }
    if (iSpecies < 0) continue;

    const double pt = fUseMCtruthParams ? part->Pt() : track->Pt() * AliPID::ParticleCharge(iSpecies);
    const double eta = fUseMCtruthParams ? part->Eta() : track->Eta();
    const double phi = fUseMCtruthParams ? part->Phi() : track->Phi();
    v.SetPtEtaPhiM(pt, eta, phi, AliPID::ParticleMass(iSpecies));

    bool hasFB5 = track->TestFilterBit(BIT(5));
    bool TPCpid = std::abs(pid->NumberOfSigmasTPC(track, static_cast<AliPID::EParticleType>(iSpecies))) < 3;
    bool hasTOF = HasTOF(track);
    double nSigmaTOF = pid->NumberOfSigmasTOF(track, static_cast<AliPID::EParticleType>(iSpecies));
    bool TOFpid = std::abs(nSigmaTOF) < 3;
    int TOFlabels[3];
    track->GetTOFLabel(TOFlabels);
    bool TOFmismatch = TOFlabels[0] != TMath::Abs(track->GetLabel());
    bool cuts[fNcuts] = {true, hasFB5, hasFB5 && TPCpid, hasFB5 && hasTOF, hasFB5 && hasTOF && TOFpid,
      hasFB5 && hasTOF && !TOFmismatch, hasFB5 && hasTOF && !TOFmismatch && TOFpid, hasFB5 && hasTOF && TOFmismatch};

    for (int iCut = 0; iCut < fNcuts; ++iCut) {
      if (cuts[iCut]) {
        fReconstructedYPhiPt[iSpecies][iCharge][iCut]->Fill(v.Rapidity(),phi,pt);
        fReconstructedEtaPhiPt[iSpecies][iCharge][iCut]->Fill(eta,phi,pt);
        if (iCut==5) {
          fNsigmaTOFvsPt[iSpecies][iCharge]->Fill(pt,nSigmaTOF);
        }
      }
    }
  } // End AOD track loop
  fNumberOfRecoPrimaryTracks->Fill(nPrimaries);

  //  Post output data.
  PostData(1,fOutputList);
}

/// Merge the output. Called once at the end of the query.
///
/// \return void
///
void AliAnalysisTaskLFefficiencies::Terminate(Option_t *) {
  return;
}

/// This function checks whether a track has or has not a prolongation in TOF.
///
/// \param track Track that has to be checked
///
bool AliAnalysisTaskLFefficiencies::HasTOF(AliVTrack *track) {
  const bool hasTOFout  = track->GetStatus() & AliVTrack::kTOFout;
  const bool hasTOFtime = track->GetStatus() & AliVTrack::kTIME;
  const bool hasGoodLength = track->GetIntegratedLength() > 350.;
  return hasTOFout && hasTOFtime && hasGoodLength;
}
