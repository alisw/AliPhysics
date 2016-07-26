#include "AliAnalysisTaskNucleiKine.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TChain.h>
#include <TRandom3.h>

#include "AliAnalysisManager.h"
#include "AliGenEventHeader.h"
#include "AliHeader.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"

#include "AliGenEventHeaderTunedPbPb.h"

ClassImp(AliAnalysisTaskNucleiKine)

AliAnalysisTaskNucleiKine::AliAnalysisTaskNucleiKine(const string name):
  AliAnalysisTaskSE(name.data()),
  fPotentialShape{"V","-0.2",0,10},
  fSpinProb{1.f},
  fSpatialDistribution{"S","0",0,10},
  fEnableDisplace{false},
  fPdgCodes {211, -211, 321, -321, 2212, -2212, 2112, -2112, 1000010020, -1000010020},
  fOutputList{nullptr},
  fEventCounter{nullptr},
  fPtSpectra{nullptr},
  fCosine{nullptr},
  fPsi2{nullptr},
  fMask{},
  fNeutrons{},
  fProtons{}
{
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskNucleiKine::~AliAnalysisTaskNucleiKine() {
  if(fOutputList) delete fOutputList;
}

void AliAnalysisTaskNucleiKine::UserCreateOutputObjects() {
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  fEventCounter = new TH1D("fEventCounter",";Centrality (%); Events / 1%",100,0,100);
  fPtSpectra = new TH2D("fPtSpectra",";Particle Species;#it{p}_{T} (GeV/#it{c})",fPdgCodes.size(),-0.5,-0.5 + fPdgCodes.size(),100,0.,10.);
  fCosine = new TH2D("fCosine",";Particle Species;cos(2(#phi-#Psi_{2}))",fPdgCodes.size(),-0.5,-0.5 + fPdgCodes.size(),100,0.,10.);
  fPsi2 = new TH1D("fPsi2",";Particle Species;#Psi_{2}",64,0,TMath::TwoPi());

  fOutputList->Add(fEventCounter);
  fOutputList->Add(fPtSpectra);
  fOutputList->Add(fCosine);
  fOutputList->Add(fPsi2);

  PostData(1,fOutputList);
}

void AliAnalysisTaskNucleiKine::UserExec(Option_t*) {
  auto to_0_2pi = [](float &x) { while (x < 0) x+= TMath::TwoPi(); };
  AliInputEventHandler* mcHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!mcHandler) return;

  AliMCEvent* mcEvent = mcHandler->MCEvent();
  if (!mcEvent) return;

  AliGenEventHeaderTunedPbPb *mcHeader = (AliGenEventHeaderTunedPbPb*)mcEvent->GenEventHeader();
  if (!mcHeader) return;
  float psi2 = mcHeader->GetPsi2();
  to_0_2pi(psi2);
  fPsi2->Fill(psi2);

  fEventCounter->Fill(1.);

  /// Clear all the vectors
  fMask.resize(mcEvent->GetNumberOfTracks(),true);
  for (int iS = 0; iS < 2; ++iS) {
    fProtons[iS].clear();
    fNeutrons[iS].clear();
  }

  /// Get the information of all the protons and neutrons to see if they could bind to form a deuteron via coalescence
  /// in this loop also a spatial coordinate is given to all the particles (to be implemented)
  for (int iTracks = 0; iTracks < mcEvent->GetNumberOfTracks(); ++iTracks) {
    AliVParticle* track = mcEvent->GetTrack(iTracks);
    if (!track) continue;
    const int pdg = track->PdgCode();
    double xyz[3]{0.,0.,0.};
    if (fEnableDisplace) {
      gRandom->Sphere(xyz[0],xyz[1],xyz[2],fSpatialDistribution.GetRandom());
    }
    if (pdg == fPdgCodes[kNeutron]) {
      fNeutrons[0].push_back({iTracks,{{xyz[0],xyz[1],xyz[2],0.},{track->Px(),track->Py(),track->Pz(),sqrt((0.939565*0.939565)+(track->P()*track->P()))}}});
    } else if(pdg == fPdgCodes[kNeutron]) {
      fNeutrons[1].push_back({iTracks,{{xyz[0],xyz[1],xyz[2],0.},{track->Px(),track->Py(),track->Pz(),sqrt((0.939565*0.939565)+(track->P()*track->P()))}}});
    } else if(pdg == fPdgCodes[kProton]) {
      fProtons[0].push_back({iTracks,{{xyz[0],xyz[1],xyz[2],0.},{track->Px(),track->Py(),track->Pz(),sqrt((0.938272*0.938272)+(track->P()*track->P()))}}});
    } else if(pdg == fPdgCodes[kAntiProton]) {
      fProtons[1].push_back({iTracks,{{xyz[0],xyz[1],xyz[2],0.},{track->Px(),track->Py(),track->Pz(),sqrt((0.938272*0.938272)+(track->P()*track->P()))}}});
    }
  }

  for (int iS = 0; iS < 2; ++iS) {
    vector<LorentzVector<PxPyPzE4D<double>>> deuterons;
    FirstPartner(iS,deuterons);
    for (auto &d : deuterons) {
      float phi = d.Phi();
      to_0_2pi(phi);
      const float pt = d.Pt();
      fPtSpectra->Fill(8+iS,pt);
      fCosine->Fill(8+iS,pt,cos(2. * (phi - psi2)));
    }
  }

  for (int iTracks = 0; iTracks < mcEvent->GetNumberOfTracks(); ++iTracks) {
    AliVParticle* track = mcEvent->GetTrack(iTracks);
    if (!track) continue;
    if (abs(track->Y()) > 0.5) continue;
    int pdg = track->PdgCode();
    float phi = track->Phi();
    to_0_2pi(phi);
    if (mcEvent->IsPhysicalPrimary(iTracks)) {
      for (size_t iParticle = 0; iParticle < fPdgCodes.size(); ++iParticle) {
        if (pdg == fPdgCodes[iParticle] && fMask[iTracks]) {
          fPtSpectra->Fill(iParticle, track->Pt());
          fCosine->Fill(iParticle,track->Pt(),cos(2. * (phi - psi2)));
        }
      }
    }
  }

  PostData(1, fOutputList);
}

void AliAnalysisTaskNucleiKine::FirstPartner(int iS, vector<LorentzVector<PxPyPzE4D<double>>> &deut) {
  for (auto& p : fProtons[iS]) {
    for (auto& n : fNeutrons[iS]) {
      if (fMask[n.first]) continue;
      auto& p0 = p.second.mom;
      auto& p1 = n.second.mom;
      const double pcm = GetPcm(p0,p1);
      const double r = (n.second.pos - p.second.pos).mag();
      if (fPotentialShape(r) + pcm < 0) {
        const double px = p0.px() + p1.px();
        const double py = p0.py() + p1.py();
        const double pz = p0.pz() + p1.pz();
        deut.push_back({px,py,pz,sqrt((px*px)+(py*py)+(pz*pz) + (1.875613*1.875613))});
        fMask[n.first] = true;
        fMask[p.first] = true;
        break;
      }
    }
  }
}
