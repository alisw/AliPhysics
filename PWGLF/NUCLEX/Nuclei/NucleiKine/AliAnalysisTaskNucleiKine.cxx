#include "AliAnalysisTaskNucleiKine.h"

#include <TH1D.h>
#include <TH2D.h>
#include <TList.h>
#include <TChain.h>

#include "AliAnalysisManager.h"
#include "AliGenEventHeader.h"
#include "AliHeader.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"


ClassImp(AliAnalysisTaskNucleiKine)

AliAnalysisTaskNucleiKine::AliAnalysisTaskNucleiKine(const string name):
  AliAnalysisTaskSE(name.data()),
  fPmax{0.2f},
  fRmax{1.e4f}, // Disabled for the time being
  fSpinProb{1.f},
  fOutputList{nullptr},
  fEventCounter{nullptr},
  fPtSpectra{nullptr},
  fPdgCodes {211, -211, 321, -321, 2212, -2212, 2112, -2112, 1000010020, -1000010020}

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

  fOutputList->Add(fEventCounter);
  fPtSpectra->Add(fPtSpectra);

  PostData(1,fOutputList);
}

void AliAnalysisTaskNucleiKine::UserExec(Option_t*) {

  AliInputEventHandler* mcHandler = dynamic_cast<AliInputEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!mcHandler) return;

  AliMCEvent* mcEvent = mcHandler->MCEvent();
  if (!mcEvent) return;

  AliGenEventHeader *mcHeader = (mcEvent->GenEventHeader());
  if (!mcHeader) return;

  /// Get the information of all the protons and neutrons to see if they could bind to form a deuteron via coalescence
  /// in this loop also a spatial coordinate is given to all the particles (to be implemented)
  vector<bool> mask(mcEvent->GetNumberOfTracks(),true);
  vector<Particle> protons[2],neutrons[2];
  for (int iTracks = 0; iTracks < mcEvent->GetNumberOfTracks(); ++iTracks) {
    AliVParticle* track = mcEvent->GetTrack(iTracks);
    if (!track) continue;
    const int pdg = track->PdgCode();
    if (pdg == fPdgCodes[kNeutron]) {
      neutrons[0].push_back({{0.f,0.f,0.f,0.f},{track->Px(),track->Py(),track->Pz(),sqrt((0.939565*0.939565)+(track->P()*track->P()))}});
    } else if(pdg == fPdgCodes[kNeutron]) {
      neutrons[1].push_back({{0.f,0.f,0.f,0.f},{track->Px(),track->Py(),track->Pz(),sqrt((0.939565*0.939565)+(track->P()*track->P()))}});
    } else if(pdg == fPdgCodes[kProton]) {
      protons[0].push_back({{0.f,0.f,0.f,0.f},{track->Px(),track->Py(),track->Pz(),sqrt((0.938272*0.938272)+(track->P()*track->P()))}});
    } else if(pdg == fPdgCodes[kAntiProton]) {
      protons[1].push_back({{0.f,0.f,0.f,0.f},{track->Px(),track->Py(),track->Pz(),sqrt((0.938272*0.938272)+(track->P()*track->P()))}});
    }
  }


  for (int iTracks = 0; iTracks < mcEvent->GetNumberOfTracks(); ++iTracks) {
    AliVParticle* track = mcEvent->GetTrack(iTracks);
    if (!track) continue;
    int pdg = track->PdgCode();
    if (mcEvent->IsPhysicalPrimary(iTracks)) {
      for (size_t iParticle = 0; iParticle < fPdgCodes.size(); ++iParticle) {
        if (pdg == fPdgCodes[iParticle] && mask[iTracks])
          fPtSpectra->Fill(iParticle, track->Pt());
      }
    }
  }

  PostData(1, fOutputList);
}

