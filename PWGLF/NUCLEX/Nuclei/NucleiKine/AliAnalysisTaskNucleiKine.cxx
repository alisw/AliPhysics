#include "AliAnalysisTaskNucleiKine.h"

#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TList.h>
#include <TMath.h>

#include "AliAnalysisManager.h"
#include "AliCollisionGeometry.h"
#include "AliGenEventHeader.h"
#include "AliHeader.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"

ClassImp(AliAnalysisTaskNucleiKine);

AliAnalysisTaskNucleiKine::AliAnalysisTaskNucleiKine(const char* name) :
  AliAnalysisTaskSE{name},
  fAfterburner{},
  fUseAfterburner{false},
  fPdgCodes{211, -211, 321, -321, 2212, -2212, 2112, -2112, 1000010020, -1000010020,3122,-3122,3312,-3312,3334,-3334},
  fParticleNames{"#pi^{+}", "#pi^{-}", "K^{+}", "K^{-}", "p", "#bar{p}", "n", "#bar{n}", "d", "#bar{d}",
      "#Lambda", "#bar{#Lambda}", "#Xi^{+}", "#Xi^{-}", "#Omega^{+}", "#Omega^{-}"},
  fIgnoreCentrality{true},
  fOutputList{nullptr},
  fEventCounter{nullptr},
  fPtSpectra{nullptr},
  fPtSpectraNoRapidity{nullptr}
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

  double b[] = {0.00,3.50,4.94,6.98,8.55,9.88,11.04};
  fEventCounter = new TH1D("fEventCounter",";Impact parameter (fm); Events",6,b);
  fReactionPlane = new TH1D("fReactionPlane",";#Psi (rad); Events",200,-TMath::PiOver2(),TMath::PiOver2());
  fPtSpectra = new TH3D("fPtSpectra",";Centrality bin;Particle Species;#it{p}_{T} (GeV/#it{c})",6,0.5,6.5,fPdgCodes.size(),-0.5,-0.5 + fPdgCodes.size(),100,0.,10.);
  fPtSpectraNoRapidity = new TH3D("fPtSpectraNoRapidity",";Centrality bin;Particle Species;#it{p}_{T} (GeV/#it{c})",6,0.5,6.5,fPdgCodes.size(),-0.5,-0.5 + fPdgCodes.size(),100,0.,10.);
  for (int iB = 1; iB <= fParticleNames.size(); ++iB) {
    fPtSpectra->GetYaxis()->SetBinLabel(iB,fParticleNames[iB-1].data());
    fPtSpectraNoRapidity->GetYaxis()->SetBinLabel(iB,fParticleNames[iB-1].data());
  }

  fOutputList->Add(fEventCounter);
  fOutputList->Add(fPtSpectra);
  fOutputList->Add(fPtSpectraNoRapidity);

  PostData(1,fOutputList);
}

void AliAnalysisTaskNucleiKine::UserExec(Option_t*) {
  AliMCEvent* mcEvent = MCEvent();
  if (!mcEvent)
    AliFatal("Missing MC event!");


  AliStack* stack = mcEvent->Stack();
  if (!stack)
    AliFatal("Missing stack.");

  if (fUseAfterburner) {
    fAfterburner.SetStack(stack);
    fAfterburner.Generate();
  }
  int nstack = stack->GetNtrack();

  float impact = 1.f;
  if (!fIgnoreCentrality) {
    AliCollisionGeometry* hd = dynamic_cast<AliCollisionGeometry*>(mcEvent->GenEventHeader());
    if (!hd)
      AliFatal("Missing collision geometry");
    float impact = hd->ImpactParameter();
  }
  int impact_bin = fEventCounter->GetXaxis()->FindBin(impact);
  fEventCounter->Fill(impact);

  for (int iTracks = 0; iTracks < nstack; ++iTracks) {
    TParticle* track = stack->Particle(iTracks);
    if (!track) continue;
    if (!stack->IsPhysicalPrimary(iTracks)) continue;
    const int pdg = track->GetPdgCode();
    for (size_t iS = 0; iS < fPdgCodes.size(); ++iS) {
      if (pdg == fPdgCodes[iS]) {
        if (TMath::Abs(track->Y()) <= 0.5) {
          fPtSpectra->Fill(impact_bin,iS, track->Pt());
        }
        fPtSpectraNoRapidity->Fill(impact_bin,iS, track->Pt());
        break;
      }
    }
  }
  PostData(1, fOutputList);
}

