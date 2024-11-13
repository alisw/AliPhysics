class TParticle;

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TMath.h"
#include "TParticlePDG.h"
#include <Riostream.h>

#include "AliAnalysisSigmaBarChargedGen.h"

ClassImp(AliAnalysisSigmaBarChargedGen)

    AliAnalysisSigmaBarChargedGen::AliAnalysisSigmaBarChargedGen()
    : AliAnalysisTaskSE(), fHistos_misc(nullptr) {
  // default constructor
}

AliAnalysisSigmaBarChargedGen::AliAnalysisSigmaBarChargedGen(
    const char *name, TString lExtraOptions)
    : AliAnalysisTaskSE(name), fHistos_misc(nullptr) {

  // Standard output
  DefineOutput(1, TList::Class()); // Miscellaneous Histograms
}

AliAnalysisSigmaBarChargedGen::~AliAnalysisSigmaBarChargedGen() {
  //------------------------------------------------
  // DESTRUCTOR
  //------------------------------------------------
  if (fHistos_misc) {
    delete fHistos_misc;
    fHistos_misc = 0x0;
  }
}

//________________________________________________________________________
void AliAnalysisSigmaBarChargedGen::UserCreateOutputObjects() {

  // miscellaneous histograms
  fHistos_misc = new THashList();
  fHistos_misc->SetName("Histos_misc");
  fHistos_misc->Add(new TH1F("hSelEvents", "", 14, 0, 14)); // #events histo

  // histograms spec
  fHistos_misc->Add(new TH1F("MC_AntiSigmaPlus_1", "", 800, 0, 20)); // vs nch
  fHistos_misc->Add(new TH1F("MC_AntiSigmaPlus_2", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_AntiSigmaMinus_1", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_AntiSigmaMinus_2", "", 800, 0, 20));

  fHistos_misc->Add(new TH1F("MC_Nbar_1", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_Nbar_prim_1", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_N_1", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_N_prim_1", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_Pbar_1", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_Pbar_prim_1", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_P_1", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_P_prim_1", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_Pion_1", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_Pion_prim_1", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_Kaon_1", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_Kaon_prim_1", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_Lambda_1", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_Lambda_prim_1", "", 800, 0, 20));

  fHistos_misc->Add(new TH1F("MC_Nbar_2", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_Nbar_prim_2", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_N_2", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_N_prim_2", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_Pbar_2", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_Pbar_prim_2", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_P_2", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_P_prim_2", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_Pion_2", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_Pion_prim_2", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_Kaon_2", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_Kaon_prim_2", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_Lambda_2", "", 800, 0, 20));
  fHistos_misc->Add(new TH1F("MC_Lambda_prim_2", "", 800, 0, 20));

  // Output posting
  PostData(1, fHistos_misc);

} // end UserCreateOutputObjects

//________________________________________________________________________
void AliAnalysisSigmaBarChargedGen::UserExec(Option_t *) {
  dynamic_cast<TH1F *>(fHistos_misc->FindObject("hSelEvents"))->Fill(1);
  // get MC event and fill histo
  AliMCEvent *lMCev = MCEvent();
  if (!lMCev) {
    Printf("ERROR: Could not retrieve MC event in file %s\n",
           fInputHandler->GetTree()->GetCurrentFile()->GetName());
    PostData(1, fHistos_misc);
    return;
  }
  dynamic_cast<TH1F *>(fHistos_misc->FindObject("hSelEvents"))->Fill(2);
  //
  AliStack *lMCstack = lMCev->Stack();
  if (!lMCstack) {
    Printf("ERROR: Could not retrieve MC stack \n");
    return;
  }

  dynamic_cast<TH1F *>(fHistos_misc->FindObject("hSelEvents"))->Fill(3);

  Bool_t fIsEvtINELgtZERO = false;

  // Loop over all primary MC particle
  for (int itrk = 0; itrk < lMCstack->GetNtrack(); itrk++) {
    TParticle *lPart = lMCstack->Particle(itrk);
    // selected charged primary particles
    if (lPart && lMCstack->IsPhysicalPrimary(itrk) &&
        ((TParticlePDG *)lPart->GetPDG())->Charge() != 0) {
      // check if event is INEL>0
      if (!fIsEvtINELgtZERO) {
        if (std::abs(lPart->Eta()) < 1) {
          fIsEvtINELgtZERO = true;
        }
      }
    }
  }

  if (!fIsEvtINELgtZERO) {
    return;
  }

  dynamic_cast<TH1F *>(fHistos_misc->FindObject("hSelEvents"))->Fill(10);

  for (int itrk = 0; itrk < lMCstack->GetNtrack(); itrk++) {
    TParticle *lPart = lMCstack->Particle(itrk);
    int pdg = (int)lPart->GetPdgCode();
    if (pdg == -3222) {
      if (TMath::Abs(lPart->Y()) < 0.5) {
        dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_AntiSigmaMinus_1"))
            ->Fill(lPart->Pt());
      }
      if (lPart->Y() > -0.96 && lPart->Y() < 0.04) {
        dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_AntiSigmaMinus_2"))
            ->Fill(lPart->Pt());
      }
    }
    if (pdg == -3112) {
      if (TMath::Abs(lPart->Y()) < 0.5) {
        dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_AntiSigmaPlus_1"))
            ->Fill(lPart->Pt());
      }
      if (lPart->Y() > -0.96 && lPart->Y() < 0.04) {
        dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_AntiSigmaPlus_2"))
            ->Fill(lPart->Pt());
      }
    }
    if (pdg == -2112) {
      if (TMath::Abs(lPart->Y()) < 0.5) {
        dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_Nbar_1"))
            ->Fill(lPart->Pt());
      }
      if (lPart->Y() > -0.96 && lPart->Y() < 0.04) {
        dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_Nbar_2"))
            ->Fill(lPart->Pt());
      }

      if (lMCstack->IsPhysicalPrimary(itrk)) {
        if (TMath::Abs(lPart->Y()) < 0.5) {
          dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_Nbar_prim_1"))
              ->Fill(lPart->Pt());
        }
        if (lPart->Y() > -0.96 && lPart->Y() < 0.04) {
          dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_Nbar_prim_2"))
              ->Fill(lPart->Pt());
        }
      }
    }
    if (pdg == 2112) {
      if (TMath::Abs(lPart->Y()) < 0.5) {
        dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_N_1"))
            ->Fill(lPart->Pt());
      }
      if (lPart->Y() > -0.96 && lPart->Y() < 0.04) {
        dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_N_2"))
            ->Fill(lPart->Pt());
      }

      if (lMCstack->IsPhysicalPrimary(itrk)) {
        if (TMath::Abs(lPart->Y()) < 0.5) {
          dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_N_prim_1"))
              ->Fill(lPart->Pt());
        }
        if (lPart->Y() > -0.96 && lPart->Y() < 0.04) {
          dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_N_prim_2"))
              ->Fill(lPart->Pt());
        }
      }
    }
    if (pdg == -2212) {
      if (TMath::Abs(lPart->Y()) < 0.5) {
        dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_Pbar_1"))
            ->Fill(lPart->Pt());
      }
      if (lPart->Y() > -0.96 && lPart->Y() < 0.04) {
        dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_Pbar_2"))
            ->Fill(lPart->Pt());
      }

      if (lMCstack->IsPhysicalPrimary(itrk)) {
        if (TMath::Abs(lPart->Y()) < 0.5) {
          dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_Pbar_prim_1"))
              ->Fill(lPart->Pt());
        }
        if (lPart->Y() > -0.96 && lPart->Y() < 0.04) {
          dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_Pbar_prim_2"))
              ->Fill(lPart->Pt());
        }
      }
    }
    if (pdg == 2212) {
      if (TMath::Abs(lPart->Y()) < 0.5) {
        dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_P_1"))
            ->Fill(lPart->Pt());
      }
      if (lPart->Y() > -0.96 && lPart->Y() < 0.04) {
        dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_P_2"))
            ->Fill(lPart->Pt());
      }

      if (lMCstack->IsPhysicalPrimary(itrk)) {
        if (TMath::Abs(lPart->Y()) < 0.5) {
          dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_P_prim_1"))
              ->Fill(lPart->Pt());
        }
        if (lPart->Y() > -0.96 && lPart->Y() < 0.04) {
          dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_P_prim_2"))
              ->Fill(lPart->Pt());
        }
      }
    }
    if (pdg == 211 || pdg == -211) {
      if (TMath::Abs(lPart->Y()) < 0.5) {
        dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_Pion_1"))
            ->Fill(lPart->Pt());
      }
      if (lPart->Y() > -0.96 && lPart->Y() < 0.04) {
        dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_Pion_2"))
            ->Fill(lPart->Pt());
      }

      if (lMCstack->IsPhysicalPrimary(itrk)) {
        if (TMath::Abs(lPart->Y()) < 0.5) {
          dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_Pion_prim_1"))
              ->Fill(lPart->Pt());
        }
        if (lPart->Y() > -0.96 && lPart->Y() < 0.04) {
          dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_Pion_prim_2"))
              ->Fill(lPart->Pt());
        }
      }
    }
    if (pdg == 321 || pdg == -321) {
      if (TMath::Abs(lPart->Y()) < 0.5) {
        dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_Kaon_1"))
            ->Fill(lPart->Pt());
      }
      if (lPart->Y() > -0.96 && lPart->Y() < 0.04) {
        dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_Kaon_2"))
            ->Fill(lPart->Pt());
      }

      if (lMCstack->IsPhysicalPrimary(itrk)) {
        if (TMath::Abs(lPart->Y()) < 0.5) {
          dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_Kaon_prim_1"))
              ->Fill(lPart->Pt());
        }
        if (lPart->Y() > -0.96 && lPart->Y() < 0.04) {
          dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_Kaon_prim_2"))
              ->Fill(lPart->Pt());
        }
      }
    }
    if (pdg == 3122 || pdg == -3122) {
      if (TMath::Abs(lPart->Y()) < 0.5) {
        dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_Lambda_1"))
            ->Fill(lPart->Pt());
      }
      if (lPart->Y() > -0.96 && lPart->Y() < 0.04) {
        dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_Lambda_2"))
            ->Fill(lPart->Pt());
      }

      if (lMCstack->IsPhysicalPrimary(itrk)) {
        if (TMath::Abs(lPart->Y()) < 0.5) {
          dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_Lambda_prim_1"))
              ->Fill(lPart->Pt());
        }
        if (lPart->Y() > -0.96 && lPart->Y() < 0.04) {
          dynamic_cast<TH1F *>(fHistos_misc->FindObject("MC_Lambda_prim_2"))
              ->Fill(lPart->Pt());
        }
      }
    }
  }

  PostData(1, fHistos_misc);
}

//________________________________________________________________________
void AliAnalysisSigmaBarChargedGen::Terminate(Option_t *) {}
