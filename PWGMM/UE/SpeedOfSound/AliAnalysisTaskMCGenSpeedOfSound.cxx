/**************************************************************************
 * Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
/*
 *  AliAnalysisTaskMCGenSpeedOfSound.cxx
 */

//_____ ROOT headers
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TF1.h>
#include <TFile.h>
#include <THnSparse.h>
#include <TList.h>
#include <TMath.h>
#include <TRandom.h>
#include <TTree.h>
#include <TTreeStream.h>

#include "TChain.h"
#include "TObjArray.h"

//_____ ALIROOT headers
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAODInputHandler.h"
#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisFilter.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliGenPythiaEventHeader.h"
#include "AliHeader.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "AliVParticle.h"

//_____ Additional includes
#include "AliGenEventHeader.h"
#include "AliVEvent.h"
// #include "AliAnalysisUtils.h"

//_____ AnalysisTask headers
#include "AliAnalysisTaskMCGenSpeedOfSound.h"
#include "AliAnalysisTaskSE.h"

//_____ STL includes
#include <iostream>

using namespace std;

ClassImp(AliAnalysisTaskMCGenSpeedOfSound)

    //_____________________________________________________________________________

    AliAnalysisTaskMCGenSpeedOfSound::AliAnalysisTaskMCGenSpeedOfSound()
    : AliAnalysisTaskSE(),
      fMcEvent(0x0),
      fMcHandler(0x0),
      fStack(0),
      hNch08(0),
      hNchV0M(0),
      hNchEtaNeg(0),
      hNchEtaPos(0),
      hNchTPCEtaGap(0),
      hNchSPDEtaGap(0),
      hPtvsNch08(0),
      hPtvsV0M(0),
      hPtvsNchEtaPos(0),
      hPtvsNchEtaNeg(0),
      hPtvsNchTPCEtaGap(0),
      hPtvsNchSPDEtaGap(0),
      fListOfObjects(0) {}

//______________________________________________________________________________

AliAnalysisTaskMCGenSpeedOfSound::AliAnalysisTaskMCGenSpeedOfSound(
    const char *name)
    : AliAnalysisTaskSE(name),
      fMcEvent(0x0),
      fMcHandler(0x0),
      fStack(0),
      hNch08(0),
      hNchV0M(0),
      hNchEtaNeg(0),
      hNchEtaPos(0),
      hNchTPCEtaGap(0),
      hNchSPDEtaGap(0),
      hPtvsNch08(0),
      hPtvsV0M(0),
      hPtvsNchEtaPos(0),
      hPtvsNchEtaNeg(0),
      hPtvsNchTPCEtaGap(0),
      hPtvsNchSPDEtaGap(0),
      fListOfObjects(0) {
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());  // Basic output slot
}

//_____________________________________________________________________________

AliAnalysisTaskMCGenSpeedOfSound::~AliAnalysisTaskMCGenSpeedOfSound() {
  // Destructor
  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor
  if (fListOfObjects) {
    delete fListOfObjects;
    fListOfObjects = 0x0;
  }
}

//______________________________________________________________________________

void AliAnalysisTaskMCGenSpeedOfSound::UserCreateOutputObjects() {
  // ### Analysis output
  fListOfObjects = new TList();
  fListOfObjects->SetOwner(kTRUE);

  constexpr int pt_Nbins{200};
  double pt_bins[pt_Nbins + 1] = {0};
  for (int i = 0; i <= pt_Nbins; ++i) {
    pt_bins[i] = i * 0.05;
  }

  constexpr int nch_Nbins{6000};
  double nch_bins[nch_Nbins + 1] = {0};
  for (int i = 0; i <= nch_Nbins; ++i) {
    nch_bins[i] = -0.5 + (double)i;
  }

  constexpr double v0mAmp_width{15.0};
  constexpr int v0mAmp_Nbins{1400};
  double v0mAmp_bins[v0mAmp_Nbins + 1] = {0};
  for (int i = 0; i <= v0mAmp_Nbins; ++i) {
    v0mAmp_bins[i] = 0.0 + i * v0mAmp_width;
  }

  constexpr int spd_Nbins{6000};
  double spd_bins[spd_Nbins + 1] = {0};
  for (int i = 0; i <= spd_Nbins; ++i) {
    spd_bins[i] = -0.5 + (double)i;
  }

  hNch08 = new TH1D("hNch08", "; #it{N}_{ch} (|#eta|<0.8); Entries", nch_Nbins,
                    nch_bins);
  fListOfObjects->Add(hNch08);

  hNchV0M = new TH1D("hNchV0M", ";#it{N}_{ch} in the V0M; Entries",
                     v0mAmp_Nbins, v0mAmp_bins);
  fListOfObjects->Add(hNchV0M);

  hNchEtaNeg =
      new TH1D("hNchEtaNeg", "; #it{N}_{ch} (-0.8#leq#eta#leq0); Entries ",
               nch_Nbins, nch_bins);
  fListOfObjects->Add(hNchEtaNeg);

  hNchEtaPos =
      new TH1D("hNchEtaPos", "; #it{N}_{ch} (0#leq#eta#leq0.8); Entries ",
               nch_Nbins, nch_bins);
  fListOfObjects->Add(hNchEtaPos);

  hNchTPCEtaGap = new TH1D("hNchTPCEtaGap",
                           "; #it{N}_{ch} (0.5#leq|#eta|#leq0.8); Entries ",
                           nch_Nbins, nch_bins);
  fListOfObjects->Add(hNchTPCEtaGap);

  hNchSPDEtaGap = new TH1D("hNchSPDEtaGap",
                           "; #it{N}_{ch} (0.7#leq|#eta|#leq1.4); Entries ",
                           spd_Nbins, spd_bins);
  fListOfObjects->Add(hNchSPDEtaGap);

  hPtvsNch08 = new TH2D(
      "hPtvsNch08",
      "; #it{N}_{ch} (|#eta|<0.8); #it{p}_{T} (|#eta|#leq0.8, GeV/#it{c})",
      nch_Nbins, nch_bins, pt_Nbins, pt_bins);
  fListOfObjects->Add(hPtvsNch08);

  hPtvsV0M = new TH2D(
      "hPtvsV0M",
      "; #it{N}_{ch} in the V0M; #it{p}_{T} (|#eta|#leq0.8, GeV/#it{c})",
      v0mAmp_Nbins, v0mAmp_bins, pt_Nbins, pt_bins);
  fListOfObjects->Add(hPtvsV0M);

  hPtvsNchEtaPos = new TH2D("hPtvsNchEtaPos",
                            "; #it{N}_{ch} (0#leq#eta#leq0.8); #it{p}_{T} "
                            "(-0.8#leq#eta#leq0, GeV/#it{c})",
                            nch_Nbins, nch_bins, pt_Nbins, pt_bins);
  fListOfObjects->Add(hPtvsNchEtaPos);

  hPtvsNchEtaNeg = new TH2D("hPtvsNchEtaNeg",
                            "; #it{N}_{ch} (-0.8#leq#eta#leq0); #it{p}_{T} "
                            "(0#leq#eta#leq0.8, GeV/#it{c})",
                            nch_Nbins, nch_bins, pt_Nbins, pt_bins);
  fListOfObjects->Add(hPtvsNchEtaNeg);

  hPtvsNchTPCEtaGap =
      new TH2D("hPtvsNchTPCEtaGap",
               "; #it{N}_{ch} (0.5#leq|#eta|#leq0.8); #it{p}_{T} "
               "(|#eta|#leq0.3, GeV/#it{c})",
               nch_Nbins, nch_bins, pt_Nbins, pt_bins);
  fListOfObjects->Add(hPtvsNchTPCEtaGap);

  hPtvsNchSPDEtaGap =
      new TH2D("hPtvsNchSPDEtaGap",
               "; #it{N}_{ch} (0.7#leq|#eta|#leq1.4); #it{p}_{T} "
               "(|#eta#|leq0.4, GeV/#it{c})",
               spd_Nbins, spd_bins, pt_Nbins, pt_bins);
  fListOfObjects->Add(hPtvsNchSPDEtaGap);

  // ### List of outputs
  PostData(1, fListOfObjects);
}
//______________________________________________________________________________
void AliAnalysisTaskMCGenSpeedOfSound::Init() {
  //
  fMcHandler = dynamic_cast<AliInputEventHandler *>(
      AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
}
//______________________________________________________________________________
void AliAnalysisTaskMCGenSpeedOfSound::UserExec(Option_t *) {
  // ### Initialize
  Init();

  // ### MC handler
  if (fMcHandler)
    fMcEvent = fMcHandler->MCEvent();
  else {
    if (fDebug > 1)
      printf("AliAnalysisTaskMCGenSpeedOfSound::Handler() fMcHandler = NULL\n");
    return;
  }

  // ### MC event
  if (!fMcEvent) {
    if (fDebug > 1)
      printf("AliAnalysisTaskMCGenSpeedOfSound::UserExec() fMcEvent = NULL \n");
    return;
  }

  fStack = ((AliMCEvent *)fMcEvent)->Stack();
  if (!fStack) {
    Printf("ERROR: Could not retrieve MC stack \n");
    cout << "Name of the file with pb :"
         << fInputHandler->GetTree()->GetCurrentFile()->GetName() << endl;
    return;
  }

  // ### MC event selection
  Bool_t isEventMCSelected = IsMCEventSelected(fMcEvent);
  if (!isEventMCSelected) return;

  GetMultipliciy();

  PostData(1, fListOfObjects);

  return;
}
//______________________________________________________________________________
Bool_t AliAnalysisTaskMCGenSpeedOfSound::IsMCEventSelected(TObject *obj) {
  Bool_t isSelected = kTRUE;

  AliMCEvent *event = 0x0;
  event = dynamic_cast<AliMCEvent *>(obj);
  if (!event) isSelected = kFALSE;

  return isSelected;
}

//______________________________________________________________________
void AliAnalysisTaskMCGenSpeedOfSound::GetMultipliciy() {
  bool isPhysPrim{false};
  double qPart{0.0};
  double etaPart{999.0};
  double ptPart{999.0};

  int eta08{0};
  int v0m{0};
  int etapos{0};
  int etaneg{0};
  int tpcetagap{0};
  int spdetagap{0};

  // ### particle loop
  for (Int_t ipart = 0; ipart < fMcEvent->GetNumberOfTracks(); ++ipart) {
    TParticle *mcPart = 0x0;
    mcPart = (TParticle *)fMcEvent->Particle(ipart);
    if (!mcPart) continue;
    // selection of primary charged particles
    if (!(mcPart->GetPDG())) continue;
    qPart = mcPart->GetPDG()->Charge() / 3.;
    if (TMath::Abs(qPart) < 0.001) continue;
    isPhysPrim = fMcEvent->IsPhysicalPrimary(ipart);
    if (!isPhysPrim) continue;

    etaPart = mcPart->Eta();
    if (TMath::Abs(etaPart) <= 0.8) {
      eta08++;
    }
    if ((2.8 < etaPart && etaPart < 5.1) ||
        (-3.7 < etaPart && etaPart < -1.7)) {
      v0m++;
    }
    if (etaPart >= -0.8 && etaPart < 0.0) {
      etaneg++;
    }
    if (etaPart >= 0.0 && etaPart <= 0.8) {
      etapos++;
    }
    if (TMath::Abs(etaPart) >= 0.5 && TMath::Abs(etaPart) <= 0.8) {
      tpcetagap++;
    }
    if (TMath::Abs(etaPart) >= 0.7 && TMath::Abs(etaPart) <= 1.4) {
      spdetagap++;
    }
  }  // particle loop

  // ### particle loop 2
  for (Int_t ipart = 0; ipart < fMcEvent->GetNumberOfTracks(); ++ipart) {
    TParticle *mcPart = 0x0;
    mcPart = (TParticle *)fMcEvent->Particle(ipart);
    if (!mcPart) continue;
    // selection of primary charged particles
    if (!(mcPart->GetPDG())) continue;
    qPart = mcPart->GetPDG()->Charge() / 3.;
    if (TMath::Abs(qPart) < 0.001) continue;
    isPhysPrim = fMcEvent->IsPhysicalPrimary(ipart);
    if (!isPhysPrim) continue;

    etaPart = mcPart->Eta();
    ptPart = mcPart->Pt();
    if (TMath::Abs(etaPart) <= 0.8) {
      hPtvsNch08->Fill(eta08, ptPart);
      hPtvsV0M->Fill(v0m, ptPart);
    }
    if (etaPart >= -0.8 && etaPart < 0.0) {
      hPtvsNchEtaPos->Fill(etapos, ptPart);
    }
    if (etaPart >= 0.0 && etaPart <= 0.8) {
      hPtvsNchEtaNeg->Fill(etaneg, ptPart);
    }
    if (TMath::Abs(etaPart) <= 0.3) {
      hPtvsNchTPCEtaGap->Fill(tpcetagap, ptPart);
    }
    if (TMath::Abs(etaPart) <= 0.4) {
      hPtvsNchSPDEtaGap->Fill(spdetagap, ptPart);
    }
  }  // particle loop

  hNch08->Fill(eta08);
  hNchV0M->Fill(v0m);
  hNchEtaNeg->Fill(etaneg);
  hNchEtaPos->Fill(etapos);
  hNchTPCEtaGap->Fill(tpcetagap);
  hNchSPDEtaGap->Fill(spdetagap);
}

//______________________________________________________________________________

void AliAnalysisTaskMCGenSpeedOfSound::Terminate(Option_t *) {
  fListOfObjects = dynamic_cast<TList *>(GetOutputData(1));
  if (!fListOfObjects) {
    Printf("ERROR: Output list not available");
    return;
  }

  return;
}

//_____________________________________________________________________________

Int_t AliAnalysisTaskMCGenSpeedOfSound::GetPidCode(Int_t pdgCode) const {
  Int_t pidCode = 999;

  switch (TMath::Abs(pdgCode)) {
    case 211:
      pidCode = 0;  // pion
      break;
    case 321:
      pidCode = 1;  // kaon
      break;
    case 2212:
      pidCode = 2;  // proton
      break;
    case 310:
      pidCode = 3;  // K0s
      break;
    case 3122:
      pidCode = 4;  // Lambda
      break;
    case 3312:
      pidCode = 5;  // Xi-
      break;
    case 3334:
      pidCode = 6;  // Omega-
      break;
    case 333:
      pidCode = 7;  // phi(1020)
      break;
    case 313:
      pidCode = 8;  // K*(892)0
      break;
    case 323:
      pidCode = 9;  // K*(892) +-
      break;
    case 3212:
      pidCode = 10;  // Sigma 0
      break;
    default:
      break;
  };

  return pidCode;
}
