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
#include "AliVVertex.h"

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
      hNchFull(0),
      hNchV0M(0),
      hTPCEtaGap(0),
      hSPDEtaGap(0),
      hNchEtaGap(0),
      hSPDFull(0),
      hSPDEtaAdj(0),
      hSPDEtaGapW(0),
      hTPCFull(0),
      hEtFull(0),
      hEtEtaGap(0),
      hPtvsNchFull(0),
      hPtvsV0M(0),
      hPtvsTPCEtaGap(0),
      hPtvsSPDEtaGap(0),
      hPtvsNchEtaGap(0),
      hPtvsSPDFull(0),
      hPtvsSPDEtaAdj(0),
      hPtvsSPDEtaGapW(0),
      hPtvsTPCFull(0),
      hPtvsEtFull(0),
      hPtvsEtEtaGap(0),
      fListOfObjects(0) {}

//______________________________________________________________________________

AliAnalysisTaskMCGenSpeedOfSound::AliAnalysisTaskMCGenSpeedOfSound(
    const char *name)
    : AliAnalysisTaskSE(name),
      fMcEvent(0x0),
      fMcHandler(0x0),
      fStack(0),
      hNchFull(0),
      hNchV0M(0),
      hTPCEtaGap(0),
      hSPDEtaGap(0),
      hNchEtaGap(0),
      hSPDFull(0),
      hSPDEtaAdj(0),
      hSPDEtaGapW(0),
      hTPCFull(0),
      hEtFull(0),
      hEtEtaGap(0),
      hPtvsNchFull(0),
      hPtvsV0M(0),
      hPtvsTPCEtaGap(0),
      hPtvsSPDEtaGap(0),
      hPtvsNchEtaGap(0),
      hPtvsSPDFull(0),
      hPtvsSPDEtaAdj(0),
      hPtvsSPDEtaGapW(0),
      hPtvsTPCFull(0),
      hPtvsEtFull(0),
      hPtvsEtEtaGap(0),
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

  constexpr int nchFull_Nbins{4500};
  double nchFull_bins[nchFull_Nbins + 1] = {0};
  for (int i = 0; i <= nchFull_Nbins; ++i) {
    nchFull_bins[i] = 0.0 + (2.0 * i);
  }

  constexpr int nchEtaGap_Nbins{4000};
  double nchEtaGap_bins[nchEtaGap_Nbins + 1] = {0};
  for (int i = 0; i <= nchEtaGap_Nbins; ++i) {
    nchEtaGap_bins[i] = 0.0 + (2.0 * i);
  }

  constexpr double v0mAmp_width{25.0};
  constexpr int v0mAmp_Nbins{2000};
  double v0mAmp_bins[v0mAmp_Nbins + 1] = {0};
  for (int i = 0; i <= v0mAmp_Nbins; ++i) {
    v0mAmp_bins[i] = 0.0 + i * v0mAmp_width;
  }

  hNchFull = new TH1D("hNchFull", "; #it{N}_{ch} (|#eta|<0.8); Entries",
                      nchFull_Nbins, nchFull_bins);

  hNchV0M = new TH1D("hNchV0M", ";#it{N}_{ch} in the V0M; Entries",
                     v0mAmp_Nbins, v0mAmp_bins);

  hNchEtaGap = new TH1D("hNchEtaGap", ";#it{N};Entries", nchEtaGap_Nbins,
                        nchEtaGap_bins);

  hTPCFull = new TH1D("hTPCFull", ";#it{N}_{ch} (|#eta|<0.8);Entries",
                      nchFull_Nbins, nchFull_bins);

  hTPCEtaGap =
      new TH1D("hTPCEtaGap", ";#it{N}_{ch} (0.5#leq|#eta|#leq0.8);Entries",
               nchEtaGap_Nbins, nchEtaGap_bins);

  hSPDFull = new TH1D("hSPDFull", ";#it{N}_{ch} (|#eta|<0.8);Entries",
                      nchFull_Nbins, nchFull_bins);

  hSPDEtaGap =
      new TH1D("hSPDEtaGap", ";#it{N}_{ch} (0.5#leq|#eta|#leq0.8);Entries",
               nchEtaGap_Nbins, nchEtaGap_bins);

  hSPDEtaAdj =
      new TH1D("hSPDEtaAdj", ";#it{N}_{ch} (0.3<|#eta|#leq0.6);Entries",
               nchEtaGap_Nbins, nchEtaGap_bins);

  hSPDEtaGapW =
      new TH1D("hSPDEtaGapW", ";#it{N}_{ch} (0.7#leq|#eta|#leq1);Entries",
               nchEtaGap_Nbins, nchEtaGap_bins);

  hEtFull = new TH1D("hEtFull", ";#it{E}_{T} (|#eta|<0.8);Entries",
                     nchFull_Nbins, nchFull_bins);

  hEtEtaGap =
      new TH1D("hEtEtaGap", ";#it{N}_{ch} (0.5#leq|#eta|#leq0.8);Entries",
               nchEtaGap_Nbins, nchEtaGap_bins);

  hPtvsV0M = new TH2D(
      "hPtvsV0M",
      "; #it{N}_{ch} in the V0M; #it{p}_{T} (|#eta|#leq0.8, GeV/#it{c})",
      v0mAmp_Nbins, v0mAmp_bins, pt_Nbins, pt_bins);

  hPtvsNchFull = new TH2D(
      "hPtvsNchFull",
      ";#it{N}_{ch} (|#eta|<0.8);#it{p}_{T} (|#eta|#leq0.8, GeV/#it{c})",
      nchFull_Nbins, nchFull_bins, pt_Nbins, pt_bins);

  hPtvsNchEtaGap = new TH2D("hPtvsNchEtaGap",
                            ";#it{N}_{ch} (0.5#leq|#eta|#leq0.8);#it{p}_{T} "
                            "(|#eta|#leq0.3, GeV/#it{c})",
                            nchEtaGap_Nbins, nchEtaGap_bins, pt_Nbins, pt_bins);

  hPtvsTPCFull = new TH2D(
      "hPtvsTPCFull",
      ";#it{N}_{ch} (|#eta|<0.8);#it{p}_{T} (|#eta|#leq0.8, GeV/#it{c})",
      nchFull_Nbins, nchFull_bins, pt_Nbins, pt_bins);

  hPtvsTPCEtaGap = new TH2D("hPtvsTPCEtaGap",
                            ";#it{N}_{ch} (0.5#leq|#eta|#leq0.8); #it{p}_{T} "
                            "(|#eta|#leq0.3, GeV/#it{c})",
                            nchEtaGap_Nbins, nchEtaGap_bins, pt_Nbins, pt_bins);

  hPtvsSPDFull = new TH2D(
      "hPtvsSPDFull",
      ";#it{N}_{ch} (|#eta|<0.8);#it{p}_{T} (|#eta|#leq0.8, GeV/#it{c})",
      nchFull_Nbins, nchFull_bins, pt_Nbins, pt_bins);

  hPtvsSPDEtaGap = new TH2D("hPtvsSPDEtaGap",
                            "; #it{N}_{ch} (0.5#leq|#eta|#leq0.8); #it{p}_{T} "
                            "(|#eta#|leq0.3, GeV/#it{c})",
                            nchEtaGap_Nbins, nchEtaGap_bins, pt_Nbins, pt_bins);

  hPtvsSPDEtaAdj = new TH2D("hPtvsSPDEtaAdj",
                            "; #it{N}_{ch} (0.3<|#eta|#leq0.6); #it{p}_{T} "
                            "(|#eta#|leq0.3, GeV/#it{c})",
                            nchEtaGap_Nbins, nchEtaGap_bins, pt_Nbins, pt_bins);

  hPtvsSPDEtaGapW =
      new TH2D("hPtvsSPDEtaGapW",
               "; #it{N}_{ch} (0.7#leq|#eta|#leq1); #it{p}_{T} "
               "(|#eta#|leq0.4, GeV/#it{c})",
               nchEtaGap_Nbins, nchEtaGap_bins, pt_Nbins, pt_bins);

  hPtvsEtFull = new TH2D(
      "hPtvsEtFull",
      ";#it{E}_{T} (|#eta|<0.8);#it{p}_{T} (|#eta|#leq0.8, GeV/#it{c})",
      nchFull_Nbins, nchFull_bins, pt_Nbins, pt_bins);

  hPtvsEtEtaGap = new TH2D("hPtvsEtEtaGap",
                           ";#it{E}_{T} (0.5#leq|#eta|#leq0.8);#it{p}_{T} "
                           "(|#eta|#leq0.3, GeV/#it{c})",
                           nchEtaGap_Nbins, nchEtaGap_bins, pt_Nbins, pt_bins);

  fListOfObjects->Add(hNchV0M);
  fListOfObjects->Add(hPtvsV0M);
  fListOfObjects->Add(hNchFull);
  fListOfObjects->Add(hPtvsNchFull);
  fListOfObjects->Add(hNchEtaGap);
  fListOfObjects->Add(hPtvsNchEtaGap);
  fListOfObjects->Add(hTPCFull);
  fListOfObjects->Add(hPtvsTPCFull);
  fListOfObjects->Add(hTPCEtaGap);
  fListOfObjects->Add(hPtvsTPCEtaGap);
  fListOfObjects->Add(hSPDFull);
  fListOfObjects->Add(hPtvsSPDFull);
  fListOfObjects->Add(hSPDEtaGap);
  fListOfObjects->Add(hPtvsSPDEtaGap);
  fListOfObjects->Add(hSPDEtaAdj);
  fListOfObjects->Add(hPtvsSPDEtaAdj);
  fListOfObjects->Add(hSPDEtaGapW);
  fListOfObjects->Add(hPtvsSPDEtaGapW);
  fListOfObjects->Add(hEtFull);
  fListOfObjects->Add(hPtvsEtFull);
  fListOfObjects->Add(hEtEtaGap);
  fListOfObjects->Add(hPtvsEtEtaGap);

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

  bool isgoodvtx{IsGoodVertex()};
  if (!isgoodvtx) {
    return;
  }

  GetMultipliciy();

  PostData(1, fListOfObjects);

  return;
}
//______________________________________________________________________________
bool AliAnalysisTaskMCGenSpeedOfSound::IsGoodVertex() const {
  bool goodvetex{false};
  const AliVVertex *vtx = (AliVVertex *)fMcEvent->GetPrimaryVertex();
  if (!vtx) {
    return goodvetex;
  }

  double zvtx{vtx->GetZ()};
  if (zvtx <= 10.0) {
    goodvetex = true;
  }
  return goodvetex;
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
  double fSPDFull{0.0};
  double fSPDEtaAdj{0.0};
  double fSPDEtaGap{0.0};
  double fSPDEtaGapW{0.0};
  double fTPCFull{0.0};
  double fTPCEtaGap{0.0};
  double fEtFull{0.0};
  double fEtEtaGap{0.0};
  double fV0{0.0};
  double fNchFull{0.0};
  double fNchEtaGap{0.0};

  // ### particle loop
  for (Int_t ipart = 0; ipart < fMcEvent->GetNumberOfTracks(); ++ipart) {
    TParticle *mcPart = 0x0;
    mcPart = (TParticle *)fMcEvent->Particle(ipart);
    if (!mcPart) continue;
    // selection of primary charged particles
    if (!(mcPart->GetPDG())) continue;
    double qPart{mcPart->GetPDG()->Charge() / 3.};
    if (TMath::Abs(qPart) < 0.001) continue;
    bool isPhysPrim{fMcEvent->IsPhysicalPrimary(ipart)};
    if (!isPhysPrim) continue;

    double etaPart{mcPart->Eta()};
    double ptPart{mcPart->Pt()};
    double mass{mcPart->GetMass()};

    if (TMath::Abs(etaPart) <= 0.8) {
      if (ptPart >= 0.0) {
        fNchFull++;
      }
      if (ptPart >= 0.03) {
        fSPDFull++;
      }
      if (ptPart >= 0.15) {
        fTPCFull++;
        fEtFull += TMath::Sqrt(ptPart * ptPart + mass * mass);
      }
    }
    if (TMath::Abs(etaPart) > 0.3 && TMath::Abs(etaPart) <= 0.6) {
      if (ptPart >= 0.03) {
        fSPDEtaAdj++;
      }
    }
    if (TMath::Abs(etaPart) >= 0.7 && TMath::Abs(etaPart) <= 1.0) {
      if (ptPart >= 0.03) {
        fSPDEtaGapW++;
      }
    }
    if (TMath::Abs(etaPart) >= 0.5 && TMath::Abs(etaPart) <= 0.8) {
      if (ptPart >= 0.0) {
        fNchEtaGap++;
      }
      if (ptPart >= 0.03) {
        fSPDEtaGap++;
      }
      if (ptPart >= 0.15) {
        fTPCEtaGap++;
        fEtEtaGap += TMath::Sqrt(ptPart * ptPart + mass * mass);
      }
    }
    if ((2.8 < etaPart && etaPart < 5.1) ||
        (-3.7 < etaPart && etaPart < -1.7)) {
      fV0++;
    }
  }  // particle loop

  // ### particle loop 2
  for (Int_t ipart = 0; ipart < fMcEvent->GetNumberOfTracks(); ++ipart) {
    TParticle *mcPart = 0x0;
    mcPart = (TParticle *)fMcEvent->Particle(ipart);
    if (!mcPart) continue;
    // selection of primary charged particles
    if (!(mcPart->GetPDG())) continue;
    double qPart{mcPart->GetPDG()->Charge() / 3.};
    if (TMath::Abs(qPart) < 0.001) continue;
    bool isPhysPrim{fMcEvent->IsPhysicalPrimary(ipart)};
    if (!isPhysPrim) continue;

    double etaPart{mcPart->Eta()};
    double ptPart{mcPart->Pt()};

    if (TMath::Abs(etaPart) <= 0.8) {
      hPtvsV0M->Fill(fV0, ptPart);
      hPtvsNchFull->Fill(fNchFull, ptPart);
      hPtvsTPCFull->Fill(fTPCFull, ptPart);
      hPtvsSPDFull->Fill(fSPDFull, ptPart);
      hPtvsEtFull->Fill(fEtFull, ptPart);
    }
    if (TMath::Abs(etaPart) <= 0.3) {
      hPtvsNchEtaGap->Fill(fNchEtaGap, ptPart);
      hPtvsTPCEtaGap->Fill(fTPCEtaGap, ptPart);
      hPtvsSPDEtaAdj->Fill(fSPDEtaAdj, ptPart);
      hPtvsSPDEtaGap->Fill(fSPDEtaGap, ptPart);
      hPtvsSPDEtaGapW->Fill(fSPDEtaGapW, ptPart);
      hPtvsEtEtaGap->Fill(fEtEtaGap, ptPart);
    }
  }  // particle loop

  hNchV0M->Fill(fV0);
  hNchFull->Fill(fNchFull);
  hNchEtaGap->Fill(fNchEtaGap);
  hTPCFull->Fill(fTPCFull);
  hTPCEtaGap->Fill(fTPCEtaGap);
  hSPDFull->Fill(fSPDFull);
  hSPDEtaAdj->Fill(fSPDEtaAdj);
  hSPDEtaGap->Fill(fSPDEtaGap);
  hSPDEtaGapW->Fill(fSPDEtaGapW);
  hEtFull->Fill(fEtFull);
  hEtEtaGap->Fill(fEtEtaGap);
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
