/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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
 * Author:    Antonio Ortiz (antonio.ortiz@nucleares.unam.mx)             *
 *            Last modification: 28/04/2023                               *
 **************************************************************************/

//_____ ROOT headers
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3D.h"
#include "TObjArray.h"
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

//_____ AnalysisTask headers
#include "AliAnalysisTaskGenMcFlattenicity.h"
#include "AliAnalysisTaskSE.h"
//_____ STL includes
#include <iostream>
using namespace std;

float min_flat = -0.01;
float max_flat = 1.01;
int nbins_flat = 1020;

// V0A
const int nRingsV0A = 4;
const int nSectorsV0A = 8;
float maxEtaV0A[nRingsV0A] = {-3.2, -2.7, -2.2, -1.7};
float minEtaV0A[nRingsV0A] = {-3.7, -3.2, -2.7, -2.2};
float maxPhiV0A[nSectorsV0A] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
float minPhiV0A[nSectorsV0A] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};
// V0C
const int nRingsV0C = 4;
const int nSectorsV0C = 8;
float maxEtaV0C[nRingsV0C] = {5.1, 4.5, 3.9, 3.4};
float minEtaV0C[nRingsV0C] = {4.5, 3.9, 3.4, 2.8};
float maxPhiV0C[nSectorsV0C] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
float minPhiV0C[nSectorsV0C] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0};

ClassImp(AliAnalysisTaskGenMcFlattenicity)

    AliAnalysisTaskGenMcFlattenicity::AliAnalysisTaskGenMcFlattenicity()
    : AliAnalysisTaskSE(), fMC(0x0), fMcHandler(0x0), fMCStack(0x0),
      fGenerator(0), fEtaCut(0.8), fIsPP(kTRUE), fPtMin(0.15), fIsINEL0(0),
      fNchv0(-1), fNchv0a(-1), fNchv0c(-1), fFlattV0(-1), hFlatt(0x0),
      hflatVsNchV0(0x0), hMultVsFlat(0x0), fOutputList(0) {
  for (Int_t pid = 0; pid < 4; ++pid) {
    hMultVsFlatVsPt[pid] = 0;
  }
  for (Int_t pid = 0; pid < 4; ++pid) {
    hFlatVsPt[pid] = 0;
  }
}
//_____________________________________________________________________________
AliAnalysisTaskGenMcFlattenicity::AliAnalysisTaskGenMcFlattenicity(
    const char *name)
    : AliAnalysisTaskSE(name), fMC(0x0), fMcHandler(0x0), fMCStack(0x0),
      fGenerator(0), fEtaCut(0.8), fIsPP(kTRUE), fPtMin(0.15), fIsINEL0(0),
      fNchv0(-1), fNchv0a(-1), fNchv0c(-1), fFlattV0(-1), hFlatt(0x0),
      hflatVsNchV0(0x0), hMultVsFlat(0x0), fOutputList(0) {
  for (Int_t pid = 0; pid < 4; ++pid) {
    hMultVsFlatVsPt[pid] = 0;
  }
  for (Int_t pid = 0; pid < 4; ++pid) {
    hFlatVsPt[pid] = 0;
  }
  // constructor
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskGenMcFlattenicity::~AliAnalysisTaskGenMcFlattenicity() {
  // destructor
  if (fOutputList) {
    delete fOutputList;
    fOutputList = 0x0;
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskGenMcFlattenicity::UserCreateOutputObjects() {
  // ### Analysis output
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);
  // create output objects

  // flattenicit binning
  const int nFlatA = 1020;
  double FlatA[nFlatA + 1];
  for (int i = 0; i < 1020; ++i) {
    FlatA[i] = 0;
    FlatA[i] = -0.01 + i * 0.001;
  }
  FlatA[1020] = 1.01;
  // multiplicity binning
  const int NnchBin = 400;
  double nchBin[NnchBin + 1];
  for (int i = 0; i < 400; ++i) {
    nchBin[i] = 0;
    nchBin[i] = i - 0.5;
  }
  nchBin[400] = 399.5;

  const int nPtBins = 43;
  double PtBins[nPtBins + 1] = {
      0.,  0.15, 0.25, 0.3,  0.35, 0.4,  0.45, 0.5, 0.55, 0.6,  0.65,
      0.7, 0.75, 0.8,  0.85, 0.9,  0.95, 1.0,  1.1, 1.2,  1.3,  1.4,
      1.5, 1.6,  1.7,  1.8,  1.9,  2.0,  2.2,  2.4, 2.6,  2.8,  3.0,
      3.2, 3.4,  3.6,  3.8,  4.0,  4.5,  5.0,  6.0, 8.0,  10.0, 20.0};
  const char *Pid[4] = {"Ch", "Pion", "Kaon", "Proton"};

  hFlatt = 0;
  hFlatt = new TH1D("hFlatt", ";flattenicity; n. ev.", nbins_flat, min_flat,
                    max_flat);
  fOutputList->Add(hFlatt);
  for (int pid = 0; pid < 4; ++pid) {
    hFlatVsPt[pid] = 0;
    hFlatVsPt[pid] = new TH2D(Form("hFlatVsPt_%s", Pid[pid]),
                              "; Flatenicity; #it{p}_{T} (GeV/#it{c})",
                              nbins_flat, min_flat, max_flat, nPtBins, PtBins);
    fOutputList->Add(hFlatVsPt[pid]);

    hMultVsFlatVsPt[pid] = 0;
    hMultVsFlatVsPt[pid] =
        new TH3D(Form("hMultVsFlatVsPt_%s", Pid[pid]),
                 "; true multiplicity; true flattenicity", NnchBin, nchBin,
                 nFlatA, FlatA, nPtBins, PtBins);
    fOutputList->Add(hMultVsFlatVsPt[pid]);
  }

  hMultVsFlat = 0;
  hMultVsFlat =
      new TH2D("hMultVsFlat", "; true multiplicity; true flattenicity", 400,
               -0.5, 399.5, nbins_flat, min_flat, max_flat);
  fOutputList->Add(hMultVsFlat);

  hflatVsNchV0 = 0;
  hflatVsNchV0 =
      new TH2D("hflatVsNchV0", "; #it{N}_{ch} (V0);1-flattenicity (V0)", 500,
               -0.5, 499.5, nbins_flat, min_flat, max_flat);
  fOutputList->Add(hflatVsNchV0);

  PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskGenMcFlattenicity::UserExec(Option_t *) {
  // ### Initialize

  fMcHandler = dynamic_cast<AliInputEventHandler *>(
      AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());

  // ### MC handler

  if (fMcHandler)
    fMC = fMcHandler->MCEvent();

  else {
    if (fDebug > 1)
      printf("AliAnalysisTaskGenUeNchTS::Handler() fMcHandler = NULL\n");
    return;
  }

  // ### MC event

  if (!fMC) {
    Printf("%s:%d MCEvent not found in Input Manager", (char *)__FILE__,
           __LINE__);
    this->Dump();
    return;
  }
  fMCStack = ((AliMCEvent *)fMC)->Stack();

  if (!fMCStack) {
    Printf("ERROR: Could not retrieve MC stack \n");
    cout << "Name of the file with pb :" << endl;
    return;
  }

  // ### MC event selection

  bool isEventMCSelected = IsMCEventSelected(fMC);
  if (!isEventMCSelected)
    return;

  AliHeader *headerMC = fMC->Header();

  bool isGoodVtxPosMC = kFALSE;

  AliGenEventHeader *genHeader = headerMC->GenEventHeader();
  TArrayF vtxMC(3); // primary vertex  MC
  vtxMC[0] = 9999;
  vtxMC[1] = 9999;
  vtxMC[2] = 9999; // initialize with dummy
  if (genHeader) {
    genHeader->PrimaryVertex(vtxMC);
  }
  if (TMath::Abs(vtxMC[2]) <= 10)
    isGoodVtxPosMC = true;
  if (!isGoodVtxPosMC)
    return;
  // Before trigger selection
  GetFlattenicity(); // leading particle at gen level
  hFlatt->Fill(fFlattV0);
  if (fNchv0a > 0 && fNchv0c > 0 && fIsINEL0) {
    hMultVsFlat->Fill(fNchv0, fFlattV0);
    hflatVsNchV0->Fill(fNchv0, 1. - fFlattV0);
    MakeAnalysis();
  }

  PostData(1, fOutputList);
  return;
}

//_____________________________________________________________________________
bool AliAnalysisTaskGenMcFlattenicity::IsMCEventSelected(TObject *obj) {

  bool isSelected = kTRUE;

  AliMCEvent *event = 0x0;
  event = dynamic_cast<AliMCEvent *>(obj);
  if (!event)
    isSelected = kFALSE;

  return isSelected;
}
//-------------------------------------------------
void AliAnalysisTaskGenMcFlattenicity::GetFlattenicity() {

  double rho = -1;
  int nCells;
  int nCells1 = nRingsV0A * nSectorsV0A;
  int nCells2 = nRingsV0C * nSectorsV0C;
  nCells = nCells1 + nCells2;

  double RhoLattice[nCells];

  for (int iCh = 0; iCh < nCells; ++iCh) {
    RhoLattice[iCh] = 0;
  }

  int nchtotal_tmp = 0;
  fNchv0a = 0;
  fNchv0c = 0;
  float eta_a, phi_a, pt_a, q_a;
  int mult_inel0;
  mult_inel0 = 0;

  for (int i = 0; i < fMC->GetNumberOfTracks(); i++) {

    AliMCParticle *particle = (AliMCParticle *)fMC->GetTrack(i);
    if (!particle)
      continue;
    if (!fMC->IsPhysicalPrimary(i))
      continue;

    eta_a = particle->Eta();
    pt_a = particle->Pt();
    phi_a = particle->Phi();

    q_a = particle->Charge();

    if (TMath::Abs(q_a) < 0.01)
      continue;
    if (pt_a <= 0)
      continue;

    // check mult eta less1
    if (TMath::Abs(eta_a) < 1.) {
      mult_inel0++;
    }
    // first part of the detector
    int i_ch = 0;
    for (int ir = 0; ir < nRingsV0A; ++ir) {
      for (int is = 0; is < nSectorsV0A; ++is) {
        if (eta_a >= minEtaV0A[ir] && eta_a < maxEtaV0A[ir] &&
            phi_a >= minPhiV0A[is] * 2.0 * M_PI / (1.0 * nSectorsV0A) &&
            phi_a < maxPhiV0A[is] * 2.0 * M_PI / (1.0 * nSectorsV0A)) {
          nchtotal_tmp++;
          fNchv0a++;
          RhoLattice[i_ch]++;
        }
        i_ch++;
      }
    }

    // second part of the detector
    for (int ir = 0; ir < nRingsV0C; ++ir) {
      for (int is = 0; is < nSectorsV0C; ++is) {
        if (eta_a >= minEtaV0C[ir] && eta_a < maxEtaV0C[ir] &&
            phi_a >= minPhiV0C[is] * 2.0 * M_PI / (1.0 * nSectorsV0C) &&
            phi_a < maxPhiV0C[is] * 2.0 * M_PI / (1.0 * nSectorsV0C)) {
          nchtotal_tmp++;
          fNchv0c++;
          RhoLattice[i_ch]++;
        }
        i_ch++;
      }
    }
  }

  // normalization to Delta_eta
  int i_ch = 0;
  for (int ir = 0; ir < nRingsV0A; ++ir) {
    double deta = TMath::Abs(maxEtaV0A[ir] - minEtaV0A[ir]);
    RhoLattice[i_ch] /= deta;
    i_ch++;
  }
  for (int ir = 0; ir < nRingsV0C; ++ir) {
    double deta = TMath::Abs(maxEtaV0C[ir] - minEtaV0C[ir]);
    RhoLattice[i_ch] /= deta;
    i_ch++;
  }

  double mRho = 0;
  for (int iCh = 0; iCh < nCells; ++iCh) {
    mRho += RhoLattice[iCh];
  }
  mRho /= (1.0 * nCells);
  // get sigma
  double sRho_tmp = 0;

  for (int iCh = 0; iCh < nCells; ++iCh) {
    sRho_tmp += TMath::Power(RhoLattice[iCh] - mRho, 2);
  }
  sRho_tmp /= (1.0 * nCells * nCells);
  double sRho = TMath::Sqrt(sRho_tmp);

  if (mRho > 0) {
    rho = sRho / mRho;
  } else {
    rho = 9999;
  }
  fNchv0 = -1;
  fNchv0 = nchtotal_tmp;
  fFlattV0 = rho;
  fIsINEL0 = false;
  if (mult_inel0 > 0) {
    fIsINEL0 = true;
  }
}
//----------------------
void AliAnalysisTaskGenMcFlattenicity::MakeAnalysis() {

  float y_a, eta_a, pt_a, q_a;
  int pdg;

  for (int i = 0; i < fMC->GetNumberOfTracks(); i++) {

    AliMCParticle *particle = (AliMCParticle *)fMC->GetTrack(i);
    if (!particle)
      continue;
    if (!fMC->IsPhysicalPrimary(i))
      continue;

    eta_a = particle->Eta();
    y_a = particle->Y();
    pt_a = particle->Pt();

    q_a = particle->Charge();
    pdg = GetPidCode(particle->PdgCode());

    if (TMath::Abs(q_a) < 0.01)
      continue;

    if (pdg >= 0 && pdg < 3) {
      if (TMath::Abs(y_a) <= fEtaCut) {
        hFlatVsPt[pdg + 1]->Fill(fFlattV0, pt_a);
        hMultVsFlatVsPt[pdg + 1]->Fill(fNchv0, fFlattV0, pt_a);
      }
    }

    // check mult eta less1
    if (TMath::Abs(eta_a) > fEtaCut) {
      continue;
    }

    hFlatVsPt[0]->Fill(fFlattV0, pt_a);
    hMultVsFlatVsPt[0]->Fill(fNchv0, fFlattV0, pt_a);
  }
}
//_______________________________________________________
int AliAnalysisTaskGenMcFlattenicity::GetPidCode(int pdgCode) {

  int pidCode = 3;

  switch (TMath::Abs(pdgCode)) {
  case 211:
    pidCode = 0; // pion
    break;
  case 321:
    pidCode = 1; // kaon
    break;
  case 2212:
    pidCode = 2; // proton
    break;
    // default:
    //	3;
  };

  return pidCode;
}
//______________________________________________________________________________
void AliAnalysisTaskGenMcFlattenicity::Terminate(Option_t *) {
  fOutputList = dynamic_cast<TList *>(GetOutputData(1));
  if (!fOutputList) {
    Printf("ERROR: Output list not available");
    return;
  }

  return;
}
