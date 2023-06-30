/****************************************************************************
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
 *                                                                        *
 * Author: Antonio Ortiz (antonio.ortiz@nucleares.unam.mx)                *
 * Anatask to compute flatenicity: PRD 107 (2023) 7, 076012               *
 **************************************************************************/

class TTree;

class AliESDtrackCuts;

#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisTaskESDfilter.h"
#include "AliAnalysisTaskSE.h"
#include "AliCentrality.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDUtils.h"
#include "AliESDVZERO.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliEventCuts.h"
#include "AliGenCocktailEventHeader.h"
#include "AliGenEventHeader.h"
#include "AliInputEventHandler.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliMCParticle.h"

#include "AliMultEstimator.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliMultVariable.h"
#include "AliMultiplicity.h"
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliStack.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TLegend.h"
#include "TList.h"
#include "TMath.h"
#include "TParticle.h"
#include "TProfile.h"
#include "TVector3.h"
#include <AliAnalysisFilter.h>
#include <AliESDVertex.h>
#include <AliHeader.h>
#include <AliMultiplicity.h>
#include <Riostream.h>
#include <TDirectory.h>
#include <TMath.h>
#include <TTree.h>
using std::cout;
using std::endl;

#include "AliAnalysisTaskFlatenicityMCpred.h"

const int nCent = 9;
double v0mClass[nCent + 1] = {0.0,  1.0,  5.0,  10.0, 20.0,
                              30.0, 40.0, 50.0, 70.0, 100.0};

using namespace std; // std namespace: so you can do things like 'cout' etc

ClassImp(AliAnalysisTaskFlatenicityMCpred) // classimp: necessary for root

    AliAnalysisTaskFlatenicityMCpred::AliAnalysisTaskFlatenicityMCpred()
    : AliAnalysisTaskSE(), fESD(0), fEventCuts(0x0), fMCStack(0), fMC(0),
      fVtxz(-1), fminrho(0.), fmaxrho(1.), fNchBin(600), fV0Mindex(-1),
      fv0multalice(-1), flatenicity_m(-1), fmultV0A_m(-1), fmultV0C_m(-1),
      flatenicity_t(-1), fmultV0A_t(-1), fmultV0C_t(-1), fnchAll(-1),
      fOutputList(0), fEtaCut(0.8), fv0mpercentile(0), fMultSelection(0x0),
      hCounter(0), hV0MBadruns(0), hMultVsFlat1(0x0), hflatVsNchV01(0x0) {
  for (int pid = 0; pid < 4; ++pid) { // add
    hMultVsFlatVsPt1[pid] = 0;
  }
  for (int pid = 0; pid < 4; ++pid) { // add
    hFlatVsPt1[pid] = 0;
  }
  for (int i_c = 0; i_c < nCent; ++i_c) {
    hMultV0MPerc[i_c] = 0;
  }
  for (int i_c = 0; i_c < nCent; ++i_c) {
    hFlatV0MPerc[i_c] = 0;
  }
  for (int pid = 0; pid < 4; ++pid) { // add
    for (int i_c = 0; i_c < nCent; ++i_c) {
      hFlatVsPt1V0MPerc[pid][i_c] = 0;
    }
  }
}

//_____________________________________________________________________________
AliAnalysisTaskFlatenicityMCpred::AliAnalysisTaskFlatenicityMCpred(
    const char *name)
    : AliAnalysisTaskSE(name), fESD(0), fEventCuts(0x0), fMCStack(0), fMC(0),
      fVtxz(-1), fminrho(0.), fmaxrho(1.), fNchBin(600), fV0Mindex(-1),
      fv0multalice(-1), flatenicity_m(-1), fmultV0A_m(-1), fmultV0C_m(-1),
      flatenicity_t(-1), fmultV0A_t(-1), fmultV0C_t(-1), fnchAll(-1),
      fOutputList(0), fEtaCut(0.8), fv0mpercentile(0), fMultSelection(0x0),
      hCounter(0), hV0MBadruns(0), hMultVsFlat1(0x0), hflatVsNchV01(0x0)

{

  for (int pid = 0; pid < 4; ++pid) {
    hMultVsFlatVsPt1[pid] = 0;
  }
  for (int pid = 0; pid < 4; ++pid) {
    hFlatVsPt1[pid] = 0;
  }
  for (int i_c = 0; i_c < nCent; ++i_c) {
    hMultV0MPerc[i_c] = 0;
  }
  for (int i_c = 0; i_c < nCent; ++i_c) {
    hFlatV0MPerc[i_c] = 0;
  }
  for (int pid = 0; pid < 4; ++pid) { // add
    for (int i_c = 0; i_c < nCent; ++i_c) {
      hFlatVsPt1V0MPerc[pid][i_c] = 0;
    }
  }

  DefineInput(0, TChain::Class()); // define the input of the analysis: in this
                                   // case you take a 'chain' of events
  DefineOutput(1, TList::Class()); // define the ouptut of the analysis: in this
                                   // case it's a list of histograms
}

//_____________________________________________________________________________
AliAnalysisTaskFlatenicityMCpred::~AliAnalysisTaskFlatenicityMCpred() {
  // destructor
  if (fOutputList) {
    delete fOutputList; // at the end of your task, it is deleted from memory by
                        // calling this function
    fOutputList = 0x0;
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskFlatenicityMCpred::UserCreateOutputObjects() {

  // flattenicit binning
  // const int nFlatA = 1020;
  double binsizerho = 0.001;
  int nFlatA = int((fmaxrho - fminrho) / binsizerho);
  double FlatA[nFlatA + 1];
  for (int i = 0; i < nFlatA; ++i) {
    FlatA[i] = 0;
    FlatA[i] = fminrho + i * binsizerho;
  }
  FlatA[nFlatA] = fmaxrho;
  // multiplicity binning
  const int NnchBin = fNchBin;
  double nchBin[NnchBin + 1];
  for (int i = 0; i < NnchBin; ++i) {
    nchBin[i] = 0;
    nchBin[i] = i - 0.5;
  }
  nchBin[NnchBin] = 1. * NnchBin - 0.5;

  const int nPtBins = 43;
  double PtBins[nPtBins + 1] = {
      0.,  0.15, 0.25, 0.3,  0.35, 0.4,  0.45, 0.5, 0.55, 0.6,  0.65,
      0.7, 0.75, 0.8,  0.85, 0.9,  0.95, 1.0,  1.1, 1.2,  1.3,  1.4,
      1.5, 1.6,  1.7,  1.8,  1.9,  2.0,  2.2,  2.4, 2.6,  2.8,  3.0,
      3.2, 3.4,  3.6,  3.8,  4.0,  4.5,  5.0,  6.0, 8.0,  10.0, 20.0};
  const char *Pid[4] = {"Ch", "Pion", "Kaon", "Proton"};

  OpenFile(1);
  fOutputList =
      new TList(); // this is a list which will contain all of your histograms
  fOutputList->SetOwner(kTRUE); // memory stuff: the list is owner of all

  hMultVsFlat1 = 0;
  hMultVsFlat1 =
      new TH2D("hMultVsFlat1", "; V0M amplitude; measured flattenicity (V0)",
               NnchBin, nchBin, nFlatA, FlatA);
  fOutputList->Add(hMultVsFlat1);

  hflatVsNchV01 = 0;
  hflatVsNchV01 =
      new TH2D("hflatVsNchV01", "; V0M amplitude;1-flattenicity (V0)", NnchBin,
               nchBin, nFlatA, FlatA);
  fOutputList->Add(hflatVsNchV01);

  for (int pid = 0; pid < 4; ++pid) {
    hFlatVsPt1[pid] = 0;
    hFlatVsPt1[pid] =
        new TH2D(Form("hFlatVsPt1_%s", Pid[pid]),
                 "; measured flattenicity (V0); #it{p}_{T} (GeV/#it{c})",
                 nFlatA, FlatA, nPtBins, PtBins);
    fOutputList->Add(hFlatVsPt1[pid]);

    hMultVsFlatVsPt1[pid] = 0;
    hMultVsFlatVsPt1[pid] =
        new TH3D(Form("hMultVsFlatVsPt1_%s", Pid[pid]),
                 "; V0M amplitude; measured flattenicity (V0M)", NnchBin,
                 nchBin, nFlatA, FlatA, nPtBins, PtBins);
    fOutputList->Add(hMultVsFlatVsPt1[pid]);
  }

  for (int i_c = 0; i_c < nCent; ++i_c) {
    hMultV0MPerc[i_c] = 0;
    hMultV0MPerc[i_c] =
        new TH1D(Form("hMultV0MPerc_%d", i_c), "", NnchBin, nchBin);
    fOutputList->Add(hMultV0MPerc[i_c]);
  }

  for (int i_c = 0; i_c < nCent; ++i_c) {
    hFlatV0MPerc[i_c] = 0;
    hFlatV0MPerc[i_c] =
        new TH1D(Form("hFlatV0MPerc_%d", i_c), "", nFlatA, FlatA);
    fOutputList->Add(hFlatV0MPerc[i_c]);
  }

  for (int pid = 0; pid < 4; ++pid) { // add
    for (int i_c = 0; i_c < nCent; ++i_c) {
      hFlatVsPt1V0MPerc[pid][i_c] = 0;
      hFlatVsPt1V0MPerc[pid][i_c] = 0;
      hFlatVsPt1V0MPerc[pid][i_c] =
          new TH2D(Form("hFlatVsPt1_%s_V0MPerc_%d", Pid[pid], i_c),
                   "; measured flattenicity (V0); #it{p}_{T} (GeV/#it{c})",
                   nFlatA, FlatA, nPtBins, PtBins);
      fOutputList->Add(hFlatVsPt1V0MPerc[pid][i_c]);
    }
  }

  hCounter = new TH1D("hCounter", "counter", 15, -0.5, 14.5);
  fOutputList->Add(hCounter);

  hV0MBadruns = new TH1D("hV0MBadruns", "", NnchBin, nchBin);
  fOutputList->Add(hV0MBadruns);

  fEventCuts.AddQAplotsToList(fOutputList);
  PostData(1, fOutputList); // postdata will notify the analysis manager of
                            // changes / updates to the
}

//_____________________________________________________________________________
void AliAnalysisTaskFlatenicityMCpred::UserExec(Option_t *) {

  AliVEvent *event = InputEvent();
  if (!event) {
    Error("UserExec", "Could not retrieve event");
    return;
  }

  fESD = dynamic_cast<AliESDEvent *>(event);

  if (!fESD) {
    Printf("%s:%d ESDEvent not found in Input Manager", (char *)__FILE__,
           __LINE__);
    this->Dump();
    return;
  }

  hCounter->Fill(0.0);

  //      E S D
  fMC = dynamic_cast<AliMCEvent *>(MCEvent());
  if (!fMC) {
    Printf("%s:%d MCEvent not found in Input Manager", (char *)__FILE__,
           __LINE__);
    this->Dump();
    return;
  }
  fMCStack = fMC->Stack();

  AliHeader *headerMC;
  // bool isGoodVtxPosMC = kFALSE;

  headerMC = fMC->Header();
  AliGenEventHeader *genHeader = headerMC->GenEventHeader();
  TArrayF vtxMC(3); // primary vertex  MC
  vtxMC[0] = 9999;
  vtxMC[1] = 9999;
  vtxMC[2] = 9999; // initialize with dummy
  if (genHeader)
    genHeader->PrimaryVertex(vtxMC);

  // if (TMath::Abs(vtxMC[2]) <= 10)
  //   isGoodVtxPosMC = kTRUE;

  /*
  HERE PURE MC TRUE
  */
  flatenicity_t = -1;
  flatenicity_t = GetFlatenicityMC();
  ExtractMultiplicitiesMC(); // here I am getting: fmultV0A_t and fmultV0C_t
  vector<double> ptMC;
  vector<double> yMC;
  vector<double> etaMC;
  vector<int> idMC;
  fnchAll = -1;
  fnchAll = FillMCarray(ptMC, yMC, etaMC, idMC);

  fMultSelection = (AliMultSelection *)fESD->FindListObject("MultSelection");
  if (!fMultSelection)
    cout << "------- No AliMultSelection Object Found --------"
         << fMultSelection << endl;

  fv0mpercentile = -999;
  fv0mpercentile = fMultSelection->GetMultiplicityPercentile("V0M");
  fv0multalice = -1;
  fv0multalice = fMultSelection->GetEstimator("V0M")->GetValue();

  // Trigger selection
  UInt_t fSelectMask = fInputHandler->IsEventSelected();
  bool isINT7selected = fSelectMask & AliVEvent::kINT7;
  if (!isINT7selected)
    return;
  hCounter->Fill(1.0);

  // Good vertex
  bool hasRecVertex = kFALSE;
  hasRecVertex = HasRecVertex();
  if (!hasRecVertex)
    return;

  hCounter->Fill(3.0);

  // Good events
  if (!fEventCuts.AcceptEvent(event)) {
    PostData(1, fOutputList);
    return;
  }
  flatenicity_m = -1;
  flatenicity_m = GetFlatenicityV0();
  ExtractMultiplicities(); // here fmultV0C and fmultV0C are obtained

  hCounter->Fill(2.0);

  fV0Mindex = -1;
  for (int i_c = 0; i_c < nCent; ++i_c) {
    if (fv0mpercentile >= v0mClass[i_c] && fv0mpercentile < v0mClass[i_c + 1]) {
      fV0Mindex = i_c;
    } else {
      continue;
    }
  }

  if (fV0Mindex >= 0 && flatenicity_m >= 0.) {
    if ((fV0Mindex == nCent - 1) && (fv0multalice > 400)) {
      hV0MBadruns->Fill(fv0multalice);
    } else {
      hMultV0MPerc[fV0Mindex]->Fill(fv0multalice);
      hMultVsFlat1->Fill(fv0multalice, flatenicity_m);
      hflatVsNchV01->Fill(fv0multalice, 1. - flatenicity_m);
      FillHistos1(fnchAll, ptMC, yMC, etaMC, idMC);
    }
  }

  PostData(1, fOutputList); // stream the result of this event to the output
                            // manager which will write it to a file
}

//______________________________________________________________________________
void AliAnalysisTaskFlatenicityMCpred::Terminate(Option_t *) {}

//______________________________________________________________________________
double AliAnalysisTaskFlatenicityMCpred::GetFlatenicityV0() {

  AliVVZERO *lVV0 = 0x0;
  AliVEvent *lVevent = 0x0;
  lVevent = dynamic_cast<AliVEvent *>(InputEvent());
  if (!lVevent) {
    AliWarning("ERROR: ESD / AOD event not available \n");
    return -1;
  }
  // Get VZERO Information for multiplicity later
  lVV0 = lVevent->GetVZEROData();
  if (!lVV0) {
    AliError("AliVVZERO not available");
    return 9999;
  }

  // Flatenicity calculation
  const int nRings = 4;
  const int nSectors = 8;
  double minEtaV0C[nRings] = {-3.7, -3.2, -2.7, -2.2};
  double maxEtaV0C[nRings] = {-3.2, -2.7, -2.2, -1.7};
  double maxEtaV0A[nRings] = {5.1, 4.5, 3.9, 3.4};
  double minEtaV0A[nRings] = {4.5, 3.9, 3.4, 2.8};
  // Grid
  const int nCells = nRings * 2 * nSectors;
  double RhoLattice[nCells];
  for (int iCh = 0; iCh < nCells; iCh++) {
    RhoLattice[iCh] = 0.0;
  }

  // before calibration
  for (int iCh = 0; iCh < nCells; iCh++) {
    double mult = lVV0->GetMultiplicity(iCh);
    RhoLattice[iCh] = mult;
  }

  int nringA = 0;
  int nringC = 0;
  for (int iCh = 0; iCh < nCells; iCh++) {
    double detaV0 = -1;
    if (iCh < 32) { // V0C
      if (iCh < 8) {
        nringC = 0;
      } else if (iCh >= 8 && iCh < 16) {
        nringC = 1;
      } else if (iCh >= 16 && iCh < 24) {
        nringC = 2;
      } else {
        nringC = 3;
      }
      detaV0 = maxEtaV0C[nringC] - minEtaV0C[nringC];
    } else { // V0A
      if (iCh < 40) {
        nringA = 0;
      } else if (iCh >= 40 && iCh < 48) {
        nringA = 1;
      } else if (iCh >= 48 && iCh < 56) {
        nringA = 2;
      } else {
        nringA = 3;
      }
      detaV0 = maxEtaV0A[nringA] - minEtaV0A[nringA];
    }
    // consider the different eta coverage
    RhoLattice[iCh] /= detaV0;
  }
  double mRho = 0;
  double flatenicity = -1;
  for (int iCh = 0; iCh < nCells; iCh++) {
    mRho += RhoLattice[iCh];
  }
  // double multiplicityV0M = mRho;
  //  average activity per cell
  mRho /= (1.0 * nCells);
  // get sigma
  double sRho_tmp = 0;
  for (int iCh = 0; iCh < nCells; iCh++) {
    sRho_tmp += TMath::Power(1.0 * RhoLattice[iCh] - mRho, 2);
  }
  sRho_tmp /= (1.0 * nCells * nCells);
  double sRho = TMath::Sqrt(sRho_tmp);
  if (mRho > 0) {
    flatenicity = sRho / mRho;
  } else {
    flatenicity = 9999;
  }
  return flatenicity;
}

//______________________________________________________________________________
double AliAnalysisTaskFlatenicityMCpred::GetFlatenicityMC() {

  // Flatenicity calculation
  const int nRings = 8;
  double maxEta[nRings] = {-3.2, -2.7, -2.2, -1.7, 5.1, 4.5, 3.9, 3.4};
  double minEta[nRings] = {-3.7, -3.2, -2.7, -2.2, 4.5, 3.9, 3.4, 2.8};

  const int nSectors = 8;
  double PhiBins[nSectors + 1];
  double deltaPhi = (2.0 * TMath::Pi()) / (1.0 * nSectors);
  for (int i_phi = 0; i_phi < nSectors + 1; ++i_phi) {
    PhiBins[i_phi] = 0;
    if (i_phi < nSectors) {
      PhiBins[i_phi] = i_phi * deltaPhi;
    } else {
      PhiBins[i_phi] = 2.0 * TMath::Pi();
    }
  }

  // Grid
  const int nCells = nRings * nSectors;
  double RhoLattice[nCells];
  for (int iCh = 0; iCh < nCells; iCh++) {
    RhoLattice[iCh] = 0.0;
  }

  int nMult = 0;
  for (int i = 0; i < fMC->GetNumberOfTracks(); ++i) {

    AliMCParticle *particle = (AliMCParticle *)fMC->GetTrack(i);
    if (!particle)
      continue;
    if (!fMC->IsPhysicalPrimary(i))
      continue;
    if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i, fMC))
      continue;
    if (particle->Pt() <= 0.0)
      continue;
    if (TMath::Abs(particle->Charge()) < 0.1)
      continue;
    double phi = particle->Phi();
    double eta = particle->Eta();

    int i_segment = 0;
    for (int i_eta = 0; i_eta < nRings; ++i_eta) {

      for (int i_phi = 0; i_phi < nSectors; ++i_phi) {

        if (eta >= minEta[i_eta] && eta < maxEta[i_eta] &&
            phi >= PhiBins[i_phi] && phi < PhiBins[i_phi + 1]) {
          nMult++;
          RhoLattice[i_segment] += 1.0;
        }
        i_segment++;
      }
    }
  }

  int i_segment = 0;
  for (int i_eta = 0; i_eta < nRings; ++i_eta) {
    for (int i_phi = 0; i_phi < nSectors; ++i_phi) {
      double deltaEta = TMath::Abs(maxEta[i_eta] - minEta[i_eta]);
      RhoLattice[i_segment] /= deltaEta;
      i_segment++;
    }
  }

  double mRho = 0;
  double flatenicity = -1;
  for (int iCh = 0; iCh < nCells; iCh++) {
    mRho += RhoLattice[iCh];
  }
  // average activity per cell
  mRho /= (1.0 * nCells);
  // get sigma
  double sRho_tmp = 0;
  for (int iCh = 0; iCh < nCells; iCh++) {
    sRho_tmp += TMath::Power(1.0 * RhoLattice[iCh] - mRho, 2);
  }
  sRho_tmp /= (1.0 * nCells * nCells);
  double sRho = TMath::Sqrt(sRho_tmp);

  if (mRho > 0) {
    flatenicity = sRho / mRho;
  } else {
    sRho = 9999;
  }
  return flatenicity;
}

//______________________________________________________________________________
void AliAnalysisTaskFlatenicityMCpred::ExtractMultiplicities() {

  AliVVZERO *lVV0 = 0x0;
  AliVEvent *lVevent = 0x0;
  lVevent = dynamic_cast<AliVEvent *>(InputEvent());
  if (!lVevent) {
    AliWarning("ERROR: ESD / AOD event not available \n");
    return;
  }
  // Get VZERO Information for multiplicity later
  lVV0 = lVevent->GetVZEROData();
  if (!lVV0) {
    AliError("AliVVZERO not available");
    return;
  }

  const int nChannels = 64;
  double fmultV0C_tmp = 0;
  double fmultV0A_tmp = 0;
  for (int iCh = 0; iCh < nChannels; iCh++) {
    double mult = lVV0->GetMultiplicity(iCh);
    if (iCh < 32) { // V0C
      fmultV0C_tmp += mult;
    } else { // V0A
      fmultV0A_tmp += mult;
    }
  }

  fmultV0C_m = fmultV0C_tmp;
  fmultV0A_m = fmultV0A_tmp;
}

//______________________________________________________________________________
void AliAnalysisTaskFlatenicityMCpred::ExtractMultiplicitiesMC() {

  fmultV0A_t = 0;
  fmultV0C_t = 0;

  for (int i = 0; i < fMC->GetNumberOfTracks(); ++i) {

    AliMCParticle *particle = (AliMCParticle *)fMC->GetTrack(i);
    if (!particle)
      continue;
    if (!fMC->IsPhysicalPrimary(i))
      continue;
    if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i, fMC))
      continue;
    if (particle->Pt() <= 0.0)
      continue;
    if (TMath::Abs(particle->Charge()) < 0.1)
      continue;
    double eta_a = particle->Eta();
    if (eta_a >= 2.8 && eta_a < 5.1) { // v0a acceptance (excluding first ring)
      fmultV0A_t++;
    }
    if (eta_a >= -3.7 && eta_a < -1.7) { // v0c
      fmultV0C_t++;
    }
  }
}

//______________________________________________________________________
bool AliAnalysisTaskFlatenicityMCpred::HasRecVertex() {

  double fMaxDeltaSpdTrackAbsolute = 0.5f;
  double fMaxDeltaSpdTrackNsigmaSPD = 1.e14f;
  double fMaxDeltaSpdTrackNsigmaTrack = 1.e14;
  double fMaxResolutionSPDvertex = 0.25f;
  double fMaxDispersionSPDvertex = 1.e14f;

  bool fRequireTrackVertex = true;
  unsigned long fFlag;
  fFlag = BIT(AliEventCuts::kNoCuts);

  const AliVVertex *vtTrc = fESD->GetPrimaryVertex();
  bool isTrackV = true;
  if (vtTrc->IsFromVertexer3D() || vtTrc->IsFromVertexerZ())
    isTrackV = false;
  const AliVVertex *vtSPD = fESD->GetPrimaryVertexSPD();

  if (vtSPD->GetNContributors() > 0)
    fFlag |= BIT(AliEventCuts::kVertexSPD);

  if (vtTrc->GetNContributors() > 1 && isTrackV)
    fFlag |= BIT(AliEventCuts::kVertexTracks);

  if (((fFlag & BIT(AliEventCuts::kVertexTracks)) || !fRequireTrackVertex) &&
      (fFlag & BIT(AliEventCuts::kVertexSPD)))
    fFlag |= BIT(AliEventCuts::kVertex);

  const AliVVertex *&vtx =
      bool(fFlag & BIT(AliEventCuts::kVertexTracks)) ? vtTrc : vtSPD;
  AliVVertex *fPrimaryVertex = const_cast<AliVVertex *>(vtx);
  if (!fPrimaryVertex)
    return kFALSE;

  /// Vertex quality cuts
  double covTrc[6], covSPD[6];
  vtTrc->GetCovarianceMatrix(covTrc);
  vtSPD->GetCovarianceMatrix(covSPD);
  double dz = bool(fFlag & AliEventCuts::kVertexSPD) &&
                      bool(fFlag & AliEventCuts::kVertexTracks)
                  ? vtTrc->GetZ() - vtSPD->GetZ()
                  : 0.; /// If one of the two vertices is not available this cut
                        /// is always passed.
  double errTot = TMath::Sqrt(covTrc[5] + covSPD[5]);
  double errTrc =
      bool(fFlag & AliEventCuts::kVertexTracks) ? TMath::Sqrt(covTrc[5]) : 1.;
  double nsigTot = TMath::Abs(dz) / errTot, nsigTrc = TMath::Abs(dz) / errTrc;
  /// vertex dispersion for run1, only for ESD, AOD code to be added here
  const AliESDVertex *vtSPDESD = dynamic_cast<const AliESDVertex *>(vtSPD);
  double vtSPDdispersion = vtSPDESD ? vtSPDESD->GetDispersion() : 0;
  if ((TMath::Abs(dz) <= fMaxDeltaSpdTrackAbsolute &&
       nsigTot <= fMaxDeltaSpdTrackNsigmaSPD &&
       nsigTrc <=
           fMaxDeltaSpdTrackNsigmaTrack) && // discrepancy track-SPD vertex
      (!vtSPD->IsFromVertexerZ() ||
       TMath::Sqrt(covSPD[5]) <= fMaxResolutionSPDvertex) &&
      (!vtSPD->IsFromVertexerZ() ||
       vtSPDdispersion <= fMaxDispersionSPDvertex) /// vertex dispersion cut for
                                                   /// run1, only for ESD
      ) // quality cut on vertexer SPD z
    fFlag |= BIT(AliEventCuts::kVertexQuality);

  bool hasVtx = (TESTBIT(fFlag, AliEventCuts::kVertex)) &&
                (TESTBIT(fFlag, AliEventCuts::kVertexQuality));
  fVtxz = vtSPD->GetZ();
  return hasVtx;
}
int AliAnalysisTaskFlatenicityMCpred::FillMCarray(vector<double> &ptArray,
                                                  vector<double> &yArray,
                                                  vector<double> &etaArray,
                                                  vector<int> &idArray) {

  ptArray.clear();
  yArray.clear();
  etaArray.clear();
  idArray.clear();

  int nParticlesChg = 0;

  for (int i = 0; i < fMC->GetNumberOfTracks(); ++i) {
    AliMCParticle *particle = (AliMCParticle *)fMC->GetTrack(i);
    if (!particle)
      continue;
    if (!fMC->IsPhysicalPrimary(i))
      continue;
    if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i, fMC))
      continue;
    if (TMath::Abs(particle->Eta()) > 2.0)
      continue;
    if (particle->Pt() <= 0)
      continue;
    if (TMath::Abs(particle->Charge()) < 0.1)
      continue;

    int idPart = -1;
    int partPDG = TMath::Abs(particle->PdgCode());
    if (partPDG == 211)
      idPart = 0; // pions
    else if (partPDG == 321)
      idPart = 1; // kaons
    else if (partPDG == 2212)
      idPart = 2; // protons
    else
      idPart = 3; // rest of the charged particles

    ptArray.push_back(particle->Pt());
    yArray.push_back(particle->Y());
    etaArray.push_back(particle->Eta());
    idArray.push_back(idPart);

    nParticlesChg++;
  }

  return nParticlesChg;
}
//___________________________________________________________________________
void AliAnalysisTaskFlatenicityMCpred::FillHistos1(int multGen,
                                                   const vector<double> &ptGen,
                                                   const vector<double> &yGen,
                                                   const vector<double> &etaGen,
                                                   const vector<int> &idGen) {
  hFlatV0MPerc[fV0Mindex]->Fill(flatenicity_m);
  if (multGen < 1)
    return;
  bool isINEL0 = false;
  int nPart = 0;
  for (int i = 0; i < multGen; ++i) {
    if (abs(etaGen[i]) < 1) {
      nPart++;
    }
  }

  if (nPart > 0) {
    isINEL0 = true;
  }
  if (isINEL0) {
    for (int i = 0; i < multGen; ++i) {

      if (idGen[i] >= 0 && idGen[i] < 3) {
        if (TMath::Abs(yGen[i]) <= fEtaCut) {
          hFlatVsPt1V0MPerc[idGen[i] + 1][fV0Mindex]->Fill(flatenicity_m,
                                                           ptGen[i]);
          hFlatVsPt1[idGen[i] + 1]->Fill(flatenicity_m, ptGen[i]);
          hMultVsFlatVsPt1[idGen[i] + 1]->Fill(fv0multalice, flatenicity_m,
                                               ptGen[i]);
        }
      }

      if (TMath::Abs(etaGen[i]) > fEtaCut) {
        continue;
      }

      hFlatVsPt1V0MPerc[0][fV0Mindex]->Fill(flatenicity_m, ptGen[i]);
      hFlatVsPt1[0]->Fill(flatenicity_m, ptGen[i]);
      hMultVsFlatVsPt1[0]->Fill(fv0multalice, flatenicity_m, ptGen[i]);
    }
  }
}
//________________________________________________
