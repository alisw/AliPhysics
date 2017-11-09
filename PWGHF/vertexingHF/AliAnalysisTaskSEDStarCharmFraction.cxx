/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

//*************************************************************************
// Class AliAnalysisTaskSEDStarCharmFraction
// AliAnalysisTask for the extraction of the fraction of prompt charm for D*
// using the impact parameter to the primary vertex
// 
// Author: Jasper van der Maarel <J.vanderMaarel@uu.nl>
//*************************************************************************

#include <TSystem.h>
#include <TParticle.h>
#include <TH1D.h>
#include <THnSparse.h>
#include <TTree.h>
#include "TROOT.h"
#include <TDatabasePDG.h>
#include <TParameter.h>
#include <AliAnalysisDataSlot.h>
#include <AliAnalysisDataContainer.h>
#include "AliRDHFCutsDStartoKpipi.h"
#include "AliMCEvent.h"
#include "AliAnalysisManager.h"
#include "AliAODMCHeader.h"
#include "AliAODHandler.h"
#include "AliLog.h"
#include "AliAODVertex.h"
#include "AliAODRecoDecay.h"
#include "AliAODRecoDecayHF.h"
#include "AliAODRecoCascadeHF.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAnalysisVertexingHF.h"
#include "AliESDtrack.h"
#include "AliAODMCParticle.h"
#include "AliNormalizationCounter.h"
#include "AliAODEvent.h"
#include "AliVertexerTracks.h"
#include "AliNeutralTrackParam.h"
#include "AliAnalysisTaskSEDStarCharmFraction.h"

ClassImp(AliAnalysisTaskSEDStarCharmFraction)

AliAnalysisTaskSEDStarCharmFraction::AliAnalysisTaskSEDStarCharmFraction() : 
  AliAnalysisTaskSE(),
  fCuts(0),
  fCounter(0),
  fReadMC(kFALSE),
  fSkipHijing(kTRUE),
  fSingleSideband(kFALSE),
  fImpParCut(0),
  fPDGMDStarD0(0),
  fNPtBins(0),
  fNEvents(0),
  fTreeCandidate(0),
  fListCandidate(0),
  fListSignal(0),
  fListSignalPrompt(0),
  fListSignalFromB(0),
  fListBackground(0),
  fIsSideband(kFALSE),
  fIsPeak(kFALSE),
  fIsSignal(kFALSE),
  fIsSignalPrompt(kFALSE),
  fIsSignalFromB(kFALSE),
  fIsBackground(kFALSE),
  fNewPrimVtx(0),
  fDStarVtx(0),
  fMagneticField(0),
  fPtForTree(0),
  fInvMassForTree(0),
  fImpParForTree(0),
  fTrueImpParForTree(0),
  fTypeForTree(0),
  fSignalTypeForTree(0)
{ // Default constructor
}

AliAnalysisTaskSEDStarCharmFraction::AliAnalysisTaskSEDStarCharmFraction(const char *name, AliRDHFCutsDStartoKpipi *cuts) : 
  AliAnalysisTaskSE(name),
  fCuts(0),
  fCounter(0),
  fReadMC(kFALSE),
  fSkipHijing(kTRUE),
  fSingleSideband(kFALSE),
  fImpParCut(0),
  fPDGMDStarD0(0),
  fNPtBins(0),
  fNEvents(0),
  fTreeCandidate(0),
  fListCandidate(0),
  fListSignal(0),
  fListSignalPrompt(0),
  fListSignalFromB(0),
  fListBackground(0),
  fIsSideband(kFALSE),
  fIsPeak(kFALSE),
  fIsSignal(kFALSE),
  fIsSignalPrompt(kFALSE),
  fIsSignalFromB(kFALSE),
  fIsBackground(kFALSE),
  fNewPrimVtx(0),
  fDStarVtx(0),
  fMagneticField(0),
  fPtForTree(0),
  fInvMassForTree(0),
  fImpParForTree(0),
  fTrueImpParForTree(0),
  fTypeForTree(0),
  fSignalTypeForTree(0)
{ // Constructor

  if (fCuts) {
    delete fCuts;
    fCuts = 0;
  }

  fCuts = new AliRDHFCutsDStartoKpipi(*cuts);
  fNPtBins = fCuts->GetNPtBins();

  DefineOutput(1, TH1D::Class());                     // Event counter
  DefineOutput(2, TList::Class());                    // Candidate
  DefineOutput(3, TList::Class());                    // Signal
  DefineOutput(4, TList::Class());                    // SignalPrompt
  DefineOutput(5, TList::Class());                    // SignalFromB
  DefineOutput(6, TList::Class());                    // Background
  DefineOutput(7, AliRDHFCutsDStartoKpipi::Class());  // Cuts
  DefineOutput(8, AliNormalizationCounter::Class());  // Normalization
  DefineOutput(9, TTree::Class());                    // Tree for candidates

  for (Int_t i=0;i<30;i++) {
    fPeakCut[i] = 0; fSidebandCut[i] = 0; fSidebandWindow[i] = 0;
  }

  //DefineOutput(5, AliNormalizationCounter::Class());
}

AliAnalysisTaskSEDStarCharmFraction::~AliAnalysisTaskSEDStarCharmFraction()
{ // Destructor 
  if (fCuts) {
    delete fCuts;
    fCuts = 0;
  }

  if (fNewPrimVtx) {
    delete fNewPrimVtx;
    fNewPrimVtx = 0;
  }

  if (fDStarVtx) {
    delete fDStarVtx;
    fDStarVtx = 0;
  }
}

void AliAnalysisTaskSEDStarCharmFraction::Init()
{ // Initialization

  AliRDHFCutsDStartoKpipi *copyCuts = new AliRDHFCutsDStartoKpipi(*fCuts);
  const char* name = GetOutputSlot(7)->GetContainer()->GetName();
  copyCuts->SetName(name);

  PostData(7, copyCuts);

  fPDGMDStarD0 = TDatabasePDG::Instance()->GetParticle(413)->Mass() - TDatabasePDG::Instance()->GetParticle(421)->Mass();
}

void AliAnalysisTaskSEDStarCharmFraction::UserCreateOutputObjects()
{ // Create histograms
  fNEvents = new TH1D(GetOutputSlot(1)->GetContainer()->GetName(), "Counter", 21, -0.5, 20.5);
  fNEvents->GetXaxis()->SetBinLabel(1, "Events");
  fNEvents->GetXaxis()->SetBinLabel(2, "Magnetic field");
  fNEvents->GetXaxis()->SetBinLabel(3, "Selected");
  fNEvents->GetXaxis()->SetBinLabel(4, "Primary vertex");
  fNEvents->GetXaxis()->SetBinLabel(5, "Candidates");
  fNEvents->GetXaxis()->SetBinLabel(6, "Passed track cuts");
  fNEvents->GetXaxis()->SetBinLabel(7, "Passed fiducial acc cuts");
  fNEvents->GetXaxis()->SetBinLabel(8, "Passed cand cuts");
  fNEvents->GetXaxis()->SetBinLabel(9, "Recalc prim vtx");
  fNEvents->GetXaxis()->SetBinLabel(10, "Calc D* vtx");
  fNEvents->GetXaxis()->SetBinLabel(11, "Pass daught imppar cut");
  fNEvents->GetXaxis()->SetBinLabel(12, "Real D* (MC)");
  fNEvents->GetXaxis()->SetBinLabel(13, "Skip from Hijing (MC)");

  fListCandidate = new TList();
  fListCandidate->SetOwner();
  fListCandidate->SetName("listCandidate");
  SetUpList(fListCandidate);

  fListSignal = new TList();
  fListSignal->SetOwner();
  fListSignal->SetName("listSignal");
  if (fReadMC) {
    SetUpList(fListSignal);
  }

  fListSignalPrompt = new TList();
  fListSignalPrompt->SetOwner();
  fListSignalPrompt->SetName("listSignalPrompt");
  if (fReadMC) {
    SetUpList(fListSignalPrompt);
  }

  fListSignalFromB = new TList();
  fListSignalFromB->SetOwner();
  fListSignalFromB->SetName("listSignalFromB");
  if (fReadMC) {
    SetUpList(fListSignalFromB);
  }

  fListBackground = new TList();
  fListBackground->SetOwner();
  fListBackground->SetName("listBackground");
  if (fReadMC) {
    SetUpList(fListBackground);
  }

  fCounter = new AliNormalizationCounter(GetOutputSlot(8)->GetContainer()->GetName());
  fCounter->Init();

  fTreeCandidate = new TTree(GetOutputSlot(9)->GetContainer()->GetName(), GetOutputSlot(9)->GetContainer()->GetName());
  fTreeCandidate->Branch("pt", &fPtForTree, "pt/D");
  fTreeCandidate->Branch("M", &fInvMassForTree, "M/D");
  fTreeCandidate->Branch("d0", &fImpParForTree, "d0/D");
  fTreeCandidate->Branch("trued0", &fTrueImpParForTree, "trued0/D");
  fTreeCandidate->Branch("type", &fTypeForTree, "type/S");
  fTreeCandidate->Branch("stype", &fSignalTypeForTree, "stype/S");

  PostData(1, fNEvents);
  PostData(2, fListCandidate);
  PostData(3, fListSignal);
  PostData(4, fListSignalPrompt);
  PostData(5, fListSignalFromB);
  PostData(6, fListBackground);
  PostData(8, fCounter);
  PostData(9, fTreeCandidate);
}

void AliAnalysisTaskSEDStarCharmFraction::SetUpList(TList *list)
{ // Fill a TList with histograms

  TString listName = list->GetName();
  listName.ReplaceAll("list", "");
  listName.ToLower();
  if (listName.EqualTo("signalfromb")) {
    listName = "signal from B";
  }
  if (listName.EqualTo("signalprompt")) {
    listName = "signal prompt";
  }

  TString regions[3] = {"", "Peak", "Sideband"}; // Invariant mass regions: all, peak region, sideband region

  // DeltaInvMass
  TH1D* hDeltaInvMass = new TH1D("DeltaInvMass", Form("D*-D^{0} invariant mass, %s; #DeltaM [GeV/c^{2}]; Entries", listName.Data()), 700, 0.13, 0.2);
  hDeltaInvMass->Sumw2();
  hDeltaInvMass->SetLineColor(4);
  hDeltaInvMass->SetMarkerColor(4);
  hDeltaInvMass->SetMarkerStyle(20);
  hDeltaInvMass->SetMarkerSize(0.6);
  list->Add(hDeltaInvMass);

  // DeltaInvMassPre
  TH1D* hDeltaInvMassPre = new TH1D("DeltaInvMassPre", Form("D*-D^{0} invariant mass pre vtx, %s; #DeltaM [GeV/c^{2}]; Entries", listName.Data()), 700, 0.13, 0.2);
  hDeltaInvMassPre->Sumw2();
  hDeltaInvMassPre->SetLineColor(4);
  hDeltaInvMassPre->SetMarkerColor(4);
  hDeltaInvMassPre->SetMarkerStyle(20);
  hDeltaInvMassPre->SetMarkerSize(0.6);
  list->Add(hDeltaInvMassPre);

  // InvMassD0
  TH1D* hInvMassD0 = new TH1D("InvMassD0", Form("D^{0} invariant mass, %s; M [GeV/c^{2}]; Entries", listName.Data()), 600, 1.6, 2.2);
  hInvMassD0->Sumw2();
  hInvMassD0->SetLineColor(4);
  hInvMassD0->SetMarkerColor(4);
  hInvMassD0->SetMarkerStyle(20);
  hInvMassD0->SetMarkerSize(0.6);
  list->Add(hInvMassD0);

  for (Int_t r=0;r<3;r++) {
    TString regionTitle(regions[r]);
    if (!regionTitle.EqualTo("")) {
      regionTitle.ToLower();
      regionTitle.Prepend(", ");
      regionTitle.Append(" region");
    }

    // ImpPar
    TH1D* hImpPar = new TH1D(Form("ImpPar%s", regions[r].Data()), Form("D* impact parameter%s, %s; Impact parameter [#mum]; Entries", regionTitle.Data(), listName.Data()), 1000, -1000., 1000.);
    hImpPar->Sumw2();
    hImpPar->SetLineColor(4);
    hImpPar->SetMarkerColor(4);
    hImpPar->SetMarkerStyle(20);
    hImpPar->SetMarkerSize(0.6);
    list->Add(hImpPar);

    // ImpParSoftPi
    TH1D* hImpParSoftPi = new TH1D(Form("ImpParSoftPi%s", regions[r].Data()), Form("#pi_{s} impact parameter%s, %s; Impact parameter [#mum]; Entries", regionTitle.Data(), listName.Data()), 1000, -1000., 1000.);
    hImpParSoftPi->Sumw2();
    hImpParSoftPi->SetLineColor(4);
    hImpParSoftPi->SetMarkerColor(4);
    hImpParSoftPi->SetMarkerStyle(20);
    hImpParSoftPi->SetMarkerSize(0.6);
    list->Add(hImpParSoftPi);

    // ImpParD0
    TH1D* hImpParD0 = new TH1D(Form("ImpParD0%s", regions[r].Data()), Form("D^{0} impact parameter%s, %s; Impact parameter [#mum]; Entries", regionTitle.Data(), listName.Data()), 1000, -1000., 1000.);
    hImpParD0->Sumw2();
    hImpParD0->SetLineColor(4);
    hImpParD0->SetMarkerColor(4);
    hImpParD0->SetMarkerStyle(20);
    hImpParD0->SetMarkerSize(0.6);
    list->Add(hImpParD0);

    //ImpParD0K
    TH1D* hImpParD0K = new TH1D(Form("ImpParD0K%s", regions[r].Data()), Form("K (from D^{0}) impact parameter%s, %s; Impact parameter [#mum]; Entries", regionTitle.Data(), listName.Data()), 1000, -1000., 1000.);
    hImpParD0K->Sumw2();
    hImpParD0K->SetLineColor(4);
    hImpParD0K->SetMarkerColor(4);
    hImpParD0K->SetMarkerStyle(20);
    hImpParD0K->SetMarkerSize(0.6);
    list->Add(hImpParD0K);

    //ImpParD0K2
    TH1D* hImpParD0K2 = new TH1D(Form("ImpParD0K2%s", regions[r].Data()), Form("K (from D^{0}) impact parameter (no recalc prim vtx)%s, %s; Impact parameter [#mum]; Entries", regionTitle.Data(), listName.Data()), 1000, -1000., 1000.);
    hImpParD0K2->Sumw2();
    hImpParD0K2->SetLineColor(4);
    hImpParD0K2->SetMarkerColor(4);
    hImpParD0K2->SetMarkerStyle(20);
    hImpParD0K2->SetMarkerSize(0.6);
    list->Add(hImpParD0K2);

    //ImpParD0Pi
    TH1D* hImpParD0Pi = new TH1D(Form("ImpParD0Pi%s", regions[r].Data()), Form("#pi (from D^{0}) impact parameter%s, %s; Impact parameter [#mum]; Entries", regionTitle.Data(), listName.Data()), 1000, -1000., 1000.);
    hImpParD0Pi->Sumw2();
    hImpParD0Pi->SetLineColor(4);
    hImpParD0Pi->SetMarkerColor(4);
    hImpParD0Pi->SetMarkerStyle(20);
    hImpParD0Pi->SetMarkerSize(0.6);
    list->Add(hImpParD0Pi);

    //ImpParD0Pi2
    TH1D* hImpParD0Pi2 = new TH1D(Form("ImpParD0Pi2%s", regions[r].Data()), Form("#pi (from D^{0}) impact parameter (no recalc prim vtx)%s, %s; Impact parameter [#mum]; Entries", regionTitle.Data(), listName.Data()), 1000, -1000., 1000.);
    hImpParD0Pi2->Sumw2();
    hImpParD0Pi2->SetLineColor(4);
    hImpParD0Pi2->SetMarkerColor(4);
    hImpParD0Pi2->SetMarkerStyle(20);
    hImpParD0Pi2->SetMarkerSize(0.6);
    list->Add(hImpParD0Pi2);
  }

  Float_t *ptLimits = fCuts->GetPtBinLimits();

  for (Int_t i=0;i<fNPtBins;i++) {
    Float_t ptLow = ptLimits[i];
    Float_t ptHigh = ptLimits[(i+1)];

    TH1D* hDeltaInvMassLoop = new TH1D(Form("DeltaInvMass_%d", i), Form("D*-D^{0} invariant mass, %s, %.1f < p_{T} < %.1f GeV/c; #DeltaM [GeV/c^{2}]; Entries", listName.Data(), ptLow, ptHigh), 700, 0.13, 0.2);
    hDeltaInvMassLoop->Sumw2();
    hDeltaInvMassLoop->SetLineColor(4);
    hDeltaInvMassLoop->SetMarkerColor(4);
    hDeltaInvMassLoop->SetMarkerStyle(20);
    hDeltaInvMassLoop->SetMarkerSize(0.6);
    list->Add(hDeltaInvMassLoop);

    TH1D* hDeltaInvMassPreLoop = new TH1D(Form("DeltaInvMassPre_%d", i), Form("D*-D^{0} invariant mass pre vtx, %s, %.1f < p_{T} < %.1f GeV/c; #DeltaM [GeV/c^{2}]; Entries", listName.Data(), ptLow, ptHigh), 700, 0.13, 0.2);
    hDeltaInvMassPreLoop->Sumw2();
    hDeltaInvMassPreLoop->SetLineColor(4);
    hDeltaInvMassPreLoop->SetMarkerColor(4);
    hDeltaInvMassPreLoop->SetMarkerStyle(20);
    hDeltaInvMassPreLoop->SetMarkerSize(0.6);
    list->Add(hDeltaInvMassPreLoop);

    TH1D* hInvMassD0Loop = new TH1D(Form("InvMassD0_%d", i), Form("D^{0} invariant mass, %s, %.1f < p_{T} < %.1f GeV/c; M [GeV/c^{2}]; Entries", listName.Data(), ptLow, ptHigh), 600, 1.6, 2.2);
    hInvMassD0Loop->Sumw2();
    hInvMassD0Loop->SetLineColor(4);
    hInvMassD0Loop->SetMarkerColor(4);
    hInvMassD0Loop->SetMarkerStyle(20);
    hInvMassD0Loop->SetMarkerSize(0.6);
    list->Add(hInvMassD0Loop);

    for (Int_t r=0;r<3;r++) {
      TString regionTitle(regions[r]);
      if (!regionTitle.EqualTo("")) {
        regionTitle.ToLower();
        regionTitle.Prepend(", ");
        regionTitle.Append(" region");
      }

      //ImpPar
      TH1D* hImpPar = new TH1D(Form("ImpPar%s_%d", regions[r].Data(), i), Form("D* impact parameter%s, %s, %.1f < p_{T} < %.1f GeV/c; Impact parameter [#mum]; Entries", regionTitle.Data(), listName.Data(), ptLow, ptHigh), 1000, -1000., 1000.);
      hImpPar->Sumw2();
      hImpPar->SetLineColor(4);
      hImpPar->SetMarkerColor(4);
      hImpPar->SetMarkerStyle(20);
      hImpPar->SetMarkerSize(0.6);
      list->Add(hImpPar);

      //ImpParSoftPi
      TH1D* hImpParSoftPi = new TH1D(Form("ImpParSoftPi%s_%d", regions[r].Data(), i), Form("#pi_{s} impact parameter%s, %s, %.1f < p_{T} < %.1f GeV/c; Impact parameter [#mum]; Entries", regionTitle.Data(), listName.Data(), ptLow, ptHigh), 1000, -1000., 1000.);
      hImpParSoftPi->Sumw2();
      hImpParSoftPi->SetLineColor(4);
      hImpParSoftPi->SetMarkerColor(4);
      hImpParSoftPi->SetMarkerStyle(20);
      hImpParSoftPi->SetMarkerSize(0.6);
      list->Add(hImpParSoftPi);

      //ImpParD0
      TH1D* hImpParD0 = new TH1D(Form("ImpParD0%s_%d", regions[r].Data(), i), Form("D^{0} impact parameter%s, %s, %.1f < p_{T} < %.1f GeV/c; Impact parameter [#mum]; Entries", regionTitle.Data(), listName.Data(), ptLow, ptHigh), 1000, -1000., 1000.);
      hImpParD0->Sumw2();
      hImpParD0->SetLineColor(4);
      hImpParD0->SetMarkerColor(4);
      hImpParD0->SetMarkerStyle(20);
      hImpParD0->SetMarkerSize(0.6);
      list->Add(hImpParD0);

      //ImpParD0K
      TH1D* hImpParD0K = new TH1D(Form("ImpParD0K%s_%d", regions[r].Data(), i), Form("K (from D^{0}) impact parameter%s, %s, %.1f < p_{T} < %.1f GeV/c; Impact parameter [#mum]; Entries", regionTitle.Data(), listName.Data(), ptLow, ptHigh), 1000, -1000., 1000.);
      hImpParD0K->Sumw2();
      hImpParD0K->SetLineColor(4);
      hImpParD0K->SetMarkerColor(4);
      hImpParD0K->SetMarkerStyle(20);
      hImpParD0K->SetMarkerSize(0.6);
      list->Add(hImpParD0K);

      //ImpParD0K2
      TH1D* hImpParD0K2 = new TH1D(Form("ImpParD0K2%s_%d", regions[r].Data(), i), Form("K (from D^{0}) impact parameter (no recalc prim vtx)%s, %s, %.1f < p_{T} < %.1f GeV/c; Impact parameter [#mum]; Entries", regionTitle.Data(), listName.Data(), ptLow, ptHigh), 1000, -1000., 1000.);
      hImpParD0K2->Sumw2();
      hImpParD0K2->SetLineColor(4);
      hImpParD0K2->SetMarkerColor(4);
      hImpParD0K2->SetMarkerStyle(20);
      hImpParD0K2->SetMarkerSize(0.6);
      list->Add(hImpParD0K2);

      //ImpParD0Pi
      TH1D* hImpParD0Pi = new TH1D(Form("ImpParD0Pi%s_%d", regions[r].Data(), i), Form("#pi (from D^{0}) impact parameter%s, %s, %.1f < p_{T} < %.1f GeV/c; Impact parameter [#mum]; Entries", regionTitle.Data(), listName.Data(), ptLow, ptHigh), 1000, -1000., 1000.);
      hImpParD0Pi->Sumw2();
      hImpParD0Pi->SetLineColor(4);
      hImpParD0Pi->SetMarkerColor(4);
      hImpParD0Pi->SetMarkerStyle(20);
      hImpParD0Pi->SetMarkerSize(0.6);
      list->Add(hImpParD0Pi);

      //ImpParD0Pi2
      TH1D* hImpParD0Pi2 = new TH1D(Form("ImpParD0Pi2%s_%d", regions[r].Data(), i), Form("#pi (from D^{0}) impact parameter (no recalc prim vtx)%s, %s, %.1f < p_{T} < %.1f GeV/c; Impact parameter [#mum]; Entries", regionTitle.Data(), listName.Data(), ptLow, ptHigh), 1000, -1000., 1000.);
      hImpParD0Pi2->Sumw2();
      hImpParD0Pi2->SetLineColor(4);
      hImpParD0Pi2->SetMarkerColor(4);
      hImpParD0Pi2->SetMarkerStyle(20);
      hImpParD0Pi2->SetMarkerSize(0.6);
      list->Add(hImpParD0Pi2);
    }
  }

  if (list == fListSignalPrompt || list == fListSignalFromB || list == fListBackground) {
    for (Int_t i=0;i<fNPtBins;i++) {
      TH1D* hProdd0d0All = new TH1D(Form("Prodd0d0All_%d", i), Form("d0d0, no cut, %s; d0d0 [cm^{2}]; Entries", listName.Data()), 2000, -0.0002, 0.0002);
      hProdd0d0All->Sumw2();
      hProdd0d0All->SetLineColor(4);
      hProdd0d0All->SetMarkerColor(4);
      hProdd0d0All->SetMarkerStyle(20);
      hProdd0d0All->SetMarkerSize(0.6);
      list->Add(hProdd0d0All);

      TH1D* hCosPntAngleAll = new TH1D(Form("CosPntAngleAll_%d", i), Form("cos #theta_{point}, no cut, %s; cos #theta_{point}; Entries", listName.Data()), 2000, -1, 1);
      hCosPntAngleAll->Sumw2();
      hCosPntAngleAll->SetLineColor(4);
      hCosPntAngleAll->SetMarkerColor(4);
      hCosPntAngleAll->SetMarkerStyle(20);
      hCosPntAngleAll->SetMarkerSize(0.6);
      list->Add(hCosPntAngleAll);

      TH1D* hNormDecLenAll = new TH1D(Form("NormDecLenAll_%d", i), Form("Normalized decay length, no cut, %s; Normalized decay length; Entries", listName.Data()), 2000, 0, 20);
      hNormDecLenAll->Sumw2();
      hNormDecLenAll->SetLineColor(4);
      hNormDecLenAll->SetMarkerColor(4);
      hNormDecLenAll->SetMarkerStyle(20);
      hNormDecLenAll->SetMarkerSize(0.6);
      list->Add(hNormDecLenAll);
    }
  }

  // Mother and grandmother particle IDs for D* daughters, background, peak region
  if (list == fListBackground) {
    Int_t sparseBins[3] = {1000, 1000, 1000}; Double_t sparseLow[3] = {-0.5, -0.5, -0.5}; Double_t sparseHigh[3] = {999.5, 999.5, 999.5};
    THnSparseF* hPDG = new THnSparseF("PDGMotherPeak", "PDG mother, peak region, background", 3, sparseBins, sparseLow, sparseHigh);
    list->Add(hPDG);
  }

  // Store true impact parameter distribution
  if (list == fListSignalFromB) {
    for (Int_t i=0;i<fNPtBins;i++) {
      Float_t ptLow = ptLimits[i];
      Float_t ptHigh = ptLimits[(i+1)];

      TH1D* hImpParTrue = new TH1D(Form("ImpParTrue_%d", i), Form("True D* impact parameter, %s, %.1f < p_{T} < %.1f GeV/c; Impact parameter [#mum]; Entries", listName.Data(), ptLow, ptHigh), 1000, -1000., 1000.);
      hImpParTrue->Sumw2();
      hImpParTrue->SetLineColor(4);
      hImpParTrue->SetMarkerColor(4);
      hImpParTrue->SetMarkerStyle(20);
      hImpParTrue->SetMarkerSize(0.6);
      list->Add(hImpParTrue);
    }
  }

  char mergeMode(102); // f

  TParameter<bool> *pSingleSideband = new TParameter<bool>("SingleSideband", fSingleSideband, mergeMode);
  list->Add(pSingleSideband);

  for (Int_t i=0;i<fNPtBins;i++) {
    TParameter<double> *pPeakCut = new TParameter<double>(Form("PeakCut_%d", i), fPeakCut[i], mergeMode);
    list->Add(pPeakCut);

    TParameter<double> *pSidebandCut = new TParameter<double>(Form("SidebandCut_%d", i), fSidebandCut[i], mergeMode);
    list->Add(pSidebandCut);

    TParameter<double> *pSidebandWindow = new TParameter<double>(Form("SidebandWindow_%d", i), fSidebandWindow[i], mergeMode);
    list->Add(pSidebandWindow);
  }
}

void AliAnalysisTaskSEDStarCharmFraction::UserExec(Option_t */*option*/)
{ // Execute analysis for current event
  if (!fInputEvent) {
    Error("UserExec", "NO EVENT FOUND!");
    return;
  }

  AliAODEvent *aodEvent = dynamic_cast<AliAODEvent*>(fInputEvent);
  TClonesArray *arrayCand = 0;

  fNEvents->Fill(0);

  if (!aodEvent && AODEvent() && IsStandardAOD()) {
    aodEvent = dynamic_cast<AliAODEvent*> (AODEvent());
    AliAODHandler* aodHandler = (AliAODHandler*) ((AliAnalysisManager::GetAnalysisManager())->GetOutputEventHandler());
    if(aodHandler->GetExtensions()) {
      AliAODExtension *ext = (AliAODExtension*) aodHandler->GetExtensions()->FindObject("AliAOD.VertexingHF.root");
      AliAODEvent *aodFromExt = ext->GetAOD();
      arrayCand = (TClonesArray*) aodFromExt->GetList()->FindObject("Dstar");
    }
  }
  else {
    arrayCand = (TClonesArray*) aodEvent->GetList()->FindObject("Dstar");
  }

  Int_t pdgDgDStartoD0pi[2] = {421, 211};
  Int_t pdgDgD0toKpi[2] = {321, 211};

  TClonesArray *mcArray = 0;
  AliAODMCHeader *mcHeader = 0;
  if (fReadMC) {
    mcArray = dynamic_cast<TClonesArray*>(aodEvent->FindListObject(AliAODMCParticle::StdBranchName()));
    if (!mcArray) {
      AliError("Could not find Monte-Carlo in AOD");
      return;
    }
    // Load MC header
    mcHeader = (AliAODMCHeader*) aodEvent->GetList()->FindObject(AliAODMCHeader::StdBranchName());
    if (!mcHeader) {
      printf("AliAnalysisTaskSEDplus::UserExec: MC header branch not found!\n");
      return;
    }
  }

  if (!aodEvent->GetPrimaryVertex() || TMath::Abs(aodEvent->GetMagneticField()) < 0.001) {
    return;
  }
  fNEvents->Fill(1);

  fCounter->StoreEvent(aodEvent, fCuts, fReadMC);

  fMagneticField = aodEvent->GetMagneticField();

  if (!fCuts->IsEventSelected(aodEvent)) {
    return;
  }
  fNEvents->Fill(2);
  
  AliAODVertex *vtx1 = (AliAODVertex*) aodEvent->GetPrimaryVertex();
  if (!vtx1 || vtx1->GetNContributors() < 1) {
    return;
  }
  fNEvents->Fill(3);

  if (!arrayCand) {
    AliInfo("Could not find array of HF vertices, skipping the event");
    return;
  }
  else {
    AliDebug(2, Form("Found %d vertices", arrayCand->GetEntriesFast()));
  }

  Int_t ptBin;
  Double_t pt;

  Int_t nSelectedProd = 0;
  Int_t nSelectedAna = 0;

  Int_t nCand = arrayCand->GetEntriesFast();
  for (Int_t iCand = 0; iCand<nCand; iCand++) {
    // AliAODVertex *origVtx;

    fIsSignal = kFALSE;
    fIsSignalPrompt = kFALSE;
    fIsSignalFromB = kFALSE;
    fIsBackground = kFALSE;

    AliAODRecoCascadeHF* cand = (AliAODRecoCascadeHF*) arrayCand->At(iCand);
    AliAODRecoDecayHF2Prong *candD0 = cand->Get2Prong();
    AliAODTrack *candSoftPi = cand->GetBachelor();

    ptBin = fCuts->PtBin(cand->Pt());

    fNEvents->Fill(4);

    if (ptBin < 0 || ptBin >= fNPtBins) {
      continue;
    }

    // Save some calculation effort by putting D* invariant mass cut
    Double_t deltaInvMass = cand->DeltaInvMass();
    if (TMath::Abs(deltaInvMass-fPDGMDStarD0) > 0.03) {
      continue;
    }

    Int_t isTkSelected = fCuts->IsSelected(cand, AliRDHFCuts::kTracks);
    if (!isTkSelected) {
      continue;
    }
    fNEvents->Fill(5);

    if (!fCuts->IsInFiducialAcceptance(cand->Pt(), cand->YDstar())) {
      continue;
    }
    fNEvents->Fill(6);

    nSelectedProd++;
    nSelectedAna++;

    // Match candidate to MC
    AliAODMCParticle *partDSt = 0;
    if (fReadMC) {
      Int_t mcLabel = cand->MatchToMC(413, 421, pdgDgDStartoD0pi, pdgDgD0toKpi, mcArray);
      if (mcLabel >= 0) {
        fIsSignal = kTRUE;

        fNEvents->Fill(11);

        partDSt = (AliAODMCParticle*) mcArray->At(mcLabel);

        if (fSkipHijing && IsFromHijing(mcArray, partDSt)) {
          fNEvents->Fill(12);
          continue;
        }

        if (IsFromB(mcArray, partDSt)) {
          fIsSignalFromB = kTRUE;
        }
        else {
          fIsSignalPrompt = kTRUE;
        }
      }
      else {
        fIsBackground = kTRUE;
      }
    }

    // Fill histograms for topological variables before applying cuts
    if (fReadMC) {
      Double_t prodd0d0 = candD0->Prodd0d0();
      Double_t cosPntAngle = candD0->CosPointingAngle();
      Double_t normDecLen = candD0->NormalizedDecayLengthXY()*(candD0->P()/candD0->Pt());

      if (fIsSignalPrompt) {
        ((TH1D*) fListSignalPrompt->FindObject(Form("Prodd0d0All_%d", ptBin)))->Fill(prodd0d0);
        ((TH1D*) fListSignalPrompt->FindObject(Form("CosPntAngleAll_%d", ptBin)))->Fill(cosPntAngle);
        ((TH1D*) fListSignalPrompt->FindObject(Form("NormDecLenAll_%d", ptBin)))->Fill(normDecLen);
      }

      if (fIsSignalFromB) {
        ((TH1D*) fListSignalFromB->FindObject(Form("Prodd0d0All_%d", ptBin)))->Fill(prodd0d0);
        ((TH1D*) fListSignalFromB->FindObject(Form("CosPntAngleAll_%d", ptBin)))->Fill(cosPntAngle);
        ((TH1D*) fListSignalFromB->FindObject(Form("NormDecLenAll_%d", ptBin)))->Fill(normDecLen);
      }

      if (fIsBackground) {
        ((TH1D*) fListBackground->FindObject(Form("Prodd0d0All_%d", ptBin)))->Fill(prodd0d0);
        ((TH1D*) fListBackground->FindObject(Form("CosPntAngleAll_%d", ptBin)))->Fill(cosPntAngle);
        ((TH1D*) fListBackground->FindObject(Form("NormDecLenAll_%d", ptBin)))->Fill(normDecLen);
      }
    }

    Int_t isSelected = fCuts->IsSelected(cand, AliRDHFCuts::kCandidate, aodEvent);
    if (!isSelected) {
      continue;
    }
    fNEvents->Fill(7);

    //Fill invariant mass pre vtx calculations
    FillHistogram("DeltaInvMassPre", deltaInvMass);
    FillHistogram(Form("DeltaInvMassPre_%d", ptBin), deltaInvMass);

    // Calculate new primary vertex without D* daughters
    fNewPrimVtx = RemoveDaughtersFromPrimaryVtx(aodEvent, cand);
    if (!fNewPrimVtx) {
      delete fNewPrimVtx; fNewPrimVtx = 0;
      continue;
    }

    /*if (cand->GetOwnPrimaryVtx()) {
      origVtx = new AliAODVertex(*cand->GetOwnPrimaryVtx());
    }
    cand->SetOwnPrimaryVtx(fNewPrimVtx);*/

    fNEvents->Fill(8);

    //D* vertexing
    fDStarVtx = ReconstructDStarVtx(cand);
    if (!fDStarVtx) {
      delete fNewPrimVtx; fNewPrimVtx = 0;
      continue;
    }

    fNEvents->Fill(9);

    Double_t d0, d0Err;
    if (CalculateImpactParameter(candSoftPi, d0, d0Err)) {
      d0 = TMath::Abs(d0*1.e4);
      if (d0 < fImpParCut) {
        delete fNewPrimVtx; fNewPrimVtx = 0;
        delete fDStarVtx; fDStarVtx = 0;
        continue;
      }
    }
    if (CalculateImpactParameter((AliAODTrack*) candD0->GetDaughter(0), d0, d0Err)) {
      d0 = TMath::Abs(d0*1.e4);
      if (d0 < fImpParCut) {
        delete fNewPrimVtx; fNewPrimVtx = 0;
        delete fDStarVtx; fDStarVtx = 0;
        continue;
      }
    }
    if (CalculateImpactParameter((AliAODTrack*) candD0->GetDaughter(1), d0, d0Err)) {
      d0 = TMath::Abs(d0*1.e4);
      if (d0 < fImpParCut) {
        delete fNewPrimVtx; fNewPrimVtx = 0;
        delete fDStarVtx; fDStarVtx = 0;
        continue;
      }
    }

    fNEvents->Fill(10);

    fTrueImpParForTree = fReadMC && fIsSignal ? CalculateTrueImpactParameterDStar(mcHeader, mcArray, cand) : -9998.;

    CheckInvMassDStar(cand);
    FillHistograms(cand);

    if (fReadMC && fIsBackground && fIsPeak) {
      AliAODTrack *daughters[3]; //0=soft pion, 1=D0 pion, 2=D0 kaon
      daughters[0] = candSoftPi;

      if (((AliAODTrack*) candD0->GetDaughter(0))->Charge()*cand->Charge() > 0) {
        // D* charge and first daughter charge are equal in sign <=> first daughter = pion
        daughters[1] = (AliAODTrack*) candD0->GetDaughter(0);
        daughters[2] = (AliAODTrack*) candD0->GetDaughter(1);
      }
      else {
        daughters[1] = (AliAODTrack*) candD0->GetDaughter(1);
        daughters[2] = (AliAODTrack*) candD0->GetDaughter(0);
      }

      Double_t PDGCodes[3] = {0., 0., 0.};

      for (Int_t iD=0; iD<3; iD++) {
        AliAODTrack *daughter = daughters[iD];

        Int_t daughterLabel = daughter->GetLabel();
        if (daughterLabel >= 0) {
          AliAODMCParticle *daughterMC = (AliAODMCParticle*) mcArray->At(daughterLabel);
          Int_t motherLabel = daughterMC->GetMother();
          if (motherLabel >= 0) {
            AliAODMCParticle *motherMC = (AliAODMCParticle*) mcArray->At(motherLabel);
            Int_t motherPDG = TMath::Abs(motherMC->GetPdgCode());

            if (iD == 0) {
              //Soft pion: store ID of mother
              PDGCodes[iD] = motherPDG;
            }
            else {
              //D0 daughters: look at grandmother

              Int_t grandmotherLabel = motherMC->GetMother();
              if (grandmotherLabel >= 0) {
                AliAODMCParticle *grandmotherMC = (AliAODMCParticle*) mcArray->At(grandmotherLabel);
                Int_t grandmotherPDG = TMath::Abs(grandmotherMC->GetPdgCode());
                PDGCodes[iD] = grandmotherPDG;
              }
            }
          }
        }
      }

      ((THnSparseF*) fListBackground->FindObject("PDGMotherPeak"))->Fill(PDGCodes);
    }

    if (fReadMC && fIsSignalFromB) {
      FillTrueImpactParameter(cand);
    }


    /*cand->UnsetOwnPrimaryVtx();
    if (origVtx) {
      cand->SetOwnPrimaryVtx(origVtx);
    }*/

    delete fNewPrimVtx; fNewPrimVtx = 0;
    delete fDStarVtx; fDStarVtx = 0;
    //delete origVtx; origVtx = 0;
  }

  fCounter->StoreCandidates(aodEvent, nSelectedProd, kTRUE);  
  fCounter->StoreCandidates(aodEvent, nSelectedAna, kFALSE); 

  PostData(1, fNEvents);
  PostData(2, fListCandidate);
  PostData(3, fListSignal);
  PostData(4, fListSignalPrompt);
  PostData(5, fListSignalFromB);
  PostData(6, fListBackground);
  PostData(8, fCounter);
  PostData(9, fTreeCandidate);
}

AliAODVertex *AliAnalysisTaskSEDStarCharmFraction::RemoveDaughtersFromPrimaryVtx(AliAODEvent *aod, AliAODRecoCascadeHF *cand)
{ // Determine primary vertex without D* daughters
  AliAODVertex *vtxAOD = aod->GetPrimaryVertex();
  if (!vtxAOD) {
    return 0;
  }

  TString title = vtxAOD->GetTitle();
  if (!title.Contains("VertexerTracks")) {
    return 0;
  }

  AliVertexerTracks *vertexer = new AliVertexerTracks(aod->GetMagneticField());

  vertexer->SetITSMode();
  vertexer->SetMinClusters(3);
  vertexer->SetConstraintOff();

  if (title.Contains("WithConstraint")) {
    Float_t diamondcovxy[3];
    aod->GetDiamondCovXY(diamondcovxy);
    Double_t pos[3] = {aod->GetDiamondX(), aod->GetDiamondY(), 0.};
    Double_t cov[6] = {diamondcovxy[0], diamondcovxy[1], diamondcovxy[2], 0., 0., 10.*10.};
    AliESDVertex *diamond = new AliESDVertex(pos, cov, 1., 1);
    vertexer->SetVtxStart(diamond);
    delete diamond; diamond=NULL;
  }

  Int_t skipped[10]; for(Int_t i=0; i<10; i++) { skipped[i] = -1; }
  Int_t nTrksToSkip = 0, id;
  AliAODTrack *t = 0;

  for (Int_t i=0; i<3; i++) {
    if (i == 0) {
      t = cand->GetBachelor();
    }
    else {
      AliAODRecoDecayHF2Prong *twoProng = cand->Get2Prong();
      t = (AliAODTrack*) twoProng->GetDaughter(i-1);
    }

    id = (Int_t) t->GetID();
    if (id < 0) {
      continue;
    }
    skipped[nTrksToSkip++] = id;
  }

  vertexer->SetSkipTracks(nTrksToSkip, skipped);
  AliESDVertex *vtxESDNew = vertexer->FindPrimaryVertex(aod);

  delete vertexer; vertexer = NULL;

  if (!vtxESDNew) {
    return 0;
  }

  if (vtxESDNew->GetNContributors() <= 0) { 
    delete vtxESDNew; vtxESDNew = NULL;
    return 0;
  }

  // convert to AliAODVertex
  Double_t pos[3], cov[6], chi2perNDF;
  vtxESDNew->GetXYZ(pos); // position
  vtxESDNew->GetCovMatrix(cov); //covariance matrix
  chi2perNDF = vtxESDNew->GetChi2toNDF();
  delete vtxESDNew; vtxESDNew = NULL;

  AliAODVertex *vtxAODNew = new AliAODVertex(pos, cov, chi2perNDF);

  return vtxAODNew;
}

AliAODVertex *AliAnalysisTaskSEDStarCharmFraction::ReconstructDStarVtx(AliAODRecoCascadeHF *cand)
{ // Determine the D* vertex

  // Double_t dcadummy, xdummy, ydummy;
  AliAODRecoDecayHF2Prong *candD0 = cand->Get2Prong();
  AliAODTrack *candSoftPi = cand->GetBachelor();
  AliNeutralTrackParam *ntpD0 = new AliNeutralTrackParam(candD0);
  // ntpD0->PropagateToDCA(fNewPrimVtx, fMagneticField, kVeryBig);
  //AliExternalTrackParam etpSoftPi; etpSoftPi.CopyFromVTrack(candSoftPi);
  AliESDtrack *candSoftPiESD = new AliESDtrack(candSoftPi);
  // dcadummy = candSoftPiESD->GetDCA(ntpD0, fMagneticField, xdummy, ydummy);

  TObjArray *tracks = new TObjArray(2);
  tracks->AddAt(candSoftPiESD, 0);
  tracks->AddAt(ntpD0, 1);
  // tracks->Add(&etpSoftPi);

  Double_t pos[3], cov[6];
  fNewPrimVtx->GetXYZ(pos);
  fNewPrimVtx->GetCovarianceMatrix(cov);
  AliESDVertex *newPrimVtxESD = new AliESDVertex(pos, cov, 100., 100, fNewPrimVtx->GetName());

  AliVertexerTracks *vertexerTracks = new AliVertexerTracks(fMagneticField);
  vertexerTracks->SetVtxStart(newPrimVtxESD);
  AliESDVertex *vertexESD = (AliESDVertex*) vertexerTracks->VertexForSelectedESDTracks(tracks);
  if (!vertexESD) {
    return 0;
  }

  if (vertexESD->GetNContributors() != 2) {
    delete vertexESD; vertexESD = 0;
    return 0;
  }

  Double_t vertRadius2 = vertexESD->GetX()*vertexESD->GetX() + vertexESD->GetY()*vertexESD->GetY();
  if (vertRadius2 > 8.) {
    // vertex outside beam pipe, reject candidate to avoid propagation through material
    delete vertexESD; vertexESD = 0;
    return 0;
  }

  Double_t pos2[3], cov2[6], chi2perNDF;
  vertexESD->GetXYZ(pos2); // position
  vertexESD->GetCovMatrix(cov2); //covariance matrix
  chi2perNDF = vertexESD->GetChi2toNDF();

  delete vertexESD; vertexESD = 0;

  AliAODVertex *vertexAOD = new AliAODVertex(pos2, cov2, chi2perNDF, 0x0, -1, AliAODVertex::kUndef, 2);

  return vertexAOD;
}

void AliAnalysisTaskSEDStarCharmFraction::CheckInvMassDStar(AliAODRecoCascadeHF *cand)
{ // Check if D* candidate falls within peak or sideband region
  Double_t invMassDiff = cand->DeltaInvMass()-fPDGMDStarD0;
  Double_t invMassDiffAbs = TMath::Abs(invMassDiff);

  fIsPeak = kFALSE;
  fIsSideband = kFALSE;

  Int_t ptBin = fCuts->PtBin(cand->Pt());
  if (ptBin < 0) {
    return;
  }

  Double_t peakCut = fPeakCut[ptBin];
  Double_t sidebandCut = fSidebandCut[ptBin];
  Double_t sidebandWindow = fSidebandWindow[ptBin];

  if (invMassDiffAbs < peakCut) {
    fIsPeak = kTRUE;
    return;
  }

  if (!fSingleSideband) {
    invMassDiff = invMassDiffAbs;
  }

  if (invMassDiff > sidebandCut && invMassDiff <= sidebandCut+sidebandWindow) {
    fIsSideband = kTRUE;
  }
}

Bool_t AliAnalysisTaskSEDStarCharmFraction::IsFromB(TClonesArray *arrayMC, const AliAODMCParticle *mcPartCandidate)
{ // Check if the MC particle comes from a B
  Int_t mother = mcPartCandidate->GetMother();
  while (mother > 0) {
    AliAODMCParticle *mcGranma = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother));
    if (mcGranma) {
      Int_t pdgGranma = mcGranma->GetPdgCode();
      Int_t abspdgGranma = TMath::Abs(pdgGranma);
      if ((abspdgGranma > 500 && abspdgGranma < 600) || (abspdgGranma > 5000 && abspdgGranma < 6000)) {
        return kTRUE;
      }
      mother = mcGranma->GetMother();
    }
    else{
      AliError("Failed casting the mother particle!");
      break;
    }
  }

  return kFALSE;
}

Bool_t AliAnalysisTaskSEDStarCharmFraction::IsFromHijing(TClonesArray *arrayMC, const AliAODMCParticle *mcPartCandidate)
{ // Check if the MC particle is from Hijing 
  Int_t mother = mcPartCandidate->GetMother();
  while (mother > 0) {
    AliAODMCParticle *mcGranma = dynamic_cast<AliAODMCParticle*>(arrayMC->At(mother));
    if (mcGranma) {
      Int_t pdgGranma = mcGranma->GetPdgCode();
      Int_t abspdgGranma = TMath::Abs(pdgGranma);
      if (abspdgGranma == 4 || abspdgGranma == 5) {
        return kFALSE;
      }
      mother = mcGranma->GetMother();
    }
    else{
      AliError("Failed casting the mother particle!");
      break;
    }
  }

  return kTRUE;
}

Bool_t AliAnalysisTaskSEDStarCharmFraction::CalculateImpactParameter(AliAODTrack *track, Double_t &d0, Double_t &d0Err)
{ // Calculate impact parameter for a track
  AliExternalTrackParam etp;
  etp.CopyFromVTrack(track);
  d0 = 0; d0Err = 0;
  Double_t dz[2], covdz[3];
  if (etp.PropagateToDCA(fNewPrimVtx, fMagneticField, 3., dz, covdz)) {
    d0 = dz[0];
    d0Err = TMath::Sqrt(covdz[0]);
    return kTRUE;
  }

  return kFALSE;
}

Double_t AliAnalysisTaskSEDStarCharmFraction::CalculateImpactParameterDStar(AliAODRecoCascadeHF *cand)
{ // Calculate impact parameter of the D*
  // Calculation from AliAODRecoDecay::ImpParXY()

  Double_t p[3] = {cand->Px(), cand->Py(), cand->Pz()};
  Double_t ptsq = p[0]*p[0] + p[1]*p[1];

  Double_t prim[3] = {fNewPrimVtx->GetX(), fNewPrimVtx->GetY(), fNewPrimVtx->GetZ()};
  Double_t sec[3] = {fDStarVtx->GetX(), fDStarVtx->GetY(), fDStarVtx->GetZ()};

  Double_t k = - (sec[0]-prim[0])*p[0] - (sec[1]-prim[1])*p[1];
  k /= ptsq;
  Double_t dx = sec[0] - prim[0] + k*p[0];
  Double_t dy = sec[1] - prim[1] + k*p[1];
  Double_t absImpPar = TMath::Sqrt(dx*dx+dy*dy);

  TVector3 mom(p[0], p[1], p[2]);
  TVector3 fline(sec[0]-prim[0], sec[1]-prim[1], sec[2]-prim[2]);
  TVector3 cross = mom.Cross(fline);

  return (cross.Z() > 0. ? absImpPar : -absImpPar);
}

Double_t AliAnalysisTaskSEDStarCharmFraction::CalculateTrueImpactParameterDStar(AliAODMCHeader *headerMC, TClonesArray *arrayMC, AliAODRecoCascadeHF *cand)
{ // Calculate true impact parameter of the D*
  AliAODRecoDecayHF2Prong *candD0 = cand->Get2Prong();
  AliAODTrack *tr1 = (AliAODTrack*) candD0->GetDaughter(0);
  AliAODTrack *tr2 = (AliAODTrack*) candD0->GetDaughter(1);
  AliAODTrack *tr3 = (AliAODTrack*) cand->GetBachelor();

  Int_t label1 = tr1->GetLabel();
  Int_t label2 = tr2->GetLabel();
  Int_t label3 = tr3->GetLabel();

  if (label1 < 0 || label2 < 0 || label3 < 0) {
    return -9999.;
  }

  AliAODMCParticle *part1 = (AliAODMCParticle*) arrayMC->At(label1);
  AliAODMCParticle *part2 = (AliAODMCParticle*) arrayMC->At(label2);
  AliAODMCParticle *part3 = (AliAODMCParticle*) arrayMC->At(label3);

  if (part1 == 0 || part2 == 0 || part3 == 0) {
    return -9999.;
  }

  Double_t PxTotal = part1->Px()+part2->Px()+part3->Px();
  Double_t PyTotal = part1->Py()+part2->Py()+part3->Py();
  Double_t PzTotal = part1->Pz()+part2->Pz()+part3->Pz();


  // True primary+secondary vertices and momentum
  Double_t prim[3] = {headerMC->GetVtxX(), headerMC->GetVtxY(), headerMC->GetVtxZ()};
  Double_t sec[3] = {part3->Xv(), part3->Yv(), part3->Zv()}; // Put soft pion production vertex as secondary vertex
  Double_t p[3] = {PxTotal, PyTotal, PzTotal};


  Double_t ptsq = p[0]*p[0] + p[1]*p[1];
  Double_t k = - (sec[0]-prim[0])*p[0] - (sec[1]-prim[1])*p[1];
  k /= ptsq;
  Double_t dx = sec[0] - prim[0] + k*p[0];
  Double_t dy = sec[1] - prim[1] + k*p[1];
  Double_t absImpPar = TMath::Sqrt(dx*dx+dy*dy)*1.e4;

  TVector3 mom(p[0], p[1], p[2]);
  TVector3 fline(sec[0]-prim[0], sec[1]-prim[1], sec[2]-prim[2]);
  TVector3 cross = mom.Cross(fline);

  return cross.Z() > 0. ? absImpPar : -absImpPar;
}

void AliAnalysisTaskSEDStarCharmFraction::FillHistograms(AliAODRecoCascadeHF *cand)
{ // Fill histograms for a D* candidate
  Double_t Pt = cand->Pt();
  Int_t ptBin = fCuts->PtBin(Pt);
  Double_t deltaInvMass = cand->DeltaInvMass();
  Double_t impPar = CalculateImpactParameterDStar(cand)*1.e4;
  Double_t invMassD0 = cand->InvMassD0();

  FillHistogram("DeltaInvMass", deltaInvMass);
  if (ptBin >= 0) {
    FillHistogram(Form("DeltaInvMass_%d", ptBin), deltaInvMass);
  }

  FillHistogram("InvMassD0", invMassD0);
  if (ptBin >= 0) {
    FillHistogram(Form("InvMassD0_%d", ptBin), invMassD0);
  }

  FillRegionHistogram("ImpPar$r", impPar);
  if (ptBin >= 0) {
    FillRegionHistogram(Form("ImpPar$r_%d", ptBin), impPar);
  }

  // Store candidate in tree
  fPtForTree = Pt;
  fInvMassForTree = deltaInvMass;
  fImpParForTree = impPar;
  // fTrueImpParForTree was already updated before this
  fTypeForTree = fIsPeak ? 1 : (fIsSideband ? (deltaInvMass < fPDGMDStarD0 ? 2 : 3) : 0);
  fSignalTypeForTree = fReadMC ? (fIsSignalPrompt ? 1 : (fIsSignalFromB ? 2 : 0)) : -1;
  fTreeCandidate->Fill();

  AliAODTrack *softPi = cand->GetBachelor();
  Double_t d0, d0Err;
  if (CalculateImpactParameter(softPi, d0, d0Err)) {
    FillRegionHistogram("ImpParSoftPi$r", d0*1.e4);
    if (ptBin >= 0) {
      FillRegionHistogram(Form("ImpParSoftPi$r_%d", ptBin), d0*1.e4);
    }
  }

  AliAODRecoDecayHF2Prong *candD0 = cand->Get2Prong();
  AliAODTrack *pi; AliAODTrack *K; Int_t piId, KId;
  if (((AliAODTrack*) candD0->GetDaughter(0))->Charge()*cand->Charge() > 0) {
    // D* and first daughter charge are equal in sign <=> first daughter = pion
    piId = 0; KId = 1;
  }
  else {
    piId = 1; KId = 0;
  }

  pi = (AliAODTrack*) candD0->GetDaughter(piId);
  K = (AliAODTrack*) candD0->GetDaughter(KId);

  //pi (from D0) impact parameter
  if (CalculateImpactParameter(pi, d0, d0Err)) {
    FillRegionHistogram("ImpParD0Pi$r", d0*1.e4);
    if (ptBin >= 0) {
      FillRegionHistogram(Form("ImpParD0Pi$r_%d", ptBin), d0*1.e4);
    }
  }

  //pi (from D0) impact parameter without recalculation of primary vertex
  d0 = candD0->Getd0Prong(piId);
  FillRegionHistogram("ImpParD0Pi2$r", d0*1.e4);
  if (ptBin >= 0) {
    FillRegionHistogram(Form("ImpParD0Pi2$r_%d", ptBin), d0*1.e4);
  }

  //K (from D0) impact parameter
  if (CalculateImpactParameter(K, d0, d0Err)) {
    FillRegionHistogram("ImpParD0K$r", d0*1.e4);
    if (ptBin >= 0) {
      FillRegionHistogram(Form("ImpParD0K$r_%d", ptBin), d0*1.e4);
    }
  }

  //K (from D0) impact parameter without recalculation of primary vertex
  d0 = candD0->Getd0Prong(KId);
  FillRegionHistogram("ImpParD0K2$r", d0*1.e4);
  if (ptBin >= 0) {
    FillRegionHistogram(Form("ImpParD0K2$r_%d", ptBin), d0*1.e4);
  }

  AliAODVertex *origVtx = 0;
  if (candD0->GetOwnPrimaryVtx()) {
    origVtx = new AliAODVertex(*candD0->GetOwnPrimaryVtx());
  }
  candD0->SetOwnPrimaryVtx(fNewPrimVtx);

  Double_t impParD0 = candD0->ImpParXY()*1.e4;
  FillRegionHistogram("ImpParD0$r", impParD0);
  if (ptBin >= 0) {
    FillRegionHistogram(Form("ImpParD0$r_%d", ptBin), impParD0);
  }

  candD0->UnsetOwnPrimaryVtx();
  if (origVtx) {
    candD0->SetOwnPrimaryVtx(origVtx);
  }
  delete origVtx;
}

void AliAnalysisTaskSEDStarCharmFraction::FillHistogram(const char *name, Double_t value)
{ // Fill a specific histogram in multiple lists
  ((TH1D*) fListCandidate->FindObject(name))->Fill(value);

  if (fIsSignal) {
    ((TH1D*) fListSignal->FindObject(name))->Fill(value);
  }

  if (fIsSignalPrompt) {
    ((TH1D*) fListSignalPrompt->FindObject(name))->Fill(value);
  }

  if (fIsSignalFromB) {
    ((TH1D*) fListSignalFromB->FindObject(name))->Fill(value);
  }

  if (fIsBackground) {
    ((TH1D*) fListBackground->FindObject(name))->Fill(value);
  }
}

void AliAnalysisTaskSEDStarCharmFraction::FillRegionHistogram(const char *name, Double_t value)
{ // Fill a specific histogram in multiple lists, in the approprate regions (all, peak region, sideband region)
  FillHistogram(TString(name).ReplaceAll("$r", "").Data(), value);

  if (fIsPeak) {
    FillHistogram(TString(name).ReplaceAll("$r", "Peak").Data(), value);
  }
  
  if (fIsSideband) {
    FillHistogram(TString(name).ReplaceAll("$r", "Sideband").Data(), value);
  }
}

void AliAnalysisTaskSEDStarCharmFraction::FillTrueImpactParameter(AliAODRecoCascadeHF* cand)
{ // Fill histogram with true impact parameter distribution for D from B
  Double_t Pt = cand->Pt();
  Int_t ptBin = fCuts->PtBin(Pt);

  if (ptBin < 0) {
    return;
  }

  if (fTrueImpParForTree == -9999.) {
    return;
  }

  TH1D* hImpParTrue = (TH1D*) fListSignalFromB->FindObject(Form("ImpParTrue_%d", ptBin));
  hImpParTrue->Fill(fTrueImpParForTree);
}

void AliAnalysisTaskSEDStarCharmFraction::Terminate(const Option_t*)
{ // Terminate
}
