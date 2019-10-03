/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
#ifndef ALIPWG4HIGHPTTRACKQA_CXX
#define ALIPWG4HIGHPTTRACKQA_CXX

#include <iostream>
#include <memory>
#include <vector>
#include "AliPWG4HighPtTrackQA.h"

#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TProfile.h"
#include "TList.h"
#include "TFile.h"
#include "TChain.h"
#include "TH3F.h"
#include "TKey.h"
#include "TSystem.h"
#include "TBits.h"

#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliMCEvent.h"
#include "AliMCEventHandler.h"
#include "AliStack.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliExternalTrackParam.h"
#include "AliLog.h"
#include "AliGenPythiaEventHeader.h"
#include "AliGenCocktailEventHeader.h"
#include "AliCentrality.h"
#include "AliAODVertex.h"
#include "AliAODEvent.h"
#include "AliAnalysisUtils.h"

ClassImp(AliPWG4HighPtTrackQA);

/**
 * @brief Constructor
 */
AliPWG4HighPtTrackQA::AliPWG4HighPtTrackQA()
    : AliAnalysisTaskSE(),
      fDataType(kESD),
      fEvent(nullptr),
      fESD(nullptr),
      fVtx(nullptr),
      fVtxAOD(nullptr),
      fTrackCuts(nullptr),
      fTrackCutsITSLoose(nullptr),
      fTrackCutsTPConly(nullptr),
      fTrackType(0),
      fFilterMask(0),
      fIncludeNoITS(kFALSE),
      fSigmaConstrainedMax(-1.),
      fPtMax(100.),
      fIsPbPb(0),
      fCentClass(10),
      fEmcalTrigger(),
      fInit(0),
      fAliAnalysisUtils(nullptr),
      fNVertCont(0),
      fNVertSPDCont(0),
      fTklVsClusSPDCut(kFALSE),
      fZvertexDiff(0.5),
      fNVariables(27),
      fVariables(nullptr),
      fITSClusterMap(0),
      fAvgTrials(1),
      fNEventAll(0),
      fNEventSel(0),
      fNEventReject(0),
      fhEvMult(0),
      fhTrackletsMult(0),
      fh1Centrality(nullptr),
      fh1Xsec(0),
      fh1Trials(0),
      fh1PtHard(0),
      fh1PtHardTrials(0),
      fh1NTracksAll(nullptr),
      fh1NTracksReject(nullptr),
      fh1NTracksSel(nullptr),
      fPtAll(0),
      fPtSel(0),
      fPtPhi(nullptr),
      fPtEta(nullptr),
      fPtEtaPhi(nullptr),
      fPtDCA2D(nullptr),
      fPtDCAZ(nullptr),
      fPtNClustersTPC(nullptr),
      fPtNClustersTPCPhi(nullptr),
      fPtNClustersTPCIter1(nullptr),
      fPtNClustersTPCIter1Phi(nullptr),
      fPtNClustersTPCShared(nullptr),
      fPtNClustersTPCSharedFrac(nullptr),
      fPtNPointITS(nullptr),
      fPtNPointITSPhi(nullptr),
      fPtChi2C(nullptr),
      fPtNSigmaToVertex(nullptr),
      fPtRelUncertainty1Pt(nullptr),
      fPtRelUncertainty1PtNClus(nullptr),
      fPtRelUncertainty1PtNClusIter1(nullptr),
      fPtRelUncertainty1PtNPointITS(nullptr),
      fPtRelUncertainty1PtITSClusterMap(nullptr),
      fPtRelUncertainty1PtChi2(nullptr),
      fPtRelUncertainty1PtChi2Iter1(nullptr),
      fPtRelUncertainty1PtPhi(nullptr),
      fPtChi2PerClusterTPC(nullptr),
      fPtChi2PerClusterTPCIter1(nullptr),
      fPtNCrossedRows(nullptr),
      fPtNCrossedRowsPhi(nullptr),
      fPtNCrossedRowsNClusFPhi(nullptr),
      fPtNCrRNCrRNClusF(nullptr),
      fPtNCrossedRowsFit(nullptr),
      fPtNCrossedRowsFitPhi(nullptr),
      fPtNCrossedRowsNClusFFitPhi(nullptr),
      fNCrossedRowsNCrossedRowsFit(nullptr),
      fNClustersNCrossedRows(nullptr),
      fNClustersNCrossedRowsFit(nullptr),
      fPtNClustersNClustersFitMap(nullptr),
      fPtTPCSignalN(nullptr),
      fPtRelUncertainty1PtNCrossedRows(nullptr),
      fPtRelUncertainty1PtNCrossedRowsFit(nullptr),
      fPtChi2Gold(nullptr),
      fPtChi2GGC(nullptr),
      fPtChi2GoldPhi(nullptr),
      fPtChi2GGCPhi(nullptr),
      fChi2GoldChi2GGC(nullptr),
      fPtChi2ITSPhi(nullptr),
      fPtSigmaY2(nullptr),
      fPtSigmaZ2(nullptr),
      fPtSigmaSnp2(nullptr),
      fPtSigmaTgl2(nullptr),
      fPtSigma1Pt2(nullptr),
      fProfPtSigmaY2(nullptr),
      fProfPtSigmaZ2(nullptr),
      fProfPtSigmaSnp2(nullptr),
      fProfPtSigmaTgl2(nullptr),
      fProfPtSigma1Pt2(nullptr),
      fProfPtSigma1Pt(nullptr),
      fProfPtPtSigma1Pt(nullptr),
      fHistList(0)
{
  SetNVariables(27);

  fPtBinEdges[0][0] = 10.;
  fPtBinEdges[0][1] = 1.;
  fPtBinEdges[1][0] = 20.;
  fPtBinEdges[1][1] = 2.;
  fPtBinEdges[2][0] = 100.;
  fPtBinEdges[2][1] = 5.;
  fPtBinEdges[3][0] = 200.;
  fPtBinEdges[3][1] = 10.;

  memset(fVertex, 0, sizeof(double) *3);
  memset(fVertex, 0, sizeof(double) * 3);
}

/**
 * @brief Constructor
 * 
 * Initialization of Inputs and Outputs
 */
AliPWG4HighPtTrackQA::AliPWG4HighPtTrackQA(const char *name) : AliAnalysisTaskSE(name),
                                                               fDataType(kESD),
                                                               fEvent(nullptr),
                                                               fESD(nullptr),
                                                               fVtx(nullptr),
                                                               fVtxAOD(nullptr),
                                                               fTrackCuts(nullptr),
                                                               fTrackCutsITSLoose(nullptr),
                                                               fTrackCutsTPConly(nullptr),
                                                               fTrackType(0),
                                                               fFilterMask(0),
                                                               fIncludeNoITS(kFALSE),
                                                               fSigmaConstrainedMax(-1.),
                                                               fPtMax(100.),
                                                               fIsPbPb(0),
                                                               fCentClass(10),
                                                               fEmcalTrigger(),
                                                               fInit(0),
                                                               fAliAnalysisUtils(nullptr),
                                                               fNVertCont(0),
                                                               fNVertSPDCont(0),
                                                               fTklVsClusSPDCut(kFALSE),
                                                               fZvertexDiff(0.5),
                                                               fNVariables(27),
                                                               fVariables(nullptr),
                                                               fITSClusterMap(0),
                                                               fAvgTrials(1),
                                                               fNEventAll(0),
                                                               fNEventSel(0),
                                                               fNEventReject(0),
                                                               fhEvMult(0),
                                                               fhTrackletsMult(0),
                                                               fh1Centrality(nullptr),
                                                               fh1Xsec(0),
                                                               fh1Trials(0),
                                                               fh1PtHard(0),
                                                               fh1PtHardTrials(0),
                                                               fh1NTracksAll(nullptr),
                                                               fh1NTracksReject(nullptr),
                                                               fh1NTracksSel(nullptr),
                                                               fPtAll(0),
                                                               fPtSel(0),
                                                               fPtPhi(nullptr),
                                                               fPtEta(nullptr),
                                                               fPtEtaPhi(nullptr),
                                                               fPtDCA2D(nullptr),
                                                               fPtDCAZ(nullptr),
                                                               fPtNClustersTPC(nullptr),
                                                               fPtNClustersTPCPhi(nullptr),
                                                               fPtNClustersTPCIter1(nullptr),
                                                               fPtNClustersTPCIter1Phi(nullptr),
                                                               fPtNClustersTPCShared(nullptr),
                                                               fPtNClustersTPCSharedFrac(nullptr),
                                                               fPtNPointITS(nullptr),
                                                               fPtNPointITSPhi(nullptr),
                                                               fPtChi2C(nullptr),
                                                               fPtNSigmaToVertex(nullptr),
                                                               fPtRelUncertainty1Pt(nullptr),
                                                               fPtRelUncertainty1PtNClus(nullptr),
                                                               fPtRelUncertainty1PtNClusIter1(nullptr),
                                                               fPtRelUncertainty1PtNPointITS(nullptr),
                                                               fPtRelUncertainty1PtITSClusterMap(nullptr),
                                                               fPtRelUncertainty1PtChi2(nullptr),
                                                               fPtRelUncertainty1PtChi2Iter1(nullptr),
                                                               fPtRelUncertainty1PtPhi(nullptr),
                                                               fPtChi2PerClusterTPC(nullptr),
                                                               fPtChi2PerClusterTPCIter1(nullptr),
                                                               fPtNCrossedRows(nullptr),
                                                               fPtNCrossedRowsPhi(nullptr),
                                                               fPtNCrossedRowsNClusFPhi(nullptr),
                                                               fPtNCrRNCrRNClusF(nullptr),
                                                               fPtNCrossedRowsFit(nullptr),
                                                               fPtNCrossedRowsFitPhi(nullptr),
                                                               fPtNCrossedRowsNClusFFitPhi(nullptr),
                                                               fNCrossedRowsNCrossedRowsFit(nullptr),
                                                               fNClustersNCrossedRows(nullptr),
                                                               fNClustersNCrossedRowsFit(nullptr),
                                                               fPtNClustersNClustersFitMap(nullptr),
                                                               fPtTPCSignalN(nullptr),
                                                               fPtRelUncertainty1PtNCrossedRows(nullptr),
                                                               fPtRelUncertainty1PtNCrossedRowsFit(nullptr),
                                                               fPtChi2Gold(nullptr),
                                                               fPtChi2GGC(nullptr),
                                                               fPtChi2GoldPhi(nullptr),
                                                               fPtChi2GGCPhi(nullptr),
                                                               fChi2GoldChi2GGC(nullptr),
                                                               fPtChi2ITSPhi(nullptr),
                                                               fPtSigmaY2(nullptr),
                                                               fPtSigmaZ2(nullptr),
                                                               fPtSigmaSnp2(nullptr),
                                                               fPtSigmaTgl2(nullptr),
                                                               fPtSigma1Pt2(nullptr),
                                                               fProfPtSigmaY2(nullptr),
                                                               fProfPtSigmaZ2(nullptr),
                                                               fProfPtSigmaSnp2(nullptr),
                                                               fProfPtSigmaTgl2(nullptr),
                                                               fProfPtSigma1Pt2(nullptr),
                                                               fProfPtSigma1Pt(nullptr),
                                                               fProfPtPtSigma1Pt(nullptr),
                                                               fHistList(0)
{
  AliDebug(2, Form("AliPWG4HighPtTrackQA Calling Constructor"));

  SetNVariables(27);

  fPtBinEdges[0][0] = 10.;
  fPtBinEdges[0][1] = 1.;
  fPtBinEdges[1][0] = 20.;
  fPtBinEdges[1][1] = 2.;
  fPtBinEdges[2][0] = 100.;
  fPtBinEdges[2][1] = 5.;
  fPtBinEdges[3][0] = 200.;
  fPtBinEdges[3][1] = 10.;

  memset(fVertex, 0, sizeof(double) * 3);
  memset(fVertexSPD, 0, sizeof(double) * 3);

  // Input slot #0 works with a TChain ESD
  DefineInput(0, TChain::Class());
  // Output slot #1 write into a TList
  DefineOutput(1, TList::Class());
}

/**
 * @brief Set variable bin sizes for pT axis in histos
 * 
 * @param region 
 * @param ptmax min. pt
 * @param ptBinWidth Width of the bin
 */
void AliPWG4HighPtTrackQA::SetPtBinEdges(Int_t region, Double_t ptmax, Double_t ptBinWidth)
{
  if (region < 3)
  {
    fPtBinEdges[region][0] = ptmax;
    fPtBinEdges[region][1] = ptBinWidth;
  }
  else
  {
    AliError("Only 3 regions alowed. Use region 0/1/2\n");
    return;
  }
}

/**
 * @brief Create output objects
 */
void AliPWG4HighPtTrackQA::UserCreateOutputObjects()
{
  AliDebug(2, Form(">> AliPWG4HighPtTrackQA::UserCreateOutputObjects \n"));

  Bool_t oldStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  OpenFile(1);
  fHistList = new TList();
  fHistList->SetOwner(kTRUE);

  Float_t fgkPtMin = 0.;
  //  Float_t fgkPtMax = fPtMax;

  //fPtBinEdges[region][0] = ptmax of region ; fPtBinEdges[region][1] = binWidth of region
  const Float_t &ptmin1 = fgkPtMin;
  const Float_t &ptmax1 = fPtBinEdges[0][0];
  const Float_t &ptmin2 = ptmax1;
  const Float_t &ptmax2 = fPtBinEdges[1][0];
  const Float_t &ptmin3 = ptmax2;
  const Float_t &ptmax3 = fPtBinEdges[2][0]; //fgkPtMax;
  const Float_t &ptmin4 = ptmax3;
  const Float_t &ptmax4 = fPtBinEdges[3][0]; //fgkPtMax;
  const Int_t nbin11 = (int)((ptmax1 - ptmin1) / fPtBinEdges[0][1]);
  const Int_t nbin12 = (int)((ptmax2 - ptmin2) / fPtBinEdges[1][1]) + nbin11;
  const Int_t nbin13 = (int)((ptmax3 - ptmin3) / fPtBinEdges[2][1]) + nbin12;
  const Int_t nbin14 = (int)((ptmax4 - ptmin4) / fPtBinEdges[3][1]) + nbin13;
  Int_t fgkNPtBins = nbin14;
  //Create array with low edges of each bin
  std::vector<double> binsPt(fgkNPtBins + 1);
  for (Int_t i = 0; i <= fgkNPtBins; i++)
  {
    if (i <= nbin11)
      binsPt[i] = (Double_t)ptmin1 + (ptmax1 - ptmin1) / nbin11 * (Double_t)i;
    if (i <= nbin12 && i > nbin11)
      binsPt[i] = (Double_t)ptmin2 + (ptmax2 - ptmin2) / (nbin12 - nbin11) * ((Double_t)i - (Double_t)nbin11);
    if (i <= nbin13 && i > nbin12)
      binsPt[i] = (Double_t)ptmin3 + (ptmax3 - ptmin3) / (nbin13 - nbin12) * ((Double_t)i - (Double_t)nbin12);
    if (i <= nbin14 && i > nbin13)
      binsPt[i] = (Double_t)ptmin4 + (ptmax4 - ptmin4) / (nbin14 - nbin13) * ((Double_t)i - (Double_t)nbin13);
  }

  Int_t fgkNPhiBins = 18 * 6;
  Float_t kMinPhi = 0.;
  Float_t kMaxPhi = 2. * TMath::Pi();
  std::vector<double> binsPhi(fgkNPhiBins + 1);
  for (Int_t i = 0; i <= fgkNPhiBins; i++)
    binsPhi[i] = (Double_t)kMinPhi + (kMaxPhi - kMinPhi) / fgkNPhiBins * (Double_t)i;

  Int_t fgkNEtaBins = 20;
  Float_t fgkEtaMin = -1.;
  Float_t fgkEtaMax = 1.;
  std::vector<double> binsEta(fgkNEtaBins + 1);
  for (Int_t i = 0; i <= fgkNEtaBins; i++)
    binsEta[i] = (Double_t)fgkEtaMin + (fgkEtaMax - fgkEtaMin) / fgkNEtaBins * (Double_t)i;

  Int_t fgkNNClustersTPCBins = 80;
  Float_t fgkNClustersTPCMin = 0.5;
  Float_t fgkNClustersTPCMax = 160.5;
  std::vector<double> binsNClustersTPC(fgkNNClustersTPCBins + 1);
  for (Int_t i = 0; i <= fgkNNClustersTPCBins; i++)
    binsNClustersTPC[i] = (Double_t)fgkNClustersTPCMin + (fgkNClustersTPCMax - fgkNClustersTPCMin) / fgkNNClustersTPCBins * (Double_t)i;

  Int_t fgkNDCA2DBins = 80;
  Float_t fgkDCA2DMin = -0.2;
  Float_t fgkDCA2DMax = 0.2;
  if (fTrackType == 1 || fTrackType == 2 || fTrackType == 4 || fTrackType == 7)
  {
    fgkDCA2DMin = -2.;
    fgkDCA2DMax = 2.;
  }
  std::vector<double> binsDCA2D(fgkNDCA2DBins + 1);
  for (Int_t i = 0; i <= fgkNDCA2DBins; i++)
    binsDCA2D[i] = (Double_t)fgkDCA2DMin + (fgkDCA2DMax - fgkDCA2DMin) / fgkNDCA2DBins * (Double_t)i;

  Int_t fgkNDCAZBins = 80;
  Float_t fgkDCAZMin = -2.;
  Float_t fgkDCAZMax = 2.;
  if (fTrackType == 1 || fTrackType == 2 || fTrackType == 4)
  {
    fgkDCAZMin = -5.;
    fgkDCAZMax = 5.;
  }
  std::vector<double> binsDCAZ(fgkNDCAZBins + 1);
  for (Int_t i = 0; i <= fgkNDCAZBins; i++)
    binsDCAZ[i] = (Double_t)fgkDCAZMin + (fgkDCAZMax - fgkDCAZMin) / fgkNDCAZBins * (Double_t)i;

  Int_t fgkNNPointITSBins = 9;
  Float_t fgkNPointITSMin = -0.5;
  Float_t fgkNPointITSMax = 8.5;
  std::vector<double> binsNPointITS(fgkNNPointITSBins + 1);
  for (Int_t i = 0; i <= fgkNNPointITSBins; i++)
    binsNPointITS[i] = (Double_t)fgkNPointITSMin + (fgkNPointITSMax - fgkNPointITSMin) / fgkNNPointITSBins * (Double_t)i;

  Int_t fgkNITSClusterMapBins = 65;
  Float_t fgkITSClusterMapMin = -0.5;
  Float_t fgkITSClusterMapMax = 64.5;
  std::vector<double> binsITSClusterMap(fgkNITSClusterMapBins + 1);
  for (Int_t i = 0; i <= fgkNITSClusterMapBins; i++)
    binsITSClusterMap[i] = (Double_t)fgkITSClusterMapMin + (fgkITSClusterMapMax - fgkITSClusterMapMin) / fgkNITSClusterMapBins * (Double_t)i;

  Int_t fgkNNSigmaToVertexBins = 9;
  Float_t fgkNSigmaToVertexMin = 0.;
  Float_t fgkNSigmaToVertexMax = 9.;
  std::vector<double> binsNSigmaToVertex(fgkNNSigmaToVertexBins + 1);
  for (Int_t i = 0; i <= fgkNNSigmaToVertexBins; i++)
    binsNSigmaToVertex[i] = (Double_t)fgkNSigmaToVertexMin + (fgkNSigmaToVertexMax - fgkNSigmaToVertexMin) / fgkNNSigmaToVertexBins * (Double_t)i;

  Int_t fgkNChi2CBins = 10;
  //  Float_t fgkChi2CMin = 0.;
  //  Float_t fgkChi2CMax = 100.; //10 sigma
  std::vector<double> binsChi2C(fgkNChi2CBins + 1);
  for (Int_t i = 0; i <= fgkNChi2CBins; i++)
    binsChi2C[i] = (Double_t)i * (Double_t)i;

  Float_t fgkRel1PtUncertaintyMin = 0.;
  Float_t fgkRel1PtUncertaintyMax = 1.;
  Float_t binEdgeRel1PtUncertainty1 = 0.3;
  Int_t fgkNRel1PtUncertaintyBins1 = 45;
  Float_t binWidthRel1PtUncertainty1 = (binEdgeRel1PtUncertainty1 - fgkRel1PtUncertaintyMin) / ((Float_t)fgkNRel1PtUncertaintyBins1);
  Int_t fgkNRel1PtUncertaintyBins2 = 35;
  Float_t binWidthRel1PtUncertainty2 = (fgkRel1PtUncertaintyMax - binEdgeRel1PtUncertainty1) / ((Float_t)fgkNRel1PtUncertaintyBins2);
  Int_t fgkNRel1PtUncertaintyBins = fgkNRel1PtUncertaintyBins1 + fgkNRel1PtUncertaintyBins2;

  std::vector<double> binsRel1PtUncertainty(fgkNRel1PtUncertaintyBins + 1);
  for (Int_t i = 0; i <= fgkNRel1PtUncertaintyBins; i++)
  {
    if (i <= fgkNRel1PtUncertaintyBins1)
      binsRel1PtUncertainty[i] = (Double_t)fgkRel1PtUncertaintyMin + (Double_t)binWidthRel1PtUncertainty1 * (Double_t)i;
    if (i <= fgkNRel1PtUncertaintyBins && i > fgkNRel1PtUncertaintyBins1)
      binsRel1PtUncertainty[i] = (Double_t)binEdgeRel1PtUncertainty1 + (Double_t)binWidthRel1PtUncertainty2 * (Double_t)(i - fgkNRel1PtUncertaintyBins1);
  }

  Int_t fgkNUncertainty1PtBins = 30;
  Float_t fgkUncertainty1PtMin = 0.;
  Float_t fgkUncertainty1PtMax = 0.1;
  if (fTrackType == 1 || fTrackType == 2 || fTrackType == 4)
    fgkUncertainty1PtMax = 0.2;
  std::vector<double> binsUncertainty1Pt(fgkNUncertainty1PtBins + 1);
  for (Int_t i = 0; i <= fgkNUncertainty1PtBins; i++)
    binsUncertainty1Pt[i] = (Double_t)fgkUncertainty1PtMin + (fgkUncertainty1PtMax - fgkUncertainty1PtMin) / fgkNUncertainty1PtBins * (Double_t)i;

  Float_t fgkChi2PerClusMin = 0.;
  Float_t fgkChi2PerClusMax = 4.;
  Int_t fgkNChi2PerClusBins = (int)(fgkChi2PerClusMax * 10.);
  std::vector<double> binsChi2PerClus(fgkNChi2PerClusBins + 1);
  for (Int_t i = 0; i <= fgkNChi2PerClusBins; i++)
    binsChi2PerClus[i] = (Double_t)fgkChi2PerClusMin + (fgkChi2PerClusMax - fgkChi2PerClusMin) / fgkNChi2PerClusBins * (Double_t)i;

  Int_t fgkNCrossedRowsNClusFBins = 45;
  Float_t fgkNCrossedRowsNClusFMin = 0.;
  Float_t fgkNCrossedRowsNClusFMax = 1.5;
  std::vector<double> binsNCrossedRowsNClusF(fgkNCrossedRowsNClusFBins + 1);
  for (Int_t i = 0; i <= fgkNCrossedRowsNClusFBins; i++)
    binsNCrossedRowsNClusF[i] = (Double_t)fgkNCrossedRowsNClusFMin + (fgkNCrossedRowsNClusFMax - fgkNCrossedRowsNClusFMin) / fgkNCrossedRowsNClusFBins * (Double_t)i;

  Float_t fgk1PtMin = 0.;
  Float_t fgk1PtMax = 6.;
  Float_t binEdge1Pt1 = 1.;
  Float_t binWidth1Pt1 = 0.05;
  Int_t fgkN1PtBins1 = (int)((binEdge1Pt1 - fgk1PtMin) / binWidth1Pt1);
  Float_t binWidth1Pt2 = 0.1;
  Int_t fgkN1PtBins2 = (int)((fgk1PtMax - binEdge1Pt1) / binWidth1Pt2);
  Int_t fgkN1PtBins = fgkN1PtBins1 + fgkN1PtBins2;
  std::vector<double> bins1Pt(fgkN1PtBins + 1);

  for (Int_t i = 0; i <= fgkN1PtBins; i++)
  {
    if (i <= fgkN1PtBins1)
      bins1Pt[i] = (Double_t)fgk1PtMin + (Double_t)(binEdge1Pt1 - fgk1PtMin) / (Double_t)fgkN1PtBins1 * (Double_t)i;
    if (i <= fgkN1PtBins && i > fgkN1PtBins1)
      bins1Pt[i] = (Double_t)binEdge1Pt1 + (Double_t)(fgk1PtMax - binEdge1Pt1) / (Double_t)fgkN1PtBins2 * (Double_t)(i - fgkN1PtBins1);
  }

  Int_t fgkNSigmaY2Bins = 50;
  Float_t fgkSigmaY2Min = 0.;
  Float_t fgkSigmaY2Max = 1.;
  if (fTrackType == 1)
    fgkSigmaY2Max = 4.;
  if (fTrackType == 2 || fTrackType == 4)
    fgkSigmaY2Max = 0.1;
  std::vector<double> binsSigmaY2(fgkNSigmaY2Bins + 1);
  for (Int_t i = 0; i <= fgkNSigmaY2Bins; i++)
    binsSigmaY2[i] = (Double_t)fgkSigmaY2Min + (fgkSigmaY2Max - fgkSigmaY2Min) / fgkNSigmaY2Bins * (Double_t)i;

  Int_t fgkNSigmaZ2Bins = 50;
  Float_t fgkSigmaZ2Min = 0.;
  Float_t fgkSigmaZ2Max = 0.4;
  std::vector<double> binsSigmaZ2(fgkNSigmaZ2Bins + 1);
  for (Int_t i = 0; i <= fgkNSigmaZ2Bins; i++)
    binsSigmaZ2[i] = (Double_t)fgkSigmaZ2Min + (fgkSigmaZ2Max - fgkSigmaZ2Min) / fgkNSigmaZ2Bins * (Double_t)i;

  Int_t fgkNSigmaSnp2Bins = 50;
  Float_t fgkSigmaSnp2Min = 0.;
  Float_t fgkSigmaSnp2Max = 0.05;
  if (fTrackType == 1)
    fgkSigmaSnp2Max = 0.2;
  if (fTrackType == 2 || fTrackType == 4)
    fgkSigmaSnp2Max = 0.1;
  std::vector<double> binsSigmaSnp2(fgkNSigmaSnp2Bins + 1);
  for (Int_t i = 0; i <= fgkNSigmaSnp2Bins; i++)
    binsSigmaSnp2[i] = (Double_t)fgkSigmaSnp2Min + (fgkSigmaSnp2Max - fgkSigmaSnp2Min) / fgkNSigmaSnp2Bins * (Double_t)i;

  Int_t fgkNSigmaTgl2Bins = 50;
  Float_t fgkSigmaTgl2Min = 0.;
  Float_t fgkSigmaTgl2Max = 0.1;
  if (fTrackType == 1)
    fgkSigmaTgl2Max = 0.2;
  if (fTrackType == 2 || fTrackType == 4)
    fgkSigmaTgl2Max = 0.1;
  std::vector<double> binsSigmaTgl2(fgkNSigmaTgl2Bins + 1);
  for (Int_t i = 0; i <= fgkNSigmaTgl2Bins; i++)
    binsSigmaTgl2[i] = (Double_t)fgkSigmaTgl2Min + (fgkSigmaTgl2Max - fgkSigmaTgl2Min) / fgkNSigmaTgl2Bins * (Double_t)i;

  Int_t fgkNSigma1Pt2Bins = 50;
  Float_t fgkSigma1Pt2Min = 0.;
  Float_t fgkSigma1Pt2Max = 1.;
  std::vector<double> binsSigma1Pt2(fgkNSigma1Pt2Bins + 1);
  for (Int_t i = 0; i <= fgkNSigma1Pt2Bins; i++)
    binsSigma1Pt2[i] = (Double_t)fgkSigma1Pt2Min + (fgkSigma1Pt2Max - fgkSigma1Pt2Min) / fgkNSigma1Pt2Bins * (Double_t)i;

  Int_t fgkTPCsignalNBins = 80;
  Float_t fgkTPCsignalNMin = 0.5;
  Float_t fgkTPCsignalNMax = 160.5;
  std::vector<double> binsTPCsignalN(fgkTPCsignalNBins + 1);
  for (Int_t i = 0; i <= fgkTPCsignalNBins; i++)
    binsTPCsignalN[i] = (Double_t)fgkTPCsignalNMin + (fgkTPCsignalNMax - fgkTPCsignalNMin) / fgkTPCsignalNBins * (Double_t)i;

  fNEventAll = new TH1F("fNEventAll", "NEventAll", 1, -0.5, 0.5);
  fHistList->Add(fNEventAll);
  fNEventSel = new TH1F("fNEventSel", "NEvent Selected for analysis", 1, -0.5, 0.5);
  fHistList->Add(fNEventSel);
  fNEventReject = new TH1F("fNEventReject", "Reason events are rejectected for analysis", 20, 0, 20);
  //Set labels
  fNEventReject->Fill("noESD", 0);
  fNEventReject->Fill("Trigger", 0);
  fNEventReject->Fill("NTracks<2", 0);
  fNEventReject->Fill("noVTX", 0);
  fNEventReject->Fill("VtxStatus", 0);
  fNEventReject->Fill("NCont<2", 0);
  fNEventReject->Fill("ZVTX>10", 0);
  fNEventReject->Fill("PileupEvent", 0);
  fNEventReject->Fill("Bkg evt", 0);
  fNEventReject->Fill("VzSPD", 0);
  fNEventReject->Fill("cent", 0);
  fNEventReject->Fill("cent>90", 0);
  fHistList->Add(fNEventReject);

  fhTrackletsMult = new TH1F("fhTrackletsMult", "fhTrackletsMult;Tracklets Mult;#events)", 300, 0, 600);
  fHistList->Add(fhTrackletsMult);

  fhEvMult = new TH1F("fhEvMult", "fhEvMult;Ref Mult;#events)", 300, 0, 600);
  fHistList->Add(fhEvMult);

  fh1Centrality = new TH1F("fh1Centrality", "fh1Centrality; Centrality %", 100, 0, 100);
  fHistList->Add(fh1Centrality);

  fh1Xsec = new TProfile("fh1Xsec", "xsec from pyxsec.root", 1, 0, 1);
  fh1Xsec->GetXaxis()->SetBinLabel(1, "<#sigma>");
  fHistList->Add(fh1Xsec);

  fh1Trials = new TH1F("fh1Trials", "trials root file", 1, 0, 1);
  fh1Trials->GetXaxis()->SetBinLabel(1, "#sum{ntrials}");
  fHistList->Add(fh1Trials);

  fh1PtHard = new TH1F("fh1PtHard", "PYTHIA Pt hard;p_{T,hard}", 350, -.5, 349.5);
  fHistList->Add(fh1PtHard);
  fh1PtHardTrials = new TH1F("fh1PtHardTrials", "PYTHIA Pt hard weight with trials;p_{T,hard}", 350, -.5, 349.5);
  fHistList->Add(fh1PtHardTrials);

  fh1NTracksAll = new TH1F("fh1NTracksAll", "fh1NTracksAll; N of ESD tracks", 1, -0.5, 0.5);
  fHistList->Add(fh1NTracksAll);

  fh1NTracksReject = new TH1F("fh1NTracksReject", "fh1NTracksReject", 1, -0.5, 0.5);
  fh1NTracksReject->Fill("noHybridTrack", 0);
  fh1NTracksReject->Fill("noITSrefit", 0);
  fh1NTracksReject->Fill("noESDtrack", 0);
  fh1NTracksReject->Fill("noTPCInner", 0);
  fh1NTracksReject->Fill("FillTPC", 0);
  fh1NTracksReject->Fill("noTPConly", 0);
  fh1NTracksReject->Fill("relate", 0);
  fh1NTracksReject->Fill("trackCuts", 0);
  fh1NTracksReject->Fill("laser", 0);
  fh1NTracksReject->Fill("chi2", 0);
  fHistList->Add(fh1NTracksReject);

  fh1NTracksSel = new TH1F("fh1NTracksSel", "fh1NTracksSel;N of Selected ESD tracks", 1, -0.5, 0.5);
  fHistList->Add(fh1NTracksSel);

  fPtAll = new TH1F("fPtAll", "PtAll; #it{p}_{T, track}", fgkNPtBins, binsPt.data());
  fHistList->Add(fPtAll);
  fPtSel = new TH1F("fPtSel", "PtSel; #it{p}_{T, track}", fgkNPtBins, binsPt.data());
  fHistList->Add(fPtSel);

  fPtPhi = new TH2F("fPtPhi", "fPtPhi; #it{p}_{T, track}; #varphi", fgkNPtBins, binsPt.data(), fgkNPhiBins, binsPhi.data());
  fHistList->Add(fPtPhi);

  fPtEta = new TH2F("fPtEta", "fPtEta; #it{p}_{T, track}; #eta", fgkNPtBins, binsPt.data(), fgkNEtaBins, binsEta.data());
  fHistList->Add(fPtEta);

  fPtEtaPhi = new TH3F("fPtEtaPhi", "fPtEtaPhi; #it{p}_{T, track}; #eta; #varphi", fgkNPtBins, binsPt.data(), fgkNEtaBins, binsEta.data(), fgkNPhiBins, binsPhi.data());
  fHistList->Add(fPtEtaPhi);

  fPtDCA2D = new TH2F("fPtDCA2D", "fPtDCA2D; #it{p}_{T, track}; DCA_{xy}", fgkNPtBins, binsPt.data(), fgkNDCA2DBins, binsDCA2D.data());
  fHistList->Add(fPtDCA2D);

  fPtDCAZ = new TH2F("fPtDCAZ", "fPtDCAZ; #it{p}_{T, track}; DCA_{z}", fgkNPtBins, binsPt.data(), fgkNDCAZBins, binsDCAZ.data());
  fHistList->Add(fPtDCAZ);

  fPtNClustersTPC = new TH2F("fPtNClustersTPC", "fPtNClustersTPC; #it{p}_{T, track}; N Clusters TPC", fgkNPtBins, binsPt.data(), fgkNNClustersTPCBins, binsNClustersTPC.data());
  fHistList->Add(fPtNClustersTPC);

  fPtNClustersTPCPhi = new TH2F("fPtNClustersTPCPhi", "fPtNClustersTPCPhi; #it{p}_{T, track}; #varphi; N Clusters TPC (TPC only track)", fgkNPhiBins, binsPhi.data(), fgkNNClustersTPCBins, binsNClustersTPC.data());
  fHistList->Add(fPtNClustersTPCPhi);

  fPtNClustersTPCIter1 = new TH2F("fPtNClustersTPCIter1", "fPtNClustersTPCIter1; #it{p}_{T, track}; N Clusters TPC (TPC only track)", fgkNPtBins, binsPt.data(), fgkNNClustersTPCBins, binsNClustersTPC.data());
  fHistList->Add(fPtNClustersTPCIter1);

  fPtNClustersTPCIter1Phi = new TH3F("fPtNClustersTPCIter1Phi", "fPtNClustersTPCIter1Phi; #it{p}_{T, track}; N Clusters TPC (TPC only track); #varphi", fgkNPtBins, binsPt.data(), fgkNNClustersTPCBins, binsNClustersTPC.data(), fgkNPhiBins, binsPhi.data());
  fHistList->Add(fPtNClustersTPCIter1Phi);

  fPtNClustersTPCShared = new TH2F("fPtNClustersTPCShared", "fPtNClustersTPCShared; #it{p}_{T, track}; N Clusters TPC Shared", fgkNPtBins, binsPt.data(), fgkNNClustersTPCBins, binsNClustersTPC.data());
  fHistList->Add(fPtNClustersTPCShared);

  fPtNClustersTPCSharedFrac = new TH2F("fPtNClustersTPCSharedFrac", "fPtNClustersTPCSharedFrac; #it{p}_{T, track}; TPC Sharef Fraction", fgkNPtBins, binsPt.data(), fgkNSigma1Pt2Bins, binsSigma1Pt2.data());
  fHistList->Add(fPtNClustersTPCSharedFrac);

  fPtNPointITS = new TH2F("fPtNPointITS", "fPtNPointITS; #it{p}_{T, track}; N points ITS", fgkNPtBins, binsPt.data(), fgkNNPointITSBins, binsNPointITS.data());
  fHistList->Add(fPtNPointITS);

  fPtNPointITSPhi = new TH3F("fPtNPointITSPhi", "fPtNPointITSPhi; #it{p}_{T, track}; N points ITS; #varphi", fgkNPtBins, binsPt.data(), fgkNNPointITSBins, binsNPointITS.data(), fgkNPhiBins, binsPhi.data());
  fHistList->Add(fPtNPointITSPhi);

  fPtChi2C = new TH2F("fPtChi2C", "fPtChi2C; #it{p}_{T, track}; #chi^{2} Constrained TPC", fgkNPtBins, binsPt.data(), fgkNChi2CBins, binsChi2C.data());
  fHistList->Add(fPtChi2C);

  fPtNSigmaToVertex = new TH2F("fPtNSigmaToVertex", "fPtNSigmaToVertex; #it{p}_{T, track}; N #sigma to Vtx ", fgkNPtBins, binsPt.data(), fgkNNSigmaToVertexBins, binsNSigmaToVertex.data());
  fHistList->Add(fPtNSigmaToVertex);

  fPtRelUncertainty1Pt = new TH2F("fPtRelUncertainty1Pt", "fPtRelUncertainty1Pt; #it{p}_{T, track}; #sigma(#it{p}_{T, track})/#it{p}_{T, track} ", fgkNPtBins, binsPt.data(), fgkNRel1PtUncertaintyBins, binsRel1PtUncertainty.data());
  fHistList->Add(fPtRelUncertainty1Pt);

  fPtRelUncertainty1PtNClus = new TH3F("fPtRelUncertainty1PtNClus", "fPtRelUncertainty1PtNClus; #it{p}_{T, track}; #sigma(#it{p}_{T, track})/#it{p}_{T, track}", fgkNPtBins, binsPt.data(), fgkNRel1PtUncertaintyBins, binsRel1PtUncertainty.data(), fgkNNClustersTPCBins, binsNClustersTPC.data());
  fHistList->Add(fPtRelUncertainty1PtNClus);

  fPtRelUncertainty1PtNClusIter1 = new TH3F("fPtRelUncertainty1PtNClusIter1", "fPtRelUncertainty1PtNClusIter1; #it{p}_{T, track};  #sigma(#it{p}_{T, track})/#it{p}_{T, track}", fgkNPtBins, binsPt.data(), fgkNRel1PtUncertaintyBins, binsRel1PtUncertainty.data(), fgkNNClustersTPCBins, binsNClustersTPC.data());
  fHistList->Add(fPtRelUncertainty1PtNClusIter1);

  fPtRelUncertainty1PtNPointITS = new TH3F("fPtRelUncertainty1PtNPointITS", "fPtRelUncertainty1PtNPointITS; #it{p}_{T, track};  #sigma(#it{p}_{T, track})/#it{p}_{T, track}", fgkNPtBins, binsPt.data(), fgkNRel1PtUncertaintyBins, binsRel1PtUncertainty.data(), fgkNNPointITSBins, binsNPointITS.data());
  fHistList->Add(fPtRelUncertainty1PtNPointITS);

  fPtRelUncertainty1PtITSClusterMap = new TH3F("fPtRelUncertainty1PtITSClusterMap", "fPtRelUncertainty1PtITSClusterMap; #it{p}_{T, track}; #sigma(#it{p}_{T, track})/#it{p}_{T, track}", fgkNPtBins, binsPt.data(), fgkNRel1PtUncertaintyBins, binsRel1PtUncertainty.data(), fgkNITSClusterMapBins, binsITSClusterMap.data());
  fHistList->Add(fPtRelUncertainty1PtITSClusterMap);

  fPtRelUncertainty1PtChi2 = new TH3F("fPtRelUncertainty1PtChi2", "fPtRelUncertainty1PtChi2; #it{p}_{T, track}; #sigma(#it{p}_{T, track})/#it{p}_{T, track}", fgkNPtBins, binsPt.data(), fgkNRel1PtUncertaintyBins, binsRel1PtUncertainty.data(), fgkNChi2PerClusBins, binsChi2PerClus.data());
  fHistList->Add(fPtRelUncertainty1PtChi2);

  fPtRelUncertainty1PtChi2Iter1 = new TH3F("fPtRelUncertainty1PtChi2Iter1", "fPtRelUncertainty1PtChi2Iter1; #it{p}_{T, track}; #sigma(#it{p}_{T, track})/#it{p}_{T, track}", fgkNPtBins, binsPt.data(), fgkNRel1PtUncertaintyBins, binsRel1PtUncertainty.data(), fgkNChi2PerClusBins, binsChi2PerClus.data());
  
  fHistList->Add(fPtRelUncertainty1PtChi2Iter1);

  fPtRelUncertainty1PtPhi = new TH3F("fPtRelUncertainty1PtPhi", "fPtRelUncertainty1PtPhi; #it{p}_{T, track};  #sigma(#it{p}_{T, track})/#it{p}_{T, track}", fgkNPtBins, binsPt.data(), fgkNRel1PtUncertaintyBins, binsRel1PtUncertainty.data(), fgkNPhiBins, binsPhi.data());
  fHistList->Add(fPtRelUncertainty1PtPhi);

  fPtChi2PerClusterTPC = new TH2F("fPtChi2PerClusterTPC", "fPtChi2PerClusterTPC; #it{p}_{T, track}; ", fgkNPtBins, binsPt.data(), fgkNChi2PerClusBins, binsChi2PerClus.data());
  fHistList->Add(fPtChi2PerClusterTPC);

  fPtChi2PerClusterTPCIter1 = new TH2F("fPtChi2PerClusterTPCIter1", "fPtChi2PerClusterTPCIter1; #it{p}_{T, track}; #chi^{2}", fgkNPtBins, binsPt.data(), fgkNChi2PerClusBins, binsChi2PerClus.data());
  fHistList->Add(fPtChi2PerClusterTPCIter1);

  fPtNCrossedRows = new TH2F("fPtNCrossedRows", "fPtNCrossedRows; #it{p}_{T, track}; N crossed rows", fgkNPtBins, binsPt.data(), fgkNNClustersTPCBins, binsNClustersTPC.data());
  fHistList->Add(fPtNCrossedRows);

  fPtNCrossedRowsPhi = new TH3F("fPtNCrossedRowsPhi", "fPtNCrossedRowsPhi; #it{p}_{T, track}; N crossed rows ; #varphi ", fgkNPtBins, binsPt.data(), fgkNNClustersTPCBins, binsNClustersTPC.data(), fgkNPhiBins, binsPhi.data());
  fHistList->Add(fPtNCrossedRowsPhi);

  fPtNCrossedRowsNClusFPhi = new TH3F("fPtNCrossedRowsNClusFPhi", "fPtNCrossedRowsNClusFPhi; #it{p}_{T, track}; N crossed raw / N findable clusters; #varphi", fgkNPtBins, binsPt.data(), fgkNCrossedRowsNClusFBins, binsNCrossedRowsNClusF.data(), fgkNPhiBins, binsPhi.data());
  fHistList->Add(fPtNCrossedRowsNClusFPhi);

  fPtNCrRNCrRNClusF = new TH3F("fPtNCrRNCrRNClusF", "fPtNCrRNCrRNClusF; #it{p}_{T, track}; N Clusters TPC; N crossed raw / N findable clusters; ", fgkNPtBins, binsPt.data(), fgkNNClustersTPCBins, binsNClustersTPC.data(), fgkNCrossedRowsNClusFBins, binsNCrossedRowsNClusF.data());
  fHistList->Add(fPtNCrRNCrRNClusF);

  fPtNCrossedRowsFit = new TH2F("fPtNCrossedRowsFit", "fPtNCrossedRowsFit; #it{p}_{T, track}; N crossed rows from fit map", fgkNPtBins, binsPt.data(), fgkNNClustersTPCBins, binsNClustersTPC.data());
  fHistList->Add(fPtNCrossedRowsFit);

  fPtNCrossedRowsFitPhi = new TH3F("fPtNCrossedRowsFitPhi", "fPtNCrossedRowsFitPhi; #it{p}_{T, track}; N crossed rows from fit map; #varphi ", fgkNPtBins, binsPt.data(), fgkNNClustersTPCBins, binsNClustersTPC.data(), fgkNPhiBins, binsPhi.data());
  fHistList->Add(fPtNCrossedRowsFitPhi);

  fPtNCrossedRowsNClusFFitPhi = new TH3F("fPtNCrossedRowsNClusFFitPhi", "fPtNCrossedRowsNClusFFitPhi; #it{p}_{T, track}; ", fgkNPtBins, binsPt.data(), fgkNCrossedRowsNClusFBins, binsNCrossedRowsNClusF.data(), fgkNPhiBins, binsPhi.data());
  fHistList->Add(fPtNCrossedRowsNClusFFitPhi);

  fNCrossedRowsNCrossedRowsFit = new TH2F("fNCrossedRowsNCrossedRowsFit", "fNCrossedRowsNCrossedRowsFit; N crossed rows; N crossed rows from fit map", fgkNNClustersTPCBins, binsNClustersTPC.data(), fgkNNClustersTPCBins, binsNClustersTPC.data());
  fHistList->Add(fNCrossedRowsNCrossedRowsFit);

  fNClustersNCrossedRows = new TH2F("fNClustersNCrossedRows", "fNClustersNCrossedRows; N clusters; N crossed rows", fgkNNClustersTPCBins, binsNClustersTPC.data(), fgkNNClustersTPCBins, binsNClustersTPC.data());
  fHistList->Add(fNClustersNCrossedRows);

  fNClustersNCrossedRowsFit = new TH2F("fNClustersNCrossedRowsFit", "fNClustersNCrossedRowsFit; N clusters; N crossed rows from fit ", fgkNNClustersTPCBins, binsNClustersTPC.data(), fgkNNClustersTPCBins, binsNClustersTPC.data());
  fHistList->Add(fNClustersNCrossedRowsFit);

  fPtNClustersNClustersFitMap = new TH3F("fPtNClustersNClustersFitMap", "fPtNClustersNClustersFitMap; #it{p}_{T, track}; N_{cls};N_{cls}^{fit map}", fgkNPtBins, binsPt.data(), fgkNNClustersTPCBins, binsNClustersTPC.data(), fgkNNClustersTPCBins, binsNClustersTPC.data());
  fHistList->Add(fPtNClustersNClustersFitMap);

  fPtRelUncertainty1PtNCrossedRows = new TH3F("fPtRelUncertainty1PtNCrossedRows", "fPtRelUncertainty1PtNCrossedRows; #it{p}_{T, track}; N crossed rows", fgkNPtBins, binsPt.data(), fgkNRel1PtUncertaintyBins, binsRel1PtUncertainty.data(), fgkNNClustersTPCBins, binsNClustersTPC.data());
  fHistList->Add(fPtRelUncertainty1PtNCrossedRows);

  fPtRelUncertainty1PtNCrossedRowsFit = new TH3F("fPtRelUncertainty1PtNCrossedRowsFit", "fPtRelUncertainty1PtNCrossedRowsFit; #it{p}_{T, track}; #sigma(1./#it{p}_{T, track}); N crossed rows from fit", fgkNPtBins, binsPt.data(), fgkNRel1PtUncertaintyBins, binsRel1PtUncertainty.data(), fgkNNClustersTPCBins, binsNClustersTPC.data());
  fHistList->Add(fPtRelUncertainty1PtNCrossedRowsFit);

  fPtChi2Gold = new TH2F("fPtChi2Gold", "fPtChi2Gold; #it{p}_{T, track}; #chi^{2} gold ", fgkNPtBins, binsPt.data(), fgkNChi2CBins, binsChi2C.data());
  fHistList->Add(fPtChi2Gold);

  fPtChi2GGC = new TH2F("fPtChi2GGC", "fPtChi2GGC; #it{p}_{T, track}; #chi^{2}_{ggc} ", fgkNPtBins, binsPt.data(), fgkNChi2CBins, binsChi2C.data());
  fHistList->Add(fPtChi2GGC);

  fPtChi2GoldPhi = new TH3F("fPtChi2GoldPhi", "fPtChi2GoldPhi; #it{p}_{T, track}; #chi^{2} gold ; #varphi", fgkNPtBins, binsPt.data(), fgkNChi2CBins, binsChi2C.data(), fgkNPhiBins, binsPhi.data());
  fHistList->Add(fPtChi2GoldPhi);

  fPtChi2GGCPhi = new TH3F("fPtChi2GGCPhi", "fPtChi2GGCPhi; #it{p}_{T, track};#chi^{2}_{ggc}; #varphi ", fgkNPtBins, binsPt.data(), fgkNChi2CBins, binsChi2C.data(), fgkNPhiBins, binsPhi.data());
  fHistList->Add(fPtChi2GGCPhi);

  fChi2GoldChi2GGC = new TH2F("fChi2GoldChi2GGC", "fChi2GoldChi2GGC;#chi^{2}_{gold};#chi^{2}_{ggc}", fgkNChi2CBins, binsChi2C.data(), fgkNChi2CBins, binsChi2C.data());
  fHistList->Add(fChi2GoldChi2GGC);

  fPtChi2ITSPhi = new TH3F("fPtChi2ITSPhi", "fPtChi2ITSPhi; #it{p}_{T, track}; #chi^{2}_{ITS};#varphi", fgkNPtBins, binsPt.data(), fgkNChi2CBins, binsChi2C.data(), fgkNPhiBins, binsPhi.data());
  fHistList->Add(fPtChi2ITSPhi);

  fPtTPCSignalN = new TH2F("fPtTPCSignalN", "fPtTPCSignalN ; #it{p}_{T, track}; TPCsignalN", fgkNPtBins, binsPt.data(), fgkTPCsignalNBins, binsTPCsignalN.data());
  fHistList->Add(fPtTPCSignalN);

  fPtSigmaY2 = new TH2F("fPtSigmaY2", "fPtSigmaY2; 1./#it{p}_{T, track}; #sigma y2 ", fgkN1PtBins, bins1Pt.data(), fgkNSigmaY2Bins, binsSigmaY2.data());
  fHistList->Add(fPtSigmaY2);

  fPtSigmaZ2 = new TH2F("fPtSigmaZ2", "fPtSigmaZ2; #it{p}_{T, track}; #sigma z2 ", fgkN1PtBins, bins1Pt.data(), fgkNSigmaZ2Bins, binsSigmaZ2.data());
  fHistList->Add(fPtSigmaZ2);

  fPtSigmaSnp2 = new TH2F("fPtSigmaSnp2", "fPtSigmaSnp2; 1./#it{p}_{T, track};  #sigma Snp2", fgkN1PtBins, bins1Pt.data(), fgkNSigmaSnp2Bins, binsSigmaSnp2.data());
  fHistList->Add(fPtSigmaSnp2);

  fPtSigmaTgl2 = new TH2F("fPtSigmaTgl2", "fPtSigmaTgl2; 1./#it{p}_{T, track}; #sigma Tgl2 ", fgkN1PtBins, bins1Pt.data(), fgkNSigmaTgl2Bins, binsSigmaTgl2.data());
  fHistList->Add(fPtSigmaTgl2);

  fPtSigma1Pt2 = new TH2F("fPtSigma1Pt2", "fPtSigma1Pt2; 1./#it{p}_{T, track}; #sigma(1./#it{p}_{T, track}); ", fgkN1PtBins, bins1Pt.data(), fgkNSigma1Pt2Bins, binsSigma1Pt2.data());
  fHistList->Add(fPtSigma1Pt2);

  fProfPtSigmaY2 = new TProfile("fProfPtSigmaY2", "fProfPtSigmaY2; 1./#it{p}_{T, track}; #sigma y2", fgkN1PtBins, bins1Pt.data());
  fHistList->Add(fProfPtSigmaY2);

  fProfPtSigmaZ2 = new TProfile("fProfPtSigmaZ2", "fProfPtSigmaZ2; 1./#it{p}_{T, track}; #sigma z2", fgkN1PtBins, bins1Pt.data());
  fHistList->Add(fProfPtSigmaZ2);

  fProfPtSigmaSnp2 = new TProfile("fProfPtSigmaSnp2", "fProfPtSigmaSnp2; 1./#it{p}_{T, track}; #sigma Snp2", fgkN1PtBins, bins1Pt.data());
  fHistList->Add(fProfPtSigmaSnp2);

  fProfPtSigmaTgl2 = new TProfile("fProfPtSigmaTgl2", "fProfPtSigmaTgl2; 1./#it{p}_{T, track}; #sigma Tgl2 ", fgkN1PtBins, bins1Pt.data());
  fHistList->Add(fProfPtSigmaTgl2);

  fProfPtSigma1Pt2 = new TProfile("fProfPtSigma1Pt2", "fProfPtSigma1Pt2", fgkN1PtBins, bins1Pt.data());
  fHistList->Add(fProfPtSigma1Pt2);

  fProfPtSigma1Pt = new TProfile("fProfPtSigma1Pt", "fProfPtSigma1Pt;p_{T};#sigma(1/p_{T}); #sigma(1./#it{p}_{T, track});", fgkNPtBins, binsPt.data());
  fHistList->Add(fProfPtSigma1Pt);

  fProfPtPtSigma1Pt = new TProfile("fProfPtPtSigma1Pt", "fProfPtPtSigma1Pt;p_{T};p_{T}#sigma(1/p_{T})", fgkNPtBins, binsPt.data());
  fHistList->Add(fProfPtPtSigma1Pt);

  TH1::AddDirectory(oldStatus);

  PostData(1, fHistList);
}

/**
 * @brief Decide if event should be selected for analysis
 * 
 * Checks following requirements:
 *  - fEvent available
 *  - trigger info from AliPhysicsSelection
 *  - MCevent available
 *  - number of reconstructed tracks > 1
 *  - primary vertex reconstructed
 *  - z-vertex < 10 cm
 *  - centrality in case of PbPb
 *
 * @return Bool_t True if the event is selected, false otherwise
 */
Bool_t AliPWG4HighPtTrackQA::SelectEvent()
{

  Bool_t selectEvent = kTRUE;

  //fEvent object available?
  if (!fEvent)
  {
    AliDebug(2, Form("ERROR: fInputEvent not available\n"));
    fNEventReject->Fill("noAliVEvent", 1);
    selectEvent = kFALSE;
    return selectEvent;
  }

  if(fEmcalTrigger.Length() && !fEvent->GetFiredTriggerClasses().Contains(fEmcalTrigger)) {
    AliDebug(2, Form("ERROR: Event not having the correct trigger string\n"));
    fNEventReject->Fill("Trigger", 1);
    selectEvent = kFALSE;
    return selectEvent;
  }

  //Check if number of reconstructed tracks is larger than 1
  if (!fEvent->GetNumberOfTracks() || fEvent->GetNumberOfTracks() < 2)
  {
    fNEventReject->Fill("NTracks<2", 1);
    selectEvent = kFALSE;
    return selectEvent;
  }

  //Check if vertex is reconstructed
  if (fDataType == kESD && dynamic_cast<AliESDEvent *>(fEvent))
  {
    fVtx = ((AliESDEvent *)fEvent)->GetPrimaryVertexTracks();

    if (!fVtx || !fVtx->GetStatus())
      fVtx = ((AliESDEvent *)fEvent)->GetPrimaryVertexSPD();

    if (!fVtx)
    {
      fNEventReject->Fill("noVTX", 1);
      selectEvent = kFALSE;
      return selectEvent;
    }

    if (!fVtx->GetStatus())
    {
      fNEventReject->Fill("VtxStatus", 1);
      selectEvent = kFALSE;
      return selectEvent;
    }

    // Need vertex cut
    if (fVtx->GetNContributors() < 2)
    {
      fNEventReject->Fill("NCont<2", 1);
      selectEvent = kFALSE;
      return selectEvent;
    }

    //Check if z-vertex < 10 cm
    double primVtx[3];
    fVtx->GetXYZ(primVtx);
    if (TMath::Sqrt(primVtx[0] * primVtx[0] + primVtx[1] * primVtx[1]) > 1. || TMath::Abs(primVtx[2] > 10.))
    {
      fNEventReject->Fill("ZVTX>10", 1);
      selectEvent = kFALSE;
      return selectEvent;
    }
  }
  else if (fDataType == kAOD && dynamic_cast<AliAODEvent *>(fEvent))
  {
    fVtxAOD = ((AliAODEvent *)fEvent)->GetPrimaryVertex();
    if (!fVtxAOD)
    {
      fNEventReject->Fill("noVTX", 1);
      selectEvent = kFALSE;
      return selectEvent;
    }

    // Need vertex cut
    if (fVtxAOD->GetNContributors() < 2)
    {
      fNEventReject->Fill("NCont<2", 1);
      selectEvent = kFALSE;
      return selectEvent;
    }

    //Check if z-vertex < 10 cm
    double primVtx[3];
    fVtxAOD->GetXYZ(primVtx);
    if (TMath::Sqrt(primVtx[0] * primVtx[0] + primVtx[1] * primVtx[1]) > 1. || TMath::Abs(primVtx[2] > 10.))
    {
      fNEventReject->Fill("ZVTX>10", 1);
      selectEvent = kFALSE;
      return selectEvent;
    }

    fAliAnalysisUtils = new AliAnalysisUtils();
    fAliAnalysisUtils->SetMinVtxContr(2);

    if (fAliAnalysisUtils->IsPileUpEvent(InputEvent()))
    {
      fNEventReject->Fill("PileupEvent", 1);
      return selectEvent;
    }

    if (fTklVsClusSPDCut && fAliAnalysisUtils->IsSPDClusterVsTrackletBG(InputEvent()))
    {
      fNEventReject->Fill("Bkg evt", 1);
      return selectEvent;
    }

    const AliVVertex *vertSPD = InputEvent()->GetPrimaryVertexSPD();
    if (vertSPD)
    {
      vertSPD->GetXYZ(fVertexSPD);
      fNVertSPDCont = vertSPD->GetNContributors();
    }

    if (fNVertSPDCont > 0 && fZvertexDiff < 999)
    {
      Double_t vzSPD = fVertexSPD[2];
      Double_t dvertex = TMath::Abs(primVtx[2] - vzSPD);
      //if difference larger than fZvertexDiff
      if (dvertex > fZvertexDiff)
      {
        fNEventReject->Fill("VzSPD", 1);
        return selectEvent;
      }
    }

    AliVMultiplicity *multiplicity = (AliVMultiplicity *)InputEvent()->GetMultiplicity();
    if (!multiplicity)
      return selectEvent;

    Int_t trackletsMult = multiplicity->GetNumberOfTracklets();
    fhTrackletsMult->Fill(trackletsMult);

    Double_t refMult = 0.;
    if (fDataType == kAOD)
    {
      if (((AliVAODHeader *)dynamic_cast<AliAODEvent *>(fEvent)->GetHeader())->GetRefMultiplicity())
        refMult = ((AliVAODHeader *)dynamic_cast<AliAODEvent *>(fEvent)->GetHeader())->GetRefMultiplicity();
    }
    else
      (refMult = 0.);
    fhEvMult->Fill(refMult);
  }
  //Centrality selection should only be done in case of PbPb
  if (IsPbPb())
  {
    Float_t cent = 0.;
    if (fCentClass != CalculateCentrality(fEvent) && fCentClass != 10)
    {
      fNEventReject->Fill("cent", 1);
      selectEvent = kFALSE;
      return selectEvent;
    }
    else
    {
      if (fDataType == kESD)
      {
        if (dynamic_cast<AliESDEvent *>(fEvent)->GetCentrality())
        {
          cent = dynamic_cast<AliESDEvent *>(fEvent)->GetCentrality()->GetCentralityPercentile("V0M");
        }
      }
      else if (fDataType == kAOD)
      {
        if (((AliVAODHeader *)dynamic_cast<AliAODEvent *>(fEvent)->GetHeader())->GetCentrality())
          cent = ((AliVAODHeader *)dynamic_cast<AliAODEvent *>(fEvent)->GetHeader())->GetCentrality();
      }
      if (cent > 90.)
      {
        fNEventReject->Fill("cent>90", 1);
        selectEvent = kFALSE;
        return selectEvent;
      }
      fh1Centrality->Fill(cent);
    }
  }

  return selectEvent;
}

/**
 * @brief Get centrality from ESD or AOD
 * 
 * @param ev Event to check
 * @return Int_t Centrality class
 */
Int_t AliPWG4HighPtTrackQA::CalculateCentrality(AliVEvent *ev)
{
  if (fDataType == kESD)
    return CalculateCentrality(dynamic_cast<AliESDEvent *>(ev));
  else if (fDataType == kAOD)
    return CalculateCentrality(dynamic_cast<AliAODEvent *>(ev));
  else
    return 5;
}

/**
 * @brief Get centrality from ESD
 */
Int_t AliPWG4HighPtTrackQA::CalculateCentrality(AliESDEvent *esd)
{
  Float_t cent = -1;

  if (esd)
  {
    if (esd->GetCentrality())
    {
      cent = esd->GetCentrality()->GetCentralityPercentile("V0M");
      if (fDebug > 3)
        printf("centrality: %f\n", cent);
    }
  }

  return GetCentralityClass(cent);
}

/**
 * @brief Get centrality from AOD
 */
Int_t AliPWG4HighPtTrackQA::CalculateCentrality(const AliAODEvent *aod)
{
  if (!aod)
    return 5;
  Float_t cent = ((AliVAODHeader *)aod->GetHeader())->GetCentrality();
  if (fDebug > 3)
    printf("centrality: %f\n", cent);

  return GetCentralityClass(cent);
}

/**
 * @brief Get centrality class
 */
Int_t AliPWG4HighPtTrackQA::GetCentralityClass(Float_t cent) const
{
  if (cent < 0)
    return 5; // OB - cent sometimes negative
  if (cent > 80)
    return 4;
  if (cent > 50)
    return 3;
  if (cent > 30)
    return 2;
  if (cent > 10)
    return 1;
  return 0;
}

void AliPWG4HighPtTrackQA::Init()
{
  if (!fInit && fDataType == kESD)
  {
    Printf("Init magnetic field ---------");
    fESD = dynamic_cast<AliESDEvent *>(InputEvent());
    fESD->InitMagneticField();
    fInit = kTRUE;
  }
}

/**
 * @brief Main loop
 * 
 * Called for each event
 * 
 * @param Option_t Not used
 */
void AliPWG4HighPtTrackQA::UserExec(Option_t *)
{
  AliDebug(2, Form(">> AliPWG4HighPtTrackQA::UserExec \n"));
  Init();
  fEvent = InputEvent();
  fESD = dynamic_cast<AliESDEvent *>(InputEvent());

  // All events without selection
  fNEventAll->Fill(0.);

  if (!SelectEvent())
  {
    // Post output data
    PostData(1, fHistList);
    return;
  }

  //Need to keep track of selected events
  fNEventSel->Fill(0.);

  fVariables = new TArrayF(fNVariables);

  if (fDataType == kESD)
    DoAnalysisESD();
  if (fDataType == kAOD)
    DoAnalysisAOD();

  //Delete old fVariables
  if (fVariables)
    delete fVariables;

  // Post output data
  PostData(1, fHistList);
}

/**
 * @brief Run analysis on ESD
 */
void AliPWG4HighPtTrackQA::DoAnalysisESD()
{

  if (!fESD)
  {
    PostData(1, fHistList);
    return;
  }

  // ---- Get MC Header information (for MC productions in pThard bins) ----
  Double_t ptHard = 0.;
  Double_t nTrials = 1; // trials for MC trigger weight for real data

  AliMCEventHandler *eventHandlerMC = dynamic_cast<AliMCEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (eventHandlerMC)
  {

    if (eventHandlerMC->MCEvent())
    {
      AliGenPythiaEventHeader *pythiaGenHeader = GetPythiaEventHeader(eventHandlerMC->MCEvent());
      if (pythiaGenHeader)
      {
        nTrials = pythiaGenHeader->Trials();
        ptHard = pythiaGenHeader->GetPtHard();

        fh1PtHard->Fill(ptHard);
        fh1PtHardTrials->Fill(ptHard, nTrials);

        fh1Trials->Fill("#sum{ntrials}", fAvgTrials);
      }
    }
  }

  Int_t nTracks = fESD->GetNumberOfTracks();
  AliDebug(2, Form("nTracks ESD%d", nTracks));

  /*
    Variables to be put in fVariables
    0: pt
    1: phi
    2: eta
    3: dca2D
    4: dcaZ 
    5: nClustersTPC
    6: nPointITS   
    7: chi2C       
    8: nSigmaToVertex
    9: trackLengthTPC
    10: chi2PerClusterTPC
    11: #crossed rows
    12: (#crossed rows)/(#findable clusters)
    13: SigmaY2
    14: SigmaZ2
    15: SigmaSnp2
    16: SigmaTgl2
    17: Sigma1Pt2
    18: NClustersTPCIter1 (Iter1: first iteration of the tracking, from the outer wall of the TPC till the inner wall is reached. TPCNclsIter1 it the number of clusters of a TPConly track without the use of ITS information)
    19: Chi2TPCIter1
    20: nClustersTPCShared
    21: Golden Chi2 - global vs TPC constrained
    22: Chi2 between global and global constrained
    23: #crossed rows from fit map
    24: (#crossed rows)/(#findable clusters) from fit map
    25: chi2ITS
    26: TPCsignalN
  */

  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++)
  {
    fh1NTracksAll->Fill(0.);

    //Get track for analysis
    AliESDtrack *track = nullptr;
    AliESDtrack *esdtrack = fESD->GetTrack(iTrack);
    if (!esdtrack)
    {
      fh1NTracksReject->Fill("noESDtrack", 1);
      continue;
    }
    std::unique_ptr<AliESDtrack> origtrack(new AliESDtrack(*esdtrack));
    if (!origtrack)
      continue;

    if (fTrackType == 4)
    {
      if (!(fTrackCuts->AcceptTrack(esdtrack)))
      {
        fh1NTracksReject->Fill("trackCuts", 1);
        continue;
      }
    }

    if (fTrackType == 1)
      track = AliESDtrackCuts::GetTPCOnlyTrack(fESD, esdtrack->GetID());
    else if (fTrackType == 2 || fTrackType == 4)
    {
      track = AliESDtrackCuts::GetTPCOnlyTrack(const_cast<AliESDEvent *>(fESD), esdtrack->GetID());
      if (!track)
      {
        fh1NTracksReject->Fill("noTPConly", 1);
        continue;
      }
      AliExternalTrackParam exParam;
      Bool_t relate = track->RelateToVertexTPC(fVtx, fESD->GetMagneticField(), kVeryBig, &exParam);
      if (!relate)
      {
        fh1NTracksReject->Fill("relate", 1);
        if (track)
          delete track;
        continue;
      }
      track->Set(exParam.GetX(), exParam.GetAlpha(), exParam.GetParameter(), exParam.GetCovariance());
    }
    else if (fTrackType == 5 || fTrackType == 6)
    {
      if (fTrackCuts->AcceptTrack(esdtrack))
      {
        continue;
      }
      else
      {
        if (!(fTrackCutsITSLoose->AcceptTrack(esdtrack)) && fTrackCutsTPConly->AcceptTrack(esdtrack))
        {

          if (fTrackType == 5)
          {
            //use TPConly constrained track
            track = AliESDtrackCuts::GetTPCOnlyTrack(fESD, esdtrack->GetID());
            if (!track)
            {
              fh1NTracksReject->Fill("noTPConly", 1);
              continue;
            }
            AliExternalTrackParam exParam;
            Bool_t relate = track->RelateToVertexTPC(fVtx, fESD->GetMagneticField(), kVeryBig, &exParam);
            if (!relate)
            {
              fh1NTracksReject->Fill("relate", 1);
              if (track)
                delete track;
              continue;
            }
            track->Set(exParam.GetX(), exParam.GetAlpha(), exParam.GetParameter(), exParam.GetCovariance());
          }
          else if (fTrackType == 6)
          {
            //use global constrained track
            track = new AliESDtrack(*esdtrack);
            track->Set(esdtrack->GetConstrainedParam()->GetX(), esdtrack->GetConstrainedParam()->GetAlpha(), esdtrack->GetConstrainedParam()->GetParameter(), esdtrack->GetConstrainedParam()->GetCovariance());
          }
        }
      }
    }
    else if (fTrackType == 7)
    {
      //use global constrained track
      track = new AliESDtrack(*esdtrack);
    }
    else
      track = esdtrack;

    if (!track)
    {
      continue;
    }

    if (fTrackType == 2 || fTrackType == 4 || fTrackType == 5)
    {
      //Cut on chi2 of constrained fit
      if (track->GetConstrainedChi2TPC() > fSigmaConstrainedMax * fSigmaConstrainedMax && fSigmaConstrainedMax > 0.)
      {
        fh1NTracksReject->Fill("chi2", 1);
        if (track)
          delete track;
        continue;
      }
    }

    fPtAll->Fill(track->Pt());

    if (!(fTrackCuts->AcceptTrack(track)) && fTrackType != 4 && fTrackType != 5 && fTrackType != 6)
    {
      fh1NTracksReject->Fill("trackCuts", 1);
      if (fTrackType == 1 || fTrackType == 2 || fTrackType == 7)
      {
        if (track)
          delete track;
      }
      continue;
    }

    if (fTrackType == 7)
    {
      if (fTrackCutsITSLoose)
      {
        if (fTrackCutsITSLoose->AcceptTrack(track))
        {
          if (track)
            delete track;
          continue;
        }
      }

      if (esdtrack->GetConstrainedParam())
        track->Set(esdtrack->GetConstrainedParam()->GetX(), esdtrack->GetConstrainedParam()->GetAlpha(), esdtrack->GetConstrainedParam()->GetParameter(), esdtrack->GetConstrainedParam()->GetCovariance());
    }

    if (!track)
    {
      if (fTrackType == 1 || fTrackType == 2 || fTrackType == 4 || fTrackType == 5 || fTrackType == 6 || fTrackType == 7)
      {
        if (track)
          delete track;
      }
      continue;
    }

    fh1NTracksSel->Fill(0.);

    TArrayF &variables = *fVariables;
    variables.Reset(0.);

    variables[kVarPT] = track->Pt();
    variables[kVarPhi] = track->Phi();
    variables[kVarEta] = track->Eta();

    Float_t dca2D = 0.;
    Float_t dcaz = 0.;

    if (fTrackType == 1 || fTrackType == 2 || fTrackType == 4)
    {
      track->GetImpactParametersTPC(dca2D, dcaz); //TPConly
    }
    else
      track->GetImpactParameters(dca2D, dcaz); //Global

    variables[kVarDCAXY] = dca2D;
    variables[kVarDCAZ] = dcaz;

    variables[kVarClustersTPC] = (float)track->GetTPCNcls();

    Int_t nPointITS = 0;
    fITSClusterMap = track->GetITSClusterMap();
    UChar_t itsMap = track->GetITSClusterMap();
    for (Int_t i = 0; i < 6; i++)
    {
      if (itsMap & (1 << i))
        nPointITS++;
    }
    variables[kVarPointsITS] = (float)nPointITS;
    Float_t chi2C = (float)track->GetConstrainedChi2();
    if (fTrackType == 1 || fTrackType == 2 || fTrackType == 4)
      chi2C = (float)track->GetConstrainedChi2TPC();
    variables[kVarChi2C] = chi2C;
    variables[kVarNsigVtx] = fTrackCuts->GetSigmaToVertex(track); // Calculates the number of sigma to the vertex for a track.

    variables[kVarRelUncertaintyInvPt] = GetTrackLengthTPC(track);

    if (variables[kVarClustersTPC] > 0.)
      variables[kVarChi2ClusterTPC] = track->GetTPCchi2() / variables[kVarClustersTPC];

    //fVariables->SetAt(track->GetTPCClusterInfo(2,1),11); //#crossed rows
    variables[kVarCrossedRowsTPC] = track->GetTPCCrossedRows(); //#crossed rows

    Float_t crossedRowsTPCNClsF = 1.; //track->GetTPCClusterInfo(2,0);
    if (track->GetTPCNclsF() > 0.)
      crossedRowsTPCNClsF = variables[11] / track->GetTPCNclsF();
    variables[kVarCrossedRowsFindTPC] = crossedRowsTPCNClsF; //(#crossed rows)/(#findable clusters)
    variables[kVarSigmaY2] = track->GetSigmaY2();
    variables[kVarSigmaZ2] = track->GetSigmaZ2();
    variables[kVarSigmaSnp2] = track->GetSigmaSnp2();
    variables[kVarSigmaTgl2] = track->GetSigmaTgl2();
    variables[kVarSigmaInvPt2] = track->GetSigma1Pt2();

    variables[kVarClustersTPCIter1] = track->GetTPCNclsIter1();
    variables[kVarChi2TPCIter1] = track->GetTPCchi2Iter1();

    variables[kVarClustersTPCShared] = track->GetTPCnclsS();

    Float_t chi2Gold = origtrack->GetChi2TPCConstrainedVsGlobal(fVtx); //GetGoldenChi2(origtrack);
    Float_t chi2GGC = GetGGCChi2(origtrack.get());

    variables[kVarChi2Gold] = chi2Gold;
    variables[kVarChi2ConstrainedGlob] = chi2GGC;

    variables[kVarNCrossedRowsFit] = GetTPCClusterInfoFitMap(track, 2, 1);
    Float_t crossedRowsTPCNClsFFit = 1.;
    if (track->GetTPCNclsF() > 0.)
      crossedRowsTPCNClsFFit = variables[kVarNCrossedRowsFit] / track->GetTPCNclsF();
    variables[kVarCrossedRowsFindFit] = crossedRowsTPCNClsFFit;

    variables[kVarChi2ITS] = track->GetITSchi2();

    variables[kVarTPCclustersdEdx] = track->GetTPCsignalN();

    TBits fitmap = track->GetTPCFitMap();
    fPtNClustersNClustersFitMap->Fill(track->Pt(), track->GetTPCNcls(), (float)fitmap.CountBits());

    FillHistograms();

    //      int mult = fTrackCuts->CountAcceptedTracks(fESD);

    if (fTrackType == 1 || fTrackType == 2 || fTrackType == 4 || fTrackType == 5 || fTrackType == 6 || fTrackType == 7)
    {
      if (track)
        delete track;
    }

  } //track loop
}

 /**
  * @brief Do QA on AOD input
  */
void AliPWG4HighPtTrackQA::DoAnalysisAOD()
{
  AliAODEvent *aod = dynamic_cast<AliAODEvent *>(fEvent);
  if (!aod)
    return;

  // ---- Get MC Header information (for MC productions in pThard bins) ----
  Double_t ptHard = 0.;
  Double_t nTrials = 1; // trials for MC trigger weight for real data

  if (MCEvent())
  {
    AliGenPythiaEventHeader *pythiaGenHeader = GetPythiaEventHeader(MCEvent());
    if (pythiaGenHeader)
    {
      nTrials = pythiaGenHeader->Trials();
      ptHard = pythiaGenHeader->GetPtHard();

      fh1PtHard->Fill(ptHard);
      fh1PtHardTrials->Fill(ptHard, nTrials);

      fh1Trials->Fill("#sum{ntrials}", fAvgTrials);
    }
  }

  AliExternalTrackParam exParam;
  for (Int_t iTrack = 0; iTrack < fEvent->GetNumberOfTracks(); iTrack++)
  {

    AliAODTrack *aodtrack = dynamic_cast<AliAODTrack *>(aod->GetTrack(iTrack));
    if (!aodtrack)
      AliFatal("Not a standard AOD");
    AliDebugStream(2) << GetName() << ": Applying filter bits " << fFilterMask << std::endl;
    if (!aodtrack->TestFilterBit(fFilterMask))
    {
      fh1NTracksReject->Fill("noHybridTrack", 1);
      continue;
    }

    if (!fIncludeNoITS)
    {
      if ((aodtrack->GetStatus() & AliESDtrack::kITSrefit) == 0)
      {
        fh1NTracksReject->Fill("noITSrefit", 1);
        continue;
      }
    }

    TArrayF & variables = *fVariables;
    variables.Reset(0.);

    variables[kVarPT] = aodtrack->Pt();
    variables[kVarPhi] = aodtrack->Phi();
    variables[kVarEta] = aodtrack->Eta();

    Double_t dca[2] = {0., 0.};
    if (aodtrack->IsGlobalConstrained())
    {
      dca[0] = aodtrack->DCA();
      dca[1] = aodtrack->ZAtDCA();
    }
    else
    {
      Double_t v[3] = {0};
      Double_t pos[3] = {0};
      fVtxAOD->GetXYZ(v);
      aodtrack->GetXYZ(pos);
      dca[0] = pos[0] - v[0];
      dca[1] = pos[1] - v[1];
    }
    variables[kVarDCAXY] = dca[0];
    variables[kVarDCAZ] = dca[1];
    variables[kVarClustersTPC] = (float)aodtrack->GetTPCNcls();
    variables[kVarPointsITS] = (float)aodtrack->GetITSNcls();
    variables[kVarChi2C] = 0.; //ConstrainedChi2TPC -> not available in AOD
    variables[kVarNsigVtx] = 0.;
    variables[kVarRelUncertaintyInvPt] = GetTrackLengthTPC(aodtrack);
    Float_t chi2pndf = aodtrack->Chi2perNDF();
    //if(fVariables->At(5)>0.) chi2pndf = aodtrack->GetTPCchi2()/fVariables->At(5);
    variables[kVarChi2ClusterTPC] = chi2pndf;
    variables[kVarCrossedRowsTPC] = GetTPCClusterInfo(aodtrack, 2, 1, 0, 159, kFALSE);
    Float_t crossedRowsTPCNClsF = 0.;
    if (aodtrack->GetTPCNclsF() > 0.)
      crossedRowsTPCNClsF = variables[11] / aodtrack->GetTPCNclsF();
    variables[kVarCrossedRowsFindTPC] = crossedRowsTPCNClsF;

    //get covariance matrix
    Double_t cov[21] = {
        0,
    };
    aodtrack->GetCovMatrix(cov);
    Double_t pxpypz[3] = {
        0,
    };
    aodtrack->PxPyPz(pxpypz);
    Double_t xyz[3] = {
        0,
    };
    aodtrack->GetXYZ(xyz);
    Short_t sign = aodtrack->Charge();
    exParam.Set(xyz, pxpypz, cov, sign);

    variables[kVarSigmaY2] = exParam.GetSigmaY2();
    variables[kVarSigmaZ2] = exParam.GetSigmaZ2();
    variables[kVarSigmaSnp2] = exParam.GetSigmaSnp2();
    variables[kVarSigmaTgl2] = exParam.GetSigmaTgl2();
    variables[kVarSigmaInvPt2] = exParam.GetSigma1Pt2();

    variables[kVarClustersTPCIter1] = 0.; //NClustersTPCIter1
    variables[kVarChi2TPCIter1] = 0.; //Chi2TPCIter1

    TBits sharedClusterMap = aodtrack->GetTPCSharedMap();
    variables[kVarClustersTPCShared] = sharedClusterMap.CountBits();

    variables[kVarChi2Gold] = 0.; //not available in AOD golden chi2
    variables[kVarChi2ConstrainedGlob] = 0.; //not available in AOD  Chi2 between global and global constrained

    variables[kVarNCrossedRowsFit] = GetTPCClusterInfo(aodtrack, 2, 1, 0, 159, kTRUE); //not available in AOD #crossed rows from fit map
    Float_t crossedRowsTPCNClsFFit = 0.;
    if (aodtrack->GetTPCNclsF() > 0.)
      crossedRowsTPCNClsFFit = variables[kVarNCrossedRowsFit] / aodtrack->GetTPCNclsF();
    variables[kVarCrossedRowsFindFit] = crossedRowsTPCNClsFFit; //(#crossed rows)/(#findable clusters) from fit map

    variables[kVarChi2ITS] = 0.;

    variables[kVarTPCclustersdEdx] = aodtrack->GetTPCsignalN();

    fPtAll->Fill(variables[0]);

    FillHistograms();
  }
}

/**
 * @brief Fill all QA histograms
 */
void AliPWG4HighPtTrackQA::FillHistograms()
{
  // declare references to array entries
  TArrayF &variables = *fVariables;
  const Float_t &pt = variables[kVarPT], &phi = variables[kVarPhi], &eta = variables[kVarEta], &dcaXY = variables[kVarDCAXY], 
                &dcaZ = variables[kVarDCAZ], &nclustersTPC = variables[kVarClustersTPC], &npointsITS = variables[kVarPointsITS], 
                &chi2constrainedTPC = variables[kVarChi2C], &sigToVtx = variables[kVarNsigVtx], 
                &relUncertaintyInvPt = variables[kVarRelUncertaintyInvPt], &chi2clusterTPC = variables[kVarChi2ClusterTPC], 
                &ncrossedrowsTPC = variables[kVarCrossedRowsTPC], &crossedrowsoverfindable = variables[kVarCrossedRowsFindTPC], 
                &sigmaY2 = variables[kVarSigmaY2], &sigmaZ2 = variables[kVarSigmaZ2], &sigmaSnp2 = variables[kVarSigmaSnp2],
                &sigmaTgl2 = variables[kVarSigmaTgl2], &sigmaInvPt = variables[kVarSigmaInvPt2], &nclustersTPCIter1 = variables[kVarClustersTPCIter1], 
                &chi2TPCIter1 = variables[kVarChi2TPCIter1], &nclustersTPCShared = variables[kVarClustersTPCShared], &chi2Gold = variables[kVarChi2Gold], 
                &chi2GlobalConstrained = variables[kVarChi2ConstrainedGlob], &ncrossedrowsFit = variables[kVarNCrossedRowsFit], 
                &ncrossedrowsFindableFit = variables[kVarCrossedRowsFindFit], &chi2ITS = variables[kVarChi2ITS], &nclustersdEdxTPC = variables[kVarTPCclustersdEdx];

  // Fill the histograms using the references
  fPtSel->Fill(pt);
  fPtPhi->Fill(pt, phi);
  fPtEta->Fill(pt, eta);
  fPtEtaPhi->Fill(pt, eta, phi);
  fPtDCA2D->Fill(pt, dcaXY);
  fPtDCAZ->Fill(pt, dcaZ);
  fPtNClustersTPC->Fill(pt, nclustersTPC);
  fPtNClustersTPCPhi->Fill(phi, nclustersTPC);
  fPtNPointITS->Fill(pt, npointsITS);
  fPtNPointITSPhi->Fill(pt, npointsITS, phi);

  fPtNClustersTPCIter1->Fill(pt, nclustersTPCIter1);
  fPtNClustersTPCIter1Phi->Fill(pt, nclustersTPCIter1, phi);
  fPtNClustersTPCShared->Fill(pt, nclustersTPCShared);
  if (nclustersTPC > 0.)
    fPtNClustersTPCSharedFrac->Fill(pt, nclustersTPCShared / nclustersTPC);

  if (nclustersTPCIter1 > 0.)
    fPtChi2PerClusterTPCIter1->Fill(pt, chi2TPCIter1 / nclustersTPCIter1);

  fPtChi2C->Fill(pt, chi2constrainedTPC);
  fPtNSigmaToVertex->Fill(pt, sigToVtx);
  fPtRelUncertainty1Pt->Fill(pt, pt * TMath::Sqrt(sigmaInvPt));
  fPtRelUncertainty1PtNClus->Fill(pt, pt * TMath::Sqrt(sigmaInvPt), nclustersTPC);
  fPtRelUncertainty1PtNClusIter1->Fill(pt, pt * TMath::Sqrt(sigmaInvPt), nclustersTPCIter1);
  fPtRelUncertainty1PtNPointITS->Fill(pt, pt * TMath::Sqrt(sigmaInvPt), npointsITS);

  fPtRelUncertainty1PtITSClusterMap->Fill(pt, pt * TMath::Sqrt(sigmaInvPt), (int)fITSClusterMap);

  fPtRelUncertainty1PtChi2->Fill(pt, pt * TMath::Sqrt(sigmaInvPt), chi2clusterTPC);
  if (nclustersTPCIter1 > 0.)
    fPtRelUncertainty1PtChi2Iter1->Fill(pt, pt * TMath::Sqrt(sigmaInvPt), chi2TPCIter1 / nclustersTPCIter1);
  fPtRelUncertainty1PtPhi->Fill(pt, pt * TMath::Sqrt(sigmaInvPt), phi);

  fPtSigmaY2->Fill(1. / pt, TMath::Sqrt(sigmaY2));
  fPtSigmaZ2->Fill(1. / pt, TMath::Sqrt(sigmaZ2));
  fPtSigmaSnp2->Fill(1. / pt, TMath::Sqrt(sigmaSnp2));
  fPtSigmaTgl2->Fill(1. / pt, TMath::Sqrt(sigmaTgl2));
  fPtSigma1Pt2->Fill(1. / pt, TMath::Sqrt(sigmaInvPt));

  fProfPtSigmaY2->Fill(1. / pt, TMath::Sqrt(sigmaY2));
  fProfPtSigmaZ2->Fill(1. / pt, TMath::Sqrt(sigmaZ2));
  fProfPtSigmaSnp2->Fill(1. / pt, TMath::Sqrt(sigmaSnp2));
  fProfPtSigmaTgl2->Fill(1. / pt, TMath::Sqrt(sigmaTgl2));
  fProfPtSigma1Pt2->Fill(1. / pt, TMath::Sqrt(sigmaInvPt));
  fProfPtSigma1Pt->Fill(pt, TMath::Sqrt(sigmaInvPt));
  fProfPtPtSigma1Pt->Fill(pt, pt * TMath::Sqrt(sigmaInvPt));

  fPtChi2PerClusterTPC->Fill(pt, chi2clusterTPC);
  fPtNCrossedRows->Fill(pt, ncrossedrowsTPC);
  fPtNCrossedRowsPhi->Fill(pt, ncrossedrowsTPC, phi);
  fPtNCrossedRowsNClusFPhi->Fill(pt, crossedrowsoverfindable, phi);
  fPtNCrRNCrRNClusF->Fill(pt, ncrossedrowsTPC, crossedrowsoverfindable);

  fPtChi2Gold->Fill(pt, chi2Gold);
  fPtChi2GGC->Fill(pt, chi2GlobalConstrained);

  fPtChi2GoldPhi->Fill(pt, chi2Gold, phi);
  fPtChi2GGCPhi->Fill(pt, chi2GlobalConstrained, phi);

  fChi2GoldChi2GGC->Fill(chi2Gold, chi2GlobalConstrained);

  fPtNCrossedRowsFit->Fill(pt, ncrossedrowsFit);
  fPtNCrossedRowsFitPhi->Fill(pt, ncrossedrowsFit, phi);
  fPtNCrossedRowsNClusFFitPhi->Fill(pt, ncrossedrowsFindableFit, phi);
  fNCrossedRowsNCrossedRowsFit->Fill(ncrossedrowsTPC, ncrossedrowsFit);

  fNClustersNCrossedRows->Fill(nclustersTPC, ncrossedrowsTPC);
  fNClustersNCrossedRowsFit->Fill(nclustersTPC, ncrossedrowsFit);

  fPtRelUncertainty1PtNCrossedRows->Fill(pt, pt * TMath::Sqrt(sigmaInvPt), ncrossedrowsTPC);
  fPtRelUncertainty1PtNCrossedRowsFit->Fill(pt, pt * TMath::Sqrt(sigmaInvPt), ncrossedrowsFit);

  if (npointsITS > 0.)
    fPtChi2ITSPhi->Fill(pt, variables[25] / npointsITS, phi);

  fPtTPCSignalN->Fill(pt, variables[26]);
}

/**
 * @brief Get the cross section and the trails either from pyxsec.root or from pysec_hists.root
 * 
 * This is to called in Notify and should provide the path to the AOD/ESD file
 * Copied from AliAnalysisTaskJetSpectrum2
 */
Bool_t AliPWG4HighPtTrackQA::PythiaInfoFromFile(const char *currFile, Float_t &fXsec, Float_t &fTrials)
{
  TString file(currFile);
  fXsec = 0;
  fTrials = 1;

  if (file.Contains("root_archive.zip#"))
  {
    Ssiz_t pos1 = file.Index("root_archive", 12, TString::kExact);
    Ssiz_t pos = file.Index("#", 1, pos1, TString::kExact);
    file.Replace(pos + 1, 20, "");
  }
  else if (file.Contains("aod_archive.zip#"))
  {
    Ssiz_t pos1 = file.Index("aod_archive", 11, TString::kExact);
    Ssiz_t pos = file.Index("#", 1, pos1, TString::kExact);
    file.Replace(pos + 1, 20, "");
  }
  else
  {
    // not an archive take the basename....
    file.ReplaceAll(gSystem->BaseName(file.Data()), "");
  }

  std::unique_ptr<TFile> fxsec(TFile::Open(Form("%s%s", file.Data(), "pyxsec.root"))); // problem that we cannot really test the existance of a file in a archive so we have to lvie with open error message from root
  if (!fxsec)
  {
    // next trial fetch the histgram file
    fxsec = std::unique_ptr<TFile>(TFile::Open(Form("%s%s", file.Data(), "pyxsec_hists.root")));
    if (!fxsec)
    {
      // not a severe condition but inciate that we have no information
      return kFALSE;
    }
    else
    {
      // find the tlist we want to be independtent of the name so use the Tkey
      TKey *key = (TKey *)fxsec->GetListOfKeys()->At(0);
      if (!key)
      {
        fxsec->Close();
        return kFALSE;
      }
      TList *list = dynamic_cast<TList *>(key->ReadObj());
      if (!list)
      {
        fxsec->Close();
        return kFALSE;
      }
      fXsec = ((TProfile *)list->FindObject("h1Xsec"))->GetBinContent(1);
      fTrials = ((TH1F *)list->FindObject("h1Trials"))->GetBinContent(1);
      fxsec->Close();
    }
  } // no tree pyxsec.root
  else
  {
    TTree *xtree = (TTree *)fxsec->Get("Xsection");
    if (!xtree)
    {
      fxsec->Close();
      return kFALSE;
    }
    UInt_t ntrials = 0;
    Double_t xsection = 0;
    xtree->SetBranchAddress("xsection", &xsection);
    xtree->SetBranchAddress("ntrials", &ntrials);
    xtree->GetEntry(0);
    fTrials = ntrials;
    fXsec = xsection;
    fxsec->Close();
  }
  return kTRUE;
}

/**
 * @brief Notify function called when file is changed
 *
 * Implemented Notify() to read the cross sections
 * and number of trials from pyxsec.root
 * Copied from AliAnalysisTaskJetSpectrum2
 * 
 * @return Bool_t false if failure, otherwise true
 */
Bool_t AliPWG4HighPtTrackQA::Notify()
{
  TTree *tree = AliAnalysisManager::GetAnalysisManager()->GetTree();
  Float_t xsection = 0;
  Float_t ftrials = 1;

  fAvgTrials = 1;
  if (tree)
  {
    TFile *curfile = tree->GetCurrentFile();
    if (!curfile)
    {
      Error("Notify", "No current file");
      return kFALSE;
    }
    if (!fh1Xsec || !fh1Trials)
    {
      //      Printf("%s%d No Histogram fh1Xsec",(char*)__FILE__,__LINE__);
      return kFALSE;
    }
    PythiaInfoFromFile(curfile->GetName(), xsection, ftrials);
    fh1Xsec->Fill("<#sigma>", xsection);
    // construct a poor man average trials
    Float_t nEntries = (Float_t)tree->GetTree()->GetEntries();
    if (ftrials >= nEntries && nEntries > 0.)
      fAvgTrials = ftrials / nEntries;
  }
  return kTRUE;
}

AliGenPythiaEventHeader *AliPWG4HighPtTrackQA::GetPythiaEventHeader(const AliMCEvent *mcEvent)
{

  if (!mcEvent)
    return 0;
  AliGenEventHeader *genHeader = mcEvent->GenEventHeader();
  AliGenPythiaEventHeader *pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader *>(genHeader);
  if (!pythiaGenHeader)
  {
    // cocktail ??
    AliGenCocktailEventHeader *genCocktailHeader = dynamic_cast<AliGenCocktailEventHeader *>(genHeader);

    if (!genCocktailHeader)
    {
      //      AliWarningGeneral(Form(" %s:%d",(char*)__FILE__,__LINE__),"Unknown header type (not Pythia or Cocktail)");
      //      AliWarning(Form("%s %d: Unknown header type (not Pythia or Cocktail)",(char*)__FILE__,__LINE__));
      return 0;
    }
    TList *headerList = genCocktailHeader->GetHeaders();
    for (Int_t i = 0; i < headerList->GetEntries(); i++)
    {
      pythiaGenHeader = dynamic_cast<AliGenPythiaEventHeader *>(headerList->At(i));
      if (pythiaGenHeader)
        break;
    }
    if (!pythiaGenHeader)
    {
      AliWarningGeneral(Form(" %s:%d", (char *)__FILE__, __LINE__), "Pythia event header not found");
      return 0;
    }
  }
  return pythiaGenHeader;
}

/**
 * @brief Get TPC cluster information
 *
 * TPC cluster information
 *   type 0: get fraction of found/findable clusters with neighbourhood definition
 *        1: findable clusters with neighbourhood definition
 *        2: found clusters
 *
 *  definition of findable clusters:
 *             a cluster is defined as findable if there is another cluster
 *             within +- nNeighbours pad rows. The idea is to overcome threshold
 *             effects with a very simple algorithm.
 * 
 * MV: copied from AliESDtrack since method is not available in AliAODTrack
 */
Float_t AliPWG4HighPtTrackQA::GetTPCClusterInfo(const AliAODTrack *tr, Int_t nNeighbours /*=3*/, Int_t type /*=0*/, Int_t row0, Int_t row1, Bool_t useFitMap) const
{
  TBits fTPCClusterMap = 0;
  if (useFitMap)
    fTPCClusterMap = tr->GetTPCFitMap();
  else
    fTPCClusterMap = tr->GetTPCClusterMap();

  if (type == 2)
    return fTPCClusterMap.CountBits();

  Int_t found = 0;
  Int_t findable = 0;
  Int_t last = -nNeighbours;

  for (Int_t i = row0; i < row1; ++i)
  {
    //look to current row
    if (fTPCClusterMap[i])
    {
      last = i;
      ++found;
      ++findable;
      continue;
    }
    //look to nNeighbours before
    if ((i - last) <= nNeighbours)
    {
      ++findable;
      continue;
    }
    //look to nNeighbours after
    for (Int_t j = i + 1; j < i + 1 + nNeighbours; ++j)
    {
      if (fTPCClusterMap[j])
      {
        ++findable;
        break;
      }
    }
  }
  if (type == 1)
    return findable;

  if (type == 0)
  {
    Float_t fraction = 0;
    if (findable > 0)
      fraction = (Float_t)found / (Float_t)findable;
    else
      fraction = 0;
    return fraction;
  }
  return 0; // undefined type - default value
}

/**
 * @brief Get TPC cluster information from fit map
 * 
 * TPC cluster information from fit map
 * type 0: get fraction of found/findable clusters with neighbourhood definition
 *      1: findable clusters with neighbourhood definition
 *      2: found clusters
 * 
 * definition of findable clusters:
 *            a cluster is defined as findable if there is another cluster
 *            within +- nNeighbours pad rows. The idea is to overcome threshold
 *            effects with a very simple algorithm.
 * 
 * @param tr Track to check
 * @param nNeighbours Number of neighbors
 * @param type Track type
 * @param row0 first row
 * @param row1 last row
 * @return Number of clusters from the fit 
 */
Float_t AliPWG4HighPtTrackQA::GetTPCClusterInfoFitMap(const AliESDtrack *tr, Int_t nNeighbours /*=3*/, Int_t type /*=0*/, Int_t row0, Int_t row1) const
{
  TBits fTPCFitMap = tr->GetTPCFitMap();
  if (type == 2)
    return fTPCFitMap.CountBits();

  Int_t found = 0;
  Int_t findable = 0;
  Int_t last = -nNeighbours;

  for (Int_t i = row0; i < row1; ++i)
  {
    //look to current row
    if (fTPCFitMap[i])
    {
      last = i;
      ++found;
      ++findable;
      continue;
    }
    //look to nNeighbours before
    if ((i - last) <= nNeighbours)
    {
      ++findable;
      continue;
    }
    //look to nNeighbours after
    for (Int_t j = i + 1; j < i + 1 + nNeighbours; ++j)
    {
      if (fTPCFitMap[j])
      {
        ++findable;
        break;
      }
    }
  }
  if (type == 1)
    return findable;

  if (type == 0)
  {
    Float_t fraction = 0;
    if (findable > 0)
      fraction = (Float_t)found / (Float_t)findable;
    else
      fraction = 0;
    return fraction;
  }
  return 0; // undefined type - default value
}

/**
 * @brief  Get track length
 * 
 * returns distance between 1st and last hit in TPC
 * distance given in number of padrows
 * 
 * @param track Track to check
 * @return Length of the track
 */
Int_t AliPWG4HighPtTrackQA::GetTrackLengthTPC(const AliESDtrack *track) const
{
  

  TBits fTPCClusterMap = track->GetTPCClusterMap();
  int firstHit = 0;
  int lastHit = 0;

  for (int i = 0; i <= 159; i++)
  {
    if (fTPCClusterMap[i] > 0)
      firstHit = i;
  }
  for (int i = 159; i >= 0; i--)
  {
    if (fTPCClusterMap[i] > 0)
      lastHit = i;
  }

  Int_t trackLength = lastHit - firstHit;

  return trackLength;
}

/**
 * @brief Get the TPC track length
 * 
 * returns distance between 1st and last hit in TPC
 * distance given in number of padrows
 *
 * @param track Track to check
 * @return Int_t Length of the track inside TPC
 */
Int_t AliPWG4HighPtTrackQA::GetTrackLengthTPC(const AliAODTrack *track) const
{

  TBits fTPCClusterMap = track->GetTPCClusterMap();
  int firstHit = 0;
  int lastHit = 0;

  for (int i = 0; i <= 159; i++)
  {
    if (fTPCClusterMap[i] > 0)
      firstHit = i;
  }
  for (int i = 159; i >= 0; i--)
  {
    if (fTPCClusterMap[i] > 0)
      lastHit = i;
  }

  Int_t trackLength = lastHit - firstHit;

  return trackLength;
}

/**
 * @brief Get the golden chi2 for the given track
 * 
 * Return chi2 between global and TPC constrained track
 * track should be the global unconstrained track
 *
 * @param origtrack Track to check
 * @return Float_t Golden chi2
 */
Float_t AliPWG4HighPtTrackQA::GetGoldenChi2(AliESDtrack *origtrack)
{
  Float_t chi2Gold = 0.;

  std::unique_ptr<AliESDtrack> tpcTrack(AliESDtrackCuts::GetTPCOnlyTrack(fESD, origtrack->GetID()));
  if (tpcTrack)
  {
    AliExternalTrackParam exParam;
    Bool_t relate = tpcTrack->RelateToVertexTPC(fVtx, fESD->GetMagneticField(), kVeryBig, &exParam);
    if (relate)
    {
      tpcTrack->Set(exParam.GetX(), exParam.GetAlpha(), exParam.GetParameter(), exParam.GetCovariance());
      //	  Double_t pTPC[2],covTPC[3];	  tpcTrack->PropagateToDCA(fVtx, fESD->GetMagneticField(), 10000,  pTPC, covTPC);
    }

    tpcTrack->Propagate(origtrack->GetAlpha(), origtrack->GetX(), fESD->GetMagneticField());
    chi2Gold = (Float_t)origtrack->GetPredictedChi2(tpcTrack.get());
  }
  return chi2Gold;
}

/**
 * @brief Get the global constrained chi2 of the track
 * 
 * Return chi2 between global and global constrained track
 * track should be the global unconstrained track
 *
 * @param origtrack Track to check 
 * @return Float_t Global constrained chi2
 */
Float_t AliPWG4HighPtTrackQA::GetGGCChi2(AliESDtrack *origtrack)
{
  Float_t chi2GGC = 0.;

  std::unique_ptr<AliESDtrack> esdtrackC(new AliESDtrack(*origtrack));
  if (esdtrackC)
  {
    if (origtrack->GetConstrainedParam())
    {
      esdtrackC->Set(origtrack->GetConstrainedParam()->GetX(), origtrack->GetConstrainedParam()->GetAlpha(), origtrack->GetConstrainedParam()->GetParameter(), origtrack->GetConstrainedParam()->GetCovariance());
      chi2GGC = (Float_t)origtrack->GetPredictedChi2(esdtrackC.get());
    }
  }

  return chi2GGC;
}

/**
 * @brief Terminate function
 * 
 * The Terminate() function is the last function to be called during
 * a query. It always runs on the client, it can be used to present
 * the results graphically or save the results to file.
 * 
 * @param Option_t Unused
 */
void AliPWG4HighPtTrackQA::Terminate(Option_t *)
{
}

#endif
