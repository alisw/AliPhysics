/************************************************************************************
 * Copyright (C) 2019, Copyright Holders of the ALICE Collaboration                 *
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
#include <iostream>
#include <memory>
#include <vector>

#include <TArrayD.h>
#include <TBinning.h>
#include <TCustomBinning.h>
#include <TH1.h>
#include <TH2.h>
#include <TLinearBinning.h>
#include <TLorentzVector.h>
#include <TRandom.h>

#include <fastjet/ClusterSequence.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/contrib/SoftDrop.hh>
#include <fastjet/config.h>
#if FASJET_VERSION_NUMBER >= 30302
#include <fastjet/tools/Recluster.hh>
#else
#include <fastjet/contrib/Recluster.hh>
#endif

#include <RooUnfoldResponse.h>

#include <AliAODEvent.h>
#include <AliAODInputHandler.h>
#include <AliAnalysisManager.h>
#include "AliAnalysisTaskEmcalSoftDropResponse.h"
#include <AliClusterContainer.h>
#include <AliEmcalJet.h>
#include <AliEmcalAnalysisFactory.h>
#include <AliJetContainer.h>
#include <AliMCParticleContainer.h>
#include <AliTrackContainer.h>
#include <AliVCluster.h>
#include <AliVParticle.h>
#include <AliVTrack.h>

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalSoftDropResponse)

    using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalSoftDropResponse::AliAnalysisTaskEmcalSoftDropResponse() : AliAnalysisTaskEmcalJet(),
                                                                               fBinningMode(kSDModeINT7),
                                                                               fFractionResponseClosure(0.5),
                                                                               fZcut(0.1),
                                                                               fBeta(0.),
                                                                               fReclusterizer(kCAAlgo),
                                                                               fSampleFraction(1.),
                                                                               fMinFractionShared(0),
                                                                               fHasResponseMatrixSparse(false),
                                                                               fHasResponseMatrixRooUnfold(true),
                                                                               fUseChargedConstituents(true),
                                                                               fUseNeutralConstituents(true),
                                                                               fUseStandardOutlierRejection(false),
                                                                               fJetTypeOutliers(kOutlierPartJet),
                                                                               fNameMCParticles("mcparticles"),
                                                                               fSampleSplitter(nullptr),
                                                                               fSampleTrimmer(nullptr),
                                                                               fPartLevelPtBinning(nullptr),
                                                                               fDetLevelPtBinning(nullptr),
                                                                               fIsEmbeddedEvent(false),
                                                                               fNamePartLevelJetContainer(""),
                                                                               fNameDetLevelJetContainer(""),
                                                                               fNameUnSubLevelJetContainer(""),
                                                                               fZgResponse(),
                                                                               fZgResponseClosure(),
                                                                               fRgResponse(),
                                                                               fRgResponseClosure(),
                                                                               fNsdResponse(),
                                                                               fNsdResponseClosure(),
                                                                               fThetagResponse(),
                                                                               fThetagResponseClosure(),
                                                                               fHistManager("AliAnalysisTaskSoftDropResponse")
{
}

AliAnalysisTaskEmcalSoftDropResponse::AliAnalysisTaskEmcalSoftDropResponse(const char *name) : AliAnalysisTaskEmcalJet(name, kTRUE),
                                                                                               fBinningMode(kSDModeINT7),
                                                                                               fFractionResponseClosure(0.5),
                                                                                               fZcut(0.1),
                                                                                               fBeta(0.),
                                                                                               fReclusterizer(kCAAlgo),
                                                                                               fSampleFraction(1.),
                                                                                               fMinFractionShared(0),
                                                                                               fHasResponseMatrixSparse(false),
                                                                                               fHasResponseMatrixRooUnfold(true),
                                                                                               fUseChargedConstituents(true),
                                                                                               fUseNeutralConstituents(true),
                                                                                               fUseStandardOutlierRejection(false),
                                                                                               fJetTypeOutliers(kOutlierPartJet),
                                                                                               fNameMCParticles("mcparticles"),
                                                                                               fSampleSplitter(nullptr),
                                                                                               fSampleTrimmer(nullptr),
                                                                                               fPartLevelPtBinning(nullptr),
                                                                                               fDetLevelPtBinning(nullptr),
                                                                                               fNamePartLevelJetContainer(""),
                                                                                               fNameDetLevelJetContainer(""),
                                                                                               fNameUnSubLevelJetContainer(""),
                                                                                               fIsEmbeddedEvent(false),
                                                                                               fZgResponse(),
                                                                                               fZgResponseClosure(),
                                                                                               fRgResponse(),
                                                                                               fRgResponseClosure(),
                                                                                               fNsdResponse(),
                                                                                               fNsdResponseClosure(),
                                                                                               fThetagResponse(),
                                                                                               fThetagResponseClosure(),
                                                                                               fHistManager(name)
{
  SetMakeGeneralHistograms(true);
}

AliAnalysisTaskEmcalSoftDropResponse::~AliAnalysisTaskEmcalSoftDropResponse()
{
  if (fPartLevelPtBinning)
    delete fPartLevelPtBinning;
  if (fDetLevelPtBinning)
    delete fDetLevelPtBinning;
  if (fSampleSplitter)
    delete fSampleSplitter;
  if (fSampleTrimmer)
    delete fSampleTrimmer;
}

void AliAnalysisTaskEmcalSoftDropResponse::UserCreateOutputObjects()
{
  AliAnalysisTaskEmcalJet::UserCreateOutputObjects();
  double R = GetJetContainer(fNameDetLevelJetContainer.Data())->GetJetRadius();

  fSampleSplitter = new TRandom;
  if (fSampleFraction < 1.)
    fSampleTrimmer = new TRandom;

  if (!fPartLevelPtBinning)
    fPartLevelPtBinning = GetDefaultPartLevelPtBinning();
  if (!fDetLevelPtBinning)
    fDetLevelPtBinning = GetDefaultDetLevelPtBinning();
  std::unique_ptr<TBinning> zgbinning(GetZgBinning()),
                            rgbinning(GetRgBinning(R)),
                            nsdbinning(new TLinearBinning(22, -1.5, 20.5)),       // Negative bins are for untagged jets
                            thetagbinning(new TLinearBinning(11, -0.1, 1.)),
                            ptbinningFine(new TLinearBinning(500, 0., 500.));
  TArrayD binEdgesZg, binEdgesRg, binEdgesNsd, binEdgesThetag, binEdgesPtPart, binEdgesPtDet, binEdgesPtFine;
  zgbinning->CreateBinEdges(binEdgesZg);
  rgbinning->CreateBinEdges(binEdgesRg);
  nsdbinning->CreateBinEdges(binEdgesNsd);
  thetagbinning->CreateBinEdges(binEdgesThetag);
  ptbinningFine->CreateBinEdges(binEdgesPtFine);
  fPartLevelPtBinning->CreateBinEdges(binEdgesPtPart);
  fDetLevelPtBinning->CreateBinEdges(binEdgesPtDet);

  const TBinning *sparsebinningZg[4] = {zgbinning.get(), ptbinningFine.get(), zgbinning.get(), ptbinningFine.get()},
                 *sparsebinningRg[4] = {rgbinning.get(), ptbinningFine.get(), rgbinning.get(), ptbinningFine.get()},
                 *sparsebinningNsd[4] = {nsdbinning.get(), ptbinningFine.get(), nsdbinning.get(), ptbinningFine.get()},
                 *sparsebinningThetag[4] = {thetagbinning.get(), ptbinningFine.get(), thetagbinning.get(), ptbinningFine.get()};

  //Need to do centrality bins in the histograms if it is not pp
  if (fForceBeamType != kpp)
  {
    for (Int_t cent = 0; cent < fNcentBins; cent++)
    {
      if(fHasResponseMatrixRooUnfold){
        fHistManager.CreateTH2(Form("hZgDetLevel_%d", cent), Form("Zg response at detector level, %d centrality bin", cent), binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
        fHistManager.CreateTH2(Form("hZgPartLevel_%d", cent), Form("Zg response at particle level, %d centrality bin", cent), binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
        fHistManager.CreateTH2(Form("hZgPartLevelTruncated_%d", cent), Form("Zg response at particle level (truncated), %d centrality bin", cent), binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
        fHistManager.CreateTH2(Form("hRgDetLevel_%d", cent), Form("Rg response at detector level, %d centrality bin", cent), binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
        fHistManager.CreateTH2(Form("hRgPartLevel_%d", cent), Form("Rg response at particle level, %d centrality bin", cent), binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
        fHistManager.CreateTH2(Form("hRgPartLevelTruncated_%d", cent), Form("Rg response at particle level (truncated), %d centrality bin", cent), binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
        fHistManager.CreateTH2(Form("hNsdDetLevel_%d", cent), Form("Nsd response at detector level, %d centrality bin", cent), binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
        fHistManager.CreateTH2(Form("hNsdPartLevel_%d", cent), Form("Nsd response at particle level, %d centrality bin", cent), binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
        fHistManager.CreateTH2(Form("hNsdPartLevelTruncated_%d", cent), Form("Nsd response at particle level (truncated), %d centrality bin", cent), binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
        fHistManager.CreateTH2(Form("hThetagDetLevel_%d", cent), Form("Thetag response at detector level, %d centrality bin", cent), binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
        fHistManager.CreateTH2(Form("hThetagPartLevel_%d", cent), Form("Thetag response at particle level, %d centrality bin", cent), binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
        fHistManager.CreateTH2(Form("hThetagPartLevelTruncated_%d", cent), Form("Thetag response at particle level (truncated), %d centrality bin", cent), binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());

        // For closure test
        fHistManager.CreateTH2(Form("hZgPartLevelClosureNoResp_%d", cent), Form("Zg response at particle level (closure test, jets not used for the response matrix), %d centrality bin", cent), binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
        fHistManager.CreateTH2(Form("hZgDetLevelClosureNoResp_%d", cent), Form("Zg response at detector level (closure test, jets not used for the response matrix), %d centrality bin", cent), binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
        fHistManager.CreateTH2(Form("hZgPartLevelClosureResp_%d", cent), Form("Zg response at particle level (closure test, jets used for the response matrix), %d centrality bin", cent), binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
        fHistManager.CreateTH2(Form("hZgDetLevelClosureResp_%d", cent), Form("Zg response at detector level (closure test, jets used for the response matrix), %d centrality bin", cent), binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
        fHistManager.CreateTH2(Form("hRgPartLevelClosureNoResp_%d", cent), Form("Rg response at particle level (closure test, jets not used for the response matrix), %d centrality bin", cent), binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
        fHistManager.CreateTH2(Form("hRgDetLevelClosureNoResp_%d", cent), Form("Rg response at detector level (closure test, jets not used for the response matrix), %d centrality bin", cent), binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
        fHistManager.CreateTH2(Form("hRgPartLevelClosureResp_%d", cent), Form("Rg response at particle level (closure test, jets used for the response matrix), %d centrality bin", cent), binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
        fHistManager.CreateTH2(Form("hRgDetLevelClosureResp_%d", cent), Form("Rg response at detector level (closure test, jets used for the response matrix), %d centrality bin", cent), binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
        fHistManager.CreateTH2(Form("hNsdPartLevelClosureNoResp_%d", cent), Form("Nsd response at particle level (closure test, jets not used for the response matrix), %d centrality bin", cent), binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
        fHistManager.CreateTH2(Form("hNsdDetLevelClosureNoResp_%d", cent), Form("Nsd response at detector level (closure test, jets not used for the response matrix), %d centrality bin", cent), binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
        fHistManager.CreateTH2(Form("hNsdPartLevelClosureResp_%d", cent), Form("Nsd response at particle level (closure test, jets used for the response matrix), %d centrality bin", cent), binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
        fHistManager.CreateTH2(Form("hNsdDetLevelClosureResp_%d", cent), Form("Nsd response at detector level (closure test, jets used for the response matrix), %d centrality bin", cent), binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
        fHistManager.CreateTH2(Form("hThetagPartLevelClosureNoResp_%d", cent), Form("Thetag response at particle level (closure test, jets not used for the response matrix), %d centrality bin", cent), binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
        fHistManager.CreateTH2(Form("hThetagDetLevelClosureNoResp_%d", cent), Form("Thetag response at detector level (closure test, jets not used for the response matrix), %d centrality bin", cent), binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
        fHistManager.CreateTH2(Form("hThetagPartLevelClosureResp_%d", cent), Form("Thetag response at particle level (closure test, jets used for the response matrix), %d centrality bin", cent), binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
        fHistManager.CreateTH2(Form("hThetagDetLevelClosureResp_%d", cent), Form("Thetag response at detector level (closure test, jets used for the response matrix), %d centrality bin", cent), binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());

        RooUnfoldResponse *r_zg = new RooUnfoldResponse(Form("hZgResponse_%d", cent), Form("z_{g} response matrix, %d centrality bin", cent)),
                          *r_rg = new RooUnfoldResponse(Form("hRgResponse_%d", cent), Form("r_{g} response matrix, %d centrality bin", cent)),
                          *r_nsd = new RooUnfoldResponse(Form("hNsdResponse_%d", cent), Form("n_{SD} response matrix, %d centrality bin", cent)),
                          *r_thetag = new RooUnfoldResponse(Form("hThetagResponse_%d", cent), Form("#Theta_{g} response matrix, %d centrality bin", cent)),
                          *r_zg_closure = new RooUnfoldResponse(Form("hZgResponseClosure_%d", cent), Form("z_{g} response matrix for the closure test, %d centrality bin", cent)),
                          *r_rg_closure = new RooUnfoldResponse(Form("hRgResponseClosure_%d", cent), Form("r_{g} response matrix for the closure test, %d centrality bin", cent)),
                          *r_nsd_closure = new RooUnfoldResponse(Form("hNsdResponseClosure_%d", cent), Form("n_{SD} response matrix for the closure test, %d centrality bin", cent)),
                          *r_thetag_closure = new RooUnfoldResponse(Form("hThetagResponseClosure_%d", cent), Form("#Theta_{g} response matrix for the closure test, %d centrality bin", cent));
        TString nameZgDetLevel = TString::Format("hZgDetLevel_%d", cent);
        TString nameZgPartLevel = TString::Format("hZgPartLevel_%d", cent);
        r_zg->Setup((TH1 *)fHistManager.FindObject(nameZgDetLevel), (TH1 *)fHistManager.FindObject(nameZgPartLevel));
        r_zg_closure->Setup((TH1 *)fHistManager.FindObject(nameZgDetLevel), (TH1 *)fHistManager.FindObject(nameZgPartLevel));
        TString nameRgDetLevel = TString::Format("hRgDetLevel_%d", cent);
        TString nameRgPartLevel = TString::Format("hRgPartLevel_%d", cent);
        r_rg->Setup((TH1 *)fHistManager.FindObject(nameRgDetLevel), (TH1 *)fHistManager.FindObject(nameRgPartLevel));
        r_rg_closure->Setup((TH1 *)fHistManager.FindObject(nameRgDetLevel), (TH1 *)fHistManager.FindObject(nameRgPartLevel));
        TString nameNsdDetLevel = TString::Format("hNsdDetLevel_%d", cent);
        TString nameNsdPartLevel = TString::Format("hNsdPartLevel_%d", cent);
        r_nsd->Setup((TH1 *)fHistManager.FindObject(nameNsdDetLevel), (TH1 *)fHistManager.FindObject(nameNsdPartLevel));
        r_nsd_closure->Setup((TH1 *)fHistManager.FindObject(nameNsdDetLevel), (TH1 *)fHistManager.FindObject(nameNsdPartLevel));
        TString nameThetagDetLevel = TString::Format("hThetagDetLevel_%d", cent);
        TString nameThetagPartLevel = TString::Format("hThetagPartLevel_%d", cent);
        r_thetag->Setup((TH1 *)fHistManager.FindObject(nameThetagDetLevel), (TH1 *)fHistManager.FindObject(nameThetagPartLevel));
        r_thetag_closure->Setup((TH1 *)fHistManager.FindObject(nameThetagDetLevel), (TH1 *)fHistManager.FindObject(nameThetagPartLevel));
        fZgResponse.push_back(r_zg);
        fZgResponseClosure.push_back(r_zg_closure);
        fRgResponse.push_back(r_rg);
        fRgResponseClosure.push_back(r_rg_closure);
        fNsdResponse.push_back(r_nsd);
        fNsdResponseClosure.push_back(r_nsd_closure);
        fThetagResponse.push_back(r_thetag);
        fThetagResponseClosure.push_back(r_thetag_closure);
      }
      if(fHasResponseMatrixSparse) {
        fHistManager.CreateTHnSparse(Form("hZgResponseSparse_%d", cent), Form("z_{g} response matrix, %d centrality bin", cent), 4, sparsebinningZg);
        fHistManager.CreateTHnSparse(Form("hRgResponseSparse_%d", cent), Form("z_{g} response matrix, %d centrality bin", cent), 4, sparsebinningRg);
        fHistManager.CreateTHnSparse(Form("hNsdResponseSparse_%d", cent), Form("z_{g} response matrix, %d centrality bin", cent), 4, sparsebinningNsd);
        fHistManager.CreateTHnSparse(Form("hThetagResponseSparse_%d", cent), Form("z_{g} response matrix, %d centrality bin", cent), 4, sparsebinningThetag);
        fHistManager.CreateTHnSparse(Form("hZgResponseClosureSparse_%d", cent), Form("z_{g} response matrix for closure test, %d centrality bin", cent), 4, sparsebinningZg);
        fHistManager.CreateTHnSparse(Form("hRgResponseClosureSparse_%d", cent), Form("z_{g} response matrix for closure test, %d centrality bin", cent), 4, sparsebinningRg);
        fHistManager.CreateTHnSparse(Form("hNsdResponseClosureSparse_%d", cent), Form("z_{g} response matrix for closure test, %d centrality bin", cent), 4, sparsebinningNsd);
        fHistManager.CreateTHnSparse(Form("hThetagResponseClosureSparse_%d", cent), Form("z_{g} response matrix for closure test, %d centrality bin", cent), 4, sparsebinningThetag);

        fHistManager.CreateTHnSparse(Form("hZgResponseClosureNoRespSparse_%d", cent), Form("z_{g} response matrix for closure test, %d centrality bin", cent), 4, sparsebinningZg);
        fHistManager.CreateTHnSparse(Form("hRgResponseClosureNoRespSparse_%d", cent), Form("z_{g} response matrix for closure test, %d centrality bin", cent), 4, sparsebinningRg);
        fHistManager.CreateTHnSparse(Form("hNsdResponseClosureNoRespSparse_%d", cent), Form("z_{g} response matrix for closure test, %d centrality bin", cent), 4, sparsebinningNsd);
        fHistManager.CreateTHnSparse(Form("hThetagResponseClosureNoRespSparse_%d", cent), Form("z_{g} response matrix for closure test, %d centrality bin", cent), 4, sparsebinningThetag);
        fHistManager.CreateTH2(Form("hZgPartLevelClosureNoRespFine_%d", cent), Form("Zg response at particle level (closure test, jets not used for the response matrix), %d centrality bin", cent), binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtFine.GetSize() - 1, binEdgesPtFine.GetArray());
        fHistManager.CreateTH2(Form("hZgDetLevelClosureNoRespFine_%d", cent), Form("Zg response at detector level (closure test, jets not used for the response matrix), %d centrality bin", cent), binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtFine.GetSize() - 1, binEdgesPtFine.GetArray());
        fHistManager.CreateTH2(Form("hRgPartLevelClosureNoRespFine_%d", cent), Form("Rg response at particle level (closure test, jets not used for the response matrix), %d centrality bin", cent), binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtFine.GetSize() - 1, binEdgesPtFine.GetArray());
        fHistManager.CreateTH2(Form("hRgDetLevelClosureNoRespFine_%d", cent), Form("Rg response at detector level (closure test, jets not used for the response matrix), %d centrality bin", cent), binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtFine.GetSize() - 1, binEdgesPtFine.GetArray());
        fHistManager.CreateTH2(Form("hNsdPartLevelClosureNoRespFine_%d", cent), Form("Nsd response at particle level (closure test, jets not used for the response matrix), %d centrality bin", cent), binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtFine.GetSize() - 1, binEdgesPtFine.GetArray());
        fHistManager.CreateTH2(Form("hNsdDetLevelClosureNoRespFine_%d", cent), Form("Nsd response at detector level (closure test, jets not used for the response matrix), %d centrality bin", cent), binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtFine.GetSize() - 1, binEdgesPtFine.GetArray());
        fHistManager.CreateTH2(Form("hThetagPartLevelClosureNoRespFine_%d", cent), Form("Thetag response at particle level (closure test, jets not used for the response matrix), %d centrality bin", cent), binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtFine.GetSize() - 1, binEdgesPtFine.GetArray());
        fHistManager.CreateTH2(Form("hThetagDetLevelClosureNoRespFine_%d", cent), Form("Thetag response at detector level (closure test, jets not used for the response matrix), %d centrality bin", cent), binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtFine.GetSize() - 1, binEdgesPtFine.GetArray());
      }

      // Residuals vs. pt,part
      fHistManager.CreateTH2(Form("hZgResiduals_%d", cent), Form("z_{g} residuals (%d centrality bin); p_{t,part} (GeV/c); z_{g, det} - z_{g, part}", cent), 350, 0., 350., 100, -1., 1.);
      fHistManager.CreateTH2(Form("hRgResiduals_%d", cent), Form("R_{g} residuals (%d centrality bin); p_{t,part} (GeV/c); R_{g, det} - R_{g, part}", cent), 350, 0., 350., 100, -1., 1.);
      fHistManager.CreateTH2(Form("hThetagResiduals_%d", cent), Form("#Theta_{g} residuals (%d centrality bin); p_{t,part} (GeV/c); #Theta_{g, det} - #Theta_{g, part}", cent), 350, 0., 350., 200, -2., 2.);
      fHistManager.CreateTH2(Form("hNsdResiduals_%d", cent), Form("n_{SD} residuals (%d centrality bin); p_{t,part} (GeV/c); n_{SD, det} - n_{SD, part}", cent), 350, 0., 350., 80, -40., 40.);
      fHistManager.CreateTH2(Form("hZgResidualsNormalized_%d", cent), Form("z_{g} residuals (normalized, %d centrality bin); p_{t,part} (GeV/c); (z_{g, det} - z_{g, part})/z_{g, part}", cent), 350, 0., 350., 100, -1., 1.);
      fHistManager.CreateTH2(Form("hRgResidualsNormalized_%d", cent), Form("R_{g} residuals (normalized, %d centrality bin); p_{t,part} (GeV/c); (R_{g, det} - R_{g, part})/R_{g, part}", cent), 350, 0., 350., 100, -1., 1.);
      fHistManager.CreateTH2(Form("hThetagResidualsNormalized_%d", cent), Form("#Theta_{g} residuals (normalized, %d centrality bin); p_{t,part} (GeV/c); (#Theta_{g, det} - #Theta_{g, part})/#Theta_{g, part}", cent), 350, 0., 350., 100, -1., 1.);
      fHistManager.CreateTH2(Form("hNsdResidualsNormalized_%d", cent), Form("n_{SD} residuals (normalized, %d centrality bin); p_{t,part} (GeV/c); (n_{SD, det} - n_{SD, part})/n_{SD, part}", cent), 350, 0., 350., 100, -10., 10.);
      // Residuals vs. Rg
      fHistManager.CreateTH2(Form("hZgResidualsRg_%d", cent), Form("z_{g} residuals vs. R_{g} (%d centrality bin); R_{g}; z_{g,det} - z_{g,part}", cent), 100, 0., 1., 100, -1., 1.);
      fHistManager.CreateTH2(Form("hRgResidualsRg_%d", cent), Form("R_{g} residuals vs. R_{g} (%d centrality bin); R_{g}; R_{g,det} - R_{g,part}", cent), 100, 0., 1., 100, -1., 1.);
      fHistManager.CreateTH2(Form("hThetagResidualsRg_%d", cent), Form("#Theta_{g} residuals vs. R_{g} (%d centrality bin); R_{g}; #Theta_{g,det} - #Theta_{g,part}", cent), 100, 0., 1., 100, -1., 1.);
      fHistManager.CreateTH2(Form("hNsdResidualsRg_%d", cent), Form("z_{g} residuals vs. R_{g} (%d centrality bin); R_{g}; z_{g,det} - z_{g,part}", cent), 100, 0., 1., 80, -40., 40.);
      fHistManager.CreateTH2(Form("hZgResidualsRgNormalized_%d", cent), Form("z_{g} residuals (normalized) vs. R_{g} (%d centrality bin); R_{g}; z_{g,det} - z_{g,part}", cent), 100, 0., 1., 100, -1., 1.);
      fHistManager.CreateTH2(Form("hRgResidualsRgNormalized_%d", cent), Form("R_{g} residuals (normalized) vs. R_{g} (%d centrality bin); R_{g}; R_{g,det} - R_{g,part}", cent), 100, 0., 1., 100, -1., 1.);
      fHistManager.CreateTH2(Form("hThetagResidualsRgNormalized_%d", cent), Form("#Theta_{g} residuals (normalized) vs. R_{g} (%d centrality bin); R_{g}; #Theta_{g,det} - #Theta_{g,part}", cent), 100, 0., 1., 100, -1., 1.);
      fHistManager.CreateTH2(Form("hNsdResidualsRgNormalized_%d", cent), Form("z_{g} residuals (normalized) vs. R_{g} (%d centrality bin); R_{g}; z_{g,det} - z_{g,part}", cent), 100, 0., 1., 100, -10., 10.); 
      // Residuals vs. Thetag
      fHistManager.CreateTH2(Form("hZgResidualsThetag_%d", cent), Form("z_{g} residuals vs. #Theta_{g} (%d centrality bin); #Theta_{g}; z_{g,det} - z_{g,part}", cent), 100, 0., 1., 100, -1., 1.);
      fHistManager.CreateTH2(Form("hRgResidualsThetag_%d", cent), Form("R_{g} residuals vs. #Theta_{g} (%d centrality bin); #Theta_{g}; R_{g,det} - R_{g,part}", cent), 100, 0., 1., 100, -1., 1.);
      fHistManager.CreateTH2(Form("hThetagResidualsThetag_%d", cent), Form("#Theta_{g} residuals vs. #Theta_{g} (%d centrality bin); #Theta_{g}; #Theta_{g,det} - #Theta_{g,part}", cent), 100, 0., 1., 100, -1., 1.);
      fHistManager.CreateTH2(Form("hNsdResidualsThetag_%d", cent), Form("z_{g} residuals vs. #Theta_{g} (%d centrality bin); #Theta_{g}; z_{g,det} - z_{g,part}", cent), 100, 0., 1., 80, -40., 40.);
      fHistManager.CreateTH2(Form("hZgResidualsThetagNormalized_%d", cent), Form("z_{g} residuals (normalized) vs. #Theta_{g} (%d centrality bin); #Theta_{g}; z_{g,det} - z_{g,part}", cent), 100, 0., 1., 100, -1., 1.);
      fHistManager.CreateTH2(Form("hRgResidualsThetagNormalized_%d", cent), Form("R_{g} residuals (normalized) vs. #Theta_{g} (%d centrality bin); #Theta_{g}; R_{g,det} - R_{g,part}", cent), 100, 0., 1., 100, -1., 1.);
      fHistManager.CreateTH2(Form("hThetagResidualsThetagNormalized_%d", cent), Form("#Theta_{g} residuals (normalized) vs. #Theta_{g} (%d centrality bin); #Theta_{g}; #Theta_{g,det} - #Theta_{g,part}", cent), 100, 0., 1., 100, -1., 1.);
      fHistManager.CreateTH2(Form("hNsdResidualsThetagNormalized_%d", cent), Form("z_{g} residuals (normalized) vs. #Theta_{g} (%d centrality bin); #Theta_{g}; z_{g,det} - z_{g,part}", cent), 100, 0., 1., 100, -10., 10.); 

      // a bit of QA stuff
      fHistManager.CreateTH1(Form("hSkippedJetsPart_%d", cent), Form("Skipped jets at part. level, %d centrality bin", cent), 350, 0., 350.);
      fHistManager.CreateTH1(Form("hSkippedJetsDet_%d", cent), Form("Skipped jets at det. level, %d centrality bin", cent), 350, 0., 350.);
      fHistManager.CreateTH2(Form("hQAEtaPhiPart_%d", cent), Form("#eta vs. #phi for selected part. level jets (%d centrality bin); #eta; #phi", cent), 100, -1., 1., 100, 0., 7.);
      fHistManager.CreateTH2(Form("hQAEtaPhiDet_%d", cent), Form("#eta vs. #phi for selected det. level jets (%d centrality bin); #eta; #phi", cent), 100, -1., 1., 100, 0., 7.);
      fHistManager.CreateTH2(Form("hQANEFPtPart_%d", cent), Form("Neutral energy fraction at part. level, %d centrality bin; p_{t} (GeV/c); NEF", cent), 350, 0., 350., 100, 0., 1.);
      fHistManager.CreateTH2(Form("hQANEFPtDet_%d", cent), Form("Neutral energy fraction at det. level, %d centrality bin; p_{t} (GeV/c); NEF", cent), 350, 0., 350., 100, 0., 1.);
      fHistManager.CreateTH2(Form("hQAZchPtPart_%d", cent), "z_{ch,max} at part. level; p_{t} (GeV/c); z_{ch,max}", 350, 0., 350., 100, 0., 1.);
      fHistManager.CreateTH2(Form("hQAZchPtDet_%d", cent), Form("z_{ch,max} at det. level, %d centrality bin; p_{t} (GeV/c); z_{ch,max}", cent), 350, 0., 350., 100, 0., 1.);
      fHistManager.CreateTH2(Form("hQAZnePtPart_%d", cent), Form("z_{ne,max} at part. level, %d centrality bin; p_{t} (GeV/c); z_{ne,max}", cent), 350, 0., 350., 100, 0., 1.);
      fHistManager.CreateTH2(Form("hQAZnePtDet_%d", cent), Form("z_{ne,max} at det. level, %d centrality bin; p_{t} (GeV/c); z_{ne,max}", cent), 350, 0., 350., 100, 0., 1.);
      fHistManager.CreateTH2(Form("hQANChPtPart_%d", cent), Form("Number of charged constituents at part. level, %d centrality bin; p_{t} (GeV/c); N_{ch}", cent), 350, 0., 350., 100, 0., 100.);
      fHistManager.CreateTH2(Form("hQANChPtDet_%d", cent), Form("Number of charged constituents at det. level, %d centrality bin; p_{t} (GeV/c); N_{ch}", cent), 350, 0., 350., 100, 0., 100.);
      fHistManager.CreateTH2(Form("hQANnePtPart_%d", cent), Form("Number of neutral constituents at part. level, %d centrality bin; p_{t} (GeV/c); N_{ne}", cent), 350, 0., 350., 100, 0., 100.);
      fHistManager.CreateTH2(Form("hQANnePtDet_%d", cent), Form("Number of neutral constituents at det. level, %d centrality bin; p_{t} (GeV/c); N_{ne}", cent), 350, 0., 350., 100, 0., 100.);
      fHistManager.CreateTH2(Form("hQAJetAreaVsJetPtPart_%d", cent), Form("Jet area vs. jet pt at particle level (%d centrality bin); p_{t} (GeV/c); Area", cent), 350, 0., 350., 200, 0., 2.);
      fHistManager.CreateTH2(Form("hQAJetAreaVsJetPtDet_%d", cent), Form("Jet area vs. jet pt at detector level (%d centrality bin); p_{t} (GeV/c); Area", cent), 350, 0., 350., 200, 0., 2.);
      fHistManager.CreateTH2(Form("hQAJetAreaVsNEFPart_%d", cent), Form("Jet area vs. NEF at particle level (%d centrality bin); NEF; Area", cent), 100, 0., 1., 200, 0., 2.);
      fHistManager.CreateTH2(Form("hQAJetAreaVsNEFDet_%d", cent), Form("Jet area vs. NEF at detector level (%d centrality bin); NEF; Area", cent), 100, 0., 1., 200, 0., 2.);
      fHistManager.CreateTH2(Form("hQAJetAreaVsNConstPart_%d", cent), Form("Jet area vs. number of consituents at particle level (%d centrality bin); Number of constituents; Area", cent), 101, -0.5, 100.5, 200, 0., 2.);
      fHistManager.CreateTH2(Form("hQAJetAreaVsNConstDet_%d", cent), Form("Jet area vs. number of consituents at detector level (%d centrality bin); Number of constituents; Area", cent), 101, -0.5, 100.5, 200, 0., 2.);
      fHistManager.CreateTH1(Form("hQAMatchingDRAbs_%d", cent), Form("Distance between part. level jet and  det. level jet (%d centrality bin)", cent), 100, 0., 1.);
      fHistManager.CreateTH1(Form("hQAMatchingDRAbs_%d", cent), Form("Distance between part. level jet and  det. level jet (%d centrality bin)", cent), 100, 0., 1.);
    }
  }
  else
  {
    if(fHasResponseMatrixRooUnfold){
      fHistManager.CreateTH2("hZgDetLevel", "Zg response at detector level", binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
      fHistManager.CreateTH2("hZgPartLevel", "Zg response at particle level", binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
      fHistManager.CreateTH2("hZgPartLevelTruncated", "Zg response at particle level (truncated)", binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
      fHistManager.CreateTH2("hRgDetLevel", "Rg response at detector level", binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
      fHistManager.CreateTH2("hRgPartLevel", "Rg response at particle level", binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
      fHistManager.CreateTH2("hRgPartLevelTruncated", "Rg response at particle level (truncated)", binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
      fHistManager.CreateTH2("hNsdDetLevel", "Nsd response at detector level", binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
      fHistManager.CreateTH2("hNsdPartLevel", "Nsd response at particle level", binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
      fHistManager.CreateTH2("hNsdPartLevelTruncated", "Nsd response at particle level (truncated)", binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
      fHistManager.CreateTH2("hThetagDetLevel", "Thetag response at detector level", binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
      fHistManager.CreateTH2("hThetagPartLevel", "Thetag response at particle level", binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
      fHistManager.CreateTH2("hThetagPartLevelTruncated", "Thetag response at particle level (truncated)", binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());

      // For closure test
      fHistManager.CreateTH2("hZgPartLevelClosureNoResp", "Zg response at particle level (closure test, jets not used for the response matrix)", binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
      fHistManager.CreateTH2("hZgDetLevelClosureNoResp", "Zg response at detector level (closure test, jets not used for the response matrix)", binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
      fHistManager.CreateTH2("hZgPartLevelClosureResp", "Zg response at particle level (closure test, jets used for the response matrix)", binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
      fHistManager.CreateTH2("hZgDetLevelClosureResp", "Zg response at detector level (closure test, jets used for the response matrix)", binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
      fHistManager.CreateTH2("hRgPartLevelClosureNoResp", "Rg response at particle level (closure test, jets not used for the response matrix)", binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
      fHistManager.CreateTH2("hRgDetLevelClosureNoResp", "Rg response at detector level (closure test, jets not used for the response matrix)", binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
      fHistManager.CreateTH2("hRgPartLevelClosureResp", "Rg response at particle level (closure test, jets used for the response matrix)", binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
      fHistManager.CreateTH2("hRgDetLevelClosureResp", "Rg response at detector level (closure test, jets used for the response matrix)", binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
      fHistManager.CreateTH2("hNsdPartLevelClosureNoResp", "Nsd response at particle level (closure test, jets not used for the response matrix)", binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
      fHistManager.CreateTH2("hNsdDetLevelClosureNoResp", "Nsd response at detector level (closure test, jets not used for the response matrix)", binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
      fHistManager.CreateTH2("hNsdPartLevelClosureResp", "Nsd response at particle level (closure test, jets used for the response matrix)", binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
      fHistManager.CreateTH2("hNsdDetLevelClosureResp", "Nsd response at detector level (closure test, jets used for the response matrix)", binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
      fHistManager.CreateTH2("hThetagPartLevelClosureNoResp", "Thetag response at particle level (closure test, jets not used for the response matrix)", binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
      fHistManager.CreateTH2("hThetagDetLevelClosureNoResp", "Thetag response at detector level (closure test, jets not used for the response matrix)", binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());
      fHistManager.CreateTH2("hThetagPartLevelClosureResp", "Thetag response at particle level (closure test, jets used for the response matrix)", binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtPart.GetSize() - 1, binEdgesPtPart.GetArray());
      fHistManager.CreateTH2("hThetagDetLevelClosureResp", "Thetag response at detector level (closure test, jets used for the response matrix)", binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtDet.GetSize() - 1, binEdgesPtDet.GetArray());

      RooUnfoldResponse *r_zg = new RooUnfoldResponse("hZgResponse", "z_{g} response matrix"),
                        *r_rg = new RooUnfoldResponse("hRgResponse", "r_{g} response matrix"),
                        *r_nsd = new RooUnfoldResponse("hNsdResponse", "n_{SD} response matrix"),
                        *r_thetag = new RooUnfoldResponse("hThetagResponse", "#Theta_{g} response matrix");
      r_zg->Setup((TH1 *)fHistManager.FindObject("hZgDetLevel"), (TH1 *)fHistManager.FindObject("hZgPartLevel"));
      r_rg->Setup((TH1 *)fHistManager.FindObject("hRgDetLevel"), (TH1 *)fHistManager.FindObject("hRgPartLevel"));
      r_nsd->Setup((TH1 *)fHistManager.FindObject("hNsdDetLevel"), (TH1 *)fHistManager.FindObject("hNsdPartLevel"));
      r_thetag->Setup((TH1 *)fHistManager.FindObject("hThetagDetLevel"), (TH1 *)fHistManager.FindObject("hThetagPartLevel"));
      fZgResponse.push_back(r_zg);
      fRgResponse.push_back(r_rg);
      fNsdResponse.push_back(r_nsd);
      fThetagResponse.push_back(r_thetag);
      // Create RooUnfold response from THnSparse
      RooUnfoldResponse *r_zg_closure = new RooUnfoldResponse("hZgResponseClosure", "z_{g} response matrix for the closure test"),
                        *r_rg_closure = new RooUnfoldResponse("hRgResponseClosure", "z_{g} response matrix for the closure test"),
                        *r_nsd_closure = new RooUnfoldResponse("hNsdResponseClosure", "z_{g} response matrix for the closure test"),
                        *r_thetag_closure = new RooUnfoldResponse("hThetagResponseClosure", "#Theta_{g} response matrix for the closure test");
      r_zg_closure->Setup((TH1 *)fHistManager.FindObject("hZgDetLevel"), (TH1 *)fHistManager.FindObject("hZgPartLevel"));
      r_rg_closure->Setup((TH1 *)fHistManager.FindObject("hRgDetLevel"), (TH1 *)fHistManager.FindObject("hRgPartLevel"));
      r_nsd_closure->Setup((TH1 *)fHistManager.FindObject("hNsdDetLevel"), (TH1 *)fHistManager.FindObject("hNsdPartLevel"));
      r_thetag_closure->Setup((TH1 *)fHistManager.FindObject("hThetagDetLevel"), (TH1 *)fHistManager.FindObject("hThetagPartLevel"));
      fZgResponseClosure.push_back(r_zg_closure);
      fRgResponseClosure.push_back(r_rg_closure);
      fNsdResponseClosure.push_back(r_nsd_closure);
      fThetagResponseClosure.push_back(r_thetag_closure);
    }

    if(fHasResponseMatrixSparse) {
      fHistManager.CreateTHnSparse("hZgResponseSparse", "z_{g} response matrix", 4, sparsebinningZg);
      fHistManager.CreateTHnSparse("hRgResponseSparse", "z_{g} response matrix", 4, sparsebinningRg);
      fHistManager.CreateTHnSparse("hNsdResponseSparse", "z_{g} response matrix", 4, sparsebinningNsd);
      fHistManager.CreateTHnSparse("hThetagResponseSparse", "z_{g} response matrix", 4, sparsebinningThetag);
      fHistManager.CreateTHnSparse("hZgResponseClosureSparse", "z_{g} response matrix for closure test", 4, sparsebinningZg);
      fHistManager.CreateTHnSparse("hRgResponseClosureSparse", "z_{g} response matrix for closure test", 4, sparsebinningRg);
      fHistManager.CreateTHnSparse("hNsdResponseClosureSparse", "z_{g} response matrix for closure test", 4, sparsebinningNsd);
      fHistManager.CreateTHnSparse("hThetagResponseClosureSparse", "z_{g} response matrix for closure test", 4, sparsebinningThetag);

      fHistManager.CreateTHnSparse("hZgResponseClosureNoRespSparse", "z_{g} response matrix for closure test pseudo data", 4, sparsebinningZg);
      fHistManager.CreateTHnSparse("hRgResponseClosureNoRespSparse", "z_{g} response matrix for closure test pseudo data", 4, sparsebinningRg);
      fHistManager.CreateTHnSparse("hNsdResponseClosureNoRespSparse", "z_{g} response matrix for closure test pseudo data", 4, sparsebinningNsd);
      fHistManager.CreateTHnSparse("hThetagResponseClosureNoRespSparse", "z_{g} response matrix for closure test pseudo data", 4, sparsebinningThetag);
      fHistManager.CreateTH2("hZgPartLevelClosureNoRespFine", "Zg response at particle level (closure test, jets not used for the response matrix)", binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtFine.GetSize() - 1, binEdgesPtFine.GetArray());
      fHistManager.CreateTH2("hZgDetLevelClosureNoRespFine", "Zg response at detector level (closure test, jets not used for the response matrix)", binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtFine.GetSize() - 1, binEdgesPtFine.GetArray());
      fHistManager.CreateTH2("hRgPartLevelClosureNoRespFine", "Rg response at particle level (closure test, jets not used for the response matrix)", binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtFine.GetSize() - 1, binEdgesPtFine.GetArray());
      fHistManager.CreateTH2("hRgDetLevelClosureNoRespFine", "Rg response at detector level (closure test, jets not used for the response matrix)", binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtFine.GetSize() - 1, binEdgesPtFine.GetArray());
      fHistManager.CreateTH2("hNsdPartLevelClosureNoRespFine", "Nsd response at particle level (closure test, jets not used for the response matrix)", binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtFine.GetSize() - 1, binEdgesPtFine.GetArray());
      fHistManager.CreateTH2("hNsdDetLevelClosureNoRespFine", "Nsd response at detector level (closure test, jets not used for the response matrix)", binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtFine.GetSize() - 1, binEdgesPtFine.GetArray());
      fHistManager.CreateTH2("hThetagPartLevelClosureNoRespFine", "Thetag response at particle level (closure test, jets not used for the response matrix)", binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtFine.GetSize() - 1, binEdgesPtFine.GetArray());
      fHistManager.CreateTH2("hThetagDetLevelClosureNoRespFine", "Thetag response at detector level (closure test, jets not used for the response matrix)", binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtFine.GetSize() - 1, binEdgesPtFine.GetArray());
    }

    // Residuals vs. pt,part
    fHistManager.CreateTH2("hZgResiduals", "z_{g} residuals vs. p_{t,part}; p_{t,part} (GeV/c); z_{g, det} - z_{g, part}", 350, 0., 350., 100, -1., 1.);
    fHistManager.CreateTH2("hRgResiduals", "R_{g} residuals vs. p_{t,part}; p_{t,part} (GeV/c); R_{g, det} - R_{g, part}", 350, 0., 350., 100, -1., 1.);
    fHistManager.CreateTH2("hThetagResiduals", "#Theta_{g} residuals vs. p_{t,part}; p_{t,part} (GeV/c); #Theta_{g, det} - #Theta_{g, part}", 350, 0., 350., 200, -2., 2.);
    fHistManager.CreateTH2("hNsdResiduals", "n_{SD} residuals vs. p_{t,part}; p_{t,part} (GeV/c); n_{SD, det} - n_{SD, part}", 350, 0., 350., 80, -40., 40.);
    fHistManager.CreateTH2("hZgResidualsNormalized", "z_{g} residuals (normalized) vs. p_{t,part}; p_{t,part} (GeV/c); (z_{g, det} - z_{g, part})/z_{g, part}", 350, 0., 350., 100, -1., 1.);
    fHistManager.CreateTH2("hRgResidualsNormalized", "R_{g} residuals (normalized) vs. p_{t,part}; p_{t,part} (GeV/c); (R_{g, det} - R_{g, part})/R_{g, part}", 350, 0., 350., 100, -1., 1.);
    fHistManager.CreateTH2("hThetagResidualsNormalized", "#Theta_{g} residuals (normalized) vs. p_{t,part}; p_{t,part} (GeV/c); (#Theta_{g, det} - #Theta_{g, part})/#Theta_{g, part}", 350, 0., 350., 100, -1., 1.);
    fHistManager.CreateTH2("hNsdResidualsNormalized", "n_{SD} residuals (normalized) vs. p_{t,part}; p_{t,part} (GeV/c); (n_{SD, det} - n_{SD, part})/n_{SD, part}", 350, 0., 350., 100, -10., 10.);
    // Residuals vs. Rg
    fHistManager.CreateTH2("hZgResidualsRg", "z_{g} residuals vs. R_{g}; R_{g}; z_{g,det} - z_{g,part}", 100, 0., 1., 100, -1., 1.);
    fHistManager.CreateTH2("hRgResidualsRg", "R_{g} residuals vs. R_{g}; R_{g}; R_{g,det} - R_{g,part}", 100, 0., 1., 100, -1., 1.);
    fHistManager.CreateTH2("hThetagResidualsRg", "#Theta_{g} residuals vs. R_{g}; R_{g}; #Theta_{g,det} - #Theta_{g,part}", 100, 0., 1., 200, -2., 2.);
    fHistManager.CreateTH2("hNsdResidualsRg", "z_{g} residuals vs. R_{g}; R_{g}; z_{g,det} - z_{g,part}", 100, 0., 1., 80, -40., 40.);
    fHistManager.CreateTH2("hZgResidualsRgNormalized", "z_{g} residuals (normalized) vs. R_{g}; R_{g}; z_{g,det} - z_{g,part}", 100, 0., 1., 100, -1., 1.);
    fHistManager.CreateTH2("hRgResidualsRgNormalized", "R_{g} residuals (normalized) vs. R_{g}; R_{g}; R_{g,det} - R_{g,part}", 100, 0., 1., 100, -1., 1.);
    fHistManager.CreateTH2("hThetagResidualsRgNormalized", "#Theta_{g} residuals (normalized) vs. R_{g}; R_{g}; #Theta_{g,det} - #Theta_{g,part}", 100, 0., 1., 100, -1., 1.);
    fHistManager.CreateTH2("hNsdResidualsRgNormalized", "z_{g} residuals (normalized) vs. R_{g}; R_{g}; z_{g,det} - z_{g,part}", 100, 0., 1., 100, -10., 10.);
    // Residuals vs. Thetag
    fHistManager.CreateTH2("hZgResidualsThetag", "z_{g} residuals vs. #Theta_{g}; #Theta_{g}; z_{g,det} - z_{g,part}", 100, 0., 1., 100, -1., 1.);
    fHistManager.CreateTH2("hRgResidualsThetag", "R_{g} residuals vs. #Theta_{g}; #Theta_{g}; R_{g,det} - R_{g,part}", 100, 0., 1., 100, -1., 1.);
    fHistManager.CreateTH2("hThetagResidualsThetag", "#Theta_{g} residuals vs. #Theta_{g}; #Theta_{g}; #Theta_{g,det} - #Theta_{g,part}", 100, 0., 1., 200, -2., 2.);
    fHistManager.CreateTH2("hNsdResidualsThetag", "z_{g} residuals vs. #Theta_{g}; #Theta_{g}; z_{g,det} - z_{g,part}", 100, 0., 1., 80, -40., 40.);
    fHistManager.CreateTH2("hZgResidualsThetagNormalized", "z_{g} residuals (normalized) vs. #Theta_{g}; #Theta_{g}; z_{g,det} - z_{g,part}", 100, 0., 1., 100, -1., 1.);
    fHistManager.CreateTH2("hRgResidualsThetagNormalized", "R_{g} residuals (normalized) vs. #Theta_{g}; #Theta_{g}; R_{g,det} - R_{g,part}", 100, 0., 1., 100, -1., 1.);
    fHistManager.CreateTH2("hThetagResidualsThetagNormalized", "#Theta_{g} residuals (normalized) vs. #Theta_{g}; #Theta_{g}; #Theta_{g,det} - #Theta_{g,part}", 100, -1., 1., 100, -1., 1.);
    fHistManager.CreateTH2("hNsdResidualsThetagNormalized", "z_{g} residuals (normalized) vs. #Theta_{g}; #Theta_{g}; z_{g,det} - z_{g,part}", 100, 0., 1., 100, -10., 10.);

    // a bit of QA stuff
    fHistManager.CreateTH1("hSkippedJetsPart", "Skipped jets at part. level", 350, 0., 350.);
    fHistManager.CreateTH1("hSkippedJetsDet", "Skipped jets at det. level", 350, 0., 350.);
    fHistManager.CreateTH2("hQAEtaPhiPart", "#eta vs. #phi for selected part. level jets; #eta; #phi", 100, -1., 1., 100, 0., 7.);
    fHistManager.CreateTH2("hQAEtaPhiDet", "#eta vs. #phi for selected det. level jets; #eta; #phi", 100, -1., 1., 100, 0., 7.);
    fHistManager.CreateTH2("hQANEFPtPart", "Neutral energy fraction at part. level; p_{t} (GeV/c); NEF", 350, 0., 350., 100, 0., 1.);
    fHistManager.CreateTH2("hQANEFPtDet", "Neutral energy fraction at det. level; p_{t} (GeV/c); NEF", 350, 0., 350., 100, 0., 1.);
    fHistManager.CreateTH2("hQAZchPtPart", "z_{ch,max} at part. level; p_{t} (GeV/c); z_{ch,max}", 350, 0., 350., 100, 0., 1.);
    fHistManager.CreateTH2("hQAZchPtDet", "z_{ch,max} at det. level; p_{t} (GeV/c); z_{ch,max}", 350, 0., 350., 100, 0., 1.);
    fHistManager.CreateTH2("hQAZnePtPart", "z_{ne,max} at part. level; p_{t} (GeV/c); z_{ne,max}", 350, 0., 350., 100, 0., 1.);
    fHistManager.CreateTH2("hQAZnePtDet", "z_{ne,max} at det. level; p_{t} (GeV/c); z_{ne,max}", 350, 0., 350., 100, 0., 1.);
    fHistManager.CreateTH2("hQANChPtPart", "Number of charged constituents at part. level; p_{t} (GeV/c); N_{ch}", 350, 0., 350., 100, 0., 100.);
    fHistManager.CreateTH2("hQANChPtDet", "Number of charged constituents at det. level; p_{t} (GeV/c); N_{ch}", 350, 0., 350., 100, 0., 100.);
    fHistManager.CreateTH2("hQANnePtPart", "Number of neutral constituents at part. level; p_{t} (GeV/c); N_{ne}", 350, 0., 350., 100, 0., 100.);
    fHistManager.CreateTH2("hQANnePtDet", "Number of neutral constituents at det. level; p_{t} (GeV/c); N_{ne}", 350, 0., 350., 100, 0., 100.);
    fHistManager.CreateTH2("hQAJetAreaVsJetPtPart", "Jet area vs. jet pt at particle level; p_{t} (GeV/c); Area", 350, 0., 350., 200, 0., 2.);
    fHistManager.CreateTH2("hQAJetAreaVsJetPtDet", "Jet area vs. jet pt at detector level; p_{t} (GeV/c); Area", 350, 0., 350., 200, 0., 2.);
    fHistManager.CreateTH2("hQAJetAreaVsNEFPart", "Jet area vs. NEF at particle level; NEF; Area", 100, 0., 1., 200, 0., 2.);
    fHistManager.CreateTH2("hQAJetAreaVsNEFDet", "Jet area vs. NEF at detector level; NEF; Area", 100, 0., 1., 200, 0., 2.);
    fHistManager.CreateTH2("hQAJetAreaVsNConstPart", "Jet area vs. number of consituents at particle level; Number of constituents; Area", 101, -0.5, 100.5, 200, 0., 2.);
    fHistManager.CreateTH2("hQAJetAreaVsNConstDet", "Jet area vs. number of consituents at detector level; Number of constituents; Area", 101, -0.5, 100.5, 200, 0., 2.);
    fHistManager.CreateTH1("hQAMatchingDRAbs", "Distance between part. level jet and  det. level jet", 100, 0., 1.);
    fHistManager.CreateTH1("hQAMatchingDRel", "Distance between part. level jet and  det. level jet", 100, 0., 1.);
  }

  // a bit of QA stuff
  fHistManager.CreateTH2("hSDUsedChargedPtjvPtcPart", "p_{t,j} vs. p_{t,const} for tracks used in SD (part. level); p_{t,j} (GeV/c); p_{t,ch} (GeV/c)", 350., 0., 350., 350, 0., 350.);
  fHistManager.CreateTH2("hSDUsedChargedPtjvPtcDet", "p_{t,j} vs. p_{t,const} for tracks used in SD (det. level); p_{t,j} (GeV/c); p_{t,ch} (GeV/c)", 350., 0., 350., 350, 0., 350.);
  fHistManager.CreateTH2("hSDUsedNeutralPtjvPtcPart", "p_{t,j} vs. p_{t,const} for clusters used in SD (part. level); p_{t,j} (GeV/c); p_{t,ne} (GeV/c)", 350., 0., 350., 350, 0., 350.);
  fHistManager.CreateTH2("hSDUsedNeutralPtjvPtcDet", "p_{t,j} vs. p_{t,const} for clusters used in SD (det. level); p_{t,j} (GeV/c); p_{t,ne} (GeV/c)", 350., 0., 350., 350, 0., 350.);
  fHistManager.CreateTH2("hSDUsedChargedPtjvPtcMaxPart", "p_{t,j} vs. p_{t,const} for max tracks used in SD (part. level); p_{t,j} (GeV/c); p_{t,ch} (GeV/c)", 350, 0., 350., 350., 0., 350.);
  fHistManager.CreateTH2("hSDUsedChargedPtjvPtcMaxDet", "p_{t,j} vs. p_{t,const} for max tracks used in SD (det. level); p_{t,j} (GeV/c); p_{t,ch} (GeV/c)", 350, 0., 350., 350., 0., 350.);
  fHistManager.CreateTH2("hSDUsedNeutralPtjvPcMaxPart", "p_{t,j} vs. p_{t,const} for max clusters used in SD (part. level); p_{t,j} (GeV/c); p_{t,ne} (GeV/c)", 350, 0., 350., 350., 0., 350.);
  fHistManager.CreateTH2("hSDUsedNeutralPtjvPcMaxDet", "p_{t,j} vs. p_{t,const} for max clusters used in SD (det. level); p_{t,j} (GeV/c); p_{t,ne} (GeV/c)", 350, 0., 350., 350., 0., 350.);
  fHistManager.CreateTH2("hSDUsedChargedEtaPhiPart", "#eta-phi for tracks used in SD (part. level); #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistManager.CreateTH2("hSDUsedChargedEtaPhiDet", "#eta-phi for tracks used in SD (det. level); #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistManager.CreateTH2("hSDUsedNeutralEtaPhiPart", "#eta vs. #phi for clusters used in SD (part. level); #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistManager.CreateTH2("hSDUsedNeutralEtaPhiDet", "#eta vs. #phi for clusters used in SD (det. level); #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistManager.CreateTH2("hSDUsedChargedEtaPhiMaxPart", "#eta-phi for tracks used in SD (part. level); #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistManager.CreateTH2("hSDUsedChargedEtaPhiMaxDet", "#eta-phi for tracks used in SD (det. level); #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistManager.CreateTH2("hSDUsedNeutralEtaPhiMaxPart", "#eta vs. #phi for clusters used in SD (part. level); #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistManager.CreateTH2("hSDUsedNeutralEtaPhiMaxDet", "#eta vs. #phi for clusters used in SD (det. level); #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistManager.CreateTH2("hSDUsedChargedDRPart", "#DeltaR vs. p_{t,jet} for tracks used in SD (part. level); p_{t,jet}; #DeltaR", 350, 0., 350., 100, 0., 1.);
  fHistManager.CreateTH2("hSDUsedChargedDRDet", "#DeltaR vs. p_{t,jet} for tracks used in SD (det. level); p_{t,jet}; #DeltaR", 350, 0., 350., 100, 0., 1.);
  fHistManager.CreateTH2("hSDUsedNeutralDRPart", "#DeltaR vs. p_{t,jet} for clusters used in SD (part. level); p_{t,jet}; #DeltaR", 350, 0., 350., 100, 0., 1.);
  fHistManager.CreateTH2("hSDUsedNeutralDRDet", "#DeltaR vs. p_{t,jet} for clusters used in SD (det. level); p_{t,jet}; #DeltaR", 350, 0., 350., 100, 0., 1.);
  fHistManager.CreateTH2("hSDUsedChargedDRMaxPart", "#DeltaR vs. p_{t,jet} for tracks used in SD (part. level); p_{t,jet}; #DeltaR", 350, 0., 350., 100, 0., 1.);
  fHistManager.CreateTH2("hSDUsedChargedDRMaxDet", "#DeltaR vs. p_{t,jet} for tracks used in SD (det. level); p_{t,jet}; #DeltaR", 350, 0., 350., 100, 0., 1.);
  fHistManager.CreateTH2("hSDUsedNeutralDRMaxPart", "#DeltaR vs. p_{t,jet} for clusters used in SD (part. level); p_{t,jet}; #DeltaR", 350, 0., 350., 100, 0., 1.);
  fHistManager.CreateTH2("hSDUsedNeutralDRMaxDet", "#DeltaR vs. p_{t,jet} for clusters used in SD (det. level); p_{t,jet}; #DeltaR", 350, 0., 350., 100, 0., 1.);
  // Cluster constituent QA
  fHistManager.CreateTH2("hSDUsedClusterTimeVsE", "Cluster time vs. energy; time (ns); E (GeV)", 1200, -600, 600, 200, 0., 200);
  fHistManager.CreateTH2("hSDUsedClusterTimeVsEFine", "Cluster time vs. energy (main region); time (ns); E (GeV)", 1000, -100, 100, 200, 0., 200);
  fHistManager.CreateTH2("hSDUsedClusterNCellVsE", "Cluster number of cells vs. energy; Number of cells; E (GeV)", 201, -0.5, 200.5, 200, 0., 200.);
  fHistManager.CreateTH2("hSDUsedlusterM02VsE", "Cluster M02 vs energy; M02; E (GeV)", 150, 0., 1.5, 200, 0., 200.);
  fHistManager.CreateTH2("hSDUsedClusterFracLeadingVsE", "Cluster frac leading cell vs energy; E (GeV); Frac. leading cell", 200, 0., 200., 110, 0., 1.1);
  fHistManager.CreateTH2("hSDUsedClusterFracLeadingVsNcell", "Cluster frac leading cell vs number of cells; Number of cells; Frac. leading cell", 201, -0.5, 200.5, 110, 0., 1.1);
  fHistManager.CreateTH1("hFracPtHardPart", "Part. level jet Pt relative to the Pt-hard of the event", 100, 0., 10.);
  fHistManager.CreateTH1("hFracPtHardDet", "Det. level jet Pt relative to the Pt-hard of the event", 100, 0., 10.);

  for (auto h : *fHistManager.GetListOfHistograms())
    fOutput->Add(h);
  for (auto r : fZgResponse)
    fOutput->Add(r);
  for (auto r : fZgResponseClosure)
    fOutput->Add(r);
  for (auto r : fRgResponse)
    fOutput->Add(r);
  for (auto r : fRgResponseClosure)
    fOutput->Add(r);
  for (auto r : fNsdResponse)
    fOutput->Add(r);
  for (auto r : fNsdResponseClosure)
    fOutput->Add(r);
  for (auto r : fThetagResponse)
    fOutput->Add(r);
  for (auto r : fThetagResponseClosure)
    fOutput->Add(r);

  PostData(1, fOutput);
}

Bool_t AliAnalysisTaskEmcalSoftDropResponse::CheckMCOutliers()
{
  if (!fMCRejectFilter)
    return true;
  if (!(fIsPythia || fIsHerwig))
    return true; // Only relevant for pt-hard production
  if(fUseStandardOutlierRejection) 
    return AliAnalysisTaskEmcal::CheckMCOutliers();
  AliDebugStream(1) << "Using custom MC outlier rejection" << std::endl;
  AliJetContainer *outlierjets(nullptr);
  switch (fJetTypeOutliers)
  {
  case kOutlierPartJet:
    outlierjets = GetJetContainer(fNamePartLevelJetContainer);
    break;
  case kOutlierDetJet:
    outlierjets = GetJetContainer(fNameDetLevelJetContainer);
    break;
  
  default:
    break;
  }
  if (!outlierjets)
    return true;

  // Check whether there is at least one particle level jet with pt above n * event pt-hard
  auto jetiter = outlierjets->accepted();
  auto max = std::max_element(jetiter.begin(), jetiter.end(), [](const AliEmcalJet *lhs, const AliEmcalJet *rhs) { return lhs->Pt() < rhs->Pt(); });
  if (max != jetiter.end())
  {
    // At least one jet found with pt > n * pt-hard
    AliDebugStream(1) << "Found max jet with pt " << (*max)->Pt() << " GeV/c" << std::endl;
    if ((*max)->Pt() > fPtHardAndJetPtFactor * fPtHard)
      return false;
  }
  return true;
}

bool AliAnalysisTaskEmcalSoftDropResponse::Run()
{
  enum EPointSD_t {
    kIndSDDet = 0,
    kIndPtDet = 1,
    kIndSDPart = 2,
    kIndPtPart = 3
  };
  AliJetContainer *partLevelJets = this->GetJetContainer(fNamePartLevelJetContainer),
                  *detLevelJets = GetJetContainer(fNameDetLevelJetContainer);
  AliClusterContainer *clusters = GetClusterContainer(EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::ClusterContainerNameFactory(fInputEvent->IsA() == AliAODEvent::Class()));
  AliTrackContainer *tracks = GetTrackContainer(EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::TrackContainerNameFactory(fInputEvent->IsA() == AliAODEvent::Class()));
  AliParticleContainer *particles = GetParticleContainer(fNameMCParticles.Data());
  double Rjet = detLevelJets->GetJetRadius();
  if (!(partLevelJets || detLevelJets))
  {
    AliErrorStream() << "Either of the jet containers not found" << std::endl;
    return kFALSE;
  }

  if (fSampleFraction < 1.)
  {
    if (fSampleTrimmer->Uniform() > fSampleFraction)
      return false;
  }

  // get truncations at detector level
  TString histname;
  if (fForceBeamType != kpp)
  {
    histname = TString::Format("hZgDetLevel_%d", fCentBin);
  }
  else
  {
    histname = "hZgDetLevel";
  }
  TH2D *fZgDetLevel = (TH2D *)fHistManager.FindObject(histname);
  auto ptmindet = fHasResponseMatrixRooUnfold ? fZgDetLevel->GetYaxis()->GetBinLowEdge(1) : 0.,
       ptmaxdet = fHasResponseMatrixRooUnfold ? fZgDetLevel->GetYaxis()->GetBinUpEdge(fZgDetLevel->GetYaxis()->GetNbins()) : 1000.;

  //when  embedding and doing the constituent subtraction there is an additional step to the detector (or hybrid) to particle level matching because the detector jet (or hybrid) is the constituent subtracted jet which is not matched so we need to find the unsubtracted jet that it corresponds to and get the matched jets from there
  AliJetContainer *jetContUS = nullptr;
  if (fIsEmbeddedEvent)
    jetContUS = GetJetContainer(fNameUnSubLevelJetContainer);

  for (auto detjet : detLevelJets->accepted())
  {
    AliEmcalJet *partjet = nullptr;
    //variables for embedded pbpb data
    AliEmcalJet *jetUS = nullptr;
    Int_t ilab = -1;
    //for embedding, find the unsubtracted jet and get it's matched detector level jet
    if (fIsEmbeddedEvent)
    {
      for (Int_t i = 0; i < jetContUS->GetNJets(); i++)
      {
        jetUS = jetContUS->GetJet(i);
        if (jetUS->GetLabel() == detjet->GetLabel())
        {
          ilab = i;
          break;
        }
      }
      if (ilab == -1)
        continue;
      jetUS = jetContUS->GetJet(ilab);
      partjet = jetUS->ClosestJet();
    }
    //if we aren't embedding then just find the matched jet
    else
      partjet = detjet->ClosestJet();
    if (!partjet)
      continue;
    //one extra level of matching needed for embedding to go from detector to particle level
    if (fIsEmbeddedEvent)
    {
      partjet = partjet->ClosestJet();
      if (!partjet)
        continue;
    }

    //cut on the shared pt fraction, when embedding the unsubtracted jet should be used
    Double_t fraction = 0;
    if (fIsEmbeddedEvent)
      fraction = jetContUS->GetFractionSharedPt(jetUS);
    else
      fraction = detLevelJets->GetFractionSharedPt(detjet);
    if (fraction < fMinFractionShared)
      continue;

    // sample splitting (for closure test)
    bool closureUseResponse = (fSampleSplitter->Uniform() < fFractionResponseClosure);

    // Get the softdrop response
    std::vector<double> softdropDet, softdropPart;

    // For QA histograms
    double znepart = 0., znedet = 0., zchpart = 0., zchdet = 0., nchpart = 0., nchdet = 0., nnepart = 0., nnedet = 0.; 
    if(clusters){
      auto leadcluster = detjet->GetLeadingCluster(clusters->GetArray());
      nnedet = detjet->GetNumberOfClusters();
      if(leadcluster){
        TLorentzVector ptvec;
        leadcluster->GetMomentum(ptvec, fVertex, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy());
        znedet = detjet->GetZ(ptvec.Px(), ptvec.Py(), ptvec.Pz());
      }
    }
    if(tracks){
      nchdet = detjet->GetNumberOfTracks();
      auto leadingtrack = detjet->GetLeadingTrack(tracks->GetArray());
      if(leadingtrack) zchdet = detjet->GetZ(leadingtrack->Px(), leadingtrack->Py(), leadingtrack->Pz());
    }

    try
    {
      softdropDet = MakeSoftdrop(*detjet, detLevelJets->GetJetRadius(), tracks, clusters);
      softdropPart = MakeSoftdrop(*partjet, partLevelJets->GetJetRadius(), particles, nullptr);
      auto deltaR = TMath::Abs(partjet->DeltaR(detjet));
      bool untaggedDet = softdropDet[0] < fZcut,
           untaggedPart = softdropPart[0] < fZcut;
      Double_t pointZg[4] = {softdropDet[0], detjet->Pt(), softdropPart[0], partjet->Pt()},
               pointRg[4] = {untaggedDet ? -0.01 : softdropDet[2], detjet->Pt(), untaggedPart ? -0.01 : softdropPart[2], partjet->Pt()},
               pointNsd[4] = {untaggedDet ? -1. : softdropDet[5], detjet->Pt(), untaggedPart ? -1. : softdropPart[5], partjet->Pt()},
               pointThetag[4] = {untaggedDet ? -0.05 : softdropDet[2]/Rjet, detjet->Pt(), untaggedPart ? -0.05 : softdropPart[2]/Rjet, partjet->Pt()};
      Double_t resZg = pointZg[kIndSDDet] - pointZg[kIndSDPart],
               resRg = pointRg[kIndSDDet] - pointRg[kIndSDPart],
               resThetag = pointThetag[kIndSDDet] - pointThetag[kIndSDPart],
               resNsd = pointNsd[kIndSDDet] - pointNsd[kIndSDPart];
      if (fForceBeamType != kpp)
      {
        // fill QA histograms
        fHistManager.FillTH2(Form("hQANEFPtDet_%d", fCentBin), detjet->Pt(), detjet->NEF());
        fHistManager.FillTH2(Form("hQANEFPtPart_%d", fCentBin), partjet->Pt(), partjet->NEF());
        fHistManager.FillTH2(Form("hQAEtaPhiPart_%d", fCentBin), partjet->Eta(), TVector2::Phi_0_2pi(partjet->Phi()));
        fHistManager.FillTH2(Form("hQAEtaPhiDet_%d", fCentBin), detjet->Eta(), TVector2::Phi_0_2pi(detjet->Phi()));
        fHistManager.FillTH2(Form("hQAJetAreaVsJetPtPart_%d", fCentBin), partjet->Pt(), partjet->Area());
        fHistManager.FillTH2(Form("hQAJetAreaVsJetPtDet_%d", fCentBin), detjet->Pt(), detjet->Area());
        fHistManager.FillTH2(Form("hQAJetAreaVsNEFPart_%d", fCentBin), partjet->NEF(), partjet->Area());
        fHistManager.FillTH2(Form("hQAJetAreaVsNEFDet_%d", fCentBin), detjet->NEF(), detjet->Area());
        fHistManager.FillTH2(Form("hQAJetAreaVsNConstPart_%d", fCentBin), partjet->GetNumberOfTracks(), partjet->Area());
        fHistManager.FillTH2(Form("hQAJetAreaVsNConstDet_%d", fCentBin), detjet->GetNumberOfClusters() + detjet->GetNumberOfTracks(), detjet->Area());
        fHistManager.FillTH1(Form("hQAMatchingDRAbs_%d", fCentBin), deltaR);
        fHistManager.FillTH1(Form("hQAMatchingDRel_%d", fCentBin), deltaR / detLevelJets->GetJetRadius());
        if(fUseChargedConstituents) {
          fHistManager.FillTH2(Form("hQAZchPtDet_%d", fCentBin), detjet->Pt(), zchdet);
          fHistManager.FillTH2(Form("hQAZchPtPart_%d", fCentBin), detjet->Pt(), zchpart);
          fHistManager.FillTH2(Form("hQANChPtDet_%d", fCentBin), detjet->Pt(), nchdet);
          fHistManager.FillTH2(Form("hQANChPtPart_%d", fCentBin), partjet->Pt(), nchpart);

        }
        if(fUseNeutralConstituents) {
          fHistManager.FillTH2(Form("hQAZnePtDet_%d", fCentBin), detjet->Pt(), znedet);
          fHistManager.FillTH2(Form("hQAZnePtPart_%d", fCentBin), partjet->Pt(), znepart);
          fHistManager.FillTH2(Form("hQANnePtDet_%d", fCentBin), detjet->Pt(), nnedet);
          fHistManager.FillTH2(Form("hQANnePtPart_%d", fCentBin), partjet->Pt(), nnepart);
        }

        if(fHasResponseMatrixRooUnfold){
          fHistManager.FillTH1(Form("hZgPartLevel_%d", fCentBin), pointZg[kIndSDPart], pointZg[kIndPtPart]);
          fHistManager.FillTH1(Form("hRgPartLevel_%d", fCentBin), pointRg[kIndSDPart], pointRg[kIndPtPart]);
          fHistManager.FillTH1(Form("hNsdPartLevel_%d", fCentBin), pointNsd[kIndSDPart], pointNsd[kIndPtPart]);
          fHistManager.FillTH1(Form("hThetagPartLevel_%d", fCentBin), pointThetag[kIndSDPart], pointThetag[kIndPtPart]);
        }

        if(!untaggedDet && !untaggedPart) {
          // Fill residuals (vs pt, Rg, Thetag)
          fHistManager.FillTH2(Form("hZgResiduals_%d", fCentBin), pointZg[kIndPtPart], resZg);
          fHistManager.FillTH2(Form("hRgResiduals_%d", fCentBin), pointRg[kIndPtPart], resRg);
          fHistManager.FillTH2(Form("hThetagResiduals_%d", fCentBin),  pointThetag[kIndPtPart], resThetag);
          fHistManager.FillTH2(Form("hNsdResiduals_%d", fCentBin),  pointNsd[kIndPtPart], resNsd);
          fHistManager.FillTH2(Form("hZgResidualsNormalized_%d", fCentBin), pointZg[kIndPtPart], resZg/pointZg[kIndSDPart]);
          fHistManager.FillTH2(Form("hRgResidualsNormalized_%d", fCentBin), pointRg[kIndPtPart], resRg/pointRg[kIndSDPart]);
          fHistManager.FillTH2(Form("hThetagResidualsNormalized_%d", fCentBin), pointThetag[kIndPtPart], resThetag/pointThetag[kIndSDPart]);
          fHistManager.FillTH2(Form("hNsdResidualsNormalized_%d", fCentBin), pointNsd[kIndPtPart], resNsd/pointNsd[kIndSDPart]);
          fHistManager.FillTH2(Form("hZgResidualsRg_%d", fCentBin), pointRg[kIndSDPart], resZg);
          fHistManager.FillTH2(Form("hRgResidualsRg_%d", fCentBin), pointRg[kIndSDPart], resRg);
          fHistManager.FillTH2(Form("hThetagResidualsRg_%d", fCentBin), pointRg[kIndSDPart], resThetag);
          fHistManager.FillTH2(Form("hNsdResidualsRg_%d", fCentBin), pointRg[kIndSDPart], resNsd);
          fHistManager.FillTH2(Form("hZgResidualsRgNormalized_%d", fCentBin), pointRg[kIndSDPart], resZg/pointZg[kIndSDPart]);
          fHistManager.FillTH2(Form("hRgResidualsRgNormalized_%d", fCentBin), pointRg[kIndSDPart], resRg/pointRg[kIndSDPart]);
          fHistManager.FillTH2(Form("hThetagResidualsRgNormalized_%d", fCentBin), pointRg[kIndSDPart], resThetag/pointThetag[kIndSDPart]);
          fHistManager.FillTH2(Form("hNsdResidualsRgNormalized_%d", fCentBin), pointRg[kIndSDPart], resNsd/pointNsd[kIndSDPart]);
          fHistManager.FillTH2(Form("hZgResidualsThetag_%d", fCentBin), pointThetag[kIndSDPart], resZg);
          fHistManager.FillTH2(Form("hRgResidualsThetag_%d", fCentBin), pointThetag[kIndSDPart], resRg);
          fHistManager.FillTH2(Form("hThetagResidualsThetag_%d", fCentBin), pointThetag[kIndSDPart], resThetag);
          fHistManager.FillTH2(Form("hNsdResidualsThetag_%d", fCentBin), pointThetag[kIndSDPart], resNsd);
          fHistManager.FillTH2(Form("hZgResidualsThetagNormalized_%d", fCentBin), pointThetag[kIndSDPart], resZg/pointZg[kIndSDPart]);
          fHistManager.FillTH2(Form("hRgResidualsThetagNormalized_%d", fCentBin), pointThetag[kIndSDPart], resRg/pointRg[kIndSDPart]);
          fHistManager.FillTH2(Form("hThetagResidualsThetagNormalized_%d", fCentBin), pointThetag[kIndSDPart], resThetag/pointThetag[kIndSDPart]);
          fHistManager.FillTH2(Form("hNsdResidualsThetagNormalized_%d", fCentBin), pointThetag[kIndSDPart], resNsd/pointNsd[kIndSDPart]);
        }
      }
      else
      {
        // fill QA histograms
        auto stat = GetStatisticsConstituentsPart(*partjet, particles);
        nchpart = stat[0];
        zchpart = stat[1];
        nnepart = stat[2];
        znepart = stat[3];
        fHistManager.FillTH2("hQANEFPtDet", detjet->Pt(), detjet->NEF());
        fHistManager.FillTH2("hQANEFPtPart", partjet->Pt(), partjet->NEF());
        fHistManager.FillTH2("hQAEtaPhiPart", partjet->Eta(), TVector2::Phi_0_2pi(partjet->Phi()));
        fHistManager.FillTH2("hQAEtaPhiDet", detjet->Eta(), TVector2::Phi_0_2pi(detjet->Phi()));
        fHistManager.FillTH2("hQAJetAreaVsJetPtPart", partjet->Pt(), partjet->Area());
        fHistManager.FillTH2("hQAJetAreaVsJetPtDet", detjet->Pt(), detjet->Area());
        fHistManager.FillTH2("hQAJetAreaVsNEFPart", partjet->NEF(), partjet->Area());
        fHistManager.FillTH2("hQAJetAreaVsNEFDet", detjet->NEF(), detjet->Area());
        fHistManager.FillTH2("hQAJetAreaVsNConstPart", partjet->GetNumberOfTracks(), partjet->Area());
        fHistManager.FillTH2("hQAJetAreaVsNConstDet", detjet->GetNumberOfClusters() + detjet->GetNumberOfTracks(), detjet->Area());
        fHistManager.FillTH1("hQAMatchingDRAbs", deltaR);
        fHistManager.FillTH1("hQAMatchingDRel", deltaR / detLevelJets->GetJetRadius());
        if(fUseChargedConstituents) {
          fHistManager.FillTH2("hQAZchPtDet", detjet->Pt(), zchdet);
          fHistManager.FillTH2("hQAZchPtPart", detjet->Pt(), zchpart);
          fHistManager.FillTH2("hQANChPtDet", detjet->Pt(), nchdet);
          fHistManager.FillTH2("hQANChPtPart", partjet->Pt(), nchpart);

        }
        if(fUseNeutralConstituents) {
          fHistManager.FillTH2("hQAZnePtDet", detjet->Pt(), znedet);
          fHistManager.FillTH2("hQAZnePtPart", partjet->Pt(), znepart);
          fHistManager.FillTH2("hQANnePtDet", detjet->Pt(), nnedet);
          fHistManager.FillTH2("hQANnePtPart", partjet->Pt(), nnepart);
        }
        // Monitor jet pt relative to pt-hard of the event
        if(fPtHard > 0.){
          fHistManager.FillTH1("hFracPtHardDet", detjet->Pt()/fPtHard);
          fHistManager.FillTH1("hFracPtHardPart", partjet->Pt() / fPtHard);
        } 

        if(fHasResponseMatrixRooUnfold){
          fHistManager.FillTH1("hZgPartLevel", pointZg[kIndSDPart], pointZg[kIndPtPart]);
          fHistManager.FillTH1("hRgPartLevel", pointRg[kIndSDPart], pointRg[kIndPtPart]);
          fHistManager.FillTH1("hNsdPartLevel", pointNsd[kIndSDPart], pointNsd[kIndPtPart]);
          fHistManager.FillTH1("hThetagPartLevel", pointThetag[kIndSDPart], pointThetag[kIndPtPart]);
        }
        if(!untaggedDet && !untaggedPart) {
          // Fill residuals
          fHistManager.FillTH2("hZgResiduals", pointZg[kIndPtPart], resZg);
          fHistManager.FillTH2("hRgResiduals", pointRg[kIndPtPart], resRg);
          fHistManager.FillTH2("hThetagResiduals",  pointThetag[kIndPtPart], resThetag);
          fHistManager.FillTH2("hNsdResiduals",  pointNsd[kIndPtPart], resNsd);
          fHistManager.FillTH2("hZgResidualsNormalized", pointZg[kIndPtPart], resZg/pointZg[kIndSDPart]);
          fHistManager.FillTH2("hRgResidualsNormalized", pointRg[kIndPtPart], resRg/pointRg[kIndSDPart]);
          fHistManager.FillTH2("hThetagResidualsNormalized", pointThetag[kIndPtPart], resThetag/pointThetag[kIndSDPart]);
          fHistManager.FillTH2("hNsdResidualsNormalized", pointNsd[kIndPtPart], resNsd/pointNsd[kIndSDPart]);
          fHistManager.FillTH2("hZgResidualsRg", pointRg[kIndSDPart], resZg);
          fHistManager.FillTH2("hRgResidualsRg", pointRg[kIndSDPart], resRg);
          fHistManager.FillTH2("hThetagResidualsRg", pointRg[kIndSDPart], resThetag);
          fHistManager.FillTH2("hNsdResidualsRg", pointRg[kIndSDPart], resNsd);
          fHistManager.FillTH2("hZgResidualsRgNormalized", pointRg[kIndSDPart], resZg/pointZg[kIndSDPart]);
          fHistManager.FillTH2("hRgResidualsRgNormalized", pointRg[kIndSDPart], resRg/pointRg[kIndSDPart]);
          fHistManager.FillTH2("hThetagResidualsRgNormalized", pointRg[kIndSDPart], resThetag/pointThetag[kIndSDPart]);
          fHistManager.FillTH2("hNsdResidualsRgNormalized", pointRg[kIndSDPart], resNsd/pointNsd[kIndSDPart]);
          fHistManager.FillTH2("hZgResidualsThetag", pointThetag[kIndSDPart], resZg);
          fHistManager.FillTH2("hRgResidualsThetag", pointThetag[kIndSDPart], resRg);
          fHistManager.FillTH2("hThetagResidualsThetag", pointThetag[kIndSDPart], resThetag);
          fHistManager.FillTH2("hNsdResidualsThetag", pointThetag[kIndSDPart], resNsd);
          fHistManager.FillTH2("hZgResidualsThetagNormalized", pointThetag[kIndSDPart], resZg/pointZg[kIndSDPart]);
          fHistManager.FillTH2("hRgResidualsThetagNormalized", pointThetag[kIndSDPart], resRg/pointRg[kIndSDPart]);
          fHistManager.FillTH2("hThetagResidualsThetagNormalized", pointThetag[kIndSDPart], resThetag/pointThetag[kIndSDPart]);
          fHistManager.FillTH2("hNsdResidualsThetagNormalized", pointThetag[kIndSDPart], resNsd/pointNsd[kIndSDPart]);
        }
      }
      if (detjet->Pt() >= ptmindet && detjet->Pt() <= ptmaxdet)
      {
        if (fForceBeamType != kpp)
        {
          if(fHasResponseMatrixRooUnfold){
            fHistManager.FillTH2(Form("hZgPartLevelTruncated_%d", fCentBin), pointZg[kIndSDPart], pointZg[kIndPtPart]);
            fHistManager.FillTH2(Form("hZgDetLevel_%d", fCentBin), pointZg[kIndSDDet], pointZg[kIndPtDet]);
            fHistManager.FillTH2(Form("hRgPartLevelTruncated_%d", fCentBin), pointRg[kIndSDPart], pointRg[kIndPtPart]);
            fHistManager.FillTH2(Form("hRgDetLevel_%d", fCentBin), pointRg[kIndSDDet], pointRg[kIndPtDet]);
            fHistManager.FillTH2(Form("hNsdPartLevelTruncated_%d", fCentBin), pointNsd[kIndSDPart], pointNsd[kIndPtPart]);
            fHistManager.FillTH2(Form("hNsdDetLevel_%d", fCentBin), pointNsd[kIndSDDet], pointNsd[kIndPtDet]);
            fHistManager.FillTH2(Form("hThetagPartLevelTruncated_%d", fCentBin), pointThetag[kIndSDPart], pointThetag[kIndPtPart]);
            fHistManager.FillTH2(Form("hThetagDetLevel_%d", fCentBin), pointThetag[kIndSDDet], pointThetag[kIndPtDet]);
            fZgResponse[fCentBin]->Fill(pointZg[kIndSDDet], pointZg[kIndPtDet], pointZg[kIndSDPart], pointZg[kIndPtPart]);
            fRgResponse[fCentBin]->Fill(pointRg[kIndSDDet], pointRg[kIndPtDet], pointRg[kIndSDPart], pointRg[kIndPtPart]);
            fNsdResponse[fCentBin]->Fill(pointNsd[kIndSDDet], pointNsd[kIndPtDet], pointNsd[kIndSDPart], pointNsd[kIndPtPart]);
            fThetagResponse[fCentBin]->Fill(pointThetag[kIndSDDet], pointThetag[kIndPtDet], pointThetag[kIndSDPart], pointThetag[kIndPtPart]);
          }
          if(fHasResponseMatrixSparse){
            fHistManager.FillTHnSparse(Form("hZgResponseSparse_%d", fCentBin), pointZg);
            fHistManager.FillTHnSparse(Form("hRgResponseSparse_%d", fCentBin), pointRg);
            fHistManager.FillTHnSparse(Form("hNsdResponseSparse_%d", fCentBin), pointNsd);
            fHistManager.FillTHnSparse(Form("hThetagResponseSparse_%d", fCentBin), pointThetag);
          }
        }
        else
        {
          if(fHasResponseMatrixRooUnfold){
            fHistManager.FillTH2("hZgPartLevelTruncated", pointZg[kIndSDPart], pointZg[kIndPtPart]);
            fHistManager.FillTH2("hZgDetLevel", pointZg[kIndSDDet], pointZg[kIndPtDet]);
            fHistManager.FillTH2("hRgPartLevelTruncated", pointRg[kIndSDPart], pointRg[kIndPtPart]);
            fHistManager.FillTH2("hRgDetLevel", pointRg[kIndSDDet], pointRg[kIndPtDet]);
            fHistManager.FillTH2("hNsdPartLevelTruncated", pointNsd[kIndSDPart], pointNsd[kIndPtPart]);
            fHistManager.FillTH2("hNsdDetLevel", pointNsd[kIndSDDet], pointNsd[kIndPtDet]);
            fHistManager.FillTH2("hThetagPartLevelTruncated", pointThetag[kIndSDPart], pointThetag[kIndPtPart]);
            fHistManager.FillTH2("hThetagDetLevel", pointThetag[kIndSDDet], pointThetag[kIndPtDet]);
            fZgResponse[0]->Fill(pointZg[kIndSDDet], pointZg[kIndPtDet], pointZg[kIndSDPart], pointZg[kIndPtPart]);
            fRgResponse[0]->Fill(pointRg[kIndSDDet], pointRg[kIndPtDet], pointRg[kIndSDPart], pointRg[kIndPtPart]);
            fNsdResponse[0]->Fill(pointNsd[kIndSDDet], pointNsd[kIndPtDet], pointNsd[kIndSDPart], pointNsd[kIndPtPart]);
            fThetagResponse[0]->Fill(pointThetag[kIndSDDet], pointThetag[kIndPtDet], pointThetag[kIndSDPart], pointThetag[kIndPtPart]);
          }
          if(fHasResponseMatrixSparse){
            fHistManager.FillTHnSparse("hZgResponseSparse", pointZg);
            fHistManager.FillTHnSparse("hRgResponseSparse", pointRg);
            fHistManager.FillTHnSparse("hNsdResponseSparse", pointNsd);
            fHistManager.FillTHnSparse("hThetagResponseSparse", pointThetag); 
          }
        }
        if (closureUseResponse)
        {
          if (fForceBeamType != kpp)
          {
            if(fHasResponseMatrixRooUnfold){
              fZgResponseClosure[fCentBin]->Fill(pointZg[kIndSDDet], pointZg[kIndPtDet], pointZg[kIndSDPart], pointZg[kIndPtPart]);
              fRgResponseClosure[fCentBin]->Fill(pointRg[kIndSDDet], pointRg[kIndPtDet], pointRg[kIndSDPart], pointRg[kIndPtPart]);
              fNsdResponseClosure[fCentBin]->Fill(pointNsd[kIndSDDet], pointNsd[kIndPtDet], pointNsd[kIndSDPart], pointNsd[kIndPtPart]);
              fThetagResponseClosure[fCentBin]->Fill(pointThetag[kIndSDDet], pointThetag[kIndPtDet], pointThetag[kIndSDPart], pointThetag[kIndPtPart]);
              fHistManager.FillTH2(Form("hZgDetLevelClosureResp_%d", fCentBin), pointZg[kIndSDDet], pointZg[kIndPtDet]);
              fHistManager.FillTH2(Form("hZgPartLevelClosureResp_%d", fCentBin), pointZg[kIndSDPart], pointZg[kIndPtPart]);
              fHistManager.FillTH2(Form("hRgDetLevelClosureResp_%d", fCentBin), pointRg[kIndSDDet], pointRg[kIndPtDet]);
              fHistManager.FillTH2(Form("hRgPartLevelClosureResp_%d", fCentBin), pointRg[kIndSDPart], pointRg[kIndPtPart]);
              fHistManager.FillTH2(Form("hNsdDetLevelClosureResp_%d", fCentBin), pointNsd[kIndSDDet], pointNsd[kIndPtDet]);
              fHistManager.FillTH2(Form("hNsdPartLevelClosureResp_%d", fCentBin), pointNsd[kIndSDPart], pointNsd[kIndPtPart]);
              fHistManager.FillTH2(Form("hThetagDetLevelClosureResp_%d", fCentBin), pointThetag[kIndSDDet], pointThetag[kIndPtDet]);
              fHistManager.FillTH2(Form("hThetagPartLevelClosureResp_%d", fCentBin), pointThetag[kIndSDPart], pointThetag[kIndPtPart]);
            }
            if(fHasResponseMatrixSparse){
              fHistManager.FillTHnSparse(Form("hZgResponseClosureSparse_%d", fCentBin), pointZg);
              fHistManager.FillTHnSparse(Form("hRgResponseClosureSparse_%d", fCentBin), pointRg);
              fHistManager.FillTHnSparse(Form("hNsdResponseClosureSparse_%d", fCentBin), pointNsd);
              fHistManager.FillTHnSparse(Form("hThetagResponseClosureSparse_%d", fCentBin), pointThetag);
            }
          }
          else
          {
            if(fHasResponseMatrixRooUnfold){
              fZgResponseClosure[0]->Fill(pointZg[kIndSDDet], pointZg[kIndPtDet], pointZg[kIndSDPart], pointZg[kIndPtPart]);
              fRgResponseClosure[0]->Fill(pointRg[kIndSDDet], pointRg[kIndPtDet], pointRg[kIndSDPart], pointRg[kIndPtPart]);
              fNsdResponseClosure[0]->Fill(pointNsd[kIndSDDet], pointNsd[kIndPtDet], pointNsd[kIndSDPart], pointNsd[kIndPtPart]);
              fThetagResponseClosure[0]->Fill(pointThetag[kIndSDDet], pointThetag[kIndPtDet], pointThetag[kIndSDPart], pointThetag[kIndPtPart]);
              fHistManager.FillTH2("hZgDetLevelClosureResp", pointZg[kIndSDDet], pointZg[kIndPtDet]);
              fHistManager.FillTH2("hZgPartLevelClosureResp", pointZg[kIndSDPart], pointZg[kIndPtPart]);
              fHistManager.FillTH2("hRgDetLevelClosureResp", pointRg[kIndSDDet], pointRg[kIndPtDet]);
              fHistManager.FillTH2("hRgPartLevelClosureResp", pointRg[kIndSDPart], pointRg[kIndPtPart]);
              fHistManager.FillTH2("hNsdDetLevelClosureResp", pointNsd[kIndSDDet], pointNsd[kIndPtDet]);
              fHistManager.FillTH2("hNsdPartLevelClosureResp", pointNsd[kIndSDPart], pointNsd[kIndPtPart]);
              fHistManager.FillTH2("hThetagDetLevelClosureResp", pointThetag[kIndSDDet], pointThetag[kIndPtDet]);
              fHistManager.FillTH2("hThetagPartLevelClosureResp", pointThetag[kIndSDPart], pointThetag[kIndPtPart]);
            }
            if(fHasResponseMatrixSparse){
              fHistManager.FillTHnSparse("hZgResponseClosureSparse", pointZg);
              fHistManager.FillTHnSparse("hRgResponseClosureSparse", pointRg);
              fHistManager.FillTHnSparse("hNsdResponseClosureSparse", pointNsd);
              fHistManager.FillTHnSparse("hThetagResponseClosureSparse", pointThetag); 
            }
          }
        }
        else
        {
          if (fForceBeamType != kpp)
          {
            if(fHasResponseMatrixSparse){
              fHistManager.FillTH2(Form("hZgPartLevelClosureNoResp_%d", fCentBin), pointZg[kIndSDPart], pointZg[kIndPtPart]);
              fHistManager.FillTH2(Form("hZgDetLevelClosureNoResp_%d", fCentBin), pointZg[kIndSDDet], pointZg[kIndPtDet]);
              fHistManager.FillTH2(Form("hRgPartLevelClosureNoResp_%d", fCentBin), pointRg[kIndSDPart], pointRg[kIndPtPart]);
              fHistManager.FillTH2(Form("hRgDetLevelClosureNoResp_%d", fCentBin), pointRg[kIndSDDet], pointRg[kIndPtDet]);
              fHistManager.FillTH2(Form("hNsdPartLevelClosureNoResp_%d", fCentBin), pointNsd[kIndSDPart], pointNsd[kIndPtPart]);
              fHistManager.FillTH2(Form("hNsdDetLevelClosureNoResp_%d", fCentBin), pointNsd[kIndSDDet], pointNsd[kIndPtDet]);
              fHistManager.FillTH2(Form("hThetagPartLevelClosureNoResp_%d", fCentBin), pointThetag[kIndSDPart], pointThetag[kIndPtPart]);
              fHistManager.FillTH2(Form("hThetagDetLevelClosureNoResp_%d", fCentBin), pointThetag[kIndSDDet], pointThetag[kIndPtDet]);
            }
            if(fHasResponseMatrixSparse){
              fHistManager.FillTHnSparse(Form("hZgResponseClosureNoRespSparse_%d", fCentBin), pointZg);
              fHistManager.FillTHnSparse(Form("hRgResponseClosureNoRespSparse_%d", fCentBin), pointRg);
              fHistManager.FillTHnSparse(Form("hNsdResponseClosureNoRespSparse_%d", fCentBin), pointNsd);
              fHistManager.FillTHnSparse(Form("hThetagResponseClosureNoRespSparse_%d", fCentBin), pointThetag); 
              fHistManager.FillTH2(Form("hZgPartLevelClosureNoRespFine_%d", fCentBin), pointZg[kIndSDPart], pointZg[kIndPtPart]);
              fHistManager.FillTH2(Form("hZgDetLevelClosureNoRespFine_%d", fCentBin), pointZg[kIndSDDet], pointZg[kIndPtDet]);
              fHistManager.FillTH2(Form("hRgPartLevelClosureNoRespFine_%d", fCentBin), pointRg[kIndSDPart], pointRg[kIndPtPart]);
              fHistManager.FillTH2(Form("hRgDetLevelClosureNoRespFine_%d", fCentBin), pointRg[kIndSDDet], pointRg[kIndPtDet]);
              fHistManager.FillTH2(Form("hNsdPartLevelClosureNoRespFine_%d", fCentBin), pointNsd[kIndSDPart], pointNsd[kIndPtPart]);
              fHistManager.FillTH2(Form("hNsdDetLevelClosureNoRespFine_%d", fCentBin), pointNsd[kIndSDDet], pointNsd[kIndPtDet]);
              fHistManager.FillTH2(Form("hThetagPartLevelClosureNoRespFine_%d", fCentBin), pointThetag[kIndSDPart], pointThetag[kIndPtPart]);
              fHistManager.FillTH2(Form("hThetagDetLevelClosureNoRespFine_%d", fCentBin), pointThetag[kIndSDDet], pointThetag[kIndPtDet]);
            }
          }
          else
          {
            if(fHasResponseMatrixRooUnfold){
              fHistManager.FillTH2("hZgDetLevelClosureNoResp", pointZg[kIndSDDet], pointZg[kIndPtDet]);
              fHistManager.FillTH2("hZgPartLevelClosureNoResp", pointZg[kIndSDPart], pointZg[kIndPtPart]);
              fHistManager.FillTH2("hRgDetLevelClosureNoResp", pointRg[kIndSDDet], pointRg[kIndPtDet]);
              fHistManager.FillTH2("hRgPartLevelClosureNoResp", pointRg[kIndSDPart], pointRg[kIndPtPart]);
              fHistManager.FillTH2("hNsdDetLevelClosureNoResp", pointNsd[kIndSDDet], pointNsd[kIndPtDet]);
              fHistManager.FillTH2("hNsdPartLevelClosureNoResp", pointNsd[kIndSDPart], pointNsd[kIndPtPart]);
              fHistManager.FillTH2("hThetagDetLevelClosureNoResp", pointThetag[kIndSDDet], pointThetag[kIndPtDet]);
              fHistManager.FillTH2("hThetagPartLevelClosureNoResp", pointThetag[kIndSDPart], pointThetag[kIndPtPart]);
            } 
            if(fHasResponseMatrixSparse) {
              fHistManager.FillTHnSparse("hZgResponseClosureNoRespSparse", pointZg);
              fHistManager.FillTHnSparse("hRgResponseClosureNoRespSparse", pointRg);
              fHistManager.FillTHnSparse("hNsdResponseClosureNoRespSparse", pointNsd);
              fHistManager.FillTHnSparse("hThetagResponseClosureNoRespSparse", pointThetag); 
              fHistManager.FillTH2("hZgDetLevelClosureNoRespFine", pointZg[kIndSDDet], pointZg[kIndPtDet]);
              fHistManager.FillTH2("hZgPartLevelClosureNoRespFine", pointZg[kIndSDPart], pointZg[kIndPtPart]);
              fHistManager.FillTH2("hRgDetLevelClosureNoRespFine", pointRg[kIndSDDet], pointRg[kIndPtDet]);
              fHistManager.FillTH2("hRgPartLevelClosureNoRespFine", pointRg[kIndSDPart], pointRg[kIndPtPart]);
              fHistManager.FillTH2("hNsdDetLevelClosureNoRespFine", pointNsd[kIndSDDet], pointNsd[kIndPtDet]);
              fHistManager.FillTH2("hNsdPartLevelClosureNoRespFine", pointNsd[kIndSDPart], pointNsd[kIndPtPart]);
              fHistManager.FillTH2("hThetagDetLevelClosureNoRespFine", pointThetag[kIndSDDet], pointThetag[kIndPtDet]);
              fHistManager.FillTH2("hThetagPartLevelClosureNoRespFine", pointThetag[kIndSDPart], pointThetag[kIndPtPart]);
            }
          }
        }
      }
    }
    catch (...)
    {
      if(fUseChargedConstituents && fUseNeutralConstituents) AliErrorStream() << "Error in softdrop evaluation - jet will be ignored" << std::endl;
      if(fForceBeamType != kpp) {
        fHistManager.FillTH1(Form("hSkippedJetsPart_%d", fCentBin), partjet->Pt());
        fHistManager.FillTH1(Form("hSkippedJetsDet_%d", fCentBin), detjet->Pt());
      } else {
        fHistManager.FillTH1("hSkippedJetsPart", partjet->Pt());
        fHistManager.FillTH1("hSkippedJetsDet", detjet->Pt());
      }
      continue;
    }
  }
  return kTRUE;
}

std::vector<double> AliAnalysisTaskEmcalSoftDropResponse::MakeSoftdrop(const AliEmcalJet &jet, double jetradius, const AliParticleContainer *tracks, const AliClusterContainer *clusters) 
{
  const int kClusterOffset = 30000; // In order to handle tracks and clusters in the same index space the cluster index needs and offset, large enough so that there is no overlap with track indices
  std::vector<fastjet::PseudoJet> constituents;
  bool isMC = dynamic_cast<const AliMCParticleContainer *>(tracks);
  AliDebugStream(2) << "Make new jet substrucutre for " << (isMC ? "MC" : "data") << " jet: Number of tracks " << jet.GetNumberOfTracks() << ", clusters " << jet.GetNumberOfClusters() << std::endl;
  fastjet::PseudoJet *maxcharged(nullptr), *maxneutral(nullptr);
  fastjet::PseudoJet inputjet(jet.Px(), jet.Py(), jet.Pz(), jet.E());
  if (tracks && (fUseChargedConstituents || isMC))
  { // Neutral particles part of particle container in case of MC
    AliDebugStream(1) << "Jet substructure: Using charged constituents" << std::endl;
    for (int itrk = 0; itrk < jet.GetNumberOfTracks(); itrk++)
    {
      auto track = jet.Track(itrk);
      if (!track->Charge() && !fUseNeutralConstituents)
        continue; // Reject neutral constituents in case of using only charged consituents
      if (track->Charge() && !fUseChargedConstituents)
        continue; // Reject charged constituents in case of using only neutral consituents
      fastjet::PseudoJet constituentTrack(track->Px(), track->Py(), track->Pz(), track->E());
      constituentTrack.set_user_index(jet.TrackAt(itrk));
      constituents.push_back(constituentTrack);
      auto &currentconstituent = constituents.back();
      if(isMC) {
        if(track->Charge()) {
          fHistManager.FillTH2("hSDUsedChargedPtjvPtcPart", jet.Pt(), constituentTrack.pt());
          fHistManager.FillTH2("hSDUsedChargedEtaPhiPart", constituentTrack.eta(), TVector2::Phi_0_2pi(constituentTrack.phi()));
          fHistManager.FillTH2("hSDUsedChargedDRPart", inputjet.pt(), inputjet.delta_R(constituentTrack));
          if(!maxcharged) {
            maxcharged = &currentconstituent;
          } else {
            if(currentconstituent.pt() > maxcharged->pt())
            maxcharged = &currentconstituent;
          }
        } else {
          fHistManager.FillTH2("hSDUsedNeutralPtjvPtcPart", jet.Pt(), constituentTrack.pt());
          fHistManager.FillTH2("hSDUsedNeutralEtaPhiPart", constituentTrack.eta(), TVector2::Phi_0_2pi(constituentTrack.phi()));
          fHistManager.FillTH2("hSDUsedNeutralDRPart", inputjet.pt(), inputjet.delta_R(constituentTrack));
          if(!maxneutral) {
            maxneutral = &currentconstituent;
          } else {
            if(currentconstituent.pt() > maxneutral->pt())
            maxneutral = &currentconstituent;
          }
        }
      } else {
        fHistManager.FillTH2("hSDUsedChargedPtjvPtcDet", jet.Pt(), constituentTrack.pt());
        fHistManager.FillTH2("hSDUsedChargedEtaPhiDet", constituentTrack.eta(), TVector2::Phi_0_2pi(constituentTrack.phi()));
        fHistManager.FillTH2("hSDUsedChargedDRDet", inputjet.pt(), inputjet.delta_R(constituentTrack));
        if(!maxcharged) {
          maxcharged = &currentconstituent;
        } else {
          if(currentconstituent.pt() > maxcharged->pt())
          maxcharged = &currentconstituent;
        }
      }
    }
  }

  if (clusters && fUseNeutralConstituents)
  {
    AliDebugStream(1) << "Jet substructure: Using neutral constituents" << std::endl;
    for (int icl = 0; icl < jet.GetNumberOfClusters(); icl++)
    {
      auto cluster = jet.ClusterAt(icl, clusters->GetArray());
      TLorentzVector clustervec;
      cluster->GetMomentum(clustervec, fVertex, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy());
      fastjet::PseudoJet constituentCluster(clustervec.Px(), clustervec.Py(), clustervec.Pz(), cluster->GetHadCorrEnergy());
      constituentCluster.set_user_index(jet.ClusterAt(icl) + kClusterOffset);
      constituents.push_back(constituentCluster);
      fHistManager.FillTH2("hSDUsedNeutralPtjvPtcDet", jet.Pt(), constituentCluster.pt());
      fHistManager.FillTH2("hSDUsedNeutralEtaPhiDet", constituentCluster.eta(), TVector2::Phi_0_2pi(constituentCluster.phi()));
      fHistManager.FillTH2("hSDUsedNeutralDRDet", inputjet.pt(), inputjet.delta_R(constituentCluster));
      fHistManager.FillTH2("hSDUsedClusterTimeVsE", cluster->GetTOF() * 1e9 - 600, clustervec.E());       // time in ns., apply 600 ns time shift
      fHistManager.FillTH2("hSDUsedClusterTimeVsEFine", cluster->GetTOF() * 1e9 - 600, clustervec.E());   // time in ns., apply 600 ns time shift
      fHistManager.FillTH2("hSDUsedClusterNCellVsE", cluster->GetNCells(), clustervec.E());
      fHistManager.FillTH2("hSDUsedlusterM02VsE", cluster->GetM02(), clustervec.E());
      double maxamplitude = 0.;
      for(int icell = 0; icell < cluster->GetNCells(); icell++) {
        double amplitude = fInputEvent->GetEMCALCells()->GetAmplitude(fInputEvent->GetEMCALCells()->GetCellPosition(cluster->GetCellAbsId(icell)));
        if(amplitude > maxamplitude) maxamplitude = amplitude;
      }
      fHistManager.FillTH2("hSDUsedClusterFracLeadingVsE", clustervec.E(), maxamplitude/cluster->E());
      fHistManager.FillTH2("hSDUsedClusterFracLeadingVsNcell", cluster->GetNCells(), maxamplitude/cluster->E());

      auto &currentconstituent = constituents.back();
      if(!maxneutral) {
        maxneutral = &currentconstituent;
      } else {
        if(currentconstituent.pt() > maxneutral->pt())
          maxneutral = &currentconstituent;
      }
    }
  }

  if(maxcharged) {
    if(isMC) {
      fHistManager.FillTH2("hSDUsedChargedPtjvPtcMaxPart", jet.Pt(), maxcharged->pt());
      fHistManager.FillTH2("hSDUsedChargedEtaPhiMaxPart", maxcharged->eta(), TVector2::Phi_0_2pi(maxcharged->phi()));
      fHistManager.FillTH2("hSDUsedChargedDRMaxPart", inputjet.pt(), inputjet.delta_R(*maxcharged));
    } else {
      fHistManager.FillTH2("hSDUsedChargedPtjvPtcMaxDet", jet.Pt(), maxcharged->pt());
      fHistManager.FillTH2("hSDUsedChargedEtaPhiMaxDet", maxcharged->eta(), TVector2::Phi_0_2pi(maxcharged->phi()));
      fHistManager.FillTH2("hSDUsedChargedDRMaxDet", inputjet.pt(), inputjet.delta_R(*maxcharged));
    }
  }

  if(maxneutral) {
    if(isMC) {
      fHistManager.FillTH2("hSDUsedNeutralPtjvPcMaxPart", jet.Pt(), maxneutral->pt());
      fHistManager.FillTH2("hSDUsedNeutralEtaPhiMaxPart", maxneutral->eta(), TVector2::Phi_0_2pi(maxneutral->phi()));
      fHistManager.FillTH2("hSDUsedNeutralDRMaxPart", inputjet.pt(), inputjet.delta_R(*maxneutral));
    } else {
      fHistManager.FillTH2("hSDUsedNeutralPtjvPcMaxDet", jet.Pt(), maxneutral->pt());
      fHistManager.FillTH2("hSDUsedNeutralEtaPhiMaxDet", maxneutral->eta(), TVector2::Phi_0_2pi(maxneutral->phi()));
      fHistManager.FillTH2("hSDUsedNeutralDRMaxDet", inputjet.pt(), inputjet.delta_R(*maxneutral));
    }
  }


  AliDebugStream(3) << "Found " << constituents.size() << " constituents for jet with pt=" << jet.Pt() << " GeV/c" << std::endl;
  if (!constituents.size())
  {
    if(fUseChargedConstituents && fUseNeutralConstituents) AliErrorStream() << "Jet has 0 constituents." << std::endl;
    throw 1;
  }
  // Redo jet finding on constituents with a
  fastjet::JetDefinition jetdef(fastjet::antikt_algorithm, jetradius * 2, static_cast<fastjet::RecombinationScheme>(0), fastjet::BestFJ30);
  fastjet::ClusterSequence jetfinder(constituents, jetdef);
  std::vector<fastjet::PseudoJet> outputjets = jetfinder.inclusive_jets(0);
  auto sdjet = outputjets[0];
  fastjet::contrib::SoftDrop softdropAlgorithm(fBeta, fZcut);
  softdropAlgorithm.set_verbose_structure(kTRUE);
  fastjet::JetAlgorithm reclusterizingAlgorithm;
  switch (fReclusterizer)
  {
  case kCAAlgo:
    reclusterizingAlgorithm = fastjet::cambridge_aachen_algorithm;
    break;
  case kKTAlgo:
    reclusterizingAlgorithm = fastjet::kt_algorithm;
    break;
  case kAKTAlgo:
    reclusterizingAlgorithm = fastjet::antikt_algorithm;
    break;
  };
#if FASTJET_VERSION_NUMBER >= 30302
  fastjet::Recluster reclusterizer(reclusterizingAlgorithm, 1, fastjet::Recluster::keep_only_hardest);
#else
  fastjet::contrib::Recluster reclusterizer(reclusterizingAlgorithm, 1, true);
#endif
  softdropAlgorithm.set_reclustering(kTRUE, &reclusterizer);
  AliDebugStream(4) << "Jet has " << sdjet.constituents().size() << " constituents" << std::endl;
  auto groomed = softdropAlgorithm(sdjet);
  auto softdropstruct = groomed.structure_of<fastjet::contrib::SoftDrop>();

  std::vector<double> result = {softdropstruct.symmetry(),
                                groomed.m(),
                                softdropstruct.delta_R(),
                                groomed.perp(),
                                softdropstruct.mu(),
                                static_cast<double>(softdropstruct.dropped_count())};
  return result;
}

std::vector<double> AliAnalysisTaskEmcalSoftDropResponse::GetStatisticsConstituentsPart(const AliEmcalJet &jet, const AliParticleContainer *particles) const {
  AliVParticle *leadingcharged = nullptr, *leadingneutral = nullptr;
  int ncharged = 0, nneutral = 0;
  for(auto ipart = 0; ipart < jet.GetNumberOfTracks(); ipart++){
    auto part = jet.TrackAt(ipart, particles->GetArray());
    if(part->Charge()) {
      ncharged++;
      if(!leadingcharged) leadingcharged = part;
      else if(part->E() > leadingcharged->E()) leadingcharged = part;
    } else {
      nneutral++;
      if(!leadingneutral) leadingneutral = part;
      else if(part->E() > leadingneutral->E()) leadingneutral = part;
    }
  }

  // Calculate z
  Double_t zch = 0, zne = 0;
  if(leadingcharged) zch = jet.GetZ(leadingcharged);
  if(leadingneutral) zne = jet.GetZ(leadingneutral);
  return {static_cast<double>(ncharged), zch, static_cast<double>(nneutral), zne};
}

TBinning *AliAnalysisTaskEmcalSoftDropResponse::GetDefaultPartLevelPtBinning() const
{
  auto binning = new TCustomBinning;
  binning->SetMinimum(0);
  switch (fBinningMode)
  {
  case kSDModeINT7:
  {
    binning->AddStep(20., 20.);
    binning->AddStep(40., 10.);
    binning->AddStep(80., 20.);
    binning->AddStep(120., 40.);
    binning->AddStep(240., 120.);
    break;
  }
  case kSDModeEJ1:
  {
    binning->AddStep(80., 80.);
    binning->AddStep(140., 10.);
    binning->AddStep(200., 20.);
    binning->AddStep(240., 40.);
    binning->AddStep(400., 160.);
    break;
  }
  case kSDModeEJ2:
  {
    binning->AddStep(70., 70.);
    binning->AddStep(100., 10.);
    binning->AddStep(140., 20.);
    binning->AddStep(400., 260.);
    break;
  }
  };
  return binning;
}

TBinning *AliAnalysisTaskEmcalSoftDropResponse::GetDefaultDetLevelPtBinning() const
{
  auto binning = new TCustomBinning;
  switch (fBinningMode)
  {
  case kSDModeINT7:
  {
    binning->SetMinimum(20);
    binning->AddStep(40., 5.);
    binning->AddStep(60., 10.);
    binning->AddStep(80., 20.);
    binning->AddStep(120., 40.);
    break;
  }
  case kSDModeEJ1:
  {
    binning->SetMinimum(80.);
    binning->AddStep(120., 5.);
    binning->AddStep(160., 10.);
    binning->AddStep(200., 20.);
    break;
  }
  case kSDModeEJ2:
  {
    binning->SetMinimum(70.);
    binning->AddStep(100., 5.);
    binning->AddStep(120., 10.);
    binning->AddStep(140., 20.);
    break;
  }
  };
  return binning;
}

TBinning *AliAnalysisTaskEmcalSoftDropResponse::GetZgBinning() const
{
  auto binning = new TCustomBinning;
  binning->SetMinimum(0.);
  binning->AddStep(0.1, 0.1);
  binning->AddStep(0.5, 0.05);
  return binning;
}

TBinning *AliAnalysisTaskEmcalSoftDropResponse::GetRgBinning(double R) const {
  auto binning = new TCustomBinning;
  binning->SetMinimum(-0.05);    // Negative bins are for untagged jets
  binning->AddStep(R, 0.05);
  return binning;
}

AliAnalysisTaskEmcalSoftDropResponse *AliAnalysisTaskEmcalSoftDropResponse::AddTaskEmcalSoftDropResponse(Double_t jetradius, AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recombinationScheme, bool ifembed, const char *namepartcont, const char *trigger)
{
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();

  Bool_t isAOD(kFALSE);
  AliInputEventHandler *inputhandler = static_cast<AliInputEventHandler *>(mgr->GetInputEventHandler());
  if (inputhandler)
  {
    if (inputhandler->IsA() == AliAODInputHandler::Class())
    {
      std::cout << "Analysing AOD events\n";
      isAOD = kTRUE;
    }
    else
    {
      std::cout << "Analysing ESD events\n";
    }
  }

  std::stringstream taskname;
  taskname << "SoftdropResponsemaker_R" << std::setw(2) << std::setfill('0') << int(jetradius * 10) << trigger;
  AliAnalysisTaskEmcalSoftDropResponse *responsemaker = new AliAnalysisTaskEmcalSoftDropResponse(taskname.str().data());
  responsemaker->SelectCollisionCandidates(AliVEvent::kINT7);
  responsemaker->SetIsEmbeddedEvent(ifembed);
  mgr->AddTask(responsemaker);

  TString partcontname(namepartcont);
  if (partcontname == "usedefault")
    partcontname = "mcparticles";
  AliParticleContainer *particles = responsemaker->AddMCParticleContainer(partcontname.Data());
  //if embedding need to specify that the particles are embedded
  if (ifembed)
    particles->SetIsEmbedding(true);
  particles->SetMinPt(0.);
  responsemaker->SetNameMCParticleContainer(partcontname.Data());

  TString partLevelTag("Jet");
  //if embedding need to specify the tag of the particle level jet from the embedding framework
  if (ifembed)
    partLevelTag = "partLevelJets";
  AliJetContainer *mcjets = responsemaker->AddJetContainer(
      jettype,
      AliJetContainer::antikt_algorithm,
      recombinationScheme,
      jetradius,
      ((jettype == AliJetContainer::kFullJet) || (jettype == AliJetContainer::kNeutralJet)) ? AliEmcalJet::kEMCALfid : AliEmcalJet::kTPC,
      particles, nullptr, partLevelTag);
  mcjets->SetJetPtCut(0.);
  mcjets->SetMaxTrackPt(1000.);
  mcjets->SetName("partjets");
  responsemaker->SetNamePartLevelJetContainer(mcjets->GetName());

  AliTrackContainer *tracks(nullptr);
  if ((jettype == AliJetContainer::kChargedJet) || (jettype == AliJetContainer::kFullJet))
  {
    tracks = responsemaker->AddTrackContainer(EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::TrackContainerNameFactory(isAOD));
    std::cout << "Track container name: " << tracks->GetName() << std::endl;
    //if embedding need to specify that the tracks are embedded
    if (ifembed)
      tracks->SetIsEmbedding(true);
    tracks->SetMinPt(0.15);
  }
  AliClusterContainer *clusters(nullptr);
  if ((jettype == AliJetContainer::kFullJet) || (jettype == AliJetContainer::kNeutralJet))
  {
    std::cout << "Using full or neutral jets ..." << std::endl;
    clusters = responsemaker->AddClusterContainer(EMCalTriggerPtAnalysis::AliEmcalAnalysisFactory::ClusterContainerNameFactory(isAOD));
    std::cout << "Cluster container name: " << clusters->GetName() << std::endl;
    clusters->SetClusHadCorrEnergyCut(0.3); // 300 MeV E-cut
    clusters->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
  }
  else
  {
    std::cout << "Using charged jets ... " << std::endl;
  }

  TString detLevelTag("Jet");
  //if embedding need to specify the tag of the detector (or hybrid) level jet from the embedding framework
  if (ifembed)
    detLevelTag = "hybridLevelJets";
  AliJetContainer *datajets = responsemaker->AddJetContainer(
      jettype,
      AliJetContainer::antikt_algorithm,
      recombinationScheme,
      jetradius,
      ((jettype == AliJetContainer::kFullJet) || (jettype == AliJetContainer::kNeutralJet)) ? AliEmcalJet::kEMCALfid : AliEmcalJet::kTPCfid,
      tracks, clusters, detLevelTag);
  datajets->SetJetPtCut(0.);
  datajets->SetMaxTrackPt(1000.);
  datajets->SetName("detjets");
  //if embedding then this jet is the unsubtracted jet and the subtracted jet is added in the run macro, if not then it is the detector level jet
  if (!ifembed)
    responsemaker->SetNameDetLevelJetContainer(datajets->GetName());
  else
    responsemaker->SetNameUnSubLevelJetContainer(datajets->GetName());

  std::string jettypestring;
  switch (jettype)
  {
  case AliJetContainer::kFullJet:
    jettypestring = "FullJets";
    break;
  case AliJetContainer::kChargedJet:
    jettypestring = "ChargedJets";
    break;
  case AliJetContainer::kNeutralJet:
    jettypestring = "NeutralJets";
    break;
  default:
    jettypestring = "Undef";
  };

  EBinningMode_t binmode(kSDModeINT7);
  std::string triggerstring(trigger);
  if (triggerstring == "EJ1")
    binmode = kSDModeEJ1;
  else if (triggerstring == "EJ2")
    binmode = kSDModeEJ2;
  responsemaker->SetBinningMode(binmode);

  // Connecting containers
  std::stringstream outputfile, histname;
  outputfile << mgr->GetCommonFileName() << ":SoftDropResponse_" << jettypestring << "_R" << std::setw(2) << std::setfill('0') << int(jetradius * 10.) << "_" << trigger;
  histname << "SoftDropResponseHistos_" << jettypestring << "_R" << std::setw(2) << std::setfill('0') << int(jetradius * 10.) << "_" << trigger;
  mgr->ConnectInput(responsemaker, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(responsemaker, 1, mgr->CreateContainer(histname.str().data(), AliEmcalList::Class(), AliAnalysisManager::kOutputContainer, outputfile.str().data()));

  return responsemaker;
}
