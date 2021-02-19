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
#include <sstream>

#include <TArrayD.h>
#include <TBinning.h>
#include <TCustomBinning.h>
#include <TH1.h>
#include <TH2.h>
#include <TLinearBinning.h>
#include <TLorentzVector.h>
#include <TRandom.h>


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
                                                                               AliAnalysisEmcalSoftdropHelperImpl(),
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
                                                                               fDropMass0Jets(false),
                                                                               fJetTypeOutliers(kOutlierPartJet),
                                                                               fNameMCParticles("mcparticles"),
                                                                               fSampleSplitter(nullptr),
                                                                               fSampleTrimmer(nullptr),
                                                                               fPartLevelPtBinning(nullptr),
                                                                               fDetLevelPtBinning(nullptr),
                                                                               fIsEmbeddedEvent(false),
                                                                               fRequirePartJetInAcceptance(false),
                                                                               fFillPlotsResiduals(true),
                                                                               fFillPlotsQAGeneral(true),
                                                                               fFillPlotsQAConstituents(true),
                                                                               fFillPlotsQAOutliers(true),
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
                                                                                               AliAnalysisEmcalSoftdropHelperImpl(),
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
                                                                                               fDropMass0Jets(false),
                                                                                               fJetTypeOutliers(kOutlierPartJet),
                                                                                               fNameMCParticles("mcparticles"),
                                                                                               fSampleSplitter(nullptr),
                                                                                               fSampleTrimmer(nullptr),
                                                                                               fPartLevelPtBinning(nullptr),
                                                                                               fDetLevelPtBinning(nullptr),
                                                                                               fIsEmbeddedEvent(false),
                                                                                               fRequirePartJetInAcceptance(false),
                                                                                               fFillPlotsResiduals(true),
                                                                                               fFillPlotsQAGeneral(true),
                                                                                               fFillPlotsQAConstituents(true),
                                                                                               fFillPlotsQAOutliers(true),
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
  double R = double(int(GetJetContainer(fNameDetLevelJetContainer.Data())->GetJetRadius() * 1000.))/1000.;  // Save cast from float to double truncating after 3rd decimal digit

  fSampleSplitter = new TRandom;
  if (fSampleFraction < 1.)
    fSampleTrimmer = new TRandom;

  if (!fPartLevelPtBinning)
    fPartLevelPtBinning = GetDefaultPartLevelPtBinning(fBinningMode);
  if (!fDetLevelPtBinning)
    fDetLevelPtBinning = GetDefaultDetLevelPtBinning(fBinningMode);
  std::unique_ptr<TBinning> zgbinning(GetZgBinning(fZcut)),
                            rgbinning(GetRgBinning(R)),
                            nsdbinning(new TLinearBinning(22, -1.5, 20.5)),       // Negative bins are for untagged jets
                            thetagbinning(new TLinearBinning(11, -0.1, 1.)),
                            ptbinningFine(new TLinearBinning(500, 0., 500.)),
                            residualsbinning(new TLinearBinning(1000, -10., 10.)),
                            rgbinningtruefine(new TLinearBinning(100, 0., 1.)),
                            tagbinning(new TLinearBinning(4, -0.5, 3.5));
  TArrayD binEdgesZg, binEdgesRg, binEdgesNsd, binEdgesThetag, binEdgesPtPart, binEdgesPtDet, binEdgesPtFine, binEdgesTag;
  zgbinning->CreateBinEdges(binEdgesZg);
  rgbinning->CreateBinEdges(binEdgesRg);
  nsdbinning->CreateBinEdges(binEdgesNsd);
  thetagbinning->CreateBinEdges(binEdgesThetag);
  ptbinningFine->CreateBinEdges(binEdgesPtFine);
  tagbinning->CreateBinEdges(binEdgesTag);
  fPartLevelPtBinning->CreateBinEdges(binEdgesPtPart);
  fDetLevelPtBinning->CreateBinEdges(binEdgesPtDet);

  const TBinning *sparsebinningZg[4] = {zgbinning.get(), ptbinningFine.get(), zgbinning.get(), ptbinningFine.get()},
                 *sparsebinningRg[4] = {rgbinning.get(), ptbinningFine.get(), rgbinning.get(), ptbinningFine.get()},
                 *sparsebinningNsd[4] = {nsdbinning.get(), ptbinningFine.get(), nsdbinning.get(), ptbinningFine.get()},
                 *sparsebinningThetag[4] = {thetagbinning.get(), ptbinningFine.get(), thetagbinning.get(), ptbinningFine.get()},
                 *rgsparsebinning[3] = {ptbinningFine.get(), rgbinningtruefine.get(), residualsbinning.get()},
                 *thetagsparsebinning[3] = {ptbinningFine.get(), rgbinningtruefine.get(), residualsbinning.get()},
                 *sparsebinningJEffPureZg[3] = {zgbinning.get(), ptbinningFine.get(), tagbinning.get()},
                 *sparsebinningJEffPureRg[3] = {rgbinning.get(), ptbinningFine.get(), tagbinning.get()},
                 *sparsebinningJEffPureThetag[3] = {thetagbinning.get(), ptbinningFine.get(), tagbinning.get()},
                 *sparsebinningJEffPureNsd[3] = {nsdbinning.get(), ptbinningFine.get(), tagbinning.get()};

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
        // Part. level distributions in the fiducial acceptance
        fHistManager.CreateTH2(Form("hZgPartLevelFine_%d", cent), Form("Zg dist at particle level_%d", cent), binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtFine.GetSize() - 1, binEdgesPtFine.GetArray());
        fHistManager.CreateTH2(Form("hRgPartLevelFine_%d", cent), Form("Rg dist at particle level_%d", cent), binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtFine.GetSize() - 1, binEdgesPtFine.GetArray());
        fHistManager.CreateTH2(Form("hNsdPartLevelFine_%d", cent), Form("Nsd dist at particle level_%d", cent), binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtFine.GetSize() - 1, binEdgesPtFine.GetArray());
        fHistManager.CreateTH2(Form("hThetagPartLevelFine_%d", cent), Form("Thetag dist at particle level_%d", cent), binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtFine.GetSize() - 1, binEdgesPtFine.GetArray());

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

      // Jet finding efficiency
      fHistManager.CreateTHnSparse(Form("hZgJetFindingEfficiency_%d", cent), Form("Jet finding efficiency as function of Zg and pt, %d centrality bin", cent), 3, sparsebinningJEffPureZg);
      fHistManager.CreateTHnSparse(Form("hRgJetFindingEfficiency_%d", cent), Form("Jet finding efficiency as function of Rg and pt, %d centrality bin", cent), 3, sparsebinningJEffPureRg);
      fHistManager.CreateTHnSparse(Form("hThetagJetFindingEfficiency_%d", cent), Form("Jet finding efficiency as function of Thetag and pt, %d centrality bin", cent), 3, sparsebinningJEffPureThetag);
      fHistManager.CreateTHnSparse(Form("hNsdJetFindingEfficiency_%d", cent), Form("Jet finding efficiency as function of Nsd and pt, %d centrality bin", cent), 3, sparsebinningJEffPureNsd);
      // split in response and no response sample for closure test
      fHistManager.CreateTHnSparse(Form("hZgJetFindingEfficiencyClosureResp_%d", cent), Form("Jet finding efficiency as function of Zg and pt (for closure test, response sample), %d centrality bin", cent), 3, sparsebinningJEffPureZg);
      fHistManager.CreateTHnSparse(Form("hRgJetFindingEfficiencyClosureResp_%d", cent), Form("Jet finding efficiency as function of Rg and pt (for closure test, response sample), %d centrality bin", cent), 3, sparsebinningJEffPureRg);
      fHistManager.CreateTHnSparse(Form("hThetagJetFindingEfficiencyClosureResp_%d", cent), Form("Jet finding efficiency as function of Thetag and pt (for closure test, response sample), %d centrality bin", cent), 3, sparsebinningJEffPureThetag);
      fHistManager.CreateTHnSparse(Form("hNsdJetFindingEfficiencyClosureResp_%d", cent), Form("Jet finding efficiency as function of Nsd and pt (for closure test, response sample), %d centrality bin", cent), 3, sparsebinningJEffPureNsd);
      fHistManager.CreateTHnSparse(Form("hZgJetFindingEfficiencyClosureNoResp_%d", cent), Form("Jet finding efficiency as function of Zg and pt (for closure test, no-response sample), %d centrality bin", cent), 3, sparsebinningJEffPureZg);
      fHistManager.CreateTHnSparse(Form("hRgJetFindingEfficiencyClosureNoResp_%d", cent), Form("Jet finding efficiency as function of Rg and pt (for closure test, no-response sample), %d centrality bin", cent), 3, sparsebinningJEffPureRg);
      fHistManager.CreateTHnSparse(Form("hThetagJetFindingEfficiencyClosureNoResp_%d", cent), Form("Jet finding efficiency as function of Thetag and pt (for closure test, no-response sample), %d centrality bin", cent), 3, sparsebinningJEffPureThetag);
      fHistManager.CreateTHnSparse(Form("hNsdJetFindingEfficiencyClosureNoResp_%d", cent), Form("Jet finding efficiency as function of Nsd and pt (for closure test, no-response sample), %d centrality bin", cent), 3, sparsebinningJEffPureNsd);
      // Jet finding purity
      fHistManager.CreateTHnSparse(Form("hZgJetFindingPurity_%d", cent), Form("Jet finding efficiency as function of Zg and pt, %d centrality bin", cent), 3, sparsebinningJEffPureZg);
      fHistManager.CreateTHnSparse(Form("hRgJetFindingPurity_%d", cent), Form("Jet finding efficiency as function of Rg and pt, %d centrality bin", cent), 3, sparsebinningJEffPureRg);
      fHistManager.CreateTHnSparse(Form("hThetagJetFindingPurity_%d", cent), Form("Jet finding efficiency as function of Thetag and pt, %d centrality bin", cent), 3, sparsebinningJEffPureThetag);
      fHistManager.CreateTHnSparse(Form("hNsdJetFindingPurity_%d", cent), Form("Jet finding efficiency as function of Nsd and pt, %d centrality bin", cent), 3, sparsebinningJEffPureNsd);
      // split in response and no response sample for closure test
      fHistManager.CreateTHnSparse(Form("hZgJetFindingPurityClosureResp_%d", cent), Form("Jet finding efficiency as function of Zg and pt (for closure test, response sample), %d centrality bin", cent), 3, sparsebinningJEffPureZg);
      fHistManager.CreateTHnSparse(Form("hRgJetFindingPurityClosureResp_%d", cent), Form("Jet finding efficiency as function of Rg and pt (for closure test, response sample), %d centrality bin", cent), 3, sparsebinningJEffPureRg);
      fHistManager.CreateTHnSparse(Form("hThetagJetFindingPurityClosureResp_%d", cent), Form("Jet finding efficiency as function of Thetag and pt (for closure test, response sample), %d centrality bin", cent), 3, sparsebinningJEffPureThetag);
      fHistManager.CreateTHnSparse(Form("hNsdJetFindingPurityClosureResp_%d", cent), Form("Jet finding efficiency as function of Nsd and pt (for closure test, response sample), %d centrality bin", cent), 3, sparsebinningJEffPureNsd);
      fHistManager.CreateTHnSparse(Form("hZgJetFindingPurityClosureNoResp_%d", cent), Form("Jet finding efficiency as function of Zg and pt (for closure test, no-response sample), %d centrality bin", cent), 3, sparsebinningJEffPureZg);
      fHistManager.CreateTHnSparse(Form("hRgJetFindingPurityClosureNoResp_%d", cent), Form("Jet finding efficiency as function of Rg and pt (for closure test, no-response sample), %d centrality bin", cent), 3, sparsebinningJEffPureRg);
      fHistManager.CreateTHnSparse(Form("hThetagJetFindingPurityClosureNoResp_%d", cent), Form("Jet finding efficiency as function of Thetag and pt (for closure test, no-response sample), %d centrality bin", cent), 3, sparsebinningJEffPureThetag);
      fHistManager.CreateTHnSparse(Form("hNsdJetFindingPurityClosureNoResp_%d", cent), Form("Jet finding efficiency as function of Nsd and pt (for closure test, no-response sample), %d centrality bin", cent), 3, sparsebinningJEffPureNsd);
      
      // Closure test (before part/det matching requirement)
      fHistManager.CreateTH2(Form("hZgClosureNoRespDetAllFine_%d", cent), Form("Zg at det. level of all jets for closure test (no-response sample), %d centrality bin", cent), binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
      fHistManager.CreateTH2(Form("hZgClosureNoRespPartAllFine_%d", cent), Form("Zg at part. level of all jets for closure test (no-response sample), %d centrality bin", cent), binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
      fHistManager.CreateTH2(Form("hZgClosureRespDetAllFine_%d", cent), Form("Zg at det. level of all jets for closure test (response sample), %d centrality bin", cent), binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
      fHistManager.CreateTH2(Form("hZgClosureRespPartAllFine_%d", cent), Form("Zg at part. level of all jets for closure test (response sample), %d centrality bin", cent), binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
      fHistManager.CreateTH2(Form("hRgClosureNoRespDetAllFine_%d", cent), Form("Rg at det. level of all jets for closure test (no-response sample), %d centrality bin", cent), binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
      fHistManager.CreateTH2(Form("hRgClosureNoRespPartAllFine_%d", cent), Form("Rg at part. level of all jets for closure test (no-response sample), %d centrality bin", cent), binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
      fHistManager.CreateTH2(Form("hRgClosureRespDetAllFine_%d", cent), Form("Rg at det. level of all jets for closure test (response sample), %d centrality bin", cent), binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
      fHistManager.CreateTH2(Form("hRgClosureRespPartAllFine_%d", cent), Form("Rg at part. level of all jets for closure test (response sample), %d centrality bin", cent), binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
      fHistManager.CreateTH2(Form("hThetagClosureNoRespDetAllFine_%d", cent), Form("Thetag at det. level of all jets for closure test (no-response sample), %d centrality bin", cent), binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
      fHistManager.CreateTH2(Form("hThetagClosureNoRespPartAllFine_%d", cent), Form("Thetag at part. level of all jets for closure test (no-response sample), %d centrality bin", cent), binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
      fHistManager.CreateTH2(Form("hThetagClosureRespDetAllFine_%d", cent), Form("Thetag at det. level of all jets for closure test (response sample), %d centrality bin", cent), binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
      fHistManager.CreateTH2(Form("hThetagClosureRespPartAllFine_%d", cent), Form("Thetag at part. level of all jets for closure test (response sample), %d centrality bin", cent), binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
      fHistManager.CreateTH2(Form("hNsdClosureNoRespDetAllFine_%d", cent), Form("Nsd at det. level of all jets for closure test (no-response sample), %d centrality bin", cent), binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
      fHistManager.CreateTH2(Form("hNsdClosureNoRespPartAllFine_%d", cent), Form("Nsd at part. level of all jets for closure test (no-response sample), %d centrality bin", cent), binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
      fHistManager.CreateTH2(Form("hNsdClosureRespDetAllFine_%d", cent), Form("Nsd at det. level of all jets for closure test (response sample), %d centrality bin", cent), binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
      fHistManager.CreateTH2(Form("hNsdClosureRespPartAllFine_%d", cent), Form("Nsd at part. level of all jets for closure test (response sample), %d centrality bin", cent), binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
      fHistManager.CreateTH2(Form("hZgClosureTruthNoRespPartAllFine_%d", cent), Form("Zg truth at part. level of all jets for closure test (no-response sample), %d centrality bin", cent), binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
      fHistManager.CreateTH2(Form("hZgClosureTruthRespPartAllFine_%d", cent), Form("Zg truth at part. level of all jets for closure test (response sample), %d centrality bin", cent), binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
      fHistManager.CreateTH2(Form("hRgClosureTruthNoRespPartAllFine_%d", cent), Form("Rg truth at part. level of all jets for closure test (no-response sample), %d centrality bin", cent), binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
      fHistManager.CreateTH2(Form("hRgClosureTruthRespPartAllFine_%d", cent), Form("Rg truth at part. level of all jets for closure test (response sample), %d centrality bin", cent), binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
      fHistManager.CreateTH2(Form("hThetagClosureTruthNoRespPartAllFine_%d", cent), Form("Thetag truth at part. level of all jets for closure test (no-response sample), %d centrality bin", cent), binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
      fHistManager.CreateTH2(Form("hThetagClosureTruthRespPartAllFine_%d", cent), Form("Thetag truth at part. level of all jets for closure test (response sample), %d centrality bin", cent), binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
      fHistManager.CreateTH2(Form("hNsdClosureTruthNoRespPartAllFine_%d", cent), Form("Nsd truth at part. level of all jets for closure test (no-response sample), %d centrality bin", cent), binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
      fHistManager.CreateTH2(Form("hNsdClosureTruthRespPartAllFine_%d", cent), Form("Nsd truth at part. level of all jets for closure test (response sample), %d centrality bin", cent), binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());

      if(fFillPlotsResiduals) {
        // Residuals vs. pt,part
        fHistManager.CreateTH2(Form("hZgResiduals_%d", cent), Form("z_{g} residuals (%d centrality bin); p_{t,part} (GeV/c); z_{g, det} - z_{g, part}", cent), 350, 0., 350., 100, -1., 1.);
        fHistManager.CreateTH2(Form("hRgResiduals_%d", cent), Form("R_{g} residuals (%d centrality bin); p_{t,part} (GeV/c); R_{g, det} - R_{g, part}", cent), 350, 0., 350., 100, -1., 1.);
        fHistManager.CreateTH2(Form("hThetagResiduals_%d", cent), Form("#Theta_{g} residuals (%d centrality bin); p_{t,part} (GeV/c); #Theta_{g, det} - #Theta_{g, part}", cent), 350, 0., 350., 200, -2., 2.);
        fHistManager.CreateTH2(Form("hNsdResiduals_%d", cent), Form("n_{SD} residuals (%d centrality bin); p_{t,part} (GeV/c); n_{SD, det} - n_{SD, part}", cent), 350, 0., 350., 80, -40., 40.);
        fHistManager.CreateTH2(Form("hZgResidualsNormalized_%d", cent), Form("z_{g} residuals (normalized, %d centrality bin); p_{t,part} (GeV/c); (z_{g, det} - z_{g, part})/z_{g, part}", cent), 350, 0., 350., 100, -1., 1.);
        fHistManager.CreateTH2(Form("hRgResidualsNormalized_%d", cent), Form("R_{g} residuals (normalized, %d centrality bin); p_{t,part} (GeV/c); (R_{g, det} - R_{g, part})/R_{g, part}", cent), 350, 0., 350., 100, -1., 1.);
        fHistManager.CreateTH2(Form("hThetagResidualsNormalized_%d", cent), Form("#Theta_{g} residuals (normalized, %d centrality bin); p_{t,part} (GeV/c); (#Theta_{g, det} - #Theta_{g, part})/#Theta_{g, part}", cent), 350, 0., 350., 100, -1., 1.);
        fHistManager.CreateTH2(Form("hNsdResidualsNormalized_%d", cent), Form("n_{SD} residuals (normalized, %d centrality bin); p_{t,part} (GeV/c); (n_{SD, det} - n_{SD, part})/n_{SD, part}", cent), 350, 0., 350., 100, -10., 10.);
        // Residuals vs. Rg/Theatg
        fHistManager.CreateTHnSparse(Form("hResidualsRg_%d", cent), Form("ResidualsRg (%d centrality bin)", cent), 3, rgsparsebinning);
        fHistManager.CreateTHnSparse(Form("hResidualsRgNormalized_%d", cent), Form("ResidualsRg (%d centrality bin)", cent), 3, rgsparsebinning);
        fHistManager.CreateTHnSparse(Form("hResidualsThetag_%d", cent), Form("ResidualsThetag (%d centrality bin)", cent), 3, thetagsparsebinning);
        fHistManager.CreateTHnSparse(Form("hResidualsThetagNormalized_%d", cent), Form("ResidualThetaRg (%d centrality bin)", cent), 3, thetagsparsebinning);
      }

      // Jets failed in soft drop
      fHistManager.CreateTH1(Form("hSkippedJetsPart_%d", cent), Form("Skipped jets at part. level, %d centrality bin", cent), 350, 0., 350.);
      fHistManager.CreateTH2(Form("hSkippedJetsPartCorr_%d", cent), Form("Pt of pairs with failed SD at part. level, %d centrality bin", cent), binEdgesPtFine, binEdgesPtFine); 
      fHistManager.CreateTH1(Form("hSkippedJetsDet_%d", cent), Form("Skipped jets at det. level, %d centrality bin", cent), 350, 0., 350.);
      fHistManager.CreateTH2(Form("hSkippedJetsDetCorr_%d", cent), Form("Pt of pairs with failed SD at det. level, %d centrality bin", cent), binEdgesPtFine, binEdgesPtFine); 
      if(fFillPlotsQAGeneral) {
        // a bit of QA stuff
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
      }
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
      // Part. level distributions in the fiducial acceptance
      fHistManager.CreateTH2("hZgPartLevelFine", "Zg dist at particle level", binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtFine.GetSize() - 1, binEdgesPtFine.GetArray());
      fHistManager.CreateTH2("hRgPartLevelFine", "Rg dist at particle level", binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtFine.GetSize() - 1, binEdgesPtFine.GetArray());
      fHistManager.CreateTH2("hNsdPartLevelFine", "Nsd dist at particle level", binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtFine.GetSize() - 1, binEdgesPtFine.GetArray());
      fHistManager.CreateTH2("hThetagPartLevelFine", "Thetag dist at particle level", binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtFine.GetSize() - 1, binEdgesPtFine.GetArray());

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

    // Jet finding efficiency
    fHistManager.CreateTHnSparse("hZgJetFindingEfficiency", "Jet finding efficiency as function of Zg and pt", 3, sparsebinningJEffPureZg);
    fHistManager.CreateTHnSparse("hRgJetFindingEfficiency", "Jet finding efficiency as function of Rg and pt", 3, sparsebinningJEffPureRg);
    fHistManager.CreateTHnSparse("hThetagJetFindingEfficiency", "Jet finding efficiency as function of Thetag and pt", 3, sparsebinningJEffPureThetag);
    fHistManager.CreateTHnSparse("hNsdJetFindingEfficiency", "Jet finding efficiency as function of Nsd and pt", 3, sparsebinningJEffPureNsd);
    // split in response and no response sample for closure test
    fHistManager.CreateTHnSparse("hZgJetFindingEfficiencyClosureResp", "Jet finding efficiency as function of Zg and pt (response sample, for closure test)", 3, sparsebinningJEffPureZg);
    fHistManager.CreateTHnSparse("hRgJetFindingEfficiencyClosureResp", "Jet finding efficiency as function of Rg and pt (response sample, for closure test)", 3, sparsebinningJEffPureRg);
    fHistManager.CreateTHnSparse("hThetagJetFindingEfficiencyClosureResp", "Jet finding efficiency as function of Thetag and pt (response sample, for closure test)", 3, sparsebinningJEffPureThetag);
    fHistManager.CreateTHnSparse("hNsdJetFindingEfficiencyClosureResp", "Jet finding efficiency as function of Nsd and pt (response sample, for closure test)", 3, sparsebinningJEffPureNsd);
    fHistManager.CreateTHnSparse("hZgJetFindingEfficiencyClosureNoResp", "Jet finding efficiency as function of Zg and pt (no-response sample, for closure test)", 3, sparsebinningJEffPureZg);
    fHistManager.CreateTHnSparse("hRgJetFindingEfficiencyClosureNoResp", "Jet finding efficiency as function of Rg and pt (no-response sample, for closure test)", 3, sparsebinningJEffPureRg);
    fHistManager.CreateTHnSparse("hThetagJetFindingEfficiencyClosureNoResp", "Jet finding efficiency as function of Thetag and pt (no-response sample, for closure test)", 3, sparsebinningJEffPureThetag);
    fHistManager.CreateTHnSparse("hNsdJetFindingEfficiencyClosureNoResp", "Jet finding efficiency as function of Nsd and pt (no-response sample, for closure test)", 3, sparsebinningJEffPureNsd);
    // Jet finding purity
    fHistManager.CreateTHnSparse("hZgJetFindingPurity", "Jet finding efficiency as function of Zg and pt", 3, sparsebinningJEffPureZg);
    fHistManager.CreateTHnSparse("hRgJetFindingPurity", "Jet finding efficiency as function of Rg and pt", 3, sparsebinningJEffPureRg);
    fHistManager.CreateTHnSparse("hThetagJetFindingPurity", "Jet finding efficiency as function of Thetag and pt", 3, sparsebinningJEffPureThetag);
    fHistManager.CreateTHnSparse("hNsdJetFindingPurity", "Jet finding efficiency as function of Nsd and pt", 3, sparsebinningJEffPureNsd);
    // split in response and no response sample for closure test
    fHistManager.CreateTHnSparse("hZgJetFindingPurityClosureResp", "Jet finding efficiency as function of Zg and pt (response sample, for closure test)", 3, sparsebinningJEffPureZg);
    fHistManager.CreateTHnSparse("hRgJetFindingPurityClosureResp", "Jet finding efficiency as function of Rg and pt (response sample, for closure test)", 3, sparsebinningJEffPureRg);
    fHistManager.CreateTHnSparse("hThetagJetFindingPurityClosureResp", "Jet finding efficiency as function of Thetag and pt (response sample, for closure test)", 3, sparsebinningJEffPureThetag);
    fHistManager.CreateTHnSparse("hNsdJetFindingPurityClosureResp", "Jet finding efficiency as function of Nsd and pt (response sample, for closure test)", 3, sparsebinningJEffPureNsd);
    fHistManager.CreateTHnSparse("hZgJetFindingPurityClosureNoResp", "Jet finding efficiency as function of Zg and pt (no-response sample, for closure test)", 3, sparsebinningJEffPureZg);
    fHistManager.CreateTHnSparse("hRgJetFindingPurityClosureNoResp", "Jet finding efficiency as function of Rg and pt (no-response sample, for closure test)", 3, sparsebinningJEffPureRg);
    fHistManager.CreateTHnSparse("hThetagJetFindingPurityClosureNoResp", "Jet finding efficiency as function of Thetag and pt (no-response sample, for closure test)", 3, sparsebinningJEffPureThetag);
    fHistManager.CreateTHnSparse("hNsdJetFindingPurityClosureNoResp", "Jet finding efficiency as function of Nsd and pt (no-response sample, for closure test)", 3, sparsebinningJEffPureNsd);

    // Closure test (before part/det matching requirement)
    fHistManager.CreateTH2("hZgClosureNoRespDetAllFine", "Zg at det. level of all jets for closure test (no-response sample)", binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
    fHistManager.CreateTH2("hZgClosureNoRespPartAllFine", "Zg at part. level of all jets for closure test (no-response sample)", binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
    fHistManager.CreateTH2("hZgClosureRespDetAllFine", "Zg at det. level of all jets for closure test (response sample)", binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
    fHistManager.CreateTH2("hZgClosureRespPartAllFine", "Zg at part. level of all jets for closure test (response sample)", binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
    fHistManager.CreateTH2("hRgClosureNoRespDetAllFine", "Rg at det. level of all jets for closure test (no-response sample)", binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
    fHistManager.CreateTH2("hRgClosureNoRespPartAllFine", "Rg at part. level of all jets for closure test (no-response sample)", binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
    fHistManager.CreateTH2("hRgClosureRespDetAllFine", "Rg at det. level of all jets for closure test (response sample)", binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
    fHistManager.CreateTH2("hRgClosureRespPartAllFine", "Rg at part. level of all jets for closure test (response sample)", binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
    fHistManager.CreateTH2("hThetagClosureNoRespDetAllFine", "Thetag at det. level of all jets for closure test (no-response sample)", binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
    fHistManager.CreateTH2("hThetagClosureNoRespPartAllFine", "Thetag at part. level of all jets for closure test (no-response sample)", binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
    fHistManager.CreateTH2("hThetagClosureRespDetAllFine", "Thetag at det. level of all jets for closure test (response sample)", binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
    fHistManager.CreateTH2("hThetagClosureRespPartAllFine", "Thetag at part. level of all jets for closure test (response sample)", binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
    fHistManager.CreateTH2("hNsdClosureNoRespDetAllFine", "Nsd at det. level of all jets for closure test (no-response sample)", binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
    fHistManager.CreateTH2("hNsdClosureNoRespPartAllFine", "Nsd at part. level of all jets for closure test (no-response sample)", binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
    fHistManager.CreateTH2("hNsdClosureRespDetAllFine", "Nsd at det. level of all jets for closure test (response sample)", binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
    fHistManager.CreateTH2("hNsdClosureRespPartAllFine", "Nsd at part. level of all jets for closure test (response sample)", binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
    fHistManager.CreateTH2("hZgClosureTruthNoRespPartAllFine", "Zg truth at part. level of all jets for closure test (no-response sample)", binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
    fHistManager.CreateTH2("hZgClosureTruthRespPartAllFine", "Zg truth at part. level of all jets for closure test (response sample)", binEdgesZg.GetSize() - 1, binEdgesZg.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
    fHistManager.CreateTH2("hRgClosureTruthNoRespPartAllFine", "Rg truth at part. level of all jets for closure test (no-response sample)", binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
    fHistManager.CreateTH2("hRgClosureTruthRespPartAllFine", "Rg truth at part. level of all jets for closure test (response sample)", binEdgesRg.GetSize() - 1, binEdgesRg.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
    fHistManager.CreateTH2("hThetagClosureTruthNoRespPartAllFine", "Thetag truth at part. level of all jets for closure test (no-response sample)", binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
    fHistManager.CreateTH2("hThetagClosureTruthRespPartAllFine", "Thetag truth at part. level of all jets for closure test (response sample)", binEdgesThetag.GetSize() - 1, binEdgesThetag.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
    fHistManager.CreateTH2("hNsdClosureTruthNoRespPartAllFine", "Nsd truth at part. level of all jets for closure test (no-response sample)", binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());
    fHistManager.CreateTH2("hNsdClosureTruthRespPartAllFine", "Nsd truth at part. level of all jets for closure test (response sample)", binEdgesNsd.GetSize() - 1, binEdgesNsd.GetArray(), binEdgesPtFine.GetSize() -1 , binEdgesPtFine.GetArray());

    if(fFillPlotsResiduals){
      // Residuals vs. pt,part
      fHistManager.CreateTH2("hZgResiduals", "z_{g} residuals vs. p_{t,part}; p_{t,part} (GeV/c); z_{g, det} - z_{g, part}", 350, 0., 350., 100, -1., 1.);
      fHistManager.CreateTH2("hRgResiduals", "R_{g} residuals vs. p_{t,part}; p_{t,part} (GeV/c); R_{g, det} - R_{g, part}", 350, 0., 350., 100, -1., 1.);
      fHistManager.CreateTH2("hThetagResiduals", "#Theta_{g} residuals vs. p_{t,part}; p_{t,part} (GeV/c); #Theta_{g, det} - #Theta_{g, part}", 350, 0., 350., 200, -2., 2.);
      fHistManager.CreateTH2("hNsdResiduals", "n_{SD} residuals vs. p_{t,part}; p_{t,part} (GeV/c); n_{SD, det} - n_{SD, part}", 350, 0., 350., 80, -40., 40.);
      fHistManager.CreateTH2("hZgResidualsNormalized", "z_{g} residuals (normalized) vs. p_{t,part}; p_{t,part} (GeV/c); (z_{g, det} - z_{g, part})/z_{g, part}", 350, 0., 350., 100, -1., 1.);
      fHistManager.CreateTH2("hRgResidualsNormalized", "R_{g} residuals (normalized) vs. p_{t,part}; p_{t,part} (GeV/c); (R_{g, det} - R_{g, part})/R_{g, part}", 350, 0., 350., 100, -1., 1.);
      fHistManager.CreateTH2("hThetagResidualsNormalized", "#Theta_{g} residuals (normalized) vs. p_{t,part}; p_{t,part} (GeV/c); (#Theta_{g, det} - #Theta_{g, part})/#Theta_{g, part}", 350, 0., 350., 100, -1., 1.);
      fHistManager.CreateTH2("hNsdResidualsNormalized", "n_{SD} residuals (normalized) vs. p_{t,part}; p_{t,part} (GeV/c); (n_{SD, det} - n_{SD, part})/n_{SD, part}", 350, 0., 350., 100, -10., 10.);
      // Residuals vs. Rg/Theatg
      fHistManager.CreateTHnSparse("hResidualsRg", "ResidualsRg", 3, rgsparsebinning);
      fHistManager.CreateTHnSparse("hResidualsRgNormalized", "ResidualsRg", 3, rgsparsebinning);
      fHistManager.CreateTHnSparse("hResidualsThetag", "ResidualsThetag", 3, thetagsparsebinning);
      fHistManager.CreateTHnSparse("hResidualsThetagNormalized", "ResidualsThetag", 3, thetagsparsebinning);
      fHistManager.CreateTH2("hResualsRgLowRg", "Rg residuals (Rg < 0.02)", 350, 0., 350, 200, -100, 100);
    }

    // Jets failed in soft drop
    fHistManager.CreateTH1("hSkippedJetsPart", "Skipped jets at part. level", 350, 0., 350.);
    fHistManager.CreateTH2("hSkippedJetsPartCorr", "Pt of pairs with failed SD at part. level", binEdgesPtFine, binEdgesPtFine); 
    fHistManager.CreateTH1("hSkippedJetsDet", "Skipped jets at det. level", 350, 0., 350.);
    fHistManager.CreateTH2("hSkippedJetsDetCorr", "Pt of pairs with failed SD at det. level", binEdgesPtFine, binEdgesPtFine); 
    if(fFillPlotsQAGeneral) {
      // a bit of QA stuff
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
  }

  if(fFillPlotsQAConstituents) {
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
  }
  if(fFillPlotsQAOutliers){
    fHistManager.CreateTH1("hFracPtHardPart", "Part. level jet Pt relative to the Pt-hard of the event", 100, 0., 10.);
    fHistManager.CreateTH1("hFracPtHardDet", "Det. level jet Pt relative to the Pt-hard of the event", 100, 0., 10.);
  }

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
  if (!(fIsPythia || fIsHerwig || fIsHepMC))
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
  enum TagStatus_t {
    kNoMatchedJet = 0,
    kMatchedJetNoSoftDrop = 1,
    kMatchedJetNoAcceptance = 2,
    kPairAccepted = 3
  };
  AliJetContainer *partLevelJets = this->GetJetContainer(fNamePartLevelJetContainer),
                  *detLevelJets = GetJetContainer(fNameDetLevelJetContainer);
  AliClusterContainer *clusters = GetClusterContainer(AliEmcalAnalysisFactory::ClusterContainerNameFactory(fInputEvent->IsA() == AliAODEvent::Class()));
  AliTrackContainer *tracks = GetTrackContainer(AliEmcalAnalysisFactory::TrackContainerNameFactory(fInputEvent->IsA() == AliAODEvent::Class()));
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

  // sample splitting (for closure test)
  // must be done at event level in order to have statistically independent
  // samples for both effciency and purity, which are different loops
  bool closureUseResponse = (fSampleSplitter->Uniform() < fFractionResponseClosure);

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

  SoftdropParams sdsettings{fReclusterizer, fBeta, fZcut,fUseChargedConstituents, fUseNeutralConstituents};
  for (auto detjet : detLevelJets->accepted()) {
    AliEmcalJet *partjet = nullptr,
                *jetUS = nullptr;  //variables for embedded pbpb data
    Int_t ilab = -1;

    //for embedding, find the unsubtracted jet and get it's matched detector level jet
    if (fIsEmbeddedEvent) {
      for (Int_t i = 0; i < jetContUS->GetNJets(); i++) {
        jetUS = jetContUS->GetJet(i);
        if (jetUS->GetLabel() == detjet->GetLabel()) {
          ilab = i;
          break;
        }
      }
      if (ilab == -1)
        continue;
      jetUS = jetContUS->GetJet(ilab);
      partjet = jetUS->ClosestJet();
    } else {
      //if we aren't embedding then just find the matched jet
      partjet = detjet->ClosestJet();
    }

    //cut on the shared pt fraction, when embedding the unsubtracted jet should be used
    Double_t fraction = 0;
    if (fIsEmbeddedEvent)
      fraction = jetContUS->GetFractionSharedPt(jetUS);
    else
      fraction = detLevelJets->GetFractionSharedPt(detjet);
    if (fraction < fMinFractionShared){
      continue;
    }

    Bool_t hasSomeError = false;
    TagStatus_t tagStatus = kNoMatchedJet;
    if (partjet) {
      //one extra level of matching needed for embedding to go from detector to particle level
      if (fIsEmbeddedEvent) {
        partjet = partjet->ClosestJet();
        if(partjet)
          tagStatus = kMatchedJetNoAcceptance;
      }
      if(partjet) {
        tagStatus = kMatchedJetNoAcceptance;
      }
    } else {
      hasSomeError = true;
      AliDebugStream(1) << "No part jet" << std::endl;
    }

    SoftdropResults softdropDet, softdropPart;
    std::vector<SoftdropResults> splittingsDet, splittingsPart;
    try {
      softdropDet = MakeSoftdrop(*detjet, detLevelJets->GetJetRadius(),false, sdsettings, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy(), fVertex, fDropMass0Jets);
      splittingsDet = IterativeDecluster(*detjet, detLevelJets->GetJetRadius(), false, sdsettings, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy(), fVertex, fDropMass0Jets);
    } catch(...) {
      // Failed SoftDrop for det. level jet.
      if(fForceBeamType != kpp) {
        fHistManager.FillTH1(Form("hSkippedJetsDet_%d", fCentBin), detjet->Pt());
        if(partjet) fHistManager.FillTH2(Form("hSkippedJetsDetCorr_%d", fCentBin), partjet->Pt(), detjet->Pt());
      } else {
        fHistManager.FillTH1("hSkippedJetsDet", detjet->Pt());
        if(partjet) fHistManager.FillTH2("hSkippedJetsDetCorr", partjet->Pt(), detjet->Pt());
      }
      continue;
    }

    // Checked requiring the presence of a part. level jet
    // - Part. jet in same acceptance
    // - Part. jet fulfilling SoftDrop
    bool partJetInAcceptance = false,
         hasMCSoftDrop = false;
    if(partjet) {
      // Get the SoftDrop params for part. level jet
      // Handle failed SoftDrop for part. level jet as impurity
      // Check condition first in case the the same acceptance is not required for the analysis
      try {
        softdropPart = MakeSoftdrop(*partjet, partLevelJets->GetJetRadius(), true, sdsettings, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy(), fVertex, fDropMass0Jets);
        splittingsPart = IterativeDecluster(*partjet, partLevelJets->GetJetRadius(),  true, sdsettings, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy(), fVertex, fDropMass0Jets);
        hasMCSoftDrop = true;
        tagStatus = kMatchedJetNoAcceptance;

        if(partjet->GetJetAcceptanceType() & detLevelJets->GetAcceptanceType()) {
          tagStatus = kPairAccepted;
          partJetInAcceptance = true;
        } else {
          AliDebugStream(1) << "Part level jet not in acceptance" << std::endl;
          hasSomeError = true;
        }
      } catch (...) {
        AliDebugStream(1) << "Part. level jet failed soft drop" << std::endl;
        hasSomeError = true;
        // Failed SoftDrop for det. level jet.
        if(fForceBeamType != kpp) {
          fHistManager.FillTH1(Form("hSkippedJetsPart_%d", fCentBin), partjet->Pt());
          fHistManager.FillTH2(Form("hSkippedJetsPartCorr_%d", fCentBin), partjet->Pt(), detjet->Pt());
        } else {
          fHistManager.FillTH1("hSkippedJetsPart", partjet->Pt());
          fHistManager.FillTH2("hSkippedJetsPartCorr", partjet->Pt(), detjet->Pt());
        }
      }

    }

    // Point Purity
    // 0 - no matched jet
    // 1 - matched jet not in the same acceptance
    // 2 - matched jet in the same acceptance
    // 3 - matched jet in same acceptance and SoftDrop of part level jet successful
    bool untaggedPurity = softdropDet.fZg < fZcut;
    double zgdet = softdropDet.fZg,
           rgdet = untaggedPurity ? -0.01 : softdropDet.fRg,
           thetagdet = untaggedPurity ? -0.05 : softdropDet.fRg/Rjet,
           nsddet = untaggedPurity ? -1. : double(splittingsDet.size());
    double pointPurityZg[3] = {zgdet, detjet->Pt(), static_cast<double>(tagStatus)},
           pointPurityRg[3] = {rgdet, detjet->Pt(), static_cast<double>(tagStatus)},
           pointPurityThetag[3] = {thetagdet, detjet->Pt(), static_cast<double>(tagStatus)},
           pointPurityNsd[3] = {nsddet, detjet->Pt(), static_cast<double>(tagStatus)};
    if(hasSomeError) {
      AliDebugStream(1) << "Was in error state, tag status " << pointPurityZg[2] << std::endl;
    }

    // Fill purity after acceptance check
    if (fForceBeamType != kpp) {
      fHistManager.FillTHnSparse(Form("hZgJetFindingPurity_%d", fCentBin), pointPurityZg);
      fHistManager.FillTHnSparse(Form("hRgJetFindingPurity_%d", fCentBin), pointPurityRg);
      fHistManager.FillTHnSparse(Form("hThetagJetFindingPurity_%d", fCentBin), pointPurityThetag);
      fHistManager.FillTHnSparse(Form("hNsdJetFindingPurity_%d", fCentBin), pointPurityNsd);
      // Split into response and no-response sample for closure test
      if(closureUseResponse) {
        fHistManager.FillTHnSparse(Form("hZgJetFindingPurityClosureResp_%d", fCentBin), pointPurityZg);
        fHistManager.FillTHnSparse(Form("hRgJetFindingPurityClosureResp_%d", fCentBin), pointPurityRg);
        fHistManager.FillTHnSparse(Form("hThetagJetFindingPurityClosureResp_%d", fCentBin), pointPurityThetag);
        fHistManager.FillTHnSparse(Form("hNsdJetFindingPurityClosureResp_%d", fCentBin), pointPurityNsd);
      } else {
        fHistManager.FillTHnSparse(Form("hZgJetFindingPurityClosureNoResp_%d", fCentBin), pointPurityZg);
        fHistManager.FillTHnSparse(Form("hRgJetFindingPurityClosureNoResp_%d", fCentBin), pointPurityRg);
        fHistManager.FillTHnSparse(Form("hThetagJetFindingPurityClosureNoResp_%d", fCentBin), pointPurityThetag);
        fHistManager.FillTHnSparse(Form("hNsdJetFindingPurityClosureNoResp_%d", fCentBin), pointPurityNsd);
      }
    } else {
      fHistManager.FillTHnSparse("hZgJetFindingPurity", pointPurityZg);
      fHistManager.FillTHnSparse("hRgJetFindingPurity", pointPurityRg);
      fHistManager.FillTHnSparse("hThetagJetFindingPurity", pointPurityThetag);
      fHistManager.FillTHnSparse("hNsdJetFindingPurity", pointPurityNsd);
      // Split into response and no-response sample for closure test
      if(closureUseResponse) {
        fHistManager.FillTHnSparse("hZgJetFindingPurityClosureResp", pointPurityZg);
        fHistManager.FillTHnSparse("hRgJetFindingPurityClosureResp", pointPurityRg);
        fHistManager.FillTHnSparse("hThetagJetFindingPurityClosureResp", pointPurityThetag);
        fHistManager.FillTHnSparse("hNsdJetFindingPurityClosureResp", pointPurityNsd);
      } else {
        fHistManager.FillTHnSparse("hZgJetFindingPurityClosureNoResp", pointPurityZg);
        fHistManager.FillTHnSparse("hRgJetFindingPurityClosureNoResp", pointPurityRg);
        fHistManager.FillTHnSparse("hThetagJetFindingPurityClosureNoResp", pointPurityThetag);
        fHistManager.FillTHnSparse("hNsdJetFindingPurityClosureNoResp", pointPurityNsd);
      }
    }
      
    // Fill det level histogram (for full closure test)
    // Distribtions before purity subtraction needed in order to test
    // the full correction chain (including purity correction and efficiency)
    // correction
    if(fForceBeamType != kpp) {
      if(closureUseResponse){
        fHistManager.FillTH2(Form("hZgClosureRespDetAllFine_%d", fCentBin), zgdet, detjet->Pt());
        fHistManager.FillTH2(Form("hRgClosureRespDetAllFine_%d", fCentBin), rgdet, detjet->Pt());
        fHistManager.FillTH2(Form("hThetagClosureRespDetAllFine_%d", fCentBin), thetagdet, detjet->Pt());
        fHistManager.FillTH2(Form("hNsdClosureRespDetAllFine_%d", fCentBin), nsddet, detjet->Pt());
      } else {
        fHistManager.FillTH2(Form("hZgClosureNoRespDetAllFine_%d", fCentBin), zgdet, detjet->Pt());
        fHistManager.FillTH2(Form("hRgClosureNoRespDetAllFine_%d", fCentBin), rgdet, detjet->Pt());
        fHistManager.FillTH2(Form("hThetagClosureNoRespDetAllFine_%d", fCentBin), thetagdet, detjet->Pt());
        fHistManager.FillTH2(Form("hNsdClosureNoRespDetAllFine_%d", fCentBin), nsddet, detjet->Pt());
      }
    } else {
      if(closureUseResponse){
        fHistManager.FillTH2("hZgClosureRespDetAllFine", zgdet, detjet->Pt());
        fHistManager.FillTH2("hRgClosureRespDetAllFine", rgdet, detjet->Pt());
        fHistManager.FillTH2("hThetagClosureRespDetAllFine", thetagdet, detjet->Pt());
        fHistManager.FillTH2("hNsdClosureRespDetAllFine", nsddet, detjet->Pt());
      } else {
        fHistManager.FillTH2("hZgClosureNoRespDetAllFine", zgdet, detjet->Pt());
        fHistManager.FillTH2("hRgClosureNoRespDetAllFine", rgdet, detjet->Pt());
        fHistManager.FillTH2("hThetagClosureNoRespDetAllFine", thetagdet, detjet->Pt());
        fHistManager.FillTH2("hNsdClosureNoRespDetAllFine", nsddet, detjet->Pt());
      }
    }

    if(!hasMCSoftDrop)
      continue;

    // point needed at that stage for part level closure all
    // Distribtions before purity subtraction needed in order to test
    // the full correction chain (including purity correction and efficiency)
    // correction
    bool untaggedDet = softdropDet.fZg < fZcut,
         untaggedPart = softdropPart.fZg < fZcut;
    Double_t pointZg[4] = {softdropDet.fZg, detjet->Pt(), softdropPart.fZg, partjet->Pt()},
             pointRg[4] = {untaggedDet ? -0.01 : softdropDet.fRg, detjet->Pt(), untaggedPart ? -0.01 : softdropPart.fRg, partjet->Pt()},
             pointNsd[4] = {untaggedDet ? -1. : double(splittingsDet.size()), detjet->Pt(), untaggedPart ? -1. : double(splittingsPart.size()), partjet->Pt()},
             pointThetag[4] = {untaggedDet ? -0.05 : softdropDet.fRg/Rjet, detjet->Pt(), untaggedPart ? -0.05 : softdropPart.fRg/Rjet, partjet->Pt()};

    // Fill truth level histogram (for full closure test)
    if(fForceBeamType != kpp) {
      if(closureUseResponse) {
        fHistManager.FillTH2(Form("hZgClosureRespPartAllFine_%d", fCentBin), pointZg[kIndSDPart], pointZg[kIndPtPart]);
        fHistManager.FillTH2(Form("hRgClosureRespPartAllFine_%d", fCentBin), pointRg[kIndSDPart], pointRg[kIndPtPart]);
        fHistManager.FillTH2(Form("hThetagClosureRespPartAllFine_%d", fCentBin), pointThetag[kIndSDPart], pointThetag[kIndPtPart]);
        fHistManager.FillTH2(Form("hNsdClosureRespPartAllFine_%d", fCentBin), pointNsd[kIndSDPart], pointNsd[kIndPtPart]);
      } else {
        fHistManager.FillTH2(Form("hZgClosureNoRespPartAllFine_%d", fCentBin), pointZg[kIndSDPart], pointZg[kIndPtPart]);
        fHistManager.FillTH2(Form("hRgClosureNoRespPartAllFine_%d", fCentBin), pointRg[kIndSDPart], pointRg[kIndPtPart]);
        fHistManager.FillTH2(Form("hThetagClosureNoRespPartAllFine_%d", fCentBin), pointThetag[kIndSDPart], pointThetag[kIndPtPart]);
        fHistManager.FillTH2(Form("hNsdClosureNoRespPartAllFine_%d", fCentBin), pointNsd[kIndSDPart], pointNsd[kIndPtPart]);
      }
    } else {
      if(closureUseResponse) {
        fHistManager.FillTH2("hZgClosureRespPartAllFine", pointZg[kIndSDPart], pointZg[kIndPtPart]);
        fHistManager.FillTH2("hRgClosureRespPartAllFine", pointRg[kIndSDPart], pointRg[kIndPtPart]);
        fHistManager.FillTH2("hThetagClosureRespPartAllFine", pointThetag[kIndSDPart], pointThetag[kIndPtPart]);
        fHistManager.FillTH2("hNsdClosureRespPartAllFine", pointNsd[kIndSDPart], pointNsd[kIndPtPart]);
      } else {
        fHistManager.FillTH2("hZgClosureNoRespPartAllFine", pointZg[kIndSDPart], pointZg[kIndPtPart]);
        fHistManager.FillTH2("hRgClosureNoRespPartAllFine", pointRg[kIndSDPart], pointRg[kIndPtPart]);
        fHistManager.FillTH2("hThetagClosureNoRespPartAllFine", pointThetag[kIndSDPart], pointThetag[kIndPtPart]);
        fHistManager.FillTH2("hNsdClosureNoRespPartAllFine", pointNsd[kIndSDPart], pointNsd[kIndPtPart]);
      }
    }

    if(fRequirePartJetInAcceptance && ! partJetInAcceptance)
      continue;

    // Pair of jets for response found - fill QA
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

    if(fFillPlotsQAConstituents) {
      FillJetQA(*partjet, true, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy());
      FillJetQA(*detjet, false, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy());
    }

    auto deltaR = TMath::Abs(partjet->DeltaR(detjet));
    Double_t resZg = pointZg[kIndSDDet] - pointZg[kIndSDPart],
             resRg = pointRg[kIndSDDet] - pointRg[kIndSDPart],
             resThetag = pointThetag[kIndSDDet] - pointThetag[kIndSDPart],
             resNsd = pointNsd[kIndSDDet] - pointNsd[kIndSDPart];
    Double_t pointResRg[3] = {pointRg[kIndPtPart], pointRg[kIndSDPart], resRg},
             pointResThetag[3] = {pointThetag[kIndPtPart], pointThetag[kIndSDPart], resThetag},
             pointResRgNormalized[3] = {pointRg[kIndPtPart], pointRg[kIndSDPart],  resRg/pointRg[kIndSDPart]},
             pointResThetagNormalized[3] = {pointThetag[kIndPtPart], pointThetag[kIndSDPart], resThetag/pointThetag[kIndSDPart]};
    if (fForceBeamType != kpp) {
      if(fFillPlotsQAGeneral) {
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
        fHistManager.FillTHnSparse(Form("hResidualsRg_%d", fCentBin), pointResRg);
        fHistManager.FillTHnSparse(Form("hResidualsRgNormalized_%d", fCentBin), pointResRgNormalized);
        fHistManager.FillTHnSparse(Form("hResidualsThetag_%d", fCentBin), pointResThetag);
        fHistManager.FillTHnSparse(Form("hResidualsThetagNormalized_%d", fCentBin), pointResThetagNormalized);
      }
    } else {
      if(fFillPlotsQAGeneral) {
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
      }

      if(fFillPlotsQAOutliers) {
        // Monitor jet pt relative to pt-hard of the event
        if(fPtHard > 0.){
          fHistManager.FillTH1("hFracPtHardDet", detjet->Pt()/fPtHard);
          fHistManager.FillTH1("hFracPtHardPart", partjet->Pt() / fPtHard);
        } 
      }

      if(fHasResponseMatrixRooUnfold){
        fHistManager.FillTH1("hZgPartLevel", pointZg[kIndSDPart], pointZg[kIndPtPart]);
        fHistManager.FillTH1("hRgPartLevel", pointRg[kIndSDPart], pointRg[kIndPtPart]);
        fHistManager.FillTH1("hNsdPartLevel", pointNsd[kIndSDPart], pointNsd[kIndPtPart]);
        fHistManager.FillTH1("hThetagPartLevel", pointThetag[kIndSDPart], pointThetag[kIndPtPart]);
      }
      if(!untaggedDet && !untaggedPart) {
        // Fill residuals
        if(fFillPlotsResiduals) {
          fHistManager.FillTH2("hZgResiduals", pointZg[kIndPtPart], resZg);
          fHistManager.FillTH2("hRgResiduals", pointRg[kIndPtPart], resRg);
          fHistManager.FillTH2("hThetagResiduals",  pointThetag[kIndPtPart], resThetag);
          fHistManager.FillTH2("hNsdResiduals",  pointNsd[kIndPtPart], resNsd);
          fHistManager.FillTH2("hZgResidualsNormalized", pointZg[kIndPtPart], resZg/pointZg[kIndSDPart]);
          fHistManager.FillTH2("hRgResidualsNormalized", pointRg[kIndPtPart], resRg/pointRg[kIndSDPart]);
          fHistManager.FillTH2("hThetagResidualsNormalized", pointThetag[kIndPtPart], resThetag/pointThetag[kIndSDPart]);
          fHistManager.FillTH2("hNsdResidualsNormalized", pointNsd[kIndPtPart], resNsd/pointNsd[kIndSDPart]);
          fHistManager.FillTHnSparse("hResidualsRg", pointResRg);
          fHistManager.FillTHnSparse("hResidualsRgNormalized", pointResRgNormalized);
          fHistManager.FillTHnSparse("hResidualsThetag", pointResThetag);
          fHistManager.FillTHnSparse("hResidualsThetagNormalized", pointResThetagNormalized);
          if(pointRg[kIndSDPart] < 0.02) fHistManager.FillTH2("hResualsRgLowRg", pointRg[kIndPtPart], resRg/pointRg[kIndSDPart]);
        }
      }
    }
    if (detjet->Pt() >= ptmindet && detjet->Pt() <= ptmaxdet) {
      if (fForceBeamType != kpp) {
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
      } else {
        if(fHasResponseMatrixRooUnfold) {
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
      if (closureUseResponse) {
        if (fForceBeamType != kpp) {
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
        } else {
          if(fHasResponseMatrixRooUnfold) {
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
      } else {
        if (fForceBeamType != kpp) {
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
        } else {
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
  

  // second loop: Loop over all part. level jet and check whether it has a corresponding detctor level jet
  // this is of relevance for the jet finding efficiency
  if(fForceBeamType == kpp){
    for(auto partjet : partLevelJets->accepted()){
      SoftdropResults softdropPart, softdropDet;
      std::vector<SoftdropResults> splittingsPart, splittingsDet;
      try{
        softdropPart = MakeSoftdrop(*partjet, partLevelJets->GetJetRadius(), true, sdsettings, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy(), fVertex, fDropMass0Jets); 
        splittingsPart = IterativeDecluster(*partjet, partLevelJets->GetJetRadius(), true, sdsettings, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy(), fVertex, fDropMass0Jets);
      } catch (...) {
        continue;
      }
      bool untaggedEfficinency = softdropPart.fZg < fZcut;
      double zgpart = softdropPart.fZg,
             rgpart = untaggedEfficinency ? -0.01 : softdropPart.fRg,
             thetagpart = untaggedEfficinency ? -0.05 : softdropPart.fRg/Rjet,
             nsdpart = untaggedEfficinency ? -1. : double(splittingsPart.size());
      // Fill 2D part. level distributions
      fHistManager.FillTH2("hZgPartLevelFine", zgpart, partjet->Pt());
      fHistManager.FillTH2("hRgPartLevelFine", rgpart , partjet->Pt());
      fHistManager.FillTH2("hNsdPartLevelFine", nsdpart, partjet->Pt());
      fHistManager.FillTH2("hThetagPartLevelFine", thetagpart, partjet->Pt());

      // Fill histograms with fully efficient truth before matching (for closure test)
      if(closureUseResponse) {
        fHistManager.FillTH2("hZgClosureTruthRespPartAllFine", zgpart, partjet->Pt());
        fHistManager.FillTH2("hRgClosureTruthRespPartAllFine", rgpart, partjet->Pt());
        fHistManager.FillTH2("hThetagClosureTruthRespPartAllFine", thetagpart, partjet->Pt());
        fHistManager.FillTH2("hNsdClosureTruthRespPartAllFine", nsdpart, partjet->Pt());
      } else {
        fHistManager.FillTH2("hZgClosureTruthNoRespPartAllFine", zgpart, partjet->Pt());
        fHistManager.FillTH2("hRgClosureTruthNoRespPartAllFine", rgpart, partjet->Pt());
        fHistManager.FillTH2("hThetagClosureTruthNoRespPartAllFine", thetagpart, partjet->Pt());
        fHistManager.FillTH2("hNsdClosureTruthNoRespPartAllFine", nsdpart, partjet->Pt());
      }

      // Handle jet finding efficiency
      // Point efficiency
      // 0 - no matched jet
      // 1 - matched jet not having SoftDrop
      // 2 - matched jet not in the same acceptance
      // 3 - matched jet in the same acceptance and has SoftDrop 
      TagStatus_t tag = kNoMatchedJet;
      auto detjet = partjet->ClosestJet();
      if(detjet) {
        tag = kMatchedJetNoSoftDrop;
        // check if for the matched det. level jet we can determine the SoftDrop
        // check this condition first in case not the same acceptance type is required
        try {
          softdropDet = MakeSoftdrop(*detjet, detLevelJets->GetJetRadius(),false, sdsettings, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy(), fVertex, fDropMass0Jets);
          splittingsDet = IterativeDecluster(*detjet, detLevelJets->GetJetRadius(),false, sdsettings, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy(), fVertex, fDropMass0Jets);
          tag = kMatchedJetNoAcceptance;
          if(detjet->GetJetAcceptanceType() & partLevelJets->GetAcceptanceType()){
            tag = kPairAccepted;
          } else {
            AliDebugStream(1) << "det level jet not in part. level acceptance" << std::endl;
          }
        } catch (...) {
          // Nothing to do
          AliDebugStream(1) << "SoftDrop falied" << std::endl;
        }
      } else {
        AliDebugStream(1) << "No det level jet found for part level jet" << std::endl;
      }
      double pointEfficiencyZg[3] = {zgpart, partjet->Pt(), static_cast<double>(tag)},
               pointEfficiencyRg[3] = {rgpart, partjet->Pt(), static_cast<double>(tag)},
               pointEfficiencyThetag[3] = {thetagpart, partjet->Pt(), static_cast<double>(tag)},
               pointEfficiencyNsd[3] = {nsdpart, partjet->Pt(), static_cast<double>(tag)};
      fHistManager.FillTHnSparse("hZgJetFindingEfficiency", pointEfficiencyZg);
      fHistManager.FillTHnSparse("hRgJetFindingEfficiency", pointEfficiencyRg);
      fHistManager.FillTHnSparse("hThetagJetFindingEfficiency", pointEfficiencyThetag);
      fHistManager.FillTHnSparse("hNsdJetFindingEfficiency", pointEfficiencyNsd);
      // Split jet finding efficiency in response and no-response sample for closure test
      if(closureUseResponse) {
        fHistManager.FillTHnSparse("hZgJetFindingEfficiencyClosureResp", pointEfficiencyZg);
        fHistManager.FillTHnSparse("hRgJetFindingEfficiencyClosureResp", pointEfficiencyRg);
        fHistManager.FillTHnSparse("hThetagJetFindingEfficiencyClosureResp", pointEfficiencyThetag);
        fHistManager.FillTHnSparse("hNsdJetFindingEfficiencyClosureResp", pointEfficiencyNsd);
      } else {
        fHistManager.FillTHnSparse("hZgJetFindingEfficiencyClosureNoResp", pointEfficiencyZg);
        fHistManager.FillTHnSparse("hRgJetFindingEfficiencyClosureNoResp", pointEfficiencyRg);
        fHistManager.FillTHnSparse("hThetagJetFindingEfficiencyClosureNoResp", pointEfficiencyThetag);
        fHistManager.FillTHnSparse("hNsdJetFindingEfficiencyClosureNoResp", pointEfficiencyNsd);
      }
    }
  }
  // efficiency loop for PbPb to be implemented when of relevance

  return kTRUE;
}

void AliAnalysisTaskEmcalSoftDropResponse::FillJetQA(const AliEmcalJet &jet, bool isPartLevel, AliVCluster::VCluUserDefEnergy_t energydef) {
  std::string tag = isPartLevel ? "Part" : "Det";
  TVector3 maxcharged, maxneutral, jetvec(jet.Px(), jet.Py(), jet.Pz());
  bool hasMaxCharged = false, hasMaxNeutral = false;
  if (fUseChargedConstituents || isPartLevel)
  { // Neutral particles part of particle container in case of MC
    AliDebugStream(1) << "Jet substructure: Using charged constituents" << std::endl;
    for (int itrk = 0; itrk < jet.GetNumberOfTracks(); itrk++)
    {
      auto track = jet.Track(itrk);
      if (!track->Charge() && !fUseNeutralConstituents)
        continue; // Reject neutral constituents in case of using only charged consituents
      if (track->Charge() && !fUseChargedConstituents)
        continue; // Reject charged constituents in case of using only neutral consituents
      TVector3 trackvec(track->Px(), track->Py(), track->Pz());
      if(track->Charge()) {
        fHistManager.FillTH2(Form("hSDUsedChargedPtjvPtc%s", tag.data()), jet.Pt(), track->Pt());
        fHistManager.FillTH2(Form("hSDUsedChargedEtaPhi%s", tag.data()), track->Eta(), TVector2::Phi_0_2pi(track->Phi()));
        fHistManager.FillTH2(Form("hSDUsedChargedDR%s", tag.data()), jet.Pt(), jetvec.DeltaR(trackvec));
        if(!hasMaxCharged) {
          maxcharged = trackvec;
          hasMaxCharged = true;
        } else {
          if(trackvec.Pt() > maxcharged.Pt())
            maxcharged = trackvec;
        }
      } else {
        fHistManager.FillTH2("hSDUsedNeutralPtjvPtcPart", jet.Pt(), track->Pt());
        fHistManager.FillTH2("hSDUsedNeutralEtaPhiPart", track->Eta(), TVector2::Phi_0_2pi(track->Phi()));
        fHistManager.FillTH2("hSDUsedNeutralDRPart", jet.Pt(), jetvec.DeltaR(trackvec));
        if(!hasMaxNeutral) {
            maxneutral = trackvec;
            hasMaxNeutral = true;
        } else {
            if(trackvec.Pt() > maxneutral.Pt())
              maxneutral = trackvec;
        }
      }
    }
  }

  if (fUseNeutralConstituents)
  {
    AliDebugStream(1) << "Jet substructure: Using neutral constituents" << std::endl;
    for (int icl = 0; icl < jet.GetNumberOfClusters(); icl++)
    {
      auto cluster = jet.Cluster(icl);
      TLorentzVector clustervec;
      cluster->GetMomentum(clustervec, fVertex, energydef);
      TVector3 clustervec3(clustervec.Pt(), clustervec.Eta(), clustervec.Phi());
      fHistManager.FillTH2("hSDUsedNeutralPtjvPtcDet", jet.Pt(), clustervec.Pt());
      fHistManager.FillTH2("hSDUsedNeutralEtaPhiDet", clustervec.Eta(), TVector2::Phi_0_2pi(clustervec.Phi()));
      fHistManager.FillTH2("hSDUsedNeutralDRDet", jet.Pt(), jetvec.DeltaR(clustervec3));
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
      if(!hasMaxNeutral) {
        maxneutral = clustervec3;
        hasMaxNeutral = true;
      } else {
          if(clustervec3.Pt() > maxneutral.Pt())
            maxneutral = clustervec3;
      }
    }
  }

  if(hasMaxCharged) {
    fHistManager.FillTH2(Form("hSDUsedChargedPtjvPtcMax%s", tag.data()), jet.Pt(), maxcharged.Pt());
    fHistManager.FillTH2(Form("hSDUsedChargedEtaPhiMax%s", tag.data()), maxcharged.Eta(), TVector2::Phi_0_2pi(maxcharged.Phi()));
    fHistManager.FillTH2(Form("hSDUsedChargedDRMax%s", tag.data()), jet.Pt(), jetvec.DeltaR(maxcharged));
  }

  if(hasMaxNeutral) {
    fHistManager.FillTH2("hSDUsedNeutralPtjvPcMaxPart", jet.Pt(), maxneutral.Pt());
    fHistManager.FillTH2("hSDUsedNeutralEtaPhiMaxPart", maxneutral.Eta(), TVector2::Phi_0_2pi(maxneutral.Phi()));
    fHistManager.FillTH2("hSDUsedNeutralDRMaxPart", jet.Pt(), jetvec.DeltaR(maxneutral));
  }
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

void AliAnalysisTaskEmcalSoftDropResponse::ConfigurePtHard(MCProductionType_t mcprodtype, const TArrayI &pthardbinning, Bool_t doMCFilter, Double_t jetptcut) {
  SetMCProductionType(mcprodtype);
  SetUsePtHardBinScaling(true);
  SetUserPtHardBinning(pthardbinning);
  if(doMCFilter) {
    SetMCFilter();
    SetJetPtFactor(jetptcut);
  }
}

void AliAnalysisTaskEmcalSoftDropResponse::ConfigureMinBias(MCProductionType_t mcprodtype){
  if(!(mcprodtype == kMCPythiaMB || mcprodtype == kMCHepMCMB)) {
    AliErrorStream() << "MC prod type not compatible with min. bias production" << std::endl;
  }
  SetMCProductionType(mcprodtype);
}

void AliAnalysisTaskEmcalSoftDropResponse::ConfigureJetSelection(Double_t minJetPtPart, Double_t minJetPtDet, Double_t maxTrackPtPart, Double_t maxTrackPtDet, Double_t maxClusterPt, Double_t minAreaPerc) {
  auto partjets = GetPartLevelJetContainer(),
       detjets = GetDetLevelJetContainer();
  
  partjets->SetJetPtCut(minJetPtPart);
  partjets->SetMaxTrackPt(maxTrackPtPart);
  detjets->SetJetPtCut(minJetPtDet);
  if(detjets->GetJetType() == AliJetContainer::kFullJet || detjets->GetJetType() == AliJetContainer::kChargedJet) {
    detjets->SetMaxTrackPt(maxTrackPtDet);
  }
  if(detjets->GetJetType() == AliJetContainer::kFullJet || detjets->GetJetType() == AliJetContainer::kNeutralJet) {
    detjets->SetMaxClusterPt(maxClusterPt);
  }
  if(minAreaPerc >= 0.) {
    detjets->SetPercAreaCut(minAreaPerc);
  }
}

AliAnalysisTaskEmcalSoftDropResponse *AliAnalysisTaskEmcalSoftDropResponse::AddTaskEmcalSoftDropResponse(Double_t jetradius, AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recombinationScheme, AliVCluster::VCluUserDefEnergy_t energydef, bool ifembed, const char *namepartcont, const char *trigger)
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
    tracks = responsemaker->AddTrackContainer(AliEmcalAnalysisFactory::TrackContainerNameFactory(isAOD));
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
    clusters = responsemaker->AddClusterContainer(AliEmcalAnalysisFactory::ClusterContainerNameFactory(isAOD));
    std::cout << "Cluster container name: " << clusters->GetName() << std::endl;
    // 300 MeV E-cut
    switch(energydef) {
      case AliVCluster::kHadCorr: 
        clusters->SetClusHadCorrEnergyCut(0.3); 
        clusters->SetDefaultClusterEnergy(AliVCluster::kHadCorr);
        break;
      case AliVCluster::kNonLinCorr:
        clusters->SetClusNonLinCorrEnergyCut(0.3); 
        clusters->SetDefaultClusterEnergy(AliVCluster::kNonLinCorr);
        break;
      case AliVCluster::kUserDefEnergy1:
      case AliVCluster::kUserDefEnergy2:
      default:
        AliErrorGeneralStream("AliAnalysisTaskEmcalSoftDropResponse::AddTaskEmcalSoftDropResponse") << "Cluster energy type not supported" << std::endl;
    }
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
