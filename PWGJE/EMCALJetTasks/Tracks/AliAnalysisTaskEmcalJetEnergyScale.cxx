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
#include <algorithm>
#include <array>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <TFile.h>
#include <TH1.h>
#include <THistManager.h>
#include <TCustomBinning.h>
#include <TLinearBinning.h>
#include <TCustomBinning.h>
#include <TPDGCode.h>
#include <TRandom.h>

#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskEmcalJetEnergyScale.h"
#include "AliEmcalMCPartonInfo.h"
#include "AliEmcalTriggerDecisionContainer.h"
#include "AliEmcalAnalysisFactory.h"
#include "AliLog.h"
#include "AliVEventHandler.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetEnergyScale)
ClassImp(PWGJE::EMCALJetTasks::AngularityHandler)
ClassImp(PWGJE::EMCALJetTasks::AngularityHandler::AngularityBin)

using namespace PWGJE::EMCALJetTasks;

AliAnalysisTaskEmcalJetEnergyScale::AliAnalysisTaskEmcalJetEnergyScale():
  AliAnalysisTaskEmcalJet(),
  fHistos(nullptr),
  fAngularityHandler(nullptr),
  fNameTracks(),
  fNameClusters(),
  fNameMCParticles(),
  fNameDetectorJets(),
  fNameParticleJets(),
  fTriggerSelectionString(),
  fNameTriggerDecisionContainer("EmcalTriggerDecision"),
  fFractionResponseClosure(0.8),
  fFillHSparse(false),
  fScaleShift(0.),
  fRequireSameAcceptance(kFALSE),
  fUseStandardOutlierRejection(kFALSE),
  fDebugMaxJetOutliers(kFALSE),
  fJetTypeOutliers(kOutlierPartJet),
  fSampleSplitter(nullptr)
{
}

AliAnalysisTaskEmcalJetEnergyScale::AliAnalysisTaskEmcalJetEnergyScale(const char *name):
  AliAnalysisTaskEmcalJet(name, true),
  fHistos(nullptr),
  fAngularityHandler(nullptr),
  fNameTracks(),
  fNameClusters(),
  fNameMCParticles(),
  fNameDetectorJets(),
  fNameParticleJets(),
  fTriggerSelectionString(),
  fNameTriggerDecisionContainer("EmcalTriggerDecision"),
  fFractionResponseClosure(0.8),
  fFillHSparse(false),
  fScaleShift(0.),
  fRequireSameAcceptance(kTRUE),
  fUseStandardOutlierRejection(kFALSE),
  fDebugMaxJetOutliers(kFALSE),
  fJetTypeOutliers(kOutlierPartJet),
  fSampleSplitter(nullptr)
{
  SetMakeGeneralHistograms(true);
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskEmcalJetEnergyScale::~AliAnalysisTaskEmcalJetEnergyScale() {
  if(fHistos) delete fHistos;
  if(fSampleSplitter) delete fSampleSplitter;
}

void AliAnalysisTaskEmcalJetEnergyScale::UserCreateOutputObjects(){
  AliAnalysisTaskEmcal::UserCreateOutputObjects();

  TCustomBinning jetPtBinningCoarseDet, jetPtBinningCoarsePart;
  jetPtBinningCoarseDet.SetMinimum(20.);
  jetPtBinningCoarseDet.AddStep(40., 2.);
  jetPtBinningCoarseDet.AddStep(60., 5.);
  jetPtBinningCoarseDet.AddStep(120., 10.);
  jetPtBinningCoarseDet.AddStep(200., 20.);
  jetPtBinningCoarsePart.SetMinimum(0);
  jetPtBinningCoarsePart.AddStep(20., 20.);
  jetPtBinningCoarsePart.AddStep(80., 10.);
  jetPtBinningCoarsePart.AddStep(200., 20.);
  jetPtBinningCoarsePart.AddStep(280., 40.);
  jetPtBinningCoarsePart.AddStep(500., 220.);

  const double kPtDetMax = 500.,
               kPtPartMax = 800.;
  const int kNPtBinsDet = 500,
            kNPtBinsPart = 800;

  fHistos = new THistManager("energyScaleHistos");
  fHistos->CreateTH1("hEventCounter", "Event counter", 1, 0.5, 1.5);
  fHistos->CreateTH1("hHardestParton", "Pt of the hardest parton", 1000, 0., 1000.);
  fHistos->CreateTH1("hHardestQuark", "Pt of the hardest parton in case it is a quark", 1000, 0., 1000.);
  fHistos->CreateTH1("hHardestGluon", "Pt of the hardest parton in case it is a gluon", 1000, 0., 1000.);
  fHistos->CreateTH1("hAllPartons", "Pt of the hardest parton", 1000, 0., 1000.);
  fHistos->CreateTH1("hAllQuarks", "Pt of the hardest parton in case it is a quark", 1000, 0., 1000.);
  fHistos->CreateTH1("hAllGluons", "Pt of the hardest parton in case it is a gluon", 1000, 0., 1000.);
  fHistos->CreateTH2("hJetEnergyScale", "Jet Energy scale; p_{t,part} (GeV/c); (p_{t,det} - p_{t,part})/p_{t,part}" , 400, 0., 400., 200, -1., 1.);
  fHistos->CreateTH2("hJetEnergyScaleCharged", "Jet Energy scale (charged part); p_{t,part} (GeV/c); (p_{t,det} - p_{t,part})/p_{t,part}" , 400, 0., 400., 200, -1., 1.);
  fHistos->CreateTH2("hJetEnergyScaleNeutral", "Jet Energy scale (neutral part); p_{t,part} (GeV/c); (p_{t,det} - p_{t,part})/p_{t,part}" , 400, 0., 400., 200, -1., 1.);
  fHistos->CreateTH2("hJetEnergyScaleChargedVsFull", "Jet Energy scale (charged part, vs. full jet pt); p_{t,part} (GeV/c); (p_{t,det} - p_{t,part})/p_{t,part}" , 400, 0., 400., 200, -1., 1.);
  fHistos->CreateTH2("hJetEnergyScaleNeutralVsFull", "Jet Energy scale (neutral part, vs. full jet pt); p_{t,part} (GeV/c); (p_{t,det} - p_{t,part})/p_{t,part}" , 400, 0., 400., 200, -1., 1.);
  fHistos->CreateTH2("hJetEnergyScaleDet", "Jet Energy scale (det); p_{t,det} (GeV/c); (p_{t,det} - p_{t,part})/p_{t,part}" , 400, 0., 400., 200, -1., 1.);
  fHistos->CreateTH2("hJetEnergyScaleDetCharged", "Jet Energy scale (det, charged part); p_{t,det} (GeV/c); (p_{t,det} - p_{t,part})/p_{t,part}" , 400, 0., 400., 200, -1., 1.);
  fHistos->CreateTH2("hJetEnergyScaleDetNeutral", "Jet Energy scale (det, neutral part); p_{t,det} (GeV/c); (p_{t,det} - p_{t,part})/p_{t,part}" , 400, 0., 400., 200, -1., 1.);
  fHistos->CreateTH2("hJetEnergyScaleDetChargedVsFull", "Jet Energy scale (det, charged part, vs. full jet pt); p_{t,det} (GeV/c); (p_{t,det} - p_{t,part})/p_{t,part}" , 400, 0., 400., 200, -1., 1.);
  fHistos->CreateTH2("hJetEnergyScaleDetNeutralVsFull", "Jet Energy scale (det, neutral part, vs. full jet pt); p_{t,det} (GeV/c); (p_{t,det} - p_{t,part})/p_{t,part}" , 400, 0., 400., 200, -1., 1.);
  fHistos->CreateTH2("hJetResponseFine", "Response matrix, fine binning", kNPtBinsDet, 0., kPtDetMax, kNPtBinsPart, 0., kPtPartMax);
  fHistos->CreateTH2("hJetResponseFineCharged", "Response matrix, fine binning", kNPtBinsDet, 0., kPtDetMax, kNPtBinsPart, 0., kPtPartMax);
  fHistos->CreateTH2("hJetResponseFineNeutral", "Response matrix, fine binning", kNPtBinsDet, 0., kPtDetMax, kNPtBinsPart, 0., kPtPartMax);
  fHistos->CreateTH2("hJetResponseFineClosure", "Response matrix, fine binning, for closure test", kNPtBinsDet, 0., kPtDetMax, kNPtBinsPart, 0., kPtPartMax);
  fHistos->CreateTH2("hJetResponseFineNoClosure", "Response matrix, fine binning, for closure test", kNPtBinsDet, 0., kPtDetMax, kNPtBinsPart, 0., kPtPartMax);
  if(fAngularityHandler) {
    fHistos->CreateTH2("hJetResponseFineAngularityLow", "Response matrix, fine binning (low angularity)", kNPtBinsDet, 0., kPtDetMax, kNPtBinsPart, 0., kPtPartMax);
    fHistos->CreateTH2("hJetResponseFineAngularityHigh", "Response matrix, fine binning (high angularity)", kNPtBinsDet, 0., kPtDetMax, kNPtBinsPart, 0., kPtPartMax);
    fHistos->CreateTH2("hJetResponseFineAngularityLowClosure", "Response matrix, fine binning for closure test(low angularity)", kNPtBinsDet, 0., kPtDetMax, kNPtBinsPart, 0., kPtPartMax);
    fHistos->CreateTH2("hJetResponseFineAngularityLowNoClosure", "Response matrix, fine binning for closure test(low angularity)", kNPtBinsDet, 0., kPtDetMax, kNPtBinsPart, 0., kPtPartMax);
    fHistos->CreateTH2("hJetResponseFineAngularityHighClosure", "Response matrix, fine binning for closure test(high angularity)", kNPtBinsDet, 0., kPtDetMax, kNPtBinsPart, 0., kPtPartMax);
    fHistos->CreateTH2("hJetResponseFineAngularityHighNoClosure", "Response matrix, fine binning for closure test(high angularity)", kNPtBinsDet, 0., kPtDetMax, kNPtBinsPart, 0., kPtPartMax);
  }
  fHistos->CreateTH1("hJetSpectrumPartAll", "Part level jet pt spectrum ", kNPtBinsPart, 0., kPtPartMax);
  fHistos->CreateTH1("hJetSpectrumPartAllClosure", "Part level jet pt spectrum for closure test", kNPtBinsPart, 0., kPtPartMax);
  fHistos->CreateTH1("hJetSpectrumPartAllNoClosure", "Part level jet pt spectrum (no-closure sample)", kNPtBinsPart, 0., kPtPartMax);
  fHistos->CreateTH1("hJetSpectrumPartAllNoEvSel", "Part level jet pt spectrum (no det-level event selection)", kNPtBinsPart, 0., kPtPartMax);
  fHistos->CreateTH2("hPurityDet", "Det. level purity", kNPtBinsDet, 0., kPtDetMax, 3, -0.5, 2.5);
  fHistos->CreateTH2("hEtaPhiPartAllNoEvSel", "#eta-#phi distibution of pure part. level jets (no event selection)", 100, -1., 1., 100, 0, TMath::TwoPi());
  fHistos->CreateTH2("hPurityDetClosure", "Det. level purity for closure test", kNPtBinsDet, 0., kPtDetMax, 3, -0.5, 2.5);
  fHistos->CreateTH2("hPurityDetNoClosure", "Det. level purity (no-closure sample)", kNPtBinsDet, 0., kPtDetMax, 3, -0.5, 2.5);
  fHistos->CreateTH2("hJetfindingEfficiencyCore", "Part. level efficiency", kNPtBinsPart, 0., kPtPartMax, 3, -0.5, 2.5);
  fHistos->CreateTH2("hJetfindingEfficiencyCoreClosure", "Part. level efficiency for closure test", kNPtBinsPart, 0., kPtPartMax, 3, -0.5, 2.5);
  fHistos->CreateTH2("hJetfindingEfficiencyCoreNoClosure", "Part. level efficiency (no-closure sample)", kNPtBinsPart, 0., kPtPartMax, 3, -0.5, 2.5);
  fHistos->CreateTH2("hAngularityPart", "Angularity at partilce level", kNPtBinsPart, 0., kPtPartMax, 100, 0., 1.);
  fHistos->CreateTH2("hAngularityDet", "Angularity at detector level", kNPtBinsDet, 0., kPtDetMax, 100, 0., 1.);
  fHistos->CreateTH1("hAngularityErrorsPart", "Number of errors calculating the angularity (part level)", kNPtBinsPart, 0., kPtPartMax);
  fHistos->CreateTH1("hAngularityErrorsDet", "Number of errors calculating the angularity (det level)", kNPtBinsDet, 0., kPtDetMax);
  fHistos->CreateTH1("hAngularityErrorsMatrix", "Number of errors filling the response according to the angularity", kNPtBinsPart, 0., kPtPartMax);
  if(fAngularityHandler) {
    fHistos->CreateTH2("hAngularityPartLow", "Angularity at partilce level", kNPtBinsPart, 0., kPtPartMax, 100, 0., 1.);
    fHistos->CreateTH2("hAngularityPartHigh", "Angularity at partilce level", kNPtBinsPart, 0., kPtPartMax, 100, 0., 1.);
  }
  if(fFillHSparse){
    TLinearBinning jetPtBinningDet(kNPtBinsDet, 0., kPtDetMax), jetPtBinningPart(600, 0., 600), nefbinning(105, 0., 1.05), ptdiffbinning(200, -1., 1.), jetEtaBinning(100, -0.9, 0.9), jetPhiBinning(100, 0., TMath::TwoPi()),
                   subsampleBinning(2, -0.5, 1.5), deltaRbinning(20, 0., 1.), statusbinningEff(3, -0.5, 2.5);

    const TBinning *diffbinning[3] = {&jetPtBinningPart, &nefbinning, &ptdiffbinning},
                   *corrbinning[6] = {&jetPtBinningPart, &jetPtBinningDet, &nefbinning, &deltaRbinning,&subsampleBinning,&subsampleBinning},
                   *effbinning[3] = {&jetPtBinningPart, &jetPtBinningDet, &statusbinningEff};

    fHistos->CreateTHnSparse("hPtDiff", "pt diff det/part", 3, diffbinning, "s");
    fHistos->CreateTHnSparse("hPtCorr", "Correlation det pt / part pt", 6, corrbinning, "s");
    fHistos->CreateTHnSparse("hJetfindingEfficiency", "Jet finding efficiency", 3, effbinning, "s");
  }

  // Debugging the JES
  std::array<double, 17> ptbins = {0., 10., 20., 30., 40., 50., 60., 80., 100., 120., 140., 160., 180., 200., 240., 280., 320.};
  for(std::size_t iptbin = 0; iptbin < ptbins.size() -1; iptbin++){
    int ptbminI = ptbins[iptbin],
        ptbmaxI = ptbins[iptbin+1];
    fHistos->CreateTH2(Form("hJESVsNEFdet_%d_%d", ptbminI, ptbmaxI), Form("JES vs. NEF_{det} for jets with %d GeV/c < p_{t,part} < %d GeV/c", ptbminI, ptbmaxI), 105, 0.,  1.05, 200, -1., 1.);
    fHistos->CreateTH2(Form("hJESVsNEFpart_%d_%d", ptbminI, ptbmaxI), Form("JES vs. NEF_{part} for jets with %d GeV/c < p_{t,part} < %d GeV/c", ptbminI, ptbmaxI), 105, 0.,  1.05, 200, -1., 1.);
    fHistos->CreateTH2(Form("hJESVsNConstDet_%d_%d", ptbminI, ptbmaxI), Form("JES vs. N_{const,det} for jets with %d GeV/c < p_{t,part} < %d GeV/c", ptbminI, ptbmaxI), 101, -0.5,  100.5, 200, -1., 1.);
    fHistos->CreateTH2(Form("hJESVsNChargedDet_%d_%d", ptbminI, ptbmaxI), Form("JES vs. N_{charged,det} for jets with %d GeV/c < p_{t,part} < %d GeV/c", ptbminI, ptbmaxI), 101, -0.5,  100.5, 200, -1., 1.);
    fHistos->CreateTH2(Form("hJESVsNNeutralDet_%d_%d", ptbminI, ptbmaxI), Form("JES vs. N_{neutral,det} for jets with %d GeV/c < p_{t,part} < %d GeV/c", ptbminI, ptbmaxI), 101, -0.5,  100.5, 200, -1., 1.);
    fHistos->CreateTH2(Form("hJESVsNConstPart_%d_%d", ptbminI, ptbmaxI), Form("JES vs. N_{const,part} for jets with %d GeV/c < p_{t,part} < %d GeV/c", ptbminI, ptbmaxI), 101, -0.5,  100.5, 200, -1., 1.);
    fHistos->CreateTH2(Form("hJESVsNChargedPart_%d_%d", ptbminI, ptbmaxI), Form("JES vs. N_{charged,part} for jets with %d GeV/c < p_{t,part} < %d GeV/c", ptbminI, ptbmaxI), 101, -0.5,  100.5, 200, -1., 1.);
    fHistos->CreateTH2(Form("hJESVsNNeutralPart_%d_%d", ptbminI, ptbmaxI), Form("JES vs. N_{neutral,part} for jets with %d GeV/c < p_{t,part} < %d GeV/c", ptbminI, ptbmaxI), 101, -0.5,  100.5, 200, -1., 1.);
    fHistos->CreateTH2(Form("hCompNEF_%d_%d", ptbminI, ptbmaxI), Form("Comparisons NEF for part. and det. jets for jets with %d GeV/c < p_{t,part} < %d GeV/c", ptbminI, ptbmaxI), 100, 0., 1.05, 100, 0., 1.);
    fHistos->CreateTH2(Form("hCompNConst_%d_%d", ptbminI, ptbmaxI), Form("Comparisons Nconst for part. and det. jets for jets with %d GeV/c < p_{t,part} < %d GeV/c", ptbminI, ptbmaxI), 101, -0.5, 100.5, 101, -0.5, 100.5);
    fHistos->CreateTH2(Form("hCompNCharged_%d_%d", ptbminI, ptbmaxI), Form("Comparisons Nch for part. and det. jets for jets with %d GeV/c < p_{t,part} < %d GeV/c", ptbminI, ptbmaxI), 101, -0.5, 100.5, 101, -0.5, 100.5);
    fHistos->CreateTH2(Form("hCompNNeutral_%d_%d", ptbminI, ptbmaxI), Form("Comparisons Nne for part. and det. jets for jets with %d GeV/c < p_{t,part} < %d GeV/c", ptbminI, ptbmaxI), 101, -0.5, 100.5, 101, -0.5, 100.5);
  }

  // A bit of QA stuff
  fHistos->CreateTH2("hQANEFPtPart", "Neutral energy fraction at part. level; p_{t} (GeV/c); NEF", 350, 0., 350., 105, 0., 1.05);
  fHistos->CreateTH2("hQANEFPtDet", "Neutral energy fraction at det. level; p_{t} (GeV/c); NEF", 350, 0., 350., 105, 0., 1.05);
  fHistos->CreateTH2("hQAZchPtPart", "z_{ch,max} at part. level; p_{t} (GeV/c); z_{ch,max}", 350, 0., 350., 100, 0., 1.);
  fHistos->CreateTH2("hQAZchPtDet", "z_{ch,max} at det. level; p_{t} (GeV/c); z_{ch,max}", 350, 0., 350., 100, 0., 1.);
  fHistos->CreateTH2("hQAZnePtPart", "z_{ne,max} at part. level; p_{t} (GeV/c); z_{ne,max}", 350, 0., 350., 100, 0., 1.);
  fHistos->CreateTH2("hQAZnePtDet", "z_{ne,max} at det. level; p_{t} (GeV/c); z_{ne,max}", 350, 0., 350., 100, 0., 1.);
  fHistos->CreateTH2("hQANChPtPart", "Number of charged constituents at part. level; p_{t} (GeV/c); N_{ch}", 350, 0., 350., 100, 0., 100.);
  fHistos->CreateTH2("hQANChPtDet", "Number of charged constituents at det. level; p_{t} (GeV/c); N_{ch}", 350, 0., 350., 100, 0., 100.);
  fHistos->CreateTH2("hQANnePtPart", "Number of neutral constituents at part. level; p_{t} (GeV/c); N_{ne}", 350, 0., 350., 100, 0., 100.);
  fHistos->CreateTH2("hQANnePtDet", "Number of neutral constituents at det. level; p_{t} (GeV/c); N_{ne}", 350, 0., 350., 100, 0., 100.);
  fHistos->CreateTH2("hQAConstPtChPart", "p_{t} of charged constituents (part. level); p_{t,j} (GeV/c); p_{t,ch} (GeV/c", 350, 0., 350., 350., 0., 350.);
  fHistos->CreateTH2("hQAConstPtChDet", "p_{t} of charged constituents (det. level); p_{t,j} (GeV/c); p_{t,ch} (GeV/c", 350, 0., 350., 350., 0., 350.);
  fHistos->CreateTH2("hQAConstPtNePart", "p_{t} of neutral constituents (part. level); p_{t,j} (GeV/c); p_{t,ne} (GeV/c)", 350, 0., 350., 350., 0., 350.);
  fHistos->CreateTH2("hQAConstPtNeDet", "p_{t} of neutral constituents (det. level); p_{t,j} (GeV/c); p_{t,ne} (GeV/c)", 350, 0., 350., 350., 0., 350.);
  fHistos->CreateTH2("hQAConstPtChMaxPart", "p_{t} of max charged constituents (part. level); p_{t,j} (GeV/c); p_{t,ch} (GeV/c", 350, 0., 350., 350., 0., 350.);
  fHistos->CreateTH2("hQAConstPtChMaxDet", "p_{t} of max charged constituents (det. level); p_{t,j} (GeV/c); p_{t,ch} (GeV/c", 350, 0., 350., 350., 0., 350.);
  fHistos->CreateTH2("hQAConstPtNeMaxPart", "p_{t} of max neutral constituents (part. level); p_{t,j} (GeV/c); p_{t,ne} (GeV/c)", 350, 0., 350., 350., 0., 350.);
  fHistos->CreateTH2("hQAConstPtNeMaxDet", "p_{t} of max neutral constituents (det. level); p_{t,j} (GeV/c); p_{t,ne} (GeV/c)", 350, 0., 350., 350., 0., 350.);
  fHistos->CreateTH2("hQAEtaPhiPart", "#eta vs. #phi for selected part. level jets; #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hQAEtaPhiDet", "#eta vs. #phi for selected det. level jets; #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hQAEtaPhiConstChPart", "#eta vs. #phi for charged constituents (part. level); #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hQAEtaPhiConstChDet", "#eta vs. #phi for charged constituents (det. level); #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hQAEtaPhiConstNePart", "#eta vs. #phi for neutral constituents (part. level); #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hQAEtaPhiConstNeDet", "#eta vs. #phi for neutral constituents (det. level); #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hQAEtaPhiConstMaxChPart", "#eta vs. #phi for max charged constituents (part. level); #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hQAEtaPhiConstMaxChDet", "#eta vs. #phi for max charged constituents (det. level); #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hQAEtaPhiConstMaxNePart", "#eta vs. #phi for max neutral constituents (part. level); #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hQAEtaPhiConstMaxNeDet", "#eta vs. #phi for max neutral constituents (det. level); #eta; #phi", 100, -1., 1., 100, 0., 7.);
  fHistos->CreateTH2("hQADeltaRChargedPart", "#DeltaR vs. p_{t,jet} of charged constituents (part. level); p_{t, jet} (GeV/c); #DeltaR", 350., 0., 350, 100, 0., 1.);
  fHistos->CreateTH2("hQADeltaRChargedDet", "#DeltaR vs. p_{t,jet} of charged constituents (det. level); p_{t, jet} (GeV/c); #DeltaR", 350., 0., 350, 100, 0., 1.);
  fHistos->CreateTH2("hQADeltaRNeutralPart", "#DeltaR vs. p_{t,jet} of neutral constituents (part. level); p_{t, jet} (GeV/c); #DeltaR", 350., 0., 350, 100, 0., 1);
  fHistos->CreateTH2("hQADeltaRNeutralDet", "#DeltaR vs. p_{t,jet} of neutral constituents (det. level); p_{t, jet} (GeV/c); #DeltaR", 350., 0., 350, 100, 0., 1);
  fHistos->CreateTH2("hQADeltaRMaxChargedPart", "#DeltaR vs. p_{t,jet} of charged constituents (part. level); p_{t, jet} (GeV/c); #DeltaR", 350., 0., 350, 100, 0., 1.);
  fHistos->CreateTH2("hQADeltaRMaxChargedDet", "#DeltaR vs. p_{t,jet} of charged constituents (det. level); p_{t, jet} (GeV/c); #DeltaR", 350., 0., 350, 100, 0., 1.);
  fHistos->CreateTH2("hQADeltaRMaxNeutralPart", "#DeltaR vs. p_{t,jet} of neutral constituents (part. level); p_{t, jet} (GeV/c); #DeltaR", 350., 0., 350, 100, 0., 1);
  fHistos->CreateTH2("hQADeltaRMaxNeutralDet", "#DeltaR vs. p_{t,jet} of neutral constituents (det. level); p_{t, jet} (GeV/c); #DeltaR", 350., 0., 350, 100, 0., 1);
  fHistos->CreateTH2("hQAJetAreaVsJetPtPart", "Jet area vs. jet pt at particle level; p_{t} (GeV/c); Area", 350, 0., 350., 200, 0., 2.);
  fHistos->CreateTH2("hQAJetAreaVsJetPtDet", "Jet area vs. jet pt at detector level; p_{t} (GeV/c); Area", 350, 0., 350., 200, 0., 2.);
  fHistos->CreateTH2("hQAJetAreaVsNEFPart", "Jet area vs. NEF at particle level; NEF; Area", 105, 0., 1.05, 200, 0.,2.);
  fHistos->CreateTH2("hQAJetAreaVsNEFDet", "Jet area vs. NEF at detector level; NEF; Area", 105, 0., 1.05, 200, 0., 2.);
  fHistos->CreateTH2("hQAJetAreaVsNConstPart", "Jet area vs. number of consituents at particle level; Number of constituents; Area", 101, -0.5, 100.5, 200, 0., 2.);
  fHistos->CreateTH2("hQAJetAreaVsNConstDet", "Jet area vs. number of consituents at detector level; Number of constituents; Area", 101, -0.5, 100.5, 200, 0., 2.);
  fHistos->CreateTH1("hQAMatchingDRAbs", "Distance between part. level jet and  det. level jet", 100, 0., 1.);
  fHistos->CreateTH1("hQAMatchingDRel", "Distance between part. level jet and  det. level jet", 100, 0., 1.);
    // Cluster constituent QA
  fHistos->CreateTH2("hQAClusterTimeVsE", "Cluster time vs. energy; time (ns); E (GeV)", 1200, -600, 600, 200, 0., 200);
  fHistos->CreateTH2("hQAClusterTimeVsEFine", "Cluster time vs. energy (main region); time (ns); E (GeV)", 1000, -100, 100, 200, 0., 200);
  fHistos->CreateTH2("hQAClusterNCellVsE", "Cluster number of cells vs. energy; Number of cells; E (GeV)", 201, -0.5, 200.5, 200, 0., 200.);
  fHistos->CreateTH2("hQAClusterM02VsE", "Cluster M02 vs energy; M02; E (GeV)", 150, 0., 1.5, 200, 0., 200.);
  fHistos->CreateTH2("hQAClusterFracLeadingVsE", "Cluster frac leading cell vs energy; E (GeV); Frac. leading cell", 200, 0., 200., 110, 0., 1.1);
  fHistos->CreateTH2("hQAClusterFracLeadingVsNcell", "Cluster frac leading cell vs number of cells; Number of cells; Frac. leading cell", 201, -0.5, 200.5, 110, 0., 1.1);
  fHistos->CreateTH1("hPartConstPi0", "Particle-level consitutent spectrum of pi0 constituents", 200, 0., 200.);
  fHistos->CreateTH1("hPartConstK0s", "Particle-level consitutent spectrum of K0 constituents", 200, 0., 200.);
  fHistos->CreateTH1("hPartConstPhoton", "Particle-level consitutent spectrum of photon constituents", 200, 0., 200.);
  fHistos->CreateTH1("hFracPtHardPart", "Part. level jet Pt relative to the Pt-hard of the event", 100, 0., 10.);
  fHistos->CreateTH1("hFracPtHardDet", "Det. level jet Pt relative to the Pt-hard of the event", 100, 0., 10.);
  if(fDebugMaxJetOutliers) {
    fHistos->CreateTH2("hDebugMaxJetPt", "Debug Detection of max. part jet; p_{t,stl}; p_{t,man}", 350, 0., 350, 350., 0., 350.);
    fHistos->CreateTH1("hDebugMaxJetError", "Number of cases where the max. jet pt differs between methods", 1., 0.5, 1.5);
  }
  for(auto h : *(fHistos->GetListOfHistograms())) fOutput->Add(h);

  fSampleSplitter = new TRandom(0);

  if(fAngularityHandler) {
    AliDebugStream(1) << *fAngularityHandler << std::endl;
  }

  PostData(1, fOutput);
}

Bool_t AliAnalysisTaskEmcalJetEnergyScale::CheckMCOutliers() {
  if(!fMCRejectFilter) return true;
  if(!(fIsPythia || fIsHerwig || fIsHepMC)) return true;    // Only relevant for pt-hard production
  if(fUseStandardOutlierRejection) return AliAnalysisTaskEmcal::CheckMCOutliers();
  AliDebugStream(2) << "Using custom MC outlier rejection" << std::endl;
  AliJetContainer *outlierjets(nullptr);
  switch(fJetTypeOutliers) {
    case kOutlierPartJet: outlierjets = GetJetContainer(fNameParticleJets); break;
    case kOutlierDetJet: outlierjets = GetJetContainer(fNameDetectorJets); break;
  };
  if(!outlierjets) return true;

  // Check whether there is at least one particle level jet with pt above n * event pt-hard
  auto jetiter = outlierjets->accepted();
  auto max = std::max_element(jetiter.begin(), jetiter.end(), [](const AliEmcalJet *lhs, const AliEmcalJet *rhs ) { return lhs->Pt() < rhs->Pt(); });
  if(max != jetiter.end())  {
    // At least one jet found with pt > n * pt-hard
    AliDebugStream(2) << "Found max jet with pt " << (*max)->Pt() << " GeV/c" << std::endl;
    if(fDebugMaxJetOutliers) {
      // cross check whether implemenation using stl gives the same result as a trivial manual iteration
      // over all jets
      // for debug purposes
      auto jetiter1 = outlierjets->accepted();
      AliEmcalJet *debugmax(nullptr);
      for(auto testjet : jetiter1) {
        if(!debugmax) debugmax = testjet;
        else {
          if(testjet->Pt() > debugmax->Pt()) debugmax = testjet;
        }
      }
      fHistos->FillTH2("hDebugMaxJetPt", (*max)->Pt(), debugmax->Pt());
      if(*max != debugmax) fHistos->FillTH1("hDebugMaxJetError", 1.);
    }
    if((*max)->Pt() > fPtHardAndJetPtFactor * fPtHard) return false;
  }
  return true;
}

void AliAnalysisTaskEmcalJetEnergyScale::UserRunBeforeEventSelection(){
  auto partjets = GetPartLevelJetContainer();
   // efficiency x acceptance: Add histos for all accepted and reconstucted accepted jets
  for(auto partjet : partjets->accepted()){
    Double_t partjetpt = partjet->Pt();
    if (fDoBkgSub && partjets->GetRhoParameter()){
      partjetpt = partjetpt - partjets->GetRhoVal() * partjet->Area();
    }
    fHistos->FillTH1("hJetSpectrumPartAllNoEvSel", partjetpt);
    fHistos->FillTH1("hEtaPhiPartAllNoEvSel", partjet->Eta(), TVector2::Phi_0_2pi(partjet->Phi()));
  }
}

Bool_t AliAnalysisTaskEmcalJetEnergyScale::Run(){
  AliDebugStream(2) << "Next event" << std::endl;
  if(!(fInputHandler->IsEventSelected() & AliVEvent::kINT7)) return false;
  if(IsSelectEmcalTriggers(fTriggerSelectionString.Data())){
    auto mctrigger = static_cast<PWG::EMCAL::AliEmcalTriggerDecisionContainer *>(fInputEvent->FindListObject(fNameTriggerDecisionContainer));
    AliDebugStream(3) << "Found trigger decision object: " << (mctrigger ? "yes" : "no") << std::endl;
    if(!mctrigger){
      AliErrorStream() <<  "Trigger decision container with name " << fNameTriggerDecisionContainer << " not found in event - not possible to select EMCAL triggers" << std::endl;
      return false;
    }
    if(!mctrigger->IsEventSelected(fTriggerSelectionString)) return false;
  }
  AliDebugStream(3) << "event selected" << std::endl;
  fHistos->FillTH1("hEventCounter", 1);

  if(fMCPartonInfo) {
    auto hardest = fMCPartonInfo->GetHardestParton();
    if(hardest) {
      fHistos->FillTH1("hHardestParton", hardest->GetMomentum().Pt());
      if(hardest->GetPdg() == 21) fHistos->FillTH1("hHardestGluon", hardest->GetMomentum().Pt());
      if(hardest->GetPdg() < 7) fHistos->FillTH1("hHardestQuark", hardest->GetMomentum().Pt());
    }
    for(auto partonObj : fMCPartonInfo->GetListOfDirectPartons()) {
      auto parton = static_cast<PWG::EMCAL::AliEmcalPartonData *>(partonObj);
      if(parton) {
        fHistos->FillTH1("hAllPartons", parton->GetMomentum().Pt());
        if(parton->GetPdg() == 21) fHistos->FillTH1("hAllGluons", parton->GetMomentum().Pt());
        if(parton->GetPdg() < 7) fHistos->FillTH1("hAllQuarks", parton->GetMomentum().Pt());
      }
    }
  }

  // SampleSplitting for closure test on event level
  bool isClosure = fSampleSplitter->Uniform() < fFractionResponseClosure;

  auto detjets = GetDetLevelJetContainer(),
       partjets = GetPartLevelJetContainer();
  AliDebugStream(4) << "Detjets has array name : " << detjets->GetArrayName() << std::endl;
  AliDebugStream(4) << "Partjets has array name : " << partjets->GetArrayName() << std::endl;
  if(!detjets || !partjets) {
    AliErrorStream() << "At least one jet container missing, exiting ..." << std::endl;
    return false;
  }
  AliClusterContainer *clusters(detjets->GetClusterContainer());
  AliTrackContainer *tracks(static_cast<AliTrackContainer *>(detjets->GetParticleContainer()));
  AliParticleContainer *particles(partjets->GetParticleContainer());
  AliDebugStream(3) << "Have both jet containers: part(" << partjets->GetNAcceptedJets() << "|" << partjets->GetNJets() << "), det(" << detjets->GetNAcceptedJets() << "|" << detjets->GetNJets() << ")" << std::endl;

  std::array<double, 17> ptbinsDebug = {0., 10., 20., 30., 40., 50., 60., 80., 100., 120., 140., 160., 180., 200., 240., 280., 320.};
  std::vector<AliEmcalJet *> acceptedjets;
  for(auto detjet : detjets->accepted()){
    AliDebugStream(4) << "Next jet" << std::endl;
    auto partjet = detjet->ClosestJet();
    if(!partjet) {
      AliDebugStream(4) << "No tagged jet" << std::endl;
      fHistos->FillTH2("hPurityDet", detjet->Pt(), 0);
      if(isClosure) {
        fHistos->FillTH2("hPurityDetClosure", detjet->Pt(), 0);
      } else {
        fHistos->FillTH2("hPurityDetNoClosure", detjet->Pt(), 0);
      }
      continue;
    } else {
      if(!(partjet->GetJetAcceptanceType() & detjets->GetAcceptanceType())){
        // Acceptance not matching
        fHistos->FillTH2("hPurityDet", detjet->Pt(), 1);
        if(isClosure){
          fHistos->FillTH2("hPurityDetClosure", detjet->Pt(), 1);
        } else {
          fHistos->FillTH2("hPurityDetNoClosure", detjet->Pt(), 1);
        }
        if(fRequireSameAcceptance) continue;
      } else {
        fHistos->FillTH2("hPurityDet", detjet->Pt(), 2);
        if(isClosure) {
          fHistos->FillTH2("hPurityDetClosure", detjet->Pt(), 2);
        } else {
          fHistos->FillTH2("hPurityDetNoClosure", detjet->Pt(), 2);
        }
      }
    }
    double angularityPart = -1., angularityDet = -1.;
    try {
      angularityPart = MakeAngularity(*partjet, AliVCluster::VCluUserDefEnergy_t::kUserDefEnergy1),
      fHistos->FillTH2("hAngularityPart", partjet->Pt(), angularityPart);
    } catch(AngularityException &e) {
      fHistos->FillTH1("hAngularityErrorsPart", partjet->Pt());
    }
    try {
      AliVCluster::VCluUserDefEnergy_t energydef = AliVCluster::VCluUserDefEnergy_t::kUserDefEnergy1;
      auto detjetscusters = detjets->GetClusterContainer();
      if(detjetscusters){
       energydef = (AliVCluster::VCluUserDefEnergy_t)detjetscusters->GetDefaultClusterEnergy();
      }
      angularityDet = MakeAngularity(*detjet, energydef);
      fHistos->FillTH2("hAngularityDet", detjet->Pt(), angularityDet);
    } catch(AngularityException &e) {
      fHistos->FillTH1("hAngularityErrorsDet", partjet->Pt());
    }
    acceptedjets.push_back(detjet);
    Double_t detjetpt = detjet->Pt();
    Double_t partjetpt = partjet->Pt();
    if(TMath::Abs(fScaleShift) > DBL_EPSILON){
      detjetpt += fScaleShift * detjetpt;
    }

    if (fDoBkgSub && detjets->GetRhoParameter()){
      detjetpt = detjetpt - detjets->GetRhoVal() * detjet->Area();
      partjetpt = partjetpt - partjets->GetRhoVal() * partjet->Area();
    }

    if(fFillHSparse) {
      Bool_t acceptancematch = false;
      if (partjet->GetJetAcceptanceType() & detjets->GetAcceptanceType()) acceptancematch = true;
      TVector3 basevec, tagvec;
      basevec.SetPtEtaPhi(detjetpt, detjet->Eta(), detjet->Phi());
      tagvec.SetPtEtaPhi(partjetpt, partjet->Eta(), partjet->Phi());
      double pointCorr[6] = {partjetpt, detjetpt, detjet->NEF(), basevec.DeltaR(tagvec), acceptancematch ? 1. : 0.,  isClosure ? 0. : 1.},
             pointDiff[3] = {partjetpt, detjet->NEF(), (detjetpt-partjetpt)/partjetpt};
      fHistos->FillTHnSparse("hPtDiff", pointDiff);
      fHistos->FillTHnSparse("hPtCorr", pointCorr);
    }

    double detptcharged = detjetpt * (1 - detjet->NEF()), detptneutral = detjetpt * detjet->NEF(),
           partptcharged = partjetpt * (1 - partjet->NEF()),  partptneutral = partjetpt * partjet->NEF();
    fHistos->FillTH2("hJetResponseFine", detjetpt, partjetpt);
    if(partptcharged) fHistos->FillTH1("hJetResponseFineCharged", detptcharged, partptcharged);
    if(partptneutral) fHistos->FillTH1("hJetResponseFineNeutral", detptneutral, partptneutral);
    fHistos->FillTH1("hJetEnergyScale", partjetpt, (detjetpt - partjetpt)/partjetpt);
    if(partptcharged) fHistos->FillTH1("hJetEnergyScaleCharged", partptcharged, (detptcharged - partptcharged)/partptcharged);
    if(partptneutral) fHistos->FillTH1("hJetEnergyScaleNeutral", partptneutral, (detptneutral - partptneutral)/partptneutral);
    if(partptcharged) fHistos->FillTH1("hJetEnergyScaleChargedVsFull", partjetpt, (detptcharged - partptcharged)/partptcharged);
    if(partptneutral) fHistos->FillTH1("hJetEnergyScaleNeutralVsFull", partjetpt, (detptneutral - partptneutral)/partptneutral);
    fHistos->FillTH1("hJetEnergyScaleDet", detjetpt, (detjetpt - partjetpt)/partjetpt);
    if(partptcharged) fHistos->FillTH1("hJetEnergyScaleDetCharged", detptcharged, (detptcharged - partptcharged)/partptcharged);
    if(partptneutral) fHistos->FillTH1("hJetEnergyScaleDetNeutral", detptneutral, (detptneutral - partptneutral)/partptneutral);
    if(partptcharged) fHistos->FillTH1("hJetEnergyScaleDetChargedVsFull", detjetpt, (detptcharged - partptcharged)/partptcharged);
    if(partptneutral) fHistos->FillTH1("hJetEnergyScaleDetNeutralVsFull", detjetpt, (detptneutral - partptneutral)/partptneutral);
    // splitting for closure test
    if(isClosure) {
      fHistos->FillTH2("hJetResponseFineClosure", detjetpt, partjetpt);
    } else {
      fHistos->FillTH2("hJetResponseFineNoClosure", detjetpt, partjetpt);
    }
    // splitting response in low / high angularity
    if(fAngularityHandler && angularityPart > -1.) {
      try{
        bool isAngularityHigh = fAngularityHandler->IsHigherAngularity(partjetpt, angularityPart);
        if(isAngularityHigh) {
          fHistos->FillTH2("hAngularityPartHigh", partjetpt, angularityPart);
          fHistos->FillTH2("hJetResponseFineAngularityHigh", detjetpt, partjetpt);
          if(isClosure) {
            fHistos->FillTH2("hJetResponseFineAngularityHighClosure", detjetpt, partjetpt);
          } else {
            fHistos->FillTH2("hJetResponseFineAngularityHighNoClosure", detjetpt, partjetpt);
          }
        } else {
          fHistos->FillTH2("hAngularityPartLow", partjetpt, angularityPart);
          fHistos->FillTH2("hJetResponseFineAngularityLow", detjetpt, partjetpt);
          if(isClosure) {
            fHistos->FillTH2("hJetResponseFineAngularityLowClosure", detjetpt, partjetpt);
          } else {
            fHistos->FillTH2("hJetResponseFineAngularityLowNoClosure", detjetpt, partjetpt);
          }
        }
      } catch(AngularityHandler::BinNotFoundException &e) {
         fHistos->FillTH1("hAngularityErrorsMatrix", partjetpt);
      }
    }

    // Fill histograms for JES debugging
    int ptbminI = -1,
        ptbmaxI = -1;
    for(std::size_t iptbin = 0; iptbin < ptbinsDebug.size() - 1; iptbin++){
      if(partjetpt >= ptbinsDebug[iptbin] && partjetpt < ptbinsDebug[iptbin+1]) {
        ptbminI = ptbinsDebug[iptbin];
        ptbmaxI = ptbinsDebug[iptbin+1];
        break;
      }
    }
    if(ptbminI > -1 && ptbmaxI > -1){
      double jes = (detjetpt - partjetpt)/partjetpt;
      fHistos->FillTH2(Form("hJESVsNEFdet_%d_%d", ptbminI, ptbmaxI), detjet->NEF(), jes);
      fHistos->FillTH2(Form("hJESVsNEFpart_%d_%d", ptbminI, ptbmaxI), partjet->NEF(), jes);
      fHistos->FillTH2(Form("hJESVsNConstDet_%d_%d", ptbminI, ptbmaxI), detjet->N(), jes);
      fHistos->FillTH2(Form("hJESVsNChargedDet_%d_%d", ptbminI, ptbmaxI), detjet->Nch(), jes);
      fHistos->FillTH2(Form("hJESVsNNeutralDet_%d_%d", ptbminI, ptbmaxI), detjet->Nn(), jes);
      fHistos->FillTH2(Form("hJESVsNConstPart_%d_%d", ptbminI, ptbmaxI), partjet->N(), jes);
      fHistos->FillTH2(Form("hJESVsNChargedPart_%d_%d", ptbminI, ptbmaxI), partjet->Nch(), jes);
      fHistos->FillTH2(Form("hJESVsNNeutralPart_%d_%d", ptbminI, ptbmaxI), partjet->Nn(), jes);
      // Add plots correlating the NEF and number of constituents between part. and det. level jets
      fHistos->FillTH2(Form("hCompNEF_%d_%d", ptbminI, ptbmaxI), partjet->NEF(), detjet->NEF());
      fHistos->FillTH2(Form("hCompNConst_%d_%d", ptbminI, ptbmaxI), partjet->N(), detjet->N());
      fHistos->FillTH2(Form("hCompNCharged_%d_%d", ptbminI, ptbmaxI), partjet->Nch(), detjet->Nch());
      fHistos->FillTH2(Form("hCompNNeutral_%d_%d", ptbminI, ptbmaxI), partjet->Nn(), detjet->Nn());
    }

    // Fill QA histograms
    fHistos->FillTH2("hQANEFPtPart", partjetpt, partjet->NEF());
    fHistos->FillTH2("hQANEFPtDet", detjetpt, detjet->NEF());
    fHistos->FillTH2("hQAEtaPhiPart", partjet->Eta(), TVector2::Phi_0_2pi(partjet->Phi()));
    fHistos->FillTH2("hQAEtaPhiDet", detjet->Eta(), TVector2::Phi_0_2pi(detjet->Phi()));
    auto deltaR = TMath::Abs(partjet->DeltaR(detjet));
    fHistos->FillTH1("hQAMatchingDRAbs", deltaR);
    fHistos->FillTH1("hQAMatchingDRel", deltaR/partjets->GetJetRadius());
    TVector3 jetvecDet(detjet->Px(), detjet->Py(), detjet->Px());
    if(clusters){
      auto leadcluster = detjet->GetLeadingCluster(clusters->GetArray());
      fHistos->FillTH2("hQANnePtDet", detjetpt, detjet->GetNumberOfClusters());

      for(auto iclust = 0; iclust < detjet->GetNumberOfClusters(); iclust++) {
        auto cluster = detjet->Cluster(iclust);
        TLorentzVector clustervec;
        cluster->GetMomentum(clustervec, fVertex, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy());
        fHistos->FillTH2("hQAConstPtNeDet", detjetpt, clustervec.Pt());
        fHistos->FillTH2("hQAEtaPhiConstNeDet", clustervec.Eta(), TVector2::Phi_0_2pi(clustervec.Phi()));
        fHistos->FillTH2("hQADeltaRNeutralDet", detjetpt, jetvecDet.DeltaR(clustervec.Vect()));
        fHistos->FillTH2("hQAClusterTimeVsE", cluster->GetTOF() * 1e9 - 600, clustervec.E());         // time in ns., apply 600 ns time shift
        fHistos->FillTH2("hQAClusterTimeVsEFine", cluster->GetTOF() * 1e9 - 600, clustervec.E());     // time in ns., apply 600 ns time shift
        fHistos->FillTH2("hQAClusterNCellVsE", cluster->GetNCells(), clustervec.E());
        fHistos->FillTH2("hQAClusterM02VsE", cluster->GetM02(), clustervec.E());
        double maxamplitude = 0.;
        for(int icell = 0; icell < cluster->GetNCells(); icell++) {
          double amplitude = fInputEvent->GetEMCALCells()->GetAmplitude(fInputEvent->GetEMCALCells()->GetCellPosition(cluster->GetCellAbsId(icell)));
          if(amplitude > maxamplitude) maxamplitude = amplitude;
        }
        fHistos->FillTH2("hQAClusterFracLeadingVsE", clustervec.E(), maxamplitude/cluster->E());
        fHistos->FillTH2("hQAClusterFracLeadingVsNcell", cluster->GetNCells(), maxamplitude/cluster->E());
      }

      if(leadcluster){
        TLorentzVector ptvec;
        leadcluster->GetMomentum(ptvec, fVertex, (AliVCluster::VCluUserDefEnergy_t)clusters->GetDefaultClusterEnergy());
        fHistos->FillTH2("hQAZnePtDet", detjetpt, detjet->GetZ(ptvec.Px(), ptvec.Py(), ptvec.Pz()));
        fHistos->FillTH2("hQAConstPtNeMaxDet", detjetpt, ptvec.Pt());
        fHistos->FillTH2("hQAEtaPhiConstMaxNeDet", ptvec.Eta(), TVector2::Phi_0_2pi(ptvec.Phi()));
        fHistos->FillTH2("hQADeltaRMaxNeutralDet", detjetpt, jetvecDet.DeltaR(ptvec.Vect()));
      }
    }
    if(tracks){
      fHistos->FillTH2("hQANChPtDet", detjetpt,  detjet->GetNumberOfTracks());
      auto leadingtrack = detjet->GetLeadingTrack(tracks->GetArray());

      for(int itrk = 0; itrk < detjet->GetNumberOfTracks(); itrk++) {
        auto trk = detjet->Track(itrk);
        fHistos->FillTH2("hQAConstPtChDet", detjetpt, trk->Pt());
        fHistos->FillTH2("hQAEtaPhiConstChDet", trk->Eta(), TVector2::Phi_0_2pi(trk->Phi()));
        fHistos->FillTH2("hQADeltaRChargedDet", detjetpt, detjet->DeltaR(trk));
      }

      if(leadingtrack){
        fHistos->FillTH2("hQAZchPtDet", detjetpt, detjet->GetZ(leadingtrack->Px(), leadingtrack->Py(), leadingtrack->Pz()));
        fHistos->FillTH2("hQAConstPtChMaxDet", detjetpt, leadingtrack->Pt());
        fHistos->FillTH2("hQAEtaPhiConstMaxChDet", leadingtrack->Eta(), leadingtrack->Phi());
        fHistos->FillTH2("hQADeltaRMaxChargedDet", detjetpt, detjet->DeltaR(leadingtrack));
      }
    }
    if(particles){
      AliVParticle *leadingcharged(nullptr), *leadingneutral(nullptr);
      int ncharged(0), nneutral(0);
      for(int ipart = 0; ipart < partjet->GetNumberOfTracks(); ipart++) {
        auto particle = partjet->Track(ipart);
        if(particle->Charge()) {
          ncharged++;
          fHistos->FillTH2("hQAConstPtChPart", partjetpt, particle->Pt());
          fHistos->FillTH2("hQAEtaPhiConstChPart", particle->Eta(), TVector2::Phi_0_2pi(particle->Phi()));
          fHistos->FillTH2("hQADeltaRChargedPart", partjetpt, partjet->DeltaR(particle));
          if(!leadingcharged) leadingcharged = particle;
          else {
            if(particle->E() > leadingcharged->E()) leadingcharged = particle;
          }
        } else {
          nneutral++;
          fHistos->FillTH2("hQAConstPtNePart", partjetpt, particle->Pt());
          fHistos->FillTH2("hQAEtaPhiConstNePart", particle->Eta(), TVector2::Phi_0_2pi(particle->Phi()));
          fHistos->FillTH2("hQADeltaRNeutralPart", partjetpt, partjet->DeltaR(particle));
          if(!leadingneutral) leadingneutral = particle;
          else {
            if(particle->E() > leadingneutral->E()) leadingneutral = particle;
          }
        }
      }
      if(leadingcharged) {
        fHistos->FillTH2("hQAConstPtChMaxPart", partjetpt, leadingcharged->Pt());
        fHistos->FillTH2("hQAEtaPhiConstMaxChPart", leadingcharged->Eta(), TVector2::Phi_0_2pi(leadingcharged->Phi()));
        fHistos->FillTH2("hQAZchPtPart", partjetpt, partjet->GetZ(leadingcharged));
        fHistos->FillTH2("hQADeltaRMaxChargedPart", partjetpt, partjet->DeltaR(leadingcharged));
      }
      if(leadingneutral) {
        fHistos->FillTH2("hQAConstPtNeMaxPart", partjetpt, leadingneutral->Pt());
        fHistos->FillTH2("hQAEtaPhiConstMaxNePart", leadingneutral->Eta(), TVector2::Phi_0_2pi(leadingneutral->Phi()));
        fHistos->FillTH2("hQAZnePtPart", partjetpt, partjet->GetZ(leadingneutral));
        fHistos->FillTH2("hQADeltaRMaxNeutralPart", partjetpt, partjet->DeltaR(leadingneutral));
      }
      fHistos->FillTH2("hQANChPtPart", partjetpt, ncharged);
      fHistos->FillTH2("hQANnePtPart", partjetpt, nneutral);
    }
    fHistos->FillTH2("hQAJetAreaVsJetPtPart", partjetpt, partjet->Area());
    fHistos->FillTH2("hQAJetAreaVsJetPtDet", detjetpt, detjet->Area());
    fHistos->FillTH2("hQAJetAreaVsNEFPart", partjet->NEF(), partjet->Area());
    fHistos->FillTH2("hQAJetAreaVsNEFDet", detjet->NEF(), detjet->Area());
    fHistos->FillTH2("hQAJetAreaVsNConstPart", partjet->GetNumberOfTracks(), partjet->Area());
    fHistos->FillTH2("hQAJetAreaVsNConstDet", detjet->GetNumberOfClusters() + detjet->GetNumberOfTracks(), detjet->Area());
    fHistos->FillTH1("hFracPtHardPart", partjetpt/fPtHard);
    fHistos->FillTH1("hFracPtHardDet", detjetpt/fPtHard);
  }

  // efficiency x acceptance: Add histos for all accepted and reconstucted accepted jets
  for(auto partjet : partjets->accepted()){
    Double_t partjetpt = partjet->Pt();
    if (fDoBkgSub && partjets->GetRhoParameter()){
      partjetpt = partjetpt - partjets->GetRhoVal() * partjet->Area();
    }

    fHistos->FillTH1("hJetSpectrumPartAll", partjetpt);
    if(isClosure) {
      fHistos->FillTH1("hJetSpectrumPartAllClosure", partjetpt);
    } else {
      fHistos->FillTH1("hJetSpectrumPartAllNoClosure", partjetpt);
    }
    for(auto itrk = 0; itrk < partjet->GetNumberOfTracks(); itrk++) {
      auto constituent = partjet->Track(itrk);
      switch(TMath::Abs(constituent->PdgCode())) {
      case kPi0: fHistos->FillTH1("hPartConstPi0", constituent->Pt()); break;
      case kK0Short: fHistos->FillTH1("hPartConstK0s", constituent->Pt()); break;
      case kGamma: fHistos->FillTH1("hPartConstPhoton", constituent->Pt()); break;
      default: break;
      };
    }
    auto detjet = partjet->ClosestJet();
    int tagstatus = 0;
    if(detjet) {
      // check whether the matched det. level jet is in the part. level acceptance
      if(detjet->GetJetAcceptanceType() & partjets->GetAcceptanceType()) tagstatus = 2;
      else tagstatus = 1;
    }
    fHistos->FillTH2("hJetfindingEfficiencyCore", partjetpt, tagstatus);
    if(isClosure) {
      fHistos->FillTH2("hJetfindingEfficiencyCoreClosure", partjetpt, tagstatus);
    } else {
      fHistos->FillTH2("hJetfindingEfficiencyCoreNoClosure", partjetpt, tagstatus);
    }
    if(fFillHSparse){
      double effvec[3] = {partjetpt, 0.,static_cast<double>(tagstatus)};
      if(detjet) {
        // Found a match
        Double_t detjetpt = detjet->Pt();
        if (fDoBkgSub && detjets->GetRhoParameter()){
          detjetpt = detjetpt - detjets->GetRhoVal() * detjet->Area();
        }

        effvec[1] = detjetpt;
        if(TMath::Abs(fScaleShift) > DBL_EPSILON){
          effvec[1] += fScaleShift * effvec[1];
        }
      }
      fHistos->FillTHnSparse("hJetfindingEfficiency", effvec);
    }
  }
  return true;
}

Double_t AliAnalysisTaskEmcalJetEnergyScale::MakeAngularity(const AliEmcalJet &jet, AliVCluster::VCluUserDefEnergy_t energydef) const {
  if(!(jet.GetNumberOfTracks() || jet.GetNumberOfClusters()))
    throw AngularityException();
  TVector3 jetvec(jet.Px(), jet.Py(), jet.Pz());
  Double_t den(0.), num(0.);
  for(UInt_t itrk = 0; itrk < jet.GetNumberOfTracks(); itrk++) {
    auto track = jet.Track(itrk);
    if(!track){
      AliErrorStream() << "Associated constituent particle / track not found\n";
      continue;
    }
    TVector3 trackvec(track->Px(), track->Py(), track->Pz());

    num +=  track->Pt() * trackvec.DrEtaPhi(jetvec);
    den += +track->Pt();
  }
  for(UInt_t icl = 0; icl < jet.GetNumberOfClusters(); icl++){
    auto clust = jet.Cluster(icl);
    if(!clust) {
      AliErrorStream() << "Associated constituent cluster not found\n";
      continue;
    }
    TLorentzVector clusterp;
    clust->GetMomentum(clusterp, fVertex, energydef);

    num += clusterp.Pt() * clusterp.Vect().DrEtaPhi(jetvec);
    den += clusterp.Pt();
  }
  return num/den;
}

bool AliAnalysisTaskEmcalJetEnergyScale::IsSelectEmcalTriggers(const TString &triggerstring) const {
  const std::array<TString, 10> kEMCALTriggers = {
    "EMC7","EJE", "EJ1", "EJ2", "DJ1", "DJ2", "EG1", "EG2", "DG1", "DG2"
  };
  bool isEMCAL = false;
  for(auto emcaltrg : kEMCALTriggers) {
    if(triggerstring.Contains(emcaltrg)) {
      isEMCAL = true;
      break;
    }
  }
  return isEMCAL;
}

void AliAnalysisTaskEmcalJetEnergyScale::SetAngularitySpitting(bool doSplit) {
  if(doSplit) {
    auto partcont = GetPartLevelJetContainer();
    double radius = -1;
    if(partcont) radius = partcont->GetJetRadius();
    if(!fAngularityHandler) {
      if(radius > -1.) {
        AliInfoStream() << "Initializing angularity handler with R=" << radius << std::endl;
        fAngularityHandler = new AngularityHandler(radius);
      } else {
        AliInfoStream() << "Default initialization of angularity without radius" << std::endl;
        fAngularityHandler = new AngularityHandler();
      }
    }
  } else {
    if(fAngularityHandler) delete fAngularityHandler;
    fAngularityHandler = nullptr;
  }
}

void AliAnalysisTaskEmcalJetEnergyScale::ConfigurePtHard(MCProductionType_t mcprodtype, const TArrayI &pthardbinning, Bool_t doMCFilter, Double_t jetptcut) {
  SetMCProductionType(mcprodtype);
  SetUsePtHardBinScaling(true);
  SetUserPtHardBinning(pthardbinning);
  if(doMCFilter) {
    SetMCFilter();
    SetJetPtFactor(jetptcut);
  }
}

void AliAnalysisTaskEmcalJetEnergyScale::ConfigureMinBias(MCProductionType_t mcprodtype){
  if(!(mcprodtype == kMCPythiaMB || mcprodtype == kMCHepMCMB)) {
    AliErrorStream() << "MC prod type not compatible with min. bias production" << std::endl;
  }
  SetMCProductionType(mcprodtype);
}

void AliAnalysisTaskEmcalJetEnergyScale::ConfigureJetSelection(Double_t minJetPtPart, Double_t minJetPtDet, Double_t maxTrackPtPart, Double_t maxTrackPtDet, Double_t maxClusterPt, Double_t minAreaPerc) {
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

AliAnalysisTaskEmcalJetEnergyScale *AliAnalysisTaskEmcalJetEnergyScale::AddTaskJetEnergyScaleBase(AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recoscheme, AliVCluster::VCluUserDefEnergy_t energydef, Double_t jetradius, Bool_t useDCAL, const char *namepartcont, const char *trigger, const char *nametrackcont, const char *nameclustercont, const char *suffix) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    ::Error("PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetEnergyScale::AddTaskJetEnergyScale", "No analysis manager available");
    return nullptr;
  }

  auto inputhandler = mgr->GetInputEventHandler();
  auto isAOD = inputhandler->IsA() == AliAODInputHandler::Class();

  std::string jettypename;
  AliJetContainer::JetAcceptanceType acceptance(AliJetContainer::kTPCfid);
  AliJetContainer::EJetType_t mcjettype(jettype);
  bool addClusterContainer(false), addTrackContainer(false);
  switch(jettype){
    case AliJetContainer::kFullJet:
        jettypename = "FullJet";
        acceptance = useDCAL ? AliJetContainer::kDCALfid : AliJetContainer::kEMCALfid;
        addClusterContainer = addTrackContainer = true;
        break;
    case AliJetContainer::kChargedJet:
        jettypename = "ChargedJet";
        acceptance = AliJetContainer::kTPCfid;
        addTrackContainer = true;
        mcjettype = AliJetContainer::kChargedJet;
        break;
    case AliJetContainer::kNeutralJet:
        jettypename = "NeutralJet";
        acceptance = useDCAL ? AliJetContainer::kDCALfid : AliJetContainer::kEMCALfid;
        addClusterContainer = true;
        break;
    case AliJetContainer::kUndefinedJetType:
        break;
  };

  std::stringstream taskname, tag;
  tag << jettypename << "_R" << std::setw(2) << std::setfill('0') << int(jetradius * 10.) << "_" << trigger;
  if(strlen(suffix)) tag << "_" << suffix;
  taskname << "EnergyScaleTask_" << tag.str();
  AliAnalysisTaskEmcalJetEnergyScale *energyscaletask = new AliAnalysisTaskEmcalJetEnergyScale(taskname.str().data());
  mgr->AddTask(energyscaletask);
  energyscaletask->SetTriggerName(trigger);

  TString partcontname(namepartcont);
  if(partcontname == "usedefault") partcontname = "mcparticles";
  auto partcont = energyscaletask->AddMCParticleContainer(partcontname.Data());
  partcont->SetMinPt(0.);
  const std::string kNameMCParticles = "MCParticles";
  partcont->SetName(kNameMCParticles.data());
  energyscaletask->SetNameMCParticles(kNameMCParticles.data());

  AliClusterContainer *clusters(nullptr);
  TString clustercontname(nameclustercont);
  if(addClusterContainer) {
    const std::string kNameClusterContainer = "EMCALClusters";
    if(clustercontname == "usedefault") clusters = energyscaletask->AddClusterContainer(AliEmcalAnalysisFactory::ClusterContainerNameFactory(isAOD));
    else clusters = energyscaletask->AddClusterContainer(clustercontname.Data());
    clusters->SetDefaultClusterEnergy(energydef);
    clusters->SetClusUserDefEnergyCut(energydef, 0.3);
    clusters->SetName(kNameClusterContainer.data());
    energyscaletask->SetNameClusters(kNameClusterContainer.data());
  }
  AliTrackContainer *tracks(nullptr);
  TString trackcontname(nametrackcont);
  if(addTrackContainer) {
    const std::string kNameTrackContainer = "Tracks";
    if(trackcontname == "usedefault") tracks = energyscaletask->AddTrackContainer(AliEmcalAnalysisFactory::TrackContainerNameFactory(isAOD));
    else tracks = energyscaletask->AddTrackContainer(trackcontname.Data());
    tracks->SetName(kNameTrackContainer.data());
    energyscaletask->SetNameTracks(kNameTrackContainer.data());
  }

  const std::string kNameJetsPart = "particleLevelJets",
                    kNameJetsDet = "detectorLevelJets";

  auto contpartjet = energyscaletask->AddJetContainer(mcjettype, AliJetContainer::antikt_algorithm, recoscheme, jetradius,
                                                      acceptance, partcont, nullptr);
  contpartjet->SetName(kNameJetsPart.data());
  energyscaletask->SetNamePartJetContainer(kNameJetsPart.data());
  std::cout << "Adding particle-level jet container with underling array: " << contpartjet->GetArrayName() << std::endl;

  auto contdetjet = energyscaletask->AddJetContainer(jettype, AliJetContainer::antikt_algorithm, recoscheme, jetradius,
                                                     acceptance, tracks, clusters);
  contdetjet->SetName(kNameJetsDet.data());
  energyscaletask->SetNameDetJetContainer(kNameJetsDet.data());
  std::cout << "Adding detector-level jet container with underling array: " << contdetjet->GetArrayName() << std::endl;

  std::stringstream outnamebuilder, listnamebuilder;
  listnamebuilder << "EnergyScaleHists_" << tag.str();
  outnamebuilder << mgr->GetCommonFileName() << ":EnergyScaleResults_" << tag.str();

  mgr->ConnectInput(energyscaletask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(energyscaletask, 1, mgr->CreateContainer(listnamebuilder.str().data(), TList::Class(), AliAnalysisManager::kOutputContainer, outnamebuilder.str().data()));
  return energyscaletask;
}

AliAnalysisTaskEmcalJetEnergyScale *AliAnalysisTaskEmcalJetEnergyScale::AddTaskJetEnergyScaleBkgSubBase(AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recoscheme, AliVCluster::VCluUserDefEnergy_t energydef, Double_t jetradius, Bool_t useDCAL, const char *namepartcont, const char *nRho, const char *nRhoMC, const char *trigger, const char *nametrackcont, const char *nameclustercont, const char *suffix) {
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    ::Error("PWGJE::EMCALJetTasks::AliAnalysisTaskEmcalJetEnergyScale::AddTaskJetEnergyScale", "No analysis manager available");
    return nullptr;
  }

  auto inputhandler = mgr->GetInputEventHandler();
  auto isAOD = inputhandler->IsA() == AliAODInputHandler::Class();

  std::string jettypename;
  AliJetContainer::JetAcceptanceType acceptance(AliJetContainer::kTPCfid);
  AliJetContainer::EJetType_t mcjettype(jettype);
  bool addClusterContainer(false), addTrackContainer(false);
  switch(jettype){
    case AliJetContainer::kFullJet:
        jettypename = "FullJet";
        acceptance = useDCAL ? AliJetContainer::kDCALfid : AliJetContainer::kEMCALfid;
        addClusterContainer = addTrackContainer = true;
        break;
    case AliJetContainer::kChargedJet:
        jettypename = "ChargedJet";
        acceptance = AliJetContainer::kTPCfid;
        addTrackContainer = true;
        mcjettype = AliJetContainer::kChargedJet;
        break;
    case AliJetContainer::kNeutralJet:
        jettypename = "NeutralJet";
        acceptance = useDCAL ? AliJetContainer::kDCALfid : AliJetContainer::kEMCALfid;
        addClusterContainer = true;
        break;
    case AliJetContainer::kUndefinedJetType:
        break;
  };

  std::stringstream taskname, tag;
  tag << jettypename << "_R" << std::setw(2) << std::setfill('0') << int(jetradius * 10.) << "_" << trigger;
  if(strlen(suffix)) tag << "_" << suffix;
  taskname << "EnergyScaleTask_" << tag.str();
  AliAnalysisTaskEmcalJetEnergyScale *energyscaletask = new AliAnalysisTaskEmcalJetEnergyScale(taskname.str().data());
  mgr->AddTask(energyscaletask);
  energyscaletask->SetTriggerName(trigger);

  TString partcontname(namepartcont);
  if(partcontname == "usedefault") partcontname = "mcparticles";
  auto partcont = energyscaletask->AddMCParticleContainer(partcontname.Data());
  partcont->SetMinPt(0.);
  const std::string kNameMCParticles = "MCParticles";
  partcont->SetName(kNameMCParticles.data());
  energyscaletask->SetNameMCParticles(kNameMCParticles.data());

  AliClusterContainer *clusters(nullptr);
  TString clustercontname(nameclustercont);
  if(addClusterContainer) {
    const std::string kNameClusterContainer = "EMCALClusters";
    if(clustercontname == "usedefault") clusters = energyscaletask->AddClusterContainer(AliEmcalAnalysisFactory::ClusterContainerNameFactory(isAOD));
    else clusters = energyscaletask->AddClusterContainer(clustercontname.Data());
    clusters->SetDefaultClusterEnergy(energydef);
    clusters->SetClusUserDefEnergyCut(energydef, 0.3);
    clusters->SetName(kNameClusterContainer.data());
    energyscaletask->SetNameClusters(kNameClusterContainer.data());
  }
  AliTrackContainer *tracks(nullptr);
  TString trackcontname(nametrackcont);
  if(addTrackContainer) {
    const std::string kNameTrackContainer = "Tracks";
    if(trackcontname == "usedefault") tracks = energyscaletask->AddTrackContainer(AliEmcalAnalysisFactory::TrackContainerNameFactory(isAOD));
    else tracks = energyscaletask->AddTrackContainer(trackcontname.Data());
    tracks->SetName(kNameTrackContainer.data());
    energyscaletask->SetNameTracks(kNameTrackContainer.data());
  }

  const std::string kNameJetsPart = "particleLevelJets",
                    kNameJetsDet = "detectorLevelJets";

  auto contpartjet = energyscaletask->AddJetContainer(mcjettype, AliJetContainer::antikt_algorithm, recoscheme, jetradius,
                                                      acceptance, partcont, nullptr);
  contpartjet->SetName(kNameJetsPart.data());
  contpartjet->SetRhoName(nRhoMC);
  energyscaletask->SetNamePartJetContainer(kNameJetsPart.data());
  std::cout << "Adding particle-level jet container with underling array: " << contpartjet->GetArrayName() << std::endl;

  auto contdetjet = energyscaletask->AddJetContainer(jettype, AliJetContainer::antikt_algorithm, recoscheme, jetradius,
                                                     acceptance, tracks, clusters);
  contdetjet->SetName(kNameJetsDet.data());
  contdetjet->SetRhoName(nRho);
  energyscaletask->SetNameDetJetContainer(kNameJetsDet.data());
  std::cout << "Adding detector-level jet container with underling array: " << contdetjet->GetArrayName() << std::endl;

  std::stringstream outnamebuilder, listnamebuilder;
  listnamebuilder << "EnergyScaleHists_" << tag.str();
  outnamebuilder << mgr->GetCommonFileName() << ":EnergyScaleResults_" << tag.str();

  mgr->ConnectInput(energyscaletask, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(energyscaletask, 1, mgr->CreateContainer(listnamebuilder.str().data(), TList::Class(), AliAnalysisManager::kOutputContainer, outnamebuilder.str().data()));
  return energyscaletask;
}


void AngularityHandler::InitFromFile(const char *filename) {
  fBins.Clear();
  std::unique_ptr<TFile> reader(TFile::Open(filename, "READ"));
  std::string rstring=Form("R%02d", int(fRadius*10.));
  auto hist = dynamic_cast<TH1 *>(reader->Get(Form("meanAngularity%s", rstring.data())));
  if(!hist) {
    AliErrorStream() << "Mean angularity dist not found for R = " << fRadius << std::endl;
    return;
  } else {
    AliDebugStream(1) << "Using histogram " << hist->GetName() << " for R=" << fRadius << std::endl;
  }
  double currentptmax = DBL_MIN, currentval = DBL_MIN; // For overflow bin (using same value as of the highest bin)
  for(int ib = 0; ib < hist->GetXaxis()->GetNbins(); ib++){
    auto min = hist->GetXaxis()->GetBinLowEdge(ib+1),
         max = hist->GetXaxis()->GetBinUpEdge(ib+1),
         value = hist->GetBinContent(ib+1);
    SetBin(min, max, value);
    if(max > currentptmax) {
      currentval = value;
      currentptmax = max;
    }
  }
  // Set overflow bin
  SetBin(currentptmax, DBL_MAX, currentval);
}

void AngularityHandler::SetBin(double min, double max, double value) {
  bool hasOverlap = false;
  for(auto b : fBins) {
    AngularityBin *nextbin = dynamic_cast<AngularityBin *>(b);
    if(!nextbin) continue;
    if(min >= nextbin->Min() && min < nextbin->Max()) hasOverlap = true;
    if(max >= nextbin->Min() && max < nextbin->Max()) hasOverlap = true;
    if(hasOverlap) break;
  }
  if(hasOverlap) {
    AliErrorStream() << "Cannot add bin with min " << min << " and max " << max << " due to overlap";
    return;
  }
  fBins.Add(new AngularityBin(min, max, value));
}

double AngularityHandler::GetValue(double pt) const {
  auto result = FindBin(pt);
  return result.Value();
}

bool AngularityHandler::IsHigherAngularity(double pt, double angularity) const {
  auto result = FindBin(pt);
  if(!result.IsOK()) throw BinNotFoundException(pt);
  return angularity > result.Value();
}

AngularityHandler::AngularityBin AngularityHandler::FindBin(double pt) const {
  AngularityBin result;
  for(auto b : fBins) {
    AngularityBin *nextbin = dynamic_cast<AngularityBin *>(b);
    if(!nextbin) continue;
    if(nextbin->IsInRange(pt)) {
      result = *nextbin;
      AliDebugStream(3) << "Found bin with minb " << result.Min() << " and maxb " << result.Max() << " for pt = " << pt << std::endl;
      break;
    }
  }
  return result;
}

void AngularityHandler::PrintStream(std::ostream &stream) const {
  stream << "Angularity Handler for R = " << fRadius << std::endl;
  stream << "===========================================" << std::endl;
  for(const auto binObject : fBins) {
    auto bin = static_cast<AngularityBin *>(binObject);
    if(bin) {
      stream << *bin << std::endl;
    }
  }
}

void AngularityHandler::AngularityBin::PrintStream(std::ostream &stream) const {
  stream << "min " << fMin << ", max " << fMax << ", value " << fValue;
}

void AngularityHandler::Print(Option_t *) const {
  std::cout << *this << std::endl;
}

std::ostream &operator<<(std::ostream &stream, const PWGJE::EMCALJetTasks::AngularityHandler &handler) {
  handler.PrintStream(stream);
  return stream;
}

std::ostream &operator<<(std::ostream &stream, const PWGJE::EMCALJetTasks::AngularityHandler::AngularityBin &bin) {
  bin.PrintStream(stream);
  return stream;
}

AliAnalysisTaskEmcalJetEnergyScale *AliAnalysisTaskEmcalJetEnergyScale::AddTaskJetEnergyScale(AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recoscheme, AliVCluster::VCluUserDefEnergy_t energydef, Double_t jetradius, Bool_t useDCAL, const char *namepartcont, const char *trigger, const char *suffix) {
  return AliAnalysisTaskEmcalJetEnergyScale::AddTaskJetEnergyScaleBase(jettype, recoscheme, energydef, jetradius, useDCAL, namepartcont, trigger, "usedefault", "usedefault", suffix);
}

AliAnalysisTaskEmcalJetEnergyScale *AliAnalysisTaskEmcalJetEnergyScale::AddTaskJetEnergyScale(AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recoscheme, AliVCluster::VCluUserDefEnergy_t energydef, Double_t jetradius, Bool_t useDCAL, const char *namepartcont, const char *trigger, const char *nametrackcont, const char *nameclustercont, const char *suffix) {
  return AliAnalysisTaskEmcalJetEnergyScale::AddTaskJetEnergyScaleBase(jettype, recoscheme, energydef, jetradius, useDCAL, namepartcont, trigger, nametrackcont, nameclustercont, suffix);
}

AliAnalysisTaskEmcalJetEnergyScale *AliAnalysisTaskEmcalJetEnergyScale::AddTaskJetEnergyScaleBkgSub(AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recoscheme, AliVCluster::VCluUserDefEnergy_t energydef, Double_t jetradius, Bool_t useDCAL, const char *namepartcont, const char *nRho, const char *nRhoMC, const char *trigger, const char *suffix) {
  return AliAnalysisTaskEmcalJetEnergyScale::AddTaskJetEnergyScaleBkgSubBase(jettype, recoscheme, energydef, jetradius, useDCAL, namepartcont, nRho, nRhoMC, trigger, "usedefault", "usedefault", suffix);
}

AliAnalysisTaskEmcalJetEnergyScale *AliAnalysisTaskEmcalJetEnergyScale::AddTaskJetEnergyScaleBkgSub(AliJetContainer::EJetType_t jettype, AliJetContainer::ERecoScheme_t recoscheme, AliVCluster::VCluUserDefEnergy_t energydef, Double_t jetradius, Bool_t useDCAL, const char *namepartcont, const char *nRho, const char *nRhoMC, const char *trigger, const char *nametrackcont, const char *nameclustercont, const char *suffix) {
  return AliAnalysisTaskEmcalJetEnergyScale::AddTaskJetEnergyScaleBkgSubBase(jettype, recoscheme, energydef, jetradius, useDCAL, namepartcont, nRho, nRhoMC, trigger, nametrackcont, nameclustercont, suffix);
}
