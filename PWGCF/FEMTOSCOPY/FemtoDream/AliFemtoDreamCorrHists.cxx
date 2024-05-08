/*
 * AliFemtoDreamCorrHists.cxx
 *
 *  Created on: Sep 12, 2017
 *      Author: gu74req
 */
#include <iostream>
#include <vector>
#include "AliFemtoDreamCorrHists.h"
#include "TMath.h"
#include "AliLog.h"
#include "TDatabasePDG.h"

ClassImp(AliFemtoDreamCorrHists)
    AliFemtoDreamCorrHists::AliFemtoDreamCorrHists()
    : fQA(nullptr),
      fResults(nullptr),
      fPairs(nullptr),
      fPairQA(nullptr),
      fMinimalBooking(false),
      fMomentumResolution(false),
      fPhiEtaPlots(false),
      fRelKThreshold(0.2),
      fSameEventDist(nullptr),
      fSameEventMinvKtandRelativeKDist(nullptr),
      fSameEventMultDist(nullptr),
      fSameEventCentDist(nullptr),
      fSameEventmTDist(nullptr),
      fSameEventmTvsMultDist(nullptr),
      fSameEventkTDist(nullptr),
      fSameEventkTandMultDist(nullptr),
      fSameEventkTandMultPtDist(nullptr),
      fSameEventkTandMultMCTrueDist(nullptr),
      fSameEventkTCentDist(nullptr),
      fSameEventmTMultDist(nullptr),
      fPtQADist(nullptr),
      fPtQADistSEPartOne(nullptr),
      fPtQADistSEPartTwo(nullptr),
      fPtQADistMEPartOne(nullptr),
      fPtQADistMEPartTwo(nullptr),
      fKstarPtQADistSEPartOne(nullptr),
      fKstarPtQADistSEPartTwo(nullptr),
      fKstarPtQADistMEPartOne(nullptr),
      fKstarPtQADistMEPartTwo(nullptr),
      fMassQADistPart1(nullptr),
      fMassQADistPart2(nullptr),
      fPairInvMassQAD(nullptr),
      fPairInvMassKstarQAD(nullptr),
      fMEMassQADistPart1(nullptr),
      fMEMassQADistPart2(nullptr),
      fPairInvMEMassQAD(nullptr),
      fPairInvMEMassKstarQAD(nullptr),
      fPairCounterSE(nullptr),
      fMixedEventDist(nullptr),
      fMixedEventMinvKtandRelativeKDist(nullptr),
      fMixedEventMultDist(nullptr),
      fMixedEventCentDist(nullptr),
      fMixedEventmTDist(nullptr),
      fMixedEventmTvsMultDist(nullptr),
      fMixedEventkTDist(nullptr),
      fMixedEventkTandMultDist(nullptr),
      fMixedEventkTandMultMCTrueDist(nullptr),
      fMixedEventkTCentDist(nullptr),
      fMixedEventmTMultDist(nullptr),
      fPairCounterME(nullptr),
      fMomResolutionSE(nullptr),
      fMomResolutionSEAll(nullptr),
      fMomResolutionME(nullptr),
      fMomResolutionMEAll(nullptr),
      fMomResolutionDist(nullptr),
      fRadiiEtaPhiSE(nullptr),
      fRadiiEtaPhiME(nullptr),
      fIntRadiiQAEtaPhiSEBefore(nullptr),
      fIntRadiiQAEtaPhiMEBefore(nullptr),
      fIntRadiiQAEtaPhiSEAfter(nullptr),
      fIntRadiiQAEtaPhiMEAfter(nullptr),
      fRadiiEtaPhiSEsmallK(nullptr),
      fRadiiEtaPhiMEsmallK(nullptr),
      fdEtadPhiSE(nullptr),
      fdEtadPhiME(nullptr),
      fdEtadPhiSEmT(nullptr),
      fdEtadPhiMEmT(nullptr),
      fEffMixingDepth(nullptr),
      fSameEventDistCommon(nullptr),
      fSameEventDistNonCommon(nullptr),
      fdEtadPhiSECommon(nullptr),
      fdEtadPhiSENonCommon(nullptr),
      fSameEventMultDistCommon(nullptr),
      fSameEventMultDistNonCommon(nullptr),
      fSameEventmTDistCommon(nullptr),
      fSameEventmTDistNonCommon(nullptr),
      fSameEventmTMultDistCommon(nullptr),
      fSameEventmTMultDistNonCommon(nullptr),
      fSameEventpTOnepTTwokStar(nullptr),
      fMixedEventpTOnepTTwokStar(nullptr),
      fDoMultBinning(false),
      fDoCentBinning(false),
      fDokTandMultBinning(false),
      fDokTandMultPtBinning(false),
      fDokTandMultMCTrueBinning(false),
      fDokTBinning(false),
      fDomTBinning(false),
      fmTMultPlots(false),
      fPtQA(false),
      fMassQA(false),
      fDokTCentralityBins(false),
      fdPhidEtaPlots(false),
      fPhiEtaPlotsSmallK(false),
      fmTDetaDPhi(false),
      fAncestors(false),
      fRemoveAncestorsResonances(false),
      fpTOnepTTwokStarPlotsmT(false),
      fpTOnepTTwokStarCutOff(3.),
      fPDGCode(),
      fmTBins(),
      fWhichPairs(),
      fCentBins(0)
{
}

AliFemtoDreamCorrHists::AliFemtoDreamCorrHists(
    const AliFemtoDreamCorrHists &hists)
    : fQA(hists.fQA),
      fResults(hists.fResults),
      fPairs(hists.fPairs),
      fPairQA(hists.fPairQA),
      fMinimalBooking(hists.fMinimalBooking),
      fMomentumResolution(hists.fMomResolutionME),
      fPhiEtaPlots(hists.fPhiEtaPlots),
      fRelKThreshold(hists.fRelKThreshold),
      fSameEventDist(hists.fSameEventDist),
      fSameEventMinvKtandRelativeKDist(hists.fSameEventMinvKtandRelativeKDist),
      fSameEventMultDist(hists.fSameEventMultDist),
      fSameEventCentDist(hists.fSameEventCentDist),
      fSameEventmTDist(hists.fSameEventmTDist),
      fSameEventmTvsMultDist(hists.fSameEventmTvsMultDist),
      fSameEventkTDist(hists.fSameEventkTDist),
      fSameEventkTandMultDist(hists.fSameEventkTandMultDist),
      fSameEventkTandMultPtDist(hists.fSameEventkTandMultPtDist),
      fSameEventkTandMultMCTrueDist(hists.fSameEventkTandMultMCTrueDist),
      fSameEventkTCentDist(hists.fSameEventkTCentDist),
      fSameEventmTMultDist(hists.fSameEventmTMultDist),
      fPtQADist(hists.fPtQADist),
      fPtQADistSEPartOne(hists.fPtQADistSEPartOne),
      fPtQADistSEPartTwo(hists.fPtQADistSEPartTwo),
      fPtQADistMEPartOne(hists.fPtQADistMEPartOne),
      fPtQADistMEPartTwo(hists.fPtQADistMEPartTwo),
      fKstarPtQADistSEPartOne(hists.fKstarPtQADistSEPartOne),
      fKstarPtQADistSEPartTwo(hists.fKstarPtQADistSEPartTwo),
      fKstarPtQADistMEPartOne(hists.fKstarPtQADistMEPartOne),
      fKstarPtQADistMEPartTwo(hists.fKstarPtQADistMEPartTwo),
      fMassQADistPart1(hists.fMassQADistPart1),
      fMassQADistPart2(hists.fMassQADistPart2),
      fPairInvMassQAD(hists.fPairInvMassQAD),
      fPairInvMassKstarQAD(hists.fPairInvMassKstarQAD),
      fMEMassQADistPart1(hists.fMEMassQADistPart1),
      fMEMassQADistPart2(hists.fMEMassQADistPart2),
      fPairInvMEMassQAD(hists.fPairInvMEMassQAD),
      fPairInvMEMassKstarQAD(hists.fPairInvMEMassKstarQAD),
      fPairCounterSE(hists.fPairCounterSE),
      fMixedEventDist(hists.fMixedEventDist),
      fMixedEventMinvKtandRelativeKDist(hists.fMixedEventMinvKtandRelativeKDist),
      fMixedEventMultDist(hists.fMixedEventMultDist),
      fMixedEventCentDist(hists.fMixedEventCentDist),
      fMixedEventmTDist(hists.fMixedEventmTDist),
      fMixedEventmTvsMultDist(hists.fMixedEventmTvsMultDist),
      fMixedEventkTDist(hists.fMixedEventkTDist),
      fMixedEventkTandMultDist(hists.fMixedEventkTandMultDist),
      fMixedEventkTandMultMCTrueDist(hists.fMixedEventkTandMultMCTrueDist),
      fMixedEventkTCentDist(hists.fMixedEventkTCentDist),
      fMixedEventmTMultDist(hists.fMixedEventmTMultDist),
      fPairCounterME(hists.fPairCounterME),
      fMomResolutionSE(hists.fMomResolutionSE),
      fMomResolutionSEAll(hists.fMomResolutionSEAll),
      fMomResolutionME(hists.fMomResolutionME),
      fMomResolutionMEAll(hists.fMomResolutionMEAll),
      fMomResolutionDist(hists.fMomResolutionDist),
      fRadiiEtaPhiSE(hists.fRadiiEtaPhiSE),
      fRadiiEtaPhiME(hists.fRadiiEtaPhiME),
      fIntRadiiQAEtaPhiSEBefore(hists.fIntRadiiQAEtaPhiSEBefore),
      fIntRadiiQAEtaPhiMEBefore(hists.fIntRadiiQAEtaPhiMEBefore),
      fIntRadiiQAEtaPhiSEAfter(hists.fIntRadiiQAEtaPhiSEAfter),
      fIntRadiiQAEtaPhiMEAfter(hists.fIntRadiiQAEtaPhiMEAfter),
      fRadiiEtaPhiSEsmallK(hists.fRadiiEtaPhiSEsmallK),
      fRadiiEtaPhiMEsmallK(hists.fRadiiEtaPhiMEsmallK),
      fdEtadPhiSE(hists.fdEtadPhiSE),
      fdEtadPhiME(hists.fdEtadPhiME),
      fdEtadPhiSEmT(hists.fdEtadPhiSEmT),
      fdEtadPhiMEmT(hists.fdEtadPhiMEmT),
      fEffMixingDepth(hists.fEffMixingDepth),
      fSameEventDistCommon(hists.fSameEventDistCommon),
      fSameEventDistNonCommon(hists.fSameEventDistCommon),
      fdEtadPhiSECommon(hists.fdEtadPhiSECommon),
      fdEtadPhiSENonCommon(hists.fdEtadPhiSENonCommon),
      fSameEventMultDistCommon(hists.fSameEventMultDistCommon),
      fSameEventMultDistNonCommon(hists.fSameEventMultDistNonCommon),
      fSameEventmTDistCommon(hists.fSameEventmTDistCommon),
      fSameEventmTDistNonCommon(hists.fSameEventmTDistNonCommon),
      fSameEventmTMultDistCommon(hists.fSameEventmTMultDistCommon),
      fSameEventmTMultDistNonCommon(hists.fSameEventmTMultDistNonCommon),
      fSameEventpTOnepTTwokStar(hists.fSameEventpTOnepTTwokStar),
      fMixedEventpTOnepTTwokStar(hists.fMixedEventpTOnepTTwokStar),
      fDoMultBinning(hists.fDoMultBinning),
      fDoCentBinning(hists.fDoCentBinning),
      fDokTandMultBinning(hists.fDokTandMultBinning),
      fDokTandMultPtBinning(hists.fDokTandMultPtBinning),
      fDokTandMultMCTrueBinning(hists.fDokTandMultMCTrueBinning),
      fDokTBinning(hists.fDokTBinning),
      fDomTBinning(hists.fDomTBinning),
      fmTMultPlots(hists.fmTMultPlots),
      fPtQA(hists.fPtQA),
      fMassQA(hists.fMassQA),
      fDokTCentralityBins(hists.fDokTCentralityBins),
      fdPhidEtaPlots(hists.fdPhidEtaPlots),
      fPhiEtaPlotsSmallK(hists.fPhiEtaPlotsSmallK),
      fmTDetaDPhi(hists.fmTDetaDPhi),
      fAncestors(hists.fAncestors),
      fRemoveAncestorsResonances(hists.fRemoveAncestorsResonances),
      fpTOnepTTwokStarPlotsmT(hists.fpTOnepTTwokStarPlotsmT),
      fpTOnepTTwokStarCutOff(hists.fpTOnepTTwokStarCutOff),
      fPDGCode(hists.fPDGCode),
      fmTBins(hists.fmTBins),
      fWhichPairs(hists.fWhichPairs),
      fCentBins(hists.fCentBins)
{
}

AliFemtoDreamCorrHists::AliFemtoDreamCorrHists(AliFemtoDreamCollConfig *conf,
                                               bool MinimalBooking)
    : fQA(nullptr),
      fResults(nullptr),
      fPairs(nullptr),
      fPairQA(nullptr),
      fMinimalBooking(MinimalBooking),
      fMomentumResolution(false),
      fPhiEtaPlots(false),
      fRelKThreshold(0.2),
      fSameEventDist(nullptr),
      fSameEventMinvKtandRelativeKDist(nullptr),
      fSameEventMultDist(nullptr),
      fSameEventCentDist(nullptr),
      fSameEventmTDist(nullptr),
      fSameEventmTvsMultDist(nullptr),
      fSameEventkTDist(nullptr),
      fSameEventkTandMultDist(nullptr),
      fSameEventkTandMultPtDist(nullptr),
      fSameEventkTandMultMCTrueDist(nullptr),
      fSameEventkTCentDist(nullptr),
      fSameEventmTMultDist(nullptr),
      fPtQADist(nullptr),
      fPtQADistSEPartOne(nullptr),
      fPtQADistSEPartTwo(nullptr),
      fPtQADistMEPartOne(nullptr),
      fPtQADistMEPartTwo(nullptr),
      fKstarPtQADistSEPartOne(nullptr),
      fKstarPtQADistSEPartTwo(nullptr),
      fKstarPtQADistMEPartOne(nullptr),
      fKstarPtQADistMEPartTwo(nullptr),
      fMassQADistPart1(nullptr),
      fMassQADistPart2(nullptr),
      fPairInvMassQAD(nullptr),
      fPairInvMassKstarQAD(nullptr),
      fMEMassQADistPart1(nullptr),
      fMEMassQADistPart2(nullptr),
      fPairInvMEMassQAD(nullptr),
      fPairInvMEMassKstarQAD(nullptr),
      fPairCounterSE(nullptr),
      fMixedEventDist(nullptr),
      fMixedEventMinvKtandRelativeKDist(nullptr),
      fMixedEventMultDist(nullptr),
      fMixedEventCentDist(nullptr),
      fMixedEventmTDist(nullptr),
      fMixedEventmTvsMultDist(nullptr),
      fMixedEventkTDist(nullptr),
      fMixedEventkTandMultDist(nullptr),
      fMixedEventkTandMultMCTrueDist(nullptr),
      fMixedEventkTCentDist(nullptr),
      fMixedEventmTMultDist(nullptr),
      fPairCounterME(nullptr),
      fMomResolutionSE(nullptr),
      fMomResolutionSEAll(nullptr),
      fMomResolutionME(nullptr),
      fMomResolutionMEAll(nullptr),
      fMomResolutionDist(nullptr),
      fRadiiEtaPhiSE(nullptr),
      fRadiiEtaPhiME(nullptr),
      fIntRadiiQAEtaPhiSEBefore(nullptr),
      fIntRadiiQAEtaPhiMEBefore(nullptr),
      fIntRadiiQAEtaPhiSEAfter(nullptr),
      fIntRadiiQAEtaPhiMEAfter(nullptr),
      fRadiiEtaPhiSEsmallK(nullptr),
      fRadiiEtaPhiMEsmallK(nullptr),
      fdEtadPhiSE(nullptr),
      fdEtadPhiME(nullptr),
      fdEtadPhiSEmT(nullptr),
      fdEtadPhiMEmT(nullptr),
      fEffMixingDepth(nullptr),
      fSameEventDistCommon(nullptr),
      fSameEventDistNonCommon(nullptr),
      fdEtadPhiSECommon(nullptr),
      fdEtadPhiSENonCommon(nullptr),
      fSameEventMultDistCommon(nullptr),
      fSameEventMultDistNonCommon(nullptr),
      fSameEventmTDistCommon(nullptr),
      fSameEventmTDistNonCommon(nullptr),
      fSameEventmTMultDistCommon(nullptr),
      fSameEventmTMultDistNonCommon(nullptr),
      fSameEventpTOnepTTwokStar(nullptr),
      fMixedEventpTOnepTTwokStar(nullptr),
      fDoMultBinning(false),
      fDoCentBinning(false),
      fDokTandMultBinning(false),
      fDokTandMultPtBinning(false),
      fDokTandMultMCTrueBinning(false),
      fDokTBinning(false),
      fDomTBinning(false),
      fmTMultPlots(false),
      fPtQA(false),
      fMassQA(false),
      fDokTCentralityBins(false),
      fdPhidEtaPlots(false),
      fPhiEtaPlotsSmallK(false),
      fmTDetaDPhi(false),
      fAncestors(false),
      fRemoveAncestorsResonances(false),
      fpTOnepTTwokStarPlotsmT(false),
      fpTOnepTTwokStarCutOff(3.),
      fPDGCode(),
      fmTBins(),
      fWhichPairs(),
      fCentBins(0)
{
  //  fMinimalBooking=MinimalBooking;
  fMomentumResolution = conf->GetDoMomResolution();
  fDoMultBinning = conf->GetDoMultBinning();
  fDoCentBinning = conf->GetDoCentBinning();
  fDokTBinning = conf->GetDokTBinning();
  fDokTCentralityBins = conf->GetDokTCentralityBinning();
  fDokTandMultBinning = conf->GetDokTandMultBinning();
  fDokTandMultPtBinning = conf->GetDokTandMultPtBinning();
  fDokTandMultMCTrueBinning = conf->GetDokTandMultMCTrueBinning();
  fDoMinvKtandRelativeKBinning = conf->GetMinvKtandRelativeKBinning();
  fDomTBinning = conf->GetDomTBinning();
  fmTMultPlots = conf->GetmTMultBinning();
  fPtQA = conf->GetDoPtQA();
  fMassQA = conf->GetDoMassQA();
  fPDGCode = conf->GetPDGCodes();
  fPhiEtaPlots = conf->GetDoPhiEtaBinning();
  fdPhidEtaPlots = conf->GetdPhidEtaPlots();
  fPhiEtaPlotsSmallK = conf->GetdPhidEtaPlotsSmallK();
  fmTDetaDPhi = conf->GetdPhidEtamTPlots();
  fAncestors = conf->GetDoAncestorsPlots();
  fRemoveAncestorsResonances = conf->GetRemoveAncestorResonances();
  fpTOnepTTwokStarPlotsmT = conf->GetDopTOnepTTwokStarPlotsmT();
  fpTOnepTTwokStarCutOff = conf->GetDopTOnepTTwokStarCutOff();

  if (fDokTCentralityBins && !fDokTBinning)
  {
    AliWarning(
        "Doing the Centrality binning without the kT Binning wont work!\n");
  }
  if (fDokTCentralityBins)
  {
    fCentBins = conf->GetCentBins();
    if (fCentBins.size() == 0)
    {
      AliWarning("Did you set the Centrality bins, since their size is 0?\n");
    }
  }
  int multbins = conf->GetNMultBins();
  int nParticles = conf->GetNParticles();
  const unsigned int nHists = conf->GetNParticleCombinations();
  std::vector<int> NBinsHist = conf->GetNBinsHist();
  std::vector<int>::iterator itNBins = NBinsHist.begin();
  std::vector<float> kRelMin = conf->GetMinKRel();
  std::vector<float>::iterator itKMin = kRelMin.begin();
  std::vector<float> kRelMax = conf->GetMaxKRel();
  std::vector<float>::iterator itKMax = kRelMax.begin();
  fmTBins = conf->GetmTBins();
  const unsigned int nmTBins = fmTBins.size();
  if (fDokTandMultBinning)
  {
    if (multbins == 0)
    {
      AliWarning("Did you set the Multiplicity bins, since their size is 0?\n");
    }
  }
  if (fDokTandMultMCTrueBinning)
  {
    if (multbins == 0)
    {
      AliWarning("Did you set the Multiplicity bins, since their size is 0?\n");
    }
  }
  if (nHists != NBinsHist.size() || nHists != kRelMin.size() || nHists != kRelMax.size())
  {
    // Todo: Replace by AliFatal!
    std::cout << "Initialzing the bin sizes, something went horribly wrong!"
              << std::endl;
  }
  // The way the histograms are assigned later is going to be for example for
  // 4 different particle species X1,X2,X3,X4:
  //     X1  X2  X3  X4
  // X1  1   2   3   4
  // X2      5   6   7
  // X3          8   9
  // X4              10<-----Number of the Histogram=Position in input vector
  // X1 corresponds the first particle array in your vector that you hand over
  // in the AliFemtoDreamPartCollection::SetEvent Method, X2 to the second and
  // so on.
  // in case we only book the most basic things we don't need this

  TDatabasePDG::Instance()->AddParticle("deuteron", "deuteron", 1.8756134,
                                        kTRUE, 0.0, 1, "Nucleus", 1000010020);
  TDatabasePDG::Instance()->AddAntiParticle("anti-deuteron", -1000010020);

  fWhichPairs = conf->GetWhichPairs();
  fResults = new TList();
  fResults->SetName("Results");
  fResults->SetOwner();

  if (!fMinimalBooking)
  {
    fQA = new TList();
    fQA->SetName("PairQA");
    fQA->SetOwner();

    fPairQA = new TList *[nHists];
    fPairCounterSE = new TH2F *[nHists];
    fPairCounterME = new TH2F *[nHists];
    fEffMixingDepth = new TH1F *[nHists];
    if (fMomentumResolution)
    {
      fMomResolutionSE = new TH2F *[nHists];
      fMomResolutionSEAll = new TH2F *[nHists];
      fMomResolutionME = new TH2F *[nHists];
      fMomResolutionMEAll = new TH2F *[nHists];
      fMomResolutionDist = new TH2F *[nHists];
    }
    else
    {
      fMomResolutionSE = nullptr;
      fMomResolutionSEAll = nullptr;
      fMomResolutionME = nullptr;
      fMomResolutionMEAll = nullptr;
      fMomResolutionDist = nullptr;
    }
    if (fPhiEtaPlots)
    {
      fRadiiEtaPhiSE = new TH2F ***[nHists];
      fRadiiEtaPhiME = new TH2F ***[nHists];
      fIntRadiiQAEtaPhiSEBefore = new TH2F **[nHists];
      fIntRadiiQAEtaPhiMEBefore = new TH2F **[nHists];
      fIntRadiiQAEtaPhiSEAfter = new TH2F **[nHists];
      fIntRadiiQAEtaPhiMEAfter = new TH2F **[nHists];
      if (fPhiEtaPlotsSmallK)
      {
        fRadiiEtaPhiSEsmallK = new TH2F ***[nHists];
        fRadiiEtaPhiMEsmallK = new TH2F ***[nHists];
      }
    }
    else
    {
      fRadiiEtaPhiSE = nullptr;
      fRadiiEtaPhiME = nullptr;
      fIntRadiiQAEtaPhiSEBefore = nullptr;
      fIntRadiiQAEtaPhiMEBefore = nullptr;
      fIntRadiiQAEtaPhiSEAfter = nullptr;
      fIntRadiiQAEtaPhiMEAfter = nullptr;
      fRadiiEtaPhiSEsmallK = nullptr;
      fRadiiEtaPhiMEsmallK = nullptr;
    }
  }
  else
  {
    fQA = nullptr;
    fPairQA = nullptr;
    fPairCounterSE = nullptr;
    fPairCounterME = nullptr;
    fEffMixingDepth = nullptr;
    fMomResolutionSE = nullptr;
    fMomResolutionSEAll = nullptr;
    fMomResolutionME = nullptr;
    fMomResolutionMEAll = nullptr;
    fMomResolutionDist = nullptr;
    fRadiiEtaPhiSE = nullptr;
    fRadiiEtaPhiME = nullptr;
    fIntRadiiQAEtaPhiSEBefore = nullptr;
    fIntRadiiQAEtaPhiMEBefore = nullptr;
    fIntRadiiQAEtaPhiSEAfter = nullptr;
    fIntRadiiQAEtaPhiMEAfter = nullptr;
    fRadiiEtaPhiSEsmallK = nullptr;
    fRadiiEtaPhiMEsmallK = nullptr;
  }
  // we always want to do this, regardless of the booking type!
  fPairs = new TList *[nHists];
  fSameEventDist = new TH1F *[nHists];
  fMixedEventDist = new TH1F *[nHists];
  if (fDoMinvKtandRelativeKBinning)
  {
    fSameEventMinvKtandRelativeKDist = new TH2F *[nHists];
    fMixedEventMinvKtandRelativeKDist = new TH2F *[nHists];
  }
  else
  {
    fSameEventMinvKtandRelativeKDist = nullptr;
    fMixedEventMinvKtandRelativeKDist = nullptr;
  }
  if (fDoMultBinning)
  {
    fSameEventMultDist = new TH2F *[nHists];
    fMixedEventMultDist = new TH2F *[nHists];
  }
  else
  {
    fSameEventMultDist = nullptr;
    fMixedEventMultDist = nullptr;
  }
  if (fDokTandMultBinning)
  {
    fSameEventkTandMultDist = new TH2F **[nHists];
    fMixedEventkTandMultDist = new TH2F **[nHists];
    if (fDokTandMultPtBinning)
    {
      fSameEventkTandMultPtDist = new TH2F **[nHists];
    }
  }
  else
  {
    fSameEventkTandMultDist = nullptr;
    fSameEventkTandMultPtDist = nullptr;
    fMixedEventkTandMultDist = nullptr;
  }
  if (fDokTandMultMCTrueBinning)
  {
    fSameEventkTandMultMCTrueDist = new TH2F **[nHists];
    fMixedEventkTandMultMCTrueDist = new TH2F **[nHists];
  }
  else
  {
    fSameEventkTandMultMCTrueDist = nullptr;
    fMixedEventkTandMultMCTrueDist = nullptr;
  }
  if (fDoCentBinning)
  {
    fSameEventCentDist = new TH2F *[nHists];
    fMixedEventCentDist = new TH2F *[nHists];
  }
  else
  {
    fSameEventCentDist = nullptr;
    fMixedEventCentDist = nullptr;
  }
  if (fDokTBinning)
  {
    fSameEventkTDist = new TH2F *[nHists];
    fMixedEventkTDist = new TH2F *[nHists];
  }
  else
  {
    fSameEventkTDist = nullptr;
    fMixedEventkTDist = nullptr;
  }
  if (fDokTCentralityBins)
  {
    fSameEventkTCentDist = new TH2F **[nHists];
    fMixedEventkTCentDist = new TH2F **[nHists];
  }
  else
  {
    fSameEventkTCentDist = nullptr;
    fMixedEventkTCentDist = nullptr;
  }
  if (fDomTBinning)
  {
    fSameEventmTDist = new TH2F *[nHists];
    fMixedEventmTDist = new TH2F *[nHists];
  }
  else
  {
    fSameEventmTDist = nullptr;
    fMixedEventmTDist = nullptr;
  }
  if (fmTMultPlots)
  {
    fSameEventmTMultDist = new TH2F **[nHists];
    fMixedEventmTMultDist = new TH2F **[nHists];
    fSameEventmTvsMultDist = new TH2F *[nHists];
    fMixedEventmTvsMultDist = new TH2F *[nHists];
  }
  else
  {
    fSameEventmTMultDist = nullptr;
    fMixedEventmTMultDist = nullptr;
    fSameEventmTvsMultDist = nullptr;
    fMixedEventmTvsMultDist = nullptr;
  }
  if (fPtQA)
  {
    fPtQADist = new TH2F *[nHists];
    fPtQADistSEPartOne = new TH2F *[nHists];
    fPtQADistSEPartTwo = new TH2F *[nHists];
    fPtQADistMEPartOne = new TH2F *[nHists];
    fPtQADistMEPartTwo = new TH2F *[nHists];

    fKstarPtQADistSEPartOne = new TH2F *[nHists];
    fKstarPtQADistSEPartTwo = new TH2F *[nHists];
    fKstarPtQADistMEPartOne = new TH2F *[nHists];
    fKstarPtQADistMEPartTwo = new TH2F *[nHists];
  }
  else
  {
    fPtQADist = nullptr;
    fPtQADistSEPartOne = nullptr;
    fPtQADistSEPartTwo = nullptr;
    fPtQADistMEPartOne = nullptr;
    fPtQADistMEPartTwo = nullptr;

    fKstarPtQADistSEPartOne = nullptr;
    fKstarPtQADistSEPartTwo = nullptr;
    fKstarPtQADistMEPartOne = nullptr;
    fKstarPtQADistMEPartTwo = nullptr;
  }
  if (fMassQA)
  {
    fMassQADistPart1 = new TH2F *[nHists];
    fMassQADistPart2 = new TH2F *[nHists];
    fPairInvMassQAD = new TH1F *[nHists];
    fPairInvMassKstarQAD = new TH2F *[nHists];

    fMEMassQADistPart1 = new TH2F *[nHists];
    fMEMassQADistPart2 = new TH2F *[nHists];
    fPairInvMEMassQAD = new TH1F *[nHists];
    fPairInvMEMassKstarQAD = new TH2F *[nHists];
  }
  else
  {
    fMassQADistPart1 = nullptr;
    fMassQADistPart2 = nullptr;
    fPairInvMassQAD = nullptr;
    fPairInvMassKstarQAD = nullptr;

    fMEMassQADistPart1 = nullptr;
    fMEMassQADistPart2 = nullptr;
    fPairInvMEMassQAD = nullptr;
    fPairInvMEMassKstarQAD = nullptr;
  }

  if (fdPhidEtaPlots)
  {
    if (!fmTDetaDPhi)
    {
      fdEtadPhiSE = new TH2F *[nHists];
      fdEtadPhiME = new TH2F *[nHists];
      fdEtadPhiSEmT = nullptr;
      fdEtadPhiMEmT = nullptr;
    }
    else
    {
      fdEtadPhiSEmT = new TH2F **[nHists];
      fdEtadPhiMEmT = new TH2F **[nHists];
      fdEtadPhiSE = nullptr;
      fdEtadPhiME = nullptr;
    }
  }
  else
  {
    fdEtadPhiSE = nullptr;
    fdEtadPhiME = nullptr;
    fdEtadPhiSEmT = nullptr;
    fdEtadPhiMEmT = nullptr;
  }
  if (fAncestors)
  {
    fSameEventDistCommon = new TH1F *[nHists];
    fSameEventDistNonCommon = new TH1F *[nHists];
    if (fdPhidEtaPlots)
    {
      fdEtadPhiSECommon = new TH2F *[nHists];
      fdEtadPhiSENonCommon = new TH2F *[nHists];
    }
    if (fDoMultBinning)
    {
      fSameEventMultDistCommon = new TH2F *[nHists];
      fSameEventMultDistNonCommon = new TH2F *[nHists];
    }
    if (fDomTBinning)
    {
      fSameEventmTDistCommon = new TH2F *[nHists];
      fSameEventmTDistNonCommon = new TH2F *[nHists];
    }
    if (fmTMultPlots)
    {
      fSameEventmTMultDistCommon = new TH2F **[nHists];
      fSameEventmTMultDistNonCommon = new TH2F **[nHists];
    }
  }
  else
  {
    fSameEventDistCommon = nullptr;
    fSameEventDistNonCommon = nullptr;
    fdEtadPhiSECommon = nullptr;
    fdEtadPhiSENonCommon = nullptr;
    fSameEventMultDistCommon = nullptr;
    fSameEventMultDistNonCommon = nullptr;
    fSameEventmTDistCommon = nullptr;
    fSameEventmTDistNonCommon = nullptr;
    fSameEventmTMultDistCommon = nullptr;
    fSameEventmTMultDistNonCommon = nullptr;
  }

  if (fpTOnepTTwokStarPlotsmT)
  {
    fSameEventpTOnepTTwokStar = new TH2F **[nHists];
    fMixedEventpTOnepTTwokStar = new TH2F **[nHists];
  }
  else
  {
    fSameEventpTOnepTTwokStar = nullptr;
    fMixedEventpTOnepTTwokStar = nullptr;
  }

  int Counter = 0;
  for (int iPar1 = 0; iPar1 < nParticles; ++iPar1)
  {
    for (int iPar2 = iPar1; iPar2 < nParticles; ++iPar2)
    {

      fPairs[Counter] = new TList();
      TString PairFolderName = TString::Format("Particle%d_Particle%d", iPar1, iPar2);
      fPairs[Counter]->SetName(PairFolderName.Data());
      fPairs[Counter]->SetOwner();
      fResults->Add(fPairs[Counter]);

      TString SameEventName = TString::Format("SEDist_Particle%d_Particle%d", iPar1,
                                              iPar2);
      fSameEventDist[Counter] = new TH1F(SameEventName.Data(),
                                         SameEventName.Data(), *itNBins,
                                         *itKMin, *itKMax);
      fPairs[Counter]->Add(fSameEventDist[Counter]);

      TString MixedEventName = TString::Format("MEDist_Particle%d_Particle%d", iPar1,
                                               iPar2);
      fMixedEventDist[Counter] = new TH1F(MixedEventName.Data(),
                                          MixedEventName.Data(), *itNBins,
                                          *itKMin, *itKMax);
      fPairs[Counter]->Add(fMixedEventDist[Counter]);
      // different multbins
      if (fDoMultBinning)
      {
        TString SameMultEventName = TString::Format("SEMultDist_Particle%d_Particle%d",
                                                    iPar1, iPar2);
        fSameEventMultDist[Counter] = new TH2F(SameMultEventName.Data(),
                                               SameMultEventName.Data(),
                                               *itNBins, *itKMin, *itKMax,
                                               multbins, 1, multbins + 1);
        fPairs[Counter]->Add(fSameEventMultDist[Counter]);

        TString MixedMultEventName = TString::Format("MEMultDist_Particle%d_Particle%d",
                                                     iPar1, iPar2);
        fMixedEventMultDist[Counter] = new TH2F(MixedMultEventName.Data(),
                                                MixedMultEventName.Data(),
                                                *itNBins, *itKMin, *itKMax,
                                                multbins, 1, multbins + 1);
        fPairs[Counter]->Add(fMixedEventMultDist[Counter]);
      }
      unsigned int DoThisPair = fWhichPairs.at(Counter);
      bool fillHists = DoThisPair > 0 ? true : false;
      if (fillHists && fDoCentBinning)
      {
        //        pp 13 TeV
        //        Mult (%)  Class η dNch/dη sys AN  Paper
        //        0 - 1 INEL>0  |η|<0.5 26.22 ±0.28 link  in preparation
        //        0 - 0.01      36.75 ±0.48
        //        0.01 - 0.05     32.43 ±0.37
        //        0.05 - 0.1      30.50 ±0.34
        //        0.1 - 1     25.59 ±0.27
        //        1 - 5     20.05 ±0.21
        //        5 - 10      16.18 ±0.17
        //        10 - 15     13.78 ±015
        //        15 - 20     12.02 ±0.13
        //        20 - 30     10.02 ±0.11
        //        30 - 40     7.93  ±0.08
        //        40 - 50     6.29  ±0.07
        //        50 - 70     4.45  ±0.05
        //        70 - 100      2.42  ±0.03
        Double_t centBins[14] = {0, 0.01, 0.05, 0.1, 1, 5, 10, 15, 20, 30, 40,
                                 50, 70, 100};
        TString SameCentEventName = TString::Format("SECentDist_Particle%d_Particle%d",
                                                    iPar1, iPar2);
        fSameEventCentDist[Counter] = new TH2F(SameCentEventName.Data(),
                                               SameCentEventName.Data(),
                                               *itNBins, *itKMin, *itKMax, 13,
                                               centBins);
        fPairs[Counter]->Add(fSameEventCentDist[Counter]);

        TString MixedCentEventName = TString::Format("MECentDist_Particle%d_Particle%d",
                                                     iPar1, iPar2);
        fMixedEventCentDist[Counter] = new TH2F(MixedCentEventName.Data(),
                                                MixedCentEventName.Data(),
                                                *itNBins, *itKMin, *itKMax, 13,
                                                centBins);
        fPairs[Counter]->Add(fMixedEventCentDist[Counter]);
      }
      // different relativeK Minv and kT bins
      if (fillHists && fDoMinvKtandRelativeKBinning)
      {
        TString SameEventMinvKtandRelativeKName = TString::Format("SERelKMinvDist_Particle%d_Particle%d",
                                                                  iPar1, iPar2);

        fSameEventMinvKtandRelativeKDist[Counter] = new TH2F(SameEventMinvKtandRelativeKName.Data(),
                                                             SameEventMinvKtandRelativeKName.Data(),
                                                             250, 0, 1,
                                                             1250, 0, 2.5);
        std::cout << "Object histogram is made for " << Counter << "  /n" << std::endl;

        fPairs[Counter]->Add(fSameEventMinvKtandRelativeKDist[Counter]);

        TString MixedEventMinvKtandRelativeKName = TString::Format("MERelKMinvDist_Particle%d_Particle%d",
                                                                   iPar1, iPar2);

        fMixedEventMinvKtandRelativeKDist[Counter] = new TH2F(MixedEventMinvKtandRelativeKName.Data(),
                                                              MixedEventMinvKtandRelativeKName.Data(),
                                                              250, 0, 1,
                                                              1250, 0, 2.5);
        fPairs[Counter]->Add(fMixedEventMinvKtandRelativeKDist[Counter]);
      }
      // kT Binning for different multbins
      if (fillHists && fDokTandMultBinning)
      {
        fSameEventkTandMultDist[Counter] = new TH2F *[multbins];
        fMixedEventkTandMultDist[Counter] = new TH2F *[multbins];
        if (fDokTandMultBinning)
        {
          fSameEventkTandMultPtDist[Counter] = new TH2F *[multbins];
        }
        for (int iMult = 1; iMult < multbins + 1; ++iMult)
        {
          TString SamekTandMultEventName = TString::Format(
              "SEkTandMultDist_Mult%i_Particle%d_Particle%d", iMult, iPar1,
              iPar2);
          fSameEventkTandMultDist[Counter][iMult] = new TH2F(
              SamekTandMultEventName.Data(), SamekTandMultEventName.Data(),
              *itNBins, *itKMin, *itKMax, *itNBins / 10, *itKMin,
              *itKMax * 1.5);
          fPairs[Counter]->Add(fSameEventkTandMultDist[Counter][iMult]);

          TString MixedkTandMultEventName = TString::Format(
              "MEkTandMultDist_Mult%i_Particle%d_Particle%d", iMult, iPar1,
              iPar2);
          fMixedEventkTandMultDist[Counter][iMult] = new TH2F(
              MixedkTandMultEventName.Data(), MixedkTandMultEventName.Data(),
              *itNBins, *itKMin, *itKMax, *itNBins / 10, *itKMin,
              *itKMax * 1.5);
          fPairs[Counter]->Add(fMixedEventkTandMultDist[Counter][iMult]);

          if (fDokTandMultPtBinning)
          {
            TString SamekTandMultPtEventName = TString::Format(
                "SEkTandMultPtDist_Mult%i_Particle%d_Particle%d", iMult, iPar1,
                iPar2);
            fSameEventkTandMultPtDist[Counter][iMult] = new TH2F(
                SamekTandMultPtEventName.Data(), SamekTandMultPtEventName.Data(),
                *itNBins, *itKMin, *itKMax, *itNBins / 10, *itKMin,
                *itKMax * 1.5);
            fPairs[Counter]->Add(fSameEventkTandMultPtDist[Counter][iMult]);
          }
        }
      }
      // kT Binning for different multbins for MC
      if (fillHists && fDokTandMultMCTrueBinning)
      {
        fSameEventkTandMultMCTrueDist[Counter] = new TH2F *[multbins];
        fMixedEventkTandMultMCTrueDist[Counter] = new TH2F *[multbins];
        for (int iMult = 1; iMult < multbins + 1; ++iMult)
        {
          TString SameEventkTandMultMCTrueName = TString::Format(
              "SEkTandMultMCTrueDist_Mult%i_Particle%d_Particle%d", iMult, iPar1,
              iPar2);
          fSameEventkTandMultMCTrueDist[Counter][iMult] = new TH2F(
              // Form("SEkTandMultMCTrueDist_Mult%i_Particle%d_Particle%d",iMult, iPar1, iPar2), Form("SEkTandMultMCTrueDist_Mult%i_Particle%d_Particle%d",iMult, iPar1, iPar2),
              SameEventkTandMultMCTrueName.Data(), SameEventkTandMultMCTrueName.Data(),
              *itNBins, *itKMin, *itKMax, *itNBins / 10, *itKMin,
              *itKMax * 1.5);
          fPairs[Counter]->Add(fSameEventkTandMultMCTrueDist[Counter][iMult]);
          TString MixedEventkTandMultMCTrueName = TString::Format(
              "MEkTandMultMCTrueDist_Mult%i_Particle%d_Particle%d", iMult, iPar1,
              iPar2);
          fMixedEventkTandMultMCTrueDist[Counter][iMult] = new TH2F(
              // Form("MEkTandMultMCTrueDist_Mult%i_Particle%d_Particle%d",iMult, iPar1, iPar2), Form("MEkTandMultMCTrueDist_Mult%i_Particle%d_Particle%d",iMult, iPar1, iPar2),
              MixedEventkTandMultMCTrueName.Data(), MixedEventkTandMultMCTrueName.Data(),
              *itNBins, *itKMin, *itKMax, *itNBins / 10, *itKMin,
              *itKMax * 1.5);
          fPairs[Counter]->Add(fMixedEventkTandMultMCTrueDist[Counter][iMult]);
        }
      }
      if (fillHists && fDokTBinning)
      {
        TString SamekTEventName = TString::Format("SEkTDist_Particle%d_Particle%d", iPar1,
                                                  iPar2);
        fSameEventkTDist[Counter] = new TH2F(SamekTEventName.Data(),
                                             SamekTEventName.Data(), *itNBins,
                                             *itKMin, *itKMax, *itNBins / 10,
                                             *itKMin, *itKMax * 1.5);
        fPairs[Counter]->Add(fSameEventkTDist[Counter]);

        TString MixedkTEventName = TString::Format("MEkTDist_Particle%d_Particle%d", iPar1,
                                                   iPar2);
        fMixedEventkTDist[Counter] = new TH2F(MixedkTEventName.Data(),
                                              MixedkTEventName.Data(), *itNBins,
                                              *itKMin, *itKMax, *itNBins / 10,
                                              *itKMin, *itKMax * 1.5);
        fPairs[Counter]->Add(fMixedEventkTDist[Counter]);
      }
      if (fillHists && fDokTCentralityBins)
      {
        const int nCentBins = fCentBins.size();
        fSameEventkTCentDist[Counter] = new TH2F *[nCentBins];
        fMixedEventkTCentDist[Counter] = new TH2F *[nCentBins];
        for (int iCent = 0; iCent < nCentBins; ++iCent)
        {
          TString SamekTCentEventName = TString::Format(
              "SEkTCentDist_Cent_%d_Particle%d_Particle%d", fCentBins[iCent],
              iPar1, iPar2);
          fSameEventkTCentDist[Counter][iCent] = new TH2F(
              SamekTCentEventName.Data(), SamekTCentEventName.Data(), *itNBins,
              *itKMin, *itKMax, *itNBins / 10, *itKMin, *itKMax * 1.5);
          fPairs[Counter]->Add(fSameEventkTCentDist[Counter][iCent]);

          TString MixedkTCentEventName = TString::Format(
              "MEkTCentDistCent%d_Particle%d_Particle%d", fCentBins[iCent],
              iPar1, iPar2);
          fMixedEventkTCentDist[Counter][iCent] = new TH2F(
              MixedkTCentEventName.Data(), MixedkTCentEventName.Data(),
              *itNBins, *itKMin, *itKMax, *itNBins / 10, *itKMin,
              *itKMax * 1.5);
          fPairs[Counter]->Add(fMixedEventkTCentDist[Counter][iCent]);
        }
      }
      if (fillHists && fmTMultPlots)
      {
        fSameEventmTMultDist[Counter] = new TH2F *[nmTBins];
        fMixedEventmTMultDist[Counter] = new TH2F *[nmTBins];
        for (unsigned int imT = 0; imT < nmTBins; ++imT)
        {
          TString SameEventmTMultName =
              TString::Format("SEmTMult_%d_Particle%d_Particle%d", imT, iPar1, iPar2);
          fSameEventmTMultDist[Counter][imT] = new TH2F(SameEventmTMultName.Data(),
                                                        SameEventmTMultName.Data(),
                                                        *itNBins, *itKMin, *itKMax,
                                                        multbins, 1, multbins + 1);
          fPairs[Counter]->Add(fSameEventmTMultDist[Counter][imT]);

          TString MixedEventmTMultName =
              TString::Format("MEmTMult_%d_Particle%d_Particle%d", imT, iPar1, iPar2);
          fMixedEventmTMultDist[Counter][imT] = new TH2F(MixedEventmTMultName.Data(),
                                                         MixedEventmTMultName.Data(),
                                                         *itNBins, *itKMin, *itKMax,
                                                         multbins, 1, multbins + 1);
          fPairs[Counter]->Add(fMixedEventmTMultDist[Counter][imT]);
        }

        TString SamemTvsMultEventName = TString::Format("SEmTvsMultDist_Particle%d_Particle%d", iPar1,
                                                        iPar2);
        fSameEventmTvsMultDist[Counter] = new TH2F(SamemTvsMultEventName.Data(),
                                                   SamemTvsMultEventName.Data(),
                                                   225, 0, 7.5, multbins, 1,
                                                   multbins + 1);
        fSameEventmTvsMultDist[Counter]->GetXaxis()->SetTitle("#it{k}* (GeV/#it{c})");
        fSameEventmTvsMultDist[Counter]->GetYaxis()->SetTitle("Mult. bin");
        fPairs[Counter]->Add(fSameEventmTvsMultDist[Counter]);

        TString MixedmTvsMultEventName = TString::Format("MEmTvsMultDist_Particle%d_Particle%d", iPar1,
                                                         iPar2);
        fMixedEventmTvsMultDist[Counter] = new TH2F(MixedmTvsMultEventName.Data(),
                                                    MixedmTvsMultEventName.Data(),
                                                    225, 0, 7.5, multbins, 1,
                                                    multbins + 1);
        fMixedEventmTvsMultDist[Counter]->GetXaxis()->SetTitle("#it{k}* (GeV/#it{c})");
        fMixedEventmTvsMultDist[Counter]->GetYaxis()->SetTitle("Mult. bin");
        fPairs[Counter]->Add(fMixedEventmTvsMultDist[Counter]);
      }

      if (fillHists && fDomTBinning)
      {
        TString SamemTEventName = TString::Format("SEmTDist_Particle%d_Particle%d", iPar1,
                                                  iPar2);
        fSameEventmTDist[Counter] = new TH2F(SamemTEventName.Data(),
                                             SamemTEventName.Data(), *itNBins,
                                             *itKMin, *itKMax, 225, 0, 7.5);
        fSameEventmTDist[Counter]->GetXaxis()->SetTitle("#it{k}* (GeV/#it{c})");
        fSameEventmTDist[Counter]->GetYaxis()->SetTitle("#it{m}_{T} (GeV/#it{c}^{2})");
        fPairs[Counter]->Add(fSameEventmTDist[Counter]);

        TString MixedmTEventName = TString::Format("MEmTDist_Particle%d_Particle%d", iPar1,
                                                   iPar2);
        fMixedEventmTDist[Counter] = new TH2F(MixedmTEventName.Data(),
                                              MixedmTEventName.Data(), *itNBins,
                                              *itKMin, *itKMax, 225, 0, 7.5);
        fMixedEventmTDist[Counter]->GetXaxis()->SetTitle("#it{k}* (GeV/#it{c})");
        fMixedEventmTDist[Counter]->GetYaxis()->SetTitle("#it{m}_{T} (GeV/#it{c}^{2})");
        fPairs[Counter]->Add(fMixedEventmTDist[Counter]);
      }

      if (fillHists && fdPhidEtaPlots)
      {
        if (!fmTDetaDPhi)
        {
          TString SameEventdPhidEtaName = TString::Format(
              "SEdPhidEtaDist_Particle%d_Particle%d", iPar1, iPar2);
          fdEtadPhiSE[Counter] = new TH2F(SameEventdPhidEtaName.Data(),
                                          SameEventdPhidEtaName.Data(), 80, -2.,
                                          2., 84, -2 * TMath::Pi() / 3,
                                          2 * TMath::Pi());
          fdEtadPhiSE[Counter]->GetXaxis()->SetTitle("#Delta#eta");
          fdEtadPhiSE[Counter]->GetYaxis()->SetTitle("#Delta#phi");
          fPairs[Counter]->Add(fdEtadPhiSE[Counter]);

          TString MixedEventdPhidEtaName = TString::Format(
              "MEdPhidEtaDist_Particle%d_Particle%d", iPar1, iPar2);
          fdEtadPhiME[Counter] = new TH2F(MixedEventdPhidEtaName.Data(),
                                          MixedEventdPhidEtaName.Data(), 80,
                                          -2., 2., 84, -2 * TMath::Pi() / 3,
                                          2 * TMath::Pi());
          fdEtadPhiME[Counter]->GetXaxis()->SetTitle("#Delta#eta");
          fdEtadPhiME[Counter]->GetYaxis()->SetTitle("#Delta#phi");
          fPairs[Counter]->Add(fdEtadPhiME[Counter]);
        }
        else
        {
          fdEtadPhiSEmT[Counter] = new TH2F *[nmTBins];
          fdEtadPhiMEmT[Counter] = new TH2F *[nmTBins];
          for (unsigned int imT = 0; imT < nmTBins; ++imT)
          {
            TString SameEventdPhidEtaName = TString::Format(
                "SEdPhidEtaDist_Particle%d_Particle%d_imT%.2f", iPar1, iPar2,
                fmTBins[imT]);
            fdEtadPhiSEmT[Counter][imT] = new TH2F(SameEventdPhidEtaName.Data(),
                                                   SameEventdPhidEtaName.Data(),
                                                   80, -2., 2., 84,
                                                   -2 * TMath::Pi() / 3,
                                                   2 * TMath::Pi());
            fdEtadPhiSEmT[Counter][imT]->GetXaxis()->SetTitle("#Delta#eta");
            fdEtadPhiSEmT[Counter][imT]->GetYaxis()->SetTitle("#Delta#phi");
            fPairs[Counter]->Add(fdEtadPhiSEmT[Counter][imT]);

            TString MixedEventdPhidEtaName = TString::Format(
                "MEdPhidEtaDist_Particle%d_Particle%d_imT%.2f", iPar1, iPar2,
                fmTBins[imT]);
            fdEtadPhiMEmT[Counter][imT] = new TH2F(
                MixedEventdPhidEtaName.Data(), MixedEventdPhidEtaName.Data(),
                80, -2., 2., 84, -2 * TMath::Pi() / 3, 2 * TMath::Pi());
            fdEtadPhiMEmT[Counter][imT]->GetXaxis()->SetTitle("#Delta#eta");
            fdEtadPhiMEmT[Counter][imT]->GetYaxis()->SetTitle("#Delta#phi");
            fPairs[Counter]->Add(fdEtadPhiMEmT[Counter][imT]);
          }
        }
      }
      // //For Common/Non Common Ancestors
      if (fillHists && fAncestors)
      {
        TString SameEventNameCommon = TString::Format("SEDistCommon_Particle%d_Particle%d", iPar1,
                                                      iPar2);
        fSameEventDistCommon[Counter] = new TH1F(SameEventNameCommon.Data(),
                                                 SameEventNameCommon.Data(), *itNBins,
                                                 *itKMin, *itKMax);
        fPairs[Counter]->Add(fSameEventDistCommon[Counter]);

        TString SameEventNameNonCommon = TString::Format("SEDistNonCommon_Particle%d_Particle%d", iPar1,
                                                         iPar2);
        fSameEventDistNonCommon[Counter] = new TH1F(SameEventNameNonCommon.Data(),
                                                    SameEventNameNonCommon.Data(), *itNBins,
                                                    *itKMin, *itKMax);
        fPairs[Counter]->Add(fSameEventDistNonCommon[Counter]);

        if (fDoMultBinning)
        {
          TString SameMultEventNameCommon = TString::Format("SEMultDistCommon_Particle%d_Particle%d",
                                                            iPar1, iPar2);
          fSameEventMultDistCommon[Counter] = new TH2F(SameMultEventNameCommon.Data(),
                                                       SameMultEventNameCommon.Data(),
                                                       *itNBins, *itKMin, *itKMax,
                                                       multbins, 1, multbins + 1);
          fPairs[Counter]->Add(fSameEventMultDistCommon[Counter]);
          TString SameMultEventNameNonCommon = TString::Format("SEMultDistNonCommon_Particle%d_Particle%d",
                                                               iPar1, iPar2);
          fSameEventMultDistNonCommon[Counter] = new TH2F(SameMultEventNameNonCommon.Data(),
                                                          SameMultEventNameNonCommon.Data(),
                                                          *itNBins, *itKMin, *itKMax,
                                                          multbins, 1, multbins + 1);
          fPairs[Counter]->Add(fSameEventMultDistNonCommon[Counter]);
        }
        if (fDomTBinning)
        {
          TString SamemTEventNameCommon = TString::Format("SEmTDistCommon_Particle%d_Particle%d", iPar1,
                                                          iPar2);
          fSameEventmTDistCommon[Counter] = new TH2F(SamemTEventNameCommon.Data(),
                                                     SamemTEventNameCommon.Data(), *itNBins,
                                                     *itKMin, *itKMax, 225, 0, 7.5);
          fPairs[Counter]->Add(fSameEventmTDistCommon[Counter]);

          TString SamemTEventNameNonCommon = TString::Format("SEmTDistNonCommon_Particle%d_Particle%d", iPar1,
                                                             iPar2);
          fSameEventmTDistNonCommon[Counter] = new TH2F(SamemTEventNameNonCommon.Data(),
                                                        SamemTEventNameNonCommon.Data(), *itNBins,
                                                        *itKMin, *itKMax, 225, 0, 7.5);
          fPairs[Counter]->Add(fSameEventmTDistNonCommon[Counter]);
        }

        if (fdPhidEtaPlots)
        {
          TString SameEventdPhidEtaNameCommon = TString::Format(
              "SEdPhidEtaDistCommon_Particle%d_Particle%d", iPar1, iPar2);
          fdEtadPhiSECommon[Counter] = new TH2F(SameEventdPhidEtaNameCommon.Data(),
                                                SameEventdPhidEtaNameCommon.Data(), 80, -2.,
                                                2., 84, -2 * TMath::Pi() / 3,
                                                2 * TMath::Pi());
          fdEtadPhiSECommon[Counter]->GetXaxis()->SetTitle("#Delta#eta");
          fdEtadPhiSECommon[Counter]->GetYaxis()->SetTitle("#Delta#phi");
          fPairs[Counter]->Add(fdEtadPhiSECommon[Counter]);

          TString SameEventdPhidEtaNameNonCommon = TString::Format(
              "SEdPhidEtaDistNonCommon_Particle%d_Particle%d", iPar1, iPar2);
          fdEtadPhiSENonCommon[Counter] = new TH2F(SameEventdPhidEtaNameNonCommon.Data(),
                                                   SameEventdPhidEtaNameNonCommon.Data(), 80, -2.,
                                                   2., 84, -2 * TMath::Pi() / 3,
                                                   2 * TMath::Pi());
          fdEtadPhiSENonCommon[Counter]->GetXaxis()->SetTitle("#Delta#eta");
          fdEtadPhiSENonCommon[Counter]->GetYaxis()->SetTitle("#Delta#phi");
          fPairs[Counter]->Add(fdEtadPhiSENonCommon[Counter]);
        }
        if (fmTMultPlots)
        {
          fSameEventmTMultDistCommon[Counter] = new TH2F *[nmTBins];
          fSameEventmTMultDistNonCommon[Counter] = new TH2F *[nmTBins];

          for (unsigned int imT = 0; imT < nmTBins; ++imT)
          {
            TString SameEventmTMultNameCommon =
                TString::Format("SEmTMultCommon_%d_Particle%d_Particle%d", imT, iPar1, iPar2);
            fSameEventmTMultDistCommon[Counter][imT] = new TH2F(SameEventmTMultNameCommon.Data(),
                                                                SameEventmTMultNameCommon.Data(),
                                                                *itNBins, *itKMin, *itKMax,
                                                                multbins, 1, multbins + 1);
            fPairs[Counter]->Add(fSameEventmTMultDistCommon[Counter][imT]);

            TString SameEventmTMultNameNonCommon =
                TString::Format("SEmTMultNonCommon_%d_Particle%d_Particle%d", imT, iPar1, iPar2);
            fSameEventmTMultDistNonCommon[Counter][imT] = new TH2F(SameEventmTMultNameNonCommon.Data(),
                                                                   SameEventmTMultNameNonCommon.Data(),
                                                                   *itNBins, *itKMin, *itKMax,
                                                                   multbins, 1, multbins + 1);
            fPairs[Counter]->Add(fSameEventmTMultDistNonCommon[Counter][imT]);
          }
        }
      }

      if (fillHists && fpTOnepTTwokStarPlotsmT)
      {

        fSameEventpTOnepTTwokStar[Counter] = new TH2F *[nmTBins];
        fMixedEventpTOnepTTwokStar[Counter] = new TH2F *[nmTBins];

        // Int_t bins[2] = {200,200};
        // Double_t xmin[2] =  {0.,0.};
        // Double_t xmax[2] = {5.0,5.0};

        for (unsigned int imT = 0; imT < nmTBins; ++imT)
        {
          TString SameEventName =
              TString::Format("SEmT_%d_pT_Particle%d_pT_Particle%d_kStarBelow%.1f", imT, iPar1, iPar2, fpTOnepTTwokStarCutOff);
          // fSameEventpTOnepTTwokStar[Counter][imT] = new THnSparseF(SameEventName.Data(),
          //                                              SameEventName.Data(),
          //                                               2,  bins, xmin, xmax);
          fSameEventpTOnepTTwokStar[Counter][imT] = new TH2F(SameEventName.Data(),
                                                             SameEventName.Data(),
                                                             200, 0., 5., 200, 0., 5.);
          fSameEventpTOnepTTwokStar[Counter][imT]->Sumw2();
          fSameEventpTOnepTTwokStar[Counter][imT]->GetXaxis()->SetTitle(Form("p_{T} Particle %d (GeV/c)", iPar1));
          fSameEventpTOnepTTwokStar[Counter][imT]->GetYaxis()->SetTitle(Form("p_{T} Particle %d (GeV/c)", iPar2));
          fPairs[Counter]->Add(fSameEventpTOnepTTwokStar[Counter][imT]);

          TString MixedEventName =
              TString::Format("MEmT_%d_pT_Particle%d_pT_Particle%d_kStar%.1f", imT, iPar1, iPar2, fpTOnepTTwokStarCutOff);
          // fMixedEventpTOnepTTwokStar[Counter][imT] = new THnSparseF(MixedEventName.Data(),
          //                                                MixedEventName.Data(),
          //                                                2,  bins, xmin, xmax);
          fMixedEventpTOnepTTwokStar[Counter][imT] = new TH2F(MixedEventName.Data(),
                                                              MixedEventName.Data(),
                                                              200, 0., 5., 200, 0., 5.);

          fMixedEventpTOnepTTwokStar[Counter][imT]->Sumw2();
          fMixedEventpTOnepTTwokStar[Counter][imT]->GetXaxis()->SetTitle(Form("p_{T} Particle %d (GeV/c)", iPar1));
          fMixedEventpTOnepTTwokStar[Counter][imT]->GetYaxis()->SetTitle(Form("p_{T} Particle %d (GeV/c)", iPar2));
          fPairs[Counter]->Add(fMixedEventpTOnepTTwokStar[Counter][imT]);
        }
      }

      if (!fMinimalBooking)
      {
        fPairQA[Counter] = new TList();
        TString PairQAName = TString::Format("QA_Particle%d_Particle%d", iPar1, iPar2);
        fPairQA[Counter]->SetName(PairQAName.Data());
        fPairQA[Counter]->SetOwner();
        fQA->Add(fPairQA[Counter]);
        TString PairCounterSEName = TString::Format("SEPairs_Particle%d_Particle%d", iPar1,
                                                    iPar2);
        fPairCounterSE[Counter] = new TH2F(PairCounterSEName.Data(),
                                           PairCounterSEName.Data(), 20, 0, 20,
                                           20, 0, 20);
        fPairCounterSE[Counter]->GetXaxis()->SetTitle(
            TString::Format("Particle%d", iPar1));
        fPairCounterSE[Counter]->GetYaxis()->SetTitle(
            TString::Format("Particle%d", iPar2));
        fPairQA[Counter]->Add(fPairCounterSE[Counter]);

        TString PairCounterMEName = TString::Format("MEPairs_Particle%d_Particle%d", iPar1,
                                                    iPar2);
        fPairCounterME[Counter] = new TH2F(PairCounterMEName.Data(),
                                           PairCounterMEName.Data(), 20, 0, 20,
                                           20, 0, 20);
        fPairCounterME[Counter]->GetXaxis()->SetTitle(
            TString::Format("Particle%d", iPar1));
        fPairCounterME[Counter]->GetYaxis()->SetTitle(
            TString::Format("Particle%d", iPar2));
        fPairQA[Counter]->Add(fPairCounterME[Counter]);

        TString EffMixingDepthName = TString::Format(
            "EffMixingDepth_Particle%d_Particle%d", iPar1, iPar2);
        int MixingDepth = conf->GetMixingDepth();
        ++MixingDepth; //+1 since counting starts at 0
        fEffMixingDepth[Counter] = new TH1F(EffMixingDepthName.Data(),
                                            EffMixingDepthName.Data(),
                                            MixingDepth, -0.5,
                                            MixingDepth - 0.5);
        fEffMixingDepth[Counter]->GetXaxis()->SetTitle("MixingDepth");
        fPairQA[Counter]->Add(fEffMixingDepth[Counter]);

        if (fillHists && fPtQA)
        {
          TString PtQAName = TString::Format("PtQA_Particle%d_Particle%d", iPar1, iPar2);
          fPtQADist[Counter] = new TH2F(PtQAName.Data(), PtQAName.Data(), 1000, 0,
                                        10, 1000, 0, 10);
          fPtQADist[Counter]->GetXaxis()->SetTitle(
              TString::Format("#it{p}_{T} Particle %d (GeV/#it{c})", iPar1));
          fPtQADist[Counter]->GetYaxis()->SetTitle(
              TString::Format("#it{p}_{T} Particle %d (GeV/#it{c})", iPar2));
          fPairQA[Counter]->Add(fPtQADist[Counter]);

          TString PtQASEPartOneName = TString::Format("PtSEPartOne_Particle%d_Particle%d",
                                                      iPar1, iPar2);
          fPtQADistSEPartOne[Counter] = new TH2F(PtQASEPartOneName.Data(),
                                                 PtQASEPartOneName.Data(), 375,
                                                 0, 7.5, multbins, 1,
                                                 multbins + 1);
          fPtQADistSEPartOne[Counter]->GetXaxis()->SetTitle(
              TString::Format("#it{p}_{T} Particle  %d (GeV/#it{c})", iPar1));
          fPtQADistSEPartOne[Counter]->GetYaxis()->SetTitle("Multiplicity");
          fPairQA[Counter]->Add(fPtQADistSEPartOne[Counter]);

          TString PtQASEPartTwoName = TString::Format("PtSEPartTwo_Particle%d_Particle%d",
                                                      iPar1, iPar2);
          fPtQADistSEPartTwo[Counter] = new TH2F(PtQASEPartTwoName.Data(),
                                                 PtQASEPartTwoName.Data(), 375,
                                                 0, 7.5, multbins, 1,
                                                 multbins + 1);
          fPtQADistSEPartTwo[Counter]->GetXaxis()->SetTitle(
              TString::Format("#it{p}_{T} Particle  %d (GeV/#it{c})", iPar2));
          fPtQADistSEPartTwo[Counter]->GetYaxis()->SetTitle("Multiplicity");
          fPairQA[Counter]->Add(fPtQADistSEPartTwo[Counter]);

          TString PtQAMEPartOneName = TString::Format("PtMEPartOne_Particle%d_Particle%d",
                                                      iPar1, iPar2);
          fPtQADistMEPartOne[Counter] = new TH2F(PtQAMEPartOneName.Data(),
                                                 PtQAMEPartOneName.Data(), 375,
                                                 0, 7.5, multbins, 1,
                                                 multbins + 1);
          fPtQADistMEPartOne[Counter]->GetXaxis()->SetTitle(
              TString::Format("#it{p}_{T} Particle  %d (GeV/#it{c})", iPar1));
          fPtQADistMEPartOne[Counter]->GetYaxis()->SetTitle("Multiplicity");
          fPairQA[Counter]->Add(fPtQADistMEPartOne[Counter]);

          TString PtQAMEPartTwoName = TString::Format("PtMEPartTwo_Particle%d_Particle%d",
                                                      iPar1, iPar2);
          fPtQADistMEPartTwo[Counter] = new TH2F(PtQAMEPartTwoName.Data(),
                                                 PtQAMEPartTwoName.Data(), 375,
                                                 0, 7.5, multbins, 1,
                                                 multbins + 1);
          fPtQADistMEPartTwo[Counter]->GetXaxis()->SetTitle(
              TString::Format("#it{p}_{T} Particle  %d (GeV/#it{c})", iPar2));
          fPtQADistMEPartTwo[Counter]->GetYaxis()->SetTitle("Multiplicity");
          fPairQA[Counter]->Add(fPtQADistMEPartTwo[Counter]);

          //------// the kstar vs pt histos //------//
          TString KstarPtQASEPartOneName = TString::Format("KstarPtSEPartOne_Particle%d_Particle%d",
                                                           iPar1, iPar2);
          fKstarPtQADistSEPartOne[Counter] = new TH2F(KstarPtQASEPartOneName.Data(),
                                                      KstarPtQASEPartOneName.Data(),
                                                      *itNBins, *itKMin, *itKMax,
                                                      375, 0, 7.5);
          fKstarPtQADistSEPartOne[Counter]->GetXaxis()->SetTitle("#it{k}* (GeV/#it{c})");
          fKstarPtQADistSEPartOne[Counter]->GetYaxis()->SetTitle(
              TString::Format("#it{p}_{T} Particle  %d (GeV/#it{c})", iPar1));
          fPairQA[Counter]->Add(fKstarPtQADistSEPartOne[Counter]);

          TString KstarPtQASEPartTwoName = TString::Format("KstarPtSEPartTwo_Particle%d_Particle%d",
                                                           iPar1, iPar2);
          fKstarPtQADistSEPartTwo[Counter] = new TH2F(KstarPtQASEPartTwoName.Data(),
                                                      KstarPtQASEPartTwoName.Data(),
                                                      *itNBins, *itKMin, *itKMax,
                                                      375, 0, 7.5);
          fKstarPtQADistSEPartTwo[Counter]->GetXaxis()->SetTitle("#it{k}* (GeV/#it{c})");
          fKstarPtQADistSEPartTwo[Counter]->GetYaxis()->SetTitle(
              TString::Format("#it{p}_{T} Particle  %d (GeV/#it{c})", iPar2));
          fPairQA[Counter]->Add(fKstarPtQADistSEPartTwo[Counter]);

          TString KstarPtQAMEPartOneName = TString::Format("KstarPtMEPartOne_Particle%d_Particle%d",
                                                           iPar1, iPar2);
          fKstarPtQADistMEPartOne[Counter] = new TH2F(KstarPtQAMEPartOneName.Data(),
                                                      KstarPtQAMEPartOneName.Data(),
                                                      *itNBins, *itKMin, *itKMax,
                                                      375, 0, 7.5);
          fKstarPtQADistMEPartOne[Counter]->GetXaxis()->SetTitle("#it{k}* (GeV/#it{c})");
          fKstarPtQADistMEPartOne[Counter]->GetYaxis()->SetTitle(
              TString::Format("#it{p}_{T} Particle  %d (GeV/#it{c})", iPar1));
          fPairQA[Counter]->Add(fKstarPtQADistMEPartOne[Counter]);

          TString KstarPtQAMEPartTwoName = TString::Format("KstarPtMEPartTwo_Particle%d_Particle%d",
                                                           iPar1, iPar2);
          fKstarPtQADistMEPartTwo[Counter] = new TH2F(KstarPtQAMEPartTwoName.Data(),
                                                      KstarPtQAMEPartTwoName.Data(),
                                                      *itNBins, *itKMin, *itKMax,
                                                      375, 0, 7.5);
          fKstarPtQADistMEPartTwo[Counter]->GetXaxis()->SetTitle("#it{k}* (GeV/#it{c})");
          fKstarPtQADistMEPartTwo[Counter]->GetYaxis()->SetTitle(
              TString::Format("#it{p}_{T} Particle  %d (GeV/#it{c})", iPar2));
          fPairQA[Counter]->Add(fKstarPtQADistMEPartTwo[Counter]);
        }

        if (fillHists && fMassQA)
        {
          TString MassQANamePart1 = TString::Format("MassQA_Particle%d_1", iPar1);
          TString MassQANamePart2 = TString::Format("MassQA_Particle%d_2", iPar2);
          TString MassQANamePart3 = TString::Format("InvMassQA_Particle%d_Particle%d",
                                                    iPar1, iPar2); //???????
          TString MassQANamePart4 = TString::Format("PDGInvMassKstarQA_Particle%d_Particle%d",
                                                    iPar1, iPar2); //???????

          const float massPart1 = TDatabasePDG::Instance()->GetParticle(
                                                              fPDGCode[iPar1])
                                      ->Mass();
          const float massPart2 = TDatabasePDG::Instance()->GetParticle(
                                                              fPDGCode[iPar2])
                                      ->Mass();

          fMassQADistPart1[Counter] = new TH2F(MassQANamePart1.Data(),
                                               MassQANamePart1.Data(), 512,
                                               massPart1 - 0.04,
                                               massPart1 + 0.04, *itNBins,
                                               *itKMin, *itKMax);
          fMassQADistPart1[Counter]->GetXaxis()->SetTitle(
              TString::Format("M_{Particle %d} (GeV/#it{c}^{2})", iPar1));
          fMassQADistPart1[Counter]->GetYaxis()->SetTitle(
              "#it{k}* (GeV/#it{c})");
          fPairQA[Counter]->Add(fMassQADistPart1[Counter]);

          fMassQADistPart2[Counter] = new TH2F(MassQANamePart2.Data(),
                                               MassQANamePart2.Data(), 512,
                                               massPart2 - 0.04,
                                               massPart2 + 0.04, *itNBins,
                                               *itKMin, *itKMax);
          fMassQADistPart2[Counter]->GetXaxis()->SetTitle(
              TString::Format("M_{Particle %d} (GeV/#it{c}^{2})", iPar2));
          fMassQADistPart2[Counter]->GetYaxis()->SetTitle(
              "#it{k}* (GeV/#it{c})");
          fPairQA[Counter]->Add(fMassQADistPart2[Counter]);

          fPairInvMassQAD[Counter] = new TH1F(MassQANamePart3.Data(),
                                              MassQANamePart3.Data(), 100,
                                              massPart1 + massPart2,
                                              4 * (massPart1 + massPart2));
          fPairInvMassQAD[Counter]->GetXaxis()->SetTitle(
              TString::Format("InvMass_{Particle %d_1}_{Particle %d_2} (GeV/#it{c}^{2})",
                              iPar1, iPar2));
          fPairQA[Counter]->Add(fPairInvMassQAD[Counter]);

          fPairInvMassKstarQAD[Counter] = new TH2F(MassQANamePart4.Data(),
                                                   MassQANamePart4.Data(), *itNBins * 2,
                                                   massPart1 + massPart2,
                                                   1.5 * (massPart1 + massPart2),
                                                   *itNBins / 2, *itKMin, *itKMax / 2.);
          fPairInvMassKstarQAD[Counter]->GetXaxis()->SetTitle(
              TString::Format("InvMass_{Particle %d_1}_{Particle %d_2} (GeV/#it{c}^{2})",
                              iPar1, iPar2));
          fPairInvMassKstarQAD[Counter]->GetYaxis()->SetTitle(
              "#it{k}* (GeV/#it{c})");
          fPairQA[Counter]->Add(fPairInvMassKstarQAD[Counter]);

          // MIXED EVENTS
          TString MEMassQANamePart1 = TString::Format("MEMassQA_Particle%d_1", iPar1);
          TString MEMassQANamePart2 = TString::Format("MEMassQA_Particle%d_2", iPar2);
          TString MEMassQANamePart3 = TString::Format("InvMEMassQA_Particle%d_Particle%d",
                                                      iPar1, iPar2); //???????
          TString MEMassQANamePart4 = TString::Format("PDGInvMEMassKstarQA_Particle%d_Particle%d",
                                                      iPar1, iPar2); //???????

          fMEMassQADistPart1[Counter] = new TH2F(MEMassQANamePart1.Data(),
                                                 MEMassQANamePart1.Data(), 512,
                                                 massPart1 - 0.04,
                                                 massPart1 + 0.04, *itNBins,
                                                 *itKMin, *itKMax);
          fMEMassQADistPart1[Counter]->GetXaxis()->SetTitle(
              TString::Format("M_{Particle %d} (GeV/#it{c}^{2})", iPar1));
          fMEMassQADistPart1[Counter]->GetYaxis()->SetTitle(
              "#it{k}* (GeV/#it{c})");
          fPairQA[Counter]->Add(fMEMassQADistPart1[Counter]);

          fMEMassQADistPart2[Counter] = new TH2F(MEMassQANamePart2.Data(),
                                                 MEMassQANamePart2.Data(), 512,
                                                 massPart2 - 0.04,
                                                 massPart2 + 0.04, *itNBins,
                                                 *itKMin, *itKMax);
          fMEMassQADistPart2[Counter]->GetXaxis()->SetTitle(
              TString::Format("M_{Particle %d} (GeV/#it{c}^{2})", iPar2));
          fMEMassQADistPart2[Counter]->GetYaxis()->SetTitle(
              "#it{k}* (GeV/#it{c})");
          fPairQA[Counter]->Add(fMEMassQADistPart2[Counter]);

          fPairInvMEMassQAD[Counter] = new TH1F(MEMassQANamePart3.Data(),
                                                MEMassQANamePart3.Data(), 100,
                                                massPart1 + massPart2,
                                                4 * (massPart1 + massPart2));
          fPairInvMEMassQAD[Counter]->GetXaxis()->SetTitle(
              TString::Format("InvMass_{Particle %d_1}_{Particle %d_2} (GeV/#it{c}^{2})",
                              iPar1, iPar2));
          fPairQA[Counter]->Add(fPairInvMEMassQAD[Counter]);

          fPairInvMEMassKstarQAD[Counter] = new TH2F(MEMassQANamePart4.Data(),
                                                     MEMassQANamePart4.Data(), *itNBins * 2,
                                                     massPart1 + massPart2,
                                                     1.5 * (massPart1 + massPart2),
                                                     *itNBins / 2, *itKMin, *itKMax / 2.);
          fPairInvMEMassKstarQAD[Counter]->GetXaxis()->SetTitle(
              TString::Format("InvMass_{Particle %d_1}_{Particle %d_2} (GeV/#it{c}^{2})",
                              iPar1, iPar2));
          fPairInvMEMassKstarQAD[Counter]->GetYaxis()->SetTitle(
              "#it{k}* (GeV/#it{c})");
          fPairQA[Counter]->Add(fPairInvMEMassKstarQAD[Counter]);
        }

        if (fillHists && fMomentumResolution)
        {
          TString MomResoSEName = TString::Format(
              "MomentumResolutionSE_Particle%d_Particle%d", iPar1, iPar2);
          fMomResolutionSE[Counter] = new TH2F(MomResoSEName.Data(),
                                               MomResoSEName.Data(), *itNBins, *itKMin, *itKMax,
                                               *itNBins, *itKMin, *itKMax);
          fMomResolutionSE[Counter]->GetXaxis()->SetTitle("k_{Generated}");
          fMomResolutionSE[Counter]->GetYaxis()->SetTitle("k_{Reco}");
          fPairQA[Counter]->Add(fMomResolutionSE[Counter]);

          TString MomResoSEAllName = TString::Format(
              "MomentumResolutionSEAll_Particle%d_Particle%d", iPar1, iPar2);
          fMomResolutionSEAll[Counter] = new TH2F(MomResoSEAllName.Data(),
                                                  MomResoSEAllName.Data(), *itNBins,
                                                  *itKMin, *itKMax, *itNBins, *itKMin, *itKMax);
          fMomResolutionSEAll[Counter]->GetXaxis()->SetTitle("k_{Generated}");
          fMomResolutionSEAll[Counter]->GetYaxis()->SetTitle("k_{Reco}");
          fPairQA[Counter]->Add(fMomResolutionSEAll[Counter]);

          TString MomResoMEName = TString::Format(
              "MomentumResolutionME_Particle%d_Particle%d", iPar1, iPar2);
          fMomResolutionME[Counter] = new TH2F(MomResoMEName.Data(),
                                               MomResoMEName.Data(), *itNBins, *itKMin, *itKMax,
                                               *itNBins, *itKMin, *itKMax);
          fMomResolutionME[Counter]->GetXaxis()->SetTitle("k_{Generated}");
          fMomResolutionME[Counter]->GetYaxis()->SetTitle("k_{Reco}");
          fPairQA[Counter]->Add(fMomResolutionME[Counter]);

          TString MomResoMEAllName = TString::Format(
              "MomentumResolutionMEAll_Particle%d_Particle%d", iPar1, iPar2);
          fMomResolutionMEAll[Counter] = new TH2F(MomResoMEAllName.Data(),
                                                  MomResoMEAllName.Data(), *itNBins, *itKMin, *itKMax,
                                                  *itNBins, *itKMin, *itKMax);
          fMomResolutionMEAll[Counter]->GetXaxis()->SetTitle("k_{Generated}");
          fMomResolutionMEAll[Counter]->GetYaxis()->SetTitle("k_{Reco}");
          fPairQA[Counter]->Add(fMomResolutionMEAll[Counter]);

          TString MomResoDistName = TString::Format(
              "MomentumResolutionDist_Particle%d_Particle%d", iPar1, iPar2);
          fMomResolutionDist[Counter] = new TH2F(MomResoDistName.Data(),
                                                 MomResoDistName.Data(), 500,
                                                 -0.3, 0.3, *itNBins, *itKMin, *itKMax);
          fMomResolutionDist[Counter]->GetXaxis()->SetTitle(
              "k_{Reco}-k_{Generated}");
          fMomResolutionDist[Counter]->GetYaxis()->SetTitle("k_{Generated}");
          fPairQA[Counter]->Add(fMomResolutionDist[Counter]);
        }
        if (fillHists && fPhiEtaPlots)
        {
          TList *RadList = new TList();
          RadList->SetOwner(true);
          RadList->SetName("PhiAtRad");
          fPairQA[Counter]->Add(RadList);
          const unsigned int nDaug1 = (unsigned int)DoThisPair / 10;
          if (nDaug1 > 9)
          {
            AliWarning(
                "you are doing something wrong, maximum of 3 Daughters supported \n");
          }
          const unsigned int nDaug2 = (unsigned int)DoThisPair % 10;
          const int nDaugComb = 81;
          fRadiiEtaPhiSE[Counter] = new TH2F **[nDaugComb]; // maximum of 81 combinations
          fRadiiEtaPhiME[Counter] = new TH2F **[nDaugComb];

          fIntRadiiQAEtaPhiSEBefore[Counter] = new TH2F *[nDaugComb]; // maximum of 81 combinations
          fIntRadiiQAEtaPhiMEBefore[Counter] = new TH2F *[nDaugComb];

          fIntRadiiQAEtaPhiSEAfter[Counter] = new TH2F *[nDaugComb]; // maximum of 81 combinations
          fIntRadiiQAEtaPhiMEAfter[Counter] = new TH2F *[nDaugComb];

          if (fPhiEtaPlotsSmallK)
          {
            fRadiiEtaPhiSEsmallK[Counter] = new TH2F **[nDaugComb];
            fRadiiEtaPhiMEsmallK[Counter] = new TH2F **[nDaugComb];
          }

          const int nRad = conf->GetNRadii();
          for (unsigned int iDaug1 = 0; iDaug1 < nDaug1; ++iDaug1)
          {
            for (unsigned int iDaug2 = 0; iDaug2 < nDaug2; ++iDaug2)
            {
              int DaugIndex = iDaug1 * 9 + iDaug2;
              fRadiiEtaPhiSE[Counter][DaugIndex] = new TH2F *[nRad];
              fRadiiEtaPhiME[Counter][DaugIndex] = new TH2F *[nRad];
              if (fPhiEtaPlotsSmallK)
              {
                fRadiiEtaPhiSEsmallK[Counter][DaugIndex] = new TH2F *[nRad];
                fRadiiEtaPhiMEsmallK[Counter][DaugIndex] = new TH2F *[nRad];
              }

              TString RadIntNameSE_Before = TString::Format(
                  "SERadQA_Before_Particle%d_Particle%d_DaugMix%d", iPar1,
                  iPar2, DaugIndex);
              TString RadIntNameME_Before = TString::Format(
                  "MERadQA_Before_Particle%d_Particle%d_DaugMix%d", iPar1,
                  iPar2, DaugIndex);
              fIntRadiiQAEtaPhiSEBefore[Counter][DaugIndex] = new TH2F(
                  RadIntNameSE_Before.Data(), RadIntNameSE_Before.Data(), 300,
                  -0.15, 0.15, 400, -0.2, 0.2);
              fIntRadiiQAEtaPhiSEBefore[Counter][DaugIndex]->GetXaxis()->SetTitle("#Delta#eta");
              fIntRadiiQAEtaPhiSEBefore[Counter][DaugIndex]->GetYaxis()->SetTitle("#Delta#phi");
              RadList->Add(
                  fIntRadiiQAEtaPhiSEBefore[Counter][DaugIndex]);
              fIntRadiiQAEtaPhiMEBefore[Counter][DaugIndex] = new TH2F(
                  RadIntNameME_Before.Data(), RadIntNameME_Before.Data(), 300,
                  -0.15, 0.15, 400, -0.2, 0.2);
              fIntRadiiQAEtaPhiMEBefore[Counter][DaugIndex]->GetXaxis()->SetTitle("#Delta#eta");
              fIntRadiiQAEtaPhiMEBefore[Counter][DaugIndex]->GetYaxis()->SetTitle("#Delta#phi");
              RadList->Add(
                  fIntRadiiQAEtaPhiMEBefore[Counter][DaugIndex]);

              TString RadIntNameSE_After = TString::Format(
                  "SERadQA_After_Particle%d_Particle%d_DaugMix%d", iPar1, iPar2,
                  DaugIndex);
              TString RadIntNameME_After = TString::Format(
                  "MERadQA_after_Particle%d_Particle%d_DaugMix%d", iPar1, iPar2,
                  DaugIndex);
              fIntRadiiQAEtaPhiSEAfter[Counter][DaugIndex] = new TH2F(
                  RadIntNameSE_After.Data(), RadIntNameSE_After.Data(), 300,
                  -0.15, 0.15, 400, -0.2, 0.2);
              fIntRadiiQAEtaPhiSEAfter[Counter][DaugIndex]->GetXaxis()->SetTitle(
                  "#Delta#eta");
              fIntRadiiQAEtaPhiSEAfter[Counter][DaugIndex]->GetYaxis()->SetTitle(
                  "#Delta#phi");
              RadList->Add(
                  fIntRadiiQAEtaPhiSEAfter[Counter][DaugIndex]);
              fIntRadiiQAEtaPhiMEAfter[Counter][DaugIndex] = new TH2F(
                  RadIntNameME_After.Data(), RadIntNameME_After.Data(), 300,
                  -0.15, 0.15, 400, -0.2, 0.2);
              fIntRadiiQAEtaPhiMEAfter[Counter][DaugIndex]->GetXaxis()->SetTitle(
                  "#Delta#eta");
              fIntRadiiQAEtaPhiMEAfter[Counter][DaugIndex]->GetYaxis()->SetTitle(
                  "#Delta#phi");
              RadList->Add(
                  fIntRadiiQAEtaPhiMEAfter[Counter][DaugIndex]);

              for (int iRad = 0; iRad < nRad; ++iRad)
              {
                TString RadNameSE = TString::Format(
                    "SERad_%i_Particle%d_Particle%d_DaugMix%d", iRad, iPar1,
                    iPar2, DaugIndex);
                TString RadNameME = TString::Format(
                    "MERad_%i_Particle%d_Particle%d_DaugMix%d", iRad, iPar1,
                    iPar2, DaugIndex);
                fRadiiEtaPhiSE[Counter][DaugIndex][iRad] = new TH2F(
                    RadNameSE.Data(), RadNameSE.Data(), 300, -0.15, 0.15, 400,
                    -0.2, 0.2);
                fRadiiEtaPhiSE[Counter][DaugIndex][iRad]->GetXaxis()->SetTitle(
                    "#Delta#eta");
                fRadiiEtaPhiSE[Counter][DaugIndex][iRad]->GetYaxis()->SetTitle(
                    "#Delta#phi");
                RadList->Add(fRadiiEtaPhiSE[Counter][DaugIndex][iRad]);
                fRadiiEtaPhiME[Counter][DaugIndex][iRad] = new TH2F(
                    RadNameME.Data(), RadNameME.Data(), 300, -0.15, 0.15, 400,
                    -0.2, 0.2);
                fRadiiEtaPhiME[Counter][DaugIndex][iRad]->GetXaxis()->SetTitle(
                    "#Delta#eta");
                fRadiiEtaPhiME[Counter][DaugIndex][iRad]->GetYaxis()->SetTitle(
                    "#Delta#phi");
                RadList->Add(fRadiiEtaPhiME[Counter][DaugIndex][iRad]);

                if (fPhiEtaPlotsSmallK)
                {
                  RadNameSE += "_smallK";
                  RadNameME += "_smallK";

                  fRadiiEtaPhiSEsmallK[Counter][DaugIndex][iRad] = new TH2F(
                      RadNameSE.Data(), RadNameSE.Data(), 300, -0.15, 0.15, 400,
                      -0.2, 0.2);
                  fRadiiEtaPhiSEsmallK[Counter][DaugIndex][iRad]->GetXaxis()->SetTitle("#Delta#eta");
                  fRadiiEtaPhiSEsmallK[Counter][DaugIndex][iRad]->GetYaxis()->SetTitle("#Delta#phi");
                  RadList->Add(
                      fRadiiEtaPhiSEsmallK[Counter][DaugIndex][iRad]);
                  fRadiiEtaPhiMEsmallK[Counter][DaugIndex][iRad] = new TH2F(
                      RadNameME.Data(), RadNameME.Data(), 300, -0.15, 0.15, 400,
                      -0.2, 0.2);
                  fRadiiEtaPhiMEsmallK[Counter][DaugIndex][iRad]->GetXaxis()->SetTitle("#Delta#eta");
                  fRadiiEtaPhiMEsmallK[Counter][DaugIndex][iRad]->GetYaxis()->SetTitle("#Delta#phi");
                  RadList->Add(
                      fRadiiEtaPhiMEsmallK[Counter][DaugIndex][iRad]);
                }
              }
            }
          }
        }
      }
      ++Counter;
      ++itNBins;
      ++itKMin;
      ++itKMax;
    }
  }
}
AliFemtoDreamCorrHists &AliFemtoDreamCorrHists::operator=(
    const AliFemtoDreamCorrHists &hists)
{
  if (this != &hists)
  {
    this->fQA = hists.fQA;
    this->fResults = hists.fResults;
    this->fPairs = hists.fPairs;
    this->fPairQA = hists.fPairQA;
    this->fMinimalBooking = hists.fMinimalBooking;
    this->fMomentumResolution = hists.fMomResolutionME;
    this->fPhiEtaPlots = hists.fPhiEtaPlots;
    this->fRelKThreshold = hists.fRelKThreshold;
    this->fSameEventDist = hists.fSameEventDist;
    this->fSameEventMinvKtandRelativeKDist = hists.fSameEventMinvKtandRelativeKDist;
    this->fSameEventMultDist = hists.fSameEventMultDist;
    this->fSameEventCentDist = hists.fSameEventCentDist;
    this->fSameEventmTDist = hists.fSameEventmTDist;
    this->fSameEventkTDist = hists.fSameEventkTDist;
    this->fSameEventmTvsMultDist = hists.fSameEventmTvsMultDist;
    this->fSameEventkTandMultDist = hists.fSameEventkTandMultDist;
    this->fSameEventkTandMultPtDist = hists.fSameEventkTandMultPtDist;
    this->fSameEventkTandMultMCTrueDist = hists.fSameEventkTandMultMCTrueDist;
    this->fSameEventkTCentDist = hists.fSameEventkTCentDist;
    this->fSameEventmTMultDist = hists.fSameEventmTMultDist;
    this->fPtQADist = hists.fPtQADist;
    this->fMassQADistPart1 = hists.fMassQADistPart1;
    this->fMassQADistPart2 = hists.fMassQADistPart2;
    this->fPairInvMassQAD = hists.fPairInvMassQAD;
    this->fPairInvMassKstarQAD = hists.fPairInvMassKstarQAD;
    this->fMEMassQADistPart1 = hists.fMEMassQADistPart1;
    this->fMEMassQADistPart2 = hists.fMEMassQADistPart2;
    this->fPairInvMEMassQAD = hists.fPairInvMEMassQAD;
    this->fPairInvMEMassKstarQAD = hists.fPairInvMEMassKstarQAD;
    this->fPairCounterSE = hists.fPairCounterSE;
    this->fMixedEventDist = hists.fMixedEventDist;
    this->fMixedEventMinvKtandRelativeKDist = hists.fMixedEventMinvKtandRelativeKDist;
    this->fMixedEventMultDist = hists.fMixedEventMultDist;
    this->fMixedEventCentDist = hists.fMixedEventCentDist;
    this->fMixedEventmTDist = hists.fMixedEventmTDist;
    this->fMixedEventkTDist = hists.fMixedEventkTDist;
    this->fMixedEventmTvsMultDist = hists.fMixedEventmTvsMultDist;
    this->fMixedEventkTandMultDist = hists.fMixedEventkTandMultDist;
    this->fMixedEventkTandMultMCTrueDist = hists.fMixedEventkTandMultMCTrueDist;
    this->fMixedEventkTCentDist = hists.fMixedEventkTCentDist;
    this->fMixedEventmTMultDist = hists.fMixedEventmTMultDist;
    this->fPairCounterME = hists.fPairCounterME;
    this->fMomResolutionME = hists.fMomResolutionME;
    this->fMomResolutionDist = hists.fMomResolutionDist;
    this->fRadiiEtaPhiSE = hists.fRadiiEtaPhiSE;
    this->fRadiiEtaPhiME = hists.fRadiiEtaPhiME;
    this->fIntRadiiQAEtaPhiSEAfter = hists.fIntRadiiQAEtaPhiSEAfter;
    this->fIntRadiiQAEtaPhiMEAfter = hists.fIntRadiiQAEtaPhiMEAfter;
    this->fRadiiEtaPhiSEsmallK = hists.fRadiiEtaPhiSEsmallK;
    this->fRadiiEtaPhiMEsmallK = hists.fRadiiEtaPhiMEsmallK;
    this->fdEtadPhiSE = hists.fdEtadPhiSE;
    this->fdEtadPhiME = hists.fdEtadPhiME;
    this->fdEtadPhiSEmT = hists.fdEtadPhiSEmT;
    this->fdEtadPhiMEmT = hists.fdEtadPhiMEmT;
    this->fEffMixingDepth = hists.fEffMixingDepth;
    this->fSameEventDistCommon = hists.fSameEventDistCommon;
    this->fSameEventDistNonCommon = hists.fSameEventDistNonCommon;
    this->fdEtadPhiSECommon = hists.fdEtadPhiSECommon;
    this->fdEtadPhiSENonCommon = hists.fdEtadPhiSENonCommon;
    this->fSameEventMultDistCommon = hists.fSameEventMultDistCommon;
    this->fSameEventMultDistNonCommon = hists.fSameEventMultDistNonCommon;
    this->fSameEventmTDistCommon = hists.fSameEventmTDistCommon;
    this->fSameEventmTDistNonCommon = hists.fSameEventmTDistNonCommon;
    this->fSameEventmTMultDistCommon = hists.fSameEventmTMultDistCommon;
    this->fSameEventmTMultDistNonCommon = hists.fSameEventmTMultDistNonCommon;
    this->fDoMultBinning = hists.fDoMultBinning;
    this->fDoCentBinning = hists.fDoCentBinning;
    this->fDokTBinning = hists.fDokTBinning;
    this->fDomTBinning = hists.fDomTBinning;
    this->fPtQA = hists.fPtQA;
    this->fDokTCentralityBins = hists.fDokTCentralityBins;
    this->fdPhidEtaPlots = hists.fdPhidEtaPlots;
    this->fCentBins = hists.fCentBins;
    this->fAncestors = hists.fAncestors;
    this->fRemoveAncestorsResonances = hists.fRemoveAncestorsResonances;
    this->fpTOnepTTwokStarPlotsmT = hists.fpTOnepTTwokStarPlotsmT;
    this->fpTOnepTTwokStarCutOff = hists.fpTOnepTTwokStarCutOff;
  }
  return *this;
}
AliFemtoDreamCorrHists::~AliFemtoDreamCorrHists()
{
  if (fPairs)
  {
    delete[] fPairs;
    delete fPairs;
  }
  if (fSameEventDist)
  {
    delete[] fSameEventDist;
    delete fSameEventDist;
  }
  if (fSameEventMinvKtandRelativeKDist)
  {
    delete[] fSameEventMinvKtandRelativeKDist;
    delete fSameEventMinvKtandRelativeKDist;
  }
  if (fSameEventMultDist)
  {
    delete[] fSameEventMultDist;
    delete fSameEventMultDist;
  }
  if (fSameEventmTDist)
  {
    delete[] fSameEventmTDist;
    delete fSameEventmTDist;
  }
  if (fSameEventmTvsMultDist)
  {
    delete[] fSameEventmTvsMultDist;
    delete fSameEventmTvsMultDist;
  }
  if (fPtQADist)
  {
    delete[] fPtQADist;
  }
  if (fMassQADistPart1)
  {
    delete[] fMassQADistPart1;
  }
  if (fMassQADistPart2)
  {
    delete[] fMassQADistPart2;
  }
  if (fPairInvMassQAD)
  {
    delete[] fPairInvMassQAD;
  }
  if (fPairInvMassKstarQAD)
  {
    delete[] fPairInvMassKstarQAD;
  }
  if (fMEMassQADistPart1)
  {
    delete[] fMEMassQADistPart1;
  }
  if (fMEMassQADistPart2)
  {
    delete[] fMEMassQADistPart2;
  }
  if (fPairInvMEMassQAD)
  {
    delete[] fPairInvMEMassQAD;
  }
  if (fPairInvMEMassKstarQAD)
  {
    delete[] fPairInvMEMassKstarQAD;
  }
  if (fSameEventkTandMultDist)
  {
    delete[] fSameEventkTandMultDist;
    delete fSameEventkTandMultDist;
  }
  if (fSameEventkTandMultPtDist)
  {
    delete[] fSameEventkTandMultPtDist;
    delete fSameEventkTandMultPtDist;
  }
  if (fSameEventkTandMultMCTrueDist)
  {
    delete[] fSameEventkTandMultMCTrueDist;
    delete fSameEventkTandMultMCTrueDist;
  }
  if (fSameEventkTDist)
  {
    delete[] fSameEventkTDist;
    delete fSameEventDist;
  }
  if (fSameEventkTCentDist)
  {
    delete[] fSameEventkTCentDist;
    delete fSameEventkTCentDist;
  }
  if (fMixedEventDist)
  {
    delete[] fMixedEventDist;
    delete fMixedEventDist;
  }
  if (fMixedEventMinvKtandRelativeKDist)
  {
    delete[] fMixedEventMinvKtandRelativeKDist;
    delete fMixedEventMinvKtandRelativeKDist;
  }
  if (fMixedEventMultDist)
  {
    delete[] fMixedEventMultDist;
    delete fMixedEventMultDist;
  }
  if (fMixedEventmTDist)
  {
    delete[] fMixedEventmTDist;
    delete fMixedEventmTDist;
  }
  if (fMixedEventkTDist)
  {
    delete[] fMixedEventkTDist;
    delete fMixedEventDist;
  }
  if (fMixedEventmTvsMultDist)
  {
    delete[] fMixedEventmTvsMultDist;
    delete fMixedEventmTvsMultDist;
  }
  if (fMixedEventkTandMultDist)
  {
    delete[] fMixedEventkTandMultDist;
    delete fMixedEventkTandMultDist;
  }
  if (fMixedEventkTandMultMCTrueDist)
  {
    delete[] fMixedEventkTandMultMCTrueDist;
    delete fMixedEventkTandMultMCTrueDist;
  }
  if (fMixedEventkTCentDist)
  {
    delete[] fMixedEventkTCentDist;
    delete fMixedEventkTCentDist;
  }
  if (fSameEventDistCommon)
  {
    delete[] fSameEventDistCommon;
    delete fSameEventDistCommon;
  }
  if (fSameEventDistNonCommon)
  {
    delete[] fSameEventDistNonCommon;
    delete fSameEventDistNonCommon;
  }
  if (fdEtadPhiSECommon)
  {
    delete[] fdEtadPhiSECommon;
    delete fdEtadPhiSECommon;
  }
  if (fdEtadPhiSENonCommon)
  {
    delete[] fdEtadPhiSENonCommon;
    delete fdEtadPhiSENonCommon;
  }
  if (fSameEventMultDistCommon)
  {
    delete[] fSameEventMultDistCommon;
    delete fSameEventMultDistCommon;
  }
  if (fSameEventMultDistNonCommon)
  {
    delete[] fSameEventMultDistNonCommon;
    delete fSameEventMultDistNonCommon;
  }
}

void AliFemtoDreamCorrHists::FillSameEventkTCentDist(int i, float kT,
                                                     float RelK, float cent)
{
  int centBin = -1;
  if (0 < cent)
  {
    for (unsigned int iCent = 0; iCent < fCentBins.size(); ++iCent)
    {
      if (cent < fCentBins[iCent])
      {
        centBin = iCent;
        break;
      }
    }
    if (centBin != -1)
    {
      fSameEventkTCentDist[i][centBin]->Fill(RelK, kT);
    }
  }
}

void AliFemtoDreamCorrHists::FillMixedEventkTCentDist(int i, float kT,
                                                      float RelK, float cent)
{
  int centBin = -1;
  if (0 < cent)
  {
    for (unsigned int iCent = 0; iCent < fCentBins.size(); ++iCent)
    {
      if (cent < fCentBins[iCent])
      {
        centBin = iCent;
        break;
      }
    }
    if (centBin != -1)
    {
      fMixedEventkTCentDist[i][centBin]->Fill(RelK, kT);
    }
  }
}

void AliFemtoDreamCorrHists::FilldPhidEtaSE(int iHist, float dPhi, float dEta,
                                            float mT)
{
  if (fmTDetaDPhi && fdPhidEtaPlots)
  {
    unsigned int pos = 0;
    for (auto it : fmTBins)
    {
      if (it > TMath::Abs(mT))
      {
        break;
      }
      else
      {
        pos++;
      }
    }
    if (pos >= fmTBins.size())
    {
      TString WarnMe = TString::Format("mT Bin for %.2f not found", mT);
      AliWarning(WarnMe.Data());
    }
    else
    {
      fdEtadPhiSEmT[iHist][pos]->Fill(dEta, dPhi);
    }
  }
  else if (fdPhidEtaPlots)
  {
    fdEtadPhiSE[iHist]->Fill(dEta, dPhi);
  }
}
void AliFemtoDreamCorrHists::FilldPhidEtaME(int iHist, float dPhi, float dEta,
                                            float mT)
{
  if (fmTDetaDPhi && fdPhidEtaPlots)
  {
    unsigned int pos = 0;
    for (auto it : fmTBins)
    {
      if (it > TMath::Abs(mT))
      {
        break;
      }
      else
      {
        pos++;
      }
    }
    if (pos >= fmTBins.size())
    {
      TString WarnMe = TString::Format("mT Bin for %.2f not found", mT);
      AliWarning(WarnMe.Data());
    }
    else
    {
      fdEtadPhiMEmT[iHist][pos]->Fill(dEta, dPhi);
    }
  }
  else if (fdPhidEtaPlots)
  {
    fdEtadPhiME[iHist]->Fill(dEta, dPhi);
  }
}

void AliFemtoDreamCorrHists::FilldPhidEtaSECommon(int iHist, float dPhi, float dEta,
                                                  float mT)
{
  if (fdPhidEtaPlots)
  {
    fdEtadPhiSECommon[iHist]->Fill(dEta, dPhi);
  }
}

void AliFemtoDreamCorrHists::FilldPhidEtaSENonCommon(int iHist, float dPhi, float dEta,
                                                     float mT)
{
  if (fdPhidEtaPlots)
  {
    fdEtadPhiSENonCommon[iHist]->Fill(dEta, dPhi);
  }
}

void AliFemtoDreamCorrHists::FillSameEventmTMultDist(int iHist, float mT, int iMult, float RelK)
{
  if (fmTMultPlots)
  {
    unsigned int pos = 0;
    for (auto it : fmTBins)
    {
      if (it > TMath::Abs(mT))
      {
        break;
      }
      else
      {
        pos++;
      }
    }
    if (pos >= fmTBins.size())
    {
      TString WarnMe = TString::Format("mT Bin for %.2f not found", mT);
      AliWarning(WarnMe.Data());
    }
    else
    {
      fSameEventmTMultDist[iHist][pos]->Fill(RelK, iMult);
      if (fSameEventmTvsMultDist[iHist] && RelK <= 0.2)
      {
        fSameEventmTvsMultDist[iHist]->Fill(mT, iMult);
      }
    }
  }
}
void AliFemtoDreamCorrHists::FillMixedEventmTMultDist(int iHist, float mT, int iMult, float RelK)
{
  if (fmTMultPlots)
  {
    unsigned int pos = 0;
    for (auto it : fmTBins)
    {
      if (it > TMath::Abs(mT))
      {
        break;
      }
      else
      {
        pos++;
      }
    }
    if (pos >= fmTBins.size())
    {
      TString WarnMe = TString::Format("mT Bin for %.2f not found", mT);
      AliWarning(WarnMe.Data());
    }
    else
    {
      fMixedEventmTMultDist[iHist][pos]->Fill(RelK, iMult);
      if (fMixedEventmTvsMultDist[iHist] && RelK <= 0.2)
      {
        fMixedEventmTvsMultDist[iHist]->Fill(mT, iMult);
      }
    }
  }
}

void AliFemtoDreamCorrHists::FillSameEventpTOnepTTwokStar(int iHist, float mT, float pTOne, float pTTwo, float RelK)
{
  if (fpTOnepTTwokStarPlotsmT)
  {
    unsigned int pos = 0;
    for (auto it : fmTBins)
    {
      if (it > TMath::Abs(mT))
      {
        break;
      }
      else
      {
        pos++;
      }
    }
    if (pos >= fmTBins.size())
    {
      TString WarnMe = TString::Format("mT Bin for %.2f not found", mT);
      AliWarning(WarnMe.Data());
    }
    else
    {
      // Double_t values[3] = {pTOne, pTTwo, RelK};
      // fSameEventpTOnepTTwokStar[iHist][pos]->Fill(values);
      if (RelK <= fpTOnepTTwokStarCutOff)
      {
        fSameEventpTOnepTTwokStar[iHist][pos]->Fill(pTOne, pTTwo);
      }
    }
  }
}

void AliFemtoDreamCorrHists::FillMixedEventpTOnepTTwokStar(int iHist, float mT, float pTOne, float pTTwo, float RelK)
{
  if (fpTOnepTTwokStarPlotsmT)
  {
    unsigned int pos = 0;
    for (auto it : fmTBins)
    {
      if (it > TMath::Abs(mT))
      {
        break;
      }
      else
      {
        pos++;
      }
    }
    if (pos >= fmTBins.size())
    {
      TString WarnMe = TString::Format("mT Bin for %.2f not found", mT);
      AliWarning(WarnMe.Data());
    }
    else
    {
      // Double_t values[3] = {pTOne, pTTwo, RelK};
      // fMixedEventpTOnepTTwokStar[iHist][pos]->Fill(values);
      if (RelK <= fpTOnepTTwokStarCutOff)
      {
        fMixedEventpTOnepTTwokStar[iHist][pos]->Fill(pTOne, pTTwo);
      }
    }
  }
}

void AliFemtoDreamCorrHists::FillSameEventmTMultDistCommon(int iHist, float mT, int iMult, float RelK)
{
  if (fmTMultPlots)
  {
    unsigned int pos = 0;
    for (auto it : fmTBins)
    {
      if (it > TMath::Abs(mT))
      {
        break;
      }
      else
      {
        pos++;
      }
    }
    if (pos >= fmTBins.size())
    {
      TString WarnMe = TString::Format("mT Bin for %.2f not found", mT);
      AliWarning(WarnMe.Data());
    }
    else
    {
      fSameEventmTMultDistCommon[iHist][pos]->Fill(RelK, iMult);
    }
  }
}

void AliFemtoDreamCorrHists::FillSameEventmTMultDistNonCommon(int iHist, float mT, int iMult, float RelK)
{
  if (fmTMultPlots)
  {
    unsigned int pos = 0;
    for (auto it : fmTBins)
    {
      if (it > TMath::Abs(mT))
      {
        break;
      }
      else
      {
        pos++;
      }
    }
    if (pos >= fmTBins.size())
    {
      TString WarnMe = TString::Format("mT Bin for %.2f not found", mT);
      AliWarning(WarnMe.Data());
    }
    else
    {
      fSameEventmTMultDistNonCommon[iHist][pos]->Fill(RelK, iMult);
    }
  }
}
