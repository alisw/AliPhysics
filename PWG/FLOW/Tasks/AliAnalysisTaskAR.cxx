/**
 * File              : AliAnalysisTaskAR.cxx
 * Author            : Anton Riedel <anton.riedel@tum.de>
 * Date              : 07.05.2021
 * Last Modified Date: 15.09.2022
 * Last Modified By  : Anton Riedel <anton.riedel@tum.de>
 */

/*************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

#include "AliAnalysisTaskAR.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliLog.h"
#include "AliMCEvent.h"
#include "AliMultSelection.h"
#include "AliVEvent.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include <TColor.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TList.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TString.h>
#include <TSystem.h>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <numeric>
#include <vector>

ClassImp(AliAnalysisTaskAR)

    AliAnalysisTaskAR::AliAnalysisTaskAR(const char *name)
    : AliAnalysisTaskSE(name),
      // Constructor
      // Base list for all output objects
      fHistList(nullptr), fHistListName("outputAnalysis"),
      // list holding all QA histograms
      fQAHistogramsList(nullptr), fQAHistogramsListName("QAHistograms"),
      fFillQAHistograms(kFALSE), fFillQACorHistogramsOnly(kFALSE),
      // sublist holding centrality estimator correlation QA histograms
      fCenCorQAHistogramsList(nullptr),
      fCenCorQAHistogramsListName("CenCorQAHistograms"),
      // sublist holding multiplicity correlation QA histograms
      fMulCorQAHistogramsList(nullptr),
      fMulCorQAHistogramsListName("MulCorQAHistograms"),
      // sublist holding centrality-multiplicity correlation QA histograms
      fCenMulCorQAHistogramsList(nullptr),
      fCenMulCorQAHistogramsListName("CenMulCorQAHistograms"),
      // sublist holding filterbit scan QA histograms
      fFBScanQAHistogramsList(nullptr),
      fFBScanQAHistogramsListName("FBScanQAHistograms"),
      // sublist holding self correlation QA histograms
      fSelfCorQAHistogramsList(nullptr),
      fSelfCorQAHistogramsListName("SelfCorQAHistograms"),
      // list holding all control histograms
      fControlHistogramsList(nullptr),
      fControlHistogramsListName("ControlHistograms"),
      // sublists holding all track control histograms
      fTrackControlHistogramsList(nullptr),
      fTrackControlHistogramsListName("TrackControlHistograms"),
      // sublists holding all event control histograms
      fEventControlHistogramsList(nullptr),
      fEventControlHistogramsListName("EventControlHistograms"),
      // cuts
      fTrackCutsCounterCumulativeName("TrackCutsCounterCumulative"),
      // fTrackCutsCounterCumulative(nullptr),
      fTrackCutsValuesName("TrackCutValues"), fTrackCutsValues(nullptr),
      fEventCutsCounterCumulativeName("EventCutsCounterCumulative"),
      // fEventCutsCounterCumulative(nullptr),
      fEventCutsValuesName("EventCutValues"), fEventCutsValues(nullptr),
      fFilterbit(128), fUseFilterbit(kFALSE), fChargedOnly(kFALSE),
      fPrimaryOnly(kFALSE), fMCPrimaryDef(kMCPhysicalPrim),
      fUseFakeTracks(kFALSE), fGlobalTracksOnly(kFALSE),
      fCentralityEstimator(kV0M), fUseCenCorCuts(kFALSE),
      fUseMulCorCuts(kFALSE), fUseCenFlatten(kFALSE), fCenFlattenHist(nullptr),
      // final results
      fFinalResultsList(nullptr), fFinalResultsListName("FinalResults"),
      fFinalResultHistogramsList(nullptr),
      fFinalResultHistogramsListName("FinalResultHistograms"),
      fFinalResultCorrelatorsList(nullptr),
      fFinalResultCorrelatorsListName("FinalResultCorrelator"),
      fFinalResultSymmetricCumulantsList(nullptr),
      fFinalResultSymmetricCumulantsListName("FinalResultSymmetricCumulant"),
      fFinalResultNormalizedSymmetricCumulantsList(nullptr),
      fFinalResultNormalizedSymmetricCumulantsListName(
          "FinalResultNormalizedSymmetricCumulant"),
      fFillControlHistogramsOnly(kFALSE), fUseNestedLoops(kFALSE),
      // flags for MC analysis
      fUseCustomSeed(kFALSE), fSeed(0), fMCOnTheFly(kFALSE), fMCClosure(kFALSE),
      fMCMultiplicity(nullptr), fLookUpTable({}), fUseFisherYates(kFALSE),
      fRandomizedTrackIndices({}), fUseFixedMultplicity(kFALSE),
      fFixedMultiplicy(12),
      // qvectors
      fUseWeights(kFALSE), fCorrelators({}), fSymmetricCumulants({}),
      fMapSCtoCor({}), fMapCorToIndex({}) {
  AliDebug(2, "AliAnalysisTaskAR::AliAnalysisTaskAR(const "
              "char *name");
  // create base list
  fHistList = new TList();
  fHistList->SetName(fHistListName);
  fHistList->SetOwner(kTRUE);
  // initialize all arrays
  InitializeArrays();
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskAR::AliAnalysisTaskAR()
    : AliAnalysisTaskSE(),
      // Dummy constructor
      // Base list for all output objects
      fHistList(nullptr), fHistListName("outputAnalysis"),
      // list holding all QA histograms
      fQAHistogramsList(nullptr), fQAHistogramsListName("QAHistograms"),
      fFillQAHistograms(kFALSE), fFillQACorHistogramsOnly(kFALSE),
      // sublist holding centrality estimator correlation QA histograms
      fCenCorQAHistogramsList(nullptr),
      fCenCorQAHistogramsListName("CenCorQAHistograms"),
      // sublist holding multiplicity correlation QA histograms
      fMulCorQAHistogramsList(nullptr),
      fMulCorQAHistogramsListName("MulCorQAHistograms"),
      // sublist holding centrality-multiplicity correlation QA histograms
      fCenMulCorQAHistogramsList(nullptr),
      fCenMulCorQAHistogramsListName("CenMulCorQAHistograms"),
      // sublist holding filterbit scan QA histograms
      fFBScanQAHistogramsList(nullptr),
      fFBScanQAHistogramsListName("FBScanQAHistograms"),
      // sublist holding self correlation QA histograms
      fSelfCorQAHistogramsList(nullptr),
      fSelfCorQAHistogramsListName("SelfCorQAHistograms"),
      // list holding all control histograms
      fControlHistogramsList(nullptr),
      fControlHistogramsListName("ControlHistograms"),
      // sublists holding all track control histograms
      fTrackControlHistogramsList(nullptr),
      fTrackControlHistogramsListName("TrackControlHistograms"),
      // sublists holding all event control histograms
      fEventControlHistogramsList(nullptr),
      fEventControlHistogramsListName("EventControlHistograms"),
      // cuts
      fTrackCutsCounterCumulativeName("TrackCutsCounterCumulative"),
      // fTrackCutsCounterCumulative(nullptr),
      fTrackCutsValuesName("TrackCutValues"), fTrackCutsValues(nullptr),
      fEventCutsCounterCumulativeName("EventCutsCounterCumulative"),
      // fEventCutsCounterCumulative(nullptr),
      fEventCutsValuesName("EventCutValues"), fEventCutsValues(nullptr),
      fFilterbit(128), fUseFilterbit(kFALSE), fChargedOnly(kFALSE),
      fPrimaryOnly(kFALSE), fMCPrimaryDef(kMCPhysicalPrim),
      fUseFakeTracks(kFALSE), fGlobalTracksOnly(kFALSE),
      fCentralityEstimator(kV0M), fUseCenCorCuts(kFALSE),
      fUseMulCorCuts(kFALSE), fUseCenFlatten(kFALSE), fCenFlattenHist(nullptr),
      // final results
      fFinalResultsList(nullptr), fFinalResultsListName("FinalResults"),
      fFinalResultHistogramsList(nullptr),
      fFinalResultHistogramsListName("FinalResultHistograms"),
      fFinalResultCorrelatorsList(nullptr),
      fFinalResultCorrelatorsListName("FinalResultCorrelator"),
      fFinalResultSymmetricCumulantsList(nullptr),
      fFinalResultSymmetricCumulantsListName("FinalResultSymmetricCumulant"),
      fFinalResultNormalizedSymmetricCumulantsList(nullptr),
      fFinalResultNormalizedSymmetricCumulantsListName(
          "FinalResultNormalizedSymmetricCumulant"),
      fFillControlHistogramsOnly(kFALSE), fUseNestedLoops(kFALSE),
      // flags for MC analysis
      fUseCustomSeed(kFALSE), fSeed(0), fMCOnTheFly(kFALSE), fMCClosure(kFALSE),
      fMCMultiplicity(nullptr), fLookUpTable({}), fUseFisherYates(kFALSE),
      fRandomizedTrackIndices({}), fUseFixedMultplicity(kFALSE),
      fFixedMultiplicy(12),
      // qvectors
      fUseWeights(kFALSE), fCorrelators({}), fSymmetricCumulants({}),
      fMapSCtoCor({}), fMapCorToIndex({}) {
  // initialize arrays
  InitializeArrays();
  AliDebug(2, "AliAnalysisTaskAR::AliAnalysisTaskAR()");
}

AliAnalysisTaskAR::~AliAnalysisTaskAR() {
  // Destructor

  // fHistlist owns (almost) all other data members, if we delete it, we will
  // recursively delete all other objects associative with this object
  if (fHistList) {
    delete fHistList;
  }

  // delete RNG
  if (gRandom) {
    delete gRandom;
  }
};

void AliAnalysisTaskAR::UserCreateOutputObjects() {
  // Called at every worker node to initialize objects

  // 1) Trick to avoid name clashes, part 1
  // 2) Book and nest all lists
  // 3) Book all objects
  // *) Trick to avoid name clashes, part 2
  Bool_t oldHistAddStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  // 2) Book and nest all lists
  this->BookAndNestAllLists();

  // 3) Book all objects
  if (fFillQAHistograms || fFillQACorHistogramsOnly) {
    this->BookQAHistograms();
  }
  this->BookControlHistograms();
  this->BookFinalResultHistograms();
  this->BookTrackBinHistograms();
  this->BookFinalResultSymmetricCumulants();
  this->BookFinalResultCorrelators();

  // seed RNG
  delete gRandom;
  fUseCustomSeed ? gRandom = new TRandom3(fSeed) : gRandom = new TRandom3(0);

  // *) Trick to avoid name clashes, part
  TH1::AddDirectory(oldHistAddStatus);

  PostData(1, fHistList);
}

void AliAnalysisTaskAR::Terminate(Option_t *) {
  // Accessing the merged output list for final compution or for off-line
  // computations (i.e. after merging)

  // fHistList = (TList *)GetOutputData(1);
  if (!fHistList) {
    std::cout << __LINE__ << ": Did not get " << fHistListName << std::endl;
    Fatal("Terminate", "Invalid Pointer to fHistList");
  }

  // get average azimuthal angle
  // should be pi, if distribution is flat, nice for trending
  fFinalResultHistograms[kAVGPHI]->SetBinContent(
      1, fTrackControlHistograms[kRECO][kPHI][kAFTER]->GetMean());
  // get average centrality
  // should be the average of the minimal and maximal centrality, if the
  // distribution is flat, nice for trending
  fFinalResultHistograms[kAVGCEN]->SetBinContent(
      1, fEventControlHistograms[kRECO][kCEN][kAFTER]->GetMean());
  // get minimum multiplicity
  // needed when running with fixed multiplicity and fischer-yates in later
  // analysis
  Int_t MinMulBin = 1;
  for (int i = 1;
       i < fEventControlHistograms[kRECO][kMULQ][kAFTER]->GetNbinsX(); i++) {
    if (fEventControlHistograms[kRECO][kMULQ][kAFTER]->GetBinContent(i) != 0) {
      MinMulBin = i;
      break;
    }
  }
  fFinalResultHistograms[kMINMUL]->SetBinContent(
      1,
      fEventControlHistograms[kRECO][kMULQ][kAFTER]->GetBinLowEdge(MinMulBin));
  // count number of events, needed for statistics in the end
  fFinalResultHistograms[kNUMBEROFEVENTS]->SetBinContent(
      1, fEventControlHistograms[kRECO][kMULQ][kAFTER]->GetEntries());
  // count number of tracks, needed for statistics in the end
  fFinalResultHistograms[kNUMBEROFTRACKS]->SetBinContent(
      1, fTrackControlHistograms[kRECO][kPHI][kAFTER]->GetEntries());

  // compute symmetric cumulant
  if (!fSymmetricCumulants.empty()) {
    FillSymmetricCumulant();
  }
}

void AliAnalysisTaskAR::InitializeArrays() {
  // Initialize all data members which are arrays in this method
  InitializeArraysForQAHistograms();
  InitializeArraysForTrackControlHistograms();
  InitializeArraysForEventControlHistograms();
  InitializeArraysForCuts();
  InitializeArraysForWeights();
  InitializeArraysForQvectors();
  InitializeArraysForFinalResultHistograms();
  // InitializeArraysForSymmetricCumulants();
  InitializeArraysForMCAnalysis();
}

void AliAnalysisTaskAR::InitializeArraysForQAHistograms() {
  // initialize array of QA histograms for the correlation between centrality
  // estimators
  // there are N(N-1)/2 such correlators, i.e. the number of elemets above/below
  // the diagonal of a square matrix
  for (int cen = 0; cen < LAST_ECENESTIMATORS * (LAST_ECENESTIMATORS - 1) / 2;
       ++cen) {
    for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
      fCenCorQAHistograms[cen][ba] = nullptr;
    }
  }

  // initialize names in a loop
  for (int cen1 = 0; cen1 < LAST_ECENESTIMATORS; ++cen1) {
    for (int cen2 = cen1 + 1; cen2 < LAST_ECENESTIMATORS; ++cen2) {
      for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
        fCenCorQAHistogramNames[IndexCorHistograms(
            cen1, cen2, LAST_ECENESTIMATORS)][ba][kNAME] =
            "fCorCenEstimatorQAHistograms[" + kCenEstimatorNames[cen1] + "+" +
            kCenEstimatorNames[cen2] + "]" + kBAName[ba];
        fCenCorQAHistogramNames[IndexCorHistograms(
            cen1, cen2, LAST_ECENESTIMATORS)][ba][kTITLE] =
            kCenEstimatorNames[cen1] + " vs " + kCenEstimatorNames[cen2] +
            kBAName[ba];
        fCenCorQAHistogramNames[IndexCorHistograms(
            cen1, cen2, LAST_ECENESTIMATORS)][ba][kXAXIS] =
            kCenEstimatorNames[cen1];
        fCenCorQAHistogramNames[IndexCorHistograms(
            cen1, cen2, LAST_ECENESTIMATORS)][ba][kYAXIS] =
            kCenEstimatorNames[cen2];
      }
    }
  }

  // set default bins
  Double_t CorCenEstimatorQAHistogramBins[2 * LAST_EBINS] = {
      // kBIN kLEDGE kUEDGE kBIN+LAST_EBINS kLEDGE+LAST_EBINS
      // kUEDGE+LAST_EBINS
      100., 0., 100., 100., 0., 100.};
  // initialize default bins
  for (int cen = 0; cen < LAST_ECENESTIMATORS * (LAST_ECENESTIMATORS - 1) / 2;
       ++cen) {
    for (int bin = 0; bin < 2 * LAST_EBINS; ++bin) {
      fCenCorQAHistogramBins[cen][bin] = CorCenEstimatorQAHistogramBins[bin];
    }
  }

  // initialize arrays for multiplicity correlation histograms
  // there are also N(N-1)/2 such correlators
  for (int mul = 0; mul < kMulEstimators * (kMulEstimators - 1) / 2; ++mul) {
    for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
      fMulCorQAHistograms[mul][ba] = nullptr;
    }
  }

  // initalize names in a loop
  for (int mul1 = 0; mul1 < kMulEstimators; ++mul1) {
    for (int mul2 = mul1 + 1; mul2 < kMulEstimators; ++mul2) {
      for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
        fMulCorQAHistogramNames[IndexCorHistograms(mul1, mul2, kMulEstimators)]
                               [ba][kNAME] =
                                   "fMulCorQAHistograms[" +
                                   kMulEstimatorNames[mul1] + "+" +
                                   kMulEstimatorNames[mul2] + "]" + kBAName[ba];
        fMulCorQAHistogramNames[IndexCorHistograms(mul1, mul2, kMulEstimators)]
                               [ba][kTITLE] =
                                   kMulEstimatorNames[mul1] + " vs " +
                                   kMulEstimatorNames[mul2] + kBAName[ba];
        fMulCorQAHistogramNames[IndexCorHistograms(mul1, mul2, kMulEstimators)]
                               [ba][kXAXIS] = kMulEstimatorNames[mul1];
        fMulCorQAHistogramNames[IndexCorHistograms(mul1, mul2, kMulEstimators)]
                               [ba][kYAXIS] = kMulEstimatorNames[mul2];
      }
    }
  }

  // set default bins
  Double_t MulCorQAHistogramBins[2 * LAST_EBINS] = {
      // kBIN kLEDGE kUEDGE kBIN+LAST_EBINS kLEDGE+LAST_EBINS
      // kUEDGE+LAST_EBINS
      1000., 0., 5000., 1000., 0., 5000.};
  // initialize default bins
  for (int mul = 0; mul < kMulEstimators * (kMulEstimators - 1) / 2; ++mul) {
    for (int bin = 0; bin < 2 * LAST_EBINS; ++bin) {
      fMulCorQAHistogramBins[mul][bin] = MulCorQAHistogramBins[bin];
    }
  }

  // initialize arrays for centrality-multiplicity correlation histograms
  for (int mul = 0; mul < kMulEstimators; mul++) {
    for (int ba = 0; ba < LAST_EBEFOREAFTER; ba++) {
      fCenMulCorQAHistograms[mul][ba] = nullptr;
    }
  }

  // initalize names in a loop
  for (int mul = 0; mul < kMulEstimators; ++mul) {
    for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
      fCenMulCorQAHistogramNames[mul][ba][kNAME] =
          "fCenMulCorQAHistogram[kCENESTIMATOR+" + kMulEstimatorNames[mul] +
          "]" + kBAName[ba];
      fCenMulCorQAHistogramNames[mul][ba][kTITLE] =
          "kCENESTIMATOR vs " + kMulEstimatorNames[mul] + kBAName[ba];
      fCenMulCorQAHistogramNames[mul][ba][kXAXIS] = "";
      fCenMulCorQAHistogramNames[mul][ba][kYAXIS] = kMulEstimatorNames[mul];
    }
  }

  // set default bins
  Double_t CenMulCorQAHistogramBins[2 * LAST_EBINS] = {
      // kBIN kLEDGE kUEDGE kBIN+LAST_EBINS kLEDGE+LAST_EBINS
      // kUEDGE+LAST_EBINS
      100., 0., 100., 1000., 0., 5000.};

  // initialize default bins
  for (int mul = 0; mul < kMulEstimators; ++mul) {
    for (int bin = 0; bin < 2 * LAST_EBINS; ++bin) {
      fCenMulCorQAHistogramBins[mul][bin] = CenMulCorQAHistogramBins[bin];
    }
  }

  // initialize array for filterbit scan QA histograms
  // i.e. count the filterbits associated with all tracks
  fFBScanQAHistogram = nullptr;
  // set name
  TString FBScanQAHistogramName[LAST_ENAME] = // NAME, TITLE, XAXIS, YAXIS
      {"fFBScanQAHistograms", "Filterbit Scan", "Filterbit", ""};
  // initialize names
  for (int name = 0; name < LAST_ENAME; ++name) {
    fFBScanQAHistogramName[name] = FBScanQAHistogramName[name];
  }

  // set bins
  Double_t FBScanQAHistogramBin[LAST_EBINS] = {kMaxFilterbit, 0.,
                                               kMaxFilterbit};
  // initialize bins
  for (int bin = 0; bin < LAST_EBINS; ++bin) {
    fFBScanQAHistogramBin[bin] = FBScanQAHistogramBin[bin];
  }

  // initialize array of track scan filterbit scan QA histograms
  // i.e. look at the spectra of kinematic varibles according to some
  // filterbit which filterbit are used are hardcode in header
  for (int track = 0; track < LAST_ETRACK; ++track) {
    for (int fb = 0; fb < kNumberofTestFilterBit; ++fb) {
      fFBTrackScanQAHistograms[track][fb] = nullptr;
    }
  }

  // set names
  TString FBTrackScanQAHistogramNames[LAST_ETRACK][LAST_ENAME] = {
      // NAME, TITLE, XAXIS, YAXIS
      {"fFBTrackScanQAHistogram[kPT]", "Filterbitscan p_{T}", "p_{t}", ""},
      {"fFBTrackScanQAHistogram[kETA]", "Filterbitscan #eta", "#eta", ""},
      {"fFBTrackScanQAHistogram[kPHI]", "Filterbitscan #varphi", "#varphi", ""},
      {"fFBTrackScanQAHistogram[kCHARGE]", "Filterbitscan Charge", "Q", ""},
      {"fFBTrackScanQAHistogram[kTPCNCLS]",
       "Filterbitscan number of TPC clusters", "", ""},
      {"fFBTrackScanQAHistogram[kTPCCROSSEDROWS]",
       "Filterbitscan number of rows crossed in TPC", "N_{TPCCROSSEDROWS}"},
      {"fFBTrackScanQAHistogram[kTPCNCLSFRACTIONSHARED]",
       "Filterbitscan number of clusters shared with TPC",
       "N_{TPCNCLSFRACTIONSHARED}"},
      {"fFBTrackScanQAHistogram[kTPCCHI2PERNDF]",
       "Filterbitscan #chi^{2}/NDF of TPC track", "N_{TPCCHI2PERNDF}"},
      {"fFBTrackScanQAHistogram[kITSNCLS]",
       "Filterbitscan number of ITS clusters", "", ""},
      {"fFBTrackScanQAHistogram[kCHI2PERNDF]", "Filterbitscan #chi^{2}/NDF", "",
       ""},
      {"fFBTrackScanQAHistogram[kDCAZ]", "Filterbitscan DCA", "", ""},
      {"fFBTrackScanQAHistogram[kDCAXY]", "Filterbitscan DCA in xy", "", ""},
  };

  // initialize names
  for (int track = 0; track < LAST_ETRACK; ++track) {
    for (int fb = 0; fb < kNumberofTestFilterBit; ++fb) {
      for (int name = 0; name < LAST_ENAME; ++name) {
        if (name == kNAME || name == kTITLE) {
          fFBTrackScanQAHistogramNames[track][fb][name] =
              FBTrackScanQAHistogramNames[track][name] +
              Form(" (%d) ", kTestFilterbit[fb]);
        } else
          fFBTrackScanQAHistogramNames[track][fb][name] =
              FBTrackScanQAHistogramNames[track][name];
      }
    }
  }
  // set default bins
  Double_t FBTrackScanHistogramBins[LAST_ETRACK][LAST_EBINS] = {
      // kBIN kLEDGE kUEDGE
      {1000., 0., 10.},           // kPT
      {1000., -2., 2.},           // kETA
      {360., 0., TMath::TwoPi()}, // kPHI
      {11., -5.5, 5.5},           // kCHARGE
      {200., 0., 200.},           // kTPCNCLS
      {200., 0., 200.},           // kTPCCROSSEDROWS
      {200., 0., 2.},             // kTPCNCLSFRACTIONSHARED
      {1000., 0., 10.},           // kTPCCHI2PERNDF
      {10., 0., 10.},             // kITSNCLS
      {1000., 0., 10.},           // kCHI2PERNDF
      {1000., -10., 10.},         // kDCAZ
      {1000, -10., 10.},          // kDCAXY
  };
  // initialize default bins
  for (int var = 0; var < LAST_ETRACK; ++var) {
    for (int bin = 0; bin < LAST_EBINS; ++bin) {
      fFBTrackScanQAHistogramBins[var][bin] =
          FBTrackScanHistogramBins[var][bin];
    }
  }

  // initialize arrays for self correlation QA histograms
  // i.e. compute k_i - k_j for i !=j for all kinematic variables
  // of all tracks in each event
  // if there are no self correlations, we expect a flat spectrum
  for (int var = 0; var < kKinematic; ++var) {
    for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
      fSelfCorQAHistograms[var][ba] = nullptr;
    }
  }

  // set names
  TString SelfCorQAHistogramNames[kKinematic][LAST_ENAME] = {
      // NAME, TITLE, XAXIS, YAXIS
      {"fSelfCorQAHistograms[kPT]", "p_{T}^{1}-p_{T}^{2}", "#Delta p_{T}", ""},
      {"fSelfCorQAHistograms[kETA]", "#eta_{1}-#eta_{2}", "#Delta #eta", ""},
      {"fSelfCorQAHistograms[kPHI]", "#varphi_{1}-#varphi_{2}",
       "#Delta #varphi", ""},
  };

  // initialize names
  for (int var = 0; var < kKinematic; ++var) {
    for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
      for (int name = 0; name < LAST_ENAME; ++name) {
        if (name == kNAME || name == kTITLE) {
          fSelfCorQAHistogramNames[var][ba][name] =
              SelfCorQAHistogramNames[var][name] + kBAName[ba];
        } else {
          fSelfCorQAHistogramNames[var][ba][name] =
              SelfCorQAHistogramNames[var][name];
        }
      }
    }
  }

  // set default bins
  Double_t SelfCorQAHistogramBins[kKinematic][LAST_EBINS] = {
      // kBIN kLEDGE kUEDGE
      {100., -0.1, 0.1}, // kPT
      {100., -0.1, 0.1}, // kPHI
      {100., -0.1, 0.1}, // kETA
  };
  // initialize default bins
  for (int var = 0; var < kKinematic; ++var) {
    for (int bin = 0; bin < LAST_EBINS; ++bin) {
      fSelfCorQAHistogramBins[var][bin] = SelfCorQAHistogramBins[var][bin];
    }
  }
}

void AliAnalysisTaskAR::InitializeArraysForTrackControlHistograms() {
  // initialize array of track control histograms
  for (int mode = 0; mode < LAST_EMODE; ++mode) {
    for (int var = 0; var < LAST_ETRACK; ++var) {
      for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
        fTrackControlHistograms[mode][var][ba] = nullptr;
      }
    }
  }

  // set names
  TString TrackControlHistogramNames[LAST_ETRACK][LAST_ENAME] = {
      // NAME, TITLE, XAXIS, YAXIS
      {"fTrackControlHistograms[kPT]", "p_{T}", "p_{T} [GeV]", ""},
      {"fTrackControlHistograms[kETA]", "#eta", "#eta", ""},
      {"fTrackControlHistograms[kPHI]", "#varphi", "#varphi", ""},
      {"fTrackControlHistograms[kCHARGE]", "Charge", "Q [e]", ""},
      {"fTrackControlHistograms[kTPCNCLS]", "Number of clusters in TPC",
       "N_{TPCNCLS}"},
      {"fTrackControlHistograms[kTPCCROSSEDROWS]",
       "Number of rows crossed in TPC", "N_{TPCCROSSEDROWS}"},
      {"fTrackControlHistograms[kTPCNCLSFRACTIONSHARED]",
       "Number of clusters shared with TPC", "N_{TPCNCLSFRACTIONSHARED}"},
      {"fTrackControlHistograms[kTPCCHI2PERNDF]", "#chi^{2}/NDF of TPC track",
       "N_{TPCCHI2PERNDF}"},
      {"fTrackControlHistograms[kITSNCLS]", "Number of clusters in ITS",
       "N_{ITSNCLS}"},
      {"fTrackControlHistograms[kCHI2PERNDF]", "#chi^{2}/NDF of track",
       "#chi^{2}/NDF", ""},
      {"fTrackControlHistograms[kDCAZ]", "DCA in Z", "DCA_{Z} [cm]"},
      {"fTrackControlHistograms[kDCAXY]", "DCA in XY", "DCA_{XY} [cm]"},
  };
  // initialize names
  for (int mode = 0; mode < LAST_EMODE; ++mode) {
    for (int var = 0; var < LAST_ETRACK; ++var) {
      for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
        for (int name = 0; name < LAST_ENAME; ++name) {
          if (name == kNAME || name == kTITLE) {
            fTrackControlHistogramNames[mode][var][ba][name] =
                kModeName[mode] + TrackControlHistogramNames[var][name] +
                kBAName[ba];
          } else {
            fTrackControlHistogramNames[mode][var][ba][name] =
                TrackControlHistogramNames[var][name];
          }
        }
      }
    }
  }

  // set default bins
  Double_t BinsTrackControlHistogramDefaults[LAST_EBINS] = {100., 0., 1000.};
  // initialize default bins
  for (int var = 0; var < LAST_ETRACK; ++var) {
    for (int bin = 0; bin < LAST_EBINS; ++bin) {
      fTrackControlHistogramBins[var][bin] =
          BinsTrackControlHistogramDefaults[bin];
    }
  }
  // initialize track cuts counter histograms
  for (int mode = 0; mode < LAST_EMODE; ++mode) {
    fTrackCutsCounter[mode] = nullptr;
  }
  // initialize name
  for (int mode = 0; mode < LAST_EMODE; ++mode) {
    fTrackCutsCounterNames[mode] = kModeName[mode] + "fTrackCutsCounter";
  }
  // set bin names of track cuts counter histogram
  TString TrackCutsCounterBinNames[LAST_ETRACK] = {
      "kPT",
      "kETA",
      "kPHI",
      "kCHARGE",
      "kTPCNCLS",
      "kTPCCROSSEDROWS",
      "kTPCNCLSFRACTIONSHARED",
      "kTPCCHI2PERNDF",
      "kITSNCLS",
      "kCHI2PERNDF",
      "kDCAZ",
      "kDCAXY",
  };
  // initialize bin names of track cuts counter histogram
  for (int name = 0; name < LAST_ETRACK; name++) {
    for (int mm = 0; mm < LAST_EMINMAX; ++mm) {
      fTrackCutsCounterBinNames[name][mm] =
          TrackCutsCounterBinNames[name] + kMMName[mm];
    }
  }
}

void AliAnalysisTaskAR::InitializeArraysForEventControlHistograms() {
  // initialize array of event control histograms
  for (int mode = 0; mode < LAST_EMODE; ++mode) {
    for (int var = 0; var < LAST_EEVENT; ++var) {
      for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
        fEventControlHistograms[mode][var][ba] = nullptr;
      }
    }
  }

  // set name
  TString EventControlHistogramNames[LAST_EEVENT][LAST_ENAME] = {
      // NAME, TITLE, XAXIS, YAXIS
      {"fEventControlHistograms[kMUL]", "Multiplicity (without track cuts)",
       "M", ""},
      {"fEventControlHistograms[kMULQ]", "Multiplicity (with track cuts)", "M",
       ""},
      {"fEventControlHistograms[kMULW]", "Multiplicity (computed from weights)",
       "M", ""},
      {"fEventControlHistograms[kMULREF]", "Reference Multipliticy", "M_{ref}",
       ""},
      {"fEventControlHistograms[kNCONTRIB]", "Number of Contributors",
       "M_{contrib}", ""},
      {"fEventControlHistograms[kCEN]", "Centrality", "Centrality Percentile",
       ""},
      {"fEventControlHistograms[kX]", "Primary Vertex X", "#bf{V}_{X} [cm]",
       ""},
      {"fEventControlHistograms[kY]", "Primary Vertex Y", "#bf{V}_{Y} [cm]",
       ""},
      {"fEventControlHistograms[kZ]", "Primary Vertex Z", "#bf{V}_{Z} [cm]",
       ""},
      {"fEventControlHistograms[kVPOS]", "Vertex Position", "|#bf{V}| [cm]",
       ""},
  };

  // initialize names
  for (int mode = 0; mode < LAST_EMODE; ++mode) {
    for (int var = 0; var < LAST_EEVENT; ++var) {
      for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
        for (int name = 0; name < LAST_ENAME; ++name) {
          if (name == kNAME || name == kTITLE) {
            fEventControlHistogramNames[mode][var][ba][name] =
                kModeName[mode] + EventControlHistogramNames[var][name] +
                kBAName[ba];
          } else
            fEventControlHistogramNames[mode][var][ba][name] =
                EventControlHistogramNames[var][name];
        }
      }
    }
  }

  // set default bins
  Double_t BinsEventControlHistogramDefaults[LAST_EBINS] = {100., 0., 10000.};
  // initialize default bins
  for (int var = 0; var < LAST_EEVENT; ++var) {
    for (int bin = 0; bin < LAST_EBINS; ++bin) {
      fEventControlHistogramBins[var][bin] =
          BinsEventControlHistogramDefaults[bin];
    }
  }
  // initialize event cuts counter histogram
  for (int mode = 0; mode < LAST_EMODE; ++mode) {
    fEventCutsCounter[mode] = nullptr;
  }
  // initialize names
  for (int mode = 0; mode < LAST_EMODE; ++mode) {
    fEventCutsCounterNames[mode] = kModeName[mode] + "fEventCutsCounter";
  }
  // initialize bin names of event cuts counter histogram
  TString EventCutsCounterBinNames[LAST_EEVENT] = {
      "kMUL", "kMULQ", "kMULW", "kMULREF", "kNCONTRIB",
      "kCEN", "kX",    "kY",    "kZ",      "kVPOS",
  };
  for (int name = 0; name < LAST_EEVENT; name++) {
    for (int mm = 0; mm < LAST_EMINMAX; ++mm) {
      fEventCutsCounterBinNames[name][mm] =
          EventCutsCounterBinNames[name] + kMMName[mm];
    }
  }
}

void AliAnalysisTaskAR::InitializeArraysForCuts() {
  // initialize all arrays for cuts

  // initialize array for track cuts
  for (int var = 0; var < LAST_ETRACK; ++var) {
    fUseTrackCuts[var] = kFALSE;
    for (int mm = 0; mm < LAST_EMINMAX; ++mm) {
      fTrackCuts[var][mm] = -99.;
    }
  }

  // initialize array for event cuts
  for (int var = 0; var < LAST_EEVENT; ++var) {
    fUseEventCuts[var] = kFALSE;
    for (int mm = 0; mm < LAST_EMINMAX; ++mm) {
      fEventCuts[var][mm] = -999;
    }
  }

  for (int i = 0; i < 2; ++i) {
    fCenCorCut[i] = -999.;
  }
  // initialize array for centrality estimators
  for (int cen = 0; cen < LAST_ECENESTIMATORS; ++cen) {
    fCentrality[cen] = 0;
  }

  for (int i = 0; i < 2; ++i) {
    fMulCorCut[i] = -999.;
  }
  // initialize array for multiplicity estimators
  for (int mul = 0; mul < kMulEstimators; ++mul) {
    fMultiplicity[mul] = 0;
  }
}

void AliAnalysisTaskAR::InitializeArraysForFinalResultHistograms() {
  // initialize array for final result histograms
  for (int var = 0; var < LAST_EFINALHIST; ++var) {
    fFinalResultHistograms[var] = nullptr;
  }

  // set name
  TString FinalResultHistogramNames[LAST_EFINALHIST][LAST_ENAME] = {
      // NAME, TITLE, XAXIS
      {"fFinalResultHistograms[kAVGPHI]", "Average #varphi", "#varphi",
       ""}, // kAVGPHI
      {"fFinalResultHistograms[kAVGCEN]", "Average centrality",
       "Centrality Percentile", ""}, // kAVGCEN
      {"fFinalResultHistograms[kMINMUL]", "Minimal multiplicity [kMULQ]", "M",
       ""}, // kMINMUL
      {"fFinalResultHistograms[kNUMBEROFEVENTS]",
       "Number of surviving of events", "", ""}, // kNUMBEROFEVENTS
      {"fFinalResultHistograms[kNUMBEROFTRACKS]",
       "Number of surviving of tracks", "", ""}, // kNUMBEROFTRACKS
  };

  // initialize name
  for (int var = 0; var < LAST_EFINALHIST; ++var) {
    for (int name = 0; name < LAST_ENAME; ++name) {
      fFinalResultHistogramNames[var][name] =
          FinalResultHistogramNames[var][name];
    }
  }

  // default bins
  Double_t BinsFinalResultHistogramDefaults[LAST_EBINS] = {// BIN LEDGE UEDGE
                                                           1., 0., 1.};
  // initialize default bins
  for (int var = 0; var < LAST_EFINALHIST; ++var) {
    for (int bin = 0; bin < LAST_EBINS; ++bin) {
      fFinalResultHistogramBins[var][bin] =
          BinsFinalResultHistogramDefaults[bin];
    }
  }
}

void AliAnalysisTaskAR::InitializeArraysForWeights() {
  // initialize all necessary components for weight and acceptance histograms
  for (int k = 0; k < kKinematic; ++k) {
    fWeightHistogram[k] = nullptr;
  }
}

void AliAnalysisTaskAR::InitializeArraysForQvectors() {
  // initalize arrays for Q-vectors
  for (Int_t h = 0; h < kMaxHarmonic; h++) {
    for (Int_t p = 0; p < kMaxPower; p++) {
      fQvector[h][p] = TComplex(0., 0.);
    }
  }
  for (int i = 0; i < kKinematic - 1; i++) {
    fTrackBinsHistogram[i] = nullptr;
  }

  for (int i = 0; i < kKinematic; i++) {
    fKinematics[i] = {};
    fKinematicWeights[i] = {};
  }
}

void AliAnalysisTaskAR::InitializeArraysForMCAnalysis() {
  // initialize arrays for MC analysis
  for (int var = 0; var < kKinematic; var++) {
    fMCKinematicVariables[var] = 0;
    fMCKinematicPDFs[var] = nullptr;
    fAcceptanceHistogram[var] = nullptr;
  }
}

void AliAnalysisTaskAR::BookAndNestAllLists() {
  // Book and nest all lists nested in the base list fHistList

  // 1. Book and nest list for QA histograms
  // 2. Book and nest list for control histograms
  // 3. Book and nest list for final results

  if (!fHistList) {
    std::cout << __LINE__ << ": Did not get " << fHistListName << std::endl;
    Fatal("BookAndNestAllLists", "Invalid Pointer");
  }
  // 1. Book and nest lists for QA histograms
  if (fFillQAHistograms || fFillQACorHistogramsOnly) {
    fQAHistogramsList = new TList();
    fQAHistogramsList->SetName(fQAHistogramsListName);
    fQAHistogramsList->SetOwner(kTRUE);
    fHistList->Add(fQAHistogramsList);

    // centrality correlation QA histograms
    fCenCorQAHistogramsList = new TList();
    fCenCorQAHistogramsList->SetName(fCenCorQAHistogramsListName);
    fCenCorQAHistogramsList->SetOwner(kTRUE);
    fQAHistogramsList->Add(fCenCorQAHistogramsList);

    // multiplicity correlation QA histograms
    fMulCorQAHistogramsList = new TList();
    fMulCorQAHistogramsList->SetName(fMulCorQAHistogramsListName);
    fMulCorQAHistogramsList->SetOwner(kTRUE);
    fQAHistogramsList->Add(fMulCorQAHistogramsList);

    // centrality-multiplicity correlation QA histograms
    fCenMulCorQAHistogramsList = new TList();
    fCenMulCorQAHistogramsList->SetName(fCenMulCorQAHistogramsListName);
    fCenMulCorQAHistogramsList->SetOwner(kTRUE);
    fQAHistogramsList->Add(fCenMulCorQAHistogramsList);

    if (!fFillQACorHistogramsOnly) {
      // filterbit QA histograms
      fFBScanQAHistogramsList = new TList();
      fFBScanQAHistogramsList->SetName(fFBScanQAHistogramsListName);
      fFBScanQAHistogramsList->SetOwner(kTRUE);
      fQAHistogramsList->Add(fFBScanQAHistogramsList);

      // self correlation QA histograms
      fSelfCorQAHistogramsList = new TList();
      fSelfCorQAHistogramsList->SetName(fSelfCorQAHistogramsListName);
      fSelfCorQAHistogramsList->SetOwner(kTRUE);
      fQAHistogramsList->Add(fSelfCorQAHistogramsList);
    }
  }

  // 2. Book and nest lists for control histograms
  fControlHistogramsList = new TList();
  fControlHistogramsList->SetName(fControlHistogramsListName);
  fControlHistogramsList->SetOwner(kTRUE);
  fHistList->Add(fControlHistogramsList);

  // track control histograms
  fTrackControlHistogramsList = new TList();
  fTrackControlHistogramsList->SetName(fTrackControlHistogramsListName);
  fTrackControlHistogramsList->SetOwner(kTRUE);
  fControlHistogramsList->Add(fTrackControlHistogramsList);

  // event control histograms
  fEventControlHistogramsList = new TList();
  fEventControlHistogramsList->SetName(fEventControlHistogramsListName);
  fEventControlHistogramsList->SetOwner(kTRUE);
  fControlHistogramsList->Add(fEventControlHistogramsList);

  // 3. Book and nest lists for final results
  fFinalResultsList = new TList();
  fFinalResultsList->SetName(fFinalResultsListName);
  fFinalResultsList->SetOwner(kTRUE);
  fHistList->Add(fFinalResultsList);

  // final result histograms
  fFinalResultHistogramsList = new TList();
  fFinalResultHistogramsList->SetName(fFinalResultHistogramsListName);
  fFinalResultHistogramsList->SetOwner(kTRUE);
  fFinalResultsList->Add(fFinalResultHistogramsList);

  // final result profiles of correlators
  fFinalResultCorrelatorsList = new TList();
  fFinalResultCorrelatorsList->SetName(fFinalResultCorrelatorsListName);
  fFinalResultCorrelatorsList->SetOwner(kTRUE);
  fFinalResultsList->Add(fFinalResultCorrelatorsList);

  // final result histograms of symmetric cummulants
  fFinalResultSymmetricCumulantsList = new TList();
  fFinalResultSymmetricCumulantsList->SetName(
      fFinalResultSymmetricCumulantsListName);
  fFinalResultSymmetricCumulantsList->SetOwner(kTRUE);
  fFinalResultsList->Add(fFinalResultSymmetricCumulantsList);

  // final result histograms for normalized symmetric cumulants
  fFinalResultNormalizedSymmetricCumulantsList = new TList();
  fFinalResultNormalizedSymmetricCumulantsList->SetName(
      fFinalResultNormalizedSymmetricCumulantsListName);
  fFinalResultNormalizedSymmetricCumulantsList->SetOwner(kTRUE);
  fFinalResultsList->Add(fFinalResultNormalizedSymmetricCumulantsList);
}

void AliAnalysisTaskAR::BookQAHistograms() {
  // Book all QA histograms

  // book centrality estimator correlation QA histograms
  for (int cen = 0; cen < LAST_ECENESTIMATORS * (LAST_ECENESTIMATORS - 1) / 2;
       ++cen) {
    for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
      fCenCorQAHistograms[cen][ba] =
          new TH2D(fCenCorQAHistogramNames[cen][ba][kNAME],
                   fCenCorQAHistogramNames[cen][ba][kTITLE],
                   fCenCorQAHistogramBins[cen][kBIN],
                   fCenCorQAHistogramBins[cen][kLEDGE],
                   fCenCorQAHistogramBins[cen][kUEDGE],
                   fCenCorQAHistogramBins[cen][kBIN + LAST_EBINS],
                   fCenCorQAHistogramBins[cen][kLEDGE + LAST_EBINS],
                   fCenCorQAHistogramBins[cen][kUEDGE + LAST_EBINS]);
      fCenCorQAHistograms[cen][ba]->SetOption("colz");
      fCenCorQAHistograms[cen][ba]->GetXaxis()->SetTitle(
          fCenCorQAHistogramNames[cen][ba][kXAXIS]);
      fCenCorQAHistograms[cen][ba]->GetYaxis()->SetTitle(
          fCenCorQAHistogramNames[cen][ba][kYAXIS]);
      fCenCorQAHistogramsList->Add(fCenCorQAHistograms[cen][ba]);
    }
  }

  // book multiplicity estimator correlation QA histograms
  for (int mul = 0; mul < kMulEstimators * (kMulEstimators - 1) / 2; ++mul) {
    for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
      fMulCorQAHistograms[mul][ba] =
          new TH2D(fMulCorQAHistogramNames[mul][ba][kNAME],
                   fMulCorQAHistogramNames[mul][ba][kTITLE],
                   fMulCorQAHistogramBins[mul][kBIN],
                   fMulCorQAHistogramBins[mul][kLEDGE],
                   fMulCorQAHistogramBins[mul][kUEDGE],
                   fMulCorQAHistogramBins[mul][kBIN + LAST_EBINS],
                   fMulCorQAHistogramBins[mul][kLEDGE + LAST_EBINS],
                   fMulCorQAHistogramBins[mul][kUEDGE + LAST_EBINS]);
      fMulCorQAHistograms[mul][ba]->SetOption("colz");
      fMulCorQAHistograms[mul][ba]->GetXaxis()->SetTitle(
          fMulCorQAHistogramNames[mul][ba][kXAXIS]);
      fMulCorQAHistograms[mul][ba]->GetYaxis()->SetTitle(
          fMulCorQAHistogramNames[mul][ba][kYAXIS]);
      fMulCorQAHistogramsList->Add(fMulCorQAHistograms[mul][ba]);
    }
  }

  // book centrality-multiplicity correlation QA histograms
  for (int mul = 0; mul < kMulEstimators; ++mul) {
    for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
      fCenMulCorQAHistograms[mul][ba] =
          new TH2D(fCenMulCorQAHistogramNames[mul][ba][kNAME],
                   fCenMulCorQAHistogramNames[mul][ba][kTITLE],
                   fCenMulCorQAHistogramBins[mul][kBIN],
                   fCenMulCorQAHistogramBins[mul][kLEDGE],
                   fCenMulCorQAHistogramBins[mul][kUEDGE],
                   fCenMulCorQAHistogramBins[mul][kBIN + LAST_EBINS],
                   fCenMulCorQAHistogramBins[mul][kLEDGE + LAST_EBINS],
                   fCenMulCorQAHistogramBins[mul][kUEDGE + LAST_EBINS]);
      fCenMulCorQAHistograms[mul][ba]->SetOption("colz");
      fCenMulCorQAHistograms[mul][ba]->GetXaxis()->SetTitle(
          kCenEstimatorNames[fCentralityEstimator]);
      fCenMulCorQAHistograms[mul][ba]->GetYaxis()->SetTitle(
          fCenMulCorQAHistogramNames[mul][ba][kYAXIS]);
      fCenMulCorQAHistogramsList->Add(fCenMulCorQAHistograms[mul][ba]);
    }
  }

  if (fFillQACorHistogramsOnly) {
    return;
  }

  // book filter bit scan QA histogram
  fFBScanQAHistogram =
      new TH1D(fFBScanQAHistogramName[kNAME], fFBScanQAHistogramName[kTITLE],
               fFBScanQAHistogramBin[kBIN], fFBScanQAHistogramBin[kLEDGE],
               fFBScanQAHistogramBin[kUEDGE]);

  // set labels of filter bit scan QA histograms
  // filterbits are powers of 2, i.e 1,2,4,...
  // label the bins accordingly up to the hardcoded maximum filter bit
  int fb = 1;
  for (int i = 0; i < kMaxFilterbit; ++i) {
    fFBScanQAHistogram->GetXaxis()->SetBinLabel(i + 1, Form("%d", fb));
    fb *= 2;
  }
  fFBScanQAHistogram->SetFillColor(kFillColor[kAFTER]);
  fFBScanQAHistogramsList->Add(fFBScanQAHistogram);

  // book track scan filterbit QA histograms
  for (int track = 0; track < LAST_ETRACK; ++track) {
    for (int fb = 0; fb < kNumberofTestFilterBit; ++fb) {
      fFBTrackScanQAHistograms[track][fb] =
          new TH1D(fFBTrackScanQAHistogramNames[track][fb][kNAME],
                   fFBTrackScanQAHistogramNames[track][fb][kTITLE],
                   fFBTrackScanQAHistogramBins[track][kBIN],
                   fFBTrackScanQAHistogramBins[track][kLEDGE],
                   fFBTrackScanQAHistogramBins[track][kUEDGE]);
      fFBTrackScanQAHistograms[track][fb]->SetFillColor(kFillColor[kAFTER]);
      fFBTrackScanQAHistograms[track][fb]->GetXaxis()->SetTitle(
          fFBTrackScanQAHistogramNames[track][fb][kXAXIS]);
      fFBTrackScanQAHistograms[track][fb]->GetYaxis()->SetTitle(
          fFBTrackScanQAHistogramNames[track][fb][kYAXIS]);
      fFBScanQAHistogramsList->Add(fFBTrackScanQAHistograms[track][fb]);
    }
  }

  // book self correlation QA histograms
  for (int var = 0; var < kKinematic; ++var) {
    for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
      fSelfCorQAHistograms[var][ba] =
          new TH1D(fSelfCorQAHistogramNames[var][ba][kNAME],
                   fSelfCorQAHistogramNames[var][ba][kTITLE],
                   fSelfCorQAHistogramBins[var][kBIN],
                   fSelfCorQAHistogramBins[var][kLEDGE],
                   fSelfCorQAHistogramBins[var][kUEDGE]);
      fSelfCorQAHistograms[var][ba]->SetFillColor(kFillColor[ba]);
      fSelfCorQAHistograms[var][ba]->SetMinimum(0.1);
      fSelfCorQAHistograms[var][ba]->GetXaxis()->SetTitle(
          fSelfCorQAHistogramNames[var][ba][kXAXIS]);
      fSelfCorQAHistograms[var][ba]->GetYaxis()->SetTitle(
          fSelfCorQAHistogramNames[var][ba][kYAXIS]);
      fSelfCorQAHistogramsList->Add(fSelfCorQAHistograms[var][ba]);
    }
  }
}

void AliAnalysisTaskAR::BookControlHistograms() {
  // Book all control histograms

  // book histogram for counting trackcuts
  // add 5 bins manually for filterbit, charged only, primary only, global
  // only and MCClosure tracks only cut
  Int_t AddBins = 5;
  for (int mode = 0; mode < LAST_EMODE; ++mode) {
    fTrackCutsCounter[mode] =
        new TH1D(fTrackCutsCounterNames[mode], fTrackCutsCounterNames[mode],
                 2 * LAST_ETRACK + AddBins, 0, 2 * LAST_ETRACK + AddBins);
    fTrackCutsCounter[mode]->SetFillColor(kFillColor[kAFTER]);
    for (int bin = 0; bin < LAST_ETRACK; ++bin) {
      for (int mm = 0; mm < LAST_EMINMAX; ++mm) {
        fTrackCutsCounter[mode]->GetXaxis()->SetBinLabel(
            2 * bin + mm + 1, fTrackCutsCounterBinNames[bin][mm]);
      }
    }
    fTrackCutsCounter[mode]->GetXaxis()->SetBinLabel(2 * LAST_ETRACK + 1,
                                                     "Filterbit");
    fTrackCutsCounter[mode]->GetXaxis()->SetBinLabel(2 * LAST_ETRACK + 2,
                                                     "ChargedOnly");
    fTrackCutsCounter[mode]->GetXaxis()->SetBinLabel(2 * LAST_ETRACK + 3,
                                                     "PrimaryOnly");
    fTrackCutsCounter[mode]->GetXaxis()->SetBinLabel(2 * LAST_ETRACK + 4,
                                                     "GlobalTracksOnly");
    fTrackCutsCounter[mode]->GetXaxis()->SetBinLabel(2 * LAST_ETRACK + 5,
                                                     "MC Closure");
    fTrackControlHistogramsList->Add(fTrackCutsCounter[mode]);
  }
  // book histogram holding cumulative track cuts
  // Double_t ctcxmin[1] = {-0.5};
  // Double_t ctcxmax[1] = {std::pow(2, LAST_ETRACK + 3) - 0.5};
  // Int_t ctcbins[1] = {static_cast<Int_t>(::pow(2, LAST_ETRACK + 3))};
  // fTrackCutsCounterCumulative = new
  // THnSparseD(fTrackCutsCounterCumulativeName,
  //                                              fTrackCutsCounterCumulativeName,
  //                                              1, ctcbins, ctcxmin,
  //                                              ctcxmax);
  // fTrackControlHistogramsList->Add(fTrackCutsCounterCumulative);

  // book histogram holding values of all track cuts
  fTrackCutsValues =
      new TProfile(fTrackCutsValuesName, fTrackCutsValuesName,
                   2 * LAST_ETRACK + AddBins, 0, 2 * LAST_ETRACK + AddBins);
  fTrackCutsValues->SetFillColor(kFillColor[kAFTER]);

  for (int bin = 0; bin < LAST_ETRACK; ++bin) {
    for (int mm = 0; mm < LAST_EMINMAX; ++mm) {
      fTrackCutsValues->Fill(2 * bin + mm, fTrackCuts[bin][mm]);
      fTrackCutsValues->GetXaxis()->SetBinLabel(
          2 * bin + mm + 1, fTrackCutsCounterBinNames[bin][mm]);
    }
  }
  if (fUseFilterbit) {
    fTrackCutsValues->Fill(2 * LAST_ETRACK, fFilterbit);
  } else {
    fTrackCutsValues->Fill(2 * LAST_ETRACK, -99);
  }
  fTrackCutsValues->GetXaxis()->SetBinLabel(2 * LAST_ETRACK + 1, "Filterbit");
  if (fChargedOnly) {
    fTrackCutsValues->Fill(2 * LAST_ETRACK + 1, 99);
  } else {
    fTrackCutsValues->Fill(2 * LAST_ETRACK + 1, -99);
  }
  fTrackCutsValues->GetXaxis()->SetBinLabel(2 * LAST_ETRACK + 2, "ChargedOnly");
  if (fPrimaryOnly) {
    fTrackCutsValues->Fill(2 * LAST_ETRACK + 2, 99);
  } else {
    fTrackCutsValues->Fill(2 * LAST_ETRACK + 2, -99);
  }
  fTrackCutsValues->GetXaxis()->SetBinLabel(2 * LAST_ETRACK + 3, "PrimaryOnly");
  if (fGlobalTracksOnly) {
    fTrackCutsValues->Fill(2 * LAST_ETRACK + 3, 99);
  } else {
    fTrackCutsValues->Fill(2 * LAST_ETRACK + 3, -99);
  }
  fTrackCutsValues->GetXaxis()->SetBinLabel(2 * LAST_ETRACK + 4,
                                            "GlobalTracksOnly");
  if (fMCClosure) {
    fTrackCutsValues->Fill(2 * LAST_ETRACK + 4, 99);
  } else {
    fTrackCutsValues->Fill(2 * LAST_ETRACK + 4, -99);
  }
  fTrackCutsValues->GetXaxis()->SetBinLabel(2 * LAST_ETRACK + 5, "MC Closure");
  fTrackControlHistogramsList->Add(fTrackCutsValues);

  // book track control histograms
  for (int mode = 0; mode < LAST_EMODE; ++mode) {
    for (int var = 0; var < LAST_ETRACK; ++var) {
      for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
        fTrackControlHistograms[mode][var][ba] =
            new TH1D(fTrackControlHistogramNames[mode][var][ba][kNAME],
                     fTrackControlHistogramNames[mode][var][ba][kTITLE],
                     fTrackControlHistogramBins[var][kBIN],
                     fTrackControlHistogramBins[var][kLEDGE],
                     fTrackControlHistogramBins[var][kUEDGE]);
        fTrackControlHistograms[mode][var][ba]->SetFillColor(kFillColor[ba]);
        fTrackControlHistograms[mode][var][ba]->SetMinimum(0.1);
        fTrackControlHistograms[mode][var][ba]->GetXaxis()->SetTitle(
            fTrackControlHistogramNames[mode][var][ba][kXAXIS]);
        fTrackControlHistograms[mode][var][ba]->GetYaxis()->SetTitle(
            fTrackControlHistogramNames[mode][var][ba][kYAXIS]);
        fTrackControlHistogramsList->Add(
            fTrackControlHistograms[mode][var][ba]);
      }
    }
  }

  // book histogram for counting event cuts
  // add 5 bins by hand for centrality/multiplicity correlation cuts and
  // centrality flattening
  for (int mode = 0; mode < LAST_EMODE; ++mode) {
    fEventCutsCounter[mode] =
        new TH1D(fEventCutsCounterNames[mode], fEventCutsCounterNames[mode],
                 2 * (LAST_EEVENT + 2) + 1, 0, 2 * (LAST_EEVENT + 2) + 1);
    fEventCutsCounter[mode]->SetFillColor(kFillColor[kAFTER]);
    for (int bin = 0; bin < LAST_EEVENT; ++bin) {
      for (int mm = 0; mm < LAST_EMINMAX; ++mm) {
        fEventCutsCounter[mode]->GetXaxis()->SetBinLabel(
            2 * bin + mm + 1, fEventCutsCounterBinNames[bin][mm]);
      }
    }
    fEventCutsCounter[mode]->GetXaxis()->SetBinLabel(2 * LAST_EEVENT + 1,
                                                     "CenCorCut[kMIN]");
    fEventCutsCounter[mode]->GetXaxis()->SetBinLabel(2 * LAST_EEVENT + 2,
                                                     "CenCorCut[kMAX]");
    fEventCutsCounter[mode]->GetXaxis()->SetBinLabel(2 * (LAST_EEVENT + 1) + 1,
                                                     "MulCorCut[kMIN]");
    fEventCutsCounter[mode]->GetXaxis()->SetBinLabel(2 * (LAST_EEVENT + 1) + 2,
                                                     "MulCorCut[kMAX]");
    fEventCutsCounter[mode]->GetXaxis()->SetBinLabel(2 * (LAST_EEVENT + 2) + 1,
                                                     "fCenFlatten");
    fEventControlHistogramsList->Add(fEventCutsCounter[mode]);
  }
  // book histogram holding cumulative event cuts
  // Double_t cecxmin[1] = {-0.5};
  // Double_t cecxmax[1] = {std::pow(2, LAST_EEVENT + 2) - 0.5};
  // Int_t cecbins[1] = {static_cast<Int_t>(TMath::Power(2, LAST_EEVENT +
  // 2))}; fEventCutsCounterCumulative = new
  // THnSparseD(fEventCutsCounterCumulativeName,
  //                                              fEventCutsCounterCumulativeName,
  //                                              1, cecbins, cecxmin,
  //                                              cecxmax);
  // fEventControlHistogramsList->Add(fEventCutsCounterCumulative);
  // book histogram holding values of all event cuts
  fEventCutsValues = new TProfile(fEventCutsValuesName, fEventCutsValuesName,
                                  2 * LAST_EEVENT + 6, 0, 2 * LAST_EEVENT + 6);
  fEventCutsValues->SetFillColor(kFillColor[kAFTER]);

  for (int bin = 0; bin < LAST_EEVENT; ++bin) {
    for (int mm = 0; mm < LAST_EMINMAX; ++mm) {
      fEventCutsValues->Fill(2 * bin + mm, fEventCuts[bin][mm]);
      fEventCutsValues->GetXaxis()->SetBinLabel(
          2 * bin + mm + 1, fEventCutsCounterBinNames[bin][mm]);
    }
  }
  if (fUseCenCorCuts) {
    fEventCutsValues->Fill(2 * LAST_EEVENT + 0, fCenCorCut[0]);
    fEventCutsValues->Fill(2 * LAST_EEVENT + 1, fCenCorCut[1]);
  } else {
    fEventCutsValues->Fill(2 * LAST_EEVENT + 0, -999);
    fEventCutsValues->Fill(2 * LAST_EEVENT + 1, -999);
  }
  fEventCutsValues->GetXaxis()->SetBinLabel(2 * LAST_EEVENT + 1, "m_{CEN}");
  fEventCutsValues->GetXaxis()->SetBinLabel(2 * LAST_EEVENT + 2, "t_{CEN}");
  if (fUseMulCorCuts) {
    fEventCutsValues->Fill(2 * (LAST_EEVENT + 1) + 0, fMulCorCut[0]);
    fEventCutsValues->Fill(2 * (LAST_EEVENT + 1) + 1, fMulCorCut[1]);
  } else {
    fEventCutsValues->Fill(2 * (LAST_EEVENT + 1) + 0, -999);
    fEventCutsValues->Fill(2 * (LAST_EEVENT + 1) + 1, -999);
  }
  fEventCutsValues->GetXaxis()->SetBinLabel(2 * (LAST_EEVENT + 1) + 1,
                                            "m_{MUL}");
  fEventCutsValues->GetXaxis()->SetBinLabel(2 * (LAST_EEVENT + 1) + 2,
                                            "t_{MUL}");
  fEventCutsValues->GetXaxis()->SetBinLabel(2 * (LAST_EEVENT + 2) + 1,
                                            "CentralityEstimator");
  fEventCutsValues->Fill(2 * (LAST_EEVENT + 2) + 0, fCentralityEstimator);
  if (fUseCenFlatten) {
    fEventCutsValues->Fill(2 * (LAST_EEVENT + 2) + 1, 999);
  } else {
    fEventCutsValues->Fill(2 * (LAST_EEVENT + 2) + 1, -999);
  }
  fEventCutsValues->GetXaxis()->SetBinLabel(2 * (LAST_EEVENT + 2) + 2,
                                            "fUseCenFlatten");
  fEventControlHistogramsList->Add(fEventCutsValues);

  // book event control histograms
  for (int mode = 0; mode < LAST_EMODE; ++mode) {
    for (int var = 0; var < LAST_EEVENT; ++var) {
      for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
        fEventControlHistograms[mode][var][ba] =
            new TH1D(fEventControlHistogramNames[mode][var][ba][kNAME],
                     fEventControlHistogramNames[mode][var][ba][kTITLE],
                     fEventControlHistogramBins[var][kBIN],
                     fEventControlHistogramBins[var][kLEDGE],
                     fEventControlHistogramBins[var][kUEDGE]);
        fEventControlHistograms[mode][var][ba]->SetFillColor(kFillColor[ba]);
        fEventControlHistograms[mode][var][ba]->SetMinimum(0.1);
        fEventControlHistograms[mode][var][ba]->GetXaxis()->SetTitle(
            fEventControlHistogramNames[mode][var][ba][kXAXIS]);
        fEventControlHistograms[mode][var][ba]->GetYaxis()->SetTitle(
            fEventControlHistogramNames[mode][var][ba][kYAXIS]);
        fEventControlHistogramsList->Add(
            fEventControlHistograms[mode][var][ba]);
      }
    }
  }
}

void AliAnalysisTaskAR::BookFinalResultHistograms() {
  // Book final result histograms

  // book final result histograms
  for (int var = 0; var < LAST_EFINALHIST; ++var) {
    fFinalResultHistograms[var] =
        new TH1D(fFinalResultHistogramNames[var][kNAME],
                 fFinalResultHistogramNames[var][kTITLE],
                 fFinalResultHistogramBins[var][kBIN],
                 fFinalResultHistogramBins[var][kLEDGE],
                 fFinalResultHistogramBins[var][kUEDGE]);
    fFinalResultHistograms[var]->SetFillColor(kcolorFinalResult);
    fFinalResultHistograms[var]->GetXaxis()->SetTitle(
        fFinalResultHistogramNames[var][2]);
    fFinalResultHistogramsList->Add(fFinalResultHistograms[var]);
  }
}

void AliAnalysisTaskAR::BookFinalResultCorrelators() {
  // Book final result profiles holding correlators
  // 3 profiles for each correlator
  //  - integrated
  //  - as a function of centrality
  //  - as a function of multiplicity

  TList *corList;
  TProfile *corProfile[LAST_EFINALRESULTPROFILE];
  Double_t bins[LAST_EBINS][LAST_EFINALRESULTPROFILE] = {
      {1, 0, 1},
      {fEventControlHistogramBins[kCEN][kBIN],
       fEventControlHistogramBins[kCEN][kLEDGE],
       fEventControlHistogramBins[kCEN][kUEDGE]},
      {fEventControlHistogramBins[kMULQ][kBIN],
       fEventControlHistogramBins[kMULQ][kLEDGE],
       fEventControlHistogramBins[kMULQ][kUEDGE]},
  };
  TString Names[LAST_EFINALRESULTPROFILE] = {
      "[kINTEGRATED]", "[kCENDEP]", "[kMULDEP]", "[kPTDEP]", "[kETADEP]"};
  TString xaxis[LAST_EFINALRESULTPROFILE] = {"", "Centrality Percentile",
                                             "Multiplicity", "p_{T}", "#eta"};
  TString corListName;
  TString corName;
  for (std::size_t i = 0; i < fCorrelators.size(); i++) {

    corListName = "v_{";
    for (std::size_t j = 0; j < fCorrelators.at(i).size(); j++) {
      if (j != fCorrelators.at(i).size() - 1) {
        corListName += Form("%d,", fCorrelators.at(i).at(j));
      } else {
        corListName += Form("%d", fCorrelators.at(i).at(j));
      }
    }
    corListName += "}";
    corList = new TList();
    corList->SetName(corListName);
    corList->SetOwner(kTRUE);
    fFinalResultCorrelatorsList->Add(corList);

    for (int i = 0; i < LAST_EFINALRESULTPROFILE; i++) {
      if (i == kPTDEP) {
        corProfile[i] = new TProfile(
            corListName + Names[i], corListName + TString(" ") + Names[i],
            fTrackBins[kPT].size() - 1, fTrackBins[kPT].data());
      } else if (i == kETADEP) {
        corProfile[i] = new TProfile(
            corListName + Names[i], corListName + TString(" ") + Names[i],
            fTrackBins[kETA].size() - 1, fTrackBins[kETA].data());
      } else
        corProfile[i] = new TProfile(corListName + Names[i],
                                     corListName + TString(" ") + Names[i],
                                     bins[i][0], bins[i][1], bins[i][2]);
      corProfile[i]->GetXaxis()->SetTitle(xaxis[i]);
      corList->Add(corProfile[i]);
    }
  }
}
void AliAnalysisTaskAR::BookTrackBinHistograms() {
  // book track bin histograms
  TString Names[kKinematic - 1] = {"TrackBinHistogram[kPT]",
                                   "TrackBinHistogram[kETA]"};

  for (int i = 0; i < kKinematic - 1; i++) {
    if (fTrackBins[i].empty()) {
      fTrackBins[i].push_back(fTrackCuts[i][kMIN]);
      fTrackBins[i].push_back(fTrackCuts[i][kMAX]);
    }
    fTrackBinsHistogram[i] = new TH1D(
        Names[i], Names[i], fTrackBins[i].size() - 1, fTrackBins[i].data());

    for (std::size_t j = 0; j < fTrackBins[i].size() - 1; j++) {
      fKinematics[i].push_back({});
      fKinematicWeights[i].push_back({});
    }
  }

  fTrackBins[kPHI].push_back(0);
  fTrackBins[kPHI].push_back(1);
  fKinematics[kPHI].push_back({});
  fKinematicWeights[kPHI].push_back({});
}

void AliAnalysisTaskAR::BookFinalResultSymmetricCumulants() {
  // book histograms holding symmetric cumulants

  TH1D *Hist[LAST_EFINALRESULTPROFILE];
  TH1D *NHist[LAST_EFINALRESULTPROFILE];
  Double_t bins[LAST_EBINS][LAST_EFINALRESULTPROFILE] = {
      {1, 0, 1},
      {fEventControlHistogramBins[kCEN][kBIN],
       fEventControlHistogramBins[kCEN][kLEDGE],
       fEventControlHistogramBins[kCEN][kUEDGE]},
      {fEventControlHistogramBins[kMULQ][kBIN],
       fEventControlHistogramBins[kMULQ][kLEDGE],
       fEventControlHistogramBins[kMULQ][kUEDGE]},
  };
  TString Names[LAST_EFINALRESULTPROFILE] = {
      "[kINTEGRATED]", "[kCENDEP]", "[kMULDEP]", "[kPTDEP]", "[kETADEP]"};
  TString xaxis[LAST_EFINALRESULTPROFILE] = {"", "Centrality Percentile",
                                             "Multiplicity"};
  TString Name;
  TList *List, *ListNormalized;
  std::vector<std::vector<Int_t>> correlators;

  Int_t Index = 0;
  for (std::size_t j = 0; j < fSymmetricCumulants.size(); j++) {
    Name = "SC(";
    for (std::size_t i = 0; i < fSymmetricCumulants.at(j).size(); i++) {
      Name += Form("%d", fSymmetricCumulants.at(j).at(i));
      if (i < fSymmetricCumulants.at(j).size() - 1) {
        Name += ",";
      }
    }
    Name += ")";

    List = new TList();
    List->SetName(Name);
    List->SetOwner(kTRUE);
    fFinalResultSymmetricCumulantsList->Add(List);

    ListNormalized = new TList();
    ListNormalized->SetName(Form("N%s", Name.Data()));
    ListNormalized->SetOwner(kTRUE);
    fFinalResultNormalizedSymmetricCumulantsList->Add(ListNormalized);

    for (int i = 0; i < LAST_EFINALRESULTPROFILE; i++) {

      if (i == kPTDEP) {
        Hist[i] = new TH1D(Name + Names[i], Name + TString(" ") + Names[i],
                           fTrackBins[kPT].size() - 1, fTrackBins[kPT].data());
      } else if (i == kETADEP) {
        Hist[i] =
            new TH1D(Name + Names[i], Name + TString(" ") + Names[i],
                     fTrackBins[kETA].size() - 1, fTrackBins[kETA].data());
      } else {
        Hist[i] = new TH1D(Name + Names[i], Name + TString(" ") + Names[i],
                           bins[i][0], bins[i][1], bins[i][2]);
      }
      Hist[i]->GetXaxis()->SetTitle(xaxis[i]);
      Hist[i]->SetFillColor(kcolorFinalResult);
      List->Add(Hist[i]);

      NHist[i] =
          dynamic_cast<TH1D *>(Hist[i]->Clone(Form("N%s", Hist[i]->GetName())));
      NHist[i]->SetTitle(Form("N%s", Hist[i]->GetTitle()));
      ListNormalized->Add(NHist[i]);
    }

    correlators = MapSCToCor(fSymmetricCumulants.at(j));
    fMapSCtoCor.insert({fSymmetricCumulants.at(j), correlators});
    for (auto cor : correlators) {
      if (std::find(fCorrelators.begin(), fCorrelators.end(), cor) !=
          fCorrelators.end()) {
        continue;
      } else {
        fCorrelators.push_back(cor);
        fMapCorToIndex.insert({cor, Index});
        Index++;
      }
    }
  }
}

std::vector<std::vector<Int_t>>
AliAnalysisTaskAR::MapSCToCor(std::vector<Int_t> sc) {
  // map symmetric cumulant to the correlators needed for its computation

  std::vector<std::vector<Int_t>> correlators;
  switch (sc.size()) {
  case 2:
    correlators = {
        {-sc.at(0), -sc.at(1), sc.at(1), sc.at(0)}, // <4>_{-l,-k,k,l}
        {-sc.at(0), sc.at(0)},                      // <2> _{-k, k}
        {-sc.at(1), sc.at(1)}                       // <2> _{-l, l}
    };
    break;
  case 3:
    correlators = {
        {-sc.at(0), -sc.at(1), -sc.at(2), sc.at(2), sc.at(1),
         sc.at(0)},                                 // <6>_{-k,-l,-n,n,l,k}
        {-sc.at(0), -sc.at(1), sc.at(1), sc.at(0)}, // <4>_{-k,-l,l,k}
        {-sc.at(0), -sc.at(2), sc.at(2), sc.at(0)}, // <4>_{-k,-n,n,k}
        {-sc.at(1), -sc.at(2), sc.at(2), sc.at(1)}, // <4>_{-l,-n,n,l}
        {-sc.at(0), sc.at(0)},                      // <2> _{-k, k}
        {-sc.at(1), sc.at(1)},                      // <2> _{-l, l}
        {-sc.at(2), sc.at(2)}};                     // <2> _ { -n, n }
    break;
  default:
    std::cout << "higher order symmetric cumulants are not implemented (yet)"
              << std::endl;
  }
  return correlators;
}

void AliAnalysisTaskAR::SetDefaultConfiguration() {
  // set default configuration

  // do not fill QA histograms by default
  fFillQAHistograms = kFALSE;
  // not even the correlation histograms only
  fFillQACorHistogramsOnly = kFALSE;
  // use nested loops only for testing and with fixed multiplicity
  fUseNestedLoops = kFALSE;
  // only use custom seed if it is supplied with dedicated setter
  fUseCustomSeed = kFALSE;
  // do not generate data on the fly by default
  fMCOnTheFly = kFALSE;
  // only do a MC closure for validation
  fMCClosure = kFALSE;
  // do not use fischer yates by default
  fUseFisherYates = kFALSE;
  // needed in combination with fischer yates
  fUseFixedMultplicity = kFALSE;
  // use fake tracks for weight computation
  // i.e. misidentified tracks
  fUseFakeTracks = kTRUE;
}

void AliAnalysisTaskAR::SetDefaultBinning() {
  // set default binning for all histograms

  // default binning of track histograms
  Double_t DefaultTrackBins[LAST_ETRACK][LAST_EBINS] = {
      // kBIN kLEDGE kUEDGE
      {1000., 0., 20},            // kPT
      {200., -1., 1.},            // kETA
      {360., 0., TMath::TwoPi()}, // kPHI
      {21., -10.5, 10.5},         // kCHARGE
      {200., 0, 200},             // kTPCNCLS
      {200., 0, 200},             // kTPCCROSSEDROWS
      {2000., 0., 2.},            // kTPCNCLSFRACTIONSHARED
      {1000., 0., 10.},           // kTPCCHI2PERNDF
      {10., 0., 10.},             // kITSNCLS
      {1000., 0., 10.},           // kCHI2PERNDF
      {1000., -4., 4},            // kDCAZ
      {1000., -4., 4},            // kDCAXY
  };

  for (int track = 0; track < LAST_ETRACK; track++) {
    fTrackControlHistogramBins[track][kBIN] = DefaultTrackBins[track][kBIN];
    fTrackControlHistogramBins[track][kLEDGE] = DefaultTrackBins[track][kLEDGE];
    fTrackControlHistogramBins[track][kUEDGE] = DefaultTrackBins[track][kUEDGE];
  }

  // default binning of event histograms
  Double_t DefaultEventBins[LAST_EEVENT][LAST_EBINS] = {
      // kBIN kLEDGE kUEDGE
      {1000., 0., 20000}, // kMUL
      {1000., 0., 5000},  // kMULQ
      {1000., 0., 5000},  // kMULW
      {1000., 0., 5000},  // kMULREF
      {1000., 0., 5000},  // kNCONTRIB
      {100., 0., 100},    // kCEN
      {1000., -20., 20},  // kX
      {1000., -20., 20.}, // kY
      {1000., -20., 20.}, // kZ
      {1000., 0., 20.},   // kVPOS
  };

  for (int event = 0; event < LAST_EEVENT; event++) {
    fEventControlHistogramBins[event][kBIN] = DefaultEventBins[event][kBIN];
    fEventControlHistogramBins[event][kLEDGE] = DefaultEventBins[event][kLEDGE];
    fEventControlHistogramBins[event][kUEDGE] = DefaultEventBins[event][kUEDGE];
  }
}

void AliAnalysisTaskAR::SetDefaultCuts(Int_t Filterbit, Double_t cenMin,
                                       Double_t cenMax) {
  // set default track cuts depending on the filterbit

  // set filterbit
  fUseFilterbit = kTRUE;
  switch (Filterbit) {
  case 1:
    fFilterbit = Filterbit;
  case 92:
    fFilterbit = Filterbit;
  case 128:
    fFilterbit = Filterbit;

    // there are no cuts on dca with filterbit 128
    fTrackCuts[kDCAXY][kMIN] = -2.4;
    fTrackCuts[kDCAXY][kMAX] = 2.4;
    fUseTrackCuts[kDCAXY] = kTRUE;

    fTrackCuts[kDCAZ][kMIN] = -3.2;
    fTrackCuts[kDCAZ][kMAX] = 3.2;
    fUseTrackCuts[kDCAZ] = kTRUE;

    break;
  case 768:
    fFilterbit = Filterbit;
  default:
    std::cout << "Filterbit " << Filterbit << "not supported" << std::endl;
  }

  // default track cuts independent of filterbit
  fPrimaryOnly = kTRUE;
  fChargedOnly = kTRUE;
  // never use global tracks in combination with filterbit
  fGlobalTracksOnly = kFALSE;
  // only need for MC
  fMCPrimaryDef = kMCPhysicalPrim;

  fTrackCuts[kPT][kMIN] = 0.2;
  fTrackCuts[kPT][kMAX] = 5.;
  fUseTrackCuts[kPT] = kTRUE;

  fTrackCuts[kETA][kMIN] = -0.8;
  fTrackCuts[kETA][kMAX] = 0.8;
  fUseTrackCuts[kETA] = kTRUE;

  fTrackCuts[kPHI][kMIN] = 0.;
  fTrackCuts[kPHI][kMAX] = TMath::TwoPi();
  fUseTrackCuts[kPHI] = kTRUE;

  fTrackCuts[kCHARGE][kMIN] = -1.5;
  fTrackCuts[kCHARGE][kMAX] = 1.5;
  fUseTrackCuts[kCHARGE] = kTRUE;

  fTrackCuts[kTPCNCLS][kMIN] = 70;
  fTrackCuts[kTPCNCLS][kMAX] = 160;
  fUseTrackCuts[kTPCNCLS] = kTRUE;

  fTrackCuts[kTPCNCLS][kMIN] = 70;
  fTrackCuts[kTPCNCLS][kMAX] = 160;
  fUseTrackCuts[kTPCNCLS] = kTRUE;

  fTrackCuts[kTPCNCLSFRACTIONSHARED][kMIN] = 0;
  fTrackCuts[kTPCNCLSFRACTIONSHARED][kMAX] = 0.4;
  fUseTrackCuts[kTPCNCLSFRACTIONSHARED] = kTRUE;

  fTrackCuts[kTPCCHI2PERNDF][kMIN] = 0.;
  fTrackCuts[kTPCCHI2PERNDF][kMAX] = 4.;
  fUseTrackCuts[kTPCCHI2PERNDF] = kTRUE;

  // event cuts
  // cut only on multiplicity estimated by number of tracks in the Q-vector
  fEventCuts[kMULQ][kMIN] = 12.;
  fEventCuts[kMULQ][kMAX] = 3000.;
  fUseEventCuts[kMULQ] = kTRUE;

  fEventCuts[kNCONTRIB][kMIN] = 2.;
  fEventCuts[kNCONTRIB][kMAX] = 1e6;
  fUseEventCuts[kNCONTRIB] = kTRUE;

  // cut on multiplicity correlation
  fUseMulCorCuts = kTRUE;
  fMulCorCut[0] = 1.4;
  fMulCorCut[1] = 300.;

  // cut on vertex position
  fEventCuts[kX][kMIN] = -1.;
  fEventCuts[kX][kMAX] = 1.;
  fUseEventCuts[kX] = kTRUE;

  fEventCuts[kY][kMIN] = -1.;
  fEventCuts[kY][kMAX] = 1.;
  fUseEventCuts[kY] = kTRUE;

  fEventCuts[kZ][kMIN] = -10.;
  fEventCuts[kZ][kMAX] = 10.;
  fUseEventCuts[kZ] = kTRUE;

  fEventCuts[kVPOS][kMIN] = 1e-6;
  fEventCuts[kVPOS][kMAX] = 15.;
  fUseEventCuts[kVPOS] = kTRUE;

  // cut on centrality
  // should be overwritten in most cases, but keep the cut open by default
  fCentralityEstimator = kV0M;
  fEventCuts[kCEN][kMIN] = cenMin;
  fEventCuts[kCEN][kMAX] = cenMax;
  fUseEventCuts[kCEN] = kTRUE;

  // cut on centrality correlation
  fUseCenCorCuts = kTRUE;
  fCenCorCut[0] = 1.;
  fCenCorCut[1] = 10.;
}

void AliAnalysisTaskAR::UserExec(Option_t *) {

  // general strategy
  // Get pointer(s) to event: reconstructed, simulated or none
  // Fill event objects
  // Check event cut
  // Start Analysis
  // -> over AOD only or
  // -> over AOD and MC (TBI over MC only) or
  // -> generate monte carlo data on the fly
  // PostData

  // get pointer to AOD event
  AliAODEvent *aAOD = nullptr;
  aAOD = dynamic_cast<AliAODEvent *>(InputEvent());
  // get pointer to MC event
  AliMCEvent *aMC = nullptr;
  aMC = MCEvent();

  // fill event control and QA histograms when running over data
  if (!fMCOnTheFly) {
    // compute event objects
    // this requries an inital loop over all tracks in the event
    FillEventObjects(aAOD, aMC);
    // fill event histograms before cut
    if (fFillQAHistograms || fFillQACorHistogramsOnly) {
      FillEventQAHistograms(kBEFORE, aAOD, aMC);
    }
    FillEventControlHistograms(kBEFORE, aAOD);
    FillEventControlHistograms(kBEFORE, aMC);

    // check if event survives event cut
    if (!SurviveEventCut(aAOD)) {
      return;
    }

    // fill event histograms after cut
    FillEventControlHistograms(kAFTER, aAOD);
    FillEventControlHistograms(kAFTER, aMC);
    if (fFillQAHistograms || fFillQACorHistogramsOnly) {
      FillEventQAHistograms(kAFTER, aAOD, aMC);
    }
  }

  // start analysis

  // clear vectors holding kinematics and weights
  ClearVectors();

  // get number of all tracks in current event
  Int_t nTracks = 0;
  if (aAOD && !aMC) { // only running over AOD
    nTracks = aAOD->GetNumberOfTracks();
  } else if (aAOD && aMC) { // running over AOD and MC data
    nTracks = aMC->GetNumberOfTracks();
  } else if (fMCOnTheFly) { // running over on the fly generated MC data
    nTracks = TMath::Ceil(fMCMultiplicity->GetRandom());
  } else {
    std::cout << __LINE__ << ": did not get number of tracks" << std::endl;
    Fatal("UserExec", "did not get number of tracks in the event");
  }

  AliAODTrack *track = nullptr;
  AliAODMCParticle *MCParticle = nullptr;
  Int_t Counter = 0;

  for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) {

    // break loop if we hit fixed multiplicity
    if (fUseFixedMultplicity) {
      if (Counter >= fFixedMultiplicy) {
        break;
      }
    }

    // if we have AOD and MC data
    if (aAOD && aMC) {
      // get a pointer to a AliAODMCparticle
      MCParticle = dynamic_cast<AliAODMCParticle *>(aMC->GetTrack(iTrack));
      if (!MCParticle) {
        continue;
      }
      // and get corresponding AODTrack, if it exists
      track = dynamic_cast<AliAODTrack *>(
          aAOD->GetTrack(fLookUpTable[MCParticle->GetLabel()]));
      // maybe track was misidentify
      if (!track && fUseFakeTracks) {
        track = dynamic_cast<AliAODTrack *>(
            aAOD->GetTrack(fLookUpTable[TMath::Abs(MCParticle->GetLabel())]));
      }
      // if running over AOD only
    } else if (aAOD && !aMC) {
      // get AODtrack directly
      // randomize tracks if necessary
      if (fUseFisherYates) {
        track = dynamic_cast<AliAODTrack *>(
            aAOD->GetTrack(fRandomizedTrackIndices.at(iTrack)));
      } else {
        track = dynamic_cast<AliAODTrack *>(aAOD->GetTrack(iTrack));
      }
      // run over on the fly generated
    } else if (fMCOnTheFly) {
      if (fMCKinematicPDFs[kPT]) {
        fMCKinematicVariables[kPT] = fMCKinematicPDFs[kPT]->GetRandom();
      }
      if (fMCKinematicPDFs[kETA]) {
        fMCKinematicVariables[kETA] = fMCKinematicPDFs[kETA]->GetRandom();
      }
      if (fMCKinematicPDFs[kPHI]) {
        fMCKinematicVariables[kPHI] = fMCKinematicPDFs[kPHI]->GetRandom();
      }
    } else {
      std::cout << __LINE__ << ": did not get kinematic variables" << std::endl;
      Fatal("UserExec", "did not get kinematic variables");
    }

    // run over data
    if (!fMCOnTheFly) {

      // fill QA track scan histograms
      if (fFillQAHistograms && !fFillQACorHistogramsOnly) {
        FillFBScanQAHistograms(track);
      }

      // fill control histogram before cutting
      FillTrackControlHistograms(kBEFORE, MCParticle);
      FillTrackControlHistograms(kBEFORE, track);

      // cut on monte carlo data, if we have any
      if (!SurviveTrackCut(MCParticle, kTRUE)) {
        continue;
      }

      // fill track control histogram after track cut on MC particle
      FillTrackControlHistograms(kAFTER, MCParticle);

      // cut on reconstructed track
      if (!SurviveTrackCut(track, kTRUE)) {
        continue;
      }

      // fill track control histogram after track cut on reconstructed track
      FillTrackControlHistograms(kAFTER, track);

      // run over on the fly generated data
    } else {

      FillTrackControlHistograms(kBEFORE, nullptr);

      if (!SurviveTrackCut(nullptr, kFALSE)) {
        continue;
      }

      FillTrackControlHistograms(kAFTER, nullptr);
    }

    // fill kinematic variables and weights into track objects
    FillTrackObjects(track);
    // increase counter, in case we want to fix multiplicity
    Counter++;
  }

  // bail out if we only want to fill control histograms
  if (fFillControlHistogramsOnly) {
    return;
  }

  // fill event control histograms for on the fly generated data
  // if we run over real data, we can do this step in the very beginning
  if (fMCOnTheFly) {
    fMultiplicity[kMUL] = nTracks;
    fMultiplicity[kMULQ] = fKinematics[kPHI].at(0).size();
    fMultiplicity[kMULW] =
        std::accumulate(fKinematicWeights[kPHI].at(0).begin(),
                        fKinematicWeights[kPHI].at(0).end(), 0);
    FillEventControlHistograms(kBEFORE, nullptr);
    FillEventControlHistograms(kAFTER, nullptr);
  }

  // fill final result profile
  FillFinalResultCorrelators();

  PostData(1, fHistList);
}

void AliAnalysisTaskAR::ClearVectors() {
  // clear vectors holding kinematics and weights of an event
  for (int k = 0; k < kKinematic - 1; ++k) {
    for (std::size_t i = 0; i < fTrackBins[k].size() - 1; i++) {
      fKinematics[k].at(i).clear();
      fKinematicWeights[k].at(i).clear();
    }
  }
  fKinematics[kPHI].at(0).clear();
  fKinematicWeights[kPHI].at(0).clear();
}

void AliAnalysisTaskAR::FillTrackObjects(AliVParticle *avp) {
  // fill kinematic variables and weights into event objects

  Double_t weight = 1.;
  // AOD track
  AliAODTrack *track = dynamic_cast<AliAODTrack *>(avp);
  if (track) {
    if (fUseWeights) {
      if (fWeightHistogram[kPT]) {
        weight *= fWeightHistogram[kPT]->GetBinContent(
            fWeightHistogram[kPT]->FindBin(track->Pt()));
      }
      if (fWeightHistogram[kETA]) {
        weight *= fWeightHistogram[kETA]->GetBinContent(
            fWeightHistogram[kETA]->FindBin(track->Eta()));
      }
      if (fWeightHistogram[kPHI]) {
        weight *= fWeightHistogram[kPHI]->GetBinContent(
            fWeightHistogram[kPHI]->FindBin(track->Phi()));
      }
    }
    fKinematics[kPHI].at(0).push_back(track->Phi());
    fKinematicWeights[kPHI].at(0).push_back(weight);

    Int_t ptBin = fTrackBinsHistogram[kPT]->FindBin(track->Pt());
    if (ptBin != 0 && ptBin <= fTrackBinsHistogram[kPT]->GetNbinsX()) {
      fKinematics[kPT].at(ptBin - 1).push_back(track->Phi());
      fKinematicWeights[kPT].at(ptBin - 1).push_back(weight);
    }
    Int_t etaBin = fTrackBinsHistogram[kETA]->FindBin(track->Eta());
    if (etaBin != 0 && etaBin <= fTrackBinsHistogram[kETA]->GetNbinsX()) {
      fKinematics[kETA].at(etaBin - 1).push_back(track->Phi());
      fKinematicWeights[kETA].at(etaBin - 1).push_back(weight);
    }
  }

  // MC particles
  AliAODMCParticle *MCparticle = dynamic_cast<AliAODMCParticle *>(avp);
  if (MCparticle) {
    if (fUseWeights) {
      if (fWeightHistogram[kPT]) {
        weight *= fWeightHistogram[kPT]->GetBinContent(
            fWeightHistogram[kPT]->FindBin(MCparticle->Pt()));
      }
      if (fWeightHistogram[kETA]) {
        weight *= fWeightHistogram[kETA]->GetBinContent(
            fWeightHistogram[kETA]->FindBin(MCparticle->Eta()));
      }
      if (fWeightHistogram[kPHI]) {
        weight *= fWeightHistogram[kPHI]->GetBinContent(
            fWeightHistogram[kPHI]->FindBin(MCparticle->Phi()));
      }
    }
    fKinematics[kPHI].at(0).push_back(MCparticle->Phi());
    fKinematicWeights[kPHI].at(0).push_back(weight);

    Int_t ptBin = fTrackBinsHistogram[kPT]->FindBin(MCparticle->Pt());
    if (ptBin != 0 && ptBin <= fTrackBinsHistogram[kPT]->GetNbinsX()) {
      fKinematics[kPT].at(ptBin - 1).push_back(MCparticle->Phi());
      fKinematicWeights[kPT].at(ptBin - 1).push_back(weight);
    }
    Int_t etaBin = fTrackBinsHistogram[kETA]->FindBin(MCparticle->Eta());
    if (etaBin != 0 && etaBin <= fTrackBinsHistogram[kETA]->GetNbinsX()) {
      fKinematics[kETA].at(etaBin - 1).push_back(MCparticle->Phi());
      fKinematicWeights[kETA].at(etaBin - 1).push_back(weight);
    }
  }

  // if neither, we run over on the fly generated data
  if (!track && !MCparticle && fMCOnTheFly) {
    if (fUseWeights) {
      if (fWeightHistogram[kPT]) {
        weight *= fWeightHistogram[kPT]->GetBinContent(
            fWeightHistogram[kPT]->FindBin(fMCKinematicVariables[kPT]));
      }
      if (fWeightHistogram[kETA]) {
        weight *= fWeightHistogram[kETA]->GetBinContent(
            fWeightHistogram[kETA]->FindBin(fMCKinematicVariables[kETA]));
      }
      if (fWeightHistogram[kPHI]) {
        weight *= fWeightHistogram[kPHI]->GetBinContent(
            fWeightHistogram[kPHI]->FindBin(fMCKinematicVariables[kPHI]));
      }
    }
    fKinematics[kPHI].at(0).push_back(fMCKinematicVariables[kPHI]);
    fKinematicWeights[kPHI].at(0).push_back(weight);

    Int_t ptBin = fTrackBinsHistogram[kPT]->FindBin(fMCKinematicVariables[kPT]);
    if (ptBin > 0 && ptBin <= fTrackBinsHistogram[kPT]->GetNbinsX()) {
      fKinematics[kPT].at(ptBin - 1).push_back(fMCKinematicVariables[kPHI]);
      fKinematicWeights[kPT].at(ptBin - 1).push_back(weight);
    }
    Int_t etaBin =
        fTrackBinsHistogram[kETA]->FindBin(fMCKinematicVariables[kETA]);
    if (etaBin > 0 && etaBin <= fTrackBinsHistogram[kETA]->GetNbinsX()) {
      fKinematics[kETA].at(etaBin - 1).push_back(fMCKinematicVariables[kPHI]);
      fKinematicWeights[kETA].at(etaBin - 1).push_back(weight);
    }
  }
}

Int_t AliAnalysisTaskAR::IndexCorHistograms(Int_t i, Int_t j, Int_t N) {
  // helper function for computing index of correlation histograms
  // this function projects 2D indeces of a matrix above the diagonal
  // to a 1D index
  //
  // example with N=4
  //    i->
  // j( 00 01 02 03)
  // |( 10 11 12 13)
  // v( 20 21 22 23)
  //  ( 30 31 32 33)
  //
  // Entry 01: IndexCorHistograms(0,1,4) => 0
  // Entry 02: IndexCorHistograms(0,2,4) => 1
  // Entry 03: IndexCorHistograms(0,3,4) => 2
  // Entry 12: IndexCorHistograms(1,2,4) => 3
  // Entry 13: IndexCorHistograms(1,3,4) => 4
  // Entry 23: IndexCorHistograms(2,3,4) => 5
  // 4*(4-1)/2=6 entries

  Int_t Index = 0;
  for (int k = 0; k < i; ++k) {
    Index += N - (k + 1);
  }
  Index += j - i - 1;
  return Index;
}

void AliAnalysisTaskAR::FillEventControlHistograms(kBeforeAfter BA,
                                                   AliVEvent *Event) {
  // fill event control histograms

  // AOD event
  AliAODEvent *AODEvent = dynamic_cast<AliAODEvent *>(Event);
  if (AODEvent) {

    // get primary vertex object
    AliAODVertex *PrimaryVertex = AODEvent->GetPrimaryVertex();

    fEventControlHistograms[kRECO][kMUL][BA]->Fill(fMultiplicity[kMUL]);
    fEventControlHistograms[kRECO][kMULQ][BA]->Fill(fMultiplicity[kMULQ]);
    fEventControlHistograms[kRECO][kMULW][BA]->Fill(fMultiplicity[kMULW]);
    fEventControlHistograms[kRECO][kMULREF][BA]->Fill(fMultiplicity[kMULREF]);
    fEventControlHistograms[kRECO][kNCONTRIB][BA]->Fill(
        fMultiplicity[kNCONTRIB]);
    fEventControlHistograms[kRECO][kCEN][BA]->Fill(
        fCentrality[fCentralityEstimator]);
    fEventControlHistograms[kRECO][kX][BA]->Fill(PrimaryVertex->GetX());
    fEventControlHistograms[kRECO][kY][BA]->Fill(PrimaryVertex->GetY());
    fEventControlHistograms[kRECO][kZ][BA]->Fill(PrimaryVertex->GetZ());
    fEventControlHistograms[kRECO][kVPOS][BA]->Fill(
        std::sqrt(PrimaryVertex->GetX() * PrimaryVertex->GetX() +
                  PrimaryVertex->GetY() * PrimaryVertex->GetY() +
                  PrimaryVertex->GetZ() * PrimaryVertex->GetZ()));
  }

  // MC event
  AliMCEvent *MCEvent = dynamic_cast<AliMCEvent *>(Event);
  if (MCEvent) {
    fEventControlHistograms[kSIM][kMUL][BA]->Fill(MCEvent->GetNumberOfTracks());
  }

  if (!AODEvent && !MCEvent && fMCOnTheFly) {
    fEventControlHistograms[kSIM][kMUL][BA]->Fill(fMultiplicity[kMUL]);
    fEventControlHistograms[kSIM][kMULQ][BA]->Fill(fMultiplicity[kMULQ]);
    fEventControlHistograms[kSIM][kMULW][BA]->Fill(fMultiplicity[kMULW]);
  }
}

void AliAnalysisTaskAR::FillTrackControlHistograms(kBeforeAfter BA,
                                                   AliVParticle *avp) {
  // fill track control histograms

  // aod track
  AliAODTrack *track = dynamic_cast<AliAODTrack *>(avp);
  if (track) {
    fTrackControlHistograms[kRECO][kPT][BA]->Fill(track->Pt());
    fTrackControlHistograms[kRECO][kETA][BA]->Fill(track->Eta());
    fTrackControlHistograms[kRECO][kPHI][BA]->Fill(track->Phi());
    fTrackControlHistograms[kRECO][kCHARGE][BA]->Fill(track->Charge());
    fTrackControlHistograms[kRECO][kTPCNCLS][BA]->Fill(track->GetTPCNcls());
    fTrackControlHistograms[kRECO][kTPCCROSSEDROWS][BA]->Fill(
        track->GetTPCCrossedRows());
    if (track->GetTPCNcls() != 0) {
      fTrackControlHistograms[kRECO][kTPCNCLSFRACTIONSHARED][BA]->Fill(
          (Double_t)track->GetTPCnclsS() / (Double_t)track->GetTPCNcls());
      fTrackControlHistograms[kRECO][kTPCCHI2PERNDF][BA]->Fill(
          track->GetTPCchi2() / track->GetTPCNcls());
    }
    fTrackControlHistograms[kRECO][kITSNCLS][BA]->Fill(track->GetITSNcls());
    fTrackControlHistograms[kRECO][kCHI2PERNDF][BA]->Fill(track->Chi2perNDF());
    fTrackControlHistograms[kRECO][kDCAZ][BA]->Fill(track->ZAtDCA());
    fTrackControlHistograms[kRECO][kDCAXY][BA]->Fill(track->DCA());
  }

  // MC particle
  AliAODMCParticle *MCParticle = dynamic_cast<AliAODMCParticle *>(avp);
  if (MCParticle) {
    fTrackControlHistograms[kSIM][kPT][BA]->Fill(MCParticle->Pt());
    fTrackControlHistograms[kSIM][kETA][BA]->Fill(MCParticle->Eta());
    fTrackControlHistograms[kSIM][kPHI][BA]->Fill(MCParticle->Phi());
    fTrackControlHistograms[kSIM][kCHARGE][BA]->Fill(MCParticle->Charge() / 3.);
  }

  // MC on the fly
  if (!avp && fMCOnTheFly) {
    fTrackControlHistograms[kSIM][kPT][BA]->Fill(fMCKinematicVariables[kPT]);
    fTrackControlHistograms[kSIM][kETA][BA]->Fill(fMCKinematicVariables[kETA]);
    fTrackControlHistograms[kSIM][kPHI][BA]->Fill(fMCKinematicVariables[kPHI]);
  }
}

void AliAnalysisTaskAR::FillEventQAHistograms(kBeforeAfter BA,
                                              AliAODEvent *AODEvent,
                                              AliMCEvent *MCEvent) {
  // fill event QA control histograms

  if (AODEvent) {
    // fill centrality estimator correlation histograms
    for (int i = 0; i < LAST_ECENESTIMATORS; ++i) {
      for (int j = i + 1; j < LAST_ECENESTIMATORS; ++j) {
        fCenCorQAHistograms[IndexCorHistograms(i, j, LAST_ECENESTIMATORS)][BA]
            ->Fill(fCentrality[i], fCentrality[j]);
      }
    }

    // fill multiplicity correlation histograms
    for (int i = 0; i < kMulEstimators; ++i) {
      for (int j = i + 1; j < kMulEstimators; ++j) {
        fMulCorQAHistograms[IndexCorHistograms(i, j, kMulEstimators)][BA]->Fill(
            fMultiplicity[i], fMultiplicity[j]);
      }
    }

    // fill centrality-multiplicity correlation histograms
    for (int i = 0; i < kMulEstimators; ++i) {
      fCenMulCorQAHistograms[i][BA]->Fill(fCentrality[fCentralityEstimator],
                                          fMultiplicity[i]);
    }

    if (fFillQACorHistogramsOnly) {
      return;
    }

    // search for self correlations with nested loop
    Int_t nTracks = AODEvent->GetNumberOfTracks();
    AliAODTrack *aTrack1 = nullptr;
    AliAODTrack *aTrack2 = nullptr;
    // starting a loop over the first track
    for (Int_t iTrack1 = 0; iTrack1 < nTracks; iTrack1++) {
      aTrack1 = dynamic_cast<AliAODTrack *>(AODEvent->GetTrack(iTrack1));
      if (!aTrack1 || !SurviveTrackCut(aTrack1, kFALSE)) {
        continue;
      }
      // starting a loop over the second track
      for (Int_t iTrack2 = iTrack1 + 1; iTrack2 < nTracks; iTrack2++) {
        aTrack2 = dynamic_cast<AliAODTrack *>(AODEvent->GetTrack(iTrack2));
        if (!aTrack2 || !SurviveTrackCut(aTrack2, kFALSE)) {
          continue;
        }
        // compute differences
        fSelfCorQAHistograms[kPT][BA]->Fill(aTrack1->Pt() - aTrack2->Pt());
        fSelfCorQAHistograms[kETA][BA]->Fill(aTrack1->Eta() - aTrack2->Eta());
        fSelfCorQAHistograms[kPHI][BA]->Fill(aTrack1->Phi() - aTrack2->Phi());
      }
    }
  }
  if (MCEvent) {
    // TBI
  }
}

void AliAnalysisTaskAR::FillFBScanQAHistograms(AliAODTrack *track) {
  // fill filter bit scan QA histograms

  // check for filterbits of the track
  // filterbits are powers of 2, i.e. 1,2,4,8,...
  int fb = 1;
  for (int i = 0; i < kMaxFilterbit; ++i) {
    if (track->TestFilterBit(fb)) {
      fFBScanQAHistogram->Fill(i);
    }
    fb *= 2;
  }

  // scan kinematic variables for different filterbits
  for (int fb = 0; fb < kNumberofTestFilterBit; ++fb) {
    if (track->TestFilterBit(kTestFilterbit[fb])) {
      fFBTrackScanQAHistograms[kPT][fb]->Fill(track->Pt());
      fFBTrackScanQAHistograms[kETA][fb]->Fill(track->Eta());
      fFBTrackScanQAHistograms[kPHI][fb]->Fill(track->Phi());
      fFBTrackScanQAHistograms[kCHARGE][fb]->Fill(track->Charge());
      fFBTrackScanQAHistograms[kTPCNCLS][fb]->Fill(track->GetTPCNcls());
      fFBTrackScanQAHistograms[kTPCCROSSEDROWS][fb]->Fill(
          track->GetTPCCrossedRows());
      if (track->GetTPCNcls() != 0) {
        fFBTrackScanQAHistograms[kTPCNCLSFRACTIONSHARED][fb]->Fill(
            (Double_t)track->GetTPCnclsS() / (Double_t)track->GetTPCNcls());
        fFBTrackScanQAHistograms[kTPCCHI2PERNDF][fb]->Fill(track->GetTPCchi2() /
                                                           track->GetTPCNcls());
      }
      fFBTrackScanQAHistograms[kITSNCLS][fb]->Fill(track->GetITSNcls());
      fFBTrackScanQAHistograms[kCHI2PERNDF][fb]->Fill(track->Chi2perNDF());
      fFBTrackScanQAHistograms[kDCAZ][fb]->Fill(track->ZAtDCA());
      fFBTrackScanQAHistograms[kDCAXY][fb]->Fill(track->DCA());
    }
  }
}

Bool_t AliAnalysisTaskAR::SurviveEventCut(AliAODEvent *aAOD) {
  // Check if the current event survives event cuts
  // return flag at the end, if one cut is not passed, set it to kFALSE
  Bool_t Flag = kTRUE;
  // Double_t CutBit = 0;

  // only cut if we pass a valid pointer
  if (aAOD) {

    // cut on multiplicity
    // number of total tracks of the event
    if (fUseEventCuts[kMUL]) {
      if (fMultiplicity[kMUL] < fEventCuts[kMUL][kMIN]) {
        fEventCutsCounter[kRECO]->Fill(2 * kMUL + kMIN + 0.5);
        // CutBit += TMath::Power(2, kMUL);
        Flag = kFALSE;
      }
      if (fMultiplicity[kMUL] > fEventCuts[kMUL][kMAX]) {
        fEventCutsCounter[kRECO]->Fill(2 * kMUL + kMAX + 0.5);
        // CutBit += TMath::Power(2, kMUL);
        Flag = kFALSE;
      }
    }
    // number of tracks that survive track cuts
    if (fUseEventCuts[kMULQ]) {
      if (fMultiplicity[kMULQ] < fEventCuts[kMULQ][kMIN]) {
        fEventCutsCounter[kRECO]->Fill(2 * kMULQ + kMIN + 0.5);
        // CutBit += TMath::Power(2, kMULQ);
        Flag = kFALSE;
      }
      if (fMultiplicity[kMULQ] > fEventCuts[kMULQ][kMAX]) {
        fEventCutsCounter[kRECO]->Fill(2 * kMULQ + kMAX + 0.5);
        // CutBit += TMath::Power(2, kMULQ);
        Flag = kFALSE;
      }
    }
    if (fUseEventCuts[kMULW]) {
      // sum of weighted surviving tracks
      if (fMultiplicity[kMULW] < fEventCuts[kMULW][kMIN]) {
        fEventCutsCounter[kRECO]->Fill(2 * kMULW + kMIN + 0.5);
        // CutBit += TMath::Power(2, kMULW);
        Flag = kFALSE;
      }
      if (fMultiplicity[kMULW] > fEventCuts[kMULW][kMAX]) {
        fEventCutsCounter[kRECO]->Fill(2 * kMULW + kMAX + 0.5);
        // CutBit += TMath::Power(2, kMULW);
        Flag = kFALSE;
      }
    }
    if (fUseEventCuts[kNCONTRIB]) {
      // numbers of contriubters to the vertex
      if (fMultiplicity[kNCONTRIB] < fEventCuts[kNCONTRIB][kMIN]) {
        fEventCutsCounter[kRECO]->Fill(2 * kNCONTRIB + kMIN + 0.5);
        // CutBit += TMath::Power(2, kNCONTRIB);
        Flag = kFALSE;
      }
      if (fMultiplicity[kNCONTRIB] > fEventCuts[kNCONTRIB][kMAX]) {
        fEventCutsCounter[kRECO]->Fill(2 * kNCONTRIB + kMAX + 0.5);
        // CutBit += TMath::Power(2, kNCONTRIB);
        Flag = kFALSE;
      }
    }
    if (fUseEventCuts[kMULREF]) {
      // cut event if it is not within the reference centrality percentile
      if (fMultiplicity[kMULREF] < fEventCuts[kMULREF][kMIN]) {
        fEventCutsCounter[kRECO]->Fill(2 * kMULREF + kMIN + 0.5);
        // CutBit += TMath::Power(2, kMULREF);
        Flag = kFALSE;
      }
      if (fMultiplicity[kMULREF] > fEventCuts[kMULREF][kMAX]) {
        fEventCutsCounter[kRECO]->Fill(2 * kMULREF + kMAX + 0.5);
        // CutBit += TMath::Power(2, kMULREF);
        Flag = kFALSE;
      }
    }
    if (fUseEventCuts[kCEN]) {
      // cut event if it is not within the centrality percentile
      if (fCentrality[fCentralityEstimator] < fEventCuts[kCEN][kMIN]) {
        fEventCutsCounter[kRECO]->Fill(2 * kCEN + kMIN + 0.5);
        // CutBit += TMath::Power(2, kCEN);
        Flag = kFALSE;
      }
      if (fCentrality[fCentralityEstimator] > fEventCuts[kCEN][kMAX]) {
        fEventCutsCounter[kRECO]->Fill(2 * kCEN + kMAX + 0.5);
        // CutBit += TMath::Power(2, kCEN);
        Flag = kFALSE;
      }
    }
    // Get primary vertex
    AliAODVertex *PrimaryVertex = aAOD->GetPrimaryVertex();
    if (!PrimaryVertex) {
      return kFALSE;
    }
    if (fUseEventCuts[kX]) {
      // cut event if primary vertex is too out of center
      if (PrimaryVertex->GetX() < fEventCuts[kX][kMIN]) {
        fEventCutsCounter[kRECO]->Fill(2 * kX + kMIN + 0.5);
        // CutBit += TMath::Power(2, kX);
        Flag = kFALSE;
      }
      if (PrimaryVertex->GetX() > fEventCuts[kX][kMAX]) {
        fEventCutsCounter[kRECO]->Fill(2 * kX + kMAX + 0.5);
        // CutBit += TMath::Power(2, kX);
        Flag = kFALSE;
      }
    }
    if (fUseEventCuts[kY]) {
      if (PrimaryVertex->GetY() < fEventCuts[kY][kMIN]) {
        fEventCutsCounter[kRECO]->Fill(2 * kY + kMIN + 0.5);
        // CutBit += TMath::Power(2, kY);
        Flag = kFALSE;
      }
      if (PrimaryVertex->GetY() > fEventCuts[kY][kMAX]) {
        fEventCutsCounter[kRECO]->Fill(2 * kY + kMAX + 0.5);
        // CutBit += TMath::Power(2, kY);
        Flag = kFALSE;
      }
    }
    if (fUseEventCuts[kZ]) {
      if (PrimaryVertex->GetZ() < fEventCuts[kZ][kMIN]) {
        fEventCutsCounter[kRECO]->Fill(2 * kZ + kMIN + 0.5);
        // CutBit += TMath::Power(2, kZ);
        Flag = kFALSE;
      }
      if (PrimaryVertex->GetZ() > fEventCuts[kZ][kMAX]) {
        fEventCutsCounter[kRECO]->Fill(2 * kZ + kMAX + 0.5);
        // CutBit += TMath::Power(2, kZ);
        Flag = kFALSE;
      }
    }
    if (fUseEventCuts[kVPOS]) {
      // additionally cut on absolute value of the vertex postion
      // there are suspicous events with |r_v|=0 that we do not trust
      Double_t VPos = std::sqrt(PrimaryVertex->GetX() * PrimaryVertex->GetX() +
                                PrimaryVertex->GetY() * PrimaryVertex->GetY() +
                                PrimaryVertex->GetZ() * PrimaryVertex->GetZ());
      if (VPos < fEventCuts[kVPOS][kMIN]) {
        fEventCutsCounter[kRECO]->Fill(2 * kVPOS + kMIN + 0.5);
        // CutBit += TMath::Power(2, kVPOS);
        Flag = kFALSE;
      }
      if (VPos > fEventCuts[kVPOS][kMAX]) {
        fEventCutsCounter[kRECO]->Fill(2 * kVPOS + kMAX + 0.5);
        // CutBit += TMath::Power(2, kVPOS);
        Flag = kFALSE;
      }
    }
    if (fUseCenCorCuts) {
      // cut on centrality estimator correlation
      // ugly! cut on fundamental observerables instead but there are some
      // really weird events we need to get rid off
      // cut away all events that are above the line
      // y=mx+t
      // and below
      // y=(x-t)/m
      // this gives a nice and symmetric cone around the diagonal y=x
      // set m>1 such that the cone gets wider for larger centralities
      Double_t m_cen = fCenCorCut[0];
      Double_t t_cen = fCenCorCut[1];
      for (int i = 0; i < LAST_ECENESTIMATORS; ++i) {
        for (int j = i + 1; j < LAST_ECENESTIMATORS; ++j) {
          if (fCentrality[j] > m_cen * fCentrality[i] + t_cen) {
            fEventCutsCounter[kRECO]->Fill(2 * LAST_EEVENT + kMAX + 0.5);
            // CutBit += TMath::Power(2, LAST_EEVENT + 0);
            Flag = kFALSE;
          }
          if (fCentrality[j] < (fCentrality[i] - t_cen) / m_cen) {
            fEventCutsCounter[kRECO]->Fill(2 * LAST_EEVENT + kMIN + 0.5);
            // CutBit += TMath::Power(2, LAST_EEVENT + 0);
            Flag = kFALSE;
          }
        }
      }
    }
    if (fUseMulCorCuts) {
      // cut on multiplicity correlation
      // ugly! cut on fundamental observerables instead but there are some
      // really weird events we need to get rid off
      // logic is same as above
      Double_t m_mul = fMulCorCut[0];
      Double_t t_mul = fMulCorCut[1];
      for (int i = 0; i < kMulEstimators; ++i) {
        for (int j = i + 1; j < kMulEstimators; ++j) {
          // skip kMUL since it is a bad multiplicity estimate
          // skip kMULW since it will differ greatly from kMULQ when we have
          // to use large weights
          if (i == kMUL || j == kMUL || i == kMULW || j == kMULW) {
            continue;
            ;
          }
          if (fMultiplicity[j] > m_mul * fMultiplicity[i] + t_mul) {
            fEventCutsCounter[kRECO]->Fill(2 * (LAST_EEVENT + 1) + kMAX + 0.5);
            // CutBit += TMath::Power(2, LAST_EEVENT + 1);
            Flag = kFALSE;
          }
          if (fMultiplicity[j] < (fMultiplicity[i] - t_mul) / m_mul) {
            fEventCutsCounter[kRECO]->Fill(2 * (LAST_EEVENT + 1) + kMIN + 0.5);
            // CutBit += TMath::Power(2, LAST_EEVENT + 1);
            Flag = kFALSE;
          }
        }
      }
    }

    if (fUseCenFlatten) {
      // flatten centrality
      // find acceptance probability as a function of centrality in
      // fCenFlattenHist
      Double_t CenProb = fCenFlattenHist->GetBinContent(
          fCenFlattenHist->FindBin(fCentrality[fCentralityEstimator]));
      if (gRandom->Uniform() > CenProb) {
        fEventCutsCounter[kRECO]->Fill(2 * (LAST_EEVENT + 2) + 0.5);
        // CutBit += TMath::Power(2, LAST_EEVENT + 2);
        Flag = kFALSE;
      }
    }
    // fEventCutsCounterCumulative->Fill(&CutBit);
  }
  return Flag;
}

Bool_t AliAnalysisTaskAR::SurviveTrackCut(AliVParticle *avp,
                                          Bool_t FillCounter) {
  // check if current track survives track cut
  // return flag at the end, if one cut fails, set it to false
  Bool_t Flag = kTRUE;
  // Double_t CutBit = 0;

  // only cut if we get valid pointer
  AliAODTrack *aTrack = dynamic_cast<AliAODTrack *>(avp);
  if (aTrack) {

    if (fUseTrackCuts[kPT]) {
      // cut PT
      if (aTrack->Pt() < fTrackCuts[kPT][kMIN]) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * kPT + kMIN + 0.5);
        }
        // CutBit += TMath::Power(2, kPT);
        Flag = kFALSE;
      }
      if (aTrack->Pt() > fTrackCuts[kPT][kMAX]) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * kPT + kMAX + 0.5);
        }
        // CutBit += TMath::Power(2, kPT);
        Flag = kFALSE;
      }
    }
    if (fUseTrackCuts[kETA]) {
      // cut ETA
      if (aTrack->Eta() < fTrackCuts[kETA][kMIN]) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * kETA + kMIN + 0.5);
        }
        // CutBit += TMath::Power(2, kETA);
        Flag = kFALSE;
      }
      if (aTrack->Eta() > fTrackCuts[kETA][kMAX]) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * kETA + kMAX + 0.5);
        }
        // CutBit += TMath::Power(2, kETA);
        Flag = kFALSE;
      }
    }
    if (fUseTrackCuts[kPHI]) {
      // cut PHI
      if (aTrack->Phi() < fTrackCuts[kPHI][kMIN]) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * kPHI + kMIN + 0.5);
        }
        // CutBit += TMath::Power(2, kPHI);
        Flag = kFALSE;
      }
      if (aTrack->Phi() > fTrackCuts[kPHI][kMAX]) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * kPHI + kMAX + 0.5);
        }
        // CutBit += TMath::Power(2, kPHI);
        Flag = kFALSE;
      }
    }
    if (fUseTrackCuts[kCHARGE]) {
      // cut on CHARGE
      if (aTrack->Charge() < fTrackCuts[kCHARGE][kMIN]) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * kCHARGE + kMIN + 0.5);
        }
        // CutBit += TMath::Power(2, kCHARGE);
        Flag = kFALSE;
      }
      if (aTrack->Charge() > fTrackCuts[kCHARGE][kMAX]) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * kCHARGE + kMAX + 0.5);
        }
        // CutBit += TMath::Power(2, kCHARGE);
        Flag = kFALSE;
      }
    }
    if (fUseTrackCuts[kTPCNCLS]) {
      // cut on number of clusters in the TPC
      if (aTrack->GetTPCNcls() < fTrackCuts[kTPCNCLS][kMIN]) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * kTPCNCLS + kMIN + 0.5);
        }
        // CutBit += TMath::Power(2, kTPCNCLS);
        Flag = kFALSE;
      }
      if (aTrack->GetTPCNcls() > fTrackCuts[kTPCNCLS][kMAX]) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * kTPCNCLS + kMAX + 0.5);
        }
        // CutBit += TMath::Power(2, kTPCNCLS);
        Flag = kFALSE;
      }
    }
    if (fUseTrackCuts[kTPCCROSSEDROWS]) {
      // cut on crossed rows in the TPC
      if (aTrack->GetTPCNCrossedRows() < fTrackCuts[kTPCCROSSEDROWS][kMIN]) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * kTPCCROSSEDROWS + kMIN + 0.5);
        }
        // CutBit += TMath::Power(2, kTPCCROSSEDROWS);
        Flag = kFALSE;
      }
      if (aTrack->GetTPCNCrossedRows() > fTrackCuts[kTPCCROSSEDROWS][kMAX]) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * kTPCCROSSEDROWS + kMAX + 0.5);
        }
        // CutBit += TMath::Power(2, kTPCCROSSEDROWS);
        Flag = kFALSE;
      }
    }

    if (fUseTrackCuts[kTPCNCLSFRACTIONSHARED]) {
      // cut on ratio of shared clusters in the TPC
      Double_t tpcnclsfractionshared;
      if (aTrack->GetTPCNcls() != 0) {
        tpcnclsfractionshared =
            (Double_t)aTrack->GetTPCnclsS() / (Double_t)aTrack->GetTPCNcls();
      } else {
        tpcnclsfractionshared = 1.;
      }
      if (tpcnclsfractionshared < fTrackCuts[kTPCNCLSFRACTIONSHARED][kMIN]) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * kTPCNCLSFRACTIONSHARED + kMIN +
                                         0.5);
        }
        // CutBit += TMath::Power(2, kTPCNCLSFRACTIONSHARED);
        Flag = kFALSE;
      }
      if (tpcnclsfractionshared > fTrackCuts[kTPCNCLSFRACTIONSHARED][kMAX]) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * kTPCNCLSFRACTIONSHARED + kMAX +
                                         0.5);
        }
        // CutBit += TMath::Power(2, kTPCNCLSFRACTIONSHARED);
        Flag = kFALSE;
      }
    }
    if (fUseTrackCuts[kTPCCHI2PERNDF]) {
      // cut on chi^2/NDF of the tracks in the TPC
      Double_t tpcchi2perndf;
      if (aTrack->GetTPCNcls() != 0) {
        tpcchi2perndf =
            (Double_t)aTrack->GetTPCchi2() / (Double_t)aTrack->GetTPCNcls();
      } else {
        tpcchi2perndf = 5;
      }
      if (tpcchi2perndf < fTrackCuts[kTPCCHI2PERNDF][kMIN]) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * kTPCCHI2PERNDF + kMIN + 0.5);
        }
        // CutBit += TMath::Power(2, kTPCCHI2PERNDF);
        Flag = kFALSE;
      }
      if (tpcchi2perndf > fTrackCuts[kTPCCHI2PERNDF][kMAX]) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * kTPCCHI2PERNDF + kMAX + 0.5);
        }
        // CutBit += TMath::Power(2, kTPCCHI2PERNDF);
        Flag = kFALSE;
      }
    }
    if (fUseTrackCuts[kITSNCLS]) {
      // cut on number of clusters in the ITS
      if (aTrack->GetITSNcls() < fTrackCuts[kITSNCLS][kMIN]) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * kITSNCLS + kMIN + 0.5);
        }
        // CutBit += TMath::Power(2, kITSNCLS);
        Flag = kFALSE;
      }
      if (aTrack->GetITSNcls() > fTrackCuts[kITSNCLS][kMAX]) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * kITSNCLS + kMAX + 0.5);
        }
        // CutBit += TMath::Power(2, kITSNCLS);
        Flag = kFALSE;
      }
    }
    if (fUseTrackCuts[kCHI2PERNDF]) {
      // cut on chi2 / NDF of the track fit
      if (aTrack->Chi2perNDF() < fTrackCuts[kCHI2PERNDF][kMIN]) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * kCHI2PERNDF + kMIN + 0.5);
        }
        // CutBit += TMath::Power(2, kCHI2PERNDF);
        Flag = kFALSE;
      }
      if (aTrack->Chi2perNDF() > fTrackCuts[kCHI2PERNDF][kMAX]) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * kCHI2PERNDF + kMAX + 0.5);
        }
        // CutBit += TMath::Power(2, kCHI2PERNDF);
        Flag = kFALSE;
      }
    }
    if (fUseTrackCuts[kDCAZ]) {
      // cut DCA in z direction
      if (aTrack->ZAtDCA() < fTrackCuts[kDCAZ][kMIN]) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * kDCAZ + kMIN + 0.5);
          // if track is not constrained it returns dummy value -999
          // makes the counter blow up
        }
        // CutBit += TMath::Power(2, kDCAZ);
        Flag = kFALSE;
      }
      if (aTrack->ZAtDCA() > fTrackCuts[kDCAZ][kMAX]) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * kDCAZ + kMAX + 0.5);
          // if track is not constrained it returns dummy value -999
          // makes the counter blow up
        }
        // CutBit += TMath::Power(2, kDCAZ);
        Flag = kFALSE;
      }
    }
    if (fUseTrackCuts[kDCAXY]) {
      // cut DCA in xy plane
      if (aTrack->DCA() < fTrackCuts[kDCAXY][kMIN]) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * kDCAXY + kMIN + 0.5);
        }
        // CutBit += TMath::Power(2, kDCAXY);
        Flag = kFALSE;
      }
      if (aTrack->DCA() > fTrackCuts[kDCAXY][kMAX]) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * kDCAXY + kMAX + 0.5);
        }
        // CutBit += TMath::Power(2, kDCAXY);
        Flag = kFALSE;
      }
    }
    if (fUseFilterbit) {
      // cut with filtertbit
      // filter bit 128 denotes TPC-only tracks, use only them for the
      // analysis, for hybrid tracks use filterbit 782
      // for more information about filterbits see the online wiki
      // the filterbits can change from run to run
      if (!aTrack->TestFilterBit(fFilterbit)) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * LAST_ETRACK + 0.5);
        }
        // CutBit += TMath::Power(2, LAST_ETRACK + 0);
        Flag = kFALSE;
      }
    }
    // if set, cut all neutral tracks away
    if (fChargedOnly) {
      if (aTrack->Charge() == 0) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * LAST_ETRACK + 1.5);
        }
        // CutBit += TMath::Power(2, LAST_ETRACK + 1);
        Flag = kFALSE;
      }
    }

    // if set, cut all non-primary tracks away
    if (fPrimaryOnly) {
      if (aTrack->GetType() != AliAODTrack::kPrimary) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * LAST_ETRACK + 2.5);
        }
        // CutBit += TMath::Power(2, LAST_ETRACK + 2);
        Flag = kFALSE;
      }
    }

    // if track id is negative, it is not a global track
    if (fGlobalTracksOnly) {
      if (aTrack->GetID() < 0) {
        if (FillCounter) {
          fTrackCutsCounter[kRECO]->Fill(2 * LAST_ETRACK + 3.5);
        }
        // CutBit += TMath::Power(2, LAST_ETRACK + 3);
        Flag = kFALSE;
      }
    }
    // fTrackCutsCounterCumulative->Fill(&CutBit);

    // when doing monte carlo closure, check if we accept the track
    // also check if we will the counter, i.e. we run in the main loop
    // we do not want to loose events because we get unlucky here
    if (fMCClosure && FillCounter) {
      if (fAcceptanceHistogram[kPT]) {
        if (fAcceptanceHistogram[kPT]->GetBinContent(
                fAcceptanceHistogram[kPT]->FindBin(aTrack->Pt())) <
            gRandom->Uniform()) {
          if (FillCounter) {
            fTrackCutsCounter[kRECO]->Fill(2 * LAST_ETRACK + 4.5);
          }
          Flag = kFALSE;
        }
      }
      if (fAcceptanceHistogram[kETA]) {
        if (fAcceptanceHistogram[kETA]->GetBinContent(
                fAcceptanceHistogram[kETA]->FindBin(aTrack->Eta())) <
            gRandom->Uniform()) {
          if (FillCounter) {
            fTrackCutsCounter[kRECO]->Fill(2 * LAST_ETRACK + 4.5);
          }
          Flag = kFALSE;
        }
      }
      if (fAcceptanceHistogram[kPHI]) {
        if (fAcceptanceHistogram[kPHI]->GetBinContent(
                fAcceptanceHistogram[kPHI]->FindBin(aTrack->Phi())) <
            gRandom->Uniform()) {
          if (FillCounter) {
            fTrackCutsCounter[kRECO]->Fill(2 * LAST_ETRACK + 4.5);
          }
          Flag = kFALSE;
        }
      }
    }
  }

  // check MC particle
  AliAODMCParticle *MCParticle = dynamic_cast<AliAODMCParticle *>(avp);
  if (MCParticle) {
    if (fUseTrackCuts[kPT]) {
      // cut PT
      if (MCParticle->Pt() < fTrackCuts[kPT][kMIN]) {
        if (FillCounter) {
          fTrackCutsCounter[kSIM]->Fill(2 * kPT + kMIN + 0.5);
        }
        Flag = kFALSE;
      }
      if (MCParticle->Pt() > fTrackCuts[kPT][kMAX]) {
        if (FillCounter) {
          fTrackCutsCounter[kSIM]->Fill(2 * kPT + kMAX + 0.5);
        }
        Flag = kFALSE;
      }
    }
    if (fUseTrackCuts[kETA]) {
      // cut ETA
      if (MCParticle->Eta() < fTrackCuts[kETA][kMIN]) {
        if (FillCounter) {
          fTrackCutsCounter[kSIM]->Fill(2 * kETA + kMIN + 0.5);
        }
        Flag = kFALSE;
      }
      if (MCParticle->Eta() > fTrackCuts[kETA][kMAX]) {
        if (FillCounter) {
          fTrackCutsCounter[kSIM]->Fill(2 * kETA + kMAX + 0.5);
        }
        Flag = kFALSE;
      }
    }
    if (fUseTrackCuts[kPHI]) {
      // cut PHI
      if (MCParticle->Phi() < fTrackCuts[kPHI][kMIN]) {
        if (FillCounter) {
          fTrackCutsCounter[kSIM]->Fill(2 * kPHI + kMIN + 0.5);
        }
        Flag = kFALSE;
      }
      if (MCParticle->Phi() > fTrackCuts[kPHI][kMAX]) {
        if (FillCounter) {
          fTrackCutsCounter[kSIM]->Fill(2 * kPHI + kMAX + 0.5);
        }
        Flag = kFALSE;
      }
    }
    if (fUseTrackCuts[kCHARGE]) {
      // cut CHARGE
      if ((MCParticle->Charge() / 3.) < fTrackCuts[kCHARGE][kMIN]) {
        if (FillCounter) {
          fTrackCutsCounter[kSIM]->Fill(2 * kCHARGE + kMIN + 0.5);
        }
        Flag = kFALSE;
      }
      if ((MCParticle->Charge() / 3.) > fTrackCuts[kCHARGE][kMAX]) {
        if (FillCounter) {
          fTrackCutsCounter[kSIM]->Fill(2 * kCHARGE + kMAX + 0.5);
        }
        Flag = kFALSE;
      }
    }
    // if set, cut all neutral particles away
    if (fChargedOnly) {
      if (MCParticle->Charge() == 0) {
        if (FillCounter) {
          fTrackCutsCounter[kSIM]->Fill(2 * LAST_ETRACK + 1.5);
        }
        Flag = kFALSE;
      }
    }
    // if set, cut all non-primary particles away
    if (fPrimaryOnly) {
      if (fMCPrimaryDef == kMCPrim && !MCParticle->IsPrimary()) {
        fTrackCutsCounter[kSIM]->Fill(2 * LAST_ETRACK + 2.5);
        Flag = kFALSE;
      }
    } else if (fMCPrimaryDef == kMCPhysicalPrim &&
               !MCParticle->IsPhysicalPrimary()) {

      fTrackCutsCounter[kSIM]->Fill(2 * LAST_ETRACK + 2.5);
      Flag = kFALSE;
    }
  }

  // if both aTrack and MCParticle are null, check if we generated monte carlo
  // data on the fly
  if (!aTrack && !MCParticle && fMCOnTheFly) {

    if (fAcceptanceHistogram[kPT]) {
      if (fAcceptanceHistogram[kPT]->GetBinContent(
              fAcceptanceHistogram[kPT]->FindBin(fMCKinematicVariables[kPT])) <
          gRandom->Uniform()) {
        Flag = kFALSE;
      }
    }
    if (fAcceptanceHistogram[kETA]) {
      if (fAcceptanceHistogram[kETA]->GetBinContent(
              fAcceptanceHistogram[kETA]->FindBin(
                  fMCKinematicVariables[kETA])) < gRandom->Uniform()) {
        Flag = kFALSE;
      }
    }
    if (fAcceptanceHistogram[kPHI]) {
      if (fAcceptanceHistogram[kPHI]->GetBinContent(
              fAcceptanceHistogram[kPHI]->FindBin(
                  fMCKinematicVariables[kPHI])) < gRandom->Uniform()) {
        Flag = kFALSE;
      }
    }
  }
  return Flag;
}

void AliAnalysisTaskAR::FillEventObjects(AliAODEvent *aAOD, AliMCEvent *aMC) {
  // get/compute event variables and fill them into data members

  // if we get two null pointers, we are generating monte carlo data on our
  // own
  if (!aAOD && !aMC) {
    return;
  }

  // get centralities
  AliMultSelection *aMS =
      dynamic_cast<AliMultSelection *>(aAOD->FindListObject("MultSelection"));
  for (int cen = 0; cen < LAST_ECENESTIMATORS; ++cen) {
    fCentrality[cen] = aMS->GetMultiplicityPercentile(kCenEstimatorNames[cen]);
  }

  // multiplicity as number of tracks
  fMultiplicity[kMUL] = 0.;

  // multiplicity as number of contributors to the primary vertex
  AliAODVertex *PrimaryVertex = aAOD->GetPrimaryVertex();
  fMultiplicity[kNCONTRIB] = PrimaryVertex->GetNContributors();

  // reference multiplicity
  // combined reference multiplicity (tracklets + ITSTPC) in |eta|<0.8
  AliAODHeader *Header = dynamic_cast<AliAODHeader *>(aAOD->GetHeader());
  fMultiplicity[kMULREF] = Header->GetRefMultiplicityComb08();

  // multiplicity as number of tracks that survive track cuts
  fMultiplicity[kMULQ] = 0.;
  // multiplicity as the weighted sum of all surviving tracks
  fMultiplicity[kMULW] = 0.;
  Double_t w = 1.;

  // get number of tracks in the event
  Int_t nTracks = aAOD->GetNumberOfTracks();

  // clear randomized indices
  if (fUseFisherYates) {
    fRandomizedTrackIndices.clear();
    for (int i = 0; i < nTracks; i++) {
      fRandomizedTrackIndices.push_back(i);
    }
  }
  Int_t jTrack = 0;
  Int_t tmp = 0;

  // reset lookuptable
  fLookUpTable.clear();

  AliAODTrack *aTrack = nullptr;

  for (Int_t iTrack = 0; iTrack < nTracks; ++iTrack) {

    // compute randomized index if necessary using Fisher-Yates shuffel
    // pseudo-code taken from
    // https://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle#Examples
    // To shuffle an array a of n elements (indices 0..n-1):
    // for i from 0 to n−2 do
    // j ← random integer such that i ≤ j < n
    // exchange a[i] and a[j]
    if (fUseFisherYates && iTrack < nTracks - 1) {
      jTrack = iTrack + gRandom->Integer(nTracks - iTrack);
      tmp = fRandomizedTrackIndices.at(iTrack);
      fRandomizedTrackIndices.at(iTrack) = fRandomizedTrackIndices.at(jTrack);
      fRandomizedTrackIndices.at(jTrack) = tmp;
    }

    // getting pointer to a track
    aTrack = dynamic_cast<AliAODTrack *>(aAOD->GetTrack(iTrack));

    // protect against invalid pointers
    if (!aTrack) {
      continue;
    }

    if (aTrack->GetID() >= 0) {
      fMultiplicity[kMUL] += 1.;
    }

    if (!SurviveTrackCut(aTrack, kFALSE)) {
      continue;
    }
    fMultiplicity[kMULQ] += 1.;

    w = 1.;
    if (fUseWeights) {
      if (fWeightHistogram[kPT]) {
        w *= fWeightHistogram[kPT]->GetBinContent(
            fWeightHistogram[kPT]->FindBin(aTrack->Pt()));
      }
      if (fWeightHistogram[kPHI]) {
        w *= fWeightHistogram[kPHI]->GetBinContent(
            fWeightHistogram[kPHI]->FindBin(aTrack->Phi()));
      }
      if (fWeightHistogram[kETA]) {
        w *= fWeightHistogram[kETA]->GetBinContent(
            fWeightHistogram[kETA]->FindBin(aTrack->Eta()));
      }
    }
    fMultiplicity[kMULW] += w;

    // since we are already looping over the events, create a look up table if
    // we also have a monte carlo event
    if (aMC) {
      // "key" = label, "value" = iTrack
      fLookUpTable.insert({aTrack->GetLabel(), iTrack});
    }
  }
}

void AliAnalysisTaskAR::CalculateQvectors(std::vector<Double_t> angles,
                                          std::vector<Double_t> weights) {
  // Calculate all Q-vectors

  // Make sure all Q-vectors are initially zero
  for (Int_t h = 0; h < kMaxHarmonic; h++) {
    for (Int_t p = 0; p < kMaxPower; p++) {
      fQvector[h][p] = TComplex(0., 0.);
    }
  }

  // Calculate Q-vectors for available angles and weights
  Double_t dPhi = 0.;
  Double_t wPhi = 1.;         // particle weight
  Double_t wPhiToPowerP = 1.; // particle weight raised to power p
  for (std::size_t i = 0; i < angles.size(); i++) {
    dPhi = angles.at(i);
    wPhi = weights.at(i);
    for (Int_t h = 0; h < kMaxHarmonic; h++) {
      for (Int_t p = 0; p < kMaxPower; p++) {
        wPhiToPowerP = pow(wPhi, p);
        fQvector[h][p] += TComplex(wPhiToPowerP * TMath::Cos(h * dPhi),
                                   wPhiToPowerP * TMath::Sin(h * dPhi));
      }
    }
  }
}

void AliAnalysisTaskAR::FillFinalResultCorrelators() {
  // fill final result profiles

  Double_t corr = 0.0;
  Double_t weight = 1.0;

  for (Int_t kin = 0; kin < kKinematic; kin++) {
    for (std::size_t Bin = 0; Bin < fTrackBins[kin].size() - 1; Bin++) {
      CalculateQvectors(fKinematics[kin].at(Bin),
                        fKinematicWeights[kin].at(Bin));
      corr = 0.;
      weight = 1.;

      // loop over all correlators
      for (std::size_t i = 0; i < fCorrelators.size(); i++) {
        // protect against insufficient amount of statistics i.e. number of
        // particles is lower then the order of correlator due to track cuts
        if (fKinematics[kin].at(Bin).size() <= fCorrelators.at(i).size()) {
          std::cout
              << "Not enough tracks in this event to compute the correlator"
              << std::endl
              << "Number of Tracks: " << fKinematics[kin].at(Bin).size()
              << std::endl
              << "Correlator: v_{";
          for (auto e : fCorrelators.at(i)) {
            std ::cout << e << ",";
          }
          std::cout << "}" << std::endl
                    << "TrackBin: " << kin << std::endl
                    << "Bin: " << Bin << std::endl;
          continue;
        }

        // compute correlator
        if (fUseNestedLoops) {
          // using nested loops
          switch (static_cast<int>(fCorrelators.at(i).size())) {
          case 2:
            corr = TwoNestedLoops(
                       fCorrelators.at(i).at(0), fCorrelators.at(i).at(1),
                       fKinematics[kin].at(Bin), fKinematicWeights[kin].at(Bin))
                       .Re();
            weight = TwoNestedLoops(0, 0, fKinematics[kin].at(Bin),
                                    fKinematicWeights[kin].at(Bin))
                         .Re();
            break;
          case 3:
            corr = ThreeNestedLoops(
                       fCorrelators.at(i).at(0), fCorrelators.at(i).at(1),
                       fCorrelators.at(i).at(2), fKinematics[kin].at(Bin),
                       fKinematicWeights[kin].at(Bin))
                       .Re();
            weight = ThreeNestedLoops(0, 0, 0, fKinematics[kin].at(Bin),
                                      fKinematicWeights[kin].at(Bin))
                         .Re();
            break;
          case 4:
            corr = FourNestedLoops(
                       fCorrelators.at(i).at(0), fCorrelators.at(i).at(1),
                       fCorrelators.at(i).at(2), fCorrelators.at(i).at(3),
                       fKinematics[kin].at(Bin), fKinematicWeights[kin].at(Bin))
                       .Re();
            weight = FourNestedLoops(0, 0, 0, 0, fKinematics[kin].at(Bin),
                                     fKinematicWeights[kin].at(Bin))
                         .Re();
            break;
          case 5:
            corr = FiveNestedLoops(
                       fCorrelators.at(i).at(0), fCorrelators.at(i).at(1),
                       fCorrelators.at(i).at(2), fCorrelators.at(i).at(3),
                       fCorrelators.at(i).at(4), fKinematics[kin].at(Bin),
                       fKinematicWeights[kin].at(Bin))
                       .Re();
            weight = FiveNestedLoops(0, 0, 0, 0, 0, fKinematics[kin].at(Bin),
                                     fKinematicWeights[kin].at(Bin))
                         .Re();
            break;
          case 6:
            corr = SixNestedLoops(
                       fCorrelators.at(i).at(0), fCorrelators.at(i).at(1),
                       fCorrelators.at(i).at(2), fCorrelators.at(i).at(3),
                       fCorrelators.at(i).at(4), fCorrelators.at(i).at(5),
                       fKinematics[kin].at(Bin), fKinematicWeights[kin].at(Bin))
                       .Re();
            weight = SixNestedLoops(0, 0, 0, 0, 0, 0, fKinematics[kin].at(Bin),
                                    fKinematicWeights[kin].at(Bin))
                         .Re();
            break;
          default:
            std::cout << "Correlators of order >6 are not implemented with "
                         "nested loops"
                      << std::endl;
          }
        } else {
          // using Qvectors
          switch (static_cast<int>(fCorrelators.at(i).size())) {
          case 2:
            corr = Two(fCorrelators.at(i).at(0), fCorrelators.at(i).at(1)).Re();
            weight = Two(0, 0).Re();
            break;
          case 3:
            corr = Three(fCorrelators.at(i).at(0), fCorrelators.at(i).at(1),
                         fCorrelators.at(i).at(2))
                       .Re();
            weight = Three(0, 0, 0).Re();
            break;
          case 4:
            corr = Four(fCorrelators.at(i).at(0), fCorrelators.at(i).at(1),
                        fCorrelators.at(i).at(2), fCorrelators.at(i).at(3))
                       .Re();
            weight = Four(0, 0, 0, 0).Re();
            break;
          case 5:
            corr = Five(fCorrelators.at(i).at(0), fCorrelators.at(i).at(1),
                        fCorrelators.at(i).at(2), fCorrelators.at(i).at(3),
                        fCorrelators.at(i).at(4))
                       .Re();
            weight = Five(0, 0, 0, 0, 0).Re();
            break;
          case 6:
            corr = Six(fCorrelators.at(i).at(0), fCorrelators.at(i).at(1),
                       fCorrelators.at(i).at(2), fCorrelators.at(i).at(3),
                       fCorrelators.at(i).at(4), fCorrelators.at(i).at(5))
                       .Re();
            weight = Six(0, 0, 0, 0, 0, 0).Re();
            break;
          default:
            corr =
                Recursion(fCorrelators.at(i).size(), fCorrelators.at(i).data())
                    .Re();
            weight =
                Recursion(
                    fCorrelators.at(i).size(),
                    std::vector<Int_t>(fCorrelators.at(i).size(), 0).data())
                    .Re();
          }
        }

        // correlators are not normalized yet
        corr /= weight;

        // fill final result profiles
        // integrated correlator
        if (kin == kPT) {
          dynamic_cast<TProfile *>(
              dynamic_cast<TList *>(fFinalResultCorrelatorsList->At(i))
                  ->At(kPTDEP))
              ->Fill(fTrackBinsHistogram[kPT]->GetBinCenter(Bin + 1), corr,
                     weight);
        } else if (kin == kETA) {
          dynamic_cast<TProfile *>(
              dynamic_cast<TList *>(fFinalResultCorrelatorsList->At(i))
                  ->At(kETADEP))
              ->Fill(fTrackBinsHistogram[kETA]->GetBinCenter(Bin + 1), corr,
                     weight);
        } else {
          dynamic_cast<TProfile *>(
              dynamic_cast<TList *>(fFinalResultCorrelatorsList->At(i))
                  ->At(kINTEGRATED))
              ->Fill(0.5, corr, weight);
          // correlator as a function of centrality
          dynamic_cast<TProfile *>(
              dynamic_cast<TList *>(fFinalResultCorrelatorsList->At(i))
                  ->At(kCENDEP))
              ->Fill(fCentrality[fCentralityEstimator], corr, weight);
          // correlator as a function of multiplicity
          dynamic_cast<TProfile *>(
              dynamic_cast<TList *>(fFinalResultCorrelatorsList->At(i))
                  ->At(kMULDEP))
              ->Fill(fMultiplicity[kMULQ], corr, weight);
        }
      }
    }
  }
}

void AliAnalysisTaskAR::FillSymmetricCumulant() {
  // fill symmetric cumulant

  for (std::size_t i = 0; i < fSymmetricCumulants.size(); i++) {

    switch (fSymmetricCumulants.at(i).size()) {
    case 2:
      SC2(fSymmetricCumulants.at(i), i);
      break;
    case 3:
      SC3(fSymmetricCumulants.at(i), i);
      break;
    default:
      std::cout
          << "Symmetric cumulants of order 4 and higher are not implemented yet"
          << std::endl;
      break;
    }
  }
}

void AliAnalysisTaskAR::SC2(std::vector<Int_t> sc, Int_t index) {

  TList *listSC_kl =
      dynamic_cast<TList *>(fFinalResultSymmetricCumulantsList->At(index));
  TList *listNSC_kl = dynamic_cast<TList *>(
      fFinalResultNormalizedSymmetricCumulantsList->At(index));

  std::vector<std::vector<Int_t>> correlators = fMapSCtoCor.at(sc);

  TList *listC_kl = dynamic_cast<TList *>(
      fFinalResultCorrelatorsList->At(fMapCorToIndex.at(correlators.at(0))));
  TList *listC_k = dynamic_cast<TList *>(
      fFinalResultCorrelatorsList->At(fMapCorToIndex.at(correlators.at(1))));
  TList *listC_l = dynamic_cast<TList *>(
      fFinalResultCorrelatorsList->At(fMapCorToIndex.at(correlators.at(2))));

  TH1D *sc_kl, *nsc_kl;
  TProfile *c_kl, *c_k, *c_l;
  Double_t sc_value, norm;

  for (Int_t i = 0; i < LAST_EFINALRESULTPROFILE; i++) {

    sc_kl = dynamic_cast<TH1D *>(listSC_kl->At(i));
    nsc_kl = dynamic_cast<TH1D *>(listNSC_kl->At(i));
    c_kl = dynamic_cast<TProfile *>(listC_kl->At(i));
    c_k = dynamic_cast<TProfile *>(listC_k->At(i));
    c_l = dynamic_cast<TProfile *>(listC_l->At(i));

    for (Int_t bin = 1; bin <= sc_kl->GetNbinsX(); bin++) {
      sc_value = c_kl->GetBinContent(bin) -
                 c_k->GetBinContent(bin) * c_l->GetBinContent(bin);
      sc_kl->SetBinContent(bin, sc_value);

      norm = c_k->GetBinContent(bin) * c_l->GetBinContent(bin);

      if (std::abs(norm) > 1e-25) {
        nsc_kl->SetBinContent(bin, sc_value / norm);
      }
    }
  }
}

void AliAnalysisTaskAR::SC3(std::vector<Int_t> sc, Int_t index) {

  TList *listSC_kln =
      dynamic_cast<TList *>(fFinalResultSymmetricCumulantsList->At(index));
  TList *listNSC_kln = dynamic_cast<TList *>(
      fFinalResultNormalizedSymmetricCumulantsList->At(index));
  std::vector<std::vector<Int_t>> correlators = fMapSCtoCor.at(sc);

  TList *listC_kln = dynamic_cast<TList *>(
      fFinalResultCorrelatorsList->At(fMapCorToIndex.at(correlators.at(0))));
  TList *listC_kl = dynamic_cast<TList *>(
      fFinalResultCorrelatorsList->At(fMapCorToIndex.at(correlators.at(1))));
  TList *listC_kn = dynamic_cast<TList *>(
      fFinalResultCorrelatorsList->At(fMapCorToIndex.at(correlators.at(2))));
  TList *listC_ln = dynamic_cast<TList *>(
      fFinalResultCorrelatorsList->At(fMapCorToIndex.at(correlators.at(3))));
  TList *listC_k = dynamic_cast<TList *>(
      fFinalResultCorrelatorsList->At(fMapCorToIndex.at(correlators.at(4))));
  TList *listC_l = dynamic_cast<TList *>(
      fFinalResultCorrelatorsList->At(fMapCorToIndex.at(correlators.at(5))));
  TList *listC_n = dynamic_cast<TList *>(
      fFinalResultCorrelatorsList->At(fMapCorToIndex.at(correlators.at(6))));

  TH1D *sc_kln, *nsc_kln;
  TProfile *c_kln, *c_kl, *c_kn, *c_ln, *c_k, *c_l, *c_n;
  Double_t sc_value, norm;

  for (Int_t i = 0; i < LAST_EFINALRESULTPROFILE; i++) {

    sc_kln = dynamic_cast<TH1D *>(listSC_kln->At(i));
    nsc_kln = dynamic_cast<TH1D *>(listNSC_kln->At(i));
    c_kln = dynamic_cast<TProfile *>(listC_kln->At(i));
    c_kl = dynamic_cast<TProfile *>(listC_kl->At(i));
    c_kn = dynamic_cast<TProfile *>(listC_kn->At(i));
    c_ln = dynamic_cast<TProfile *>(listC_ln->At(i));
    c_k = dynamic_cast<TProfile *>(listC_k->At(i));
    c_l = dynamic_cast<TProfile *>(listC_l->At(i));
    c_n = dynamic_cast<TProfile *>(listC_n->At(i));

    for (Int_t bin = 1; bin <= sc_kln->GetNbinsX(); bin++) {
      sc_value = c_kln->GetBinContent(bin) -
                 c_kl->GetBinContent(bin) * c_n->GetBinContent(bin) -
                 c_kn->GetBinContent(bin) * c_l->GetBinContent(bin) -
                 c_ln->GetBinContent(bin) * c_k->GetBinContent(bin) +
                 2 * c_k->GetBinContent(bin) * c_l->GetBinContent(bin) *
                     c_n->GetBinContent(bin);
      sc_kln->SetBinContent(bin, sc_value);

      norm = c_k->GetBinContent(bin) * c_l->GetBinContent(bin) *
             c_n->GetBinContent(bin);

      if (std::abs(norm) > 1e-25) {
        nsc_kln->SetBinContent(bin, sc_value / norm);
      }
    }
  }
}

TComplex AliAnalysisTaskAR::Q(Int_t n, Int_t p) {
  // return Qvector from fQvector array

  if (n > kMaxHarmonic || p > kMaxPower) {
    std::cout << __LINE__ << ": running out of bounds" << std::endl;
    Fatal("Q", "Running out of bounds in fQvector");
  }
  if (n >= 0) {
    return fQvector[n][p];
  }
  return TComplex::Conjugate(fQvector[-n][p]);
}

TComplex AliAnalysisTaskAR::Two(Int_t n1, Int_t n2) {
  // Generic two-particle correlation <exp[i(n1*phi1+n2*phi2)]>.
  TComplex two = Q(n1, 1) * Q(n2, 1) - Q(n1 + n2, 2);
  return two;
}

TComplex AliAnalysisTaskAR::Three(Int_t n1, Int_t n2, Int_t n3) {
  // Generic three-particle correlation <exp[i(n1*phi1+n2*phi2+n3*phi3)]>.
  TComplex three = Q(n1, 1) * Q(n2, 1) * Q(n3, 1) - Q(n1 + n2, 2) * Q(n3, 1) -
                   Q(n2, 1) * Q(n1 + n3, 2) - Q(n1, 1) * Q(n2 + n3, 2) +
                   2. * Q(n1 + n2 + n3, 3);
  return three;
}

TComplex AliAnalysisTaskAR::Four(Int_t n1, Int_t n2, Int_t n3, Int_t n4) {
  // Generic four-particle correlation
  // <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4)]>.
  TComplex four =
      Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4, 1) -
      Q(n1 + n2, 2) * Q(n3, 1) * Q(n4, 1) -
      Q(n2, 1) * Q(n1 + n3, 2) * Q(n4, 1) -
      Q(n1, 1) * Q(n2 + n3, 2) * Q(n4, 1) + 2. * Q(n1 + n2 + n3, 3) * Q(n4, 1) -
      Q(n2, 1) * Q(n3, 1) * Q(n1 + n4, 2) + Q(n2 + n3, 2) * Q(n1 + n4, 2) -
      Q(n1, 1) * Q(n3, 1) * Q(n2 + n4, 2) + Q(n1 + n3, 2) * Q(n2 + n4, 2) +
      2. * Q(n3, 1) * Q(n1 + n2 + n4, 3) - Q(n1, 1) * Q(n2, 1) * Q(n3 + n4, 2) +
      Q(n1 + n2, 2) * Q(n3 + n4, 2) + 2. * Q(n2, 1) * Q(n1 + n3 + n4, 3) +
      2. * Q(n1, 1) * Q(n2 + n3 + n4, 3) - 6. * Q(n1 + n2 + n3 + n4, 4);

  return four;
}

TComplex AliAnalysisTaskAR::Five(Int_t n1, Int_t n2, Int_t n3, Int_t n4,
                                 Int_t n5) {
  // Generic five-particle correlation
  // <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5)]>.
  TComplex five = Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) -
                  Q(n1 + n2, 2) * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) -
                  Q(n2, 1) * Q(n1 + n3, 2) * Q(n4, 1) * Q(n5, 1) -
                  Q(n1, 1) * Q(n2 + n3, 2) * Q(n4, 1) * Q(n5, 1) +
                  2. * Q(n1 + n2 + n3, 3) * Q(n4, 1) * Q(n5, 1) -
                  Q(n2, 1) * Q(n3, 1) * Q(n1 + n4, 2) * Q(n5, 1) +
                  Q(n2 + n3, 2) * Q(n1 + n4, 2) * Q(n5, 1) -
                  Q(n1, 1) * Q(n3, 1) * Q(n2 + n4, 2) * Q(n5, 1) +
                  Q(n1 + n3, 2) * Q(n2 + n4, 2) * Q(n5, 1) +
                  2. * Q(n3, 1) * Q(n1 + n2 + n4, 3) * Q(n5, 1) -
                  Q(n1, 1) * Q(n2, 1) * Q(n3 + n4, 2) * Q(n5, 1) +
                  Q(n1 + n2, 2) * Q(n3 + n4, 2) * Q(n5, 1) +
                  2. * Q(n2, 1) * Q(n1 + n3 + n4, 3) * Q(n5, 1) +
                  2. * Q(n1, 1) * Q(n2 + n3 + n4, 3) * Q(n5, 1) -
                  6. * Q(n1 + n2 + n3 + n4, 4) * Q(n5, 1) -
                  Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n1 + n5, 2) +
                  Q(n2 + n3, 2) * Q(n4, 1) * Q(n1 + n5, 2) +
                  Q(n3, 1) * Q(n2 + n4, 2) * Q(n1 + n5, 2) +
                  Q(n2, 1) * Q(n3 + n4, 2) * Q(n1 + n5, 2) -
                  2. * Q(n2 + n3 + n4, 3) * Q(n1 + n5, 2) -
                  Q(n1, 1) * Q(n3, 1) * Q(n4, 1) * Q(n2 + n5, 2) +
                  Q(n1 + n3, 2) * Q(n4, 1) * Q(n2 + n5, 2) +
                  Q(n3, 1) * Q(n1 + n4, 2) * Q(n2 + n5, 2) +
                  Q(n1, 1) * Q(n3 + n4, 2) * Q(n2 + n5, 2) -
                  2. * Q(n1 + n3 + n4, 3) * Q(n2 + n5, 2) +
                  2. * Q(n3, 1) * Q(n4, 1) * Q(n1 + n2 + n5, 3) -
                  2. * Q(n3 + n4, 2) * Q(n1 + n2 + n5, 3) -
                  Q(n1, 1) * Q(n2, 1) * Q(n4, 1) * Q(n3 + n5, 2) +
                  Q(n1 + n2, 2) * Q(n4, 1) * Q(n3 + n5, 2) +
                  Q(n2, 1) * Q(n1 + n4, 2) * Q(n3 + n5, 2) +
                  Q(n1, 1) * Q(n2 + n4, 2) * Q(n3 + n5, 2) -
                  2. * Q(n1 + n2 + n4, 3) * Q(n3 + n5, 2) +
                  2. * Q(n2, 1) * Q(n4, 1) * Q(n1 + n3 + n5, 3) -
                  2. * Q(n2 + n4, 2) * Q(n1 + n3 + n5, 3) +
                  2. * Q(n1, 1) * Q(n4, 1) * Q(n2 + n3 + n5, 3) -
                  2. * Q(n1 + n4, 2) * Q(n2 + n3 + n5, 3) -
                  6. * Q(n4, 1) * Q(n1 + n2 + n3 + n5, 4) -
                  Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4 + n5, 2) +
                  Q(n1 + n2, 2) * Q(n3, 1) * Q(n4 + n5, 2) +
                  Q(n2, 1) * Q(n1 + n3, 2) * Q(n4 + n5, 2) +
                  Q(n1, 1) * Q(n2 + n3, 2) * Q(n4 + n5, 2) -
                  2. * Q(n1 + n2 + n3, 3) * Q(n4 + n5, 2) +
                  2. * Q(n2, 1) * Q(n3, 1) * Q(n1 + n4 + n5, 3) -
                  2. * Q(n2 + n3, 2) * Q(n1 + n4 + n5, 3) +
                  2. * Q(n1, 1) * Q(n3, 1) * Q(n2 + n4 + n5, 3) -
                  2. * Q(n1 + n3, 2) * Q(n2 + n4 + n5, 3) -
                  6. * Q(n3, 1) * Q(n1 + n2 + n4 + n5, 4) +
                  2. * Q(n1, 1) * Q(n2, 1) * Q(n3 + n4 + n5, 3) -
                  2. * Q(n1 + n2, 2) * Q(n3 + n4 + n5, 3) -
                  6. * Q(n2, 1) * Q(n1 + n3 + n4 + n5, 4) -
                  6. * Q(n1, 1) * Q(n2 + n3 + n4 + n5, 4) +
                  24. * Q(n1 + n2 + n3 + n4 + n5, 5);
  return five;
}

TComplex AliAnalysisTaskAR::Six(Int_t n1, Int_t n2, Int_t n3, Int_t n4,
                                Int_t n5, Int_t n6) {
  // Generic six-particle correlation
  // <exp[i(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6)]>.
  TComplex six =
      Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) * Q(n6, 1) -
      Q(n1 + n2, 2) * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) * Q(n6, 1) -
      Q(n2, 1) * Q(n1 + n3, 2) * Q(n4, 1) * Q(n5, 1) * Q(n6, 1) -
      Q(n1, 1) * Q(n2 + n3, 2) * Q(n4, 1) * Q(n5, 1) * Q(n6, 1) +
      2. * Q(n1 + n2 + n3, 3) * Q(n4, 1) * Q(n5, 1) * Q(n6, 1) -
      Q(n2, 1) * Q(n3, 1) * Q(n1 + n4, 2) * Q(n5, 1) * Q(n6, 1) +
      Q(n2 + n3, 2) * Q(n1 + n4, 2) * Q(n5, 1) * Q(n6, 1) -
      Q(n1, 1) * Q(n3, 1) * Q(n2 + n4, 2) * Q(n5, 1) * Q(n6, 1) +
      Q(n1 + n3, 2) * Q(n2 + n4, 2) * Q(n5, 1) * Q(n6, 1) +
      2. * Q(n3, 1) * Q(n1 + n2 + n4, 3) * Q(n5, 1) * Q(n6, 1) -
      Q(n1, 1) * Q(n2, 1) * Q(n3 + n4, 2) * Q(n5, 1) * Q(n6, 1) +
      Q(n1 + n2, 2) * Q(n3 + n4, 2) * Q(n5, 1) * Q(n6, 1) +
      2. * Q(n2, 1) * Q(n1 + n3 + n4, 3) * Q(n5, 1) * Q(n6, 1) +
      2. * Q(n1, 1) * Q(n2 + n3 + n4, 3) * Q(n5, 1) * Q(n6, 1) -
      6. * Q(n1 + n2 + n3 + n4, 4) * Q(n5, 1) * Q(n6, 1) -
      Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n1 + n5, 2) * Q(n6, 1) +
      Q(n2 + n3, 2) * Q(n4, 1) * Q(n1 + n5, 2) * Q(n6, 1) +
      Q(n3, 1) * Q(n2 + n4, 2) * Q(n1 + n5, 2) * Q(n6, 1) +
      Q(n2, 1) * Q(n3 + n4, 2) * Q(n1 + n5, 2) * Q(n6, 1) -
      2. * Q(n2 + n3 + n4, 3) * Q(n1 + n5, 2) * Q(n6, 1) -
      Q(n1, 1) * Q(n3, 1) * Q(n4, 1) * Q(n2 + n5, 2) * Q(n6, 1) +
      Q(n1 + n3, 2) * Q(n4, 1) * Q(n2 + n5, 2) * Q(n6, 1) +
      Q(n3, 1) * Q(n1 + n4, 2) * Q(n2 + n5, 2) * Q(n6, 1) +
      Q(n1, 1) * Q(n3 + n4, 2) * Q(n2 + n5, 2) * Q(n6, 1) -
      2. * Q(n1 + n3 + n4, 3) * Q(n2 + n5, 2) * Q(n6, 1) +
      2. * Q(n3, 1) * Q(n4, 1) * Q(n1 + n2 + n5, 3) * Q(n6, 1) -
      2. * Q(n3 + n4, 2) * Q(n1 + n2 + n5, 3) * Q(n6, 1) -
      Q(n1, 1) * Q(n2, 1) * Q(n4, 1) * Q(n3 + n5, 2) * Q(n6, 1) +
      Q(n1 + n2, 2) * Q(n4, 1) * Q(n3 + n5, 2) * Q(n6, 1) +
      Q(n2, 1) * Q(n1 + n4, 2) * Q(n3 + n5, 2) * Q(n6, 1) +
      Q(n1, 1) * Q(n2 + n4, 2) * Q(n3 + n5, 2) * Q(n6, 1) -
      2. * Q(n1 + n2 + n4, 3) * Q(n3 + n5, 2) * Q(n6, 1) +
      2. * Q(n2, 1) * Q(n4, 1) * Q(n1 + n3 + n5, 3) * Q(n6, 1) -
      2. * Q(n2 + n4, 2) * Q(n1 + n3 + n5, 3) * Q(n6, 1) +
      2. * Q(n1, 1) * Q(n4, 1) * Q(n2 + n3 + n5, 3) * Q(n6, 1) -
      2. * Q(n1 + n4, 2) * Q(n2 + n3 + n5, 3) * Q(n6, 1) -
      6. * Q(n4, 1) * Q(n1 + n2 + n3 + n5, 4) * Q(n6, 1) -
      Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4 + n5, 2) * Q(n6, 1) +
      Q(n1 + n2, 2) * Q(n3, 1) * Q(n4 + n5, 2) * Q(n6, 1) +
      Q(n2, 1) * Q(n1 + n3, 2) * Q(n4 + n5, 2) * Q(n6, 1) +
      Q(n1, 1) * Q(n2 + n3, 2) * Q(n4 + n5, 2) * Q(n6, 1) -
      2. * Q(n1 + n2 + n3, 3) * Q(n4 + n5, 2) * Q(n6, 1) +
      2. * Q(n2, 1) * Q(n3, 1) * Q(n1 + n4 + n5, 3) * Q(n6, 1) -
      2. * Q(n2 + n3, 2) * Q(n1 + n4 + n5, 3) * Q(n6, 1) +
      2. * Q(n1, 1) * Q(n3, 1) * Q(n2 + n4 + n5, 3) * Q(n6, 1) -
      2. * Q(n1 + n3, 2) * Q(n2 + n4 + n5, 3) * Q(n6, 1) -
      6. * Q(n3, 1) * Q(n1 + n2 + n4 + n5, 4) * Q(n6, 1) +
      2. * Q(n1, 1) * Q(n2, 1) * Q(n3 + n4 + n5, 3) * Q(n6, 1) -
      2. * Q(n1 + n2, 2) * Q(n3 + n4 + n5, 3) * Q(n6, 1) -
      6. * Q(n2, 1) * Q(n1 + n3 + n4 + n5, 4) * Q(n6, 1) -
      6. * Q(n1, 1) * Q(n2 + n3 + n4 + n5, 4) * Q(n6, 1) +
      24. * Q(n1 + n2 + n3 + n4 + n5, 5) * Q(n6, 1) -
      Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) * Q(n1 + n6, 2) +
      Q(n2 + n3, 2) * Q(n4, 1) * Q(n5, 1) * Q(n1 + n6, 2) +
      Q(n3, 1) * Q(n2 + n4, 2) * Q(n5, 1) * Q(n1 + n6, 2) +
      Q(n2, 1) * Q(n3 + n4, 2) * Q(n5, 1) * Q(n1 + n6, 2) -
      2. * Q(n2 + n3 + n4, 3) * Q(n5, 1) * Q(n1 + n6, 2) +
      Q(n3, 1) * Q(n4, 1) * Q(n2 + n5, 2) * Q(n1 + n6, 2) -
      Q(n3 + n4, 2) * Q(n2 + n5, 2) * Q(n1 + n6, 2) +
      Q(n2, 1) * Q(n4, 1) * Q(n3 + n5, 2) * Q(n1 + n6, 2) -
      Q(n2 + n4, 2) * Q(n3 + n5, 2) * Q(n1 + n6, 2) -
      2. * Q(n4, 1) * Q(n2 + n3 + n5, 3) * Q(n1 + n6, 2) +
      Q(n2, 1) * Q(n3, 1) * Q(n4 + n5, 2) * Q(n1 + n6, 2) -
      Q(n2 + n3, 2) * Q(n4 + n5, 2) * Q(n1 + n6, 2) -
      2. * Q(n3, 1) * Q(n2 + n4 + n5, 3) * Q(n1 + n6, 2) -
      2. * Q(n2, 1) * Q(n3 + n4 + n5, 3) * Q(n1 + n6, 2) +
      6. * Q(n2 + n3 + n4 + n5, 4) * Q(n1 + n6, 2) -
      Q(n1, 1) * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) * Q(n2 + n6, 2) +
      Q(n1 + n3, 2) * Q(n4, 1) * Q(n5, 1) * Q(n2 + n6, 2) +
      Q(n3, 1) * Q(n1 + n4, 2) * Q(n5, 1) * Q(n2 + n6, 2) +
      Q(n1, 1) * Q(n3 + n4, 2) * Q(n5, 1) * Q(n2 + n6, 2) -
      2. * Q(n1 + n3 + n4, 3) * Q(n5, 1) * Q(n2 + n6, 2) +
      Q(n3, 1) * Q(n4, 1) * Q(n1 + n5, 2) * Q(n2 + n6, 2) -
      Q(n3 + n4, 2) * Q(n1 + n5, 2) * Q(n2 + n6, 2) +
      Q(n1, 1) * Q(n4, 1) * Q(n3 + n5, 2) * Q(n2 + n6, 2) -
      Q(n1 + n4, 2) * Q(n3 + n5, 2) * Q(n2 + n6, 2) -
      2. * Q(n4, 1) * Q(n1 + n3 + n5, 3) * Q(n2 + n6, 2) +
      Q(n1, 1) * Q(n3, 1) * Q(n4 + n5, 2) * Q(n2 + n6, 2) -
      Q(n1 + n3, 2) * Q(n4 + n5, 2) * Q(n2 + n6, 2) -
      2. * Q(n3, 1) * Q(n1 + n4 + n5, 3) * Q(n2 + n6, 2) -
      2. * Q(n1, 1) * Q(n3 + n4 + n5, 3) * Q(n2 + n6, 2) +
      6. * Q(n1 + n3 + n4 + n5, 4) * Q(n2 + n6, 2) +
      2. * Q(n3, 1) * Q(n4, 1) * Q(n5, 1) * Q(n1 + n2 + n6, 3) -
      2. * Q(n3 + n4, 2) * Q(n5, 1) * Q(n1 + n2 + n6, 3) -
      2. * Q(n4, 1) * Q(n3 + n5, 2) * Q(n1 + n2 + n6, 3) -
      2. * Q(n3, 1) * Q(n4 + n5, 2) * Q(n1 + n2 + n6, 3) +
      4. * Q(n3 + n4 + n5, 3) * Q(n1 + n2 + n6, 3) -
      Q(n1, 1) * Q(n2, 1) * Q(n4, 1) * Q(n5, 1) * Q(n3 + n6, 2) +
      Q(n1 + n2, 2) * Q(n4, 1) * Q(n5, 1) * Q(n3 + n6, 2) +
      Q(n2, 1) * Q(n1 + n4, 2) * Q(n5, 1) * Q(n3 + n6, 2) +
      Q(n1, 1) * Q(n2 + n4, 2) * Q(n5, 1) * Q(n3 + n6, 2) -
      2. * Q(n1 + n2 + n4, 3) * Q(n5, 1) * Q(n3 + n6, 2) +
      Q(n2, 1) * Q(n4, 1) * Q(n1 + n5, 2) * Q(n3 + n6, 2) -
      Q(n2 + n4, 2) * Q(n1 + n5, 2) * Q(n3 + n6, 2) +
      Q(n1, 1) * Q(n4, 1) * Q(n2 + n5, 2) * Q(n3 + n6, 2) -
      Q(n1 + n4, 2) * Q(n2 + n5, 2) * Q(n3 + n6, 2) -
      2. * Q(n4, 1) * Q(n1 + n2 + n5, 3) * Q(n3 + n6, 2) +
      Q(n1, 1) * Q(n2, 1) * Q(n4 + n5, 2) * Q(n3 + n6, 2) -
      Q(n1 + n2, 2) * Q(n4 + n5, 2) * Q(n3 + n6, 2) -
      2. * Q(n2, 1) * Q(n1 + n4 + n5, 3) * Q(n3 + n6, 2) -
      2. * Q(n1, 1) * Q(n2 + n4 + n5, 3) * Q(n3 + n6, 2) +
      6. * Q(n1 + n2 + n4 + n5, 4) * Q(n3 + n6, 2) +
      2. * Q(n2, 1) * Q(n4, 1) * Q(n5, 1) * Q(n1 + n3 + n6, 3) -
      2. * Q(n2 + n4, 2) * Q(n5, 1) * Q(n1 + n3 + n6, 3) -
      2. * Q(n4, 1) * Q(n2 + n5, 2) * Q(n1 + n3 + n6, 3) -
      2. * Q(n2, 1) * Q(n4 + n5, 2) * Q(n1 + n3 + n6, 3) +
      4. * Q(n2 + n4 + n5, 3) * Q(n1 + n3 + n6, 3) +
      2. * Q(n1, 1) * Q(n4, 1) * Q(n5, 1) * Q(n2 + n3 + n6, 3) -
      2. * Q(n1 + n4, 2) * Q(n5, 1) * Q(n2 + n3 + n6, 3) -
      2. * Q(n4, 1) * Q(n1 + n5, 2) * Q(n2 + n3 + n6, 3) -
      2. * Q(n1, 1) * Q(n4 + n5, 2) * Q(n2 + n3 + n6, 3) +
      4. * Q(n1 + n4 + n5, 3) * Q(n2 + n3 + n6, 3) -
      6. * Q(n4, 1) * Q(n5, 1) * Q(n1 + n2 + n3 + n6, 4) +
      6. * Q(n4 + n5, 2) * Q(n1 + n2 + n3 + n6, 4) -
      Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n5, 1) * Q(n4 + n6, 2) +
      Q(n1 + n2, 2) * Q(n3, 1) * Q(n5, 1) * Q(n4 + n6, 2) +
      Q(n2, 1) * Q(n1 + n3, 2) * Q(n5, 1) * Q(n4 + n6, 2) +
      Q(n1, 1) * Q(n2 + n3, 2) * Q(n5, 1) * Q(n4 + n6, 2) -
      2. * Q(n1 + n2 + n3, 3) * Q(n5, 1) * Q(n4 + n6, 2) +
      Q(n2, 1) * Q(n3, 1) * Q(n1 + n5, 2) * Q(n4 + n6, 2) -
      Q(n2 + n3, 2) * Q(n1 + n5, 2) * Q(n4 + n6, 2) +
      Q(n1, 1) * Q(n3, 1) * Q(n2 + n5, 2) * Q(n4 + n6, 2) -
      Q(n1 + n3, 2) * Q(n2 + n5, 2) * Q(n4 + n6, 2) -
      2. * Q(n3, 1) * Q(n1 + n2 + n5, 3) * Q(n4 + n6, 2) +
      Q(n1, 1) * Q(n2, 1) * Q(n3 + n5, 2) * Q(n4 + n6, 2) -
      Q(n1 + n2, 2) * Q(n3 + n5, 2) * Q(n4 + n6, 2) -
      2. * Q(n2, 1) * Q(n1 + n3 + n5, 3) * Q(n4 + n6, 2) -
      2. * Q(n1, 1) * Q(n2 + n3 + n5, 3) * Q(n4 + n6, 2) +
      6. * Q(n1 + n2 + n3 + n5, 4) * Q(n4 + n6, 2) +
      2. * Q(n2, 1) * Q(n3, 1) * Q(n5, 1) * Q(n1 + n4 + n6, 3) -
      2. * Q(n2 + n3, 2) * Q(n5, 1) * Q(n1 + n4 + n6, 3) -
      2. * Q(n3, 1) * Q(n2 + n5, 2) * Q(n1 + n4 + n6, 3) -
      2. * Q(n2, 1) * Q(n3 + n5, 2) * Q(n1 + n4 + n6, 3) +
      4. * Q(n2 + n3 + n5, 3) * Q(n1 + n4 + n6, 3) +
      2. * Q(n1, 1) * Q(n3, 1) * Q(n5, 1) * Q(n2 + n4 + n6, 3) -
      2. * Q(n1 + n3, 2) * Q(n5, 1) * Q(n2 + n4 + n6, 3) -
      2. * Q(n3, 1) * Q(n1 + n5, 2) * Q(n2 + n4 + n6, 3) -
      2. * Q(n1, 1) * Q(n3 + n5, 2) * Q(n2 + n4 + n6, 3) +
      4. * Q(n1 + n3 + n5, 3) * Q(n2 + n4 + n6, 3) -
      6. * Q(n3, 1) * Q(n5, 1) * Q(n1 + n2 + n4 + n6, 4) +
      6. * Q(n3 + n5, 2) * Q(n1 + n2 + n4 + n6, 4) +
      2. * Q(n1, 1) * Q(n2, 1) * Q(n5, 1) * Q(n3 + n4 + n6, 3) -
      2. * Q(n1 + n2, 2) * Q(n5, 1) * Q(n3 + n4 + n6, 3) -
      2. * Q(n2, 1) * Q(n1 + n5, 2) * Q(n3 + n4 + n6, 3) -
      2. * Q(n1, 1) * Q(n2 + n5, 2) * Q(n3 + n4 + n6, 3) +
      4. * Q(n1 + n2 + n5, 3) * Q(n3 + n4 + n6, 3) -
      6. * Q(n2, 1) * Q(n5, 1) * Q(n1 + n3 + n4 + n6, 4) +
      6. * Q(n2 + n5, 2) * Q(n1 + n3 + n4 + n6, 4) -
      6. * Q(n1, 1) * Q(n5, 1) * Q(n2 + n3 + n4 + n6, 4) +
      6. * Q(n1 + n5, 2) * Q(n2 + n3 + n4 + n6, 4) +
      24. * Q(n5, 1) * Q(n1 + n2 + n3 + n4 + n6, 5) -
      Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n5 + n6, 2) +
      Q(n1 + n2, 2) * Q(n3, 1) * Q(n4, 1) * Q(n5 + n6, 2) +
      Q(n2, 1) * Q(n1 + n3, 2) * Q(n4, 1) * Q(n5 + n6, 2) +
      Q(n1, 1) * Q(n2 + n3, 2) * Q(n4, 1) * Q(n5 + n6, 2) -
      2. * Q(n1 + n2 + n3, 3) * Q(n4, 1) * Q(n5 + n6, 2) +
      Q(n2, 1) * Q(n3, 1) * Q(n1 + n4, 2) * Q(n5 + n6, 2) -
      Q(n2 + n3, 2) * Q(n1 + n4, 2) * Q(n5 + n6, 2) +
      Q(n1, 1) * Q(n3, 1) * Q(n2 + n4, 2) * Q(n5 + n6, 2) -
      Q(n1 + n3, 2) * Q(n2 + n4, 2) * Q(n5 + n6, 2) -
      2. * Q(n3, 1) * Q(n1 + n2 + n4, 3) * Q(n5 + n6, 2) +
      Q(n1, 1) * Q(n2, 1) * Q(n3 + n4, 2) * Q(n5 + n6, 2) -
      Q(n1 + n2, 2) * Q(n3 + n4, 2) * Q(n5 + n6, 2) -
      2. * Q(n2, 1) * Q(n1 + n3 + n4, 3) * Q(n5 + n6, 2) -
      2. * Q(n1, 1) * Q(n2 + n3 + n4, 3) * Q(n5 + n6, 2) +
      6. * Q(n1 + n2 + n3 + n4, 4) * Q(n5 + n6, 2) +
      2. * Q(n2, 1) * Q(n3, 1) * Q(n4, 1) * Q(n1 + n5 + n6, 3) -
      2. * Q(n2 + n3, 2) * Q(n4, 1) * Q(n1 + n5 + n6, 3) -
      2. * Q(n3, 1) * Q(n2 + n4, 2) * Q(n1 + n5 + n6, 3) -
      2. * Q(n2, 1) * Q(n3 + n4, 2) * Q(n1 + n5 + n6, 3) +
      4. * Q(n2 + n3 + n4, 3) * Q(n1 + n5 + n6, 3) +
      2. * Q(n1, 1) * Q(n3, 1) * Q(n4, 1) * Q(n2 + n5 + n6, 3) -
      2. * Q(n1 + n3, 2) * Q(n4, 1) * Q(n2 + n5 + n6, 3) -
      2. * Q(n3, 1) * Q(n1 + n4, 2) * Q(n2 + n5 + n6, 3) -
      2. * Q(n1, 1) * Q(n3 + n4, 2) * Q(n2 + n5 + n6, 3) +
      4. * Q(n1 + n3 + n4, 3) * Q(n2 + n5 + n6, 3) -
      6. * Q(n3, 1) * Q(n4, 1) * Q(n1 + n2 + n5 + n6, 4) +
      6. * Q(n3 + n4, 2) * Q(n1 + n2 + n5 + n6, 4) +
      2. * Q(n1, 1) * Q(n2, 1) * Q(n4, 1) * Q(n3 + n5 + n6, 3) -
      2. * Q(n1 + n2, 2) * Q(n4, 1) * Q(n3 + n5 + n6, 3) -
      2. * Q(n2, 1) * Q(n1 + n4, 2) * Q(n3 + n5 + n6, 3) -
      2. * Q(n1, 1) * Q(n2 + n4, 2) * Q(n3 + n5 + n6, 3) +
      4. * Q(n1 + n2 + n4, 3) * Q(n3 + n5 + n6, 3) -
      6. * Q(n2, 1) * Q(n4, 1) * Q(n1 + n3 + n5 + n6, 4) +
      6. * Q(n2 + n4, 2) * Q(n1 + n3 + n5 + n6, 4) -
      6. * Q(n1, 1) * Q(n4, 1) * Q(n2 + n3 + n5 + n6, 4) +
      6. * Q(n1 + n4, 2) * Q(n2 + n3 + n5 + n6, 4) +
      24. * Q(n4, 1) * Q(n1 + n2 + n3 + n5 + n6, 5) +
      2. * Q(n1, 1) * Q(n2, 1) * Q(n3, 1) * Q(n4 + n5 + n6, 3) -
      2. * Q(n1 + n2, 2) * Q(n3, 1) * Q(n4 + n5 + n6, 3) -
      2. * Q(n2, 1) * Q(n1 + n3, 2) * Q(n4 + n5 + n6, 3) -
      2. * Q(n1, 1) * Q(n2 + n3, 2) * Q(n4 + n5 + n6, 3) +
      4. * Q(n1 + n2 + n3, 3) * Q(n4 + n5 + n6, 3) -
      6. * Q(n2, 1) * Q(n3, 1) * Q(n1 + n4 + n5 + n6, 4) +
      6. * Q(n2 + n3, 2) * Q(n1 + n4 + n5 + n6, 4) -
      6. * Q(n1, 1) * Q(n3, 1) * Q(n2 + n4 + n5 + n6, 4) +
      6. * Q(n1 + n3, 2) * Q(n2 + n4 + n5 + n6, 4) +
      24. * Q(n3, 1) * Q(n1 + n2 + n4 + n5 + n6, 5) -
      6. * Q(n1, 1) * Q(n2, 1) * Q(n3 + n4 + n5 + n6, 4) +
      6. * Q(n1 + n2, 2) * Q(n3 + n4 + n5 + n6, 4) +
      24. * Q(n2, 1) * Q(n1 + n3 + n4 + n5 + n6, 5) +
      24. * Q(n1, 1) * Q(n2 + n3 + n4 + n5 + n6, 5) -
      120. * Q(n1 + n2 + n3 + n4 + n5 + n6, 6);
  return six;
}

TComplex AliAnalysisTaskAR::Recursion(Int_t n, Int_t *harmonic,
                                      Int_t mult /* = 1*/, Int_t skip /*= 0*/) {
  // Calculate multi-particle correlators by using recursion (an improved
  // faster version) originally developed by Kristjan Gulbrandsen
  // (gulbrand@nbi.dk).

  Int_t nm1 = n - 1;
  TComplex c(Q(harmonic[nm1], mult));
  if (nm1 == 0)
    return c;
  c *= Recursion(nm1, harmonic);
  if (nm1 == skip)
    return c;

  Int_t multp1 = mult + 1;
  Int_t nm2 = n - 2;
  Int_t counter1 = 0;
  Int_t hhold = harmonic[counter1];
  harmonic[counter1] = harmonic[nm2];
  harmonic[nm2] = hhold + harmonic[nm1];
  TComplex c2(Recursion(nm1, harmonic, multp1, nm2));
  Int_t counter2 = n - 3;
  while (counter2 >= skip) {
    harmonic[nm2] = harmonic[counter1];
    harmonic[counter1] = hhold;
    ++counter1;
    hhold = harmonic[counter1];
    harmonic[counter1] = harmonic[nm2];
    harmonic[nm2] = hhold + harmonic[nm1];
    c2 += Recursion(nm1, harmonic, multp1, counter2);
    --counter2;
  }
  harmonic[nm2] = harmonic[counter1];
  harmonic[counter1] = hhold;

  if (mult == 1)
    return c - c2;
  return c - Double_t(mult) * c2;
}

TComplex AliAnalysisTaskAR::TwoNestedLoops(Int_t n1, Int_t n2,
                                           std::vector<Double_t> angles,
                                           std::vector<Double_t> weights) {
  // Calculation of <cos(n1*phi1+n2*phi2)> and <sin(n1*phi1+n2*phi2)>
  // with two nested loops

  TComplex Two(0., 0.);

  Double_t phi1 = 0., phi2 = 0.; // particle angle
  Double_t w1 = 1., w2 = 1.;     // particle weight
  for (std::size_t i1 = 0; i1 < angles.size(); i1++) {
    phi1 = angles.at(i1);
    w1 = weights.at(i1);
    for (std::size_t i2 = 0; i2 < angles.size(); i2++) {
      if (i2 == i1) {
        continue;
      } // Get rid of autocorrelations
      phi2 = angles.at(i2);
      w2 = weights.at(i2);
      Two += TComplex(w1 * w2 * TMath::Cos(n1 * phi1 + n2 * phi2),
                      w1 * w2 * TMath::Sin(n1 * phi1 + n2 * phi2));
    }
  }
  return Two;
}

TComplex AliAnalysisTaskAR::ThreeNestedLoops(Int_t n1, Int_t n2, Int_t n3,
                                             std::vector<Double_t> angles,
                                             std::vector<Double_t> weights) {
  // Calculation of <cos(n1*phi1+n2*phi2+n3*phi3)> and
  // <sin(n1*phi1+n2*phi2+n3*phi3)> with three nested loops.

  TComplex Q(0., 0.);
  Double_t phi1 = 0., phi2 = 0., phi3 = 0.; // particle angle
  Double_t w1 = 1., w2 = 1., w3 = 1.;       // particle weight
  for (std::size_t i1 = 0; i1 < angles.size(); i1++) {
    phi1 = angles.at(i1);
    w1 = weights.at(i1);
    for (std::size_t i2 = 0; i2 < angles.size(); i2++) {
      if (i2 == i1) {
        continue;
      } // Get rid of autocorrelations
      phi2 = angles.at(i2);
      w2 = weights.at(i2);
      for (std::size_t i3 = 0; i3 < angles.size(); i3++) {
        if (i3 == i1 || i3 == i2) {
          continue;
        } // Get rid of autocorrelations
        phi3 = angles.at(i3);
        w3 = weights.at(i3);
        Q += TComplex(
            w1 * w2 * w3 * TMath::Cos(n1 * phi1 + n2 * phi2 + n3 * phi3),
            w1 * w2 * w3 * TMath::Sin(n1 * phi1 + n2 * phi2 + n3 * phi3));
      }
    }
  }
  return Q;
}

TComplex AliAnalysisTaskAR::FourNestedLoops(Int_t n1, Int_t n2, Int_t n3,
                                            Int_t n4,
                                            std::vector<Double_t> angles,
                                            std::vector<Double_t> weights) {
  // Calculation of <cos(n1*phi1+n2*phi2+n3*phi3+n4*phi4)> and
  // <sin(n1*phi1+n2*phi2+n3*phi3+n4*phi4)> with four nested loops.

  TComplex Q(0., 0.);
  Double_t phi1 = 0., phi2 = 0., phi3 = 0., phi4 = 0.; // particle angle
  Double_t w1 = 1., w2 = 1., w3 = 1., w4 = 1.;         // particle weight
  for (std::size_t i1 = 0; i1 < angles.size(); i1++) {
    phi1 = angles.at(i1);
    w1 = weights.at(i1);
    for (std::size_t i2 = 0; i2 < angles.size(); i2++) {
      if (i2 == i1) {
        continue;
      } // Get rid of autocorrelations
      phi2 = angles.at(i2);
      w2 = weights.at(i2);
      for (std::size_t i3 = 0; i3 < angles.size(); i3++) {
        if (i3 == i1 || i3 == i2) {
          continue;
        } // Get rid of autocorrelations
        phi3 = angles.at(i3);
        w3 = weights.at(i3);
        for (std::size_t i4 = 0; i4 < angles.size(); i4++) {
          if (i4 == i1 || i4 == i2 || i4 == i3) {
            continue;
          } // Get rid of autocorrelations
          phi4 = angles.at(i4);
          w4 = weights.at(i4);
          Q += TComplex(
              w1 * w2 * w3 * w4 *
                  TMath::Cos(n1 * phi1 + n2 * phi2 + n3 * phi3 + n4 * phi4),
              w1 * w2 * w3 * w4 *
                  TMath::Sin(n1 * phi1 + n2 * phi2 + n3 * phi3 + n4 * phi4));
        }
      }
    }
  }
  return Q;
}

TComplex AliAnalysisTaskAR::FiveNestedLoops(Int_t n1, Int_t n2, Int_t n3,
                                            Int_t n4, Int_t n5,
                                            std::vector<Double_t> angles,
                                            std::vector<Double_t> weights) {
  // Calculation of <cos(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5)> and
  // <sin(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5)> with four nested loops.

  TComplex Q(0., 0.);
  Double_t phi1 = 0., phi2 = 0., phi3 = 0., phi4 = 0.,
           phi5 = 0.;                                   // particle angle
  Double_t w1 = 1., w2 = 1., w3 = 1., w4 = 1., w5 = 1.; // particle weight
  for (std::size_t i1 = 0; i1 < angles.size(); i1++) {
    phi1 = angles.at(i1);
    w1 = weights.at(i1);
    for (std::size_t i2 = 0; i2 < angles.size(); i2++) {
      if (i2 == i1) {
        continue;
      } // Get rid of autocorrelations
      phi2 = angles.at(i2);
      w2 = weights.at(i2);
      for (std::size_t i3 = 0; i3 < angles.size(); i3++) {
        if (i3 == i1 || i3 == i2) {
          continue;
        } // Get rid of autocorrelations
        phi3 = angles.at(i3);
        w3 = weights.at(i3);
        for (std::size_t i4 = 0; i4 < angles.size(); i4++) {
          if (i4 == i1 || i4 == i2 || i4 == i3) {
            continue;
          } // Get rid of autocorrelations
          phi4 = angles.at(i4);
          w4 = weights.at(i4);
          for (std::size_t i5 = 0; i5 < angles.size(); i5++) {
            if (i5 == i1 || i5 == i2 || i5 == i3 || i5 == i4) {
              continue;
            } // Get rid of autocorrelations
            phi5 = angles.at(i5);
            w5 = weights.at(i5);
            Q += TComplex(w1 * w2 * w3 * w4 * w5 *
                              TMath::Cos(n1 * phi1 + n2 * phi2 + n3 * phi3 +
                                         n4 * phi4 + n5 * phi5),
                          w1 * w2 * w3 * w4 * w5 *
                              TMath::Sin(n1 * phi1 + n2 * phi2 + n3 * phi3 +
                                         n4 * phi4 + n5 * phi5));
          }
        }
      }
    }
  }
  return Q;
}

TComplex AliAnalysisTaskAR::SixNestedLoops(Int_t n1, Int_t n2, Int_t n3,
                                           Int_t n4, Int_t n5, Int_t n6,
                                           std::vector<Double_t> angles,
                                           std::vector<Double_t> weights) {
  // Calculation of <cos(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5*n6*phi6)>
  // and <sin(n1*phi1+n2*phi2+n3*phi3+n4*phi4+n5*phi5+n6*phi6)> with four
  // nested loops.

  TComplex Q(0., 0.);
  Double_t phi1 = 0., phi2 = 0., phi3 = 0., phi4 = 0., phi5 = 0.,
           phi6 = 0.; // particle angle
  Double_t w1 = 1., w2 = 1., w3 = 1., w4 = 1., w5 = 1.,
           w6 = 1.; // particle weight
  for (std::size_t i1 = 0; i1 < angles.size(); i1++) {
    phi1 = angles.at(i1);
    w1 = weights.at(i1);
    for (std::size_t i2 = 0; i2 < angles.size(); i2++) {
      if (i2 == i1) {
        continue;
      } // Get rid of autocorrelations
      phi2 = angles.at(i2);
      w2 = weights.at(i2);
      for (std::size_t i3 = 0; i3 < angles.size(); i3++) {
        if (i3 == i1 || i3 == i2) {
          continue;
        } // Get rid of autocorrelations
        phi3 = angles.at(i3);
        w3 = weights.at(i3);
        for (std::size_t i4 = 0; i4 < angles.size(); i4++) {
          if (i4 == i1 || i4 == i2 || i4 == i3) {
            continue;
          } // Get rid of autocorrelations
          phi4 = angles.at(i4);
          w4 = weights.at(i4);
          for (std::size_t i5 = 0; i5 < angles.size(); i5++) {
            if (i5 == i1 || i5 == i2 || i5 == i3 || i5 == i4) {
              continue;
            } // Get rid of autocorrelations
            phi5 = angles.at(i5);
            w5 = weights.at(i5);
            for (std::size_t i6 = 0; i6 < angles.size(); i6++) {
              if (i6 == i1 || i6 == i2 || i6 == i3 || i6 == i4 || i6 == i5) {
                continue;
              } // Get rid of autocorrelations
              phi6 = angles.at(i6);
              w6 = weights.at(i6);
              Q += TComplex(w1 * w2 * w3 * w4 * w5 * w6 *
                                TMath::Cos(n1 * phi1 + n2 * phi2 + n3 * phi3 +
                                           n4 * phi4 + n5 * phi5 + n6 * phi6),
                            w1 * w2 * w3 * w4 * w5 * w6 *
                                TMath::Sin(n1 * phi1 + n2 * phi2 + n3 * phi3 +
                                           n4 * phi4 + n5 * phi5 + n6 * phi6));
            }
          }
        }
      }
    }
  }
  return Q;
}

void AliAnalysisTaskAR::SetCenCorQAHistogramBinning(
    Int_t cen1, Int_t xnbins, Double_t xlowerEdge, Double_t xupperEdge,
    Int_t cen2, Int_t ynbins, Double_t ylowerEdge, Double_t yupperEdge) {
  if (cen1 >= LAST_ECENESTIMATORS || cen2 >= LAST_ECENESTIMATORS) {
    std::cout << __LINE__ << ": running out of bounds" << std::endl;
    Fatal("SetCenCorQAHistogramBinning",
          "Running out of bounds in SetCenCorQAHistogramBinning");
  }
  if (xupperEdge < xlowerEdge && yupperEdge < ylowerEdge) {
    std::cout << __LINE__ << ": upper edge has to be larger than the lower edge"
              << std::endl;
    Fatal("SetCenCorQAHistogramBinning",
          ": upper edge has to be larger than the lower edge");
  }
  this->fCenCorQAHistogramBins[IndexCorHistograms(
      cen1, cen2, LAST_ECENESTIMATORS)][kBIN] = xnbins;
  this->fCenCorQAHistogramBins[IndexCorHistograms(
      cen1, cen2, LAST_ECENESTIMATORS)][kLEDGE] = xlowerEdge;
  this->fCenCorQAHistogramBins[IndexCorHistograms(
      cen1, cen2, LAST_ECENESTIMATORS)][kUEDGE] = xupperEdge;
  this->fCenCorQAHistogramBins[IndexCorHistograms(
      cen1, cen2, LAST_ECENESTIMATORS)][kBIN + LAST_EBINS] = ynbins;
  this->fCenCorQAHistogramBins[IndexCorHistograms(
      cen1, cen2, LAST_ECENESTIMATORS)][kLEDGE + LAST_EBINS] = ylowerEdge;
  this->fCenCorQAHistogramBins[IndexCorHistograms(
      cen1, cen2, LAST_ECENESTIMATORS)][kUEDGE + LAST_EBINS] = yupperEdge;
}

void AliAnalysisTaskAR::SetMulCorQAHistogramBinning(
    Int_t mul1, Int_t xnbins, Double_t xlowerEdge, Double_t xupperEdge,
    Int_t mul2, Int_t ynbins, Double_t ylowerEdge, Double_t yupperEdge) {
  if (mul1 >= kMulEstimators || mul2 >= kMulEstimators) {
    std::cout << __LINE__ << ": running out of bounds" << std::endl;
    Fatal("SetMulCorQAHistogramBinning",
          "Running out of bounds in SetMulCorQAHistogramBinning");
  }
  if (xupperEdge < xlowerEdge && yupperEdge < ylowerEdge) {
    std::cout << __LINE__ << ": upper edge has to be larger than the lower edge"
              << std::endl;
    Fatal("SetMulCorQAHistogramBinning",
          ": upper edge has to be larger than the lower edge");
  }
  this->fMulCorQAHistogramBins[IndexCorHistograms(mul1, mul2, kMulEstimators)]
                              [kBIN] = xnbins;
  this->fMulCorQAHistogramBins[IndexCorHistograms(mul1, mul2, kMulEstimators)]
                              [kLEDGE] = xlowerEdge;
  this->fMulCorQAHistogramBins[IndexCorHistograms(mul1, mul2, kMulEstimators)]
                              [kUEDGE] = xupperEdge;
  this->fMulCorQAHistogramBins[IndexCorHistograms(mul1, mul2, kMulEstimators)]
                              [kBIN + LAST_EBINS] = ynbins;
  this->fMulCorQAHistogramBins[IndexCorHistograms(mul1, mul2, kMulEstimators)]
                              [kLEDGE + LAST_EBINS] = ylowerEdge;
  this->fMulCorQAHistogramBins[IndexCorHistograms(mul1, mul2, kMulEstimators)]
                              [kUEDGE + LAST_EBINS] = yupperEdge;
}

void AliAnalysisTaskAR::SetCenMulCorQAHistogramBinning(
    Int_t xnbins, Double_t xlowerEdge, Double_t xupperEdge, Int_t mul,
    Int_t ynbins, Double_t ylowerEdge, Double_t yupperEdge) {
  if (mul >= kMulEstimators) {
    std::cout << __LINE__ << ": running out of bounds" << std::endl;
    Fatal("SetCenMulCorQAHistogramBinning",
          "Running out of bounds in SetMulCorQAHistogramBinning");
  }
  if (xupperEdge < xlowerEdge && yupperEdge < ylowerEdge) {
    std::cout << __LINE__ << ": upper edge has to be larger than the lower edge"
              << std::endl;
    Fatal("SetCenMulCorQAHistogramBinning",
          ": upper edge has to be larger than the lower edge");
  }
  this->fCenMulCorQAHistogramBins[mul][kBIN] = xnbins;
  this->fCenMulCorQAHistogramBins[mul][kLEDGE] = xlowerEdge;
  this->fCenMulCorQAHistogramBins[mul][kUEDGE] = xupperEdge;
  this->fCenMulCorQAHistogramBins[mul][kBIN + LAST_EBINS] = ynbins;
  this->fCenMulCorQAHistogramBins[mul][kLEDGE + LAST_EBINS] = ylowerEdge;
  this->fCenMulCorQAHistogramBins[mul][kUEDGE + LAST_EBINS] = yupperEdge;
}

void AliAnalysisTaskAR::SetAcceptanceHistogram(kTrack kinematic,
                                               const char *Filename,
                                               const char *Histname) {
  // set a acceptance histograms

  // check if index is out of range
  if (kinematic > kKinematic) {
    std::cout << __LINE__ << ": Out of range" << std::endl;
    Fatal("SetAccpetanceHistogram", "Out of range");
  }
  // check if file exists
  if (gSystem->AccessPathName(Filename, kFileExists)) {
    std::cout << __LINE__ << ": File does not exist" << std::endl;
    Fatal("SetAcceptanceHistogram", "Invalid file name");
  }
  TFile *file = new TFile(Filename, "READ");
  if (!file) {
    std::cout << __LINE__ << ": Cannot open file" << std::endl;
    Fatal("SetAcceptanceHistogram", "ROOT file cannot be read");
  }
  this->fAcceptanceHistogram[kinematic] =
      dynamic_cast<TH1D *>(file->Get(Histname));
  if (!fAcceptanceHistogram[kinematic]) {
    std::cout << __LINE__ << ": No acceptance histogram" << std::endl;
    Fatal("SetAcceptanceHistogram", "Cannot get acceptance histogram");
  }
  // keeps the histogram in memory after we close the file
  this->fAcceptanceHistogram[kinematic]->SetDirectory(0);
  file->Close();
}

void AliAnalysisTaskAR::SetWeightHistogram(kTrack kinematic,
                                           const char *Filename,
                                           const char *Histname) {
  // set weight histogram

  // check if index is out of range
  if (kinematic > kKinematic) {
    std::cout << __LINE__ << ": Out of range" << std::endl;
    Fatal("SetAccpetanceHistogram", "Out of range");
  }
  // check if file exists
  if (gSystem->AccessPathName(Filename, kFileExists)) {
    std::cout << __LINE__ << ": File does not exist" << std::endl;
    Fatal("SetWeightHistogram", "Invalid file name");
  }
  TFile *file = new TFile(Filename, "READ");
  if (!file) {
    std::cout << __LINE__ << ": Cannot open file" << std::endl;
    Fatal("SetWeightHistogram", "ROOT file cannot be read");
  }
  this->fWeightHistogram[kinematic] = dynamic_cast<TH1D *>(file->Get(Histname));
  if (!fWeightHistogram[kinematic]) {
    std::cout << __LINE__ << ": No acceptance histogram" << std::endl;
    Fatal("SetWeightHistogram", "Cannot get weight histogram");
  }
  // keeps the histogram in memory after we close the file
  this->fWeightHistogram[kinematic]->SetDirectory(0);
  file->Close();
  this->fUseWeights = kTRUE;
}

void AliAnalysisTaskAR::SetCenFlattenHist(const char *Filename,
                                          const char *Histname) {
  // get histogram for centrality flattening
  // check if file exists
  if (gSystem->AccessPathName(Filename, kFileExists)) {
    std::cout << __LINE__ << ": File does not exist" << std::endl;
    Fatal("SetCenFlattenHist", "Invalid file name");
  }
  TFile *file = new TFile(Filename, "READ");
  if (!file) {
    std::cout << __LINE__ << ": Cannot open file" << std::endl;
    Fatal("SetCenFlattenHist", "ROOT file cannot be read");
  }
  this->fCenFlattenHist = dynamic_cast<TH1D *>(file->Get(Histname));
  if (!fCenFlattenHist) {
    std::cout << __LINE__ << ": No histogram" << std::endl;
    Fatal("SetCenFlattenHist", "Cannot get histogram");
  }
  // keeps the histogram in memory after we close the file
  this->fCenFlattenHist->SetDirectory(0);
  file->Close();
  this->fUseCenFlatten = kTRUE;
}

void AliAnalysisTaskAR::GetPointers(TList *histList) {
  // Initialize pointer for base list fHistList so we can initialize all of
  // objects and call terminate off-line

  fHistList = histList;
  if (!fHistList) {
    std::cout << __LINE__ << ": Did not get " << fHistListName << std::endl;
    Fatal("GetPointers", "Invalid Pointer");
  }

  // initialize all other objects
  this->GetPointersForControlHistograms();
  this->GetPointersForQAHistograms();
  this->GetPointersForFinalResults();
}

void AliAnalysisTaskAR::GetPointersForQAHistograms() {
  // get pointers for QA Histograms

  // get pointer for fControlHistograms
  fQAHistogramsList =
      dynamic_cast<TList *>(fHistList->FindObject(fQAHistogramsListName));
  // if the pointer is null, then there was no QA
  if (!fQAHistogramsList) {
    return;
  }

  // get pointer for fCenCorQAHistogramsList, if it is there
  fCenCorQAHistogramsList = dynamic_cast<TList *>(
      fQAHistogramsList->FindObject(fCenCorQAHistogramsListName));

  if (fCenCorQAHistogramsList) {
    // get pointers for centrality correlation histograms
    for (int cen = 0; cen < LAST_ECENESTIMATORS * (LAST_ECENESTIMATORS - 1) / 2;
         ++cen) {
      for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
        fCenCorQAHistograms[cen][ba] =
            dynamic_cast<TH2D *>(fCenCorQAHistogramsList->FindObject(
                fCenCorQAHistogramNames[cen][ba][kNAME]));
        if (!fCenCorQAHistograms[cen][ba]) {
          std::cout << __LINE__ << ": Did not get "
                    << fCenCorQAHistogramNames[cen][ba][kNAME] << std::endl;
          Fatal("GetPointersForQAHistograms", "Invalid Pointer");
        }
      }
    }
  }

  // get pointer for fMulCorQAHistogramsList, if it is there
  fMulCorQAHistogramsList = dynamic_cast<TList *>(
      fQAHistogramsList->FindObject(fMulCorQAHistogramsListName));

  if (fMulCorQAHistogramsList) {
    // get pointers for multiplicity correlation histograms
    for (int mul = 0; mul < kMulEstimators; ++mul) {
      for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
        fMulCorQAHistograms[mul][ba] =
            dynamic_cast<TH2D *>(fMulCorQAHistogramsList->FindObject(
                fMulCorQAHistogramNames[mul][ba][kNAME]));
        if (!fMulCorQAHistograms[mul][ba]) {
          std::cout << __LINE__ << ": Did not get "
                    << fMulCorQAHistogramNames[mul][ba][kNAME] << std::endl;
          Fatal("GetPointersForQAHistograms", "Invalid Pointer");
        }
      }
    }
  }

  // get pointer for fCenMulCorQAHistogramsList, if it is there
  fCenMulCorQAHistogramsList = dynamic_cast<TList *>(
      fQAHistogramsList->FindObject(fCenMulCorQAHistogramsListName));

  if (fCenMulCorQAHistogramsList) {
    // get pointers for multiplicity correlation histograms
    for (int mul = 0; mul < kMulEstimators; ++mul) {
      for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
        fCenMulCorQAHistograms[mul][ba] =
            dynamic_cast<TH2D *>(fCenMulCorQAHistogramsList->FindObject(
                fCenMulCorQAHistogramNames[mul][ba][kNAME]));
        if (!fCenMulCorQAHistograms[mul][ba]) {
          std::cout << __LINE__ << ": Did not get "
                    << fCenMulCorQAHistogramNames[mul][ba][kNAME] << std::endl;
          Fatal("GetPointersForQAHistograms", "Invalid Pointer");
        }
      }
    }
  }

  // get pointer for fFBScanQAHistogramsList, if it is there
  fFBScanQAHistogramsList = dynamic_cast<TList *>(
      fQAHistogramsList->FindObject(fFBScanQAHistogramsListName));

  if (fFBScanQAHistogramsList) {
    // get pointer for filter bit scan histogram
    fFBScanQAHistogram = dynamic_cast<TH1D *>(
        fFBScanQAHistogramsList->FindObject(fFBScanQAHistogramName[kNAME]));

    // get pointer track scan filterbit QA histograms
    for (int track = 0; track < LAST_ETRACK; ++track) {
      for (int fb = 0; fb < kNumberofTestFilterBit; ++fb) {
        fFBTrackScanQAHistograms[track][fb] =
            dynamic_cast<TH1D *>(fFBScanQAHistogramsList->FindObject(
                fFBTrackScanQAHistogramNames[track][fb][kNAME]));
        if (!fFBTrackScanQAHistograms[track][fb]) {
          std::cout << __LINE__ << ": Did not get "
                    << fFBTrackScanQAHistogramNames[track][fb][kNAME]
                    << std::endl;
          Fatal("GetPointersForQAHistograms", "Invalid Pointer");
        }
      }
    }
  }

  // get pointer for fSelfCorQAHistogramsList, if it is there
  fSelfCorQAHistogramsList = dynamic_cast<TList *>(
      fQAHistogramsList->FindObject(fSelfCorQAHistogramsListName));

  if (fSelfCorQAHistogramsList) {
    // get pointers for self correlation QA histograms
    for (int var = 0; var < kKinematic; ++var) {
      for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
        fSelfCorQAHistograms[var][ba] =
            dynamic_cast<TH1D *>(fSelfCorQAHistogramsList->FindObject(
                fSelfCorQAHistogramNames[var][ba][kNAME]));
        if (!fSelfCorQAHistograms[var][ba]) {
          std::cout << __LINE__ << ": Did not get "
                    << fSelfCorQAHistogramNames[var][ba][kNAME] << std::endl;
          Fatal("GetPointersForQAHistograms", "Invalid Pointer");
        }
      }
    }
  }
}

void AliAnalysisTaskAR::GetPointersForControlHistograms() {
  // get pointers for Control Histograms

  // get pointer for fControlHistograms
  fControlHistogramsList =
      dynamic_cast<TList *>(fHistList->FindObject(fControlHistogramsListName));
  if (!fControlHistogramsList) {
    std::cout << __LINE__ << ": Did not get " << fControlHistogramsListName
              << std::endl;
    Fatal("GetPointersForControlHistograms", "Invalid Pointer");
  }

  // get pointer for fTrackControlHistogramsList
  fTrackControlHistogramsList = dynamic_cast<TList *>(
      fControlHistogramsList->FindObject(fTrackControlHistogramsListName));
  if (!fTrackControlHistogramsList) {
    std::cout << __LINE__ << ": Did not get " << fTrackControlHistogramsListName
              << std::endl;
    Fatal("GetPointersForControlHistograms", "Invalid Pointer");
  }

  // get pointers for track cut counter histograms
  for (int mode = 0; mode < LAST_EMODE; ++mode) {
    fTrackCutsCounter[mode] = dynamic_cast<TH1D *>(
        fTrackControlHistogramsList->FindObject(fTrackCutsCounterNames[mode]));
  }

  fTrackCutsValues = dynamic_cast<TProfile *>(
      fTrackControlHistogramsList->FindObject(fTrackCutsValuesName));

  // fTrackCutsCounterCumulative = dynamic_cast<THnSparseD *>(
  //     fTrackControlHistogramsList->FindObject(fTrackCutsCounterCumulativeName));

  // get all pointers for track control histograms
  for (int mode = 0; mode < LAST_EMODE; ++mode) {
    for (int var = 0; var < LAST_ETRACK; ++var) {
      for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
        fTrackControlHistograms[mode][var][ba] =
            dynamic_cast<TH1D *>(fTrackControlHistogramsList->FindObject(
                fTrackControlHistogramNames[mode][var][ba][kNAME]));
        if (!fTrackControlHistograms[mode][var][ba]) {
          std::cout << __LINE__ << ": Did not get "
                    << fTrackControlHistogramNames[mode][var][ba][kNAME]
                    << std::endl;
          Fatal("GetPointersForControlHistograms", "Invalid Pointer");
        }
      }
    }
  }

  // get pointer for fEventControlHistogramsList
  fEventControlHistogramsList = dynamic_cast<TList *>(
      fControlHistogramsList->FindObject(fEventControlHistogramsListName));
  if (!fEventControlHistogramsList) {
    std::cout << __LINE__ << ": Did not get " << fEventControlHistogramsListName
              << std::endl;
    Fatal("GetPointersForControlHistograms", "Invalid Pointer");
  }

  // get pointers for event cut counter histograms
  for (int mode = 0; mode < LAST_EMODE; ++mode) {
    fEventCutsCounter[mode] = dynamic_cast<TH1D *>(
        fEventControlHistogramsList->FindObject(fEventCutsCounterNames[mode]));
  }

  fEventCutsValues = dynamic_cast<TProfile *>(
      fEventControlHistogramsList->FindObject(fEventCutsValuesName));

  // fEventCutsCounterCumulative = dynamic_cast<THnSparseD *>(
  //     fEventControlHistogramsList->FindObject(fEventCutsCounterCumulativeName));

  // get all pointers for event control histograms
  for (int mode = 0; mode < LAST_EMODE; ++mode) {
    for (int var = 0; var < LAST_EEVENT; ++var) {
      for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
        fEventControlHistograms[mode][var][ba] =
            dynamic_cast<TH1D *>(fEventControlHistogramsList->FindObject(
                fEventControlHistogramNames[mode][var][ba][kNAME]));
        if (!fEventControlHistograms[mode][var][ba]) {
          std::cout << __LINE__ << ": Did not get "
                    << fEventControlHistogramNames[mode][var][ba][kNAME]
                    << std::endl;
          Fatal("GetPointersForControlHistograms", "Invalid Pointer");
        }
      }
    }
  }
}

void AliAnalysisTaskAR::GetPointersForFinalResults() {
  // Get pointers for all final result

  // Get pointer for fFinalResultsList
  fFinalResultsList =
      dynamic_cast<TList *>(fHistList->FindObject(fFinalResultsListName));
  if (!fFinalResultsList) {
    std::cout << __LINE__ << ": Did not get " << fFinalResultsListName
              << std::endl;
    Fatal("GetPointersForOutputHistograms", "Invalid Pointer");
  }

  // get pointers for fFinalResultHistogramsList-
  fFinalResultHistogramsList = dynamic_cast<TList *>(
      fFinalResultsList->FindObject(fFinalResultHistogramsListName));
  if (!fFinalResultHistogramsList) {
    std::cout << __LINE__ << ": Did not get " << fFinalResultHistogramsListName
              << std::endl;
    Fatal("GetPointersForFinalResults", "Invalid Pointer");
  }
  // get pointers for all final result histograms
  for (int var = 0; var < LAST_EFINALHIST; ++var) {
    fFinalResultHistograms[var] =
        dynamic_cast<TH1D *>(fFinalResultHistogramsList->FindObject(
            fFinalResultHistogramNames[var][kNAME]));
    if (!fFinalResultHistograms[var]) {
      std::cout << __LINE__ << ": Did not get "
                << fFinalResultHistogramNames[var][kNAME] << std::endl;
      Fatal("GetPointersForFinalResults", "Invalid Pointer");
    }
  }

  // get pointers for fFinalResultCorrelatorsList
  fFinalResultCorrelatorsList = dynamic_cast<TList *>(
      fFinalResultsList->FindObject(fFinalResultCorrelatorsListName));
  if (!fFinalResultCorrelatorsList) {
    std::cout << __LINE__ << ": Did not get " << fFinalResultCorrelatorsListName
              << std::endl;
    Fatal("GetPointersForFinalResults", "Invalid Pointer");
  }

  // get pointers for fFinalResultCorrelatorsList
  fFinalResultSymmetricCumulantsList = dynamic_cast<TList *>(
      fFinalResultsList->FindObject(fFinalResultSymmetricCumulantsListName));
  if (!fFinalResultSymmetricCumulantsList) {
    std::cout << __LINE__ << ": Did not get "
              << fFinalResultSymmetricCumulantsListName << std::endl;
    Fatal("GetPointersForFinalResults", "Invalid Pointer");
  }

  // initialize vectors
  TString name;
  std::vector<Int_t> sc;
  fSymmetricCumulants.clear();
  fCorrelators.clear();
  fMapSCtoCor.clear();

  // get pointers for fFinalResultSymmetricCumulantsList
  for (auto list : *fFinalResultSymmetricCumulantsList) {
    sc.clear();
    name = TString(list->GetName());
    name.ReplaceAll("SC(", "");
    name.ReplaceAll(")", "");
    name.ReplaceAll(",", "");
    for (Int_t i = 0; i < name.Length(); i++) {
      sc.push_back(TString(name[i]).Atoi());
    }
    fSymmetricCumulants.push_back(sc);
  }

  // get pointers for fFinalResultNormalizedSymmetricCumulantsList
  fFinalResultNormalizedSymmetricCumulantsList =
      dynamic_cast<TList *>(fFinalResultsList->FindObject(
          fFinalResultNormalizedSymmetricCumulantsListName));
  if (!fFinalResultNormalizedSymmetricCumulantsList) {
    std::cout << __LINE__ << ": Did not get "
              << fFinalResultNormalizedSymmetricCumulantsListName << std::endl;
    Fatal("GetPointersForFinalResults", "Invalid Pointer");
  }

  Int_t Index = 0;
  std::vector<std::vector<Int_t>> correlators;
  for (std::size_t j = 0; j < fSymmetricCumulants.size(); j++) {
    correlators = MapSCToCor(fSymmetricCumulants.at(j));
    fMapSCtoCor.insert({fSymmetricCumulants.at(j), correlators});
    for (auto cor : correlators) {
      if (std::find(fCorrelators.begin(), fCorrelators.end(), cor) !=
          fCorrelators.end()) {
        continue;
      } else {
        fCorrelators.push_back(cor);
        fMapCorToIndex.insert({cor, Index});
        Index++;
      }
    }
  }
}
