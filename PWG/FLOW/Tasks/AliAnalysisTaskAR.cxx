/**
 * File              : AliAnalysisTaskAR.cxx
 * Author            : Anton Riedel <anton.riedel@tum.de>
 * Date              : 07.05.2021
 * Last Modified Date: 29.07.2021
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

#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskAR.h"
#include "AliLog.h"
#include "AliMultSelection.h"
#include <TColor.h>
#include <TFile.h>
#include <TH1.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <cstdlib>
#include <iostream>

ClassImp(AliAnalysisTaskAR)

    AliAnalysisTaskAR::AliAnalysisTaskAR(const char *name,
                                         Bool_t useParticleWeights)
    : AliAnalysisTaskSE(name),
      // Base list for all output objects
      fHistList(nullptr), fHistListName("outputStudentAnalysis"),
      // list holding all QA histograms
      fQAHistogramsList(nullptr), fQAHistogramsListName("QAHistograms"),
      fFillQAHistograms(kFALSE),
      // list holding centrality estimator correlation histograms
      fCenCorQAHistogramsList(nullptr),
      fCenCorQAHistogramsListName("CenCorQAHistograms"),
      // list holding filterbit scan histograms
      fFBScanQAHistogramsList(nullptr),
      fFBScanQAHistogramsListName("FBScanQAHistograms"),
      // list holding self correlation QA histograms
      fSelfCorQAHistogramsList(nullptr),
      fSelfCorQAHistogramsListName("SelfCorQAHistograms"),
      // list holding all control histograms
      fControlHistogramsList(nullptr),
      fControlHistogramsListName("ControlHistograms"),
      // sublists for track control histograms
      fTrackControlHistogramsList(nullptr),
      fTrackControlHistogramsListName("TrackControlHistograms"),
      // sublists for event control histograms
      fEventControlHistogramsList(nullptr),
      fEventControlHistogramsListName("EventControlHistograms"),
      // cuts
      fCentralitySelCriterion("V0M"), fTrackCutsCounter(nullptr),
      fEventCutsCounter(nullptr), fFilterbit(128), fPrimaryOnly(kFALSE),
      // Final results
      fFinalResultsList(nullptr), fFinalResultsListName("FinalResults"),
      // flags for MC analysis
      fMCAnalysisList(nullptr), fMCAnalysisListName("MCAnalysis"),
      fMCAnalaysis(kFALSE), fMCClosure(kFALSE), fSeed(0),
      fUseCustomSeed(kFALSE), fMCPdf(nullptr), fMCPdfName("pdf"),
      fMCFlowHarmonics(nullptr),
      fMCNumberOfParticlesPerEventFluctuations(kFALSE),
      fMCNumberOfParticlesPerEvent(500),
      // qvectors
      fQvectorList(nullptr), fPhi({}), fWeights({}),
      fAcceptanceHistogram(nullptr), fWeightHistogram(nullptr),
      fUseWeights(kFALSE), fResetWeights(kFALSE), fCorrelators({}) {
  // Constructor

  AliDebug(2, "AliAnalysisTaskAR::AliAnalysisTaskAR(const "
              "char *name, Bool_t useParticleWeights)");

  // Base list
  fHistList = new TList();
  fHistList->SetName(fHistListName);
  fHistList->SetOwner(kTRUE);

  // Initialize all arrays
  this->InitializeArrays();

  // Define input and output slots here
  // Input slot #0 works with an AliFlowEventSimple
  // DefineInput(0, AliFlowEventSimple::Class());
  // Input slot #1 is needed for the weights input file:
  // if(useParticleWeights)
  //{
  // DefineInput(1, TList::Class());
  //}
  // Output slot #0 is reserved
  // Output slot #1 writes into a TList container

  DefineOutput(1, TList::Class());

  // if (useParticleWeights) {
  // TBI
  // }
}

AliAnalysisTaskAR::AliAnalysisTaskAR()
    : AliAnalysisTaskSE(),
      // Dummy constructor
      // Base list for all output objects
      fHistList(nullptr), fHistListName("outputStudentAnalysis"),
      // list holding all QA histograms
      fQAHistogramsList(nullptr), fQAHistogramsListName("QAHistograms"),
      fFillQAHistograms(kFALSE),
      // list holding centrality estimator correlation histograms
      fCenCorQAHistogramsList(nullptr),
      fCenCorQAHistogramsListName("CenCorQAHistograms"),
      // list holding filterbit scan histograms
      fFBScanQAHistogramsList(nullptr),
      fFBScanQAHistogramsListName("FBScanQAHistograms"),
      // list holding self correlation QA histograms
      fSelfCorQAHistogramsList(nullptr),
      fSelfCorQAHistogramsListName("SelfCorQAHistograms"),
      // list holding all control histograms
      fControlHistogramsList(nullptr),
      fControlHistogramsListName("ControlHistograms"),
      // sublists for track control histograms
      fTrackControlHistogramsList(nullptr),
      fTrackControlHistogramsListName("TrackControlHistograms"),
      // sublists for event control histograms
      fEventControlHistogramsList(nullptr),
      fEventControlHistogramsListName("EventControlHistograms"),
      // cuts
      fCentralitySelCriterion("V0M"), fTrackCutsCounter(nullptr),
      fEventCutsCounter(nullptr), fFilterbit(128), fPrimaryOnly(kFALSE),
      // Final results
      fFinalResultsList(nullptr), fFinalResultsListName("FinalResults"),
      // flags for MC analysis
      fMCAnalysisList(nullptr), fMCAnalysisListName("MCAnalysis"),
      fMCAnalaysis(kFALSE), fMCClosure(kFALSE), fSeed(0),
      fUseCustomSeed(kFALSE), fMCPdf(nullptr), fMCPdfName("pdf"),
      fMCFlowHarmonics(nullptr),
      fMCNumberOfParticlesPerEventFluctuations(kFALSE),
      fMCNumberOfParticlesPerEvent(500),
      // qvectors
      fQvectorList(nullptr), fPhi({}), fWeights({}),
      fAcceptanceHistogram(nullptr), fWeightHistogram(nullptr),
      fUseWeights(kFALSE), fResetWeights(kFALSE), fCorrelators({}) {
  // initialize arrays in dummy constructor !!!!
  this->InitializeArrays();

  AliDebug(2, "AliAnalysisTaskAR::AliAnalysisTaskAR()");
}

AliAnalysisTaskAR::~AliAnalysisTaskAR() {
  // Destructor

  // fHlist owns all other data members, if we delete it, we will recursively
  // delete all other objects associative with this object
  if (fHistList) {
    delete fHistList;
  }

  if (fMCAnalaysis || fMCClosure) {
    delete gRandom;
  }
  if (fMCAnalaysis) {
    delete fMCPdf;
  }
};

void AliAnalysisTaskAR::UserCreateOutputObjects() {
  // Called at every worker node to initialize.

  // 1) Trick to avoid name clashes, part 1;
  // 2) Book and nest all lists;
  // 3) Book all objects;
  // *) Trick to avoid name clashes, part 2.

  // 1) Trick to avoid name clashes, part 1
  Bool_t oldHistAddStatus = TH1::AddDirectoryStatus();
  TH1::AddDirectory(kFALSE);

  // 2) Book and nest all lists
  this->BookAndNestAllLists();

  // 3) Book all objects
  if (fFillQAHistograms) {
    this->BookQAHistograms();
  }
  this->BookControlHistograms();
  this->BookFinalResultHistograms();
  this->BookFinalResultProfiles();
  if (fMCAnalaysis) {
    this->BookMCObjects();
  }
  if (fMCAnalaysis || fMCClosure) {
    delete gRandom;
    fUseCustomSeed ? gRandom = new TRandom3(fSeed) : gRandom = new TRandom3(0);
  }

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

  // Do some calculation in offline mode here:

  // ... your code for offline calculations ...

  /* get average value of phi and write it into its own histogram */
  fFinalResultHistograms[kPHIAVG]->SetBinContent(
      1, fTrackControlHistograms[kPHI][kAFTER]->GetMean());

  // for local MC Analysis
  // compute analytical values for correlators and fill them into profile
  if (fMCAnalaysis) {
    Double_t theoryValue = 1.;
    for (auto V : fCorrelators) {
      theoryValue = 1.;
      for (auto i : V) {
        theoryValue *= fMCFlowHarmonics->GetAt(abs(i) - 1);
      }
      fFinalResultProfiles[kHARTHEO]->Fill(V.size() - 1.5, theoryValue);
    }
  }
}

void AliAnalysisTaskAR::InitializeArrays() {
  // Initialize all data members which are arrays in this method
  InitializeArraysForQAHistograms();
  InitializeArraysForTrackControlHistograms();
  InitializeArraysForEventControlHistograms();
  InitializeArraysForCuts();
  InitializeArraysForQvectors();
  InitializeArraysForFinalResultHistograms();
  InitializeArraysForFinalResultProfiles();
  InitializeArraysForMCAnalysis();
}

void AliAnalysisTaskAR::InitializeArraysForQAHistograms() {
  // initialize array of QA histograms for the correlation between centrality
  // estimators
  // there are N(N-1)/2 such correlators, i.e. the number of elemets above the
  // diagonal of a square matrix
  for (int cen = 0; cen < LAST_ECENESTIMATORS * (LAST_ECENESTIMATORS - 1) / 2;
       ++cen) {
    for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
      fCenCorQAHistograms[cen][ba] = nullptr;
    }
  }

  // names for QA histograms for the correlation between centrality estimators
  TString CenCorQAHistogramNames[LAST_ECENESTIMATORS *
                                 (LAST_ECENESTIMATORS - 1) / 2][LAST_ENAME] = {
      // NAME, TITLE, XAXIS, YAXIS
      {"fCorCenEstimatorQAHistograms[kV0M+kCL0]", "V0M vs CL0", "V0M", "CL0"},
      {"fCorCenEstimatorQAHistograms[kV0M+kCL1]", "V0M vs CL1", "V0M", "CL1"},
      {"fCorCenEstimatorQAHistograms[kV0M+kSPDTRACKLETS]",
       "V0M vs SPDTracklets", "V0M", "SPDTracklets"},
      {"fCorCenEstimatorQAHistograms[kCL0+kCL1]", "CL0 vs CL1", "CL0", "CL1"},
      {"fCorCenEstimatorQAHistograms[kCL0+kSPDTRACKLETS]",
       "CL0 vs kSPDTRACKLETS", "CL0", "SPDTracklets"},
      {"fCorCenEstimatorQAHistograms[kCL1+kSPDTRACKLETS]",
       "CL1 vs kSPDTRACKLETS", "CL1", "SPDTracklets"},
  };

  // initialize names of QA histograms for the correlation between centrality
  // estimators
  for (int cen = 0; cen < LAST_ECENESTIMATORS * (LAST_ECENESTIMATORS - 1) / 2;
       ++cen) {
    for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
      for (int name = 0; name < LAST_ENAME; ++name) {
        if (name == kNAME || name == kTITLE) {
          fCenCorQAHistogramNames[cen][ba][name] =
              CenCorQAHistogramNames[cen][name] + kBAName[ba];
        } else
          fCenCorQAHistogramNames[cen][ba][name] =
              CenCorQAHistogramNames[cen][name];
      }
    }
  }

  // default bins for QA histograms for the correlation between centrality
  // estimators
  Double_t CorCenEstimatorQAHistogramBins[LAST_ECENESTIMATORS *
                                          (LAST_ECENESTIMATORS - 1) /
                                          2][2 * LAST_EBINS] = {
      // kBIN kLEDGE kUEDGE kBIN kLEDGE kUEDGE
      {50., 0., 100., 50., 0., 100.}, // kV0M + kCL0
      {50., 0., 100., 50., 0., 100.}, // kV0M + kCL1
      {50., 0., 100., 50., 0., 100.}, // kV0M + kSPDTRACKLETS
      {50., 0., 100., 50., 0., 100.}, // CL0 + kCL1
      {50., 0., 100., 50., 0., 100.}, // CL0 + kSPDTRACKLETS
      {50., 0., 100., 50., 0., 100.}, // CL1 + kSPDTRACKLETS
  };
  // initialize array of bins and edges for QA histograms for the correlation
  // between centrality
  for (int cen = 0; cen < LAST_ECENESTIMATORS * (LAST_ECENESTIMATORS - 1) / 2;
       ++cen) {
    for (int bin = 0; bin < 2 * LAST_EBINS; ++bin) {
      fCenCorQAHistogramBins[cen][bin] =
          CorCenEstimatorQAHistogramBins[cen][bin];
    }
  }

  // initialize array for filterbit scan QA histograms
  fFBScanQAHistogram = nullptr;
  // names for filterbit scan QA histograms
  TString FBScanQAHistogramName[LAST_ENAME] = // NAME, TITLE, XAXIS, YAXIS
      {"fFBScanQAHistograms", "Filterbit Scan", "Filterbit", ""};
  // initialize names for track control histograms
  for (int name = 0; name < LAST_ENAME; ++name) {
    fFBScanQAHistogramName[name] = FBScanQAHistogramName[name];
  }

  // default bins for filterbit scan histograms
  Double_t FBScanQAHistogramBin[LAST_EBINS] = {kMaxFilterbit, 0.,
                                               kMaxFilterbit};
  // initialize array of bins and edges for filterbit scan histogram
  for (int bin = 0; bin < LAST_EBINS; ++bin) {
    fFBScanQAHistogramBin[bin] = FBScanQAHistogramBin[bin];
  }

  // initialize array of track filterbit scan QA histograms
  for (int track = 0; track < LAST_ETRACK; ++track) {
    for (int fb = 0; fb < kNumberofTestFilterBit; ++fb) {
      fFBTrackScanQAHistograms[track][fb] = nullptr;
    }
  }

  // names for track filterbit scan QA histograms
  TString FBTrackScanQAHistogramNames[LAST_ETRACK][LAST_ENAME] = {
      // NAME, TITLE, XAXIS, YAXIS
      {"fFBTrackScanQAHistogram[kPT]", "Filterbitscan p_{T}", "p_{t}", ""},
      {"fFBTrackScanQAHistogram[kPHI]", "Filterbitscan #varphi", "#varphi", ""},
      {"fFBTrackScanQAHistogram[kETA]", "Filterbitscan #eta", "#eta", ""},
      {"fFBTrackScanQAHistogram[kCHARGE]", "Filterbitscan Charge", "Q", ""},
      {"fFBTrackScanQAHistogram[kTPCNCLS]",
       "Filterbitscan number of TPC clusters", "", ""},
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
  // default bins
  Double_t FBTrackScanHistogramBins[LAST_ETRACK][LAST_EBINS] = {
      // kBIN kLEDGE kUEDGE
      {100., 0., 10.},            // kPT
      {360., 0., TMath::TwoPi()}, // kPHI
      {200., -2., 2.},            // kETA
      {5., -2.5, 2.5},            // kCHARGE
      {160., 0., 160.},           // kTPCNCLS
      {10., 0., 10.},             // kITSNCLS
      {100., 0., 10.},            // kCHI2PERNDF
      {100., -10., 10.},          // kDCAZ
      {100, -10., 10.},           // kDCAXY
  };
  // initialize array of bins and edges
  for (int var = 0; var < LAST_ETRACK; ++var) {
    for (int bin = 0; bin < LAST_EBINS; ++bin) {
      fFBTrackScanQAHistogramBins[var][bin] =
          FBTrackScanHistogramBins[var][bin];
    }
  }

  // initialize arrays for self correlation QA histograms
  for (int var = 0; var < kKinematic; ++var) {
    for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
      fSelfCorQAHistograms[var][ba] = nullptr;
    }
  }

  // names for self correlation QA histograms
  TString SelfCorQAHistogramNames[kKinematic][LAST_ENAME] = {
      // NAME, TITLE, XAXIS, YAXIS
      {"fSelfCorQAHistograms[kPT]", "p_{T}^{1}-p_{T}^{2}", "#Delta p_{T}", ""},
      {"fSelfCorQAHistograms[kPHI]", "#varphi_{1}-#varphi_{2}",
       "#Delta #varphi", ""},
      {"fSelfCorQAHistograms[kETA]", "#eta_{1}-#eta_{2}", "#Delta #eta", ""},
  };

  // initialize names for self correlation QA histograms
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

  Double_t SelfCorQAHistogramBins[kKinematic][LAST_EBINS] = {
      // kBIN kLEDGE kUEDGE
      {100., -0.1, 0.1}, // kPT
      {100., -0.1, 0.1}, // kPHI
      {100., -0.1, 0.1}, // kETA
  };
  // initialize array of bins and edges for track control histograms
  for (int var = 0; var < kKinematic; ++var) {
    for (int bin = 0; bin < LAST_EBINS; ++bin) {
      fSelfCorQAHistogramBins[var][bin] = SelfCorQAHistogramBins[var][bin];
    }
  }
}

void AliAnalysisTaskAR::InitializeArraysForTrackControlHistograms() {
  // initialize array of track control histograms
  for (int var = 0; var < LAST_ETRACK; ++var) {
    for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
      fTrackControlHistograms[var][ba] = nullptr;
    }
  }

  // names for track control histograms
  TString TrackControlHistogramNames[LAST_ETRACK][LAST_ENAME] = {
      // NAME, TITLE, XAXIS, YAXIS
      {"fTrackControlHistograms[kPT]", "p_{T}", "p_{T}", ""},
      {"fTrackControlHistograms[kPHI]", "#varphi", "#varphi", ""},
      {"fTrackControlHistograms[kETA]", "#eta", "#eta", ""},
      {"fTrackControlHistograms[kCHARGE]", "Charge", "Q", ""},
      {"fTrackControlHistograms[kTPCNCLS]", "Number of clusters in TPC", ""},
      {"fTrackControlHistograms[kITSNCLS]", "Number of clusters in ITS", ""},
      {"fTrackControlHistograms[kCHI2PERNDF]", "CHI2PERNDF of track", "", ""},
      {"fTrackControlHistograms[kDCAZ]", "DCA in Z", ""},
      {"fTrackControlHistograms[kDCAXY]", "DCA in XY", ""}, // kBEFORE
  };
  // initialize names for track control histograms
  for (int var = 0; var < LAST_ETRACK; ++var) {
    for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
      for (int name = 0; name < LAST_ENAME; ++name) {
        if (name == kNAME || name == kTITLE) {
          fTrackControlHistogramNames[var][ba][name] =
              TrackControlHistogramNames[var][name] + kBAName[ba];
        } else {
          fTrackControlHistogramNames[var][ba][name] =
              TrackControlHistogramNames[var][name];
        }
      }
    }
  }

  // default bins for track control histograms
  Double_t BinsTrackControlHistogramDefaults[LAST_ETRACK][LAST_EBINS] = {
      // kBIN kLEDGE kUEDGE
      {100., 0., 10.},            // kPT
      {360., 0., TMath::TwoPi()}, // kPHI
      {200., -2., 2.},            // kETA
      {7., -3.5, 3.5},            // kCHARGE
      {160., 0., 160.},           // kTPCNCLS
      {1000., 0., 1000.},         // kITSNCLS
      {100., 0., 10.},            // kCHI2PERNDF
      {100., -10., 10.},          // kDCAZ
      {100, -10., 10.},           // kDCAXY
  };
  // initialize array of bins and edges for track control histograms
  for (int var = 0; var < LAST_ETRACK; ++var) {
    for (int bin = 0; bin < LAST_EBINS; ++bin) {
      fTrackControlHistogramBins[var][bin] =
          BinsTrackControlHistogramDefaults[var][bin];
    }
  }
  // initialize bin names of track cuts counter histogram
  TString TrackCutsCounterBinNames[LAST_ETRACK] = {
      "kPT",      "kPHI",        "kETA",  "kCHARGE", "kTPCNCLS",
      "kITSNCLS", "kCHI2PERNDF", "kDCAZ", "kDCAXY",
  };
  for (int name = 0; name < LAST_ETRACK; name++) {
    fTrackCutsCounterBinNames[name] = TrackCutsCounterBinNames[name];
  }
}

void AliAnalysisTaskAR::InitializeArraysForEventControlHistograms() {
  // initialize array of event control histograms
  for (int var = 0; var < LAST_EEVENT; ++var) {
    for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
      fEventControlHistograms[var][ba] = nullptr;
    }
  }

  // name of event control histograms
  TString EventControlHistogramNames[LAST_EEVENT][LAST_ENAME] = {
      // NAME, TITLE, XAXIS, YAXIS
      {"fEventControlHistograms[kX]", "Primary vertex X", "X", ""},
      {"fEventControlHistograms[kY]", "Primary vertex Y", "Y", ""},
      {"fEventControlHistograms[kZ]", "Primary vertex Z", "Z", ""},
      {"fEventControlHistograms[kCEN]", "centrality", "Centrality Percentile",
       ""},
      {"fEventControlHistograms[kMUL]", "multiplicity", "M", ""},
      {"fEventControlHistograms[kNCONTRIB]", "Number of Contributers",
       "#Contributors", ""},
  };
  // initialize names for event control histograms
  for (int var = 0; var < LAST_EEVENT; ++var) {
    for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
      for (int name = 0; name < LAST_ENAME; ++name) {
        if (name == kNAME || name == kTITLE) {
          fEventControlHistogramNames[var][ba][name] =
              EventControlHistogramNames[var][name] + kBAName[ba];
        } else
          fEventControlHistogramNames[var][ba][name] =
              EventControlHistogramNames[var][name];
      }
    }
  }

  // default bins for event control histograms
  Double_t BinsEventControlHistogramDefaults[LAST_EEVENT][LAST_EBINS] = {
      // kBIN kLEDGE kUEDGE
      {40., -20., 20.},   // kX
      {40., -20., 20.},   // kY
      {40., -20., 20.},   // kZ
      {10., 0., 100},     // kCEN
      {200., 0., 20000.}, // kMUL
      {100., 0., 5000.},  // kNCONTRIB
  };
  // initialize array of bins and edges for track control histograms
  for (int var = 0; var < LAST_EEVENT; ++var) {
    for (int bin = 0; bin < LAST_EBINS; ++bin) {
      fEventControlHistogramBins[var][bin] =
          BinsEventControlHistogramDefaults[var][bin];
    }
  }
  // initialize bin names of event cuts counter histogram
  TString EventCutsCounterBinNames[LAST_EEVENT] = {
      "kX", "kY", "kZ", "kCEN", "kMUL", "kNCONTRIB",
  };
  for (int name = 0; name < LAST_EEVENT; name++) {
    fEventCutsCounterBinNames[name] = EventCutsCounterBinNames[name];
  }
}

void AliAnalysisTaskAR::InitializeArraysForCuts() {
  // initialize all arrays for cuts

  // default track cuts
  Double_t TrackCutDefaults[LAST_ETRACK][LAST_EMINMAX] = {
      // MIN MAX
      {0., 5.},             // kPT
      {0., TMath::TwoPi()}, // kPHI
      {-3., 3.},            // kETA
      {0., 160.},           // kTPCNCLS
      {0., 10.},            // kITSNCLS
      {0., 10.},            // kCHI2PERNDF
      {-10., 10},           // kDCAZ
      {-10., 10},           // kDCAXY
  };
  // initialize array for track cuts
  for (int var = 0; var < LAST_ETRACK; ++var) {
    for (int mm = 0; mm < LAST_EMINMAX; ++mm) {
      fTrackCuts[var][mm] = TrackCutDefaults[var][mm];
    }
  }

  // default event cuts
  Double_t EventCutDefaults[LAST_EEVENT][LAST_EMINMAX]{
      // MIN MAX
      {-20., 20.}, // kX
      {-20., 20.}, // kY
      {-20., 20.}, // kZ
      {0., 100.},  // kCEN
      {0., 20000}, // kMUL
      {0., 1e6},   // kNCONTRIB
  };
  // initialize array for event cuts
  for (int var = 0; var < LAST_EEVENT; ++var) {
    for (int mm = 0; mm < LAST_EMINMAX; ++mm) {
      fEventCuts[var][mm] = EventCutDefaults[var][mm];
    }
  }
}

void AliAnalysisTaskAR::InitializeArraysForFinalResultHistograms() {
  // initialize array for final result histograms
  for (int var = 0; var < LAST_EFINALHIST; ++var) {
    fFinalResultHistograms[var] = nullptr;
  }

  TString FinalResultHistogramNames[LAST_EFINALHIST][LAST_ENAME] = {
      // NAME, TITLE, XAXIS
      {"fFinalResultHistograms[PHIAVG]", "Average #varphi", "#varphi",
       ""}, // PHIAVG
  };

  // initialize names for final result histograms
  for (int var = 0; var < LAST_EFINALHIST; ++var) {
    for (int name = 0; name < LAST_ENAME; ++name) {
      fFinalResultHistogramNames[var][name] =
          FinalResultHistogramNames[var][name];
    }
  }

  // default bins for final result histograms
  Double_t BinsFinalResultHistogramDefaults[LAST_EFINALHIST][LAST_EBINS] = {
      // BIN LEDGE UEDGE
      {1., 0., 1.}, // AVGPHI
  };
  // initialize array of bins and edges for track control histograms
  for (int var = 0; var < LAST_EFINALHIST; ++var) {
    for (int bin = 0; bin < LAST_EBINS; ++bin) {
      fFinalResultHistogramBins[var][bin] =
          BinsFinalResultHistogramDefaults[var][bin];
    }
  }
}

void AliAnalysisTaskAR::InitializeArraysForQvectors() {
  // Make sure all Q-vectors are initially zero
  for (Int_t h = 0; h < kMaxHarmonic; h++) {
    for (Int_t p = 0; p < kMaxPower; p++) {
      fQvector[h][p] = TComplex(0., 0.);
    }
  }
}

void AliAnalysisTaskAR::InitializeArraysForFinalResultProfiles() {
  // initialize array for final result profiles
  for (int var = 0; var < LAST_EFINALPROFILE; ++var) {
    fFinalResultProfiles[var] = nullptr;
  }

  TString FinalResultProfileNames[LAST_EFINALPROFILE][LAST_ENAME] = {
      // kNAME, kTITLE, kXAXIS
      {"fFinalResultProfiles[kHARDATA]", "Flow Harmonics (Data)", "",
       ""}, // kHARDATA
      {"fFinalResultProfiles[kHARDATARESET]",
       "Flow Harmonics (Data, weights reset)", ""}, // kHARDATARESET
      {"fFinalResultProfiles[kHARTHEO]", "Flow Harmonics (Theory)", "",
       ""}, // kHARTHEO
  };

  // initialize names for final result profiles
  for (int var = 0; var < LAST_EFINALPROFILE; ++var) {
    for (int name = 0; name < LAST_ENAME; ++name) {
      fFinalResultProfileNames[var][name] = FinalResultProfileNames[var][name];
    }
  }

  // default bins for final result histograms
  Double_t BinsFinalResultProfileDefaults[LAST_EFINALPROFILE][LAST_EBINS] = {
      // kBIN kLEDGE kUEDGE
      {1., 0., 1.}, // kHARDATA
      {1., 0., 1.}, // kHARDATARESET
      {1., 0., 1.}, // kHARTHEO
  };
  // initialize array of bins and edges for final result profiles
  for (int var = 0; var < LAST_EFINALPROFILE; ++var) {
    for (int bin = 0; bin < LAST_EBINS; ++bin) {
      fFinalResultProfileBins[var][bin] =
          BinsFinalResultProfileDefaults[var][bin];
    }
  }
}

void AliAnalysisTaskAR::InitializeArraysForMCAnalysis() {
  // initialize arrays for MC analysis

  // range of pdf
  Double_t MCPdfRangeDefaults[LAST_EMINMAX] = {0.0, TMath::TwoPi()};
  for (int i = 0; i < LAST_EMINMAX; ++i) {
    fMCPdfRange[i] = MCPdfRangeDefaults[i];
  }

  // range of fluctuations of number of particles produced per event
  Int_t MCNumberOfParticlesPerEventRangeDefaults[LAST_EMINMAX] = {500, 1000};
  for (int i = 0; i < LAST_EMINMAX; ++i) {
    fMCNumberOfParticlesPerEventRange[i] =
        MCNumberOfParticlesPerEventRangeDefaults[i];
  }
}

void AliAnalysisTaskAR::BookAndNestAllLists() {
  // Book and nest all lists nested in the base list fHistList

  // 1. Book and nest lists for QA histograms
  // 2. Book and nest lists for control histograms
  // 3. Book and nest lists for final results

  if (!fHistList) {
    std::cout << __LINE__ << ": Did not get " << fHistListName << std::endl;
    Fatal("BookAndNestAllLists", "Invalid Pointer");
  }
  // 1. Book and nest lists for QA histograms
  if (fFillQAHistograms) {
    fQAHistogramsList = new TList();
    fQAHistogramsList->SetName(fQAHistogramsListName);
    fQAHistogramsList->SetOwner(kTRUE);
    fHistList->Add(fQAHistogramsList);

    // centrality correlation QA histograms
    fCenCorQAHistogramsList = new TList();
    fCenCorQAHistogramsList->SetName(fCenCorQAHistogramsListName);
    fCenCorQAHistogramsList->SetOwner(kTRUE);
    fQAHistogramsList->Add(fCenCorQAHistogramsList);

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

  // 2. Book and nest lists for control histograms:
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

  // 3. Book and nest lists for final results:
  fFinalResultsList = new TList();
  fFinalResultsList->SetName(fFinalResultsListName);
  fFinalResultsList->SetOwner(kTRUE);
  fHistList->Add(fFinalResultsList);

  // 4. Book and nest lists for MC Analsysis
  if (fMCAnalaysis) {
    fMCAnalysisList = new TList();
    fMCAnalysisList->SetName(fMCAnalysisListName);
    fMCAnalysisList->SetOwner(kTRUE);
    fHistList->Add(fMCAnalysisList);
  }
}

void AliAnalysisTaskAR::BookQAHistograms() {
  // Book all QA histograms

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
      fCenCorQAHistograms[cen][ba]->SetStats(kFALSE);
      fCenCorQAHistograms[cen][ba]->SetOption("colz");
      fCenCorQAHistograms[cen][ba]->GetXaxis()->SetTitle(
          fCenCorQAHistogramNames[cen][ba][kXAXIS]);
      fCenCorQAHistograms[cen][ba]->GetYaxis()->SetTitle(
          fCenCorQAHistogramNames[cen][ba][kYAXIS]);
      fCenCorQAHistogramsList->Add(fCenCorQAHistograms[cen][ba]);
    }
  }

  fFBScanQAHistogram =
      new TH1D(fFBScanQAHistogramName[kNAME], fFBScanQAHistogramName[kTITLE],
               fFBScanQAHistogramBin[kBIN], fFBScanQAHistogramBin[kLEDGE],
               fFBScanQAHistogramBin[kUEDGE]);
  int fb = 1;
  for (int i = 0; i < kMaxFilterbit; ++i) {
    fFBScanQAHistogram->GetXaxis()->SetBinLabel(i + 1, Form("%d", fb));
    fb *= 2;
  }
  fFBScanQAHistogram->SetStats(kFALSE);
  fFBScanQAHistogram->SetFillColor(kFillColor[kAFTER]);
  fFBScanQAHistogramsList->Add(fFBScanQAHistogram);

  for (int track = 0; track < LAST_ETRACK; ++track) {
    for (int fb = 0; fb < kNumberofTestFilterBit; ++fb) {
      fFBTrackScanQAHistograms[track][fb] =
          new TH1D(fFBTrackScanQAHistogramNames[track][fb][kNAME],
                   fFBTrackScanQAHistogramNames[track][fb][kTITLE],
                   fFBTrackScanQAHistogramBins[track][kBIN],
                   fFBTrackScanQAHistogramBins[track][kLEDGE],
                   fFBTrackScanQAHistogramBins[track][kUEDGE]);
      fFBTrackScanQAHistograms[track][fb]->SetStats(kFALSE);
      fFBTrackScanQAHistograms[track][fb]->SetFillColor(kFillColor[kAFTER]);
      fFBTrackScanQAHistograms[track][fb]->GetXaxis()->SetTitle(
          fFBTrackScanQAHistogramNames[track][fb][kXAXIS]);
      fFBTrackScanQAHistograms[track][fb]->GetYaxis()->SetTitle(
          fFBTrackScanQAHistogramNames[track][fb][kYAXIS]);
      fFBScanQAHistogramsList->Add(fFBTrackScanQAHistograms[track][fb]);
    }
  }

  for (int var = 0; var < kKinematic; ++var) {
    for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
      fSelfCorQAHistograms[var][ba] =
          new TH1D(fSelfCorQAHistogramNames[var][ba][kNAME],
                   fSelfCorQAHistogramNames[var][ba][kTITLE],
                   fSelfCorQAHistogramBins[var][kBIN],
                   fSelfCorQAHistogramBins[var][kLEDGE],
                   fSelfCorQAHistogramBins[var][kUEDGE]);
      fSelfCorQAHistograms[var][ba]->SetStats(kFALSE);
      fSelfCorQAHistograms[var][ba]->SetFillColor(kFillColor[ba]);
      fSelfCorQAHistograms[var][ba]->SetMinimum(0.0);
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
  fTrackCutsCounter = new TH1D("fTrackCutsCounter", "Count track cuts",
                               LAST_ETRACK, 0, LAST_ETRACK);
  fTrackCutsCounter->SetStats(kFALSE);
  fTrackCutsCounter->SetFillColor(kFillColor[kAFTER]);
  for (int bin = 0; bin < fTrackCutsCounter->GetNbinsX(); ++bin) {
    fTrackCutsCounter->GetXaxis()->SetBinLabel(bin + 1,
                                               fTrackCutsCounterBinNames[bin]);
  }
  fTrackControlHistogramsList->Add(fTrackCutsCounter);
  // book track control histograms
  for (int var = 0; var < LAST_ETRACK; ++var) {
    for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
      fTrackControlHistograms[var][ba] =
          new TH1D(fTrackControlHistogramNames[var][ba][kNAME],
                   fTrackControlHistogramNames[var][ba][kTITLE],
                   fTrackControlHistogramBins[var][kBIN],
                   fTrackControlHistogramBins[var][kLEDGE],
                   fTrackControlHistogramBins[var][kUEDGE]);
      fTrackControlHistograms[var][ba]->SetStats(kFALSE);
      fTrackControlHistograms[var][ba]->SetFillColor(kFillColor[ba]);
      fTrackControlHistograms[var][ba]->SetMinimum(0.0);
      fTrackControlHistograms[var][ba]->GetXaxis()->SetTitle(
          fTrackControlHistogramNames[var][ba][kXAXIS]);
      fTrackControlHistograms[var][ba]->GetYaxis()->SetTitle(
          fTrackControlHistogramNames[var][ba][kYAXIS]);
      fTrackControlHistogramsList->Add(fTrackControlHistograms[var][ba]);
    }
  }

  // book histogram for counting event cuts
  fEventCutsCounter = new TH1D("fEventCutsCounter", "Count track cuts",
                               LAST_EEVENT, 0, LAST_EEVENT);
  fEventCutsCounter->SetStats(kFALSE);
  fEventCutsCounter->SetFillColor(kFillColor[kAFTER]);
  for (int bin = 0; bin < fEventCutsCounter->GetNbinsX(); ++bin) {
    fEventCutsCounter->GetXaxis()->SetBinLabel(bin + 1,
                                               fEventCutsCounterBinNames[bin]);
  }
  fEventControlHistogramsList->Add(fEventCutsCounter);
  // book event control histograms
  for (int var = 0; var < LAST_EEVENT; ++var) {
    for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
      fEventControlHistograms[var][ba] =
          new TH1D(fEventControlHistogramNames[var][ba][kNAME],
                   fEventControlHistogramNames[var][ba][kTITLE],
                   fEventControlHistogramBins[var][kBIN],
                   fEventControlHistogramBins[var][kLEDGE],
                   fEventControlHistogramBins[var][kUEDGE]);
      fEventControlHistograms[var][ba]->SetStats(kFALSE);
      fEventControlHistograms[var][ba]->SetFillColor(kFillColor[ba]);
      fEventControlHistograms[var][ba]->SetMinimum(0.0);
      fEventControlHistograms[var][ba]->GetXaxis()->SetTitle(
          fEventControlHistogramNames[var][ba][kXAXIS]);
      fEventControlHistograms[var][ba]->GetYaxis()->SetTitle(
          fEventControlHistogramNames[var][ba][kYAXIS]);
      fEventControlHistogramsList->Add(fEventControlHistograms[var][ba]);
    }
  }
}

void AliAnalysisTaskAR::BookFinalResultHistograms() {
  // Book all histograms to hold the final results

  Color_t colorFinalResult = kBlue - 10;

  // book event control histograms
  for (int var = 0; var < LAST_EFINALHIST; ++var) {
    fFinalResultHistograms[var] = new TH1D(
        fFinalResultHistogramNames[var][0], fFinalResultHistogramNames[var][1],
        fFinalResultHistogramBins[var][kBIN],
        fFinalResultHistogramBins[var][kLEDGE],
        fFinalResultHistogramBins[var][kUEDGE]);
    fFinalResultHistograms[var]->SetStats(kFALSE);
    fFinalResultHistograms[var]->SetFillColor(colorFinalResult);
    fFinalResultHistograms[var]->GetXaxis()->SetTitle(
        fFinalResultHistogramNames[var][2]);
    fFinalResultsList->Add(fFinalResultHistograms[var]);
  }
}

void AliAnalysisTaskAR::BookFinalResultProfiles() {
  // Book all profiles to hold the final results

  // book final result profiles
  for (int var = 0; var < LAST_EFINALPROFILE; ++var) {
    fFinalResultProfiles[var] = new TProfile(
        fFinalResultProfileNames[var][0], fFinalResultProfileNames[var][1],
        fFinalResultProfileBins[var][kBIN],
        fFinalResultProfileBins[var][kLEDGE],
        fFinalResultProfileBins[var][kUEDGE], nullptr);
    fFinalResultProfiles[var]->SetStats(kFALSE);
    fFinalResultProfiles[var]->GetXaxis()->SetTitle(
        fFinalResultProfileNames[var][2]);
    fFinalResultProfiles[var]->Sumw2();
    fFinalResultsList->Add(fFinalResultProfiles[var]);
  }
}

void AliAnalysisTaskAR::BookMCObjects() {
  // book objects need for MC analysis

  // protect at some point if fMCFlowHarmonics is empty
  if (!fMCFlowHarmonics) {
    std::cout << __LINE__ << ": no flow harmonics defined" << std::endl;
    Fatal("BookMCObjects", "Invalid Pointer");
  }

  // base setup for the pdf for MC analysis with flow harmonics
  // 1. generate formula, i.e. fourier series
  // 2. set flow harmonics as parameters as given by fMCFlowHarmonics
  // 3. leave symmetry planes and set them later on a event by event basis

  // generate formula
  TString Formula = "1+";
  for (int i = 1; i <= fMCFlowHarmonics->GetSize(); ++i) {
    Formula += Form("2*[%d]*TMath::Cos(%d*(x-[%d]))", 2 * i - 1, i, 2 * i);
    if (i < fMCFlowHarmonics->GetSize()) {
      Formula += "+";
    }
  }
  Formula = "(" + Formula + ")/TMath::TwoPi()";
  // create TF1 object
  fMCPdf = new TF1(fMCPdfName, Formula, 0., TMath::TwoPi());
  fMCAnalysisList->Add(fMCPdf);

  // set flow harmonics
  // flow harmonics are parameters with odd index
  for (int i = 0; i < fMCFlowHarmonics->GetSize(); ++i) {
    fMCPdf->SetParameter(2 * i + 1, fMCFlowHarmonics->GetAt(i));
  }
}

void AliAnalysisTaskAR::UserExec(Option_t *) {
  // if you do MC analysis locally, call MCOnTheFlyExec and bail out
  if (fMCAnalaysis) {
    MCOnTheFlyExec();
    return;
  }

  // general strategy for real data
  // 1. Reset event-by-event objects
  // 2. Get pointer to AOD event
  // 3. Start analysis over AODs, i.e. fill fPhi
  // 4. Fill event objects
  // 5. PostData

  // 1. Reset event-by-event objects
  fPhi.clear();
  fWeights.clear();

  // 2. Get pointer to AOD event
  AliAODEvent *aAOD = dynamic_cast<AliAODEvent *>(InputEvent()); // from TaskSE
  if (!aAOD) {
    return;
  }

  // fillhistograms before cut
  if (fFillQAHistograms) {
    FillEventQAHistograms(kBEFORE, aAOD);
  }
  FillEventControlHistograms(kBEFORE, aAOD);

  // cut event
  if (!SurviveEventCut(aAOD)) {
    return;
  }

  // fill histogram after event cut
  if (fFillQAHistograms) {
    FillEventQAHistograms(kAFTER, aAOD);
  }
  FillEventControlHistograms(kAFTER, aAOD);

  // 3. Start analysis over AODs:

  //  number of all tracks in current event
  Int_t nTracks = aAOD->GetNumberOfTracks();

  // count number of valid tracks before and after cutting for computing
  // multiplicity
  Int_t nTracks_beforeCut = 0;
  Int_t nTracks_afterCut = 0;

  // loop over all tracks in the event
  for (Int_t iTrack = 0; iTrack < nTracks; iTrack++) {

    // getting a pointer to a track
    AliAODTrack *aTrack = dynamic_cast<AliAODTrack *>(aAOD->GetTrack(iTrack));

    // protect against invalid pointers
    if (!aTrack) {
      continue;
    }

    // fill filter bit scan QA histogram
    if (fFillQAHistograms) {
      FillFBScanQAHistograms(aTrack);
    }

    // get kinematic variables of the track
    // Double_t pt = aTrack->Pt();
    Double_t phi = aTrack->Phi();
    // Double_t eta = aTrack->Eta();

    // fill track control histograms before track cut
    FillTrackControlHistograms(kBEFORE, aTrack);
    nTracks_beforeCut++;

    // cut track
    if (!SurviveTrackCut(aTrack)) {
      continue;
    }

    // fill track control histograms after track cut
    FillTrackControlHistograms(kAFTER, aTrack);
    nTracks_afterCut++;

    // finally, fill azimuthal angles into vector
    if (!fMCClosure) {
      fPhi.push_back(phi);
    } else {
      Int_t AcceptanceBin = fAcceptanceHistogram->FindBin(phi);
      if (gRandom->Uniform() <=
          fAcceptanceHistogram->GetBinContent(AcceptanceBin)) {
        fPhi.push_back(phi);
        fWeights.push_back(fWeightHistogram->GetBinContent(AcceptanceBin));
      }
    }
  }

  // fill control histogram for Multiplicity after counting all tracks
  // fEventControlHistograms[kMUL][kBEFORE]->Fill(nTracks_beforeCut);
  // fEventControlHistograms[kMUL][kAFTER]->Fill(nTracks_afterCut);

  // 4. Fill event objects
  // reset weights if required

  // calculate all qvectors
  CalculateQvectors();

  // fill final result profile
  FillFinalResultProfile(kHARDATA);

  // if option is given, repeat with weights resetted
  if (fResetWeights) {
    std::fill(fWeights.begin(), fWeights.end(), 1.);
    CalculateQvectors();
    FillFinalResultProfile(kHARDATARESET);
  }

  // d) PostData:
  PostData(1, fHistList);
}

void AliAnalysisTaskAR::MCOnTheFlyExec() {
  // call this method for local monte carlo analysis

  // reset angles and weights
  fPhi.clear();
  fWeights.clear();

  // set symmetry planes for MC analysis
  MCPdfSymmetryPlanesSetup();
  // loop over all particles in an event
  Double_t Phi = 0.0;
  Int_t AcceptanceBin = 0;
  for (int i = 0; i < GetMCNumberOfParticlesPerEvent(); ++i) {
    Phi = fMCPdf->GetRandom();

    if (fUseWeights) {
      AcceptanceBin = fAcceptanceHistogram->FindBin(Phi);
      if (gRandom->Uniform() <=
          fAcceptanceHistogram->GetBinContent(AcceptanceBin)) {
        fPhi.push_back(Phi);
        fWeights.push_back(fWeightHistogram->GetBinContent(AcceptanceBin));
      }
    } else {
      fPhi.push_back(Phi);
    }
  }

  // compute Q-vectors
  CalculateQvectors();
  // fill data into final result profile
  FillFinalResultProfile(kHARDATA);
  // if option is given, reset weight and recompute Q-vectors
  if (fResetWeights) {
    std::fill(fWeights.begin(), fWeights.end(), 1.);
    CalculateQvectors();
    FillFinalResultProfile(kHARDATARESET);
  }
}

void AliAnalysisTaskAR::FillEventControlHistograms(kBeforeAfter BA,
                                                   AliAODEvent *event) {
  // fill event control histograms

  // get centrality percentile
  Double_t centralityPercentile =
      dynamic_cast<AliMultSelection *>(event->FindListObject("MultSelection"))
          ->GetMultiplicityPercentile(fCentralitySelCriterion);

  // get primary vertex object
  AliAODVertex *PrimaryVertex = event->GetPrimaryVertex();

  // fill control histograms
  fEventControlHistograms[kMUL][BA]->Fill(event->GetNumberOfTracks());
  fEventControlHistograms[kCEN][BA]->Fill(centralityPercentile);
  fEventControlHistograms[kNCONTRIB][BA]->Fill(
      PrimaryVertex->GetNContributors());
  fEventControlHistograms[kX][BA]->Fill(PrimaryVertex->GetX());
  fEventControlHistograms[kY][BA]->Fill(PrimaryVertex->GetY());
  fEventControlHistograms[kZ][BA]->Fill(PrimaryVertex->GetZ());
}

void AliAnalysisTaskAR::FillTrackControlHistograms(kBeforeAfter BA,
                                                   AliAODTrack *track) {
  // fill track control histograms
  fTrackControlHistograms[kPT][BA]->Fill(track->Pt());
  fTrackControlHistograms[kPHI][BA]->Fill(track->Phi());
  fTrackControlHistograms[kETA][BA]->Fill(track->Eta());
  fTrackControlHistograms[kCHARGE][BA]->Fill(track->Charge());
  fTrackControlHistograms[kTPCNCLS][BA]->Fill(track->GetTPCNcls());
  fTrackControlHistograms[kITSNCLS][BA]->Fill(track->GetITSNcls());
  fTrackControlHistograms[kCHI2PERNDF][BA]->Fill(track->Chi2perNDF());
  fTrackControlHistograms[kDCAZ][BA]->Fill(track->ZAtDCA());
  fTrackControlHistograms[kDCAXY][BA]->Fill(track->DCA());
}

void AliAnalysisTaskAR::FillEventQAHistograms(kBeforeAfter BA,
                                              AliAODEvent *event) {
  // fill QA control histograms

  // get all centrality percentiles
  Double_t centralityPercentile[LAST_ECENESTIMATORS];

  centralityPercentile[kV0M] =
      dynamic_cast<AliMultSelection *>(event->FindListObject("MultSelection"))
          ->GetMultiplicityPercentile("V0M");
  centralityPercentile[kCL0] =
      dynamic_cast<AliMultSelection *>(event->FindListObject("MultSelection"))
          ->GetMultiplicityPercentile("CL0");
  centralityPercentile[kCL1] =
      dynamic_cast<AliMultSelection *>(event->FindListObject("MultSelection"))
          ->GetMultiplicityPercentile("CL1");
  centralityPercentile[kSPDTRACKLETS] =
      dynamic_cast<AliMultSelection *>(event->FindListObject("MultSelection"))
          ->GetMultiplicityPercentile("SPDTracklets");

  fCenCorQAHistograms[0][BA]->Fill(centralityPercentile[kV0M],
                                   centralityPercentile[kCL0]);
  fCenCorQAHistograms[1][BA]->Fill(centralityPercentile[kV0M],
                                   centralityPercentile[kCL1]);
  fCenCorQAHistograms[2][BA]->Fill(centralityPercentile[kV0M],
                                   centralityPercentile[kSPDTRACKLETS]);
  fCenCorQAHistograms[3][BA]->Fill(centralityPercentile[kV0M],
                                   centralityPercentile[kCL0]);
  fCenCorQAHistograms[4][BA]->Fill(centralityPercentile[kV0M],
                                   centralityPercentile[kCL1]);
  fCenCorQAHistograms[5][BA]->Fill(centralityPercentile[kV0M],
                                   centralityPercentile[kCL1]);

  // search for self correlations with nested loop
  Int_t nTracks = event->GetNumberOfTracks();
  // starting a loop over the first track
  for (Int_t iTrack1 = 0; iTrack1 < nTracks; iTrack1++) {
    AliAODTrack *aTrack1 =
        dynamic_cast<AliAODTrack *>(event->GetTrack(iTrack1));
    if (!aTrack1 || !SurviveTrackCut(aTrack1)) {
      continue;
    }
    // starting a loop over the second track
    for (Int_t iTrack2 = iTrack1 + 1; iTrack2 < nTracks; iTrack2++) {
      AliAODTrack *aTrack2 =
          dynamic_cast<AliAODTrack *>(event->GetTrack(iTrack2));
      if (!aTrack2 || !SurviveTrackCut(aTrack2)) {
        continue;
      }
      fSelfCorQAHistograms[kPHI][BA]->Fill(aTrack1->Phi() - aTrack2->Phi());
      fSelfCorQAHistograms[kPT][BA]->Fill(aTrack1->Pt() - aTrack2->Pt());
      fSelfCorQAHistograms[kETA][BA]->Fill(aTrack1->Eta() - aTrack2->Eta());
    }
  }
}

void AliAnalysisTaskAR::FillFBScanQAHistograms(AliAODTrack *track) {
  // fill track QA histograms
  int fb = 1;
  for (int i = 0; i < kMaxFilterbit; ++i) {
    if (track->TestFilterBit(fb)) {
      fFBScanQAHistogram->Fill(i);
    }
    fb *= 2;
  }
  for (int fb = 0; fb < kNumberofTestFilterBit; ++fb) {
    if (track->TestFilterBit(kTestFilterbit[fb])) {
      fFBTrackScanQAHistograms[kPT][fb]->Fill(track->Pt());
      fFBTrackScanQAHistograms[kPHI][fb]->Fill(track->Phi());
      fFBTrackScanQAHistograms[kETA][fb]->Fill(track->Eta());
      fFBTrackScanQAHistograms[kCHARGE][fb]->Fill(track->Charge());
      fFBTrackScanQAHistograms[kTPCNCLS][fb]->Fill(track->GetTPCNcls());
      fFBTrackScanQAHistograms[kITSNCLS][fb]->Fill(track->GetITSNcls());
      fFBTrackScanQAHistograms[kCHI2PERNDF][fb]->Fill(track->Chi2perNDF());
      fFBTrackScanQAHistograms[kDCAZ][fb]->Fill(track->ZAtDCA());
      fFBTrackScanQAHistograms[kDCAXY][fb]->Fill(track->DCA());
    }
  }
}

Bool_t AliAnalysisTaskAR::SurviveEventCut(AliVEvent *ave) {

  // Check if the current event survives event cuts

  // Determine Ali{MC,ESD,AOD}Event:
  // AliMCEvent *aMC = dynamic_cast<AliMCEvent*>(ave);
  // AliESDEvent *aESD = dynamic_cast<AliESDEvent*>(ave);
  // get object for determining centrality

  // return flag at the event, if one cut is not passed, set it to false
  Bool_t Flag = kTRUE;

  // cast into AOD event
  AliAODEvent *aAOD = dynamic_cast<AliAODEvent *>(ave);
  if (!aAOD) {
    return kFALSE;
  }

  // get centrality percentile
  AliMultSelection *ams =
      (AliMultSelection *)aAOD->FindListObject("MultSelection");
  if (!ams) {
    return kFALSE;
  }
  Double_t MultiplicityPercentile =
      ams->GetMultiplicityPercentile(fCentralitySelCriterion);

  // cut event if it is not within the centrality percentile
  if ((MultiplicityPercentile < fEventCuts[kCEN][kMIN]) ||
      (MultiplicityPercentile > fEventCuts[kCEN][kMAX])) {
    fEventCutsCounter->Fill(kCEN + 0.5);
    Flag = kFALSE;
  }

  // Get primary vertex
  AliAODVertex *PrimaryVertex = aAOD->GetPrimaryVertex();
  if (!PrimaryVertex) {
    return kFALSE;
  }

  // cut event if primary vertex is too out of center
  if ((PrimaryVertex->GetX() < fEventCuts[kX][kMIN]) ||
      (PrimaryVertex->GetX() > fEventCuts[kX][kMAX])) {
    fEventCutsCounter->Fill(kX + 0.5);
    Flag = kFALSE;
  }
  if ((PrimaryVertex->GetY() < fEventCuts[kY][kMIN]) ||
      (PrimaryVertex->GetY() > fEventCuts[kY][kMAX])) {
    fEventCutsCounter->Fill(kY + 0.5);
    Flag = kFALSE;
  }
  if ((PrimaryVertex->GetZ() < fEventCuts[kZ][kMIN]) ||
      (PrimaryVertex->GetZ() > fEventCuts[kZ][kMAX])) {
    fEventCutsCounter->Fill(kZ + 0.5);
    Flag = kFALSE;
  }

  // cut on multiplicity
  if ((aAOD->GetNumberOfTracks() < fEventCuts[kMUL][kMIN]) ||
      (aAOD->GetNumberOfTracks() > fEventCuts[kMUL][kMAX])) {
    fEventCutsCounter->Fill(kMUL + 0.5);
    Flag = kFALSE;
  }

  // cut on numbers of contriubters to the vertex
  if ((PrimaryVertex->GetNContributors() < fEventCuts[kNCONTRIB][kMIN]) ||
      (PrimaryVertex->GetNContributors() > fEventCuts[kNCONTRIB][kMAX])) {
    fEventCutsCounter->Fill(kNCONTRIB + 0.5);
    Flag = kFALSE;
  }

  return Flag;
}

Bool_t AliAnalysisTaskAR::SurviveTrackCut(AliAODTrack *aTrack) {
  // check if current track survives track cut

  // return flag at the end, if one cut fails, set it to false
  Bool_t Flag = kTRUE;

  // if set, cut all non-primary particles away
  if (fPrimaryOnly) {
    if (aTrack->GetType() != AliAODTrack::kPrimary) {
      Flag = kFALSE;
    }
  }
  // cut PT
  if ((aTrack->Pt() < fTrackCuts[kPT][kMIN]) ||
      (aTrack->Pt() > fTrackCuts[kPT][kMAX])) {
    fTrackCutsCounter->Fill(kPT + 0.5);
    Flag = kFALSE;
  }
  // cut PHI
  if ((aTrack->Phi() < fTrackCuts[kPHI][kMIN]) ||
      (aTrack->Phi() > fTrackCuts[kPHI][kMAX])) {
    fTrackCutsCounter->Fill(kPHI + 0.5);
    Flag = kFALSE;
  }
  // cut ETA
  if ((aTrack->Eta() < fTrackCuts[kETA][kMIN]) ||
      (aTrack->Eta() > fTrackCuts[kETA][kMAX])) {
    fTrackCutsCounter->Fill(kETA + 0.5);
    Flag = kFALSE;
  }
  // cut on number of clusters in the TPC
  if ((aTrack->GetTPCNcls() < fTrackCuts[kTPCNCLS][kMIN]) ||
      (aTrack->GetTPCNcls() > fTrackCuts[kTPCNCLS][kMAX])) {
    fTrackCutsCounter->Fill(kTPCNCLS + 0.5);
    Flag = kFALSE;
  }
  // cut on number of clusters in the ITS
  if ((aTrack->GetITSNcls() < fTrackCuts[kITSNCLS][kMIN]) ||
      (aTrack->GetITSNcls() > fTrackCuts[kITSNCLS][kMAX])) {
    fTrackCutsCounter->Fill(kITSNCLS + 0.5);
    Flag = kFALSE;
  }
  // cut on chi2 / NDF of the track fit
  if ((aTrack->Chi2perNDF() < fTrackCuts[kCHI2PERNDF][kMIN]) ||
      (aTrack->Chi2perNDF() > fTrackCuts[kCHI2PERNDF][kMAX])) {
    fTrackCutsCounter->Fill(kCHI2PERNDF + 0.5);
    Flag = kFALSE;
  }
  // cut DCA in z direction
  if ((aTrack->ZAtDCA() < fTrackCuts[kDCAZ][kMIN]) ||
      (aTrack->ZAtDCA() > fTrackCuts[kDCAZ][kMAX])) {
    fTrackCutsCounter->Fill(kDCAZ + 0.5);
    Flag = kFALSE;
  }
  // cut DCA in xy plane
  if ((aTrack->DCA() < fTrackCuts[kDCAXY][kMIN]) ||
      (aTrack->DCA() > fTrackCuts[kDCAXY][kMAX])) {
    fTrackCutsCounter->Fill(kDCAXY + 0.5);
    Flag = kFALSE;
  }

  // cut with filtertbit
  // filter bit 128 denotes TPC-only tracks, use only them for the
  // analysis
  // for hybrid tracks use filterbit 782
  // for more information about filterbits see the online week
  // the filterbits can change from run to run
  // fill control histograms
  if (!aTrack->TestFilterBit(fFilterbit)) {
    Flag = kFALSE;
  }
  return Flag;
}

void AliAnalysisTaskAR::MCPdfSymmetryPlanesSetup() {
  // set symmetry planes randomly on a event by event basis
  // Double_t Psi = 0;
  Double_t Psi = gRandom->Uniform(0., TMath::TwoPi());
  for (int i = 0; i < fMCFlowHarmonics->GetSize(); ++i) {
    fMCPdf->SetParameter(2 * (i + 1), Psi);
  }
}

Int_t AliAnalysisTaskAR::GetMCNumberOfParticlesPerEvent() {

  if (!fMCNumberOfParticlesPerEventFluctuations) {
    return fMCNumberOfParticlesPerEvent;
  }

  return gRandom->Uniform(fMCNumberOfParticlesPerEventRange[kMIN],
                          fMCNumberOfParticlesPerEventRange[kMAX]);
};

void AliAnalysisTaskAR::CalculateQvectors() {
  // Calculate Q-vectors

  // 1) Make sure all Q-vectors are initially zero
  for (Int_t h = 0; h < kMaxHarmonic; h++) {
    for (Int_t p = 0; p < kMaxPower; p++) {
      fQvector[h][p] = TComplex(0., 0.);
    }
  }

  // 2) Calculate Q-vectors for available angles and weights:
  Double_t dPhi = 0.;
  Double_t wPhi = 1.;         // particle weight
  Double_t wPhiToPowerP = 1.; // particle weight raised to power p
  for (std::size_t i = 0; i < fPhi.size(); i++) {
    dPhi = fPhi[i];
    if (fUseWeights) {
      wPhi = fWeights[i];
    }
    for (Int_t h = 0; h < kMaxHarmonic; h++) {
      for (Int_t p = 0; p < kMaxPower; p++) {
        if (fUseWeights) {
          wPhiToPowerP = pow(wPhi, p);
        }
        fQvector[h][p] += TComplex(wPhiToPowerP * TMath::Cos(h * dPhi),
                                   wPhiToPowerP * TMath::Sin(h * dPhi));
      }
    }
  }
}

void AliAnalysisTaskAR::FillFinalResultProfile(kFinalProfile fp) {
  // fill final result profiles

  Double_t corr = 0.0;
  Double_t weight = 0.0;
  Int_t index = 0;

  // loop over all correlators
  for (auto Corr : fCorrelators) {
    // protect against insufficient amount of statistics
    // i.e. number of paritcles is lower then the order of correlator due to
    // track cuts
    if (fPhi.size() < Corr.size()) {
      return;
    }

    // compute correlator
    switch (static_cast<int>(Corr.size())) {
    case 2:
      corr = Two(Corr.at(0), Corr.at(1)).Re();
      weight = Two(0, 0).Re();
      break;
    case 3:
      corr = Three(Corr.at(0), Corr.at(1), Corr.at(2)).Re();
      weight = Three(0, 0, 0).Re();
      break;
    case 4:
      corr = Four(Corr.at(0), Corr.at(1), Corr.at(2), Corr.at(3)).Re();
      weight = Four(0, 0, 0, 0).Re();
      break;
    case 5:
      corr =
          Five(Corr.at(0), Corr.at(1), Corr.at(2), Corr.at(3), Corr.at(4)).Re();
      weight = Five(0, 0, 0, 0, 0).Re();
      break;
    case 6:
      corr = Six(Corr.at(0), Corr.at(1), Corr.at(2), Corr.at(3), Corr.at(4),
                 Corr.at(5))
                 .Re();
      weight = Six(0, 0, 0, 0, 0, 0).Re();
      break;
    default:
      corr = Recursion(Corr.size(), Corr.data()).Re();
      weight = Recursion(Corr.size(), std::vector<Int_t>(Corr.size(), 0).data())
                   .Re();
    }

    index++;
    // fill final resutl profile
    fFinalResultProfiles[fp]->Fill(index - 0.5, corr / weight, weight);
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

Double_t AliAnalysisTaskAR::CombinatorialWeight(Int_t n) {
  // calculate combinatrial weight for Qvectors
  if (n >= static_cast<Int_t>(fPhi.size())) {
    std::cout << __LINE__ << ": Two few particles for this correlator"
              << std::endl;
    Fatal("Combinatorial weight",
          "order of correlator is larger then number of particles");
  }
  Double_t w = 1.0;
  for (int i = 0; i < n; ++i) {
    w *= (fPhi.size() - i);
  }
  return w;
}

TComplex AliAnalysisTaskAR::TwoNestedLoops(Int_t n1, Int_t n2) {
  // Calculation of <cos(n1*phi1+n2*phi2)> and <sin(n1*phi1+n2*phi2)>
  // with two nested loops

  TComplex Two(0., 0.);

  Double_t phi1 = 0., phi2 = 0.; // particle angle
  Double_t w1 = 1., w2 = 1.;     // particle weight
  for (std::size_t i1 = 0; i1 < fPhi.size(); i1++) {
    phi1 = fPhi[i1];
    if (fUseWeights) {
      w1 = fWeights[i1];
    }
    for (std::size_t i2 = 0; i2 < fPhi.size(); i2++) {
      if (i2 == i1) {
        continue;
      } // Get rid of autocorrelations
      phi2 = fPhi[i2];
      if (fUseWeights) {
        w2 = fWeights[i2];
      }
      Two += TComplex(TMath::Cos(n1 * w1 * phi1 + n2 * w2 * phi2),
                      TMath::Sin(n1 * w1 * phi1 + n2 * w2 * phi2));
    }
  }
  return Two / CombinatorialWeight(2);
}

TComplex AliAnalysisTaskAR::ThreeNestedLoops(Int_t n1, Int_t n2, Int_t n3) {
  // Calculation of <cos(n1*phi1+n2*phi2+n3*phi3)> and
  // <sin(n1*phi1+n2*phi2+n3*phi3)> with three nested loops.

  TComplex Q(0., 0.);
  Double_t phi1 = 0., phi2 = 0., phi3 = 0.; // particle angle
  Double_t w1 = 1., w2 = 1., w3 = 1.;       // particle weight
  for (std::size_t i1 = 0; i1 < fPhi.size(); i1++) {
    phi1 = fPhi[i1];
    if (fUseWeights) {
      w1 = fWeights[i1];
    }
    for (std::size_t i2 = 0; i2 < fPhi.size(); i2++) {
      if (i2 == i1) {
        continue;
      } // Get rid of autocorrelations
      phi2 = fPhi[i2];
      if (fUseWeights) {
        w2 = fWeights[i2];
      }
      for (std::size_t i3 = 0; i3 < fPhi.size(); i3++) {
        if (i3 == i1 || i3 == i2) {
          continue;
        } // Get rid of autocorrelations
        phi3 = fPhi[i3];
        if (fUseWeights) {
          w3 = fWeights[i3];
        }
        Q += TComplex(
            TMath::Cos(n1 * w1 * phi1 + n2 * w2 * phi2 + n3 * w3 * phi3),
            TMath::Sin(n1 * w1 * phi1 + n2 * w2 * phi2 + n3 * w3 * phi3));
      }
    }
  }
  return Q / CombinatorialWeight(3);
}

TComplex AliAnalysisTaskAR::FourNestedLoops(Int_t n1, Int_t n2, Int_t n3,
                                            Int_t n4) {
  // Calculation of <cos(n1*phi1+n2*phi2+n3*phi3+n4*phi4)> and
  // <sin(n1*phi1+n2*phi2+n3*phi3+n4*phi4)> with four nested loops.

  TComplex Q(0., 0.);
  Double_t phi1 = 0., phi2 = 0., phi3 = 0., phi4 = 0.; // particle angle
  Double_t w1 = 1., w2 = 1., w3 = 1., w4 = 1.;         // particle weight
  for (std::size_t i1 = 0; i1 < fPhi.size(); i1++) {
    phi1 = fPhi[i1];
    if (fUseWeights) {
      w1 = fWeights[i1];
    }
    for (std::size_t i2 = 0; i2 < fPhi.size(); i2++) {
      if (i2 == i1) {
        continue;
      } // Get rid of autocorrelations
      phi2 = fPhi[i2];
      if (fUseWeights) {
        w2 = fWeights[i2];
      }
      for (std::size_t i3 = 0; i3 < fPhi.size(); i3++) {
        if (i3 == i1 || i3 == i2) {
          continue;
        } // Get rid of autocorrelations
        phi3 = fPhi[i3];
        if (fUseWeights) {
          w3 = fWeights[i3];
        }
        for (std::size_t i4 = 0; i4 < fPhi.size(); i4++) {
          if (i4 == i1 || i4 == i2 || i4 == i3) {
            continue;
          } // Get rid of autocorrelations
          phi4 = fPhi[i4];
          if (fUseWeights) {
            w4 = fWeights[i4];
          }
          Q += TComplex(TMath::Cos(n1 * w1 * phi1 + n2 * w2 * phi2 +
                                   n3 * w3 * phi3 + n4 * w4 * phi4),
                        TMath::Sin(n1 * w1 * phi1 + n2 * w2 * phi2 +
                                   n3 * w3 * phi3 + n4 * w4 * phi4));
        }
      }
    }
  }
  return Q / CombinatorialWeight(4);
}

void AliAnalysisTaskAR::SetAcceptanceHistogram(const char *Filename,
                                               const char *Histname) {
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
  this->fAcceptanceHistogram = dynamic_cast<TH1D *>(file->Get(Histname));
  if (!fAcceptanceHistogram) {
    std::cout << __LINE__ << ": No acceptance histogram" << std::endl;
    Fatal("SetAcceptanceHistogram", "Cannot get acceptance histogram");
  }
  // keeps the histogram in memory after we close the file
  this->fAcceptanceHistogram->SetDirectory(0);
  file->Close();
}

void AliAnalysisTaskAR::SetWeightHistogram(const char *Filename,
                                           const char *Histname) {
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
  this->fWeightHistogram = dynamic_cast<TH1D *>(file->Get(Histname));
  if (!fWeightHistogram) {
    std::cout << __LINE__ << ": No acceptance histogram" << std::endl;
    Fatal("SetWeightHistogram", "Cannot get weight histogram");
  }
  // keeps the histogram in memory after we close the file
  this->fWeightHistogram->SetDirectory(0);
  file->Close();
}

void AliAnalysisTaskAR::GetPointers(TList *histList) {
  // Initialize pointer for base list fHistList so we can initialize all ot
  // bjects and call terminate off-lin

  fHistList = histList;
  if (!fHistList) {
    std::cout << __LINE__ << ": Did not get " << fHistListName << std::endl;
    Fatal("GetPointers", "Invalid Pointer");
  }

  // initialize all other objects
  this->GetPointersForControlHistograms();
  this->GetPointersForFinalResultHistograms();
  this->GetPointersForFinalResultProfiles();
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
      fHistList->FindObject(fTrackControlHistogramsListName));
  if (!fTrackControlHistogramsList) {
    std::cout << __LINE__ << ": Did not get " << fTrackControlHistogramsListName
              << std::endl;
    Fatal("GetPointersForControlHistograms", "Invalid Pointer");
  }

  // get all pointers for track control histograms
  for (int var = 0; var < LAST_ETRACK; ++var) {
    for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
      fTrackControlHistograms[var][ba] =
          dynamic_cast<TH1D *>(fTrackControlHistogramsList->FindObject(
              fTrackControlHistogramNames[var][ba][0]));
      if (!fTrackControlHistograms[var][ba]) {
        std::cout << __LINE__ << ": Did not get "
                  << fTrackControlHistogramNames[var][ba][0] << std::endl;
        Fatal("GetPointersForControlHistograms", "Invalid Pointer");
      }
    }
  }

  // get pointer for fEventControlHistogramsList
  fEventControlHistogramsList = dynamic_cast<TList *>(
      fHistList->FindObject(fEventControlHistogramsListName));
  if (!fEventControlHistogramsList) {
    std::cout << __LINE__ << ": Did not get " << fEventControlHistogramsListName
              << std::endl;
    Fatal("GetPointersForControlHistograms", "Invalid Pointer");
  }

  // get all pointers for event control histograms
  for (int var = 0; var < LAST_EEVENT; ++var) {
    for (int ba = 0; ba < LAST_EBEFOREAFTER; ++ba) {
      fEventControlHistograms[var][ba] =
          dynamic_cast<TH1D *>(fEventControlHistogramsList->FindObject(
              fEventControlHistogramNames[var][ba][0]));
      if (!fEventControlHistograms[var][ba]) {
        std::cout << __LINE__ << ": Did not get "
                  << fEventControlHistogramNames[var][ba][0] << std::endl;
        Fatal("GetPointersForControlHistograms", "Invalid Pointer");
      }
    }
  }
}

void AliAnalysisTaskAR::GetPointersForFinalResultHistograms() {
  // Get pointers for final result Histograms

  // Get pointer for fFinalResultsList
  fFinalResultsList =
      dynamic_cast<TList *>(fHistList->FindObject(fFinalResultsListName));
  if (!fFinalResultsList) {
    std::cout << __LINE__ << ": Did not get " << fFinalResultsListName
              << std::endl;
    Fatal("GetPointersForOutputHistograms", "Invalid Pointer");
  }

  // get all pointers for final result histograms
  for (int var = 0; var < LAST_EFINALHIST; ++var) {
    fFinalResultHistograms[var] = dynamic_cast<TH1D *>(
        fFinalResultsList->FindObject(fFinalResultHistogramNames[var][0]));
    if (!fFinalResultHistograms[var]) {
      std::cout << __LINE__ << ": Did not get "
                << fFinalResultHistogramNames[var][0] << std::endl;
      Fatal("GetPointersForOutputHistograms", "Invalid Pointer");
    }
  }

  // Set again all flags:
  // fFillBuffers = (Bool_t)fBuffersFlagsPro->GetBinContent(1);
  // fMaxBuffer = fBuffersFlagsPro->GetBinContent(2);
}

void AliAnalysisTaskAR::GetPointersForFinalResultProfiles() {
  // Get pointers for final result Histograms

  // Get pointer for fFinalResultsList
  fFinalResultsList =
      dynamic_cast<TList *>(fHistList->FindObject(fFinalResultsListName));
  if (!fFinalResultsList) {
    std::cout << __LINE__ << ": Did not get " << fFinalResultsListName
              << std::endl;
    Fatal("GetPointersForOutputHistograms", "Invalid Pointer");
  }

  // get all pointers for final result histograms
  for (int var = 0; var < LAST_EFINALPROFILE; ++var) {
    fFinalResultProfiles[var] = dynamic_cast<TProfile *>(
        fFinalResultsList->FindObject(fFinalResultProfileNames[var][0]));
    if (!fFinalResultProfiles[var]) {
      std::cout << __LINE__ << ": Did not get "
                << fFinalResultProfileNames[var][0] << std::endl;
      Fatal("GetPointersForOutputProfiles", "Invalid Pointer");
    }
  }

  // Set again all flags:
  // fFillBuffers = (Bool_t)fBuffersFlagsPro->GetBinContent(1);
  // fMaxBuffer = fBuffersFlagsPro->GetBinContent(2);
}
