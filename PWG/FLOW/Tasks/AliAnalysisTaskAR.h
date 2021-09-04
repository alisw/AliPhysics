/**
 * File              : AliAnalysisTaskAR.h
 * Author            : Anton Riedel <anton.riedel@tum.de>
 * Date              : 07.05.2021
 * Last Modified Date: 04.09.2021
 * Last Modified By  : Anton Riedel <anton.riedel@tum.de>
 */

/*
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.
 * See cxx source for full Copyright notice
 * $Id$
 */

#ifndef ALIANALYSISTASKAR_H
#define ALIANALYSISTASKAR_H

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAnalysisTaskSE.h"
#include "AliVEvent.h"
#include <Riostream.h>
#include <TComplex.h>
#include <TDataType.h>
#include <TExMap.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TString.h>
#include <vector>

// global constants for qvectors
const Int_t kMaxHarmonic = 20;
const Int_t kMaxCorrelator = 20;
const Int_t kMaxPower = 20;
// global constants for QA filterbit scan
const Int_t kMaxFilterbit = 15; // 2^(15-1)=16384
const Int_t kNumberofTestFilterBit = 5;
const Int_t kTestFilterbit[5] = {1, 128, 256, 512, 768};
// centrality estimators
enum kCenEstimators { kV0M, kCL0, kCL1, kSPDTRACKLETS, LAST_ECENESTIMATORS };
const TString kCenEstimatorNames[LAST_ECENESTIMATORS] = {"V0M", "CL0", "CL1",
                                                         "SPDTracklets"};
// event variables
enum kEvent {
  kMUL,
  kMULQ,
  kMULW,
  kMULREF,
  kNCONTRIB,
  kCEN,
  kX,
  kY,
  kZ,
  kVPOS,
  LAST_EEVENT
};
// multiplicity estimators
const Int_t kMulEstimators = kNCONTRIB + 1;
const TString kMulEstimatorNames[kMulEstimators] = {"kMUL", "kMULQ", "kMULW",
                                                    "kMULREF", "kNCONTRIB"};
// track variables
enum kTrack {
  kPT,
  kPHI,
  kETA,
  kCHARGE,
  kTPCNCLS,
  kITSNCLS,
  kCHI2PERNDF,
  kDCAZ,
  kDCAXY,
  LAST_ETRACK
};
// kinematic variables
const Int_t kKinematic = kETA + 1;
// various gloabl objects
enum kFinalHist { kPHIAVG, LAST_EFINALHIST };
enum kFinalProfile { kHARDATA, kHARDATARESET, kHARTHEO, LAST_EFINALPROFILE };
enum kBins { kBIN, kLEDGE, kUEDGE, LAST_EBINS };
enum kName { kNAME, kTITLE, kXAXIS, kYAXIS, LAST_ENAME };
enum kMinMax { kMIN, kMAX, LAST_EMINMAX };
const TString kMMName[LAST_EMINMAX] = {"[kMIN]", "[kMAX]"};
enum kBeforeAfter { kBEFORE, kAFTER, LAST_EBEFOREAFTER };
const TString kBAName[LAST_EBEFOREAFTER] = {"[kBEFORE]", "[kAFTER]"};
const Color_t kFillColor[LAST_EBEFOREAFTER] = {kRed - 10, kGreen - 10};
enum kMode { kRECO, kSIM, LAST_EMODE };
const TString kModeName[LAST_EMODE] = {"[kRECO]", "[kSIM]"};

class AliAnalysisTaskAR : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskAR();
  AliAnalysisTaskAR(const char *name, Bool_t useParticleWeights = kFALSE);
  virtual ~AliAnalysisTaskAR();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);

  // methods called in the constructor
  virtual void InitializeArrays();
  virtual void InitializeArraysForQAHistograms();
  virtual void InitializeArraysForTrackControlHistograms();
  virtual void InitializeArraysForEventControlHistograms();
  virtual void InitializeArraysForCuts();
  virtual void InitializeArraysForWeights();
  virtual void InitializeArraysForQvectors();
  virtual void InitializeArraysForFinalResultHistograms();
  virtual void InitializeArraysForFinalResultProfiles();
  virtual void InitializeArraysForMCAnalysis();

  // methods called in UserCreateOutputObjects()
  virtual void BookAndNestAllLists();
  virtual void BookQAHistograms();
  virtual void BookControlHistograms();
  virtual void BookFinalResultHistograms();
  virtual void BookFinalResultProfiles();
  virtual void BookMCOnTheFlyObjects();

  // functions called in UserExec()
  virtual void MCOnTheFlyExec();
  virtual void FillEventQAHistograms(kBeforeAfter BA, AliVEvent *event);
  virtual void FillFBScanQAHistograms(AliAODTrack *track);
  virtual void FillEventControlHistograms(kBeforeAfter BA, AliVEvent *event);
  virtual void FillTrackControlHistograms(kBeforeAfter BA, AliVParticle *avp);
  virtual void FillFinalResultProfile(kFinalProfile fp);
  virtual Bool_t SurviveEventCut(AliVEvent *ave);
  virtual Bool_t SurviveTrackCut(AliVParticle *aTrack, Bool_t FillCounter);
  virtual void FillEventObjects(AliAODEvent *aAOD, AliMCEvent *aMC);
  virtual void ClearVectors();
  virtual void FillTrackObjects(AliAODTrack *track);
  virtual void AggregateWeights();
  virtual void ResetWeights();
  virtual Int_t IndexCorHistograms(Int_t i, Int_t j, Int_t N);

  // methods called in MCOnTheFlyExec()
  virtual void MCPdfSymmetryPlanesSetup();
  virtual Int_t GetMCNumberOfParticlesPerEvent();

  // methods for computing qvectors
  void CalculateQvectors();
  TComplex Q(Int_t n, Int_t p);
  TComplex Two(Int_t n1, Int_t n2);
  TComplex Three(Int_t n1, Int_t n2, Int_t n3);
  TComplex Four(Int_t n1, Int_t n2, Int_t n3, Int_t n4);
  TComplex Five(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5);
  TComplex Six(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6);
  TComplex Recursion(Int_t n, Int_t *harmonic, Int_t mult = 1, Int_t skip = 0);
  TComplex TwoNestedLoops(Int_t n1, Int_t n2);
  TComplex ThreeNestedLoops(Int_t n1, Int_t n2, Int_t n3);
  TComplex FourNestedLoops(Int_t n1, Int_t n2, Int_t n3, Int_t n4);
  // TComplex FiveNestedLoops(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5);
  // TComplex SixNestedLoops(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5,
  // Int_t n6);
  Double_t CombinatorialWeight(Int_t n);

  // GetPointers Methods in case we need to manually trigger Terminate()
  virtual void GetPointers(TList *list);
  virtual void GetPointersForQAHistograms();
  virtual void GetPointersForControlHistograms();
  virtual void GetPointersForFinalResultHistograms();
  virtual void GetPointersForFinalResultProfiles();

  // setters and getters for list objects
  void SetControlHistogramsList(TList *const chl) {
    this->fControlHistogramsList = chl;
  };
  TList *GetControlHistogramsList() const {
    return this->fControlHistogramsList;
  }
  void SetFinalResultsList(TList *const frl) { this->fFinalResultsList = frl; };
  TList *GetFinalResultsList() const { return this->fFinalResultsList; }

  // setters for QA histograms
  void SetFillQAHistograms(Bool_t option) { fFillQAHistograms = option; }
  // generic setter for centrality correlation QA histogram binning
  void SetCenCorQAHistogramBinning(Int_t cen1, Int_t xnbins,
                                   Double_t xlowerEdge, Double_t xupperEdge,
                                   Int_t cen2, Int_t ynbins,
                                   Double_t ylowerEdge, Double_t yupperEdge);
  // generic setter for multiplicity correlation QA histogram binning
  void SetMulCorQAHistogramBinning(Int_t mul1, Int_t xnbins,
                                   Double_t xlowerEdge, Double_t xupperEdge,
                                   Int_t mul2, Int_t ynbins,
                                   Double_t ylowerEdge, Double_t yupperEdge);
  // generic setter for track scan QA histograms
  void SetFBTrackScanQAHistogramBinning(kTrack Track, Int_t nbins,
                                        Double_t lowerEdge,
                                        Double_t upperEdge) {
    if (Track >= LAST_ETRACK) {
      std::cout << __LINE__ << ": running out of bounds" << std::endl;
      Fatal("SetFBTrackScanQAHistogramBinning",
            "Running out of bounds in SetFBTrackScanQAHistogramBinning");
    }
    if (upperEdge < lowerEdge) {
      std::cout << __LINE__
                << ": upper edge has to be larger than the lower edge"
                << std::endl;
      Fatal("SetFBTrackScanQAHistogramBinning",
            ": upper edge has to be larger than the lower edge");
    }
    this->fFBTrackScanQAHistogramBins[Track][kBIN] = nbins;
    this->fFBTrackScanQAHistogramBins[Track][kLEDGE] = lowerEdge;
    this->fFBTrackScanQAHistogramBins[Track][kUEDGE] = upperEdge;
  }
  // generic setter for self correlation QA histogram binning
  void SetSelfCorQAHistogramBinning(kTrack Track, Int_t nbins,
                                    Double_t lowerEdge, Double_t upperEdge) {
    if (Track >= kKinematic) {
      std::cout << __LINE__ << ": running out of bounds" << std::endl;
      Fatal("SetSelfCorQAHistogramBinning",
            "Running out of bounds in SetSelfCorQAHistogramBinning");
    }
    if (upperEdge < lowerEdge) {
      std::cout << __LINE__
                << ": upper edge has to be larger than the lower edge"
                << std::endl;
      Fatal("SetSelfCorQAHistogramBinning",
            ": upper edge has to be larger than the lower edge");
    }
    this->fSelfCorQAHistogramBins[Track][kBIN] = nbins;
    this->fSelfCorQAHistogramBins[Track][kLEDGE] = lowerEdge;
    this->fSelfCorQAHistogramBins[Track][kUEDGE] = upperEdge;
  }
  // generic setter for track control histogram binning
  void SetTrackControlHistogramBinning(kTrack Track, Int_t nbins,
                                       Double_t lowerEdge, Double_t upperEdge) {
    if (Track >= LAST_ETRACK) {
      std::cout << __LINE__ << ": running out of bounds" << std::endl;
      Fatal("SetTrackControlHistogramBinning",
            "Running out of bounds in SetTrackControlHistogramBinning");
    }
    if (upperEdge < lowerEdge) {
      std::cout << __LINE__
                << ": upper edge has to be larger than the lower edge"
                << std::endl;
      Fatal("SetEventControlHistogramBinning",
            ": upper edge has to be larger than the lower edge");
    }
    this->fTrackControlHistogramBins[Track][kBIN] = nbins;
    this->fTrackControlHistogramBins[Track][kLEDGE] = lowerEdge;
    this->fTrackControlHistogramBins[Track][kUEDGE] = upperEdge;
  }
  // generic setter for event histogram binning
  void SetEventControlHistogramBinning(kEvent Event, Int_t nbins,
                                       Double_t lowerEdge, Double_t upperEdge) {
    if (Event >= LAST_EEVENT) {
      std::cout << __LINE__ << ": running out of bounds" << std::endl;
      Fatal("SetEventControlHistogramBinning",
            "Running out of bounds in SetEventControlHistogramBinning");
    }
    if (upperEdge < lowerEdge) {
      std::cout << __LINE__
                << ": upper edge has to be larger than the lower edge"
                << std::endl;
      Fatal("SetEventControlHistogramBinning",
            ": upper edge has to be larger than the lower edge");
    }
    this->fEventControlHistogramBins[Event][kBIN] = nbins;
    this->fEventControlHistogramBins[Event][kLEDGE] = lowerEdge;
    this->fEventControlHistogramBins[Event][kUEDGE] = upperEdge;
  }

  // setters for cuts
  // centrality selection criterion
  void SetCentralityEstimator(kCenEstimators CentralityEstimator) {
    if (CentralityEstimator >= LAST_ECENESTIMATORS) {
      std::cout << __LINE__ << ": running out of bounds" << std::endl;
      Fatal("SetCentralityEstimator",
            "Running out of bounds in SetCentralityEstimator");
    }
    this->fCentralityEstimator = CentralityEstimator;
  }
  // generic setter for track cuts
  void SetTrackCuts(kTrack Track, Double_t min, Double_t max) {
    if (Track >= LAST_ETRACK) {
      std::cout << __LINE__ << ": running out of bounds" << std::endl;
      Fatal("SetTrackCuts", "Running out of bounds in SetTrackCuts");
    }
    if (max < min) {
      std::cout << __LINE__ << ": maximum has to be larger than the minimum"
                << std::endl;
      Fatal("SetTrackCuts", ": maximum has to be larger than the minimum");
    }
    this->fTrackCuts[Track][kMIN] = min;
    this->fTrackCuts[Track][kMAX] = max;
  }
  // generic setter for event cuts
  void SetEventCuts(kEvent Event, Double_t min, Double_t max) {
    if (Event >= LAST_EEVENT) {
      std::cout << __LINE__ << ": running out of bounds" << std::endl;
      Fatal("SetEventCuts", "Running out of bounds in SetEventCuts");
    }
    if (max < min) {
      std::cout << __LINE__ << ": maximum has to be larger than the minimum"
                << std::endl;
      Fatal("SetEventCuts", ": maximum has to be larger than the minimum");
    }
    this->fEventCuts[Event][kMIN] = min;
    this->fEventCuts[Event][kMAX] = max;
  }
  // setter for centrality correlation cut
  // void SetCenCorCut(Double_t m, Double_t t) {
  //   if (m < 1.) {
  //     std::cout << __LINE__ << ": slope too small" << std::endl;
  //     Fatal("SetCenCorCut", "slope too small");
  //   }
  //   if (t < 1.) {
  //     std::cout << __LINE__ << ": offset too small" << std::endl;
  //     Fatal("SetCenCorCut", "offset too small");
  //   }
  //   this->fCenCorCut[0] = m;
  //   this->fCenCorCut[1] = t;
  void SetCenCorCut(Double_t m, Double_t t) {
    this->fCenCorCut[0] = m;
    this->fCenCorCut[1] = t;
  }
  // setter for multiplicity correlation cut
  void SetMulCorCut(Double_t m, Double_t t) {
    this->fMulCorCut[0] = m;
    this->fMulCorCut[1] = t;
  }
  // filterbit
  // depends strongly on the data set
  // typical choices are 1,128,256,768
  void SetFilterbit(Int_t Filterbit) { this->fFilterbit = Filterbit; }
  // cut all non-primary particles away
  void SetPrimaryOnlyCut(Bool_t option) { this->fPrimaryOnly = option; }

  // setters for MC on the fly analysis
  void SetMCOnTheFly(Bool_t option) { this->fMCOnTheFly = option; }
  void SetMCClosure(Bool_t option) { this->fMCClosure = option; }
  void SetUseWeights(kTrack kinematic, Bool_t option) {
    if (kinematic >= kKinematic) {
      std::cout << __LINE__ << ": Out of range" << std::endl;
      Fatal("SetUseWeights", "Out of range");
    }
    this->fUseWeights[kinematic] = option;
  }
  // reset weights and redo the analysis
  void SetResetWeights(kTrack kinematic, Bool_t option) {
    if (kinematic >= kKinematic) {
      std::cout << __LINE__ << ": Out of range" << std::endl;
      Fatal("SetResetWeights", "Out of range");
    }
    this->fResetWeights[kinematic] = option;
  }
  void SetUseCustomSeed(const UInt_t seed) {
    this->fSeed = seed;
    this->fUseCustomSeed = kTRUE;
  }
  // set flow harmonics for pdf
  void SetMCFlowHarmonics(std::vector<Double_t> FlowHarmonics) {
    if (FlowHarmonics.size() >= kMaxHarmonic) {
      std::cout << __LINE__ << ": Vector exceeds maximum allowed harmonic"
                << std::endl;
      Fatal("SetFlowHarmonics", "Too many harmonics");
    }
    fMCFlowHarmonics = FlowHarmonics;
  }
  void SetMCPdfRange(Double_t min, Double_t max) {
    fMCPdfRange[kMIN] = min;
    fMCPdfRange[kMAX] = max;
  }
  void SetMCNumberOfParticlesPerEvent(Int_t n) {
    fMCNumberOfParticlesPerEvent = n;
  }
  void SetMCNumberOfParticlesPerEventRange(Int_t min, Int_t max) {
    fMCNumberOfParticlesPerEventFluctuations = kTRUE;
    fMCNumberOfParticlesPerEventRange[kMIN] = min;
    fMCNumberOfParticlesPerEventRange[kMAX] = max;
  }

  // setters for acceptance and weight histograms for monte carlo closure
  void SetAcceptanceHistogram(kTrack kinematic, TH1D *AcceptanceHistogram) {
    if (!AcceptanceHistogram) {
      std::cout << __LINE__ << ": Did not get acceptance histogram"
                << std::endl;
      Fatal("SetAccpetanceHistogram", "Invalid pointer");
    }
    if (kinematic >= kKinematic) {
      std::cout << __LINE__ << ": Out of range" << std::endl;
      Fatal("SetAccpetanceHistogram", "Out of range");
    }
    this->fAcceptanceHistogram[kinematic] = AcceptanceHistogram;
  }
  void SetAcceptanceHistogram(kTrack kinematic, const char *Filename,
                              const char *Histname);
  void SetWeightHistogram(kTrack kinematic, TH1D *WeightHistogram) {
    if (!WeightHistogram) {
      std::cout << __LINE__ << ": Did not get weight histogram" << std::endl;
      Fatal("SetWeightHistogram", "Invalid pointer");
    }
    if (kinematic >= kKinematic) {
      std::cout << __LINE__ << ": Out of range" << std::endl;
      Fatal("SetWeightHistogram", "Out of range");
    }
    this->fUseWeights[kinematic] = kTRUE;
    this->fWeightHistogram[kinematic] = WeightHistogram;
  }
  void SetWeightHistogram(kTrack kinematic, const char *Filename,
                          const char *Histname);

  // set correlators to be computed
  void SetCorrelators(std::vector<std::vector<Int_t>> correlators) {
    this->fCorrelators = correlators;
    for (int i = 0; i < LAST_EFINALPROFILE; ++i) {
      fFinalResultProfileBins[i][kBIN] = fCorrelators.size();
      fFinalResultProfileBins[i][kLEDGE] = 0;
      fFinalResultProfileBins[i][kUEDGE] = fCorrelators.size();
    }
  }

private:
  AliAnalysisTaskAR(const AliAnalysisTaskAR &aatmpf);
  AliAnalysisTaskAR &operator=(const AliAnalysisTaskAR &aatmpf);

  // base list
  TList *fHistList;
  TString fHistListName;

  // QA histograms
  TList *fQAHistogramsList;
  TString fQAHistogramsListName;
  Bool_t fFillQAHistograms;
  // array holding all centrality estimates
  Double_t fCentrality[LAST_ECENESTIMATORS];
  // centrality correlation histograms
  TList *fCenCorQAHistogramsList;
  TString fCenCorQAHistogramsListName;
  TH2D *fCenCorQAHistograms[LAST_ECENESTIMATORS * (LAST_ECENESTIMATORS - 1) / 2]
                           [LAST_EBEFOREAFTER];
  TString fCenCorQAHistogramNames[LAST_ECENESTIMATORS *
                                  (LAST_ECENESTIMATORS - 1) /
                                  2][LAST_EBEFOREAFTER][LAST_ENAME];
  Double_t fCenCorQAHistogramBins[LAST_ECENESTIMATORS *
                                  (LAST_ECENESTIMATORS - 1) /
                                  2][2 * LAST_EBINS];
  // array holding all multiplicity estimates
  Double_t fMultiplicity[kMulEstimators];
  // multiplicity correlation histograms
  TList *fMulCorQAHistogramsList;
  TString fMulCorQAHistogramsListName;
  TH2D *fMulCorQAHistograms[kMulEstimators * (kMulEstimators - 1) / 2]
                           [LAST_EBEFOREAFTER];
  TString fMulCorQAHistogramNames[kMulEstimators * (kMulEstimators - 1) / 2]
                                 [LAST_EBEFOREAFTER][LAST_ENAME];
  Double_t fMulCorQAHistogramBins[kMulEstimators * (kMulEstimators - 1) / 2]
                                 [2 * LAST_EBINS];
  // filterbit scans histograms
  TList *fFBScanQAHistogramsList;
  TString fFBScanQAHistogramsListName;
  TH1D *fFBScanQAHistogram;
  TString fFBScanQAHistogramName[LAST_ENAME];
  Double_t fFBScanQAHistogramBin[LAST_EBINS];
  TH1D *fFBTrackScanQAHistograms[LAST_ETRACK][kNumberofTestFilterBit];
  TString fFBTrackScanQAHistogramNames[LAST_ETRACK][kNumberofTestFilterBit]
                                      [LAST_ENAME];
  Double_t fFBTrackScanQAHistogramBins[LAST_ETRACK][LAST_EBINS];
  // self correlation histograms
  TList *fSelfCorQAHistogramsList;
  TString fSelfCorQAHistogramsListName;
  TH1D *fSelfCorQAHistograms[kKinematic][LAST_EBEFOREAFTER];
  TString fSelfCorQAHistogramNames[kKinematic][LAST_EBEFOREAFTER][LAST_ENAME];
  Double_t fSelfCorQAHistogramBins[kKinematic][LAST_EBINS];

  // control histograms
  TList *fControlHistogramsList;
  TString fControlHistogramsListName;
  // track control histograms
  TList *fTrackControlHistogramsList;
  TString fTrackControlHistogramsListName;
  TH1D *fTrackControlHistograms[LAST_EMODE][LAST_ETRACK][LAST_EBEFOREAFTER];
  TString fTrackControlHistogramNames[LAST_EMODE][LAST_ETRACK]
                                     [LAST_EBEFOREAFTER][LAST_ENAME];
  Double_t fTrackControlHistogramBins[LAST_ETRACK][LAST_EBINS];
  // event control historams
  TList *fEventControlHistogramsList;
  TString fEventControlHistogramsListName;
  TH1D *fEventControlHistograms[LAST_EMODE][LAST_EEVENT][LAST_EBEFOREAFTER];
  TString fEventControlHistogramNames[LAST_EMODE][LAST_EEVENT]
                                     [LAST_EBEFOREAFTER][LAST_ENAME];
  Double_t fEventControlHistogramBins[LAST_EEVENT][LAST_EBINS];

  // cuts
  Double_t fTrackCuts[LAST_ETRACK][LAST_EMINMAX];
  TH1D *fTrackCutsCounter[LAST_EMODE];
  TString fTrackCutsCounterNames[LAST_EMODE];
  TString fTrackCutsCounterBinNames[LAST_ETRACK][LAST_EMINMAX];
  Double_t fEventCuts[LAST_EEVENT][LAST_EMINMAX];
  TH1D *fEventCutsCounter[LAST_EMODE];
  TString fEventCutsCounterNames[LAST_EMODE];
  TString fEventCutsCounterBinNames[LAST_EEVENT][LAST_EMINMAX];
  Int_t fFilterbit;
  Bool_t fPrimaryOnly;
  kCenEstimators fCentralityEstimator;
  Double_t fCenCorCut[2];
  Double_t fMulCorCut[2];

  // Final results
  TList *fFinalResultsList;
  TString fFinalResultsListName;
  // array holding final result histograms
  TH1D *fFinalResultHistograms[LAST_EFINALHIST];
  TString fFinalResultHistogramNames[LAST_EFINALHIST][LAST_ENAME];
  Double_t fFinalResultHistogramBins[LAST_EFINALHIST][LAST_EBINS];
  // array holding final result profiles
  TProfile *fFinalResultProfiles[LAST_EFINALPROFILE];
  TString fFinalResultProfileNames[LAST_EFINALPROFILE][LAST_ENAME];
  Double_t fFinalResultProfileBins[LAST_EFINALPROFILE][LAST_EBINS];

  // Monte Carlo on the fly/closure
  Bool_t fMCOnTheFly;
  Bool_t fMCClosure;
  UInt_t fSeed;
  Bool_t fUseCustomSeed;
  TF1 *fMCPdf;
  TString fMCPdfName;
  Double_t fMCPdfRange[LAST_EMINMAX];
  std::vector<Double_t> fMCFlowHarmonics;
  Bool_t fMCNumberOfParticlesPerEventFluctuations;
  Int_t fMCNumberOfParticlesPerEvent;
  Int_t fMCNumberOfParticlesPerEventRange[LAST_EMINMAX];

  // Look up tabel between MC and data particles
  TExMap *fLookUpTable;

  // qvectors
  TComplex fQvector[kMaxHarmonic][kMaxPower];
  std::vector<Double_t> fKinematics[kKinematic];
  std::vector<Double_t> fKinematicWeights[kKinematic];
  std::vector<Double_t> fWeightsAggregated;
  TH1D *fAcceptanceHistogram[kKinematic];
  TH1D *fWeightHistogram[kKinematic];
  Bool_t fUseWeights[kKinematic];
  Bool_t fUseWeightsAggregated;
  Bool_t fResetWeights[kKinematic];
  Bool_t fResetWeightsAggregated;
  std::vector<std::vector<Int_t>> fCorrelators;

  // increase this counter in each new version
  ClassDef(AliAnalysisTaskAR, 10);
};

#endif
