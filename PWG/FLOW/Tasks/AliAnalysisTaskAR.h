/**
 * File              : AliAnalysisTaskAR.h
 * Author            : Anton Riedel <anton.riedel@tum.de>
 * Date              : 07.05.2021
 * Last Modified Date: 18.10.2021
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
#include "AliAODMCParticle.h"
#include "AliAODTrack.h"
#include "AliAnalysisTaskSE.h"
#include "AliVEvent.h"
#include "AliVParticle.h"
#include <Riostream.h>
#include <TComplex.h>
#include <TDataType.h>
#include <TExMap.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
// #include <THnSparse.h>
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
const Int_t kNumberofTestFilterBit = 6;
const Int_t kTestFilterbit[kNumberofTestFilterBit] = {1,   92,  128,
                                                      256, 512, 768};
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
  kTPCCROSSEDROWS,
  kTPCNCLSFRACTIONSHARED,
  kTPCCHI2PERNDF,
  kITSNCLS,
  kCHI2PERNDF,
  kDCAZ,
  kDCAXY,
  LAST_ETRACK
};
// kinematic variables
const Int_t kKinematic = kETA + 1;
// final result histograms
enum kFinalResultHist {
  kAVGPHI,
  kAVGCEN,
  kMINMUL,
  kNUMBEROFEVENTS,
  kNUMBEROFTRACKS,
  LAST_EFINALHIST
};
enum kFinalResultProfile {
  kINTEGRATED,
  kCENDEP,
  kMULDEP,
  LAST_EFINALRESULTPROFILE
};
// various gloabl objects
enum kBins { kBIN, kLEDGE, kUEDGE, LAST_EBINS };
enum kName { kNAME, kTITLE, kXAXIS, kYAXIS, LAST_ENAME };
enum kMinMax { kMIN, kMAX, LAST_EMINMAX };
const TString kMMName[LAST_EMINMAX] = {"[kMIN]", "[kMAX]"};
enum kBeforeAfter { kBEFORE, kAFTER, LAST_EBEFOREAFTER };
const TString kBAName[LAST_EBEFOREAFTER] = {"[kBEFORE]", "[kAFTER]"};
const Color_t kFillColor[LAST_EBEFOREAFTER] = {kRed - 10, kGreen - 10};
const Color_t kcolorFinalResult = kBlue - 10;
enum kMode { kRECO, kSIM, LAST_EMODE };
const TString kModeName[LAST_EMODE] = {"[kRECO]", "[kSIM]"};
enum kMCPrimaryDef { kMCPrim, kMCPhysicalPrim };

class AliAnalysisTaskAR : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskAR();
  AliAnalysisTaskAR(const char *name);
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
  virtual void InitializeArraysForMCAnalysis();

  // methods called in UserCreateOutputObjects()
  virtual void BookAndNestAllLists();
  virtual void BookQAHistograms();
  virtual void BookControlHistograms();
  virtual void BookFinalResultHistograms();
  virtual void BookFinalResultCorrelators();
  virtual void BookFinalResultSymmetricCumulants();

  // functions for setting default values for binning/cuts
  void SetDefaultConfiguration();
  void SetDefaultBinning();
  void SetDefaultCuts(Int_t Filterbit = 128, Double_t cenMin = 0,
                      Double_t cenMax = 100);

  // functions called in UserExec()
  virtual void FillEventQAHistograms(kBeforeAfter BA, AliAODEvent *AODEvent,
                                     AliMCEvent *MCEvent);
  virtual void FillFBScanQAHistograms(AliAODTrack *track);
  virtual void FillEventControlHistograms(kBeforeAfter BA, AliVEvent *Event);
  virtual void FillTrackControlHistograms(kBeforeAfter BA, AliVParticle *track);
  virtual void FillFinalResultCorrelators();
  virtual void FillSymmetricCumulant();
  virtual Bool_t SurviveEventCut(AliAODEvent *aAOD);
  virtual Bool_t SurviveTrackCut(AliVParticle *aTrack, Bool_t FillCounter);
  virtual void FillEventObjects(AliAODEvent *aAOD, AliMCEvent *aMC);
  virtual void ClearVectors();
  virtual void FillTrackObjects(AliVParticle *avp);
  virtual void AggregateWeights();
  virtual Int_t IndexCorHistograms(Int_t i, Int_t j, Int_t N);

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
  TComplex FiveNestedLoops(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5);
  TComplex SixNestedLoops(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5,
                          Int_t n6);

  // methods for calculating symmetric cumulants
  void SC2(kFinalResultProfile rp);
  void SC3(kFinalResultProfile rp);

  // GetPointers Methods in case we need to manually trigger Terminate()
  virtual void GetPointers(TList *list);
  virtual void GetPointersForQAHistograms();
  virtual void GetPointersForControlHistograms();
  virtual void GetPointersForFinalResults();

  // setters and getters for list objects
  TList *GetCenCorQAHistogramList() const { return fCenCorQAHistogramsList; }
  TObject *GetCenCorQAHistogram(kBeforeAfter ba, kCenEstimators cen1,
                                kCenEstimators cen2) {
    return fCenCorQAHistogramsList->FindObject(
        fCenCorQAHistogramNames[IndexCorHistograms(
            cen1, cen2, LAST_ECENESTIMATORS)][ba][0]);
  }
  TList *GetMulCorQAHistogramList() const { return fMulCorQAHistogramsList; }
  TObject *GetMulCorQAHistogram(kBeforeAfter ba, Int_t mul1, Int_t mul2) {
    return fMulCorQAHistogramsList->FindObject(
        fMulCorQAHistogramNames[IndexCorHistograms(mul1, mul2, kMulEstimators)]
                               [ba][0]);
  }
  TList *GetSelfCorQAHistogramList() const { return fSelfCorQAHistogramsList; }
  void SetControlHistogramsList(TList *const chl) {
    this->fControlHistogramsList = chl;
  };
  TList *GetControlHistogramsList() const {
    return this->fControlHistogramsList;
  }
  TObject *GetEventControlHistogram(kMode mode, kEvent event,
                                    kBeforeAfter ba) const {
    return fEventControlHistogramsList->FindObject(
        fEventControlHistogramNames[mode][event][ba][0]);
  }
  TObject *GetTrackControlHistogram(kMode mode, kTrack track,
                                    kBeforeAfter ba) const {
    return fTrackControlHistogramsList->FindObject(
        fTrackControlHistogramNames[mode][track][ba][0]);
  }
  void SetFinalResultsList(TList *const frl) { this->fFinalResultsList = frl; };
  TList *GetFinalResultsList() const { return this->fFinalResultsList; }

  void SetFinalResultHistogramsList(TList *const frl) {
    this->fFinalResultHistogramsList = frl;
  };
  TList *GetFinalResultHistogramssList() const {
    return this->fFinalResultHistogramsList;
  }
  void SetFinalResultProfilesList(TList *const frl) {
    this->fFinalResultCorrelatorsList = frl;
  };
  TList *GetFinalResultProfilesList() const {
    return this->fFinalResultCorrelatorsList;
  }
  void GetCorrelatorValues(std::vector<Double_t> *cor) {
    cor->clear();
    TList *list;
    for (auto List : *fFinalResultCorrelatorsList) {
      list = dynamic_cast<TList *>(List);
      cor->push_back(dynamic_cast<TProfile *>(list->At(0))->GetBinContent(1));
    }
  }

  // setters for QA histograms
  void SetFillQAHistograms(Bool_t option) { fFillQAHistograms = option; }
  void SetFillQACorHistogramsOnly(Bool_t option) {
    fFillQACorHistogramsOnly = option;
  }
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
  // generic setter for multiplicity correlation QA histogram binning
  void SetCenMulCorQAHistogramBinning(Int_t xnbins, Double_t xlowerEdge,
                                      Double_t xupperEdge, Int_t mul,
                                      Int_t ynbins, Double_t ylowerEdge,
                                      Double_t yupperEdge);
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
  // only fill control histograms
  void SetFillControlHistogramsOnly(Bool_t option) {
    this->fFillControlHistogramsOnly = option;
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
    this->fUseTrackCuts[Track] = kTRUE;
  }
  void SetTrackCuts(kTrack Track, Bool_t option) {
    this->fUseTrackCuts[Track] = option;
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
    this->fUseEventCuts[Event] = kTRUE;
  }
  void SetEventCuts(kEvent Event, Bool_t option) {
    this->fUseEventCuts[Event] = option;
  }
  void SetCenCorCut(Double_t m, Double_t t) {
    this->fCenCorCut[0] = m;
    this->fCenCorCut[1] = t;
    this->fUseCenCorCuts = kTRUE;
  }
  void SetCenCorCut(Bool_t option) { this->fUseCenCorCuts = option; }
  // setter for multiplicity correlation cut
  void SetMulCorCut(Double_t m, Double_t t) {
    this->fMulCorCut[0] = m;
    this->fMulCorCut[1] = t;
    this->fUseMulCorCuts = kTRUE;
  }
  void SetMulCorCut(Bool_t option) { this->fUseMulCorCuts = option; }
  // filterbit
  // depends strongly on the data set
  // typical choices are 1,128,256,768
  void SetFilterbit(Int_t Filterbit) {
    this->fFilterbit = Filterbit;
    this->fUseFilterbit = kTRUE;
  }
  // cut all neutral particles away
  void SetChargedOnlyCut(Bool_t option) { this->fChargedOnly = option; }
  // cut all non-primary particles away
  void SetPrimaryOnlyCut(Bool_t option) { this->fPrimaryOnly = option; }
  void SetPrimaryDefinitionInMC(kMCPrimaryDef def) {
    this->fMCPrimaryDef = def;
  }
  // cut all non-global track away
  void SetGlobalTracksOnlyCut(Bool_t option) {
    this->fGlobalTracksOnly = option;
  }
  // flatten centrality if necessary
  void SetCenFlattenHist(TH1D *hist) {
    if (!hist) {
      std::cout << __LINE__ << ": Did not get centrality flattening histogram"
                << std::endl;
      Fatal("SetCenFlattenHist", "Invalid pointer");
    }
    this->fUseCenFlatten = kTRUE;
    this->fCenFlattenHist = hist;
  };
  void SetCenFlattenHist(const char *Filename, const char *Histname);

  // setters for MC on the fly/closure analysis
  void SetMCClosure(Bool_t option) { this->fMCClosure = option; }
  void SetMCOnTheFly(Bool_t option) { this->fMCOnTheFly = option; }
  void SetCustomSeed(const UInt_t seed) {
    this->fSeed = seed;
    this->fUseCustomSeed = kTRUE;
  }
  void SetMCMultiplicityPdf(TF1 *pdf) { this->fMCMultiplicity = pdf; }
  void SetMCKinematicPdf(Int_t kinematic, TF1 *pdf) {
    if (kinematic >= kKinematic) {
      std::cout << __LINE__ << ": Out of range" << std::endl;
      Fatal("SetMCKinematicPdf", "Out of range");
    }
    this->fMCKinematicPDFs[kinematic] = pdf;
  }
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

  // setters for weight histograms
  void SetUseWeights(kTrack kinematic, Bool_t option) {
    if (kinematic >= kKinematic) {
      std::cout << __LINE__ << ": Out of range" << std::endl;
      Fatal("SetUseWeights", "Out of range");
    }
    this->fUseWeights[kinematic] = option;
  }
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
  }

  // set symmetric cumulant to be computed
  void SetSymmetricCumulant(std::vector<Int_t> SC) {
    this->fSymmetricCumulants = SC;
  }

  // use nested loops for computation of correlators
  void SetUseNestedLoops(Bool_t option) { this->fUseNestedLoops = option; }

  // use fischer-yates for indices randomization
  void SetFisherYates(Bool_t option) { this->fUseFisherYates = option; }
  // use a fixed multiplicity
  void SetFixedMultiplicity(Int_t FixedMultiplicity) {
    this->fUseFixedMultplicity = kTRUE;
    this->fFixedMultiplicy = FixedMultiplicity;
  }
  void SetUseFixedMultiplicity(Bool_t option) {
    this->fUseFixedMultplicity = option;
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
  // only fill correlation histograms
  Bool_t fFillQACorHistogramsOnly;
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
  // multiplicity-centrality correlation histograms
  TList *fCenMulCorQAHistogramsList;
  TString fCenMulCorQAHistogramsListName;
  TH2D *fCenMulCorQAHistograms[kMulEstimators][LAST_EBEFOREAFTER];
  TString fCenMulCorQAHistogramNames[kMulEstimators][LAST_EBEFOREAFTER]
                                    [LAST_ENAME];
  Double_t fCenMulCorQAHistogramBins[kMulEstimators][2 * LAST_EBINS];
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
  Bool_t fUseTrackCuts[LAST_ETRACK];
  TH1D *fTrackCutsCounter[LAST_EMODE];
  TString fTrackCutsCounterCumulativeName;
  // THnSparseD *fTrackCutsCounterCumulative;
  TString fTrackCutsValuesName;
  TH1D *fTrackCutsValues;
  TString fTrackCutsCounterNames[LAST_EMODE];
  TString fTrackCutsCounterBinNames[LAST_ETRACK][LAST_EMINMAX];
  Double_t fEventCuts[LAST_EEVENT][LAST_EMINMAX];
  Bool_t fUseEventCuts[LAST_EEVENT];
  TH1D *fEventCutsCounter[LAST_EMODE];
  TString fEventCutsCounterCumulativeName;
  // THnSparseD *fEventCutsCounterCumulative;
  TString fEventCutsValuesName;
  TH1D *fEventCutsValues;
  TString fEventCutsCounterNames[LAST_EMODE];
  TString fEventCutsCounterBinNames[LAST_EEVENT][LAST_EMINMAX];
  Int_t fFilterbit;
  Bool_t fUseFilterbit;
  Bool_t fChargedOnly;
  Bool_t fPrimaryOnly;
  kMCPrimaryDef fMCPrimaryDef;
  Bool_t fGlobalTracksOnly;
  kCenEstimators fCentralityEstimator;
  Double_t fCenCorCut[2];
  Bool_t fUseCenCorCuts;
  Double_t fMulCorCut[2];
  Bool_t fUseMulCorCuts;
  Bool_t fUseCenFlatten;
  TH1D *fCenFlattenHist;

  // Final results
  TList *fFinalResultsList;
  TString fFinalResultsListName;
  // array holding final result histograms
  TList *fFinalResultHistogramsList;
  TString fFinalResultHistogramsListName;
  TH1D *fFinalResultHistograms[LAST_EFINALHIST];
  TString fFinalResultHistogramNames[LAST_EFINALHIST][LAST_ENAME];
  Double_t fFinalResultHistogramBins[LAST_EFINALHIST][LAST_EBINS];
  // final result correlators
  // will be generated by CreateUserObjects depending on the correlators we have
  TList *fFinalResultCorrelatorsList;
  TString fFinalResultCorrelatorsListName;
  TList *fFinalResultSymmetricCumulantsList;
  TString fFinalResultSymmetricCumulantsListName;
  // only fill control histograms
  Bool_t fFillControlHistogramsOnly;
  Bool_t fUseNestedLoops;

  // Seed for RNG
  Bool_t fUseCustomSeed;
  UInt_t fSeed;

  // Monte Carlo on the fly/Closure
  Bool_t fMCOnTheFly;
  Bool_t fMCClosure;
  TF1 *fMCKinematicPDFs[kKinematic];
  Double_t fMCKinematicVariables[kKinematic];
  TH1D *fAcceptanceHistogram[kKinematic];
  TF1 *fMCMultiplicity;

  // Look up tabel between MC and data particles
  TExMap *fLookUpTable;

  // use Fisher-Yates algorithm to randomize tracks
  Bool_t fUseFisherYates;
  std::vector<Int_t> fRandomizedTrackIndices;
  // needed when sampling a fixed number of tracks per event
  Bool_t fUseFixedMultplicity;
  Int_t fFixedMultiplicy;

  // qvectors
  TComplex fQvector[kMaxHarmonic][kMaxPower];
  std::vector<Double_t> fKinematics[kKinematic];
  std::vector<Double_t> fKinematicWeights[kKinematic];
  std::vector<Double_t> fWeightsAggregated;
  TH1D *fWeightHistogram[kKinematic];
  Bool_t fUseWeights[kKinematic];
  Bool_t fUseWeightsAggregated;
  std::vector<std::vector<Int_t>> fCorrelators;
  std::vector<Int_t> fSymmetricCumulants;

  // increase this counter in each new version
  ClassDef(AliAnalysisTaskAR, 16);
};

#endif
