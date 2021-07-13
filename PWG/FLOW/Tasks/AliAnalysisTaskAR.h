/**
 * File              : AliAnalysisTaskAR.h
 * Author            : Anton Riedel <anton.riedel@tum.de>
 * Date              : 07.05.2021
 * Last Modified Date: 13.07.2021
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
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TString.h>

// global constants
const Int_t kMaxHarmonic = 20;
const Int_t kMaxCorrelator = 20;
const Int_t kMaxPower = 20;

// enumerations
enum kEvent { kCEN, kMUL, kNCONTRIB, LAST_EEVENT };
enum kTrack { kPT, kPHI, kETA, kDCAZ, kDCAXY, LAST_ETRACK };
enum kXYZ { kX, kY, kZ, LAST_EXYZ };
enum kFinalHist { kPHIAVG, LAST_EFINALHIST };
enum kFinalProfile { kHARDATA, kHARDATARESET, kHARTHEO, LAST_EFINALPROFILE };
enum kBins { kBIN, kLEDGE, kUEDGE, LAST_EBINS };
enum kName { kNAME, kTITLE, kXAXIS, LAST_ENAME };
enum kBeforeAfter { kBEFORE, kAFTER, LAST_EBEFOREAFTER };
enum kMinMax { kMIN, kMAX, LAST_EMINMAX };

class AliAnalysisTaskAR : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskAR();
  AliAnalysisTaskAR(const char *name, Bool_t useParticleWeights = kFALSE);
  virtual ~AliAnalysisTaskAR();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);

  // Methods called in the constructor:
  virtual void InitializeArrays();
  virtual void InitializeArraysForTrackControlHistograms();
  virtual void InitializeArraysForEventControlHistograms();
  virtual void InitializeArraysForCuts();
  virtual void InitializeArraysForQvectors();
  virtual void InitializeArraysForFinalResultHistograms();
  virtual void InitializeArraysForFinalResultProfiles();
  virtual void InitializeArraysForMCAnalysis();

  // Methods called in UserCreateOutputObjects():
  virtual void BookAndNestAllLists();
  virtual void BookControlHistograms();
  virtual void BookFinalResultHistograms();
  virtual void BookFinalResultProfiles();
  virtual void BookMCObjects();

  // split calls in UserExec() depending on flags
  virtual void MCOnTheFlyExec();
  virtual void FillFinalResultProfile(kFinalProfile fp);

  // methods called in AODExec():
  virtual Bool_t SurviveEventCut(AliVEvent *ave);
  virtual Bool_t SurviveTrackCut(AliAODTrack *aTrack);

  // methods called MCOnTheFlyExec()
  virtual void MCPdfSymmetryPlanesSetup();
  virtual Int_t GetMCNumberOfParticlesPerEvent();

  // Methods for computing qvectors
  void CalculateQvectors();
  TComplex Q(Int_t n, Int_t p);
  TComplex Two(Int_t n1, Int_t n2);
  TComplex Three(Int_t n1, Int_t n2, Int_t n3);
  TComplex Four(Int_t n1, Int_t n2, Int_t n3, Int_t n4);
  TComplex Five(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5);
  TComplex Six(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6);
  TComplex Recursion(Int_t n, Int_t *harmonic, Int_t mult = 1, Int_t skip = 0);

  // methods for computing nested loops
  TComplex TwoNestedLoops(Int_t n1, Int_t n2);
  TComplex ThreeNestedLoops(Int_t n1, Int_t n2, Int_t n3);
  TComplex FourNestedLoops(Int_t n1, Int_t n2, Int_t n3, Int_t n4);
  // TComplex FiveNestedLoops(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5);
  // TComplex SixNestedLoops(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5,
  // Int_t n6);
  Double_t CombinatorialWeight(Int_t n);

  // GetPointers Methods in case we need to manually trigger Terminate()
  virtual void GetPointers(TList *list);
  virtual void GetPointersForControlHistograms();
  virtual void GetPointersForFinalResultHistograms();
  virtual void GetPointersForFinalResultProfiles();

  // Setters and getters
  void SetControlHistogramsList(TList *const chl) {
    this->fControlHistogramsList = chl;
  };
  TList *GetControlHistogramsList() const {
    return this->fControlHistogramsList;
  }
  void SetFinalResultsList(TList *const frl) { this->fFinalResultsList = frl; };
  TList *GetFinalResultsList() const { return this->fFinalResultsList; }

  // generic setter for track control histogram binning
  void SetTrackControlHistogramBinning(kTrack Track, Int_t nbins, Double_t min,
                                       Double_t max) {
    if (Track > LAST_ETRACK) {
      std::cout << __LINE__ << ": running out of bounds" << std::endl;
      Fatal("SetTrackControlHistogramBinning",
            "Running out of bounds in SetTrackControlHistogramBinning");
    }
    this->fBinsTrackControlHistograms[Track][kBIN] = nbins;
    this->fBinsTrackControlHistograms[Track][kLEDGE] = min;
    this->fBinsTrackControlHistograms[Track][kUEDGE] = max;
  }
  // generic setter for event histogram binning
  void SetEventControlHistogramBinning(kEvent Event, Int_t nbins, Double_t min,
                                       Double_t max) {
    if (Event > LAST_EEVENT) {
      std::cout << __LINE__ << ": running out of bounds" << std::endl;
      Fatal("SetEventControlHistogramBinning",
            "Running out of bounds in SetEventControlHistogramBinning");
    }
    this->fBinsEventControlHistograms[Event][kBIN] = nbins;
    this->fBinsEventControlHistograms[Event][kLEDGE] = min;
    this->fBinsEventControlHistograms[Event][kUEDGE] = max;
  }

  // setters for cuts
  // centrality selection criterion
  // only use V0M, CL0/1, SPDTracklets
  void SetCentralitySelCriterion(TString SelCriterion) {
    this->fCentralitySelCriterion = SelCriterion;
  }
  // generic setter for track cuts
  void SetTrackCuts(kTrack Track, Double_t min, Double_t max) {
    if (Track > LAST_ETRACK) {
      std::cout << __LINE__ << ": running out of bounds" << std::endl;
      Fatal("SetTrackCuts", "Running out of bounds in SetTrackCuts");
    }
    this->fTrackCuts[Track][kMIN] = min;
    this->fTrackCuts[Track][kMAX] = max;
  }
  // generic setter for event cuts
  void SetEventCuts(kEvent Event, Double_t min, Double_t max) {
    if (Event > LAST_EEVENT) {
      std::cout << __LINE__ << ": running out of bounds" << std::endl;
      Fatal("SetEventCuts", "Running out of bounds in SetEventCuts");
    }
    this->fEventCuts[Event][kMIN] = min;
    this->fEventCuts[Event][kMAX] = max;
  }
  // generic setter primary vertex cut
  void SetPrimaryVertexCuts(kXYZ xyz, Double_t min, Double_t max) {
    if (xyz > LAST_EXYZ) {
      std::cout << __LINE__ << ": running out of bounds" << std::endl;
      Fatal("SetPrimaryVertexCuts",
            "Running out of bounds in SetPrimaryVertexCuts");
    }
    this->fPrimaryVertexCuts[xyz][kMIN] = min;
    this->fPrimaryVertexCuts[xyz][kMAX] = max;
  }
  // filterbit
  // depends strongly on the data set
  // typical choices are 1,128,256,768
  void SetFilterbit(Int_t Filterbit) { this->fFilterbit = Filterbit; }
  // cut all non-primary particles
  void SetPrimaryOnlyCut(Bool_t option) { this->fPrimaryOnly = option; }

  // setters for MC analsys
  void SetMCAnalysis(Bool_t option) { this->fMCAnalaysis = option; }
  void SetMCClosure(Bool_t option) { this->fMCClosure = option; }
  void SetUseWeights(Bool_t option) { this->fUseWeights = option; }
  void SetResetWeights(Bool_t option) { this->fResetWeights = option; }
  void SetUseCustomSeed(const UInt_t seed) {
    this->fSeed = seed;
    this->fUseCustomSeed = kTRUE;
  }
  void SetMCFlowHarmonics(TArrayD *array) {
    if (array->GetSize() > kMaxHarmonic) {
      std::cout << __LINE__ << ": Array exceeds maximum allowed harmonic"
                << std::endl;
      Fatal("SetFlowHarmonics", "Too many harmonics");
    }
    fMCFlowHarmonics = array;
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
  void SetAcceptanceHistogram(TH1F *AcceptanceHistogram) {
    if (!AcceptanceHistogram) {
      std::cout << __LINE__ << ": Did not get acceptance histogram"
                << std::endl;
      Fatal("SetAccpetanceHistogram", "Invalid pointer");
    }
    this->fAcceptanceHistogram = AcceptanceHistogram;
  }
  void SetAcceptanceHistogram(const char *Filename, const char *Histname);
  void SetWeightHistogram(TH1F *WeightHistogram) {
    if (!WeightHistogram) {
      std::cout << __LINE__ << ": Did not get weight histogram" << std::endl;
      Fatal("SetWeightHistogram", "Invalid pointer");
    }
    this->fWeightHistogram = WeightHistogram;
  }
  void SetWeightHistogram(const char *Filename, const char *Histname);

  // set correlators we want to compute
  void SetCorrelators(std::vector<std::vector<Int_t>> correlators) {
    this->fCorrelators = correlators;
    for (int i = 0; i < LAST_EFINALPROFILE; ++i) {
      fBinsFinalResultProfiles[i][kBIN] = fCorrelators.size();
      fBinsFinalResultProfiles[i][kLEDGE] = 0;
      fBinsFinalResultProfiles[i][kUEDGE] = fCorrelators.size();
    }
  }

private:
  AliAnalysisTaskAR(const AliAnalysisTaskAR &aatmpf);
  AliAnalysisTaskAR &operator=(const AliAnalysisTaskAR &aatmpf);

  // base list holding all output object (a.k.a. grandmother of all lists)
  TList *fHistList;
  TString fHistListName;

  // control histograms
  TList *fControlHistogramsList;
  TString fControlHistogramsListName;
  TH1F *fTrackControlHistograms[LAST_ETRACK][LAST_EBEFOREAFTER];
  TString fTrackControlHistogramNames[LAST_ETRACK][LAST_EBEFOREAFTER]
                                     [LAST_ENAME];
  // array holding bins of track control histograms
  Double_t fBinsTrackControlHistograms[LAST_ETRACK][LAST_EBINS];
  // array holding event control histograms
  TH1F *fEventControlHistograms[LAST_EEVENT][LAST_EBEFOREAFTER];
  TString fEventControlHistogramNames[LAST_EEVENT][LAST_EBEFOREAFTER]
                                     [LAST_ENAME];
  Double_t fBinsEventControlHistograms[LAST_EEVENT][LAST_EBINS];

  // cuts
  TString fCentralitySelCriterion;
  Double_t fTrackCuts[LAST_ETRACK][LAST_EMINMAX];
  Double_t fEventCuts[LAST_EEVENT][LAST_EMINMAX];
  Double_t fPrimaryVertexCuts[LAST_EXYZ][LAST_EMINMAX];
  Int_t fFilterbit;
  Bool_t fPrimaryOnly;

  // Final results
  TList *fFinalResultsList;
  TString fFinalResultsListName;
  // array holding final result histograms
  TH1F *fFinalResultHistograms[LAST_EFINALHIST];
  TString fFinalResultHistogramNames[LAST_EFINALHIST][LAST_ENAME];
  Double_t fBinsFinalResultHistograms[LAST_EFINALHIST][LAST_EBINS];
  // arayy holding final resutl profiles
  TProfile *fFinalResultProfiles[LAST_EFINALPROFILE];
  TString fFinalResultProfileNames[LAST_EFINALPROFILE][LAST_ENAME];
  Double_t fBinsFinalResultProfiles[LAST_EFINALPROFILE][LAST_EBINS];

  // Monte Carlo analysis
  TList *fMCAnalysisList;
  TString fMCAnalysisListName;
  Bool_t fMCAnalaysis;
  Bool_t fMCClosure;
  UInt_t fSeed;
  Bool_t fUseCustomSeed;
  TF1 *fMCPdf;
  TString fMCPdfName;
  Double_t fMCPdfRange[LAST_EMINMAX];
  TArrayD *fMCFlowHarmonics;
  Bool_t fMCNumberOfParticlesPerEventFluctuations;
  Int_t fMCNumberOfParticlesPerEvent;
  Int_t fMCNumberOfParticlesPerEventRange[LAST_EMINMAX];

  // qvectors
  TList *fQvectorList;
  TComplex fQvector[kMaxHarmonic][kMaxPower];
  std::vector<Double_t> fPhi;
  std::vector<Double_t> fWeights;
  TH1F *fAcceptanceHistogram;
  TH1F *fWeightHistogram;
  Bool_t fUseWeights;
  Bool_t fResetWeights;
  std::vector<std::vector<Int_t>> fCorrelators;

  // Increase this counter in each new version:
  ClassDef(AliAnalysisTaskAR, 6);
};

#endif
