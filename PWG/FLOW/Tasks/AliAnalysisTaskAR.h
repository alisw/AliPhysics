/**
 * File              : AliAnalysisTaskAR.h
 * Author            : Anton Riedel <anton.riedel@tum.de>
 * Date              : 07.05.2021
 * Last Modified Date: 05.07.2021
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
#include "Riostream.h"
#include "TComplex.h"
#include "TDataType.h"
#include "TF1.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TRandom3.h"
#include "TString.h"

/* global constants */
const Int_t kMaxHarmonic = 20;
const Int_t kMaxCorrelator = 20;
const Int_t kMaxPower = 20;

/* enumerations */
enum kEvent { kCEN, kMUL, LAST_EEVENT };
enum kTrack { kPT, kPHI, kETA, LAST_ETRACK };
enum kFinalHist { kPHIAVG, LAST_EFINALHIST };
enum kFinalProfile { kHARDATA, kHARDATARESET, kHARTHEO, LAST_EFINALPROFILE };
enum kBins { kBIN, kLEDGE, kUEDGE, LAST_EBINS };
enum kName { kNAME, kTITLE, kXAXIS, LAST_ENAME };
enum kBeforeAfter { kBEFORE, kAFTER, LAST_EBEFOREAFTER };
enum kMinMax { kMIN, kMAX, LAST_EMINMAX };
enum kXYZ { kX, kY, kZ, LAST_EXYZ };

class AliAnalysisTaskAR : public AliAnalysisTaskSE {
public:
  AliAnalysisTaskAR();
  AliAnalysisTaskAR(const char *name, Bool_t useParticleWeights = kFALSE);
  virtual ~AliAnalysisTaskAR();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);

  /* Methods called in the constructor: */
  virtual void InitializeArrays();
  virtual void InitializeArraysForTrackControlHistograms();
  virtual void InitializeArraysForEventControlHistograms();
  virtual void InitializeArraysForCuts();
  virtual void InitializeArraysForQvectors();
  virtual void InitializeArraysForFinalResultHistograms();
  virtual void InitializeArraysForFinalResultProfiles();
  virtual void InitializeArraysForMCAnalysis();

  /* Methods called in UserCreateOutputObjects(): */
  virtual void BookAndNestAllLists();
  virtual void BookControlHistograms();
  virtual void BookFinalResultHistograms();
  virtual void BookFinalResultProfiles();
  virtual void BookMCObjects();

  /* split calls in UserExec() depending on flags */
  virtual void AODExec();
  virtual void MCOnTheFlyExec();
  virtual void FillFinalResultProfiles(kFinalProfile fp);

  /* methods called in AODExec(): */
  virtual Bool_t SurviveEventCut(AliVEvent *ave);
  virtual Bool_t SurviveTrackCut(AliAODTrack *aTrack);

  /* methods called MCOnTheFlyExec() */
  virtual void MCPdfSymmetryPlanesSetup();
  virtual Int_t GetMCNumberOfParticlesPerEvent();

  /* Methods for computing qvectors */
  void CalculateQvectors();
  TComplex Q(Int_t n, Int_t p);
  TComplex Two(Int_t n1, Int_t n2);
  TComplex Three(Int_t n1, Int_t n2, Int_t n3);
  TComplex Four(Int_t n1, Int_t n2, Int_t n3, Int_t n4);
  TComplex Five(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5);
  TComplex Six(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5, Int_t n6);
  TComplex Recursion(Int_t n, Int_t *harmonic, Int_t mult = 1, Int_t skip = 0);

  /* methods for computing nested loops */
  TComplex TwoNestedLoops(Int_t n1, Int_t n2);
  TComplex ThreeNestedLoops(Int_t n1, Int_t n2, Int_t n3);
  TComplex FourNestedLoops(Int_t n1, Int_t n2, Int_t n3, Int_t n4);
  /* TComplex FiveNestedLoops(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t
   * n5);
   */
  /* TComplex SixNestedLoops(Int_t n1, Int_t n2, Int_t n3, Int_t n4, Int_t n5,
   */
  /*                         Int_t n6); */
  Double_t CombinatorialWeight(Int_t n);

  /* GetPointers Methods in case we need to manually trigger Terminate()
   */
  virtual void GetPointers(TList *list);
  virtual void GetPointersForControlHistograms();
  virtual void GetPointersForFinalResultHistograms();
  virtual void GetPointersForFinalResultProfiles();

  /* Setters and getters for data analysis */
  void SetControlHistogramsList(TList *const chl) {
    this->fControlHistogramsList = chl;
  };
  TList *GetControlHistogramsList() const {
    return this->fControlHistogramsList;
  }
  void SetFinalResultsList(TList *const frl) { this->fFinalResultsList = frl; };
  TList *GetFinalResultsList() const { return this->fFinalResultsList; }

  void SetPtBinning(Int_t nbins, Double_t min, Double_t max) {
    this->fBinsTrackControlHistograms[kPT][kBIN] = nbins;
    this->fBinsTrackControlHistograms[kPT][kLEDGE] = min;
    this->fBinsTrackControlHistograms[kPT][kUEDGE] = max;
  };
  void SetPhiBinning(Int_t nbins, Double_t min, Double_t max) {
    this->fBinsTrackControlHistograms[kPHI][kBIN] = nbins;
    this->fBinsTrackControlHistograms[kPHI][kLEDGE] = min;
    this->fBinsTrackControlHistograms[kPHI][kUEDGE] = max;
  };
  void SetEtaBinning(Int_t nbins, Double_t min, Double_t max) {
    this->fBinsTrackControlHistograms[kETA][kBIN] = nbins;
    this->fBinsTrackControlHistograms[kETA][kLEDGE] = min;
    this->fBinsTrackControlHistograms[kETA][kUEDGE] = max;
  };

  void SetCenBinning(Int_t nbins, Double_t min, Double_t max) {
    this->fBinsEventControlHistograms[kCEN][kBIN] = nbins;
    this->fBinsEventControlHistograms[kCEN][kLEDGE] = min;
    this->fBinsEventControlHistograms[kCEN][kUEDGE] = max;
  };

  void SetMulBinning(Int_t nbins, Double_t min, Double_t max) {
    this->fBinsEventControlHistograms[kMUL][kBIN] = nbins;
    this->fBinsEventControlHistograms[kMUL][kLEDGE] = min;
    this->fBinsEventControlHistograms[kMUL][kUEDGE] = max;
  };

  void SetCentralitySelCriterion(TString SelCriterion) {
    this->fCentralitySelCriterion = SelCriterion;
  }

  void SetPtCuts(Double_t min, Double_t max) {
    this->fTrackCuts[kPT][kMIN] = min;
    this->fTrackCuts[kPT][kMAX] = max;
  }

  void SetPhiCuts(Double_t min, Double_t max) {
    this->fTrackCuts[kPHI][kMIN] = min;
    this->fTrackCuts[kPHI][kMAX] = max;
  }

  void SetEtaCuts(Double_t min, Double_t max) {
    this->fTrackCuts[kETA][kMIN] = min;
    this->fTrackCuts[kETA][kMAX] = max;
  }

  void SetPrimaryVertexXCuts(Double_t min, Double_t max) {
    this->fPrimaryVertexCuts[kX][kMIN] = min;
    this->fPrimaryVertexCuts[kX][kMAX] = max;
  }

  void SetPrimaryVertexYCuts(Double_t min, Double_t max) {
    this->fPrimaryVertexCuts[kY][kMIN] = min;
    this->fPrimaryVertexCuts[kY][kMAX] = max;
  }

  void SetPrimaryVertexZCuts(Double_t min, Double_t max) {
    this->fPrimaryVertexCuts[kZ][kMIN] = min;
    this->fPrimaryVertexCuts[kZ][kMAX] = max;
  }

  void SetFilterbit(Int_t Filterbit) { this->fFilterbit = Filterbit; }

  /* setters and getters for MC analsys */
  void SetMCAnalysis(Bool_t option) { this->fMCAnalaysis = option; }
  void SetMCClosure(Bool_t option) { this->fMCClosure = option; }
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
    this->fAcceptanceHistogram = AcceptanceHistogram;
  }
  void SetWeightHistogram(TH1F *WeightHistogram) {
    this->fWeightHistogram = WeightHistogram;
  }
  void SetUseWeights(Bool_t option) { this->fUseWeights = option; }
  void SetResetWeights(Bool_t option) { this->fResetWeights = option; }

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

  /* base list holding all output object (a.k.a. grandmother of all lists)
   */
  TList *fHistList;
  TString fHistListName;

  /* control histograms */
  TList *fControlHistogramsList;
  TString fControlHistogramsListName;
  TH1F *fTrackControlHistograms[LAST_ETRACK][LAST_EBEFOREAFTER];
  TString fTrackControlHistogramNames[LAST_ETRACK][LAST_EBEFOREAFTER]
                                     [LAST_ENAME];
  /* array holding bins of track control histograms */
  Double_t fBinsTrackControlHistograms[LAST_ETRACK][LAST_EBINS];
  /* array holding event control histograms */
  TH1F *fEventControlHistograms[LAST_EEVENT][LAST_EBEFOREAFTER];
  TString fEventControlHistogramNames[LAST_EEVENT][LAST_EBEFOREAFTER]
                                     [LAST_ENAME];
  Double_t fBinsEventControlHistograms[LAST_EEVENT][LAST_EBINS];

  /* cuts */
  TString fCentralitySelCriterion;
  Double_t fTrackCuts[LAST_ETRACK][LAST_EMINMAX];
  Double_t fPrimaryVertexCuts[LAST_EXYZ][LAST_EMINMAX];
  Int_t fFilterbit;

  /* Final results */
  TList *fFinalResultsList;
  TString fFinalResultsListName;
  /* array holding final result histograms */
  TH1F *fFinalResultHistograms[LAST_EFINALHIST];
  TString fFinalResultHistogramNames[LAST_EFINALHIST][LAST_ENAME];
  Double_t fBinsFinalResultHistograms[LAST_EFINALHIST][LAST_EBINS];
  TProfile *fFinalResultProfiles[LAST_EFINALPROFILE];
  TString fFinalResultProfileNames[LAST_EFINALPROFILE][LAST_ENAME];
  Double_t fBinsFinalResultProfiles[LAST_EFINALPROFILE][LAST_EBINS];

  /* Monte Carlo analysis */
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

  /* qvectors */
  TList *fQvectorList;
  TComplex fQvector[kMaxHarmonic][kMaxPower];
  std::vector<Double_t> fPhi;
  std::vector<Double_t> fWeights;
  TH1F *fAcceptanceHistogram;
  TH1F *fWeightHistogram;
  Bool_t fUseWeights;
  Bool_t fResetWeights;
  std::vector<std::vector<Int_t>> fCorrelators;

  /* Increase this counter in each new version: */
  ClassDef(AliAnalysisTaskAR, 6);
};

#endif
