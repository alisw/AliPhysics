/* -------------------------------------------------------------------------- /
/ Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.      /
/ See cxx source for full Copyright notice                                    /
/ --------------------------------------------------------------------------- /
/ Analysis task for the computation of the numerator and denominators for the /
/ various SPC measured for Run2 Pb-Pb data                                    /
/                                                                             /
/ Authors: Cindy Mordasini (cindy.mordasini@cern.ch)                          /
/          Maxim Virta                                                        /
/ -------------------------------------------------------------------------- */
#ifndef ALIANALYSISSPCRUN2_H
#define ALIANALYSISSPCRUN2_H

#include "AliHeader.h"
#include "TComplex.h"
#include "TList.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TProfile.h"
#include "TSystem.h"
#include "AliJBaseTrack.h"
#include "TClonesArray.h"

#include "AliJEfficiency.h"

class TClonesArray;
class AliAnalysisSPCRun2 {
 public:
  // Methods inherited from AliAnalysisTaskSE.
  AliAnalysisSPCRun2();
  AliAnalysisSPCRun2(const char *name);
  virtual ~AliAnalysisSPCRun2();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  // General methods specific to this analysis class.
  virtual void InitializeArrays();
  virtual void BookAndNestAllLists();
  virtual void BookControlHistograms();
  virtual void BookFinalResultsHistograms();
  virtual void CalculateQvectors(Int_t c_nParticles, Double_t* c_angles, Double_t* c_weights);
  TComplex Q(Int_t n, Int_t p);
  TComplex Recursion(Int_t n, Int_t* harmonic, Int_t mult, Int_t skip);
  virtual void MainTask(Int_t centBin, Int_t mult, Double_t* m_angles, Double_t* m_weights);
  virtual void Correlation(Int_t c_nPart, Int_t c_nHarmo, Int_t* harmo, Double_t *correlData);
  virtual void ComputeTPCWithEtaGaps(Int_t centBin, Int_t mult, Double_t* m_angles, Double_t* m_weights, Double_t* m_pseudo);

  /// General setters/getters.
  void SetInputList(TClonesArray *inputarray) {fInputList = inputarray;}
  TClonesArray *GetInputList() const {return fInputList;}
  TList* GetMainList() const{return fHistList;}
  void SetDebugLevel(Int_t debuglevel) {
    fDebugLevel = debuglevel;
    cout << "Setting Debug Level = " << fDebugLevel << endl;
  }
  void SetSaveAllQA(Bool_t SaveQA) {bSaveAllQA = SaveQA;}

  /// Centrality-related methods.
  void SetEventCentrality(Double_t cent) {fCentrality = cent;}
  Int_t SelectCentrality(Double_t centValue);
  void SetCentrality(Float_t cen0, Float_t cen1, Float_t cen2, Float_t cen3, Float_t cen4,
      Float_t cen5, Float_t cen6, Float_t cen7, Float_t cen8, Float_t cen9, Float_t cen10,
      Float_t cen11, Float_t cen12, Float_t cen13, Float_t cen14, Float_t cen15, Float_t cen16) {
    fcent_0 = cen0; fcent_1 = cen1; fcent_2 = cen2; fcent_3 = cen3; fcent_4 = cen4;
    fcent_5 = cen5; fcent_6 = cen6; fcent_7 = cen7; fcent_8 = cen8; fcent_9 = cen9;
    fcent_10 = cen10; fcent_11 = cen11; fcent_12 = cen12; fcent_13 = cen13; fcent_14 = cen14;
    fcent_15 = cen15; fcent_16 = cen16;
  }
  virtual void SetInitializeCentralityArray();

  void SetMinNuPar(Int_t top) {fMinNumberPart = top;} 
  Int_t GetMinNuPar() const {return fMinNumberPart;}
  void SetEtaGaps(Bool_t ComputeEtaGap, Float_t EtaGap) {
    bComputeEtaGap = ComputeEtaGap;
    fEtaGap = EtaGap;
  }

  /// Correlator-related methods.
  void SetCorrSet(Int_t obsInd, Int_t harmo[8]) {
    for (int i = 0; i < 8; i++) {fHarmosArray[obsInd][i] = harmo[i];}
  }

  // LOKI: Added in HIJING investigations.
 void SetUseWeights(Bool_t WeightsNUE, Bool_t WeightsNUA){this->bUseWeightsNUE = WeightsNUE; this->bUseWeightsNUA = WeightsNUA;}

 private:
  AliAnalysisSPCRun2(const AliAnalysisSPCRun2& aat);
  AliAnalysisSPCRun2& operator=(const AliAnalysisSPCRun2& aat);

  TClonesArray *fInputList;         // Input tracks selected in the JCatalyst.
  TList *fHistList;                 // Base list to hold all output objects.
  Int_t fDebugLevel;                // Verbosity of the class in the terminal.

  Double_t fCentrality;             // Centrality of the current event.
  Float_t fcent_0, fcent_1, fcent_2, fcent_3, fcent_4, fcent_5, fcent_6, fcent_7, fcent_8,
    fcent_9, fcent_10, fcent_11, fcent_12, fcent_13, fcent_14, fcent_15, fcent_16;
      // fcent_i holds the edge of a centrality bin.
  Int_t fCentralityBins;            //! Number of centrality bins (Size(array)-1).
  Int_t fMinNumberPart;             // Minimum number of particles to get valid correlators.

  Bool_t bUseWeightsNUE;            // kTRUE: Use non-unit particle weights for NUE corrections.
  Bool_t bUseWeightsNUA;            // kTRUE: Use non-unit particle weights for NUA corrections.
  Bool_t bComputeEtaGap;            // kTRUE: Calculate 2-particle eta gaps (default: kFALSE).
  Float_t fEtaGap;                  // Value of the eta gap itself.

  Bool_t bSaveAllQA;                // kTRUE: Save the standard QA histograms (default: kTRUE).
  TH1F *fCounterHistogram;          //! Counters for some QA checks.
  TProfile *fProfileTrackCuts;      //! Storage of the values used in some cuts.

  TComplex fQvector[113][15];       // All combinations of Q-vectors.
  Int_t fHarmosArray[12][8];        // Array of combinations of harmonics for the SPC.
  
  TList *fCentralityList[16];       //! Results per centrality bins. Up to 16 possible bins.
    // Size: [fMaxHarmonic*fMaxCorrelator+1][fMaxCorrelator+1]
    // Can deal with maximum 10 different SPC.
  TList *fControlHistogramsList[16];  //! List to hold all control histograms for a centrality bin.
  TList *fFinalResultsList[16];     //! List to hold all histograms with final results.

  TH1F *fCentralityHistogram[16];   //! Centrality of the final tracks.
  TH1F *fMultHistogram[16];         //! Number of the final tracks per event.
  TH1F *fPTHistogram[16];           //! Transverse momentum of the final tracks.
  TH1F *fPhiHistogram[16];          //! Azimuthal angles of the final tracks.
  TH1F *fEtaHistogram[16];          //! Pseudorapidity of the final tracks.
  TH1I *fChargeHistogram[16];       //! Electric charge after the track selection.
  TProfile *fPhiWeightProfile[16];  //! QA for the NUA weights.

  TProfile *fResults[16];           //! Final numerators and denominators.
  TProfile *fResultsAlternativeError[16]; //! 'fResults' with different error options
  TProfile *fCovResults[16];        //! Storage of the terms needed for the covariance.
  TProfile *fJoinedCovResults[16];  //! Storage of the joined covariance term calculated as one
                                    //  correlator <z> instead of product of two correlators <x*y>.
  TProfile *fProfileTPCEta[16];     //! Profile for 2-particle eta gap computation.
  Float_t fCentralityArray[17];     //! Edges for the centrality division. (0-80% with 5% width).

  ClassDef(AliAnalysisSPCRun2, 1); 
};

#endif  // ALIANALYSISSPCRUN2
