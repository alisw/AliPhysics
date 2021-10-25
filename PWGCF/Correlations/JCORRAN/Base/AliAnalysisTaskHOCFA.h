#ifndef ALIANALYSISTASKHOCFA_H
#define ALIANALYSISTASKHOCFA_H

/* -------------------------------------------------------------------------- /
/ Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.      /
/ See cxx source for full Copyright notice                                    /
/ --------------------------------------------------------------------------- /
/ Analysis task for the computation of the correlators needed for the CFA.    /
/                                                                             /
/ Author: Cindy Mordasini (cindy.mordasini@cern.ch)                           /
/ -------------------------------------------------------------------------- */

#include "TSystem.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TList.h"
#include "TComplex.h"
#include "AliJEfficiency.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TProfile.h"
#include <TExMap.h>
#include <sstream>

class TClonesArray;

class AliAnalysisTaskHOCFA {
public:
// Methods inherited from AliAnalysisTaskSE.
  AliAnalysisTaskHOCFA();
  AliAnalysisTaskHOCFA(const char *name, Bool_t useWeights=kFALSE);
  virtual ~AliAnalysisTaskHOCFA();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

// Methods specific for this class.
  void SetInputList(TClonesArray *inputArray) {fInputList = inputArray;}
  TClonesArray *GetInputList() const {return fInputList;}
  void SetDebugLevel(Int_t debug) {fDebugLevel = debug;}
  TList *GetCFAlist() const {return fMainList;}
  virtual void InitialiseArrayMembers();

  void SetCentralityBinning(Int_t nBins) {fNCentralityBins = nBins;}
  virtual void SetCentralityArray(TString values);
  void SetEventCentrality(Float_t cent) {fCentrality = cent;}
  Int_t GetCentralityBin(Float_t cent);

  virtual void BookFinalResults();

  void SetMinMultiplicity(Int_t minMult) {fMultiplicityMin = minMult;}
  void SetParticleWeights(Bool_t weightsNUE, Bool_t weightsNUA) {fUseWeightsNUE = weightsNUE; fUseWeightsNUA = weightsNUA;}
  void SetNumberCombi(Int_t combi) {fNCombi = combi;}
  virtual void SetHarmoArray(TString combiString);

  virtual void CalculateQvectors(Long64_t multiplicity, Double_t angles[], Double_t pWeights[]);
  TComplex Q(Int_t n, Int_t p);
  TComplex CalculateRecursion(Int_t n, Int_t *harmonic, Int_t mult=1, Int_t skip=0);
  virtual void ComputeAllTerms(); // TBC: Do I need to pass the angles and weights?
  virtual void CalculateCorrelator(Int_t combi, Int_t bin, Int_t nParticles, Int_t harmonics[], Double_t *errorTerms);


private:
  AliAnalysisTaskHOCFA(const AliAnalysisTaskHOCFA& aat);
  AliAnalysisTaskHOCFA& operator=(const AliAnalysisTaskHOCFA& aat);
  TClonesArray *fInputList; // List of selected tracks.
  Int_t fDebugLevel;  // Choose the quantity of terminal outputs for debugging.
  TList *fMainList; // Main TList where all results for this analysis are saved.
  TList *fCentralityList[16]; //! TLists for the results per each centrality bin.

  Int_t fNCentralityBins; //! Number of centrality bins in the division (Size(array)-1).
  Float_t fCentralityArray[17]; //! Edges for the centrality division (0-80% with 5% width).
  Float_t fCentrality;  // Centrality of the current event.
  Int_t fCentralityBin; //! Corresponding centrality bin.

  Long64_t fMultiplicity; //! Multiplicity of the event after full selection.
  Int_t fMultiplicityMin; // Minimum multiplicity to calculate the correlators.
  Bool_t fUseWeightsNUE; // kTrue: Enable the non-unit NUE corrections.
  Bool_t fUseWeightsNUA; // kTrue: Enable the non-unit NUA corrections.
  Int_t fNCombi;  // Number of combinations of harmonics (max 6).
  Int_t fHarmoArray[6][3];  // Combinations of harmonics for the SCs/ACs.
    // 6: max number of possible observables, 3: number of harmonics.
  Int_t fPowers[15][3]; // List of terms by their power.
    // 15 terms, 3 harmonics.
  TComplex fQvectors[81][11]; // Needed combinations of Q-vectors.
    // Size: [(v8*10part)+1][10part+1].

  TH1F *fHistoCent[16]; //! Centrality distribution of the trimmed tracks.
  TH1I *fHistoMulti[16];  //! Multiplicity distribution of the trimmed tracks.
  TH1D *fHistoPt[16]; //! pT distribution of the trimmed tracks.
  TH1D *fHistoEta[16];  //! eta distribution of the trimmed tracks.
  TH1D *fHistoPhi[16];  //! phi distribution of the trimmed tracks.
  TH1I *fHistoCharge[16]; //! Charge distribution of the trimmed tracks.
  TProfile *fCorrelTerms[6][16];  //! Combinations of correlators for SCs/ACs.
    // 6: Max number of SCs/ACs, 16: Max number of centrality bins.
  TProfile *fErrorTerms[6][16]; //! Terms for the error propagation for NSC(k,l,m).

  ClassDef(AliAnalysisTaskHOCFA, 3);
};

#endif
