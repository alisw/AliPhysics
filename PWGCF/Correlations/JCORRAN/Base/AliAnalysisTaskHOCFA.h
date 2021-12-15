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
#include "TList.h"
#include "TString.h"
#include "TComplex.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TProfile.h"
#include "AliJEfficiency.h"
#include <sstream>

class TClonesArray;

class AliAnalysisTaskHOCFA {
 public:
// Methods inherited from AliAnalysisTaskSE.
  AliAnalysisTaskHOCFA();
  AliAnalysisTaskHOCFA(const char *name);
  virtual ~AliAnalysisTaskHOCFA();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *);

// Methods specific to this class.
  void SetInputList(TClonesArray *inputArray) {fInputList = inputArray;}
  TClonesArray *GetInputList() const {return fInputList;}
  void SetDebugLevel(int debug) {fDebugLevel = debug;}
  TList *GetCFAlist() const {return fMainList;}
  virtual void InitialiseArrayMembers();

  void SetCentralityBinning(int nBins) {fNCentralityBins = nBins;}
  virtual void SetCentralityArray(TString values);
  void SetEventCentrality(float cent) {fCentrality = cent;}
  int GetCentralityBin(float cent);

  virtual void BookFinalResults();

  void SetMinMultiplicity(int minMult) {fMultiplicityMin = minMult;}
  void SetPtRange(double minPt, double maxPt) {fPtMin = minPt; fPtMax = maxPt;}
  void SetParticleWeights(bool weightsNUE, bool weightsNUA) {
    fUseWeightsNUE = weightsNUE; fUseWeightsNUA = weightsNUA;
  }
  void SetObservable(bool observ) {fGetSC3h = observ;}
  void SetNumberCombi(int combi) {fNCombi = combi;}
  virtual void SetHarmoArray(TString combiString);
  void SetEtaGaps(bool etaGap, float myGap) {
    fGetEtaGap = etaGap; fEtaGap = myGap;
  }

  virtual void CalculateQvectors(Long64_t multiplicity, double angles[], double pWeights[]);
  TComplex Q(int n, int p);
  TComplex CalculateRecursion(int n, int *harmonic, int mult=1, int skip=0);
  virtual void ComputeAllTerms();
  virtual void CalculateCorrelator(int combi, int bin, int nParticles, int harmonics[], double *errorTerms);
  virtual void ComputeEtaGaps(Long64_t multiplicity, double angles[], double pWeights[], double pseudorapidity[]);

private:
  AliAnalysisTaskHOCFA(const AliAnalysisTaskHOCFA& aat);
  AliAnalysisTaskHOCFA& operator=(const AliAnalysisTaskHOCFA& aat);
  TClonesArray *fInputList;     // List of selected tracks.
  int fDebugLevel;              // Select how much is printed in the terminal.
  TList *fMainList;             // Main list for all the results.
  TList *fCentralityList[16];   //! Results per each centrality bin.

  int fNCentralityBins;         //! Number of centrality bins (Size(array)-1).
  float fCentralityArray[17];   //! Edges for the centrality division.
    // (0-80% with 5% width).
  float fCentrality;            // Centrality of the current event.
  int fCentralityBin;           //! Corresponding centrality bin.

  Long64_t fMultiplicity;       //! Multiplicity after full selection.
  int fMultiplicityMin;         // Minimum multiplicity to have valid events.
  double fPtMin;                // Minimum transverse momentum.
  double fPtMax;                // Maximum transverse momentum.
  bool fUseWeightsNUE;          // kTRUE: Enable the non-unit NUE corrections.
  bool fUseWeightsNUA;          // kTRUE: Enable the non-unit NUA corrections.
  bool fGetSC3h;                // kTRUE: Calculate SC(k,l,m), else AC(m,n).
  bool fGetEtaGap;              // kTRUE: Get the 2p correlators with an eta gap.
  float fEtaGap;                // Value of the gap (default: 0.).
  int fNCombi;                  // Number of combinations of harmonics (max 6).
  int fHarmoArray[6][3];        // Combinations of harmonics for the CFA.
    // 6: max number of possible observables, 3: number of harmonics.
  int fPowers[15][3];           // List of terms by their power.
    // 15 terms, 3 harmonics.
  TComplex fQvectors[81][11];   // Needed combinations of Q-vectors.
    // Size: [(v8*10part)+1][10part+1].

  TProfile *fHistoConfig;       //! Configuration of the analysis.
  TH1F *fHistoCent[16];         //! Centrality distribution of trimmed tracks.
  TH1I *fHistoMulti[16];        //! Multiplicity of the trimmed tracks.
  TH1D *fHistoPt[16];           //! pT distribution of the trimmed tracks.
  TH1D *fHistoEta[16];          //! eta distribution of the trimmed tracks.
  TH1D *fHistoPhi[16];          //! phi distribution of the trimmed tracks.
  TH1I *fHistoCharge[16];       //! Charge distribution of the trimmed tracks.
  TProfile *fCorrelTerms[6][16];      //! Combinations of correlators for CFA.
    // 6: Max number of SCs/ACs, 16: Max number of centrality bins.
  TProfile *fCorrelEtaGap[16];        //! 2-particle correlators for v_1-v_8 with eta gap.
  TProfile *fErrorTermsSC3h[6][16];   //! Error propagation (SC(k,l,m).
  TProfile *fErrorTermsAC41[6][16];   //! Error propagation (AC_41(m,n)).

  ClassDef(AliAnalysisTaskHOCFA, 5);
};

#endif