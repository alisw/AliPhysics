/* -------------------------------------------------------------------------- /
/ Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.      /
/ See cxx source for full Copyright notice                                    /
/ --------------------------------------------------------------------------- /
/ Analysis task for the computation of the correlators needed for the CFA.    /
/                                                                             /
/ Author: Cindy Mordasini (cindy.mordasini@cern.ch)                           /
/ -------------------------------------------------------------------------- */
#ifndef ALIANALYSISTASKHOCFA_H
#define ALIANALYSISTASKHOCFA_H

#include <sstream>

#include "TComplex.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TList.h"
#include "TProfile.h"
#include "TString.h"
#include "TSystem.h"

#include "AliJEfficiency.h"

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

  // General methods specific to this class.  //TODO.
  virtual void InitialiseArrayMembers();
  virtual void BookFinalResults();
  virtual void CalculateQvectors(Long64_t myMulti,
    double myAngles[], double myWeights[]);
  TComplex Q(int n, int p);
  TComplex CalculateRecursion(int n, int *harmonic, int mult=1, int skip=0);
  virtual void ComputeAllTerms(float myCentWeight);
  virtual void CalculateCorrelator(int myMulti, int myHarmos[],
    TProfile *myProfile, int myBin, int myPowers[], float myCentWeight);
  //virtual void CalculateCorrelator(int combi, int bin, int nParticles, int harmonics[], double *errorTerms);
  virtual void ComputeEtaGaps(Long64_t multiplicity, double angles[], double pWeights[],
    double pseudorapidity[], float myCentWeight);  

  // Setters/getters specific to this class.
  void SetInputList(TClonesArray *inputArray) {fInputList = inputArray;}
  TClonesArray *GetInputList() const {return fInputList;}
  TList *GetCFAlist() const {return fMainList;}
  void SetDebugLevel(int debug) {fDebugLevel = debug;}

  /// Centrality and multiplicity.
  virtual void SetCentralityArray(TString values);
  void SetEventCentrality(float cent) {fCentrality = cent;}
  void SetCentralityBinning(int nBins) {fNCentralityBins = nBins;}
  int GetCentralityBin(float myCent);
  void SetMinMultiplicity(int minMult) {fMultiplicityMin = minMult;}

  /// Kinematic ranges and corrections.
  void SetPtRange(double minPt, double maxPt) {fPtMin = minPt; fPtMax = maxPt;}
  void SetEtaGap(bool etaGap, float myGap) {fApplyEtaGap = etaGap; fEtaGap = myGap;}
  void SetParticleWeights(bool weightsNUE, bool weightsNUA) {
    fUseWeightsNUE = weightsNUE; fUseWeightsNUA = weightsNUA;
  }
  void SetCentralityWeights(bool weightsCent) {fUseWeightsCent = weightsCent;}

  /// Analysis observables.
  void SetObservable(bool thisObs, bool thisOrder) {fGetSC = thisObs; fGetLowerHarmos = thisOrder;
    if (thisObs) {fNCombi = 1;}       // All 13 combinations are obtained in one go for the SC.
    else {
      if (thisOrder) {fNCombi = 7;}   // We get the AC only for the 7 lowest 2-h combinations.
      else {fNCombi = 6;}             // We get the AC only for the 6 higher 2-h combinations.
    }
  }

 private:
  AliAnalysisTaskHOCFA(const AliAnalysisTaskHOCFA& aat);
  AliAnalysisTaskHOCFA& operator=(const AliAnalysisTaskHOCFA& aat);

  TClonesArray *fInputList;     // Input tracks selected in the catalyst.
  TList *fMainList;             // Mother list for all the outputs.
  int fDebugLevel;              // Verbosity of the class (0: default, 5: debug, 10: full verbose).

  TList *fCentralityList[16];   //! Result list per each centrality bin.
  Long64_t fMultiplicity;       //! Multiplicity after full selection.  
  float fCentralityArray[17];   //! Edges for the centrality division. (0-80% with 5% width).
  float fCentrality;            // Centrality of the current event.
  int fNCentralityBins;         //! Number of centrality bins (Size(array)-1).
  int fCentralityBin;           //! Centrality bin of the current event.
  int fMultiplicityMin;         // Minimum multiplicity to have a valid event weight.

  double fPtMin;                // Minimum transverse momentum.
  double fPtMax;                // Maximum transverse momentum.
  float fEtaGap;                // Value of the gap (default: 0.).
  bool fApplyEtaGap;            // kTRUE: Get the 2-p correlators with an eta gap.
  bool fUseWeightsNUE;          // kTRUE: Enable the non-unit NUE corrections.
  bool fUseWeightsNUA;          // kTRUE: Enable the non-unit NUA corrections.
  bool fUseWeightsCent;         // kTRUE: Enable the non-unit centrality correction for LHC15o.

  int fHarmoArray2h[13][2];     // Combinations of 2-harmonics for AC/SC.
    // 13: max number of combinations of 2 harmonics.
  int fHarmoArray3h[9][3];      // Combinations of 3-harmonics for SC only.
  int fPowers[13][2];           // List of 2-h terms by their power for AC only.
  TComplex fQvectors[81][11];   // All the needed combinations of Q-vectors.
    // Size: [(v8*10part)+1][10part+1].
  int fNCombi;                  // Number of combinations of harmonics in one wagon (max 7).
  bool fGetSC;                  // kTRUE: Measure 2-h and 3-h SC, else 2-h AC (default: true).
  bool fGetLowerHarmos;         // kTRUE: Measure the terms for the lower harmonics (default: true).

  TProfile *fHistoConfig;       //! Configuration of the analysis (8 bins).
  TH1F *fHistoCent[16];         //! Centrality distribution of selected tracks.
  TH1F *fHistoCentCorrect[16];  //! Corrected centrality distribution of the selected tracks.
  TH1I *fHistoMulti[16];        //! Multiplicity of the selected tracks.
  TH1D *fHistoPt[16];           //! pT distribution of the selected tracks.
  TH1D *fHistoEta[16];          //! eta distribution of the selected tracks.
  TH1D *fHistoPhi[16];          //! phi distribution of the selected tracks.
  TH1I *fHistoCharge[16];       //! Charge distribution of the selected tracks.

  TProfile *fCorrel2p[16];        //! 2-p terms for v1-v8 without eta gap (8 bins).
  TProfile *fCorrel2pEtaGap[16];  //! 2-p terms for v1-v8 with fixed eta gap (8 bins).
  TProfile *fCorrel2h[7][16];     //! 2-harmonic terms for SC/AC (13 bins).
    // 7: max number of combi in one wagon (all in 1 for SC), 16: max number of centrality bins.
  TProfile *fCorrel3h[16];        //! 3-harmonic terms for SC (9 bins, not enabled for AC).

/*
  TProfile *fErrorTermsSC3h[6][16];   //! Error propagation (SC(k,l,m).
  TProfile *fErrorTermsAC41[6][16];   //! Error propagation (AC_41(m,n)).
*/

  ClassDef(AliAnalysisTaskHOCFA, 7);
};

#endif  // ALIANALYSISTASKHOCFA_H