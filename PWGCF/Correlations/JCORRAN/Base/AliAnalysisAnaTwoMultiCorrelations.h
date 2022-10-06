
/* -------------------------------------------------------------------------- /
/ Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.      /
/ See cxx source for full Copyright notice                                    /
/ $Id$                                                                        /
/                                                                             /
/ Analysis task computing the 2-, 4- and 6-particle correlators for different /
/ combinations of flow amplitudes up to v_6. The script can take as an input  /
/ Monte Carlo simulations at reco and kine levels (e.g. HIJING) as well as    /
/ experimental Pb-Pb data.                                                    /
/                                                                             /
/ Author: Cindy Mordasini (cindy.mordasini@cern.ch)                           /
/ Version 20 from the 16.09.2020.                                             /
/ -------------------------------------------------------------------------- */

#ifndef AliAnalysisAnaTwoMultiCorrelations_H
#define AliAnalysisAnaTwoMultiCorrelations_H

#include "TSystem.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TList.h"
#include "TComplex.h"
#include "AliJEfficiency.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TProfile.h"
#include <TExMap.h>
class TClonesArray;

class AliAnalysisAnaTwoMultiCorrelations{// : public AliAnalysisTaskSE

public:
/* Mandatory functions for the class to work properly within the framework. ---------------- */
  AliAnalysisAnaTwoMultiCorrelations();
  AliAnalysisAnaTwoMultiCorrelations(const char *name, Bool_t useParticleWeights=kFALSE);
  virtual ~AliAnalysisAnaTwoMultiCorrelations();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);

/* Setters and getters for the data members.                                                 */
// 1. Configuration of the general parameters of the analysis.
  void SetAnalysisParameters(Bool_t getAC,
    Int_t harmoK, Int_t harmoL, Int_t harmoM)
  {
    this->fComputeACs     = getAC;
    this->fACHarmoOne     = harmoK;
    this->fACHarmoTwo     = harmoL;
    this->fACHarmoThree   = harmoM;
  } // End: void SetAnalysisParameters(Bool_t, Bool_t, Bool_t, Int_t, Int_t, Int_t).
  
  void SetInputList(TClonesArray *inputarray){fInputList = inputarray;}
  TClonesArray * GetInputList() const{return fInputList;}
  void SetEventCentrality( float cent ){fCentrality = cent;}
  Int_t GetCentralityBin(Float_t fFstCentrality);
  AliJEfficiency* GetAliJEfficiency() const{return fEfficiency;}
  TList* GetTList() const{return fMPCList;} // get the list for external task
// 2. Configuration of the centrality parameters.
  void SetCentrality(Int_t const nBins, Float_t minCentrality, Float_t maxCentrality,
    Int_t totalBins)
  {
    this->fCentralityMin        = minCentrality;
    this->fCentralityMax        = maxCentrality;
    this->fTotalCentralityBin   = totalBins;
  } // End: void SetCentrality(Float_t, Float_t, Bool_t, Bool_t).

// 6. Configuration of the use of the weights.
  void SetParticleWeights(Bool_t useTable, Bool_t useNonUnitWeights,
    Bool_t usePt, Bool_t usePhi, Bool_t useEta)
  {
    this->fUseParticleWeights   = useNonUnitWeights;
    this->fUsePtWeights         = usePt;
    this->fUsePhiWeights        = usePhi;
    this->fUseEtaWeights        = useEta;
  } // End: void SetParticleWeights(Bool_t, Bool_t, Bool_t, Bool_t, Bool_t).

  void SetInputParticleWeights(TString fileWeight);

  void SetJWeights(Bool_t useJEfficiency, Int_t indexFilter)
  {
    this->fUseJEfficiency   = useJEfficiency;
    this->fFilterbitIndex   = indexFilter;
  } // End: void SetJWeights(Bool_t, Int_t).

// 7. Configuration of the histograms and results.
  void SetReducedQvectors(Int_t maxHarmo, Int_t maxCorrel, Int_t powerK, Bool_t doEtaGaps)
  {
    this->fHighestHarmonic      = maxHarmo;
    this->fLargestCorrelators   = maxCorrel;
    this->fReducedQPower        = powerK;
    this->fComputeEtaGaps       = doEtaGaps;
  } // End: void SetReducedQvectors(Int_t, Int_t, Int_t).

  void SetBinningEvents(Int_t nBinsMulti)
  {
    this->fNumberBinsMulti    = nBinsMulti;
  } // End: void SetBinningEvents(Int_t, Int_t, Float_t).

/* Methods called in the constructors. */
  virtual void  InitialiseArraysOfDataMembers();

/* Methods called in "UserExec". */
  virtual void  AnalyseRecoEvent(); // Do the normal analysis at reconstructed level.
  virtual void  CalculateQvectors(long long numberOfParticles, Float_t angles[], Float_t pWeights[]);
  TComplex      Q(Int_t n, Int_t p);
  virtual void  ComputeReducedQvectors(long long numberOfParticles);
  TComplex      CalculateRecursion(Int_t n, Int_t *harmonic, Int_t mult=1, Int_t skip=0);
  virtual void  ComputeSCsCorrelators(long long numberOfParticles, Float_t angles[], Float_t pWeights[]);
  virtual void  ComputeACsCorrelators(long long numberOfParticles, Float_t angles[], Float_t pWeights[]);
  virtual void  ComputeTPCWithEtaGaps(long long numberOfParticles, Float_t angles[], Float_t pWeights[], Float_t pseudorapidity[]);
  void SetDebugLevel(int debuglevel){
    fDebugLevel = debuglevel; cout <<"setting Debug Level = " << fDebugLevel << endl;}
/* Methods called in "Terminate"............................................................ */

/* Methods called in "UserCreateOutputObjects".............................................. */
  virtual void  BookAllLists();
  virtual void  BookMPCList();
/* ----------------------------------------------------------------------------------------- */
private:
  AliAnalysisAnaTwoMultiCorrelations(const AliAnalysisAnaTwoMultiCorrelations& aattmc);
  AliAnalysisAnaTwoMultiCorrelations& operator=(const AliAnalysisAnaTwoMultiCorrelations& aattmc);
  TClonesArray *fInputList;
// 1. General parameters for the configuration of the analysis.
  Int_t fACHarmoOne;        // First harmonic in ACs with two and three harmonics.
  Int_t fACHarmoTwo;        // Second harmonic in ACs with two and three harmonics.
  Int_t fACHarmoThree;      // Third harmonic in ACs with three harmonics.
    // kFALSE: compute the multi-particle correlators at reco level.
  Bool_t fComputeACs;       // kTRUE: get the correlators needed for ACs.
  Bool_t fWriteMinimum;

// 2. Parameters related to the centrality.
  Int_t fNumberBinsMulti;         // Number of bins for fHistoMultiplicity.
  Int_t fTotalCentralityBin;      // Total number of centrality bins in the analysis range.
  Float_t fCentrality;
  Int_t fCentralityBin; 

  long long fInitialMultiplicity; //! Initial number of tracks in event.
  Float_t fCentralityMin;         // Minimum of the centrality analysis range.
  Float_t fCentralityMax;         // Maximum of the centrality analysis range.  

// 4. Parameters related to the HMOs selection.
  Int_t fMultiplicityMin;       // Minimum multiplicity needed for the event weight.
    // (Minimum non-strict -> The event can still have this multiplicity and be selected.)
// 5. Parameters related to the track selection

// 6. Parameters related to the efficiency and acceptance weights.
  Bool_t fUseParticleWeights; // kTRUE: use non-unit particle weight.
  Bool_t fUsePtWeights;       // kTRUE: use pT weights for NUE.
  Bool_t fUsePhiWeights;      // kTRUE: use phi weights for NUA.
  Bool_t fUseEtaWeights;      // kTRUE: use eta weights for NUA.

  AliJEfficiency *fEfficiency;  // Used to apply NUE to the data.
  Bool_t fFirstEvent;           ///< True if this is the first event analyzed.
  Bool_t fUseJEfficiency;       // Use JEfficiency code to get the pT-efficiency?
  Int_t fFilterbitIndex;        // Index used for the efficiency correction, must correspond to the main filter.
    // 0: TPCOnly; 6: hybrid (This work for AOD86, I am not sure if it is work for new AOD).
  TH1F *fHistoEfficiency;   //! Distribution of the efficiency correction.
  TH1F *fHistoEffInverse;   //! Distribution of the inverse of the efficiency correction.

// 7. Parameters related to the multi-particle correlations.
  Int_t fHighestHarmonic;       // Largest order of flow amplitude to compute (default: 8).
  Int_t fLargestCorrelators;    // Maximum number of particles in the correlators (default: 10).
  Int_t fReducedQPower;         // Power k for the reduced Q-vectors (default: 0).
  TComplex fQvectors[81][11];    // All the needed combinations of Q-vectors.
    // Size: [(fHighestHarmonic*fLargestCorrelators)+1][fLargestCorrelators+1].
  TList *fMPCList;              //! Daughter list for the multi-particle correlations techniques.
  TH1F *fHistoReducedQvectors[9][8];     //! Modulus of the reduced Q-vectors distributions for a given k.
  TProfile *fProfileTwoPartCorrel[9];    //! 2-particle correlators (SC: 8 bins, AC: 3 bins).
  TProfile *fProfileFourPartCorrel[9];   //! 4-particle correlators (SC: 36 bins, AC: 4 bins).
  TProfile *fProfileFourPartCorrelCheck[9]; //! <4>_{j,k,-j,-k} for j,k = 1..8 cross-check (28 bins).
  TProfile *fProfileSixPartCorrel[9];    //! 6-particle correlators (SC: 56 bins, AC: 4 bins).
  TProfile *fProfileEightPartCorrel[9];  //! 8-particle correlators (No SC, AC: 3 bins).
  TProfile *fProfileTenPartCorrel[9];    //! 10-particle correlators (No SC, AC: 1 bin).

// 8. Parameters related to the 2-particle correlations with eta gaps.
  Bool_t fComputeEtaGaps;   // kTRUE: compute the 2-particle correlator with eta gaps.
  TList *fTPCEtaList;       //! Daughter list with the 2-p correlators calculated with eta gaps.
  TProfile *fProfileTPCEta[9][11]; //! <2>_{j,-j} for j = 1..8 with eta gaps (8 bins per profile).
  int fDebugLevel; //
/* ----------------------------------------------------------------------------------------- */
/* Version number to handle the objects through the iterations of the code.                  */
  ClassDef(AliAnalysisAnaTwoMultiCorrelations, 20);
};  // End of the class.

#endif


