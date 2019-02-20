
//--------------------------------------------------------------------------------------// 
// Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.               //
// See cxx source for full Copyright notice                                             //
// $Id$                                                                                 //
//--------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------//
// Analysis task for the data analysis of the correlations between the anisotropic flow //
// harmonics v_n with the Pb-Pb data taken by the ALICE experiment.                     //
// The current script computes the multiparticle correlators using the method of the    //
// Q-vectors for a maximum of 6 different harmonics and 8 particles).                   //
//                                                                                      //
// Author: Cindy Mordasini (cindy.mordasini@cern.ch)                                    //
//--------------------------------------------------------------------------------------//

#ifndef ALIANALYSISTASKTWOMULTICORRELATIONS_H
#define ALIANALYSISTASKTWOMULTICORRELATIONS_H

#include "AliAnalysisTaskSE.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliAODMCParticle.h"
#include "AliVEvent.h"
#include "TProfile.h"
#include "TComplex.h"
#include "TH1D.h"

//======================================================================================//

class AliAnalysisTaskTwoMultiCorrelations : public AliAnalysisTaskSE{
public:
// The six following functions are mandatory for AliAnalysisTaskSE to use the class properly.
  AliAnalysisTaskTwoMultiCorrelations();
  AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights=kFALSE);
  virtual ~AliAnalysisTaskTwoMultiCorrelations();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);

//--------------------------------------------------------------------------------------//
// Setters and getters for the data members.
  void SetAnalysisType(Bool_t bothAnalysis, Bool_t mcAnalysis, Bool_t aodAnalysis)
  {
    this->fProcessBothKineAndReco = bothAnalysis;
    this->fProcessOnlyKine = mcAnalysis;
    this->fProcessOnlyReco = aodAnalysis;
  } // End: SetAnalysisType(Bool_t, Bool_t, Bool_t)

  void SetGeneralParameters(Int_t maxNumberCorrelations, Int_t highestHarmonic, Int_t nbHarmonics, Bool_t useParticleWeights, Bool_t computeNestedLoops)
  {
    this->fMaxNumberCorrelations = maxNumberCorrelations;
    this->fMaxFlowHarmonic = highestHarmonic;
    this->fNumberHarmonicsInSC = nbHarmonics;
    this->fUseParticleWeights = useParticleWeights;
    this->fComputeNestedLoops = computeNestedLoops;
  };  // End: void SetGeneralParameters(Int_t, Bool_t, Bool_t, Bool_t)

  void SetHarmonics(Int_t nOne, Int_t nTwo, Int_t nThree, Int_t nFour, Int_t nFive, Int_t nSix, Int_t nSeven, Int_t nEight)
  {
    this->fHarmonicOne = nOne;
    this->fHarmonicTwo = nTwo;
    this->fHarmonicThree = nThree;
    this->fHarmonicFour = nFour;
    this->fHarmonicFive = nFive;
    this->fHarmonicSix = nSix;
    this->fHarmonicSeven = nSeven;
    this->fHarmonicEight = nEight;
  };  // End: void SetHarmonics(Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t, Int_t)

  void SetCentralityRange(Int_t const nBins, Float_t minCentrality, Float_t maxCentrality)
  {
    this->fCentralityMin = minCentrality;
    this->fCentralityMax = maxCentrality;
  };  // End: void SetCentralityRange(Int_t const, Float_t, Float_t)

  void SetCutsRange(Float_t minPtCut, Float_t maxPtCut, Float_t minEtaCut, Float_t maxEtaCut)
  {
    this->fPtMin = minPtCut;
    this->fPtMax = maxPtCut;
    this->fEtaMin = minEtaCut;
    this->fEtaMax = maxEtaCut;
  };  // End: void SetCutsRange(Float_t, Float_t, Float_t, Float_t)

  void SetPreCutControlList(TList* const pccl) {this->fControlListPreCuts = pccl;};
  TList* GetPreCutControlList() const {return this->fControlListPreCuts;}
  void SetPostCutControlList(TList* const poccl) {this->fControlListPostCuts = poccl;};
  TList* GetPostCutControlList() const {return this->fControlListPostCuts;}
  void SetCorrelationResultsList(TList* const crl) {this->fFinalList = crl;};
  TList* GetCorrelationResultsList() const {return this->fFinalList;}

// Methods called in the constructor.
  virtual void InitialiseArraysOfQvectors();
  virtual void InitialiseArraysOfTProfiles();

// Methods called in UserCreateOutputObjects().
  virtual void BookAndNestAllLists();
  virtual void BookControlListPreCuts();
  virtual void BookControlListPostCuts();
  virtual void BookFinalList();

// Methods called in UserExec(Option_t *).
  virtual void AODanalysis(AliAODEvent *aAODevent);
  virtual void MCanalysis(AliMCEvent *aMCevent);
  virtual void CalculateQvectors(long long nParticles, Double_t angles[], Double_t weights[]);
  virtual void GSCfullAnalysis(long long nParticles, Double_t angles[], Double_t weights[]);
  TComplex Q(Int_t n, Int_t p);
  TComplex CalculateRecursion(Int_t n, Int_t *harmonic, Int_t mult=1, Int_t skip=0);
  Double_t ComputeTwoNestedLoops(long long nParticles, Int_t *harmonic, Double_t angles[], Double_t weights[], TProfile *profile);
  Double_t ComputeThreeNestedLoops(long long nParticles, Int_t *harmonic, Double_t angles[], Double_t weights[], TProfile *profile);
  Double_t ComputeFourNestedLoops(long long nParticles, Int_t *harmonic, Double_t angles[], Double_t weights[], TProfile *profile);

// Methods called in Terminate(Option_t *).


//--------------------------------------------------------------------------------------//
private:
  AliAnalysisTaskTwoMultiCorrelations(const AliAnalysisTaskTwoMultiCorrelations& aattmc);
  AliAnalysisTaskTwoMultiCorrelations& operator=(const AliAnalysisTaskTwoMultiCorrelations& aattmc);

// General parameters.
  Int_t fMaxNumberCorrelations; // Maximum number of particles in the correlator (Not going higher than 8-p correlations).
  Int_t fMaxFlowHarmonic; // Maximum harmonic n for v_n (Not going further than v_6).
  TComplex fQvectors[49][9];  // All needed combinations of Q-vectors (size: [fMaxFlowHarmonic*fMaxNumberCorrelations+1][fMaxNumberCorrelations+1]).
  Int_t fNumberHarmonicsInSC;  // Number of harmonics in the GSC.

  TString *fAnalysisType; //! Type of analysis: MC or AOD.
  Bool_t fProcessBothKineAndReco; // Process both MC and AOD.
  Bool_t fProcessOnlyKine;  // Process only MC.
  Bool_t fProcessOnlyReco;  // Process only AOD.

  Bool_t fUseParticleWeights; // Use non-unit particle weights.
  Bool_t fComputeNestedLoops; // Computate the nested loops for cross-check.

// Structure of the output file.
  TList* fMainList; // Main output list.
  TList* fControlListPreCuts; // Secondary list with the control histograms before the cuts.
  TList* fControlListPostCuts;  // Secondary list with the control histograms after the cuts.
  TList* fFinalList;  // Secondary list with all the results of the correlators.
  TList *fListTwoParticles; // Tertiary list with the 2-p correlations.
  TList *fListThreeParticles; // Tertiary list with the 3-p correlations.
  TList *fListFourParticles;  // Tertiary list with the 4-p correlations.
  TList *fListSixParticles; // Tertiary list with the 6- and 8-p correlations.

// Ranges for the centrality and cuts.
  Double_t fCentralityMin;  // Minimum of the centrality.
  Double_t fCentralityMax;  // Maximum of the centrality.
  Double_t fPtMin;  // Minimum of the transverse momentum.
  Double_t fPtMax;  // Maximum of the transverse momentum.
  Double_t fEtaMin; // Minimum of the pseudorapidity.
  Double_t fEtaMax;  // Maximum of the pseudorapidity.

// Harmonics.
  Int_t fHarmonicOne; // Harmonic n_1.
  Int_t fHarmonicTwo; // Harmonic n_2.
  Int_t fHarmonicThree; // Harmonic n_3.
  Int_t fHarmonicFour;  // Harmonic n_4.
  Int_t fHarmonicFive;  // Harmonic n_5.
  Int_t fHarmonicSix; // Harmonic n_6.
  Int_t fHarmonicSeven; // Harmonic n_7.
  Int_t fHarmonicEight; // Harmonic n_8.

// Control histograms before the application of the cuts.
  TH1D *fHistoCentralityPreCuts;  //! Centrality distribution.
  TH1D *fHistoPtPreCuts;  //! Transverse momentum distribution.
  TH1D *fHistoEtaPreCuts; //! Pseudorapidity distribution.
  TH1D *fHistoPhiPreCuts; //! Azimuthal angles distribution.
  TH1D *fHistoMultiplicityPreCuts;  //! Multiplicity distribution.

// Control histograms after the application of the cuts.
  TH1D *fHistoMultiplicityPostCuts; //! Multiplicity distribution.
  TH1D *fHistoPtPostCuts; //! Transverse momentum distribution.
  TH1D *fHistoEtaPostCuts;  //! Pseudorapidity distribution.
  TH1D *fHistoPhiPostCuts;  //! Azimuthal angles distribution.

// TProfiles for the final results.
  TProfile *fProfileCosineTwoParticles[6];  //! <2>_{j,-j}, j: k,l,m,n,(k+l),(k-l).
  TProfile *fProfileCosineTwoNestedLoops[6];  //! <2>_{j,-j} with nested loops.
  TProfile *fProfileTwoCosine[5]; //! <<2>_{i,-i}<2>_{j,-j}> for ij: kl,km,lm, and <<2>_{i,-i}<2>_{j,-j}<2>_{h,-h}>, ijh: klm and <<4>_{k,l,-k,-l} <2>_{m,-m}>.
  TProfile *fProfileTwoCosineNestedLoops[5];  //! N<<2>_{i,-i}<2>_{j,-j}> with nested loops.
  TProfile *fProfileCosineThreeParticles[4];  //! <3>_{h,i,j}, h: (k+l),(k-l) (first two), j: (l-k),(-k-l) (last two).
  TProfile *fProfileCosineThreeNestedLoops[4];  //! <3>_{h,i,j} with nested loops.
  TProfile *fProfileCosineFourParticles[6]; //! <4>_{i,j,-i,-j}, ij: kl,km,lm,kn,ln,mn.
  TProfile *fProfileCosineFourNestedLoops[6]; //! <4>_{i,j,-i,-j} with nested loops.
  TProfile *fProfileCosineSixParticles[4];  //! <6>_{h,i,j,-h,-i,-j}, hij: klm,kln,kmn,lmn.
  TProfile *fProfileCosineEightParticles; //! <8>_{k,l,m,n,-k,-l,-m,-n}.

// Version counter for the submissions on Grid.
/// Increase the counter by one when the latest version changes the structure
/// of the output file (version of the 2018-11-19, counter: 5)
  ClassDef(AliAnalysisTaskTwoMultiCorrelations,5);
}; // End of the definition of the class AliAnalysisTaskTwoMultiCorrelations

#endif

