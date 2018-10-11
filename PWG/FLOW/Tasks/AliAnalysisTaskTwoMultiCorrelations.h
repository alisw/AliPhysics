
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
#include "AliVEvent.h"
#include "TProfile.h"
#include "TComplex.h"
#include "TH1D.h"

//======================================================================================//

class AliAnalysisTaskTwoMultiCorrelations : public AliAnalysisTaskSE{
public:
// The six following functions are mandatory for AliAnalysisTaskSE to use the class properly
  AliAnalysisTaskTwoMultiCorrelations();
  AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights=kFALSE);
  virtual ~AliAnalysisTaskTwoMultiCorrelations();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);


// Organisation of the methods
/// Setters and getters for the data members
/// Methods called in the constructor
/// Methods called in UserCreateOutputObjects()
/// Methods called in UserExec(Option_t *)
/// Methods called in Terminate(Option_t *)

// Setters and getters for the data members
  void SetGeneralParameters(Int_t maxNumberCorrelations, Int_t highestHarmonic, Int_t nbHarmonics, Bool_t useParticleWeights, Bool_t computeNestedLoops, Bool_t computeSine)
  {
    this->fMaxNumberCorrelations = maxNumberCorrelations;
    this->fHighHarmonic = highestHarmonic;
    this->fNumberDifferentHarmonics = nbHarmonics;
    this->fUseParticleWeights = useParticleWeights;
    this->fComputeNestedLoops = computeNestedLoops;
    this->fComputeSine = computeSine;
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
    this->fMinCentrality = minCentrality;
    this->fMaxCentrality = maxCentrality;
  };  // End: void SetCentralityRange(Int_t const, Float_t, Float_t)

  void SetCutsRange(Float_t minPtCut, Float_t maxPtCut, Float_t minEtaCut, Float_t maxEtaCut)
  {
    this->fMinPtCut = minPtCut;
    this->fMaxPtCut = maxPtCut;
    this->fMinEtaCut = minEtaCut;
    this->fMaxEtaCut = maxEtaCut;
  };  // End: void SetCutsRange(Float_t, Float_t, Float_t, Float_t)

  void SetPreCutControlList(TList* const pccl) {this->fPreCutControlList = pccl;};
  TList* GetPreCutControlList() const {return this->fPreCutControlList;}
  void SetPostCutControlList(TList* const poccl) {this->fPostCutControlList = poccl;};
  TList* GetPostCutControlList() const {return this->fPostCutControlList;}
  void SetCorrelationResultsList(TList* const crl) {this->fCorrelationResultsList = crl;};
  TList* GetCorrelationResultsList() const {return this->fCorrelationResultsList;}

// Methods called in the constructor
  virtual void InitialiseArraysOfQvectors();
  virtual void InitialiseArraysOfTProfiles();

// Methods called in UserCreateOutputObjects()
  virtual void BookAndNestAllLists();
  virtual void BookPreCutControlList();
  virtual void BookPostCutControlList();
  virtual void BookCorrelationResultsList();

// Methods called in UserExec(Option_t *)
  virtual void CalculateQvectors(Int_t nParticles, Double_t phi[], Double_t particleWeight[]);
  TComplex Q(Int_t n, Int_t p);
  TComplex CalculateRecursionWithQvectors(Int_t nPartCorr, Int_t harmonics[], Int_t p=1, Int_t skip=0);
  virtual void ComputeTwoParticleCorrelationWithQvectors(Int_t numeratorFirstTwoParticleHarmonics[], Int_t numeratorLastTwoParticleHarmonics[], Int_t indexTProfile);
  virtual void ComputeFourParticleCorrelationWithQvectors(Int_t numeratorFourParticleHarmonics[], Int_t indexTProfile);
  virtual void ComputeSixParticleCorrelationWithQvectors(Int_t numeratorSixParticleHarmonics[], Int_t indexTProfile);
  virtual void ComputeEightParticleCorrelationWithQvectors(Int_t numeratorEightParticleHarmonics[]);
  virtual void ComputeCorrelationsWithTwoNestedLoops(Int_t numeratorTwoParticleHarmonics[], Int_t numeratorTwoLastParticleHarmonics[], Int_t indexTProfile, Int_t nParticles, Double_t phi[], Double_t particleWeight[]);
  virtual void ComputeCorrelationsWithFourNestedLoops(Int_t numeratorFourParticleHarmonics[], Int_t indexTProfile, Int_t nParticles, Double_t phi[], Double_t particleWeight[]);

// Methods called in Terminate(Option_t *)


private:
  AliAnalysisTaskTwoMultiCorrelations(const AliAnalysisTaskTwoMultiCorrelations& aattmc);
  AliAnalysisTaskTwoMultiCorrelations& operator=(const AliAnalysisTaskTwoMultiCorrelations& aattmc);

// Organisation of the data members
/// General parameters for the analysis
/// Harmonics
/// Ranges for the centrality and cuts
/// Structure of the output
/// Histograms and TProfiles

// General parameters
  Int_t fMaxNumberCorrelations;  // Maximum number of particles in the correlator (not going than 4-harmonic 8-particle symmetric cumulants)
  Int_t fHighHarmonic; // Maximum harmonic n for v_n (not going higher than v6)
  TComplex fQvectors[49][9];  // All needed combinations of Q-vectors (size: [fHighHarmonic*fMaxNumberCorrelations+1][fMaxNumberCorrelations+1])
  Int_t fNumberDifferentHarmonics;  // Number of different harmonics (2-4)
  Bool_t fUseParticleWeights; // kTRUE: enable the use of non-unit particle weights
  Bool_t fComputeNestedLoops; // kTRUE: enable the computation of the nested loops for validation
  Bool_t fComputeSine;  // kTRUE: enable the computation and save of the sine component

// Harmonics
  Int_t fHarmonicOne; // Harmonic n_1
  Int_t fHarmonicTwo; // Harmonic n_2
  Int_t fHarmonicThree; // Harmonic n_3
  Int_t fHarmonicFour;  // Harmonic n_4
  Int_t fHarmonicFive;  // Harmonic n_5
  Int_t fHarmonicSix; // Harmonic n_6
  Int_t fHarmonicSeven; // Harmonic n_7
  Int_t fHarmonicEight; // Harmonic n_8

// Ranges for the centrality and cuts
  Float_t fMinCentrality; // Minimum of the centrality
  Float_t fMaxCentrality; // Maximum of the centrality
  Float_t fMinPtCut;  // Minimum of the range for the transverse momentum
  Float_t fMaxPtCut;  // Maximum of the range for the transverse momentum
  Float_t fMinEtaCut; // Minimum of the range for the pseudorapidity
  Float_t fMaxEtaCut; // Maximum of the range for the pseudorapidity

// Structure of the output file
  TList* fOutputMainList; // Main list holding all the output objects
  TList* fPreCutControlList;  // List holding all the control histograms before the cuts
  TList* fPostCutControlList; // List holding all the control histograms after the cuts
  TList* fCorrelationResultsList;  // List holding all the multiparticle correlations

// Control histograms before the application of the cuts
  TH1D *fCentralityPreCutHisto; //! Centrality distribution
  TH1D *fMultiplicityPreCutHisto;  //! Multiplicity distribution
  TH1D *fPtPreCutControlHisto;  //! Transverse momentum distribution
  TH1D *fEtaPreCutControlHisto; //! Pseudorapidity distribution
  TH1D *fPhiPreCutHisto; //! Azimuthal angles distribution

// Control histograms after the application of the cuts
  TH1D *fMultiplicityPostCutHisto;  //! Multiplicity distribution
  TH1D *fPtPostCutControlHisto;  //! Transverse momentum distribution
  TH1D *fEtaPostCutControlHisto; //! Pseudorapidity distribution
  TH1D *fPhiPostCutHisto; //! Azimuthal angles distribution

// Results histograms for the multiparticle correlations
  TProfile *fTwoParticleCorrelationProfile[2][4]; //! 2-p correlations ([cos,sin][m,n,p,q]) with Q-vectors
  TProfile *fFourParticleCorrelationProfile[2][6];  //! 4-p correlations ([cos,sin][max number of combinations of 2 elements out of 4]) with Q-vectors
  TProfile *fSixParticleCorrelationProfile[2][4];  //! 6-p correlations ([cos,sin][max number of combinations of 3 elements out of 4]) with Q-vectors
  TProfile *fEightParticleCorrelationProfile[2];  //! 8-p correlations ([cos,sin]) with Q-vectors
  TProfile *fTwoCosineAverageProfile[2][4];  //! [<<cos(m(phi1-phi2))><cos(n(phi1-phi2))>>,<<sin(m(phi1-phi2))><sin(n(phi1-phi2))>>][6] with Q-vectors
  TProfile *fTwoNestedCorrelationProfile[2][4];  //! 2-p correlations ([cos,sin][m,n,p,q]) with nested loops
  TProfile *fFourNestedCorrelationProfile[2][6]; //! 4-p correlations ([cos,sin][max number of combinations of 2 elements out of 4]) with nested loops
  TProfile *fTwoCosineAverageNestedProfile[2][4];  //! [<<cos(m(phi1-phi2))><cos(n(phi1-phi2))>>,<<sin(m(phi1-phi2))><sin(n(phi1-phi2))>>][4] with nested loops

// Version counter for the submissions on Grid
/// Increase the counter by one when the latest version changes the structure
/// of the output file (version of the 2018-10-09, counter: 4)
  ClassDef(AliAnalysisTaskTwoMultiCorrelations,4);

}; // End of the definition of the class AliAnalysisTaskTwoMultiCorrelations

#endif

