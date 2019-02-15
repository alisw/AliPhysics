
//--------------------------------------------------------------------------------------// 
// Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.               //
// See cxx source for full Copyright notice                                             //
// $Id$                                                                                 //
//--------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------//
// Analysis task for the computation of the multiparticle correlations with different   //
// flow harmonics v_n. This script can takes the Monte Carlo simulations data (e.g.     //
// HIJING), as well as the experimental Pb-Pb data taken by the ALICE experiment.       //
// The current script computes the multiparticle correlators using the method of the    //
// Q-vectors for a maximum of 6 different harmonics and 8 particles).                   //
//                                                                                      //
// Author: Cindy Mordasini (cindy.mordasini@cern.ch)                                    //
// Version: 12.02.2019                                                                  //
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
// Definition of the class.
class AliAnalysisTaskTwoMultiCorrelations : public AliAnalysisTaskSE{
public:
/* The six following functions are mandatory for AliAnalysisTaskSE to use the class properly. */
  AliAnalysisTaskTwoMultiCorrelations();
  AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights=kFALSE);
  virtual ~AliAnalysisTaskTwoMultiCorrelations();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);

//--------------------------------------------------------------------------------------//
// Setters and getters for all the data members.
  void SetGeneralParameters(Int_t maxNumberCorrelations, Int_t highestHarmonic, Int_t nbHarmonics, Bool_t useParticleWeights, Bool_t computeNestedLoops)
  {
    this->fMaxNumberCorrelations = maxNumberCorrelations;
    this->fMaxFlowHarmonic = highestHarmonic;
    this->fNumberHarmonicsInSC = nbHarmonics;
    this->fUseParticleWeights = useParticleWeights;
    this->fComputeNestedLoops = computeNestedLoops;
  } // End: void SetGeneralParameters().

  void SetControlListEventCuts(TList* const sclec) {this->fControlListEventCuts = sclec;};
  TList* GetControlListEventCuts() const {return this->fControlListEventCuts;}
  void SetControlListTrackCuts(TList* const scltc) {this->fControlListTrackCuts = scltc;};
  TList* GetControlListTrackCuts() const {return this->fControlListTrackCuts;}
  void SetListCorrelations(TList* const slc) {this->fListCorrelations = slc;};
  TList* GetListCorrelations() const {return this->fListCorrelations;}

  void SetAnalysisType(Bool_t bothAnalysis, Bool_t mcAnalysis, Bool_t aodAnalysis)
  {
    this->fProcessBothKineAndReco = bothAnalysis;
    this->fProcessOnlyKine = mcAnalysis;
    this->fProcessOnlyReco = aodAnalysis;
  } // End: SetAnalysisType().

  void SetCentralityEstimation(Bool_t useSPD, Bool_t useVzero, Int_t const nBins, Float_t minCentrality, Float_t maxCentrality)
  {
    this->fUseSPDForCentrality = useSPD;
    this->fUseVZeroForCentrality = useVzero;
    this->fCentralityMin = minCentrality;
    this->fCentralityMax = maxCentrality;
  } // End: void SetCentralityEstimation().

  void SetEventSelection(Bool_t cutOnVertexX, Float_t minVertexX, Float_t maxVertexX, Bool_t cutOnVertexY, Float_t minVertexY, Float_t maxVertexY, Bool_t cutOnVertexZ, Float_t minVertexZ, Float_t maxVertexZ)
  {
    this->fCutOnVertexX = cutOnVertexX;
    this->fVertexMinX = minVertexX;
    this->fVertexMaxX = maxVertexX;
    this->fCutOnVertexY = cutOnVertexY;
    this->fVertexMinY = minVertexY;
    this->fVertexMaxY = maxVertexY;
    this->fCutOnVertexZ = cutOnVertexZ;
    this->fVertexMinZ = minVertexZ;
    this->fVertexMaxZ = maxVertexZ;
  } // End: void SetEventSelection().

  void SetTrackSelection(Float_t minPt, Float_t maxPt, Float_t minEta, Float_t maxEta, Int_t minNumberOfClustersTPC, Float_t maxDCAxy, Float_t maxDCAz)
  {
    this->fPtMin = minPt;
    this->fPtMax = maxPt;
    this->fEtaMin = minEta;
    this->fEtaMax = maxEta;
    this->fNumberOfTPCMin = minNumberOfClustersTPC;
    this->fDCAxyMax = maxDCAxy;
    this->fDCAzMax = maxDCAz;
  } // End: void SetTrackSelection().

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
  } // End: void SetHarmonics().

// Methods called in the constructor.
  virtual void InitialiseArraysOfQvectors();
  virtual void InitialiseArraysOfTProfiles();

// Methods called in UserCreateOutputObjects().
  virtual void BookAllLists();
  virtual void BookControlListEventCuts();
  virtual void BookControlListTrackCuts();
  virtual void BookListCorrelations();

// Methods called in UserExec(Option_t *).
  virtual void AODanalysis(AliAODEvent *aAODevent);
  virtual void MCanalysis(AliMCEvent *aMCevent);
  Bool_t CreateTrackSelection(Double_t currentPt, Double_t currentEta, Int_t currentNumberOfTPC, Double_t currentDCAXY, Double_t currentDCAZ);
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
  Int_t fMaxNumberCorrelations; // Maximum number of particles in the correlator (default: 8).
  Int_t fMaxFlowHarmonic; // Maximum harmonic n for v_n (default: 6).
  TComplex fQvectors[49][9];  // All needed combinations of Q-vectors (size: [fMaxFlowHarmonic*fMaxNumberCorrelations+1][fMaxNumberCorrelations+1]).
  Int_t fNumberHarmonicsInSC; // Number of harmonics in the GSC (default: 3).
  Bool_t fUseParticleWeights; // Use non-unit particle weights (default: kFALSE).
  Bool_t fComputeNestedLoops; // Compute the nested loops for cross-check (default: kFALSE).

// Structure of the output file.
  TList *fMainList; // Main output list.
  TList *fControlListEventCuts; // Secondary list with the observables for the event cuts.
  TList *fControlListTrackCuts; // Secondary list with the observables for the track cuts.
  TList *fListCorrelations; // Secondary list with the results for all the correlations.
  TList *fListTwoParticles; // Tertiary list with the 2-p correlations.
  TList *fListThreeParticles; // Tertiary list with the 3-p correlations.
  TList *fListFourParticles;  // Tertiary list with the 4-p correlations.
  TList *fListSixParticles; // Tertiary list with the 6- and 8-p correlations.

// Control histograms for the distribution of observables for the event selection.
  TH1D *fHistoCentrality; //! Centrality distribution.
  TH1D *fHistoNumberOfTracksBeforeCuts; //! Number of tracks in the event before any selection.
  TH1D *fHistoNumberOfTracksAfterEventCuts; //! Number of tracks remaining after the event selection.
  TH1D *fHistoVertexXBeforeCuts;  //! x-position of the PV before the event selection.
  TH1D *fHistoVertexYBeforeCuts;  //! y-position of the PV before the event selection.
  TH1D *fHistoVertexZBeforeCuts;  //! z-position of the PV before the event selection.
  TH1D *fHistoVertexXAfterEventCuts;  //! x-position of the PV after the event selection.
  TH1D *fHistoVertexYAfterEventCuts;  //! y-position of the PV after the event selection.
  TH1D *fHistoVertexZAfterEventCuts;  //! z-position of the PV after the event selection.

// Control histrograms for the distribution of observables for the track selection.
  TH1D *fHistoNumberOfTracksAfterAllCuts; //! Number of tracks remaining after both the event and the track selection.
  TH1D *fHistoPtBeforeCuts; //! Transverse momentum before the track selection.
  TH1D *fHistoPtAfterCuts;  //! Transverse momentum distribution after the track selection.  
  TH1D *fHistoEtaBeforeCuts;  //! Pseudorapidity distribution before the track selection.
  TH1D *fHistoEtaAfterCuts; //! Pseudorapidity distribution after the track selection.
  TH1D *fHistoPhiBeforeCuts;  //! Azimuthal angles distribution before the track selection.
  TH1D *fHistoPhiAfterCuts; //! Azimuthal angles distribution after the track selection.
  TH1D *fHistoTPCClustersBeforeCuts;  //! Number of TPC clusters before the track selection.
  TH1D *fHistoTPCClustersAfterCuts; //! Number of TPC clusters after the track selection.
  TH1D *fHistoDCAXYBeforeCuts; //! xy-plane of the DCA before the track selection.
  TH1D *fHistoDCAZBeforeCuts; //! z-coordinate of the DCA before the track selection.
  TH1D *fHistoDCAXYAfterCuts;  //! xy-plane of the DCA after the track selection.
  TH1D *fHistoDCAZAfterCuts;  //! z-coordinate of the DCA after the track selection.

// TProfiles with the final multiparticle correlations.
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

// Type of analysis.
  TString *fAnalysisType; //! Type of analysis: MC or AOD.
  Bool_t fProcessBothKineAndReco; // Process both MC and AOD (default: kFALSE).
  Bool_t fProcessOnlyKine;  // Process only MC (default: kFALSE).
  Bool_t fProcessOnlyReco;  // Process only AOD (default: kFALSE).

// Determination of the centrality.
  TString *fCentralitySelection;  //! Detector for the centrality estimation: SPD (CL1) or V0 (V0M).
  Bool_t fUseSPDForCentrality;  // Use the SPD detector (default: kFALSE).
  Bool_t fUseVZeroForCentrality;  // Use the V0 detector (default: kFALSE).
  Double_t fCentralityMin;  // Minimum value for the centrality percentile (default: 0).
  Double_t fCentralityMax;  // Maximum value for the centrality percentile (default: 100).

// Event selection.
  Bool_t fCutOnVertexX; // Apply the cuts on the x-position of the PV (default: kFALSE).
  Double_t fVertexMinX; // Minimum of the x-position of the PV (default: -44).
  Double_t fVertexMaxX; // Maximum of the x-position of the PV (default: -44).

  Bool_t fCutOnVertexY; // Apply the cuts on the y-position of the PV (default: kFALSE).
  Double_t fVertexMinY; // Minimum of the y-position of the PV (default: -44).
  Double_t fVertexMaxY; // Maximum of the y-position of the PV (default: -44).

  Bool_t fCutOnVertexZ; // Apply the cuts on the y-position of the PV (default: kFALSE).
  Double_t fVertexMinZ; // Minimum of the z-position of the PV (default: -10 cm).
  Double_t fVertexMaxZ; // Maximum of the z-position of the PV (default: 10 cm).

// Track selection.
  Double_t fPtMin;  // Minimum value of the transverse momentum (default: 0.2 GeV).
  Double_t fPtMax;  // Maximum value of the transverse momentum (default: 5 GeV).

  Double_t fEtaMin; // Minimum value of the pseudorapidity (default: -0.8).
  Double_t fEtaMax; // Maximum value of the pseudorapidity (default: 0.8).

  Double_t fNumberOfTPCMin; // Minimum number of TPC clusters (default: 70?).

  Double_t fDCAxyMax;  // Maximum value for the xy-coordinate of the DCA (default: 3.2 cm).
  Double_t fDCAzMax;  // Maximum value for the z-coordinate of the DCA (default: 2.4 cm).

// Harmonics.
  Int_t fHarmonicOne; // Harmonic n_1 (default: 2).
  Int_t fHarmonicTwo; // Harmonic n_2 (default: -2).
  Int_t fHarmonicThree; // Harmonic n_3 (default: 3).
  Int_t fHarmonicFour;  // Harmonic n_4 (default: -3).
  Int_t fHarmonicFive;  // Harmonic n_5 (default: 4).
  Int_t fHarmonicSix; // Harmonic n_6 (default: -4).
  Int_t fHarmonicSeven; // Harmonic n_7 (default: 0).
  Int_t fHarmonicEight; // Harmonic n_8 (default: 0).

// Version counter for the submissions on Grid.
/// Increase the counter by one when the latest version changes the structure
/// of the output file (version of the 2019-02-08, counter: 6).
  ClassDef(AliAnalysisTaskTwoMultiCorrelations,6);

};  // End: class AliAnalysisTaskTwoMultiCorrelations().

#endif

