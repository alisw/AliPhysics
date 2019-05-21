
//--------------------------------------------------------------------------------------// 
// Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.               //
// See cxx source for full Copyright notice                                             //
// $Id$                                                                                 //
//--------------------------------------------------------------------------------------//

//--------------------------------------------------------------------------------------//
// Analysis task for the computation of the multiparticle correlations for the flow     //
// harmonics v_1 to v_6. This version of the script compute the 2-, 4- and 6- particle  //
// correlations for all the useful combinations of these six harmonics. It can take     //
// Monte Carlo simulations data (e.g. HIJING), as well as the experimental Pb-Pb data   //
// taken by the ALICE experiment.                                                       //
// The method used to compute the multiparticle correlations is the Generic Framework   //
// based on Q-vectors. A setter lets open the possibility to cross-check the results    //
// with nested loops.                                                                   //
//                                                                                      //
// Author: Cindy Mordasini (cindy.mordasini@cern.ch)                                    //
// Version: 27.02.2019                                                                  //
//--------------------------------------------------------------------------------------//

#ifndef ALIANALYSISTASKTWOMULTICORRELATIONS_H
#define ALIANALYSISTASKTWOMULTICORRELATIONS_H

#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliVEvent.h"
#include "AliAODTrack.h"
#include "TList.h"
#include "TComplex.h"
#include "TH2I.h"
#include "TH1D.h"
#include "TH1I.h"
#include "TProfile.h"

//######################################################################################//
// Definition of the class.
//======================================================================================//
class AliAnalysisTaskTwoMultiCorrelations : public AliAnalysisTaskSE
{
public:
/* These six functions are mandatory for the class to work properly. */
  AliAnalysisTaskTwoMultiCorrelations();
  AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights=kFALSE);
  virtual ~AliAnalysisTaskTwoMultiCorrelations();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);

//--------------------------------------------------------------------------------------//
// Setters and getters for the data members.
  void SetMultiplicityList(TList* const sml) {this->fMultiplicityList = sml;};
  TList* GetMultiplcityList() const {return this->fMultiplicityList;}
  void SetEventSelectionList(TList* const sesl) {this->fEventSelectionList = sesl;};
  TList* GetEventSelectionList() const {return this->fEventSelectionList;}
  void SetTrackSelectionList(TList* const stsl) {this->fTrackSelectionList = stsl;};
  TList* GetTrackSelectionList() const {return this->fTrackSelectionList;}
  void SetMultiParticleCorrelationsList(TList* const smpcl) {this->fMultiParticleCorrelationsList = smpcl;};
  TList* GetMultiParticleCorrelationsList() const {return this->fMultiParticleCorrelationsList;}
  void SetTwoParticleCorrelationsWithEtaGapsList(TList* const stpcl) {this->fTwoParticleCorrelationsWithEtaGapsList = stpcl;};
  TList* GetTwoParticleCorrelationsWithEtaGapsList() const {return this->fTwoParticleCorrelationsWithEtaGapsList;}

/// Configure the general parameters of the analysis.
  void SetAnalysisDimensions(Int_t maxFlowHarmonic, Int_t maxCorrelatedParticles)
  {
    this->fHighestFlowHarmonic = maxFlowHarmonic;
    this->fMaxNumberOfParticlesInCorrelations = maxCorrelatedParticles;
  }

  void SetTypeOfFiles(Bool_t aodFiles, Bool_t mcFiles, Bool_t bothFiles)
  {
    this->fProcessOnlyAOD = aodFiles;
    this->fProcessOnlyMC = mcFiles;
    this->fProcessBothMCandAOD = bothFiles;
  }

  void SetGeneralParameters( Bool_t computeEtaGaps, Bool_t crosscheckWithNestedLoops, Bool_t doCorrelationPlot, Bool_t useNonUnitWeights)
  {
    this->fComputeEtaGaps = computeEtaGaps;
    this->fCrosscheckWithNestedLoops = crosscheckWithNestedLoops;
    this->fDoTDCorrelationHisto = doCorrelationPlot;
    this->fUseParticleWeights = useNonUnitWeights;
  }

/// Configure the parameters related to the number of tracks.
  void SetCentralityDetermination(Bool_t useVZero, Bool_t useSPD, Int_t const nBins, Double_t minCentrality, Double_t maxCentrality)
  {
    this->fCentralityFromVZero = useVZero;
    this->fCentralityFromSPD = useSPD;
    this->fCentralityMin = minCentrality;
    this->fCentralityMax = maxCentrality;
  }

  void SetHistoParameters(Int_t HNOTbins, Double_t HNOTmax, Int_t HFCbins, Double_t HFCmax)
  {
    this->fNumberOfBinsHNOT = HNOTbins;
    this->fMaxBinHNOT = HNOTmax;
    this->fNumberOfBinsHFC = HFCbins;
    this->fMaxBinHFC = HFCmax;
  }

/// Configure the parameters related to the event selection criteria.
  void SetPVSelection(Bool_t cutOnVertexX, Double_t minVertexX, Double_t maxVertexX, Bool_t cutOnVertexY, Double_t minVertexY, Double_t maxVertexY, Bool_t cutOnVertexZ, Double_t minVertexZ, Double_t maxVertexZ)
  {
    this->fCutOnPVX = cutOnVertexX;
    this->fPVXMin = minVertexX;
    this->fPVXMax = maxVertexX;
    this->fCutOnPVY = cutOnVertexY;
    this->fPVYMin = minVertexY;
    this->fPVYMax = maxVertexY;
    this->fCutOnPVZ = cutOnVertexZ;
    this->fPVZMin = minVertexZ;
    this->fPVZMax = maxVertexZ;
  }

// Configure the parameters related to the track selection criteria.
  void SetFiltersSelection(Int_t firstFilter, Int_t globalFilter)
  {
    this->fMainFilter = firstFilter;
    this->fGlobalFilter = globalFilter;
  }

  void SetCorrelationSelection(Int_t minNumberOfTracks, Bool_t removeHMoutliers, Double_t minA, Double_t minB, Double_t maxA, Double_t maxB)
  {
    this->fMultiplicityMin = minNumberOfTracks;
    this->fCutOnTDCorrelations = removeHMoutliers;
    this->fMultiplicityMinA = minA;
    this->fMultiplicityMinB = minB;
    this->fMultiplicityMaxA = maxA;
    this->fMultiplicityMaxB = maxB;
  }

  void SetPtSelection(Bool_t cutOnPt, Double_t minPt, Double_t maxPt)
  {
    this->fCutOnPt = cutOnPt;
    this->fPtMin = minPt;
    this->fPtMax = maxPt;
  }

  void SetEtaSelection(Bool_t cutOnEta, Double_t minEta, Double_t maxEta)
  {
    this->fCutOnEta = cutOnEta;
    this->fEtaMin = minEta;
    this->fEtaMax = maxEta;
  }

  void SetTPCSelection(Bool_t cutOnTPC, Int_t minNumberOfClustersTPC, Bool_t cutOnChiSquareTPC, Double_t minChiSquareTPC, Double_t maxChiSquareTPC)
  {
    this->fCutOnNumberOfTPC = cutOnTPC;
    this->fNumberOfTPCMin = minNumberOfClustersTPC;
    this->fCutOnChiSquarePInTPC = cutOnChiSquareTPC;
    this->fChiSquarePInTPCMin = minChiSquareTPC;
    this->fChiSquarePInTPCMax = maxChiSquareTPC;
  }

  void SetDCASelection(Bool_t cutOnDCAxy, Double_t maxDCAxy, Bool_t cutOnDCAz, Double_t maxDCAz)
  {
    this->fCutOnDCAxy = cutOnDCAxy;
    this->fDCAxyMax = maxDCAxy;
    this->fCutOnDCAz = cutOnDCAz;
    this->fDCAzMax = maxDCAz;
  }

  void SetChargeSelection(Bool_t cutOnElectricCharge, Int_t electricCharge)
  {
    this->fCutOnCharge = cutOnElectricCharge;
    this->fCharge = electricCharge;
  }

  void SetITSSelection(Bool_t cutOnITS, Int_t minNumberOfClustersITS)
  {
    this->fCutOnNumberOfITS = cutOnITS;
    this->fNumberOfITSMin = minNumberOfClustersITS;
  }

  void SetEventTracksSelection(Bool_t cutMaxNumberOfTracks, Int_t maxTracksZero, Int_t maxTracksFive, Int_t maxTracksTen, Int_t maxTracksTwenty, Int_t maxTracksThirty, Int_t maxTracksForty, Int_t maxTracksFifty, Int_t maxTracksSixty, Int_t maxTracksSeventy)
  {
    this->fCutOnTracksMax = cutMaxNumberOfTracks;
    this->fNumberOfTracksMaxZero = maxTracksZero;
    this->fNumberOfTracksMaxFive = maxTracksFive;
    this->fNumberOfTracksMaxTen = maxTracksTen;
    this->fNumberOfTracksMaxTwenty = maxTracksTwenty;
    this->fNumberOfTracksMaxThirty = maxTracksThirty;
    this->fNumberOfTracksMaxForty = maxTracksForty;
    this->fNumberOfTracksMaxFifty = maxTracksFifty;
    this->fNumberOfTracksMaxSixty = maxTracksSixty;
    this->fNumberOfTracksMaxSeventy = maxTracksSeventy;
  } 

// Configure the parameters related to the multi-particle correlations.
// Configure the parameters related to the 2-particle correlations with eta gaps.

//--------------------------------------------------------------------------------------//
// Methods called in the constructors.
  virtual void InitialiseArraysOfQvectors();
  virtual void InitialiseArraysOfTProfiles();

//--------------------------------------------------------------------------------------//
// Methods called in 'UserCreateOutputObjects'.
  virtual void BookAllLists();
  virtual void BookMultiplicityList();
  virtual void BookEventSelectionList();
  virtual void BookTrackSelectionList();
  virtual void BookMultiParticleCorrelationsList();
  virtual void BookTwoParticleCorrelationsWithEtaGapsList();

//--------------------------------------------------------------------------------------//
// Methods called in 'UserExec'.
  virtual void AnalyseMCevent(AliMCEvent *aMCevent);
  virtual void AnalyseAODevent(AliAODEvent *aAODevent);
  Bool_t ApplyTrackSelection(AliAODTrack *aAODtrack);
  virtual void CalculateQvectors(long long numberOfParticles, Double_t angles[], Double_t pWeights[]);
  TComplex Q(Int_t n, Int_t p);
  virtual void ComputeMultiparticleCorrelations(long long numberOfParticles, Double_t angles[], Double_t pWeights[]);
  TComplex CalculateRecursion(Int_t n, Int_t *harmonic, Int_t mult=1, Int_t skip=0);
  virtual void ComputeTwoNestedLoops(long long nParticles, Int_t *harmonic, Double_t aAngles[], Double_t weights[], TProfile *profile, Double_t middleBin);
  virtual void ComputeFourNestedLoops(long long nParticles, Int_t *harmonic, Double_t aAngles[], Double_t weights[], TProfile *profile, Double_t middleBin);
  virtual void ComputeTwoParticleEtaGaps(long long nParticles, Double_t angles[], Double_t pWeights[], Double_t pseudorapidity[]);

//--------------------------------------------------------------------------------------//
// Methods called in 'Terminate'.


//======================================================================================//
// Data members.
//--------------------------------------------------------------------------------------//
private:
  AliAnalysisTaskTwoMultiCorrelations(const AliAnalysisTaskTwoMultiCorrelations& aattmc);
  AliAnalysisTaskTwoMultiCorrelations& operator=(const AliAnalysisTaskTwoMultiCorrelations& aattmc);

// General parameters of the analysis.
  TList *fMainList; // Mother list in the output file.

  Int_t fHighestFlowHarmonic; // Highest flow harmonic taken into account.
  Int_t fMaxNumberOfParticlesInCorrelations;  // Maximum number of particles used in the correlators.
  TComplex fQvectors[49][9];  // All needed combinations of Q-vectors. (size: [fHighestFlowHarmonic*fMaxNumberCorrelations+1][fMaxNumberCorrelations+1])

  Bool_t fProcessOnlyAOD; // Process only AOD files (or Reco)?
  Bool_t fProcessOnlyMC;  // Process only MC files (or Kine)?
  Bool_t fProcessBothMCandAOD;  // Process both MC and AOD files?

  Bool_t fComputeEtaGaps; // Compute the eta gaps method for the 2-particle correlations?
  Bool_t fCrosscheckWithNestedLoops;  // Crosscheck the results with the nested loops for the 2- and 4-particle correlations?
  Bool_t fDoTDCorrelationHisto; // Fill of the 2D correlation histogram?
  Bool_t fUseParticleWeights; // Use non-unit particle weights? [TBA: full gestion of the non-unit weights in the task.]

// Parameters related to the number of tracks.
  TList *fMultiplicityList; // Daughter list for the histograms containing the number of tracks.
  TH1D *fHistoCentrality; //! Distribution of the centrality of the events.
  TH1I *fHistoInitialNumberOfTracks;  //! Distribution of the initial number of tracks.
  TH1I *fHistoIntermediateNumberOfTracks; //! Distribution of the number of tracks before the track selection.
  TH1I *fHistoFinalNumberOfTracks;  //! Distribution of the final number of tracks.
  TH2I *fHistoFilterCorrelations; //! 2D correlation histogram with multiplicity for two filter bits before the track selection.
  TH2I *fHistoFinalFilterCorrelations; //! 2D correlation histogram with multiplicity for two filter bits after the track selection.

  Bool_t fCentralityFromVZero;  // Use the V0 detector to estimate the centrality of the events?
  Bool_t fCentralityFromSPD;  // Use the SPD detector to estimate the centrality of the events? 
  Double_t fCentralityMin;  // Minimum value of the centrality percentile.
  Double_t fCentralityMax;  // Maximum value of the centrality percentile.

  Int_t fNumberOfBinsHNOT;  // Number of bins for the TH1I with the number of tracks.
  Double_t fMaxBinHNOT; // Last bin for the TH1I with the number of tracks.
  Int_t fNumberOfBinsHFC; // Number of bins for the TH2I with the correlations between the filters.
  Double_t fMaxBinHFC;  // Last bin for the TH2I with the correlations between the filters.

// Parameters related to the event selection criteria.
  TList *fEventSelectionList; // Daughter list for the histograms containing the event selection criteria.
  TH1D *fHistoInitialPVX; //! Initial distribution of the x-position of the PV.
  TH1D *fHistoFinalPVX; //! Final distribution of the x-position of the PV.
  TH1D *fHistoInitialPVY; //! Initial distribution of the y-position of the PV.
  TH1D *fHistoFinalPVY; //! Final distribution of the y-position of the PV.
  TH1D *fHistoInitialPVZ; //! Initial distribution of the z-position of the PV.
  TH1D *fHistoFinalPVZ; //! Final distribution of the z-position of the PV.

  Bool_t fCutOnPVX; // Apply the cuts on the x-position of the PV?
  Bool_t fCutOnPVY; // Apply the cuts on the y-position of the PV?
  Bool_t fCutOnPVZ; // Apply the cuts on the z-position of the PV?

  Double_t fPVXMin; // Minimum cut on the x-position of the PV.
  Double_t fPVXMax; // Maximum cut on the x-position of the PV.
  Double_t fPVYMin; // Minimum cut on the y-position of the PV.
  Double_t fPVYMax; // Maximum cut on the y-position of the PV.
  Double_t fPVZMin; // Minimum cut on the z-position of the PV.
  Double_t fPVZMax; // Maximum cut on the z-position of the PV.

// Parameters related to the track selection criteria.
  TList *fTrackSelectionList; // Daughter list for the histograms containing the track selection criteria.
  TH1D *fHistoIntermediatePt; //! Distribution of the transverse momentum before the track selection.
  TH1D *fHistoFinalPt;  //! Distribution of the transverse momentum after the full selection.
  TH1D *fHistoIntermediateEta;  //! Distribution of the pseudorapidity before the track selection.
  TH1D *fHistoFinalEta; //! Distribution of the pseudorapidity after the full selection.
  TH1D *fHistoIntermediatePhi;  //! Distribution of the azimuthal angles before the track selection.
  TH1D *fHistoFinalPhi; //! Distribution of the azimuthal angles after the full selection.
  TH1I *fHistoIntermediateNumberOfTPC;  //! Distribution of the number of TPC clusters before the track selection.
  TH1I *fHistoFinalNumberOfTPC; //! Distribution of the number of TPC clusters after the full selection.
  TH1D *fHistoIntermediateChiSquare;  //! Distribution of the chi^2 of the track momentum in TPC before the track selection.
  TH1D *fHistoFinalChiSquare; //! Distribution of the chi^2 of the track momentum in TPC after the full selection.
  TH1D *fHistoIntermediateDCAxy;  //! Distribution of the xy-coordinate of the DCA before the track selection.
  TH1D *fHistoFinalDCAxy; //! Distribution of the xy-coordinate of the DCA after the full selection.
  TH1D *fHistoIntermediateDCAz; //! Distribution of the z-coordinate of the DCA before the track selection.
  TH1D *fHistoFinalDCAz;  //! Distribution of the z-coordinate of the DCA after the full selection.
  TH1I *fHistoIntermediateCharge; //! Distribution of the electric charge of the tracks before the track selection.
  TH1I *fHistoFinalCharge;  //! Distribution of the electric charge of the tracks after the full selection.
  TH1I *fHistoIntermediateNumberOfITS;  //! Distribution of the number of clusters in the ITS before the track selection.
  TH1I *fHistoFinalNumberOfITS; //! Distribution of the number of clusters in the ITS after the track selection.

  Bool_t fCutOnTDCorrelations;  // Apply the cuts on the number of tracks to remove high multiplicity outliers?
  Bool_t fCutOnPt;  // Apply the cuts on the transverse momentum?
  Bool_t fCutOnEta; // Apply the cuts on the pseudorapidity?
  Bool_t fCutOnNumberOfTPC; // Apply the cut on the number of TPC clusters?
  Bool_t fCutOnChiSquarePInTPC; // Apply the cuts on chi^2 of the track momentum in TPC?
  Bool_t fCutOnDCAxy; // Apply the cut on the xy-coordinate of the DCA?
  Bool_t fCutOnDCAz;  // Apply the cut on the z-coordinate of the DCA?
  Bool_t fCutOnCharge;  // Apply the cut on the electric charge of the tracks?
  Bool_t fCutOnNumberOfITS; // Apply the cut on the number of ITS clusters?

  Int_t fMainFilter; // Main filter bit used in the analysis.
  Int_t fGlobalFilter;  // Second filter bit used to remove high multiplicity outliers.
  Int_t fMultiplicityMin;  // Strict minimum of the number of tracks required in an event to have an event weight which makes sense.
  Double_t fMultiplicityMinA; // a in 'a(global multiplicity) + b' for the minimum boundary.
  Double_t fMultiplicityMinB; // b in 'a(global multiplicity) + b' for the minimum boundary.
  Double_t fMultiplicityMaxA; // a in 'a(global multiplicity) + b' for the maximum boundary.
  Double_t fMultiplicityMaxB; // b in 'a(global multiplicity) + b' for the maximum boundary.

  Double_t fPtMin;  // Minimum cut on the transverse momentum.
  Double_t fPtMax;  // Maximum cut on the transverse momentum.
  Double_t fEtaMin; // Minimum cut on the pseudorapidity.
  Double_t fEtaMax; // Maximum cut on the pseudorapidity.
  Int_t fNumberOfTPCMin;  // Minimum number of TPC clusters.
  Double_t fChiSquarePInTPCMin; // Minimum cut on chi^2 of the track momentum in TPC.
  Double_t fChiSquarePInTPCMax; // Maximum cut on chi^2 of the track momentum in TPC.
  Double_t fDCAxyMax; // Maximum cut on the xy-coordinate of the DCA.
  Double_t fDCAzMax;  // Maximum cut on the z-coordinate of the DCA.
  Int_t fCharge;  // Electric charge of all the tracks.
  Int_t fNumberOfITSMin;  // Minimum number of ITS clusters.

/// Brute-force method to remove the HM outliers.
  Bool_t fCutOnTracksMax; // Apply maximum limits on the multiplicity to remove high multiplicity outliers?
  Int_t fNumberOfTracksMaxZero; // Strict maximum number of tracks needed in one event for 0-5% centrality.
  Int_t fNumberOfTracksMaxFive; // Strict maximum number of tracks needed in one event for 5-10% centrality.
  Int_t fNumberOfTracksMaxTen;  // Strict maximum number of tracks needed in one event for 10-20% centrality.
  Int_t fNumberOfTracksMaxTwenty; // Strict maximum number of tracks needed in one event for 20-30% centrality.
  Int_t fNumberOfTracksMaxThirty; // Strict maximum number of tracks needed in one event for 30-40% centrality.
  Int_t fNumberOfTracksMaxForty;  // Strict maximum number of tracks needed in one event for 40-50% centrality.
  Int_t fNumberOfTracksMaxFifty;  // Strict maximum number of tracks needed in one event for 50-60% centrality.
  Int_t fNumberOfTracksMaxSixty;  // Strict maximum number of tracks needed in one event for 60-70% centrality.
  Int_t fNumberOfTracksMaxSeventy;  // Strict maximum number of tracks needed in one event for 70-80% centrality.

// Parameters related to the multi-particle correlations.
  TList *fMultiParticleCorrelationsList; // Daughter list with the multi-particle correlations.
  TProfile *fProfileTwoParticleCorrelations;  //! <2>_{j,-j} for j = 1..6. (6 bins)
  TProfile *fProfileFourParticleCorrelations; //! <4>_{j,k,-j,-k} for j = 1..6, k = j..6. (21 bins)
  TProfile *fProfileFourParticleCorrelationsCrossCheck; //! <4>_{j,k,-j,-k} for j = 2..6, k = 1..j-1. (15 bins)
  TProfile *fProfileSixParticleCorrelations;  //! <6>_{j,k,l,-j,-k,-l} for j = 1..4, k = 2..5 (k > j), l = 3..6 (l > k). (20 bins)
  TProfile *fProfileTwoParticleCorrelationsNestedLoops; //! <2>_{j,-j} for j = 1..6 with nested loops. (6 bins)
  TProfile *fProfileFourParticleCorrelationsNestedLoops;  //! <4>_{j,k,-j,-k} for j = 1..6, k = j..6 with nested loops. (21 bins)

// Parameters related to the 2-particle correlations with eta gaps.
  TList *fTwoParticleCorrelationsWithEtaGapsList; // Daughter list with the 2-p correlations calculated with eta gaps.
  TProfile *fProfileTwoParticleCorrelationsWithEtaGaps[6];  //! <2>_{j,-j} for j = 1..6 with the eta gaps method. (11 bins per profile).

//--------------------------------------------------------------------------------------//
// Version number to handle properly objects written before and after the changes.
// Version 11, date: 2019-05-20.
  ClassDef(AliAnalysisTaskTwoMultiCorrelations, 11);

};  // End: class AliAnalysisTaskTwoMultiCorrelations.

#endif
