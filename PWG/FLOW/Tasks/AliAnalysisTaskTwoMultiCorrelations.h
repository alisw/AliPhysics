
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
  void SetQAListBeforeSelection(TList* const slbs) {this->fQAListBeforeSelection = slbs;};
  TList* GetQAListBeforeSelection() const {return this->fQAListBeforeSelection;}
  void SetQAListAfterSelection(TList* const slas) {this->fQAListAfterSelection = slas;};
  TList* GetQAListAfterSelection() const {return this->fQAListAfterSelection;}
  void SetListCorrelations(TList* const slc) {this->fListCorrelations = slc;};
  TList* GetListCorrelations() const {return this->fListCorrelations;}

  void SetGeneralParameters(Int_t maxParticlesInCorrelations, Int_t maxFlowHarmonic, Bool_t useNonUnitWeights, Bool_t crossCheckFourParticle, Bool_t crossCheckNestedLoops)
  {
    this->fMaxNumberOfParticlesInCorrelations = maxParticlesInCorrelations;
    this->fHighestFlowHarmonic = maxFlowHarmonic;
    this->fUseParticleWeights = useNonUnitWeights;
    this->fCrossCheckFourParticleCorrelations = crossCheckFourParticle;
    this->fCrossCheckWithNestedLoops = crossCheckNestedLoops;
  }

  void SetAnalysisType(Bool_t aodFiles, Bool_t mcFiles, Bool_t bothFiles)
  {
    this->fProcessOnlyAOD = aodFiles;
    this->fProcessOnlyMC = mcFiles;
    this->fProcessBothMCandAOD = bothFiles;
  }

  void SetCentralityEstimation(Bool_t useVZero, Bool_t useSPD, Int_t const nBins, Double_t minCentrality, Double_t maxCentrality)
  {
    this->fCentralityFromVZero = useVZero;
    this->fCentralityFromSPD = useSPD;
    this->fCentralityMin = minCentrality;
    this->fCentralityMax = maxCentrality;
  }

  void SetVertexSelection(Bool_t cutOnVertexX, Double_t minVertexX, Double_t maxVertexX, Bool_t cutOnVertexY, Double_t minVertexY, Double_t maxVertexY, Bool_t cutOnVertexZ, Double_t minVertexZ, Double_t maxVertexZ)
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
  }

  void SetEventTracksSelection(Int_t minNumberOfTracks, Bool_t cutMaxNumberOfTracks, Int_t maxTracksZero, Int_t maxTracksFive, Int_t maxTracksTen, Int_t maxTracksTwenty, Int_t maxTracksThirty, Int_t maxTracksForty, Int_t maxTracksFifty, Int_t maxTracksSixty, Int_t maxTracksSeventy)
  {
    this->fNumberOfTracksMin = minNumberOfTracks;
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

/// Setters of the differents cuts of the track selection.
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

  void SetTPCSelection(Int_t trackFilter, Bool_t cutOnTPC, Int_t minNumberOfClustersTPC, Bool_t cutOnChiSquareTPC, Double_t minChiSquareTPC, Double_t maxChiSquareTPC)
  {
    this->fFilter = trackFilter;
    this->fCutOnNumberOfTPC = cutOnTPC;
    this->fNumberOfTPCMin = minNumberOfClustersTPC;
    this->fCutOnChiSquarePInTPC = cutOnChiSquareTPC;
    this->fChiSquarePInTPCMin = minChiSquareTPC;
    this->fChiSquarePInTPCMax = maxChiSquareTPC;
  }

  void SetDCASelection(Bool_t cutOnDCA, Double_t maxDCAxy, Double_t maxDCAz)
  {
    this->fCutOnDCA = cutOnDCA;
    this->fDCAxyMax = maxDCAxy;
    this->fDCAzMax = maxDCAz;
  }

  void SetChargeSelection(Bool_t cutOnElectricCharge, Int_t electricCharge)
  {
    this->fCutOnCharge = cutOnElectricCharge;
    this->fCharge = electricCharge;
  }

  void SetQAHistoForEventSelection(Int_t HNOTbins, Double_t HNOTmax)
  {
    this->fHNOTNumberOfBins = HNOTbins;
    this->fHNOTMax = HNOTmax;
  }

//--------------------------------------------------------------------------------------//
// Methods called in the constructors.
  virtual void InitialiseArraysOfQvectors();

//--------------------------------------------------------------------------------------//
// Methods called in 'UserCreateOutputObjects'.
  virtual void BookAllLists();
  virtual void BookQAListBeforeSelection();
  virtual void BookQAListAfterSelection();
  virtual void BookListCorrelations();

//--------------------------------------------------------------------------------------//
// Methods called in 'UserExec'.
  virtual void AnalyseAODevent(AliAODEvent *aAODevent);
  virtual void AnalyseMCevent(AliMCEvent *aMCevent);
  Bool_t ApplyTrackSelection(Double_t momentum, Double_t pseudorapidity, Int_t NclustersInTPC, Double_t TPCchiSquare, Double_t xyDCA, Double_t zDCA, Int_t eCharge);
  virtual void CalculateQvectors(long long numberOfParticles, Double_t angles[], Double_t pWeights[]);
  TComplex Q(Int_t n, Int_t p);
  virtual void ComputeMultiparticleCorrelations(long long numberOfParticles, Double_t angles[], Double_t pWeights[]);
  TComplex CalculateRecursion(Int_t n, Int_t *harmonic, Int_t mult=1, Int_t skip=0);
  virtual void ComputeTwoNestedLoops(long long nParticles, Int_t *harmonic, Double_t aAngles[], Double_t weights[], TProfile *profile, Double_t middleBin);
  virtual void ComputeFourNestedLoops(long long nParticles, Int_t *harmonic, Double_t aAngles[], Double_t weights[], TProfile *profile, Double_t middleBin);

//--------------------------------------------------------------------------------------//
// Methods called in 'Terminate'.


//======================================================================================//
// Data members.
//--------------------------------------------------------------------------------------//
private:
  AliAnalysisTaskTwoMultiCorrelations(const AliAnalysisTaskTwoMultiCorrelations& aattmc);
  AliAnalysisTaskTwoMultiCorrelations& operator=(const AliAnalysisTaskTwoMultiCorrelations& aattmc);

// Structure of the output file.
  TList *fMainList; // Mother list inside the output file.
  TList *fQAListBeforeSelection; // Daughter list with the observables before selection (event or track).
  TList *fQAListAfterSelection; // Daughter list with the observables after the full selection.
  TList *fListCorrelations; // Daughter list with the multiparticle correlations.

// General parameters.
  Int_t fMaxNumberOfParticlesInCorrelations;  // Maximum number of particles in the correlations. (default: 8)
  Int_t fHighestFlowHarmonic; // Highest flow harmonic taken into account. (default: v_6).
  TComplex fQvectors[49][9];  // All needed combinations of Q-vectors. (size: [fHighestFlowHarmonic*fMaxNumberCorrelations+1][fMaxNumberCorrelations+1])
  Bool_t fUseParticleWeights; // Use non-unit particle weights. (default: kFALSE)
  Bool_t fCrossCheckFourParticleCorrelations; // Compute the doubled combinations of harmonics for <4> to cross-check the results. (default: kFALSE)
  Bool_t fCrossCheckWithNestedLoops;  // Compute the nested loops for the 2- and 4-particle correlations to cross-check the results. (default: kFALSE)

// Type of files used in the analysis.
  Bool_t fProcessOnlyAOD; // Process only AOD files (or Reco). (default: kFALSE)
  Bool_t fProcessOnlyMC;  // Process only MC files (or Kine). (default: kFALSE)
  Bool_t fProcessBothMCandAOD;  // Process both MC and AOD files. (default: kFALSE)

// Determination of the centrality.
  Bool_t fCentralityFromVZero;  // Use the V0 detector to estimate the centrality of the events. (default: kFALSE)
  Bool_t fCentralityFromSPD;  // Use the SPD detector to estimate the centrality of the events. (default: kFALSE)
  Double_t fCentralityMin;  // Minimum value for the centrality percentile. (default: 0)
  Double_t fCentralityMax;  // Maximum value for the centrality percentile. (default: 100)

// Event selection.
  Bool_t fCutOnVertexX; // Apply the cuts on the x-position of the PV? (default: kFALSE)
  Double_t fVertexMinX; // Minimum of the x-position of the PV. (default: -44)
  Double_t fVertexMaxX; // Maximum of the x-position of the PV. (default: -44)

  Bool_t fCutOnVertexY; // Apply the cuts on the y-position of the PV? (default: kFALSE)
  Double_t fVertexMinY; // Minimum of the y-position of the PV. (default: -44)
  Double_t fVertexMaxY; // Maximum of the y-position of the PV. (default: -44)

  Bool_t fCutOnVertexZ; // Apply the cuts on the y-position of the PV? (default: kFALSE)
  Double_t fVertexMinZ; // Minimum of the z-position of the PV. (default: -10 cm)
  Double_t fVertexMaxZ; // Maximum of the z-position of the PV. (default: 10 cm)

  Int_t fNumberOfTracksMin; // Strict minimum number of tracks needed in an event to have an event weight which makes sense. (default: 6)
  Bool_t fCutOnTracksMax; // Apply maximum limits on the multiplicity to remove high multiplicity outliers? (default: kFALSE)
  Int_t fNumberOfTracksMaxZero; // Strict maximum number of tracks needed in one event for 0-5% centrality. (default: 0)
  Int_t fNumberOfTracksMaxFive; // Strict maximum number of tracks needed in one event for 5-10% centrality. (default: 0)
  Int_t fNumberOfTracksMaxTen;  // Strict maximum number of tracks needed in one event for 10-20% centrality. (default: 0)
  Int_t fNumberOfTracksMaxTwenty; // Strict maximum number of tracks needed in one event for 20-30% centrality. (default: 0)
  Int_t fNumberOfTracksMaxThirty; // Strict maximum number of tracks needed in one event for 30-40% centrality. (default: 0)
  Int_t fNumberOfTracksMaxForty;  // Strict maximum number of tracks needed in one event for 40-50% centrality. (default: 0)
  Int_t fNumberOfTracksMaxFifty;  // Strict maximum number of tracks needed in one event for 50-60% centrality. (default: 0)
  Int_t fNumberOfTracksMaxSixty;  // Strict maximum number of tracks needed in one event for 60-70% centrality. (default: 0)
  Int_t fNumberOfTracksMaxSeventy;  // Strict maximum number of tracks needed in one event for 70-80% centrality. (default: 0)

// Track selection.
  Bool_t fCutOnPt;  // Apply the cuts on the transverse momentum? (default: kFALSE)
  Double_t fPtMin;  // Minimum value of the transverse momentum. (default: 0.2 GeV)
  Double_t fPtMax;  // Maximum value of the transverse momentum. (default: 5 GeV)

  Bool_t fCutOnEta; // Apply the cuts on the pseudorapidity? (default: kFALSE)
  Double_t fEtaMin; // Minimum value of the pseudorapidity. (default: -0.8)
  Double_t fEtaMax; // Maximum value of the pseudorapidity. (default: 0.8)

  Int_t fFilter;  // Filter bit used on the tracks. (default: 128)
  Bool_t fCutOnNumberOfTPC; // Apply the cut on the number of TPC clusters? (default: kFALSE)
  Int_t fNumberOfTPCMin; // Minimum number of TPC clusters. (default: 70)
  Bool_t fCutOnChiSquarePInTPC; // Apply the cuts on chi^2 of the track momentum in TPC? (default: kFALSE)
  Double_t fChiSquarePInTPCMin;  // Minimum value of chi^2 of the track momentum in TPC. (default: 0.1)
  Double_t fChiSquarePInTPCMax;  // Maximum value of chi^2 of the track momentum in TPC. (default: 4.)

  Bool_t fCutOnDCA; // Apply the cuts on the DCA coordinates? (default: kFALSE)  
  Double_t fDCAxyMax;  // Maximum value for the xy-coordinate of the DCA. (default: 3.2 cm)
  Double_t fDCAzMax;  // Maximum value for the z-coordinate of the DCA. (default: 2.4 cm)

  Bool_t fCutOnCharge;  // Apply the cuts on the electric charge of the tracks? (default: kFALSE)
  Int_t fCharge;  // Charge of the tracks. (default: 0)

// TH1D with the observables for the event selection.
  TH1D *fHistoCentrality; //! Distribution of the centrality of the events.
  TH1I *fHistoInitialNumberOfTracks;  //! Distribution of the initial number of tracks.
  TH1I *fHistoNumberOfTracksBeforeTrackSelection; //! Distribution of the number of tracks before the track selection.
  TH1I *fHistoFinalNumberOfTracks;  //! Final number of tracks.
  Int_t fHNOTNumberOfBins;  // Number of bins for 'fHisto*NumberOfTracks'. (default: 30000)
  Double_t fHNOTMax;  // Maximum value for 'fHisto*NumberOfTracks'. (default: 30000)
  TH1D *fHistoVertexXBeforeSelection; //! Distribution of the initial PV x-position.
  TH1D *fHistoVertexXAfterSelection;  //! Distribution of the PV x-position after the full selection.
  TH1D *fHistoVertexYBeforeSelection; //! Distribution of the initial PV y-position.
  TH1D *fHistoVertexYAfterSelection;  //! Distribution of the PV y-position after the full selection.
  TH1D *fHistoVertexZBeforeSelection; //! Distribution of the initial PV z-position.
  TH1D *fHistoVertexZAfterSelection;  //! Distribution of the PV z-position after the full selection.

// TH1D with the observables for the track selection.
  TH1D *fHistoPtBeforeSelection;  //! Distribution of the transverse momentum before the track selection.
  TH1D *fHistoPtAfterSelection; //! Distribution of the transverse momentum after the full selection.
  TH1D *fHistoEtaBeforeSelection; //! Distribution of the pseudorapidity before the track selection.
  TH1D *fHistoEtaAfterSelection;  //! Distribution of the pseudorapidity after the full selection.
  TH1D *fHistoPhiBeforeSelection; //! Distribution of the azimuthal angles before the track selection.
  TH1D *fHistoPhiAfterSelection;  //! Distribution of the azimuthal angles after the full selection.
  TH1I *fHistoTPCClustersBeforeSelection; //! Distribution of the number of TPC clusters before the track selection.
  TH1I *fHistoTPCClustersAfterSelection;  //! Distribution of the number of TPC clusters after the full selection.
  TH1D *fHistoTPCChiSquareBeforeSelection;  //! Distribution of the chi square of the track momentum in the TPC before the track selection.
  TH1D *fHistoTPCChiSquareAfterSelection; //! Distribution of the chi square of the track momentum in the TPC after the full selection.
  TH1D *fHistoDCAXYBeforeSelection; //! Distribution of the xy-plane of the DCA before the track selection.
  TH1D *fHistoDCAXYAfterSelection;  //! Distribution of the xy-plane of the DCA after the full selection.
  TH1D *fHistoDCAZBeforeSelection;  //! Distribution of the z-coordinate of the DCA before the track selection.
  TH1D *fHistoDCAZAfterSelection; //! Distribution of the z-coordinate of the DCA after the full selection.
  TH1I *fHistoChargeBeforeSelection;  //! Distribution of the electric charge of the tracks before the track selection.
  TH1I *fHistoChargeAfterSelection; //! Distribution of the electric charge of the tracks after the full selection.

// TProfiles with the final multiparticle correlations.
  TProfile *fProfileTwoParticleCorrelations;  //! <2>_{j,-j} for j = 1..6. (6 bins)
  TProfile *fProfileFourParticleCorrelations; //! <4>_{j,k,-j,-k} for j = 1..6, k = j..6. (21 bins)
  TProfile *fProfileFourParticleCorrelationsCrossCheck; //! <4>_{j,k,-j,-k} for j = 2..6, k = 1..j-1. (15 bins, cross-check purpose)
  TProfile *fProfileSixParticleCorrelations;  //! <6>_{j,k,l,-j,-k,-l} for j = 1..4, k = 2..5 (k > j), l = 3..6 (l > k). (20 bins)
  TProfile *fProfileTwoParticleCorrelationsNestedLoops; //! <2>_{j,-j} for j = 1..6 with nested loops. (6 bins)
  TProfile *fProfileFourParticleCorrelationsNestedLoops;  //! <4>_{j,k,-j,-k} for j = 1..6, k = j..6 with nested loops. (21 bins)

//--------------------------------------------------------------------------------------//
// Version number to handle properly objects written before and after the changes.
// Version 9, date: 2019-04-15.
  ClassDef(AliAnalysisTaskTwoMultiCorrelations, 9);

};  // End: class AliAnalysisTaskTwoMultiCorrelations.

#endif
