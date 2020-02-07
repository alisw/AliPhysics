
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
/ Version 17 from the 07.02.2020.                                             /
/ -------------------------------------------------------------------------- */

#ifndef ALIANALYSISTASKTWOMULTICORRELATIONS_H
#define ALIANALYSISTASKTWOMULTICORRELATIONS_H

#include "AliAnalysisTaskSE.h"
#include "TSystem.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TList.h"
#include "TComplex.h"
#include "AliJEfficiency.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliMultSelection.h"
#include "AliAODTrack.h"
#include "AliAODMCParticle.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2I.h"
#include "TProfile.h"
#include <TExMap.h>

class AliAnalysisTaskTwoMultiCorrelations : public AliAnalysisTaskSE
{
public:
/* Mandatory functions for the class to work properly within the framework. */
  AliAnalysisTaskTwoMultiCorrelations();
  AliAnalysisTaskTwoMultiCorrelations(const char *name, Bool_t useParticleWeights=kFALSE);
  virtual ~AliAnalysisTaskTwoMultiCorrelations();
  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *);
  virtual void Terminate(Option_t *);

/* Setters and getters for the data members. */
  void    SetMultiplicityList(TList* const sml) {this->fMultiplicityList = sml;};
  TList*  GetMultiplcityList() const {return this->fMultiplicityList;}
  void    SetEventSelectionList(TList* const sesl) {this->fEventSelectionList = sesl;};
  TList*  GetEventSelectionList() const {return this->fEventSelectionList;}
  void    SetTrackSelectionList(TList* const stsl) {this->fTrackSelectionList = stsl;};
  TList*  GetTrackSelectionList() const {return this->fTrackSelectionList;}
  void    SetMPCList(TList* const smpcl) {this->fMPCList = smpcl;};
  TList*  GetMPCList() const {return this->fMPCList;}
  void    SetTPCEtaList(TList* const stpcl) {this->fTPCEtaList = stpcl;};
  TList*  GetTPCEtaList() const {return this->fTPCEtaList;}

// 1. Configuration of the general parameters.
  void SetAnalysisParameters(Int_t maxHarmo, Int_t maxCorrel,
    Bool_t doKine, Bool_t drawTriangles, Bool_t doEtaGaps)
  {
    this->fHighestHarmonic        = maxHarmo;
    this->fLargestCorrelators     = maxCorrel;
    this->fDoKineAnalysis         = doKine;
    this->fDoCorrelationsHisto    = drawTriangles;
    this->fComputeEtaGaps         = doEtaGaps;
  } // End: void SetAnalysisParameters(Int_t, Int_t, Bool_t, Bool_t, Bool_t).

// 2. Configuration of the centrality parameters.
  void SetCentrality(Int_t const nBins, Float_t minCentrality, Float_t maxCentrality,
    Bool_t useSPD, Bool_t useVZero)
  {
    this->fCentralityMin          = minCentrality;
    this->fCentralityMax          = maxCentrality;
    this->fCentralityFromSPD      = useSPD;
    this->fCentralityFromVZero    = useVZero;
  } // End: void SetCentrality(Float_t, Float_t, Bool_t, Bool_t).

// 4. Configuration of the physics event selection.
  void SetPVxSelection(Bool_t cutPVx, Float_t minPVx, Float_t maxPVx)
  {
    this->fCutOnPVx   = cutPVx;
    this->fPVxMin     = minPVx;
    this->fPVxMax     = maxPVx;
  } // End: void SetPVxSelection(Bool_t, Float_t, Float_t).

  void SetPVySelection(Bool_t cutPVy, Float_t minPVy, Float_t maxPVy)
  {
    this->fCutOnPVy   = cutPVy;
    this->fPVyMin     = minPVy;
    this->fPVyMax     = maxPVy;
  } // End: void SetPVySelection(Bool_t, Float_t, Float_t).

  void SetPVzSelection(Bool_t cutPVz, Float_t minPVz, Float_t maxPVz)
  {
    this->fCutOnPVz   = cutPVz;
    this->fPVzMin     = minPVz;
    this->fPVzMax     = maxPVz;
  } // End: void SetPVzSelection(Bool_t, Float_t, Float_t).

// 5. Configuration of the HMOs selection.
  void SetHMOsSelection(Int_t minMultiplicity, Int_t mainFilter, Int_t secondFilter,
    Bool_t cutHMOs, Float_t minA, Float_t minB, Float_t maxA, Float_t maxB)
  {
    this->fMultiplicityMin    = minMultiplicity;
    this->fMainFilter         = mainFilter;
    this->fGlobalFilter       = secondFilter;
    this->fCutOnHMOs          = cutHMOs;
    this->fMultiplicityMinA   = minA;
    this->fMultiplicityMinB   = minB;
    this->fMultiplicityMaxA   = maxA;
    this->fMultiplicityMaxB   = maxB;
  } // End: void SetHMOsSelection(Int_t, Int_t, Int_t, Bool_t, Float_t, Float_t, Float_t, Float_t).

// 6. Configuration of the track selection.
  void SetPtSelection(Bool_t cutPt, Float_t minPt, Float_t maxPt)
  {
    this->fCutOnPt    = cutPt;
    this->fPtMin      = minPt;
    this->fPtMax      = maxPt;
  } // End: void SetPtSelection(Bool_t, Float_t, Float_t).

  void SetEtaSelection(Bool_t cutEta, Float_t minEta, Float_t maxEta)
  {
    this->fCutOnEta   = cutEta;
    this->fEtaMin     = minEta;
    this->fEtaMax     = maxEta;
  } // End: void SetEtaSelection(Bool_t, Float_t, Float_t).

  void SetNTPCSelection(Bool_t cutNTPC, Int_t minNTPC)
  {
    this->fCutOnNTPC    = cutNTPC;
    this->fNTPCMin      = minNTPC;
  } // End: void SetNTPCSelection(Bool_t, Int_t).

  void SetChiSelection(Bool_t cutChi, Float_t minChi, Float_t maxChi)
  {
    this->fCutOnChi   = cutChi;
    this->fChiMin     = minChi;
    this->fChiMax     = maxChi;
  } // End: void SetChiSelection(Bool_t, Float_t, Float_t).

  void SetNITSSelection(Bool_t cutNITS, Int_t minNITS)
  {
    this->fCutOnNITS    = cutNITS;
    this->fNITSMin      = minNITS;
  } // End: void SetNITSSelection(Bool_t, Int_t).

  void SetDCAxySelection(Bool_t cutDCAxy, Float_t maxDCAxy)
  {
    this->fCutOnDCAxy   = cutDCAxy;
    this->fDCAxyMax     = maxDCAxy;
  } // End: void SetDCAxySelection(Bool_t, Float_t).

  void SetDCAzSelection(Bool_t cutDCAz, Float_t maxDCAz)
  {
    this->fCutOnDCAz    = cutDCAz;
    this->fDCAzMax      = maxDCAz;
  } // End: void SetDCAzSelection(Bool_t, Float_t).

  void SetChargeSelection(Bool_t cutCharge, Bool_t keepPositive)
  {
    this->fCutOnCharge                = cutCharge;
    this->fKeepPositiveCharges        = keepPositive;
  } // End: void SetChargeSelection(Bool_t, Bool_t).

// 7. Configuration of the use of the weights.
  void SetParticleWeights(Bool_t useTable, Bool_t useNonUnitWeights,
    Bool_t usePt, Bool_t usePhi, Bool_t useEta, Bool_t useLocal,
    TString path, TString period)
  {
    this->fUseKineRecoTable           = useTable;
    this->fUseNonUnitParticleWeights  = useNonUnitWeights;
    this->fUsePtWeights               = usePt;
    this->fUsePhiWeights              = usePhi;
    this->fUseEtaWeights              = useEta;
    this->fUseLocalFiles              = useLocal;
    this->fPathToWeights              = path;
    this->fPeriodUsedForWeight        = period;
  } // End: void SetParticleWeights(Bool_t, Bool_t, Bool_t, Bool_t, Bool_t, Bool_t, TString).

  void SetJWeights(Bool_t useJEfficiency, Int_t indexFilter)
  {
    this->fUseJEfficiency   = useJEfficiency;
    this->fFilterbitIndex   = indexFilter;
  } // End: void SetJWeights(Bool_t, Int_t).

// 8. Configuration of the histograms and results.
  void SetReducedQvectors(Int_t powerK)
  {
    this->fReducedQPower    = powerK;
  } // End: void SetReducedQvectors(Int_t).

  void SetBinning(Int_t nBinsPt, Int_t nBinsEta, Int_t nBinsPhi)
  {
    this->fNumberBinsPt     = nBinsPt;
    this->fNumberBinsEta    = nBinsEta;
    this->fNumberBinsPhi    = nBinsPhi;
  } // End: void SetBinning(Int_t, Int_t, Int_t).

/* Methods called in the constructors. */
  virtual void  InitialiseArraysOfQvectors();
  virtual void  InitialiseArraysOfHistos();
  virtual void  InitialiseArraysOfTProfiles();

/* Methods called in "UserExec". */
  virtual void  AnalyseRecoEvent(); // Do the normal analysis at reconstructed level.
  virtual void  GetRatioDistributions();  // Get the distributions for the weights.
  Bool_t        ApplyCentralitySelection();
  Bool_t        ApplyEventSelection(Bool_t isRecoEvent);
  Bool_t        RemoveHMOs();
  Bool_t        ApplyTrackSelection(AliAODTrack *aAODtrack);
  Bool_t        ApplyTrackSelection(AliAODMCParticle *aMCtrack);
  virtual void  CalculateWeight(Int_t RunNumber, long long numberOfParticles, Float_t* pWeights, Float_t* angles, Float_t* pt, Float_t* eta);
  TH1F*         GetHistogramWithWeights(Int_t RunNumber, const char *variable);
  virtual void  CalculateQvectors(long long numberOfParticles, Float_t angles[], Float_t pWeights[]);
  TComplex      Q(Int_t n, Int_t p);
  virtual void  ComputeReducedQvectors(long long numberOfParticles);
  TComplex      CalculateRecursion(Int_t n, Int_t *harmonic, Int_t mult=1, Int_t skip=0);
  virtual void  ComputeAllCorrelators(long long numberOfParticles, Float_t angles[], Float_t pWeights[]);
  virtual void  ComputeTPCWithEtaGaps(long long numberOfParticles, Float_t angles[], Float_t pWeights[], Float_t pseudorapidity[]);

/* Methods called in "Terminate". */

/* Methods called in "UserCreateOutputObjects". */
  virtual void  BookAllLists();
  virtual void  BookMultiplicityList();
  virtual void  BookEventSelectionList();
  virtual void  BookTrackSelectionList();
  virtual void  BookMPCList();
  virtual void  BookTPCEtaList();

/* ----------------------------------------------------------------------------------------- */
private:
  AliAnalysisTaskTwoMultiCorrelations(const AliAnalysisTaskTwoMultiCorrelations& aattmc);
  AliAnalysisTaskTwoMultiCorrelations& operator=(const AliAnalysisTaskTwoMultiCorrelations& aattmc);

// 1. General parameters for the configuration of the analysis.
  TList     *fMainList;             // Mother list for the output file.
  Int_t     fHighestHarmonic;       // Highest flow harmonic to use (default: v_8).
  Int_t     fLargestCorrelators;    // Maximum number of particles in the correlators (default: 8 particles).
  Bool_t    fDoKineAnalysis;        // kTRUE: get the distributions at kine level for the weights,
    // kFALSE: compute the multi-particle correlators at reco level.
  Bool_t    fDoCorrelationsHisto;   // kTRUE: fill the 2d filter correlations histograms.
  Bool_t    fComputeEtaGaps;        // kTRUE: compute the 2-particle correlator with eta gaps.
  TComplex  fQvectors[65][9];       // All the needed combinations of Q-vectors.
    // Size: [(fHighestHarmonic*fLargestCorrelators)+1][fLargestCorrelators+1].
  AliAODEvent *fRecoEvent;          //! Input AOD event at reco level for ALICE and MC-generator AOD files.
  AliMCEvent  *fKineEvent;          //! Input MC event at kine level for MC-generator AOD files.

// 2. Parameters related to the centrality.
  TList     *fMultiplicityList;     // Daughter list for the multiplicity histograms.
  Float_t   fCentralityMin;         // Minimum of the centrality range (default: 0.0).
  Float_t   fCentralityMax;         // Maximum of the centrality range (default: 100.0).
  Bool_t    fCentralityFromSPD;     // kTRUE: use centrality estimations from SPD clusters (CL1).
  Bool_t    fCentralityFromVZero;   // kTRUE: use centrality estimations from V0A + V0C (V0M).
  Float_t   fCentrality;            //! Centrality of the current AOD event (default: -1).
  TH1I      *fHistoEvents;          //! Number of events at each selection step (1 bin per step).
  TH1F      *fHistoCentrality[2];   //! Centrality distributions [before/after the centrality selection].

// 3. Parameters related to the multiplicity.
//    [before/after the ** selection].
  long long   fInitialMultiplicity;           //! Initial number of tracks in an event (default: 0).
  TH1I        *fHistoMultiplicity[2];         //! Multiplicity distributions (** the event/track selection).
  TH1I        *fHistoMultiplicityMain[2];     //! Multiplicity distributions for the main filter bit (** HMOs).
  TH1I        *fHistoMultiplicityGlobal[2];   //! Multiplicity distributions for the global filter bit (** HMOS).
  TH2I        *fHistoCorrelatedFilters[2];    //! 2D distribution of the two filters (** HMOs).

// 4. Parameters related to the physics event selection.
//    [before/after the event selection.]
  TList     *fEventSelectionList;   // Daughter list for the event QA histograms.
  Bool_t    fCutOnPVx;              // kTRUE: apply the cuts on PV_x.
  Float_t   fPVxMin;                // Minimum PV_x [cm] (default: -44.0).
  Float_t   fPVxMax;                // Maximum PV_x [cm] (default: 44.0).
  Bool_t    fCutOnPVy;              // kTRUE: apply the cuts on PV_y.
  Float_t   fPVyMin;                // Minimum PV_y [cm] (default: -44.0).
  Float_t   fPVyMax;                // Maximum PV_y [cm] (default: 44.0).
  Bool_t    fCutOnPVz;              // kTRUE: apply the cuts on PV_z.
  Float_t   fPVzMin;                // Minimum PV_z [cm] (default: -10.0).
  Float_t   fPVzMax;                // Maximum PV_z [cm] (default: 10.0).
  TH1F      *fHistoPVx[2];          //! PV_x distributions.
  TH1F      *fHistoPVy[2];          //! PV_y distributions.
  TH1F      *fHistoPVz[2];          //! PV_z distributions.

// 5. Parameters related to the HMOs selection.
  Int_t     fMultiplicityMin;     // Minimum multiplicity needed for the event weight.
    // (Minimum non-strict -> The event can still have this multiplicity and be selected.)
  Int_t     fMainFilter;          // Main filter bit used in the analysis (default: 128).
  Int_t     fGlobalFilter;        // Second filter bit used to remove HMOs (default: 256).
  Bool_t    fCutOnHMOs;           // kTRUE: apply the cuts to remove the HMOs.
  Float_t   fMultiplicityMinA;    // a in 'a(global multiplicity) + b' for the minimum boundary.
  Float_t   fMultiplicityMinB;    // b in 'a(global multiplicity) + b' for the minimum boundary.
  Float_t   fMultiplicityMaxA;    // a in 'a(global multiplicity) + b' for the maximum boundary.
  Float_t   fMultiplicityMaxB;    // b in 'a(global multiplicity) + b' for the maximum boundary.

// 6. Parameters related to the track selection.
//    [before/after the track selecton.]
  TList     *fTrackSelectionList;   // Daughter list for the track QA histograms.
  Bool_t    fCutOnPt;               // kTRUE: apply the cuts on the transverse momentum.
  Float_t   fPtMin;                 // Minimum p_T [GeV/c] (default: 0.2).
  Float_t   fPtMax;                 // Maximum p_T [GeV/c] (default: 5.0).
  Bool_t    fCutOnEta;              // kTRUE: apply the cuts on the pseudorapidity.
  Float_t   fEtaMin;                // Minimum eta (default: -0.8).
  Float_t   fEtaMax;                // Maximum eta (default: 0.8)
  Bool_t    fCutOnNTPC;             // kTRUE: apply the cut on the number of TPC clusters.
  Int_t     fNTPCMin;               // Minimum N_TPC (default: 70).
  Bool_t    fCutOnChi;              // kTRUE: apply the cuts on chi^2 of the track momentum in TPC.
  Float_t   fChiMin;                // Minimum chi^2 (default: 0.1).
  Float_t   fChiMax;                // Maximum chi^2 (default: 4.0).
  Bool_t    fCutOnNITS;             // kTRUE: apply the cut on the number of ITS clusters.
  Int_t     fNITSMin;               // Minimum N_ITS (default: 2).
  Bool_t    fCutOnDCAxy;            // kTRUE: apply the cut on the xy-distance of DCA.
  Float_t   fDCAxyMax;              // Maximum DCA_xy [cm] (default: 2.4).
  Bool_t    fCutOnDCAz;             // kTRUE: apply the cut on the z-distance of DCA.
  Float_t   fDCAzMax;               // Maximum DCA_z [cm] (default: 3.2).
  Bool_t    fCutOnCharge;           // kTRUE: select only one type of charge.
  Bool_t    fKeepPositiveCharges;   // kTRUE: select only positive charge (fCutOnCharge must be kTRUE).
  Int_t     fNumberBinsPt;          // Binning for the pT distributions (default: 1000).
  Int_t     fNumberBinsEta;         // Binning for the eta distributions (default: 1000).
  Int_t     fNumberBinsPhi;         // Binning for the phi distributions (default: 720).
  TH1F      *fHistoPt[2];           //! pT distributions.
  TH1F      *fHistoEta[2];          //! eta distributions.
  TH1F      *fHistoPhi[2];          //! phi distributions.
  TH1I      *fHistoNTPC[2];         //! N_TPC distributions.
  TH1F      *fHistoChiSquare[2];    //! chi^2/nDF distributions.
  TH1I      *fHistoNITS[2];         //! N_ITS distributions.
  TH1F      *fHistoDCAxy[2];        //! DCA_xy distributions.
  TH1F      *fHistoDCAz[2];         //! DCA_z distributions.
  TH1I      *fHistoCharge[2];       //! Charge distributions.

// 7. Parameters related to the efficiency and acceptance weights.
  Bool_t    fUseKineRecoTable;      // kTRUE: use the kine-reco mapping table (for kine level).
  Bool_t    fUseNonUnitParticleWeights; // kTRUE: use non-unit particle weight.
  Bool_t    fUsePtWeights;          // kTRUE: use pT weights for NUE.
  Bool_t    fUsePhiWeights;         // kTRUE: use phi weights for NUA.
  Bool_t    fUseEtaWeights;         // kTRUE: use eta weights for NUA.
  Bool_t    fUseLocalFiles;         // kTRUE: use local files, kFALSE: use GRID files.
  TString   fPathToWeights;         // Path to the "Weights" directory.
    // (To use on AliEn: "alien:///alice/cern.ch/user/c/cimordas/Weights")
    // (To use locally: "/home/cindy/TestAliAnalysisTask/Test_UseWeights")
  TString   fPeriodUsedForWeight;   // Name of the period to use for the weights.
  TH1F      *fHistoWeights[3];      //! Histograms with the weights for phi, pT, eta.

  AliJEfficiency *fEfficiency;      // Used to apply NUE to the data.
  Bool_t    fFirstEvent;            ///< True if this is the first event analyzed.
  Bool_t    fUseJEfficiency;        // Use JEfficiency code to get the pT-efficiency?
  Int_t     fFilterbitIndex;        // Index used for the efficiency correction, must correspond to the main filter.
    // 0: TPCOnly; 6: hybrid (This work for AOD86, I am not sure if it is work for new AOD).
  TH1F      *fHistoEfficiency;      //! Distribution of the efficiency correction.
  TH1F      *fHistoEffInverse;      //! Distribution of the inverse of the efficiency correction.

// 8. Parameters related to the multi-particle correlations.
  TList       *fMPCList;                  // Daughter list for the multi-particle correlations techniques.
  Int_t       fReducedQPower;             // Power k for the reduced Q-vectors (default: 0).
  TH1F        *fHistoReducedQvectors[8];  //! Modulus of the reduced Q-vectors distributions for a given k.
  TProfile    *fProfileTwoPartCorrel;     //! <2>_{j,-j} for j = 1..8 (8 bins).
  TProfile    *fProfileFourPartCorrel;    //! <4>_{j,k,-j,-k} for j,k = 1..8 (36 bins).
  TProfile    *fProfileFourPartCorrelCheck; //! <4>_{j,k,-j,-k} for j,k = 1..8 for cross-check (28 bins).
  TProfile    *fProfileSixPartCorrel;     //! <6>_{j,k,l,-j,-k,-l} for j,k,l = 1..8 (56 bins).

// 9. Parameters related to the 2-particle correlations with eta gaps.
  TList       *fTPCEtaList;           // Daughter list with the 2-p correlators calculated with eta gaps.
  TProfile    *fProfileTPCEta[11];    //! <2>_{j,-j} for j = 1..8 with eta gaps (8 bins per profile).

/* ----------------------------------------------------------------------------------------- */
/* Version number to handle the objects through the iterations of the code.                  */
  ClassDef(AliAnalysisTaskTwoMultiCorrelations, 17);
};  // End of the class.

#endif


