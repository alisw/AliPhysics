
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
/ Version 15 from the 13.11.2019.                                             /
/ -------------------------------------------------------------------------- */

#ifndef ALIANALYSISTASKTWOMULTICORRELATIONS_H
#define ALIANALYSISTASKTWOMULTICORRELATIONS_H

#include "AliAnalysisTaskSE.h"
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

  void SetAnalysisParameters(Int_t maxHarmo, Int_t maxCorrelator, Bool_t doEtaGaps,
    Bool_t doTrianglePlots, Bool_t useParticleWeights)
  {
    this->fHighestHarmonic            = maxHarmo;
    this->fLargestCorrelators         = maxCorrelator;
    this->fComputeEtaGaps             = doEtaGaps;
    this->fDoCorrelationsHisto        = doTrianglePlots;
    this->fUseNonUnitParticleWeights  = useParticleWeights;
  }

  void SetTypeOfFiles(Bool_t useAOD, Bool_t useMC, Bool_t useBoth)
  {
    this->fAnalyseOnlyAOD             = useAOD;
    this->fAnalyseOnlyMC              = useMC;
    this->fAnalyseBothMCandAOD        = useBoth;
  }

  void SetCentrality(Bool_t useVZero, Bool_t useSPD, Int_t const nBins,
    Float_t minCentrality, Float_t maxCentrality)
  {
    this->fCentralityFromVZero        = useVZero;
    this->fCentralityFromSPD          = useSPD;
    this->fCentralityMin              = minCentrality;
    this->fCentralityMax              = maxCentrality;
  }

  void SetPVxyzSelection(Bool_t cutPVx, Float_t minPVx, Float_t maxPVx,
    Bool_t cutPVy, Float_t minPVy, Float_t maxPVy,
    Bool_t cutPVz, Float_t minPVz, Float_t maxPVz)
  {
    this->fCutOnPVx                   = cutPVx;
    this->fPVxMin                     = minPVx;
    this->fPVxMax                     = maxPVx;
    this->fCutOnPVy                   = cutPVy;
    this->fPVyMin                     = minPVy;
    this->fPVyMax                     = maxPVy;
    this->fCutOnPVz                   = cutPVz;
    this->fPVzMin                     = minPVz;
    this->fPVzMax                     = maxPVz;
  }

  void SetHMOsSelection(Int_t minMultiplicity, Int_t mainFilter, Int_t secondFilter,
    Int_t indexFilter, Bool_t cutHMOs, Float_t minA, Float_t minB,
    Float_t maxA, Float_t maxB)
  {
    this->fMultiplicityMin            = minMultiplicity;
    this->fMainFilter                 = mainFilter;
    this->fGlobalFilter               = secondFilter;
    this->fFilterbitIndex             = indexFilter;
    this->fCutOnHMOs                  = cutHMOs;
    this->fMultiplicityMinA           = minA;
    this->fMultiplicityMinB           = minB;
    this->fMultiplicityMaxA           = maxA;
    this->fMultiplicityMaxB           = maxB;
  }

  void SetPtSelection(Bool_t cutPt, Float_t minPt, Float_t maxPt)
  {
    this->fCutOnPt                    = cutPt;
    this->fPtMin                      = minPt;
    this->fPtMax                     = maxPt;
  }

  void SetEtaSelection(Bool_t cutEta, Float_t minEta, Float_t maxEta)
  {
    this->fCutOnEta                   = cutEta;
    this->fEtaMin                     = minEta;
    this->fEtaMax                     = maxEta;
  }

  void SetNTPCSelection(Bool_t cutNTPC, Int_t minNTPC)
  {
    this->fCutOnNTPC                  = cutNTPC;
    this->fNTPCMin                    = minNTPC;
  }

  void SetChiSelection(Bool_t cutChi, Float_t minChi, Float_t maxChi)
  {
    this->fCutOnChi                   = cutChi;
    this->fChiMin                     = minChi;
    this->fChiMax                     = maxChi;
  }

  void SetNITSSelection(Bool_t cutNITS, Int_t minNITS)
  {
    this->fCutOnNITS                  = cutNITS;
    this->fNITSMin                    = minNITS;
  }

  void SetDCAxyzSelection(Bool_t cutDCAxy, Float_t maxDCAxy,
    Bool_t cutDCAz, Float_t maxDCAz)
  {
    this->fCutOnDCAxy                 = cutDCAxy;
    this->fDCAxyMax                   = maxDCAxy;
    this->fCutOnDCAz                  = cutDCAz;
    this->fDCAzMax                    = maxDCAz;
  }

  void SetChargeSelection(Bool_t cutCharge, Bool_t keepPositive)
  {
    this->fCutOnCharge                = cutCharge;
    this->fKeepPositiveCharges        = keepPositive;
  }

  void SetReducedQvectors(Int_t powerK)
  {
    this->fReducedQPower              = powerK;
  }

/* Methods called in the constructors. */
  virtual void  InitialiseArraysOfQvectors();
  virtual void  InitialiseArraysOfHistos();
  virtual void  InitialiseArraysOfTProfiles();

/* Methods called in "UserExec". */
  virtual void  AnalyseAODevent();
  virtual void  AnalyseMCevent();
  Bool_t        ApplyPhysicsEventSelection();
  Bool_t        RemoveHMOs();
  Bool_t        ApplyTrackSelectionAOD(AliAODTrack *aAODtrack);
  Bool_t        ApplyTrackSelectionMC(AliAODMCParticle *aMCtrack);
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

/* ------------------------------------------------------------------------- */
private:
  AliAnalysisTaskTwoMultiCorrelations(const AliAnalysisTaskTwoMultiCorrelations& aattmc);
  AliAnalysisTaskTwoMultiCorrelations& operator=(const AliAnalysisTaskTwoMultiCorrelations& aattmc);

// General parameters for the analysis.
  TList       *fMainList;                   // Mother list for the output file.
  Int_t       fHighestHarmonic;             // Highest measured flow harmonic (default: v_8).
  Int_t       fLargestCorrelators;          // Maximum number of particles in the correlators (default: 8 particles).
  TComplex    fQvectors[65][9];             // All needed combinations of Q-vectors.
    // Size: [(fHighestHarmonic*fLargestCorrelators)+1][fLargestCorrelators+1].
  Bool_t      fComputeEtaGaps;              // Compute the 2-particle correlations with the eta gaps?
  Bool_t      fDoCorrelationsHisto;         // Fill the 2D correlations histograms?
  Bool_t      fUseNonUnitParticleWeights;   // Use non-unit particle weights?
  AliJEfficiency *fEfficiency;              // Used to apply NUE to the data.
  Bool_t      fFirstEvent;                  ///< True if this is the first event analyzed.

// Type of files used in the analysis.
  Bool_t      fAnalyseOnlyAOD;              // Analyse only AOD files (or Reco)?
  Bool_t      fAnalyseOnlyMC;               // Analyse only MC files (or Kine)?
  Bool_t      fAnalyseBothMCandAOD;         // Analyse both MC and AOD files?
  AliAODEvent *fAODevent;                   //! Input AOD event.
  AliMCEvent  *fMCevent;                    //! Input MC event.

// Parameters and histograms related to the centrality.
  TList       *fMultiplicityList;           // Daughter list for the multiplicity histograms.
  Bool_t      fCentralityFromVZero;         // Use V0M estimations?
  Bool_t      fCentralityFromSPD;           // Use SPD clusters estimations? 
  Float_t     fCentralityMin;               // Minimum of the centrality range.
  Float_t     fCentralityMax;               // Maximum of the centrality range.
  AliMultSelection *fMultSelection;         //! Class for the centrality estimations.
  TH1F        *fHistoCentrality[2];         //! Centrality distributions (before/after the event selection).

// Parameters and histograms related to the number of tracks.
  long long   fInitialMultiplicity;         //! Initial number of tracks in an event.
  TH1I        *fHistoMultiplicity[2];       //! Multiplicity distributions (before the event selection/after the track selection).
  TH1I        *fHistoMultiplicityMain[2];   //! Multiplicity distributions for the main filter bit (before/after the HMOs).
  TH1I        *fHistoMultiplicityGlobal[2];   //! Multiplicity distributions for the global filter bit (before/after the HMOs).
  TH2I        *fHistoCorrelatedFilters[2];    //! 2D distribution of the two filters (before/after the HMOs).

// Parameters and histograms related to the event selection.
  TList       *fEventSelectionList;         // Daughter list for the event QA histograms.
  Bool_t      fCutOnPVx;                    // Apply the cuts on PV_x?
  Float_t     fPVxMin;                      // Minimum PV_x [cm].
  Float_t     fPVxMax;                      // Maximum PV_x [cm].
  Bool_t      fCutOnPVy;                    // Apply the cuts on PV_y?
  Float_t     fPVyMin;                      // Minimum PV_y [cm].
  Float_t     fPVyMax;                      // Maximum PV_y [cm].
  Bool_t      fCutOnPVz;                    // Apply the cuts on PV_z?
  Float_t     fPVzMin;                      // Minimum PV_z [cm].
  Float_t     fPVzMax;                      // Maximum PV_z [cm].
  Int_t       fMultiplicityMin;             // Minimum multiplicity needed for the event weight.
    // (Minimum non-strict -> The event can still have this multiplicity and be selected.)
  Int_t       fMainFilter;                  // Main filter bit used in the analysis.
  Int_t       fGlobalFilter;                // Second filter bit used to remove HMOs.
  Bool_t      fCutOnHMOs;                   // Apply the cuts on the HMOs?
  Float_t     fMultiplicityMinA;            // a in 'a(global multiplicity) + b' for the minimum boundary.
  Float_t     fMultiplicityMinB;            // b in 'a(global multiplicity) + b' for the minimum boundary.
  Float_t     fMultiplicityMaxA;            // a in 'a(global multiplicity) + b' for the maximum boundary.
  Float_t     fMultiplicityMaxB;            // b in 'a(global multiplicity) + b' for the maximum boundary.
  TH1F        *fHistoPVx[2];                //! PV_x distributions (before/after the event selection).
  TH1F        *fHistoPVy[2];                //! PV_y distributions (before/after the event selection).
  TH1F        *fHistoPVz[2];                //! PV_z distributions (before/after the event selection).

// Parameters and histograms related to the track selection.
  TList       *fTrackSelectionList;         // Daughter list for the track QA histograms.
  Int_t       fFilterbitIndex;              // Index used for the efficiency correction, must correspond to the main filter.
    // 0: TPCOnly; 6: hybrid (This work for AOD86, I am not sure if it is work for new AOD).
  Bool_t      fCutOnPt;                     // Apply the cuts on p_T?
  Float_t     fPtMin;                       // Minimum p_T [GeV/c].
  Float_t     fPtMax;                       // Maximum p_T [GeV/c].
  TH1F        *fHistoEfficiency;            //! Distribution of the efficiency correction.
  TH1F        *fHistoEffInverse;            //! Distribution of the inverse of the efficiency correction.
  Bool_t      fCutOnEta;                    // Apply the cuts on eta?
  Float_t     fEtaMin;                      // Minimum eta.
  Float_t     fEtaMax;                      // Maximum eta.
  Bool_t      fCutOnNTPC;                   // Apply the cuts on N_TPC?
  Int_t       fNTPCMin;                     // Minimum N_TPC.
  Bool_t      fCutOnChi;                    // Apply the cuts on chi^2 of the track momentum in TPC?
  Float_t     fChiMin;                      // Minimum chi^2.
  Float_t     fChiMax;                      // Maximum chi^2.
  Bool_t      fCutOnNITS;                   // Apply the cut on N_ITS?
  Int_t       fNITSMin;                     // Minimum N_ITS.
  Bool_t      fCutOnDCAxy;                  // Apply the cut on DCA_xy?
  Float_t     fDCAxyMax;                    // Maximum DCA_xy [cm].
  Bool_t      fCutOnDCAz;                   // Apply the cut on DCA_z?
  Float_t     fDCAzMax;                     // Maximum DCA_z [cm].
  Bool_t      fCutOnCharge;                 // Select only one type of charges?
  Bool_t      fKeepPositiveCharges;         // Keep only the positive charges? (fCutOnCharge must be kTRUE).
  TH1F        *fHistoPt[2];                 //! pT distributions (before/after the track selection).
  TH1F        *fHistoEta[2];                //! eta distributions (before/after the track selection).
  TH1F        *fHistoPhi[2];                //! phi distributions (before/after the track selection).
  TH1I        *fHistoNTPC[2];               //! N_TPC distributions (before/after the track selection).
  TH1F        *fHistoChiSquare[2];          //! chi^2/nDF distributions (before/after the track selection).
  TH1I        *fHistoNITS[2];               //! N_ITS distributions (before/after the track selection).
  TH1F        *fHistoDCAxy[2];              //! DCA_xy distributions (before/after the track selection).
  TH1F        *fHistoDCAz[2];               //! DCA_z distributions (before/after the track selection).
  TH1I        *fHistoCharge[2];             //! Charge distributions (before/after the track selection).

// Parameters related to the multi-particle correlations.
  TList       *fMPCList;                    // Daughter list for the multi-particle correlations techniques.
  Int_t       fReducedQPower;               // Power k for the reduced Q-vectors (default: 0).
  TH1F        *fHistoReducedQvectors[8];    //! Modulus of the reduced Q-vectors distributions for a given k.
  TProfile    *fProfileTwoPartCorrel;       //! <2>_{j,-j} for j = 1..8 (8 bins).
  TProfile    *fProfileFourPartCorrel;      //! <4>_{j,k,-j,-k} for j,k = 1..8 (36 bins).
  TProfile    *fProfileFourPartCorrelCheck; //! <4>_{j,k,-j,-k} for j,k = 1..8 for cross-check (28 bins).
  TProfile    *fProfileSixPartCorrel;       //! <6>_{j,k,l,-j,-k,-l} for j,k,l = 1..8 (56 bins).

// Parameters related to the 2-particle correlations with eta gaps.
  TList       *fTPCEtaList;                 // Daughter list with the 2-p correlators calculated with eta gaps.
  TProfile    *fProfileTPCEta[11];          //! <2>_{j,-j} for j = 1..8 with eta gaps (8 bins per profile).

/* ------------------------------------------------------------------------- */
/* Version number to handle the objects through the iterations of the code.  */
  ClassDef(AliAnalysisTaskTwoMultiCorrelations, 15);
};    // End of the class.

#endif

