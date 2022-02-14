#ifndef AliCascadeModule_H
#define AliCascadeModule_H
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/***********************************************

  Lambda Analysis Module - Header
  -------------------------------

This version: 27th April 2012 

--- David Dobrigkeit Chinellato
    daviddc@ifi.unicamp.br

***********************************************/

class AliCascadeModule{
public: 
  //Constructor
  AliCascadeModule(); 
  AliCascadeModule(TString ParticleType);

  //Set Files to Use
  void SetRealDataFile      ( TString RealDataFilename     );
  void SetMCDataFile        ( TString MCDataFilename       );
  void SetOutputFile        ( TString OutputFilename       );

  //Set Pt Bin Limits
  void SetPtBinLimits(Long_t got_ptbinnumb, const Double_t *got_ptbinlimits);

  //Set Rapidity Window
  void SetRapidityWindow(Double_t got_RapidityBoundaryLower, Double_t got_RapidityBoundaryUpper);
  
  //Type (either CMS or LAB)
  void SetRapidityType(TString lRapidityType = "CMS");
  
  //yCMS Shift (for pA analysis)
  void SetRapidityShift( Double_t lRapShift = 0.465);

  //Set Normalization Strategy
  void SetNormalizationStrategy ( TString lNormStrat );
  
  //Set Centrality Estimator
  void SetCentralityEstimator ( TString lCentralityEstimator );
  
  //Set CINT1B/INEL to normalize to yield
  void SetCINT1BoverINEL(Double_t got_CINT1BoverINEL);

  //Set Cuts - topological
  void SetCutV0Radius                 (Double_t cut);
  void SetCutCascRadius               (Double_t cut);
  void SetCutDCANegToPV               (Double_t cut);
  void SetCutDCAPosToPV               (Double_t cut);
  void SetCutDCABachToPV              (Double_t cut);
  void SetCutDCAV0ToPV                (Double_t cut);
  void SetCutDCACascToPV              (Double_t cut);
  void SetCutDCAV0Daughters           (Double_t cut);
  void SetCutDCACascDaughters         (Double_t cut);
  void SetCutV0CosPA                  (Double_t cut);
  void SetCutCascCosPA                (Double_t cut);
  void SetCutV0Mass                   (Double_t cut);
  void SetCutDCAxyCascToPV            (Double_t cut);
  void SetCutDCAzCascToPV             (Double_t cut);
  void SetCutDCAzPosToPV              (Double_t cut);
  void SetCutDCAzNegToPV              (Double_t cut);
  void SetCutDCAzBachToPV             (Double_t cut);
  void SetCutDCABachToBaryon          (Double_t cut);
  void SetCutBBCosPA                  (Double_t cut);
  //
  void SetCutV0Radius                 (TH1F* hCut);
  void SetCutCascRadius               (TH1F* hCut);
  void SetCutDCANegToPV               (TH1F* hCut);
  void SetCutDCAPosToPV               (TH1F* hCut);
  void SetCutDCABachToPV              (TH1F* hCut);
  void SetCutDCAV0ToPV                (TH1F* hCut);
  void SetCutDCACascToPV              (TH1F* hCut);
  void SetCutDCAV0Daughters           (TH1F* hCut);
  void SetCutDCACascDaughters         (TH1F* hCut);
  void SetCutV0CosPA                  (TH1F* hCut);
  void SetCutCascCosPA                (TH1F* hCut);
  void SetCutV0Mass                   (TH1F* hCut);
  void SetCutDCAxyCascToPV            (TH1F* hCut);
  void SetCutDCAzCascToPV             (TH1F* hCut);
  void SetCutDCAzPosToPV              (TH1F* hCut);
  void SetCutDCAzNegToPV              (TH1F* hCut);
  void SetCutDCAzBachToPV             (TH1F* hCut);
  void SetCutDCABachToBaryon          (TH1F* hCut);
  void SetCutBBCosPA                  (TH1F* hCut);

  //Set Cuts - other
  void SetCutProperLifetime           (Double_t cut);
  void SetCutCompetingSpecies         (Double_t cut);
  void SetCutTPCPIDNSigmas            (Double_t cut);
  void SetCutSigmaForSignalExtraction (Double_t cut);
  void SetCutLeastNumberOfClusters    (Double_t cut);
  void SetCutMinTrackLength           (Double_t cut);
  void SetCutMaxChi2PerCluster        (Double_t cut);
  void SetCutDaughterEta              (Double_t cut);
  void SetCutCausality                (Double_t cut);
  void SetCutTOFExpTDiffNeg           (Double_t min, Double_t max);
  void SetCutTOFExpTDiffPos           (Double_t min, Double_t max);
  void SetCutTOFExpTDiffBach          (Double_t min, Double_t max);
  void SetCutTOFSignalNeg             (Double_t min, Double_t max);
  void SetCutTOFSignalPos             (Double_t min, Double_t max);
  void SetCutTOFSignalBach            (Double_t min, Double_t max);
  //
  void SetCutProperLifetime           (TH1F* hCut);
  void SetCutCompetingSpecies         (TH1F* hCut);
  void SetCutTPCPIDNSigmas            (TH1F* hCut);
  void SetCutSigmaForSignalExtraction (TH1F* hCut);
  void SetCutLeastNumberOfClusters    (TH1F* hCut);
  void SetCutMinTrackLength           (TH1F* hCut);
  void SetCutMaxChi2PerCluster        (TH1F* hCut);
  void SetCutDaughterEta              (TH1F* hCut);
  void SetCutCausality                (TH1F* hCut);
  void SetCutTOFExpTDiffNeg           (TH1F* hCutMin, TH1F* hCutMax);
  void SetCutTOFExpTDiffPos           (TH1F* hCutMin, TH1F* hCutMax);
  void SetCutTOFExpTDiffBach          (TH1F* hCutMin, TH1F* hCutMax);
  void SetCutTOFSignalNeg             (TH1F* hCutMin, TH1F* hCutMax);
  void SetCutTOFSignalPos             (TH1F* hCutMin, TH1F* hCutMax);
  void SetCutTOFSignalBach            (TH1F* hCutMin, TH1F* hCutMax);

  //Set Use TOF match for daughter tracks
  void SetUseTOFmatchBach (Bool_t usetof = kTRUE, Double_t ptmin = 0., Double_t ptmax = 100.);
  void SetUseTOFmatchNeg  (Bool_t usetof = kTRUE, Double_t ptmin = 0., Double_t ptmax = 100.);
  void SetUseTOFmatchPos  (Bool_t usetof = kTRUE, Double_t ptmin = 0., Double_t ptmax = 100.);
  void SetUseTOFmatchOne  (Bool_t usetof = kTRUE, Double_t ptmin = 0., Double_t ptmax = 100.);

  //Set Use ITS refit for daughter tracks
  void SetUseITSrefitBach (Bool_t useits = kTRUE, Double_t ptmin = 0., Double_t ptmax = 100.);
  void SetUseITSrefitNeg  (Bool_t useits = kTRUE, Double_t ptmin = 0., Double_t ptmax = 100.);
  void SetUseITSrefitPos  (Bool_t useits = kTRUE, Double_t ptmin = 0., Double_t ptmax = 100.);
  void SetUseITSrefitOne  (Bool_t useits = kTRUE, Double_t ptmin = 0., Double_t ptmax = 100.);
  
  //Set Use ITS refit or TOF match for daughter tracks
  void SetUseITSTOFBach (Bool_t useitstof = kTRUE, Double_t ptmin = 0., Double_t ptmax = 100.);
  void SetUseITSTOFNeg  (Bool_t useitstof = kTRUE, Double_t ptmin = 0., Double_t ptmax = 100.);
  void SetUseITSTOFPos  (Bool_t useitstof = kTRUE, Double_t ptmin = 0., Double_t ptmax = 100.);
  void SetUseITSTOFOne  (Bool_t useitstof = kTRUE, Double_t ptmin = 0., Double_t ptmax = 100.);

  //Set Fit Background or not 
  void SetFitBackground ( Bool_t fitBgSwitch );

  //Set Histos with Peak Position and Width
  void SetPeakPositionAndWidthFromFit(TF1* fMeanFit, TF1* fSigmaFit);
  void SetUsePeakPositionAndWidthFromFit( Bool_t flag );

  //Use Template
  void SetUseMCBackgroundTemplate        ( Bool_t fUseMCBackgroundTemplateReceived        );
  void SetUseMinBiasMCBackgroundTemplate ( Bool_t fUseMinBiasMCBackgroundTemplateReceived );

  //Multiplicity Study Setters
  void SetPerformMultiplicityStudy ( Bool_t lPerformMultStudy );
  void SetLowMultValue             ( Double_t lLoMultBound    );
  void SetHighMultValue            ( Double_t lHiMultBound    );
  void SetMultRange                ( Double_t lLoMultBound , Double_t lHiMultBound );
  void SetUseIntegratedEfficiency  ( Bool_t lUseIntegratedEfficiency );

  //Set 
  void SetGeantFlukaCorrection( TF1* func );  

  //Set MV Pileup rejection
  void SetMVPileupRejection( Bool_t lMVPileSwitch );

  //Reject events with non empty neighboring BCs
  void SetMinDistToClosestNonEmptyBC(Int_t lNBCs); 

  //Save histograms for selection variables
  void SetSaveVarHistos( Bool_t lSaveSwitch );

  //Get Cuts - topological
  Double_t GetCutV0Radius                 (Double_t pt);
  Double_t GetCutCascRadius               (Double_t pt);
  Double_t GetCutDCANegToPV               (Double_t pt);
  Double_t GetCutDCAPosToPV               (Double_t pt);
  Double_t GetCutDCABachToPV              (Double_t pt);
  Double_t GetCutDCAV0ToPV                (Double_t pt);
  Double_t GetCutDCACascToPV              (Double_t pt);
  Double_t GetCutDCAV0Daughters           (Double_t pt);
  Double_t GetCutDCACascDaughters         (Double_t pt);
  Double_t GetCutV0CosPA                  (Double_t pt);
  Double_t GetCutCascCosPA                (Double_t pt);
  Double_t GetCutV0Mass                   (Double_t pt);
  Double_t GetCutDCAxyCascToPV            (Double_t pt);
  Double_t GetCutDCAzCascToPV             (Double_t pt);
  Double_t GetCutDCAzPosToPV              (Double_t pt);
  Double_t GetCutDCAzNegToPV              (Double_t pt);
  Double_t GetCutDCAzBachToPV             (Double_t pt);
  Double_t GetCutDCABachToBaryon          (Double_t pt);
  Double_t GetCutBBCosPA                  (Double_t pt);

  //Get Cuts - other
  Double_t GetCutProperLifetime           (Double_t pt);
  Double_t GetCutCompetingSpecies         (Double_t pt);
  Double_t GetCutTPCPIDNSigmas            (Double_t pt);
  Double_t GetCutSigmaForSignalExtraction (Double_t pt);
  Double_t GetCutLeastNumberOfClusters    (Double_t pt);
  Double_t GetCutMinTrackLength           (Double_t pt);
  Double_t GetCutMaxChi2PerCluster        (Double_t pt);
  Double_t GetCutDaughterEta              (Double_t pt);
  Double_t GetCutCausality                (Double_t pt);
  Double_t GetCutMinTOFExpTDiffNeg        (Double_t pt);
  Double_t GetCutMaxTOFExpTDiffNeg        (Double_t pt);
  Double_t GetCutMinTOFExpTDiffPos        (Double_t pt);
  Double_t GetCutMaxTOFExpTDiffPos        (Double_t pt);
  Double_t GetCutMinTOFExpTDiffBach       (Double_t pt);
  Double_t GetCutMaxTOFExpTDiffBach       (Double_t pt);
  Double_t GetCutMinTOFSignalNeg          (Double_t pt);
  Double_t GetCutMaxTOFSignalNeg          (Double_t pt);
  Double_t GetCutMinTOFSignalPos          (Double_t pt);
  Double_t GetCutMaxTOFSignalPos          (Double_t pt);
  Double_t GetCutMinTOFSignalBach         (Double_t pt);
  Double_t GetCutMaxTOFSignalBach         (Double_t pt);
   
  //Check TOF match 
  Bool_t CheckTOFmatchBach( Double_t pt, Double_t tdiff );
  Bool_t CheckTOFmatchNeg ( Double_t pt, Double_t tdiff );
  Bool_t CheckTOFmatchPos ( Double_t pt, Double_t tdiff );
  Bool_t CheckTOFmatchOne ( Double_t pt, Double_t bachtdiff, Double_t negtdiff, Double_t postdiff );

  //Check ITS refit
  Bool_t CheckITSrefitBach( Double_t pt, ULong64_t trackstatus );
  Bool_t CheckITSrefitNeg ( Double_t pt, ULong64_t trackstatus );
  Bool_t CheckITSrefitPos ( Double_t pt, ULong64_t trackstatus );
  Bool_t CheckITSrefitOne ( Double_t pt, ULong64_t btrackstatus, ULong64_t ntrackstatus, ULong64_t ptrackstatus );

  //Check ITSTOF
  Bool_t CheckITSTOFBach( Double_t pt, ULong64_t trackstatus, Double_t tdiff );
  Bool_t CheckITSTOFNeg ( Double_t pt, ULong64_t trackstatus, Double_t tdiff );
  Bool_t CheckITSTOFPos ( Double_t pt, ULong64_t trackstatus, Double_t tdiff );
  Bool_t CheckITSTOFOne ( Double_t pt, ULong64_t btrackstatus, ULong64_t ntrackstatus, ULong64_t ptrackstatus, Double_t bachtdiff, Double_t negtdiff, Double_t postdiff );


  //Get pT range for TOF match 
  Double_t GetPtMinForTOFmatchBach();
  Double_t GetPtMaxForTOFmatchBach();
  Double_t GetPtMinForTOFmatchNeg();
  Double_t GetPtMaxForTOFmatchNeg();
  Double_t GetPtMinForTOFmatchPos();
  Double_t GetPtMaxForTOFmatchPos();
  Double_t GetPtMinForTOFmatchOne();
  Double_t GetPtMaxForTOFmatchOne();

  //Get pT range for ITS refit 
  Double_t GetPtMinForITSrefitBach();
  Double_t GetPtMaxForITSrefitBach();
  Double_t GetPtMinForITSrefitNeg();
  Double_t GetPtMaxForITSrefitNeg();
  Double_t GetPtMinForITSrefitPos();
  Double_t GetPtMaxForITSrefitPos();
  Double_t GetPtMinForITSrefitOne();
  Double_t GetPtMaxForITSrefitOne();

  //Get pT range for ITSTOF 
  Double_t GetPtMinForITSTOFBach();
  Double_t GetPtMaxForITSTOFBach();
  Double_t GetPtMinForITSTOFNeg();
  Double_t GetPtMaxForITSTOFNeg();
  Double_t GetPtMinForITSTOFPos();
  Double_t GetPtMaxForITSTOFPos();
  Double_t GetPtMinForITSTOFOne();
  Double_t GetPtMaxForITSTOFOne();

  //Get Use TOF match
  Bool_t GetUseTOFmatchBach();
  Bool_t GetUseTOFmatchNeg();
  Bool_t GetUseTOFmatchPos();
  Bool_t GetUseTOFmatchOne();

  //Get Use ITS refit 
  Bool_t GetUseITSrefitBach();
  Bool_t GetUseITSrefitNeg();
  Bool_t GetUseITSrefitPos();
  Bool_t GetUseITSrefitOne();

  //Get Use ITSTOF
  Bool_t GetUseITSTOFBach();
  Bool_t GetUseITSTOFNeg();
  Bool_t GetUseITSTOFPos();
  Bool_t GetUseITSTOFOne();

  //Do Analysis
  void DoAnalysis();

  //Set Default Cuts
  void SetDefaultCuts();

  //Print Configuration
  void PrintConfiguration();

  //Auxiliary Functions 
  TString IntToString(int input);
  TString DoubleToString(double input);
  Double_t ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr );
  Double_t MyGeant3FlukaCorrectionForProtons(const Double_t *x, const Double_t *par);
  Double_t MyGeant3FlukaCorrectionForAntiProtons(const Double_t *x, const Double_t *par);
  Double_t MyLevyPtXi(const Double_t *pt, const Double_t *par);
  Double_t MyBgPol1(const Double_t *x, const Double_t *par);
  Double_t MyBgPolToEval1(const Double_t *x, const Double_t *par);
  Double_t RoundToThousandth( const Double_t lToRound );

  void SetPerformResolutionTests( Bool_t lPerformResolutionTests ); 

private:
  //root file names
  TString fRealDataFile;
  TString fMCDataFile;
  TString fOutputDataFile;
  TString fRapidityType;

  //Store Pt Bin Limits
  //Max Number of Pt Bins Set here (100)
  Double_t fptbinlimits[200];
  Double_t fptX[200];
  Long_t fptbinnumb;

  //Rapidity Range Window
  Double_t fRapidityBoundaryLower;
  Double_t fRapidityBoundaryUpper;
  Double_t fRapidityShift;

  //CINT1B/INEL for normalization
  Double_t fCINT1BoverINELratio; 
  
  //Main Analysis Parameters
  //--- 11 Topological Selections
  Double_t fCutV0Radius;   //1
  Double_t fCutCascRadius; //2
  Double_t fCutDCANegToPV; //3
  Double_t fCutDCAPosToPV; //4
  Double_t fCutDCABachToPV;//5
  Double_t fCutDCAV0Daughters; //6
  Double_t fCutV0CosPA;    //7
  Double_t fCutCascCosPA;  //8
  Double_t fCutV0Mass;     //9
  Double_t fCutDCAV0ToPV;  //10
  Double_t fCutDCACascToPV; 
  Double_t fCutDCACascDaughters; //11
  Double_t fCutDCAxyCascToPV;
  Double_t fCutDCAzCascToPV;
  Double_t fCutDCAzNegToPV;//12
  Double_t fCutDCAzPosToPV;//13
  Double_t fCutDCAzBachToPV;//14
  Double_t fCutDCABachToBaryon;//15
  Double_t fCutBBCosPA;//16
  //--- additional parameters
  Double_t fCutNegInnerP;
  Double_t fCutPosInnerP;
  Double_t fCutBachInnerP;

  //--- Causality
  Double_t fCutCausality; 
  //--- Proper Lifetime   
  Double_t fCutProperLifetime;
  //--- Competing Species Rejection; 
  Double_t fCutCompetingSpecies; 
  //--- TPC dE/dx N-sigmas
  Double_t fCutTPCPIDNSigmas;
  //--- Sigmas for signal extraction
  Double_t fCutNSigmasForSignalExtraction;
  //--- Smallest Number of Crossed Rows in TPC accepted for tracks
  Double_t fCutLeastNumberOfClusters;
  //--- Minimum Track Length
  Double_t fCutMinTrackLength;
  //--- Maximum chi2/cluster
  Double_t fCutMaxChi2PerCluster;
  //--- Daughter Track eta cut 
  Double_t fCutDaughterEta;

  //--- TOF match 
  Double_t fPtMinTOFmatchBach;
  Double_t fPtMaxTOFmatchBach;
  Double_t fPtMinTOFmatchNeg;
  Double_t fPtMaxTOFmatchNeg;
  Double_t fPtMinTOFmatchPos;
  Double_t fPtMaxTOFmatchPos;
  Double_t fPtMinTOFmatchOne;
  Double_t fPtMaxTOFmatchOne;
  Double_t fCutMinTOFExpTDiffNeg;
  Double_t fCutMaxTOFExpTDiffNeg;
  Double_t fCutMinTOFExpTDiffPos;
  Double_t fCutMaxTOFExpTDiffPos;
  Double_t fCutMinTOFExpTDiffBach;
  Double_t fCutMaxTOFExpTDiffBach;
  Double_t fCutMinTOFSignalNeg;
  Double_t fCutMaxTOFSignalNeg;
  Double_t fCutMinTOFSignalPos;
  Double_t fCutMaxTOFSignalPos;
  Double_t fCutMinTOFSignalBach;
  Double_t fCutMaxTOFSignalBach;

  //--- ITS refit 
  Double_t fPtMinITSrefitBach;
  Double_t fPtMaxITSrefitBach;
  Double_t fPtMinITSrefitNeg;
  Double_t fPtMaxITSrefitNeg;
  Double_t fPtMinITSrefitPos;
  Double_t fPtMaxITSrefitPos;
  Double_t fPtMinITSrefitOne;
  Double_t fPtMaxITSrefitOne;

  //--- ITSTOF
  Double_t fPtMinITSTOFBach;
  Double_t fPtMaxITSTOFBach;
  Double_t fPtMinITSTOFNeg;
  Double_t fPtMaxITSTOFNeg;
  Double_t fPtMinITSTOFPos;
  Double_t fPtMaxITSTOFPos;
  Double_t fPtMinITSTOFOne;
  Double_t fPtMaxITSTOFOne;

  //--- TOF match
  Bool_t fUseTOFmatchBach;
  Bool_t fUseTOFmatchNeg;
  Bool_t fUseTOFmatchPos;
  Bool_t fUseTOFmatchOne;

  //--- ITS refit
  Bool_t fUseITSrefitBach;
  Bool_t fUseITSrefitNeg;
  Bool_t fUseITSrefitPos;
  Bool_t fUseITSrefitOne;

  //--- ITSTOF
  Bool_t fUseITSTOFBach;
  Bool_t fUseITSTOFNeg;
  Bool_t fUseITSTOFPos;
  Bool_t fUseITSTOFOne;

  //--- 11 Topological Selections
  TList* fListOfPtDepCuts;
  TH1F* fHistCutV0Radius;   //1
  TH1F* fHistCutCascRadius; //2
  TH1F* fHistCutDCANegToPV; //3
  TH1F* fHistCutDCAPosToPV; //4
  TH1F* fHistCutDCABachToPV;//5
  TH1F* fHistCutDCAV0Daughters; //6
  TH1F* fHistCutV0CosPA;    //7
  TH1F* fHistCutCascCosPA;  //8
  TH1F* fHistCutV0Mass;     //9
  TH1F* fHistCutDCAV0ToPV;  //10
  TH1F* fHistCutDCACascToPV;
  TH1F* fHistCutDCACascDaughters; //11
  TH1F* fHistCutDCAxyCascToPV;
  TH1F* fHistCutDCAzCascToPV;
  TH1F* fHistCutDCAzNegToPV;//12
  TH1F* fHistCutDCAzPosToPV;//13
  TH1F* fHistCutDCAzBachToPV;//14
  TH1F* fHistCutDCABachToBaryon;//15
  TH1F* fHistCutBBCosPA;//16
  //--- additional parameters
  TH1F* fHistCutNegInnerP;
  TH1F* fHistCutPosInnerP;
  TH1F* fHistCutBachInnerP;

  //--- Causality
  TH1F* fHistCutCausality; 
  //--- Proper Lifetime   
  TH1F* fHistCutProperLifetime;
  //--- Competing Species Rejection; 
  TH1F* fHistCutCompetingSpecies; 
  //--- TPC dE/dx N-sigmas
  TH1F* fHistCutTPCPIDNSigmas;
  //--- Sigmas for signal extraction
  TH1F* fHistCutNSigmasForSignalExtraction;
  //--- Smallest Number of Crossed Rows in TPC accepted for tracks
  TH1F* fHistCutLeastNumberOfClusters;
  //--- Minimum Track Length
  TH1F* fHistCutMinTrackLength;
  //--- Maximum chi2/cluster
  TH1F* fHistCutMaxChi2PerCluster;
  //--- Daughter Track eta cut 
  TH1F* fHistCutDaughterEta;

  //--- TOF cuts 
  TH1F* fHistCutMinTOFExpTDiffNeg;
  TH1F* fHistCutMaxTOFExpTDiffNeg;
  TH1F* fHistCutMinTOFExpTDiffPos;
  TH1F* fHistCutMaxTOFExpTDiffPos;
  TH1F* fHistCutMinTOFExpTDiffBach;
  TH1F* fHistCutMaxTOFExpTDiffBach;
  TH1F* fHistCutMinTOFSignalNeg;
  TH1F* fHistCutMaxTOFSignalNeg;
  TH1F* fHistCutMinTOFSignalPos;
  TH1F* fHistCutMaxTOFSignalPos;
  TH1F* fHistCutMinTOFSignalBach;
  TH1F* fHistCutMaxTOFSignalBach;


  //WhichParticle switch: "Lambda", "AntiLambda" or "K0Short"
  TString fWhichParticle;

  //Do Fitting instead of bin counting switch
  Bool_t fFitBackgroundSwitch;

  //User provided functions for peak position and width
  TF1* fPeakPositionFit;
  TF1* fPeakWidthFit;
  Bool_t fUsePeakPositionAndWidthFromFit;

  //(Takes Precedence) Use MC background template option! 
  Bool_t fUseMCBackgroundTemplate; 
  Bool_t fUseMinBiasMCBackgroundTemplate; 

  //AntiProton Geant/Fluka correction
  TF1* fFuncGeantFlukaCorr;

  //Centrality Estimator (default V0M)
  TString fWhichEstimator;
  
  //MVPileup rejection
  Bool_t fMVPileupSwitch;

  //Minimum distance to closest non empty BC
  Int_t fMinDistToClosestNonEmptyBC;

  //Save histos switch
  Bool_t fSaveVarHistosSwitch;

  //Perform multiplicity selection 
  // --- (mult estimator in pp, centrality in PbPb)
  Bool_t fPerformMultiplicityStudy;
  Bool_t fUseIntegratedEfficiency;
  Double_t fLoMultBound;
  Double_t fHiMultBound;

  //Perform Resolution Tests (skip to save time) 
  Bool_t fPerformResolutionTests; 
  
  //Normalization Strategy Switch:
  // --- NoSelection..........: Usual method
  // --- FinalAnalysis........: Final event selection criteria
  TString fNormalizationSwitch;

};

inline void AliCascadeModule::SetUseTOFmatchBach(Bool_t usetof, Double_t ptmin, Double_t ptmax) { fUseTOFmatchBach = usetof; fPtMinTOFmatchBach = ptmin; fPtMaxTOFmatchBach = ptmax; }
inline void AliCascadeModule::SetUseTOFmatchNeg (Bool_t usetof, Double_t ptmin, Double_t ptmax) { fUseTOFmatchNeg  = usetof; fPtMinTOFmatchNeg  = ptmin; fPtMaxTOFmatchNeg  = ptmax; }
inline void AliCascadeModule::SetUseTOFmatchPos (Bool_t usetof, Double_t ptmin, Double_t ptmax) { fUseTOFmatchPos  = usetof; fPtMinTOFmatchPos  = ptmin; fPtMaxTOFmatchPos  = ptmax; }
inline void AliCascadeModule::SetUseTOFmatchOne (Bool_t usetof, Double_t ptmin, Double_t ptmax) { fUseTOFmatchOne  = usetof; fPtMinTOFmatchOne  = ptmin; fPtMaxTOFmatchOne  = ptmax; }

inline void AliCascadeModule::SetUseITSrefitBach(Bool_t useits, Double_t ptmin, Double_t ptmax) { fUseITSrefitBach = useits; fPtMinITSrefitBach = ptmin; fPtMaxITSrefitBach = ptmax; }
inline void AliCascadeModule::SetUseITSrefitNeg (Bool_t useits, Double_t ptmin, Double_t ptmax) { fUseITSrefitNeg  = useits; fPtMinITSrefitNeg  = ptmin; fPtMaxITSrefitNeg  = ptmax; }
inline void AliCascadeModule::SetUseITSrefitPos (Bool_t useits, Double_t ptmin, Double_t ptmax) { fUseITSrefitPos  = useits; fPtMinITSrefitPos  = ptmin; fPtMaxITSrefitPos  = ptmax; }
inline void AliCascadeModule::SetUseITSrefitOne (Bool_t useits, Double_t ptmin, Double_t ptmax) { fUseITSrefitOne  = useits; fPtMinITSrefitOne  = ptmin; fPtMaxITSrefitOne  = ptmax; }

inline void AliCascadeModule::SetUseITSTOFBach(Bool_t useitstof, Double_t ptmin, Double_t ptmax) { fUseITSTOFBach = useitstof; fPtMinITSTOFBach = ptmin; fPtMaxITSTOFBach = ptmax; }
inline void AliCascadeModule::SetUseITSTOFNeg (Bool_t useitstof, Double_t ptmin, Double_t ptmax) { fUseITSTOFNeg  = useitstof; fPtMinITSTOFNeg  = ptmin; fPtMaxITSTOFNeg  = ptmax; }
inline void AliCascadeModule::SetUseITSTOFPos (Bool_t useitstof, Double_t ptmin, Double_t ptmax) { fUseITSTOFPos  = useitstof; fPtMinITSTOFPos  = ptmin; fPtMaxITSTOFPos  = ptmax; }
inline void AliCascadeModule::SetUseITSTOFOne (Bool_t useitstof, Double_t ptmin, Double_t ptmax) { fUseITSTOFOne  = useitstof; fPtMinITSTOFOne  = ptmin; fPtMaxITSTOFOne  = ptmax; }

inline void AliCascadeModule::SetGeantFlukaCorrection(TF1* func) { fFuncGeantFlukaCorr = func; }

inline Bool_t   AliCascadeModule::GetUseTOFmatchBach()      { return fUseTOFmatchBach;   }
inline Bool_t   AliCascadeModule::GetUseTOFmatchNeg()       { return fUseTOFmatchNeg;    }
inline Bool_t   AliCascadeModule::GetUseTOFmatchPos()       { return fUseTOFmatchPos;    }
inline Bool_t   AliCascadeModule::GetUseTOFmatchOne()       { return fUseTOFmatchOne;    }

inline Bool_t   AliCascadeModule::GetUseITSrefitBach()      { return fUseITSrefitBach;   }
inline Bool_t   AliCascadeModule::GetUseITSrefitNeg()       { return fUseITSrefitNeg;    }
inline Bool_t   AliCascadeModule::GetUseITSrefitPos()       { return fUseITSrefitPos;    }
inline Bool_t   AliCascadeModule::GetUseITSrefitOne()       { return fUseITSrefitOne;    }

inline Bool_t   AliCascadeModule::GetUseITSTOFBach()        { return fUseITSTOFBach;   }
inline Bool_t   AliCascadeModule::GetUseITSTOFNeg()         { return fUseITSTOFNeg;    }
inline Bool_t   AliCascadeModule::GetUseITSTOFPos()         { return fUseITSTOFPos;    }
inline Bool_t   AliCascadeModule::GetUseITSTOFOne()         { return fUseITSTOFOne;    }

inline Double_t AliCascadeModule::GetPtMinForTOFmatchBach() { return fPtMinTOFmatchBach; } 
inline Double_t AliCascadeModule::GetPtMaxForTOFmatchBach() { return fPtMaxTOFmatchBach; } 
inline Double_t AliCascadeModule::GetPtMinForTOFmatchNeg()  { return fPtMinTOFmatchNeg;  } 
inline Double_t AliCascadeModule::GetPtMaxForTOFmatchNeg()  { return fPtMaxTOFmatchNeg;  } 
inline Double_t AliCascadeModule::GetPtMinForTOFmatchPos()  { return fPtMinTOFmatchPos;  } 
inline Double_t AliCascadeModule::GetPtMaxForTOFmatchPos()  { return fPtMaxTOFmatchPos;  } 
inline Double_t AliCascadeModule::GetPtMinForTOFmatchOne()  { return fPtMinTOFmatchOne;  } 
inline Double_t AliCascadeModule::GetPtMaxForTOFmatchOne()  { return fPtMaxTOFmatchOne;  } 

inline Double_t AliCascadeModule::GetPtMinForITSrefitBach() { return fPtMinITSrefitBach; } 
inline Double_t AliCascadeModule::GetPtMaxForITSrefitBach() { return fPtMaxITSrefitBach; } 
inline Double_t AliCascadeModule::GetPtMinForITSrefitNeg()  { return fPtMinITSrefitNeg;  } 
inline Double_t AliCascadeModule::GetPtMaxForITSrefitNeg()  { return fPtMaxITSrefitNeg;  } 
inline Double_t AliCascadeModule::GetPtMinForITSrefitPos()  { return fPtMinITSrefitPos;  } 
inline Double_t AliCascadeModule::GetPtMaxForITSrefitPos()  { return fPtMaxITSrefitPos;  } 
inline Double_t AliCascadeModule::GetPtMinForITSrefitOne()  { return fPtMinITSrefitOne;  } 
inline Double_t AliCascadeModule::GetPtMaxForITSrefitOne()  { return fPtMaxITSrefitOne;  } 

inline Double_t AliCascadeModule::GetPtMinForITSTOFBach()   { return fPtMinITSTOFBach; } 
inline Double_t AliCascadeModule::GetPtMaxForITSTOFBach()   { return fPtMaxITSTOFBach; } 
inline Double_t AliCascadeModule::GetPtMinForITSTOFNeg()    { return fPtMinITSTOFNeg;  } 
inline Double_t AliCascadeModule::GetPtMaxForITSTOFNeg()    { return fPtMaxITSTOFNeg;  } 
inline Double_t AliCascadeModule::GetPtMinForITSTOFPos()    { return fPtMinITSTOFPos;  } 
inline Double_t AliCascadeModule::GetPtMaxForITSTOFPos()    { return fPtMaxITSTOFPos;  } 
inline Double_t AliCascadeModule::GetPtMinForITSTOFOne()    { return fPtMinITSTOFOne;  } 
inline Double_t AliCascadeModule::GetPtMaxForITSTOFOne()    { return fPtMaxITSTOFOne;  } 

#endif


