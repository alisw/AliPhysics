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

  Cascade Analysis Module
  ----------------------

  This Class enables the analysis of output files
  created with the

     AliAnalysisTaskStrangenessVsMultiplicity
      and
     AliAnalysisTaskStrangenessVsMultiplicityMC
 
  grid tasks.

  It constructs corrected Xi, AntiXi, Omega and
  AntiOmega spectra

  This version: late 2018

  To be compiled as a library
  (.L AliCascadeModule.cxx++ or something to that effect)

 --- Rafael Derradi de Souza
 --- David Dobrigkeit Chinellato
     daviddc@ifi.unicamp.br

 ***********************************************/

//--- For C++ ----
#include <TApplication.h>
#include <TMatrix.h>
#include <TMatrixD.h>
#include <TROOT.h>
#include <TMath.h>
#include <TF1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TLine.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TString.h>
#include <THashList.h>
#include <TDirectoryFile.h>
#include <cstdlib>
using namespace std;

//--- For ROOT ---
#include "AliVTrack.h"
#include "AliCascadeModule.h"

AliCascadeModule::AliCascadeModule()
{
    //Dummy Constructor that sets some default values
    //so that nothing too horrible will happen (hopefully)
    fWhichParticle           = "XiMinus"; //Default
    fCINT1BoverINELratio     = 1;
    fRapidityBoundaryLower   = -0.5;
    fRapidityBoundaryUpper   = 0.5;
    fRapidityType            = "LAB";
    fRapidityShift           = 0.465;
    fRealDataFile            = "";
    fMCDataFile              = "";
    fOutputDataFile          = "";

    fWhichEstimator          = "fRefMultEta5";

    //Selections
    fCutDaughterEta                  = 0.8;
    fCutV0Radius                     = -1;
    fCutCascRadius                   = -1;
    fCutV0CosPA                      = -2;
    fCutCascCosPA                    = -1;
    fCutDCANegToPV                   = -1;
    fCutDCAPosToPV                   = -1;
    fCutDCABachToPV                  = -1;
    fCutDCAV0ToPV                    = -1;
    fCutDCACascToPV                  = 1.e6;
    fCutDCAV0Daughters               = 1.e6;
    fCutDCACascDaughters             = 1.e6;
    fCutDCAxyCascToPV                = 1.e6;
    fCutDCAzCascToPV                 = 1.e6;
    fCutDCAzNegToPV                  = 1.e6;
    fCutDCAzPosToPV                  = 1.e6;
    fCutDCAzBachToPV                 = 1.e6;
    fCutV0Mass                       = -1.;
    fCutDCABachToBaryon              = -1.;
    fCutBBCosPA                      = 1.e6;
    fCutProperLifetime               = 1e+6;
    fCutTPCPIDNSigmas                = 1e+6;
    fCutLeastNumberOfClusters        = 70;
    fCutMinTrackLength               = -1.;
    fCutMaxChi2PerCluster            = 1.e9;
    fCutCompetingSpecies             = -1;
    fCutNegInnerP                    = -1.;
    fCutPosInnerP                    = -1.;
    fCutBachInnerP                   = -1.;
    fCutMinTOFExpTDiffNeg            = -99.;
    fCutMaxTOFExpTDiffNeg            =  99.;
    fCutMinTOFExpTDiffPos            = -99.;
    fCutMaxTOFExpTDiffPos            =  99.;
    fCutMinTOFExpTDiffBach           = -99.;
    fCutMaxTOFExpTDiffBach           =  99.;
    fCutMinTOFSignalNeg              = -99.;
    fCutMaxTOFSignalNeg              =  99.;
    fCutMinTOFSignalPos              = -99.;
    fCutMaxTOFSignalPos              =  99.;
    fCutMinTOFSignalBach             = -99.;
    fCutMaxTOFSignalBach             =  99.;
    fUseTOFmatchBach                 = kFALSE;
    fUseTOFmatchNeg                  = kFALSE;
    fUseTOFmatchPos                  = kFALSE;
    fUseTOFmatchOne                  = kFALSE;
    fPtMinTOFmatchBach               =   0.;
    fPtMaxTOFmatchBach               = 100.;
    fPtMinTOFmatchNeg                =   0.;
    fPtMaxTOFmatchNeg                = 100.;
    fPtMinTOFmatchPos                =   0.;
    fPtMaxTOFmatchPos                = 100.;
    fPtMinTOFmatchOne                =   0.;
    fPtMaxTOFmatchOne                = 100.;
    fUseITSrefitBach                 = kFALSE;
    fUseITSrefitNeg                  = kFALSE;
    fUseITSrefitPos                  = kFALSE;
    fUseITSrefitOne                  = kFALSE;
    fPtMinITSrefitBach               =   0.;
    fPtMaxITSrefitBach               = 100.;
    fPtMinITSrefitNeg                =   0.;
    fPtMaxITSrefitNeg                = 100.;
    fPtMinITSrefitPos                =   0.;
    fPtMaxITSrefitPos                = 100.;
    fPtMinITSrefitOne                =   0.;
    fPtMaxITSrefitOne                = 100.;
    fUseITSTOFBach                   = kFALSE;
    fUseITSTOFNeg                    = kFALSE;
    fUseITSTOFPos                    = kFALSE;
    fUseITSTOFOne                    = kFALSE;
    fPtMinITSTOFBach                 =   0.;
    fPtMaxITSTOFBach                 = 100.;
    fPtMinITSTOFNeg                  =   0.;
    fPtMaxITSTOFNeg                  = 100.;
    fPtMinITSTOFPos                  =   0.;
    fPtMaxITSTOFPos                  = 100.;
    fPtMinITSTOFOne                  =   0.;
    fPtMaxITSTOFOne                  = 100.;
    fCutNSigmasForSignalExtraction   = 4;

    fListOfPtDepCuts = new TList();
    fHistCutV0Radius                    = 0x0; 
    fHistCutCascRadius                  = 0x0; 
    fHistCutV0CosPA                     = 0x0; 
    fHistCutCascCosPA                   = 0x0; 
    fHistCutDCANegToPV                  = 0x0; 
    fHistCutDCAPosToPV                  = 0x0; 
    fHistCutDCABachToPV                 = 0x0; 
    fHistCutDCAV0ToPV                   = 0x0; 
    fHistCutDCACascToPV                 = 0x0;
    fHistCutDCAV0Daughters              = 0x0; 
    fHistCutDCACascDaughters            = 0x0; 
    fHistCutDCAxyCascToPV               = 0x0;
    fHistCutDCAzCascToPV                = 0x0;
    fHistCutDCAzNegToPV                 = 0x0; 
    fHistCutDCAzPosToPV                 = 0x0; 
    fHistCutDCAzBachToPV                = 0x0; 
    fHistCutV0Mass                      = 0x0; 
    fHistCutDCABachToBaryon             = 0x0; 
    fHistCutBBCosPA                     = 0x0; 
    fHistCutProperLifetime              = 0x0; 
    fHistCutTPCPIDNSigmas               = 0x0; 
    fHistCutLeastNumberOfClusters       = 0x0; 
    fHistCutMinTrackLength              = 0x0; 
    fHistCutMaxChi2PerCluster           = 0x0; 
    fHistCutCompetingSpecies            = 0x0; 
    fHistCutNegInnerP                   = 0x0; 
    fHistCutPosInnerP                   = 0x0; 
    fHistCutBachInnerP                  = 0x0; 
    fHistCutMinTOFExpTDiffNeg           = 0x0; 
    fHistCutMaxTOFExpTDiffNeg           = 0x0; 
    fHistCutMinTOFExpTDiffPos           = 0x0; 
    fHistCutMaxTOFExpTDiffPos           = 0x0; 
    fHistCutMinTOFExpTDiffBach          = 0x0; 
    fHistCutMaxTOFExpTDiffBach          = 0x0; 
    fHistCutMinTOFSignalNeg             = 0x0; 
    fHistCutMaxTOFSignalNeg             = 0x0; 
    fHistCutMinTOFSignalPos             = 0x0; 
    fHistCutMaxTOFSignalPos             = 0x0; 
    fHistCutMinTOFSignalBach            = 0x0; 
    fHistCutMaxTOFSignalBach            = 0x0; 
    fHistCutNSigmasForSignalExtraction  = 0x0;

    fPeakPositionFit                    = 0x0;
    fPeakWidthFit                       = 0x0;

    //Default: Use bin counting
    fFitBackgroundSwitch = kFALSE;

    fUseMCBackgroundTemplate = kFALSE;
    fUseMinBiasMCBackgroundTemplate = kTRUE;

    //Default: recompute peak positions and widths
    fUsePeakPositionAndWidthFromFit = kFALSE;

    //Default: No Causality
    fCutCausality = -1;

    //Default: Min-Bias
    fPerformMultiplicityStudy = kFALSE;
    fLoMultBound       = -1;
    fHiMultBound       = 10000;

    //Multiplicity selection for MC
    fUseIntegratedEfficiency = kTRUE;

    //Geant/Fluka correction
    fFuncGeantFlukaCorr = 0x0;

    //Default: do not check for MVPileup
    fMVPileupSwitch    = kFALSE;

    //Default: do not reject events with non empty neighbors
    fMinDistToClosestNonEmptyBC = 0;

    //Default: save histos for selection variables
    fSaveVarHistosSwitch = kTRUE;

    //Pt Bins: undefined
    fptbinnumb = -1;
    for( Int_t ipt = 0; ipt<200; ipt++) fptbinlimits[ipt] = 0;

    //Flag for resolution tests
    fPerformResolutionTests = kFALSE;
}

AliCascadeModule::AliCascadeModule(TString fParticleType)
{
    // Allows definition of Particle Type in analysis.
    // Possible Options are "XiMinus", "XiPlus", "OmegaMinus" and "OmegaPlus".
    // If some other string is given, this constructor will
    // default to "XiMinus".
    fWhichParticle = fParticleType;
    if(fWhichParticle!="XiMinus"&&fWhichParticle!="XiPlus"&&fWhichParticle!="OmegaMinus"&&fWhichParticle!="OmegaPlus"){
        cout<<"Particle Type "<<fParticleType<<" unknown. Set to XiMinus."<<endl;
        fWhichParticle = "XiMinus";
    }
    fCINT1BoverINELratio     = 1;
    fRapidityBoundaryLower   = -0.5;
    fRapidityBoundaryUpper   = 0.5;
    fRapidityType            = "LAB";
    fRapidityShift           = 0.465;
    fRealDataFile            = "";
    fMCDataFile              = "";
    fOutputDataFile          = "";

    fWhichEstimator          = "fRefMultEta5";

    //Selections
    fCutDaughterEta                  = 0.8;
    fCutV0Radius                     = -1;
    fCutCascRadius                   = -1;
    fCutV0CosPA                      = -2;
    fCutCascCosPA                    = -1;
    fCutDCANegToPV                   = -1;
    fCutDCAPosToPV                   = -1;
    fCutDCABachToPV                  = -1;
    fCutDCAV0ToPV                    = -1;
    fCutDCACascToPV                  = 1.e6;
    fCutDCAV0Daughters               = 1.e6;
    fCutDCACascDaughters             = 1.e6;
    fCutDCAxyCascToPV                = 1.e6;
    fCutDCAzCascToPV                 = 1.e6;
    fCutDCAzNegToPV                  = 1.e6;
    fCutDCAzPosToPV                  = 1.e6;
    fCutDCAzBachToPV                 = 1.e6;
    fCutV0Mass                       = -1.;
    fCutDCABachToBaryon              = -1.;
    fCutBBCosPA                      = 1.e6;
    fCutProperLifetime               = 1e+6;
    fCutTPCPIDNSigmas                = 1e+6;
    fCutLeastNumberOfClusters        = 70;
    fCutMinTrackLength               = -1.;
    fCutMaxChi2PerCluster            = 1.e9;
    fCutCompetingSpecies             = -1;
    fCutNegInnerP                    = -1.;
    fCutPosInnerP                    = -1.;
    fCutBachInnerP                   = -1.;
    fCutMinTOFExpTDiffNeg            = -99.;
    fCutMaxTOFExpTDiffNeg            =  99.;
    fCutMinTOFExpTDiffPos            = -99.;
    fCutMaxTOFExpTDiffPos            =  99.;
    fCutMinTOFExpTDiffBach           = -99.;
    fCutMaxTOFExpTDiffBach           =  99.;
    fCutMinTOFSignalNeg              = -99.;
    fCutMaxTOFSignalNeg              =  99.;
    fCutMinTOFSignalPos              = -99.;
    fCutMaxTOFSignalPos              =  99.;
    fCutMinTOFSignalBach             = -99.;
    fCutMaxTOFSignalBach             =  99.;
    fUseTOFmatchBach                 = kFALSE;
    fUseTOFmatchNeg                  = kFALSE;
    fUseTOFmatchPos                  = kFALSE;
    fUseTOFmatchOne                  = kFALSE;
    fPtMinTOFmatchBach               =   0.;
    fPtMaxTOFmatchBach               = 100.;
    fPtMinTOFmatchNeg                =   0.;
    fPtMaxTOFmatchNeg                = 100.;
    fPtMinTOFmatchPos                =   0.;
    fPtMaxTOFmatchPos                = 100.;
    fPtMinTOFmatchOne                =   0.;
    fPtMaxTOFmatchOne                = 100.;
    fUseITSrefitBach                 = kFALSE;
    fUseITSrefitNeg                  = kFALSE;
    fUseITSrefitPos                  = kFALSE;
    fUseITSrefitOne                  = kFALSE;
    fPtMinITSrefitBach               =   0.;
    fPtMaxITSrefitBach               = 100.;
    fPtMinITSrefitNeg                =   0.;
    fPtMaxITSrefitNeg                = 100.;
    fPtMinITSrefitPos                =   0.;
    fPtMaxITSrefitPos                = 100.;
    fPtMinITSrefitOne                =   0.;
    fPtMaxITSrefitOne                = 100.;
    fUseITSTOFBach                   = kFALSE;
    fUseITSTOFNeg                    = kFALSE;
    fUseITSTOFPos                    = kFALSE;
    fUseITSTOFOne                    = kFALSE;
    fPtMinITSTOFBach                 =   0.;
    fPtMaxITSTOFBach                 = 100.;
    fPtMinITSTOFNeg                  =   0.;
    fPtMaxITSTOFNeg                  = 100.;
    fPtMinITSTOFPos                  =   0.;
    fPtMaxITSTOFPos                  = 100.;
    fPtMinITSTOFOne                  =   0.;
    fPtMaxITSTOFOne                  = 100.;
    fCutNSigmasForSignalExtraction   = 4;

    fListOfPtDepCuts = new TList();
    fHistCutV0Radius                    = 0x0; 
    fHistCutCascRadius                  = 0x0; 
    fHistCutV0CosPA                     = 0x0; 
    fHistCutCascCosPA                   = 0x0; 
    fHistCutDCANegToPV                  = 0x0; 
    fHistCutDCAPosToPV                  = 0x0; 
    fHistCutDCABachToPV                 = 0x0; 
    fHistCutDCAV0ToPV                   = 0x0; 
    fHistCutDCACascToPV                 = 0x0;
    fHistCutDCAV0Daughters              = 0x0; 
    fHistCutDCACascDaughters            = 0x0; 
    fHistCutDCAxyCascToPV               = 0x0;
    fHistCutDCAzCascToPV                = 0x0;
    fHistCutDCAzNegToPV                 = 0x0; 
    fHistCutDCAzPosToPV                 = 0x0; 
    fHistCutDCAzBachToPV                = 0x0; 
    fHistCutV0Mass                      = 0x0; 
    fHistCutDCABachToBaryon             = 0x0; 
    fHistCutBBCosPA                     = 0x0; 
    fHistCutProperLifetime              = 0x0; 
    fHistCutTPCPIDNSigmas               = 0x0; 
    fHistCutLeastNumberOfClusters       = 0x0; 
    fHistCutMinTrackLength              = 0x0; 
    fHistCutMaxChi2PerCluster           = 0x0; 
    fHistCutCompetingSpecies            = 0x0; 
    fHistCutNegInnerP                   = 0x0; 
    fHistCutPosInnerP                   = 0x0; 
    fHistCutBachInnerP                  = 0x0; 
    fHistCutMinTOFExpTDiffNeg           = 0x0; 
    fHistCutMaxTOFExpTDiffNeg           = 0x0; 
    fHistCutMinTOFExpTDiffPos           = 0x0; 
    fHistCutMaxTOFExpTDiffPos           = 0x0; 
    fHistCutMinTOFExpTDiffBach          = 0x0; 
    fHistCutMaxTOFExpTDiffBach          = 0x0; 
    fHistCutMinTOFSignalNeg             = 0x0; 
    fHistCutMaxTOFSignalNeg             = 0x0; 
    fHistCutMinTOFSignalPos             = 0x0; 
    fHistCutMaxTOFSignalPos             = 0x0; 
    fHistCutMinTOFSignalBach            = 0x0; 
    fHistCutMaxTOFSignalBach            = 0x0; 
    fHistCutNSigmasForSignalExtraction  = 0x0;

    fPeakPositionFit                    = 0x0;
    fPeakWidthFit                       = 0x0;

    //Default: Use bin counting
    fFitBackgroundSwitch = kFALSE;

    fUseMCBackgroundTemplate = kFALSE;
    fUseMinBiasMCBackgroundTemplate = kTRUE;

    //Default: recompute peak positions and widths
    fUsePeakPositionAndWidthFromFit = kFALSE;

    //Default: No Causality
    fCutCausality = -1;

    //Default: Min-Bias
    fPerformMultiplicityStudy = kFALSE;
    fLoMultBound       = -1;
    fHiMultBound       = 10000;

    //Multiplicity selection for MC
    fUseIntegratedEfficiency = kTRUE;

    //Geant/Fluka correction
    fFuncGeantFlukaCorr = 0x0;

    //Default: do not check for MVPileup
    fMVPileupSwitch    = kFALSE;

    //Default: do not reject events with non empty neighbors
    fMinDistToClosestNonEmptyBC = 0;

    //Default: save histos for selection variables
    fSaveVarHistosSwitch = kTRUE;

    //Pt Bins: undefined
    fptbinnumb = -1;
    for( Int_t ipt = 0; ipt<200; ipt++) fptbinlimits[ipt] = 0;

    //Flag for resolution tests
    fPerformResolutionTests = kFALSE;
}

/***********************************************
  --- Setters For Configuration ---
 ***********************************************/

// Filename Setters
void AliCascadeModule::SetRealDataFile    ( TString RealDataFilename     ){
    //Set root file containing real data candidates.
    fRealDataFile = RealDataFilename;
}
void AliCascadeModule::SetMCDataFile      ( TString MCDataFilename       ){
    //Set root file containing Monte Carlo data (for efficiency computation).
    fMCDataFile   = MCDataFilename;
}
void AliCascadeModule::SetOutputFile      ( TString OutputFilename       ){
    //Set root filename for the analysis output.
    fOutputDataFile = OutputFilename;
    cout<<"[AliCascadeModule] Set output file to \'"<<OutputFilename<<"\'."<<endl;
}

// Bin Limit Setter
void AliCascadeModule::SetPtBinLimits(Long_t got_fptbinnumb,const Double_t *got_fptbinlimits){
    //Function to set pt binning. First argument is the number of pt bins, second is
    //an array with bin limits.
    fptbinnumb = got_fptbinnumb;
    for(int ix = 0;ix<fptbinnumb+1;ix++){
        fptbinlimits[ix] = got_fptbinlimits[ix];
    }
    for(int ix = 0;ix<fptbinnumb;ix++){
        fptX[ix] = (fptbinlimits[ix+1] + fptbinlimits[ix])/2.;
    }
    cout<<"[AliCascadeModule] Received "<<fptbinnumb<<" pt bins, set accordingly."<<endl;
}

// Rapidity Window Setter
void AliCascadeModule::SetRapidityWindow(Double_t got_RapidityBoundaryLower, Double_t got_RapidityBoundaryUpper){
    //Set Rapidity Boundary used in analysis.
    //Value provided will be used as upper limit for the modulus of y (|y|< given value)
    fRapidityBoundaryLower = got_RapidityBoundaryLower;
    fRapidityBoundaryUpper = got_RapidityBoundaryUpper;
    cout<<"[AliCascadeModule] Received "<<got_RapidityBoundaryLower<<" to "<<got_RapidityBoundaryUpper<<" as rapidity limits, set accordingly."<<endl;
}

// Rapidity Type Setter
void AliCascadeModule::SetRapidityType(TString lRapidityType){
    //Set Rapidity Type used: either
    // - CMS (center of mass, *not* the experiment)
    // - LAB (Laboratory reference)
    fRapidityType = lRapidityType;
    cout<<"[AliCascadeModule] Received rapidity type: Y("<<lRapidityType<<") for analysis, set accordingly."<<endl;
}

// Rapidity Shift Setter
void AliCascadeModule::SetRapidityShift( Double_t lRapShift){
    //Set Rapidity shift used for yCMS calculation
    fRapidityShift = lRapShift;
    cout<<"[AliCascadeModule] Received rapidity shift: "<<fRapidityShift<<" for analysis, set accordingly."<<endl;
}

void AliCascadeModule::SetNormalizationStrategy ( TString lNormStrat ){
    //Set method used to compute normalization
    fNormalizationSwitch = lNormStrat;
    cout<<"[AliCascadeModule] Received Normalization treatment method: "<<lNormStrat<<endl;
}

void AliCascadeModule::SetCentralityEstimator ( TString lCentralityEstimator ){
    //Set centrality estimator (one of CL1, TRK, V0M, V0A, ZNA)
    fWhichEstimator = lCentralityEstimator;
    cout<<"[AliCascadeModule] Received Centrality Estimator: "<<lCentralityEstimator<<endl;
}

// CINT1B/INEL Setter for normalization to yields
void AliCascadeModule::SetCINT1BoverINEL(Double_t got_fCINT1BoverINEL){
    //Set CINT1B/INEL ratio (determined by collaboration).
    fCINT1BoverINELratio = got_fCINT1BoverINEL;
    cout<<"[AliCascadeModule] Received CINT1B/INEL = "<<got_fCINT1BoverINEL<<" to normalize with, set accordingly."<<endl;
}

// Topological Selection Setters //////////////////////////////////////////////
//1
void AliCascadeModule::SetCutV0Radius(Double_t cut){
    //Set minimum decay radius for the V0 in centimeters.
    //Note: this is (R_{2D}, cylindrical, centered around ALICE detector)
    fCutV0Radius = cut;
    cout<<"[AliCascadeModule] Received V0 Radius (min value) = "<<cut<<endl;
}
//2
void AliCascadeModule::SetCutCascRadius(Double_t cut){
    //Set minimum decay radius for the V0 in centimeters.
    //Note: this is (R_{2D}, cylindrical, centered around ALICE detector)
    fCutCascRadius = cut;
    cout<<"[AliCascadeModule] Received Casc Radius (min value) = "<<cut<<endl;
}
//3
void AliCascadeModule::SetCutDCANegToPV(Double_t cut){
    //Set minimum distance of closest approach between V0 negative daughter
    //track and primary vertex (in centimeters).
    fCutDCANegToPV = cut;
    cout<<"[AliCascadeModule] Received DCA Negative track to PV (min value) = "<<cut<<endl;
}
//4
void AliCascadeModule::SetCutDCAPosToPV(Double_t cut){
    //Set minimum distance of closest approach between V0 positive daughter
    //track and primary vertex (in centimeters).
    fCutDCAPosToPV = cut;
    cout<<"[AliCascadeModule] Received DCA Positive track to PV (min value) = "<<cut<<endl;
}
//5
void AliCascadeModule::SetCutDCABachToPV(Double_t cut){
    //Set minimum distance of closest approach between V0 negative daughter
    //track and primary vertex (in centimeters).
    fCutDCABachToPV = cut;
    cout<<"[AliCascadeModule] Received DCA Bachelor track to PV (min value) = "<<cut<<endl;
}
//6
void AliCascadeModule::SetCutDCAV0Daughters(Double_t cut){
    //Set minimum distance of closest approach between V0 daughter
    //tracks, in sigmas. This is not in centimeters because if the
    //tracks have been determined with ITS refit the resolution is
    //greatly improved; thus, the cut may be tighter in that case.
    //Using sigmas will therefore take the tracking method in
    //consideration.
    fCutDCAV0Daughters = cut;
    cout<<"[AliCascadeModule] Received DCA V0 Daughters (max value) = "<<cut<<endl;
}
//7
void AliCascadeModule::SetCutDCACascDaughters(Double_t cut){
    //Set minimum distance of closest approach between Cascade daughters.
    fCutDCACascDaughters = cut;
    cout<<"[AliCascadeModule] Received DCA Casc Daughters (max value) = "<<cut<<endl;
}
//8
void AliCascadeModule::SetCutV0CosPA(Double_t cut){
    //Set minimum value for the cosine of pointing angle of the V0.
    fCutV0CosPA = cut;
    cout<<"[AliCascadeModule] Received V0 Cosine of Pointing Angle (min value) = "<<cut<<endl;
}
//9
void AliCascadeModule::SetCutCascCosPA(Double_t cut){
    //Set minimum value for the cosine of pointing angle of the cascade.
    fCutCascCosPA = cut;
    cout<<"[AliCascadeModule] Received Casc Cosine of Pointing Angle (min value) = "<<cut<<endl;
}
//10
void AliCascadeModule::SetCutV0Mass(Double_t cut){
    //Set mass window width for Lambdas (in GeV/c^2)
    fCutV0Mass = cut;
    cout<<"[AliCascadeModule] Received Lambda invariant mass window width = "<<cut<<endl;
}
//11
void AliCascadeModule::SetCutDCAV0ToPV(Double_t cut){
    //Set minimum distance of closest approach between V0 and PV
    fCutDCAV0ToPV = cut;
    cout<<"[AliCascadeModule] Received DCA V0 to PV (min value) = "<<cut<<endl;
}
//
void AliCascadeModule::SetCutDCACascToPV(Double_t cut){
    //Set maximum distance of closest approach between V0 and PV
    fCutDCACascToPV = cut;
    cout<<"[AliCascadeModule] Received DCA Cascade to PV (max value) = "<<cut<<endl;
}
//
void AliCascadeModule::SetCutDCAxyCascToPV(Double_t cut){
    //Set maximum distance of closest approach between V0 and PV in the xy plane
    fCutDCAxyCascToPV = cut;
    cout<<"[AliCascadeModule] Received DCAxy Cascade to PV (max value) = "<<cut<<endl;
}
//
void AliCascadeModule::SetCutDCAzCascToPV(Double_t cut){
    //Set maximum distance of closest approach between V0 and PV along the z direction
    fCutDCAzCascToPV = cut;
    cout<<"[AliCascadeModule] Received DCAz Cascade to PV (max value) = "<<cut<<endl;
}
//12
void AliCascadeModule::SetCutDCAzNegToPV(Double_t cut) {
    //Set maximum distance of closest approach along z direction between V0 negative daughter
    //track and primary vertex (in centimeters).
    fCutDCAzNegToPV = cut;
    cout<<"[AliCascadeModule] Received DCAz Negative track to PV (max value) = "<<cut<<endl;
}
//13
void AliCascadeModule::SetCutDCAzPosToPV(Double_t cut) {
    //Set maximum distance of closest approach along z direction between V0 positive daughter
    //track and primary vertex (in centimeters).
    fCutDCAzPosToPV = cut;
    cout<<"[AliCascadeModule] Received DCAz Positive track to PV (max value) = "<<cut<<endl;
}
//14
void AliCascadeModule::SetCutDCAzBachToPV(Double_t cut){
    //Set maximum distance of closest approach along z direction between bachelor
    //track and primary vertex (in centimeters).
    fCutDCAzBachToPV = cut;
    cout<<"[AliCascadeModule] Received DCAz Bachelor track to PV (max value) = "<<cut<<endl;
}
//15
void AliCascadeModule::SetCutDCABachToBaryon(Double_t cut){
    //Set minimum distance between bachelor and baryon daughter tracks
    //This cut is particularly useful to cure the bump structure in Xi analysis
    fCutDCABachToBaryon = cut;
    cout<<"[AliCascadeModule] Received DCA Bachelor to Baryon (min value) = "<<cut<<endl;
}
//16
void AliCascadeModule::SetCutBBCosPA(Double_t cut){
    //Set maximum value for the cosine of pointing angle from bachelor-baryon (wrong vertex)
    //This cut is particularly useful to cure the bump structure in Xi analysis
    fCutBBCosPA = cut;
    cout<<"[AliCascadeModule] Received Bachelor-Baryon Cosine of Pointing Angle (max value) = "<<cut<<endl;
}
///////////////////////////////////////////////////////////////////////////////

// pT-dep Topological Selection Setters ///////////////////////////////////////
//1
void AliCascadeModule::SetCutV0Radius(TH1F* hCut){
    //Set minimum decay radius for the V0 in centimeters.
    //Note: this is (R_{2D}, cylindrical, centered around ALICE detector)
    fHistCutV0Radius = hCut;
    fListOfPtDepCuts->Add(fHistCutV0Radius);
    cout<<"[AliCascadeModule] Received pt-dep V0 Radius (min value)"<<endl;
}
//2
void AliCascadeModule::SetCutCascRadius(TH1F* hCut){
    //Set minimum decay radius for the V0 in centimeters.
    //Note: this is (R_{2D}, cylindrical, centered around ALICE detector)
    fHistCutCascRadius = hCut;
    fListOfPtDepCuts->Add(fHistCutCascRadius);
    cout<<"[AliCascadeModule] Received pt-dep Casc Radius (min value)"<<endl;
}
//3
void AliCascadeModule::SetCutDCANegToPV(TH1F* hCut){
    //Set minimum distance of closest approach between V0 negative daughter
    //track and primary vertex (in centimeters).
    fHistCutDCANegToPV = hCut;
    fListOfPtDepCuts->Add(fHistCutDCANegToPV);
    cout<<"[AliCascadeModule] Received pt-dep DCA Negative track to PV (min value)"<<endl;
}
//4
void AliCascadeModule::SetCutDCAPosToPV(TH1F* hCut){
    //Set minimum distance of closest approach between V0 positive daughter
    //track and primary vertex (in centimeters).
    fHistCutDCAPosToPV = hCut;
    fListOfPtDepCuts->Add(fHistCutDCAPosToPV);
    cout<<"[AliCascadeModule] Received pt-dep DCA Positive track to PV (min value)"<<endl;
}
//5
void AliCascadeModule::SetCutDCABachToPV(TH1F* hCut){
    //Set minimum distance of closest approach between V0 negative daughter
    //track and primary vertex (in centimeters).
    fHistCutDCABachToPV = hCut;
    fListOfPtDepCuts->Add(fHistCutDCABachToPV);
    cout<<"[AliCascadeModule] Received pt-dep DCA Bachelor track to PV (min value)"<<endl;
}
//6
void AliCascadeModule::SetCutDCAV0Daughters(TH1F* hCut){
    //Set minimum distance of closest approach between V0 daughter
    //tracks, in sigmas. This is not in centimeters because if the
    //tracks have been determined with ITS refit the resolution is
    //greatly improved; thus, the cut may be tighter in that case.
    //Using sigmas will therefore take the tracking method in
    //consideration.
    fHistCutDCAV0Daughters = hCut;
    fListOfPtDepCuts->Add(fHistCutDCAV0Daughters);
    cout<<"[AliCascadeModule] Received pt-dep DCA V0 Daughters (max value)"<<endl;
}
//7
void AliCascadeModule::SetCutDCACascDaughters(TH1F* hCut){
    //Set minimum distance of closest approach between Cascade daughters.
    fHistCutDCACascDaughters = hCut;
    fListOfPtDepCuts->Add(fHistCutDCACascDaughters);
    cout<<"[AliCascadeModule] Received pt-dep DCA Casc Daughters (max value)"<<endl;
}
//8
void AliCascadeModule::SetCutV0CosPA(TH1F* hCut){
    //Set minimum value for the cosine of pointing angle of the V0.
    fHistCutV0CosPA = hCut;
    fListOfPtDepCuts->Add(fHistCutV0CosPA);
    cout<<"[AliCascadeModule] Received pt-dep V0 Cosine of Pointing Angle (min value)"<<endl;
}
//9
void AliCascadeModule::SetCutCascCosPA(TH1F* hCut){
    //Set minimum value for the cosine of pointing angle of the cascade.
    fHistCutCascCosPA = hCut;
    fListOfPtDepCuts->Add(fHistCutCascCosPA);
    cout<<"[AliCascadeModule] Received pt-dep Casc Cosine of Pointing Angle (min value)"<<endl;
}
//10
void AliCascadeModule::SetCutV0Mass(TH1F* hCut){
    //Set mass window width for Lambdas (in GeV/c^2)
    fHistCutV0Mass = hCut;
    fListOfPtDepCuts->Add(fHistCutV0Mass);
    cout<<"[AliCascadeModule] Received pt-dep Lambda invariant mass window width"<<endl;
}
//11
void AliCascadeModule::SetCutDCAV0ToPV(TH1F* hCut){
    //Set minimum distance of closest approach between V0 and PV
    fHistCutDCAV0ToPV = hCut;
    fListOfPtDepCuts->Add(fHistCutDCAV0ToPV);
    cout<<"[AliCascadeModule] Received pt-dep DCA V0 to PV (min value)"<<endl;
}
//
void AliCascadeModule::SetCutDCACascToPV(TH1F* hCut){
    //Set maximum distance of closest approach between Cascade and PV
    fHistCutDCACascToPV = hCut;
    fListOfPtDepCuts->Add(fHistCutDCACascToPV);
    cout<<"[AliCascadeModule] Received pt-dep DCA Cascade to PV (max value)"<<endl;
}
//
void AliCascadeModule::SetCutDCAxyCascToPV(TH1F* hCut){
    //Set maximum distance of closest approach between Cascade and PV in the xy plane
    fHistCutDCAxyCascToPV = hCut;
    fListOfPtDepCuts->Add(fHistCutDCAxyCascToPV);
    cout<<"[AliCascadeModule] Received pt-dep DCAxy Cascade to PV (max value)"<<endl;
}
//
void AliCascadeModule::SetCutDCAzCascToPV(TH1F* hCut){
    //Set maximum distance of closest approach between Cascade and PV along the z direction
    fHistCutDCAzCascToPV = hCut;
    fListOfPtDepCuts->Add(fHistCutDCAzCascToPV);
    cout<<"[AliCascadeModule] Received pt-dep DCAz Cascade to PV (max value)"<<endl;
}
//12
void AliCascadeModule::SetCutDCAzNegToPV(TH1F* hCut) {
    //Set maximum distance of closest approach along z direction between V0 negative daughter
    //track and primary vertex (in centimeters).
    fHistCutDCAzNegToPV = hCut;
    fListOfPtDepCuts->Add(fHistCutDCAzNegToPV);
    cout<<"[AliCascadeModule] Received pt-dep DCAz Negative track to PV (max value)"<<endl;
}
//13
void AliCascadeModule::SetCutDCAzPosToPV(TH1F* hCut) {
    //Set maximum distance of closest approach along z direction between V0 positive daughter
    //track and primary vertex (in centimeters).
    fHistCutDCAzPosToPV = hCut;
    fListOfPtDepCuts->Add(fHistCutDCAzPosToPV);
    cout<<"[AliCascadeModule] Received pt-dep DCAz Positive track to PV (max value)"<<endl;
}
//14
void AliCascadeModule::SetCutDCAzBachToPV(TH1F* hCut){
    //Set maximum distance of closest approach along z direction between bachelor
    //track and primary vertex (in centimeters).
    fHistCutDCAzBachToPV = hCut;
    fListOfPtDepCuts->Add(fHistCutDCAzBachToPV);
    cout<<"[AliCascadeModule] Received pt-dep DCAz Bachelor track to PV (max value)"<<endl;
}
//15
void AliCascadeModule::SetCutDCABachToBaryon(TH1F* hCut){
    //Set minimum distance between bachelor and baryon daughter tracks
    //This cut is particularly useful to cure the bump structure in Xi analysis
    fHistCutDCABachToBaryon = hCut;
    fListOfPtDepCuts->Add(fHistCutDCABachToBaryon);
    cout<<"[AliCascadeModule] Received pt-dep DCA Bachelor to Baryon (min value)"<<endl;
}
//16
void AliCascadeModule::SetCutBBCosPA(TH1F* hCut){
    //Set maximum value for the cosine of pointing angle from bachelor-baryon (wrong vertex)
    //This cut is particularly useful to cure the bump structure in Xi analysis
    fHistCutBBCosPA = hCut;
    fListOfPtDepCuts->Add(fHistCutBBCosPA);
    cout<<"[AliCascadeModule] Received pt-dep Bachelor-Baryon Cosine of Pointing Angle (max value)"<<endl;
}
///////////////////////////////////////////////////////////////////////////////

// Other Selection Setters ////////////////////////////////////////////////////
void AliCascadeModule::SetCutProperLifetime(Double_t cut){
    //Set maximum value for m*L/p variable for the V0.
    //This is the "proper lifetime" selection and is usually called a
    //"c*tau cut". Should be set to a value larger than the c*tau for
    //the V0 considered.
    fCutProperLifetime = cut;
    cout<<"[AliCascadeModule] Received proper lifetime cut (max value) = "<<cut<<endl;
}

void AliCascadeModule::SetCutCompetingSpecies(Double_t cut){
    //Set minimum distance in xi invariant mass space for Omega Analysis.
    fCutCompetingSpecies = cut;
    cout<<"[AliCascadeModule] Received competing species cut (min value) = "<<cut<<endl;
}
void AliCascadeModule::SetCutTPCPIDNSigmas(Double_t cut){
    //Set maximum deviation from the expected energy loss in
    //the TPC, in multiples of sigmas as computed from the AliPIDResponse
    //object. Selection is only used in real data and should thus
    //be very loose to ensure negligible signal loss.
    fCutTPCPIDNSigmas = cut;
    cout<<"[AliCascadeModule] Received TPC N-sigmas selection (max dist from BB curve) = "<<cut<<endl;
}
void AliCascadeModule::SetCutSigmaForSignalExtraction(Double_t cut){
    //Set number of sigmas for the signal extraction method. The value
    //provided is the number of sigmas from the peak position used for
    //the peak sampling, i.e. the peak will be from [<m>-cut*sigma,<m>+cut*sigma]
    //while the background sampling regions will be from
    //[<m>-2*cut*sigma,<m>+2*cut*sigma.
    fCutNSigmasForSignalExtraction = cut;
    cout<<"[AliCascadeModule] Received N-sigmas for sig. ext.: peak is (-"<<cut<<",+"<<cut<<") in sigmas"<<endl;
}
void AliCascadeModule::SetCutLeastNumberOfClusters(Double_t cut){
    //Set smallest allowed number of TPC clusters for the V0 daughter tracks.
    fCutLeastNumberOfClusters = cut;
    cout<<"[AliCascadeModule] Received Least Nbr of crossed rows (min value) = "<<cut<<endl;
}
void AliCascadeModule::SetCutMinTrackLength(Double_t cut) {
    //Set minimum lenght allowed for Cascade daughter tracks
    fCutMinTrackLength = cut;
    cout<<"[AliCascadeModule] Received Daughter track length cut (min value) = "<<cut<<endl;
}
void AliCascadeModule::SetCutMaxChi2PerCluster(Double_t cut) {
    //Set maximum chi2/cluster for the track
    fCutMaxChi2PerCluster = cut;
    cout<<"[AliCascadeModule] Received chi2/cluster cut (max value) = "<<cut<<endl;
}
void AliCascadeModule::SetCutDaughterEta(Double_t cut){
    //Set maximum eta value allowed for the V0 daughter tracks (|eta|<cut).
    fCutDaughterEta = cut;
    cout<<"[AliCascadeModule] Received Daughter |eta| cut (max value) = "<<cut<<endl;
}

void AliCascadeModule::SetCutCausality(Double_t cut){
    //Set maximum eta value allowed for the V0 daughter tracks (|eta|<cut).
    fCutCausality = cut;
    cout<<"[AliCascadeModule] Received causality cut (min value) = "<<cut<<endl;
}
///////////////////////////////////////////////////////////////////////////////

// pT-dep - Other Selection Setters ///////////////////////////////////////////
void AliCascadeModule::SetCutProperLifetime(TH1F* hCut){
    //Set maximum value for m*L/p variable for the V0.
    //This is the "proper lifetime" selection and is usually called a
    //"c*tau cut". Should be set to a value larger than the c*tau for
    //the V0 considered.
    fHistCutProperLifetime = hCut;
    fListOfPtDepCuts->Add(fHistCutProperLifetime);
    cout<<"[AliCascadeModule] Received pt-dep proper lifetime cut (max value)"<<endl;
}

void AliCascadeModule::SetCutCompetingSpecies(TH1F* hCut){
    //Set minimum distance in xi invariant mass space for Omega Analysis.
    fHistCutCompetingSpecies = hCut;
    fListOfPtDepCuts->Add(fHistCutCompetingSpecies);
    cout<<"[AliCascadeModule] Received pt-dep competing species cut (min value)"<<endl;
}
void AliCascadeModule::SetCutTPCPIDNSigmas(TH1F* hCut){
    //Set maximum deviation from the expected energy loss in
    //the TPC, in multiples of sigmas as computed from the AliPIDResponse
    //object. Selection is only used in real data and should thus
    //be very loose to ensure negligible signal loss.
    fHistCutTPCPIDNSigmas = hCut;
    fListOfPtDepCuts->Add(fHistCutTPCPIDNSigmas);
    cout<<"[AliCascadeModule] Received pt-dep TPC N-sigmas selection (max dist from BB curve)"<<endl;
}
void AliCascadeModule::SetCutSigmaForSignalExtraction(TH1F* hCut){
    //Set number of sigmas for the signal extraction method. The value
    //provided is the number of sigmas from the peak position used for
    //the peak sampling, i.e. the peak will be from [<m>-cut*sigma,<m>+cut*sigma]
    //while the background sampling regions will be from
    //[<m>-2*cut*sigma,<m>+2*cut*sigma.
    fHistCutNSigmasForSignalExtraction = hCut;
    fListOfPtDepCuts->Add(fHistCutNSigmasForSignalExtraction);
    cout<<"[AliCascadeModule] Received pt-dep N-sigmas for sig. ext."<<endl;
}
void AliCascadeModule::SetCutLeastNumberOfClusters(TH1F* hCut){
    //Set smallest allowed number of TPC clusters for the V0 daughter tracks.
    fHistCutLeastNumberOfClusters = hCut;
    fListOfPtDepCuts->Add(fHistCutLeastNumberOfClusters);
    cout<<"[AliCascadeModule] Received pt-dep Least Nbr of crossed rows (min value)"<<endl;
}
void AliCascadeModule::SetCutMinTrackLength(TH1F* hCut) {
    //Set pT-dep minimum track length allowed for Cascade daughter tracks.
    fHistCutMinTrackLength = hCut;
    fHistCutMinTrackLength->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fHistCutMinTrackLength->GetYaxis()->SetTitle("MinTrackLength");
    cout<<"[AliCascadeModule] Received pT-dep track length cut (min value)"<<endl;
}
void AliCascadeModule::SetCutMaxChi2PerCluster(TH1F* hCut) {
    //Set pT-dep maximum chi2/cluster for Cascade daughter tracks.
    fHistCutMaxChi2PerCluster = hCut;
    fHistCutMaxChi2PerCluster->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    fHistCutMaxChi2PerCluster->GetYaxis()->SetTitle("MaxChi2PerCluster");
    cout<<"[AliCascadeModule] Received pT-dep chi2/cluster cut (max value)"<<endl;
}
void AliCascadeModule::SetCutDaughterEta(TH1F* hCut){
    //Set maximum eta value allowed for the V0 daughter tracks (|eta|<cut).
    fHistCutDaughterEta = hCut;
    fListOfPtDepCuts->Add(fHistCutDaughterEta);
    cout<<"[AliCascadeModule] Received pt-dep Daughter |eta| cut (max value)"<<endl;
}

void AliCascadeModule::SetCutCausality(TH1F* hCut){
    //Set maximum eta value allowed for the V0 daughter tracks (|eta|<cut).
    fHistCutCausality = hCut;
    fListOfPtDepCuts->Add(fHistCutCausality);
    cout<<"[AliCascadeModule] Received pt-dep causality cut (min value)"<<endl;
}

void AliCascadeModule::SetCutTOFExpTDiffNeg( TH1F* hCutMin, TH1F* hCutMax ) {
    //Set a pt-dep cut to limit the
    //time difference between expected and measured time-of-flight
    fHistCutMinTOFExpTDiffNeg = hCutMin;
    fHistCutMaxTOFExpTDiffNeg = hCutMax;
    fListOfPtDepCuts->Add(fHistCutMinTOFExpTDiffNeg);
    fListOfPtDepCuts->Add(fHistCutMaxTOFExpTDiffNeg);
    cout<<"[AliCascadeModule] Received pt-dep TOFExpTDiff cut for negative track (min, max values)"<<endl;
}
void AliCascadeModule::SetCutTOFExpTDiffPos( TH1F* hCutMin, TH1F* hCutMax ) {
    //Set a pt-dep cut to limit the
    //time difference between expected and measured time-of-flight
    fHistCutMinTOFExpTDiffPos = hCutMin;
    fHistCutMaxTOFExpTDiffPos = hCutMax;
    fListOfPtDepCuts->Add(fHistCutMinTOFExpTDiffPos);
    fListOfPtDepCuts->Add(fHistCutMaxTOFExpTDiffPos);
    cout<<"[AliCascadeModule] Received pt-dep TOFExpTDiff cut for positive track (min, max values)"<<endl;
}
void AliCascadeModule::SetCutTOFExpTDiffBach( TH1F* hCutMin, TH1F* hCutMax ) {
    //Set a pt-dep cut to limit the
    //time difference between expected and measured time-of-flight
    fHistCutMinTOFExpTDiffBach = hCutMin;
    fHistCutMaxTOFExpTDiffBach = hCutMax;
    fListOfPtDepCuts->Add(fHistCutMinTOFExpTDiffBach);
    fListOfPtDepCuts->Add(fHistCutMaxTOFExpTDiffBach);
    cout<<"[AliCascadeModule] Received pt-dep TOFExpTDiff cut for bachelor track (min, max values)"<<endl;
}

void AliCascadeModule::SetCutTOFSignalNeg( TH1F* hCutMin, TH1F* hCutMax ) {
    //Set a pt-dep cut to limit the measured time-of-flight
    fHistCutMinTOFSignalNeg = hCutMin;
    fHistCutMaxTOFSignalNeg = hCutMax;
    fListOfPtDepCuts->Add(fHistCutMinTOFSignalNeg);
    fListOfPtDepCuts->Add(fHistCutMaxTOFSignalNeg);
    cout<<"[AliCascadeModule] Received pt-dep TOF cut for negative track (min, max values)"<<endl;
}
void AliCascadeModule::SetCutTOFSignalPos( TH1F* hCutMin, TH1F* hCutMax ) {
    //Set a pt-dep cut to limit the measured time-of-flight
    fHistCutMinTOFSignalPos = hCutMin;
    fHistCutMaxTOFSignalPos = hCutMax;
    fListOfPtDepCuts->Add(fHistCutMinTOFSignalPos);
    fListOfPtDepCuts->Add(fHistCutMaxTOFSignalPos);
    cout<<"[AliCascadeModule] Received pt-dep TOF cut for positive track (min, max values)"<<endl;
}
void AliCascadeModule::SetCutTOFSignalBach( TH1F* hCutMin, TH1F* hCutMax ) {
    //Set a pt-dep cut to limit the measured time-of-flight
    fHistCutMinTOFSignalBach = hCutMin;
    fHistCutMaxTOFSignalBach = hCutMax;
    fListOfPtDepCuts->Add(fHistCutMinTOFSignalBach);
    fListOfPtDepCuts->Add(fHistCutMaxTOFSignalBach);
    cout<<"[AliCascadeModule] Received pt-dep TOF cut for bachelor track (min, max values)"<<endl;
}
///////////////////////////////////////////////////////////////////////////////

void AliCascadeModule::SetFitBackground ( Bool_t fitBgSwitch ){
    //Turns on background fitting for signal extraction instead of pure
    //bin counting. Useful, among other things, for systematics.
    fFitBackgroundSwitch = fitBgSwitch;
}

void AliCascadeModule::SetPeakPositionAndWidthFromFit(TF1* fMeanFit, TF1* fSigmaFit) {
    //Set user defined peak position for signal extraction
    fPeakPositionFit = fMeanFit;
    fPeakPositionFit->SetName("fPeakPositionFit");
    fPeakWidthFit = fSigmaFit;
    fPeakWidthFit->SetName("fPeakWidthFit");
    fUsePeakPositionAndWidthFromFit = kTRUE;
    cout<<"[AliCascadeModule] Switched to use provided functions for peak position and width"<<endl;
}

void AliCascadeModule::SetUsePeakPositionAndWidthFromFit( Bool_t flag ) {
    //Switch to use peak position and width from user provided function
    if( !(fPeakPositionFit && fPeakWidthFit) ) {
        cout<<"[AliCascadeModule] Please, provide parametrized functions for peak postion and width"<<endl;
        cout<<"                   using AliCascadeModule::SetPeakPositionAndWidthFromFit( fPosition, fWidth )"<<endl;
        fUsePeakPositionAndWidthFromFit = kFALSE;
        exit(-1);
    }
    else fUsePeakPositionAndWidthFromFit = flag;
}

void AliCascadeModule::SetUseMCBackgroundTemplate ( Bool_t fUseMCBackgroundTemplateReceived ){
    //Turns on MC template. Warning: takes precendence over Fitting!
    fUseMCBackgroundTemplate = fUseMCBackgroundTemplateReceived;
}

void AliCascadeModule::SetUseMinBiasMCBackgroundTemplate ( Bool_t fUseMinBiasMCBackgroundTemplateReceived ){
    //Turns on MC template. Warning: takes precendence over Fitting!
    fUseMinBiasMCBackgroundTemplate = fUseMinBiasMCBackgroundTemplateReceived;
}

void AliCascadeModule::SetPerformMultiplicityStudy ( Bool_t lPerformMultStudy ){
    //Turns on selection according to multiplicity of the event.
    //Boundaries are set in charged track multiplicity (pp) or in
    //centrality percentiles (PbPb). This is a requirement for studying
    //PbPb.
    fPerformMultiplicityStudy = lPerformMultStudy;
}

void AliCascadeModule::SetUseIntegratedEfficiency ( Bool_t lUseIntegratedEfficiency ){
    //To not select in multiplicity for efficiency
    fUseIntegratedEfficiency = lUseIntegratedEfficiency;
}

void AliCascadeModule::SetLowMultValue ( Double_t lLoMultBound       ){
    //Lower boundary (inclusive) in integer number for mult selection.
    //Note: If in PbPb and you want, say, 10-20% centrality, set this
    //to 10.
    fLoMultBound = lLoMultBound;
}

void AliCascadeModule::SetHighMultValue ( Double_t lHiMultBound       ){
    //Lower boundary (inclusive) in integer number for mult selection.
    //Note: If in PbPb and you want, say, 10-20% centrality, set this
    //to 10.
    fHiMultBound = lHiMultBound;
}

void AliCascadeModule::SetMultRange     ( Double_t lLoMultBound , Double_t lHiMultBound ){
    //Lower and Upper boundary (inclusive) in integer number for mult selection.
    //Sets both boundaries simultaneously for ease of use.
    fLoMultBound = lLoMultBound;
    fHiMultBound = lHiMultBound;
}

void AliCascadeModule::SetPerformResolutionTests ( Bool_t lPerformResolutionTests ){
    //Turns on resolution tests. Will consume a bit more time, but will store results
    //on resolution vs pT.
    fPerformResolutionTests = lPerformResolutionTests;
}

void AliCascadeModule::SetMVPileupRejection( Bool_t lMVPileupSwitch ) {
    //Reject events tagged with multiple vertices
    fMVPileupSwitch = lMVPileupSwitch;
}

void AliCascadeModule::SetMinDistToClosestNonEmptyBC( Int_t lNBCs ) {
    //Reject events with non empty neighboring BCs
    fMinDistToClosestNonEmptyBC = lNBCs;
}

void AliCascadeModule::SetCutTOFExpTDiffBach( Double_t min, Double_t max ) {
    //Set cut on the measured time-of-flight of the bachelor daughter
    fUseTOFmatchBach = kTRUE;
    fCutMinTOFExpTDiffBach = min;
    fCutMaxTOFExpTDiffBach = max;
    cout<<"[AliCascadeModule] Received TOF cut for bachelor track (min, max values) = ("<<min<<", "<<max<<")"<<endl;
}

void AliCascadeModule::SetCutTOFExpTDiffNeg( Double_t min, Double_t max ) {
    //Set cut on the measured time-of-flight of the negative daughter
    fUseTOFmatchNeg = kTRUE;
    fCutMinTOFExpTDiffNeg = min;
    fCutMaxTOFExpTDiffNeg = max;
    cout<<"[AliCascadeModule] Received TOF cut for negative track (min, max values) = ("<<min<<", "<<max<<")"<<endl;
}

void AliCascadeModule::SetCutTOFExpTDiffPos( Double_t min, Double_t max ) {
    //Set cut on the measured time-of-flight of the positive daughter
    fUseTOFmatchPos = kTRUE;
    fCutMinTOFExpTDiffPos = min;
    fCutMaxTOFExpTDiffPos = max;
    cout<<"[AliCascadeModule] Received TOF cut for positive track (min, max values) = ("<<min<<", "<<max<<")"<<endl;
}

void AliCascadeModule::SetCutTOFSignalNeg( Double_t min, Double_t max ) {
    //Set cut on the  measured time-of-flight of the negative daughter
    fUseTOFmatchNeg = kTRUE;
    fCutMinTOFSignalNeg = min;
    fCutMaxTOFSignalNeg = max;
    cout<<"[AliCascadeModule] Received TOF cut for negative track (min, max values) = ("<<min<<", "<<max<<")"<<endl;
}

void AliCascadeModule::SetCutTOFSignalPos( Double_t min, Double_t max ) {
    //Set cut on the measured time-of-flight of the positive daughter
    fUseTOFmatchPos = kTRUE;
    fCutMinTOFSignalPos = min;
    fCutMaxTOFSignalPos = max;
    cout<<"[AliCascadeModule] Received TOF cut for positive track (min, max values) = ("<<min<<", "<<max<<")"<<endl;
}

void AliCascadeModule::SetCutTOFSignalBach( Double_t min, Double_t max ) {
    //Set cut on the measured time-of-flight of the bachelor daughter
    fUseTOFmatchBach = kTRUE;
    fCutMinTOFSignalBach = min;
    fCutMaxTOFSignalBach = max;
    cout<<"[AliCascadeModule] Received TOF cut for bachelor track (min, max values) = ("<<min<<", "<<max<<")"<<endl;
}

void AliCascadeModule::SetSaveVarHistos( Bool_t lSaveSwitch ) {
    //Save histograms with selection variable distributions
    fSaveVarHistosSwitch = lSaveSwitch;
}

void AliCascadeModule::SetDefaultCuts(){
    //Sets Default cuts for analysis. (adjusted for adequate pp analysis)
    cout<<"[AliCascadeModule] Setting default cuts for particle species: "<<fWhichParticle<<endl;
    //Set Cuts - topological
    if ( fWhichParticle ==    "XiMinus" || fWhichParticle ==    "XiPlus" ) SetCutV0Radius(1.2); //1
    if ( fWhichParticle == "OmegaMinus" || fWhichParticle == "OmegaPlus" ) SetCutV0Radius(1.1); //1

    //These two depend on Species!
    if ( fWhichParticle == "XiMinus" || fWhichParticle == "OmegaMinus" ){
        //If particle then positive track is the baryon (proton) -> bends less
        SetCutDCANegToPV(0.04);
        SetCutDCAPosToPV(0.03);
    }

    if ( fWhichParticle == "XiPlus" || fWhichParticle == "OmegaPlus" ){
        //If anti-particle then positive track is the meson (pion) -> bends more
        SetCutDCANegToPV(0.03);
        SetCutDCAPosToPV(0.04);
    }

    SetCutDCABachToPV     (0.040); //4
    SetCutDCAV0Daughters  (  1.5); //5

    //Super Loose Cuts (check)
    SetCutV0CosPA         (0.97); //6
    SetCutCascCosPA       (0.97); //7

    SetCutDCAV0ToPV       ( 0.06); //8

    if ( fWhichParticle ==    "XiMinus" || fWhichParticle ==    "XiPlus" ) SetCutCascRadius(0.6); //1
    if ( fWhichParticle == "OmegaMinus" || fWhichParticle == "OmegaPlus" ) SetCutCascRadius(0.5); //1

    SetCutV0Mass          (0.008); //10
    SetCutDCACascDaughters(  1.3); //1

    SetCutDCAzNegToPV    (1.e+6); // do not apply by default
    SetCutDCAzPosToPV    (1.e+6); // do not apply by default
    SetCutDCAzBachToPV   (1.e+6); // do not apply by default

    if ( fWhichParticle == "XiMinus" || fWhichParticle == "XiPlus" ){
        SetCutProperLifetime    (3*4.91); // 3 times c*tau
        SetCutCompetingSpecies  (    -1); // Do not use.
    }
    if ( fWhichParticle == "OmegaMinus" || fWhichParticle == "OmegaPlus" ){
        SetCutProperLifetime    (3*2.46); // 3 times c*tau
        SetCutCompetingSpecies  ( 0.008); // 8 MeV/c^2 from Xi mass.
    }

    SetCutTPCPIDNSigmas                           (    4);
    SetCutSigmaForSignalExtraction                (    4);
    SetCutLeastNumberOfClusters                   (   70);
    SetCutMinTrackLength                          (  -1.);
    SetCutMaxChi2PerCluster                       ( 1.e9);
    SetCutDaughterEta                             (  0.8);

    SetCutDCABachToBaryon ( -1. );
    SetCutBBCosPA         ( 1.e6 );
}

// Selection Getters //////////////////////////////////////////////////////////
Double_t AliCascadeModule::GetCutV0Radius(Double_t pt){
    //Get minimum decay radius for the V0 in centimeters.
    //Note: this is (R_{2D}, cylindrical, centered around ALICE detector)
    if(fHistCutV0Radius) {
        return (Double_t)fHistCutV0Radius->GetBinContent(fHistCutV0Radius->FindBin(pt));
    }
    else return (Double_t)fCutV0Radius;
}
Double_t AliCascadeModule::GetCutCascRadius(Double_t pt){
    //Get minimum decay radius for the V0 in centimeters.
    //Note: this is (R_{2D}, cylindrical, centered around ALICE detector)
    if(fHistCutCascRadius) {
        return (Double_t)fHistCutCascRadius->GetBinContent(fHistCutCascRadius->FindBin(pt));
    }
    else return (Double_t)fCutCascRadius;
}
Double_t AliCascadeModule::GetCutDCANegToPV(Double_t pt){
    //Get minimum distance of closest approach between V0 negative daughter
    //track and primary vertex (in centimeters).
    if(fHistCutDCANegToPV) {
        return (Double_t)fHistCutDCANegToPV->GetBinContent(fHistCutDCANegToPV->FindBin(pt));
    }
    else return (Double_t)fCutDCANegToPV;
}
Double_t AliCascadeModule::GetCutDCAPosToPV(Double_t pt){
    //Get minimum distance of closest approach between V0 positive daughter
    //track and primary vertex (in centimeters).
    if(fHistCutDCAPosToPV) {
        return (Double_t)fHistCutDCAPosToPV->GetBinContent(fHistCutDCAPosToPV->FindBin(pt));
    }
    else return (Double_t)fCutDCAPosToPV;
}
Double_t AliCascadeModule::GetCutDCABachToPV(Double_t pt){
    //Get minimum distance of closest approach between V0 negative daughter
    //track and primary vertex (in centimeters).
    if(fHistCutDCABachToPV) {
        return (Double_t)fHistCutDCABachToPV->GetBinContent(fHistCutDCABachToPV->FindBin(pt));
    }
    else return (Double_t)fCutDCABachToPV;
}
Double_t AliCascadeModule::GetCutDCAV0Daughters(Double_t pt){
    //Get minimum distance of closest approach between V0 daughter
    //tracks, in sigmas. This is not in centimeters because if the
    //tracks have been determined with ITS refit the resolution is
    //greatly improved; thus, the cut may be tighter in that case.
    //Using sigmas will therefore take the tracking method in
    //consideration.
    if(fHistCutDCAV0Daughters) {
        return (Double_t)fHistCutDCAV0Daughters->GetBinContent(fHistCutDCAV0Daughters->FindBin(pt));
    }
    else return (Double_t)fCutDCAV0Daughters;
}
Double_t AliCascadeModule::GetCutDCACascDaughters(Double_t pt){
    //Get minimum distance of closest approach between Cascade daughters.
    if(fHistCutDCACascDaughters) {
        return (Double_t)fHistCutDCACascDaughters->GetBinContent(fHistCutDCACascDaughters->FindBin(pt));
    }
    else return (Double_t)fCutDCACascDaughters;
}
Double_t AliCascadeModule::GetCutV0CosPA(Double_t pt){
    //Get minimum value for the cosine of pointing angle of the V0.
    if(fHistCutV0CosPA) {
        return (Double_t)fHistCutV0CosPA->GetBinContent(fHistCutV0CosPA->FindBin(pt));
    }
    else return (Double_t)fCutV0CosPA;
}
Double_t AliCascadeModule::GetCutCascCosPA(Double_t pt){
    //Get minimum value for the cosine of pointing angle of the cascade.
    if(fHistCutCascCosPA) {
        return (Double_t)fHistCutCascCosPA->GetBinContent(fHistCutCascCosPA->FindBin(pt));
    }
    else return (Double_t)fCutCascCosPA;
}
Double_t AliCascadeModule::GetCutV0Mass(Double_t pt){
    //Get mass window width for Lambdas (in GeV/c^2)
    if(fHistCutV0Mass) {
        return (Double_t)fHistCutV0Mass->GetBinContent(fHistCutV0Mass->FindBin(pt));
    }
    else return (Double_t)fCutV0Mass;
}
Double_t AliCascadeModule::GetCutDCAV0ToPV(Double_t pt){
    //Get minimum distance of closest approach between V0 and PV
    if(fHistCutDCAV0ToPV) {
        return (Double_t)fHistCutDCAV0ToPV->GetBinContent(fHistCutDCAV0ToPV->FindBin(pt));
    }
    else return (Double_t)fCutDCAV0ToPV;
}
Double_t AliCascadeModule::GetCutDCACascToPV(Double_t pt){
    //Get maximum distance of closest approach between Cascade and PV
    if(fHistCutDCACascToPV) {
        return (Double_t)fHistCutDCACascToPV->GetBinContent(fHistCutDCACascToPV->FindBin(pt));
    }
    else return (Double_t)fCutDCACascToPV;
}
Double_t AliCascadeModule::GetCutDCAxyCascToPV(Double_t pt){
    //Get maximum distance of closest approach between Cascade and PV in the xy plane
    if(fHistCutDCAxyCascToPV) {
        return (Double_t)fHistCutDCAxyCascToPV->GetBinContent(fHistCutDCAxyCascToPV->FindBin(pt));
    }
    else return (Double_t)fCutDCAxyCascToPV;
}
Double_t AliCascadeModule::GetCutDCAzCascToPV(Double_t pt){
    //Get maximum distance of closest approach between Cascade and PV along the z direction
    if(fHistCutDCAzCascToPV) {
        return (Double_t)fHistCutDCAzCascToPV->GetBinContent(fHistCutDCAzCascToPV->FindBin(pt));
    }
    else return (Double_t)fCutDCAzCascToPV;
}
Double_t AliCascadeModule::GetCutDCAzNegToPV(Double_t pt) {
    //Get maximum distance of closest approach along z direction between V0 negative daughter
    //track and primary vertex (in centimeters).
    if(fHistCutDCAzNegToPV) {
        return (Double_t)fHistCutDCAzNegToPV->GetBinContent(fHistCutDCAzNegToPV->FindBin(pt));
    }
    else return (Double_t)fCutDCAzNegToPV;
}
Double_t AliCascadeModule::GetCutDCAzPosToPV(Double_t pt) {
    //Get maximum distance of closest approach along z direction between V0 positive daughter
    //track and primary vertex (in centimeters).
    if(fHistCutDCAzPosToPV) {
        return (Double_t)fHistCutDCAzPosToPV->GetBinContent(fHistCutDCAzPosToPV->FindBin(pt));
    }
    else return (Double_t)fCutDCAzPosToPV;
}
Double_t AliCascadeModule::GetCutDCAzBachToPV(Double_t pt){
    //Get maximum distance of closest approach along z direction between bachelor
    //track and primary vertex (in centimeters).
    if(fHistCutDCAzBachToPV) {
        return (Double_t)fHistCutDCAzBachToPV->GetBinContent(fHistCutDCAzBachToPV->FindBin(pt));
    }
    else return (Double_t)fCutDCAzBachToPV;
}
Double_t AliCascadeModule::GetCutDCABachToBaryon(Double_t pt){
    //Get minimum distance between bachelor and baryon daughter tracks
    if(fHistCutDCABachToBaryon) {
        return (Double_t)fHistCutDCABachToBaryon->GetBinContent(fHistCutDCABachToBaryon->FindBin(pt));
    }
    else return (Double_t)fCutDCABachToBaryon;
}
Double_t AliCascadeModule::GetCutBBCosPA(Double_t pt){
    //Get maximum value for the cosine of pointing angle from bachelor-baryon (wrong vertex)
    if(fHistCutBBCosPA) {
        return (Double_t)fHistCutBBCosPA->GetBinContent(fHistCutBBCosPA->FindBin(pt));
    }
    else return (Double_t)fCutBBCosPA;
}
///////////////////////////////////////////////////////////////////////////////


// Other Selection Getters ////////////////////////////////////////////////////
Double_t AliCascadeModule::GetCutProperLifetime(Double_t pt){
    //Get maximum value for m*L/p variable for the V0.
    //This is the "proper lifetime" selection and is usually called a
    //"c*tau cut". Should be set to a value larger than the c*tau for
    //the V0 considered.
    if(fHistCutProperLifetime) {
        return (Double_t)fHistCutProperLifetime->GetBinContent(fHistCutProperLifetime->FindBin(pt));
    }
    else return (Double_t)fCutProperLifetime;
}

Double_t AliCascadeModule::GetCutCompetingSpecies(Double_t pt){
    //Get minimum distance in xi invariant mass space for Omega Analysis.
    if(fHistCutCompetingSpecies) {
        return (Double_t)fHistCutCompetingSpecies->GetBinContent(fHistCutCompetingSpecies->FindBin(pt));
    }
    else return (Double_t)fCutCompetingSpecies;
}
Double_t AliCascadeModule::GetCutTPCPIDNSigmas(Double_t pt){
    //Get maximum deviation from the expected energy loss in
    //the TPC, in multiples of sigmas as computed from the AliPIDResponse
    //object. Selection is only used in real data and should thus
    //be very loose to ensure negligible signal loss.
    if(fHistCutTPCPIDNSigmas) {
        return (Double_t)fHistCutTPCPIDNSigmas->GetBinContent(fHistCutTPCPIDNSigmas->FindBin(pt));
    }
    else return (Double_t)fCutTPCPIDNSigmas;
}
Double_t AliCascadeModule::GetCutSigmaForSignalExtraction(Double_t pt){
    //Get number of sigmas for the signal extraction method. The value
    //provided is the number of sigmas from the peak position used for
    //the peak sampling, i.e. the peak will be from [<m>-cut*sigma,<m>+cut*sigma]
    //while the background sampling regions will be from
    //[<m>-2*cut*sigma,<m>+2*cut*sigma.
    if(fHistCutNSigmasForSignalExtraction) {
        return (Double_t)fHistCutNSigmasForSignalExtraction->GetBinContent(fHistCutNSigmasForSignalExtraction->FindBin(pt));
    }
    else return (Double_t)fCutNSigmasForSignalExtraction;
}
Double_t AliCascadeModule::GetCutLeastNumberOfClusters(Double_t pt){
    //Get smallest allowed number of TPC clusters for the V0 daughter tracks.
    if(fHistCutLeastNumberOfClusters) {
        return (Double_t)fHistCutLeastNumberOfClusters->GetBinContent(fHistCutLeastNumberOfClusters->FindBin(pt));
    }
    else return (Double_t)fCutLeastNumberOfClusters;
}
Double_t AliCascadeModule::GetCutMinTrackLength(Double_t pt) {
    //Get minimum track length allowed for the Cascade daughter tracks (|eta|<cut) for this pt.
    if(fHistCutMinTrackLength) {
        return (Double_t)fHistCutMinTrackLength->GetBinContent(fHistCutMinTrackLength->FindBin(pt));
    }
    else return (Double_t)fCutMinTrackLength;
}
Double_t AliCascadeModule::GetCutMaxChi2PerCluster(Double_t pt) {
    //Get maximum chi2/cluster for the Cascade daughter tracks for this pt.
    if(fHistCutMaxChi2PerCluster) {
        return (Double_t)fHistCutMaxChi2PerCluster->GetBinContent(fHistCutMaxChi2PerCluster->FindBin(pt));
    }
    else return (Double_t)fCutMaxChi2PerCluster;
}
Double_t AliCascadeModule::GetCutDaughterEta(Double_t pt){
    //Get maximum eta value allowed for the V0 daughter tracks (|eta|<cut).
    if(fHistCutDaughterEta) {
        return (Double_t)fHistCutDaughterEta->GetBinContent(fHistCutDaughterEta->FindBin(pt));
    }
    else return (Double_t)fCutDaughterEta;
}

Double_t AliCascadeModule::GetCutCausality(Double_t pt){
    //Get maximum eta value allowed for the V0 daughter tracks (|eta|<cut).
    if(fHistCutCausality) {
        return (Double_t)fHistCutCausality->GetBinContent(fHistCutCausality->FindBin(pt));
    }
    else return (Double_t)fCutCausality;
}

Double_t AliCascadeModule::GetCutMinTOFExpTDiffNeg(Double_t pt){
    //Get minimum tdiff value allowed for the negative daughter tracks.
    if(fHistCutMinTOFExpTDiffNeg) {
        return (Double_t)fHistCutMinTOFExpTDiffNeg->GetBinContent(fHistCutMinTOFExpTDiffNeg->FindBin(pt));
    }
    else return (Double_t)fCutMinTOFExpTDiffNeg;
}
Double_t AliCascadeModule::GetCutMaxTOFExpTDiffNeg(Double_t pt){
    //Get maximum tdiff value allowed for the negative daughter tracks.
    if(fHistCutMaxTOFExpTDiffNeg) {
        return (Double_t)fHistCutMaxTOFExpTDiffNeg->GetBinContent(fHistCutMaxTOFExpTDiffNeg->FindBin(pt));
    }
    else return (Double_t)fCutMaxTOFExpTDiffNeg;
}

Double_t AliCascadeModule::GetCutMinTOFExpTDiffPos(Double_t pt){
    //Get minimum tdiff value allowed for the positive daughter tracks.
    if(fHistCutMinTOFExpTDiffPos) {
        return (Double_t)fHistCutMinTOFExpTDiffPos->GetBinContent(fHistCutMinTOFExpTDiffPos->FindBin(pt));
    }
    else return (Double_t)fCutMinTOFExpTDiffPos;
}
Double_t AliCascadeModule::GetCutMaxTOFExpTDiffPos(Double_t pt){
    //Get maximum tdiff value allowed for the positive daughter tracks.
    if(fHistCutMaxTOFExpTDiffPos) {
        return (Double_t)fHistCutMaxTOFExpTDiffPos->GetBinContent(fHistCutMaxTOFExpTDiffPos->FindBin(pt));
    }
    else return (Double_t)fCutMaxTOFExpTDiffPos;
}

Double_t AliCascadeModule::GetCutMinTOFExpTDiffBach(Double_t pt){
    //Get minimum tdiff value allowed for the bachelor daughter tracks.
    if(fHistCutMinTOFExpTDiffBach) {
        return (Double_t)fHistCutMinTOFExpTDiffBach->GetBinContent(fHistCutMinTOFExpTDiffBach->FindBin(pt));
    }
    else return (Double_t)fCutMinTOFExpTDiffBach;
}
Double_t AliCascadeModule::GetCutMaxTOFExpTDiffBach(Double_t pt){
    //Get maximum tdiff value allowed for the bachelor daughter tracks.
    if(fHistCutMaxTOFExpTDiffBach) {
        return (Double_t)fHistCutMaxTOFExpTDiffBach->GetBinContent(fHistCutMaxTOFExpTDiffBach->FindBin(pt));
    }
    else return (Double_t)fCutMaxTOFExpTDiffBach;
}

Double_t AliCascadeModule::GetCutMinTOFSignalNeg(Double_t pt){
    //Get minimum tdiff value allowed for the negative daughter tracks.
    if(fHistCutMinTOFSignalNeg) {
        return (Double_t)fHistCutMinTOFSignalNeg->GetBinContent(fHistCutMinTOFSignalNeg->FindBin(pt));
    }
    else return (Double_t)fCutMinTOFSignalNeg;
}
Double_t AliCascadeModule::GetCutMaxTOFSignalNeg(Double_t pt){
    //Get maximum tdiff value allowed for the negative daughter tracks.
    if(fHistCutMaxTOFSignalNeg) {
        return (Double_t)fHistCutMaxTOFSignalNeg->GetBinContent(fHistCutMaxTOFSignalNeg->FindBin(pt));
    }
    else return (Double_t)fCutMaxTOFSignalNeg;
}

Double_t AliCascadeModule::GetCutMinTOFSignalPos(Double_t pt){
    //Get minimum tdiff value allowed for the positive daughter tracks.
    if(fHistCutMinTOFSignalPos) {
        return (Double_t)fHistCutMinTOFSignalPos->GetBinContent(fHistCutMinTOFSignalPos->FindBin(pt));
    }
    else return (Double_t)fCutMinTOFSignalPos;
}
Double_t AliCascadeModule::GetCutMaxTOFSignalPos(Double_t pt){
    //Get maximum tdiff value allowed for the positive daughter tracks.
    if(fHistCutMaxTOFSignalPos) {
        return (Double_t)fHistCutMaxTOFSignalPos->GetBinContent(fHistCutMaxTOFSignalPos->FindBin(pt));
    }
    else return (Double_t)fCutMaxTOFSignalPos;
}

Double_t AliCascadeModule::GetCutMinTOFSignalBach(Double_t pt){
    //Get minimum tdiff value allowed for the bachelor daughter tracks.
    if(fHistCutMinTOFSignalBach) {
        return (Double_t)fHistCutMinTOFSignalBach->GetBinContent(fHistCutMinTOFSignalBach->FindBin(pt));
    }
    else return (Double_t)fCutMinTOFSignalBach;
}
Double_t AliCascadeModule::GetCutMaxTOFSignalBach(Double_t pt){
    //Get maximum tdiff value allowed for the bachelor daughter tracks.
    if(fHistCutMaxTOFSignalBach) {
        return (Double_t)fHistCutMaxTOFSignalBach->GetBinContent(fHistCutMaxTOFSignalBach->FindBin(pt));
    }
    else return (Double_t)fCutMaxTOFSignalBach;
}

//
Bool_t AliCascadeModule::CheckTOFmatchBach(Double_t pt, Double_t tdiff) {
    //Check time difference between expected and measure tof for bachelor daughter
    if( !fUseTOFmatchBach ) return kTRUE; // don't care about TOF
    if( (TMath::Abs(tdiff+2500)<1.e-3) && (pt>fPtMinTOFmatchBach) && (pt<fPtMaxTOFmatchBach) ) return kTRUE; // no TOF info for this track
    Double_t min = GetCutMinTOFExpTDiffBach( pt );
    Double_t max = GetCutMaxTOFExpTDiffBach( pt );
    if(tdiff>min && tdiff<max) return kTRUE;
    else return kFALSE;
}
Bool_t AliCascadeModule::CheckTOFmatchNeg(Double_t pt, Double_t tdiff) {
    //Check time difference between expected and measure tof for bachelor daughter
    if( !fUseTOFmatchNeg ) return kTRUE; // don't care about TOF
    if( (TMath::Abs(tdiff+2500)<1.e-3) && (pt>fPtMinTOFmatchNeg) && (pt<fPtMaxTOFmatchNeg) ) return kTRUE; // no TOF info for this track
    Double_t min = GetCutMinTOFExpTDiffNeg( pt );
    Double_t max = GetCutMaxTOFExpTDiffNeg( pt );
    if(tdiff>min && tdiff<max) return kTRUE;
    else return kFALSE;
}
Bool_t AliCascadeModule::CheckTOFmatchPos(Double_t pt, Double_t tdiff) {
    //Check time difference between expected and measure tof for bachelor daughter
    if( !fUseTOFmatchPos ) return kTRUE; // don't care about TOF
    if( (TMath::Abs(tdiff+2500)<1.e-3) && (pt>fPtMinTOFmatchPos) && (pt<fPtMaxTOFmatchPos) ) return kTRUE; // no TOF info for this track
    Double_t min = GetCutMinTOFExpTDiffPos( pt );
    Double_t max = GetCutMaxTOFExpTDiffPos( pt );
    if(tdiff>min && tdiff<max) return kTRUE;
    else return kFALSE;
}
Bool_t AliCascadeModule::CheckTOFmatchOne(Double_t pt, Double_t bachtdiff, Double_t negtdiff, Double_t postdiff) {
    //Check time difference between expected and measure tof for one of the daughters
    if( !fUseTOFmatchOne || (pt<fPtMinTOFmatchOne) || (pt>fPtMaxTOFmatchOne) ) return kTRUE; // don't care about TOF
    Double_t bachmin = GetCutMinTOFExpTDiffBach( pt );
    Double_t bachmax = GetCutMaxTOFExpTDiffBach( pt );
    Double_t negmin  = GetCutMinTOFExpTDiffNeg ( pt );
    Double_t negmax  = GetCutMaxTOFExpTDiffNeg ( pt );
    Double_t posmin  = GetCutMinTOFExpTDiffPos ( pt );
    Double_t posmax  = GetCutMaxTOFExpTDiffPos ( pt );
    if( (bachtdiff>bachmin && bachtdiff<bachmax) || (negtdiff>negmin && negtdiff<negmax) || (postdiff>posmin && postdiff<posmax) ) return kTRUE;
    else return kFALSE;
}

//
Bool_t AliCascadeModule::CheckITSrefitBach(Double_t pt, ULong64_t trackstatus) {
    //Check ITS refit status for bachelor daughter
    if( !fUseITSrefitBach || (pt<fPtMinITSrefitBach) || (pt>fPtMaxITSrefitBach) ) return kTRUE; // don't care about ITS
    if( trackstatus & AliVTrack::kITSrefit ) return kTRUE; 
    return kFALSE;
}
Bool_t AliCascadeModule::CheckITSrefitNeg(Double_t pt, ULong64_t trackstatus) {
    //Check ITS refit status for negative daughter
    if( !fUseITSrefitNeg || (pt<fPtMinITSrefitNeg) || (pt>fPtMaxITSrefitNeg) ) return kTRUE; // don't care about ITS
    if( trackstatus & AliVTrack::kITSrefit ) return kTRUE; 
    return kFALSE;
}
Bool_t AliCascadeModule::CheckITSrefitPos(Double_t pt, ULong64_t trackstatus) {
    //Check ITS refit status for positive daughter
    if( !fUseITSrefitPos || (pt<fPtMinITSrefitPos) || (pt>fPtMaxITSrefitPos) ) return kTRUE; // don't care about ITS
    if( trackstatus & AliVTrack::kITSrefit ) return kTRUE; 
    return kFALSE;
}
Bool_t AliCascadeModule::CheckITSrefitOne(Double_t pt, ULong64_t btrackstatus, ULong64_t ntrackstatus, ULong64_t ptrackstatus) {
    //Check ITS refit status for one of the daughters 
    if( !fUseITSrefitOne || (pt<fPtMinITSrefitOne) || (pt>fPtMaxITSrefitOne) ) return kTRUE; // don't care about ITS
    if( (btrackstatus & AliVTrack::kITSrefit) | (ntrackstatus & AliVTrack::kITSrefit) | (ptrackstatus & AliVTrack::kITSrefit) ) return kTRUE; 
    return kFALSE;
}

//
Bool_t AliCascadeModule::CheckITSTOFBach(Double_t pt, ULong64_t btrackstatus, 
                                                     Double_t     bachtdiff) {
    //Check whether ITS or TOF is available for the bachelor track 
    if( !fUseITSTOFBach  || (pt<fPtMinITSTOFBach) || (pt>fPtMaxITSTOFBach) ) return kTRUE;

    Bool_t isITSrefitBach = kFALSE;
    if( (btrackstatus & AliVTrack::kITSrefit) ) isITSrefitBach = kTRUE; 

    Bool_t isTOFmatchBach = kFALSE;
    Double_t bachmin = GetCutMinTOFExpTDiffBach( pt );
    Double_t bachmax = GetCutMaxTOFExpTDiffBach( pt );
    if( (bachtdiff>bachmin && bachtdiff<bachmax) ) isTOFmatchBach = kTRUE;

    if( isITSrefitBach || isTOFmatchBach ) return kTRUE;
    else return kFALSE;
}
Bool_t AliCascadeModule::CheckITSTOFNeg(Double_t pt, ULong64_t ntrackstatus, 
                                                     Double_t     negtdiff) {
    //Check whether ITS or TOF is available for the negative track 
    if( !fUseITSTOFNeg  || (pt<fPtMinITSTOFNeg) || (pt>fPtMaxITSTOFNeg) ) return kTRUE;

    Bool_t isITSrefitNeg = kFALSE;
    if( (ntrackstatus & AliVTrack::kITSrefit) ) isITSrefitNeg = kTRUE; 

    Bool_t isTOFmatchNeg = kFALSE;
    Double_t negmin = GetCutMinTOFExpTDiffNeg( pt );
    Double_t negmax = GetCutMaxTOFExpTDiffNeg( pt );
    if( (negtdiff>negmin && negtdiff<negmax) ) isTOFmatchNeg = kTRUE;

    if( isITSrefitNeg || isTOFmatchNeg ) return kTRUE;
    else return kFALSE;
}
Bool_t AliCascadeModule::CheckITSTOFPos(Double_t pt, ULong64_t ptrackstatus, 
                                                     Double_t     postdiff) {
    //Check whether ITS or TOF is available for the positive track 
    if( !fUseITSTOFPos  || (pt<fPtMinITSTOFPos) || (pt>fPtMaxITSTOFPos) ) return kTRUE;

    Bool_t isITSrefitPos = kFALSE;
    if( (ptrackstatus & AliVTrack::kITSrefit) ) isITSrefitPos = kTRUE; 

    Bool_t isTOFmatchPos = kFALSE;
    Double_t posmin = GetCutMinTOFExpTDiffPos( pt );
    Double_t posmax = GetCutMaxTOFExpTDiffPos( pt );
    if( (postdiff>posmin && postdiff<posmax) ) isTOFmatchPos = kTRUE;

    if( isITSrefitPos || isTOFmatchPos ) return kTRUE;
    else return kFALSE;
}
Bool_t AliCascadeModule::CheckITSTOFOne(Double_t pt, ULong64_t btrackstatus, ULong64_t ntrackstatus, ULong64_t ptrackstatus, 
                                                     Double_t     bachtdiff, Double_t      negtdiff, Double_t      postdiff) {
    //Check whether ITS or TOF is available for one of the daughters 
    if( !fUseITSTOFOne  || (pt<fPtMinITSTOFOne) || (pt>fPtMaxITSTOFOne) ) return kTRUE;

    Bool_t isITSrefitOne = kFALSE;
    if( (btrackstatus & AliVTrack::kITSrefit) | (ntrackstatus & AliVTrack::kITSrefit) | (ptrackstatus & AliVTrack::kITSrefit) ) isITSrefitOne = kTRUE; 

    Bool_t isTOFmatchOne = kFALSE;
    Double_t bachmin = GetCutMinTOFExpTDiffBach( pt );
    Double_t bachmax = GetCutMaxTOFExpTDiffBach( pt );
    Double_t negmin  = GetCutMinTOFExpTDiffNeg ( pt );
    Double_t negmax  = GetCutMaxTOFExpTDiffNeg ( pt );
    Double_t posmin  = GetCutMinTOFExpTDiffPos ( pt );
    Double_t posmax  = GetCutMaxTOFExpTDiffPos ( pt );
    if( (bachtdiff>bachmin && bachtdiff<bachmax) || (negtdiff>negmin && negtdiff<negmax) || (postdiff>posmin && postdiff<posmax) ) isTOFmatchOne = kTRUE;

    if( isITSrefitOne || isTOFmatchOne ) return kTRUE;
    else return kFALSE;
}
///////////////////////////////////////////////////////////////////////////////

void AliCascadeModule::PrintConfiguration() {
    //Print current analysis configuration
    Double_t lParticleMass = 1.32171;
    if( fWhichParticle == "OmegaMinus" || fWhichParticle == "OmegaPlus" ) lParticleMass = 1.67245;
    cout<<"--------------- Configuration --------------------------"<<endl;
    cout<<" Analysed Particle.............: "<<fWhichParticle<<endl;
    cout<<" This Particle Mass............: "<<lParticleMass<<endl;
    cout<<" Rapidity Window, Lower Limit..: "<<fRapidityBoundaryLower<<endl;
    cout<<" Rapidity Window, Upper Limit..: "<<fRapidityBoundaryUpper<<endl;
    cout<<" Rapidity Type.................: y("<<fRapidityType<<")"<<endl;
    if( fRapidityType == "CMS")
        cout<<" Rapidity Shift................: "<<fRapidityShift<<endl;
    cout<<" CINT1B/INEL Ratio used........: "<<fCINT1BoverINELratio<<endl;
    if(fFuncGeantFlukaCorr) cout<<" Geant-Fluka correction........: "<<fFuncGeantFlukaCorr->GetName()<<endl;
    cout<<" -------------- Topological Selections -----------------"<<endl;
    cout<<" V0 Decay Radius...............: "<<Form("%s", fHistCutV0Radius         ? "pT-dep" : Form("%g", fCutV0Radius        ))<<endl; //1
    cout<<" Cascade Decay Radius..........: "<<Form("%s", fHistCutCascRadius       ? "pT-dep" : Form("%g", fCutCascRadius      ))<<endl; //2
    cout<<" DCA Negative track to PV......: "<<Form("%s", fHistCutDCANegToPV       ? "pT-dep" : Form("%g", fCutDCANegToPV      ))<<endl; //3
    cout<<" DCA Positive track to PV......: "<<Form("%s", fHistCutDCAPosToPV       ? "pT-dep" : Form("%g", fCutDCAPosToPV      ))<<endl; //4
    cout<<" DCA Bachelor track to PV......: "<<Form("%s", fHistCutDCABachToPV      ? "pT-dep" : Form("%g", fCutDCABachToPV     ))<<endl; //5
    cout<<" DCA V0 to PV..................: "<<Form("%s", fHistCutDCAV0ToPV        ? "pT-dep" : Form("%g", fCutDCAV0ToPV       ))<<endl; //6
    cout<<" DCA Cascade to PV.............: "<<Form("%s", fHistCutDCACascToPV      ? "pT-dep" : Form("%g", fCutDCACascToPV     ))<<endl;
    cout<<" DCA V0 Daughters..............: "<<Form("%s", fHistCutDCAV0Daughters   ? "pT-dep" : Form("%g", fCutDCAV0Daughters  ))<<endl; //7
    cout<<" DCA Cascade Daughters.........: "<<Form("%s", fHistCutDCACascDaughters ? "pT-dep" : Form("%g", fCutDCACascDaughters))<<endl; //8
    cout<<" Cosine of Pointing Angle V0...: "<<Form("%s", fHistCutV0CosPA          ? "pT-dep" : Form("%g", fCutV0CosPA         ))<<endl; //9
    cout<<" Cosine of Pointing Angle Casc.: "<<Form("%s", fHistCutCascCosPA        ? "pT-dep" : Form("%g", fCutCascCosPA       ))<<endl; //10
    cout<<" V0 Mass window (GeV/c^2)......: "<<Form("%s", fHistCutV0Mass           ? "pT-dep" : Form("%g", fCutV0Mass          ))<<endl; //11
    cout<<" DCAxy Cascade to PV...........: "<<Form("%s", fHistCutDCAxyCascToPV    ? "pT-dep" : Form("%g", fCutDCAxyCascToPV   ))<<endl;
    cout<<" DCAz Cascade to PV............: "<<Form("%s", fHistCutDCAzCascToPV     ? "pT-dep" : Form("%g", fCutDCAzCascToPV    ))<<endl;
    cout<<" DCAz Negative track to PV.....: "<<Form("%s", fHistCutDCAzNegToPV      ? "pT-dep" : Form("%g", fCutDCAzNegToPV     ))<<endl; //12
    cout<<" DCAz Positive track to PV.....: "<<Form("%s", fHistCutDCAzPosToPV      ? "pT-dep" : Form("%g", fCutDCAzPosToPV     ))<<endl; //13
    cout<<" DCAz Bachelor track to PV.....: "<<Form("%s", fHistCutDCAzBachToPV     ? "pT-dep" : Form("%g", fCutDCAzBachToPV    ))<<endl; //14
    cout<<" DCA Bachelor to Baryon........: "<<Form("%s", fHistCutDCABachToBaryon  ? "pT-dep" : Form("%g", fCutDCABachToBaryon ))<<endl; //15
    cout<<" Bachelor-Baryon CosPA.........: "<<Form("%s", fHistCutBBCosPA          ? "pT-dep" : Form("%g", fCutBBCosPA         ))<<endl; //16
    cout<<" -------------- Other Selections -----------------------"<<endl;
    cout<<" Proper Lifetime cut (cm)......: "<<Form("%s", fHistCutProperLifetime             ? "pT-dep" : Form("%g", fCutProperLifetime            ))<<endl;
    cout<<" TPC dE/dx sigmas cut (Real)...: "<<Form("%s", fHistCutTPCPIDNSigmas              ? "pT-dep" : Form("%g", fCutTPCPIDNSigmas             ))<<endl;
    cout<<" CutNSigmasForSignalExtraction.: "<<Form("%s", fHistCutNSigmasForSignalExtraction ? "pT-dep" : Form("%g", fCutNSigmasForSignalExtraction))<<endl;
    cout<<" Least # of Clusters...........: "<<Form("%s", fHistCutLeastNumberOfClusters      ? "pT-dep" : Form("%g", fCutLeastNumberOfClusters     ))<<endl;
    cout<<" Minimum track length..........: "<<Form("%s", fHistCutMinTrackLength             ? "pT-dep" : Form("%g", fCutMinTrackLength            ))<<endl;
    cout<<" Maximum chi2/cluster..........: "<<Form("%s", fHistCutMaxChi2PerCluster          ? "pT-dep" : Form("%g", fCutMaxChi2PerCluster         ))<<endl;
    cout<<" Daughter Track |eta| < .......: "<<Form("%s", fHistCutDaughterEta                ? "pT-dep" : Form("%g", fCutDaughterEta               ))<<endl;
    cout<<" Competing Species Rejection...: "<<Form("%s", fHistCutCompetingSpecies           ? "pT-dep" : Form("%g", fCutCompetingSpecies          ))<<endl;
    cout<<" Check MV Pileup...............: "<<Form("%s", (fMVPileupSwitch ? "YES" : "NO"))<<endl;
    cout<<" Min Distance to closest BC....: "<<Form("%d BCs", fMinDistToClosestNonEmptyBC)<<endl;
    cout<<" -------------- TOF match ------------------------------"<<endl;
    if(fUseTOFmatchOne || fUseTOFmatchNeg  || fUseITSTOFOne || fUseITSTOFNeg ) cout<<" TOFExpTDiff Negative Track....: "<<Form("%s", (fHistCutMinTOFExpTDiffNeg  && fHistCutMaxTOFExpTDiffNeg ) ? "pT-dep" : Form("(%g, %g) ns", fCutMinTOFExpTDiffNeg,  fCutMaxTOFExpTDiffNeg))<<endl;
    if(fUseTOFmatchOne || fUseTOFmatchPos  || fUseITSTOFOne || fUseITSTOFPos ) cout<<" TOFExpTDiff Positive Track....: "<<Form("%s", (fHistCutMinTOFExpTDiffPos  && fHistCutMaxTOFExpTDiffPos ) ? "pT-dep" : Form("(%g, %g) ns", fCutMinTOFExpTDiffPos,  fCutMaxTOFExpTDiffPos))<<endl;
    if(fUseTOFmatchOne || fUseTOFmatchBach || fUseITSTOFOne || fUseITSTOFBach) cout<<" TOFExpTDiff Bachelor Track....: "<<Form("%s", (fHistCutMinTOFExpTDiffBach && fHistCutMaxTOFExpTDiffBach) ? "pT-dep" : Form("(%g, %g) ns", fCutMinTOFExpTDiffBach, fCutMaxTOFExpTDiffBach))<<endl;
    if(fUseTOFmatchOne || fUseTOFmatchNeg  || fUseITSTOFOne || fUseITSTOFNeg ) cout<<" TOFSignal Negative Track......: "<<Form("%s", (fHistCutMinTOFSignalNeg  && fHistCutMaxTOFSignalNeg ) ? "pT-dep" : Form("(%g, %g) ns", fCutMinTOFSignalNeg,  fCutMaxTOFSignalNeg))<<endl;
    if(fUseTOFmatchOne || fUseTOFmatchPos  || fUseITSTOFOne || fUseITSTOFPos ) cout<<" TOFSignal Positive Track......: "<<Form("%s", (fHistCutMinTOFSignalPos  && fHistCutMaxTOFSignalPos ) ? "pT-dep" : Form("(%g, %g) ns", fCutMinTOFSignalPos,  fCutMaxTOFSignalPos))<<endl;
    if(fUseTOFmatchOne || fUseTOFmatchBach || fUseITSTOFOne || fUseITSTOFBach) cout<<" TOFSignal Bachelor Track......: "<<Form("%s", (fHistCutMinTOFSignalBach && fHistCutMaxTOFSignalBach) ? "pT-dep" : Form("(%g, %g) ns", fCutMinTOFSignalBach, fCutMaxTOFSignalBach))<<endl;
    if(fUseTOFmatchOne || fUseTOFmatchNeg  || fUseITSTOFOne || fUseITSTOFNeg ) cout<<" Neg.  Track pT range for match: ("<<fPtMinTOFmatchNeg<<", "<<fPtMaxTOFmatchNeg<<") GeV/c"<<endl;
    if(fUseTOFmatchOne || fUseTOFmatchPos  || fUseITSTOFOne || fUseITSTOFPos ) cout<<" Pos.  Track pT range for match: ("<<fPtMinTOFmatchPos<<", "<<fPtMaxTOFmatchPos<<") GeV/c"<<endl;
    if(fUseTOFmatchOne || fUseTOFmatchBach || fUseITSTOFOne || fUseITSTOFBach) cout<<" Bach. Track pT range for match: ("<<fPtMinTOFmatchBach<<", "<<fPtMaxTOFmatchBach<<") GeV/c"<<endl;
    cout<<" -------------- ITS refit ------------------------------"<<endl;
    if(fUseITSrefitOne || fUseITSrefitNeg  || fUseITSTOFOne || fUseITSTOFNeg ) cout<<" Neg.  Track pT range for refit: ("<<fPtMinITSrefitNeg<<", "<<fPtMaxITSrefitNeg<<") GeV/c"<<endl;
    if(fUseITSrefitOne || fUseITSrefitPos  || fUseITSTOFOne || fUseITSTOFPos ) cout<<" Pos.  Track pT range for refit: ("<<fPtMinITSrefitPos<<", "<<fPtMaxITSrefitPos<<") GeV/c"<<endl;
    if(fUseITSrefitOne || fUseITSrefitBach || fUseITSTOFOne || fUseITSTOFBach) cout<<" Bach. Track pT range for refit: ("<<fPtMinITSrefitBach<<", "<<fPtMaxITSrefitBach<<") GeV/c"<<endl;
    cout<<"--------------- File Names -----------------------------"<<endl;
    cout<<" Real Data File................: "<<fRealDataFile<<endl;
    cout<<" MC File.......................: "<<fMCDataFile<<endl;
    cout<<" Analysis output file..........: "<<fOutputDataFile<<endl;
    cout<<"--------------------------------------------------------"<<endl;
    cout<<endl;
}

TString AliCascadeModule::IntToString(int input){
    //Integer to TString Converter
    char dummyChar[50];
    sprintf(dummyChar, "%d", input);
    TString outputstring = dummyChar;
    return outputstring;
}

TString AliCascadeModule::DoubleToString(double input){
    //Double to TString Converter
    char dummyChar[50];
    sprintf(dummyChar, "%.3f", input);
    TString outputstring = dummyChar;
    return outputstring;
}

Double_t AliCascadeModule::ErrorInRatio ( Double_t A, Double_t Aerr, Double_t B, Double_t Berr ){
    //Error in a Ratio
    if(B!=0){
        Double_t errorfromtop = Aerr*Aerr / (B*B) ;
        Double_t errorfrombottom = ((A*A)/(B*B*B*B)) * Berr * Berr;
        return TMath::Sqrt( errorfromtop + errorfrombottom );
    }
    return 1;
}

Double_t AliCascadeModule::MyGeant3FlukaCorrectionForProtons(const Double_t *x, const Double_t *par){
    //Parametrization Used for Geant3/Fluka Correction for protons
    //Credit: Antonin Maire
    return 1 - par[0]*TMath::Exp(par[1]*x[0]) + par[2];
}


Double_t AliCascadeModule::MyGeant3FlukaCorrectionForAntiProtons(const Double_t *x, const Double_t *par){
    // Parametrization Used for Geant3/Fluka Correction for antiprotons
    // Credit: Antonin Maire
    //
    // La fonction A*TMath::Erf( B.x ) ne marche pas à bas pt.
    // Différentes fonctions ont été testées.
    // La fonction suivante semble donner de bons résultats
    // On peut jouer sur la puissance n du terme "1/x^n" pour repousser le pt auquel la fonction prend la valeur 1.0
    // Ici, pour n = 0.2, la fonction prend la valeur 0.9990 en pt = 10 GeV/c

    return 1 - par[0]*TMath::Exp(par[1]*x[0]) + par[2] + par[3]*1/TMath::Power(x[0], 0.2)*TMath::Log(x[0]);
}

Double_t AliCascadeModule::MyLevyPtXi(const Double_t *pt, const Double_t *par)
{
    //Levy Fit Function
    Double_t lMass  = 1.32171; //Xi Mass
    Double_t ldNdy  = par[0];
    Double_t lTemp = par[1];
    Double_t lPower = par[2];

    Double_t lBigCoef = ((lPower-1)*(lPower-2)) / (lPower*lTemp*(lPower*lTemp+lMass*(lPower-2)));
    Double_t lInPower = 1 + (TMath::Sqrt(pt[0]*pt[0]+lMass*lMass)-lMass) / (lPower*lTemp);

    return ldNdy * pt[0] * lBigCoef * TMath::Power(lInPower,(-1)*lPower);
    //return ldNdy * lBigCoef * TMath::Power(lInPower,(-1)*lPower);
}

Double_t AliCascadeModule::MyBgPol1(const Double_t *x, const Double_t *par)
{
    //Function for background fitting, rejects peak region
    if ( x[0] > par[2] && x[0] < par[3]) {
        TF1::RejectPoint();
        return 0;
    }
    return par[0] + par[1]*x[0];
}

Double_t AliCascadeModule::MyBgPolToEval1(const Double_t *x, const Double_t *par)
{
    //Just a plain linear function.
    return par[0] + par[1]*x[0];
}

Double_t AliCascadeModule::RoundToThousandth( const Double_t lToRound ){
    //Round any number to a hundredth...
    //Well, within machine precision...
    return TMath::Nint( 1000 * lToRound ) * 0.001;
}

void AliCascadeModule::DoAnalysis(){
    //---------------------------
    // Cascade Analysis Function
    //---------------------------
    //
    //
    //Consists of the following steps:
    //
    // (1) Loop over Real data candidates, acquire peak position and widths.
    // (2) Loop over Real data candidates, extract signal with variable extraction
    //     areas. Two loops: totally avoids binning data, allowing for sub-MeV/c^2
    //     granularity in signal extraction areas.
    // (3) Loop over MC data reconstructed candidates, find associated-to-MC primary
    //     candidates for efficiency numerator.
    // (4) Perform Geant3/Fluka correction.
    // (5) Get generated primary Cascade histograms from MC file for efficiency denominator.
    // (6) Perform detection efficiency correction and compute final corrected spectra.
    //
    //
    // Normalization:
    //  --- Number of Inelastic Events estimated to be:
    //      <Number of CINT1B triggers> / <CINT1B/INEL ratio>
    //
    //  --- Efficiency denominator: filled at pre-physics selection stage. All signal
    //      losses are therefore estimated from Monte Carlo.
    //
    //  --- Pileup: not included in MC. Fraction of events removed by IsPileupFromSPD
    //      is artificially removed from Number of Inelastic Events as well.
    //
    // Output: written to specified file (with SetOutputFile(...)).
    //   --- includes spectra and efficiencies in base directory.
    //   --- includes a number of different subdirectories storing plots such as
    //       such as invariant mass histos, signal extraction, pt resolution, Geant-Fluka
    //       related distributions, feeddown related distributions, and so on.

    //Set Batch Mode: Ignore all canvases in this section
    gROOT->SetBatch (kTRUE);

    cout<<"======================================"<<endl;
    cout<<" -------- Spectra Extraction -------- "<<endl;
    cout<<"======================================"<<endl;
    cout<<endl;
    if(fptbinnumb == -1){
        cout<<"[AliCascadeModule] It's unclear what kind of pt binning you want."<<endl;
        cout<<"[AliCascadeModule] Most likely you forgot to set it using SetPtBinLimits..."<<endl;
        cout<<"[AliCascadeModule] Analysis will NOT be done. Returning."<<endl;
        return;
    }
    //Mass: Needed for decay length
    Double_t lParticleMass = 1.32171;
    if( fWhichParticle == "OmegaMinus" || fWhichParticle == "OmegaPlus" ) lParticleMass = 1.67245;

    // Set variable names and corresponding enum for consistent usage in the code
    enum                 {  V0RADIUS ,  CASCRADIUS ,  DCANEGTOPV ,  DCAPOSTOPV ,  DCABACHTOPV ,  DCAV0TOPV ,  DCACASCTOPV ,  DCAV0DAUGHTERS ,  DCACASCDAUGHTERS ,  DCAXYCASCTOPV ,  DCAZCASCTOPV ,  DCAZNEGTOPV ,  DCAZPOSTOPV ,  DCAZBACHTOPV ,  DCABACHTOBARYON ,  V0PA ,  CASCPA ,  BBPA ,  V0MASS ,  PROPERLIFETIME ,  COMPETINGSPECIES ,  TPCNEGDEDX ,  TPCPOSDEDX ,  TPCBACHDEDX ,  TPCNEGPIDNSIGMAS ,  TPCPOSPIDNSIGMAS ,    TPCBACHPIDNSIGMAS ,  TPCNCLUSTERS ,  MAXCHI2PERCLUSTER ,  MINTRACKLENGTH ,  NEGDAUGHTERETA ,  POSDAUGHTERETA ,  BACHDAUGHTERETA ,  NEGINNERP ,  POSINNERP ,  BACHINNERP , NEGTOFEXPTDIFF ,  POSTOFEXPTDIFF ,  BACHTOFEXPTDIFF , NEGTOFSIGNAL ,  POSTOFSIGNAL ,  BACHTOFSIGNAL  ,  NVARS };
    TString lVarName[] = { "V0Radius", "CascRadius", "DCANegToPV", "DCAPosToPV", "DCABachToPV", "DCAV0ToPV", "DCACascToPV", "DCAV0Daughters", "DCACascDaughters", "DCAxyCascToPV", "DCAzCascToPV", "DCAzNegToPV", "DCAzPosToPV", "DCAzBachToPV", "DCABachToBaryon", "V0PA", "CascPA", "BBPA", "V0Mass", "ProperLifetime", "CompetingSpecies", "TPCNegdEdx", "TPCPosdEdx", "TPCBachdEdx", "TPCNegPIDNSigmas", "TPCPosPIDNSigmas", "TPCBachPIDNSigmas", "TPCNClusters", "MaxChi2PerCluster", "MinTrackLength", "NegDaughterEta", "PosDaughterEta", "BachDaughterEta", "NegInnerP", "PosInnerP", "BachInnerP", "NegTOFExpTDiff", "PosTOFExpTDiff", "BachTOFExpTDiff", "NegTOFSignal", "PosTOFSignal", "BachTOFSignal" };

    PrintConfiguration();

    //save current configuration
    TH1F* fHistConfig = new TH1F("fHistConfig","Analysis configuration", 50, 0, 50);
    //
    // General info
    fHistConfig->Fill( "ParticleMass"                       , lParticleMass);
    fHistConfig->Fill( "RapidityBoundaryLower"              , fRapidityBoundaryLower);
    fHistConfig->Fill( "RapidityBoundaryUpper"              , fRapidityBoundaryUpper);
    fHistConfig->Fill( "RapidityShift"                      , fRapidityShift);
    fHistConfig->Fill( lVarName[ NEGDAUGHTERETA     ].Data(), fCutDaughterEta);
    fHistConfig->Fill( lVarName[ POSDAUGHTERETA     ].Data(), fCutDaughterEta);
    fHistConfig->Fill( lVarName[ BACHDAUGHTERETA    ].Data(), fCutDaughterEta);
    //
    // Event selections
    fHistConfig->Fill( "MVPileupRejection"                  , fMVPileupSwitch);
    fHistConfig->Fill( "MinDistToClosestNonEmptyBC"         , fMinDistToClosestNonEmptyBC);
    //
    // Topological cuts
    fHistConfig->Fill( lVarName[ V0RADIUS           ].Data(), fCutV0Radius);
    fHistConfig->Fill( lVarName[ CASCRADIUS         ].Data(), fCutCascRadius);
    fHistConfig->Fill( lVarName[ V0PA               ].Data(), TMath::ACos(fCutV0CosPA));
    fHistConfig->Fill( lVarName[ CASCPA             ].Data(), TMath::ACos(fCutCascCosPA));
    fHistConfig->Fill( lVarName[ DCANEGTOPV         ].Data(), fCutDCANegToPV);
    fHistConfig->Fill( lVarName[ DCAPOSTOPV         ].Data(), fCutDCAPosToPV);
    fHistConfig->Fill( lVarName[ DCABACHTOPV        ].Data(), fCutDCABachToPV);
    fHistConfig->Fill( lVarName[ DCAV0TOPV          ].Data(), fCutDCAV0ToPV);
    fHistConfig->Fill( lVarName[ DCACASCTOPV        ].Data(), fCutDCACascToPV);
    fHistConfig->Fill( lVarName[ DCAV0DAUGHTERS     ].Data(), fCutDCAV0Daughters);
    fHistConfig->Fill( lVarName[ DCACASCDAUGHTERS   ].Data(), fCutDCACascDaughters);
    fHistConfig->Fill( lVarName[ DCAXYCASCTOPV      ].Data(), fCutDCAxyCascToPV);
    fHistConfig->Fill( lVarName[ DCAZCASCTOPV       ].Data(), fCutDCAzCascToPV);
    fHistConfig->Fill( lVarName[ DCAZNEGTOPV        ].Data(), fCutDCAzNegToPV);
    fHistConfig->Fill( lVarName[ DCAZPOSTOPV        ].Data(), fCutDCAzPosToPV);
    fHistConfig->Fill( lVarName[ DCAZBACHTOPV       ].Data(), fCutDCAzBachToPV);
    fHistConfig->Fill( lVarName[ V0MASS             ].Data(), fCutV0Mass);
    fHistConfig->Fill( lVarName[ DCABACHTOBARYON    ].Data(), fCutDCABachToBaryon);
    fHistConfig->Fill( lVarName[ BBPA               ].Data(), TMath::ACos(fCutBBCosPA));
    //
    // Other cuts
    fHistConfig->Fill( lVarName[ PROPERLIFETIME     ].Data(), fCutProperLifetime);
    fHistConfig->Fill( lVarName[ TPCNEGPIDNSIGMAS   ].Data(), fCutTPCPIDNSigmas);
    fHistConfig->Fill( lVarName[ TPCPOSPIDNSIGMAS   ].Data(), fCutTPCPIDNSigmas);
    fHistConfig->Fill( lVarName[ TPCBACHPIDNSIGMAS  ].Data(), fCutTPCPIDNSigmas);
    fHistConfig->Fill( lVarName[ TPCNCLUSTERS       ].Data(), fCutLeastNumberOfClusters);
    fHistConfig->Fill( lVarName[ MINTRACKLENGTH     ].Data(), fCutMinTrackLength);
    fHistConfig->Fill( lVarName[ MAXCHI2PERCLUSTER  ].Data(), fCutMaxChi2PerCluster);
    fHistConfig->Fill( lVarName[ COMPETINGSPECIES   ].Data(), fCutCompetingSpecies);
    fHistConfig->Fill( lVarName[ NEGINNERP          ].Data(), fCutNegInnerP);
    fHistConfig->Fill( lVarName[ POSINNERP          ].Data(), fCutPosInnerP);
    fHistConfig->Fill( lVarName[ BACHINNERP         ].Data(), fCutBachInnerP);
    fHistConfig->Fill( "MinNegTOFExpTDiff"                  , fCutMinTOFExpTDiffNeg);
    fHistConfig->Fill( "MaxNegTOFExpTDiff"                  , fCutMaxTOFExpTDiffNeg);
    fHistConfig->Fill( "MinPosTOFExpTDiff"                  , fCutMinTOFExpTDiffPos);
    fHistConfig->Fill( "MaxPosTOFExpTDiff"                  , fCutMaxTOFExpTDiffPos);
    fHistConfig->Fill( "MinBachTOFExpTDiff"                 , fCutMinTOFExpTDiffBach);
    fHistConfig->Fill( "MaxBachTOFExpTDiff"                 , fCutMaxTOFExpTDiffBach);
    fHistConfig->Fill( "MinNegTOFSignal"                    , fCutMinTOFSignalNeg);
    fHistConfig->Fill( "MaxNegTOFSignal"                    , fCutMaxTOFSignalNeg);
    fHistConfig->Fill( "MinPosTOFSignal"                    , fCutMinTOFSignalPos);
    fHistConfig->Fill( "MaxPosTOFSignal"                    , fCutMaxTOFSignalPos);
    fHistConfig->Fill( "MinBachTOFSignal"                   , fCutMinTOFSignalBach);
    fHistConfig->Fill( "MaxBachTOFSignal"                   , fCutMaxTOFSignalBach);
    fHistConfig->Fill( "SigmasForSignalExtraction"          , fCutNSigmasForSignalExtraction);

    // save number of events processed
    TH1F* fHistEventCounter = new TH1F("fHistEventCounter","Number of Events Analysed", 2, 0, 2);

    //=== Real data loop 1: Acquire peak positions, widths ===
    //--- Preparing... ---
    //defining helping histogram - only used to FindBin Index!========
    TH1F* fHistPt          = new TH1F("fHistPt"         , "Dummy;p_{T} (GeV/c);Counts", fptbinnumb, fptbinlimits);
    TH1F* fHistPtRaw       = new TH1F("fHistPtRaw"      , "Dummy;p_{T} (GeV/c);Counts", fptbinnumb, fptbinlimits);
    TH1F* fHistPtRawMC	   = new TH1F("fHistPtRawMC"    , "Dummy;p_{T} (GeV/c);Counts", fptbinnumb, fptbinlimits);
    TH1F* fHistPtSignal    = new TH1F("fHistPtSignal"   , "Dummy;p_{T} (GeV/c);Counts", fptbinnumb, fptbinlimits);
    TH1F* fHistPtSignalMC  = new TH1F("fHistPtSignalMC" , "Dummy;p_{T} (GeV/c);Counts", fptbinnumb, fptbinlimits);
    TH1F* fHistPtGenerated = new TH1F("fHistPtGenerated", "Dummy;p_{T} (GeV/c);Counts", 200, 0, 20);
    TH1F* fHistGenPerEvent = new TH1F("fHistGenPerEvent", "Dummy;p_{T} (GeV/c);Counts",   1, 0,  1);
    //Peak Position Histograms
    TH1F* fHistPeakPosition	  = new TH1F("fHistPeakPosition"  , "Peak Position (Real);p_{T} (GeV/c);Peak Position"   , fptbinnumb, fptbinlimits);
    TH1F* fHistPeakPositionMC = new TH1F("fHistPeakPositionMC", "Peak Position (MC);p_{T} (GeV/c);Peak Position"     , fptbinnumb, fptbinlimits);
    //Signal to Noise Histograms
    TH1F* fHistSigToNoise     = new TH1F("fHistSigToNoise"    , "Signal to Noise Ratio (real);p_{T} (GeV/c);Sig / Bg", fptbinnumb, fptbinlimits);
    TH1F* fHistSigToNoiseMC	  = new TH1F("fHistSigToNoiseMC"  , "Signal to Noise Ratio (MC);p_{T} (GeV/c);Sig / Bg"  , fptbinnumb, fptbinlimits);

    //Signal Extraction Range Histogram
    TH1F* fHistSignalExtractionRange = new TH1F("fHistSignalExtractionRange","Sig. Ext. Range;p_{T} (GeV/c);Range", fptbinnumb, fptbinlimits);

    //Resolution Histogram (filled with MC)
    TH2F* f2dHistPtResolution = new TH2F("f2dHistPtResolution", "p_{t} Resolution;p_{t} (reco);p_{t} (mc)"       , fptbinnumb, fptbinlimits, fptbinnumb, fptbinlimits);
    TH2F* f2dHistPtBlur       = new TH2F("f2dHistPtBlur"      , "p_{t} blurring matrix;p_{t} (reco);p_{t} (mc)"  , fptbinnumb, fptbinlimits, fptbinnumb, fptbinlimits);
    TH2F* f2dHistPtSharpen    = new TH2F("f2dHistPtSharpen"   , "p_{t} sharpening matrix;p_{t} (reco);p_{t} (mc)", fptbinnumb, fptbinlimits, fptbinnumb, fptbinlimits);

    //Efficiency Denominator + Numerator
    TH1F* fHistEffNumerator   = new TH1F("fHistEffNumerator"  , "Efficiency Numerator;p_{T} (GeV/c);Reconst.+Assoc.", fptbinnumb, fptbinlimits);
    TH1F* fHistEffDenominator = new TH1F("fHistEffDenominator", "Efficiency Denominator;p_{T} (GeV/c);Generated"    , fptbinnumb, fptbinlimits);

    /********************************************************

      ---> Setting up histogram boundaries

     ********************************************************/

    Int_t lWeAreAtBin = 0;
    Double_t lHistoLowerBoundary = -1;
    Double_t lHistoUpperBoundary = -1;
    Long_t lHistoNBins = -1;
    if ( fWhichParticle == "XiMinus" || fWhichParticle == "XiPlus" ){
        lHistoLowerBoundary = 1.322 - 0.1;
        lHistoUpperBoundary = 1.322 + 0.1;
        lHistoNBins = 200; //1MeV/c^2 binning
    }
    if ( fWhichParticle == "OmegaMinus" || fWhichParticle == "OmegaPlus" ){
        lHistoLowerBoundary = 1.672 - 0.1;
        lHistoUpperBoundary = 1.672 + 0.1;
        lHistoNBins = 200; //1MeV/c^2 binning
    }
    //Setting Up Base Histograms to use===============================
    TH1F* lHistoMBCasc[100];
    TH1F* lHistoSelectedCascMC[100];
    TH1F* lHistoBgTemplateMC[100];

    TCanvas* lCanvasHistoMBCasc[100];
    TCanvas* lCanvasHistoSelectedCascMC[100];
    //Selection-Level (as opposed to Min-Bias)
    TH1F* lHistoSelectedCasc[100];
    TCanvas* lCanvasHistoSelectedCasc[100];
    //TH1F* lHistoSelectedCascMC[100]; //OPTIONAL


    TH1F *lHistResolution[100];
    TF1 *lHistResolutionGaussian[100];
    char lHistResolutionName[100];
    TH1F* fHistResolutionVsPt 		= new TH1F("fHistResolutionVsPt","Resolution vs p_{t};p_{t} (GeV/c);#Delta(p_{t}) = p_{t}-p_{t}^{true} (GeV/c)",fptbinnumb,fptbinlimits);
    TH1F* fHistResolutionVsPtDivByBinWidth 		= new TH1F("fHistResolutionVsPtDivByBinWidth","Resolution vs p_{t};p_{t} (GeV/c);#Delta(p_{t})/#Delta^{bin}(p_{t}) = (p_{t}-p_{t}^{true})/#Delta^{bin}(p_{t}) (GeV/c)",fptbinnumb,fptbinlimits);
    TH1F* fHistResolutionVsPtWithGaussians 		= new TH1F("fHistResolutionVsPtWithGaussians","Resolution vs p_{t};p_{t} (GeV/c);#Delta(p_{t}) = p_{t}-p_{t}^{true} (GeV/c)",fptbinnumb,fptbinlimits);
    TH1F* fHistResolutionVsPtDivByBinWidthWithGaussians 		= new TH1F("fHistResolutionVsPtDivByBinWidthWithGaussians","Resolution vs p_{t};p_{t} (GeV/c);#Delta(p_{t})/#Delta^{bin}(p_{t}) = (p_{t}-p_{t}^{true})/#Delta^{bin}(p_{t}) (GeV/c)",fptbinnumb,fptbinlimits);
    char histname[80]; 	TString bindescription = "";
    for(Int_t ihist=0; ihist<100; ++ihist) {
        //__________________________________________________________________________________
        //Histo For Real Data
        sprintf(histname,"lHistoMBCasc%i",ihist);
        if(fWhichParticle == "XiMinus")     bindescription="#Xi^{-}, bin #";
        if(fWhichParticle == "XiPlus")      bindescription="#bar{#Xi}^{+}, bin #";
        if(fWhichParticle == "OmegaMinus")  bindescription="#Omega^{-}, bin #";
        if(fWhichParticle == "OmegaPlus")   bindescription="#bar{#Omega}^{+}, bin #";
        bindescription.Append(IntToString( ihist ));
        if ( ihist < fptbinnumb ){
            bindescription.Append(", ");
            bindescription.Append(DoubleToString(fptbinlimits[ ihist ]));
            bindescription.Append("-");
            bindescription.Append(DoubleToString(fptbinlimits[ ihist+1 ]));
            bindescription.Append("GeV/c");
        }
        lHistoMBCasc[ihist]   =	new TH1F(histname,"Candidates;M(appropriate) (GeV/c^{2});Counts", lHistoNBins,lHistoLowerBoundary,lHistoUpperBoundary);
        lHistoMBCasc[ihist]->SetTitle(bindescription);
        sprintf(histname,"lCanvasHistoMBCasc%i",ihist);
        lCanvasHistoMBCasc[ihist] = new TCanvas(histname, bindescription, 800, 600);
        lCanvasHistoMBCasc[ihist] -> SetFillColor(kWhite);
        lCanvasHistoMBCasc[ihist] -> SetLeftMargin( 0.175 );
        lCanvasHistoMBCasc[ihist] -> SetBottomMargin( 0.175 );

        //__________________________________________________________________________________
        //Histo for Real Data at Selection Level
        sprintf(histname,"lHistoSelectedCasc%i",ihist);
        if(fWhichParticle == "XiMinus")     bindescription="#Xi^{-}, bin #";
        if(fWhichParticle == "XiPlus")      bindescription="#bar{#Xi}^{+}, bin #";
        if(fWhichParticle == "OmegaMinus")  bindescription="#Omega^{-}, bin #";
        if(fWhichParticle == "OmegaPlus")   bindescription="#bar{#Omega}^{+}, bin #";
        bindescription.Append(IntToString( ihist ));
        if ( ihist < fptbinnumb ){
            bindescription.Append(", ");
            bindescription.Append(DoubleToString(fptbinlimits[ ihist ]));
            bindescription.Append("-");
            bindescription.Append(DoubleToString(fptbinlimits[ ihist+1 ]));
            bindescription.Append("GeV/c");
        }
        lHistoSelectedCasc[ihist]   =	new TH1F(histname,"Candidates;M(appropriate) (GeV/c^{2});Counts", lHistoNBins,lHistoLowerBoundary,lHistoUpperBoundary);
        lHistoSelectedCasc[ihist]->SetTitle(bindescription);
        sprintf(histname,"lCanvasHistoSelectedCasc%i",ihist);
        lCanvasHistoSelectedCasc[ihist] = new TCanvas(histname, bindescription, 800, 600);
        lCanvasHistoSelectedCasc[ihist] -> SetFillColor(kWhite);
        lCanvasHistoSelectedCasc[ihist] -> SetLeftMargin( 0.175 );
        lCanvasHistoSelectedCasc[ihist] -> SetBottomMargin( 0.175 );

        //__________________________________________________________________________________
        //Histo for MC
        sprintf(histname,"lHistoSelectedCascMC%i",ihist);
        if(fWhichParticle == "XiMinus")     bindescription="#Xi^{-}, bin #";
        if(fWhichParticle == "XiPlus")      bindescription="#bar{#Xi}^{+}, bin #";
        if(fWhichParticle == "OmegaMinus")  bindescription="#Omega^{-}, bin #";
        if(fWhichParticle == "OmegaPlus")   bindescription="#bar{#Omega}^{+}, bin #";
        bindescription.Append(IntToString( ihist ));
        if ( ihist < fptbinnumb ){
            bindescription.Append(", ");
            bindescription.Append(DoubleToString(fptbinlimits[ ihist ]));
            bindescription.Append("-");
            bindescription.Append(DoubleToString(fptbinlimits[ ihist+1 ]));
            bindescription.Append("GeV/c");
        }
        lHistoSelectedCascMC[ihist]   =	new TH1F(histname,"Candidates;M(appropriate) (GeV/c^{2});Counts", lHistoNBins,lHistoLowerBoundary,lHistoUpperBoundary);
        lHistoSelectedCascMC[ihist]->SetTitle(bindescription);
        sprintf(histname,"lCanvasHistoSelectedCascMC%i",ihist);
        lCanvasHistoSelectedCascMC[ihist] = new TCanvas(histname, bindescription, 800, 600);
        lCanvasHistoSelectedCascMC[ihist] -> SetFillColor(kWhite);
        lCanvasHistoSelectedCascMC[ihist] -> SetLeftMargin( 0.175 );
        lCanvasHistoSelectedCascMC[ihist] -> SetBottomMargin( 0.175 );

        //__________________________________________________________________________________
        //Histo for MC at Selection Level - Will be used for MC template method
        //lHistoBgTemplateMC
        sprintf(histname,"lHistoBgTemplateMC%i",ihist);
        if(fWhichParticle == "XiMinus")     bindescription="#Xi^{-}, bin #";
        if(fWhichParticle == "XiPlus")      bindescription="#bar{#Xi}^{+}, bin #";
        if(fWhichParticle == "OmegaMinus")  bindescription="#Omega^{-}, bin #";
        if(fWhichParticle == "OmegaPlus")   bindescription="#bar{#Omega}^{+}, bin #";
        bindescription.Append(IntToString( ihist ));
        if ( ihist < fptbinnumb ){
            bindescription.Append(", ");
            bindescription.Append(DoubleToString(fptbinlimits[ ihist ]));
            bindescription.Append("-");
            bindescription.Append(DoubleToString(fptbinlimits[ ihist+1 ]));
            bindescription.Append("GeV/c");
        }
        lHistoBgTemplateMC[ihist]   =	new TH1F(histname,"Candidates;M(appropriate) (GeV/c^{2});Counts", lHistoNBins,lHistoLowerBoundary,lHistoUpperBoundary);
        lHistoBgTemplateMC[ihist]->SetTitle(bindescription);

        //__________________________________________________________________________________
        //Histo for resolution tests
        sprintf(lHistResolutionName,"lHistResolution%i",ihist);
        Long_t lNumberOfBinsResolution = 1000;
        if ( ihist < fptbinnumb ) {
            if ( (fptbinlimits[ihist+1]+fptbinlimits[ihist])*0.5 > 5 )  lNumberOfBinsResolution = 100;
            if ( (fptbinlimits[ihist+1]+fptbinlimits[ihist])*0.5 > 10 ) lNumberOfBinsResolution = 10;
        }

        //cout << ihist << endl;//" ---> " << lNumberOfBinsResolution << endl;
        lHistResolution[ihist] = new TH1F ( lHistResolutionName, bindescription, 1000, -5., 5.); //histo with 5MeV/c precision!
        //lHistResolution[ihist] = new TH1F ( lHistResolutionName, bindescription, lNumberOfBinsResolution, -5., 5.); //histo with 5MeV/c precision!
        sprintf(lHistResolutionName, "lHistResolutionGaussian%i", ihist);
        lHistResolutionGaussian[ihist] = new TF1(lHistResolutionName, "[0]*TMath::Gaus(x,[1],[2])",-5,5);

    }
    //================================================================
    cout<<endl;

    // save distributions for each variable in the ttree (signal/bkg) -- rafael
    Int_t lVarNbins[NVARS];                     Float_t lVarBinLo[NVARS];                     Float_t lVarBinHi[NVARS];
    lVarNbins[ V0RADIUS          ] =   400;     lVarBinLo[ V0RADIUS          ] =     0.0;     lVarBinHi[ V0RADIUS          ] =    40.0;
    lVarNbins[ CASCRADIUS        ] =   400;     lVarBinLo[ CASCRADIUS        ] =     0.0;     lVarBinHi[ CASCRADIUS        ] =    40.0;
    lVarNbins[ DCANEGTOPV        ] =   500;     lVarBinLo[ DCANEGTOPV        ] =     0.0;     lVarBinHi[ DCANEGTOPV        ] =    10.0;
    lVarNbins[ DCAPOSTOPV        ] =   500;     lVarBinLo[ DCAPOSTOPV        ] =     0.0;     lVarBinHi[ DCAPOSTOPV        ] =    10.0;
    lVarNbins[ DCABACHTOPV       ] =   500;     lVarBinLo[ DCABACHTOPV       ] =     0.0;     lVarBinHi[ DCABACHTOPV       ] =    10.0;
    lVarNbins[ DCAV0TOPV         ] =   500;     lVarBinLo[ DCAV0TOPV         ] =     0.0;     lVarBinHi[ DCAV0TOPV         ] =    10.0;
    lVarNbins[ DCACASCTOPV       ] =   500;     lVarBinLo[ DCACASCTOPV       ] =     0.0;     lVarBinHi[ DCACASCTOPV       ] =    10.0;
    lVarNbins[ DCAXYCASCTOPV     ] =   500;     lVarBinLo[ DCAXYCASCTOPV     ] =     0.0;     lVarBinHi[ DCAXYCASCTOPV     ] =    10.0;
    lVarNbins[ DCAZCASCTOPV      ] =   500;     lVarBinLo[ DCAZCASCTOPV      ] =     0.0;     lVarBinHi[ DCAZCASCTOPV      ] =    10.0;
    lVarNbins[ DCAV0DAUGHTERS    ] =   400;     lVarBinLo[ DCAV0DAUGHTERS    ] =     0.0;     lVarBinHi[ DCAV0DAUGHTERS    ] =     2.0;
    lVarNbins[ DCACASCDAUGHTERS  ] =   400;     lVarBinLo[ DCACASCDAUGHTERS  ] =     0.0;     lVarBinHi[ DCACASCDAUGHTERS  ] =     2.0;
    lVarNbins[ DCAZNEGTOPV       ] =  1000;     lVarBinLo[ DCAZNEGTOPV       ] =  -100.0;     lVarBinHi[ DCAZNEGTOPV       ] =   100.0;
    lVarNbins[ DCAZPOSTOPV       ] =  1000;     lVarBinLo[ DCAZPOSTOPV       ] =  -100.0;     lVarBinHi[ DCAZPOSTOPV       ] =   100.0;
    lVarNbins[ DCAZBACHTOPV      ] =  1000;     lVarBinLo[ DCAZBACHTOPV      ] =  -100.0;     lVarBinHi[ DCAZBACHTOPV      ] =   100.0;
    lVarNbins[ DCABACHTOBARYON   ] =   500;     lVarBinLo[ DCABACHTOBARYON   ] =     0.0;     lVarBinHi[ DCABACHTOBARYON   ] =    10.0;
    lVarNbins[ V0PA              ] =  1000;     lVarBinLo[ V0PA              ] =     0.0;     lVarBinHi[ V0PA              ] =     1.0;
    lVarNbins[ CASCPA            ] =  1000;     lVarBinLo[ CASCPA            ] =     0.0;     lVarBinHi[ CASCPA            ] =     1.0;
    lVarNbins[ BBPA              ] =  1000;     lVarBinLo[ BBPA              ] =     0.0;     lVarBinHi[ BBPA              ] =     1.0;
    lVarNbins[ V0MASS            ] =   400;     lVarBinLo[ V0MASS            ] =     1.0;     lVarBinHi[ V0MASS            ] =     1.2;
    lVarNbins[ PROPERLIFETIME    ] =   500;     lVarBinLo[ PROPERLIFETIME    ] =     0.0;     lVarBinHi[ PROPERLIFETIME    ] =    50.0;
    lVarNbins[ COMPETINGSPECIES  ] =  1000;     lVarBinLo[ COMPETINGSPECIES  ] =     1.2;     lVarBinHi[ COMPETINGSPECIES  ] =     2.2;
    lVarNbins[ TPCNEGDEDX        ] =  1000;     lVarBinLo[ TPCNEGDEDX        ] =     0.0;     lVarBinHi[ TPCNEGDEDX        ] =   500.0;
    lVarNbins[ TPCPOSDEDX        ] =  1000;     lVarBinLo[ TPCPOSDEDX        ] =     0.0;     lVarBinHi[ TPCPOSDEDX        ] =   500.0;
    lVarNbins[ TPCBACHDEDX       ] =  1000;     lVarBinLo[ TPCBACHDEDX       ] =     0.0;     lVarBinHi[ TPCBACHDEDX       ] =   500.0;
    lVarNbins[ TPCNEGPIDNSIGMAS  ] =   200;     lVarBinLo[ TPCNEGPIDNSIGMAS  ] =   -10.0;     lVarBinHi[ TPCNEGPIDNSIGMAS  ] =    10.0;
    lVarNbins[ TPCPOSPIDNSIGMAS  ] =   200;     lVarBinLo[ TPCPOSPIDNSIGMAS  ] =   -10.0;     lVarBinHi[ TPCPOSPIDNSIGMAS  ] =    10.0;
    lVarNbins[ TPCBACHPIDNSIGMAS ] =   200;     lVarBinLo[ TPCBACHPIDNSIGMAS ] =   -10.0;     lVarBinHi[ TPCBACHPIDNSIGMAS ] =    10.0;
    lVarNbins[ TPCNCLUSTERS      ] =   160;     lVarBinLo[ TPCNCLUSTERS      ] =     0.0;     lVarBinHi[ TPCNCLUSTERS      ] =   160.0;
    lVarNbins[ MAXCHI2PERCLUSTER ] =  1000;     lVarBinLo[ MAXCHI2PERCLUSTER ] =     0.0;     lVarBinHi[ MAXCHI2PERCLUSTER ] =    10.0;
    lVarNbins[ MINTRACKLENGTH    ] =   200;     lVarBinLo[ MINTRACKLENGTH    ] =     0.0;     lVarBinHi[ MINTRACKLENGTH    ] =   200.0;
    lVarNbins[ NEGDAUGHTERETA    ] =   200;     lVarBinLo[ NEGDAUGHTERETA    ] =    -1.0;     lVarBinHi[ NEGDAUGHTERETA    ] =     1.0;
    lVarNbins[ POSDAUGHTERETA    ] =   200;     lVarBinLo[ POSDAUGHTERETA    ] =    -1.0;     lVarBinHi[ POSDAUGHTERETA    ] =     1.0;
    lVarNbins[ BACHDAUGHTERETA   ] =   200;     lVarBinLo[ BACHDAUGHTERETA   ] =    -1.0;     lVarBinHi[ BACHDAUGHTERETA   ] =     1.0;
    lVarNbins[ NEGINNERP         ] =  1000;     lVarBinLo[ NEGINNERP         ] =     0.0;     lVarBinHi[ NEGINNERP         ] =    10.0;
    lVarNbins[ POSINNERP         ] =  1000;     lVarBinLo[ POSINNERP         ] =     0.0;     lVarBinHi[ POSINNERP         ] =    10.0;
    lVarNbins[ BACHINNERP        ] =  1000;     lVarBinLo[ BACHINNERP        ] =     0.0;     lVarBinHi[ BACHINNERP        ] =    10.0;
    lVarNbins[ NEGTOFEXPTDIFF    ] =  1000;     lVarBinLo[ NEGTOFEXPTDIFF    ] =   -50.0;     lVarBinHi[ NEGTOFEXPTDIFF    ] =    50.0;
    lVarNbins[ POSTOFEXPTDIFF    ] =  1000;     lVarBinLo[ POSTOFEXPTDIFF    ] =   -50.0;     lVarBinHi[ POSTOFEXPTDIFF    ] =    50.0;
    lVarNbins[ BACHTOFEXPTDIFF   ] =  1000;     lVarBinLo[ BACHTOFEXPTDIFF   ] =   -50.0;     lVarBinHi[ BACHTOFEXPTDIFF   ] =    50.0;
    lVarNbins[ NEGTOFSIGNAL      ] =  2000;     lVarBinLo[ NEGTOFSIGNAL      ] =  -100.0;     lVarBinHi[ NEGTOFSIGNAL      ] =   100.0;
    lVarNbins[ POSTOFSIGNAL      ] =  2000;     lVarBinLo[ POSTOFSIGNAL      ] =  -100.0;     lVarBinHi[ POSTOFSIGNAL      ] =   100.0;
    lVarNbins[ BACHTOFSIGNAL     ] =  2000;     lVarBinLo[ BACHTOFSIGNAL     ] =  -100.0;     lVarBinHi[ BACHTOFSIGNAL     ] =   100.0;

    TList* lListOfVarHistosAll = new TList();
    TList* lListOfVarHistosSel = new TList();
    TList* lListOfVarHistosMCAll = new TList();
    TList* lListOfVarHistosMCSel = new TList();
    TList* lListOfVarHistosMCAssAll = new TList();
    TList* lListOfVarHistosMCAssSel = new TList();
    //
    TH2F** lHistSigPlusCenterBgAll = new TH2F*[NVARS];
    TH2F** lHistLeftPlusRightBgAll = new TH2F*[NVARS];
    TH2F** lHistSigPlusCenterBgSel = new TH2F*[NVARS];
    TH2F** lHistLeftPlusRightBgSel = new TH2F*[NVARS];
    //
    TH2F** lHistSigPlusCenterBgMCAll = new TH2F*[NVARS];
    TH2F** lHistLeftPlusRightBgMCAll = new TH2F*[NVARS];
    TH2F** lHistSigPlusCenterBgMCSel = new TH2F*[NVARS];
    TH2F** lHistLeftPlusRightBgMCSel = new TH2F*[NVARS];
    //
    TH2F** lHistSigPlusCenterBgMCAssAll = new TH2F*[NVARS];
    TH2F** lHistLeftPlusRightBgMCAssAll = new TH2F*[NVARS];
    TH2F** lHistSigPlusCenterBgMCAssSel = new TH2F*[NVARS];
    TH2F** lHistLeftPlusRightBgMCAssSel = new TH2F*[NVARS];

    for(Int_t ivar=0; ivar<NVARS; ++ivar) {
        if(fSaveVarHistosSwitch) {
            // data
            lHistSigPlusCenterBgAll[ivar] = new TH2F(Form("lHistSigPlusCenterBgDataAll%s", lVarName[ivar].Data()), Form("Signal+Center Bg (Data);#it{p}_{T} (GeV/c);%s", lVarName[ivar].Data()), fptbinnumb, fptbinlimits, lVarNbins[ivar], lVarBinLo[ivar], lVarBinHi[ivar]);
            lHistLeftPlusRightBgAll[ivar] = new TH2F(Form("lHistLeftPlusRightBgDataAll%s", lVarName[ivar].Data()), Form("Left+Right Bg (Data);#it{p}_{T} (GeV/c);%s", lVarName[ivar].Data()), fptbinnumb, fptbinlimits, lVarNbins[ivar], lVarBinLo[ivar], lVarBinHi[ivar]);
            lHistSigPlusCenterBgSel[ivar] = new TH2F(Form("lHistSigPlusCenterBgDataSel%s", lVarName[ivar].Data()), Form("Signal+Center Bg (Data);#it{p}_{T} (GeV/c);%s", lVarName[ivar].Data()), fptbinnumb, fptbinlimits, lVarNbins[ivar], lVarBinLo[ivar], lVarBinHi[ivar]);
            lHistLeftPlusRightBgSel[ivar] = new TH2F(Form("lHistLeftPlusRightBgDataSel%s", lVarName[ivar].Data()), Form("Left+Right Bg (Data);#it{p}_{T} (GeV/c);%s", lVarName[ivar].Data()), fptbinnumb, fptbinlimits, lVarNbins[ivar], lVarBinLo[ivar], lVarBinHi[ivar]);
            lListOfVarHistosAll->Add(lHistSigPlusCenterBgAll[ivar]);
            lListOfVarHistosAll->Add(lHistLeftPlusRightBgAll[ivar]);
            lListOfVarHistosSel->Add(lHistSigPlusCenterBgSel[ivar]);
            lListOfVarHistosSel->Add(lHistLeftPlusRightBgSel[ivar]);
            // mc
            lHistSigPlusCenterBgMCAll[ivar] = new TH2F(Form("lHistSigPlusCenterBgMCAll%s", lVarName[ivar].Data()), Form("Signal+Center Bg (MC);#it{p}_{T} (GeV/c);%s", lVarName[ivar].Data()), fptbinnumb, fptbinlimits, lVarNbins[ivar], lVarBinLo[ivar], lVarBinHi[ivar]);
            lHistLeftPlusRightBgMCAll[ivar] = new TH2F(Form("lHistLeftPlusRightBgMCAll%s", lVarName[ivar].Data()), Form("Left+Right Bg (MC);#it{p}_{T} (GeV/c);%s", lVarName[ivar].Data()), fptbinnumb, fptbinlimits, lVarNbins[ivar], lVarBinLo[ivar], lVarBinHi[ivar]);
            lHistSigPlusCenterBgMCSel[ivar] = new TH2F(Form("lHistSigPlusCenterBgMCSel%s", lVarName[ivar].Data()), Form("Signal+Center Bg (MC);#it{p}_{T} (GeV/c);%s", lVarName[ivar].Data()), fptbinnumb, fptbinlimits, lVarNbins[ivar], lVarBinLo[ivar], lVarBinHi[ivar]);
            lHistLeftPlusRightBgMCSel[ivar] = new TH2F(Form("lHistLeftPlusRightBgMCSel%s", lVarName[ivar].Data()), Form("Left+Right Bg (MC);#it{p}_{T} (GeV/c);%s", lVarName[ivar].Data()), fptbinnumb, fptbinlimits, lVarNbins[ivar], lVarBinLo[ivar], lVarBinHi[ivar]);
            lListOfVarHistosMCAll->Add(lHistSigPlusCenterBgMCAll[ivar]);
            lListOfVarHistosMCAll->Add(lHistLeftPlusRightBgMCAll[ivar]);
            lListOfVarHistosMCSel->Add(lHistSigPlusCenterBgMCSel[ivar]);
            lListOfVarHistosMCSel->Add(lHistLeftPlusRightBgMCSel[ivar]);
            // mc assoc
            lHistSigPlusCenterBgMCAssAll[ivar] = new TH2F(Form("lHistSigPlusCenterBgMCAssAll%s", lVarName[ivar].Data()), Form("Signal+Center Bg (MC);#it{p}_{T} (GeV/c);%s", lVarName[ivar].Data()), fptbinnumb, fptbinlimits, lVarNbins[ivar], lVarBinLo[ivar], lVarBinHi[ivar]);
            lHistLeftPlusRightBgMCAssAll[ivar] = new TH2F(Form("lHistLeftPlusRightBgMCAssAll%s", lVarName[ivar].Data()), Form("Left+Right Bg (MC);#it{p}_{T} (GeV/c);%s", lVarName[ivar].Data()), fptbinnumb, fptbinlimits, lVarNbins[ivar], lVarBinLo[ivar], lVarBinHi[ivar]);
            lHistSigPlusCenterBgMCAssSel[ivar] = new TH2F(Form("lHistSigPlusCenterBgMCAssSel%s", lVarName[ivar].Data()), Form("Signal+Center Bg (MC);#it{p}_{T} (GeV/c);%s", lVarName[ivar].Data()), fptbinnumb, fptbinlimits, lVarNbins[ivar], lVarBinLo[ivar], lVarBinHi[ivar]);
            lHistLeftPlusRightBgMCAssSel[ivar] = new TH2F(Form("lHistLeftPlusRightBgMCAssSel%s", lVarName[ivar].Data()), Form("Left+Right Bg (MC);#it{p}_{T} (GeV/c);%s", lVarName[ivar].Data()), fptbinnumb, fptbinlimits, lVarNbins[ivar], lVarBinLo[ivar], lVarBinHi[ivar]);
            lListOfVarHistosMCAssAll->Add(lHistSigPlusCenterBgMCAssAll[ivar]);
            lListOfVarHistosMCAssAll->Add(lHistLeftPlusRightBgMCAssAll[ivar]);
            lListOfVarHistosMCAssSel->Add(lHistSigPlusCenterBgMCAssSel[ivar]);
            lListOfVarHistosMCAssSel->Add(lHistLeftPlusRightBgMCAssSel[ivar]);
            //
            // rename x-axis of dE/dx histos
            if( ivar == TPCNEGDEDX || ivar == TPCPOSDEDX ) {
                lHistSigPlusCenterBgAll[ivar]->GetXaxis()->SetTitle("#it{p} (GeV/c)");
                lHistLeftPlusRightBgAll[ivar]->GetXaxis()->SetTitle("#it{p} (GeV/c)");
                lHistSigPlusCenterBgSel[ivar]->GetXaxis()->SetTitle("#it{p} (GeV/c)");
                lHistLeftPlusRightBgSel[ivar]->GetXaxis()->SetTitle("#it{p} (GeV/c)");
                lHistSigPlusCenterBgMCAll[ivar]->GetXaxis()->SetTitle("#it{p} (GeV/c)");
                lHistLeftPlusRightBgMCAll[ivar]->GetXaxis()->SetTitle("#it{p} (GeV/c)");
                lHistSigPlusCenterBgMCSel[ivar]->GetXaxis()->SetTitle("#it{p} (GeV/c)");
                lHistLeftPlusRightBgMCSel[ivar]->GetXaxis()->SetTitle("#it{p} (GeV/c)");
                lHistSigPlusCenterBgMCAssAll[ivar]->GetXaxis()->SetTitle("#it{p} (GeV/c)");
                lHistLeftPlusRightBgMCAssAll[ivar]->GetXaxis()->SetTitle("#it{p} (GeV/c)");
                lHistSigPlusCenterBgMCAssSel[ivar]->GetXaxis()->SetTitle("#it{p} (GeV/c)");
                lHistLeftPlusRightBgMCAssSel[ivar]->GetXaxis()->SetTitle("#it{p} (GeV/c)");
            }
        }
        else {
            lHistSigPlusCenterBgAll[ivar] = 0x0;
            lHistLeftPlusRightBgAll[ivar] = 0x0;
            lHistSigPlusCenterBgSel[ivar] = 0x0;
            lHistLeftPlusRightBgSel[ivar] = 0x0;
            lHistSigPlusCenterBgMCAll[ivar] = 0x0;
            lHistLeftPlusRightBgMCAll[ivar] = 0x0;
            lHistSigPlusCenterBgMCSel[ivar] = 0x0;
            lHistLeftPlusRightBgMCSel[ivar] = 0x0;
            lHistSigPlusCenterBgMCAssAll[ivar] = 0x0;
            lHistLeftPlusRightBgMCAssAll[ivar] = 0x0;
            lHistSigPlusCenterBgMCAssSel[ivar] = 0x0;
            lHistLeftPlusRightBgMCAssSel[ivar] = 0x0;
        }
    }

    //================================================================
    cout<<endl;





    cout<<"--------------- Open Real Data File --------------------"<<endl;
    TFile* file = TFile::Open(fRealDataFile, "READ");
    TList* clist      = (TList*)file->Get("PWGLF_StrVsMult/cList");
    TTree* lTreeEvent = (TTree*)file->Get("PWGLF_StrVsMult/fTreeEvent");
    TTree* lTree      = (TTree*)file->Get("PWGLF_StrVsMult/fTreeCascade");

    Double_t lLoMultBound =   0.;
    Double_t lHiMultBound = 100.;
    if(fPerformMultiplicityStudy) {
        lLoMultBound = fLoMultBound;
        lHiMultBound = fHiMultBound;
    }

    Float_t fCentrality = 0.;
    Bool_t  fMVPileupFlag = 0;
    Int_t   fClosestNonEmptyBC = 0;

    lTreeEvent->SetBranchAddress("fCentrality", &fCentrality);
    lTreeEvent->SetBranchAddress("fMVPileupFlag", &fMVPileupFlag);
    lTreeEvent->SetBranchAddress("fClosestNonEmptyBC", &fClosestNonEmptyBC);

    TH1F* fHistV0MultiplicityForTrigEvt;
    TH1F* fHistV0MultiplicityForSelEvtNoTPCOnly;
    TH1F* fHistV0MultiplicityForSelEvtNoTPCOnlyNoPileup;

    cout<<"--------------------------------------------------------"<<endl;
    cout<<" Will now loop over events, please wait..."<<endl;
    Long_t lNEvents = 0;
    for(Long_t iEv = 0; iEv<lTreeEvent->GetEntries(); iEv++) {
        lTreeEvent->GetEntry(iEv);
        if( iEv % ( lTreeEvent->GetEntries() / 10 ) == 0 ) cout<<" At Event "<<iEv<<" out of "<<lTreeEvent->GetEntries()<<endl;
        // check MV Pileup rejection
        if( fMVPileupSwitch && !fMVPileupFlag ) continue;
        // check distance to closest non empty BC
        if( TMath::Abs( fClosestNonEmptyBC ) < fMinDistToClosestNonEmptyBC ) continue; 
        //Multiplicity Switch
        if( fPerformMultiplicityStudy == kTRUE &&  //inside mult bin
            fCentrality>fLoMultBound &&
            fCentrality<fHiMultBound
          ) lNEvents++;
    }
    cout<<" Number of events, this multiplicity....: "<<lNEvents <<endl;
    cout<<"--------------------------------------------------------"<<endl;



    //Variable Definition=============================================
    //Kinematic
    Float_t lPt, lRap, lPtMC, lNegEta, lPosEta, lBachEta;
    Float_t lPosPx, lPosPy, lPosPz;
    Float_t lNegPx, lNegPy, lNegPz;
    Float_t lBachPx, lBachPy, lBachPz;
    //Invariant Masses
    Float_t lInvariantMass; //wildcard
    Float_t lCompetingParticleMass = -1; //Competing Species rejection
    //Topological variable definitions
    Float_t lDcaV0Daughters; //1
    Float_t lDcaPosToPrimVertex,  lDcaNegToPrimVertex, lDcaBachToPrimVertex; //2, 3, 4
    Float_t lDCAzPosToPrimVertex, lDCAzNegToPrimVertex, lDCAzBachToPrimVertex; // 12, 13, 14;
    Float_t lDCABachToBaryon;//15
    Float_t lDCAxyCascToPV, lDCAzCascToPV, lDCACascToPV;
    //Cosine of Pointing Angle variable
    Float_t lV0CosinePointingAngle, lCascCosinePointingAngle; //5, 6
    Float_t lV0CosinePointingAngleSpecial;
    Float_t lBBCosPA;
    //Decay Radius and distance over total momentum
    Float_t lV0Radius, lCascRadius, lDistOverTotMom; //7, 8
    Float_t lV0Mass, lDcaCascDaughters; //9, 10
    Float_t lDcaV0ToPV; //11
    //Least Number of TPC Clusters
    Int_t lLeastNbrClusters;
    Float_t lMinTrackLength;
    Float_t lMaxChi2PerCluster;
    //TPC dE/dx acquired with AliPIDResponse class
    Float_t lNSigmasPosProton,lNSigmasNegProton,lNSigmasPosPion,lNSigmasNegPion,
            lNSigmasBachPion, lNSigmasBachKaon;
    Float_t lPosdEdx, lNegdEdx, lBachdEdx;
    Float_t lPosInnerP, lNegInnerP, lBachInnerP;
    //ITS
    ULong64_t lNegTrackStatus, lPosTrackStatus, lBachTrackStatus;
    //TOF
    Float_t lNegTOFExpTDiff, lPosTOFExpTDiff, lBachTOFExpTDiff;
    Float_t lNegTOFSignal, lPosTOFSignal, lBachTOFSignal;
    //Charge
    Int_t lCharge = 0;
    //Multiplicity Variable
    Float_t lMultiplicity = -1.;
    //================================================================

    //Linking to Tree=================================================
    //--- Base Variables ----------------------------------------------
    lTree->SetBranchAddress("fTreeCascVarCharge"  ,&lCharge );
    lTree->SetBranchAddress("fTreeCascVarPosEta"  ,&lPosEta );
    lTree->SetBranchAddress("fTreeCascVarNegEta"  ,&lNegEta );
    lTree->SetBranchAddress("fTreeCascVarBachEta" ,&lBachEta);
    lTree->SetBranchAddress("fTreeCascVarPt",&lPt);
    if ( fWhichParticle == "XiMinus"      )  lTree->SetBranchAddress("fTreeCascVarMassAsXi"   ,&lInvariantMass);
    if ( fWhichParticle == "XiPlus"       )  lTree->SetBranchAddress("fTreeCascVarMassAsXi"   ,&lInvariantMass);
    if ( fWhichParticle == "OmegaMinus" || fWhichParticle == "OmegaPlus" ){
        lTree->SetBranchAddress("fTreeCascVarMassAsOmega",&lInvariantMass);
        lTree->SetBranchAddress("fTreeCascVarMassAsXi"   ,&lCompetingParticleMass); //For Competing Rejection
    }
    if ( fWhichParticle == "XiMinus"    || fWhichParticle == "XiPlus"    )
        lTree->SetBranchAddress("fTreeCascVarRapXi"   ,&lRap);
    if ( fWhichParticle == "OmegaMinus" || fWhichParticle == "OmegaPlus" )
        lTree->SetBranchAddress("fTreeCascVarRapOmega",&lRap);
    lTree->SetBranchAddress("fTreeCascVarLeastNbrClusters",&lLeastNbrClusters);
    lTree->SetBranchAddress("fTreeCascVarPosPx"     , &lPosPx);
    lTree->SetBranchAddress("fTreeCascVarPosPy"     , &lPosPy);
    lTree->SetBranchAddress("fTreeCascVarPosPz"     , &lPosPz);
    lTree->SetBranchAddress("fTreeCascVarNegPx"     , &lNegPx);
    lTree->SetBranchAddress("fTreeCascVarNegPy"     , &lNegPy);
    lTree->SetBranchAddress("fTreeCascVarNegPz"     , &lNegPz);
    lTree->SetBranchAddress("fTreeCascVarBachPx"    , &lBachPx);
    lTree->SetBranchAddress("fTreeCascVarBachPy"    , &lBachPy);
    lTree->SetBranchAddress("fTreeCascVarBachPz"    , &lBachPz);
    lTree->SetBranchAddress("fTreeCascVarPosInnerP",  &lPosInnerP);
    lTree->SetBranchAddress("fTreeCascVarNegInnerP",  &lNegInnerP);
    lTree->SetBranchAddress("fTreeCascVarBachInnerP", &lBachInnerP);
    //--- TPC Variables -----------------------------------------------
    lTree->SetBranchAddress("fTreeCascVarPosNSigmaProton",&lNSigmasPosProton);
    lTree->SetBranchAddress("fTreeCascVarNegNSigmaProton",&lNSigmasNegProton);
    lTree->SetBranchAddress("fTreeCascVarPosNSigmaPion",&lNSigmasPosPion);
    lTree->SetBranchAddress("fTreeCascVarNegNSigmaPion",&lNSigmasNegPion);
    lTree->SetBranchAddress("fTreeCascVarBachNSigmaPion",&lNSigmasBachPion);
    lTree->SetBranchAddress("fTreeCascVarBachNSigmaKaon",&lNSigmasBachKaon);
    lTree->SetBranchAddress("fTreeCascVarPosdEdx",  &lPosdEdx);
    lTree->SetBranchAddress("fTreeCascVarNegdEdx",  &lNegdEdx);
    lTree->SetBranchAddress("fTreeCascVarBachdEdx", &lBachdEdx);
    lTree->SetBranchAddress("fTreeCascVarMinTrackLength", &lMinTrackLength);
    lTree->SetBranchAddress("fTreeCascVarMaxChi2PerCluster", &lMaxChi2PerCluster);
    //--- Topological selection variables -----------------------------
    lTree->SetBranchAddress("fTreeCascVarV0Radius",&lV0Radius); //1
    lTree->SetBranchAddress("fTreeCascVarCascRadius",&lCascRadius); //2
    lTree->SetBranchAddress("fTreeCascVarV0Mass",&lV0Mass); //3
    lTree->SetBranchAddress("fTreeCascVarV0CosPointingAngle",&lV0CosinePointingAngle); //4
    lTree->SetBranchAddress("fTreeCascVarV0CosPointingAngleSpecial",&lV0CosinePointingAngleSpecial);
    lTree->SetBranchAddress("fTreeCascVarCascCosPointingAngle",&lCascCosinePointingAngle); //5
    lTree->SetBranchAddress("fTreeCascVarDCANegToPrimVtx",&lDcaNegToPrimVertex); //6
    lTree->SetBranchAddress("fTreeCascVarDCAPosToPrimVtx",&lDcaPosToPrimVertex); //7
    lTree->SetBranchAddress("fTreeCascVarDCABachToPrimVtx",&lDcaBachToPrimVertex); //8
    lTree->SetBranchAddress("fTreeCascVarDCAV0ToPrimVtx",&lDcaV0ToPV); //9
    lTree->SetBranchAddress("fTreeCascVarCascDCAtoPVxy", &lDCAxyCascToPV);
    lTree->SetBranchAddress("fTreeCascVarCascDCAtoPVz", &lDCAzCascToPV);
    lTree->SetBranchAddress("fTreeCascVarDCAV0Daughters",&lDcaV0Daughters); //10
    lTree->SetBranchAddress("fTreeCascVarDCACascDaughters",&lDcaCascDaughters); //11
    lTree->SetBranchAddress("fTreeCascVarDistOverTotMom",&lDistOverTotMom); //11
    lTree->SetBranchAddress("fTreeCascVarNegDCAz",  &lDCAzNegToPrimVertex);
    lTree->SetBranchAddress("fTreeCascVarPosDCAz",  &lDCAzPosToPrimVertex);
    lTree->SetBranchAddress("fTreeCascVarBachDCAz", &lDCAzBachToPrimVertex);
    lTree->SetBranchAddress("fTreeCascVarDCABachToBaryon", &lDCABachToBaryon);
    lTree->SetBranchAddress("fTreeCascVarWrongCosPA", &lBBCosPA);
    //--- ITS flag -----------------------------------------------------
    lTree->SetBranchAddress("fTreeCascVarPosTrackStatus", &lPosTrackStatus);
    lTree->SetBranchAddress("fTreeCascVarNegTrackStatus", &lNegTrackStatus);
    lTree->SetBranchAddress("fTreeCascVarBachTrackStatus", &lBachTrackStatus);
    //--- TOF info -----------------------------------------------------
    lTree->SetBranchAddress("fTreeCascVarNegTOFExpTDiff",  &lNegTOFExpTDiff);
    lTree->SetBranchAddress("fTreeCascVarPosTOFExpTDiff",  &lPosTOFExpTDiff);
    lTree->SetBranchAddress("fTreeCascVarBachTOFExpTDiff", &lBachTOFExpTDiff);
    lTree->SetBranchAddress("fTreeCascVarNegTOFSignal",  &lNegTOFSignal);
    lTree->SetBranchAddress("fTreeCascVarPosTOFSignal",  &lPosTOFSignal);
    lTree->SetBranchAddress("fTreeCascVarBachTOFSignal", &lBachTOFSignal);
    //--- Multiplicity Variable ----------------------------------------
    lTree->SetBranchAddress("fTreeCascVarCentrality",&fCentrality);
    //--- MV pileup flag -----------------------------------------------
    lTree->SetBranchAddress("fTreeCascVarMVPileupFlag", &fMVPileupFlag);
    lTree->SetBranchAddress("fTreeCascVarClosestNonEmptyBC", &fClosestNonEmptyBC);
    //================================================================

    //== Variables for holding peak position, width ==============
    Double_t lPeakPosition[100];
    Double_t lPeakWidth[100];
    Double_t lLeftBgLeftLimit[100];
    Double_t lLeftBgRightLimit[100];
    Double_t lPeakLeftLimit[100];
    Double_t lPeakRightLimit[100];
    Double_t lRightBgLeftLimit[100];
    Double_t lRightBgRightLimit[100];
    //May be needed if bg areas are different in the future
    //Double_t lScaleFactor[100];

    TLine *lLineLeftMost[100];
    TLine *lLineLeft[100];
    TLine *lLineRight[100];
    TLine *lLineRightMost[100];

    TLine *lLineLeftMostMC[100];
    TLine *lLineLeftMC[100];
    TLine *lLineRightMC[100];
    TLine *lLineRightMostMC[100];

    char fgausname[100];
    TF1 *fgausPt[100];

    Long_t lNCandidates = lTree->GetEntries();
    Long_t lOneTenthOfNCandidates = ((double)(lNCandidates) / 10. );

    if(!fUsePeakPositionAndWidthFromFit) {
        cout<<"--------------- Real Data File Loop 1 ------------------"<<endl;
        for(Long_t icand = 0;icand<lNCandidates;icand++){
            lTree->GetEntry(icand);

            // check MV Pileup rejection
            if( fMVPileupSwitch && !fMVPileupFlag ) continue;

            // check distance to closest non empty BC
            if( TMath::Abs( fClosestNonEmptyBC ) < fMinDistToClosestNonEmptyBC ) continue; 

            //Multiplicity Switch -- use integrated sample for peak finding
            lMultiplicity = (Double_t)fCentrality;
            if( fPerformMultiplicityStudy && (lMultiplicity<0. || lMultiplicity>100.) ) continue;

            if( icand % lOneTenthOfNCandidates == 0 )
                cout<<" Currently at candidate........: "<<icand<<" / "<<lNCandidates<<" ( "<<(long)(((double)(icand)/(double)(lNCandidates))*(100.+1e-3))<<"% )"<<endl;

            //CMS Shift: apply before!
            if(fRapidityType == "CMS") lRap = lRap + fRapidityShift;

            lMultiplicity = (Double_t)fCentrality;

            //Compute 3D DCA Cascade to PV
            lDCACascToPV = TMath::Sqrt( lDCAxyCascToPV*lDCAxyCascToPV + lDCAzCascToPV*lDCAzCascToPV );

            //Now check validity
            if( lRap<fRapidityBoundaryUpper && lRap>fRapidityBoundaryLower &&
                (//charge condition (x-check)
                  (fWhichParticle == "XiMinus"    && lCharge == -1) ||
                  (fWhichParticle == "XiPlus"     && lCharge ==  1) ||
                  (fWhichParticle == "OmegaMinus" && lCharge == -1) ||
                  (fWhichParticle == "OmegaPlus"  && lCharge ==  1)
                ) &&
                TMath::Abs(lNegEta)       <= GetCutDaughterEta      ( lPt ) &&
                TMath::Abs(lPosEta)       <= GetCutDaughterEta      ( lPt ) &&
                TMath::Abs(lBachEta)      <= GetCutDaughterEta      ( lPt ) &&
                //Topological Selections
                lV0Radius                 >= GetCutV0Radius         ( lPt ) &&
                lCascRadius               >= GetCutCascRadius       ( lPt ) &&
                TMath::Abs(lV0Mass-1.116) <= GetCutV0Mass           ( lPt ) &&
                lV0CosinePointingAngle    >= GetCutV0CosPA          ( lPt ) &&
                lCascCosinePointingAngle  >= GetCutCascCosPA        ( lPt ) &&
                lDcaNegToPrimVertex       >= GetCutDCANegToPV       ( lPt ) &&
                lDcaPosToPrimVertex       >= GetCutDCAPosToPV       ( lPt ) &&
                lDcaBachToPrimVertex      >= GetCutDCABachToPV      ( lPt ) &&
                lDcaV0Daughters           <= GetCutDCAV0Daughters   ( lPt ) &&
                lDcaCascDaughters         <= GetCutDCACascDaughters ( lPt ) &&
                lDcaV0ToPV                >= GetCutDCAV0ToPV        ( lPt ) &&
                lDCACascToPV              <= GetCutDCACascToPV      ( lPt ) &&
                lDCAxyCascToPV            <= GetCutDCAxyCascToPV    ( lPt ) &&
                lDCAzCascToPV             <= GetCutDCAzCascToPV     ( lPt ) &&
                lParticleMass*lDistOverTotMom     <= GetCutProperLifetime   ( lPt ) &&
                lDCAzNegToPrimVertex      <= GetCutDCAzNegToPV      ( lPt ) &&
                lDCAzPosToPrimVertex      <= GetCutDCAzPosToPV      ( lPt ) &&
                lDCAzBachToPrimVertex     <= GetCutDCAzBachToPV     ( lPt ) &&
                lDCABachToBaryon          >= GetCutDCABachToBaryon  ( lPt ) &&
                lBBCosPA                  <= GetCutBBCosPA          ( lPt ) &&
                lMinTrackLength           >= GetCutMinTrackLength   ( lPt ) &&

                //Competing Species Rejection (only for Omegas)
                TMath::Abs( lCompetingParticleMass - 1.32171 ) >= fCutCompetingSpecies &&

                //Causality cut
                TMath::Abs( lV0Radius - lCascRadius ) >= GetCutCausality( lPt ) &&

                //Nclusters cut
                lLeastNbrClusters >= GetCutLeastNumberOfClusters( lPt ) &&

                ( //official response code
                  ( fWhichParticle == "XiMinus"
                    && TMath::Abs(lNSigmasNegPion)   <= GetCutTPCPIDNSigmas( lPt )
                    && TMath::Abs(lNSigmasPosProton) <= GetCutTPCPIDNSigmas( lPt )
                    && TMath::Abs(lNSigmasBachPion)  <= GetCutTPCPIDNSigmas( lPt )) ||
                  ( fWhichParticle == "XiPlus"
                    && TMath::Abs(lNSigmasPosPion)   <= GetCutTPCPIDNSigmas( lPt )
                    && TMath::Abs(lNSigmasNegProton) <= GetCutTPCPIDNSigmas( lPt )
                    && TMath::Abs(lNSigmasBachPion)  <= GetCutTPCPIDNSigmas( lPt )) ||
                  ( fWhichParticle == "OmegaMinus"
                    && TMath::Abs(lNSigmasNegPion)   <= GetCutTPCPIDNSigmas( lPt )
                    && TMath::Abs(lNSigmasPosProton) <= GetCutTPCPIDNSigmas( lPt )
                    && TMath::Abs(lNSigmasBachKaon)  <= GetCutTPCPIDNSigmas( lPt )) ||
                  ( fWhichParticle == "OmegaPlus"
                    && TMath::Abs(lNSigmasPosPion)   <= GetCutTPCPIDNSigmas( lPt )
                    && TMath::Abs(lNSigmasNegProton) <= GetCutTPCPIDNSigmas( lPt )
                    && TMath::Abs(lNSigmasBachKaon)  <= GetCutTPCPIDNSigmas( lPt ))
                ) &&

                ( //TOF match
                  CheckTOFmatchBach( lPt, lBachTOFExpTDiff ) &&
                  CheckTOFmatchNeg(  lPt, lNegTOFExpTDiff  ) &&
                  CheckTOFmatchPos(  lPt, lPosTOFExpTDiff  ) &&
                  CheckTOFmatchOne(  lPt, lBachTOFExpTDiff, lNegTOFExpTDiff, lPosTOFExpTDiff )
                ) &&
                 
                ( //ITS refit 
                  CheckITSrefitBach( lPt, lBachTrackStatus ) &&
                  CheckITSrefitNeg(  lPt, lNegTrackStatus  ) &&
                  CheckITSrefitPos(  lPt, lPosTrackStatus  ) &&
                  CheckITSrefitOne(  lPt, lBachTrackStatus, lNegTrackStatus, lPosTrackStatus )
                ) &&
                
                ( //ITSTOF 
                  CheckITSTOFBach( lPt, lBachTrackStatus, lBachTOFExpTDiff ) &&
                  CheckITSTOFNeg(  lPt, lNegTrackStatus,  lNegTOFExpTDiff  ) &&
                  CheckITSTOFPos(  lPt, lPosTrackStatus,  lPosTOFExpTDiff  ) &&
                  CheckITSTOFOne(  lPt, lBachTrackStatus, lNegTrackStatus, lPosTrackStatus, lBachTOFExpTDiff, lNegTOFExpTDiff, lPosTOFExpTDiff )
                )

              ) { // Start Entry Loop
                  lWeAreAtBin = fHistPt->FindBin(lPt)-1;
                  if(lWeAreAtBin == -1) lWeAreAtBin = 99; //UnderFlow, special treatment
                  lHistoMBCasc[lWeAreAtBin]->Fill(lInvariantMass); //fill with specific inv mass
            }
        }
        cout<<"--------------- Loop Completed -------------------------"<<endl;
        cout<<endl;

        cout<<"--------------- Peak Finding (gauss+linear) ------------"<<endl;
        // Define temp histos to store peak position and width
        TH1F *lHistPeakPosition = new TH1F(Form("PeakPosition_%s", fWhichParticle.Data()), "Peak Position", fptbinnumb, fptbinlimits);
        TH1F *lHistPeakWidth    = new TH1F(Form("PeakWidth_%s", fWhichParticle.Data()),    "Peak Width",    fptbinnumb, fptbinlimits);

        for(Int_t ibin = 0; ibin<fptbinnumb; ibin++){
            cout<<"---> Peak Finding, bin #"<<ibin<<"..."<<endl;
            sprintf(fgausname,"fGausPt%i",ibin);
            if ( fWhichParticle == "XiMinus" || fWhichParticle == "XiPlus" ){
                fgausPt[ibin]= new TF1(fgausname,"[0]*TMath::Gaus(x,[1],[2])+[3]*x+[4]", 1.322-0.03, 1.322+0.03 );
                fgausPt[ibin]->SetParameter(1,1.322);
                fgausPt[ibin]->SetParameter(2,0.0025);
                fgausPt[ibin]->SetParLimits(2,0.001,0.01);
            }
            if ( fWhichParticle == "OmegaMinus" || fWhichParticle == "OmegaPlus"){
                fgausPt[ibin]= new TF1(fgausname,"[0]*TMath::Gaus(x,[1],[2])+[3]*x+[4]", 1.672-0.03, 1.672+0.03 );
                fgausPt[ibin]->SetParameter(1,1.672);
                fgausPt[ibin]->SetParameter(2,0.0025);
                fgausPt[ibin]->SetParLimits(2,0.001,0.01);
            }
            fgausPt[ibin]->SetParameter(0,lHistoMBCasc[ibin]->GetMaximum() * 0.9);
            fgausPt[ibin]->SetParameter(3,0);
            fgausPt[ibin]->SetParameter(4,lHistoMBCasc[ibin]->GetMaximum() * 0.1);
            lHistoMBCasc[ibin]->Fit(fgausname,"QREM0");
            lPeakPosition[ibin] = fgausPt[ibin]->GetParameter(1);
            lPeakWidth[ibin] = TMath::Abs( fgausPt[ibin]->GetParameter(2) );
            cout<<"---> ["<<fptbinlimits[ibin]<<" - "<<fptbinlimits[ibin+1]<<" GeV/c]\tPeak at: "<<lPeakPosition[ibin]<<", sigma = "<<lPeakWidth[ibin]<<endl;
            //Find Corresponding Limits In this bin
            lLeftBgLeftLimit[ibin]  = lPeakPosition[ibin] - 2.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin];
            lLeftBgRightLimit[ibin] = lPeakPosition[ibin] - 1.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin];
            lPeakLeftLimit[ibin]    = lPeakPosition[ibin] - 1.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin];
            lPeakRightLimit[ibin]   = lPeakPosition[ibin] + 1.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin];
            lRightBgLeftLimit[ibin] = lPeakPosition[ibin] + 1.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin];
            lRightBgRightLimit[ibin]= lPeakPosition[ibin] + 2.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin];
            //lScaleFactor[ibin] = (lPeakRightLimit[ibin] - lPeakLeftLimit[ibin]) / ( (lLeftBgRightLimit[ibin]-lLeftBgLeftLimit[ibin]) + (lRightBgRightLimit[ibin] - lRightBgLeftLimit[ibin]) );
            fHistPeakPosition->SetBinContent(ibin+1, lPeakPosition[ibin]);
            fHistPeakPosition->SetBinError(ibin+1, lPeakWidth[ibin]);
            //
            lHistPeakPosition->SetBinContent(ibin+1, lPeakPosition[ibin]);
            lHistPeakPosition->SetBinError(ibin+1, fgausPt[ibin]->GetParError(1));
            lHistPeakWidth->SetBinContent(ibin+1, lPeakWidth[ibin]);
            lHistPeakWidth->SetBinError(ibin+1, fgausPt[ibin]->GetParError(2));
            //
            //Create Signal Extraction Range Histogram
            fHistSignalExtractionRange->SetBinContent(ibin+1, lPeakPosition[ibin]);
            fHistSignalExtractionRange->SetBinError(ibin+1, 2.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin] );
            //Create appropriate TLine Objects for Canvases
            lLineLeftMost[ibin]  = new TLine( lLeftBgLeftLimit[ibin],   0, lLeftBgLeftLimit[ibin],   lHistoMBCasc[ibin]->GetMaximum() * 0.95 );
            lLineLeft[ibin]      = new TLine( lLeftBgRightLimit[ibin],  0, lLeftBgRightLimit[ibin],  lHistoMBCasc[ibin]->GetMaximum() * 0.95 );
            lLineRight[ibin]     = new TLine( lRightBgLeftLimit[ibin],  0, lRightBgLeftLimit[ibin],  lHistoMBCasc[ibin]->GetMaximum() * 0.95 );
            lLineRightMost[ibin] = new TLine( lRightBgRightLimit[ibin], 0, lRightBgRightLimit[ibin], lHistoMBCasc[ibin]->GetMaximum() * 0.95 );

            //Preparing Canvas for storing...
            lCanvasHistoMBCasc[ibin]->cd();
            lHistoMBCasc[ibin]->Draw();
            lLineLeftMost[ibin]->Draw();
            lLineLeft[ibin]->Draw();
            lLineRight[ibin]->Draw();
            lLineRightMost[ibin]->Draw();
            fgausPt[ibin]->Draw("same");

            //Preparing Canvas for storing, selection level
            lCanvasHistoSelectedCasc[ibin]->cd();
            lHistoSelectedCasc[ibin]->Draw();
            lLineLeftMost[ibin]->Draw();
            lLineLeft[ibin]->Draw();
            lLineRight[ibin]->Draw();
            lLineRightMost[ibin]->Draw();

        }
        // Save peak position and signal extraction range
        TString fOutputSigExtParamsName = fOutputDataFile;
        fOutputSigExtParamsName.ReplaceAll("Results", "SigExtParams");
        TFile *fOutputSigExtParams = new TFile(fOutputSigExtParamsName.Data(), "RECREATE");
        lHistPeakPosition->Write();
        lHistPeakWidth->Write();
        fOutputSigExtParams->Close();
        delete lHistPeakPosition; lHistPeakPosition = 0;
        delete lHistPeakWidth;    lHistPeakWidth    = 0;
        cout<<"--------------- Peak Finding Finished ------------------"<<endl;
    }
    else {
        // Use peak position and width provided by the user
        cout<<"---> Peak Position and Width (provided by the user):"<<endl;
        if( !(fPeakPositionFit && fPeakWidthFit) ) {
            cout<<"[AliCascadeModule] Error: functions for peak position and width not found!"<<endl;
            exit(-1);
        }
        cout<<"     Peak Position Fit: "<<fPeakPositionFit->GetName()<<endl;
        cout<<"     Peak Width Fit...: "<<fPeakWidthFit->GetName()<<endl;
        fHistPeakPosition->SetTitle(Form("%s (provided by the user)", fHistPeakPosition->GetTitle()));
        fHistSignalExtractionRange->SetTitle(Form("%s (provided by the user)", fHistSignalExtractionRange->GetTitle()));
        for(Int_t ibin = 0; ibin<fptbinnumb; ibin++){
            //cout<<"---> pT bin #"<<ibin<<"..."<<endl;
            Double_t lThisPtBinCenter = (fptbinlimits[ibin]+fptbinlimits[ibin+1])/2.;
            lPeakPosition[ibin] = fPeakPositionFit->Eval( lThisPtBinCenter );
            lPeakWidth[ibin]    = fPeakWidthFit->Eval( lThisPtBinCenter );
            cout<<"---> ["<<fptbinlimits[ibin]<<" - "<<fptbinlimits[ibin+1]<<" GeV/c]\tPeak at: "<<lPeakPosition[ibin]<<", sigma = "<<lPeakWidth[ibin]<<endl;
            //Find Corresponding Limits In this bin
            lLeftBgLeftLimit[ibin]  = lPeakPosition[ibin] - 2.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin];
            lLeftBgRightLimit[ibin] = lPeakPosition[ibin] - 1.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin];
            lPeakLeftLimit[ibin]    = lPeakPosition[ibin] - 1.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin];
            lPeakRightLimit[ibin]   = lPeakPosition[ibin] + 1.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin];
            lRightBgLeftLimit[ibin] = lPeakPosition[ibin] + 1.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin];
            lRightBgRightLimit[ibin]= lPeakPosition[ibin] + 2.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin];
            //Create appropriate TLine Objects for Canvases
            lLineLeftMost[ibin]  = new TLine( lLeftBgLeftLimit[ibin],   0, lLeftBgLeftLimit[ibin],   lHistoMBCasc[ibin]->GetMaximum() * 0.95 );
            lLineLeft[ibin]      = new TLine( lLeftBgRightLimit[ibin],  0, lLeftBgRightLimit[ibin],  lHistoMBCasc[ibin]->GetMaximum() * 0.95 );
            lLineRight[ibin]     = new TLine( lRightBgLeftLimit[ibin],  0, lRightBgLeftLimit[ibin],  lHistoMBCasc[ibin]->GetMaximum() * 0.95 );
            lLineRightMost[ibin] = new TLine( lRightBgRightLimit[ibin], 0, lRightBgRightLimit[ibin], lHistoMBCasc[ibin]->GetMaximum() * 0.95 );
            //Copy to the default histos to be saved
            fHistPeakPosition->SetBinContent(ibin+1, lPeakPosition[ibin]);
            fHistPeakPosition->SetBinError(ibin+1, lPeakWidth[ibin]);
            //Create Signal Extraction Range Histogram
            fHistSignalExtractionRange->SetBinContent(ibin+1, lPeakPosition[ibin]);
            fHistSignalExtractionRange->SetBinError(ibin+1, 2.*fCutNSigmasForSignalExtraction*lPeakWidth[ibin] );
            //Preparing Canvas for storing...
            lCanvasHistoMBCasc[ibin]->cd();
            lHistoMBCasc[ibin]->Draw();
            lLineLeftMost[ibin]->Draw();
            lLineLeft[ibin]->Draw();
            lLineRight[ibin]->Draw();
            lLineRightMost[ibin]->Draw();
            //fgausPt[ibin]->Draw("same");
            //Preparing Canvas for storing, selection level
            lCanvasHistoSelectedCasc[ibin]->cd();
            lHistoSelectedCasc[ibin]->Draw();
            lLineLeftMost[ibin]->Draw();
            lLineLeft[ibin]->Draw();
            lLineRight[ibin]->Draw();
            lLineRightMost[ibin]->Draw();
        }
    }

    //Defining Signal holding variables===============================
    Double_t lSigRealCasc[100];
    Double_t lSigErrRealCasc[100];
    Double_t lSigMCCasc[100];
    Double_t lSigErrMCCasc[100];
    //================================================================
    Double_t lLeftPlusRightBgCasc[100];
    Double_t lSigPlusCenterBgCasc[100];
    Double_t lLeftPlusRightBgCascMC[100];
    Double_t lSigPlusCenterBgCascMC[100];
    Double_t lLeftPlusRightBgTemplateCascMC[100];
    Double_t lCenterBgTemplateCascMC[100];
    for(Long_t n1 = 0; n1<100; n1++){
        lLeftPlusRightBgCasc[n1]=0;
        lSigPlusCenterBgCasc[n1]=0;
        lLeftPlusRightBgCascMC[n1]=0;
        lSigPlusCenterBgCascMC[n1]=0;
        lLeftPlusRightBgTemplateCascMC[n1]=0;
        lCenterBgTemplateCascMC[n1]=0;
    }
    //================================================================
    cout<<endl;
    cout<<"--------------- Real Data File Loop 2 ------------------"<<endl;
    for(Long_t icand = 0;icand<lNCandidates;icand++){
        lTree->GetEntry(icand);

        // check MV Pileup rejection
        if( fMVPileupSwitch && !fMVPileupFlag ) continue;

        // check distance to closest non empty BC
        if( TMath::Abs( fClosestNonEmptyBC ) < fMinDistToClosestNonEmptyBC ) continue; 

        //Multiplicity Switch
        lMultiplicity = (Double_t)fCentrality;
        if( fPerformMultiplicityStudy && (lMultiplicity<lLoMultBound || lMultiplicity>lHiMultBound) ) continue;

        if( icand % lOneTenthOfNCandidates == 0 )
            cout<<" Currently at candidate........: "<<icand<<" / "<<lNCandidates<<" ( "<<(long)(((double)(icand)/(double)(lNCandidates))*(100.+1e-3))<<"% )"<<endl;

        //CMS Shift: apply before!
        if(fRapidityType == "CMS") lRap = lRap + fRapidityShift; //DON'T CARE, this is pp

        //Compute 3D DCA Cascade to PV
        lDCACascToPV = TMath::Sqrt( lDCAxyCascToPV*lDCAxyCascToPV + lDCAzCascToPV*lDCAzCascToPV );

        //--- TPC dE/dx QA
        Float_t lPosP  = TMath::Sqrt( lPosPx*lPosPx + lPosPy*lPosPy + lPosPz*lPosPz );
        Float_t lNegP  = TMath::Sqrt( lNegPx*lNegPx + lNegPy*lNegPy + lNegPz*lNegPz );
        Float_t lBachP = TMath::Sqrt( lBachPx*lBachPx + lBachPy*lBachPy + lBachPz*lBachPz );
        Float_t lPosNSigmas  = -999;
        Float_t lNegNSigmas  = -999;
        Float_t lBachNSigmas = -999;
        if(fWhichParticle=="XiMinus")    { lPosNSigmas = lNSigmasPosProton; lNegNSigmas = lNSigmasNegPion;   lBachNSigmas = lNSigmasBachPion;  }
        if(fWhichParticle=="XiPlus")     { lPosNSigmas = lNSigmasPosPion;   lNegNSigmas = lNSigmasNegProton; lBachNSigmas = lNSigmasBachPion;  }
        if(fWhichParticle=="OmegaMinus") { lPosNSigmas = lNSigmasPosProton; lNegNSigmas = lNSigmasNegPion;   lBachNSigmas = lNSigmasBachKaon;  }
        if(fWhichParticle=="OmegaPlus")  { lPosNSigmas = lNSigmasPosPion;   lNegNSigmas = lNSigmasNegProton; lBachNSigmas = lNSigmasBachKaon;  }

        lWeAreAtBin = fHistPt->FindBin(lPt)-1;
        if(lWeAreAtBin == -1) lWeAreAtBin = 99; //UnderFlow, special treatment

        if(fSaveVarHistosSwitch) {
            //--- Fill histograms for signal+center background candidates
            if(lInvariantMass>lPeakLeftLimit[lWeAreAtBin] && lInvariantMass<lPeakRightLimit[lWeAreAtBin]) {
                lHistSigPlusCenterBgAll[ V0RADIUS          ] -> Fill( lPt, lV0Radius );
                lHistSigPlusCenterBgAll[ CASCRADIUS        ] -> Fill( lPt, lCascRadius );
                lHistSigPlusCenterBgAll[ DCANEGTOPV        ] -> Fill( lPt, lDcaNegToPrimVertex );
                lHistSigPlusCenterBgAll[ DCAPOSTOPV        ] -> Fill( lPt, lDcaPosToPrimVertex );
                lHistSigPlusCenterBgAll[ DCABACHTOPV       ] -> Fill( lPt, lDcaBachToPrimVertex );
                lHistSigPlusCenterBgAll[ DCAV0TOPV         ] -> Fill( lPt, lDcaV0ToPV );
                lHistSigPlusCenterBgAll[ DCACASCTOPV       ] -> Fill( lPt, lDCACascToPV );
                lHistSigPlusCenterBgAll[ DCAXYCASCTOPV     ] -> Fill( lPt, lDCAxyCascToPV );
                lHistSigPlusCenterBgAll[ DCAZCASCTOPV      ] -> Fill( lPt, lDCAzCascToPV );
                lHistSigPlusCenterBgAll[ DCAV0DAUGHTERS    ] -> Fill( lPt, lDcaV0Daughters );
                lHistSigPlusCenterBgAll[ DCACASCDAUGHTERS  ] -> Fill( lPt, lDcaCascDaughters );
                lHistSigPlusCenterBgAll[ DCAZNEGTOPV       ] -> Fill( lPt, lDCAzNegToPrimVertex );
                lHistSigPlusCenterBgAll[ DCAZPOSTOPV       ] -> Fill( lPt, lDCAzPosToPrimVertex );
                lHistSigPlusCenterBgAll[ DCAZBACHTOPV      ] -> Fill( lPt, lDCAzBachToPrimVertex );
                lHistSigPlusCenterBgAll[ DCABACHTOBARYON   ] -> Fill( lPt, lDCABachToBaryon );
                lHistSigPlusCenterBgAll[ V0PA              ] -> Fill( lPt, TMath::ACos( lV0CosinePointingAngle ) );
                lHistSigPlusCenterBgAll[ CASCPA            ] -> Fill( lPt, TMath::ACos( lCascCosinePointingAngle ) );
                lHistSigPlusCenterBgAll[ BBPA              ] -> Fill( lPt, TMath::ACos( lBBCosPA ) );
                lHistSigPlusCenterBgAll[ V0MASS            ] -> Fill( lPt, lV0Mass );
                lHistSigPlusCenterBgAll[ PROPERLIFETIME    ] -> Fill( lPt, lDistOverTotMom*lParticleMass );
                lHistSigPlusCenterBgAll[ COMPETINGSPECIES  ] -> Fill( lPt, lCompetingParticleMass );
                lHistSigPlusCenterBgAll[ TPCNEGDEDX        ] -> Fill( lNegP,  lNegdEdx );
                lHistSigPlusCenterBgAll[ TPCPOSDEDX        ] -> Fill( lPosP,  lPosdEdx );
                lHistSigPlusCenterBgAll[ TPCBACHDEDX       ] -> Fill( lBachP, lBachdEdx );
                lHistSigPlusCenterBgAll[ TPCNEGPIDNSIGMAS  ] -> Fill( lPt, lNegNSigmas );
                lHistSigPlusCenterBgAll[ TPCPOSPIDNSIGMAS  ] -> Fill( lPt, lPosNSigmas );
                lHistSigPlusCenterBgAll[ TPCBACHPIDNSIGMAS ] -> Fill( lPt, lBachNSigmas );
                lHistSigPlusCenterBgAll[ TPCNCLUSTERS      ] -> Fill( lPt, lLeastNbrClusters );
                lHistSigPlusCenterBgAll[ MAXCHI2PERCLUSTER ] -> Fill( lPt, lMaxChi2PerCluster );
                lHistSigPlusCenterBgAll[ MINTRACKLENGTH    ] -> Fill( lPt, lMinTrackLength );
                lHistSigPlusCenterBgAll[ NEGDAUGHTERETA    ] -> Fill( lPt, lNegEta );
                lHistSigPlusCenterBgAll[ POSDAUGHTERETA    ] -> Fill( lPt, lPosEta );
                lHistSigPlusCenterBgAll[ BACHDAUGHTERETA   ] -> Fill( lPt, lBachEta );
                lHistSigPlusCenterBgAll[ NEGINNERP         ] -> Fill( lPt, lNegInnerP );
                lHistSigPlusCenterBgAll[ POSINNERP         ] -> Fill( lPt, lPosInnerP );
                lHistSigPlusCenterBgAll[ BACHINNERP        ] -> Fill( lPt, lBachInnerP );
                lHistSigPlusCenterBgAll[ NEGTOFEXPTDIFF    ] -> Fill( lPt, lNegTOFExpTDiff );
                lHistSigPlusCenterBgAll[ POSTOFEXPTDIFF    ] -> Fill( lPt, lPosTOFExpTDiff );
                lHistSigPlusCenterBgAll[ BACHTOFEXPTDIFF   ] -> Fill( lPt, lBachTOFExpTDiff );
                lHistSigPlusCenterBgAll[ NEGTOFSIGNAL      ] -> Fill( lPt, lNegTOFSignal );
                lHistSigPlusCenterBgAll[ POSTOFSIGNAL      ] -> Fill( lPt, lPosTOFSignal );
                lHistSigPlusCenterBgAll[ BACHTOFSIGNAL     ] -> Fill( lPt, lBachTOFSignal );
            }

            //--- Fill histograms for left+right background candidates
            if( (lInvariantMass>lLeftBgLeftLimit[lWeAreAtBin]  && lInvariantMass<lLeftBgRightLimit[lWeAreAtBin]  ) ||
                (lInvariantMass>lRightBgLeftLimit[lWeAreAtBin] && lInvariantMass<lRightBgRightLimit[lWeAreAtBin] ) ) {
                lHistLeftPlusRightBgAll[ V0RADIUS          ] -> Fill( lPt, lV0Radius );
                lHistLeftPlusRightBgAll[ CASCRADIUS        ] -> Fill( lPt, lCascRadius );
                lHistLeftPlusRightBgAll[ DCANEGTOPV        ] -> Fill( lPt, lDcaNegToPrimVertex );
                lHistLeftPlusRightBgAll[ DCAPOSTOPV        ] -> Fill( lPt, lDcaPosToPrimVertex );
                lHistLeftPlusRightBgAll[ DCABACHTOPV       ] -> Fill( lPt, lDcaBachToPrimVertex );
                lHistLeftPlusRightBgAll[ DCAV0TOPV         ] -> Fill( lPt, lDcaV0ToPV );
                lHistLeftPlusRightBgAll[ DCACASCTOPV       ] -> Fill( lPt, lDCACascToPV );
                lHistLeftPlusRightBgAll[ DCAXYCASCTOPV     ] -> Fill( lPt, lDCAxyCascToPV );
                lHistLeftPlusRightBgAll[ DCAZCASCTOPV      ] -> Fill( lPt, lDCAzCascToPV );
                lHistLeftPlusRightBgAll[ DCAV0DAUGHTERS    ] -> Fill( lPt, lDcaV0Daughters );
                lHistLeftPlusRightBgAll[ DCACASCDAUGHTERS  ] -> Fill( lPt, lDcaCascDaughters );
                lHistLeftPlusRightBgAll[ DCAZNEGTOPV       ] -> Fill( lPt, lDCAzNegToPrimVertex );
                lHistLeftPlusRightBgAll[ DCAZPOSTOPV       ] -> Fill( lPt, lDCAzPosToPrimVertex );
                lHistLeftPlusRightBgAll[ DCAZBACHTOPV      ] -> Fill( lPt, lDCAzBachToPrimVertex );
                lHistLeftPlusRightBgAll[ DCABACHTOBARYON   ] -> Fill( lPt, lDCABachToBaryon );
                lHistLeftPlusRightBgAll[ V0PA              ] -> Fill( lPt, TMath::ACos( lV0CosinePointingAngle ) );
                lHistLeftPlusRightBgAll[ CASCPA            ] -> Fill( lPt, TMath::ACos( lCascCosinePointingAngle ) );
                lHistLeftPlusRightBgAll[ BBPA              ] -> Fill( lPt, TMath::ACos( lBBCosPA ) );
                lHistLeftPlusRightBgAll[ V0MASS            ] -> Fill( lPt, lV0Mass );
                lHistLeftPlusRightBgAll[ PROPERLIFETIME    ] -> Fill( lPt, lDistOverTotMom*lParticleMass );
                lHistLeftPlusRightBgAll[ COMPETINGSPECIES  ] -> Fill( lPt, lCompetingParticleMass );
                lHistLeftPlusRightBgAll[ TPCNEGDEDX        ] -> Fill( lNegP,  lNegdEdx );
                lHistLeftPlusRightBgAll[ TPCPOSDEDX        ] -> Fill( lPosP,  lPosdEdx );
                lHistLeftPlusRightBgAll[ TPCBACHDEDX       ] -> Fill( lBachP, lBachdEdx );
                lHistLeftPlusRightBgAll[ TPCNEGPIDNSIGMAS  ] -> Fill( lPt, lNegNSigmas );
                lHistLeftPlusRightBgAll[ TPCPOSPIDNSIGMAS  ] -> Fill( lPt, lPosNSigmas );
                lHistLeftPlusRightBgAll[ TPCBACHPIDNSIGMAS ] -> Fill( lPt, lBachNSigmas );
                lHistLeftPlusRightBgAll[ TPCNCLUSTERS      ] -> Fill( lPt, lLeastNbrClusters );
                lHistLeftPlusRightBgAll[ MAXCHI2PERCLUSTER ] -> Fill( lPt, lMaxChi2PerCluster );
                lHistLeftPlusRightBgAll[ MINTRACKLENGTH    ] -> Fill( lPt, lMinTrackLength );
                lHistLeftPlusRightBgAll[ NEGDAUGHTERETA    ] -> Fill( lPt, lNegEta );
                lHistLeftPlusRightBgAll[ POSDAUGHTERETA    ] -> Fill( lPt, lPosEta );
                lHistLeftPlusRightBgAll[ BACHDAUGHTERETA   ] -> Fill( lPt, lBachEta );
                lHistLeftPlusRightBgAll[ NEGINNERP         ] -> Fill( lPt, lNegInnerP );
                lHistLeftPlusRightBgAll[ POSINNERP         ] -> Fill( lPt, lPosInnerP );
                lHistLeftPlusRightBgAll[ BACHINNERP        ] -> Fill( lPt, lBachInnerP );
                lHistLeftPlusRightBgAll[ NEGTOFEXPTDIFF    ] -> Fill( lPt, lNegTOFExpTDiff );
                lHistLeftPlusRightBgAll[ POSTOFEXPTDIFF    ] -> Fill( lPt, lPosTOFExpTDiff );
                lHistLeftPlusRightBgAll[ BACHTOFEXPTDIFF   ] -> Fill( lPt, lBachTOFExpTDiff );
                lHistLeftPlusRightBgAll[ NEGTOFSIGNAL    ] -> Fill( lPt, lNegTOFSignal );
                lHistLeftPlusRightBgAll[ POSTOFSIGNAL    ] -> Fill( lPt, lPosTOFSignal );
                lHistLeftPlusRightBgAll[ BACHTOFSIGNAL   ] -> Fill( lPt, lBachTOFSignal );
            }
        }

        //Now check validity
        if( lRap<fRapidityBoundaryUpper && lRap>fRapidityBoundaryLower &&
            (//charge condition (x-check)
              (fWhichParticle == "XiMinus"    && lCharge == -1) ||
              (fWhichParticle == "XiPlus"     && lCharge ==  1) ||
              (fWhichParticle == "OmegaMinus" && lCharge == -1) ||
              (fWhichParticle == "OmegaPlus"  && lCharge ==  1)
            ) &&
            TMath::Abs(lNegEta)       <= GetCutDaughterEta      ( lPt ) &&
            TMath::Abs(lPosEta)       <= GetCutDaughterEta      ( lPt ) &&
            TMath::Abs(lBachEta)      <= GetCutDaughterEta      ( lPt ) &&
            //Topological Selections
            lV0Radius                 >= GetCutV0Radius         ( lPt ) &&
            lCascRadius               >= GetCutCascRadius       ( lPt ) &&
            TMath::Abs(lV0Mass-1.116) <= GetCutV0Mass           ( lPt ) &&
            lV0CosinePointingAngle    >= GetCutV0CosPA          ( lPt ) &&
            lCascCosinePointingAngle  >= GetCutCascCosPA        ( lPt ) &&
            lDcaNegToPrimVertex       >= GetCutDCANegToPV       ( lPt ) &&
            lDcaPosToPrimVertex       >= GetCutDCAPosToPV       ( lPt ) &&
            lDcaBachToPrimVertex      >= GetCutDCABachToPV      ( lPt ) &&
            lDcaV0Daughters           <= GetCutDCAV0Daughters   ( lPt ) &&
            lDcaCascDaughters         <= GetCutDCACascDaughters ( lPt ) &&
            lDcaV0ToPV                >= GetCutDCAV0ToPV        ( lPt ) &&
            lDCACascToPV              <= GetCutDCACascToPV      ( lPt ) &&
            lDCAxyCascToPV            <= GetCutDCAxyCascToPV    ( lPt ) &&
            lDCAzCascToPV             <= GetCutDCAzCascToPV     ( lPt ) &&
            lParticleMass*lDistOverTotMom     <= GetCutProperLifetime   ( lPt ) &&
            lDCAzNegToPrimVertex      <= GetCutDCAzNegToPV      ( lPt ) &&
            lDCAzPosToPrimVertex      <= GetCutDCAzPosToPV      ( lPt ) &&
            lDCAzBachToPrimVertex     <= GetCutDCAzBachToPV     ( lPt ) &&
            lDCABachToBaryon          >= GetCutDCABachToBaryon  ( lPt ) &&
            lBBCosPA                  <= GetCutBBCosPA          ( lPt ) &&
            lMinTrackLength           >= GetCutMinTrackLength   ( lPt ) &&

            //Competing Species Rejection (only for Omegas)
            TMath::Abs( lCompetingParticleMass - 1.32171 ) >= GetCutCompetingSpecies( lPt ) &&

            //Causality Cut
            TMath::Abs( lV0Radius - lCascRadius ) >= GetCutCausality( lPt ) &&

            //Nclusters cut
            lLeastNbrClusters >= GetCutLeastNumberOfClusters( lPt ) &&

            ( //official response code
              ( fWhichParticle == "XiMinus"
                && TMath::Abs(lNSigmasNegPion)   <= GetCutTPCPIDNSigmas( lPt )
                && TMath::Abs(lNSigmasPosProton) <= GetCutTPCPIDNSigmas( lPt )
                && TMath::Abs(lNSigmasBachPion)  <= GetCutTPCPIDNSigmas( lPt )) ||
              ( fWhichParticle == "XiPlus"
                && TMath::Abs(lNSigmasPosPion)   <= GetCutTPCPIDNSigmas( lPt )
                && TMath::Abs(lNSigmasNegProton) <= GetCutTPCPIDNSigmas( lPt )
                && TMath::Abs(lNSigmasBachPion)  <= GetCutTPCPIDNSigmas( lPt )) ||
              ( fWhichParticle == "OmegaMinus"
                && TMath::Abs(lNSigmasNegPion)   <= GetCutTPCPIDNSigmas( lPt )
                && TMath::Abs(lNSigmasPosProton) <= GetCutTPCPIDNSigmas( lPt )
                && TMath::Abs(lNSigmasBachKaon)  <= GetCutTPCPIDNSigmas( lPt )) ||
              ( fWhichParticle == "OmegaPlus"
                && TMath::Abs(lNSigmasPosPion)   <= GetCutTPCPIDNSigmas( lPt )
                && TMath::Abs(lNSigmasNegProton) <= GetCutTPCPIDNSigmas( lPt )
                && TMath::Abs(lNSigmasBachKaon)  <= GetCutTPCPIDNSigmas( lPt ))
            ) &&

            ( //TOF match 
              CheckTOFmatchBach( lPt, lBachTOFExpTDiff ) &&
              CheckTOFmatchNeg(  lPt, lNegTOFExpTDiff  ) &&
              CheckTOFmatchPos(  lPt, lPosTOFExpTDiff  ) &&
              CheckTOFmatchOne(  lPt, lBachTOFExpTDiff, lNegTOFExpTDiff, lPosTOFExpTDiff )
            ) &&

            ( //ITS refit 
              CheckITSrefitBach( lPt, lBachTrackStatus ) &&
              CheckITSrefitNeg(  lPt, lNegTrackStatus  ) &&
              CheckITSrefitPos(  lPt, lPosTrackStatus  ) &&
              CheckITSrefitOne(  lPt, lBachTrackStatus, lNegTrackStatus, lPosTrackStatus )
            ) &&
                
            ( //ITSTOF 
              CheckITSTOFBach( lPt, lBachTrackStatus, lBachTOFExpTDiff ) &&
              CheckITSTOFNeg(  lPt, lNegTrackStatus,  lNegTOFExpTDiff  ) &&
              CheckITSTOFPos(  lPt, lPosTrackStatus,  lPosTOFExpTDiff  ) &&
              CheckITSTOFOne(  lPt, lBachTrackStatus, lNegTrackStatus, lPosTrackStatus, lBachTOFExpTDiff, lNegTOFExpTDiff, lPosTOFExpTDiff )
            )

          ) { // Start Entry Loop
            lHistoSelectedCasc[lWeAreAtBin]->Fill(lInvariantMass); //fill with specific inv mass

            //--- Peak Region
            if( lInvariantMass>lPeakLeftLimit[lWeAreAtBin] && lInvariantMass<lPeakRightLimit[lWeAreAtBin] ) {
                lSigPlusCenterBgCasc[lWeAreAtBin]++;

                //--- Fill histograms for signal+center background candidates
                if(fSaveVarHistosSwitch) {
                    lHistSigPlusCenterBgSel[ V0RADIUS          ] -> Fill( lPt, lV0Radius );
                    lHistSigPlusCenterBgSel[ CASCRADIUS        ] -> Fill( lPt, lCascRadius );
                    lHistSigPlusCenterBgSel[ DCANEGTOPV        ] -> Fill( lPt, lDcaNegToPrimVertex );
                    lHistSigPlusCenterBgSel[ DCAPOSTOPV        ] -> Fill( lPt, lDcaPosToPrimVertex );
                    lHistSigPlusCenterBgSel[ DCABACHTOPV       ] -> Fill( lPt, lDcaBachToPrimVertex );
                    lHistSigPlusCenterBgSel[ DCAV0TOPV         ] -> Fill( lPt, lDcaV0ToPV );
                    lHistSigPlusCenterBgSel[ DCACASCTOPV       ] -> Fill( lPt, lDCACascToPV );
                    lHistSigPlusCenterBgSel[ DCAXYCASCTOPV     ] -> Fill( lPt, lDCAxyCascToPV );
                    lHistSigPlusCenterBgSel[ DCAZCASCTOPV      ] -> Fill( lPt, lDCAzCascToPV );
                    lHistSigPlusCenterBgSel[ DCAV0DAUGHTERS    ] -> Fill( lPt, lDcaV0Daughters );
                    lHistSigPlusCenterBgSel[ DCACASCDAUGHTERS  ] -> Fill( lPt, lDcaCascDaughters );
                    lHistSigPlusCenterBgSel[ DCAZNEGTOPV       ] -> Fill( lPt, lDCAzNegToPrimVertex );
                    lHistSigPlusCenterBgSel[ DCAZPOSTOPV       ] -> Fill( lPt, lDCAzPosToPrimVertex );
                    lHistSigPlusCenterBgSel[ DCAZBACHTOPV      ] -> Fill( lPt, lDCAzBachToPrimVertex );
                    lHistSigPlusCenterBgSel[ DCABACHTOBARYON   ] -> Fill( lPt, lDCABachToBaryon );
                    lHistSigPlusCenterBgSel[ V0PA              ] -> Fill( lPt, TMath::ACos( lV0CosinePointingAngle ) );
                    lHistSigPlusCenterBgSel[ CASCPA            ] -> Fill( lPt, TMath::ACos( lCascCosinePointingAngle ) );
                    lHistSigPlusCenterBgSel[ BBPA              ] -> Fill( lPt, TMath::ACos( lBBCosPA ) );
                    lHistSigPlusCenterBgSel[ V0MASS            ] -> Fill( lPt, lV0Mass );
                    lHistSigPlusCenterBgSel[ PROPERLIFETIME    ] -> Fill( lPt, lDistOverTotMom*lParticleMass );
                    lHistSigPlusCenterBgSel[ COMPETINGSPECIES  ] -> Fill( lPt, lCompetingParticleMass );
                    lHistSigPlusCenterBgSel[ TPCNEGDEDX        ] -> Fill( lNegP,  lNegdEdx );
                    lHistSigPlusCenterBgSel[ TPCPOSDEDX        ] -> Fill( lPosP,  lPosdEdx );
                    lHistSigPlusCenterBgSel[ TPCBACHDEDX       ] -> Fill( lBachP, lBachdEdx );
                    lHistSigPlusCenterBgSel[ TPCNEGPIDNSIGMAS  ] -> Fill( lPt, lNegNSigmas );
                    lHistSigPlusCenterBgSel[ TPCPOSPIDNSIGMAS  ] -> Fill( lPt, lPosNSigmas );
                    lHistSigPlusCenterBgSel[ TPCBACHPIDNSIGMAS ] -> Fill( lPt, lBachNSigmas );
                    lHistSigPlusCenterBgSel[ TPCNCLUSTERS      ] -> Fill( lPt, lLeastNbrClusters );
                    lHistSigPlusCenterBgSel[ MAXCHI2PERCLUSTER ] -> Fill( lPt, lMaxChi2PerCluster );
                    lHistSigPlusCenterBgSel[ MINTRACKLENGTH    ] -> Fill( lPt, lMinTrackLength );
                    lHistSigPlusCenterBgSel[ NEGDAUGHTERETA    ] -> Fill( lPt, lNegEta );
                    lHistSigPlusCenterBgSel[ POSDAUGHTERETA    ] -> Fill( lPt, lPosEta );
                    lHistSigPlusCenterBgSel[ BACHDAUGHTERETA   ] -> Fill( lPt, lBachEta );
                    lHistSigPlusCenterBgSel[ NEGINNERP         ] -> Fill( lPt, lNegInnerP );
                    lHistSigPlusCenterBgSel[ POSINNERP         ] -> Fill( lPt, lPosInnerP );
                    lHistSigPlusCenterBgSel[ BACHINNERP        ] -> Fill( lPt, lBachInnerP );
                    lHistSigPlusCenterBgSel[ NEGTOFEXPTDIFF    ] -> Fill( lPt, lNegTOFExpTDiff );
                    lHistSigPlusCenterBgSel[ POSTOFEXPTDIFF    ] -> Fill( lPt, lPosTOFExpTDiff );
                    lHistSigPlusCenterBgSel[ BACHTOFEXPTDIFF   ] -> Fill( lPt, lBachTOFExpTDiff );
                    lHistSigPlusCenterBgSel[ NEGTOFSIGNAL      ] -> Fill( lPt, lNegTOFSignal );
                    lHistSigPlusCenterBgSel[ POSTOFSIGNAL      ] -> Fill( lPt, lPosTOFSignal );
                    lHistSigPlusCenterBgSel[ BACHTOFSIGNAL     ] -> Fill( lPt, lBachTOFSignal );
                }
            }

            //--- Left and Right Background Region
            if( (lInvariantMass>lLeftBgLeftLimit[lWeAreAtBin]  && lInvariantMass<lLeftBgRightLimit[lWeAreAtBin]  ) ||
                (lInvariantMass>lRightBgLeftLimit[lWeAreAtBin] && lInvariantMass<lRightBgRightLimit[lWeAreAtBin] ) ) {
                lLeftPlusRightBgCasc[lWeAreAtBin]++;

                //--- Fill histograms for left+right background candidates
                if(fSaveVarHistosSwitch) {
                    lHistLeftPlusRightBgSel[ V0RADIUS          ] -> Fill( lPt, lV0Radius );
                    lHistLeftPlusRightBgSel[ CASCRADIUS        ] -> Fill( lPt, lCascRadius );
                    lHistLeftPlusRightBgSel[ DCANEGTOPV        ] -> Fill( lPt, lDcaNegToPrimVertex );
                    lHistLeftPlusRightBgSel[ DCAPOSTOPV        ] -> Fill( lPt, lDcaPosToPrimVertex );
                    lHistLeftPlusRightBgSel[ DCABACHTOPV       ] -> Fill( lPt, lDcaBachToPrimVertex );
                    lHistLeftPlusRightBgSel[ DCAV0TOPV         ] -> Fill( lPt, lDcaV0ToPV );
                    lHistLeftPlusRightBgSel[ DCACASCTOPV       ] -> Fill( lPt, lDCACascToPV );
                    lHistLeftPlusRightBgSel[ DCAXYCASCTOPV     ] -> Fill( lPt, lDCAxyCascToPV );
                    lHistLeftPlusRightBgSel[ DCAZCASCTOPV      ] -> Fill( lPt, lDCAzCascToPV );
                    lHistLeftPlusRightBgSel[ DCAV0DAUGHTERS    ] -> Fill( lPt, lDcaV0Daughters );
                    lHistLeftPlusRightBgSel[ DCACASCDAUGHTERS  ] -> Fill( lPt, lDcaCascDaughters );
                    lHistLeftPlusRightBgSel[ DCAZNEGTOPV       ] -> Fill( lPt, lDCAzNegToPrimVertex );
                    lHistLeftPlusRightBgSel[ DCAZPOSTOPV       ] -> Fill( lPt, lDCAzPosToPrimVertex );
                    lHistLeftPlusRightBgSel[ DCAZBACHTOPV      ] -> Fill( lPt, lDCAzBachToPrimVertex );
                    lHistLeftPlusRightBgSel[ DCABACHTOBARYON   ] -> Fill( lPt, lDCABachToBaryon );
                    lHistLeftPlusRightBgSel[ V0PA              ] -> Fill( lPt, TMath::ACos( lV0CosinePointingAngle ) );
                    lHistLeftPlusRightBgSel[ CASCPA            ] -> Fill( lPt, TMath::ACos( lCascCosinePointingAngle ) );
                    lHistLeftPlusRightBgSel[ BBPA              ] -> Fill( lPt, TMath::ACos( lBBCosPA ) );
                    lHistLeftPlusRightBgSel[ V0MASS            ] -> Fill( lPt, lV0Mass );
                    lHistLeftPlusRightBgSel[ PROPERLIFETIME    ] -> Fill( lPt, lDistOverTotMom*lParticleMass );
                    lHistLeftPlusRightBgSel[ COMPETINGSPECIES  ] -> Fill( lPt, lCompetingParticleMass );
                    lHistLeftPlusRightBgSel[ TPCNEGDEDX        ] -> Fill( lNegP,  lNegdEdx );
                    lHistLeftPlusRightBgSel[ TPCPOSDEDX        ] -> Fill( lPosP,  lPosdEdx );
                    lHistLeftPlusRightBgSel[ TPCBACHDEDX       ] -> Fill( lBachP, lBachdEdx );
                    lHistLeftPlusRightBgSel[ TPCNEGPIDNSIGMAS  ] -> Fill( lPt, lNegNSigmas );
                    lHistLeftPlusRightBgSel[ TPCPOSPIDNSIGMAS  ] -> Fill( lPt, lPosNSigmas );
                    lHistLeftPlusRightBgSel[ TPCBACHPIDNSIGMAS ] -> Fill( lPt, lBachNSigmas );
                    lHistLeftPlusRightBgSel[ TPCNCLUSTERS      ] -> Fill( lPt, lLeastNbrClusters );
                    lHistLeftPlusRightBgSel[ MAXCHI2PERCLUSTER ] -> Fill( lPt, lMaxChi2PerCluster );
                    lHistLeftPlusRightBgSel[ MINTRACKLENGTH    ] -> Fill( lPt, lMinTrackLength );
                    lHistLeftPlusRightBgSel[ NEGDAUGHTERETA    ] -> Fill( lPt, lNegEta );
                    lHistLeftPlusRightBgSel[ POSDAUGHTERETA    ] -> Fill( lPt, lPosEta );
                    lHistLeftPlusRightBgSel[ BACHDAUGHTERETA   ] -> Fill( lPt, lBachEta );
                    lHistLeftPlusRightBgSel[ NEGINNERP         ] -> Fill( lPt, lNegInnerP );
                    lHistLeftPlusRightBgSel[ POSINNERP         ] -> Fill( lPt, lPosInnerP );
                    lHistLeftPlusRightBgSel[ BACHINNERP        ] -> Fill( lPt, lBachInnerP );
                    lHistLeftPlusRightBgSel[ NEGTOFEXPTDIFF    ] -> Fill( lPt, lNegTOFExpTDiff );
                    lHistLeftPlusRightBgSel[ POSTOFEXPTDIFF    ] -> Fill( lPt, lPosTOFExpTDiff );
                    lHistLeftPlusRightBgSel[ BACHTOFEXPTDIFF   ] -> Fill( lPt, lBachTOFExpTDiff );
                    lHistLeftPlusRightBgSel[ NEGTOFSIGNAL      ] -> Fill( lPt, lNegTOFSignal );
                    lHistLeftPlusRightBgSel[ POSTOFSIGNAL      ] -> Fill( lPt, lPosTOFSignal );
                    lHistLeftPlusRightBgSel[ BACHTOFSIGNAL     ] -> Fill( lPt, lBachTOFSignal );
                }
            }

        } // End Entry Loop
    }
    cout<<"--------------- Loop Completed -------------------------"<<endl;
    cout<<endl;

    cout<<"--------------- Memory Cleanup -------------------------"<<endl;
    if ( clist ) {
        clist->Delete();
        delete clist;
    }

    file->Close("R");
    file->Delete();
    delete file;
    cout<<endl;

    TF1 *lfitNoise[50];
    TF1 *lSampleNoise[50];
    char lFitNameOne[50];

    if( fFitBackgroundSwitch ){
        cout<<"--------------- Backgrounds: fitting with linear -------"<<endl;
        for(long i=0; i<fptbinnumb; i++){
            //Define Function to Fit Background
            sprintf(lFitNameOne,"lfitNoise%i",((int)(i)));
            //cout<<"creating fitnoise, named "<<lFitNameOne<<endl;
            lfitNoise[i] = new TF1(lFitNameOne, this, &AliCascadeModule::MyBgPol1,
                    RoundToThousandth ( lLeftBgLeftLimit[i]   ),
                    RoundToThousandth ( lRightBgRightLimit[i] ), 4 , "AliCascadeModule", "MyBgPol1");
            lfitNoise[i] -> FixParameter(2,RoundToThousandth ( lLeftBgRightLimit[i] ) );
            lfitNoise[i] -> FixParameter(3,RoundToThousandth ( lRightBgLeftLimit[i] ) );
            lfitNoise[i] -> SetParameter(0,lLeftPlusRightBgCasc[i] * lHistoMBCasc[i]->GetBinWidth(5) / (lRightBgLeftLimit[i]-lLeftBgRightLimit[i] + 1e-6 ) );
            lfitNoise[i] -> SetParameter(1,0);
            cout<<"Guessed Parameter 0 for "<<lFitNameOne<<" to be "<<lLeftPlusRightBgCasc[i] * lHistoMBCasc[i]->GetBinWidth(5) / (lRightBgLeftLimit[i]-lLeftBgRightLimit[i] + 1e-6 )<<endl;
            sprintf(lFitNameOne,"lSampleNoise%i",((int)(i)));

            //Define Function to Sample Background
            //cout<<"creating sample "<<i<<endl;
            lSampleNoise[i] = new TF1(lFitNameOne, this, &AliCascadeModule::MyBgPolToEval1,
                    RoundToThousandth ( lLeftBgLeftLimit[i]   ),
                    RoundToThousandth ( lRightBgRightLimit[i] ), 2 , "AliCacadeModule", "MyBgPolToEval1");
        }
        for(long i=0; i<fptbinnumb; i++){
            //cout<<"Fitting function for bin "<<i<<", get name = "<<lfitNoise[i]->GetName()<<endl;
            sprintf(lFitNameOne,"lfitNoise%i",((int)(i)));
            lHistoMBCasc[i] -> Fit(lFitNameOne,"LLrie+0");
            lSampleNoise[i]->SetParameter(0, lfitNoise[i]->GetParameter(0) );
            lSampleNoise[i]->SetParameter(1, lfitNoise[i]->GetParameter(1) );
        }
        for(long i=0; i<fptbinnumb; i++){
            cout<<"Overriding Background info: Was "<<lLeftPlusRightBgCasc[i]<<", is now "<<lSampleNoise[i]->Integral( lPeakLeftLimit[i], lPeakRightLimit[i] )/lHistoMBCasc[i]->GetBinWidth(5)<<endl;
            lLeftPlusRightBgCasc[i] = lSampleNoise[i]->Integral( lPeakLeftLimit[i], lPeakRightLimit[i] )/lHistoMBCasc[i]->GetBinWidth(5);
        }
        cout<<"--------------- Fitting Finished! ----------------------"<<endl;
    }

    //=============================================================
    // Compute Signal + Sig to noise
    for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++) {

        //Signal computation
        lSigRealCasc   [ipoint] = lSigPlusCenterBgCasc[ipoint] - lLeftPlusRightBgCasc[ipoint];
        lSigErrRealCasc[ipoint] = TMath::Sqrt( lSigPlusCenterBgCasc[ipoint] + lLeftPlusRightBgCasc[ipoint] );
        //Error Formula: Equivalent to Sqrt(S+B+B) = Sqrt(S+2B)

        //Attribute Raw Counts
        fHistPtRaw->SetBinContent(ipoint+1, lSigRealCasc   [ipoint] / ((fptbinlimits[ipoint+1] - fptbinlimits[ipoint]) * lNEvents * (fRapidityBoundaryUpper - fRapidityBoundaryLower)) );
        fHistPtRaw->SetBinError  (ipoint+1, lSigErrRealCasc[ipoint] / ((fptbinlimits[ipoint+1] - fptbinlimits[ipoint]) * lNEvents * (fRapidityBoundaryUpper - fRapidityBoundaryLower)) );

        //Attribute Raw Counts
        fHistPtSignal->SetBinContent(ipoint+1, lSigRealCasc   [ipoint] );
        fHistPtSignal->SetBinError  (ipoint+1, lSigErrRealCasc[ipoint] );

        //Signal-to-noise computation
        if( lLeftPlusRightBgCasc[ipoint] != 0 ) {
            fHistSigToNoise->SetBinContent(ipoint+1, (lSigPlusCenterBgCasc[ipoint] - lLeftPlusRightBgCasc[ipoint]) / lLeftPlusRightBgCasc[ipoint] );
        } else {
            fHistSigToNoise->SetBinContent(ipoint+1, -1) ; //-1 means: no background
        }
    }
    //=============================================================

    cout<<"--------------- Extracted Signal -----------------------"<<endl;
    for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++){
        cout<<"---> ["<<fptbinlimits[ipoint]<<" - "<<fptbinlimits[ipoint+1]<<" GeV/c]\tSignal: "<<lSigRealCasc[ipoint]<<" +/- "<<lSigErrRealCasc[ipoint]<<endl;
    }
    cout<<"--------------------------------------------------------"<<endl;

    //=============================================================
    // Preparations for MC loop
    //Extra info available only in MC
    Int_t lPID = 0;
    Int_t lPIDBachelor = 0;
    Int_t lPIDNegative = 0;
    Int_t lPIDPositive = 0;
    Int_t lPrimaryStatus = 0;
    Float_t lRapMC = 0;
    //=============================================================
    cout<<endl;

    //Prepare containers for Geant-Fluka correction
    TH1F *lProtonMomentum[100];
    char lNameOne[100];
    for(Int_t ibin=0; ibin<100; ibin++) {
        sprintf(lNameOne, "lProtonMomentumBin%i", ibin);
        lProtonMomentum[ibin] = new TH1F(lNameOne, "", 800, 0., 20.);
        if(fWhichParticle == "XiMinus")    bindescription = "#Xi^{-}, bin #";
        if(fWhichParticle == "XiPlus")     bindescription = "#bar{#Xi}^{+}, bin #";
        if(fWhichParticle == "OmegaMinus") bindescription = "#Omega^{-}, bin #";
        if(fWhichParticle == "OmegaPlus")  bindescription = "#bar{#Omega}^{+}, bin #";
        bindescription.Append(IntToString( ibin ));
        if( ibin < fptbinnumb ) {
            bindescription.Append(", ");
            bindescription.Append(DoubleToString(fptbinlimits[ ibin ]));
            bindescription.Append("-");
            bindescription.Append(DoubleToString(fptbinlimits[ ibin+1 ]));
            bindescription.Append("GeV/c");
        }
        lProtonMomentum[ibin]->SetTitle(bindescription);
    }

    //MC File Acquisition=============================================
    cout<<"--------------- Opening MC file ------------------------"<<endl;
    TFile* fileMC = TFile::Open(fMCDataFile, "READ");
    TList* clistMC      = (TList*)fileMC->Get("PWGLF_StrVsMult_MC/cList");
    TTree *lTreeMC      = (TTree*)fileMC->Get("PWGLF_StrVsMult_MC/fTreeCascade");
    TTree* lTreeEventMC = (TTree*)fileMC->Get("PWGLF_StrVsMult_MC/fTreeEvent");
    TH1F* fHistMultiplicityBeforeTrigSelMC;

    Double_t lLoMultBoundMC = fLoMultBound;
    Double_t lHiMultBoundMC = fHiMultBound;
    if(fUseIntegratedEfficiency) {
        cout << " Using integrated multiplicity for efficiency computation." << endl;
        lLoMultBoundMC =   0.;
        lHiMultBoundMC = 100.;
    }

    lTreeEventMC->SetBranchAddress("fCentrality",   &fCentrality);
    lTreeEventMC->SetBranchAddress("fMVPileupFlag", &fMVPileupFlag);

    cout<<"--------------------------------------------------------"<<endl;
    cout<<" Will now loop over events, please wait..."<<endl;
    Long_t lNEventsMC = 0;
    for(Long_t iEv = 0; iEv<lTreeEventMC->GetEntries(); iEv++) {
        lTreeEventMC->GetEntry(iEv);
        if( iEv % ( lTreeEventMC->GetEntries() / 10 ) == 0 ) cout<<" At Event "<<iEv<<" out of "<<lTreeEventMC->GetEntries()<<endl;
        // check MV Pileup rejection
        if( fMVPileupSwitch && !fMVPileupFlag ) continue;
        // Multiplicity Switch
        if( fPerformMultiplicityStudy == kTRUE &&  //inside mult bin
            fCentrality>lLoMultBoundMC &&
            fCentrality<lHiMultBoundMC
          ) lNEventsMC++;
    }
    cout<<" Number of events (MC), this multiplicity....: "<<lNEventsMC <<endl;
    cout<<"--------------------------------------------------------"<<endl;

    //================================================================

    //Linking to Tree=================================================
    //--- Base Variables ----------------------------------------------
    lTreeMC->SetBranchAddress("fTreeCascVarCharge"  , &lCharge  );
    lTreeMC->SetBranchAddress("fTreeCascVarPosEta"  , &lPosEta  );
    lTreeMC->SetBranchAddress("fTreeCascVarNegEta"  , &lNegEta  );
    lTreeMC->SetBranchAddress("fTreeCascVarBachEta" , &lBachEta );
    lTreeMC->SetBranchAddress("fTreeCascVarPt"      , &lPt      );
    lTreeMC->SetBranchAddress("fTreeCascVarPtMC"    , &lPtMC    );
    if ( fWhichParticle == "XiMinus" )  lTreeMC->SetBranchAddress("fTreeCascVarMassAsXi", &lInvariantMass);
    if ( fWhichParticle == "XiPlus"  )  lTreeMC->SetBranchAddress("fTreeCascVarMassAsXi", &lInvariantMass);
    if ( fWhichParticle == "OmegaMinus" || fWhichParticle == "OmegaPlus" ) {
        lTreeMC->SetBranchAddress("fTreeCascVarMassAsOmega", &lInvariantMass);
        lTreeMC->SetBranchAddress("fTreeCascVarMassAsXi"   , &lCompetingParticleMass); //For Competing Rejection
    }
    if ( fWhichParticle == "XiMinus"    || fWhichParticle == "XiPlus"    )
        lTreeMC->SetBranchAddress("fTreeCascVarRapXi"   , &lRap);
    if ( fWhichParticle == "OmegaMinus" || fWhichParticle == "OmegaPlus" )
        lTreeMC->SetBranchAddress("fTreeCascVarRapOmega", &lRap);
    lTreeMC->SetBranchAddress("fTreeCascVarRapMC", &lRapMC);
    lTreeMC->SetBranchAddress("fTreeCascVarLeastNbrClusters", &lLeastNbrClusters);
    lTreeMC->SetBranchAddress("fTreeCascVarPosPx"           , &lPosPx);
    lTreeMC->SetBranchAddress("fTreeCascVarPosPy"           , &lPosPy);
    lTreeMC->SetBranchAddress("fTreeCascVarPosPz"           , &lPosPz);
    lTreeMC->SetBranchAddress("fTreeCascVarNegPx"           , &lNegPx);
    lTreeMC->SetBranchAddress("fTreeCascVarNegPy"           , &lNegPy);
    lTreeMC->SetBranchAddress("fTreeCascVarNegPz"           , &lNegPz);
    lTreeMC->SetBranchAddress("fTreeCascVarBachPx"          , &lBachPx);
    lTreeMC->SetBranchAddress("fTreeCascVarBachPy"          , &lBachPy);
    lTreeMC->SetBranchAddress("fTreeCascVarBachPz"          , &lBachPz);
    lTreeMC->SetBranchAddress("fTreeCascVarPosInnerP"       , &lPosInnerP);
    lTreeMC->SetBranchAddress("fTreeCascVarNegInnerP"       , &lNegInnerP);
    lTreeMC->SetBranchAddress("fTreeCascVarBachInnerP"      , &lBachInnerP);
    //--- TPC Variables -----------------------------------------------
    lTreeMC->SetBranchAddress("fTreeCascVarPosNSigmaProton", &lNSigmasPosProton);
    lTreeMC->SetBranchAddress("fTreeCascVarNegNSigmaProton", &lNSigmasNegProton);
    lTreeMC->SetBranchAddress("fTreeCascVarPosNSigmaPion"  , &lNSigmasPosPion);
    lTreeMC->SetBranchAddress("fTreeCascVarNegNSigmaPion"  , &lNSigmasNegPion);
    lTreeMC->SetBranchAddress("fTreeCascVarBachNSigmaPion" , &lNSigmasBachPion);
    lTreeMC->SetBranchAddress("fTreeCascVarBachNSigmaKaon" , &lNSigmasBachKaon);
    lTreeMC->SetBranchAddress("fTreeCascVarPosdEdx"        , &lPosdEdx);
    lTreeMC->SetBranchAddress("fTreeCascVarNegdEdx"        , &lNegdEdx);
    lTreeMC->SetBranchAddress("fTreeCascVarBachdEdx"       , &lBachdEdx);
    lTreeMC->SetBranchAddress("fTreeCascVarMinTrackLength" , &lMinTrackLength);
    lTreeMC->SetBranchAddress("fTreeCascVarMaxChi2PerCluster", &lMaxChi2PerCluster);
    //--- Topological selection variables -----------------------------
    lTreeMC->SetBranchAddress("fTreeCascVarV0Radius"                 , &lV0Radius);
    lTreeMC->SetBranchAddress("fTreeCascVarCascRadius"               , &lCascRadius);
    lTreeMC->SetBranchAddress("fTreeCascVarV0Mass"                   , &lV0Mass);
    lTreeMC->SetBranchAddress("fTreeCascVarV0CosPointingAngle"       , &lV0CosinePointingAngle);
    lTreeMC->SetBranchAddress("fTreeCascVarV0CosPointingAngleSpecial", &lV0CosinePointingAngleSpecial);
    lTreeMC->SetBranchAddress("fTreeCascVarCascCosPointingAngle"     , &lCascCosinePointingAngle);
    lTreeMC->SetBranchAddress("fTreeCascVarDCANegToPrimVtx"          , &lDcaNegToPrimVertex);
    lTreeMC->SetBranchAddress("fTreeCascVarDCAPosToPrimVtx"          , &lDcaPosToPrimVertex);
    lTreeMC->SetBranchAddress("fTreeCascVarDCABachToPrimVtx"         , &lDcaBachToPrimVertex);
    lTreeMC->SetBranchAddress("fTreeCascVarDCAV0ToPrimVtx"           , &lDcaV0ToPV);
    lTreeMC->SetBranchAddress("fTreeCascVarCascDCAtoPVxy"            , &lDCAxyCascToPV);
    lTreeMC->SetBranchAddress("fTreeCascVarCascDCAtoPVz"             , &lDCAzCascToPV);
    lTreeMC->SetBranchAddress("fTreeCascVarDCAV0Daughters"           , &lDcaV0Daughters);
    lTreeMC->SetBranchAddress("fTreeCascVarDCACascDaughters"         , &lDcaCascDaughters);
    lTreeMC->SetBranchAddress("fTreeCascVarDistOverTotMom"           , &lDistOverTotMom);
    lTreeMC->SetBranchAddress("fTreeCascVarNegDCAz"                  , &lDCAzNegToPrimVertex);
    lTreeMC->SetBranchAddress("fTreeCascVarPosDCAz"                  , &lDCAzPosToPrimVertex);
    lTreeMC->SetBranchAddress("fTreeCascVarBachDCAz"                 , &lDCAzBachToPrimVertex);
    lTreeMC->SetBranchAddress("fTreeCascVarDCABachToBaryon"          , &lDCABachToBaryon);
    lTreeMC->SetBranchAddress("fTreeCascVarWrongCosPA"               , &lBBCosPA);
    //--- MC Truth Info ------------------------------------------------
    lTreeMC->SetBranchAddress("fTreeCascVarPID"              , &lPID);
    lTreeMC->SetBranchAddress("fTreeCascVarPIDBachelor"      , &lPIDBachelor);
    lTreeMC->SetBranchAddress("fTreeCascVarPIDNegative"      , &lPIDNegative);
    lTreeMC->SetBranchAddress("fTreeCascVarPIDPositive"      , &lPIDPositive);
    lTreeMC->SetBranchAddress("fTreeCascVarIsPhysicalPrimary", &lPrimaryStatus);
    //--- ITS flag -----------------------------------------------------
    lTreeMC->SetBranchAddress("fTreeCascVarPosTrackStatus" , &lPosTrackStatus);
    lTreeMC->SetBranchAddress("fTreeCascVarNegTrackStatus" , &lNegTrackStatus);
    lTreeMC->SetBranchAddress("fTreeCascVarBachTrackStatus", &lBachTrackStatus);
    //--- TOF info -----------------------------------------------------
    lTreeMC->SetBranchAddress("fTreeCascVarNegTOFExpTDiff" , &lNegTOFExpTDiff);
    lTreeMC->SetBranchAddress("fTreeCascVarPosTOFExpTDiff" , &lPosTOFExpTDiff);
    lTreeMC->SetBranchAddress("fTreeCascVarBachTOFExpTDiff", &lBachTOFExpTDiff);
    lTreeMC->SetBranchAddress("fTreeCascVarNegTOFSignal" , &lNegTOFSignal);
    lTreeMC->SetBranchAddress("fTreeCascVarPosTOFSignal" , &lPosTOFSignal);
    lTreeMC->SetBranchAddress("fTreeCascVarBachTOFSignal", &lBachTOFSignal);
    //--- Multiplicity Variable ----------------------------------------
    lTreeMC->SetBranchAddress("fTreeCascVarCentrality", &fCentrality);
    //--- MV pileup flag -----------------------------------------------
    lTreeMC->SetBranchAddress("fTreeCascVarMVPileupFlag", &fMVPileupFlag);
    //================================================================

    Long_t lNCandidatesMC = lTreeMC->GetEntries();

    //================================================================
    cout<<endl;
    Long_t lOneTenthOfNCandidatesMC = ((double)(lNCandidatesMC) / 10. );
    cout<<"--------------- MC Data File Loop  ---------------------"<<endl;
    for(Long_t icand = 0;icand<lNCandidatesMC;icand++){
        lTreeMC->GetEntry(icand);

        // check MV Pileup rejection
        if( fMVPileupSwitch && !fMVPileupFlag ) continue;

        //Multiplicity Switch
        lMultiplicity = (Double_t)fCentrality;
        if( fPerformMultiplicityStudy && (lMultiplicity<lLoMultBoundMC || lMultiplicity>lHiMultBoundMC) ) continue;

        if( icand % lOneTenthOfNCandidatesMC == 0 )
            cout<<" Currently at candidate........: "<<icand<<" / "<<lNCandidatesMC<<" ( "<<(long)(((double)(icand)/(double)(lNCandidatesMC))*(100.+1e-3))<<"% )"<<endl;

        // NOTE: pT_reco (and rap_reco) used for efficiency computation -- see details
        //       in Fiorella's presentation here: https://indico.cern.ch/event/675315

        //CMS Shift: apply before!
        if(fRapidityType == "CMS") lRap   = lRap   + fRapidityShift;
        if(fRapidityType == "CMS") lRapMC = lRapMC + fRapidityShift;

        //Compute 3D DCA Cascade to PV
        lDCACascToPV = TMath::Sqrt( lDCAxyCascToPV*lDCAxyCascToPV + lDCAzCascToPV*lDCAzCascToPV );

        //--- TPC dE/dx QA
        Float_t lPosP  = TMath::Sqrt( lPosPx*lPosPx + lPosPy*lPosPy + lPosPz*lPosPz );
        Float_t lNegP  = TMath::Sqrt( lNegPx*lNegPx + lNegPy*lNegPy + lNegPz*lNegPz );
        Float_t lBachP = TMath::Sqrt( lBachPx*lBachPx + lBachPy*lBachPy + lBachPz*lBachPz );
        Float_t lPosNSigmas  = -999;
        Float_t lNegNSigmas  = -999;
        Float_t lBachNSigmas = -999;
        if(fWhichParticle=="XiMinus")    { lPosNSigmas = lNSigmasPosProton; lNegNSigmas = lNSigmasNegPion;   lBachNSigmas = lNSigmasBachPion;  }
        if(fWhichParticle=="XiPlus")     { lPosNSigmas = lNSigmasPosPion;   lNegNSigmas = lNSigmasNegProton; lBachNSigmas = lNSigmasBachPion;  }
        if(fWhichParticle=="OmegaMinus") { lPosNSigmas = lNSigmasPosProton; lNegNSigmas = lNSigmasNegPion;   lBachNSigmas = lNSigmasBachKaon;  }
        if(fWhichParticle=="OmegaPlus")  { lPosNSigmas = lNSigmasPosPion;   lNegNSigmas = lNSigmasNegProton; lBachNSigmas = lNSigmasBachKaon;  }

        lWeAreAtBin = fHistPt->FindBin( lPt ) - 1;
        if(lWeAreAtBin == -1) lWeAreAtBin = 99; //UnderFlow, special treatment

        if(fSaveVarHistosSwitch) {
            //--- Fill histograms for signal+center background candidates
            if(lInvariantMass>lPeakLeftLimit[lWeAreAtBin] && lInvariantMass<lPeakRightLimit[lWeAreAtBin]) {
                lHistSigPlusCenterBgMCAll[ V0RADIUS          ] -> Fill( lPt, lV0Radius );
                lHistSigPlusCenterBgMCAll[ CASCRADIUS        ] -> Fill( lPt, lCascRadius );
                lHistSigPlusCenterBgMCAll[ DCANEGTOPV        ] -> Fill( lPt, lDcaNegToPrimVertex );
                lHistSigPlusCenterBgMCAll[ DCAPOSTOPV        ] -> Fill( lPt, lDcaPosToPrimVertex );
                lHistSigPlusCenterBgMCAll[ DCABACHTOPV       ] -> Fill( lPt, lDcaBachToPrimVertex );
                lHistSigPlusCenterBgMCAll[ DCAV0TOPV         ] -> Fill( lPt, lDcaV0ToPV );
                lHistSigPlusCenterBgMCAll[ DCACASCTOPV       ] -> Fill( lPt, lDCACascToPV );
                lHistSigPlusCenterBgMCAll[ DCAXYCASCTOPV     ] -> Fill( lPt, lDCAxyCascToPV );
                lHistSigPlusCenterBgMCAll[ DCAZCASCTOPV      ] -> Fill( lPt, lDCAzCascToPV );
                lHistSigPlusCenterBgMCAll[ DCAV0DAUGHTERS    ] -> Fill( lPt, lDcaV0Daughters );
                lHistSigPlusCenterBgMCAll[ DCACASCDAUGHTERS  ] -> Fill( lPt, lDcaCascDaughters );
                lHistSigPlusCenterBgMCAll[ DCAZNEGTOPV       ] -> Fill( lPt, lDCAzNegToPrimVertex );
                lHistSigPlusCenterBgMCAll[ DCAZPOSTOPV       ] -> Fill( lPt, lDCAzPosToPrimVertex );
                lHistSigPlusCenterBgMCAll[ DCAZBACHTOPV      ] -> Fill( lPt, lDCAzBachToPrimVertex );
                lHistSigPlusCenterBgMCAll[ DCABACHTOBARYON   ] -> Fill( lPt, lDCABachToBaryon );
                lHistSigPlusCenterBgMCAll[ V0PA              ] -> Fill( lPt, TMath::ACos( lV0CosinePointingAngle ) );
                lHistSigPlusCenterBgMCAll[ CASCPA            ] -> Fill( lPt, TMath::ACos( lCascCosinePointingAngle ) );
                lHistSigPlusCenterBgMCAll[ BBPA              ] -> Fill( lPt, TMath::ACos( lBBCosPA ) );
                lHistSigPlusCenterBgMCAll[ V0MASS            ] -> Fill( lPt, lV0Mass );
                lHistSigPlusCenterBgMCAll[ PROPERLIFETIME    ] -> Fill( lPt, lDistOverTotMom*lParticleMass );
                lHistSigPlusCenterBgMCAll[ COMPETINGSPECIES  ] -> Fill( lPt, lCompetingParticleMass );
                lHistSigPlusCenterBgMCAll[ TPCNEGDEDX        ] -> Fill( lNegP,  lNegdEdx );
                lHistSigPlusCenterBgMCAll[ TPCPOSDEDX        ] -> Fill( lPosP,  lPosdEdx );
                lHistSigPlusCenterBgMCAll[ TPCBACHDEDX       ] -> Fill( lBachP, lBachdEdx );
                lHistSigPlusCenterBgMCAll[ TPCNEGPIDNSIGMAS  ] -> Fill( lPt, lNegNSigmas );
                lHistSigPlusCenterBgMCAll[ TPCPOSPIDNSIGMAS  ] -> Fill( lPt, lPosNSigmas );
                lHistSigPlusCenterBgMCAll[ TPCBACHPIDNSIGMAS ] -> Fill( lPt, lBachNSigmas );
                lHistSigPlusCenterBgMCAll[ TPCNCLUSTERS      ] -> Fill( lPt, lLeastNbrClusters );
                lHistSigPlusCenterBgMCAll[ MAXCHI2PERCLUSTER ] -> Fill( lPt, lMaxChi2PerCluster );
                lHistSigPlusCenterBgMCAll[ MINTRACKLENGTH    ] -> Fill( lPt, lMinTrackLength );
                lHistSigPlusCenterBgMCAll[ NEGDAUGHTERETA    ] -> Fill( lPt, lNegEta );
                lHistSigPlusCenterBgMCAll[ POSDAUGHTERETA    ] -> Fill( lPt, lPosEta );
                lHistSigPlusCenterBgMCAll[ BACHDAUGHTERETA   ] -> Fill( lPt, lBachEta );
                lHistSigPlusCenterBgMCAll[ NEGINNERP         ] -> Fill( lPt, lNegInnerP );
                lHistSigPlusCenterBgMCAll[ POSINNERP         ] -> Fill( lPt, lPosInnerP );
                lHistSigPlusCenterBgMCAll[ BACHINNERP        ] -> Fill( lPt, lBachInnerP );
                lHistSigPlusCenterBgMCAll[ NEGTOFEXPTDIFF    ] -> Fill( lPt, lNegTOFExpTDiff );
                lHistSigPlusCenterBgMCAll[ POSTOFEXPTDIFF    ] -> Fill( lPt, lPosTOFExpTDiff );
                lHistSigPlusCenterBgMCAll[ BACHTOFEXPTDIFF   ] -> Fill( lPt, lBachTOFExpTDiff );
                lHistSigPlusCenterBgMCAll[ NEGTOFSIGNAL      ] -> Fill( lPt, lNegTOFSignal );
                lHistSigPlusCenterBgMCAll[ POSTOFSIGNAL      ] -> Fill( lPt, lPosTOFSignal );
                lHistSigPlusCenterBgMCAll[ BACHTOFSIGNAL     ] -> Fill( lPt, lBachTOFSignal );
            }

            //--- Fill histograms for left+right background candidates
            if( (lInvariantMass>lLeftBgLeftLimit[lWeAreAtBin]  && lInvariantMass<lLeftBgRightLimit[lWeAreAtBin]  ) ||
                (lInvariantMass>lRightBgLeftLimit[lWeAreAtBin] && lInvariantMass<lRightBgRightLimit[lWeAreAtBin] ) ) {
                lHistLeftPlusRightBgMCAll[ V0RADIUS          ] -> Fill( lPt, lV0Radius );
                lHistLeftPlusRightBgMCAll[ CASCRADIUS        ] -> Fill( lPt, lCascRadius );
                lHistLeftPlusRightBgMCAll[ DCANEGTOPV        ] -> Fill( lPt, lDcaNegToPrimVertex );
                lHistLeftPlusRightBgMCAll[ DCAPOSTOPV        ] -> Fill( lPt, lDcaPosToPrimVertex );
                lHistLeftPlusRightBgMCAll[ DCABACHTOPV       ] -> Fill( lPt, lDcaBachToPrimVertex );
                lHistLeftPlusRightBgMCAll[ DCAV0TOPV         ] -> Fill( lPt, lDcaV0ToPV );
                lHistLeftPlusRightBgMCAll[ DCACASCTOPV       ] -> Fill( lPt, lDCACascToPV );
                lHistLeftPlusRightBgMCAll[ DCAXYCASCTOPV     ] -> Fill( lPt, lDCAxyCascToPV );
                lHistLeftPlusRightBgMCAll[ DCAZCASCTOPV      ] -> Fill( lPt, lDCAzCascToPV );
                lHistLeftPlusRightBgMCAll[ DCAV0DAUGHTERS    ] -> Fill( lPt, lDcaV0Daughters );
                lHistLeftPlusRightBgMCAll[ DCACASCDAUGHTERS  ] -> Fill( lPt, lDcaCascDaughters );
                lHistLeftPlusRightBgMCAll[ DCAZNEGTOPV       ] -> Fill( lPt, lDCAzNegToPrimVertex );
                lHistLeftPlusRightBgMCAll[ DCAZPOSTOPV       ] -> Fill( lPt, lDCAzPosToPrimVertex );
                lHistLeftPlusRightBgMCAll[ DCAZBACHTOPV      ] -> Fill( lPt, lDCAzBachToPrimVertex );
                lHistLeftPlusRightBgMCAll[ DCABACHTOBARYON   ] -> Fill( lPt, lDCABachToBaryon );
                lHistLeftPlusRightBgMCAll[ V0PA              ] -> Fill( lPt, TMath::ACos( lV0CosinePointingAngle ) );
                lHistLeftPlusRightBgMCAll[ CASCPA            ] -> Fill( lPt, TMath::ACos( lCascCosinePointingAngle ) );
                lHistLeftPlusRightBgMCAll[ BBPA              ] -> Fill( lPt, TMath::ACos( lBBCosPA ) );
                lHistLeftPlusRightBgMCAll[ V0MASS            ] -> Fill( lPt, lV0Mass );
                lHistLeftPlusRightBgMCAll[ PROPERLIFETIME    ] -> Fill( lPt, lDistOverTotMom*lParticleMass );
                lHistLeftPlusRightBgMCAll[ COMPETINGSPECIES  ] -> Fill( lPt, lCompetingParticleMass );
                lHistLeftPlusRightBgMCAll[ TPCNEGDEDX        ] -> Fill( lNegP,  lNegdEdx );
                lHistLeftPlusRightBgMCAll[ TPCPOSDEDX        ] -> Fill( lPosP,  lPosdEdx );
                lHistLeftPlusRightBgMCAll[ TPCBACHDEDX       ] -> Fill( lBachP, lBachdEdx );
                lHistLeftPlusRightBgMCAll[ TPCNEGPIDNSIGMAS  ] -> Fill( lPt, lNegNSigmas );
                lHistLeftPlusRightBgMCAll[ TPCPOSPIDNSIGMAS  ] -> Fill( lPt, lPosNSigmas );
                lHistLeftPlusRightBgMCAll[ TPCBACHPIDNSIGMAS ] -> Fill( lPt, lBachNSigmas );
                lHistLeftPlusRightBgMCAll[ TPCNCLUSTERS      ] -> Fill( lPt, lLeastNbrClusters );
                lHistLeftPlusRightBgMCAll[ MAXCHI2PERCLUSTER ] -> Fill( lPt, lMaxChi2PerCluster );
                lHistLeftPlusRightBgMCAll[ MINTRACKLENGTH    ] -> Fill( lPt, lMinTrackLength );
                lHistLeftPlusRightBgMCAll[ NEGDAUGHTERETA    ] -> Fill( lPt, lNegEta );
                lHistLeftPlusRightBgMCAll[ POSDAUGHTERETA    ] -> Fill( lPt, lPosEta );
                lHistLeftPlusRightBgMCAll[ BACHDAUGHTERETA   ] -> Fill( lPt, lBachEta );
                lHistLeftPlusRightBgMCAll[ NEGINNERP         ] -> Fill( lPt, lNegInnerP );
                lHistLeftPlusRightBgMCAll[ POSINNERP         ] -> Fill( lPt, lPosInnerP );
                lHistLeftPlusRightBgMCAll[ BACHINNERP        ] -> Fill( lPt, lBachInnerP );
                lHistLeftPlusRightBgMCAll[ NEGTOFEXPTDIFF    ] -> Fill( lPt, lNegTOFExpTDiff );
                lHistLeftPlusRightBgMCAll[ POSTOFEXPTDIFF    ] -> Fill( lPt, lPosTOFExpTDiff );
                lHistLeftPlusRightBgMCAll[ BACHTOFEXPTDIFF   ] -> Fill( lPt, lBachTOFExpTDiff );
                lHistLeftPlusRightBgMCAll[ NEGTOFSIGNAL      ] -> Fill( lPt, lNegTOFSignal );
                lHistLeftPlusRightBgMCAll[ POSTOFSIGNAL      ] -> Fill( lPt, lPosTOFSignal );
                lHistLeftPlusRightBgMCAll[ BACHTOFSIGNAL     ] -> Fill( lPt, lBachTOFSignal );
            }

            //perfect PID association, IsPhysicalPrimary association
            if( ( fWhichParticle    == "XiMinus"
                  && lPID           ==  3312 // Candidate is a XiMinus
                  && lPIDPositive   ==  2212 // Pos Daughter is p
                  && lPIDNegative   ==  -211 // Neg Daughter is pi-
                  && lPIDBachelor   ==  -211 // Bach Daughter is pi-
                  && lPrimaryStatus ==     1 // Physical Primary Criterion
                ) ||
                ( fWhichParticle    == "XiPlus"
                  && lPID           == -3312 // Candidate is a XiMinus
                  && lPIDPositive   ==   211 // Pos Daughter is pi+
                  && lPIDNegative   == -2212 // Neg Daughter is pbar
                  && lPIDBachelor   ==   211 // Bach Daughter is pi+
                  && lPrimaryStatus ==     1 // Physical Primary Criterion
                ) ||
                ( fWhichParticle    == "OmegaMinus"
                  && lPID           ==  3334 // Candidate is a XiMinus
                  && lPIDPositive   ==  2212 // Pos Daughter is p
                  && lPIDNegative   ==  -211 // Neg Daughter is pi-
                  && lPIDBachelor   ==  -321 // Bach Daughter is K-
                  && lPrimaryStatus ==     1 // Physical Primary Criterion
                ) ||
                ( fWhichParticle    == "OmegaPlus"
                  && lPID           == -3334 // Candidate is a XiMinus
                  && lPIDPositive   ==   211 // Pos Daughter is pi+
                  && lPIDNegative   == -2212 // Neg Daughter is pbar
                  && lPIDBachelor   ==   321 // Bach Daughter is K-
                  && lPrimaryStatus ==     1 // Physical Primary Criterion
                )
              ) {

                //--- Fill histograms for signal+center background candidates
                if(lInvariantMass>lPeakLeftLimit[lWeAreAtBin] && lInvariantMass<lPeakRightLimit[lWeAreAtBin]) {
                    lHistSigPlusCenterBgMCAssAll[ V0RADIUS          ] -> Fill( lPt, lV0Radius );
                    lHistSigPlusCenterBgMCAssAll[ CASCRADIUS        ] -> Fill( lPt, lCascRadius );
                    lHistSigPlusCenterBgMCAssAll[ DCANEGTOPV        ] -> Fill( lPt, lDcaNegToPrimVertex );
                    lHistSigPlusCenterBgMCAssAll[ DCAPOSTOPV        ] -> Fill( lPt, lDcaPosToPrimVertex );
                    lHistSigPlusCenterBgMCAssAll[ DCABACHTOPV       ] -> Fill( lPt, lDcaBachToPrimVertex );
                    lHistSigPlusCenterBgMCAssAll[ DCAV0TOPV         ] -> Fill( lPt, lDcaV0ToPV );
                    lHistSigPlusCenterBgMCAssAll[ DCACASCTOPV       ] -> Fill( lPt, lDCACascToPV );
                    lHistSigPlusCenterBgMCAssAll[ DCAXYCASCTOPV     ] -> Fill( lPt, lDCAxyCascToPV );
                    lHistSigPlusCenterBgMCAssAll[ DCAZCASCTOPV      ] -> Fill( lPt, lDCAzCascToPV );
                    lHistSigPlusCenterBgMCAssAll[ DCAV0DAUGHTERS    ] -> Fill( lPt, lDcaV0Daughters );
                    lHistSigPlusCenterBgMCAssAll[ DCACASCDAUGHTERS  ] -> Fill( lPt, lDcaCascDaughters );
                    lHistSigPlusCenterBgMCAssAll[ DCAZNEGTOPV       ] -> Fill( lPt, lDCAzNegToPrimVertex );
                    lHistSigPlusCenterBgMCAssAll[ DCAZPOSTOPV       ] -> Fill( lPt, lDCAzPosToPrimVertex );
                    lHistSigPlusCenterBgMCAssAll[ DCAZBACHTOPV      ] -> Fill( lPt, lDCAzBachToPrimVertex );
                    lHistSigPlusCenterBgMCAssAll[ DCABACHTOBARYON   ] -> Fill( lPt, lDCABachToBaryon );
                    lHistSigPlusCenterBgMCAssAll[ V0PA              ] -> Fill( lPt, TMath::ACos( lV0CosinePointingAngle ) );
                    lHistSigPlusCenterBgMCAssAll[ CASCPA            ] -> Fill( lPt, TMath::ACos( lCascCosinePointingAngle ) );
                    lHistSigPlusCenterBgMCAssAll[ BBPA              ] -> Fill( lPt, TMath::ACos( lBBCosPA ) );
                    lHistSigPlusCenterBgMCAssAll[ V0MASS            ] -> Fill( lPt, lV0Mass );
                    lHistSigPlusCenterBgMCAssAll[ PROPERLIFETIME    ] -> Fill( lPt, lDistOverTotMom*lParticleMass );
                    lHistSigPlusCenterBgMCAssAll[ COMPETINGSPECIES  ] -> Fill( lPt, lCompetingParticleMass );
                    lHistSigPlusCenterBgMCAssAll[ TPCNEGDEDX        ] -> Fill( lNegP,  lNegdEdx );
                    lHistSigPlusCenterBgMCAssAll[ TPCPOSDEDX        ] -> Fill( lPosP,  lPosdEdx );
                    lHistSigPlusCenterBgMCAssAll[ TPCBACHDEDX       ] -> Fill( lBachP, lBachdEdx );
                    lHistSigPlusCenterBgMCAssAll[ TPCNEGPIDNSIGMAS  ] -> Fill( lPt, lNegNSigmas );
                    lHistSigPlusCenterBgMCAssAll[ TPCPOSPIDNSIGMAS  ] -> Fill( lPt, lPosNSigmas );
                    lHistSigPlusCenterBgMCAssAll[ TPCBACHPIDNSIGMAS ] -> Fill( lPt, lBachNSigmas );
                    lHistSigPlusCenterBgMCAssAll[ TPCNCLUSTERS      ] -> Fill( lPt, lLeastNbrClusters );
                    lHistSigPlusCenterBgMCAssAll[ MAXCHI2PERCLUSTER ] -> Fill( lPt, lMaxChi2PerCluster );
                    lHistSigPlusCenterBgMCAssAll[ MINTRACKLENGTH    ] -> Fill( lPt, lMinTrackLength );
                    lHistSigPlusCenterBgMCAssAll[ NEGDAUGHTERETA    ] -> Fill( lPt, lNegEta );
                    lHistSigPlusCenterBgMCAssAll[ POSDAUGHTERETA    ] -> Fill( lPt, lPosEta );
                    lHistSigPlusCenterBgMCAssAll[ BACHDAUGHTERETA   ] -> Fill( lPt, lBachEta );
                    lHistSigPlusCenterBgMCAssAll[ NEGINNERP         ] -> Fill( lPt, lNegInnerP );
                    lHistSigPlusCenterBgMCAssAll[ POSINNERP         ] -> Fill( lPt, lPosInnerP );
                    lHistSigPlusCenterBgMCAssAll[ BACHINNERP        ] -> Fill( lPt, lBachInnerP );
                    lHistSigPlusCenterBgMCAssAll[ NEGTOFEXPTDIFF    ] -> Fill( lPt, lNegTOFExpTDiff );
                    lHistSigPlusCenterBgMCAssAll[ POSTOFEXPTDIFF    ] -> Fill( lPt, lPosTOFExpTDiff );
                    lHistSigPlusCenterBgMCAssAll[ BACHTOFEXPTDIFF   ] -> Fill( lPt, lBachTOFExpTDiff );
                    lHistSigPlusCenterBgMCAssAll[ NEGTOFSIGNAL      ] -> Fill( lPt, lNegTOFSignal );
                    lHistSigPlusCenterBgMCAssAll[ POSTOFSIGNAL      ] -> Fill( lPt, lPosTOFSignal );
                    lHistSigPlusCenterBgMCAssAll[ BACHTOFSIGNAL     ] -> Fill( lPt, lBachTOFSignal );
                }

                //--- Fill histograms for left+right background candidates
                if( (lInvariantMass>lLeftBgLeftLimit[lWeAreAtBin]  && lInvariantMass<lLeftBgRightLimit[lWeAreAtBin]  ) ||
                    (lInvariantMass>lRightBgLeftLimit[lWeAreAtBin] && lInvariantMass<lRightBgRightLimit[lWeAreAtBin] ) ) {
                    lHistLeftPlusRightBgMCAssAll[ V0RADIUS          ] -> Fill( lPt, lV0Radius );
                    lHistLeftPlusRightBgMCAssAll[ CASCRADIUS        ] -> Fill( lPt, lCascRadius );
                    lHistLeftPlusRightBgMCAssAll[ DCANEGTOPV        ] -> Fill( lPt, lDcaNegToPrimVertex );
                    lHistLeftPlusRightBgMCAssAll[ DCAPOSTOPV        ] -> Fill( lPt, lDcaPosToPrimVertex );
                    lHistLeftPlusRightBgMCAssAll[ DCABACHTOPV       ] -> Fill( lPt, lDcaBachToPrimVertex );
                    lHistLeftPlusRightBgMCAssAll[ DCAV0TOPV         ] -> Fill( lPt, lDcaV0ToPV );
                    lHistLeftPlusRightBgMCAssAll[ DCACASCTOPV       ] -> Fill( lPt, lDCACascToPV );
                    lHistLeftPlusRightBgMCAssAll[ DCAXYCASCTOPV     ] -> Fill( lPt, lDCAxyCascToPV );
                    lHistLeftPlusRightBgMCAssAll[ DCAZCASCTOPV      ] -> Fill( lPt, lDCAzCascToPV );
                    lHistLeftPlusRightBgMCAssAll[ DCAV0DAUGHTERS    ] -> Fill( lPt, lDcaV0Daughters );
                    lHistLeftPlusRightBgMCAssAll[ DCACASCDAUGHTERS  ] -> Fill( lPt, lDcaCascDaughters );
                    lHistLeftPlusRightBgMCAssAll[ DCAZNEGTOPV       ] -> Fill( lPt, lDCAzNegToPrimVertex );
                    lHistLeftPlusRightBgMCAssAll[ DCAZPOSTOPV       ] -> Fill( lPt, lDCAzPosToPrimVertex );
                    lHistLeftPlusRightBgMCAssAll[ DCAZBACHTOPV      ] -> Fill( lPt, lDCAzBachToPrimVertex );
                    lHistLeftPlusRightBgMCAssAll[ DCABACHTOBARYON   ] -> Fill( lPt, lDCABachToBaryon );
                    lHistLeftPlusRightBgMCAssAll[ V0PA              ] -> Fill( lPt, TMath::ACos( lV0CosinePointingAngle ) );
                    lHistLeftPlusRightBgMCAssAll[ CASCPA            ] -> Fill( lPt, TMath::ACos( lCascCosinePointingAngle ) );
                    lHistLeftPlusRightBgMCAssAll[ BBPA              ] -> Fill( lPt, TMath::ACos( lBBCosPA ) );
                    lHistLeftPlusRightBgMCAssAll[ V0MASS            ] -> Fill( lPt, lV0Mass );
                    lHistLeftPlusRightBgMCAssAll[ PROPERLIFETIME    ] -> Fill( lPt, lDistOverTotMom*lParticleMass );
                    lHistLeftPlusRightBgMCAssAll[ COMPETINGSPECIES  ] -> Fill( lPt, lCompetingParticleMass );
                    lHistLeftPlusRightBgMCAssAll[ TPCNEGDEDX        ] -> Fill( lNegP,  lNegdEdx );
                    lHistLeftPlusRightBgMCAssAll[ TPCPOSDEDX        ] -> Fill( lPosP,  lPosdEdx );
                    lHistLeftPlusRightBgMCAssAll[ TPCBACHDEDX       ] -> Fill( lBachP, lBachdEdx );
                    lHistLeftPlusRightBgMCAssAll[ TPCNEGPIDNSIGMAS  ] -> Fill( lPt, lNegNSigmas );
                    lHistLeftPlusRightBgMCAssAll[ TPCPOSPIDNSIGMAS  ] -> Fill( lPt, lPosNSigmas );
                    lHistLeftPlusRightBgMCAssAll[ TPCBACHPIDNSIGMAS ] -> Fill( lPt, lBachNSigmas );
                    lHistLeftPlusRightBgMCAssAll[ TPCNCLUSTERS      ] -> Fill( lPt, lLeastNbrClusters );
                    lHistLeftPlusRightBgMCAssAll[ MAXCHI2PERCLUSTER ] -> Fill( lPt, lMaxChi2PerCluster );
                    lHistLeftPlusRightBgMCAssAll[ MINTRACKLENGTH    ] -> Fill( lPt, lMinTrackLength );
                    lHistLeftPlusRightBgMCAssAll[ NEGDAUGHTERETA    ] -> Fill( lPt, lNegEta );
                    lHistLeftPlusRightBgMCAssAll[ POSDAUGHTERETA    ] -> Fill( lPt, lPosEta );
                    lHistLeftPlusRightBgMCAssAll[ BACHDAUGHTERETA   ] -> Fill( lPt, lBachEta );
                    lHistLeftPlusRightBgMCAssAll[ NEGINNERP         ] -> Fill( lPt, lNegInnerP );
                    lHistLeftPlusRightBgMCAssAll[ POSINNERP         ] -> Fill( lPt, lPosInnerP );
                    lHistLeftPlusRightBgMCAssAll[ BACHINNERP        ] -> Fill( lPt, lBachInnerP );
                    lHistLeftPlusRightBgMCAssAll[ NEGTOFEXPTDIFF    ] -> Fill( lPt, lNegTOFExpTDiff );
                    lHistLeftPlusRightBgMCAssAll[ POSTOFEXPTDIFF    ] -> Fill( lPt, lPosTOFExpTDiff );
                    lHistLeftPlusRightBgMCAssAll[ BACHTOFEXPTDIFF   ] -> Fill( lPt, lBachTOFExpTDiff );
                    lHistLeftPlusRightBgMCAssAll[ NEGTOFSIGNAL      ] -> Fill( lPt, lNegTOFSignal );
                    lHistLeftPlusRightBgMCAssAll[ POSTOFSIGNAL      ] -> Fill( lPt, lPosTOFSignal );
                    lHistLeftPlusRightBgMCAssAll[ BACHTOFSIGNAL     ] -> Fill( lPt, lBachTOFSignal );
                }

            }

        }

        //Now check validity
        if( lRap<fRapidityBoundaryUpper && lRap>fRapidityBoundaryLower &&
            (//charge condition (x-check)
              (fWhichParticle == "XiMinus"    && lCharge == -1) ||
              (fWhichParticle == "XiPlus"     && lCharge ==  1) ||
              (fWhichParticle == "OmegaMinus" && lCharge == -1) ||
              (fWhichParticle == "OmegaPlus"  && lCharge ==  1)
            ) && //Daughter Eta
            TMath::Abs(lNegEta)       <= GetCutDaughterEta      ( lPt ) &&
            TMath::Abs(lPosEta)       <= GetCutDaughterEta      ( lPt ) &&
            TMath::Abs(lBachEta)      <= GetCutDaughterEta      ( lPt ) &&
            //Topological Selections
            lV0Radius                 >= GetCutV0Radius         ( lPt ) &&
            lCascRadius               >= GetCutCascRadius       ( lPt ) &&
            TMath::Abs(lV0Mass-1.116) <= GetCutV0Mass           ( lPt ) &&
            lV0CosinePointingAngle    >= GetCutV0CosPA          ( lPt ) &&
            lCascCosinePointingAngle  >= GetCutCascCosPA        ( lPt ) &&
            lDcaNegToPrimVertex       >= GetCutDCANegToPV       ( lPt ) &&
            lDcaPosToPrimVertex       >= GetCutDCAPosToPV       ( lPt ) &&
            lDcaBachToPrimVertex      >= GetCutDCABachToPV      ( lPt ) &&
            lDcaV0Daughters           <= GetCutDCAV0Daughters   ( lPt ) &&
            lDcaCascDaughters         <= GetCutDCACascDaughters ( lPt ) &&
            lDcaV0ToPV                >= GetCutDCAV0ToPV        ( lPt ) &&
            lDCACascToPV              <= GetCutDCACascToPV      ( lPt ) &&
            lDCAxyCascToPV            <= GetCutDCAxyCascToPV    ( lPt ) &&
            lDCAzCascToPV             <= GetCutDCAzCascToPV     ( lPt ) &&
            lParticleMass*lDistOverTotMom     <= GetCutProperLifetime   ( lPt ) &&
            lDCAzNegToPrimVertex      <= GetCutDCAzNegToPV      ( lPt ) &&
            lDCAzPosToPrimVertex      <= GetCutDCAzPosToPV      ( lPt ) &&
            lDCAzBachToPrimVertex     <= GetCutDCAzBachToPV     ( lPt ) &&
            lDCABachToBaryon          >= GetCutDCABachToBaryon  ( lPt ) &&
            lBBCosPA                  <= GetCutBBCosPA          ( lPt ) &&
            lMinTrackLength           >= GetCutMinTrackLength   ( lPt ) &&

            //Competing Species Rejection (only for Omegas)
            TMath::Abs( lCompetingParticleMass - 1.32171 ) >= GetCutCompetingSpecies ( lPt ) &&

            //Causality Cut
            TMath::Abs( lV0Radius - lCascRadius ) >= GetCutCausality( lPt ) &&

            //Nclusters cut
            lLeastNbrClusters >= GetCutLeastNumberOfClusters( lPt ) &&

            //( //official response code -- NOT APPLIED TO MC
            //  ( fWhichParticle == "XiMinus"
            //    && TMath::Abs(lNSigmasNegPion)   <= GetCutTPCPIDNSigmas( lPt )
            //    && TMath::Abs(lNSigmasPosProton) <= GetCutTPCPIDNSigmas( lPt )
            //    && TMath::Abs(lNSigmasBachPion)  <= GetCutTPCPIDNSigmas( lPt )) ||
            //  ( fWhichParticle == "XiPlus"
            //    && TMath::Abs(lNSigmasPosPion)   <= GetCutTPCPIDNSigmas( lPt )
            //    && TMath::Abs(lNSigmasNegProton) <= GetCutTPCPIDNSigmas( lPt )
            //    && TMath::Abs(lNSigmasBachPion)  <= GetCutTPCPIDNSigmas( lPt )) ||
            //  ( fWhichParticle == "OmegaMinus"
            //    && TMath::Abs(lNSigmasNegPion)   <= GetCutTPCPIDNSigmas( lPt )
            //    && TMath::Abs(lNSigmasPosProton) <= GetCutTPCPIDNSigmas( lPt )
            //    && TMath::Abs(lNSigmasBachKaon)  <= GetCutTPCPIDNSigmas( lPt )) ||
            //  ( fWhichParticle == "OmegaPlus"
            //    && TMath::Abs(lNSigmasPosPion)   <= GetCutTPCPIDNSigmas( lPt )
            //    && TMath::Abs(lNSigmasNegProton) <= GetCutTPCPIDNSigmas( lPt )
            //    && TMath::Abs(lNSigmasBachKaon)  <= GetCutTPCPIDNSigmas( lPt ))
            //) &&


            ( //TOF match 
              CheckTOFmatchBach( lPt, lBachTOFExpTDiff ) &&
              CheckTOFmatchNeg(  lPt, lNegTOFExpTDiff  ) &&
              CheckTOFmatchPos(  lPt, lPosTOFExpTDiff  ) &&
              CheckTOFmatchOne(  lPt, lBachTOFExpTDiff, lNegTOFExpTDiff, lPosTOFExpTDiff )
            ) &&

            ( //ITS refit 
              CheckITSrefitBach( lPt, lBachTrackStatus ) &&
              CheckITSrefitNeg(  lPt, lNegTrackStatus  ) &&
              CheckITSrefitPos(  lPt, lPosTrackStatus  ) &&
              CheckITSrefitOne(  lPt, lBachTrackStatus, lNegTrackStatus, lPosTrackStatus )
            ) &&
                
            ( //ITSTOF 
              CheckITSTOFBach( lPt, lBachTrackStatus, lBachTOFExpTDiff ) &&
              CheckITSTOFNeg(  lPt, lNegTrackStatus,  lNegTOFExpTDiff  ) &&
              CheckITSTOFPos(  lPt, lPosTrackStatus,  lPosTOFExpTDiff  ) &&
              CheckITSTOFOne(  lPt, lBachTrackStatus, lNegTrackStatus, lPosTrackStatus, lBachTOFExpTDiff, lNegTOFExpTDiff, lPosTOFExpTDiff )
            )

          ) { // Start Entry Loop

            if(fSaveVarHistosSwitch) {
                //--- Fill histograms for signal+center background candidates
                if( lInvariantMass>lPeakLeftLimit[lWeAreAtBin] && lInvariantMass<lPeakRightLimit[lWeAreAtBin] ) {
                    lHistSigPlusCenterBgMCSel[ V0RADIUS          ] -> Fill( lPt, lV0Radius );
                    lHistSigPlusCenterBgMCSel[ CASCRADIUS        ] -> Fill( lPt, lCascRadius );
                    lHistSigPlusCenterBgMCSel[ DCANEGTOPV        ] -> Fill( lPt, lDcaNegToPrimVertex );
                    lHistSigPlusCenterBgMCSel[ DCAPOSTOPV        ] -> Fill( lPt, lDcaPosToPrimVertex );
                    lHistSigPlusCenterBgMCSel[ DCABACHTOPV       ] -> Fill( lPt, lDcaBachToPrimVertex );
                    lHistSigPlusCenterBgMCSel[ DCAV0TOPV         ] -> Fill( lPt, lDcaV0ToPV );
                    lHistSigPlusCenterBgMCSel[ DCACASCTOPV       ] -> Fill( lPt, lDCACascToPV );
                    lHistSigPlusCenterBgMCSel[ DCAXYCASCTOPV     ] -> Fill( lPt, lDCAxyCascToPV );
                    lHistSigPlusCenterBgMCSel[ DCAZCASCTOPV      ] -> Fill( lPt, lDCAzCascToPV );
                    lHistSigPlusCenterBgMCSel[ DCAV0DAUGHTERS    ] -> Fill( lPt, lDcaV0Daughters );
                    lHistSigPlusCenterBgMCSel[ DCACASCDAUGHTERS  ] -> Fill( lPt, lDcaCascDaughters );
                    lHistSigPlusCenterBgMCSel[ DCAZNEGTOPV       ] -> Fill( lPt, lDCAzNegToPrimVertex );
                    lHistSigPlusCenterBgMCSel[ DCAZPOSTOPV       ] -> Fill( lPt, lDCAzPosToPrimVertex );
                    lHistSigPlusCenterBgMCSel[ DCAZBACHTOPV      ] -> Fill( lPt, lDCAzBachToPrimVertex );
                    lHistSigPlusCenterBgMCSel[ DCABACHTOBARYON   ] -> Fill( lPt, lDCABachToBaryon );
                    lHistSigPlusCenterBgMCSel[ V0PA              ] -> Fill( lPt, TMath::ACos( lV0CosinePointingAngle ) );
                    lHistSigPlusCenterBgMCSel[ CASCPA            ] -> Fill( lPt, TMath::ACos( lCascCosinePointingAngle ) );
                    lHistSigPlusCenterBgMCSel[ BBPA              ] -> Fill( lPt, TMath::ACos( lBBCosPA ) );
                    lHistSigPlusCenterBgMCSel[ V0MASS            ] -> Fill( lPt, lV0Mass );
                    lHistSigPlusCenterBgMCSel[ PROPERLIFETIME    ] -> Fill( lPt, lDistOverTotMom*lParticleMass );
                    lHistSigPlusCenterBgMCSel[ COMPETINGSPECIES  ] -> Fill( lPt, lCompetingParticleMass );
                    lHistSigPlusCenterBgMCSel[ TPCNEGDEDX        ] -> Fill( lNegP,  lNegdEdx );
                    lHistSigPlusCenterBgMCSel[ TPCPOSDEDX        ] -> Fill( lPosP,  lPosdEdx );
                    lHistSigPlusCenterBgMCSel[ TPCBACHDEDX       ] -> Fill( lBachP, lBachdEdx );
                    lHistSigPlusCenterBgMCSel[ TPCNEGPIDNSIGMAS  ] -> Fill( lPt, lNegNSigmas );
                    lHistSigPlusCenterBgMCSel[ TPCPOSPIDNSIGMAS  ] -> Fill( lPt, lPosNSigmas );
                    lHistSigPlusCenterBgMCSel[ TPCBACHPIDNSIGMAS ] -> Fill( lPt, lBachNSigmas );
                    lHistSigPlusCenterBgMCSel[ TPCNCLUSTERS      ] -> Fill( lPt, lLeastNbrClusters );
                    lHistSigPlusCenterBgMCSel[ MAXCHI2PERCLUSTER ] -> Fill( lPt, lMaxChi2PerCluster );
                    lHistSigPlusCenterBgMCSel[ MINTRACKLENGTH    ] -> Fill( lPt, lMinTrackLength );
                    lHistSigPlusCenterBgMCSel[ NEGDAUGHTERETA    ] -> Fill( lPt, lNegEta );
                    lHistSigPlusCenterBgMCSel[ POSDAUGHTERETA    ] -> Fill( lPt, lPosEta );
                    lHistSigPlusCenterBgMCSel[ BACHDAUGHTERETA   ] -> Fill( lPt, lBachEta );
                    lHistSigPlusCenterBgMCSel[ NEGINNERP         ] -> Fill( lPt, lNegInnerP );
                    lHistSigPlusCenterBgMCSel[ POSINNERP         ] -> Fill( lPt, lPosInnerP );
                    lHistSigPlusCenterBgMCSel[ BACHINNERP        ] -> Fill( lPt, lBachInnerP );
                    lHistSigPlusCenterBgMCSel[ NEGTOFEXPTDIFF    ] -> Fill( lPt, lNegTOFExpTDiff );
                    lHistSigPlusCenterBgMCSel[ POSTOFEXPTDIFF    ] -> Fill( lPt, lPosTOFExpTDiff );
                    lHistSigPlusCenterBgMCSel[ BACHTOFEXPTDIFF   ] -> Fill( lPt, lBachTOFExpTDiff );
                    lHistSigPlusCenterBgMCSel[ NEGTOFSIGNAL      ] -> Fill( lPt, lNegTOFSignal );
                    lHistSigPlusCenterBgMCSel[ POSTOFSIGNAL      ] -> Fill( lPt, lPosTOFSignal );
                    lHistSigPlusCenterBgMCSel[ BACHTOFSIGNAL     ] -> Fill( lPt, lBachTOFSignal );
                }

                //--- Fill histograms for left+right background candidates
                if( (lInvariantMass>lLeftBgLeftLimit[lWeAreAtBin]  && lInvariantMass<lLeftBgRightLimit[lWeAreAtBin]  ) ||
                    (lInvariantMass>lRightBgLeftLimit[lWeAreAtBin] && lInvariantMass<lRightBgRightLimit[lWeAreAtBin] ) ) {
                    lHistLeftPlusRightBgMCSel[ V0RADIUS          ] -> Fill( lPt, lV0Radius );
                    lHistLeftPlusRightBgMCSel[ CASCRADIUS        ] -> Fill( lPt, lCascRadius );
                    lHistLeftPlusRightBgMCSel[ DCANEGTOPV        ] -> Fill( lPt, lDcaNegToPrimVertex );
                    lHistLeftPlusRightBgMCSel[ DCAPOSTOPV        ] -> Fill( lPt, lDcaPosToPrimVertex );
                    lHistLeftPlusRightBgMCSel[ DCABACHTOPV       ] -> Fill( lPt, lDcaBachToPrimVertex );
                    lHistLeftPlusRightBgMCSel[ DCAV0TOPV         ] -> Fill( lPt, lDcaV0ToPV );
                    lHistLeftPlusRightBgMCSel[ DCACASCTOPV       ] -> Fill( lPt, lDCACascToPV );
                    lHistLeftPlusRightBgMCSel[ DCAXYCASCTOPV     ] -> Fill( lPt, lDCAxyCascToPV );
                    lHistLeftPlusRightBgMCSel[ DCAZCASCTOPV      ] -> Fill( lPt, lDCAzCascToPV );
                    lHistLeftPlusRightBgMCSel[ DCAV0DAUGHTERS    ] -> Fill( lPt, lDcaV0Daughters );
                    lHistLeftPlusRightBgMCSel[ DCACASCDAUGHTERS  ] -> Fill( lPt, lDcaCascDaughters );
                    lHistLeftPlusRightBgMCSel[ DCAZNEGTOPV       ] -> Fill( lPt, lDCAzNegToPrimVertex );
                    lHistLeftPlusRightBgMCSel[ DCAZPOSTOPV       ] -> Fill( lPt, lDCAzPosToPrimVertex );
                    lHistLeftPlusRightBgMCSel[ DCAZBACHTOPV      ] -> Fill( lPt, lDCAzBachToPrimVertex );
                    lHistLeftPlusRightBgMCSel[ DCABACHTOBARYON   ] -> Fill( lPt, lDCABachToBaryon );
                    lHistLeftPlusRightBgMCSel[ V0PA              ] -> Fill( lPt, TMath::ACos( lV0CosinePointingAngle ) );
                    lHistLeftPlusRightBgMCSel[ CASCPA            ] -> Fill( lPt, TMath::ACos( lCascCosinePointingAngle ) );
                    lHistLeftPlusRightBgMCSel[ BBPA              ] -> Fill( lPt, TMath::ACos( lBBCosPA ) );
                    lHistLeftPlusRightBgMCSel[ V0MASS            ] -> Fill( lPt, lV0Mass );
                    lHistLeftPlusRightBgMCSel[ PROPERLIFETIME    ] -> Fill( lPt, lDistOverTotMom*lParticleMass );
                    lHistLeftPlusRightBgMCSel[ COMPETINGSPECIES  ] -> Fill( lPt, lCompetingParticleMass );
                    lHistLeftPlusRightBgMCSel[ TPCNEGDEDX        ] -> Fill( lNegP,  lNegdEdx );
                    lHistLeftPlusRightBgMCSel[ TPCPOSDEDX        ] -> Fill( lPosP,  lPosdEdx );
                    lHistLeftPlusRightBgMCSel[ TPCBACHDEDX       ] -> Fill( lBachP, lBachdEdx );
                    lHistLeftPlusRightBgMCSel[ TPCNEGPIDNSIGMAS  ] -> Fill( lPt, lNegNSigmas );
                    lHistLeftPlusRightBgMCSel[ TPCPOSPIDNSIGMAS  ] -> Fill( lPt, lPosNSigmas );
                    lHistLeftPlusRightBgMCSel[ TPCBACHPIDNSIGMAS ] -> Fill( lPt, lBachNSigmas );
                    lHistLeftPlusRightBgMCSel[ TPCNCLUSTERS      ] -> Fill( lPt, lLeastNbrClusters );
                    lHistLeftPlusRightBgMCSel[ MAXCHI2PERCLUSTER ] -> Fill( lPt, lMaxChi2PerCluster );
                    lHistLeftPlusRightBgMCSel[ MINTRACKLENGTH    ] -> Fill( lPt, lMinTrackLength );
                    lHistLeftPlusRightBgMCSel[ NEGDAUGHTERETA    ] -> Fill( lPt, lNegEta );
                    lHistLeftPlusRightBgMCSel[ POSDAUGHTERETA    ] -> Fill( lPt, lPosEta );
                    lHistLeftPlusRightBgMCSel[ BACHDAUGHTERETA   ] -> Fill( lPt, lBachEta );
                    lHistLeftPlusRightBgMCSel[ NEGINNERP         ] -> Fill( lPt, lNegInnerP );
                    lHistLeftPlusRightBgMCSel[ POSINNERP         ] -> Fill( lPt, lPosInnerP );
                    lHistLeftPlusRightBgMCSel[ BACHINNERP        ] -> Fill( lPt, lBachInnerP );
                    lHistLeftPlusRightBgMCSel[ NEGTOFEXPTDIFF    ] -> Fill( lPt, lNegTOFExpTDiff );
                    lHistLeftPlusRightBgMCSel[ POSTOFEXPTDIFF    ] -> Fill( lPt, lPosTOFExpTDiff );
                    lHistLeftPlusRightBgMCSel[ BACHTOFEXPTDIFF   ] -> Fill( lPt, lBachTOFExpTDiff );
                    lHistLeftPlusRightBgMCSel[ NEGTOFSIGNAL      ] -> Fill( lPt, lNegTOFSignal );
                    lHistLeftPlusRightBgMCSel[ POSTOFSIGNAL      ] -> Fill( lPt, lPosTOFSignal );
                    lHistLeftPlusRightBgMCSel[ BACHTOFSIGNAL     ] -> Fill( lPt, lBachTOFSignal );
                }

            }

            //perfect PID association, IsPhysicalPrimary association
            if( ( fWhichParticle    == "XiMinus"
                  && lPID           ==  3312 // Candidate is a XiMinus
                  && lPIDPositive   ==  2212 // Pos Daughter is p
                  && lPIDNegative   ==  -211 // Neg Daughter is pi-
                  && lPIDBachelor   ==  -211 // Bach Daughter is pi-
                  && lPrimaryStatus ==     1 // Physical Primary Criterion
                ) ||
                ( fWhichParticle    == "XiPlus"
                  && lPID           == -3312 // Candidate is a XiMinus
                  && lPIDPositive   ==   211 // Pos Daughter is p
                  && lPIDNegative   == -2212 // Neg Daughter is pbar
                  && lPIDBachelor   ==   211 // Bach Daughter is pi+
                  && lPrimaryStatus ==     1 // Physical Primary Criterion
                ) ||
                ( fWhichParticle    == "OmegaMinus"
                  && lPID           ==  3334 // Candidate is a XiMinus
                  && lPIDPositive   ==  2212 // Pos Daughter is p
                  && lPIDNegative   ==  -211 // Neg Daughter is pi-
                  && lPIDBachelor   ==  -321 // Bach Daughter is K-
                  && lPrimaryStatus ==     1 // Physical Primary Criterion
                ) ||
                ( fWhichParticle    == "OmegaPlus"
                  && lPID           == -3334 // Candidate is a XiMinus
                  && lPIDPositive   ==   211 // Pos Daughter is pi+
                  && lPIDNegative   == -2212 // Neg Daughter is pbar
                  && lPIDBachelor   ==   321 // Bach Daughter is K+
                  && lPrimaryStatus ==     1 // Physical Primary Criterion
                )
              ) {
                  lHistoSelectedCascMC[lWeAreAtBin]->Fill(lInvariantMass); //fill with specific inv mass
                  f2dHistPtResolution->Fill(lPt, lPtMC);

                  //--- Peak Region
                  if(lInvariantMass>lPeakLeftLimit[lWeAreAtBin]      && lInvariantMass<lPeakRightLimit[lWeAreAtBin]     ) {
                      lSigPlusCenterBgCascMC[lWeAreAtBin]++;
                  }

                //--- Left+Right Background Region
                if( (lInvariantMass>lLeftBgLeftLimit[lWeAreAtBin]    && lInvariantMass<lLeftBgRightLimit[lWeAreAtBin]) ||
                    (lInvariantMass>lRightBgLeftLimit[lWeAreAtBin]   && lInvariantMass<lRightBgRightLimit[lWeAreAtBin]) ) {
                    lLeftPlusRightBgCascMC[lWeAreAtBin]++;
                }

                //--- Info needed for Geant-Fluka correction
                if( fWhichParticle.Contains("Minus") ) lProtonMomentum[lWeAreAtBin]->Fill( TMath::Sqrt( lPosPx*lPosPx + lPosPy*lPosPy ) );
                if( fWhichParticle.Contains("Plus")  ) lProtonMomentum[lWeAreAtBin]->Fill( TMath::Sqrt( lNegPx*lNegPx + lNegPy*lNegPy ) );

                //--- Resolution tests
                lHistResolution[lWeAreAtBin]->Fill(lPt - lPtMC);

                if(fSaveVarHistosSwitch) {
                    //--- Fill histograms for signal+center background candidates
                    if(lInvariantMass>lPeakLeftLimit[lWeAreAtBin] && lInvariantMass<lPeakRightLimit[lWeAreAtBin]) {
                        lHistSigPlusCenterBgMCAssSel[ V0RADIUS          ] -> Fill( lPt, lV0Radius );
                        lHistSigPlusCenterBgMCAssSel[ CASCRADIUS        ] -> Fill( lPt, lCascRadius );
                        lHistSigPlusCenterBgMCAssSel[ DCANEGTOPV        ] -> Fill( lPt, lDcaNegToPrimVertex );
                        lHistSigPlusCenterBgMCAssSel[ DCAPOSTOPV        ] -> Fill( lPt, lDcaPosToPrimVertex );
                        lHistSigPlusCenterBgMCAssSel[ DCABACHTOPV       ] -> Fill( lPt, lDcaBachToPrimVertex );
                        lHistSigPlusCenterBgMCAssSel[ DCAV0TOPV         ] -> Fill( lPt, lDcaV0ToPV );
                        lHistSigPlusCenterBgMCAssSel[ DCACASCTOPV       ] -> Fill( lPt, lDCACascToPV );
                        lHistSigPlusCenterBgMCAssSel[ DCAXYCASCTOPV     ] -> Fill( lPt, lDCAxyCascToPV );
                        lHistSigPlusCenterBgMCAssSel[ DCAZCASCTOPV      ] -> Fill( lPt, lDCAzCascToPV );
                        lHistSigPlusCenterBgMCAssSel[ DCAV0DAUGHTERS    ] -> Fill( lPt, lDcaV0Daughters );
                        lHistSigPlusCenterBgMCAssSel[ DCACASCDAUGHTERS  ] -> Fill( lPt, lDcaCascDaughters );
                        lHistSigPlusCenterBgMCAssSel[ DCAZNEGTOPV       ] -> Fill( lPt, lDCAzNegToPrimVertex );
                        lHistSigPlusCenterBgMCAssSel[ DCAZPOSTOPV       ] -> Fill( lPt, lDCAzPosToPrimVertex );
                        lHistSigPlusCenterBgMCAssSel[ DCAZBACHTOPV      ] -> Fill( lPt, lDCAzBachToPrimVertex );
                        lHistSigPlusCenterBgMCAssSel[ DCABACHTOBARYON   ] -> Fill( lPt, lDCABachToBaryon );
                        lHistSigPlusCenterBgMCAssSel[ V0PA              ] -> Fill( lPt, TMath::ACos( lV0CosinePointingAngle ) );
                        lHistSigPlusCenterBgMCAssSel[ CASCPA            ] -> Fill( lPt, TMath::ACos( lCascCosinePointingAngle ) );
                        lHistSigPlusCenterBgMCAssSel[ BBPA              ] -> Fill( lPt, TMath::ACos( lBBCosPA ) );
                        lHistSigPlusCenterBgMCAssSel[ V0MASS            ] -> Fill( lPt, lV0Mass );
                        lHistSigPlusCenterBgMCAssSel[ PROPERLIFETIME    ] -> Fill( lPt, lDistOverTotMom*lParticleMass );
                        lHistSigPlusCenterBgMCAssSel[ COMPETINGSPECIES  ] -> Fill( lPt, lCompetingParticleMass );
                        lHistSigPlusCenterBgMCAssSel[ TPCNEGDEDX        ] -> Fill( lNegP,  lNegdEdx );
                        lHistSigPlusCenterBgMCAssSel[ TPCPOSDEDX        ] -> Fill( lPosP,  lPosdEdx );
                        lHistSigPlusCenterBgMCAssSel[ TPCBACHDEDX       ] -> Fill( lBachP, lBachdEdx );
                        lHistSigPlusCenterBgMCAssSel[ TPCNEGPIDNSIGMAS  ] -> Fill( lPt, lNegNSigmas );
                        lHistSigPlusCenterBgMCAssSel[ TPCPOSPIDNSIGMAS  ] -> Fill( lPt, lPosNSigmas );
                        lHistSigPlusCenterBgMCAssSel[ TPCBACHPIDNSIGMAS ] -> Fill( lPt, lBachNSigmas );
                        lHistSigPlusCenterBgMCAssSel[ TPCNCLUSTERS      ] -> Fill( lPt, lLeastNbrClusters );
                        lHistSigPlusCenterBgMCAssSel[ MAXCHI2PERCLUSTER ] -> Fill( lPt, lMaxChi2PerCluster );
                        lHistSigPlusCenterBgMCAssSel[ MINTRACKLENGTH    ] -> Fill( lPt, lMinTrackLength );
                        lHistSigPlusCenterBgMCAssSel[ NEGDAUGHTERETA    ] -> Fill( lPt, lNegEta );
                        lHistSigPlusCenterBgMCAssSel[ POSDAUGHTERETA    ] -> Fill( lPt, lPosEta );
                        lHistSigPlusCenterBgMCAssSel[ BACHDAUGHTERETA   ] -> Fill( lPt, lBachEta );
                        lHistSigPlusCenterBgMCAssSel[ NEGINNERP         ] -> Fill( lPt, lNegInnerP );
                        lHistSigPlusCenterBgMCAssSel[ POSINNERP         ] -> Fill( lPt, lPosInnerP );
                        lHistSigPlusCenterBgMCAssSel[ BACHINNERP        ] -> Fill( lPt, lBachInnerP );
                        lHistSigPlusCenterBgMCAssSel[ NEGTOFEXPTDIFF    ] -> Fill( lPt, lNegTOFExpTDiff );
                        lHistSigPlusCenterBgMCAssSel[ POSTOFEXPTDIFF    ] -> Fill( lPt, lPosTOFExpTDiff );
                        lHistSigPlusCenterBgMCAssSel[ BACHTOFEXPTDIFF   ] -> Fill( lPt, lBachTOFExpTDiff );
                        lHistSigPlusCenterBgMCAssSel[ NEGTOFSIGNAL      ] -> Fill( lPt, lNegTOFSignal );
                        lHistSigPlusCenterBgMCAssSel[ POSTOFSIGNAL      ] -> Fill( lPt, lPosTOFSignal );
                        lHistSigPlusCenterBgMCAssSel[ BACHTOFSIGNAL     ] -> Fill( lPt, lBachTOFSignal );
                    }

                    //--- Fill histograms for left+right background candidates
                    if( (lInvariantMass>lLeftBgLeftLimit[lWeAreAtBin]  && lInvariantMass<lLeftBgRightLimit[lWeAreAtBin]  ) ||
                        (lInvariantMass>lRightBgLeftLimit[lWeAreAtBin] && lInvariantMass<lRightBgRightLimit[lWeAreAtBin] ) ) {
                        lHistLeftPlusRightBgMCAssSel[ V0RADIUS          ] -> Fill( lPt, lV0Radius );
                        lHistLeftPlusRightBgMCAssSel[ CASCRADIUS        ] -> Fill( lPt, lCascRadius );
                        lHistLeftPlusRightBgMCAssSel[ DCANEGTOPV        ] -> Fill( lPt, lDcaNegToPrimVertex );
                        lHistLeftPlusRightBgMCAssSel[ DCAPOSTOPV        ] -> Fill( lPt, lDcaPosToPrimVertex );
                        lHistLeftPlusRightBgMCAssSel[ DCABACHTOPV       ] -> Fill( lPt, lDcaBachToPrimVertex );
                        lHistLeftPlusRightBgMCAssSel[ DCAV0TOPV         ] -> Fill( lPt, lDcaV0ToPV );
                        lHistLeftPlusRightBgMCAssSel[ DCACASCTOPV       ] -> Fill( lPt, lDCACascToPV );
                        lHistLeftPlusRightBgMCAssSel[ DCAXYCASCTOPV     ] -> Fill( lPt, lDCAxyCascToPV );
                        lHistLeftPlusRightBgMCAssSel[ DCAZCASCTOPV      ] -> Fill( lPt, lDCAzCascToPV );
                        lHistLeftPlusRightBgMCAssSel[ DCAV0DAUGHTERS    ] -> Fill( lPt, lDcaV0Daughters );
                        lHistLeftPlusRightBgMCAssSel[ DCACASCDAUGHTERS  ] -> Fill( lPt, lDcaCascDaughters );
                        lHistLeftPlusRightBgMCAssSel[ DCAZNEGTOPV       ] -> Fill( lPt, lDCAzNegToPrimVertex );
                        lHistLeftPlusRightBgMCAssSel[ DCAZPOSTOPV       ] -> Fill( lPt, lDCAzPosToPrimVertex );
                        lHistLeftPlusRightBgMCAssSel[ DCAZBACHTOPV      ] -> Fill( lPt, lDCAzBachToPrimVertex );
                        lHistLeftPlusRightBgMCAssSel[ DCABACHTOBARYON   ] -> Fill( lPt, lDCABachToBaryon );
                        lHistLeftPlusRightBgMCAssSel[ V0PA              ] -> Fill( lPt, TMath::ACos( lV0CosinePointingAngle ) );
                        lHistLeftPlusRightBgMCAssSel[ CASCPA            ] -> Fill( lPt, TMath::ACos( lCascCosinePointingAngle ) );
                        lHistLeftPlusRightBgMCAssSel[ BBPA              ] -> Fill( lPt, TMath::ACos( lBBCosPA ) );
                        lHistLeftPlusRightBgMCAssSel[ V0MASS            ] -> Fill( lPt, lV0Mass );
                        lHistLeftPlusRightBgMCAssSel[ PROPERLIFETIME    ] -> Fill( lPt, lDistOverTotMom*lParticleMass );
                        lHistLeftPlusRightBgMCAssSel[ COMPETINGSPECIES  ] -> Fill( lPt, lCompetingParticleMass );
                        lHistLeftPlusRightBgMCAssSel[ TPCNEGDEDX        ] -> Fill( lNegP,  lNegdEdx );
                        lHistLeftPlusRightBgMCAssSel[ TPCPOSDEDX        ] -> Fill( lPosP,  lPosdEdx );
                        lHistLeftPlusRightBgMCAssSel[ TPCBACHDEDX       ] -> Fill( lBachP, lBachdEdx );
                        lHistLeftPlusRightBgMCAssSel[ TPCNEGPIDNSIGMAS  ] -> Fill( lPt, lNegNSigmas );
                        lHistLeftPlusRightBgMCAssSel[ TPCPOSPIDNSIGMAS  ] -> Fill( lPt, lPosNSigmas );
                        lHistLeftPlusRightBgMCAssSel[ TPCBACHPIDNSIGMAS ] -> Fill( lPt, lBachNSigmas );
                        lHistLeftPlusRightBgMCAssSel[ TPCNCLUSTERS      ] -> Fill( lPt, lLeastNbrClusters );
                        lHistLeftPlusRightBgMCAssSel[ MAXCHI2PERCLUSTER ] -> Fill( lPt, lMaxChi2PerCluster );
                        lHistLeftPlusRightBgMCAssSel[ MINTRACKLENGTH    ] -> Fill( lPt, lMinTrackLength );
                        lHistLeftPlusRightBgMCAssSel[ NEGDAUGHTERETA    ] -> Fill( lPt, lNegEta );
                        lHistLeftPlusRightBgMCAssSel[ POSDAUGHTERETA    ] -> Fill( lPt, lPosEta );
                        lHistLeftPlusRightBgMCAssSel[ BACHDAUGHTERETA   ] -> Fill( lPt, lBachEta );
                        lHistLeftPlusRightBgMCAssSel[ NEGINNERP         ] -> Fill( lPt, lNegInnerP );
                        lHistLeftPlusRightBgMCAssSel[ POSINNERP         ] -> Fill( lPt, lPosInnerP );
                        lHistLeftPlusRightBgMCAssSel[ BACHINNERP        ] -> Fill( lPt, lBachInnerP );
                        lHistLeftPlusRightBgMCAssSel[ NEGTOFEXPTDIFF    ] -> Fill( lPt, lNegTOFExpTDiff );
                        lHistLeftPlusRightBgMCAssSel[ POSTOFEXPTDIFF    ] -> Fill( lPt, lPosTOFExpTDiff );
                        lHistLeftPlusRightBgMCAssSel[ BACHTOFEXPTDIFF   ] -> Fill( lPt, lBachTOFExpTDiff );
                        lHistLeftPlusRightBgMCAssSel[ NEGTOFSIGNAL      ] -> Fill( lPt, lNegTOFSignal );
                        lHistLeftPlusRightBgMCAssSel[ POSTOFSIGNAL      ] -> Fill( lPt, lPosTOFSignal );
                        lHistLeftPlusRightBgMCAssSel[ BACHTOFSIGNAL     ] -> Fill( lPt, lBachTOFSignal );
                    }
                }
            }
            //____________________________________________________________________________________

            //____________________________________________________________________________________
            //This is the MC template acquisition, if requested
            if( lRap<fRapidityBoundaryUpper && lRap>fRapidityBoundaryLower &&
                ( (lMultiplicity>=fLoMultBound && lMultiplicity<=fHiMultBound) || fUseMinBiasMCBackgroundTemplate ) &&
                ( //Perfect MC NON-association (for background estimates)
                  (fWhichParticle == "XiMinus"    && lPID!= 3312 ) ||
                  (fWhichParticle == "XiPlus"     && lPID!=-3312 ) ||
                  (fWhichParticle == "OmegaMinus" && lPID!= 3334 ) ||
                  (fWhichParticle == "OmegaPlus"  && lPID!=-3334 )
                ) &&
                ( //Remove "Pi->Mu"-decayed cascades: may affect result particularly if special MC used
                  (fWhichParticle == "XiMinus"    && ( lPIDNegative != 13 && lPIDBachelor != 13 ) ) ||
                  (fWhichParticle == "XiPlus"     && ( lPIDPositive !=-13 && lPIDBachelor !=-13 ) ) ||
                  (fWhichParticle == "OmegaMinus" && ( lPIDNegative != 13 && lPIDBachelor !=-13 ) ) ||
                  (fWhichParticle == "OmegaPlus"  && ( lPIDPositive !=-13 && lPIDBachelor !=-13 ) )
                )
              ) {
                lWeAreAtBin = fHistPt->FindBin(lPt)-1;
                if(lWeAreAtBin == -1) lWeAreAtBin = 99; //UnderFlow, special treatment
                lHistoBgTemplateMC[lWeAreAtBin]->Fill(lInvariantMass); //fill with specific inv mass

                //--- Peak Region
                if(lInvariantMass>lPeakLeftLimit[lWeAreAtBin]      && lInvariantMass<lPeakRightLimit[lWeAreAtBin]     ) {
                    lCenterBgTemplateCascMC[lWeAreAtBin]++;
                }

                //--- Left+Right Background Region
                if( (lInvariantMass>lLeftBgLeftLimit[lWeAreAtBin]    && lInvariantMass<lLeftBgRightLimit[lWeAreAtBin]) ||
                    (lInvariantMass>lRightBgLeftLimit[lWeAreAtBin]   && lInvariantMass<lRightBgRightLimit[lWeAreAtBin]) ) {
                    lLeftPlusRightBgTemplateCascMC[lWeAreAtBin]++;
                }
            }

        } // End Entry Loop
    }
    cout<<"--------------- Loop Completed -------------------------"<<endl;
    cout<<endl;


    cout<<"--------------- X-Check Peak Finding (gauss+linear) ----"<<endl;
    //== Variables for holding peak position, width ==============
    Double_t lPeakPositionMC[100];
    Double_t lPeakWidthMC[100];

    char fgausnameMC[100];
    TF1 *fgausPtMC[100];
    for(Int_t ibin = 0; ibin<fptbinnumb; ibin++){
        cout<<"---> Peak Finding, bin #"<<ibin<<" (perfect count = "<<lHistoSelectedCascMC[ibin]->GetEntries()<<")..."<<endl;
        sprintf(fgausnameMC,"fGausPtMC%i",ibin);
        if ( fWhichParticle == "XiMinus" || fWhichParticle == "XiPlus" ){
            fgausPtMC[ibin]= new TF1(fgausnameMC,"[0]*TMath::Gaus(x,[1],[2])+[3]*x+[4]", 1.322-0.03, 1.322+0.03 );
            fgausPtMC[ibin]->SetParameter(1,1.322);
            fgausPtMC[ibin]->SetParameter(2,0.0025);
            fgausPtMC[ibin]->SetParLimits(2,0.001,0.01);
        }
        if ( fWhichParticle == "OmegaMinus" || fWhichParticle == "OmegaPlus"){
            fgausPtMC[ibin]= new TF1(fgausnameMC,"[0]*TMath::Gaus(x,[1],[2])+[3]*x+[4]", 1.672-0.03, 1.672+0.03 );
            fgausPtMC[ibin]->SetParameter(1,1.672);
            fgausPtMC[ibin]->SetParameter(2,0.0025);
            fgausPtMC[ibin]->SetParLimits(2,0.001,0.01);
        }
        fgausPtMC[ibin]->SetParameter(0,lHistoSelectedCascMC[ibin]->GetMaximum() * 0.9);
        //fgausPtMC[ibin]->SetParameter(3,0);
        //fgausPtMC[ibin]->SetParameter(4,HistoMBCascMC[ibin]->GetMaximum() * 0.1);
        lHistoSelectedCascMC[ibin]->Fit(fgausnameMC,"QREM0");
        lPeakPositionMC[ibin] = TMath::Abs( fgausPtMC[ibin]->GetParameter(1) );
        lPeakWidthMC[ibin] = fgausPtMC[ibin]->GetParameter(2);
        cout<<"---> ["<<fptbinlimits[ibin]<<" - "<<fptbinlimits[ibin+1]<<" GeV/c]\tPeak at: "<<lPeakPositionMC[ibin]<<", sigma = "<<lPeakWidthMC[ibin]<<endl;
        fHistPeakPositionMC->SetBinContent(ibin+1, lPeakPositionMC[ibin]);
        fHistPeakPositionMC->SetBinError(ibin+1, lPeakWidthMC[ibin]);

        //Create appropriate TLine Objects for Canvases
        lLineLeftMostMC[ibin]  = new TLine( lLeftBgLeftLimit[ibin],   0, lLeftBgLeftLimit[ibin],   lHistoSelectedCascMC[ibin]->GetMaximum() * 0.95 );
        lLineLeftMC[ibin]      = new TLine( lLeftBgRightLimit[ibin],  0, lLeftBgRightLimit[ibin],  lHistoSelectedCascMC[ibin]->GetMaximum() * 0.95 );
        lLineRightMC[ibin]     = new TLine( lRightBgLeftLimit[ibin],  0, lRightBgLeftLimit[ibin],  lHistoSelectedCascMC[ibin]->GetMaximum() * 0.95 );
        lLineRightMostMC[ibin] = new TLine( lRightBgRightLimit[ibin], 0, lRightBgRightLimit[ibin], lHistoSelectedCascMC[ibin]->GetMaximum() * 0.95 );

        //Preparing Canvas for storing...
        lCanvasHistoSelectedCascMC[ibin]->cd();
        lHistoSelectedCascMC[ibin]->Draw();
        lLineLeftMostMC[ibin]->Draw();
        lLineLeftMC[ibin]->Draw();
        lLineRightMC[ibin]->Draw();
        lLineRightMostMC[ibin]->Draw();
        fgausPtMC[ibin]->Draw("same");

    }
    cout<<"--------------- Peak Finding Finished (in MC) ----------"<<endl;
    cout<<endl;

    TF1 *lfitNoiseMC[50];
    TF1 *lSampleNoiseMC[50];
    char lFitNameOneMC[50];

    //=============================================================
    if( fFitBackgroundSwitch ){
        cout<<"--------------- Backgrounds: fitting with linear -------"<<endl;
        for(long i=0; i<fptbinnumb; i++){
            //Define Function to Fit Background
            sprintf(lFitNameOneMC,"lfitNoiseMC%i",(int)i);
            cout<<"creating fitnoise, named "<<lFitNameOneMC<<endl;
            lfitNoiseMC[i] = new TF1(lFitNameOneMC, this, &AliCascadeModule::MyBgPol1,
                    RoundToThousandth ( lLeftBgLeftLimit[i]   ),
                    RoundToThousandth ( lRightBgRightLimit[i] ), 4 , "AliCascadeModule", "MyBgPol1");
            lfitNoiseMC[i] -> FixParameter(2,RoundToThousandth ( lLeftBgRightLimit[i] ) );
            lfitNoiseMC[i] -> FixParameter(3,RoundToThousandth ( lRightBgLeftLimit[i] ) );
            lfitNoiseMC[i] -> SetParameter(0,lLeftPlusRightBgCascMC[i]*lHistoSelectedCascMC[i]->GetMaximum() / (lSigPlusCenterBgCascMC[i]+1e-6) );
            lfitNoiseMC[i] -> SetParameter(1,0);
            sprintf(lFitNameOneMC,"lSampleNoiseMC%i",(int)i);

            //Define Function to Sample Background
            cout<<"creating sample "<<i<<endl;
            lSampleNoiseMC[i] = new TF1(lFitNameOneMC, this, &AliCascadeModule::MyBgPolToEval1,
                    RoundToThousandth ( lLeftBgLeftLimit[i]   ),
                    RoundToThousandth ( lRightBgRightLimit[i] ), 2, "AliCascadeModule", "MyBgPolToEval1" );
        }
        for(long i=0; i<fptbinnumb; i++){
            cout<<"Fitting function for bin "<<i<<", get name = "<<lfitNoiseMC[i]->GetName()<<endl;
            sprintf(lFitNameOneMC,"lfitNoiseMC%i",(int)i);
            lHistoSelectedCascMC[i] -> Fit(lFitNameOneMC,"LLrie+0");
            lSampleNoiseMC[i]->SetParameter(0, lfitNoiseMC[i]->GetParameter(0) );
            lSampleNoiseMC[i]->SetParameter(1, lfitNoiseMC[i]->GetParameter(1) );
        }
        for(long i=0; i<fptbinnumb; i++){
            cout<<"Overriding Background info: Was "<<lLeftPlusRightBgCascMC[i]<<", is now "<<lSampleNoiseMC[i]->Integral( lPeakLeftLimit[i], lPeakRightLimit[i] )/lHistoSelectedCascMC[i]->GetBinWidth(5)<<endl;
            lLeftPlusRightBgCascMC[i] = lSampleNoiseMC[i]->Integral( lPeakLeftLimit[i], lPeakRightLimit[i] )/lHistoSelectedCascMC[i]->GetBinWidth(5);
        }
        cout<<"--------------- Fitting Finished! ----------------------"<<endl;
    }
    //=============================================================

    //=============================================================
    //Optional: Compute background from MC template
    // Scaling rule:
    //
    // (Real data BG estimate in peak) = (MC BG estimate in peak) * (Real data BG in side-bands) / (MC BG in side-bands)
    //
    // i.e., set:
    //
    // lLeftPlusRightBgCasc[i] = lCenterBgTemplateCascMC[i] * lLeftPlusRightBgCasc[i] / lLeftPlusRightBgTemplateCascMC[i]
    //
    // ---> Error propagation: all counters above are uncorrelated, i.e. simple statistical propagation is sufficient.

    /*
       cout<<"Using MC as background template for real data. WARNING: EXPERIMENTAL"<<endl;
       Double_t lScalingFactor         = 1;
       Double_t lScalingFactorError    = 1;
       Double_t lNewValue              = 1;
       Double_t lNewValueError         = 1;
       for(long i=0; i<fptbinnumb; i++){

       cout<<"bin #"<<i<<", lCenterBgTemplateCascMC[i] = "<<lCenterBgTemplateCascMC[i]<<", leftright = "<<lLeftPlusRightBgTemplateCascMC[i]<<endl;
       cout<<"Entries in histo: "<<lHistoBgTemplateMC[i]->GetEntries()<<endl;
       lScalingFactor      = lCenterBgTemplateCascMC[i] / lLeftPlusRightBgTemplateCascMC[i];
       lScalingFactorError = ErrorInRatio(lCenterBgTemplateCascMC[i],        TMath::Sqrt(lCenterBgTemplateCascMC[i]),
       lLeftPlusRightBgTemplateCascMC[i], TMath::Sqrt(lLeftPlusRightBgTemplateCascMC[i]) );
       lNewValue       = lLeftPlusRightBgCasc[i] * lScalingFactor;
       lNewValueError  = TMath::Sqrt( lLeftPlusRightBgCasc[i]*lScalingFactor*lScalingFactor + TMath::Power( lScalingFactorError*lLeftPlusRightBgCasc[i],2 ) );

       cout<<"Bin #"<<i<<", Replaced BG = "<<lLeftPlusRightBgCasc[i]<<"\t+/-\t"<<TMath::Sqrt(lLeftPlusRightBgCasc[i])<<" with "<<lNewValue<<"\t+/-\t"<<lNewValueError<<endl;
       lHistoBgTemplateMC[i]->Scale( lScalingFactor ) ; //Scaling histogram

       if ( fUseMCBackgroundTemplate ) {
    //Signal computation
    lSigRealCasc[i] = lSigPlusCenterBgCasc[i] - lNewValue;
    lSigErrRealCasc[i] = TMath::Sqrt(lSigPlusCenterBgCasc[i]+lNewValueError*lNewValueError);
    //Error Formula: New, has to be revised.

    //Attribute Raw Counts
    fHistPtRaw->SetBinContent(i+1, lSigRealCasc[i]/((fptbinlimits[i+1]-fptbinlimits[i])*lNEvents*(fRapidityBoundaryUpper-fRapidityBoundaryLower)) );
    fHistPtRaw->SetBinError(i+1, lSigErrRealCasc[i]/((fptbinlimits[i+1]-fptbinlimits[i])*lNEvents*(fRapidityBoundaryUpper-fRapidityBoundaryLower)) );

    //Attribute Raw Signal
    fHistPtSignal->SetBinContent(i+1, lSigRealCasc[i] );
    fHistPtSignal->SetBinError(i+1, lSigErrRealCasc[i] );

    //Signal-to-noise computation
    if( lLeftPlusRightBgCasc[i] != 0 ){
    fHistSigToNoise->SetBinContent(i+1, (lSigPlusCenterBgCasc[i] - lNewValue)/lNewValue );
    }else{
    fHistSigToNoise->SetBinContent(i+1, -1) ; //-1 means: no background
    }
    }

    }
    //=============================================================
     */

    //=============================================================
    // Compute Signal
    for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++){
        lSigMCCasc[ipoint] = lSigPlusCenterBgCascMC[ipoint] - lLeftPlusRightBgCascMC[ipoint];
        lSigErrMCCasc[ipoint] = TMath::Sqrt(lSigPlusCenterBgCascMC[ipoint]+lLeftPlusRightBgCascMC[ipoint]);
        //Error Formula: Equivalent to Sqrt(S+B+B) = Sqrt(S+2B)

        //Attribute Raw Counts
        fHistPtRawMC->SetBinContent(ipoint+1, lSigMCCasc[ipoint] / ((fptbinlimits[ipoint+1] - fptbinlimits[ipoint]) * lNEventsMC * (fRapidityBoundaryUpper - fRapidityBoundaryLower)) );
        fHistPtRawMC->SetBinError(ipoint+1, lSigErrMCCasc[ipoint] / ((fptbinlimits[ipoint+1] - fptbinlimits[ipoint]) * lNEventsMC * (fRapidityBoundaryUpper - fRapidityBoundaryLower)) );

        //Attribute Raw Counts
        fHistPtSignalMC->SetBinContent(ipoint+1, lSigMCCasc[ipoint] );
        fHistPtSignalMC->SetBinError(ipoint+1, lSigErrMCCasc[ipoint] );

        //Signal-to-noise computation
        if( lLeftPlusRightBgCascMC[ipoint] != 0 ){
            fHistSigToNoiseMC->SetBinContent(ipoint+1, (lSigPlusCenterBgCascMC[ipoint] - lLeftPlusRightBgCascMC[ipoint]) / lLeftPlusRightBgCascMC[ipoint] );
        }else{
            fHistSigToNoiseMC->SetBinContent(ipoint+1, -1) ; //-1 means: no background
        }
    }
    //=============================================================
    cout<<"--------------- Extracted Signal (MC) ------------------"<<endl;
    for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++){
        cout<<"---> ["<<fptbinlimits[ipoint]<<" - "<<fptbinlimits[ipoint+1]<<" GeV/c]\tSignal: "<<lSigMCCasc[ipoint]<<" +/- "<<lSigErrMCCasc[ipoint]<<endl;
    }
    cout<<"--------------------------------------------------------"<<endl;
    cout<<endl;


    cout<<"--------------- Process Generated Cascades -------------"<<endl;
    //Filling histogram with original MC particles====================
    TH1D* fHistGen = 0x0;
    cout<<"Getting Histogram..."<<endl;
    //------- OLD VERSION ----------------------------------------
    //if(fWhichParticle == "XiMinus")
    //  fHistGen 			= (TH1F*)clistMC->FindObject("fHistPt_GenXiMinus");
    //if(fWhichParticle == "XiPlus")
    //  fHistGen 			= (TH1F*)clistMC->FindObject("fHistPt_GenXiPlus");
    //if(fWhichParticle == "OmegaMinus")
    //  fHistGen 			= (TH1F*)clistMC->FindObject("fHistPt_GenOmegaMinus");
    //if(fWhichParticle == "OmegaPlus")
    //  fHistGen 			= (TH1F*)clistMC->FindObject("fHistPt_GenOmegaPlus");
    //
    //------- NEW VERSION ----------------------------------------
    TH3D* h3DGenerated = (TH3D*)clistMC->FindObject(Form("fHistGeneratedPtVsYVsCentrality%s", fWhichParticle.Data()));
    //
    // find projection bins in rapidity
    Int_t rapBin_min = h3DGenerated->GetYaxis()->FindBin( fRapidityBoundaryLower+1.e-6 );
    Int_t rapBin_max = h3DGenerated->GetYaxis()->FindBin( fRapidityBoundaryUpper-1.e-6 );
    //
    // find projection bins in multiplicity
    Int_t multBin_min = h3DGenerated->GetZaxis()->FindBin( 0.+1.e-6 );
    Int_t multBin_max = h3DGenerated->GetZaxis()->FindBin( 100.-1.e-6 );
    //
    // get the projection in pt
    fHistGen = (TH1D*)h3DGenerated->ProjectionX(Form("fHistPt_Gen%s", fWhichParticle.Data()), rapBin_min, rapBin_max, multBin_min, multBin_max);

    TH1F *fHistMCCountbyptCasc	= new TH1F("fHistMCCountbyptCasc","Cascade MC count;p_{T} (GeV/c);Counts", fptbinnumb, fptbinlimits);

    Double_t temppt;
    for(long i = 1; i<fHistGen->GetNbinsX()+1;i++){
        temppt = fHistGen->GetXaxis()->GetBinCenter(i);
        for(long filling = 0; filling<fHistGen->GetBinContent(i); filling++){
            fHistMCCountbyptCasc->Fill(temppt);
        }
    }
    cout<<"--------------- Generated Casc Dump --------------------"<<endl;
    for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++){
        cout<<"---> ["<<fptbinlimits[ipoint]<<" - "<<fptbinlimits[ipoint+1]<<" GeV/c]\tSignal: "<<fHistMCCountbyptCasc->GetBinContent(ipoint+1)<<" +/- "<<TMath::Sqrt(fHistMCCountbyptCasc->GetBinContent(ipoint+1))<<endl;
    }
    cout<<"--------------------------------------------------------"<<endl;
    cout<<endl;

    //=============================================================
    // Resolution tests
    if( fPerformResolutionTests ) {
        cout<<"--------------------------------------------------------"<<endl;
        cout<<" Resolution tests"<<endl;
        cout<<"--------------------------------------------------------"<<endl;
        cout<<" Usual tests: (Pt - PtMC) histograms "<<endl;
        cout<<" Fitting gaussians..."<<endl;
        for(Int_t ibin=0;ibin<fptbinnumb;ibin++){
            cout<<"---> Fitting Gaussian to bin #"<<ibin<<endl;
            lHistResolutionGaussian[ibin]->SetParameter(0, lHistResolution[ibin]->GetMaximum() );
            lHistResolutionGaussian[ibin]->SetParameter(1, lHistResolution[ibin]->GetMean()    );
            lHistResolutionGaussian[ibin]->SetParameter(2, lHistResolution[ibin]->GetRMS()     );
            sprintf(lHistResolutionName,"lHistResolutionGaussian%i",ibin);
            lHistResolution[ibin]->Fit(lHistResolutionName,"IREM0");
        }
    }

    //=============================================================
    // Compute Efficiency
    Double_t lEfficiency[100];
    Double_t lEfficiencyError[100];
    for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++){
        Double_t lEffNumerator = lSigMCCasc[ipoint];
        Double_t lEffDenominator = fHistMCCountbyptCasc->GetBinContent(ipoint+1);
        fHistEffNumerator->SetBinContent(ipoint+1, lEffNumerator);
        fHistEffDenominator->SetBinContent(ipoint+1, lEffDenominator);
        lEfficiency[ipoint] = lEffNumerator / lEffDenominator;
        lEfficiencyError[ipoint] = ErrorInRatio(
                lSigMCCasc[ipoint],
                lSigErrMCCasc[ipoint],
                fHistMCCountbyptCasc->GetBinContent(ipoint+1),
                TMath::Sqrt(fHistMCCountbyptCasc->GetBinContent(ipoint+1))
                );
    }
    //=============================================================
    cout<<endl;
    cout<<"--------------- Memory Cleanup -------------------------"<<endl;
    //fHistDummyCasc->Delete();
    fHistMCCountbyptCasc->Delete();
    clistMC->Delete();
    delete clistMC;
    fileMC->Close("R");
    fileMC->Delete();
    delete fileMC;
    cout<<endl;
    cout<<"--------------- Pure Efficiency Numbers ----------------"<<endl;
    TH1F* fHistPureEfficiency 		= new TH1F("fHistPureEfficiency", "Pure Efficiency;p_{T} (GeV/c);Pure Efficiency", fptbinnumb, fptbinlimits);
    TH1F* fHistEfficiency 		    = new TH1F("fHistEfficiency"    , "Efficiency;p_{T} (GeV/c);Efficiency"          , fptbinnumb, fptbinlimits);
    if( fWhichParticle == "XiMinus"      ) fHistPureEfficiency->SetTitle("#Xi^{-} Efficiency (no Geant-Fluka corr.)");
    if( fWhichParticle == "XiPlus"       ) fHistPureEfficiency->SetTitle("#bar{#Xi}^{+} Efficiency (no Geant-Fluka corr.)");
    if( fWhichParticle == "OmegaMinus"   ) fHistPureEfficiency->SetTitle("#Omega^{-} Efficiency (no Geant-Fluka corr.)");
    if( fWhichParticle == "OmegaPlus"    ) fHistPureEfficiency->SetTitle("#bar{#Omega}^{+} Efficiency (no Geant-Fluka corr.)");

    if( fWhichParticle == "XiMinus"      ) fHistEfficiency->SetTitle("#Xi^{-} Efficiency (with Geant-Fluka corr.)");
    if( fWhichParticle == "XiPlus"       ) fHistEfficiency->SetTitle("#bar{#Xi}^{+} Efficiency (with Geant-Fluka corr.)");
    if( fWhichParticle == "OmegaMinus"   ) fHistEfficiency->SetTitle("#Omega^{-} Efficiency (with Geant-Fluka corr.)");
    if( fWhichParticle == "OmegaPlus"    ) fHistEfficiency->SetTitle("#bar{#Omega}^{+} Efficiency (with Geant-Fluka corr.)");

    for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++){
        cout<<"---> ["<<fptbinlimits[ipoint]<<" - "<<fptbinlimits[ipoint+1]<<" GeV/c]\tEff.: "<<lEfficiency[ipoint]<<" +/- "<<lEfficiencyError[ipoint]<<endl;
        fHistPureEfficiency->SetBinContent(ipoint+1, lEfficiency[ipoint]);
        fHistPureEfficiency->SetBinError(ipoint+1, lEfficiencyError[ipoint]);
    }
    cout<<"--------------------------------------------------------"<<endl;
    cout<<endl;

    //do Geant-Fluka correction
    if( fFuncGeantFlukaCorr ) {
        cout<<"--------------- Geant-Fluka Correction ----------------"<<endl;
        cout<<"Correction received: "<<fFuncGeantFlukaCorr->GetName()<<endl;
        cout<<"Embedding Geant-Fluka in Efficiencies"<<endl;
        for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++) {
            printf("---> Xi pT bin: [%lf - %lf]\tProton pT: %.4lf GeV/c\n", fptbinlimits[ipoint], fptbinlimits[ipoint+1], lProtonMomentum[ipoint]->GetMean());
            lEfficiency[ipoint]      *= fFuncGeantFlukaCorr->Eval( lProtonMomentum[ipoint]->GetMean() ) ;
            lEfficiencyError[ipoint] *= fFuncGeantFlukaCorr->Eval( lProtonMomentum[ipoint]->GetMean() ) ;
        }
        cout<<"Geant-Fluka Corrected Efficiencies:"<<endl;
        for(Int_t ipoint = 0; ipoint<fptbinnumb; ipoint++) {
            cout<<"---> ["<<fptbinlimits[ipoint]<<" - "<<fptbinlimits[ipoint+1]<<" GeV/c]\tEff.: "<<lEfficiency[ipoint]<<" +/- "<<lEfficiencyError[ipoint]<<endl;
            fHistEfficiency->SetBinContent(ipoint+1, lEfficiency[ipoint]);
            fHistEfficiency->SetBinError(ipoint+1, lEfficiencyError[ipoint]);
        }
        cout<<"--------------------------------------------------------"<<endl;
        cout<<endl;
    } 
    else {
        cout<<"--------------------------------------------------------"<<endl;
        cout<<" WARNING - Geant-Fluka correction not being applied!"<<endl;
        cout<<"--------------------------------------------------------"<<endl;
        cout<<endl;
    }//end Geant-Fluka correction if

    //=========================================================================
    //---> At this stage, everthing's just ready for the actual spectrum
    //---> computation to occur!
    Double_t lSpectrum[100];
    Double_t lSpectrumError[100];

    for(Int_t ibin=0;ibin<fptbinnumb;ibin++){
        lSpectrum[ibin]      = (TMath::Abs(lEfficiency[ibin]>1.e-9)) ? (lSigRealCasc[ibin] / lEfficiency [ibin]) : 0.;
        lSpectrumError[ibin] = (TMath::Abs(lEfficiency[ibin]>1.e-9)) ? ErrorInRatio(
                lSigRealCasc    [ibin],
                lSigErrRealCasc [ibin],
                lEfficiency      [ibin],
                lEfficiencyError [ibin]
                ) : 0.;
    }

    //Divide by: Bin Width, Rapidity Window, N_{events}
    for(Int_t ibin=0;ibin<fptbinnumb;ibin++){
        lSpectrum[ibin] /= (fptbinlimits[ibin+1]-fptbinlimits[ibin]);
        lSpectrum[ibin] /= lNEvents;
        lSpectrum[ibin] /= (fRapidityBoundaryUpper-fRapidityBoundaryLower);
        lSpectrumError[ibin] /= (fptbinlimits[ibin+1]-fptbinlimits[ibin]);
        lSpectrumError[ibin] /= lNEvents;
        lSpectrumError[ibin] /= (fRapidityBoundaryUpper-fRapidityBoundaryLower);
    }

    TH1F* fHistPtXiMinus = new TH1F("fHistPtXiMinus","#Xi^{-} Corrected Spectrum;p_{T} (GeV/c);1/N #frac{d^{2}N}{dydp_{t}}",fptbinnumb,fptbinlimits);
    TH1F* fHistPtXiPlus = new TH1F("fHistPtXiPlus","#bar{#Xi}^{+} Corrected Spectrum;p_{T} (GeV/c);1/N #frac{d^{2}N}{dydp_{t}}",fptbinnumb,fptbinlimits);
    TH1F* fHistPtOmegaMinus = new TH1F("fHistPtOmegaMinus","#Omega^{-} Corrected Spectrum;p_{T} (GeV/c);1/N #frac{d^{2}N}{dydp_{t}}",fptbinnumb,fptbinlimits);
    TH1F* fHistPtOmegaPlus = new TH1F("fHistPtOmegaPlus","#bar{#Omega}^{+} Corrected Spectrum;p_{T} (GeV/c);1/N #frac{d^{2}N}{dydp_{t}}",fptbinnumb,fptbinlimits);

    //Copy to Histogram
    for(Int_t ibin=0;ibin<fptbinnumb;ibin++){
        if(fWhichParticle == "XiMinus" ){
            fHistPtXiMinus->SetBinContent( ibin+1, lSpectrum[ibin]      );
            fHistPtXiMinus->SetBinError  ( ibin+1, lSpectrumError[ibin] );
        }
        if(fWhichParticle == "XiPlus" ){
            fHistPtXiPlus->SetBinContent( ibin+1, lSpectrum[ibin]      );
            fHistPtXiPlus->SetBinError  ( ibin+1, lSpectrumError[ibin] );
        }
        if(fWhichParticle == "OmegaMinus" ){
            fHistPtOmegaMinus->SetBinContent( ibin+1, lSpectrum[ibin]      );
            fHistPtOmegaMinus->SetBinError  ( ibin+1, lSpectrumError[ibin] );
        }
        if(fWhichParticle == "OmegaPlus" ){
            fHistPtOmegaPlus->SetBinContent( ibin+1, lSpectrum[ibin]      );
            fHistPtOmegaPlus->SetBinError  ( ibin+1, lSpectrumError[ibin] );
        }
    }


    //=========================================================================

    cout<<"--------------- Result Output --------------------------"<<endl;
    cout<<" ---> Writing information to "<<fOutputDataFile<<endl;
    // Open an output file
    TFile* lResultsFile = TFile::Open(fOutputDataFile, "RECREATE");
    if (!lResultsFile || !lResultsFile->IsOpen()){
        cout<<"Error! Couldn't open file!"<<endl;
        return;
    }

    //Preparing Signal Extraction Range Canvas
    TCanvas *cSigExtRange = new TCanvas("cSigExtRange","Extraction Range",900,600);
    cSigExtRange->SetFillColor(kWhite);
    cSigExtRange->SetLeftMargin(0.17);
    cSigExtRange->SetRightMargin(0.17);
    cSigExtRange->cd();
    fHistSignalExtractionRange->SetFillColor(18);
    fHistSignalExtractionRange->SetMarkerStyle(20);
    if( fWhichParticle == "XiMinus" || fWhichParticle == "XiPlus" ){
        fHistSignalExtractionRange->GetYaxis()->SetRangeUser(1.322-0.08,1.322+0.08);
    }
    if( fWhichParticle == "OmegaMinus" || fWhichParticle == "OmegaPlus" ){
        fHistSignalExtractionRange->GetYaxis()->SetRangeUser(1.672-0.08,1.672+0.08);
    }
    fHistSignalExtractionRange->Draw("E2");

    //Saving Invariant Mass Plots (real data)
    lResultsFile->cd();
    TDirectoryFile *lInvMassReal = new TDirectoryFile("lInvMassReal","Invariant Mass Plots (Real Data)");
    lInvMassReal->cd();
    cSigExtRange->Write();
    if(!fUsePeakPositionAndWidthFromFit) {
        for(Int_t ibin = 0; ibin<fptbinnumb; ibin++) {
            lCanvasHistoMBCasc[ibin] -> Write();
        }
        for(Int_t ibin = 0; ibin<fptbinnumb; ibin++) {
            lCanvasHistoSelectedCasc[ibin] -> Write();
        }
    }

    TDirectoryFile *lInvMassRealRawData = new TDirectoryFile("lInvMassRealRawData","Objects for Inv Mass Plots (Real Data)");
    lInvMassRealRawData->cd();
    fHistSignalExtractionRange->Write();
    fHistPeakPosition->Write();
    fHistSigToNoise->Write();
    fHistPtRaw->Write();
    fHistPtSignal->Write();
    if(!fUsePeakPositionAndWidthFromFit) {
        for(Int_t ibin = 0; ibin<fptbinnumb; ibin++) {
            lHistoMBCasc[ibin] -> Write();
            fgausPt[ibin]      -> Write();
            if( fFitBackgroundSwitch ) lfitNoise[ibin] -> Write();
        }
    }
    for(Int_t ibin = 0; ibin<fptbinnumb; ibin++) {
        lHistoSelectedCasc[ibin] -> Write();
    }

    //Saving Invariant Mass Plots (MC)
    lResultsFile->cd();
    TDirectoryFile *lInvMassMC = new TDirectoryFile("lInvMassMC","Invariant Mass Plots (Monte Carlo)");
    lInvMassMC->cd();

    for(Int_t ibin = 0; ibin<fptbinnumb; ibin++) {
        lCanvasHistoSelectedCascMC[ibin] -> Write();
    }

    TDirectoryFile *lInvMassMCRawData = new TDirectoryFile("lInvMassMCRawData","Objects for Inv Mass Plots (MC)");
    lInvMassMCRawData->cd();
    fHistPtRawMC->Write();
    fHistPtSignalMC->Write();
    fHistPeakPositionMC->Write();
    fHistSigToNoiseMC->Write();
    f2dHistPtResolution->Write();
    f2dHistPtBlur->Write();
    f2dHistPtSharpen->Write();
    fHistPtGenerated->Write();
    fHistGenPerEvent->Write();

    for(Int_t ibin = 0; ibin<fptbinnumb; ibin++) {
        lHistoSelectedCascMC[ibin] ->Write();
        fgausPtMC[ibin]         ->Write();
        if( fFitBackgroundSwitch ) lfitNoiseMC[ibin] -> Write();
    }
    for(Int_t ibin = 0; ibin<fptbinnumb; ibin++) {
        lHistoBgTemplateMC[ibin]->Write();
    }

    //Saving Geant-Fluka Correction Data (MC)
    if( fFuncGeantFlukaCorr ) {
        lResultsFile->cd();
        TDirectoryFile *lGFCorrection = new TDirectoryFile("lGFCorrection","Geant-Fluka Correction Histograms");
        lGFCorrection->cd();

        fFuncGeantFlukaCorr->Write();
        for(Int_t ibin = 0; ibin<fptbinnumb; ibin++) {
            lProtonMomentum[ibin]->Write();
        }
    }

    fHistEffNumerator->Write();
    fHistEffDenominator->Write();

    //Saving pT-dep. cuts
    lResultsFile->cd();
    TDirectoryFile *lPtDepCuts = new TDirectoryFile("lPtDepCuts", "Objects for pT-dependent cuts");
    lPtDepCuts->cd();
    fListOfPtDepCuts->Write();

    //Saving Resolution Information
    //preparing...
    for(Int_t ibin=0; ibin<fptbinnumb; ibin++){
        fHistResolutionVsPt->SetBinContent(ibin+1,lHistResolution[ibin]->GetMean());
        fHistResolutionVsPt->SetBinError(ibin+1,lHistResolution[ibin]->GetRMS());
        fHistResolutionVsPtDivByBinWidth->SetBinContent(ibin+1,lHistResolution[ibin]->GetMean() / (fptbinlimits[ibin+1]-fptbinlimits[ibin]) );
        fHistResolutionVsPtDivByBinWidth->SetBinError(ibin+1,lHistResolution[ibin]->GetRMS() / (fptbinlimits[ibin+1]-fptbinlimits[ibin]) );

        fHistResolutionVsPtWithGaussians->SetBinContent(ibin+1,lHistResolutionGaussian[ibin]->GetParameter(1));
        fHistResolutionVsPtWithGaussians->SetBinError(ibin+1,lHistResolutionGaussian[ibin]->GetParameter(2));
        fHistResolutionVsPtDivByBinWidthWithGaussians->SetBinContent(ibin+1,lHistResolutionGaussian[ibin]->GetParameter(1) / (fptbinlimits[ibin+1]-fptbinlimits[ibin]) );
        fHistResolutionVsPtDivByBinWidthWithGaussians->SetBinError(ibin+1,lHistResolutionGaussian[ibin]->GetParameter(2) / (fptbinlimits[ibin+1]-fptbinlimits[ibin]) );
    }

    lResultsFile->cd();
    TDirectoryFile *lDirResolution = new TDirectoryFile("lInfoResolution","Resolution Information");
    lDirResolution->cd();
    fHistResolutionVsPt->Write();
    fHistResolutionVsPtDivByBinWidth->Write();
    fHistResolutionVsPtWithGaussians->Write();
    fHistResolutionVsPtDivByBinWidthWithGaussians->Write();
    for(Int_t ibin = 0; ibin<fptbinnumb; ibin++) {
        lHistResolution[ibin] ->Write();
        lHistResolutionGaussian[ibin]->Write();
    }

    //Saving distributions of selection variables
    lResultsFile->cd();
    TDirectoryFile *lDirVariables = new TDirectoryFile("lSelectionVariables", "Selection Variables");
    lDirVariables->cd();
    if(fSaveVarHistosSwitch) {
        lListOfVarHistosAll->Write();
        lListOfVarHistosSel->Write();
        lListOfVarHistosMCAll->Write();
        lListOfVarHistosMCSel->Write();
        lListOfVarHistosMCAssAll->Write();
        lListOfVarHistosMCAssSel->Write();
    }

    //save configuration histo
    lResultsFile->cd();
    fHistConfig->Write();

    //save histo with number of events analysed
    fHistEventCounter->Fill("# Events (Data)", lNEvents);
    fHistEventCounter->Fill("# Events (MC)",   lNEventsMC);
    fHistEventCounter->Write();

    fHistPureEfficiency->Write();
    fHistEfficiency->Write();

    if(fWhichParticle == "XiMinus"   ) fHistPtXiMinus     ->  Write();
    if(fWhichParticle == "XiPlus"    ) fHistPtXiPlus      ->  Write();
    if(fWhichParticle == "OmegaMinus") fHistPtOmegaMinus  ->  Write();
    if(fWhichParticle == "OmegaPlus" ) fHistPtOmegaPlus   ->  Write();

    lResultsFile->Write();
    lResultsFile->Close();
    delete lResultsFile;

    //================================================
    //Manual Cleanup of all created Histograms
    //================================================
    // needed if you want to re-run the whole thing without
    // memory leaks (systematics, etc)

    //switch on if you want large amounts of printout
    Bool_t lDebugCleaningProcess = kFALSE;

    if (lDebugCleaningProcess) cout<<"fHistConfig->Delete()"<<endl;
    fHistConfig->Delete();
    if (lDebugCleaningProcess) cout<<"fHistEventCounter->Delete()"<<endl;
    fHistEventCounter->Delete();
    if (lDebugCleaningProcess) cout<<"fHistPt->Delete()"<<endl;
    fHistPt->Delete();
    fHistPtRaw->Delete();
    fHistPtRawMC->Delete();
    fHistPtSignal->Delete();
    fHistPtSignalMC->Delete();
    fHistPtGenerated->Delete();
    fHistGenPerEvent->Delete();
    if (lDebugCleaningProcess) cout<<"fHistPeakPosition->Delete()"<<endl;
    fHistPeakPosition->Delete();
    if (lDebugCleaningProcess) cout<<"fHistPeakPositionMC->Delete()"<<endl;
    fHistPeakPositionMC->Delete();
    if (lDebugCleaningProcess) cout<<"fHistSigToNoise->Delete()"<<endl;
    fHistSigToNoise->Delete();
    if (lDebugCleaningProcess) cout<<"fHistSigToNoiseMC->Delete()"<<endl;
    fHistSigToNoiseMC->Delete();
    if (lDebugCleaningProcess) cout<<"fHistSignalExtractionRange->Delete()"<<endl;
    fHistSignalExtractionRange->Delete();
    if(lDebugCleaningProcess) cout<<"f2dHistPtResolution->Delete()"<<endl;
    f2dHistPtResolution->Delete();

    if(lDebugCleaningProcess) cout<<"lfitNoise*[*]->Delete(); lSampleNoise*[*]->Delete()"<<endl;
    if( fFitBackgroundSwitch ){
        for(long i=0; i<fptbinnumb; i++){
            lfitNoise[i]    -> Delete();
            lSampleNoise[i] -> Delete();
            lfitNoiseMC[i]    -> Delete();
            lSampleNoiseMC[i] -> Delete();
        }
    }

    //pt-by-pt histos
    fHistResolutionVsPt->Delete();
    fHistResolutionVsPtDivByBinWidth->Delete();
    fHistResolutionVsPtWithGaussians->Delete();
    fHistResolutionVsPtDivByBinWidthWithGaussians->Delete();
    f2dHistPtBlur->Delete();
    f2dHistPtSharpen->Delete();
    if (lDebugCleaningProcess) cout<<"lHistoMBCasc*[*]->Delete()"<<endl;
    for(Int_t ihist=0;ihist<100;ihist++){
        lHistoMBCasc[ihist]->Delete();
        lHistoSelectedCasc[ihist]->Delete();
        lHistoSelectedCascMC[ihist]->Delete();
        lHistoBgTemplateMC[ihist]->Delete();
        lHistResolution[ihist]->Delete();
        lHistResolutionGaussian[ihist]->Delete();
    }

    if (lDebugCleaningProcess) cout<<"lLine*[*]->Delete()"<<endl;
    //gaussian fit functions, drawing lines
    for(Int_t ibin = 0; ibin<fptbinnumb; ibin++){
        lLineLeftMost[ibin] ->Delete();
        lLineLeft[ibin]     ->Delete();
        lLineRight[ibin]    ->Delete();
        lLineRightMost[ibin]->Delete();
        if(!fUsePeakPositionAndWidthFromFit) fgausPt[ibin]->Delete();
        //
        lLineLeftMostMC[ibin] ->Delete();
        lLineLeftMC[ibin]     ->Delete();
        lLineRightMC[ibin]    ->Delete();
        lLineRightMostMC[ibin]->Delete();
        fgausPtMC[ibin]->Delete();
    }
    if (lDebugCleaningProcess) cout<<"lCanvasHistoMBCasc*[*]->Delete()"<<endl;
    for(Int_t ihist=0;ihist<100;ihist++){
        lCanvasHistoMBCasc[ihist]->Close();
        lCanvasHistoSelectedCasc[ihist]->Close();
        lCanvasHistoSelectedCascMC[ihist]->Close();
        delete lCanvasHistoMBCasc[ihist];
        delete lCanvasHistoSelectedCasc[ihist];
        delete lCanvasHistoSelectedCascMC[ihist];
    }

    if (lDebugCleaningProcess) cout<<"cSigExtRange->Delete()"<<endl;
    cSigExtRange->Close();
    delete cSigExtRange;

    if (lDebugCleaningProcess) cout<<"fHist*Efficiency->Delete()"<<endl;
    //MC info, efficiencies
    fHistPureEfficiency->Delete();
    fHistEfficiency->Delete();

    if (lDebugCleaningProcess) cout<<"fHistPt*->Delete()"<<endl;
    //Corrected Spectra Histograms
    fHistPtXiMinus->Delete();
    fHistPtXiPlus->Delete();
    fHistPtOmegaMinus->Delete();
    fHistPtOmegaPlus->Delete();


    //Exit Batch Mode
    gROOT->SetBatch (kFALSE);

    cout<<"--------------------------------------------------------"<<endl;
    cout<<" There, done! "<<endl;
    cout<<endl;




}

