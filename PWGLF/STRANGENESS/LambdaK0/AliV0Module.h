#ifndef ALIV0MODULE_H
#define ALIV0MODULE_H
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

class AliV0Module{
public: 
  //Constructor
  AliV0Module(); 
  AliV0Module(TString ParticleType);

  //Set Files to Use
  void SetRealDataFile      ( TString RealDataFilename     );
  void SetMCDataFile        ( TString MCDataFilename       );
  void SetFeedDownDataFile  ( TString FeedDownDataFilename );
  void SetOutputFile        ( TString OutputFilename       );

  //Set Pt Bin Limits
  void SetPtBinLimits(Long_t got_ptbinnumb, const Double_t *got_ptbinlimits);

  //Set Rapidity Window
  void SetRapidityWindow(Double_t got_RapidityBoundary);

  //Set CINT1B/INEL to normalize to yield
  void SetCINT1BoverINEL(Double_t got_CINT1BoverINEL);

  //Set Cuts - topological
  void SetCutV0Radius       (Double_t cut);
  void SetCutDCANegToPV     (Double_t cut);
  void SetCutDCAPosToPV     (Double_t cut);
  void SetCutDCAV0Daughters (Double_t cut);
  void SetCutV0CosPA        (Double_t cut);

  //Set Cuts - other
  void SetCutProperLifetime                         (Double_t cut);
  void SetCutTPCPIDNSigmas                          (Double_t cut);
  void SetCutSigmaForSignalExtraction               (Double_t cut);
  void SetCutLeastNumberOfCrossedRows               (Double_t cut);
  void SetCutLeastNumberOfCrossedRowsOverFindable   (Double_t cut);
  void SetCutDaughterEta                            (Double_t cut);
  void SetCutCompetingV0Rejection                   (Double_t cut);

  //Set Feeddown treatment
  void SetFeeddownTreatment ( TString FDMethod );

  //Set Fit Background or not 
  void SetFitBackground ( Bool_t fitBgSwitch );

  //Do Analysis
  void DoAnalysis();

  //Set Default Cuts
  void SetDefaultCuts();

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

private:
  //root file names
  TString fRealDataFile;
  TString fMCDataFile;
  TString fFeedDownDataFile;
  TString fOutputDataFile;

  //Store Pt Bin Limits
  //Max Number of Pt Bins Set here (100)
  Double_t fptbinlimits[100];
  Double_t fptX[100];
  Long_t fptbinnumb;

  //Rapidity Range Window
  Double_t fRapidityBoundary;

  //CINT1B/INEL for normalization
  Double_t fCINT1BoverINELratio; 
  
  //Main Analysis Parameters
  //--- 5 Topological Selections
  Double_t fCutV0Radius;
  Double_t fCutDCANegToPV;
  Double_t fCutDCAPosToPV;
  Double_t fCutDCAV0Daughters;
  Double_t fCutV0CosPA;
  //--- Proper Lifetime
  Double_t fCutProperLifetime;
  //--- TPC dE/dx N-sigmas
  Double_t fCutTPCPIDNSigmas;
  //--- Sigmas for signal extraction
  Double_t fCutNSigmasForSignalExtraction;
  //--- Smallest Number of Crossed Rows in TPC accepted for tracks
  Double_t fCutLeastNumberOfCrossedRows;
  Double_t fCutLeastNumberOfCrossedRowsOverFindable;
  //--- Daughter Track eta cut 
  Double_t fCutDaughterEta;
  //--- Competing V0 Species Rejection
  Double_t fCutCompetingV0Rejection;

  //WhichParticle switch: "Lambda", "AntiLambda" or "K0Short"
  TString fWhichParticle;

  //Do Fitting instead of bin counting switch
  Bool_t fFitBackgroundSwitch;

  //FeedDownTreatment switch:  
  // --- NoFD.............: Doesn't feeddown subtract at all 
  // --- DoubleChargedXi..: Multiply Charged Xi subtraction by 2
  // --- UseMCRatio.......: Fill FD Matrix with Xi- and Xi0 
  //  (effectively use Xi0/Xi- from MC, should not be used with 
  //   Xi- and Xi+ enhanced datasets!!)
  TString fFDSwitch;
};
#endif
