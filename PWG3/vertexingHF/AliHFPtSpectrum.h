#ifndef ALIHFPTSPECTRUM_H
#define ALIHFPTSPECTRUM_H

/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//***********************************************************************
// Class AliHFPtSpectrum
// Base class for feed-down corrections on heavy-flavour decays
// computes the cross-section via one of the three implemented methods:
//   0) Consider no feed-down prediction 
//   1) Subtract the feed-down with the "fc" method 
//       Yield = Reco * fc;  where fc = 1 / ( 1 + (eff_b/eff_c)*(N_b/N_c) ) ;
//   2) Subtract the feed-down with the "Nb" method
//       Yield = Reco - Feed-down (exact formula on the function implementation)
//
// Author: Z.Conesa, zconesa@in2p3.fr
//***********************************************************************

#include "TNamed.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"

class AliHFPtSpectrum: public TNamed
{

 public:
  
  // Constructor
  AliHFPtSpectrum(const char* name="AliHFPtSpectrum", const char* title="HF feed down correction class", Int_t option=1);
  // Copy constructor
  AliHFPtSpectrum(const AliHFPtSpectrum &rhs);
  // Assignment operator
  AliHFPtSpectrum& operator=(const AliHFPtSpectrum &source);
  // Destructor
  virtual ~AliHFPtSpectrum();

  //
  // Setters
  //
  // Set the theoretical direct & feeddown pt spectrum
  void SetMCptSpectra(TH1 *hDirect, TH1 *hFeedDown);
  // Set the theoretical feeddown pt spectrum
  void SetFeedDownMCptSpectra(TH1 *hFeedDown);
  // Set the theoretical direct & feeddown pt spectrum upper and lower bounds
  void SetMCptDistributionsBounds(TH1 *hDirectMax, TH1 *hDirectMin, TH1 *hFeedDownMax, TH1 *hFeedDownMin);
  // Set the theoretical feeddown pt spectrum upper and lower bounds
  void SetFeedDownMCptDistributionsBounds(TH1 *hFeedDownMax, TH1 *hFeedDownMin);
  // Set the acceptance and efficiency corrections for direct  
  void SetDirectAccEffCorrection(TH1 *hDirectEff);
  // Set the acceptance and efficiency corrections for direct & feeddown 
  void SetAccEffCorrection(TH1 *hDirectEff, TH1 *hFeedDownEff);
  // Set the reconstructed spectrum
  void SetReconstructedSpectrum(TH1 *hRec);
  // Set the calculation option flag for feed-down correction: 0=none, 1=fc , 2=Nb 
  void SetFeedDownCalculationOption(Int_t option){ fFeedDownOption = option; }
  // Set if the calculation has to consider asymmetric uncertaInt_ties or not
  void SetComputeAsymmetricUncertainties(bool flag){ fAsymUncertainties = flag; }
  // Set the luminosity and its uncertainty
  void SetLuminosity(Double_t luminosity, Double_t unc){
    fLuminosity[0]=luminosity;  fLuminosity[1]=unc;
  }
  // Set the trigger efficiency and its uncertainty
  void SetTriggerEfficiency(Double_t efficiency, Double_t unc){
    fTrigEfficiency[0]=efficiency; fTrigEfficiency[1]=unc;
  }

  //
  // Getters
  //
  // Return the TGraphAsymmErrors of the feed-down correction
  TGraphAsymmErrors * GetFeedDownCorrection_fc() { return (fgFc ?  fgFc : NULL); }
  // Return the histogram of the feed-down correction
  TH1 * GetHistoFeedDownCorrection_fc() { return (fhFc ?  fhFc : NULL); }
  // Return the histograms of the feed-down correction bounds
  TH1 * GetHistoUpperLimitFeedDownCorrection_fc() { return (fhFc_max ? fhFc_max : NULL); }
  TH1 * GetHistoLowerLimitFeedDownCorrection_fc() { return (fhFc_min ? fhFc_min : NULL); }
  // Return the TGraphAsymmErrors of the yield after feed-down correction 
  TGraphAsymmErrors * GetFeedDownCorrectedSpectrum() { return (fgYieldCorr ? fgYieldCorr : NULL); }
  // Return the histogram of the yield after feed-down correction 
  TH1 * GetHistoFeedDownCorrectedSpectrum() { return (fhYieldCorr ? fhYieldCorr : NULL); }
  // Return the histogram of the yield after feed-down correction bounds
  TH1 * GetHistoUpperLimitFeedDownCorrectedSpectrum() { return (fhYieldCorr_max ? fhYieldCorr_max : NULL); }
  TH1 * GetHistoLowerLimitFeedDownCorrectedSpectrum() { return (fhYieldCorr_min ? fhYieldCorr_min : NULL); }
  // Return the equivalent invariant cross-section TGraphAsymmErrors
  TGraphAsymmErrors * GetCrossSectionFromYieldSpectrum() { return (fgSigmaCorr ? fgSigmaCorr : NULL); }
  // Return the equivalent invariant cross-section histogram
  TH1 * GetHistoCrossSectionFromYieldSpectrum() { return (fhSigmaCorr ? fhSigmaCorr : NULL); }
  // Return the equivalent invariant cross-section histogram bounds
  TH1 * GetHistoUpperLimitCrossSectionFromYieldSpectrum() { return (fhSigmaCorr_max ? fhSigmaCorr_max : NULL); }
  TH1 * GetHistoLowerLimitCrossSectionFromYieldSpectrum() { return (fhSigmaCorr_min ? fhSigmaCorr_min : NULL); }

  //
  // Main function:
  //    Compute the invariant cross-section from the yield (correct it)
  void ComputeHFPtSpectrum(Double_t delta_y=1.0, Double_t BR_c=1.0, Double_t BR_b_decay=1.0);

  //
  // Functions to  reweight histograms for testing purposes: 
  //   to reweight the simulation: hToReweight is reweighted as hReference/hToReweight
  TH1 * ReweightHisto(TH1 *hToReweight, TH1 *hReference);
  //   to reweight the reco-histos: hRecToReweight is reweighted as hReference/hMCToReweight
  TH1 * ReweightRecHisto(TH1 *hRecToReweight, TH1 *hMCToReweight, TH1 *hMCReference);


 protected:

  // Initialization 
  Bool_t Initialize();
  
  // Basic functions
  //
  // Compute the feed-down correction via fc-method
  void CalculateFeedDownCorrection_fc(); 
  // Correct the yield for feed-down correction via fc-method
  void CalculateFeedDownCorrectedSpectrum_fc(); 
  // Correct the yield for feed-down correction via Nb-method
  void CalculateFeedDownCorrectedSpectrum_Nb(Float_t delta_y, Double_t BR_b_decay); 

  // Check histograms consistency function
  Bool_t CheckHistosConsistency(TH1 *h1, TH1 *h2);

  //
  // Input spectra
  //
  TH1 *fhDirectMCpt;            // Input MC c-->D spectra
  TH1 *fhFeedDownMCpt;          // Input MC b-->D spectra
  TH1 *fhDirectMCpt_max;        // Input MC maximum c-->D spectra
  TH1 *fhDirectMCpt_min;        // Input MC minimum c-->D spectra
  TH1 *fhFeedDownMCpt_max;      // Input MC maximum b-->D spectra
  TH1 *fhFeedDownMCpt_min;      // Input MC minimum b-->D spectra
  TH1 *fhDirectEffpt;           // c-->D Acceptance and efficiency correction
  TH1 *fhFeedDownEffpt;         // b-->D Acceptance and efficiency correction
  TH1 *fhRECpt;                 // all reconstructed D
  //
  // Normalization factors
  Double_t fLuminosity[2];           // analyzed luminosity & uncertainty
  Double_t fTrigEfficiency[2];       // trigger efficiency & uncertainty

  //
  // Output spectra
  //
  TH1 *fhFc;                            // Correction histo fc = 1 / ( 1 + (eff_b/eff_c)*(N_b/N_c) ) 
  TH1 *fhFc_max;                        // Maximum fc histo
  TH1 *fhFc_min;                        // Minimum fc histo
  TGraphAsymmErrors * fgFc;             // Correction as TGraphAsymmErrors
  TH1 *fhYieldCorr;                     // Corrected yield  
  TH1 *fhYieldCorr_max;                 // Maximum corrected yield  
  TH1 *fhYieldCorr_min;                 // Minimum corrected yield  
  TGraphAsymmErrors * fgYieldCorr;      // Corrected yield as TGraphAsymmErrors
  TH1 *fhSigmaCorr;                     // Corrected cross-section  
  TH1 *fhSigmaCorr_max;                 // Maximum corrected cross-section  
  TH1 *fhSigmaCorr_min;                 // Minimum corrected cross-section
  TGraphAsymmErrors * fgSigmaCorr;      // Corrected cross-section as TGraphAsymmErrors

  //
  Int_t fFeedDownOption;            // feed-down correction flag: 0=none, 1=fc, 2=Nb 
  Bool_t fAsymUncertainties;        // flag: asymmetric uncertainties are (1) or not (0) considered


  ClassDef(AliHFPtSpectrum,1) // Class for Heavy Flavor spectra corrections
};

#endif
