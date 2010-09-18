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
  void SetComputeAsymmetricUncertainties(Bool_t flag){ fAsymUncertainties = flag; }
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
  TGraphAsymmErrors * GetFeedDownCorrectionFc() { return (fgFc ?  fgFc : NULL); }
  // Return the histogram of the feed-down correction
  TH1 * GetHistoFeedDownCorrectionFc() { return (fhFc ?  fhFc : NULL); }
  // Return the histograms of the feed-down correction bounds
  TH1 * GetHistoUpperLimitFeedDownCorrectionFc() { return (fhFcMax ? fhFcMax : NULL); }
  TH1 * GetHistoLowerLimitFeedDownCorrectionFc() { return (fhFcMin ? fhFcMin : NULL); }
  // Return the TGraphAsymmErrors of the yield after feed-down correction 
  TGraphAsymmErrors * GetFeedDownCorrectedSpectrum() { return (fgYieldCorr ? fgYieldCorr : NULL); }
  // Return the histogram of the yield after feed-down correction 
  TH1 * GetHistoFeedDownCorrectedSpectrum() { return (fhYieldCorr ? fhYieldCorr : NULL); }
  // Return the histogram of the yield after feed-down correction bounds
  TH1 * GetHistoUpperLimitFeedDownCorrectedSpectrum() { return (fhYieldCorrMax ? fhYieldCorrMax : NULL); }
  TH1 * GetHistoLowerLimitFeedDownCorrectedSpectrum() { return (fhYieldCorrMin ? fhYieldCorrMin : NULL); }
  // Return the equivalent invariant cross-section TGraphAsymmErrors
  TGraphAsymmErrors * GetCrossSectionFromYieldSpectrum() { return (fgSigmaCorr ? fgSigmaCorr : NULL); }
  // Return the equivalent invariant cross-section histogram
  TH1 * GetHistoCrossSectionFromYieldSpectrum() { return (fhSigmaCorr ? fhSigmaCorr : NULL); }
  // Return the equivalent invariant cross-section histogram bounds
  TH1 * GetHistoUpperLimitCrossSectionFromYieldSpectrum() { return (fhSigmaCorrMax ? fhSigmaCorrMax : NULL); }
  TH1 * GetHistoLowerLimitCrossSectionFromYieldSpectrum() { return (fhSigmaCorrMin ? fhSigmaCorrMin : NULL); }

  //
  // Main function:
  //    Compute the invariant cross-section from the yield (correct it)
  // variables : analysed delta_y, BR for the final correction, BR b --> decay (relative to the input theoretical prediction)
  void ComputeHFPtSpectrum(Double_t deltaY=1.0, Double_t branchingRatioC=1.0, Double_t branchingRatioBintoFinalDecay=1.0);

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
  void CalculateFeedDownCorrectionFc(); 
  // Correct the yield for feed-down correction via fc-method
  void CalculateFeedDownCorrectedSpectrumFc(); 
  // Correct the yield for feed-down correction via Nb-method
  void CalculateFeedDownCorrectedSpectrumNb(Double_t deltaY, Double_t branchingRatioBintoFinalDecay); 

  // Check histograms consistency function
  Bool_t CheckHistosConsistency(TH1 *h1, TH1 *h2);

  //
  // Input spectra
  //
  TH1 *fhDirectMCpt;            // Input MC c-->D spectra
  TH1 *fhFeedDownMCpt;          // Input MC b-->D spectra
  TH1 *fhDirectMCptMax;         // Input MC maximum c-->D spectra
  TH1 *fhDirectMCptMin;         // Input MC minimum c-->D spectra
  TH1 *fhFeedDownMCptMax;       // Input MC maximum b-->D spectra
  TH1 *fhFeedDownMCptMin;       // Input MC minimum b-->D spectra
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
  TH1 *fhFcMax;                         // Maximum fc histo
  TH1 *fhFcMin;                         // Minimum fc histo
  TGraphAsymmErrors * fgFc;             // Correction as TGraphAsymmErrors
  TH1 *fhYieldCorr;                     // Corrected yield  
  TH1 *fhYieldCorrMax;                  // Maximum corrected yield  
  TH1 *fhYieldCorrMin;                  // Minimum corrected yield  
  TGraphAsymmErrors * fgYieldCorr;      // Corrected yield as TGraphAsymmErrors
  TH1 *fhSigmaCorr;                     // Corrected cross-section  
  TH1 *fhSigmaCorrMax;                  // Maximum corrected cross-section  
  TH1 *fhSigmaCorrMin;                  // Minimum corrected cross-section
  TGraphAsymmErrors * fgSigmaCorr;      // Corrected cross-section as TGraphAsymmErrors

  //
  Int_t fFeedDownOption;            // feed-down correction flag: 0=none, 1=fc, 2=Nb 
  Bool_t fAsymUncertainties;        // flag: asymmetric uncertainties are (1) or not (0) considered


  ClassDef(AliHFPtSpectrum,1) // Class for Heavy Flavor spectra corrections
};

#endif
