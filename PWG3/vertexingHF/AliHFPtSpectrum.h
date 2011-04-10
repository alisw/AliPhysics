#ifndef ALIHFPTSPECTRUM_H
#define ALIHFPTSPECTRUM_H

/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

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
//  (the corrected yields per bin are divided by the bin-width)
//
// Author: Z.Conesa, zconesa@in2p3.fr
//***********************************************************************

#include "TNamed.h"
#include "TMath.h"

class TH1;
class TGraphAsymmErrors;
class AliHFSystErr;


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
  void SetMCptSpectra(TH1D *hDirect, TH1D *hFeedDown);
  // Set the theoretical feeddown pt spectrum
  void SetFeedDownMCptSpectra(TH1D *hFeedDown);
  // Set the theoretical direct & feeddown pt spectrum upper and lower bounds
  void SetMCptDistributionsBounds(TH1D *hDirectMax, TH1D *hDirectMin, TH1D *hFeedDownMax, TH1D *hFeedDownMin);
  // Set the theoretical feeddown pt spectrum upper and lower bounds
  void SetFeedDownMCptDistributionsBounds(TH1D *hFeedDownMax, TH1D *hFeedDownMin);
  // Set the acceptance and efficiency corrections for direct  
  void SetDirectAccEffCorrection(TH1D *hDirectEff);
  // Set the acceptance and efficiency corrections for direct & feeddown 
  void SetAccEffCorrection(TH1D *hDirectEff, TH1D *hFeedDownEff);
  // Set the reconstructed spectrum
  void SetReconstructedSpectrum(TH1D *hRec);
  void SetReconstructedSpectrumSystematics(TGraphAsymmErrors *gRec); 
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
  // Set global acceptance x efficiency correction uncertainty (in percentages)
  void SetAccEffPercentageUncertainty(Double_t globalEffUnc, Double_t globalBCEffRatioUnc){
    fGlobalEfficiencyUncertainties[0] = globalEffUnc;
    fGlobalEfficiencyUncertainties[1] = globalBCEffRatioUnc;
  }
  // Set the normalization factors
  void SetNormalization(Double_t normalization){
    fLuminosity[0]=normalization; fTrigEfficiency[0]=1.0;
  }
  void SetNormalization(Double_t nevents, Double_t sigma){
    fLuminosity[0]=nevents/sigma; fTrigEfficiency[0]=1.0;
  }
  void SetNormalization(Double_t nevents, Double_t sigma, Double_t sigmaunc){
    fLuminosity[0] = nevents/sigma; 
    fTrigEfficiency[0] = 1.0;
    fLuminosity[1] = fLuminosity[0] * TMath::Sqrt( (1/nevents) + (sigmaunc/sigma)*(sigmaunc/sigma) );
  }

  //
  // Getters
  //
  // Return the theoretical predictions used for the calculation (rebinned if needed)
  TH1D * GetDirectTheoreticalSpectrum() const { return (fhDirectMCpt ? (TH1D*)fhDirectMCpt : NULL); }
  TH1D * GetDirectTheoreticalUpperLimitSpectrum() const { return (fhDirectMCptMax ? (TH1D*)fhDirectMCptMax : NULL); }
  TH1D * GetDirectTheoreticalLowerLimitSpectrum() const { return (fhDirectMCptMin ? (TH1D*)fhDirectMCptMin : NULL); }
  TH1D * GetFeedDownTheoreticalSpectrum() const { return (fhFeedDownMCpt ? (TH1D*)fhFeedDownMCpt : NULL); }
  TH1D * GetFeedDownTheoreticalUpperLimitSpectrum() const { return (fhFeedDownMCptMax ? (TH1D*)fhFeedDownMCptMax : NULL); }
  TH1D * GetFeedDownTheoreticalLowerLimitSpectrum() const { return (fhFeedDownMCptMin ? (TH1D*)fhFeedDownMCptMin : NULL); }
  // Return the acceptance and efficiency corrections (rebinned if needed)
  TH1D * GetDirectAccEffCorrection() const { return (fhDirectEffpt ? (TH1D*)fhDirectEffpt : NULL); }
  TH1D * GetFeedDownAccEffCorrection() const { return (fhFeedDownEffpt ? (TH1D*)fhFeedDownEffpt : NULL); }
  // Return the TGraphAsymmErrors of the feed-down correction (extreme systematics)
  TGraphAsymmErrors * GetFeedDownCorrectionFcExtreme() const { return (fgFcExtreme ?  fgFcExtreme : NULL); }
  // Return the TGraphAsymmErrors of the feed-down correction (conservative systematics)
  TGraphAsymmErrors * GetFeedDownCorrectionFcConservative() const { return (fgFcConservative ?  fgFcConservative : NULL); }
  // Return the histogram of the feed-down correction
  TH1D * GetHistoFeedDownCorrectionFc() const { return (fhFc ?  (TH1D*)fhFc : NULL); }
  // Return the histograms of the feed-down correction bounds
  TH1D * GetHistoUpperLimitFeedDownCorrectionFc() const { return (fhFcMax ? (TH1D*)fhFcMax : NULL); }
  TH1D * GetHistoLowerLimitFeedDownCorrectionFc() const { return (fhFcMin ? (TH1D*)fhFcMin : NULL); }
  // Return the TGraphAsymmErrors of the yield after feed-down correction (systematics but feed-down) 
  TGraphAsymmErrors * GetFeedDownCorrectedSpectrum() const { return (fgYieldCorr ? fgYieldCorr : NULL); }
  // Return the TGraphAsymmErrors of the yield after feed-down correction (feed-down extreme systematics)
  TGraphAsymmErrors * GetFeedDownCorrectedSpectrumExtreme() const { return (fgYieldCorrExtreme ? fgYieldCorrExtreme : NULL); }
  // Return the TGraphAsymmErrors of the yield after feed-down correction (feed-down conservative systematics)
  TGraphAsymmErrors * GetFeedDownCorrectedSpectrumConservative() const { return (fgYieldCorrConservative ? fgYieldCorrConservative : NULL); }
  // Return the histogram of the yield after feed-down correction 
  TH1D * GetHistoFeedDownCorrectedSpectrum() const { return (fhYieldCorr ? (TH1D*)fhYieldCorr : NULL); }
  // Return the histogram of the yield after feed-down correction bounds
  TH1D * GetHistoUpperLimitFeedDownCorrectedSpectrum() const { return (fhYieldCorrMax ? (TH1D*)fhYieldCorrMax : NULL); }
  TH1D * GetHistoLowerLimitFeedDownCorrectedSpectrum() const { return (fhYieldCorrMin ? (TH1D*)fhYieldCorrMin : NULL); }
  // Return the equivalent invariant cross-section TGraphAsymmErrors (systematics but feed-down) 
  TGraphAsymmErrors * GetCrossSectionFromYieldSpectrum() const { return (fgSigmaCorr ? fgSigmaCorr : NULL); }
  // Return the equivalent invariant cross-section TGraphAsymmErrors (feed-down extreme systematics)
  TGraphAsymmErrors * GetCrossSectionFromYieldSpectrumExtreme() const { return (fgSigmaCorrExtreme ? fgSigmaCorrExtreme : NULL); }
  // Return the equivalent invariant cross-section TGraphAsymmErrors (feed-down conservative systematics)
  TGraphAsymmErrors * GetCrossSectionFromYieldSpectrumConservative() const { return (fgSigmaCorrConservative ? fgSigmaCorrConservative : NULL); }
  // Return the equivalent invariant cross-section histogram
  TH1D * GetHistoCrossSectionFromYieldSpectrum() const { return (fhSigmaCorr ? (TH1D*)fhSigmaCorr : NULL); }
  // Return the equivalent invariant cross-section histogram bounds
  TH1D * GetHistoUpperLimitCrossSectionFromYieldSpectrum() const { return (fhSigmaCorrMax ? (TH1D*)fhSigmaCorrMax : NULL); }
  TH1D * GetHistoLowerLimitCrossSectionFromYieldSpectrum() const { return (fhSigmaCorrMin ? (TH1D*)fhSigmaCorrMin : NULL); }

  //
  // Main function:
  //    Compute the invariant cross-section from the yield (correct it)
  // variables : analysed delta_y, BR for the final correction, BR b --> decay (relative to the input theoretical prediction)
  void ComputeHFPtSpectrum(Double_t deltaY=1.0, Double_t branchingRatioC=1.0, Double_t branchingRatioBintoFinalDecay=1.0);

  // Compute the systematic uncertainties
  //   taking as input the AliHFSystErr uncertainties
  void ComputeSystUncertainties(AliHFSystErr *systematics, Bool_t combineFeedDown);
  //
  // Drawing the corrected spectrum comparing to theoretical prediction
  void DrawSpectrum(TGraphAsymmErrors *gPrediction);

  //
  // Basic functions
  // 
  void EstimateAndSetDirectEfficiencyRecoBin(TH1D *hSimu, TH1D *hReco);
  void EstimateAndSetFeedDownEfficiencyRecoBin(TH1D *hSimu, TH1D *hReco);
  
  //
  // Functions to  reweight histograms for testing purposes: 
  //   to reweight the simulation: hToReweight is reweighted as hReference/hToReweight
  TH1D * ReweightHisto(TH1D *hToReweight, TH1D *hReference);
  //   to reweight the reco-histos: hRecToReweight is reweighted as hReference/hMCToReweight
  TH1D * ReweightRecHisto(TH1D *hRecToReweight, TH1D *hMCToReweight, TH1D *hMCReference);


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
  Bool_t CheckHistosConsistency(TH1D *h1, TH1D *h2);
  // Function to rebin the theoretical spectra in the data-reconstructed spectra binning
  TH1D * RebinTheoreticalSpectra(TH1D *hTheory, const char *name);
  // Function to estimate the efficiency in the data-reconstructed spectra binning
  TH1D * EstimateEfficiencyRecoBin(TH1D *hSimu, TH1D *hReco, const char *name);


  //
  // Input spectra
  //
  TH1D *fhDirectMCpt;            // Input MC c-->D spectra
  TH1D *fhFeedDownMCpt;          // Input MC b-->D spectra
  TH1D *fhDirectMCptMax;         // Input MC maximum c-->D spectra
  TH1D *fhDirectMCptMin;         // Input MC minimum c-->D spectra
  TH1D *fhFeedDownMCptMax;       // Input MC maximum b-->D spectra
  TH1D *fhFeedDownMCptMin;       // Input MC minimum b-->D spectra
  TH1D *fhDirectEffpt;           // c-->D Acceptance and efficiency correction
  TH1D *fhFeedDownEffpt;         // b-->D Acceptance and efficiency correction
  TH1D *fhRECpt;                 // all reconstructed D
  //
  TGraphAsymmErrors *fgRECSystematics; // all reconstructed D Systematic uncertainties
  //
  // Normalization factors
  Double_t fLuminosity[2];           // analyzed luminosity & uncertainty
  Double_t fTrigEfficiency[2];       // trigger efficiency & uncertainty
  Double_t fGlobalEfficiencyUncertainties[2]; // uncertainties on the efficiency [0]=c, b, [1]=b/c

  //
  // Output spectra
  //
  TH1D *fhFc;                            // Correction histo fc = 1 / ( 1 + (eff_b/eff_c)*(N_b/N_c) ) 
  TH1D *fhFcMax;                         // Maximum fc histo
  TH1D *fhFcMin;                         // Minimum fc histo
  TGraphAsymmErrors * fgFcExtreme;       // Extreme correction as TGraphAsymmErrors
  TGraphAsymmErrors * fgFcConservative;  // Extreme correction as TGraphAsymmErrors
  TH1D *fhYieldCorr;                     // Corrected yield (stat unc. only)
  TH1D *fhYieldCorrMax;                  // Maximum corrected yield  
  TH1D *fhYieldCorrMin;                  // Minimum corrected yield  
  TGraphAsymmErrors * fgYieldCorr;              // Corrected yield as TGraphAsymmErrors  (syst but feed-down)
  TGraphAsymmErrors * fgYieldCorrExtreme;       // Extreme corrected yield as TGraphAsymmErrors  (syst from feed-down)
  TGraphAsymmErrors * fgYieldCorrConservative;  // Conservative corrected yield as TGraphAsymmErrors  (syst from feed-down) 
  TH1D *fhSigmaCorr;                     // Corrected cross-section (stat unc. only)
  TH1D *fhSigmaCorrMax;                  // Maximum corrected cross-section  
  TH1D *fhSigmaCorrMin;                  // Minimum corrected cross-section
  TGraphAsymmErrors * fgSigmaCorr;              // Corrected cross-section as TGraphAsymmErrors (syst but feed-down)
  TGraphAsymmErrors * fgSigmaCorrExtreme;       // Extreme corrected cross-section as TGraphAsymmErrors (syst from feed-down)
  TGraphAsymmErrors * fgSigmaCorrConservative;  // Conservative corrected cross-section as TGraphAsymmErrors  (syst from feed-down)

  //
  Int_t fFeedDownOption;            // feed-down correction flag: 0=none, 1=fc, 2=Nb 
  Bool_t fAsymUncertainties;        // flag: asymmetric uncertainties are (1) or not (0) considered


  ClassDef(AliHFPtSpectrum,1) // Class for Heavy Flavor spectra corrections
};

#endif
