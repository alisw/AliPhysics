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
//
//  In HIC you can also evaluate how the feed-down correction is influenced by an energy loss hypothesis: 
//      Raa(c-->D) / Raa(b-->D) defined here as Rcb for the "fc" method
//      Raa(b-->D) defined here as Rb for the "Nb" method
//
// Author: Z.Conesa, zconesa@in2p3.fr
//***********************************************************************

#include "TNamed.h"
#include "TMath.h"

class TH1;
class TH2;
class TNtuple;
class TGraphAsymmErrors;


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
  // Set if the calculation has to consider Ratio(c/b eloss) hypothesis 
  void SetComputeElossHypothesis(Bool_t flag){ fPbPbElossHypothesis = flag; }
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
    fLuminosity[0]=normalization;
  }
  void SetNormalization(Int_t nevents, Double_t sigma){
    fLuminosity[0]=nevents/sigma;
    fNevts = nevents;
  }
  void SetNormalization(Int_t nevents, Double_t sigma, Double_t sigmaunc){
    fLuminosity[0] = nevents/sigma; 
    fLuminosity[1] = fLuminosity[0] * TMath::Sqrt( (1/nevents) + (sigmaunc/sigma)*(sigmaunc/sigma) );
    fNevts = nevents;
  }
  //
  // Set the Tab parameter and its uncertainty
  void SetTabParameter(Double_t tabvalue, Double_t uncertainty){
    fTab[0] = tabvalue;
    fTab[1] = uncertainty;
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
  // Return whether the Ratio(c/b eloss) hypothesis has been considered
  Bool_t IsElossHypothesisCalculated(){ return fPbPbElossHypothesis; }
  // Return the TGraphAsymmErrors of the feed-down correction (extreme systematics)
  TGraphAsymmErrors * GetFeedDownCorrectionFcExtreme() const { return (fgFcExtreme ?  fgFcExtreme : NULL); }
  // Return the TGraphAsymmErrors of the feed-down correction (conservative systematics)
  TGraphAsymmErrors * GetFeedDownCorrectionFcConservative() const { return (fgFcConservative ?  fgFcConservative : NULL); }
  // Return the histogram of the feed-down correction
  TH1D * GetHistoFeedDownCorrectionFc() const { return (fhFc ?  (TH1D*)fhFc : NULL); }
  // Return the histograms of the feed-down correction bounds
  TH1D * GetHistoUpperLimitFeedDownCorrectionFc() const { return (fhFcMax ? (TH1D*)fhFcMax : NULL); }
  TH1D * GetHistoLowerLimitFeedDownCorrectionFc() const { return (fhFcMin ? (TH1D*)fhFcMin : NULL); }
  // Return the histogram of the feed-down correction times the Ratio(c/b eloss)
  TH2D * GetHistoFeedDownCorrectionFcVsEloss() const { return (fhFcRcb ?  (TH2D*)fhFcRcb : NULL); }
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
  // Return the histogram of the yield after feed-down correction vs the Ratio(c/b eloss)
  TH2D * GetHistoFeedDownCorrectedSpectrumVsEloss() const { return (fhYieldCorrRcb ? (TH2D*)fhYieldCorrRcb : NULL); }
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
  // Return the cross section systematics from data systematics
  TH1D * GetHistoCrossSectionDataSystematics() const { return (fhSigmaCorrDataSyst ? (TH1D*)fhSigmaCorrDataSyst : NULL); }
  //
  // PbPb special calculations 
  // Return the equivalent invariant cross-section histogram vs the Ratio(c/b eloss)
  TH2D * GetHistoCrossSectionFromYieldSpectrumVsEloss() const { return (fhSigmaCorrRcb ? (TH2D*)fhSigmaCorrRcb : NULL); }
  // Return the ntuple of the calculation vs the Ratio(c/b eloss)
  TNtuple * GetNtupleCrossSectionVsEloss() { return (fnSigma ? (TNtuple*)fnSigma : NULL); }
  //
  //
  // Histograms to keep track of the influence of the efficiencies statistical uncertainty on the cross-section
  TH1D * GetDirectStatEffUncOnSigma() const { return (TH1D*)fhStatUncEffcSigma; }
  TH1D * GetFeedDownStatEffUncOnSigma() const { return (TH1D*)fhStatUncEffbSigma; }
  // Histograms to keep track of the influence of the efficiencies statistical uncertainty on the feed-down correction factor
  TH1D * GetDirectStatEffUncOnFc() const { return (TH1D*)fhStatUncEffcFD; }
  TH1D * GetFeedDownStatEffUncOnFc() const { return (TH1D*)fhStatUncEffbFD; }


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
  // Functionality to find the y-axis bin of a TH2 for a given y-value
  Int_t FindTH2YBin(TH2D *histo, Float_t yvalue);


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
  Int_t fNevts;                      // nb of analyzed events
  Double_t fLuminosity[2];           // analyzed luminosity & uncertainty
  Double_t fTrigEfficiency[2];       // trigger efficiency & uncertainty
  Double_t fGlobalEfficiencyUncertainties[2]; // uncertainties on the efficiency [0]=c, b, [1]=b/c
  Double_t fTab[2];                   // Tab parameter and its uncertainty

  //
  // Output spectra
  //
  TH1D *fhFc;                            // Correction histo fc = 1 / ( 1 + (eff_b/eff_c)*(N_b/N_c) ) 
  TH1D *fhFcMax;                         // Maximum fc histo
  TH1D *fhFcMin;                         // Minimum fc histo
  TH2D *fhFcRcb;                         // Correction histo fc vs the Ratio(c/b eloss)
  TGraphAsymmErrors * fgFcExtreme;       // Extreme correction as TGraphAsymmErrors
  TGraphAsymmErrors * fgFcConservative;  // Extreme correction as TGraphAsymmErrors
  TH1D *fhYieldCorr;                     // Corrected yield (stat unc. only)
  TH1D *fhYieldCorrMax;                  // Maximum corrected yield  
  TH1D *fhYieldCorrMin;                  // Minimum corrected yield  
  TH2D *fhYieldCorrRcb;                  // Corrected yield (stat unc. only) vs the Ratio(c/b eloss)
  TGraphAsymmErrors * fgYieldCorr;              // Corrected yield as TGraphAsymmErrors  (syst but feed-down)
  TGraphAsymmErrors * fgYieldCorrExtreme;       // Extreme corrected yield as TGraphAsymmErrors  (syst from feed-down)
  TGraphAsymmErrors * fgYieldCorrConservative;  // Conservative corrected yield as TGraphAsymmErrors  (syst from feed-down) 
  TH1D *fhSigmaCorr;                     // Corrected cross-section (stat unc. only)
  TH1D *fhSigmaCorrMax;                  // Maximum corrected cross-section  
  TH1D *fhSigmaCorrMin;                  // Minimum corrected cross-section
  TH1D *fhSigmaCorrDataSyst;             // Corrected cross-section (syst. unc. from data only)
  TH2D *fhSigmaCorrRcb;                  // Corrected cross-section (stat unc. only) vs the Ratio(c/b eloss)
  TGraphAsymmErrors * fgSigmaCorr;              // Corrected cross-section as TGraphAsymmErrors (syst but feed-down)
  TGraphAsymmErrors * fgSigmaCorrExtreme;       // Extreme corrected cross-section as TGraphAsymmErrors (syst from feed-down)
  TGraphAsymmErrors * fgSigmaCorrConservative;  // Conservative corrected cross-section as TGraphAsymmErrors  (syst from feed-down)
  //
  TNtuple *fnSigma;     // Ntuple of the calculation vs the Ratio(c/b eloss)

  //
  Int_t fFeedDownOption;            // feed-down correction flag: 0=none, 1=fc, 2=Nb 
  Bool_t fAsymUncertainties;        // flag: asymmetric uncertainties are (1) or not (0) considered
  Bool_t fPbPbElossHypothesis;      // flag: whether to do estimates vs Ratio(c/b eloss) hypothesis

  //
  TH1D *fhStatUncEffcSigma;       // Uncertainty on the cross-section due to the prompt efficiency statistical uncertainty
  TH1D *fhStatUncEffbSigma;       // Uncertainty on the cross-section due to the feed-down efficiency statistical uncertainty
  TH1D *fhStatUncEffcFD;          // Uncertainty on the feed-down correction due to the prompt efficiency statistical uncertainty
  TH1D *fhStatUncEffbFD;          // Uncertainty on the feed-down correction due to the feed-down efficiency statistical uncertainty

  ClassDef(AliHFPtSpectrum,2) // Class for Heavy Flavor spectra corrections
};

#endif
