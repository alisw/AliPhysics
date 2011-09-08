/**************************************************************************
 * Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
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

#include <Riostream.h>

#include "TMath.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2.h"
#include "TH2D.h"
#include "TNtuple.h"
#include "TGraphAsymmErrors.h"
#include "TNamed.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "AliLog.h"
#include "AliHFSystErr.h"
#include "AliHFPtSpectrum.h"

ClassImp(AliHFPtSpectrum)

//_________________________________________________________________________________________________________
AliHFPtSpectrum::AliHFPtSpectrum(const char* name, const char* title, Int_t option):
  TNamed(name,title),
  fhDirectMCpt(NULL),
  fhFeedDownMCpt(NULL),
  fhDirectMCptMax(NULL),
  fhDirectMCptMin(NULL),
  fhFeedDownMCptMax(NULL),
  fhFeedDownMCptMin(NULL),
  fhDirectEffpt(NULL),
  fhFeedDownEffpt(NULL),
  fhRECpt(NULL),
  fgRECSystematics(NULL),
  fNevts(1),
  fLuminosity(),
  fTrigEfficiency(),
  fGlobalEfficiencyUncertainties(),
  fTab(),
  fhFc(NULL),
  fhFcMax(NULL),
  fhFcMin(NULL),
  fhFcRcb(NULL),
  fgFcExtreme(NULL),
  fgFcConservative(NULL),
  fhYieldCorr(NULL),
  fhYieldCorrMax(NULL),
  fhYieldCorrMin(NULL),
  fhYieldCorrRcb(NULL),
  fgYieldCorr(NULL),
  fgYieldCorrExtreme(NULL),
  fgYieldCorrConservative(NULL),
  fhSigmaCorr(NULL),
  fhSigmaCorrMax(NULL),
  fhSigmaCorrMin(NULL),
  fhSigmaCorrDataSyst(NULL),
  fhSigmaCorrRcb(NULL),
  fgSigmaCorr(NULL),
  fgSigmaCorrExtreme(NULL),
  fgSigmaCorrConservative(NULL),
  fnSigma(NULL),
  fFeedDownOption(option),
  fAsymUncertainties(kTRUE),
  fPbPbElossHypothesis(kFALSE),
  fhStatUncEffcSigma(NULL),
  fhStatUncEffbSigma(NULL),
  fhStatUncEffcFD(NULL),
  fhStatUncEffbFD(NULL)
{
  //
  // Default constructor
  //

  fLuminosity[0]=1.;  fLuminosity[1]=0.;  
  fTrigEfficiency[0]=1.; fTrigEfficiency[1]=0.; 
  fGlobalEfficiencyUncertainties[0]=0.; fGlobalEfficiencyUncertainties[1]=0.;
  fTab[0]=1.;  fTab[1]=0.;

}

//_________________________________________________________________________________________________________
AliHFPtSpectrum::AliHFPtSpectrum(const AliHFPtSpectrum &rhs):
  TNamed(rhs),
  fhDirectMCpt(rhs.fhDirectMCpt),
  fhFeedDownMCpt(rhs.fhFeedDownMCpt),
  fhDirectMCptMax(rhs.fhDirectMCptMax),
  fhDirectMCptMin(rhs.fhDirectMCptMin),
  fhFeedDownMCptMax(rhs.fhFeedDownMCptMax),
  fhFeedDownMCptMin(rhs.fhFeedDownMCptMin),
  fhDirectEffpt(rhs.fhDirectEffpt),
  fhFeedDownEffpt(rhs.fhFeedDownEffpt),
  fhRECpt(rhs.fhRECpt),
  fgRECSystematics(rhs.fgRECSystematics),
  fNevts(rhs.fNevts),
  fLuminosity(),
  fTrigEfficiency(),
  fGlobalEfficiencyUncertainties(),
  fTab(),
  fhFc(rhs.fhFc),
  fhFcMax(rhs.fhFcMax),
  fhFcMin(rhs.fhFcMin),
  fhFcRcb(rhs.fhFcRcb),
  fgFcExtreme(rhs.fgFcExtreme),
  fgFcConservative(rhs.fgFcConservative),
  fhYieldCorr(rhs.fhYieldCorr),
  fhYieldCorrMax(rhs.fhYieldCorrMax),
  fhYieldCorrMin(rhs.fhYieldCorrMin),
  fhYieldCorrRcb(rhs.fhYieldCorrRcb),
  fgYieldCorr(rhs.fgYieldCorr),
  fgYieldCorrExtreme(rhs.fgYieldCorrExtreme),
  fgYieldCorrConservative(rhs.fgYieldCorrConservative),
  fhSigmaCorr(rhs.fhSigmaCorr),
  fhSigmaCorrMax(rhs.fhSigmaCorrMax),
  fhSigmaCorrMin(rhs.fhSigmaCorrMin),
  fhSigmaCorrDataSyst(rhs.fhSigmaCorrDataSyst),
  fhSigmaCorrRcb(rhs.fhSigmaCorrRcb),
  fgSigmaCorr(rhs.fgSigmaCorr),
  fgSigmaCorrExtreme(rhs.fgSigmaCorrExtreme),
  fgSigmaCorrConservative(rhs.fgSigmaCorrConservative),
  fnSigma(rhs.fnSigma),
  fFeedDownOption(rhs.fFeedDownOption),
  fAsymUncertainties(rhs.fAsymUncertainties),
  fPbPbElossHypothesis(rhs.fPbPbElossHypothesis),
  fhStatUncEffcSigma(NULL),
  fhStatUncEffbSigma(NULL),
  fhStatUncEffcFD(NULL),
  fhStatUncEffbFD(NULL)
{
  //
  // Copy constructor
  //

  for(Int_t i=0; i<2; i++){
    fLuminosity[i] = rhs.fLuminosity[i];
    fTrigEfficiency[i] = rhs.fTrigEfficiency[i];
    fGlobalEfficiencyUncertainties[i] = rhs.fGlobalEfficiencyUncertainties[i];
    fTab[i] = rhs.fTab[i];
  }

}

//_________________________________________________________________________________________________________
AliHFPtSpectrum &AliHFPtSpectrum::operator=(const AliHFPtSpectrum &source){
  //
  // Assignment operator
  //

  if (&source == this) return *this;
  
  fhDirectMCpt = source.fhDirectMCpt;
  fhFeedDownMCpt = source.fhFeedDownMCpt;
  fhDirectMCptMax = source.fhDirectMCptMax;
  fhDirectMCptMin = source.fhDirectMCptMin;
  fhFeedDownMCptMax = source.fhFeedDownMCptMax;
  fhFeedDownMCptMin = source.fhFeedDownMCptMin;
  fhDirectEffpt = source.fhDirectEffpt;
  fhFeedDownEffpt = source.fhFeedDownEffpt;
  fhRECpt = source.fhRECpt;
  fgRECSystematics = source.fgRECSystematics;
  fNevts = source.fNevts;
  fhFc = source.fhFc;
  fhFcMax = source.fhFcMax;
  fhFcMin = source.fhFcMin;
  fhFcRcb = source.fhFcRcb;
  fgFcExtreme = source.fgFcExtreme;
  fgFcConservative = source.fgFcConservative;
  fhYieldCorr = source.fhYieldCorr;
  fhYieldCorrMax = source.fhYieldCorrMax;
  fhYieldCorrMin = source.fhYieldCorrMin;
  fhYieldCorrRcb = source.fhYieldCorrRcb;
  fgYieldCorr = source.fgYieldCorr;
  fgYieldCorrExtreme = source.fgYieldCorrExtreme;
  fgYieldCorrConservative = source.fgYieldCorrConservative;
  fhSigmaCorr = source.fhSigmaCorr;
  fhSigmaCorrMax = source.fhSigmaCorrMax;
  fhSigmaCorrMin = source.fhSigmaCorrMin;
  fhSigmaCorrDataSyst = source.fhSigmaCorrDataSyst;
  fhSigmaCorrRcb = source.fhSigmaCorrRcb;
  fgSigmaCorr = source.fgSigmaCorr;
  fgSigmaCorrExtreme = source.fgSigmaCorrExtreme;
  fgSigmaCorrConservative = source.fgSigmaCorrConservative;
  fnSigma = source.fnSigma;
  fFeedDownOption = source.fFeedDownOption;
  fAsymUncertainties = source.fAsymUncertainties;
  fPbPbElossHypothesis = source.fPbPbElossHypothesis;
  
  for(Int_t i=0; i<2; i++){
    fLuminosity[i] = source.fLuminosity[i];
    fTrigEfficiency[i] = source.fTrigEfficiency[i];
    fGlobalEfficiencyUncertainties[i] = source.fGlobalEfficiencyUncertainties[i];
    fTab[i] = source.fTab[i];
  }

  return *this;
}

//_________________________________________________________________________________________________________
AliHFPtSpectrum::~AliHFPtSpectrum(){
  //
  // Destructor
  //
  if (fhDirectMCpt) delete fhDirectMCpt;    
  if (fhFeedDownMCpt) delete fhFeedDownMCpt;  
  if (fhDirectMCptMax) delete fhDirectMCptMax; 
  if (fhDirectMCptMin) delete fhDirectMCptMin; 
  if (fhFeedDownMCptMax) delete fhFeedDownMCptMax;
  if (fhFeedDownMCptMin) delete fhFeedDownMCptMin;
  if (fhDirectEffpt) delete fhDirectEffpt;    
  if (fhFeedDownEffpt) delete fhFeedDownEffpt;  
  if (fhRECpt) delete fhRECpt;    
  if (fgRECSystematics) delete fgRECSystematics;
  if (fhFc) delete fhFc;
  if (fhFcMax) delete fhFcMax;
  if (fhFcMin) delete fhFcMin;
  if (fhFcRcb) delete fhFcRcb;
  if (fgFcExtreme) delete fgFcExtreme;  
  if (fgFcConservative) delete fgFcConservative;
  if (fhYieldCorr) delete fhYieldCorr;                
  if (fhYieldCorrMax) delete fhYieldCorrMax;             
  if (fhYieldCorrMin) delete fhYieldCorrMin;    
  if (fhYieldCorrRcb) delete fhYieldCorrRcb;
  if (fgYieldCorr) delete fgYieldCorr;  
  if (fgYieldCorrExtreme) delete fgYieldCorrExtreme;
  if (fgYieldCorrConservative) delete fgYieldCorrConservative;
  if (fhSigmaCorr) delete fhSigmaCorr;                  
  if (fhSigmaCorrMax) delete fhSigmaCorrMax;               
  if (fhSigmaCorrMin) delete fhSigmaCorrMin; 
  if (fhSigmaCorrDataSyst) delete fhSigmaCorrDataSyst;
  if (fgSigmaCorr) delete fgSigmaCorr;    
  if (fgSigmaCorrExtreme) delete fgSigmaCorrExtreme;
  if (fgSigmaCorrConservative) delete fgSigmaCorrConservative;
  if (fnSigma) delete fnSigma;
}
  

//_________________________________________________________________________________________________________
TH1D * AliHFPtSpectrum::RebinTheoreticalSpectra(TH1D *hTheory, const char *name) {
  //
  // Function to rebin the theoretical spectrum 
  //  with respect to the real-data reconstructed spectrum binning 
  //
  
  if (!hTheory || !fhRECpt) {
    AliError("Feed-down or reconstructed spectra don't exist");
    return NULL;
  }

  //
  // Get the reconstructed spectra bins & limits
  Int_t nbins = fhRECpt->GetNbinsX();
  Int_t nbinsMC = hTheory->GetNbinsX();
  Double_t *limits = new Double_t[nbins+1];
  Double_t xlow=0., binwidth=0.;
  for (Int_t i=1; i<=nbins; i++) {
    binwidth = fhRECpt->GetBinWidth(i);
    xlow = fhRECpt->GetBinLowEdge(i);
    limits[i-1] = xlow;
  }
  limits[nbins] = xlow + binwidth;

  // Check that the reconstructed spectra binning 
  // is larger than the theoretical one
  Double_t thbinwidth = hTheory->GetBinWidth(1);
  for (Int_t i=1; i<=nbins; i++) {
    binwidth = fhRECpt->GetBinWidth(i);
    if ( thbinwidth > binwidth ) {
      AliInfo(" Beware it seems that the reconstructed spectra has a smaller binning than the theoretical predictions !! ");
    }
  }

  //
  // Define a new histogram with the real-data reconstructed spectrum binning 
  TH1D * hTheoryRebin = new TH1D(name," theoretical rebinned prediction",nbins,limits);

  Double_t sum[nbins], items[nbins];
  for (Int_t ibin=0; ibin<nbins; ibin++) {
    sum[ibin]=0.; items[ibin]=0.;
  }
  for (Int_t ibin=0; ibin<=nbinsMC; ibin++){
    
    for (Int_t ibinrec=0; ibinrec<nbins; ibinrec++){
      if (hTheory->GetBinCenter(ibin)>limits[ibinrec] && 
	  hTheory->GetBinCenter(ibin)<limits[ibinrec+1]){
	sum[ibinrec]+=hTheory->GetBinContent(ibin);
	items[ibinrec]+=1.;
      }
    }
    
  }

  // set the theoretical rebinned spectra to ( sum-bins / n-bins ) per new bin
  for (Int_t ibinrec=0; ibinrec<nbins; ibinrec++) {
    hTheoryRebin->SetBinContent(ibinrec+1,sum[ibinrec]/items[ibinrec]);
  }
  
  return (TH1D*)hTheoryRebin;
}

//_________________________________________________________________________________________________________
void AliHFPtSpectrum::SetMCptSpectra(TH1D *hDirect, TH1D *hFeedDown){
  //
  // Set the MonteCarlo or Theoretical spectra
  //  both for direct and feed-down contributions
  //
  
  if (!hDirect || !hFeedDown || !fhRECpt) {
    AliError("One or both (direct, feed-down) spectra or the reconstructed spectra don't exist");
    return;
  }

  Bool_t areconsistent = kTRUE;
  areconsistent = CheckHistosConsistency(hDirect,hFeedDown);
  if (!areconsistent) {
    AliInfo("Histograms are not consistent (bin width, bounds)"); 
    return;
  }

  //
  // Rebin the theoretical predictions to the reconstructed spectra binning 
  //
  fhDirectMCpt = RebinTheoreticalSpectra(hDirect,"fhDirectMCpt");
  fhDirectMCpt->SetNameTitle("fhDirectMCpt"," direct theoretical prediction");
  fhFeedDownMCpt = RebinTheoreticalSpectra(hFeedDown,"fhFeedDownMCpt");
  fhFeedDownMCpt->SetNameTitle("fhFeedDownMCpt"," feed-down theoretical prediction");

}

//_________________________________________________________________________________________________________
void AliHFPtSpectrum::SetFeedDownMCptSpectra(TH1D *hFeedDown){
  //
  // Set the MonteCarlo or Theoretical spectra
  //  for feed-down contribution
  //
  
  if (!hFeedDown || !fhRECpt) {
    AliError("Feed-down or reconstructed spectra don't exist");
    return;
  }

  //
  // Rebin the theoretical predictions to the reconstructed spectra binning 
  //
  fhFeedDownMCpt = RebinTheoreticalSpectra(hFeedDown,"fhFeedDownMCpt");
  fhFeedDownMCpt->SetNameTitle("fhFeedDownMCpt"," feed-down theoretical prediction");

}

//_________________________________________________________________________________________________________
void AliHFPtSpectrum::SetMCptDistributionsBounds(TH1D *hDirectMax, TH1D *hDirectMin, TH1D *hFeedDownMax, TH1D *hFeedDownMin){
  //
  // Set the maximum and minimum MonteCarlo or Theoretical spectra
  //  both for direct and feed-down contributions
  // used in case uncertainties are asymmetric and ca not be on the "basic histograms"
  //

  if (!hDirectMax || !hDirectMin || !hFeedDownMax|| !hFeedDownMin || !fhRECpt) {
    AliError("One or all of the max/min direct/feed-down or the reconstructed spectra don't exist");
    return;
  }

  Bool_t areconsistent = kTRUE; 
  areconsistent &= CheckHistosConsistency(hDirectMax,hDirectMin);
  areconsistent &= CheckHistosConsistency(hFeedDownMax,hFeedDownMin);
  areconsistent &= CheckHistosConsistency(hDirectMax,hFeedDownMax);
  if (!areconsistent) {
    AliInfo("Histograms are not consistent (bin width, bounds)"); 
    return;
  }


  //
  // Rebin the theoretical predictions to the reconstructed spectra binning 
  //
  fhDirectMCptMax = RebinTheoreticalSpectra(hDirectMax,"fhDirectMCptMax");
  fhDirectMCptMax->SetNameTitle("fhDirectMCptMax"," maximum direct theoretical prediction");
  fhDirectMCptMin = RebinTheoreticalSpectra(hDirectMin,"fhDirectMCptMin");
  fhDirectMCptMin->SetNameTitle("fhDirectMCptMin"," minimum direct theoretical prediction");
  fhFeedDownMCptMax = RebinTheoreticalSpectra(hFeedDownMax,"fhFeedDownMCptMax");
  fhFeedDownMCptMax->SetNameTitle("fhFeedDownMCptMax"," maximum feed-down theoretical prediction");
  fhFeedDownMCptMin = RebinTheoreticalSpectra(hFeedDownMin,"fhFeedDownMCptMin");
  fhFeedDownMCptMin->SetNameTitle("fhFeedDownMCptMin"," minimum feed-down theoretical prediction");

}

//_________________________________________________________________________________________________________
void AliHFPtSpectrum::SetFeedDownMCptDistributionsBounds(TH1D *hFeedDownMax, TH1D *hFeedDownMin){
  //
  // Set the maximum and minimum MonteCarlo or Theoretical spectra
  //   for feed-down contributions
  // used in case uncertainties are asymmetric and can not be on the "basic histogram"
  //

  if (!hFeedDownMax || !hFeedDownMin || !fhRECpt) {
    AliError("One or all of the max/min direct/feed-down spectra don't exist");
    return;
  }

  Bool_t areconsistent = kTRUE; 
  areconsistent &= CheckHistosConsistency(hFeedDownMax,hFeedDownMin);
  if (!areconsistent) {
    AliInfo("Histograms are not consistent (bin width, bounds)"); 
    return;
  }


  //
  // Rebin the theoretical predictions to the reconstructed spectra binning 
  //
  fhFeedDownMCptMax = RebinTheoreticalSpectra(hFeedDownMax,"fhFeedDownMCptMax");
  fhFeedDownMCptMax->SetNameTitle("fhFeedDownMCptMax"," maximum feed-down theoretical prediction");
  fhFeedDownMCptMin = RebinTheoreticalSpectra(hFeedDownMin,"fhFeedDownMCptMin");
  fhFeedDownMCptMin->SetNameTitle("fhFeedDownMCptMin"," minimum feed-down theoretical prediction");

}

//_________________________________________________________________________________________________________
void AliHFPtSpectrum::SetDirectAccEffCorrection(TH1D *hDirectEff){
  //
  // Set the Acceptance and Efficiency corrections 
  //   for the direct contribution
  //
  
  if (!hDirectEff) {
    AliError("The direct acceptance and efficiency corrections doesn't exist");
    return;
  }

  fhDirectEffpt = (TH1D*)hDirectEff->Clone();
  fhDirectEffpt->SetNameTitle("fhDirectEffpt"," direct acceptance x efficiency correction");
}

//_________________________________________________________________________________________________________
void AliHFPtSpectrum::SetAccEffCorrection(TH1D *hDirectEff, TH1D *hFeedDownEff){
  //
  // Set the Acceptance and Efficiency corrections 
  //  both for direct and feed-down contributions
  //
  
  if (!hDirectEff || !hFeedDownEff) {
    AliError("One or both (direct, feed-down) acceptance and efficiency corrections don't exist");
    return;
  }

  Bool_t areconsistent=kTRUE;
  areconsistent = CheckHistosConsistency(hDirectEff,hFeedDownEff);
  if (!areconsistent) {
    AliInfo("Histograms are not consistent (bin width, bounds)"); 
    return;
  }

  fhDirectEffpt = (TH1D*)hDirectEff->Clone();
  fhFeedDownEffpt = (TH1D*)hFeedDownEff->Clone();
  fhDirectEffpt->SetNameTitle("fhDirectEffpt"," direct acceptance x efficiency correction");
  fhFeedDownEffpt->SetNameTitle("fhFeedDownEffpt"," feed-down acceptance x efficiency correction");
}

//_________________________________________________________________________________________________________
void AliHFPtSpectrum::SetReconstructedSpectrum(TH1D *hRec) {
  //
  // Set the reconstructed spectrum
  //
  
  if (!hRec) {
    AliError("The reconstructed spectrum doesn't exist");
    return;
  }

  fhRECpt = (TH1D*)hRec->Clone();
  fhRECpt->SetNameTitle("fhRECpt"," reconstructed spectrum");
}

//_________________________________________________________________________________________________________
void AliHFPtSpectrum::SetReconstructedSpectrumSystematics(TGraphAsymmErrors *gRec) {
  //
  // Set the reconstructed spectrum (uncorrected yield) systematic uncertainties
  // 

  // Check the compatibility with the reconstructed spectrum
  Double_t gbinwidth = gRec->GetErrorXlow(1) + gRec->GetErrorXhigh(1) ;
  Double_t hbinwidth = fhRECpt->GetBinWidth(1);
  Double_t gxbincenter=0., gybincenter=0.;
  gRec->GetPoint(1,gxbincenter,gybincenter);
  Double_t hbincenter = fhRECpt->GetBinCenter(1);
  if ( (gbinwidth != hbinwidth) || (gxbincenter!=hbincenter) ) {
    AliError(" The reconstructed spectrum and its systematics don't seem compatible");
    return;
  }
  
  fgRECSystematics = gRec; 
}

//_________________________________________________________________________________________________________
void AliHFPtSpectrum::ComputeHFPtSpectrum(Double_t deltaY, Double_t branchingRatioC, Double_t branchingRatioBintoFinalDecay) {
  //
  // Main function to compute the corrected cross-section:
  // variables : analysed delta_y, BR for the final correction,
  //             BR b --> D --> decay (relative to the input theoretical prediction)
  //
  //   Sigma = ( 1. / (lumi * delta_y * BR_c * eff_trig * eff_c ) ) * spectra (corrected for feed-down)
  //
  // Uncertainties: (stat) delta_sigma = sigma * sqrt ( (delta_spectra/spectra)^2 )
  //  (syst but feed-down) delta_sigma = sigma * sqrt ( (delta_spectra_syst/spectra)^2 + (delta_lumi/lumi)^2 + (delta_eff_trig/eff_trig)^2 + (delta_eff/eff)^2 )
  //      (feed-down syst) delta_sigma = sigma * sqrt ( (delta_spectra_fd/spectra_fd)^2 )
  //
  //  In HIC the feed-down correction varies with an energy loss hypothesis:
  //      Raa(c-->D) / Raa(b-->D) for the "fc" method, Raa(b-->D) for the "Nb" method (see exact formulas in the functions)
  //

  //
  // First: Initialization
  //
  Bool_t areHistosOk = Initialize();
  if (!areHistosOk) {
    AliInfo(" Histos not properly initialized. Check : inconsistent binning ? missing histos ?");
    return;
  }

  //
  // Second: Correct for feed-down
  //
  if (fFeedDownOption==1) {
    // Compute the feed-down correction via fc-method
    CalculateFeedDownCorrectionFc(); 
    // Correct the yield for feed-down correction via fc-method
    CalculateFeedDownCorrectedSpectrumFc(); 
  }
  else if (fFeedDownOption==2) {
    // Correct the yield for feed-down correction via Nb-method
    CalculateFeedDownCorrectedSpectrumNb(deltaY,branchingRatioBintoFinalDecay); 
  }
  else if (fFeedDownOption==0) { 
    // If there is no need for feed-down correction,
    //    the "corrected" yield is equal to the raw yield
    fhYieldCorr = (TH1D*)fhRECpt->Clone();
    fhYieldCorr->SetNameTitle("fhYieldCorr","un-corrected yield");
    fhYieldCorrMax = (TH1D*)fhRECpt->Clone();
    fhYieldCorrMin = (TH1D*)fhRECpt->Clone();
    fhYieldCorrMax->SetNameTitle("fhYieldCorrMax","un-corrected yield");
    fhYieldCorrMin->SetNameTitle("fhYieldCorrMin","un-corrected yield");
    fAsymUncertainties=kFALSE;
  }
  else { 
    AliInfo(" Are you sure the feed-down correction option is right ?"); 
  }

  // Print out information
  printf("\n\n     Correcting the spectra with : \n   luminosity = %2.2e +- %2.2e, trigger efficiency = %2.2e +- %2.2e, \n    delta_y = %2.2f, BR_c = %2.2e, BR_b_decay = %2.2e \n    %2.2f percent uncertainty on the efficiencies, and %2.2f percent uncertainty on the b/c efficiencies ratio \n\n",fLuminosity[0],fLuminosity[1],fTrigEfficiency[0],fTrigEfficiency[1],deltaY,branchingRatioC,branchingRatioBintoFinalDecay,fGlobalEfficiencyUncertainties[0],fGlobalEfficiencyUncertainties[1]);
  if (fPbPbElossHypothesis)  printf("\n\n     The considered Tab is  %4.2e +- %2.2e \n\n",fTab[0],fTab[1]);

  //
  // Finally: Correct from yields to cross-section
  //
  Int_t nbins = fhRECpt->GetNbinsX();
  Double_t binwidth = fhRECpt->GetBinWidth(1);
  Double_t *limits = new Double_t[nbins+1]; 
  Double_t *binwidths = new Double_t[nbins]; 
  Double_t xlow=0.;
  for (Int_t i=1; i<=nbins; i++) {
    binwidth = fhRECpt->GetBinWidth(i);
    xlow = fhRECpt->GetBinLowEdge(i);
    limits[i-1] = xlow;
    binwidths[i-1] = binwidth;
  }
  limits[nbins] = xlow + binwidth;

  
  // declare the output histograms
  fhSigmaCorr = new TH1D("fhSigmaCorr","corrected sigma",nbins,limits);
  fhSigmaCorrMax = new TH1D("fhSigmaCorrMax","max corrected sigma",nbins,limits);
  fhSigmaCorrMin = new TH1D("fhSigmaCorrMin","min corrected sigma",nbins,limits);
  fhSigmaCorrDataSyst = new TH1D("fhSigmaCorrDataSyst","data syst uncertainties on the corrected sigma",nbins,limits);
  if (fPbPbElossHypothesis && fFeedDownOption==1) {
    fhSigmaCorrRcb = new TH2D("fhSigmaCorrRcb","corrected sigma vs Rcb Eloss hypothesis; p_{T} [GeV/c] ; Rcb Eloss hypothesis ; #sigma",nbins,limits,800,0.,4.);
    fnSigma = new TNtuple("fnSigma"," Sigma ntuple calculation","pt:Signal:Rcb:fc:Yield:Sigma");
  }
  if (fPbPbElossHypothesis && fFeedDownOption==2) {
    fhSigmaCorrRcb = new TH2D("fhSigmaCorrRcb","corrected sigma vs Rb Eloss hypothesis; p_{T} [GeV/c] ; Rb Eloss hypothesis ; #sigma",nbins,limits,800,0.,4.);
    fnSigma = new TNtuple("fnSigma"," Sigma ntuple calculation","pt:Signal:Rb:fc:Yield:Sigma");
  }
  // and the output TGraphAsymmErrors
  if (fAsymUncertainties){
    fgSigmaCorr = new TGraphAsymmErrors(nbins+1);
    fgSigmaCorrExtreme = new TGraphAsymmErrors(nbins+1);
    fgSigmaCorrConservative = new TGraphAsymmErrors(nbins+1);
  }
  fhStatUncEffcSigma = new TH1D("fhStatUncEffcSigma","direct charm stat unc on the cross section",nbins,limits);
  fhStatUncEffbSigma = new TH1D("fhStatUncEffbSigma","secondary charm stat unc on the cross section",nbins,limits);


  // protect against null denominator
  if (deltaY==0. || fLuminosity[0]==0. || fTrigEfficiency[0]==0. || branchingRatioC==0.) {
    AliError(" Hey you ! Why luminosity or trigger-efficiency or the c-BR or delta_y are set to zero ?! ");
    delete [] limits;
    delete [] binwidths;
    return ;
  }

  Double_t value=0, errValue=0, errvalueMax=0., errvalueMin=0.;
  Double_t errvalueExtremeMax=0., errvalueExtremeMin=0.;
  Double_t errvalueConservativeMax=0., errvalueConservativeMin=0.;
  Double_t errvalueStatUncEffc=0., errvalueStatUncEffb=0.;
  for(Int_t ibin=1; ibin<=nbins; ibin++){
    
    // Variables initialization
    value=0.; errValue=0.; errvalueMax=0.; errvalueMin=0.;
    errvalueExtremeMax=0.; errvalueExtremeMin=0.;
    errvalueConservativeMax=0.; errvalueConservativeMin=0.;
    errvalueStatUncEffc=0.; errvalueStatUncEffb=0.;

    // Sigma calculation
    //   Sigma = ( 1. / (lumi * delta_y * BR_c * eff_trig * eff_c ) ) * spectra (corrected for feed-down)
    value = (fhDirectEffpt->GetBinContent(ibin) && fhDirectEffpt->GetBinContent(ibin)!=0. && fhRECpt->GetBinContent(ibin)>0.) ? 
      ( fhYieldCorr->GetBinContent(ibin) / ( deltaY * branchingRatioC * fLuminosity[0] * fTrigEfficiency[0] * fhDirectEffpt->GetBinContent(ibin) ) )
      : 0. ;

    // Sigma statistical uncertainty:
    //   delta_sigma = sigma * sqrt ( (delta_spectra/spectra)^2 )
    errValue = (value!=0.) ?  value * (fhYieldCorr->GetBinError(ibin)/fhYieldCorr->GetBinContent(ibin))  : 0. ;

    //    cout<< " x "<< fhRECpt->GetBinCenter(ibin) << " sigma " << value << " +- "<< errValue << " (stat)"<<endl;

    //
    // Sigma systematic uncertainties
    //
    if (fAsymUncertainties && value>0.) {
      
      //  (syst but feed-down) delta_sigma = sigma * sqrt ( (delta_spectra_syst/spectra)^2 +
      //                                     (delta_lumi/lumi)^2 + (delta_eff_trig/eff_trig)^2 + (delta_eff/eff)^2  + (global_eff)^2 )
      errvalueMax = value * TMath::Sqrt( (fgYieldCorr->GetErrorYhigh(ibin)/fhYieldCorr->GetBinContent(ibin))*(fgYieldCorr->GetErrorYhigh(ibin)/fhYieldCorr->GetBinContent(ibin)) + 
					 (fLuminosity[1]/fLuminosity[0])*(fLuminosity[1]/fLuminosity[0]) + 
					 (fTrigEfficiency[1]/fTrigEfficiency[0])*(fTrigEfficiency[1]/fTrigEfficiency[0])  +
					 (fhDirectEffpt->GetBinError(ibin)/fhDirectEffpt->GetBinContent(ibin))*(fhDirectEffpt->GetBinError(ibin)/fhDirectEffpt->GetBinContent(ibin)) +
					 fGlobalEfficiencyUncertainties[0]*fGlobalEfficiencyUncertainties[0] );
      errvalueMin = value * TMath::Sqrt( (fgYieldCorr->GetErrorYlow(ibin)/fhYieldCorr->GetBinContent(ibin))*(fgYieldCorr->GetErrorYlow(ibin)/fhYieldCorr->GetBinContent(ibin)) + 
					 (fLuminosity[1]/fLuminosity[0])*(fLuminosity[1]/fLuminosity[0]) + 
					 (fTrigEfficiency[1]/fTrigEfficiency[0])*(fTrigEfficiency[1]/fTrigEfficiency[0])  +
					 (fhDirectEffpt->GetBinError(ibin)/fhDirectEffpt->GetBinContent(ibin))*(fhDirectEffpt->GetBinError(ibin)/fhDirectEffpt->GetBinContent(ibin))  +
					 fGlobalEfficiencyUncertainties[0]*fGlobalEfficiencyUncertainties[0] );

      // Uncertainties from feed-down
      //      (feed-down syst) delta_sigma = sigma * sqrt ( (delta_spectra_fd/spectra_fd)^2 )
      //   extreme case
      errvalueExtremeMax = value * (fgYieldCorrExtreme->GetErrorYhigh(ibin)/fhYieldCorr->GetBinContent(ibin));
      errvalueExtremeMin =  value * (fgYieldCorrExtreme->GetErrorYlow(ibin)/fhYieldCorr->GetBinContent(ibin));
      //
      //   conservative case
      errvalueConservativeMax = value * (fgYieldCorrConservative->GetErrorYhigh(ibin)/fhYieldCorr->GetBinContent(ibin));
      errvalueConservativeMin =  value * (fgYieldCorrConservative->GetErrorYlow(ibin)/fhYieldCorr->GetBinContent(ibin));


      // stat unc of the efficiencies, separately
      errvalueStatUncEffc = value * (fhDirectEffpt->GetBinError(ibin)/fhDirectEffpt->GetBinContent(ibin)) ;
      errvalueStatUncEffb = 0.;

    }
    else {
      // protect against null denominator
      errvalueMax = (value!=0.) ?
	value * TMath::Sqrt( (fLuminosity[1]/fLuminosity[0])*(fLuminosity[1]/fLuminosity[0]) + 
			     (fTrigEfficiency[1]/fTrigEfficiency[0])*(fTrigEfficiency[1]/fTrigEfficiency[0])  +
			     (fhDirectEffpt->GetBinError(ibin)/fhDirectEffpt->GetBinContent(ibin))*(fhDirectEffpt->GetBinError(ibin)/fhDirectEffpt->GetBinContent(ibin))  +
			     fGlobalEfficiencyUncertainties[0]*fGlobalEfficiencyUncertainties[0] )
	: 0. ;
      errvalueMin = errvalueMax;
    }
    
    //
    // Fill the histograms
    //
    fhSigmaCorr->SetBinContent(ibin,value);
    fhSigmaCorr->SetBinError(ibin,errValue);
    //
    // Fill the histos and ntuple vs the Eloss hypothesis
    //
    if (fPbPbElossHypothesis) {
      // Loop over the Eloss hypothesis
      for (Float_t rval=0.0025; rval<4.0; rval+=0.005) {
	Int_t rbin = FindTH2YBin(fhYieldCorrRcb,rval);
	Double_t yieldRcbvalue = fhYieldCorrRcb->GetBinContent(ibin,rbin);
	// Sigma calculation
	//   Sigma = ( 1. / (lumi * delta_y * BR_c * eff_trig * eff_c ) ) * spectra (corrected for feed-down)
	Double_t sigmaRcbvalue = (fhDirectEffpt->GetBinContent(ibin) && fhDirectEffpt->GetBinContent(ibin)!=0.) ? 
	  ( yieldRcbvalue / ( deltaY * branchingRatioC * fLuminosity[0] * fTrigEfficiency[0] * fhDirectEffpt->GetBinContent(ibin) ) )
	  : 0. ;
	fhSigmaCorrRcb->Fill( fhSigmaCorr->GetBinCenter(ibin) , rval, sigmaRcbvalue );
	// 	if(ibin==3) 
	// 	  cout << " pt "<< fhRECpt->GetBinCenter(ibin) <<" bin "<< ibin<<" rval="<<rval<<", rbin="<<rbin<<" fc-value="<< fhFcRcb->GetBinContent(ibin,rbin) <<", yield-fcRbvalue="<<yieldRcbvalue<<", sigma-fcRbvalue="<<sigmaRcbvalue<<endl;
	fnSigma->Fill(fhRECpt->GetBinCenter(ibin), fhRECpt->GetBinContent(ibin),
		      rval, fhFcRcb->GetBinContent(ibin,rbin),
		      yieldRcbvalue, sigmaRcbvalue );
      }
    }
    //
    // Fill the TGraphAsymmErrors
    if (fAsymUncertainties) {
      Double_t x = fhYieldCorr->GetBinCenter(ibin);
      fgSigmaCorr->SetPoint(ibin,x,value); // i,x,y
      fgSigmaCorr->SetPointError(ibin,(binwidths[ibin-1]/2.),(binwidths[ibin-1]/2.),errvalueMin,errvalueMax); // i,xl,xh,yl,yh
      fhSigmaCorrMax->SetBinContent(ibin,value+errvalueMax);
      fhSigmaCorrMin->SetBinContent(ibin,value-errvalueMin);
      fgSigmaCorrExtreme->SetPoint(ibin,x,value); // i,x,y
      fgSigmaCorrExtreme->SetPointError(ibin,(binwidths[ibin-1]/2.),(binwidths[ibin-1]/2.),errvalueExtremeMin,errvalueExtremeMax); // i,xl,xh,yl,yh
      fgSigmaCorrConservative->SetPoint(ibin,x,value); // i,x,y
      fgSigmaCorrConservative->SetPointError(ibin,(binwidths[ibin-1]/2.),(binwidths[ibin-1]/2.),errvalueConservativeMin,errvalueConservativeMax); // i,xl,xh,yl,yh

      fhStatUncEffcSigma->SetBinContent(ibin,0.); 
      if(value>0.) fhStatUncEffcSigma->SetBinError(ibin,((errvalueStatUncEffc/value)*100.));
      fhStatUncEffbSigma->SetBinContent(ibin,0.); fhStatUncEffbSigma->SetBinError(ibin,0.);
      //      cout << " pt "<< fhRECpt->GetBinCenter(ibin) <<" bin "<< ibin<<" stat-unc-c-sigma "<< errvalueStatUncEffc/value << endl;
    }
    
  }
  delete [] binwidths;
  delete [] limits;

}

//_________________________________________________________________________________________________________
TH1D * AliHFPtSpectrum::EstimateEfficiencyRecoBin(TH1D *hSimu, TH1D *hReco, const char *name) {
  //
  // Function that computes the acceptance and efficiency correction
  //  based on the simulated and reconstructed spectra
  //  and using the reconstructed spectra bin width
  //
  //  eff = reco/sim ; err_eff = sqrt( eff*(1-eff) )/ sqrt( sim )
  // 

  if(!fhRECpt){
    AliInfo("Hey, the reconstructed histogram was not set yet !"); 
    return NULL;
  }

  Int_t nbins = fhRECpt->GetNbinsX();
  Double_t *limits = new Double_t[nbins+1];
  Double_t xlow=0.,binwidth=0.;
  for (Int_t i=1; i<=nbins; i++) {
    binwidth = fhRECpt->GetBinWidth(i);
    xlow = fhRECpt->GetBinLowEdge(i);
    limits[i-1] = xlow;
  }
  limits[nbins] = xlow + binwidth;

  TH1D * hEfficiency = new TH1D(name," acceptance #times efficiency",nbins,limits);
  
  Double_t *sumSimu=new Double_t[nbins];
  Double_t *sumReco=new Double_t[nbins];
  for (Int_t ibin=0; ibin<nbins; ibin++){
    sumSimu[ibin]=0.;  sumReco[ibin]=0.;
  }
  for (Int_t ibin=0; ibin<=hSimu->GetNbinsX(); ibin++){

    for (Int_t ibinrec=0; ibinrec<nbins; ibinrec++){
      if ( hSimu->GetBinCenter(ibin)>limits[ibinrec] && 
	   hSimu->GetBinCenter(ibin)<limits[ibinrec+1] ) {
	sumSimu[ibinrec]+=hSimu->GetBinContent(ibin);
      }
      if ( hReco->GetBinCenter(ibin)>limits[ibinrec] && 
	   hReco->GetBinCenter(ibin)<limits[ibinrec+1] ) {
	sumReco[ibinrec]+=hReco->GetBinContent(ibin);
      }
    }
    
  }


  // the efficiency is computed as reco/sim (in each bin)
  //  its uncertainty is err_eff = sqrt( eff*(1-eff) )/ sqrt( sim )
  Double_t eff=0., erreff=0.;
  for (Int_t ibinrec=0; ibinrec<nbins; ibinrec++) {
    if (sumSimu[ibinrec]!= 0. && sumReco[ibinrec]!=0.) {
      eff = sumReco[ibinrec] / sumSimu[ibinrec] ;
      // protection in case eff > 1.0
      // test calculation (make the argument of the sqrt positive)
      erreff = TMath::Sqrt( eff * TMath::Abs(1.0 - eff) ) / TMath::Sqrt( sumSimu[ibinrec] );
    }
    else { eff=0.0; erreff=0.; }
    hEfficiency->SetBinContent(ibinrec+1,eff);
    hEfficiency->SetBinError(ibinrec+1,erreff);
  }

  delete [] sumSimu;
  delete [] sumReco;

  return (TH1D*)hEfficiency;
}

//_________________________________________________________________________________________________________
void AliHFPtSpectrum::EstimateAndSetDirectEfficiencyRecoBin(TH1D *hSimu, TH1D *hReco) {
  //
  // Function that computes the Direct  acceptance and efficiency correction
  //  based on the simulated and reconstructed spectra
  //  and using the reconstructed spectra bin width
  //
  //  eff = reco/sim ; err_eff = sqrt( eff*(1-eff) )/ sqrt( sim )
  // 

  if(!fhRECpt || !hSimu || !hReco){
    AliError("Hey, the reconstructed histogram was not set yet !"); 
    return;
  }

  fhDirectEffpt = EstimateEfficiencyRecoBin(hSimu,hReco,"fhDirectEffpt");
  fhDirectEffpt->SetNameTitle("fhDirectEffpt"," direct acceptance #times efficiency");

}

//_________________________________________________________________________________________________________
void AliHFPtSpectrum::EstimateAndSetFeedDownEfficiencyRecoBin(TH1D *hSimu, TH1D *hReco) {
  //
  // Function that computes the Feed-Down acceptance and efficiency correction
  //  based on the simulated and reconstructed spectra
  //  and using the reconstructed spectra bin width
  //
  //  eff = reco/sim ; err_eff = sqrt( eff*(1-eff) )/ sqrt( sim )
  // 
  
  if(!fhRECpt || !hSimu || !hReco){
    AliError("Hey, the reconstructed histogram was not set yet !"); 
    return;
  }

  fhFeedDownEffpt = EstimateEfficiencyRecoBin(hSimu,hReco,"fhFeedDownEffpt");
  fhFeedDownEffpt->SetNameTitle("fhFeedDownEffpt"," feed-down acceptance #times efficiency");

}

//_________________________________________________________________________________________________________
Bool_t AliHFPtSpectrum::Initialize(){
  //
  // Initialization of the variables (histograms)
  //

  if (fFeedDownOption==0) { 
    AliInfo("Getting ready for the corrections without feed-down consideration");
  } else if (fFeedDownOption==1) { 
    AliInfo("Getting ready for the fc feed-down correction calculation");
  } else if (fFeedDownOption==2) {
    AliInfo("Getting ready for the Nb feed-down correction calculation");
  } else { AliError("The calculation option must be <=2"); return kFALSE; }

  // Start checking the input histograms consistency
  Bool_t areconsistent=kTRUE;

  // General checks 
  if (!fhDirectEffpt || !fhRECpt) {
    AliError(" Reconstructed spectra and/or the Nc efficiency distributions are not defined");
    return kFALSE;
  }
  areconsistent &= CheckHistosConsistency(fhRECpt,fhDirectEffpt);
  if (!areconsistent) {
    AliInfo("Histograms required for Nb correction are not consistent (bin width, bounds)"); 
    return kFALSE;
  }
  if (fFeedDownOption==0) return kTRUE;

  //
  // Common checks for options 1 (fc) & 2(Nb)
  if (!fhFeedDownMCpt || !fhFeedDownEffpt) {
    AliError(" Theoretical Nb and/or the Nb efficiency distributions are not defined");
    return kFALSE;
  }
  areconsistent &= CheckHistosConsistency(fhRECpt,fhFeedDownMCpt);
  areconsistent &= CheckHistosConsistency(fhFeedDownMCpt,fhFeedDownEffpt);
  if (fAsymUncertainties) {
    if (!fhFeedDownMCptMax || !fhFeedDownMCptMin) {
      AliError(" Max/Min theoretical Nb distributions are not defined");
      return kFALSE;
    }
    areconsistent &= CheckHistosConsistency(fhFeedDownMCpt,fhFeedDownMCptMax);
  }
  if (!areconsistent) {
    AliInfo("Histograms required for Nb correction are not consistent (bin width, bounds)"); 
    return kFALSE;
  }
  if (fFeedDownOption>1) return kTRUE;  

  //
  // Now checks for option 1 (fc correction) 
  if (!fhDirectMCpt) {
    AliError("Theoretical Nc distributions is not defined");
    return kFALSE;
  }
  areconsistent &= CheckHistosConsistency(fhDirectMCpt,fhFeedDownMCpt);
  areconsistent &= CheckHistosConsistency(fhDirectMCpt,fhDirectEffpt);
  if (fAsymUncertainties) {
    if (!fhDirectMCptMax || !fhDirectMCptMin) {
      AliError(" Max/Min theoretical Nc distributions are not defined");
      return kFALSE;
    }
    areconsistent &= CheckHistosConsistency(fhDirectMCpt,fhDirectMCptMax);
  }
  if (!areconsistent) {
    AliInfo("Histograms required for fc correction are not consistent (bin width, bounds)"); 
    return kFALSE;
  }

  return kTRUE;
}

//_________________________________________________________________________________________________________
Bool_t AliHFPtSpectrum::CheckHistosConsistency(TH1D *h1, TH1D *h2){
  //
  // Check the histograms consistency (bins, limits)
  //

  if (!h1 || !h2) {
    AliError("One or both histograms don't exist");
    return kFALSE;
  }

  Double_t binwidth1 = h1->GetBinWidth(1);
  Double_t binwidth2 = h2->GetBinWidth(1);
  Double_t min1 = h1->GetBinCenter(1) - (binwidth1/2.) ; 
//   Double_t max1 = h1->GetBinCenter(nbins1) + (binwidth1/2.) ;
  Double_t min2 = h2->GetBinCenter(1) - (binwidth2/2.) ;
//   Double_t max2 = h2->GetBinCenter(nbins2) + (binwidth2/2.) ;

  if (binwidth1!=binwidth2) {
    AliInfo(" histograms with different bin width");
    return kFALSE;
  }
  if (min1!=min2) {
    AliInfo(" histograms with different minimum");
    return kFALSE;
  }
//   if (max1!=max2) {
//     AliInfo(" histograms with different maximum");
//     return kFALSE;
//   }

  return kTRUE;
}

//_________________________________________________________________________________________________________
void AliHFPtSpectrum::CalculateFeedDownCorrectionFc(){ 
  //
  // Compute fc factor and its uncertainties bin by bin
  //   fc = 1 / ( 1 + (eff_b/eff_c)*(N_b/N_c) ) 
  //
  // uncertainties: (conservative) combine the upper/lower N_b & N_c predictions together
  //                (extreme) combine the upper N_b predictions with the lower N_c predictions & viceversa
  //                systematic uncertainty on the acceptance x efficiency b/c ratio are included 
  //
  //  In addition, in HIC the feed-down correction varies with an energy loss hypothesis: Raa(c-->D) / Raa(b-->D) = Rcb
  //	       fc (Rcb) = ( 1. / ( 1 + (eff_b/eff_c)*(N_b/N_c)* (1/Rcb) ) );
  //
  AliInfo("Calculating the feed-down correction factor (fc method)");
  
  // define the variables
  Int_t nbins = fhRECpt->GetNbinsX();
  Double_t binwidth = fhRECpt->GetBinWidth(1);
  Double_t *limits = new Double_t[nbins+1];
  Double_t *binwidths = new Double_t[nbins];
  Double_t xlow=0.;
  for (Int_t i=1; i<=nbins; i++) {
    binwidth = fhRECpt->GetBinWidth(i);
    xlow = fhRECpt->GetBinLowEdge(i);
    limits[i-1] = xlow;
    binwidths[i-1] = binwidth;
  }
  limits[nbins] = xlow + binwidth;

  Double_t correction=1.;
  Double_t theoryRatio=1.;
  Double_t effRatio=1.; 
  Double_t correctionExtremeA=1., correctionExtremeB=1.;
  Double_t theoryRatioExtremeA=1., theoryRatioExtremeB=1.;
  Double_t correctionConservativeA=1., correctionConservativeB=1.;
  Double_t theoryRatioConservativeA=1., theoryRatioConservativeB=1.;
  Double_t correctionUnc=0.;
  Double_t correctionExtremeAUnc=0., correctionExtremeBUnc=0.;
  Double_t correctionConservativeAUnc=0., correctionConservativeBUnc=0.;

  // declare the output histograms
  fhFc = new TH1D("fhFc","fc correction factor",nbins,limits);
  fhFcMax = new TH1D("fhFcMax","max fc correction factor",nbins,limits);
  fhFcMin = new TH1D("fhFcMin","min fc correction factor",nbins,limits);
  if(fPbPbElossHypothesis) fhFcRcb = new TH2D("fhFcRcb","fc correction factor vs Rcb Eloss hypothesis; p_{T} [GeV/c] ; Rcb Eloss hypothesis ; fc correction",nbins,limits,800,0.,4.);
  // two local control histograms
  TH1D *hTheoryRatio = new TH1D("hTheoryRatio","Theoretical B-->D over c-->D (feed-down/direct) ratio",nbins,limits);
  TH1D *hEffRatio = new TH1D("hEffRatio","Efficiency B-->D over c-->D (feed-down/direct) ratio",nbins,limits);
  // and the output TGraphAsymmErrors
  if (fAsymUncertainties) {
    fgFcExtreme = new TGraphAsymmErrors(nbins+1);
    fgFcExtreme->SetNameTitle("fgFcExtreme","fgFcExtreme");
    fgFcConservative = new TGraphAsymmErrors(nbins+1);
    fgFcConservative->SetNameTitle("fgFcConservative","fgFcConservative");
  }

  fhStatUncEffcFD = new TH1D("fhStatUncEffcFD","direct charm stat unc on the feed-down correction",nbins,limits);
  fhStatUncEffbFD = new TH1D("fhStatUncEffbFD","secondary charm stat unc on the feed-down correction",nbins,limits);
  Double_t correctionConservativeAUncStatEffc=0., correctionConservativeBUncStatEffc=0.;
  Double_t correctionConservativeAUncStatEffb=0., correctionConservativeBUncStatEffb=0.;

  //
  // Compute fc
  //
  for (Int_t ibin=1; ibin<=nbins; ibin++) {

    //  theory_ratio = (N_b/N_c) 
    theoryRatio = (fhDirectMCpt->GetBinContent(ibin)>0. && fhFeedDownMCpt->GetBinContent(ibin)>0.) ? 
      fhFeedDownMCpt->GetBinContent(ibin) / fhDirectMCpt->GetBinContent(ibin) : 1.0 ;

    //
    // Calculate the uncertainty [ considering only the theoretical uncertainties on Nb & Nc for now !!! ]
    //
    // extreme A = direct-max, feed-down-min
    theoryRatioExtremeA = (fhDirectMCptMax->GetBinContent(ibin)>0. && fhFeedDownMCptMin->GetBinContent(ibin)>0.) ? 
      fhFeedDownMCptMin->GetBinContent(ibin) / fhDirectMCptMax->GetBinContent(ibin) : 1.0 ;
    // extreme B = direct-min, feed-down-max
    theoryRatioExtremeB = (fhDirectMCptMin->GetBinContent(ibin)>0. && fhDirectMCptMax->GetBinContent(ibin)>0.) ? 
      fhFeedDownMCptMax->GetBinContent(ibin) / fhDirectMCptMin->GetBinContent(ibin) : 1.0 ;
    // conservative A = direct-max, feed-down-max
    theoryRatioConservativeA = (fhDirectMCptMax->GetBinContent(ibin)>0. && fhFeedDownMCptMin->GetBinContent(ibin)>0.) ? 
      fhFeedDownMCptMax->GetBinContent(ibin) / fhDirectMCptMax->GetBinContent(ibin) : 1.0 ;
    // conservative B = direct-min, feed-down-min
    theoryRatioConservativeB = (fhDirectMCptMin->GetBinContent(ibin)>0. && fhDirectMCptMax->GetBinContent(ibin)>0.) ? 
      fhFeedDownMCptMin->GetBinContent(ibin) / fhDirectMCptMin->GetBinContent(ibin) : 1.0 ;

    //  eff_ratio = (eff_b/eff_c)
    effRatio = (fhDirectEffpt->GetBinContent(ibin) && fhDirectEffpt->GetBinContent(ibin)!=0.) ? 
      fhFeedDownEffpt->GetBinContent(ibin) / fhDirectEffpt->GetBinContent(ibin) : 1.0 ;

    //   fc = 1 / ( 1 + (eff_b/eff_c)*(N_b/N_c) ) 
    if( TMath::Abs(effRatio - 1.0)<0.01 || TMath::Abs(theoryRatio - 1.0)<0.01 ) {
      correction = 1.0;
      correctionExtremeA = 1.0;
      correctionExtremeB = 1.0;
      correctionConservativeA = 1.0; 
      correctionConservativeB = 1.0;
    }
    else {
      correction = ( 1. / ( 1 + ( effRatio * theoryRatio ) ) );
      correctionExtremeA = ( 1. / ( 1 + ( effRatio * theoryRatioExtremeA ) ) );
      correctionExtremeB = ( 1. / ( 1 + ( effRatio * theoryRatioExtremeB ) ) );
      correctionConservativeA = ( 1. / ( 1 + ( effRatio * theoryRatioConservativeA ) ) );
      correctionConservativeB = ( 1. / ( 1 + ( effRatio * theoryRatioConservativeB ) ) );
    }


    // fc uncertainty from (eff_b/eff_c) = fc^2 * (N_b/N_c) * delta(eff_b/eff_c)
    //  delta(eff_b/eff_c) is a percentage = effRatio * sqrt( fGlobalEfficiencyUncertainties[1]^2 + unc_eff_c ^2 + unc_eff_b ^2 ) 
    correctionUnc = correction*correction * theoryRatio * effRatio *
      TMath::Sqrt( fGlobalEfficiencyUncertainties[1]*fGlobalEfficiencyUncertainties[1] + 
		   (fhFeedDownEffpt->GetBinError(ibin)/fhFeedDownEffpt->GetBinContent(ibin))*(fhFeedDownEffpt->GetBinError(ibin)/fhFeedDownEffpt->GetBinContent(ibin)) +
		   (fhDirectEffpt->GetBinError(ibin)/fhDirectEffpt->GetBinContent(ibin))*(fhDirectEffpt->GetBinError(ibin)/fhDirectEffpt->GetBinContent(ibin)) 
		   );
    correctionExtremeAUnc = correctionExtremeA*correctionExtremeA * theoryRatioExtremeA  * effRatio *
      TMath::Sqrt( fGlobalEfficiencyUncertainties[1]*fGlobalEfficiencyUncertainties[1] + 
		   (fhFeedDownEffpt->GetBinError(ibin)/fhFeedDownEffpt->GetBinContent(ibin))*(fhFeedDownEffpt->GetBinError(ibin)/fhFeedDownEffpt->GetBinContent(ibin)) +
		   (fhDirectEffpt->GetBinError(ibin)/fhDirectEffpt->GetBinContent(ibin))*(fhDirectEffpt->GetBinError(ibin)/fhDirectEffpt->GetBinContent(ibin)) 
		   );
    correctionExtremeBUnc = correctionExtremeB*correctionExtremeB * theoryRatioExtremeB  * effRatio *
      TMath::Sqrt( fGlobalEfficiencyUncertainties[1]*fGlobalEfficiencyUncertainties[1] + 
		   (fhFeedDownEffpt->GetBinError(ibin)/fhFeedDownEffpt->GetBinContent(ibin))*(fhFeedDownEffpt->GetBinError(ibin)/fhFeedDownEffpt->GetBinContent(ibin)) +
		   (fhDirectEffpt->GetBinError(ibin)/fhDirectEffpt->GetBinContent(ibin))*(fhDirectEffpt->GetBinError(ibin)/fhDirectEffpt->GetBinContent(ibin)) 
		   );
    correctionConservativeAUnc = correctionConservativeA*correctionConservativeA * theoryRatioConservativeA  *effRatio *
      TMath::Sqrt( fGlobalEfficiencyUncertainties[1]*fGlobalEfficiencyUncertainties[1] + 
		   (fhFeedDownEffpt->GetBinError(ibin)/fhFeedDownEffpt->GetBinContent(ibin))*(fhFeedDownEffpt->GetBinError(ibin)/fhFeedDownEffpt->GetBinContent(ibin)) +
		   (fhDirectEffpt->GetBinError(ibin)/fhDirectEffpt->GetBinContent(ibin))*(fhDirectEffpt->GetBinError(ibin)/fhDirectEffpt->GetBinContent(ibin)) 
		   );
    //
    correctionConservativeAUncStatEffc = correctionConservativeA*correctionConservativeA * theoryRatioConservativeA  *effRatio * 
      (fhDirectEffpt->GetBinError(ibin)/fhDirectEffpt->GetBinContent(ibin));
    correctionConservativeAUncStatEffb = correctionConservativeA*correctionConservativeA * theoryRatioConservativeA  *effRatio * 
      (fhFeedDownEffpt->GetBinError(ibin)/fhFeedDownEffpt->GetBinContent(ibin));

    correctionConservativeBUnc = correctionConservativeB*correctionConservativeB * theoryRatioConservativeB  *effRatio *
      TMath::Sqrt( fGlobalEfficiencyUncertainties[1]*fGlobalEfficiencyUncertainties[1] + 
		   (fhFeedDownEffpt->GetBinError(ibin)/fhFeedDownEffpt->GetBinContent(ibin))*(fhFeedDownEffpt->GetBinError(ibin)/fhFeedDownEffpt->GetBinContent(ibin)) +
		   (fhDirectEffpt->GetBinError(ibin)/fhDirectEffpt->GetBinContent(ibin))*(fhDirectEffpt->GetBinError(ibin)/fhDirectEffpt->GetBinContent(ibin)) 
		   );
    correctionConservativeBUncStatEffb = correctionConservativeB*correctionConservativeB * theoryRatioConservativeB  *effRatio * 
      (fhFeedDownEffpt->GetBinError(ibin)/fhFeedDownEffpt->GetBinContent(ibin));
    correctionConservativeBUncStatEffc = correctionConservativeB*correctionConservativeB * theoryRatioConservativeB  *effRatio * 
      (fhDirectEffpt->GetBinError(ibin)/fhDirectEffpt->GetBinContent(ibin));


    // Fill in the histograms
    hTheoryRatio->SetBinContent(ibin,theoryRatio);
    hEffRatio->SetBinContent(ibin,effRatio);
    fhFc->SetBinContent(ibin,correction);
    //
    // Estimate how the result varies vs charm/beauty Eloss hypothesis
    //
    if ( TMath::Abs(correction-1.0)>0.01 && fPbPbElossHypothesis){
      // Loop over the Eloss hypothesis
      //      Int_t rbin=0;
      for (Float_t rval=0.0025; rval<4.0; rval+=0.005){
	Double_t correctionRcb = ( 1. / ( 1 + ( effRatio * theoryRatio * (1/rval) ) ) );
	fhFcRcb->Fill( fhFc->GetBinCenter(ibin) , rval, correctionRcb );
	// 	if(ibin==3){
	// 	  cout << " pt "<< fhFc->GetBinCenter(ibin) <<" bin "<< ibin<<" rval="<<rval<<", rbin="<<rbin<<", fc-Rcb-value="<<correctionRcb<<endl;
	// 	  rbin++;
	// 	}
      }
    }
    //
    // Fill the rest of (asymmetric) histograms
    //
    if (fAsymUncertainties) {
      Double_t x = fhDirectMCpt->GetBinCenter(ibin);
      Double_t val[4] = { correctionExtremeA + correctionExtremeAUnc, correctionExtremeA - correctionExtremeAUnc, 
			  correctionExtremeB + correctionExtremeBUnc, correctionExtremeB - correctionExtremeBUnc };
      Double_t uncExtremeMin = correction - TMath::MinElement(4,val);
      Double_t uncExtremeMax = TMath::MaxElement(4,val) - correction;
      fgFcExtreme->SetPoint(ibin,x,correction); // i,x,y
      fgFcExtreme->SetPointError(ibin,(binwidths[ibin-1]/2.),(binwidths[ibin-1]/2.),uncExtremeMin,uncExtremeMax); // i,xl,xh,yl,yh
      fhFcMax->SetBinContent(ibin,correction+uncExtremeMax);
      fhFcMin->SetBinContent(ibin,correction-uncExtremeMin);
      Double_t consval[4] = { correctionConservativeA - correctionConservativeAUnc, correctionConservativeA + correctionConservativeAUnc,
			      correctionConservativeB - correctionConservativeBUnc, correctionConservativeB + correctionConservativeBUnc};
      Double_t uncConservativeMin = correction - TMath::MinElement(4,consval);
      Double_t uncConservativeMax =  TMath::MaxElement(4,consval) - correction;
      fgFcConservative->SetPoint(ibin,x,correction); // i,x,y
      fgFcConservative->SetPointError(ibin,(binwidths[ibin-1]/2.),(binwidths[ibin-1]/2.),uncConservativeMin,uncConservativeMax); // i,xl,xh,yl,yh
      if( !(correction>0.) ){
	fgFcExtreme->SetPoint(ibin,x,0.); // i,x,y
	fgFcExtreme->SetPointError(ibin,(binwidths[ibin-1]/2.),(binwidths[ibin-1]/2.),0.,0.); // i,xl,xh,yl,yh
	fgFcConservative->SetPoint(ibin,x,0.); // i,x,y
	fgFcConservative->SetPointError(ibin,(binwidths[ibin-1]/2.),(binwidths[ibin-1]/2.),0.,0.); // i,xl,xh,yl,yh
      }

      Double_t valStatEffc[2] = { correctionConservativeAUncStatEffc/correctionConservativeA, 
				  correctionConservativeBUncStatEffc/correctionConservativeB };
      Double_t valStatEffb[2] = { correctionConservativeAUncStatEffb/correctionConservativeA, 
				  correctionConservativeBUncStatEffb/correctionConservativeB };
      Double_t uncConservativeStatEffc = TMath::MaxElement(2,valStatEffc);
      Double_t uncConservativeStatEffb = TMath::MaxElement(2,valStatEffb);
      fhStatUncEffcFD->SetBinContent(ibin,0.); fhStatUncEffcFD->SetBinError(ibin,uncConservativeStatEffc*100.);
      fhStatUncEffbFD->SetBinContent(ibin,0.); fhStatUncEffbFD->SetBinError(ibin,uncConservativeStatEffb*100.);
      //      cout << " pt "<< fhStatUncEffcFD->GetBinCenter(ibin) <<" bin "<< ibin<<" fc-stat-c ="<<uncConservativeStatEffc<<" fc-stat-b ="<<uncConservativeStatEffb<<endl;
      
    }

  }
  delete [] binwidths;
  delete [] limits;

}

//_________________________________________________________________________________________________________
void AliHFPtSpectrum::CalculateFeedDownCorrectedSpectrumFc(){
  //
  // Compute the feed-down corrected spectrum if feed-down correction is done via fc factor (bin by bin)
  //    physics = reco * fc / bin-width
  //
  //    uncertainty:             (stat) delta_physics = physics * sqrt ( (delta_reco/reco)^2 )
  //               (syst but feed-down) delta_physics = physics * sqrt ( (delta_reco_syst/reco)^2 )
  //                   (feed-down syst) delta_physics = physics * sqrt ( (delta_fc/fc)^2 )
  //
  //    ( Calculation done bin by bin )
  //
  //  In addition, in HIC the feed-down correction varies with an energy loss hypothesis: Raa(c-->D) / Raa(b-->D) = Rcb

  AliInfo(" Calculating the feed-down corrected spectrum (fc method)");

  if (!fhFc || !fhRECpt) {
    AliError(" Reconstructed or fc distributions are not defined");
    return;
  }

  Int_t nbins = fhRECpt->GetNbinsX();
  Double_t value = 0., errvalue = 0., errvalueMax= 0., errvalueMin= 0.;
  Double_t valueExtremeMax= 0., valueExtremeMin= 0.;
  Double_t valueConservativeMax= 0., valueConservativeMin= 0.;
  Double_t binwidth = fhRECpt->GetBinWidth(1);
  Double_t *limits = new Double_t[nbins+1];
  Double_t *binwidths = new Double_t[nbins];
  Double_t xlow=0.;
  for (Int_t i=1; i<=nbins; i++) {
    binwidth = fhRECpt->GetBinWidth(i);
    xlow = fhRECpt->GetBinLowEdge(i);
    limits[i-1] = xlow;
    binwidths[i-1] = binwidth;
  }
  limits[nbins] = xlow + binwidth;
  
  // declare the output histograms
  fhYieldCorr = new TH1D("fhYieldCorr","corrected yield (by fc)",nbins,limits);
  fhYieldCorrMax = new TH1D("fhYieldCorrMax","max corrected yield (by fc)",nbins,limits);
  fhYieldCorrMin = new TH1D("fhYieldCorrMin","min corrected yield (by fc)",nbins,limits);  
  if(fPbPbElossHypothesis) fhYieldCorrRcb = new TH2D("fhYieldCorrRcb","corrected yield (by fc) vs Rcb Eloss hypothesis; p_{T} [GeV/c] ; Rcb Eloss hypothesis ; corrected yield",nbins,limits,800,0.,4.);
  // and the output TGraphAsymmErrors
  if (fAsymUncertainties){
    fgYieldCorr = new TGraphAsymmErrors(nbins+1);
    fgYieldCorrExtreme = new TGraphAsymmErrors(nbins+1);
    fgYieldCorrConservative = new TGraphAsymmErrors(nbins+1);
  }
  
  //
  // Do the calculation
  // 
  for (Int_t ibin=1; ibin<=nbins; ibin++) {

    // calculate the value 
    //    physics = reco * fc / bin-width
    value = (fhRECpt->GetBinContent(ibin) && fhFc->GetBinContent(ibin)) ? 
      fhRECpt->GetBinContent(ibin) * fhFc->GetBinContent(ibin) : 0. ;
    value /= fhRECpt->GetBinWidth(ibin) ;
    
    // Statistical uncertainty 
    //    (stat) delta_physics = physics * sqrt ( (delta_reco/reco)^2 )
    errvalue = (value!=0. && fhRECpt->GetBinContent(ibin) && fhRECpt->GetBinContent(ibin)!=0.) ?
      value * (fhRECpt->GetBinError(ibin)/fhRECpt->GetBinContent(ibin)) : 0. ; 

    // Calculate the systematic uncertainties
    //    (syst but feed-down) delta_physics = physics * sqrt ( (delta_reco_syst/reco)^2 )
    //        (feed-down syst) delta_physics = physics * sqrt ( (delta_fc/fc)^2 )
    //
    // Protect against null denominator. If so, define uncertainty as null
    if (fhRECpt->GetBinContent(ibin) && fhRECpt->GetBinContent(ibin)!=0.) {

      if (fAsymUncertainties) {

	// Systematics but feed-down
	if (fgRECSystematics) { 
	  errvalueMax = value * ( fgRECSystematics->GetErrorYhigh(ibin) / fhRECpt->GetBinContent(ibin) );
	  errvalueMin = value * ( fgRECSystematics->GetErrorYlow(ibin) / fhRECpt->GetBinContent(ibin) );
	}
	else { errvalueMax = 0.; errvalueMin = 0.; }
	
	// Extreme feed-down systematics
	valueExtremeMax = fhRECpt->GetBinContent(ibin) * ( fhFc->GetBinContent(ibin) + fgFcExtreme->GetErrorYhigh(ibin) ) / fhRECpt->GetBinWidth(ibin) ;
	valueExtremeMin = fhRECpt->GetBinContent(ibin) * ( fhFc->GetBinContent(ibin) - fgFcExtreme->GetErrorYlow(ibin) ) / fhRECpt->GetBinWidth(ibin) ;
	
	// Conservative feed-down systematics
	valueConservativeMax = fhRECpt->GetBinContent(ibin) * ( fhFc->GetBinContent(ibin) + fgFcConservative->GetErrorYhigh(ibin) ) / fhRECpt->GetBinWidth(ibin) ;
	valueConservativeMin = fhRECpt->GetBinContent(ibin) * ( fhFc->GetBinContent(ibin) - fgFcConservative->GetErrorYlow(ibin) ) / fhRECpt->GetBinWidth(ibin) ;

      }

    }
    else { errvalueMax = 0.; errvalueMin = 0.; }
    
    //
    // Fill in the histograms
    //
    fhYieldCorr->SetBinContent(ibin,value);
    fhYieldCorr->SetBinError(ibin,errvalue);
    //
    // Fill the histos and ntuple vs the Eloss hypothesis
    //
    if (fPbPbElossHypothesis) {
      // Loop over the Eloss hypothesis
      for (Float_t rval=0.0025; rval<4.0; rval+=0.005){
	Int_t rbin = FindTH2YBin(fhYieldCorrRcb,rval);
	Double_t fcRcbvalue = fhFcRcb->GetBinContent(ibin,rbin);
	//    physics = reco * fcRcb / bin-width
	Double_t Rcbvalue = (fhRECpt->GetBinContent(ibin) && fcRcbvalue) ? 
	  fhRECpt->GetBinContent(ibin) * fcRcbvalue : 0. ;
	Rcbvalue /= fhRECpt->GetBinWidth(ibin) ;
	fhYieldCorrRcb->Fill( fhYieldCorr->GetBinCenter(ibin) , rval, Rcbvalue );
	// 	  cout << " pt "<< fhRECpt->GetBinCenter(ibin) <<" bin "<< ibin<<" rval="<<rval<<", rbin="<<rbin<<" fc-fcRbvalue="<<fcRcbvalue<<", yield="<<Rcbvalue<<endl;
      }
    }
    if (fAsymUncertainties) {
      Double_t center = fhYieldCorr->GetBinCenter(ibin);
      fgYieldCorr->SetPoint(ibin,center,value); // i,x,y
      fgYieldCorr->SetPointError(ibin,(binwidths[ibin-1]/2.),(binwidths[ibin-1]/2.),errvalueMin,errvalueMax); // i,xl,xh,yl,yh
      fhYieldCorrMax->SetBinContent(ibin,value+errvalueMax); 
      fhYieldCorrMin->SetBinContent(ibin,value-errvalueMin);
      fgYieldCorrExtreme->SetPoint(ibin,center,value);
      fgYieldCorrExtreme->SetPointError(ibin,(binwidths[ibin-1]/2.),(binwidths[ibin-1]/2.),value-valueExtremeMin,valueExtremeMax-value);
      fgYieldCorrConservative->SetPoint(ibin,center,value);
      fgYieldCorrConservative->SetPointError(ibin,(binwidths[ibin-1]/2.),(binwidths[ibin-1]/2.),value-valueConservativeMin,valueConservativeMax-value);
    }

  }
  delete [] binwidths;
  delete [] limits;
  
}

//_________________________________________________________________________________________________________
void AliHFPtSpectrum::CalculateFeedDownCorrectedSpectrumNb(Double_t deltaY, Double_t branchingRatioBintoFinalDecay) {
  //
  // Compute the feed-down corrected spectrum if feed-down correction is done via Nb (bin by bin)
  //    physics =  [ reco  - (lumi * delta_y * BR_b * eff_trig * eff_b * Nb_th) ] / bin-width
  //
  //    uncertainty:   (stat)  delta_physics = sqrt ( (delta_reco)^2 )  / bin-width
  //     (syst but feed-down)  delta_physics = sqrt ( (delta_reco_syst)^2 )  / bin-width
  //         (feed-down syst)  delta_physics = sqrt ( (k*delta_lumi/lumi)^2 + (k*delta_eff_trig/eff_trig)^2 
  //                                                   + (k*delta_Nb/Nb)^2 + (k*delta_eff/eff)^2  + (k*global_eff_ratio)^2 ) / bin-width
  //                    where k = lumi * delta_y * BR_b * eff_trig * eff_b * Nb_th
  //
  //  In addition, in HIC the feed-down correction varies with an energy loss hypothesis: Raa(b-->D) = Rb
  //    physics =  [ reco  - ( Tab * Nevt * delta_y * BR_b * eff_trig * eff_b * Nb_th * Rb ) ] / bin-width
  //
  AliInfo("Calculating the feed-down correction factor and spectrum (Nb method)");

  Int_t nbins = fhRECpt->GetNbinsX();
  Double_t binwidth = fhRECpt->GetBinWidth(1);
  Double_t value = 0., errvalue = 0., errvalueMax = 0., errvalueMin = 0., kfactor = 0.;
  Double_t errvalueExtremeMax = 0., errvalueExtremeMin = 0.;
  Double_t *limits = new Double_t[nbins+1];
  Double_t *binwidths = new Double_t[nbins];
  Double_t xlow=0.;
  for (Int_t i=1; i<=nbins; i++) {
    binwidth = fhRECpt->GetBinWidth(i);
    xlow = fhRECpt->GetBinLowEdge(i);
    limits[i-1] = xlow;
    binwidths[i-1] = binwidth;
  }
  limits[nbins] = xlow + binwidth;
  
  // declare the output histograms
  fhYieldCorr = new TH1D("fhYieldCorr","corrected yield (by Nb)",nbins,limits);
  fhYieldCorrMax = new TH1D("fhYieldCorrMax","max corrected yield (by Nb)",nbins,limits);
  fhYieldCorrMin = new TH1D("fhYieldCorrMin","min corrected yield (by Nb)",nbins,limits);
  if(fPbPbElossHypothesis) {
    fhFcRcb = new TH2D("fhFcRcb","fc correction factor (Nb method) vs Rb Eloss hypothesis; p_{T} [GeV/c] ; Rb Eloss hypothesis ; fc correction",nbins,limits,800,0.,4.);
    fhYieldCorrRcb = new TH2D("fhYieldCorrRcb","corrected yield (by Nb) vs Rb Eloss hypothesis; p_{T} [GeV/c] ; Rb Eloss hypothesis ; corrected yield",nbins,limits,800,0.,4.);
  }
  // and the output TGraphAsymmErrors
  if (fAsymUncertainties){
    fgYieldCorr = new TGraphAsymmErrors(nbins+1);
    fgYieldCorrExtreme = new TGraphAsymmErrors(nbins+1);
    fgYieldCorrConservative = new TGraphAsymmErrors(nbins+1);
    // Define fc-conservative 
    fgFcConservative = new TGraphAsymmErrors(nbins+1);
    AliInfo(" Beware the conservative & extreme uncertainties are equal by definition !");
  }

  // variables to define fc-conservative 
  double correction=0, correctionMax=0., correctionMin=0.;

  fhStatUncEffcFD = new TH1D("fhStatUncEffcFD","direct charm stat unc on the feed-down correction",nbins,limits);
  fhStatUncEffbFD = new TH1D("fhStatUncEffbFD","secondary charm stat unc on the feed-down correction",nbins,limits);
  Double_t correctionUncStatEffc=0.;
  Double_t correctionUncStatEffb=0.;


  //
  // Do the calculation
  // 
  for (Int_t ibin=1; ibin<=nbins; ibin++) {

    // Calculate the value
    //    physics =  [ reco  - (lumi * delta_y * BR_b * eff_trig * eff_b * Nb_th) ] / bin-width
    //  In HIC :   physics =  [ reco  - ( Tab * Nevt * delta_y * BR_b * eff_trig * eff_b * Nb_th * Rb ) ] / bin-width
    //
    //
    Double_t frac = 1.0, errfrac =0.;
    if(fPbPbElossHypothesis) {
      frac = fTab[0]*fNevts; 
      errfrac = frac * TMath::Sqrt( (fTab[1]/fTab[0])*(fTab[1]/fTab[0]) + (1/fNevts) );
    } else {
      frac = fLuminosity[0]; 
      errfrac = fLuminosity[1];
    }
    
    value = ( fhRECpt->GetBinContent(ibin)>0. && fhRECpt->GetBinContent(ibin)!=0. && 
	      fhFeedDownMCpt->GetBinContent(ibin)>0. && fhFeedDownEffpt->GetBinContent(ibin)>0. ) ?
      fhRECpt->GetBinContent(ibin) - frac*(deltaY*branchingRatioBintoFinalDecay*fTrigEfficiency[0]*fhFeedDownEffpt->GetBinContent(ibin)*fhFeedDownMCpt->GetBinContent(ibin) * fhRECpt->GetBinWidth(ibin) ) 
      : 0. ;
    value /= fhRECpt->GetBinWidth(ibin);
    if (value<0.) value =0.;

    //  Statistical uncertainty:   delta_physics = sqrt ( (delta_reco)^2 )  / bin-width
    errvalue = (value!=0. && fhRECpt->GetBinError(ibin) && fhRECpt->GetBinError(ibin)!=0.)  ? 
      fhRECpt->GetBinError(ibin) : 0.;
    errvalue /= fhRECpt->GetBinWidth(ibin);

    // Correction (fc) : Estimate of the relative amount feed-down subtracted
    // correction =  [ 1  - (lumi * delta_y * BR_b * eff_trig * eff_b * Nb_th)/reco ] 
    // in HIC: correction =  [ 1  - ( Tab * Nevt * delta_y * BR_b * eff_trig * eff_b * Nb_th)/reco ]
    correction = (value>0.) ? 
      1 - (frac*deltaY*branchingRatioBintoFinalDecay*fTrigEfficiency[0]*fhFeedDownEffpt->GetBinContent(ibin)*fhFeedDownMCpt->GetBinContent(ibin) * fhRECpt->GetBinWidth(ibin) ) / fhRECpt->GetBinContent(ibin)  : 0. ;
    if (correction<0.) correction = 0.;

    // Systematic uncertainties
    //     (syst but feed-down)  delta_physics = sqrt ( (delta_reco_syst)^2 )  / bin-width
    //         (feed-down syst)  delta_physics = sqrt ( (k*delta_lumi/lumi)^2 + (k*delta_eff_trig/eff_trig)^2 
    //                                                   + (k*delta_Nb/Nb)^2 + (k*delta_eff/eff)^2 + (k*global_eff_ratio)^2 ) / bin-width
    //                    where k = lumi * delta_y * BR_b * eff_trig * eff_b * Nb_th * bin-width
    kfactor = frac*deltaY*branchingRatioBintoFinalDecay*fTrigEfficiency[0]*fhFeedDownEffpt->GetBinContent(ibin)*fhFeedDownMCpt->GetBinContent(ibin) * fhRECpt->GetBinWidth(ibin) ;
    //
    if (fAsymUncertainties && value>0.) {
      Double_t nb =  fhFeedDownMCpt->GetBinContent(ibin);
      Double_t nbDmax = fhFeedDownMCptMax->GetBinContent(ibin) - fhFeedDownMCpt->GetBinContent(ibin);
      Double_t nbDmin = fhFeedDownMCpt->GetBinContent(ibin) - fhFeedDownMCptMin->GetBinContent(ibin);

      // Systematics but feed-down
      if (fgRECSystematics){
	errvalueMax = fgRECSystematics->GetErrorYhigh(ibin) / fhRECpt->GetBinWidth(ibin) ;
	errvalueMin = fgRECSystematics->GetErrorYlow(ibin) / fhRECpt->GetBinWidth(ibin);
      }
      else { errvalueMax = 0.; errvalueMin = 0.; }
  
      // Feed-down systematics
      // min value with the maximum Nb
      errvalueExtremeMin = TMath::Sqrt( ( (kfactor*errfrac/frac)*(kfactor*errfrac/frac) ) +
					( (kfactor*fTrigEfficiency[1]/fTrigEfficiency[0])*(kfactor*fTrigEfficiency[1]/fTrigEfficiency[0]) ) +
					( (kfactor*nbDmax/nb)*(kfactor*nbDmax/nb) )  +
					( (kfactor*fhFeedDownEffpt->GetBinError(ibin)/fhFeedDownEffpt->GetBinContent(ibin))*(kfactor*fhFeedDownEffpt->GetBinError(ibin)/fhFeedDownEffpt->GetBinContent(ibin)) ) +
					( (kfactor*fGlobalEfficiencyUncertainties[1])*(kfactor*fGlobalEfficiencyUncertainties[1]) ) 
					) / fhRECpt->GetBinWidth(ibin);
      // max value with the minimum Nb
      errvalueExtremeMax =  TMath::Sqrt( ( (kfactor*errfrac/frac)*(kfactor*errfrac/frac) ) +
					 ( (kfactor*fTrigEfficiency[1]/fTrigEfficiency[0])*(kfactor*fTrigEfficiency[1]/fTrigEfficiency[0]) ) +
					 ( (kfactor*nbDmin/nb)*(kfactor*nbDmin/nb) )  +
					 ( (kfactor*fhFeedDownEffpt->GetBinError(ibin)/fhFeedDownEffpt->GetBinContent(ibin))*(kfactor*fhFeedDownEffpt->GetBinError(ibin)/fhFeedDownEffpt->GetBinContent(ibin))	) +
					 ( (kfactor*fGlobalEfficiencyUncertainties[1])*(kfactor*fGlobalEfficiencyUncertainties[1]) )
					 ) / fhRECpt->GetBinWidth(ibin);

      // Correction systematics (fc)
      // min value with the maximum Nb
      correctionMin = TMath::Sqrt( ( (kfactor*errfrac/frac)*(kfactor*errfrac/frac) ) + 
				   ( (kfactor*fTrigEfficiency[1]/fTrigEfficiency[0])*(kfactor*fTrigEfficiency[1]/fTrigEfficiency[0]) ) +
				   ( (kfactor*nbDmax/nb)*(kfactor*nbDmax/nb) )  +
				   ( (kfactor*fhFeedDownEffpt->GetBinError(ibin)/fhFeedDownEffpt->GetBinContent(ibin))*(kfactor*fhFeedDownEffpt->GetBinError(ibin)/fhFeedDownEffpt->GetBinContent(ibin)) ) +
				   ( (kfactor*fGlobalEfficiencyUncertainties[1])*(kfactor*fGlobalEfficiencyUncertainties[1]) ) 
				   ) / fhRECpt->GetBinContent(ibin) ;
      // max value with the minimum Nb
      correctionMax =  TMath::Sqrt( ( (kfactor*errfrac/frac)*(kfactor*errfrac/frac) ) + 
				    ( (kfactor*fTrigEfficiency[1]/fTrigEfficiency[0])*(kfactor*fTrigEfficiency[1]/fTrigEfficiency[0]) ) +
				    ( (kfactor*nbDmin/nb)*(kfactor*nbDmin/nb) )  +
				    ( (kfactor*fhFeedDownEffpt->GetBinError(ibin)/fhFeedDownEffpt->GetBinContent(ibin))*(kfactor*fhFeedDownEffpt->GetBinError(ibin)/fhFeedDownEffpt->GetBinContent(ibin))	) +
				    ( (kfactor*fGlobalEfficiencyUncertainties[1])*(kfactor*fGlobalEfficiencyUncertainties[1]) )
				    ) / fhRECpt->GetBinContent(ibin) ;
      //
      correctionUncStatEffb = TMath::Sqrt(  ( (kfactor*fhFeedDownEffpt->GetBinError(ibin)/fhFeedDownEffpt->GetBinContent(ibin))*(kfactor*fhFeedDownEffpt->GetBinError(ibin)/fhFeedDownEffpt->GetBinContent(ibin))	)
					    ) / fhRECpt->GetBinContent(ibin) ;
      correctionUncStatEffc = 0.;
    }
    else{ // Don't consider Nb uncertainty in this case [ to be tested!!! ]
      errvalueExtremeMax =  TMath::Sqrt( ( (kfactor*errfrac/frac)*(kfactor*errfrac/frac) ) +
					 ( (kfactor*fTrigEfficiency[1]/fTrigEfficiency[0])*(kfactor*fTrigEfficiency[1]/fTrigEfficiency[0]) )  +
					 ( (kfactor*fhFeedDownEffpt->GetBinError(ibin)/fhFeedDownEffpt->GetBinContent(ibin))*(kfactor*fhFeedDownEffpt->GetBinError(ibin)/fhFeedDownEffpt->GetBinContent(ibin))	)  +
					 ( (kfactor*fGlobalEfficiencyUncertainties[1])*(kfactor*fGlobalEfficiencyUncertainties[1]) )
					 ) / fhRECpt->GetBinWidth(ibin);
      errvalueExtremeMin =  errvalueExtremeMax ;
    }
    

    // fill in histograms
    fhYieldCorr->SetBinContent(ibin,value);
    fhYieldCorr->SetBinError(ibin,errvalue);    
    //
    // Estimate how the result varies vs charm/beauty Eloss hypothesis
    //
    if ( correction>0.0001 && fPbPbElossHypothesis){
      // Loop over the Eloss hypothesis
      //      Int_t rbin=0;
      for (Float_t rval=0.0025; rval<4.0; rval+=0.005){
	// correction =  [ 1  - (Tab *Nevt * delta_y * BR_b * eff_trig * eff_b * Nb_th *binwidth )* (rval) /reco ] 
	Double_t fcRcbvalue = 1 - (fTab[0]*fNevts*deltaY*branchingRatioBintoFinalDecay*fTrigEfficiency[0]*fhFeedDownEffpt->GetBinContent(ibin)*fhFeedDownMCpt->GetBinContent(ibin)*fhRECpt->GetBinWidth(ibin) * rval ) / fhRECpt->GetBinContent(ibin) ;
	if(fcRcbvalue<0.) fcRcbvalue=0.;
	fhFcRcb->Fill( fhRECpt->GetBinCenter(ibin) , rval, fcRcbvalue );
	//    physics = reco * fcRcb / bin-width
	Double_t Rcbvalue = (fhRECpt->GetBinContent(ibin) && fcRcbvalue) ? 
	  fhRECpt->GetBinContent(ibin) * fcRcbvalue : 0. ;
	Rcbvalue /= fhRECpt->GetBinWidth(ibin) ;
	fhYieldCorrRcb->Fill( fhYieldCorr->GetBinCenter(ibin) , rval, Rcbvalue );
	// 	if(ibin==3){
	// 	  cout << " pt "<< fhFcRcb->GetBinCenter(ibin) <<" bin "<< ibin<<" rval="<<rval<<", rbin="<<rbin<<", fc-Rb-value="<< fcRcbvalue << ", yield-Rb-value="<< Rcbvalue <<endl;
	//	cout << " pt "<< fhFcRcb->GetBinCenter(ibin) <<" bin "<< ibin<<" rval="<<rval<<", fc-Rb-value="<< fcRcbvalue << ", yield-Rb-value="<< Rcbvalue <<endl;
	// 	  rbin++;
	// 	}
      }
    }
    //
    // Fill the rest of (asymmetric) histograms
    //
    if (fAsymUncertainties) {
      Double_t x = fhYieldCorr->GetBinCenter(ibin);
      fgYieldCorr->SetPoint(ibin,x,value); // i,x,y
      fgYieldCorr->SetPointError(ibin,(binwidths[ibin-1]/2.),(binwidths[ibin-1]/2.),errvalueMin,errvalueMax); // i,xl,xh,yl,yh
      fhYieldCorrMax->SetBinContent(ibin,value+errvalueMax); 
      fhYieldCorrMin->SetBinContent(ibin,value-errvalueMin);
      fgYieldCorrExtreme->SetPoint(ibin,x,value); // i,x,y
      fgYieldCorrExtreme->SetPointError(ibin,(binwidths[ibin-1]/2.),(binwidths[ibin-1]/2.),errvalueExtremeMin,errvalueExtremeMax); // i,xl,xh,yl,yh
      fgYieldCorrConservative->SetPoint(ibin,x,value); // i,x,y
      fgYieldCorrConservative->SetPointError(ibin,(binwidths[ibin-1]/2.),(binwidths[ibin-1]/2.),errvalueExtremeMin,errvalueExtremeMax); // i,xl,xh,yl,yh
      //      cout << " bin " << ibin << ", correction " << correction << ", min correction unc " << correctionMin << ", max correction unc " << correctionMax << endl;
      if(correction>0.){
	fgFcConservative->SetPoint(ibin,x,correction);
	fgFcConservative->SetPointError(ibin,(binwidths[ibin-1]/2.),(binwidths[ibin-1]/2.),correctionMin,correctionMax);
	
	fhStatUncEffbFD->SetBinContent(ibin,0.); fhStatUncEffbFD->SetBinError(ibin,correctionUncStatEffb/correction*100.);
	fhStatUncEffcFD->SetBinContent(ibin,0.); fhStatUncEffcFD->SetBinError(ibin,correctionUncStatEffc/correction*100.);
	//	cout << " pt "<< fhStatUncEffcFD->GetBinCenter(ibin) <<" bin "<< ibin<<" fc-stat-c ="<< correctionUncStatEffc/correction <<" fc-stat-b ="<< correctionUncStatEffb/correction <<endl;
      }
      else{
	fgFcConservative->SetPoint(ibin,x,0.);
	fgFcConservative->SetPointError(ibin,(binwidths[ibin-1]/2.),(binwidths[ibin-1]/2.),0.,0.);
      }
    }

  }
  delete [] binwidths;
  delete [] limits;
  
}


//_________________________________________________________________________________________________________
void AliHFPtSpectrum::ComputeSystUncertainties(AliHFSystErr *systematics, Bool_t combineFeedDown) {
  //
  // Function that re-calculates the global systematic uncertainties
  //   by calling the class AliHFSystErr and combining those
  //   (in quadrature) with the feed-down subtraction uncertainties
  //

  // Estimate the feed-down uncertainty in percentage
  Int_t nentries = fgSigmaCorrConservative->GetN();
  TGraphAsymmErrors *grErrFeeddown = new TGraphAsymmErrors(nentries);
  Double_t x=0., y=0., errx=0., erryl=0., erryh=0;
  for(Int_t i=0; i<nentries; i++) {
    x=0.; y=0.; errx=0.; erryl=0.; erryh=0.;
    fgSigmaCorrConservative->GetPoint(i,x,y);
    if(y>0.){
      errx = fgSigmaCorrConservative->GetErrorXlow(i) ;
      erryl = fgSigmaCorrConservative->GetErrorYlow(i) / y ;
      erryh = fgSigmaCorrConservative->GetErrorYhigh(i) / y ;
    }
    //    cout << " x "<< x << " +- "<<errx<<" , y "<<y<<" + "<<erryh<<" - "<<erryl<<endl; 
    grErrFeeddown->SetPoint(i,x,0.);
    grErrFeeddown->SetPointError(i,errx,errx,erryl,erryh); //i, xl, xh, yl, yh
  }

  // Draw all the systematics independently
  systematics->DrawErrors(grErrFeeddown);

  // Set the sigma systematic uncertainties
  // possibly combine with the feed-down uncertainties 
  Double_t errylcomb=0., erryhcomb=0;
  for(Int_t i=1; i<nentries; i++) {
    fgSigmaCorr->GetPoint(i,x,y);
    errx = grErrFeeddown->GetErrorXlow(i) ;
    erryl = grErrFeeddown->GetErrorYlow(i);
    erryh = grErrFeeddown->GetErrorYhigh(i);
    if (combineFeedDown) {
      errylcomb = systematics->GetTotalSystErr(x,erryl) * y ;
      erryhcomb = systematics->GetTotalSystErr(x,erryh) * y ;
    } else {
      errylcomb = systematics->GetTotalSystErr(x) * y ;
      erryhcomb = systematics->GetTotalSystErr(x) * y ;
    }
    fgSigmaCorr->SetPointError(i,errx,errx,errylcomb,erryhcomb);
    //
    fhSigmaCorrDataSyst->SetBinContent(i,y);
    erryl = systematics->GetTotalSystErr(x) * y ;
    fhSigmaCorrDataSyst->SetBinError(i,erryl);
  }

}


//_________________________________________________________________________________________________________
void AliHFPtSpectrum::DrawSpectrum(TGraphAsymmErrors *gPrediction) {
  //
  // Example method to draw the corrected spectrum & the theoretical prediction
  //

  TCanvas *csigma = new TCanvas("csigma","Draw the corrected cross-section & the prediction");
  csigma->SetFillColor(0);
  gPrediction->GetXaxis()->SetTitleSize(0.05);
  gPrediction->GetXaxis()->SetTitleOffset(0.95);
  gPrediction->GetYaxis()->SetTitleSize(0.05);
  gPrediction->GetYaxis()->SetTitleOffset(0.95);
  gPrediction->GetXaxis()->SetTitle("p_{T}  [GeV]");
  gPrediction->GetYaxis()->SetTitle("BR #times #frac{d#sigma}{dp_{T}} |_{|y|<0.5}   [pb/GeV]");
  gPrediction->SetLineColor(kGreen+2);
  gPrediction->SetLineWidth(3);
  gPrediction->SetFillColor(kGreen+1);
  gPrediction->Draw("3CA");
  fgSigmaCorr->SetLineColor(kRed);
  fgSigmaCorr->SetLineWidth(1);
  fgSigmaCorr->SetFillColor(kRed);
  fgSigmaCorr->SetFillStyle(0);
  fgSigmaCorr->Draw("2");
  fhSigmaCorr->SetMarkerColor(kRed);
  fhSigmaCorr->Draw("esame");
  csigma->SetLogy();
  TLegend * leg = new TLegend(0.7,0.75,0.87,0.5);
  leg->SetBorderSize(0);
  leg->SetLineColor(0);
  leg->SetFillColor(0);
  leg->SetTextFont(42);
  leg->AddEntry(gPrediction,"FONLL ","fl");
  leg->AddEntry(fhSigmaCorr,"data stat. unc.","pl");
  leg->AddEntry(fgSigmaCorr,"data syst. unc.","f");
  leg->Draw();
  csigma->Draw();

}

//_________________________________________________________________________________________________________
TH1D * AliHFPtSpectrum::ReweightHisto(TH1D *hToReweight, TH1D *hReference){
  //
  // Function to  reweight histograms for testing purposes: 
  // This function takes the histo hToReweight and reweights 
  //  it (its pt shape) with respect to hReference 
  // 

  // check histograms consistency
  Bool_t areconsistent=kTRUE;
  areconsistent &= CheckHistosConsistency(hToReweight,hReference);
  if (!areconsistent) {
    AliInfo("the histograms to reweight are not consistent (bin width, bounds)"); 
    return NULL;
  }

  // define a new empty histogram
  TH1D *hReweighted = (TH1D*)hToReweight->Clone("hReweighted");
  hReweighted->Reset();
  Double_t weight=1.0;
  Double_t yvalue=1.0; 
  Double_t integralRef = hReference->Integral();
  Double_t integralH = hToReweight->Integral();

  // now reweight the spectra
  //
  // the weight is the relative probability of the given pt bin in the reference histo
  //  divided by its relative probability (to normalize it) on the histo to re-weight
  for (Int_t i=0; i<=hToReweight->GetNbinsX(); i++) {
    weight = (hReference->GetBinContent(i)/integralRef) / (hToReweight->GetBinContent(i)/integralH) ;
    yvalue = hToReweight->GetBinContent(i);
    hReweighted->SetBinContent(i,yvalue*weight);
  }

  return (TH1D*)hReweighted;
}

//_________________________________________________________________________________________________________
TH1D * AliHFPtSpectrum::ReweightRecHisto(TH1D *hRecToReweight, TH1D *hMCToReweight, TH1D *hMCReference){
  //
  // Function to  reweight histograms for testing purposes: 
  // This function takes the histo hToReweight and reweights 
  //  it (its pt shape) with respect to hReference /hMCToReweight
  // 

  // check histograms consistency
  Bool_t areconsistent=kTRUE;
  areconsistent &= CheckHistosConsistency(hMCToReweight,hMCReference);
  areconsistent &= CheckHistosConsistency(hRecToReweight,hMCReference);
  if (!areconsistent) {
    AliInfo("the histograms to reweight are not consistent (bin width, bounds)"); 
    return NULL;
  }

  // define a new empty histogram
  TH1D *hReweighted = (TH1D*)hMCToReweight->Clone("hReweighted");
  hReweighted->Reset();
  TH1D *hRecReweighted = (TH1D*)hRecToReweight->Clone("hRecReweighted");
  hRecReweighted->Reset();
  Double_t weight=1.0;
  Double_t yvalue=1.0, yrecvalue=1.0; 
  Double_t integralRef = hMCReference->Integral();
  Double_t integralH = hMCToReweight->Integral();

  // now reweight the spectra
  //
  // the weight is the relative probability of the given pt bin 
  //  that should be applied in the MC histo to get the reference histo shape
  //  Probabilities are properly normalized.
  for (Int_t i=0; i<=hMCToReweight->GetNbinsX(); i++) {
    weight = (hMCReference->GetBinContent(i)/integralRef) / (hMCToReweight->GetBinContent(i)/integralH) ;
    yvalue = hMCToReweight->GetBinContent(i);
    hReweighted->SetBinContent(i,yvalue*weight);
    yrecvalue = hRecToReweight->GetBinContent(i);
    hRecReweighted->SetBinContent(i,yrecvalue*weight);
  }

  return (TH1D*)hRecReweighted;
}



//_________________________________________________________________________________________________________
Int_t AliHFPtSpectrum::FindTH2YBin(TH2D *histo, Float_t yvalue){
  //
  // Function to find the y-axis bin of a TH2 for a given y-value
  //
  
  Int_t nbins = histo->GetNbinsY();
  Int_t ybin=0;
  for (int j=0; j<=nbins; j++) {
    Float_t value = histo->GetYaxis()->GetBinCenter(j);
    Float_t width = histo->GetYaxis()->GetBinWidth(j);
    //    if( TMath::Abs(yvalue-value)<= width/2. ) {
    if( TMath::Abs(yvalue-value)<= width ) {
      ybin =j;
      //      cout <<" value "<<value << ", yval "<< yvalue<<", bin width "<<width/2.<< " y ="<<ybin<<endl;
      break;
    }
  }
  
  return ybin;
}
