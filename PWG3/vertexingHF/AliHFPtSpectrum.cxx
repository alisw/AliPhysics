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

#include <Riostream.h>

#include "TMath.h"
#include "TH1.h"
#include "TH1D.h"
#include "TGraphAsymmErrors.h"

#include "AliLog.h"
#include "AliHFPtSpectrum.h"

ClassImp(AliHFPtSpectrum)

//_________________________________________________________________________________________________________
AliHFPtSpectrum::AliHFPtSpectrum(const char* name, const char* title, Int_t option):
  TNamed(name,title),
  fhDirectMCpt(),
  fhFeedDownMCpt(),
  fhDirectMCptMax(),
  fhDirectMCptMin(),
  fhFeedDownMCptMax(),
  fhFeedDownMCptMin(),
  fhDirectEffpt(),
  fhFeedDownEffpt(),
  fhRECpt(),
  fLuminosity(),
  fTrigEfficiency(),
  fhFc(),
  fhFcMax(),
  fhFcMin(),
  fgFc(),
  fhYieldCorr(),
  fhYieldCorrMax(),
  fhYieldCorrMin(),
  fgYieldCorr(),
  fhSigmaCorr(),
  fhSigmaCorrMax(),
  fhSigmaCorrMin(),
  fgSigmaCorr(),
  fFeedDownOption(option),
  fAsymUncertainties(kTRUE)
{
  //
  // Default constructor
  //

  fLuminosity[0]=1.;  fLuminosity[1]=0.;  
  fTrigEfficiency[0]=1.; fTrigEfficiency[1]=0.; 

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
  fLuminosity(),
  fTrigEfficiency(),
  fhFc(rhs.fhFc),
  fhFcMax(rhs.fhFcMax),
  fhFcMin(rhs.fhFcMin),
  fgFc(rhs.fgFc),
  fhYieldCorr(rhs.fhYieldCorr),
  fhYieldCorrMax(rhs.fhYieldCorrMax),
  fhYieldCorrMin(rhs.fhYieldCorrMin),
  fgYieldCorr(rhs.fgYieldCorr),
  fhSigmaCorr(rhs.fhSigmaCorr),
  fhSigmaCorrMax(rhs.fhSigmaCorrMax),
  fhSigmaCorrMin(rhs.fhSigmaCorrMin),
  fgSigmaCorr(rhs.fgSigmaCorr),
  fFeedDownOption(rhs.fFeedDownOption),
  fAsymUncertainties(rhs.fAsymUncertainties)
{
  //
  // Copy constructor
  //

  for(Int_t i=0; i<2; i++){
    fLuminosity[i] = rhs.fLuminosity[i];
    fTrigEfficiency[i] = rhs.fTrigEfficiency[i];
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
  fhFc = source.fhFc;
  fhFcMax = source.fhFcMax;
  fhFcMin = source.fhFcMin;
  fgFc = source.fgFc;
  fhYieldCorr = source.fhYieldCorr;
  fhYieldCorrMax = source.fhYieldCorrMax;
  fhYieldCorrMin = source.fhYieldCorrMin;
  fgYieldCorr = source.fgYieldCorr;
  fhSigmaCorr = source.fhSigmaCorr;
  fhSigmaCorrMax = source.fhSigmaCorrMax;
  fhSigmaCorrMin = source.fhSigmaCorrMin;
  fgSigmaCorr = source.fgSigmaCorr;
  fFeedDownOption = source.fFeedDownOption;
  fAsymUncertainties = source.fAsymUncertainties;
  
  for(Int_t i=0; i<2; i++){
    fLuminosity[i] = source.fLuminosity[i];
    fTrigEfficiency[i] = source.fTrigEfficiency[i];
  }

  return *this;
}

//_________________________________________________________________________________________________________
AliHFPtSpectrum::~AliHFPtSpectrum(){
  //
  // Destructor
  //
  ;
}
  

//_________________________________________________________________________________________________________
void AliHFPtSpectrum::SetMCptSpectra(TH1 *hDirect, TH1 *hFeedDown){
  //
  // Set the MonteCarlo or Theoretical spectra
  //  both for direct and feed-down contributions
  //
  
  if (!hDirect || !hFeedDown) {
    AliError("One or both (direct, feed-down) spectra don't exist");
    return;
  }

  Bool_t areconsistent = kTRUE;
  areconsistent = CheckHistosConsistency(hDirect,hFeedDown);
  if (!areconsistent) {
    AliInfo("Histograms are not consistent (bin width, bounds)"); 
    return;
  }

  fhDirectMCpt = hDirect;
  fhFeedDownMCpt = hFeedDown;
}

//_________________________________________________________________________________________________________
void AliHFPtSpectrum::SetFeedDownMCptSpectra(TH1 *hFeedDown){
  //
  // Set the MonteCarlo or Theoretical spectra
  //  for feed-down contribution
  //
  
  if (!hFeedDown) {
    AliError("Feed-down spectra don't exist");
    return;
  }
  fhFeedDownMCpt = hFeedDown;
}

//_________________________________________________________________________________________________________
void AliHFPtSpectrum::SetMCptDistributionsBounds(TH1 *hDirectMax, TH1 *hDirectMin, TH1 *hFeedDownMax, TH1 *hFeedDownMin){
  //
  // Set the maximum and minimum MonteCarlo or Theoretical spectra
  //  both for direct and feed-down contributions
  // used in case uncertainties are asymmetric and ca not be on the "basic histograms"
  //

  if (!hDirectMax || !hDirectMin || !hFeedDownMax|| !hFeedDownMin) {
    AliError("One or all of the max/min direct/feed-down spectra don't exist");
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

  fhDirectMCptMax = hDirectMax;
  fhDirectMCptMin = hDirectMin;
  fhFeedDownMCptMax = hFeedDownMax;
  fhFeedDownMCptMin = hFeedDownMin;
}

//_________________________________________________________________________________________________________
void AliHFPtSpectrum::SetFeedDownMCptDistributionsBounds(TH1 *hFeedDownMax, TH1 *hFeedDownMin){
  //
  // Set the maximum and minimum MonteCarlo or Theoretical spectra
  //   for feed-down contributions
  // used in case uncertainties are asymmetric and can not be on the "basic histogram"
  //

  if (!hFeedDownMax || !hFeedDownMin) {
    AliError("One or all of the max/min direct/feed-down spectra don't exist");
    return;
  }

  Bool_t areconsistent = kTRUE; 
  areconsistent &= CheckHistosConsistency(hFeedDownMax,hFeedDownMin);
  if (!areconsistent) {
    AliInfo("Histograms are not consistent (bin width, bounds)"); 
    return;
  }

  fhFeedDownMCptMax = hFeedDownMax;
  fhFeedDownMCptMin = hFeedDownMin;
}

//_________________________________________________________________________________________________________
void AliHFPtSpectrum::SetDirectAccEffCorrection(TH1 *hDirectEff){
  //
  // Set the Acceptance and Efficiency corrections 
  //   for the direct contribution
  //
  
  if (!hDirectEff) {
    AliError("The direct acceptance and efficiency corrections doesn't exist");
    return;
  }

  fhDirectEffpt = hDirectEff;
}

//_________________________________________________________________________________________________________
void AliHFPtSpectrum::SetAccEffCorrection(TH1 *hDirectEff, TH1 *hFeedDownEff){
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

  fhDirectEffpt = hDirectEff;
  fhFeedDownEffpt = hFeedDownEff;
}

//_________________________________________________________________________________________________________
void AliHFPtSpectrum::SetReconstructedSpectrum(TH1 *hRec) {
  //
  // Set the reconstructed spectrum
  //
  
  if (!hRec) {
    AliError("The reconstructed spectrum doesn't exist");
    return;
  }

  fhRECpt = hRec;
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
  // Uncertainties: delta_sigma = sigma * sqrt ( (delta_reco/reco)^2 + (delta_lumi/lumi)^2 + (delta_eff_trig/eff_trig)^2  )

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
    fhYieldCorr = fhRECpt;
    fhYieldCorr->SetNameTitle("fhYieldCorr","un-corrected yield");
    fhYieldCorrMax = fhRECpt;
    fhYieldCorrMin = fhRECpt;
    fhYieldCorrMax->SetNameTitle("fhYieldCorrMax","un-corrected yield");
    fhYieldCorrMin->SetNameTitle("fhYieldCorrMin","un-corrected yield");
    fAsymUncertainties=kFALSE;
  }
  else { 
    AliInfo(" Are you sure the feed-down correction option is right ?"); 
  }

  // Print out information
  printf("\n\n     Correcting the spectra with : \n   luminosity = %2.2e +- %2.2e, trigger efficiency = %2.2e +- %2.2e, \n    delta_y = %2.2f, BR_c = %2.2e, BR_b_decay = %2.2e \n\n",fLuminosity[0],fLuminosity[1],fTrigEfficiency[0],fTrigEfficiency[1],deltaY,branchingRatioC,branchingRatioBintoFinalDecay);

  //
  // Finally: Correct from yields to cross-section
  //
  Int_t nbins = fhRECpt->GetNbinsX();
  Double_t binwidth = fhRECpt->GetBinWidth(1);
  Double_t xmin = fhRECpt->GetBinCenter(1) - (binwidth/2.) ; 
  Double_t xmax = fhRECpt->GetBinCenter(nbins) + (binwidth/2.) ; 
  
  // declare the output histograms
  TH1D *hSigmaCorr = new TH1D("hSigmaCorr","corrected sigma",nbins,xmin,xmax);
  TH1D *hSigmaCorrMax = new TH1D("hSigmaCorrMax","max corrected sigma",nbins,xmin,xmax);
  TH1D *hSigmaCorrMin = new TH1D("hSigmaCorrMin","min corrected sigma",nbins,xmin,xmax);
  // and the output TGraphAsymmErrors
  if (fAsymUncertainties & !fgSigmaCorr) fgSigmaCorr = new TGraphAsymmErrors(nbins);

  // protect against null denominator
  if (deltaY==0. || fLuminosity[0]==0. || fTrigEfficiency[0]==0. || branchingRatioC==0.) {
    AliError(" Hey you ! Why luminosity or trigger-efficiency or the c-BR or delta_y are set to zero ?! ");
    return ;
  }

  Double_t value=0, valueMax=0., valueMin=0.;
  for(Int_t ibin=0; ibin<=nbins; ibin++){

    //   Sigma = ( 1. / (lumi * delta_y * BR_c * eff_trig * eff_c ) ) * spectra (corrected for feed-down)
    value = (fhDirectEffpt->GetBinContent(ibin) && fhDirectEffpt->GetBinContent(ibin)!=0.) ? 
      ( fhYieldCorr->GetBinContent(ibin) / ( deltaY * branchingRatioC * fLuminosity[0] * fTrigEfficiency[0] * fhDirectEffpt->GetBinContent(ibin) ) )
      : 0. ;
    
    // Uncertainties: delta_sigma = sigma * sqrt ( (delta_reco/reco)^2 + (delta_lumi/lumi)^2 + (delta_eff_trig/eff_trig)^2  )
    if (fAsymUncertainties) {
      valueMax = value * TMath::Sqrt( (fgYieldCorr->GetErrorYhigh(ibin)/fhYieldCorr->GetBinContent(ibin))* (fgYieldCorr->GetErrorYhigh(ibin)/fhYieldCorr->GetBinContent(ibin)) + 
				       (fLuminosity[1]/fLuminosity[0])*(fLuminosity[1]/fLuminosity[0]) + 
				       (fTrigEfficiency[1]/fTrigEfficiency[0])*(fTrigEfficiency[1]/fTrigEfficiency[0])  );
      valueMin = value * TMath::Sqrt( (fgYieldCorr->GetErrorYlow(ibin)/fhYieldCorr->GetBinContent(ibin))* (fgYieldCorr->GetErrorYlow(ibin)/fhYieldCorr->GetBinContent(ibin)) + 
				       (fLuminosity[1]/fLuminosity[0])*(fLuminosity[1]/fLuminosity[0]) + 
				       (fTrigEfficiency[1]/fTrigEfficiency[0])*(fTrigEfficiency[1]/fTrigEfficiency[0])  );
    }
    else {
      // protect against null denominator
      valueMax = (value!=0.) ?
	value * TMath::Sqrt( (fhYieldCorr->GetBinError(ibin)/fhYieldCorr->GetBinContent(ibin))* (fhYieldCorr->GetBinError(ibin)/fhYieldCorr->GetBinContent(ibin)) + 
			     (fLuminosity[1]/fLuminosity[0])*(fLuminosity[1]/fLuminosity[0]) + 
			     (fTrigEfficiency[1]/fTrigEfficiency[0])*(fTrigEfficiency[1]/fTrigEfficiency[0])  )
	: 0. ;
      valueMin = valueMax;
    }
    
    // Fill the histograms
    hSigmaCorr->SetBinContent(ibin,value);
    hSigmaCorr->SetBinError(ibin,valueMax);
    hSigmaCorrMax->SetBinContent(ibin,valueMax);
    hSigmaCorrMin->SetBinContent(ibin,valueMin);
    // Fill the TGraphAsymmErrors
    if (fAsymUncertainties) {
      Double_t x = fhYieldCorr->GetBinCenter(ibin);
      fgSigmaCorr->SetPoint(ibin,x,value); // i,x,y
      fgSigmaCorr->SetPointError(ibin,(binwidth/2.),(binwidth/2.),valueMin,valueMax); // i,xl,xh,yl,yh
    }
    
  }

  fhSigmaCorr = hSigmaCorr ;
  fhSigmaCorrMax = hSigmaCorrMax ;
  fhSigmaCorrMin = hSigmaCorrMin ;
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
Bool_t AliHFPtSpectrum::CheckHistosConsistency(TH1 *h1, TH1 *h2){
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
  
  // define the variables
  Int_t nbins = fhRECpt->GetNbinsX();
  Double_t binwidth = fhRECpt->GetBinWidth(1);
  Double_t xmin = fhRECpt->GetBinCenter(1) - (binwidth/2.) ; 
  Double_t xmax = fhRECpt->GetBinCenter(nbins) + (binwidth/2.) ; 
  Double_t correction=1.;
  Double_t correctionMax=1., correctionMin=1.;
  Double_t theoryRatio=1.;
  Double_t effRatio=1.; 
  
  // declare the output histograms
  TH1D *hfc = new TH1D("hfc","fc correction factor",nbins,xmin,xmax);
  TH1D *hfcMax = new TH1D("hfcMax","max fc correction factor",nbins,xmin,xmax);
  TH1D *hfcMin = new TH1D("hfcMin","min fc correction factor",nbins,xmin,xmax);
  // two local control histograms
  TH1D *hTheoryRatio = new TH1D("hTheoryRatio","Theoretical B-->D over c-->D (feed-down/direct) ratio",nbins,xmin,xmax);
  TH1D *hEffRatio = new TH1D("hEffRatio","Efficiency B-->D over c-->D (feed-down/direct) ratio",nbins,xmin,xmax);
  // and the output TGraphAsymmErrors
  if (fAsymUncertainties & !fgFc) fgFc = new TGraphAsymmErrors(nbins);

  //
  // Compute fc
  //
  for (Int_t ibin=0; ibin<=nbins; ibin++) {

    //  theory_ratio = (N_b/N_c) 
    theoryRatio = (fhDirectMCpt->GetBinContent(ibin) && fhDirectMCpt->GetBinContent(ibin)!=0.) ? fhFeedDownMCpt->GetBinContent(ibin) / fhDirectMCpt->GetBinContent(ibin) : 1.0 ;
    //  eff_ratio = (eff_b/eff_c)
    effRatio = (fhDirectEffpt->GetBinContent(ibin) && fhDirectEffpt->GetBinContent(ibin)!=0.) ? fhFeedDownEffpt->GetBinContent(ibin) / fhDirectEffpt->GetBinContent(ibin) : 1.0 ;
    //   fc = 1 / ( 1 + (eff_b/eff_c)*(N_b/N_c) ) 
    correction = (effRatio && theoryRatio) ? ( 1. / ( 1 + ( effRatio * theoryRatio ) ) ) : 1.0 ;

    // Calculate the uncertainty [ considering only the theoretical uncertainties on Nb & Nc for now !!! ]
    // delta_fc = fc^2 * (Nb/Nc) * sqrt ( (delta_Nb/Nb)^2 + (delta_Nc/Nc)^2 ) 
    Double_t deltaNbMax = fhFeedDownMCptMax->GetBinContent(ibin) - fhFeedDownMCpt->GetBinContent(ibin) ;
    Double_t deltaNbMin = fhFeedDownMCpt->GetBinContent(ibin) - fhFeedDownMCptMin->GetBinContent(ibin) ;
    Double_t deltaNcMax = fhDirectMCptMax->GetBinContent(ibin) - fhDirectMCpt->GetBinContent(ibin) ;
    Double_t deltaNcMin = fhDirectMCpt->GetBinContent(ibin) - fhDirectMCptMin->GetBinContent(ibin) ;

    // Protect against null denominator. If so, define uncertainty as null
    if (fhFeedDownMCpt->GetBinContent(ibin) && fhFeedDownMCpt->GetBinContent(ibin)!=0. && 
	fhDirectMCpt->GetBinContent(ibin) && fhDirectMCpt->GetBinContent(ibin)!=0. ) {
      correctionMax = correction*correction*theoryRatio * 
	TMath::Sqrt( 
		    (deltaNbMax/fhFeedDownMCpt->GetBinContent(ibin))*(deltaNbMax/fhFeedDownMCpt->GetBinContent(ibin)) + 
		    (deltaNcMax/fhDirectMCpt->GetBinContent(ibin))*(deltaNcMax/fhDirectMCpt->GetBinContent(ibin)) 
		     );
      correctionMin = correction*correction*theoryRatio * 
	TMath::Sqrt( 
		    (deltaNbMin/fhFeedDownMCpt->GetBinContent(ibin))*(deltaNbMin/fhFeedDownMCpt->GetBinContent(ibin)) + 
		    (deltaNcMin/fhDirectMCpt->GetBinContent(ibin))*(deltaNcMin/fhDirectMCpt->GetBinContent(ibin)) 
		     );
    }
    else { correctionMax = 0.; correctionMin = 0.; }


    // Fill in the histograms
    hTheoryRatio->SetBinContent(ibin,theoryRatio);
    hEffRatio->SetBinContent(ibin,effRatio);
    hfc->SetBinContent(ibin,correction);
    hfcMax->SetBinContent(ibin,correction+correctionMax);
    hfcMin->SetBinContent(ibin,correction-correctionMin);
    if (fAsymUncertainties) {
      Double_t x = fhDirectMCpt->GetBinCenter(ibin);
      fgFc->SetPoint(ibin,x,correction); // i,x,y
      fgFc->SetPointError(ibin,(binwidth/2.),(binwidth/2.),correctionMin,correctionMax); // i,xl,xh,yl,yh
    }

  }

  fhFc = hfc;
  fhFcMax = hfcMax;
  fhFcMin = hfcMin;
}

//_________________________________________________________________________________________________________
void AliHFPtSpectrum::CalculateFeedDownCorrectedSpectrumFc(){
  //
  // Compute the feed-down corrected spectrum if feed-down correction is done via fc factor (bin by bin)
  //    physics = reco * fc 
  //
  //    uncertainty: delta_physics = physics * sqrt ( (delta_reco/reco)^2 + (delta_fc/fc)^2 )
  //
  //    ( Calculation done bin by bin )

  if (!fhFc || !fhRECpt) {
    AliError(" Reconstructed or fc distributions are not defined");
    return;
  }

  Int_t nbins = fhRECpt->GetNbinsX();
  Double_t value = 0., valueDmax= 0., valueDmin= 0.;
  Double_t binwidth = fhRECpt->GetBinWidth(1);
  Double_t xmin = fhRECpt->GetBinCenter(1) - (binwidth/2.) ; 
  Double_t xmax = fhRECpt->GetBinCenter(nbins) + (binwidth/2.) ; 
  
  // declare the output histograms
  TH1D *hYield = new TH1D("hYield","corrected yield (by fc)",nbins,xmin,xmax);
  TH1D *hYieldMax = new TH1D("hYieldMax","max corrected yield (by fc)",nbins,xmin,xmax);
  TH1D *hYieldMin = new TH1D("hYieldMin","min corrected yield (by fc)",nbins,xmin,xmax);
  // and the output TGraphAsymmErrors
  if (fAsymUncertainties & !fgYieldCorr) fgYieldCorr = new TGraphAsymmErrors(nbins);
  
  //
  // Do the calculation
  // 
  for (Int_t ibin=0; ibin<=nbins; ibin++) {

    // calculate the value 
    value = fhRECpt->GetBinContent(ibin) * fhFc->GetBinContent(ibin) ;

    // calculate the value uncertainty
    // Protect against null denominator. If so, define uncertainty as null
    if (fhRECpt->GetBinContent(ibin) && fhRECpt->GetBinContent(ibin)!=0.) {

      if (fAsymUncertainties) {

	if (fhFc->GetBinContent(ibin) && fhFc->GetBinContent(ibin)!=0.) {
	  valueDmax = value * TMath::Sqrt( ( (fhRECpt->GetBinError(ibin)/fhRECpt->GetBinContent(ibin))*(fhRECpt->GetBinError(ibin)/fhRECpt->GetBinContent(ibin)) ) + ( (fgFc->GetErrorYhigh(ibin)/fhFc->GetBinContent(ibin))*(fgFc->GetErrorYhigh(ibin)/fhFc->GetBinContent(ibin)) )  );
	  valueDmin = value * TMath::Sqrt( ( (fhRECpt->GetBinError(ibin)/fhRECpt->GetBinContent(ibin))*(fhRECpt->GetBinError(ibin)/fhRECpt->GetBinContent(ibin)) ) + ( (fgFc->GetErrorYlow(ibin)/fhFc->GetBinContent(ibin))*(fgFc->GetErrorYlow(ibin)/fhFc->GetBinContent(ibin)) ) );
	}
	else { valueDmax = 0.; valueDmin = 0.; }

      }
      else { // Don't consider fc uncertainty in this case [ to be tested!!! ]
	valueDmax = value * (fhRECpt->GetBinError(ibin)/fhRECpt->GetBinContent(ibin)) ; 
	valueDmin = value * (fhRECpt->GetBinError(ibin)/fhRECpt->GetBinContent(ibin)) ;
      }

    }
    else { valueDmax = 0.; valueDmin = 0.; }
    
    // fill in the histograms
    hYield->SetBinContent(ibin,value);
    hYield->SetBinError(ibin,valueDmax);
    hYieldMax->SetBinContent(ibin,value+valueDmax); 
    hYieldMin->SetBinContent(ibin,value-valueDmin);
    if (fAsymUncertainties) {
      Double_t center = hYield->GetBinCenter(ibin);
      fgYieldCorr->SetPoint(ibin,center,value); // i,x,y
      fgYieldCorr->SetPointError(ibin,(binwidth/2.),(binwidth/2.),valueDmin,valueDmax); // i,xl,xh,yl,yh
    }

  }
  
  fhYieldCorr =  hYield;
  fhYieldCorrMax = hYieldMax; 
  fhYieldCorrMin = hYieldMin; 
}

//_________________________________________________________________________________________________________
void AliHFPtSpectrum::CalculateFeedDownCorrectedSpectrumNb(Double_t deltaY, Double_t branchingRatioBintoFinalDecay) {
  //
  // Compute the feed-down corrected spectrum if feed-down correction is done via Nb (bin by bin)
  //    physics =  reco  - (lumi * delta_y * BR_b * eff_trig * eff_b * Nb_th) 
  //
  //    uncertainty: delta_physics = sqrt ( (delta_reco)^2 + (k*delta_lumi/lumi)^2 + 
  //                                        (k*delta_eff_trig/eff_trig)^2 + (k*delta_Nb/Nb)^2 )
  //                    where k = lumi * delta_y * BR_b * eff_trig * eff_b * Nb_th
  //

  Int_t nbins = fhRECpt->GetNbinsX();
  Double_t binwidth = fhRECpt->GetBinWidth(1);
  Double_t value = 0., valueDmax = 0., valueDmin = 0., kfactor = 0.;
  Double_t xmin = fhRECpt->GetBinCenter(1) - (binwidth/2.) ; 
  Double_t xmax = fhRECpt->GetBinCenter(nbins) + (binwidth/2.) ; 
  
  // declare the output histograms
  TH1D *hYield = new TH1D("hYield","corrected yield (by Nb)",nbins,xmin,xmax);
  TH1D *hYieldMax = new TH1D("hYieldMax","max corrected yield (by Nb)",nbins,xmin,xmax);
  TH1D *hYieldMin = new TH1D("hYieldMin","min corrected yield (by Nb)",nbins,xmin,xmax);
  // and the output TGraphAsymmErrors
  if (fAsymUncertainties & !fgYieldCorr) fgYieldCorr = new TGraphAsymmErrors(nbins);

  //
  // Do the calculation
  // 
  for (Int_t ibin=0; ibin<=nbins; ibin++) {
    
    // calculate the value
    value = fhRECpt->GetBinContent(ibin) - (deltaY*branchingRatioBintoFinalDecay*fLuminosity[0]*fTrigEfficiency[0]*fhFeedDownEffpt->GetBinContent(ibin)*fhFeedDownMCpt->GetBinContent(ibin) );

    kfactor = deltaY*branchingRatioBintoFinalDecay*fLuminosity[0]*fTrigEfficiency[0]*fhFeedDownEffpt->GetBinContent(ibin)*fhFeedDownMCpt->GetBinContent(ibin) ;

    // calculate the value uncertainty
    if (fAsymUncertainties) {
      Double_t Nb =  fhFeedDownMCpt->GetBinContent(ibin);
      Double_t NbDmax = fhFeedDownMCptMax->GetBinContent(ibin) - fhFeedDownMCpt->GetBinContent(ibin);
      Double_t NbDmin = fhFeedDownMCpt->GetBinContent(ibin) - fhFeedDownMCptMin->GetBinContent(ibin);
      valueDmax = TMath::Sqrt( ( fhRECpt->GetBinError(ibin)*fhRECpt->GetBinError(ibin) ) +
			       ( (kfactor*fLuminosity[1]/fLuminosity[0])*(kfactor*fLuminosity[1]/fLuminosity[0]) ) +
			       ( (kfactor*fTrigEfficiency[1]/fTrigEfficiency[0])*(kfactor*fTrigEfficiency[1]/fTrigEfficiency[0]) ) +
			       ( (kfactor*NbDmax/Nb)*(kfactor*NbDmax/Nb) ) 	);
      valueDmin =  TMath::Sqrt( ( fhRECpt->GetBinError(ibin)*fhRECpt->GetBinError(ibin) ) +
				( (kfactor*fLuminosity[1]/fLuminosity[0])*(kfactor*fLuminosity[1]/fLuminosity[0]) ) +
				( (kfactor*fTrigEfficiency[1]/fTrigEfficiency[0])*(kfactor*fTrigEfficiency[1]/fTrigEfficiency[0]) ) +
				( (kfactor*NbDmin/Nb)*(kfactor*NbDmin/Nb) ) 	);
    }
    else{ // Don't consider Nb uncertainty in this case [ to be tested!!! ]
      valueDmax =  TMath::Sqrt( ( fhRECpt->GetBinError(ibin)*fhRECpt->GetBinError(ibin) ) +
				( (kfactor*fLuminosity[1]/fLuminosity[0])*(kfactor*fLuminosity[1]/fLuminosity[0]) ) +
				( (kfactor*fTrigEfficiency[1]/fTrigEfficiency[0])*(kfactor*fTrigEfficiency[1]/fTrigEfficiency[0]) ) 	);
      valueDmin =  valueDmax ;
    }
    
    // fill in histograms
    hYield->SetBinContent(ibin,value);
    hYield->SetBinError(ibin,valueDmax);
    hYieldMax->SetBinContent(ibin,value+valueDmax); 
    hYieldMin->SetBinContent(ibin,value-valueDmin);
    if (fAsymUncertainties) {
      Double_t x = hYield->GetBinCenter(ibin);
      fgYieldCorr->SetPoint(ibin,x,value); // i,x,y
      fgYieldCorr->SetPointError(ibin,(binwidth/2.),(binwidth/2.),valueDmin,valueDmax); // i,xl,xh,yl,yh
    }

  }
  
  fhYieldCorr =  hYield;
  fhYieldCorrMax = hYieldMax; 
  fhYieldCorrMin = hYieldMin; 
}


//_________________________________________________________________________________________________________
TH1 * AliHFPtSpectrum::ReweightHisto(TH1 *hToReweight, TH1 *hReference){
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
  TH1 *hReweighted = (TH1*)hToReweight->Clone("hReweighted");
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

  return (TH1*)hReweighted;
}

//_________________________________________________________________________________________________________
TH1 * AliHFPtSpectrum::ReweightRecHisto(TH1 *hRecToReweight, TH1 *hMCToReweight, TH1 *hMCReference){
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
  TH1 *hReweighted = (TH1*)hMCToReweight->Clone("hReweighted");
  hReweighted->Reset();
  TH1 *hRecReweighted = (TH1*)hRecToReweight->Clone("hRecReweighted");
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

  return (TH1*)hRecReweighted;
}

