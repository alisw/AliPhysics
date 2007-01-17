/**************************************************************************
 * Copyright(c) 2007, ALICE Experiment at CERN, All rights reserved.      *
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

// The class which simulates the pulse shape from the PHOS FEE shaper,
// make sampled amplitudes, digitize them.
// Use case:
//   AliPHOSPulseGenerator *pulse = new AliPHOSPulseGenerator(energy,time);
//   Int_t *adcHG = new Int_t[pulse->GetRawFormatTimeBins()];
//   Int_t *adcLG= new Int_t[pulse->GetRawFormatTimeBins()];
//   pulse->MakeSamples();
//   pulse->GetSamples(adcHG, adcHG) ; 
//   pulse->Print();
//   pulse->Draw();
//
// Author: Yuri Kharlov, after Yves Schutz and Per Thomas Hille

// --- ROOT system ---
#include "TMath.h" 
#include "TF1.h" 
#include "TGraph.h" 
#include "TCanvas.h" 

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliPHOSPulseGenerator.h"

// --- Standard library ---
#include <cmath>
#include <iostream>

using std::cout;
using std::endl; 

ClassImp(AliPHOSPulseGenerator) 

Double_t AliPHOSPulseGenerator::fgCapa        = 1.;        // 1pF 
Int_t    AliPHOSPulseGenerator::fgOrder       = 2 ;        // order of the Gamma function
Double_t AliPHOSPulseGenerator::fgTimeMax     = 2.56E-5 ;  // each sample is over 100 ns fTimeMax/fTimeBins
Double_t AliPHOSPulseGenerator::fgTimePeak    = 4.1E-6 ;   // 4 micro seconds
Double_t AliPHOSPulseGenerator::fgTimeTrigger = 100E-9 ;   // 100ns, just for a reference
Double_t AliPHOSPulseGenerator::fgHighCharge  = 8.2;       // adjusted for a high gain range of 5.12 GeV (10 bits)
Double_t AliPHOSPulseGenerator::fgHighGain    = 6.64;
Double_t AliPHOSPulseGenerator::fgHighLowGainFactor = 16.; // adjusted for a low gain range of 82 GeV (10 bits) 

//-----------------------------------------------------------------------------
AliPHOSPulseGenerator::AliPHOSPulseGenerator(Double_t a, Double_t t0)
  : TObject(), fAmplitude(a), fTZero(t0), fDataHG(0), fDataLG(0), fDigitize(kTRUE)
{
  // Contruct a pulsegenrator object and initializes all necessary parameters
  // @param a digit amplitude in GeV
  // @param t0 time delay in nanoseconds of signal relative the first sample. 
  // This value should be between 0 and Ts, where Ts is the sample interval

  fDataHG = new Double_t[fkTimeBins];
  fDataLG = new Double_t[fkTimeBins];
}

//-----------------------------------------------------------------------------
AliPHOSPulseGenerator::AliPHOSPulseGenerator(const AliPHOSPulseGenerator & pulse)
  : TObject(), fAmplitude(pulse.fAmplitude), fTZero(pulse.fTZero), fDataHG(0), fDataLG(0), fDigitize(kTRUE)
{
  fDataHG = new Double_t[pulse.fkTimeBins];
  fDataLG = new Double_t[pulse.fkTimeBins];
  for (Int_t i=0; i<pulse.fkTimeBins; i++) {
    fDataHG[i] = pulse.fDataHG[i];
    fDataLG[i] = pulse.fDataHG[i];
  }
}

//-----------------------------------------------------------------------------
AliPHOSPulseGenerator::~AliPHOSPulseGenerator()
{
  // Destructor: delete arrays of samples

  delete [] fDataHG;
  fDataHG=0;
  delete [] fDataLG;
  fDataLG=0;
}

//-----------------------------------------------------------------------------
void AliPHOSPulseGenerator::AddBaseline(Double_t baselineLevel)
{
  // Adds a baseline offset to the signal
  // @param baselineLevel The basline level to add
  for (Int_t i=0; i<fkTimeBins; i++) {
    fDataHG[i] += baselineLevel;
    fDataLG[i] += baselineLevel;
  }
}

//-----------------------------------------------------------------------------
void AliPHOSPulseGenerator::AddNoise(Double_t * /* sigma */)
{
  // Adds Gaussian white noise to the sample array given by *dataPtr.
  // @param sigma the noise amplitude in entities of ADC levels  
  
  AliError("not implemented yet");
}


//-----------------------------------------------------------------------------
void AliPHOSPulseGenerator::AddNoise(Double_t * /* sigma */, Double_t /* cutoff */)
{
  //Adds correlated Gaussian noise with cutof frequency "cutoff"
  // @param sigma noise amplitude in entities of ADC levels
  // @param -30DB cutoff frequency of the noise in entities of sampling frequency

  AliError("not implemented yet");
}

//-----------------------------------------------------------------------------
void AliPHOSPulseGenerator::AddPretriggerSamples(Int_t nPresamples)
{
  // Adds pretrigger samples to the sample arrays and replace them
  // with concatinated and truncated arrays

  Double_t *tmpDataHG = new Double_t[fkTimeBins];
  Double_t *tmpDataLG = new Double_t[fkTimeBins];
  Int_t i;
  for (i=0; i<fkTimeBins; i++) {
    tmpDataHG[i] = fDataHG[i];
    tmpDataLG[i] = fDataLG[i];
  }
  for (i=0; i<fkTimeBins; i++) {
    if (i<nPresamples) {
      tmpDataHG[i] = 0.;
      tmpDataLG[i] = 0.;
    }
    else {
      tmpDataHG[i] = fDataHG[i-nPresamples];
      tmpDataLG[i] = fDataLG[i-nPresamples];
    }
  }
  delete [] tmpDataHG;
  delete [] tmpDataLG;
}

//-----------------------------------------------------------------------------
void AliPHOSPulseGenerator::Digitize()
{
  // Emulates ADC: rounds down to nearest integer value all amplitudes
  for (Int_t i=0; i<fkTimeBins; i++) {
    fDataHG[i] = (Double_t) ((Int_t)(fDataHG[i] + 0.5));
    fDataLG[i] = (Double_t) ((Int_t)(fDataLG[i] + 0.5));
  }
}

//-----------------------------------------------------------------------------
Double_t AliPHOSPulseGenerator::RawResponseFunction(Double_t *x, Double_t *par) 
{
  // Shape of the electronics raw reponse:
  // It is a semi-gaussian, 2nd order Gamma function of the general form
  // v(t) = n**n * Q * A**n / C *(t/tp)**n * exp(-n * t/tp) with 
  // tp : peaking time par[0]
  // n  : order of the function
  // C  : integrating capacitor in the preamplifier
  // A  : open loop gain of the preamplifier
  // Q  : the total APD charge to be measured Q = C * energy
  
  Double_t signal ;
  Double_t xx = x[0] - ( fgTimeTrigger + par[3] ) ; 

  if (xx < 0 || xx > fgTimeMax) 
    signal = 0. ;  
  else {
    Double_t fac = par[0] * TMath::Power(fgOrder, fgOrder) * TMath::Power(par[1], fgOrder)/fgCapa ; 
    signal = fac * par[2] * TMath::Power(xx/fgTimePeak, fgOrder) * TMath::Exp(-fgOrder*(xx/fgTimePeak)) ;
  }
  return signal ;  
}

//__________________________________________________________________
Double_t AliPHOSPulseGenerator::RawResponseFunctionMax(Double_t charge, Double_t gain) 
{
  // Maximum value of the shaper response function
  return ( charge * TMath::Power(fgOrder, fgOrder) * TMath::Power(gain, fgOrder) 
	   / ( fgCapa * TMath::Exp(fgOrder) ) );  
}

//__________________________________________________________________
Bool_t AliPHOSPulseGenerator::MakeSamples()
{
  // for a start time fTZero and an amplitude fAmplitude given by digit, 
  // calculates the raw sampled response AliPHOSPulseGenerator::RawResponseFunction

  const Int_t kRawSignalOverflow = 0x3FF ; 
  Bool_t lowGain = kFALSE ; 

  TF1 signalF("signal", RawResponseFunction, 0, GetRawFormatTimeMax(), 4);

  for (Int_t iTime = 0; iTime < GetRawFormatTimeBins(); iTime++) {
    signalF.SetParameter(0, fgHighCharge) ; 
    signalF.SetParameter(1, fgHighGain) ; 
    signalF.SetParameter(2, fAmplitude) ; 
    signalF.SetParameter(3, fTZero) ; 
    Double_t time = iTime * GetRawFormatTimeMax() / GetRawFormatTimeBins() ;
    Double_t signal = signalF.Eval(time) ;     
    if ( static_cast<Int_t>(signal+0.5) > kRawSignalOverflow ){  // larger than 10 bits 
      signal = kRawSignalOverflow ;
      lowGain = kTRUE ; 
    }
    fDataHG[iTime] = signal;

    signalF.SetParameter(0, GetRawFormatLowCharge() ) ;     
    signalF.SetParameter(1, GetRawFormatLowGain() ) ; 
    signal = signalF.Eval(time) ;  
    if ( static_cast<Int_t>(signal+0.5) > kRawSignalOverflow)  // larger than 10 bits 
      signal = kRawSignalOverflow ;
    fDataLG[iTime] = signal;

  }
  return lowGain ; 
}

//__________________________________________________________________
void AliPHOSPulseGenerator::GetSamples(Int_t *adcHG, Int_t *adcLG) const
{
  // Return integer sample arrays adcHG and adcLG
  for (Int_t iTime = 0; iTime < GetRawFormatTimeBins(); iTime++) {
    adcHG[iTime] = static_cast<Int_t>(fDataHG[iTime]) ;
    adcLG[iTime] = static_cast<Int_t>(fDataLG[iTime]) ;
  }
}

//__________________________________________________________________
void AliPHOSPulseGenerator::Print(Option_t*) const
{
  // Prints sampled amplitudes to stdout
  Int_t i;
  cout << "High gain: ";
  for (i=0; i<fkTimeBins; i++)
    cout << (Int_t)fDataHG[i] << " ";
  cout << endl;

  cout << "Low  gain: ";
  for (i=0; i<fkTimeBins; i++)
    cout << (Int_t)fDataLG[i] << " ";
  cout << endl;
}

//__________________________________________________________________
void AliPHOSPulseGenerator::Draw(Option_t*)
{
  // Draw graphs with high and low gain samples

  Double_t *time = new Double_t[fkTimeBins];
  for (Int_t iTime = 0; iTime < GetRawFormatTimeBins(); iTime++) {
    time[iTime] = iTime * GetRawFormatTimeMax() / GetRawFormatTimeBins() ;
  }
  Int_t nPoints = fkTimeBins;
  TGraph *graphHG = new TGraph(nPoints,time,fDataHG);
  TGraph *graphLG = new TGraph(nPoints,time,fDataLG);
  graphHG->SetMarkerStyle(20);
  graphLG->SetMarkerStyle(24);
  graphHG->SetMarkerSize(0.3);
  graphLG->SetMarkerSize(0.3);
  graphHG->SetTitle("High gain samples");
  graphLG->SetTitle("Low gain samples");

  TCanvas *c1 = new TCanvas("c1","Raw ALTRO samples",10,10,700,500);
  c1->Divide(2,1);
  c1->cd(1);
  gPad->SetLeftMargin(0.15);
  graphHG->Draw("AP");
  graphHG->GetHistogram()->SetTitleOffset(1.6,"Y");
  graphHG->GetHistogram()->SetXTitle("time, #musec");
  graphHG->GetHistogram()->SetYTitle("Amplitude, ADC counts");
  c1->cd(2);
  gPad->SetLeftMargin(0.15);
  graphLG->Draw("AP");
  graphLG->GetHistogram()->SetTitleOffset(1.6,"Y");
  graphLG->GetHistogram()->SetXTitle("time, #musec");
  graphLG->GetHistogram()->SetYTitle("Amplitude, ADC counts");
  c1->Update();
}
