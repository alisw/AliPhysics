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

/* $Id: $ */

// This class extracts the signal parameters (energy, time, quality)
// from ALTRO samples. Energy is in ADC counts, time is in time bin units.
// Class uses FastFitting algorithm to fit sample and extract time and Amplitude 
// and evaluate sample quality = (chi^2/NDF)/some parameterization providing 
// efficiency close to 100%
// 
// Typical use case:
//     AliPHOSRawFitter *fitter=new AliPHOSRawFitter();
//     fitter->SetChannelGeo(module,cellX,cellZ,caloFlag);
//     fitter->SetCalibData(fgCalibData) ;
//     fitter->Eval(sig,sigStart,sigLength);
//     Double_t amplitude = fitter.GetEnergy();
//     Double_t time      = fitter.GetTime();
//     Bool_t   isLowGain = fitter.GetCaloFlag()==0;

// Author: Dmitri Peressounko (after A.Pavlinov - see RAW/AliCaloFastAltroFitv0.cxx)

// --- ROOT system ---
#include "TArrayI.h"
#include "TList.h"
#include "TMath.h"
#include "TH1I.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TROOT.h"

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliPHOSCalibData.h"
#include "AliPHOSRawFitterv3.h"
#include "AliPHOSPulseGenerator.h"

ClassImp(AliPHOSRawFitterv3)

//-----------------------------------------------------------------------------
AliPHOSRawFitterv3::AliPHOSRawFitterv3():
  AliPHOSRawFitterv0()
{
  //Default constructor.
}

//-----------------------------------------------------------------------------
AliPHOSRawFitterv3::~AliPHOSRawFitterv3()
{
  //Destructor.
}

//-----------------------------------------------------------------------------
AliPHOSRawFitterv3::AliPHOSRawFitterv3(const AliPHOSRawFitterv3 &phosFitter ):
  AliPHOSRawFitterv0(phosFitter) 
{
  //Copy constructor.
}

//-----------------------------------------------------------------------------
AliPHOSRawFitterv3& AliPHOSRawFitterv3::operator = (const AliPHOSRawFitterv3 & /*phosFitter*/)
{
  //Assignment operator.
  return *this;
}

//-----------------------------------------------------------------------------
Bool_t AliPHOSRawFitterv3::Eval(const UShort_t *signal, Int_t sigStart, Int_t sigLength)
{
  //Extract an energy deposited in the crystal,
  //crystal' position (module,column,row),
  //time and gain (high or low).
  //First collects sample, then evaluates it and if it has
  //reasonable shape, fits it with Gamma2 function and extracts 
  //energy and time.

  if (fCaloFlag == 2 || fNBunches > 1) {
    fQuality = 1000;
    return kTRUE;
  }
  if(fCaloFlag!=0 && fCaloFlag!=1){//Corrupted sample
    fQuality=2000;
    fEnergy=0 ;
    return kTRUE;
  }

  Float_t pedMean = 0;
  Float_t pedRMS  = 0;
  Int_t   nPed    = 0;
  const Float_t kBaseLine   = 1.0;
  const Int_t   kPreSamples = 10;
  

  //We tryed individual taus for each channel, but
  //this approach seems to be unstable. Much better results are obtaned with
  //fixed decay time for all channels.
  const Double_t tau=22.18 ;

  TArrayD samples(sigLength); // array of sample amplitudes
  TArrayD times(sigLength); // array of sample time stamps
  for (Int_t i=sigLength-kPreSamples; i<sigLength; i++) {
    nPed++;
    pedMean += signal[i];
    pedRMS  += signal[i]*signal[i] ;
  }

  fEnergy = -111;
  fQuality= 999. ;
  const Float_t maxEtoFit=5 ; //fit only samples above this energy, accept all samples (with good aRMS) below it
  Double_t pedestal = 0;

  if (fPedSubtract) {
    if (nPed > 0) {
      fPedestalRMS=(pedRMS - pedMean*pedMean/nPed)/nPed ;
      if(fPedestalRMS > 0.) 
	fPedestalRMS = TMath::Sqrt(fPedestalRMS) ;
      pedestal = (Double_t)(pedMean/nPed); // pedestal subtraction
      fEnergy -= pedestal; // pedestal subtraction
    }
    else
      return kFALSE;
  }
  else {
    //take pedestals from DB
    pedestal = (Double_t) fAmpOffset ;
    if (fCalibData) {
      Float_t truePed       = fCalibData->GetADCpedestalEmc(fModule, fCellZ, fCellX) ;
      Int_t   altroSettings = fCalibData->GetAltroOffsetEmc(fModule, fCellZ, fCellX) ;
      pedestal += truePed - altroSettings ;
    }
    else{
//      AliWarning(Form("Can not read data from OCDB")) ;
    }
    fEnergy-=pedestal ;
  }

  if (fEnergy < kBaseLine) fEnergy = 0;
  //Evaluate time
   fTime = sigStart-sigLength-3;
  for (Int_t i=0; i<sigLength; i++) {
    samples.AddAt(signal[i]-pedestal,sigLength-i-1);
    times.AddAt(i/tau ,i);
  }
  
  //calculate time and energy
  Int_t    maxBin=0 ;
  Int_t    maxAmp=0 ;
  Int_t    nMax = 0 ; //number of points in plato
  Double_t aMean =0. ;
  Double_t aRMS  =0. ;
  Double_t wts   =0 ;
  for (Int_t i=0; i<sigLength; i++){
    if(signal[i] > pedestal){
      Double_t de = signal[i] - pedestal ;
      if(de > 1.) {
	aMean += de*i ;
	aRMS  += de*i*i ;
	wts   += de; 
      }
      if(signal[i] >  maxAmp){
        maxAmp = signal[i]; 
        nMax=0;
	maxBin = i ;
      }
      if(signal[i] == maxAmp){
        nMax++;
      }
    }
  }

  if (maxBin==sigLength-1){//bad "rising" sample
    fEnergy =    0. ;
    fTime   = -999. ;
    fQuality=  999. ;
    return kTRUE ;
  }

  fEnergy=Double_t(maxAmp)-pedestal ;
  fOverflow =0 ;  //look for plato on the top of sample
  if (fEnergy>500 &&  //this is not fluctuation of soft sample
     nMax>2){ //and there is a plato
    fOverflow = kTRUE ;
  }
  
  if (wts > 0) {
    aMean /= wts; 
    aRMS   = aRMS/wts - aMean*aMean;
  }

  //do not take too small energies
  if (fEnergy < kBaseLine) 
    fEnergy = 0;
  
  //do not test quality of too soft samples
  if (fEnergy < maxEtoFit){
    if (aRMS < 2.) //sigle peak
      fQuality = 999. ;
    else
      fQuality =   0. ;
    //Evaluate time of signal arriving
    return kTRUE ;
  }
      
  // if sample has reasonable mean and RMS, try to fit it with gamma2
  //This method can not analyse overflow samples
  if(fOverflow){
    fQuality = 99. ;
    return kTRUE ;
  }
  // First estimate t0
  Double_t a=0,b=0,c=0 ;
  Int_t minI=0 ;
  if (fPedSubtract) 
    minI=kPreSamples ;
  for(Int_t i=minI; i<sigLength; i++){
    if(samples.At(i)<=0.)
      continue ;
    Double_t t= times.At(i) ;
    Double_t f02 = TMath::Exp(-2.*t);
    Double_t f12 = t*f02;
    Double_t f22 = t*f12;
    // Derivatives
    Double_t f02d = -2.*f02;
    Double_t f12d = f02 - 2.*f12;
    Double_t f22d = 2.*(f12 - f22);
    a += f02d * samples.At(i) ;
    b -= f12d * samples.At(i) ;
    c += f22d * samples.At(i) ;
  }
  
  //Find roots
  if(a==0.){
    if(b==0.){ //no roots
      fQuality = 2000. ;
      if(AliLog::GetDebugLevel("PHOS","AliPHOSRawFitterv3")>3){
        printf(" a=%f, b=%f, c=%f \n",a,b,c) ;
        goto plot ;
      }
      return kTRUE ;
    }
    Double_t t1=-c/b ;
    Double_t amp=0.,den=0.; ;
    for(Int_t i=minI; i<sigLength; i++){
      if(samples.At(i)<=0.)
        continue ;
      Double_t dt = times.At(i) - t1;
      Double_t f = dt*dt*TMath::Exp(-2.*dt);
      amp += f*samples.At(i);
      den += f*f;
    }
    if(den>0.0) amp /= den;
    // chi2 calculation
    fQuality=0.;
    for(Int_t i=minI; i<sigLength; i++){
      if(samples.At(i)<=0.)
        continue ;
      Double_t t = times.At(i)- t1;
      Double_t dy = samples.At(i)- amp*t*t*TMath::Exp(-2.*t) ;
      fQuality += dy*dy;
    }
    fTime+=t1*tau ;
    fEnergy = amp*TMath::Exp(-2.);
    fQuality/= sigLength ; //If we have overflow the number of actually fitted points is smaller, but chi2 in this case is not important.
  }
  else{
    Double_t det = b*b - a*c;
    if(det>=1.e-6 && det<0.0) {
      det = 0.0; //rounding error
    }
    if(det<0.){ //Problem
      fQuality = 1500. ;
      if(AliLog::GetDebugLevel("PHOS","AliPHOSRawFitterv3")>3){
        printf(" det=%e \n",det) ;
        goto plot ;
      }
      return kTRUE ;
    }

    det = TMath::Sqrt(det);
    Double_t t1 = (-b + det) / a;
//    Double_t t2 = (-b - det) / a; //second root is wrong one
    Double_t amp1=0., den1=0. ;
    for(Int_t i=minI; i<sigLength; i++){
      if(samples.At(i)<=0.)
        continue ;
      Double_t dt1 = times.At(i) - t1;
      Double_t f01 = dt1*dt1*TMath::Exp(-2.*dt1);
      amp1 += f01*samples.At(i);
      den1 += f01*f01;
    }
    if(den1>0.0) amp1 /= den1;
    Double_t chi1=0.; // chi2=0. ;
    for(Int_t i=minI; i<sigLength; i++){
      if(samples.At(i)<=0.)
        continue ;
      Double_t dt1 = times.At(i)- t1;
      Double_t dy1 = samples.At(i)- amp1*dt1*dt1*TMath::Exp(-2.*dt1) ;
      chi1 += dy1*dy1;
    }
    fEnergy=amp1*TMath::Exp(-2.) ; ; 
    fTime+=t1*tau ;
    fQuality=chi1/sigLength ;
  } 

  //Impose cut on quality
  fQuality/=2.+0.004*fEnergy*fEnergy ;

  //Draw corrupted samples
  if(AliLog::GetDebugLevel("PHOS","AliPHOSRawFitterv3")>3){
    plot:
    if(fEnergy > 30. && fQuality >1. && !fOverflow ){ //Draw only bad samples
//    if(!fOverflow ){ //Draw only bad samples
      printf("Sample par: amp=%f,  t0=%f, Quality=%f \n",fEnergy,fTime,fQuality) ;
      TH1I * h = (TH1I*)gROOT->FindObjectAny("hSamples") ;
      if(!h) h = new TH1I("hSamples","Samples",65,0.,65.) ;
      h->Reset() ;
      for (Int_t i=0; i<sigLength; i++) {
        h->SetBinContent(i+1,samples.At(i)+pedestal) ;
      }
      TF1 * fffit = new TF1("fffit","[0]+[1]*((x-[2])/[3])^2*exp(2.-2.*(x-[2])/[3])",0.,200.) ;
      fffit->SetParameters(pedestal,fEnergy,fTime,tau) ;
      fffit->SetLineColor(2) ;
      TCanvas * can = (TCanvas*)gROOT->FindObjectAny("cSamples") ;
      if(!can){
        can = new TCanvas("cSamples","cSamples",10,10,600,600) ;
        can->SetFillColor(0) ;
        can->SetFillStyle(0) ;
        can->Range(0,0,1,1);
        can->SetBorderSize(0);
      }
      can->cd() ;
  
      TPad * spectrum_1 = new TPad("spectrum_1", "spectrum_1",0.001,0.32,0.99,0.99);
      spectrum_1->Draw();
      spectrum_1->cd();
      spectrum_1->Range(0,0,1,1);
      spectrum_1->SetFillColor(0);
      spectrum_1->SetFillStyle(4000);
      spectrum_1->SetBorderSize(1);
      spectrum_1->SetBottomMargin(0.012);
      spectrum_1->SetTopMargin(0.03);
      spectrum_1->SetLeftMargin(0.10);
      spectrum_1->SetRightMargin(0.05);

      char title[155] ;
      snprintf(title,155,"Sample, mod=%d, x=%d, z=%d, Quality=%5.1f",fModule,fCellX,fCellZ,fQuality) ;
      h->SetTitle(title) ;
      h->Draw() ;
      fffit->Draw("same") ;

      snprintf(title,155,"mod%d_x%d_z%d_HG_qu%4.1f",fModule,fCellX,fCellZ,fQuality) ;
      TFile fout("samples_bad.root","update") ;
      h->Write(title);
      fout.Close() ;

      can->cd() ;
      TPad *spectrum_2 = new TPad("spectrum_2", "spectrum_2",0.001,0.01,0.99,0.32);
      spectrum_2->SetFillColor(0) ;
      spectrum_2->SetFillStyle(0) ;
      spectrum_2->SetGridy() ;
      spectrum_2->Draw();
      spectrum_2->Range(0,0,1,1);
      spectrum_2->SetFillColor(0);
      spectrum_2->SetBorderSize(1);
      spectrum_2->SetTopMargin(0.01);
      spectrum_2->SetBottomMargin(0.25);
      spectrum_2->SetLeftMargin(0.10);
      spectrum_2->SetRightMargin(0.05);
      spectrum_2->cd() ;

      TH1I * hd = (TH1I*)gROOT->FindObjectAny("hSamplesDif") ;
      if(!hd) hd = new TH1I("hd","Samples",65,0.,65.) ;
      hd->Reset() ;
      for (Int_t i=0; i<sigLength; i++) {
        hd->SetBinContent(i+1,TMath::Max(-1023.,TMath::Min(1023.,samples.At(i)+pedestal-fffit->Eval(i)))) ;
      }
      hd->Draw();
/* 
      can->Update() ;
      printf("Press <enter> to continue\n") ;
      getchar();
*/

      delete fffit ;
      delete spectrum_1 ;
      delete spectrum_2 ;
    }
  }
  
  return kTRUE;
}
