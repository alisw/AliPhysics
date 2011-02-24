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
#include "AliPHOSRawFitterv2.h"
#include "AliPHOSPulseGenerator.h"

ClassImp(AliPHOSRawFitterv2)

//-----------------------------------------------------------------------------
AliPHOSRawFitterv2::AliPHOSRawFitterv2():
  AliPHOSRawFitterv0(),
  fAlpha(0.1),fBeta(0.035),fMax(0) 
{
  //Default constructor.
}

//-----------------------------------------------------------------------------
AliPHOSRawFitterv2::~AliPHOSRawFitterv2()
{
  //Destructor.
}

//-----------------------------------------------------------------------------
AliPHOSRawFitterv2::AliPHOSRawFitterv2(const AliPHOSRawFitterv2 &phosFitter ):
  AliPHOSRawFitterv0(phosFitter),
  fAlpha(0.1),fBeta(0.035),fMax(0)
{
  //Copy constructor.
}

//-----------------------------------------------------------------------------
AliPHOSRawFitterv2& AliPHOSRawFitterv2::operator = (const AliPHOSRawFitterv2 & /*phosFitter*/)
{
  //Assignment operator.
  return *this;
}

//-----------------------------------------------------------------------------
Bool_t AliPHOSRawFitterv2::Eval(const UShort_t *signal, Int_t sigStart, Int_t sigLength)
{
  //Extract an energy deposited in the crystal,
  //crystal' position (module,column,row),
  //time and gain (high or low).
  //First collects sample, then evaluates it and if it has
  //reasonable shape, fits it with Gamma2 function and extracts 
  //energy and time.

  const Float_t maxEtoFit=5 ; //fit only samples above this energy, accept all samples (with good aRMS) below it
  const Float_t kBaseLine   = 1.0;
  const Int_t   kPreSamples = 10;

  fOverflow = kFALSE ;  
  fEnergy=0 ;
  if (fCaloFlag == 2 || fNBunches > 1) {
    fQuality = 150;
    return kTRUE;
  }
  if(fCaloFlag!=0 && fCaloFlag!=1){//Corrupted sample
    fQuality=200;
    fEnergy=0 ;
    return kTRUE;
  }

  //Evaluate pedestals 
  Float_t pedMean = 0;
  Float_t pedRMS  = 0;
  Int_t   nPed    = 0;
  for (Int_t i=sigLength-kPreSamples; i<sigLength; i++) {
    nPed++;
    pedMean += signal[i];
    pedRMS  += signal[i]*signal[i] ;
  }

  fEnergy = -111;
  fQuality= 999. ;
  Double_t pedestal = 0;

  if (fPedSubtract) {
    if (nPed > 0) {
      fPedestalRMS=(pedRMS - pedMean*pedMean/nPed)/nPed ;
      if(fPedestalRMS > 0.) 
	fPedestalRMS = TMath::Sqrt(fPedestalRMS) ;
      pedestal = (Double_t)(pedMean/nPed); // pedestal subtraction
    }
    else
      return kFALSE;
  }
  else {
    //take pedestals from DB
    pedestal = (Double_t) fAmpOffset ;
  }

  
  //calculate rough quality of the sample and check for overflow
  Int_t    maxBin=0 ;
  Int_t    maxAmp=0 ;
  Int_t    minAmp= signal[0] ;
  Int_t    nMax = 0 ; //number of points in plato
  Double_t aMean =0. ;
  Double_t aRMS  =0. ;
  Double_t wts   =0 ;
  Bool_t falling = kTRUE ; //Bad monotoneusly falling sample
  Bool_t rising = kTRUE ; //Bad monotoneusly riging sample
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
      if(signal[i] <  minAmp)
        minAmp=signal[i] ;
      if(falling && i>0 && signal[i]<signal[i-1])
        falling=kFALSE ;
      if(rising && i>0 && signal[i]>signal[i-1])
        rising=kFALSE ;
    }
  }

  if(rising || falling){//bad "rising" or falling  sample
    fEnergy =    0. ;
    fTime   = 0. ; //-999. ;
    fQuality=  250. ;
    return kTRUE ;
  }
  if(maxAmp-minAmp<3 && maxAmp>7 && sigLength>20){ //bad flat sample
    fEnergy =    0. ;
    fTime   = 0; //-999. ;
    fQuality=  260. ;
    return kTRUE ;
  }

  fEnergy=Double_t(maxAmp)-pedestal ;
  if (fEnergy < kBaseLine) fEnergy = 0;
  fTime = sigStart-sigLength-3;

  //do not test quality of too soft samples
  if (wts > 0) {
    aMean /= wts; 
    aRMS   = aRMS/wts - aMean*aMean;
  }
  if (fEnergy <= maxEtoFit){
    if (aRMS < 2.) //sigle peak
      fQuality = 299. ;
    else
      fQuality =   0. ;
    //Evaluate time of signal arriving
    return kTRUE ;
  }

  //look for plato on the top of sample
  if (fEnergy>500 &&  //this is not fluctuation of soft sample
     nMax>2){ //and there is a plato
    fOverflow = kTRUE ;
  }
  

  //do not fit High Gain samples with overflow
  if(fCaloFlag==1 && fOverflow){
    fQuality = 99. ;
    return kTRUE;

  }

  //----Now fit sample with reasonable shape------
  TArrayD samples(sigLength); // array of sample amplitudes
  TArrayD times(sigLength); // array of sample time stamps
  for (Int_t i=0; i<sigLength; i++) {
    samples.AddAt(signal[i]-pedestal,sigLength-i-1);
    times.AddAt(double(i),i);
  }
      
  if(fMax==0)
    FindMax() ;
  if(!FindAmpT(samples,times)){
    if(AliLog::GetDebugLevel("PHOS","AliPHOSRawFitterv2")>3){
      goto plot ;
    }
    else{
      return kFALSE ;
    }
  }
  fEnergy*=fMax ;
  fTime += sigStart-sigLength-3;


  //Impose cut on quality
//  fQuality/=4. ;
  fQuality/=1.+0.005*fEnergy ;

  //Draw corrupted samples
  if(AliLog::GetDebugLevel("PHOS","AliPHOSRawFitterv2")>3){
    if(fEnergy > 50. ){
    plot:
      printf("Sample par: amp=%f,  t0=%f, Quality=%f \n",fEnergy,fTime,fQuality) ;
      TH1I * h = (TH1I*)gROOT->FindObjectAny("hSamples") ;
      if(!h) h = new TH1I("hSamples","Samples",65,0.,65.) ;
      h->Reset() ;
      for (Int_t i=0; i<sigLength; i++) {
        h->SetBinContent(i+1,float(samples.At(i))) ;
      }
//      TF1 * fffit = new TF1("fffit","[0]+[1]*((x-[2])/[3])^2*exp(2.-2.*(x-[2])/[3])",0.,200.) ;
      TF1 * fffit = new TF1("fffit","[0]*((x-[1])*(x-[1])*exp(-[2]*(x-[1]))+(x-[1])*exp(-[3]*(x-[1])))",0.,60.) ;
      fffit->SetParameters(fEnergy/fMax,fTime-(sigStart-sigLength-3),fAlpha,fBeta) ;
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
//      h->Fit(fffit,"","",0.,51.) ;
      h->Draw() ;
      fffit->Draw("same") ;
/*
      sprintf(title,"mod%d_x%d_z%d_HG_qu%4.1f",fModule,fCellX,fCellZ,fQuality) ;
      TFile fout("samples_bad.root","update") ;
      h->Write(title);
      fout.Close() ;
*/
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
 
      can->Update() ;
      printf("Press <enter> to continue\n") ;
      getchar();


      delete fffit ;
      delete spectrum_1 ;
      delete spectrum_2 ;
    }
  }
  
  return kTRUE;
}
//------------------------------------------------------------------
Bool_t AliPHOSRawFitterv2::FindAmpT(TArrayD samples, TArrayD times){
// makes fit

  const Int_t nMaxIter=50 ;   //Maximal number of iterations
  const Double_t epsdt = 1.e-3 ; //expected precision of t0 calculation

  Double_t dTime=times.At(0)-0.5 ; //Most probable Initial approximation
//printf(" start fit... \n") ;

  Int_t nPoints = samples.GetSize() ;
  Double_t dea=TMath::Exp(-fAlpha) ;
  Double_t deb=TMath::Exp(-fBeta) ;
  Double_t dt=1.,timeOld=dTime,dfOld=0. ; 
  for(Int_t iter=0; iter<nMaxIter; iter++){
    Double_t yy=0.;
    Double_t yf=0. ;
    Double_t ydf=0. ;
    Double_t yddf=0. ;
    Double_t ff=0. ;
    Double_t fdf=0. ;
    Double_t dfdf=0. ;
    Double_t fddf=0. ;
    Int_t nfit=0 ;
    Double_t aexp=TMath::Exp(-fAlpha*(times.At(0)-1.-dTime)) ;
    Double_t bexp=TMath::Exp(-fBeta*(times.At(0)-1.-dTime)) ;
    for(Int_t i=0; i<nPoints; i++){
      Double_t t= times.At(i)-dTime ;
      aexp*=dea ;
      bexp*=deb ;
      if(t<0.) continue ;
      Double_t y=samples.At(i) ;
      if(y<=fAmpThreshold)
        continue ;
      nfit++ ;
      Double_t at=fAlpha*t ;
      Double_t bt = fBeta*t ;
      Double_t phi=t*(t*aexp+bexp) ;
      Double_t dphi=t*aexp*(2.-at)+(1.-bt)*bexp ;
      Double_t ddphi=aexp*(2.-4.*at+at*at)+bexp*fBeta*(bt-2.) ;
      yy+=y*y ;
      yf+=y*phi ;
      ydf+=y*dphi ;
      yddf+=y*ddphi ;
      ff+=phi*phi ;
      fdf+=phi*dphi ;
      dfdf+=dphi*dphi ;
      fddf+=phi*ddphi ;
    }

    if(ff==0.||nfit==0. ){
      fQuality=199 ;
      return kFALSE ;
    }
    Double_t f=ydf*ff-yf*fdf ;     //d(chi2)/dt
    Double_t df=yf*(dfdf+fddf)-yddf*ff-ydf*fdf;
    if(df<=0.){ //we are too far from the root. In the wicinity of root df>0
      if(iter!=0 && dfOld>0.){//If at previous step df was OK, just reduce step size
        dt*=0.5 ;
        dTime=timeOld+dt ;  
        continue ;
      }
      if(f<0){ //f<0 => dTime is too small and we still do not know root region
        dTime+=2. ;
        continue ;
      }
      else{ //dTime is too large, we are beyond the root region
        dTime-=2. ;
        continue ;
      }
    }
    dt=-f/df ; 
    if(TMath::Abs(dt)<epsdt){
      fQuality=(yy-yf*yf/ff)/nfit ;
      fEnergy=yf/ff ;  //ff!=0 already tested
      fTime=dTime ;
      return kTRUE ;
    }
    //In some cases time steps are huge (derivative ~0)
    if(dt>10.) dt=10. ;   //restrict step size
    if(dt<-10.) dt=-5.3 ; //this restriction should be asimmetric to avoid jumping from one point to another
    timeOld=dTime ;  //remember current position for the case
    dfOld=df ;       //of reduction of dt step size
    dTime+=dt ;

    if(dTime>100. || dTime<-30.){ //this is corrupted sample, do not spend time improving accuracy.
      fQuality=(yy-yf*yf/ff)/nfit ;
      fEnergy=yf/ff ;  //ff!=0 already tested
      fTime=dTime ;
      return kFALSE ;
    }

  }
  //failed to find a root, too many iterations
  fQuality=99.;
  fEnergy=0 ; 
  return kFALSE ;
}
//_________________________________________
void AliPHOSRawFitterv2::FindMax(){
  //Finds maxumum of currecnt parameterization
  Double_t t=2./fAlpha ;
  fMax = t*t*TMath::Exp(-fAlpha*t)+t*TMath::Exp(-fBeta*t) ;
  Double_t dt=15 ;
  while(dt>0.01){
     Double_t dfdt=(2.*t-fAlpha*t*t)*TMath::Exp(-fAlpha*t)+(1.-fBeta*t)*TMath::Exp(-fBeta*t) ;
     if(dfdt>0.)
        t+=dt ;
     else
       t-=dt ;
     Double_t maxNew = t*t*TMath::Exp(-fAlpha*t)+t*TMath::Exp(-fBeta*t) ;
     if(maxNew>fMax)
        fMax=maxNew ;
     else{
       dt/=2 ;
       if(dfdt<0.)
         t+=dt ;
       else
        t-=dt ;
     }
  }   
}
 
