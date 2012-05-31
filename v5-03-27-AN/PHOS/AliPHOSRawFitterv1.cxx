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
// A fitting algorithm evaluates the energy and the time from Minuit minimization
// 
// Typical use case:
//     AliPHOSRawFitter *fitter=new AliPHOSRawFitter();
//     fitter->SetChannelGeo(module,cellX,cellZ,caloFlag);
//     fitter->SetCalibData(fgCalibData) ;
//     fitter->Eval(sig,sigStart,sigLength);
//     Double_t amplitude = fitter.GetEnergy();
//     Double_t time      = fitter.GetTime();
//     Bool_t   isLowGain = fitter.GetCaloFlag()==0;

// Author: Dmitri Peressounko (Oct.2008)
// Modified: Yuri Kharlov (Jul.2009)

// --- ROOT system ---
#include "TArrayI.h"
#include "TList.h"
#include "TMath.h"
#include "TMinuit.h"

// --- AliRoot header files ---
#include "AliLog.h"
#include "AliPHOSCalibData.h"
#include "AliPHOSRawFitterv1.h"
#include "AliPHOSPulseGenerator.h"

ClassImp(AliPHOSRawFitterv1)

//-----------------------------------------------------------------------------
AliPHOSRawFitterv1::AliPHOSRawFitterv1():
  AliPHOSRawFitterv0(),
  fSampleParamsLow(0x0),
  fSampleParamsHigh(0x0),
  fToFit(0x0)
{
  //Default constructor.
  if(!gMinuit) 
    gMinuit = new TMinuit(100);
  fSampleParamsHigh =new TArrayD(7) ;
  fSampleParamsHigh->AddAt(2.174,0) ;
  fSampleParamsHigh->AddAt(0.106,1) ;
  fSampleParamsHigh->AddAt(0.173,2) ;
  fSampleParamsHigh->AddAt(0.06106,3) ;
  //last two parameters are pedestal and overflow
  fSampleParamsLow=new TArrayD(7) ;
  fSampleParamsLow->AddAt(2.456,0) ;
  fSampleParamsLow->AddAt(0.137,1) ;
  fSampleParamsLow->AddAt(2.276,2) ;
  fSampleParamsLow->AddAt(0.08246,3) ;
  fToFit = new TList() ;
}

//-----------------------------------------------------------------------------
AliPHOSRawFitterv1::~AliPHOSRawFitterv1()
{
  //Destructor.
  //Destructor
  if(fSampleParamsLow){
    delete fSampleParamsLow ; 
    fSampleParamsLow=0 ;
  }
  if(fSampleParamsHigh){
    delete fSampleParamsHigh ;
    fSampleParamsHigh=0;
  }
  if(fToFit){
    delete fToFit ;
    fToFit=0 ;
  }
}

//-----------------------------------------------------------------------------
AliPHOSRawFitterv1::AliPHOSRawFitterv1(const AliPHOSRawFitterv1 &phosFitter ):
  AliPHOSRawFitterv0(phosFitter),
  fSampleParamsLow(0x0),
  fSampleParamsHigh(0x0),
  fToFit(0x0)
{
  //Copy constructor.
  fToFit = new TList() ;
  fSampleParamsLow =new TArrayD(*(phosFitter.fSampleParamsLow)) ;
  fSampleParamsHigh=new TArrayD(*(phosFitter.fSampleParamsHigh)) ;
}

//-----------------------------------------------------------------------------
AliPHOSRawFitterv1& AliPHOSRawFitterv1::operator = (const AliPHOSRawFitterv1 &phosFitter)
{
  //Assignment operator.
  if(this != &phosFitter) {
    fToFit = new TList() ;
    if(fSampleParamsLow){
      fSampleParamsLow = phosFitter.fSampleParamsLow ;
      fSampleParamsHigh= phosFitter.fSampleParamsHigh ;
    }
    else{
      fSampleParamsLow =new TArrayD(*(phosFitter.fSampleParamsLow)) ; 
      fSampleParamsHigh=new TArrayD(*(phosFitter.fSampleParamsHigh)) ;
    }
  }
  return *this;
}

//-----------------------------------------------------------------------------
Bool_t AliPHOSRawFitterv1::Eval(const UShort_t *signal, Int_t sigStart, Int_t sigLength)
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

  Float_t pedMean = 0;
  Float_t pedRMS  = 0;
  Int_t   nPed    = 0;
  const Float_t kBaseLine   = 1.0;
  const Int_t   kPreSamples = 10;
  
  TArrayI *fSamples = new TArrayI(sigLength); // array of sample amplitudes
  TArrayI *fTimes   = new TArrayI(sigLength); // array of sample time stamps
  for (Int_t i=0; i<sigLength; i++) {
    if (i<kPreSamples) {
      nPed++;
      pedMean += signal[i];
      pedRMS  += signal[i]*signal[i] ;
    }
    fSamples->AddAt(signal[i],sigLength-i-1);
    fTimes  ->AddAt(i ,i);
  }

  fEnergy = -111;
  fQuality= 999. ;
  const Float_t sampleMaxHG=102.332 ;  //maximal height of HG sample with given parameterization
  const Float_t sampleMaxLG=277.196 ;  //maximal height of LG sample with given parameterization
  const Float_t maxEtoFit=5 ; //fit only samples above this energy, accept all samples (with good aRMS) below it
  Double_t pedestal = 0;

  if (fPedSubtract) {
    if (nPed > 0) {
      fPedestalRMS=(pedRMS - pedMean*pedMean/nPed)/nPed ;
      if(fPedestalRMS > 0.) 
	fPedestalRMS = TMath::Sqrt(fPedestalRMS) ;
      fEnergy -= (Double_t)(pedMean/nPed); // pedestal subtraction
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
      AliWarning(Form("Can not read data from OCDB")) ;
    }
    fEnergy-=pedestal ;
  }

  if (fEnergy < kBaseLine) fEnergy = 0;
  //Evaluate time
  Int_t iStart = 0;
  while(iStart<sigLength && fSamples->At(iStart)-pedestal <kBaseLine) iStart++ ;
  fTime = sigStart-sigLength+iStart; 
  
  //calculate time and energy
  Int_t    maxBin=0 ;
  Int_t    maxAmp=0 ;
  Double_t aMean =0. ;
  Double_t aRMS  =0. ;
  Double_t wts   =0 ;
  Int_t    tStart=0 ;

  for (Int_t i=0; i<sigLength; i++){
    if(signal[i] > pedestal){
      Double_t de = signal[i] - pedestal ;
      if(de > 1.) {
	aMean += de*i ;
	aRMS  += de*i*i ;
	wts   += de; 
      }
      if(de > 2 && tStart==0) 
	tStart = i ;
      if(maxAmp < signal[i]){
	maxBin = i ;
	maxAmp = signal[i] ;
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
     maxBin<sigLength-1 && fSamples->At(maxBin+1)==maxAmp){ //and there is a plato
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
    fTime = tStart;
    if (aRMS < 2.) //sigle peak
      fQuality = 999. ;
    else
      fQuality =   0. ;
    return kTRUE ;
  }
      
  // if sample has reasonable mean and RMS, try to fit it with gamma2
  
  gMinuit->mncler();                     // Reset Minuit's list of paramters
  gMinuit->SetPrintLevel(-1) ;           // No Printout
  gMinuit->SetFCN(AliPHOSRawFitterv1::UnfoldingChiSquare) ;  

  // To set the address of the minimization function 
  
  fToFit->Clear("nodelete") ;
  Double_t b=0,bmin=0,bmax=0 ;
  if      (fCaloFlag == 0){ // Low gain
    fSampleParamsLow->AddAt(pedestal,4) ;
    if (fOverflow)
      fSampleParamsLow->AddAt(double(maxAmp),5) ;
    else
      fSampleParamsLow->AddAt(double(1023),5) ;
    fSampleParamsLow->AddAt(double(iStart),6) ;
    fToFit->AddFirst((TObject*)fSampleParamsLow) ; 
    b=fSampleParamsLow->At(2) ;
    bmin=0.5 ;
    bmax=10. ;
  }
  else if (fCaloFlag == 1){ // High gain
    fSampleParamsHigh->AddAt(pedestal,4) ;
    if (fOverflow)
      fSampleParamsHigh->AddAt(double(maxAmp),5) ;
    else
      fSampleParamsHigh->AddAt(double(1023),5);
    fSampleParamsHigh->AddAt(double(iStart),6);
    fToFit->AddFirst((TObject*)fSampleParamsHigh) ; 
    b=fSampleParamsHigh->At(2) ;
    bmin=0.05 ;
    bmax=0.4 ;
  }
  fToFit->AddLast((TObject*)fSamples) ;
  fToFit->AddLast((TObject*)fTimes) ;
  
  gMinuit->SetObjectFit((TObject*)fToFit) ;         // To tranfer pointer to UnfoldingChiSquare
  Int_t ierflg ;
  gMinuit->mnparm(0, "t0",  1.*tStart, 0.01, -500., 500., ierflg) ;
  if(ierflg != 0){
    //	  AliWarning(Form("Unable to set initial value for fit procedure : t0=%e\n",1.*tStart) ) ;
   fEnergy =   0. ;
    fTime   =-999. ;
    fQuality= 999. ;
    return kTRUE ; //will scan further
  }
  Double_t amp0=0; 
  if      (fCaloFlag == 0) // Low gain
    amp0 = fEnergy/sampleMaxLG;
  else if (fCaloFlag == 1) // High gain
    amp0 = fEnergy/sampleMaxHG;
  
  gMinuit->mnparm(1, "Energy", amp0 , 0.01*amp0, 0, 0, ierflg) ;
  if(ierflg != 0){
    //	  AliWarning(Form("Unable to set initial value for fit procedure : E=%e\n", amp0)) ;
    fEnergy =   0. ;
    fTime   =-999. ;
    fQuality= 999. ;
    return kTRUE ; //will scan further
  }
  
  gMinuit->mnparm(2, "p2", b, 0.01*b, bmin, bmax, ierflg) ;
  if(ierflg != 0){                                         
    //        AliWarning(Form("Unable to set initial value for fit procedure : E=%e\n", amp0)) ;  
    fEnergy =   0. ;
    fTime   =-999. ;
    fQuality= 999. ;
    return kTRUE ; //will scan further  
  }             
  
  Double_t p0 = 0.0001 ; // "Tolerance" Evaluation stops when EDM = 0.0001*p0 ; The number of function call slightly
  //  depends on it. 
  Double_t p1 = 1.0 ;
  Double_t p2 = 0.0 ;
  gMinuit->mnexcm("SET STR", &p2, 0, ierflg) ;   // force TMinuit to reduce function calls  
  gMinuit->mnexcm("SET GRA", &p1, 1, ierflg) ;   // force TMinuit to use my gradient  
  //	gMinuit->SetMaxIterations(100);
  gMinuit->mnexcm("SET NOW", &p2 , 0, ierflg) ;  // No Warnings
  
  gMinuit->mnexcm("MIGRAD", &p0, 0, ierflg) ;    // minimize 
  
  Double_t err=0.,t0err=0. ;
  Double_t t0=0.,efit=0. ;
  gMinuit->GetParameter(0,t0, t0err) ;
  gMinuit->GetParameter(1,efit, err) ;
  
  Double_t bfit=0., berr=0. ;
  gMinuit->GetParameter(2,bfit,berr) ;
  
  //Calculate total energy
  //this is parameterization of dependence of pulse height on parameter b
  if(fCaloFlag == 0) // Low gain
    efit *= 99.54910 + 78.65038*bfit ;
  else if(fCaloFlag == 1) // High gain
    efit *= 80.33109 + 128.6433*bfit ;
  
  if(efit < 0. || efit > 10000.){
    //set energy to previously found max
    fTime   =-999.;
    fQuality= 999 ;
    return kTRUE;
  }                                                                             
  
  //evaluate fit quality
  Double_t fmin=0.,fedm=0.,errdef=0. ;
  Int_t npari,nparx,istat;
  gMinuit->mnstat(fmin,fedm,errdef,npari,nparx,istat) ;
  fQuality = fmin/sigLength ;
  //compare quality with some parameterization
  if      (fCaloFlag == 0) // Low gain
    fQuality /= 2.00 + 0.0020*fEnergy ;
  else if (fCaloFlag == 1) // High gain
    fQuality /= 0.75 + 0.0025*fEnergy ;
  
  fEnergy = efit ;
  fTime  += t0 - 4.024*bfit ; //-10.402*bfit+4.669*bfit*bfit ; //Correction for 70 samples
//  fTime  += sigStart;
  
  delete fSamples ;
  delete fTimes ;
  return kTRUE;
}
//-----------------------------------------------------------------------------
Double_t AliPHOSRawFitterv1::Gamma2(Double_t dt,Double_t en,Double_t b,TArrayD * params){  //Function for fitting samples
  //parameters:
  //dt-time after start
  //en-amplutude
  //function parameters
  
  Double_t ped=params->At(4) ;
  if(dt<0.)
    return ped ; //pedestal
  else
    return ped+en*(TMath::Power(dt,params->At(0))*TMath::Exp(-dt*params->At(1))+b*dt*dt*TMath::Exp(-dt*params->At(3))) ;
}
//_____________________________________________________________________________
void AliPHOSRawFitterv1::UnfoldingChiSquare(Int_t & /*nPar*/, Double_t * Grad, Double_t & fret, Double_t * x, Int_t iflag)
{
  // Number of parameters, Gradient, Chi squared, parameters, what to do

  TList * toFit= (TList*)gMinuit->GetObjectFit() ;
  TArrayD * params=(TArrayD*)toFit->At(0) ; 
  TArrayI * samples = (TArrayI*)toFit->At(1) ;
  TArrayI * times = (TArrayI*)toFit->At(2) ;

  fret = 0. ;     
  if(iflag == 2)
    for(Int_t iparam = 0 ; iparam < 3 ; iparam++)    
      Grad[iparam] = 0 ; // Will evaluate gradient
  
  Double_t t0=x[0] ;
  Double_t en=x[1] ;
  Double_t b=x[2] ;
  Double_t n=params->At(0) ;
  Double_t alpha=params->At(1) ;
  Double_t beta=params->At(3) ;
  //  Double_t ped=params->At(4) ;

  Double_t overflow=params->At(5) ;
  Int_t iBin = (Int_t) params->At(6) ;
  Int_t nSamples=TMath::Min(iBin+70,samples->GetSize()) ; //Here we set number of points to fit (70)
  // iBin - first non-zero sample 
  Int_t tStep=times->At(iBin+1)-times->At(iBin) ;
  Double_t ddt=times->At(iBin)-t0-tStep ;
  Double_t exp1=TMath::Exp(-alpha*ddt) ;
  Double_t exp2=TMath::Exp(-beta*ddt) ;
  Double_t dexp1=TMath::Exp(-alpha*tStep) ;
  Double_t dexp2=TMath::Exp(-beta*tStep) ;
  for(Int_t i = iBin; i<nSamples ; i++) {
    Double_t dt=double(times->At(i))-t0 ;
    Double_t fsample = double(samples->At(i)) ;
    Double_t diff=0. ;
    exp1*=dexp1 ;
    exp2*=dexp2 ;
//    if(fsample>=overflow)
//      continue ;    
    if(dt<=0.){
      diff=fsample ; 
      fret += diff*diff ;
      continue ;
    }
    
    Double_t dtn=TMath::Power(dt,n) ;
    Double_t dtnE=dtn*exp1 ;
    Double_t dt2E=dt*dt*exp2 ;
    Double_t fit=en*(dtnE + b*dt2E) ;
    if(fsample>=overflow && fit>=overflow)
      continue ;    

    diff = fsample - fit ;
    fret += diff*diff ;
    if(iflag == 2){  // calculate gradient
      Grad[0] += en*diff*(dtnE*(n/dt-alpha)+b*dt2E*(2./dt-beta))  ; //derivative over t0
      Grad[1] -= diff*(dtnE+b*dt2E) ;
      Grad[2] -= en*diff*dt2E ;
    }
    
  }
  if(iflag == 2)
    for(Int_t iparam = 0 ; iparam < 3 ; iparam++)    
      Grad[iparam] *= 2. ; 
}
