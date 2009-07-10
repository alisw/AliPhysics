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

// This class plots samples and qualify their quality according to their shape
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
#include "TList.h"
#include "TMath.h"
#include "TMinuit.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
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
  fNtimeSamples(25),
  fRMScut(11.)
{
  //Default constructor.
  fLGpar[0]=0.971 ;
  fLGpar[1]=0.0465;
  fLGpar[2]=1.56  ;
  fHGpar[0]=0.941 ; 
  fHGpar[1]=0.0436;
  fHGpar[2]=1.65  ;
}

//-----------------------------------------------------------------------------
AliPHOSRawFitterv2::~AliPHOSRawFitterv2()
{
  //Destructor.
}

//-----------------------------------------------------------------------------
AliPHOSRawFitterv2::AliPHOSRawFitterv2(const AliPHOSRawFitterv2 &phosFitter ):
  AliPHOSRawFitterv0(phosFitter), 
  fNtimeSamples(25),
  fRMScut(11.)
{
  //Copy constructor.
  fNtimeSamples=phosFitter.fNtimeSamples ;
  for(Int_t i=0; i<3;i++){
    fLGpar[i]=phosFitter.fLGpar[i] ;
    fHGpar[i]=phosFitter.fHGpar[i] ;
  }
  fRMScut=phosFitter.fRMScut ;
}

//-----------------------------------------------------------------------------
AliPHOSRawFitterv2& AliPHOSRawFitterv2::operator = (const AliPHOSRawFitterv2 &phosFitter)
{
  //Assignment operator.

  fNtimeSamples=phosFitter.fNtimeSamples ;
  for(Int_t i=0; i<3;i++){
    fLGpar[i]=phosFitter.fLGpar[i] ;
    fHGpar[i]=phosFitter.fHGpar[i] ;
  }
  fRMScut=phosFitter.fRMScut ;
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

  if (fNBunches > 1) {
    fQuality = 1000;
    return kTRUE;
  }

  const Float_t kBaseLine   = 1.0;
  const Int_t   kPreSamples = 10;

  Float_t  pedMean   = 0;
  Float_t  pedRMS    = 0;
  Int_t    nPed      = 0;
  UShort_t maxSample = 0;
  Int_t    nMax      = 0;

  TCanvas * cs = (TCanvas*)gROOT->FindObjectAny("CSample") ;
  if(!cs)
    cs = new TCanvas("CSample","CSample") ;

  TH1D * h = (TH1D*)gROOT->FindObjectAny("hSample") ;
  if(!h) h = new TH1D("hSample","",200,0.,200.) ;

  Double_t pedestal = 0;
  for (Int_t i=0; i<sigLength; i++) {
    if (i<kPreSamples) {
      nPed++;
      pedMean += signal[i];
      pedRMS  += signal[i]*signal[i] ;
    }
    if(signal[i] > maxSample) maxSample = signal[i];
    if(signal[i] == maxSample) nMax++;
    h->SetBinContent(i+1,signal[i]) ;
  }
  fEnergy = (Double_t)maxSample;
  if (maxSample > 900 && nMax > 2) fOverflow = kTRUE;

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

  if (fEnergy < kBaseLine) {
    fEnergy = 0;
    return kTRUE;
  }
  
  // calculate time
  fTime=0. ;
  Double_t tRMS = 0. ;
  Double_t tW = 0. ;
  Int_t cnts=0 ;
  Double_t a=0,b=0,c=0 ;
  if(fCaloFlag == 0){ // Low gain
    a=fLGpar[0] ; 
    b=fLGpar[1] ; 
    c=fLGpar[2] ; 
  }
  else if(fCaloFlag == 1){ // High gain
    a=fHGpar[0] ; 
    b=fHGpar[1] ; 
    c=fHGpar[2] ; 
  }


  fQuality = 0. ;
  
  for(Int_t i=1; i<sigLength && cnts<fNtimeSamples; i++){
    if(signal[i] < pedestal)
      continue ;
    Double_t de = signal[i]   - pedestal ;
    Double_t av = signal[i-1] - pedestal + de;
    if(av<=0.) //this is fluctuation around pedestal, skip it
      continue ;
    Double_t ds = signal[i] - signal[i-1] ;
    Double_t ti = ds/av ;       // calculate log. derivative
    ti = a/(ti+b)-c*ti ;        // and compare with parameterization
    ti = i - ti ; 
    Double_t wi = TMath::Abs(ds) ;
    fTime += ti*wi ;
    tW    += wi;
    tRMS  += ti*ti*wi ;
    cnts++ ;
  } 

  if(tW>0.){
    fTime/=tW ;
    fQuality = tRMS/tW-fTime*fTime ;
    fTime+=sigStart;
  }
  else{
    fTime=-999. ;
    fQuality=999. ;
  }

  Bool_t isBad = 0 ;
  for(Int_t i=1; i<sigLength-1&&!isBad; i++){
    if(signal[i] > signal[i-1]+5 && signal[i] > signal[i+1]+5) { //single jump
      isBad=1 ;
    }
  }
  if(pedestal < 10.)
    isBad=1 ;
  
  if(fPedestalRMS > 0.1)
    isBad=1 ;
  
  for(Int_t i=1; i<sigLength-1&&!isBad; i++){
    if(signal[i] < pedestal-1)
      isBad=1 ;
  }
  
  if(fEnergy>10. && !isBad ){
    printf("fE=%f, ped=%f, fQuality=%f, pedRMS=%f \n",fEnergy,pedestal,fQuality,pedRMS) ;
    if(fOverflow)printf(" Overflow \n") ;
    if(isBad)printf("bad") ;
    cs->cd() ;
    h->Draw() ;
    cs->Update() ;
    getchar() ;
  }

  return kTRUE;
}
