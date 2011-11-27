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


// This class extracts the signal parameters (energy, time, quality)
// from ALTRO samples. Energy is in ADC counts, time is in time bin units.
// If sample is not in saturation, a coarse algorithm is used (a-la AliPHOSRawFitterv0)
// If sample is in saturation, the unsaturated part of the sample is fit a-la AliPHOSRawFitterv1
// 
//     AliPHOSRawFitterv4 *fitterv4=new AliPHOSRawFitterv4();
//     fitterv4->SetChannelGeo(module,cellX,cellZ,caloFlag);
//     fitterv4->SetCalibData(fgCalibData) ;
//     fitterv4->Eval(sig,sigStart,sigLength);
//     Double_t amplitude = fitterv4.GetEnergy();
//     Double_t time      = fitterv4.GetTime();
//     Bool_t   isLowGain = fitterv4.GetCaloFlag()==0;

// Author: Dmitri Peressounko

// --- ROOT system ---
#include "TArrayI.h"
#include "TMath.h"
#include "TObject.h"
#include "TArrayD.h"
#include "TList.h"
#include "TMinuit.h"


//Used for debug
//#include "TROOT.h"
//#include "TH1.h"
//#include "TCanvas.h"
//#include "TPad.h"
//#include "TF1.h"



// --- AliRoot header files ---
#include "AliPHOSRawFitterv4.h"
#include "AliPHOSCalibData.h"
#include "AliLog.h"

ClassImp(AliPHOSRawFitterv4)

//-----------------------------------------------------------------------------
AliPHOSRawFitterv4::AliPHOSRawFitterv4():
  AliPHOSRawFitterv1(),
  fFitHighGain(0) 
{
  //Default constructor
}

//-----------------------------------------------------------------------------
AliPHOSRawFitterv4::~AliPHOSRawFitterv4()
{
}

//-----------------------------------------------------------------------------
AliPHOSRawFitterv4::AliPHOSRawFitterv4(const AliPHOSRawFitterv4 &phosFitter ):
  AliPHOSRawFitterv1((AliPHOSRawFitterv1)phosFitter),
  fFitHighGain(0) 
{
  //Copy constructor
  fFitHighGain=phosFitter.fFitHighGain;
}


//-----------------------------------------------------------------------------

Bool_t AliPHOSRawFitterv4::Eval(const UShort_t *signal, Int_t sigStart, Int_t sigLength)
{
  // Calculate signal parameters (energy, time, quality) from array of samples
  // Energy is a maximum sample minus pedestal 9
  // Time is the first time bin
  // Signal overflows is there are at least 3 samples of the same amplitude above 900

  fOverflow= kFALSE ;
  fEnergy  = 0;
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

  for (Int_t i=0; i<sigLength; i++) {
    if (i>sigLength-kPreSamples) { //inverse signal time order
      nPed++;
      pedMean += signal[i];
      pedRMS  += signal[i]*signal[i] ;
    }
    if(signal[i] > maxSample ){ maxSample = signal[i]; nMax=0;}
    if(signal[i] == maxSample) nMax++;

  }

  fEnergy = (Double_t)maxSample;
  if (maxSample > 850 && nMax > 2) fOverflow = kTRUE;
  
  Double_t pedestal = 0 ;
  if (fPedSubtract) {
    if (nPed > 0) {
      fPedestalRMS=(pedRMS - pedMean*pedMean/nPed)/nPed ;
      if(fPedestalRMS > 0.) 
	fPedestalRMS = TMath::Sqrt(fPedestalRMS) ;
      pedestal = (Double_t)(pedMean/nPed);
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
      AliDebug(2,Form("Pedestal and offset are not read from OCDB. Use 0 for their values.")) ;
    }
  }
  fEnergy-=pedestal ;
  if (fEnergy < kBaseLine) fEnergy = 0;
  
  
  if(fOverflow && ((fCaloFlag == 0) || fFitHighGain)){ //by default fit LowGain only
    TArrayI *samples = new TArrayI(sigLength); // array of sample amplitudes
    TArrayI *times   = new TArrayI(sigLength); // array of sample time stamps
    Bool_t result = kTRUE ;
    for (Int_t i=0; i<sigLength; i++) {
      samples->AddAt(signal[i]-pedestal,sigLength-i-1);
      times  ->AddAt(i ,i);
    }
    //Prepare fit parameters
    //Evaluate time
    Int_t iStart = 0;
    while(iStart<sigLength && samples->At(iStart) <kBaseLine) iStart++ ;
    if      (fCaloFlag == 0){ // Low gain
      fSampleParamsLow->AddAt(pedestal,4) ;
      fSampleParamsLow->AddAt(double(maxSample),5) ;
      fSampleParamsLow->AddAt(double(iStart),6) ;
      fToFit->AddFirst((TObject*)fSampleParamsLow) ; 
    }
    else if (fCaloFlag == 1){ // High gain
      fSampleParamsHigh->AddAt(pedestal,4) ;
      fSampleParamsHigh->AddAt(double(maxSample),5) ;
      fSampleParamsHigh->AddAt(double(iStart),6) ;
      fToFit->AddFirst((TObject*)fSampleParamsHigh) ; 
    }
    result=EvalWithFitting(samples,times); 
    delete samples ;
    delete times ;
    
    return result; 
  
  }
  
  
  //Sample withour overflow - Evaluate time
  fTime = sigStart-sigLength-3; 
  const Int_t nLine= 6 ;        //Parameters of fitting
  const Float_t eMinTOF = 10. ; //Choosed from beam-test and cosmic analyis
  const Float_t kAmp=0.35 ;     //Result slightly depends on them, so no getters
  // Avoid too low peak:
  if(fEnergy < eMinTOF){
     return kTRUE;
  }
  // Find index posK (kLevel is a level of "timestamp" point Tk):
  Int_t posK =sigLength-1 ; //last point before crossing k-level
  Double_t levelK = pedestal + kAmp*fEnergy;
  while(signal[posK] <= levelK && posK>=0){
     posK-- ;
  }
  posK++ ;

  if(posK == 0 || posK==sigLength-1){
    return kTRUE; 
  }

  // Find crosing point by solving linear equation (least squares method)
  Int_t np = 0;
  Int_t iup=posK-1 ;
  Int_t idn=posK ;
  Double_t sx = 0., sy = 0., sxx = 0., sxy = 0.;
  Double_t x,y ;

  while(np<nLine){
    //point above crossing point
    if(iup>=0){
      x = sigLength-iup-1;
      y = signal[iup];
      sx += x;
      sy += y;
      sxx += (x*x);
      sxy += (x*y);
      np++ ;
      iup-- ;
    }
    //Point below crossing point
    if(idn<sigLength){
      if(signal[idn]<pedestal){
        idn=sigLength-1 ; //do not scan further
	idn++ ;
        continue ;
      }
      x = sigLength-idn-1;
      y = signal[idn];
      sx += x;
      sy += y;
      sxx += (x*x);
      sxy += (x*y);
      np++;
      idn++ ;
    }
    if(idn>=sigLength && iup<0){
      break ; //can not fit futher
    }
  }

  Double_t det = np*sxx - sx*sx;
  if(det == 0){
    return kTRUE;
  }
  Double_t c1 = (np*sxy - sx*sy)/det;  //slope
  Double_t c0 = (sy-c1*sx)/np; //offset
  if(c1 == 0){
    return kTRUE;
  }

  // Find where the line cross kLevel:
  fTime += (levelK - c0)/c1-5. ; //5: mean offset between k-Level and start times
  return kTRUE;

}
//===================================================
Bool_t AliPHOSRawFitterv4::EvalWithFitting(TArrayI*samples, TArrayI * times){

  // if sample has reasonable mean and RMS, try to fit it with gamma2
  const Float_t sampleMaxHG=102.332 ; 
  const Float_t sampleMaxLG=277.196 ; 
  
  gMinuit->mncler();                     // Reset Minuit's list of paramters
  gMinuit->SetPrintLevel(-1) ;           // No Printout
  gMinuit->SetFCN(AliPHOSRawFitterv1::UnfoldingChiSquare) ;  
  // To set the address of the minimization function 
  
  fToFit->Clear("nodelete") ;
  Double_t b=0,bmin=0,bmax=0 ;
  if      (fCaloFlag == 0){ // Low gain
    b=fSampleParamsLow->At(2) ;
    bmin=0.5 ;
    bmax=10. ;
    fToFit->AddFirst((TObject*)fSampleParamsLow) ; 
  }
  else if (fCaloFlag == 1){ // High gain
    b=fSampleParamsHigh->At(2) ;
    bmin=0.05 ;
    bmax=0.4 ;
    fToFit->AddFirst((TObject*)fSampleParamsHigh) ; 
  }
  fToFit->AddLast((TObject*)samples) ;
  fToFit->AddLast((TObject*)times) ;
  

  gMinuit->SetObjectFit((TObject*)fToFit) ;         // To tranfer pointer to UnfoldingChiSquare
  Int_t ierflg ;
  gMinuit->mnparm(0, "t0",  1., 0.01, -50., 50., ierflg) ;
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
  //////	gMinuit->SetMaxIterations(100);
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
    fEnergy = efit*(99.54910 + 78.65038*bfit) ;
  else if(fCaloFlag == 1) // High gain
    fEnergy=efit*(80.33109 + 128.6433*bfit) ;
  
  if(fEnergy < 0. || fEnergy > 10000.){
    //set energy to previously found max
    fTime   =-999.;
    fQuality= 999 ;
    fEnergy=0. ;
    return kTRUE;
  }                                                                             
  
  //evaluate fit quality
  Double_t fmin=0.,fedm=0.,errdef=0. ;
  Int_t npari,nparx,istat;
  gMinuit->mnstat(fmin,fedm,errdef,npari,nparx,istat) ;
  fQuality = fmin/samples->GetSize() ;
  //compare quality with some parameterization
  if      (fCaloFlag == 0) // Low gain
    fQuality /= 2.00 + 0.0020*fEnergy ;
  else if (fCaloFlag == 1) // High gain
    fQuality /= 0.75 + 0.0025*fEnergy ;
  
  fTime  += t0 - 4.024*bfit ;

  if(fQuality==0.){//no points to fit)
    fTime   =-999.;
    fQuality= 1999 ;
    fEnergy=0. ;  
  }
  
/*  
  if(1){
    TH1I * h = (TH1I*)gROOT->FindObjectAny("hSamples") ;
    if(!h) h = new TH1I("hSamples","Samples",65,0.,65.) ;
    h->Reset() ;
    for (Int_t i=0; i<samples->GetSize(); i++) {
        h->SetBinContent(i+1,float(samples->At(i))) ;
    }
    TF1 * fffit = new TF1("fffit","[0]*(abs(x-[1])^[3]*exp(-[2]*(x-[1]))+[4]*(x-[1])*(x-[1])*exp(-[5]*(x-[1])))",0.,60.) ;
    TArrayD * params=(TArrayD*)fToFit->At(0) ; 
    Double_t n=params->At(0) ;
    Double_t alpha=params->At(1) ;
    Double_t beta=params->At(3) ;
    fffit->SetParameters(efit,t0,alpha,n,bfit,beta) ;
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
  
//      TPad * spectrum_1 = new TPad("spectrum_1", "spectrum_1",0.001,0.32,0.99,0.99);
      TPad * spectrum_1 = new TPad("spectrum_1", "spectrum_1",0.001,0.001,0.999,0.999);
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


//       can->cd() ;
//       TPad *spectrum_2 = new TPad("spectrum_2", "spectrum_2",0.001,0.01,0.99,0.32);
//       spectrum_2->SetFillColor(0) ;
//       spectrum_2->SetFillStyle(0) ;
//       spectrum_2->SetGridy() ;
//       spectrum_2->Draw();
//       spectrum_2->Range(0,0,1,1);
//       spectrum_2->SetFillColor(0);
//       spectrum_2->SetBorderSize(1);
//       spectrum_2->SetTopMargin(0.01);
//       spectrum_2->SetBottomMargin(0.25);
//       spectrum_2->SetLeftMargin(0.10);
//       spectrum_2->SetRightMargin(0.05);
//       spectrum_2->cd() ;
// 
//       TH1I * hd = (TH1I*)gROOT->FindObjectAny("hSamplesDif") ;
//       if(!hd) hd = new TH1I("hd","Samples",65,0.,65.) ;
//       hd->Reset() ;
//       for (Int_t i=0; i<samples->GetSize(); i++) {
//         hd->SetBinContent(i+1,TMath::Max(-1023.,TMath::Min(1023.,samples->At(i)-fffit->Eval(i)))) ;
//       }
//       hd->Draw();

      can->Update() ;
      printf("Press <enter> to continue\n") ;
      getchar();


      delete fffit ;
      delete spectrum_1 ;
//      delete spectrum_2 ;
    }
*/
  
  return kTRUE;

}
