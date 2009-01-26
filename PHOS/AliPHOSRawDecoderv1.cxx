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

// This class decodes the stream of ALTRO samples to extract
// the PHOS "digits" of current event. Uses fitting procedure
// to separate reasonable samples
// 
// Typical use case:
//     AliRawReader* rf = new AliRawReaderDate("2006run2211.raw");
//     AliPHOSRawDecoder dc(rf);
//     while (rf->NextEvent()) {
//       dc.SubtractPedestals(kTRUE);
//       while ( dc.NextDigit() ) {
//         Int_t module = dc.GetModule();
//         Int_t column = dc.GetColumn();
//         Int_t row = dc.GetRow();
//         Double_t amplitude = dc.GetEnergy();
//         Double_t time = dc.GetTime();
//         Bool_t IsLowGain = dc.IsLowGain();
//            ..........
//       }
//     }

// Author: Dmitri Peressounko

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
#include "AliPHOSCalibData.h"
#include "AliPHOSRawDecoderv1.h"
#include "AliPHOSPulseGenerator.h"

ClassImp(AliPHOSRawDecoderv1)

//-----------------------------------------------------------------------------
  AliPHOSRawDecoderv1::AliPHOSRawDecoderv1():AliPHOSRawDecoder(),
fSampleParamsLow(0x0),fSampleParamsHigh(0x0),fToFit(0x0)
{
  //Default constructor.
}

//-----------------------------------------------------------------------------
AliPHOSRawDecoderv1::AliPHOSRawDecoderv1(AliRawReader* rawReader,  AliAltroMapping **mapping):
  AliPHOSRawDecoder(rawReader,mapping),
  fSampleParamsLow(0x0),fSampleParamsHigh(0x0),fToFit(0x0)
{
  //Construct a decoder object.
  //Is is user responsibility to provide next raw event 
  //using AliRawReader::NextEvent().

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
AliPHOSRawDecoderv1::~AliPHOSRawDecoderv1()
{
  //Destructor.
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
AliPHOSRawDecoderv1::AliPHOSRawDecoderv1(const AliPHOSRawDecoderv1 &phosDecoder ):
  AliPHOSRawDecoder(phosDecoder), 
  fSampleParamsLow(0x0),fSampleParamsHigh(0x0),fToFit(0x0)
{
  //Copy constructor.
  fToFit = new TList() ;
  fSampleParamsLow =new TArrayD(*(phosDecoder.fSampleParamsLow)) ;
  fSampleParamsHigh=new TArrayD(*(phosDecoder.fSampleParamsHigh)) ;
}

//-----------------------------------------------------------------------------
AliPHOSRawDecoderv1& AliPHOSRawDecoderv1::operator = (const AliPHOSRawDecoderv1 &phosDecoder)
{
  //Assignment operator.

  //  if(this != &phosDecoder) {
  //  }
  fToFit = new TList() ;
  if(fSampleParamsLow){
    fSampleParamsLow = phosDecoder.fSampleParamsLow ;
    fSampleParamsHigh= phosDecoder.fSampleParamsHigh ;
  }
  else{
    fSampleParamsLow =new TArrayD(*(phosDecoder.fSampleParamsLow)) ; 
    fSampleParamsHigh=new TArrayD(*(phosDecoder.fSampleParamsHigh)) ;
  }
  return *this;
}

//-----------------------------------------------------------------------------
Bool_t AliPHOSRawDecoderv1::NextDigit()
{
  //Extract an energy deposited in the crystal,
  //crystal' position (module,column,row),
  //time and gain (high or low).
  //First collects sample, then evaluates it and if it has
  //reasonable shape, fits it with Gamma2 function and extracts 
  //energy and time.

//Debug=====================
//  TCanvas * c = 0; //(TCanvas*)gROOT->FindObjectAny("CSample") ;
//  if(!c)
//    c = new TCanvas("CSample","CSample") ;
// 
//  TH1D * h = 0 ; //(TH1D*)gROOT->FindObjectAny("hSample") ;
//  if(!h)
//    h=new TH1D("hSample","",200,0.,200.) ;
// 
//  TF1 * fff = 0 ; //(TF1*)gROOT->FindObjectAny("fff") ;
//  if(!fff)
//    fff = new TF1("fff","[0]+[1]*((abs(x-[2]))^[3]*exp(-(x-[2])*[4])+[5]*(x-[2])*(x-[2])*exp(-(x-[2])*[6]))",0.,1000.) ;
//End debug===========
  
  AliCaloRawStream* in = fCaloStream;
  
  Int_t    iBin     = fSamples->GetSize() ;
  Int_t    tLength  = 0;
  fEnergy = -111;
  Float_t pedMean = 0;
  Float_t pedRMS = 0;
  Int_t   nPed = 0;
  Float_t baseLine = 1.0;
  const Float_t nPreSamples = 10;
  fQuality= 999. ;
  const Float_t sampleMaxHG=102.332 ;  //maximal height of HG sample with given parameterization
  const Float_t sampleMaxLG=277.196 ;  //maximal height of HG sample with given parameterization
  const Float_t maxEtoFit=5 ; //fit only samples above this energy, accept all samples (with good aRMS) below it
  fSamples->Reset();
  fTimes  ->Reset();

  while ( in->Next() ) { 

    if(!tLength) {
      tLength = in->GetTimeLength();
      if(tLength!=fSamples->GetSize()) {
	delete fSamples ;
	delete fTimes ;
	fSamples = new TArrayI(tLength);
	fTimes = new TArrayI(tLength);
        iBin= fSamples->GetSize() ;
      }
      else{
        fSamples->Reset() ;
      }
    }
    
    // Fit the full sample
    if((in->IsNewHWAddress() && iBin != fSamples->GetSize()) //new HW address
       ||(iBin<=0)) {  //or new signal in same address

      //First remember new sample
      fNewLowGainFlag = in->IsLowGain();
      fNewModule      = in->GetModule()+1;
      fNewRow         = in->GetRow()   +1;
      fNewColumn      = in->GetColumn()+1;
      fNewAmp         = in->GetSignal() ;
      fNewTime        = in->GetTime() ;  
  
      //now handle already collected 
      Double_t pedestal =0. ;
      fPedestalRMS=0. ;
      if(fPedSubtract){ 
	if (nPed > 0){
	  pedestal = (Double_t)(pedMean/nPed); 
          fPedestalRMS=pedRMS/nPed-pedestal*pedestal ;
          if(fPedestalRMS>0.) fPedestalRMS=TMath::Sqrt(fPedestalRMS) ;
        }
	else
	  return kFALSE;
      }
      else{
        //take pedestals from DB
        pedestal = fAmpOffset ;
        if(fCalibData){
           Float_t truePed = fCalibData->GetADCpedestalEmc(fModule, fColumn, fRow) ;
           Int_t   altroSettings = fCalibData->GetAltroOffsetEmc(fModule, fColumn, fRow) ;
           pedestal += truePed - altroSettings ;
         }
         else{
//           printf("AliPHOSRawDecoderv1::NextDigit(): Can not read data from OCDB \n") ;
         }
      }

      //calculate time and energy
      Int_t maxBin=0 ;
      Int_t maxAmp=0 ;
      Double_t aMean=0. ;
      Double_t aRMS=0. ;
      Double_t wts=0 ;
      Int_t tStart = 0 ;
      for(Int_t i=iBin; i<fSamples->GetSize(); i++){
        if(fSamples->At(i)>pedestal){
          Double_t de=fSamples->At(i)-pedestal ;
          if(de>1.){
            aMean+=de*i ;
            aRMS+=de*i*i ;
            wts+=de; 
          }
          if(de>2 && tStart==0) 
            tStart=i ;
          if(maxAmp<fSamples->At(i)){
            maxBin=i ;
            maxAmp=fSamples->At(i) ;
          }
        }
      }
      if(maxBin==fSamples->GetSize()-1){//bad "rising" sample
        fEnergy=0. ;
        fTime=-999.;
        fQuality= 999. ;
        return kTRUE ;
      }
      fEnergy=Double_t(maxAmp)-pedestal ;
      fOverflow =0 ;  //look for plato on the top of sample
      if(fEnergy>500 &&  //this is not fluctuation of soft sample
         maxBin<fSamples->GetSize()-1 && fSamples->At(maxBin+1)==maxAmp){ //and there is a plato
         fOverflow = kTRUE ;
      }
      
      if(wts>0){
	aMean/=wts; 
	aRMS=aRMS/wts-aMean*aMean;
      }

      //do not take too small energies
      if(fEnergy < baseLine) 
         fEnergy = 0;

      //do not test quality of too soft samples
      if(fEnergy<maxEtoFit){
        fTime=fTimes->At(tStart);
        if(aRMS<2.) //sigle peak
          fQuality=999. ;
        else
          fQuality= 0. ;
        return kTRUE ;
      }

      
//Debug:=====Draw sample
//if(fEnergy>pedestal+10.){
//if(fLowGainFlag && fEnergy>2){
//  if(!c)
//    if(!fLowGainFlag && fRow==32 && fColumn==18){
//    TCanvas *c = new TCanvas("CSample","CSample") ;
//    c->cd() ;
//    h->Draw() ;
//    c->Update() ;
// printf("fEnergy=%f, aRMS=%f \n",fEnergy,aRMS) ;   
//getchar() ;
//}
//======================

      //IF sample has reasonable mean and RMS, try to fit it with gamma2
	
	gMinuit->mncler();                     // Reset Minuit's list of paramters
	gMinuit->SetPrintLevel(-1) ;           // No Printout
	gMinuit->SetFCN(AliPHOSRawDecoderv1::UnfoldingChiSquare) ;  
	// To set the address of the minimization function 
 	
        fToFit->Clear("nodelete") ;
	Double_t b,bmin,bmax ;
	if(fLowGainFlag){
	  fSampleParamsLow->AddAt(pedestal,4) ;
	  if(fOverflow)
	    fSampleParamsLow->AddAt(double(maxAmp),5) ;
	  else
	    fSampleParamsLow->AddAt(double(1023),5) ;
	  fSampleParamsLow->AddAt(double(iBin),6) ;
	  fToFit->AddFirst((TObject*)fSampleParamsLow) ; 
	  b=fSampleParamsLow->At(2) ;
	  bmin=0.5 ;
	  bmax=10. ;
	}
	else{
	  fSampleParamsHigh->AddAt(pedestal,4) ;
	  if(fOverflow)
	    fSampleParamsHigh->AddAt(double(maxAmp),5) ;
	  else
	    fSampleParamsHigh->AddAt(double(1023),5);
	  fSampleParamsHigh->AddAt(double(iBin),6);
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
	  fEnergy=0. ;
	  fTime=-999. ;
          fQuality=999 ;
	  return kTRUE ; //will scan further
	}
        Double_t amp0; 
        if(fLowGainFlag)
          amp0=fEnergy/sampleMaxLG;
        else
          amp0=fEnergy/sampleMaxHG;

	gMinuit->mnparm(1, "Energy", amp0 , 0.01*amp0, 0, 0, ierflg) ;
	if(ierflg != 0){
//	  AliWarning(Form("Unable to set initial value for fit procedure : E=%e\n", amp0)) ;
	  fEnergy=0. ;
	  fTime=-999. ;
          fQuality=999 ;
	  return kTRUE ; //will scan further
	}

        gMinuit->mnparm(2, "p2", b, 0.01*b, bmin, bmax, ierflg) ;
        if(ierflg != 0){                                         
//        AliWarning(Form("Unable to set initial value for fit procedure : E=%e\n", amp0)) ;  
          fEnergy=0. ;           
          fTime=-999. ;         
          fQuality=999 ;       
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
	
	Double_t err,t0err ;
	Double_t t0,efit ;
	gMinuit->GetParameter(0,t0, t0err) ;    
	gMinuit->GetParameter(1,efit, err) ;    

        Double_t bfit, berr ;
        gMinuit->GetParameter(2,bfit,berr) ;

        //Calculate total energy
        //this isparameterization of depetendence of pulse height on parameter b
        if(fLowGainFlag)
          efit*=99.54910 + 78.65038*bfit ;
        else
          efit*=80.33109+128.6433*bfit ;

        if(efit<0. || efit > 10000.){
//set energy to previously found max
//          fEnergy=0 ; //bad sample                                                    
          fTime=-999.;                                                                
          fQuality=999 ;                                                              
          return kTRUE;
        }                                                                             
 
        //evaluate fit quality
	Double_t fmin,fedm,errdef ;
	Int_t npari,nparx,istat;
	gMinuit->mnstat(fmin,fedm,errdef,npari,nparx,istat) ;
        fQuality=fmin/(fSamples->GetSize()-iBin) ;
        //compare quality with some parameterization
        if(fLowGainFlag){
          fQuality/=2.+0.002*fEnergy ;
        }
        else{
          fQuality/=0.75+0.0025*fEnergy ;
        }

//Debug================
//        Double_t n,alpha,beta ;
//        Double_t en ;
//       if(fLowGainFlag){
//          n=fSampleParamsLow->At(0) ;
//          alpha=fSampleParamsLow->At(1) ;
//          beta=fSampleParamsLow->At(3) ;
//          en=efit/(99.54910 + 78.65038*bfit) ;
//        }
//        else{
//          n=fSampleParamsHigh->At(0) ;
//          alpha=fSampleParamsHigh->At(1) ;
//          beta=fSampleParamsHigh->At(3) ;
//          en=efit/(80.33109+128.6433*bfit) ;
//        }
//
////    if( fQuality > 1 && fEnergy > 20. && !fOverflow){
////    if(!fLowGainFlag && fRow==32 && fColumn==18){
//{
//    printf("Col=%d, row=%d, qual=%f, E=%f, t0=%f, b=%f\n",fColumn,fRow,fQuality,efit,t0,bfit) ;
//    printf("    Energy = %f \n",fEnergy) ;
//    TCanvas * c = new TCanvas("samp") ;
//    c->cd() ;
//    h->Draw() ;
//    if(fLowGainFlag){
//      fff->SetParameters(pedestal,en,t0,n,alpha,bfit,beta) ;
//    }
//    else{
//     fff->SetParameters(pedestal,en,t0,n,alpha,bfit,beta) ;
//    }
//////    for(Int_t i=1;i<=h->GetNbinsX(); i++){
////       Double_t x=h->GetBinCenter(i) ;
////       h->SetBinContent(i,h->GetBinContent(i)-fff->Eval(x)) ;
////    }
////    h->SetMinimum(-15.) ;
////    h->SetMaximum(15.) ;
//    h->Draw() ;
//    fff->Draw("same") ;
//    c->Update();
//    getchar() ;
//    }
//====================

      fEnergy=efit ;
      fTime=t0-4.024*bfit ; //-10.402*bfit+4.669*bfit*bfit ; //Correction for 70 samples
//      fTime=t0+2.8*bfit ; //-10.402*bfit+4.669*bfit*bfit ; //Correction for 50 samples
//      fQuality = bfit ;
      return kTRUE;
    }
    
    fLowGainFlag = in->IsLowGain();
    fModule = in->GetModule()+1;
    fRow    = in->GetRow()   +1;
    fColumn = in->GetColumn()+1;

    //add previouly taken if coincides
    if(fLowGainFlag==fNewLowGainFlag && fModule==fNewModule &&
       fRow==fNewRow && fColumn==fNewColumn){
       iBin--;
       if(fPedSubtract && fNewTime < nPreSamples) {
         pedMean += in->GetSignal();
         pedRMS += in->GetSignal()*in->GetSignal() ;
         nPed++;
       }
       fSamples->AddAt(fNewAmp,iBin);
       fTimes->AddAt(fNewTime,iBin);
    
       //Mark that we already take it
       fNewModule=-1 ;
    }
    
    // Fill array with samples
    iBin--;
    if(fPedSubtract && (in->GetTime() < nPreSamples)) {
      pedMean += in->GetSignal();
      pedRMS += in->GetSignal()*in->GetSignal() ;
      nPed++;
    }
    fSamples->AddAt(in->GetSignal(),iBin);
    fTimes->AddAt(in->GetTime(),iBin);
 
//Debug==============
//    h->SetBinContent(in->GetTime(),in->GetSignal()) ;
//EndDebug==============
    
  } // in.Next()
  
  return kFALSE;
}
//_____________________________________________________________________________
void AliPHOSRawDecoderv1::UnfoldingChiSquare(Int_t & /*nPar*/, Double_t * Grad, Double_t & fret, Double_t * x, Int_t iflag)
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
  Double_t ped=params->At(4) ;

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
    if(fsample>=overflow)
      continue ;
    Double_t diff ;
    exp1*=dexp1 ;
    exp2*=dexp2 ;
    if(dt<=0.){
      diff=fsample - ped ; 
      fret += diff*diff ;
      continue ;
    }
    Double_t dtn=TMath::Power(dt,n) ;
    Double_t dtnE=dtn*exp1 ;
    Double_t dt2E=dt*dt*exp2 ;
    Double_t fit=ped+en*(dtnE + b*dt2E) ;
//    if(fit>=overflow){
//      diff=fsample-overflow ;
//      fret += diff*diff ;
//      //zero gradient here
//    }
//    else{
      diff = fsample - fit ;
      fret += diff*diff ;
      if(iflag == 2){  // calculate gradient
        Grad[0] += en*diff*(dtnE*(n/dt-alpha)+b*dt2E*(2./dt-beta))  ; //derivative over t0
        Grad[1] -= diff*(dtnE+b*dt2E) ;
        Grad[2] -= en*diff*dt2E ;
      }
//    }
  }
  if(iflag == 2)
    for(Int_t iparam = 0 ; iparam < 3 ; iparam++)    
      Grad[iparam] *= 2. ; 
}
//-----------------------------------------------------------------------------
Double_t AliPHOSRawDecoderv1::Gamma2(Double_t dt,Double_t en,Double_t b,TArrayD * params){  //Function for fitting samples
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

