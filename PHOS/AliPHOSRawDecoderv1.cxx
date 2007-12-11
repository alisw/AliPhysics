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
//       dc.SetOldRCUFormat(kTRUE);
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
//#include "AliLog.h"
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
  fSampleParamsHigh =new TArrayD(5) ;
  fSampleParamsHigh->AddAt(4.25,0) ;
  fSampleParamsHigh->AddAt(0.094,1) ;
  fSampleParamsHigh->AddAt(0.0151,2) ;
  fSampleParamsHigh->AddAt(0.0384,3) ;
  fSampleParamsLow=new TArrayD(5) ;
  fSampleParamsLow->AddAt(5.14,0) ;
  fSampleParamsLow->AddAt(0.0970,1) ;
  fSampleParamsLow->AddAt(0.0088,2) ;
  fSampleParamsLow->AddAt(0.0346,3) ;
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

//  TCanvas * c  = (TCanvas *)gROOT->FindObjectAny("canvMy") ;
//  TH1S * h = new TH1S("s","",200,0.5,200.5) ;
//  TF1 * fff = new TF1("fff","[0]+[1]*((x-[2])+[3]*(x-[2])*(x-[2]))*(exp(-(x-[2])*[4])+[5]*exp(-(x-[2])*[6]))",0.,1000.) ;
  
  AliCaloRawStream* in = fCaloStream;
  
  Int_t    iBin     = 0;
  Int_t    tLength  = 0;
  fEnergy = -111;
  Float_t pedMean = 0;
  Int_t   nPed = 0;
  Float_t baseLine = 1.0;
  const Float_t nPreSamples = 10;
  
  while ( in->Next() ) { 
    
    if(!tLength) {
      tLength = in->GetTimeLength();
      if(tLength!=fSamples->GetSize()) {
	delete fSamples ;
	fSamples = new TArrayI(tLength);
      }
      else{
	for(Int_t i=0; i<fSamples->GetSize(); i++){
	  fSamples->AddAt(0,i) ;
	}
      }
    }
    
    // Fit the full sample
    if(in->IsNewHWAddress() && iBin>0) {
      
      Double_t pedestal =0. ;
      if(fPedSubtract){ 
	if (nPed > 0)
	  pedestal = (Double_t)(pedMean/nPed); 
	else
	  return kFALSE;
      }
      
      //calculate energy
      //first estimate if this sample looks like gamma2 function
      Double_t aMean=0. ;
      Double_t aRMS=0. ;
      Int_t tStart = 0 ;
      Int_t cnts=0 ;
      for(Int_t i=0; i<fSamples->GetSize(); i++){
	if(fSamples->At(i)>0){
	  Double_t de=fSamples->At(i)-pedestal ;
	  aMean+=de ;
	  aRMS+=de*de ;
	  cnts++;
	  if(de>2 && tStart==0)
	    tStart=i ;
	  if(fSamples->At(i)>fEnergy)
	    fEnergy=fSamples->At(i) ;
	}
      }
      if(cnts>0){
	aMean/=cnts; 
	aRMS=aRMS/cnts-aMean*aMean;
      }
      
      //IF sample has reasonable mean and RMS, try to fit it with gamma2
      if(fEnergy>2.&& cnts >20 && aMean>0. && aRMS>2.){ //more or less reasonable sample
	
	gMinuit->mncler();                     // Reset Minuit's list of paramters
	gMinuit->SetPrintLevel(-1) ;           // No Printout
	gMinuit->SetFCN(AliPHOSRawDecoderv1::UnfoldingChiSquare) ;  
	// To set the address of the minimization function 
 	
       fToFit->Clear() ;
       if(fLowGainFlag){
         fSampleParamsLow->AddAt(pedestal,4) ;
         fToFit->AddFirst((TObject*)fSampleParamsLow) ; 
        }
        else{
         fSampleParamsHigh->AddAt(pedestal,4) ;
         fToFit->AddFirst((TObject*)fSampleParamsHigh) ; 
        }
        fToFit->AddLast((TObject*)fSamples) ;

	gMinuit->SetObjectFit((TObject*)fToFit) ;         // To tranfer pointer to UnfoldingChiSquare
	Int_t ierflg ;
	gMinuit->mnparm(0, "t0",  1.*tStart, 0.1, 0, 0, ierflg) ;
	if(ierflg != 0){
//	  AliWarning(Form("Unable to set initial value for fit procedure : t0=%e\n",1.*tStart) ) ;
	  fEnergy=0. ;
	  fTime=-999. ;
	  return kTRUE ; //will scan further
	}
        Double_t amp0=(fEnergy-pedestal)*0.0032;

	gMinuit->mnparm(1, "Energy", amp0 , 0.001*amp0, 0, 0, ierflg) ;
	if(ierflg != 0){
//	  AliWarning(Form("Unable to set initial value for fit procedure : E=%e\n", amp0)) ;
	  fEnergy=0. ;
	  fTime=-999. ;
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
	
        Double_t a,alpha ;
        if(fLowGainFlag){
          a=fSampleParamsLow->At(0) ;
          alpha=fSampleParamsLow->At(1) ;
        }
        else{
          a=fSampleParamsHigh->At(0) ;
          alpha=fSampleParamsHigh->At(1) ;
        }

//    c->cd() ;
//    h->Draw() ;
//    if(fLowGainFlag){
//      fff->SetParameters(pedestal,efit,t0,a,alpha,fSampleParamsLow->At(2),fSampleParamsLow->At(3)) ;
//    }
//    else{
//      fff->SetParameters(pedestal,efit,t0,a,alpha,fSampleParamsHigh->At(2),fSampleParamsHigh->At(3)) ;
//    }
//    fff->Draw("same") ;
//    c->Update();
          
        efit*=(2.*a+TMath::Sqrt(4.*a*a+alpha*alpha))/alpha/alpha*TMath::Exp(-1.+(alpha-TMath::Sqrt(4.*a*a+alpha*alpha))/2./a) ;
//printf("efit=%f, t0=%e +- %e, ped=%f \n",efit,t0,t0err,pedestal) ;
	Double_t fmin,fedm,errdef ;
	Int_t npari,nparx,istat;
	gMinuit->mnstat(fmin,fedm,errdef,npari,nparx,istat) ;
  
//if(fLowGainFlag)
// printf("LowGain \n") ;
//else
// printf("highGain \n") ;

//printf("fmin=%e \n",fmin) ;
//getchar() ;	

	if(1){ //fmin < 3.+0.3*efit ){ //Chi^2 of a good sample
	  if(efit>0.){
	    fEnergy=efit ;
	    fTime=t0 ;
	  }
	}
	else{
	  fEnergy=0 ; //bad sample
	  fTime=-999.;
	}
	if(fLowGainFlag)
	  fEnergy *= fPulseGenerator->GetRawFormatHighLowGainFactor(); // *16 
	
        fTime*=fPulseGenerator->GetRawFormatTimeTrigger() ; 

	if (fEnergy < baseLine) fEnergy = 0;
      }
      else{ //bad sample
	fEnergy=0. ;
	fTime=-999. ;
      }
      
      return kTRUE;
    }
    
    fLowGainFlag = in->IsLowGain();
    fTime = fPulseGenerator->GetRawFormatTimeTrigger() * in->GetTime();
    fModule = in->GetModule()+1;
    fRow    = in->GetRow()   +1;
    fColumn = in->GetColumn()+1;
    
    
    // Fill array with samples
    iBin++;                                                             
    if(tLength-iBin < nPreSamples) {
      pedMean += in->GetSignal();
      nPed++;
    }
    fSamples->AddAt(in->GetSignal(),tLength-iBin);
//    h->SetBinContent(tLength-iBin+1,in->GetSignal()) ;
    
  } // in.Next()
  
  return kFALSE;
}
//_____________________________________________________________________________
void AliPHOSRawDecoderv1::UnfoldingChiSquare(Int_t & /*nPar*/, Double_t * Grad, Double_t & fret, Double_t * x, Int_t iflag)
{
  // Calculates the Chi square for the samples minimization
  // Number of parameters, Gradient, Chi squared, parameters, what to do

  TList * toFit= (TList*)gMinuit->GetObjectFit() ;
  TArrayD * params=(TArrayD*)toFit->At(0) ; 
  TArrayI * samples = (TArrayI*)toFit->At(1) ;

  fret = 0. ;     
  if(iflag == 2)
    for(Int_t iparam = 0 ; iparam < 2 ; iparam++)    
      Grad[iparam] = 0 ; // Will evaluate gradient
  
  Int_t nSamples=samples->GetSize() ; //Math::Min(70,samples->GetSize()) ;
  Double_t t0=x[0] ;
  Double_t en=x[1] ;
  Double_t a=params->At(0) ;
  Double_t alpha=params->At(1) ;
  Double_t b=params->At(2) ;
  Double_t beta=params->At(3) ;
  
  for(Int_t i = 0 ; i < nSamples ; i++) {
    if(samples->At(i)==0 || samples->At(i)==1023) //zero or overflow
      continue ;
    Double_t dt=i*1.-t0 ;
    Double_t diff=float(samples->At(i))-Gamma2(dt,en,params) ;
    Double_t w=0.1+0.005*i ; //Mean Pedestal RMS + rising modulation
    //    if(w==0)w=1. ;
    diff/=w ;
    if(iflag == 2){  // calculate gradient
      if(dt>=0.){
	Grad[0] += -2.*en*diff*((alpha*dt*(1.+a*dt)-1.-2.*a*dt)*TMath::Exp(-alpha*dt)+
                              b*(beta*dt*(1.+a*dt)-1.-2.*a*dt)*TMath::Exp(-beta*dt)) /w ; //derivative over t0
	Grad[1] += -2.*diff*(dt+a*dt*dt)*(TMath::Exp(-alpha*dt)+
                     b*TMath::Exp(-dt*beta))/w ;
      }
    }
    fret += diff*diff ;
  }
  if(nSamples){
    fret/=nSamples ;
    if(iflag == 2){
      for(Int_t iparam = 0 ; iparam < 2 ; iparam++)    
	Grad[iparam] /= nSamples ;
    }
  }
  
}
//-----------------------------------------------------------------------------
Double_t AliPHOSRawDecoderv1::Gamma2(Double_t dt,Double_t en,TArrayD * params){  //Function for fitting samples
  //parameters:
  //dt-time after start
  //en-amplutude
  //function parameters
  
  Double_t ped=params->At(4) ;
  if(dt<0.)
    return ped ; //pedestal
  else
    return ped+en*(dt+params->At(0)*dt*dt)*(TMath::Exp(-dt*params->At(1))+params->At(2)*TMath::Exp(-dt*params->At(3))) ;
}

