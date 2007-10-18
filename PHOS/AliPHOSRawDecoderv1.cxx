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
#include "TMath.h"
#include "TMinuit.h"

// --- AliRoot header files ---
#include "AliPHOSRawDecoderv1.h"
#include "AliPHOSPulseGenerator.h"


ClassImp(AliPHOSRawDecoderv1)

//-----------------------------------------------------------------------------
  AliPHOSRawDecoderv1::AliPHOSRawDecoderv1():AliPHOSRawDecoder()
{
  //Default constructor.
}

//-----------------------------------------------------------------------------
AliPHOSRawDecoderv1::AliPHOSRawDecoderv1(AliRawReader* rawReader,  AliAltroMapping **mapping):
  AliPHOSRawDecoder(rawReader,mapping)
{
  //Construct a decoder object.
  //Is is user responsibility to provide next raw event 
  //using AliRawReader::NextEvent().

  if(!gMinuit) 
    gMinuit = new TMinuit(100);

}

//-----------------------------------------------------------------------------
AliPHOSRawDecoderv1::~AliPHOSRawDecoderv1()
{
  //Destructor.
}

//-----------------------------------------------------------------------------
AliPHOSRawDecoderv1::AliPHOSRawDecoderv1(const AliPHOSRawDecoderv1 &phosDecoder ):
  AliPHOSRawDecoder(phosDecoder)
{
  //Copy constructor.
}

//-----------------------------------------------------------------------------
AliPHOSRawDecoderv1& AliPHOSRawDecoderv1::operator = (const AliPHOSRawDecoderv1 &phosDecode)
{
  //Assignment operator.

  //  if(this != &phosDecode) {
  //  }

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
      
      if(fPedSubtract){ 
	fEnergy-=pedestal ;
      }
      
      //IF sample has reasonable mean and RMS, try to fit it with gamma2
      if(fEnergy>2.&& cnts >20 && aMean>0. && aRMS>2.){ //more or less reasonable sample
	
	gMinuit->mncler();                     // Reset Minuit's list of paramters
	gMinuit->SetPrintLevel(-1) ;           // No Printout
	gMinuit->SetFCN(AliPHOSRawDecoderv1::UnfoldingChiSquare) ;  
	// To set the address of the minimization function 
	
	gMinuit->SetObjectFit((TObject*)fSamples) ;         // To tranfer pointer to UnfoldingChiSquare
	Int_t ierflg ;
	gMinuit->mnparm(0, "p",  pedestal, 0.1, 0, 0, ierflg) ;
	if(ierflg != 0){ 
	  //	  Warning("NextDigit", "Uunable to set initial value for fit procedure : ped = %f\n", pedestal ) ;
	  fEnergy=0. ;
	  fTime=-999. ;
	  return kTRUE ; //will scan further
	}
	gMinuit->mnparm(1, "t0",  1.*tStart, 0.1, 0, 0, ierflg) ;
	if(ierflg != 0){
	  //	  Warning("NextDigit", "Unable to set initial value for fit procedure : t0\n" ) ;
	  fEnergy=0. ;
	  fTime=-999. ;
	  return kTRUE ; //will scan further
	}
	gMinuit->mnparm(2, "Energy", fEnergy*0.018 , 0.001*fEnergy, 0, 0, ierflg) ;
	if(ierflg != 0){
	  //	  Warning("NextDigit", "Unable to set initial value for fit procedure : E=%f\n", fEnergy*0.018) ;
	  fEnergy=0. ;
	  fTime=-999. ;
	  return kTRUE ; //will scan further
	}
	gMinuit->mnparm(3, "Slope", 0.09 , 0.001, 0.001, 1., ierflg) ;
	if(ierflg != 0){
	  //	  Warning("NextDigit", "Unable to set initial value for fit procedure : Slope\n") ;
	  fEnergy=0. ;
	  fTime=-999. ;
	  return kTRUE ; //will scan further
	}
	
	Double_t p0 = 1. ; // "Tolerance" Evaluation stops when EDM = 0.0001*p0 ; The number of function call slightly
	//  depends on it. 
	Double_t p1 = 1.0 ;
	Double_t p2 = 0.0 ;
	gMinuit->mnexcm("SET STR", &p2, 0, ierflg) ;   // force TMinuit to reduce function calls  
	gMinuit->mnexcm("SET GRA", &p1, 1, ierflg) ;   // force TMinuit to use my gradient  
	//	gMinuit->SetMaxIterations(100);
	gMinuit->mnexcm("SET NOW", &p2 , 0, ierflg) ;  // No Warnings
	
	gMinuit->mnexcm("MIGRAD", &p0, 0, ierflg) ;    // minimize 
	
	Double_t err ;
	Double_t p,t0,efit,slope ;
	gMinuit->GetParameter(0,p, err) ;    
	gMinuit->GetParameter(1,t0, err) ;    
	gMinuit->GetParameter(2,efit, err) ;    
	gMinuit->GetParameter(3,slope, err) ;    
	
	efit=efit*4.*TMath::Exp(-2.)/slope/slope ; //slope can not be zero - par limits
	
	Double_t fmin,fedm,errdef ;
	Int_t npari,nparx,istat;
	gMinuit->mnstat(fmin,fedm,errdef,npari,nparx,istat) ;
	
	if(fmin < cnts*(efit/10.) ){ //Chi^2 of a good sample
	  if(efit<1050.&& efit>0.){
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
    
  } // in.Next()
  
  return kFALSE;
}
//_____________________________________________________________________________
void AliPHOSRawDecoderv1::UnfoldingChiSquare(Int_t & /*nPar*/, Double_t * Grad, Double_t & fret, Double_t * x, Int_t iflag)
{
  // Calculates the Chi square for the samples minimization
  // Number of parameters, Gradient, Chi squared, parameters, what to do

  TArrayI * samples = (TArrayI*)gMinuit->GetObjectFit() ;

  fret = 0. ;     
  if(iflag == 2)
    for(Int_t iparam = 0 ; iparam < 4 ; iparam++)    
      Grad[iparam] = 0 ; // Will evaluate gradient
  
  Int_t nSamples=samples->GetSize() ;
  Double_t p =x[0] ;
  Double_t t0=x[1] ;
  Double_t en=x[2] ;
  Double_t a =x[3] ;
  
  for(Int_t i = 0 ; i < nSamples ; i++) {
    if(samples->At(i)==0)
      continue ;
    Double_t dt=i*1.-t0 ;
    Double_t diff=samples->At(i)-Gamma2(dt,p,en,a) ;
    Double_t w=1. ; //TMath::Ceil(TMath::Abs(samples->At(i)-ped)/10.) ;
    //    if(w==0)w=1. ;
    diff/=w ;
    if(iflag == 2){  // calculate gradient
      Grad[0] += -2.*diff/w ; //derivative over pedestal
      if(dt>=0.){
	Grad[1] += -2.*en*diff*(a*dt-2.*dt)*dt*TMath::Exp(-a*dt)/w ; //derivative over t0
	Grad[2] += -2.*diff*dt*dt*TMath::Exp(-a*dt)/w ;
	Grad[3] +=  2.*en*diff*dt*dt*dt*TMath::Exp(-a*dt)/w ;
      }
    }
    fret += diff*diff ;
  }
  /*
  if(nSamples){
    fret/=nSamples ;
    if(iflag == 2){
      for(Int_t iparam = 0 ; iparam < 4 ; iparam++)    
	Grad[iparam] /= nSamples ;
    }
  }
  */
  
  
}
//-----------------------------------------------------------------------------
Double_t AliPHOSRawDecoderv1::Gamma2(Double_t dt,Double_t p,Double_t en,Double_t a){
  //Function for fitting samples
  //parameters:
  //p-pedestal
  //en-amplutude
  //a-decay time

  if(dt<0.)
    return p ; //pedestal
  else
    return p+en*dt*dt*TMath::Exp(-dt*a) ;
}

