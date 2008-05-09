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

// Author: Dmitri Peressounko using idea of Y.Kucheryaev

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
#include "AliPHOSRawDecoderv2.h"
#include "AliPHOSPulseGenerator.h"


ClassImp(AliPHOSRawDecoderv2)

//-----------------------------------------------------------------------------
  AliPHOSRawDecoderv2::AliPHOSRawDecoderv2():AliPHOSRawDecoder(),
fNtimeSamples(25),fRMScut(11.)
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
AliPHOSRawDecoderv2::AliPHOSRawDecoderv2(AliRawReader* rawReader,  AliAltroMapping **mapping):
  AliPHOSRawDecoder(rawReader,mapping),
fNtimeSamples(25),fRMScut(11.)
{
  //Construct a decoder object.
  //Is is user responsibility to provide next raw event 
  //using AliRawReader::NextEvent().
  fLGpar[0]=0.971 ;
  fLGpar[1]=0.0465;
  fLGpar[2]=1.56  ;
  fHGpar[0]=0.941 ; 
  fHGpar[1]=0.0436;
  fHGpar[2]=1.65  ;
}

//-----------------------------------------------------------------------------
AliPHOSRawDecoderv2::~AliPHOSRawDecoderv2()
{
  //Destructor.
  //Nothing to delete
}

//-----------------------------------------------------------------------------
AliPHOSRawDecoderv2::AliPHOSRawDecoderv2(const AliPHOSRawDecoderv2 &phosDecoder ):
  AliPHOSRawDecoder(phosDecoder), 
fNtimeSamples(25),fRMScut(11.)
{
  //Copy constructor.
  fNtimeSamples=phosDecoder.fNtimeSamples ;
  for(Int_t i=0; i<3;i++){
    fLGpar[i]=phosDecoder.fLGpar[i] ;
    fHGpar[i]=phosDecoder.fHGpar[i] ;
  }
  fRMScut=phosDecoder.fRMScut ;
}

//-----------------------------------------------------------------------------
AliPHOSRawDecoderv2& AliPHOSRawDecoderv2::operator = (const AliPHOSRawDecoderv2 &phosDecoder)
{
  //Assignment operator.

  fNtimeSamples=phosDecoder.fNtimeSamples ;
  for(Int_t i=0; i<3;i++){
    fLGpar[i]=phosDecoder.fLGpar[i] ;
    fHGpar[i]=phosDecoder.fHGpar[i] ;
  }
  fRMScut=phosDecoder.fRMScut ;
  return *this;
}

//-----------------------------------------------------------------------------
Bool_t AliPHOSRawDecoderv2::NextDigit()
{
  //Extract an energy deposited in the crystal,
  //crystal' position (module,column,row),
  //time and gain (high or low).
  //First collects sample, then evaluates it and if it has
  //reasonable shape, fits it with Gamma2 function and extracts 
  //energy and time.

  TCanvas * cs = (TCanvas*)gROOT->FindObjectAny("CSample") ;
  if(!cs)
    cs = new TCanvas("CSample","CSample") ;

  TH1D * h = (TH1D*)gROOT->FindObjectAny("hSample") ;
  if(!h)
    h=new TH1D("hSample","",200,0.,200.) ;


  AliCaloRawStream* in = fCaloStream;
  
  Int_t    iBin     = fSamples->GetSize();
  fEnergy = 0;
  Double_t pedMean = 0;
  Double_t pedRMS = 0;
  Int_t   nPed = 0;
  Double_t baseLine = 1.0;
  const Int_t nPreSamples = 10;
  fQuality = 0. ;
  
  while ( in->Next() ) { 
    
    // Fit the full sample
    if(in->IsNewHWAddress() && iBin!=fSamples->GetSize()) {
      
      Double_t pedestal =0. ;
      if(fPedSubtract){ 
	if (nPed > 0)
	  pedestal = (Double_t)(pedMean/nPed); 
	else
	  return kFALSE;
      }
      for(Int_t i=0; i<fSamples->GetSize(); i++){
        h->SetBinContent(i+1,fSamples->At(i)) ;
      }      

      //calculate time and energy
      Int_t maxBin=0 ;
      Int_t maxAmp=0 ; 
      for(Int_t i=iBin; i<fSamples->GetSize(); i++){
        if(maxAmp<fSamples->At(i)){
          maxBin=i ;
          maxAmp=fSamples->At(i) ;
        }
      }
      if(maxBin==fSamples->GetSize()-1){//bad sample 
        fEnergy=0. ;                                                                                                                       
        fTime=-999.;                                                                                                                       
        return kTRUE ;                                                                                                                     
      } 
      fEnergy=Double_t(maxAmp)-pedestal ;
      fOverflow =0 ;  //look for plato on the top of sample
      if(fEnergy>500 &&  //this is not fluctuation of soft sample
         maxBin<fSamples->GetSize()-1 && fSamples->At(maxBin+1)==maxAmp){ //and there is a plato
         fOverflow = kTRUE ;
      }

//    if(fEnergy>500.){
// if(fRow==54 && fColumn==24){
//    printf("fE=%f, ped=%f, row=%d, col=%d \n",fEnergy,pedestal,fRow,fColumn) ;
//    if(fOverflow)printf(" Overflow \n") ;
//    else printf("iBin=%d, maxBin=%d, maxAmp=%d,Amp(+1)=%d,Amp(-1)=%d  \n",iBin,maxBin,maxAmp,fSamples->At(maxBin+1),fSamples->At(maxBin-1)) ;
//    cs->cd() ;
//    h->Draw() ;
//    cs->Update() ;
//    getchar() ;
// }

      if(fOverflow)
        return kTRUE ; //do not calculate energy and time for overflowed channels

      if(fEnergy<baseLine){ //do not evaluate time, drop this sample
        fEnergy=0. ;
        fTime=-999.;
        return kTRUE ;
      }

      //else calculate time
      fTime=0. ;
      Double_t tRMS = 0. ;
      Double_t tW = 0. ;
      Int_t cnts=0 ;
      Double_t a,b,c ;
      if(fLowGainFlag){
        a=fLGpar[0] ; 
        b=fLGpar[1] ; 
        c=fLGpar[2] ; 
      }
      else{
        a=fHGpar[0] ; 
        b=fHGpar[1] ; 
        c=fHGpar[2] ; 
      }

      for(Int_t i=iBin+1; i<fSamples->GetSize()&& cnts<fNtimeSamples; i++){
	if(fSamples->At(i)<pedestal)
          continue ;
//Presently we do not use channels with overflow
//        if(fOverflow && (fSamples->At(i)==maxAmp ||fSamples->At(i-1)==maxAmp)) //can not calculate time in overflow bin
//          continue ;
        if(fTimes->At(i)-fTimes->At(i-1)!=1) //do not use samples with non-consequtive points
          continue ;
	Double_t de=fSamples->At(i)-pedestal ;
        Double_t av = de+fSamples->At(i-1)-pedestal ;
        if(av<=0.) //this is fluctuation around pedestal, scip
          continue ;
        Double_t ds = fSamples->At(i)-fSamples->At(i-1) ;
        Double_t ti = ds/av ;     // calculate log. derivative
        ti=a/(ti+b)-c*ti ;        // and compare with parameterization
        ti=Double_t(fTimes->At(i))-ti ; 
        Double_t wi = TMath::Abs(ds) ;
        fTime+=ti*wi ;
        tW+=wi;
        tRMS+=ti*ti*wi ;
        cnts++ ;
      } 
 
      if(tW>0.){
        fTime/=tW ;
        fQuality = tRMS/tW-fTime*fTime ;
        //Normalize quality
//printf("t0=%f, RMS=%f, cut=%f \n",fTime,tRMS,fRMScut) ;
//        if(tRMS>=fRMScut){ //bad sample
//          fTime=-999. ;
//          fEnergy=0. ;
//        }
      }
      else{
        fTime=-999. ;
        fQuality=999. ;
      }

      Bool_t isBad = 0 ;
      for(Int_t i=iBin+1; i<fSamples->GetSize()-1&&!isBad; i++){
        if(fSamples->At(i)>fSamples->At(i-1)+5 && fSamples->At(i)>fSamples->At(i+1)+5) { //single jump
          isBad=1 ;
        }
      }
      if(pedestal<10.)
        isBad=1 ;

      pedRMS=pedRMS/nPed-pedestal*pedestal ;
      if(pedRMS>0.1)
        isBad=1 ;

      for(Int_t i=iBin+1; i<fSamples->GetSize()-1&&!isBad; i++){                                                                           
         if(fSamples->At(i)<pedestal-1)
           isBad=1 ;
      }

      //two maxima


    if(fEnergy>10. && !isBad ){
    printf("fE=%f, ped=%f, fQuality=%f, pedRMS=%f \n",fEnergy,pedestal,fQuality,pedRMS) ;
    if(fOverflow)printf(" Overflow \n") ;
    if(isBad)printf("bad") ;
//    else printf("iBin=%d, maxBin=%d, maxAmp=%d,Amp(+1)=%d,Amp(-1)=%d  \n",iBin,maxBin,maxAmp,fSamples->At(maxBin+1),fSamples->At(maxBin-1)) ;
    cs->cd() ;
    h->Draw() ;
    cs->Update() ;
    getchar() ;
 }


      return kTRUE ; //will scan further
    }
    
    fLowGainFlag = in->IsLowGain();
    fModule = in->GetModule()+1;
    fRow    = in->GetRow()   +1;
    fColumn = in->GetColumn()+1;
    
    
    // Fill array with samples
    iBin--;                                                             
    if(fPedSubtract && (in->GetTime() < nPreSamples)) {
      pedMean += in->GetSignal();
      pedRMS += in->GetSignal()*in->GetSignal();
      nPed++;
    }
    fSamples->AddAt(in->GetSignal(),iBin);
    fTimes->AddAt(in->GetTime(),iBin);
  } // in.Next()
  
  return kFALSE;
}
