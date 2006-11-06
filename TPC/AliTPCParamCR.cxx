/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

///////////////////////////////////////////////////////////////////////
//  Manager and of geomety  classes for set: TPC                     //
//                                                                   //
//  !sectors are numbered from  0                                     //
//  !pad rows are numbered from 0                                     //
//  
//  27.7.   - AliTPCPaaramSr object for TPC 
//            TPC with straight pad rows 
//  Origin:  Marian Ivanov, Uni. of Bratislava, ivanov@fmph.uniba.sk // 
//                                                                   //  
///////////////////////////////////////////////////////////////////////


#include <Riostream.h>
#include <TMath.h>
#include <TClass.h>
#include <AliTPCParamCR.h>
#include "AliTPCPRF2D.h"
#include "AliTPCRF1D.h"



ClassImp(AliTPCParamCR)
static const  Int_t kMaxRows=600;
static const  Float_t  kEdgeSectorSpace = 2.5;

AliTPCParamCR::AliTPCParamCR()
              :AliTPCParam(),
	       fInnerPRF(0),
	       fOuter1PRF(0),
	       fOuter2PRF(0),
               fTimeRF(0),
	       fFacSigma(0.)
{   
  //
  //constructor set the default parameters

  fFacSigma = Float_t(2.);
  SetDefault();
  Update();
}
AliTPCParamCR::AliTPCParamCR(const AliTPCParamCR &param)
              :AliTPCParam(),
	       fInnerPRF(0),
	       fOuter1PRF(0),
	       fOuter2PRF(0),
               fTimeRF(0),
	       fFacSigma(0.)
{
  //
  // copy constructor - dummy
  //
  fFacSigma= param.fFacSigma;
}
AliTPCParamCR & AliTPCParamCR::operator =(const AliTPCParamCR & param)
{
  //
  // assignment operator - dummy
  //
   fFacSigma= param.fFacSigma;
  return (*this); 
}
AliTPCParamCR::~AliTPCParamCR()
{
  //
  //destructor destroy some dynmicaly alocated variables
  if (fInnerPRF != 0) delete fInnerPRF;
  if (fOuter1PRF != 0) delete fOuter1PRF;
  if (fOuter2PRF != 0) delete fOuter2PRF;
  if (fTimeRF != 0) delete fTimeRF;
}

void AliTPCParamCR::SetDefault()
{
  //set default TPC param   
  fbStatus = kFALSE;
  AliTPCParam::SetDefault();  
}  

Int_t  AliTPCParamCR::CalcResponse(Float_t* xyz, Int_t * index, Int_t /*dummy*/)
{
  //
  //calculate bin response as function of the input position -x 
  //return number of valid response bin
  //
  //we suppose that coordinata is expressed in float digits 
  // it's mean coordinate system 8
  //xyz[0] - float padrow xyz[1] is float pad  (center pad is number 0) and xyz[2] is float time bin
  if ( (fInnerPRF==0)||(fOuter1PRF==0) ||(fOuter2PRF==0)||(fTimeRF==0) ){ 
    Error("AliTPCParamCR", "response function were not adjusted");
    return -1;
  }
  
  Float_t sfpadrow;   // sigma of response function
  Float_t sfpad;      // sigma  of 
  Float_t sftime= fFacSigma*fTimeRF->GetSigma()/fZWidth;     //3 sigma of time response
  if (index[1]<fNInnerSector){
    sfpadrow =fFacSigma*fInnerPRF->GetSigmaY()/fInnerPadPitchLength;
    sfpad    =fFacSigma*fInnerPRF->GetSigmaX()/fInnerPadPitchWidth;
  }else{
    if(index[2]<fNRowUp1){
      sfpadrow =fFacSigma*fOuter1PRF->GetSigmaY()/fOuter1PadPitchLength;
      sfpad    =fFacSigma*fOuter1PRF->GetSigmaX()/fOuterPadPitchWidth;}
    else{ sfpadrow =fFacSigma*fOuter2PRF->GetSigmaY()/fOuter2PadPitchLength;
    sfpad    =fFacSigma*fOuter2PRF->GetSigmaX()/fOuterPadPitchWidth;}
  }

  Int_t fpadrow = TMath::Nint(xyz[0]-sfpadrow);  //"first" padrow
  Int_t fpad    = TMath::Nint(xyz[1]-sfpad);     //first pad
  Int_t ftime   = TMath::Nint(xyz[2]+fTimeRF->GetOffset()-sftime);    // first time
  Int_t lpadrow = TMath::Min(TMath::Nint(xyz[0]+sfpadrow),fpadrow+19);  //"last" padrow
  Int_t lpad    = TMath::Min(TMath::Nint(xyz[1]+sfpad),fpad+19);     //last pad
  Int_t ltime   = TMath::Min(TMath::Nint(xyz[2]+fTimeRF->GetOffset()+sftime),ftime+19);    // last time
   
  Float_t  padres[20][20];  //I don't expect bigger number of bins
  Float_t  timeres[20];     
  //calculate padresponse function 
  Int_t padrow; 
  for (padrow = fpadrow;padrow<=lpadrow;padrow++)
    for (Int_t pad = fpad;pad<=lpad;pad++){
      Float_t dy = (xyz[0]-Float_t(padrow));
      Float_t dx = (xyz[1]-Float_t(pad));
      if (index[1]<fNInnerSector)
	padres[padrow-fpadrow][pad-fpad]=fInnerPRF->GetPRF(dx*fInnerPadPitchWidth,dy*fInnerPadPitchLength);
      else{
	if(index[2]<fNRowUp1){
	  padres[padrow-fpadrow][pad-fpad]=fOuter1PRF->GetPRF(dx*fOuterPadPitchWidth,dy*fOuter1PadPitchLength);}
	else{ padres[padrow-fpadrow][pad-fpad]=fOuter2PRF->GetPRF(dx*fOuterPadPitchWidth,dy*fOuter2PadPitchLength);}}}     
    
  //calculate time response function

  Int_t time;
  for (time = ftime;time<=ltime;time++) timeres[time-ftime]= fTimeRF->GetRF((xyz[2]-Float_t(time))*fZWidth); 
    
  //write over threshold values to stack
  Int_t cindex3=-1;
  Int_t cindex=0;
  Float_t cweight = 0;
  for (padrow = fpadrow;padrow<=lpadrow;padrow++)
    for (Int_t pad = fpad;pad<=lpad;pad++)
      for (time = ftime;time<=ltime;time++){
	cweight = timeres[time-ftime]*padres[padrow-fpadrow][pad-fpad];
	if (cweight>fResponseThreshold) {
	  fResponseBin[++cindex3]=padrow;
	  fResponseBin[++cindex3]=pad;
	  fResponseBin[++cindex3]=time;
	  fResponseWeight[++cindex]=cweight;
	}
      }
  fCurrentMax=cindex;	
  return fCurrentMax;
}

void AliTPCParamCR::CRXYZtoXYZ(Float_t *xyz,
	       const Int_t &sector, const Int_t & padrow, Int_t option) const  
{  
  //transform relative coordinates to absolute
  Bool_t rel = ( (option&2)!=0);
  Int_t index[2]={sector,padrow};
  if (rel==kTRUE)      Transform4to3(xyz,index);//if the position is relative to pad row  
  Transform2to1(xyz,index);
}

void AliTPCParamCR::XYZtoCRXYZ(Float_t *xyz,
			     Int_t &sector, Int_t & padrow, Int_t option) const
{
   //transform global position to the position relative to the sector padrow
  //if option=0  X calculate absolute            calculate sector
  //if option=1  X           absolute            use input sector
  //if option=2  X           relative to pad row calculate sector
  //if option=3  X           relative            use input sector
  //!!!!!!!!! WE start to calculate rows from row = 0
  Int_t index[2];
  Bool_t rel = ( (option&2)!=0);  

  //option 0 and 2  means that we don't have information about sector
  if ((option&1)==0)   Transform0to1(xyz,index);  //we calculate sector number 
  else
    index[0]=sector;
  Transform1to2(xyz,index);
  Transform2to3(xyz,index);
  //if we store relative position calculate position relative to pad row
  if (rel==kTRUE) Transform3to4(xyz,index);
  sector = index[0];
  padrow = index[1];
}


         
Bool_t AliTPCParamCR::Update()
{
  Int_t i;
  if (AliTPCParam::Update()==kFALSE) return kFALSE;
  fbStatus = kFALSE;

 Float_t firstrow = fInnerRadiusLow + 2.225 ;   
 for( i= 0;i<fNRowLow;i++)
   {
     Float_t x = firstrow + fInnerPadPitchLength*(Float_t)i;  
     fPadRowLow[i]=x;
     // number of pads per row
/*Float_t y = (x-0.5*fInnerPadPitchLength)*tan(fInnerAngle/2.)-fInnerWireMount-
	fInnerPadPitchWidth/2.;*/
     Float_t y = x*tan(fInnerAngle/2.)-fInnerWireMount;
     fYInner[i]=y;
     fNPadsLow[i] = 1+2*(Int_t)(y/fInnerPadPitchWidth) ;

   }
 firstrow = fOuterRadiusLow + 1.6;
 for(i=0;i<fNRowUp;i++)
   {
     if(i<fNRowUp1){
       Float_t x = firstrow + fOuter1PadPitchLength*(Float_t)i; 
       fPadRowUp[i]=x;
/*Float_t y =(x-0.5*fOuter1PadPitchLength)*tan(fOuterAngle/2.)-fOuterWireMount-
  fOuterPadPitchWidth/2.;*/
       Float_t y = x*tan(fInnerAngle/2.)-fInnerWireMount;
       fNPadsUp[i] = 1+2*(Int_t)(y/fOuterPadPitchWidth) ;
       fYOuter[i] = y;   
       if(i==fNRowUp1-1) {
           fLastWireUp1=fPadRowUp[i] +0.375;
           firstrow = fPadRowUp[i] + 0.5*(fOuter1PadPitchLength+fOuter2PadPitchLength);
       }
     }
     else
       {
	 Float_t x = firstrow + fOuter2PadPitchLength*(Float_t)(i-64);
/*Float_t y =(x-0.5*fOuter2PadPitchLength)*tan(fOuterAngle/2.)-fOuterWireMount-
  fOuterPadPitchWidth/2.;*/
	 Float_t y = x*tan(fInnerAngle/2.)-fInnerWireMount;
         fNPadsUp[i] = 1+2*(Int_t)(y/fOuterPadPitchWidth) ;
	 fYOuter[i] = y;
       }
   }   
     
  fNtRows = fNInnerSector*fNRowLow+fNOuterSector*fNRowUp;
  fbStatus = kTRUE;
  return kTRUE;
}



void AliTPCParamCR::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliTPC.

   if (R__b.IsReading()) {
      AliTPCParamCR::Class()->ReadBuffer(R__b, this);
      Update();
   } else {
      AliTPCParamCR::Class()->WriteBuffer(R__b, this);
   }
}










