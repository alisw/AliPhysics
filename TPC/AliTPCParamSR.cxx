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

/*
$Log$
Revision 1.1.4.2  2000/04/10 11:36:13  kowal2

New Detector parameters handling class

*/

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


#include <iostream.h>
#include <TMath.h>
#include <TObject.h>
#include <AliTPCParamSR.h>



ClassImp(AliTPCParamSR)
const static  Int_t kMaxRows=600;
const static  Float_t  kEdgeSectorSpace = 2.5;
const static Float_t kFacSigmaPadRow=2.;
const static Float_t kFacSigmaPad=3.;
const static Float_t kFacSigmaTime=3.;


AliTPCParamSR::AliTPCParamSR()
{   
  //
  //constructor set the default parameters
  fInnerPRF=0;
  fOuterPRF=0;
  fTimeRF = 0;
  fFacSigmaPadRow = Float_t(kFacSigmaPadRow);
  fFacSigmaPad = Float_t(kFacSigmaPad);
  fFacSigmaTime = Float_t(kFacSigmaTime);


  SetDefault();
  Update();
}

AliTPCParamSR::~AliTPCParamSR()
{
  //
  //destructor destroy some dynmicaly alocated variables
  if (fInnerPRF != 0) delete fInnerPRF;
  if (fOuterPRF != 0) delete fOuterPRF;
  if (fTimeRF != 0) delete fTimeRF;
}

void AliTPCParamSR::SetDefault()
{
  //set default TPC param   
  fbStatus = kFALSE;
  AliTPCParam::SetDefault();  
}  

Int_t  AliTPCParamSR::CalcResponse(Float_t* xyz, Int_t * index, Int_t row)
{
  //
  //calculate bin response as function of the input position -x 
  //return number of valid response bin
  //
  //we suppose that coordinate is expressed in float digits 
  // it's mean coordinate system 8
  //xyz[0] - float padrow xyz[1] is float pad  (center pad is number 0) and xyz[2] is float time bin
  if ( (fInnerPRF==0)||(fOuterPRF==0)||(fTimeRF==0) ){ 
    Error("AliTPCParamSR", "response function was not adjusted");
    return -1;
  }
  
  Float_t sfpadrow;   // sigma of response function
  Float_t sfpad;      // sigma  of 
  Float_t sftime= fFacSigmaTime*fTimeRF->GetSigma()/fZWidth;     //3 sigma of time response
  if (index[1]<fNInnerSector){
    sfpadrow =fFacSigmaPadRow*fInnerPRF->GetSigmaY()/fInnerPadPitchLength;
    sfpad    =fFacSigmaPad*fInnerPRF->GetSigmaX()/fInnerPadPitchWidth;
  }else{
    sfpadrow =fFacSigmaPadRow*fOuterPRF->GetSigmaY()/fOuterPadPitchLength;
    sfpad    =fFacSigmaPad*fOuterPRF->GetSigmaX()/fOuterPadPitchWidth;
  }

  Int_t fpadrow = TMath::Max(TMath::Nint(index[2]+xyz[0]-sfpadrow),0);  //"first" padrow
  Int_t fpad    = TMath::Nint(xyz[1]-sfpad);     //first pad
  Int_t ftime   = TMath::Max(TMath::Nint(xyz[2]+GetZOffset()/GetZWidth()-sftime),0);  // first time
  Int_t lpadrow = TMath::Min(TMath::Nint(index[2]+xyz[0]+sfpadrow),fpadrow+19);  //"last" padrow
  lpadrow       = TMath::Min(GetNRow(index[1])-1,lpadrow);
  Int_t lpad    = TMath::Min(TMath::Nint(xyz[1]+sfpad),fpad+19);     //last pad
  Int_t ltime   = TMath::Min(TMath::Nint(xyz[2]+GetZOffset()/GetZWidth()+sftime),ftime+19);    // last time
  ltime         = TMath::Min(ltime,GetMaxTBin()-1); 
 
  if (row>=0) { //if we are interesting about given pad row
    if (fpadrow<=row) fpadrow =row;
    else 
      return 0;
    if (lpadrow>=row) lpadrow = row;
    else 
      return 0;
  }

 
  Float_t  padres[20][20];  //I don't expect bigger number of bins
  Float_t  timeres[20];     
  Int_t cindex3=0;
  Int_t cindex=0;
  Float_t cweight = 0;
  if (fpadrow>=0) {
  //calculate padresponse function    
  Int_t padrow, pad;
  for (padrow = fpadrow;padrow<=lpadrow;padrow++)
    for (pad = fpad;pad<=lpad;pad++){
      Float_t dy = (-xyz[0]+Float_t(index[2]-padrow));
      Float_t dx = (-xyz[1]+Float_t(pad));
      if (index[1]<fNInnerSector)
	padres[padrow-fpadrow][pad-fpad]=fInnerPRF->GetPRF(dx*fInnerPadPitchWidth,dy*fInnerPadPitchLength);
      else
	padres[padrow-fpadrow][pad-fpad]=fOuterPRF->GetPRF(dx*fOuterPadPitchWidth,dy*fOuterPadPitchLength);          }
  //calculate time response function
  Int_t time;
  for (time = ftime;time<=ltime;time++) 
    timeres[time-ftime]= fTimeRF->GetRF((-xyz[2]+Float_t(time))*fZWidth);     
  //write over threshold values to stack
  for (padrow = fpadrow;padrow<=lpadrow;padrow++)
    for (pad = fpad;pad<=lpad;pad++)
      for (time = ftime;time<=ltime;time++){
	cweight = timeres[time-ftime]*padres[padrow-fpadrow][pad-fpad];
	if (cweight>fResponseThreshold) {
	  fResponseBin[cindex3]=padrow;
	  fResponseBin[cindex3+1]=pad;
	  fResponseBin[cindex3+2]=time;
	  cindex3+=3;  
	  fResponseWeight[cindex]=cweight;
	  cindex++;
	}
      }
  }
  fCurrentMax=cindex;	
  return fCurrentMax;
}

void AliTPCParamSR::TransformTo8(Float_t *xyz, Int_t *index) const
{
  //
  // transformate point to digit coordinate
  //
  if (index[0]==0) Transform0to1(xyz,index);
  if (index[0]==1) Transform1to2(xyz,index);
  if (index[0]==2) Transform2to3(xyz,index);
  if (index[0]==3) Transform3to4(xyz,index);
  if (index[0]==4) Transform4to8(xyz,index);
}

void AliTPCParamSR::TransformTo2(Float_t *xyz, Int_t *index) const
{
  //
  //transformate point to rotated coordinate
  //
  //we suppose that   
  if (index[0]==0) Transform0to1(xyz,index);
  if (index[0]==1) Transform1to2(xyz,index);
  if (index[0]==4) Transform4to3(xyz,index);
  if (index[0]==8) {  //if we are in digit coordinate system transform to global
    Transform8to4(xyz,index);
    Transform4to3(xyz,index);  
  }
}

void AliTPCParamSR::CRXYZtoXYZ(Float_t *xyz,
	       const Int_t &sector, const Int_t & padrow, Int_t option) const  
{  
  //transform relative coordinates to absolute
  Bool_t rel = ( (option&2)!=0);
  Int_t index[2]={sector,padrow};
  if (rel==kTRUE)      Transform4to3(xyz,index);//if the position is relative to pad row  
  Transform2to1(xyz,index);
}

void AliTPCParamSR::XYZtoCRXYZ(Float_t *xyz,
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

Float_t AliTPCParamSR::GetPrimaryLoss(Float_t *x, Int_t *index, Float_t *angle)
{
  //
  //
  Float_t padlength=GetPadPitchLength(index[1]);
  Float_t a1=TMath::Sin(angle[0]);
  a1*=a1;
  Float_t a2=TMath::Sin(angle[1]);
  a2*=a2;
  Float_t length =padlength*TMath::Sqrt(1+a1+a2);
  return length*fNPrimLoss;
}

Float_t AliTPCParamSR::GetTotalLoss(Float_t *x, Int_t *index, Float_t *angle)
{
  //
  //
  Float_t padlength=GetPadPitchLength(index[1]);
  Float_t a1=TMath::Sin(angle[0]);
  a1*=a1;
  Float_t a2=TMath::Sin(angle[1]);
  a2*=a2;
  Float_t length =padlength*TMath::Sqrt(1+a1+a2);
  return length*fNTotalLoss;
  
}


void AliTPCParamSR::GetClusterSize(Float_t *x, Int_t *index, Float_t *angle, Int_t mode, Float_t *sigma)
{
  //
  //return cluster sigma2 (x,y) for particle at position x
  // in this case x coordinata is in drift direction
  //and y in pad row direction
  //we suppose that input coordinate system is digit system
   
  Float_t  xx;
  Float_t lx[3] = {x[0],x[1],x[2]};
  Int_t   li[3] = {index[0],index[1],index[2]};
  TransformTo2(lx,li);
  //  Float_t  sigmadiff;
  sigma[0]=0;
  sigma[1]=0;
  
  xx = lx[2];  //calculate drift length in cm
  if (xx>0) {
    sigma[0]+= xx*GetDiffL()*GetDiffL();
    sigma[1]+= xx*GetDiffT()*GetDiffT(); 
  }


  //sigma[0]=sigma[1]=0;
  if (GetTimeRF()!=0) sigma[0]+=GetTimeRF()->GetSigma()*GetTimeRF()->GetSigma();
  if ( (index[1]<fNInnerSector) &&(GetInnerPRF()!=0))   
    sigma[1]+=GetInnerPRF()->GetSigmaX()*GetInnerPRF()->GetSigmaX();
  if ( (index[1]>=fNInnerSector) && (GetOuterPRF()!=0))
    sigma[1]+=GetOuterPRF()->GetSigmaX()*GetOuterPRF()->GetSigmaX();


  sigma[0]/= GetZWidth()*GetZWidth();
  sigma[1]/=GetPadPitchWidth(index[0])*GetPadPitchWidth(index[0]);
}




void AliTPCParamSR::GetSpaceResolution(Float_t *x, Int_t *index, Float_t *angle, 
				       Float_t amplitude, Int_t mode, Float_t *sigma)
{
  //
  //
  //
  
}
Float_t  AliTPCParamSR::GetAmp(Float_t *x, Int_t *index, Float_t *angle)
{
  //
  //
  //
  return 0;
}

Float_t * AliTPCParamSR::GetAnglesAccMomentum(Float_t *x, Int_t * index, Float_t* momentum, Float_t *angle)
{
  //
  //calculate angle of track to padrow at given position
  // for given magnetic field and momentum of the particle
  //

  TransformTo2(x,index);
  AliDetectorParam::GetAnglesAccMomentum(x,index,momentum,angle);    
  Float_t addangle = TMath::ASin(x[1]/GetPadRowRadii(index[1],index[2]));
  angle[1] +=addangle;
  return angle;				 
}

         
Bool_t AliTPCParamSR::Update()
{
  
  //
  // update some calculated parameter which must be updated after changing "base"
  // parameters 
  // for example we can change size of pads and according this recalculate number
  // of pad rows, number of of pads in given row ....
  Int_t i;
  if (AliTPCParam::Update()==kFALSE) return kFALSE;
  fbStatus = kFALSE;

  // adjust lower sectors pad row positions and pad numbers 
  fNRowLow   =  (Int_t(1.001+((fRInnerLastWire-fRInnerFirstWire)/fInnerWWPitch))
	       -2*fInnerDummyWire)/fNInnerWiresPerPad;  
  if ( kMaxRows<fNRowLow) fNRowUp = kMaxRows;
  if (1>fNRowLow) return kFALSE;
 
  //Float_t firstpad = fRInnerFirstWire+(fInnerDummyWire-0.5)*fInnerWWPitch
  //    +fInnerPadPitchLength/2.;
  Float_t lastpad = fRInnerLastWire-(fInnerDummyWire-0.5)*fInnerWWPitch
    -fInnerPadPitchLength/2.;
  Float_t firstpad = lastpad-Float_t(fNRowLow-1)*fInnerPadPitchLength;
  
  for (i = 0;i<fNRowLow;i++) 
    {
       Float_t x  = firstpad +fInnerPadPitchLength*(Float_t)i;       
       Float_t y = (x-0.5*fInnerPadPitchLength)*tan(fInnerAngle/2.)-fInnerWireMount-
	            fInnerPadPitchWidth/2.;
       fPadRowLow[i] = x;
       fNPadsLow[i] = 1+2*(Int_t)(y/fInnerPadPitchWidth) ;
       }

  // adjust upper sectors pad row positions and pad numbers
  fNRowUp   = (Int_t(1.001+((fROuterLastWire-fROuterFirstWire)/fOuterWWPitch))
	       -2*fOuterDummyWire)/fNOuterWiresPerPad; 
  if ( kMaxRows<fNRowUp) fNRowUp = kMaxRows;
  if (1>fNRowUp) return kFALSE;
  firstpad = fROuterFirstWire+(fOuterDummyWire-0.5)*fOuterWWPitch
    +fOuterPadPitchLength/2.;
 
  for (i = 0;i<fNRowUp;i++) 
    {
       Float_t x  = firstpad + fOuterPadPitchLength*(Float_t)i;      
       Float_t y = (x-0.5*fOuterPadPitchLength)*tan(fOuterAngle/2.)-fOuterWireMount-
	            fInnerPadPitchWidth/2.;
       fPadRowUp[i] = x;
       fNPadsUp[i] = 1+2*(Int_t)(y/fOuterPadPitchWidth) ;
    }
  fNtRows = fNInnerSector*fNRowLow+fNOuterSector*fNRowUp;
  return kTRUE;
}



void AliTPCParamSR::Streamer(TBuffer &R__b)
{
   // Stream an object of class AliTPC.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      //      TObject::Streamer(R__b);
      AliTPCParam::Streamer(R__b);
      //      if (R__v < 2) return;
       Update();
   } else {
      R__b.WriteVersion(AliTPCParamSR::IsA());
      //TObject::Streamer(R__b);  
      AliTPCParam::Streamer(R__b);    
   }
}




