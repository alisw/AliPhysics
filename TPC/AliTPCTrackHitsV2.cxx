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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Projection Chamber  track hits object                                //
//
//  Origin: Marian Ivanov , GSI Darmstadt
//
// AliTPCTrackHitsV2
//   Container for Track Hits - based on standard TClonesArray -
//   fArray of AliTPCTrackHitsParamV2 
//   In AliTPCTrackHitsParamV2 - parameterization of the track segment  is stored 
//   for each of the track segment - relative position ( distance between  hits) and
//   charge of the hits is stored - comparing to classical TClonesArray of AliTPChit -
//   comperssion factor of 5-7 (depending on the required precision) -
//   In future release AliTPCTrackHitsV2 - will replace old AliTPCTrackHits - which were not
//   based on standard ROOT containers
//   Basic function:
//      // during building Container
//   AddHitKartez(Int_t volumeID, Int_t trackID, Double_t x, Double_t y, Double_t z,Int_t q)
//   void SetHitPrecision(Double_t prec) {fPrecision=prec;}
//   void SetStepPrecision(Double_t prec) {fStep=prec;}   
//   Bool_t  FlushHitStack(Bool_t force=kTRUE);    
//      //at the end necessary to have Container in consistent state
//    
//     // looping over Container
//   Bool_t  First(), Bool_t Next() - iterators - return status of the operation
//   AliTPChit * GetHit(); - return current hit   


//Begin_Html
/*
<img src="gif/AliTPCTrackHitsV2.gif">
*/
//End_Html
//                                                                           //
//                                                                          //
///////////////////////////////////////////////////////////////////////////////
//

#include "AliTPCTrackHitsV2.h"
#include "TClonesArray.h"    
#include "AliTPC.h"



ClassImp(AliTPCTrackHitsV2) 
ClassImp(AliTrackHitsParamV2)  

  //
Int_t AliTrackHitsParamV2::fgCounter1 =0;
Int_t AliTrackHitsParamV2::fgCounter2 =0;
//
Int_t AliTPCTrackHitsV2::fgCounter1 =0;
Int_t AliTPCTrackHitsV2::fgCounter2 =0;
//
const Double_t AliTPCTrackHitsV2::fgkPrecision=1e-6;  //precision 
const Double_t AliTPCTrackHitsV2::fgkPrecision2=1e-20;  //precision




struct  AliTPCTempHitInfoV2 {
  friend class AliTPCTrackHitsV2;
protected:
  enum    { kStackSize = 10000};
  AliTPCTempHitInfoV2();   
  AliTPCTempHitInfoV2(const AliTPCTempHitInfoV2 &hit)
    {hit.Copy(*this);}
  AliTPCTempHitInfoV2& operator = (const AliTPCTempHitInfoV2 &hit)
     {hit.Copy(*this); return (*this);}
  void     NewParam(Double_t r, Double_t z, Double_t fi, Int_t q);
  void     SetHit(Double_t r, Double_t z, Double_t fi, Int_t q);
  Double_t * GetPosition(Int_t index){return &fPositionStack[index*3];}
  void    UpdateParam(Double_t maxdelta); //recal
  void   Fit2(Double_t fSumY, Double_t fSumYX, Double_t fSumYX2,
	    Double_t fSumX,  Double_t fSumX2, Double_t fSumX3, 
	    Double_t fSumX4, Int_t n,
	      Double_t &a, Double_t &b, Double_t &c);
  void  Fit(AliTrackHitsParamV2 * param);
  Double_t fSumDr;    // Sum of Dr
  Double_t fSumDr2;   // Square of sum of Dr
  Double_t fSumDr3;   // Cube of sum of Dr
  Double_t fSumDr4;   // Fourth power of sum of Dr
  Double_t fSumDFi;  //  Sum of DFi
  Double_t fSumDFiDr; //  Sum of DFiDr
  Double_t fSumDFiDr2;//  Sum of square of DFiDr
  Double_t fSumDZ;     // Sum of DZ
  Double_t fSumDZDr;  //  Sum of DZDr
  Double_t fSumDZDr2;  // Sum of square of DZDr
  Double_t fOldR;     //previos r
  Double_t fPositionStack[3*kStackSize];  //position stack 
  UInt_t   fQStack[kStackSize];           //Q stack
  UInt_t fStackIndex;   //current stack index 
  //  UInt_t fInfoIndex;    //current track info index
  UInt_t fParamIndex;   //current track parameters index
  //  AliTrackHitsInfo  * fInfo; //current track info
  AliTrackHitsParamV2 * fParam; //current track param
private:
  void Copy(AliTPCTempHitInfoV2 &) const
  {printf("Not Implemented\n"); exit(1);}
};


AliTPCTempHitInfoV2::AliTPCTempHitInfoV2()
{
  //
  // Standard constructor
  // set to default value
  //
  fSumDr=fSumDr2=fSumDr3=fSumDr4=
    fSumDFi=fSumDFiDr=fSumDFiDr2=
    fSumDZ=fSumDZDr=fSumDZDr2=0;  
  fStackIndex = 0;
  //  fInfoIndex  = 0;
  fParamIndex = 0;
}


void AliTPCTempHitInfoV2::NewParam(Double_t r, Double_t z, Double_t fi, Int_t q)
{
  //
  //reset stack and sum parameters
  //store line initial point
  //
  fSumDr=fSumDr2=fSumDr3=fSumDr4=
    fSumDFi=fSumDFiDr=fSumDFiDr2=
    fSumDZ=fSumDZDr=fSumDZDr2=0;  
  fStackIndex=0;
  fParam->fR = r;
  fOldR = r;
  fParam->fZ = z;
  fParam->fFi = fi;
  fParam->fAn = 0.;
  fParam->fAd = 0.;
  fParam->fTheta =0.;
  fParam->fThetaD =0.;
  SetHit(r,z,fi,q);
}

void AliTPCTempHitInfoV2::SetHit(Double_t r, Double_t z, Double_t fi, Int_t q)
{
  //
  //add hit to the stack
  //recalculate new estimete of line parameters
  Double_t *f = GetPosition(fStackIndex);  
  f[0] = r;
  f[1] = z;
  f[2] = fi;
  fQStack[fStackIndex]=q;
  if (fStackIndex==0) return;
  Double_t dr  = (r-fParam->fR);
  if (TMath::Abs(dr)<AliTPCTrackHitsV2::fgkPrecision) dr =AliTPCTrackHitsV2::fgkPrecision;
  Double_t dfi = fi-fParam->fFi;
  Double_t dz  = z -fParam->fZ; 
  Double_t dr2 =dr*dr;
  Double_t dr3 =dr2*dr;
  Double_t dr4 =dr3*dr;
  fSumDr +=dr;
  fSumDr2+=dr2;
  fSumDr3+=dr3;
  fSumDr4+=dr4;
  fSumDFi +=dfi;
  fSumDFiDr+=dfi*dr;
  fSumDFiDr2+=dfi*dr2;
  fSumDZ +=dz;
  fSumDZDr+=dz*dr;
  fSumDZDr2+=dz*dr2;
  
  //update fit parameters
  //
  Double_t det = fSumDr2*fSumDr4-fSumDr3*fSumDr3;
  if (TMath::Abs(det)<AliTPCTrackHitsV2::fgkPrecision2) return;
  if ( ( fStackIndex>1 )  ){
    fParam->fAn = (fSumDr4*fSumDFiDr-fSumDr3*fSumDFiDr2)/det;
    fParam->fAd = (fSumDr2*fSumDFiDr2-fSumDr3*fSumDFiDr)/det;
  }
  else
    fParam->fAn = fSumDFiDr/fSumDr2;
  if ( ( fStackIndex>1 )  ){
    fParam->fTheta = (fSumDr4*fSumDZDr-fSumDr3*fSumDZDr2)/det;
    fParam->fThetaD= (fSumDr2*fSumDZDr2-fSumDr3*fSumDZDr)/det;
  }
  else
    fParam->fTheta = fSumDZDr/fSumDr2; 
}


void   AliTPCTempHitInfoV2::UpdateParam(Double_t maxdelta)
{
  //
  // recalc parameters not fixing origin point
  //
  if (fStackIndex>5){ 
    Double_t a,b,c;
    a=b=c=0;
    Fit2(fSumDFi, fSumDFiDr, fSumDFiDr2, fSumDr,fSumDr2,fSumDr3,fSumDr4,
	 fStackIndex, a,b,c);
    if (TMath::Abs(a)<maxdelta){
      fParam->fFi +=a/fParam->fR;    
      fParam->fAn = b;    
      fParam->fAd = c;                  
    }
    Fit2(fSumDZ, fSumDZDr, fSumDZDr2, fSumDr,fSumDr2,fSumDr3,fSumDr4,
	 fStackIndex, a,b,c) ;   
    if (TMath::Abs(a)<maxdelta){
      fParam->fZ +=a;    
      fParam->fTheta = b;    
      fParam->fThetaD = c;   
    }                         
  }
      
}

void   AliTPCTempHitInfoV2::Fit2(Double_t fSumY, Double_t fSumYX, Double_t fSumYX2,
	    Double_t fSumX,  Double_t fSumX2, Double_t fSumX3, 
	    Double_t fSumX4, Int_t n,
	    Double_t &a, Double_t &b, Double_t &c)
{
  //
  // fit of second order
  //
  Double_t det = 
    n* (fSumX2*fSumX4-fSumX3*fSumX3) -
    fSumX*      (fSumX*fSumX4-fSumX3*fSumX2)+
    fSumX2*     (fSumX*fSumX3-fSumX2*fSumX2);
    
  if (TMath::Abs(det)> AliTPCTrackHitsV2::fgkPrecision) {    
    a = 
      (fSumY * (fSumX2*fSumX4-fSumX3*fSumX3)-
       fSumX *(fSumYX*fSumX4-fSumYX2*fSumX3)+
       fSumX2*(fSumYX*fSumX3-fSumYX2*fSumX2))/det; 
    b=
      (n*(fSumYX*fSumX4-fSumX3*fSumYX2)-
      fSumY*(fSumX*fSumX4-fSumX3*fSumX2)+
      fSumX2*(fSumX*fSumYX2-fSumYX*fSumX2))/det;
    c=
      (n*(fSumX2*fSumYX2-fSumYX*fSumX3)-
       fSumX*(fSumX*fSumYX2-fSumYX*fSumX2)+
       fSumY*(fSumX*fSumX3-fSumX2*fSumX2))/det;  
  }
}

void   AliTPCTempHitInfoV2::Fit(AliTrackHitsParamV2 * param)
{
  //
  // fit fixing first and the last point 
  // result stored in new param
  //
  Double_t dx2  = (GetPosition(fStackIndex))[0]-fParam->fR;
  Double_t det = fSumDr4+dx2*fSumDr2-2*dx2*fSumDr3;
  if ( (TMath::Abs(det)> AliTPCTrackHitsV2::fgkPrecision) &&
       ((TMath::Abs(dx2)> AliTPCTrackHitsV2::fgkPrecision))){
    Double_t dfi2 = (GetPosition(fStackIndex))[1]-fParam->fFi;
    param->fAd = (fSumDFiDr2+dfi2*fSumDr-dx2*fSumDFiDr-dfi2*fSumDr3/dx2)/det;
    param->fAn  = (dfi2-param->fAd*dx2*dx2)/dx2;
    
    Double_t dz2 = (GetPosition(fStackIndex))[1]-fParam->fZ;
    param->fTheta = (fSumDZDr2+dz2*fSumDr-dx2*fSumDZDr-dz2*fSumDr3/dx2)/det;
    param->fTheta  = (dz2-param->fAd*dx2*dx2)/dx2;
  }
  
}

AliTrackHitsParamV2::AliTrackHitsParamV2()
{
  //
  // default constructor
  //
  fgCounter1++;
  fgCounter2++;
  fHitDistance=0;
  fCharge=0;
  fNHits=0;
}

AliTrackHitsParamV2::~AliTrackHitsParamV2()
{
  //
  // Standard destructor
  //
  fgCounter1--;
  if (fHitDistance) {
    delete[]fHitDistance;  
    fHitDistance=0;
  }
  if (fCharge){
    delete[]fCharge;  
    fCharge =0;
  }
}


AliTPCTrackHitsV2::AliTPCTrackHitsV2()
{
  //
  //default constructor
  //
  const Float_t kHitPrecision=0.002; //default precision for hit position in cm
  const Float_t kStep =0.003;  //30 mum step 
  const UShort_t kMaxDistance =100;  //maximum distance 100  

  fPrecision=kHitPrecision; //precision in cm
  fStep = kStep; //step size
  fMaxDistance = kMaxDistance; //maximum distance
  fTempInfo =0;
  fSize=0;
  //fTrackHitsInfo = new AliObjectArray("AliTrackHitsInfo"); 
  //fTrackHitsParam = new AliObjectArray("AliTrackHitsParamV2");
  //fHitsPosAndQ = new TArrayOfArrayVStack("AliHitInfo");
  fArray  = new TClonesArray("AliTrackHitsParamV2");
  fCurrentHit = new AliTPCCurrentHitV2;
  fVolumes =0;
  fNVolumes =0;
  fHit =0;
  fgCounter1++;
  fgCounter2++;

} 

AliTPCTrackHitsV2::~AliTPCTrackHitsV2()
{
  //
  //default destructor
  //
  //  if (fTrackHitsInfo) delete fTrackHitsInfo;
  if (fArray) {
    delete fArray;
    fArray =0;
  }
  //if (fHitsPosAndQ) delete fHitsPosAndQ;
  if (fCurrentHit) delete fCurrentHit;
  if (fTempInfo) delete fTempInfo;
  if (fVolumes) {
    delete [] fVolumes;
    fVolumes =0;
    fNVolumes=0;
  }
  if (fHit){
    delete fHit;
    fHit=0;
  }
  fgCounter1--;
}

void AliTPCTrackHitsV2::Clear()
{
  //
  // clear object  
  //
  fSize = 0;
  if (fArray){
    for (Int_t i=0;i<fArray->GetEntriesFast();i++){
      AliTrackHitsParamV2 * par = (AliTrackHitsParamV2 *)fArray->UncheckedAt(i);
      par->~AliTrackHitsParamV2();  // delete object
    }
    fArray->Clear();  
  }
  if (fTempInfo){
    delete fTempInfo; 
    delete fHit;
    fHit =0;
    fTempInfo =0;
  } 
  if (fVolumes){
    delete [] fVolumes;
    fVolumes=0;
    fNVolumes=0;
  }
}


void AliTPCTrackHitsV2::AddHitKartez(Int_t volumeID, Int_t trackID, Double_t x, 
	      Double_t y, Double_t z,Int_t q)
{
  //
  // add hit to the container - it add hit at the end - input in global coordinata
  //
  Double_t r = TMath::Sqrt(x*x+y*y);
  Double_t fi = TMath::ACos(x/r);
  if (y<0) fi*=-1.;
    AddHit(volumeID,trackID,r,z,fi,q);
}


void AliTPCTrackHitsV2::AddHit(Int_t volumeID, Int_t trackID, 
			     Double_t r, Double_t z, Double_t fi, Int_t q)
{
  //
  // Adding one hit
  //
  fSize++;
  Bool_t diff=kFALSE;
  if (!fTempInfo) { //initialisation of track  - initialisation of parameters
    fTempInfo = new AliTPCTempHitInfoV2;
    fTempInfo->fParam = new((*fArray)[0]) AliTrackHitsParamV2;
    fTempInfo->fParam->fVolumeID = volumeID;
    fTempInfo->fParam->fTrackID = trackID;
    AddVolume(volumeID);
    //
    fTempInfo->fParamIndex = 0;
    fTempInfo->NewParam(r,z,fi,q);
    return;
  }
    
  // if new volume or new trackID  
  if ( (volumeID!=fTempInfo->fParam->fVolumeID) || 
       (trackID!=fTempInfo->fParam->fTrackID)){
    if (volumeID!=fTempInfo->fParam->fVolumeID) AddVolume(volumeID);
    diff=kTRUE;
    FlushHitStack(kTRUE);        

    fTempInfo->fParamIndex++;   
    fTempInfo->fParam =  new((*fArray)[fTempInfo->fParamIndex]) AliTrackHitsParamV2;   
    fTempInfo->fParam->fVolumeID = volumeID;
    fTempInfo->fParam->fTrackID = trackID;   
    fTempInfo->NewParam(r,z,fi,q);
    return;
  }
     
  //calculate current fit precission to next point
  AliTrackHitsParamV2 &param = *(fTempInfo->fParam);
  Double_t dd=0;
  Double_t dl=0;
  Double_t ratio=0;
  Double_t dr,dz,dfi,ddz,ddfi;
  Double_t drhit,ddl;
  dr=dz=dfi=ddz=ddfi=0;
  drhit = r-fTempInfo->fOldR;
  { 
    //Double_t dfi2 = param.fAn+2*param.fAd*(r-param.fR); 
    Double_t dfi2 = param.fAn;
    dfi2*=dfi2*fTempInfo->fOldR*fTempInfo->fOldR;
    //Double_t ddz2 =  param.fTheta+2*param.fThetaD*(r-param.fR);
    Double_t ddz2 =  param.fTheta;
    ddz2*=ddz2;
    ratio = TMath::Sqrt(1.+ dfi2+ ddz2);  
  }
  //
  //  dl = fStep * Short_t(TMath::Nint(drhit*ratio/fStep));   // MI change - range check
  dl = drhit*ratio/fStep;
  if (TMath::Abs(dl)>32765) dl =0;
  dl = fStep * Short_t(TMath::Nint(dl));
  //
  ddl = dl - drhit*ratio; 
  fTempInfo->fOldR += dl/ratio; 

  if (fTempInfo->fStackIndex>2){     
    dr = r-param.fR;        
    dz =  z-param.fZ;  
    dfi = fi-param.fFi;
    ddz = dr*param.fTheta+dr*dr*param.fThetaD-dz;
    ddfi= dr*param.fAn+dr*dr*param.fAd-dfi;    
    dd  = TMath::Sqrt(ddz*ddz+r*r*ddfi*ddfi+ddl*ddl); 
    //
  }        
  //safety factor 1.25
  if ( ( (dd*1.25>fPrecision) ) ||  
       (fTempInfo->fStackIndex+4>fTempInfo->kStackSize) || 
       (TMath::Abs(dl/fStep)>fMaxDistance)  ) 
    diff=kTRUE;
  else{  // if precision OK
    fTempInfo->fStackIndex++;   
    fTempInfo->SetHit(r,z,fi,q);
    return;
  }  


  //if parameter changed 
  if (FlushHitStack(kFALSE)){   //if full buffer flushed
    fTempInfo->fParamIndex++;
    fTempInfo->fParam =  new((*fArray)[fTempInfo->fParamIndex]) AliTrackHitsParamV2;   
    fTempInfo->fParam->fVolumeID = volumeID;
    fTempInfo->fParam->fTrackID = trackID;   
    fTempInfo->NewParam(r,z,fi,q);
  }
  else{
    fTempInfo->fStackIndex++;
    fTempInfo->SetHit(r,z,fi,q);              
  }
}   

Bool_t AliTPCTrackHitsV2::FlushHitStack(Bool_t force)
{
  //
  // write fHitsPosAndQ information from the stack to te arrays
  //
  if (!fTempInfo) return kFALSE; 
 
  AliTrackHitsParamV2 & param = *(fTempInfo->fParam);
  //recalculate track parameter not fixing first point
  fTempInfo->UpdateParam(fStep/4.);
  //fTempInfo->Fit(fTempInfo->fParam);  //- fixing the first and the last point

  Double_t oldr = param.fR; 
  UInt_t i;
  Double_t dd;
  param.fNHits = fTempInfo->fStackIndex+1;
  if (param.fHitDistance) delete []param.fHitDistance;
  if (param.fCharge) delete []param.fCharge;
  param.fHitDistance = new Short_t[param.fNHits];
  param.fCharge = new Short_t[param.fNHits];

   
  for (i=0; i <= fTempInfo->fStackIndex; i++){
    Double_t * position = fTempInfo->GetPosition(i);
    Double_t   dr = position[0]-oldr;
    Double_t   ratio; 
    { 
      //Double_t dfi2 = param.fAn+2*param.fAd*(position[0]-param.fR);
      Double_t dfi2 = param.fAn;
      dfi2*=dfi2*oldr*oldr;
      //Double_t ddz2 =  param.fTheta+2*param.fThetaD*(position[0]-param.fR);
      Double_t ddz2 =  param.fTheta;
      ddz2*=ddz2;
      ratio = TMath::Sqrt(1.+ dfi2+ ddz2);  
    }

    //    Double_t dl = fStep*(Short_t)TMath::Nint(dr*ratio/fStep);   //MI change 
    Double_t dl = dr*ratio/fStep;
    if (TMath::Abs(dl)>32765) dl =0;
    dl = fStep * Short_t(TMath::Nint(dl));

    dr = dl/ratio; 
    oldr+=dr;
    //calculate precission
    AliTrackHitsParamV2 &param = *(fTempInfo->fParam);    
    //real deltas
    Double_t dr1=  position[0]-param.fR;
    Double_t dz =  position[1]-param.fZ;
    Double_t dfi = position[2]-param.fFi;
    //extrapolated deltas
    Double_t dr2 = oldr-param.fR; 
    Double_t ddr = dr2-dr1;
    Double_t ddz = dr2*param.fTheta+dr2*dr2*param.fThetaD-dz;
    Double_t ddfi= dr2*param.fAn+dr2*dr2*param.fAd-dfi;    
    dd = TMath::Sqrt(ddz*ddz+oldr*oldr*ddfi*ddfi+ddr*ddr); 


    if ( (dd>fPrecision) ){ 
      //if ( (dd<0) ){ 
      if (i==0){
	param.fAn = 0;
	param.fAd = 0;
	param.fTheta =0;
        param.fThetaD =0;
	Double_t ddz = dr2*param.fTheta+dr2*dr2*param.fThetaD-dz;
	Double_t ddfi= dr2*param.fAn+dr2*dr2*param.fAd-dfi;    
	dl = 0;
	dd = TMath::Sqrt(ddz*ddz+oldr*oldr*ddfi*ddfi+ddr*ddr); 
      }
      else
     	break;
    }

    param.fHitDistance[i]= Short_t(TMath::Nint(dl/fStep));
    param.fCharge[i]= Short_t(fTempInfo->fQStack[i]);
  }    
  
  if (i<=fTempInfo->fStackIndex){ //if previous iteration not succesfull 
    Short_t * charge = new Short_t[i];
    Short_t * hitDistance= new Short_t[i];
    memcpy(charge, param.fCharge,sizeof(Short_t)*i);
    memcpy(hitDistance, param.fHitDistance,sizeof(Short_t)*i);
    delete [] param.fCharge;
    delete [] param.fHitDistance;
    param.fNHits= i;
    param.fCharge = charge;
    param.fHitDistance = hitDistance;
    //
    Int_t volumeID = fTempInfo->fParam->fVolumeID;
    Int_t  trackID =fTempInfo->fParam->fTrackID;   
    fTempInfo->fParamIndex++;
    fTempInfo->fParam = new((*fArray)[fTempInfo->fParamIndex]) AliTrackHitsParamV2; 
    Double_t * p = fTempInfo->GetPosition(i);
    UInt_t index2 = fTempInfo->fStackIndex;
    fTempInfo->NewParam(p[0],p[1],p[2],fTempInfo->fQStack[i]);
    fTempInfo->fParam->fVolumeID= volumeID;
    fTempInfo->fParam->fTrackID= trackID;
    if (i+1<=index2) FlushHitStack2(i+1,index2);

    if (force) return      FlushHitStack(kTRUE);      
    return kFALSE;
  }  
  return kTRUE;
} 
 

void AliTPCTrackHitsV2::FlushHitStack2(Int_t index1, Int_t index2)
{
  //
  // second iteration flush stack
  // call only for hits where first iteration were not succesfully interpolated
  //
  Double_t * positionstack = new Double_t[3*(index2-index1+1)];
  UInt_t   * qstack        = new UInt_t[index2-index1+1];
  memcpy(positionstack, &fTempInfo->fPositionStack[3*index1],
	 (3*(index2-index1+1))*sizeof(Double_t));
  memcpy(qstack, &fTempInfo->fQStack[index1],(index2-index1+1)*sizeof(UInt_t));
  Double_t *p = positionstack;
  for (Int_t j=0; j<=index2-index1;j++){ 
    fTempInfo->fStackIndex++;
    fTempInfo->SetHit(p[3*j+0],p[3*j+1],p[3*j+2],qstack[j]);
  }  
  delete []positionstack;
  delete []qstack;
}


void AliTPCTrackHitsV2::AddVolume(Int_t volume)
{
  //
  //add volumes to tthe list of volumes
  //
  Int_t * volumes = new Int_t[fNVolumes+1];
  if (fVolumes) memcpy(volumes,fVolumes,(fNVolumes)*sizeof(Int_t));
  volumes[fNVolumes]=volume;
  fNVolumes++;
  if (fVolumes) delete []fVolumes;
  fVolumes = volumes;  
}


Bool_t AliTPCTrackHitsV2::First()
{
  //
  //set Current hit for the first hit
  //

  if (fArray->GetSize()<=0) {
    fCurrentHit->fStatus = kFALSE;
    return kFALSE;
  }

  AliTrackHitsParamV2 *param = (AliTrackHitsParamV2 *)fArray->At(0);
  if (!fHit) fHit = new AliTPChit;
  if (!(param) ) {
    fCurrentHit->fStatus = kFALSE;
    return kFALSE;
  }
  //
  fCurrentHit->fParamIndex = 0;
  fCurrentHit->fStackIndex = 0;
  //
  //
  ((AliTPChit*)fHit)->fSector = param->fVolumeID;
  ((AliTPChit*)fHit)->SetTrack(param->fTrackID);
  ((AliTPChit*)fHit)->SetX(param->fR*TMath::Cos(param->fFi));
  ((AliTPChit*)fHit)->SetY(param->fR*TMath::Sin(param->fFi));
  ((AliTPChit*)fHit)->SetZ(param->fZ); 
  ((AliTPChit*)fHit)->fQ = param->fCharge[0];     
  /*
    fCurrentHit->fHit.fSector = param->fVolumeID;
    fCurrentHit->fHit.SetTrack(param->fTrackID);
    fCurrentHit->fHit.SetX(param->fR*TMath::Cos(param->fFi));
    fCurrentHit->fHit.SetY(param->fR*TMath::Sin(param->fFi));
    fCurrentHit->fHit.SetZ(param->fZ); 
    fCurrentHit->fHit.fQ = param->fCharge[0];   
  */
  fCurrentHit->fR = param->fR;
  
  return fCurrentHit->fStatus = kTRUE;
}

Bool_t AliTPCTrackHitsV2::Next()
{
  //
  // Hit iterator  
  //
  if (!(fCurrentHit->fStatus)) 
    return kFALSE;

  fCurrentHit->fStackIndex++;

  AliTrackHitsParamV2 *param =  (AliTrackHitsParamV2 *)fArray->At(fCurrentHit->fParamIndex);
  if (fCurrentHit->fStackIndex>=param->fNHits){
    fCurrentHit->fParamIndex++;
    if (fCurrentHit->fParamIndex>=fArray->GetEntriesFast()){
      fCurrentHit->fStatus=kFALSE;
      return kFALSE;
    }
    param =  (AliTrackHitsParamV2 *)fArray->At(fCurrentHit->fParamIndex);
    fCurrentHit->fStackIndex=0; 
    fCurrentHit->fR = param->fR;
  }



  Double_t ratio;
  { 
    //    Double_t dfi2 = param->fAn+2*param->fAd*(fCurrentHit->fR-param->fR);
    Double_t dfi2 = param->fAn;
    dfi2*=dfi2*fCurrentHit->fR*fCurrentHit->fR;
    //    Double_t ddz2 = param->fTheta+2*param->fThetaD*(fCurrentHit->fR-param->fR);
    Double_t ddz2 =  param->fTheta;
    ddz2*=ddz2;
    ratio = TMath::Sqrt(1.+ dfi2+ ddz2);  
  }

  fCurrentHit->fR += fStep*param->fHitDistance[fCurrentHit->fStackIndex]/ratio;

  Double_t dR = fCurrentHit->fR - param->fR;
  Double_t fi = param->fFi + (param->fAn*dR+param->fAd*dR*dR);
  Double_t z  = param->fZ + (param->fTheta*dR+param->fThetaD*dR*dR);
  /*
  fCurrentHit->fHit.fQ = param->fCharge[fCurrentHit->fStackIndex];  
  fCurrentHit->fHit.SetX(fCurrentHit->fR*TMath::Cos(fi));
  fCurrentHit->fHit.SetY(fCurrentHit->fR*TMath::Sin(fi));
  fCurrentHit->fHit.SetZ(z);   
  fCurrentHit->fHit.fSector = param->fVolumeID;
  fCurrentHit->fHit.SetTrack(param->fTrackID);
  */
  ((AliTPChit*)fHit)->fQ = param->fCharge[fCurrentHit->fStackIndex];  
  ((AliTPChit*)fHit)->SetX(fCurrentHit->fR*TMath::Cos(fi));
  ((AliTPChit*)fHit)->SetY(fCurrentHit->fR*TMath::Sin(fi));
  ((AliTPChit*)fHit)->SetZ(z);   
  ((AliTPChit*)fHit)->fSector = param->fVolumeID;
  ((AliTPChit*)fHit)->SetTrack(param->fTrackID);

  return kTRUE;
}
  
AliHit * AliTPCTrackHitsV2::GetHit() const
{
  //
  // Return one hit
  //
   return (fCurrentHit->fStatus)? fHit:0;
  //return &fCurrentHit->fHit;

} 
 
AliTrackHitsParamV2 * AliTPCTrackHitsV2::GetParam()
{
  //
  // Return current parameters
  //
  return (fCurrentHit->fStatus)? 
    (AliTrackHitsParamV2 *)fArray->At(fCurrentHit->fParamIndex):0;
}

