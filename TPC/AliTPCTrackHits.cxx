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
//  Class for storing simulated AliTPCHits  for given track                  //
//             -average compression comparing to classical ClonesArray is    //
//              around 5-7 (depending on the required hit precision)       //
//                                                                           //
//Begin_Html
/*
<img src="gif/AliTPCTrackHits.gif">
*/
//End_Html
//                                                                           //
//                                                                          //
///////////////////////////////////////////////////////////////////////////////

#include <TError.h>

#include "AliTPC.h"
#include "AliTPCTrackHits.h"
//           Interface classes
#include "AliTPCTrackHitsInterfaces.h"

ClassImp(AliTPCTrackHits) 
LClassImp(AliTrackHitsInfo) 
LClassImp(AliTrackHitsParam)  
LClassImp(AliHitInfo)

Int_t AliTrackHitsInfo::fgCounter1 =0;
Int_t AliTrackHitsInfo::fgCounter2 =0;
Int_t AliTrackHitsParam::fgCounter1 =0;
Int_t AliTrackHitsParam::fgCounter2 =0;
Int_t AliHitInfo::fgCounter1 =0;
Int_t AliHitInfo::fgCounter2 =0;
Int_t AliTPCTrackHits::fgCounter1 =0;
Int_t AliTPCTrackHits::fgCounter2 =0;
const Double_t AliTPCTrackHits::fgkPrecision=1e-6;  //precision 
const Double_t AliTPCTrackHits::fgkPrecision2=1e-20;  //precision
const Double_t AliTPCTrackHits::fgkTimePrecision=20.e-9;  //hit time precision 


class AliTPCCurrentHit {
private:
  AliTPChit fHit;     //   - hit in "standard" representation
  UInt_t   fInfoIndex;//   - current info pointer 
  UInt_t   fParamIndex;//  - current param pointer
  UInt_t   fStackIndex; // - current hit stack index
  Double_t fR;   //current Radius
  Bool_t  fStatus; //current status    
};   


class  AliTPCTempHitInfo {
private:
  enum    { kStackSize = 100};
  AliTPCTempHitInfo(); 
  AliTPCTempHitInfo(const AliTPCTempHitInfo &)
    {::Fatal("copy ctor","Not implemented\n");}
  AliTPCTempHitInfo & operator = (const AliTPCTempHitInfo &)
    {::Fatal("= operator","Not implemented\n");return *this;}

  void     NewParam(Double_t r, Double_t z, Double_t fi, Int_t q, Float_t time);
  void     SetHit(Double_t r, Double_t z, Double_t fi, Int_t q, Float_t time);
  Double_t * GetPosition(Int_t index){return &fPositionStack[index*3];}
  void    UpdateParam(Double_t maxdelta); //recal
  void   Fit2(Double_t fSumY, Double_t fSumYX, Double_t fSumYX2,
	    Double_t fSumX,  Double_t fSumX2, Double_t fSumX3, 
	    Double_t fSumX4, Int_t n,
	      Double_t &a, Double_t &b, Double_t &c);
  void  Fit(AliTrackHitsParam * param);
  Double_t fSumDr;    //fSumDr
  Double_t fSumDr2;   //fSumDr2
  Double_t fSumDr3;   //  fSumDr3
  Double_t fSumDr4;   //fSumDr4
  Double_t fSumDFi;  //fSumDFi
  Double_t fSumDFiDr; // fSumDFiDr 
  Double_t fSumDFiDr2;//fSumDFiDr2
  Double_t fSumDZ;     //fSumDZ
  Double_t fSumDZDr;  //fSumDZDr
  Double_t fSumDZDr2;  //fSumDZDr2
  Double_t fOldR;     //previos r
  Double_t fPositionStack[3*kStackSize];  //position stack 
  UInt_t   fQStack[kStackSize];           //Q stack
  Float_t  fTimeStack[kStackSize];        //time stack
  UInt_t fStackIndex;   //current stack index 
  UInt_t fInfoIndex;    //current track info index
  UInt_t fParamIndex;   //current track parameters index
  AliTrackHitsInfo  * fInfo; //current track info
  AliTrackHitsParam * fParam; //current track param
};


AliTPCTempHitInfo::AliTPCTempHitInfo()
{
  //
  //set to default value
  fSumDr=fSumDr2=fSumDr3=fSumDr4=
    fSumDFi=fSumDFiDr=fSumDFiDr2=
    fSumDZ=fSumDZDr=fSumDZDr2=0;  
  fStackIndex = 0;
  fInfoIndex  = 0;
  fParamIndex = 0;
}


void AliTPCTempHitInfo::NewParam(Double_t r, Double_t z, Double_t fi, Int_t q, Float_t time)
{
  //
  //reset stack and sum parameters
  //store line initial point
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
  SetHit(r,z,fi,q,time);
}

void AliTPCTempHitInfo::SetHit(Double_t r, Double_t z, Double_t fi, Int_t q, Float_t time)
{
  //
  //add hit to the stack
  //recalculate new estimete of line parameters
  Double_t *f = GetPosition(fStackIndex);  
  f[0] = r;
  f[1] = z;
  f[2] = fi;
  fQStack[fStackIndex]=q;
  fTimeStack[fStackIndex]=time;
  if (fStackIndex==0) return;
  Double_t dr  = (r-fParam->fR);
  if (TMath::Abs(dr)<AliTPCTrackHits::fgkPrecision) dr =AliTPCTrackHits::fgkPrecision;
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
  if (TMath::Abs(det)<AliTPCTrackHits::fgkPrecision2) return;
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


void   AliTPCTempHitInfo::UpdateParam(Double_t maxdelta)
{
  //recalc parameters not fixing origin point
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
void   AliTPCTempHitInfo::Fit2(Double_t fSumY, Double_t fSumYX, Double_t fSumYX2,
	    Double_t fSumX,  Double_t fSumX2, Double_t fSumX3, 
	    Double_t fSumX4, Int_t n,
	    Double_t &a, Double_t &b, Double_t &c)
{
  //fit of second order
  Double_t det = 
    n* (fSumX2*fSumX4-fSumX3*fSumX3) -
    fSumX*      (fSumX*fSumX4-fSumX3*fSumX2)+
    fSumX2*     (fSumX*fSumX3-fSumX2*fSumX2);
    
  if (TMath::Abs(det)> AliTPCTrackHits::fgkPrecision) {    
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

void   AliTPCTempHitInfo::Fit(AliTrackHitsParam * param)
{
  // fit fixing first and the last point 
  //result stored in new param
  Double_t dx2  = (GetPosition(fStackIndex))[0]-fParam->fR;
  Double_t det = fSumDr4+dx2*fSumDr2-2*dx2*fSumDr3;
  if ( (TMath::Abs(det)> AliTPCTrackHits::fgkPrecision) &&
       ((TMath::Abs(dx2)> AliTPCTrackHits::fgkPrecision))){
    Double_t dfi2 = (GetPosition(fStackIndex))[1]-fParam->fFi;
    param->fAd = (fSumDFiDr2+dfi2*fSumDr-dx2*fSumDFiDr-dfi2*fSumDr3/dx2)/det;
    param->fAn  = (dfi2-param->fAd*dx2*dx2)/dx2;
    
    Double_t dz2 = (GetPosition(fStackIndex))[1]-fParam->fZ;
    param->fTheta = (fSumDZDr2+dz2*fSumDr-dx2*fSumDZDr-dz2*fSumDr3/dx2)/det;
    param->fTheta  = (dz2-param->fAd*dx2*dx2)/dx2;
  }
  
}

//______________________________________________________________________
AliTrackHitsInfo::AliTrackHitsInfo() : 
  fTrackID(0),
  fVolumeID(0),
  fHitParamIndex(0)
{
  //
  // Default constructor
  //
  fgCounter1++;
  fgCounter2++;
}

//______________________________________________________________________
AliTrackHitsParam::AliTrackHitsParam() :
  fR(0),
  fZ(0),
  fFi(0),
  fAn(0),
  fAd(0),
  fTheta(0),
  fThetaD(0)
{
  //
  // Default constructor
  //
  fgCounter1++;
  fgCounter2++;
}


AliTPCTrackHits::AliTPCTrackHits()
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
  fTrackHitsInfo = new AliObjectArray("AliTrackHitsInfo"); 
  fTrackHitsParam = new AliObjectArray("AliTrackHitsParam");
  fHitsPosAndQ = new TArrayOfArrayVStack("AliHitInfo");
  fCurrentHit = new AliTPCCurrentHit;
  fgCounter1++;
  fgCounter2++;

} 

AliTPCTrackHits::AliTPCTrackHits(const AliTPCTrackHits& r) : TObject(r)
{
  //dummy
}
AliTPCTrackHits &AliTPCTrackHits::operator=(const AliTPCTrackHits& /* r */)
{
  //dummy
  return *this;
}

AliTPCTrackHits::~AliTPCTrackHits()
{
  //
  //default destructor
  //
  if (fTrackHitsInfo) delete fTrackHitsInfo;
  if (fTrackHitsParam) delete fTrackHitsParam;
  if (fHitsPosAndQ) delete fHitsPosAndQ;
  if (fCurrentHit) delete fCurrentHit;
  if (fTempInfo) delete fTempInfo;
  fgCounter1--;
}

void AliTPCTrackHits::Clear()
{
  //
  //clear object 
  fTrackHitsInfo->Clear();
  fTrackHitsParam->Clear();  
  //fTrackHitsInfo->Resize(0);
  //fTrackHitsParam->Resize(0);
  fHitsPosAndQ->Clear();

  if (fTempInfo){
    delete fTempInfo; 
    fTempInfo =0;
  } 
}


void AliTPCTrackHits::AddHitKartez(Int_t volumeID, Int_t trackID, Double_t x, 
	      Double_t y, Double_t z,Int_t q, Float_t time)
{
  //add hits (cartesian)
  Double_t r = TMath::Sqrt(x*x+y*y);
  Double_t fi = TMath::ACos(x/r);
  if (y<0) fi*=-1.;
    AddHit(volumeID,trackID,r,z,fi,q,time);
}

void AliTPCTrackHits::AddHit(Int_t volumeID, Int_t trackID, 
			     Double_t r, Double_t z, Double_t fi, Int_t q, Float_t time)
{
  //
  Bool_t diff=kFALSE;
  if (!fTempInfo) { //initialsation of track
    fTempInfo = new AliTPCTempHitInfo;
    //
    if (fTrackHitsInfo->GetCapacity()<10) fTrackHitsInfo->Reserve(10);
    fTrackHitsInfo->Resize(1);
    fTempInfo->fInfoIndex =0;
    if (fTrackHitsParam->GetCapacity()<100) fTrackHitsParam->Reserve(100);    
    fTrackHitsParam->Resize(1);
    //
    fTempInfo->fInfo = 
      (AliTrackHitsInfo*) (fTrackHitsInfo->At(0));
    fTempInfo->fInfo->fVolumeID = volumeID;
    fTempInfo->fInfo->fTrackID = trackID;
    fTempInfo->fInfo->fHitParamIndex =0;     
    fTempInfo->fInfoIndex = 0;
    //
    fTempInfo->fParam = 
      (AliTrackHitsParam*) (fTrackHitsParam->At(0));
    fTempInfo->fParamIndex = 0;
    fTempInfo->NewParam(r,z,fi,q,time);
    return;
  }
  
  Int_t size = fHitsPosAndQ->GetSize();
  if (size>(Int_t)fTempInfo->fParamIndex) {
    fTempInfo->fParamIndex++;
    if (fTempInfo->fParamIndex+1>fTrackHitsParam->GetSize()) 
      fTrackHitsParam->Resize(fTempInfo->fParamIndex+1);    
    fTempInfo->fParam = 
      (AliTrackHitsParam*) (fTrackHitsParam->At(fTempInfo->fParamIndex));  
    fTempInfo->NewParam(r,z,fi,q,time);
    return;
  }
  

  // if new volume or new trackID  
  if ( (volumeID!=fTempInfo->fInfo->fVolumeID) || 
       (trackID!=fTempInfo->fInfo->fTrackID)){
    diff=kTRUE;

    FlushHitStack(kTRUE);
      
    fTempInfo->fInfoIndex++;
    if (fTempInfo->fInfoIndex+1>fTrackHitsInfo->GetSize()) 
      fTrackHitsInfo->Resize(fTempInfo->fInfoIndex+1);
    fTempInfo->fInfo = 
      (AliTrackHitsInfo*) (fTrackHitsInfo->At(fTempInfo->fInfoIndex));      
    fTempInfo->fInfo->fVolumeID = volumeID;
    fTempInfo->fInfo->fTrackID = trackID;
    fTempInfo->fInfo->fHitParamIndex =fTempInfo->fParamIndex+1;  
    // FlushHitStack(kTRUE);

    fTempInfo->fParamIndex++;
    if (fTempInfo->fParamIndex+1>fTrackHitsParam->GetSize()) 
      fTrackHitsParam->Resize(fTempInfo->fParamIndex+1);    
    fTempInfo->fParam = 
      (AliTrackHitsParam*) (fTrackHitsParam->At(fTempInfo->fParamIndex));  
    fTempInfo->NewParam(r,z,fi,q,time);
    return;
  }
     
  //calculate current fit precission to next point
  AliTrackHitsParam &param = *(fTempInfo->fParam);
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
 

  dl = (TMath::Abs(drhit*ratio/fStep)<32000) ?  fStep * Short_t(TMath::Nint(drhit*ratio/fStep)):0;
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
  else{
    fTempInfo->fStackIndex++;   
    fTempInfo->SetHit(r,z,fi,q,time);
    return;
  }  
  //if parameter changed 
  if (FlushHitStack(kFALSE)){   //if full buffer flushed
    fTempInfo->fParamIndex++;
    if (fTempInfo->fParamIndex+1>fTrackHitsParam->GetSize()) 
      fTrackHitsParam->Resize(fTempInfo->fParamIndex+1);    
    fTempInfo->fParam = 
      (AliTrackHitsParam*) (fTrackHitsParam->At(fTempInfo->fParamIndex));  
    fTempInfo->NewParam(r,z,fi,q,time);
  }
  else{
    fTempInfo->fStackIndex++;
    fTempInfo->SetHit(r,z,fi,q,time);              
  }
}   

Bool_t AliTPCTrackHits::FlushHitStack(Bool_t force)
{
  //
  //write fHitsPosAndQ information from the stack to te arrays
  if (!fTempInfo) return kFALSE; 
  Int_t size = fHitsPosAndQ->GetSize();

  if ( (size>0)&&(size!=(Int_t)fTempInfo->fParamIndex)) return kFALSE;  

  if (fHitsPosAndQ->Push(fTempInfo->fStackIndex+1)!=fTempInfo->fParamIndex){
    cout<<"internal error - contact MI\n";
    return kFALSE;
  }
  AliHitInfo * info;
 
  AliTrackHitsParam & param = *(fTempInfo->fParam);
  //recalculate track parameter not fixing first point
  fTempInfo->UpdateParam(fStep/4.);
  //fTempInfo->Fit(fTempInfo->fParam);  //- fixing the first and the last point

  Double_t oldr = param.fR; 
  //cout<<"C3"<<fTempInfo->fStackIndex<<"\n"<<flush;
  UInt_t i;
  Double_t dd;
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

    Double_t dl = (TMath::Abs(dr*ratio/fStep)<32000) ? fStep*(Short_t)TMath::Nint(dr*ratio/fStep):0;  
    dr = dl/ratio; 
    oldr+=dr;
    //calculate precission
    AliTrackHitsParam &param = *(fTempInfo->fParam);    
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

    info = (AliHitInfo*)(fHitsPosAndQ->At(fTempInfo->fParamIndex,i));
    info->fHitDistance = (TMath::Abs(dl/fStep)<32000) ?Short_t(TMath::Nint(dl/fStep)):0;
    info->fCharge = Short_t(fTempInfo->fQStack[i]);
    info->fTime = TMath::Nint(fTempInfo->fTimeStack[i]/AliTPCTrackHits::fgkTimePrecision);
    /*
    cout<<"C2";
    cout<<" "<<fTempInfo->fStackIndex<<" \t";
    cout<<" "<<i<<" \t";
    cout<<position[0]<<"\t";
    cout<<position[1]<<"\t"; 
    cout<<position[2]<<"\t";
    cout<<param.fAn<<"\t";
    cout<<param.fTheta<<"\t";
    cout<<dr1<<"\t"<<ddr<<"\t"<<ddz<<"\t"<<ddfi<<"\t"<<dd<<"\n"<<flush;    
    */
  }    
  
  if (i<=fTempInfo->fStackIndex){ //if previous iteration not succesfull 
    fHitsPosAndQ->Resize(fTempInfo->fParamIndex,i);
    //
    fTempInfo->fParamIndex++;
    if (fTempInfo->fParamIndex+1>fTrackHitsParam->GetSize()) 
      fTrackHitsParam->Resize(fTempInfo->fParamIndex+1);    
    fTempInfo->fParam = 
	(AliTrackHitsParam*) (fTrackHitsParam->At(fTempInfo->fParamIndex));    
    Double_t * p = fTempInfo->GetPosition(i);
    UInt_t index2 = fTempInfo->fStackIndex;
    fTempInfo->NewParam(p[0],p[1],p[2],fTempInfo->fQStack[i],fTempInfo->fTimeStack[i]);        
    if (i+1<=index2) FlushHitStack2(i+1,index2);

    if (force) return      FlushHitStack(kTRUE);      
    return kFALSE;
  }  
  return kTRUE;
} 
 

void AliTPCTrackHits::FlushHitStack2(Int_t index1, Int_t index2)
{
  //
  // second iteration flush stack
  // call only for hits where first iteration were not succesfully interpolated  
  Double_t * positionstack = new Double_t[3*(index2-index1+1)];
  UInt_t   * qstack        = new UInt_t[index2-index1+1];
  Float_t  * timestack        = new Float_t[index2-index1+1];
  memcpy(positionstack, &fTempInfo->fPositionStack[3*index1],
	 (3*(index2-index1+1))*sizeof(Double_t));
  memcpy(qstack, &fTempInfo->fQStack[index1],(index2-index1+1)*sizeof(UInt_t));
  memcpy(timestack, &fTempInfo->fTimeStack[index1],(index2-index1+1)*sizeof(Float_t));
  Double_t *p = positionstack;
  for (Int_t j=0; j<=index2-index1;j++){ 
    fTempInfo->fStackIndex++;
    fTempInfo->SetHit(p[3*j+0],p[3*j+1],p[3*j+2],qstack[j],timestack[j]);
  }  
  delete []positionstack;
  delete []qstack;
  delete []timestack;
}


   




  

Bool_t AliTPCTrackHits::First()
{
  //
  //set Current hit for the first hit
  //
  AliTrackHitsInfo *info = (AliTrackHitsInfo *)fTrackHitsInfo->At(0);
  AliTrackHitsParam *param = (AliTrackHitsParam *)fTrackHitsParam->At(0);
  AliHitInfo * hinfo = (AliHitInfo *)fHitsPosAndQ->At(0,0);

  if (!(info) || !(param) || !(hinfo) ) {
    fCurrentHit->fStatus = kFALSE;
    return kFALSE;
  }

  fCurrentHit->fInfoIndex  = 0;
  fCurrentHit->fParamIndex = 0;
  fCurrentHit->fStackIndex = 0;

  fCurrentHit->fHit.fSector = info->fVolumeID;
  fCurrentHit->fHit.SetTrack(info->fTrackID);
  fCurrentHit->fHit.SetX(param->fR*TMath::Cos(param->fFi));
  fCurrentHit->fHit.SetY(param->fR*TMath::Sin(param->fFi));
  fCurrentHit->fHit.SetZ(param->fZ); 
  fCurrentHit->fHit.fQ = (Float_t)(hinfo->fCharge*AliTPCTrackHits::fgkTimePrecision);
  fCurrentHit->fHit.fTime = hinfo->fTime;
   
  fCurrentHit->fR = param->fR;
  
  return fCurrentHit->fStatus = kTRUE;
}


/*
Bool_t AliTPCTrackHits::Next()
{
  //
  //  
  if (!(fCurrentHit->fStatus)) 
    return kFALSE;

  fCurrentHit->fStackIndex++;
  AliHitInfo * hinfo = (AliHitInfo *)fHitsPosAndQ->At(fCurrentHit->fParamIndex,
						      fCurrentHit->fStackIndex);
  AliTrackHitsInfo *info = (AliTrackHitsInfo *)fTrackHitsInfo->At(fCurrentHit->fInfoIndex);
  AliTrackHitsParam *param =  (AliTrackHitsParam *)fTrackHitsParam->At(fCurrentHit->fParamIndex);

  if (!hinfo) {
    hinfo = (AliHitInfo *)fHitsPosAndQ->At(fCurrentHit->fParamIndex+1, 0);
    if (!hinfo) 
      return fCurrentHit->fStatus = kFALSE;
    if (hinfo){ 
      fCurrentHit->fParamIndex++;
      fCurrentHit->fStackIndex = 0;
      param = (AliTrackHitsParam *)fTrackHitsParam->At(fCurrentHit->fParamIndex);
      if (!param) 
	return fCurrentHit->fStatus = kFALSE;     
      fCurrentHit->fR = param->fR;

      if ((fCurrentHit->fInfoIndex+1<fTrackHitsInfo->GetSize())
	&&((info+1)->fHitParamIndex<=fCurrentHit->fParamIndex)){
	fCurrentHit->fInfoIndex++;
	info = (AliTrackHitsInfo *)fTrackHitsInfo->At(fCurrentHit->fInfoIndex);
	if (!info) 
	  return fCurrentHit->fStatus = kFALSE;
	fCurrentHit->fHit.fSector = info->fVolumeID;
	fCurrentHit->fHit.SetTrack(info->fTrackID);
      }
    }  
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

  fCurrentHit->fHit.fQ = hinfo->fCharge;
  fCurrentHit->fHit.fTime = (Float_t)(hinfo->fTime*AliTPCTrackHits::fgkTimePrecision);
  fCurrentHit->fR += fStep*hinfo->fHitDistance/ratio;
  Double_t dR = fCurrentHit->fR - param->fR;
  //Double_t dR =0;
  Double_t fi = param->fFi + (param->fAn*dR+param->fAd*dR*dR);
  Double_t z  = param->fZ + (param->fTheta*dR+param->fThetaD*dR*dR);
  
  fCurrentHit->fHit.SetX(fCurrentHit->fR*TMath::Cos(fi));
  fCurrentHit->fHit.SetY(fCurrentHit->fR*TMath::Sin(fi));
  fCurrentHit->fHit.SetZ(z);  
  return kTRUE;
}

*/  
AliTPChit * AliTPCTrackHits::GetHit() const
{
  //
   return (fCurrentHit->fStatus)? &fCurrentHit->fHit:0;
  //return &fCurrentHit->fHit;

}  



Bool_t AliTPCTrackHits::Next(Int_t id)
{
  //
  //  
  if (!(fCurrentHit->fStatus)) 
    return kFALSE;

  //  fCurrentHit->fStackIndex++;
  AliHitInfo * hinfo = (AliHitInfo *)fHitsPosAndQ->At(fCurrentHit->fParamIndex,
						      fCurrentHit->fStackIndex);
  AliTrackHitsInfo *info = (AliTrackHitsInfo *)fTrackHitsInfo->At(fCurrentHit->fInfoIndex);
  if (!info) {
    fCurrentHit->fStatus = kFALSE;
    return kFALSE; 
  }
  AliTrackHitsParam *param =  (AliTrackHitsParam *)fTrackHitsParam->At(fCurrentHit->fParamIndex);

  if ( (id>=0) && (info!=0) && (info->fVolumeID!=id)){
    fCurrentHit->fInfoIndex++;
    info = (AliTrackHitsInfo *)fTrackHitsInfo->At(fCurrentHit->fInfoIndex);
    if (!info) {
      fCurrentHit->fStatus = kFALSE;
      return kFALSE;
    }
    fCurrentHit->fParamIndex = info->fHitParamIndex;
    param =  (AliTrackHitsParam *)fTrackHitsParam->At(fCurrentHit->fParamIndex);
    fCurrentHit->fStackIndex =0;    
    fCurrentHit->fR = param->fR;
    return Next(id);
  }
  if (!info) {
    fCurrentHit->fStatus = kFALSE;
    return kFALSE; 
  }
  if (!hinfo) {
    hinfo = (AliHitInfo *)fHitsPosAndQ->At(fCurrentHit->fParamIndex+1, 0);
    if (!hinfo){ 
      fCurrentHit->fStatus = kFALSE;
      return kFALSE;
    }
    if (hinfo){ 
      fCurrentHit->fParamIndex++;
      fCurrentHit->fStackIndex = 0;
      param = (AliTrackHitsParam *)fTrackHitsParam->At(fCurrentHit->fParamIndex);
      if (!param){ 
	fCurrentHit->fStatus = kFALSE;
	return kFALSE;     
      }
      fCurrentHit->fR = param->fR;

      if ((fCurrentHit->fInfoIndex+1<fTrackHitsInfo->GetSize())
	&&((info+1)->fHitParamIndex<=fCurrentHit->fParamIndex)){
	fCurrentHit->fInfoIndex++;
	info = (AliTrackHitsInfo *)fTrackHitsInfo->At(fCurrentHit->fInfoIndex);
	if (!info){ 
	  fCurrentHit->fStatus = kFALSE;
	  return kFALSE;
	}
	if ( (id>=0) && (info!=0) && (info->fVolumeID!=id)){
	  return Next(id);
	}
	fCurrentHit->fHit.fSector = info->fVolumeID;
	fCurrentHit->fHit.SetTrack(info->fTrackID);
      }
    }  
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

  fCurrentHit->fHit.fQ = hinfo->fCharge;
  fCurrentHit->fHit.fTime = (Float_t)(hinfo->fTime*AliTPCTrackHits::fgkTimePrecision);
  fCurrentHit->fR += fStep*hinfo->fHitDistance/ratio;
  Double_t dR = fCurrentHit->fR - param->fR;
  //Double_t dR =0;
  Double_t fi = param->fFi + (param->fAn*dR+param->fAd*dR*dR);
  Double_t z  = param->fZ + (param->fTheta*dR+param->fThetaD*dR*dR);
  
  fCurrentHit->fHit.SetX(fCurrentHit->fR*TMath::Cos(fi));
  fCurrentHit->fHit.SetY(fCurrentHit->fR*TMath::Sin(fi));
  fCurrentHit->fHit.SetZ(z);  
  fCurrentHit->fHit.fSector = info->fVolumeID;
  fCurrentHit->fHit.SetTrack(info->fTrackID);
  //
  fCurrentHit->fStatus = kTRUE;
  fCurrentHit->fStackIndex++;
  return kTRUE;
}
  

AliTrackHitsParam * AliTPCTrackHits::GetParam()
{
  //
  return (fCurrentHit->fStatus)? (AliTrackHitsParam *)fTrackHitsParam->At(fCurrentHit->fParamIndex) :0;
} 

AliHitInfo * AliTPCTrackHits::GetHitInfo()
{
  //
  return (fCurrentHit->fStatus)? 
    (AliHitInfo *)fHitsPosAndQ->At(fCurrentHit->fParamIndex,fCurrentHit->fStackIndex) :0;
} 



