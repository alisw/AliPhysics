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
//------------------------------------------------------------------
// this class is used for a simple
// track fitting using the polynom
// bending in both directions are allowed
//-----------------------------------------------------------------


#include "AliTPCpolyTrack.h"
#include "TMath.h"

ClassImp(AliTPCpolyTrack)


AliTPCpolyTrack::AliTPCpolyTrack()
  :TObject(),
   fA(0), 
   fB(0), 
   fC(0), 
   fD(0), 
   fE(0), 
   fF(0), 
   fMaxX(0), 
   fMinX(0), 
   fSumW(0),    
   fSumX(0),    
   fSumX2(0),   
   fSumX3(0),    
   fSumX4(0),   
   fSumY(0),    
   fSumYX(0),   
   fSumYX2(0),  
   fSumZ(0),     
   fSumZX(0),   
   fSumZX2(0), 
   fNPoints(0)
{
  Reset();
}

void   AliTPCpolyTrack::Reset()
{
  //
  // reset track
  fSumX = fSumX2= fSumX3=fSumX4 = fSumY=fSumYX=fSumYX2=fSumZ=fSumZX=fSumZX2=fSumW =0;
  fNPoints = 0;
}

void AliTPCpolyTrack::AddPoint(Double_t x, Double_t y, Double_t z,Double_t sy, Double_t sz)
{
  //
  //
  if (fNPoints==0){
    fMaxX = x;
    fMinX = x;
  }else{
    if (x>fMaxX) fMaxX=x;
    if (x<fMinX) fMinX=x;
  }

  Double_t x2 = x*x; 
  Double_t w = 2./(sy+sz);
  fSumW += w;
  //
  fSumX       += x*w;
  fSumX2      += x2*w;
  fSumX3      += x2*x*w;
  fSumX4      += x2*x2*w;
  //
  fSumY       +=y*w;
  fSumYX      +=y*x*w;
  fSumYX2     +=y*x2*w;
  //
  fSumZ       +=z*w;
  fSumZX      +=z*x*w;
  fSumZX2     +=z*x2*w;
  //
  fX[fNPoints] = x;
  fY[fNPoints] = y;
  fZ[fNPoints] = z;  
  fSY[fNPoints] = sy;
  fSZ[fNPoints] = sz;  

  fNPoints++;
  
}

void  AliTPCpolyTrack::UpdateParameters()
{
  //
  //
  //Update fit parameters
  if (fNPoints>4){
    Fit2(fSumY,fSumYX,fSumYX2,fSumX,fSumX2,fSumX3,fSumX4,fSumW,fA,fB,fC);
    //    Fit2(fSumZ,fSumZX,fSumZX2,fSumX,fSumX2,fSumX3,fSumX4,fNPoints,fD,fE,fF);
    Fit1(fSumZ,fSumZX,fSumX,fSumX2,fSumW,fD,fE,fF);
  }
  else
    {
      Fit1(fSumY,fSumYX,fSumX,fSumX2,fSumW,fA,fB,fC);
      Fit1(fSumZ,fSumZX,fSumX,fSumX2,fSumW,fD,fE,fF);
    }
}

void  AliTPCpolyTrack::UpdateParameters(Int_t ny, Int_t nz)
{
  //
  //
  //Update fit parameters
  if (ny>1){
    Fit2(fSumY,fSumYX,fSumYX2,fSumX,fSumX2,fSumX3,fSumX4,fSumW,fA,fB,fC);
  }
  else{ 
    Fit1(fSumY,fSumYX,fSumX,fSumX2,fSumW,fA,fB,fC);
  }
  if (nz>1){
    Fit2(fSumZ,fSumZX,fSumZX2,fSumX,fSumX2,fSumX3,fSumX4,fNPoints,fD,fE,fF);
  }
  else{
    Fit1(fSumZ,fSumZX,fSumX,fSumX2,fSumW,fD,fE,fF);
  }
}



void  AliTPCpolyTrack::Fit2(Double_t fSumY, Double_t fSumYX, Double_t fSumYX2,
	    Double_t fSumX,  Double_t fSumX2, Double_t fSumX3, 
	    Double_t fSumX4, Double_t fSumW,
	    Double_t &a, Double_t &b, Double_t &c)
{
  //fit of second order
  Double_t det = 
    fSumW* (fSumX2*fSumX4-fSumX3*fSumX3) -
    fSumX*      (fSumX*fSumX4-fSumX3*fSumX2)+
    fSumX2*     (fSumX*fSumX3-fSumX2*fSumX2);
    
  if (TMath::Abs(det)> 0.000000000000001) {    
    a = 
      (fSumY * (fSumX2*fSumX4-fSumX3*fSumX3)-
       fSumX *(fSumYX*fSumX4-fSumYX2*fSumX3)+
       fSumX2*(fSumYX*fSumX3-fSumYX2*fSumX2))/det; 
    b=
      (fSumW*(fSumYX*fSumX4-fSumX3*fSumYX2)-
      fSumY*(fSumX*fSumX4-fSumX3*fSumX2)+
      fSumX2*(fSumX*fSumYX2-fSumYX*fSumX2))/det;
    c=
      (fSumW*(fSumX2*fSumYX2-fSumYX*fSumX3)-
       fSumX*(fSumX*fSumYX2-fSumYX*fSumX2)+
       fSumY*(fSumX*fSumX3-fSumX2*fSumX2))/det;  
  }
}

void  AliTPCpolyTrack::Fit1(Double_t fSumY, Double_t fSumYX, 
	      Double_t fSumX,  Double_t fSumX2, 
	      Double_t fSumW, Double_t &a, Double_t &b, Double_t &c)
{
  //
  //
  //
  Double_t det = fSumW*fSumX2-fSumX*fSumX;
  if (TMath::Abs(det)> 0.000000000000001) { 
    b = (fSumW*fSumYX-fSumX*fSumY)/det;
    a = (fSumX2*fSumY-fSumX*fSumYX)/det;
    c = 0;
  }else{
    a =fSumYX/fSumX;
    b =0;
    c =0;
  }

}

void AliTPCpolyTrack::Refit(AliTPCpolyTrack &track, Double_t deltay, Double_t deltaz)
{
  //
  // refit with cut on distortion
  //
  track.Reset();
  //first refit to temporary
  AliTPCpolyTrack track0;
  track0.Reset();
  for (Int_t i=0;i<fNPoints;i++){
    Double_t y,z;
    GetFitPoint(fX[i],y,z);
    if ( (TMath::Abs(y-fY[i])<deltay)&&(TMath::Abs(z-fZ[i])<deltaz)){
      track0.AddPoint(fX[i],y,z);
    }
  }
  if (track0.GetN()>2) 
    track0.UpdateParameters();
  else 
    return;
  //
  for (Int_t i=0;i<fNPoints;i++){
    Double_t y,z;
    track0.GetFitPoint(fX[i],y,z);
    if ( (TMath::Abs(y-fY[i])<deltay)&&(TMath::Abs(z-fZ[i])<deltaz)){
      track.AddPoint(fX[i],y,z);
    }
  }
  if (track.GetN()>2) 
    track.UpdateParameters();

}

void AliTPCpolyTrack::Refit(AliTPCpolyTrack &track, Double_t deltay, Double_t deltaz, Int_t nfirst, Int_t ny, Int_t nz)
{
  //
  // refit with cut on distortion
  //
  track.Reset();
  //first refit to temporary
  AliTPCpolyTrack track0;
  track0.Reset();
  for (Int_t i=0;i<fNPoints;i++){
    Double_t y,z;
    GetFitPoint(fX[i],y,z);
    if ( (TMath::Abs(y-fY[i])<deltay)&&(TMath::Abs(z-fZ[i])<deltaz)){
      track0.AddPoint(fX[i],y,z);
    }    
  }
  if (track0.GetN()>2){ 
    if (track0.GetN()>nfirst)
      track0.UpdateParameters(ny,nz);
    else 
      track0.UpdateParameters(1,1);
  }
  else 
    return;
  //
  for (Int_t i=0;i<fNPoints;i++){
    Double_t y,z;
    track0.GetFitPoint(fX[i],y,z);
    if ( (TMath::Abs(y-fY[i])<deltay)&&(TMath::Abs(z-fZ[i])<deltaz)){
      track.AddPoint(fX[i],y,z);
    }
  }
  if (track.GetN()>2) 
    track.UpdateParameters(ny,nz);

}

