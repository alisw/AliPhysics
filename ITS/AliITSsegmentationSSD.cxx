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

#include <TMath.h>
#include <TF1.h>
#include <iostream.h>
#include "AliITSsegmentationSSD.h"
#include "AliITSgeom.h"


ClassImp(AliITSsegmentationSSD)
AliITSsegmentationSSD::AliITSsegmentationSSD(){
  // default constructor
   fGeom=0;
   fCorr=0;
   SetDetSize();
   SetPadSize();
   SetNPads();
   SetAngles();
   fLayer =0;
}
//------------------------------
AliITSsegmentationSSD::AliITSsegmentationSSD(AliITSgeom *geom){
  // constuctor
   fGeom=geom;
   fCorr=0;
   SetDetSize();
   SetPadSize();
   SetNPads();
   SetAngles();
   fLayer =0;

}
//____________________________________________________________________________
AliITSsegmentationSSD& AliITSsegmentationSSD::operator=(AliITSsegmentationSSD &source){
// Operator =
     if(this==&source) return *this;
     this->fNstrips = source.fNstrips;
     this->fStereoP = source.fStereoP;
     this->fStereoN = source.fStereoN;
     this->fStereoPl5 = source.fStereoPl5;
     this->fStereoNl5 = source.fStereoNl5;
     this->fStereoPl6 = source.fStereoPl6;
     this->fStereoNl6 = source.fStereoNl6;
     this->fLayer   = source.fLayer;
     this->fPitch   = source.fPitch;
     this->fDz      = source.fDz;
     this->fDx      = source.fDx;
     this->fDy      = source.fDy;
     this->fLayer   = source.fLayer;
     this->fGeom    = source.fGeom; // copy only the pointer
     this->fCorr    = new TF1(*(source.fCorr)); // make a proper copy
     return *this;
     
}
//____________________________________________________________________________
AliITSsegmentationSSD::AliITSsegmentationSSD(AliITSsegmentationSSD &source){
  // copy constructor
  *this = source;
}
//------------------------------
void AliITSsegmentationSSD::Init(){
  // standard initalizer

    SetPadSize();
    SetNPads();
    SetAngles();

}
//-------------------------------------------------------
void AliITSsegmentationSSD::Angles(Float_t &aP,Float_t &aN)
     {
	 if (fLayer == 5)
	 {
         aP = fStereoPl5;
         aN = fStereoNl5;
	 }
   
	 if (fLayer == 6)
	 {
         aP = fStereoPl6;
         aN = fStereoNl6;
	 }
     }

     void AliITSsegmentationSSD::SetLayer(Int_t l)
     {
     if (l==5) fLayer =5;
     if (l==6) fLayer =6;
     }

//-------------------------------------------------------
void AliITSsegmentationSSD::GetPadIxz(Float_t x,Float_t z,Int_t &iP,Int_t &iN)
{
  // returns P and N sided strip numbers for a given location.

    // expects x, z in microns

    Float_t StereoP, StereoN;
    Angles(StereoP,StereoN);
    Float_t tanP=TMath::Tan(StereoP);
    Float_t tanN=TMath::Tan(StereoN);
    Float_t x1=x,z1=z;
    x1 += fDx/2;
    z1 += fDz/2;
    Float_t  ldX = x1 - z1*tanP;          // distance from left-down edge 
    iP = (Int_t)(ldX/fPitch);
    iP = (iP<0)? -1: iP;      
    iP = (iP>fNstrips)? -1: iP;

    ldX = x1 - tanN*(fDz - z1);
    iN = (Int_t)(ldX/fPitch);
    iN = (iN<0)? -1: iN;
    iN = (iN>fNstrips)? -1: iN;

}
//-------------------------------------------------------
void AliITSsegmentationSSD::GetPadCxz(Int_t iP,Int_t iN,Float_t &x,Float_t &z)
{
    // actually this is the GetCrossing(Float_t &,Float_t &) 

    // returns local x, z  in microns !

  Float_t Dx = fDx; // detector size in x direction, microns
  Float_t Dz = fDz; // detector size in z direction, microns
  Float_t xP; // x coordinate in the P side from the first P strip
  Float_t xN; // x coordinate in the N side from the first N strip
  Float_t StereoP, StereoN;
  Angles(StereoP,StereoN);
  Float_t kP=TMath::Tan(StereoP);
  Float_t kN=TMath::Tan(StereoN);

    xP=iP*fPitch;
    xN=iN*fPitch; 
    x = xP + kP*(Dz*kN-xP+xN)/(kP+kN);
    z = (Dz*kN-xP+xN)/(kP+kN); 
    x -= Dx/2;
    z -= Dz/2;
    //if(TMath::Abs(z) > Dz/2) cout<<"Warning, wrong z local ="<<z<<endl; 
    // Check that zL is inside the detector for the 
    // correspondent xP and xN coordinates

    return;   
}


