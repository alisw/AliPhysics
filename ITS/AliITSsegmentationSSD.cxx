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

#include "AliITSsegmentationSSD.h"
#include "AliITSgeom.h"


ClassImp(AliITSsegmentationSSD)
AliITSsegmentationSSD::AliITSsegmentationSSD(){
  // default constructor
   fGeom=0;
   fCorr=0;
   SetDetSize();
   SetCellSize();
   SetNCells();
   SetAngles();
}
//------------------------------
AliITSsegmentationSSD::AliITSsegmentationSSD(AliITSgeom *geom){
  // constuctor
   fGeom=geom;
   fCorr=0;
   SetDetSize();
   SetCellSize();
   SetNCells();
   SetAngles();
   //Init(); 
}
//____________________________________________________________________________
AliITSsegmentationSSD& AliITSsegmentationSSD::operator=(AliITSsegmentationSSD &source){
// Operator =
     if(this==&source) return *this;
     this->fNstrips = source.fNstrips;
     this->fStereoP = source.fStereoP;
     this->fStereoN = source.fStereoN;
     this->fPitch   = source.fPitch;
     this->fDz      = source.fDz;
     this->fDx      = source.fDx;
     this->fDy      = source.fDy;
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

    AliITSgeomSSD *gssd = (AliITSgeomSSD *) (fGeom->GetShape(5,1,1));
    const Float_t kconv=10000.;
    fDx = 2.*kconv*gssd->GetDx();
    fDz = 2.*kconv*gssd->GetDz();
    fDy = 2.*kconv*gssd->GetDy();
    SetCellSize();
    SetNCells();
    SetAngles();

}
//-------------------------------------------------------
void AliITSsegmentationSSD::GetCellIxz(Float_t &x,Float_t &z,Int_t &iP,Int_t &iN)
{
  // returns P and N sided strip numbers for a given location.

    // expects x, z in microns

    Float_t tanP=TMath::Tan(fStereoP);
    Float_t tanN=TMath::Tan(fStereoN);

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
void AliITSsegmentationSSD::GetCellCxz(Int_t iP,Int_t iN,Float_t &x,Float_t &z)
{
    // actually this is the GetCrossing(Float_t &,Float_t &) 

    // returns x, z  in microns !

    Float_t flag=2*fDx;

    Float_t tanP=TMath::Tan(fStereoP);
    Float_t tanN=TMath::Tan(fStereoN);

    Float_t dx = 0.1;
    x = iP*fPitch;
    z = iN*fPitch; 

    if(tanP + tanN  == 0) {x=z=flag; return ;}

    z = (z - x + tanN * fDz) / (tanP + tanN);    
    x = x + tanP * z;                         

    x -= fDx/2;
    z -= fDz/2;

    if ( ( z < -(fDz/2+dx) ) || ( z > (fDz/2+dx) ) ) {x=z=flag; return ;}
    if ( ( x < -(fDx/2+dx) ) || ( x > (fDx/2+dx) ) ) {x=z=flag; return ;}

    return;   
}
