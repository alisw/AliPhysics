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

#include "AliITSsegmentationSDD.h"
#include "AliITS.h"
#include "AliITSgeom.h"
#include "AliRun.h"
#include "AliITSresponse.h"


ClassImp(AliITSsegmentationSDD)
//------------------------------
AliITSsegmentationSDD::AliITSsegmentationSDD(AliITSgeom* geom, AliITSresponse *resp){
  // constructor
   fGeom=geom;
   fResponse=resp;
   fCorr=0;
   SetDetSize();
   SetPadSize();
   SetNPads();

}
//_____________________________________________________________________________
AliITSsegmentationSDD::AliITSsegmentationSDD(){
  // standard constructor
   fGeom=0;
   fResponse=0;  
   fCorr=0;
   SetDetSize();
   SetPadSize();
   SetNPads();

}
//_____________________________________________________________________________
AliITSsegmentationSDD& AliITSsegmentationSDD::operator=(AliITSsegmentationSDD &source){
  // Operator =
  if(this==&source) return *this;
  this->fNsamples = source.fNsamples;
  this->fNanodes  = source.fNanodes;
  this->fPitch    = source.fPitch;
  this->fTimeStep = source.fTimeStep;
  this->fDx       = source.fDx;
  this->fDz       = source.fDz;
  this->fDy       = source.fDy;
  this->fCorr     = new TF1(*(source.fCorr));
  this->fGeom     = source.fGeom; // Just copy the pointer
  this->fResponse = source.fResponse; //Just copy the pointer
  return *this;
}
//___________________________________________________________________________
AliITSsegmentationSDD::AliITSsegmentationSDD(AliITSsegmentationSDD &source){
  // Copy constructor
   *this = source;
}
//------------------------------
void AliITSsegmentationSDD::Init(){
  // Standard initilisation routine

   if(!fGeom) {
     return;
     //fGeom = ((AliITS*)gAlice->GetModule("ITS"))->GetITSgeom();
   }
   AliITSgeomSDD *gsdd = (AliITSgeomSDD *) (fGeom->GetShape(3,1,1));

   const Float_t kconv=10000.;
   fDz = 2.*kconv*gsdd->GetDz();
   fDx = kconv*gsdd->GetDx();
   fDy = 2.*kconv*gsdd->GetDy();
}

//------------------------------
void AliITSsegmentationSDD::
Neighbours(Int_t iX, Int_t iZ, Int_t* Nlist, Int_t Xlist[8], Int_t Zlist[8]){
  // returns neighbours for use in Cluster Finder routines and the like

    if(iX >= fNanodes) printf("iX > fNanodes %d %d\n",iX,fNanodes);
    if(iZ >= fNsamples) printf("iZ > fNsamples %d %d\n",iZ,fNsamples);
    *Nlist=4;
    Xlist[0]=Xlist[1]=iX;
    if(iX && (iX != fNanodes/2)) Xlist[2]=iX-1;
    else Xlist[2]=iX;
    if ((iX !=fNanodes/2 -1) && (iX != fNanodes)) Xlist[3]=iX+1;
    else Xlist[3]=iX;
    if(iZ) Zlist[0]=iZ-1;
    else Zlist[0]=iZ;
    if (iZ < fNsamples) Zlist[1]=iZ+1;
    else Zlist[1]=iZ;
    Zlist[2]=Zlist[3]=iZ;
}
//------------------------------
void AliITSsegmentationSDD::GetPadIxz(Float_t x,Float_t z,Int_t &timebin,Int_t &anode){
//  Returns cell coordinates (time sample,anode) for given real local coordinates (x,z)

    // expects x, z in cm

    const Float_t kconv=10000;  // cm->um

    Float_t speed=fResponse->DriftSpeed();
    Int_t na = fNanodes/2;
    Float_t driftpath=fDx-TMath::Abs(kconv*x);
    timebin=(Int_t)(driftpath/speed/fTimeStep);
    anode=(Int_t)(kconv*z/fPitch) + na/2;
    if (x > 0) anode += na;

    timebin+=1;
    anode+=1;

}

//------------------------------
void AliITSsegmentationSDD::GetPadCxz(Int_t timebin,Int_t anode,Float_t &x ,Float_t &z){
    // Transform from cell to real local coordinates
  
    // returns x, z in cm

    const Float_t kconv=10000;  // um->cm

    Float_t speed=fResponse->DriftSpeed();
    Int_t na = fNanodes/2;
    Float_t driftpath=(timebin+1)*fTimeStep*speed;
    if (anode >= na) x=(fDx-driftpath)/kconv;
    else x = -(fDx-driftpath)/kconv;
    if (anode >= na) anode-=na;
    z=((anode+1)*fPitch-fDz/2)/kconv;

}

//------------------------------
void AliITSsegmentationSDD::GetPadTxz(Float_t &x,Float_t &z){
    // Get anode and time bucket as floats - numbering from 0

    // expects x, z in cm

    const Float_t kconv=10000;  // cm->um

    //Float_t x0=x;
    Float_t speed=fResponse->DriftSpeed();
    //Int_t na = fNanodes/2;
    Float_t driftpath=fDx-TMath::Abs(kconv*x);
    x=driftpath/speed/fTimeStep;
    z=kconv*z/fPitch;
    // z=kconv*z/fPitch + (float)na/2;
    //if (x0 > 0) z += (float)na;

}
//------------------------------
void AliITSsegmentationSDD::GetLocal(Int_t module,Float_t *g ,Float_t *l){
  // returns local coordinates from global
    if(!fGeom) {
        return;
        //fGeom = ((AliITS*)gAlice->GetModule("ITS"))->GetITSgeom();
    }
    fGeom->GtoL(module,g,l);
}
//------------------------------
void AliITSsegmentationSDD::GetGlobal(Int_t module,Float_t *l ,Float_t *g){
  // return global coordinates from local
    if(!fGeom) {
        return;
        //fGeom = ((AliITS*)gAlice->GetModule("ITS"))->GetITSgeom();
    }

    fGeom->LtoG(module,l,g);

}
