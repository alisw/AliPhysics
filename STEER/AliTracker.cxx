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

//-------------------------------------------------------------------------
//               Implementation of the AliTracker class
//  that is the base for AliTPCtracker, AliITStrackerV2 and AliTRDtracker    
//        Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-------------------------------------------------------------------------
#include <TClass.h>
#include <TMath.h>

#include "AliMagF.h"
#include "AliTracker.h"
#include "AliCluster.h"
#include "AliKalmanTrack.h"

Bool_t AliTracker::fgUniformField=kTRUE;
Double_t AliTracker::fgBz=kAlmost0Field;
const AliMagF *AliTracker::fgkFieldMap=0;

ClassImp(AliTracker)

AliTracker::AliTracker():
  TObject(),
  fX(0),
  fY(0),
  fZ(0),
  fSigmaX(0.005),
  fSigmaY(0.005),
  fSigmaZ(0.010)
{
  //--------------------------------------------------------------------
  // The default constructor.
  //--------------------------------------------------------------------
  if (!fgkFieldMap) AliWarning("Field map is not set. Call AliTracker::SetFieldMap before creating a tracker!");
}

//__________________________________________________________________________
AliTracker::AliTracker(const AliTracker &atr):
  TObject(atr),
  fX(atr.fX),
  fY(atr.fY),
  fZ(atr.fZ),
  fSigmaX(atr.fSigmaX),
  fSigmaY(atr.fSigmaY),
  fSigmaZ(atr.fSigmaZ)
{
  //--------------------------------------------------------------------
  // The default constructor.
  //--------------------------------------------------------------------
  if (!fgkFieldMap) AliWarning("Field map is not set. Call AliTracker::SetFieldMap before creating a tracker!");
}

//__________________________________________________________________________
void AliTracker::SetFieldMap(const AliMagF* map, Bool_t uni) {
  //--------------------------------------------------------------------
  //This passes the field map to the reconstruction.
  //--------------------------------------------------------------------
  if (map==0) AliFatalClass("Can't access the field map !");

  if (fgkFieldMap) {
     AliWarningClass("The magnetic field map has been already set !");
     return;
  }

  fgUniformField=uni;
  fgkFieldMap=map;

  //Float_t r[3]={0.,0.,0.},b[3]; map->Field(r,b);
  //Double_t bz=-b[2];
 
  Double_t bz=-map->SolenoidField();
  fgBz=TMath::Sign(kAlmost0Field,bz) + bz;

}

//__________________________________________________________________________
void AliTracker::CookLabel(AliKalmanTrack *t, Float_t wrong) const {
  //--------------------------------------------------------------------
  //This function "cooks" a track label. If label<0, this track is fake.
  //--------------------------------------------------------------------
  Int_t noc=t->GetNumberOfClusters();
  Int_t *lb=new Int_t[noc];
  Int_t *mx=new Int_t[noc];
  AliCluster **clusters=new AliCluster*[noc];

  Int_t i;
  for (i=0; i<noc; i++) {
     lb[i]=mx[i]=0;
     Int_t index=t->GetClusterIndex(i);
     clusters[i]=GetCluster(index);
  }

  Int_t lab=123456789;
  for (i=0; i<noc; i++) {
    AliCluster *c=clusters[i];
    lab=TMath::Abs(c->GetLabel(0));
    Int_t j;
    for (j=0; j<noc; j++) if (lb[j]==lab || mx[j]==0) break;
    lb[j]=lab;
    (mx[j])++;
  }

  Int_t max=0;
  for (i=0; i<noc; i++) if (mx[i]>max) {max=mx[i]; lab=lb[i];}
    
  for (i=0; i<noc; i++) {
    AliCluster *c=clusters[i];
    //if (TMath::Abs(c->GetLabel(1)) == lab ||
    //    TMath::Abs(c->GetLabel(2)) == lab ) max++;
    if (TMath::Abs(c->GetLabel(0)!=lab))
	if (TMath::Abs(c->GetLabel(1)) == lab ||
	    TMath::Abs(c->GetLabel(2)) == lab ) max++;
  }

  if ((1.- Float_t(max)/noc) > wrong) lab=-lab;
  t->SetFakeRatio((1.- Float_t(max)/noc));
  t->SetLabel(lab);

  delete[] lb;
  delete[] mx;
  delete[] clusters;
}

//____________________________________________________________________________
void AliTracker::UseClusters(const AliKalmanTrack *t, Int_t from) const {
  //------------------------------------------------------------------
  //This function marks clusters associated with the track.
  //------------------------------------------------------------------
  Int_t noc=t->GetNumberOfClusters();
  for (Int_t i=from; i<noc; i++) {
     Int_t index=t->GetClusterIndex(i);
     AliCluster *c=GetCluster(index); 
     c->Use();   
  }
}

Double_t AliTracker::GetBz(Float_t *r) {
  //------------------------------------------------------------------
  // Returns Bz (kG) at the point "r" .
  //------------------------------------------------------------------
    Float_t b[3]; fgkFieldMap->Field(r,b);
    Double_t bz=-Double_t(b[2]);
    return  (TMath::Sign(kAlmost0Field,bz) + bz);
}
