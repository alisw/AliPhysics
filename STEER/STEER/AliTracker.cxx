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

#include <TMath.h>
#include <TH1F.h>
#include <TGeoMatrix.h>

#include "AliTracker.h"
#include "AliGeomManager.h"
#include "AliCluster.h"
#include "AliKalmanTrack.h"
#include "AliGlobalQADataMaker.h"

Bool_t AliTracker::fFillResiduals=kFALSE;
TObjArray **AliTracker::fResiduals=NULL;
AliRecoParam::EventSpecie_t AliTracker::fEventSpecie=AliRecoParam::kDefault;

ClassImp(AliTracker)

AliTracker::AliTracker():
  AliTrackerBase(),
  fEventInfo(NULL)
{
  //--------------------------------------------------------------------
  // The default constructor.
  //--------------------------------------------------------------------
}

//__________________________________________________________________________
AliTracker::AliTracker(const AliTracker &atr):
  AliTrackerBase(atr),
  fEventInfo(atr.fEventInfo)
{
  //--------------------------------------------------------------------
  // The default constructor.
  //--------------------------------------------------------------------
}

//__________________________________________________________________________
void AliTracker::FillClusterArray(TObjArray* /*array*/) const
{
  // Publishes all pointers to clusters known to the tracker into the
  // passed object array.
  // The ownership is not transfered - the caller is not expected to delete
  // the clusters.

  AliWarning("should be overriden by a sub-class.");
}

//__________________________________________________________________________
void AliTracker::CookLabel(AliKalmanTrack *t, Float_t wrong) const {
  //--------------------------------------------------------------------
  //This function "cooks" a track label. If label<0, this track is fake.
  //--------------------------------------------------------------------
  Int_t noc=t->GetNumberOfClusters();
  if (noc<1) return;
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
    if (j<noc) {
       lb[j]=lab;
       (mx[j])++;
    }
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

void AliTracker::FillResiduals(const AliExternalTrackParam *t,
			      Double_t *p, Double_t *cov, 
                              UShort_t id, Bool_t updated) {
  //
  // This function fills the histograms of residuals 
  // The array of these histos is external for this AliTracker class.
  // Normally, this array belong to AliGlobalQADataMaker class.  
  // 
  if (!fFillResiduals) return; 
  if (!fResiduals) return; 

  const Double_t *residuals=t->GetResiduals(p,cov,updated);
  if (!residuals) return;

  TH1F *h=0;
  Int_t esIndex = AliRecoParam::AConvert(fEventSpecie) ; 
  AliGeomManager::ELayerID layer=AliGeomManager::VolUIDToLayer(id);
  h=(TH1F*)fResiduals[esIndex]->At(2*layer-2);
  if (h) h->Fill(residuals[0]);
  h=(TH1F*)fResiduals[esIndex]->At(2*layer-1);
  if (h) h->Fill(residuals[1]);

  if (layer==5) {
    if (p[1]<0) {  // SSD1 absolute residuals
      h = (TH1F*)fResiduals[esIndex]->At(40);
      if (h) h->Fill(t->GetY()-p[0]); //C side
      h = (TH1F*)fResiduals[esIndex]->At(41);
      if (h) h->Fill(t->GetZ()-p[1]);
    } else {             
      h = (TH1F*)fResiduals[esIndex]->At(42);
      if (h) h->Fill(t->GetY()-p[0]); //A side
      h = (TH1F*)fResiduals[esIndex]->At(43);
      if (h) h->Fill(t->GetZ()-p[1]);
    }           
  }
  if (layer==6) {  // SSD2 absolute residuals
    if (p[1]<0) {
      h = (TH1F*)fResiduals[esIndex]->At(44);
      if (h) h->Fill(t->GetY()-p[0]); //C side
      h = (TH1F*)fResiduals[esIndex]->At(45);
      if (h) h->Fill(t->GetZ()-p[1]);
    } else {
      h = (TH1F*)fResiduals[esIndex]->At(46);
      if (h) h->Fill(t->GetY()-p[0]); //A side
      h = (TH1F*)fResiduals[esIndex]->At(47);
      if (h) h->Fill(t->GetZ()-p[1]);
    }
  }

}

void AliTracker::FillResiduals(const AliExternalTrackParam *t,
                               const AliCluster *c, Bool_t /*updated*/) {
  //
  // This function fills the histograms of residuals 
  // The array of these histos is external for this AliTracker class.
  // Normally, this array belong to AliGlobalQADataMaker class.  
  // 
  // For the moment, the residuals are absolute !
  //

  if (!fFillResiduals) return; 
  if (!fResiduals) return; 

  UShort_t id=c->GetVolumeId();
  const TGeoHMatrix *matrixT2L=AliGeomManager::GetTracking2LocalMatrix(id);

  // Position of the cluster in the tracking c.s.
  Double_t clsTrk[3]={c->GetX(), c->GetY(), c->GetZ()};
  // Position of the cluster in the local module c.s.
  Double_t clsLoc[3]={0.,0.,0.};
  matrixT2L->LocalToMaster(clsTrk,clsLoc);


  // Position of the intersection point in the tracking c.s.
  Double_t trkTrk[3]={t->GetX(),t->GetY(),t->GetZ()};
  // Position of the intersection point in the local module c.s.
  Double_t trkLoc[3]={0.,0.,0.};
  matrixT2L->LocalToMaster(trkTrk,trkLoc);

  Double_t residuals[2]={trkLoc[0]-clsLoc[0], trkLoc[2]-clsLoc[2]};

  TH1F *h=0;
  Int_t esIndex = AliRecoParam::AConvert(fEventSpecie) ; 
  AliGeomManager::ELayerID layer=AliGeomManager::VolUIDToLayer(id);
  h=(TH1F*)fResiduals[esIndex]->At(2*layer-2);
  if (h) h->Fill(residuals[0]);
  h=(TH1F*)fResiduals[esIndex]->At(2*layer-1);
  if (h) h->Fill(residuals[1]);

}

