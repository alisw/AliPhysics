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
//                             Origin
//               Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-------------------------------------------------------------------------

#include <TMath.h>

#include "AliTracker.h"
#include "AliCluster.h"
#include "AliKalmanTrack.h"
#include "AliRun.h"
#include "AliMagF.h"

#include "TFile.h"
#include "Riostream.h"


extern AliRun* gAlice;

ClassImp(AliTracker)

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
    if (TMath::Abs(c->GetLabel(1)) == lab ||
        TMath::Abs(c->GetLabel(2)) == lab ) max++;
  }

  if ((1.- Float_t(max)/noc) > wrong) lab=-lab;
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

////////////////////////////////////////////////////////////////////////
Int_t AliTracker::SetFieldFactor() {
//
// Utility class to set the value of the magnetic field in the barrel
// It supposes that the correct object gAlice is in the memory
//
   AliKalmanTrack::SetConvConst(1000/0.299792458/gAlice->Field()->SolenoidField());
   cout<<"Magnetic field in kGauss: "<<gAlice->Field()->SolenoidField()<<endl;
   return 0;
}
////////////////////////////////////////////////////////////////////////
Int_t AliTracker::SetFieldFactor(TFile *file, Bool_t deletegAlice) {
//
// Utility class to set the value of the magnetic field in the barrel
// gAlice object is read from the file, and optionally deleted
// 
  if (!(gAlice=(AliRun*)file->Get("gAlice"))) {
    cerr<<"gAlice has not been found in file "<<file->GetName();
    return 1;
  }   
  Int_t rc = SetFieldFactor();
  if (deletegAlice) {
    delete gAlice;  
    gAlice = 0;
  }
  return rc;
}
////////////////////////////////////////////////////////////////////////
Int_t AliTracker::SetFieldFactor(Char_t* fileName, Bool_t closeFile) {
//
// Utility class to set the value of the magnetic field in the barrel
// gAlice object is read from the file, the file is optionally closed
// 
   TFile *file=TFile::Open(fileName);
   if (!file->IsOpen()) {cerr<<"Cannnot open "<<fileName<<" !\n"; return 1;}
   Int_t rc = SetFieldFactor(file, closeFile) ;
   if (closeFile) file->Close();
   return rc;
}
////////////////////////////////////////////////////////////////////////

