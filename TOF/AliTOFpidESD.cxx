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

//-----------------------------------------------------------------
//           Implementation of the TOF PID class
// Very naive one... Should be made better by the detector experts...
//      Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
//-----------------------------------------------------------------
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TError.h"

#include "AliTOFpidESD.h"
#include "AliESD.h"
#include "AliESDtrack.h"
#include "AliTOFdigit.h"
#include "AliTOFGeometry.h"

#include <stdlib.h>


ClassImp(AliTOFpidESD)

//_________________________________________________________________________
AliTOFpidESD::AliTOFpidESD(Double_t *param) {
  //
  //  The main constructor
  //
  fR=378.; 
  fDy=AliTOFGeometry::XPad(); fDz=AliTOFGeometry::ZPad(); 
  fN=0; fEventN=0;

  fSigma=param[0];
  fRange=param[1];

}

//_________________________________________________________________________
Int_t AliTOFpidESD::LoadClusters(TTree *dTree, AliTOFGeometry *geom) {
  //--------------------------------------------------------------------
  //This function loads the TOF clusters
  //--------------------------------------------------------------------
  TBranch *branch=dTree->GetBranch("TOF");
  if (!branch) { 
    Error("LoadClusters"," can't get the branch with the TOF digits !\n");
    return 1;
  }

  TClonesArray dummy("AliTOFdigit",10000), *digits=&dummy;
  branch->SetAddress(&digits);

  dTree->GetEvent(0);
  Int_t nd=digits->GetEntriesFast();
  Info("LoadClusters","number of digits: %d",nd);

  for (Int_t i=0; i<nd; i++) {
    AliTOFdigit *d=(AliTOFdigit*)digits->UncheckedAt(i);
    Int_t dig[5]; Float_t g[3];
    dig[0]=d->GetSector();
    dig[1]=d->GetPlate();
    dig[2]=d->GetStrip();
    dig[3]=d->GetPadz();
    dig[4]=d->GetPadx();

    geom->GetPos(dig,g);

    Double_t h[5];
    h[0]=TMath::Sqrt(g[0]*g[0]+g[1]*g[1]);
    h[1]=TMath::ATan2(g[1],g[0]); h[2]=g[2]; 
    h[3]=d->GetTdc(); h[4]=d->GetAdc();

    AliTOFcluster *cl=new AliTOFcluster(h,d->GetTracks(),i);
    InsertCluster(cl);
  }  

  return 0;
}

//_________________________________________________________________________
void AliTOFpidESD::UnloadClusters() {
  //--------------------------------------------------------------------
  //This function unloads TOF clusters
  //--------------------------------------------------------------------
  for (Int_t i=0; i<fN; i++) delete fClusters[i];
  fN=0;
}

//_________________________________________________________________________
Int_t AliTOFpidESD::InsertCluster(AliTOFcluster *c) {
  //--------------------------------------------------------------------
  //This function adds a cluster to the array of clusters sorted in Z
  //--------------------------------------------------------------------
  if (fN==kMaxCluster) {
    Error("InsertCluster","Too many clusters !\n");
    return 1;
  }

  if (fN==0) {fClusters[fN++]=c; return 0;}
  Int_t i=FindClusterIndex(c->GetZ());
  memmove(fClusters+i+1 ,fClusters+i,(fN-i)*sizeof(AliTOFcluster*));
  fClusters[i]=c; fN++;

  return 0;
}

//_________________________________________________________________________
Int_t AliTOFpidESD::FindClusterIndex(Double_t z) const {
  //--------------------------------------------------------------------
  // This function returns the index of the nearest cluster 
  //--------------------------------------------------------------------
  if (fN==0) return 0;
  if (z <= fClusters[0]->GetZ()) return 0;
  if (z > fClusters[fN-1]->GetZ()) return fN;
  Int_t b=0, e=fN-1, m=(b+e)/2;
  for (; b<e; m=(b+e)/2) {
    if (z > fClusters[m]->GetZ()) b=m+1;
    else e=m; 
  }
  return m;
}

static int cmp(const void *p1, const void *p2) {
  AliESDtrack *t1=*((AliESDtrack**)p1);
  AliESDtrack *t2=*((AliESDtrack**)p2);
  Double_t c1[15]; t1->GetExternalCovariance(c1);
  Double_t c2[15]; t2->GetExternalCovariance(c2);
  if (c1[0]*c1[2] <c2[0]*c2[2]) return -1;
  if (c1[0]*c1[2]==c2[0]*c2[2]) return 0;
  return 1;
}

//_________________________________________________________________________
Int_t AliTOFpidESD::MakePID(AliESD *event)
{
  //
  //  This function calculates the "detector response" PID probabilities
  //                Just for a bare hint... 

  static const Double_t kMasses[]={
    0.000511, 0.105658, 0.139570, 0.493677, 0.938272, 1.875613
  };

  Int_t ntrk=event->GetNumberOfTracks();
  AliESDtrack **tracks=new AliESDtrack*[ntrk];

  Int_t i;
  for (i=0; i<ntrk; i++) {
    AliESDtrack *t=event->GetTrack(i);
    tracks[i]=t;
  }
  qsort(tracks,ntrk,sizeof(AliESDtrack*),cmp);

  Int_t nmatch=0;
  for (i=0; i<ntrk; i++) {
    AliESDtrack *t=tracks[i];

    if ((t->GetStatus()&AliESDtrack::kTRDout)==0) continue;
    if ((t->GetStatus()&AliESDtrack::kTRDStop)!=0) continue;

    Double_t x,par[5]; t->GetExternalParameters(x,par);
    Double_t cov[15]; t->GetExternalCovariance(cov);

Double_t dphi=(5*TMath::Sqrt(cov[0]) + 0.5*fDy + 2.5*TMath::Abs(par[2]))/fR; 
Double_t dz=5*TMath::Sqrt(cov[2]) + 0.5*fDz + 2.5*TMath::Abs(par[3]);

    Double_t phi=TMath::ATan2(par[0],x) + t->GetAlpha();
    if (phi<-TMath::Pi()) phi+=2*TMath::Pi();
    if (phi>=TMath::Pi()) phi-=2*TMath::Pi();
    Double_t z=par[1];

    Double_t d2max=1000.;
    Int_t index=-1;
    for (Int_t k=FindClusterIndex(z-dz); k<fN; k++) {
      AliTOFcluster *c=fClusters[k];

      if (c->GetZ() > z+dz) break;
      if (c->IsUsed()) continue;

      Double_t dph=TMath::Abs(c->GetPhi()-phi);
      if (dph>TMath::Pi()) dph-=2*TMath::Pi();
      if (dph>dphi) continue;

      Double_t d2=dph*dph*fR*fR + (c->GetZ()-z)*(c->GetZ()-z);
      if (d2 > d2max) continue;

      d2max=d2;
      index=k;
    }

    if (index<0) {
      //Info("MakePID","matching failed ! %d",TMath::Abs(t->GetLabel()));
       continue;
    }

    nmatch++;

    AliTOFcluster *c=fClusters[index];
    c->Use();

    Double_t tof=50*c->GetTDC()+32; // in ps
    t->SetTOFsignal(tof);
    t->SetTOFcluster(c->GetIndex());

    if ((t->GetStatus()&AliESDtrack::kTIME)==0) continue;

    Double_t time[10]; t->GetIntegratedTimes(time);

    //track length correction
    Double_t rc=TMath::Sqrt(c->GetR()*c->GetR() + c->GetZ()*c->GetZ());
    Double_t rt=TMath::Sqrt(x*x + par[0]*par[0] + par[1]*par[1]);
    Double_t dlt=rc-rt;

    Double_t p[10];
    Double_t mom=t->GetP();
    for (Int_t j=0; j<AliESDtrack::kSPECIES; j++) {
      Double_t mass=kMasses[j];
      Double_t dpp=0.01;      //mean relative pt resolution;
      if (mom>0.5) dpp=0.01*mom;
      Double_t sigma=dpp*time[j]/(1.+ mom*mom/(mass*mass));
      sigma=TMath::Sqrt(sigma*sigma + fSigma*fSigma);

      time[j]+=dlt/3e-2*TMath::Sqrt(mom*mom+mass*mass)/mom;

      if (TMath::Abs(tof-time[j]) > fRange*sigma) {
	p[j]=TMath::Exp(-0.5*fRange*fRange)/sigma;
        continue;
      }
      p[j]=TMath::Exp(-0.5*(tof-time[j])*(tof-time[j])/(sigma*sigma))/sigma;
    }
    /*
    if (c->GetLabel(0)!=TMath::Abs(t->GetLabel())) {
       cerr<<"Wrong matching: "<<t->GetLabel()<<endl;
       continue;
    }
    */
    t->SetTOFpid(p);

  }

  Info("MakePID","Number of matched ESD track: %d",nmatch);

  delete[] tracks;

  return 0;
}

