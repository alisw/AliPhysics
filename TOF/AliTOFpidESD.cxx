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

#include <TGeant3.h>              // For DtoM and InitGeo, should  
#include "../STRUCT/AliBODY.h"    // be removed after we switch to 
#include "../STRUCT/AliFRAMEv2.h" // AliTOGeometry::GetPos
#include "AliTOFv4T0.h"

#include <stdlib.h>


ClassImp(AliTOFpidESD)

//_________________________________________________________________________
static Int_t InitGeo() {
  //gSystem->Load("libgeant321");
  //new TGeant3("C++ Interface to Geant3");

  AliBODY *body = new AliBODY("BODY", "Alice envelop");
  body->CreateGeometry();

  AliFRAMEv2 *frame = new AliFRAMEv2("FRAME", "Space Frame");
  frame->SetHoles(1);
  frame->CreateGeometry();

  AliTOF *tof = new AliTOFv4T0("TOF", "TOF with Holes");
  tof->CreateGeometry();

  return 0;
} 

//_________________________________________________________________________
static Int_t DtoM(Int_t *dig, Float_t *g) {
  const Int_t kMAX=13;
  Int_t lnam[kMAX],lnum[kMAX];

  TGeant3 *geant3=(TGeant3*)gMC;
  if (!geant3) {::Info("DtoM","no geant3 found !"); return 1;}

  strncpy((Char_t*)(lnam+0),"ALIC",4); lnum[0]=1;
  strncpy((Char_t*)(lnam+1),"B077",4); lnum[1]=1;

//4 padx  if z<=0 then ...
if ((dig[1]==4)||(dig[1]==3)) dig[4]=AliTOFGeometry::NpadX()-dig[4];
else if (dig[1]==2) {
   if (dig[2]>7) dig[4]=AliTOFGeometry::NpadX()-dig[4];
   else if (dig[2]==7) {
      if (dig[3]==1) dig[4]=AliTOFGeometry::NpadX()-dig[4];
      else dig[4]+=1; 
   } else dig[4]+=1; 
} else dig[4]+=1;

//3 padz
if ((dig[1]==3)||(dig[1]==4)) dig[3]=AliTOFGeometry::NpadZ()-dig[3];
else dig[3]+=1;

//2 strip
if (dig[1]==0) dig[2]=AliTOFGeometry::NStripC()-dig[2]; 
else if (dig[1]==1) dig[2]=AliTOFGeometry::NStripB()-dig[2];
else dig[2]+=1; 

//1 module
dig[1]=AliTOFGeometry::NPlates()-dig[1];

//0 sector
dig[0]+=(AliTOFGeometry::NSectors()/2); dig[0]%=AliTOFGeometry::NSectors();
dig[0]+=1;

  //sector
  switch (dig[0]) {
     case 1:
         strncpy((Char_t*)(lnam+2),"B071",4); lnum[2]=9;
         strncpy((Char_t*)(lnam+3),"BTO1",4); lnum[3]=1;
         break;
     case 2:
         strncpy((Char_t*)(lnam+2),"B071",4); lnum[2]=10;
         strncpy((Char_t*)(lnam+3),"BTO1",4); lnum[3]=1;
         break;
     case 3:
         strncpy((Char_t*)(lnam+2),"B074",4); lnum[2]=4;
         strncpy((Char_t*)(lnam+3),"BTO2",4); lnum[3]=1;
         if (dig[1]>3) lnum[3]=2;
         break;
     case 4:
         strncpy((Char_t*)(lnam+2),"B074",4); lnum[2]=5;
         strncpy((Char_t*)(lnam+3),"BTO2",4); lnum[3]=1;
         if (dig[1]>3) lnum[3]=2;
         break;
     case 5:
         strncpy((Char_t*)(lnam+2),"B074",4); lnum[2]=1;
         strncpy((Char_t*)(lnam+3),"BTO2",4); lnum[3]=1;
         if (dig[1]>3) lnum[3]=2;
         break;
     case 6:
         strncpy((Char_t*)(lnam+2),"B074",4); lnum[2]=2;
         strncpy((Char_t*)(lnam+3),"BTO2",4); lnum[3]=1;
         if (dig[1]>3) lnum[3]=2;
         break;
     case 7:
         strncpy((Char_t*)(lnam+2),"B074",4); lnum[2]=3;
         strncpy((Char_t*)(lnam+3),"BTO2",4); lnum[3]=1;
         if (dig[1]>3) lnum[3]=2;
         break;
     case 8:
         strncpy((Char_t*)(lnam+2),"B071",4); lnum[2]=1;
         strncpy((Char_t*)(lnam+3),"BTO1",4); lnum[3]=1;
         break;
     case 9:
         strncpy((Char_t*)(lnam+2),"B071",4); lnum[2]=2;
         strncpy((Char_t*)(lnam+3),"BTO1",4); lnum[3]=1;
         break;
     case 10:
         strncpy((Char_t*)(lnam+2),"B075",4); lnum[2]=1;
         strncpy((Char_t*)(lnam+3),"BTO3",4); lnum[3]=1;
         if (dig[1]>4) lnum[3]=2;
         break;
     case 11:
         strncpy((Char_t*)(lnam+2),"B075",4); lnum[2]=2;
         strncpy((Char_t*)(lnam+3),"BTO3",4); lnum[3]=1;
         if (dig[1]>4) lnum[3]=2;
         break;
     case 12:
         strncpy((Char_t*)(lnam+2),"B075",4); lnum[2]=3;
         strncpy((Char_t*)(lnam+3),"BTO3",4); lnum[3]=1;
         if (dig[1]>4) lnum[3]=2;
         break;
     case 13:
         strncpy((Char_t*)(lnam+2),"B071",4); lnum[2]=3;
         strncpy((Char_t*)(lnam+3),"BTO1",4); lnum[3]=1;
         break;
     case 14:
         strncpy((Char_t*)(lnam+2),"B071",4); lnum[2]=4;
         strncpy((Char_t*)(lnam+3),"BTO1",4); lnum[3]=1;
         break;
     case 15:
         strncpy((Char_t*)(lnam+2),"B071",4); lnum[2]=5;
         strncpy((Char_t*)(lnam+3),"BTO1",4); lnum[3]=1;
         break;
     case 16:
         strncpy((Char_t*)(lnam+2),"B071",4); lnum[2]=6;
         strncpy((Char_t*)(lnam+3),"BTO1",4); lnum[3]=1;
         break;
     case 17:
         strncpy((Char_t*)(lnam+2),"B071",4); lnum[2]=7;
         strncpy((Char_t*)(lnam+3),"BTO1",4); lnum[3]=1;
         break;
     case 18:
         strncpy((Char_t*)(lnam+2),"B071",4); lnum[2]=8;
         strncpy((Char_t*)(lnam+3),"BTO1",4); lnum[3]=1;
         break;
     default:
         ::Info("DtoM","Wrong sector number : %d !",dig[0]);
         return 2;
  }
  //module
  switch (dig[1]) {
     case 1:
        strncpy((Char_t*)(lnam+4),"FTOC",4); lnum[4]=1;
        if (dig[0]==1  || dig[0]==2 ) lnum[4]=2;
        if (dig[0]>=8  && dig[0]<=9) lnum[4]=2;
        if (dig[0]>=13 && dig[0]<=18) lnum[4]=2;
        strncpy((Char_t*)(lnam+5),"FLTC",4); lnum[5]=0;
        break;
     case 2:
        strncpy((Char_t*)(lnam+4),"FTOB",4); lnum[4]=1;
        if (dig[0]==1  || dig[0]==2 ) lnum[4]=2;
        if (dig[0]>=8  && dig[0]<=9) lnum[4]=2;
        if (dig[0]>=13 && dig[0]<=18) lnum[4]=2;
        strncpy((Char_t*)(lnam+5),"FLTB",4); lnum[5]=0;
        break;
     case 3:
        strncpy((Char_t*)(lnam+4),"FTOA",4); lnum[4]=0;
        strncpy((Char_t*)(lnam+5),"FLTA",4); lnum[5]=0;
        break;
     case 4:
        strncpy((Char_t*)(lnam+4),"FTOB",4); lnum[4]=1;
        strncpy((Char_t*)(lnam+5),"FLTB",4); lnum[5]=0;
        break;
     case 5:
        strncpy((Char_t*)(lnam+4),"FTOC",4); lnum[4]=1;
        strncpy((Char_t*)(lnam+5),"FLTC",4); lnum[5]=0;
        break;
     default:
        ::Info("DtoM","Wrong module number : %d !",dig[1]);
        return 3;
  }

  //strip
  if (dig[2]<0/*1*/ || dig[2]>19) {
     ::Info("DtoM","Wrong strip number : %d !",dig[2]);
     return 3;
  }
  strncpy((Char_t*)(lnam+6),"FSTR",4); lnum[6]=dig[2];
  strncpy((Char_t*)(lnam+7),"FSEN",4); lnum[7]=0;

  //divz
  if (dig[3]<1 || dig[3]>2) {
     ::Info("DtoM","Wrong z-division number : %d !",dig[3]);
     return 3;
  }
  strncpy((Char_t*)(lnam+8),"FSEZ",4); lnum[8]=dig[3];

  //divx
  if (dig[4]<1 || dig[4]>48) {
     ::Info("DtoM","Wrong x-division number : %d !",dig[4]);
     return 3;
  }
  strncpy((Char_t*)(lnam+9),"FSEX",4); lnum[9]=dig[4];
  strncpy((Char_t*)(lnam+10),"FPAD",4); lnum[10]=0;

  Gcvolu_t *gcvolu=geant3->Gcvolu(); gcvolu->nlevel=0;
  Int_t err=geant3->Glvolu(11,lnam,lnum);  //11-th level
  if (err) {
    ::Info("DtoM","Wrong volume !"); return 3;
  }
  Float_t l[3]={0.,0.,0}; geant3->Gdtom(l,g,1);
  return 0;
}


//_________________________________________________________________________
AliTOFpidESD::AliTOFpidESD(Double_t *param) throw (const Char_t *) {
  //
  //  The main constructor
  //
  if (InitGeo()) throw "AliTOFpidESD: can not initialize the geometry !\n";


  if (fTOFGeometry) delete fTOFGeometry;
  fTOFGeometry = new AliTOFGeometry();

  fR=378.; 
  fDy=fTOFGeometry->XPad(); fDz=fTOFGeometry->ZPad(); 
  fN=0; fEventN=0;

  fSigma=param[0];
  fRange=param[1];

}


//_________________________________________________________________________

//_________________________________________________________________________
Int_t AliTOFpidESD::LoadClusters(const TFile *df) {

  //--------------------------------------------------------------------
  //This function loads the TOF clusters
  //--------------------------------------------------------------------

  if (!((TFile *)df)->IsOpen()) {
    Error("LoadClusters","file with the TOF digits has not been open !\n");
    return 1;
  }

  Char_t   name[100]; 
  sprintf(name,"TreeD%d",GetEventNumber());
  TTree *dTree=(TTree*)((TFile *)df)->Get(name);
  if (!dTree) { 
    Error("LoadClusters"," can't get the tree with digits !\n");
    return 1;
  }
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

    // fTOFGeometry->GetPos(dig,g);   // uncomment this 
    DtoM(dig,g);

    Double_t h[5];
    h[0]=TMath::Sqrt(g[0]*g[0]+g[1]*g[1]);
    h[1]=TMath::ATan2(g[1],g[0]); h[2]=g[2]; 
    h[3]=d->GetTdc(); h[4]=d->GetAdc();

    AliTOFcluster *cl=new AliTOFcluster(h,d->GetTracks(),i);
    InsertCluster(cl);
  }  

  delete dTree;
  return 0;
}

//_________________________________________________________________________
Int_t AliTOFpidESD::LoadClusters(TTree *dTree) {
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

    //fTOFGeometry->GetPos(dig,g);
    DtoM(dig,g);

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
      if (dph>TMath::Pi()) dph-=TMath::Pi();
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
	p[j]=TMath::Exp(-0.5*fRange*fRange);
        continue;
      }
      p[j]=TMath::Exp(-0.5*(tof-time[j])*(tof-time[j])/(sigma*sigma));
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

