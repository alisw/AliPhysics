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

//-------------------------------------------------------
//          Implementation of the TPC tracker
//
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------

#include <TObjArray.h>
#include <TFile.h>
#include <TTree.h>
#include <AliRunLoader.h>
#include <AliLoader.h>
#include <Riostream.h>

#include "AliTPCtracker.h"
#include "AliTPCcluster.h"
#include "AliTPCParam.h"
#include "AliTPCClustersRow.h"
#include "AliTPCcluster.h"

//_____________________________________________________________________________
AliTPCtracker::AliTPCtracker(const AliTPCParam *par, Int_t eventn, const char* evfoldname):
AliTracker(), fkNIS(par->GetNInnerSector()/2), fkNOS(par->GetNOuterSector()/2),
fEvFolderName(evfoldname)
{
  //---------------------------------------------------------------------
  // The main TPC tracker constructor
  //---------------------------------------------------------------------
  cout<<"fkNIS = "<<fkNIS<<endl;
  cout<<"fkNOS = "<<fkNOS<<endl;

  fInnerSec=new AliTPCSector[fkNIS];         
  fOuterSec=new AliTPCSector[fkNOS];

  Int_t i;
  
  for (i=0; i<fkNIS; i++) fInnerSec[i].Setup(par,0);
  for (i=0; i<fkNOS; i++) fOuterSec[i].Setup(par,1);

  fN=0;  fSectors=0;


  fEventN = eventn;
  AliRunLoader* rl = AliRunLoader::GetRunLoader(fEvFolderName);
  if (rl == 0x0)
   {
     Error("AliTPCtracker","Can not get RL from specified folder %s",fEvFolderName.Data());
     return;
   }
  rl->GetEvent(fEventN);

  AliLoader* tpcl = rl->GetLoader("TPCLoader");
  if (tpcl == 0x0)
   {
     Error("AliTPCtracker","Can not get TPC Laoder from Run Loader");
     return;
   }
  
  if (tpcl->TreeR() == 0x0) tpcl->LoadRecPoints("read");
  TTree* treeR = tpcl->TreeR();
  if (treeR == 0x0)
   {
     cout<<"Error: Can not get TreeR\n";
   }
  fClustersArray.Setup(par);
  fClustersArray.SetClusterType("AliTPCcluster");
  fClustersArray.ConnectTree(treeR);

  fSeeds=0;
}

//_____________________________________________________________________________
AliTPCtracker::~AliTPCtracker() {
  //------------------------------------------------------------------
  // TPC tracker destructor
  //------------------------------------------------------------------
  delete[] fInnerSec;
  delete[] fOuterSec;
  if (fSeeds) {
    fSeeds->Delete(); 
    delete fSeeds;
  }
}

//_____________________________________________________________________________
Double_t SigmaY2(Double_t r, Double_t tgl, Double_t pt)
{
  //
  // Parametrised error of the cluster reconstruction (pad direction)   
  //
  // Sigma rphi
  const Float_t kArphi=0.41818e-2;
  const Float_t kBrphi=0.17460e-4;
  const Float_t kCrphi=0.30993e-2;
  const Float_t kDrphi=0.41061e-3;
  
  pt=TMath::Abs(pt)*1000.;
  Double_t x=r/pt;
  tgl=TMath::Abs(tgl);
  Double_t s=kArphi - kBrphi*r*tgl + kCrphi*x*x + kDrphi*x;
  if (s<0.4e-3) s=0.4e-3;
  s*=1.3; //Iouri Belikov

  return s;
}

//_____________________________________________________________________________
Double_t SigmaZ2(Double_t r, Double_t tgl) 
{
  //
  // Parametrised error of the cluster reconstruction (drift direction)
  //
  // Sigma z
  const Float_t kAz=0.39614e-2;
  const Float_t kBz=0.22443e-4;
  const Float_t kCz=0.51504e-1;
  

  tgl=TMath::Abs(tgl);
  Double_t s=kAz - kBz*r*tgl + kCz*tgl*tgl;
  if (s<0.4e-3) s=0.4e-3;
  s*=1.3; //Iouri Belikov

  return s;
}

//_____________________________________________________________________________
Double_t f1(Double_t x1,Double_t y1,
                   Double_t x2,Double_t y2,
                   Double_t x3,Double_t y3) 
{
  //-----------------------------------------------------------------
  // Initial approximation of the track curvature
  //-----------------------------------------------------------------
  Double_t d=(x2-x1)*(y3-y2)-(x3-x2)*(y2-y1);
  Double_t a=0.5*((y3-y2)*(y2*y2-y1*y1+x2*x2-x1*x1)-
                  (y2-y1)*(y3*y3-y2*y2+x3*x3-x2*x2));
  Double_t b=0.5*((x2-x1)*(y3*y3-y2*y2+x3*x3-x2*x2)-
                  (x3-x2)*(y2*y2-y1*y1+x2*x2-x1*x1));

  Double_t xr=TMath::Abs(d/(d*x1-a)), yr=d/(d*y1-b);

  return -xr*yr/sqrt(xr*xr+yr*yr); 
}


//_____________________________________________________________________________
Double_t f2(Double_t x1,Double_t y1,
                   Double_t x2,Double_t y2,
                   Double_t x3,Double_t y3) 
{
  //-----------------------------------------------------------------
  // Initial approximation of the track curvature times center of curvature
  //-----------------------------------------------------------------
  Double_t d=(x2-x1)*(y3-y2)-(x3-x2)*(y2-y1);
  Double_t a=0.5*((y3-y2)*(y2*y2-y1*y1+x2*x2-x1*x1)-
                  (y2-y1)*(y3*y3-y2*y2+x3*x3-x2*x2));
  Double_t b=0.5*((x2-x1)*(y3*y3-y2*y2+x3*x3-x2*x2)-
                  (x3-x2)*(y2*y2-y1*y1+x2*x2-x1*x1));

  Double_t xr=TMath::Abs(d/(d*x1-a)), yr=d/(d*y1-b);
  
  return -a/(d*y1-b)*xr/sqrt(xr*xr+yr*yr);
}

//_____________________________________________________________________________
Double_t f3(Double_t x1,Double_t y1, 
                   Double_t x2,Double_t y2,
                   Double_t z1,Double_t z2) 
{
  //-----------------------------------------------------------------
  // Initial approximation of the tangent of the track dip angle
  //-----------------------------------------------------------------
  return (z1 - z2)/sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

//_____________________________________________________________________________
void AliTPCtracker::LoadOuterSectors() {
  //-----------------------------------------------------------------
  // This function fills outer TPC sectors with clusters.
  //-----------------------------------------------------------------
  UInt_t index;
  TTree* tree = fClustersArray.GetTree();
  if (tree == 0x0)
   {
    Error("LoadOuterSectors","Can not get tree from fClustersArray");
    return;
   }
  Int_t j=Int_t(tree->GetEntries());
  cout<<"fClustersArray.GetTree()->GetEntries() = "<<j<<endl;
  for (Int_t i=0; i<j; i++) {
      AliSegmentID *s=fClustersArray.LoadEntry(i);
      Int_t sec,row;
      AliTPCParam *par=(AliTPCParam*)fClustersArray.GetParam();
      par->AdjustSectorRow(s->GetID(),sec,row);
      if (sec<fkNIS*2) continue;
      AliTPCClustersRow *clrow=fClustersArray.GetRow(sec,row);
      Int_t ncl=clrow->GetArray()->GetEntriesFast();
      while (ncl--) {
         AliTPCcluster *c=(AliTPCcluster*)(*clrow)[ncl];
         index=(((sec<<8)+row)<<16)+ncl;
         fOuterSec[(sec-fkNIS*2)%fkNOS][row].InsertCluster(c,index);
      }
  }

  fN=fkNOS;
  fSectors=fOuterSec;
}

//_____________________________________________________________________________
void AliTPCtracker::UnloadOuterSectors() {
  //-----------------------------------------------------------------
  // This function clears outer TPC sectors.
  //-----------------------------------------------------------------
  Int_t nup=fOuterSec->GetNRows();
  for (Int_t i=0; i<fkNOS; i++) {
    for (Int_t j=0; j<nup; j++) {
      if (fClustersArray.GetRow(i+fkNIS*2,j)) 
          fClustersArray.ClearRow(i+fkNIS*2,j);
      if (fClustersArray.GetRow(i+fkNIS*2+fkNOS,j)) 
          fClustersArray.ClearRow(i+fkNIS*2+fkNOS,j);
    }
  }

  fN=0;
  fSectors=0;
}

//_____________________________________________________________________________
void AliTPCtracker::LoadInnerSectors() {
  //-----------------------------------------------------------------
  // This function fills inner TPC sectors with clusters.
  //-----------------------------------------------------------------
  UInt_t index;
  Int_t j=Int_t(fClustersArray.GetTree()->GetEntries());
  for (Int_t i=0; i<j; i++) {
      AliSegmentID *s=fClustersArray.LoadEntry(i);
      Int_t sec,row;
      AliTPCParam *par=(AliTPCParam*)fClustersArray.GetParam();
      par->AdjustSectorRow(s->GetID(),sec,row);
      if (sec>=fkNIS*2) continue;
      AliTPCClustersRow *clrow=fClustersArray.GetRow(sec,row);
      Int_t ncl=clrow->GetArray()->GetEntriesFast();
      while (ncl--) {
         AliTPCcluster *c=(AliTPCcluster*)(*clrow)[ncl];
         index=(((sec<<8)+row)<<16)+ncl;
         fInnerSec[sec%fkNIS][row].InsertCluster(c,index);
      }
  }

  fN=fkNIS;
  fSectors=fInnerSec;
}

//_____________________________________________________________________________
void AliTPCtracker::UnloadInnerSectors() {
  //-----------------------------------------------------------------
  // This function clears inner TPC sectors.
  //-----------------------------------------------------------------
  Int_t nlow=fInnerSec->GetNRows();
  for (Int_t i=0; i<fkNIS; i++) {
    for (Int_t j=0; j<nlow; j++) {
      if (fClustersArray.GetRow(i,j)) fClustersArray.ClearRow(i,j);
      if (fClustersArray.GetRow(i+fkNIS,j)) fClustersArray.ClearRow(i+fkNIS,j);
    }
  }

  fN=0;
  fSectors=0;
}

//_____________________________________________________________________________
Int_t AliTPCtracker::FollowProlongation(AliTPCseed& t, Int_t rf) {
  //-----------------------------------------------------------------
  // This function tries to find a track prolongation.
  //-----------------------------------------------------------------
  Double_t xt=t.GetX();
  const Int_t kSKIP=(t.GetNumberOfClusters()<10) ? kRowsToSkip : 
                    Int_t(0.5*fSectors->GetNRows());
  Int_t tryAgain=kSKIP;

  Double_t alpha=t.GetAlpha() - fSectors->GetAlphaShift();
  if (alpha > 2.*TMath::Pi()) alpha -= 2.*TMath::Pi();  
  if (alpha < 0.            ) alpha += 2.*TMath::Pi();  
  Int_t s=Int_t(alpha/fSectors->GetAlpha())%fN;

  Int_t nrows=fSectors->GetRowNumber(xt)-1;
  for (Int_t nr=nrows; nr>=rf; nr--) {
    Double_t x=fSectors->GetX(nr), ymax=fSectors->GetMaxY(nr);
    if (!t.PropagateTo(x)) return 0;

    AliTPCcluster *cl=0;
    UInt_t index=0;
    Double_t maxchi2=kMaxCHI2;
    const AliTPCRow &krow=fSectors[s][nr];
    Double_t pt=t.GetConvConst()/(100/0.299792458/0.2)/t.Get1Pt();
    Double_t sy2=SigmaY2(t.GetX(),t.GetTgl(),pt);
    Double_t sz2=SigmaZ2(t.GetX(),t.GetTgl());
    Double_t road=4.*sqrt(t.GetSigmaY2() + sy2), y=t.GetY(), z=t.GetZ();

    if (road>kMaxROAD) {
      if (t.GetNumberOfClusters()>4) 
        cerr<<t.GetNumberOfClusters()
        <<"FindProlongation warning: Too broad road !\n"; 
      return 0;
    }

    if (krow) {
      for (Int_t i=krow.Find(y-road); i<krow; i++) {
       AliTPCcluster *c=(AliTPCcluster*)(krow[i]);
       if (c->GetY() > y+road) break;
       if (c->IsUsed()) continue;
       if ((c->GetZ()-z)*(c->GetZ()-z) > 16.*(t.GetSigmaZ2()+sz2)) continue;
       Double_t chi2=t.GetPredictedChi2(c);
       if (chi2 > maxchi2) continue;
       maxchi2=chi2;
       cl=c;
        index=krow.GetIndex(i);       
      }
    }
    if (cl) {
      Float_t l=fSectors->GetPadPitchWidth();
      Float_t corr=1.; if (nr>63) corr=0.67; // new (third) pad response !
      t.SetSampledEdx(cl->GetQ()/l*corr,t.GetNumberOfClusters());
      if (!t.Update(cl,maxchi2,index)) {
         if (!tryAgain--) return 0;
      } else tryAgain=kSKIP;
    } else {
      if (tryAgain==0) break;
      if (y > ymax) {
         s = (s+1) % fN;
         if (!t.Rotate(fSectors->GetAlpha())) return 0;
      } else if (y <-ymax) {
         s = (s-1+fN) % fN;
         if (!t.Rotate(-fSectors->GetAlpha())) return 0;
      }
      tryAgain--;
    }
  }

  return 1;
}

//_____________________________________________________________________________
Int_t AliTPCtracker::FollowBackProlongation
(AliTPCseed& seed, const AliTPCtrack &track) {
  //-----------------------------------------------------------------
  // This function propagates tracks back through the TPC
  //-----------------------------------------------------------------
  Double_t alpha=seed.GetAlpha() - fSectors->GetAlphaShift();
  if (alpha > 2.*TMath::Pi()) alpha -= 2.*TMath::Pi();  
  if (alpha < 0.            ) alpha += 2.*TMath::Pi();  
  Int_t s=Int_t(alpha/fSectors->GetAlpha())%fN;

  Int_t idx=-1, sec=-1, row=-1;
  Int_t nc=seed.GetLabel(); //index of the cluster to start with
  if (nc--) {
     idx=track.GetClusterIndex(nc);
     sec=(idx&0xff000000)>>24; row=(idx&0x00ff0000)>>16;
  }
  if (fSectors==fInnerSec) { if (sec >= 2*fkNIS) row=-1; } 
  else { if (sec <  2*fkNIS) row=-1; }   

  Int_t nr=fSectors->GetNRows();
  for (Int_t i=0; i<nr; i++) {
    Double_t x=fSectors->GetX(i), ymax=fSectors->GetMaxY(i);

    if (!seed.PropagateTo(x)) return 0;

    Double_t y=seed.GetY();
    if (y > ymax) {
       s = (s+1) % fN;
       if (!seed.Rotate(fSectors->GetAlpha())) return 0;
    } else if (y <-ymax) {
       s = (s-1+fN) % fN;
       if (!seed.Rotate(-fSectors->GetAlpha())) return 0;
    }

    AliTPCcluster *cl=0;
    Int_t index=0;
    Double_t maxchi2=kMaxCHI2;
    Double_t pt=seed.GetConvConst()/(100/0.299792458/0.2)/seed.Get1Pt();
    Double_t sy2=SigmaY2(seed.GetX(),seed.GetTgl(),pt);
    Double_t sz2=SigmaZ2(seed.GetX(),seed.GetTgl());
    Double_t road=4.*sqrt(seed.GetSigmaY2() + sy2), z=seed.GetZ();
    if (road>kMaxROAD) {
      cerr<<seed.GetNumberOfClusters()
          <<"AliTPCtracker::FollowBackProlongation: Too broad road !\n"; 
      return 0;
    }


    Int_t accepted=seed.GetNumberOfClusters();
    if (row==i) {
       //try to accept already found cluster
       AliTPCcluster *c=(AliTPCcluster*)GetCluster(idx);
       Double_t chi2;
       if ((chi2=seed.GetPredictedChi2(c))<maxchi2 || accepted<27) { 
          index=idx; cl=c; maxchi2=chi2;
       } //else cerr<<"AliTPCtracker::FollowBackProlongation: oulier !\n";
          
       if (nc--) {
          idx=track.GetClusterIndex(nc); 
          sec=(idx&0xff000000)>>24; row=(idx&0x00ff0000)>>16;
       } 
       if (fSectors==fInnerSec) { if (sec >= 2*fkNIS) row=-1; }
       else { if (sec <  2*fkNIS) row=-1; }   

    }
    if (!cl) {
       //try to fill the gap
       const AliTPCRow &krow=fSectors[s][i];
       if (accepted>27)
       if (krow) {
          for (Int_t i=krow.Find(y-road); i<krow; i++) {
           AliTPCcluster *c=(AliTPCcluster*)(krow[i]);
           if (c->GetY() > y+road) break;
           if (c->IsUsed()) continue;
        if ((c->GetZ()-z)*(c->GetZ()-z)>16.*(seed.GetSigmaZ2()+sz2)) continue;
           Double_t chi2=seed.GetPredictedChi2(c);
           if (chi2 > maxchi2) continue;
           maxchi2=chi2;
           cl=c;
            index=krow.GetIndex(i);
          }
       }
    }

    if (cl) {
      Float_t l=fSectors->GetPadPitchWidth();
      Float_t corr=1.; if (i>63) corr=0.67; // new (third) pad response !
      seed.SetSampledEdx(cl->GetQ()/l*corr,seed.GetNumberOfClusters());
      seed.Update(cl,maxchi2,index);
    }

  }

  seed.SetLabel(nc);

  return 1;
}

//_____________________________________________________________________________
void AliTPCtracker::MakeSeeds(Int_t i1, Int_t i2) {
  //-----------------------------------------------------------------
  // This function creates track seeds.
  //-----------------------------------------------------------------
  cout<<" Making Seeds i1="<<i1<<"  i2="<<i2<<"\n";
  if (fSeeds==0) fSeeds=new TObjArray(15000);

  Double_t x[5], c[15];

  Double_t alpha=fOuterSec->GetAlpha(), shift=fOuterSec->GetAlphaShift();
  Double_t cs=cos(alpha), sn=sin(alpha);

  Double_t x1 =fOuterSec->GetX(i1);
  Double_t xx2=fOuterSec->GetX(i2);

  for (Int_t ns=0; ns<fkNOS; ns++) {
    Int_t nl=fOuterSec[(ns-1+fkNOS)%fkNOS][i2];
    Int_t nm=fOuterSec[ns][i2];
    Int_t nu=fOuterSec[(ns+1)%fkNOS][i2];
    const AliTPCRow& kr1=fOuterSec[ns][i1];
    for (Int_t is=0; is < kr1; is++) {
      Double_t y1=kr1[is]->GetY(), z1=kr1[is]->GetZ();
      for (Int_t js=0; js < nl+nm+nu; js++) {
       const AliTPCcluster *kcl;
        Double_t x2,   y2,   z2;
        Double_t x3=GetX(), y3=GetY(), z3=GetZ();

       if (js<nl) {
         const AliTPCRow& kr2=fOuterSec[(ns-1+fkNOS)%fkNOS][i2];
         kcl=kr2[js];
          y2=kcl->GetY(); z2=kcl->GetZ();
          x2= xx2*cs+y2*sn;
          y2=-xx2*sn+y2*cs;
       } else 
         if (js<nl+nm) {
           const AliTPCRow& kr2=fOuterSec[ns][i2];
           kcl=kr2[js-nl];
            x2=xx2; y2=kcl->GetY(); z2=kcl->GetZ();
         } else {
           const AliTPCRow& kr2=fOuterSec[(ns+1)%fkNOS][i2];
           kcl=kr2[js-nl-nm];
            y2=kcl->GetY(); z2=kcl->GetZ();
            x2=xx2*cs-y2*sn;
            y2=xx2*sn+y2*cs;
         }

        Double_t zz=z1 - (z1-z3)/(x1-x3)*(x1-x2); 
        if (TMath::Abs(zz-z2)>5.) continue;

        Double_t d=(x2-x1)*(0.-y2)-(0.-x2)*(y2-y1);
        if (d==0.) {cerr<<"MakeSeeds warning: Straight seed !\n"; continue;}

       x[0]=y1;
       x[1]=z1;
       x[4]=f1(x1,y1,x2,y2,x3,y3);
       if (TMath::Abs(x[4]) >= 0.0066) continue;
       x[2]=f2(x1,y1,x2,y2,x3,y3);
       //if (TMath::Abs(x[4]*x1-x[2]) >= 0.99999) continue;
       x[3]=f3(x1,y1,x2,y2,z1,z2);
       if (TMath::Abs(x[3]) > 1.2) continue;
       Double_t a=asin(x[2]);
       Double_t zv=z1 - x[3]/x[4]*(a+asin(x[4]*x1-x[2]));
       if (TMath::Abs(zv-z3)>10.) continue; 

        Double_t sy1=kr1[is]->GetSigmaY2(), sz1=kr1[is]->GetSigmaZ2();
        Double_t sy2=kcl->GetSigmaY2(),     sz2=kcl->GetSigmaZ2();
	//Double_t sy3=400*3./12., sy=0.1, sz=0.1;
        Double_t sy3=25000*x[4]*x[4]+0.1, sy=0.1, sz=0.1;

       Double_t f40=(f1(x1,y1+sy,x2,y2,x3,y3)-x[4])/sy;
       Double_t f42=(f1(x1,y1,x2,y2+sy,x3,y3)-x[4])/sy;
       Double_t f43=(f1(x1,y1,x2,y2,x3,y3+sy)-x[4])/sy;
       Double_t f20=(f2(x1,y1+sy,x2,y2,x3,y3)-x[2])/sy;
       Double_t f22=(f2(x1,y1,x2,y2+sy,x3,y3)-x[2])/sy;
       Double_t f23=(f2(x1,y1,x2,y2,x3,y3+sy)-x[2])/sy;
       Double_t f30=(f3(x1,y1+sy,x2,y2,z1,z2)-x[3])/sy;
       Double_t f31=(f3(x1,y1,x2,y2,z1+sz,z2)-x[3])/sz;
       Double_t f32=(f3(x1,y1,x2,y2+sy,z1,z2)-x[3])/sy;
       Double_t f34=(f3(x1,y1,x2,y2,z1,z2+sz)-x[3])/sz;

        c[0]=sy1;
        c[1]=0.;       c[2]=sz1;
        c[3]=f20*sy1;  c[4]=0.;       c[5]=f20*sy1*f20+f22*sy2*f22+f23*sy3*f23;
        c[6]=f30*sy1;  c[7]=f31*sz1;  c[8]=f30*sy1*f20+f32*sy2*f22;
                       c[9]=f30*sy1*f30+f31*sz1*f31+f32*sy2*f32+f34*sz2*f34;
        c[10]=f40*sy1; c[11]=0.; c[12]=f40*sy1*f20+f42*sy2*f22+f43*sy3*f23;
        c[13]=f30*sy1*f40+f32*sy2*f42;
        c[14]=f40*sy1*f40+f42*sy2*f42+f43*sy3*f43;

        UInt_t index=kr1.GetIndex(is);
       AliTPCseed *track=new AliTPCseed(index, x, c, x1, ns*alpha+shift);
        Float_t l=fOuterSec->GetPadPitchWidth();
        track->SetSampledEdx(kr1[is]->GetQ()/l,0);

        Int_t rc=FollowProlongation(*track, i2);
        if (rc==0 || track->GetNumberOfClusters()<(i1-i2)/2) delete track;
        else 
         {
           fSeeds->AddLast(track); 
           cout<<"Adding seed   "<<fSeeds->GetEntries()<<"\r";
         }
      }
    }
  }
}

//_____________________________________________________________________________
Int_t AliTPCtracker::ReadSeeds(const TFile *inp) {
  //-----------------------------------------------------------------
  // This function reades track seeds.
  //-----------------------------------------------------------------
  TDirectory *savedir=gDirectory; 

  TFile *in=(TFile*)inp;
  if (!in->IsOpen()) {
     cerr<<"AliTPCtracker::ReadSeeds(): input file is not open !\n";
     return 1;
  }

  in->cd();
  TTree *seedTree=(TTree*)in->Get("Seeds");
  if (!seedTree) {
     cerr<<"AliTPCtracker::ReadSeeds(): ";
     cerr<<"can't get a tree with track seeds !\n";
     return 2;
  }
  AliTPCtrack *seed=new AliTPCtrack; 
  seedTree->SetBranchAddress("tracks",&seed);
  
  if (fSeeds==0) fSeeds=new TObjArray(15000);

  Int_t n=(Int_t)seedTree->GetEntries();
  for (Int_t i=0; i<n; i++) {
     seedTree->GetEvent(i);
     fSeeds->AddLast(new AliTPCseed(*seed,seed->GetAlpha()));
  }
  
  delete seed;

  delete seedTree; //Thanks to Mariana Bondila

  savedir->cd();

  return 0;
}

//_____________________________________________________________________________
Int_t AliTPCtracker::Clusters2Tracks() 
{
  //-----------------------------------------------------------------
  // This is a track finder.
  //-----------------------------------------------------------------
  Int_t retval = 0;

  AliRunLoader* rl = AliRunLoader::GetRunLoader(fEvFolderName);
  if (rl == 0x0)
   {
     Error("Clusters2Tracks","Can not get RL from specified folder %s",fEvFolderName.Data());
     return 1;
   }
  
  retval = rl->GetEvent(fEventN);
  if (retval)
   {
     Error("Clusters2Tracks","Error while getting event %d",fEventN);
     return 1;
   }
  
  AliLoader* tpcl = rl->GetLoader("TPCLoader");
  if (tpcl == 0x0)
   {
     Error("Clusters2Tracks","Can not get TPC Laoder from Run Loader");
     return 1;
   }
 
  if ( tpcl->TreeR() == 0x0)
   { 
    retval = tpcl->LoadRecPoints("READ");
    if (retval)
     {
       Error("Clusters2Tracks","Error while loading Reconstructed Points");
       return 1;
     }
   }
   
  if (tpcl->TreeT() == 0x0) tpcl->MakeTree("T");

  TTree &tracktree = *(tpcl->TreeT());
  
  TBranch* br= tracktree.GetBranch("tracks");
  if (br)
   {
     Error("Clusters2Tracks","Branch \"tracks\" already exists in TreeT for TPC");
     return 1;
   }
  
  AliTPCtrack *iotrack=0;
  tracktree.Branch("tracks","AliTPCtrack",&iotrack,32000,0);

  LoadOuterSectors();

  //find track seeds
  Int_t nup=fOuterSec->GetNRows(), nlow=fInnerSec->GetNRows();
  Int_t nrows=nlow+nup;
  if (fSeeds==0) {
     Int_t gap=Int_t(0.125*nrows), shift=Int_t(0.5*gap);
     MakeSeeds(nup-1, nup-1-gap);
     MakeSeeds(nup-1-shift, nup-1-shift-gap);
  }    
  fSeeds->Sort();

  //tracking in outer sectors
  Int_t nseed=fSeeds->GetEntriesFast();
  cout<<"nseed="<<nseed<<endl;
  Int_t i;
  for (i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)fSeeds->UncheckedAt(i), &t=*pt;
    if (FollowProlongation(t)) {
       UseClusters(&t);
       continue;
    }
    delete fSeeds->RemoveAt(i);
  }  
  //UnloadOuterSectors();

  //tracking in inner sectors
  LoadInnerSectors();
  Int_t found=0;
  for (i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)fSeeds->UncheckedAt(i), &t=*pt;
    if (!pt) continue;
    Int_t nc=t.GetNumberOfClusters();

    Double_t alpha=t.GetAlpha() - fInnerSec->GetAlphaShift();
    if (alpha > 2.*TMath::Pi()) alpha -= 2.*TMath::Pi();
    if (alpha < 0.            ) alpha += 2.*TMath::Pi();
    Int_t ns=Int_t(alpha/fInnerSec->GetAlpha())%fkNIS;

    alpha=ns*fInnerSec->GetAlpha() + fInnerSec->GetAlphaShift() - t.GetAlpha();

    if (t.Rotate(alpha)) {
       if (FollowProlongation(t)) {
          if (t.GetNumberOfClusters() >= Int_t(0.4*nrows)) {
             t.CookdEdx();
             CookLabel(pt,0.1); //For comparison only
             iotrack=pt;
             tracktree.Fill();
             UseClusters(&t,nc);
             cerr<<found++<<'\r';
          }
       }
    }
    delete fSeeds->RemoveAt(i); 
  }  
  UnloadInnerSectors();
  UnloadOuterSectors();

  tpcl->WriteTracks("OVERWRITE");
  
  cerr<<"Number of found tracks : "<<found<<endl;

  return 0;
}

//_____________________________________________________________________________
Int_t AliTPCtracker::PropagateBack() 
 {
  //-----------------------------------------------------------------
  // This function propagates tracks back through the TPC.
  //-----------------------------------------------------------------
  
  cout<<"This method is not converted to NewIO yet\n";
  return 1;

  fSeeds=new TObjArray(15000);
  Int_t retval = 0;

  AliRunLoader* rl = AliRunLoader::GetRunLoader(fEvFolderName);
  if (rl == 0x0)
   {
     Error("Clusters2Tracks","Can not get RL from specified folder %s",fEvFolderName.Data());
     return 1;
   }
  
  retval = rl->GetEvent(fEventN);
  if (retval)
   {
     Error("Clusters2Tracks","Error while getting event %d",fEventN);
     return 1;
   }
  
  AliLoader* tpcl = rl->GetLoader("TPCLoader");
  if (tpcl == 0x0)
   {
     Error("Clusters2Tracks","Can not get TPC Laoder from Run Loader");
     return 1;
   }

  AliLoader* itsl = rl->GetLoader("ITSLoader");
  if (tpcl == 0x0)
   {
     Error("Clusters2Tracks","Can not get ITS Laoder from Run Loader");
     return 1;
   }
  
  if (itsl->TreeT() == 0x0) itsl->LoadTracks();

  TTree *bckTree=itsl->TreeT();
  if (!bckTree) {
     cerr<<"AliTPCtracker::PropagateBack() ";
     cerr<<"can't get a tree with back propagated ITS tracks !\n";
     return 3;
  }
  
  AliTPCtrack *bckTrack=new AliTPCtrack; 
  bckTree->SetBranchAddress("tracks",&bckTrack);


  TFile* in = 0x0;
  TFile* out = 0x0;
  cout<<"And NOW there will be a segmentation violation!!!!\n";
  bckTree=(TTree*)in->Get("TreeT_ITSb_0");
  if (!bckTree) {
     cerr<<"AliTPCtracker::PropagateBack() ";
     cerr<<"can't get a tree with back propagated ITS tracks !\n";
     return 3;
  }


  TTree *tpcTree=(TTree*)in->Get("TreeT_TPC_0");
  if (!tpcTree) {
     cerr<<"AliTPCtracker::PropagateBack() ";
     cerr<<"can't get a tree with TPC tracks !\n";
     return 4;
  }
  AliTPCtrack *tpcTrack=new AliTPCtrack; 
  tpcTree->SetBranchAddress("tracks",&tpcTrack);

//*** Prepare an array of tracks to be back propagated
  Int_t nup=fOuterSec->GetNRows(), nlow=fInnerSec->GetNRows();
  Int_t nrows=nlow+nup;

  TObjArray tracks(15000);
  Int_t i=0,j=0;
  Int_t tpcN=(Int_t)tpcTree->GetEntries();
  Int_t bckN=(Int_t)bckTree->GetEntries();
  if (j<bckN) bckTree->GetEvent(j++);
  for (i=0; i<tpcN; i++) {
     tpcTree->GetEvent(i);
     Double_t alpha=tpcTrack->GetAlpha();
     if (j<bckN &&
     TMath::Abs(tpcTrack->GetLabel())==TMath::Abs(bckTrack->GetLabel())) {
        if (!bckTrack->Rotate(alpha-bckTrack->GetAlpha())) continue;
        fSeeds->AddLast(new AliTPCseed(*bckTrack,bckTrack->GetAlpha()));
        bckTree->GetEvent(j++);
     } else {
        tpcTrack->ResetCovariance();
        fSeeds->AddLast(new AliTPCseed(*tpcTrack,alpha));
     }
     tracks.AddLast(new AliTPCtrack(*tpcTrack));
  }

  out->cd();
  TTree backTree("TreeT_TPCb_0","Tree with back propagated TPC tracks");
  AliTPCtrack *otrack=0;
  backTree.Branch("tracks","AliTPCtrack",&otrack,32000,0);

//*** Back propagation through inner sectors
  LoadInnerSectors();
  Int_t nseed=fSeeds->GetEntriesFast();
  for (i=0; i<nseed; i++) {
     AliTPCseed *ps=(AliTPCseed*)fSeeds->UncheckedAt(i), &s=*ps;
     const AliTPCtrack *pt=(AliTPCtrack*)tracks.UncheckedAt(i), &t=*pt;

     Int_t nc=t.GetNumberOfClusters();
     s.SetLabel(nc-1); //set number of the cluster to start with

     if (FollowBackProlongation(s,t)) {
       UseClusters(&s);
        continue;
     }
     delete fSeeds->RemoveAt(i);
  }  
  UnloadInnerSectors();

//*** Back propagation through outer sectors
  LoadOuterSectors();
  Int_t found=0;
  for (i=0; i<nseed; i++) {
     AliTPCseed *ps=(AliTPCseed*)fSeeds->UncheckedAt(i), &s=*ps;
     if (!ps) continue;
     Int_t nc=s.GetNumberOfClusters();
     const AliTPCtrack *pt=(AliTPCtrack*)tracks.UncheckedAt(i), &t=*pt;

     Double_t alpha=s.GetAlpha() - fSectors->GetAlphaShift();
     if (alpha > 2.*TMath::Pi()) alpha -= 2.*TMath::Pi();
     if (alpha < 0.            ) alpha += 2.*TMath::Pi();
     Int_t ns=Int_t(alpha/fSectors->GetAlpha())%fN;

     alpha =ns*fSectors->GetAlpha() + fSectors->GetAlphaShift();
     alpha-=s.GetAlpha();

     if (s.Rotate(alpha)) {
        if (FollowBackProlongation(s,t)) {
          if (s.GetNumberOfClusters() >= Int_t(0.4*nrows)) {
             s.CookdEdx();
              s.SetLabel(t.GetLabel());
              UseClusters(&s,nc);
              otrack=ps;
              backTree.Fill();
              cerr<<found++<<'\r';
              continue;
          }
       }
     }
     delete fSeeds->RemoveAt(i);
  }  
  UnloadOuterSectors();

  backTree.Write();

  cerr<<"Number of seeds: "<<nseed<<endl;
  cerr<<"Number of back propagated ITS tracks: "<<bckN<<endl;
  cerr<<"Number of back propagated TPC tracks: "<<found<<endl;

  delete bckTrack;
  delete tpcTrack;

  delete bckTree; //Thanks to Mariana Bondila
  delete tpcTree; //Thanks to Mariana Bondila

  return 0;
}

//_________________________________________________________________________
AliCluster *AliTPCtracker::GetCluster(Int_t index) const {
  //--------------------------------------------------------------------
  //       Return pointer to a given cluster
  //--------------------------------------------------------------------
  Int_t sec=(index&0xff000000)>>24; 
  Int_t row=(index&0x00ff0000)>>16; 
  Int_t ncl=(index&0x0000ffff)>>00;

  AliTPCClustersRow *clrow=((AliTPCtracker *) this)->fClustersArray.GetRow(sec,row);
  return (AliCluster*)(*clrow)[ncl];      
}

//__________________________________________________________________________
void AliTPCtracker::CookLabel(AliKalmanTrack *t, Float_t wrong) const {
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

  else {
     Int_t tail=Int_t(0.10*noc);
     max=0;
     for (i=1; i<=tail; i++) {
       AliCluster *c=clusters[noc-i];
       if (lab == TMath::Abs(c->GetLabel(0)) ||
           lab == TMath::Abs(c->GetLabel(1)) ||
           lab == TMath::Abs(c->GetLabel(2))) max++;
     }
     if (max < Int_t(0.5*tail)) lab=-lab;
  }

  t->SetLabel(lab);

  delete[] lb;
  delete[] mx;
  delete[] clusters;
}

//_________________________________________________________________________
void AliTPCtracker::AliTPCSector::Setup(const AliTPCParam *par, Int_t f) {
  //-----------------------------------------------------------------------
  // Setup inner sector
  //-----------------------------------------------------------------------
  if (f==0) {
     fAlpha=par->GetInnerAngle();
     fAlphaShift=par->GetInnerAngleShift();
     fPadPitchWidth=par->GetInnerPadPitchWidth();
     f1PadPitchLength=par->GetInnerPadPitchLength();
     f2PadPitchLength=f1PadPitchLength;

     fN=par->GetNRowLow();
//     cout<<"par->GetNRowLow() = "<<par->GetNRowLow()<<"  fN = "<<fN<<endl;
     fRow=new AliTPCRow[fN];
     for (Int_t i=0; i<fN; i++) fRow[i].SetX(par->GetPadRowRadiiLow(i));
  } else {
     fAlpha=par->GetOuterAngle();
     fAlphaShift=par->GetOuterAngleShift();
     fPadPitchWidth=par->GetOuterPadPitchWidth();
     f1PadPitchLength=par->GetOuter1PadPitchLength();
     f2PadPitchLength=par->GetOuter2PadPitchLength();

     fN=par->GetNRowUp();
     fRow=new AliTPCRow[fN];
     for (Int_t i=0; i<fN; i++){ 
     fRow[i].SetX(par->GetPadRowRadiiUp(i));
     }
  } 
}

//_________________________________________________________________________
void 
AliTPCtracker::AliTPCRow::InsertCluster(const AliTPCcluster* c, UInt_t index) {
  //-----------------------------------------------------------------------
  // Insert a cluster into this pad row in accordence with its y-coordinate
  //-----------------------------------------------------------------------
  if (fN==kMaxClusterPerRow) {
    cerr<<"AliTPCRow::InsertCluster(): Too many clusters !\n"; return;
  }
  if (fN==0) {fIndex[0]=index; fClusters[fN++]=c; return;}
  Int_t i=Find(c->GetY());
  memmove(fClusters+i+1 ,fClusters+i,(fN-i)*sizeof(AliTPCcluster*));
  memmove(fIndex   +i+1 ,fIndex   +i,(fN-i)*sizeof(UInt_t));
  fIndex[i]=index; fClusters[i]=c; fN++;
}

//___________________________________________________________________
Int_t AliTPCtracker::AliTPCRow::Find(Double_t y) const {
  //-----------------------------------------------------------------------
  // Return the index of the nearest cluster 
  //-----------------------------------------------------------------------
  if(fN<=0) return 0;
  if (y <= fClusters[0]->GetY()) return 0;
  if (y > fClusters[fN-1]->GetY()) return fN;
  Int_t b=0, e=fN-1, m=(b+e)/2;
  for (; b<e; m=(b+e)/2) {
    if (y > fClusters[m]->GetY()) b=m+1;
    else e=m; 
  }
  return m;
}

//_____________________________________________________________________________
void AliTPCtracker::AliTPCseed::CookdEdx(Double_t low, Double_t up) {
  //-----------------------------------------------------------------
  // This funtion calculates dE/dX within the "low" and "up" cuts.
  //-----------------------------------------------------------------
  Int_t i;
  Int_t nc=GetNumberOfClusters();
  Int_t * index = new Int_t[nc];
  TMath::Sort(nc, fdEdxSample,index,kFALSE);

  Int_t nl=Int_t(low*nc), nu=Int_t(up*nc);
  Float_t dedx=0;
  for (i=nl; i<=nu; i++) dedx += fdEdxSample[index[i]];
  dedx /= (nu-nl+1);
  SetdEdx(dedx);

  delete [] index;

  //Very rough PID
  Double_t p=TMath::Sqrt((1.+ GetTgl()*GetTgl())/(Get1Pt()*Get1Pt()));

  Double_t log1=TMath::Log(p+0.45), log2=TMath::Log(p+0.12);
  if (p<0.6) {
    if (dedx < 34 + 30/(p+0.45)/(p+0.45) + 24*log1) {SetMass(0.13957); return;}
    if (dedx < 34 + 30/(p+0.12)/(p+0.12) + 24*log2) {SetMass(0.49368); return;}
    SetMass(0.93827); return;
  }

  if (p<1.2) {
    if (dedx < 34 + 30/(p+0.12)/(p+0.12) + 24*log2) {SetMass(0.13957); return;}
    SetMass(0.93827); return;
  }

  SetMass(0.13957); return;

}



