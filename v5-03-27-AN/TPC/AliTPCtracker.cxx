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
#include <TError.h>
#include <TFile.h>
#include <TTree.h>

#include "AliESDEvent.h"
#include "AliESDtrack.h"

#include "AliTPCtracker.h"
#include "AliTPCcluster.h"
#include "AliTPCParam.h"
#include "AliClusters.h"

ClassImp(AliTPCtracker)

//_____________________________________________________________________________
AliTPCtracker::AliTPCtracker(): 
  AliTracker(), 
  fkNIS(0), 
  fInnerSec(0),         
  fkNOS(0),
  fOuterSec(0),
  fN(0),
  fSectors(0),
  fParam(0),
  fSeeds(0)
{
  //
  // The default TPC tracker constructor
  //
}

//_____________________________________________________________________________
AliTPCtracker::AliTPCtracker(const AliTPCParam *par): 
  AliTracker(), 
  fkNIS(par->GetNInnerSector()/2), 
  fInnerSec(new AliTPCSector[fkNIS]),         
  fkNOS(par->GetNOuterSector()/2),
  fOuterSec(new AliTPCSector[fkNOS]),
  fN(0),
  fSectors(0),
  fParam((AliTPCParam*) par),
  fSeeds(0)
{
  //---------------------------------------------------------------------
  // The main TPC tracker constructor
  //---------------------------------------------------------------------

  Int_t i;
  for (i=0; i<fkNIS; i++) fInnerSec[i].Setup(par,0);
  for (i=0; i<fkNOS; i++) fOuterSec[i].Setup(par,1);

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
Int_t AliTPCtracker::LoadClusters(TTree *cTree) {
  //-----------------------------------------------------------------
  // This function loads TPC clusters.
  //-----------------------------------------------------------------
  TBranch *branch=cTree->GetBranch("Segment");
  if (!branch) {
     Error("LoadClusters","Can't get the branch !");
     return 1;
  }

  AliClusters carray, *addr=&carray;
  addr = new AliClusters("AliTPCcluster");
  branch->SetAddress(&addr);

  Int_t nentr=(Int_t)cTree->GetEntries();

  for (Int_t i=0; i<nentr; i++) {
      cTree->GetEvent(i);

      Int_t ncl=carray.GetArray()->GetEntriesFast();

      Int_t nir=fInnerSec->GetNRows(), nor=fOuterSec->GetNRows();
      Int_t id=carray.GetID();
      if ((id<0) || (id>2*(fkNIS*nir + fkNOS*nor))) {
	 Fatal("LoadClusters","Wrong index !");
      }
      Int_t outindex = 2*fkNIS*nir;
      if (id<outindex) {
         Int_t sec = id/nir;
         Int_t row = id - sec*nir;
         sec %= fkNIS;
         AliTPCRow &padrow=fInnerSec[sec][row];
         while (ncl--) {
           AliTPCcluster *c=(AliTPCcluster*)carray[ncl];
           padrow.InsertCluster(c,sec,row);
         }
      } else {
         id -= outindex;
         Int_t sec = id/nor;
         Int_t row = id - sec*nor;
         sec %= fkNOS;
         AliTPCRow &padrow=fOuterSec[sec][row];
         while (ncl--) {
           AliTPCcluster *c=(AliTPCcluster*)carray[ncl];
           padrow.InsertCluster(c,sec+fkNIS,row);
         }
      }
      carray.GetArray()->Clear();
  }

  return 0;
}

//_____________________________________________________________________________
void AliTPCtracker::UnloadClusters() {
  //-----------------------------------------------------------------
  // This function unloads TPC clusters.
  //-----------------------------------------------------------------
  Int_t i;
  for (i=0; i<fkNIS; i++) {
    Int_t nr=fInnerSec->GetNRows();
    for (Int_t n=0; n<nr; n++) fInnerSec[i][n].ResetClusters();
  }
  for (i=0; i<fkNOS; i++) {
    Int_t nr=fOuterSec->GetNRows();
    for (Int_t n=0; n<nr; n++) fOuterSec[i][n].ResetClusters();
  }
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
    Double_t pt=t.GetSignedPt();
    Double_t sy2=AliTPCcluster::SigmaY2(t.GetX(),t.GetTgl(),pt);
    Double_t sz2=AliTPCcluster::SigmaZ2(t.GetX(),t.GetTgl());
    Double_t road=4.*sqrt(t.GetSigmaY2() + sy2), y=t.GetY(), z=t.GetZ();

    if (road>kMaxROAD) {
      if (t.GetNumberOfClusters()>4) 
         Warning("FindProlongation","Too broad road !"); 
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

Int_t AliTPCtracker::FollowRefitInward(AliTPCseed *seed, AliTPCtrack *track) {
  //
  // This function propagates seed inward TPC using old clusters
  // from the track.
  // 
  // Sylwester Radomski, GSI
  // 26.02.2003
  //

  // loop over rows

  Int_t nRows = fSectors->GetNRows();
  for (Int_t iRow = nRows-1; iRow >= 0; iRow--) {

    Double_t x = fSectors->GetX(iRow);
    if (!seed->PropagateTo(x)) return 0;

    // try to find an assigned cluster in this row

    AliTPCcluster* cluster = NULL;
    Int_t idx = -1;
    Int_t sec = -1;
    for (Int_t iCluster = 0; iCluster < track->GetNumberOfClusters(); iCluster++) {
      idx = track->GetClusterIndex(iCluster); 
      sec = (idx&0xff000000)>>24; 
      Int_t row = (idx&0x00ff0000)>>16;
      if (((fSectors == fInnerSec) && (sec >= fkNIS)) ||
	  ((fSectors == fOuterSec) && (sec < fkNIS))) continue;
      if (row == iRow) {
	cluster = (AliTPCcluster*) GetCluster(idx);
	break;
      }
    }

    // update the track seed with the found cluster

    if (cluster) {
      Double_t dAlpha = fParam->GetAngle(sec) - seed->GetAlpha();
      if (TMath::Abs(dAlpha) > 0.0001) {
	if (!seed->Rotate(dAlpha)) return 0;
	if (!seed->PropagateTo(x)) return 0;
      }

      seed->Update(cluster, seed->GetPredictedChi2(cluster), idx);
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
  Int_t nc=track.GetNumberOfClusters();

  if (nc--) {
     idx=track.GetClusterIndex(nc);
     sec=(idx&0xff000000)>>24; row=(idx&0x00ff0000)>>16;
  }
  if (fSectors==fInnerSec) { if (sec >= fkNIS) row=-1; } 
  else { if (sec <  fkNIS) row=-1; }   

  Int_t nr=fSectors->GetNRows();
  for (Int_t i=0; i<nr; i++) {
    Double_t x=fSectors->GetX(i), ymax=fSectors->GetMaxY(i);
    Double_t y;
    if (!seed.GetYAt(x,GetBz(),y)) return 0;
 
    if (y > ymax) {
       s = (s+1) % fN;
       if (!seed.Rotate(fSectors->GetAlpha())) return 0;
    } else if (y <-ymax) {
       s = (s-1+fN) % fN;
       if (!seed.Rotate(-fSectors->GetAlpha())) return 0;
    }

    if (!seed.PropagateTo(x)) return 0;

    AliTPCcluster *cl=0;
    Int_t index=0;
    Double_t maxchi2=kMaxCHI2;
    Double_t pt=seed.GetSignedPt();
    Double_t sy2=AliTPCcluster::SigmaY2(seed.GetX(),seed.GetTgl(),pt);
    Double_t sz2=AliTPCcluster::SigmaZ2(seed.GetX(),seed.GetTgl());
    Double_t road=4.*sqrt(seed.GetSigmaY2() + sy2), z=seed.GetZ();
    if (road>kMaxROAD) {
      Warning("FollowBackProlongation","Too broad road !"); 
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
       if (fSectors==fInnerSec) { if (sec >= fkNIS) row=-1; }
       else { if (sec < fkNIS) row=-1; }   

    }
    if (!cl) {
       //try to fill the gap
       const AliTPCRow &krow=fSectors[s][i];
       if (accepted>27)
       if (krow) {
          for (Int_t icl=krow.Find(y-road); icl<krow; icl++) {
	    AliTPCcluster *c=(AliTPCcluster*)(krow[icl]);
	    if (c->GetY() > y+road) break;
	    if (c->IsUsed()) continue;
	 if ((c->GetZ()-z)*(c->GetZ()-z)>16.*(seed.GetSigmaZ2()+sz2)) continue;
	    Double_t chi2=seed.GetPredictedChi2(c);
	    if (chi2 > maxchi2) continue;
	    maxchi2=chi2;
	    cl=c;
            index=krow.GetIndex(icl);
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

  return 1;
}

//_____________________________________________________________________________
void AliTPCtracker::MakeSeeds(Int_t i1, Int_t i2) {
  //-----------------------------------------------------------------
  // This function creates track seeds.
  //-----------------------------------------------------------------
  Double_t x[5], c[15];

  Double_t alpha=fSectors->GetAlpha(), shift=fSectors->GetAlphaShift();
  Double_t cs=cos(alpha), sn=sin(alpha);

  Double_t x1 =fSectors->GetX(i1);
  Double_t xx2=fSectors->GetX(i2);

  for (Int_t ns=0; ns<fN; ns++) {
    Int_t nl=fSectors[(ns-1+fN)%fN][i2];
    Int_t nm=fSectors[ns][i2];
    Int_t nu=fSectors[(ns+1)%fN][i2];
    const AliTPCRow& kr1=fSectors[ns][i1];
    for (Int_t is=0; is < kr1; is++) {
      Double_t y1=kr1[is]->GetY(), z1=kr1[is]->GetZ();
      for (Int_t js=0; js < nl+nm+nu; js++) {
	const AliTPCcluster *kcl;
        Double_t x2,   y2,   z2;
        Double_t x3=GetX(), y3=GetY(), z3=GetZ();

	if (js<nl) {
	  const AliTPCRow& kr2=fSectors[(ns-1+fN)%fN][i2];
	  kcl=kr2[js];
          y2=kcl->GetY(); z2=kcl->GetZ();
          x2= xx2*cs+y2*sn;
          y2=-xx2*sn+y2*cs;
	} else 
	  if (js<nl+nm) {
	    const AliTPCRow& kr2=fSectors[ns][i2];
	    kcl=kr2[js-nl];
            x2=xx2; y2=kcl->GetY(); z2=kcl->GetZ();
	  } else {
	    const AliTPCRow& kr2=fSectors[(ns+1)%fN][i2];
	    kcl=kr2[js-nl-nm];
            y2=kcl->GetY(); z2=kcl->GetZ();
            x2=xx2*cs-y2*sn;
            y2=xx2*sn+y2*cs;
	  }

        Double_t zz=z1 - (z1-z3)/(x1-x3)*(x1-x2); 
        if (TMath::Abs(zz-z2)>5.) continue;

        Double_t d=(x2-x1)*(0.-y2)-(0.-x2)*(y2-y1);
        if (d==0.) {
           Warning("MakeSeeds","Straight seed !"); 
           continue;
        }
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

        Int_t index=kr1.GetIndex(is);
	AliTPCseed *track=new AliTPCseed(x1, ns*alpha+shift, x, c, index);
        Float_t l=fSectors->GetPadPitchWidth();
        track->SetSampledEdx(kr1[is]->GetQ()/l,0);

        Int_t rc=FollowProlongation(*track, i2);
        if (rc==0 || track->GetNumberOfClusters()<(i1-i2)/2) delete track;
        else fSeeds->AddLast(track); 
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
     Error("ReadSeeds","Input file has not been open !");
     return 1;
  }

  in->cd();
  TTree *seedTree=(TTree*)in->Get("Seeds");
  if (!seedTree) {
     Error("ReadSeeds","Can't get a tree with track seeds !");
     return 2;
  }
  AliTPCtrack *seed=new AliTPCtrack; 
  seedTree->SetBranchAddress("tracks",&seed);
  
  if (fSeeds==0) fSeeds=new TObjArray(15000);

  Int_t n=(Int_t)seedTree->GetEntries();
  for (Int_t i=0; i<n; i++) {
     seedTree->GetEvent(i);
     seed->ResetClusters();
     fSeeds->AddLast(new AliTPCseed(*seed));
  }
  
  delete seed;

  delete seedTree; //Thanks to Mariana Bondila

  savedir->cd();

  return 0;
}

//_____________________________________________________________________________
Int_t AliTPCtracker::Clusters2Tracks(AliESDEvent *event) {
  //-----------------------------------------------------------------
  // This is a track finder.
  // The clusters must be already loaded ! 
  //-----------------------------------------------------------------

  //find track seeds
  Int_t nup=fOuterSec->GetNRows(), nlow=fInnerSec->GetNRows();
  Int_t nrows=nlow+nup;
  if (fSeeds==0) {
     Int_t gap=Int_t(0.125*nrows), shift=Int_t(0.5*gap);
     fSectors=fOuterSec; fN=fkNOS;
     fSeeds=new TObjArray(15000);
     MakeSeeds(nup-1, nup-1-gap);
     MakeSeeds(nup-1-shift, nup-1-shift-gap);
  }
  fSeeds->Sort();

  Int_t nseed=fSeeds->GetEntriesFast();
  for (Int_t i=0; i<nseed; i++) {
    //tracking in the outer sectors
    fSectors=fOuterSec; fN=fkNOS;

    AliTPCseed *pt=(AliTPCseed*)fSeeds->UncheckedAt(i), &t=*pt;
    if (!FollowProlongation(t)) {
       delete fSeeds->RemoveAt(i);
       continue;
    }

    //tracking in the inner sectors
    fSectors=fInnerSec; fN=fkNIS;

    Double_t alpha=t.GetAlpha() - fInnerSec->GetAlphaShift();
    if (alpha > 2.*TMath::Pi()) alpha -= 2.*TMath::Pi();
    if (alpha < 0.            ) alpha += 2.*TMath::Pi();
    Int_t ns=Int_t(alpha/fInnerSec->GetAlpha())%fkNIS;

    alpha=ns*fInnerSec->GetAlpha()+fInnerSec->GetAlphaShift()-t.GetAlpha();

    if (t.Rotate(alpha)) {
      if (FollowProlongation(t)) {
        if (t.GetNumberOfClusters() >= Int_t(0.4*nrows)) {
          t.CookdEdx();
          CookLabel(pt,0.1); //For comparison only
          pt->PropagateTo(fParam->GetInnerRadiusLow());
          AliESDtrack iotrack;
          iotrack.UpdateTrackParams(pt,AliESDtrack::kTPCin);

          event->AddTrack(&iotrack);

          UseClusters(&t);
        }
      }
    }
    delete fSeeds->RemoveAt(i);
  }

  Info("Clusters2Tracks","Number of found tracks : %d",
       event->GetNumberOfTracks());

  fSeeds->Clear(); delete fSeeds; fSeeds=0;

  return 0;
}

//_____________________________________________________________________________
Int_t AliTPCtracker::RefitInward(AliESDEvent* event) {
  //
  // The function propagates tracks throught TPC inward
  // using already associated clusters.
  // The clusters must be already loaded !
  //

  Int_t nTracks = event->GetNumberOfTracks();
  Int_t nRefited = 0;

  for (Int_t i = 0; i < nTracks; i++) {
    AliESDtrack* track = event->GetTrack(i);
    ULong_t status = track->GetStatus();

    if ( (status & AliESDtrack::kTPCrefit) != 0 ) continue;    
    if ( (status & AliESDtrack::kTPCout ) == 0 ) continue;

    if ( (status & AliESDtrack::kTRDout ) != 0 ) 
      if ( (status & AliESDtrack::kTRDrefit ) == 0 ) continue;

    AliTPCtrack* tpcTrack = new AliTPCtrack(*track);
    AliTPCseed* seed=new AliTPCseed(*tpcTrack); seed->ResetClusters(); 

    if ( (status & AliESDtrack::kTRDrefit) == 0 ) seed->ResetCovariance(10.);

    fSectors = fOuterSec;

    Int_t res = FollowRefitInward(seed, tpcTrack);
    UseClusters(seed);
    Int_t nc = seed->GetNumberOfClusters();

    fSectors = fInnerSec;

    res = FollowRefitInward(seed, tpcTrack);
    UseClusters(seed, nc);

    if (res) {
      seed->PropagateTo(fParam->GetInnerRadiusLow());
      seed->SetLabel(tpcTrack->GetLabel());
      seed->SetdEdx(tpcTrack->GetdEdx());
      track->UpdateTrackParams(seed, AliESDtrack::kTPCrefit);
      nRefited++;
    }

    delete seed;
    delete tpcTrack;
  }

  Info("RefitInward","Number of refitted tracks : %d",nRefited);

  return 0;
}

Int_t AliTPCtracker::PropagateBack(AliESDEvent *event) {
  //-----------------------------------------------------------------
  // This function propagates tracks back through the TPC.
  // The clusters must be already loaded !
  //-----------------------------------------------------------------
  Int_t nentr=event->GetNumberOfTracks();
  Info("PropagateBack", "Number of ESD tracks: %d\n", nentr);

  Int_t ntrk=0;
  for (Int_t i=0; i<nentr; i++) {
    AliESDtrack *esd=event->GetTrack(i);
    ULong_t status=esd->GetStatus();

    if ( (status & AliESDtrack::kTPCin ) == 0 ) continue;
    if ( (status & AliESDtrack::kTPCout) != 0 ) continue;
    if ( (status & AliESDtrack::kITSin) != 0 )
       if ( (status & AliESDtrack::kITSout) == 0 ) continue;

    const AliTPCtrack t(*esd);
    AliTPCseed s(t); s.ResetClusters();

    if ( (status & AliESDtrack::kITSout) == 0 ) s.ResetCovariance(10.);

    //inner sectors
    fSectors=fInnerSec; fN=fkNIS;

    Double_t alpha=s.GetAlpha() - fSectors->GetAlphaShift();
    if (alpha > 2.*TMath::Pi()) alpha -= 2.*TMath::Pi();
    if (alpha < 0.            ) alpha += 2.*TMath::Pi();
    Int_t ns=Int_t(alpha/fSectors->GetAlpha())%fN;
    alpha =ns*fSectors->GetAlpha() + fSectors->GetAlphaShift();
    alpha-=s.GetAlpha();

    if (!s.Rotate(alpha)) continue;
    if (!FollowBackProlongation(s,t)) continue;

    UseClusters(&s);

    //outer sectors
    fSectors=fOuterSec; fN=fkNOS;

    Int_t nc=s.GetNumberOfClusters();

    alpha=s.GetAlpha() - fSectors->GetAlphaShift();
    if (alpha > 2.*TMath::Pi()) alpha -= 2.*TMath::Pi();
    if (alpha < 0.            ) alpha += 2.*TMath::Pi();
    ns=Int_t(alpha/fSectors->GetAlpha())%fN;

    alpha =ns*fSectors->GetAlpha() + fSectors->GetAlphaShift();
    alpha-=s.GetAlpha();

    if (!s.Rotate(alpha)) continue;
    if (!FollowBackProlongation(s,t)) continue;
    {
    Int_t nrows=fOuterSec->GetNRows()+fInnerSec->GetNRows();
    if (s.GetNumberOfClusters() < Int_t(0.4*nrows)) continue;
    }
    s.PropagateTo(fParam->GetOuterRadiusUp());
    s.CookdEdx();
    CookLabel(&s,0.1); //For comparison only
    UseClusters(&s,nc);
    esd->UpdateTrackParams(&s,AliESDtrack::kTPCout);
    ntrk++;
  }
  Info("PropagateBack","Number of back propagated tracks: %d",ntrk);

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

  const AliTPCcluster *cl=0;

  if (sec<fkNIS) {
    cl=fInnerSec[sec][row].GetUnsortedCluster(ncl); 
  } else {
    sec-=fkNIS;
    cl=fOuterSec[sec][row].GetUnsortedCluster(ncl);
  }

  return (AliCluster*)cl;      
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
    if (j<noc) {
       lb[j]=lab;
       (mx[j])++;
    }
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
void AliTPCtracker::
AliTPCRow::InsertCluster(const AliTPCcluster* c, Int_t sec, Int_t row) {
  //-----------------------------------------------------------------------
  // Insert a cluster into this pad row in accordence with its y-coordinate
  //-----------------------------------------------------------------------
  if (fN==kMaxClusterPerRow) {
     ::Error("InsertCluster","Too many clusters !"); 
     return;
  }

  Int_t index=(((sec<<8)+row)<<16)+fN;

  if (fN==0) {
     fSize=kMaxClusterPerRow/8;
     fClusterArray=new AliTPCcluster[fSize];
     fIndex[0]=index;
     fClusterArray[0]=*c; fClusters[fN++]=fClusterArray; 
     return;
  }

  if (fN==fSize) {
     Int_t size=fSize*2;
     AliTPCcluster *buff=new AliTPCcluster[size];
     memcpy(buff,fClusterArray,fSize*sizeof(AliTPCcluster));
     for (Int_t i=0; i<fN; i++) 
        fClusters[i]=buff+(fClusters[i]-fClusterArray);
     delete[] fClusterArray;
     fClusterArray=buff;
     fSize=size;
  }

  Int_t i=Find(c->GetY());
  memmove(fClusters+i+1 ,fClusters+i,(fN-i)*sizeof(AliTPCcluster*));
  memmove(fIndex   +i+1 ,fIndex   +i,(fN-i)*sizeof(UInt_t));
  fIndex[i]=index; 
  fClusters[i]=fClusterArray+fN; fClusterArray[fN++]=*c;
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

  Int_t swap;//stupid sorting
  do {
    swap=0;
    for (i=0; i<nc-1; i++) {
      if (fdEdxSample[i]<=fdEdxSample[i+1]) continue;
      Float_t tmp=fdEdxSample[i];
      fdEdxSample[i]=fdEdxSample[i+1]; fdEdxSample[i+1]=tmp;
      swap++;
    }
  } while (swap);

  Int_t nl=Int_t(low*nc), nu=Int_t(up*nc);
  Float_t dedx=0;
  for (i=nl; i<=nu; i++) dedx += fdEdxSample[i];
  dedx /= (nu-nl+1);
  SetdEdx(dedx);

  return;
}



