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

/*
$Log$
Revision 1.5  2000/12/20 07:51:59  kowal2
Changes suggested by Alessandra and Paolo to avoid overlapped
data fields in encapsulated classes.

Revision 1.4  2000/11/02 07:27:16  kowal2
code corrections

Revision 1.2  2000/06/30 12:07:50  kowal2
Updated from the TPC-PreRelease branch

Revision 1.1.2.1  2000/06/25 08:53:55  kowal2
Splitted from AliTPCtracking

*/

//-------------------------------------------------------
//          Implementation of the TPC tracker
//
//   Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------

#include <TObjArray.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream.h>

#include "AliTPCtracker.h"
#include "AliTPCcluster.h"
#include "AliTPCClustersArray.h"
#include "AliTPCClustersRow.h"

const AliTPCParam *AliTPCtracker::AliTPCSector::fgParam;

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
inline Double_t f1(Double_t x1,Double_t y1,
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
inline Double_t f2(Double_t x1,Double_t y1,
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
inline Double_t f3(Double_t x1,Double_t y1, 
                   Double_t x2,Double_t y2,
                   Double_t z1,Double_t z2) 
{
  //-----------------------------------------------------------------
  // Initial approximation of the tangent of the track dip angle
  //-----------------------------------------------------------------
  return (z1 - z2)/sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

//_____________________________________________________________________________
Int_t AliTPCtracker::FindProlongation(AliTPCseed& t, const AliTPCSector *sec,
Int_t s, Int_t rf) 
{
  //-----------------------------------------------------------------
  // This function tries to find a track prolongation.
  //-----------------------------------------------------------------
  const Int_t kSKIP=(t.GetNumberOfClusters()<10) ? 10 : 
                    Int_t(0.5*sec->GetNRows());
  Int_t tryAgain=kSKIP;
  Double_t alpha=sec->GetAlpha();
  Int_t ns=Int_t(2*TMath::Pi()/alpha+0.5);

  for (Int_t nr=sec->GetRowNumber(t.GetX())-1; nr>=rf; nr--) {
    Double_t x=sec->GetX(nr), ymax=sec->GetMaxY(nr);
    if (!t.PropagateTo(x)) return 0;

    AliTPCcluster *cl=0;
    UInt_t index=0;
    Double_t maxchi2=12.;
    const AliTPCRow &krow=sec[s][nr];
    Double_t sy2=SigmaY2(t.GetX(),t.GetTgl(),1./t.Get1Pt());
    Double_t sz2=SigmaZ2(t.GetX(),t.GetTgl());
    Double_t road=4.*sqrt(t.GetSigmaY2() + sy2), y=t.GetY(), z=t.GetZ();

    if (road>30) {
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
      Float_t l=sec->GetPadPitchWidth();
      t.SetSampledEdx(cl->GetQ()/l,t.GetNumberOfClusters());
      if (!t.Update(cl,maxchi2,index)) {
         if (!tryAgain--) return 0;
      } else tryAgain=kSKIP;
    } else {
      if (tryAgain==0) break;
      if (y > ymax) {
         s = (s+1) % ns;
         if (!t.Rotate(alpha)) return 0;
      } else if (y <-ymax) {
         s = (s-1+ns) % ns;
         if (!t.Rotate(-alpha)) return 0;
      }
      tryAgain--;
    }
  }

  return 1;

}

//_____________________________________________________________________________
void 
AliTPCtracker::MakeSeeds(TObjArray& seeds,const AliTPCSector *sec, Int_t max,
Int_t i1, Int_t i2)
{
  //-----------------------------------------------------------------
  // This function creates track seeds.
  //-----------------------------------------------------------------
  Double_t x[5], c[15];

  Double_t alpha=sec->GetAlpha(), shift=sec->GetAlphaShift();
  Double_t cs=cos(alpha), sn=sin(alpha);

  Double_t x1 =sec->GetX(i1);
  Double_t xx2=sec->GetX(i2);

  for (Int_t ns=0; ns<max; ns++) {
    Int_t nl=sec[(ns-1+max)%max][i2];
    Int_t nm=sec[ns][i2];
    Int_t nu=sec[(ns+1)%max][i2];
    const AliTPCRow& kr1=sec[ns][i1];
    for (Int_t is=0; is < kr1; is++) {
      Double_t y1=kr1[is]->GetY(), z1=kr1[is]->GetZ();
      for (Int_t js=0; js < nl+nm+nu; js++) {
	const AliTPCcluster *kcl;
        Double_t x2,   y2,   z2;
        Double_t x3=0.,y3=0.;

	if (js<nl) {
	  const AliTPCRow& kr2=sec[(ns-1+max)%max][i2];
	  kcl=kr2[js];
          y2=kcl->GetY(); z2=kcl->GetZ();
          x2= xx2*cs+y2*sn;
          y2=-xx2*sn+y2*cs;
	} else 
	  if (js<nl+nm) {
	    const AliTPCRow& kr2=sec[ns][i2];
	    kcl=kr2[js-nl];
            x2=xx2; y2=kcl->GetY(); z2=kcl->GetZ();
	  } else {
	    const AliTPCRow& kr2=sec[(ns+1)%max][i2];
	    kcl=kr2[js-nl-nm];
            y2=kcl->GetY(); z2=kcl->GetZ();
            x2=xx2*cs-y2*sn;
            y2=xx2*sn+y2*cs;
	  }

        Double_t zz=z1 - z1/x1*(x1-x2); 
        if (TMath::Abs(zz-z2)>5.) continue;

        Double_t d=(x2-x1)*(0.-y2)-(0.-x2)*(y2-y1);
        if (d==0.) {cerr<<"MakeSeeds warning: Straight seed !\n"; continue;}

	x[0]=y1;
	x[1]=z1;
	x[3]=f1(x1,y1,x2,y2,x3,y3);
	if (TMath::Abs(x[3]) >= 0.0066) continue;
	x[2]=f2(x1,y1,x2,y2,x3,y3);
	//if (TMath::Abs(x[3]*x1-x[2]) >= 0.99999) continue;
	x[4]=f3(x1,y1,x2,y2,z1,z2);
	if (TMath::Abs(x[4]) > 1.2) continue;
	Double_t a=asin(x[2]);
	Double_t zv=z1 - x[4]/x[3]*(a+asin(x[3]*x1-x[2]));
	if (TMath::Abs(zv)>10.) continue; 

        Double_t sy1=kr1[is]->GetSigmaY2(), sz1=kr1[is]->GetSigmaZ2();
        Double_t sy2=kcl->GetSigmaY2(),     sz2=kcl->GetSigmaZ2();
	Double_t sy3=100*0.025, sy=0.1, sz=0.1;

	Double_t f30=(f1(x1,y1+sy,x2,y2,x3,y3)-x[3])/sy;
	Double_t f32=(f1(x1,y1,x2,y2+sy,x3,y3)-x[3])/sy;
	Double_t f34=(f1(x1,y1,x2,y2,x3,y3+sy)-x[3])/sy;
	Double_t f20=(f2(x1,y1+sy,x2,y2,x3,y3)-x[2])/sy;
	Double_t f22=(f2(x1,y1,x2,y2+sy,x3,y3)-x[2])/sy;
	Double_t f24=(f2(x1,y1,x2,y2,x3,y3+sy)-x[2])/sy;
	Double_t f40=(f3(x1,y1+sy,x2,y2,z1,z2)-x[4])/sy;
	Double_t f41=(f3(x1,y1,x2,y2,z1+sz,z2)-x[4])/sz;
	Double_t f42=(f3(x1,y1,x2,y2+sy,z1,z2)-x[4])/sy;
	Double_t f43=(f3(x1,y1,x2,y2,z1,z2+sz)-x[4])/sz;

        c[0]=sy1;
        c[1]=0.;       c[2]=sz1;
        c[3]=f20*sy1;  c[4]=0.;       c[5]=f20*sy1*f20+f22*sy2*f22+f24*sy3*f24;
        c[6]=f30*sy1;  c[7]=0.;       c[8]=f30*sy1*f20+f32*sy2*f22+f34*sy3*f24;
                                      c[9]=f30*sy1*f30+f32*sy2*f32+f34*sy3*f34;
        c[10]=f40*sy1; c[11]=f41*sz1; c[12]=f40*sy1*f20+f42*sy2*f22;
        c[13]=f40*sy1*f30+f42*sy2*f32;
        c[14]=f40*sy1*f40+f41*sz1*f41+f42*sy2*f42+f43*sz2*f43;

        UInt_t index=kr1.GetIndex(is);
	AliTPCseed *track=new AliTPCseed(index, x, c, x1, ns*alpha+shift);
        Float_t l=sec->GetPadPitchWidth();
        track->SetSampledEdx(kr1[is]->GetQ()/l,0);

        Int_t rc=FindProlongation(*track,sec,ns,i2);
        if (rc==0 || track->GetNumberOfClusters()<(i1-i2)/2) delete track;
        else seeds.AddLast(track); 
      }
    }
  }
}

//_____________________________________________________________________________
Int_t AliTPCtracker::Clusters2Tracks(const AliTPCParam *par, TFile *of) {
  //-----------------------------------------------------------------
  // This is a track finder.
  //-----------------------------------------------------------------
  TDirectory *savedir=gDirectory; 

  if (!of->IsOpen()) {
     cerr<<"AliTPCtracker::Clusters2Tracks(): output file not open !\n";
     return 1;
  }

  AliTPCClustersArray carray;
  carray.Setup(par);
  carray.SetClusterType("AliTPCcluster");
  carray.ConnectTree("Segment Tree");

  of->cd();
  TTree tracktree("TreeT","Tree with TPC tracks");
  AliTPCtrack *iotrack=0;
  tracktree.Branch("tracks","AliTPCtrack",&iotrack,32000,0);

  AliTPCSector::SetParam(par);

  const Int_t kNIS=par->GetNInnerSector()/2;
  AliTPCSSector *ssec=new AliTPCSSector[kNIS];         
  Int_t nlow=ssec->GetNRows();     

  const Int_t kNOS=par->GetNOuterSector()/2;
  AliTPCLSector *lsec=new AliTPCLSector[kNOS];
  Int_t nup=lsec->GetNRows();
    
  //Load outer sectors
  UInt_t index;
  Int_t i,j;

  j=Int_t(carray.GetTree()->GetEntries());
  for (i=0; i<j; i++) {
      AliSegmentID *s=carray.LoadEntry(i);
      Int_t sec,row;
      par->AdjustSectorRow(s->GetID(),sec,row);
      if (sec<kNIS*2) continue;
      AliTPCClustersRow *clrow=carray.GetRow(sec,row);
      Int_t ncl=clrow->GetArray()->GetEntriesFast();
      while (ncl--) {
         AliTPCcluster *c=(AliTPCcluster*)(*clrow)[ncl];
         index=(((sec<<8)+row)<<16)+ncl;
         lsec[(sec-kNIS*2)%kNOS][row].InsertCluster(c,index);
      }
  }

  //find track seeds
  TObjArray seeds(20000);
  Int_t nrows=nlow+nup;
  Int_t gap=Int_t(0.125*nrows), shift=Int_t(0.5*gap);
  MakeSeeds(seeds, lsec, kNOS, nup-1, nup-1-gap);
  MakeSeeds(seeds, lsec, kNOS, nup-1-shift, nup-1-shift-gap);    
  seeds.Sort();

  //tracking in outer sectors
  Int_t nseed=seeds.GetEntriesFast();
  for (i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)seeds.UncheckedAt(i), &t=*pt;
    Double_t alpha=t.GetAlpha();
    if (alpha > 2.*TMath::Pi()) alpha -= 2.*TMath::Pi();  
    if (alpha < 0.            ) alpha += 2.*TMath::Pi();  
    Int_t ns=Int_t(alpha/lsec->GetAlpha())%kNOS;

    if (FindProlongation(t,lsec,ns)) {
       t.UseClusters(&carray);
       continue;
    }
    delete seeds.RemoveAt(i);
  }  

  //unload outer sectors
  for (i=0; i<kNOS; i++) {
      for (j=0; j<nup; j++) {
          if (carray.GetRow(i+kNIS*2,j)) carray.ClearRow(i+kNIS*2,j);
          if (carray.GetRow(i+kNIS*2+kNOS,j)) carray.ClearRow(i+kNIS*2+kNOS,j);
      }
  }

  //load inner sectors
  j=Int_t(carray.GetTree()->GetEntries());
  for (i=0; i<j; i++) {
      AliSegmentID *s=carray.LoadEntry(i);
      Int_t sec,row;
      par->AdjustSectorRow(s->GetID(),sec,row);
      if (sec>=kNIS*2) continue;
      AliTPCClustersRow *clrow=carray.GetRow(sec,row);
      Int_t ncl=clrow->GetArray()->GetEntriesFast();
      while (ncl--) {
         AliTPCcluster *c=(AliTPCcluster*)(*clrow)[ncl];
         index=(((sec<<8)+row)<<16)+ncl;
         ssec[sec%kNIS][row].InsertCluster(c,index);
      }
  }

  //tracking in inner sectors
  Int_t found=0;
  for (i=0; i<nseed; i++) {
    AliTPCseed *pt=(AliTPCseed*)seeds.UncheckedAt(i), &t=*pt;
    if (!pt) continue;
    Int_t nc=t.GetNumberOfClusters();

    Double_t alpha=t.GetAlpha() - ssec->GetAlphaShift();
    if (alpha > 2.*TMath::Pi()) alpha -= 2.*TMath::Pi();
    if (alpha < 0.            ) alpha += 2.*TMath::Pi();
    Int_t ns=Int_t(alpha/ssec->GetAlpha())%kNIS;

    alpha=ns*ssec->GetAlpha() + ssec->GetAlphaShift() - t.GetAlpha();

    if (t.Rotate(alpha)) {
       if (FindProlongation(t,ssec,ns)) {
          if (t.GetNumberOfClusters() >= Int_t(0.4*nrows)) {
             t.CookdEdx();
             //t.CookLabel(&carray);
             iotrack=pt;
             tracktree.Fill();
             t.UseClusters(&carray,nc);
             cerr<<found++<<'\r';
          }
       }
    }
    delete seeds.RemoveAt(i); 
  }  
  tracktree.Write();

  //unload inner sectors
  for (i=0; i<kNIS; i++) {
      for (j=0; j<nlow; j++) {
          if (carray.GetRow(i,j)) carray.ClearRow(i,j);
          if (carray.GetRow(i+kNIS,j)) carray.ClearRow(i+kNIS,j);
      }
  }

  cerr<<"Number of found tracks : "<<found<<endl;

  delete[] ssec;
  delete[] lsec;

  savedir->cd();

  return 0;
}

//_________________________________________________________________________
void 
AliTPCtracker::AliTPCRow::InsertCluster(const AliTPCcluster* c, UInt_t index) {
  //-----------------------------------------------------------------------
  // Insert a cluster into this pad row in accordence with its y-coordinate
  //-----------------------------------------------------------------------
  if (fN==kMAXCLUSTER) {
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
}

//_____________________________________________________________________________
void AliTPCtracker::AliTPCseed::UseClusters(AliTPCClustersArray *ca, Int_t n) {
  //-----------------------------------------------------------------
  // This function marks clusters associated with this track.
  //-----------------------------------------------------------------
  Int_t nc=GetNumberOfClusters();
  Int_t sec,row,ncl;

  for (Int_t i=n; i<nc; i++) {
     GetCluster(i,sec,row,ncl);
     AliTPCClustersRow *clrow=ca->GetRow(sec,row);
     AliTPCcluster *c=(AliTPCcluster*)(*clrow)[ncl]; 
     c->Use();   
  }
}


