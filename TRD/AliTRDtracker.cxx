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
Revision 1.7  2000/12/08 16:07:02  cblume
Update of the tracking by Sergei

Revision 1.6  2000/11/30 17:38:08  cblume
Changes to get in line with new STEER and EVGEN

Revision 1.5  2000/11/14 14:40:27  cblume
Correction for the Sun compiler (kTRUE and kFALSE)

Revision 1.4  2000/11/10 14:57:52  cblume
Changes in the geometry constants for the DEC compiler

Revision 1.3  2000/10/15 23:40:01  cblume
Remove AliTRDconst

Revision 1.2  2000/10/06 16:49:46  cblume
Made Getters const

Revision 1.1.2.2  2000/10/04 16:34:58  cblume
Replace include files by forward declarations

Revision 1.1.2.1  2000/09/22 14:47:52  cblume
Add the tracking code

*/   

#include <iostream.h>

#include <TFile.h>
#include <TROOT.h>
#include <TBranch.h>
#include <TTree.h>

#include "AliRun.h"
#include "AliTRD.h"
#include "AliTRDgeometry.h"
#include "AliTRDrecPoint.h" 
#include "AliTRDcluster.h" 
#include "AliTRDtrack.h"
#include "AliTRDtrackingSector.h"
#include "AliTRDtimeBin.h"

#include "AliTRDtracker.h"

ClassImp(AliTRDtracker) 

  const  Int_t     AliTRDtracker::fSeedGap            = 35;  
  const  Int_t     AliTRDtracker::fSeedStep           = 5;   


  const  Float_t   AliTRDtracker::fMinClustersInTrack = 0.5;  
  const  Float_t   AliTRDtracker::fMinClustersInSeed  = 0.5;  
  const  Float_t   AliTRDtracker::fSeedDepth          = 0.5; 
  const  Float_t   AliTRDtracker::fSkipDepth          = 0.2;
  const  Float_t   AliTRDtracker::fMaxSeedDeltaZ      = 30.;  
  const  Float_t   AliTRDtracker::fMaxSeedC           = 0.01; 
  const  Float_t   AliTRDtracker::fMaxSeedTan         = 1.2;  
  const  Float_t   AliTRDtracker::fMaxSeedVertexZ     = 200.; 
  const  Float_t   AliTRDtracker::fLabelFraction      = 0.5;  
  const  Float_t   AliTRDtracker::fWideRoad           = 30.;

  const  Double_t  AliTRDtracker::fMaxChi2            = 12.; 
  const  Double_t  AliTRDtracker::fSeedErrorSY        = 0.1;
  const  Double_t  AliTRDtracker::fSeedErrorSY3       = 2.5;
  const  Double_t  AliTRDtracker::fSeedErrorSZ        = 0.1;

//____________________________________________________________________
AliTRDtracker::AliTRDtracker()
{
  //
  // Default constructor
  //   

  fEvent     = 0;
  fGeom      = NULL;

  fNclusters = 0;
  fClusters  = NULL; 
  fNseeds    = 0;
  fSeeds     = NULL;
  fNtracks   = 0;
  fTracks    = NULL;

}   

//____________________________________________________________________
AliTRDtracker::AliTRDtracker(const Text_t* name, const Text_t* title)
                  :TNamed(name, title)
{
  fEvent     = 0;
  fGeom      = NULL;

  fNclusters = 0;
  fClusters  = new TObjArray(2000); 
  fNseeds    = 0;
  fSeeds     = new TObjArray(20000);
  fNtracks   = 0;
  fTracks    = new TObjArray(10000);

}   

//___________________________________________________________________
AliTRDtracker::~AliTRDtracker()
{
  delete fClusters;
  delete fTracks;
  delete fSeeds;
  delete fGeom;
}   

//___________________________________________________________________
void AliTRDtracker::Clusters2Tracks()
{
  Int_t inner, outer;
  Int_t fTotalNofTimeBins = fGeom->GetTimeMax() * AliTRDgeometry::Nplan();
  Int_t nSteps = (Int_t) (fTotalNofTimeBins * fSeedDepth)/fSeedStep;

  for(Int_t i=0; i<nSteps; i++) {
    printf("step %d out of %d \n", i+1, nSteps);
    outer=fTotalNofTimeBins-1-i*fSeedStep; inner=outer-fSeedGap;
    MakeSeeds(inner,outer);
    FindTracks();
  } 
}          

//_____________________________________________________________________
Double_t AliTRDtracker::ExpectedSigmaY2(Double_t r, Double_t tgl, Double_t pt)
{
  // Parametrised "expected" error of the cluster reconstruction in Y 

  Double_t s = 0.2;    
  return s;
}

//_____________________________________________________________________
Double_t AliTRDtracker::ExpectedSigmaZ2(Double_t r, Double_t tgl)
{
  // Parametrised "expected" error of the cluster reconstruction in Z 

  Double_t s, pad = fGeom->GetRowPadSize();
  s = pad * pad /12.;  
  return s;
}                  

//_____________________________________________________________________
inline Double_t f1trd(Double_t x1,Double_t y1,
                      Double_t x2,Double_t y2,
                      Double_t x3,Double_t y3)
{
  // Initial approximation of the track curvature
  // Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch

  Double_t d=(x2-x1)*(y3-y2)-(x3-x2)*(y2-y1);
  Double_t a=0.5*((y3-y2)*(y2*y2-y1*y1+x2*x2-x1*x1)-
                  (y2-y1)*(y3*y3-y2*y2+x3*x3-x2*x2));
  Double_t b=0.5*((x2-x1)*(y3*y3-y2*y2+x3*x3-x2*x2)-
                  (x3-x2)*(y2*y2-y1*y1+x2*x2-x1*x1));

  Double_t xr=TMath::Abs(d/(d*x1-a)), yr=d/(d*y1-b);

  return -xr*yr/sqrt(xr*xr+yr*yr);
}          

//_____________________________________________________________________
inline Double_t f2trd(Double_t x1,Double_t y1,
                      Double_t x2,Double_t y2,
                      Double_t x3,Double_t y3)
{
  // Initial approximation of the track curvature times center of curvature
  // Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch

  Double_t d=(x2-x1)*(y3-y2)-(x3-x2)*(y2-y1);
  Double_t a=0.5*((y3-y2)*(y2*y2-y1*y1+x2*x2-x1*x1)-
                  (y2-y1)*(y3*y3-y2*y2+x3*x3-x2*x2));
  Double_t b=0.5*((x2-x1)*(y3*y3-y2*y2+x3*x3-x2*x2)-
                  (x3-x2)*(y2*y2-y1*y1+x2*x2-x1*x1));

  Double_t xr=TMath::Abs(d/(d*x1-a)), yr=d/(d*y1-b);

  return -a/(d*y1-b)*xr/sqrt(xr*xr+yr*yr);
}          

//_____________________________________________________________________
inline Double_t f3trd(Double_t x1,Double_t y1,
                      Double_t x2,Double_t y2,
                      Double_t z1,Double_t z2)
{
  // Initial approximation of the tangent of the track dip angle
  // Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch

  return (z1 - z2)/sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}            


//___________________________________________________________________

Int_t AliTRDtracker::FindProlongation(AliTRDtrack& t, AliTRDtrackingSector *sec,
                            Int_t s, Int_t rf)
{
  // Starting from current position on track=t this function tries 
  // to extrapolate the track up to timeBin=rf and to confirm prolongation
  // if a close cluster is found. *sec is a pointer to allocated
  // array of sectors, in which the initial sector has index=s. 

  const Int_t TIME_BINS_TO_SKIP=Int_t(fSkipDepth*sec->GetNtimeBins());
  Int_t try_again=TIME_BINS_TO_SKIP;

  Double_t alpha=AliTRDgeometry::GetAlpha();

  Int_t ns=Int_t(2*TMath::Pi()/alpha+0.5);

  for (Int_t nr=sec->GetTimeBinNumber(t.GetX())-1; nr>=rf; nr--) {
    Double_t x=sec->GetX(nr), ymax=sec->GetMaxY(nr);
    if (!t.PropagateTo(x)) {
      cerr<<"Can't propagate to x = "<<x<<endl;
      return 0;
    }

    AliTRDcluster *cl=0;
    UInt_t index=0;

    Double_t max_chi2=fMaxChi2;

    AliTRDtimeBin& time_bin=sec[s][nr];
    Double_t sy2=ExpectedSigmaY2(t.GetX(),t.GetTgl(),t.GetPt());
    Double_t sz2=ExpectedSigmaZ2(t.GetX(),t.GetTgl());
    Double_t road=5.*sqrt(t.GetSigmaY2() + sy2), y=t.GetY(), z=t.GetZ();

    if (road>fWideRoad) {
      if (t.GetNclusters() > 4) {
	cerr<<t.GetNclusters()<<" FindProlongation: Road is too wide !\n";
      }
      return 0;
    }       

    if (time_bin) {
      for (Int_t i=time_bin.Find(y-road); i<time_bin; i++) {
        AliTRDcluster* c=(AliTRDcluster*)(time_bin[i]);
        if (c->GetY() > y+road) break;
        if (c->IsUsed() > 0) continue;

	if((c->GetZ()-z)*(c->GetZ()-z) > 25. + sz2) continue;

        Double_t chi2=t.GetPredictedChi2(c);

        if (chi2 > max_chi2) continue;
        max_chi2=chi2;
        cl=c;
        index=time_bin.GetIndex(i);
      }   
    }
    if (cl) {

      //      Float_t l=sec->GetPitch();
      //      t.SetSampledEdx(cl->fQ/l,Int_t(t));  
     
      t.Update(cl,max_chi2,index);

      try_again=TIME_BINS_TO_SKIP;
    } else {
      if (try_again==0) break;
      if (y > ymax) {
	cerr<<"y > ymax: "<<y<<" > "<<ymax<<endl;
         s = (s+1) % ns;
         if (!t.Rotate(alpha)) {
	   cerr<<"Failed to rotate, alpha = "<<alpha<<endl;
	   return 0;
	 }
      } else if (y <-ymax) {
	cerr<<"y < -ymax: "<<y<<" < "<<-ymax<<endl;
         s = (s-1+ns) % ns;
         if (!t.Rotate(-alpha)) {
	   cerr<<"Failed to rotate, alpha = "<<alpha<<endl;
	   return 0;
	 }
      }
      try_again--;
    }
  }

  return 1;
}          



//_____________________________________________________________________________
void AliTRDtracker::GetEvent(const Char_t *hitfile, const Char_t *clusterfile)
{
  // Opens a ROOT-file with TRD-clusters and reads in the cluster-tree

  ReadClusters(fClusters, clusterfile);

  // get geometry from the file with hits

  TFile *fInputFile = (TFile*) gROOT->GetListOfFiles()->FindObject(hitfile);
  if (!fInputFile) {
    printf("AliTRDtracker::Open -- ");
    printf("Open the ALIROOT-file %s.\n",hitfile);
    fInputFile = new TFile(hitfile);
  }
  else {
    printf("AliTRDtracker::Open -- ");
    printf("%s is already open.\n",hitfile);
  }

  // Get AliRun object from file or create it if not on file

  gAlice = (AliRun*) fInputFile->Get("gAlice");
  if (gAlice) {
    printf("AliTRDtracker::GetEvent -- ");
    printf("AliRun object found on file.\n");
  }
  else {
    printf("AliTRDtracker::GetEvent -- ");
    printf("Could not find AliRun object.\n");
  }

  // Import the Trees for the event nEvent in the file
  Int_t nparticles = gAlice->GetEvent(fEvent);
  cerr<<"nparticles = "<<nparticles<<endl;

  if (nparticles <= 0) {
    printf("AliTRDtracker::GetEvent -- ");
    printf("No entries in the trees for event %d.\n",fEvent);
  }

  AliTRD *TRD = (AliTRD*) gAlice->GetDetector("TRD");
  fGeom = TRD->GetGeometry();

}     


//_____________________________________________________________________________
void AliTRDtracker::SetUpSectors(AliTRDtrackingSector *sec)
{
  // Fills clusters into TRD tracking_sectors 
  // Note that the numbering scheme for the TRD tracking_sectors 
  // differs from that of TRD sectors

  for (Int_t i=0; i<AliTRDgeometry::Nsect(); i++) sec[i].SetUp();

  //  Sort clusters into AliTRDtimeBin's within AliTRDtrackSector's 

  cerr<<"MakeSeeds: sorting clusters"<<endl;
              
  Int_t ncl=fClusters->GetEntriesFast();
  UInt_t index;
  while (ncl--) {
    printf("\r %d left",ncl); 
    AliTRDcluster *c=(AliTRDcluster*)fClusters->UncheckedAt(ncl);
    Int_t detector=c->GetDetector(), local_time_bin=c->GetLocalTimeBin();
    Int_t sector=fGeom->GetSector(detector);

    Int_t tracking_sector=sector;
    if(sector > 0) tracking_sector = AliTRDgeometry::kNsect - sector;

    Int_t tb=sec[sector].GetTimeBin(detector,local_time_bin); 
    index=ncl;
    sec[tracking_sector][tb].InsertCluster(c,index);
  }    
  printf("\r\n");
}


//_____________________________________________________________________________
void AliTRDtracker::MakeSeeds(Int_t inner, Int_t outer)
{
  // Creates track seeds using clusters in timeBins=i1,i2

  Int_t i2 = inner, i1 = outer; 

  if (!fClusters) return; 

  AliTRDtrackingSector fTrSec[AliTRDgeometry::kNsect];
  SetUpSectors(fTrSec);

  // find seeds

  Double_t x[5], c[15];
  Int_t max_sec=AliTRDgeometry::kNsect;

  Double_t alpha=AliTRDgeometry::GetAlpha();
  Double_t shift=AliTRDgeometry::GetAlpha()/2.;
  Double_t cs=cos(alpha), sn=sin(alpha);  

  Double_t x1 =fTrSec[0].GetX(i1);
  Double_t xx2=fTrSec[0].GetX(i2);  

  for (Int_t ns=0; ns<max_sec; ns++) {

    Int_t nl=fTrSec[(ns-1+max_sec)%max_sec][i2];
    Int_t nm=fTrSec[ns][i2];
    Int_t nu=fTrSec[(ns+1)%max_sec][i2];

    AliTRDtimeBin& r1=fTrSec[ns][i1];

    for (Int_t is=0; is < r1; is++) {
      Double_t y1=r1[is]->GetY(), z1=r1[is]->GetZ();

      for (Int_t js=0; js < nl+nm+nu; js++) {
	const AliTRDcluster *cl;
        Double_t x2,   y2,   z2;
        Double_t x3=0.,y3=0.;  

        if (js<nl) {
	  AliTRDtimeBin& r2=fTrSec[(ns-1+max_sec)%max_sec][i2];
	  cl=r2[js]; 
	  y2=cl->GetY(); z2=cl->GetZ();

          x2= xx2*cs+y2*sn;
          y2=-xx2*sn+y2*cs;          

        } else
          if (js<nl+nm) {
	    AliTRDtimeBin& r2=fTrSec[ns][i2];
	    cl=r2[js-nl];
	    x2=xx2; y2=cl->GetY(); z2=cl->GetZ();
          } else {
	    AliTRDtimeBin& r2=fTrSec[(ns+1)%max_sec][i2];
            cl=r2[js-nl-nm];
	    y2=cl->GetY(); z2=cl->GetZ();

            x2=xx2*cs-y2*sn;
            y2=xx2*sn+y2*cs;   

          }

        Double_t zz=z1 - z1/x1*(x1-x2);
	
        if (TMath::Abs(zz-z2)>fMaxSeedDeltaZ) continue;   

        Double_t d=(x2-x1)*(0.-y2)-(0.-x2)*(y2-y1);
        if (d==0.) {cerr<<"TRD MakeSeeds: Straight seed !\n"; continue;}

        x[0]=y1;
        x[1]=z1;
        x[2]=f1trd(x1,y1,x2,y2,x3,y3);

        if (TMath::Abs(x[2]) >= fMaxSeedC) continue;

        x[3]=f2trd(x1,y1,x2,y2,x3,y3);

        if (TMath::Abs(x[2]*x1-x[3]) >= 0.99999) continue;

        x[4]=f3trd(x1,y1,x2,y2,z1,z2);

        if (TMath::Abs(x[4]) > fMaxSeedTan) continue;

        Double_t a=asin(x[3]);
        Double_t zv=z1 - x[4]/x[2]*(a+asin(x[2]*x1-x[3]));
        if (TMath::Abs(zv)>fMaxSeedVertexZ) continue;    

        Double_t sy1=r1[is]->GetSigmaY2(), sz1=r1[is]->GetSigmaZ2();
        Double_t sy2=cl->GetSigmaY2(),     sz2=cl->GetSigmaZ2();
        Double_t sy3=fSeedErrorSY3, sy=fSeedErrorSY, sz=fSeedErrorSZ;

        Double_t f20=(f1trd(x1,y1+sy,x2,y2,x3,y3)-x[2])/sy;
        Double_t f22=(f1trd(x1,y1,x2,y2+sy,x3,y3)-x[2])/sy;
        Double_t f24=(f1trd(x1,y1,x2,y2,x3,y3+sy)-x[2])/sy;
        Double_t f30=(f2trd(x1,y1+sy,x2,y2,x3,y3)-x[3])/sy;
        Double_t f32=(f2trd(x1,y1,x2,y2+sy,x3,y3)-x[3])/sy;
        Double_t f34=(f2trd(x1,y1,x2,y2,x3,y3+sy)-x[3])/sy;
        Double_t f40=(f3trd(x1,y1+sy,x2,y2,z1,z2)-x[4])/sy;
        Double_t f41=(f3trd(x1,y1,x2,y2,z1+sz,z2)-x[4])/sz;
        Double_t f42=(f3trd(x1,y1,x2,y2+sy,z1,z2)-x[4])/sy;
        Double_t f43=(f3trd(x1,y1,x2,y2,z1,z2+sz)-x[4])/sz;

        c[0]=sy1;
        c[1]=0.;       c[2]=sz1;
        c[3]=f20*sy1;  c[4]=0.;       c[5]=f20*sy1*f20+f22*sy2*f22+f24*sy3*f24;
        c[6]=f30*sy1;  c[7]=0.;       c[8]=f30*sy1*f20+f32*sy2*f22+f34*sy3*f24;
                                      c[9]=f30*sy1*f30+f32*sy2*f32+f34*sy3*f34;
        c[10]=f40*sy1; c[11]=f41*sz1; c[12]=f40*sy1*f20+f42*sy2*f22;
        c[13]=f40*sy1*f30+f42*sy2*f32;
        c[14]=f40*sy1*f40+f41*sz1*f41+f42*sy2*f42+f43*sz2*f43;   

        UInt_t index=r1.GetIndex(is);
        AliTRDtrack *track=new AliTRDtrack(index, x, c, x1, ns*alpha+shift); 

        Int_t rc=FindProlongation(*track,fTrSec,ns,i2);
       
        if ((rc < 0) || 
            (track->GetNclusters() < (i1-i2)*fMinClustersInSeed)) delete track;
        else { 
	  fSeeds->AddLast(track); fNseeds++; 
	  cerr<<"found seed "<<fNseeds<<endl;
	}
      }
    }
  }

  fSeeds->Sort();
}          

//___________________________________________________________________
void AliTRDtracker::FindTracks() 
{
  if (!fClusters) return; 

  AliTRDtrackingSector fTrSec[AliTRDgeometry::kNsect];
  SetUpSectors(fTrSec);

  // find tracks

  Int_t num_of_time_bins = fTrSec[0].GetNtimeBins(); 
  Int_t nseed=fSeeds->GetEntriesFast();

  Int_t nSeedClusters;
  for (Int_t i=0; i<nseed; i++) {
    cerr<<"FindTracks: seed "<<i+1<<" out of "<<nseed<<endl;

    AliTRDtrack& t=*((AliTRDtrack*)fSeeds->UncheckedAt(i));

    nSeedClusters = t.GetNclusters();
    Double_t alpha=t.GetAlpha();

    if (alpha > 2.*TMath::Pi()) alpha -= 2.*TMath::Pi();
    if (alpha < 0.            ) alpha += 2.*TMath::Pi();
    Int_t ns=Int_t(alpha/AliTRDgeometry::GetAlpha())%AliTRDgeometry::kNsect;

    if (FindProlongation(t,fTrSec,ns)) {
      cerr<<"No of clusters in the track = "<<t.GetNclusters()<<endl; 
      if (t.GetNclusters() >= Int_t(fMinClustersInTrack*num_of_time_bins)) {
	Int_t label = GetTrackLabel(t);
	t.SetLabel(label);
	UseClusters(t);

        AliTRDtrack *pt = new AliTRDtrack(t);
        fTracks->AddLast(pt); fNtracks++;     

	cerr<<"found track "<<fNtracks<<endl;
      }                         
    }     
    delete fSeeds->RemoveAt(i);  
  }            
}

//__________________________________________________________________
void AliTRDtracker::UseClusters(AliTRDtrack t) {
  Int_t ncl=t.GetNclusters();
  for (Int_t i=0; i<ncl; i++) {
    Int_t index = t.GetClusterIndex(i);
    AliTRDcluster *c=(AliTRDcluster*)fClusters->UncheckedAt(index);
    c->Use();
  }
}

//__________________________________________________________________
Int_t AliTRDtracker::GetTrackLabel(AliTRDtrack t) {

  Int_t label=123456789, index, i, j;
  Int_t ncl=t.GetNclusters();
  const Int_t range = AliTRDgeometry::kNplan * fGeom->GetTimeMax();
  Bool_t label_added;

  //  Int_t s[range][2];
  Int_t **s = new Int_t* [range];
  for (i=0; i<range; i++) {
    s[i] = new Int_t[2];
  }
  for (i=0; i<range; i++) {
    s[i][0]=-1;
    s[i][1]=0;
  }

  Int_t t0,t1,t2;
  for (i=0; i<ncl; i++) {
    index=t.GetClusterIndex(i);
    AliTRDcluster *c=(AliTRDcluster*)fClusters->UncheckedAt(index);
    t0=c->GetTrackIndex(0);
    t1=c->GetTrackIndex(1);
    t2=c->GetTrackIndex(2);
  }

  for (i=0; i<ncl; i++) {
    index=t.GetClusterIndex(i);
    AliTRDcluster *c=(AliTRDcluster*)fClusters->UncheckedAt(index);
    for (Int_t k=0; k<3; k++) { 
      label=c->GetTrackIndex(k);
      label_added=kFALSE; j=0;
      if (label >= 0) {
	while ( (!label_added) && ( j < range ) ) {
	  if (s[j][0]==label || s[j][1]==0) {
	    s[j][0]=label; 
	    s[j][1]=s[j][1]+1; 
	    label_added=kTRUE;
	  }
	  j++;
	}
      }
    }
  }

  Int_t max=0;
  label = -123456789;

  for (i=0; i<range; i++) {
    if (s[i][1]>max) {
      max=s[i][1]; label=s[i][0];
    }
  }
  delete []s;
  if(max > ncl*fLabelFraction) return label;
  else return -1;
}

//___________________________________________________________________

Int_t AliTRDtracker::WriteTracks(const Char_t *filename) {

  TDirectory *savedir=gDirectory;   

  TFile *out=TFile::Open(filename,"RECREATE");

  TTree tracktree("TreeT","Tree with TRD tracks");

  AliTRDtrack *iotrack=0;
  tracktree.Branch("tracks","AliTRDtrack",&iotrack,32000,0);  

  Int_t ntracks=fTracks->GetEntriesFast();

  for (Int_t i=0; i<ntracks; i++) {
    AliTRDtrack *pt=(AliTRDtrack*)fTracks->UncheckedAt(i);
    iotrack=pt;
    tracktree.Fill(); 
    cerr<<"WriteTracks: put track "<<i<<" in the tree"<<endl;
  }

  tracktree.Write(); 
  out->Close(); 

  savedir->cd();  

  cerr<<"WriteTracks: done"<<endl;
  return 0;
}

//_____________________________________________________________________________
void AliTRDtracker::ReadClusters(TObjArray *array, const Char_t *filename, 
Int_t option) 
{
  //
  // Reads AliTRDclusters (option >= 0) or AliTRDrecPoints (option < 0) 
  // from the file. The names of the cluster tree and branches 
  // should match the ones used in AliTRDclusterizer::WriteClusters()
  //

  TDirectory *savedir=gDirectory; 

  TFile *file = TFile::Open(filename);
  if (!file->IsOpen()) {printf("Can't open file %s !\n",filename); return;} 

  TTree *tree = (TTree*)file->Get("ClusterTree");
  Int_t nentr = (Int_t) tree->GetEntries();
  printf("found %d entries in %s.\n",nentr,tree->GetName());

  TBranch *branch;
  TObjArray *ioArray = new TObjArray(400);

  if( option < 0 ) {
    branch = tree->GetBranch("RecPoints");

    for (Int_t i=0; i<nentr; i++) {
      branch->SetAddress(&ioArray);
      tree->GetEvent(i);
      Int_t npoints = ioArray->GetEntriesFast();
      printf("Read %d rec. points from entry %d \n", npoints, i);

      for(Int_t j=0; j<npoints; j++) {
	AliTRDrecPoint *p=(AliTRDrecPoint*)ioArray->UncheckedAt(j);
	array->AddLast(p);
	ioArray->RemoveAt(j);
      }
    }
  }
  else {
    branch = tree->GetBranch("Clusters");

    for (Int_t i=0; i<nentr; i++) {
      branch->SetAddress(&ioArray);
      tree->GetEvent(i);
      Int_t npoints = ioArray->GetEntriesFast();
      printf("Read %d clusters from entry %d \n", npoints, i);

      for(Int_t j=0; j<npoints; j++) {
	AliTRDcluster *c=(AliTRDcluster*)ioArray->UncheckedAt(j);
	array->AddLast(c);
	ioArray->RemoveAt(j);
      }
    }
  }

  file->Close();                   
  savedir->cd(); 
  
}

