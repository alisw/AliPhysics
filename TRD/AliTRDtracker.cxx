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
Revision 1.14  2001/11/14 10:50:46  cblume
Changes in digits IO. Add merging of summable digits

Revision 1.13  2001/05/30 12:17:47  hristov
Loop variables declared once

Revision 1.12  2001/05/28 17:07:58  hristov
Last minute changes; ExB correction in AliTRDclusterizerV1; taking into account of material in G10 TEC frames and material between TEC planes (C.Blume,S.Sedykh)

Revision 1.8  2000/12/20 13:00:44  cblume
Modifications for the HP-compiler

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

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  The TRD tracker                                                           //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

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

  const  Float_t     AliTRDtracker::fgkSeedDepth          = 0.5; 
  const  Float_t     AliTRDtracker::fgkSeedStep           = 0.05;   
  const  Float_t     AliTRDtracker::fgkSeedGap            = 0.25;  

  const  Float_t     AliTRDtracker::fgkMaxSeedDeltaZ12    = 40.;  
  const  Float_t     AliTRDtracker::fgkMaxSeedDeltaZ      = 25.;  
  const  Float_t     AliTRDtracker::fgkMaxSeedC           = 0.0052; 
  const  Float_t     AliTRDtracker::fgkMaxSeedTan         = 1.2;  
  const  Float_t     AliTRDtracker::fgkMaxSeedVertexZ     = 150.; 

  const  Double_t    AliTRDtracker::fgkSeedErrorSY        = 0.2;
  const  Double_t    AliTRDtracker::fgkSeedErrorSY3       = 2.5;
  const  Double_t    AliTRDtracker::fgkSeedErrorSZ        = 0.1;

  const  Float_t     AliTRDtracker::fgkMinClustersInSeed  = 0.7;  

  const  Float_t     AliTRDtracker::fgkMinClustersInTrack = 0.5;  
  const  Float_t     AliTRDtracker::fgkSkipDepth          = 0.05;
  const  Float_t     AliTRDtracker::fgkLabelFraction      = 0.5;  
  const  Float_t     AliTRDtracker::fgkWideRoad           = 20.;

  const  Double_t    AliTRDtracker::fgkMaxChi2            = 24.; 

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

  fSY2corr = 0.025;
}   

//____________________________________________________________________
AliTRDtracker::AliTRDtracker(const Text_t* name, const Text_t* title)
                  :TNamed(name, title)
{
  //
  // TRD tracker contructor
  //

  fEvent     = 0;
  fGeom      = NULL;

  fNclusters = 0;
  fClusters  = new TObjArray(2000); 
  fNseeds    = 0;
  fSeeds     = new TObjArray(20000);
  fNtracks   = 0;
  fTracks    = new TObjArray(10000);

  fSY2corr = 0.025;
}   

//___________________________________________________________________
AliTRDtracker::~AliTRDtracker()
{
  //
  // Destructor
  //

  delete fClusters;
  delete fTracks;
  delete fSeeds;
  delete fGeom;

}   

//___________________________________________________________________
void AliTRDtracker::Clusters2Tracks(TH1F *hs, TH1F *hd)
{
  //
  // Do the trackfinding
  //

  Int_t i;

  Int_t inner, outer;
  Int_t fTotalNofTimeBins = fGeom->GetTimeMax() * AliTRDgeometry::Nplan();
  Int_t nSteps = (Int_t) (fgkSeedDepth / fgkSeedStep);
  Int_t gap = (Int_t) (fTotalNofTimeBins * fgkSeedGap);
  Int_t step = (Int_t) (fTotalNofTimeBins * fgkSeedStep);


  //  nSteps = 1;

  if (!fClusters) return; 

  AliTRDtrackingSector fTrSec[AliTRDgeometry::kNsect];
  SetUpSectors(fTrSec);

  // make a first turn looking for seed ends in the same (n,n) 
  // and in the adjacent sectors (n,n+1)

  for(i=0; i<nSteps; i++) {
    printf("step %d out of %d \n", i+1, nSteps);
    outer=fTotalNofTimeBins-1-i*step; inner=outer-gap;
    MakeSeeds(inner,outer, fTrSec, 1, hs, hd);
    FindTracks(fTrSec, hs, hd);
  } 

  // make a second turn looking for seed ends in next-to-adjacent 
  // sectors (n,n+2)

  for(i=0; i<nSteps; i++) {
    printf("step %d out of %d \n", i+1, nSteps);
    outer=fTotalNofTimeBins-1-i*step; inner=outer-gap;
    MakeSeeds(inner, outer, fTrSec, 2, hs, hd);
    FindTracks(fTrSec,hs,hd);
  } 

}          

//_____________________________________________________________________
Double_t AliTRDtracker::ExpectedSigmaY2(Double_t r, Double_t tgl, Double_t pt) const
{
  //
  // Parametrised "expected" error of the cluster reconstruction in Y 
  //

  Double_t s = 0.08 * 0.08;    
  return s;

}

//_____________________________________________________________________
Double_t AliTRDtracker::ExpectedSigmaZ2(Double_t r, Double_t tgl) const
{
  //
  // Parametrised "expected" error of the cluster reconstruction in Z 
  //

  Double_t s = 6 * 6 /12.;  
  return s;

}                  

//_____________________________________________________________________
Double_t f1trd(Double_t x1,Double_t y1,
	       Double_t x2,Double_t y2,
	       Double_t x3,Double_t y3)
{
  //
  // Initial approximation of the track curvature
  // Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
  //

  Double_t d=(x2-x1)*(y3-y2)-(x3-x2)*(y2-y1);
  Double_t a=0.5*((y3-y2)*(y2*y2-y1*y1+x2*x2-x1*x1)-
                  (y2-y1)*(y3*y3-y2*y2+x3*x3-x2*x2));
  Double_t b=0.5*((x2-x1)*(y3*y3-y2*y2+x3*x3-x2*x2)-
                  (x3-x2)*(y2*y2-y1*y1+x2*x2-x1*x1));

  Double_t xr=TMath::Abs(d/(d*x1-a)), yr=d/(d*y1-b);

  return -xr*yr/sqrt(xr*xr+yr*yr);

}          

//_____________________________________________________________________
Double_t f2trd(Double_t x1,Double_t y1,
	       Double_t x2,Double_t y2,
	       Double_t x3,Double_t y3)
{
  //
  // Initial approximation of the track curvature times center of curvature
  // Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
  //

  Double_t d=(x2-x1)*(y3-y2)-(x3-x2)*(y2-y1);
  Double_t a=0.5*((y3-y2)*(y2*y2-y1*y1+x2*x2-x1*x1)-
                  (y2-y1)*(y3*y3-y2*y2+x3*x3-x2*x2));
  Double_t b=0.5*((x2-x1)*(y3*y3-y2*y2+x3*x3-x2*x2)-
                  (x3-x2)*(y2*y2-y1*y1+x2*x2-x1*x1));

  Double_t xr=TMath::Abs(d/(d*x1-a)), yr=d/(d*y1-b);

  return -a/(d*y1-b)*xr/sqrt(xr*xr+yr*yr);

}          

//_____________________________________________________________________
Double_t f3trd(Double_t x1,Double_t y1,
	       Double_t x2,Double_t y2,
	       Double_t z1,Double_t z2)
{
  //
  // Initial approximation of the tangent of the track dip angle
  // Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch
  //

  return (z1 - z2)/sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));

}            


//___________________________________________________________________

Int_t AliTRDtracker::FindProlongation(AliTRDtrack& t, AliTRDtrackingSector *sec,
                            Int_t s, Int_t rf, Int_t matchedIndex, 
				      TH1F *hs, TH1F *hd)
{
  //
  // Starting from current position on track=t this function tries 
  // to extrapolate the track up to timeBin=rf and to confirm prolongation
  // if a close cluster is found. *sec is a pointer to allocated
  // array of sectors, in which the initial sector has index=s. 
  //

  //  TH1F *hsame = hs;     
  //  TH1F *hdiff = hd;   

  //  Bool_t good_match;

  const Int_t kTimeBinsToSkip=Int_t(fgkSkipDepth*sec->GetNtimeBins());
  Int_t tryAgain=kTimeBinsToSkip;

  Double_t alpha=AliTRDgeometry::GetAlpha();

  Int_t ns=Int_t(2*TMath::Pi()/alpha+0.5);

  Double_t x0, rho;

  for (Int_t nr=sec->GetTimeBinNumber(t.GetX())-1; nr>=rf; nr--) {

    Double_t x=sec->GetX(nr);
    Double_t ymax=x*TMath::Tan(0.5*alpha);

    rho = 0.00295; x0 = 11.09;  // TEC
    if(sec->TECframe(nr,t.GetY(),t.GetZ())) { 
      rho = 1.7; x0 = 33.0;     // G10 frame 
    } 
    if(TMath::Abs(x - t.GetX()) > 3) { 
      rho = 0.0559; x0 = 55.6;  // radiator
    }
    if (!t.PropagateTo(x,x0,rho,0.139)) {
      cerr<<"Can't propagate to x = "<<x<<endl;
      return 0;
    }

    AliTRDtimeBin& timeBin=sec[s][nr];
    Double_t sy2=ExpectedSigmaY2(t.GetX(),t.GetTgl(),t.GetPt());
    Double_t sz2=ExpectedSigmaZ2(t.GetX(),t.GetTgl());
    Double_t road=25.*sqrt(t.GetSigmaY2() + sy2), y=t.GetY(), z=t.GetZ();

    if (road>fgkWideRoad) {
      if (t.GetNclusters() > 4) {
	cerr<<t.GetNclusters()<<" FindProlongation: Road is too wide !\n";
      }
      return 0;
    }       

    AliTRDcluster *cl=0;
    UInt_t index=0;
    //    Int_t ncl = 0;

    Double_t maxChi2=fgkMaxChi2;

    if (timeBin) {

      for (Int_t i=timeBin.Find(y-road); i<timeBin; i++) {
        AliTRDcluster* c=(AliTRDcluster*)(timeBin[i]);

	//	good_match = kFALSE;
	//	if((c->GetTrackIndex(0) == matchedIndex) ||
	//   (c->GetTrackIndex(1) == matchedIndex) ||
	//   (c->GetTrackIndex(2) == matchedIndex)) good_match = kTRUE;
	//	  if(hsame) hsame->Fill(TMath::Abs(c->GetY()-y)/road);
	//	  if(hdiff) hdiff->Fill(road);

        if (c->GetY() > y+road) break;
        if (c->IsUsed() > 0) continue;

	//	if(good_match) hsame->Fill(TMath::Abs(c->GetZ()-z));
	//	else hdiff->Fill(TMath::Abs(c->GetZ()-z));

	//	if(!good_match) continue;

	if((c->GetZ()-z)*(c->GetZ()-z) > 3 * 12 * sz2) continue;

        Double_t chi2=t.GetPredictedChi2(c);

	//	if((c->GetTrackIndex(0) == matchedIndex) ||
	//	   (c->GetTrackIndex(1) == matchedIndex) ||
	//	   (c->GetTrackIndex(2) == matchedIndex))
	//	  hdiff->Fill(chi2);

	//	ncl++;

        if (chi2 > maxChi2) continue;
        maxChi2=chi2;
        cl=c;
        index=timeBin.GetIndex(i);
      }   
      
      if(!cl) {

	for (Int_t i=timeBin.Find(y-road); i<timeBin; i++) {
	  AliTRDcluster* c=(AliTRDcluster*)(timeBin[i]);

	  if (c->GetY() > y+road) break;
	  if (c->IsUsed() > 0) continue;	  
	  if((c->GetZ()-z)*(c->GetZ()-z) > 3.25 * 12 * sz2) continue;
	  
	  Double_t chi2=t.GetPredictedChi2(c);
	  
	  //	  ncl++;

	  if (chi2 > maxChi2) continue;
	  maxChi2=chi2;
	  cl=c;
	  index=timeBin.GetIndex(i);
	}   
      }
      
    }

    if (cl) {

      t.Update(cl,maxChi2,index);

      tryAgain=kTimeBinsToSkip;
    } else {
      if (tryAgain==0) break;
      if (y > ymax) {
	//	cerr<<"y > ymax: "<<y<<" > "<<ymax<<endl;
         s = (s+1) % ns;
         if (!t.Rotate(alpha)) {
	   cerr<<"Failed to rotate, alpha = "<<alpha<<endl;
	   return 0;
	 }
      } else if (y <-ymax) {
	//	cerr<<"y < -ymax: "<<y<<" < "<<-ymax<<endl;
         s = (s-1+ns) % ns;
         if (!t.Rotate(-alpha)) {
	   cerr<<"Failed to rotate, alpha = "<<alpha<<endl;
	   return 0;
	 }
      }
      if(!sec->TECframe(nr,y,z)) tryAgain--;
    }
  }

  return 1;
}          



//_____________________________________________________________________________
void AliTRDtracker::GetEvent(const Char_t *hitfile, const Char_t *clusterfile)
{
  //
  // Opens a ROOT-file with TRD-clusters and reads the cluster-tree in
  //

  ReadClusters(fClusters, clusterfile);

  // get geometry from the file with hits

  TFile *fInputFile = (TFile*) gROOT->GetListOfFiles()->FindObject(hitfile);
  if (!fInputFile) {
    printf("AliTRDtracker::Open -- ");
    printf("Open the ALIROOT-file %s.\n",hitfile);
    fInputFile = new TFile(hitfile,"UPDATE");
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

  /*  
  // Import the Trees for the event nEvent in the file
  Int_t nparticles = gAlice->GetEvent(fEvent);
  cerr<<"nparticles = "<<nparticles<<endl;

  if (nparticles <= 0) {
    printf("AliTRDtracker::GetEvent -- ");
    printf("No entries in the trees for event %d.\n",fEvent);
  }
  */  

  AliTRD *trd = (AliTRD*) gAlice->GetDetector("TRD");
  fGeom = trd->GetGeometry();

}     


//_____________________________________________________________________________
void AliTRDtracker::SetUpSectors(AliTRDtrackingSector *sec)
{
  //
  // Fills clusters into TRD tracking_sectors 
  // Note that the numbering scheme for the TRD tracking_sectors 
  // differs from that of TRD sectors
  //

  for (Int_t i=0; i<AliTRDgeometry::Nsect(); i++) sec[i].SetUp();

  //  Sort clusters into AliTRDtimeBin's within AliTRDtrackSector's 

  cerr<<"SetUpSectors: sorting clusters"<<endl;
              
  Int_t ncl=fClusters->GetEntriesFast();
  UInt_t index;
  while (ncl--) {
    printf("\r %d left  ",ncl); 
    AliTRDcluster *c=(AliTRDcluster*)fClusters->UncheckedAt(ncl);
    Int_t detector=c->GetDetector(), localTimeBin=c->GetLocalTimeBin();
    Int_t sector=fGeom->GetSector(detector);

    Int_t trackingSector = AliTRDgeometry::kNsect - sector - 1;

    Int_t tb=sec[sector].GetTimeBin(detector,localTimeBin); 
    index=ncl;
    sec[trackingSector][tb].InsertCluster(c,index);

  }    
  printf("\r\n");
}


//_____________________________________________________________________________
void AliTRDtracker::MakeSeeds(Int_t inner, Int_t outer, 
                              AliTRDtrackingSector* fTrSec, Int_t turn,
			      TH1F *hs, TH1F *hd)
{
  //
  // Creates track seeds using clusters in timeBins=i1,i2
  //

  Int_t i2 = inner, i1 = outer; 
  Int_t ti[3], to[3];
  Int_t nprim = 85600/2;


  TH1F *hsame = hs;
  TH1F *hdiff = hd;   
  Bool_t match = false;
  Int_t matchedIndex;

  // find seeds

  Double_t x[5], c[15];
  Int_t maxSec=AliTRDgeometry::kNsect;

  Double_t alpha=AliTRDgeometry::GetAlpha();
  Double_t shift=AliTRDgeometry::GetAlpha()/2.;
  Double_t cs=cos(alpha), sn=sin(alpha);  
  Double_t cs2=cos(2.*alpha), sn2=sin(2.*alpha);  

  Double_t x1 =fTrSec[0].GetX(i1);
  Double_t xx2=fTrSec[0].GetX(i2);  


  printf("\n");

  if((turn != 1)&&(turn != 2)) {
    printf("*** Error in MakeSeeds: unexpected turn = %d  \n", turn);
    return;
  }


  for (Int_t ns=0; ns<maxSec; ns++) {

    printf("\n MakeSeeds: sector %d \n", ns); 

    Int_t nl2=fTrSec[(ns-2+maxSec)%maxSec][i2]; 
    Int_t nl=fTrSec[(ns-1+maxSec)%maxSec][i2];
    Int_t nm=fTrSec[ns][i2];
    Int_t nu=fTrSec[(ns+1)%maxSec][i2];
    Int_t nu2=fTrSec[(ns+2)%maxSec][i2]; 

    AliTRDtimeBin& r1=fTrSec[ns][i1];

    for (Int_t is=0; is < r1; is++) {
      Double_t y1=r1[is]->GetY(), z1=r1[is]->GetZ();
      for(Int_t ii=0; ii<3; ii++) to[ii] = r1[is]->GetTrackIndex(ii); 

      for (Int_t js=0; js < nl2+nl+nm+nu+nu2; js++) {
	  
	const AliTRDcluster *cl;
	Double_t x2,   y2,   z2;
	Double_t x3=0., y3=0.;  

	if (js<nl2) {
	  if(turn != 2) continue;
	  AliTRDtimeBin& r2=fTrSec[(ns-2+maxSec)%maxSec][i2];
	  cl=r2[js];
	  y2=cl->GetY(); z2=cl->GetZ();
	  for(Int_t ii=0; ii<3; ii++) ti[ii] = cl->GetTrackIndex(ii);

	  x2= xx2*cs2+y2*sn2;
	  y2=-xx2*sn2+y2*cs2;
	}        
	else if (js<nl2+nl) {
	  if(turn != 1) continue;
	  AliTRDtimeBin& r2=fTrSec[(ns-1+maxSec)%maxSec][i2];
	  cl=r2[js-nl2];
	  y2=cl->GetY(); z2=cl->GetZ();
	  for(Int_t ii=0; ii<3; ii++) ti[ii] = cl->GetTrackIndex(ii);

	  x2= xx2*cs+y2*sn;
	  y2=-xx2*sn+y2*cs;

	}
	else if (js<nl2+nl+nm) {
	  if(turn != 1) continue;
	  AliTRDtimeBin& r2=fTrSec[ns][i2];
	  cl=r2[js-nl2-nl];
	  x2=xx2; y2=cl->GetY(); z2=cl->GetZ();
	  for(Int_t ii=0; ii<3; ii++) ti[ii] = cl->GetTrackIndex(ii);
	}
	else if (js<nl2+nl+nm+nu) {
	  if(turn != 1) continue;
	  AliTRDtimeBin& r2=fTrSec[(ns+1)%maxSec][i2];
	  cl=r2[js-nl2-nl-nm];
	  y2=cl->GetY(); z2=cl->GetZ();
	  for(Int_t ii=0; ii<3; ii++) ti[ii] = cl->GetTrackIndex(ii);

	  x2=xx2*cs-y2*sn;
	  y2=xx2*sn+y2*cs;

	}
	else {
	  if(turn != 2) continue;
	  AliTRDtimeBin& r2=fTrSec[(ns+2)%maxSec][i2];
	  cl=r2[js-nl2-nl-nm-nu];
	  y2=cl->GetY(); z2=cl->GetZ();
	  for(Int_t ii=0; ii<3; ii++) ti[ii] = cl->GetTrackIndex(ii);
	  
	  x2=xx2*cs2-y2*sn2;
	  y2=xx2*sn2+y2*cs2;
	}         
	

	match = false;
	matchedIndex = -1;
	for (Int_t ii=0; ii<3; ii++) {
	  // cerr<<"ti["<<ii<<"] = "<<ti[ii]<<"; to["<<ii<<"] = "<<to[ii]<<endl;
	  if(ti[ii] < 0) continue;
	  if(ti[ii] >= nprim) continue;
	  for (Int_t kk=0; kk<3; kk++) {
	    if(to[kk] < 0) continue;
	    if(to[kk] >= nprim) continue;
	    if(ti[ii] == to[kk]) {
	      //cerr<<"ti["<<ii<<"] = "<<ti[ii]<<" = "<<to[kk]<<" = to["<<kk<<"]"<<endl;
	      matchedIndex = ti[ii];
	      match = true;
	    }
	  }
	}                 
	
	if(TMath::Abs(z1-z2) > fgkMaxSeedDeltaZ12) continue;

        Double_t zz=z1 - z1/x1*(x1-x2);

        if (TMath::Abs(zz-z2)>fgkMaxSeedDeltaZ) continue;   

        Double_t d=(x2-x1)*(0.-y2)-(0.-x2)*(y2-y1);
        if (d==0.) {cerr<<"TRD MakeSeeds: Straight seed !\n"; continue;}

        x[0]=y1;
        x[1]=z1;
        x[2]=f1trd(x1,y1,x2,y2,x3,y3);

        if (TMath::Abs(x[2]) > fgkMaxSeedC) continue;

        x[3]=f2trd(x1,y1,x2,y2,x3,y3);

        if (TMath::Abs(x[2]*x1-x[3]) >= 0.99999) continue;

        x[4]=f3trd(x1,y1,x2,y2,z1,z2);

        if (TMath::Abs(x[4]) > fgkMaxSeedTan) continue;

        Double_t a=asin(x[3]);
        Double_t zv=z1 - x[4]/x[2]*(a+asin(x[2]*x1-x[3]));

        if (TMath::Abs(zv)>fgkMaxSeedVertexZ) continue;    

        Double_t sy1=r1[is]->GetSigmaY2(), sz1=r1[is]->GetSigmaZ2();
        Double_t sy2=cl->GetSigmaY2(),     sz2=cl->GetSigmaZ2();
        Double_t sy3=fgkSeedErrorSY3, sy=fgkSeedErrorSY, sz=fgkSeedErrorSZ;

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
	
        AliTRDtrack *track=new AliTRDtrack(r1[is],index,x,c,x1,ns*alpha+shift); 

        Int_t rc=FindProlongation(*track,fTrSec,ns,i2,matchedIndex,hsame,hdiff);

	//	if (match) hsame->Fill((Float_t) track->GetNclusters());
	//	else hdiff->Fill((Float_t) track->GetNclusters());  
	//	delete track;
	//	continue;

        if ((rc < 0) || 
            (track->GetNclusters() < (i1-i2)*fgkMinClustersInSeed)) delete track;
        else { 
	  fSeeds->AddLast(track); fNseeds++; 
	  printf("\r found seed %d  ", fNseeds);
	}
      }
    }
  }

  fSeeds->Sort();
}          

//_____________________________________________________________________________
void AliTRDtracker::ReadClusters(TObjArray *array, const Char_t *filename) 
{
  //
  // Reads AliTRDclusters (option >= 0) or AliTRDrecPoints (option < 0) 
  // from the file. The names of the cluster tree and branches 
  // should match the ones used in AliTRDclusterizer::WriteClusters()
  //

  TDirectory *savedir=gDirectory; 

  TFile *file = TFile::Open(filename);
  if (!file->IsOpen()) {printf("Can't open file %s !\n",filename); return;} 

  Char_t treeName[12];
  sprintf(treeName,"TreeR%d_TRD",fEvent);
  TTree *clusterTree = (TTree*) file->Get(treeName);

  TObjArray *clusterArray = new TObjArray(400); 
 
  clusterTree->GetBranch("TRDcluster")->SetAddress(&clusterArray); 
  
  Int_t nEntries = (Int_t) clusterTree->GetEntries();
  printf("found %d entries in %s.\n",nEntries,clusterTree->GetName());

  // Loop through all entries in the tree
  Int_t nbytes;
  AliTRDcluster *c = 0;

  for (Int_t iEntry = 0; iEntry < nEntries; iEntry++) {    
    
    // Import the tree
    nbytes += clusterTree->GetEvent(iEntry);  

    // Get the number of points in the detector
    Int_t nCluster = clusterArray->GetEntriesFast();  
    printf("Read %d clusters from entry %d \n", nCluster, iEntry);

    // Loop through all TRD digits
    for (Int_t iCluster = 0; iCluster < nCluster; iCluster++) { 
      c = (AliTRDcluster*)clusterArray->UncheckedAt(iCluster);
      AliTRDcluster *co = new AliTRDcluster(*c);
      co->SetSigmaY2(c->GetSigmaY2() * fSY2corr);
      array->AddLast(co);
      delete clusterArray->RemoveAt(iCluster); 
    }
  }

  file->Close();                   
  delete clusterArray;
  savedir->cd(); 
  
}

//___________________________________________________________________
void AliTRDtracker::FindTracks(AliTRDtrackingSector* fTrSec, TH1F *hs, TH1F *hd) 
{
  //
  // Finds tracks in TRD
  //

  TH1F *hsame = hs;
  TH1F *hdiff = hd;   

  Int_t numOfTimeBins = fTrSec[0].GetNtimeBins(); 
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

    Int_t label = GetTrackLabel(t);

    if (FindProlongation(t,fTrSec,ns,0,label,hsame,hdiff)) {
      cerr<<"No of clusters in the track = "<<t.GetNclusters()<<endl; 
      if (t.GetNclusters() >= Int_t(fgkMinClustersInTrack*numOfTimeBins)) {
	Int_t label = GetTrackLabel(t);
	t.SetLabel(label);
	t.CookdEdx();
	UseClusters(t);

        AliTRDtrack *pt = new AliTRDtrack(t);
        fTracks->AddLast(pt); fNtracks++;     

	cerr<<"found track "<<fNtracks<<endl;
      }                         
    }     
    delete fSeeds->RemoveAt(i);  
    fNseeds--;
  }            
}

//__________________________________________________________________
void AliTRDtracker::UseClusters(AliTRDtrack t) 
{
  //
  // Mark used cluster
  //

  Int_t ncl=t.GetNclusters();
  for (Int_t i=0; i<ncl; i++) {
    Int_t index = t.GetClusterIndex(i);
    AliTRDcluster *c=(AliTRDcluster*)fClusters->UncheckedAt(index);
    c->Use();
  }

}

//__________________________________________________________________
Int_t AliTRDtracker::GetTrackLabel(AliTRDtrack t) 
{
  //
  // Get MC label
  //

  Int_t label=123456789, index, i, j;
  Int_t ncl=t.GetNclusters();
  const Int_t kRange = AliTRDgeometry::kNplan * fGeom->GetTimeMax();
  Bool_t labelAdded;

  //  Int_t s[kRange][2];
  Int_t **s = new Int_t* [kRange];
  for (i=0; i<kRange; i++) {
    s[i] = new Int_t[2];
  }
  for (i=0; i<kRange; i++) {
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
      labelAdded=kFALSE; j=0;
      if (label >= 0) {
	while ( (!labelAdded) && ( j < kRange ) ) {
	  if (s[j][0]==label || s[j][1]==0) {
	    s[j][0]=label; 
	    s[j][1]=s[j][1]+1; 
	    labelAdded=kTRUE;
	  }
	  j++;
	}
      }
    }
  }

  Int_t max=0;
  label = -123456789;

  for (i=0; i<kRange; i++) {
    if (s[i][1]>max) {
      max=s[i][1]; label=s[i][0];
    }
  }
  delete []s;
  if(max > ncl*fgkLabelFraction) return label;
  else return -1;
}

//___________________________________________________________________
Int_t AliTRDtracker::WriteTracks(const Char_t *filename) 
{
  //
  // Write the tracks to the output file
  //

  TDirectory *savedir=gDirectory;   

  //TFile *out=TFile::Open(filename,"RECREATE");
  TFile *out = (TFile*) gROOT->GetListOfFiles()->FindObject(filename);
  if (!out) {
    printf("AliTRDtracker::Open -- ");
    printf("Open the ALIROOT-file %s.\n",filename);
    out = new TFile(filename,"RECREATE");
  }
  else {
    printf("AliTRDtracker::Open -- ");
    printf("%s is already open.\n",filename);
  }

  Char_t treeName[12];
  sprintf(treeName,"TreeT%d_TRD",fEvent);
  TTree tracktree(treeName,"Tree with TRD tracks");

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
  return 1;

}

