///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  Time Projection Chamber                                                  //
//  This class contains the basic functions for the Time Projection Chamber  //
//  detector. Functions specific to one particular geometry are              //
//  contained in the derived classes                                         //
//                                                                           //
//Begin_Html
/*
<img src="picts/AliTPCClass.gif">
*/
//End_Html
//                                                                           //
//                                                                          //
///////////////////////////////////////////////////////////////////////////////

#include <TMath.h>
#include <TRandom.h>
#include <TVector.h>
#include <TMatrix.h>
#include <TGeometry.h>
#include <TNode.h>
#include <TTUBS.h>
#include <TObjectTable.h>
#include "TParticle.h"
#include "AliTPC.h"
#include "AliRun.h"
#include <iostream.h>
#include <fstream.h>
#include "AliMC.h"

//MI change
#include "AliTPCParam.h"
#include "AliTPCD.h"
#include "AliTPCPRF2D.h"
#include "AliTPCRF1D.h"


ClassImp(AliTPC) 

//_____________________________________________________________________________
AliTPC::AliTPC()
{
  //
  // Default constructor
  //
  fIshunt   = 0;
  fClusters = 0;
  fHits     = 0;
  fDigits   = 0;
  fTracks   = 0;
  fNsectors = 0;
  fNtracks  = 0;
  fNclusters= 0;
  //MI changes
  fDigParam= new AliTPCD();
  fDigits = fDigParam->GetArray();
}
 
//_____________________________________________________________________________
AliTPC::AliTPC(const char *name, const char *title)
      : AliDetector(name,title)
{
  //
  // Standard constructor
  //

  //
  // Initialise arrays of hits and digits 
  fHits     = new TClonesArray("AliTPChit",  176);
  //  fDigits   = new TClonesArray("AliTPCdigit",10000);
  //MI change
  fDigParam= new AliTPCD;
  fDigits = fDigParam->GetArray();
  //
  // Initialise counters
  fClusters = 0;
  fTracks   = 0;
  fNsectors = 72;
  fNtracks  = 0;
  fNclusters= 0;
  fDigitsIndex = new Int_t[fNsectors+1];
  fClustersIndex = new Int_t[fNsectors+1];
  //
  fIshunt     =  0;
  //
  // Initialise color attributes
  SetMarkerColor(kYellow);
}

//_____________________________________________________________________________
AliTPC::~AliTPC()
{
  //
  // TPC destructor
  //
  fIshunt   = 0;
  delete fHits;
  delete fDigits;
  delete fClusters;
  delete fTracks;
  delete fDigParam;
  if (fDigitsIndex)   delete [] fDigitsIndex;
  if (fClustersIndex) delete [] fClustersIndex;
}

//_____________________________________________________________________________
void AliTPC::AddCluster(Float_t *hits, Int_t *tracks)
{
  //
  // Add a simulated cluster to the list
  //
  if(!fClusters) fClusters=new TClonesArray("AliTPCcluster",10000);
  TClonesArray &lclusters = *fClusters;
  new(lclusters[fNclusters++]) AliTPCcluster(hits,tracks);
}
 
//_____________________________________________________________________________
void AliTPC::AddCluster(const AliTPCcluster &c)
{
  //
  // Add a simulated cluster copy to the list
  //
  if(!fClusters) fClusters=new TClonesArray("AliTPCcluster",10000);
  TClonesArray &lclusters = *fClusters;
  new(lclusters[fNclusters++]) AliTPCcluster(c);
}
 
//_____________________________________________________________________________
void AliTPC::AddDigit(Int_t *tracks, Int_t *digits)
{
  //
  // Add a TPC digit to the list
  //
  //  TClonesArray &ldigits = *fDigits;
  //MI change 
  TClonesArray &ldigits = *fDigParam->GetArray();
  new(ldigits[fNdigits++]) AliTPCdigit(tracks,digits);
}
 
//_____________________________________________________________________________
void AliTPC::AddHit(Int_t track, Int_t *vol, Float_t *hits)
{
  //
  // Add a hit to the list
  //
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliTPChit(fIshunt,track,vol,hits);
}
 
//_____________________________________________________________________________
void AliTPC::AddTrack(Float_t *hits)
{
  //
  // Add a track to the list of tracks
  //
  TClonesArray &ltracks = *fTracks;
  new(ltracks[fNtracks++]) AliTPCtrack(hits);
}

//_____________________________________________________________________________
void AliTPC::AddTrack(const AliTPCtrack& t)
{
  //
  // Add a track copy to the list of tracks
  //
  if(!fTracks) fTracks=new TClonesArray("AliTPCtrack",10000);
  TClonesArray &ltracks = *fTracks;
  new(ltracks[fNtracks++]) AliTPCtrack(t);
}

//_____________________________________________________________________________
void AliTPC::BuildGeometry()
{
  //
  // Build TPC ROOT TNode geometry for the event display
  //
  TNode *Node, *Top;
  TTUBS *tubs;
  Int_t i;
  const int kColorTPC=19;
  char name[5], title[20];
  const Double_t kDegrad=TMath::Pi()/180;
  const Double_t loAng=30;
  const Double_t hiAng=15;
  const Int_t nLo = Int_t (360/loAng+0.5);
  const Int_t nHi = Int_t (360/hiAng+0.5);
  const Double_t loCorr = 1/TMath::Cos(0.5*loAng*kDegrad);
  const Double_t hiCorr = 1/TMath::Cos(0.5*hiAng*kDegrad);
  //
  // Get ALICE top node
  Top=gAlice->GetGeometry()->GetNode("alice");
  //
  // Inner sectors
  for(i=0;i<nLo;i++) {
    sprintf(name,"LS%2.2d",i);
    sprintf(title,"TPC low sector %d",i);
    tubs = new TTUBS(name,title,"void",88*loCorr,136*loCorr,250,loAng*(i-0.5),loAng*(i+0.5));
    tubs->SetNumberOfDivisions(1);
    Top->cd();
    Node = new TNode(name,title,name,0,0,0,"");
    Node->SetLineColor(kColorTPC);
    fNodes->Add(Node);
  }
  // Outer sectors
  for(i=0;i<nHi;i++) {
    sprintf(name,"US%2.2d",i);
    sprintf(title,"TPC upper sector %d",i);
    tubs = new TTUBS(name,title,"void",142*hiCorr,250*hiCorr,250,hiAng*(i-0.5),hiAng*(i+0.5));
    tubs->SetNumberOfDivisions(1);
    Top->cd();
    Node = new TNode(name,title,name,0,0,0,"");
    Node->SetLineColor(kColorTPC);
    fNodes->Add(Node);
  }
}
//_____________________________________________________________________________
Int_t AliTPC::DistancetoPrimitive(Int_t , Int_t )
{
  //
  // Calculate distance from TPC to mouse on the display
  // Dummy procedure
  //
  return 9999;
}

//_____________________________________________________________________________
//const int MAX_CLUSTER=nrow_low+nrow_up; 
const int MAX_CLUSTER=200; 
const int S_MAXSEC=24;
const int L_MAXSEC=48;
const int ROWS_TO_SKIP=21;
const Float_t MAX_CHI2=12.;


//_____________________________________________________________________________
static Double_t SigmaY2(Double_t r, Double_t tgl, Double_t pt)
{
  //
  // Calculate spread in Y
  //
  pt=TMath::Abs(pt)*1000.;
  Double_t x=r/pt;
  tgl=TMath::Abs(tgl);
  Double_t s=a_rphi - b_rphi*r*tgl + c_rphi*x*x + d_rphi*x;
  if (s<0.4e-3) s=0.4e-3;
  return s;
}

//_____________________________________________________________________________
static Double_t SigmaZ2(Double_t r, Double_t tgl) 
{
  //
  // Calculate spread in Z
  //
  tgl=TMath::Abs(tgl);
  Double_t s=a_z - b_z*r*tgl + c_z*tgl*tgl;
  if (s<0.4e-3) s=0.4e-3;
  return s;
}

//_____________________________________________________________________________
inline Double_t f1(Double_t x1,Double_t y1,   //C
                   Double_t x2,Double_t y2,
                   Double_t x3,Double_t y3) 
{
  //
  // Function f1
  //
  Double_t d=(x2-x1)*(y3-y2)-(x3-x2)*(y2-y1);
  Double_t a=0.5*((y3-y2)*(y2*y2-y1*y1+x2*x2-x1*x1)-
                  (y2-y1)*(y3*y3-y2*y2+x3*x3-x2*x2));
  Double_t b=0.5*((x2-x1)*(y3*y3-y2*y2+x3*x3-x2*x2)-
                  (x3-x2)*(y2*y2-y1*y1+x2*x2-x1*x1));

  Double_t xr=TMath::Abs(d/(d*x1-a)), yr=d/(d*y1-b);

  return -xr*yr/sqrt(xr*xr+yr*yr); 
}


//_____________________________________________________________________________
inline Double_t f2(Double_t x1,Double_t y1,  //eta=C*x0
                   Double_t x2,Double_t y2,
                   Double_t x3,Double_t y3) 
{
  //
  // Function f2
  //
  Double_t d=(x2-x1)*(y3-y2)-(x3-x2)*(y2-y1);
  Double_t a=0.5*((y3-y2)*(y2*y2-y1*y1+x2*x2-x1*x1)-
                  (y2-y1)*(y3*y3-y2*y2+x3*x3-x2*x2));
  Double_t b=0.5*((x2-x1)*(y3*y3-y2*y2+x3*x3-x2*x2)-
                  (x3-x2)*(y2*y2-y1*y1+x2*x2-x1*x1));

  Double_t xr=TMath::Abs(d/(d*x1-a)), yr=d/(d*y1-b);
  
  return -a/(d*y1-b)*xr/sqrt(xr*xr+yr*yr);
}

//_____________________________________________________________________________
inline Double_t f3(Double_t x1,Double_t y1,  //tgl
                   Double_t x2,Double_t y2,
                   Double_t z1,Double_t z2) 
{
  //
  // Function f3
  //
  return (z1 - z2)/sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
}

//_____________________________________________________________________________
static int FindProlongation(AliTPCtrack& t, const AliTPCSector *sec, 
			    int s, int ri, int rf=0) 
{
  //
  // Propagate track
  //
  int try_again=ROWS_TO_SKIP;
  Double_t alpha=sec->GetAlpha();
  int ns=int(2*TMath::Pi()/alpha)+1;

  for (int nr=ri; nr>=rf; nr--) {
    Double_t x=sec[s].GetX(nr), ymax=sec[s].GetMaxY(nr);
    if (!t.PropagateTo(x)) return -1;

    const AliTPCcluster *cl=0;
    Double_t max_chi2=MAX_CHI2;
    const AliTPCRow& row=sec[s][nr];
    Double_t sy2=SigmaY2(t.GetX(),t.GetTgl(),t.GetPt());
    Double_t sz2=SigmaZ2(t.GetX(),t.GetTgl());
    Double_t road=3.*sqrt(t.GetSigmaY2() + 4*sy2), y=t.GetY(), z=t.GetZ();

    if (road>30) {
      if (t>3) cerr<<t<<" AliTPCtrack warning: Too broad road !\n"; 
      return -1;
    }

    if (row) {
      for (int i=row.Find(y-road); i<row; i++) {
	AliTPCcluster* c=(AliTPCcluster*)(row[i]);
	if (c->fY > y+road) break;
	if (c->IsUsed()) continue;
	if ((c->fZ - z)*(c->fZ - z) > 9.*(t.GetSigmaZ2() + 4*sz2)) continue;
	Double_t chi2=t.GetPredictedChi2(c);
	if (chi2 > max_chi2) continue;
	max_chi2=chi2;
	cl=c;       
      }
    }
    if (cl) {
      t.Update(cl,max_chi2);
      try_again=ROWS_TO_SKIP;
    } else {
      if (try_again==0) break;
      if (y > ymax) {
         s = (s+1) % ns;
         if (!t.Rotate(alpha)) return -1;
      } else if (y <-ymax) {
           s = (s-1+ns) % ns;
           if (!t.Rotate(-alpha)) return -1;
      }
      try_again--;
    }
  }

  return s;
}


//_____________________________________________________________________________
static void MakeSeeds(TObjArray& seeds,const AliTPCSector* sec,int i1,int i2,
const AliTPCParam *p)
{
  //
  // Find seed for tracking
  //
  TMatrix C(5,5); TVector x(5);
  int max_sec=L_MAXSEC/2;
  for (int ns=0; ns<max_sec; ns++) {
    int nl=sec[(ns-1+max_sec)%max_sec][i2];
    int nm=sec[ns][i2];
    int nu=sec[(ns+1)%max_sec][i2];
    Double_t alpha=sec[ns].GetAlpha();
    const AliTPCRow& r1=sec[ns][i1];
    for (int is=0; is < r1; is++) {
      Double_t x1=sec[ns].GetX(i1), y1=r1[is]->fY, z1=r1[is]->fZ;
      for (int js=0; js < nl+nm+nu; js++) {
	const AliTPCcluster *cl;
	Double_t cs,sn;
	int ks;
	
	if (js<nl) {
	  ks=(ns-1+max_sec)%max_sec;
	  const AliTPCRow& r2=sec[(ns-1+max_sec)%max_sec][i2];
	  cl=r2[js];
	  cs=cos(alpha); sn=sin(alpha);
	} else 
	  if (js<nl+nm) {
	    ks=ns;
	    const AliTPCRow& r2=sec[ns][i2];
	    cl=r2[js-nl];
	    cs=1; sn=0.;
	  } else {
	    ks=(ns+1)%max_sec;
	    const AliTPCRow& r2=sec[(ns+1)%max_sec][i2];
	    cl=r2[js-nl-nm];
	    cs=cos(alpha); sn=-sin(alpha);
	  }
	Double_t x2=sec[ns].GetX(i2), y2=cl->fY, z2=cl->fZ;
	//if (z1*z2 < 0) continue;
	//if (TMath::Abs(z1) < TMath::Abs(z2)) continue;
	
	Double_t tmp= x2*cs+y2*sn;
	y2 =-x2*sn+y2*cs;
	x2=tmp;     
	
	x(0)=y1;
	x(1)=z1;
	x(2)=f1(x1,y1,x2,y2,0.,0.);
	x(3)=f2(x1,y1,x2,y2,0.,0.);
	x(4)=f3(x1,y1,x2,y2,z1,z2);
	
	if (TMath::Abs(x(2)*x1-x(3)) >= 0.999) continue;
	
	if (TMath::Abs(x(4)) > 1.2) continue;

	Double_t a=asin(x(3));
	Double_t zv=z1 - x(4)/x(2)*(a+asin(x(2)*x1-x(3)));
	if (TMath::Abs(zv)>33.) continue; 

	TMatrix X(6,6); X=0.; 
	X(0,0)=r1[is]->fSigmaY2; X(1,1)=r1[is]->fSigmaZ2;
	X(2,2)=cl->fSigmaY2;     X(3,3)=cl->fSigmaZ2;
	X(4,4)=3./12.; X(5,5)=3./12.;
	TMatrix F(5,6); F.UnitMatrix();
	Double_t sy=sqrt(X(0,0)), sz=sqrt(X(1,1));
	F(2,0)=(f1(x1,y1+sy,x2,y2,0.,0.)-x(2))/sy;
	F(2,2)=(f1(x1,y1,x2,y2+sy,0.,0.)-x(2))/sy;
	F(2,4)=(f1(x1,y1,x2,y2,0.,0.+sy)-x(2))/sy;
	F(3,0)=(f2(x1,y1+sy,x2,y2,0.,0.)-x(3))/sy;
	F(3,2)=(f2(x1,y1,x2,y2+sy,0.,0.)-x(3))/sy;
	F(3,4)=(f2(x1,y1,x2,y2,0.,0.+sy)-x(3))/sy;
	F(4,0)=(f3(x1,y1+sy,x2,y2,z1,z2)-x(4))/sy;
	F(4,1)=(f3(x1,y1,x2,y2,z1+sz,z2)-x(4))/sz;
	F(4,2)=(f3(x1,y1,x2,y2+sy,z1,z2)-x(4))/sy;
	F(4,3)=(f3(x1,y1,x2,y2,z1,z2+sz)-x(4))/sz;
	F(4,4)=0;
	F(3,3)=0;
	
	TMatrix t(F,TMatrix::kMult,X);
	C.Mult(t,TMatrix(TMatrix::kTransposed,F));

	TrackSeed *track=new TrackSeed(*(r1[is]),x,C,p);
	int rc=FindProlongation(*track,sec,ns,i1-1,i2);
        if (rc<0 || *track<(i1-i2)/2) delete track;
        else seeds.AddLast(track); 
      }
    }
  }
}

//_____________________________________________________________________________
void AliTPC::Clusters2Tracks()
{
  //
  // TPC Track finder from clusters.
  //
  if (!fClusters) return;

  AliTPCParam *p=&fDigParam->GetParam();
  Int_t nrow_low=p->GetNRowLow();
  Int_t nrow_up=p->GetNRowUp();

  AliTPCSSector ssec[S_MAXSEC/2];
  for (int i=0; i<S_MAXSEC/2; i++) ssec[i].SetUp(p);

  AliTPCLSector lsec[L_MAXSEC/2];
  for (int j=0; j<L_MAXSEC/2; j++) lsec[j].SetUp(p);

  int ncl=fClusters->GetEntriesFast();
  while (ncl--) {
    AliTPCcluster *c=(AliTPCcluster*)fClusters->UncheckedAt(ncl);

    int sec=int(c->fSector)-1, row=int(c->fPadRow)-1;
    
    if (sec<24) {
      if (row<0 || row>nrow_low) {cerr<<"low !!!"<<row<<endl; continue;}
      ssec[sec%12][row].InsertCluster(c);
    } else {
      if (row<0 || row>nrow_up ) {cerr<<"up  !!!"<<row<<endl; continue;}
      sec -= 24;
      lsec[sec%24][row].InsertCluster(c);
    }
  }
  
  
  TObjArray seeds(20000);
  MakeSeeds(seeds,lsec,nrow_up-1,nrow_up-1-8,p);
  MakeSeeds(seeds,lsec,nrow_up-1-4,nrow_up-1-4-8,p);
    
  seeds.Sort();
  
  int found=0;
  int nseed=seeds.GetEntriesFast();
  
  for (int s=0; s<nseed; s++) {
    AliTPCtrack& t=*((AliTPCtrack*)seeds.UncheckedAt(s));
    Double_t alpha=t.GetAlpha();
    if (alpha > 2.*TMath::Pi()) alpha -= 2.*TMath::Pi();  
    if (alpha < 0.            ) alpha += 2.*TMath::Pi();  
    int ns=int(alpha/lsec->GetAlpha() + 0.5);
    
    Double_t x=t.GetX();
    int nr;
    if (x<p->GetPadRowRadiiUp(nrow_up-1-4-7)) nr=nrow_up-1-4-8;
    else if (x<p->GetPadRowRadiiUp(nrow_up-1-7)) nr=nrow_up-1-8;
    else {cerr<<x<<" =x !!!\n"; continue;}

    int ls=FindProlongation(t,lsec,ns,nr-1);
    if (ls<0) continue;
    x=t.GetX(); alpha=lsec[ls].GetAlpha();          //
    Double_t phi=ls*alpha + atan(t.GetY()/x);       // Find S-sector
    int ss=int(0.5*(phi/alpha+1));                  //
    alpha *= 2*(ss-0.5*ls);                         // and rotation angle
    if (!t.Rotate(alpha)) continue;                 //
    ss %= (S_MAXSEC/2);                             //
    
    if (FindProlongation(t,ssec,ss,nrow_low-1)<0) continue;
    if (t < 30) continue;
    
    AddTrack(t);
    t.UseClusters();
    cerr<<found++<<'\r';
  }  
}

//_____________________________________________________________________________
void AliTPC::CreateMaterials()
{
  //-----------------------------------------------
  // Create Materials for for TPC
  //-----------------------------------------------

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------
  
  AliMC* pMC = AliMC::GetMC();

  Int_t ISXFLD=gAlice->Field()->Integ();
  Float_t SXMGMX=gAlice->Field()->Max();
  
  Float_t absl, radl, a, d, z;
  Float_t dg;
  Float_t x0ne;
  Float_t buf[1];
  Int_t nbuf;
  
  // --- Methane (CH4) --- 
  Float_t am[2] = { 12.,1. };
  Float_t zm[2] = { 6.,1. };
  Float_t wm[2] = { 1.,4. };
  Float_t dm    = 7.17e-4;
  // --- The Neon CO2 90/10 mixture --- 
  Float_t ag[2] = { 20.18 };
  Float_t zg[2] = { 10. };
  Float_t wg[2] = { .8,.2 };
  Float_t dne   = 9e-4;        // --- Neon density in g/cm3 ---
  
  // --- Mylar (C5H4O2) --- 
  Float_t amy[3] = { 12.,1.,16. };
  Float_t zmy[3] = { 6.,1.,8. };
  Float_t wmy[3] = { 5.,4.,2. };
  Float_t dmy    = 1.39;
  // --- CO2 --- 
  Float_t ac[2] = { 12.,16. };
  Float_t zc[2] = { 6.,8. };
  Float_t wc[2] = { 1.,2. };
  Float_t dc    = .001977;
  // --- Carbon density and radiation length --- 
  Float_t densc = 2.265;
  Float_t radlc = 18.8;
  // --- Silicon --- 
  Float_t asi   = 28.09;
  Float_t zsi   = 14.;
  Float_t desi  = 2.33;
  Float_t radsi = 9.36;
  
  // --- Define the various materials for GEANT --- 
  AliMaterial(0, "Al $", 26.98, 13., 2.7, 8.9, 37.2);
  x0ne = 28.94 / dne;
  AliMaterial(1, "Ne $", 20.18, 10., dne, x0ne, 999.);
  
  // --  Methane, defined by the proportions of atoms 
  
  AliMixture(2, "Methane$", am, zm, dm, -2, wm);
  
  // --- CO2, defined by the proportion of atoms 
  
  AliMixture(7, "CO2$", ac, zc, dc, -2, wc);
  
  // --  Get A,Z etc. for CO2 
  
  char namate[21];
  pMC->Gfmate((*fIdmate)[7], namate, a, z, d, radl, absl, buf, nbuf);
  ag[1] = a;
  zg[1] = z;
  dg = dne * .9 + dc * .1;
  
  // --  Create Ne/CO2 90/10 mixture 
  
  AliMixture(3, "Gas-mixt $", ag, zg, dg, 2, wg);
  AliMixture(4, "Gas-mixt $", ag, zg, dg, 2, wg);
  
  AliMaterial(5, "G10$", 20., 10., 1.7, 19.4, 999.);
  AliMixture(6, "Mylar$", amy, zmy, dmy, -3, wmy);
  
  a = ac[0];
  z = zc[0];
  AliMaterial(8, "Carbon", a, z, densc, radlc, 999.);
  
  AliMaterial(9, "Silicon", asi, zsi, desi, radsi, 999.);
  AliMaterial(99, "Air$", 14.61, 7.3, .001205, 30420., 67500.);
  
  AliMedium(400, "Al wall$",  0, 0, ISXFLD, SXMGMX, 10., .1, .1, .1,   .1);
  AliMedium(402, "Gas mix1$", 3, 0, ISXFLD, SXMGMX, 10., .01,.1, .001, .01);
  AliMedium(403, "Gas mix2$", 3, 0, ISXFLD, SXMGMX, 10., .01,.1, .001, .01);
  AliMedium(404, "Gas mix3$", 4, 1, ISXFLD, SXMGMX, 10., .01,.1, .001, .01);
  AliMedium(405, "G10 pln$",  5, 0, ISXFLD, SXMGMX, 10., .1, .1, .1,   .1 );
  AliMedium(406, "Mylar  $",  6, 0, ISXFLD, SXMGMX, 10., .01,.1, .001, .01);
  AliMedium(407, "CO2    $",  7, 0, ISXFLD, SXMGMX, 10., .01,.1, .01,  .01);
  AliMedium(408, "Carbon $",  8, 0, ISXFLD, SXMGMX, 10., .1, .1, .1,   .1 );
  AliMedium(409, "Silicon$",  9, 0, ISXFLD, SXMGMX, 10., .1, .1, .1,   .1 );
  AliMedium(499, "Air gap$", 99, 0, ISXFLD, SXMGMX, 10., .1, .1, .1,   .1 );
}

//_____________________________________________________________________________
struct Bin {
   const AliTPCdigit *dig;
   int idx;
   Bin() {dig=0; idx=-1;}
};

struct PreCluster : public AliTPCcluster {
   const AliTPCdigit* summit;
   int idx;
   int cut;
   int npeaks;
   PreCluster();
};
PreCluster::PreCluster() : AliTPCcluster() {cut=npeaks=0;}


//_____________________________________________________________________________
static void FindCluster(int i, int j, Bin bins[][MAXTPCTBK+2], PreCluster &c) 
{
  //
  // Find clusters
  //
  Bin& b=bins[i][j];
  double q=double(b.dig->fSignal);

  if (q<0) { // digit is at the edge of the pad row
    q=-q;
    c.cut=1;
  } 
  if (b.idx >= 0 && b.idx != c.idx) {
    c.idx=b.idx;
    c.npeaks++;
  }
  
  if (q > TMath::Abs(c.summit->fSignal)) c.summit=b.dig;
  
  c.fY += i*q;
  c.fZ += j*q;
  c.fSigmaY2 += i*i*q;
  c.fSigmaZ2 += j*j*q;
  c.fQ += q;

  b.dig = 0;  b.idx = c.idx;
  
  if (bins[i-1][j].dig) FindCluster(i-1,j,bins,c);
  if (bins[i][j-1].dig) FindCluster(i,j-1,bins,c);
  if (bins[i+1][j].dig) FindCluster(i+1,j,bins,c);
  if (bins[i][j+1].dig) FindCluster(i,j+1,bins,c);
  
}

//_____________________________________________________________________________
void AliTPC::Digits2Clusters()
{
  //
  // simple TPC cluster finder from digits.
  //
  //
  AliTPCParam * fTPCParam = &(fDigParam->GetParam());
  
  const Int_t MAX_PAD=200+2, MAX_BUCKET=MAXTPCTBK+2;
  const Int_t Q_min=60;
  const Int_t THRESHOLD=20;
  
  TTree *t=(TTree*)gDirectory->Get("TreeD0_Param1");
  t->GetBranch("Digits")->SetAddress(&fDigits);
  Int_t sectors_by_rows=(Int_t)t->GetEntries();
  
  int ncls=0;
  
  for (Int_t n=0; n<sectors_by_rows; n++) {
    if (!t->GetEvent(n)) continue;
    Bin bins[MAX_PAD][MAX_BUCKET];
    AliTPCdigit *dig=(AliTPCdigit*)fDigits->UncheckedAt(0);
    Int_t nsec=dig->fSector, nrow=dig->fPadRow;
    Int_t ndigits=fDigits->GetEntriesFast();
    
    int npads;  int sign_z;
    if (nsec<25) {
      sign_z=(nsec<13) ? 1 : -1;
      npads=fTPCParam->GetNPadsLow(nrow-1);
    } else {
      sign_z=(nsec<49) ? 1 : -1;
      npads=fTPCParam->GetNPadsUp(nrow-1);
    }
    
    int ndig;
    for (ndig=0; ndig<ndigits; ndig++) {
      dig=(AliTPCdigit*)fDigits->UncheckedAt(ndig);
      int i=dig->fPad, j=dig->fTime;
      if (dig->fSignal >= THRESHOLD) bins[i][j].dig=dig;
      if (i==1 || i==npads || j==1 || j==MAXTPCTBK) dig->fSignal*=-1;
    }
    
    int ncl=0;
    int i,j;
    
    for (i=1; i<MAX_PAD-1; i++) {
      for (j=1; j<MAX_BUCKET-1; j++) {
	if (bins[i][j].dig == 0) continue;
	PreCluster c; c.summit=bins[i][j].dig; c.idx=ncls;
	FindCluster(i, j, bins, c);
	c.fY /= c.fQ;
	c.fZ /= c.fQ;

        double s2 = c.fSigmaY2/c.fQ - c.fY*c.fY;
        c.fSigmaY2 = s2 + 1./12.;
        c.fSigmaY2 *= fTPCParam->GetPadPitchWidth()*
                      fTPCParam->GetPadPitchWidth();
        if (s2 != 0.) c.fSigmaY2 *= 0.022*8*4;

        s2 = c.fSigmaZ2/c.fQ - c.fZ*c.fZ;
        c.fSigmaZ2 = s2 + 1./12.;
        c.fSigmaZ2 *= fTPCParam->GetZWidth()*fTPCParam->GetZWidth();
        if (s2 != 0.) c.fSigmaZ2 *= 0.068*4*4;

	c.fY = (c.fY - 0.5 - 0.5*npads)*fTPCParam->GetPadPitchWidth();
	c.fZ = fTPCParam->GetZWidth()*(c.fZ+1); 
	c.fZ -= 3.*fTPCParam->GetZSigma(); // PASA delay 
	c.fZ = sign_z*(z_end - c.fZ);
	//c.fZ += 0.023;
	c.fSector=nsec;
	c.fPadRow=nrow;
	c.fTracks[0]=c.summit->fTracks[0];
	c.fTracks[1]=c.summit->fTracks[1];
	c.fTracks[2]=c.summit->fTracks[2];

	if (c.cut) {
	  c.fSigmaY2 *= 25.;
	  c.fSigmaZ2 *= 4.;
	}
	
	AddCluster(c); ncls++; ncl++;
      }
    }
    
    for (ndig=0; ndig<ndigits; ndig++) {
      dig=(AliTPCdigit*)fDigits->UncheckedAt(ndig);
      if (TMath::Abs(dig->fSignal) >= 0) 
	bins[dig->fPad][dig->fTime].dig=dig;
    }
    
    for (i=1; i<MAX_PAD-1; i++) {
      for (j=1; j<MAX_BUCKET-1; j++) {
	if (bins[i][j].dig == 0) continue;
	PreCluster c; c.summit=bins[i][j].dig; c.idx=ncls;
	FindCluster(i, j, bins, c);
	if (c.fQ <= Q_min) continue; //noise cluster
        if (c.npeaks>1) continue;    //overlapped cluster
	c.fY /= c.fQ;
	c.fZ /= c.fQ;

        double s2 = c.fSigmaY2/c.fQ - c.fY*c.fY;
        c.fSigmaY2 = s2 + 1./12.;
        c.fSigmaY2 *= fTPCParam->GetPadPitchWidth()*
                      fTPCParam->GetPadPitchWidth();
        if (s2 != 0.) c.fSigmaY2 *= 0.022*4*0.6*4;

        s2 = c.fSigmaZ2/c.fQ - c.fZ*c.fZ;
        c.fSigmaZ2 = s2 + 1./12.;
        c.fSigmaZ2 *= fTPCParam->GetZWidth()*fTPCParam->GetZWidth();
        if (s2 != 0.) c.fSigmaZ2 *= 0.068*4*0.4;

	c.fY = (c.fY - 0.5 - 0.5*npads)*fTPCParam->GetPadPitchWidth();
	c.fZ = fTPCParam->GetZWidth()*(c.fZ+1); 
	c.fZ -= 3.*fTPCParam->GetZSigma(); // PASA delay 
	c.fZ = sign_z*(z_end - c.fZ);
	//c.fZ += 0.023;
	c.fSector=nsec;
	c.fPadRow=nrow;
	c.fTracks[0]=c.summit->fTracks[0];
	c.fTracks[1]=c.summit->fTracks[1];
	c.fTracks[2]=c.summit->fTracks[2];
	
	if (c.cut) {
	  c.fSigmaY2 *= 25.;
	  c.fSigmaZ2 *= 4.;
	}
	
	if (c.npeaks==0) {AddCluster(c); ncls++; ncl++;}
	else {
	  new ((*fClusters)[c.idx]) AliTPCcluster(c);
	}
      }
    }
    
    cerr<<"sector, row, digits, clusters: "
	<<nsec<<' '<<nrow<<' '<<ndigits<<' '<<ncl<<"                  \r";
    
    fDigits->Clear();
    
  }
}

//_____________________________________________________________________________
void AliTPC::ElDiff(Float_t *xyz)
{
  //--------------------------------------------------
  // calculates the diffusion of a single electron
  //--------------------------------------------------

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------
  AliTPCParam * fTPCParam = &(fDigParam->GetParam());
  Float_t driftl;
  
  Float_t z0=xyz[2];

  driftl=z_end-TMath::Abs(xyz[2]);

  if(driftl<0.01) driftl=0.01;

  // check the attachment

  driftl=TMath::Sqrt(driftl);

  //  Float_t sig_t = driftl*diff_t;
  //Float_t sig_l = driftl*diff_l;
  Float_t sig_t = driftl*fTPCParam->GetDiffT();
  Float_t sig_l = driftl*fTPCParam->GetDiffL();
  xyz[0]=gRandom->Gaus(xyz[0],sig_t);
  xyz[1]=gRandom->Gaus(xyz[1],sig_t);
  xyz[2]=gRandom->Gaus(xyz[2],sig_l);
  
  if (TMath::Abs(xyz[2])>z_end){
    xyz[2]=TMath::Sign(z_end,z0);
  }
  if(xyz[2]*z0 < 0.){
    Float_t eps = 0.0001;
    xyz[2]=TMath::Sign(eps,z0);
  } 
}

//_____________________________________________________________________________
void AliTPC::Hits2Clusters()
{
  //--------------------------------------------------------
  // TPC simple cluster generator from hits
  // obtained from the TPC Fast Simulator
  // The point errors are taken from the parametrization
  //--------------------------------------------------------

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------
  AliTPCParam * fTPCParam = &(fDigParam->GetParam());
  Float_t sigma_rphi,sigma_z,cl_rphi,cl_z;
  //
  TParticle *particle; // pointer to a given particle
  AliTPChit *tpcHit; // pointer to a sigle TPC hit
  TClonesArray *Particles; //pointer to the particle list
  Int_t sector,nhits;
  Int_t ipart;
  Float_t xyz[5];
  Float_t pl,pt,tanth,rpad,ratio;
  Float_t cph,sph;
  
  //---------------------------------------------------------------
  //  Get the access to the tracks 
  //---------------------------------------------------------------
  
  TTree *TH = gAlice->TreeH();
  Stat_t ntracks = TH->GetEntries();
  
  //------------------------------------------------------------
  // Loop over all sectors (72 sectors)
  // Sectors 1-24 are lower sectors, 1-12 z>0, 13-24 z<0
  // Sectors 25-72 are upper sectors, 25-48 z>0, 49-72 z<0
  //
  // First cluster for sector 1 starts at "0"
  //------------------------------------------------------------
  
  
  fClustersIndex[0] = 0;
  
  //
  for(Int_t isec=1;isec<fNsectors+1;isec++){
    //MI change
    fTPCParam->AdjustAngles(isec,cph,sph);
    
    //------------------------------------------------------------
    // Loop over tracks
    //------------------------------------------------------------
    
    for(Int_t track=0;track<ntracks;track++){
      ResetHits();
      TH->GetEvent(track);
      //
      //  Get number of the TPC hits and a pointer
      //  to the particles
      //
      nhits=fHits->GetEntriesFast();
      Particles=gAlice->Particles();
      //
      // Loop over hits
      //
      for(Int_t hit=0;hit<nhits;hit++){
	tpcHit=(AliTPChit*)fHits->UncheckedAt(hit);
	sector=tpcHit->fSector; // sector number
	if(sector != isec) continue; //terminate iteration
	ipart=tpcHit->fTrack;
	particle=(TParticle*)Particles->UncheckedAt(ipart);
	pl=particle->Pz();
	pt=particle->Pt();
	if(pt < 1.e-9) pt=1.e-9;
	tanth=pl/pt;
	tanth = TMath::Abs(tanth);
	rpad=TMath::Sqrt(tpcHit->fX*tpcHit->fX + tpcHit->fY*tpcHit->fY);
	ratio=0.001*rpad/pt; // pt must be in MeV/c - historical reason
	
	//   space-point resolutions
	
	sigma_rphi=SigmaY2(rpad,tanth,pt);
	sigma_z   =SigmaZ2(rpad,tanth   );
	
	//   cluster widths
	
	cl_rphi=ac_rphi-bc_rphi*rpad*tanth+cc_rphi*ratio*ratio;
	cl_z=ac_z-bc_z*rpad*tanth+cc_z*tanth*tanth;
	
	// temporary protection
	
	if(sigma_rphi < 0.) sigma_rphi=0.4e-3;
	if(sigma_z < 0.) sigma_z=0.4e-3;
	if(cl_rphi < 0.) cl_rphi=2.5e-3;
	if(cl_z < 0.) cl_z=2.5e-5;
	
	//
	
	//
	// smearing --> rotate to the 1 (13) or to the 25 (49) sector,
	// then the inaccuracy in a X-Y plane is only along Y (pad row)!
	//
	Float_t xprim= tpcHit->fX*cph + tpcHit->fY*sph;
	Float_t yprim=-tpcHit->fX*sph + tpcHit->fY*cph;
	xyz[0]=gRandom->Gaus(yprim,TMath::Sqrt(sigma_rphi));   // y
        Double_t alpha=(sector<25) ? alpha_low : alpha_up;
	if (TMath::Abs(xyz[0]/xprim) > TMath::Tan(0.5*alpha)) xyz[0]=yprim;
	xyz[1]=gRandom->Gaus(tpcHit->fZ,TMath::Sqrt(sigma_z)); // z 
        if (TMath::Abs(xyz[1]) > 250) xyz[1]=tpcHit->fZ;
	xyz[2]=tpcHit->fQ;                                     // q
	xyz[3]=sigma_rphi;                                     // fSigmaY2
	xyz[4]=sigma_z;                                        // fSigmaZ2
	
	//find row number
	//MI we must change
	Int_t row = fTPCParam->GetPadRow(sector,xprim) ;	
	// and finally add the cluster
	Int_t tracks[5]={tpcHit->fTrack, -1, -1, sector, row+1};
	AddCluster(xyz,tracks);
	
      } // end of loop over hits
    }   // end of loop over tracks 
    
    fClustersIndex[isec] = fNclusters; // update clusters index
    
  } // end of loop over sectors
  
  fClustersIndex[fNsectors]--; // set end of the clusters buffer
  
}


void AliTPC::Hits2Digits()  
{ 

 //----------------------------------------------------
  // Loop over all sectors (72 sectors)
  // Sectors 1-24 are lower sectors, 1-12 z>0, 13-24 z<0
  // Sectors 25-72 are upper sectors, 25-48 z>0, 49-72 z<0
  //----
  for(Int_t isec=1;isec<fNsectors+1;isec++)  Hits2DigitsSector(isec);
}


//_____________________________________________________________________________
void AliTPC::Hits2DigitsSector(Int_t isec)
{
  //-------------------------------------------------------------------
  // TPC conversion from hits to digits.
  //------------------------------------------------------------------- 

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------

  //-------------------------------------------------------
  //  Get the access to the track hits
  //-------------------------------------------------------

  AliTPCParam * fTPCParam = &(fDigParam->GetParam());
  TTree *TH = gAlice->TreeH(); // pointer to the hits tree
  Stat_t ntracks = TH->GetEntries();

  if( ntracks > 0){

  //------------------------------------------- 
  //  Only if there are any tracks...
  //-------------------------------------------

    
    // TObjArrays for three neighouring pad-rows

    TObjArray **rowTriplet = new TObjArray* [3]; 
    
    // TObjArray-s for each pad-row

    TObjArray **row;
      
    for (Int_t trip=0;trip<3;trip++){  
      rowTriplet[trip]=new TObjArray;
    }


    
      printf("*** Processing sector number %d ***\n",isec);

      Int_t nrows =fTPCParam->GetNRow(isec);

      row= new TObjArray* [nrows];
    
      MakeSector(isec,nrows,TH,ntracks,row);

      //--------------------------------------------------------
      //   Digitize this sector, row by row
      //   row[i] is the pointer to the TObjArray of TVectors,
      //   each one containing electrons accepted on this
      //   row, assigned into tracks
      //--------------------------------------------------------

      Int_t i;

      for (i=0;i<nrows;i++){

	// Triplets for i = 0 and i=1 are identical!
	// The same for i = nrows-1 and nrows!

        if(i != 1 && i != nrows-1){
	   MakeTriplet(i,rowTriplet,row);
	 }

	    DigitizeRow(i,isec,rowTriplet);

	    fDigParam->Fill();

            Int_t ndig=fDigParam->GetArray()->GetEntriesFast();

            printf("*** Sector, row, digits %d %d %d ***\n",isec,i,ndig);

	    ResetDigits(); // reset digits for this row after storing them
             
       } // end of the sector digitization
     
       // delete the last triplet

       for (i=0;i<3;i++) rowTriplet[i]->Delete();
          
       delete [] row; // delete the array of pointers to TObjArray-s
        
  } // ntracks >0
} // end of Hits2Digits
//_____________________________________________________________________________
void AliTPC::MakeTriplet(Int_t row,
                         TObjArray **rowTriplet, TObjArray **prow)
{
  //------------------------------------------------------------------
  //  Makes the "triplet" of the neighbouring pad-row for the
  //  digitization including the cross-talk between the pad-rows
  //------------------------------------------------------------------

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------

  AliTPCParam * fTPCParam = &(fDigParam->GetParam());
  Float_t gasgain = fTPCParam->GetGasGain();
  Int_t nTracks[3];

  Int_t nElements,nElectrons;

  TVector *pv;
  TVector *tv;

  //-------------------------------------------------------------------
  // pv is an "old" track, i.e. label + triplets of (x,y,z) 
  // for each electron
  //
  //-------------------------------------------------------------------


  Int_t i1,i2;
  Int_t nel,nt;

  if(row == 0 || row == 1){

  // create entire triplet for the first AND the second row

    nTracks[0] = prow[0]->GetEntries();
    nTracks[1] = prow[1]->GetEntries();
    nTracks[2] = prow[2]->GetEntries();

    for(i1=0;i1<3;i1++){
      nt = nTracks[i1]; // number of tracks for this row

      for(i2=0;i2<nt;i2++){
        pv = (TVector*)prow[i1]->At(i2);
        TVector &v1 = *pv;
        nElements = pv->GetNrows(); 
        nElectrons = (nElements-1)/3;

        tv = new TVector(4*nElectrons+1); // create TVector for a modified track
        TVector &v2 = *tv;
        v2(0)=v1(0); //track label

        for(nel=0;nel<nElectrons;nel++){        
          Int_t idx1 = nel*3;
          Int_t idx2 = nel*4;
       	  // Avalanche, including fluctuations
          Int_t aval = (Int_t) (-gasgain*TMath::Log(gRandom->Rndm()));
          v2(idx2+1)= v1(idx1+1);
          v2(idx2+2)= v1(idx1+2);
          v2(idx2+3)= v1(idx1+3);
          v2(idx2+4)= (Float_t)aval; // in number of electrons!        
	} // end of loop over electrons
	//
	//  Add this track to a row 
	//

        rowTriplet[i1]->Add(tv); 


      } // end of loop over tracks for this row

      prow[i1]->Delete(); // remove "old" tracks
      delete prow[i1]; // delete  TObjArray for this row
      prow[i1]=0; // set pointer to NULL

    } // end of loop over row triplets


  }
  else{
   
    rowTriplet[0]->Delete(); // remove old lower row

    nTracks[0]=rowTriplet[1]->GetEntries(); // previous middle row
    nTracks[1]=rowTriplet[2]->GetEntries(); // previous upper row
    nTracks[2]=prow[row+1]->GetEntries(); // next row
    

    //------------------------------------------- 
    //  shift new tracks downwards
    //-------------------------------------------

    for(i1=0;i1<nTracks[0];i1++){
      pv=(TVector*)rowTriplet[1]->At(i1);
      rowTriplet[0]->Add(pv); 
    }       

    rowTriplet[1]->Clear(); // leave tracks on the heap!!!

    for(i1=0;i1<nTracks[1];i1++){
      pv=(TVector*)rowTriplet[2]->At(i1);
      rowTriplet[1]->Add(pv);
    }

    rowTriplet[2]->Clear(); // leave tracks on the heap!!!

    //---------------------------------------------
    //  Create new upper row
    //---------------------------------------------

    

    for(i1=0;i1<nTracks[2];i1++){
        pv = (TVector*)prow[row+1]->At(i1);
        TVector &v1 = *pv;
        nElements = pv->GetNrows();
        nElectrons = (nElements-1)/3;

        tv = new TVector(4*nElectrons+1); // create TVector for a modified track
        TVector &v2 = *tv;
        v2(0)=v1(0); //track label

        for(nel=0;nel<nElectrons;nel++){
        
          Int_t idx1 = nel*3;
          Int_t idx2 = nel*4;
	  // Avalanche, including fluctuations
          Int_t aval = (Int_t) 
	    (-gasgain*TMath::Log(gRandom->Rndm()));
          
          v2(idx2+1)= v1(idx1+1); 
          v2(idx2+2)= v1(idx1+2);
          v2(idx2+3)= v1(idx1+3);
          v2(idx2+4)= (Float_t)aval; // in number of electrons!        
	} // end of loop over electrons

          rowTriplet[2]->Add(tv);
     
    } // end of loop over tracks
         
    prow[row+1]->Delete(); // delete tracks for this row
    delete prow[row+1];  // delete TObjArray for this row
    prow[row+1]=0; // set a pointer to NULL

  }

}  // end of MakeTriplet
//_____________________________________________________________________________
void AliTPC::ExB(Float_t *xyz)
{
  //------------------------------------------------
  //  ExB at the wires and wire number calulation
  //------------------------------------------------
  
  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------
   AliTPCParam * fTPCParam = &(fDigParam->GetParam());

   Float_t x1=xyz[0];
   fTPCParam->GetWire(x1);        //calculate nearest wire position
   Float_t dx=xyz[0]-x1;
   xyz[1]+=dx*fTPCParam->GetOmegaTau();

} // end of ExB
//_____________________________________________________________________________
void AliTPC::DigitizeRow(Int_t irow,Int_t isec,TObjArray **rowTriplet)
{
  //-----------------------------------------------------------
  // Single row digitization, coupling from the neighbouring
  // rows taken into account
  //-----------------------------------------------------------

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------
 

  AliTPCParam * fTPCParam = &(fDigParam->GetParam());
  Float_t chipgain= fTPCParam->GetChipGain();
  Float_t zerosup = fTPCParam->GetZeroSup();
  Int_t nrows =fTPCParam->GetNRow(isec);
  
  Int_t nTracks[3];
  Int_t n_of_pads[3];
  Int_t IndexRange[4];
  Int_t i1;
  Int_t iFlag; 

  //
  // iFlag = 0 -> inner row, iFlag = 1 -> the middle one, iFlag = 2 -> the outer one
  //

  nTracks[0]=rowTriplet[0]->GetEntries(); // lower row
  nTracks[1]=rowTriplet[1]->GetEntries(); // middle row
  nTracks[2]=rowTriplet[2]->GetEntries(); // upper row

    
  if(irow == 0){
    iFlag=0;
    n_of_pads[0]=fTPCParam->GetNPads(isec,0);
    n_of_pads[1]=fTPCParam->GetNPads(isec,1);
  }
  else if(irow == nrows-1){
     iFlag=2;
     n_of_pads[1]=fTPCParam->GetNPads(isec,irow-1);
     n_of_pads[2]=fTPCParam->GetNPads(isec,irow);
  }
  else {
    iFlag=1;
    for(i1=0;i1<3;i1++){
       n_of_pads[i1]=fTPCParam->GetNPads(isec,irow-1+i1);
    }
  }
 
  //
  //  Integrated signal for this row
  //  and a single track signal
  // 
   
  TMatrix *m1   = new TMatrix(1,n_of_pads[iFlag],1,MAXTPCTBK); // integrated
  TMatrix *m2   = new TMatrix(1,n_of_pads[iFlag],1,MAXTPCTBK); // single

  //

  TMatrix &Total  = *m1;

  //  Array of pointers to the label-signal list

  Int_t NofDigits = n_of_pads[iFlag]*MAXTPCTBK; // number of digits for this row

  Float_t  **pList = new Float_t* [NofDigits]; 

  Int_t lp;

  for(lp=0;lp<NofDigits;lp++)pList[lp]=0; // set all pointers to NULL

  //
  //  Straight signal and cross-talk, cross-talk is integrated over all
  //  tracks and added to the signal at the very end
  //
   

  for (i1=0;i1<nTracks[iFlag];i1++){

    m2->Zero(); // clear single track signal matrix
  
    Float_t TrackLabel = 
      GetSignal(rowTriplet[iFlag],i1,n_of_pads[iFlag],m2,m1,IndexRange); 

    GetList(TrackLabel,n_of_pads[iFlag],m2,IndexRange,pList);

  }

  // 
  //  Cross talk from the neighbouring pad-rows
  //

  TMatrix *m3 =  new TMatrix(1,n_of_pads[iFlag],1,MAXTPCTBK); // cross-talk

  TMatrix &Cross = *m3;

  if(iFlag == 0){

    // cross-talk from the outer row only (first pad row)

    GetCrossTalk(0,rowTriplet[1],nTracks[1],n_of_pads,m3);

  }
  else if(iFlag == 2){

    // cross-talk from the inner row only (last pad row)

    GetCrossTalk(2,rowTriplet[1],nTracks[1],n_of_pads,m3);

  }
  else{

    // cross-talk from both inner and outer rows

    GetCrossTalk(3,rowTriplet[0],nTracks[0],n_of_pads,m3); // inner
    GetCrossTalk(4,rowTriplet[2],nTracks[2],n_of_pads,m3); //outer
  }

  Total += Cross; // add the cross-talk

  //
  //  Convert analog signal to ADC counts
  //
   
  Int_t tracks[3];
  Int_t digits[5];


  for(Int_t ip=1;ip<n_of_pads[iFlag]+1;ip++){
    for(Int_t it=1;it<MAXTPCTBK+1;it++){

      Float_t q = Total(ip,it);

      Int_t gi =(it-1)*n_of_pads[iFlag]+ip-1; // global index

      q = gRandom->Gaus(q,fTPCParam->GetNoise()); // apply noise
      q *= (q_el*1.e15); // convert to fC
      q *= chipgain; // convert to mV   
      q *= (adc_sat/dyn_range); // convert to ADC counts  

      if(q <zerosup) continue; // do not fill zeros
      if(q > adc_sat) q = adc_sat;  // saturation

      //
      //  "real" signal or electronic noise (list = -1)?
      //    

      for(Int_t j1=0;j1<3;j1++){
        tracks[j1] = (pList[gi]) ?(Int_t)(*(pList[gi]+j1)) : -1;
      }

      digits[0]=isec;
      digits[1]=irow+1;
      digits[2]=ip;
      digits[3]=it;
      digits[4]= (Int_t)q;

      //  Add this digit

      AddDigit(tracks,digits);
    
    } // end of loop over time buckets
  }  // end of lop over pads 

  //
  //  This row has been digitized, delete nonused stuff
  //

  for(lp=0;lp<NofDigits;lp++){
    if(pList[lp]) delete [] pList[lp];
  }
  
  delete [] pList;

  delete m1;
  delete m2;
  delete m3;

} // end of DigitizeRow
//_____________________________________________________________________________
Float_t AliTPC::GetSignal(TObjArray *p1, Int_t ntr, Int_t np, TMatrix *m1, TMatrix *m2,
                          Int_t *IndexRange)
{

  //---------------------------------------------------------------
  //  Calculates 2-D signal (pad,time) for a single track,
  //  returns a pointer to the signal matrix and the track label 
  //  No digitization is performed at this level!!!
  //---------------------------------------------------------------

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------

  TVector *tv;
  AliTPCParam * fTPCParam = &(fDigParam->GetParam());
  AliTPCPRF2D * fPRF2D = &(fDigParam->GetPRF2D());
  AliTPCRF1D  * fRF    = &(fDigParam->GetRF()); 
 
  //to make the code faster we put parameters  to the stack

  Float_t zwidth  = fTPCParam->GetZWidth();
  Float_t zwidthm1  =1./zwidth;

  tv = (TVector*)p1->At(ntr); // pointer to a track
  TVector &v = *tv;
  
  Float_t label = v(0);

  Int_t CentralPad = (np+1)/2;
  Int_t PadNumber;
  Int_t nElectrons = (tv->GetNrows()-1)/4;
  Float_t range=((np-1)/2 + 0.5)*fTPCParam->GetPadPitchWidth(); // pad range
  range -= 0.5; // dead zone, 5mm from the edge, according to H.G. Fischer

  Float_t IneffFactor = 0.5; // inefficiency in the gain close to the edge, as above


  Float_t PadSignal[7]; // signal from a single electron

  TMatrix &signal = *m1;
  TMatrix &total = *m2;


  IndexRange[0]=9999; // min pad
  IndexRange[1]=-1; // max pad
  IndexRange[2]=9999; //min time
  IndexRange[3]=-1; // max time

  //
  //  Loop over all electrons
  //

  for(Int_t nel=0; nel<nElectrons; nel++){
   Int_t idx=nel*4;
   Float_t xwire = v(idx+1);
   Float_t y = v(idx+2);
   Float_t z = v(idx+3);


   Float_t absy=TMath::Abs(y);
   
   if(absy < 0.5*fTPCParam->GetPadPitchWidth()){
     PadNumber=CentralPad;
   }
   else if (absy < range){
     PadNumber=(Int_t) ((absy-0.5*fTPCParam->GetPadPitchWidth())/fTPCParam->GetPadPitchWidth() +1.);
     PadNumber=(Int_t) (TMath::Sign((Float_t)PadNumber, y)+CentralPad);
   }
   else continue; // electron out of pad-range , lost at the sector edge
    
   Float_t aval = (absy<range-0.5) ? v(idx+4):v(idx+4)*IneffFactor;
   

   Float_t dist = y - (Float_t)(PadNumber-CentralPad)*fTPCParam->GetPadPitchWidth();
   for (Int_t i=0;i<7;i++){
     PadSignal[i]=fPRF2D->GetPRF(dist+(i-3)*fTPCParam->GetPadPitchWidth(),xwire)*aval;
     PadSignal[i] *= fTPCParam->GetPadCoupling();
   }

   Int_t  LeftPad = TMath::Max(1,PadNumber-3);
   Int_t  RightPad = TMath::Min(np,PadNumber+3);

   Int_t pmin=LeftPad-PadNumber+3; // lower index of the pad_signal vector
   Int_t pmax=RightPad-PadNumber+3; // upper index     
   
   Float_t z_drift = (z_end-z)*zwidthm1;
   Float_t z_offset = z_drift-(Int_t)z_drift;
  //distance to the centre of nearest time bin (in time bin units)
   Int_t FirstBucket = (Int_t)z_drift+1; 


   // loop over time bins (4 bins is enough - 3 sigma truncated Gaussian)
   for (Int_t i2=0;i2<4;i2++){          
     Int_t TrueTime = FirstBucket+i2; // current time bucket
     Float_t dz   = (Float_t(i2)+z_offset)*zwidth; 
     Float_t ampl = fRF->GetRF(dz); 
     if( (TrueTime>MAXTPCTBK) ) break; // beyond the time range
     
     IndexRange[2]=TMath::Min(IndexRange[2],TrueTime); // min time
     IndexRange[3]=TMath::Max(IndexRange[3],TrueTime); // max time

     // loop over pads, from pmin to pmax
     for(Int_t i3=pmin;i3<=pmax;i3++){
       Int_t TruePad = LeftPad+i3-pmin;
       IndexRange[0]=TMath::Min(IndexRange[0],TruePad); // min pad
       IndexRange[1]=TMath::Max(IndexRange[1],TruePad); // max pad
       signal(TruePad,TrueTime)+=(PadSignal[i3]*ampl); // not converted to charge!!!
       total(TruePad,TrueTime)+=(PadSignal[i3]*ampl); // not converted to charge!!!
     } // end of pads loop
   }  // end of time loop
  } // end of loop over electrons

  return label; // returns track label when finished
}

//_____________________________________________________________________________
void AliTPC::GetList(Float_t label,Int_t np,TMatrix *m,Int_t *IndexRange,
                     Float_t **pList)
{
  //----------------------------------------------------------------------
  //  Updates the list of tracks contributing to digits for a given row
  //----------------------------------------------------------------------

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------

  TMatrix &signal = *m;

  // lop over nonzero digits

  for(Int_t it=IndexRange[2];it<IndexRange[3]+1;it++){
    for(Int_t ip=IndexRange[0];ip<IndexRange[1]+1;ip++){


        Int_t GlobalIndex = (it-1)*np+ip-1; // GlobalIndex starts from 0!
        
        if(!pList[GlobalIndex]){
        
          // 
	  // Create new list (6 elements - 3 signals and 3 labels),
	  // but only if the signal is > 0. 
	  //

          if(signal(ip,it)>0.){

          pList[GlobalIndex] = new Float_t [6];

	  // set list to -1 

          *pList[GlobalIndex] = -1.;
          *(pList[GlobalIndex]+1) = -1.;
          *(pList[GlobalIndex]+2) = -1.;
          *(pList[GlobalIndex]+3) = -1.;
          *(pList[GlobalIndex]+4) = -1.;
          *(pList[GlobalIndex]+5) = -1.;


          *pList[GlobalIndex] = label;
          *(pList[GlobalIndex]+3) = signal(ip,it);}
        }
        else{

	  // check the signal magnitude

          Float_t highest = *(pList[GlobalIndex]+3);
          Float_t middle = *(pList[GlobalIndex]+4);
          Float_t lowest = *(pList[GlobalIndex]+5);

	  //
	  //  compare the new signal with already existing list
	  //

          if(signal(ip,it)<lowest) continue; // neglect this track

	  //

          if (signal(ip,it)>highest){
            *(pList[GlobalIndex]+5) = middle;
            *(pList[GlobalIndex]+4) = highest;
            *(pList[GlobalIndex]+3) = signal(ip,it);

            *(pList[GlobalIndex]+2) = *(pList[GlobalIndex]+1);
            *(pList[GlobalIndex]+1) = *pList[GlobalIndex];
            *pList[GlobalIndex] = label;
	  }
          else if (signal(ip,it)>middle){
            *(pList[GlobalIndex]+5) = middle;
            *(pList[GlobalIndex]+4) = signal(ip,it);

            *(pList[GlobalIndex]+2) = *(pList[GlobalIndex]+1);
            *(pList[GlobalIndex]+1) = label;
	  }
          else{
            *(pList[GlobalIndex]+5) = signal(ip,it);
            *(pList[GlobalIndex]+2) = label;
	  }
        }

    } // end of loop over pads
  } // end of loop over time bins




}//end of GetList
//___________________________________________________________________
void AliTPC::MakeSector(Int_t isec,Int_t nrows,TTree *TH,
                        Stat_t ntracks,TObjArray **row)
{

  //-----------------------------------------------------------------
  // Prepares the sector digitization, creates the vectors of
  // tracks for each row of this sector. The track vector
  // contains the track label and the position of electrons.
  //-----------------------------------------------------------------

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------

  AliTPCParam * fTPCParam = &(fDigParam->GetParam());
  Int_t i;
  Float_t xyz[3]; 

  AliTPChit *tpcHit; // pointer to a sigle TPC hit    
 
  //----------------------------------------------
  // Create TObjArray-s, one for each row,
  // each TObjArray will store the TVectors
  // of electrons, one TVector per each track.
  //---------------------------------------------- 
    
  for(i=0; i<nrows; i++){
    row[i] = new TObjArray;
  }
  Int_t *n_of_electrons = new Int_t [nrows]; // electron counter for each row
  TVector **tr = new TVector* [nrows]; //pointers to the track vectors

  //--------------------------------------------------------------------
  //  Loop over tracks, the "track" contains the full history
  //--------------------------------------------------------------------

  Int_t previousTrack,currentTrack;
  previousTrack = -1; // nothing to store so far!

  for(Int_t track=0;track<ntracks;track++){

    ResetHits();

    TH->GetEvent(track); // get next track
    Int_t nhits = fHits->GetEntriesFast(); // get number of hits for this track

    if(nhits == 0) continue; // no hits in the TPC for this track

    //--------------------------------------------------------------
    //  Loop over hits
    //--------------------------------------------------------------

    for(Int_t hit=0;hit<nhits;hit++){

      tpcHit = (AliTPChit*)fHits->UncheckedAt(hit); // get a pointer to a hit
      
      Int_t sector=tpcHit->fSector; // sector number
      if(sector != isec) continue; //terminate iteration

	currentTrack = tpcHit->fTrack; // track number
        if(currentTrack != previousTrack){
                          
           // store already filled fTrack
              
	   for(i=0;i<nrows;i++){
             if(previousTrack != -1){
	       if(n_of_electrons[i]>0){
	         TVector &v = *tr[i];
		 v(0) = previousTrack;
                 tr[i]->ResizeTo(3*n_of_electrons[i]+1); // shrink if necessary
	         row[i]->Add(tr[i]);                     
	       }
               else{
                 delete tr[i]; // delete empty TVector
                 tr[i]=0;
	       }
	     }

             n_of_electrons[i]=0;
             tr[i] = new TVector(361); // TVectors for the next fTrack

	   } // end of loop over rows
	       
           previousTrack=currentTrack; // update track label 
	}
	   
	Int_t QI = (Int_t) (tpcHit->fQ); // energy loss (number of electrons)

       //---------------------------------------------------
       //  Calculate the electron attachment probability
       //---------------------------------------------------

        Float_t time = 1.e6*(z_end-TMath::Abs(tpcHit->fZ))/fTPCParam->GetDriftV(); 
	// in microseconds!	
	Float_t AttProb = fTPCParam->GetAttCoef()*
	  fTPCParam->GetOxyCont()*time; //  fraction! 
   
	//-----------------------------------------------
	//  Loop over electrons
	//-----------------------------------------------

        for(Int_t nel=0;nel<QI;nel++){
          // skip if electron lost due to the attachment
          if((gRandom->Rndm(0)) < AttProb) continue; // electron lost!
	  xyz[0]=tpcHit->fX;
	  xyz[1]=tpcHit->fY;
	  xyz[2]=tpcHit->fZ;
	  ElDiff(xyz); // Appply the diffusion
	  Int_t row_number;
	  fTPCParam->XYZtoCRXYZ(xyz,isec,row_number,3);

	  //transform position to local coordinates
	  //option 3 means that we calculate x position relative to 
	  //nearest pad row 

	  if ((row_number<0)||row_number>=fTPCParam->GetNRow(isec)) continue;
	  ExB(xyz); // ExB effect at the sense wires
	  n_of_electrons[row_number]++;	  
	  //----------------------------------
	  // Expand vector if necessary
	  //----------------------------------
	  if(n_of_electrons[row_number]>120){
	    Int_t range = tr[row_number]->GetNrows();
	    if(n_of_electrons[row_number] > (range-1)/3){
	      tr[row_number]->ResizeTo(range+150); // Add 50 electrons
	    }
	  }
	  
	  TVector &v = *tr[row_number];
	  Int_t idx = 3*n_of_electrons[row_number]-2;

	  v(idx)=  xyz[0];   // X
	  v(idx+1)=xyz[1];   // Y (along the pad-row)
          v(idx+2)=xyz[2];   // Z
	    
	} // end of loop over electrons
        
      } // end of loop over hits
    } // end of loop over tracks

    //
    //   store remaining track (the last one) if not empty
    //

     for(i=0;i<nrows;i++){
       if(n_of_electrons[i]>0){
          TVector &v = *tr[i];
	  v(0) = previousTrack;
          tr[i]->ResizeTo(3*n_of_electrons[i]+1); // shrink if necessary
	  row[i]->Add(tr[i]);  
	}
	else{
          delete tr[i];
          tr[i]=0;
	}  
      }  

          delete [] tr;
          delete [] n_of_electrons; 

} // end of MakeSector
//_____________________________________________________________________________
void AliTPC::GetCrossTalk (Int_t iFlag,TObjArray *p,Int_t ntracks,Int_t *npads,
                           TMatrix *m)
{

  //-------------------------------------------------------------
  // Calculates the cross-talk from one row (inner or outer)
  //-------------------------------------------------------------

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------

  //
  // iFlag=2 & 3 --> cross-talk from the inner row
  // iFlag=0 & 4 --> cross-talk from the outer row
  //
  AliTPCParam * fTPCParam = &(fDigParam->GetParam());
  AliTPCPRF2D * fPRF2D = &(fDigParam->GetPRF2D());
  AliTPCRF1D  * fRF    = &(fDigParam->GetRF()); 
 
  //to make code faster

  Float_t zwidth  = fTPCParam->GetZWidth();
  Float_t zwidthm1  =1/fTPCParam->GetZWidth();

 Int_t PadNumber;
 Float_t xwire;

 Int_t nPadsSignal; // for this pads the signal is calculated
 Float_t range;     // sense wire range
 Int_t nPadsDiff;

 Float_t IneffFactor=0.5; // gain inefficiency close to the edges

 if(iFlag == 0){   
   // 1-->0
   nPadsSignal = npads[1];
   range = ((npads[1]-1)/2 + 0.5)*fTPCParam->GetPadPitchWidth();  
   nPadsDiff = (npads[1]-npads[0])/2;
 }  
 else if (iFlag == 2){
   // 1-->2
   nPadsSignal = npads[2];
   range = ((npads[1]-1)/2 + 0.5)*fTPCParam->GetPadPitchWidth();
   nPadsDiff = 0;
 }
 else if (iFlag == 3){
   // 0-->1
   nPadsSignal = npads[1];
   range = ((npads[0]-1)/2 + 0.5)*fTPCParam->GetPadPitchWidth();
   nPadsDiff = 0;
 }
 else{
   // 2-->1
   nPadsSignal = npads[2];
   range = ((npads[2]-1)/2 + 0.5)*fTPCParam->GetPadPitchWidth();
   nPadsDiff = (npads[2]-npads[1])/2;
 }

 range-=0.5; // dead zone close to the edges

 TVector *tv; 
 TMatrix &signal = *m;

 Int_t CentralPad = (nPadsSignal+1)/2;
 Float_t PadSignal[7]; // signal from a single electron
 // Loop over tracks
 for(Int_t nt=0;nt<ntracks;nt++){
   tv=(TVector*)p->At(nt); // pointer to a track
   TVector &v = *tv;
   Int_t nElectrons = (tv->GetNrows()-1)/4;
   // Loop over electrons
   for(Int_t nel=0; nel<nElectrons; nel++){
     Int_t idx=nel*4;
     xwire=v(idx+1);
 
     if (iFlag==0) xwire+=fTPCParam->GetPadPitchLength();
     if (iFlag==2)  xwire-=fTPCParam->GetPadPitchLength();
     if (iFlag==3)  xwire-=fTPCParam->GetPadPitchLength();
     if (iFlag==4)  xwire+=fTPCParam->GetPadPitchLength();  
   
     //  electron acceptance for the cross-talk (only the closest wire)  

     Float_t dxMax = fTPCParam->GetPadPitchLength()*0.5+fTPCParam->GetWWPitch();
     if(TMath::Abs(xwire)>dxMax) continue;

     Float_t y = v(idx+2);
     Float_t z = v(idx+3);
     Float_t absy=TMath::Abs(y);

     if(absy < 0.5*fTPCParam->GetPadPitchWidth()){
       PadNumber=CentralPad;
     }
     else if (absy < range){
       PadNumber=(Int_t) ((absy-0.5*fTPCParam->GetPadPitchWidth())/fTPCParam->GetPadPitchWidth() +1.);
       PadNumber=(Int_t) (TMath::Sign((Float_t)PadNumber, y)+CentralPad);
     }
     else continue; // electron out of sense wire range, lost at the sector edge

     Float_t aval = (absy<range-0.5) ? v(idx+4):v(idx+4)*IneffFactor;

     Float_t dist = y - (Float_t)(PadNumber-CentralPad)*fTPCParam->GetPadPitchWidth();
       
     for (Int_t i=0;i<7;i++){
       PadSignal[i]=fPRF2D->GetPRF(dist+(3-i)*fTPCParam->GetPadPitchWidth(),xwire)*aval;

       PadSignal[i] *= fTPCParam->GetPadCoupling();
     }
     // real pad range

     Int_t  LeftPad = TMath::Max(1,PadNumber-3);
     Int_t  RightPad = TMath::Min(nPadsSignal,PadNumber+3);

     Int_t pmin=LeftPad-PadNumber+3; // lower index of the pad_signal vector
     Int_t pmax=RightPad-PadNumber+3; // upper index  


     Float_t z_drift = (z_end-z)*zwidthm1;
     Float_t z_offset = z_drift-(Int_t)z_drift;
     //distance to the centre of nearest time bin (in time bin units)
     Int_t FirstBucket = (Int_t)z_drift+1; 
     // MI check it --time offset
     for (Int_t i2=0;i2<4;i2++){     
       Int_t TrueTime = FirstBucket+i2; // current time bucket
       Float_t dz   = (Float_t(i2)+z_offset)*zwidth; 
       Float_t ampl = fRF->GetRF(dz); 
       if((TrueTime>MAXTPCTBK)) break; // beyond the time range


       // loop over pads, from pmin to pmax

       for(Int_t i3=pmin;i3<pmax+1;i3++){
         Int_t TruePad = LeftPad+i3-pmin;

         if(TruePad<nPadsDiff+1 || TruePad > nPadsSignal-nPadsDiff) continue;

         TruePad -= nPadsDiff;
         signal(TruePad,TrueTime)+=(PadSignal[i3]*ampl); // not converted to charge!

       } // end of loop over pads
     } // end of loop over time bins

   } // end of loop over electrons 

 } // end of loop over tracks

} // end of GetCrossTalk



//_____________________________________________________________________________
void AliTPC::Init()
{
  //
  // Initialise TPC detector after definition of geometry
  //
  Int_t i;
  //
  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" TPC_INIT ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");
  //
  for(i=0;i<80;i++) printf("*");
  printf("\n");
}

//_____________________________________________________________________________
void AliTPC::MakeBranch(Option_t* option)
{
  //
  // Create Tree branches for the TPC.
  //
  Int_t buffersize = 4000;
  char branchname[10];
  sprintf(branchname,"%s",GetName());

  AliDetector::MakeBranch(option);

  char *D = strstr(option,"D");

  if (fDigits   && gAlice->TreeD() && D) {
    gAlice->TreeD()->Branch(branchname,&fDigits, buffersize);
    printf("Making Branch %s for digits\n",branchname);
  }	

  char *R = strstr(option,"R");

  if (fClusters && gAlice->TreeR() && R) {
    gAlice->TreeR()->Branch(branchname,&fClusters, buffersize);
    printf("Making Branch %s for Clusters\n",branchname);
  }	
}
 
//_____________________________________________________________________________
void AliTPC::ResetDigits()
{
  //
  // Reset number of digits and the digits array for this detector
  // reset clusters
  //
  fNdigits   = 0;
  //  if (fDigits)   fDigits->Clear();
  if (fDigParam->GetArray()!=0)  fDigParam->GetArray()->Clear();
  fNclusters = 0;
  if (fClusters) fClusters->Clear();
}

//_____________________________________________________________________________
void AliTPC::SetSecAL(Int_t sec)
{
  //---------------------------------------------------
  // Activate/deactivate selection for lower sectors
  //---------------------------------------------------

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------

  fSecAL = sec;
}

//_____________________________________________________________________________
void AliTPC::SetSecAU(Int_t sec)
{
  //----------------------------------------------------
  // Activate/deactivate selection for upper sectors
  //---------------------------------------------------

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------

  fSecAU = sec;
}

//_____________________________________________________________________________
void AliTPC::SetSecLows(Int_t s1,Int_t s2,Int_t s3,Int_t s4,Int_t s5, Int_t s6)
{
  //----------------------------------------
  // Select active lower sectors
  //----------------------------------------

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------

  fSecLows[0] = s1;
  fSecLows[1] = s2;
  fSecLows[2] = s3;
  fSecLows[3] = s4;
  fSecLows[4] = s5;
  fSecLows[5] = s6;
}

//_____________________________________________________________________________
void AliTPC::SetSecUps(Int_t s1,Int_t s2,Int_t s3,Int_t s4,Int_t s5, Int_t s6,
                       Int_t s7, Int_t s8 ,Int_t s9 ,Int_t s10, 
                       Int_t s11 , Int_t s12)
{
  //--------------------------------
  // Select active upper sectors
  //--------------------------------

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------

  fSecUps[0] = s1;
  fSecUps[1] = s2;
  fSecUps[2] = s3;
  fSecUps[3] = s4;
  fSecUps[4] = s5;
  fSecUps[5] = s6;
  fSecUps[6] = s7;
  fSecUps[7] = s8;
  fSecUps[8] = s9;
  fSecUps[9] = s10;
  fSecUps[10] = s11;
  fSecUps[11] = s12;
}

//_____________________________________________________________________________
void AliTPC::SetSens(Int_t sens)
{

  //-------------------------------------------------------------
  // Activates/deactivates the sensitive strips at the center of
  // the pad row -- this is for the space-point resolution calculations
  //-------------------------------------------------------------

  //-----------------------------------------------------------------
  // Origin: Marek Kowalski  IFJ, Krakow, Marek.Kowalski@ifj.edu.pl
  //-----------------------------------------------------------------

  fSens = sens;
}

//_____________________________________________________________________________
void AliTPC::Streamer(TBuffer &R__b)
{
  //
  // Stream an object of class AliTPC.
  //
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      AliDetector::Streamer(R__b);
      if (R__v < 2) return;
      R__b >> fNsectors;
      R__b >> fNclusters;
      R__b >> fNtracks;
      fClustersIndex = new Int_t[fNsectors+1];
      fDigitsIndex   = new Int_t[fNsectors+1];
   } else {
      R__b.WriteVersion(AliTPC::IsA());
      AliDetector::Streamer(R__b);
      R__b << fNsectors;
      R__b << fNclusters;
      R__b << fNtracks;
   }
}
  
ClassImp(AliTPCcluster)
 
//_____________________________________________________________________________
AliTPCcluster::AliTPCcluster(Float_t *hits, Int_t *lab)
{
  //
  // Creates a simulated cluster for the TPC
  //
  fTracks[0]  = lab[0];
  fTracks[1]  = lab[1];
  fTracks[2]  = lab[2];
  fSector     = lab[3];
  fPadRow     = lab[4];
  fY          = hits[0];
  fZ          = hits[1];
  fQ          = hits[2];
  fSigmaY2    = hits[3];
  fSigmaZ2    = hits[4];
}
 
//_____________________________________________________________________________
void AliTPCcluster::GetXYZ(Float_t *x, const AliTPCParam *par) const 
{
  //
  // Transformation from local to global coordinate system
  //
  x[0]=par->GetPadRowRadii(fSector,fPadRow-1);
  x[1]=fY;
  x[2]=fZ;
  par->CRXYZtoXYZ(x,fSector,fPadRow-1,1);
}
 
//_____________________________________________________________________________
Int_t AliTPCcluster::Compare(TObject * o)
{
  //
  // compare two clusters according y coordinata
  AliTPCcluster *cl= (AliTPCcluster *)o;
  if (fY<cl->fY) return -1;
  if (fY==cl->fY) return 0;
  return 1;  
}

Bool_t AliTPCcluster::IsSortable() const
{
  //
  //make AliTPCcluster sortabale
  return kTRUE; 
}



ClassImp(AliTPCdigit)
 
//_____________________________________________________________________________
AliTPCdigit::AliTPCdigit(Int_t *tracks, Int_t *digits):
  AliDigit(tracks)
{
  //
  // Creates a TPC digit object
  //
  fSector     = digits[0];
  fPadRow     = digits[1];
  fPad        = digits[2];
  fTime       = digits[3];
  fSignal     = digits[4];
}

 
ClassImp(AliTPChit)
 
//_____________________________________________________________________________
AliTPChit::AliTPChit(Int_t shunt, Int_t track, Int_t *vol, Float_t *hits):
AliHit(shunt,track)
{
  //
  // Creates a TPC hit object
  //
  fSector     = vol[0];
  fPadRow     = vol[1];
  fX          = hits[0];
  fY          = hits[1];
  fZ          = hits[2];
  fQ          = hits[3];
}
 
 
ClassImp(AliTPCtrack)
 
//_____________________________________________________________________________
AliTPCtrack::AliTPCtrack(Float_t *hits)
{
  //
  // Default creator for a TPC reconstructed track object
  //
  ref=hits[0]; // This is dummy code !
}

AliTPCtrack::AliTPCtrack(const AliTPCcluster& c,const TVector& xx,
			 const TMatrix& CC, const AliTPCParam *p):
  x(xx),C(CC),clusters(MAX_CLUSTER)
{
  //
  // Standard creator for a TPC reconstructed track object
  //
  chi2=0.;
  int sec=c.fSector-1, row=c.fPadRow-1;
  ref=p->GetPadRowRadii(sec+1,row);

  if (sec<24) { 
    fAlpha=(sec%12)*alpha_low;
  } else { 
    fAlpha=(sec%24)*alpha_up;
  }
  clusters.AddLast((AliTPCcluster*)(&c));
}

//_____________________________________________________________________________
AliTPCtrack::AliTPCtrack(const AliTPCtrack& t) : x(t.x), C(t.C),
  clusters(t.clusters.GetEntriesFast()) 
{
  //
  // Copy creator for a TPC reconstructed track
  //
  ref=t.ref;
  chi2=t.chi2;
  fAlpha=t.fAlpha;
  int n=t.clusters.GetEntriesFast();
  for (int i=0; i<n; i++) clusters.AddLast(t.clusters.UncheckedAt(i));
}

//_____________________________________________________________________________
Double_t AliTPCtrack::GetY(Double_t xk) const 
{
  //
  //
  //
  Double_t c2=x(2)*xk - x(3);
  if (TMath::Abs(c2) >= 0.999) {
    if (*this>10) cerr<<*this<<" AliTPCtrack warning: No y for this x !\n";
    return 0.;
  }
  Double_t c1=x(2)*ref - x(3);
  Double_t r1=sqrt(1.-c1*c1), r2=sqrt(1.-c2*c2);
  Double_t dx=xk-ref;
  return x(0) + dx*(c1+c2)/(r1+r2);
}

//_____________________________________________________________________________
int AliTPCtrack::PropagateTo(Double_t xk,Double_t x0,Double_t rho,Double_t pm)
{
  //
  // Propagate a TPC reconstructed track
  //
  if (TMath::Abs(x(2)*xk - x(3)) >= 0.999) {
    if (*this>3) cerr<<*this<<" AliTPCtrack warning: Propagation failed !\n";
    return 0;
  }

  Double_t x1=ref, x2=x1+0.5*(xk-x1), dx=x2-x1, y1=x(0), z1=x(1);
  Double_t c1=x(2)*x1 - x(3), r1=sqrt(1.- c1*c1);
  Double_t c2=x(2)*x2 - x(3), r2=sqrt(1.- c2*c2);
  
  x(0) += dx*(c1+c2)/(r1+r2);
  x(1) += dx*(c1+c2)/(c1*r2 + c2*r1)*x(4);

  TMatrix F(5,5); F.UnitMatrix();
  Double_t rr=r1+r2, cc=c1+c2, xx=x1+x2;
  F(0,2)= dx*(rr*xx + cc*(c1*x1/r1+c2*x2/r2))/(rr*rr);
  F(0,3)=-dx*(2*rr + cc*(c1/r1 + c2/r2))/(rr*rr);
  Double_t cr=c1*r2+c2*r1;
  F(1,2)= dx*x(4)*(cr*xx-cc*(r1*x2-c2*c1*x1/r1+r2*x1-c1*c2*x2/r2))/(cr*cr);
  F(1,3)=-dx*x(4)*(2*cr + cc*(c2*c1/r1-r1 + c1*c2/r2-r2))/(cr*cr);
  F(1,4)= dx*cc/cr; 
  TMatrix tmp(F,TMatrix::kMult,C);
  C.Mult(tmp,TMatrix(TMatrix::kTransposed,F));
  
  ref=x2;
  
  //Multiple scattering******************
  Double_t ey=x(2)*ref - x(3);
  Double_t ex=sqrt(1-ey*ey);
  Double_t ez=x(4);
  TMatrix Q(5,5); Q=0.;
  Q(2,2)=ez*ez+ey*ey;   Q(2,3)=-ex*ey;       Q(2,4)=-ex*ez;
  Q(3,2)=Q(2,3);        Q(3,3)= ez*ez+ex*ex; Q(3,4)=-ey*ez;
  Q(4,2)=Q(2,4);        Q(4,3)= Q(3,4);      Q(4,4)=1.;
  
  F=0;
  F(2,2)=-x(2)*ex;          F(2,3)=-x(2)*ey;
  F(3,2)=-ex*(x(2)*ref-ey); F(3,3)=-(1.+ x(2)*ref*ey - ey*ey);
  F(4,2)=-ez*ex;            F(4,3)=-ez*ey;           F(4,4)=1.;
  
  tmp.Mult(F,Q);
  Q.Mult(tmp,TMatrix(TMatrix::kTransposed,F));
  
  Double_t p2=GetPt()*GetPt()*(1.+x(4)*x(4));
  Double_t beta2=p2/(p2 + pm*pm);
  Double_t d=sqrt((x1-ref)*(x1-ref)+(y1-x(0))*(y1-x(0))+(z1-x(1))*(z1-x(1)));
  d*=2.;
  Double_t theta2=14.1*14.1/(beta2*p2*1e6)*d/x0*rho;
  Q*=theta2;
  C+=Q;
  
  //Energy losses************************
  Double_t dE=0.153e-3/beta2*(log(5940*beta2/(1-beta2)) - beta2)*d*rho;
  if (x1 < x2) dE=-dE;
  x(2)*=(1.- sqrt(p2+pm*pm)/p2*dE);
  //x(3)*=(1.- sqrt(p2+pm*pm)/p2*dE);
  
  x1=ref; x2=xk; y1=x(0); z1=x(1);
  c1=x(2)*x1 - x(3); r1=sqrt(1.- c1*c1);
  c2=x(2)*x2 - x(3); r2=sqrt(1.- c2*c2);
  
  x(0) += dx*(c1+c2)/(r1+r2);
  x(1) += dx*(c1+c2)/(c1*r2 + c2*r1)*x(4);
  
  F.UnitMatrix();
  rr=r1+r2; cc=c1+c2; xx=x1+x2;
  F(0,2)= dx*(rr*xx + cc*(c1*x1/r1+c2*x2/r2))/(rr*rr);
  F(0,3)=-dx*(2*rr + cc*(c1/r1 + c2/r2))/(rr*rr);
  cr=c1*r2+c2*r1;
  F(1,2)= dx*x(4)*(cr*xx-cc*(r1*x2-c2*c1*x1/r1+r2*x1-c1*c2*x2/r2))/(cr*cr);
  F(1,3)=-dx*x(4)*(2*cr + cc*(c2*c1/r1-r1 + c1*c2/r2-r2))/(cr*cr);
  F(1,4)= dx*cc/cr; 
  tmp.Mult(F,C);
  C.Mult(tmp,TMatrix(TMatrix::kTransposed,F));
  
  ref=x2;
  
  return 1;
}

//_____________________________________________________________________________
void AliTPCtrack::PropagateToVertex(Double_t x0,Double_t rho,Double_t pm) 
{
  //
  // Propagate a reconstructed track from the vertex
  //
  Double_t c=x(2)*ref - x(3);
  Double_t tgf=-x(3)/(x(2)*x(0) + sqrt(1-c*c));
  Double_t snf=tgf/sqrt(1.+ tgf*tgf);
  Double_t xv=(x(3)+snf)/x(2);
  PropagateTo(xv,x0,rho,pm);
}

//_____________________________________________________________________________
void AliTPCtrack::Update(const AliTPCcluster *c, Double_t chisq)
{
  //
  // Update statistics for a reconstructed TPC track
  //
  TMatrix H(2,5); H.UnitMatrix();
  TMatrix Ht(TMatrix::kTransposed,H);
  TVector m(2);   m(0)=c->fY; m(1)=c->fZ; 
  TMatrix V(2,2); V(0,0)=c->fSigmaY2; V(0,1)=0.; V(1,0)=0.; V(1,1)=c->fSigmaZ2;

  TMatrix tmp(H,TMatrix::kMult,C);
  TMatrix R(tmp,TMatrix::kMult,Ht); R+=V;
  
  Double_t det=(Double_t)R(0,0)*R(1,1) - (Double_t)R(0,1)*R(1,0);
  R(0,1)=R(0,0); R(0,0)=R(1,1); R(1,1)=R(0,1); 
  R(1,0)*=-1; R(0,1)=R(1,0);
  R*=1./det;
  
  //R.Invert();
  
  TMatrix K(C,TMatrix::kMult,Ht); K*=R;
  
  TVector savex=x;
  x*=H; x-=m; x*=-1; x*=K; x+=savex;
  if (TMath::Abs(x(2)*ref-x(3)) >= 0.999) {
    if (*this>3) cerr<<*this<<" AliTPCtrack warning: Filtering failed !\n";
    x=savex;
    return;
  }
  
  TMatrix saveC=C;
  C.Mult(K,tmp); C-=saveC; C*=-1;
  
  clusters.AddLast((AliTPCcluster*)c);
  chi2 += chisq;
}

//_____________________________________________________________________________
int AliTPCtrack::Rotate(Double_t alpha)
{
  //
  // Rotate a reconstructed TPC track
  //
  fAlpha += alpha;
  
  Double_t x1=ref, y1=x(0);
  Double_t ca=cos(alpha), sa=sin(alpha);
  Double_t r1=x(2)*ref - x(3);
  
  ref = x1*ca + y1*sa;
  x(0)=-x1*sa + y1*ca;
  x(3)=x(3)*ca + (x(2)*y1 + sqrt(1.- r1*r1))*sa;
  
  Double_t r2=x(2)*ref - x(3);
  if (TMath::Abs(r2) >= 0.999) {
    if (*this>3) cerr<<*this<<" AliTPCtrack warning: Rotation failed !\n";
    return 0;
  }
  
  Double_t y0=x(0) + sqrt(1.- r2*r2)/x(2);
  if ((x(0)-y0)*x(2) >= 0.) {
    if (*this>3) cerr<<*this<<" AliTPCtrack warning: Rotation failed !!!\n";
    return 0;
  }
  
  TMatrix F(5,5); F.UnitMatrix();
  F(0,0)=ca;
  F(3,0)=x(2)*sa; 
  F(3,2)=(y1 - r1*x1/sqrt(1.- r1*r1))*sa; 
  F(3,3)= ca + sa*r1/sqrt(1.- r1*r1);
  TMatrix tmp(F,TMatrix::kMult,C); 
  // Double_t dy2=C(0,0);
  C.Mult(tmp,TMatrix(TMatrix::kTransposed,F));
  // C(0,0)+=dy2*sa*sa*r1*r1/(1.- r1*r1);
  // C(1,1)+=dy2*sa*sa*x(4)*x(4)/(1.- r1*r1);
  
  return 1;
}

//_____________________________________________________________________________
void AliTPCtrack::UseClusters() const 
{
  //
  //
  //
  int num_of_clusters=clusters.GetEntriesFast();
  for (int i=0; i<num_of_clusters; i++) {
    //if (i<=14) continue;
    AliTPCcluster *c=(AliTPCcluster*)clusters.UncheckedAt(i);
    c->Use();   
  }
}

//_____________________________________________________________________________
Double_t AliTPCtrack::GetPredictedChi2(const AliTPCcluster *c) const 
{
  //
  // Calculate chi2 for a reconstructed TPC track
  //
  TMatrix H(2,5); H.UnitMatrix();
  TVector m(2);   m(0)=c->fY; m(1)=c->fZ; 
  TMatrix V(2,2); V(0,0)=c->fSigmaY2; V(0,1)=0.; V(1,0)=0.; V(1,1)=c->fSigmaZ2;
  TVector res=x;  res*=H; res-=m; //res*=-1; 
  TMatrix tmp(H,TMatrix::kMult,C);
  TMatrix R(tmp,TMatrix::kMult,TMatrix(TMatrix::kTransposed,H)); R+=V;
  
  Double_t det=(Double_t)R(0,0)*R(1,1) - (Double_t)R(0,1)*R(1,0);
  if (TMath::Abs(det) < 1.e-10) {
    if (*this>3) cerr<<*this<<" AliTPCtrack warning: Singular matrix !\n";
    return 1e10;
  }
  R(0,1)=R(0,0); R(0,0)=R(1,1); R(1,1)=R(0,1); 
  R(1,0)*=-1; R(0,1)=R(1,0);
  R*=1./det;
  
  //R.Invert();
  
  TVector r=res;
  res*=R;
  return r*res;
}

//_____________________________________________________________________________
int AliTPCtrack::GetLab() const 
{
  //
  //
  //
  int lab=123456789;
  struct {
    int lab;
    int max;
  } s[MAX_CLUSTER]={{0,0}};
  
  int i;
  int num_of_clusters=clusters.GetEntriesFast();
  for (i=0; i<num_of_clusters; i++) {
    AliTPCcluster *c=(AliTPCcluster*)clusters.UncheckedAt(i);
    lab=TMath::Abs(c->fTracks[0]);
    int j;
    for (j=0; j<MAX_CLUSTER; j++)
      if (s[j].lab==lab || s[j].max==0) break;
    s[j].lab=lab;
    s[j].max++;
  }
  
  int max=0;
  for (i=0; i<num_of_clusters; i++) 
    if (s[i].max>max) {max=s[i].max; lab=s[i].lab;}

  for (i=0; i<num_of_clusters; i++) {
    AliTPCcluster *c=(AliTPCcluster*)clusters.UncheckedAt(i);
    if (TMath::Abs(c->fTracks[1]) == lab ||
        TMath::Abs(c->fTracks[2]) == lab ) max++;
  }
  
  if (1.-float(max)/num_of_clusters > 0.10) return -lab;
  
  if (num_of_clusters < 6) return lab;
  
  max=0;
  for (i=1; i<=6; i++) {
    AliTPCcluster *c=(AliTPCcluster*)clusters.UncheckedAt(num_of_clusters-i);
    if (lab == TMath::Abs(c->fTracks[0]) ||
        lab == TMath::Abs(c->fTracks[1]) ||
        lab == TMath::Abs(c->fTracks[2])) max++;
  }
  if (max<3) return -lab;
  
  return lab;
}

//_____________________________________________________________________________
void AliTPCtrack::GetPxPyPz(Double_t& px, Double_t& py, Double_t& pz) const 
{
  //
  // Get reconstructed TPC track momentum
  //
  Double_t pt=0.3*FIELD/TMath::Abs(x(2))/100; // GeV/c
  Double_t r=x(2)*ref-x(3);
  Double_t y0=x(0) + sqrt(1.- r*r)/x(2);
  px=-pt*(x(0)-y0)*x(2);    //cos(phi);
  py=-pt*(x(3)-ref*x(2));   //sin(phi);
  pz=pt*x(4);
  Double_t tmp=px*TMath::Cos(fAlpha) - py*TMath::Sin(fAlpha);
  py=px*TMath::Sin(fAlpha) + py*TMath::Cos(fAlpha);
  px=tmp;  
}

//_____________________________________________________________________________
//
//     Classes for internal tracking use
//

//_____________________________________________________________________________
void AliTPCRow::InsertCluster(const AliTPCcluster* c) 
{
  //
  // Insert a cluster in the list
  //
  if (num_of_clusters==MAX_CLUSTER_PER_ROW) {
    cerr<<"AliTPCRow::InsertCluster(): Too many clusters !\n"; return;
  }
  if (num_of_clusters==0) {clusters[num_of_clusters++]=c; return;}
  int i=Find(c->fY);
  memmove(clusters+i+1 ,clusters+i,(num_of_clusters-i)*sizeof(AliTPCcluster*));
  clusters[i]=c; num_of_clusters++;
}

//_____________________________________________________________________________
int AliTPCRow::Find(Double_t y) const 
{
  //
  //
  //
  if (y <= clusters[0]->fY) return 0;
  if (y > clusters[num_of_clusters-1]->fY) return num_of_clusters;
  int b=0, e=num_of_clusters-1, m=(b+e)/2;
  for (; b<e; m=(b+e)/2) {
    if (y > clusters[m]->fY) b=m+1;
    else e=m; 
  }
  return m;
}

