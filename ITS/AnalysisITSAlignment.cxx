//////////////////////////////////////////////////////////////////////////
//  Alice ITS first detector alignment program.                         //
//                                                                      //
// version: 0.0.0 Draft.                                                //
// Date: April 18 1999                                                  //
// By: Bjorn S. Nilsen                                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include <fstream.h>
#include <stdio.h>
#include <time.h>
#include "AliITS.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TProfile.h"
#include "TH2.h"
#include "TArray.h"
#include "TRandom.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TParticle.h"
#include "AliRun.h"
#include "AliITSgeom.h"
#include "AnalysisITSAlignment.h"
#include "AliITSstatistics2.h"
#include "AliITSstatistics.h"
#include "AliITSAlignmentTrack.h"
#include "AliITSAlignmentModule.h"

#define PRIMARYONLY 0  // if PIMARY_ONLY ==0 then all tracks
#define MOMENTUMCUT 1.0 // Sets the momentum cut for the tracks.

//Double_t resz[]={0.00680,0.00650,0.00265,0.00265,0.0830,0.0830};
//Double_t resx[]={0.00105,0.00115,0.00340,0.00340,0.0020,0.0020};
// From table 1.3 of ITS TDR. Spatial precisions and Cell size
Double_t  resz[]={0.00700,0.00700,0.00280,0.00280,0.08300,0.08300};
Double_t  resx[]={0.00120,0.00120,0.00380,0.00380,0.00200,0.00200};
Double_t cellx[]={0.00500,0.00500,0.01500,0.01500,0.00950,0.00950};
Double_t celly[]={0.00150,0.00150,0.00300,0.00300,0.00300,0.00300};
Double_t cellz[]={0.03000,0.03000,0.03000,0.03000,4.00000,4.00000};

//_________________________________________________________________________
void HitsTo(ClustAl_tl *trk,Int_t &ntrk,Int_t nt,TTree *TH,AliITS *ITS,
	    Float_t nsigmaT1,Float_t nsigmaT2,Float_t nsigmaT3,
	    Float_t nsigmaR1,Float_t nsigmaR2,Float_t nsigmaR3){

   // Local variables
   Int_t        nClustMin = 3;
   const Int_t  jmax = 50;
   Int_t        t,nb,nh,h,i,ly=-1,ld=-1,dt=-1,tt=-1,j=0;
   ClustAl_sl   ht[jmax];
   Int_t        id[3],ie[3],nht[jmax],track;
   Float_t      gp[3],lp[3],p[3],tof,gp2[3];
   Float_t      ns;
   Float_t      px[jmax],py[jmax],pz[jmax];
   TClonesArray *ITShits = ITS->Hits();
   AliITShit    *itsHit;
   // initilize and creat AliITSgeom data
   AliITSgeom   *gm  = ITS->GetITSgeom();
   AliITSgeom   &gm2 = *gm;

//   nt = 100; // lets speed things up.

   Float_t tran[]={90.e-4/3.,30.e-4/3.,780.e-4/3.}; // r,rphi,z
   Float_t  rot[]={0.4e-3/3.,3.0e-3/3.,100.e-3/3.}; // about r, rphi, z
   tran[0] *= nsigmaT1;
   tran[1] *= nsigmaT2;
   tran[2] *= nsigmaT3;
   rot[0]  *= nsigmaR1;
   rot[1]  *= nsigmaR2;
   rot[2]  *= nsigmaR3;

   printf("Randomize by tran=%fcm %fcm %fcm\n",tran[0],tran[1],tran[2]);
   printf("Randomize by rot=%frad %frad %frad\n",rot[0],rot[1],rot[2]);

   gm2.RandomCylindericalChange(tran,rot);

   for(i=0;i<jmax;i++){// zero for for first track.
       ht[i].xl = ht[i].xg = 0.0;
       ht[i].yl = ht[i].yg = 0.0;
       ht[i].zl = ht[i].zg = 0.0;
       ht[i].lay = ht[i].lad = ht[i].det = 0;
       px[i]    = 0.0;
       py[i]    = 0.0;
       pz[i]    = 0.0;
       nht[i]   = 0;
   } // end init to zero.

   j = 0;
   for(t=0;t<nt;t++){ // Loop over tracks.
      gAlice->ResetHits();
      nb = TH->GetEvent(t);
      nh = ITShits->GetEntriesFast();
      for(h=0;h<nh;h++){ // hits
         itsHit = (AliITShit *)ITShits->UncheckedAt(h);
	 track  = itsHit->fTrack;
	 if(track != tt){ // if new track
	     if(j>nClustMin)if(tt<nt||PRIMARYONLY==0){ 
                 // if there are enought points save data
		 trk[ntrk].track  = tt;
		 trk[ntrk].nclust = j+1;
		 ns               = 1./(Float_t (nht[0]));
		 trk[ntrk].px     = px[0]*ns;
		 trk[ntrk].py     = py[0]*ns;
		 trk[ntrk].pz     = pz[0]*ns;
		 trk[ntrk].pt     = sqrt(px[0]*px[0]+py[0]*py[0])*ns;
		 trk[ntrk].p      = sqrt(px[0]*px[0]+py[0]*py[0]+
					 pz[0]*pz[0]            )*ns;
		 trk[ntrk].clust  = new ClustAl_sl[j+1];
		 for(i=0;i<=j;i++){ // loop over all detector hit.
		     ns                     = 1./(Float_t(nht[i]));
		     trk[ntrk].clust[i].lay = ie[0] = ht[i].lay;
		     trk[ntrk].clust[i].lad = ie[1] = ht[i].lad;
		     trk[ntrk].clust[i].det = ie[2] = ht[i].det;
		     lp[0]                  = ht[i].xl*ns;
		     lp[1]                  = ht[i].yl*ns;
		     lp[2]                  = ht[i].zl*ns;
		     // set in detector resolution on local coordiante system
		     // z->zl rphi->xl yl-> thickness leave it alone
		     lp[0] = lp[0]/resx[ie[0]-1];
		     lp[0] = (Float_t)((Int_t)lp[0]);
		     lp[0] = lp[0]*resx[ie[0]-1] + 0.5*resx[ie[0]-1];
		     lp[2] = lp[2]/resz[ie[0]-1];
		     lp[2] = (Float_t)((Int_t)lp[2]);
		     lp[2] = lp[2]*resz[ie[0]-1] + 0.5*resz[ie[0]-1];
		     //
		     trk[ntrk].clust[i].xl  = lp[0];
		     trk[ntrk].clust[i].yl  = lp[1];
		     trk[ntrk].clust[i].zl  = lp[2];
		     gm2.LtoG(ie,lp,gp);
		     gm->LtoG(ie,lp,gp2);
		     trk[ntrk].clust[i].xg  = gp[0];
		     trk[ntrk].clust[i].yg  = gp[1];
		     trk[ntrk].clust[i].zg  = gp[2];
		 } // end for i: detectors clust loop
		 ntrk++; // set up for next track
	     } // end if j>nClustMin
	     for(i=0;i<jmax;i++){ // zero out structure
		 ht[i].xl = ht[i].xg = 0.0;
		 ht[i].yl = ht[i].yg = 0.0;
		 ht[i].zl = ht[i].zg = 0.0;
		 ht[i].lay = ht[i].lad = ht[i].det = 0;
		 px[i]    = 0.0;
		 py[i]    = 0.0;
		 pz[i]    = 0.0;
		 nht[i]   = 0;
	     } // end for i
	     j  = -1; // since ly=ld=dt=0 forces a j++ call.
	     ly = 0; // zero old detector values
	     ld = 0;
	     dt = 0;
	     tt = track;
	 } // end if track != tt
         itsHit->GetPositionG(gp[0],gp[1],gp[2],tof);
         itsHit->GetDetectorID(id[0],id[1],id[2]);
         gm->GtoL(id,gp,lp);
	 itsHit->GetMomentumG(p[0],p[1],p[2]);
         if(!(id[0]==ly && id[1]==ld && id[2]==dt)){//if Not the same detector
	     j++; // incriment detector counter j.
	     ly = id[0];
	     ld = id[1];
	     dt = id[2];
	 } // end if id==idold
	 ht[j].lay = id[0];
	 ht[j].lad = id[1];
	 ht[j].det = id[2];
	 ht[j].xl += lp[0];
	 ht[j].yl += lp[1];
	 ht[j].zl += lp[2];
	 ht[j].xg  = gp[0];
	 ht[j].yg  = gp[1];
	 ht[j].zg  = gp[2];
	 px[j]    += p[0];
	 py[j]    += p[1];
	 pz[j]    += p[2];
	 nht[j]++;
      } // end for h
  } // end for t: track loop
  return;
}
//______________________________________________________________________
void HitsToClustAl(ClustAl_tl *trk,Int_t &ntrk,Int_t nt,TTree *TH,AliITS *ITS,
		   Float_t fraction){

   // Local variables
   Int_t        nClustMin = 3,Iseed;
//   Int_t        icount=0,icountMAX=100;
   const Int_t  jmax = 50;
   Int_t        t,nb,nh,h,i,ly=-1,ld=-1,dt=-1,tt=-1,j=0;
   ClustAl_sl   ht[jmax];
   Int_t        id[3],ie[3],nht[jmax],track;
   Float_t      gp[3],lp[3],tof;
   Float_t      ns;
   TClonesArray *ITShits = ITS->Hits(),*Prt = gAlice->Particles();
   AliITShit    *itsHit;
   TParticle    *prt;
   // initilize and creat AliITSgeom data
   AliITSgeom   *gm  = ITS->GetITSgeom();

   for(i=0;i<jmax;i++){// zero for for first track.
       ht[i].xl  = 0.0;
       ht[i].yl  = 0.0;
       ht[i].zl  = 0.0;
       ht[i].lay = ht[i].lad = ht[i].det = 0;
       nht[i]    = 0;
   } // end init to zero.

   j = 0;
   Iseed = time(0);
   printf("HitsToclustAl: Iseed=%d ",Iseed);
   gRandom->SetSeed(Iseed);
   printf("gRandom->Rndm(1)=%f\n",gRandom->Rndm(1));
   for(t=0;t<nt;t++){ // Loop over tracks.
      if(fraction<gRandom->Rndm(1)) continue; // skip some tracks
      gAlice->ResetHits();
      nb = TH->GetEvent(t);
      nh = ITShits->GetEntriesFast();
      for(h=0;h<nh;h++){ // hits
         itsHit = (AliITShit *)ITShits->UncheckedAt(h);
	 track  = itsHit->fTrack;
         prt= (TParticle *)Prt->UncheckedAt(track);
         if(prt->P()<MOMENTUMCUT) continue;
	 if(track != tt){ // if new track
	     if(j>nClustMin)if(tt<nt||PRIMARYONLY==0){ 
//		 if(icount>icountMAX) return;
//		 icount++;
                 // if there are enought points save data
		 prt              = (TParticle *)Prt->UncheckedAt(tt);
		 trk[ntrk].track  = tt;
		 trk[ntrk].nclust = j+1;
		 ns               = 1./(Float_t (nht[0]));
		 trk[ntrk].px     = prt->Px();
		 trk[ntrk].py     = prt->Py();
		 trk[ntrk].pz     = prt->Pz();
		 trk[ntrk].pt     = prt->Pt();
		 trk[ntrk].p      = prt->P();
		 trk[ntrk].clust  = new ClustAl_sl[j+1];
		 for(i=0;i<=j;i++){ // loop over all detector hit.
		     ns                     = 1./(Float_t(nht[i]));
		     trk[ntrk].clust[i].lay = ie[0] = ht[i].lay;
		     trk[ntrk].clust[i].lad = ie[1] = ht[i].lad;
		     trk[ntrk].clust[i].det = ie[2] = ht[i].det;
		     lp[0]                  = ht[i].xl*ns;
		     lp[1]                  = ht[i].yl*ns;
		     lp[2]                  = ht[i].zl*ns;
		     // set in detector resolution on local coordiante system
		     // z->zl rphi->xl yl-> thickness leave it alone
		     lp[0] = lp[0]/resx[ie[0]-1];
		     lp[0] = (Float_t)((Int_t)lp[0]);
		     lp[0] = lp[0]*resx[ie[0]-1] + 0.5*resx[ie[0]-1];
		     lp[2] = lp[2]/resz[ie[0]-1];
		     lp[2] = (Float_t)((Int_t)lp[2]);
		     lp[2] = lp[2]*resz[ie[0]-1] + 0.5*resz[ie[0]-1];
		     //
		     trk[ntrk].clust[i].xl  = lp[0];
		     trk[ntrk].clust[i].yl  = lp[1];
		     trk[ntrk].clust[i].zl  = lp[2];
		 } // end for i: detectors clust loop
		 ntrk++; // set up for next track
	     } // end if j>nClustMin
	     for(i=0;i<jmax;i++){ // zero out structure
		 ht[i].xl = ht[i].xg = 0.0;
		 ht[i].yl = ht[i].yg = 0.0;
		 ht[i].zl = ht[i].zg = 0.0;
		 ht[i].lay = ht[i].lad = ht[i].det = 0;
		 nht[i]   = 0;
	     } // end for i
	     j  = -1; // since ly=ld=dt=0 forces a j++ call.
	     ly = 0; // zero old detector values
	     ld = 0;
	     dt = 0;
	     tt = track;
	 } // end if track != tt
         itsHit->GetPositionG(gp[0],gp[1],gp[2],tof);
         itsHit->GetDetectorID(id[0],id[1],id[2]);
         gm->GtoL(id,gp,lp);
         if(!(id[0]==ly && id[1]==ld && id[2]==dt)){//if Not the same detector
	     j++; // incriment detector counter j.
	     ly = id[0];
	     ld = id[1];
	     dt = id[2];
	 } // end if id==idold
	 ht[j].lay = id[0];
	 ht[j].lad = id[1];
	 ht[j].det = id[2];
	 ht[j].xl += lp[0];
	 ht[j].yl += lp[1];
	 ht[j].zl += lp[2];
	 nht[j]++;
      } // end for h
  } // end for t: track loop
  return;
}
//______________________________________________________________________
void SetDetectorResolusion(Int_t l,Float_t *xl){
    //set in detector resolution on local coordiante system
    // z->zl rphi->xl yl-> thickness leave it alone
    // changed from res to cell for detector simulations.

    l--;
    xl[0] = xl[0]/cellx[l];
    xl[0] = (Float_t)((Int_t)xl[0]);
    xl[0] = xl[0]*cellx[l] + 0.5*cellx[l];
    xl[2] = xl[2]/cellz[l];
    xl[2] = (Float_t)((Int_t)xl[2]);
    xl[2] = xl[2]*cellz[l] + 0.5*cellz[l];
    xl[1] = 0.0;
}
//______________________________________________________________________
void FillAliITSAlignmentTrack(AliITSAlignmentTrack *trk,Int_t &ntrk,Int_t nt,
			      TTree *TH,AliITS *ITS,Float_t fraction){
    const Int_t iMAX=20,jMAX=50,nClustMin=3;
    Int_t   ist[iMAX],imax;
    Int_t   t,i,h,j,nh;
    Int_t   lay,lad,det,index,track,Iseed;
    Int_t   indexl[jMAX];
    Float_t tof,xl[3],exl[3][3];
    AliITSstatistics *sxl[jMAX],*szl[jMAX];
    TParticle    *prt;
    AliITShit    *itshit;
    TClonesArray *ITShits = ITS->Hits();
    AliITSgeom   *gm = ITS->GetITSgeom();
    TClonesArray *Prt = gAlice->Particles();

    printf("Entered FillAliITSAlignmentTrack\n");
    for(i=0;i<jMAX;i++){
	sxl[i] = new AliITSstatistics(2);
	szl[i] = new AliITSstatistics(2);
    }  // end for i
    for(i=0;i<3;i++)for(j=0;j<3;j++) exl[i][j] = 0.0;
    Iseed = time(0);
    gRandom->SetSeed(Iseed);
    printf("Entering track loop. Iseed=%d\n",Iseed);
    for(t=0;t<nt;t++) if(fraction>=gRandom->Rndm(1)){
//	gAlice->ResetHits();
	TH->GetEvent(t);
	nh = ITShits->GetEntriesFast();
	i = 0;
	ist[0] = 0;
	printf("Entering hit loop for track=%d\n",t);
	for(h=1;h<nh&&i<iMAX-1;h++) 
	    if(((AliITShit *)ITShits->UncheckedAt(ist[i]))->GetTrack() !=
	       ((AliITShit *)ITShits->UncheckedAt(h))->GetTrack()) ist[++i]=h;
	ist[++i]=nh;
	imax = i+1;
	printf("entering loop i nh=%d imax=%d\n",nh,imax);
	for(i=1;i<imax;i++){ // loop over tracks from primary track
	    printf("Getting hit ist[%d-1]=%d\n",i,ist[i-1]);
	    itshit = (AliITShit *)ITShits->UncheckedAt(ist[i-1]);
	    track  = itshit->GetTrack();
	    prt    = (TParticle *)Prt->UncheckedAt(track);
	    if(prt->P()<MOMENTUMCUT) continue;
	    printf("pass P cut P=%f\n",prt->P());
	    for(j=0;j<jMAX;j++){
		sxl[j]->Reset();
		szl[j]->Reset();
	    } // end for j
	    printf("Exiting the reset loop, ist[%d-1]=%d\n",i,ist[i-1]);
	    j = 0;
	    indexl[0]=((AliITShit *)ITShits->UncheckedAt(ist[i-1]))->
		                                                GetDetector();
	    printf("looping over this track's hits ist[%d-1]=%d to ist[%d]=%d\n",i,ist[i-1],i,ist[i]);
	    for(h=ist[i-1];h<ist[i];h++){ // loop over hits
		itshit = (AliITShit *)ITShits->UncheckedAt(h);
		index = itshit->GetDetector();
		if(indexl[j]!=index) j++;
		indexl[j] = index;
		itshit->GetPositionL(xl[0],xl[1],xl[2],tof);
		sxl[j]->AddValue((Double_t) xl[0],1.0);
		szl[j]->AddValue((Double_t) xl[2],1.0);
	    } // end for h
	    if(j<nClustMin) continue; // get next track
	    printf("Setting up tracks ntrk=%d #digits for track=%d\n",ntrk,j+1);
	    trk[ntrk].CreatePoints(j+1);
	    printf("CreatePoints with ntrk=%d and j+1=%d\n",ntrk,j+1);
	    trk[ntrk].SetTParticle(prt);
	    printf("SetTParticle\n");
	    trk[ntrk].SetTrackNumber(track);
	    printf("SetTrackNumber=%d\n",track);
	    for(h=0;h<=j;h++){
		xl[0] = sxl[h]->GetMean();
		xl[1] = 0.0;
		xl[2] = szl[h]->GetMean();
		gm->GetModuleId(indexl[h],lay,lad,det);
		printf("setting detector resolusion for index[%d]=%d (%d,%d,%d)\n",h,indexl[h],lay,lad,det);
		SetDetectorResolusion(lay,xl);
		printf("Detector resolusion set\n");
		exl[0][0] = 1./(resx[lay-1]*resx[lay-1]);
		exl[1][1] = 12./(celly[lay-1]*celly[lay-1]);
		exl[2][2] = 1./(resz[lay-1]*resz[lay-1]);
		printf("Adding point to trk[%d] indexl[%d]=%d\n",ntrk,h,indexl[h]);
		trk[ntrk].AddPointLastL(indexl[h],(Float_t *)xl,
					                     (Float_t **)exl);
	    } // end for h
	    ntrk++; // set up for next track
	} // end for i
    } // end for t // end if(fraction>=gRandom->Rndm(1))
    for(i=0;i<jMAX;i++){
	delete sxl[i];
	delete szl[i];
    }  // end for i
}
//______________________________________________________________________
void FillGlobalPositions(ClustAl_tl *trk,Int_t ntrk,AliITSgeom *g){
    Int_t   i,j,id[3];
    Float_t lx[3],gx[3];

    for(i=0;i<ntrk;i++) for(j=0;j<trk[i].nclust;j++){
	lx[0] = trk[i].clust[j].xl;
	lx[1] = trk[i].clust[j].yl;
	lx[2] = trk[i].clust[j].zl;
	id[0] = trk[i].clust[j].lay;
	id[1] = trk[i].clust[j].lad;
	id[2] = trk[i].clust[j].det;
	g->LtoG(id,lx,gx);
	trk[i].clust[j].xg = gx[0];
	trk[i].clust[j].yg = gx[1];
	trk[i].clust[j].zg = gx[2];
    } // end for i,j
    return;
}
//______________________________________________________________________
void PlotGeomChanges(AliITSgeom *gt,AliITSgeom *gc,TFile *Hfile,Float_t *Rdta){
    Int_t    ly,ld,dt;
    Float_t  x,y,z;
    Double_t xd,yd,zd;
    Double_t rt,rc,phit,phic,zt,zc;
    Float_t  A0,B0,C0,A1,B1,C1;
    Double_t dr,drphi,dz,dphi;
    Double_t drmi,drma,drphimi,drphima,dzmi,dzma;

    AliITSstatistics *R    = new AliITSstatistics(2);
    AliITSstatistics *RPhi = new AliITSstatistics(2);
    AliITSstatistics *Z    = new AliITSstatistics(2);
    AliITSstatistics *A    = new AliITSstatistics(2);
    AliITSstatistics *B    = new AliITSstatistics(2);
    AliITSstatistics *C    = new AliITSstatistics(2);
    TH1F *Gr  = new TH1F("Gr","Radial Displacement (cm)",500,-1.0,1.0);
    drmi = -1.0;
    drma = +1.0;
    Gr->Sumw2();
    Gr->SetXTitle("Displacement (cm)");
    TH1F *Grphi  = new TH1F("Grphi","RPhi Displacement (cm)",500,-1.0,1.0);
    drphimi = -1.0;
    drphima = +1.0;
    Grphi->Sumw2();
    Grphi->SetXTitle("Displacement (cm)");
    TH1F *Gz  = new TH1F("Gz","Z Displacement (cm)",500,-1.0,1.0);
    dzmi = -1.0;
    dzma = +1.0;
    Gz->Sumw2();
    Gz->SetXTitle("Displacement (cm)");

    Float_t weight=1.0;
    for(ly=1;ly<=gt->GetNlayers();ly++)for(ld=1;ld<=gt->GetNladders(ly);ld++)
    for(dt=1;dt<=gt->GetNdetectors(ly);dt++){
	gt->GetTrans(ly,ld,dt,x,y,z);xd=x;yd=y;zd=z;
	gt->GetAngles(ly,ld,dt,A0,B0,C0);
	rt    = TMath::Hypot(yd,xd);
	phit  = TMath::ATan2(yd,xd);
	if(phit<0.0) phit += 2.0*TMath::Pi();
	zt    = zd;
	gc->GetTrans(ly,ld,dt,x,y,z);xd=x;yd=y;zd=z;
	gc->GetAngles(ly,ld,dt,A1,B1,C1);
	rc    = TMath::Hypot(yd,xd);
	phic  = TMath::ATan2(yd,xd);
	if(phic<0.0) phic += 2.0*TMath::Pi();
	zc    = zd;
	dr    = rt-rc;
	dphi  = phit - phic;
	if(dphi>2.0*TMath::Pi()) dphi -= 2.0*TMath::Pi();
	drphi = 0.5*(rt + rc)*dphi; // change in phi as measured in rphi
	dz    = zt - zc;
	Gr->Fill(dr,1.0);
	Grphi->Fill(drphi,1.0);
	Gz->Fill(dz,1.0);
	if(dr>=drmi&&dr<drma) R->AddValue(dr,1.0);
	if(drphi>=drphimi&&drphi<drphima) RPhi->AddValue(drphi,1.0);
	if(dz>=dzmi&&dz<dzma) Z->AddValue(dz,1.0);
	A->AddValue(A0-A1,weight);
	B->AddValue(B0-B1,weight);
	C->AddValue(C0-C1,weight);
    } // end for ly,ld,dt
    Hfile->Write();
    Rdta[0] = R->GetRMS();
    Rdta[1] = RPhi->GetRMS();
    Rdta[2] = Z->GetRMS();
    Rdta[3] = A->GetRMS();
    Rdta[4] = B->GetRMS();
    Rdta[5] = C->GetRMS();
    printf("PlotGeomChanges: RMS(r) = %f RMS(rphi) = %f RMS(z) = %f\n",
	   R->GetRMS(),RPhi->GetRMS(),Z->GetRMS());
    delete Gr;
    delete Grphi;
    delete Gz;
    delete R;
    delete RPhi;
    delete Z;
    delete A;
    delete B;
    delete C;
    return;
}
//______________________________________________________________________
Int_t FitTrackToLine(ClustAl_tl &trk){
   // Local Variables
   Int_t   i,l,j=0;
   Double_t x,y,z,wx,wy,dx,dy;
   Double_t a,b,c,d;
   Double_t xb,sxb,zxb;
   Double_t yb,syb,zyb;
   Double_t xzb,yzb,z2xb,z2yb;
   Double_t sum,phi=0.0,b0,d0;

   trk.qual = -1.0;
   if(trk.nclust<3) return -1;

   b  = (trk.clust[0].xg-trk.clust[trk.nclust].xg)/
        (trk.clust[0].zg-trk.clust[trk.nclust].zg);
   d  = (trk.clust[0].yg-trk.clust[trk.nclust].yg)/
        (trk.clust[0].zg-trk.clust[trk.nclust].zg);
   do{
     b0 = b;
     d0 = d;
     xb  = 0.0; sxb = 0.0; zxb  = 0.0;
     yb  = 0.0; syb = 0.0; zyb  = 0.0;
     xzb = 0.0; yzb = 0.0; z2xb = 0.0; z2yb = 0.0;
     for(i=0;i<trk.nclust;i++){
      l    = trk.clust[i].lay - 1;  // zero based
      x    = trk.clust[i].xg;
      y    = trk.clust[i].yg;
      z    = trk.clust[i].zg;
      phi  = atan2(y,x);
      wx   = 1.0/(resx[l]*resx[l]*cos(phi)*cos(phi)+b0*b0*resz[l]*resz[l]);// 1.0/rms^2
      wy   = 1.0/(resx[l]*resx[l]*sin(phi)*sin(phi)+d0*d0*resz[l]*resz[l]);// 1.0/rms^2
      xb   += x*wx;
      sxb  += wx;
      zxb  += z*wx;
      yb   += y*wy;
      syb  += wy;
      zyb  += z*wy;
      xzb  += x*z*wx;
      yzb  += y*z*wy;
      z2xb += z*z*wx;
      z2yb += z*z*wy;
     } // end for i
     dx = zxb*zxb - z2xb*sxb;
     if(dx!=0.0){
	 a  = zxb*xzb - xb*z2xb;
	 a  = a/dx;
	 b  = xb*zxb  - sxb*xzb;
	 b  = b/dx;
     }else{
	  printf("FitTrackToLine: Error dx=0 track=%d zxb=%f z2xb=%f sxb=%f\n",
		 trk.track,zxb,z2xb,sxb);
	  a = 0.0;
	  b = 1E10;
	  j++;
	  if(j>100) return(0);
     } // end if dx!=0.0
     dy = zyb*zyb - z2yb*syb;
     if(dy!=0.0){
	 c  = zyb*yzb - yb*z2yb;
	 c  = c/dy;
	 d  = yb*zyb  - syb*yzb;
	 d  = d/dy;
     }else{
	  printf("FitTrackToLine: Error dy=0 track=%d zyb=%f z2yb=%f syb=%f\n",
		 trk.track,zyb,z2yb,syb);
	  c = 0.0;
	  d = 1.E10;
	  j++;
	  if(j>100) return(0);
     } // end if dy!=0.0
   } while(fabs(b0-b)<1.E-5 && fabs(d0-d)<1.E-5);
   trk.a = a;
   trk.b = b;
   trk.c = c;
   trk.d = d;
   sum   = 0.0;
   for(i=0;i<trk.nclust;i++){
      l  = trk.clust[i].lay - 1;  // zero based
      x  = trk.clust[i].xg;
      y  = trk.clust[i].yg;
      z  = trk.clust[i].zg;
      x  = x - (a + b*z);  // x=a+bz
      y  = y - (c + d*z);  // y=c+dz
      wx =  resx[l]*resx[l]*cos(phi)*cos(phi) + b*b * resz[l]*resz[l];
      wy =  resx[l]*resx[l]*sin(phi)*sin(phi) + d*d * resz[l]*resz[l];
      if(wx==0.0||wy==0.0) {
	  printf("FitTrackToLine: error wx=%f wy=%f trk.clust[%d].lay=%d\n",
                 wx,wy,i,trk.clust[i].lay);
      }// end if sxz or syz == 0.
      sum += x*x/wx + y*y/wy;
   } // end for i #2
   trk.qual = sum;
   if(trk.nclust>2) trk.qual /= Double_t (trk.nclust-2);
   // per degree of freedom 2 in xz 2 in yz.
   return 0;
}
//______________________________________________________________________
void LtoLline(const Int_t *id1,const Int_t *id2,
	      Double_t a1,Double_t b1,Double_t c1,Double_t d1,
	      Double_t &a2,Double_t &b2,Double_t &c2,Double_t &d2,
	      AliITSgeom *gm){
    Double_t x11[3],x12[3],x21[3],x22[3],h;

    x11[0] = a1;
    x11[1] = 0.0;
    x11[2] = c1;
    x12[0] = a1+b1;
    x12[1] = 1.0;
    x12[2] = c1+d1;
    gm->LtoL(id1,id2,x11,x21);
    gm->LtoL(id1,id2,x12,x22);
    h = x21[1] - x22[1];
    if(h!=0.0){
	b2 = (x21[0] - x22[0])/h;
	d2 = (x21[2] - x22[2])/h;
	a2 = x21[0] - b2*x21[1];
	c2 = x21[2] - d2*x21[1];
    }else{
	printf("LtoLline: error line in plane of detector\n");
    } // end if h!=0
    return;
}
//______________________________________________________________________
Int_t FitTrackToLineL(ClustAl_tl &trk,AliITSgeom *gm){
   // Local Variables
   Int_t   i,id0[3],id[3];
   Double_t wx/*,wy*/,wz;
   Double_t a,b,c,d;
   Double_t xg[3],xl[3],x2g[3],x2l[3];
   AliITSstatistics2 *Fx  = new AliITSstatistics2(2);
   AliITSstatistics2 *Fz  = new AliITSstatistics2(2);

   trk.qual = -1.0;
   if(trk.nclust<3) return -1;

   Int_t Npts = trk.nclust;
   id0[0] = trk.clust[0].lay;
   id0[1] = trk.clust[0].lad;
   id0[2] = trk.clust[0].det;
   for(i=0;i<Npts;i++){
       id[0] = trk.clust[i].lay;
       id[1] = trk.clust[i].lad;
       id[2] = trk.clust[i].det;
       xg[0] = trk.clust[i].xg;
       xg[1] = trk.clust[i].yg;
       xg[2] = trk.clust[i].zg;
       gm->GtoL(id0,xg,xl);
       wx = 1.0/(resx[id[0]-1]*resx[id[0]-1]);
       wz = 1.0/(resz[id[0]-1]*resz[id[0]-1]);
       Fx->AddValue(xl[0],xl[1],wx);
       Fz->AddValue(xl[2],xl[1],wz);
   } // end for i
   trk.qual  = Fx->FitToLine(a,b);
   trk.qual += Fz->FitToLine(c,d);
   trk.a0 = a;
   trk.b0 = b;
   trk.c0 = c;
   trk.d0 = d;
   xl[0]  = a;
   xl[1]  = 0.0;
   xl[2]  = c;
   x2l[0] = a+b;
   x2l[1] = 1.0;
   x2l[2] = c+d;
   gm->LtoG(id0,xl,xg);
   gm->LtoG(id0,x2l,x2g);
   c = xg[2] - x2g[2];
   if(c!=0.0){
       b = (xg[0] - x2g[0])/c;
       d = (xg[1] - x2g[1])/c;
       a = xg[0] - b*xg[2];
       c = xg[1] - d*xg[2];
       trk.a = a;
       trk.b = b;
       trk.c = c;
       trk.d = d;
   }else{
       return -1;
   }// end if c!=0.0
   return 0;
}
//______________________________________________________________________
Int_t FindCircleCenter(Double_t &x0,Double_t &y0,Double_t x1,Double_t y1,
		       Double_t  x2,Double_t  y2,Double_t x3,Double_t y3){
////////////////////////////////////////////////////////////////////////
//     This was derived as folows. Given three non-linear points find
// the circle that is therefor defined by those three non-linar points.
// Assume that the circle is centers at x0,y0 and has a radous R. Then
// (1) R^2 = (x1-x0)^2 + (y1-y0)^2
// (2) R^2 = (x2-x0)^2 + (y2-y0)^2
// (3) R^2 = (x3-x0)^2 + (y3-y0)^2.
// Now consider the two equations derived from the above
// (1) - (2) = x1^2 - x2^2 -2x0(x1-x2) + y1^2 - y2Y2 -2Y0(y1-y2) = 0
// (1) - (3) = x1^2 - x3^2 -2x0(x1-x3) + y1^2 - y3Y2 -2Y0(y1-y3) = 0
// solving these two equations for x0 and y0 gives
// x0 = +{(y1-y2)*(y1-y3)*(y1-y3)+x1*x1*(y1-y3)+x2*x2*(y3-y1)+x3*x3*(y1-y2)}/2d
// y0 = -{(x1-x2)*(x1-x3)*(x1-x3)+y1*y1*(x1-x3)+y2*y2*(x3-x1)+y3*y3*(x1-x2)}/2d
// with d = (x1-x2)*(y1-y3) - (x1-x3)*(y1-y2)
////////////////////////////////////////////////////////////////////////
    Double_t d;

    d = (x1-x2)*(y1-y3) - (x1-x3)*(y1-y2);
    if(d==0.0) return 0;  // fits to a line!

    x0 = (y1-y2)*(y1-y3)*(y1-y3)+x1*x1*(y1-y3)+x2*x2*(y3-y1)+x3*x3*(y1-y2);
    y0 = (x1-x2)*(x1-x3)*(x1-x3)+y1*y1*(x1-x3)+y2*y2*(x3-x1)+y3*y3*(x1-x2);
    x0 = +0.5*x0/d;
    y0 = -0.5*y0/d;

    return 1;
}
//______________________________________________________________________
void FitAllTracks(ClustAl_tl *trk,Int_t ntrk,Float_t *v0,AliITSgeom *gm,
		  const char *sfile,TFile *Hfile,Float_t *Fdta,Int_t *Ndta){
   // Local Variables
   Int_t    i,j,k,id0[3],id[3];
   Double_t xm,ym,zm,dt,ad/*,rh,rp,rphih,rphip*/;
   Double_t tp,tx,rx,rz;
   Double_t a,b,c,d,qual;
   Double_t xl[3],xg[3];
   Double_t trkqualMAX = 50.0;
   Bool_t   printit = kFALSE;
   char     filename[80],hid[10];
   Int_t    Nqualmax = 10,NqualgtMAX=0;
   Int_t    Nqualless[10];
// Distrobution statisitics
   AliITSstatistics *Sad    = new AliITSstatistics(4);
   AliITSstatistics *SrxzL  = new AliITSstatistics(4);
   AliITSstatistics *SrxL[6];
   AliITSstatistics *SrzL[6];
   for(i=0;i<6;i++) {
       SrxL[i] = new AliITSstatistics(4);
       SrzL[i] = new AliITSstatistics(4);
   } // end for i

   for(i=0;i<Nqualmax;i++) Nqualless[i] = 0;

// Setup histograms

   TH2F *Ttqp  = new TH2F("Ttqp","Track quality vs. momentum",
			 500,0.0,50.0,500,0.0,10.0);
   Ttqp->Sumw2();
   Ttqp->SetXTitle("Chi^2/degree freedom");
   Ttqp->SetYTitle("Momentum GeV/c");
   TH2F *Tdttq = new TH2F("Tdttq",
			  "Distance from true vertex vs. track quality",
			  500,0.0,0.10,100,0.0,50.0);
   Tdttq->Sumw2();
   Tdttq->SetXTitle("Distanct to true vertex (cm)");
   Tdttq->SetYTitle("Chi^2/degree freedom");
   TH2F *Tadtq = new TH2F("Tadtq",
                     "atan(sqrt(b*b+d*d))-atan(pt/|pz|) vs. track quality",
			  500,-0.03,0.03,500,0.0,50.0);
   Tadtq->Sumw2();
   Tadtq->SetXTitle("theta fitted line - atan(pt/p) (rad)");
   Tadtq->SetYTitle("Chi^2/degree freedom");
   TH2F *Taddt = new TH2F("Taddt",
         "atan(sqrt(b*b+d*d))-atan(pt/|pz|) vs. Distance from true vertex",
			  500,-0.03,0.03,200,0.0,0.10);
   Taddt->Sumw2();
   Taddt->SetXTitle("theta fitted line - atan(pt/p) (rad)");
   Taddt->SetYTitle("distance to true vertex (cm)");
   TH2F *Tadp = new TH2F("Tadp",
         "atan(sqrt(b*b+d*d))-atan(pt/|pz|) vs. momentum",
			  500,-0.03,0.03,200,0.0,10.0);
   Tadp->Sumw2();
   Tadp->SetXTitle("theta fitted line - atan(pt/p) (rad)");
   Tadp->SetYTitle("Momentum (GeV/c)");

   TH2F *TrxrzL[7];
//   TH2F *TrrprzG[7];
   for(i=0;i<7;i++){
       sprintf(filename,"Layer %1.1d: xi-x(i) vs. zi-z(i) local",i+1);
       if(i==6) sprintf(filename,"Sum of all layers: xi-x(i) "
			"vs. zi-z(i) local");
       sprintf(hid,"TrxrzL%1.1d",i+1);
       TrxrzL[i] = new TH2F(hid,filename,500,-1.0,1.0,500,-1.0,1.0);
       TrxrzL[i]->Sumw2();
       TrxrzL[i]->SetXTitle("xi-x(i) (cm) local");
       TrxrzL[i]->SetYTitle("zi-z(i) (cm) local");
//
//       sprintf(filename,"Layer %1.1d: rPhii-rPhi(i) vs. zi-z(i) local",i+1);
//       if(i==6) sprintf(filename,"Sum of all layers: rPhii-rPhi(i) "
//			"vs. zi-z(i) local");
//       sprintf(hid,"TrrprzG%1.1d",i+1);
//       TrrprzG[i] = new TH2F(hid,filename,500,-1.0,1.0,500,-1.0,1.0);
//       TrrprzG[i]->Sumw2();
//       TrrprzG[i]->SetXTitle("rPhii-rPhi(i) (cm) local");
//       TrrprzG[i]->SetYTitle("zi-z(i) (cm) local");
   } // end for i

   TH1D *Tptqp[10];
   TH1D *Ttqall = new TH1D("Ttqall","Track quality for all tracks fit",
			   500,0.0,1000.0);
   Ttqall->Sumw2();
   Ttqall->SetXTitle("Track quality: Chi squared per degree freedom");

// Fit each track and fill histograms and the like.

   for(i=0;i<ntrk;i++){
       if(trk[i].p < MOMENTUMCUT) continue;
      if(FitTrackToLineL(trk[i],gm)!=0) continue;  // fit track to a line
      qual = trk[i].qual;
      if(qual<0.0) continue;
      // Fill histograms and statistics before cource chi squared cut
      Ttqall->Fill(qual,1.0);
      for(j=0;j<Nqualmax;j++) if(qual<(Double_t) j+1) Nqualless[j]++;
//      if(qual>trkqualMAX) {
//	  NqualgtMAX++;
//	  continue;
//      } // end if iqual>trkqualMAX=50.0
      a = trk[i].a;
      b = trk[i].b;
      c = trk[i].c;
      d = trk[i].d;
      zm  = b*(v0[0]-a) + d*(v0[1]-c);
      dt  = (b*b + d*d);
      if(dt!=0.0) zm /= dt;
      else{
	  printf("FitAllTracks: trk[%d].b=%f trk[%d].d=%f\n",i,b,i,d);
	  zm=0.0;
      }// end if else dt!=0.0
      xm  = a + b*zm;
      ym  = c + d*zm;
      dt   = sqrt((xm-v0[0])*(xm-v0[0]) + (ym-v0[1])*(ym-v0[1]) +
		 (zm-v0[2])*(zm-v0[2]));
      tx  = atan(sqrt(b*b + d*d));
      tp  = atan2(trk[i].pt,fabs(trk[i].pz));
      ad  = tx-tp;
      // Fill histograms and statistics
      Tadtq->Fill(ad,qual,1.0);
      Sad->AddValue(ad,1.0);
      Taddt->Fill(ad,dt,1.0);
      Tadp->Fill(ad,trk[i].p,1.0);
      Tdttq->Fill(dt,qual,1.0);
      Ttqp->Fill(qual,trk[i].p,1.0);
      //
      id0[0] = trk[i].clust[0].lay;
      id0[1] = trk[i].clust[0].lad;
      id0[2] = trk[i].clust[0].det;
      xg[0] = trk[i].clust[0].xg;
      xg[1] = trk[i].clust[0].yg;
      xg[2] = trk[i].clust[0].zg;
      gm->GtoL(id0,xg,xl);
      rx = xl[0] - trk[i].a0 - trk[i].b0 * xl[1];
      rz = xl[2] - trk[i].c0 - trk[i].d0 * xl[1];
      SrxL[id0[0]-1]->AddValue(rx,1.0);
      SrzL[id0[0]-1]->AddValue(rz,1.0);
      SrxzL->AddValue(rx,1.0);
      SrxzL->AddValue(rz,1.0);
      TrxrzL[6]->Fill(rx,rz,1.0);
      TrxrzL[id0[0]-1]->Fill(rx,rz,1.0);
      for(j=1;j<trk[i].nclust;j++){
	  id[0] = trk[i].clust[j].lay;
	  id[1] = trk[i].clust[j].lad;
	  id[2] = trk[i].clust[j].det;
	  LtoLline(id0,id,trk[i].a0,trk[i].b0,trk[i].c0,trk[i].d0,a,b,c,d,gm);
	  xg[0] = trk[i].clust[j].xg;
	  xg[1] = trk[i].clust[j].yg;
	  xg[2] = trk[i].clust[j].zg;
	  gm->GtoL(id,xg,xl);
	  rx = xl[0] - a - b * xl[1];
	  rz = xl[2] - c - d * xl[1];
	  SrxL[id[0]-1]->AddValue(rx,1.0);
	  SrzL[id[0]-1]->AddValue(rz,1.0);
	  SrxzL->AddValue(rx,1.0);
	  SrxzL->AddValue(rz,1.0);
	  TrxrzL[6]->Fill(rx,rz,1.0);
	  TrxrzL[id[0]-1]->Fill(rx,rz,1.0);
      } // end for j
   } // end for i

// Write out information

   for(i=0;i<6;i++){
       Fdta[4*i+0] = SrxL[i]->GetRMS();
       Fdta[4*i+1] = SrxL[i]->GetErrorRMS();
       Fdta[4*i+2] = SrzL[i]->GetRMS();
       Fdta[4*i+3] = SrzL[i]->GetErrorRMS();
   } // end for i
   for(i=0;i<10;i++) Ndta[i] = Nqualless[i];
   Ndta[10] = Sad->GetN();
   Ndta[11] = NqualgtMAX;

   printf("FitAllTracks: %d tracks cut leaving %d, with a chi squared <%f\n",
	  NqualgtMAX,Sad->GetN(),trkqualMAX);
   printf("The number of tracks with chi squared <1");
   for(i=1;i<Nqualmax;i++) printf(",%d",i+1);
   printf("\nFitAllTracks: %d",Nqualless[0]);
   for(i=1;i<Nqualmax;i++) printf(",%d",Nqualless[i]);
   printf("\n");
   printf("FitAllTracks: RMSs of ad rxz = %e+-%e %e+-%e\n",
	  Sad->GetRMS(),Sad->GetErrorRMS(),
	  SrxzL->GetRMS(),SrxzL->GetErrorRMS());

   printf("FitAllTracks: Residuals by layer x=-rphi z ");
   for(i=0;i<6;i++) printf(":%d:%e+-%e %e+-%e ",i+1,
			   SrxL[i]->GetRMS(),SrxL[i]->GetErrorRMS(),
			   SrzL[i]->GetRMS(),SrzL[i]->GetErrorRMS());
   printf("\n");

// Setup and fill projections

   for(i=0;i<10;i++){
       xm = 0.5*Double_t(i);
       ym = xm+0.5;
       j  = Ttqp->GetYaxis()->FindBin(xm);
       k  = Ttqp->GetYaxis()->FindBin(ym);
       sprintf(filename,"Track Quality for %3.1f<p<%3.1f",xm,ym);
       Tptqp[i] = Ttqp->ProjectionX(filename,j,k,"E");
       Tptqp[i]->SetXTitle("Track Quality (Chi squared/df)");
   } // end for i

   TH1D *Ttq  = Ttqp->ProjectionX("Ttq",0,Ttqp->GetNbinsY()+1,"E");
   Ttq->SetXTitle("Chi^2/degree freedom");
   TH1D *Tp   = Ttqp->ProjectionY("Tp", 0,Ttqp->GetNbinsX()+1,"E");
   Tp->SetXTitle("Momentum GeV/c");
   TH1D *Tdt  = Tdttq->ProjectionX("Tdt",0,Tdttq->GetNbinsY()+1,"E");
   Tdt->SetXTitle("Distanct to true vertex (cm)");
   TH1D *Tad  = Tadtq->ProjectionX("Tad",0,Tadtq->GetNbinsY()+1,"E");
   Tad->SetXTitle("theta fitted line - atan(pt/p) (rad)");
   TH1D *TrxL[7],*TrzL[7];
//   TH1D *TrrphiG[7],*TrzG[7];
   for(i=0;i<7;i++){
       sprintf(hid,"TrxL%1.1d",i+1);
       TrxL[i] = TrxrzL[i]->ProjectionX(hid,0,TrxrzL[i]->GetNbinsY()+1,"E");
       TrxL[i]->SetXTitle("xi-x(i) (cm)");
       sprintf(hid,"TrzL%1.1d",i+1);
       TrzL[i] = TrxrzL[i]->ProjectionY(hid,0,TrxrzL[i]->GetNbinsX()+1,"E");
       TrzL[i]->SetXTitle("zi-z(i) (cm)");
//
//       sprintf(hid,"TrrphiG%1.1d",i+1);
//       TrrphiG[i] = TrrprzG[i]->ProjectionX(hid,0,
//					      TrrprzG[i]->GetNbinsY()+1,"E");
//       TrrphiG[i]->SetXTitle("rPhii-rPhi(i) (cm)");
//       sprintf(hid,"TrzG%1.1d",i+1);
//       TrzG[i] = TrrprzG[i]->ProjectionY(hid,0,
//					   TrrprzG[i]->GetNbinsX()+1,"E");
//       TrzG[i]->SetXTitle("zGi-zG(i) (cm)");
   } // end for i


   Hfile->Write();

   if(printit){
       TCanvas *c0 = new TCanvas("c0","Track quality distribution",
				 500,100,600,700);
       Ttqall->Draw();
       sprintf(filename,"%s_T_tqall.ps",sfile);
       if(printit) c0->Print(filename);
       Ttqp->Draw("COL");
       sprintf(filename,"%s_T_tq_p.ps",sfile);
       if(printit) c0->Print(filename);
       Ttq->Draw();
       sprintf(filename,"%s_T_tq.ps",sfile);
       if(printit) c0->Print(filename);
       Tp->Draw();
       sprintf(filename,"%s_T_p.ps",sfile);
       if(printit) c0->Print(filename);
       Tdttq->Draw("COL");
       sprintf(filename,"%s_T_dt_tq.ps",sfile);
       if(printit) c0->Print(filename);
       Tdt->Draw();
       sprintf(filename,"%s_T_dt.ps",sfile);
       if(printit) c0->Print(filename);
       Tadtq->Draw("COL");
       sprintf(filename,"%s_T_ad_tq.ps",sfile);
       if(printit) c0->Print(filename);
       Tad->Draw();
       sprintf(filename,"%s_T_ad.ps",sfile);
       if(printit) c0->Print(filename);
       Tadp->Draw("COL");
       sprintf(filename,"%s_T_ad_p.ps",sfile);
       if(printit) c0->Print(filename);
       Taddt->Draw("COL");
       sprintf(filename,"%s_T_ad_dt.ps",sfile);
       if(printit) c0->Print(filename);
       for(i=0;i<7;i++){
	   TrxrzL[i]->Draw("COL");
	   sprintf(filename,"%s_T%1.1d_rx_rz.ps",sfile,i+1);
	   if(printit) c0->Print(filename);
	   TrxL[i]->Draw();
	   sprintf(filename,"%s_T%1.1d_rx.ps",sfile,i+1);
	   if(printit) c0->Print(filename);
	   TrzL[i]->Draw();
	   sprintf(filename,"%s_T%1.1d_rz.ps",sfile,i+1);
	   if(printit) c0->Print(filename);
       } // end for i
       for(i=0;i<10;i++){
	   Tptqp[i]->Draw();
	   sprintf(filename,"%s_T_tq_p%1.1d.ps",sfile,i);
	   if(printit) c0->Print(filename);
       } // end for i
   } // end if printit

// Delet allocated stuff.
   for(i=0;i<7;i++) {
       delete TrxL[i];
       delete TrzL[i];
//       delete TrrphiG[i];
//       delete TrzG[i];
   } // end for i
   for(i=0;i<10;i++) delete Tptqp[i];
   delete Ttqp;
   delete Tdttq;
   delete Tadtq;
   delete Taddt;
   delete Tadp;
   for(i=0;i<7;i++) delete TrxrzL[i];
//   for(i=0;i<7;i++) delete TrrprzG[i];
   delete Ttqall;
   delete Ttq;
   delete Tp;
   delete Tdt;
   delete Tad;
//   printf("finished with track fitting\n");
   delete Sad;
   delete SrxzL;
   for(i=0;i<6;i++) delete SrxL[i];
   for(i=0;i<6;i++) delete SrzL[i];
   return;
}
//______________________________________________________________________
void FitAllTracksG(ClustAl_tl *trk,Int_t ntrk,Float_t *v0,AliITSgeom *gm,
		  const char *sfile,TFile *Hfile){
   // Local Variables
   Int_t    i,j,k,lay,lad,det;
   Double_t xm,ym,zm,dt,ad,rh,rp,rphih,rphip;
   Double_t tp,tx,rx,rz/*,rr*/,rrphi;
   Double_t xg[3],xl[3],xp,yp,xh,yh,a,b,c,d,x1[3],x2[3],qual;
   Double_t trkqualMAX = 50.0;
   Bool_t   printit = kFALSE;
   char     filename[80],hid[10];
   Int_t    Nqualmax = 10,NqualgtMAX=0;
   Int_t    Nqualless[10];
// Distrobution statisitics
   AliITSstatistics *Sad    = new AliITSstatistics(4);
   AliITSstatistics *SrxzL  = new AliITSstatistics(4);
   AliITSstatistics *SrxL[6];
   AliITSstatistics *SrzL[6];
   for(i=0;i<6;i++) {
       SrxL[i] = new AliITSstatistics(4);
       SrzL[i] = new AliITSstatistics(4);
   } // end for i

   for(i=0;i<Nqualmax;i++) Nqualless[i] = 0;

// Setup histograms

   TH2F *Ttqp  = new TH2F("Ttqp","Track quality vs. momentum",
			 500,0.0,50.0,500,0.0,10.0);
   Ttqp->Sumw2();
   Ttqp->SetXTitle("Chi^2/degree freedom");
   Ttqp->SetYTitle("Momentum GeV/c");
   TH2F *Tdttq = new TH2F("Tdttq",
			  "Distance from true vertex vs. track quality",
			  500,0.0,0.10,100,0.0,50.0);
   Tdttq->Sumw2();
   Tdttq->SetXTitle("Distanct to true vertex (cm)");
   Tdttq->SetYTitle("Chi^2/degree freedom");
   TH2F *Tadtq = new TH2F("Tadtq",
                     "atan(sqrt(b*b+d*d))-atan(pt/|pz|) vs. track quality",
			  500,-0.03,0.03,500,0.0,50.0);
   Tadtq->Sumw2();
   Tadtq->SetXTitle("theta fitted line - atan(pt/p) (rad)");
   Tadtq->SetYTitle("Chi^2/degree freedom");
   TH2F *Taddt = new TH2F("Taddt",
         "atan(sqrt(b*b+d*d))-atan(pt/|pz|) vs. Distance from true vertex",
			  500,-0.03,0.03,200,0.0,0.10);
   Taddt->Sumw2();
   Taddt->SetXTitle("theta fitted line - atan(pt/p) (rad)");
   Taddt->SetYTitle("distance to true vertex (cm)");
   TH2F *Tadp = new TH2F("Tadp",
         "atan(sqrt(b*b+d*d))-atan(pt/|pz|) vs. momentum",
			  500,-0.03,0.03,200,0.0,10.0);
   Tadp->Sumw2();
   Tadp->SetXTitle("theta fitted line - atan(pt/p) (rad)");
   Tadp->SetYTitle("Momentum (GeV/c)");

   TH2F *TrxrzL[7],*TrrprzG[7];
   for(i=0;i<7;i++){
       sprintf(filename,"Layer %1.1d: xi-x(i) vs. zi-z(i) local",i+1);
       if(i==6) sprintf(filename,"Sum of all layers: xi-x(i) "
			"vs. zi-z(i) local");
       sprintf(hid,"TrxrzL%1.1d",i+1);
       TrxrzL[i] = new TH2F(hid,filename,500,-1.0,1.0,500,-1.0,1.0);
       TrxrzL[i]->Sumw2();
       TrxrzL[i]->SetXTitle("xi-x(i) (cm) local");
       TrxrzL[i]->SetYTitle("zi-z(i) (cm) local");
//
       sprintf(filename,"Layer %1.1d: rPhii-rPhi(i) vs. zi-z(i) local",i+1);
       if(i==6) sprintf(filename,"Sum of all layers: rPhii-rPhi(i) "
			"vs. zi-z(i) local");
       sprintf(hid,"TrrprzG%1.1d",i+1);
       TrrprzG[i] = new TH2F(hid,filename,500,-1.0,1.0,500,-1.0,1.0);
       TrrprzG[i]->Sumw2();
       TrrprzG[i]->SetXTitle("rPhii-rPhi(i) (cm) local");
       TrrprzG[i]->SetYTitle("zi-z(i) (cm) local");
   } // end for i

   TH1D *Tptqp[10];
   TH1D *Ttqall = new TH1D("Ttqall","Track quality for all tracks fit",
			   500,0.0,1000.0);
   Ttqall->Sumw2();
   Ttqall->SetXTitle("Track quality: Chi squared per degree freedom");

// Fit each track and fill histograms and the like.

   for(i=0;i<ntrk;i++){
       if(FitTrackToLine(trk[i])!=0) continue;  // fit track to a line
      qual = trk[i].qual;
      if(qual<0.0) continue;
      // Fill histograms and statistics before cource chi squared cut
      Ttqall->Fill(qual,1.0);
      for(j=0;j<Nqualmax;j++) if(qual<(Double_t) j+1) Nqualless[j]++;
      if(qual>trkqualMAX) {
	  NqualgtMAX++;
	  continue;
      } // end if iqual>trkqualMAX=50.0
      a = trk[i].a;
      b = trk[i].b;
      c = trk[i].c;
      d = trk[i].d;
      zm  = b*(v0[0]-a) + d*(v0[1]-c);
      dt  = (b*b + d*d);
      if(dt!=0.0) zm /= dt;
      else{
	  printf("FitAllTracks: trk[%d].b=%f trk[%d].d=%f\n",i,b,i,d);
	  zm=0.0;
      }// end if else dt!=0.0
      xm  = a + b*zm;
      ym  = c + d*zm;
      dt   = sqrt((xm-v0[0])*(xm-v0[0]) + (ym-v0[1])*(ym-v0[1]) +
		 (zm-v0[2])*(zm-v0[2]));
      tx  = atan(sqrt(b*b + d*d));
      tp  = atan2(trk[i].pt,fabs(trk[i].pz));
      ad  = tx-tp;
      // Fill histograms and statistics
      Tadtq->Fill(ad,qual,1.0);
      Sad->AddValue(ad,1.0);
      Taddt->Fill(ad,dt,1.0);
      Tadp->Fill(ad,trk[i].p,1.0);
      Tdttq->Fill(dt,qual,1.0);
      Ttqp->Fill(qual,trk[i].p,1.0);
      for(j=0;j<trk[i].nclust;j++){
	  lay = trk[i].clust[j].lay;
	  lad = trk[i].clust[j].lad;
	  det = trk[i].clust[j].det;
	  xh = trk[i].clust[j].xg;
	  yh = trk[i].clust[j].yg;
	  rh    = hypot(xh,yh);
	  rphih = atan2(yh,xh);
	  xp = a+b*trk[i].clust[j].zg;
	  yp = c+d*trk[i].clust[j].zg;
	  rp    = hypot(xp,yp);
	  rphip = atan2(yp,xp);
	  //
	  if(fabs(rphih-rphip)>TMath::Pi()){
	      if(rphih<rphip) rphih += 2.0*TMath::Pi();
	      else rphip = 2.0*TMath::Pi();
	  } //
	  rrphi = rh*rphih-rp*rphip;
	  //
	  x1[0] = a; x1[1] = b; x1[2] = 0.0;
	  gm->GtoL(lay,lad,det,x1,xl);
	  x1[0] = xl[0]; x1[1] = xl[1]; x1[2] = xl[2];
	  if(trk[i].clust[j].zg!=0.0){
	      x2[0] = xp; x2[1] = yp; x2[2] = trk[i].clust[j].zg;
	  }else{
	      x2[0] = a+b; x2[1] = c+d; x2[2] = 1.0;
	  } // end if trk[i].clst[j].zg!=0
	  gm->GtoL(lay,lad,det,x2,xl);
	  x2[0] = xl[0]; x2[1] = xl[1]; x2[2] = xl[2];
	  xg[0] = trk[i].clust[j].xg;
	  xg[1] = trk[i].clust[j].yg;
	  xg[2] = trk[i].clust[j].zg;
	  xl[0] = trk[i].clust[j].xl;
	  xl[1] = trk[i].clust[j].yl;
	  xl[2] = trk[i].clust[j].zl;
	  if(x1[1]!=x2[1]){
	      rx = x1[0] - x1[1]*(x1[0]-x2[0])/(x1[1]-x2[1]);
	      rz = x1[2] - x1[1]*(x1[2]-x2[2])/(x1[1]-x2[1]);
	  }else{
	      continue;
	  } // end if x1[1]!=x2[1]
	  rx = rx - xl[0];
	  rz = rz - xl[2];
	  // Fill histograms and statistics
	  TrxrzL[lay-1]->Fill(rx,rz,1.0);
	  TrxrzL[6]->Fill(rx,rz,1.0);
	  TrrprzG[lay-1]->Fill(rrphi,rz,1.0);
	  TrrprzG[6]->Fill(rrphi,rz,1.0);
	  SrxL[lay-1]->AddValue(rx,1.0);
	  SrzL[lay-1]->AddValue(rz,1.0);
	  SrxzL->AddValue(rx,1.0);
	  SrxzL->AddValue(rz,1.0);
      } // end for j
   } // end for i

// Write out information

   printf("FitAllTracks: %d tracks cut leaving %d, with a chi squared >%f\n",
	  NqualgtMAX,Sad->GetN(),trkqualMAX);
   printf("The number of tracks with chi squared <1");
   for(i=1;i<Nqualmax;i++) printf(",%d",i+1);
   printf("\nFitAllTracks: %d",Nqualless[0]);
   for(i=1;i<Nqualmax;i++) printf(",%d",Nqualless[i]);
   printf("\n");
   printf("FitAllTracks: RMSs of ad, rxz = %f+-%f %f+-%f\n",
	  Sad->GetRMS(),Sad->GetErrorRMS(),
	  SrxzL->GetRMS(),SrxzL->GetErrorRMS());

   printf("FitAllTracks: Residuals by layer x=-rphi z ");
   for(i=0;i<6;i++) printf("%f+-%f %f+-%f ",
			   SrxL[i]->GetRMS(),SrxL[i]->GetErrorRMS(),
			   SrzL[i]->GetRMS(),SrzL[i]->GetErrorRMS());
   printf("\n");

// Setup and fill projections

   for(i=0;i<10;i++){
       xm = 0.5*Double_t(i);
       ym = xm+0.5;
       j  = Ttqp->GetYaxis()->FindBin(xm);
       k  = Ttqp->GetYaxis()->FindBin(ym);
       sprintf(filename,"Track Quality for %3.1f<p<%3.1f",xm,ym);
       Tptqp[i] = Ttqp->ProjectionX(filename,j,k,"E");
       Tptqp[i]->SetXTitle("Track Quality (Chi squared/df)");
   } // end for i

   TH1D *Ttq  = Ttqp->ProjectionX("Ttq",0,Ttqp->GetNbinsY()+1,"E");
   Ttq->SetXTitle("Chi^2/degree freedom");
   TH1D *Tp   = Ttqp->ProjectionY("Tp", 0,Ttqp->GetNbinsX()+1,"E");
   Tp->SetXTitle("Momentum GeV/c");
   TH1D *Tdt  = Tdttq->ProjectionX("Tdt",0,Tdttq->GetNbinsY()+1,"E");
   Tdt->SetXTitle("Distanct to true vertex (cm)");
   TH1D *Tad  = Tadtq->ProjectionX("Tad",0,Tadtq->GetNbinsY()+1,"E");
   Tad->SetXTitle("theta fitted line - atan(pt/p) (rad)");
   TH1D *TrxL[7],*TrzL[7],*TrrphiG[7],*TrzG[7];
   for(i=0;i<7;i++){
       sprintf(hid,"TrxL%1.1d",i+1);
       TrxL[i] = TrxrzL[i]->ProjectionX(hid,0,TrxrzL[i]->GetNbinsY()+1,"E");
       TrxL[i]->SetXTitle("xi-x(i) (cm)");
       sprintf(hid,"TrzL%1.1d",i+1);
       TrzL[i] = TrxrzL[i]->ProjectionY(hid,0,TrxrzL[i]->GetNbinsX()+1,"E");
       TrzL[i]->SetXTitle("zi-z(i) (cm)");
//
       sprintf(hid,"TrrphiG%1.1d",i+1);
       TrrphiG[i] = TrrprzG[i]->ProjectionX(hid,0,
					      TrrprzG[i]->GetNbinsY()+1,"E");
       TrrphiG[i]->SetXTitle("rPhii-rPhi(i) (cm)");
       sprintf(hid,"TrzG%1.1d",i+1);
       TrzG[i] = TrrprzG[i]->ProjectionY(hid,0,
					   TrrprzG[i]->GetNbinsX()+1,"E");
       TrzG[i]->SetXTitle("zGi-zG(i) (cm)");
   } // end for i


   Hfile->Write();

   if(printit){
       TCanvas *c0 = new TCanvas("c0","Track quality distribution",
				 500,100,600,700);
       Ttqall->Draw();
       sprintf(filename,"%s_T_tqall.ps",sfile);
       if(printit) c0->Print(filename);
       Ttqp->Draw("COL");
       sprintf(filename,"%s_T_tq_p.ps",sfile);
       if(printit) c0->Print(filename);
       Ttq->Draw();
       sprintf(filename,"%s_T_tq.ps",sfile);
       if(printit) c0->Print(filename);
       Tp->Draw();
       sprintf(filename,"%s_T_p.ps",sfile);
       if(printit) c0->Print(filename);
       Tdttq->Draw("COL");
       sprintf(filename,"%s_T_dt_tq.ps",sfile);
       if(printit) c0->Print(filename);
       Tdt->Draw();
       sprintf(filename,"%s_T_dt.ps",sfile);
       if(printit) c0->Print(filename);
       Tadtq->Draw("COL");
       sprintf(filename,"%s_T_ad_tq.ps",sfile);
       if(printit) c0->Print(filename);
       Tad->Draw();
       sprintf(filename,"%s_T_ad.ps",sfile);
       if(printit) c0->Print(filename);
       Tadp->Draw("COL");
       sprintf(filename,"%s_T_ad_p.ps",sfile);
       if(printit) c0->Print(filename);
       Taddt->Draw("COL");
       sprintf(filename,"%s_T_ad_dt.ps",sfile);
       if(printit) c0->Print(filename);
       for(i=0;i<7;i++){
	   TrxrzL[i]->Draw("COL");
	   sprintf(filename,"%s_T%1.1d_rx_rz.ps",sfile,i+1);
	   if(printit) c0->Print(filename);
	   TrxL[i]->Draw();
	   sprintf(filename,"%s_T%1.1d_rx.ps",sfile,i+1);
	   if(printit) c0->Print(filename);
	   TrzL[i]->Draw();
	   sprintf(filename,"%s_T%1.1d_rz.ps",sfile,i+1);
	   if(printit) c0->Print(filename);
       } // end for i
       for(i=0;i<10;i++){
	   Tptqp[i]->Draw();
	   sprintf(filename,"%s_T_tq_p%1.1d.ps",sfile,i);
	   if(printit) c0->Print(filename);
       } // end for i
   } // end if printit

// Delet allocated stuff.
   for(i=0;i<7;i++) {
       delete TrxL[i];
       delete TrzL[i];
       delete TrrphiG[i];
       delete TrzG[i];
   } // end for i
   for(i=0;i<10;i++) delete Tptqp[i];
   delete Ttqp;
   delete Tdttq;
   delete Tadtq;
   delete Taddt;
   delete Tadp;
   for(i=0;i<7;i++) delete TrxrzL[i];
   for(i=0;i<7;i++) delete TrrprzG[i];
   delete Ttqall;
   delete Ttq;
   delete Tp;
   delete Tdt;
   delete Tad;
//   printf("finished with track fitting\n");
   delete Sad;
   delete SrxzL;
   for(i=0;i<6;i++) delete SrxL[i];
   for(i=0;i<6;i++) delete SrzL[i];
   return;
}
//______________________________________________________________________
void FindVertex2(ClustAl_tl &trk1,ClustAl_tl &trk2,Double_t *vt,Double_t &d){
   // Local Variables
   Double_t x1,y1,x2,y2;
   Double_t z1,z2;
   Double_t a1,b1,c1,d1;
   Double_t a2,b2,c2,d2;
   Double_t da,db,dc,dd,dbd;
   Double_t den,num1,num2;

   a1 = trk1.a; b1 = trk1.b; c1 = trk1.c; d1 = trk1.c;
   a2 = trk2.a; b2 = trk2.b; c2 = trk2.c; d2 = trk2.c;

   // Find z1 and z2 of points of closest approch.
   da = a1-a2;
   db = b1-b2;
   dc = c1-c2;
   dd = d1-d2;
   dbd = b1*d2-d1*b2;
   den = -db*db - dd*dd - dbd*dbd;
   num1 = da*(db+d1*dbd) + dc*(dd-b1*dbd);
   num2 = da*(db+d2*dbd) + dc*(dd-b2*dbd);
   if(den!=0.0){
       z1 = num1/den;
       z2 = num2/den;
   }else{ // parallel lines
       z1 = 0.0;
       z2 = 0.0;
   } // end if den!-0.0

   // find coordinate of closest approch and distance between.
   x1    = a1+b1*z1; y1 = c1+d1*z1;
   x2    = a2+b2*z2; y2 = c2*d2*z2;
   vt[0] = 0.5*(x1+x2);
   vt[1] = 0.5*(y1+y2);
   vt[2] = 0.5*(z1+z2);
   d     = sqrt( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) );

   return;
}
//______________________________________________________________________
void FitVertexAll(ClustAl_tl *trk,Int_t ntrk,const char *sfile,TFile *Hfile){
   // Local Variables
   Int_t    i,j,k;
   Double_t d,vz,r,q;
   Double_t vt[3];
   Bool_t   printit = kFALSE;
   char     filename[80],hid[10];
   Double_t trkqualMAX = 50.0;
   const Int_t    IntQualMax=10;
   Int_t    QualInt[IntQualMax];
   AliITSstatistics *Svz   = new AliITSstatistics(4);
   AliITSstatistics *Svx   = new AliITSstatistics(4);
   AliITSstatistics *Svy   = new AliITSstatistics(4);
   //   Double_t vzmi=-0.20,vzma=0.20,vxmi=-1.0,vxma=1.0,vymi=-1.0,vyma=1.0;

   for(i=0;i<IntQualMax;i++) QualInt[i] = 0;

   TH1F *Vvzbtq[IntQualMax];
   for(i=0;i<IntQualMax;i++){
       sprintf(hid,"Vvzbta%2.2d",i+1);
       sprintf(filename,"Z of vertex of pairs for tracks with chi squared <%d",
	       i+1);
       Vvzbtq[i] = new TH1F(hid,filename,100,-1.0,1.0);
       Vvzbtq[i]->Sumw2();
       Vvzbtq[i]->SetXTitle("Z of vertex (cm)");
   } // end for i

   TH2F *Vvztq  = new TH2F("Vvztq","Z of vertex of pairs vs. track quality cut",
			   500,-0.20,0.20,200,0.0,20.0);
   Vvztq->Sumw2();
   Vvztq->SetXTitle("Z of vertex (cm)");
   Vvztq->SetYTitle("Chi^2/degree freedom");
   TH2F *Vvzvr  = new TH2F("Vvzvr","Z vs. R of vertex of pairs",
			                       200,-0.20,0.20,400,0.0,5.0);
   Vvzvr->Sumw2();
   Vvzvr->SetXTitle("Z of vertex (cm)");
   Vvzvr->SetYTitle("R of vertex (cm)");
   TH2F *Vdtq   = new TH2F("Vdtq","Distance between lines vs. track quality cut",
			   500,0.0,0.2,200,0.0,50.0);
   Vdtq->Sumw2();
   Vdtq->SetXTitle("minimum distance between lines (cm)");
   Vdtq->SetYTitle("Chi^2/degree freedom");
   TH2F *Vvzd   = new TH2F("Vvzd","Z vertex vs. Dist between lines",
                                               500,-5.0,5.0,200,0.0,5.0);
   Vvzd->Sumw2();
   Vvzd->SetXTitle("Z of vertex (cm)");
   Vvzd->SetYTitle("minimmum distance between lined (cm)");
   TH2F *Vvxvy  = new TH2F("Vvxvy","X vertex vs. Y vertex",
                                               200,-10.0,10.0,200,-10.0,10.0);
   Vvxvy->Sumw2();
   Vvxvy->SetXTitle("X of vertex (cm)");
   Vvxvy->SetYTitle("Y of vertex (cm)");
   TH2F *Vvxvz  = new TH2F("Vvxvz","X vertex vs. Z vertex",
                                               200,-10.0,10.0,200,-10.0,10.0);
   Vvxvz->Sumw2();
   Vvxvz->SetXTitle("X of vertex (cm)");
   Vvxvz->SetYTitle("Z of vertex (cm)");
   TH2F *Vvyvz  = new TH2F("Vvyvz","Y vertex vs. Z vertex",
                                               200,-10.0,10.0,200,-10.0,10.0);
   Vvyvz->Sumw2();
   Vvyvz->SetXTitle("Y of vertex (cm)");
   Vvyvz->SetYTitle("Z of vertex (cm)");

   for(i=0;i<ntrk-1;i++){
       if(trk[i].qual>trkqualMAX||trk[i].qual<0.0) continue;
       for(j=i+1;j<ntrk;j++){
	   if(trk[j].qual>trkqualMAX||trk[i].qual<0.0) continue;
	   FindVertex2(trk[i],trk[j],vt,d);
	   q  = TMath::Max(trk[i].qual,trk[j].qual);
	   r  = vt[0]*vt[0] + vt[1]*vt[1];
	   r  = sqrt(r);
	   vz = vt[2];
           Vvxvy->Fill(vt[0],vt[1],1.0);
           Vvxvz->Fill(vt[0],vt[2],1.0);
           Vvyvz->Fill(vt[1],vt[2],1.0);
	   Svx->AddValue(vt[0],1.0);
	   Svy->AddValue(vt[1],1.0);
	   Svz->AddValue(vt[2],1.0);
	   Vvztq->Fill(vz,q,1.0);
	   Vvzvr->Fill(vz,r,1.0);
	   Vdtq->Fill(d,q,1.0);
	   Vvzd->Fill(vz,d,1.0);
	   for(k=0;k<IntQualMax;k++)if(q<(Double_t)k+1.)QualInt[k]++;
	   for(k=0;k<IntQualMax;k++)if(q<(Double_t)k+1.)Vvzbtq[k]->Fill(vz,1.);
       } // end for j
   } // end for i

   TH1D *Vd  = Vvzd->ProjectionY ("Vd", 0,Vvzd->GetNbinsX()+1,"E");
   Vd->SetXTitle("minimum distance between lines (cm)");
   TH1D *Vvr = Vvzvr->ProjectionY("Vvr",0,Vvzvr->GetNbinsX()+1,"E");
   Vvr->SetXTitle("R of vertex (cm)");
   TH1D *Vtq = Vdtq->ProjectionY ("Vtq",0,Vdtq->GetNbinsX()+1,"E");
   Vtq->SetXTitle("Chi^2/degree freedom");
//
   TH1D *Vvx = Vvxvy->ProjectionX("Vvx",0,Vvxvy->GetNbinsY()+1,"E");
   Vvx->SetXTitle("X of vertex (cm)");
   TH1D *Vvy = Vvxvy->ProjectionY("Vvy",0,Vvxvy->GetNbinsX()+1,"E");
   Vvy->SetXTitle("Y of vertex (cm)");
   TH1D *Vvz = Vvxvz->ProjectionX ("Vvz",0,Vvxvz->GetNbinsY()+1,"E");
   Vvz->SetXTitle("Z of vertex (cm)");

   printf("FitVertexAll: N(track pairs Qual <=1");
   for(k=1;k<IntQualMax;k++) printf(",%d",k+1);
   printf("\nFitVertexAll: %d",QualInt[0]);
   for(k=1;k<IntQualMax;k++) printf(",%d",QualInt[k]);
   printf("FitVertexAll: RMS of vx vy vz=%e+-%e %e+-%e %e+-%e\n",
	  Svx->GetRMS(),Svx->GetErrorRMS(),
	  Svy->GetRMS(),Svy->GetErrorRMS(),
	  Svz->GetRMS(),Svy->GetErrorRMS());
   Hfile->Write();

   if(printit){
       TCanvas *c1 = new TCanvas("c1","Vertex info",500,100,600,700);
       for(k=0;k<IntQualMax;k++){
	   Vvzbtq[k]->Draw();
	   sprintf(filename,"%s_V_vz_tq%2.2d.ps",sfile,k+1);
	   if(printit) c1->Print(filename);
       } // end for k
       Vvzd->Draw("COL");
       sprintf(filename,"%s_V_vz_d.ps",sfile);
       if(printit) c1->Print(filename);
       Vvzvr->Draw("COL");
       sprintf(filename,"%s_V_vz_vr.ps",sfile);
       if(printit) c1->Print(filename);
       Vdtq->Draw("COL");
       sprintf(filename,"%s_V_d_tq.ps",sfile);
       if(printit) c1->Print(filename);
       Vvztq->Draw("COL");
       sprintf(filename,"%s_V_vz_tq.ps",sfile);
       if(printit) c1->Print(filename);
       Vvz->Draw();
       sprintf(filename,"%s_V_vz.ps",sfile);
       if(printit) c1->Print(filename);
       Vd->Draw();
       sprintf(filename,"%s_V_d.ps",sfile);
       if(printit) c1->Print(filename);
       Vvr->Draw();
       sprintf(filename,"%s_V_vr.ps",sfile);
       if(printit) c1->Print(filename);
       Vtq->Draw();
       sprintf(filename,"%s_V_tq.ps",sfile);
       if(printit) c1->Print(filename);
       Vvxvy->Draw("COL");
       sprintf(filename,"%s_V_vx_vy.ps",sfile);
       if(printit) c1->Print(filename);
       Vvxvz->Draw("COL");
       sprintf(filename,"%s_V_vx_vz.ps",sfile);
       if(printit) c1->Print(filename);
       Vvyvz->Draw("COL");
       sprintf(filename,"%s_V_vy_vz.ps",sfile);
       if(printit) c1->Print(filename);
       Vvx->Draw();
       sprintf(filename,"%s_V_vx.ps",sfile);
       if(printit) c1->Print(filename);
       Vvy->Draw();
       sprintf(filename,"%s_V_vy.ps",sfile);
       if(printit) c1->Print(filename);
   } // end if printit
   for(k=0;k<IntQualMax;k++) delete Vvzbtq[k];
   delete Vvzd;
   delete Vvzvr;
   delete Vdtq;
   delete Vvztq;
   delete Vvz;
   delete Vd;
   delete Vvr;
   delete Vtq;
   delete Vvxvy;
   delete Vvxvz;
   delete Vvyvz;
   delete Vvx;
   delete Vvy;
   delete Svz;
   delete Svx;
   delete Svy;

   return;
}
//______________________________________________________________________
void OnlyOneGeometry(char *filename,AliITSgeom *gm,AliITSgeom &gm2,
		     Float_t trans[],Float_t rot[]){

    ifstream f_in(filename);
    if(f_in){
	printf("reading from geometry file %s\n",filename);
	gm2.ReadGeom(f_in);
	f_in.close();
    }else{
	f_in.close();
	gm2 = *gm;
	gm2.RandomCylindericalChange(trans,rot);
	printf("writting to geometry file %s\n",filename);
	ofstream f_out(filename);
	gm2.PrintGeom(f_out);
	f_out.close();
    } // end if
    return;
}
//______________________________________________________________________
void deleteClustAl(ClustAl_tl *trk,Int_t ntrk){

    Int_t i;

    for(i=0;i<ntrk;i++) delete[] trk[i].clust;
    return;
}
//______________________________________________________________________
void RunAlignment(Int_t evnt,Float_t fraction=1.0){
// define some variables for later use.
//
    printf("gAlice=%p and evnt=%d\n",gAlice,evnt);
//
    Int_t nparticles = gAlice->GetEvent(evnt);
    printf("nparticles %d\n",nparticles);
    if (nparticles <= 0) return; /* get next event */

// Get pointers to Alice detectors and Clusts containers
    AliITS *ITS = (AliITS*)gAlice->GetDetector("ITS");
    if(!ITS) return;          /* error no ITS data exit */
    TTree *TH   = gAlice->TreeH();
    Int_t Ntrkp =  (Int_t) TH->GetEntries(),ntrk;
    AliITSgeom gm2,gm3,*gm = ITS->GetITSgeom();
    Int_t Nmods = gm->GetIndexMax();
    Float_t   trans[15] ={0.0E-0,1.0E-4,4.0E-4,7.0E-4,1.0E-3,
			  2.0E-3,4.0E-3,6.0E-3,8.0E-3,1.0E-2,
			  2.0E-2,3.0E-2,5.0E-2,7.5E-2,1.0E-1}; // cm
    Float_t   rots[15] ={0.0E-0,1.0E-4,4.0E-4,7.0E-4,1.0E-3,
			 2.0E-3,4.0E-3,6.0E-3,8.0E-3,1.0E-2,
			 2.0E-2,3.0E-2,5.0E-2,7.5E-2,1.0E-1}; // rad
   Float_t     tran[3] = {0.0,0.0,0.0},rot[3] = {0.0,0.0,0.0};
   TFile *Hfile;
   Float_t     Rdta[6];
   Int_t Itimes=0,i,j,badmod;
   Double_t Chi2b;

// Array (stucture) of clusts for the first and second layer
// this should be replaced with either clusters or digits
// when they are proporly defined.
    AliITSAlignmentTrack *trk = new AliITSAlignmentTrack[Ntrkp];
    TObjArray *mods = new TObjArray(Nmods);
    AliITSAlignmentModule *mod;

    printf("Ntrkp=%d\n",Ntrkp);

    FillAliITSAlignmentTrack(trk,ntrk,Ntrkp,TH,ITS,fraction);

    for(Int_t Isigmas=0;Isigmas<1;Isigmas++){
//
//      tran[0] = sigma1;
//	tran[1] = sigma2;
//	tran[2] = sigma3;
	if(Itimes==0){ tran[0] = trans[Isigmas];
	}else tran[0] = 0.0;
	if(Itimes==1){ tran[1] = trans[Isigmas];
	}else tran[1] = 0.0;
	if(Itimes==2){ tran[2] = trans[Isigmas];
	}else tran[2] = 0.0;
	if(Itimes==3){ rot[0] = rots[Isigmas];
	}else rot[0] = 0.0;
	if(Itimes==4){ rot[1] = rots[Isigmas];
	}else rot[1] = 0.0;
	if(Itimes==5){ rot[2] = rots[Isigmas];
	}else rot[2] = 0.0;
	printf("tran= %e %e %e (cm), rot=%e %e %e (rad)\n",
	       tran[0],tran[1],tran[2],rot[0],rot[1],rot[2]);
//
	gm2 = *gm;
	gm2.RandomCylindericalChange(tran,rot);
	gm3 = gm2;
	Hfile = new TFile("Alignment_geom_0_2.root","RECREATE",
			  "comparison of geometry before refitting");
        PlotGeomChanges(gm,&gm2,Hfile,Rdta);
	Hfile -> Close();

	for(i=0;i<Nmods;i++) mods->AddAt(new AliITSAlignmentModule(i,&gm3,
                                                                  ntrk,trk),i);
//
	j = 0;
	do{
	    for(i=0;i<ntrk;i++) trk[i].SetGlobalPosition(&gm2);
	    //
	    for(i=0;i<ntrk;i++){
		trk[i].FitToFunction(2,&gm3);
	    } // end for i
	    // find detector with the worst Chi2.
	    j++;
	    Chi2b  = 0.0;
	    badmod = 0;
	    for(i=0;i<Nmods;i++){
		if(((AliITSAlignmentModule *)(mods->At(i)))->GetChi2()>Chi2b){
		    Chi2b =((AliITSAlignmentModule *)(mods->At(i)))->GetChi2();
		    badmod = i;
		} // end if
	    } // end for i
	    mod = (AliITSAlignmentModule *)(mods->At(badmod));
	    mod->AlignModule();
	    gm3.SetTrans(badmod,mod->GetTranslationVector());
	    gm3.SetByAngles(badmod,mod->GetRotationAngles());
	}while(j<Nmods);
    } // end for Isigmas
//
    Hfile = new TFile("Alignment_geom_0_3.root","RECREATE",
		      "comparison of geometry after refitting");
    PlotGeomChanges(gm,&gm3,Hfile,Rdta);
    Hfile->Close();
    printf("Event %d done\n",evnt);
//
    delete[] trk;            // now delet memory allocated above.
}
//______________________________________________________________________
