#include <stdio.h>
#include <math.h>
#include <time.h>
#include "TMinuit.h"
#include "TRandom.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TProfile.h"
#include "TH2.h"
#include "TArray.h"
#include "TCanvas.h"
#include "TString.h"
#include "AliRun.h"
#include "AliITS.h"
#include "AliITSgeom.h"
#include "TParticle.h"
#include "vertex.h"

// Global in this file only. Needed for Minuit
Int_t gstrt,gend;  // starting and ending bin values
TH1D  *gH1;        // histogram to be fit.

void FindVertexs(Int_t evnt,Float_t frac=1.0,Float_t len=16.82){

      // Get pointers to Alice detectors and Hits containers
      AliDetector  *ITS       = gAlice->GetDetector("ITS");
      if(!ITS) return;          // error no ITS data exit
//      TClonesArray *Particles = gAlice->Particles();
      TClonesArray *ITShits   = ITS->Hits();
      TTree        *TH        = gAlice->TreeH();
      Int_t        ntracks    =  (Int_t)TH->GetEntries();
//      ntracks = 10;

      printf(" %d primary tracks in the file\n",ntracks);
      // Array (stucture) of hits for the first and second layer
      // this should be replaced with either clusters or digits
      // when they are proporly defined.
      Hit_tl *spd1 = new Hit_tl[ntracks];
      Hit_tl *spd2 = new Hit_tl[ntracks];
      Hit_tl **spdi = new Hit_tl*[ntracks];
      Hit_tl **spdo = new Hit_tl*[ntracks];
      for(Int_t i=0;i<ntracks;i++){
	  spdi[i] = &(spd1[i]);
//	  printf("spdi[%d]=%p spd1[%d]=%p  ",i,spdi[i],i,&(spd1[i]));
	  spdo[i] = &(spd2[i]);
//	  printf("spdo[%d]=%p spd2[%d]=%p\n",i,spdo[i],i,&(spd2[i]));
      } // end for i
      Int_t   i1max,i2max;
      Float_t vz;

      HitsToV(spdi,i1max,spdo,i2max,ntracks,TH,ITShits,frac,len);
      printf("back in Macro i1max=%d i2max=%d\n",i1max,i2max);

//      Float_t r1=0.0,r2=0.0;
//      for(i=0;i<TMath::Max(i1max,i2max);i++){
//	  if(i<i1max) r1 += spdi[i]->r;
//	  if(i<i2max) r2 += spdo[i]->r;
//      } // end for i
//      printf("<r1>=%f i1max=%d\t<r2>=%f i2max=%d\n",
//	     r1/Float_t(i1max),i1max,r2/Float_t(i2max),i2max);

      vz = vertexSlow(spdi,i1max,spdo,i2max);

      printf("Slow Sorted event=%d Zvertex=%f\n",evnt,vz);

//      vz = vertex(spdi,i1max,spdo,i2max);

//      printf("Phi sorted event=%d Zvertex=%f\n",evnt,vz);

//      vz = vertexEta(spdi,i1max,spdo,i2max);

//      printf("Eta sorted event=%d Zvertex=%f\n",evnt,vz);

      return;
}


Bool_t L_SortPhi(const Hit_tl *s1,const Hit_tl *s2){
  // Phi sorting function for qsort.
   Float_t a;

   a = s1->phir - s2->phir;
   if(a<0.0) return kTRUE;
   if(a>0.0) return kFALSE;
   a = s1->etar - s2->etar;
   if(a<0.0) return kTRUE;
   if(a>0.0) return kFALSE;
   return kFALSE;
}

Bool_t L_SortEta(const Hit_tl *s1,const Hit_tl *s2){
  // Eta sorting function for qsort.
   Float_t a;

   a = s1->etar - s2->etar;
   if(a<0.0) return kTRUE;
   if(a>0.0) return kFALSE;
   a = s1->phir - s2->phir;
   if(a<0.0) return kTRUE;
   if(a>0.0) return kFALSE;
   return kFALSE;
}

void hpsortPhi(Hit_tl **ra,Int_t n){
   Int_t i,ir,j,l;
   Hit_tl *rra;

   if(n<2) return;

   l  = ((n-1) >> 1) +1; // devide 2 + 1
   ir = n-1;
   for(;;){
     if(l>0){
        rra = ra[--l];  // decrament first
     }else{
        rra    = ra[ir];
        ra[ir] = ra[0];
        if(--ir == 0){  // decrament first
           ra[0] = rra;
           break;
        } // if --ra == 0 
     } // end l>0 
     i = l;
     j = l+1;
     while(j<=ir){
        if( j<ir && L_SortPhi(ra[j],ra[j+1]) ) j++;
        if( L_SortPhi(rra,ra[j]) ){
           ra[i] = ra[j];
           i = j;
           j <<= 1; // time 2.
        }else{
           break;
        } // end if func()
     } // end while
     ra[i] = rra;
   } // end for ever 
}


void hpsortEta(Hit_tl **ra,Int_t n){
   Int_t i,ir,j,l;
   Hit_tl *rra;

   if(n<2) return;

   l  = ((n-1) >> 1) +1; // devide 2 + 1
   ir = n-1;
   for(;;){
     if(l>0){
        rra = ra[--l];  // decrament first
     }else{
        rra    = ra[ir];
        ra[ir] = ra[0];
        if(--ir == 0){  // decrament first
           ra[0] = rra;
           break;
        } // if --ra == 0 
     } // end l>0 
     i = l;
     j = l+1;
     while(j<=ir){
        if( j<ir && L_SortEta(ra[j],ra[j+1]) ) j++;
        if( L_SortEta(rra,ra[j]) ){
           ra[i] = ra[j];
           i = j;
           j <<= 1; // time 2.
        }else{
           break;
        } // end if func() 
     } // end while 
     ra[i] = rra;
   } // end for ever 
}

void fillStructure(Hit_tl **spd,Int_t &nspd,
		   Float_t *xv,Int_t *id,Int_t track,Int_t nhits){
    Float_t PI2 = 2.0*TMath::Pi();
//    Float_t PI  = TMath::Pi();
    Int_t i;
    Float_t x,y,z,zr,r,phi,eta,rphi;
    Float_t dz    = 0.0300; // 300 microns
    Float_t drphi = 0.0050; //  50 microns

    i = nspd;
//   if(nspd<5) printf("fill: spd=%p spd[%d]=%p i=%d spd[i]=%p id=%d %d %d\n",
//		      spd,nspd,spd[nspd],i,spd[i],id[0],id[1],id[2]);
    x = xv[0];
    y = xv[1];
    z = xv[2];
    r = sqrt(x*x+y*y);
    phi = atan2(y,x);
    if(phi<0.0) phi += PI2;
    rphi = r*phi;
    eta  = -0.5*tan(0.5*atan2(r,z));

    spd[i]->track = track;
    spd[i]->n     = nhits;
    spd[i]->lad   = id[1];
    spd[i]->det   = id[2];
    spd[i]->x     = x;
    spd[i]->y     = y;
    spd[i]->z     = z;
    spd[i]->r     = r;
    spd[i]->phi   = phi;
    spd[i]->eta   = eta;

    zr   = dz    * Float_t(Int_t(z   /   dz)) + 0.5*dz;
    rphi = drphi * Float_t(Int_t(rphi/drphi)) + 0.5*drphi;

    spd[i]->xr   =  cos(phi)*rphi/r;
    spd[i]->yr   = -sin(phi)*rphi/r;
    spd[i]->zr   = zr;
    spd[i]->phir = rphi/r;
    spd[i]->etar = -0.5*tan(0.5*atan2(r,zr));

    nspd++;
    return;
}

void HitsToV(Hit_tl **spdi,Int_t &nspdi,Hit_tl **spdo,Int_t &nspdo,
	    Int_t ntracks,TTree *TH,TClonesArray *ITShits,
	     Float_t fraction=1.0,Float_t len=16.82){
    Int_t     i,t,h,n,nb,nhits,id[4],idold[4],Iseed,ieta=0,itrk=0,ieta2=0;
    Float_t   x[3],xb[3],xl[3];
    AliITShit *itsHit;
    TClonesArray *Part = gAlice->Particles();
    TParticle *part;

    nspdi = nspdo = 0;

    idold[0] = idold[1] = idold[2] = idold[3] = 0;
    xb[0]    = xb[1]    = xb[2]    = 0.0;
    n = 0;
    Iseed = time(0);
    printf("HitsToV: Iseed=%d fraction=%f Pixel length=%fcm\n",
	   Iseed,fraction,len);
    gRandom->SetSeed(Iseed);
//    printf("HitsToV: gRandom->Rndm(1)=%f\n",gRandom->Rndm(1));
    for(t=0;t<ntracks;t++){
	if(fraction<gRandom->Rndm(1)) continue; // skip some tracks
	itrk++;
	gAlice->ResetHits();
	nb       = TH->GetEvent(t);
	nhits    = ITShits->GetEntriesFast();
	for(h=0;h<nhits;h++){
	    itsHit = (AliITShit *) ITShits->UncheckedAt(h);
	    itsHit->GetDetectorID(id[0],id[1],id[2]); id[3]=t;
	    itsHit->GetPositionG(x[0],x[1],x[2]);
	    itsHit->GetPositionL(xl[0],xl[1],xl[2]);
	    if(h==0){
		part = (TParticle *) (Part->UncheckedAt(itsHit->GetTrack()));
		if(TMath::Abs(part->Eta())<=1.0) ieta++;
		if(TMath::Abs(part->Eta())<=0.5) ieta2++;
	    } // end if h==0
	    if(TMath::Abs(x[2]/len) >= 1.0) continue;
	    if(x[0]==0.0&&x[1]==0.0){
		printf("Hitsto: t=%d h=%d/%d id=%d,%d,%d x=%f,%f,%f\n",
		       t,h,nhits,id[0],id[1],id[2],x[0],x[1],x[2]); 
		continue;
	    } // end if x[0]==x[1]==0.0
	    if(!(id[0]==idold[0]&&id[1]==idold[1]&&
		 id[2]==idold[2]&&id[3]==idold[3])){
		if(!(n<=0 || (xb[0]==0.0&&xb[1]==0.0))){
		    for(i=0;i<3;i++) xb[i]   /= Float_t(n);
		    if(idold[0]==1){
			fillStructure(spdi,nspdi,xb,idold,t,n);
		    }
		    if(idold[0]==2){
			fillStructure(spdo,nspdo,xb,idold,t,n);
		    }
		    if(nspdi>ntracks || nspdo>ntracks){
			printf("Hitsto: fill error,"
			       " nspdi=%d nspdo=%d ntracks=%d\n",
			       nspdi,nspdo,ntracks);
		    } // end if fill error
		} // end if n>0
		for(i=0;i<4;i++) idold[i] = id[i];
		for(i=0;i<3;i++) xb[i]    = 0.0;
		n = 0;
	    } // end if id != idold
	    for(i=0;i<3;i++) xb[i] += x[i];
	    n++;
	} // end for h
    } // end for t
    printf("exiting HitsToV: %d primary tracks in eta=+-1\n",ieta);
    printf("exiting HitsToV: %d primary tracks #2 in eta=+-0.5\n",ieta2);
    printf("exiting HitsToV: %d primary tracks in file used\n",itrk);
    return;
}

Float_t vertex(Hit_tl **spdi,Int_t i1max,Hit_tl **spdo,Int_t i2max){
   Float_t r1    = 3.910078; // radius at which hit is from beam axis for spd1
   Float_t r2    = 6.955933; // radius at which hit is from beam axis for spd2
   Float_t DPhi  = 0.005;    // maximum allowed difference in phi angle
   Float_t DZ    = 12.5;      // maximum allowed difference in z position
   Float_t avt   = 0.0;
   Int_t   nbinx = 2000;
   Int_t   start = 0;
   Bool_t  mod   = kFALSE;
   Int_t   i,j;
   Float_t zv,av,su,i0,i1,i2,x0,x1,dphi,dz;
   Float_t PI2 = 2.0*TMath::Pi();
   Float_t PI  = TMath::Pi();

   // sort according to phi.
   hpsortPhi(spdi,i1max);
   hpsortPhi(spdo,i2max);

   // find best vertex allong z.
   TH2S *Pvzphi   = new TH2S("Pvzphi","Phi: Posible Z vertecies vs. Phi",
                              nbinx,-1.0,1.0,32,0.0,PI2);
   Pvzphi->Sumw2(); // collect good statitics
   TH1F *Pvzfl   = new TH1F("Pvzfl","Phi: Posible Z vertecies flattened",
			 nbinx,-1.0,1.0);
   Pvzfl->Sumw2();
   TH1F *Pvztr  = new TH1F("Pvztr","Phi: Z Vertex found for True Tracks",
			 nbinx,-1.0,1.0);
   Pvztr->Sumw2(); // collect good statitics

   for(i=0;i<i1max;i++){
       if(spdi[i]->r==0.0) {printf("spdi[%d]->r=0.0\n",i);continue;}
       for(;spdo[start]->phir<spdi[i]->phir-1.5*DPhi;start++);
       for(j=start;(spdo[j]->phir < spdi[i]->phir+DPhi) && (j<i2max);j++){
	   dphi = fabs(spdo[j]->phir - spdi[i]->phir); if(dphi>PI) dphi -= PI;
	   dz   = fabs(spdo[j]->zr   - spdi[i]->zr);
	   if(dphi>DPhi) continue;
           if(dz>DZ)     continue; // If outside dz range skip it
	   r1   = spdi[i]->r;
	   r2   = spdo[j]->r;
	   if(r2-r1!=0.0) zv   = (r2*spdi[i]->zr - r1*spdo[j]->zr)/(r2-r1);
	   else{
	       printf("vertex: spdi[%d]->r=%f = spdo[%d]->r=%f\n",i,r1,j,r2);
	       continue;
	   } // end if else r2-r1!=0.0
	   Pvzphi->Fill(zv,spdi[i]->phir);
	   Pvzfl->Fill(zv);
	   if(spdi[i]->track == spdo[j]->track) Pvztr->Fill(zv);
       } // end for j
   } // end for i

   TH1D *Pvzdef  = Pvzphi->ProjectionX("Phi: Posible Z vertecies",1,nbinx,"E");

   i0 = Pvzfl->GetBinContent(1);
   x0 = Float_t(1);
   i1 = Pvzfl->GetBinContent(nbinx);
   x1 = Float_t(nbinx);
   su = x0-x1; if(su==0.0) return -200.0;
   su = (i0-i1)/su;
   x0 = x0-su*i0;
   x1 = su;
   for(i=1;i<=nbinx;i++){
      su = x1*Float_t(i) + x0;
      Pvzfl->AddBinContent(i,-su);
   } // end for i

   printf("Phi:            mean=%f RMS=%f w2=%f\n",
               Pvzdef->GetMean(),Pvzdef->GetRMS(),Pvzdef->GetSumOfWeights());
   printf("Phi: Flattened mean=%f RMS=%f w2=%f\n",
	          Pvzfl->GetMean(),Pvzfl->GetRMS(),Pvzfl->GetSumOfWeights());

   av = 3.0 * Pvzfl->Integral(1,nbinx)/Float_t(nbinx);
   printf("Phi: Flattened av=%f Pvzfl->Max=%f nbinx=%d\n",
	  av,Pvzfl->GetMaximum(),nbinx);

   su  = 0.0;
   avt = 0.0;
   i1  = Pvzfl->GetBinContent(1);
   i2  = Pvzfl->GetBinContent(2);
   for(i=2;i<nbinx;i++){
       i0 = i1; // old i-1 value
       i1 = i2; // old i value
       i2 = Pvzfl->GetBinContent(i+1); // new i+1 value
       if(i1 > av && (i0 > av || i2 > av ) ){
	   if(!mod) mod = kTRUE;
       } else if(mod) break;
       if(mod){ // inside peak
	   su  += i1;
	   avt += i1 * Pvzfl->GetXaxis()->GetBinCenter(i);
       } // end if mod
   } // end for i

   if(su!=0.0) zv = avt/su;
   else        zv = -100.0; // an unphysical value

   TCanvas *c0 = new TCanvas("c0","Alice ITS vertex finder", 400,10,600,700);
   Pvzphi->Draw("col");
   c0->Print("vertex5_P_vz_phi.ps");
   Pvzdef->Draw();
   c0->Print("vertex5_P_vz_def.ps");
   Pvzfl->Draw();
   c0->Print("vertex5_P_vz_flat.ps");
   Pvztr->Draw();
   c0->Print("vertex5_P_vz_true.ps");

   return zv;
}

void Chi2Gauss(Int_t &npar,Double_t *gin,Double_t &f,
	       Double_t *par,Int_t iflag){
    Double_t chi2 = 0.0;
    Double_t delta,h,x,eh;

    for(Int_t i=gstrt;i<gend;i++){
	h = gH1->GetBinContent(i);
	eh = TMath::Sqrt(h);
	if(eh <= 0.0) eh = 1.0;
	x = gH1->GetXaxis()->GetBinCenter(i);
	delta = (h - par[0]*TMath::Gaus(x,par[1],par[2]) - par[3])/eh;
	chi2 += delta*delta;
    } // end for i
    f = chi2;
    return;
}

Float_t vertexSlow(Hit_tl **spdi,Int_t i1max,Hit_tl **spdo,Int_t i2max){
   Float_t r1    = 3.910078; // radius at which hit is from beam axis for spd1
   Float_t r2    = 6.955933; // radius at which hit is from beam axis for spd2
   Float_t DPhi  = 0.005;    // maximum allowed difference in phi angle
   Float_t BDphi;
   Float_t DZ    = 12.5;      // maximum allowed difference in z position
   Float_t avt   = 0.0;
   Int_t   nbinx = 2000;
   Int_t   nbst;
   Bool_t  mod   = kFALSE;
   Int_t   i,j;
   Float_t zv,av,su,i0,i1,i2,x0,x1,dphi,dz;
   Float_t PI2 = 2.0*TMath::Pi();
   Float_t PI  = TMath::Pi();

   if(i1max<=0||i2max<=0) return -1.0E10;

   // find best vertex allong z.
   TH2S *Svzphi  = new TH2S("Svzphi","Slow: Posible Z vertecies vs. Phi",
			    nbinx,+0.0,10.0,32,0.0,PI2);
   Svzphi->Sumw2(); // collect good statitics
   TH2S *Svzdphi  = new TH2S("Svzdpii","Slow: Posible Z vertecies vs. DPhi",
			     200,+4.0,6.0,20,0.0,10*DPhi);
   Svzdphi->Sumw2(); // collect good statitics
   TH1F *Svzfl  = new TH1F("Svzfl","Slow: Posible Z vertecies Flattened",
			    nbinx,+4.0,6.0);
   Svzfl->Sumw2();
   TH1F *Svztr = new TH1F("Svztr","Slow: Z Vertex found by True Tracks",
			  nbinx,+4.0,6.0);
   Svztr->Sumw2(); // collect good statitics

   printf("Svertex: i1max=%d i2max=%d\n",i1max,i2max);

   for(i=0;i<i1max;i++) for(j=0;j<i2max;j++) {
       dphi = fabs(spdo[j]->phir - spdi[i]->phir); if(dphi>PI) dphi -= PI;
       dz   = fabs(spdo[j]->zr   - spdi[i]->zr);
       r1   = spdi[i]->r;
       r2   = spdo[j]->r;
       if(r2-r1!=0.0) zv   = (r2*spdi[i]->zr - r1*spdo[j]->zr)/(r2-r1);
       else{
	   printf("vertex_slow: spdi[%d]->r=%f = spdo[%d]->r=%f\n",i,r1,j,r2);
	   continue;
       } // end if else r1-r2!=0.0
//       if(j<10&&i<10) printf("zv=%e dphi=%e,r1=%e,r2=%e,dz=%e\n",
//			     zv,dphi,r1,r2,dz);
       Svzdphi->Fill(zv,dphi);
       if(dphi>DPhi) continue; // If outside DPhi (momentum) range, skip it.
       if(dz>DZ)     continue; // If outside dz range, skip it.
       Svzphi ->Fill(zv,spdi[i]->phir);
       Svzfl->Fill(zv);
       if(spdi[i]->track == spdo[j]->track) Svztr->Fill(zv);
   } // end for i,j

   TH1D *Svzdef = Svzphi->ProjectionX("Slow: Posible Z vertecies",0,nbinx,"E");

   i0 = Svzfl->GetBinContent(1);
   x0 = Float_t(1);
   i1 = Svzfl->GetBinContent(nbinx);
   x1 = Float_t(nbinx);
   su = x0-x1; if(su==0.0) return -200.0;
   su = (i0-i1)/su;
   x0 = x0-su*i0;
   x1 = su;
   for(i=1;i<=nbinx;i++){
      su = x1*Float_t(i) + x0;
      Svzfl->AddBinContent(i,-su);
   } // end for i

   printf("Slow:           mean=%f RMS=%f w2=%f\n",
	  Svzdef->GetMean(),Svzdef->GetRMS(),Svzdef->GetSumOfWeights());
   printf("Slow: Flattened mean=%f RMS=%f w2=%f\n",
	  Svzfl->GetMean(),Svzfl->GetRMS(),Svzfl->GetSumOfWeights());

   av = 3.0 * Svzfl->Integral(1,nbinx)/Float_t(nbinx);
   printf("Slow: Flattened av=%f Tvxps->Max=%f nbinx=%d\n",
	  av,Svzfl->GetMaximum(),nbinx);

   {// find best dPhi that masimizes the signal/noise ration.
       Float_t sn=0.0,sig=0.0,nois=0.0,ns;
       Int_t   iend = Svzdphi->GetYaxis()->GetNbins();
       Int_t   jend = Svzdphi->GetXaxis()->GetNbins();
       Int_t   ipeak = Svzdphi->GetXaxis()->FindBin(Svzdef->GetMean());
       nbst = 0;
       for(i=0;i<=iend;i++){
	   sig += Svzdphi->GetCellContent(ipeak,i);
	   ns   = 0.0;
	   for(j=1;j<6;j++) ns += Svzdphi->GetCellContent(j,i);
	   for(j=jend-6;j<jend;j++) ns += Svzdphi->GetCellContent(j,i);
	   nois += 0.1*ns;
	   if(nois<=0.0) continue;
	   if((sig-nois)/nois>sn){
	       sn = (sig-nois)/nois;
	       nbst = i;
	   } // end if
//	   printf("Svertex: bin=%d signal/noise=%e\n",i,(sig-nois)/nois);
       } // end for i
   } // end find best dPhi

   if(nbst<=0) nbst = 1; // must have some data.
   BDphi = Svzdphi->GetYaxis()->GetBinUpEdge(nbst);
   TH1D *Svzbst = Svzdphi->ProjectionX("Slow: Best, Z vertecies",0,nbst,"E");

   { // Start Muinuit fitting
       // initilize Minuit package
       gH1   = Svzbst; // histogram to be fit
       gstrt = Svzbst->GetXaxis()->FindBin(
        Svzdef->GetMean()-2.0*(Svzdef->GetRMS()));//histogram start bin value
       gend  = Svzbst->GetXaxis()->FindBin(
        Svzdef->GetMean()+2.0*(Svzdef->GetRMS()));//histogram end   bin value
       TMinuit *gMinuit = new TMinuit(4); // init Minuit
       gMinuit->SetFCN(Chi2Gauss);  // chi^2 function with Gaussian
       Double_t arglist[10];  // Munuit parameter array
       Int_t ierflg = 0; // Muniut error flag
       //arglist[0] = 1.;
       //gMinuit->mnexcm("SET ERR",arglist,1,ierglg);
       // Set starting values and step size for parameters
       Double_t vstart[4],step[4];
       { // find background
	   Float_t ns   = 0.0;
	   Int_t jend = Svzbst->GetXaxis()->GetNbins();
	   for(j=1;j<6;j++) ns += Svzbst->GetBinContent(j);
	   for(j=jend-6;j<jend;j++) ns += Svzbst->GetBinContent(j);
	   vstart[3] = 0.1*ns;
       } // end find backgrount
       vstart[1] = Svzbst->GetMean();
       vstart[2] = Svzbst->GetRMS();
       if(vstart[2] == 0.0) vstart[2] = 0.04;
       vstart[0] = (Svzbst->GetEntries() - (Svzbst->GetNbinsX())*vstart[3])*
	   vstart[2]/TMath::Sqrt(TMath::Pi());
       if(vstart[0]<=0.0) vstart[0] = 1.0;
       for(i=0;i<4;i++) step[i] = 0.05*vstart[i];
       if(vstart[3] <= 0.0) step[3] = 0.1;
       step[1] = 0.01; // mean expected at about zero set step my hand
       gMinuit->mnparm(0,"Const",vstart[0],step[0],0,0,ierflg);
       gMinuit->mnparm(1,"Mean" ,vstart[1],step[1],0,0,ierflg);
       gMinuit->mnparm(2,"Sigma",vstart[2],step[2],0,0,ierflg);
       gMinuit->mnparm(3,"Offst",vstart[3],step[3],0,0,ierflg);
       // Now ready for minimization step
       //arglist[0] = 500.; // Maximum number of calls
       //arglist[1] = 1.;   // Tolorance
       gMinuit->mnexcm("MIGRAD",arglist,0,ierflg); // do minimization
       gMinuit->mnexcm("MINO",arglist,0,ierflg);   // find best errors
       // Get results
       Double_t parmin,edm,errdef;
       Int_t nvpar,nparx,icstat;
       gMinuit->mnstat(parmin,edm,errdef,nvpar,nparx,icstat);
       printf("Svertex: chi2gauss=%e edist=%e istat=%d dPhi=%e\n",
	      parmin,edm,icstat,BDphi);
       TString chnam = TString(10);
       Double_t par[4],epar[4],empar[4],eppar[4];
       for(i=0;i<4;i++){ 
	   gMinuit->mnpout(i,chnam,par[i],epar[i],empar[i],eppar[i],ierflg);
	   gMinuit->mnerrs(i,eppar[i],empar[i],epar[i],arglist[i]);
	   printf("Svertex: %s = %e +- %e (%e,%e)\n",
		  chnam.Data(),par[i],epar[i],empar[i],eppar[i]);
       } // end for i
   } // End Muinuit fitting
   su  = 0.0;
   avt = 0.0;
   i1  = Svzfl->GetBinContent(1);
   i2  = Svzfl->GetBinContent(2);
   for(i=2;i<nbinx;i++){
       i0 = i1; // old i-1 value
       i1 = i2; // old i value
       i2 = Svzfl->GetBinContent(i+1); // new i+1 value
       if(i1 > av && (i0 > av || i2 > av ) ){
	   if(!mod) mod = kTRUE;
       } else if(mod) break;
       if(mod){ // inside peak
	   su  += i1;
	   avt += i1 * Svzfl->GetXaxis()->GetBinCenter(i);
       } // end if mod
   } // end for i

   if(su!=0.0) zv = avt/su;
   else        zv = -100.0; // an unphysical value

   TCanvas *c1 = new TCanvas("c1","Slow Alice ITS vertex finder",
        	             450,10,600,700);
   Svzphi->Draw("col");
   c1->Print("vertex5_S_vz_phi.eps");
   Svzdphi->Draw();
   c1->Print("vertex5_S_vz_dphi.eps");
   Svzdef->Draw();
   c1->Print("vertex5_S_vz_def.eps");
   Svzbst->Draw();
   c1->Print("vertex5_S_vz_bst.eps");
   Svztr->Draw();
   c1->Print("vertex5_S_vz_true.eps");

   return zv;
}

Float_t vertexEta(Hit_tl **spdi,Int_t i1max,Hit_tl **spdo,Int_t i2max){
   Float_t r1    = 3.910078;// radius at which hit is from beam axis for spd1
   Float_t r2    = 6.955933;// radius at which hit is from beam axis for spd2
   Float_t DPhi  = 0.005;   // maximum allowed difference in phi angle
   Float_t DEta  = 0.100;   // maximum allowed difference in eta/pseudorapidity
   Float_t DZ    = 12.5;     // maximum allowed difference in z position
   Float_t avt   = 0.0;
   Int_t   nbinx = 2000;
   Int_t   start = 0;
   Bool_t  mod   = kFALSE;
   Int_t   i,j;
   Float_t zv=0.0,av,su,i0,i1,i2,x0,x1,dphi,dz;
   Float_t PI2 = 2.0*TMath::Pi();
   Float_t PI  = TMath::Pi();

   // sort according to phi.
   hpsortEta(spdi,i1max);
   hpsortEta(spdo,i2max);

   // find best vertex allong z.
   TH2S *Evzphi   = new TH2S("Evzphi","Eta: Posible Z vertecies vs. Phi",
			     nbinx,-5.0,5.0,32,0.0,PI2);
   Evzphi->Sumw2(); // collect good statitics
   TH2S *Evzdphi  = new TH2S("Evzdphi","Eta: Posible Z vertecies vs. DPhi",
			     200,-1.0,1.0,20,0.0,10*DPhi);
   Evzdphi->Sumw2(); // collect good statitics
   TH1F *Evzfl   = new TH1F("Evzfl","Eta: Posible Z vertecies Flattened",
			  nbinx,-5.0,5.0);
   Evzfl->Sumw2();
   TH1F *Evztr  = new TH1F("Evztr","Eta: Z Vertex found by True Tracks",
			   nbinx,-5.0,5.0);
   Evztr->Sumw2(); // collect good statitics

   for(i=0;i<i1max;i++){
       for(;spdo[start]->etar < spdi[i]->etar-1.5*DEta;start++);
       for(j=start;(spdo[j]->etar < spdi[i]->etar+DEta) && (j < i2max);j++){
	   dphi = fabs(spdi[i]->phir - spdo[j]->phir); if(dphi>PI) dphi -= PI;
	   dz   = fabs(spdi[i]->zr   - spdo[j]->zr);
	   r1   = spdi[i]->r;
	   r2   = spdo[j]->r;
	   if(r2-r1!=0.0) zv   = (r2*spdi[i]->zr - r1*spdo[j]->zr)/(r2-r1);
	   else printf("vertex_Eta: spdi[%d]->r=%f = spdo[%d]->r=%f\n",
		       i,r1,j,r2);
	   Evzdphi->Fill(zv,dphi);
	   if(dphi>DPhi) continue;
	   if(dz>DZ)     continue;
	   Evzphi->Fill(zv,spdi[i]->phir);
	   Evzfl->Fill(zv);
	   if(spdi[i]->track == spdo[j]->track) Evztr ->Fill(zv);
       } // end for j
   } // end for i

   TH1D *Evzdef  = Evzphi->ProjectionX("Eta: Z vertecies Eta",1,nbinx,"E");

   i0 = Evzfl->GetBinContent(1);
   x0 = Float_t(1);
   i1 = Evzfl->GetBinContent(nbinx);
   x1 = Float_t(nbinx);
   su = x0-x1; if(su==0.0) return -200.0;
   su = (i0-i1)/su;
   x0 = x0-su*i0;
   x1 = su;
   for(i=1;i<=nbinx;i++){
      su = x1*Float_t(i) + x0;
      Evzfl->AddBinContent(i,-su);
   } // end for i

   printf("Eta:           mean=%f RMS=%f w2=%f\n",
	  Evzdef->GetMean(),Evzdef->GetRMS(),Evzdef->GetSumOfWeights());
   printf("Eta: Flattened mean=%f RMS=%f w2=%f\n",
	  Evzfl->GetMean(),Evzfl->GetRMS(),Evzfl->GetSumOfWeights());

   av = 3.0 * Evzfl->Integral(1,nbinx)/Float_t(nbinx);
   printf("Eta: Flattened av=%f TvxpE->Max=%f nbinx=%d\n",
	  av,Evzfl->GetMaximum(),nbinx);

   su  = 0.0;
   avt = 0.0;
   i1  = Evzfl->GetBinContent(1);
   i2  = Evzfl->GetBinContent(2);
   for(i=2;i<nbinx;i++){
       i0 = i1; // old i-1 value
       i1 = i2; // old i value
       i2 = Evzfl->GetBinContent(i+1); // new i+1 value
       if(i1 > av && (i0 > av || i2 > av ) ){
	   if(!mod) mod = kTRUE;
       } else if(mod) break;
       if(mod){ // inside peak
	   su  += i1;
	   avt += i1 * Evzfl->GetXaxis()->GetBinCenter(i);
       } // end if mod
   } // end for i

   if(su!=0.0) zv = avt/su;
   else        zv = -100.0; // an unphysical value

   TCanvas *c2 = new TCanvas("c2","Alice ITS vertex finder Eta", 
			     500,10,600,700);
   Evzphi->Draw("col");
   c2->Print("vertex5_E_vz_phi.ps");
   Evzdphi->Draw();
   c2->Print("vertex5_E_vz_dphi.ps");
   Evzdef->Draw();
   c2->Print("vertex5_E_vz_def.ps");
   Evzfl->Draw();
   c2->Print("vertex5_E_vz_flat.ps");
   Evztr->Draw();
   c2->Print("vertex5_E_vz_true.ps");

   return zv;
}
