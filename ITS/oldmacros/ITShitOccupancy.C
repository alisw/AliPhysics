/*
// Some time it is nice to compile the code to find all of the typos.
#include <stdio.h>
#include <math.h>
#include "TROOT.h"
#include "TSystem.h"
#include "TClassTable.h"
#include "TFile.h"
#include "AliITS.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TArray.h"
#include "TCanvas.h"
#include "AliRun.h"
#include "AliITSgeom.h"
#include "TParticlet.h"

extern TSystem     *gSystem;
extern TROOT       *gROOT;
extern TClassTable *gClassTable;
extern AliRun      *gAlice;
*/

Double_t dEtadZ(Double_t z,Double_t r){
    Float_t a;

    a = TMath::Sqrt(z*z+r*r);
    if(a!=0.0) a = 1./a;
    return(a);
}

Double_t PdEtadZ(Double_t *z,Double_t *r){
    return(dEtadZ(*z,*r));
}

Double_t ZfromEtaR(Double_t eta,Double_t r){
    Double_t z;

    z   = 2.0*TMath::ATan(TMath::Exp(-eta));
    z   = r/TMath::Tan(z);
    return(z);
}

Double_t EtafromZr(Double_t z,Double_t r){
    Double_t a,eta;

    a   = 0.5*TMath::ATan2(r,z);    // Tan2(y,x) = y/x
    eta = -TMath::Log(TMath::Abs(TMath::Tan(a)));
    if(a<0.0) eta = -eta;
    return (eta);
}

Double_t PEtafromZr(Double_t *z,Double_t *r){
    return (EtafromZr(*z,*r));
}

Double_t Weight(Double_t *a,Double_t *b){
    Double_t eta,z,Etap,rt,dEta,dEtap,weight,REtap;

    eta    = a[0];
    Etap   = a[1];
    rt     = b[0];
    dEta   = b[1];
    dEtap  = b[2];
    REtap  = b[3];
    z      = ZfromEtaR(eta,rt);
    weight = TMath::Abs(dEtadZ(z,rt));
    weight = weight*REtap/(rt*2.0*TMath::Pi()*dEta*dEtap);
    return weight;
}

void ITShitOccupancy (const char *filename="galice.root",Int_t evNumber=0){
/////////////////////////////////////////////////////////////////////////
//
//   Written by Roberto Barbera
//   Modified by Bjorn S. Nilsen May 19 1999
//
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of aliroot
//   and fill some histograms. The Occupancy values produce are only
//   approximate and are known to have problems. No warranty of any kind
//   may be give or assumed. 
//
//     Root > .L ITShitOccupancy.C//this loads the macro in memory
//     Root > ITShitOccupancy();  //by default process first event   
//     Root > ITShitOccupancy("galiceA.root");  //process file galiceA.root
//     Root > ITShitOccupancy("galiceA.root",2);//process file galiceA.root 
//                                                 third event
//
//     Root > .x ITShitOccupancy.C(); //by default process first event   
//     Root > .x ITShitOccupancy.C("galiceA.root"); //process file galiceA.root
//     Root > .x ITShitOccupancy.C("galiceA.root",2);//process file 
//                                                     galiceA.root third event
//     The test figures for this macro are kept in 
// $ALICE_ROOT/ITS/figures/ITShitOccupancy.figures as a series of post-
// script files and one ITShitOccupancy.root file. This is because ther
// are so many figures produced.
//
////////////////////////////////////////////////////////////////////////
//
// The code computes the Occupancy of HITS in the different ITS detectors
// as a function of z. This does not represent the true occupancy expected
// in the ITS since there will be more than one hit per charge cluster in
// the different detector of the ITS. To do a proper occupancy calculation
// a full simulation of all of the detectors in the ITS will be needed. In
// addition a proper background needs to be used.
//
//  This code has been modified to be as compatible with the standard (May 26
//  1999) version of the official released ITS code as possible. Functions 
//  introduce in the unofficial release of the ITS code (March 199) have
//  been commented out and replaced by the variables defined in the ITShits
//  class. As soon as these functions are officially part of the aliroot 
//  release. These functions should be used instead to maximise future 
//  portability.
//
//  As with all ITS software produced under the deadline of the ITS TDR, there
//  is no real guarantee that this code will function, produce meaning full
//  results, solve in any way any possible problem or question, or be supported
//  in any kind. Questions or comments can be addressed to 
//  Nilsen@mps.ohio-state.edu (for as long as I am there) or 
//  Roberto.Barbera@ct.infn.it. A response may never come.
//
//  This exclamation should be copy written and thereby not used my any other
//  software writer, especially those writing code for ALICE. We expect 
//  everyone else's code to work perfectly and do the job as advertised, either
//  explicitly, by word of mouth, or by any other mystical means.
//
////////////////////////////////////////////////////////////////////////
// Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   } // end if gClassTable...
//
// Connect the Root Galice file containing Geometry, Kine and Hits
   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
   if (!file) file = new TFile(filename);
   printf("Reading from root file %s\n",filename);
//
// Get AliRun object from file or create it if not on file
   if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   } // end if !gAlice
//
// Get pointers to Alice detectors and Hits containers
   Int_t  nparticles  = gAlice->GetEvent(evNumber);
   if (nparticles <= 0) return;
   AliITS       *ITS       = gAlice->GetDetector("ITS");
   if(!ITS) {printf("ITS detector not found. Exiting\n"); return;}
   TClonesArray *Particles = gAlice->Particles();
   TClonesArray *ITShits   = ITS->Hits();
   TTree        *TH        = gAlice->TreeH();
   Int_t        ntracks    = TH->GetEntries();
   AliITSgeom   *gm        = ITS->GetITSgeom();
   Int_t        Nlayers    = gm->GetNlayers();
//
   printf("%d primary tracks and their secondary tracks read in\n",ntracks);
//
// Import the Kine and Hits Trees for the event evNumber in the file
   Float_t x,y,z,mass,e,r,rt;
   Float_t px,py,pz,pt;
   Float_t destep;
   Float_t theta,phi;
   Float_t eta,etap;
   Float_t dzeta,weight;
   Int_t idpart;
//
   Float_t ztot[]  = {33.52,// total length of layer 1
                      33.52,// total length of layer 2
                      46.10,// total length of layer 3
                      60.77,// total length of layer 4
                      90.42,// total length of layer 5
                     101.95};// total length of layer 6
   Float_t raprx[] = {4.0,7.0,14.9,23.8,39.1,43.6}; // aproximate layer radii
//
   Float_t totdet[]={256*256*4*20,//total number of channels layer 1 =  5242880
                     256*256*4*40,//total number of channels layer 2 = 10485760
                   2*256*256*6*14,//total number of channels layer 3 = 11010048
                   2*256*256*8*22,//total number of channels layer 4 = 23068672
                        768*23*34,//total number of channels layer 5 =  1201152
                        768*26*38};//total number of channels layer 6 = 1517568
//
   Float_t clusize,clusize_z,clusize_rphi;   
//
   Float_t threshold[]={8.0e-6,8.0e-6, // noise= 200 e- ; thr = 2000 e- ;
                        3.0e-6,3.0e-6, // noise= 250 e- ; thr= 3 sigma noise
                        6.0e-6,6.0e-6};// noise= 500 e- ; thr= 3 sigma noise
//
   Int_t     nhitlayer[6];
   Int_t     ntrklayer[6];
   Int_t     layer,ladder,det;
   Int_t     nbytes = 0;
   Int_t     i,j,jj,hit,ipart,ipartold,iparent;
   Int_t     nhits;
   Int_t     sector,plane;
   TParticle *particle,*pparticle;
   AliITShit *itsHit;

// Create histograms

   Float_t etamin = -1.;
   Float_t etamax = 1.;
   Float_t ptmin  = 0.;
   Float_t ptmax  = 3.;
   Int_t   nbin   = 40;
   Int_t   nbinz  = 200;
   Float_t deta   = (etamax-etamin)/nbin;
   Float_t EtapMax=2.436246,EtapMin=-2.436246; // from Config.C file
   TVector etabin(nbin+1);
   TVector zbin(nbin+1);

// Open a root file to save histograms in
   char *Hfilename = "ITShitOccupancy.root";
   TFile *Hfile    = (TFile *)gROOT->GetListOfFiles()->FindObject(Hfilename);
   if(Hfile) Hfile->Close();
   Hfile = new TFile(Hfilename,"RECREATE","Histograms for ITS Occupancy");
   printf("Writting histograms to root file %s\n",Hfilename);
//
   TH2F *hzetap[6];
   TH1D *hzetap_z[6],*hzetap_etap[6];
   TH2F *hzetapden[6];
   TH1D *hzetapden_z[6],*hzetapden_etap[6];
   TH2F *hEtaEtapden[6];
   TH1D *hEtaEtapden_Eta[6],*hEtaEtapden_Etap[6];
   TH2F *hEtaEtapdenWeight[6];
   TH1D *hEtaEtapdenWeight_Eta[6];
   TF2  *fweight;
   TF1  *fdEtadZ,*fEta;
   TH1F *hdEtadZ[6],*hEta[6];

   { // "Book" histograms
       char histid[7],htit[40];
       for(i=0;i<Nlayers;i++){
	   sprintf(histid,"hzetap%1.1d",i+1);
	   sprintf(htit,"Hits for Layer %1.1d",i+1);
	   hzetap[i] = new TH2F(histid,htit,nbinz,-0.6*ztot[i],0.6*ztot[i],
				20,floor(EtapMin)-1.0,ceil(EtapMax)+1.0);
	   hzetap[i]->SetXTitle("Z (cm)");
	   hzetap[i]->SetYTitle("Beam particle pseudorapidity");
	   hzetap[i]->Sumw2();
	   //
	   sprintf(histid,"hzetapden%1.1d",i+1);
	   sprintf(htit,"Density of Hits for Layer %1.1d",i+1);
	   hzetapden[i] = new TH2F(histid,htit,nbinz,-0.6*ztot[i],0.6*ztot[i],
			       20,floor(EtapMin)-1.0,ceil(EtapMax)+1.0);
	   hzetapden[i]->Sumw2();
	   hzetapden[i]->SetXTitle("Z (cm)");
	   hzetapden[i]->SetYTitle("Beam particle pseudorapidity");
	   nhitlayer[i] = 0.0;
	   ntrklayer[i] = 0.0;
	   //
	   sprintf(histid,"hEtaEtapden%1.1d",i+1);
	   sprintf(htit,"Hits for Layer %1.1d by Eta",i+1);
	   hEtaEtapden[i] = new TH2F(histid,htit,nbinz,
				  -EtafromZr(0.6*ztot[i],raprx[i]),
				   EtafromZr(0.6*ztot[i],raprx[i]),
				20,floor(EtapMin)-1.0,ceil(EtapMax)+1.0);
	   hEtaEtapden[i]->SetXTitle("Track Pseudorapidity");
	   hEtaEtapden[i]->SetYTitle("Beam particle pseudorapidity");
	   hEtaEtapden[i]->Sumw2();
       } // end for i

       sprintf(histid,"fweight");
       fweight = new TF2(histid,Weight,-EtafromZr(0.6*ztot[0],raprx[0]),
				   EtafromZr(0.6*ztot[0],raprx[0]),
			 floor(EtapMin)-1.0,ceil(EtapMax)+1.0,4);
       Double_t params[4];
       sprintf(histid,"fdEtadZ");
       fdEtadZ = new TF1(histid,PdEtadZ,floor(EtapMin)-1.0,
			 ceil(EtapMax)+1.0,1);
       sprintf(histid,"fEta");
       fEta = new TF1(histid,PEtafromZr,floor(EtapMin)-1.0,
		      ceil(EtapMax)+1.0,1);
       for(i=0;i<Nlayers;i++){
	   params[0] = raprx[i];
	   params[1] = hEtaEtapden[i]->GetXaxis()->GetBinWidth(1);
	   params[2] = hEtaEtapden[i]->GetYaxis()->GetBinWidth(1);
	   params[3] = EtapMax-EtapMin;
	   fweight->SetParameters(params);
	   sprintf(histid,"hEtaEtapdenWeight%1.1d",i+1);
	   sprintf(htit,"Weight for Layer %1.1d by Eta",i+1);
	   hEtaEtapdenWeight[i] = new TH2F(histid,htit,nbinz,
				  -EtafromZr(0.6*ztot[i],raprx[i]),
				   EtafromZr(0.6*ztot[i],raprx[i]),
				20,floor(EtapMin)-1.0,ceil(EtapMax)+1.0);
	   hEtaEtapdenWeight[i]->SetXTitle("Track Pseudorapidity");
	   hEtaEtapdenWeight[i]->SetYTitle("Beam particle pseudorapidity");
	   hEtaEtapdenWeight[i]->Sumw2();
	   hEtaEtapdenWeight[i]->Eval(fweight);
//
	   fdEtadZ->SetParameters(params);
	   sprintf(histid,"hdEtadZ%1.1d",i+1);
	   sprintf(htit,"dEtadZ(z) for Layer %1.1d",i+1);
	   hdEtadZ[i] = new TH1F(histid,htit,nbinz,
				  -0.6*ztot[i],0.6*ztot[i]);
	   hdEtadZ[i]->SetXTitle("Z (cm)");
	   hdEtadZ[i]->Sumw2();
	   hdEtadZ[i]->Eval(fdEtadZ);
//
	   fEta->SetParameters(params);
	   sprintf(histid,"hEta%1.1d",i+1);
	   sprintf(htit,"Eta(z) for Layer %1.1d",i+1);
	   hEta[i] = new TH1F(histid,htit,nbinz,
				  -0.6*ztot[i],0.6*ztot[i]);
	   hEta[i]->SetXTitle("Z (cm)");
	   hEta[i]->Sumw2();
	   hEta[i]->Eval(fEta);
       } // end for i

   } //  end "Book" histograms
// Start loop on tracks in the hits container

   Double_t xa[2],par[4];

   for (Int_t track=0; track<ntracks;track++) {
     gAlice->ResetHits();
     nbytes += TH->GetEvent(track);
     nhits   = ITShits->GetEntriesFast();
     for (hit=0;hit<nhits;hit++) {
	 itsHit   = (AliITShit*)ITShits->UncheckedAt(hit);
	 destep   = itsHit->GetIonization();
	 // With this new version, to be able to do proper detector
	 // simulations, the requirment that a "hit" leave some
	 // energy deposited has been removed.
	 if(destep<=0.0) continue; // skip hits without energy loss.
	 ipart    = itsHit->fTrack;
	 particle = (TParticle*)Particles->UncheckedAt(ipart);
	 iparent  = particle->GetFirstMother();
	 pparticle= (TParticle*)Particles->UncheckedAt(ntracks-track);
	 etap     = TMath::ATan2(pparticle->Pt(),pparticle->Pz());
	 etap     = -TMath::Log(TMath::Tan(0.5*etap));
	 itsHit->GetPositionG(x,y,z);
	 x -= pparticle->Vx();  // correct for primary vertex position
	 y -= pparticle->Vy();  // correct for primary vertex position
	 z -= pparticle->Vz();  // correct for primary vertex position
	 itsHit->GetMomentumG(px,py,pz);
	 itsHit->GetDetectorID(layer,ladder,det);
	 r        = sqrt(x*x+y*y+z*z);
	 rt       = sqrt(x*x+y*y);
	 theta    = TMath::ACos(z/r)*180./TMath::Pi();
	 phi      = TMath::ATan2(y,x); if(phi<0.0) phi += TMath::Pi();
	 eta      = EtafromZr(z,rt);
	 pt       = TMath::Sqrt(px*px+py*py);
	 if(ipart!=ipartold)ntrklayer[layer-1] += 1;
	 ipartold = ipart;
	 nhitlayer[layer-1] += 1;

	 switch (layer) {
        // Calculate average cluster sizes
	 case 1:
	     clusize_rphi = 1.4;
	     if(TMath::Abs(eta)>=0.0 && TMath::Abs(eta)< 0.2) clusize_z=1.2;
	     if(TMath::Abs(eta)>=0.2 && TMath::Abs(eta)< 0.4) clusize_z=1.2;
	     if(TMath::Abs(eta)>=0.4 && TMath::Abs(eta)< 0.6) clusize_z=1.3;
	     if(TMath::Abs(eta)>=0.6 && TMath::Abs(eta)< 0.8) clusize_z=1.4;
	     if(TMath::Abs(eta)>=0.8 && TMath::Abs(eta)< 1.0) clusize_z=1.5;
	     clusize = clusize_z*clusize_rphi;
	 case 2:
	     clusize_rphi = 1.4;
	     if(TMath::Abs(eta)>=0.0 && TMath::Abs(eta)< 0.2) clusize_z=1.2;
	     if(TMath::Abs(eta)>=0.2 && TMath::Abs(eta)< 0.4) clusize_z=1.2;
	     if(TMath::Abs(eta)>=0.4 && TMath::Abs(eta)< 0.6) clusize_z=1.3;
	     if(TMath::Abs(eta)>=0.6 && TMath::Abs(eta)< 0.8) clusize_z=1.4;
	     if(TMath::Abs(eta)>=0.8 && TMath::Abs(eta)< 1.0) clusize_z=1.5;
	     clusize = clusize_z*clusize_rphi;
	 case 3:
	     clusize_rphi = 8.0;
	     clusize_z    = 3.0;
	     clusize      = clusize_z*clusize_rphi;
	 case 4:
	     clusize_rphi = 8.0;
	     clusize_z    = 3.0;
	     clusize      = clusize_z*clusize_rphi;
	 case 5:
	     clusize_z = 1.0;
	     if(pt>=0.0 && pt< 0.2) clusize_rphi=2.0;
	     if(pt>=0.2 && pt< 0.4) clusize_rphi=1.9;
	     if(pt>=0.4 && pt< 0.6) clusize_rphi=1.6;
	     if(pt>=0.6 && pt< 0.8) clusize_rphi=1.6;
	     if(pt>=0.8 && pt< 1.0) clusize_rphi=1.7;
	     if(pt>=1.0 && pt< 1.5) clusize_rphi=1.7;
	     if(pt>=1.5 && pt< 2.5) clusize_rphi=1.7;
	     if(pt>=2.5           ) clusize_rphi=1.7;
	     clusize = clusize_z*clusize_rphi;
	 case 6:
	     clusize_z = 1.0;
	     if(pt>=0.0 && pt< 0.2) clusize_rphi=2.0;
	     if(pt>=0.2 && pt< 0.4) clusize_rphi=1.9;
	     if(pt>=0.4 && pt< 0.6) clusize_rphi=1.6;
	     if(pt>=0.6 && pt< 0.8) clusize_rphi=1.6;
	     if(pt>=0.8 && pt< 1.0) clusize_rphi=1.7;
	     if(pt>=1.0 && pt< 1.5) clusize_rphi=1.7;
	     if(pt>=1.5 && pt< 2.5) clusize_rphi=1.7;
	     if(pt>=2.5           ) clusize_rphi=1.7;
	     clusize = clusize_z*clusize_rphi;
	 } // end switch layer
//       weightE = clusize*ztot[layer-1]/totdet[layer-1];
	 xa[0]  = eta; xa[1]  = etap;
	 par[0] = rt;
	 par[1] = hEtaEtapden[layer-1]->GetXaxis()->GetBinWidth(1);
	 par[2] = hEtaEtapden[layer-1]->GetYaxis()->GetBinWidth(1);
	 par[3] = EtapMax-EtapMin;
	 weight = Weight(xa,par);
	 hEtaEtapden[layer-1]->Fill(eta,etap,weight);
	 weight = 1.;
	 hzetap[layer-1]->Fill(z,etap,weight);
	 weight = par[3]/(rt*2.0*TMath::Pi()*
		       hzetapden[layer-1]->GetXaxis()->GetBinWidth(1)*
	               hzetapden[layer-1]->GetYaxis()->GetBinWidth(1));
	 hzetapden[layer-1]->Fill(z,etap,weight);
       } // end for hit
   } // end for track

   for(i=0;i<Nlayers;i++)printf("No. of hits   in layer %d = %d\n",i+1,nhitlayer[i]);
   for(i=0;i<Nlayers;i++)printf("No. of tracks in layer %d = %d\n",i+1,ntrklayer[i]);

//Create canvases, set the view ranges and show histograms

     TCanvas *canvas = new TCanvas("canvas","ITShitOccupancy",1);

     Bool_t printit = kTRUE;
     char psfilename[40];
     for(i=0;i<Nlayers;i++){
         hzetap[i]->Draw("colz");
	 sprintf(psfilename,"ITShitOccupancy_%1.1d_z_etap.ps",i+1);
	 if(printit) canvas->Print(psfilename);
	 sprintf(psfilename,"hzetap_z%1.1d",i+1);
	 hzetap_z[i]    = hzetap[i]->ProjectionX(psfilename,0,
						 hzetap[i]->GetNbinsY()+1,"E");
	 hzetap_z[i]->SetXTitle("Z (cm)");
	 hzetap_z[i]->SetYTitle("Number of Hits per cm");
	 sprintf(psfilename,"hzetap_etap%1.1d",i+1);
	 hzetap_etap[i] = hzetap[i]->ProjectionY(psfilename,0,
						 hzetap[i]->GetNbinsX()+1,"E");
	 hzetap_etap[i]->SetXTitle("Beam particle pseudorapidity");
	 hzetap_etap[i]->SetYTitle("Number of Hits");
	 hzetap_z[i]->Draw();
	 sprintf(psfilename,"ITShitOccupancy_%1.1d_z.ps",i+1);
	 if(printit) canvas->Print(psfilename);
	 hzetap_etap[i]->Draw();
	 sprintf(psfilename,"ITShitOccupancy_%1.1d_etap.ps",i+1);
	 if(printit) canvas->Print(psfilename);
//
         hzetapden[i]->Draw("colz");
	 sprintf(psfilename,"ITSHitDen_%1.1d_z_etap.ps",i+1);
	 if(printit) canvas->Print(psfilename);
	 sprintf(psfilename,"hzetapden_z%1.1d",i+1);
	 hzetapden_z[i]    = hzetapden[i]->ProjectionX(psfilename,0,
                                        hzetapden[i]->GetNbinsY()+1,"E");
	 hzetapden_z[i]->SetXTitle("Z (cm)");
	 hzetapden_z[i]->SetYTitle("Number of Hits (cm^-2)");
	 sprintf(psfilename,"hzetapden_etap%1.1d",i+1);
	 hzetapden_etap[i] = hzetapden[i]->ProjectionY(psfilename,0,
                                      hzetapden[i]->GetNbinsX()+1,"E");
	 hzetapden_etap[i]->SetXTitle("Beam particle pseudorapidity");
	 hzetapden_etap[i]->SetYTitle("Number of Hits (cm^-2)");
	 hzetapden_z[i]->Draw();
	 sprintf(psfilename,"ITSHitDen_%1.1d_z.ps",i+1);
	 if(printit) canvas->Print(psfilename);
	 hzetapden_etap[i]->Draw();
	 sprintf(psfilename,"ITSHitDen_%1.1d_etap.ps",i+1);
	 if(printit) canvas->Print(psfilename);
//
         hEtaEtapden[i]->Draw("colz");
	 sprintf(psfilename,"ITShitDensity_%1.1d_Eta_Etap.ps",i+1);
	 if(printit) canvas->Print(psfilename);
	 sprintf(psfilename,"hEtaEtapden_Eta%1.1d",i+1);
	 hEtaEtapden_Eta[i]    = hEtaEtapden[i]->ProjectionX(psfilename,0,
                                         hEtaEtapden[i]->GetNbinsY()+1,"E");
	 hEtaEtapden_Eta[i]->SetXTitle("Pseudorapidity");
	 hEtaEtapden_Eta[i]->SetYTitle("Number of Hits (cm^-2)");
	 sprintf(psfilename,"hEtaEtapden_Etap%1.1d",i+1);
	 hEtaEtapden_Etap[i] = hEtaEtapden[i]->ProjectionY(psfilename,0,
                                         hEtaEtapden[i]->GetNbinsX()+1,"E");
	 hEtaEtapden_Etap[i]->SetXTitle("Beam particle pseudorapidity");
	 hEtaEtapden_Etap[i]->SetYTitle("Number of Hits (cm^-2)");
	 hEtaEtapden_Eta[i]->Draw();
	 sprintf(psfilename,"ITShitDensity_%1.1d_Eta.ps",i+1);
	 if(printit) canvas->Print(psfilename);
	 hEtaEtapden_Etap[i]->Draw();
	 sprintf(psfilename,"ITShitDensity_%1.1d_Etap.ps",i+1);
	 if(printit) canvas->Print(psfilename);
//
         hEtaEtapdenWeight[i]->Draw("colz");
	 sprintf(psfilename,"ITShitDensityWeight_%1.1d_Eta_Etap.ps",i+1);
	 if(printit) canvas->Print(psfilename);
	 sprintf(psfilename,"hEtaEtapdenWeight_Eta%1.1d",i+1);
	 hEtaEtapdenWeight_Eta[i] =
                               hEtaEtapdenWeight[i]->ProjectionX(psfilename,
                             0,hEtaEtapdenWeight[i]->GetNbinsY()+1,"E");
	 hEtaEtapdenWeight_Eta[i]->SetXTitle("Pseudorapidity");
	 hEtaEtapdenWeight_Eta[i]->SetYTitle("Weight (cm^-2)");
	 hEtaEtapdenWeight_Eta[i]->Draw();
	 sprintf(psfilename,"ITShitDensityWeight_%1.1d_Eta.ps",i+1);
	 if(printit) canvas->Print(psfilename);
     } // end for i
//
     Hfile->Write();
//
     return;
}
