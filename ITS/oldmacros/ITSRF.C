#define __COMPILED__
#ifdef __COMPILED__
#include <iostream.h>
#include <AliITS.h>
#include <AliITSRiemannFit.h>
#include <AliRun.h>
#include <TROOT.h>
#include <TClonesArray.h>
#include <TFitter.h>
#include <TMinuit.h>
#include <TF2.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TObjArray.h>
#include <TParticle.h>
#include <TTree.h>
#include <TVector3.h>
#include <TVirtualFitter.h>
#endif

//    Writes the data of helices and the position of the true II-ary, to perform fit 
//    of II-ary and calculation of resolution (in another macro)
//
void ITSRF(const char *filename="galice.root", 
	   const char *outfile="RiemannFit.dat",
	   Int_t evnt=0) {
  
  // Connect the Root Galice file containing Geometry, Kine and Hits
  TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
  if (!file){ file = new TFile(filename);printf("Opening new file\n");}
  printf("Root input file is %s\n",filename);
  
  // Get AliRun object from file or create it if not on file
  if (!gAlice) {
    gAlice = (AliRun*)file->Get("gAlice");
    if (gAlice) printf("AliRun object found on file\n");
    if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
  }
  
  // Get pointers to Alice detectors and Hits containers
  AliITS  *ITS       = (AliITS*)gAlice->GetModule("ITS");
  if(!ITS) {cout<<"not ITS"<<endl; return;}
  Int_t nparticles = gAlice->GetEvent(evnt);   // number of particles generated (including daughters)

  //
  //  Treat ITS information per module
  //
  
  Int_t nmodules;
  ITS->InitModules(-1,nmodules); 
  cout<<" Filling modules..."<<endl;
  ITS->FillModules(evnt,-1,nmodules," "," ");

  //
  //  Get the recontruction tree
  //
  TTree *TR = gAlice->TreeR();
  Int_t nent = (Int_t)TR->GetEntries();
  
  //
  // Initialize Riemann-fit class
  //
  AliITSRiemannFit *rfit = new AliITSRiemannFit();
  rfit->InitPoints(evnt,nent,ITS,TR,nparticles);
  
  const TVector3 Bfield(0.0,0.0,0.2); // ALICE magnetic field
  //
  //   Setting the variables
  //
  
  TParticle *MPart;
  
  Double_t r0[3]={0.0,0.0,0.0}, w, chi;
  Double_t Pt, Px, Py, Pz, charge, Ptot;
  Double_t x0, y0, rho, fd0, fphi, z0, vpar, chisql;
  Double_t CorrLin, ResSum; 
  Int_t status[4]={0,0,0,0}, gfit_status, PdgPart;
  
  FILE *outdat=fopen(outfile,"w");

  cout<<endl<<" Looping on "<< nparticles <<" particles..."<<endl;

  Int_t nAcc=0;
  for(Int_t part = 0; part < nparticles ; part++) {
    if (!(part%100)) cout<<part<<"\r";
    MPart = (TParticle*)gAlice->Particle(part);
    PdgPart = MPart->GetPdgCode();

    //
    //  get track fit
    // 
    Pt = MPart->Pt();
    Px = MPart->Px();
    Py = MPart->Py();
    Pz = MPart->Pz();
    Ptot = MPart->P();
    charge = PdgPart/TMath::Abs(PdgPart);  // +/- 1
    gfit_status = (Int_t)rfit->FitHelix(part, charge, Px, Py, Pz, fd0, fphi, x0, y0, 
					rho, w, z0, vpar, chisql, CorrLin, ResSum);
    // With this offset chi = 0  at the 
    // distance of closest approach in bending plane
    //
    chi=-TMath::Pi()-fphi;  
    if(gfit_status > 0) {
      switch(gfit_status) {
      case 1:  status[0]++;break;
      case 2:  status[1]++;break;
      case 11:  status[2]++;break;
      case 12:  status[3]++; break;
      default: cout<<"CASO NON PREVISTO IN K:   "<<gfit_status<<endl;break;
      }
      continue; 
    }

    //
    //
    //
    Int_t written=0;
    written=fprintf(outdat,"%d ",   part); 
    written=fprintf(outdat,"%f ",   r0[0]); 
    written=fprintf(outdat,"%f ",   r0[1]);
    written=fprintf(outdat,"%f \n", r0[2]);
    written=fprintf(outdat,"%f \n", ResSum); 
    written=fprintf(outdat," %f ",  rho);
    written=fprintf(outdat,"%f ",   w);
    written=fprintf(outdat,"%f ",   fphi); 
    written=fprintf(outdat,"%f \n", vpar); 
    nAcc++;

  } // loop on particles
  fclose(outdat);

  printf(" missed fits because of (resp.):\n no six points, no solution to cubic equation, \n z-direction inversion, tangent limit \n");
  printf(" %d %d %d %d %d \n",status[0],status[1],status[2],status[3]);


  return;
} 


