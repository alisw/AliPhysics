#include "iostream.h"

void ITSreadTest (Int_t evNumber1=0,Int_t evNumber2=0) 
{
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and do some analysis.
//   
/////////////////////////////////////////////////////////////////////////

// Dynamically link some shared libs

   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   }

// Connect the Root Galice file containing Geometry, Kine and Hits

   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   if (!file) file = new TFile("galice.root");
   file->ls();

// Get AliRun object from file or create it if not on file

   if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }
 
//
//   Loop over events 
//
   Int_t Nh=0;
   Int_t Nh1=0;
   for (int nev=0; nev<= evNumber2; nev++) {
     Int_t nparticles = gAlice->GetEvent(nev);
     cout << "nev         " << nev <<endl;
     cout << "nparticles  " << nparticles <<endl;
     if (nev < evNumber1) continue;
     if (nparticles <= 0) return;

     TTree *TH = gAlice->TreeH();
     Int_t ntracks = TH->GetEntries();
     cout<<"ntracks "<<ntracks<<endl;

   Int_t nbytes = 0;

   AliITSRawClusterSDD  *ITSclust;
   AliITSRecPoint  *ITSrecp;

// Get pointers to Alice detectors and Digits containers
   AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
   if(!ITS) return;
   TClonesArray *Particles = gAlice->Particles();

     // fill modules with sorted by module hits
     Int_t nmodules;
     ITS->InitModules(-1,nmodules); 
     ITS->FillModules(nev,-1,nmodules," "," ");
     //get pointer to modules array
     TObjArray *ITSmodules = ITS->GetModules();
     AliITShit *itsHit;

     TTree *TR=gAlice->TreeR();
     Int_t nentr=TR->GetEntries();
     printf("Found %d entries in the TreeR (must be one per module per event!)\n",nentr);
     TClonesArray *ITSrecpoints  = ITS->RecPoints();

     // get the Tree for clusters
     ITS->GetTreeC(nev);
     TTree *TC=ITS->TreeC();
     Int_t nent=TC->GetEntries();
     printf("Found %d entries in the tree (must be one per module per event!)\n",nent);
   
     for (Int_t idettype=0;idettype<3;idettype++) {

       TClonesArray *ITSclusters  = ITS->ClustersAddress(idettype);
       //printf ("ITSclusters %p \n",ITSclusters);

          if (idettype != 1) continue;
	  for (Int_t mod=0; mod<nent; mod++) {
	      AliITSmodule *itsModule = (AliITSmodule*)ITSmodules->At(mod);
	      Int_t nhits = itsModule->GetNhits();
              printf("module nhits %d %d\n",mod,nhits);
	      for (Int_t ihit=0;ihit<nhits;ihit++) {
		itsHit=(AliITShit*)itsModule->GetHit(ihit);
                printf("ihit x y z %d %f %f %f\n",ihit,itsHit->GetXG(),itsHit->GetYG(),itsHit->GetZG());
	      }
              
       
              ITS->ResetClusters();
              TC->GetEvent(mod);
	      Int_t nclust = ITSclusters->GetEntries();
	      if (nclust) printf("Found %d clust for module %d in det type %d \n",nclust,mod,idettype);
	      if (!nclust) continue;

              ITS->ResetRecPoints();
              TR->GetEvent(mod);
	      Int_t nrecp = ITSrecpoints->GetEntries();
	      if (nrecp) printf("Found %d recp for module %d  \n",nrecp,mod);

	      for (Int_t clust=0;clust<nclust;clust++) {
		ITSrecp   = (AliITSRecPoint*)ITSrecpoints->UncheckedAt(clust);
		printf("recp %d %f %f %f\n",clust,ITSrecp->GetX(),ITSrecp->GetZ(),ITSrecp->GetQ());
		ITSclust   = (AliITSRawClusterSDD*)ITSclusters->UncheckedAt(clust);
		printf("clust %d %f %f %f\n",clust,ITSclust->X(),ITSclust->Z(),ITSclust->Q());
	      }
	  }        
     }

     ITS->ClearModules();

   }   // event loop 

     cout<<"END  test for clusters and hits "<<endl;

     file->Close();   
}



