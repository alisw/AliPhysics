void RICH (Int_t evNumber=0) 
{
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   for RICH
//   
//     Root > .L RICH.C   //this loads the macro in memory
//     Root > RICH();     //by default process first event   
//     Root > RICH(2);    //process third event
/////////////////////////////////////////////////////////////////////////


// Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
      gSystem->Load("libGeant3Dummy.so");   // a dummy version of Geant3
      gSystem->Load("PHOS/libPHOSdummy.so");        // the standard Alice classes 
      gSystem->Load("libgalice.so");        // the standard Alice classes 
   }
      
// Connect the Root Galice file containing Geometry, Kine and Hits
   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   if (!file) file = new TFile("galice.root");

// Get AliRun object from file or create it if not on file
   if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }
   
// Import the Kine and Hits Trees for the event evNumber in the file
   Int_t nparticles = gAlice->GetEvent(evNumber);
   if (nparticles <= 0) return;
   Float_t x,y,z,mass,e;
   Int_t nbytes = 0;
   Int_t j,hit,ipart;
   Int_t nhits;
   Int_t sector,plane;
   GParticle *particle;
   AliTPChit  *tpcHit;
   AliTRDhit  *trdHit;

// Get pointers to Alice detectors and Hits containers
   AliRICH *RICH  = gAlice->GetDetector("RICH");
   TClonesArray *Particles = gAlice->Particles();
   if (RICH) {
     TClonesArray *RICHhits   = RICH->Hits();
     TClonesArray *Mips       = RICH->Mips();
     TClonesArray *Ckovs      = RICH->Ckovs();
     TClonesArray *Padhits    = RICH->Padhits();
   }
   printf("RICHhits, Mips, Ckovs, Pad %d, %d, %d, %d\n",RICHhits,Mips,Ckovs,Padhits);

   TTree *TH = gAlice->TreeH();
   Int_t ntracks    = TH->GetEntries();

   // Start loop on tracks in the hits containers
   for (Int_t track=0; track<ntracks;track++) {
     gAlice->ResetHits();
     nbytes += TH->GetEvent(track);
     Int_t nrich=RICHhits->GetEntriesFast();
     if(nrich) {
       printf("Number of hit entries %d\n",nrich);
       printf("Number of Mips entries %d\n",Mips->GetEntriesFast());
       printf("Number of Ckovs entries %d\n",Ckovs->GetEntriesFast());
       printf("Number of Padhits entries %d\n",Padhits->GetEntriesFast());
     }
   }
}
