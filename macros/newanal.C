void newanal (Int_t evNumber=0) 
{
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and fill some histograms.
//   
//     Root > .L anal.C   //this loads the macro in memory
//     Root > anal();     //by default process first event   
//     Root > anal(2);    //process third event
//Begin_Html
/*
<img src="gif/anal.gif">
*/
//End_Html
/////////////////////////////////////////////////////////////////////////


// Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
     gROOT->LoadMacro("loadlibs.C");
     loadlibs();
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

// Get pointers to Alice detectors and Hits containers
   AliDetector *TPC  = gAlice->GetDetector("TPC");
   AliDetector *TRD  = gAlice->GetDetector("TRD");
   TClonesArray *Particles = gAlice->Particles();

   Int_t ntracks    = gAlice->TreeH()->GetEntries();

   // Create histograms
   TH1F *hSectors = new TH1F("hSectors","Number of tracks hits per sector",75,1,76);
   TH1F *hTRD     = new TH1F("hTRD","Number of planes crossed per track in TRD",6,1,7);
   TH1F *hTRDprim = new TH1F("hTRDprim","Number of planes crossed per primary track in TRD",6,1,7);
   
// Start loop on tracks in the hits containers
   for (Int_t track=0; track<ntracks;track++) {
     if(TPC) {
     // ======>Histogram TPC
     // histogram number of tracks per sector
       for(AliTPChit* tpcHit=(AliTPChit*)TPC->FirstHit(track); tpcHit; tpcHit=(AliTPChit*)TPC->NextHit()) {
         sector = tpcHit->fPadRow;
         hSectors->Fill(sector);
       }
     }        
     // =======>Histogram TRD
     // histogram number of planes crossed per track
     // same for primary tracks only
     if (TRD) {
       for(AliTRDhit* trdHit=(AliTRDhit*)TRD->FirstHit(-1); trdHit; trdHit=(AliTRDhit*)TRD->NextHit()) {
         ipart    = trdHit->fTrack;
         particle = (GParticle*)Particles->UncheckedAt(ipart);
         plane = trdHit->fPlane;
         hTRD->Fill(plane);
         if (particle->GetParent() < 0) hTRDprim->Fill(plane);
       }
     }        
   }

//Create a canvas, set the view range, show histograms
   TCanvas *c1 = new TCanvas("c1","Alice TPC and ITS hits",400,10,600,700);
   c1->Divide(1,2);
   c1->cd(1);
   gPad->SetFillColor(33);
   hSectors->SetFillColor(42);
   hSectors->Fit("pol1");
   c1->cd(2);
   gPad->SetFillColor(33);
   hTRD->SetFillColor(42);
   hTRD->Draw();
   hTRDprim->SetFillColor(46);
   hTRDprim->Draw("same");
   c1->Print("anal.ps");
}
