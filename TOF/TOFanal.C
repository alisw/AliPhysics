void tofanal (Int_t evNumber=0) 
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
<img src="picts/tofanal.gif">
*/
//End_Html
/////////////////////////////////////////////////////////////////////////


// Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
     gROOT->LoadMacro("loadlibs.C");
     loadlibs();
   } else {
      delete gAlice;
      gAlice = 0;
   }
      
// Connect the Root Galice file containing Geometry, Kine and Hits
   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   if (!file) file = new TFile("galice.root");

// Get AliRun object from file or create it if not on file
   if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","TOF test program");
   }
      
// Import the Kine and Hits Trees for the event evNumber in the file
   gAlice->GetEvent(evNumber);
   Float_t x,y,z,mass,e;
   Int_t nbytes = 0;
   Int_t j,hit,ipart;
   Int_t nhits;
   Float_t tof;
   TParticle *particle;

// Get pointers to Alice detectors and Hits containers
   AliDetector *TOF  = gAlice->GetDetector("TOF");
   TClonesArray *Particles = gAlice->Particles();

   Int_t ntracks    = gAlice->TreeH()->GetEntries();

   // Create histograms
   TH1F *hTOF = new TH1F("TOF","Time-of-flight distribution",100,0,10e-8);
   TH1F *hTOFprim = new TH1F("TOFprim","Time-of-flight distribution of primaries",100,0,10e-8);
// Start loop on tracks in the hits containers
   for (Int_t track=0; track<ntracks;track++) {
     if(TOF) {
     // ======>Histogram TOF
       for(AliTOFhit* tofHit=(AliTOFhit*)TOF->FirstHit(track); tofHit; tofHit=(AliTOFhit*)TOF->NextHit()) {
         tof = tofHit->fTof;
         hTOF->Fill(tof);
         ipart    = tofHit->fTrack;
         particle = (TParticle*)Particles->UncheckedAt(ipart);
         if (particle->GetFirstMother() < 0) hTOFprim->Fill(tof);
       }
     }        
   }

//Create a canvas, set the view range, show histograms
   TCanvas *c1 = new TCanvas("c1","Alice TOF hits",400,10,600,700);
   c1->Divide(1,2);
   c1->cd(1);
   gPad->SetFillColor(33);
   gPad->SetLogy();
   hTOF->SetFillColor(42);
   hTOF->Draw();
   //   hSectors->Fit("pol1");
   c1->cd(2);
   gPad->SetFillColor(33);
   gPad->SetLogy();
   hTOFprim->SetFillColor(42);
   hTOFprim->Draw();
   //   hTOFprim->Draw("same");
   c1->Print("tofanal.ps");
}
