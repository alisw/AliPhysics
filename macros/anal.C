void anal (Int_t evNumber=0) 
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
<img src="picts/anal.gif">
*/
//End_Html
/////////////////////////////////////////////////////////////////////////


  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  }
      
  // Connect the Root Galice file containing Geometry, Kine and Hits
  TFile *file = new TFile("galice.root","read");
  if (!file) {
    cerr<<"anal.C: No galice.root file. Please run the simulation\n";
    return 1;
  }

  // Get AliRun object from file
  if (gAlice) delete gAlice;
  gAlice = (AliRun*)file->Get("gAlice");
  if (gAlice) 
    cerr<<"anal.C: AliRun object found on file\n";
  else {
    cerr<<"anal.C: AliRun object not found on file\n";
    return 2;
  }
      
  // Import the Kine and Hits Trees for the event evNumber in the file
  Int_t nparticles = gAlice->GetEvent(evNumber);
  if (nparticles <= 0) return;

  // Get pointers to Alice detectors and Hits containers
  AliTPC *TPC  = gAlice->GetDetector("TPC");
  AliTRD *TRD  = gAlice->GetDetector("TRD");
  AliTRDgeometry *TRDgeometry   = TRD->GetGeometry();


  TTree *TH = gAlice->TreeH();
  Int_t ntracks    = TH->GetEntries();

  // Create histograms
  TH1F *hSectors = new TH1F("hSectors","Number of tracks hits per sector",72,-0.5,71.5);
  TH1F *hTRD     = new TH1F("hTRD","Number of planes crossed per track in TRD",6,-0.5,5.5);
  TH1F *hTRDprim = new TH1F("hTRDprim","Number of planes crossed per primary track in TRD",6,-0.5,5.5);
   
  // Start loop on tracks in the hits containers
  for (Int_t track=0; track<ntracks;track++) {

    // ======>Histogram TPC
    // histogram the number of tracks per sector
    AliTPChit  *tpcHit;
     if (TPC) {
       for (tpcHit=(AliTPChit*)TPC->FirstHit(track);tpcHit;
            tpcHit=(AliTPChit*)TPC->NextHit()) {
         Int_t sector = tpcHit->fSector;
         hSectors->Fill(sector);
       }
     }        
     // =======>Histogram TRD
     // histogram number of planes crossed per track
     // same for primary tracks only
     AliTRDhit  *trdHit;
     if (TRD) {
       for (trdHit=(AliTRDhit*)TRD->FirstHit(track);trdHit;
            trdHit=(AliTRDhit*)TRD->NextHit()) {
         Int_t ipart    = trdHit->Track();
         TParticle * particle = (TParticle*)gAlice->Particle(ipart);
         UShort_t detector = trdHit->GetDetector();
         Int_t plane = TRDgeometry->GetPlane(detector);
         hTRD->Fill(plane);
         if (particle->GetFirstMother() < 0) hTRDprim->Fill(plane);
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











