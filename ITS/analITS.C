void analITS (Int_t evNumber=0) 
{
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and fill some histograms.
//   
//     Root > .L analITS.C   //this loads the macro in memory
//     Root > analITS();     //by default process first event   
//     Root > analITS(2);    //process third event
//Begin_Html
/*
<img src="figures/analITS_ref.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>A reference plot produced by analITS.C.
</font>
<pre>
 */
//End_Html
/////////////////////////////////////////////////////////////////////////


// Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   } // end if gClassTable...
      
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
   Float_t x,y,z,mass,e,r;
   Int_t nbytes = 0;
   Int_t j,hit,ipart;
   Int_t nhits;
   Int_t sector,plane;
   TParticle *particle;
   AliTPChit  *tpcHit;
   AliITShit  *itsHit;

// Get pointers to Alice detectors and Hits containers
   AliDetector *TPC  = gAlice->GetDetector("TPC");
   AliDetector *ITS  = gAlice->GetDetector("ITS");
   TClonesArray *Particles = gAlice->Particles();
   if (TPC) TClonesArray *TPChits   = TPC->Hits();
   if (ITS) TClonesArray *ITShits   = ITS->Hits();

   TTree *TH = gAlice->TreeH();
   Int_t ntracks    = TH->GetEntries();

   // Create histograms
   TH1F *hITSZ    = new TH1F("hITSZ","Z of hits in ITS",100,-50.,50.);
   TH1F *hITSR    = new TH1F("hITSR","R of hits in ITS",100,0.,50.);
   TH1F *hITSlayer= new TH1F("hITSlayer","JLAY of hits in ITS",6, 0., 6.);
   TH1F *hITSTr   = new TH1F("hITSTr","Track number for hits in ITS",100,0.,50000.);

// Start loop on tracks in the hits containers
   for (Int_t track=0; track<ntracks;track++) {
     gAlice->ResetHits();
     nbytes += TH->GetEvent(track);
     Int_t i=0;
     if (ITS) {
       nhits = ITShits->GetEntriesFast();
       for (hit=0;hit<nhits;hit++) {
	 itsHit   = (AliITShit*)ITShits->UncheckedAt(hit);
	 // With this new version, to be able to do proper detector
	 // simulations, the requirment that a "hit" leave some
	 // energy deposited has been removed.
	 if(itsHit->GetIonization()<=0.0) continue;
         ipart    = itsHit->fTrack;
         particle = (TParticle*)Particles->UncheckedAt(ipart);
         z = itsHit->GetZG();
         hITSZ->Fill(z);
         r = sqrt(itsHit->GetXG()*itsHit->GetXG() + 
		  itsHit->GetYG()*itsHit->GetYG());
         hITSR->Fill(r);
         hITSlayer->Fill((Float_t)itsHit->GetLayer());
         hITSTr->Fill(itsHit->fTrack);
         i++;
       }
     }        
   }

//Create a canvas, set the view range, show histograms
   TCanvas *c1 = new TCanvas("c1","Alice TPC and ITS hits",400,10,600,700);
   c1->Divide(2,2);
   c1->cd(1);
   gPad->SetFillColor(33);
   hITSZ->SetFillColor(33);
   hITSZ->Draw();
   c1->cd(2);
   gPad->SetFillColor(33);
   hITSR->SetFillColor(46);
   hITSR->Draw();
   c1->cd(3);
   gPad->SetFillColor(33);
   hITSlayer->SetFillColor(46);
   hITSlayer->Draw();
   c1->cd(4);
   gPad->SetFillColor(33);
   hITSTr->SetFillColor(46);
   hITSTr->Draw();
   c1->Print("analITS.gif");
}
