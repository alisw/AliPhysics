void analITSgeom (const char *filename="galice.root",Int_t evNumber=0){
/////////////////////////////////////////////////////////////////////////
//   This macro is a small example of a ROOT macro
//   illustrating how to read the output of GALICE
//   and fill some histograms.
//   
//     Root > .L analITSgeom.C   //this loads the macro in memory
//     Root > analITSgeom();     //by default process first event   
//     Root > analITgeomS(2);    //process third event
//Begin_Html
/*
<img src="figures/analITSgeom_ref.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>A reference plot produced by analITSgeom.C.
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
   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
   if (!file) file = new TFile(filename);

// Get AliRun object from file or create it if not on file
   if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }
      
// Import the Kine and Hits Trees for the event evNumber in the file
   Int_t nparticles = gAlice->GetEvent(evNumber);
   if (nparticles <= 0) return;
   Float_t x,y,z,mass,e,r,phi,flad;
   Int_t nbytes = 0;
   Int_t j,hit,ipart,lay;
   Int_t nhits;
   Int_t sector,plane;
   TParticle *particle;
   AliITShit  *itsHit;

// Get pointers to Alice detectors and Hits containers
   AliDetector *ITS  = gAlice->GetDetector("ITS");
   if(!ITS) return;
   TClonesArray *Particles = gAlice->Particles();
   TClonesArray *ITShits   = ITS->Hits();
   AliITSgeom   *gm        = ((AliITS *)ITS)->GetITSgeom();

   TTree *TH        = gAlice->TreeH();
   Int_t ntracks    = TH->GetEntries();

   // Create histograms
   Float_t pi2 = 2.0*TMath::Pi();
   Int_t   Nlad = gm->GetNladders(1);
   TH2F *hITS1 = new TH2F("hITS1","Ladder# vs. angle for Layer 1",
			  Nlad,0.5,(Float_t)Nlad+0.5,100,0.0,pi2);
   Int_t   Nlad = gm->GetNladders(2);
   TH2F *hITS2 = new TH2F("hITS2","Ladder# vs. angle for Layer 2",
			  Nlad,0.5,(Float_t)Nlad+0.5,100,0.0,pi2);
   Int_t   Nlad = gm->GetNladders(3);
   TH2F *hITS3 = new TH2F("hITS3","Ladder# vs. angle for Layer 3",
			  Nlad,0.5,(Float_t)Nlad+0.5,100,0.0,pi2);
   Int_t   Nlad = gm->GetNladders(4);
   TH2F *hITS4 = new TH2F("hITS4","Ladder# vs. angle for Layer 4",
			  Nlad,0.5,(Float_t)Nlad+0.5,100,0.0,pi2);
   Int_t   Nlad = gm->GetNladders(5);
   TH2F *hITS5 = new TH2F("hITS5","Ladder# vs. angle for Layer 5",
			  Nlad,0.5,(Float_t)Nlad+0.5,100,0.0,pi2);
   Int_t   Nlad = gm->GetNladders(6);
   TH2F *hITS6 = new TH2F("hITS6","Ladder# vs. angle for Layer 6",
			  Nlad,0.5,(Float_t)Nlad+0.5,100,0.0,pi2);

// Start loop on tracks in the hits containers
   for (Int_t track=0; track<ntracks;track++) {
     gAlice->ResetHits();
     nbytes += TH->GetEvent(track);
     Int_t i=0;
       nhits = ITShits->GetEntriesFast();
       for (hit=0;hit<nhits;hit++) {
	 itsHit   = (AliITShit*)ITShits->UncheckedAt(hit);
	 // With this new version, to be able to do proper detector
	 // simulations, the requirment that a "hit" leave some
	 // energy deposited has been removed.
	 if(itsHit->GetIonization()<=0.0) continue;
	 phi = TMath::ATan2(itsHit->GetYG(),itsHit->GetXG());
	 if(phi<0.0) phi += pi2;
	 lay = itsHit->GetLayer();
	 flad = ((Float_t)itsHit->GetLadder()) + 0.5;
         if(lay==1) hITS1->Fill(flad,phi,1.0);
         if(lay==2) hITS2->Fill(flad,phi,1.0);
         if(lay==3) hITS3->Fill(flad,phi,1.0);
         if(lay==4) hITS4->Fill(flad,phi,1.0);
         if(lay==5) hITS5->Fill(flad,phi,1.0);
         if(lay==6) hITS6->Fill(flad,phi,1.0);
       } // end for hit
   }  // end for track

//Create a canvas, set the view range, show histograms
   TCanvas *c1 = new TCanvas("c1"," ITS geometry test",400,10,600,700);
   c1->Divide(2,3);
   c1->cd(1);
   hITS1->Draw();
   c1->cd(2);
   hITS2->Draw();
   c1->cd(3);
   hITS3->Draw();
   c1->cd(4);
   hITS4->Draw();
   c1->cd(5);
   hITS5->Draw();
   c1->cd(6);
   hITS6->Draw();
   c1->Print("analITSgeom.gif");
}
