#ifndef __CINT__
#include "iostream.h"
#endif

Bool_t TPCSortTracks(Int_t event = 0)
{
   TFile *fileTracks   = TFile::Open("AliTPCtracks.root");
   TFile *fileClusters = TFile::Open("AliTPCclusters.root");
   TFile *fileEvent    = TFile::Open("galice.root");

   // get TPC parameterization
   AliTPCParam *param=(AliTPCParam *)fileEvent->Get("75x40_100x60_150x60");
   if (!param) {
		cerr << "(TPCSortTracks) ERROR: TPC parameters have not been found!" << endl;
		return kFALSE;
	}

   // read and sort tracks
   Int_t i;
	TSortedList tracks_list;
   AliTPCtrack *iotrack = 0;
   TTree *tracktree = (TTree*)fileTracks->Get(Form("TreeT_TPC_%d", event));
   Int_t nentr = (Int_t)tracktree->GetEntries();
   for (i = 0; i < nentr; i++) {
     iotrack = new AliTPCtrack;
     tracktree->SetBranchAddress("tracks", &iotrack);
     tracktree->GetEvent(i);
     tracks_list.Add(iotrack);
   }   
	delete tracktree;
	
   // assign to each track its GEANT label
   fileClusters->cd();
	AliTPCtracker *tracker = new AliTPCtracker(param, event);
   tracker->LoadInnerSectors();
   tracker->LoadOuterSectors();
	TListIter iter(&tracks_list);
   for (i = 0; i < nentr; i++) {
     iotrack = (AliTPCtrack*)iter.Next();
	  if (!iotrack) {
	     cerr << "(TPCSortTracks) WARNING: Track no. " << i << " is NULL!!!" << endl;
		  continue;  
	  }
     tracker->CookLabel(iotrack, 0.1);
   }   
   delete tracker;
	
   // create the new TTree of TPC tracks sorted w.r. to Pt
	tracktree = new TTree(Form("TreeT_TPC_%d", event),"Tree with TPC tracks sorted w.r to pt");
   tracktree->Branch("tracks", "AliTPCtrack", &iotrack, 32000, 0);
	iter.Reset();
   for (i = 0; i < nentr; i++) {
     iotrack = (AliTPCtrack*)iter.Next();
	  if (!iotrack) {
	     cerr << "(TPCSortTracks) WARNING: Track no. " << i << " is NULL!!!" << endl;
		  continue;  
	  }
     tracktree->Fill();
   }
	
	// save the new tree into new file
	TFile *fileOutput = TFile::Open("AliTPCtracksSorted.root","recreate");
   tracktree->Write();
   fileOutput->Close();
   fileEvent->Close();
   fileClusters->Close();
	fileTracks->Close();
	
	return kTRUE;
}


void AliITSTrackingV1(Int_t evNumber1=0,Int_t evNumber2=0, Int_t min_t=-1, Int_t max_t=0,Bool_t flagvert=1, Bool_t realmass=0, const char *filename="galice.root") {

  ///////////////// Dynamically link some shared libs ////////////////////////////////
  
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  } else {
    delete gAlice;
    gAlice=0;
  }
  
  // Connect the Root Galice file containing Geometry, Kine and Hits
   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
   if (!file) file = new TFile("galice.root","UPDATE");
   //if (!file) file = new TFile(filename);

// Get AliRun object from file or create it if not on file
   if (!gAlice) {
      gAlice = (AliRun*)file->Get("gAlice");
      if (gAlice) printf("AliRun object found on file\n");
      if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }
   
   AliKalmanTrack::SetMagneticField(gAlice->Field()->SolenoidField() / 10.0);

  
  cout << "Sorting TPC tracks w.r. to transverse momentum...";
  Bool_t success_sorting = TPCSortTracks();
  if (success_sorting) {
     cout << "DONE!" << endl;
  }
  else {
	cout << "Some error occurred..." << endl;
	return 1;
  }

  AliITS* IITTSS =(AliITS *)gAlice->GetDetector("ITS");        
  if (!IITTSS) return;                                           

//
//   Loop over events 
//
   Int_t Nh=0;
   Int_t Nh1=0;
   for (Int_t nev=0; nev<= evNumber2; nev++) {
	AliITSTrackerV1 *ITStracker = new AliITSTrackerV1(IITTSS,nev,flagvert);
     Int_t nparticles = gAlice->GetEvent(nev);
     cout << "nev         " << nev <<endl;
     cout << "nparticles  " << nparticles <<endl;
     if (nev < evNumber1) continue;
     if (nparticles <= 0) return;

     TTree *TR=gAlice->TreeR();
     Int_t nent=TR->GetEntries();
     //printf("Found %d entries in the TreeR (must be one per module per event!)\n",nent);


     TStopwatch timer;
	  
	  timer.Start();
     ITStracker->DoTracking(nev,min_t,max_t,file,realmass);    
     timer.Stop(); timer.Print();
	 AliITSgeom *g1 = IITTSS->GetITSgeom();
    Int_t NumOfModules = g1->GetIndexMax();	   
	  ITStracker->DelMatrix(NumOfModules);
	  delete ITStracker;  
   }   // event loop  
   file->Close();   
}

