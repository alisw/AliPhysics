void FMDMerger (Int_t evNumber1=0, Int_t evNumber2=0, Int_t ibg=0, Int_t bgr=10) 
{
// Dynamically link some shared libs

   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   }

// Connect the Root Galice file containing Geometry, Kine and Hits

   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("pipe.root");
   if (file) file->Close(); 
   file = new TFile("galice.root","UPDATE");
   cout<<" Just starting... file "<<file<<endl;

// Get AliRun object from file or create it if not on file

   if (!gAlice) {
       gAlice = (AliRun*)file->Get("gAlice");
       if (gAlice) printf("AliRun object found on file\n");
       if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
   }
   AliFMD *FMD  = (AliFMD*) gAlice->GetModule("FMD");

   if (FMD) {
// creation
       AliFMDMerger* merger = new AliFMDMerger();
       cout<<" merger "<<merger<<endl;
       /*
       // granularity
	merger->SetRingsSi1(128);
        merger->SetRingsSi2(64);
        merger->SetSectorsSi1(20);
        merger->SetSectorsSi2(24);
       */
// configuration
       if (ibg) {
	 merger->SetMode(ibg);
	 merger->SetBackgroundFileName("bg.root");
	 cout<<" background"<<endl;
      }
       // pass
       FMD->SetMerger(merger);
   }
// Action !
//
//   Loop over events              
//
    for (int nev=evNumber1; nev<= evNumber2; nev++) {
	Int_t nparticles = gAlice->GetEvent(nev);
	cout << "nev         " << nev <<endl;
	cout << "nparticles  " << nparticles <<endl;
	if (nev < evNumber1) continue;
	if (nparticles <= 0) return;
	Int_t nbgr_ev = Int_t(nev*bgr/(evNumber2+1));
	cout<<" nbgr_ev "<<nbgr_ev<<endl;
	
	if (ibg) {
	    merger->SetBackgroundEventNumber(nbgr_ev);
	}

	gAlice->SDigits2Digits("FMD");

	char hname[30];
	sprintf(hname,"TreeD%d",nev);
	//	gAlice->TreeD()->Write(hname);
	//	cout<<hname<<" was written in file"<<file<<endl;
	//	gAlice->TreeD()->Print();
	//reset tree
	gAlice->TreeD()->Reset();

    }   // event loop 
}














