RICHpatrec (Int_t evNumber1=0,Int_t evNumber2=0) {
// Dynamically link some shared libs
 
    if (gClassTable->GetID("AliRun") < 0) {
	gROOT->LoadMacro("loadlibs.C");
	loadlibs();
    }
    else {
      //delete gAlice;
      gAlice = 0;
    }
// Connect the Root Galice file containing Geometry, Kine and Hits
    
    TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
    if (file) file->Close(); 
    file = new TFile("galice.root","UPDATE");
    file->ls();
    
    printf ("I'm after Map \n");
    
// Get AliRun object from file or create it if not on file
    
    if (!gAlice) {
	gAlice = (AliRun*)file->Get("gAlice");
	if (gAlice) printf("AliRun object found on file\n");
	if (!gAlice) gAlice = new AliRun("gAlice","RICH reconstruction program");
    }else {
      delete gAlice;
      gAlice = (AliRun*)file->Get("gAlice");
      	if (gAlice) printf("AliRun object found on file\n");
	if (!gAlice) gAlice = new AliRun("gAlice","Alice test program");
    }
    
    printf ("I'm after gAlice \n");
    
    AliRICH *RICH  = (AliRICH*) gAlice->GetDetector("RICH"); 
    
    // Create Recontruction algorithm object
    AliRICHPatRec *detect = new AliRICHPatRec("RICH patrec algorithm","Reconstruction");
    
// Reconstruct
    // Event Loop
    //
    for (int nev=0; nev<= evNumber2; nev++) {
      Int_t nparticles = gAlice->GetEvent(nev);
      cout <<endl<< "Processing event:" <<nev<<endl;
      cout << "Particles       :" <<nparticles<<endl;
      if (nev < evNumber1) continue;
      if (nparticles <= 0) return;
      if (RICH) detect->PatRec();
      char hname[30];
      sprintf(hname,"TreeR%d",nev);
      gAlice->TreeR()->Write(hname);
      gAlice->TreeR()->Reset();
   } // event loop  
   
    printf("\nEnd of Macro*************************************\n");
}



