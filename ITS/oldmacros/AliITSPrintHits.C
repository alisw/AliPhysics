void AliITSPrintHits(Int_t track, Int_t hit=-1, Int_t evNumber=0
		     ,const char *filename="galice.root"){
/*
  This macro will print the a specific hit or all of the hits from a give 
track.
*/
  if(gAlice){
    delete gAlice;
    gAlice=0;
  }else{
    // Dynamically link some shared libs
    if(gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
    } // end if 
  }
// Connect the Root Galice file containing Geometry, Kine and Hits
    TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
    if(!file) file = new TFile(filename);
 
// Get AliRun object from file or create it if not on file
    if(!gAlice) {
        gAlice = (AliRun*)file->Get("gAlice");
        if(gAlice) printf("AliRun object found on file\n");
        if(!gAlice) gAlice = new AliRun("gAlice","Alice test program");
    } // end if !gAlice
 
// Set event pointer to this event
    Int_t nparticles = gAlice->GetEvent(evNumber);
    if (nparticles <= 0){
        cout << "No particles found for event " << evNumber;
        cout << " in file " << filename << endl;
        return;
    } // end if nparticles <=0 
// Pointer to specific detector hits.
    AliITShit    *itsHit; 
// Get pointers to ALL Alice detectors and Hits containers
    AliITS    *ITS    = (AliITS*)    gAlice->GetDetector("ITS"); 
    Int_t t,h,ntracks = gAlice->TreeH()->GetEntries();  
// Start loop on tracks in the hits containers
    for(t=0; t<ntracks;t++){
	gAlice->ResetHits();
	gAlice->TreeH()->GetEvent(t);
	h = 0;
	for(itsHit=(AliITShit*)ITS->FirstHit(-1);itsHit;
	    itsHit=(AliITShit*)ITS->NextHit()){
	    if((track==itsHit->GetTrack())&&(h==hit||hit<0)){
		cout << h << " ";
		itsHit->Print((ostream*)cout);
		cout << endl;
	    } // end if
	    h++;
	} // end for itsHit
    } // end for t
}
