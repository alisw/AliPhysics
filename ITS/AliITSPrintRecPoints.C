void AliITSPrintRecPoints(TString rfn="galice.root",Int_t mod=-1,
			  Int_t evnt=-1){
  // Macro to print out the recpoints for all or a specific module

  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  } 
  else {
    if(gAlice){
      delete gAlice->GetRunLoader();
      delete gAlice;
      gAlice=0;
    }
  }

  gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSstandard.C");

  AliRunLoader *rl = AccessFile(rfn); // Set up to read in Data
  Int_t retval = rl->LoadHeader();
  if (retval){
    cerr<<"AliITSPrintRecPoints.C : LoadHeader returned error"<<endl;
    return;
  }

  AliITSLoader* ITSloader =  (AliITSLoader*) rl->GetLoader("ITSLoader");

  if(!ITSloader){
    cerr<<"AliITSPrintRecPoints.C :  ITS loader not found"<<endl;
    return;
  }

  ITSloader->LoadHits("read");
  ITSloader->LoadDigits("read");
  ITSloader->LoadRecPoints("read");
  AliITS *ITS = (AliITS*)gAlice->GetDetector("ITS");
  if(!ITS){
    cout << "Error: no ITS found. Aborting"<<endl;
    return;
  } // end if !ITS

  Int_t evNumber1 = 0;
  Int_t evNumber2 = gAlice->GetEventsPerRun();
  if(evnt>=0){
    evNumber1 = evnt;
    evNumber2 = evnt+1;
  } // end if evnt>=0
  Int_t mod1 = 0;
  Int_t mod2 = ITS->GetITSgeom()->GetIndexMax();
  if(mod>=0){
    mod1 = mod;
    mod2 = mod+1;
  } // end if mod>=0
  TClonesArray *rpa;
  AliITSRecPoint *rp = 0;
  AliITSDetTypeRec* rec = new AliITSDetTypeRec();
  rec->SetLoader(ITSloader);
  rec->SetITSgeom(ITS->GetITSgeom());
  rec->SetDefaults();

  Int_t event,m,i,i2;
  for(event = evNumber1; event < evNumber2; event++){
    rl->GetEvent(event);
    rec->SetTreeAddress();
    for(m=mod1;m<mod2;m++){
      rec->ResetRecPoints();
      TTree *TR = ITSloader->TreeR();
      TR->GetEvent(m);
      rpa = rec->RecPoints();
      i2 = rpa->GetEntriesFast();
      cout <<  "Event=" << event << " module=" << m <<
	" Number of Recpoints=" << i2 <<endl;
      for(i=0;i<i2;i++){
	rp = (AliITSRecPoint*)(rpa->At(i));
	cout << i << " ";
	rp->Print((ostream*)cout);
	cout << endl;
      } // end for i
    } // end for m
  } // end for event

}
