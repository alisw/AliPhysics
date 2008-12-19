void AliITSPrintHits(TString hfn="galice.root",Int_t mod=-1,
		     Int_t evnt=-1){
  // Macro to print out the recpoints for all or a specific module

  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  } 
  else {
    if(gAlice){
      delete AliRunLoader::GetRunLoader();
      delete gAlice;
      gAlice=0;
    }
  }
  gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSstandard.C");

  AliRunLoader *rl = AccessFile(hfn); // Set up to read in Data
  Int_t retval = rl->LoadHeader();
  if (retval){
    cerr<<"AliITSPrintHits.C : LoadHeader returned error"<<endl;
    return;
  }

  AliITSLoader* ITSloader =  (AliITSLoader*) rl->GetLoader("ITSLoader");

  if(!ITSloader){
    cerr<<"AliITSPrintHits.C :  ITS loader not found"<<endl;
    return;
  }

  ITSloader->LoadHits("read");
  AliITS *ITS = (AliITS*)gAlice->GetDetector("ITS");
  if(!ITS){
    cout << "Error: no ITS found. Aborting"<<endl;
    return;
  } // end if !ITS

  Int_t evNumber1 = 0;
  Int_t evNumber2 = AliRunLoader::GetNumberOfEvents();
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
  AliITShit *hp = 0;

  Int_t nmodules,size=-1;
  Int_t event,m,i,i2,hit,trk;
  for(event = evNumber1; event < evNumber2; event++){
    cout<<"Processing event "<<event<<endl;
    rl->GetEvent(event);
    ITS->InitModules(size,nmodules);
    ITS->FillModules(event,0,-1," "," ");
    for(m=mod1;m<mod2;m++){
      i2 = (ITS->GetModule(m))->GetNhits();
      cout <<  "Event=" << event << " module=" << m <<
	" Number of Hits=" << i2 <<endl;
      for(i=0;i<i2;i++){
	trk = (ITS->GetModule(m))->GetHitTrackIndex(i);
	hit = (ITS->GetModule(m))->GetHitHitIndex(i);
	hp  = (ITS->GetModule(m))->GetHit(i);
	cout << i << " trk#="<<trk<<" hit#="<< hit << " ";
	hp->Print((ostream*)cout);
	cout << endl;
      } // end for i
    } // end for m
    ITS->ClearModules();
  } // end for event
}
