void ReadACORDEHits (Int_t numberOfEvents=-1,
		     const char *filename="galice.root")

  //  produces some plots from acorde hits

{
  // Dynamically link some shared libs
  if (gClassTable->GetID("AliRun") < 0) 
    {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
    }
  if (gAlice)
    {
      delete gAlice->GetRunLoader();
      delete gAlice;
      gAlice = 0x0;
    }

  // get loaders
  AliRunLoader *rl =
    AliRunLoader::Open("galice.root",
		       AliConfig::GetDefaultEventFolderName(),"read");
  if (!rl)
    {
      cerr<<"Can't load RunLoader from file! \n";      return 0x0;
    }
  
  rl->LoadgAlice();
  gAlice = rl->GetAliRun();
  if (!gAlice)
    {
      cerr << " AliRun object not found on file \n";     return 0x0; 
    }
  rl->LoadHeader();
  

  // Get the pointer to the ACORDE detector
  AliLoader *acordel = rl->GetLoader("ACORDELoader");
  AliACORDE *ACORDE = (AliACORDE *) gAlice->GetDetector("ACORDE");
  if (ACORDE == 0x0 || acordel == 0x0) {
    cerr << " Can not find ACORDE or ACORDELoader \n";    return 0x0;  
  }

  // get number of events 
  if (numberOfEvents<0) numberOfEvents=(Int_t)(rl->GetNumberOfEvents());
  
  //create histograms
  TH1F *elossH = new TH1F("eloss" ,"Energy Loss ",1000,0.,10);
  TH2F *geoH = new TH2F("geo"," ACORDE geometry seen by hits",
			250, -750,750, 250, -500,500);

  for (Int_t ievent=0; ievent<numberOfEvents; ievent++) {
    if ((ievent%10) == 0)  printf ("Processing event %d \n", ievent);
    rl->GetEvent(ievent);
    
    // Get the pointer Hit tree
    acordel->LoadHits();
    TTree *hitTree = acordel->TreeH();
    ACORDE->SetTreeAddress();
    if (!hitTree) {
      cout << " No TreeH found" << endl;      return 0x0; //rc;
    }
    
    rl->LoadKinematics();
    Int_t nTrack = (Int_t) hitTree->GetEntries();
    
    // Start loop on tracks in the hits containers
    for(Int_t iTrack=0; iTrack<nTrack;iTrack++){
      ACORDE->ResetHits();
      hitTree->GetEvent(iTrack);          
      if(ACORDE) {	 
	for(acordeHit=(AliACORDEhit*)ACORDE->FirstHit(-1);acordeHit;
	    acordeHit=(AliACORDEhit*)ACORDE->NextHit()) {
	  elossH->Fill( (Float_t)(acordeHit->Eloss()*1000.0));
	  geoH->Fill(acordeHit->X(),acordeHit->Z());
	} // end for acordeHits
      } // end if ACORDE
    } // end for iTrack  
  }
 // save histos in a root file
  TFile *fout = new TFile("ACORDE_hits.root","RECREATE");
  elossH->Write();
  geoH->Write();
  fout->Close();
}
