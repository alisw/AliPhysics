void sim(Int_t nev=10) {
  AliSimulation simu;
  //simu.SetMakeSDigits("TRD TOF PHOS HMPID EMCAL FMD PMD T0 ZDC VZERO MUON HLT");
  //simu.SetMakeDigits ("TRD TOF PHOS HMPID EMCAL FMD PMD T0 ZDC VZERO MUON HLT");
  //simu.SetMakeDigitsFromHits("ITS TPC");
  //simu.SetWriteRawData("ALL","raw.root",kTRUE);
  simu.SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/?cacheFold=/tmp/CDBCache?operateDisconnected=kFALSE");

  // No write access to the OCDB => specific storage
  simu.SetSpecificStorage("GRP/GRP/Data",
                          Form("local://%s",gSystem->pwd()));

  //simu.SetRunQA(":") ; 
  //simu.SetRunHLT("");
  //simu.SetWriteGRPEntry(kFALSE);

  //Merge with HIJING
  //Load necessary libraries and connect to the GRID
  gSystem->Load("libNetx.so") ;
  gSystem->Load("libRAliEn.so");
  TGrid::Connect("alien://") ;
    
  //Feed Grid with collection file
  TGridCollection * collection = (TGridCollection*) TAlienCollection::Open("collection.xml");   
  if (! collection) {
    AliError(Form("collection: %s not found", kXML)) ;
    return kFALSE ;
  }
  TGridResult* result = collection->GetGridResult("",0 ,0);
  TString alienURL = result->GetKey(0, "turl") ;//split all files in collection, only one is sent per job
  cout << "================== " << alienURL << endl ;
  simu.MergeWith(alienURL);

  TStopwatch timer;
  timer.Start();
  simu.Run(nev);
  WriteXsection();
  timer.Stop();
  timer.Print();
}

WriteXsection()
{
  TPythia6 *pythia = TPythia6::Instance();
  pythia->Pystat(1);
  Double_t xsection = pythia->GetPARI(1);
  Int_t    ntrials  = pythia->GetMSTI(5);

  TTree   *tree   = new TTree("Xsection","Pythia cross section");
  TBranch *branch = tree->Branch("xsection", &xsection, "X/D");
  TBranch *branch = tree->Branch("ntrials" , &ntrials , "X/i");
  tree->Fill();

  TFile *file = new TFile("pyxsec.root","recreate");
  tree->Write();
  file->Close();

  cout << "Pythia cross section: " << xsection 
       << ", number of trials: " << ntrials << endl;
}
