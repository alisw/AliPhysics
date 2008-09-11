void sim(Int_t nev=4) {
  const char * kYear = "08" ; 
  AliSimulation simu;
  simu.SetMakeSDigits("TRD TOF PHOS HMPID  EMCAL MUON FMD PMD T0 ZDC VZERO");
  simu.SetMakeDigits ("TRD TOF PHOS HMPID  EMCAL MUON FMD PMD T0 ZDC VZERO");
  simu.SetMakeDigitsFromHits("ITS TPC");
  simu.SetWriteRawData("ALL","raw.root",kTRUE);
//  simu.SetDefaultStorage("alien://Folder=/alice/data/2008/LHC08d/OCDB/");
  simu.SetDefaultStorage("local://$ALICE_ROOT");
  simu.SetSpecificStorage("EMCAL/*","local://DB");

  simu.SetRunQA("ALL:ALL") ; 
 // AliQA::SetQARefStorage(Form("%s%s/", AliQA::GetQARefDefaultStorage(), kYear)) ;
  AliQA::SetQARefStorage("local://$ALICE_ROOT") ;
  // AliQA::SetQARefDataDirName(AliQA::kMONTECARLO) ; //RUN_TYPE
  for (Int_t det = 0 ; det < AliQA::kNDET ; det++) {
    simu.SetQACycles(det, 2) ;
  }
  
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

