void sim(Int_t nev=150) {
  AliSimulation simu;
  simu.SetMakeSDigits("TRD TOF PHOS HMPID EMCAL MUON FMD ZDC T0 VZERO");
  simu.SetMakeDigitsFromHits("ITS TPC");

  simu.SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal/");
  
  // No write access to the OCDB => local specific storage
  simu.SetSpecificStorage("GRP/GRP/Data",
			  Form("local://%s/",gSystem->pwd()));

  Printf("%s:%d %d events",(char*)__FILE__,__LINE__,nev);

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
 UInt_t    ntrials  = pythia->GetMSTI(5);

 TFile *file = new TFile("pyxsec.root","recreate");
 TTree   *tree   = new TTree("Xsection","Pythia cross section");
 TBranch *branch = tree->Branch("xsection", &xsection, "X/D");
 TBranch *branch = tree->Branch("ntrials" , &ntrials , "X/i");
 tree->Fill();



 tree->Write();
 file->Close();

 cout << "Pythia cross section: " << xsection
      << ", number of trials: " << ntrials << endl;
}

