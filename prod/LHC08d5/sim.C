void sim(Int_t nev=10) {

  AliSimulation simu;
  simu.SetMakeSDigits("TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0 VZERO");
  simu.SetMakeDigitsFromHits("ITS TPC");

  simu.SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Ideal");

  // No write access to the OCDB => specific storage
  simu.SetSpecificStorage("GRP/GRP/Data",
			       Form("local://%s",gSystem->pwd()));

  TStopwatch timer;
  timer.Start();
  simu.Run(nev);
  timer.Stop();
  timer.Print();
}
