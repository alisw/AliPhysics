void sim(Int_t nev=100) {
  AliSimulation simulator;
  simulator.SetMakeSDigits("TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD START VZERO");
  simulator.SetMakeDigitsFromHits("ITS TPC");
  simulator.SetWriteRawData("ALL","raw.root",kTRUE);
  simulator.SetDefaultStorage("alien://Folder=/alice/data/2006/LHC06a/CDB");

  TGrid::Connect("alien://");
  if (!gGrid) {
 	fprintf(stderr,"Error: cannot authenticate to the API service - exit!\n");
	exit(-2);
  }	

  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
