// Simulation steering macro

// mikael.mieskolainen@cern.ch, 22.9.2015

void sim(Int_t nev = 3) {

  if (gSystem->Getenv("EVENT")) {
    nev = atoi(gSystem->Getenv("EVENT"));
  }

  AliSimulation simulator;
  simulator.SetMakeSDigits("TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0 VZERO");
  simulator.SetMakeDigitsFromHits("ITS TPC");
  simulator.SetWriteRawData("ALL", "raw.root", kTRUE);

  simulator.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  simulator.SetSpecificStorage("GRP/GRP/Data", Form("local://%s", gSystem->pwd()));

  // Quality Assurance
  simulator.SetRunQA("ALL:ALL");
  simulator.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;

  for (Int_t det = 0; det < AliQA::kNDET; ++det) {
    simulator.SetQACycles((AliQAv1::DETECTORINDEX_t)det, nev+1) ;
  }

  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
