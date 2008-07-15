void sim(Int_t nev=20) {
  gSystem->Load("libProof");
  gSystem->Load("libGui");
  if (!strcmp(gSystem->GetBuildArch(),"macosx")) gSystem->Load("libf95");
  gROOT->Macro("loadlibssim.C");
  new AliRun("gAlice","The ALICE Off-line Simulation Framework");

  AliSimulation simulator;
  simulator.SetMakeSDigits("TRD TOF PHOS HMPID EMCAL MUON FMD ZDC PMD T0 VZERO");
  simulator.SetMakeDigitsFromHits("ITS TPC");
/*   simulator.SetWriteRawData("ALL","raw.root",kTRUE); */
  simulator.SetRunHLT("");
  simulator.SetQA(kFALSE);
 
  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
