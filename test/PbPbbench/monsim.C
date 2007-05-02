void monsim(Int_t nev=1){ 
  // MonaLisa monitoring
  gSystem->Load("libNet.so");
  gSystem->Load("libMonaLisa.so");

  new TMonaLisaWriter(gSystem->Getenv("TEST_PLATFORMID"),"Simulation PbPb","aliendb3.cern.ch");

  gROOT->LoadMacro("sim.C");
  sim(nev);
  gMonitoringWriter->SendProcessingProgress(1,1,kTRUE);  

  // Send the size of the raw.root file

  FileStat_t buf;
  gSystem->GetPathInfo("./raw.root",buf);

  TList *valuelist = new TList();
  valuelist->SetOwner(kTRUE);

  TMonaLisaValue* valdouble = new TMonaLisaValue("raw.root size",buf.fSize);
  valuelist->Add(valdouble);

  gMonitoringWriter->SendParameters(valuelist);
  delete valuelist;

  printf("#Test finished successfully#\n");
} 
