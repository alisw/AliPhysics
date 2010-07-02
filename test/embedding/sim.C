void sim(Int_t embrun) 
{
  AliSimulation sim;
  if (embrun == 4) {
    AliCDBManager *cdbm = AliCDBManager::Instance();
    cdbm->SetRun(atoi(gSystem->Getenv("DC_RUN")));
    cdbm->SetDefaultStorage("local://$ALICE_ROOT/OCDB");     
    cdbm->SetSpecificStorage("GRP/GRP/Data",Form("local://%s",gSystem->pwd()));
    sim.SetMakeSDigits("ITS TPC TRD TOF VZERO");  

    sim.ConvertRaw2SDigits("raw.root","AliESDs.root");
    return;
  }
  
  if (embrun == 2) {
    sim.SetRunGeneration(kFALSE);
    sim.SetMakeSDigits("");
    sim.SetMakeDigitsFromHits("");
  }
  else {
    sim.SetRunGeneration(kTRUE);
    sim.SetMakeSDigits("ITS TPC TRD TOF VZERO");
  }
  sim.SetRunSimulation(kTRUE);
  sim.SetMakeDigits("ITS TPC TRD TOF VZERO");
  sim.SetWriteRawData("ITS TPC TRD TOF VZERO","raw.root",kTRUE);
  if (embrun == 1)
    sim.MergeWith("../BackgroundSDigits/galice.root",1);

  sim.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  sim.SetSpecificStorage("GRP/GRP/Data",
			       Form("local://%s",gSystem->pwd()));
  sim.SetRunQA(":") ; 
  AliQA::SetQARefStorage("local://$ALICE_ROOT/OCDB") ;
  
  for (Int_t det = 0 ; det < AliQA::kNDET ; det++) {
    sim.SetQACycles(det, 1) ;
  }

//   sim.SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-15-Release/Full/");
//   sim.SetRunHLT("");
//   sim.SetQA(kFALSE);

  sim.Run(1);
}
