void sim(Int_t embrun) 
{
  AliSimulation sim;
  if (embrun == 2) {
    sim.SetRunGeneration(kFALSE);
    sim.SetMakeSDigits("");
  }
  else {
    sim.SetRunGeneration(kTRUE);
    sim.SetMakeSDigits("ITS TPC TRD TOF");
  }
  sim.SetRunSimulation(kTRUE);
  sim.SetMakeDigits("ITS TPC TRD TOF");
  sim.SetWriteRawData("ITS TPC TRD TOF","raw.root",kTRUE);
  if (embrun == 1)
    sim.MergeWith("../Background/galice.root",1);
  sim.SetDefaultStorage("alien://Folder=/alice/simulation/2008/v4-12-Release/Full/");
  sim.SetRunHLT("");
  sim.SetQA(kFALSE);

  sim.Run(1);
}
