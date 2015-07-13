/* 

Simple macro to test EMCAL Simulation

J.L. Klay
LLNL
 
*/

void TestEMCALSimulation(Int_t nev =10, Bool_t raw = kFALSE) {

  AliSimulation simulator;
  simulator.SetConfigFile("Config.C");
  simulator.SetMakeSDigits("EMCAL");
  simulator.SetMakeDigits ("EMCAL");

  //simulator.SetRunGeneration(kFALSE); // Generate or not particles
  //simulator.SetRunSimulation(kFALSE); // Generate or not HITS (detector response) or not, start from SDigits

  if(raw)  simulator.SetWriteRawData("EMCAL","raw.root",kTRUE);

  //OCDB settings
  simulator.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  simulator.SetSpecificStorage("GRP/GRP/Data",
                               Form("local://%s",gSystem->pwd()));

  // In case of anchoring MC, comment previous OCDB lines
  // select the appropriate year
  //simulator.SetDefaultStorage("alien://Folder=/alice/data/2011/OCDB");
  //simulator.UseVertexFromCDB();
  //simulator.UseMagFieldFromGRP();

  //Avoid the HLT to run
  simulator.SetRunHLT("");

  //Avoid QA
  simulator.SetRunQA(":");
  
  TStopwatch timer;
  timer.Start();

//  simulator.SetRunNumber(159582); // LHC11d run
  
  simulator.Run(nev);

  timer.Stop();
  timer.Print();

}
