/* 

Simple macro to test EMCAL Simulation

J.L. Klay
LLNL
 
*/

void TestEMCALSimulation(Int_t nev =10, Bool_t raw = kFALSE) {

  AliSimulation simulator;
  simulator.SetConfigFile("Config.C");
  simulator.SetMakeSDigits("EMCAL");
  simulator.SetMakeDigits("EMCAL");
  if(raw)  simulator.SetWriteRawData("EMCAL","raw.root",kTRUE);

  //OCDB settings
  simulator.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  simulator.SetSpecificStorage("GRP/GRP/Data",
                               Form("local://%s",gSystem->pwd()));

  //Avoid the HLT to run
  simulator.SetRunHLT("");

  TStopwatch timer;
  timer.Start();

//  simulator.SetRunNumber(140234);
  
  simulator.Run(nev);

  timer.Stop();
  timer.Print();

}
