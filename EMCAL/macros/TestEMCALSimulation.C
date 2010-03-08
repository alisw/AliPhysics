/* 

Simple macro to test EMCAL Simulation

J.L. Klay
LLNL
 
*/

void TestEMCALSimulation(Int_t nev =10) {

  AliSimulation simulator;
  simulator.SetConfigFile("Config.C");
  simulator.SetMakeSDigits("EMCAL");
  simulator.SetMakeDigits("EMCAL");
  simulator.SetWriteRawData("EMCAL","raw.root",kTRUE);

  //OCDB settings
  simulator.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  simulator.SetSpecificStorage("GRP/GRP/Data",
                               Form("local://%s",gSystem->pwd()));

  //Avoid the HLT to run
  simulator.SetRunHLT("");

  TStopwatch timer;
  timer.Start();

  simulator.Run(nev);

  timer.Stop();
  timer.Print();

}
