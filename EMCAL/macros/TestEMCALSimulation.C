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


  TStopwatch timer;
  timer.Start();

  simulator.Run(nev);

  timer.Stop();
  timer.Print();

}
