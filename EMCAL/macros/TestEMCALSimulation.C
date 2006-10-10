/* 

Simple macro to test EMCAL Simulation

J.L. Klay
LLNL
 
*/

void TestEMCALSimulation(Int_t nev =1) {

  AliSimulation simulator;

  TStopwatch timer;
  timer.Start();
  simulator.SetConfigFile("Config.C");
  simulator.SetMakeSDigits("EMCAL");
  simulator.SetMakeDigits("EMCAL");
  simulator.Run(nev);
  timer.Stop();
  timer.Print();

}
