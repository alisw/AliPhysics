// $Id$
//
// Macro for running lego simulation in test/vmctest/lego.

void sim(Int_t nev=1, const TString& config, const TString& det)  {
  AliSimulation simulator(config);
  TString configFunction("Config(\"");
  configFunction.Append(det.Data());
  configFunction.Append(TString("\");"));
  cout << configFunction.Data() << endl;
  gAlice->SetConfigFunction(configFunction.Data());
  TStopwatch timer;
  timer.Start();
  simulator.RunLego(config.Data());
  timer.Stop();
  timer.Print();
}
