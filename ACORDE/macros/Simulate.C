void Simulate()
{
  AliSimulation sim;
  sim.SetConfigFile("./Config.C");
  //sim.SetMakeSDigits(" ");
  // sim.SetMakeDigits("");
  //sim.SetWriteRawData(" ");
  TStopwatch w;
  w.Start();
  sim.Run(10);
  w.Stop();
  w.Print();
}
