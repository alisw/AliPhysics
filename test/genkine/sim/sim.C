void sim(Int_t nev=1) {
  AliSimulation simulator;

  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
