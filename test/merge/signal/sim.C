void sim(Int_t nev=1) {
  AliSimulation simulator;
  simulator.MergeWith("../backgr/galice.root",3);
  simulator.SetWriteRawData("ALL","raw.root",kTRUE);

  TStopwatch timer;
  timer.Start();
  simulator.Run(nev);
  timer.Stop();
  timer.Print();
}
