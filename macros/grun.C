void grun (Int_t nevent=1, const char *config="Config.C")
{
  //
  // Simple macro to run aliroot in a batch mode
  //
  gAlice->Init(config);
  TStopwatch timer;
  timer.Start();
  gAlice->Run(nevent);
  timer.Stop();
  timer.Print();
}
