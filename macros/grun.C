void grun (Int_t nevent=1, const char *config="Config.C")
{
  //
  // Simple macro to run aliroot in a batch mode
  //
  gAlice->Init(config);
  TStopwatch timer;
  timer.Start();
//
//  If nevent is negative it is assumed that in config the 
//  global variable eventsPerRun has been set.
//
  if (nevent < 0) 
      nevent = eventsPerRun;
  gAlice->Run(nevent);
  timer.Stop();
  timer.Print();
}
