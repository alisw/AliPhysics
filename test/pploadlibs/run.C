void run() 
{
  TStopwatch timer;
  timer.Start();
  gSystem->Exec("root.exe -q sim.C > sim.log 2>&1");
  gSystem->Exec("root.exe -q rec.C > rec.log 2>&1");
  timer.Stop();
  timer.Print();
}
