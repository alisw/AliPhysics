void rec() {
  
  //  AliLog::SetClassDebugLevel("AliReconstruction",1);
  AliLog::SetClassDebugLevel("AliITSUReconstructor",1);

  TDatime t;

  gSystem->Load("libITSUpgradeBase");
  gSystem->Load("libITSUpgradeRec");

  gSystem->Exec("rm -rf *RecPoints* AliESD*");

  // Set ITS upgrade reconstructor
  gPluginMgr->AddHandler("AliReconstructor", "*",
			 "AliITSUReconstructor","ITS", "AliITSUReconstructor()");
  
  AliReconstruction rec;

  rec.SetRunReconstruction("ITS"); // run cluster finder
  rec.SetRunTracking(""); // Turn on with ITS when tracker is implemented

  rec.SetRunVertexFinder(kFALSE); // to be implemented - CreateVertexer
  rec.SetRunMultFinder(kFALSE);   // to be implemented - CreateMultFinder
  rec.SetRunPlaneEff(kFALSE);     // to be implemented - CreateTrackleter

  //  rec.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  rec.SetSpecificStorage("GRP/GRP/Data",
			 Form("local://%s",gSystem->pwd()));
  rec.SetSpecificStorage("ITS/Align/Data",
			 Form("local://%s",gSystem->pwd()));
  

  rec.SetRunQA(":");
  rec.SetRunGlobalQA(0);
  AliLog::Flush();

  TStopwatch timer;
  timer.Start();
  rec.Run();
  timer.Stop();
  timer.Print();

  printf("\n\n\n TDatime \n");

  t.Print();
  t.Set();
  t.Print();
}
