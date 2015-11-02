
void rec() {
  
  AliLog::SetClassDebugLevel("AliITSUReconstructor",1);

  // Set ITS upgrade reconstructor
  gPluginMgr->AddHandler("AliReconstructor", "ITS",
			 "AliITSUReconstructor","ITS", "AliITSUReconstructor()");
  
  AliReconstruction rec;
  rec.SetRunReconstruction("ITS"); // run cluster finder
  rec.SetRunTracking("ITS"); 
  
  
  rec.SetRunVertexFinder(kTRUE);  // to be implemented - CreateVertexer
  rec.SetRunMultFinder(kFALSE);   // to be implemented - CreateMultFinder
  rec.SetRunPlaneEff(kFALSE);     // to be implemented - CreateTrackleter

  rec.SetSpecificStorage("GRP/GRP/Data",
			 Form("local://%s",gSystem->pwd()));
  rec.SetSpecificStorage("ITS/Align/Data",
			 Form("local://%s",gSystem->pwd()));
  rec.SetSpecificStorage("ITS/Calib/RecoParam",
			 Form("local://%s",gSystem->pwd()));
  

  rec.SetRunQA(":");
  rec.SetRunGlobalQA(0);
  AliLog::Flush();

  AliITSURecoParam *par=AliITSURecoParam::GetHighFluxParam();
  par->SetTracker(2);    // 1 is the Cooked Matrix tracker, 2 is the Cellular Automaton tracker  
  par->SetSAonly(kTRUE); // kFALSE is the TPC+ITS mode
  rec.SetRecoParam("ITS",par);

  TStopwatch timer;
  timer.Start();
  //
  rec.Run();
  timer.Stop();
  timer.Print();
}

