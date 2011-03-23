void rec() {
  //AliLog::SetClassDebugLevel("AliITStrackerUpgrade",5);

  TDatime t;

  gSystem->Load("libITSUpgradeBase.so");
  gSystem->Load("libITSUpgradeRec.so");

  gSystem->Exec("rm -rf *RecPoints* AliESD*");
  AliITSRecoParam *p = AliITSRecoParam::GetLowFluxParam();
  p->SetTrackerSAOnly();
  p->SetInwardFindingSA();

  AliReconstruction rec;
  rec.SetRecoParam("ITS",p);
  rec.SetRunReconstruction("ITS");
  rec.SetUpgradeModule("ITS");
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunTracking("ITS");
  rec.SetFillESD("ITS");
  rec.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  rec.SetRunMultFinder(kFALSE);
  rec.SetSpecificStorage("GRP/GRP/Data",
			 Form("local://%s",gSystem->pwd()));
  rec.SetRunPlaneEff(kFALSE);
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
