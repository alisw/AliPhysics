void Rec(TString file="/Users/jklay/Projects/LHC/alice/work/beamtest07/Period_LHC07a_EMCAL.Run_000000518.Host_001.Seq_10.root")
{
  // Reconstruction of RAW data from the input file raw.root
  // Boris Polichtchouk, 31 Aug 2007

  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT");
  AliCDBManager::Instance()->SetSpecificStorage("EMCAL/Calib/Data","local://");

  //AliLog::SetGlobalDebugLevel(2);

  AliReconstruction rec ;
  rec.SetOption("EMCAL","OldRCUFormat");
  rec.SetRunTracking("") ;
  rec.SetRunVertexFinder(kFALSE) ; 
  rec.SetRunLocalReconstruction("EMCAL") ;
  rec.SetFillESD("EMCAL") ;

  rec.SetInput(file.Data());  // read RAW data

  rec.SetEventRange(0,10);
  rec.SetNumberOfEventsPerFile(-1);

  rec.Run();


}
