void Rec(TString file="/scratch/alicehp2/commun/testbeam07/LHC07a_EMCAL/000000190/raw/07000000190001.10.root")
{

  // Modified for aliroot v4-13
  // Cynthia Hadjidakis, 30 May 2008
  // Reconstruction of RAW data from the input file raw.root
  // Boris Polichtchouk, 31 Aug 2007


  AliCDBManager::Instance()->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  AliCDBManager::Instance()->SetSpecificStorage("EMCAL/Calib/Data","local://$ALICE_ROOT/OCDB/EMCAL/beamtest07/");

  //AliLog::SetGlobalDebugLevel(2);

  AliReconstruction rec ;
  rec.SetOption("EMCAL","OldRCUFormat");
  rec.SetRunQA(kFALSE);       // bug with QA
  rec.SetRunTracking("") ;
  rec.SetRunVertexFinder(kFALSE) ;
  rec.SetRunLocalReconstruction("EMCAL") ;
  rec.SetFillESD("EMCAL") ;
  
  rec.SetInput(file.Data());  // read RAW data

  rec.SetEventRange(0,-1);
  rec.SetNumberOfEventsPerFile(-1);
  rec.Run();

  delete rec;
  gObjectTable->Print();

}


