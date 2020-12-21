void rec() {

  AliReconstruction reco;

    reco.SetRunVertexFinder(kFALSE) ;
    reco.SetRunV0Finder(kFALSE) ;
    reco.SetRunMultFinder(kFALSE) ;

   AliQAManager::QAManager(AliQAv1::kSIMMODE) ;


  reco.SetCDBSnapshotMode("Sim/OCDBrec.root");
  reco.SetRunReconstruction("PHOS") ;
  reco.SetRunQA(":") ;
  reco.SetRunGlobalQA(kFALSE) ;

  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
