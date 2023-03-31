void rec() {

  AliReconstruction reco;

    reco.SetRunVertexFinder(kFALSE) ;
    reco.SetRunV0Finder(kFALSE) ;
    reco.SetRunMultFinder(kFALSE) ;

   AliQAManager::QAManager(AliQAv1::kSIMMODE) ;


  reco.SetDefaultStorage("alien://Folder=/alice/data/2018/OCDB");
  reco.SetRunReconstruction("PHOS") ;
  reco.SetRunQA(":") ;
  reco.SetRunGlobalQA(kFALSE) ;

 reco.SetSpecificStorage("GRP/GRP/Data",
                          Form("local://%s",gSystem->pwd()));
  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
