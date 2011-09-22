void rec() {

  AliReconstruction reco;

//
// switch off cleanESD, write ESDfriends and Alignment data
  
//  reco.SetCleanESD(kFALSE);
//  reco.SetWriteESDfriend();
//  reco.SetWriteAlignmentData();
    reco.SetRunVertexFinder(kFALSE) ;
    reco.SetRunV0Finder(kFALSE) ;
    reco.SetRunMultFinder(kFALSE) ;

   AliQAManager::QAManager(AliQAv1::kSIMMODE) ;

//
// ITS Efficiency and tracking errors

//
// Residual OCDB

  reco.SetDefaultStorage("alien://folder=/alice/data/2010/OCDB");
//  reco.SetDefaultStorage("local://./OCDB");
  reco.SetSpecificStorage("PHOS/*/*","local://./OCDB");


//  reco.SetRunReconstruction("PHOS EMCAL") ;
  reco.SetRunReconstruction("PHOS") ;
  reco.SetRunQA(":") ;
  reco.SetRunGlobalQA(kFALSE) ;
//
// GRP from local OCDB

 reco.SetSpecificStorage("GRP/GRP/Data",
                          Form("local://%s",gSystem->pwd()));
//  reco.SetSpecificStorage("GRP/GRP/Data", "alien://Folder=/alice/data/2010/OCDB");


//
  TStopwatch timer;
  timer.Start();
  reco.Run();
  timer.Stop();
  timer.Print();
}
