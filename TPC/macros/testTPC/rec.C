void rec(const char *filename="data.root")
{

  char *ocdbpath = gSystem->Getenv("OCDB_PATH");
  if (ocdbpath==0){
    ocdbpath="alien://folder=/alice/data/2007/LHC07w/OCDB/";
  }
  printf("OCDB PATH = %s\n",ocdbpath);
   

  //gSystem->Load("libXrdClient.so");
  //gSystem->Load("libNetx.so"); 
  AliLog::SetClassDebugLevel("AliTPCRawStream",-5);
  AliLog::SetClassDebugLevel("AliRawReaderDate",-5);
  AliLog::SetClassDebugLevel("AliTPCAltroMapping",-5);
  AliLog::SetModuleDebugLevel("RAW",-5);
  AliLog::SetGlobalLogLevel(3);

  //
  // First version of the reconstruction
  // script for the FDR'07

  // Set the CDB storage location
  // AliLog::SetModuleDebugLevel("STEER",2);
  AliCDBManager * man = AliCDBManager::Instance();
  //man->SetDefaultStorage("alien://folder=/alice/data/2007/LHC07w/OCDB/");
  man->SetDefaultStorage(ocdbpath);
  //  man->SetSpecificStorage("TPC/Calib/Parameters","local:///data/test2007/");
  // man->SetSpecificStorage("TPC/Calib/PadNoise","local:///data/test2007/");
  //   man->SetSpecificStorage("ITS/Calib/DDLMapSDD","local://$ALICE_ROOT/OCDB");
  //   man->SetSpecificStorage("MUON/Calib/Mapping","local://$ALICE_ROOT/OCDB");
  //   man->SetSpecificStorage("MUON/Calib/DDLStore","local://$ALICE_ROOT/OCDB");

  // TPC settings
  AliLog::SetClassDebugLevel("AliTPCclustererMI",2);
  AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetCosmicTestParam(kFALSE);
  tpcRecoParam->SetTimeInterval(60,940);
  tpcRecoParam->Dump();
  AliTPCReconstructor::SetRecoParam(tpcRecoParam);
  AliTPCReconstructor::SetStreamLevel(1);

  // PHOS settings
  AliPHOSRecoParam* recEmc = new AliPHOSRecoParamEmc();
  recEmc->SetSubtractPedestals(kTRUE);
  recEmc->SetMinE(0.05);
  recEmc->SetClusteringThreshold(0.10);
  AliPHOSReconstructor::SetRecoParamEmc(recEmc);

  // T0 settings
  AliLog::SetModuleDebugLevel("T0", 10);

  // MUON settings
  AliLog::SetClassDebugLevel("AliMUONRawStreamTracker",3);
 
  // AliReconstruction settings
  AliReconstruction rec;
  rec.SetUniformFieldTracking(kFALSE);  
  rec.SetWriteESDfriend(kTRUE);
  rec.SetWriteAlignmentData();
  rec.SetInput(filename);
  //
  //rec.SetRunLocalReconstruction("");
  rec.SetRunReconstruction("TPC");
  rec.SetFillESD("TPC");
  rec.SetRunV0Finder(kFALSE);  
  rec.SetRunVertexFinder(kFALSE);

  rec.SetRunQA(kTRUE);
  
  // Detector options if any
  rec.SetOption("MUON","SAVEDIGITS");
  rec.SetOption("TPC","OldRCUFormat");
  rec.SetOption("PHOS","OldRCUFormat");

  // To be enabled when CTP readout starts
  rec.SetFillTriggerESD(kFALSE);

  // all events in one single file
  rec.SetNumberOfEventsPerFile(-1);

  // switch off cleanESD
  rec.SetCleanESD(kFALSE);

  //AliLog::SetGlobalDebugLevel(2);
  rec.Run();

  cout << "-----------------------------------------------------------------" << endl;
  cout << "-----------------------------------------------------------------" << endl;
  cout << "--------- Reconstruction Completed. Start merging QAs -----------" << endl;
  cout << "-----------------------------------------------------------------" << endl;
  cout << "-----------------------------------------------------------------" << endl;
  //
  cout <<" EXITING RECONSTRUNCTION SESSION\n";
  //
  exit();
  AliTPCcalibDB::Instance()->GetParameters()->Dump();
  //  AliQADataMakerSteer qas;
  // qas.Merge();
}
