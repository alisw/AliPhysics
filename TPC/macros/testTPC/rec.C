void rec(const char *filename="data.root")
{
  gSystem->Load("libXrdClient.so");
  gSystem->Load("libNetX.so");
  //
  // First version of the reconstruction
  // script for the FDR'07

  // Set the CDB storage location
  // AliLog::SetModuleDebugLevel("STEER",2);
  AliCDBManager * man = AliCDBManager::Instance();
  //man->SetDefaultStorage("alien://folder=/alice/data/2007/LHC07w/OCDB/");
  man->SetDefaultStorage("/data/test2007/OCDB");
  man->SetSpecificStorage("GRP/GRP/Data","local://$ALICE_ROOT");
  man->SetSpecificStorage("ITS/Calib/DDLMapSDD","local://$ALICE_ROOT");
  man->SetSpecificStorage("MUON/Calib/Mapping","local://$ALICE_ROOT");
  man->SetSpecificStorage("MUON/Calib/DDLStore","local://$ALICE_ROOT");

  // TPC settings
  AliLog::SetClassDebugLevel("AliTPCclustererMI",2);
  AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetCosmicTestParam(kTRUE);
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
 
  // Tracking settings
  //  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., 1);
  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 0., 10., 2);
  AliTracker::SetFieldMap(field,1);

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
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunQA(kFALSE);
  //
  //rec.SetEventRange(1,5);

  // In case some detectors have to be switched off...
  //  rec.SetRunLocalReconstruction("ALL");
  //  rec.SetRunTracking("ALL");
  //  rec.SetFillESD("ALL");
  // Disable vertex finder for the moment
  //  rec.SetRunVertexFinder(kFALSE);

  // To be enabled if some equipment IDs are not set correctly by DAQ
  //  rec.SetEquipmentIdMap("EquipmentIdMap.data");

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
  AliQADataMakerSteer qas;
  qas.Merge();
}
