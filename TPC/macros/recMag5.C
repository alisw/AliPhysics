void recMag5(const char *filename="data.root")
{
  /////////////////////////////////////////////////////////////////////////////////////////
  //
  // First version of the reconstruction
  // script for the FDR'08

  // Set the CDB storage location
  // AliLog::SetModuleDebugLevel("STEER",2);
  AliCDBManager * man = AliCDBManager::Instance();
  //  man->SetDefaultStorage("local://LocalCDB");
  man->SetDefaultStorage("alien://folder=/alice/data/2008/LHC08b/OCDB/"); 
  //man->SetDefaultStorage("local:///data/OCDB");
  
  // Files that we can not read from alien...solved
  //  man->SetSpecificStorage("ITS/Calib/MapsAnodeSDD","local://$ALICE_ROOT/OCDB");
  //  man->SetSpecificStorage("ITS/Calib/MapsTimeSDD","local://$ALICE_ROOT/OCDB");
  //  man->SetSpecificStorage("TPC/Calib/ExB","local://$ALICE_ROOT/OCDB");

  // Objects not found if using LHC07w database...solved
  //  man->SetSpecificStorage("ITS/Calib/MapsAnodeSDD","local:///afs/cern.ch/user/c/cheshkov/public/OCDB");
  // man->SetSpecificStorage("GRP/GRP/Data","local://$ALICE_ROOT/OCDB");
  // man->SetSpecificStorage("ITS/Calib/DDLMapSDD","local://$ALICE_ROOT/OCDB");
  // man->SetSpecificStorage("MUON/Calib/Mapping","local://$ALICE_ROOT/OCDB");
  // man->SetSpecificStorage("MUON/Calib/DDLStore","local://$ALICE_ROOT/OCDB");
  
  AliITSRecoParam * itsRecoParam =  AliITSRecoParam::GetCosmicTestParam();
  itsRecoParam->SetFactorSAWindowSizes(20);
  itsRecoParam->SetClusterErrorsParam(2);
  itsRecoParam->SetFindV0s(kFALSE);
  itsRecoParam->SetAddVirtualClustersInDeadZone(kFALSE);
  itsRecoParam->SetUseAmplitudeInfo(kFALSE);
  // In case we want to switch off a layer
  //  itsRecoParam->SetLayerToSkip(<N>);
  // itsRecoParam->SetLayerToSkip(4);
  // itsRecoParam->SetLayerToSkip(5);
 itsRecoParam->SetLayerToSkip(2);
 itsRecoParam->SetLayerToSkip(3);
 itsRecoParam->SetClusterMisalError(1.0); // [cm]
 itsRecoParam->SetSAUseAllClusters();
 AliITSReconstructor::SetRecoParam(itsRecoParam); 
 
  // TPC settings
  AliLog::SetClassDebugLevel("AliTPCclusterer",2);
  AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetCosmicTestParam(kFALSE);
  Double_t sysError[5]={0.3,3, 0.3/150., 3./150.,0.3/(150*150.)};
  tpcRecoParam->SetSystematicError(sysError);
  tpcRecoParam->SetTimeInterval(60,940);
  tpcRecoParam->Dump();
  AliTPCReconstructor::SetRecoParam(tpcRecoParam);
  AliTPCReconstructor::SetStreamLevel(10);


  // TRD setting
  AliTRDrawStreamBase::SetRawStreamVersion("TB");

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
  AliMUONRecoParam *muonRecoParam = AliMUONRecoParam::GetLowFluxParam();
  muonRecoParam->CombineClusterTrackReco(kTRUE);
  muonRecoParam->SetCalibrationMode("NOGAIN");
  //muonRecoParam->SetClusteringMode("PEAKFIT");
  //muonRecoParam->SetClusteringMode("PEAKCOG");
  muonRecoParam->Print("FULL");
  AliRecoParam::Instance()->RegisterRecoParam(muonRecoParam);
 
  // Tracking settings
  Double_t mostProbPt=0.35;
  AliExternalTrackParam::SetMostProbablePt(mostProbPt);

  // AliReconstruction settings
  AliReconstruction rec;
  rec.SetUniformFieldTracking(kFALSE);
  rec.SetWriteESDfriend(kTRUE);
  rec.SetWriteAlignmentData();
  rec.SetInput(filename);
  rec.SetRunReconstruction("ALL");
  rec.SetUseTrackingErrorsForAlignment("ITS");

  // In case some detectors have to be switched off...
  //  rec.SetRunLocalReconstruction("ALL");
  //  rec.SetRunTracking("ALL");
  //  rec.SetFillESD("ALL");
  // Enable vertex finder - it is needed for cosmic track reco
  rec.SetRunVertexFinder(kFALSE);

  // To be enabled if some equipment IDs are not set correctly by DAQ
  //  rec.SetEquipmentIdMap("EquipmentIdMap.data");

  // Detector options if any
  //  rec.SetOption("ITS","cosmics,onlyITS");
  //rec.SetOption("ITS","cosmics,MI");
  rec.SetOption("ITS","cosmics");
  rec.SetOption("MUON","SAVEDIGITS");
  rec.SetOption("TPC","OldRCUFormat");
  rec.SetOption("PHOS","OldRCUFormat");
  rec.SetOption("T0","cosmic");

  // To be enabled when CTP readout starts
  rec.SetFillTriggerESD(kFALSE);

  // all events in one single file
  rec.SetNumberOfEventsPerFile(-1);

  // switch off cleanESD
  rec.SetCleanESD(kFALSE);

  // AliLog::SetGlobalDebugLevel(2);

  rec.SetRunQA(kFALSE);
  AliLog::Flush();
  //TPC addition 
  rec.SetRunReconstruction("ITS TPC");
  rec.SetFillESD("ITS TPC"); 
  //rec.SetEventRange(0,40);
  // filter logs
  AliLog::SetClassDebugLevel("AliTPCRawStream",-5);
  AliLog::SetClassDebugLevel("AliRawReaderDate",-5);
  AliLog::SetClassDebugLevel("AliTPCAltroMapping",-5);
  AliLog::SetModuleDebugLevel("RAW",-5);  
  // 
  rec.Run();

  //cout << "-----------------------------------------------------------------" << endl;
  //cout << "-----------------------------------------------------------------" << endl;
  //cout << "--------- Reconstruction Completed. Start merging QAs -----------" << endl;
  //cout << "-----------------------------------------------------------------" << endl;
  //cout << "-----------------------------------------------------------------" << endl;
  //AliQADataMakerSteer qas;
  //qas.Merge();
  //  exit(10);
}
