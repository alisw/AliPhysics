void rec(Int_t runNumber = 0, const char* year = "08", const char *localFileName = NULL)
{
  // Offline shifter reconstruction macro

  TString filename;

  if (!localFileName) {

    cout << "Going to run the reconstruction for run: " << runNumber << endl;

    // connect to the grid 
    TGrid * grid = 0x0 ; 
    grid = TGrid::Connect("alien://") ; 
		
    // make the file name pattern year and run number
    TString pattern;
    pattern.Form("%9d",runNumber);
    pattern.ReplaceAll(" ", "0") ; 
    pattern.Prepend(year);
    pattern.Append("*0.root");

    // find the files associated to this run
    // get the list of files from AliEn directly 
    TString baseDir; 
    baseDir.Form("/alice/data/20%s/",year);

    cout << "Looking for raw-data files with pattern " << pattern << " in folder " << baseDir << endl;

    TGridResult *result = grid->Query(baseDir, pattern);

    TList *fileList = result->GetFileInfoList();

    cout << fileList->GetEntries() << " raw-data files found" << endl;
    if ( fileList->GetEntries() == 0) {
      cout << "Exiting..." << endl;
      return;
    }

    // Take the first (or last?) file...
    TFileInfo *fi =  (TFileInfo *)fileList->At(0); 
    //  TFileInfo *fi =  (TFileInfo *)fileList->At(fileList->GetEntries()-1); 

    cout << "Getting the file:" << fi->GetCurrentUrl()->GetUrl() << endl;
    fi->Dump();

    filename = fi->GetCurrentUrl()->GetUrl();
  }
  else {
    // In case of local raw-data file...
    filename = localFileName;
  }

  AliLog::Flush();

  /////////////////////////////////////////////////////////////////////////////////////////
  //
  // First version of the reconstruction
  // script for the FDR'08

  // Set the CDB storage location
  // AliLog::SetModuleDebugLevel("STEER",2);
  AliCDBManager * man = AliCDBManager::Instance();
  //  man->SetDefaultStorage("local://LocalCDB");
  man->SetDefaultStorage("alien://folder=/alice/data/2008/LHC08a/OCDB/");
  
  // Files that we can not read from alien...solved
  //  man->SetSpecificStorage("ITS/Calib/MapsAnodeSDD","local://$ALICE_ROOT");
  //  man->SetSpecificStorage("ITS/Calib/MapsTimeSDD","local://$ALICE_ROOT");
  //  man->SetSpecificStorage("TPC/Calib/ExB","local://$ALICE_ROOT");

  // Objects not found if using LHC07w database...solved
  //  man->SetSpecificStorage("ITS/Calib/MapsAnodeSDD","local:///afs/cern.ch/user/c/cheshkov/public/OCDB");
  // man->SetSpecificStorage("GRP/GRP/Data","local://$ALICE_ROOT");
  // man->SetSpecificStorage("ITS/Calib/DDLMapSDD","local://$ALICE_ROOT");
  // man->SetSpecificStorage("MUON/Calib/Mapping","local://$ALICE_ROOT");
  // man->SetSpecificStorage("MUON/Calib/DDLStore","local://$ALICE_ROOT");

  // ITS settings
  AliITSRecoParam * itsRecoParam = AliITSRecoParam::GetCosmicTestParam();
  itsRecoParam->SetClusterErrorsParam(2);
  itsRecoParam->SetFindV0s(kFALSE);
  itsRecoParam->SetAddVirtualClustersInDeadZone(kFALSE);
  itsRecoParam->SetUseAmplitudeInfo(kFALSE);
  // In case we want to switch off a layer
  //  itsRecoParam->SetLayerToSkip(<N>);
  itsRecoParam->SetLayerToSkip(4);
  itsRecoParam->SetLayerToSkip(5);
  itsRecoParam->SetLayerToSkip(2);
  itsRecoParam->SetLayerToSkip(3);
  AliITSReconstructor::SetRecoParam(itsRecoParam);

  // TPC settings
  AliLog::SetClassDebugLevel("AliTPCclustererMI",2);
  AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetCosmicTestParam(kTRUE);
  tpcRecoParam->SetTimeInterval(60,940);
  tpcRecoParam->Dump();
  AliTPCReconstructor::SetRecoParam(tpcRecoParam);
  AliTPCReconstructor::SetStreamLevel(1);

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
  //  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., 1);
  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 0., 10., 2);
  AliTracker::SetFieldMap(field,1);
  Double_t mostProbPt=0.35;
  AliExternalTrackParam::SetMostProbablePt(mostProbPt);

  // AliReconstruction settings
  AliReconstruction rec;
  rec.SetUniformFieldTracking(kFALSE);
  rec.SetWriteESDfriend(kTRUE);
  rec.SetWriteAlignmentData();
  rec.SetInput(filename.Data());
  rec.SetRunReconstruction("ALL");
  rec.SetUseTrackingErrorsForAlignment("ITS");

  // In case some detectors have to be switched off...
  //  rec.SetRunLocalReconstruction("ALL");
  //  rec.SetRunTracking("ALL");
  //  rec.SetFillESD("ALL");
  // Enable vertex finder - it is needed for cosmic track reco
  rec.SetRunVertexFinder(kTRUE);

  // To be enabled if some equipment IDs are not set correctly by DAQ
  //  rec.SetEquipmentIdMap("EquipmentIdMap.data");

  // Detector options if any
  rec.SetOption("ITS","cosmics,onlyITS");
  rec.SetOption("MUON","SAVEDIGITS");
  rec.SetOption("TPC","OldRCUFormat");
  rec.SetOption("PHOS","OldRCUFormat");

  // To be enabled when CTP readout starts
  rec.SetFillTriggerESD(kFALSE);

  // all events in one single file
  rec.SetNumberOfEventsPerFile(-1);

  // switch off cleanESD
  rec.SetCleanESD(kFALSE);

  // rec.SetEventRange(0,15);
  // AliLog::SetGlobalDebugLevel(2);

  rec.SetRunQA(kFALSE);
  AliLog::Flush();
  rec.Run();

  //cout << "-----------------------------------------------------------------" << endl;
  //cout << "-----------------------------------------------------------------" << endl;
  //cout << "--------- Reconstruction Completed. Start merging QAs -----------" << endl;
  //cout << "-----------------------------------------------------------------" << endl;
  //cout << "-----------------------------------------------------------------" << endl;
  //AliQADataMakerSteer qas;
  //qas.Merge();
}
