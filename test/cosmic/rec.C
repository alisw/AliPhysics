void rec(const char *filename="raw.root", const Int_t mfield=1)
{
  /////////////////////////////////////////////////////////////////////////////////////////
  //
  // Second version of the reconstruction
  // script for the 2008 cosmic data (LHC08b) 
  //
  /////////////////////////////////////////////////////////////////////////////////////////
  //AliLog::SetGlobalLogLevel(AliLog::kWarning);
  AliLog::SetGlobalLogLevel(AliLog::kError);

  gSystem->Load("libRAliEn.so");
  gSystem->Load("libNet.so");
  gSystem->Load("libMonaLisa.so");
  new TMonaLisaWriter(0, "GridAliRoot-rec.C", 0, 0, "global");
  gSystem->Setenv("APMON_INTERVAL", "120");

  // Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();
    man->SetDefaultStorage("local://$ALICE_ROOT");
  //man->SetDefaultStorage("alien://folder=/alice/data/2008/LHC08a/OCDB/");
  man->SetSpecificStorage("ITS/Calib/*","local://$ALICE_ROOT");
  
  // Example in case a specific CDB storage is needed
  //  man->SetSpecificStorage("ITS/Calib/MapsAnodeSDD","local://$ALICE_ROOT");

  // Reconstruction settings
  AliReconstruction rec;

  // ITS settings
  AliITSRecoParam * itsRecoParam = AliITSRecoParam::GetCosmicTestParam();
  itsRecoParam->SetFactorSAWindowSizes(20);
  itsRecoParam->SetClusterErrorsParam(2);
  itsRecoParam->SetFindV0s(kFALSE);
  itsRecoParam->SetAddVirtualClustersInDeadZone(kFALSE);
  itsRecoParam->SetUseAmplitudeInfo(kFALSE);
  // In case we want to switch off a layer
  //  itsRecoParam->SetLayerToSkip(<N>);
  //  itsRecoParam->SetLayerToSkip(4);
  //  itsRecoParam->SetLayerToSkip(5);
  itsRecoParam->SetLayerToSkip(2);
  itsRecoParam->SetLayerToSkip(3);
  //itsRecoParam->SetSAOnePointTracks();
  itsRecoParam->SetClusterMisalError(0.1); // [cm]
  itsRecoParam->SetSAUseAllClusters();
  rec.SetRecoParam("ITS",itsRecoParam);

  // TPC settings
  //AliLog::SetClassDebugLevel("AliTPCclustererMI",2);
  AliTPCRecoParam * tpcRecoParam = AliTPCRecoParam::GetCosmicTestParam(kFALSE);
  tpcRecoParam->SetTimeInterval(60,940);
  Double_t sysError[5]={0.3,1, 0.3/150., 1./150.,0.3/(150*150.)};
  tpcRecoParam->SetSystematicError(sysError);
  tpcRecoParam->SetMinMaxCutAbs(4.);
  tpcRecoParam->SetMinLeftRightCutAbs(6.);
  tpcRecoParam->SetMinUpDownCutAbs(6.);
  //  tpcRecoParam->Dump();
  rec.SetRecoParam("TPC",tpcRecoParam);
  AliTPCReconstructor::SetStreamLevel(1);

  // TRD setting
  // Settings for the TRD Raw Reader
  AliTRDrawStreamBase::SetRawStreamVersion("TB");
  AliTRDrawStreamTB::SetNoErrorWarning();
  AliTRDrawStreamTB::AllowCorruptedData();
  AliTRDrawStreamTB::DisableStackNumberChecker();
  AliTRDrawStreamTB::DisableStackLinkNumberChecker();
  AliTRDrawStreamTB::SetSubtractBaseline(10);
  
  // TRD reconstruction params
  AliTRDrecoParam *fTRDrecoParam = AliTRDrecoParam::GetCosmicTestParam();
  rec.SetRecoParam("TRD",fTRDrecoParam);
  AliTRDtrackerV1::SetNTimeBins(30);

  // PHOS settings
  AliPHOSRecoParam* recPHOS = new AliPHOSRecoParam();
  recPHOS->SetEMCSubtractPedestals(kTRUE);
  recPHOS->SetEMCMinE(0.05);
  recPHOS->SetEMCClusteringThreshold(0.10);
  rec.SetRecoParam("PHOS",recPHOS);

  // T0 settings
  //AliLog::SetModuleDebugLevel("T0", 10);

  // MUON settings
  //AliLog::SetClassDebugLevel("AliMUONRawStreamTracker",3);
  AliMUONRecoParam *muonRecoParam = AliMUONRecoParam::GetLowFluxParam();
  muonRecoParam->CombineClusterTrackReco(kTRUE);
  muonRecoParam->SetCalibrationMode("NOGAIN");
  //muonRecoParam->SetClusteringMode("PEAKFIT");
  //muonRecoParam->SetClusteringMode("PEAKCOG");
  muonRecoParam->Print("FULL");
  rec.SetRecoParam("MUON",muonRecoParam);
 
  // Tracking settings
  AliMagWrapCheb* field;
  if (mfield)
     field = new AliMagWrapCheb("Maps","Maps", 2, 1., 10., AliMagWrapCheb::k5kG);
   else
     field = new AliMagWrapCheb("Maps","Maps", 2, 0., 10., AliMagWrapCheb::k2kG);
  AliTracker::SetFieldMap(field,1);
  Double_t mostProbPt=0.35;
  AliExternalTrackParam::SetMostProbablePt(mostProbPt);

  // AliReconstruction settings
  rec.SetUniformFieldTracking(kFALSE);
  rec.SetWriteESDfriend(kTRUE);
  rec.SetWriteAlignmentData();
  rec.SetInput(filename);
  //  rec.SetRunReconstruction("ALL");
  rec.SetUseTrackingErrorsForAlignment("ITS");

  // In case some detectors have to be switched off...
  rec.SetRunReconstruction("ITS TPC TRD TOF HMPID PHOS MUON FMD PMD T0 VZERO ZDC ACORDE");

  // Enable vertex finder - it is needed for cosmic track reco
  rec.SetRunVertexFinder(kTRUE);

  // To be enabled if some equipment IDs are not set correctly by DAQ
  //  rec.SetEquipmentIdMap("EquipmentIdMap.data");

  // Detector options if any
  rec.SetOption("ITS","cosmics");
  rec.SetOption("MUON","SAVEDIGITS");
  rec.SetOption("T0","cosmic");

  // Enabled when CTP readout starts
  rec.SetFillTriggerESD(kTRUE);

  // all events in one single file
  rec.SetNumberOfEventsPerFile(-1);

  // switch off cleanESD
  rec.SetCleanESD(kFALSE);

  //rec.SetEventRange(0,15);

  rec.SetRunQA("ITS TPC:ESD RECPOINT");
  rec.SetRunGlobalQA(kTRUE);
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
