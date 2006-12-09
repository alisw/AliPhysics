//
// example macro for reconstruction of the TPC raw data
//
// The path to the Calibration parameters is for the moment hard-wired in the code
// Taken from /afs/
//
//

void recTPC(Int_t type, const char *filename="data.root")
{
  //
  // Set path to calibration data
  //
  // type variable = 0 - cosmic test
  //               = 1 - laser test   
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT");
  //man->SetRun(0);
  //man->SetSpecificStorage("TPC/*/*","local:///afs/cern.ch/user/m/mivanov/public/Calib");
  //
  // Set reconstruction parameters
  //
  AliLog::SetClassDebugLevel("AliTPCclustererMI",2);
  AliTPCRecoParam * tpcRecoParam = 0;
  if (type==0)  tpcRecoParam = AliTPCRecoParam::GetCosmicTestParam(kTRUE);
  if (type>0)  tpcRecoParam = AliTPCRecoParam::GetLaserTestParam(kTRUE);
  tpcRecoParam->Dump();
  AliTPCReconstructor::SetRecoParam(tpcRecoParam);
  AliTPCReconstructor::SetStreamLevel(1);
  //
  //
  //
  AliReconstruction rec;  
  rec.SetDefaultStorage("local://$ALICE_ROOT");
  rec.SetSpecificStorage("TPC/*/*","local:///afs/cern.ch/user/m/mivanov/public/Calib");
  rec.SetLoadAlignData("");
  rec.SetWriteESDfriend(kTRUE);
  rec.SetInput(filename);
  rec.SetEquipmentIdMap("EquipmentIdMap.data");
  rec.SetRunReconstruction("TPC");
  rec.SetOption("TPC","PedestalSubtraction OldRCUFormat");
  //  rec.SetRunLocalReconstruction("");
  //  rec.SetRunTracking("TPC");
  rec.SetFillESD("TPC");
  rec.SetFillTriggerESD(kFALSE);
  rec.SetRunVertexFinder(kFALSE);
  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., 1);
  AliTracker::SetFieldMap(field,1);
  rec.SetWriteAlignmentData(kTRUE);
  rec.Run();
}

void recTracking(Int_t type, const char *filename="data.root")
{
  //
  // Set path to calibration data
  //
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT");
  man->SetRun(0);
  man->SetSpecificStorage("TPC/*/*","local:///afs/cern.ch/user/m/mivanov/public/Calib");
  //
  // Set reconstruction parameters
  //
  AliLog::SetClassDebugLevel("AliTPCclustererMI",2);

  AliTPCRecoParam * tpcRecoParam = 0;
  if (type==0)  tpcRecoParam = AliTPCRecoParam::GetCosmicTestParam(kTRUE);
  if (type>0)  tpcRecoParam = AliTPCRecoParam::GetLaserTestParam(kTRUE);

  AliTPCReconstructor::SetRecoParam(tpcRecoParam);
  AliTPCReconstructor::SetStreamLevel(1);

  //
  //
  //
  AliReconstruction rec;
  rec.SetSpecificStorage("TPC/*/*","local:///afs/cern.ch/user/m/mivanov/public/Calib");
  rec.SetLoadAlignData("");
  rec.SetWriteESDfriend(kTRUE);
  rec.SetInput(filename);
  rec.SetEquipmentIdMap("EquipmentIdMap.data");
  //rec.SetRunReconstruction("TPC");
  rec.SetOption("TPC","PedestalSubtraction OldRCUFormat");
  rec.SetRunLocalReconstruction("");
  rec.SetRunTracking("TPC");
  rec.SetFillESD("TPC");
  rec.SetFillTriggerESD(kFALSE);
  rec.SetRunVertexFinder(kFALSE);
  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., 1);
  AliTracker::SetFieldMap(field,1);
  rec.SetWriteAlignmentData(kTRUE);
  rec.Run(0);
}


