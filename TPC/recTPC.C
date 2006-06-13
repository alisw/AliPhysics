//
// example macro for reconstruction of the TPC raw data
//
// The path to the Calibration parameters is for the moment hard-wired in the code
// Taken from /afs/
//
//

void recTPC(const char *filename="../dataroot/run385.001.root")
{
  AliLog::SetClassDebugLevel("AliTPCclustererMI",2);
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT");
  man->SetRun(0);

  man->SetSpecificStorage("TPC","local:///afs/cern.ch/user/m/mivanov/public/Calib");

  AliReconstruction rec;  
  rec.SetSpecificStorage("TPC","local:///afs/cern.ch/user/m/mivanov/public/Calib");
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
  AliTPCReconstructor::SetCtgRange(100.);
  AliTPCReconstructor::SetStreamLevel(1);
  rec.SetWriteAlignmentData(kTRUE);
  rec.Run();
}

void recTracking(const char *filename="../run439.001.root", Int_t nevents=1)
{
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT");
  man->SetRun(0);

  man->SetSpecificStorage("TPC","local:///afs/cern.ch/user/m/mivanov/public/Calib");
  AliReconstruction rec;
  //rec.SetSpecificStorage("TPC","local:///afs/cern.ch/user/m/mivanov/public/Calib");
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
  AliTPCReconstructor::SetCtgRange(100.);
  AliTPCReconstructor::SetStreamLevel(1);
  rec.SetWriteAlignmentData(kTRUE);
  rec.Run(0,nevents);
}


