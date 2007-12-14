void recTPC2007(Int_t type, const char *filename="data.root")
{

  // .x recTPC(0,"rfio:/castor/cern.ch/alice/tpc/2007/12/13/00/07000011524011.980.root")
  //
  // Set path to calibration data
  //
  // type variable = 0 - cosmic test
  //               = 1 - laser test   
  //
  // Set reconstruction parameters
  ////
  gSystem->Load("/afs/cern.ch/alice/tpctest/root/HEAD/lib/libRAliEn.so");  
  TGrid::Connect("alien://");
  //

  AliLog::SetClassDebugLevel("AliTPCclustererMI",2);
  AliTPCRecoParam * tpcRecoParam = 0;
  if (type==0)  tpcRecoParam = AliTPCRecoParam::GetCosmicTestParam(kTRUE);
  if (type>0)  tpcRecoParam = AliTPCRecoParam::GetLaserTestParam(kTRUE);
  tpcRecoParam->SetTimeInterval(60,1000);
  tpcRecoParam->Dump();
  AliTPCReconstructor::SetRecoParam(tpcRecoParam);
  AliTPCReconstructor::SetStreamLevel(1);
  AliMagFMaps* field = new AliMagFMaps("Maps","Maps", 2, 1., 10., 1);
  AliTracker::SetFieldMap(field,1);
  //
  //
  AliReconstruction rec;  
  rec.SetDefaultStorage("alien://folder=/alice/data/2007/LHC07w/OCDB/");
  rec.SetWriteESDfriend(kTRUE);
  rec.SetInput(filename);
  rec.SetRunReconstruction("TPC");
  rec.SetEquipmentIdMap("EquipmentIdMap.data");
  //
  rec.SetOption("TPC","PedestalSubtraction OldRCUFormat");
  rec.SetFillESD("TPC");
  rec.SetFillTriggerESD(kFALSE);
  rec.SetRunVertexFinder(kFALSE);
  rec.Run();
}

  
