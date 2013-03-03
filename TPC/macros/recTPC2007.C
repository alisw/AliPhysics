void recTPC2007(Int_t type, const char *filename="data.root")
{
  //  .L $ALICE_ROOT/TPC/macros/recTPC2007.C
  
  // recTPC2007(0,"rfio:/castor/cern.ch/alice/data/2007/12/16/11/07000013151011.20.root")
  //recTPC2007(0,"alien:///alice/data/2007/LHC07w_TPC/000012674/raw/07000012674011.90.root")
  //  recTPC2007(0,"/data/test2007/07000012674011/07000012674011.90.root")

  //
  // Set path to calibration data
  //
  // type variable = 0 - cosmic test
  //               = 1 - laser test   
  //
  // Set reconstruction parameters
  ////
  //gSystem->Load("/afs/cern.ch/alice/tpctest/root/HEAD/lib/libRAliEn.so");  
  //gSystem->Load("$ROOTSYS/lib/libXrdClient.so") 
  TGrid::Connect("alien://");
  //

  AliLog::SetClassDebugLevel("AliTPCclusterer",2);
  AliTPCRecoParam * tpcRecoParam = 0;
  if (type==0)  tpcRecoParam = AliTPCRecoParam::GetCosmicTestParam(kTRUE);
  if (type>0)  tpcRecoParam = AliTPCRecoParam::GetLaserTestParam(kTRUE);
  tpcRecoParam->SetTimeInterval(60,940);
  tpcRecoParam->Dump();
  //
  //
 //  tpcRecoParam->SetMinMaxCutAbs(4.);
//   tpcRecoParam->SetMinLeftRightCutAbs(6.);
//   tpcRecoParam->SetMinUpDownCutAbs(7.);
//   tpcRecoParam->SetMinMaxCutSigma(4.);
//   tpcRecoParam->SetMinLeftRightCutSigma(6.);
//   tpcRecoParam->SetMinUpDownCutSigma(7.);
  //
  //
  AliTPCReconstructor::SetRecoParam(tpcRecoParam);
  AliTPCReconstructor::SetStreamLevel(100);
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
  rec.SetCleanESD(kFALSE);
  rec.Run();
}

  
