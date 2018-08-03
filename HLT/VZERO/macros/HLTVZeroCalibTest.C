void HLTVZeroCalibTest()
{
  // set up HLT system to enable configuration registration
  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();

  AliHLTConfiguration clusterStatOrig("VZEROCalib", "VZEROOnlineCalib", "VZERO-RECO", "");
  AliHLTConfiguration rootWriter("RootWriter", "ROOTFileWriter", "VZEROCalib", "-directory testDir -datafile test.root");
}
