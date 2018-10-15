void HLTVZeroCalibTest()
{
  // set up HLT system to enable configuration registration
  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();

  AliHLTConfiguration vzeroCalib("VZEROCalib", "VZEROOnlineCalib", "VZERO-RECO ITS-SPD-vertexer", "");
  AliHLTConfiguration rootWriter("RootWriter", "ROOTFileWriter", "VZEROCalib", "-directory testDir -datafile test.root");
}
