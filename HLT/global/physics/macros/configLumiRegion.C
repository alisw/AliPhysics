void configLumiRegion(const char* parent = "GLOBAL-esd-converter")
{
  // set up HLT system to enable configuration registration
  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();

  // writer configuration
  AliHLTConfiguration lumi("lumiRegion" , "LumiRegComponent" , parent , "-pushback-period 2");

  TString writerArg(Form("-directory lumiRegionFiles -datafile test_lumiRegion.root"));

  // -- The RootFileWriter 
  AliHLTConfiguration rootWriter("RootWriter", "ROOTFileWriter", "lumiRegion", writerArg.Data() );

}
