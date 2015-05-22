void testconfigED(const char* parent = "GLOBAL-esd-converter", const char* config="")
{
  gROOT->Macro("$ALICE_PHYSICS/PWGPP/CalibMacros/CPass0/LoadLibraries.C");
  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();

  // writer configuration
  AliHLTConfiguration calib("ED" , "ZMQsink" , parent , config);

  // -- The RootFileWriter 
  //TString writerArg(Form("-directory testDir -datafile test.root"));
  //AliHLTConfiguration rootWriter("RootWriter", "ROOTFileWriter", "test", writerArg.Data() );
}
