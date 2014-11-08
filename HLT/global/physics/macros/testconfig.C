void testconfig(const char* parent = "GLOBAL-esd-converter")
{
  // set up HLT system to enable configuration registration
  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();

  //pHLT->LoadComponentLibraries("libANALYSIS.so");  
  //pHLT->LoadComponentLibraries("libANALYSISalice.so");  
  /*
  pHLT->LoadComponentLibraries("libESD.so");  
  pHLT->LoadComponentLibraries("libSTEER.so");  
  pHLT->LoadComponentLibraries("libSTEERBase.so");  
  pHLT->LoadComponentLibraries("libAOD.so");  
  pHLT->LoadComponentLibraries("libANALYSIS.so");  
  pHLT->LoadComponentLibraries("libANALYSISalice.so");  

  pHLT->LoadComponentLibraries("libHLTbase.so");
  pHLT->LoadComponentLibraries("libAliHLTUtil.so");
  pHLT->LoadComponentLibraries("libAliHLTGlobal.so");  
  */

  /*
    pHLT->LoadComponentLibraries("libAliHLTMUON.so");  
    pHLT->LoadComponentLibraries("libAliHLTTPC.so");  
    pHLT->LoadComponentLibraries("libAliHLTTRD.so");  
  */


  // writer configuration
  AliHLTConfiguration calib("test" , "AnaManagerComponent" , parent , "");

  TString writerArg(Form("-directory testDir -datafile test.root"));

  // -- The RootFileWriter 
  AliHLTConfiguration rootWriter("RootWriter", "ROOTFileWriter", "test", writerArg.Data() );

  //pHLT->BuildTaskList("RootWriter");
  //pHLT->PrintTaskList();
}
