void calibTPCProcessor_config(const char* parent = "TPC-globalmerger TPC-ClusterTransformation")
{
  // set up HLT system to enable configuration registration
  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();

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
  AliHLTConfiguration calib("TPCtest" , "TPCCalibProcessor" , parent , "");

  TString writerArg(Form("-directory analysis -datafile analysis_test_calib.root"));

  // -- The RootFileWriter 
  AliHLTConfiguration rootWriter("RootWriter", "ROOTFileWriter", "TPCtest", writerArg.Data() );

  //pHLT->BuildTaskList("RootWriter");
  //pHLT->PrintTaskList();
}
