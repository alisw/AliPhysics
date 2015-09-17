void testconfigCalib(const char* parent = "GLOBAL-esd-converter", const char* tpcConfig="AddTaskMacro=$ALICE_PHYSICS/PWGPP/CalibMacros/CPass0/AddTaskTPCCalib.C AddTaskArgs=\"TPCCalib:CalibTimeDrift\" WriteAnalysisToFile=1")
{
  //load the libraries needed by the calib code
  gROOT->Macro("$ALICE_PHYSICS/PWGPP/CalibMacros/CPass0/LoadLibraries.C");
  // set up HLT system to enable configuration registration
  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();

  //pHLT->LoadComponentLibraries("libANALYSIS.so");  
  //pHLT->LoadComponentLibraries("libANALYSISalice.so");  

  /*pHLT->LoadComponentLibraries("libESD.so");
  pHLT->LoadComponentLibraries("libSTEER.so");  
  pHLT->LoadComponentLibraries("libSTEERBase.so");  
  pHLT->LoadComponentLibraries("libAOD.so");  
  pHLT->LoadComponentLibraries("libANALYSIS.so");  
  pHLT->LoadComponentLibraries("libANALYSISalice.so");  

  pHLT->LoadComponentLibraries("libHLTbase.so");
  pHLT->LoadComponentLibraries("libAliHLTUtil.so");
  pHLT->LoadComponentLibraries("libAliHLTGlobal.so");*/

  /*pHLT->LoadComponentLibraries("libSTAT.so");
  pHLT->LoadComponentLibraries("libANALYSISalice.so");
  pHLT->LoadComponentLibraries("libANALYSIScalib.so");
  //
  // detector libraries
  //
  pHLT->LoadComponentLibraries("libTPCcalib.so");
  pHLT->LoadComponentLibraries("libTRDcalib.so");
  pHLT->LoadComponentLibraries("libT0calib.so");
  pHLT->LoadComponentLibraries("libTOFcalib.so");
  //
  // PWGPP libraries
  //
  pHLT->LoadComponentLibraries("libANALYSISalice.so");
  pHLT->LoadComponentLibraries("libANALYSIScalib.so");
  pHLT->LoadComponentLibraries("libTENDER.so");
  pHLT->LoadComponentLibraries("libPWGPP.so");*/

  /*
    pHLT->LoadComponentLibraries("libAliHLTMUON.so");  
    pHLT->LoadComponentLibraries("libAliHLTTPC.so");  
    pHLT->LoadComponentLibraries("libAliHLTTRD.so");  
  */


  // writer configuration
  AliHLTConfiguration calib("calib" , "HLTAnalysisManagerComponent" , parent , tpcConfig);

  TString writerArg(Form("-directory testDir -datafile test.root"));

  // -- The RootFileWriter 
  AliHLTConfiguration rootWriter("RootWriter", "ROOTFileWriter", "calib", writerArg.Data() );

  //pHLT->BuildTaskList("RootWriter");
  //pHLT->PrintTaskList();
}
