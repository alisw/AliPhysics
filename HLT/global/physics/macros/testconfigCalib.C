void testconfigCalib(const char* parent = "GLOBAL-flat-esd-converter", const char* tpcConfig="AddTaskMacro=$ALICE_PHYSICS/PWGPP/CalibMacros/CPass0/AddTaskTPCCalib.C(\"TPCCalib:CalibTimeDrift\") WriteAnalysisToFile=0 -pushback-period=1000")
{
  // set up HLT system to enable configuration registration
  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();

  AliHLTConfiguration calibTPC("calib" , "HLTAnalysisManagerComponent" , parent , tpcConfig);
  //AliHLTConfiguration vzeroBG("calib" , "HLTAnalysisManagerComponent" , parent , "AddTaskMacro=$ALICE_PHYSICS/PWGPP/BeamGasMonitoring/macros/AddTaskBGMonitorQA.C AddTaskArgs= -pushback-period=1000");

  AliHLTConfiguration rootWriter("RootWriterTPCcalib", "ROOTFileWriter", "calib", "-directory hltoutput -datafile hltoutput.root");

}
