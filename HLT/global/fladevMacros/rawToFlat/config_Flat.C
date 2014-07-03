void config_Flat( TString directory="outFlat", TString fileName="outFlatHLT.dat")
{

  cout<<"Now entering config_Flat"<<endl;
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
  

  // -- The RootFileWriter 
	  AliHLTConfiguration RootWriter("RootWriter", "FileWriter", "GLOBAL-flat-esd-converter", "-directory " + directory + " -datafile " + fileName );

  //pHLT->BuildTaskList("RootWriter");
  //pHLT->PrintTaskList();
}
