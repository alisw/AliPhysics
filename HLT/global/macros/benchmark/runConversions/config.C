void config( TString parent="GLOBAL-flat-esd-converter", TString directory="outFlat", TString fileName="outFlatHLT.dat")
{

  cout<<"Now entering config.C"<<endl;
  // set up HLT system to enable configuration registration
  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();



  // writer configuration
  

	// -- The RootFileWriter 
	if(fileName.EndsWith(".root")){
  cout<<"use RootFileWriter"<<endl;
		AliHLTConfiguration RootWriter("RootWriter", "ROOTFileWriter", parent, "-directory " + directory + " -datafile " + fileName + " -concatenate-blocks" + " -concatenate-events");
	
		
	}
	else{
  cout<<"use normal file writer"<<endl;
		 AliHLTConfiguration RootWriter("RootWriter", "FileWriter", parent, "-directory " + directory + " -datafile " + fileName + " -concatenate-blocks" + " -concatenate-events");
	 }
  //pHLT->BuildTaskList("RootWriter");
  //pHLT->PrintTaskList();
}
