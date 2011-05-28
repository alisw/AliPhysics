void hlt_rec(const char* input="./") {

  // For real data:
 // AliCDBManager::Instance()->SetDefaultStorage("raw://");
 
  gStyle->SetPalette(1);

  if(!gSystem->AccessPathName("galice.root")){
    cerr << "please delete the galice.root or run at different place." << endl;
    return;
  }
 
  if (!input) {
    cerr << "please specify input or run without arguments" << endl;
    return;
  }

  //////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // init the HLT system in order to define the analysis chain below
  //
  
  AliHLTSystem* gHLT=AliHLTPluginBase::GetInstance();

  gHLT->LoadComponentLibraries("libAliHLTUtil.so");
  gHLT->LoadComponentLibraries("libAliHLTRCU.so");
  gHLT->LoadComponentLibraries("libAliHLTCalo.so");
  gHLT->LoadComponentLibraries("libAliHLTEMCAL.so");
  gHLT->LoadComponentLibraries("libAliHLTPHOS.so");
  gHLT->LoadComponentLibraries("libAliHLTGlobal.so");

  //////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // define the analysis chain to be run
  //

  // Define which modules and RCUs to use
  int moduleStart = 0;
  int moduleEnd = 9;
  int rcuStart = 0;
  int rcuEnd = 1;
  
  // The number of RCUs per module in EMCAL (used for calculating data specification)

  Int_t rcusPerModule = 2;
  Int_t ddlOffset = 4608; // The DDL offset for EMCAL (for PHOS the number is 1792)
  
  TString arg, fwInput, ecInput, dm;

  // Loop over the modules/RCUs
  for (int module = moduleStart; module <= moduleEnd; module++) 
    {
      TString clInput;
      
      for(int rcu = rcuStart; rcu <= rcuEnd; rcu++) 
	{
	  TString publisher, rawanalyzer, dm;

	  // Raw data publisher
	  publisher.Form("EMC-RP_%02d_%d", module, rcu);
	  
	 arg.Form("-verbose -minid %d -datatype 'DDL_RAW ' 'EMCA'  -dataspec 0x%x ", ddlOffset + module*(rcusPerModule) + rcu, 0x1 << (module*rcusPerModule + rcu));
	  AliHLTConfiguration pubConf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());
	  cout << arg.Data() << endl;

	  // Raw analyzer
	  arg = "";
	  rawanalyzer.Form("EMC-RA_%02d_%d", module, rcu);
	  AliHLTConfiguration rawConf(rawanalyzer.Data(), "EmcalRawCrude", publisher.Data(), arg.Data());
	   
	  // If you want to have a component subscribing to all the raw analyzers this will be it's input
	  if(fwInput.Length() > 0) fwInput += " ";
	  fwInput+=rawanalyzer;
	  cout << ">>>>>>>>>> fwInput:" << fwInput << endl;
	  
          // digit maker components
          dm.Form("EMC-DM_%02d_%d", module, rcu);
          arg="";
	  // HI-GAIN ~ 15/256 ------ CHECK LOW GAIN FACTOR
          arg.Form("-highgainfactor 0.0167 -lowgainfactor 0.08");
          AliHLTConfiguration dmConf(dm.Data(), "EmcalDigitMaker", rawanalyzer.Data(), arg.Data());
          if(clInput.Length() > 0) clInput += " ";
          clInput+=dm;
	  cout << ">>>>>>>>>> clInput:" << clInput << endl;	
	
	} // end of RCU loop

   	TString cl, ca;

      	cl.Form("EMC-CL_%02d", module);
      	arg = "";
	// digitthreshold was 0.005
      	arg.Form("-digitthreshold 0.0167 -recpointthreshold 0.1");
      	AliHLTConfiguration clConf(cl.Data(), "EmcalClusterizer", clInput.Data(), arg.Data());
      	if(ecInput.Length() > 0) ecInput += " ";
      	ecInput += cl;

    } // end of module loop

  TString ec, test;

     // The call the histo maker 
  	arg.Form("");
  	arg.Form("-pushfraction 5 -beverbose 1");
	test = ecInput + " " + fwInput + " ";
  	AliHLTConfiguration hfConf("emcalHisto","EmcalRawHistoMaker",test.Data(), arg.Data());
	ecInput += " emcalHisto";
	cout << ">>>>>>>>>> test:" << test << endl;	
		
    	ec.Form("ESD-CONVERTER");
    	arg = "";
 	AliHLTConfiguration esdcconf(ec.Data(), "GlobalEsdConverter"   , ecInput.Data(), "");


  // The filewriter 
//  arg.Form("-datatype 'CHANNELT' 'EMCA' -datafile RAout.root -concatenate-blocks -concatenate-events");
//  AliHLTConfiguration fwConf("filewriter", "FileWriter", fwInput.Data(), arg.Data());
  
  // The call the histo maker 
//  arg.Form("");
//  arg.Form("-pushfraction 5 -beverbose 1");
//  AliHLTConfiguration hfConf("emcalHisto","EmcalRawHistoMaker",fwInput.Data(), arg.Data());
  
  // Write the root file 
  //AliHLTConfiguration rwConf("rootFileHisto","ROOTFileWriter", "emcalHisto ESD-CONVERTER","-datafile roothisto.root -concatenate-events -overwrite");



  
  //////////////////////////////////////////////////////////////////////
  //
  // Init and run the reconstruction
  // All but HLT reconstruction is switched off 
  //
  /////////////////////////////////////////////////////////////////////

  TString option="libAliHLTUtil.so libAliHLTRCU.so libAliHLTCalo.so libAliHLTEMCAL.so libAliHLTGlobal.so chains=";
  option+="ESD-CONVERTER";
  //option+="rootFileHisto loglevel=0x5f";

  AliReconstruction rec;
  

  // uncomment for simulation
  rec.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  rec.SetSpecificStorage("GRP/GRP/Data", Form("local://%s",gSystem->pwd()));

  rec.SetRunReconstruction(":");
  rec.SetEventRange(0,10);
  rec.SetInput(input);
  
  //rec.SetInput("raw.root");
  rec.SetRunVertexFinder(kFALSE);

  rec.SetRunReconstruction("HLT");
  //rec.SetRunTracking(":");
  rec.SetLoadAlignFromCDB(0);
  rec.SetRunQA(":");

  rec.SetOption("HLT", option);
  
  rec.Run();


}
