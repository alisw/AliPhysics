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
  // New handling of the specification: Specification will be the DDL ID
  // There will be 1 publisher per ddl, but the raw analyzers receive input from ALL DDLs
  
  Int_t ddlOffset = 4608; // The DDL offset for EMCAL (for PHOS the number is 1792)
  
  TString arg, fwInput, ecInput, clInput, rps, test;

  TString tmInput;   // Input for trigger maker
  TString tdInput;   // Input for trigger data maker
  
  for(int module = 0; module <= AliDAQ::NumberOfDdls("EMCAL"); module++) {
    TString publisher;
    // raw data publisher components
    publisher.Form("EMCAL-RP_%02d", module);
    arg.Form("-verbose -minid %d -datatype 'DDL_RAW ' 'EMCA'  -dataspec %d ", ddlOffset + module, module);
                
    if(rps.Length()) rps += " ";
    rps += publisher;
    AliHLTConfiguration pubConf(publisher.Data(), "AliRawReaderPublisher", NULL , arg.Data());
  }
                
  // Raw analyzer
  arg = "";
  AliHLTConfiguration rawConf("EMCAL-RA", "EmcalRawCrude", rps.Data(), arg.Data());
  fwInput+="EMCAL-RA";
                
  // Raw analyzer for TRU data
  arg = "";
  AliHLTConfiguration truConf("EMCAL-TRU", "EmcalTruAnalyzer", rps.Data(), arg.Data());
  if(tdInput.Length() > 0) tdInput += " ";
  tdInput += "EMCAL-TRU";
  if(fwInput.Length() > 0) fwInput += " ";
  fwInput+="EMCAL-TRU";

  // STU raw analyser
  AliHLTConfiguration stuConf("EMCAL-STU", "EmcalStuAnalyzer", rps.Data(), "");
  tdInput += " EMCAL-STU";
  if(fwInput.Length() > 0) fwInput += " ";
  fwInput+="EMCAL-STU";

  // digit maker components
  arg="";
  arg.Form("-sethighgainfactor 0.0153 -setlowgainfactor 0.2448 -setdigitthresholds 0.005 0.002");
  AliHLTConfiguration dmConf("EMCAL-DM", "EmcalDigitMaker", "EMCAL-RA", arg.Data());
  if(tmInput.Length() > 0) tmInput += " ";
  tmInput += " EMCAL-DM";
  if(clInput.Length() > 0) clInput += " ";
  clInput += " EMCAL-DM";
            
  arg = "";
  arg.Form("-digitthreshold 0.005 -recpointthreshold 0.1 -modulemode");
  AliHLTConfiguration clConf("EMCAL-CF", "EmcalClusterizer", clInput.Data(), arg.Data());
  ecInput += " EMCAL-CF";
            
  // Tigger data merger
  AliHLTConfiguration trgdata("EMCAL-TRG", "EmcalTriggerDataMaker", tdInput.Data(), "");
  tmInput += " EMCAL-TRG";
  fwInput += " EMCAL-TRG";
  ecInput += " EMCAL-TRG";

  AliHLTConfiguration trgmaker("EMCAL-TM", "EmcalTriggerMaker", tmInput.Data(), "");
  fwInput += " EMCAL-TM";

  
  // The call the histo maker 
  arg.Form("");
  arg.Form("-pushfraction 5 -beverbose 1");
  test = ecInput + " " + fwInput + " ";
  AliHLTConfiguration hfConf("emcalHisto","EmcalRawHistoMaker",test.Data(), arg.Data());
  ecInput += " emcalHisto";
  cout << ">>>>>>>>>> test:" << test << endl;	
 
  // SINK!
  AliHLTConfiguration esdcconf("ESD-CONVERTER", "GlobalEsdConverter"   , ecInput.Data(), "");
  
  
  // The filewriter 
  // SINK!
  arg.Form("-datatype 'CHANNELT' 'EMCA' -datafile RAout.root -concatenate-blocks -concatenate-events");
  AliHLTConfiguration fwConf("filewriter", "FileWriter", fwInput.Data(), arg.Data());
 
  // The call the histo maker 
  //arg.Form("");
  //arg.Form("-pushfraction 5 -beverbose 1");
  //AliHLTConfiguration hfConf("emcalHisto","EmcalRawHistoMaker",fwInput.Data(), arg.Data());
  
  // Write the root file 
  AliHLTConfiguration rwConf("rootFileHisto","ROOTFileWriter", "emcalHisto ESD-CONVERTER","-datafile roothisto.root -concatenate-events -overwrite");
  
  //////////////////////////////////////////////////////////////////////
  //
  // Init and run the reconstruction
  // All but HLT reconstruction is switched off 
  //
  /////////////////////////////////////////////////////////////////////
  
  TString option="libAliHLTUtil.so libAliHLTRCU.so libAliHLTCalo.so libAliHLTEMCAL.so libAliHLTGlobal.so chains=";
  //option+="ESD-CONVERTER";
  option+="rootFileHisto loglevel=0x5f";
  
  AliReconstruction rec;
  
  
  // uncomment for simulation
  // rec.SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  // rec.SetSpecificStorage("GRP/GRP/Data", Form("local://%s",gSystem->pwd()));
  //rec.SetDefaultStorage("local:///cvmfs/alice-ocdb.cern.ch/calibration/data/2013/OCDB");
  rec.SetDefaultStorage("raw://");
  
  rec.SetRunReconstruction(":");
  rec.SetEventRange(0,100);
  rec.SetInput(input);
  
  rec.SetInput("raw.root");
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunMultFinder(kFALSE);
  rec.SetRunVertexFinderTracks(kFALSE);
  rec.SetRunV0Finder(kFALSE);
  rec.SetRunCascadeFinder(kFALSE);
  
  rec.SetRunReconstruction("HLT");
  //rec.SetRunTracking(":");
  rec.SetLoadAlignFromCDB(0);
  rec.SetRunQA(":");
  rec.SetRunGlobalQA(kFALSE);
  
  rec.SetOption("HLT", option);
  
  rec.Run();
  
  
}
