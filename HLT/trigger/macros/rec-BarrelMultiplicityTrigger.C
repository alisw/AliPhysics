void rec_BarrelMultiplicityTrigger(const char *filename="../raw.root")
{
  if(!gSystem->AccessPathName("galice.root")){
    cerr << "You have to run in a subfolder of the original data directory. If so, please delete the galice.root." << endl;
    return;
  }

  if (!filename) {
    cerr << "please specify input or run without arguments" << endl;
    return;
  }
  // Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage("local:///opt/HLT-public/OCDB/LHC09c");

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // init the HLT system in order to define the analysis chain below
  //
  AliHLTSystem* gHLT=AliHLTPluginBase::GetInstance();
 
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // define the analysis chain to be run
  //
  AliHLTConfiguration pubconf("hltesd-publisher", "ESDMCEventPublisher", NULL , "-entrytype HLTESD -datapath ..");

  AliHLTConfiguration triggerconf("multiplicity-trigger", "BarrelMultiplicityTrigger", "hltesd-publisher" , "-max-tdca 60");
  AliHLTConfiguration globaltriggerconf("global-trigger", "HLTGlobalTrigger", "multiplicity-trigger" , "");

  // Reconstruction settings
  AliReconstruction rec;

  // QA options
  rec.SetRunQA(":") ;

  // AliReconstruction settings
  rec.SetInput(filename);
  //rec.SetEventRange(0,20);
  rec.SetRunReconstruction("HLT TPC");
  rec.SetOption("HLT", "loglevel=0x7c chains=hltesd-publisher,global-trigger ignore-hltout");

  rec.Run();

}
