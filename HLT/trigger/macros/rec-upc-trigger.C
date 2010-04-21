// $Id$
/**
 * @file rec-upc-trigger.C
 * @brief Test macro for UPC trigger
 *
 * Usage:
 * <pre>
 *   aliroot -b -q rec-upc-trigger.C | tee rec-upc-trigger.log
 * </pre>
 *
 * The macro asumes raw data to be available in the rawx folders, either
 * simulated or real data. A different input can be specified as parameter
 * <pre>
 *   aliroot -b -q rec-upc-trigger.C'("input.root")'
 * </pre>
 *
 *
 * @author Kyrre Skjerdal (kyrre.skjerdal@cern.ch)
 */

void rec_upc_trigger(const char *filename="raw.root"){

  if(!gSystem->AccessPathName("galice.root")){
    cerr << "Remove galice.root or run in a different folder." << endl;
    return;
  }

  if (!filename) {
    cerr << "please specify input or run without arguments" << endl;
    return;
  }
  
  // Set the CDB storage location
  AliCDBManager *man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB/");
  

  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // init the HLT system in order to define the analysis chain below
  //
  AliHLTSystem *gHLT = AliHLTPluginBase::GetInstance();
 
 
 
  ///////////////////////////////////////////////////////////////////////////////////////////////////
  //
  // define the analysis chain to be run
  //
 
 // AliHLTConfiguration pubconf("hltesd-publisher", "ESDMCEventPublisher", NULL , "-entrytype HLTESD -datapath ..");

  AliHLTConfiguration triggerconf("upc", "UpcTrigger", "GLOBAL-esd-converter", "");
 //AliHLTConfiguration globaltriggerconf("global-trigger", "HLTGlobalTrigger", "multiplicity-trigger" , "");

  // Reconstruction settings
  AliReconstruction rec;

  // QA options
  rec.SetRunQA(":") ;

  // AliReconstruction settings
  rec.SetInput(filename);
  rec.SetSpecificStorage("GRP/GRP/Data",Form("local://%s",gSystem->pwd()));
 
  rec.SetEventRange(0,100);
 
  rec.SetRunVertexFinder(kFALSE);
  
  rec.SetRunLocalReconstruction("HLT");
  rec.SetRunTracking("");
  rec.SetFillESD("");
  
  //rec.SetRunReconstruction("HLT");
  //rec.SetLoadAlignFromCDB(0);
  
  rec.SetOption("HLT", "loglevel=0x7c chains=upc");

  rec.Run();

}
