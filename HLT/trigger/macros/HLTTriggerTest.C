/**
 * @file HLTTriggerTest.C
 * @brief Macro for testing HLT Triggers
 *
 * This macro uses the ESDMCEventPublisher to publish AliESDEvents.
 * 
 * Results are written to disk by a rootfile writer
 *
 * @author Jochen Thaeder <jochen@thaeder.de>
 * @ingroup alihlt_trigger
 */


// ---------------------------------------------------------------------------- 


/** HLTTriggerTest test macro
 *  @param nEvents Number of events which should be processed
 */
void HLTTriggerTest(Int_t nEvents=1, 
		    const Char_t* triggername="H-Barrel_pT_Single-V0001.001",
		    const Char_t* esdpath="./") {
  
  // -- Set the CDB storage location
  // ---------------------------------
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  man->SetRun(120000);

  TString writerInput;
  TString globalInput;
  TString triggerInput;

  // -- Switch Logging
  // -------------------
  AliLog::SetGlobalLogLevel( AliLog::kError );
  AliHLTLogging log;
  log.SwitchAliLog(0);

  // -- Initialize HLT
  // -------------------
  AliHLTSystem gHLT;
  gHLT.SetGlobalLoggingLevel(0x7D);
  
  gHLT.LoadComponentLibraries("libESD.so");  
  gHLT.LoadComponentLibraries("libSTEER.so");  
  gHLT.LoadComponentLibraries("libSTEERBase.so");  
  gHLT.LoadComponentLibraries("libAOD.so");  
  gHLT.LoadComponentLibraries("libANALYSIS.so");  
  gHLT.LoadComponentLibraries("libANALYSISalice.so");  

  gHLT.LoadComponentLibraries("libHLTbase.so");
  gHLT.LoadComponentLibraries("libAliHLTUtil.so");

  gHLT.LoadComponentLibraries("libAliHLTMUON.so");  
  gHLT.LoadComponentLibraries("libAliHLTTPC.so");  
  gHLT.LoadComponentLibraries("libAliHLTTRD.so");  

  gHLT.LoadComponentLibraries("libAliHLTTrigger.so");  

  // ----------------------------//
  // -                         - //
  // -- Parameters            -- //
  // -                         - //
  // ----------------------------//

  TString triggerName(triggername);
  TString esdPath(esdpath);

  // ----------------------------//
  // -                         - //
  // -- Publisher  Components -- //
  // -                         - //
  // ----------------------------//

  // ------------------------------------------
  // -- The ESDMCEventPublisher   
  // ------------------------------------------
  TString publisherId("ESDMCEventPublisher");
  // ------------------------------------------
  TString publisherArg(Form("-entrytype HLTESD -datapath %s", esdPath.Data()));
  
  AliHLTConfiguration ESDMCEventPublisher(publisherId.Data(), publisherId.Data(), NULL, publisherArg.Data() );
  
  if (!triggerInput.IsNull()) triggerInput+=" ";
  triggerInput += publisherId;
    
  // ----------------------------//
  // -                         - //
  // -- Processing Components -- //
  // -                         - //
  // ----------------------------//

  // ------------------------------------------
  // -- Trigger
  // ------------------------------------------
  TString triggerId(triggerName);
  // ------------------------------------------

  TString triggerArg(Form("-triggername %s", triggerName.Data()));
  
  AliHLTConfiguration trigger(triggerId.Data(), "BarrelMultiplicityTrigger", 
			      triggerInput.Data(), triggerArg.Data()); 
  
  if (!globalInput.IsNull()) globalInput+=" ";
  globalInput += triggerId;

  if (!writerInput.IsNull()) writerInput+=" ";
  writerInput += triggerId;

  // ------------------------------------------
  // -- Global Trigger 
  // ------------------------------------------
  TString globalId("GlobalTrigger");
  // ------------------------------------------

  TString triggerArg(Form("-triggername %s", triggerName.Data()));
  
  AliHLTConfiguration global(globalId.Data(), "HLTGlobalTrigger", globalInput.Data(), ""); 
  
  if (!writerInput.IsNull()) writerInput+=" ";
  writerInput += globalId;

  // ----------------------------//
  // -                         - //
  // --    Sink Components    -- //
  // -                         - //
  // ----------------------------//
  
  TString writerArg(Form("-directory analysis -datafile trigger_%d_%s.root", nEvents, triggerName.Data() ));

  // -- The RootFileWriter 
  AliHLTConfiguration rootWriter("RootWriter", "ROOTFileWriter", writerInput.Data(), writerArg.Data() );
  
  // --------------------------- //
  // -                         - //
  // --         Run           -- //
  // -                         - //
  // --------------------------- //

  gHLT.BuildTaskList("RootWriter");
  gHLT.Run(nEvents);

  return;
}
