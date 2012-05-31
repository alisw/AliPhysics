/**
 * @file HLTMultiplicityCorrelationsTest.C
 * @brief Macro for testing HLT MultipicityCorrelations
 *
 * This macro uses the ESDMCEventPublisher to publish AliESDEvents.
 * 
 * Results are written to disk by a rootfile writer
 *
 * @author Jochen Thaeder <jochen@thaeder.de>
 * @ingroup alihlt_physics
 */


// ---------------------------------------------------------------------------- 

/** HLTMultiplicityCorrelationsTest test macro
 *  @param nEvents Number of events which should be processed
 */
void HLTMultiplicityCorrelationsTest(const Char_t* esdpath="./",
				     const Char_t *cdbURI="local://$ALICE_ROOT/OCDB", 
				     Int_t nEvents=1) {
  
  // -- Set the CDB storage location
  // ---------------------------------
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage(cdbURI);
  man->SetRun(120000);

  TString writerInput;
  TString analysisInput;

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

  /*
    gHLT.LoadComponentLibraries("libAliHLTMUON.so");  
    gHLT.LoadComponentLibraries("libAliHLTTPC.so");  
    gHLT.LoadComponentLibraries("libAliHLTTRD.so");  
  */

  gHLT.LoadComponentLibraries("libAliHLTGlobal.so");  

  // ----------------------------//
  // -                         - //
  // -- Parameters            -- //
  // -                         - //
  // ----------------------------//

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
  TString publisherArg(Form("-entrytype ESD -datapath %s", esdPath.Data()));
  
  AliHLTConfiguration ESDMCEventPublisher(publisherId.Data(), publisherId.Data(), NULL, publisherArg.Data() );
  
  if (!analysisInput.IsNull()) analysisInput+=" ";
  analysisInput += publisherId;
    
  // ----------------------------//
  // -                         - //
  // -- Processing Components -- //
  // -                         - //
  // ----------------------------//

  // ------------------------------------------
  // -- Analysis
  // ------------------------------------------
  TString analysisId("MultiplicityCorrelations");
  // ------------------------------------------

  TString analysisArg("");
  
  AliHLTConfiguration analysis(analysisId.Data(),analysisId.Data(),
			       analysisInput.Data(), analysisArg.Data()); 
  
  if (!writerInput.IsNull()) writerInput+=" ";
  writerInput += analysisId;

  // ----------------------------//
  // -                         - //
  // --    Sink Components    -- //
  // -                         - //
  // ----------------------------//
  
  TString writerArg(Form("-directory analysis -datafile analysis_%d_%s.root", nEvents, analysisId.Data() ));

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
