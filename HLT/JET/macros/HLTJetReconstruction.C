/**
 * @file HLTJetReconstruction.C
 * @brief Macro for testing HLT Jet Reconstruction
 *
 * This macro uses the ESDMCEventPublisher to publish AliMCEvents and AliESDEvents.
 * The processing is done be JetFinder's
 * Results are written to disk by a rootfile writer
 *
 * @author thaeder@kip.uni-heidelberg.de
 * @ingroup alihlt_jet
 */

/** HLTJetReconstruction test macro
 *  @param nEvents Number of events which should be processed
 */
void HLTJetReconstruction(Int_t nEvents=1) {

  TString writerInput;
  TString arg;

  // this is just a tool to switch the logging systems
  AliHLTLogging log;
  log.SwitchAliLog(1);

  AliHLTSystem gHLT;
  gHLT.SetGlobalLoggingLevel(0x7F);
  gHLT.LoadComponentLibraries("libAliHLTUtil.so");
  gHLT.LoadComponentLibraries("libAliHLTJET.so");  

  // -                         - //
  // -- Publisher  Components -- //
  // -                         - //

  arg.Form("-entrytype MCFAST -dataspec 0x0000001F -datapath /home/jthaeder/jet/data/v4-16-Rev-01/FastGen/kPythia6Jets104_125_14TeV/JET-ETA=-0.2,0.2_JET-ET=50,1000_R=0.4_500ev");

  // -- The ESDMCEventPublisher 
  AliHLTConfiguration ESDMCEventPublisher("ESDMCEventPublisher", "ESDMCEventPublisher", 
					  NULL, arg.Data() );
  
  if (!writerInput.IsNull()) writerInput+=" ";
  writerInput+="ESDMCEventPublisher";
  
  // -                         - //
  // -- Processing Components -- //
  // -                         - //
  AliHLTConfiguration jetFinder("JETConeJet", "JETConeJetFinder",
				"ESDMCEventPublisher","");

  if (!writerInput.IsNull()) writerInput+=" ";
  writerInput+="JETConeJet";

  // -                         - //
  // --    Sink Components    -- //
  // -                         - //

  // -- The RootFileWriter 
  AliHLTConfiguration rootWriter("RootWriter", "ROOTFileWriter", 
				 writerInput.Data(), "-datafile event");


  // -                         - //
  // --         Run           -- //
  // -                         - //

  gHLT.BuildTaskList("RootWriter");
  gHLT.Run(nEvents);

}
