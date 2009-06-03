/**
 * @file MCGeneratorPublisher.C
 * @brief Macro for testing of the generation and publishing of AliHLTMCEvent
 *
 * This macro is a testing/example macro of how to use the MCGenerator 
 * (AliHLTMCGeneratorComponent) and RootFileWriter (AliHLTRootFileWriter). 
 * It defines only two component in the chain, the publisher, which 
 * publishes the content of the root files according to the selection. 
 * Be aware there can be several root objects in one root file.
 *
 * For more descriptions, especially the used datatypes and specification:
 * @see AliHLTMCGeneratorComponent
 *
 * @author thaeder@kip.uni-heidelberg.de
 * @ingroup alihlt_tutorial
 */

// ---------------------------------------------------------------------------- 
// -- Pythia Parameter 
// ---------------------------------------------------------------------------- 

enum PprRun_t {
  kPythia6Jets20_24,   kPythia6Jets24_29,   kPythia6Jets29_35,
  kPythia6Jets35_42,   kPythia6Jets42_50,   kPythia6Jets50_60,
  kPythia6Jets60_72,   kPythia6Jets72_86,   kPythia6Jets86_104,
  kPythia6Jets104_125, kPythia6Jets125_150, kPythia6Jets150_180,
  kPyJetJet, kPyGammaJetPHOS, kRunMax
};

const Char_t* pprRunName[] = {
  "kPythia6Jets20_24",   "kPythia6Jets24_29",   "kPythia6Jets29_35",
  "kPythia6Jets35_42",   "kPythia6Jets42_50",   "kPythia6Jets50_60",
  "kPythia6Jets60_72",   "kPythia6Jets72_86",   "kPythia6Jets86_104",
  "kPythia6Jets104_125", "kPythia6Jets125_150", "kPythia6Jets150_180",
  "kPyJetJet", "kPyGammaJetPHOS"
};


/** ESDMCEventPublisher test macro
 *  @param nEvents Number of events which should be processed
 */
void MCGeneratorPublisher(Int_t nEvents=10, PprRun_t runType = kPythia6Jets104_125) {

  TString writerInput;

  // -- Switch Logging
  // -------------------
  AliLog::SetGlobalLogLevel( AliLog::kError );
  AliHLTLogging log;
  log.SwitchAliLog(0);

  // -- Initialize HLT
  // -------------------
  AliHLTSystem gHLT;
  gHLT.SetGlobalLoggingLevel(0x7D);
  
  gHLT.LoadComponentLibraries("libAliHLTUtil.so");

  gHLT.LoadComponentLibraries("liblhapdf.so");  
  gHLT.LoadComponentLibraries("libEGPythia6.so");   
  gHLT.LoadComponentLibraries("libpythia6.so");  
  gHLT.LoadComponentLibraries("libAliPythia6.so");  

  // ----------------------------//
  // -                         - //
  // -- Parameters            -- //
  // -                         - //
  // ----------------------------//

  Float_t coneRadius = 0.7;
  Float_t cutPtSeed  = 4.0;
  Float_t cutEtJet   = 7.0;

  Int_t seed = 12345;

  // ----------------------------//
  // -                         - //
  // -- Publisher  Components -- //
  // -                         - //
  // ----------------------------//
 
  TString generatorId( Form("MCGenerator_%s", pprRunName[runType]) );
  
  TString generatorArg( Form("-seed %d -nevents %d -runtype %d -coneRadius %.1f -jetCutMinEt %.1f", 
			     seed, nEvents, runType, coneRadius, cutEtJet));
  
  AliHLTConfiguration mcGenerator(generatorId.Data(), "MCGenerator", NULL, generatorArg.Data() );
  
  if (!writerInput.IsNull()) writerInput+=" ";
  writerInput += generatorId;

  // ----------------------------//
  // -                         - //
  // --    Sink Components    -- //
  // -                         - //
  // ----------------------------//
  
  TString writerArg( Form("-directory analysis -datafile analyze_%d_%s -write-all-events", nEvents, pprRunName[runType] ));

  // -- The RootFileWriter 
  AliHLTConfiguration rootWriter("RootWriter", "ROOTFileWriter", 
				 writerInput.Data(), writerArg.Data() );
  
  // --------------------------- //
  // -                         - //
  // --         Run           -- //
  // -                         - //
  // --------------------------- //

  gHLT.BuildTaskList("RootWriter");
  gHLT.Run(nEvents);

  return;
}
