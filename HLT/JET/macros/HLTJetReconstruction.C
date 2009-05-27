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

// ---------------------------------------------------------------------------- 


/** HLTJetReconstruction test macro
 *  @param nEvents Number of events which should be processed
 */
void HLTJetReconstruction(Int_t nEvents=1, Int_t idx=0, Bool_t generate=kFALSE, PprRun_t runType = kPythia6Jets104_125 ) {

  TString writerInput;
  TString analysisInput;
  TString jetInput;

  // -- Switch Logging
  // -------------------
  AliLog::SetGlobalLogLevel( AliLog::kError );
  AliHLTLogging log;
  log.SwitchAliLog(0);

  // -- Initialize HLT
  // -------------------
  AliHLTSystem gHLT;
  gHLT.SetGlobalLoggingLevel(0x7D);

  gHLT.LoadComponentLibraries("libCGAL.so");
  gHLT.LoadComponentLibraries("libfastjet.so");

  gHLT.LoadComponentLibraries("libESD.so");  
  gHLT.LoadComponentLibraries("libSTEER.so");  
  gHLT.LoadComponentLibraries("libSTEERBase.so");  
  gHLT.LoadComponentLibraries("libAOD.so");  
  gHLT.LoadComponentLibraries("libANALYSIS.so");  
  gHLT.LoadComponentLibraries("libANALYSISalice.so");  
  gHLT.LoadComponentLibraries("libJETAN.so");  

  gHLT.LoadComponentLibraries("libAliHLTUtil.so");
  gHLT.LoadComponentLibraries("libAliHLTJET.so");  

  gHLT.LoadComponentLibraries("liblhapdf.so");  
  gHLT.LoadComponentLibraries("libEGPythia6.so");   
  gHLT.LoadComponentLibraries("libpythia6.so");  
  gHLT.LoadComponentLibraries("libAliPythia6.so");  

  // ----------------------------//
  // -                         - //
  // -- Parameters            -- //
  // -                         - //
  // ----------------------------//

  Float_t aConeRadius[] = { 0.4, 0.7 };
  Float_t aCutPtSeed[]  = { 4.0, 7.0, 10.0 };
  Float_t aCutEtJet[]   = { 4.0, 7.0, 10.0, 15.0 };

  Float_t coneRadius = 0.7;
  Float_t cutPtSeed  = 4.0;
  Float_t cutEtJet   = 7.0;

  Int_t seed = 12345;

  // ----------------------------//
  // -                         - //
  // -- Publisher  Components -- //
  // -                         - //
  // ----------------------------//

  // ------------------------------------------
  // -- The ESDMCEventPublisher   
  // ------------------------------------------
  if ( ! generate ) {

    TString publisherId("ESDMCEventPublisher");

    TString publisherArg( Form("-entrytype MCFAST -dataspec 0x0000001F -datapath /home/jthaeder/jet/data/HEAD_2009-03-17/FastGen/kPythia6Jets104_125_14TeV/JET-ETA=-0.2,0.2_JET-ET=50,1000_R=0.4_10ev") );
    
    AliHLTConfiguration ESDMCEventPublisher(publisherId.Data(), "ESDMCEventPublisher", NULL, publisherArg.Data() );

    if (!analysisInput.IsNull()) analysisInput+=" ";
    analysisInput += publisherId;

    if (!jetInput.IsNull()) jetInput+=" ";
    jetInput += publisherId;
  }

  // ------------------------------------------
  // -- The MCGenerator
  // ------------------------------------------
  else {
    
    TString generatorId( Form("MCGenerator_%s", pprRunName[runType]) );

    TString generatorArg( Form("-seed %d -nevents %d -runtype %d -coneRadius %.1f -jetCutMinEt %.1f", 
			       seed, nEvents, runType, coneRadius, cutEtJet));
    
    AliHLTConfiguration mcGenerator(generatorId.Data(), "MCGenerator", NULL, generatorArg.Data() );

    if (!analysisInput.IsNull()) analysisInput+=" ";
    analysisInput += generatorId;

    if (!jetInput.IsNull()) jetInput+=" ";
    jetInput += generatorId;
  }

  // ----------------------------//
  // -                         - //
  // -- Processing Components -- //
  // -                         - //
  // ----------------------------//
#if 1
  // ------------------------------------------
  // -- ConeJetFinder
  // ------------------------------------------

  TString jetId("JETConeJet");

  TString jetArg( Form("-coneRadius %.1f -trackCutMinPt 0.0 -seedCutMinPt %.1f -jetCutMinEt %.1f",
		       coneRadius, cutPtSeed, cutEtJet) );
  
  AliHLTConfiguration jetCone(jetId.Data(), "JETConeJetFinder", jetInput.Data(), jetArg.Data()); 
  
  if (!analysisInput.IsNull()) analysisInput+=" ";
  analysisInput += jetId;
  
#else
  // ------------------------------------------
  // -- FastJetFinder
  // ------------------------------------------

  AliHLTConfiguration jetFinder("JETFastJet_Kt", "JETFastJetFinder",
				jetInput.Data(),"-finderType kt");
    
  if (!writerInput.IsNull()) writerInput+=" ";
  writerInput+="JETFastJet_Kt";
  
  
  AliHLTConfiguration jetFinder("JETFastJet_AntiKt", "JETFastJetFinder",
				jetInput.Data(),"-finderType antikt");
  
  if (!writerInput.IsNull()) writerInput+=" ";
  writerInput+="JETFastJet_AntiKt";
#endif
  
  // ------------------------------------------
  // -- Jet Analysis 
  // ------------------------------------------

  TString analysisId("JETAnalysis");
  TString analysisArg("");
  
  AliHLTConfiguration jetAnalysis(analysisId.Data(), "JETAnalysis", analysisInput.Data(), analysisArg.Data() );

  if (!writerInput.IsNull()) writerInput+=" ";
  writerInput += analysisId;
  
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
