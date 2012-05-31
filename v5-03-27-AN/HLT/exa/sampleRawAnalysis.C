// $Id$
void sampleRawAnalysis();
/**
 * @file sampleRawAnalysis.C
 * @brief Example macro to run the AliHLTSampleRawAnalysisComponent in
 * AliReconstruction.
 *
 * The component subscribes to DDL raw data published by the
 * AliHLTRawReaderPublisherComponent. The macros requires a raw data file
 * and a corresponding GRP entry.
 *
 * <pre>
 * Usage: aliroot -b -q -l \
 *     sampleRawAnalysis.C'("rawfile", "cdb", minEvent, maxEvent)'
 *
 * Examples:
 *     sampleRawAnalysis.C'("raw.root", minEvent, MaxEvent)'
 *     sampleRawAnalysis.C'("./", minEvent, MaxEvent)'
 *     sampleRawAnalysis.C'("alien:///alice/data/2010/.../raw/....root")' 
 *
 * Defaults
 *     cdb="local://$ALICE_ROOT/OCDB" -> take local OCDB from distribution
 *     minEvent=-1   -> no lower event selection
 *     maxEvent=-1   -> no upper event selection
 *
 * </pre>
 *
 * The input file can be a local raw.root file but also a file from the
 * GRID. The separate DDL files generated in simulation can be accessed
 * using AliRawReaderFile by speficying "directory/".
 *
 * Since the macro runs AliReconstruction the OCDB needs to be set up, in
 * particular the GRP entry. If testing with a local OCDB you have to
 * simulate some events and run the macro in the folder of the simulation.
 * Also HLT components configure from objects in the OCDB.
 *
 * Note: You need a valid GRID token, if you want to access files directly
 * from the Grid, use 'alien-token-init' of your alien installation.
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_tutorial
 */
void sampleRawAnalysis(const char *filename,
		       const char *cdbURI,
		       int minEvent=-1,
		       int maxEvent=-1)
{
  if(!gSystem->AccessPathName("galice.root")){
    cerr << "AliReconstruction on raw data requires to delete galice.root, or run at different place." << endl;
    cerr << "!!! DO NOT DELETE the galice.root of your simulation, but create a subfolder !!!!" << endl;
    return;
  }

  if (gSystem->AccessPathName(filename)) {
    cerr << "can not find file " << filename << endl;
    return;
  }

  // connect to the GRID if we use a file or OCDB from the GRID
  TString struri=cdbURI;
  TString strfile=filename;
  if (struri.BeginsWith("raw://") ||
      strfile.Contains("://") && !strfile.Contains("local://")) {
    TGrid::Connect("alien");
  }

  // Set the CDB storage location
  AliCDBManager * man = AliCDBManager::Instance();
  man->SetDefaultStorage(cdbURI);
  if (struri.BeginsWith("local://")) {
    // set specific storage for GRP entry
    // search in the working directory and one level above, the latter
    // follows the standard simulation setup like e.g. in test/ppbench
    if (!gSystem->AccessPathName("GRP/GRP/Data")) {
      man->SetSpecificStorage("GRP/GRP/Data", "local://$PWD");
    } else if (!gSystem->AccessPathName("../GRP/GRP/Data")) {
      man->SetSpecificStorage("GRP/GRP/Data", "local://$PWD/..");      
    } else {
      cerr << "can not find a GRP entry, please run the macro in the folder" << endl;
      cerr << "of a simulated data sample, or specify a GRID OCDB" << endl;
      sampleRawAnalysis();      
      return;
    }
  }

  //////////////////////////////////////////////////////////////////////////////////////
  //
  // Reconstruction settings
  AliReconstruction rec;

  if (minEvent>=0 || maxEvent>minEvent) {
    if (minEvent<0) minEvent=0;
    if (maxEvent<minEvent) maxEvent=minEvent;
    rec.SetEventRange(minEvent,maxEvent);
  }
  rec.SetRunReconstruction("HLT");
  rec.SetLoadAlignFromCDB(kFALSE);
  rec.SetWriteESDfriend(kFALSE);
  // due to bug ...
  // StopOnError needs to be disabled
  rec.SetStopOnError(kFALSE);
  rec.SetRunVertexFinder(kFALSE);
  rec.SetInput(filename);

  // QA options
  rec.SetRunQA(":") ;
  //rec.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;

  //////////////////////////////////////////////////////////////////////////////////////
  //
  // setup the HLT system
  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();

  // define a data publisher configuration
  // arguments:
  //  1) id of the configuartion, later used to refer to this configuration
  //  2) id of the component to run
  //  3) parents, no parents for publisher components
  //  4) optional component arguments
  //              publish link #512, the first link (numbered 0) of SSD
  AliHLTConfiguration publisher("RAW-Publisher", 
				"AliRawReaderPublisher", 
				"", 
				"-minid 512 "
				"-datatype 'DDL_RAW ' 'ISSD' "
				"-dataspec 0x01"
				);

  // define configuration of the SampleRawAnalyis component
  // arguments:
  //  1) id of the configuartion, later used to refer to this configuration
  //  2) id of the component to run
  //  3) parents, here the publisher configuration defined above
  //  4) optional component arguments
  AliHLTConfiguration rawanalysis("RAW-Analysis", 
				  "SampleRawAnalysis", 
				  "RAW-Publisher", 
				  "-mandatory1 test "
				  "-verbose"
				  );

  // set option for the HLT module in AliReconstruction
  // arguments
  //  - ignore-hltout : ignore the HLTOUT payload from the HLT DDLs
  //  - libraries to be used as plugins
  //  - loglevel=0x7c : Important, Info, Warning, Error, Fatal
  //  - chains=RAW-Analysis : chains to be run
  rec.SetOption("HLT",
		"ignore-hltout " 
		"libAliHLTUtil.so libAliHLTSample.so "
		"loglevel=0x7c "
		"chains=RAW-Analysis "
		);

  rec.SetRunPlaneEff(kFALSE);

  // switch off cleanESD
  rec.SetCleanESD(kFALSE);

  AliLog::Flush();
  rec.Run();

}

void sampleRawAnalysis(const char *filename,
		       int minEvent=-1,
		       int maxEvent=-1)
{
  sampleRawAnalysis(filename, "local://$ALICE_ROOT/OCDB", minEvent, maxEvent);
}

void sampleRawAnalysis()
{
  cout << "sampleRawAnalysis: Run HLT component 'SampleRawAnalyis' in AliReconstruction" << endl;
  cout << " Usage: aliroot -b -q -l \\" << endl;
  cout << "     sampleRawAnalysis.C'(\"file\", \"cdb\", minEvent, maxEvent)'" << endl;
  cout << "" << endl;
  cout << " Examples:" << endl;
  cout << "     sampleRawAnalysis.C'(\"raw.root\", minEvent, MaxEvent)'" << endl;
  cout << "     sampleRawAnalysis.C'(\"./\", minEvent, MaxEvent)'" << endl;
  cout << "     sampleRawAnalysis.C'(\"alien:///alice/data/2010/.../raw/....root\")' " << endl;
  cout << "" << endl;
  cout << " Defaults" << endl;
  cout << "     cdb=\"local://$ALICE_ROOT/OCDB\"  -> take local OCDB" << endl;
  cout << "     minEvent=-1   -> no lower event selection" << endl;
  cout << "     maxEvent=-1   -> no upper event selection" << endl;
}
