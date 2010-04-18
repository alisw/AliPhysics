// $Id$
/**
 * @file sampleEsdAnalysis.C
 * @brief Example macro to run the AliHLTSampleESDAnalysisComponent in
 * AliReconstruction.
 *
 * The component subscribes to the output of the default HLT reconstruction
 * chain and extracts the ESD from the input.
 *
 * <pre>
 * Usage: aliroot -b -q -l \
 *     sampleEsdAnalysis.C'("file", "cdb", minEvent, maxEvent)'
 *
 * Examples:
 *     sampleEsdAnalysis.C'("alien:///alice/data/2009/.../....root")' 
 *     sampleEsdAnalysis.C'("raw.root", "local://$PWD", minEvent, MaxEvent)'
 *
 * Defaults
 *     cdb="raw://"  -> take OCDB from GRID
 *     minEvent=-1   -> no lower event selection
 *     maxEvent=-1   -> no upper event selection
 *
 * </pre>
 *
 * The input file can be a file on the grid, indicated by the tag
 * 'alien://' indicates. By default also the OCDB is set to the GRID.
 * If either the file or the OCDB is taken from the GRID, the macros
 * connects to the Grid in the beginning.
 *
 * As for the OCDB it is always a good idea to use the OCDB from the
 * Grid as this will contain all the necessary objects and the latest
 * calibration. The special URI 'raw://' is most advisable as it selects
 * the storage automatically from the run number. Other options are e.g.
 * - "alien://folder=/alice/data/2010/OCDB"
 * - "local://$ALICE_ROOT/OCDB"
 *
 * Note: You need a valid GRID token, if you want to access files directly
 * from the Grid, use 'alien-token-init' of your alien installation.
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_tutorial
 */
void sampleEsdAnalysis(const char *filename,
		       const char *cdbURI,
		       int minEvent=-1,
		       int maxEvent=-1)
{
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
  if (struri.BeginsWith("local://") && !gSystem->AccessPathName("GRP/GRP/Data")) {
    // set specific storage for GRP entry
    man->SetSpecificStorage("GRP/GRP/Data", "local://$PWD");
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
  rec.SetRunReconstruction("HLT ITS TPC");
  rec.SetWriteESDfriend(kFALSE);
  rec.SetInput(filename);

  // QA options
  rec.SetRunQA(":") ;
  //rec.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;

  //////////////////////////////////////////////////////////////////////////////////////
  //
  // setup the HLT system
  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();

  // define a configuration for the SampleESDAnalysis component
  // arguments:
  //  1) id of the configuartion, later used to refer to this configuration
  //  2) id of the component to run
  //  3) parents, here the configuration 'GLOBAL-esd-converter' of the libAliHLTGlobal
  //  4) optional component arguments
  AliHLTConfiguration esdanalysis("ESD-Analysis", "SampleESDAnalysis", "GLOBAL-esd-converter", "");

  // set option for the HLT module in AliReconstruction
  // arguments
  //  - ignore-hltout : ignore the HLTOUT payload from the HLT DDLs
  //  - libraries to be used as plugins
  //  - loglevel=0x79 : Important, Warning, Error, Fatal
  //  - chains=ESD-Analysis : chains to be run
  rec.SetOption("HLT",
		"ignore-hltout " 
		"libAliHLTUtil.so libAliHLTGlobal.so libAliHLTSample.so "
		"loglevel=0x79 "
		"chains=ESD-Analysis "
		);

  rec.SetRunPlaneEff(kFALSE);

  // switch off cleanESD
  rec.SetCleanESD(kFALSE);

  AliLog::Flush();
  rec.Run();

}

void sampleEsdAnalysis(const char *filename,
		       int minEvent=-1,
		       int maxEvent=-1)
{
  sampleEsdAnalysis(filename, "raw://", minEvent, maxEvent, modules, hltOptions);
}

void sampleEsdAnalysis()
{
  cout << "sampleEsdAnalysis: Run HLT component 'SampleESDAnalyis' in AliReconstruction" << endl;
  cout << " Usage: aliroot -b -q -l \\" << endl;
  cout << "     sampleEsdAnalysis.C'(\"file\", \"cdb\", minEvent, maxEvent)'" << endl;
  cout << "" << endl;
  cout << " Examples:" << endl;
  cout << "     sampleEsdAnalysis.C'(\"alien:///alice/data/2009/.../....root\")' " << endl;
  cout << "     sampleEsdAnalysis.C'(\"raw.root\", \"local://$PWD\", minEvent, MaxEvent)'" << endl;
  cout << "" << endl;
  cout << " Defaults" << endl;
  cout << "     cdb=\"raw://\"  -> take OCDB from GRID" << endl;
  cout << "     minEvent=-1   -> no lower event selection" << endl;
  cout << "     maxEvent=-1   -> no upper event selection" << endl;
}
