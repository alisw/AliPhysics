// $Id$
/**
 * @file sampleCalibrationProcessor.C
 * @brief Run the SampleCalibration component (AliHLTSampleCalibrationComponent)
 *
 * <pre>
 * Usage: aliroot -b -q -l \
 *     sampleCalibrationProcessor.C'("file", "cdb", minEvent, maxEvent)'
 *
 * Examples:
 *     sampleCalibrationProcessor.C'("alien:///alice/data/2009/.../....root")' 
 *     sampleCalibrationProcessor.C'("raw.root", "local://$PWD", minEvent, MaxEvent)'
 *
 * Defaults
 *     cdb="local://$ALICE_ROOT/OCDB"  -> take local OCDB
 *     minEvent=-1   -> no lower event selection
 *     maxEvent=-1   -> no upper event selection
 *
 * </pre>
 *
 * An HLT chain is set up with one publisher for raw data which the calibration
 * component subscribes to. The output of the component is directed to a
 * ROOTFileWriter producing the file 'calibdata.root'. Another FileWriter
 * writes the data block designated for storage in the FXS. 
 * The HLT chain is run in the context of AliReconstruction.
 * 
 * The input file can be a file on the grid, indicated by the tag
 * 'alien://'.
 * If either the file or the OCDB is taken from the GRID, the macros
 * connects to the Grid in the beginning.
 * Note: You need a valid GRID token, use 'alien-token-init' of your
 * alien installation.
 *
 * This example uses the local OCDB as default, but other locations can
 * be specified, like e.g. "raw://".
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_tutorial
 */
void sampleCalibrationProcessor(const char *filename,
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

  // Reconstruction settings
  AliReconstruction rec;

  if (minEvent>=0 || maxEvent>minEvent) {
    if (minEvent<0) minEvent=0;
    if (maxEvent<minEvent) maxEvent=minEvent;
    rec.SetEventRange(minEvent,maxEvent);
  }

  rec.SetRunReconstruction("HLT");

  // QA options
  rec.SetRunQA(":") ;

  // AliReconstruction settings
  rec.SetWriteESDfriend(kFALSE);
  rec.SetInput(filename);

  rec.SetRunPlaneEff(kFALSE);
  rec.SetRunVertexFinder(kFALSE);
  rec.SetRunMultFinder(kFALSE);

  // switch off cleanESD
  rec.SetCleanESD(kFALSE);

  // setup the HLT chain
  AliHLTSystem* pHLT=AliHLTPluginBase::GetInstance();
  
  // publisher for raw data
  AliHLTConfiguration rawpub("rawpub", "AliRawReaderPublisher", "", "-minid 768 -maxid 983 -datatype 'DDL_RAW ' SMPL");

  // the configuration of the calibration component
  // send the produced object not more than once every 2 sec
  AliHLTConfiguration smplcalib("smplcalib", "SampleCalibration", "rawpub", "-pushback-period=2");

  // writer for data objects produced by the calibration component during
  // ProcessCalibration, 
  AliHLTConfiguration eventwriter("eventwriter", "ROOTFileWriter", "smplcalib", "-datafile calibdata.root -concatenate-events -overwrite -write-all-events");

  // filter the output of the smplcalib component and forward only the blocks
  // of type 'FXS_CAL ' designated for storage by the FXSSubscriber
  AliHLTConfiguration fxsfilter("fxsfilter", "BlockFilter", "smplcalib", "-typeid 'FXS_CAL ' ");
  // write the block, a file like 'EOR_*HLT:FXS_CAL' will be written containing the
  // binary data block of the FXS header and the streamed object
  AliHLTConfiguration fxswriter("fxswriter", "FileWriter", "fxsfilter", "-write-all-events");

  rec.SetOption("HLT", "loglevel=0x7c ignore-hltout libAliHLTUtil.so libAliHLTSample.so chains=fxswriter,eventwriter");

  AliLog::Flush();
  rec.Run();

}

void sampleCalibrationProcessor(const char *filename,
		  int minEvent=-1,
		  int maxEvent=-1)
{
  sampleCalibrationProcessor(filename, "local://$ALICE_ROOT/OCDB", minEvent, maxEvent);
}

void sampleCalibrationProcessor()
{
  cout << "sampleCalibrationProcessor: Run AliRoot reconstruction locally" << endl;
  cout << " Usage: aliroot -b -q -l \\" << endl;
  cout << "     sampleCalibrationProcessor.C'(\"file\", \"cdb\", minEvent, maxEvent)'" << endl;
  cout << "" << endl;
  cout << " Examples:" << endl;
  cout << "     sampleCalibrationProcessor.C'(\"alien:///alice/data/2009/.../....root\")' " << endl;
  cout << "     sampleCalibrationProcessor.C'(\"raw.root\", \"local://$PWD\", minEvent, MaxEvent)'" << endl;
  cout << "" << endl;
  cout << " Defaults" << endl;
  cout << "     cdb=\"local://$ALICE_ROOT/OCDB\"  -> take local OCDB " << endl;
  cout << "     minEvent=-1   -> no lower event selection" << endl;
  cout << "     maxEvent=-1   -> no upper event selection" << endl;
}
