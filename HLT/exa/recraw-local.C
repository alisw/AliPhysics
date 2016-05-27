// $Id$
/**
 * @file recraw-local.C
 * @brief Run reconstruction of raw data locally
 *
 * <pre>
 * Usage: aliroot -b -q -l \
 *     recraw-local.C'("file", "cdb", minEvent, maxEvent, modules)'
 *
 * Examples:
 *     recraw-local.C'("alien:///alice/data/2009/.../....root")' 
 *     recraw-local.C'("raw://run12345")'
 *     recraw-local.C'("raw://run12345", minEvent, MaxEvent)'
 *     recraw-local.C'("raw.root", "local://$PWD", minEvent, MaxEvent)'
 *     recraw-local.C'("10000135031045.10..root", "raw://", -1, -1)'
 *     recraw-local.C'("alien:///alice/data/2010/LHC10f/000135031/raw/10000135031045.10.root", "raw://", -1, -1)'
 *
 * Defaults
 *     cdb="raw://"  -> take OCDB from GRID
 *     minEvent=-1   -> no lower event selection
 *     maxEvent=-1   -> no upper event selection
 *     modules="ALL" -> all modules
 *     hltOption="loglevel=0x7c" -> logging level info and above
 *
 * </pre>
 *
 * The input file can be a file on the grid, indicated by the tag
 * 'alien://' indicates. By default also the OCDB is set to the GRID.
 * If either the file or the OCDB is taken from the GRID, the macros
 * connects to the Grid in the beginning.
 *
 * Input files can be specified via te run number when using the tag
 * 'raw://' followed by the string 'run12345' where the number needs
 * to be adjusted.
 *
 * As for the OCDB it is always a good idea to use the OCDB from the
 * Grid as this will contain all the necessary objects and the latest
 * calibration. The special URI 'raw://' is most advisable as it selects
 * the storage automatically from the run number. Other options are e.g.
 * - "alien://folder=/alice/data/2010/OCDB"
 * - "local://$ALICE_ROOT/OCDB"
 *
 * Re-running the HLT reconstruction
 * By specifying the hlt options, the HLT chain can be re-run instead
 * of just extracting the online result. E.g. the following options
 * specify to ignore the HLTOUT payload and run the two chains defined
 * in the agents. The translation of the online configuration into
 * an HLT offline chain is under development.
 * <pre>
 *   ignore-hltout chains=GLOBAL-esd-converter,TPC-clusters
 * <pre>
 *
 * Note: You need a valid GRID token, use 'alien-token-init' of your
 * alien installation.
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_qa
 */
void recraw_local(const char *filename,
    const char *cdbURI,
    int minEvent=-1,
    int maxEvent=-1,
    const char *modules="ALL",
    const char *hltOptions="loglevel=0x7c",
    const char *cdbDrain=NULL,
    const char* cdbSnapshot=NULL)
{
  if(!gSystem->AccessPathName("galice.root")){
    cerr << "AliReconstruction on raw data requires to delete galice.root, or run at different place." << endl;
    cerr << "!!! DO NOT DELETE the galice.root of your simulation, but create a subfolder !!!!" << endl;
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
  
  //if found, configure any specific storages
  //same logic as in CPassX and benchmark
  if (gROOT->LoadMacro("localOCDBaccessConfig.C")==0) {
    localOCDBaccessConfig();
  }  

  if (struri.BeginsWith("local://")) {
    // set specific storage for GRP entry
    // search in the working directory and one level above, the latter
    // follows the standard simulation setup like e.g. in test/ppbench
    if (!gSystem->AccessPathName("GRP/GRP/Data")) {
      man->SetSpecificStorage("GRP/GRP/Data", "local://$PWD");
    } else if (!gSystem->AccessPathName("../GRP/GRP/Data")) {
      man->SetSpecificStorage("GRP/GRP/Data", "local://$PWD/..");      
    }
  }
  if (cdbDrain) man->SetDrain(cdbDrain);

  if (cdbSnapshot) {
    man->SetSnapshotMode(cdbSnapshot);
  }

  // Reconstruction settings
  AliReconstruction rec;

  if (minEvent>=0 || maxEvent>minEvent) {
    if (minEvent<0) minEvent=0;
    if (maxEvent<minEvent) maxEvent=minEvent;
    rec.SetEventRange(minEvent,maxEvent);
  }

  TString strModules=modules;
  if (modules)
    rec.SetRunReconstruction(modules);
  else
    rec.SetRunReconstruction("ALL");

  // QA options
  TString qaOptions="HLT TPC";
  if (!strModules.Contains("TPC")) qaOptions.ReplaceAll("TPC", "");
  qaOptions+=":ALL";
  rec.SetRunQA(qaOptions) ;
  //rec.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;

  // AliReconstruction settings
  rec.SetWriteESDfriend(kTRUE);
  rec.SetRunVertexFinder(strModules.Contains("ITS"));
  rec.SetRunMultFinder(strModules.Contains("ITS"));
  rec.SetInput(filename);
  rec.SetOption("HLT", hltOptions);

  rec.SetRunPlaneEff(kFALSE);

  // switch off cleanESD
  rec.SetCleanESD(kFALSE);

  AliLog::Flush();
  rec.Run();

}

void recraw_local(const char *filename,
		  int minEvent=-1,
		  int maxEvent=-1,
		  const char *modules="ALL",
		  const char *hltOptions="loglevel=0x7f")
{
  recraw_local(filename, "raw://", minEvent, maxEvent, modules, hltOptions);
}

void recraw_local()
{
  cout << "recraw-local: Run AliRoot reconstruction locally" << endl;
  cout << " Usage: aliroot -b -q -l \\" << endl;
  cout << "     recraw-local.C'(\"file\", \"cdb\", minEvent, maxEvent, modules, hltOptions)'" << endl;
  cout << "" << endl;
  cout << " Examples:" << endl;
  cout << "     recraw-local.C'(\"alien:///alice/data/2009/.../....root\")' " << endl;
  cout << "     recraw-local.C'(\"raw://run12345\")'" << endl;
  cout << "     recraw-local.C'(\"raw://run12345\", minEvent, MaxEvent)'" << endl;
  cout << "     recraw-local.C'(\"raw.root\", \"local://$PWD\", minEvent, MaxEvent)'" << endl;
  cout << "" << endl;
  cout << " Defaults" << endl;
  cout << "     cdb=\"raw://\"  -> take OCDB from GRID" << endl;
  cout << "     minEvent=-1   -> no lower event selection" << endl;
  cout << "     maxEvent=-1   -> no upper event selection" << endl;
  cout << "     modules=\"ALL\" -> all modules" << endl;
  cout << "     hltOption=\"loglevel=0x7c\" -> logging level info and above" << endl;
}
