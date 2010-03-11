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
 *
 * Defaults
 *     cdb="raw://"  -> take OCDB from GRID
 *     minEvent=-1   -> no lower event selection
 *     maxEvent=-1   -> no upper event selection
 *     modules="ALL" -> all modules
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
		  const char *modules="ALL")
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

  // Reconstruction settings
  AliReconstruction rec;

  if (minEvent>=0 || maxEvent>minEvent) {
    if (minEvent<0) minEvent=0;
    if (maxEvent<minEvent) maxEvent=minEvent;
    rec.SetEventRange(minEvent,maxEvent);
  }

  if (modules)
    rec.SetRunReconstruction(modules);
  else
    rec.SetRunReconstruction("ALL");

  // QA options
  rec.SetRunQA("HLT TPC:ALL") ;
  //rec.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;

  // AliReconstruction settings
  rec.SetWriteESDfriend(kTRUE);
  rec.SetInput(filename);
  rec.SetOption("HLT","loglevel=0x3f");

  rec.SetRunPlaneEff(kFALSE);

  // switch off cleanESD
  rec.SetCleanESD(kFALSE);

  AliLog::Flush();
  rec.Run();

}

void recraw_local(const char *filename,
		  int minEvent=-1,
		  int maxEvent=-1,
		  const char *modules="ALL")
{
  recraw_local(filename, "raw://", minEvent, maxEvent, modules);
}

void recraw_local()
{
  cout << "recraw-local: Run AliRoot reconstruction locally" << endl;
  cout << " Usage: aliroot -b -q -l \\" << endl;
  cout << "     recraw-local.C'(\"file\", \"cdb\", minEvent, maxEvent, modules)'" << endl;
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
}
