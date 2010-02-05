// $Id$
/**
 * @file recraw-local.C
 * @brief Run reconstruction of raw data locally
 *
 * <pre>
 * Usage: aliroot -b -q -l \
 *     recraw-local.C'("file", "cdb", minEvent, maxEvent)'
 *
 * Examples:
 *     recraw-local.C'("alien:///alice/data/2009/.../.root")' 
 *     recraw-local.C'("raw://run12345")'
 *     recraw-local.C'("raw://run12345", minEvent, MaxEvent)'
 *     recraw-local.C'("raw.root", "local://$PWD", minEvent, MaxEvent)'
 *
 * Defaults
 *     file=raw.root"
 *     cdb="raw://"
 *     minEvent=-1
 *     maxEvent=-1
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
 * Note: You need a valid token, use 'alien-token-init' of your alien
 * installation.
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_qa
 */
void recraw_local(const char *filename,
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

  // Reconstruction settings
  AliReconstruction rec;

  if (minEvent>=0 || maxEvent>minEvent) {
    if (minEvent<0) minEvent=0;
    if (maxEvent<minEvent) maxEvent=minEvent;
    rec.SetEventRange(minEvent,maxEvent);
  }

  rec.SetRunReconstruction("ALL");

  // QA options
  rec.SetRunQA("ALL:ALL") ;
  rec.SetQARefDefaultStorage("local://$ALICE_ROOT/QAref") ;

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

void recraw_local(const char *filename="raw.root",
		  int minEvent=-1,
		  int maxEvent=-1)
{
  recraw_local(filename, "raw://", minEvent, maxEvent);
}
