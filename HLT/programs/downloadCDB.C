// $Id$
///
/// @file downloadCDB.C
/// @brief Download a snapshot of the CDB to another folder
///
/// <pre>
/// Usage: aliroot -b -q -l \
///     downloadCDB.C'(runno, "from", "to", "path")'
///
/// Examples:
///     downloadCDB.C'(144991, "alien://folder=/alice/data/2011/OCDB", "local:///tmp/144991/OCDB")'
///
/// Defaults
///     path="*/*/*"  -> download everything
///
/// </pre>
///
///
void downloadCDB(Int_t runnr,
      const char* from,
      const char* to,
      const char* path="*/*/*") 
{
  AliCDBManager* man=AliCDBManager::Instance();
  man->SetDefaultStorage(from);
  man->SetRun(runnr);
  man->SetDrain(to);
  AliCDBPath cdbpath(path);
  man->GetAll(path, runnr);

  if(TString(getenv("OCDB_SNAPSHOT_CREATE")) == TString("kTRUE"))
  {
    TString snapshotFile(getenv("OCDB_SNAPSHOT_FILENAME"));

    TString snapshotFileOut = "OCDB.root";
    if(!(snapshotFile.IsNull() || snapshotFile.IsWhitespace()))
    {
      snapshotFileOut = snapshotFile;
    }

    man->DumpToSnapshotFile("OCDB.root",kFALSE);
  }
}

void downloadCDB()
{
  cout << " Usage:" << endl;
  cout << "   aliroot -b -q -l \\" << endl;
  cout << "     downloadCDB.C'(runno, \"from\", \"to\", \"path\")'" << endl;
  cout << "" << endl;
  cout << " Examples:" << endl;
  cout << "   aliroot -b -q -l \\" << endl;
  cout << "     downloadCDB.C'(144991, \"alien://folder=/alice/data/2011/OCDB\", \"local:///tmp/144991/OCDB\")'" << endl;
  cout << "" << endl;
  cout << " Defaults" << endl;
  cout << "     path=\"*/*/*\"  -> download everything" << endl;
}
