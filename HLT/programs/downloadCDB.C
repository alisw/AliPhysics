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
  man->SetDrain(to);
  AliCDBPath cdbpath(path);
  man->GetAll(path, runnr);
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
