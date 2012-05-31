// $Id$
/**
 * @file adjustOCDBObject.C
 * @brief Tool to adjust properties of an OCDB object and write it back to OCDB.
 *
 * Usage:
 * <pre>
 *   aliroot -b -q adjustOCDBObject("file",
 *                                  targetCDB="local://$PWD",
 *                                  firstrun=-1,
 *                                  lastrun=-1,
 *                                  version=-1,
 *                                  subversion=-1,
 *                                  comment=NULL,
 *                                  responsible=NULL,
 *                                  alirootv=NULL) 
 * </pre>
 *
 * The macro opens an OCDB entry directly as a file, changes properties
 * of the entry and uses the CDBmanager to write it back to some new
 * location.
 *
 * Required parameters:
 * - filename
 * - OCDB target
 * Optional parameters, if the default values are specified the
 * corresponding key in the entry is not changed
 * - first run
 * - last run, if smaller than first run 'AliCDBRunRange::Infinity' is used
 * - version
 * - subversion
 * - comment
 * - responsible
 * - aliroot version 
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_tutorial
 */
void adjustOCDBObject(const char* file,
		      const char* targetCDB="local://$PWD",
		      int firstrun=-1,
		      int lastrun=-1,
		      int version=-1,
		      int subversion=-1,
		      const char* comment=NULL,
		      const char* responsible=NULL,
		      const char* alirootv=NULL) 
{
  TFile* origfile=new TFile(file);
  if (origfile->IsZombie()) {
    cerr << "error opening file " << file << endl;
    return;
  }
  
  AliCDBEntry* cdbEntry=NULL;
  origfile->GetObject("AliCDBEntry", cdbEntry);
  if (!cdbEntry) {
    cerr << "can not find CDB entry in file " << file << endl;
  }

  AliCDBId& cdbId=cdbEntry->GetId();
  if (version>=0) cdbId.SetVersion(version);
  if (subversion>=0) cdbId.SetSubVersion(subversion);
  if (firstrun>=0) cdbId.SetFirstRun(firstrun);
  if (lastrun>=0) cdbId.SetLastRun(lastrun);
  if (lastrun<firstrun) cdbId.SetLastRun(AliCDBRunRange::Infinity());

  AliCDBMetaData* meta=cdbEntry->GetMetaData();
  if (comment) meta->SetComment(comment);
  if (responsible) meta->SetResponsible(responsible);
  if (alirootv) meta->SetAliRootVersion(alirootv);

  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(targetCDB);
  man->Put(cdbEntry);

}

void adjustOCDBObject() 
{
  cout << "adjustOCDBObject.C   adjust properties of an OCDB entry" << endl;
  cout << "usage:" << endl;
  cout << "   aliroot -b -q adjustOCDBObject'(\"file\", " << endl;
  cout << "                                  targetCDB=\"local://$PWD\", " << endl;
  cout << "                                  firstrun=-1, " << endl;
  cout << "                                  lastrun=-1, " << endl;
  cout << "                                  version=-1, " << endl;
  cout << "                                  subversion=-1, " << endl;
  cout << "                                  comment=NULL, " << endl;
  cout << "                                  responsible=NULL, " << endl;
  cout << "                                  alirootv=NULL)'  " << endl;
  cout << " e.g." << endl;
  cout << "   aliroot -b -q adjustOCDBObject'(\"myfile.root\", \"$ALICE_ROOT/OCDB\", 42, 42)'" << endl;
}
