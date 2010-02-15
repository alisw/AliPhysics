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

  AliCDBMetaData* meta=cdbEntry->GetMetaData();
  if (comment) meta->SetComment(comment);
  if (responsible) meta->SetResponsible(responsible);
  if (alirootv) meta->SetAliRootVersion(alirootv);

  AliCDBManager* man = AliCDBManager::Instance();
  man->SetDefaultStorage(targetCDB);
  man->Put(cdbEntry);

}
