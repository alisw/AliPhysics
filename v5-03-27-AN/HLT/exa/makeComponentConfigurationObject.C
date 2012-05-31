// $Id$
/**
 * @file makeComponentConfigurationObject.C
 * @brief Creation of HLT component configuration objects in OCDB
 *
 * <pre>
 * Usage:
 *  aliroot -b -q makeComponentConfigurationObject.C'("path", "param", "uri", runmin, runmax)'
 *  aliroot -b -q makeComponentConfigurationObject.C'("path", "key", "param", "uri", runmin, runmax)'
 * </pre>
 *
 * Create an OCDB entry with a TObjString containing param.
 * Many HLT components understand configuration strings containing
 * arguments and parameters just like the command line arguments.
 * This macro facilitates the creation of an appropriate object
 * from a parameter string.
 * As another approach the TObjString parameters are stored in a TMap
 * associated to a key. A TMap object is generated if 'key' is specified.
 *
 * Parameters: <br>
 * - path           path of the entry within the OCDB
 * - param (opt)    string to be stored in the TObjSting, default empty
 * - uri   (opt)    the OCDB URI, default $ALICE_ROOT   
 * - runmin (opt)   default 0
 * - runmax (opt)   default 999999999
 *
 * Note: The configuration procedure of an HLT component is not
 * restricted to that scheme. The implementation is up to the
 * developer and more complex objects are possible.
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_tutorial
 */
void makeComponentConfigurationObject(const char* path, 
				      const char* key,
				      const char* param,
				      const char* cdbUri,
				      int runmin=0,
				      int runmax=AliCDBRunRange::Infinity(),
				      int runNo=0)
{
  AliCDBManager* man = AliCDBManager::Instance();
  if (!man) {
    cerr << "can not get AliCDBManager" << end;
    exit;
  }
  TString storage;
  if (!man->IsDefaultStorageSet()) {
    if (cdbUri) {
      storage=cdbUri;
      if (storage.Contains("://")==0) {
	storage="local://"; storage+=cdbUri;
      }
    } else {
      storage="local://$ALICE_ROOT/OCDB";
    }
    man->SetDefaultStorage(storage);
  } else {
    storage = man->GetDefaultStorage()->GetURI();
  }

  TMap* pMap=NULL;

  // load existing object and init TMap
  AliCDBEntry* pExisting=NULL;
  AliCDBStorage* pStorage=AliCDBManager::Instance()->GetDefaultStorage();
  if (key && pStorage->GetLatestVersion(path, runNo)>=0) {
    pExisting=pStorage->Get(path, runNo);
    if (pExisting->GetObject()->IsA() == TMap::Class()) {
      pMap=(TMap*)pExisting->GetObject()->Clone();
    }
  }  

  if (key && !pMap) pMap=new TMap;

  // here is the actual content of the configuration object
  TObject* obj=new TObjString(param);
  if (pMap) {
    if (pMap->FindObject(key)) {
      pMap->Remove(new TObjString(key));
    }
    pMap->Add(new TObjString(key), obj);
    obj=pMap;
  }

  AliCDBPath cdbPath(path);
  AliCDBId cdbId(cdbPath, runmin, runmax);
  AliCDBMetaData* cdbMetaData=NULL;
  if (pExisting) cdbMetaData=pExisting->GetMetaData();
  else cdbMetaData=new AliCDBMetaData;
  man->Put(obj, cdbId, cdbMetaData);
}

void makeComponentConfigurationObject(const char* path, const char* param="",
				      const char* cdbUri=NULL,
				      int runmin=0,
				      int runmax=AliCDBRunRange::Infinity())
{
  makeComponentConfigurationObject(path, NULL, param, cdbUri, runmin, runmax);
}

void makeComponentConfigurationObject(const char* path, 
				      int runNo,
				      const char* key,
				      const char* param)
{
  makeComponentConfigurationObject(path, key, param, NULL, 0, AliCDBRunRange::Infinity(), runNo);
}

void makeComponentConfigurationObject()
{
  cout << "===============================================================" << endl;
  cout << "usage: aliroot -b -q -l makeComponentConfigurationObject.C'(\"path\", \"param\", \"uri\", rangemin, rangemax)'" << endl << endl;
  cout << "  path           path of the entry within the OCDB" << endl;
  cout << "  param (opt)    string to be stored in the TObjSting, default empty" << endl;
  cout << "  uri   (opt)    the OCDB URI, default $ALICE_ROOT/OCDB   " << endl;
  cout << "  rangemin (opt) default 0" << endl;
  cout << "  rangemax (opt) default AliCDBRunRange::Infinity()" << endl;
  cout << "===============================================================" << endl;
  cout << "usage: aliroot -b -q -l makeComponentConfigurationObject.C'(\"path\", \"key\", \"param\", \"uri\", rangemin, rangemax)'" << endl << endl;
  cout << "  path           path of the entry within the OCDB" << endl;
  cout << "  param (opt)    string to be stored in the TObjSting, default empty" << endl;
  cout << "  uri   (opt)    the OCDB URI, default $ALICE_ROOT/OCDB   " << endl;
  cout << "  rangemin (opt) default 0" << endl;
  cout << "  rangemax (opt) default AliCDBRunRange::Infinity()" << endl;
  cout << "===============================================================" << endl;
}
