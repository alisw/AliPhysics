// $Id$
/**
 * @file makeComponentConfigurationObject.C
 * @brief Creation of HLT component configuration objects in OCDB
 *
 * <pre>
 * Usage: aliroot -b -q makeComponentConfigurationObject.C'("path", "param", "uri", runmin, runmax)'
 * </pre>
 *
 * Create an OCDB entry with a TObjString containing param.
 * Many HLT components understand configuration strings containing
 * arguments and parameters just like the command line arguments.
 * This macro facilitates the creation of an appropriate object
 * from a parameter string.
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
void makeComponentConfigurationObject(const char* path, const char* param="",
				      const char* cdbUri=NULL,
				      int runmin=0,
				      int runmax=999999999)
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
      if (storage.Contains("://")!=0) {
	storage="local://"; storage+=cdbUri;
      }
    } else {
      storage="local://$ALICE_ROOT/OCDB";
    }
    man->SetDefaultStorage(storage);
  } else {
    storage = man->GetDefaultStorage()->GetURI();
  }

  // here is the actual content of the configuration object
  TObjString obj=param;
  AliCDBPath cdbPath(path);
  AliCDBId cdbId(cdbPath, runmin, runmax);
  AliCDBMetaData cdbMetaData;
  man->Put(&obj, cdbId, &cdbMetaData);
  cout << "adding TObjString type OCDB object " << path << " (" << (param[0]==0?"<empty>":param) << ") [" << runmin << "," << runmax << "] in " << storage << endl;
}

void makeComponentConfigurationObject()
{
  cout << "===============================================================" << endl;
  cout << "usage: aliroot -b -q -l makeComponentConfigurationObject.C'(\"path\", \"param\", \"uri\", rangemin, rangemax)'" << endl << endl;
  cout << "  path           path of the entry within the OCDB" << endl;
  cout << "  param (opt)    string to be stored in the TObjSting, default empty" << endl;
  cout << "  uri   (opt)    the OCDB URI, default $ALICE_ROOT   " << endl;
  cout << "  rangemin (opt) default 0" << endl;
  cout << "  rangemax (opt) default 999999999" << endl;
  cout << "===============================================================" << endl;
}
