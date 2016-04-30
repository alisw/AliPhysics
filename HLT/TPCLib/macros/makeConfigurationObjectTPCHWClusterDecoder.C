//-*- Mode: C++ -*-
/**
 * @file makeConfigurationObjectTPCHWClusterDecoder.C
 * @brief Creation of HLT component configuration objects in OCDB
 *
 * <pre>
 * Usage: aliroot -b -q makeConfigurationObjectTPCHWClusterDecoder.C'("param", "uri", runMin, runMax)'
 * </pre>
 *
 * Create an OCDB entry with a TObjString containing param for
 * TPCHWClusterDecoder
 *
 * Parameters: <br>
 * - param (opt)    string to be stored in the TObjSting, default empty
 * - uri   (opt)    the OCDB URI, default $ALICE_ROOT   
 * - runMin (opt)   default 0
 * - runMax (opt)   default 999999999
 * 
 * Current Param : 
 *  - ""  <pre> aliroot -b -q makeConfigurationObjectTPCHWClusterDecoder.C </pre>
 *
 * @author Sergey Gorbunov <sergey.gorbunov@fias.uni-frankfurt.de>
 * @ingroup alihlt_tpc
 */
void makeConfigurationObjectTPCHWClusterDecoder(const Char_t* param="", const Char_t* cdbUri=NULL,
				      Int_t runMin=0, Int_t runMax=AliCDBRunRange::Infinity()) {

  // --------------------------------------
  // -- Setup CDB
  // --------------------------------------

  AliCDBManager* man = AliCDBManager::Instance();
  if (!man) {
    cerr << "Error : Can not get AliCDBManager" << end;
    exit;
  }

  TString storage;
  if (!man->IsDefaultStorageSet()) {
    if ( cdbUri ) {
      storage = cdbUri;
      if ( storage.Contains("://") == 0 ) {
	storage = "local://"; 
	storage += cdbUri;
      }
    } 
    else {
      storage="local://$ALICE_ROOT/OCDB";
    }
    man->SetDefaultStorage(storage);
  } 
  else {
    storage = man->GetDefaultStorage()->GetURI();
  }

  TString path("HLT/ConfigTPC/TPCHWClusterDecoder");

  // --------------------------------------
  // -- Create Config Object
  // --------------------------------------

  // here is the actual content of the configuration object
  TObjString configParam=param;

  TObject *configObj = static_cast<TObject*>(&configParam);
  // --------------------------------------
  // -- Fill Object
  // --------------------------------------
  
  if ( !configObj ) {
    cerr << "Error : No configuration object created" << endl;
    return;
  }
    
  AliCDBPath cdbPath(path);
  AliCDBId   cdbId(cdbPath, runMin, runMax);
  AliCDBMetaData cdbMetaData;
  man->Put(configObj, cdbId, &cdbMetaData);

  printf("Adding %s type OCDB object to %s [%d,%d] in %s \n",
	 configObj->ClassName(), 
	 path.Data(),
	 runMin, runMax, storage.Data());
}

