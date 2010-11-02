//-*- Mode: C++ -*-
// $Id: makeConfigurationObjectdNdPtAnalysis.C$
/**
 * @file makeConfigurationObjectdNdPtAnalysis.C
 * @brief Creation of HLT component configuration objects in OCDB
 *
 * <pre>
 * Usage: aliroot -b -q makeConfigurationObjectdNdPtAnalysis.C'("param", "uri", runMin, runMax)'
 * </pre>
 *
 * Create an OCDB entry with a TObjString containing param for the 
 * ZDC reconstruction.
 *
 * Parameters: <br>
 * - param (opt)    string to be stored in the TObjSting, default empty
 * - uri   (opt)    the OCDB URI, default $ALICE_ROOT   
 * - runMin (opt)   default 0
 * - runMax (opt)   default 999999999
 * 
 * Current Param : 
 *  - ""  <pre> aliroot -b -q makeConfigurationObjectdNdPtAnalysis.C </pre>
 *
 * @author Chiara Oppedisano <Chiara.Oppedisano@to.infn.it>
 *         Jochen Thaeder <jochen@thaeder.de>
 * @ingroup alihlt_zdc
 */
void makeConfigurationObjectdNdPtAnalysis(const Char_t* param="-vertexZRange 20. -meanVertexXYZ 0. 0. 0. -meanVertexXYZSigma 1. 1. 10. -etaRange 0.9 -ptRange 0.15 1.e10 -maxDCAr 3.0 -maxDCAz 30.0 -maxDCAToVertexXY 3.0 -maxDCAToVertexZ 3.0 -requireSigmaToVertex 1", const Char_t* cdbUri=NULL,
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

  TString path("HLT/ConfigAnalysis/dNdPtAnalysis");

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

