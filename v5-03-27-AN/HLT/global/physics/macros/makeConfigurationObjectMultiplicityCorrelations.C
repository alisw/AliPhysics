//-*- Mode: C++ -*-
// $Id: makeConfigurationObjectMultiplicityCorrelations.C$
/**
 * @file makeConfigurationObjectMultiplicityCorrelations.C
 * @brief Creation of HLT component configuration objects in OCDB
 *
 * <pre>
 * Usage: aliroot -b -q makeConfigurationObjectMultiplicityCorrelations.C'("param", "uri", runMin, runMax)'
 * </pre>
 *
 * Create an OCDB entry with a TObjString containing param for
 * MultiplicityCorrelations
 *
 * Parameters: <br>
 * - param (opt)    string to be stored in the TObjSting, default empty
 * - uri   (opt)    the OCDB URI, default $ALICE_ROOT   
 * - runMin (opt)   default 0
 * - runMax (opt)   default 999999999
 * - centralityFile (opt) set new centralityFile
 * 
 * Current Param : 
 *  - ""  <pre> aliroot -b -q makeConfigurationObjectMultiplicityCorrelations.C </pre>
 *
 * @author Jochen Thaeder <jochen@thaeder.de>
 * @ingroup alihlt_physics
 */
void makeConfigurationObjectMultiplicityCorrelations(const Char_t* param="-addTrigger CPBI1 -addTrigger CPBI2", const Char_t* cdbUri=NULL,
						     Int_t runMin=0, Int_t runMax=AliCDBRunRange::Infinity(),
						     Char_t *centralityFile="centrality.root") {

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

  if (param) {
    TString path0("HLT/ConfigGlobal/MultiplicityCorrelations");
    // --------------------------------------
    // -- Create Config Object 0
    // --------------------------------------
    
    // here is the actual content of the configuration object
    TObjString configParam=param;
    
    TObject *configObj = static_cast<TObject*>(&configParam);
    // --------------------------------------
    // -- Fill Object 0
    // --------------------------------------
    
    if ( !configObj ) {
      cerr << "Error : No configuration object created" << endl;
      return;
    }
    
    AliCDBPath cdbPath0(path0);
    AliCDBId   cdbId0(cdbPath0, runMin, runMax);
    AliCDBMetaData cdbMetaData0;
    man->Put(configObj, cdbId0, &cdbMetaData0);
    
    printf("Adding %s type OCDB object to %s [%d,%d] in %s \n",
	   configObj->ClassName(), 
	   path0.Data(),
	   runMin, runMax, storage.Data());
  }

  if (centralityFile) {
    TString path1("HLT/ConfigGlobal/MultiplicityCorrelationsCentrality");
    // --------------------------------------
    // -- Create Config Object 1
    // --------------------------------------
   
    // here is the actual content of the configuration object
    TFile *file = TFile::Open(centralityFile);
    TH1F *centrality = static_cast<TH1F*>(file->Get("fHOutMultV0M_percentile")); 
    
    TObject *configObjCentrality = static_cast<TObject*>(centrality);
    // --------------------------------------
    // -- Fill Object
    // --------------------------------------
    
    if ( !configObjCentrality ) {
      cerr << "Error : No centrality configuration object created" << endl;
      return;
    }
    
    AliCDBPath cdbPath1(path1);
    AliCDBId   cdbId1(cdbPath1, runMin, runMax);
    AliCDBMetaData cdbMetaData1;
    man->Put(configObjCentrality, cdbId1, &cdbMetaData1);

    printf("Adding %s type OCDB object to %s [%d,%d] in %s \n",
	   configObjCentrality->ClassName(), 
	   path1.Data(),
	   runMin, runMax, storage.Data());
  }
}

