//-*- Mode: C++ -*-
// $Id: makeTriggerConfigurationObject.C  $
/**
 * @file makeTriggerConfigurationObject.C
 * @brief Creation of HLT component configuration objects in OCDB
 *
 * <pre>
 * Usage: aliroot -b -q makeTriggerConfigurationObject.C'("triggerName", "cdbUri", runMin, runMax)'
 * </pre>
 *
 * Create an CDB entry for a certain trigger.
 * 
 * According to <Trigger-Identifier>, <Minor Version> and <Major Versions>, 
 * CDB objects are created and saved in the CDB. This objects can be any TObject
 * and have to implemented according to the trigger needs.
 *
 * Note : <br>
 * In the path in the CDB '-' is replaced by '_._'  
 *
 * Parameters: <br>
 * - triggerName    Trigger Name following the standard
 *                  H-<Trigger-Identifier>-VXXXX.YYY
                    - VXXXX being the major version number, specifying the settings 
                    - YYY   being the minor version number, specifying the alogrithm version
 * - cdbUri (opt)   the CDB URI, default $ALICE_ROOT   
 * - runMin (opt)   default 0
 * - runMax (opt)   default 999999999
 *
 * Implemented Trigger : <br>
 *  - Barrel_pT_Single  -> AliHLTESDTrackCuts object
 *
 * Usage Example : <br>
 *  aliroot -b -l -q makeTriggerConfigurationObject.C'("H-Barrel_pT_Single-V0001.001")'
 *
 * @author Jochen Thaeder <jochen@thaeder.de>
 * @ingroup alihlt_trigger
 */

// #################################################################################
void makeTriggerConfigurationObject(const Char_t* triggerName, const Char_t* cdbUri=NULL, 
				    Int_t runMin=0, Int_t runMax=AliCDBRunRange::Infinity() ) {

  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libHLTbase.so");
  gSystem->Load("libAliHLTUtil.so");

  // --------------------------------------
  // -- Parse Trigger Name
  // --------------------------------------

  TString name(triggerName);
  if ( !name.BeginsWith("H-") ) {
    cerr << "Error : " << name.Data() << " is not a valid trigger name!" << endl;
    return;
  }
  
  TObjArray* entries = name.Tokenize("-");
  if (entries->GetEntriesFast() != 3 ) {
    cerr << "Error : " << name.Data() << " is not a valid trigger name!" << endl;
    return;
  }

  TString id = (static_cast<TObjString*>(entries->UncheckedAt(1)))->GetString();

  TObjArray* version = (static_cast<TObjString*>(entries->UncheckedAt(2)))->GetString().Tokenize(".");
  if (version->GetEntriesFast() != 2 ) {
    cerr << "Error : " << name.Data() << " is not a valid trigger name!" << endl;
    return;
  }

  // -- Major version ID = Settings
  TString major = (static_cast<TObjString*>(version->UncheckedAt(0)))->GetString();

  // -- Minor version ID = Algorithm
  TString minor = (static_cast<TObjString*>(version->UncheckedAt(1)))->GetString();
  

  // -------------------------------------------------------------------------------------

  printf ("Trigger Setup for : %s \n", triggerName);

  // --------------------------------------
  // -- Create Configuration Object
  // --------------------------------------
  
  TObject* configObj = NULL;

  if ( !id.CompareTo("Barrel_pT_Single") ) { // ----------------------------------------
    
    AliHLTESDTrackCuts* esdTrackCuts = NULL;
    TString description;

    if ( !minor.CompareTo("001") ) {
      printf (" - Trigger Setup for 2010 pp data.\n");
      
      esdTrackCuts = AliHLTESDTrackCuts::GetStandardTrackCuts2010pp();
      description += esdTrackCuts->GetTitle();
      
      // -- pt cut 1 GeV - minClusters 80
      if ( !major.CompareTo("V0001") ) {
	esdTrackCuts->SetPtRange(1.);
	esdTrackCuts->SetMinNClustersTPC(80);
	description += " && p_t > 1 GeV/c && clustersTPC > 80";
      }
      // -- pt cut 2 GeV - minClusters 80
      else if ( !major.CompareTo("V0002") ) {
	esdTrackCuts->SetPtRange(2.);
	esdTrackCuts->SetMinNClustersTPC(80);
	description += " && p_t > 2 GeV/c && clustersTPC > 80";
      }
      // -- pt cut 3 GeV - minClusters 100
      else if ( !major.CompareTo("V0003") ) {
	esdTrackCuts->SetPtRange(3.);
	esdTrackCuts->SetMinNClustersTPC(100);
	description += " && p_t > 3 GeV/c && clustersTPC > 100";
      }
      else {
	cerr << "Error : Major version " << major.Data() << " is not implemented!" << endl;
	return;
      }
    }
    else {
      cerr << "Error : Minor version " << minor.Data() << " is not implemented!" << endl;
      return;
    }
    
    if ( !esdTrackCuts ) {
      cerr << "Error : No AliHLTESDTrackCuts object created" << endl;
      return;
    }
    
    esdTrackCuts->SetTitle(description);
    configObj = static_cast<TObject*>(esdTrackCuts);
  }
  else { // ------------------------------------------------------------------------------
    cerr << "Error : Trigger name " << id.Data() << " is not implemented!" << endl;
      return;
  }

  // -------------------------------------------------------------------------------------  
  printf(" - TrackCuts : %s\n", esdTrackCuts->GetTitle());

  // --------------------------------------
  // -- Setup CDB
  // --------------------------------------

  AliCDBManager* man = AliCDBManager::Instance();
  if (!man) {
    cerr << "Error : Can not get AliCDBManager" << endl;
    return;
  }
  
  TString storage;
  if ( !man->IsDefaultStorageSet() ) {
    if ( cdbUri ) {
      storage = cdbUri;
      if ( storage.Contains("://") == 0 ) {
	storage = "local://"; 
	storage += cdbUri;
      }
    } 
    else {
      storage = "local://$ALICE_ROOT/OCDB";
    }
    man->SetDefaultStorage(storage);
  } 
  else {
    storage = man->GetDefaultStorage()->GetURI();
  }

  TString path("HLT/ConfigHLT/");
  path += name.ReplaceAll("-",1,"_._",3);

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
  
  // --------------------------------------
  // -- Clean up
  // --------------------------------------
  if (version)
    delete version;
  version = NULL;

  if (entries)
    delete entries;
  entries = NULL;

  if (esdTrackCuts)
    delete esdTrackCuts;
  esdTrackCuts = NULL;
}
// #################################################################################
void makeTriggerConfigurationObject() {

  cout << "===============================================================" << endl;
  cout << "usage: aliroot -b -q -l makeComponentConfigurationObject.C'(\"triggerName\", \"cdbUri\", rangeMin, rangeMax)'" << endl << endl;
  cout << "  triggerName    Trigger Name following the standard" << endl;
  cout << "                 H-<Trigger-Identifier>-VXXXX.YYY" << endl;
  cout << "                    - VXXXX being the major version number, specifying the settings" << endl;
  cout << "                    - YYY   being the minor version number, specifying the alogrithm version" << endl;
  cout << "  cdbUri   (opt) the OCDB URI, default $ALICE_ROOT/OCDB   " << endl;
  cout << "  rangeMin (opt) default 0" << endl;
  cout << "  rangeMax (opt) default 999999999" << endl << endl;
  cout << "example: aliroot -b -l -q makeTriggerConfigurationObject.C'(\"H-Barrel_pT_Single-V0001.001\")'" << endl;
  cout << "===============================================================" << endl;
}
