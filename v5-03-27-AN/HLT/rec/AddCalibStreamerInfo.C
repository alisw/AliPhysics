// $Id$
/**
 * @file AddCalibStreamerInfo.C
 * @brief Extraction of data blocks from HLTOUT
 *
 * <pre>
 * Usage: aliroot -b -q AddCalibStreamerInfo.C'(classname, uri, version, firstRun, lastRun)'
 *     classname    name of the class to add the streamer info
 *     uri          optional URI of the OCDB, default local://$ALICE_ROOT/OCDB
 *     version      version number of the entry (optional)
 *     firstRun     first run (optional)
 *     lastRun      last run (optional)    
 * </pre>
 *
 * The macro checks whether the streamer info of the current version
 * of the class 'classname' is available in the
 * HLT/Calib/StreamerInfo entry, and adds it if not.
 * 
 * Streamer infos are needed in order to extract TObjects from the HLT raw data
 * payload.
 *
 * @author Matthias.Richter@ift.uib.no
 * @ingroup alihlt_rec
 */
void AddCalibStreamerInfo(
			  const char* className,
			  const char* cdbPath = "local://$ALICE_ROOT/OCDB",
			  Int_t version = 0,
			  Int_t firstRun = 0,
			  Int_t lastRun = AliCDBRunRange::Infinity()
			  )
{
  const char* gkCalibStreamerInfoEntry="HLT/Calib/StreamerInfo";

  // Setup the CDB default storage and run number.
  AliCDBManager* cdbManager = AliCDBManager::Instance();
  if (cdbManager == NULL) {
    cerr << "ERROR: Global CDB manager object does not exist." << endl;
    return;
  }
  AliCDBStorage* storage = cdbManager->GetStorage(cdbPath);
  if (storage == NULL) {
    cerr << "ERROR: Could not get storage for: " << cdbPath << endl;
    return;
  }

  TClass* pClass=TClass::GetClass(className);
  if (!pClass) {
    cerr << "ERROR: Can not find TClass object: " << className << endl;
    return;
  }
  Int_t classVersion=pClass->GetClassVersion();

  TObjArray* pInfos=NULL;
  AliCDBRunRange runRange(firstRun, lastRun);
  AliCDBEntry* pExistingEntry=storage->Get(gkCalibStreamerInfoEntry, runRange);
  if (pExistingEntry) pInfos=(TObjArray*)pExistingEntry->GetObject();
  else pInfos=new TObjArray;

  TStreamerInfo* pInfo=NULL;
  for (int i=0; i<pInfos->GetEntriesFast(); i++) {
    if (pInfos->At(i)==NULL) continue;

    if (pInfos->At(i)->IsA()!=TStreamerInfo::Class()) {
      cout << "skipping object " << pInfos->At(i)->GetName() << " class " << pInfos->At(i)->Class()->GetName() << endl;
      continue;
    }

    pInfo=(TStreamerInfo*)pInfos->At(i);
    TString infoname=pInfo->GetName();
    if (infoname.CompareTo(className)!=0) continue;

    if (pInfo->GetClassVersion()==classVersion) {
      cout << "nothing to be done, class " << className << " version " << classVersion << " already in the object" << endl;
      return;
    }
  }
  pInfo=new TStreamerInfo(pClass);
  pInfo->Build();
  pInfos->AddAtFree(pInfo);
  
  ///////////////////////////////////////////////////////////////////////////////////////////	
  // Write the updated object to OCDB
  AliCDBId id(gkCalibStreamerInfoEntry, firstRun, lastRun, version);
  AliCDBMetaData* metaData = new AliCDBMetaData();
  metaData->SetResponsible("HLT");
  metaData->SetComment("Streamer info for streamed objects in the HLT raw data payload.");
  storage->Put(pInfos, id, metaData);
}
