//update the list of streamers contained in the OCDB for a specified class
//do it recursively for all data member types.
//blame: Mikolaj Krzewicki, mkrzewic@cern.ch

void makeFullHLTStreamerInfo(const char* className = "AliESDEvent",
                             const char* ocdbSource = "local:///cvmfs/alice.cern.ch/calibration/data/2015/OCDB",
                             Int_t version = 0,
                             Int_t firstRun = 0,
                             Int_t lastRun = AliCDBRunRange::Infinity(),
                             const char* ocdbTarget = "local://OCDB"
                             )
{
  // Setup the CDB default storage and run number.  
  const char* gkCalibStreamerInfoEntry="HLT/Calib/StreamerInfo";

  // Setup the CDB default storage and run number.
  AliCDBManager* cdbManager = AliCDBManager::Instance();
  if (cdbManager == NULL) {
    cerr << "ERROR: Global CDB manager object does not exist." << endl;
    return;
  }
  AliCDBStorage* ocdbSourceStorage = cdbManager->GetStorage(ocdbSource);
  if (ocdbSourceStorage == NULL) {
    cerr << "ERROR: Could not get storage for: " << cdbPath << endl;
    return;
  }
  AliCDBStorage* ocdbTargetStorage = cdbManager->GetStorage(ocdbTarget);
  if (ocdbSourceStorage == NULL) {
    cerr << "ERROR: Could not get storage for: " << cdbPath << endl;
    return;
  }
  AliCDBRunRange runRange(firstRun, lastRun);
  AliCDBEntry* pExistingEntry=ocdbSourceStorage->Get(gkCalibStreamerInfoEntry, runRange);  
  
  //get the list of streamers from OCDB
  TObjArray* list = NULL;
  if (pExistingEntry) list=(TObjArray*)pExistingEntry->GetObject();
  else list = new TObjArray;  
                             
  //update the list of streamers for selected class
  TClass* pClass = TClass::GetClass(className);
  if (!pClass) {return;}
  updateStreamerInfos(pClass, list);

  // Write the updated object to target OCDB
  AliCDBId id(gkCalibStreamerInfoEntry, firstRun, lastRun, version);
  AliCDBMetaData* metaData = new AliCDBMetaData();
  metaData->SetResponsible("HLT");
  metaData->SetComment("Streamer info for streamed objects in the HLT raw data payload.");
  ocdbTargetStorage->Put(list, id, metaData);  
}

void updateStreamerInfos(TClass* pClass, TObjArray* existingStreamerObjects)
{
  //update a TObjArray of stremer infos with TStreamerInfos of the class
  // including TStreamerInfos of all the data member types recursively
  if (HaveStreamersForClass(pClass, existingStreamerObjects)) {return;}
  
  TStreamerInfo* streamerInfo = new TStreamerInfo(pClass);
  streamerInfo->Build();
  cout << "adding streamer info for class " << pClass->GetName()
    << "version " << pClass->GetClassVersion()
    << endl;
  existingStreamerObjects->AddAtFree(streamerInfo);
  
  TList* listOfDataMembers = pClass->GetListOfDataMembers();
  if (!listOfDataMembers) {return;}
  TIter nextDataMember(listOfDataMembers);
  TDataMember* member=NULL;
  while (member = dynamic_cast<TDataMember*>(nextDataMember()))
  {
    TString memberClassName = member->GetTypeName();
    cout << "  scanning member of type " << memberClassName 
         << " is basic: " << member->IsBasic()
         << " is enum: " << member->IsEnum()
         << endl;
    TClass* memberClass = TClass::GetClass(member->GetTypeName());
    if (!memberClass) {continue;}
    cout << "processing class " << memberClass->GetName() << endl;
    updateStreamerInfos(memberClass, existingStreamerObjects);
  }
}

Bool_t HaveStreamersForClass(TClass* pClass, TObjArray* pInfos)
{
  //check is the class description is already in the list
  TString className = pClass->GetName();
  Version_t classVersion = pClass->GetClassVersion();
  
  TStreamerInfo* pInfo=NULL;
  for (int i=0; i<pInfos->GetEntriesFast(); i++) {
    if (pInfos->At(i)==NULL) continue;
    if (pInfos->At(i)->IsA()!=TStreamerInfo::Class()) {
      cout << "  skipping object " << pInfos->At(i)->GetName() 
           << " class "          << pInfos->At(i)->Class()->GetName() 
           << endl;
      continue;
    }
    pInfo=dynamic_cast<TStreamerInfo*>(pInfos->At(i));
    if (!pInfo) continue;

    TString infoName=pInfo->GetName();
    if (infoName.CompareTo(className)!=0) continue;

    if (pInfo->GetClassVersion()==classVersion) {
      cout << "  nothing to be done, class " 
           << className    << " version " 
           << classVersion << " already in the list" 
           << endl;
      return kTRUE;
    }
  }
  return kFALSE;
}
