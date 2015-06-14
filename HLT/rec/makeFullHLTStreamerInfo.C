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
  TObjArray* listOld = NULL;
  if (pExistingEntry) listOld=(TObjArray*)pExistingEntry->GetObject();
  TObjArray* listNew      = new TObjArray;  

  //create a complete list of streamers for selected class
  TClass* pClass = TClass::GetClass(className);
  if (!pClass) {return;}
  createStreamerInfos(pClass, listNew);

  listNew->Print();

  //update the old list with new definitions
  updateStreamerList(listOld, listNew);

  //listOld->Print();

  // Write the updated object to target OCDB
  AliCDBId id(gkCalibStreamerInfoEntry, firstRun, lastRun, version);
  AliCDBMetaData* metaData = new AliCDBMetaData();
  metaData->SetResponsible("HLT");
  metaData->SetComment("Streamer info for streamed objects in the HLT raw data payload.");
  ocdbTargetStorage->Put(listOld, id, metaData);  
}

void createStreamerInfos(TClass* pClass, TObjArray* listOfStreamerInfos)
{
  //update a TObjArray of stremer infos with TStreamerInfos of the class
  // including TStreamerInfos of all the data member types recursively
  cout << "processing class " << pClass->GetName() << endl;
  if (HaveStreamersForClass(pClass, listOfStreamerInfos)) return;

  TStreamerInfo* streamerInfo = new TStreamerInfo(pClass);
  streamerInfo->Build();
  cout << "adding streamer info for class " << pClass->GetName()
    << "version " << pClass->GetClassVersion()
    << endl;
  listOfStreamerInfos->AddAtFree(streamerInfo);

  TList* listOfDataMembers = pClass->GetListOfDataMembers();
  if (!listOfDataMembers) {return;}
  TIter nextDataMember(listOfDataMembers);
  TDataMember* member=NULL;
  while (member = dynamic_cast<TDataMember*>(nextDataMember()))
  {
    TString memberClassName = member->GetTypeName();
    cout << "  scanning member " << member->GetName()
      << " of type " << memberClassName 
      << endl;
    TClass* memberClass = TClass::GetClass(member->GetTypeName());
    if (!memberClass) {continue;}
    createStreamerInfos(memberClass, listOfStreamerInfos);
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

void updateStreamerList(TObjArray* oldList, TObjArray* newList)
{
  //update the old list with new/unique entries from new list
  for (int j=0; j<oldList->GetEntriesFast(); j++)
  {
    TStreamerInfo* oldStreamerInfo=dynamic_cast<TStreamerInfo*>(oldList->At(j));
    if (!oldStreamerInfo) continue;
    TString oldName = oldStreamerInfo->GetName();
    int oldVersion = oldStreamerInfo->GetClassVersion();

    for (int i=0; i<newList->GetEntriesFast(); i++) {
      TStreamerInfo* newStreamerInfo=dynamic_cast<TStreamerInfo*>(newList->At(i));
      if (!newStreamerInfo) continue;
      TString newName = newStreamerInfo->GetName();
      int newVersion = newStreamerInfo->GetClassVersion();
      //update; or not
      if (! oldName.EqualTo(newName)) continue;
      if (oldVersion == newVersion) continue;

      cout << "adding " << newName << " v" << newVersion
        << " (old: " << oldName << " v" << oldVersion << ")"<<endl;
      oldList->AddAtFree(newStreamerInfo);
    }
  }
  return;
}
