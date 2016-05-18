/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
//-------------------------------------------------------------------------
//   Implementation of AliCDBManager and AliCDBParam classe
//   Author: Alberto Colla
//   e-mail: Alberto.Colla@cern.ch
//-------------------------------------------------------------------------

#include <fstream>

#include "AliCDBManager.h"
#include "AliCDBStorage.h"
#include "AliLog.h"
#include "AliCDBDump.h"
#include "AliCDBLocal.h"
#include "AliCDBGrid.h"
#include "AliCDBEntry.h"
#include "AliCDBHandler.h"

#include <TObjString.h>
#include <TSAXParser.h>
#include <TFile.h>
#include <TKey.h>
#include <TUUID.h>
#include <TGrid.h>
#include "TMessage.h"
#include "TObject.h"
#include "TRegexp.h"

ClassImp(AliCDBParam)

ClassImp(AliCDBManager)

//TODO OCDB and Reference folder should not be fully hardcoded but built from run number (or year/LHC period)
TString AliCDBManager::fgkCondUri("alien://folder=/alice/cern.ch/user/a/aliprod/testCDB/CDB?user=aliprod");
TString AliCDBManager::fgkRefUri("alien://folder=/alice/cern.ch/user/a/aliprod/testCDB/Reference?user=aliprod");
TString AliCDBManager::fgkMCIdealStorage("alien://folder=/alice/simulation/2008/v4-15-Release/Ideal");
TString AliCDBManager::fgkMCFullStorage("alien://folder=/alice/simulation/2008/v4-15-Release/Full");
TString AliCDBManager::fgkMCResidualStorage("alien://folder=/alice/simulation/2008/v4-15-Release/Residual");
TString AliCDBManager::fgkOCDBFolderXMLfile("alien:///alice/data/OCDBFoldervsRunRange.xml");
AliCDBManager* AliCDBManager::fgInstance = 0x0;

//_____________________________________________________________________________
AliCDBManager* AliCDBManager::Instance(TMap *entryCache, Int_t run) {
// returns AliCDBManager instance (singleton)

  if (!fgInstance) {
    fgInstance = new AliCDBManager();
    if (!entryCache)
      fgInstance->Init();
    else
      fgInstance->InitFromCache(entryCache,run);
  }

  return fgInstance;
}

//_____________________________________________________________________________
void AliCDBManager::Init() {
// factory registering

  RegisterFactory(new AliCDBDumpFactory());
  RegisterFactory(new AliCDBLocalFactory());
  // AliCDBGridFactory is registered only if AliEn libraries are enabled in Root
  if(!gSystem->Exec("root-config --has-alien 2>/dev/null |grep yes 2>&1 > /dev/null")){ // returns 0 if yes
    AliInfo("AliEn classes enabled in Root. AliCDBGrid factory registered.");
    RegisterFactory(new AliCDBGridFactory());
    fCondParam = CreateParameter(fgkCondUri);
    fRefParam = CreateParameter(fgkRefUri);
  }

  InitShortLived();
}

//_____________________________________________________________________________
void AliCDBManager::InitFromCache(TMap *entryCache, Int_t run) {
// initialize manager from existing cache
// used on the slaves in case of parallel reconstruction
  SetRun(run);

  TIter iter(entryCache->GetTable());
  TPair* pair = 0;

  while((pair = dynamic_cast<TPair*> (iter.Next()))){
    fEntryCache.Add(pair->Key(),pair->Value());
  }
  // fEntry is the new owner of the cache
  fEntryCache.SetOwnerKeyValue(kTRUE,kTRUE);
  entryCache->SetOwnerKeyValue(kFALSE,kFALSE);
  AliInfo(Form("%d cache entries have been loaded",fEntryCache.GetEntries()));
}

//_____________________________________________________________________________
void  AliCDBManager::DumpToSnapshotFile(const char* snapshotFileName, Bool_t singleKeys) const {
//
// If singleKeys is true, dump the entries map and the ids list to the snapshot file
// (provided mostly for historical reasons, the file is then read with InitFromSnapshot),
// otherwise write to file each AliCDBEntry separately (the is the preferred way, the file
// is then read with SetSnapshotMode).

  // open the file
  TFile *f = TFile::Open(snapshotFileName,"RECREATE");
  if (!f || f->IsZombie()){
    AliError(Form("Cannot open file %s",snapshotFileName));
    return;
  }

  AliInfo(Form("Dumping entriesMap (entries'cache) with %d entries!\n", fEntryCache.GetEntries()));
  AliInfo(Form("Dumping entriesList with %d entries!\n", fIds->GetEntries()));

  f->cd();
  if(singleKeys){
    f->WriteObject(&fEntryCache,"CDBentriesMap");
    f->WriteObject(fIds,"CDBidsList");
  }else{
    // We write the entries one by one named by their calibration path
    TIter iter(fEntryCache.GetTable());
    TPair* pair = 0;
    while((pair = dynamic_cast<TPair*> (iter.Next()))){
      TObjString *os = dynamic_cast<TObjString*>(pair->Key());
      if (!os) continue;
      TString path = os->GetString();
      AliCDBEntry *entry = dynamic_cast<AliCDBEntry*>(pair->Value());
      if (!entry) continue;
      path.ReplaceAll("/","*");
      entry->Write(path.Data());
    }
  }
  f->Close();
  delete f;
}

//_____________________________________________________________________________
void  AliCDBManager::DumpToLightSnapshotFile(const char* lightSnapshotFileName) const {
// The light snapshot does not contain the CDB objects (AliCDBEntries) but
// only the information identifying them, that is the map of storages and
// the list of AliCDBIds, as in the UserInfo of AliESDs.root

  // open the file
  TFile *f = TFile::Open(lightSnapshotFileName,"RECREATE");
  if (!f || f->IsZombie()){
    AliError(Form("Cannot open file %s",lightSnapshotFileName));
    return;
  }

  AliInfo(Form("Dumping map of storages with %d entries!\n", fStorageMap->GetEntries()));
  AliInfo(Form("Dumping entriesList with %d entries!\n", fIds->GetEntries()));
  f->WriteObject(fStorageMap,"cdbStoragesMap");
  f->WriteObject(fIds,"CDBidsList");

  f->Close();
  delete f;
}

//_____________________________________________________________________________
Bool_t AliCDBManager::InitFromSnapshot(const char* snapshotFileName, Bool_t overwrite){
// initialize manager from a CDB snapshot, that is add the entries
// to the entries map and the ids to the ids list taking them from
// the map and the list found in the input file

  // if the manager is locked it cannot initialize from a snapshot
  if(fLock) {
    AliError("Being locked I cannot initialize from the snapshot!");
    return kFALSE;
  }	

  // open the file
  TString snapshotFile(snapshotFileName);
  if(snapshotFile.BeginsWith("alien://")){
    if(!gGrid) {
      TGrid::Connect("alien://","");
      if(!gGrid) {
        AliError("Connection to alien failed!");
        return kFALSE;
      }
    }
  }

  TFile *f = TFile::Open(snapshotFileName);
  if (!f || f->IsZombie()){
    AliError(Form("Cannot open file %s",snapshotFileName));
    return kFALSE;
  }

  // retrieve entries' map from snapshot file
  TMap *entriesMap = 0;
  TIter next(f->GetListOfKeys());
  TKey *key;
  while ((key = (TKey*)next())) {
    if (strcmp(key->GetClassName(),"TMap") != 0) continue;
    entriesMap = (TMap*)key->ReadObj();
    break;
  }
  if (!entriesMap || entriesMap->GetEntries()==0){
    AliError("Cannot get valid map of CDB entries from snapshot file");
    return kFALSE;
  }

  // retrieve ids' list from snapshot file
  TList *idsList = 0;
  TIter nextKey(f->GetListOfKeys());
  TKey *keyN;
  while ((keyN = (TKey*)nextKey())) {
    if (strcmp(keyN->GetClassName(),"TList") != 0) continue;
    idsList = (TList*)keyN->ReadObj();
    break;
  }
  if (!idsList || idsList->GetEntries()==0){
    AliError("Cannot get valid list of CDB entries from snapshot file");
    return kFALSE;
  }

  // Add each (entry,id) from the snapshot to the memory: entry to the cache, id to the list of ids.
  // If "overwrite" is false: add the entry to the cache and its id to the list of ids
  // only if neither of them is already there.
  // If "overwrite" is true: write the snapshot entry,id in any case. If something
  // was already there for that calibration type, remove it and issue a warning
  TIter iterObj(entriesMap->GetTable());
  TPair* pair = 0;
  Int_t nAdded=0;
  while((pair = dynamic_cast<TPair*> (iterObj.Next()))){
    TObjString* os = (TObjString*) pair->Key();
    TString path = os->GetString();
    TIter iterId(idsList);
    AliCDBId* id=0;
    AliCDBId* correspondingId=0;
    while((id = dynamic_cast<AliCDBId*> (iterId.Next()))){
      TString idpath(id->GetPath());
      if(idpath==path){
        correspondingId=id;
        break;
      }
    }
    if(!correspondingId){
      AliError(Form("id for \"%s\" not found in the snapshot (while entry was). This entry is skipped!",path.Data()));
      break;
    }
    Bool_t cached = fEntryCache.Contains(path.Data());
    Bool_t registeredId = kFALSE;
    TIter iter(fIds);
    AliCDBId *idT = 0;
    while((idT = dynamic_cast<AliCDBId*> (iter.Next()))){
      if(idT->GetPath()==path){
        registeredId = kTRUE;
        break;
      }
    }

    if(overwrite){
      if(cached || registeredId){
        AliWarning(Form("An entry was already cached for \"%s\". Removing it before caching from snapshot",path.Data()));
        UnloadFromCache(path.Data());
      }
      fEntryCache.Add(pair->Key(),pair->Value());
      fIds->Add(id);
      nAdded++;
    }else{
      if(cached || registeredId){
        AliWarning(Form("An entry was already cached for \"%s\". Not adding this object from snapshot",path.Data()));
      }else{
        fEntryCache.Add(pair->Key(),pair->Value());
        fIds->Add(id);
        nAdded++;
      }
    }
  }

  // fEntry is the new owner of the cache
  fEntryCache.SetOwnerKeyValue(kTRUE,kTRUE);
  entriesMap->SetOwnerKeyValue(kFALSE,kFALSE);
  fIds->SetOwner(kTRUE);
  idsList->SetOwner(kFALSE);
  AliInfo(Form("%d new (entry,id) cached. Total number %d",nAdded,fEntryCache.GetEntries()));

  f->Close();
  delete f;

  return kTRUE;
}

//_____________________________________________________________________________
void AliCDBManager::Destroy() {
// delete ALCDBManager instance and active storages

  if (fgInstance) {
    //fgInstance->Delete();
    delete fgInstance;
    fgInstance = 0x0;
  }
}

//_____________________________________________________________________________
AliCDBManager::AliCDBManager():
  TObject(),
  fFactories(),
  fActiveStorages(),
  fSpecificStorages(),
  fEntryCache(),
  fPromptEntryCache(),
  fIds(0),
  fStorageMap(0),
  fShortLived(0),
  fDefaultStorage(NULL),
  fDrainStorage(NULL),
  fCondParam(0),
  fRefParam(0),
  fRun(-1),
  fCache(kTRUE),
  fPromptCache(kFALSE),
  fLock(kFALSE),
  fSnapshotMode(kFALSE),
  fSnapshotFile(0),
  fOCDBUploadMode(kFALSE),
  fRaw(kFALSE),
  fCvmfsOcdb(""),
  fStartRunLHCPeriod(-1),
  fEndRunLHCPeriod(-1),
  fLHCPeriod(""),
  fKey(0)
{
  // default constuctor
  fFactories.SetOwner(1);
  fActiveStorages.SetOwner(1);
  fSpecificStorages.SetOwner(1);
  fEntryCache.SetName("CDBEntryCache");
  fEntryCache.SetOwnerKeyValue(kTRUE,kTRUE);
  fPromptEntryCache.SetOwnerKeyValue(kTRUE,kTRUE);

  fStorageMap = new TMap();
  fStorageMap->SetOwner(1);
  fIds = new TList();
  fIds->SetOwner(1);
}

//_____________________________________________________________________________
AliCDBManager::~AliCDBManager() {
// destructor
  ClearCache();
  ClearPromptCache();
  DestroyActiveStorages();
  fFactories.Delete();
  fDrainStorage = 0x0;
  fDefaultStorage = 0x0;
  delete fStorageMap; fStorageMap = 0;
  delete fIds; fIds = 0;
  delete fCondParam;
  delete fRefParam;
  delete fShortLived; fShortLived = 0x0;
  //fSnapshotCache = 0;
  //fSnapshotIdsList = 0;
  if(fSnapshotMode){
    fSnapshotFile->Close();
    fSnapshotFile = 0;
  }
}

//_____________________________________________________________________________
void AliCDBManager::PutActiveStorage(AliCDBParam* param, AliCDBStorage* storage){
// put a storage object into the list of active storages

  fActiveStorages.Add(param, storage);
  AliDebug(1, Form("Active storages: %d", fActiveStorages.GetEntries()));
}

//_____________________________________________________________________________
void AliCDBManager::RegisterFactory(AliCDBStorageFactory* factory) {
// add a storage factory to the list of registerd factories

  if (!fFactories.Contains(factory)) {
    fFactories.Add(factory);
  }
}

//_____________________________________________________________________________
Bool_t AliCDBManager::HasStorage(const char* dbString) const {
// check if dbString is a URI valid for one of the registered factories

  TIter iter(&fFactories);

  AliCDBStorageFactory* factory=0;
  while ((factory = (AliCDBStorageFactory*) iter.Next())) {

    if (factory->Validate(dbString)) {
      return kTRUE;
    }
  }

  return kFALSE;
}

//_____________________________________________________________________________
AliCDBParam* AliCDBManager::CreateParameter(const char* dbString) const {
// create AliCDBParam object from URI string

  TString uriString(dbString);

  if ( !fCvmfsOcdb.IsNull() && uriString.BeginsWith("alien://")) {
    AlienToCvmfsUri(uriString);
  }

  TIter iter(&fFactories);

  AliCDBStorageFactory* factory=0;
  while ((factory = (AliCDBStorageFactory*) iter.Next())) {
    AliCDBParam* param = factory->CreateParameter(uriString);
    if(param) return param;
  }

  return NULL;
}

//_____________________________________________________________________________
void AliCDBManager::AlienToCvmfsUri(TString& uriString) const {
// convert alien storage uri to local:///cvmfs storage uri (called when OCDB_PATH is set)

  TObjArray *arr = uriString.Tokenize('?');
  TIter iter(arr);
  TObjString *str = 0;
  TString entryKey = "";
  TString entryValue  = "";
  TString newUriString = "";
  while((str = (TObjString*) iter.Next())){
    TString entry(str->String());
    Int_t indeq = entry.Index('=');
    entryKey = entry(0, indeq+1);
    entryValue = entry(indeq+1, entry.Length()-indeq);

    if ( entryKey.Contains("folder", TString::kIgnoreCase) )
    {

      TRegexp re_RawFolder("^/alice/data/20[0-9]+/OCDB");
      TRegexp re_MCFolder("^/alice/simulation/2008/v4-15-Release");
      TString rawFolder = entryValue(re_RawFolder);
      TString mcFolder = entryValue(re_MCFolder);
      if ( !rawFolder.IsNull() ){
        entryValue.Replace(0, 6, "/cvmfs/alice-ocdb.cern.ch/calibration");
        //entryValue.Replace(entryValue.Length()-4, entryValue.Length(), "");
      } else if ( !mcFolder.IsNull() ){
        entryValue.Replace(0,36,"/cvmfs/alice-ocdb.cern.ch/calibration/MC");
      } else {
        AliFatal(Form("Environment variable for cvmfs OCDB folder set for an invalid OCDB storage:\n   %s", entryValue.Data()));
      }
    } else {
      newUriString += entryKey;
    }
    newUriString += entryValue;
    newUriString += '?';
  }
  newUriString.Prepend("local://");
  newUriString.Remove(TString::kTrailing, '?');
  uriString = newUriString;
}

//_____________________________________________________________________________
AliCDBStorage* AliCDBManager::GetStorage(const char* dbString) {
// Get the CDB storage corresponding to the URI string passed as argument
// If "raw://" is passed, get the storage for the raw OCDB for the current run (fRun)

  TString uriString(dbString);
  if (uriString.EqualTo("raw://")) {
    if (!fLHCPeriod.IsNull() && !fLHCPeriod.IsWhitespace()) {
      return GetDefaultStorage();
    } else {
      TString lhcPeriod("");
      Int_t startRun = -1, endRun = -1;
      GetLHCPeriodAgainstAlienFile(fRun, lhcPeriod, startRun, endRun);
      return GetStorage(lhcPeriod.Data());
    }
  }

  AliCDBParam* param = CreateParameter(dbString);
  if (!param) {
    AliError(Form("Failed to activate requested storage! Check URI: %s", dbString));
    return NULL;
  }

  AliCDBStorage* aStorage = GetStorage(param);

  delete param;
  return aStorage;
}

//_____________________________________________________________________________
AliCDBStorage* AliCDBManager::GetStorage(const AliCDBParam* param) {
// get storage object from AliCDBParam object

  // if the list of active storages already contains
  // the requested storage, return it
  AliCDBStorage* aStorage = GetActiveStorage(param);
  if (aStorage) {
    return aStorage;
  }

  // if lock is ON, cannot activate more storages!
  if(fLock) {
    if (fDefaultStorage) {
      AliFatal("Lock is ON, and default storage is already set: "
          "cannot reset it or activate more storages!");
    }
  }	

  // loop on the list of registered factories
  TIter iter(&fFactories);
  AliCDBStorageFactory* factory=0;
  while ((factory = (AliCDBStorageFactory*) iter.Next())) {

    // each factory tries to create its storage from the parameter
    aStorage = factory->Create(param);
    if (aStorage) {
      PutActiveStorage(param->CloneParam(), aStorage);
      aStorage->SetURI(param->GetURI());
      if(fRun >= 0) {
        if( aStorage->GetType() == "alien" || aStorage->GetType() == "local" )
          aStorage->QueryCDB(fRun);
      }
      return aStorage;
    }
  }

  AliError(Form("Failed to activate requested storage! Check URI: %s", param->GetURI().Data()));

  return NULL;
}

//_____________________________________________________________________________
AliCDBStorage* AliCDBManager::GetActiveStorage(const AliCDBParam* param) {
// get a storage object from the list of active storages

  return dynamic_cast<AliCDBStorage*> (fActiveStorages.GetValue(param));
}

//_____________________________________________________________________________
TList* AliCDBManager::GetActiveStorages() {
// return list of active storages
// user has responsibility to delete returned object

  TList* result = new TList();

  TIter iter(fActiveStorages.GetTable());
  TPair* aPair=0;
  while ((aPair = (TPair*) iter.Next())) {
    result->Add(aPair->Value());
  }

  return result;
}

//_____________________________________________________________________________
void AliCDBManager::SetDrain(const char* dbString) {
// set drain storage from URI string

  fDrainStorage = GetStorage(dbString);	
}

//_____________________________________________________________________________
void AliCDBManager::SetDrain(const AliCDBParam* param) {
// set drain storage from AliCDBParam

  fDrainStorage = GetStorage(param);
}

//_____________________________________________________________________________
void AliCDBManager::SetDrain(AliCDBStorage* storage) {
// set drain storage from another active storage

  fDrainStorage = storage;
}

//_____________________________________________________________________________
Bool_t AliCDBManager::Drain(AliCDBEntry *entry) {
// drain retrieved object to drain storage

  AliDebug(2, "Draining into drain storage...");
  return fDrainStorage->Put(entry);
}

//____________________________________________________________________________
Bool_t AliCDBManager::SetOCDBUploadMode() {
// Set the framework in official upload mode. This tells the framework to upload
// objects to cvmfs after they have been uploaded to AliEn OCDBs.
// It return false if the executable to upload to cvmfs is not found.

  TString cvmfsUploadExecutable("$HOME/bin/ocdb-cvmfs");
  gSystem->ExpandPathName(cvmfsUploadExecutable);
  if ( gSystem->AccessPathName(cvmfsUploadExecutable) )
    return kFALSE;
  fOCDBUploadMode = kTRUE;
  return kTRUE;
}

//____________________________________________________________________________
void AliCDBManager::SetDefaultStorage(const char* storageUri) {
// sets default storage from URI string

  // if in the cvmfs case (triggered by environment variable) check for path validity
  // and modify Uri if it is "raw://"
  TString cvmfsOcdb(gSystem->Getenv("OCDB_PATH"));
  if (! cvmfsOcdb.IsNull()){
    fCvmfsOcdb = cvmfsOcdb;
    ValidateCvmfsCase();
  }

  // checking whether we are in the raw case
  TString uriTemp(storageUri);
  if (uriTemp == "raw://") {
    fRaw = kTRUE; // read then by SetRun to check if the method has to be called again with expanded uri
    AliInfo("Setting the run-number will set the corresponding OCDB for raw data reconstruction.");
    return;
  }

  AliCDBStorage* bckStorage = fDefaultStorage;

  fDefaultStorage = GetStorage(storageUri);

  if(!fDefaultStorage) return;

  if(bckStorage && (fDefaultStorage != bckStorage)){
    AliWarning("Existing default storage replaced: clearing cache!");
    ClearCache();
  }

  if (fStorageMap->Contains("default")) {
    delete fStorageMap->Remove(((TPair*)fStorageMap->FindObject("default"))->Key());
  }
  fStorageMap->Add(new TObjString("default"), new TObjString(fDefaultStorage->GetURI()));
}

//_____________________________________________________________________________
void AliCDBManager::SetDefaultStorage(const AliCDBParam* param) {
// set default storage from AliCDBParam object

  AliCDBStorage* bckStorage = fDefaultStorage;

  fDefaultStorage = GetStorage(param);

  if(!fDefaultStorage) return;

  if(bckStorage && (fDefaultStorage != bckStorage)){
    AliWarning("Existing default storage replaced: clearing cache!");
    ClearCache();
  }

  if (fStorageMap->Contains("default")) {
    delete fStorageMap->Remove(((TPair*)fStorageMap->FindObject("default"))->Key());
  }
  fStorageMap->Add(new TObjString("default"), new TObjString(fDefaultStorage->GetURI()));
}

//_____________________________________________________________________________
void AliCDBManager::SetDefaultStorage(AliCDBStorage* storage) {
// set default storage from another active storage

  // if lock is ON, cannot activate more storages!
  if(fLock) {
    if (fDefaultStorage) {
      AliFatal("Lock is ON, and default storage is already set: "
          "cannot reset it or activate more storages!");
    }
  }	

  if (!storage) {
    UnsetDefaultStorage();
    return;
  }

  AliCDBStorage* bckStorage = fDefaultStorage;

  fDefaultStorage = storage;

  if(bckStorage && (fDefaultStorage != bckStorage)){
    AliWarning("Existing default storage replaced: clearing cache!");
    ClearCache();
  }

  if (fStorageMap->Contains("default")) {
    delete fStorageMap->Remove(((TPair*)fStorageMap->FindObject("default"))->Key());
  }
  fStorageMap->Add(new TObjString("default"), new TObjString(fDefaultStorage->GetURI()));
}

//_____________________________________________________________________________
void AliCDBManager::SetDefaultStorage(const char* mcString, const char* simType) {
// sets default storage for MC data
// mcString MUST be "MC",
// simType can be "Ideal","Residual","Full"

  TString strmcString(mcString);
  TString strsimType(simType);
  TString dbString;
  if (strmcString != "MC"){
    AliFatal("Method requires first string to be MC!");
  }
  else {
    if (strsimType == "Ideal"){
      dbString = fgkMCIdealStorage;
    }
    else if (strsimType == "Full"){
      dbString = fgkMCFullStorage;
    }
    else if (strsimType == "Residual"){
      dbString = fgkMCResidualStorage;
    }
    else {
      AliFatal("Error in setting the storage for MC data, second argument MUST be either \"Ideal\" or \"Full\" or \"Residual\".");
    }

    SetDefaultStorage(dbString.Data());
    fStartRunLHCPeriod=0;
    fEndRunLHCPeriod=AliCDBRunRange::Infinity();
    if(!fDefaultStorage) AliFatal(Form("%s storage not there! Please check!",dbString.Data()));
  }
}

//_____________________________________________________________________________
void AliCDBManager::ValidateCvmfsCase() const {
  // The OCDB_PATH variable contains the path to the directory in /cvmfs/ which is
  // an AliRoot tag based snapshot of the AliEn file catalogue (e.g. 
  // /cvmfs/alice.cern.ch/x86_64-2.6-gnu-4.1.2/Packages/OCDB/v5-05-76-AN).
  // The directory has to contain:
  // 1) <data|MC>/20??.list.gz gzipped text files listing the OCDB files (seen by that AliRoot tag)
  // 2) bin/getOCDBFilesPerRun.sh   (shell+awk) script extracting from 1) the list
  //    of valid files for the given run.

    if (! fCvmfsOcdb.BeginsWith("/cvmfs"))  //!!!! to be commented out for testing
      AliFatal(Form("OCDB_PATH set to an invalid path: %s", fCvmfsOcdb.Data()));

    TString cvmfsUri(fCvmfsOcdb);
    gSystem->ExpandPathName(cvmfsUri);
    if (gSystem->AccessPathName(cvmfsUri))
      AliFatal(Form("OCDB_PATH set to an invalid path: %s", cvmfsUri.Data()));
    
    // check that we find the two scripts we need

    AliDebug(3, "OCDB_PATH envvar is set. Changing OCDB storage from alien:// to local:///cvmfs type.");
    cvmfsUri = cvmfsUri.Strip(TString::kTrailing, '/');
    cvmfsUri.Append("/bin/getOCDBFilesPerRun.sh");
    if (gSystem->AccessPathName(cvmfsUri))
      AliFatal(Form("Cannot find valid script: %s", cvmfsUri.Data()));
}

//_____________________________________________________________________________
void AliCDBManager::SetDefaultStorageFromRun(Int_t run) {
// set default storage from the run number - to be used only with raw data	

  // if lock is ON, cannot activate more storages!
  if(fLock) {
    if (fDefaultStorage) {
      AliFatal("Lock is ON, and default storage is already set: "
          "cannot activate default storage from run number");
    }
  }	

  TString lhcPeriod("");
  Int_t startRun = 0, endRun = 0;
  if (! fCvmfsOcdb.IsNull()) { // fRaw and cvmfs case: set LHC period from cvmfs file
    GetLHCPeriodAgainstCvmfsFile(run, lhcPeriod, startRun, endRun);
  } else {                     // fRaw: set LHC period from AliEn XML file
    GetLHCPeriodAgainstAlienFile(run, lhcPeriod, startRun, endRun);
  }

  fLHCPeriod = lhcPeriod;
  fStartRunLHCPeriod = startRun;
  fEndRunLHCPeriod = endRun;

  SetDefaultStorage(fLHCPeriod.Data());
  if(!fDefaultStorage) AliFatal(Form("%s storage not there! Please check!",fLHCPeriod.Data()));

}

//_____________________________________________________________________________
void AliCDBManager::GetLHCPeriodAgainstAlienFile(Int_t run, TString& lhcPeriod, Int_t& startRun, Int_t& endRun) {
// set LHC period (year + first, last run) comparing run number and AliEn XML file
 
// retrieve XML file from alien
  if(!gGrid) {
    TGrid::Connect("alien://","");
    if(!gGrid) {
      AliError("Connection to alien failed!");
      return;
    }
  }
  TUUID uuid;
  TString rndname = gSystem->TempDirectory(); // "/tmp/";
  rndname += "/";
  rndname += "OCDBFolderXML.";
  rndname += uuid.AsString();
  rndname += ".xml";
  AliDebug(2, Form("file to be copied = %s", fgkOCDBFolderXMLfile.Data()));
  if (!TFile::Cp(fgkOCDBFolderXMLfile.Data(), rndname.Data())) {
    AliFatal(Form("Cannot make a local copy of OCDBFolder xml file in %s",rndname.Data()));
  }
  AliCDBHandler* saxcdb = new AliCDBHandler();
  saxcdb->SetRun(run);
  TSAXParser *saxParser = new TSAXParser();
  saxParser->ConnectToHandler("AliCDBHandler", saxcdb);
  saxParser->ParseFile(rndname.Data());
  AliInfo(Form(" LHC folder = %s", saxcdb->GetOCDBFolder().Data()));
  AliInfo(Form(" LHC period start run = %d", saxcdb->GetStartRunRange()));
  AliInfo(Form(" LHC period end run = %d", saxcdb->GetEndRunRange()));
  lhcPeriod = saxcdb->GetOCDBFolder();
  startRun = saxcdb->GetStartRunRange();
  endRun = saxcdb->GetEndRunRange();
}

//_____________________________________________________________________________
void AliCDBManager::GetLHCPeriodAgainstCvmfsFile(Int_t run, TString& lhcPeriod, Int_t& startRun, Int_t& endRun) {
// set LHC period (year + first, last run) comparing run number and CVMFS file
// We don't want to connect to AliEn just to set the uri from the runnumber
// for that we use the script getUriFromYear.sh in the cvmfs AliRoot package

  TString getYearScript(fCvmfsOcdb);
  getYearScript = getYearScript.Strip(TString::kTrailing, '/');
  getYearScript.Append("/bin/getUriFromYear.sh");
  if (gSystem->AccessPathName(getYearScript))
    AliFatal(Form("Cannot find valid script: %s", getYearScript.Data()));
  TString inoutFile(gSystem->WorkingDirectory());
  inoutFile += "/uri_range_";
  inoutFile += TString::Itoa(run,10);
  TString command(getYearScript);
  command += ' ';
  command += TString::Itoa(run,10);
  command += Form(" > %s", inoutFile.Data());
  AliDebug(3, Form("Running command: \"%s\"",command.Data()));
  Int_t result = gSystem->Exec(command.Data());
  if(result != 0) {
    AliFatal(Form("Was not able to execute \"%s\"", command.Data()));
  }

  // now read the file with the uri and first and last run
  std::ifstream file(inoutFile.Data());
  if (!file.is_open()) {
    AliFatal(Form("Error opening file \"%s\"!", inoutFile.Data()));
  }
  TString line;
  TObjArray* oStringsArray = 0;
  while (line.ReadLine(file)){
    oStringsArray = line.Tokenize(' ');
  }
  TObjString *oStrUri = dynamic_cast<TObjString*> (oStringsArray->At(0));
  TObjString *oStrFirst = dynamic_cast<TObjString*> (oStringsArray->At(1));
  TString firstRun = oStrFirst->GetString();
  TObjString *oStrLast = dynamic_cast<TObjString*> (oStringsArray->At(2));
  TString lastRun = oStrLast->GetString();

  lhcPeriod = oStrUri->GetString();
  startRun = firstRun.Atoi();
  endRun = lastRun.Atoi();

  file.close();
}

//_____________________________________________________________________________
void AliCDBManager::UnsetDefaultStorage() {
// Unset default storage

  // if lock is ON, action is forbidden!
  if(fLock) {
    if (fDefaultStorage) {
      AliFatal("Lock is ON: cannot unset default storage!");
    }
  }

  if (fDefaultStorage) {
    AliWarning("Clearing cache!");
    ClearCache();
  }

  fRun = fStartRunLHCPeriod = fEndRunLHCPeriod = -1;
  fRaw = kFALSE;

  fDefaultStorage = 0x0;
}

//_____________________________________________________________________________
void AliCDBManager::SetSpecificStorage(const char* calibType, const char* dbString, Int_t version, Int_t subVersion) {
// sets storage specific for detector or calibration type (works with AliCDBManager::Get(...))

  AliCDBParam *aPar = CreateParameter(dbString);
  if(!aPar) return;
  SetSpecificStorage(calibType, aPar, version, subVersion);
  delete aPar;
}

//_____________________________________________________________________________
void AliCDBManager::SetSpecificStorage(const char* calibType, const AliCDBParam* param, Int_t version, Int_t subVersion) {
// sets storage specific for detector or calibration type (works with AliCDBManager::Get(...))
// Default storage should be defined prior to any specific storages, e.g.:
// AliCDBManager::instance()->SetDefaultStorage("alien://");
// AliCDBManager::instance()->SetSpecificStorage("TPC/*","local://DB_TPC");
// AliCDBManager::instance()->SetSpecificStorage("*/Align/*","local://DB_TPCAlign");
// calibType must be a valid CDB path! (3 level folder structure)
// Specific version/subversion is set in the uniqueid of the AliCDBParam value stored in the
// specific storages map

  if(!fDefaultStorage && !fRaw) {
    AliError("Please activate a default storage first!");
    return;
  }

  AliCDBPath aPath(calibType);
  if(!aPath.IsValid()){
    AliError(Form("Not a valid path: %s", calibType));
    return;
  }

  TObjString *objCalibType = new TObjString(aPath.GetPath());
  if(fSpecificStorages.Contains(objCalibType)){
    AliWarning(Form("Storage \"%s\" already activated! It will be replaced by the new one",
          calibType));
    AliCDBParam *checkPar = dynamic_cast<AliCDBParam*> (fSpecificStorages.GetValue(calibType));
    if(checkPar) delete checkPar;
    delete fSpecificStorages.Remove(objCalibType);
  }
  AliCDBStorage *aStorage = GetStorage(param);
  if(!aStorage) return;

  // Set the unique id of the AliCDBParam stored in the map to store specific version/subversion
  UInt_t uId = ((subVersion+1)<<16) + (version+1);
  AliCDBParam *specificParam = param->CloneParam();
  specificParam->SetUniqueID(uId);
  fSpecificStorages.Add(objCalibType, specificParam);

  if(fStorageMap->Contains(objCalibType)){
    delete fStorageMap->Remove(objCalibType);
  }
  fStorageMap->Add(objCalibType->Clone(), new TObjString(param->GetURI()));

}

//_____________________________________________________________________________
AliCDBStorage* AliCDBManager::GetSpecificStorage(const char* calibType) {
// get storage specific for detector or calibration type

  AliCDBPath calibPath(calibType);
  if(!calibPath.IsValid()) return NULL;

  AliCDBParam *checkPar = (AliCDBParam*) fSpecificStorages.GetValue(calibPath.GetPath());
  if(!checkPar){
    AliError(Form("%s storage not found!", calibType));
    return NULL;
  } else {
    return GetStorage(checkPar);
  }

}

//_____________________________________________________________________________
AliCDBParam* AliCDBManager::SelectSpecificStorage(const TString& path) {
// select storage valid for path from the list of specific storages

  AliCDBPath aPath(path);
  if(!aPath.IsValid()) return NULL;

  TIter iter(&fSpecificStorages);
  TObjString *aCalibType=0;
  AliCDBPath tmpPath("null/null/null");
  AliCDBParam* aPar=0;
  while((aCalibType = (TObjString*) iter.Next())){
    AliCDBPath calibTypePath(aCalibType->GetName());
    if(calibTypePath.Comprises(aPath)) {
      if(calibTypePath.Comprises(tmpPath)) continue;
      aPar = (AliCDBParam*) fSpecificStorages.GetValue(aCalibType);
      tmpPath.SetPath(calibTypePath.GetPath());
    }
  }
  return aPar;
}

//_____________________________________________________________________________
AliCDBEntry* AliCDBManager::Get(const AliCDBPath& path, Int_t runNumber,
    Int_t version, Int_t subVersion) {
  // get an AliCDBEntry object from the database

  if(runNumber < 0){
    // RunNumber is not specified. Try with fRun
    if (fRun < 0){
      AliError("Run number neither specified in query nor set in AliCDBManager! Use AliCDBManager::SetRun.");
      return NULL;
    }
    runNumber = fRun;
  }

  return Get(AliCDBId(path, runNumber, runNumber, version, subVersion));
}

//_____________________________________________________________________________
AliCDBEntry* AliCDBManager::Get(const AliCDBPath& path,
    const AliCDBRunRange& runRange, Int_t version,
    Int_t subVersion) {
  // get an AliCDBEntry object from the database!

  return Get(AliCDBId(path, runRange, version, subVersion));
}

//_____________________________________________________________________________
AliCDBEntry* AliCDBManager::Get(const AliCDBId& queryId, Bool_t forceCaching) {
// get an AliCDBEntry object from the database

  // check if queryId's path and runRange are valid
  // queryId is invalid also if version is not specified and subversion is!
  if (!queryId.IsValid()) {
    AliError(Form("Invalid query: %s", queryId.ToString().Data()));
    return NULL;
  }

  // queryId is not specified if path contains wildcard or run range= [-1,-1]
  if (!queryId.IsSpecified()) {
    AliError(Form("Unspecified query: %s",
          queryId.ToString().Data()));
    return NULL;
  }

  if(fLock && !(fRun >= queryId.GetFirstRun() && fRun <= queryId.GetLastRun()))
    AliFatal("Lock is ON: cannot use different run number than the internal one!");

  if(fCache && !(fRun >= queryId.GetFirstRun() && fRun <= queryId.GetLastRun()))
    AliWarning("Run number explicitly set in query: CDB cache temporarily disabled!");

  AliCDBEntry *entry=0;

  //first of all, look in the prompt cache, for online overrides
  if(fPromptCache && queryId.GetFirstRun() == fRun)
    entry = (AliCDBEntry*) fPromptEntryCache.GetValue(queryId.GetPath());
  if(entry) {
    AliDebug(2, Form("Object %s retrieved from PROMPT cache !!",queryId.GetPath().Data()));
    return entry;
  }

  // first look into map of cached objects
  if(fCache && queryId.GetFirstRun() == fRun)
    entry = (AliCDBEntry*) fEntryCache.GetValue(queryId.GetPath());
  if(entry) {
    AliDebug(2, Form("Object %s retrieved from cache !!",queryId.GetPath().Data()));
    return entry;
  }

  // if snapshot flag is set, try getting from the snapshot
  // but in the case a specific storage is specified for this path
  AliCDBParam *aPar=SelectSpecificStorage(queryId.GetPath());
  if(!aPar){
    if(fSnapshotMode && queryId.GetFirstRun() == fRun)
    {
      entry = GetEntryFromSnapshot(queryId.GetPath());
      if(entry) {
        AliInfo(Form("Object \"%s\" retrieved from the snapshot.",queryId.GetPath().Data()));
        if(queryId.GetFirstRun() == fRun) // no need to check fCache, fSnapshotMode not possible otherwise
          CacheEntry(queryId.GetPath(), entry);

        if(!fIds->Contains(&entry->GetId()))
          fIds->Add(entry->GetId().Clone());

        return entry;
      }
    }
  }

  // Entry is not in cache (and, in case we are in snapshot mode, not in the snapshot either)
  // => retrieve it from the storage and cache it!!
  if(!fDefaultStorage) {
    AliError("No storage set!");
    return NULL;
  }

  Int_t version = -1, subVersion = -1;
  AliCDBStorage *aStorage=0;
  if(aPar) {
    aStorage=GetStorage(aPar);
    TString str = aPar->GetURI();
    UInt_t uId = aPar->GetUniqueID();
    version = Int_t(uId&0xffff) - 1;
    subVersion = Int_t(uId>>16) - 1;
    AliDebug(2,Form("Looking into storage: %s",str.Data()));
  } else {
    aStorage=GetDefaultStorage();
    AliDebug(2,"Looking into default storage");
  }

  AliCDBId finalQueryId(queryId);
  if(version >= 0) {
    AliDebug(2,Form("Specific version set to: %d", version));
    finalQueryId.SetVersion(version);
  }
  if(subVersion >= 0) {
    AliDebug(2,Form("Specific subversion set to: %d", subVersion));
    finalQueryId.SetSubVersion(subVersion);
  }
  entry = aStorage->Get(finalQueryId);

  if(entry && fCache && (queryId.GetFirstRun()==fRun || forceCaching)){
    CacheEntry(queryId.GetPath(), entry);
  }

  if(entry && !fIds->Contains(&entry->GetId())){
    fIds->Add(entry->GetId().Clone());
  }

  return entry;
}

//_____________________________________________________________________________
AliCDBEntry* AliCDBManager::GetEntryFromSnapshot(const char* path) {
// get the entry from the open snapshot file

  TString sPath(path);
  sPath.ReplaceAll("/","*");
  if(!fSnapshotFile){
    AliError("No snapshot file is open!");
    return 0;
  }
  AliCDBEntry *entry = dynamic_cast<AliCDBEntry*>(fSnapshotFile->Get(sPath.Data()));
  if(!entry){
    AliDebug(2,Form("Cannot get a CDB entry for \"%s\" from snapshot file",path));
    return 0;
  }

  return entry;
}

//_____________________________________________________________________________
Bool_t AliCDBManager::SetSnapshotMode(const char* snapshotFileName) {
// set the manager in snapshot mode

  if(!fCache){
    AliError("Cannot set the CDB manage in snapshot mode if the cache is not active!");
    return kFALSE;
  }

  //open snapshot file
  TString snapshotFile(snapshotFileName);
  if(snapshotFile.BeginsWith("alien://")){
    if(!gGrid) {
      TGrid::Connect("alien://","");
      if(!gGrid) {
        AliError("Connection to alien failed!");
        return kFALSE;
      }
    }
  }

  fSnapshotFile = TFile::Open(snapshotFileName);
  if (!fSnapshotFile || fSnapshotFile->IsZombie()){
    AliError(Form("Cannot open file %s",snapshotFileName));
    return kFALSE;
  }

  AliInfo("The CDB manager is set in snapshot mode!");
  fSnapshotMode = kTRUE;
  return kTRUE;

}

//_____________________________________________________________________________
const char* AliCDBManager::GetURI(const char* path) {
// return the URI of the storage where to look for path

  if(!IsDefaultStorageSet()) return 0;

  AliCDBParam *aPar=SelectSpecificStorage(path);

  if(aPar) {
    return aPar->GetURI().Data();

  } else {
    return GetDefaultStorage()->GetURI().Data();
  }

  return 0;
}

//_____________________________________________________________________________
Int_t AliCDBManager::GetStartRunLHCPeriod(){
// get the first run of validity
// for the current period
// if set
  if(fStartRunLHCPeriod==-1)
    AliWarning("Run-range not yet set for the current LHC period.");
  return fStartRunLHCPeriod;
}

//_____________________________________________________________________________
Int_t AliCDBManager::GetEndRunLHCPeriod(){
// get the last run of validity
// for the current period
// if set
  if(fEndRunLHCPeriod==-1)
    AliWarning("Run-range not yet set for the current LHC period.");
  return fEndRunLHCPeriod;
}

//_____________________________________________________________________________
TString AliCDBManager::GetLHCPeriod(){
// get the current LHC period string
//
  if(fLHCPeriod.IsWhitespace() || fLHCPeriod.IsNull())
    AliWarning("LHC period (OCDB folder) not yet set");
  return fLHCPeriod;
}

//_____________________________________________________________________________
AliCDBId* AliCDBManager::GetId(const AliCDBPath& path, Int_t runNumber,
    Int_t version, Int_t subVersion) {
  // get the AliCDBId of the valid object from the database (does not retrieve the object)
  // User must delete returned object!

  if(runNumber < 0){
    // RunNumber is not specified. Try with fRun
    if (fRun < 0){
      AliError("Run number neither specified in query nor set in AliCDBManager! Use AliCDBManager::SetRun.");
      return NULL;
    }
    runNumber = fRun;
  }

  return GetId(AliCDBId(path, runNumber, runNumber, version, subVersion));
}

//_____________________________________________________________________________
AliCDBId* AliCDBManager::GetId(const AliCDBPath& path,
    const AliCDBRunRange& runRange, Int_t version,
    Int_t subVersion) {
  // get the AliCDBId of the valid object from the database (does not retrieve the object)
  // User must delete returned object!

  return GetId(AliCDBId(path, runRange, version, subVersion));
}

//_____________________________________________________________________________
AliCDBId* AliCDBManager::GetId(const AliCDBId& query) {
// get the AliCDBId of the valid object from the database (does not retrieve the object)
// User must delete returned object!

  if(!fDefaultStorage) {
    AliError("No storage set!");
    return NULL;
  }

  // check if query's path and runRange are valid
  // query is invalid also if version is not specified and subversion is!
  if (!query.IsValid()) {
    AliError(Form("Invalid query: %s", query.ToString().Data()));
    return NULL;
  }

  // query is not specified if path contains wildcard or run range= [-1,-1]
  if (!query.IsSpecified()) {
    AliError(Form("Unspecified query: %s",
          query.ToString().Data()));
    return NULL;
  }

  if(fCache && query.GetFirstRun() != fRun)
    AliWarning("Run number explicitly set in query: CDB cache temporarily disabled!");

  AliCDBEntry* entry = 0;

  // first look into map of cached objects
  if(fCache && query.GetFirstRun() == fRun)
    entry = (AliCDBEntry*) fEntryCache.GetValue(query.GetPath());

  if(entry) {
    AliDebug(2, Form("Object %s retrieved from cache !!",query.GetPath().Data()));
    return dynamic_cast<AliCDBId*> (entry->GetId().Clone());
  }

  // Entry is not in cache -> retrieve it from CDB and cache it!!
  AliCDBStorage *aStorage=0;
  AliCDBParam *aPar=SelectSpecificStorage(query.GetPath());

  if(aPar) {
    aStorage=GetStorage(aPar);
    TString str = aPar->GetURI();
    AliDebug(2,Form("Looking into storage: %s",str.Data()));

  } else {
    aStorage=GetDefaultStorage();
    AliDebug(2,"Looking into default storage");
  }

  return aStorage->GetId(query);

}

//_____________________________________________________________________________
TList* AliCDBManager::GetAll(const AliCDBPath& path, Int_t runNumber,
    Int_t version, Int_t subVersion) {
  // get multiple AliCDBEntry objects from the database

  if(runNumber < 0){
    // RunNumber is not specified. Try with fRun
    if (fRun < 0){
      AliError("Run number neither specified in query nor set in AliCDBManager! Use AliCDBManager::SetRun.");
      return NULL;
    }
    runNumber = fRun;
  }

  return GetAll(AliCDBId(path, runNumber, runNumber, version, 	
        subVersion));
}

//_____________________________________________________________________________
TList* AliCDBManager::GetAll(const AliCDBPath& path,
    const AliCDBRunRange& runRange, Int_t version, Int_t subVersion) {
  // get multiple AliCDBEntry objects from the database

  return GetAll(AliCDBId(path, runRange, version, subVersion));
}

//_____________________________________________________________________________
TList* AliCDBManager::GetAll(const AliCDBId& query) {
// get multiple AliCDBEntry objects from the database
// Warning: this method works correctly only for queries of the type "Detector/*"
// 		and not for more specific queries e.g. "Detector/Calib/*" !
// Warning #2: Entries are cached, but GetAll will keep on retrieving objects from OCDB!
// 		To get an object from cache use Get() function

  if(!fDefaultStorage) {
    AliError("No storage set!");
    return NULL;
  }

  if (!query.IsValid()) {
    AliError(Form("Invalid query: %s", query.ToString().Data()));
    return NULL;
  }

  if((fSpecificStorages.GetEntries()>0) && query.GetPath().BeginsWith('*')){
    // if specific storages are active a query with "*" is ambiguous
    AliError("Query too generic in this context!");
    return NULL;
  }

  if (query.IsAnyRange()) {
    AliError(Form("Unspecified run or runrange: %s",
          query.ToString().Data()));
    return NULL;
  }

  if(fLock && query.GetFirstRun() != fRun)
    AliFatal("Lock is ON: cannot use different run number than the internal one!");

  AliCDBParam *aPar = SelectSpecificStorage(query.GetPath());

  AliCDBStorage *aStorage;
  if(aPar) {
    aStorage=GetStorage(aPar);
    AliDebug(2,Form("Looking into storage: %s", aPar->GetURI().Data()));

  } else {
    aStorage=GetDefaultStorage();
    AliDebug(2,Form("Looking into default storage: %s", aStorage->GetURI().Data()));
  }

  TList *result = 0;
  if(aStorage) result = aStorage->GetAll(query);
  if(!result) return 0;

  // loop on result to check whether entries should be re-queried with specific storages
  if(fSpecificStorages.GetEntries()>0 && ! (fSpecificStorages.GetEntries() == 1 && aPar)) {
    AliInfo("Now look into all other specific storages...");

    TIter iter(result);
    AliCDBEntry* chkEntry=0;

    while((chkEntry = dynamic_cast<AliCDBEntry*> (iter.Next()))){
      AliCDBId& chkId = chkEntry->GetId();
      AliDebug(2, Form("Checking id %s ", chkId.GetPath().Data()));
      AliCDBParam *chkPar=SelectSpecificStorage(chkId.GetPath());
      if (!chkPar || aPar == chkPar) continue;
      AliCDBStorage *chkStorage = GetStorage(chkPar);
      AliDebug(2, Form("Found specific storage! %s", chkPar->GetURI().Data()));

      chkId.SetRunRange(query.GetFirstRun(), query.GetLastRun());
      UInt_t uId = chkPar->GetUniqueID();
      Int_t version = -1, subVersion = -1;
      version = Int_t(uId&0xffff) - 1;
      subVersion = Int_t(uId>>16) - 1;
      if(version!=-1){
        chkId.SetVersion(version);
      }else{
        chkId.SetVersion(query.GetVersion());
      }
      if(subVersion!=-1){
        chkId.SetSubVersion(subVersion);
      }else{
        chkId.SetSubVersion(query.GetSubVersion());
      }

      AliCDBEntry *newEntry=0;
      if(chkStorage) newEntry = chkStorage->Get(chkId);
      if(!newEntry) continue;

      // object is found in specific storage: replace entry in the result list!
      chkEntry->SetOwner(1);
      delete result->Remove(chkEntry);
      result->AddFirst(newEntry);
    }

    Int_t nEntries = result->GetEntries();
    AliInfo("After look into other specific storages, result list is:");
    for(int i=0; i<nEntries;i++){
      AliCDBEntry *entry = (AliCDBEntry*) result->At(i);
      AliInfo(Form("%s",entry->GetId().ToString().Data()));
    }
  }

  // caching entries
  TIter iter(result);
  AliCDBEntry* entry=0;
  while((entry = dynamic_cast<AliCDBEntry*> (iter.Next()))){

    if(!fIds->Contains(&entry->GetId())){
      fIds->Add(entry->GetId().Clone());
    }
    if(fCache && (query.GetFirstRun() == fRun)){
      CacheEntry(entry->GetId().GetPath(), entry);
    }
  }


  return result;
}

//_____________________________________________________________________________
Bool_t AliCDBManager::Put(TObject* object, const AliCDBId& id, AliCDBMetaData* metaData, const char* mirrors, DataType type){
// store an AliCDBEntry object into the database

  if (object==0x0) {
    AliError("Null Entry! No storage will be done!");
    return kFALSE;
  }

  AliCDBEntry anEntry(object, id, metaData);
  return Put(&anEntry, mirrors, type);

}


//_____________________________________________________________________________
Bool_t AliCDBManager::Put(AliCDBEntry* entry, const char* mirrors, DataType type){
// store an AliCDBEntry object into the database

  if(type == kPrivate && !fDefaultStorage) {
    AliError("No storage set!");
    return kFALSE;
  }

  if (!entry){
    AliError("No entry!");
    return kFALSE;
  }

  if (entry->GetObject()==0x0){
    AliError("No valid object in CDB entry!");
    return kFALSE;
  }

  if (!entry->GetId().IsValid()) {
    AliError(Form("Invalid entry ID: %s",
          entry->GetId().ToString().Data()));
    return kFALSE;
  }	

  if (!entry->GetId().IsSpecified()) {
    AliError(Form("Unspecified entry ID: %s",
          entry->GetId().ToString().Data()));
    return kFALSE;
  }

  AliCDBId id = entry->GetId();
  AliCDBParam *aPar = SelectSpecificStorage(id.GetPath());

  AliCDBStorage *aStorage=0;

  if(aPar) {
    aStorage=GetStorage(aPar);
  } else {
    switch(type){
      case kCondition:
        aStorage = GetStorage(fCondParam);
        break;
      case kReference:
        aStorage = GetStorage(fRefParam);
        break;
      case kPrivate:
        aStorage = GetDefaultStorage();
        break;
    }
  }

  AliDebug(2,Form("Storing object into storage: %s", aStorage->GetURI().Data()));

  TString strMirrors(mirrors);
  Bool_t result = kFALSE;
  if(!strMirrors.IsNull() && !strMirrors.IsWhitespace())
    result = aStorage->Put(entry, mirrors, type);
  else
    result = aStorage->Put(entry, "", type);

  if(fRun >= 0) QueryCDB();

  return result;


}

//_____________________________________________________________________________
void AliCDBManager::SetMirrorSEs(const char* mirrors) {
// set mirror Storage Elements for the default storage, if it is of type "alien"
  if(fDefaultStorage->GetType() != "alien"){
    AliInfo("The default storage is not of type \"alien\". Settings for Storage Elements are not taken into account!");
    return;
  }
  fDefaultStorage->SetMirrorSEs(mirrors);
}

//_____________________________________________________________________________
const char* AliCDBManager::GetMirrorSEs() const {
// get mirror Storage Elements for the default storage, if it is of type "alien"
  if(fDefaultStorage->GetType() != "alien"){
    AliInfo("The default storage is not of type \"alien\". Settings for Storage Elements are not taken into account!");
    return "";
  }
  return fDefaultStorage->GetMirrorSEs();
}

//_____________________________________________________________________________
void AliCDBManager::CacheEntry(const char* path, AliCDBEntry* entry) {
// cache AliCDBEntry. Cache is valid until run number is changed.

  AliCDBEntry *chkEntry = dynamic_cast<AliCDBEntry*> (fEntryCache.GetValue(path));

  if(chkEntry) {
    AliDebug(2, Form("Object %s already in cache !!", path));
    return;
  } else {
    AliDebug(2,Form("Caching entry %s", path));
  }

  fEntryCache.Add(new TObjString(path), entry);
  AliDebug(2,Form("Cache entries: %d", fEntryCache.GetEntries()));

}

//_____________________________________________________________________________
void AliCDBManager::Print(Option_t* /*option*/) const {
// Print list of active storages and their URIs

  TString output=Form("Run number = %d; ",fRun);
  output += "Cache is ";
  if(!fCache) output += "NOT ";
  output += Form("ACTIVE; Number of active storages: %d\n",fActiveStorages.GetEntries());

  if(fDefaultStorage) {
    output += Form("\t*** Default Storage URI: \"%s\"\n",fDefaultStorage->GetURI().Data());
    //		AliInfo(output.Data());
  }
  if(fSpecificStorages.GetEntries()>0) {
    TIter iter(fSpecificStorages.GetTable());
    TPair *aPair=0;
    Int_t i=1;
    while((aPair = (TPair*) iter.Next())){
      output += Form("\t*** Specific storage %d: Path \"%s\" -> URI \"%s\"\n",
          i++, ((TObjString*) aPair->Key())->GetName(),
          ((AliCDBParam*) aPair->Value())->GetURI().Data());
    }
  }
  if(fDrainStorage) {
    output += Form("*** Drain Storage URI: %s\n",fDrainStorage->GetURI().Data());
  }
  AliInfo(output.Data());
}

//_____________________________________________________________________________
void AliCDBManager::SetRun(Int_t run) {
// Sets current run number.
// When the run number changes the caching is cleared.

  if(fRun == run)
    return;

  if(fLock && fRun >= 0) {
    AliFatal("Lock is ON, cannot reset run number!");
  }	

  fRun = run;

  if (fRaw) {
    // here the LHCPeriod xml file is parsed; the string containing the correct period is returned; the default storage is set
    if (fStartRunLHCPeriod <= run && fEndRunLHCPeriod >= run){
      AliInfo("LHCPeriod alien folder for current run already in memory");
    }else{
      SetDefaultStorageFromRun(fRun);
      if(fEntryCache.GetEntries()!=0) ClearCache();
      return;
    }
  }
  ClearCache();
  QueryCDB();
}

//_____________________________________________________________________________
void AliCDBManager::ClearPromptCache(){
// clear AliCDBEntry prompt cache

  AliDebug(2, Form("PROMPT cache entries to be deleted: %d",fPromptEntryCache.GetEntries()));

  fPromptEntryCache.DeleteAll();

  AliDebug(2, Form("After deleting - PROMPT cache entries: %d",fPromptEntryCache.GetEntries()));
}

//_____________________________________________________________________________
void AliCDBManager::ClearCache(){
// clear AliCDBEntry cache

  AliDebug(2, Form("Cache entries to be deleted: %d",fEntryCache.GetEntries()));

  /*
  // To clean entries one by one
  TIter iter(fEntryCache.GetTable());
  TPair* pair=0;
  while((pair= dynamic_cast<TPair*> (iter.Next()))){

  TObjString* key = dynamic_cast<TObjString*> (pair->Key());
  AliCDBEntry* entry = dynamic_cast<AliCDBEntry*> (pair->Value());
  AliDebug(2, Form("Deleting entry: %s", key->GetName()));
  if (entry) delete entry;
  delete fEntryCache.Remove(key);
  }
  */
  fEntryCache.DeleteAll();
  AliDebug(2, Form("After deleting - Cache entries: %d",fEntryCache.GetEntries()));
}

//_____________________________________________________________________________
void AliCDBManager::UnloadFromCache(const char* path){
// unload cached object
// that is remove the entry from the cache and the id from the list of ids
//
  if(!fActiveStorages.GetEntries()) {
    AliDebug(2, Form("No active storages. Object \"%s\" is not unloaded from cache", path));
    return;
  }

  AliCDBPath queryPath(path);
  if(!queryPath.IsValid()) return;

  if(!queryPath.IsWildcard()) { // path is not wildcard, get it directly from the cache and unload it!
    if(fEntryCache.Contains(path)){
      AliDebug(2, Form("Unloading object \"%s\" from cache and from list of ids", path));
      TObjString pathStr(path);
      delete fEntryCache.Remove(&pathStr);
      // we do not remove from the list of Id's (it's not very coherent but we leave the
      // id for the benefit of the userinfo)
      /*
         TIter iter(fIds);
         AliCDBId *id = 0;
         while((id = dynamic_cast<AliCDBId*> (iter.Next()))){
         if(queryPath.Comprises(id->GetPath()))
         delete fIds->Remove(id);
         }*/
    } else {
      AliWarning(Form("Cache does not contain object \"%s\"!", path));
    }
    AliDebug(2, Form("Cache entries: %d",fEntryCache.GetEntries()));
    return;
  }

  // path is wildcard: loop on the cache and unload all comprised objects!
  TIter iter(fEntryCache.GetTable());
  TPair* pair = 0;
  Int_t removed=0;

  while((pair = dynamic_cast<TPair*> (iter.Next()))){
    AliCDBPath entryPath = pair->Key()->GetName();
    if(queryPath.Comprises(entryPath)) {
      AliDebug(2, Form("Unloading object \"%s\" from cache and from list of ids", entryPath.GetPath().Data()));
      TObjString pathStr(entryPath.GetPath());
      delete fEntryCache.Remove(&pathStr);
      removed++;

      // we do not remove from the list of Id's (it's not very coherent but we leave the
      // id for the benefit of the userinfo)
      /*
         TIter iterids(fIds);
         AliCDBId *anId = 0;
         while((anId = dynamic_cast<AliCDBId*> (iterids.Next()))){
         AliCDBPath aPath = anId->GetPath();
         TString aPathStr = aPath.GetPath();
         if(queryPath.Comprises(aPath)) {
         delete fIds->Remove(anId);
         }
         }*/
    }
  }
  AliDebug(2,Form("Cache entries and ids removed: %d   Remaining: %d",removed,fEntryCache.GetEntries()));
}

//_____________________________________________________________________________
void AliCDBManager::DestroyActiveStorages() {
// delete list of active storages

  fActiveStorages.DeleteAll();
  fSpecificStorages.DeleteAll();
}

//_____________________________________________________________________________
void AliCDBManager::DestroyActiveStorage(AliCDBStorage* /*storage*/) {
// destroys active storage

  /*
     TIter iter(fActiveStorages.GetTable());
     TPair* aPair;
     while ((aPair = (TPair*) iter.Next())) {
     if(storage == (AliCDBStorage*) aPair->Value())
     delete fActiveStorages.Remove(aPair->Key());
     storage->Delete(); storage=0x0;
     }
     */

}

//_____________________________________________________________________________
void AliCDBManager::QueryCDB() {
// query default and specific storages for files valid for fRun. Every storage loads the Ids into its list.

  if (fRun < 0){
    AliError("Run number not yet set! Use AliCDBManager::SetRun.");
    return;
  }
  if (!fDefaultStorage){
    AliError("Default storage is not set! Use AliCDBManager::SetDefaultStorage");
    return;
  }
  if(fDefaultStorage->GetType() == "alien" || fDefaultStorage->GetType() == "local"){
    fDefaultStorage->QueryCDB(fRun);
  //} else {
  //	AliDebug(2,"Skipping query for valid files, it used only in grid...");
  }

  TIter iter(&fSpecificStorages);
  TObjString *aCalibType=0;
  AliCDBParam* aPar=0;
  while((aCalibType = dynamic_cast<TObjString*> (iter.Next()))){
    aPar = (AliCDBParam*) fSpecificStorages.GetValue(aCalibType);
    if(aPar) {
      AliDebug(2,Form("Querying specific storage %s",aCalibType->GetName()));
      AliCDBStorage *aStorage = GetStorage(aPar);
      if(aStorage->GetType() == "alien" || aStorage->GetType() == "local"){
        aStorage->QueryCDB(fRun, aCalibType->GetName());
      } else {
        AliDebug(2,
            "Skipping query for valid files, it is used only in grid...");
      }
    }
  }
}

//______________________________________________________________________________________________
const char* AliCDBManager::GetDataTypeName(DataType type) {
// returns the name (string) of the data type

  switch (type){
    case kCondition: return "Conditions";
    case kReference: return "Reference";
    case kPrivate: return "Private";
  }
  return 0;

}

//______________________________________________________________________________________________
Bool_t AliCDBManager::DiffObjects(const char *cdbFile1, const char *cdbFile2) const {
// Compare byte-by-byte the objects contained in the CDB entry in two different files,
// whose name is passed as input
// Return value:
//   kTRUE - in case the content of the OCDB object (persistent part) is exactly the same
//   kFALSE - otherwise

  TString f1Str(cdbFile1);
  TString f2Str(cdbFile2);
  if (!gGrid && ( f1Str.BeginsWith("alien://") || f2Str.BeginsWith("alien://") ))
    TGrid::Connect("alien://");

  TFile * f1 = TFile::Open(cdbFile1);
  if (!f1){
    Printf("Cannot open file \"%s\"",cdbFile1);
    return kFALSE;
  }
  TFile * f2 = TFile::Open(cdbFile2);
  if (!f2){
    Printf("Cannot open file \"%s\"",cdbFile2);
    return kFALSE;
  }

  AliCDBEntry * entry1 = (AliCDBEntry*)f1->Get("AliCDBEntry");
  if (!entry1){
    Printf("Cannot get CDB entry from file \"%s\"",cdbFile1);
    return kFALSE;
  }
  AliCDBEntry * entry2 = (AliCDBEntry*)f2->Get("AliCDBEntry");
  if (!entry2){
    Printf("Cannot get CDB entry from file \"%s\"",cdbFile2);
    return kFALSE;
  }

  // stream the two objects in the buffer of two TMessages
  TObject* object1 = entry1->GetObject();
  TObject* object2 = entry2->GetObject();
  TMessage * file1 = new TMessage(TBuffer::kWrite);
  file1->WriteObject(object1);
  Int_t size1 = file1->Length();
  TMessage * file2 = new TMessage(TBuffer::kWrite);
  file2->WriteObject(object2);
  Int_t size2 = file2->Length();
  if (size1!=size2){
    Printf("Problem 2:  OCDB entry of different size (%d,%d)",size1,size2);
    return kFALSE;
  }

  // if the two buffers have the same size, check that they are the same byte-by-byte
  Int_t countDiff=0;
  char* buf1 = file1->Buffer();
  char* buf2 = file2->Buffer();
  //for (Int_t i=0; i<size1; i++)    if (file1->Buffer()[i]!=file2->Buffer()[i]) countDiff++;
  for(Int_t i=0; i<size1; i++)
    if (buf1[i]!=buf2[i]) countDiff++;

  if (countDiff>0){
    Printf("The CDB objects differ by %d bytes.", countDiff);
    return kFALSE;
  }

  Printf("The CDB objects are the same in the two files.");
  return kTRUE;
}

//______________________________________________________________________________________________
void AliCDBManager::InitShortLived() {
// Init the list of short-lived objects
// currently disabled

  fShortLived=0x0;

  // 	fShortLived = new TList();
  // 	fShortLived->SetOwner(1);
  //
  // 	fShortLived->Add(new TObjString("EMCAL/Calib/Data"));
  //
  // 	fShortLived->Add(new TObjString("HMPID/Calib/Nmean"));
  // 	fShortLived->Add(new TObjString("HMPID/Calib/Qthre"));
  //
  // 	fShortLived->Add(new TObjString("ITS/Calib/CalibSPD"));
  //
  // 	fShortLived->Add(new TObjString("MUON/Calib/Gains"));
  // 	fShortLived->Add(new TObjString("MUON/Calib/HV"));
  // 	fShortLived->Add(new TObjString("MUON/Calib/Pedestals"));
  //
  // 	fShortLived->Add(new TObjString("PHOS/Calib/CpvGainPedestals"));
  // 	fShortLived->Add(new TObjString("PHOS/Calib/EmcGainPedestals"));
  //
  // 	fShortLived->Add(new TObjString("PMD/Calib/Data"));
  //
  // 	fShortLived->Add(new TObjString("TRD/Calib/ChamberGainFactor"));
  // 	fShortLived->Add(new TObjString("TRD/Calib/LocalGainFactor"));
  // 	fShortLived->Add(new TObjString("TRD/Calib/ChamberT0"));
  // 	fShortLived->Add(new TObjString("TRD/Calib/LocalT0"));
  // 	fShortLived->Add(new TObjString("TRD/Calib/ChamberVdrift"));
  // 	fShortLived->Add(new TObjString("TRD/Calib/LocalVdrift"));
  //
  // 	fShortLived->Add(new TObjString("ZDC/Calib/Data"));

}

//______________________________________________________________________________________________
Bool_t AliCDBManager::IsShortLived(const char* path) {
// returns the name (string) of the data type

  if(!fShortLived) return kFALSE;

  AliCDBPath aPath(path);
  if(!aPath.IsValid()){
    AliError(Form("Not a valid path: %s", path));
    return kFALSE;
  }

  return fShortLived->Contains(path);

}

//______________________________________________________________________________________________
ULong64_t AliCDBManager::SetLock(Bool_t lock, ULong64_t key){
// To lock/unlock user must provide the key. A new key is provided after
// each successful lock. User should always backup the returned key and
// use it on next access.
  ULong64_t msk=0x0000000ffffffff;
  if (fLock == lock) return 0;  // nothing to be done
  if (lock) {
    // User wants to lock - check his identity
    if (fKey) {
      // Lock has a user - check his key
      if (fKey != key) {
        AliFatal("Wrong key provided to lock CDB. Please remove CDB lock access from your code !");
        return 0;
      }
    }
    // Provide new key
    fKey = gSystem->Now();
    fLock = kTRUE;
    // give to experts possibility to unlock the CDB
    SetUniqueID(UInt_t(fKey&msk));
    return fKey;
  }
  // User wants to unlock - check the provided key
  if (key != fKey) {
    AliFatal("Lock is ON: wrong key provided");
    return 0;
  }
  fLock = kFALSE;
  return key;
}

///////////////////////////////////////////////////////////
// AliCDBManager Parameter class                         //
// interface to specific AliCDBParameter class           //
// (AliCDBGridParam, AliCDBLocalParam, AliCDBDumpParam)  //
///////////////////////////////////////////////////////////

AliCDBParam::AliCDBParam():
  fType(),
  fURI()
{
  // constructor

}

//_____________________________________________________________________________
AliCDBParam::~AliCDBParam() {
  // destructor

}

void AliCDBManager::ExtractBaseFolder(TString& url)
{
  // TBD RS
  // remove everything but the url - 
  // Exact copy of the AliReconstuction::Rectify.... (to be removed)
  // 
  //
  TString sbs;
  if (!(sbs=url("\\?User=[^?]*")).IsNull())                url.ReplaceAll(sbs,"");
  if (!(sbs=url("\\?DBFolder=[^?]*")).IsNull())            url.ReplaceAll("?DB","");
  if (!(sbs=url("\\?SE=[^?]*")).IsNull())                  url.ReplaceAll(sbs,"");
  if (!(sbs=url("\\?CacheFolder=[^?]*")).IsNull())         url.ReplaceAll(sbs,"");
  if (!(sbs=url("\\?OperateDisconnected=[^?]*")).IsNull()) url.ReplaceAll(sbs,"");
  if (!(sbs=url("\\?CacheSize=[^?]*")).IsNull())           url.ReplaceAll(sbs,"");  
  if (!(sbs=url("\\?CleanupInterval=[^?]*")).IsNull())     url.ReplaceAll(sbs,"");  
  Bool_t slash=kFALSE,space=kFALSE;
  while ( (slash=url.EndsWith("/")) || (space=url.EndsWith(" ")) ) {
    if (slash) url = url.Strip(TString::kTrailing,'/');
    if (space) url = url.Strip(TString::kTrailing,' ');
  }
  //url.ToLower();
  //
}

//_____________________________________________________________________________
void AliCDBManager::PromptCacheEntry(const char* path, AliCDBEntry* entry)
{
  //put (or replace) a CDB entry in the prompt cache, and enable the cache
  //this is owner of object, so old entry is destroyed

  TPair* pair = static_cast<TPair*>(fEntryCache.FindObject(path));
  if (pair)
  {
    AliDebug(2, Form("Object %s already in cache, replacing!!", path));
    TObject* key = pair->Key();
    TObject* value = pair->Value();
    fEntryCache.Remove(key);
    delete key;
    //delete value; //We must not delete this here. Actually, the CDB Manager CANNOT take ownership. If some function has querried that Object, it might still have a pointer to it and we DO NOT know.
  }
  AliDebug(2,Form("Caching entry %s", path));
  fPromptEntryCache.Add(new TObjString(path), entry);
  SetPromptCacheFlag(kTRUE);
}

