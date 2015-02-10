/// \file CreateIonTailObject.C
/// 
/// CreateIonTailObject
/// Generic script in order to create a IonTail object
/// 
/// Usage: `aliroot -b -q CreateIonTailObject.C`
///
/// Read object:
///
/// ~~~{.cpp} 
/// TFile* f = TFile::Open("/tmp/ocdb/TPC/Calib/IonTail/Run0_999999999_v0_s0.root")
/// AliCDBEntry* entry = static_cast<AliCDBEntry*>(f.Get("AliCDBEntry"))
/// TObjArray * arr = entry->GetObject();
/// arr->ls();
/// ~~~

void CreateIonTailObject( const Char_t* objectPath = "/u/marsland/MAF/OCDBimp/AllResponseFunctions.root",
			  const Int_t minRun = 0, 
			  const Int_t maxRun = AliCDBRunRange::Infinity(),
			  const Char_t* newStoragePath = "local:///tmp/ocdb", 
			  const Char_t* author = "Mesut Arslandok", 
			  const Char_t *comment = "Create new IonTail object",
			  const Char_t* alirootVersion = "05-02-Rev-35") {
  
  // Get Input
  TFile* inFile = TFile::Open(objectPath);
  if (!inFile) {
    printf("File %s could not be found!\n", objectPath);
    return -1;
  }

  TObjArray* arr = static_cast<TObjArray*>(inFile->Get("arrResponse"));


  // -- Write out
  // -------------------------------------------------------------------
  
  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible(author);
  metaData->SetBeamPeriod(1);
  metaData->SetAliRootVersion(alirootVersion); 
  metaData->SetComment(comment);

  AliCDBId id("TPC/Calib/IonTail", minRun, maxRun);
  AliCDBStorage * gStorage = AliCDBManager::Instance()->GetStorage(newStoragePath);
  gStorage->Put(arr, id, metaData);    

  return;
}
