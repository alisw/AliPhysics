#ifndef ALI_CDB_MANAGER_H
#define ALI_CDB_MANAGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliCDBManager                                            //
//                                                                 //
/////////////////////////////////////////////////////////////////////

#include <TObject.h>
#include <TList.h>
#include <TMap.h>
#include <TSystem.h>
#include <TFile.h>

class AliCDBEntry;
class AliCDBId;
class AliCDBPath;
class AliCDBRunRange;
class AliCDBMetaData;
class AliCDBStorage;
class AliCDBStorageFactory;
class AliCDBParam;

class AliCDBManager: public TObject {

  public:
    enum DataType {kCondition=0, kReference, kPrivate};

    void RegisterFactory(AliCDBStorageFactory* factory);

    Bool_t HasStorage(const char* dbString) const;

    AliCDBParam* CreateParameter(const char* dbString) const;
    AliCDBParam* GetCondParam() const {return fCondParam;}
    AliCDBParam* GetRefParam() const {return fRefParam;}
    static const char* GetDataTypeName(DataType type);

    AliCDBStorage* GetStorage(const char* dbString);

    TList* GetActiveStorages();

    const TMap* GetStorageMap() const {return fStorageMap;}
    const TList* GetRetrievedIds() const {return fIds;}

    const TMap* GetSnapshotMap() const;
    const TList* GetSnapshotRetrievedIds() const;
    void CreateMapListCopy(TMap& mapCopy, TList& listCopy) const;

    void SetDefaultStorage(const char* dbString);
    void SetDefaultStorage(const AliCDBParam* param);
    void SetDefaultStorage(AliCDBStorage *storage);
    void SetDefaultStorage(const char* runType, const char* simType);
    void SetDefaultStorageFromRun(Int_t run);

    Bool_t IsDefaultStorageSet() const {return fDefaultStorage != 0;}
    AliCDBStorage* GetDefaultStorage() const {return fDefaultStorage;}
    void UnsetDefaultStorage();

    void SetSpecificStorage(const char* calibType, const char* dbString, Int_t version = -1, Int_t subVersion = -1);

    //this puts an object in the fPromptEntryCache to override some storage with
    //in-memory object (to be used online)
    void PromptCacheEntry(const char* calibType, AliCDBEntry* entry);
    void SetPromptCacheFlag(Bool_t f) {fPromptCache = f;}
    Bool_t GetPromptCacheFlag() const {return fPromptCache;}

    AliCDBStorage* GetSpecificStorage(const char* calibType);

    void SetDrain(const char* dbString);
    void SetDrain(const AliCDBParam* param);
    void SetDrain(AliCDBStorage *storage);
    void UnsetDrain(){fDrainStorage = 0x0;}
    Bool_t IsDrainSet() const {return fDrainStorage != 0;}
    Bool_t Drain(AliCDBEntry* entry);

    Bool_t SetOCDBUploadMode();
    void UnsetOCDBUploadMode() { fOCDBUploadMode=kFALSE; }
    Bool_t IsOCDBUploadMode() const { return fOCDBUploadMode; }

    AliCDBEntry* Get(const AliCDBId& query, Bool_t forceCaching=kFALSE);
    AliCDBEntry* Get(const AliCDBPath& path, Int_t runNumber=-1,
        Int_t version = -1, Int_t subVersion = -1);
    AliCDBEntry* Get(const AliCDBPath& path, const AliCDBRunRange& runRange,
        Int_t version = -1, Int_t subVersion = -1);
    AliCDBEntry* GetEntryFromSnapshot(const char* path);

    const char* GetURI(const char* path);				 

    TList* GetAll(const AliCDBId& query);
    TList* GetAll(const AliCDBPath& path, Int_t runNumber=-1,
        Int_t version = -1, Int_t subVersion = -1);
    TList* GetAll(const AliCDBPath& path, const AliCDBRunRange& runRange,
        Int_t version = -1, Int_t subVersion = -1); 

    Bool_t Put(TObject* object, const AliCDBId& id, AliCDBMetaData* metaData,
        const char* mirrors="", DataType type=kPrivate);
    Bool_t Put(AliCDBEntry* entry, const char* mirrors="", DataType type=kPrivate);

    void SetCacheFlag(Bool_t cacheFlag) {fCache=cacheFlag;}
    Bool_t GetCacheFlag() const {return fCache;}

    ULong64_t SetLock(Bool_t lockFlag=kTRUE, ULong64_t key=0);
    Bool_t GetLock() const {return fLock;}

    void SetRaw(Bool_t rawFlag){fRaw=rawFlag;}
    Bool_t GetRaw() const {return fRaw;}

    void SetRun(Int_t run);
    Int_t GetRun() const {return fRun;}

    void SetMirrorSEs(const char* mirrors);
    const char* GetMirrorSEs() const;

    void DestroyActiveStorages();
    void DestroyActiveStorage(AliCDBStorage* storage);

    void QueryCDB();

    void Print(Option_t* option="") const;

    static void Destroy();
    ~AliCDBManager();

    void ClearCache();
    void ClearPromptCache();
    void UnloadFromCache(const char* path);
    const TMap* GetEntryCache() const {return &fEntryCache;}

    Bool_t IsShortLived(const char* path);

    static AliCDBManager* Instance(TMap *entryCache = NULL, Int_t run = -1);

    void Init();
    void InitFromCache(TMap *entryCache, Int_t run);
    Bool_t InitFromSnapshot(const char* snapshotFileName, Bool_t overwrite=kTRUE);
    Bool_t SetSnapshotMode(const char* snapshotFileName="OCDB.root");
    void UnsetSnapshotMode() {fSnapshotMode=kFALSE;}
    void DumpToSnapshotFile(const char* snapshotFileName, Bool_t singleKeys) const;
    void DumpToLightSnapshotFile(const char* lightSnapshotFileName) const;

    Int_t GetStartRunLHCPeriod();
    Int_t GetEndRunLHCPeriod();
    TString GetLHCPeriod();
    TString GetCvmfsOcdbTag() const {return fCvmfsOcdb;}

    Bool_t DiffObjects(const char *cdbFile1, const char *cdbFile2) const;
    void ExtractBaseFolder(TString& url); // remove everything but the url from OCDB path


  protected:

    static TString fgkCondUri;	// URI of the Conditions data base folder
    static TString fgkRefUri;	// URI of the Reference data base folder
    static TString fgkMCIdealStorage;	// URI of the MC-Ideal Conditions data base folder form MC data
    static TString fgkMCFullStorage;	// URI of the MC-Full Conditions data base folder form MC data
    static TString fgkMCResidualStorage;	// URI of the MC-Residual Conditions data base folder form MC data
    static TString fgkOCDBFolderXMLfile;	// alien path of the XML file for OCDB folder <--> Run range correspondance

    AliCDBManager() ; 
    AliCDBManager(const AliCDBManager & source);
    AliCDBManager & operator=(const AliCDBManager & source);

    static AliCDBManager* fgInstance; // AliCDBManager instance

    AliCDBStorage* GetStorage(const AliCDBParam* param);
    AliCDBStorage* GetActiveStorage(const AliCDBParam* param);
    void PutActiveStorage(AliCDBParam* param, AliCDBStorage* storage);
    void SetSpecificStorage(const char* calibType, const AliCDBParam* param, Int_t version = -1, Int_t subVersion = -1);
    void AlienToCvmfsUri(TString& uriString) const;
    void ValidateCvmfsCase() const;
    void GetLHCPeriodAgainstAlienFile(Int_t run, TString& lhcPeriod, Int_t& startRun, Int_t& endRun);
    void GetLHCPeriodAgainstCvmfsFile(Int_t run, TString& lhcPeriod, Int_t& startRun, Int_t& endRun);

    void CacheEntry(const char* path, AliCDBEntry* entry);

    AliCDBParam* SelectSpecificStorage(const TString& path);

    AliCDBId* GetId(const AliCDBId& query);
    AliCDBId* GetId(const AliCDBPath& path, Int_t runNumber=-1,
        Int_t version = -1, Int_t subVersion = -1);
    AliCDBId* GetId(const AliCDBPath& path, const AliCDBRunRange& runRange,
        Int_t version = -1, Int_t subVersion = -1);


    //	void Init();
    void InitShortLived();
    //	void InitFromCache(TMap *entryCache, Int_t run);


    TList fFactories; 		//! list of registered storage factories
    TMap fActiveStorages;		//! list of active storages
    TMap fSpecificStorages;         //! list of detector-specific storages
    TMap fEntryCache;    	  	//! cache of the retrieved objects
    TMap fPromptEntryCache;   //! cache for in-memory objects to override objects on storage (to be used online)

    TList* fIds;           	//! List of the retrieved object Id's (to be streamed to file)
    TMap* fStorageMap;      //! list of storages (to be streamed to file)
    TList* fShortLived; 	//! List of short lived objects

    AliCDBStorage *fDefaultStorage;	//! pointer to default storage
    AliCDBStorage *fDrainStorage;	//! pointer to drain storage

    AliCDBParam* fCondParam; 	// Conditions data storage parameters
    AliCDBParam* fRefParam;		// Reference data storage parameters

    Int_t fRun;			//! The run number
    Bool_t fCache;			//! The cache flag
    Bool_t fPromptCache; //! The prompt cache flag
    Bool_t fLock; 	//! Lock flag, if ON default storage and run number cannot be reset

    Bool_t fSnapshotMode;           //! flag saying if we are in snapshot mode
    TFile *fSnapshotFile;
    Bool_t fOCDBUploadMode;         //! flag for uploads to Official CDBs (upload to cvmfs must follow upload to AliEn)

    Bool_t fRaw;   // flag to say whether we are in the raw case
    TString fCvmfsOcdb;       // set from $OCDB_PATH, points to a cvmfs AliRoot package
    Int_t fStartRunLHCPeriod; // 1st run of the LHC period set
    Int_t fEndRunLHCPeriod;   // last run of the LHC period set
    TString fLHCPeriod;       // LHC period alien folder

  private:
    ULong64_t fKey;  //! Key for locking/unlocking


    ClassDef(AliCDBManager, 0);
};


/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliCDBStorageFactory                                     //
//                                                                 //
/////////////////////////////////////////////////////////////////////

class AliCDBParam;
class AliCDBStorageFactory: public TObject {
  friend class AliCDBManager;

  public:
  virtual ~AliCDBStorageFactory(){}
  virtual Bool_t Validate(const char* dbString) = 0;
  virtual AliCDBParam* CreateParameter(const char* dbString) = 0;	

  protected:
  virtual AliCDBStorage* Create(const AliCDBParam* param) = 0;

  ClassDef(AliCDBStorageFactory, 0);
};

/////////////////////////////////////////////////////////////////////
//                                                                 //
//  class AliCDBParam                                              //
//                                                                 //
/////////////////////////////////////////////////////////////////////

class AliCDBParam: public TObject {

  public:

    AliCDBParam();
    virtual ~AliCDBParam();

    const TString& GetType() const {return fType;};
    const TString& GetURI() const {return fURI;};

    virtual AliCDBParam* CloneParam() const = 0;

  protected:

    void SetType(const char* type) {fType = type;};
    void SetURI(const char* uri) {fURI = uri;};

  private:

    TString fType; //! CDB type
    TString fURI;  //! CDB URI

    ClassDef(AliCDBParam, 0);
};

#endif
