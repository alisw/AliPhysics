
/////////////////////////////////////////////////////////////////
// Class to generate temperature sensor data base entries.
//
// Existing data base structure read at start of processsing.
// 20/12-2006 HH.
// Modification log:
/////////////////////////////////////////////////////////////////

#ifndef AliDCSGenDB_h
#define AliDCSGenDB_h

#include <TROOT.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TH1F.h>
#include <TFile.h>
#include <TObjArray.h>

#include "AliCDBMetaData.h"
#include "AliCDBManager.h"
#include "AliCDBId.h"
#include "AliCDBStorage.h"
#include "AliDCSSensorArray.h"
#include "AliLog.h"
#include "TSystem.h"

class AliDCSGenDB : public TObject {

public:

// Constructors

  AliDCSGenDB();
  AliDCSGenDB(const char* defaultStorage, const char* specificStorage);
  ~AliDCSGenDB();

// Functionality

  void            MakeCalib(const char *file, const char *fMap,
                            const TTimeStamp& startTime,
			    const TTimeStamp& endTime,
			    Int_t firstRun, Int_t lastRun, const char *calibDir);
  void            MakeConfig(const char *file, Int_t firstRun, Int_t lastRun, 
                             const char *confDir);
  AliCDBMetaData* CreateMetaObject(const char *objectClassName);
  void            StoreObject(const char* cdbPath, TObject* object, AliCDBMetaData* metaData);
  void            Init(Int_t run, const char *configDir, 
                       const char *specificDir, 
		       const char *sensorClass="AliDCSSensorArray");
  static TClonesArray *  ReadList(const char* fname, const char *title="dcsConf");
  static TTree        *  ReadListTree(const char* fname, const char *title="dcsConf");

// Getters/Setters

  void            SetFirstRun(Int_t frun){fFirstRun=frun;}
  void            SetLastRun(Int_t lrun) {fLastRun=lrun;}
  TMap*           SetGraphFile(const char* fname);
  void            SetConfTree(TTree *tree) {fConfTree=tree;}
  TTree*          GetConfTree() const {return fConfTree;}
  const TString&  GetSpecificStorage() const { return fSpecificStorage;}
  void            SetSpecificStorage (const TString& specificStorage) { fSpecificStorage=specificStorage; }
  const TString&  GetDefaultStorage() const { return fDefaultStorage;}
  void            SetDefaultStorage (const TString& defaultStorage) { fDefaultStorage=defaultStorage; }
  const AliDCSSensorArray* GetSensorArray() const {return fSensor;}
  void            SetSensorArray(AliDCSSensorArray *arr) { fSensor=arr; }


protected:
  AliDCSGenDB(const AliDCSGenDB& org);
  AliDCSGenDB& operator= (const AliDCSGenDB& org);

   Int_t          fFirstRun;        // first run in validity period
   Int_t          fLastRun;         // last run in validity period
   TString        fSpecificStorage; // specific storage for data base
   TString        fDefaultStorage;  // default storage for data base
   AliDCSSensorArray  *fSensor;     // array of DCS sensors
   AliCDBStorage  *fStorLoc;        // pointer to CDB storage
   AliCDBMetaData *fMetaData;       // data base metadata
   TTree          *fConfTree;	    // configuration tree

   ClassDef(AliDCSGenDB,1)
 };
#endif

