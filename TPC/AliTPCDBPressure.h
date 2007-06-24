/////////////////////////////////////////////////////////////////
// Class to generate temperature sensor data base entries. 
// 
// Existing data base structure read at start of processsing.
// 20/12-2006 HH.
// Modification log:
/////////////////////////////////////////////////////////////////

#ifndef AliTPCDBPressure_h
#define AliTPCDBPressure_h

#include <TROOT.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TH1F.h>
#include <TFile.h>
#include <TObjArray.h>

#include "AliTPCcalibDB.h"
#include "AliCDBMetaData.h"
#include "AliCDBManager.h"
#include "AliCDBId.h"
#include "AliCDBStorage.h"
#include "AliDCSSensorArray.h"
#include "AliLog.h"
#include "TSystem.h"

class AliTPCDBPressure : public TObject {

public:

  AliTPCDBPressure();
  AliTPCDBPressure(const AliTPCDBPressure& org);
  ~AliTPCDBPressure();
  AliTPCDBPressure& operator= (const AliTPCDBPressure& org);
  void            Copy(TObject &c) const;
  void            MakeCalib(const char *file, const char *fMap,
                            const TTimeStamp& startTime,
			    const TTimeStamp& endTime, 
			    Int_t firstRun, Int_t lastRun);
  void            MakeConfig(const char *file, Int_t firstRun, Int_t lastRun); 
  AliCDBMetaData* CreateMetaObject(const char *objectClassName);
  void            StoreObject(const char* cdbPath, TObject* object, AliCDBMetaData* metaData);
  void            Init(Int_t run);
  void            InitDB(Int_t run);
  void            SetFirstRun(Int_t frun){fFirstRun=frun;}
  void            SetLastRun(Int_t lrun) {fLastRun=lrun;}
  TMap*           SetGraphFile(const char* fname);
  void            SetConfTree(TTree *tree) {fConfTree=tree;}
  TTree*          GetConfTree() const {return fConfTree;} 
  static TClonesArray *  ReadList(const char* fname);
  static TTree        *  ReadListTree(const char* fname);

private:

   Int_t          fFirstRun;      // first run in validity period
   Int_t          fLastRun;       // last run in validity period
   AliDCSSensorArray  *fPressure; // array of pressure sensors
   AliCDBStorage  *fStorLoc;      // pointer to CDB storage
   AliTPCcalibDB  *fCalib;        // calibration object    
   AliCDBMetaData *fMetaData;     // data base metadata
   TTree          *fConfTree;	  // configuration tree

   ClassDef(AliTPCDBPressure,1)
 };
#endif

