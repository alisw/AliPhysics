/////////////////////////////////////////////////////////////////
// Class to generate temperature sensor data base entries. 
// 
// Existing data base structure read at start of processsing.
// 20/12-2006 HH.
// Modification log:
/////////////////////////////////////////////////////////////////

#ifndef AliTPCDBTemp_h
#define AliTPCDBTemp_h

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
#include "AliTPCSensorTempArray.h"
#include "AliLog.h"

class AliTPCDBTemp : public TObject {

public:

  AliTPCDBTemp();
  AliTPCDBTemp(const AliTPCDBTemp& org);
  ~AliTPCDBTemp();
  AliTPCDBTemp& operator= (const AliTPCDBTemp& org);
  void            Copy(TObject &c) const;
  void            MakeCalib(const char *file, const char *fMap,
                            const TTimeStamp& startTime,
			    const TTimeStamp& endTime, Int_t run);
  void            MakeConfig(const char *file, Int_t firstRun, Int_t lastRun); 
  AliCDBMetaData* CreateMetaObject(const char *objectClassName);
  void            StoreObject(const char* cdbPath, TObject* object, AliCDBMetaData* metaData);
  void            Init(Int_t run);
  void            InitDB(Int_t run);
  void            SetFirstRun(Int_t frun){fFirstRun=frun;}
  void            SetLastRun(Int_t lrun) {fLastRun=lrun;}
  TMap*           SetGraphFile(const char* fname);
  void            SetConfTree(TTree* tree) {fConfTree=tree;}
  TTree*          GetConfTree() const {return fConfTree;} 
  static TClonesArray *  ReadList(const char* fname);
  static TTree        *  ReadListTree(const char* fname);

private:

   Int_t          fFirstRun;	   // first run in validity period
   Int_t          fLastRun;        // last run in validity period
   AliTPCSensorTempArray  *fTemperature; // array of temperature sensors
   AliCDBStorage  *fStorLoc;	   // pointer to CDB storage
   AliTPCcalibDB  *fCalib;	   // calibration object
   AliCDBMetaData *fMetaData;	   // data base metadata
   TTree          *fConfTree;	   // configuration tree
   
   ClassDef(AliTPCDBTemp,1)
};
#endif
