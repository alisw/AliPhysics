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

class AliTPCDBTemp {

public:

  AliTPCDBTemp();
  void            MakeCalib(const char *file, const char *fMap,
                            const TTimeStamp& startTime,
			    const TTimeStamp& endTime, Int_t run);
  AliCDBMetaData* CreateMetaObject(const char *objectClassName);
  void            StoreObject(const char* cdbPath, TObject* object, AliCDBMetaData* metaData);
  void            Init(Int_t run);
  void            InitDB(Int_t run);
  void            SetFirstRun(Int_t frun){fFirstRun=frun;}
  void            SetLastRun(Int_t lrun) {fLastRun=lrun;}
  TMap*           SetGraphFile(const char* fname);

private:

   Int_t          fFirstRun;
   Int_t          fLastRun;
   AliTPCSensorTempArray  *fTemperature;
   AliCDBStorage  *fStorLoc;
   AliTPCcalibDB  *fCalib;
   AliCDBMetaData *fMetaData;
};
#endif
