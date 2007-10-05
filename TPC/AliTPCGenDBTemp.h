/////////////////////////////////////////////////////////////////
// Class to generate temperature sensor data base entries.
//
// Existing data base structure read at start of processsing.
// 20/12-2006 HH.
// Modification log:
/////////////////////////////////////////////////////////////////

#ifndef AliTPCGenDBTemp_h
#define AliTPCGenDBTemp_h

#include <TROOT.h>
#include <TFile.h>
#include <TObjArray.h>

#include "AliTPCSensorTempArray.h"
#include "AliLog.h"
#include "AliDCSGenDB.h"

class AliTPCGenDBTemp : public AliDCSGenDB {

public:

// constructors

  AliTPCGenDBTemp();
  AliTPCGenDBTemp(const char *defaultStorage, const char *specificStorage);
  ~AliTPCGenDBTemp();
  void            MakeCalib(const char *file, const char *fMap,
                            const TTimeStamp& startTime,
			    const TTimeStamp& endTime, Int_t run);

// functionality

  static TClonesArray *  ReadList(const char* fname);
  static TTree        *  ReadListTree(const char* fname);

// getters/setters


private:
  AliTPCGenDBTemp(const AliTPCGenDBTemp& org);
  AliTPCGenDBTemp& operator= (const AliTPCGenDBTemp& org);

   ClassDef(AliTPCGenDBTemp,1)
};
#endif
