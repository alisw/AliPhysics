/////////////////////////////////////////////////////////////////
// Class to generate temperature sensor data base entries.
//
// Existing data base structure read at start of processsing.
// 20/12-2006 HH.
// Modification log:
/////////////////////////////////////////////////////////////////

#ifndef AliTPCGenDBConf_h
#define AliTPCGenDBConf_h

#include <TROOT.h>
#include <TFile.h>
#include <TObjArray.h>
#include <TEnv.h>

#include "AliTPCPreprocessor.h"
#include "AliLog.h"
#include "AliDCSGenDB.h"

class AliTPCGenDBConf : public AliDCSGenDB {

public:

// constructors

  AliTPCGenDBConf();
  AliTPCGenDBConf(const char *defaultStorage, const char *specificStorage);
  ~AliTPCGenDBConf();

// functionality

  void            MakeConfig(const char *file, Int_t firstRun, Int_t lastRun, const char *confDir);

// getters/setters


private:
  AliTPCGenDBConf(const AliTPCGenDBConf& org);
  AliTPCGenDBConf& operator= (const AliTPCGenDBConf& org);

   ClassDef(AliTPCGenDBConf,1)
};
#endif
