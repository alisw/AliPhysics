#ifndef ALIEMCALSETUPTASK_H
#define ALIEMCALSETUPTASK_H

// $Id$

class TClonesArray;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"

class AliEmcalSetupTask : public AliAnalysisTaskSE {
 public:
  AliEmcalSetupTask();
  AliEmcalSetupTask(const char *name);
  virtual ~AliEmcalSetupTask();

  void UserExec(Option_t *option);
  void SetGeoPath(const char *n)  { fGeoPath  = n; }
  void SetOadbPath(const char *n) { fOadbPath = n; }
  void SetOcdbPath(const char *n) { fOcdbPath = n; }

 protected:
  TString            fOcdbPath;        // path to ocdb (def=none)
  TString            fOadbPath;        // path to oadb
  TString            fGeoPath;         // path to geometry
  Bool_t             fIsInit;          //!=true then already initialized 

 private:
  AliEmcalSetupTask(const AliEmcalSetupTask&);            // not implemented
  AliEmcalSetupTask &operator=(const AliEmcalSetupTask&); // not implemented

  ClassDef(AliEmcalSetupTask, 2); // Class to setup geometry for EMCal
};

#endif
