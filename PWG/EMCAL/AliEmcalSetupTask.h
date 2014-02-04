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

  void               SetGeoPath(const char *n)  { fGeoPath  = n; }
  void               SetOadbPath(const char *n) { fOadbPath = n; }
  void               SetOcdbPath(const char *n) { fOcdbPath = n; }
  void               SetObjs(const char *n)     { fObjs     = n; }

 protected:
  void               ConnectInputData(Option_t *option = "");
  void               UserExec(Option_t *option);
  void               Setup(Int_t runno);
  void               Terminate(Option_t *option);

  TString            fOcdbPath;        // path to ocdb (def=uselocal)
  TString            fOadbPath;        // path to oadb
  TString            fGeoPath;         // path to geometry
  TString            fObjs;            // string of objects for alignment to apply
  Bool_t             fIsInit;          //!=true then already initialized 
  TString            fLocalOcdb;       //!directory path to local ocdb
  TString            fLocalOcdbStor;   //!storage path to local ocdb

 private:
  AliEmcalSetupTask(const AliEmcalSetupTask&);            // not implemented
  AliEmcalSetupTask &operator=(const AliEmcalSetupTask&); // not implemented

  ClassDef(AliEmcalSetupTask, 5); // Class to setup geometry for EMCal
};

#endif
