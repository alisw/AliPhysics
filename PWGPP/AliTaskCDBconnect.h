#ifndef ALITASKCDBCONNECT_H
#define ALITASKCDBCONNECT_H

//==============================================================================
// TaskCDBconnect - task just allowing connection to CDB (no lock)
//==============================================================================

#ifndef ALIANALYSISTASK_H
#include "AliAnalysisTask.h"
#endif

class AliCDBManager;
class AliGRPManager;
class AliESDEvent;
class AliESDInputHandler;

class AliTaskCDBconnect : public AliAnalysisTask {
private:
  Bool_t                    fFallBackToRaw;  // allow fallback to raw if cvmfs is not mounted
  Int_t                     fRun;            // Current run
  ULong64_t                 fLock;           // CDB lock
  TString                   fStorage;        // Storage (if cvmfs                   
  TObjArray                 fSpecCDBUri;     // Array with specific CDB storages
  AliGRPManager            *fGRPManager;     //! Pointer to GRP manager
  AliTaskCDBconnect(const AliTaskCDBconnect &other);
  AliTaskCDBconnect& operator=(const AliTaskCDBconnect &other);

  void                      InitGRP();
  //
public:
  AliTaskCDBconnect();
  AliTaskCDBconnect(const char *name, const char *storage="raw://", Int_t run=0, Bool_t fallback=kFALSE);
  virtual ~AliTaskCDBconnect();
  Int_t                     GetRun()        const {return fRun;}
  AliGRPManager*            GetGRPManager() const {return (AliGRPManager*)fGRPManager;}
  virtual void              Exec(Option_t *option);
  virtual void              CreateOutputObjects();
  virtual void              ConnectInputData(Option_t *option = "");
  void                      SetSpecificStorage(const char* calibType, const char* dbString,
                                               Int_t version = -1, Int_t subVersion = -1);
  void                      SetFallBackToRaw(Bool_t v)       {fFallBackToRaw = v;}
  Bool_t                    GetFallBackToRaw()         const {return fFallBackToRaw;}
  ClassDef(AliTaskCDBconnect,5)  // Class giving CDB connectivity
};
#endif
