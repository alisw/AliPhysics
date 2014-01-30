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
  Int_t                     fRun;            // Current run
  AliGRPManager            *fGRPManager;     //! Pointer to GRP manager

  AliTaskCDBconnect(const AliTaskCDBconnect &other);
  AliTaskCDBconnect& operator=(const AliTaskCDBconnect &other);

  void                      InitGRP();
  //
public:
  AliTaskCDBconnect();
  AliTaskCDBconnect(const char *name, const char *storage="raw://", Int_t run=0);
  virtual ~AliTaskCDBconnect();
  Int_t                     GetRun()        const {return fRun;}
  AliGRPManager*            GetGRPManager() const {return (AliGRPManager*)fGRPManager;}
  virtual void              Exec(Option_t *option);
  virtual void              CreateOutputObjects();
  void                      SetSpecificStorage(const char* calibType, const char* dbString);
    
  ClassDef(AliTaskCDBconnect,2)  // Class giving CDB connectivity
};
#endif
