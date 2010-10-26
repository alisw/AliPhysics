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
  Bool_t                    fRunChanged;     //! Flag for run change.
  AliESDInputHandler       *fESDhandler;     //! Pointer to ESD input handler
  AliESDEvent              *fESD;            //! Pointer to current ESD event
  AliGRPManager            *fGRPManager;     //! Pointer to GRP manager

  AliTaskCDBconnect(const AliTaskCDBconnect &other);
  AliTaskCDBconnect& operator=(const AliTaskCDBconnect &other);

public:
  AliTaskCDBconnect();
  AliTaskCDBconnect(const char *name);
  virtual ~AliTaskCDBconnect();
  AliESDInputHandler       *GetESDhandler() const {return fESDhandler;}
  AliESDEvent              *GetEvent() const {return fESD;}
  Int_t                     GetRun() const {return fRun;}
  void                      InitGRP();
  Bool_t                    RunChanged() const {return fRunChanged;}
  void                      SetRunNumber(Int_t run) {fRun = run;}
  // Run control
  virtual void              ConnectInputData(Option_t *option = "");
  virtual void              CreateOutputObjects();
  virtual void              LocalInit();
  virtual Bool_t            Notify();
  virtual void              Exec(Option_t *option);
  virtual void              Terminate(Option_t *option);
    
  ClassDef(AliTaskCDBconnect,1)  // Class giving CDB connectivity
};
#endif
