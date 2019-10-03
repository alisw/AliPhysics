#ifndef ALIPHOSTENDERTASK_H
#define ALIPHOSTENDERTASK_H


#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class AliPHOSTenderSupply;

class AliPHOSTenderTask : public AliAnalysisTaskSE {

protected:
  AliPHOSTenderSupply      *fPHOSTender;
  
private:
  AliPHOSTenderTask(const AliPHOSTenderTask &other);
  AliPHOSTenderTask& operator=(const AliPHOSTenderTask &other);

public:  
  AliPHOSTenderTask();
  AliPHOSTenderTask(const char *name);
  virtual ~AliPHOSTenderTask();

  void                      SetPHOSTenderSupply(AliPHOSTenderSupply *supply);
  AliPHOSTenderSupply*      GetPHOSTenderSupply() {return fPHOSTender;}

  virtual void              ConnectInputData(Option_t *option);
  virtual void              UserCreateOutputObjects();
  virtual void              UserExec(Option_t*);
  virtual void              NotifyRun();
    
  ClassDef(AliPHOSTenderTask,1) // Wrapper class to hold tender supply for AOD usage
};
#endif
