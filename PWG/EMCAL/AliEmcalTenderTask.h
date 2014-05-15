#ifndef ALIEMCALTENDERTASK_H
#define ALIEMCALTENDERTASK_H

// $Id$

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

class AliEMCALTenderSupply;

class AliEmcalTenderTask : public AliAnalysisTaskSE {

protected:
  AliEMCALTenderSupply      *fEMCALTender;
  
private:
  AliEmcalTenderTask(const AliEmcalTenderTask &other);
  AliEmcalTenderTask& operator=(const AliEmcalTenderTask &other);

public:  
  AliEmcalTenderTask();
  AliEmcalTenderTask(const char *name);
  virtual ~AliEmcalTenderTask();

  void                      SetEMCALTenderSupply(AliEMCALTenderSupply *supply)       ;
  AliEMCALTenderSupply*     GetEMCALTenderSupply()                              const;

  virtual void              ConnectInputData(Option_t *option);
  virtual void              UserCreateOutputObjects();
  virtual void              UserExec(Option_t*);
    
  ClassDef(AliEmcalTenderTask,2) // Wrapper class to hold tender supply for AOD usage
};
#endif
