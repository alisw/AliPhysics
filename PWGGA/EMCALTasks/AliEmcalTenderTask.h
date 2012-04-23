#ifndef ALIEMCALTENDERTASK_H
#define ALIEMCAKTENDERTASK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliEmcalTenderTask.h  $ */


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

  void                      SetEMCALTenderSupply(AliEMCALTenderSupply *supply);

  // Run control
  virtual void              ConnectInputData(Option_t *option);
  virtual void              UserCreateOutputObjects();
  virtual void              UserExec(Option_t *option);
    
  ClassDef(AliEmcalTenderTask,1) 
};
#endif
#ifndef ALIEMCALTENDERTASK_H
#define ALIEMCAKTENDERTASK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliEmcalTenderTask.h  $ */


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

  void                      SetEMCALTenderSupply(AliEMCALTenderSupply *supply);

  // Run control
  virtual void              ConnectInputData(Option_t *option);
  virtual void              UserCreateOutputObjects();
  virtual void              UserExec(Option_t *option);
    
  ClassDef(AliEmcalTenderTask,1) 
};
#endif
