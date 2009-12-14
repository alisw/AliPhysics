/* $Id: AliEventStatsTask.h 35782 2009-10-22 11:54:31Z jgrosseo $ */

#ifndef ALIEVENTSTATSTASK_H
#define ALIEVENTSTATSTASK_H

#include "AliAnalysisTaskSE.h"

class AliPhysicsSelection;

class AliEventStatsTask : public AliAnalysisTaskSE {
  public:
    AliEventStatsTask(const char* opt = "");
    virtual ~AliEventStatsTask();

    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t*);
    virtual void   Terminate(Option_t*);

    void SetOption(const char* opt) { fOption = opt; }

 protected:
    TList* fOutput;                  //! list send on output slot 1

    TString fOption;      // option string  
    
    AliPhysicsSelection* fPhysicsSelection; //! event selection class

 private:
    AliEventStatsTask(const AliEventStatsTask&);
    AliEventStatsTask& operator=(const AliEventStatsTask&);

  ClassDef(AliEventStatsTask, 1);
};

#endif
