#ifndef ALIRSNANALYSISTASK_H
#define ALIRSNANALYSISTASK_H

#ifndef ALIANALYSISTASKSE_H
#include "AliAnalysisTaskSE.h"
#endif

#include <TObjArray.h>

class TList;
class AliMixInputEventHandler;
class AliMultiInputEventHandler;
class AliRsnLoop;

class AliRsnAnalysisTask : public AliAnalysisTaskSE {

public:

   AliRsnAnalysisTask();
   AliRsnAnalysisTask(const char *name);
   AliRsnAnalysisTask(const AliRsnAnalysisTask&);
   AliRsnAnalysisTask& operator=(const AliRsnAnalysisTask&);
   virtual ~AliRsnAnalysisTask();

   virtual void     UserCreateOutputObjects();
   virtual void     UserExec(Option_t *option);
   virtual void     UserExecMix(Option_t*);
   virtual void     Terminate(Option_t *);
   
   void             AddLoop(AliRsnLoop *object);
   void             InitInputHandlers();
   
private:

   TList                      *fOutput;        //  output list
   TObjArray                   fRsnObjects;    //  list of computation objects

   AliMultiInputEventHandler  *fInputEHMain;   //! input multi handler
   AliMixInputEventHandler    *fInputEHMix;    //! mix input handler

   ClassDef(AliRsnAnalysisTask, 1); // AliRsnAnalysisTask
};

#endif

