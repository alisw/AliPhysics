#ifndef ALIMESPIDTASK_H
#define ALIMESPIDTASK_H

////////////////////////////////////////////////////////////////////////////
//  PID task for Multiplicity and Event Shape group                      //
//  Authors:                                                              //
//    Cristi Andrei <Cristian.Andrei@cern.ch>                             //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALIMESBASETASK_H
#include "AliMESbaseTask.h"
#endif

class AliMESpidTask : public AliMESbaseTask
{
public:
  AliMESpidTask();
  AliMESpidTask(const char *name);
  virtual ~AliMESpidTask();

  Bool_t         BuildQAHistos();
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *opt);
  virtual Bool_t PostProcess();

private:
  AliMESpidTask(const AliMESpidTask&);
  AliMESpidTask& operator=(const AliMESpidTask&);

  ClassDef(AliMESpidTask, 1)            // PID task for the Multi Event Shape
};

#endif

