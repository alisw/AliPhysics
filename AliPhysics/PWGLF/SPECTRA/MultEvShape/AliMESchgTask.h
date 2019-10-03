#ifndef ALIMESCHGTASK_H
#define ALIMESCHGTASK_H

////////////////////////////////////////////////////////////////////////////
//  Chg task for Multiplicity and Event Shape group                       //
//  Authors:                                                              //
//    Andrei Herghelegiu <aherghe@niham.nipne.ro>                         //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#ifndef ALIMESBASETASK_H
#include "AliMESbaseTask.h"
#endif

class AliMESchgTask : public AliMESbaseTask
{
public:
  AliMESchgTask();
  AliMESchgTask(const char *name);
  virtual ~AliMESchgTask();

  Bool_t         BuildQAHistos();
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *opt);
  virtual Bool_t PostProcess();

private:
  AliMESchgTask(const AliMESchgTask&);
  AliMESchgTask& operator=(const AliMESchgTask&);

  ClassDef(AliMESchgTask, 1)            // Chg task for the Multi Event Shape
};

#endif

