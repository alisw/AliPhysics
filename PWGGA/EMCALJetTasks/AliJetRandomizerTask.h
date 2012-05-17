#ifndef ALIJETRANDOMIZERTASK_H
#define ALIJETRANDOMIZERTASK_H

// $Id$

class TClonesArray;

#include "AliJetModelBaseTask.h"

class AliJetRandomizerTask : public AliJetModelBaseTask {
 public:
  AliJetRandomizerTask();
  AliJetRandomizerTask(const char *name); 
  virtual ~AliJetRandomizerTask();

  virtual void         UserExec(Option_t* /*option*/);

 protected:

  virtual void         Run();               // do randomizing

 private:
  AliJetRandomizerTask(const AliJetRandomizerTask&);            // not implemented
  AliJetRandomizerTask &operator=(const AliJetRandomizerTask&); // not implemented

  ClassDef(AliJetRandomizerTask, 1) // Jet randomizer task
};
#endif
