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

  void         SetRandomizeEta(Int_t opt = 1)    { fRandomizeEta = opt; }

  void         UserExec(Option_t* /*option*/);


 protected:
  void         Run();

  Int_t        fRandomizeEta;  //0 = do not randomize eta; 1 = randomize eta uniformly; 2 = invert eta sign

 private:
  AliJetRandomizerTask(const AliJetRandomizerTask&);            // not implemented
  AliJetRandomizerTask &operator=(const AliJetRandomizerTask&); // not implemented

  ClassDef(AliJetRandomizerTask, 2) // Jet randomizer task
};
#endif
