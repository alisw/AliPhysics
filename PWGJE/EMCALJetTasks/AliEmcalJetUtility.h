#ifndef ALIEMCALJETUTILITY_H
#define ALIEMCALJETUTILITY_H

#include <TNamed.h>

#include "AliFJWrapper.h"

class AliEmcalJetTask;
class AliEmcalJet;
class AliFJWrapper;

class AliEmcalJetUtility : public TNamed
{
 public:

  AliEmcalJetUtility();
  AliEmcalJetUtility(const char* name);
  AliEmcalJetUtility(const AliEmcalJetUtility &jet);
  AliEmcalJetUtility& operator=(const AliEmcalJetUtility &jet);
  ~AliEmcalJetUtility() {;}

  virtual void Init() = 0;                                                        // Executed only once at the end of AliEmcalJetTask::DoInit()
  virtual void InitEvent(AliFJWrapper& fjw) = 0;                                  //
  virtual void Prepare(AliFJWrapper& fjw) = 0;                                    // Executed for each event at the beginning of AliEmcalJetTask::FillJetBranch()
  virtual void ProcessJet(AliEmcalJet* jet, Int_t ij, AliFJWrapper& fjw) = 0;     // Executed for each jet in the loop in AliEmcalJetTask::FillJetBranch()
  virtual void Terminate(AliFJWrapper& fjw) = 0;                                  // Executed for each event at the end of AliEmcalJetTask::FillJetBranch()

  void SetJetTask(AliEmcalJetTask* jetTask) { fJetTask = jetTask; }

 protected:

  AliEmcalJetTask       *fJetTask     ; // pointer to the main jet task
  Bool_t                 fInit        ; //! whether or not the utility has been initialized

  ClassDef(AliEmcalJetUtility, 1) // Abstract Emcal jet utility class
};
#endif
