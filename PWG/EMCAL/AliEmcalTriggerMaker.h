#ifndef ALIEMCALTRIGGERMAKER_H
#define ALIEMCALTRIGGERMAKER_H

// $Id$

class TClonesArray;
class AliEmcalTriggerSetupInfo;

#include "AliAnalysisTaskEmcal.h"

class AliEmcalTriggerMaker : public AliAnalysisTaskEmcal {
 public:
  AliEmcalTriggerMaker();
  AliEmcalTriggerMaker(const char *name);
  virtual ~AliEmcalTriggerMaker();

  void ExecOnce();
  Bool_t Run();

  void SetCaloTriggersOutName(const char *name) { fCaloTriggersOutName      = name; }
  void SetCaloTriggerSetupOutName(const char *name) { fCaloTriggerSetupOutName      = name; }

 protected:  
  TString            fCaloTriggersOutName;    // name of output track array
  TString            fCaloTriggerSetupOutName;    // name of output track array
  TClonesArray      *fCaloTriggersOut;        //!trigger array out
  AliEmcalTriggerSetupInfo  *fCaloTriggerSetupOut;        //!trigger setup

 private:
  AliEmcalTriggerMaker(const AliEmcalTriggerMaker&);            // not implemented
  AliEmcalTriggerMaker &operator=(const AliEmcalTriggerMaker&); // not implemented

  ClassDef(AliEmcalTriggerMaker, 1); // Task to make array of EMCAL particle
};
#endif
