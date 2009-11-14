/* $Id: AliTriggerTask.h 35782 2009-10-22 11:54:31Z jgrosseo $ */

#ifndef ALITRIGGERTASK_H
#define ALITRIGGERTASK_H

#include "AliAnalysisTask.h"
#include "AliPWG0Helper.h"

class TH1;
class AliESDEvent;

class AliTriggerTask : public AliAnalysisTask {
  public:
    AliTriggerTask(const char* opt = "");
    virtual ~AliTriggerTask();

    virtual void   ConnectInputData(Option_t *);
    virtual void   CreateOutputObjects();
    virtual void   Exec(Option_t*);
    virtual void   Terminate(Option_t*);

    void SetOption(const char* opt) { fOption = opt; }
    void SetTimes(UInt_t start, UInt_t end) { fStartTime = start; fEndTime = end; }

 protected:
    AliESDEvent *fESD;    //! ESD object
    TList* fOutput;                  //! list send on output slot 0

    TString fOption;      // option string
    UInt_t fStartTime;    // run start time
    UInt_t fEndTime;      // run end time

    Int_t fNTriggers;     //! number triggers
    AliPWG0Helper::Trigger* fTriggerList;  //! list of triggers
    TH1** fStats;                 //! trigger stats

 private:
    AliTriggerTask(const AliTriggerTask&);
    AliTriggerTask& operator=(const AliTriggerTask&);

  ClassDef(AliTriggerTask, 1);
};

#endif
