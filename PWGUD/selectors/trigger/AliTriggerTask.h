/* $Id: AliTriggerTask.h 35782 2009-10-22 11:54:31Z jgrosseo $ */

#ifndef ALITRIGGERTASK_H
#define ALITRIGGERTASK_H

#include "AliAnalysisTask.h"
#include "AliPWG0Helper.h"
#include "TParameter.h"

class TH1;
class AliESDEvent;
class AliPhysicsSelection;

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
    void SetUseOrbits(Bool_t flag) { fUseOrbits = flag; }
    void SetPhysicsSelection(AliPhysicsSelection* selection) { fPhysicsSelection = selection; }

 protected:
    AliESDEvent *fESD;    //! ESD object
    TList* fOutput;                  //! list send on output slot 0

    TString fOption;      // option string  
    UInt_t fStartTime;    // run start time
    UInt_t fEndTime;      // run end time
    Bool_t fUseOrbits;    // use orbits instead of time stamps on the axes
    
    TParameter<Long_t>* fFirstOrbit; // first orbit occuring
    TParameter<Long_t>* fLastOrbit; // first orbit occuring

    Int_t fNTriggers;                           //! number triggers
    AliTriggerAnalysis::Trigger* fTriggerList;  //! list of triggers
    Int_t fNTriggerClasses;                     //! number of trigger classes
    const char** fTriggerClassesList;           //! list of trigger classes
    TH1*** fStats;                              //! trigger stats
    
    AliPhysicsSelection* fPhysicsSelection; // trigger object

 private:
    AliTriggerTask(const AliTriggerTask&);
    AliTriggerTask& operator=(const AliTriggerTask&);

  ClassDef(AliTriggerTask, 1);
};

#endif
