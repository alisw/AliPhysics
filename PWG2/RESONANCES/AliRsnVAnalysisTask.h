#ifndef ALIRSNVANALYSISTASK_H
#define ALIRSNVANALYSISTASK_H

//
// Class AliRsnVAnalysisTask
//
// Virtual Class derivated from AliAnalysisTaskSE which will be base class
// for all RSN Multi tasks
//
// authors: Martin Vala (martin.vala@cern.ch)
//          Alberto Pulvirenti (alberto.pulvirenti@ct.infn.it)
//

#include "AliAnalysisTaskSE.h"
#include "AliMixInputEventHandler.h"

#include "AliRsnEvent.h"
#include "AliRsnVATProcessInfo.h"

class AliESDEvent;
class AliAODEvent;
class AliMCEvent;

class AliRsnVAnalysisTask : public AliAnalysisTaskSE {
public:

   AliRsnVAnalysisTask(const char *name = "AliRsnVAnalysisTask", Bool_t mcOnly = kFALSE);
   AliRsnVAnalysisTask(const AliRsnVAnalysisTask& copy);
   AliRsnVAnalysisTask& operator= (const AliRsnVAnalysisTask& /*copy*/) { return *this; }
   virtual ~AliRsnVAnalysisTask() {/* Does nothing*/;}

   // basic interface methods
   virtual void    LocalInit();
   virtual Bool_t  UserNotify();
   virtual void    ConnectInputData(Option_t *opt);
   virtual void    UserCreateOutputObjects();
   virtual void    UserExec(Option_t* opt);
   virtual void    UserExecMix(Option_t* option = "");
   virtual void    Terminate(Option_t* opt);

   // customized methods (to be implemented in derived classes)
   virtual void    RsnUserCreateOutputObjects();
   virtual void    RsnUserExec(Option_t*);
   virtual void    RsnUserExecMix(Option_t*);
   virtual void    RsnTerminate(Option_t*);
   virtual Bool_t  RsnEventProcess();

   // getters
   AliRsnEvent*           GetRsnEvent(Int_t i = 0) {return &fRsnEvent[i];}
   AliRsnVATProcessInfo*  GetInfo()                {return &fTaskInfo;}
   Bool_t                 IsMixing()               {return fIsMixing;}
   Bool_t                 IsUsingMixingRange()     {return fUseMixingRange;}

   // setters
   void SetMCOnly(Bool_t mcOnly = kTRUE)                           {fMCOnly = mcOnly;}
   void SetLogType(AliLog::EType_t type, const char *classes = "") {fLogType = type; fLogClassesString = classes;}
   void SetPrintInfoNumber(const Long64_t &num = 100)              {fTaskInfo.SetPrintInfoNumber(num);}
   void SetMixing(Bool_t doMix = kTRUE)                            {fIsMixing = doMix;}
   void UseMixingRange(Bool_t useMixRange = kTRUE)                 {fUseMixingRange = useMixRange;}

protected:

   AliLog::EType_t          fLogType;             //  log type
   TString                  fLogClassesString;    //  all classes string divided with ":"

   AliESDEvent             *fESDEvent[2];         //  ESD event
   AliMCEvent              *fMCEvent[2];          //  MC event
   AliAODEvent             *fAODEventIn[2];       //  AOD event from input
   AliAODEvent             *fAODEventOut[2];      //  AOD event from output from previous taks
 
   Bool_t                   fIsMixing;            //  flag is using mixing
   Bool_t                   fMCOnly;              //  use only MC information
   AliRsnEvent              fRsnEvent[2];         //  interface to event for RSN package

   TList                   *fInfoList;            //! output list for informations
   AliRsnVATProcessInfo     fTaskInfo;            //  task info

   AliMixInputEventHandler *fMixedEH;             //! mixed event hadnler
   Bool_t                   fUseMixingRange;      //  flag

   void                     SetupMixingEvents();
   void                     SetDebugForAllClasses();

   ClassDef(AliRsnVAnalysisTask, 1)
};

#endif
