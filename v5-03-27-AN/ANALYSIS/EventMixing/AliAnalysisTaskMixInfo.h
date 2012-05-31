#ifndef ALIANALYSISTASKMIXINFO_H
#define ALIANALYSISTASKMIXINFO_H

//
// Class AliAnalysisTaskMixInfo
//
// AliAnalysisTaskMixInfo is task
// for mixing info
//
// authors:
//          Martin Vala (martin.vala@cern.ch)
//

#include "AliLog.h"

#include "AliAnalysisTaskSE.h"

class AliMixInputEventHandler;
class TList;
class AliMixInfo;
class AliAnalysisTaskMixInfo : public AliAnalysisTaskSE {
public:
   AliAnalysisTaskMixInfo(const char *name = "<default name>");
   virtual ~AliAnalysisTaskMixInfo();

   virtual void    UserCreateOutputObjects();
   virtual void    UserExec(Option_t *option);
   virtual void    Terminate(Option_t *);
   virtual void    UserExecMix(Option_t *option = "");
   virtual void    FinishTaskOutput();

   void            InitInputHandlers();
   void            InitMixInfo();
   // sets log type to list of classes
   void            SetLogType(AliLog::EType_t type, TString allClasses = "");
   // sets correctly debug level to AliLog for all classes listed in fLogClassesString
   void            SetDebugForAllClasses();

   void            PrintEventInfo();

private:

   AliMultiInputEventHandler  *fInputEHMain;       //! input multi handler
   AliMixInputEventHandler    *fInputEHMix;        //! mix input handler

   TList                      *fOutputList;        //! output list
   AliMixInfo                 *fMixInfo;           //! mix info

   Long64_t                    fCurrentEntryTmp;   //! temporary current entry number

   AliLog::EType_t             fLogType;           // log type
   TString                     fLogClassesString;  // all classes string divided with ":"

   AliAnalysisTaskMixInfo(const AliAnalysisTaskMixInfo &); // not implemented
   AliAnalysisTaskMixInfo &operator=(const AliAnalysisTaskMixInfo &); // not implemented

   ClassDef(AliAnalysisTaskMixInfo, 1); // example of analysis
};

#endif
