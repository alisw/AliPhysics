#ifndef ALIRSNMINIMonitorTask_H
#define ALIRSNMINIMonitorTask_H

//
// Analysis task for 'mini' sub-package
// Contains all definitions needed for running an analysis:
// -- global event cut
// -- a list of track cuts (any number)
// -- definitions of output histograms
// -- values to be computed.
//

#include <TString.h>
#include <TClonesArray.h>

#include "AliAnalysisTaskSE.h"

#include "AliRsnEvent.h"
#include "AliRsnMiniMonitor.h"

class TList;

class AliTriggerAnalysis;
class AliRsnMiniEvent;
class AliRsnCutSet;

class AliRsnMiniMonitorTask : public AliAnalysisTaskSE {

public:

   AliRsnMiniMonitorTask();
   AliRsnMiniMonitorTask(const char *name, Bool_t isMC = kFALSE);
   AliRsnMiniMonitorTask(const AliRsnMiniMonitorTask &copy);
   AliRsnMiniMonitorTask& operator=(const AliRsnMiniMonitorTask &copy);
   virtual ~AliRsnMiniMonitorTask();

   void                UseMC(Bool_t yn = kTRUE)           {fUseMC = yn;}                     
   void                UseCentrality(const char *type)    {fUseCentrality = kTRUE; fCentralityType = type; fCentralityType.ToUpper();}
   void                UseMultiplicity(const char *type)  {fUseCentrality = kFALSE; fCentralityType = type; fCentralityType.ToUpper();}
   void                SetEventCuts(AliRsnCutSet *cuts)   {fEventCuts    = cuts;}
   Int_t               AddTrackCuts(AliRsnCutSet *cuts);
   TClonesArray       *Outputs()                          {return &fHistograms;}
   
   virtual void        UserCreateOutputObjects();
   virtual void        UserExec(Option_t*);
   virtual void        Terminate(Option_t*);
   
   AliRsnMiniMonitor  *CreateMonitor(const char *name, AliRsnMiniMonitor::EType type, Int_t cutID);
  
private:

   Char_t   CheckCurrentEvent();
   Double_t ComputeCentrality(Bool_t isESD);

   Bool_t               fUseMC;           //  use or not MC info
   Int_t                fEvNum;           //! absolute event counter
   Bool_t               fUseCentrality;   //  if true, use centrality for event, otherwise use multiplicity
   TString              fCentralityType;  //  definition used to choose what centrality or multiplicity to use
                       
   TList               *fOutput;          //  output list
   TClonesArray         fHistograms;      //  list of histogram definitions
                       
   AliRsnCutSet        *fEventCuts;       //  cuts on events
   TObjArray            fTrackCuts;       //  list of single track cuts
   AliRsnEvent          fRsnEvent;        //! interface object to the event
   Bool_t               fBigOutput;       //  flag if open file for output list

   ClassDef(AliRsnMiniMonitorTask, 1);   // AliRsnMiniMonitorTask
};

inline AliRsnMiniMonitor* AliRsnMiniMonitorTask::CreateMonitor
(const char *name, AliRsnMiniMonitor::EType type, Int_t cutID)
{
//
// Create a new histogram definition in the task,
// which is then returned to the user for its configuration
//

   Int_t n = fHistograms.GetEntries(); 
   AliRsnMiniMonitor *newDef = new (fHistograms[n]) AliRsnMiniMonitor(name, type, cutID);
   
   return newDef;
}

#endif
