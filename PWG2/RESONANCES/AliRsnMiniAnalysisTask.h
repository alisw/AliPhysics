#ifndef ALIRSNMINIANALYSISTASK_H
#define ALIRSNMINIANALYSISTASK_H

//
// Analysis task for 'mini' sub-package
// Contains all definitions needed for running an analysis:
// -- global event cut
// -- a list of track cuts (any number)
// -- definitions of output histograms
// -- values to be computed.
//

#include "AliAnalysisTaskSE.h"

#include <TString.h>
#include <TClonesArray.h>

#include "AliRsnEvent.h"
#include "AliRsnMiniValue.h"
#include "AliRsnMiniOutput.h"

class TList;
class AliTriggerAnalysis;
class AliRsnMiniEvent;
class AliRsnCutSet;

class AliRsnMiniAnalysisTask : public AliAnalysisTaskSE {

public:

   AliRsnMiniAnalysisTask();
   AliRsnMiniAnalysisTask(const char *name, Bool_t isMC = kFALSE);
   AliRsnMiniAnalysisTask(const AliRsnMiniAnalysisTask &copy);
   AliRsnMiniAnalysisTask& operator=(const AliRsnMiniAnalysisTask &copy);
   virtual ~AliRsnMiniAnalysisTask();

   virtual void        UserCreateOutputObjects();
   virtual void        UserExec(Option_t *option);
   virtual void        Terminate(Option_t *);
   virtual void        FinishTaskOutput();
  
   void                UseMC(Bool_t yn = kTRUE)           {fUseMC = yn;}                     
   void                UseCentrality(const char *type)    {fUseCentrality = kTRUE; fCentralityType = type; fCentralityType.ToUpper();}
   void                UseMultiplicity(const char *type)  {fUseCentrality = kFALSE; fCentralityType = type; fCentralityType.ToUpper();}
   void                SetNMix(Int_t nmix)                {fNMix = nmix;}
   void                SetMaxDiffMult (Double_t val)      {fMaxDiffMult  = val;}
   void                SetMaxDiffVz   (Double_t val)      {fMaxDiffVz    = val;}
   void                SetMaxDiffAngle(Double_t val)      {fMaxDiffAngle = val;}
   void                SetEventCuts(AliRsnCutSet *cuts)   {fEventCuts    = cuts;}
   Int_t               AddTrackCuts(AliRsnCutSet *cuts);
   
   TClonesArray       *Outputs()        {return &fHistograms;}
   TClonesArray       *Values()         {return &fValues;}
   
   Int_t               ValueID(AliRsnMiniValue::EType type, Bool_t useMC = kFALSE);
   Int_t               CreateValue(AliRsnMiniValue::EType type, Bool_t useMC = kFALSE); 
   AliRsnMiniOutput   *CreateOutput(const char *name, AliRsnMiniOutput::EOutputType type, AliRsnMiniOutput::EComputation src);
   AliRsnMiniOutput   *CreateOutput(const char *name, const char *outType, const char *compType);
  
private:

   Char_t   CheckCurrentEvent();
   Double_t ComputeCentrality(Bool_t isESD);
   void     FillTrueMotherESD(AliRsnMiniEvent *event);
   void     FillTrueMotherAOD(AliRsnMiniEvent *event);
   void     StoreTrueMother(AliRsnMiniPair *pair, AliRsnMiniEvent *event);
   void     ProcessEvents(AliRsnMiniEvent *evMain, AliRsnMiniEvent *evMix = 0x0);

   Bool_t               fUseMC;           //  use or not MC info
   Int_t                fEvNum;           //! absolute event counter
   Bool_t               fUseCentrality;   //  if true, use centrality for event, otherwise use multiplicity
   TString              fCentralityType;  //  definition used to choose what centrality or multiplicity to use
                       
   Int_t                fNMix;            //  mixing --> required number of mixes
   Double_t             fMaxDiffMult;     //  mixing --> max difference in multiplicity
   Double_t             fMaxDiffVz;       //  mixing --> max difference in Vz of prim vert
   Double_t             fMaxDiffAngle;    //  mixing --> max difference in reaction plane angle
                       
   TList               *fOutput;          //  output list
   TClonesArray         fHistograms;      //  list of histogram definitions
   TClonesArray         fValues;          //  list of values to be computed
   TH1F                *fHEventStat;      //  histogram of event statistics
                       
   AliRsnCutSet        *fEventCuts;       //  cuts on events
   TObjArray            fTrackCuts;       //  list of single track cuts
   AliRsnEvent          fRsnEvent;        //! interface object to the event
   TTree               *fEvBuffer;        //! mini-event buffer
   TArrayI              fNMixed;          //! array to keep trace of how many times an event was mixed
   AliTriggerAnalysis  *fTriggerAna;      //! trigger analysis
   AliESDtrackCuts     *fESDtrackCuts;    //! quality cut for ESD tracks

   ClassDef(AliRsnMiniAnalysisTask, 1); // AliRsnMiniAnalysisTask
};

inline Int_t AliRsnMiniAnalysisTask::CreateValue(AliRsnMiniValue::EType type, Bool_t useMC)
{
//
// Create a new value in the task,
// and returns its ID, which is needed for setting up histograms.
// If that value was already initialized, returns its ID and does not recreate it.
//

   Int_t valID = ValueID(type, useMC);
   if (valID >= 0 && valID < fValues.GetEntries()) {
      AliInfo(Form("Value '%s' is already created in slot #%d", AliRsnMiniValue::ValueName(type, useMC), valID));
   } else {
      valID = fValues.GetEntries();
      AliInfo(Form("Creating value '%s' in slot #%d", AliRsnMiniValue::ValueName(type, useMC), valID));
      new (fValues[valID]) AliRsnMiniValue(type, useMC);
   }
   
   return valID;
}

inline Int_t AliRsnMiniAnalysisTask::ValueID(AliRsnMiniValue::EType type, Bool_t useMC)
{
//
// Searches if a value computation is initialized
//

   const char *name = AliRsnMiniValue::ValueName(type, useMC);
   TObject *obj = fValues.FindObject(name);
   if (obj) 
      return fValues.IndexOf(obj); 
   else
      return -1;
}

inline AliRsnMiniOutput* AliRsnMiniAnalysisTask::CreateOutput
(const char *name, AliRsnMiniOutput::EOutputType type, AliRsnMiniOutput::EComputation src)
{
//
// Create a new histogram definition in the task,
// which is then returned to the user for its configuration
//

   Int_t n = fHistograms.GetEntries(); 
   AliRsnMiniOutput *newDef = new (fHistograms[n]) AliRsnMiniOutput(name, type, src);
   
   return newDef;
}

inline AliRsnMiniOutput* AliRsnMiniAnalysisTask::CreateOutput
(const char *name, const char *outType, const char *compType)
{
//
// Create a new histogram definition in the task,
// which is then returned to the user for its configuration
//

   Int_t n = fHistograms.GetEntries(); 
   AliRsnMiniOutput *newDef = new (fHistograms[n]) AliRsnMiniOutput(name, outType, compType);
   
   return newDef;
}
#endif
