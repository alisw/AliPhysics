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

#include <TString.h>
#include <TClonesArray.h>

#include "AliAnalysisTaskSE.h"

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
   AliRsnMiniAnalysisTask &operator=(const AliRsnMiniAnalysisTask &copy);
   virtual ~AliRsnMiniAnalysisTask();

   void                UseMC(Bool_t yn = kTRUE)           {fUseMC = yn;}
   void                UseESDTriggerMask(UInt_t trgMask)     {fTriggerMask = trgMask;}
   void                UseCentrality(const char *type)    {fUseCentrality = kTRUE; fCentralityType = type; fCentralityType.ToUpper();}
   void                SetUseCentralityPatch(Bool_t isAOD049) {fUseAOD049CentralityPatch = isAOD049;}
   void                UseMultiplicity(const char *type)  {fUseCentrality = kFALSE; fCentralityType = type; fCentralityType.ToUpper();}
   void                UseContinuousMix()                 {fContinuousMix = kTRUE;}
   void                UseBinnedMix()                     {fContinuousMix = kFALSE;}
   void                SetNMix(Int_t nmix)                {fNMix = nmix;}
   void                SetMaxDiffMult (Double_t val)      {fMaxDiffMult  = val;}
   void                SetMaxDiffVz   (Double_t val)      {fMaxDiffVz    = val;}
   void                SetMaxDiffAngle(Double_t val)      {fMaxDiffAngle = val;}
   void                SetEventCuts(AliRsnCutSet *cuts)   {fEventCuts    = cuts;}
   void                SetMixPrintRefresh(Int_t n)        {fMixPrintRefresh = n;}
   void                SetMaxNDaughters(Short_t n)        {fMaxNDaughters = n;}
   void                SetCheckMomentumConservation(Bool_t checkP) {fCheckP = checkP;}
   Int_t               AddTrackCuts(AliRsnCutSet *cuts);
   TClonesArray       *Outputs()                          {return &fHistograms;}
   TClonesArray       *Values()                           {return &fValues;}
   Short_t             GetMaxNDaughters()                 {return fMaxNDaughters;}
   void                SetEventQAHist(TString type,TH2F *histo);
   void                UseBigOutput(Bool_t b=kTRUE) { fBigOutput = b; }

   virtual void        UserCreateOutputObjects();
   virtual void        UserExec(Option_t *);
   virtual void        Terminate(Option_t *);
   virtual void        FinishTaskOutput();

   Int_t               ValueID(AliRsnMiniValue::EType type, Bool_t useMC = kFALSE);
   Int_t               CreateValue(AliRsnMiniValue::EType type, Bool_t useMC = kFALSE);
   AliRsnMiniOutput   *CreateOutput(const char *name, AliRsnMiniOutput::EOutputType type, AliRsnMiniOutput::EComputation src);
   AliRsnMiniOutput   *CreateOutput(const char *name, const char *outType, const char *compType);

private:

   Char_t   CheckCurrentEvent();
   void     FillMiniEvent(Char_t evType);
   Double_t ComputeAngle();
   Double_t ComputeCentrality(Bool_t isESD);
   Double_t ComputeMultiplicity(Bool_t isESD,TString type);
   Double_t ApplyCentralityPatchAOD049();
   void     FillTrueMotherESD(AliRsnMiniEvent *event);
   void     FillTrueMotherAOD(AliRsnMiniEvent *event);
   void     StoreTrueMother(AliRsnMiniPair *pair, AliRsnMiniEvent *event);
   Bool_t   EventsMatch(AliRsnMiniEvent *event1, AliRsnMiniEvent *event2);

   Bool_t               fUseMC;           //  use or not MC info
   Int_t                fEvNum;           //! absolute event counter
   UInt_t               fTriggerMask;   //trigger mask
   Bool_t               fUseCentrality;   //  if true, use centrality for event, otherwise use multiplicity
   TString              fCentralityType;  //  definition used to choose what centrality or multiplicity to use
   Bool_t               fUseAOD049CentralityPatch; //flag to enable AOD049 centrality patch

   Bool_t               fContinuousMix;   //  mixing --> technique chosen (continuous or binned)
   Int_t                fNMix;            //  mixing --> required number of mixes
   Double_t             fMaxDiffMult;     //  mixing --> max difference in multiplicity
   Double_t             fMaxDiffVz;       //  mixing --> max difference in Vz of prim vert
   Double_t             fMaxDiffAngle;    //  mixing --> max difference in reaction plane angle

   TList               *fOutput;          //  output list
   TClonesArray         fHistograms;      //  list of histogram definitions
   TClonesArray         fValues;          //  list of values to be computed
   TH1F                *fHEventStat;      //  histogram of event statistics
   TH1F                *fHAEventsVsMulti; //  histogram of event statistics
   TH2F                *fHAEventVz;       //  histogram of vertex-z vs. multiplicity/centrality
   TH2F                *fHAEventMultiCent;//  histogram of multiplicity vs. centrality
   TH2F                *fHAEventPlane;    //  histogram of event plane vs. multiplicity/centrality

   AliRsnCutSet        *fEventCuts;       //  cuts on events
   TObjArray            fTrackCuts;       //  list of single track cuts
   AliRsnEvent          fRsnEvent;        //! interface object to the event
   TTree               *fEvBuffer;        //! mini-event buffer
   AliTriggerAnalysis  *fTriggerAna;      //! trigger analysis
   AliESDtrackCuts     *fESDtrackCuts;    //! quality cut for ESD tracks
   AliRsnMiniEvent     *fMiniEvent;       //! mini-event cursor
   Bool_t               fBigOutput;       // flag if open file for output list
   Int_t                fMixPrintRefresh; // how often info in mixing part is printed
   Short_t              fMaxNDaughters;   // maximum number of allowed mother's daughter
   Bool_t               fCheckP;          // flag to set in order to check the momentum conservation for mothers

   ClassDef(AliRsnMiniAnalysisTask, 7);   // AliRsnMiniAnalysisTask
};


#endif
