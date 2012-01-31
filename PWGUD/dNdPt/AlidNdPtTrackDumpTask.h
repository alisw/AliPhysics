#ifndef ALIDNDPTTRACKDUMPTASK_H
#define ALIDNDPTTRACKDUMPTASK_H

//------------------------------------------------------------------------------
// Task to dump track information 
// TPC constrained and TP+ITS combined 
// for outliers analysis.
// 
// Author: J.Otwinowski 19/06/2011 
//------------------------------------------------------------------------------

class AliESDEvent;
class AliMCEvent;
class AlidNdPtEventCuts;
class AlidNdPtAcceptanceCuts;
class AliESDtrackCuts;
class AlidNdPt;
class AlidNdPtAnalysis;
class AlidNdPtCorrection;
class AliMagFMaps;
class AliESDEvent; 
class AliMCEvent; 
class TList;
class TTree;
class TTreeSRedirector;

#include "AliTriggerAnalysis.h"
#include "AliAnalysisTaskSE.h"
#include "dNdPt/AlidNdPtHelper.h"

class AlidNdPtTrackDumpTask : public AliAnalysisTaskSE {
 public:


  AlidNdPtTrackDumpTask(const char *name = "AlidNdPtTrackDumpTask");
  virtual ~AlidNdPtTrackDumpTask();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  virtual Bool_t Notify();
  virtual void   FinishTaskOutput();
  void SetUseMCInfo(Bool_t info)           { fUseMCInfo = info; }
  Bool_t IsUseMCInfo() const               { return fUseMCInfo; }
  
  // Process events
  virtual void Process(AliESDEvent *const esdEvent=0, AliMCEvent *const mcEvent=0, AliESDfriend *const esdFriend=0);

  void SetEventCuts(AlidNdPtEventCuts* const cuts)              { fdNdPtEventCuts = cuts; }
  void SetAcceptanceCuts(AlidNdPtAcceptanceCuts* const cuts)    { fdNdPtAcceptanceCuts = cuts; }
  void SetRecAcceptanceCuts(AlidNdPtAcceptanceCuts* const cuts) { fdNdPtRecAcceptanceCuts = cuts; }
  void SetTrackCuts(AliESDtrackCuts* const cuts)                { fEsdTrackCuts = cuts; }
  void SetTrigger(const AliTriggerAnalysis::Trigger trigger)    { fTrigger = trigger; }
  void SetAnalysisMode(const AlidNdPtHelper::AnalysisMode mode) { fAnalysisMode = mode; }

  AlidNdPtEventCuts* GetEventCuts() const                       { return fdNdPtEventCuts; }
  AlidNdPtAcceptanceCuts* GetAcceptanceCuts() const             { return fdNdPtAcceptanceCuts; }
  AlidNdPtAcceptanceCuts* GetRecAcceptanceCuts() const          { return fdNdPtRecAcceptanceCuts; }  
  AliESDtrackCuts* GetTrackCuts() const                         { return fEsdTrackCuts; }
  AliTriggerAnalysis::Trigger GetTrigger() const                { return fTrigger; }
  AlidNdPtHelper::AnalysisMode GetAnalysisMode() const          { return fAnalysisMode; }

  TString GetCentralityEstimator() const {return fCentralityEstimator; }
  void SetCentralityEstimator(TString centEst="V0M") { fCentralityEstimator = centEst; }

  Bool_t IsFromConversion(const Int_t label, AliStack *const stack);
  Bool_t IsFromMaterial(const Int_t label, AliStack *const stack);
  Bool_t IsFromStrangeness(const Int_t label, AliStack *const stack);
  TParticle *GetMother(TParticle *const particle, AliStack *const stack);

  Bool_t ConstrainTPCInner(AliExternalTrackParam *const tpcInnerC, const AliESDVertex* vtx, Double_t b[3]);
  Bool_t ConstrainTrackInner(AliExternalTrackParam *const trackInnerC, const AliESDVertex* vtx, Double_t mass, Double_t b[3]);


 private:

  AliESDEvent *fESD;    //! ESD event
  AliMCEvent *fMC;      //! MC event
  AliESDfriend *fESDfriend; //! ESDfriend event
  TList* fOutput;       //! list send on output slot 0
  TIterator *fPitList;  //! iterator over the output objetcs  

  Bool_t fUseMCInfo;        //! use MC information

  AlidNdPtEventCuts      *fdNdPtEventCuts;      //! event cuts
  AlidNdPtAcceptanceCuts *fdNdPtAcceptanceCuts; //! acceptance cuts  
  AlidNdPtAcceptanceCuts *fdNdPtRecAcceptanceCuts; //! additional recontruction acceptance cuts (not used for MC truth)
  AliESDtrackCuts *fEsdTrackCuts;          //! esd track cuts
  AliTriggerAnalysis::Trigger fTrigger;    //! trigger settings
  AlidNdPtHelper::AnalysisMode fAnalysisMode;   //! analysis mode TPC only, TPC + ITS

  TTree* fOutputSummary;      //! tree to dump output
  TTreeSRedirector* fTreeSRedirector;      //! temp tree to dump output

  TString fCentralityEstimator;     //! use centrality can be "VOM" (default), "FMD", "TRK", "TKL", "CL0", "CL1", "V0MvsFMD", "TKLvsV0M", "ZEMvsZDC"


  AlidNdPtTrackDumpTask(const AlidNdPtTrackDumpTask&); // not implemented
  AlidNdPtTrackDumpTask& operator=(const AlidNdPtTrackDumpTask&); // not implemented
  
  ClassDef(AlidNdPtTrackDumpTask, 1); // example of analysis
};

#endif
