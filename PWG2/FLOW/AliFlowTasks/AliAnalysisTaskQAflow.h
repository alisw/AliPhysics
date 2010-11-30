#ifndef ALIANALYSISTASKEFFICIENCYPRIMARIES_CXX
#define ALIANALYSISTASKEFFICIENCYPRIMARIES_CXX

class TObjArray;
class TNtuple;
class AliESDEvent;
class AliFlowEventCuts;
class AliFlowTrackCuts;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskQAflow: public AliAnalysisTaskSE
{
  public:
    AliAnalysisTaskQAflow();
    AliAnalysisTaskQAflow(const char* name);
    virtual ~AliAnalysisTaskQAflow();
    
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);


    void SetTrackCuts(AliFlowTrackCuts* trackcuts) {fTrackCuts=trackcuts;}
    void SetEventCuts(AliFlowEventCuts* eventcuts) {fEventCuts=eventcuts;}
    void SetFillNTuple(Bool_t b=kTRUE) {fFillNtuple=b;}
    void SetDoCorrelations(Bool_t b=kTRUE) {fDoCorrelations=b;}
  
  private:
    TObjArray* fOutput; //output histograms
    Bool_t fFillNtuple;  //fil; the ntuple
    Bool_t fDoCorrelations; //do the slow loopinloop correlations
    TNtuple* fNtuple; //output ntuple
    AliFlowEventCuts* fEventCuts; //AliAnalysisCuts - applied before analysis - for comparing different event classes
    AliFlowTrackCuts* fTrackCuts; //AliFlowTrackCuts go in here

    AliAnalysisTaskQAflow(const AliAnalysisTaskQAflow&); // not implemented
    AliAnalysisTaskQAflow& operator=(const AliAnalysisTaskQAflow&); // not implemented

    ClassDef(AliAnalysisTaskQAflow, 2); // example of analysis 
};

#endif

