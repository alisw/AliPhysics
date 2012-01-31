#ifndef ALIANALYSISTASKQAPMDFLOW_CXX
#define ALIANALYSISTASKQAPMDFLOW_CXX

class TObjArray;
class TNtuple;
class AliESDEvent;
class AliFlowEventCuts;
class AliFlowTrackCuts;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskQAPmdflow: public AliAnalysisTaskSE
{
  public:
    AliAnalysisTaskQAPmdflow();
    AliAnalysisTaskQAPmdflow(const char* name);
    virtual ~AliAnalysisTaskQAPmdflow();
    
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);


    void SetTrackCuts(AliFlowTrackCuts* trackcutsrp) {fRPTrackCuts=trackcutsrp;}
    void SetPOITrackCuts(AliFlowTrackCuts* trackcutspoi) {fPOITrackCuts=trackcutspoi;}
    void SetEventCuts(AliFlowEventCuts* eventcuts) {fEventCuts=eventcuts;}
  
  private:
    TObjArray* fOutput; //output histograms
    AliFlowEventCuts* fEventCuts; //AliAnalysisCuts - applied before analysis - for comparing different event classes
    AliFlowTrackCuts* fRPTrackCuts; //AliFlowTrackCuts go in here
    AliFlowTrackCuts* fPOITrackCuts; //AliFlowTrackCuts go in here

    AliAnalysisTaskQAPmdflow(const AliAnalysisTaskQAPmdflow&); // not implemented
    AliAnalysisTaskQAPmdflow& operator=(const AliAnalysisTaskQAPmdflow&); // not implemented

    ClassDef(AliAnalysisTaskQAPmdflow, 1); // example of analysis 
};

#endif
