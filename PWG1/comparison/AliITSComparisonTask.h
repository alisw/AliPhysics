#ifndef AliITSComparisonTask_cxx
#define AliITSComparisonTask_cxx

class TList;
class TH1F;
class TH2F;
class TClonesArray;

#ifdef __MAKECINT__
#pragma link C++ class AliMCComparisonTrack+;
#endif

#include "AliAnalysisTaskSE.h"

class AliITSComparisonTask: public AliAnalysisTaskSE
{
  public:
    AliITSComparisonTask();
    AliITSComparisonTask(const char* name);
    virtual ~AliITSComparisonTask() {}
    
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);
  
  private:
    TList* fListOfHistos;
    
    TH1F* hgood;
    TH1F* hfound;
    TH1F* hfake;
    TH1F* hp;
    TH1F* hl;
    TH1F* hpt;
    TH1F* htip;
    TH1F* he;
    TH2F* hep;
    TH1F* hgoodPhi;
    TH1F* hfoundPhi;
    TH1F* hlip;
    
    AliITSComparisonTask(const AliITSComparisonTask&); // not implemented
    AliITSComparisonTask& operator=(const AliITSComparisonTask&); // not implemented

    ClassDef(AliITSComparisonTask, 1); // example of analysis 
};

#endif
