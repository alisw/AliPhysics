#ifndef AliTPCComparisonTask_cxx
#define AliTPCComparisonTask_cxx

class TList;
class TH1F;
class TH2F;
class TClonesArray;

#ifdef __MAKECINT__
#pragma link C++ class AliMCComparisonTrack+;
#endif

#include "AliAnalysisTaskSE.h"

class AliTPCComparisonTask: public AliAnalysisTaskSE
{
  public:
    AliTPCComparisonTask(const char* name = "AliTPCComparisonTask");
    virtual ~AliTPCComparisonTask() {}
    
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
    TH1F* hmpt;
    TH1F* he;
    TH2F* hep;
    TH1F* hgoodPhi;
    TH1F* hfoundPhi;
    
    AliTPCComparisonTask(const AliTPCComparisonTask&); // not implemented
    AliTPCComparisonTask& operator=(const AliTPCComparisonTask&); // not implemented

    ClassDef(AliTPCComparisonTask, 1); // example of analysis 
};

#endif
