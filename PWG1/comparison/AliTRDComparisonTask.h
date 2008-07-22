#ifndef AliTRDComparisonTask_cxx
#define AliTRDComparisonTask_cxx

class TList;
class TH1F;
class TH2F;
class TClonesArray;

#ifdef __MAKECINT__
#pragma link C++ class AliMCComparisonTrack+;
#endif

#include "AliAnalysisTaskSE.h"

class AliTRDComparisonTask: public AliAnalysisTaskSE
{
  public:
    AliTRDComparisonTask();
    AliTRDComparisonTask(const char* name);
    virtual ~AliTRDComparisonTask() {}
    
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
    TH1F* hz;
    TH1F* hc;
    
    AliTRDComparisonTask(const AliTRDComparisonTask&); // not implemented
    AliTRDComparisonTask& operator=(const AliTRDComparisonTask&); // not implemented

    ClassDef(AliTRDComparisonTask, 1); // example of analysis 
};

#endif
