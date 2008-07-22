#ifndef AliTOFComparisonTask_cxx
#define AliTOFComparisonTask_cxx

class TList;
class TH1F;
class TH2F;
class TClonesArray;

#ifdef __MAKECINT__
#pragma link C++ class AliMCComparisonTrack+;
#endif

#include "AliAnalysisTaskSE.h"

class AliTOFComparisonTask: public AliAnalysisTaskSE
{
  public:
    AliTOFComparisonTask();
    AliTOFComparisonTask(const char* name);
    virtual ~AliTOFComparisonTask() {}
    
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);
  
  private:
    TList* fListOfHistos;
    
    TH1F* hgood;
    TH1F* hfound;
    TH1F* hfake;
    TH1F* hgoodPhi;
    TH1F* hfoundPhi;
    TH1F* hgoodl;
    TH1F* hfakel;
    TH1F* hfoundl;
    
    AliTOFComparisonTask(const AliTOFComparisonTask&); // not implemented
    AliTOFComparisonTask& operator=(const AliTOFComparisonTask&); // not implemented

    ClassDef(AliTOFComparisonTask, 1); // example of analysis 
};

#endif
