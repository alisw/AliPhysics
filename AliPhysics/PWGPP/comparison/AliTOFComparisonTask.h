#ifndef ALITOFCOMPARISONTASK_H
#define ALITOFCOMPARISONTASK_H

//-------------------------------------------------------------------------
//
// This is the PROOF-enabled version of TOF/AliTOFComparison.C macro.
// Origin:  Andrei.Zalite@cern.ch
//
//-------------------------------------------------------------------------

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
    TList* fListOfHistos;// the list of the output histos
    
    TH1F* fGood;         // good tracks
    TH1F* fFound;        // found tracks
    TH1F* fFake;         // fake tracks
    TH1F* fGoodPhi;      // phi for good tracks
    TH1F* fFoundPhi;     // phi for found tracks
    TH1F* fGoodl;        // tan(lambda) for good tracks
    TH1F* fFakel;        // tan(lambda) for fake tracks
    TH1F* fFoundl;       // tan(lambda) for found tracks
    
    AliTOFComparisonTask(const AliTOFComparisonTask&); // not implemented
    AliTOFComparisonTask& operator=(const AliTOFComparisonTask&); // not implemented

    ClassDef(AliTOFComparisonTask, 1); // example of analysis 
};

#endif
