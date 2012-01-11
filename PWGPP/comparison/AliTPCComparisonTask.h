#ifndef ALITPCCOMPARISONTASK_H
#define ALITPCCOMPARISONTASK_H

//-------------------------------------------------------------------------
//
// This is the PROOF-enabled version of TPC/AliTPCComparison.C macro.
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

class AliTPCComparisonTask: public AliAnalysisTaskSE
{
  public:
    AliTPCComparisonTask();
    AliTPCComparisonTask(const char* name);
    virtual ~AliTPCComparisonTask() {}
    
    virtual void UserCreateOutputObjects();
    virtual void UserExec(Option_t *option);
    virtual void Terminate(Option_t *);
  
  private:
    TList* fListOfHistos; // the list of the output histos
    
    TH1F* fGood;          // good tracks
    TH1F* fFound;         // found tracks
    TH1F* fFake;          // fake tracks
    TH1F* fP;             // phi resolution 
    TH1F* fL;             // lambda resolution
    TH1F* fPt;            // pt resolution
    TH1F* fHmpt;          // high-pt resolution
    TH1F* fE;             // dE/dx for MIP
    TH2F* fEp;            // dE/dx vs momentum
    TH1F* fGoodPhi;       // phi for good tracks
    TH1F* fFoundPhi;      // phi for found tracks
    
    AliTPCComparisonTask(const AliTPCComparisonTask&); // not implemented
    AliTPCComparisonTask& operator=(const AliTPCComparisonTask&); // not implemented

    ClassDef(AliTPCComparisonTask, 1); // example of analysis 
};

#endif
