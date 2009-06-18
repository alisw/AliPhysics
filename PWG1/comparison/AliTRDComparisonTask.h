#ifndef ALITRDCOMPARISONTASK_H
#define ALITRDCOMPARISONTASK_H

//-------------------------------------------------------------------------
//
// This is the PROOF-enabled version of TRD/Macros/AliTRDComparisonV2.C macro.
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
    TList* fListOfHistos;; // The list of output histos
    
    TH1F* fGood;           // good tracks
    TH1F* fFound;          // found tracks
    TH1F* fFake;           // fake tracks
    TH1F* fP;              // phi resolution
    TH1F* fL;              // lambda resolution 
    TH1F* fPt;             // pt resolution
    TH1F* fHmpt;           // Y and Z resolution
    TH1F* fE;              // dE/dx for MIP
    TH2F* fEp;             // dE/dx vs momentum
    TH1F* fGoodPhi;        // phi for good tracks
    TH1F* fFoundPhi;       // phi for found tracks
    TH1F* fZ;              // z resolution
    TH1F* fC;              // number of assigned clusters
    
    AliTRDComparisonTask(const AliTRDComparisonTask&); // not implemented
    AliTRDComparisonTask& operator=(const AliTRDComparisonTask&); // not implemented

    ClassDef(AliTRDComparisonTask, 1); // example of analysis 
};

#endif
