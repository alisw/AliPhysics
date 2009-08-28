#ifndef ALIANALYSISTASKLRC_H
#define ALIANALYSISTASKLRC_H

// Analysis task for Long Range Correlation (LRC) analysis using TPC data
// This includes a TList of AliLRCProcess objects that are processing LRC analysis
// for a given Eta window 


// Author : Andrey Ivanov , St.Peterburg State University
// Email: Andrey.Ivanov@cern.ch

// Version line : 3.5
// Version 3.5.5

#include "AliAnalysisTaskSE.h"
#include "AliLRCProcess.h"

class AliAnalysisTaskLRC : public AliAnalysisTaskSE {

public:
 
 
  //Constructors
  
  AliAnalysisTaskLRC(const char *name = "AliAnalysisTaskLRC",Bool_t runKine=kFALSE);
  virtual ~AliAnalysisTaskLRC() {}
  
  //AliAnalysisTaskSE overloading
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
    
  void AddLRCProcess(AliLRCProcess *newProc); //Adds new AliLRCProcess to analysis task
  
  Double_t fMaxPtLimit;  //Max Pt filter
  Double_t fMinPtLimit;  // Min Pt filter 
  Bool_t fDropKineE;     // Force to drop any e+- in Kine

private:
  
  TList fLRCproc;       //  AliLRCProcess objects list
  TList* fOutList;      //! Output data container 
  
  Bool_t fRunKine;      // ESD/AOD  - KINE switch
  
      
  AliAnalysisTaskLRC(const AliAnalysisTaskLRC&); // not implemented
  AliAnalysisTaskLRC& operator=(const AliAnalysisTaskLRC&); // not implemented
  
  ClassDef(AliAnalysisTaskLRC, 1); 
};

#endif
