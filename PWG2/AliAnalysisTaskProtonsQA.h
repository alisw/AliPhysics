#ifndef AliAnalysisTaskProtonsQA_cxx
#define AliAnalysisTaskProtonsQA_cxx

// Analysis task creating a the 2d y-p_t spectrum of p and antip
// Author: Panos Cristakoglou
class TList;
class AliESDEvent;
class AliMCEvent;
class AliProtonQAAnalysis;

#include "AliAnalysisTask.h"

class AliAnalysisTaskProtonsQA : public AliAnalysisTask {
 public:
  AliAnalysisTaskProtonsQA();
  AliAnalysisTaskProtonsQA(const char *name);
  virtual ~AliAnalysisTaskProtonsQA() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

 private:
  AliESDEvent *fESD;    //ESD object
  AliMCEvent  *fMC;     //MC object

  TList  *fList0; //TList output object
  TList  *fList1; //TList output object
  TList  *fList2; //TList output object
  
  AliProtonQAAnalysis *fAnalysis; //analysis object
 
  AliAnalysisTaskProtonsQA(const AliAnalysisTaskProtonsQA&); // not implemented
  AliAnalysisTaskProtonsQA& operator=(const AliAnalysisTaskProtonsQA&); // not implemented
  
  ClassDef(AliAnalysisTaskProtonsQA, 1); // example of analysis
};

#endif
