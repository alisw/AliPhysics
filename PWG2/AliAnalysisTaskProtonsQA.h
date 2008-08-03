#ifndef AliAnalysisTaskProtonsQA_cxx
#define AliAnalysisTaskProtonsQA_cxx

// Analysis task creating a the 2d y-p_t spectrum of p and antip
// Author: Panos Cristakoglou
class TString;
class TList;
class AliESDEvent;
class AliAODEvent;
class AliMCEvent;
class AliProtonAnalysis;
class TF1;

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

  TList  *fList; //TList output object                                                                                
  
  AliProtonAnalysis *fAnalysis; //analysis object                                                                     
  
  AliAnalysisTaskProtonsQA(const AliAnalysisTaskProtonsQA&); // not implemented                                           
  AliAnalysisTaskProtonsQA& operator=(const AliAnalysisTaskProtonsQA&); // not implemented                                 
  
  ClassDef(AliAnalysisTaskProtonsQA, 1); // example of analysis
};

#endif
