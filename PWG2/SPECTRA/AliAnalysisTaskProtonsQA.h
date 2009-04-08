#ifndef ALIANALYSISTASKPROTONSQA_H
#define ALIANALYSISTASKPROTONSQA_H

//-------------------------------------------------------------------------
//                 Class AliAnalysisTaskProtonQA
//   This is the task for the baryon (proton) QA analysis
//
//    Origin: Panos Christakoglou | Panos.Christakoglou@cern.ch
//-------------------------------------------------------------------------

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

  void SetAnalysisObject(AliProtonQAAnalysis *const analysis) {
    fProtonQAAnalysis = analysis;}

 
 private:
  AliESDEvent *fESD;    //ESD object
  AliMCEvent  *fMC;     //MC object

  TList  *fList0; //TList output object
  TList  *fList1; //TList output object
  TList  *fList2; //TList output object
  TList  *fList3; //TList output object
  TList  *fList4; //TList output object
  TList  *fList5; //TList output object
  TList  *fList6; //TList output object
  TList  *fList7; //TList output object
  
  AliProtonQAAnalysis *fProtonQAAnalysis; //analysis object
 
  AliAnalysisTaskProtonsQA(const AliAnalysisTaskProtonsQA&); // not implemented
  AliAnalysisTaskProtonsQA& operator=(const AliAnalysisTaskProtonsQA&); // not implemented
  
  ClassDef(AliAnalysisTaskProtonsQA, 1); // example of analysis
};

#endif
