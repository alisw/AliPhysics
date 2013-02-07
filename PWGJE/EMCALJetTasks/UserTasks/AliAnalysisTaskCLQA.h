#ifndef ALIANALYSISTASKCLQA_H
#define ALIANALYSISTASKCLQA_H

// $Id: AliAnalysisTaskCLQA.h 58847 2012-09-30 18:11:49Z loizides $

class TClonesArray;
class TString;
class TH1F;
class TH2F;
class TH3F;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskCLQA : public AliAnalysisTaskEmcalJet {
 public:
  AliAnalysisTaskCLQA();
  AliAnalysisTaskCLQA(const char *name);
  virtual ~AliAnalysisTaskCLQA();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

 protected:

  Bool_t                      FillHistograms()                                              ;
  Bool_t                      RetrieveEventObjects()                                        ;

 private:
  AliAnalysisTaskCLQA(const AliAnalysisTaskCLQA&);            // not implemented
  AliAnalysisTaskCLQA &operator=(const AliAnalysisTaskCLQA&); // not implemented

  ClassDef(AliAnalysisTaskCLQA, 1) // Constantin's Task
};
#endif
