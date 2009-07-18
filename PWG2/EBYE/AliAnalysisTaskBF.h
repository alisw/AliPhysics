#ifndef ALIANALYSISTASKBF_CXX
#define ALIANALYSISTASKBF_CXX

// Analysis task for the BF code
// Authors: Panos Cristakoglou@cern.ch

class AliBalance;
class AliESDEvent;

#include "AliAnalysisTask.h"

class AliAnalysisTaskBF : public AliAnalysisTask {
 public:
  AliAnalysisTaskBF(const char *name = "AliAnalysisTaskBF");
  virtual ~AliAnalysisTaskBF() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetAnalysisObject(AliBalance *const analysis) {
    fBalance = analysis;}
  
 private:
  AliESDEvent *fESD;    //ESD object

  AliBalance *fBalance; //BF object
   
  AliAnalysisTaskBF(const AliAnalysisTaskBF&); // not implemented
  AliAnalysisTaskBF& operator=(const AliAnalysisTaskBF&); // not implemented
  
  ClassDef(AliAnalysisTaskBF, 1); // example of analysis
};

#endif
