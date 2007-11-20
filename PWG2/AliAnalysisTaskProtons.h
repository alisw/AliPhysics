#ifndef AliAnalysisTaskProtons_cxx
#define AliAnalysisTaskProtons_cxx

// Analysis task creating a the 2d y-p_t spectrum of p and antip
// Author: Panos Cristakoglou

class TList;
class AliESDEvent;

#include "PWG2spectra/SPECTRA/AliProtonAnalysis.h"
#include "AliAnalysisTask.h"

class AliAnalysisTaskProtons : public AliAnalysisTask {
 public:
  AliAnalysisTaskProtons(const char *name = "AliAnalysisTaskProtons");
  virtual ~AliAnalysisTaskProtons() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:
  AliESDEvent *fESD;    //ESD object
  TList  *fList; //TList output object
  AliProtonAnalysis *fAnalysis; //analysis object
   
  AliAnalysisTaskProtons(const AliAnalysisTaskProtons&); // not implemented
  AliAnalysisTaskProtons& operator=(const AliAnalysisTaskProtons&); // not implemented
  
  ClassDef(AliAnalysisTaskProtons, 1); // example of analysis
};

#endif
