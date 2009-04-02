#ifndef AliAnalysisTaskProtons_cxx
#define AliAnalysisTaskProtons_cxx

// Analysis task to run the \bar{p}/p analysis
// Author: Panos Cristakoglou
//class TString;
class TList;
class AliESDEvent;
class AliAODEvent;
class AliMCEvent;
class AliProtonAnalysis;

#include "AliAnalysisTask.h"

class AliAnalysisTaskProtons : public AliAnalysisTask {
 public:
  AliAnalysisTaskProtons();
  AliAnalysisTaskProtons(const char *name);
  virtual ~AliAnalysisTaskProtons() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetAnalysisObject(AliProtonAnalysis *analysis) {
    fProtonAnalysis = analysis;}
  
 private:
  AliESDEvent *fESD;    //ESD object 
  AliAODEvent *fAOD;    //AOD object
  AliMCEvent  *fMC;     //MC object 
  
  TList  *fList; //TList output object 
  
  AliProtonAnalysis *fProtonAnalysis; //analysis object 
  
  AliAnalysisTaskProtons(const AliAnalysisTaskProtons&); // not implemented
  AliAnalysisTaskProtons& operator=(const AliAnalysisTaskProtons&); // not implemented
  
  ClassDef(AliAnalysisTaskProtons, 1); // example of analysis
};

#endif
