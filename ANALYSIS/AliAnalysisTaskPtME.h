// gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include");
#ifndef AliAnalysisTaskPt_cxx
#define AliAnalysisTaskPt_cxx

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

class TH1F;
class AliESDEvent;

#include "AliAnalysisTaskME.h"

class AliAnalysisTaskPtME : public AliAnalysisTaskME {
 public:
  AliAnalysisTaskPtME(const char *name = "AliAnalysisTaskPt");
  virtual ~AliAnalysisTaskPtME() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:
  TH1F        *fHistPt; //Pt spectrum
   
  AliAnalysisTaskPtME(const AliAnalysisTaskPtME&); // not implemented
  AliAnalysisTaskPtME& operator=(const AliAnalysisTaskPtME&); // not implemented
  
  ClassDef(AliAnalysisTaskPtME, 1); // example of analysis
};

#endif
