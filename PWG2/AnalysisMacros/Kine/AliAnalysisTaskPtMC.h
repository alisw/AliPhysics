#ifndef AliAnalysisTaskPtMC_cxx
#define AliAnalysisTaskPtMC_cxx

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

class TH1F;
class AliESDEvent;

#include "AliAnalysisTask.h"

class AliAnalysisTaskPtMC : public AliAnalysisTask {
 public:
  AliAnalysisTaskPtMC(const char *name = "AliAnalysisTaskPtMC");
  virtual ~AliAnalysisTaskPtMC() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:
  AliESDEvent *fESD;    //ESD object
  TH1F        *fHistPt; //Pt spectrum
   
  AliAnalysisTaskPtMC(const AliAnalysisTaskPtMC&); // not implemented
  AliAnalysisTaskPtMC& operator=(const AliAnalysisTaskPtMC&); // not implemented

  ClassDef(AliAnalysisTaskPtMC, 1); // example of analysis
};

#endif
