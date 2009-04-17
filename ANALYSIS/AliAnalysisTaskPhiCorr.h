// gSystem->SetIncludePath("-I$ROOTSYS/include -I$ALICE_ROOT/include");
#ifndef AliAnalysisTaskPt_cxx
#define AliAnalysisTaskPt_cxx

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

class TH1F;
class TList;
class AliESDEvent;

#include "AliAnalysisTaskME.h"
#include "AliMixedEvent.h"

class AliAnalysisTaskPhiCorr : public AliAnalysisTaskME {
 public:
  AliAnalysisTaskPhiCorr(const char *name = "AliAnalysisTaskPt");
  virtual ~AliAnalysisTaskPhiCorr() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
 private:
  TList       *fHists;      // List with histos
  TH1F        *fHistDphiCO; // Pt spectrum
  TH1F        *fHistDphiUC; // Pt spectrum
  AliMixedEvent fMixedEvent;
  
  AliAnalysisTaskPhiCorr(const AliAnalysisTaskPhiCorr&); // not implemented
  AliAnalysisTaskPhiCorr& operator=(const AliAnalysisTaskPhiCorr&); // not implemented
  
  ClassDef(AliAnalysisTaskPhiCorr, 1); // example of analysis
};

#endif
