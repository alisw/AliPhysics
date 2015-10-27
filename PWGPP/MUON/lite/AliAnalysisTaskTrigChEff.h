#ifndef ALIANALYSISTASKTRIGCHEFF_H
#define ALIANALYSISTASKTRIGCHEFF_H

/* $Id$ */ 

//
// Class for trigger chamber efficiency calculations
// and tests
//
// Author: Diego Stocco
//

#include "AliVAnalysisMuon.h"

class AliMuonTrackCuts;
class TList;
class TObjArray;
class TString;
class AliTrigChEffOutput;

class AliAnalysisTaskTrigChEff : public AliVAnalysisMuon {
 public:
  AliAnalysisTaskTrigChEff();
  AliAnalysisTaskTrigChEff(const char *name, const AliMuonTrackCuts& cuts);
  virtual ~AliAnalysisTaskTrigChEff();

  void Terminate(Option_t *option);
  void FinishTaskOutput();

  void MyUserCreateOutputObjects();
  void ProcessEvent(TString physSel, const TObjArray& selectTrigClasses, TString centrality);

 private:

  AliAnalysisTaskTrigChEff(const AliAnalysisTaskTrigChEff&);
  AliAnalysisTaskTrigChEff& operator=(const AliAnalysisTaskTrigChEff&);

  AliTrigChEffOutput* fAnalysisOutput; //!<! Output handler
  TList*  fList;     //!<! TList output object

  ClassDef(AliAnalysisTaskTrigChEff, 5); // Trigger chamber efficiencies
};

#endif
