#ifndef ALIANALYSISTASKTRIGGERSTUDY_H
#define ALIANALYSISTASKTRIGGERSTUDY_H

#include "AliAnalysisTaskSE.h"

//-------------------------------------------------------------------------
//                      AliAnalysisTaskTriggerStudy
// 
// 
//
//
// Author: Michele Floris, CERN
//-------------------------------------------------------------------------


class AliESDEvent;
class AliESDtrackCuts;
class AliHistoListWrapper;
class AliTriggerAnalysis;
class AliAnalysisTaskTriggerStudy : public AliAnalysisTaskSE {

public:

  AliAnalysisTaskTriggerStudy();
  AliAnalysisTaskTriggerStudy(const char * name);
  AliAnalysisTaskTriggerStudy(const AliAnalysisTaskTriggerStudy& obj) ;
  ~AliAnalysisTaskTriggerStudy();

  void SetIsMC(Bool_t flag=kTRUE) { fIsMC = flag;}
  AliHistoListWrapper * GetHistoList() { return fHistoList;}

  TH1 * GetHistoTracklets(const char * name, const char * title);

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  

private:

  //
  AliESDEvent *  fESD;    //! ESD object  AliVEvent*     fEvent;
  AliHistoListWrapper  * fHistoList; // wrapper for the list, takes care of merging + histo booking and getters  
  Bool_t fIsMC; // true if processing montecarlo
  
  AliTriggerAnalysis * fTriggerAnalysis;

  AliAnalysisTaskTriggerStudy& operator=(const AliAnalysisTaskTriggerStudy& task);
  
  ClassDef(AliAnalysisTaskTriggerStudy, 2)


};

#endif
