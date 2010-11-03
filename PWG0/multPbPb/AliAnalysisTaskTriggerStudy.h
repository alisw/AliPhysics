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

  TH1 * GetHistoTracklets   (const char * name, const char * title);
  void  FillTriggerOverlaps (const char * name, const char * title, Int_t nFastOrOffline, Bool_t v0A, Bool_t v0C, Bool_t OM2, 
			     Bool_t OM3, Bool_t cMBS2A,Bool_t cMBS2C, Bool_t cMBAC) ;

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  

private:

  //
  AliESDEvent *  fESD;    //! ESD object  AliVEvent*     fEvent;
  AliHistoListWrapper  * fHistoList; // wrapper for the list, takes care of merging + histo booking and getters  
  Bool_t fIsMC; // true if processing montecarlo
  
  AliTriggerAnalysis * fTriggerAnalysis; // trigger analysis object, to get the offline triggers
  TString fHistoSuffix; // suffix appended to all histos, set in the user exec.


  AliAnalysisTaskTriggerStudy& operator=(const AliAnalysisTaskTriggerStudy& task);
  
  ClassDef(AliAnalysisTaskTriggerStudy, 2)


};

#endif
