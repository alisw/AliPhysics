/*
 * **********************************************************
 * Virtual class for processing trees of AliReducedEventInfo
 * Authors: Ionut-Cristian Arsene, iarsene@cern.ch, i.c.arsene@fys.uio.no
 *                Jacobus Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
 * Creation date: 2015/10/01
 *********************************************************
 */


#ifndef ALIREDUCEDANALYSISTASKSE_H
#define ALIREDUCEDANALYSISTASKSE_H

#include <TObject.h>
#include <TString.h> 

#include "AliReducedVarManager.h"
#include "AliHistogramManager.h"
#include "AliReducedBaseEvent.h"

//________________________________________________________________
class AliReducedAnalysisTaskSE : public TObject {
  
public:
  AliReducedAnalysisTaskSE();
  AliReducedAnalysisTaskSE(const Char_t* name, const Char_t* title);
  virtual ~AliReducedAnalysisTaskSE();
  
  // Virtual functions, to be implemented in the inheriting analysis classes
  // initialization (typically called in AliAnalysisTask::UserCreateOutputObjects())
  virtual void Init();
  // process a given event (typically called in AliAnalysisTask::UserExec())
  virtual void Process();
  // finish, to be executed after all events were processed
  virtual void Finish();
  // add output objects;
  
  // setters
  void SetEvent(AliReducedBaseEvent* event) {fEvent = event;}
  //void SetHistogramManager(AliHistogramManager* histos) {fHistosManager = histos;}  
  
  // getters
  //AliHistogramManager* GetHistogramManager() const {return fHistosManager;}
  virtual AliHistogramManager* GetHistogramManager() const = 0;  //{return 0x0;}
  AliReducedBaseEvent* GetEvent() const {return fEvent;}
  
protected:
  AliReducedAnalysisTaskSE(const AliReducedAnalysisTaskSE& task);             
  AliReducedAnalysisTaskSE& operator=(const AliReducedAnalysisTaskSE& task);      
  
  //AliHistogramManager* fHistosManager;   // Histogram manager
  
  TString fName;             // name
  TString fTitle;                // title
    
  AliReducedBaseEvent* fEvent;           //! current event to be processed
  Float_t fValues[AliReducedVarManager::kNVars];   // array of values to hold information for histograms
  
  ClassDef(AliReducedAnalysisTaskSE, 2)
};

#endif
