#ifndef ALIDNDPTTASK_H
#define ALIDNDPTTASK_H

//------------------------------------------------------------------------------
// Task for dNdPt analysis.
// 
// Author: J.Otwinowski 04/11/2008 
//------------------------------------------------------------------------------

class AliESDEvent;
class AliMCEvent;
class AlidNdPtEventCuts;
class AlidNdPtAcceptanceCuts;
class AliESDtrackCuts;
class AlidNdPt;
class AlidNdPtAnalysis;
class AlidNdPtCorrection;
class AliMagFMaps;
class TList;

#include "dNdPt/AlidNdPtHelper.h"
#include "AliAnalysisTaskSE.h"

class AlidNdPtTask : public AliAnalysisTaskSE {
 public:
  AlidNdPtTask(const char *name = "AlidNdPtTask");
  virtual ~AlidNdPtTask();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  virtual Bool_t Notify();

  Bool_t AddAnalysisObject(AlidNdPt *pObj);
  void SetUseMCInfo(Bool_t info)           { fUseMCInfo = info; }
  
 private:

  AliESDEvent *fESD;    //! ESD event
  AliMCEvent *fMC;      //! MC event
  TList* fOutput;       //! list send on output slot 0
  TIterator *fPitList;  //! iterator over the output objetcs  
  TList *fCompList;     // list of comparison objects

  Bool_t fUseMCInfo;        // use MC information

  AlidNdPtTask(const AlidNdPtTask&); // not implemented
  AlidNdPtTask& operator=(const AlidNdPtTask&); // not implemented
  
  ClassDef(AlidNdPtTask, 3); // example of analysis
};

#endif
