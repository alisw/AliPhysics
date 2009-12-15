#ifndef AliPERFORMANCETASK_H
#define AliPERFORMANCETASK_H

//------------------------------------------------------------------------------
// Task to run reconstruction performance. 
// 
// Author: J.Otwinowski 01/04/2009 
//------------------------------------------------------------------------------

class AliESDEvent;
class AliESDfriend;
class AliMCEvent;
class AliPerformanceObject;
class AliMagF;
class TList;

#include "AliAnalysisTaskSE.h"

class AliPerformanceTask : public AliAnalysisTaskSE {
 public:
  AliPerformanceTask();
  AliPerformanceTask(const char *name, const char *title);
  virtual ~AliPerformanceTask();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  virtual Bool_t Notify();

  // Add comparison objects
  Bool_t AddPerformanceObject(AliPerformanceObject* comp);

  // Use MC
  void SetUseMCInfo(Bool_t useMCInfo = kFALSE) {fUseMCInfo = useMCInfo;}

  // Use ESD friend
  void SetUseESDfriend(Bool_t useESDFriend = kFALSE) {fUseESDfriend = useESDFriend;}

 private:
  AliESDEvent *fESD;   //! ESD event
  AliESDfriend *fESDfriend; //! ESD friend event
  AliMCEvent *fMC;    //! MC event

  TList *fOutput;             //! list send on output slot 0
  TIterator *fPitList;        //! iterator over the output objetcs  
  TList *fCompList;           // list of comparison objects

  Bool_t fUseMCInfo;          // use MC information
  Bool_t fUseESDfriend;       // use ESD friend

  AliPerformanceTask(const AliPerformanceTask&); // not implemented
  AliPerformanceTask& operator=(const AliPerformanceTask&); // not implemented
  
  ClassDef(AliPerformanceTask, 1); // example of analysis
};

#endif
