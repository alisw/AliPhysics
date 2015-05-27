#ifndef ALIGLOALGTASK_H
#define ALIGLOALGTASK_H

#include "AliAnalysisTaskSE.h"
#include <TString.h>
class AliAlgSteer;

class AliGloAlgTask : public AliAnalysisTaskSE {
 public:
  //
  AliGloAlgTask(const char *name = "AliGloAlgTask");
  virtual ~AliGloAlgTask(); 
  
  virtual void  LocalInit();
  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t *option);
  virtual void  NotifyRun();
  virtual void  Terminate(Option_t *);
  //
  UInt_t     GetTriggerSelection()                            const {return fTrigSel;}
  void       SetTriggerSelection(UInt_t sel=AliVEvent::kAny)        {fTrigSel = sel;}
  //
  void        SetConfMacroName(const char* nm="alignConf.C")        {fConfMacroName = nm;}
  const char* GetConfMacroName()                              const {return fConfMacroName.Data();}
  //
  void        SetIniParFileName(const char* nm=0)                   {fIniParFileName = nm;}
  const char* GetIniParFileName()                              const {return fIniParFileName.Data();}
  //
 protected:
  //
  TList*       fOutput;                   // output list send on output slot 1
  UInt_t       fTrigSel;                  // trigger selection
  AliAlgSteer* fAlgSteer;                 // alignment steering
  //
  TString      fIniParFileName;          // initial parameters file
  TString      fConfMacroName;           // name of alignment configuration macro
  //
 private:    
  AliGloAlgTask(const AliGloAlgTask&); // not implemented
  AliGloAlgTask& operator=(const AliGloAlgTask&); // not implemented 
  //  
  ClassDef(AliGloAlgTask, 1);  
};


#endif
