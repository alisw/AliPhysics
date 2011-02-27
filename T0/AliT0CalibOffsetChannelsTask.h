#ifndef AliT0CalibOffsetChannelsTask_cxx
#define AliT0CalibOffsetChannelsTask_cxx

// task determines mean and sigma of T0 signals  ORA, ORC, ORA-ORC, ORA+ORC/2  
// Authors: FK  

class TH1F;
class AliESDEvent;
class AliT0CalibSeasonTimeShift;

#include "AliAnalysisTaskSE.h"

class AliT0CalibOffsetChannelsTask : public AliAnalysisTaskSE {
 public:
  AliT0CalibOffsetChannelsTask();
  AliT0CalibOffsetChannelsTask(const char *name);
  virtual ~AliT0CalibOffsetChannelsTask(); 
  
  virtual void ConnectInputData(Option_t *option);
  virtual void   UserCreateOutputObjects();
  //  virtual void     Process(AliESDEvent *event);
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  TObjArray* GetOffsetHistos() {return fTzeroObject;}
  
 private:
  AliESDEvent *fESD;        //! ESD object
  TObjArray *fTzeroObject;  //array with CFDi-CFD1 and  CFDi
  TH1F        *fTimeDiff[24];   //! CFDi-CFD1 vs Npmt   
  TH1F        *fCFD[24];   //! CFDi  vs Npmt 
   int         fRunNumber;
  
   //  AliT0CalibSeasonTimeShift *fTzeroObject;
 
  AliT0CalibOffsetChannelsTask(const AliT0CalibOffsetChannelsTask&); // not implemented
  AliT0CalibOffsetChannelsTask& operator=(const AliT0CalibOffsetChannelsTask&); // not implemented
  
  ClassDef(AliT0CalibOffsetChannelsTask, 1); // example of analysis
};

#endif
