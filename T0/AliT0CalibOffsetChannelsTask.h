#ifndef AliT0CalibOffsetChannelsTask_cxx
#define AliT0CalibOffsetChannelsTask_cxx

// task determines mean and sigma of T0 signals  ORA, ORC, ORA-ORC, ORA+ORC/2  
// Authors: FK  

class TH1F;
class TObjArray; 
class AliESDEvent;

#include "AliAnalysisTaskSE.h"

class AliT0CalibOffsetChannelsTask : public AliAnalysisTaskSE {
 public:
  AliT0CalibOffsetChannelsTask();
  AliT0CalibOffsetChannelsTask(const char *name);
  virtual ~AliT0CalibOffsetChannelsTask(); 
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  TObjArray* GetOffsetHistos() {return fTzeroObject;}
  
 private:
  AliESDEvent *fESD;          //! ESD object
  TObjArray   *fTzeroObject;  // array with CFDi-CFD1 and  CFDi
  TH1F        *fTimeDiff[24]; //! CFDi-CFD1 vs Npmt   
  TH1F        *fCFD[24];      //! CFDi  vs Npmt 
  TH1F        *fTzeroORA;     //! or A spectrum    
  TH1F        *fTzeroORC;     //! or C spectrum    
  TH1F        *fResolution;   //! or A minus or C spectrum    
  TH1F        *fTzeroORAplusORC; //! ORA+ORC /2 
  int         fRunNumber;
  
 
  AliT0CalibOffsetChannelsTask(const AliT0CalibOffsetChannelsTask&); // not implemented
  AliT0CalibOffsetChannelsTask& operator=(const AliT0CalibOffsetChannelsTask&); // not implemented
  
  ClassDef(AliT0CalibOffsetChannelsTask, 1); // example of analysis
};

#endif
