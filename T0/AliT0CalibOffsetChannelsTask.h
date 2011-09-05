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

  Float_t  GetCFDvalue(Int_t channel)  const {return fCDBcfds[channel];}
  Float_t* GetCFDvalue()          const {return (float*)fCDBcfds;}
  Float_t  GetTimeEq(Int_t channel)        const {return fCDBdelays[channel];}
  Float_t* GetTimeEq()          const {return (float*)fCDBdelays;}
  void SetCFDvalue(Int_t channel, Float_t val) {fCDBcfds[channel]=val;}
  void SetTimeEq(Int_t channel, Float_t val) {fCDBdelays[channel]=val;}
  Float_t *GetT0Means() { return fCDBT0s;}
  void SetT0Means(Int_t ihist, Float_t mean ) {fCDBT0s[ihist]=mean;};

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
  Float_t fCDBdelays[24];  //time delays from OCDB
  Float_t fCDBcfds[24];    // mean CFD from OCDB
  Float_t fCDBT0s[4];      //position T0AC, T0A, T0A, resolution

  
 
  AliT0CalibOffsetChannelsTask(const AliT0CalibOffsetChannelsTask&); // not implemented
  AliT0CalibOffsetChannelsTask& operator=(const AliT0CalibOffsetChannelsTask&); // not implemented
  
  ClassDef(AliT0CalibOffsetChannelsTask, 2); // example of analysis
};

#endif
