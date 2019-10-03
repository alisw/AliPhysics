#ifndef AliFakeTrackTask_H
#define  AliFakeTrackTask_H
//#include <fstream>
class TH1F;
class TH2F;
class TH3F;
class AliESDEvent;
class AliESDtrack;
class  AliESDtrackCuts;
//class AliESDpidCuts;
class AliPIDResponse ; 
class AliESDpid;
class TGraph;
class AliStack;
class TChain;
#include "AliAnalysisTaskSE.h"
//#include "AliESDtrackCuts.h"
#include "AliESDpidCuts.h"
#include "AliESDpid.h"




class AliFakeTrackTask : public AliAnalysisTaskSE {
 public:
  AliFakeTrackTask(const char *name = "AliFakeTrackTask");
  virtual ~AliFakeTrackTask();
  void SetTrackCuts(AliESDtrackCuts* cuts){ftrackcuts=cuts;} 	
  
  //virtual void   ConnectInputData(Option_t *);
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *){}; 
  //virtual void   LocalInit();
 
 
 
 private:
 
 
 AliESDEvent *fESD;    //ESD object    
 TH3F* fptvsTPCsignalvsITSsignalAll; // pt vs. TPC signal vs. ITS signal
 TH3F* fptvsTPCsignalvsITSsignalGlobalgood; // 	good tracks pt vs. TPC signal vs. ITS signal
 TH3F* fptvsTPCsignalvsITSsignalGlobalfake; // global fake pt vs. TPC signal vs. ITS signal
 TH3F* fptvsTPCsignalvsITSsignalTPCfake; // global fake pt vs. TPC signal vs. ITS signal
 TH3F* fptvsTPCsignalvsITSsignalITSfake; // global fake pt vs. TPC signal vs. ITS signal
 TH1F* ffakestat; // types of the fake tracks 	

 AliESDtrackCuts* ftrackcuts; // track cuts for track selection 	 	
 TList* flistout; //out list
 AliPIDResponse   *fPIDResponse; // pointer to pid response   

 ClassDef(AliFakeTrackTask, 1); // example of analysis
};

#endif
