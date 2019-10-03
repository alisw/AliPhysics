#ifndef ALITPCTASKPID_H
#define ALITPCTASKPID_H

// ROOT includes
#include <TList.h>
#include <TH1.h>
#include <TH2.h>

// AliRoot includes
#include <AliAnalysisTask.h>
#include <AliESDEvent.h>
#include <AliESDfriend.h>
#include <AliESDtrack.h>
#include <AliESDfriendTrack.h>
#include <AliTPCseed.h>
#include <TString.h>
class AliMCEventHadnler;
class TParticle;
class AliMCInfo;
class AliESDRecInfo;
class AliESDEvent;
class AliMCEvent;
class THnSparse;
class TObjArray;

class AliTPCtaskPID : public AliAnalysisTask {
public:
  AliTPCtaskPID();
  AliTPCtaskPID(const char *name);
  AliTPCtaskPID(const AliTPCtaskPID& info);
  virtual ~AliTPCtaskPID();  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  virtual void   FinishTaskOutput();
  void           Init();
  //
  void           ProcessMCInfo();
  //
  THnSparse * GetTPCsignal(){return fTPCsignal;}
  THnSparse * GetTPCsignalNorm(){return fTPCsignalNorm;}
  THnSparse * GetTPCr(){return fTPCr;}
  //
  static void BinLogX(TAxis *axis);
protected:
  void RegisterDebugOutput();
  AliTPCtaskPID& operator=(const AliTPCtaskPID& /*info*/) { return *this;}
  AliMCEvent  * fMCinfo;          //! MC event handler
  AliESDEvent * fESD;             //! current esd event
  //
  //
  // 
  TObjArray  *fList; //TList output object
  THnSparse * fTPCsignal;         //raw tpc signal - dEdx
  THnSparse * fTPCsignalNorm;     //normalized TPC signal
  THnSparse * fTPCr;              //TPC pid info from ESD
  ClassDef(AliTPCtaskPID, 1); // Analysis task base class for tracks
};

#endif
