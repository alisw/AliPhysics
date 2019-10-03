#ifndef ALITPCOMPARISONPID_H
#define ALITPCOMPARISONPID_H

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
class AliGenInfoMaker;
class TTreeSRedirector;
class AliMCEventHadnler;
class TParticle;
class AliMCInfo;
class AliESDRecInfo;
class AliESDEvent;
class AliMCEvent;
class THnSparse;

class AliTPCComparisonPID : public AliAnalysisTask {
public:
  AliTPCComparisonPID();
  AliTPCComparisonPID(const char *name);
  AliTPCComparisonPID(const AliTPCComparisonPID& info);
  virtual ~AliTPCComparisonPID();  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  virtual void   FinishTaskOutput();
  void           SetDebugOuputhPath(const char * name){fDebugOutputPath=name;}
  void           Init();
  //
  void           ProcessMCInfo();
  //
  THnSparse * GetTPCsignal(){return fTPCsignal;}
  THnSparse * GetTPCsignalNorm(){return fTPCsignalNorm;}
  //
  // debug streamer part
  //
  TTreeSRedirector *GetDebugStreamer();
  void       SetStreamLevel(Int_t streamLevel){fStreamLevel=streamLevel;}
  void       SetDebugLevel(Int_t level) {fDebugLevel = level;}
  Int_t      GetStreamLevel() const {return fStreamLevel;}
  Int_t      GetDebugLevel() const {return fDebugLevel;}
  //
protected:
  void RegisterDebugOutput();
  AliTPCComparisonPID& operator=(const AliTPCComparisonPID& /*info*/) { return *this;}
  AliMCEvent  * fMCinfo;          //! MC event handler
  AliESDEvent * fESD;             //! current esd event
  //
  //
  //
  THnSparse * fTPCsignal;         //raw tpc signal - dEdx
  THnSparse * fTPCsignalNorm;     //normalized TPC signal
  //
  TTreeSRedirector *fDebugStreamer;     //! debug streamer
  Int_t  fStreamLevel;                  //  debug stream level 
  Int_t  fDebugLevel;                   //  debug level
  TString      fDebugOutputPath; // debug output path
  ClassDef(AliTPCComparisonPID, 1); // Analysis task base class for tracks
};

#endif
