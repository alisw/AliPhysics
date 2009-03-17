#ifndef ALIGENINFOTASK_H
#define ALIGENINFOTASK_H

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
class AliGenInfoMaker;
class TTreeSRedirector;
class AliMCEventHadnler;
class TParticle;
class AliMCInfo;
class AliESDRecInfo;
class AliESDEvent;
class AliESDfriend;
class AliMCEvent;
class AliComparisonObject;

class AliGenInfoTask : public AliAnalysisTask {
 public:
 AliGenInfoTask();
 AliGenInfoTask(const char *name);
  virtual ~AliGenInfoTask();
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  virtual void   FinishTaskOutput();
  //
  //
  void ProcessMCInfo();
  void ProcessESDInfo();
  void ProcessComparison();
  void DumpInfo();
  //
  //
  // debug streamer part
  //
  TTreeSRedirector *GetDebugStreamer();
  void       SetStreamLevel(Int_t streamLevel){fStreamLevel=streamLevel;}
  void       SetDebugLevel(Int_t level) {fDebugLevel = level;}
  Int_t      GetStreamLevel() const {return fStreamLevel;}
  Int_t      GetDebugLevel() const {return fDebugLevel;}
  //
  Bool_t     AcceptParticle(TParticle *part);  
  AliMCInfo *GetTrack(Int_t index, Bool_t force=kFALSE);
  AliESDRecInfo *GetRecTrack(Int_t index, Bool_t force=kFALSE);
  Bool_t     AddComparisonObject(AliComparisonObject *pObj);
  void             RegisterDebugOutput(const char *path);
  void         SetDebugOuputhPath(const char * name){fDebugOutputPath=name;}
 protected:
  AliGenInfoTask(const AliGenInfoTask& /*info*/);
  AliGenInfoTask& operator=(const AliGenInfoTask& /*info*/) { return *this;}
  AliMCEvent  * fMCinfo;     //! MC event handler
  AliESDEvent * fESD;             //! current esd event
  AliESDfriend * fESDfriend;             //! current esd event
  //
  TObjArray   *fCompList;        // comparison object list
  //
  TClonesArray *fGenTracksArray;  //clones array with filtered particles
  TClonesArray *fGenKinkArray;    //clones array with filtered Kinks
  TClonesArray *fGenV0Array;      //clones array with filtered V0s
  //
  TClonesArray *fRecTracksArray;  //clones array with filtered tracks 
  //
  //
  TTreeSRedirector *fDebugStreamer;     //! debug streamer
  Int_t  fStreamLevel;                  //  debug stream level 
  Int_t  fDebugLevel;                   //  debug level
  TString      fDebugOutputPath; // debug output path
  ClassDef(AliGenInfoTask, 1); // Analysis task base class for tracks
};

#endif
