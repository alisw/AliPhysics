#ifndef ALIMCTRACKINGTESTTASK_H
#define ALIMCTRACKINGTESTTASK_H

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
class AliComparisonObject;
class AliTrackComparison;

class AliMCTrackingTestTask : public AliAnalysisTask {
 public:
 AliMCTrackingTestTask();
 AliMCTrackingTestTask(const char *name);
  virtual ~AliMCTrackingTestTask();
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  virtual void FinishTaskOutput();
  void         SetDebugOuputhPath(const char * name){fDebugOutputPath=name;}

  //
  void           ProcessMCInfo();
  void           ProcessRefTracker(AliTrackReference* refIn, AliTrackReference* refOut, TParticle*part, Int_t type);
  
  void           FitTrackRefs(TParticle * part, TClonesArray * trefs);

  Bool_t         IsFindable(Int_t label, Float_t minTrackLength);
  Bool_t         AddComparisonObject(AliTrackComparison* comp);
  //
  // debug streamer part
  //
  TTreeSRedirector *GetDebugStreamer();
  void       SetStreamLevel(Int_t streamLevel){fStreamLevel=streamLevel;}
  void       SetDebugLevel(Int_t level) {fDebugLevel = level;}
  Int_t      GetStreamLevel() const {return fStreamLevel;}
  Int_t      GetDebugLevel() const {return fDebugLevel;}
  //
  static AliExternalTrackParam * MakeTrack(const AliTrackReference* ref, TParticle*part);
  static Bool_t  PropagateToPoint(AliExternalTrackParam *param, Double_t *xyz, Double_t mass,  Float_t step);
 protected:
  void RegisterDebugOutput();
  AliMCTrackingTestTask(const AliMCTrackingTestTask& /*info*/);
  AliMCTrackingTestTask& operator=(const AliMCTrackingTestTask& /*info*/) { return *this;}
  AliMCEvent  * fMCinfo;          //! MC event handler
  AliESDEvent * fESD;             //! current esd event
  
  Int_t         fCurrentRun;      //  current run number

  //
  //
  //
  TTreeSRedirector *fDebugStreamer;     //! debug streamer
  Int_t  fStreamLevel;                  //  debug stream level 
  Int_t  fDebugLevel;                   //  debug level
  TString      fDebugOutputPath; // debug output path

  TList* fOutList;
  TIterator *fPitList;        //! iterator over the output objetcs  
  TList *fCompList; 

  ClassDef(AliMCTrackingTestTask, 1); // Analysis task base class for tracks
};

#endif
