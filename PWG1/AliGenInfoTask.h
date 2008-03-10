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

class AliGenInfoTask : public AliAnalysisTask {
 public:
 AliGenInfoTask();
 AliGenInfoTask(const char *name);
  virtual ~AliGenInfoTask() {};
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);
  virtual void   SetDebugLevel(Int_t level) {fDebug = level;}
  virtual void   SetMaxTracks(Int_t max=10) {fMaxTracks = max;}
  
 protected:
  AliGenInfoTask(const AliGenInfoTask& /*info*/);
  AliGenInfoTask& operator=(const AliGenInfoTask& /*info*/) { return *this;}

  virtual Int_t FillTrackHistograms(Int_t nTracks, AliESDtrack* track, 
				    AliESDfriendTrack* friendTrack, 
				    AliTPCseed* seed);
  AliGenInfoMaker *fGenMaker;    // gen Maker
  Int_t         fDebug;          //  Debug flag
  AliESDEvent*  fESD;            //! ESD
  AliESDfriend* fESDfriend;      //! ESD friend
  TList*        fListOfHists;    //! Output list of histograms
  
  Int_t         fMaxTracks;      // Max tracks in histogram
  TH1F*         hESDTracks;      //! N ESD tracks
  TH1F*         hGoodTracks;     //! GOOD tracks

  ClassDef(AliGenInfoTask, 1); // Analysis task base class for TPC tracks and clusters
};

#endif
