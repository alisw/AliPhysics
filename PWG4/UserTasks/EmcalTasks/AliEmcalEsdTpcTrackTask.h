#ifndef ALIEMCALESDTPCTRACKTASK_H
#define ALIEMCALESDTPCTRACKTASK_H

// $Id$

class TClonesArray;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"

class AliEmcalEsdTpcTrackTask : public AliAnalysisTaskSE {
 public:
  AliEmcalEsdTpcTrackTask();
  AliEmcalEsdTpcTrackTask(const char *name);
  virtual ~AliEmcalEsdTpcTrackTask();

  void UserCreateOutputObjects();
  void UserExec(Option_t *option);
   
  void SetDoSpdVtxConstrain(Bool_t b)            { fDoSpdVtxCon      = b;    }
  void SetHybridTrackCuts(AliESDtrackCuts *cuts) { fHybridTrackCuts  = cuts; }
  void SetTrackCuts(AliESDtrackCuts *cuts)       { fEsdTrackCuts     = cuts; }
  void SetTracksName(const char *name)           { fTracksName       = name; }

 protected:
  AliESDtrackCuts   *fEsdTrackCuts;      // esd track cuts
  Bool_t             fDoSpdVtxCon;       // if true then do vertex constraint
  AliESDtrackCuts   *fHybridTrackCuts;   // hybrid track cuts
  TString            fTracksName;        // name of tracks 
  AliESDEvent       *fEsdEv;             //!esd event
  TClonesArray      *fTracks;            //!track array

 private:
  AliEmcalEsdTpcTrackTask(const AliEmcalEsdTpcTrackTask&);            // not implemented
  AliEmcalEsdTpcTrackTask &operator=(const AliEmcalEsdTpcTrackTask&); // not implemented

  ClassDef(AliEmcalEsdTpcTrackTask, 1); // Class to constrain TPC tracks to SPD vertex
};

#endif
