#ifndef ALIEMCALAODTRACKFILTERTASK_H
#define ALIEMCALAODTRACKFILTERTASK_H

// $Id: AliEmcalAodTrackFilterTask.h 54003 2012-01-19 16:40:42Z loizides $

class TClonesArray;
class AliAODEvent;
class AliAODTrack;

#include "AliAnalysisTaskSE.h"

class AliEmcalAodTrackFilterTask : public AliAnalysisTaskSE {
 public:
  AliEmcalAodTrackFilterTask();
  AliEmcalAodTrackFilterTask(const char *name);
  virtual ~AliEmcalAodTrackFilterTask();

  void UserCreateOutputObjects();
  void UserExec(Option_t *option);
  void SetRunPeriod(const char *p);
  Bool_t AcceptTrack(AliAODTrack *track);
  void RetrieveEventObjects();
  AliAODTrack* GetTrack(const Int_t i) const;
  Int_t GetNumberOfTracks() const;
   
  void SetAODfilterBit(Int_t b)                  { fAODfilterBit     = b;    }
  void SetTracksName(const char *name)           { fTracksOutName    = name; }
  void SetTracksIn(const char *name)             { fTracksInName     = name; }

 protected:
  Int_t              fAODfilterBit;         // if true then do vertex constraint
  TString            fTracksOutName;        // name of output tracks 
  TString            fTracksInName;         // name of input tracks
  AliAODEvent       *fAOD;                  //!aod event
  TClonesArray      *fTracksIn;             //!track array in
  TClonesArray      *fTracksOut;            //!track array out

 private:
  AliEmcalAodTrackFilterTask(const AliEmcalAodTrackFilterTask&);            // not implemented
  AliEmcalAodTrackFilterTask &operator=(const AliEmcalAodTrackFilterTask&); // not implemented

  ClassDef(AliEmcalAodTrackFilterTask, 1); // Class to filter hybrid tracks in AOD events
};

#endif
