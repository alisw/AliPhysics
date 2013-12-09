#ifndef ALIEMCALESDTPCTRACKTASK_H
#define ALIEMCALESDTPCTRACKTASK_H

// $Id$

class TClonesArray;
class AliESDEvent;
class AliESDtrack;
class AliESDtrackCuts;
class AliEMCALRecoUtils;

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
  void SetIncludeNoITS(Bool_t f)                 { fIncludeNoITS     = f;    }

  void SetDoPropagation(Bool_t b)                { fDoPropagation    = b;    }
  void SetDist(Double_t d)                       { fDist             = d;    }
  void SetMinPtProp(Double_t pt)                 { fMinPtCutProp     = pt;   }
  void SetRecoUtils(AliEMCALRecoUtils *ru)       { fRecoUtils        = ru;   }

 protected:
  void PropagateTrackToEMCal(AliESDtrack *esdTrack);

  AliESDtrackCuts   *fEsdTrackCuts;      // esd track cuts
  Bool_t             fDoSpdVtxCon;       // if true then do vertex constraint
  AliESDtrackCuts   *fHybridTrackCuts;   // hybrid track cuts
  TString            fTracksName;        // name of tracks 
  Bool_t             fIncludeNoITS;      // includes tracks with failed ITS refit
  Bool_t             fDoPropagation;     // propagate all hybrid tracks to EMCal surface
  AliEMCALRecoUtils *fRecoUtils;         // esd reco utils
  Double_t           fDist;              // distance to surface (430cm default)
  Double_t           fMinPtCutProp;      // minimum track pt cut for propagated tracks (350 MeV/c default)
  AliESDEvent       *fEsdEv;             //!esd event
  TClonesArray      *fTracks;            //!track array

 private:
  AliEmcalEsdTpcTrackTask(const AliEmcalEsdTpcTrackTask&);            // not implemented
  AliEmcalEsdTpcTrackTask &operator=(const AliEmcalEsdTpcTrackTask&); // not implemented

  ClassDef(AliEmcalEsdTpcTrackTask, 2); // Class to constrain TPC tracks to SPD vertex
};

#endif
