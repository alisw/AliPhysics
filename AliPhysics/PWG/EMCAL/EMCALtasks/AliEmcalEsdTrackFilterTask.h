#ifndef ALIEMCALESDTRACKFILTERTASK_H
#define ALIEMCALESDTRACKFILTERTASK_H

class TClonesArray;
class AliESDEvent;
class AliESDtrackCuts;

#include <TF1.h>

#include "AliAnalysisTaskSE.h"
#include "AliESDtrackCuts.h"

class AliEmcalEsdTrackFilterTask : public AliAnalysisTaskSE {
 public:
  AliEmcalEsdTrackFilterTask();
  AliEmcalEsdTrackFilterTask(const char *name);
  virtual ~AliEmcalEsdTrackFilterTask();

  void SetDist(Double_t d)                       { fDist             = d;    }
  void SetDoPropagation(Bool_t b)                { fDoPropagation    = b;    }
  void SetDoSpdVtxConstrain(Bool_t b)            { fDoSpdVtxCon      = b;    }
  void SetHybridTrackCuts(AliESDtrackCuts *cuts) { fHybridTrackCuts  = cuts; }
  void SetIncludeNoITS(Bool_t f)                 { fIncludeNoITS     = f;    }
  void SetTrackCuts(AliESDtrackCuts *cuts)       { fEsdTrackCuts     = cuts; }
  void SetTracksName(const char *name)           { fTracksName       = name; }
  void SetTrackEfficiency(Double_t eff = 0.95)   { fTrackEfficiency  = new TF1("eff", "[0]", 0, 500); fTrackEfficiency->FixParameter(0,eff); }
  void SetTrackEfficiency(TF1* eff)              { fTrackEfficiency  = eff;  }
  void SetMC(Bool_t b)                           { fIsMC             = b  ;  }

 protected:
  void UserCreateOutputObjects();
  void UserExec(Option_t *option);

  AliESDtrackCuts   *fEsdTrackCuts;      // esd track cuts
  Bool_t             fDoSpdVtxCon;       // if true then do vertex constraint
  AliESDtrackCuts   *fHybridTrackCuts;   // hybrid track cuts
  TString            fTracksName;        // name of tracks 
  Bool_t             fIncludeNoITS;      // includes tracks with failed ITS refit
  Bool_t             fDoPropagation;     // propagate all hybrid tracks to EMCal surface
  Double_t           fDist;              // distance to surface (440cm default)
  TF1               *fTrackEfficiency;   // track efficiency
  Bool_t             fIsMC;              // whether it is a MC event or not
  AliESDEvent       *fEsdEv;             //!esd event
  TClonesArray      *fTracks;            //!track array

 private:
  AliEmcalEsdTrackFilterTask(const AliEmcalEsdTrackFilterTask&);            // not implemented
  AliEmcalEsdTrackFilterTask &operator=(const AliEmcalEsdTrackFilterTask&); // not implemented

  ClassDef(AliEmcalEsdTrackFilterTask, 4); // Class to constrain TPC tracks to SPD vertex
};

#endif
