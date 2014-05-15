#ifndef ALIEMCALTRACKPROPAGATORTASK_H
#define ALIEMCALTRACKPROPAGATORTASK_H

// $Id: AliEmcalTrackPropagatorTask.h | Mon Dec 9 12:59:28 2013 +0100 | Constantin Loizides  $

class TClonesArray;

#include "AliAnalysisTaskSE.h"

class AliEmcalTrackPropagatorTask : public AliAnalysisTaskSE {
 public:
  AliEmcalTrackPropagatorTask();
  AliEmcalTrackPropagatorTask(const char *name);
  virtual ~AliEmcalTrackPropagatorTask();

  void               SetDist(Double_t d)                 { fDist               = d;    }
  void               SetOnlyIfNotSet(Bool_t b)           { fOnlyIfNotSet       = b; }
  void               SetTracksInName(const char *n)      { fTracksInName       = n; }
  void               SetTracksOutName(const char *n)     { fTracksOutName      = n; }

 protected:
  void               UserCreateOutputObjects();
  void               UserExec(Option_t *option);
   
  TString            fTracksInName;      // name of tracks in  
  TString            fTracksOutName;     // name of tracks out
  Double_t           fDist;              // distance to surface (440cm default)
  Bool_t             fOnlyIfNotSet;      // only attempt if not already at surface
  TClonesArray      *fTracksIn;          //!track array in
  TClonesArray      *fTracksOut;         //!track array out

 private:
  AliEmcalTrackPropagatorTask(const AliEmcalTrackPropagatorTask&);            // not implemented
  AliEmcalTrackPropagatorTask &operator=(const AliEmcalTrackPropagatorTask&); // not implemented

  ClassDef(AliEmcalTrackPropagatorTask, 2); // Class to propagate and store track parameters at EMCAL surface
};
#endif
