#ifndef ALIEMCALMCTRACKSELECTOR_H
#define ALIEMCALMCTRACKSELECTOR_H

// $Id$

class TClonesArray;

#include "AliAnalysisTaskSE.h"

class AliEmcalMCTrackSelector : public AliAnalysisTaskSE {
 public:
  AliEmcalMCTrackSelector();
  AliEmcalMCTrackSelector(const char *name);
  virtual ~AliEmcalMCTrackSelector();

  void UserCreateOutputObjects();
  void UserExec(Option_t *option);

  void SetTracksOutName(const char *name)            { fTracksOutName    = name ; }
  void SetRejectNK(Bool_t r = kTRUE)                 { fRejectNK         = r    ; }
  void SetChargedMC(Bool_t c = kTRUE)                { fChargedMC        = c    ; }

 protected:
  TString            fTracksOutName;        // name of output track array
  Bool_t             fRejectNK;             // true = reject k0l and neutrons
  Bool_t             fChargedMC;            // true = only charged particles
  TClonesArray      *fTracksOut;            //!track array out

 private:
  AliEmcalMCTrackSelector(const AliEmcalMCTrackSelector&);            // not implemented
  AliEmcalMCTrackSelector &operator=(const AliEmcalMCTrackSelector&); // not implemented

  ClassDef(AliEmcalMCTrackSelector, 0); // Task to select tracks in MC events
};
#endif
