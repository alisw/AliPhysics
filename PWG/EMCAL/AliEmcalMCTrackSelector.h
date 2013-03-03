#ifndef ALIEMCALMCTRACKSELECTOR_H
#define ALIEMCALMCTRACKSELECTOR_H

// $Id$

class TClonesArray;
class TString;
class TH1I;
class AliNamedArrayI;

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
  Int_t              GetNumberOfTracks() const;
  AliVParticle      *GetTrack(Int_t i);
  void               AddTrack(AliVParticle *track, Int_t nacc);

  TString            fTracksOutName;        // name of output track array
  Bool_t             fRejectNK;             // true = reject k0l and neutrons
  Bool_t             fChargedMC;            // true = only charged particles
  Bool_t             fInit;                 // true = task initialized
  TString            fTracksMapName;        // name of the track map
  Bool_t             fEsdMode;              //!switch for ESD/AOD mode
  TClonesArray      *fTracksIn;             //!track array in (AOD only)
  TClonesArray      *fTracksOut;            //!track array out
  AliNamedArrayI    *fTracksMap;            //!track mapping

 private:
  AliEmcalMCTrackSelector(const AliEmcalMCTrackSelector&);            // not implemented
  AliEmcalMCTrackSelector &operator=(const AliEmcalMCTrackSelector&); // not implemented

  ClassDef(AliEmcalMCTrackSelector, 2); // Task to select tracks in MC events
};
#endif
