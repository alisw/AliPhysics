#ifndef ALIEMCALPICOTRACKMAKER_H
#define ALIEMCALPICOTRACKMAKER_H

// $Id$

class TClonesArray;
class AliVEvent;
class AliVTrack;
class AliESDtrackCuts;

#include "AliAnalysisTaskSE.h"

class AliEmcalPicoTrackMaker : public AliAnalysisTaskSE {
 public:
  AliEmcalPicoTrackMaker();
  AliEmcalPicoTrackMaker(const char *name);
  virtual ~AliEmcalPicoTrackMaker();

  void UserCreateOutputObjects();
  void UserExec(Option_t *option);

  void SetAODfilterBits(Int_t b0 = 0, Int_t b1 = 0)  { fAODfilterBits[0] = b0; fAODfilterBits[1] = b1; }
  void SetESDtrackCuts(AliESDtrackCuts *cuts)        { fESDtrackCuts     = cuts; }
  void SetTracksInName(const char *name)             { fTracksInName     = name; }
  void SetTracksOutName(const char *name)            { fTracksOutName    = name; }
  void SetMaxTrackPt(Float_t pt)                     { fMaxTrackPt       = pt  ; }

 protected:
  Int_t              fAODfilterBits[2];     // AOD track filter bit map
  AliESDtrackCuts   *fESDtrackCuts;         // ESD track cuts
  TString            fTracksOutName;        // name of output track array
  TString            fTracksInName;         // name of input track array
  Float_t            fMaxTrackPt;           // max pt of tracks
  TClonesArray      *fTracksIn;             //!track array in
  TClonesArray      *fTracksOut;            //!track array out

 private:
  AliEmcalPicoTrackMaker(const AliEmcalPicoTrackMaker&);            // not implemented
  AliEmcalPicoTrackMaker &operator=(const AliEmcalPicoTrackMaker&); // not implemented

  ClassDef(AliEmcalPicoTrackMaker, 1); // Task to make PicoTracks in AOD/ESD events
};
#endif
