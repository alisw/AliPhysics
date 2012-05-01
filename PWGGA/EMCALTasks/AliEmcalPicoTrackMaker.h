#ifndef ALIEMCALPICOTRACKMAKER_H
#define ALIEMCALPICOTRACKMAKER_H

// $Id: AliEmcalPicoTrackMaker.h 54003 2012-01-19 16:40:42Z loizides $

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

 protected:
  Int_t              fAODfilterBits[2];     // AOD track filter bit map
  AliESDtrackCuts   *fESDtrackCuts;         // ESD track cuts
  TString            fTracksOutName;        // name of output track array
  TString            fTracksInName;         // name of input track array
  TClonesArray      *fTracksIn;             //!track array in
  TClonesArray      *fTracksOut;            //!track array out

 private:
  AliEmcalPicoTrackMaker(const AliEmcalPicoTrackMaker&);            // not implemented
  AliEmcalPicoTrackMaker &operator=(const AliEmcalPicoTrackMaker&); // not implemented

  ClassDef(AliEmcalPicoTrackMaker, 1); // Task to make PicoTracks in AOD/ESD events
};
#endif
