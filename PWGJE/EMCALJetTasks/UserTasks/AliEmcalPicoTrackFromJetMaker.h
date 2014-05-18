#ifndef ALIEMCALPICOTRACKFROMJETMAKER_H
#define ALIEMCALPICOTRACKFROMJETMAKER_H

class TClonesArray;
class AliVTrack;

#include "AliAnalysisTaskSE.h"

class AliEmcalPicoTrackFromJetMaker : public AliAnalysisTaskSE {
 public:
  AliEmcalPicoTrackFromJetMaker();
  AliEmcalPicoTrackFromJetMaker(const char *name);
  virtual ~AliEmcalPicoTrackFromJetMaker();

  void               SetJetsInName(const char *name)                   { fJetsInName        = name; }
  void               SetTracksOutName(const char *name)                { fTracksOutName     = name; }

 protected:
  void               UserCreateOutputObjects();
  void               UserExec(Option_t *option);

  TString            fTracksOutName;        // name of output track array
  TString            fJetsInName;           // name of input jet array
  TClonesArray      *fJetsIn;               //!jet array in
  TClonesArray      *fTracksOut;            //!track array out

 private:
  AliEmcalPicoTrackFromJetMaker(const AliEmcalPicoTrackFromJetMaker&);            // not implemented
  AliEmcalPicoTrackFromJetMaker &operator=(const AliEmcalPicoTrackFromJetMaker&); // not implemented

  ClassDef(AliEmcalPicoTrackFromJetMaker, 1); // Task to make PicoTracks from jet 4-vectors
};
#endif
