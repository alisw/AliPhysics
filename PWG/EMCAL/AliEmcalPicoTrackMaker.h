#ifndef ALIEMCALPICOTRACKMAKER_H
#define ALIEMCALPICOTRACKMAKER_H

// $Id: AliEmcalPicoTrackMaker.h | Fri Dec 6 10:29:04 2013 +0100 | Constantin Loizides  $

class TClonesArray;
class AliVTrack;

#include "AliAnalysisTaskSE.h"

class AliEmcalPicoTrackMaker : public AliAnalysisTaskSE {
 public:
  AliEmcalPicoTrackMaker();
  AliEmcalPicoTrackMaker(const char *name);
  virtual ~AliEmcalPicoTrackMaker();

  void               SetTrackEfficiency(Double_t eff = 0.95)           { fTrackEfficiency   = eff ; }
  void               SetTrackEtaLimits(Double_t min, Double_t max)     { fMaxTrackEta       = max ; fMinTrackEta      = min ; }
  void               SetTrackPhiLimits(Double_t min, Double_t max)     { fMaxTrackPhi       = max ; fMinTrackPhi      = min ; }
  void               SetTrackPtLimits(Double_t min, Double_t max)      { fMaxTrackPt        = max ; fMinTrackPt       = min ; }
  void               SetTracksInName(const char *name)                 { fTracksInName      = name; }
  void               SetTracksOutName(const char *name)                { fTracksOutName     = name; }

 protected:
  void               UserCreateOutputObjects();
  void               UserExec(Option_t *option);

  Int_t              fAODfilterBits[2];     // AOD track filter bit map
  TString            fTracksOutName;        // name of output track array
  TString            fTracksInName;         // name of input track array
  Double_t           fMinTrackPt;           // mix pt of tracks
  Double_t           fMaxTrackPt;           // max pt of tracks
  Double_t           fMinTrackEta;          // cut on track eta
  Double_t           fMaxTrackEta;          // cut on track eta
  Double_t           fMinTrackPhi;          // cut on track phi
  Double_t           fMaxTrackPhi;          // cut on track phi
  Double_t           fTrackEfficiency;      // track efficiency
  TClonesArray      *fTracksIn;             //!track array in
  TClonesArray      *fTracksOut;            //!track array out

 private:
  AliEmcalPicoTrackMaker(const AliEmcalPicoTrackMaker&);            // not implemented
  AliEmcalPicoTrackMaker &operator=(const AliEmcalPicoTrackMaker&); // not implemented

  ClassDef(AliEmcalPicoTrackMaker, 7); // Task to make PicoTracks in AOD/ESD events
};
#endif
