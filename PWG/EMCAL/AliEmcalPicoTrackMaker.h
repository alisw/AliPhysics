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

  void SetAODfilterBits(Int_t b0 = 0, Int_t b1 = 0)      { fAODfilterBits[0] = b0  ; fAODfilterBits[1] = b1  ; }
  void SetESDtrackCuts(AliESDtrackCuts *cuts)            { fESDtrackCuts     = cuts; }
  void SetTracksInName(const char *name)                 { fTracksInName     = name; }
  void SetTracksOutName(const char *name)                { fTracksOutName    = name; }
  void SetTrackPtLimits(Double_t min, Double_t max)      { fMaxTrackPt       = max ; fMinTrackPt       = min ; }
  void SetTrackEtaLimits(Double_t min, Double_t max)     { fMaxTrackEta      = max ; fMinTrackEta      = min ; }
  void SetTrackPhiLimits(Double_t min, Double_t max)     { fMaxTrackPhi      = max ; fMinTrackPhi      = min ; }
  void SetTrackEfficiency(Double_t eff = 0.95)           { fTrackEfficiency  = eff ; }
  void SetIncludeNoITS(Bool_t f)                         { fIncludeNoITS     = f   ; }
  void SetUseNegativeLabels(Bool_t f)                    { fUseNegativeLabels= f   ; }
  void SetMC(Bool_t a)                                   { fIsMC             = a   ; }

 protected:
  Int_t              fAODfilterBits[2];     // AOD track filter bit map
  AliESDtrackCuts   *fESDtrackCuts;         // ESD track cuts
  TString            fTracksOutName;        // name of output track array
  TString            fTracksInName;         // name of input track array
  Double_t           fMinTrackPt;           // mix pt of tracks
  Double_t           fMaxTrackPt;           // max pt of tracks
  Double_t           fMinTrackEta;          // cut on track eta
  Double_t           fMaxTrackEta;          // cut on track eta
  Double_t           fMinTrackPhi;          // cut on track phi
  Double_t           fMaxTrackPhi;          // cut on track phi
  Double_t           fTrackEfficiency;      // track efficiency
  Bool_t             fIncludeNoITS;         // includes tracks with failed ITS refit
  Bool_t             fUseNegativeLabels;    // whether or not should use negative MC labels
  Bool_t             fIsMC;                 // whether it is a MC event or not
  TClonesArray      *fTracksIn;             //!track array in
  TClonesArray      *fTracksOut;            //!track array out

 private:
  AliEmcalPicoTrackMaker(const AliEmcalPicoTrackMaker&);            // not implemented
  AliEmcalPicoTrackMaker &operator=(const AliEmcalPicoTrackMaker&); // not implemented

  ClassDef(AliEmcalPicoTrackMaker, 5); // Task to make PicoTracks in AOD/ESD events
};
#endif
