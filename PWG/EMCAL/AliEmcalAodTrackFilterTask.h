#ifndef ALIEMCALAODTRACKFILTERTASK_H
#define ALIEMCALAODTRACKFILTERTASK_H

// $Id$

class TClonesArray;
class AliVEvent;
class AliVTrack;

#include "AliAnalysisTaskSE.h"

class AliEmcalAodTrackFilterTask : public AliAnalysisTaskSE {
 public:
  AliEmcalAodTrackFilterTask();
  AliEmcalAodTrackFilterTask(const char *name);
  virtual ~AliEmcalAodTrackFilterTask();

  void UserCreateOutputObjects();
  void UserExec(Option_t *option);

  void SetAODfilterBits(Int_t b0 = 0, Int_t b1 = 0)         { fAODfilterBits[0]  = b0  ; fAODfilterBits[1] = b1  ; }
  void SetCutMaxFractionSharedTPCClusters(Double_t c = 0.4) { fCutMaxFrShTPCClus = c   ; }
  void SetDist(Double_t d)                                  { fDist              = d   ; }
  void SetDoPropagation(Bool_t b)                           { fDoPropagation     = b   ; }
  void SetIncludeNoITS(Bool_t f)                            { fIncludeNoITS      = f   ; }
  void SetMC(Bool_t a)                                      { fIsMC              = a   ; }
  void SetModifyTrack(Bool_t b)                             { fModifyTrack       = b;    }
  void SetTrackEfficiency(Double_t eff = 0.95)              { fTrackEfficiency   = eff ; }
  void SetTrackEtaLimits(Double_t min, Double_t max)        { fMaxTrackEta       = max ; fMinTrackEta      = min ; }
  void SetTrackPhiLimits(Double_t min, Double_t max)        { fMaxTrackPhi       = max ; fMinTrackPhi      = min ; }
  void SetTrackPtLimits(Double_t min, Double_t max)         { fMaxTrackPt        = max ; fMinTrackPt       = min ; }
  void SetTracksInName(const char *name)                    { fTracksInName      = name; }
  void SetTracksOutName(const char *name)                   { fTracksOutName     = name; }
  void SetUseNegativeLabels(Bool_t f)                       { fUseNegativeLabels = f   ; }

 protected:
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
  Bool_t             fIncludeNoITS;         // includes tracks with failed ITS refit
  Bool_t             fUseNegativeLabels;    // whether or not should use negative MC labels
  Bool_t             fIsMC;                 // whether it is a MC event or not
  Double_t           fCutMaxFrShTPCClus;    // max fraction of shared TPC clusters
  Bool_t             fModifyTrack;          // if true then overwrite some fields in AodTrack
  Bool_t             fDoPropagation;        // if true then propagate all hybrid tracks to EMCal surface
  Double_t           fDist;                 // distance to surface (440cm default)
  TClonesArray      *fTracksIn;             //!track array in
  TClonesArray      *fTracksOut;            //!track array out

 private:
  AliEmcalAodTrackFilterTask(const AliEmcalAodTrackFilterTask&);            // not implemented
  AliEmcalAodTrackFilterTask &operator=(const AliEmcalAodTrackFilterTask&); // not implemented

  ClassDef(AliEmcalAodTrackFilterTask, 1); // Task to filter Aod tracks
};
#endif
