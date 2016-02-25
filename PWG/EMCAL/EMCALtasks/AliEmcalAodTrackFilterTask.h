#ifndef ALIEMCALAODTRACKFILTERTASK_H
#define ALIEMCALAODTRACKFILTERTASK_H

class TClonesArray;

#include <TF1.h>

#include "AliAnalysisTaskSE.h"

class AliEmcalAodTrackFilterTask : public AliAnalysisTaskSE {
 public:
  AliEmcalAodTrackFilterTask();
  AliEmcalAodTrackFilterTask(const char *name);
  virtual ~AliEmcalAodTrackFilterTask();

  void               SetAODfilterBits(Int_t b0 = 0, Int_t b1 = 0)         { fAODfilterBits[0]  = b0  ; fAODfilterBits[1] = b1  ; }
  void               SetAttemptProp(Bool_t b)                             { fAttemptProp       = b   ; }
  void               SetAttemptPropMatch(Bool_t b)                        { fAttemptPropMatch  = b   ; }
  void               SetCutMaxFractionSharedTPCClusters(Double_t c = 0.4) { fCutMaxFrShTPCClus = c   ; }
  void               SetDist(Double_t d)                                  { fDist              = d   ; }
  void               SetDoPropagation(Bool_t b)                           { fDoPropagation     = b   ; }
  void               SetIncludeNoITS(Bool_t f)                            { fIncludeNoITS      = f   ; }
  void               SetMC(Bool_t b)                                      { fIsMC              = b   ; }
  void               SetTracksInName(const char *name)                    { fTracksInName      = name; }
  void               SetTracksOutName(const char *name)                   { fTracksOutName     = name; }
  void               SetUseNegativeLabels(Bool_t f)                       { fUseNegativeLabels = f   ; }
  void               SetTrackEfficiency(Double_t eff = 0.95)              { fTrackEfficiency  = new TF1("eff", "[0]", 0, 500); fTrackEfficiency->FixParameter(0,eff); }
  void               SetKeepInvMassTag(Bool_t f)                          { fKeepInvMassTag = f ; }
  void               SetTrackEfficiency(TF1* eff)                         { fTrackEfficiency  = eff  ; }

 protected:
  void               UserCreateOutputObjects();
  void               UserExec(Option_t *option);

  Int_t              fAODfilterBits[2];     // AOD track filter bit map
  TString            fTracksOutName;        // name of output track array
  TString            fTracksInName;         // name of input track array
  Bool_t             fIncludeNoITS;         // includes tracks with failed ITS refit
  Double_t           fCutMaxFrShTPCClus;    // max fraction of shared TPC clusters
  Bool_t             fUseNegativeLabels;    // whether or not should use negative MC labels
  Bool_t             fIsMC;                 // whether it is a MC event or not
  Bool_t             fDoPropagation;        // if true then propagate all hybrid tracks to EMCal surface
  Bool_t             fAttemptProp;          // if true then attempt to propagate if not done yet
  Bool_t             fAttemptPropMatch;     // if true then attempt to propagate if not done yet but IsEMCAL is true
  Bool_t             fKeepInvMassTag;     // if true then pass in track container labels for tagging tracks in jets
  Double_t           fDist;                 // distance to surface (440cm default)
  TF1               *fTrackEfficiency;      // track efficiency
  TClonesArray      *fTracksIn;             //!track array in
  TClonesArray      *fTracksOut;            //!track array out

 private:
  AliEmcalAodTrackFilterTask(const AliEmcalAodTrackFilterTask&);            // not implemented
  AliEmcalAodTrackFilterTask &operator=(const AliEmcalAodTrackFilterTask&); // not implemented

  ClassDef(AliEmcalAodTrackFilterTask, 4); // Task to filter Aod tracks
};
#endif
