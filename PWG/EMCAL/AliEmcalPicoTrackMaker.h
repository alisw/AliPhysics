#ifndef ALIEMCALPICOTRACKMAKER_H
#define ALIEMCALPICOTRACKMAKER_H

class TClonesArray;
class AliVParticle;
class AliNamedArrayI;

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
  void               SetMCParticlesName(const char *name)              { fMCParticlesName   = name; }
  void               SetCopyMCFlag(Bool_t c, const char* name)         { fCopyMCFlag        = c   ; fMCParticlesName  = name; }
  

 protected:
  void               UserCreateOutputObjects();
  void               UserExec(Option_t *option);

  AliVParticle*      GetMCParticle(Int_t label);

  Int_t              fAODfilterBits[2];     // AOD track filter bit map
  TString            fTracksOutName;        // name of output track array
  TString            fTracksInName;         // name of input track array
  TString            fMCParticlesName;      // name of MC particle array, used by IsHIJINGParticle
  Double_t           fMinTrackPt;           // mix pt of tracks
  Double_t           fMaxTrackPt;           // max pt of tracks
  Double_t           fMinTrackEta;          // cut on track eta
  Double_t           fMaxTrackEta;          // cut on track eta
  Double_t           fMinTrackPhi;          // cut on track phi
  Double_t           fMaxTrackPhi;          // cut on track phi
  Double_t           fTrackEfficiency;      // track efficiency
  Bool_t             fCopyMCFlag;           // copy MC flag
  TClonesArray      *fTracksIn;             //!track array in
  TClonesArray      *fTracksOut;            //!track array out
  TClonesArray      *fMCParticles;          //!MC particle array
  AliNamedArrayI    *fMCParticlesMap;       //!MC particle map
  Bool_t             fInit;                 //!true = task initialized

 private:
  AliEmcalPicoTrackMaker(const AliEmcalPicoTrackMaker&);            // not implemented
  AliEmcalPicoTrackMaker &operator=(const AliEmcalPicoTrackMaker&); // not implemented

  ClassDef(AliEmcalPicoTrackMaker, 8); // Task to make PicoTracks in AOD/ESD events
};
#endif
