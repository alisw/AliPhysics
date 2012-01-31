/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

// AliStarTrackCuts:
// A track cut class for AliStarTrack
// origin: Mikolaj Krzewicki (mikolaj.krzewicki@cern.ch)

#ifndef ALISTARTRACKCUTS_H
#define ALISTARTRACKCUTS_H

#include <float.h>
#include "TNamed.h"

class AliStarTrack;

class AliStarTrackCuts : public TNamed {

 public:
  AliStarTrackCuts();
  //AliStarTrackCuts(const AliStarTrackCuts& someCuts);
  //AliStarTrackCuts& operator=(const AliStarTrackCuts& someCuts);
  virtual  ~AliStarTrackCuts() {}
  
  Bool_t PassesCuts(const AliStarTrack* track) const;
  static AliStarTrackCuts* StandardCuts();

  void SetIDMax(const Int_t value) {fIDMax=value;fCutID=kTRUE;};
  void SetIDMin(const Int_t value) {fIDMin=value;fCutID=kTRUE;};
  void SetChargeMax(const Int_t value) {fChargeMax=value;fCutCharge=kTRUE;};
  void SetChargeMin(const Int_t value) {fChargeMin=value;fCutCharge=kTRUE;};
  void SetEtaMax(const Float_t value) {fEtaMax=value;fCutEta=kTRUE;};
  void SetEtaMin(const Float_t value) {fEtaMin=value;fCutEta=kTRUE;};
  void SetPhiMax(const Float_t value) {fPhiMax=value;fCutPhi=kTRUE;};
  void SetPhiMin(const Float_t value) {fPhiMin=value;fCutPhi=kTRUE;};
  void SetPtMax(const Float_t value) {fPtMax=value;fCutPt=kTRUE;};
  void SetPtMin(const Float_t value) {fPtMin=value;fCutPt=kTRUE;};
  void SetDCAMax(const Float_t value) {fDCAMax=value;fCutDCA=kTRUE;};
  void SetDCAMin(const Float_t value) {fDCAMin=value;fCutDCA=kTRUE;};
  void SetNHitsMax(const Int_t value) {fNHitsMax=value;fCutNHits=kTRUE;};
  void SetNHitsMin(const Int_t value) {fNHitsMin=value;fCutNHits=kTRUE;};
  void SetNHitsFitMax(const Int_t value) {fNHitsFitMax=value;fCutNHitsFit=kTRUE;};
  void SetNHitsFitMin(const Int_t value) {fNHitsFitMin=value;fCutNHitsFit=kTRUE;};
  void SetNHitsPossMax(const Int_t value) {fNHitsPossMax=value;fCutNHitsPoss=kTRUE;};
  void SetNHitsPossMin(const Int_t value) {fNHitsPossMin=value;fCutNHitsPoss=kTRUE;};
  void SetNHitsDedxMax(const Int_t value) {fNHitsDedxMax=value;fCutNHitsDedx=kTRUE;};
  void SetNHitsDedxMin(const Int_t value) {fNHitsDedxMin=value;fCutNHitsDedx=kTRUE;};
  void SetdEdxMax(const Float_t value) {fdEdxMax=value;fCutdEdx=kTRUE;};
  void SetdEdxMin(const Float_t value) {fdEdxMin=value;fCutdEdx=kTRUE;};
  void SetNSigElectMax(const Float_t value) {fNSigElectMax=value;fCutNSigElect=kTRUE;};
  void SetNSigElectMin(const Float_t value) {fNSigElectMin=value;fCutNSigElect=kTRUE;};
  void SetNSigPiMax(const Float_t value) {fNSigPiMax=value;fCutNSigPi=kTRUE;};
  void SetNSigPiMin(const Float_t value) {fNSigPiMin=value;fCutNSigPi=kTRUE;};
  void SetNSigKMax(const Float_t value) {fNSigKMax=value;fCutNSigK=kTRUE;};
  void SetNSigKMin(const Float_t value) {fNSigKMin=value;fCutNSigK=kTRUE;};
  void SetNSigProtonMax(const Float_t value) {fNSigProtonMax=value;fCutNSigProton=kTRUE;};
  void SetNSigProtonMin(const Float_t value) {fNSigProtonMin=value;fCutNSigProton=kTRUE;};
  void SetFitRatioMax(const Float_t value) {fFitRatioMax=value;fCutFitRatio=kTRUE;};
  void SetFitRatioMin(const Float_t value) {fFitRatioMin=value;fCutFitRatio=kTRUE;};

  Int_t GetIDMax() const {return fIDMax;}
  Int_t GetIDMin() const {return fIDMin;}
  Int_t GetChargeMax() const {return fChargeMax;}
  Int_t GetChargeMin() const {return fChargeMin;}
  Float_t GetEtaMax() const {return fEtaMax;}
  Float_t GetEtaMin() const {return fEtaMin;}
  Float_t GetPhiMax() const {return fPhiMax;}
  Float_t GetPhiMin() const {return fPhiMin;}
  Float_t GetPtMax() const {return fPtMax;}
  Float_t GetPtMin() const {return fPtMin;}
  Float_t GetDCAMax() const {return fDCAMax;}
  Float_t GetDCAMin() const {return fDCAMin;}
  Int_t GetNHitsMax() const {return fNHitsMax;}
  Int_t GetNHitsMin() const {return fNHitsMin;}
  Int_t GetNHitsFitMax() const {return fNHitsFitMax;}
  Int_t GetNHitsFitMin() const {return fNHitsFitMin;}
  Int_t GetNHitsPossMax() const {return fNHitsPossMax;}
  Int_t GetNHitsPossMin() const {return fNHitsPossMin;}
  Int_t GetNHitsDedxMax() const {return fNHitsDedxMax;}
  Int_t GetNHitsDedxMin() const {return fNHitsDedxMin;}
  Float_t GetdEdxMax() const {return fdEdxMax;}
  Float_t GetdEdxMin() const {return fdEdxMin;}
  Float_t GetNSigElectMax() const {return fNSigElectMax;}
  Float_t GetNSigElectMin() const {return fNSigElectMin;}
  Float_t GetNSigPiMax() const {return fNSigPiMax;}
  Float_t GetNSigPiMin() const {return fNSigPiMin;}
  Float_t GetNSigKMax() const {return fNSigKMax;}
  Float_t GetNSigKMin() const {return fNSigKMin;}
  Float_t GetNSigProtonMax() const {return fNSigProtonMax;}
  Float_t GetNSigProtonMin() const {return fNSigProtonMin;}
  Float_t GetFitRatioMax() const {return fFitRatioMax;};
  Float_t GetFitRatioMin() const {return fFitRatioMin;};

 private:
  Bool_t   fCutID;  //cut on id
  Int_t fIDMax;     //id max
  Int_t fIDMin;     //id min
  Bool_t   fCutCharge;  //cut on charge
  Int_t fChargeMax;     //max charge
  Int_t fChargeMin;     //min charge
  Bool_t   fCutEta;     //cut on eta
  Float_t fEtaMax;      //eta max
  Float_t fEtaMin;      //eta min
  Bool_t   fCutPhi;     //cut on phi
  Float_t fPhiMax;      //phi max
  Float_t fPhiMin;      //phi min
  Bool_t   fCutPt;      //cut on pt
  Float_t fPtMax;       //pt max
  Float_t fPtMin;       //pt min
  Bool_t   fCutDCA;     //cut dca
  Float_t fDCAMax;      //max dca
  Float_t fDCAMin;      //min dca
  Bool_t   fCutNHits;   //cut nhits
  Int_t fNHitsMax;      //max nhits
  Int_t fNHitsMin;      //min nhits
  Bool_t   fCutNHitsFit;//cut nhistfit
  Int_t fNHitsFitMax;   //nhistfit max
  Int_t fNHitsFitMin;   //nhistfit min
  Bool_t   fCutNHitsPoss;  //cut nhitsposs
  Int_t fNHitsPossMax;     //nhitsposs max
  Int_t fNHitsPossMin;     //nhitsposs min
  Bool_t   fCutNHitsDedx;  //cut nhitsdedx
  Int_t fNHitsDedxMax;     //nhitsdedx max
  Int_t fNHitsDedxMin;     //nhits min
  Bool_t   fCutdEdx;       //cut dedx
  Float_t fdEdxMax;        //dedx max
  Float_t fdEdxMin;        //dedx min
  Bool_t   fCutNSigElect;  //cut nsigelect
  Float_t fNSigElectMax;   //nsigelect max
  Float_t fNSigElectMin;   //nsigelect min
  Bool_t   fCutNSigPi;     //cut nsigpi
  Float_t fNSigPiMax;      //max nsigpi
  Float_t fNSigPiMin;      //min nsigpi
  Bool_t   fCutNSigK;      //cut nsigk
  Float_t fNSigKMax;       //max nsigk
  Float_t fNSigKMin;       //min nsigk
  Bool_t   fCutNSigProton; //cut nsigproton
  Float_t fNSigProtonMax;  //max nsigproton
  Float_t fNSigProtonMin;  //min nsigproton
  Bool_t  fCutFitRatio;    //cut fitratio
  Float_t fFitRatioMax;    //max fitratio
  Float_t fFitRatioMin;    //min fitratio

  ClassDef(AliStarTrackCuts,1)
};

#endif


