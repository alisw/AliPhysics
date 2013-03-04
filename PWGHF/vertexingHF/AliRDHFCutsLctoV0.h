#ifndef ALIRDHFCUTSLCTOV0_H
#define ALIRDHFCUTSLCTOV0_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//***********************************************************
// Class AliRDHFCutsLctoV0
// class for cuts on AOD reconstructed Lc-> V0 + bachelor
//***********************************************************

#include "AliAODPidHF.h"
#include "AliRDHFCuts.h"

class AliRDHFCutsLctoV0 : public AliRDHFCuts
{
 public:

  enum ELctoV0channel {
    kLcToK0Spr=0x0001,
    kLcToLBarpi=0x0002,
    kLcToLpi=0x0004
  };

 enum ELctoV0pidStrategy {
  kTOFandTPC=0,
  kTOForTPCveto=1,
  kTOFandTPCasym1=2,
  kTOFandTPCasym2=3
 };

  AliRDHFCutsLctoV0(const char* name="CutsLctoV0", Short_t v0channel=0);

  virtual ~AliRDHFCutsLctoV0();

  AliRDHFCutsLctoV0(const AliRDHFCutsLctoV0& source);
  AliRDHFCutsLctoV0& operator=(const AliRDHFCutsLctoV0& source);

  using AliRDHFCuts::GetCutVarsForOpt;
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters);

  using AliRDHFCuts::IsSelected;
  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel);

  using AliRDHFCuts::IsSelectedPID;
  virtual Int_t IsSelectedPID(AliAODRecoDecayHF* obj);

  Int_t IsSelectedSingleCut(TObject* obj, Int_t selectionLevel, Int_t cutIndex);

  Int_t CombinePIDCuts (Int_t returnvalue, Int_t returnvaluePID) const;

  Float_t GetMassCut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(0,iPtBin)] : 1.e6);}
  Float_t GetDCACut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(7,iPtBin)] : 1.e6);}

  void SetPidSelectionFlag(Int_t a) {fPidSelectionFlag=a;}
  Int_t GetPidSelectionFlag() {return fPidSelectionFlag;}

  Int_t GetV0Type();

  virtual void SetStandardCutsPP2010();
  virtual void SetStandardCutsPbPb2010();
  virtual void SetStandardCutsPbPb2011();

  virtual Bool_t IsInFiducialAcceptance(Double_t pt,Double_t y) const;

  void SetPidV0pos(AliAODPidHF* pidV0pos) {
    if (fPidHFV0pos) delete fPidHFV0pos;
    fPidHFV0pos = new AliAODPidHF(*pidV0pos);
  }
  void SetPidV0neg(AliAODPidHF* pidV0neg) {
    if (fPidHFV0neg) delete fPidHFV0neg;
    fPidHFV0neg = new AliAODPidHF(*pidV0neg);
  }

  AliAODPidHF * GetPidV0pos() { return fPidHFV0pos; }
  AliAODPidHF * GetPidV0neg() { return fPidHFV0neg; }

  void AddTrackCutsV0daughters(AliESDtrackCuts* v0daug) {
    fV0daughtersCuts = new AliESDtrackCuts(*v0daug);
  }
  virtual AliESDtrackCuts *GetTrackCutsV0daughters() const {return fV0daughtersCuts;}

  virtual void PrintAll() const;
 protected:

  void CheckPID(AliAODTrack *bachelor, AliAODTrack *v0Neg, AliAODTrack *v0Pos,
		Bool_t &isBachelorID1, Bool_t &isV0NegID2, Bool_t &isV0PosID4);

 private:

  Int_t fPidSelectionFlag;
  AliAODPidHF *fPidHFV0pos;
  AliAODPidHF *fPidHFV0neg;
  AliESDtrackCuts *fV0daughtersCuts; // cuts for v0 daughters (AOD converted to ESD on the flight!)
  Float_t     fV0Type; // V0 type -- should be defined as in AliRDHFCuts.h

  //UShort_t fV0channel;

  ClassDef(AliRDHFCutsLctoV0,4);  // class for cuts on AOD reconstructed Lc->V0+bachelor
};

#endif
