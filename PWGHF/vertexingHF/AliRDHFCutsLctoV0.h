#ifndef ALIRDHFCUTSLCTOV0_H
#define ALIRDHFCUTSLCTOV0_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//***********************************************************
/// \class Class AliRDHFCutsLctoV0
/// \brief class for cuts on AOD reconstructed Lc-> V0 + bachelor
//***********************************************************

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
  kTOFandTPCasym2=3,
  kTPClowTOFhigh=4,
  kTPClowTOFintermediateTOForTPChigh=5,
  kTPClowTOFintermediateTOForTPChighA=6,
  kTPClowTOForTPCintermediate_high=7,
  kTPClowTOForTPCintermediate_highA=8,
  kCombinedTPCTOF=9
 };

  AliRDHFCutsLctoV0(const char* name="CutsLctoV0", Short_t v0channel=0);

  virtual ~AliRDHFCutsLctoV0();

  AliRDHFCutsLctoV0(const AliRDHFCutsLctoV0& source);
  AliRDHFCutsLctoV0& operator=(const AliRDHFCutsLctoV0& source);

  using AliRDHFCuts::GetCutVarsForOpt;
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters);

  using AliRDHFCuts::IsSelected;
  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel) 
                         {return IsSelected(obj,selectionLevel,0);}
  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel,AliAODEvent* aod);

  using AliRDHFCuts::IsSelectedPID;
  virtual Int_t IsSelectedPID(AliAODRecoDecayHF* obj);


  Bool_t PreSelect(TObject* obj, AliAODv0 *v0, AliVTrack *bachelorTrack);

  Int_t IsSelectedSingleCut(TObject* obj, Int_t selectionLevel, Int_t cutIndex, AliAODEvent* aod=0x0);

  Int_t CombineCuts (Int_t returnvalueTrack, Int_t returnvalue, Int_t returnvaluePID) const;

  Float_t GetMassCut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(0,iPtBin)] : 1.e6);}
  Float_t GetDCACut(Int_t iPtBin=0) const { return (GetCuts() ? fCutsRD[GetGlobalIndex(7,iPtBin)] : 1.e6);}

  void SetPidSelectionFlag(Int_t a) {fPidSelectionFlag=a;}
  Int_t GetPidSelectionFlag() {return fPidSelectionFlag;}

  Bool_t AreLctoV0DaughtersSelected(AliAODRecoDecayHF *rd, AliAODEvent* aod=0x0) const;

  Int_t GetV0Type();

  void SetHighPtCut(Float_t highPtCut) {fHighPtCut=highPtCut;};
  Float_t GetHighPtCut() const {return fHighPtCut;};

  void SetLowPtCut(Float_t lowPtCut) {fLowPtCut=lowPtCut;};
  Float_t GetLowPtCut() const {return fLowPtCut;};

  void SetMinCombinedProbability(Int_t ptBins, Float_t *ptMin);
  const Float_t *GetMinCombinedProbability() {return fMinCombinedProbability;}

  void SetExcludedCut(Int_t excludedCut) {fExcludedCut=excludedCut;}
  Int_t GetExcludedCut(){return fExcludedCut;}

  virtual void SetStandardCutsPP2010();
  virtual void SetStandardCutsPbPb2010();
  virtual void SetStandardCutsPbPb2011();

  virtual Bool_t IsInFiducialAcceptance(Double_t pt,Double_t y) const;

  void AddTrackCutsV0daughters(AliESDtrackCuts* v0daug)
  { delete fV0daughtersCuts; fV0daughtersCuts = new AliESDtrackCuts(*v0daug); }
  virtual AliESDtrackCuts *GetTrackCutsV0daughters() const {return fV0daughtersCuts;}

  Double_t GetProtonEmissionAngleCMS(AliAODRecoDecayHF *d);
  Double_t GetReSignedd0(AliAODRecoDecayHF *d);
  void SetMagneticField(Double_t a){fBzkG = a;}

  virtual void PrintAll() const;
 protected:

  void CheckPID(Int_t candPtBin, AliAODTrack *bachelor, AliAODTrack * /*v0Neg*/, AliAODTrack * /*v0Pos*/,
		Bool_t &isBachelorID1, Bool_t &isBachelorID2, Bool_t &isBachelorID4);

 private:

  Int_t fPidSelectionFlag;
  AliESDtrackCuts *fV0daughtersCuts; /// cuts for v0 daughters (AOD converted to ESD on the flight!)
  Float_t     fV0Type; /// V0 type -- should be defined as in AliRDHFCuts.h
  Float_t fHighPtCut;  /// high pT cut separation for proton identification
  Float_t fLowPtCut;   /// low pT cut separation for proton identification
  Int_t   fExcludedCut; /// cut to be excluded (-1=none)
  Float_t *fMinCombinedProbability; //[fnPtBins] min value for combined PID probabilities
  Double_t fBzkG; //Magnetic field for propagation

  //UShort_t fV0channel;

  /// \cond CLASSIMP    
  ClassDef(AliRDHFCutsLctoV0,8);  /// class for cuts on AOD reconstructed Lc->V0+bachelor
  /// \endcond
};

#endif
