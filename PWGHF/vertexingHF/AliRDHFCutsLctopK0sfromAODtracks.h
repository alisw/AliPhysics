#ifndef ALIRDHFCUTSLCTOPK0SFROMAODTRACKS_H
#define ALIRDHFCUTSLCTOPK0SFROMAODTRACKS_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//***********************************************************
/// \class Class AliRDHFCutsXictoPLUSXiPiPifromAODtracks
/// \brief class for cuts on AOD reconstructed Lc-> pK0s
//***********************************************************

#include "AliRDHFCuts.h"

class AliRDHFCutsLctopK0sfromAODtracks : public AliRDHFCuts
{
 public:

  enum EPIDStrategy{
    kNSigmaCuts,
    kCombinedCuts
  };

  AliRDHFCutsLctopK0sfromAODtracks(const char* name="CutsLctopK0s");
  virtual ~AliRDHFCutsLctopK0sfromAODtracks();
  AliRDHFCutsLctopK0sfromAODtracks(const AliRDHFCutsLctopK0sfromAODtracks& source);
  AliRDHFCutsLctopK0sfromAODtracks& operator=(const AliRDHFCutsLctopK0sfromAODtracks& source);

  using AliRDHFCuts::GetCutVarsForOpt;
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters);

  using AliRDHFCuts::IsSelected;
  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel);
  using AliRDHFCuts::IsSelectedPID;
  virtual Int_t IsSelectedPID(AliAODRecoDecayHF* obj);
  Int_t IsSelectedCombinedPID(AliAODRecoDecayHF* obj);
  Double_t GetProtonProbabilityTPCTOF(AliAODTrack *trk);

  void SetPIDStrategy(EPIDStrategy pidStrategy){fPIDStrategy=pidStrategy;}
  EPIDStrategy GetPIDStrategy() const {return fPIDStrategy;}
  void SetCombinedPIDThreshold(Double_t a){fCombinedPIDThreshold=a;}
  Double_t GetCombinedPIDThreshold(){return fCombinedPIDThreshold;}

  void SetUseOnTheFlyV0(Bool_t a) { fUseOnTheFlyV0=a; }
  Bool_t GetUseOnTheFlyV0() { return fUseOnTheFlyV0; }

  Bool_t SingleTrkCuts(AliAODTrack *trk);
  Bool_t SingleV0Cuts(AliAODv0 *v0, AliAODVertex *vert);
  Bool_t SelectWithRoughCuts(AliAODv0 *v0, AliAODTrack *trk1);

  void SetProdTrackPtMin(Double_t a){fProdTrackPtMin=a;}
  void SetProdTrackEtaRange(Double_t a){fProdTrackEtaRange=a;}
  void SetProdUseAODFilterBit(Bool_t a){fProdUseAODFilterBit=a;}
  void SetProdV0MassTolK0s(Double_t a){fProdV0MassTolK0s=a;}
  void SetProdV0PtMin(Double_t a){fProdV0PtMin=a;}
  void SetProdV0CosPointingAngleToPrimVtxMin(Double_t a){fProdV0CosPointingAngleToPrimVtxMin=a;}
  void SetProdV0DcaDaughtersMax(Double_t a){fProdV0DcaDaughtersMax=a;}
  void SetProdV0DaughterEtaRange(Double_t a){fProdV0DaughterEtaRange=a;}
  void SetProdV0DaughterPtMin(Double_t a){fProdV0DaughterPtMin=a;}
  void SetProdV0DaughterTPCClusterMin(Double_t a){fProdV0DaughterTPCClusterMin=a;}
  
  void SetProdRoughMassTol(Double_t a){fProdRoughMassTol=a;}
  void SetProdRoughPtMin(Double_t a){fProdRoughPtMin=a;}
  
  Double_t GetProdTrackPtMin(){return fProdTrackPtMin;}
  Double_t GetProdTrackEtaRange(){return fProdTrackEtaRange;}
  Bool_t   GetProdUseAODFilterBit(){return fProdUseAODFilterBit;}
  Double_t GetProdV0MassTolK0s(){return fProdV0MassTolK0s;}
  Double_t GetProdV0PtMin(){return fProdV0PtMin;}
  Double_t GetProdV0CosPointingAngleToPrimVtxMin(){return fProdV0CosPointingAngleToPrimVtxMin;}
  Double_t GetProdV0DcaDaughtersMax(){return fProdV0DcaDaughtersMax;}
  Double_t GetProdV0DaughterEtaRange(){return fProdV0DaughterEtaRange;}
  Double_t GetProdV0DaughterPtMin(){return fProdV0DaughterPtMin;}
  Double_t GetProdV0DaughterTPCClusterMin(){return fProdV0DaughterTPCClusterMin;}
  
  Double_t GetProdRoughMassTol(){return fProdRoughMassTol;}
  Double_t GetProdRoughPtMin(){return fProdRoughPtMin;}
  
  Double_t CalculateLcCosPAXY(AliAODRecoDecayHF *obj);

 protected:
	
 private:

  EPIDStrategy fPIDStrategy;        /// PID strategy
  Double_t fCombinedPIDThreshold;   /// Threshold used in  IsSelectedCombinedPID
  Bool_t   fUseOnTheFlyV0;          /// Flag to check if we use on-the-fly v0
  
  Double_t fProdTrackPtMin;         /// Minimum Track pT used before object creation
  Double_t fProdTrackEtaRange;      /// eta range used before object creation
  Bool_t   fProdUseAODFilterBit;    /// Flag for AOD filter Bit used before object creation
  Double_t fProdV0MassTolK0s;       /// K0s mass selection  used before object creation
  Double_t fProdV0PtMin;            /// Minimum K0s pT used before object creation
  Double_t fProdV0CosPointingAngleToPrimVtxMin; /// V0 pointing angle used before object creation
  Double_t fProdV0DcaDaughtersMax;  /// Max DCA between V0 daughters used before object creation
  Double_t fProdV0DaughterEtaRange; /// V0Daughter eta range used before object creation
  Double_t fProdV0DaughterPtMin;    /// V0 Daughter pT min used before object creation
  Double_t fProdV0DaughterTPCClusterMin; /// V0 daughter Minimum TPC cluster pT used before object creation
  Double_t fProdRoughMassTol;       /// Mass cut for Lc used before object creation
  Double_t fProdRoughPtMin;         /// pT cut for Lc used before object creation
  
  /// \cond CLASSIMP
  ClassDef(AliRDHFCutsLctopK0sfromAODtracks,2);
  /// \endcond
};

#endif
