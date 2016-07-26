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

class TLorentzVector;

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
  Int_t IsSelected(TLorentzVector* t1, TLorentzVector *t2, Double_t *info, Int_t selectionLevel);
  Int_t IsSelectedCombinedPID(AliAODRecoDecayHF* obj);
  Double_t GetProtonProbabilityTPCTOF(AliAODTrack *trk);
  Bool_t IsSelectedProtonID(AliAODTrack* trk);
  Bool_t IsSelectedKaonID(AliAODTrack* trk);

  void SetPIDStrategy(EPIDStrategy pidStrategy){fPIDStrategy=pidStrategy;}
  EPIDStrategy GetPIDStrategy() const {return fPIDStrategy;}
  void SetCombinedPIDThreshold(Double_t a){fCombinedPIDThreshold=a;}
  Double_t GetCombinedPIDThreshold(){return fCombinedPIDThreshold;}

  void SetUseOnTheFlyV0(Bool_t a) { fUseOnTheFlyV0=a; }
  Bool_t GetUseOnTheFlyV0() { return fUseOnTheFlyV0; }

  Bool_t SingleTrkCuts(AliAODTrack *trk, AliAODTrack *trkpid, AliAODVertex *vtx);
  Bool_t SingleKaonCuts(AliAODTrack *trk, AliAODVertex *vtx);
  Bool_t SingleV0Cuts(AliAODv0 *v0, AliAODVertex *vert);
	Bool_t TagV0(AliAODTrack *etrk, AliAODEvent *evt, Int_t ntrk, Double_t &minmass);
	Bool_t TagV0SameSign(AliAODTrack *etrk, AliAODEvent *evt, Int_t ntrk, Double_t &minmass);
  Bool_t SelectWithRoughCuts(AliAODv0 *v0, AliAODTrack *trk1);
  Bool_t SelectWithRoughCuts(TLorentzVector *v0, TLorentzVector *trk1);
  Bool_t SelectWithRoughCutsWS(AliAODTrack *vka, AliAODTrack *trk1);
  Bool_t SelectWithRoughCutsWS(TLorentzVector *vka, TLorentzVector *trk1);

  void SetMagneticField(Double_t a){fBzkG = a;}
  void SetPrimaryVertex(Double_t *a){fPrimVert[0] = a[0];fPrimVert[1] = a[1];fPrimVert[2] = a[2];}

  void SetProdTrackTPCNclsPIDMin(Int_t a){fProdTrackTPCNclsPIDMin=a;}
  void SetProdTrackTPCNclsRatioMin(Double_t a){fProdTrackTPCNclsRatioMin=a;}
  void SetProdUseAODFilterBit(Bool_t a){fProdUseAODFilterBit=a;}
  void SetProdAODFilterBit(Int_t a){fProdAODFilterBit=a;}
  void SetProdRejectTrackWithShared(Bool_t a){fProdRejectTrackWithShared=a;}
  void SetProdV0MassTolK0s(Double_t a){fProdV0MassTolK0s=a;}
  void SetProdV0MassRejLambda(Double_t a){fProdV0MassRejLambda=a;}
  void SetProdV0MassRejPhoton(Double_t a){fProdV0MassRejPhoton=a;}
  void SetProdV0PtMin(Double_t a){fProdV0PtMin=a;}
  void SetProdV0CosPointingAngleToPrimVtxMin(Double_t a){fProdV0CosPointingAngleToPrimVtxMin=a;}
  void SetProdV0DcaDaughtersMax(Double_t a){fProdV0DcaDaughtersMax=a;}
  void SetProdV0DaughterEtaRange(Double_t a){fProdV0DaughterEtaRange=a;}
  void SetProdV0DaughterPtMin(Double_t a){fProdV0DaughterPtMin=a;}
  void SetProdV0DaughterTPCClusterMin(Double_t a){fProdV0DaughterTPCClusterMin=a;}
  void SetProdV0EtaRange(Double_t a, Double_t b){fProdV0EtaMin=a;fProdV0EtaMax=b;}
  void SetProdV0RapRange(Double_t a, Double_t b){fProdV0RapMin=a;fProdV0RapMax=b;}
  void SetProdRoughMassTol(Double_t a){fProdRoughMassTol=a;}
  void SetProdRoughPtMin(Double_t a){fProdRoughPtMin=a;}
	void SetTagV0MassTol(Double_t a){fTagV0MassTol=a;}
  
  Int_t GetProdTrackTPCNclsPIDMin(){return fProdTrackTPCNclsPIDMin;}
  Double_t GetProdTrackTPCNclsRatioMin(){return fProdTrackTPCNclsRatioMin;}
  Bool_t   GetProdUseAODFilterBit(){return fProdUseAODFilterBit;}
  Int_t   GetProdAODFilterBit(){return fProdAODFilterBit;}
  Bool_t   GetProdRejectTrackWithShared(){return fProdRejectTrackWithShared;}
  Double_t GetProdV0MassTolK0s(){return fProdV0MassTolK0s;}
  Double_t GetProdV0MassRejLambda(){return fProdV0MassRejLambda;}
  Double_t GetProdV0MassRejPhoton(){return fProdV0MassRejPhoton;}
  Double_t GetProdV0PtMin(){return fProdV0PtMin;}
  Double_t GetProdV0CosPointingAngleToPrimVtxMin(){return fProdV0CosPointingAngleToPrimVtxMin;}
  Double_t GetProdV0DcaDaughtersMax(){return fProdV0DcaDaughtersMax;}
  Double_t GetProdV0DaughterEtaRange(){return fProdV0DaughterEtaRange;}
  Double_t GetProdV0DaughterPtMin(){return fProdV0DaughterPtMin;}
  Double_t GetProdV0DaughterTPCClusterMin(){return fProdV0DaughterTPCClusterMin;}
  void GetProdV0EtaRange(Double_t &a, Double_t &b){a=fProdV0EtaMin;b=fProdV0EtaMax;}
  void GetProdV0RapRange(Double_t &a, Double_t &b){a=fProdV0RapMin;b=fProdV0RapMax;}
  Double_t GetProdRoughMassTol(){return fProdRoughMassTol;}
  Double_t GetProdRoughPtMin(){return fProdRoughPtMin;}
  Double_t GetTagV0MassTol(){return fTagV0MassTol;}
  
  Double_t CalculateLcCosPAXY(AliAODRecoDecayHF *obj);

	void SetMixingWeights(Int_t nbinpr, Double_t *bins_pr, Int_t nbink0s, Double_t *bins_k0s, Double_t *p0val, Double_t *p1val, Double_t *p2val, Double_t *p3val);
	Double_t GetMixingWeight(Double_t dphi, Double_t deta, Double_t pt_pr, Double_t pt_k0s);

 protected:
	
 private:

  EPIDStrategy fPIDStrategy;        /// PID strategy
  Double_t fCombinedPIDThreshold;   /// Threshold used in  IsSelectedCombinedPID
  Bool_t   fUseOnTheFlyV0;          /// Flag to check if we use on-the-fly v0
  Double_t fBzkG; ///B field
  Double_t fPrimVert[3];///Primary vertex
  
  Int_t fProdTrackTPCNclsPIDMin;      /// Min. Number of TPC PID cluster
  Double_t fProdTrackTPCNclsRatioMin;      /// Min. Number of TPC PID cluster
  Bool_t   fProdUseAODFilterBit;    /// Flag for AOD filter Bit used before object creation
  Int_t    fProdAODFilterBit;    /// AOD filter Bit used before object creation
  Bool_t   fProdRejectTrackWithShared;    /// Flag to Reject tracks with shared clusters
  Double_t fProdV0MassTolK0s;       /// K0s mass selection  used before object creation
  Double_t fProdV0MassRejLambda;       /// lambda mass rejection  used before object creation
  Double_t fProdV0MassRejPhoton;       /// photon mass rejection  used before object creation
  Double_t fProdV0PtMin;            /// Minimum K0s pT used before object creation
  Double_t fProdV0CosPointingAngleToPrimVtxMin; /// V0 pointing angle used before object creation
  Double_t fProdV0DcaDaughtersMax;  /// Max DCA between V0 daughters used before object creation
  Double_t fProdV0DaughterEtaRange; /// V0Daughter eta range used before object creation
  Double_t fProdV0DaughterPtMin;    /// V0 Daughter pT min used before object creation
  Double_t fProdV0DaughterTPCClusterMin; /// V0 daughter Minimum TPC cluster pT used before object creation
	Double_t fProdV0EtaMin; /// Minimum eta of cascade
	Double_t fProdV0EtaMax; /// Maximum eta of cascade
	Double_t fProdV0RapMin; /// Minimum rapidity of cascade
	Double_t fProdV0RapMax; /// Maximum rapidity of cascade
  Double_t fProdRoughMassTol;       /// Mass cut for Lc used before object creation
  Double_t fProdRoughPtMin;         /// pT cut for Lc used before object creation

	Int_t fNWeightingProtonBinLimits; ///Number of bins for proton
	Double_t *fWeightingProtonBins; //[fNWeightingProtonBinLimits] ptbin of proton weighting
	Int_t fNWeightingK0sBinLimits; ///Number of bins for k0s
	Double_t *fWeightingK0sBins; //[fNWeightingK0sBinLimits] ptbin of k0s weighting
	Int_t fNWeightingBins; ///Number of bins for mixing weight should be proton x k0s
	Double_t *fWeight_p0; //[fNWeightingBins] p0 values
	Double_t *fWeight_p1; //[fNWeightingBins] p1 values
	Double_t *fWeight_p2; //[fNWeightingBins] p2 values
	Double_t *fWeight_p3; //[fNWeightingBins] p3 values

	Double_t fTagV0MassTol; //V0 tagging tolorance
  
  /// \cond CLASSIMP
  ClassDef(AliRDHFCutsLctopK0sfromAODtracks,4);
  /// \endcond
};

#endif
