#ifndef ALIRDHFCUTSKFP_H
#define ALIRDHFCUTSKFP_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//***********************************************************
/// \class Class AliRDHFCutsKFP
/// \brief class for cuts on AOD reconstructed Xic0 -> Xi pi
/// \author Jianhui Zhu <zjh@mail.ccnu.edu.cn>, Central China Normal University & GSI Helmholtz Centre for Heavy Ion Research
/// \date Apr 26, 2019
//***********************************************************

/* $Id$ */

#include "AliRDHFCuts.h"

class AliRDHFCutsKFP : public AliRDHFCuts
{
 public:

  enum EPIDStrategy{
    kNSigmaCuts,
    kCombinedCuts
  };

  AliRDHFCutsKFP(const char* name="CutsXicZerotoXiPi");
  virtual ~AliRDHFCutsKFP();
  AliRDHFCutsKFP(const AliRDHFCutsKFP& source);
  AliRDHFCutsKFP& operator=(const AliRDHFCutsKFP& source);

  using AliRDHFCuts::GetCutVarsForOpt;
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters);

  using AliRDHFCuts::IsSelected;
  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel);
  using AliRDHFCuts::IsSelectedPID;
  virtual Int_t IsSelectedPID(AliAODRecoDecayHF* obj);
  Int_t IsSelectedCombinedPID(AliAODRecoDecayHF* obj);
  Double_t GetPionProbabilityTPCTOF(AliAODTrack *trk);

  void SetPIDStrategy(EPIDStrategy pidStrategy){fPIDStrategy=pidStrategy;}
  EPIDStrategy GetPIDStrategy() const {return fPIDStrategy;}
  void SetCombinedPIDThreshold(Double_t a){fCombinedPIDThreshold=a;}
  void SetUseXic0PID(Bool_t a){fUseXic0PID=a;}
  Double_t GetCombinedPIDThreshold(){return fCombinedPIDThreshold;}
  void SetPidPiFromXic0(AliAODPidHF* pidPion) {
      if(fPidObjPiFromXic0) delete fPidObjPiFromXic0;
      fPidObjPiFromXic0 = new AliAODPidHF(*pidPion);
      }

  void SetPidPiFromXi(AliAODPidHF* pidPion) { 
      if(fPidObjPiFromXi) delete fPidObjPiFromXi;
      fPidObjPiFromXi=new AliAODPidHF(*pidPion);
      }

  void SetPidKaFromOmega(AliAODPidHF* pidKaon) { 
      if(fPidObjKaFromOmega) delete fPidObjKaFromOmega;
      fPidObjKaFromOmega=new AliAODPidHF(*pidKaon);
      }

  void SetPidPiFromV0(AliAODPidHF* pidPion) { 
      if(fPidObjPiFromV0) delete fPidObjPiFromV0;
      fPidObjPiFromV0=new AliAODPidHF(*pidPion);
      }

  void SetPidPrFromV0(AliAODPidHF* pidProton) { 
      if(fPidObjPrFromV0) delete fPidObjPrFromV0;
      fPidObjPrFromV0=new AliAODPidHF(*pidProton);
      }

  AliAODPidHF* GetPidPiFromXic0() const {return fPidObjPiFromXic0;}
  AliAODPidHF* GetPidPiFromXi() const {return fPidObjPiFromXi;}
  AliAODPidHF* GetPidKaFromOmega() const {return fPidObjKaFromOmega;}
  AliAODPidHF* GetPidPiFromV0() const {return fPidObjPiFromV0;}
  AliAODPidHF* GetPidPrFromV0() const {return fPidObjPrFromV0;}

  Bool_t SingleTrkCuts(AliAODTrack *trk);
  Bool_t SingleCascadeCuts(AliAODcascade *casc, Double_t *vert, Bool_t anaOmegacZero);
  Bool_t SingleV0LambdaTotCuts(AliAODv0 *v0);
  Bool_t SingleCascCuts(AliAODcascade *casc);
  Bool_t LambdaPIDCuts(AliAODv0 *v0);
  Bool_t AntiLambdaPIDCuts(AliAODv0 *v0);
  Bool_t SinglePionPoolCuts(AliAODTrack *trk);
  Bool_t PassedTrackQualityCuts_PrimaryPion(AliAODTrack *trk);
  Bool_t PassedTrackQualityCuts_SecondaryPion(AliAODTrack *trk);
  Bool_t SingleCascadeCutsRef(AliAODcascade *casc, Double_t *vert, Bool_t anaOmegaZero);
  Bool_t SelectWithRoughCuts(AliAODcascade *casc, AliAODTrack *trk1);

  void SetPtMinXic0(Double_t a){fPtMinXic0=a;}
  void SetPtMinPiFromXic0(Double_t a){fPtMinPiFromXic0=a;}
  void SetPtMinPiFromXi(Double_t a){fPtMinPiFromXi=a;}
  void SetProdTrackEtaRange(Double_t a){fProdTrackEtaRange=a;}
  void SetProdUseAODFilterBit(Bool_t a){fProdUseAODFilterBit=a;}
  void SetProdMassTolLambda(Double_t a){fProdMassTolLambda=a;}
  void SetProdMassTolK0S(Double_t a){fProdMassTolK0S=a;}
  void SetProdMassTolXi(Double_t a){fProdMassTolXi=a;}
  void SetProdMassTolXic0(Double_t a){fProdMassTolXic0=a;}
  void SetProdMassTolOmega(Double_t a){fProdMassTolOmega=a;}
  void SetProdMassRejOmega(Double_t a){fProdMassRejOmega=a;}
  void SetProdRfidMinV0(Double_t a){fProdRfidMinV0=a;}
  void SetProdRfidMaxV0(Double_t a){fProdRfidMaxV0=a;}
  void SetProdRfidMinXi(Double_t a){fProdRfidMinXi=a;}
  void SetProdRfidMaxXi(Double_t a){fProdRfidMaxXi=a;}
  void SetProdCascProperDecayLengthMax(Double_t a){fProdCascProperDecayLengthMax=a;}
  void SetProdDcaXiDaughtersMax(Double_t a){fProdDcaXiDaughtersMax=a;}
  void SetProdDcaV0DaughtersMax(Double_t a){fProdDcaV0DaughtersMax=a;}
  void SetProdDcaBachToPrimVertexMin(Double_t a){fProdDcaBachToPrimVertexMin=a;}
  void SetProdDcaV0ToPrimVertexMin(Double_t a){fProdDcaV0ToPrimVertexMin=a;}
  void SetProdDcaV0PrToPrimVertexMin(Double_t a){fProdDcaV0PrToPrimVertexMin=a;}
  void SetProdDcaV0PiToPrimVertexMin(Double_t a){fProdDcaV0PiToPrimVertexMin=a;}
  void SetProdXiCosineOfPoiningAngleMin(Double_t a){fProdXiCosineOfPoiningAngleMin=a;}
  void SetProdV0CosineOfPoiningAngleXiMin(Double_t a){fProdV0CosineOfPoiningAngleXiMin=a;}
  void SetProdRoughMassTol(Double_t a){fProdRoughMassTol=a;}
  void SetProdRoughPtMin(Double_t a){fProdRoughPtMin=a;}
  void SetProdLikeSignDcaMax(Double_t a){fProdLikeSignDcaMax=a;}
  void SetProdCascNTPCClustersMin(Double_t a){fProdCascNTPCClustersMin=a;}
  void SetProdChi2TPCV0PrMax(Double_t a) {fProdChi2TPCV0PrMax=a;}
  void SetProdChi2TPCV0PiMax(Double_t a) {fProdChi2TPCV0PiMax=a;}
  void SetKFPLam_Chi2geoMax(Double_t a) {fKFPLam_Chi2geoMax=a;}
  void SetKFPLam_Chi2topoMax(Double_t a) {fKFPLam_Chi2topoMax=a;}
  void SetKFPLam_lDeltalMin(Double_t a) {fKFPLam_lDeltalMin=a;}
  void SetKFPXi_Chi2geoMax(Double_t a) {fKFPXi_Chi2geoMax=a;}
  void SetKFPXi_Chi2topoMax(Double_t a) {fKFPXi_Chi2topoMax=a;}
  void SetKFPXi_lDeltalMin(Double_t a) {fKFPXi_lDeltalMin=a;}
  void SetKFPXic0_Chi2geoMax(Double_t a) {fKFPXic0_Chi2geoMax=a;}
  void SetProdTrackTPCNCrossedRowsMin(Int_t a) {fProdTrackTPCNCrossedRowsMin=a;}
  void SetProdTrackTPCNCrossedRowsRatioMin(Double_t a) {fProdTrackTPCNCrossedRowsRatioMin=a;}
  void SetProdTrackTPCsignalNMin(Int_t a) {fProdTrackTPCsignalNMin=a;}
  void SetPriTrackChi2perNDFMax(Double_t a) {fPriTrackChi2perNDFMax=a;}
  void SetPriTrackITSNclsMin(Int_t a) {fPriTrackITSNclsMin=a;}

  Double_t GetPtMinXic0(){return fPtMinXic0;}
  Double_t GetPtMinPiFromXic0(){return fPtMinPiFromXic0;}
  Double_t GetPtMinPiFromXi(){return fPtMinPiFromXi;}
  Double_t GetProdTrackEtaRange(){return fProdTrackEtaRange;}
  Bool_t   GetProdUseAODFilterBit(){return fProdUseAODFilterBit;}
  Double_t GetProdMassTolLambda(){return fProdMassTolLambda;}
  Double_t GetProdMassTolK0S(){return fProdMassTolK0S;}
  Double_t GetProdMassTolXi(){return fProdMassTolXi;}
  Double_t GetProdMassTolXic0(){return fProdMassTolXic0;}
  Double_t GetProdMassTolOmega(){return fProdMassTolOmega;}
  Double_t GetProdMassRejOmega(){return fProdMassRejOmega;}
  Double_t GetProdRfidMinV0(){return fProdRfidMinV0;}
  Double_t GetProdRfidMaxV0(){return fProdRfidMaxV0;}
  Double_t GetProdRfidMinXi(){return fProdRfidMinXi;}
  Double_t GetProdRfidMaxXi(){return fProdRfidMaxXi;}
  Double_t GetProdCascProperDecayLengthMax(){return fProdCascProperDecayLengthMax;}
  Double_t GetProdDcaXiDaughtersMax(){return fProdDcaXiDaughtersMax;}
  Double_t GetProdDcaV0DaughtersMax(){return fProdDcaV0DaughtersMax;}
  Double_t GetProdDcaBachToPrimVertexMin(){return fProdDcaBachToPrimVertexMin;}
  Double_t GetProdDcaV0ToPrimVertexMin(){return fProdDcaV0ToPrimVertexMin;}
  Double_t GetProdDcaV0PrToPrimVertexMin(){return fProdDcaV0PrToPrimVertexMin;}
  Double_t GetProdDcaV0PiToPrimVertexMin(){return fProdDcaV0PiToPrimVertexMin;}
  Double_t GetProdXiCosineOfPoiningAngleMin(){return fProdXiCosineOfPoiningAngleMin;}
  Double_t GetProdV0CosineOfPoiningAngleXiMin(){return fProdV0CosineOfPoiningAngleXiMin;}
  Double_t GetProdRoughMassTol(){return fProdRoughMassTol;}
  Double_t GetProdRoughPtMin(){return fProdRoughPtMin;}
  Double_t GetProdLikeSignDcaMax(){return fProdLikeSignDcaMax;}
  Double_t GetProdCascNTPCClustersMin(){return fProdCascNTPCClustersMin;}
  Double_t GetProdChi2TPCV0PrMax() {return fProdChi2TPCV0PrMax;}
  Double_t GetProdChi2TPCV0PiMax() {return fProdChi2TPCV0PiMax;}
  Double_t GetKFPLam_Chi2geoMax() {return fKFPLam_Chi2geoMax;}
  Double_t GetKFPLam_Chi2topoMax() {return fKFPLam_Chi2topoMax;}
  Double_t GetKFPLam_lDeltalMin() {return fKFPLam_lDeltalMin;}
  Double_t GetKFPXi_Chi2geoMax() {return fKFPXi_Chi2geoMax;}
  Double_t GetKFPXi_Chi2topoMax() {return fKFPXi_Chi2topoMax;}
  Double_t GetKFPXi_lDeltalMin() {return fKFPXi_lDeltalMin;}
  Double_t GetKFPXic0_Chi2geoMax() {return fKFPXic0_Chi2geoMax;}
  Int_t    GetProdTrackTPCNCrossedRowsMin() {return fProdTrackTPCNCrossedRowsMin;}
  Double_t GetProdTrackTPCNCrossedRowsRatioMin() {return fProdTrackTPCNCrossedRowsRatioMin;}
  Int_t    GetProdTrackTPCsignalNMin() {return fProdTrackTPCsignalNMin;}
  Double_t GetPriTrackChi2perNDFMax() {return fPriTrackChi2perNDFMax;}
  Int_t    GetPriTrackITSNclsMin() {return fPriTrackITSNclsMin;}

  void useSetNPtBins(Int_t nptBins){SetNPtBins(nptBins);}
 protected:
	
 private:

  EPIDStrategy fPIDStrategy;        /// PID Strategy
  Double_t fCombinedPIDThreshold;   /// PID threshold used in IsSelectedCombinedPID
  Bool_t fUseXic0PID;                   /// Use PID or not
  AliAODPidHF *fPidObjPiFromXic0;       /// PID object for Xic0-pion
  AliAODPidHF *fPidObjPiFromXi;         /// PID object for cascade-pion
  AliAODPidHF *fPidObjKaFromOmega;      /// PID object for cascade-kaon
  AliAODPidHF *fPidObjPrFromV0;         /// PID object for V0-proton
  AliAODPidHF *fPidObjPiFromV0;         /// PID object for V0-pion

  Double_t fPtMinXic0;               /// Minimum pT of Xic0
  Double_t fPtMinPiFromXic0;         /// Minimum Bachelor pT of pion from Xic0 decay
  Double_t fPtMinPiFromXi;         /// Minimum Bachelor pT of pion from Xi decay
  Double_t fProdTrackEtaRange;      /// Bachelor Eta range
  Bool_t   fProdUseAODFilterBit;    /// Use AODfilterBit or not
  Double_t fProdMassTolLambda;      /// Tolerance of Lambda mass from PDG value
  Double_t fProdMassTolK0S;         /// Tolerance of K0S mass from PDG value
  Double_t fProdMassTolXi;          /// Tolerance of Xi mass from PDG value
  Double_t fProdMassTolXic0;          /// Tolerance of Xic0 mass from PDG value
  Double_t fProdMassTolOmega;          /// Tolerance of Omega mass from PDG value
  Double_t fProdMassRejOmega;          /// Rejection range of Omega mass from PDG value
  Double_t fProdRfidMinV0;          /// Minimum Decay vertex of V0
  Double_t fProdRfidMaxV0;          /// Max Decay vertex of V0
  Double_t fProdRfidMinXi;          /// Minimum Decay vertex of Xi
  Double_t fProdRfidMaxXi;          /// Max Decay vertex of Xi
  Double_t fProdCascProperDecayLengthMax;        /// mL/p of cascade
  Double_t fProdDcaXiDaughtersMax;  /// Max Dca between Xi daughters
  Double_t fProdDcaV0DaughtersMax;  /// Max Dca between V0 daughters
  Double_t fProdDcaBachToPrimVertexMin;  /// Min Dca between Bachelor and PV
  Double_t fProdDcaV0ToPrimVertexMin;  /// Min Dca between v0 and PV
  Double_t fProdDcaV0PrToPrimVertexMin;  /// Min Dca between v0-proton and PV
  Double_t fProdDcaV0PiToPrimVertexMin;  /// Min Dca between v0-pion and PV
  Double_t fProdXiCosineOfPoiningAngleMin;  /// Min Xi cos pointing angle  to PV
  Double_t fProdV0CosineOfPoiningAngleXiMin;  // /Min V0 cos pointing angle  to Xi vertex
  Double_t fProdCascNTPCClustersMin;         /// Minimum number of TPC clusters
  Double_t fProdChi2TPCV0PrMax;               /// Max. chi2 of TPC clusters for v0-proton
  Double_t fProdChi2TPCV0PiMax;              /// Max. chi2 of TPC clusters for v0-pion
  Double_t fProdLikeSignDcaMax;     /// Maximum DCA of pions
  Double_t fProdRoughMassTol;       /// Tolerance of Xic mass from PDG value
  Double_t fProdRoughPtMin;         /// Minimum pT of Xic
  Int_t    fProdTrackTPCNCrossedRowsMin; /// Minimum number of TPC crossed rows
  Double_t fProdTrackTPCNCrossedRowsRatioMin; /// Minimum ratio of TPC crossed rows / findable ratio
  Int_t    fProdTrackTPCsignalNMin; /// Minimum number of TPC clusters for dE/dx
  Double_t fPriTrackChi2perNDFMax; /// Max. chi2/NDF of momentum fit
  Int_t    fPriTrackITSNclsMin; /// Min. number of points in ITS

  Double_t fKFPLam_Chi2geoMax;          /// chi2/ndf cut of lambda reconstruction from KFParticle
  Double_t fKFPLam_Chi2topoMax;          /// chi2/ndf cut of lambda reconstruction from KFParticle
  Double_t fKFPLam_lDeltalMin;          /// l/Deltal cut of lambda reconstruction from KFParticle
  Double_t fKFPXi_Chi2geoMax;         /// chi2/ndf cut of Xi- reconstruction from KFParticle
  Double_t fKFPXi_Chi2topoMax;         /// chi2/ndf cut of Xi- reconstruction from KFParticle
  Double_t fKFPXi_lDeltalMin;          /// l/Deltal cut of lambda reconstruction from KFParticle
  Double_t fKFPXic0_Chi2geoMax;         /// chi2/ndf cut of Xic0 reconstruction from KFParticle


  /// \cond CLASSIMP
  ClassDef(AliRDHFCutsKFP, 3);
  /// \endcond
};

#endif
