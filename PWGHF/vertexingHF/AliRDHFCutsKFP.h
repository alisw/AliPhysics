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

#include "AliAODRecoCascadeHF.h"
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
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters) {
                           return GetCutVarsForOpt(d,vars,nvars,pdgdaughters,0x0); }
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters, AliAODEvent *aod);

  using AliRDHFCuts::IsSelected;
  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel)
                        {return IsSelected(obj,selectionLevel,0);}
  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel, AliAODEvent* aod);
  using AliRDHFCuts::IsSelectedPID;
  virtual Int_t IsSelectedPID(AliAODRecoDecayHF* obj);
  Int_t IsSelectedCombinedPID(AliAODRecoDecayHF* obj);
  Double_t GetPionProbabilityTPCTOF(AliAODTrack *trk);

  void SetPIDStrategy(EPIDStrategy pidStrategy){fPIDStrategy=pidStrategy;}
  EPIDStrategy GetPIDStrategy() const {return fPIDStrategy;}
  void SetCombinedPIDThreshold(Double_t a){fCombinedPIDThreshold=a;}
  void SetUseXic0PID(Bool_t a){fUseXic0PID=a;}
  void SetUseXicPlusPID(Bool_t a){fUseXicPlusPID=a;}
  void SetUseXiccPPPID(Bool_t a){fUseXiccPPPID=a;}
  void SetUseLcPID(Bool_t a){fUseLcPID=a;}
  Double_t GetCombinedPIDThreshold(){return fCombinedPIDThreshold;}
  void SetPidDaughter(AliAODPidHF* pid) {
      if(fPidObjDau) delete fPidObjDau;
      fPidObjDau = new AliAODPidHF(*pid);
      }
  void SetPidPiFromXic0(AliAODPidHF* pidPion) {
      if(fPidObjPiFromXic0) delete fPidObjPiFromXic0;
      fPidObjPiFromXic0 = new AliAODPidHF(*pidPion);
      }

  void SetPidPiFromXiccPP(AliAODPidHF* pidPion) {
      if(fPidObjPiFromXiccPP) delete fPidObjPiFromXiccPP;
      fPidObjPiFromXiccPP = new AliAODPidHF(*pidPion);
      }

  void SetPidPiFromXicPlus(AliAODPidHF* pidPion) {
      if(fPidObjPiFromXicPlus) delete fPidObjPiFromXicPlus;
      fPidObjPiFromXicPlus = new AliAODPidHF(*pidPion);
      }

  void SetPidPiFromXi(AliAODPidHF* pidPion) { 
      if(fPidObjPiFromXi) delete fPidObjPiFromXi;
      fPidObjPiFromXi=new AliAODPidHF(*pidPion);
      }

  void SetPidPiFromV0(AliAODPidHF* pidPion) { 
      if(fPidObjPiFromV0) delete fPidObjPiFromV0;
      fPidObjPiFromV0=new AliAODPidHF(*pidPion);
      }

  void SetPidPrFromV0(AliAODPidHF* pidProton) { 
      if(fPidObjPrFromV0) delete fPidObjPrFromV0;
      fPidObjPrFromV0=new AliAODPidHF(*pidProton);
      }
  void SetWeightFunction(TF1 *weight) {
      if (fWeight) delete fWeight;
      fWeight=new TF1 (*weight);
      }
  void SetWeightFunctionUp(TF1 *weight_up) {
      if (fWeight_up) delete fWeight_up;
      fWeight_up=new TF1 (*weight_up);
      }
  void SetWeightFunctionDw(TF1 *weight_dw) {
      if (fWeight_dw) delete fWeight_dw;
      fWeight_dw=new TF1 (*weight_dw);
      }

  AliAODPidHF* GetPidDaughter() const {return fPidObjDau;}
  AliAODPidHF* GetPidPiFromXic0() const {return fPidObjPiFromXic0;}
  AliAODPidHF* GetPidPiFromXiccPP() const {return fPidObjPiFromXiccPP;}
  AliAODPidHF* GetPidPiFromXicPlus() const {return fPidObjPiFromXicPlus;}
  AliAODPidHF* GetPidPiFromXi() const {return fPidObjPiFromXi;}
  AliAODPidHF* GetPidPiFromV0() const {return fPidObjPiFromV0;}
  AliAODPidHF* GetPidPrFromV0() const {return fPidObjPrFromV0;}
  TF1* GetWeightFunction()   const {return fWeight;}
  TF1* GetWeightFunctionUp() const {return fWeight_up;}
  TF1* GetWeightFunctionDw() const {return fWeight_dw;}

  Bool_t SingleTrkCuts(AliAODTrack *trk);
  Bool_t SingleCascadeCuts(AliAODcascade *casc, Double_t *vert, Bool_t anaOmegacZero);
  Bool_t SingleV0LambdaTotCuts(AliAODv0 *v0);
  Bool_t SingleCascCuts(AliAODcascade *casc, Bool_t IsAnaOmegac0);
  Bool_t LambdaPIDCuts(AliAODv0 *v0);
  Bool_t AntiLambdaPIDCuts(AliAODv0 *v0);
  Bool_t SinglePionPoolCuts(AliAODTrack *trk);
  Bool_t PassedTrackQualityCuts_PrimaryPion(AliAODTrack *trk);
  Bool_t PassedTrackQualityCuts_SecondaryPion(AliAODTrack *trk);
  Bool_t SingleCascadeCutsRef(AliAODcascade *casc, Double_t *vert, Bool_t anaOmegaZero);
  Bool_t SelectWithRoughCuts(AliAODcascade *casc, AliAODTrack *trk1);
  Bool_t PreSelForLc2pKs0(AliAODRecoCascadeHF *Lc2pKs0);
  Bool_t PreSelForLc2Lpi(AliAODRecoCascadeHF *Lc2Lpi);

  void SetPtMinLc(Double_t a){fPtMinLc=a;}
  void SetPtMinPrFromLc(Double_t a){fPtMinPrFromLc=a;}
  void SetPtMinPiFromLc(Double_t a){fPtMinPiFromLc=a;}
  void SetPtMinXic0(Double_t a){fPtMinXic0=a;}
  void SetPtMinXiccPP(Double_t a){fPtMinXiccPP=a;}
  void SetPtMinXicPlus(Double_t a){fPtMinXicPlus=a;}
  void SetPtMinPiFromPrimary(Double_t a){fPtMinPiFromPrimary=a;}
  void SetPtMinPiFromXi(Double_t a){fPtMinPiFromXi=a;}
  void SetPtMinPiFromXic0ForML(Double_t a){fPtMinPiFromXic0ForML=a;}
  void SetProdTrackEtaRange(Double_t a){fProdTrackEtaRange=a;}
  void SetProdUseAODFilterBit(Bool_t a){fProdUseAODFilterBit=a;}
  void SetProdMassTolLc(Double_t a){fProdMassTolLc=a;}
  void SetProdMassTolLambda(Double_t a){fProdMassTolLambda=a;}
  void SetProdMassRejLambda(Double_t a){fProdMassRejLambda=a;}
  void SetProdMassTolKs0(Double_t a){fProdMassTolKs0=a;}
  void SetProdMassRejKs0(Double_t a){fProdMassRejKs0=a;}
  void SetProdMassTolXi(Double_t a){fProdMassTolXi=a;}
  void SetProdMassTolXic0(Double_t a){fProdMassTolXic0=a;}
  void SetProdMassTolXiccPP(Double_t a){fProdMassTolXiccPP=a;}
  void SetProdMassTolXicPlus(Double_t a){fProdMassTolXicPlus=a;}
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
  void SetKFPKs0_Chi2geoMax(Double_t a) {fKFPKs0_Chi2geoMax=a;}
  void SetKFPKs0_lDeltalMin(Double_t a) {fKFPKs0_lDeltalMin=a;}
  void SetKFPKs0_Chi2topoMax(Double_t a) {fKFPKs0_Chi2topoMax=a;}
  void SetKFPLc_Chi2geoMax(Double_t a) {fKFPLc_Chi2geoMax=a;}
  void SetKFPLam_Chi2geoMax(Double_t a) {fKFPLam_Chi2geoMax=a;}
  void SetKFPLam_Chi2topoMin(Double_t a) {fKFPLam_Chi2topoMin=a;}
  void SetKFPLam_lDeltalMin(Double_t a) {fKFPLam_lDeltalMin=a;}
  void SetKFPXi_Chi2geoMax(Double_t a) {fKFPXi_Chi2geoMax=a;}
  void SetKFPXi_Chi2topoMax(Double_t a) {fKFPXi_Chi2topoMax=a;}
  void SetKFPXi_lDeltalMin(Double_t a) {fKFPXi_lDeltalMin=a;}
  void SetKFPXic0_Chi2geoMax(Double_t a) {fKFPXic0_Chi2geoMax=a;}
  void SetKFPXiccPP_Chi2geoMax(Double_t a) {fKFPXiccPP_Chi2geoMax=a;}
  void SetKFPXicPlus_Chi2geoMax(Double_t a) {fKFPXicPlus_Chi2geoMax=a;}
  void SetProdTrackTPCNCrossedRowsMin(Int_t a) {fProdTrackTPCNCrossedRowsMin=a;}
  void SetProdTrackTPCNCrossedRowsRatioMin(Double_t a) {fProdTrackTPCNCrossedRowsRatioMin=a;}
  void SetProdTrackTPCsignalNMin(Int_t a) {fProdTrackTPCsignalNMin=a;}
  void SetPriTrackChi2perNDFMax(Double_t a) {fPriTrackChi2perNDFMax=a;}
  void SetPriTrackITSNclsMin(Int_t a) {fPriTrackITSNclsMin=a;}
  void SetCut_nSigmaTPC_PiFromXic0_or_PiFromOmegac0(Double_t a) {fCut_nSigmaTPC_PiFromXic0_or_PiFromOmegac0=a;}
  void SetCut_nSigmaTPC_PiFromXi_or_KaFromOmega(Double_t a) {fCut_nSigmaTPC_PiFromXi_or_KaFromOmega=a;}
  void SetCut_nSigmaTPC_PrFromLam(Double_t a) {fCut_nSigmaTPC_PrFromLam=a;}
  void SetCut_nSigmaTPC_PiFromLam(Double_t a) {fCut_nSigmaTPC_PiFromLam=a;}
  void SetCut_nSigmaTOF_PiFromXic0_or_PiFromOmegac0(Double_t a) {fCut_nSigmaTOF_PiFromXic0_or_PiFromOmegac0=a;}
  void SetCut_nSigmaTOF_PiFromXi_or_KaFromOmega(Double_t a) {fCut_nSigmaTOF_PiFromXi_or_KaFromOmega=a;}
  void SetCut_nSigmaTOF_PrFromLam(Double_t a) {fCut_nSigmaTOF_PrFromLam=a;}
  void SetCut_nSigmaTOF_PiFromLam(Double_t a) {fCut_nSigmaTOF_PiFromLam=a;}
  void SetCut_Armenteros(Double_t a) {fCut_Armenteros=a;}

  Double_t GetPtMinLc(){return fPtMinLc;}
  Double_t GetPtMinPrFromLc(){return fPtMinPrFromLc;}
  Double_t GetPtMinPiFromLc(){return fPtMinPiFromLc;}
  Double_t GetPtMinXic0(){return fPtMinXic0;}
  Double_t GetPtMinXiccPP(){return fPtMinXiccPP;}
  Double_t GetPtMinXicPlus(){return fPtMinXicPlus;}
  Double_t GetPtMinPiFromPrimary(){return fPtMinPiFromPrimary;}
  Double_t GetPtMinPiFromXi(){return fPtMinPiFromXi;}
  Double_t GetPtMinPiFromXic0ForML(){return fPtMinPiFromXic0ForML;}
  Double_t GetProdTrackEtaRange(){return fProdTrackEtaRange;}
  Bool_t   GetProdUseAODFilterBit(){return fProdUseAODFilterBit;}
  Double_t GetProdMassTolLc(){return fProdMassTolLc;}
  Double_t GetProdMassTolLambda(){return fProdMassTolLambda;}
  Double_t GetProdMassRejLambda(){return fProdMassRejLambda;}
  Double_t GetProdMassTolKs0(){return fProdMassTolKs0;}
  Double_t GetProdMassRejKs0(){return fProdMassRejKs0;}
  Double_t GetProdMassTolXi(){return fProdMassTolXi;}
  Double_t GetProdMassTolXic0(){return fProdMassTolXic0;}
  Double_t GetProdMassTolXiccPP(){return fProdMassTolXiccPP;}
  Double_t GetProdMassTolXicPlus(){return fProdMassTolXicPlus;}
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
  Double_t GetKFPKs0_Chi2geoMax() {return fKFPKs0_Chi2geoMax;}
  Double_t GetKFPKs0_lDeltalMin() {return fKFPKs0_lDeltalMin;}
  Double_t GetKFPKs0_Chi2topoMax() {return fKFPKs0_Chi2topoMax;}
  Double_t GetKFPLc_Chi2geoMax() {return fKFPLc_Chi2geoMax;}
  Double_t GetKFPLam_Chi2geoMax() {return fKFPLam_Chi2geoMax;}
  Double_t GetKFPLam_Chi2topoMin() {return fKFPLam_Chi2topoMin;}
  Double_t GetKFPLam_lDeltalMin() {return fKFPLam_lDeltalMin;}
  Double_t GetKFPXi_Chi2geoMax() {return fKFPXi_Chi2geoMax;}
  Double_t GetKFPXi_Chi2topoMax() {return fKFPXi_Chi2topoMax;}
  Double_t GetKFPXi_lDeltalMin() {return fKFPXi_lDeltalMin;}
  Double_t GetKFPXic0_Chi2geoMax() {return fKFPXic0_Chi2geoMax;}
  Double_t GetKFPXiccPP_Chi2geoMax() {return fKFPXiccPP_Chi2geoMax;}
  Double_t GetKFPXicPlus_Chi2geoMax() {return fKFPXicPlus_Chi2geoMax;}
  Int_t    GetProdTrackTPCNCrossedRowsMin() {return fProdTrackTPCNCrossedRowsMin;}
  Double_t GetProdTrackTPCNCrossedRowsRatioMin() {return fProdTrackTPCNCrossedRowsRatioMin;}
  Int_t    GetProdTrackTPCsignalNMin() {return fProdTrackTPCsignalNMin;}
  Double_t GetPriTrackChi2perNDFMax() {return fPriTrackChi2perNDFMax;}
  Int_t    GetPriTrackITSNclsMin() {return fPriTrackITSNclsMin;}
  Double_t GetCut_nSigmaTPC_PiFromXic0_or_PiFromOmegac0() {return fCut_nSigmaTPC_PiFromXic0_or_PiFromOmegac0;}
  Double_t GetCut_nSigmaTPC_PiFromXi_or_KaFromOmega() {return fCut_nSigmaTPC_PiFromXi_or_KaFromOmega;}
  Double_t GetCut_nSigmaTPC_PrFromLam() {return fCut_nSigmaTPC_PrFromLam;}
  Double_t GetCut_nSigmaTPC_PiFromLam() {return fCut_nSigmaTPC_PiFromLam;}
  Double_t GetCut_nSigmaTOF_PiFromXic0_or_PiFromOmegac0() {return fCut_nSigmaTOF_PiFromXic0_or_PiFromOmegac0;}
  Double_t GetCut_nSigmaTOF_PiFromXi_or_KaFromOmega() {return fCut_nSigmaTOF_PiFromXi_or_KaFromOmega;}
  Double_t GetCut_nSigmaTOF_PrFromLam() {return fCut_nSigmaTOF_PrFromLam;}
  Double_t GetCut_nSigmaTOF_PiFromLam() {return fCut_nSigmaTOF_PiFromLam;}
  Double_t GetCut_Armenteros() {return fCut_Armenteros;}

  void useSetNPtBins(Int_t nptBins){SetNPtBins(nptBins);}
 protected:
	
 private:

  EPIDStrategy fPIDStrategy;        /// PID Strategy
  Double_t fCombinedPIDThreshold;   /// PID threshold used in IsSelectedCombinedPID
  Bool_t fUseLcPID;                   /// Use PID or not for Lc
  Bool_t fUseXic0PID;                 /// Use PID or not for Xic0
  Bool_t fUseXiccPPPID;              /// Use PID or not for Xicc++
  Bool_t fUseXicPlusPID;              /// Use PID or not for XicPlus
  AliAODPidHF *fPidObjDau;              /// PID object for all daughter tracks
  AliAODPidHF *fPidObjPiFromXiccPP;     /// PID object for Xicc++ - pion
  AliAODPidHF *fPidObjPiFromXicPlus;    /// PID object for XicPlus-pion
  AliAODPidHF *fPidObjPiFromXic0;       /// PID object for Xic0-pion
  AliAODPidHF *fPidObjPiFromXi;         /// PID object for cascade-pion
  AliAODPidHF *fPidObjPrFromV0;         /// PID object for V0-proton
  AliAODPidHF *fPidObjPiFromV0;         /// PID object for V0-pion

  Double_t fPtMinLc;                 /// Minimum pT of Lc
  Double_t fPtMinPrFromLc;           /// Minimum pT of proton from Lc decay
  Double_t fPtMinXic0;               /// Minimum pT of Xic0
  Double_t fPtMinPiFromPrimary;      /// Minimum Bachelor pT of pion from primary vtx
  Double_t fPtMinPiFromXi;         /// Minimum Bachelor pT of pion from Xi decay
  Double_t fPtMinPiFromXic0ForML;   /// Minimum Bachelor pT of pion from Xic0 decay for ML
  Double_t fPtMinPiFromLc;           /// Minimum pT of pion from Lc decay
  Double_t fPtMinXiccPP;            /// Minimum pT of Xicc++
  Double_t fPtMinXicPlus;            /// Minimum pT of XicPlus
  Double_t fProdTrackEtaRange;      /// Bachelor Eta range
  Bool_t   fProdUseAODFilterBit;    /// Use AODfilterBit or not
  Double_t fProdMassTolLc;          /// Tolerance of Lc mass from PDG value
  Double_t fProdMassTolKs0;         /// Tolerance of Ks0 mass from PDG value
  Double_t fProdMassTolLambda;      /// Tolerance of Lambda mass from PDG value
  Double_t fProdMassTolXi;          /// Tolerance of Xi mass from PDG value
  Double_t fProdMassTolXic0;        /// Tolerance of Xic0 mass from PDG value
  Double_t fProdRfidMinV0;          /// Minimum Decay vertex of V0
  Double_t fProdRfidMaxV0;          /// Max Decay vertex of V0
  Double_t fProdRfidMinXi;          /// Minimum Decay vertex of Xi
  Double_t fProdRfidMaxXi;          /// Max Decay vertex of Xi
  Double_t fProdCascProperDecayLengthMax;        /// mL/p of cascade
  Double_t fProdMassRejLambda;      /// Rejection of Lambda mass from PDG value
  Double_t fProdMassRejKs0;         /// Rejection of Ks0 mass from PDG value
  Double_t fProdMassTolXiccPP;     /// Tolerance of Xicc++ mass from PDG value
  Double_t fProdMassTolXicPlus;     /// Tolerance of XicPlus mass from PDG value
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

  Double_t fKFPKs0_Chi2geoMax;  /// chi2/ndf(geo) cut of Ks0 reconstruction from KFParticle
  Double_t fKFPKs0_lDeltalMin;  /// l/Deltal cut of Ks0 reconstruction from KFParticle
  Double_t fKFPKs0_Chi2topoMax; /// chi2/ndf(topo) of Ks0 reconstruction from KFParticle
  Double_t fKFPLc_Chi2geoMax;   /// chi2/ndf(geo) cut of Lc reconstruction from KFParticle
  Double_t fKFPLam_Chi2geoMax;  /// chi2/ndf(geo) cut of lambda reconstruction from KFParticle
  Double_t fKFPLam_Chi2topoMin; /// chi2/ndf(topo) cut of lambda reconstruction from KFParticle
  Double_t fKFPLam_lDeltalMin;  /// l/Deltal cut of lambda reconstruction from KFParticle
  Double_t fKFPXi_Chi2geoMax;   /// chi2/ndf(geo) cut of Xi- reconstruction from KFParticle
  Double_t fKFPXi_Chi2topoMax;  /// chi2/ndf(topo) cut of Xi- reconstruction from KFParticle
  Double_t fKFPXi_lDeltalMin;   /// l/Deltal cut of lambda reconstruction from KFParticle
  Double_t fKFPXic0_Chi2geoMax; /// chi2/ndf(geo) cut of Xic0 reconstruction from KFParticle
  Double_t fKFPXicPlus_Chi2geoMax; /// chi2/ndf(geo) cut of XicPlus reconstruction from KFParticle
  Double_t fKFPXiccPP_Chi2geoMax; /// chi2/ndf(geo) cut of Xicc++ reconstruction from KFParticle
  Double_t fCut_nSigmaTPC_PiFromXic0_or_PiFromOmegac0; /// nSigmaTPC_PiFromXic0 or nSigmaTPC_PiFromOmegac0 cut
  Double_t fCut_nSigmaTPC_PiFromXi_or_KaFromOmega; // nSigmaTPC_PiFromXi or nSigmaTPC_KaFromOmega cut
  Double_t fCut_nSigmaTPC_PrFromLam; /// nSigmaTPC_PrFromLam cut
  Double_t fCut_nSigmaTPC_PiFromLam; /// nSigmaTPC_PiFromLam cut
  Double_t fCut_nSigmaTOF_PiFromXic0_or_PiFromOmegac0; /// nSigmaTOF_PiFromXic0 or nSigmaTOF_PiFromOmegac0 cut
  Double_t fCut_nSigmaTOF_PiFromXi_or_KaFromOmega; /// nSigmaTOF_PiFromXi or nSigmaTOF_KaFromOmega cut
  Double_t fCut_nSigmaTOF_PrFromLam; /// nSigmaTOF_PrFromLam cut
  Double_t fCut_nSigmaTOF_PiFromLam; /// nSigmaTOF_PiFromLam cut
  Double_t fCut_Armenteros; /// Armenteros-Podolanski (qT/|alpha|) cut

  TF1 *fWeight; /// Weight
  TF1 *fWeight_up; /// Weight_up
  TF1 *fWeight_dw; /// Weight_dw

  /// \cond CLASSIMP
  ClassDef(AliRDHFCutsKFP, 11);
  /// \endcond
};

#endif
