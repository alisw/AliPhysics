#ifndef ALIRDHFCUTSXICZEROTOXIPIFROMAODTRACKS_H
#define ALIRDHFCUTSXICZEROTOXIPIFROMAODTRACKS_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//***********************************************************
/// \class Class AliRDHFCutsXicZerotoXiPifromAODtracks
/// \brief class for cuts on AOD reconstructed Xic-> Xi pi
//  \Modified by Jianhui Zhu, zjh@mail.ccnu.edu.cn
//***********************************************************

#include "AliRDHFCuts.h"

class AliRDHFCutsXicZerotoXiPifromAODtracks : public AliRDHFCuts
{
 public:

  enum EPIDStrategy{
    kNSigmaCuts,
    kCombinedCuts
  };

  AliRDHFCutsXicZerotoXiPifromAODtracks(const char* name="CutsXicZerotoXiPi");
  virtual ~AliRDHFCutsXicZerotoXiPifromAODtracks();
  AliRDHFCutsXicZerotoXiPifromAODtracks(const AliRDHFCutsXicZerotoXiPifromAODtracks& source);
  AliRDHFCutsXicZerotoXiPifromAODtracks& operator=(const AliRDHFCutsXicZerotoXiPifromAODtracks& source);

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
  void SetUseCascadePID(Bool_t a){fUseCascadePID=a;}
  Double_t GetCombinedPIDThreshold(){return fCombinedPIDThreshold;}
  void SetPidCascPi(AliAODPidHF* pidPion) { 
      if(fPidObjCascPi) delete fPidObjCascPi;
      fPidObjCascPi=new AliAODPidHF(*pidPion);
      }
  AliAODPidHF* GetPidCascPi() const {return fPidObjCascPi;}
  void SetPidCascPr(AliAODPidHF* pidProton) { 
      if(fPidObjCascPr) delete fPidObjCascPr;
      fPidObjCascPr=new AliAODPidHF(*pidProton);
      }
  AliAODPidHF* GetPidCascPr() const {return fPidObjCascPr;}
  void SetPidCascKa(AliAODPidHF* pidKaon) { 
      if(fPidObjCascKa) delete fPidObjCascKa;
      fPidObjCascKa=new AliAODPidHF(*pidKaon);
      }
  AliAODPidHF* GetPidCascKa() const {return fPidObjCascKa;}


  Bool_t SingleTrkCuts(AliAODTrack *trk);
  Bool_t SingleCascadeCuts(AliAODcascade *casc, Double_t *vert, Bool_t anaOmegacZero);
  Bool_t SingleCascadeCutsRef(AliAODcascade *casc, Double_t *vert, Bool_t anaOmegaZero);
  Bool_t SelectWithRoughCuts(AliAODcascade *casc, AliAODTrack *trk1);

  void SetProdTrackPtMin(Double_t a){fProdTrackPtMin=a;}
  void SetProdTrackEtaRange(Double_t a){fProdTrackEtaRange=a;}
  void SetProdUseAODFilterBit(Bool_t a){fProdUseAODFilterBit=a;}
  void SetProdMassTolLambda(Double_t a){fProdMassTolLambda=a;}
  void SetProdMassTolXi(Double_t a){fProdMassTolXi=a;}
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
  void SetProdCascNTPCCrossedRowsMin(Double_t a){fProdCascNTPCCrossedRowsMin=a;}
  void SetProdCascNTPCCrossedOverFindalbleRatioMin(Double_t a){fProdCascNTPCCrossedOverFindableRatioMin=a;}

  Double_t GetProdTrackPtMin(){return fProdTrackPtMin;}
  Double_t GetProdTrackEtaRange(){return fProdTrackEtaRange;}
  Bool_t   GetProdUseAODFilterBit(){return fProdUseAODFilterBit;}
  Double_t GetProdMassTolLambda(){return fProdMassTolLambda;}
  Double_t GetProdMassTolXi(){return fProdMassTolXi;}
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
  Double_t GetProdCascNTPCCrossedRowsMin(){return fProdCascNTPCCrossedRowsMin;}
  Double_t GetProdCascNTPCCrossedOverFindalbleRatioMin(){return fProdCascNTPCCrossedOverFindableRatioMin;}
  
  void useSetNPtBins(Int_t nptBins){SetNPtBins(nptBins);}
 protected:
	
 private:

  EPIDStrategy fPIDStrategy;        /// PID Strategy
  Double_t fCombinedPIDThreshold;   /// PID threshold used in IsSelectedCombinedPID
  Bool_t fUseCascadePID;            /// Use PID for cascade or not
  AliAODPidHF *fPidObjCascPi;         /// PID object for cascade-pion
  AliAODPidHF *fPidObjCascPr;         /// PID object for cascade-proton
  AliAODPidHF *fPidObjCascKa;         /// PID object for cascade-kaon

  Double_t fProdTrackPtMin;         /// Minimum Bachelor pT
  Double_t fProdTrackEtaRange;      /// Bachelor Eta range
  Bool_t   fProdUseAODFilterBit;    /// Use AODfilterBit or not
  Double_t fProdMassTolLambda;      /// Tolerance of Lambda mass from PDG value
  Double_t fProdMassTolXi;          /// Tolerance of Xi mass from PDG value
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
  Double_t fProdCascNTPCCrossedRowsMin;   
  Double_t fProdCascNTPCCrossedOverFindableRatioMin;
  Double_t fProdLikeSignDcaMax;     /// Maximum DCA of pions
  Double_t fProdRoughMassTol;       /// Tolerance of Xic mass from PDG value
  Double_t fProdRoughPtMin;         /// Minimum pT of Xic

  /// \cond CLASSIMP
  ClassDef(AliRDHFCutsXicZerotoXiPifromAODtracks, 3); 
  /// \endcond
};

#endif
