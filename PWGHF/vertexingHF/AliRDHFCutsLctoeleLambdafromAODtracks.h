#ifndef ALIRDHFCUTSLCTOELELAMBDAFROMAODTRACKS_H
#define ALIRDHFCUTSLCTOELELAMBDAFROMAODTRACKS_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//***********************************************************
/// \class Class AliRDHFCutsLctoeleLambdafromAODtracks
/// \briefclass for cuts on AOD reconstructed Lc-> ele + Lambda
//***********************************************************

#include "AliRDHFCuts.h"

class AliRDHFCutsLctoeleLambdafromAODtracks : public AliRDHFCuts
{
 public:

  enum EPIDStrategy{
    kNSigmaCuts,
    kNSigmaCustomizedCuts,
    kNSigmaCustomizedPtDepCuts,
    kCombinedCuts
  };

  AliRDHFCutsLctoeleLambdafromAODtracks(const char* name="CutsLctoeleLambda");
  virtual ~AliRDHFCutsLctoeleLambdafromAODtracks();
  AliRDHFCutsLctoeleLambdafromAODtracks(const AliRDHFCutsLctoeleLambdafromAODtracks& source);
  AliRDHFCutsLctoeleLambdafromAODtracks& operator=(const AliRDHFCutsLctoeleLambdafromAODtracks& source);

  using AliRDHFCuts::GetCutVarsForOpt;
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters);

  using AliRDHFCuts::IsSelected;
  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel);
  using AliRDHFCuts::IsSelectedPID;
  virtual Int_t IsSelectedPID(AliAODRecoDecayHF* obj);
  Int_t IsSelectedCombinedPID(AliAODRecoDecayHF* obj);
  Bool_t IsSelectedeID(AliAODTrack* trk);
  Bool_t IsSelectedCustomizedeID(AliAODTrack* trk);
  Bool_t IsSelectedCustomizedPtDepeID(AliAODTrack* trk);
  Bool_t IsSelectedCombinedeID(AliAODTrack* trk);

  void SetPIDStrategy(EPIDStrategy pidStrategy){fPIDStrategy=pidStrategy;}
  EPIDStrategy GetPIDStrategy() const {return fPIDStrategy;}
  void SetCombinedPIDThreshold(Double_t a){fCombinedPIDThreshold=a;}
  Double_t GetCombinedPIDThreshold(){return fCombinedPIDThreshold;}

  void SetUseOnTheFlyV0(Bool_t a) { fUseOnTheFlyV0=a; }
  Bool_t GetUseOnTheFlyV0() { return fUseOnTheFlyV0; }

  Bool_t SingleTrkCuts(AliAODTrack *trk, AliAODVertex *vert);
  Bool_t SingleTrkCutsNoPID(AliAODTrack *trk, AliAODVertex *vert);
  Bool_t SingleV0Cuts(AliAODv0 *v0, AliAODVertex *vert);
  Bool_t SelectWithRoughCuts(AliAODv0 *v0, AliAODTrack *trk1);

  void SetProdTrackTPCNclsPIDMin(Int_t a){fProdTrackTPCNclsPIDMin=a;}
  void SetProdTrackTPCNclsRatioMin(Double_t a){fProdTrackTPCNclsRatioMin=a;}
  void SetProdUseAODFilterBit(Bool_t a){fProdUseAODFilterBit=a;}
  void SetProdV0MassTolLambda(Double_t a){fProdV0MassTolLambda=a;}
  void SetProdV0PtMin(Double_t a){fProdV0PtMin=a;}
  void SetProdV0CosPointingAngleToPrimVtxMin(Double_t a){fProdV0CosPointingAngleToPrimVtxMin=a;}
  void SetProdV0DcaDaughtersMax(Double_t a){fProdV0DcaDaughtersMax=a;}
  void SetProdV0DaughterEtaRange(Double_t a){fProdV0DaughterEtaRange=a;}
  void SetProdV0DaughterPtMin(Double_t a){fProdV0DaughterPtMin=a;}
  void SetProdV0DaughterTPCClusterMin(Double_t a){fProdV0DaughterTPCClusterMin=a;}
  void SetProdV0DaughterTPCCrossRatioMin(Double_t a){fProdV0DaughterTPCCrossRatioMin=a;}
  void SetProdRfidMinV0(Double_t a){fProdRfidMinV0=a;}
  void SetProdRfidMaxV0(Double_t a){fProdRfidMaxV0=a;}
  void SetProdDcaV0ToPrimVertexMin(Double_t a){fProdDcaV0ToPrimVertexMin=a;}
  void SetProdDcaV0PrToPrimVertexMin(Double_t a){fProdDcaV0PrToPrimVertexMin=a;}
  void SetProdDcaV0PiToPrimVertexMin(Double_t a){fProdDcaV0PiToPrimVertexMin=a;}
  void SetProdV0ProperDecayLengthMax(Double_t a){fProdV0ProperDecayLengthMax=a;}
  void SetProdMassRejK0s(Double_t a){fProdMassRejK0s=a;}
  void SetProdV0EtaRange(Double_t a, Double_t b){fProdV0EtaMin=a;fProdV0EtaMax=b;}
  void SetProdV0RapRange(Double_t a, Double_t b){fProdV0RapMin=a;fProdV0RapMax=b;}
  
  void SetProdRoughMassTol(Double_t a){fProdRoughMassTol=a;}
  void SetProdRoughPtMin(Double_t a){fProdRoughPtMin=a;}
  
  Int_t GetProdTrackTPCNclsPIDMin(){return fProdTrackTPCNclsPIDMin;}
  Double_t GetProdTrackTPCNclsRatioMin(){return fProdTrackTPCNclsRatioMin;}
  Bool_t   GetProdUseAODFilterBit(){return fProdUseAODFilterBit;}
  Double_t GetProdV0MassTolLambda(){return fProdV0MassTolLambda;}
  Double_t GetProdV0PtMin(){return fProdV0PtMin;}
  Double_t GetProdV0CosPointingAngleToPrimVtxMin(){return fProdV0CosPointingAngleToPrimVtxMin;}
  Double_t GetProdV0DcaDaughtersMax(){return fProdV0DcaDaughtersMax;}
  Double_t GetProdV0DaughterEtaRange(){return fProdV0DaughterEtaRange;}
  Double_t GetProdV0DaughterPtMin(){return fProdV0DaughterPtMin;}
  Double_t GetProdV0DaughterTPCClusterMin(){return fProdV0DaughterTPCClusterMin;}
  Double_t GetProdV0DaughterTPCCrossRatioMin(){return fProdV0DaughterTPCCrossRatioMin;}
  Double_t GetProdRfidMinV0(){return fProdRfidMinV0;}
  Double_t GetProdRfidMaxV0(){return fProdRfidMaxV0;}
  Double_t GetProdDcaV0ToPrimVertexMin(){return fProdDcaV0ToPrimVertexMin;}
  Double_t GetProdDcaV0PrToPrimVertexMin(){return fProdDcaV0PrToPrimVertexMin;}
  Double_t GetProdDcaV0PiToPrimVertexMin(){return fProdDcaV0PiToPrimVertexMin;}
  Double_t GetProdV0ProperDecayLengthMax(){return fProdV0ProperDecayLengthMax;}
  Double_t GetProdMassRejK0s(){return fProdMassRejK0s;}
  void GetProdV0EtaRange(Double_t &a, Double_t &b){a=fProdV0EtaMin;b=fProdV0EtaMax;}
  void GetProdV0RapRange(Double_t &a, Double_t &b){a=fProdV0RapMin;b=fProdV0RapMax;}

  Double_t GetProdRoughMassTol(){return fProdRoughMassTol;}
  Double_t GetProdRoughPtMin(){return fProdRoughPtMin;}

  void SetUseLambdaPID(Bool_t a){fUseLambdaPID=a;}
  Bool_t GetUseLambdaPID(){return fUseLambdaPID;}
  void SetPidProton(AliAODPidHF* pidProton) { 
      if(fPidObjProton) delete fPidObjProton;
      fPidObjProton=new AliAODPidHF(*pidProton);
      }
  AliAODPidHF* GetPidProton() const {return fPidObjProton;}
  void SetPidPion(AliAODPidHF* pidPion) { 
      if(fPidObjPion) delete fPidObjPion;
      fPidObjPion=new AliAODPidHF(*pidPion);
      }
  AliAODPidHF* GetPidPion() const {return fPidObjPion;}
	void GetSigmaElectronTPCRange(Double_t &a,Double_t &b){a=fSigmaElectronTPCMin;b=fSigmaElectronTPCMax;}
	void GetSigmaElectronTOFRange(Double_t &a,Double_t &b){a=fSigmaElectronTOFMin;b=fSigmaElectronTOFMax;}
	void GetSigmaElectronTPCPtDepPars(Double_t &a,Double_t &b){a=fSigmaElectronTPCPtDepPar0;b=fSigmaElectronTPCPtDepPar1;}

	void SetExcludePionTPC(Bool_t a){fExcludePionTPC=a;}
	void SetExcludeProtonTPC(Bool_t a){fExcludeProtonTPC=a;}
	void SetExcludeKaonTPC(Bool_t a){fExcludeKaonTPC=a;}
	void SetExcludenSigmaPionTPC(Double_t a){fExcludenSigmaPionTPC=a;}
	void SetExcludenSigmaProtonTPC(Double_t a){fExcludenSigmaProtonTPC=a;}
	void SetExcludenSigmaKaonTPC(Double_t a){fExcludenSigmaKaonTPC=a;}
	void SetSigmaElectronTPCRange(Double_t a,Double_t b){fSigmaElectronTPCMin=a;fSigmaElectronTPCMax=b;}
	void SetSigmaElectronTOFRange(Double_t a,Double_t b){fSigmaElectronTOFMin=a;fSigmaElectronTOFMax=b;}
	void SetSigmaElectronTPCPtDepPars(Double_t a,Double_t b){fSigmaElectronTPCPtDepPar0=a;fSigmaElectronTPCPtDepPar1=b;}

 protected:
	
 private:

  EPIDStrategy fPIDStrategy;        /// PID strategy
  Double_t fCombinedPIDThreshold;   /// Threshold used in  IsSelectedCombinedPID
  Bool_t fUseLambdaPID;            /// Use PID for proton from Lc
  AliAODPidHF *fPidObjProton;         /// PID object for proton from Lc
  AliAODPidHF *fPidObjPion;         /// PID object for proton from Lc
  Bool_t   fUseOnTheFlyV0;          /// Flag to check if we use on-the-fly v0
  
  Int_t fProdTrackTPCNclsPIDMin;      /// Min. Number of TPC PID cluster
  Double_t fProdTrackTPCNclsRatioMin;      /// Min. Number of TPC PID cluster
  Bool_t   fProdUseAODFilterBit;    /// Flag for AOD filter Bit used before object creation
  Double_t fProdV0MassTolLambda;       /// Lambda mass selection  used before object creation
  Double_t fProdV0PtMin;            /// Minimum Lambda pT used before object creation
  Double_t fProdV0CosPointingAngleToPrimVtxMin;/// V0 pointing angle used before object creation
  Double_t fProdV0DcaDaughtersMax;  /// Max DCA between V0 daughters used before object creation
  Double_t fProdV0DaughterEtaRange; /// V0Daughter eta range used before object creation
  Double_t fProdV0DaughterPtMin;    /// V0 Daughter pT min used before object creation
  Double_t fProdV0DaughterTPCClusterMin;/// V0 daughter Minimum TPC cluster pT used before object creation
  Double_t fProdV0DaughterTPCCrossRatioMin;/// V0 daughter Minimum TPC cluster pT used before object creation
  Double_t fProdRfidMinV0;          /// Minimum Decay vertex of V0
  Double_t fProdRfidMaxV0;          /// Max Decay vertex of V0
  Double_t fProdDcaV0ToPrimVertexMin;  /// Min Dca between v0 and PV
  Double_t fProdDcaV0PrToPrimVertexMin;  /// Min Dca between v0-proton and PV
  Double_t fProdDcaV0PiToPrimVertexMin;  /// Min Dca between v0-pion and PV
  Double_t fProdV0ProperDecayLengthMax;        /// mL/p of cascade
  Double_t fProdMassRejK0s;          /// Rejection range of Omega mass from PDG value
	Double_t fProdV0EtaMin; /// Minimum eta of cascade
	Double_t fProdV0EtaMax; /// Maximum eta of cascade
	Double_t fProdV0RapMin; /// Minimum rapidity of cascade
	Double_t fProdV0RapMax; /// Maximum rapidity of cascade
  Double_t fProdRoughMassTol;       /// Mass cut for Lc used before object creation
  Double_t fProdRoughPtMin;         /// pT cut for Lc used before object creation

	Bool_t fExcludePionTPC;    /// Flag wheter to exlude pion band
	Bool_t fExcludeProtonTPC;  /// Flag wheter to exlude proton band
	Bool_t fExcludeKaonTPC;    /// Flag wheter to exlude proton band
	Double_t fExcludenSigmaPionTPC; /// nSigma to exclude for pion band
	Double_t fExcludenSigmaProtonTPC; /// nSigma to exclude for proton band
	Double_t fExcludenSigmaKaonTPC; /// nSigma to exclude for Kaon band
	Double_t fSigmaElectronTPCMin;  /// nSigma to exclude for Kaon band
	Double_t fSigmaElectronTPCPtDepPar0; /// nSigma electron lower limit (par0)
	Double_t fSigmaElectronTPCPtDepPar1; /// nSigma electron lower limit (par1)
	Double_t fSigmaElectronTPCMax; /// nSigma to exclude for Kaon band
	Double_t fSigmaElectronTOFMin; /// nSigma to exclude for Kaon band
	Double_t fSigmaElectronTOFMax; /// nSigma to exclude for Kaon band

  /// \cond CLASSIMP     
  ClassDef(AliRDHFCutsLctoeleLambdafromAODtracks,3);
  /// \endcond
};

#endif
