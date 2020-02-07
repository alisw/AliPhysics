#ifndef ALIRDHFCUTSXICTOELEXIFROMAODTRACKS_H
#define ALIRDHFCUTSXICTOELEXIFROMAODTRACKS_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//***********************************************************
/// \class Class AliRDHFCutsXictoeleXifromAODtracks
/// \brief class for cuts on AOD reconstructed Xic-> ele + Xi
//***********************************************************

#include "AliRDHFCuts.h"

class AliRDHFCutsXictoeleXifromAODtracks : public AliRDHFCuts
{
 public:

  enum EPIDStrategy{
    kNSigmaCuts,
    kNSigmaCustomizedCuts,
    kNSigmaCustomizedPtDepCuts,
    kCombinedCuts
  };

  AliRDHFCutsXictoeleXifromAODtracks(const char* name="CutsXictoeleXi");
  virtual ~AliRDHFCutsXictoeleXifromAODtracks();
  AliRDHFCutsXictoeleXifromAODtracks(const AliRDHFCutsXictoeleXifromAODtracks& source);
  AliRDHFCutsXictoeleXifromAODtracks& operator=(const AliRDHFCutsXictoeleXifromAODtracks& source);

  using AliRDHFCuts::GetCutVarsForOpt;
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters);

  using AliRDHFCuts::IsSelected;
  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel);
  using AliRDHFCuts::IsSelectedPID;
  virtual Int_t IsSelectedPID(AliAODRecoDecayHF* obj);
  Int_t IsSelectedCombinedPID(AliAODRecoDecayHF* obj);
  Bool_t IsSelectedeID(AliAODTrack* trk);
  Bool_t IsSelectedCustomizedeID(AliAODTrack* trk);
  Bool_t IsSelectedCustomizedPtDepeID(AliAODTrack* trk, AliAODTrack *trkpid);
  Bool_t IsSelectedCombinedeID(AliAODTrack* trk);
  Int_t IsSelected(TLorentzVector* vtrk, TLorentzVector *vcasc, Double_t *cutvars, Int_t selectionLevel);

  void SetPIDStrategy(EPIDStrategy pidStrategy){fPIDStrategy=pidStrategy;}
  EPIDStrategy GetPIDStrategy() const {return fPIDStrategy;}
  void SetCombinedPIDThreshold(Double_t a){fCombinedPIDThreshold=a;}
  Double_t GetCombinedPIDThreshold(){return fCombinedPIDThreshold;}

  void SetUseV0Topology(Int_t a) { fUseV0Topology=a; }
  Int_t GetUseV0Topology() { return fUseV0Topology; }

  Bool_t SingleTrkCuts(AliAODTrack *trk, AliAODTrack *trkpid, AliAODVertex *primvert);
  Bool_t SingleTrkCutsNoPID(AliAODTrack *trk, AliAODTrack *trkpid, AliAODVertex *primvert);
  Bool_t SingleCascadeCuts(AliAODcascade *casc, Double_t *vert);
	Bool_t TagConversions(AliAODTrack *etrk, Int_t *id2index, AliAODEvent *evt, Int_t ntrk, Double_t &minmass);
	Bool_t TagConversionsSameSign(AliAODTrack *etrk, Int_t *id2index, AliAODEvent *evt, Int_t ntrk, Double_t &minmass);
  Bool_t SelectWithRoughCuts(AliAODcascade *casc, AliAODTrack *trk1);

  void SetMagneticField(Double_t a){fBzkG = a;}
  void SetPrimaryVertex(Double_t *a){fPrimVert[0] = a[0];fPrimVert[1] = a[1];fPrimVert[2] = a[2];}

  void SetProdTrackTPCNclsPIDMin(Int_t a){fProdTrackTPCNclsPIDMin=a;}
  void SetProdTrackTPCNclsRatioMin(Double_t a){fProdTrackTPCNclsRatioMin=a;}
  void SetProdUseAODFilterBit(Bool_t a){fProdUseAODFilterBit=a;}
  void SetProdAODFilterBit(Int_t a){fProdAODFilterBit=a;}
  void SetProdRejectTrackWithShared(Bool_t a){fProdRejectTrackWithShared=a;}
  void SetProdMassTolLambda(Double_t a){fProdMassTolLambda=a;}
  void SetProdMassTolXiRough(Double_t a){fProdMassTolXiRough=a;}
  void SetProdMassTolXi(Double_t a){fProdMassTolXi=a;}
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
  void SetProdCascNTPCClustersMin(Double_t a){fProdCascNTPCClustersMin=a;}
  void SetProdCascCutMinNCrossedRowsTPC(Double_t a){fProdCascCutMinNCrossedRowsTPC=a;}
  void SetProdCascratioCrossedRowsOverFindableClusterTPC(Double_t a){fProdCascratioCrossedRowsOverFindableClusterTPC=a;}
  void SetProdCascEtaRange(Double_t a, Double_t b){fProdCascEtaMin=a;fProdCascEtaMax=b;}
  void SetProdCascRapRange(Double_t a, Double_t b){fProdCascRapMin=a;fProdCascRapMax=b;}
  void SetProdCascPtRange(Double_t a, Double_t b){fProdCascPtMin=a;fProdCascPtMax=b;}
  void SetProdRoughMassTol(Double_t a){fProdRoughMassTol=a;}
  void SetProdRoughPtMin(Double_t a){fProdRoughPtMin=a;}
  
  Int_t GetProdTrackTPCNclsPIDMin(){return fProdTrackTPCNclsPIDMin;}
  Double_t GetProdTrackTPCNclsRatioMin(){return fProdTrackTPCNclsRatioMin;}
  Bool_t   GetProdUseAODFilterBit(){return fProdUseAODFilterBit;}
  Int_t   GetProdAODFilterBit(){return fProdAODFilterBit;}
  Bool_t   GetProdRejectTrackWithShared(){return fProdRejectTrackWithShared;}
  Double_t GetProdMassTolLambda(){return fProdMassTolLambda;}
  Double_t GetProdMassTolXiRough(){return fProdMassTolXiRough;}
  Double_t GetProdMassTolXi(){return fProdMassTolXi;}
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
  Double_t GetProdCascNTPCClustersMin(){return fProdCascNTPCClustersMin;}
  Double_t GetProdCascCutMinNCrossedRowsTPC(){return fProdCascCutMinNCrossedRowsTPC;}
  Double_t GetProdCascratioCrossedRowsOverFindableClusterTPC(){return fProdCascratioCrossedRowsOverFindableClusterTPC;}
  void GetProdCascEtaRange(Double_t &a, Double_t &b){a=fProdCascEtaMin;b=fProdCascEtaMax;}
  void GetProdCascRapRange(Double_t &a, Double_t &b){a=fProdCascRapMin;b=fProdCascRapMax;}
  void GetProdCascPtRange(Double_t &a, Double_t &b){a=fProdCascPtMin;b=fProdCascPtMax;}

  Bool_t GetUseCascadePID(){return fUseCascadePID;}
  void SetUseCascadePID(Bool_t a){fUseCascadePID=a;}
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
	void GetSigmaElectronTPCRange(Double_t &a,Double_t &b){a=fSigmaElectronTPCMin;b=fSigmaElectronTPCMax;}
	void GetSigmaElectronTOFRange(Double_t &a,Double_t &b){a=fSigmaElectronTOFMin;b=fSigmaElectronTOFMax;}
	void GetSigmaElectronTPCPtDepPars(Double_t &a,Double_t &b){a=fSigmaElectronTPCPtDepPar0;b=fSigmaElectronTPCPtDepPar1;}
	void GetSigmaElectronTPCPtDepPars(Double_t &a,Double_t &b,Double_t &c){a=fSigmaElectronTPCPtDepPar0;b=fSigmaElectronTPCPtDepPar1;c=fSigmaElectronTPCPtDepPar2;}
	Double_t GetConversionMassMax(){return fConversionMassMax;}
	Double_t GetEleXiMassMax(){return fEleXiMassMax;}

	void SetExcludePionTPC(Bool_t a){fExcludePionTPC=a;}
	void SetExcludeProtonTPC(Bool_t a){fExcludeProtonTPC=a;}
	void SetExcludeKaonTPC(Bool_t a){fExcludeKaonTPC=a;}
	void SetExcludenSigmaPionTPC(Double_t a){fExcludenSigmaPionTPC=a;}
	void SetExcludenSigmaProtonTPC(Double_t a){fExcludenSigmaProtonTPC=a;}
	void SetExcludenSigmaKaonTPC(Double_t a){fExcludenSigmaKaonTPC=a;}
	void SetSigmaElectronTPCRange(Double_t a,Double_t b){fSigmaElectronTPCMin=a;fSigmaElectronTPCMax=b;}
	void SetSigmaElectronTOFRange(Double_t a,Double_t b){fSigmaElectronTOFMin=a;fSigmaElectronTOFMax=b;}
	void SetSigmaElectronTPCPtDepPars(Double_t a,Double_t b){fSigmaElectronTPCPtDepPar0=a;fSigmaElectronTPCPtDepPar1=b;}
	void SetSigmaElectronTPCPtDepPars(Double_t a,Double_t b,Double_t c){fSigmaElectronTPCPtDepPar0=a;fSigmaElectronTPCPtDepPar1=b;fSigmaElectronTPCPtDepPar2=c;}
	void SetConversionMassMax(Double_t a){fConversionMassMax=a;}
	void SetEleXiMassMax(Double_t a){fEleXiMassMax=a;}

	Bool_t IsPeakRegion(AliAODcascade *c);
	Bool_t IsPeakRegion(TLorentzVector *c);
	Bool_t IsSideBand(AliAODcascade *c);
	Bool_t IsSideBand(TLorentzVector *c);

  void SetSftPosR125(AliAODTrack *track,Double_t bfield,Double_t priVtx[3], Double_t *XSftR125);
  Double_t dEtaSR125(Double_t *postrack1,Double_t *postrack2);
  Double_t dPhiSR125(Double_t *postrack1,Double_t *postrack2);
  Double_t GetdPhiSdEtaSR125(AliAODTrack *tracke, AliAODTrack *trackp,AliAODTrack *trackn, AliAODTrack *trackb, Double_t bfield,Double_t priVtx[3], Double_t &dPhiS_ep, Double_t &dEtaS_ep,Double_t &dPhiS_en, Double_t &dEtaS_en, Double_t &dPhiS_eb, Double_t &dEtaS_eb);
  Double_t DeltaPhi(AliAODcascade *casc, AliAODTrack *trk);
  Double_t CosOpeningAngle(AliAODcascade *casc, AliAODTrack *trk);
  Double_t DeltaEta(AliAODcascade *casc, AliAODTrack *trk);

 protected:
	
 private:

  EPIDStrategy fPIDStrategy;        /// PID strategy
  Double_t fCombinedPIDThreshold;   /// Threshold used in  IsSelectedCombinedPID
  Bool_t fUseCascadePID;            /// Use PID for cascade or not
  AliAODPidHF *fPidObjCascPi;         /// PID object for cascade-pion
  AliAODPidHF *fPidObjCascPr;         /// PID object for cascade-proton
//  Bool_t   fUseOnTheFlyV0;          /// Flag to check if we use on-the-fly v0
  Int_t   fUseV0Topology;          /// 0: Cowboy+Sailor 1: Cowboy 2:Sailor
  Double_t fBzkG; ///B field
  Double_t fPrimVert[3];///Primary vertex
  
  Int_t fProdTrackTPCNclsPIDMin;      /// Min. Number of TPC PID cluster
  Double_t fProdTrackTPCNclsRatioMin;      /// Min. Number of TPC PID cluster
  Bool_t   fProdUseAODFilterBit;    /// Flag for AOD filter Bit used before object creation
  Int_t    fProdAODFilterBit;    /// AOD filter Bit used before object creation
  Bool_t   fProdRejectTrackWithShared;    /// Flag to Reject tracks with shared clusters
  Double_t fProdMassTolLambda;      /// Tolerance of Lambda mass from PDG value
  Double_t fProdMassTolXiRough;          /// Tolerance of Xi mass from PDG value (including sideband)
  Double_t fProdMassTolXi;          /// Tolerance of Xi mass from PDG value
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
  Double_t fProdV0CosineOfPoiningAngleXiMin;  /// Min V0 cos pointing angle  to Xi vertex
  Double_t fProdCascNTPCClustersMin;         /// Minimum number of TPC clusters
  Double_t fProdCascCutMinNCrossedRowsTPC;  // Minimum number of Crossed Rows TPC
  Double_t fProdCascratioCrossedRowsOverFindableClusterTPC; // Minmum ratio of CrossedRows Over Findable Cluster
  Double_t fProdCascEtaMin; /// Minimum eta of cascade
  Double_t fProdCascEtaMax; /// Maximum eta of cascade
  Double_t fProdCascRapMin; /// Minimum rapidity of cascade
  Double_t fProdCascRapMax; /// Maximum rapidity of cascade
  Double_t fProdCascPtMin; /// Minimum Pt of cascade
  Double_t fProdCascPtMax; /// Maximum Pt of cascade
  Double_t fProdRoughMassTol;       /// Mass cut for Lc used before object creation
  Double_t fProdRoughPtMin;         /// pT cut for Lc used before object creation

	Bool_t fExcludePionTPC; /// Flag wheter to exlude pion band
	Bool_t fExcludeProtonTPC; /// Flag wheter to exlude proton band
	Bool_t fExcludeKaonTPC; /// Flag wheter to exlude proton band
	Double_t fExcludenSigmaPionTPC; /// nSigma to exclude for pion band
	Double_t fExcludenSigmaProtonTPC; /// nSigma to exclude for proton band
	Double_t fExcludenSigmaKaonTPC; /// nSigma to exclude for Kaon band
	Double_t fSigmaElectronTPCMin; /// nSigma to exclude for Kaon band
	Double_t fSigmaElectronTPCPtDepPar0; /// nSigma electron lower limit (par0)
	Double_t fSigmaElectronTPCPtDepPar1; /// nSigma electron lower limit (par1)
	Double_t fSigmaElectronTPCPtDepPar2; /// nSigma electron lower limit (par2)
	Double_t fSigmaElectronTPCMax; /// nSigma to exclude for Kaon band
	Double_t fSigmaElectronTOFMin; /// nSigma to exclude for Kaon band
	Double_t fSigmaElectronTOFMax; /// nSigma to exclude for Kaon band

	Double_t fConversionMassMax; /// Conversion mass
	Double_t fEleXiMassMax; /// e-xi mass max
  
  /// \cond CLASSIMP
  ClassDef(AliRDHFCutsXictoeleXifromAODtracks,9);
  /// \endcond
};

#endif
