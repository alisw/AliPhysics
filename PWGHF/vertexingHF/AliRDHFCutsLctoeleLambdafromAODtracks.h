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
  Bool_t IsSelectedCustomizedPtDepeID(AliAODTrack* trk, AliAODTrack *trkpid);
  Bool_t IsSelectedCombinedeID(AliAODTrack* trk);
  Int_t IsSelected(TLorentzVector* vtrk, TLorentzVector *vv0, Double_t *cutvars, Int_t selectionLevel);

  void SetPIDStrategy(EPIDStrategy pidStrategy){fPIDStrategy=pidStrategy;}
  EPIDStrategy GetPIDStrategy() const {return fPIDStrategy;}
  void SetCombinedPIDThreshold(Double_t a){fCombinedPIDThreshold=a;}
  Double_t GetCombinedPIDThreshold(){return fCombinedPIDThreshold;}

  void SetUseOnTheFlyV0(Bool_t a) { fUseOnTheFlyV0=a; }
  Bool_t GetUseOnTheFlyV0() { return fUseOnTheFlyV0; }
  void SetUseV0Topology(Int_t a) { fUseV0Topology=a; }
  Int_t GetUseV0Topology() { return fUseV0Topology; }

  Bool_t SingleTrkCuts(AliAODTrack *trk, AliAODTrack *trkpid, AliAODVertex *vert);
  Bool_t SingleTrkCutsNoPID(AliAODTrack *trk, AliAODTrack *trkpid, AliAODVertex *vert);
  Bool_t SingleV0Cuts(AliAODv0 *v0, AliAODVertex *vert);
	Bool_t TagConversions(AliAODTrack *etrk, Int_t *id2index, AliAODEvent *evt, Int_t ntrk, Double_t &minmass);
	Bool_t TagConversionsSameSign(AliAODTrack *etrk, Int_t *id2index, AliAODEvent *evt, Int_t ntrk, Double_t &minmass);
  Bool_t SelectWithRoughCuts(AliAODv0 *v0, AliAODTrack *trk1);

  void SetMagneticField(Double_t a){fBzkG = a;}
  void SetPrimaryVertex(Double_t *a){fPrimVert[0] = a[0];fPrimVert[1] = a[1];fPrimVert[2] = a[2];}

  void SetProdTrackTPCNclsPIDMin(Int_t a){fProdTrackTPCNclsPIDMin=a;}
  void SetProdTrackTPCNclsRatioMin(Double_t a){fProdTrackTPCNclsRatioMin=a;}
  void SetProdUseAODFilterBit(Bool_t a){fProdUseAODFilterBit=a;}
  void SetProdAODFilterBit(Int_t a){fProdAODFilterBit=a;}
  void SetProdRejectTrackWithShared(Bool_t a){fProdRejectTrackWithShared=a;}
  void SetProdV0KinkRejection(Bool_t a){fProdV0KinkRejection=a;}
  void SetProdV0MassTolLambda(Double_t a){fProdV0MassTolLambda=a;}
  void SetProdV0MassTolLambdaRough(Double_t a){fProdV0MassTolLambdaRough=a;}
  void SetProdV0PtMin(Double_t a){fProdV0PtMin=a;}
  void SetProdV0PtMax(Double_t a){fProdV0PtMax=a;}
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
  Int_t   GetProdAODFilterBit(){return fProdAODFilterBit;}
  Bool_t   GetProdRejectTrackWithShared(){return fProdRejectTrackWithShared;}
  Bool_t   GetProdV0KinkRejection(){return fProdV0KinkRejection;}
  Double_t GetProdV0MassTolLambda(){return fProdV0MassTolLambda;}
  Double_t GetProdV0MassTolLambdaRough(){return fProdV0MassTolLambdaRough;}
  Double_t GetProdV0PtMin(){return fProdV0PtMin;}
  Double_t GetProdV0PtMax(){return fProdV0PtMax;}
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
	void GetSigmaElectronTPCPtDepPars(Double_t &a,Double_t &b,Double_t &c){a=fSigmaElectronTPCPtDepPar0;b=fSigmaElectronTPCPtDepPar1;c=fSigmaElectronTPCPtDepPar2;}
	Double_t GetConversionMassMax(){return fConversionMassMax;}

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
	Bool_t IsPeakRegion(AliAODv0 *c);
	Bool_t IsPeakRegion(TLorentzVector *c);
	Bool_t IsSideBand(AliAODv0 *c);
	Bool_t IsSideBand(TLorentzVector *c);
  void SetSftPosR125(AliAODTrack *track,Double_t bfield,Double_t priVtx[3], Double_t *XSftR125);
  void SetSftPosR(AliAODTrack *track,Double_t bfield,Double_t R, Double_t priVtx[3], Double_t *XSftR);
  Double_t dEtaSR125(Double_t *postrack1,Double_t *postrack2);
  Double_t dPhiSR125(Double_t *postrack1,Double_t *postrack2);
  Double_t GetdPhiSdEtaSR125(AliAODTrack *tracke, AliAODTrack *trackp,AliAODTrack *trackn, Double_t bfield,Double_t priVtx[3], Double_t &dPhiS_ep, Double_t &dEtaS_ep,Double_t &dPhiS_en, Double_t &dEtaS_en);
  Double_t CalculatePhotonMass(AliAODTrack *track1, AliAODTrack *track2);
  Double_t DeltaPhi(AliAODv0 *v0, AliAODTrack *trk);
  Double_t DeltaEta(AliAODv0 *v0, AliAODTrack *trk);

 protected:
	
 private:

  EPIDStrategy fPIDStrategy;        /// PID strategy
  Double_t fCombinedPIDThreshold;   /// Threshold used in  IsSelectedCombinedPID
  Bool_t fUseLambdaPID;            /// Use PID for proton from Lc
  AliAODPidHF *fPidObjProton;         /// PID object for proton from Lc
  AliAODPidHF *fPidObjPion;         /// PID object for proton from Lc
  Bool_t   fUseOnTheFlyV0;          /// Flag to check if we use on-the-fly v0
  Int_t   fUseV0Topology;          /// 0: Cowboy+Sailor 1: Cowboy 2:Sailor
  Double_t fBzkG; ///B field
  Double_t fPrimVert[3]; ///Primary vertex
  
  Int_t fProdTrackTPCNclsPIDMin;      /// Min. Number of TPC PID cluster
  Double_t fProdTrackTPCNclsRatioMin;      /// Min. Number of TPC PID cluster
  Bool_t   fProdUseAODFilterBit;    /// Flag for AOD filter Bit used before object creation
  Int_t    fProdAODFilterBit;    /// AOD filter Bit used before object creation
  Bool_t   fProdRejectTrackWithShared;    /// Flag to Reject tracks with shared clusters
  Bool_t   fProdV0KinkRejection;    /// Flag to Reject v0 kinks
  Double_t fProdV0MassTolLambda;       /// Lambda mass selection  used before object creation
  Double_t fProdV0MassTolLambdaRough;       /// Lambda mass selection  used before object creation
  Double_t fProdV0PtMin;            /// Minimum Lambda pT used before object creation
  Double_t fProdV0PtMax;            /// Max Lambda pT used before object creation
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
	Double_t fSigmaElectronTPCPtDepPar2; /// nSigma electron lower limit (par2)
	Double_t fSigmaElectronTPCMax; /// nSigma to exclude for Kaon band
	Double_t fSigmaElectronTOFMin; /// nSigma to exclude for Kaon band
	Double_t fSigmaElectronTOFMax; /// nSigma to exclude for Kaon band

	Double_t fConversionMassMax; /// Conversion mass

  /// \cond CLASSIMP     
  ClassDef(AliRDHFCutsLctoeleLambdafromAODtracks,9);
  /// \endcond
};

#endif
