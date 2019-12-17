#ifndef ALIANALYSISVERTEXINGHF_H
#define ALIANALYSISVERTEXINGHF_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//-------------------------------------------------------------------------
/// \class                     Class AliAnalysisVertexingHF
/// \brief           Reconstruction of heavy-flavour decay candidates
///
/// \author Origin: E.Bruna, G.E.Bruno, A.Dainese, F.Prino, R.Romita, X.M.Zhang
/// \author Contact: andrea.dainese@pd.infn.it
//-------------------------------------------------------------------------

#include <TNamed.h>
#include <TList.h>

#include "AliAnalysisFilter.h"
#include "AliESDtrackCuts.h"

class AliPIDResponse;
class AliESDVertex;
class AliAODRecoDecay;
class AliAODRecoDecayHF;
class AliAODRecoDecayHF2Prong;
class AliAODRecoDecayHF3Prong;
class AliAODRecoDecayHF4Prong;
class AliAODRecoCascadeHF;
class AliAnalysisFilter;
class AliRDHFCuts;
class AliRDHFCutsD0toKpi;
class AliRDHFCutsJpsitoee;
class AliRDHFCutsDplustoK0spi;
class AliRDHFCutsDplustoKpipi;
class AliRDHFCutsDstoK0sK;
class AliRDHFCutsDstoKKpi;
class AliRDHFCutsLctopKpi;
class AliRDHFCutsLctoV0;
class AliRDHFCutsD0toKpipipi;
class AliRDHFCutsDStartoKpipi;
class AliESDtrack;
class AliVEvent;
class AliAODVertex;
class AliVertexerTracks;
class AliESDv0;
class AliAODv0;

//-----------------------------------------------------------------------------
class AliAnalysisVertexingHF : public TNamed {
 public:
  //
  AliAnalysisVertexingHF();
  AliAnalysisVertexingHF(const AliAnalysisVertexingHF& source);
  AliAnalysisVertexingHF& operator=(const AliAnalysisVertexingHF& source);
  virtual ~AliAnalysisVertexingHF();

  void FindCandidates(AliVEvent *event,
		      TClonesArray *aodVerticesHFTClArr,
		      TClonesArray *aodD0toKpiTClArr,
		      TClonesArray *aodJPSItoEleTClArr,
		      TClonesArray *aodCharm3ProngTClArr,
		      TClonesArray *aodCharm4ProngTClArr,
		      TClonesArray *aodDstarTClArr,
		      TClonesArray *aodCascadesTClArr,
		      TClonesArray *aodLikeSign2ProngTClArr,
		      TClonesArray *aodLikeSign3ProngTClArr);

  TList* FillListOfCuts();
  void FixReferences(AliAODEvent *aod);
  Bool_t FillRecoCand(AliVEvent *event,AliAODRecoDecayHF3Prong *rd3);
  Bool_t FillRecoCand(AliVEvent *event,AliAODRecoDecayHF2Prong *rd2);
  Bool_t FillRecoCasc(AliVEvent *event,AliAODRecoCascadeHF *rc,Bool_t isDStar,Bool_t recoSecVtx=kFALSE);
  Bool_t RecoSecondaryVertexForCascades(AliVEvent *event, AliAODRecoCascadeHF *rc);
  void PrintStatus() const;
  void SetSecVtxWithKF() { fSecVtxWithKF=kTRUE; }
  void SetD0toKpiOn() { fD0toKpi=kTRUE; }
  void SetD0toKpiOff() { fD0toKpi=kFALSE; }
  void SetJPSItoEleOn() { fJPSItoEle=kTRUE; }
  void SetJPSItoEleOff() { fJPSItoEle=kFALSE; }
  void Set3ProngOn() { f3Prong=kTRUE; }
  void Set3ProngOff() { f3Prong=kFALSE; }
  void Set4ProngOn() { f4Prong=kTRUE; }
  void Set4ProngOff() { f4Prong=kFALSE; }
  void SetDstarOn() { fDstar=kTRUE; }
  void SetDstarOff() { fDstar=kFALSE; }
  void SetCascadesOn() { fCascades=kTRUE; }
  void SetCascadesOff() { fCascades=kFALSE; }
  void SetLikeSignOn() { fLikeSign=kTRUE; }
  void SetLikeSignOff() {fLikeSign=kFALSE; fLikeSign3prong=kFALSE;}
  void SetLikeSign3prongOn() { fLikeSign=kTRUE; fLikeSign3prong=kTRUE; }
  void SetLikeSign3prongOff() { fLikeSign3prong=kFALSE; }
  void SetMixEventOn() { fMixEvent=kTRUE; }
  void SetMixEventOff() { fMixEvent=kFALSE; }
  void SetInputAOD() { fInputAOD=kTRUE; }
  void SetMakeReducedRHF(Bool_t makeredAOD=kFALSE) { fMakeReducedRHF=makeredAOD; }
  Bool_t GetD0toKpi() const { return fD0toKpi; }
  Bool_t GetJPSItoEle() const { return fJPSItoEle; }
  Bool_t Get3Prong() const { return f3Prong; }
  Bool_t Get4Prong() const { return f4Prong; }
  Bool_t GetDstar()  const { return fDstar; }
  Bool_t GetCascades() const { return fCascades; }
  Bool_t GetLikeSign() const { return fLikeSign; }
  Bool_t GetLikeSign3prong() const { return fLikeSign3prong; }
  Bool_t GetMixEvent() const { return fMixEvent; }
  Bool_t GetInputAOD() const { return fInputAOD; }
  Bool_t GetRecoPrimVtxSkippingTrks() const {return fRecoPrimVtxSkippingTrks;}
  Bool_t GetRmTrksFromPrimVtx() const {return fRmTrksFromPrimVtx;}
  Bool_t GetMakeReducedRHF() const {return fMakeReducedRHF;}
  void SetFindVertexForDstar(Bool_t vtx=kTRUE) { fFindVertexForDstar=vtx; }
  void SetFindVertexForCascades(Bool_t vtx=kTRUE) { fFindVertexForCascades=vtx; }

  void  SetV0TypeForCascadeVertex(Int_t type) {fV0TypeForCascadeVertex = type;}
  Int_t GetV0TypeForCascadeVertex()           { return fV0TypeForCascadeVertex;}

  void SetRecoPrimVtxSkippingTrks()
    { fRecoPrimVtxSkippingTrks=kTRUE; fRmTrksFromPrimVtx=kFALSE;}
  void UnsetRecoPrimVtxSkippingTrks()
    { fRecoPrimVtxSkippingTrks=kFALSE; fRmTrksFromPrimVtx=kFALSE;}
  void SetRmTrksFromPrimVtx()
    {fRmTrksFromPrimVtx=kTRUE; fRecoPrimVtxSkippingTrks=kFALSE; }
  void SetTrackFilter(AliAnalysisFilter* trackF) {
    /// switch off the TOF selection that cannot be applied with AODTracks
    TList *l = (TList*)trackF->GetCuts();
    AliESDtrackCuts *tcuts = (AliESDtrackCuts*)l->FindObject("AliESDtrackCuts");
    if(tcuts->GetFlagCutTOFdistance()) tcuts->SetFlagCutTOFdistance(kFALSE);
    fTrackFilter = trackF;
  }
  void SetTrackFilter2prongPbCentral(Float_t maxPercentile, AliAnalysisFilter* trackF) {
    /// switch off the TOF selection that cannot be applied with AODTracks
    TList *l = (TList*)trackF->GetCuts();
    AliESDtrackCuts *tcuts = (AliESDtrackCuts*)l->FindObject("AliESDtrackCuts");
    if(tcuts->GetFlagCutTOFdistance()) tcuts->SetFlagCutTOFdistance(kFALSE);
    fTrackFilter2prongCentral = trackF;
    fMaxCentPercentileForTightCuts=maxPercentile;
  }
  void SetTrackFilter3prongPbCentral(Float_t maxPercentile, AliAnalysisFilter* trackF) {
    /// switch off the TOF selection that cannot be applied with AODTracks
    TList *l = (TList*)trackF->GetCuts();
    AliESDtrackCuts *tcuts = (AliESDtrackCuts*)l->FindObject("AliESDtrackCuts");
    if(tcuts->GetFlagCutTOFdistance()) tcuts->SetFlagCutTOFdistance(kFALSE);
    fTrackFilter3prongCentral = trackF;
    fMaxCentPercentileForTightCuts=maxPercentile;
  }
  void SetTrackFilterSoftPi(AliAnalysisFilter* trackF) {
    /// switch off the TOF selection that cannot be applied with AODTracks
    TList *l = (TList*)trackF->GetCuts();
    AliESDtrackCuts *tcuts = (AliESDtrackCuts*)l->FindObject("AliESDtrackCuts");
    if(tcuts->GetFlagCutTOFdistance()) tcuts->SetFlagCutTOFdistance(kFALSE);
    fTrackFilterSoftPi = trackF;
  }
  void SetTrackFilterBachelor(AliAnalysisFilter* trackF) {
    /// switch off the TOF selection that cannot be applied with AODTracks
    TList *l = (TList*)trackF->GetCuts();
    AliESDtrackCuts *tcuts = (AliESDtrackCuts*)l->FindObject("AliESDtrackCuts");
    if(tcuts->GetFlagCutTOFdistance()) tcuts->SetFlagCutTOFdistance(kFALSE);
    fTrackFilterBachelor = trackF;
  }
  AliAnalysisFilter* GetTrackFilter() const { return fTrackFilter; }
  AliAnalysisFilter* GetTrackFilterSoftPi() const { return fTrackFilterSoftPi; }
  void SetCutsD0toKpi(AliRDHFCutsD0toKpi* cuts) { fCutsD0toKpi = cuts; }
  AliRDHFCutsD0toKpi* GetCutsD0toKpi() const { return fCutsD0toKpi; }
  void SetCutsJpsitoee(AliRDHFCutsJpsitoee* cuts) { fCutsJpsitoee = cuts; }
  AliRDHFCutsJpsitoee* GetCutsJpsitoee() const { return fCutsJpsitoee; }
  void SetCutsDplustoK0spi(AliRDHFCutsDplustoK0spi* cuts) { fCutsDplustoK0spi = cuts; }
  AliRDHFCutsDplustoK0spi* GetCutsDplustoK0spi() const { return fCutsDplustoK0spi; }
  void SetCutsDplustoKpipi(AliRDHFCutsDplustoKpipi* cuts) { fCutsDplustoKpipi = cuts; }
  AliRDHFCutsDplustoKpipi* GetCutsDplustoKpipi() const { return fCutsDplustoKpipi; }
  void SetCutsDstoK0sK(AliRDHFCutsDstoK0sK* cuts) { fCutsDstoK0sK = cuts; }
  AliRDHFCutsDstoK0sK* GetCutsDstoK0sK() const { return fCutsDstoK0sK; }
  void SetCutsDstoKKpi(AliRDHFCutsDstoKKpi* cuts) { fCutsDstoKKpi = cuts; }
  AliRDHFCutsDstoKKpi* GetCutsDstoKKpi() const { return fCutsDstoKKpi; }
  void SetCutsLctopKpi(AliRDHFCutsLctopKpi* cuts) { fCutsLctopKpi = cuts; }
  AliRDHFCutsLctopKpi* GetCutsLctopKpi() const { return fCutsLctopKpi; }
  void SetCutsLctoV0(AliRDHFCutsLctoV0* cuts) { fCutsLctoV0 = cuts; }
  AliRDHFCutsLctoV0* GetCutsLctoV0() const { return fCutsLctoV0; }
  void SetCutsD0toKpipipi(AliRDHFCutsD0toKpipipi* cuts) { fCutsD0toKpipipi = cuts; }
  AliRDHFCutsD0toKpipipi* GetCutsD0toKpipipi() const { return fCutsD0toKpipipi; }
  void SetCutsDStartoKpipi(AliRDHFCutsDStartoKpipi* cuts) { fCutsDStartoKpipi = cuts; }
  AliRDHFCutsDStartoKpipi* GetCutsDStartoKpipi() const { return fCutsDStartoKpipi; }
  void SetMassCutBeforeVertexing(Bool_t flag) { fMassCutBeforeVertexing=flag; }

  void SetMasses();
  Bool_t CheckCutsConsistency();

  void SetUseTPCPID(Bool_t opt=kTRUE){fUseTPCPID=opt;}
  void SetUseTOFPID(Bool_t opt=kTRUE){fUseTOFPID=opt;}
  void SetUseTPCPIDOnlyIfNoTOF(Bool_t opt=kTRUE){fUseTPCPIDOnlyIfNoTOF=opt;}
  void SetMaxMomForTPCPid(Double_t mom){fMaxMomForTPCPid=mom;}
  void SetUsePidTag(Bool_t opt=kTRUE){fUsePidTag=opt;}
  void SetnSigmaTPCforPionSel(Double_t nsl, Double_t nsh){
    fnSigmaTPCPionLow=nsl; fnSigmaTPCPionHi=nsh;}
  void SetnSigmaTOFforPionSel(Double_t nsl, Double_t nsh){
    fnSigmaTOFPionLow=nsl; fnSigmaTOFPionHi=nsh;}
  void SetnSigmaTPCforKaonSel(Double_t nsl, Double_t nsh){
    fnSigmaTPCKaonLow=nsl; fnSigmaTPCKaonHi=nsh;}
  void SetnSigmaTOFforKaonSel(Double_t nsl, Double_t nsh){
    fnSigmaTOFKaonLow=nsl; fnSigmaTOFKaonHi=nsh;}
  void SetnSigmaTPCforProtonSel(Double_t nsl, Double_t nsh){
    fnSigmaTPCProtonLow=nsl; fnSigmaTPCProtonHi=nsh;}
  void SetnSigmaTOFforProtonSel(Double_t nsl, Double_t nsh){
    fnSigmaTOFProtonLow=nsl; fnSigmaTOFProtonHi=nsh;}

  void SetUseKaonPIDfor3Prong(Bool_t opt=kTRUE){fUseKaonPIDfor3Prong=opt;}
  void SetNotUseProtonPIDforLambdaC(){fUsePIDforLc=0;}
  void SetUseProtonPIDforLambdaC(){fUsePIDforLc=1;}
  void SetUseProtonAndPionPIDforLambdaC(){fUsePIDforLc=2;}
  void SetUseKaonPIDforDs(Bool_t opt=kTRUE){fUseKaonPIDforDs=opt;}

  void SetNotUseProtonPIDforLambdaC2V0(){fUsePIDforLc2V0=kFALSE;} //clm
  void SetUseProtonPIDforLambdaC2V0(){fUsePIDforLc2V0=kTRUE;}     //clm



  void GetnSigmaTOFforPionSel(Double_t& minnsigma, Double_t& maxnsigma) const {
    minnsigma=fnSigmaTOFPionLow;maxnsigma=fnSigmaTOFPionHi;
  }
  void GetnSigmaTPCforPionSel(Double_t& minnsigma, Double_t& maxnsigma) const {
    minnsigma=fnSigmaTPCPionLow;maxnsigma=fnSigmaTPCPionHi;
  }
  void GetnSigmaTOFforKaonSel(Double_t& minnsigma, Double_t& maxnsigma) const {
    minnsigma=fnSigmaTOFKaonLow;maxnsigma=fnSigmaTOFKaonHi;
  }
  void GetnSigmaTPCforKaonSel(Double_t& minnsigma, Double_t& maxnsigma) const {
    minnsigma=fnSigmaTPCKaonLow;maxnsigma=fnSigmaTPCKaonHi;
  }
  void GetnSigmaTOFforProtonSel(Double_t& minnsigma, Double_t& maxnsigma) const {
    minnsigma=fnSigmaTOFProtonLow;maxnsigma=fnSigmaTOFProtonHi;
  }
  void GetnSigmaTPCforProtonSel(Double_t& minnsigma, Double_t& maxnsigma) const {
    minnsigma=fnSigmaTPCProtonLow;maxnsigma=fnSigmaTPCProtonHi;
  }

  Bool_t GetUseTPCPID() const {return fUseTPCPID;}
  Bool_t GetUseTOFPID() const {return fUseTOFPID;}
  Bool_t GetUseTPCPIDOnlyIfNoTOF() const {return fUseTPCPIDOnlyIfNoTOF;}
  Double_t GetMaxMomForTPCPid() const {return fMaxMomForTPCPid;}

  Bool_t GetUseKaonPIDfor3Prong() const {return fUseKaonPIDfor3Prong;}
  Int_t GetUseProtonPIDforLambdaC() const {return fUsePIDforLc;}
  Bool_t GetUseKaonPIDforDs() const {return fUseKaonPIDforDs;}
  Bool_t GetUseProtonPIDforLambdaC2V0() const {return fUsePIDforLc2V0;}

  void SetPidResponse(AliPIDResponse* p){fPidResponse=p;}

  AliAODTrack *GetProng(AliVEvent *event,AliAODRecoDecayHF *rd,Int_t iprong);
  Bool_t GetTrackMomentumAtSecVert(AliESDtrack* tr, AliAODVertex* secVert, Double_t momentum[3]) const;

  //
 private:
  //
  enum { kBitDispl = 0, kBitSoftPi = 1, kBit3Prong = 2, kBitPionCompat = 3, kBitKaonCompat = 4, kBitProtonCompat = 5, kBitBachelor = 6};

  Bool_t fInputAOD; /// input from AOD (kTRUE) or ESD (kFALSE)
  Int_t fAODMapSize; /// size of fAODMap
  /// map between index and ID for AOD tracks
  Int_t *fAODMap; //[fAODMapSize] map between index and ID for AOD tracks

  AliVertexerTracks* fVertexerTracks; /// vertexer, to compute secondary vertices
  Double_t fBzkG; /// z componenent of field in kG

  Bool_t fSecVtxWithKF; /// if kTRUE use KF vertexer, else AliVertexerTracks

  Bool_t fRecoPrimVtxSkippingTrks; /// flag for primary vertex reco on the fly
                                   /// for each candidate, w/o its daughters
  Bool_t fRmTrksFromPrimVtx; /// flag for fast removal of daughters from
                             /// the primary vertex

  AliESDVertex *fV1; /// primary vertex
  AliAODVertex *fV1AOD; /// primary vertex (AOD format)

  /// flag to enable candidates production
  Bool_t fD0toKpi;   /// D0->Kpi
  Bool_t fJPSItoEle; /// Jpsi->ee
  Bool_t f3Prong;    /// D+,Ds,Lc
  Bool_t f4Prong;    /// D0->Kpipipi
  Bool_t fDstar;     /// D*->D0pi
  Bool_t fCascades;  /// cascades, Lc --> v0+track, D+ --> K0s+Pi, Ds --> K0s+K
  Bool_t fLikeSign;  /// Like-sign pairs
  Bool_t fLikeSign3prong;  /// Like-sign triplets
  Bool_t fMixEvent; /// event mixing

  AliPIDResponse* fPidResponse; /// PID response
  Bool_t fUseKaonPIDfor3Prong;  /// Kaon PID usage for 3 prongs
  Int_t  fUsePIDforLc;          /// PID for Lambdac: 0=no, 1=proton, 2=p and pi
  Bool_t fUsePIDforLc2V0;       /// PID for Lambdac 2 V0: 0=no, 1=proton,
  Bool_t fUseKaonPIDforDs;      /// Kaon PID usage for Ds
  Bool_t fUseTPCPID;            /// switch use/not use TPC PID
  Bool_t fUseTOFPID;            /// switch use/not use TOF PID
  Bool_t fUseTPCPIDOnlyIfNoTOF; /// use TPC PID only for tracks that without TOF
  Double_t fMaxMomForTPCPid;    /// upper momentum limit to apply TPC PID
  Bool_t fUsePidTag;            /// flag to control usage of PID tagging
  Double_t fnSigmaTPCPionLow;   /// Low cut value on n. of sigmas for pi TPC PID
  Double_t fnSigmaTPCPionHi;    /// High cut value on n. of sigmas for pi TPC PID
  Double_t fnSigmaTOFPionLow;   /// Low cut value on n. of sigmas for pi TOF PID
  Double_t fnSigmaTOFPionHi;    /// High cut value on n. of sigmas for pi TOF PID
  Double_t fnSigmaTPCKaonLow;   /// Low cut value on n. of sigmas for K TPC PID
  Double_t fnSigmaTPCKaonHi;    /// High cut value on n. of sigmas for K TPC PID
  Double_t fnSigmaTOFKaonLow;   /// Low cut value on n. of sigmas for K TOF PID
  Double_t fnSigmaTOFKaonHi;    /// High cut value on n. of sigmas for K TOF PID
  Double_t fnSigmaTPCProtonLow; /// Low cut value on n. of sigmas for p TPC PID
  Double_t fnSigmaTPCProtonHi;  /// High cut value on n. of sigmas for p TPC PID
  Double_t fnSigmaTOFProtonLow; /// Low cut value on n. of sigmas for p TOF PID
  Double_t fnSigmaTOFProtonHi;  /// High cut value on n. of sigmas for p TOF PID

  Float_t fMaxCentPercentileForTightCuts; /// max. centrality percentile for using tight cuts

  // single-track cuts
  AliAnalysisFilter *fTrackFilter; ///  Track Filter for displaced vertices
  AliAnalysisFilter *fTrackFilter2prongCentral; ///  Track Filter for displaced vertices in PbPb central events (tighter cuts) for 2 prong (D0->Kpi)
  AliAnalysisFilter *fTrackFilter3prongCentral; ///  Track Filter for displaced vertices in PbPb central events (tighter cuts) for 3 prong (D+, Ds, Lc)
  AliAnalysisFilter *fTrackFilterSoftPi; ///  Track Filter for D* soft pion
  AliAnalysisFilter *fTrackFilterBachelor; ///  Track Filter for bachelor
  // candidates cuts
  AliRDHFCutsD0toKpi *fCutsD0toKpi; /// D0->Kpi cuts
  AliRDHFCutsJpsitoee *fCutsJpsitoee; /// J/psi->ee cuts
  AliRDHFCutsDplustoK0spi *fCutsDplustoK0spi; /// D+->K0s+pi
  AliRDHFCutsDplustoKpipi *fCutsDplustoKpipi; /// D+->Kpipi cuts
  AliRDHFCutsDstoK0sK *fCutsDstoK0sK; /// Ds->K0s+K
  AliRDHFCutsDstoKKpi *fCutsDstoKKpi; /// Ds->KKpi cuts
  AliRDHFCutsLctopKpi *fCutsLctopKpi; /// Lc->pKpi cuts
  AliRDHFCutsLctoV0 *fCutsLctoV0; /// Lc --> v0 + bachelor cuts
  AliRDHFCutsD0toKpipipi *fCutsD0toKpipipi; /// D0->Kpipipi cuts
  AliRDHFCutsDStartoKpipi *fCutsDStartoKpipi; /// Dstar->D0pi cuts

  TList *fListOfCuts;    /// pointer to list of cuts for output file
  Bool_t fFindVertexForDstar; /// reconstruct a secondary vertex or assume it's from the primary vertex
  Bool_t fFindVertexForCascades;  /// reconstruct a secondary vertex or assume it's from the primary vertex
  Int_t  fV0TypeForCascadeVertex;  /// Select which V0 type we want to use for the cascas
  Bool_t fMassCutBeforeVertexing; /// to go faster in PbPb
  // dummies for invariant mass calculation
  AliAODRecoDecay *fMassCalc2; /// for 2 prong
  AliAODRecoDecay *fMassCalc3; /// for 3 prong
  AliAODRecoDecay *fMassCalc4; /// for 4 prong
  Double_t fMinPt3Prong; /// minimum pt for 3 prong candidates
  Bool_t fOKInvMassD0; /// pair fullfilling D0 inv mass selection
  Bool_t fOKInvMassJpsi; /// pair fullfilling Jpsi inv mass selection
  Bool_t fOKInvMassDplus; /// triplet fullfilling D+ inv mass selection
  Bool_t fOKInvMassDs; /// triplet fullfilling Ds inv mass selection
  Bool_t fOKInvMassLc; /// triplet fullfilling Lc inv mass selection
  Bool_t fOKInvMassDstar; /// combination fullfilling D* inv mass selection
  Bool_t fOKInvMassD0to4p; /// 4tracks fullfilling D0 inv mass selection
  Bool_t fOKInvMassLctoV0; /// triplet fullfilling Lc inv mass selection

  Int_t  fnTrksTotal;
  Int_t  fnSeleTrksTotal;
  Bool_t fMakeReducedRHF;// switch the reduction of dAOD size on/off

  Double_t fMassDzero;
  Double_t fMassDplus;
  Double_t fMassDs;
  Double_t fMassLambdaC;
  Double_t fMassDstar;
  Double_t fMassJpsi;
  Double_t fMassPhi;
  Double_t fMassK;

  //
  void AddRefs(AliAODVertex *v,AliAODRecoDecayHF *rd,const AliVEvent *event,
	       const TObjArray *trkArray) const;
  void AddDaughterRefs(AliAODVertex *v,const AliVEvent *event,
		       const TObjArray *trkArray) const;
  AliAODRecoDecayHF2Prong* Make2Prong(TObjArray *twoTrackArray1,AliVEvent *event,
				      AliAODVertex *secVert,Double_t dcap1n1,
				      Bool_t &okD0,Bool_t &okJPSI,Bool_t &okD0fromDstar,
				      Bool_t callFromCascade=kFALSE, Bool_t refill=kFALSE,
				      AliAODRecoDecayHF2Prong *rd=0x0);
  AliAODRecoDecayHF3Prong* Make3Prong(TObjArray *threeTrackArray,AliVEvent *event,
				      AliAODVertex *secVert,
				      Double_t dispersion,
				      const AliAODVertex *vertexp1n1,
				      TObjArray *twoTrackArray2,
				      Double_t dcap1n1,Double_t dcap2n1,Double_t dcap1p2,
				      Bool_t useForLc, Bool_t useForDs,
				      Bool_t &ok3Prong);
  AliAODRecoDecayHF3Prong* Make3Prong(TObjArray *threeTrackArray,AliVEvent *event,
                                      AliAODVertex *secVert,
                                      Double_t dispersion,
                                      Double32_t dist12,
                                      Double32_t dist23,
                                      Double_t dcap1n1,Double_t dcap2n1,Double_t dcap1p2,
                                      AliAODRecoDecayHF3Prong *rd);
  AliAODRecoDecayHF4Prong* Make4Prong(TObjArray *fourTrackArray,AliVEvent *event,
                                      AliAODVertex *secVert,
                                      const AliAODVertex *vertexp1n1,
                                      const AliAODVertex *vertexp1n1p2,
                                      Double_t dcap1n1,Double_t dcap1n2,
                                      Double_t dcap2n1,Double_t dcap2n2,
                                      Bool_t &ok4Prong);
  AliAODRecoCascadeHF* MakeCascade(TObjArray *twoTrackArray,AliVEvent *event,
				   AliAODVertex *secVert,
				   AliAODRecoDecayHF2Prong *rd2Prong,
				   Double_t dca,
				   Bool_t &okDstar);
  AliAODRecoCascadeHF* MakeCascade(TObjArray *twoTrackArray,AliVEvent *event,
				   AliAODVertex *secVert,
				   AliAODv0 *v0,
				   Double_t dca,
				   Bool_t &okCascades);

  void MapAODtracks(AliVEvent *aod);
  AliAODVertex* PrimaryVertex(const TObjArray *trkArray=0x0,AliVEvent *event=0x0) const;
  AliAODVertex* ReconstructSecondaryVertex(TObjArray *trkArray,Double_t &dispersion,Bool_t useTRefArray=kTRUE) const;

  Bool_t SelectInvMassAndPt3prong(Double_t *px,Double_t *py,Double_t *pz, Int_t pidLcStatus=3);
  Bool_t SelectInvMassAndPt4prong(Double_t *px,Double_t *py,Double_t *pz);
  Bool_t SelectInvMassAndPtD0Kpi(Double_t *px,Double_t *py,Double_t *pz);
  Bool_t SelectInvMassAndPtJpsiee(Double_t *px,Double_t *py,Double_t *pz);
  Bool_t SelectInvMassAndPtDstarD0pi(Double_t *px,Double_t *py,Double_t *pz);
  Bool_t SelectInvMassAndPtCascade(Double_t *px,Double_t *py,Double_t *pz);

  Bool_t SelectInvMassAndPt3prong(TObjArray *trkArray);
  Bool_t SelectInvMassAndPt4prong(TObjArray *trkArray);
  Bool_t SelectInvMassAndPtDstarD0pi(TObjArray *trkArray);
  Bool_t SelectInvMassAndPtCascade(TObjArray *trkArray);

  void   SelectTracksAndCopyVertex(const AliVEvent *event,Int_t trkEntries,
				   TObjArray &seleTrksArray,
				   TObjArray &tracksAtVertex,
				   Int_t &nSeleTrks,
				   UChar_t *seleFlags,Int_t *evtNumber);
  void SetParametersAtVertex(AliESDtrack* esdt, const AliExternalTrackParam* extpar) const;

  Bool_t SingleTrkCuts(AliESDtrack *trk,Float_t centralityperc, Bool_t &okDisplaced,Bool_t &okSoftPi, Bool_t &ok3prong, Bool_t &okBachelor) const;

  void   SetSelectionBitForPID(AliRDHFCuts *cuts,AliAODRecoDecayHF *rd,Int_t bit);

  AliAODv0* TransformESDv0toAODv0(AliESDv0 *esdv0,
				  TObjArray *twoTrackArrayV0);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisVertexingHF,30);  // Reconstruction of HF decay candidates
  /// \endcond
};

#endif
