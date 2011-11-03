#ifndef ALIANALYSISVERTEXINGHF_H
#define ALIANALYSISVERTEXINGHF_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//-------------------------------------------------------------------------
//                      Class AliAnalysisVertexingHF
//            Reconstruction of heavy-flavour decay candidates
//      
//  Origin: E.Bruna, G.E.Bruno, A.Dainese, F.Prino, R.Romita, X.M.Zhang
//  Contact: andrea.dainese@pd.infn.it
//-------------------------------------------------------------------------

#include <TNamed.h>

class TList;
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
class AliRDHFCutsDplustoKpipi;
class AliRDHFCutsDstoKKpi;
class AliRDHFCutsLctopKpi;
class AliRDHFCutsLctoV0;
class AliRDHFCutsD0toKpipipi;
class AliRDHFCutsDStartoKpipi;
class AliESDtrack;
class AliVEvent;
class AliAODVertex;
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
  void SetFindVertexForDstar(Bool_t vtx=kTRUE) { fFindVertexForDstar=vtx; }
  void SetFindVertexForCascades(Bool_t vtx=kTRUE) { fFindVertexForCascades=vtx; }
  void SetRecoPrimVtxSkippingTrks() 
    { fRecoPrimVtxSkippingTrks=kTRUE; fRmTrksFromPrimVtx=kFALSE;}
  void UnsetRecoPrimVtxSkippingTrks()
    { fRecoPrimVtxSkippingTrks=kFALSE; fRmTrksFromPrimVtx=kFALSE;}
  void SetRmTrksFromPrimVtx() 
    {fRmTrksFromPrimVtx=kTRUE; fRecoPrimVtxSkippingTrks=kFALSE; }
  void SetTrackFilter(AliAnalysisFilter* trackF) { fTrackFilter = trackF; }
  void SetTrackFilterSoftPi(AliAnalysisFilter* trackF) { fTrackFilterSoftPi = trackF; }
  AliAnalysisFilter* GetTrackFilter() const { return fTrackFilter; }
  AliAnalysisFilter* GetTrackFilterSoftPi() const { return fTrackFilterSoftPi; }
  void SetCutsD0toKpi(AliRDHFCutsD0toKpi* cuts) { fCutsD0toKpi = cuts; }
  AliRDHFCutsD0toKpi* GetCutsD0toKpi() const { return fCutsD0toKpi; }
  void SetCutsJpsitoee(AliRDHFCutsJpsitoee* cuts) { fCutsJpsitoee = cuts; }
  AliRDHFCutsJpsitoee* GetCutsJpsitoee() const { return fCutsJpsitoee; }
  void SetCutsDplustoKpipi(AliRDHFCutsDplustoKpipi* cuts) { fCutsDplustoKpipi = cuts; }
  AliRDHFCutsDplustoKpipi* GetCutsDplustoKpipi() const { return fCutsDplustoKpipi; }
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

  //
 private:
  //
  enum { kBitDispl = 0, kBitSoftPi = 1 };

  Bool_t fInputAOD; // input from AOD (kTRUE) or ESD (kFALSE) 
  Int_t fAODMapSize; // size of fAODMap 
  Int_t *fAODMap; //[fAODMapSize] map between index and ID for AOD tracks

  Double_t fBzkG; // z componenent of field in kG

  Bool_t fSecVtxWithKF; // if kTRUE use KF vertexer, else AliVertexerTracks

  Bool_t fRecoPrimVtxSkippingTrks; // flag for primary vertex reco on the fly
                                   // for each candidate, w/o its daughters
  Bool_t fRmTrksFromPrimVtx; // flag for fast removal of daughters from 
                             // the primary vertex

  AliESDVertex *fV1; // primary vertex

  // flag to enable candidates production
  Bool_t fD0toKpi;   // D0->Kpi 
  Bool_t fJPSItoEle; // Jpsi->ee
  Bool_t f3Prong;    // D+,Ds,Lc
  Bool_t f4Prong;    // D0->Kpipipi
  Bool_t fDstar;     // D*->D0pi
  Bool_t fCascades;  // cascades, Lc --> v0+track
  Bool_t fLikeSign;  // Like-sign pairs
  Bool_t fLikeSign3prong;  // Like-sign triplets
  Bool_t fMixEvent; // event mixing

  // single-track cuts
  AliAnalysisFilter *fTrackFilter; //  Track Filter for displaced vertices
  AliAnalysisFilter *fTrackFilterSoftPi; //  Track Filter for D* soft pion
  // candidates cuts
  AliRDHFCutsD0toKpi *fCutsD0toKpi; // D0->Kpi cuts
  AliRDHFCutsJpsitoee *fCutsJpsitoee; // J/psi->ee cuts
  AliRDHFCutsDplustoKpipi *fCutsDplustoKpipi; // D+->Kpipi cuts
  AliRDHFCutsDstoKKpi *fCutsDstoKKpi; // Ds->KKpi cuts
  AliRDHFCutsLctopKpi *fCutsLctopKpi; // Lc->pKpi cuts
  AliRDHFCutsLctoV0 *fCutsLctoV0; // Lc --> v0 + bachelor cuts
  AliRDHFCutsD0toKpipipi *fCutsD0toKpipipi; // D0->Kpipipi cuts
  AliRDHFCutsDStartoKpipi *fCutsDStartoKpipi; // Dstar->D0pi cuts

  TList *fListOfCuts;    // pointer to list of cuts for output file
  Bool_t fFindVertexForDstar; // reconstruct a secondary vertex or assume it's from the primary vertex
  Bool_t fFindVertexForCascades;  // reconstruct a secondary vertex or assume it's from the primary vertex
  Bool_t fMassCutBeforeVertexing; // to go faster in PbPb
  // dummies for invariant mass calculation
  AliAODRecoDecay *fMassCalc2; // for 2 prong
  AliAODRecoDecay *fMassCalc3; // for 3 prong
  AliAODRecoDecay *fMassCalc4; // for 4 prong

  Int_t fnTrksTotal;
  Int_t fnSeleTrksTotal;

  //
  void AddRefs(AliAODVertex *v,AliAODRecoDecayHF *rd,const AliVEvent *event,
	       const TObjArray *trkArray) const;
  void AddDaughterRefs(AliAODVertex *v,const AliVEvent *event,
		       const TObjArray *trkArray) const;
  AliAODRecoDecayHF2Prong* Make2Prong(TObjArray *twoTrackArray1,AliVEvent *event,
				      AliAODVertex *secVert,Double_t dcap1n1,
				      Bool_t &okD0,Bool_t &okJPSI,Bool_t &okD0fromDstar) const;
  AliAODRecoDecayHF3Prong* Make3Prong(TObjArray *threeTrackArray,AliVEvent *event,
				      AliAODVertex *secVert,
				      Double_t dispersion,
				      const AliAODVertex *vertexp1n1,
				      const AliAODVertex *vertexp2n1,
				      Double_t dcap1n1,Double_t dcap2n1,Double_t dcap1p2,
				      Bool_t &ok3Prong) const;
  AliAODRecoDecayHF4Prong* Make4Prong(TObjArray *fourTrackArray,AliVEvent *event,
                                      AliAODVertex *secVert,
                                      const AliAODVertex *vertexp1n1,
                                      const AliAODVertex *vertexp1n1p2,
                                      Double_t dcap1n1,Double_t dcap1n2,
                                      Double_t dcap2n1,Double_t dcap2n2,
                                      Bool_t &ok4Prong) const;
  AliAODRecoCascadeHF* MakeCascade(TObjArray *twoTrackArray,AliVEvent *event,
				   AliAODVertex *secVert,
				   AliAODRecoDecayHF2Prong *rd2Prong,
				   Double_t dca,
				   Bool_t &okDstar) const;
  AliAODRecoCascadeHF* MakeCascade(TObjArray *twoTrackArray,AliVEvent *event,
				   AliAODVertex *secVert,
				   AliAODv0 *v0,
				   Double_t dca,
				   Bool_t &okCascades) const;

  AliAODVertex* PrimaryVertex(const TObjArray *trkArray=0x0,AliVEvent *event=0x0) const;
  AliAODVertex* ReconstructSecondaryVertex(TObjArray *trkArray,Double_t &dispersion,Bool_t useTRefArray=kTRUE) const;
  Bool_t SelectInvMassAndPt(Int_t decay,Int_t nprongs,
		       Double_t *px,Double_t *py,Double_t *pz) const;
  Bool_t SelectInvMassAndPt(TObjArray *trkArray) const;
  void   SelectTracksAndCopyVertex(const AliVEvent *event,Int_t trkEntries,
				   TObjArray &seleTrksArray,
				   TObjArray &tracksAtVertex,
				   Int_t &nSeleTrks,
				   UChar_t *seleFlags,Int_t *evtNumber);
  void SetParametersAtVertex(AliESDtrack* esdt, AliExternalTrackParam* extpar) const;

  Bool_t SingleTrkCuts(AliESDtrack *trk,Bool_t &okDisplaced,Bool_t &okSoftPi) const;

  void   SetSelectionBitForPID(AliRDHFCuts *cuts,AliAODRecoDecayHF *rd,Int_t bit);

  AliAODv0* TransformESDv0toAODv0(AliESDv0 *esdv0, 
				  TObjArray *twoTrackArrayV0);

  //
  ClassDef(AliAnalysisVertexingHF,19);  // Reconstruction of HF decay candidates
};


#endif








