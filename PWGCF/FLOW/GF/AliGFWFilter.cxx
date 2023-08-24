#include "AliGFWFilter.h"
using namespace GFWFlags;
AliGFWFilter::AliGFWFilter():
  fRetFlags(0),
  flTrFlag(0),
  fEventFlag(0),
  fTrackMasks({}),
  fNTrMasks(0),
  fEventMasks({}),
  fNEvMasks(0),
  fAODTrack(0),
  fEventCuts(0),
  fPtMin(0.2),
  fPtMax(3.0),
  fEtaMin(-0.8),
  fEtaMax(0.8)
{
};
AliGFWFilter::~AliGFWFilter() {
};
void AliGFWFilter::CleanUp() {
};
Bool_t AliGFWFilter::AcceptVertex(AliAODEvent *inEv, Double_t *lvtxXYZ) {
  const AliAODVertex* vtx = dynamic_cast<const AliAODVertex*>(inEv->GetPrimaryVertex());
  if(!vtx || vtx->GetNContributors() < 1)
    return kFALSE;
  const AliAODVertex* vtxSPD = dynamic_cast<const AliAODVertex*>(inEv->GetPrimaryVertexSPD());
  Double_t dMaxResol = 0.25; // suggested from DPG
  Double_t cov[6] = {0};
  vtxSPD->GetCovarianceMatrix(cov);
  Double_t zRes = TMath::Sqrt(cov[5]);
  if ( vtxSPD->IsFromVertexerZ() && (zRes > dMaxResol)) return kFALSE;
  const Double_t aodVtxZ = vtx->GetZ();
  vtx->GetXYZ(lvtxXYZ);
  return kTRUE;
};
void AliGFWFilter::CheckEvent(AliVEvent* inEv) {
  AliAODEvent *fAOD = dynamic_cast<AliAODEvent*>(inEv);
  // AliESDEvent *fESD = dynamic_cast<AliESDEvent*>(inEv); //ESD not implemented yet
  if(!fRetFlags) {fRetFlags = new AliGFWFlags(); }
  else fRetFlags->CleanUp();
  if(fAOD) {
    fEventFlag=0;
    //If event cuts are set, then only consider events that pass the event selection
    if(fEventCuts)  { if(fEventCuts->AcceptEvent(fAOD)) AddEv(klEventCuts); else return; }
    else AddEv(klEventCuts); //Otherwise, mark all as passed (assuming ES is done outside)
    Double_t vtxXYZ[3];
    if(AcceptVertex(fAOD,vtxXYZ)) AddEv(klVtxOK);
    else return; //If vertex ain't accepted, no point to continue
    Double_t lvtxz = TMath::Abs(fAOD->GetPrimaryVertex()->GetZ());
    if(lvtxz<5)  AddEv(klVtxZ5);
    if(lvtxz<7)  AddEv(klVtxZ7);
    if(lvtxz<9)  AddEv(klVtxZ9);
    if(lvtxz<10) AddEv(klVtxZ10);
    UInt_t evfl = calculateEventFlag();
    fRetFlags->SetEventFlags(evfl);
    Double_t pt, eta, trXYZ[3], trDCAxy;
    UShort_t nTPCCLSsh, nTPCcl;
    for(Int_t i=0;i<fAOD->GetNumberOfTracks();i++) {
      fAODTrack = (AliAODTrack*)fAOD->GetTrack(i);
      if(!fAODTrack) continue;
      pt = fAODTrack->Pt();
      eta = fAODTrack->Eta();
      if(pt<fPtMin || pt>fPtMax) continue;
      if(eta<fEtaMin || eta>fEtaMax) continue;
      flTrFlag=0;
      //Testing filter bits:
      if(fAODTrack->TestFilterBit(32)) AddTr(klFB32);
      if(fAODTrack->TestFilterBit(64)) AddTr(klFB64);
      if(fAODTrack->TestFilterBit(256)) AddTr(klFB256);
      if(fAODTrack->TestFilterBit(512)) AddTr(klFB512);
      //Testing shared clusters
      nTPCCLSsh = fAODTrack->GetTPCnclsS();
      nTPCcl    = fAODTrack->GetTPCncls();
      if(nTPCCLSsh*1.0/nTPCcl <= 0.4) AddTr(klSharedClusters);
      //Testing hit on first layer SDD
      if(fAODTrack->HasPointOnITSLayer(4)) AddTr(klHitOnSDD); //2xSPD, 2xSSD, 2xSDD
      if(!fAODTrack->HasPointOnITSLayer(0) && !fAODTrack->HasPointOnITSLayer(1)) AddTr(klNoSPD); //No hits in SPD -> Necessary for FB 96<->768 overlap
      //Calculating DCAs:
      fAODTrack->GetXYZ(trXYZ);
      trXYZ[0] -= vtxXYZ[0];
      trXYZ[1] -= vtxXYZ[1];
      trXYZ[2] -= vtxXYZ[2];
      trDCAxy = TMath::Sqrt(trXYZ[0]*trXYZ[0] + trXYZ[1]*trXYZ[1]);
      //Checking DCAz:
      if(TMath::Abs(trXYZ[2])<2) AddTr(klDCAz20);
      if(TMath::Abs(trXYZ[2])<1) AddTr(klDCAz10);
      if(TMath::Abs(trXYZ[2])<0.5) AddTr(klDCAz05);
      //Checking DCAxy:
      Double_t dcaCut2010 = f_DCAxy2010(pt);
      Double_t dcaCut2011 = f_DCAxy2011(pt);
      if(trDCAxy<dcaCut2010*7) AddTr(klDCAxy2010);
      if(trDCAxy<dcaCut2011*7) AddTr(klDCAxy2011);
      if(trDCAxy<dcaCut2011*4) AddTr(klDCAxy4Sigma);
      if(trDCAxy<dcaCut2011*8) AddTr(klDCAxy8Sigma);
      if(trDCAxy<dcaCut2011*10) AddTr(klDCAxy10Sigma);
      //Checking TPC chi2 per cluster:
      Double_t tpcChi2PerCluster = fAODTrack->GetTPCchi2perCluster();
      if(tpcChi2PerCluster<=2.) AddTr(klTPCchi2PC20);
      if(tpcChi2PerCluster<=2.5) AddTr(klTPCchi2PC25);
      if(tpcChi2PerCluster<=3.0) AddTr(klTPCchi2PC30);
      if(tpcChi2PerCluster<=4.0) AddTr(klTPCchi2PC40);
      //Checking number of TPC clusters:
      UShort_t nTPCCls = fAODTrack->GetTPCNclsF();
      if(nTPCCls>=70) AddTr(klNTPCcls70);
      if(nTPCCls>=80) AddTr(klNTPCcls80);
      if(nTPCCls>=90) AddTr(klNTPCcls90);
      if(nTPCCls>=100) AddTr(klNTPCcls100);
      //Derived track cuts -- only for diff. modifications of filter bits, where OR is required
      if(TSB(klFB32)||TSB(klFB64)) AddTr(klFB96);
      if(TSB(klFB96)&&TSB(klSharedClusters)) AddTr(klFB96Tuned); //Tuned to overlap with 768 (modified)
      if(TSB(klFB256)||TSB(klFB512)) AddTr(klFB768);
      if(TSB(klFB256)||TB(klFB512+klNoSPD+klHitOnSDD)) AddTr(klFB768Tuned); //Tuned to overlap with 96 (modified). Second part is that only second part requires hit in SDD
      UInt_t flagToAdd = calculateTrackFlag();
      if(flagToAdd) fRetFlags->AddTrackFlags(i,flagToAdd); //Only count ones that pass any cuts
    };
  };
};
//klFB32, klFB64, klFB256, klFB512, klFB96, klFB96Tuned, klFB768, klFB768Tuned,
//klSharedClusters, klHitOnSDD, -- included in *Tuned fb's already
//klDCAz20, klDCAz10, klDCAz05,
//klDCAxy2010, klDCAxy2011, klDCAxy8Sigma, klDCAxy4Sigma, klDCAxy10Sigma,
//klTPCchi2PC25, klTPCchi2PC20, klTPCchi2PC30, klTPCchi2PC40
//klNTPCcls70, klNTPCcls80, klNTPCcls90, klNTPCcls100
void AliGFWFilter::CreateStandardCutMasks(kLocalTrackFlags lStandardChi2Cut) {
  //Standard cuts:
  //Nominal -- FB96:
  fTrackMasks.push_back(klFB96 + klDCAz20 + klDCAxy2011 + lStandardChi2Cut + klNTPCcls70);
  //FB768:
  fTrackMasks.push_back(klFB768 + lStandardChi2Cut + klNTPCcls70);
  //FB96, |dcaZ| < 1:
  fTrackMasks.push_back(klFB96 + klDCAz10 + klDCAxy2011 + lStandardChi2Cut + klNTPCcls70);
  //FB96, |dcaZ| < 0.5:
  fTrackMasks.push_back(klFB96 + klDCAz05 + klDCAxy2011 + lStandardChi2Cut + klNTPCcls70);
  //FB96, |DCAxy| 4 sigma:
  fTrackMasks.push_back(klFB96 + klDCAz20 + klDCAxy4Sigma + lStandardChi2Cut + klNTPCcls70);
  //FB96, |DCAxy| 10 sigma:
  fTrackMasks.push_back(klFB96 + klDCAz20 + klDCAxy10Sigma + lStandardChi2Cut + klNTPCcls70);
  //FB96, chi2 per cluster <2:
  fTrackMasks.push_back(klFB96 + klDCAz20 + klDCAxy2011 + klTPCchi2PC20 + klNTPCcls70);
  //FB96, chi2 per cluster <3:
  fTrackMasks.push_back(klFB96 + klDCAz20 + klDCAxy2011 + klTPCchi2PC30 + klNTPCcls70);
  //FB96, Ntpc clusters > 80:
  fTrackMasks.push_back(klFB96 + klDCAz20 + klDCAxy2011 + lStandardChi2Cut + klNTPCcls80);
  //FB96, Ntpc clusters > 90:
  fTrackMasks.push_back(klFB96 + klDCAz20 + klDCAxy2011 + lStandardChi2Cut + klNTPCcls90);
  //FB96, Ntpc clusters > 100:
  fTrackMasks.push_back(klFB96 + klDCAz20 + klDCAxy2011 + lStandardChi2Cut + klNTPCcls100);
  //More reasonable ones:
  //Nominal, tuned to overlap with 768
  fTrackMasks.push_back(klFB96Tuned + klDCAz20 + klDCAxy2011 + lStandardChi2Cut + klNTPCcls70);
  //FB768 tuned to overlap with 96
  fTrackMasks.push_back(klFB768Tuned + klDCAz20 + klDCAxy2011 + lStandardChi2Cut + klNTPCcls70);
  //FB768 tuned, smaller DCAz cut
  fTrackMasks.push_back(klFB768Tuned + klDCAz05 + klDCAxy2011 + lStandardChi2Cut + klNTPCcls70);
  //FB768 tuned, DCAxy 4 sigma
  fTrackMasks.push_back(klFB768Tuned + klDCAz20 + klDCAxy4Sigma + lStandardChi2Cut + klNTPCcls70);
  //FB768 tuned, DCAxy 10 sigma
  fTrackMasks.push_back(klFB768Tuned + klDCAz20 + klDCAxy10Sigma + lStandardChi2Cut + klNTPCcls70);
  //FB768tunes tuned, TPC chi2 <2
  fTrackMasks.push_back(klFB768Tuned + klDCAz20 + klDCAxy2011 + klTPCchi2PC20 + klNTPCcls70);
  //FB768tunes tuned, TPC chi2 <3
  fTrackMasks.push_back(klFB768Tuned + klDCAz20 + klDCAxy2011 + klTPCchi2PC30 + klNTPCcls70);
  //FB768tunes tuned, nTPC clusters >90
  fTrackMasks.push_back(klFB768Tuned + klDCAz20 + klDCAxy2011 + lStandardChi2Cut + klNTPCcls90);
  //FB96 with 2010 AND 2011 DCAxy
  fTrackMasks.push_back(klFB96 + klDCAz20 + klDCAxy2010 + klDCAxy2011 + lStandardChi2Cut + klNTPCcls70);
  //Standard FB96 with chi2 < 2.5 cut -- for testing LHC15o_pass1, assuming that standard is 4.0
  fTrackMasks.push_back(klFB96 + klDCAz20 + klDCAxy2011 + klTPCchi2PC25 + klNTPCcls70);
  //Same as the nominal tracks, but need to replicate for NUA rebinning
  fTrackMasks.push_back(klFB96 + klDCAz20 + klDCAxy2011 + lStandardChi2Cut + klNTPCcls70);

  //Event cuts:
  // if(!fEventMasks) fEventMasks = new UInt_t[gNEventFlags];
  //Nominal:
  fEventMasks.push_back(klVtxZ10 + klEventCuts);
  //vtx z<9
  fEventMasks.push_back(klVtxZ9 + klEventCuts);
  //vtx z<7
  fEventMasks.push_back(klVtxZ7 + klEventCuts);
  //vtx z<5
  fEventMasks.push_back(klVtxZ5 + klEventCuts);

  fNTrMasks = (Int_t)fTrackMasks.size(); //fetch number of track cuts actually added, also so that we don't need to check the size for each track
  fNEvMasks = (Int_t)fEventMasks.size(); //same as above, but for events
};
void AliGFWFilter::AddCustomCuts(Bool_t cleanFirst, UInt_t lEv, UInt_t lTr) {
  if(cleanFirst) {fEventMasks.clear(); fTrackMasks.clear(); };
  fTrackMasks.push_back(lTr);
  fEventMasks.push_back(lEv);
  fNTrMasks = (Int_t)fTrackMasks.size(); //fetch number of track cuts actually added, also so that we don't need to check the size for each track
  fNEvMasks = (Int_t)fEventMasks.size(); //same as above, but for events
}
UInt_t AliGFWFilter::calculateEventFlag() {
  UInt_t retFlag = 0;
  for(Int_t i=0;i<fNEvMasks;i++) if(EB(fEventMasks[i])) retFlag|=(1<<i);
  return retFlag;
};
UInt_t AliGFWFilter::calculateTrackFlag() {
  UInt_t retFlag = 0;
  for(Int_t i=0;i<fNTrMasks;i++) if(TB(fTrackMasks[i])) retFlag|=(1<<i);
  return retFlag;
}
