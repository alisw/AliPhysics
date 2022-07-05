#include "AliAnalysisTaskJetQ.h"

AliAnalysisTaskJetQ::AliAnalysisTaskJetQ():
 AliAnalysisTaskSE(),
 fAOD(0),
 fMCEvent(0),
 fPoolMgr(0),
 fPoolTrackArray(0),
 fTriggerType(AliVEvent::kINT7),
 fOutList(0),
 fCentAxis(0),
 fVzAxis(0),
 fPtAxis(0),
 fNormCounter(0),
 fCorrPlot(0),
 fMixCorrPlot(0),
 fPtAssocMin(2.),
 fPtAssocMax(2.5),
 fPtTriggMin(6.),
 fPtTriggMax(8.)
{}
//_____________________________________________________________________________
AliAnalysisTaskJetQ::AliAnalysisTaskJetQ(const char* name):
 AliAnalysisTaskSE(name),
 fAOD(0),
 fMCEvent(0),
 fPoolMgr(0),
 fPoolTrackArray(0),
 fTriggerType(AliVEvent::kINT7),
 fOutList(0),
 fCentAxis(0),
 fVzAxis(0),
 fPtAxis(0),
 fNormCounter(0),
 fCorrPlot(0),
 fMixCorrPlot(0),
 fPtAssocMin(2.),
 fPtAssocMax(2.5),
 fPtTriggMin(6.),
 fPtTriggMax(8.)
{
    DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskJetQ::~AliAnalysisTaskJetQ()
{
  delete fPoolTrackArray;
}
//_____________________________________________________________________________
void AliAnalysisTaskJetQ::UserCreateOutputObjects()
{
    OpenFile(1);
    //Setting up variables
    if(fCentBins.size()<2) {
      Double_t tempCentBins[] = {0,100};
      SetCentralityBins(1,tempCentBins);
    };
    if(fVzBins.size()<2) {
      Double_t vzs[]   = {-10.,10.};
      SetVtxZBins(1,vzs);
    };
    if(fPtBins.size()<2) {
      Double_t ptbs[] = {2.,2.5};
      SetPtBins(1,ptbs);
    };
    fPtAssocMin = fPtAxis->GetBinLowEdge(1);
    fPtAssocMax = fPtAxis->GetBinUpEdge(fPtAxis->GetNbins());
    fPtDif = fPtBins.size()>2; //if only one pT bin, then assoc is not pt-dif. and use TH2 instead of TH3
    if(!fEvMixPars[0]) SetEventMixingCapacity(30,5000,10,10);
    //Setting up the pool manager
    // Double_t cents[] = {0.,50.,70.};
    fPoolMgr = new AliEventPoolManager(fEvMixPars[0], fEvMixPars[1], fCentBins.size()-1, fCentBins.data(), fVzBins.size()-1, fVzBins.data());
      //fNCentBins,fCentBins.data(), fNzVtxBins, fzVtxBins.data());
    if (!fPoolMgr) { AliError("Event Pool manager not created!"); return; }
    fPoolMgr->SetTargetValues(fEvMixPars[1], fEvMixPars[2]*0.01, fEvMixPars[3]);
    fPoolTrackArray = new TObjArray();
    fPoolTrackArray->SetOwner(kTRUE);
    // creating output lists
    fOutList = new TList();
    fOutList->SetOwner(kTRUE);
    TH1D *vtzBefore = new TH1D("vtzBefore","vtzBefore",20,-10,10);
    TH1D *vtzAfter  = new TH1D("vtzAfter","vtzAfter",20,-10,10);
    fNormCounter    = new TH2D("NormCounter","NormCounter; multi/cent; index",fCentBins.size()-1, fCentBins.data(), 4, -0.5, 3.5);
    fOutList->Add(vtzBefore);
    fOutList->Add(vtzAfter);
    fOutList->Add(fNormCounter);
    //Setting up correlation plots
    fCorrPlot    = new TH1**[fCentAxis->GetNbins()];
    fMixCorrPlot = new TH1**[fCentAxis->GetNbins()];
    for(Int_t iCent=0;iCent<fCentAxis->GetNbins();iCent++) {
      fCorrPlot[iCent] = new TH1*[fPtAxis->GetNbins()];
      fMixCorrPlot[iCent] = new TH1*[fPtAxis->GetNbins()];
      for(Int_t iPt=0;iPt<fPtAxis->GetNbins();iPt++) {
        fCorrPlot[iCent][iPt]       = new TH3D(Form("fCorr_Cent%i_Pt%i",iCent,iPt),Form("fCorr_Cent%i_Pt%i; #Delta#phi (rad); #Delta#eta; z_{vtx}",iCent,iPt),100,-C_PI_HALF,C_PI_TH,200, -1.6, 1.6,fVzBins.size()-1,-10,10);
        fMixCorrPlot[iCent][iPt]    = new TH3D(Form("fMixCorr_Cent%i_Pt%i",iCent,iPt),Form("fMixCorr_Cent%i_Pt%i; #Delta#phi (rad); #Delta#eta; z_{vtx}",iCent,iPt),100,-C_PI_HALF,C_PI_TH,200, -1.6, 1.6,fVzBins.size()-1,-10,10);
        fCorrPlot[iCent][iPt]->GetZaxis()->Set(fVzBins.size()-1,fVzBins.data());
        fMixCorrPlot[iCent][iPt]->GetZaxis()->Set(fVzBins.size()-1,fVzBins.data());
        fOutList->Add(fCorrPlot[iCent][iPt]);
        fOutList->Add(fMixCorrPlot[iCent][iPt]);
      };
    };
    PostData(1, fOutList);
}
//_____________________________________________________________________________
void AliAnalysisTaskJetQ::UserExec(Option_t *)
{
  Bool_t fIsMC = kFALSE;
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD) return;
  if(fIsMC) {
    fMCEvent = dynamic_cast<AliMCEvent *>(MCEvent());
    if (!fMCEvent) return;
  }
  Double_t vz = fAOD->GetPrimaryVertex()->GetZ();
  AliMultSelection *lMultSel = (AliMultSelection*)fInputEvent->FindListObject("MultSelection");
  Double_t l_Cent = lMultSel->GetMultiplicityPercentile("V0M");
  if(!CheckTrigger(l_Cent)) return;
  Double_t vtxXYZ[] = {0.,0.,0.};
  if(!AcceptAOD(fAOD, vtxXYZ)) return;
  // Double_t vz = fAOD->GetPrimaryVertex()->GetZ();
  ((TH1D*)fOutList->At(0))->Fill(vz);
  Int_t ind = FindGivenPt(fPtTriggMin,fPtTriggMax);
  if(ind<0) return;
  Int_t i_Cent = fCentAxis->FindBin(l_Cent);
  Int_t i_vz   = fVzAxis->FindBin(vz);
  if(!i_Cent || i_Cent>fCentAxis->GetNbins()) return; //out of centrality/vz range
  if(!i_vz || i_vz>fVzAxis->GetNbins()) return; //out of centrality/vz range
  ((TH1D*)fOutList->At(1))->Fill(vz);


  fNormCounter->Fill(l_Cent,0); //Number of triggers
  Int_t nPairs = FillCorrelations(ind,i_Cent,vz);
  AliEventPool *pool = fPoolMgr->GetEventPool(i_Cent-1, i_vz-1);
  if(!pool) { printf("Could not find the event pool!\n"); return; };
  // printf("Current numbe of events in pool: %i\n",pool->GetCurrentNEvents());
  if(pool->IsReady()) { Int_t nMixPairs = FillMixedEvent(ind, pool, i_Cent, vz);  };
  pool->UpdatePool((TObjArray*)fPoolTrackArray->Clone());
  fPoolTrackArray->Clear();
  PostData(1, fOutList);
}
//_____________________________________________________________________________
void AliAnalysisTaskJetQ::Terminate(Option_t *)
{}
//_____________________________________________________________________________

Bool_t AliAnalysisTaskJetQ::CheckTrigger(Double_t l_Cent) {
  UInt_t fSelMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  //Apparently, MB trigger can also mark special triggers, leaving depleted regions in multi. To avoid this, pass true, if MB has been triggered.
  //This would fail if spec. triggers would also flag MB trigger, which seems to NOT be the case.
  if(!(fTriggerType&fSelMask)) { return kFALSE; }; //printf("Returning from the generic check\n");
  if(fSelMask&(fTriggerType&(AliVEvent::kINT7+AliVEvent::kMB))) {return kTRUE; }; //printf("Passed by MB trigger!\n");
  if((fSelMask&fTriggerType&AliVEvent::kCentral) && l_Cent>10) {return kFALSE; }; //printf("Returnning from kCent case\n");
  if((fSelMask&fTriggerType&AliVEvent::kSemiCentral) && (l_Cent<30 || l_Cent>50)) {return kFALSE; }; //printf("Returning from kSC case\n");
  return kTRUE;
}
Bool_t AliAnalysisTaskJetQ::AcceptAOD(AliAODEvent *inEv, Double_t *lvtxXYZ) {
  if(!fEventCuts.AcceptEvent(inEv)) return 0;
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
  if(TMath::Abs(aodVtxZ) > 10)
    return kFALSE;
  vtx->GetXYZ(lvtxXYZ);
  return kTRUE;
};
Int_t AliAnalysisTaskJetQ::FindGivenPt(const Double_t &ptMin, const Double_t &ptMax) {
  Int_t ind=-1;
  AliAODTrack *lTrack;
  Double_t lPtMax=0;
  for(Int_t i=0;i<fAOD->GetNumberOfTracks();i++) {
    lTrack = (AliAODTrack*)fAOD->GetTrack(i);
    if(!lTrack->TestFilterBit(96)) continue;
    lPtMax = lTrack->Pt();
    if(lPtMax > ptMin && lPtMax < ptMax) { ind=i; break; };
  };
  return ind;
}
Int_t AliAnalysisTaskJetQ::FillCorrelations(Int_t &triggerIndex, Int_t &centVal, Double_t &vzValue) {
  Int_t nPairs=0;
  AliAODTrack *lTrack;
  AliAODTrack *lTriggerTrack = (AliAODTrack*)fAOD->GetTrack(triggerIndex);
  Double_t l_TrEta = lTriggerTrack->Eta();
  Double_t l_TrPhi = lTriggerTrack->Phi();
  Int_t centBin = centVal-1;
  for(Int_t i=0;i<fAOD->GetNumberOfTracks();i++) {
    lTrack = (AliAODTrack*)fAOD->GetTrack(i);
    if(!lTrack->TestFilterBit(96)) continue;
    Int_t ptInd = fPtAxis->FindBin(lTrack->Pt());
    if(!ptInd || ptInd>fPtAxis->GetNbins()) continue;
    ptInd--;
    Double_t d_Eta = l_TrEta-lTrack->Eta();
    Double_t d_Phi = l_TrPhi-lTrack->Phi();
    fixPhi(d_Phi);
    // printf("Attempting to fill a track with %f %f %f (%s)\n",d_Phi,d_Eta,lPt,fPtDif);
    fill3DHist(fCorrPlot[centBin][ptInd],d_Phi,d_Eta,vzValue);
    nPairs++;
    fPoolTrackArray->Add(new AliBasicParticle(lTrack->Eta(),lTrack->Phi(),lTrack->Pt(), lTrack->Charge()));
  };
  return nPairs;
};
Int_t AliAnalysisTaskJetQ::FillMixedEvent(Int_t &triggerIndex, AliEventPool *l_pool, Int_t &centVal, Double_t &vzValue) {
  Int_t nPairs=0;
  AliAODTrack *lTrack;
  AliAODTrack *lTriggerTrack = (AliAODTrack*)fAOD->GetTrack(triggerIndex);
  Double_t l_TrEta = lTriggerTrack->Eta();
  Double_t l_TrPhi = lTriggerTrack->Phi();
  Double_t l_AsEta;
  Double_t l_AsPhi;
  Double_t l_AsPt;
  TObjArray *mixTracks;
  Int_t centInd = centVal-1;
  for(Int_t i=0;i<l_pool->GetCurrentNEvents();i++) {
    mixTracks = l_pool->GetEvent(i);
    TObjArrayIter *myIt = new TObjArrayIter(mixTracks);
    while(AliBasicParticle *rTr = (AliBasicParticle*)myIt->Next()) {
      Int_t ptInd = fPtAxis->FindBin(rTr->Pt());
      if(!ptInd || ptInd>fPtAxis->GetNbins()) continue;
      ptInd-=1;
      l_AsEta = l_TrEta - rTr->Eta();
      l_AsPhi = l_TrPhi - rTr->Phi();
      fixPhi(l_AsPhi);
      fill3DHist(fMixCorrPlot[centInd][ptInd],l_AsPhi,l_AsEta,vzValue);
      nPairs++;
    };
  };
  return nPairs;
};
