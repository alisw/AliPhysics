#include "AliEffFDContainer.h"
AliEffFDContainer::AliEffFDContainer():
  TNamed("",""),
  fOutList(0),
  fCutList(0),
  fChi2Cut(1e10),
  fDCAXYPtCut(0),
  fIsMC(kFALSE),
  flMCEvent(0),
  flESDEvent(0),
  flMultSel(0),
  flmcWeightsHandler(0),
  flMCSpectraWeights(0),
  fEff(0),
  fDCA(0),
  fWithinDCA(0),
  fPtBins(0),
  fNPtBins(0),
  fCentBins(0),
  fNCentBins(0),
  fCentEst("V0M"),
  fCent(-999),
  fPtMin(0.2),
  fPtMax(3.0),
  fEta(0.8),
  htmp(0)
{
  fOutList = new TList();
  fOutList->SetOwner(kTRUE);
  fCutList = new TList();
  fCutList->SetOwner(kTRUE);
};
AliEffFDContainer::AliEffFDContainer(TString lName, TString lTitle):
  TNamed(lName,lTitle),
  fOutList(0),
  fCutList(0),
  fChi2Cut(1e10),
  fDCAXYPtCut(0),
  fIsMC(kFALSE),
  flMCEvent(0),
  flESDEvent(0),
  flMultSel(0),
  flmcWeightsHandler(0),
  flMCSpectraWeights(0),
  fEff(0),
  fDCA(0),
  fWithinDCA(0),
  fPtBins(0),
  fNPtBins(0),
  fCentBins(0),
  fNCentBins(0),
  fCentEst("V0M"),
  fCent(-999),
  fPtMin(0.2),
  fPtMax(3.0),
  fEta(0.8),
  htmp(0)
{
  fOutList = new TList();
  fOutList->SetOwner(kTRUE);
  fOutList->SetName(lName.Data());
  fCutList = new TList();
  fCutList->SetOwner(kTRUE);
  fCutList->SetName(Form("%s_Cuts",lName.Data()));
  //For testing purposes
  htmp = new TH1D("h1","h1",10,0,10);
  fOutList->Add(htmp);
};
void AliEffFDContainer::NewEvent(AliVEvent &InputEventAddress, AliMCEvent &MCEventAddress) {
  flMCEvent = &MCEventAddress;
  flESDEvent = dynamic_cast<AliESDEvent*>(&InputEventAddress);
  flMultSel= dynamic_cast<AliMultSelection*>(flESDEvent->FindListObject("MultSelection"));
  flmcWeightsHandler = static_cast<AliMCSpectraWeightsHandler*>(flESDEvent->FindListObject("fMCSpectraWeights"));
  flMCSpectraWeights = (flmcWeightsHandler) ? flmcWeightsHandler->fMCSpectraWeight : nullptr;
  fCent = flMultSel->GetMultiplicityPercentile(fCentEst.Data());
};
void AliEffFDContainer::Fill(AliVEvent &InputEventAddress, AliMCEvent &MCEventAddress) {
  NewEvent(InputEventAddress,MCEventAddress);
  printf("Address comparison: %p vs %p\n",flESDEvent,(void*)&InputEventAddress);
  if(!flESDEvent) {printf("ESD event not set! Not filling...\n"); return; };
  if(!flmcWeightsHandler) {printf("MCWeightsHandler not found! Not filling...\n"); return; };
  if(!flMCSpectraWeights) {printf("fMCSpectraWeights not found! Not filling...\n"); return; };
  Int_t nPrimPart = flMCEvent->GetNumberOfTracks();
  Int_t nTracks = flESDEvent->GetNumberOfTracks();
  AliMCParticle *lPart;
  AliESDtrack *lTrack;
  Double_t pt, eta;
  Double_t CompWeight;
  //Particle loop
  for(Int_t i=0;i<nPrimPart;i++) {
    if(AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i, flMCEvent)) continue;
    lPart = (AliMCParticle*)flMCEvent->GetTrack(i);
    if(!lPart->IsPhysicalPrimary()) continue;
    eta = lPart->Eta();
    if(TMath::Abs(eta)>fEta) continue;
    pt = lPart->Pt();
    if(pt<fPtMin || pt>fPtMax) continue;
    CompWeight = flMCSpectraWeights->GetMCSpectraWeightNominal(lPart->Particle());
    //Proceed with filling
    fEff[0]->Fill(pt,fCent,CompWeight);
    fEff[2]->Fill(pt,fCent);
  };
  //Track loop
  for(Int_t i=0;i<nTracks;i++) {
    lTrack = (AliESDtrack*)flESDEvent->GetTrack(i);
    if(!lTrack) continue;
    //Check track cuts
    Bool_t passTC=kFALSE;
    Bool_t passDCA=kFALSE;
    Bool_t passChi2=kFALSE;
    Float_t dcaXY, dcaZ;
    lTrack->GetImpactParameters(dcaXY,dcaZ);
    eta = lTrack->Eta();
    if(TMath::Abs(eta)>fEta) continue;
    pt  = lTrack->Pt();
    if(pt<fPtMin || pt>fPtMax) continue;
    for(Int_t iTC=0;iTC<fCutList->GetEntries();iTC++) {
      if(((AliESDtrackCuts*)fCutList->At(iTC))->AcceptTrack(lTrack)) passTC=kTRUE;
      if(passTC) break; //as soon as one of the track cuts have passed, we're done here
    }
    if(!passTC) continue; //If track cuts not passed, then don't bother anymore
    if(fDCAXYPtCut) passDCA = TMath::Abs(dcaXY)<fDCAXYPtCut->Eval(pt); else passDCA=kTRUE;
    Double_t lChi2Constrained = GetChi2TPCConstrained(lTrack);
    if(fChi2Cut<1e9) passChi2 = (lChi2Constrained<fChi2Cut && lChi2Constrained>=0); else passChi2=kTRUE;
    //Fetch corresponding particle
    Int_t fLabel = lTrack->GetLabel();
    Int_t index = TMath::Abs(fLabel);
    if(AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(index, flMCEvent)) continue;
    if (index < 0) continue;
    lPart = (AliMCParticle*)flMCEvent->GetTrack(index);
    if(!lPart) continue;
    eta = lPart->Eta();
    if(TMath::Abs(eta)>fEta) continue;
    pt = lPart->Pt();
    if(pt<fPtMin || pt>fPtMax) continue;
    CompWeight = flMCSpectraWeights->GetMCSpectraWeightNominal(lPart->Particle());
    if(lPart->IsPhysicalPrimary()) {
      fDCA[0]->Fill(pt,0.,dcaXY,CompWeight);
      if(passChi2) fDCA[1]->Fill(pt,0.,dcaXY,CompWeight);
      if(passDCA) {
        fWithinDCA[0]->Fill(pt,0.,CompWeight);
        if(passChi2) {
          fWithinDCA[1]->Fill(pt,0.,CompWeight);
          fEff[1]->Fill(pt,fCent,CompWeight);
          fEff[3]->Fill(pt,fCent);
        };
      };
    } else if(lPart->IsSecondaryFromWeakDecay()) {
      fDCA[0]->Fill(pt,1.,dcaXY,CompWeight);
      if(passChi2) fDCA[1]->Fill(pt,1.,dcaXY,CompWeight);
      if(passDCA) {
        fWithinDCA[0]->Fill(pt,1.,CompWeight);
        if(passChi2) fWithinDCA[1]->Fill(pt,1.,CompWeight);
      };
    } else if(lPart->IsSecondaryFromMaterial()) {
      fDCA[0]->Fill(pt,2.,dcaXY,CompWeight);
      if(passChi2) fDCA[1]->Fill(pt,2.,dcaXY,CompWeight);
      if(passDCA) {
        fWithinDCA[0]->Fill(pt,2.,CompWeight);
        if(passChi2) fWithinDCA[1]->Fill(pt,2.,CompWeight);
      };
    };
  };
};
void AliEffFDContainer::StoreBins(Int_t nBins, Double_t *lBins, Int_t &tNBins, Double_t *&tBins) {
  if(tBins) delete [] tBins;
  tNBins = nBins;
  tBins = new Double_t[nBins+1];
  for(Int_t i=0;i<=nBins;i++) tBins[i] = lBins[i];
};
void AliEffFDContainer::CreateHistograms() {
  if(!fCentBins) {
    Double_t centBins[] = {0,100};
    SetCentralityBins(1,centBins);
  };
  if(!fPtBins) {
    Double_t ptBins[] = {0,1,2,3,4,5,6,7,8,9,10};
    SetPtBins(10,ptBins);
  }
  fEff = new TH2D*[4];
  fEff[0] = new TH2D("nChGen_Weighted","ChGen_Weighted; #it{p}_{T} (GeV/#it{c}); Cent.",fNPtBins,fPtBins,fNCentBins,fCentBins);
  fEff[1] = new TH2D("nChRec_Weighted","ChGen_Weighted; #it{p}_{T} (GeV/#it{c}); Cent.",fNPtBins,fPtBins,fNCentBins,fCentBins);
  fEff[2] = new TH2D("nChGen_Uneighted","ChGen_Weighted; #it{p}_{T} (GeV/#it{c}); Cent.",fNPtBins,fPtBins,fNCentBins,fCentBins);
  fEff[3] = new TH2D("nChRec_Uneighted","ChGen_Weighted; #it{p}_{T} (GeV/#it{c}); Cent.",fNPtBins,fPtBins,fNCentBins,fCentBins);
  for(Int_t i=0;i<4;i++) { fOutList->Add(fEff[i]); };
  //Setting up DCAs
  Double_t binsDCA[79] = {-3.00, -2.90, -2.80, -2.70, -2.60, -2.50, -2.40, -2.30, -2.20, -2.10, -2.00, -1.90, -1.80, -1.70, -1.60, -1.50, -1.40, -1.30, -1.20, -1.10, -1.00, -0.90, -0.80, -0.70, -0.60, -0.50, -0.40, -0.30, -0.20, -0.10, -0.09, -0.08, -0.07, -0.06, -0.05, -0.04, -0.03, -0.02, -0.01, 0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90, 3.00};
  Int_t NbinsDCA = 78;
  Double_t binsType[4] = {-0.5,0.5,1.5,2.5};
  Int_t NbinsType = 3;
  fDCA = new TH3D*[2];
  fWithinDCA = new TH2D*[2];
  fDCA[0] = new TH3D("DCAxy_noChi2","DCAxy_noChi2; pt; type; dcaxy", fNPtBins, fPtBins, NbinsType, binsType, NbinsDCA, binsDCA);
  fDCA[1] = new TH3D("DCAxy_withChi2","DCAxy_withChi2; pt; type; dcaxy", fNPtBins, fPtBins, NbinsType, binsType, NbinsDCA, binsDCA);
  fWithinDCA[0] = new TH2D("WithinDCA_noChi2","WithinDCA_noChi2; pt; type", fNPtBins, fPtBins, NbinsType, binsType);
  fWithinDCA[1] = new TH2D("WithinDCA_withChi2","WithinDCA_withChi2; pt; type", fNPtBins, fPtBins, NbinsType, binsType);
  for(Int_t i=0;i<2;i++) { fOutList->Add(fDCA[i]); };
  for(Int_t i=0;i<2;i++) { fOutList->Add(fWithinDCA[i]); };
};
void AliEffFDContainer::AddCut(AliESDtrackCuts *inCut) {
  AliESDtrackCuts *tCut = (AliESDtrackCuts*)inCut->Clone();
  fChi2Cut = tCut->GetMaxChi2TPCConstrainedGlobal();
  if(fChi2Cut<1e9) {tCut->SetMaxChi2TPCConstrainedGlobal(); }; //Same comparison as in AliESDtrackCuts
  TString fDCAcut = tCut->GetMaxDCAToVertexXYPtDep();
  if(!fDCAcut.IsNull()) {
    fDCAcut.ReplaceAll("pt","x"); //Because in AliESDtrackCuts, it's pt, and not x
    tCut->SetMaxDCAToVertexXYPtDep("3*(pt^0)");
    if(!fDCAXYPtCut) fDCAXYPtCut = new TF1("DCAxy_pt_cut",fDCAcut.Data(),0.1,100);
  };
  fCutList->Add(tCut);
}
Double_t AliEffFDContainer::GetChi2TPCConstrained(const AliESDtrack *l_Tr) {
  // get vertex
  const AliESDVertex* vertex = 0;
  Double_t chi2TPCConstrainedVsGlobal=-2; //Initial value
  vertex = (AliESDVertex*)flESDEvent->GetPrimaryVertexTracks();
  if (!vertex || !vertex->GetStatus()) vertex = (AliESDVertex*)flESDEvent->GetPrimaryVertexSPD();
  //TPC vertex not considered
  // if ((!vertex || !vertex->GetStatus()) && fCutMaxChi2TPCConstrainedVsGlobalVertexType & kVertexTPC) vertex = esdEvent->GetPrimaryVertexTPC();
  if (vertex->GetStatus()) chi2TPCConstrainedVsGlobal = l_Tr->GetChi2TPCConstrainedVsGlobal(vertex);
  return chi2TPCConstrainedVsGlobal;
}
