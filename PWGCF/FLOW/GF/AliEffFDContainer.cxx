/*
  Container to calculate charged-particle efficiencies with PCC and also record DCAxy distributions for feed-down estimation.
  The PCC framework (AliMCSpectraWeights) used here is written by Patrick Huhn.
  For postprocessing of output to get efficiencies and feed-down corrections, use the
  postprocessing macro at: https://github.com/vvislavi/Feeddown
  If used, please acknowledge the authors of AliMCSpectraWeights as well as myself
  Author: Vytautas Vislavicius
*/
#include "AliEffFDContainer.h"
AliEffFDContainer::AliEffFDContainer():
  TNamed("",""),
  fOutList(0),
  fCutList(0),
  fChi2Cut(1e10),
  fDCAXYPtCut(0),
  fIsMC(kFALSE),
  fUseGenPt(kTRUE),
  fInitialized(kFALSE),
  flMCEvent(0),
  flESDEvent(0),
  flMultSel(0),
  flmcWeightsHandler(0),
  flMCSpectraWeights(0),
  fEff(0),
  fDCA(0),
  fWithinDCA(0),
  fIdentifier(new TNamed("Identifier","")),
  fPtBins(0),
  fNPtBins(0),
  fCentBins(0),
  fNCentBins(0),
  fCentEst("V0M"),
  fCent(-999),
  fPtMin(0.2),
  fPtMax(3.0),
  fEta(0.8),
  fEtaLow(-9999)
{
  fOutList = new TList();
  fOutList->SetOwner(kTRUE);
  fCutList = new TList();
  fCutList->SetOwner(kTRUE);
};
AliEffFDContainer::AliEffFDContainer(TString lName, TString lTitle, Bool_t lIsMC):
  TNamed(lName,lTitle),
  fOutList(0),
  fCutList(0),
  fChi2Cut(1e10),
  fDCAXYPtCut(0),
  fIsMC(lIsMC),
  fUseGenPt(kTRUE),
  fInitialized(kFALSE),
  flMCEvent(0),
  flESDEvent(0),
  flMultSel(0),
  flmcWeightsHandler(0),
  flMCSpectraWeights(0),
  fEff(0),
  fDCA(0),
  fWithinDCA(0),
  fIdentifier(new TNamed("Identifier","")),
  fPtBins(0),
  fNPtBins(0),
  fCentBins(0),
  fNCentBins(0),
  fCentEst("V0M"),
  fCent(-999),
  fPtMin(0.2),
  fPtMax(3.0),
  fEta(0.8),
  fEtaLow(-9999)
{
  fOutList = new TList();
  fOutList->SetOwner(kTRUE);
  fOutList->SetName(lName.Data());
  fCutList = new TList();
  fCutList->SetOwner(kTRUE);
  fCutList->SetName(Form("%s_Cuts",lName.Data()));
};
AliEffFDContainer::~AliEffFDContainer() {
  delete fOutList;
  delete fCutList;
}
Long64_t AliEffFDContainer::Merge(TCollection *collist) {
  Long64_t nmerged=0;
  AliEffFDContainer *l_EFFD = 0;
  TIter all_EFFD(collist);
  while ((l_EFFD = ((AliEffFDContainer*) all_EFFD()))) {
    TList *tarL = l_EFFD->fOutList;
    if(!tarL) continue;
    if(!fOutList) {
      fOutList = (TList*)tarL->Clone();
    } else
      for(Int_t i=0; i<fOutList->GetEntries(); i++) ((TH1*)fOutList->At(i))->Add((TH1*)tarL->At(i));
    nmerged++;
  };
  return nmerged;
};

void AliEffFDContainer::NewEvent(AliESDEvent &inputESD) { //Data part
  flESDEvent = &inputESD;
  flMultSel= dynamic_cast<AliMultSelection*>(flESDEvent->FindListObject("MultSelection"));
  fCent = flMultSel->GetMultiplicityPercentile(fCentEst.Data());
};
void AliEffFDContainer::NewEvent(AliMCEvent &inputMC) { //MC part
  flMCEvent = &inputMC;
  flmcWeightsHandler = static_cast<AliMCSpectraWeightsHandler*>(flESDEvent->FindListObject("fMCSpectraWeights"));
  flMCSpectraWeights = (flmcWeightsHandler) ? flmcWeightsHandler->fMCSpectraWeight : nullptr;
};
void AliEffFDContainer::Fill(AliESDEvent &inputESD, AliMCEvent &inputMC) {
  if(!fIsMC) {
    printf("\n\n\n");
    printf("Hi! I see you've called AliEddDFContainer::Fill(...) for MC event while the container was set up for data! You probably forgot to set the correct MC flag in the constructor.\n");
    printf("I would love to fix this for you, but this would create more problems in the output. Unfortunatelly, I will have to crash now...\n");
    AliFatal("Please set the correct MC flag in the constructor or use the appropriate fill method!\n");
  };
  if(!fInitialized) CreateHistograms();
  NewEvent(inputESD);
  NewEvent(inputMC); //MC event _has_ to go after data setup b/c of the weight handler setup. Little bit inconvenient, in case the fIsMC is not setup correct.
  // printf("Address comparison: %p vs %p\n",flESDEvent,(void*)&inputESD);
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
    if(lPart->Charge()==0.) continue;
    eta = lPart->Eta();
    if(!CheckEta(eta)) continue;
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
    if(!CheckEta(eta)) continue;
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
    if(lPart->Charge()==0.) continue;
    eta = lPart->Eta();
    if(!CheckEta(eta)) continue;
    if(fUseGenPt) pt = lPart->Pt();
    if(pt<fPtMin || pt>fPtMax) continue;
    CompWeight = flMCSpectraWeights->GetMCSpectraWeightNominal(lPart->Particle());
    Double_t secWeight = flMCSpectraWeights->GetWeightForSecondaryParticle(lPart->Particle());
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
      if(passChi2) fDCA[1]->Fill(pt,1.,dcaXY,secWeight);
      if(passDCA) {
        fWithinDCA[0]->Fill(pt,1.,secWeight);
        if(passChi2) fWithinDCA[1]->Fill(pt,1.,secWeight);
      };
    } else if(lPart->IsSecondaryFromMaterial()) {
      fDCA[0]->Fill(pt,2.,dcaXY,1);
      if(passChi2) fDCA[1]->Fill(pt,2.,dcaXY,1);
      if(passDCA) {
        fWithinDCA[0]->Fill(pt,2.,1);
        if(passChi2) fWithinDCA[1]->Fill(pt,2.,1);
      };
    };
  };
};
void AliEffFDContainer::Fill(AliESDEvent &inputESD) {
  if(!fInitialized) CreateHistograms();
  NewEvent(inputESD);
  // printf("Address comparison: %p vs %p\n",flESDEvent,(void*)&inputESD);
  if(!flESDEvent) {printf("ESD event not set! Not filling...\n"); return; };
  Int_t nTracks = flESDEvent->GetNumberOfTracks();
  AliESDtrack *lTrack;
  Double_t pt, eta;
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
    if(!CheckEta(eta)) continue;
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
    //Second argument could be centrality here
    fDCA[0]->Fill(pt,fCent,dcaXY);
    if(passChi2) fDCA[1]->Fill(pt,fCent,dcaXY);
    if(passDCA) {
      fWithinDCA[0]->Fill(pt,fCent);
      if(passChi2) {
        fWithinDCA[1]->Fill(pt,fCent);
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
void AliEffFDContainer::CreateHistograms(Bool_t forceRecreate) {
  if(forceRecreate) {
    if(fDCA) { delete [] fDCA; fDCA = 0; };
    if(fWithinDCA) { delete [] fWithinDCA; fWithinDCA = 0; };
    if(fEff) { delete [] fEff; fEff = 0; };
  }
  if(!fCentBins) {
    Double_t centBins[] = {0,100};
    SetCentralityBins(1,centBins);
  };
  if(!fPtBins) {
    Double_t ptBins[] = {0,1,2,3,4,5,6,7,8,9,10};
    SetPtBins(10,ptBins);
  }
  //Setting up DCAs
  Double_t binsDCA[61] = {-3.00, -2.90, -2.80, -2.70, -2.60, -2.50, -2.40, -2.30, -2.20, -2.10, -2.00, -1.90, -1.80, -1.70, -1.60, -1.50, -1.40, -1.30, -1.20, -1.10, -1.00, -0.90, -0.80, -0.70, -0.60, -0.50, -0.40, -0.30, -0.20, -0.10, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00, 1.10, 1.20, 1.30, 1.40, 1.50, 1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.70, 2.80, 2.90, 3.00};
  Int_t NbinsDCA = 60;
  Double_t binsDCAFine[301] = {-1.50, -1.49, -1.48, -1.47, -1.46, -1.45, -1.44, -1.43, -1.42, -1.41, -1.40, -1.39, -1.38, -1.37, -1.36, -1.35, -1.34, -1.33, -1.32, -1.31, -1.30, -1.29, -1.28, -1.27, -1.26, -1.25, -1.24, -1.23, -1.22, -1.21, -1.20, -1.19, -1.18, -1.17, -1.16, -1.15, -1.14, -1.13, -1.12, -1.11, -1.10, -1.09, -1.08, -1.07, -1.06, -1.05, -1.04, -1.03, -1.02, -1.01, -1.00, -0.99, -0.98, -0.97, -0.96, -0.95, -0.94, -0.93, -0.92, -0.91, -0.90, -0.89, -0.88, -0.87, -0.86, -0.85, -0.84, -0.83, -0.82, -0.81, -0.80, -0.79, -0.78, -0.77, -0.76, -0.75, -0.74, -0.73, -0.72, -0.71, -0.70, -0.69, -0.68, -0.67, -0.66, -0.65, -0.64, -0.63, -0.62, -0.61, -0.60, -0.59, -0.58, -0.57, -0.56, -0.55, -0.54, -0.53, -0.52, -0.51, -0.50, -0.49, -0.48, -0.47, -0.46, -0.45, -0.44, -0.43, -0.42, -0.41, -0.40, -0.39, -0.38, -0.37, -0.36, -0.35, -0.34, -0.33, -0.32, -0.31, -0.30, -0.29, -0.28, -0.27, -0.26, -0.25, -0.24, -0.23, -0.22, -0.21, -0.20, -0.19, -0.18, -0.17, -0.16, -0.15, -0.14, -0.13, -0.12, -0.11, -0.10, -0.09, -0.08, -0.07, -0.06, -0.05, -0.04, -0.03, -0.02, -0.01, -0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25, 0.26, 0.27, 0.28, 0.29, 0.30, 0.31, 0.32, 0.33, 0.34, 0.35, 0.36, 0.37, 0.38, 0.39, 0.40, 0.41, 0.42, 0.43, 0.44, 0.45, 0.46, 0.47, 0.48, 0.49, 0.50, 0.51, 0.52, 0.53, 0.54, 0.55, 0.56, 0.57, 0.58, 0.59, 0.60, 0.61, 0.62, 0.63, 0.64, 0.65, 0.66, 0.67, 0.68, 0.69, 0.70, 0.71, 0.72, 0.73, 0.74, 0.75, 0.76, 0.77, 0.78, 0.79, 0.80, 0.81, 0.82, 0.83, 0.84, 0.85, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91, 0.92, 0.93, 0.94, 0.95, 0.96, 0.97, 0.98, 0.99, 1.00, 1.01, 1.02, 1.03, 1.04, 1.05, 1.06, 1.07, 1.08, 1.09, 1.10, 1.11, 1.12, 1.13, 1.14, 1.15, 1.16, 1.17, 1.18, 1.19, 1.20, 1.21, 1.22, 1.23, 1.24, 1.25, 1.26, 1.27, 1.28, 1.29, 1.30, 1.31, 1.32, 1.33, 1.34, 1.35, 1.36, 1.37, 1.38, 1.39, 1.40, 1.41, 1.42, 1.43, 1.44, 1.45, 1.46, 1.47, 1.48, 1.49, 1.50};
  Int_t NbinsDCAFine = 300;
  Double_t binsType[4] = {-0.5,0.5,1.5,2.5};
  Int_t NbinsType = 3;
  Double_t *lYAxis = fIsMC?binsType:fCentBins;
  Int_t     nYAxis = fIsMC?NbinsType:fNCentBins;
  fDCA = new TH3D*[2];
  fWithinDCA = new TH2D*[2];
  fDCA[0] = new TH3D(makeName("DCAxy_noChi2"),"DCAxy_noChi2; pt; type; dcaxy", fNPtBins, fPtBins, nYAxis, lYAxis, NbinsDCA, binsDCA);
  fDCA[1] = new TH3D(makeName("DCAxy_withChi2"),"DCAxy_withChi2; pt; type; dcaxy", fNPtBins, fPtBins, nYAxis, lYAxis, NbinsDCAFine, binsDCAFine);
  fWithinDCA[0] = new TH2D(makeName("WithinDCA_noChi2"),"WithinDCA_noChi2; pt; type", fNPtBins, fPtBins, nYAxis, lYAxis);
  fWithinDCA[1] = new TH2D(makeName("WithinDCA_withChi2"),"WithinDCA_withChi2; pt; type", fNPtBins, fPtBins, nYAxis, lYAxis);
  for(Int_t i=0;i<2;i++) { fDCA[i]->SetDirectory(0); fOutList->Add(fDCA[i]); };
  for(Int_t i=0;i<2;i++) { fWithinDCA[i]->SetDirectory(0); fOutList->Add(fWithinDCA[i]); };
  fInitialized=kTRUE;
  if(!fIsMC) return; //If running on data, no need to fill in efficiencies
  //Setting up efficiencies
  fEff = new TH2D*[4];
  fEff[0] = new TH2D(makeName("nChGen_Weighted"),"ChGen_Weighted; #it{p}_{T} (GeV/#it{c}); Cent.",fNPtBins,fPtBins,fNCentBins,fCentBins);
  fEff[1] = new TH2D(makeName("nChRec_Weighted"),"ChGen_Weighted; #it{p}_{T} (GeV/#it{c}); Cent.",fNPtBins,fPtBins,fNCentBins,fCentBins);
  fEff[2] = new TH2D(makeName("nChGen_Uneighted"),"ChGen_Weighted; #it{p}_{T} (GeV/#it{c}); Cent.",fNPtBins,fPtBins,fNCentBins,fCentBins);
  fEff[3] = new TH2D(makeName("nChRec_Uneighted"),"ChGen_Weighted; #it{p}_{T} (GeV/#it{c}); Cent.",fNPtBins,fPtBins,fNCentBins,fCentBins);
  for(Int_t i=0;i<4;i++) { fEff[i]->SetDirectory(0); fOutList->Add(fEff[i]); };


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
void AliEffFDContainer::AddCut(Int_t fFilterBit) {
  AliESDtrackCuts *lTC;
  Int_t nTotFlag=0;
  if(fFilterBit & 32) { //FB32
    lTC = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
    AddCut(lTC);
    delete lTC;
    nTotFlag+=32;
    printf("Adding track cuts corresponding to FB32...\n");
  };
  if(fFilterBit & 64) {//FB64
    lTC = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
    lTC->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kNone);
    lTC->SetClusterRequirementITS(AliESDtrackCuts::kSDD,AliESDtrackCuts::kFirst);
    AddCut(lTC);
    delete lTC;
    nTotFlag+=64;
    printf("Adding track cuts corresponding to FB64...\n");
  };
  if(fFilterBit & 256) {//FB256
    lTC = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
    lTC->SetMaxDCAToVertexXY(2.4);
    lTC->SetMaxDCAToVertexZ(3.2);
    lTC->SetDCAToVertex2D(kTRUE);
    lTC->SetMaxChi2TPCConstrainedGlobal(36);
    lTC->SetMaxFractionSharedTPCClusters(0.4);
    AddCut(lTC);
    delete lTC;
    nTotFlag+=256;
    printf("Adding track cuts corresponding to FB256...\n");
  };
  if(fFilterBit & 512) {//FB256
    lTC = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
    lTC->SetMaxDCAToVertexXY(2.4);
    lTC->SetMaxDCAToVertexZ(3.2);
    lTC->SetDCAToVertex2D(kTRUE);
    lTC->SetMaxChi2TPCConstrainedGlobal(36);
    lTC->SetMaxFractionSharedTPCClusters(0.4);
    lTC->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
    lTC->SetRequireITSRefit(kTRUE);
    AddCut(lTC);
    delete lTC;
    nTotFlag+=512;
    printf("Adding track cuts corresponding to FB512...\n");
  };
  if(nTotFlag!=fFilterBit) {
    printf("Warning!!! I could not add the requested (full) filter bit; instead, I've only added FB %i. Please construct the relevant remaining part as AliESDtrackCuts instead...\n", nTotFlag);
  }
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
Bool_t AliEffFDContainer::AddContainer(AliEffFDContainer *target) {
  TList *tlist = target->GetOutList();
  if(tlist->GetEntries() != fOutList->GetEntries()) { printf("Containers have different number of contents!\n"); return 0;};
  for(Int_t i=0;i<tlist->GetEntries();i++) {
    TH1 *h1 = (TH1*)fOutList->At(i);
    TH1 *h2 = (TH1*)tlist->At(i);
    h1->Add(h2);
  }
  return 1;
}
