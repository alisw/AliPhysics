
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
  fEvNomFlag(1),
  fTrNomFlag(1),
  fAddPID(kFALSE),
  fNSpecies(1),
  fPIDResponse(0),
  fBayesPID(0),
  fMinBayesProb({0.95,0.85,0.85}),
  fIsMC(kFALSE),
  fUseGenPt(kTRUE),
  fInitialized(kFALSE),
  flMCEvent(0),
  flESDEvent(0),
  flAODEvent(0),
  flMultSel(0),
  flmcWeightsHandler(0),
  flMCSpectraWeights(0),
  fMCSWeights(0),
  fEff(0),
  fDCA(0),
  fWithinDCA(0),
  fPurity(0),
  fFDvsPhi(0),
  fIdentifier(new TNamed("Identifier","")),
  fSpNames({"","pi_","ka_","pr_"}),
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
  fEvNomFlag(1),
  fTrNomFlag(1),
  fAddPID(kFALSE),
  fNSpecies(1),
  fPIDResponse(0),
  fBayesPID(0),
  fMinBayesProb({0.95,0.85,0.85}),
  fIsMC(lIsMC),
  fUseGenPt(kTRUE),
  fInitialized(kFALSE),
  flMCEvent(0),
  flESDEvent(0),
  flAODEvent(0),
  flMultSel(0),
  flmcWeightsHandler(0),
  flMCSpectraWeights(0),
  fMCSWeights(0),
  fEff(0),
  fDCA(0),
  fWithinDCA(0),
  fPurity(0),
  fFDvsPhi(0),
  fIdentifier(new TNamed("Identifier","")),
  fSpNames({"","pi_","ka_","pr_"}),
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
void AliEffFDContainer::NewAODEvent(AliAODEvent &inputAOD, AliMCEvent &inputMC) { //AOD part
  flAODEvent = &inputAOD;
  flMCEvent = &inputMC;
  flMultSel= dynamic_cast<AliMultSelection*>(flAODEvent->FindListObject("MultSelection"));
  fCent = flMultSel->GetMultiplicityPercentile(fCentEst.Data());
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
  if(fCent<fCentBins[0] || fCent>fCentBins[fNCentBins]) return;
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
    fEff[0][0]->Fill(pt,fCent,CompWeight);
    fEff[0][2]->Fill(pt,fCent);
    if(fAddPID) {
      Int_t pidIndex = GetTruePIDIndex(TMath::Abs(lPart->PdgCode()))+1;
      if(pidIndex) {
        fEff[pidIndex][0]->Fill(pt,fCent,CompWeight);
        fEff[pidIndex][2]->Fill(pt,fCent);
      }
    }
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
    // eta = lPart->Eta(); ///Eta check either on gen. or rec., not on both
    // if(!CheckEta(eta)) continue; //Eta check either on gen. or rec., not on both
    if(fUseGenPt) pt = lPart->Pt();
    if(pt<fPtMin || pt>fPtMax) continue;
    CompWeight = flMCSpectraWeights->GetMCSpectraWeightNominal(lPart->Particle());
    Double_t secWeight = flMCSpectraWeights->GetWeightForSecondaryParticle(index);
    Int_t lBayesPIDIndex=0;
    Int_t lTruePIDIndex =0;
    Bool_t IndexMatch=kFALSE;
    Bool_t fillPIDHists=kFALSE;
    if(fAddPID) {
      lBayesPIDIndex = GetBayesPIDIndex(lTrack)+1;
      lTruePIDIndex  = GetTruePIDIndex(TMath::Abs(lPart->PdgCode()))+1;
      IndexMatch = (lBayesPIDIndex == lTruePIDIndex);
      fillPIDHists = (lBayesPIDIndex > 0) && IndexMatch;
    };
    if(lBayesPIDIndex && passDCA && passChi2) {
      fPurity[lBayesPIDIndex-1][0]->Fill(pt,fCent,CompWeight);
      fPurity[lBayesPIDIndex-1][1]->Fill(pt,fCent);
      if(lPart->IsPhysicalPrimary()) {
        fPurity[lBayesPIDIndex-1][2]->Fill(pt,fCent,CompWeight);
        fPurity[lBayesPIDIndex-1][3]->Fill(pt,fCent);
      }
    };
    if(lPart->IsPhysicalPrimary()) {
      fDCA[0][0]->Fill(pt,0.,dcaXY,CompWeight);
      if(passChi2) fDCA[0][1]->Fill(pt,0.,dcaXY,CompWeight);
      if(passDCA) {
        fWithinDCA[0][0]->Fill(pt,0.,fCent,CompWeight);
        if(passChi2) {
          fWithinDCA[0][1]->Fill(pt,0.,fCent,CompWeight);
          fEff[0][1]->Fill(pt,fCent,CompWeight);
          fEff[0][3]->Fill(pt,fCent);
        };
      };
      //PID part
      if(fillPIDHists) {
        fDCA[lBayesPIDIndex][0]->Fill(pt,0.,dcaXY,CompWeight);
        if(passChi2) fDCA[lBayesPIDIndex][1]->Fill(pt,0.,dcaXY,CompWeight);
        if(passDCA) {
          fWithinDCA[lBayesPIDIndex][0]->Fill(pt,0.,fCent,CompWeight);
          if(passChi2) {
            fWithinDCA[lBayesPIDIndex][1]->Fill(pt,0.,fCent,CompWeight);
            fEff[lBayesPIDIndex][1]->Fill(pt,fCent,CompWeight);
            fEff[lBayesPIDIndex][3]->Fill(pt,fCent);
          };
        };
      };
    } else if(lPart->IsSecondaryFromWeakDecay()) {
      fDCA[0][0]->Fill(pt,1.,dcaXY,CompWeight);
      if(passChi2) fDCA[0][1]->Fill(pt,1.,dcaXY,secWeight);
      if(passDCA) {
        fWithinDCA[0][0]->Fill(pt,1.,fCent,secWeight);
        if(passChi2) fWithinDCA[0][1]->Fill(pt,1.,fCent,secWeight);
      };
      //PID part:
      if(fillPIDHists) {
        fDCA[lBayesPIDIndex][0]->Fill(pt,1.,dcaXY,CompWeight);
        if(passChi2) fDCA[lBayesPIDIndex][1]->Fill(pt,1.,dcaXY,secWeight);
        if(passDCA) {
          fWithinDCA[lBayesPIDIndex][0]->Fill(pt,1.,fCent,secWeight);
          if(passChi2) fWithinDCA[lBayesPIDIndex][1]->Fill(pt,1.,fCent,secWeight);
        };
      };
    } else if(lPart->IsSecondaryFromMaterial()) {
      fDCA[0][0]->Fill(pt,2.,dcaXY,1);
      if(passChi2) fDCA[0][1]->Fill(pt,2.,dcaXY,1);
      if(passDCA) {
        fWithinDCA[0][0]->Fill(pt,2.,fCent,1);
        if(passChi2) fWithinDCA[0][1]->Fill(pt,2.,fCent,1);
      };
      //PID part:
      if(fillPIDHists) {
        fDCA[lBayesPIDIndex][0]->Fill(pt,2.,dcaXY,1);
        if(passChi2) fDCA[lBayesPIDIndex][1]->Fill(pt,2.,dcaXY,1);
        if(passDCA) {
          fWithinDCA[lBayesPIDIndex][0]->Fill(pt,2.,fCent,1);
          if(passChi2) fWithinDCA[lBayesPIDIndex][1]->Fill(pt,2.,fCent,1);
        };
      }
    };
  };
};
void AliEffFDContainer::Fill(AliAODEvent &inputAOD, AliMCEvent &inputMC) {
  if(!fIsMC) {
    printf("\n\n\n");
    printf("Hi! I see you've called AliEddDFContainer::Fill(...) for MC event while the container was set up for data! You probably forgot to set the correct MC flag in the constructor.\n");
    printf("I would love to fix this for you, but this would create more problems in the output. Unfortunatelly, I will have to crash now...\n");
    AliFatal("Please set the correct MC flag in the constructor or use the appropriate fill method!\n");
  };
  if(!fInitialized) CreateHistograms();
  NewAODEvent(inputAOD, inputMC);
  //For testing purposes, bypassing centrality check. Local files all have centrality 199. Need to double-check why
  if(fCent<fCentBins[0] || fCent>fCentBins[fNCentBins]) return;

  if(!flAODEvent) {printf("AOD event not set! Not filling...\n"); return; };
  AliGFWFlags *lFlags = (AliGFWFlags*)flAODEvent->FindListObject("GFWFlags");
  if(!lFlags) {printf("GFWFlags were not found!\n"); return; };
  UInt_t gEventFlag = lFlags->GetEventFlags();
  if(!(gEventFlag&fEvNomFlag)) return; //If not the selected event flag, then move on

  Int_t nPrimPart = flMCEvent->GetNumberOfTracks();
  Int_t nTracks = flAODEvent->GetNumberOfTracks();
  AliMCParticle *lPart;
  AliAODTrack *lTrack;
  Double_t pt, eta;
  Double_t CompWeight;
  Int_t nMulti = CalculateMult();
  Int_t l_MultiBin=fMCSWeights->GetZaxis()->FindBin(nMulti);
  if(l_MultiBin>fMCSWeights->GetNbinsZ()) l_MultiBin--;
  else if(l_MultiBin<1) l_MultiBin=1;
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
    Int_t primIndex = GetMCSWPrimIndex(lPart)+1;
    Int_t ptBin = fMCSWeights->GetYaxis()->FindBin(pt);
    if(ptBin>fMCSWeights->GetNbinsY()) ptBin--;
    else if(ptBin<1) ptBin=1;
    CompWeight = fMCSWeights->GetBinContent(primIndex,ptBin,l_MultiBin);//flMCSpectraWeights->GetMCSpectraWeightNominal(lPart->Particle());
    //Proceed with filling
    fEff[0][0]->Fill(pt,fCent,CompWeight);
    fEff[0][2]->Fill(pt,fCent);
    if(fAddPID) {
      Int_t pidIndex = GetTruePIDIndex(TMath::Abs(lPart->PdgCode()))+1;
      if(pidIndex) {
        fEff[pidIndex][0]->Fill(pt,fCent,CompWeight);
        fEff[pidIndex][2]->Fill(pt,fCent);
      }
    }
  };
  //Track loop
  for(Int_t iTr=0;iTr<lFlags->GetNFiltered();iTr++) {
    UInt_t gTrackFlags = lFlags->GetTrackFlag(iTr);
    if(!(gTrackFlags&fTrNomFlag)) continue; //Check if we want to accept the track
    Int_t trInd = lFlags->GetTrackIndex(iTr);
    lTrack = (AliAODTrack*)flAODEvent->GetTrack(trInd);
    pt  = lTrack->Pt();
    //Tracks already prefiltered with a given eta, but we shouldn't cut on both rec. & gen. eta! So here just an extra check of rec. eta, in case of running with multiple eta bins
    eta = lTrack->Eta();
    if(!CheckEta(eta)) continue;
    //Fetch the corresponding MC particle
    Int_t fLabel = lTrack->GetLabel();
    Int_t index = TMath::Abs(fLabel);
    if(AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(index, flMCEvent)) continue;
    if (index < 0) continue;
    lPart = (AliMCParticle*)flMCEvent->GetTrack(index);
    if(!lPart) continue;
    if(lPart->Charge()==0.) continue;
    //Weights -- should be fetched from histograms
    Int_t primIndex = GetMCSWPrimIndex(lPart)+1;
    Int_t ptBin = fMCSWeights->GetYaxis()->FindBin(pt);
    if(ptBin>fMCSWeights->GetNbinsY()) ptBin--;
    else if(ptBin<1) ptBin=1;
    Double_t CompWeight = fMCSWeights->GetBinContent(primIndex,ptBin,l_MultiBin);
    Double_t secWeight=1;
    if(lPart->IsSecondaryFromWeakDecay()) {
      Double_t ptMother;
      Int_t secIndex = GetMCSWMotherIndex(lPart,ptMother)+1;
      if(!secIndex) continue;
      Int_t ptMotherBin = fMCSWeights->GetYaxis()->FindBin(ptMother);
      if(ptMotherBin>fMCSWeights->GetNbinsY()) ptMotherBin--;
      else if(ptMotherBin<1) ptMotherBin=1;
      secWeight = fMCSWeights->GetBinContent(secIndex,ptMotherBin,l_MultiBin);
    }
    if(fUseGenPt) pt = lPart->Pt();
    //PID part:
    Int_t lBayesPIDIndex=0, lTruePIDIndex=0;
    Bool_t IndexMatch=kFALSE, fillPIDHists=kFALSE;
    if(fAddPID) {
      lBayesPIDIndex = GetBayesPIDIndex(lTrack)+1;
      lTruePIDIndex  = GetTruePIDIndex(TMath::Abs(lPart->PdgCode()))+1;
      IndexMatch = (lBayesPIDIndex == lTruePIDIndex);
      fillPIDHists = (lBayesPIDIndex > 0) && IndexMatch;
    };
    if(lBayesPIDIndex) {
      fPurity[lBayesPIDIndex-1][0]->Fill(pt,fCent,CompWeight);
      fPurity[lBayesPIDIndex-1][1]->Fill(pt,fCent);
      if(lPart->IsPhysicalPrimary()) {
        fPurity[lBayesPIDIndex-1][2]->Fill(pt,fCent,CompWeight);
        fPurity[lBayesPIDIndex-1][3]->Fill(pt,fCent);
      }
    };
    //Fetch phi of the track and recenter:
    Double_t l_phi = lTrack->Phi();
    if(l_phi<0) l_phi+=TMath::TwoPi();
    //Start filling:
    //Filling DCA distributions doesn't make much sense
    //Also, everything passes the chi2 here (AOD tracks)
    if(lPart->IsPhysicalPrimary()) {
      fWithinDCA[0][0]->Fill(pt,0.,fCent,CompWeight);
      fWithinDCA[0][1]->Fill(pt,0.,fCent,1.);
      fEff[0][1]->Fill(pt,fCent,CompWeight);
      fEff[0][3]->Fill(pt,fCent);
      //Prim vs Phi. Primaries contribute to both, primaries and all
      fFDvsPhi[0][0]->Fill(pt,fCent,l_phi,CompWeight);
      fFDvsPhi[0][1]->Fill(pt,fCent,l_phi,CompWeight);
      //PID part
      if(fillPIDHists) {
        fWithinDCA[lBayesPIDIndex][0]->Fill(pt,0.,fCent,CompWeight);
        fWithinDCA[lBayesPIDIndex][1]->Fill(pt,0.,fCent,1.);
        fEff[lBayesPIDIndex][1]->Fill(pt,fCent,CompWeight);
        fEff[lBayesPIDIndex][3]->Fill(pt,fCent);
        //Prim vs Phi. Primaries contribute to both, primaries and all
        fFDvsPhi[lBayesPIDIndex][0]->Fill(pt,fCent,l_phi,CompWeight);
        fFDvsPhi[lBayesPIDIndex][1]->Fill(pt,fCent,l_phi,CompWeight);
      };
    } else if(lPart->IsSecondaryFromWeakDecay()) {
      fWithinDCA[0][0]->Fill(pt,1.,fCent,secWeight);
      fWithinDCA[0][1]->Fill(pt,1.,fCent,1.);
      //Prim vs Phi. Secondaries contribute to all
      fFDvsPhi[0][1]->Fill(pt,fCent,l_phi,CompWeight);
      //PID part:
      if(fillPIDHists) {
        fWithinDCA[lBayesPIDIndex][0]->Fill(pt,1.,fCent,secWeight);
        fWithinDCA[lBayesPIDIndex][1]->Fill(pt,1.,fCent,1);
        //Prim vs Phi. Secondaries contribute to all
        fFDvsPhi[lBayesPIDIndex][1]->Fill(pt,fCent,l_phi,CompWeight);

      }
    } else if(lPart->IsSecondaryFromMaterial()) {
      fWithinDCA[0][0]->Fill(pt,2.,fCent,1.);
      fWithinDCA[0][1]->Fill(pt,2.,fCent,1.);
      //Prim vs Phi. Secondaries contribute to all
      fFDvsPhi[0][1]->Fill(pt,fCent,l_phi,CompWeight);
      //PID part:
      if(fillPIDHists) {
        fWithinDCA[lBayesPIDIndex][0]->Fill(pt,2.,fCent,1);
        fWithinDCA[lBayesPIDIndex][1]->Fill(pt,2.,fCent,1);
        //Prim vs Phi. Secondaries contribute to all
        fFDvsPhi[lBayesPIDIndex][1]->Fill(pt,fCent,l_phi,CompWeight);
      }
    }
  }
};
void AliEffFDContainer::Fill(AliESDEvent &inputESD) {
  if(!fInitialized) CreateHistograms();
  NewEvent(inputESD);
  if(fCent<fCentBins[0] || fCent>fCentBins[fNCentBins]) return;
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
    fDCA[0][0]->Fill(pt,fCent,dcaXY);
    if(passChi2) fDCA[0][1]->Fill(pt,fCent,dcaXY);
    if(passDCA) {
      fWithinDCA[0][0]->Fill(pt,0.,fCent);
      if(passChi2) {
        fWithinDCA[0][1]->Fill(pt,0.,fCent);
      };
    };
    //PID part
    if(fAddPID) {
      Int_t lBayesPIDIndex = GetBayesPIDIndex(lTrack)+1;
      if(lBayesPIDIndex) {
        fDCA[lBayesPIDIndex][0]->Fill(pt,fCent,dcaXY);
        if(passChi2) fDCA[lBayesPIDIndex][1]->Fill(pt,fCent,dcaXY);
        if(passDCA) {
          fWithinDCA[lBayesPIDIndex][0]->Fill(pt,0.,fCent);
          if(passChi2) {
            fWithinDCA[lBayesPIDIndex][1]->Fill(pt,0.,fCent);
          };
        };
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
  //Phi bins for FD vs phi:
  const Int_t fNPhiBins=100;
  Double_t fPhiBins[fNPhiBins+1];
  for(Int_t i=0;i<=fNPhiBins;i++) fPhiBins[i] = TMath::TwoPi()/fNPhiBins * i;
  //Setting up y-axis: either centrality or prim/secondary time
  Double_t *lYAxis = fIsMC?binsType:fCentBins;
  Int_t     nYAxis = fIsMC?NbinsType:fNCentBins;
  if(fAddPID) fNSpecies=4; else fNSpecies=1;
  fDCA       = new TH3D**[fNSpecies];
  fWithinDCA = new TH3D**[fNSpecies];
  fEff       = new TH2D**[fNSpecies];
  fFDvsPhi   = new TH3D**[fNSpecies];
  if(fAddPID) fPurity = new TH2D**[fNSpecies-1]; //No need for purity for charged tracks
  for(Int_t iSpecie=0; iSpecie<fNSpecies; iSpecie++) {
    fDCA[iSpecie] = new TH3D*[2];
    fWithinDCA[iSpecie] = new TH3D*[2];
    fDCA[iSpecie][0] = new TH3D(makeName("DCAxy_noChi2",iSpecie),"DCAxy_noChi2; pt; type; dcaxy", fNPtBins, fPtBins, nYAxis, lYAxis, NbinsDCA, binsDCA);
    fDCA[iSpecie][1] = new TH3D(makeName("DCAxy_withChi2",iSpecie),"DCAxy_withChi2; pt; type; dcaxy", fNPtBins, fPtBins, nYAxis, lYAxis, NbinsDCAFine, binsDCAFine);
    fWithinDCA[iSpecie][0] = new TH3D(makeName("WithinDCA_noChi2",iSpecie),"WithinDCA_noChi2; pt; type; centrality", fNPtBins, fPtBins, NbinsType, binsType, fNCentBins, fCentBins);
    fWithinDCA[iSpecie][1] = new TH3D(makeName("WithinDCA_withChi2",iSpecie),"WithinDCA_withChi2; pt; type; centrality", fNPtBins, fPtBins, NbinsType, binsType, fNCentBins, fCentBins);
    for(Int_t i=0;i<2;i++) { fDCA[iSpecie][i]->SetDirectory(0); fOutList->Add(fDCA[iSpecie][i]); };
    for(Int_t i=0;i<2;i++) { fWithinDCA[iSpecie][i]->SetDirectory(0); fOutList->Add(fWithinDCA[iSpecie][i]); };
    if(!fIsMC) continue; //If running on data, no need to fill in efficiencies
    //Setting up efficiencies
    fEff[iSpecie] = new TH2D*[4];
    fEff[iSpecie][0] = new TH2D(makeName("nChGen_Weighted",iSpecie),"ChGen_Weighted; #it{p}_{T} (GeV/#it{c}); Cent.",fNPtBins,fPtBins,fNCentBins,fCentBins);
    fEff[iSpecie][1] = new TH2D(makeName("nChRec_Weighted",iSpecie),"ChGen_Weighted; #it{p}_{T} (GeV/#it{c}); Cent.",fNPtBins,fPtBins,fNCentBins,fCentBins);
    fEff[iSpecie][2] = new TH2D(makeName("nChGen_Uneighted",iSpecie),"ChGen_Weighted; #it{p}_{T} (GeV/#it{c}); Cent.",fNPtBins,fPtBins,fNCentBins,fCentBins);
    fEff[iSpecie][3] = new TH2D(makeName("nChRec_Uneighted",iSpecie),"ChGen_Weighted; #it{p}_{T} (GeV/#it{c}); Cent.",fNPtBins,fPtBins,fNCentBins,fCentBins);
    for(Int_t i=0;i<4;i++) { fEff[iSpecie][i]->SetDirectory(0); fOutList->Add(fEff[iSpecie][i]); };
    //Feeddown vs phi (in given eta bin)
    fFDvsPhi[iSpecie] = new TH3D*[2];
    fFDvsPhi[iSpecie][0] = new TH3D(makeName("PrimVsPhi",iSpecie),"PrimariesVsPhi; #it{p}_{T} (GeV/#it{c}); cent.; #varphi", fNPtBins, fPtBins, fNCentBins, fCentBins, fNPhiBins, fPhiBins);
    fFDvsPhi[iSpecie][1] = new TH3D(makeName("AllVsPhi",iSpecie), "AllVsPhi; #it{p}_{T} (GeV/#it{c}); cent.; #varphi", fNPtBins, fPtBins, fNCentBins, fCentBins, fNPhiBins, fPhiBins);
    for(Int_t i=0;i<2;i++) { fFDvsPhi[iSpecie][i]->SetDirectory(0); fOutList->Add(fFDvsPhi[iSpecie][i]); };
    if(fAddPID) {
      if(!iSpecie) continue; //Not adding it for charged. Then we have to be carefull when filling
      fPurity[iSpecie-1] = new TH2D*[4];
      fPurity[iSpecie-1][0] = new TH2D(makeName("PurityAll_Weighted",iSpecie),"SelectedByBayes; #it{p}_{T} (GeV/#it{c}); Cent.",fNPtBins,fPtBins,fNCentBins,fCentBins);
      fPurity[iSpecie-1][1] = new TH2D(makeName("PurityAll_Uneighted",iSpecie),"SelectedByBayes; #it{p}_{T} (GeV/#it{c}); Cent.",fNPtBins,fPtBins,fNCentBins,fCentBins);
      fPurity[iSpecie-1][2] = new TH2D(makeName("PurityPrimary_Weighted",iSpecie),"SelectedByBayes; #it{p}_{T} (GeV/#it{c}); Cent.",fNPtBins,fPtBins,fNCentBins,fCentBins);
      fPurity[iSpecie-1][3] = new TH2D(makeName("PurityPrimary_Uneighted",iSpecie),"SelectedByBayes; #it{p}_{T} (GeV/#it{c}); Cent.",fNPtBins,fPtBins,fNCentBins,fCentBins);
      for(Int_t i=0;i<4;i++) { fPurity[iSpecie-1][i]->SetDirectory(0); fOutList->Add(fPurity[iSpecie-1][i]); };
    };
  };
  fInitialized=kTRUE;
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
Int_t AliEffFDContainer::GetBayesPIDIndex(AliVTrack *l_track) {
  Double_t l_Probs[AliPID::kSPECIES];
  Bool_t l_TOFUsed = fBayesPID->ComputeProbabilities(l_track, fPIDResponse, l_Probs) & AliPIDResponse::kDetTOF;
  Int_t pidInd = 0;
  for(Int_t i=0;i<AliPID::kSPECIES; i++) pidInd=(l_Probs[i]>l_Probs[pidInd])?i:pidInd;
  Int_t retInd = pidInd-AliPID::kPion; //Not interested in e+mu, so realign to 0
  if(retInd<0 || retInd>2) return -1; //Shouldn't be larger than 2, but just to be safe
  if(l_Probs[pidInd] < fMinBayesProb[retInd]) return -1;
  //check nsigma cuts
  if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(l_track,(AliPID::EParticleType)pidInd))>3) return -1;
  if(l_TOFUsed) if(TMath::Abs(fPIDResponse->NumberOfSigmasTOF(l_track,(AliPID::EParticleType)pidInd))>3) return -1;
  return retInd;
}
Int_t AliEffFDContainer::GetTruePIDIndex(const Int_t &pdgcode) {
  if(pdgcode==211) return 0;
  else if(pdgcode==321) return 1;
  else if(pdgcode==2212) return 2;
  else return -1;
};
Int_t AliEffFDContainer::GetMCSWPrimIndex(AliMCParticle *part) {
  Int_t ipdg = TMath::Abs(part->PdgCode());
  if (ipdg == 211) return AliMCSpectraWeights::ParticleType::kPion;
  if (ipdg == 321) return AliMCSpectraWeights::ParticleType::kKaon;
  if (ipdg == 2212) return AliMCSpectraWeights::ParticleType::kProtons;
  if (ipdg == 3222) return AliMCSpectraWeights::ParticleType::kSigmaPlus;
  if (ipdg == 3112) return AliMCSpectraWeights::ParticleType::kSigmaMinus;
  if(ipdg == 3122) return AliMCSpectraWeights::ParticleType::kLambda;
  return AliMCSpectraWeights::ParticleType::kRest;
};
Int_t AliEffFDContainer::GetMCSWMotherIndex(AliMCParticle *part, Double_t &ptMother) {
  auto const absPDG = TMath::Abs(part->PdgCode());
  auto motherPartLabel = part->GetMother();
  if (motherPartLabel<0) return -1;
  auto const motherPart = (AliMCParticle*)flMCEvent->GetTrack(motherPartLabel);
  if(!motherPart) return -1;
  ptMother = motherPart->Pt();
  auto const motherPDG = TMath::Abs(motherPart->PdgCode());
  //Lambda case
  if ((motherPDG == 3122 || motherPDG == 3222 || motherPDG == 3112 || motherPDG == 3212) && motherPart->IsPhysicalPrimary()) return AliMCSpectraWeights::ParticleType::kSigmaPlus;
  //K0short case
  if((motherPDG == 310 || motherPDG == 130 || motherPDG == 311 || motherPDG == 321) && motherPart->IsPhysicalPrimary()) return AliMCSpectraWeights::ParticleType::kKaon;
  //secondary from pion
  if( motherPDG == 211 && motherPart->IsPhysicalPrimary()) return AliMCSpectraWeights::ParticleType::kPion;
  //Xi->lambda->proton
  if((motherPDG == 3122 && TMath::Abs(flMCEvent->MotherOfParticle(motherPartLabel)->PdgCode()) == 3312)) return AliMCSpectraWeights::ParticleType::kSigmaPlus;
  //Otherwise
  return -1;
};
Int_t AliEffFDContainer::CalculateMult() {
  if(!flMCEvent) {printf("MC event not set!\n"); return -1; };
  Int_t retCount=0;
  for(Int_t i=0;i<flMCEvent->GetNumberOfTracks();i++) {
    AliMCParticle *mcp = dynamic_cast<AliMCParticle*>(flMCEvent->GetTrack(i));
    if(!mcp) continue;
    if(!mcp->IsPhysicalPrimary()) continue;
    if(TMath::Abs(mcp->Charge()) < 0.01) continue;
    if(TMath::Abs(mcp->Eta())>0.5) continue;
    if(mcp->Pt()<0.05) continue;
    retCount++;
  };
  return retCount;
};
TH2 *AliEffFDContainer::get2DRatio(TString numID, TString denID, Int_t iSpecie) {
  TH2 *hNum = (TH2*)fetchObj(numID,iSpecie);
  if(!hNum) {printf("Could not find %s in %s!\n",makeName(numID,iSpecie).Data(), this->GetName()); return 0; };
  TH2 *hDen = (TH2*)fetchObj(denID,iSpecie);
  if(!hDen) {printf("Could not find %s in %s!\n",makeName(denID,iSpecie).Data(), this->GetName()); return 0; };
  hNum = (TH2*)hNum->Clone(denID.Data());
  hNum->SetDirectory(0);
  hNum->Divide(hDen);
  return hNum;
};
TH1 *AliEffFDContainer::get1DRatio(TString numID, TString denID, Int_t iSpecie, Int_t yb1, Int_t yb2) {
  TH2 *hNum = (TH2*)fetchObj(numID,iSpecie);
  if(!hNum) {printf("Could not find %s in %s!\n",makeName(numID,iSpecie).Data(), this->GetName()); return 0; };
  TH2 *hDen = (TH2*)fetchObj(denID,iSpecie);
  if(!hDen) {printf("Could not find %s in %s!\n",makeName(denID,iSpecie).Data(), this->GetName()); return 0; };
  if(yb1<1) yb1=1;
  if(yb2<yb1 || yb2>hNum->GetNbinsY()) yb2 = hNum->GetNbinsY();
  TH1 *h1Num = hNum->ProjectionX(numID.Data(),yb1,yb2);
  TH1 *h1Den = hDen->ProjectionX(denID.Data(),yb1,yb2);
  h1Num->SetDirectory(0);
  h1Num->Divide(h1Den);
  delete h1Den;
  return h1Num;
};
TH1 *AliEffFDContainer::getPureFeeddown(Int_t iSpecie, Int_t centBin1, Int_t centBin2) {
  TObject *obj = fetchObj("WithinDCA_withChi2",iSpecie);
  //Older versions have TH2 instead of TH3 (that is, no centrality. Don't want to break it.)
  TH2 *l_h2 = dynamic_cast<TH2*>(obj);
  if(!l_h2) {
    TH3 *l_h3 = dynamic_cast<TH3*>(obj);
    if(!l_h3) {
      printf("Could not find %s in %s!\n",makeName("WithinDCA_withChi2",iSpecie).Data(), this->GetName());
      return 0;
    };
    if(centBin1<1) centBin1 = 1;
    if(centBin2>l_h3->GetNbinsZ() || centBin2<1) centBin2 = l_h3->GetNbinsZ();
    l_h3->GetZaxis()->SetRange(centBin1,centBin2);
    l_h2 = (TH2*)l_h3->Project3D("yx");
    l_h3->GetZaxis()->SetRange(1,l_h3->GetNbinsZ());
  };
  return getPureFeeddown(iSpecie,l_h2);
}
TH1 *AliEffFDContainer::getPureFeeddown(Int_t iSpecie, TH2 *inh) {
  TH1 *hPrim = (TH1*)inh->ProjectionX(makeName("PrimOverAll",iSpecie),1,1);
  hPrim->SetDirectory(0);
  TH1 *hAll  = (TH1*)inh->ProjectionX("All",1,inh->GetNbinsX());
  hPrim->Divide(hAll);
  delete hAll;
  return hPrim;
}
TH2 *AliEffFDContainer::getFDvsPhi(Int_t iSpecie, Bool_t RatioToIntegrated, Int_t cBin1, Int_t cBin2) {
  TH3D *hPrim = (TH3D*)fetchObj("PrimVsPhi",iSpecie);
  if(cBin1<1) cBin1=1;
  if(cBin2<cBin1 || cBin2<1 || cBin2>hPrim->GetNbinsY()) cBin2=hPrim->GetNbinsY();
  TString hName("FDvsPhi");
  if(cBin1==1 && cBin2==hPrim->GetNbinsY()) hName.Append("_MB");
  else hName.Append(Form("_CentBin_%i_%i",cBin1,cBin2));
  TH3D *hAll = (TH3D*)fetchObj("AllVsPhi",iSpecie);
  hPrim->GetYaxis()->SetRange(cBin1,cBin2);
  hAll ->GetYaxis()->SetRange(cBin1,cBin2);
  TH2D *h2Prim = (TH2D*)hPrim->Project3D("zx");
  TH2D *h2All  = (TH2D*)hAll ->Project3D("zx");
  h2Prim->SetDirectory(0);
  h2Prim->SetName(makeName(hName,iSpecie).Data());
  // h2Prim->Divide(h2Prim,h2All,1,1,"B"); //if we want binomial error propagation
  h2Prim->Divide(h2All);
  delete h2All;
  if(RatioToIntegrated) { //Normalize by central value, if required
    TH1D *h1Prim = (TH1D*)hPrim->Project3D("x");
    TH1D *h1All  = (TH1D*)hAll ->Project3D("x");
    h1Prim->Divide(h1All);
    delete h1All;
    for(Int_t ix=1;ix<=h1Prim->GetNbinsX();ix++) {
      Double_t intVal = h1Prim->GetBinContent(ix);
      if(!intVal) continue;
      for(Int_t iy=1;iy<=h2Prim->GetNbinsY();iy++) {
        h2Prim->SetBinContent(ix,iy,h2Prim->GetBinContent(ix,iy)/intVal);
        h2Prim->SetBinError(ix,iy,h2Prim->GetBinError(ix,iy)/intVal);
      };
    };
    delete h1Prim;
  }
  hPrim->GetYaxis()->SetRange(1,hPrim->GetNbinsY());
  hAll->GetYaxis()->SetRange(1,hAll->GetNbinsY());
  return h2Prim;
}
