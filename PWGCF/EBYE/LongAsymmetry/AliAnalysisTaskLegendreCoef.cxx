////////////////////////////////////////////////////////
// AliAnalysisTaskLegendreCoef:
// Description: Analysis task computes the background 
// and extracts legendre coefficients from eta dist
// Author: Raquel Quishpe (raquel.quishpe@cern.ch)
////////////////////////////////////////////////////////
#include "TChain.h"
#include "TH3D.h"
#include "TTree.h"
#include "TProfile.h"
#include "TString.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAODTrack.h"
#include "AliMultSelection.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskLegendreCoef.h"
#include "AliAODMCParticle.h"
#include "AliEventCuts.h"
#include "AliAODMCHeader.h"
#include "TRandom2.h"

ClassImp(AliAnalysisTaskLegendreCoef)

AliAnalysisTaskLegendreCoef::AliAnalysisTaskLegendreCoef() : AliAnalysisTaskSE(),
  fAOD(0), fOutputList(0),
  fIsMC(0), fChi2DoF(3), fTPCNcls(70), fPtmin(0.2), fPtmax(2), fEtaMin(-0.8), fEtaMax(0.8), fBit(96),
  fIsPileUpCuts(0), fIsBuildBG(0), fIsBuildLG(0), 
  fPosBackgroundHist(0), fNegBackgroundHist(0), fChargedBackgroundHist(0),
  fMCPosBackgroundHist(0), fMCNegBackgroundHist(0), fMCChargedBackgroundHist(0), fNeventCentHist(0),
  fGenName("Hijing"), fPileUpLevel(2), fTPCNCrossedRows(70), fPVzMax(8.0), fPVzMin(0.0), fPVzSign(0), fNetabins(16), fIsRunFBOnly(0), fEventCuts(0)
{

}
//_____________________________________________________________________________
AliAnalysisTaskLegendreCoef::AliAnalysisTaskLegendreCoef(const char* name) : AliAnalysisTaskSE(name),
  fAOD(0), fOutputList(0),
  fIsMC(0), fChi2DoF(3), fTPCNcls(70), fPtmin(0.2), fPtmax(2), fEtaMin(-0.8), fEtaMax(0.8), fBit(96),
  fIsPileUpCuts(0), fIsBuildBG(0), fIsBuildLG(0), 
  fPosBackgroundHist(0), fNegBackgroundHist(0), fChargedBackgroundHist(0), 
  fMCPosBackgroundHist(0), fMCNegBackgroundHist(0), fMCChargedBackgroundHist(0), fNeventCentHist(0),
  fGenName("Hijing"), fPileUpLevel(2), fTPCNCrossedRows(70), fPVzMax(8.0), fPVzMin(0.0), fPVzSign(0), fNetabins(16), fIsRunFBOnly(0), fEventCuts(0)
{
  // Default constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container
  DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskLegendreCoef::~AliAnalysisTaskLegendreCoef()
{
  // destructor
  if(fOutputList) {
    delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
  }
}
//_____________________________________________________________________________

void AliAnalysisTaskLegendreCoef::UserCreateOutputObjects()
{
  // Initialize output list of containers
  if (fOutputList != NULL) {
    delete fOutputList;
    fOutputList = NULL;
  }
  if (!fOutputList) {
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
  }
  double centbins[9] = {0.01,5,10,20,30,40,50,60,70};


  if(fIsBuildBG){
    fOutputList->Add(new TH1D("NeventsCentHist", "nevents;centrality", 8,centbins));//output histogram when building background -  all tracks
    fOutputList->Add(new TH2D("PosBGHistOut", "postrack;eta;centrality", fNetabins,fEtaMin,fEtaMax,8,centbins));//output histogram when building background -  positive tracks
    fOutputList->Add(new TH2D("NegBGHistOut", "negtrack;eta;centrality", fNetabins,fEtaMin,fEtaMax,8,centbins));//output histogram when building background -  negative tracks
    fOutputList->Add(new TH2D("ChargedBGHistOut", "track;eta;centrality", fNetabins,fEtaMin,fEtaMax,8,centbins));//output histogram when building background -  charged tracks
    fOutputList->Add(new TH3D("PhiEtaCentHist", "track;phi;eta;centrality", 24, 0 ,TMath::Pi()*2.0, fNetabins,fEtaMin,fEtaMax, 100, 0.01, 70.0));//QA hist phi vs eta
    fOutputList->Add(new TH2D("PtCentHist", "track;pt;centrality", 70,fPtmin,fPtmax,8,centbins));//QA hist pt 
    if(fIsMC) {
      fOutputList->Add(new TH1D("GenNeventsCentHist", "nevents;centrality", 8,centbins));//output histogram when building background -  all tracks
      fOutputList->Add(new TH1D("RecNeventsCentHist", "nevents;centrality", 8,centbins));//output histogram when building background -  all tracks
      fOutputList->Add(new TH2D("MCPosBGHistOut", "postrack;eta;centrality", fNetabins,fEtaMin,fEtaMax,8,centbins));//output MC histogram when building background -  positive tracks
      fOutputList->Add(new TH2D("MCNegBGHistOut", "negtrack;eta;centrality", fNetabins,fEtaMin,fEtaMax,8,centbins));//output MC histogram when building background -  negative tracks
      fOutputList->Add(new TH2D("MCChargedBGHistOut", "track;eta;centrality", fNetabins,fEtaMin,fEtaMax,8,centbins));//output MC histogram when building background -  charged tracks  
      fOutputList->Add(new TH3D("GenPhiEtaCentHist", "track;phi;eta;centrality", 24, 0 ,TMath::Pi()*2.0, fNetabins,fEtaMin,fEtaMax, 100, 0.01, 70.0));//QA hist phi vs eta
      fOutputList->Add(new TH3D("RecPhiEtaCentHist", "track;phi;eta;centrality", 24, 0 ,TMath::Pi()*2.0, fNetabins,fEtaMin,fEtaMax, 100, 0.01, 70.0));//QA hist phi vs eta
      fOutputList->Add(new TH2D("GenPtCentHist", "track;pt;centrality", 70,fPtmin,fPtmax,8,centbins));//QA hist pt 
      fOutputList->Add(new TH2D("RecPtCentHist", "track;pt;centrality", 70,fPtmin,fPtmax,8,centbins));//QA hist pt
    }
  }
    
  if(fIsBuildLG) {
    //for charged
    fOutputList->Add(new TProfile("a1","a1;centrality;",8,centbins, ""));//first order legendre cofficient - direct
    fOutputList->Add(new TProfile("Ra1","Ra1;centrality;",8,centbins, ""));//first order legendre cofficient - random Fill(ran->Uniform(-fEta,fEta))
    fOutputList->Add(new TProfile("a2","a2;centrality;",8,centbins, ""));//second order legendre cofficient - direct
    fOutputList->Add(new TProfile("Ra2","Ra2;centrality;",8,centbins, ""));//second order legendre cofficient - random Fill(ran->Uniform(-fEta,fEta))
    //for positive
    fOutputList->Add(new TProfile("a1pos","a1;centrality;",8,centbins, ""));//first order legendre cofficient - direct
    fOutputList->Add(new TProfile("Ra1pos","Ra1;centrality;",8,centbins, ""));//first order legendre cofficient - random Fill(ran->Uniform(-fEta,fEta))
    fOutputList->Add(new TProfile("a2pos","a2;centrality;",8,centbins, ""));//second order legendre cofficient - direct
    fOutputList->Add(new TProfile("Ra2pos","Ra2;centrality;",8,centbins, ""));//second order legendre cofficient - random Fill(ran->Uniform(-fEta,fEta))
    //for negative
    fOutputList->Add(new TProfile("a1neg","a1;centrality;",8,centbins, ""));//first order legendre cofficient - direct
    fOutputList->Add(new TProfile("Ra1neg","Ra1;centrality;",8,centbins, ""));//first order legendre cofficient - random Fill(ran->Uniform(-fEta,fEta))
    fOutputList->Add(new TProfile("a2neg","a2;centrality;",8,centbins, ""));//second order legendre cofficient - direct
    fOutputList->Add(new TProfile("Ra2neg","Ra2;centrality;",8,centbins, ""));//second order legendre cofficient - random Fill(ran->Uniform(-fEta,fEta))
    if(fIsMC) {
      //for charged
      fOutputList->Add(new TProfile("a1MC","a1;centrality;",8,centbins, ""));//first order legendre cofficient - direct
      fOutputList->Add(new TProfile("Ra1MC","Ra1;centrality;",8,centbins, ""));//first order legendre cofficient - random Fill(ran->Uniform(-fEta,fEta))
      fOutputList->Add(new TProfile("a2MC","a2;centrality;",8,centbins, ""));//second order legendre cofficient - direct
      fOutputList->Add(new TProfile("Ra2MC","Ra2;centrality;",8,centbins, ""));//second order legendre cofficient - random Fill(ran->Uniform(-fEta,fEta))
      //for positive
      fOutputList->Add(new TProfile("a1posMC","a1;centrality;",8,centbins, ""));//first order legendre cofficient - direct
      fOutputList->Add(new TProfile("Ra1posMC","Ra1;centrality;",8,centbins, ""));//first order legendre cofficient - random Fill(ran->Uniform(-fEta,fEta))
      fOutputList->Add(new TProfile("a2posMC","a2;centrality;",8,centbins, ""));//second order legendre cofficient - direct
      fOutputList->Add(new TProfile("Ra2posMC","Ra2;centrality;",8,centbins, ""));//second order legendre cofficient - random Fill(ran->Uniform(-fEta,fEta))
      //for negative
      fOutputList->Add(new TProfile("a1negMC","a1;centrality;",8,centbins, ""));//first order legendre cofficient - direct
      fOutputList->Add(new TProfile("Ra1negMC","Ra1;centrality;",8,centbins, ""));//first order legendre cofficient - random Fill(ran->Uniform(-fEta,fEta))
      fOutputList->Add(new TProfile("a2negMC","a2;centrality;",8,centbins, ""));//second order legendre cofficient - direct
      fOutputList->Add(new TProfile("Ra2negMC","Ra2;centrality;",8,centbins, ""));//second order legendre cofficient - random Fill(ran->Uniform(-fEta,fEta))
    }
  }

  PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskLegendreCoef::BuildBackground()
{
  //printf("Building background!\n");

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD) return;

  if(fIsPileUpCuts){
    fEventCuts.fUseITSTPCCluCorrelationCut = fPileUpLevel;
    if (!fEventCuts.AcceptEvent(fAOD)) return;
  }

  TClonesArray *stack =0;
  if(fIsMC) {
    TList *lst = fAOD->GetList();
    stack = (TClonesArray*)lst->FindObject(AliAODMCParticle::StdBranchName());
    if(!stack) return;
  }

  //making a cut in pvz -8 to 8cm
  const AliAODVertex *PrimaryVertex = fAOD->GetVertex(0);
  if(!PrimaryVertex) return;
  Float_t PVz = PrimaryVertex->GetZ();
  if(fabs(PVz)>fPVzMax) return;
  if(fabs(PVz)<fPVzMin) return;
  if(fPVzSign==1 && PVz<0) return;//take only positive pvz
  if(fPVzSign==-1 && PVz>0) return;//take only negative pvz

  AliMultSelection *MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
  if(!MultSelection) return;
  
  Float_t Cent = MultSelection->GetMultiplicityPercentile("V0M"); //centrality
  if(Cent>70.0 || Cent<0.01) return;
  ((TH1D*) fOutputList->FindObject("NeventsCentHist"))->Fill(Cent);//Nevents vs centrality

  if(fIsMC){
    int nMCTracks;
    if (!stack) nMCTracks = 0;
    else nMCTracks = stack->GetEntries();
    
    if(fIsPileUpCuts){
      AliAODMCHeader *mcHeader = 0;
      mcHeader = (AliAODMCHeader*)fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName());
      if(!mcHeader) {
        printf("AliAnalysisTaskSEHFTreeCreator::UserExec: MC header branch not found!\n");
        return;
      }
      Bool_t isPileupInGeneratedEvent = kFALSE;
      isPileupInGeneratedEvent = AliAnalysisUtils::IsPileupInGeneratedEvent(mcHeader, fGenName);
      if(isPileupInGeneratedEvent) return;
    }
    ((TH1D*) fOutputList->FindObject("GenNeventsCentHist"))->Fill(Cent);//Nevents vs centrality
    for (Int_t i(0); i < nMCTracks; i++) {
      AliAODMCParticle *p1=(AliAODMCParticle*)stack->UncheckedAt(i);
      if (!p1) continue;
      if(abs(p1->Charge())<=1) continue;
      if(!p1->IsPhysicalPrimary()) continue;
      if((fabs(p1->GetPdgCode())==211)||(fabs(p1->GetPdgCode())==2212)||(fabs(p1->GetPdgCode())==321)){
        //build background
        Float_t pTMCgen   = p1->Pt();
        Float_t phiMCgen  = p1->Phi();
        Float_t etaMCgen  = p1->Eta();
        if(etaMCgen> fEtaMax || etaMCgen<fEtaMin ) continue;
        if(pTMCgen < fPtmin|| pTMCgen > fPtmax) continue;
        ((TH2D*) fOutputList->FindObject("MCChargedBGHistOut"))->Fill(etaMCgen, Cent);
        if(p1->Charge() > 0) ((TH2D*) fOutputList->FindObject("MCPosBGHistOut"))->Fill(etaMCgen, Cent);
        if(p1->Charge() < 0) ((TH2D*) fOutputList->FindObject("MCNegBGHistOut"))->Fill(etaMCgen, Cent);
        ((TH3D*) fOutputList->FindObject("GenPhiEtaCentHist"))->Fill(phiMCgen,etaMCgen, Cent);
        ((TH2D*) fOutputList->FindObject("GenPtCentHist"))->Fill(pTMCgen, Cent);
      }
    }
  }

//fill hist nevent vs cent
  for(Int_t i(0); i < fAOD->GetNumberOfTracks(); i++) {                 // loop over all these tracks
    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
    if(!track) continue;
    if(!fIsMC){
      if(track->Charge()==0)continue;//only get charged tracks
      if(track->Eta() > fEtaMax || track->Eta() < fEtaMin) continue;//eta cut
      if(track->Pt() < fPtmin|| track->Pt() > fPtmax) continue; //pt cut
      if(!fIsRunFBOnly){
        if(track->GetTPCNcls()<fTPCNcls || track->GetTPCNCrossedRows()<fTPCNCrossedRows || track->Chi2perNDF() > fChi2DoF) continue;// cut in TPC Ncls , crossed rows and chi2/dof  
      }
      if(!track->TestFilterBit(fBit)) continue;
      //build background for raw data
      //printf("filter bit is %i\n",fBit);
      ((TH2D*) fOutputList->FindObject("ChargedBGHistOut"))->Fill(track->Eta(), Cent);
      if(track->Charge() > 0) ((TH2D*) fOutputList->FindObject("PosBGHistOut"))->Fill(track->Eta(), Cent);
      if(track->Charge() < 0) ((TH2D*) fOutputList->FindObject("NegBGHistOut"))->Fill(track->Eta(), Cent);
      ((TH3D*) fOutputList->FindObject("PhiEtaCentHist"))->Fill(track->Phi(),track->Eta(), Cent);
      ((TH2D*) fOutputList->FindObject("PtCentHist"))->Fill(track->Pt(), Cent);
    }else{
      //build background for MC reconstructed
      int label = TMath::Abs(track->GetLabel());
      AliAODMCParticle* mcTrack = dynamic_cast<AliAODMCParticle*>(stack->At(label));
      if (!mcTrack) continue;
      if(abs(mcTrack->Charge())<=1)continue;//only get charged tracks
      if(!mcTrack->IsPhysicalPrimary()) continue;  
      if((fabs(mcTrack->GetPdgCode())==211)||(fabs(mcTrack->GetPdgCode())==2212)||(fabs(mcTrack->GetPdgCode())==321)){
        Float_t pTrec   = mcTrack->Pt();
        Float_t phirec  = mcTrack->Phi();
        Float_t etarec  = mcTrack->Eta();
        Float_t ncl = track->GetTPCNcls();
        Float_t crossedrows = track->GetTPCNCrossedRows();
        Float_t chi2 = track->Chi2perNDF();
        Int_t isfb = 0;
        if(track->TestFilterBit(fBit)) isfb=1;
        if(isfb==0) continue;//choose filterbit
        if(abs(mcTrack->Charge())<=1)continue;//only get charged tracks
        if(etarec > fEtaMax || etarec < fEtaMin) continue;//eta cut
        if(pTrec < fPtmin|| pTrec > fPtmax) continue; //pt cut
        if(!fIsRunFBOnly){         
          if(ncl<fTPCNcls || crossedrows<fTPCNCrossedRows || chi2 > fChi2DoF) continue;// cut in TPC Ncls , crossed rows and chi2/dof  
        }
        ((TH2D*) fOutputList->FindObject("ChargedBGHistOut"))->Fill(etarec, Cent);
        if(mcTrack->Charge() > 0) ((TH2D*) fOutputList->FindObject("PosBGHistOut"))->Fill(etarec, Cent);
        if(mcTrack->Charge() < 0) ((TH2D*) fOutputList->FindObject("NegBGHistOut"))->Fill(etarec, Cent);
        ((TH3D*) fOutputList->FindObject("RecPhiEtaCentHist"))->Fill(phirec, etarec, Cent);
        ((TH2D*) fOutputList->FindObject("RecPtCentHist"))->Fill(pTrec, Cent);     
      }
    }
  }
  return;
}

//_____________________________________________________________________________
void AliAnalysisTaskLegendreCoef::BuildSignal()
{
  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD) return;

  if(fIsPileUpCuts){
    fEventCuts.fUseITSTPCCluCorrelationCut = fPileUpLevel;
    if (!fEventCuts.AcceptEvent(fAOD)) return;
  }

  TClonesArray *stack =0;
  if(fIsMC) {
    TList *lst = fAOD->GetList();
    stack = (TClonesArray*)lst->FindObject(AliAODMCParticle::StdBranchName());
    if(!stack) return;
  }

  //making a cut in pvz -8 to 8cm
  const AliAODVertex *PrimaryVertex = fAOD->GetVertex(0);
  if(!PrimaryVertex) return;
  Float_t PVz = PrimaryVertex->GetZ();
  if(fabs(PVz)>fPVzMax) return;
  if(fabs(PVz)<fPVzMin) return;
  if(fPVzSign==1 && PVz<0) return;//take only positive pvz
  if(fPVzSign==-1 && PVz>0) return;//take only negative pvz
  
  AliMultSelection *MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
  if(!MultSelection) return;
  
  Float_t Cent = MultSelection->GetMultiplicityPercentile("V0M"); //centrality
  if(Cent>70.0 || Cent<0.01) return;
  double centbins[9] = {0.01,5,10,20,30,40,50,60,70};
  TAxis *centaxis = new TAxis(8,centbins);
  TH1D *posSignal, *negSignal, *chargedSignal;
  TH1D *MCposSignal, *MCnegSignal, *MCchargedSignal;

  posSignal = new TH1D("PosSignal", "postrack;eta;centrality", fNetabins,fEtaMin,fEtaMax);
  negSignal = new TH1D("NegSignal", "postrack;eta;centrality", fNetabins,fEtaMin,fEtaMax);
  chargedSignal = new TH1D("chargedSignal", "postrack;eta;centrality", fNetabins,fEtaMin,fEtaMax);
  MCposSignal = new TH1D("MCPosSignal", "postrack;eta;centrality", fNetabins,fEtaMin,fEtaMax);
  MCnegSignal = new TH1D("MCNegSignal", "postrack;eta;centrality", fNetabins,fEtaMin,fEtaMax);
  MCchargedSignal = new TH1D("MCchargedSignal", "postrack;eta;centrality", fNetabins,fEtaMin,fEtaMax);


  if(fIsMC){
    int nMCTracks;
    if (!stack) nMCTracks = 0;
    else nMCTracks = stack->GetEntries();
  
    if(fIsPileUpCuts){
      AliAODMCHeader *mcHeader = 0;
      mcHeader = (AliAODMCHeader*)fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName());
      if(!mcHeader) {
        printf("AliAnalysisTaskSEHFTreeCreator::UserExec: MC header branch not found!\n");
        return;
      }
      Bool_t isPileupInGeneratedEvent = kFALSE;
      isPileupInGeneratedEvent = AliAnalysisUtils::IsPileupInGeneratedEvent(mcHeader,fGenName);
      if(isPileupInGeneratedEvent) return;
    }
  
    for (Int_t i(0); i < nMCTracks; i++) {
      AliAODMCParticle *p1=(AliAODMCParticle*)stack->UncheckedAt(i);
      if (!p1) continue;
      if(abs(p1->Charge())<=1)continue;//only get charged tracks
      if(!p1->IsPhysicalPrimary()) continue;
      Float_t etaMCgen = p1->Eta();
      if(etaMCgen > fEtaMax || etaMCgen<fEtaMin ) continue;
      if(p1->Pt() < fPtmin|| p1->Pt() > fPtmax) continue;
      if((fabs(p1->GetPdgCode())==211)||(fabs(p1->GetPdgCode())==2212)||(fabs(p1->GetPdgCode())==321)){
        // calculate signal    
        MCchargedSignal->Fill(etaMCgen);
        if(p1->Charge() > 0) MCposSignal->Fill(etaMCgen);
        if(p1->Charge() < 0) MCnegSignal->Fill(etaMCgen);        
      }
    }
  }

  for(Int_t i(0); i < fAOD->GetNumberOfTracks(); i++) {                 // loop over all these tracks
    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
    if(!track) continue;
    if(!fIsMC){
    //build signal for raw data
      if(track->Charge()==0)continue;//only get charged tracks
      Float_t etaRaw = track->Eta();
      if(etaRaw > fEtaMax || etaRaw < fEtaMin) continue;//eta cut
      if(track->Pt() < fPtmin|| track->Pt() > fPtmax) continue; //pt cut
      if(!fIsRunFBOnly){
        if(track->GetTPCNcls()<fTPCNcls || track->GetTPCNCrossedRows()<fTPCNCrossedRows || track->Chi2perNDF() > fChi2DoF) continue;// cut in TPC Ncls, crossed rows and chi2/dof   
      }
      if(track->TestFilterBit(fBit)) {
        chargedSignal->Fill(etaRaw);
        if(track->Charge() > 0) posSignal->Fill(etaRaw);
        if(track->Charge() < 0) negSignal->Fill(etaRaw);
      }
    }else{
      //build signal for MC reconstructed
      int label = TMath::Abs(track->GetLabel());
      AliAODMCParticle* mcTrack = dynamic_cast<AliAODMCParticle*>(stack->At(label));
      if (!mcTrack) continue;
      if(abs(mcTrack->Charge())<=1)continue;//only get charged tracks
      if(!mcTrack->IsPhysicalPrimary()) continue;    
      Float_t etaRec = mcTrack->Eta();
      if(etaRec> fEtaMax || etaRec < fEtaMin) continue;//eta cut
      if(mcTrack->Pt() < fPtmin|| mcTrack->Pt() > fPtmax) continue; //pt cut 
      if(!fIsRunFBOnly){ 
        if(track->GetTPCNcls()<fTPCNcls || track->GetTPCNCrossedRows()<fTPCNCrossedRows || track->Chi2perNDF() > fChi2DoF) continue;// cut in TPC Ncls, crossed rows and chi2/dof       
      }
      if(!track->TestFilterBit(fBit)) continue;//only filterbit
      if((fabs(mcTrack->GetPdgCode())==211)||(fabs(mcTrack->GetPdgCode())==2212)||(fabs(mcTrack->GetPdgCode())==321)){
        chargedSignal->Fill(etaRec);
        if(mcTrack->Charge() > 0) posSignal->Fill(etaRec);
        if(mcTrack->Charge() < 0) negSignal->Fill(etaRec);
      }
    }
  }

  //calculate coefficients if histogram is full
  int ncentbin = centaxis->FindBin(Cent);
  int neventcent = fNeventCentHist->GetBinContent(ncentbin);

  if(chargedSignal->Integral()>0) {
    TH1D *chargedBG = fChargedBackgroundHist->ProjectionX("chargedbackground", ncentbin,ncentbin);
    // printf("BACKGROUND normal is %f\n",chargedBG->Integral());
    chargedBG->Scale(1.0/neventcent);
    BuildCoefficients(chargedSignal, chargedBG, Cent, "");
    if(posSignal->Integral()>0) {
      TH1D *posBG = fPosBackgroundHist->ProjectionX("posbackground", ncentbin,ncentbin);
      posBG->Scale(1.0/neventcent);
      BuildCoefficients(posSignal, posBG, Cent, "pos");   
    }
    if(negSignal->Integral()>0) {
      TH1D *negBG = fNegBackgroundHist->ProjectionX("negbackground", ncentbin,ncentbin);
      negBG->Scale(1.0/neventcent);
      BuildCoefficients(negSignal, negBG, Cent, "neg");
    }
  }


  //calculate coefficients if histogram is full
  if(fIsMC && MCchargedSignal->Integral()>0) {
    TH1D *MCchargedBG = fMCChargedBackgroundHist->ProjectionX("MCchargedbackground", ncentbin,ncentbin);
    MCchargedBG->Scale(1.0/neventcent);
    BuildCoefficients(MCchargedSignal, MCchargedBG, Cent, "MC");
    if(MCposSignal->Integral()>0) {
      TH1D *MCposBG = fMCPosBackgroundHist->ProjectionX("MCposbackground", ncentbin,ncentbin);
      MCposBG->Scale(1.0/neventcent);
      BuildCoefficients(MCposSignal, MCposBG, Cent, "posMC");   
    }
    if(MCnegSignal->Integral()>0) {
      TH1D *MCnegBG = fMCNegBackgroundHist->ProjectionX("negbackground", ncentbin,ncentbin);
      MCnegBG->Scale(1.0/neventcent);
      BuildCoefficients(MCnegSignal, MCnegBG, Cent, "negMC");
    }
  }
  
  delete (posSignal);
  delete (negSignal); 
  delete (chargedSignal);
  delete (MCposSignal);
  delete (MCnegSignal); 
  delete (MCchargedSignal);
  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskLegendreCoef::BuildCoefficients(TH1D *signal, TH1D *background, Float_t centrality, TString type)
{
  TH1D *RanHist[5];//random distribution histograms
  int ntracks = signal->Integral();
  int n; 
  char histname[50];
  //normalizing signal hist 
  signal->Divide(background);
    // printf("signal normal is %f\n",signal->Integral());

  for (int s=0; s<5; s++){//5 random histograms
    n = sprintf (histname, "RanHist%i", s+1);
    RanHist[s] = new TH1D(histname,histname, fNetabins, fEtaMin, fEtaMax);
    RanHist[s]->Sumw2();
  }

  //generating random uniform distributions
  for (int s=0; s<5; s++){
    RanHist[s]->Reset("ICE");
    for (int rn=0; rn<ntracks; rn++) {
      RanHist[s]->Fill(background->GetRandom());
    }
    RanHist[s]->Divide(background); 
  }

  //calculating the an coefficients  
  ((TProfile*) fOutputList->FindObject("a1"+type))->Fill(centrality, pow(GetSingleAnCoef(1,signal),2.0));
  ((TProfile*) fOutputList->FindObject("a2"+type))->Fill(centrality, pow(GetSingleAnCoef(2,signal),2.0));

  for (int s=0; s<5; s++){
    ((TProfile*) fOutputList->FindObject("Ra1"+type))->Fill(centrality, pow(GetSingleAnCoef(1,RanHist[s]),2.0));
    ((TProfile*) fOutputList->FindObject("Ra2"+type))->Fill(centrality, pow(GetSingleAnCoef(2,RanHist[s]),2.0));
  }  

  for (int s=0; s<5; s++){
    delete (RanHist[s]);
  }
  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskLegendreCoef::UserExec(Option_t *)
{
  if(fIsBuildBG) BuildBackground();
  if(fIsBuildLG) BuildSignal();
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskLegendreCoef::GetSingleAnCoef(int order, TH1D *hist)
{   
  double coef = 0.0;
  double etawidth = (fEtaMax-fEtaMin)/2.0;
  for(Int_t i=1;i<=hist->GetNbinsX();i++){
    coef += sqrt((1.0/etawidth)*(order+0.5))*(hist->GetBinContent(i)-1)*LegPol(order,hist->GetBinCenter(i));
  }
  coef = coef*hist->GetBinWidth(1);    
  return coef;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskLegendreCoef::LegPol(int order, Double_t x){
    double etawidth = (fEtaMax-fEtaMin)/2.0;
    if( order == 1 ) return x/etawidth;
    if( order == 2 ) return 0.5*(3.0*pow(x/etawidth,2)-1.0);
    if( order == 3 ) return 0.5*(5.0*pow(x/etawidth,3)-3.0*x/etawidth);
    if( order == 4 ) return (35.0*pow(x/etawidth,4)-30.0*pow(x/etawidth,2)+3.0)/8.0;
    else return 0.0;
}
//_____________________________________________________________________________
void AliAnalysisTaskLegendreCoef::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
