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
  fIsMC(0), fChi2DoF(4), fTPCNcls(70), fPtmin(0.2), fPtmax(2), fEta(0.8), fBit(96),
  fIsPileUpCuts(0), fIsBuildBG(0), fIsBuildLG(0), 
  fPosBackgroundHist(0), fNegBackgroundHist(0), fChargedBackgroundHist(0),
  fMCPosBackgroundHist(0), fMCNegBackgroundHist(0), fMCChargedBackgroundHist(0), fNeventCentHist(0),
  fEventCuts(0)
{

}
//_____________________________________________________________________________
AliAnalysisTaskLegendreCoef::AliAnalysisTaskLegendreCoef(const char* name) : AliAnalysisTaskSE(name),
  fAOD(0), fOutputList(0),
  fIsMC(0), fChi2DoF(4), fTPCNcls(70), fPtmin(0.2), fPtmax(2), fEta(0.8), fBit(96),
  fIsPileUpCuts(0), fIsBuildBG(0), fIsBuildLG(0), 
  fPosBackgroundHist(0), fNegBackgroundHist(0), fChargedBackgroundHist(0), 
  fMCPosBackgroundHist(0), fMCNegBackgroundHist(0), fMCChargedBackgroundHist(0), fNeventCentHist(0),
  fEventCuts(0)
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
    fOutputList->Add(new TH2D("PosBGHistOut", "postrack;eta;centrality", 16,-fEta,fEta,8,centbins));//output histogram when building background -  positive tracks
    fOutputList->Add(new TH2D("NegBGHistOut", "negtrack;eta;centrality", 16,-fEta,fEta,8,centbins));//output histogram when building background -  negative tracks
    fOutputList->Add(new TH2D("ChargedBGHistOut", "track;eta;centrality", 16,-fEta,fEta,8,centbins));//output histogram when building background -  charged tracks
    if(fIsMC) {
      fOutputList->Add(new TH2D("MCPosBGHistOut", "postrack;eta;centrality", 16,-fEta,fEta,8,centbins));//output MC histogram when building background -  positive tracks
      fOutputList->Add(new TH2D("MCNegBGHistOut", "negtrack;eta;centrality", 16,-fEta,fEta,8,centbins));//output MC histogram when building background -  negative tracks
      fOutputList->Add(new TH2D("MCChargedBGHistOut", "track;eta;centrality", 16,-fEta,fEta,8,centbins));//output MC histogram when building background -  charged tracks  
    }
  }
    
  if(fIsBuildLG) {
    //for charged
    fOutputList->Add(new TProfile("a1","a1;centrality;",8,centbins, ""));//first order legendre cofficient - direct
    fOutputList->Add(new TProfile("Ra1","Ra1;centrality;",8,centbins, ""));//first order legendre cofficient - random Fill(ran->Uniform(-fEta,fEta))
    fOutputList->Add(new TProfile("RDa1","RDa1;centrality;",8,centbins, ""));//first order legendre cofficient - random sample Fill(etaBG[cb]->GetRandom())
    fOutputList->Add(new TProfile("a2","a2;centrality;",8,centbins, ""));//second order legendre cofficient - direct
    fOutputList->Add(new TProfile("Ra2","Ra2;centrality;",8,centbins, ""));//second order legendre cofficient - random Fill(ran->Uniform(-fEta,fEta))
    fOutputList->Add(new TProfile("RDa2","RDa2;centrality;",8,centbins, ""));//second order legendre cofficient - random sample Fill(etaBG[cb]->GetRandom())
    //for positive
    fOutputList->Add(new TProfile("a1pos","a1;centrality;",8,centbins, ""));//first order legendre cofficient - direct
    fOutputList->Add(new TProfile("Ra1pos","Ra1;centrality;",8,centbins, ""));//first order legendre cofficient - random Fill(ran->Uniform(-fEta,fEta))
    fOutputList->Add(new TProfile("RDa1pos","RDa1;centrality;",8,centbins, ""));//first order legendre cofficient - random sample Fill(etaBG[cb]->GetRandom())
    fOutputList->Add(new TProfile("a2pos","a2;centrality;",8,centbins, ""));//second order legendre cofficient - direct
    fOutputList->Add(new TProfile("Ra2pos","Ra2;centrality;",8,centbins, ""));//second order legendre cofficient - random Fill(ran->Uniform(-fEta,fEta))
    fOutputList->Add(new TProfile("RDa2pos","RDa2;centrality;",8,centbins, ""));//second order legendre cofficient - random sample Fill(etaBG[cb]->GetRandom())
    //for negative
    fOutputList->Add(new TProfile("a1neg","a1;centrality;",8,centbins, ""));//first order legendre cofficient - direct
    fOutputList->Add(new TProfile("Ra1neg","Ra1;centrality;",8,centbins, ""));//first order legendre cofficient - random Fill(ran->Uniform(-fEta,fEta))
    fOutputList->Add(new TProfile("RDa1neg","RDa1;centrality;",8,centbins, ""));//first order legendre cofficient - random sample Fill(etaBG[cb]->GetRandom())
    fOutputList->Add(new TProfile("a2neg","a2;centrality;",8,centbins, ""));//second order legendre cofficient - direct
    fOutputList->Add(new TProfile("Ra2neg","Ra2;centrality;",8,centbins, ""));//second order legendre cofficient - random Fill(ran->Uniform(-fEta,fEta))
    fOutputList->Add(new TProfile("RDa2neg","RDa2;centrality;",8,centbins, ""));//second order legendre cofficient - random sample Fill(etaBG[cb]->GetRandom())
    if(fIsMC) {
      //for charged
      fOutputList->Add(new TProfile("a1MC","a1;centrality;",8,centbins, ""));//first order legendre cofficient - direct
      fOutputList->Add(new TProfile("Ra1MC","Ra1;centrality;",8,centbins, ""));//first order legendre cofficient - random Fill(ran->Uniform(-fEta,fEta))
      fOutputList->Add(new TProfile("RDa1MC","RDa1;centrality;",8,centbins, ""));//first order legendre cofficient - random sample Fill(etaBG[cb]->GetRandom())
      fOutputList->Add(new TProfile("a2MC","a2;centrality;",8,centbins, ""));//second order legendre cofficient - direct
      fOutputList->Add(new TProfile("Ra2MC","Ra2;centrality;",8,centbins, ""));//second order legendre cofficient - random Fill(ran->Uniform(-fEta,fEta))
      fOutputList->Add(new TProfile("RDa2MC","RDa2;centrality;",8,centbins, ""));//second order legendre cofficient - random sample Fill(etaBG[cb]->GetRandom())
      //for positive
      fOutputList->Add(new TProfile("a1posMC","a1;centrality;",8,centbins, ""));//first order legendre cofficient - direct
      fOutputList->Add(new TProfile("Ra1posMC","Ra1;centrality;",8,centbins, ""));//first order legendre cofficient - random Fill(ran->Uniform(-fEta,fEta))
      fOutputList->Add(new TProfile("RDa1posMC","RDa1;centrality;",8,centbins, ""));//first order legendre cofficient - random sample Fill(etaBG[cb]->GetRandom())
      fOutputList->Add(new TProfile("a2posMC","a2;centrality;",8,centbins, ""));//second order legendre cofficient - direct
      fOutputList->Add(new TProfile("Ra2posMC","Ra2;centrality;",8,centbins, ""));//second order legendre cofficient - random Fill(ran->Uniform(-fEta,fEta))
      fOutputList->Add(new TProfile("RDa2posMC","RDa2;centrality;",8,centbins, ""));//second order legendre cofficient - random sample Fill(etaBG[cb]->GetRandom())
      //for negative
      fOutputList->Add(new TProfile("a1negMC","a1;centrality;",8,centbins, ""));//first order legendre cofficient - direct
      fOutputList->Add(new TProfile("Ra1negMC","Ra1;centrality;",8,centbins, ""));//first order legendre cofficient - random Fill(ran->Uniform(-fEta,fEta))
      fOutputList->Add(new TProfile("RDa1negMC","RDa1;centrality;",8,centbins, ""));//first order legendre cofficient - random sample Fill(etaBG[cb]->GetRandom())
      fOutputList->Add(new TProfile("a2negMC","a2;centrality;",8,centbins, ""));//second order legendre cofficient - direct
      fOutputList->Add(new TProfile("Ra2negMC","Ra2;centrality;",8,centbins, ""));//second order legendre cofficient - random Fill(ran->Uniform(-fEta,fEta))
      fOutputList->Add(new TProfile("RDa2negMC","RDa2;centrality;",8,centbins, ""));//second order legendre cofficient - random sample Fill(etaBG[cb]->GetRandom())
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

  if(!fIsMC && fIsPileUpCuts){
    fEventCuts.fUseITSTPCCluCorrelationCut = true;
    if (!fEventCuts.AcceptEvent(fAOD)) return;
  }

  //making a cut in pvz -8 to 8cm
  const AliAODVertex *PrimaryVertex = fAOD->GetVertex(0);
  if(!PrimaryVertex) return;
  Float_t PVz = PrimaryVertex->GetZ();
  if(fabs(PVz)>8.0) return;
  
  AliMultSelection *MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
  if(!MultSelection) return;
  
  Float_t Cent = MultSelection->GetMultiplicityPercentile("V0M"); //centrality
  //mCentCL0 = MultSelection->GetMultiplicityPercentile("CL0");
  //mCentCL1 = MultSelection->GetMultiplicityPercentile("CL1");
  if(Cent>70.0 || Cent<0.01) return;
  ((TH1D*) fOutputList->FindObject("NeventsCentHist"))->Fill(Cent);//Nevents vs centrality

//fill hist nevent vs cent
  for(Int_t i(0); i < fAOD->GetNumberOfTracks(); i++) {                 // loop over all these tracks
    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
    if(!track) continue;
    if(fabs(track->Eta()) > fEta) continue;//eta cut
    if(track->Pt() < fPtmin|| track->Pt() > fPtmax) continue; //pt cut
    if(track->GetTPCNcls()<fTPCNcls || track->Chi2perNDF() > fChi2DoF) continue;// cut in TPC Ncls and chi2/dof   
    if(track->TestFilterBit(fBit)) {
      //build background
      //printf("filter bit is %i\n",fBit);

      ((TH2D*) fOutputList->FindObject("ChargedBGHistOut"))->Fill(track->Eta(), Cent);
      if(track->Charge() > 0) ((TH2D*) fOutputList->FindObject("PosBGHistOut"))->Fill(track->Eta(), Cent);
      if(track->Charge() < 0) ((TH2D*) fOutputList->FindObject("NegBGHistOut"))->Fill(track->Eta(), Cent);
    }
  }
  
  if(fIsMC){
    TClonesArray *stack = 0;
    TList *lst = fAOD->GetList();
    stack = (TClonesArray*)lst->FindObject(AliAODMCParticle::StdBranchName());
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
      Bool_t isParticleFromOutOfBunchPileUpEvent = kFALSE;
      isParticleFromOutOfBunchPileUpEvent = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(-1, mcHeader, stack);
      if(isParticleFromOutOfBunchPileUpEvent) return;
    }
    
    for (Int_t i(0); i < nMCTracks; i++) {
      AliAODMCParticle *p1=(AliAODMCParticle*)stack->UncheckedAt(i);
      if (!p1) continue;
      if(p1->Charge()!=-3 && p1->Charge()!=+3) continue;// x3 by convention
      if(!p1->IsPrimary()) continue;
      if(!p1->IsPhysicalPrimary()) continue;
      if(fabs(p1->Eta()) > fEta ) continue;
      if(p1->Pt() < fPtmin|| p1->Pt() > fPtmax) continue;
      if((fabs(p1->GetPdgCode())==211)||(fabs(p1->GetPdgCode())==2212)||(fabs(p1->GetPdgCode())==321)){
        //build background
        ((TH2D*) fOutputList->FindObject("MCChargedBGHistOut"))->Fill(p1->Eta(), Cent);
        if(p1->Charge() > 0) ((TH2D*) fOutputList->FindObject("MCPosBGHistOut"))->Fill(p1->Eta(), Cent);
        if(p1->Charge() < 0) ((TH2D*) fOutputList->FindObject("MCNegBGHistOut"))->Fill(p1->Eta(), Cent);
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

  if(!fIsMC && fIsPileUpCuts){
    fEventCuts.fUseITSTPCCluCorrelationCut = true;
    if (!fEventCuts.AcceptEvent(fAOD)) return;
  }

  //making a cut in pvz -8 to 8cm
  const AliAODVertex *PrimaryVertex = fAOD->GetVertex(0);
  if(!PrimaryVertex) return;
  Float_t PVz = PrimaryVertex->GetZ();
  if(fabs(PVz)>8.0) return;
  
  AliMultSelection *MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
  if(!MultSelection) return;
  
  Float_t Cent = MultSelection->GetMultiplicityPercentile("V0M"); //centrality
  //mCentCL0 = MultSelection->GetMultiplicityPercentile("CL0");
  //mCentCL1 = MultSelection->GetMultiplicityPercentile("CL1");
  if(Cent>70.0 || Cent<0.01) return;
  double centbins[9] = {0.01,5,10,20,30,40,50,60,70};
  TAxis *centaxis = new TAxis(8,centbins);
  TH1D *posSignal, *negSignal, *chargedSignal;
  TH1D *MCposSignal, *MCnegSignal, *MCchargedSignal;

  posSignal = new TH1D("PosSignal", "postrack;eta;centrality", 16,-fEta,fEta);
  negSignal = new TH1D("NegSignal", "postrack;eta;centrality", 16,-fEta,fEta);
  chargedSignal = new TH1D("chargedSignal", "postrack;eta;centrality", 16,-fEta,fEta);
  MCposSignal = new TH1D("MCPosSignal", "postrack;eta;centrality", 16,-fEta,fEta);
  MCnegSignal = new TH1D("MCNegSignal", "postrack;eta;centrality", 16,-fEta,fEta);
  MCchargedSignal = new TH1D("MCchargedSignal", "postrack;eta;centrality", 16,-fEta,fEta);

  for(Int_t i(0); i < fAOD->GetNumberOfTracks(); i++) {                 // loop over all these tracks
    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
    if(!track) continue;
    if(fabs(track->Eta()) > fEta) continue;//eta cut
    if(track->Pt() < fPtmin|| track->Pt() > fPtmax) continue; //pt cut
    if(track->GetTPCNcls()<fTPCNcls || track->Chi2perNDF() > fChi2DoF) continue;// cut in TPC Ncls and chi2/dof   
    if(track->TestFilterBit(fBit)) {
      // calculate signal    
      //printf("filter bit is %i\n",fBit);

      chargedSignal->Fill(track->Eta());
      if(track->Charge() > 0) posSignal->Fill(track->Eta());
      if(track->Charge() < 0) negSignal->Fill(track->Eta());
    }
  }
  
  //calculate coefficients if histogram is full
  int ncentbin = centaxis->FindBin(Cent);
  int neventcent = fNeventCentHist->GetBinContent(ncentbin);
  if(chargedSignal->Integral()>0) {
    TH1D *chargedBG = fChargedBackgroundHist->ProjectionX("chargedbackground", ncentbin,ncentbin);
    chargedBG->Scale(1.0/neventcent);
    BuildCoefficients(chargedSignal, chargedBG, Cent, "");
  }
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

  if(fIsMC){
    TClonesArray *stack = 0;
    TList *lst = fAOD->GetList();
    stack = (TClonesArray*)lst->FindObject(AliAODMCParticle::StdBranchName());
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
      Bool_t isParticleFromOutOfBunchPileUpEvent = kFALSE;
      isParticleFromOutOfBunchPileUpEvent = AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(-1, mcHeader, stack);
      if(isParticleFromOutOfBunchPileUpEvent) return;
    }
  
    for (Int_t i(0); i < nMCTracks; i++) {
      AliAODMCParticle *p1=(AliAODMCParticle*)stack->UncheckedAt(i);
      if (!p1) continue;
      if(p1->Charge()!=-3 && p1->Charge()!=+3) continue;// x3 by convention
      if(!p1->IsPrimary()) continue;
      if(!p1->IsPhysicalPrimary()) continue;
      if(fabs(p1->Eta()) > fEta ) continue;
      if(p1->Pt() < fPtmin|| p1->Pt() > fPtmax) continue;
      if((fabs(p1->GetPdgCode())==211)||(fabs(p1->GetPdgCode())==2212)||(fabs(p1->GetPdgCode())==321)){
        // calculate signal    
        MCchargedSignal->Fill(p1->Eta());
        if(p1->Charge() > 0) MCposSignal->Fill(p1->Eta());
        if(p1->Charge() < 0) MCnegSignal->Fill(p1->Eta());        
      }
    }
    //calculate coefficients if histogram is full
    if(MCchargedSignal->Integral()>0) {
      TH1D *MCchargedBG = fMCChargedBackgroundHist->ProjectionX("MCchargedbackground", ncentbin,ncentbin);
      MCchargedBG->Scale(1.0/neventcent);
      BuildCoefficients(MCchargedSignal, MCchargedBG, Cent, "MC");
    }
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
void AliAnalysisTaskLegendreCoef::BuildCoefficients(TH1D *signal, TH1D *background, Float_t centrality, char *type)
{
  TH1D *RanDistHist[5],*RanHist[5];//random distribution histograms
  TRandom2 *ran = new TRandom2();//random number for uniform distribution
  double ntracks = signal->Integral();
  int n; 
  char histname[50];
  for (int s=0; s<5; s++){//5 random histograms
    n = sprintf (histname, "RanHist%i", s+1);
    RanHist[s] = new TH1D(histname,histname, 16, -fEta, fEta);
    RanHist[s]->Sumw2();
    n = sprintf (histname, "RanDistHist%i", s+1);
    RanDistHist[s] = new TH1D(histname,histname, 16, -fEta, fEta);
    RanDistHist[s]->Sumw2();
  }

  //generating random uniform distributions
  for (int s=0; s<5; s++){
    RanHist[s]->Reset("ICE");
    for (int rn=0; rn<ntracks; rn++) {RanHist[s]->Fill(ran->Uniform(-fEta,fEta));}
        //printf("first ranhist normal is %f\n",RanHist[s]->Integral());
    RanHist[s]->Scale(16.0/(double)ntracks);
    //printf("ranhist normal is %f\n",RanHist[s]->Integral());

  }

  //generating random sample distributions
  for (int s=0; s<5; s++){
    RanDistHist[s]->Reset("ICE");
    for (int rn=0; rn<ntracks; rn++) {RanDistHist[s]->Fill(background->GetRandom());}
    RanDistHist[s]->Divide(background);
    //printf("RanDistHist[s] normal is %f\n",RanDistHist[s]->Integral());
  }

  //normalizing signal hist
  signal->Divide(background);
  signal->Scale(16.0/signal->Integral());

  //printf("signal normal is %f\n",signal->Integral());

  //calculating the an coefficients  
  n = sprintf (histname, "a1%s", type);
  ((TProfile*) fOutputList->FindObject(histname))->Fill(centrality, pow(GetSingleAnCoef(1,signal),2.0));
  n = sprintf (histname, "a2%s", type);
  ((TProfile*) fOutputList->FindObject(histname))->Fill(centrality, pow(GetSingleAnCoef(2,signal),2.0));
  for (int s=0; s<5; s++){
    n = sprintf (histname, "Ra1%s", type);
    ((TProfile*) fOutputList->FindObject(histname))->Fill(centrality, pow(GetSingleAnCoef(1,RanHist[s]),2.0));
    n = sprintf (histname, "RDa1%s", type);
    ((TProfile*) fOutputList->FindObject(histname))->Fill(centrality, pow(GetSingleAnCoef(1,RanDistHist[s]),2.0));
    n = sprintf (histname, "Ra2%s", type);
    ((TProfile*) fOutputList->FindObject(histname))->Fill(centrality, pow(GetSingleAnCoef(2,RanHist[s]),2.0));
    n = sprintf (histname, "RDa2%s", type);
    ((TProfile*) fOutputList->FindObject(histname))->Fill(centrality, pow(GetSingleAnCoef(2,RanDistHist[s]),2.0));
  }  

  for (int s=0; s<5; s++){
    delete (RanHist[s]);
    delete (RanDistHist[s]);    
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
  for(Int_t i=1;i<=hist->GetNbinsX();i++){
    coef += sqrt((1.0/fEta)*(order+0.5))*(hist->GetBinContent(i)-1)*LegPol(order,hist->GetBinCenter(i));
  }
  coef = coef*hist->GetBinWidth(1);    
  return coef;
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskLegendreCoef::LegPol(int order, Double_t x){
    if( order == 1 ) return x/fEta;
    if( order == 2 ) return 0.5*(3.0*pow(x/fEta,2)-1.0);
    if( order == 3 ) return 0.5*(5.0*pow(x/fEta,3)-3.0*x/fEta);
    if( order == 4 ) return (35.0*pow(x/fEta,4)-30.0*pow(x/fEta,2)+3.0)/8.0;
    else return 0.0;
}
//_____________________________________________________________________________
void AliAnalysisTaskLegendreCoef::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
