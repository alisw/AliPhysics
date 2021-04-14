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

ClassImp(AliAnalysisTaskLegendreCoef)

AliAnalysisTaskLegendreCoef::AliAnalysisTaskLegendreCoef() : AliAnalysisTaskSE(),
  fAOD(0), fOutputList(0),
  fIsMC(0), fChi2DoF(4), fTPCNcls(70), fPtmin(0.2), fPtmax(2), fEta(0.8), 
  fIsPileUpCuts(0), fIsBuildBG(0), fIsBuildLG(0), 
  fPosBackgroundHist(0), fNegBackgroundHist(0), fChargedBackgroundHist(0),
  fMCPosBackgroundHist(0), fMCNegBackgroundHist(0), fMCChargedBackgroundHist(0),
  fEventCuts(0)
{

}
//_____________________________________________________________________________
AliAnalysisTaskLegendreCoef::AliAnalysisTaskLegendreCoef(const char* name) : AliAnalysisTaskSE(name),
  fAOD(0), fOutputList(0),
  fIsMC(0), fChi2DoF(4), fTPCNcls(70), fPtmin(0.2), fPtmax(2), fEta(0.8), 
  fIsPileUpCuts(0), fIsBuildBG(0), fIsBuildLG(0), 
  fPosBackgroundHist(0), fNegBackgroundHist(0), fChargedBackgroundHist(0),
  fMCPosBackgroundHist(0), fMCNegBackgroundHist(0), fMCChargedBackgroundHist(0),
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
    fOutputList->Add(new TProfile("a1","a1;centrality;",8,centbins, ""));//first order legendre cofficient - direct
    fOutputList->Add(new TProfile("Ra1","Ra1;centrality;",8,centbins, ""));//first order legendre cofficient - random Fill(ran->Uniform(-fEta,fEta))
    fOutputList->Add(new TProfile("RDa1","RDa1;centrality;",8,centbins, ""));//first order legendre cofficient - random sample Fill(etaBG[cb]->GetRandom())
    fOutputList->Add(new TProfile("a2","a2;centrality;",8,centbins, ""));//second order legendre cofficient - direct
    fOutputList->Add(new TProfile("Ra2","Ra2;centrality;",8,centbins, ""));//second order legendre cofficient - random Fill(ran->Uniform(-fEta,fEta))
    fOutputList->Add(new TProfile("RDa2","RDa2;centrality;",8,centbins, ""));//second order legendre cofficient - random sample Fill(etaBG[cb]->GetRandom())
    if(fIsMC) {
      fOutputList->Add(new TProfile("MCa1","a1;centrality;",8,centbins, ""));//first order legendre cofficient MC - direct
      fOutputList->Add(new TProfile("MCa2","a2;centrality;",8,centbins, ""));//second order legendre cofficient MC - direct
    }
  }

  PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskLegendreCoef::UserExec(Option_t *)
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

  for(Int_t i(0); i < fAOD->GetNumberOfTracks(); i++) {                 // loop over all these tracks
    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
    if(!track) continue;
    if(fabs(track->Eta()) > fEta) continue;//eta cut
    if(track->Pt() < fPtmin|| track->Pt() > fPtmax) continue; //pt cut
    if(track->GetTPCNcls()<fTPCNcls || track->Chi2perNDF() > fChi2DoF) continue;// cut in TPC Ncls and chi2/dof   
    if(track->TestFilterBit(96)) {
      if(fIsBuildBG){
        //build background
        ((TH2D*) fOutputList->FindObject("ChargedBGHistOut"))->Fill(track->Eta(), Cent);
        if(track->Charge() > 0) ((TH2D*) fOutputList->FindObject("PosBGHistOut"))->Fill(track->Eta(), Cent);
        if(track->Charge() < 0) ((TH2D*) fOutputList->FindObject("NegBGHistOut"))->Fill(track->Eta(), Cent);
      }
      // (TODO)calculate legendre coefficients    
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
        if(fIsBuildBG){
          //build background
          ((TH2D*) fOutputList->FindObject("MCChargedBGHistOut"))->Fill(p1->Eta(), Cent);
          if(p1->Charge() > 0) ((TH2D*) fOutputList->FindObject("MCPosBGHistOut"))->Fill(p1->Eta(), Cent);
          if(p1->Charge() < 0) ((TH2D*) fOutputList->FindObject("MCNegBGHistOut"))->Fill(p1->Eta(), Cent);
        }
        // (TODO)calculate legendre coefficients    
      }
    }
  }
}

//_____________________________________________________________________________
void AliAnalysisTaskLegendreCoef::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
