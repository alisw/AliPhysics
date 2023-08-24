////////////////////////////////////////////////////////
// AliAnalysisTaskLegendreCoef_local:
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
#include "AliAnalysisTaskLegendreCoef_local.h"
#include "AliAODMCParticle.h"
#include "AliEventCuts.h"
#include "AliAODMCHeader.h"
#include "TRandom2.h"
#include <TTreeStream.h>

ClassImp(AliAnalysisTaskLegendreCoef_local)

AliAnalysisTaskLegendreCoef_local::AliAnalysisTaskLegendreCoef_local() : AliAnalysisTaskSE(),
  fAOD(0), fOutputList(0),
  fIsMC(0), fChi2DoF(3), fTPCNcls(70), fPtmin(0.2), fPtmax(2), fEtaMin(-0.8), fEtaMax(0.8), fBit(96),
  fIsPileUpCuts(0), fIsBuildBG(0), fIsBuildLG(0), 
  fPosBackgroundHist(0), fNegBackgroundHist(0), fChargedBackgroundHist(0),
  fMCPosBackgroundHist(0), fMCNegBackgroundHist(0), fMCChargedBackgroundHist(0), fNeventCentHist(0),
  fGenName("Hijing"), fPileUpLevel(2), fTPCNCrossedRows(70), fPVzMax(8.0), fPVzMin(0.0), fPVzSign(0), fNetabins(16), PVz(0), Cent(0), fTreeSRedirector(0x0),fTreeMCrec(0x0),fTreeMCgen(0x0),fTreeRaw(0x0),fEventCuts(0)
{

}
//_____________________________________________________________________________
AliAnalysisTaskLegendreCoef_local::AliAnalysisTaskLegendreCoef_local(const char* name) : AliAnalysisTaskSE(name),
  fAOD(0), fOutputList(0),
  fIsMC(0), fChi2DoF(3), fTPCNcls(70), fPtmin(0.2), fPtmax(2), fEtaMin(-0.8), fEtaMax(0.8), fBit(96),
  fIsPileUpCuts(0), fIsBuildBG(0), fIsBuildLG(0), 
  fPosBackgroundHist(0), fNegBackgroundHist(0), fChargedBackgroundHist(0), 
  fMCPosBackgroundHist(0), fMCNegBackgroundHist(0), fMCChargedBackgroundHist(0), fNeventCentHist(0),
  fGenName("Hijing"), fPileUpLevel(2), fTPCNCrossedRows(70), fPVzMax(8.0), fPVzMin(0.0), fPVzSign(0), fNetabins(16), PVz(0), Cent(0), fTreeSRedirector(0x0),fTreeMCrec(0x0),fTreeMCgen(0x0),fTreeRaw(0x0),fEventCuts(0)
{
  // Default constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  // DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TList container
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
  DefineOutput(3, TTree::Class());
  DefineOutput(4, TTree::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskLegendreCoef_local::~AliAnalysisTaskLegendreCoef_local()
{
  // destructor
  if(fOutputList) {
    delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
  }
  if (fTreeSRedirector) delete fTreeSRedirector;

}
//_____________________________________________________________________________

void AliAnalysisTaskLegendreCoef_local::UserCreateOutputObjects()
{
  // Initialize output list of containers
  OpenFile(1);
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
    fOutputList->Add(new TH2D("PtCentHist", "track;pt;centrality", 50,fPtmin,fPtmax,8,centbins));//QA hist pt 
    if(fIsMC) {
      fOutputList->Add(new TH1D("GenNeventsCentHist", "nevents;centrality", 8,centbins));//output histogram when building background -  all tracks
      fOutputList->Add(new TH1D("RecNeventsCentHist", "nevents;centrality", 8,centbins));//output histogram when building background -  all tracks
      fOutputList->Add(new TH2D("MCPosBGHistOut", "postrack;eta;centrality", fNetabins,fEtaMin,fEtaMax,8,centbins));//output MC histogram when building background -  positive tracks
      fOutputList->Add(new TH2D("MCNegBGHistOut", "negtrack;eta;centrality", fNetabins,fEtaMin,fEtaMax,8,centbins));//output MC histogram when building background -  negative tracks
      fOutputList->Add(new TH2D("MCChargedBGHistOut", "track;eta;centrality", fNetabins,fEtaMin,fEtaMax,8,centbins));//output MC histogram when building background -  charged tracks  
      fOutputList->Add(new TH3D("GenPhiEtaCentHist", "track;phi;eta;centrality", 24, 0 ,TMath::Pi()*2.0, fNetabins,fEtaMin,fEtaMax, 100, 0.01, 70.0));//QA hist phi vs eta
      fOutputList->Add(new TH3D("RecPhiEtaCentHist", "track;phi;eta;centrality", 24, 0 ,TMath::Pi()*2.0, fNetabins,fEtaMin,fEtaMax, 100, 0.01, 70.0));//QA hist phi vs eta
      fOutputList->Add(new TH2D("GenPtCentHist", "track;pt;centrality", 50,fPtmin,fPtmax,8,centbins));//QA hist pt 
      fOutputList->Add(new TH2D("RecPtCentHist", "track;pt;centrality", 50,fPtmin,fPtmax,8,centbins));//QA hist pt
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

  // ************************************************************************
  //   Trees
  // ************************************************************************
  //
  fTreeSRedirector = new TTreeSRedirector();
  fTreeMCrec    = ((*fTreeSRedirector)<<"mcRec").GetTree();
  fTreeMCgen    = ((*fTreeSRedirector)<<"mcGen").GetTree();
  fTreeRaw    = ((*fTreeSRedirector)<<"Raw").GetTree();


  PostData(1, fOutputList);
  PostData(2, fTreeMCrec);
  PostData(3, fTreeMCgen);
  PostData(4, fTreeRaw);

}
//_____________________________________________________________________________
void AliAnalysisTaskLegendreCoef_local::BuildBackground()
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
  PVz = PrimaryVertex->GetZ();
  if(fabs(PVz)>fPVzMax) return;
  if(fabs(PVz)<fPVzMin) return;
  if(fPVzSign==1){//take only positive pvz
    if(PVz<0) return;
  }
  if(fPVzSign==-1){//take only negative pvz
    if(PVz>0) return;
  }

  AliMultSelection *MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
  if(!MultSelection) return;
  
  Cent = MultSelection->GetMultiplicityPercentile("V0M"); //centrality
  //mCentCL0 = MultSelection->GetMultiplicityPercentile("CL0");
  //mCentCL1 = MultSelection->GetMultiplicityPercentile("CL1");
  if(Cent>70.0 || Cent<0.01) return;
  ((TH1D*) fOutputList->FindObject("NeventsCentHist"))->Fill(Cent);//Nevents vs centrality

  if(!fTreeSRedirector) return;

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
      DumpMCGen(fAOD, p1);
      if(fabs(p1->Charge()<1))continue;//only get charged tracks
      if(!p1->IsPhysicalPrimary()) continue;
      if(p1->Eta()> fEtaMax || p1->Eta()<fEtaMin ) continue;
      if(p1->Pt() < fPtmin|| p1->Pt() > fPtmax) continue;
      if((fabs(p1->GetPdgCode())==211)||(fabs(p1->GetPdgCode())==2212)||(fabs(p1->GetPdgCode())==321)){
        //build background
        ((TH2D*) fOutputList->FindObject("MCChargedBGHistOut"))->Fill(p1->Eta(), Cent);
        if(p1->Charge() > 0) ((TH2D*) fOutputList->FindObject("MCPosBGHistOut"))->Fill(p1->Eta(), Cent);
        if(p1->Charge() < 0) ((TH2D*) fOutputList->FindObject("MCNegBGHistOut"))->Fill(p1->Eta(), Cent);
        ((TH3D*) fOutputList->FindObject("GenPhiEtaCentHist"))->Fill(p1->Phi(),p1->Eta(), Cent);
        ((TH2D*) fOutputList->FindObject("GenPtCentHist"))->Fill(p1->Pt(), Cent);
      }
    }
  }

  //fill hist nevent vs cent
  for(Int_t i(0); i < fAOD->GetNumberOfTracks(); i++) {                 // loop over all these tracks
    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
    if(!track) continue;
    if (fIsMC) {
      int label = TMath::Abs(track->GetLabel());
      AliAODMCParticle* mcTrack = dynamic_cast<AliAODMCParticle*>(stack->At(label));
      DumpMCRec(fAOD,track,mcTrack);
    }
    if(!fIsMC) DumpRaw(fAOD,track);

    if(fabs(track->Charge())<1) continue;//only get charged tracks
    if(track->Eta() > fEtaMax || track->Eta() < fEtaMin) continue;//eta cut
    if(track->Pt() < fPtmin|| track->Pt() > fPtmax) continue; //pt cut
    if(track->GetTPCNcls()<fTPCNcls || track->GetTPCNCrossedRows()<fTPCNCrossedRows || track->Chi2perNDF() > fChi2DoF) continue;// cut in TPC Ncls , crossed rows and chi2/dof  
    if(track->TestFilterBit(fBit)) {
      if(!fIsMC){
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
        if(fabs(mcTrack->Charge())<1)continue;//only get charged tracks
        // if(!mcTrack->IsPrimary()) continue;
        if(!mcTrack->IsPhysicalPrimary()) continue;  
        //
        // count Tracks 
        if((fabs(mcTrack->GetPdgCode())==211)||(fabs(mcTrack->GetPdgCode())==2212)||(fabs(mcTrack->GetPdgCode())==321)){
          ((TH2D*) fOutputList->FindObject("ChargedBGHistOut"))->Fill(mcTrack->Eta(), Cent);
          if(mcTrack->Charge() > 0) ((TH2D*) fOutputList->FindObject("PosBGHistOut"))->Fill(mcTrack->Eta(), Cent);
          if(mcTrack->Charge() < 0) ((TH2D*) fOutputList->FindObject("NegBGHistOut"))->Fill(mcTrack->Eta(), Cent);
          ((TH3D*) fOutputList->FindObject("RecPhiEtaCentHist"))->Fill(mcTrack->Phi(),mcTrack->Eta(), Cent);
          ((TH2D*) fOutputList->FindObject("RecPtCentHist"))->Fill(mcTrack->Pt(), Cent);
        }
      }
    }
  }
  return;
}
//_____________________________________________________________________________
void AliAnalysisTaskLegendreCoef_local::DumpMCRec(AliAODEvent *fAOD, AliAODTrack *trackAOD, AliAODMCParticle* mcTrack)
{

  //
  // Dump information for the track cut tuning
  Int_t nTPCClusters = fAOD->GetNumberOfTPCClusters();
  Int_t nITSClusters = 0;
  AliVMultiplicity *multiObj = fAOD->GetMultiplicity();
  for(Int_t i=2;i<6;i++) nITSClusters += multiObj->GetNumberOfITSClusters(i);
  //
  Float_t naodtracks = fAOD->GetNumberOfTracks();
  Float_t pt = trackAOD->Pt();
  Double_t chi2perndf = trackAOD->Chi2perNDF();
  UShort_t tpcncls = trackAOD->GetTPCncls();
  Double_t tpcchi2 = trackAOD->GetTPCchi2();
  Double_t tpcchi2percl = trackAOD->GetTPCchi2perCluster();
  Double_t tpcchi2perndf = trackAOD->GetTPCchi2perNDF();
  Double_t oldtpcchi2perndf = trackAOD->GetOldTPCchi2perNDF();
  Short_t charge = trackAOD->Charge();
  Bool_t istpconly = trackAOD->IsTPCOnly();
  Bool_t istpcconst = trackAOD->IsTPCConstrained();
  Bool_t isglobalconst = trackAOD->IsGlobalConstrained();
  UShort_t tpcncrossedrows = trackAOD->GetTPCNCrossedRows();
  Double_t ptot = trackAOD->GetTPCmomentum();
  Double_t tpctgl = trackAOD->GetTPCTgl();
  Double_t tpcsignaln = trackAOD->GetTPCsignalN();
  Bool_t bit96 = trackAOD->TestFilterBit(96);
  Bool_t bit128 = trackAOD->TestFilterBit(128);
  Bool_t bit16 = trackAOD->TestFilterBit(16);
  Bool_t bit768 = trackAOD->TestFilterBit(768);
  Bool_t bit512 = trackAOD->TestFilterBit(512);
  Double_t dEdx = trackAOD->GetTPCsignal();
  Float_t pv[2],cov[3];
  trackAOD->GetImpactParameters(pv,cov);
  Double_t leta = trackAOD->Eta();
  Double_t p = trackAOD->P();
  Double_t phi = trackAOD->Phi();
  Double_t theta = trackAOD->Theta();
  Float_t misscl = trackAOD->GetTPCClusterInfo(3,0,0,159);
  Float_t itsncls = trackAOD->GetITSNcls();
  Double_t itschi2 = trackAOD->GetITSchi2();
  //
  // Variables from as you use 
  Bool_t bPdgCond = ((fabs(mcTrack->GetPdgCode())==211)||(fabs(mcTrack->GetPdgCode())==2212)||(fabs(mcTrack->GetPdgCode())==321));
  Bool_t fTrackBitDef = trackAOD->TestFilterBit(fBit);
  Bool_t bEtaAcc = (trackAOD->Eta() > fEtaMax || trackAOD->Eta() < fEtaMin) ;
  Bool_t bPtAcc = (trackAOD->Pt() < fPtmin|| trackAOD->Pt() > fPtmax);
  Float_t fChi2Dof = trackAOD->Chi2perNDF();
  Int_t fNcls = trackAOD->GetTPCNcls();
  Int_t fNcrows = trackAOD->GetTPCNCrossedRows();
  Bool_t bPrim = (mcTrack->IsPrimary());
  Bool_t bPhysicalPrim = (mcTrack->IsPhysicalPrimary());
  Int_t fPDGcode = mcTrack->GetPdgCode();
  Int_t fCharge = mcTrack->Charge();
  //
  if (fTreeSRedirector){
    (*fTreeSRedirector)<<"mcRec"<<
    "tpcclmult=" << nTPCClusters <<  // eta
    "itsclmult=" << nITSClusters <<  // eta
    "naodtracks=" << naodtracks <<  // eta
    "dEdx=" << dEdx <<  // eta
    "phi=" << phi <<  // eta
    "theta=" << theta <<  // eta
    "ptot=" << ptot <<  // eta
    "pt=" << pt <<  // eta
    "p=" << p <<  // eta
    "eta=" << leta <<  // eta
    "dcaxy=" << pv[0] <<  // eta
    "dcaz=" << pv[1] <<  // eta
    "cov0=" << cov[0] <<  // eta
    "cov1=" << cov[1] <<  // eta
    "cov2=" << cov[2] <<  // eta
    "bit96=" << bit96 <<  // eta
    "bit128=" << bit128 <<  // eta
    "bit16=" << bit16 <<  //   eta
    "bit768=" << bit768 <<  // eta
    "bit512=" << bit512 <<  // eta
    "misscl=" << misscl <<  // eta
    "chi2perndf=" << chi2perndf <<  // eta
    "tpcncls=" << tpcncls <<  // eta
    "itsncls=" << itsncls <<  // eta
    "tpcchi2=" << tpcchi2 <<  // eta
    "itschi2=" << itschi2 <<  // eta
    "tpcchi2percl=" << tpcchi2percl <<  // real tpcchi2 which is used
    "tpcchi2perndf=" << tpcchi2perndf <<  // eta
    "oldtpcchi2perndf=" << oldtpcchi2perndf <<  // eta
    "charge=" << charge <<  // eta
    "istpconly=" << istpconly <<  // eta
    "istpcconst=" << istpcconst <<  // eta
    "isglobalconst=" << isglobalconst <<  // eta
    "tpcncrossedrows=" << tpcncrossedrows <<  // eta
    "tpctgl=" << tpctgl <<  // eta
    "tpcsignaln=" << tpcsignaln <<  // eta
    //
    "rBitDef=" << fTrackBitDef << 
    "rEtaAcc=" << bEtaAcc << 
    "rPtAcc=" << bPtAcc << 
    "rPrim=" << bPrim << 
    "rPhysicalPrim=" << bPhysicalPrim <<
    "rPDGcode=" << fPDGcode << 
    "rChi2Dof=" << fChi2Dof << 
    "rNcls=" <<fNcls << 
    "rNcrows=" << fNcrows <<
    "rCharge=" << fCharge << 
    "rPdgCond=" << bPdgCond << 
    "rcent="   << Cent <<     // Centrality
    "rvZ="     << PVz <<  //PVz
    "rPdgCond=" << bPdgCond << 
    "\n";
  }

}
//_____________________________________________________________________________
void AliAnalysisTaskLegendreCoef_local::DumpMCGen(AliAODEvent *fAOD, AliAODMCParticle* mcTrack)
{
  //
  // Variables from as you use 
  Bool_t bPdgCond = ((fabs(mcTrack->GetPdgCode())==211)||(fabs(mcTrack->GetPdgCode())==2212)||(fabs(mcTrack->GetPdgCode())==321));
  Bool_t bEtaAcc = (mcTrack->Eta() > fEtaMax || mcTrack->Eta() < fEtaMin) ;
  Bool_t bPtAcc = (mcTrack->Pt() < fPtmin|| mcTrack->Pt() > fPtmax);
  Bool_t bPrim = (mcTrack->IsPrimary());
  Bool_t bPhysicalPrim = (mcTrack->IsPhysicalPrimary());
  Int_t fPDGcode = mcTrack->GetPdgCode();
  Int_t fCharge = mcTrack->Charge();
  Bool_t fChargeCond = (mcTrack->Charge()!=-3 && mcTrack->Charge()!=+3);
  Float_t pTMCgen   = mcTrack->Pt();
  Float_t phiMCgen  = mcTrack->Phi();
  Float_t etaMCgen  = mcTrack->Eta();
  Float_t chargeGen = mcTrack->Eta();
  //
  if (fTreeSRedirector){
    (*fTreeSRedirector)<<"mcGen"<<
    "rEtaAcc=" << bEtaAcc << 
    "rPtAcc=" << bPtAcc << 
    "rPrim=" << bPrim << 
    "rPhysicalPrim=" << bPhysicalPrim <<
    "rPDGcode=" << fPDGcode << 
    "rCharge=" << fCharge << 
    "rPdgCond=" << bPdgCond << 
    "rcharge="  << chargeGen<<         // transverse momentum
    "rpT="      << pTMCgen<<           // transverse momentum
    "reta="     << etaMCgen <<          // mc eta
    "rphi="     << phiMCgen <<          // mc eta
    "rcent="    << Cent <<     // Centrality
    "rvZ="      << PVz <<  //PVz
    "\n";
  }

}
//_____________________________________________________________________________
void AliAnalysisTaskLegendreCoef_local::DumpRaw(AliAODEvent *fAOD, AliAODTrack *trackAOD)
{

  //
  // Dump information for the track cut tuning
  Int_t nTPCClusters = fAOD->GetNumberOfTPCClusters();
  Int_t nITSClusters = 0;
  AliVMultiplicity *multiObj = fAOD->GetMultiplicity();
  for(Int_t i=2;i<6;i++) nITSClusters += multiObj->GetNumberOfITSClusters(i);
  //
  Float_t naodtracks = fAOD->GetNumberOfTracks();
  Float_t pt = trackAOD->Pt();
  Double_t chi2perndf = trackAOD->Chi2perNDF();
  UShort_t tpcncls = trackAOD->GetTPCncls();
  Double_t tpcchi2 = trackAOD->GetTPCchi2();
  Double_t tpcchi2percl = trackAOD->GetTPCchi2perCluster();
  Double_t tpcchi2perndf = trackAOD->GetTPCchi2perNDF();
  Double_t oldtpcchi2perndf = trackAOD->GetOldTPCchi2perNDF();
  Short_t charge = trackAOD->Charge();
  Bool_t istpconly = trackAOD->IsTPCOnly();
  Bool_t istpcconst = trackAOD->IsTPCConstrained();
  Bool_t isglobalconst = trackAOD->IsGlobalConstrained();
  UShort_t tpcncrossedrows = trackAOD->GetTPCNCrossedRows();
  Double_t ptot = trackAOD->GetTPCmomentum();
  Double_t tpctgl = trackAOD->GetTPCTgl();
  Double_t tpcsignaln = trackAOD->GetTPCsignalN();
  Bool_t bit96 = trackAOD->TestFilterBit(96);
  Bool_t bit128 = trackAOD->TestFilterBit(128);
  Bool_t bit16 = trackAOD->TestFilterBit(16);
  Bool_t bit768 = trackAOD->TestFilterBit(768);
  Bool_t bit512 = trackAOD->TestFilterBit(512);
  Double_t dEdx = trackAOD->GetTPCsignal();
  Float_t pv[2],cov[3];
  trackAOD->GetImpactParameters(pv,cov);
  Double_t leta = trackAOD->Eta();
  Double_t p = trackAOD->P();
  Double_t phi = trackAOD->Phi();
  Double_t theta = trackAOD->Theta();
  Float_t misscl = trackAOD->GetTPCClusterInfo(3,0,0,159);
  Float_t itsncls = trackAOD->GetITSNcls();
  Double_t itschi2 = trackAOD->GetITSchi2();
  //
  // Variables from as you use 
  Bool_t fTrackBitDef = trackAOD->TestFilterBit(fBit);
  Bool_t bEtaAcc = (trackAOD->Eta() > fEtaMax || trackAOD->Eta() < fEtaMin) ;
  Bool_t bPtAcc = (trackAOD->Pt() < fPtmin|| trackAOD->Pt() > fPtmax);
  Float_t fChi2Dof = trackAOD->Chi2perNDF();
  Int_t fNcls = trackAOD->GetTPCNcls();
  Int_t fNcrows = trackAOD->GetTPCNCrossedRows();
  //
  if (fTreeSRedirector){
    (*fTreeSRedirector)<<"Raw"<<
    "tpcclmult=" << nTPCClusters <<  // eta
    "itsclmult=" << nITSClusters <<  // eta
    "naodtracks=" << naodtracks <<  // eta
    "dEdx=" << dEdx <<  // eta
    "phi=" << phi <<  // eta
    "theta=" << theta <<  // eta
    "ptot=" << ptot <<  // eta
    "pt=" << pt <<  // eta
    "p=" << p <<  // eta
    "eta=" << leta <<  // eta
    "dcaxy=" << pv[0] <<  // eta
    "dcaz=" << pv[1] <<  // eta
    "cov0=" << cov[0] <<  // eta
    "cov1=" << cov[1] <<  // eta
    "cov2=" << cov[2] <<  // eta
    "bit96=" << bit96 <<  // eta
    "bit128=" << bit128 <<  // eta
    "bit16=" << bit16 <<  //   eta
    "bit768=" << bit768 <<  // eta
    "bit512=" << bit512 <<  // eta
    "misscl=" << misscl <<  // eta
    "chi2perndf=" << chi2perndf <<  // eta
    "tpcncls=" << tpcncls <<  // eta
    "itsncls=" << itsncls <<  // eta
    "tpcchi2=" << tpcchi2 <<  // eta
    "itschi2=" << itschi2 <<  // eta
    "tpcchi2percl=" << tpcchi2percl <<  // real tpcchi2 which is used
    "tpcchi2perndf=" << tpcchi2perndf <<  // eta
    "oldtpcchi2perndf=" << oldtpcchi2perndf <<  // eta
    "charge=" << charge <<  // eta
    "istpconly=" << istpconly <<  // eta
    "istpcconst=" << istpcconst <<  // eta
    "isglobalconst=" << isglobalconst <<  // eta
    "tpcncrossedrows=" << tpcncrossedrows <<  // eta
    "tpctgl=" << tpctgl <<  // eta
    "tpcsignaln=" << tpcsignaln <<  // eta
    //
    "rBitDef=" << fTrackBitDef << 
    "rEtaAcc=" << bEtaAcc << 
    "rPtAcc=" << bPtAcc << 
    "rChi2Dof=" << fChi2Dof << 
    "rNcls=" <<fNcls << 
    "rNcrows=" << fNcrows <<
    "rcent="   << Cent <<     // Centrality
    "rvZ="     << PVz <<  //PVz
    "\n";
  }

}
//_____________________________________________________________________________
void AliAnalysisTaskLegendreCoef_local::BuildSignal()
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
  if(fPVzSign==1){//take only positive pvz
    if(PVz<0) return;
  }
  if(fPVzSign==-1){//take only negative pvz
    if(PVz>0) return;
  }
  
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
      if(p1->Charge()==0)continue;//only get charged tracks
      if(p1->Charge()!=-3 && p1->Charge()!=+3) continue;// x3 by convention
      if(!p1->IsPrimary()) continue;
      if(!p1->IsPhysicalPrimary()) continue;
      if(p1->Eta() > fEtaMax || p1->Eta()<fEtaMin ) continue;
      if(p1->Pt() < fPtmin|| p1->Pt() > fPtmax) continue;
      if((fabs(p1->GetPdgCode())==211)||(fabs(p1->GetPdgCode())==2212)||(fabs(p1->GetPdgCode())==321)){
        // calculate signal    
        MCchargedSignal->Fill(p1->Eta());
        if(p1->Charge() > 0) MCposSignal->Fill(p1->Eta());
        if(p1->Charge() < 0) MCnegSignal->Fill(p1->Eta());        
      }
    }


  for(Int_t i(0); i < fAOD->GetNumberOfTracks(); i++) {                 // loop over all these tracks
    AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
    if(!track) continue;
    if(track->Charge()==0)continue;//only get charged tracks
    if(track->Eta() > fEtaMax || track->Eta() < fEtaMin) continue;//eta cut
    if(track->Pt() < fPtmin|| track->Pt() > fPtmax) continue; //pt cut
    if(track->GetTPCNcls()<fTPCNcls || track->GetTPCNCrossedRows()<fTPCNCrossedRows || track->Chi2perNDF() > fChi2DoF) continue;// cut in TPC Ncls, crossed rows and chi2/dof   
    if(track->TestFilterBit(fBit)) {
      if(!fIsMC){
        //build signal for raw data
        chargedSignal->Fill(track->Eta());
        if(track->Charge() > 0) posSignal->Fill(track->Eta());
        if(track->Charge() < 0) negSignal->Fill(track->Eta());
      }else{
        //build signal for MC reconstructed
        int label = TMath::Abs(track->GetLabel());
        AliAODMCParticle* mcTrack = dynamic_cast<AliAODMCParticle*>(stack->At(label));
        if (!mcTrack) continue;
        if(mcTrack->Charge()==0)continue;//only get charged tracks
        if(!mcTrack->IsPrimary()) continue;
        if(!mcTrack->IsPhysicalPrimary()) continue;    
        if((fabs(mcTrack->GetPdgCode())==211)||(fabs(mcTrack->GetPdgCode())==2212)||(fabs(mcTrack->GetPdgCode())==321)){
          chargedSignal->Fill(mcTrack->Eta());
          if(mcTrack->Charge() > 0) posSignal->Fill(mcTrack->Eta());
          if(mcTrack->Charge() < 0) negSignal->Fill(mcTrack->Eta());
        }
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
    if(MCchargedSignal->Integral()>0) {
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
void AliAnalysisTaskLegendreCoef_local::BuildCoefficients(TH1D *signal, TH1D *background, Float_t centrality, TString type)
{
  TH1D *RanHist[5];//random distribution histograms
  int ntracks = signal->Integral();
  int n; 
  char histname[50];
    printf("tracks is %i\n",ntracks);

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
void AliAnalysisTaskLegendreCoef_local::UserExec(Option_t *)
{
  if(fIsBuildBG) BuildBackground();
  if(fIsBuildLG) BuildSignal();
}
//_____________________________________________________________________________
Double_t AliAnalysisTaskLegendreCoef_local::GetSingleAnCoef(int order, TH1D *hist)
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
Double_t AliAnalysisTaskLegendreCoef_local::LegPol(int order, Double_t x){
    double etawidth = (fEtaMax-fEtaMin)/2.0;
    if( order == 1 ) return x/etawidth;
    if( order == 2 ) return 0.5*(3.0*pow(x/etawidth,2)-1.0);
    if( order == 3 ) return 0.5*(5.0*pow(x/etawidth,3)-3.0*x/etawidth);
    if( order == 4 ) return (35.0*pow(x/etawidth,4)-30.0*pow(x/etawidth,2)+3.0)/8.0;
    else return 0.0;
}
//_____________________________________________________________________________
void AliAnalysisTaskLegendreCoef_local::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________