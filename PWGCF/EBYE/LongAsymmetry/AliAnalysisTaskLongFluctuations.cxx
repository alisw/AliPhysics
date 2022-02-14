////////////////////////////////////////////////////////
// AliAnalysisTaskLongFluctuations:
// Description: Analysis task 
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
#include "AliAnalysisTaskLongFluctuations.h"
#include "AliAODMCParticle.h"
#include "AliEventCuts.h"
#include "AliAODMCHeader.h"


ClassImp(AliAnalysisTaskLongFluctuations)

AliAnalysisTaskLongFluctuations::AliAnalysisTaskLongFluctuations() : AliAnalysisTaskSE(),
fAOD(0), fOutputList(0), fNt(0), fIsMC(0), fChi2DoF(4), fTPCNcls(70), fPtmin(0.2), fPtmax(2), fEta(0.8), fEventCuts(0)

{
    mCentV0M = 0.;
    mCentCL0 = 0.;
    mCentCL1 = 0.;
    mPVz = 0.;
    mZDCN1 = 0.;
    mZDCN2 = 0.;
    mNTLs = 0.;
    mNGlob = 0.;
    mNHyb = 0.;
    mSpTGlob = 0.;
    mSpTHyb = 0.;
    mNMC = 0.;
    mSpTMC = 0.;
    
    for(int i=0;i<16;i++){
        aEtaPos[i] = 0.;
        aEtaNeg[i] = 0.;
        aEtaPosMC[i] = 0.;
        aEtaNegMC[i] = 0.;
    }
    
}
//_____________________________________________________________________________
AliAnalysisTaskLongFluctuations::AliAnalysisTaskLongFluctuations(const char* name) : AliAnalysisTaskSE(name),
    fAOD(0), fOutputList(0), fNt(0), fIsMC(0), fChi2DoF(4), fTPCNcls(70), fPtmin(0.2), fPtmax(2), fEta(0.8), fEventCuts(0)
{
    mCentV0M = 0.;
    mCentCL0 = 0.;
    mCentCL1 = 0.;
    mPVz = 0.;
    mZDCN1 = 0.;
    mZDCN2 = 0.;
    mNTLs = 0.;
    mNGlob = 0.;
    mNHyb = 0.;
    mSpTGlob = 0.;
    mSpTHyb = 0.;
    mNMC = 0.;
    mSpTMC = 0.;
    
    for(int i=0;i<16;i++){
        aEtaPos[i] = 0.;
        aEtaNeg[i] = 0.;
        aEtaPosMC[i] = 0.;
        aEtaNegMC[i] = 0.;
    }

    // constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskLongFluctuations::~AliAnalysisTaskLongFluctuations()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________

void AliAnalysisTaskLongFluctuations::SetMCRead(Bool_t bs){fIsMC=bs;}
void AliAnalysisTaskLongFluctuations::SetChi2DoF(Double_t Chi2DoF){fChi2DoF = Chi2DoF;}
void AliAnalysisTaskLongFluctuations::SetNclTPC(Int_t ncl){fTPCNcls = ncl;}
void AliAnalysisTaskLongFluctuations::SetPtLimits(Double_t ptmin, Double_t ptmax){fPtmin = ptmin; fPtmax = ptmax;}
void AliAnalysisTaskLongFluctuations::SetEtaLimit(Double_t etalimit){fEta = etalimit;}
void AliAnalysisTaskLongFluctuations::SetPileUpRead(Bool_t ps){fIsPileUpCuts=ps; }

void AliAnalysisTaskLongFluctuations::UserCreateOutputObjects()
{
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
    
    fNt = new TTree("nt","Long Fluc data");
    fNt->Branch("CentV0M",&mCentV0M);
    fNt->Branch("CentCL0",&mCentCL0);
    fNt->Branch("CentCL1",&mCentCL1);
    fNt->Branch("ZDCN1",&mZDCN1);
    fNt->Branch("ZDCN2",&mZDCN2);
    fNt->Branch("NTLs",&mNTLs);
    fNt->Branch("NGlob",&mNGlob);
    fNt->Branch("NHyb",&mNHyb);
    fNt->Branch("SpTGlob",&mSpTGlob);
    fNt->Branch("SpTHyb",&mSpTHyb);
    fNt->Branch("EtaPos",aEtaPos,"EtaPos[16]/F");
    fNt->Branch("EtaNeg",aEtaNeg,"EtaNeg[16]/F");
    
    if(fIsMC) {
        fNt->Branch("NMC",&mNMC);
        fNt->Branch("SpTMC",&mSpTMC);
        fNt->Branch("EtaPosMC",aEtaPosMC,"EtaPosMC[16]/F");
        fNt->Branch("EtaNegMC",aEtaNegMC,"EtaNegMC[16]/F");
    }

    fOutputList->Add(fNt);
    PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskLongFluctuations::UserExec(Option_t *)
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
    mPVz = PrimaryVertex->GetZ();
    if(fabs(mPVz)>8.0) return;
    
    AliMultSelection *MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
    if(!MultSelection) return;
    
    mCentV0M = MultSelection->GetMultiplicityPercentile("V0M"); //centrality
    mCentCL0 = MultSelection->GetMultiplicityPercentile("CL0");
    mCentCL1 = MultSelection->GetMultiplicityPercentile("CL1");
    
    mZDCN1 = fAOD->GetZDCN1Energy();
    mZDCN2 = fAOD->GetZDCN2Energy();
    mNTLs = fAOD->GetTracklets()->GetNumberOfTracklets();
    
    mNGlob = 0;
    mNHyb = 0;
    mSpTGlob = 0;
    mSpTHyb = 0;
    for(int i=0;i<16;i++){
        aEtaPos[i] = 0.;
        aEtaNeg[i] = 0.;
    }
    
    TH1D etaPos("etaPos","",16,-fEta,fEta);
    TH1D etaNeg("etaNeg","",16,-fEta,fEta);
    
    for(Int_t i(0); i < fAOD->GetNumberOfTracks(); i++) {                 // loop over all these tracks
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
        if(!track) continue;
        if(fabs(track->Eta()) > fEta) continue;//eta cut
        if(track->Pt() < fPtmin|| track->Pt() > fPtmax) continue; //pt cut
        if(track->GetTPCNcls()<fTPCNcls || track->Chi2perNDF() > fChi2DoF) continue;// cut in TPC Ncls and chi2/dof
        
        if(track->TestFilterBit(96)) {
            mNGlob++;
            mSpTGlob+=track->Pt();
            if(track->Charge() > 0) etaPos.Fill(track->Eta());
            if(track->Charge() < 0) etaNeg.Fill(track->Eta());
        }
        if(track->TestFilterBit(768)) {
            mNHyb++;
            mSpTHyb+=track->Pt();
        }
    }
    for(int i=0;i<16;i++){
        aEtaPos[i] = etaPos.GetBinContent(i+1);
        aEtaNeg[i] = etaNeg.GetBinContent(i+1);
    }
    
    if(fIsMC){
        TH1D etaPosMC("etaPosMC","",20,-fEta,fEta);
        TH1D etaNegMC("etaNegMC","",20,-fEta,fEta);
        
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



        mNMC  = 0;
        mSpTMC = 0;
        for(int i=0;i<16;i++){
            aEtaPosMC[i] = 0.;
            aEtaNegMC[i] = 0.;
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
                mNMC++;
                mSpTMC+=p1->Pt();
                if(p1->Charge() > 0) etaPosMC.Fill(p1->Eta());
                if(p1->Charge() < 0) etaNegMC.Fill(p1->Eta());
            }
        }
        for(int i=0;i<16;i++){
            aEtaPosMC[i] = etaPosMC.GetBinContent(i+1);
            aEtaNegMC[i] = etaNegMC.GetBinContent(i+1);
        }
    }
    fNt->Fill();
    PostData(1, fOutputList);
                                                      
}

//_____________________________________________________________________________
void AliAnalysisTaskLongFluctuations::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
