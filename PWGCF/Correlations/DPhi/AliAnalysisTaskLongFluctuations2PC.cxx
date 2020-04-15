////////////////////////////////////////////////////////
// AliAnalysisTaskLongFluctuations2PC:
// Description: Analysis task for 2PC for eta1, eta2
// Author: Raquel Quishpe (raquel.quishpe@cern.ch)
////////////////////////////////////////////////////////

#include "TChain.h"
#include "TH3D.h"
#include "TProfile.h"
#include "TString.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAODTrack.h"
#include "AliMultSelection.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskLongFluctuations2PC.h"
#include "AliAODMCParticle.h"

ClassImp(AliAnalysisTaskLongFluctuations2PC)

AliAnalysisTaskLongFluctuations2PC::AliAnalysisTaskLongFluctuations2PC() : AliAnalysisTaskSE(),
    fAOD(0), fOutputList(0), fMultD(0), fMultDpT(0), fMultDMC(0), fMultDpTMC(0), fMultDReMC(0), fPVzCentNevents(0), fCentMultQA(0), fAllCentQA(0), fIsMC(0), fFB(0), fCentrality("V0M"), fChi2DoF(0), fTPCNcls(0), fPtmin(0), fPtmax(0), fEta(0)
{
    for(Int_t i(0);i<9;i++) {
        fHistSigMC[i]=0; fHistBgMC[i]=0; fHistSigGapMC[i]=0; fHistBgGapMC[i]=0;
        for(Int_t j(0);j<17;j++) {fHistSig[i][j]=0; fHistBg[i][j]=0; fHistSigGap[i][j]=0; fHistBgGap[i][j]=0;}
    }
}
//_____________________________________________________________________________
AliAnalysisTaskLongFluctuations2PC::AliAnalysisTaskLongFluctuations2PC(const char* name) : AliAnalysisTaskSE(name),
    fAOD(0), fOutputList(0), fMultD(0), fMultDpT(0), fMultDMC(0), fMultDpTMC(0), fMultDReMC(0), fPVzCentNevents(0), fCentMultQA(0), fAllCentQA(0), fIsMC(0), fFB(0), fCentrality("V0M"), fChi2DoF(0), fTPCNcls(0), fPtmin(0), fPtmax(0), fEta(0)
{
    for(Int_t i(0);i<9;i++) {
        fHistSigMC[i]=0; fHistBgMC[i]=0; fHistSigGapMC[i]=0; fHistBgGapMC[i]=0;
        for(Int_t j(0);j<17;j++) {fHistSig[i][j]=0; fHistBg[i][j]=0; fHistSigGap[i][j]=0; fHistBgGap[i][j]=0;}
    }
    // constructor
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}
//_____________________________________________________________________________
AliAnalysisTaskLongFluctuations2PC::~AliAnalysisTaskLongFluctuations2PC()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
}
//_____________________________________________________________________________

void AliAnalysisTaskLongFluctuations2PC::SetMCRead(Bool_t bs){fIsMC=bs;}
void AliAnalysisTaskLongFluctuations2PC::SetFilterBit(Int_t fb){fFB = fb;}
void AliAnalysisTaskLongFluctuations2PC::SetCentrality(TString Cent){fCentrality = Cent;}
void AliAnalysisTaskLongFluctuations2PC::SetChi2DoF(Double_t Chi2DoF){fChi2DoF = Chi2DoF;}
void AliAnalysisTaskLongFluctuations2PC::SetNclTPC(Int_t ncl){fTPCNcls = ncl;}
void AliAnalysisTaskLongFluctuations2PC::SetPtLimits(Double_t ptmin, Double_t ptmax){fPtmin = ptmin; fPtmax = ptmax;}
void AliAnalysisTaskLongFluctuations2PC::SetEtaLimit(Double_t etalimit){fEta = etalimit;}


void AliAnalysisTaskLongFluctuations2PC::UserCreateOutputObjects()
{
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);

    Double_t centbins[9] = {0.,5,10,20,30,40,50,60,70};
    Double_t pvzbins[17] = {-8.,-7.,-6.,-5.,-4.,-3.,-2.,-1.,0.,1.,2.,3.,4.,5.,6.,7.,8.};
    Int_t etaBins = 40;
    char hname[20];
    
    fMultD = new TH1D("multD","",4000,-0.5,3999.5);
    fMultD->SetXTitle("N_{ch}");
    fMultDpT = new TProfile("multDpT","",4000,-0.5,3999.5);
    fMultDpT->SetXTitle("N_{ch}");
    fMultDpT->SetYTitle("#langle p_{T} #rangle");
    fPVzCentNevents = new TH2D("PVzCentNevents","",16,pvzbins,8,centbins);
    fPVzCentNevents->SetXTitle("PVz");
    fPVzCentNevents->SetYTitle("Cent");
    fCentMultQA = new TH2D("CentMultQA","",100,0.0,100.0,4000,-0.5,3999.5);
    fCentMultQA->SetXTitle("Cent V0M");
    fCentMultQA->SetYTitle("Mult");
    fAllCentQA = new TH3D("AllCentQA","",100,0,100,100,0,100,100,0,100);
    fAllCentQA->SetXTitle("Cent V0M");
    fAllCentQA->SetYTitle("Cent CL0");
    fAllCentQA->SetZTitle("Cent TPC");
    
    fOutputList->Add(fMultD);
    fOutputList->Add(fMultDpT);
    fOutputList->Add(fPVzCentNevents);
    fOutputList->Add(fCentMultQA);
    fOutputList->Add(fAllCentQA);

    if(fIsMC) {
        fMultDMC = new TH1D("multDMC","",4000,-0.5,3999.5);
        fMultDMC->SetXTitle("N_{ch}");
        fMultDpTMC = new TProfile("multDpTMC","",4000,-0.5,3999.5);
        fMultDpTMC->SetXTitle("N_{ch}");
        fMultDpTMC->SetYTitle("#langle p_{T} #rangle");
        fOutputList->Add(fMultDMC);
        fOutputList->Add(fMultDpTMC);
        fMultDReMC = new TH2D("multDReMC","",4000,-0.5,3999.5,4000,-0.5,3999.5);
        fMultDReMC->SetXTitle("Reco N_{ch}");
        fMultDReMC->SetYTitle("Gen N_{ch}");
        fOutputList->Add(fMultDReMC);

    }
    for(Int_t i(0);i<9;i++){
        for(Int_t j(0);j<17;j++){
            sprintf(hname,"histSig_%i_%i",i,j);
            fHistSig[i][j] = new TH2D(hname,"",etaBins,-0.8,0.8,etaBins, -0.8,0.8);
            fHistSig[i][j]->SetXTitle("#eta_{1}");
            fHistSig[i][j]->SetYTitle("#eta_{2}");
            fOutputList->Add(fHistSig[i][j]);

            sprintf(hname,"histBg_%i_%i",i,j);
            fHistBg[i][j] = new TH1D(hname,"",etaBins,-0.8,0.8);
            fHistBg[i][j]->SetXTitle("#eta");
            fOutputList->Add(fHistBg[i][j]);

            sprintf(hname,"histSigGap_%i_%i",i,j);
            fHistSigGap[i][j] = new TH2D(hname,"",etaBins,-0.8,0.8,etaBins,-0.8,0.8);
            fHistSigGap[i][j]->SetXTitle("#eta_{1}");
            fHistSigGap[i][j]->SetYTitle("#eta_{2}");
            fOutputList->Add(fHistSigGap[i][j]);

            sprintf(hname,"histBgGap_%i_%i",i,j);
            fHistBgGap[i][j] = new TH2D(hname,"",etaBins,-0.8,0.8,25, 0.0, 2.0*PI);
            fHistBgGap[i][j]->SetXTitle("#eta");
            fHistBgGap[i][j]->SetYTitle("#phi");           
            fOutputList->Add(fHistBgGap[i][j]);

        }
        if(fIsMC) {
            sprintf(hname,"histSigMC_%i",i);
            fHistSigMC[i] = new TH2D(hname,"",etaBins,-0.8,0.8,etaBins,-0.8,0.8);
            fHistSigMC[i]->SetXTitle("#eta_{1}");
            fHistSigMC[i]->SetYTitle("#eta_{2}");
            fOutputList->Add(fHistSigMC[i]);

            sprintf(hname,"histBgMC_%i",i);
            fHistBgMC[i] = new TH1D(hname,"",etaBins,-0.8,0.8);
            fHistBgMC[i]->SetXTitle("#eta");
            fOutputList->Add(fHistBgMC[i]);

            sprintf(hname,"histSigGapMC_%i",i);
            fHistSigGapMC[i] = new TH2D(hname,"",etaBins,-0.8,0.8,etaBins,-0.8,0.8);
            fHistSigGapMC[i]->SetXTitle("#eta_{1}");
            fHistSigGapMC[i]->SetYTitle("#eta_{2}");
            fOutputList->Add(fHistSigGapMC[i]);

            sprintf(hname,"histBgGapMC_%i",i);
            fHistBgGapMC[i] = new TH2D(hname,"",etaBins,-0.8,0.8,25, 0.0, 2.0*PI);
            fHistBgGapMC[i]->SetXTitle("#eta");
            fHistBgGapMC[i]->SetYTitle("#phi");
            fOutputList->Add(fHistBgGapMC[i]);
        }
    }
    
    PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskLongFluctuations2PC::UserExec(Option_t *)
{
        
    Double_t centbins[9] = {0.,5,10,20,30,40,50,60,70};
    Double_t pvzbins[17] = {-8.,-7.,-6.,-5.,-4.,-3.,-2.,-1.,0.,1.,2.,3.,4.,5.,6.,7.,8.};
    Float_t mPVz, mCent;
    Double_t eta1, eta2, phi1, phi2;
    TAxis *centaxis = new TAxis(8, centbins);
    TAxis *pvzaxis = new TAxis(16, pvzbins);

    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) return;
    
    //making a cut in pvz -8 to 8cm
    const AliAODVertex *PrimaryVertex = fAOD->GetVertex(0);
    if(!PrimaryVertex) return;
    mPVz = PrimaryVertex->GetZ();
    if(fabs(mPVz)<0.000001 || fabs(mPVz)>8.0) return;
    
    AliMultSelection *MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
    if(!MultSelection) return;
    mCent = MultSelection->GetMultiplicityPercentile(fCentrality); //centrality
    if(mCent < 0.000001) return;
    if(fAOD->IsPileupFromSPD(20)) return;
    
    fAllCentQA->Fill(MultSelection->GetMultiplicityPercentile("V0M"),MultSelection->GetMultiplicityPercentile("CL0"),MultSelection->GetMultiplicityPercentile("CL1"));
    fPVzCentNevents->Fill(mPVz,mCent);
    Int_t bCent(centaxis->FindBin(mCent)-1);
    Int_t bPVz(pvzaxis->FindBin(mPVz)-1);
    Int_t nTracks(fAOD->GetNumberOfTracks());           // see how many tracks there are in the event
    
    Double_t nAcc = 0.0;
    Double_t mpT = 0.0;
        
    for(Int_t i(0); i < nTracks; i++) {                 // loop over all these tracks
        AliAODTrack* track1 = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
        if(!track1 || !track1->TestFilterBit(fFB)) continue;
        if(fabs(track1->Eta()) > fEta) continue;//eta cut
        if(track1->Pt() < fPtmin|| track1->Pt() > fPtmax) continue; //pt cut
        if(track1->GetTPCNcls()<fTPCNcls || track1->Chi2perNDF() > fChi2DoF) continue;// cut in TPC Ncls and chi2/dof
        eta1 = track1->Eta();
        phi1 = track1->Phi();
        fHistBg[bCent][bPVz]->Fill(eta1);
        fHistBgGap[bCent][bPVz]->Fill(eta1, phi1);
        nAcc+=1.0;
        mpT+= track1->Pt();
        for(Int_t j(0); j < nTracks; j++) {
            AliAODTrack* track2 = static_cast<AliAODTrack*>(fAOD->GetTrack(j));
            if(!track2 || !track2->TestFilterBit(fFB)) continue;
            if(fabs(track2->Eta()) > fEta) continue;//eta cut
            if(track2->Pt() < fPtmin|| track2->Pt() > fPtmax) continue; //pt cut
            if(track2->GetTPCNcls()<fTPCNcls || track2->Chi2perNDF() > fChi2DoF) continue;//cut in TPC Ncls and chi2/dof
            eta2 = track2->Eta();
            phi2 = track2->Phi();
            fHistSig[bCent][bPVz]->Fill(eta1,eta2);
            if(fabs (phi1-phi2)>PI/2. && fabs(phi1-phi2)<3.*PI/2.){//gap signal histogram
                fHistSigGap[bCent][bPVz]->Fill(eta1, eta2);
            }
        }
    }
    if (nAcc > 0.001) {
        mpT/=nAcc;
        fMultD->Fill(nAcc);
        fMultDpT->Fill(nAcc,mpT);
        fCentMultQA->Fill(mCent,nAcc);
    }
   
    if(fIsMC){
        Double_t nAccMC = 0.0;
        mpT = 0.0;
        TClonesArray *stack = 0;
        TList *lst = fAOD->GetList();
        stack = (TClonesArray*)lst->FindObject(AliAODMCParticle::StdBranchName());
        int nMCTracks;
        if (!stack) nMCTracks = 0;
        else nMCTracks = stack->GetEntries();
        for (Int_t i(0); i < nMCTracks; i++) {
             AliAODMCParticle *p1=(AliAODMCParticle*)stack->UncheckedAt(i);
             if (!p1) continue;
             if(p1->Charge()!=-3 && p1->Charge()!=+3) continue;// x3 by convention
             if(!p1->IsPrimary()) continue;
             if(!p1->IsPhysicalPrimary()) continue;
             if(fabs(p1->Eta()) > fEta ) continue;
             if(p1->Pt() < fPtmin|| p1->Pt() > fPtmax) continue;
             if((fabs(p1->GetPdgCode())==211)||(fabs(p1->GetPdgCode())==2212)||(fabs(p1->GetPdgCode())==321)){
                 eta1 = p1->Eta();
                 phi1 = p1->Phi();
                 fHistBgMC[bCent]->Fill(eta1);
                 fHistBgGapMC[bCent]->Fill(eta1, phi1);
                 nAccMC+=1.0;
                 mpT+= p1->Pt();
                 for (Int_t j(0); j < nMCTracks; j++) {
                       AliAODMCParticle *p2=(AliAODMCParticle*)stack->UncheckedAt(j);
                       if (!p2) continue;
                       if(p2->Charge()!=-3 && p2->Charge()!=+3) continue;// x3 by convention
                       if(!p2->IsPrimary()) continue;
                       if(!p2->IsPhysicalPrimary()) continue;
                       if(fabs(p2->Eta()) > fEta ) continue;
                       if(p2->Pt() < fPtmin|| p2->Pt() > fPtmax) continue;
                       if((fabs(p2->GetPdgCode())==211)||(fabs(p2->GetPdgCode())==2212)||(fabs(p2->GetPdgCode())==321)){
                           eta2 = p2->Eta();
                           phi2 = p2->Phi();
                           fHistSigMC[bCent]->Fill(eta1,eta2);
                             if(fabs (phi1-phi2)>PI/2. && fabs(phi1-phi2)<3.*PI/2.){//gap signal histogram
                                 fHistSigGapMC[bCent]->Fill(eta1, eta2);
                             }
                 
                       }
                 }
             }
        }
        if (nAcc > 0.001 && nAccMC > 0.001) {
            mpT/=nAccMC;
            fMultDMC->Fill(nAccMC);
            fMultDpTMC->Fill(nAccMC,mpT);
            fMultDReMC->Fill(nAcc,nAccMC);
        }
    }

    PostData(1, fOutputList);
                                                      
}
//_____________________________________________________________________________
void AliAnalysisTaskLongFluctuations2PC::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
