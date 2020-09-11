////////////////////////////////////////////////////////
// AliAnalysisTaskLongFluctuations:
// Description: Analysis task 
// Author: Raquel Quishpe (raquel.quishpe@cern.ch)
////////////////////////////////////////////////////////

#include "TChain.h"
#include "TH3D.h"
#include "TNtuple.h"
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

ClassImp(AliAnalysisTaskLongFluctuations)

AliAnalysisTaskLongFluctuations::AliAnalysisTaskLongFluctuations() : AliAnalysisTaskSE(),
fAOD(0), fOutputList(0), fNt(0), fIsMC(0), fChi2DoF(4), fTPCNcls(70), fPtmin(0.2), fPtmax(2), fEta(0.8),fEtaGlobCentPVz(0), fEtaTPCCentPVz(0), fEtaHybCentPVz(0), fEtaMCCent(0)

{

}
//_____________________________________________________________________________
AliAnalysisTaskLongFluctuations::AliAnalysisTaskLongFluctuations(const char* name) : AliAnalysisTaskSE(name),
    fAOD(0), fOutputList(0), fNt(0), fIsMC(0), fChi2DoF(4), fTPCNcls(70), fPtmin(0.2), fPtmax(2), fEta(0.8),fEtaGlobCentPVz(0), fEtaTPCCentPVz(0), fEtaHybCentPVz(0), fEtaMCCent(0)
{
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

void AliAnalysisTaskLongFluctuations::UserCreateOutputObjects()
{
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
    
    fEtaGlobCentPVz = new TH3D("EtaGlobCentPVz","",20,-(Double_t)fEta,(Double_t)fEta,100,0.0,100.0,16,-8.0,8.0);
    fEtaTPCCentPVz = new TH3D("EtaTPCCentPVz","",20,-(Double_t)fEta,(Double_t)fEta,100,0,100,16,-8.0,8.0);
    fEtaHybCentPVz = new TH3D("EtaHybCentPVz","",20,-(Double_t)fEta,(Double_t)fEta,100,0.0,100.0,16,-8.0,8.0);
    fOutputList->Add(fEtaGlobCentPVz);
    fOutputList->Add(fEtaTPCCentPVz);
    fOutputList->Add(fEtaHybCentPVz);
    
    if(fIsMC) {fNt = new TNtuple("nt","LongFluctuations","CentV0M:CentCL0:CentCL1:PVz:ZDCN1:ZDCP1:ZDCN2:ZDCP2:NTLs:NGlob:NTPC:NHyb:SpTGlob:SpTTPC:SpTHyb:a1Glob:a2Glob:a3Glob:a4Glob:a1TPC:a2TPC:a3TPC:a4TPC:a1Hyb:a2Hyb:a3Hyb:a4Hyb:a1pos:a2pos:a3pos:a4pos:a1neg:a2neg:a3neg:a4neg:NMC:SpTMC:a1MC:a2MC:a3MC:a4MC:a1posMC:a2posMC:a3posMC:a4posMC:a1negMC:a2negMC:a3negMC:a4negMC");
        fEtaMCCent = new TH2D("EtaMCCent","",20,-(Double_t)fEta,(Double_t)fEta,100,0.0,100.0);
        fOutputList->Add(fEtaMCCent);
    }
    else fNt = new TNtuple("nt","LongFluctuations","CentV0M:CentCL0:CentCL1:PVz:ZDCN1:ZDCP1:ZDCN2:ZDCP2:NTLs:NGlob:NTPC:NHyb:SpTGlob:SpTTPC:SpTHyb:a1Glob:a2Glob:a3Glob:a4Glob:a1TPC:a2TPC:a3TPC:a4TPC:a1Hyb:a2Hyb:a3Hyb:a4Hyb:a1pos:a2pos:a3pos:a4pos:a1neg:a2neg:a3neg:a4neg");
    fOutputList->Add(fNt);
    PostData(1, fOutputList);
}
//_____________________________________________________________________________
void AliAnalysisTaskLongFluctuations::UserExec(Option_t *)
{
        
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) return;

    Float_t CentV0M, CentCL0, CentCL1, PVz, ZDCN1, ZDCP1, ZDCN2 ,ZDCP2, NTLs, NGlob, NTPC, NHyb, SpTGlob, SpTTPC, SpTHyb, a1Glob, a2Glob, a3Glob, a4Glob ,a1TPC, a2TPC, a3TPC, a4TPC, a1Hyb, a2Hyb, a3Hyb, a4Hyb, a1pos, a2pos, a3pos, a4pos, a1neg, a2neg ,a3neg, a4neg;
    
    //making a cut in pvz -8 to 8cm
    const AliAODVertex *PrimaryVertex = fAOD->GetVertex(0);
    if(!PrimaryVertex) return;
    PVz = PrimaryVertex->GetZ();
    if(fabs(PVz)>8.0) return;
    
    AliMultSelection *MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
    if(!MultSelection) return;
    
    CentV0M = MultSelection->GetMultiplicityPercentile("V0M"); //centrality
    CentCL0 = MultSelection->GetMultiplicityPercentile("CL0");
    CentCL1 = MultSelection->GetMultiplicityPercentile("CL1");
    
    ZDCN1 = fAOD->GetZDCN1Energy();
    ZDCP1 = fAOD->GetZDCP1Energy();
    ZDCN2 = fAOD->GetZDCN2Energy();
    ZDCP2 = fAOD->GetZDCP2Energy();
    NTLs = fAOD->GetTracklets()->GetNumberOfTracklets();
    
    NGlob = 0;
    NTPC  = 0;
    NHyb = 0;
    SpTGlob = 0;
    SpTTPC = 0;
    SpTHyb = 0;
    
    TH1D etaTPC("etaTPC","",20,-fEta,fEta);
    TH1D etaGlob("etaGlob","",20,-fEta,fEta);
    TH1D etaHyb("etaHyb","",20,-fEta,fEta);
    TH1D etaPos("etaPos","",20,-fEta,fEta);
    TH1D etaNeg("etaNeg","",20,-fEta,fEta);
    
    for(Int_t i(0); i < fAOD->GetNumberOfTracks(); i++) {                 // loop over all these tracks
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));         // get a track (type AliAODTrack) from the event
        if(!track) continue;
        if(fabs(track->Eta()) > fEta) continue;//eta cut
        if(track->Pt() < fPtmin|| track->Pt() > fPtmax) continue; //pt cut
        if(track->GetTPCNcls()<fTPCNcls || track->Chi2perNDF() > fChi2DoF) continue;// cut in TPC Ncls and chi2/dof
        
        if(track->TestFilterBit(96)) {
            NGlob++;
            SpTGlob+=track->Pt();
            etaGlob.Fill(track->Eta());
            fEtaGlobCentPVz->Fill(track->Eta(),CentV0M,PVz);
        }
        if(track->TestFilterBit(128)) {
            NTPC++;
            SpTTPC+=track->Pt();
            etaTPC.Fill(track->Eta());
            fEtaTPCCentPVz->Fill(track->Eta(),CentV0M,PVz);
        }
        if(track->TestFilterBit(768)) {
            NHyb++;
            SpTHyb+=track->Pt();
            etaHyb.Fill(track->Eta());
            if(track->Charge() > 0) etaPos.Fill(track->Eta());
            if(track->Charge() < 0) etaNeg.Fill(track->Eta());
            fEtaHybCentPVz->Fill(track->Eta(),CentV0M,PVz);
        }
    }
    
    
    if(etaGlob.Integral() < 1 || etaTPC.Integral() < 1 || etaHyb.Integral() < 1 || etaPos.Integral() < 1 || etaNeg.Integral() < 1) return;
    etaGlob.Scale(1.0/(etaGlob.Integral()/20.));
    etaTPC.Scale(1.0/(etaTPC.Integral()/20.));
    etaHyb.Scale(1.0/(etaHyb.Integral()/20.));
    etaPos.Scale(1.0/(etaPos.Integral()/20.));
    etaNeg.Scale(1.0/(etaNeg.Integral()/20.));
    
    a1Glob = GetAnCoeff(1,&etaGlob);
    a2Glob = GetAnCoeff(2,&etaGlob);
    a3Glob = GetAnCoeff(3,&etaGlob);
    a4Glob = GetAnCoeff(4,&etaGlob);

    a1TPC = GetAnCoeff(1,&etaTPC);
    a2TPC = GetAnCoeff(2,&etaTPC);
    a3TPC = GetAnCoeff(3,&etaTPC);
    a4TPC = GetAnCoeff(4,&etaTPC);

    a1Hyb = GetAnCoeff(1,&etaHyb);
    a2Hyb = GetAnCoeff(2,&etaHyb);
    a3Hyb = GetAnCoeff(3,&etaHyb);
    a4Hyb = GetAnCoeff(4,&etaHyb);

    a1pos = GetAnCoeff(1,&etaPos);
    a2pos = GetAnCoeff(2,&etaPos);
    a3pos = GetAnCoeff(3,&etaPos);
    a4pos = GetAnCoeff(4,&etaPos);
    
    a1neg = GetAnCoeff(1,&etaNeg);
    a2neg = GetAnCoeff(2,&etaNeg);
    a3neg = GetAnCoeff(3,&etaNeg);
    a4neg = GetAnCoeff(4,&etaNeg);
    
    if(fIsMC){
        TH1D etaMC("etaMC","",20,-fEta,fEta);
        TH1D etaPosMC("etaPosMC","",20,-fEta,fEta);
        TH1D etaNegMC("etaNegMC","",20,-fEta,fEta);
        
        Float_t NMC, SpTMC, a1MC,a2MC, a3MC, a4MC, a1posMC, a2posMC, a3posMC, a4posMC, a1negMC, a2negMC, a3negMC, a4negMC;
        TClonesArray *stack = 0;
        TList *lst = fAOD->GetList();
        stack = (TClonesArray*)lst->FindObject(AliAODMCParticle::StdBranchName());
        int nMCTracks;
        if (!stack) nMCTracks = 0;
        else nMCTracks = stack->GetEntries();
        NMC  = 0;
        SpTMC = 0;
        for (Int_t i(0); i < nMCTracks; i++) {
            AliAODMCParticle *p1=(AliAODMCParticle*)stack->UncheckedAt(i);
            if (!p1) continue;
            if(p1->Charge()!=-3 && p1->Charge()!=+3) continue;// x3 by convention
            if(!p1->IsPrimary()) continue;
            if(!p1->IsPhysicalPrimary()) continue;
            if(fabs(p1->Eta()) > fEta ) continue;
            if(p1->Pt() < fPtmin|| p1->Pt() > fPtmax) continue;
            if((fabs(p1->GetPdgCode())==211)||(fabs(p1->GetPdgCode())==2212)||(fabs(p1->GetPdgCode())==321)){
                NMC++;
                SpTMC+=p1->Pt();
                etaMC.Fill(p1->Eta());
                if(p1->Charge() > 0) etaPosMC.Fill(p1->Eta());
                if(p1->Charge() < 0) etaNegMC.Fill(p1->Eta());
                fEtaMCCent->Fill(p1->Eta(),CentV0M);
            }
        }
        if(etaMC.Integral() < 1 || etaPosMC.Integral() < 1 || etaNegMC.Integral() < 1 ) return;
        etaMC.Scale(1.0/(etaMC.Integral()/20.));
        etaPosMC.Scale(1.0/(etaPosMC.Integral()/20.));
        etaNegMC.Scale(1.0/(etaNegMC.Integral()/20.));

        a1MC = GetAnCoeff(1,&etaMC);
        a2MC = GetAnCoeff(2,&etaMC);
        a3MC = GetAnCoeff(3,&etaMC);
        a4MC = GetAnCoeff(4,&etaMC);

        a1posMC = GetAnCoeff(1,&etaPosMC);
        a2posMC = GetAnCoeff(2,&etaPosMC);
        a3posMC = GetAnCoeff(3,&etaPosMC);
        a4posMC = GetAnCoeff(4,&etaPosMC);

        a1negMC = GetAnCoeff(1,&etaNegMC);
        a2negMC = GetAnCoeff(2,&etaNegMC);
        a3negMC = GetAnCoeff(3,&etaNegMC);
        a4negMC = GetAnCoeff(4,&etaNegMC);


        float fillnt[49] = {CentV0M, CentCL0, CentCL1, PVz, ZDCN1, ZDCP1, ZDCN2 ,ZDCP2, NTLs, NGlob, NTPC, NHyb, SpTGlob, SpTTPC, SpTHyb,a1Glob, a2Glob, a3Glob, a4Glob ,a1TPC, a2TPC, a3TPC, a4TPC, a1Hyb, a2Hyb, a3Hyb, a4Hyb, a1pos, a2pos, a3pos, a4pos, a1neg, a2neg ,a3neg, a4neg , NMC, SpTMC, a1MC, a2MC, a3MC, a4MC, a1posMC, a2posMC, a3posMC, a4posMC, a1negMC, a2negMC, a3negMC, a4negMC};
        fNt->Fill(fillnt);
    }
    else {
            float fillnt[35] = {CentV0M, CentCL0, CentCL1, PVz, ZDCN1, ZDCP1, ZDCN2 ,ZDCP2, NTLs, NGlob, NTPC, NHyb, SpTGlob, SpTTPC, SpTHyb,a1Glob, a2Glob, a3Glob, a4Glob ,a1TPC, a2TPC, a3TPC, a4TPC, a1Hyb, a2Hyb, a3Hyb, a4Hyb, a1pos, a2pos, a3pos, a4pos, a1neg, a2neg ,a3neg, a4neg};
            fNt->Fill(fillnt);
        }
   
    PostData(1, fOutputList);
                                                      
}

Float_t AliAnalysisTaskLongFluctuations::GetAnCoeff(Int_t order, TH1D *hist){
    Float_t an = 0;
    for(Int_t i=1;i<=hist->GetNbinsX();i++){
        an += sqrt((1.0/fEta)*(order+0.5))*(hist->GetBinContent(i)-1)*LegPol(order,hist->GetBinCenter(i));
    }
    an = an*hist->GetBinWidth(1);
    return an;
}

Float_t AliAnalysisTaskLongFluctuations::LegPol(Int_t order, Double_t x){
    if( order == 1 ) return x/fEta;
    if( order == 2 ) return 0.5*(3.0*pow(x/fEta,2)-1.0);
    if( order == 3 ) return 0.5*(5.0*pow(x/fEta,3)-3.0*x/fEta);
    if( order == 4 ) return (35.0*pow(x/fEta,4)-30.0*pow(x/fEta,2)+3.0)/8.0;
    else return 0.0;
}

//_____________________________________________________________________________
void AliAnalysisTaskLongFluctuations::Terminate(Option_t *)
{
    // terminate
    // called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
