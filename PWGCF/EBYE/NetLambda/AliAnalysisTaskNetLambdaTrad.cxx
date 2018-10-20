// For: Net Lambda fluctuation analysis via traditional method
// By: Ejiro Umaka Apr 2018
// Updated Oct 19; Provisions for feeddown correction
// Parts of the code taken from:
// AliEbyEPidEfficiencyContamination.cxx
// AliAnalysisTaskStrangenessVsMultiplicityMCRun2.cxx
// AliAnalysisTaskNetLambdaIdent.cxx

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliAODTrack.h"
#include "AliESDtrack.h"
#include "AliEventCuts.h"
#include "AliExternalTrackParam.h"
#include "AliAnalysisFilter.h"
#include "AliVMultiplicity.h"
#include "AliAnalysisUtils.h"
#include "AliAODMCParticle.h"
#include "AliStack.h"
#include "AliPIDResponse.h"
#include "AliMCEventHandler.h"
#include "AliV0vertexer.h"
#include "AliESDv0Cuts.h"
#include "AliMultSelection.h"
#include "TMath.h"
#include "TFile.h"
#include "TList.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TH3F.h"
#include "THnSparse.h"


#include "AliAnalysisTaskNetLambdaTrad.h"

ClassImp(AliAnalysisTaskNetLambdaTrad)

//-----------------------------------------------------------------------------
AliAnalysisTaskNetLambdaTrad::AliAnalysisTaskNetLambdaTrad(const char* name) :
AliAnalysisTaskSE(name),
fESD(0x0),
fAOD(0x0),
fPIDResponse(0x0),
fEventCuts(0),
fListHist(0x0),
fTreeV0(0x0),
fHistEventCounter(0x0),
fHistCentrality(0x0),
f2fHistGenCentVsPtLambda(0x0),
f2fHistGenCentVsPtAntiLambda(0x0),
f2fHistV0CentVsPtLambda(0x0),
f2fHistV0CentVsPtAntiLambda(0x0),
f2fHistInvMassVsPtLambda(0x0),
f2fHistInvMassVsPtLambdaRec(0x0),
f2fHistInvMassVsPtAntiLambda(0x0),
f2fHistInvMassVsPtAntiLambdaRec(0x0),
f2fHistRecPrimariesCentVsPtLambda(0x0),
f2fHistRecPrimariesCentVsPtAntiLambda(0x0),
f1fHistmassctLambda(0x0),
f1fHistmassctAntiLambda(0x0),
f2fHistPtmassctLambda(0x0),
f2fHistPtmassctAntiLambda(0x0),
f2fHistLambdaSecFromWeakDecay(0x0),
f2fHistAntiLambdaSecFromWeakDecay(0x0),
f2fHistLambdaMaterial(0x0),
f2fHistAntiLambdaMaterial(0x0),
f2fHistLambdaMisId(0x0),
f2fHistAntiLambdaMisId(0x0),
f2fHisthXiPlus(0x0),
f2fHisthXiMinus(0x0),
f2fHisthXiZero(0x0),
f2fHisthXiZeroAnti(0x0),
f2fHistLambdafromXi(0x0),
f2fHistAntiLambdafromXi(0x0),
f2fHistLRecstat(0x0),
f2fHistARecstat(0x0),
f2fHistLGenstat(0x0),
f2fHistAGenstat(0x0),
fCentrality(-1),
fTreeVariablePID(-1),
fTreeVariablePIDMother(-1),
fTreeVariablePIDPositive(-1),
fTreeVariablePIDNegative(-1),
fTreeVariableposTrMotherPID(-1),
fTreeVariablenegTrMotherPID(-1),
fNptBins(37),
fIsMC(kTRUE),
fIsAOD(kFALSE),
fEvSel(AliVEvent::kINT7),
fTreeVariableInvMassLambda(0),
fTreeVariableInvMassAntiLambda(0),
fTreeVariableDcaV0Daughters(0),
fTreeVariableDcaV0ToPrimVertex(0),
fTreeVariableDcaPosToPrimVertex(0),
fTreeVariableDcaNegToPrimVertex(0),

fPtBinNplusNminusV0(NULL),
fPtBinNplusNminusV0masscut(NULL),
fPtBinNplusNminusRec(NULL),
fPtBinNplusNminusRecPri(NULL),
fPtBinNplusNminusGen(NULL)

{
    Info("AliAnalysisTaskNetLambdaTrad","Calling Constructor");
    
    DefineInput(0,TChain::Class());
    DefineOutput(1,TList::Class());
    DefineOutput(2,TTree::Class());
    
}
//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

void AliAnalysisTaskNetLambdaTrad::UserCreateOutputObjects()
{
    fListHist = new TList();
    fListHist->SetOwner();
    
    fHistEventCounter = new TH1D( "fHistEventCounter", ";Evt. Sel. Step;Count",2,0,2);
    fHistEventCounter->GetXaxis()->SetBinLabel(1, "Processed");
    fHistEventCounter->GetXaxis()->SetBinLabel(2, "Selected");
    fListHist->Add(fHistEventCounter);
    
    fHistCentrality = new TH1D( "fHistCentrality", "fHistCentrality",100,0,100);
    fListHist->Add(fHistCentrality);
    
    const Int_t xNbins = 100;
    Double_t xBinEdge[xNbins+1];
    for(Int_t iBin = 0 ; iBin <= xNbins; iBin++)
    {
        xBinEdge[iBin] = iBin - 0.5;
    }
    
    const Double_t LambdaPtBins[]={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.5,8.0,10.0,12.0};
    const Int_t fNptBins=sizeof(LambdaPtBins)/sizeof(Double_t) - 1;
    
    Int_t    listbins=33;
    Double_t ptMax=12.;
    Double_t XiBins[listbins+1];
    for (Int_t i=0; i<=(listbins+1); i++) XiBins[i]=i*(ptMax+2)/listbins;
    
//    Printf("nBins is %d\n", fNptBins);
    
    f2fHistLRecstat = new TH2F("f2fHistLRecstat","f2fHistLRecstat", 100, -0.5, 99.5, 1900, -0.5, 1899.5);
    fListHist->Add(f2fHistLRecstat);
    
    f2fHistARecstat = new TH2F("f2fHistARecstat","f2fHistARecstat", 100, -0.5, 99.5, 1900, -0.5, 1899.5);
    fListHist->Add(f2fHistARecstat);
    
    f2fHistV0CentVsPtLambda = new TH2F("f2fHistV0CentVsPtLambda"," Centrality Vs #Lambda Rec Pt", xNbins, xBinEdge, fNptBins, LambdaPtBins);
    fListHist->Add(f2fHistV0CentVsPtLambda);
    
    f2fHistV0CentVsPtAntiLambda = new TH2F("f2fHistV0CentVsPtAntiLambda","Centrality Vs Rec #bar{#Lambda} Pt", xNbins, xBinEdge, fNptBins, LambdaPtBins);
    fListHist->Add(f2fHistV0CentVsPtAntiLambda);
    
    f2fHistInvMassVsPtLambda = new TH2F("f2fHistInvMassVsPtLambda","Inv mass #Lambda Vs Pt",100,1.08,1.16, fNptBins, LambdaPtBins);
    fListHist->Add(f2fHistInvMassVsPtLambda);
    
    f2fHistInvMassVsPtAntiLambda = new TH2F("f2fHistInvMassVsPtAntiLambda","Inv mass #bar{#Lambda} Vs Pt",100,1.08,1.16, fNptBins, LambdaPtBins);
    fListHist->Add(f2fHistInvMassVsPtAntiLambda);
    
    f2fHistPtmassctLambda = new TH2F("f2fHistPtmassctLambda","#Lambda masscut",xNbins, xBinEdge, fNptBins, LambdaPtBins);
    fListHist->Add(f2fHistPtmassctLambda);
    
    f2fHistPtmassctAntiLambda = new TH2F("f2fHistPtmassctAntiLambda","#bar{#Lambda} masscut", xNbins, xBinEdge, fNptBins, LambdaPtBins);
    fListHist->Add(f2fHistPtmassctAntiLambda);
    
    f1fHistmassctLambda = new TH1F("f1fHistmassctLambda","#Lambda masscut 1D",100,1.1,1.14);
    fListHist->Add(f1fHistmassctLambda);
    
    f1fHistmassctAntiLambda = new TH1F("f1fHistmassctAntiLambda","#bar{#Lambda} masscut 1D",100,1.1,1.14);
    fListHist->Add(f1fHistmassctAntiLambda);
    
    const Int_t dim = 75; //37 pt bins *2 + 1 cent bin
    Int_t bin[dim]    = { 100,500, 500, 500,500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 200, 200,200, 200, 200, 200,200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200,
        500, 500, 500,500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 500, 200, 200,200, 200, 200, 200,200, 200, 200, 200, 200, 200, 200, 200, 200, 200, 200};
    
    Double_t min[dim] =  { -0.5,-0.5, -0.5, -0.5,-0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,-0.5, -0.5, -0.5, -0.5,-0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
        -0.5, -0.5, -0.5,-0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,-0.5, -0.5, -0.5, -0.5,-0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5};

    Double_t max[dim] = { 99.5,499.5, 499.5, 499.5,499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 199.5, 199.5,199.5, 199.5, 199.5, 199.5,199.5, 199.5, 199.5, 199.5, 199.5, 199.5, 199.5, 199.5, 199.5, 199.5, 199.5,
        499.5, 499.5, 499.5,499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 199.5, 199.5,199.5, 199.5, 199.5, 199.5,199.5, 199.5, 199.5, 199.5, 199.5, 199.5, 199.5, 199.5, 199.5, 199.5, 199.5};
    
    
    fPtBinNplusNminusV0 = new THnSparseI("fPtBinNplusNminusV0","cent-nlambda-nantilambda", dim, bin, min, max);
    fListHist->Add(fPtBinNplusNminusV0);
    
    fPtBinNplusNminusV0masscut = new THnSparseI("fPtBinNplusNminusV0masscut","cent-nlambda-nantilambda masscut", dim, bin, min, max);
    fListHist->Add(fPtBinNplusNminusV0masscut);
    
    
    OpenFile(2);
    fTreeV0 = new TTree("fTreeV0","V0 Candidates");
    fTreeV0->Branch("fTreeVariableDcaV0Daughters",&fTreeVariableDcaV0Daughters,"fTreeVariableDcaV0Daughters/F");
    fTreeV0->Branch("fTreeVariableDcaV0ToPrimVertex",&fTreeVariableDcaV0ToPrimVertex,"fTreeVariableDcaV0ToPrimVertex/F");
    fTreeV0->Branch("fTreeVariableDcaPosToPrimVertex",&fTreeVariableDcaPosToPrimVertex,"fTreeVariableDcaPosToPrimVertex/F");
    fTreeV0->Branch("fTreeVariableDcaNegToPrimVertex",&fTreeVariableDcaNegToPrimVertex,"fTreeVariableDcaNegToPrimVertex/F");
    fTreeV0->Branch("fTreeVariableInvMassLambda",&fTreeVariableInvMassLambda,"fTreeVariableInvMassLambda/F");
    fTreeV0->Branch("fTreeVariableInvMassAntiLambda",&fTreeVariableInvMassAntiLambda,"fTreeVariableInvMassAntiLambda/F");
    
    
    if(fIsMC)
    {
        fPtBinNplusNminusGen = new THnSparseI("fPtBinNplusNminusGen","cent-nlambda-nantilambda", dim, bin, min, max);
        fListHist->Add(fPtBinNplusNminusGen);
        
        fPtBinNplusNminusRec = new THnSparseI("fPtBinNplusNminusRec","cent-nlambda-nantilambda", dim, bin, min, max);
        fListHist->Add(fPtBinNplusNminusRec);
        
        fPtBinNplusNminusRecPri = new THnSparseI("fPtBinNplusNminusRecPri","cent-nlambda-nantilambda", dim, bin, min, max);
        fListHist->Add(fPtBinNplusNminusRecPri);
        
        f2fHistGenCentVsPtLambda = new TH2F( "f2fHistGenCentVsPtLambda", "Centrality Vs #Lambda Gen Pt", xNbins, xBinEdge, fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistGenCentVsPtLambda);
        
        f2fHistGenCentVsPtAntiLambda = new TH2F( "f2fHistGenCentVsPtAntiLambda", "Centrality Vs #bar{#Lambda} Gen Pt", xNbins, xBinEdge, fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistGenCentVsPtAntiLambda);
        
        f2fHisthXiPlus = new TH2F( "f2fHisthXiPlus", "Centrality Vs Xi+ Gen Pt", xNbins, xBinEdge, fNptBins, LambdaPtBins);
        fListHist->Add(f2fHisthXiPlus);
        
        f2fHisthXiMinus = new TH2F( "f2fHisthXiMinus", "Centrality Vs Xi- Gen Pt", xNbins, xBinEdge, fNptBins, LambdaPtBins);
        fListHist->Add(f2fHisthXiMinus);
        
        f2fHisthXiZero = new TH2F( "f2fHisthXiZero", "Centrality Vs Xi Zero Gen Pt", xNbins, xBinEdge, fNptBins, LambdaPtBins);
        fListHist->Add(f2fHisthXiZero);
        
        f2fHisthXiZeroAnti = new TH2F( "f2fHisthXiZeroAnti", "Centrality Vs Xi Zero bar Gen Pt", xNbins, xBinEdge, fNptBins, LambdaPtBins);
        fListHist->Add(f2fHisthXiZeroAnti);
        
        f2fHistLambdafromXi = new TH2F( "f2fHistLambdafromXi", "Lambda pt vs Xi pt", fNptBins, LambdaPtBins, listbins, XiBins);
        fListHist->Add(f2fHistLambdafromXi);
        
        f2fHistAntiLambdafromXi = new TH2F( "f2fHistAntiLambdafromXi", "#bar{#Lambda} pt vs Xi pt", fNptBins, LambdaPtBins, listbins, XiBins);
        fListHist->Add(f2fHistAntiLambdafromXi);
        
        f2fHistInvMassVsPtLambdaRec = new TH2F("f2fHistInvMassVsPtLambdaRec","Inv mass #Lambda Vs Pt Rec",100,1.08,1.16, fNptBins,LambdaPtBins);
        fListHist->Add(f2fHistInvMassVsPtLambdaRec);
        
        f2fHistInvMassVsPtAntiLambdaRec = new TH2F("f2fHistInvMassVsPtAntiLambdaRec","Inv mass #bar{#Lambda} Vs Pt Rec",100,1.08,1.16,fNptBins,LambdaPtBins);
        fListHist->Add(f2fHistInvMassVsPtAntiLambdaRec);
        
        f2fHistLambdaSecFromWeakDecay = new TH2F("f2fHistLambdaSecFromWeakDecay","#Lambda from weak decays", xNbins, xBinEdge, fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistLambdaSecFromWeakDecay);
        
        f2fHistAntiLambdaSecFromWeakDecay = new TH2F("f2fHistAntiLambdaSecFromWeakDecay","#bar{#Lambda} from weak decays", xNbins, xBinEdge, fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistAntiLambdaSecFromWeakDecay);
        
        f2fHistLambdaMaterial = new TH2F("f2fHistLambdaMaterial","#Lambda from material", xNbins, xBinEdge, fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistLambdaMaterial);
        
        f2fHistAntiLambdaMaterial = new TH2F("f2fHistAntiLambdaMaterial","#bar{#Lambda} from material", xNbins, xBinEdge, fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistAntiLambdaMaterial);
        
        f2fHistLambdaMisId = new TH2F("f2fHistLambdaMisId","#Lambda Mis Identified", xNbins, xBinEdge, fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistLambdaMisId);
        
        f2fHistAntiLambdaMisId = new TH2F("f2fHistAntiLambdaMisId","#bar{#Lambda}  Mis Identified", xNbins, xBinEdge, fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistAntiLambdaMisId);
   
        f2fHistRecPrimariesCentVsPtLambda = new TH2F("f2fHistRecPrimariesCentVsPtLambda","#Lambda primaries", xNbins, xBinEdge, fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistRecPrimariesCentVsPtLambda);
        
        f2fHistRecPrimariesCentVsPtAntiLambda = new TH2F("f2fHistRecPrimariesCentVsPtAntiLambda","#bar{#Lambda} primaries", xNbins, xBinEdge, fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistRecPrimariesCentVsPtAntiLambda);
        
        f2fHistLGenstat = new TH2F("f2fHistLGenstat","f2fHistLGenstat", 100, -0.5, 99.5, 1900, -0.5, 1899.5);
        fListHist->Add(f2fHistLGenstat);
        
        f2fHistAGenstat = new TH2F("f2fHistAGenstat","f2fHistAGenstat", 100, -0.5, 99.5, 1900, -0.5, 1899.5);
        fListHist->Add(f2fHistAGenstat);
        
        fTreeV0->Branch("fTreeVariablePID",&fTreeVariablePID);
        fTreeV0->Branch("fTreeVariableposTrMotherPID",&fTreeVariableposTrMotherPID);
        fTreeV0->Branch("fTreeVariablenegTrMotherPID",&fTreeVariablenegTrMotherPID);
        fTreeV0->Branch("fTreeVariablePIDMother",&fTreeVariablePIDMother);
        fTreeV0->Branch("fTreeVariablePIDPositive",&fTreeVariablePIDPositive);
        fTreeV0->Branch("fTreeVariablePIDNegative",&fTreeVariablePIDNegative);
        
    }
    
    PostData(1,fListHist);
    PostData(2,fTreeV0);
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskNetLambdaTrad::UserExec(Option_t *)
{
    
    const Int_t dim = fNptBins*2;
    Int_t ptCh[dim];
    Int_t ptChCut[dim];
    Int_t ptChMC[dim];
    Int_t ptChMCrec[dim];
    Int_t ptChMCrecPri[dim];


    for(Int_t idx = 0; idx < dim; idx++)
    {
        ptCh[idx] = 0.;
        ptChCut[idx] = 0.;
        ptChMC[idx] = 0;
        ptChMCrec[idx] = 0;
        ptChMCrecPri[idx] = 0;
    }
    
    
    if (!fInputEvent) return;
    
    if(fIsAOD)
    {
        fAOD = dynamic_cast<AliAODEvent*>(fInputEvent);
        if (!fAOD) return;
    }
    else
    {
        fESD = dynamic_cast<AliESDEvent*>(InputEvent());
        if (!fESD) return;
    }
    
    fPIDResponse = fInputHandler->GetPIDResponse();
    if(!fPIDResponse) return;
    
    AliStack *stack = 0x0;
    if(fIsMC)
    {
        if(!fIsAOD)
        {
            AliMCEventHandler *mcH = dynamic_cast<AliMCEventHandler*>((AliAnalysisManager::GetAnalysisManager())->GetMCtruthEventHandler());
            if(!mcH) return;
            fMCEvent=mcH->MCEvent();
        }
        
        if(!fMCEvent) return;
        
        if(!fIsAOD)
        {
            stack = fMCEvent->Stack();
            if(!stack) return;
        }
    }
    
    AliMultSelection *MultSelection = (AliMultSelection*) fInputEvent->FindListObject("MultSelection");
    if(!MultSelection) return;
    
    if(!(fInputHandler->IsEventSelected() & fEvSel)) return;
    
    Double_t vVtx[3];
    
    AliVVertex *vvertex = (AliVVertex*)fInputEvent->GetPrimaryVertex();
    if (!vvertex) return;
    vVtx[0] = vvertex->GetX();
    vVtx[1] = vvertex->GetY();
    vVtx[2] = vvertex->GetZ();
    
    if(vVtx[2] < -10. || vVtx[2] > 10.) return;
    
    fCentrality = MultSelection->GetMultiplicityPercentile("V0M");
    
    if( fCentrality < 0 || fCentrality >=80 ) return;
    if (!fEventCuts.AcceptEvent(fInputEvent)) return;//pileup cut
    
    fHistEventCounter->Fill(1.5);
    fHistCentrality->Fill(fCentrality);
    
    
    Int_t nGen = 0;
    Double_t lRap = 0.0;
    Double_t nRecL = 0.0;
    Double_t nRecA = 0.0;
    Double_t nGenL = 0.0;
    Double_t nGenA = 0.0;
    
    if(fIsMC)
    {
        
        if(fIsAOD) nGen = fMCEvent->GetNumberOfTracks();
        else nGen = stack->GetNtrack();
        for(Int_t iGen = 0; iGen < nGen; iGen++)
        {
            Int_t genpid = -1;
            Float_t gpt = 0.0, eta = 0.0,abseta =0.0;
            
            if(fIsAOD)
            {
                AliAODMCParticle* mctrack = (AliAODMCParticle*)fMCEvent->GetTrack(iGen);
                if(!mctrack) continue;
                if(!(mctrack->IsPhysicalPrimary())) continue;
                genpid = mctrack->PdgCode();
                gpt = mctrack->Pt();
                eta = mctrack->Eta();
            }
            else
            {
                TParticle* mctrack = stack->Particle(iGen);
                if(!mctrack) continue;
                if(!(stack->IsPhysicalPrimary(iGen))) continue;
                genpid = mctrack->GetPdgCode();
                gpt = mctrack->Pt();
                eta = mctrack->Eta();
            }
            
            abseta = TMath::Abs(eta);
            if(abseta > 0.8) continue;
            
            Int_t iptbinMC = GetPtBin(gpt);
            if( iptbinMC < 0 || iptbinMC > fNptBins-1 ) continue;
            
            if(genpid == 3122)
            {
                f2fHistGenCentVsPtLambda->Fill(fCentrality, gpt);
                nGenL += 1.;
                ptChMC[iptbinMC] += 1;
            }
            
            if(genpid == -3122)
            {
                f2fHistGenCentVsPtAntiLambda->Fill(fCentrality,gpt);
                nGenA += 1.;
                ptChMC[iptbinMC+fNptBins] += 1;
            }
            
            if(genpid == -3312)
            {
                f2fHisthXiPlus->Fill(fCentrality,gpt);
            }
            if(genpid == 3312)
            {
                f2fHisthXiMinus->Fill(fCentrality,gpt);
            }
            if(genpid == 3322)
            {
                f2fHisthXiZero->Fill(fCentrality,gpt);
            }
            if(genpid == -3322)
            {
                f2fHisthXiZeroAnti->Fill(fCentrality,gpt);
            }
            
        } // end loop over generated particles
        
        f2fHistLGenstat->Fill(fCentrality, nGenL);
        f2fHistAGenstat->Fill(fCentrality, nGenA);
        
        
        Double_t ptContainerMC[dim+1];
        ptContainerMC[0] = (Double_t) fCentrality;
        for(Int_t i = 1; i <= dim; i++)
        {
            ptContainerMC[i] = ptChMC[i-1];
        }
        fPtBinNplusNminusGen->Fill(ptContainerMC);
        
    }
    
    Int_t nV0 = 0;
    if(fIsAOD) nV0 = fAOD->GetNumberOfV0s();
    else nV0 = fESD->GetNumberOfV0s();
    AliESDv0 *esdv0 = 0x0;
    AliESDtrack *esdpTrack = 0x0;
    AliESDtrack *esdnTrack = 0x0;
    AliAODv0 *aodv0 = 0x0;
    AliAODTrack *aodpTrack = 0x0;
    AliAODTrack *aodnTrack = 0x0;
    
    for(Int_t iV0 = 0; iV0 < nV0; iV0++)
    {
        esdv0 = 0x0;
        esdpTrack = 0x0;
        esdnTrack = 0x0;
        aodv0 = 0x0;
        aodpTrack = 0x0;
        aodnTrack = 0x0;
        
        Double_t  vertx[3];
        
        Float_t invMassLambda = -999, invMassAntiLambda = -999;
        Float_t V0pt = -999, eta = -999, pmom = -999;
        Float_t ppt = -999,  peta = -999, posprnsg = -999;
        Float_t npt = -999,  neta = -999, negprnsg = -999;
        Bool_t  ontheflystat = kFALSE;
        Float_t dcaPosToVertex = -999, dcaNegToVertex = -999, dcaDaughters = -999, dcaV0ToVertex = -999, cosPointingAngle = -999;
        
        if(fIsAOD)
        {
            aodv0 = fAOD->GetV0(iV0);
            if(!aodv0) continue;
            aodpTrack = (AliAODTrack*)aodv0->GetDaughter(0);
            if(!aodpTrack) continue;
            aodnTrack = (AliAODTrack*)aodv0->GetDaughter(1);
            if(!aodnTrack) continue;
            
            aodv0->GetXYZ(vertx);
            if(aodpTrack->Charge() == aodnTrack->Charge())
            {
                continue;
            }
            
            Float_t lPosTrackCrossedRows = aodpTrack->GetTPCClusterInfo(2,1);
            Float_t lNegTrackCrossedRows = aodnTrack->GetTPCClusterInfo(2,1);
            fTreeVariableLeastNbrCrossedRows = lPosTrackCrossedRows;
            if( lNegTrackCrossedRows < fTreeVariableLeastNbrCrossedRows )
                fTreeVariableLeastNbrCrossedRows = (Int_t) lNegTrackCrossedRows;
            
            if( !(aodpTrack->GetStatus() & AliAODTrack::kTPCrefit)) continue;
            if( !(aodnTrack->GetStatus() & AliAODTrack::kTPCrefit)) continue;
            
            if ( ( ( aodpTrack->GetTPCClusterInfo(2,1) ) < 70 ) || ( ( aodnTrack->GetTPCClusterInfo(2,1) ) < 70 ) ) continue;
            
            if( aodpTrack->GetKinkIndex(0)>0 || aodnTrack->GetKinkIndex(0)>0 ) continue;
            
            if( aodpTrack->GetTPCNclsF()<=0 || aodnTrack->GetTPCNclsF()<=0 ) continue;
            
            Float_t lPosTrackCrossedRowsOverFindable = lPosTrackCrossedRows / ((double)(aodpTrack->GetTPCNclsF()));
            Float_t lNegTrackCrossedRowsOverFindable = lNegTrackCrossedRows / ((double)(aodnTrack->GetTPCNclsF()));
            
            fTreeVariableLeastRatioCrossedRowsOverFindable = lPosTrackCrossedRowsOverFindable;
            if( lNegTrackCrossedRowsOverFindable < fTreeVariableLeastRatioCrossedRowsOverFindable)
                fTreeVariableLeastRatioCrossedRowsOverFindable = lNegTrackCrossedRowsOverFindable;
            
            if (fTreeVariableLeastRatioCrossedRowsOverFindable < 0.8 ) continue;
            
            posprnsg = fPIDResponse->NumberOfSigmasTPC(aodpTrack, AliPID::kProton);
            negprnsg = fPIDResponse->NumberOfSigmasTPC(aodnTrack, AliPID::kProton);
            
            V0pt = aodv0->Pt();
            pmom = aodv0->P();
            eta = aodv0->Eta();
            ppt = aodpTrack->Pt();
            peta = aodpTrack->Eta();
            npt = aodnTrack->Pt();
            neta = aodnTrack->Eta();
            ontheflystat = aodv0->GetOnFlyStatus();
            
            
            dcaPosToVertex = aodv0->DcaPosToPrimVertex();
            dcaNegToVertex = aodv0->DcaNegToPrimVertex();
            cosPointingAngle = aodv0->CosPointingAngle(vVtx);
            dcaDaughters = aodv0->DcaV0Daughters();
            dcaV0ToVertex = aodv0->DcaV0ToPrimVertex();
            invMassLambda = aodv0->MassLambda();
            invMassAntiLambda = aodv0->MassAntiLambda();
            
        }
        else
        {
            esdv0 = fESD->GetV0(iV0);
            if(!esdv0) continue;
            esdpTrack =  (AliESDtrack*)fESD->GetTrack(TMath::Abs(esdv0->GetPindex()));
            if(!esdpTrack) continue;
            esdnTrack = (AliESDtrack*)fESD->GetTrack(TMath::Abs(esdv0->GetNindex()));
            if(!esdnTrack) continue;
            esdv0->GetXYZ(vertx[0], vertx[1], vertx[2]);
            
            if(esdpTrack->Charge() == esdnTrack->Charge())
            {
                continue;
            }
            Float_t lPosTrackCrossedRows = esdpTrack->GetTPCClusterInfo(2,1);
            Float_t lNegTrackCrossedRows = esdnTrack->GetTPCClusterInfo(2,1);
            
            fTreeVariableLeastNbrCrossedRows = (Int_t) lPosTrackCrossedRows;
            if( lNegTrackCrossedRows < fTreeVariableLeastNbrCrossedRows )
                fTreeVariableLeastNbrCrossedRows = (Int_t) lNegTrackCrossedRows;
            
            if( !(esdpTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
            if( !(esdnTrack->GetStatus() & AliESDtrack::kTPCrefit)) continue;
            
            if ( ( ( esdpTrack->GetTPCClusterInfo(2,1) ) < 70 ) || ( ( esdnTrack->GetTPCClusterInfo(2,1) ) < 70 ) ) continue;
            
            if( esdpTrack->GetKinkIndex(0)>0 || esdnTrack->GetKinkIndex(0)>0 ) continue;
            
            if( esdpTrack->GetTPCNclsF()<=0 || esdnTrack->GetTPCNclsF()<=0 ) continue;
            
            Float_t lPosTrackCrossedRowsOverFindable = lPosTrackCrossedRows / ((double)(esdpTrack->GetTPCNclsF()));
            Float_t lNegTrackCrossedRowsOverFindable = lNegTrackCrossedRows / ((double)(esdnTrack->GetTPCNclsF()));
            
            fTreeVariableLeastRatioCrossedRowsOverFindable = lPosTrackCrossedRowsOverFindable;
            if( lNegTrackCrossedRowsOverFindable < fTreeVariableLeastRatioCrossedRowsOverFindable )
                fTreeVariableLeastRatioCrossedRowsOverFindable = lNegTrackCrossedRowsOverFindable;
            
            if ( fTreeVariableLeastRatioCrossedRowsOverFindable < 0.8 ) continue;
            
            V0pt = esdv0->Pt();
            pmom = esdv0->P();
            eta = esdv0->Eta();
            ppt = esdpTrack->Pt();
            peta = esdpTrack->Eta();
            npt = esdnTrack->Pt();
            neta = esdnTrack->Eta();
            ontheflystat = esdv0->GetOnFlyStatus();
            
            posprnsg = fPIDResponse->NumberOfSigmasTPC(esdpTrack, AliPID::kProton);
            negprnsg = fPIDResponse->NumberOfSigmasTPC(esdnTrack, AliPID::kProton);
            
            Float_t vb1, vb2;
            Float_t e1, e2;
            esdpTrack->GetImpactParameters(vb1,vb2);
            dcaPosToVertex = TMath::Sqrt(vb1*vb1+vb2*vb2);
            esdnTrack->GetImpactParameters(e1,e2);
            dcaNegToVertex = TMath::Sqrt(e1*e1+e2*e2);
            cosPointingAngle = esdv0->GetV0CosineOfPointingAngle();
            dcaDaughters = esdv0->GetDcaV0Daughters();
            dcaV0ToVertex = esdv0->GetD(vVtx[0],vVtx[1],vVtx[2]);
            
            esdv0->ChangeMassHypothesis(3122);
            invMassLambda = esdv0->GetEffMass();
            esdv0->ChangeMassHypothesis(-3122);
            invMassAntiLambda = esdv0->GetEffMass();
        }
        
        Float_t v0Radius = TMath::Sqrt(vertx[0]*vertx[0]+vertx[1]*vertx[1]);
        if(TMath::Abs(peta) > 1) continue;
        if(TMath::Abs(neta) > 1) continue;
        if(cosPointingAngle < 0.999) continue;
        if(dcaDaughters > 1.5) continue;
        if(v0Radius < 5.0) continue;
        if(v0Radius > 200.) continue;
        
        fTreeVariableDcaV0ToPrimVertex = dcaV0ToVertex;
        fTreeVariableDcaV0Daughters = dcaDaughters;
        fTreeVariableDcaPosToPrimVertex = dcaPosToVertex;
        fTreeVariableDcaNegToPrimVertex = dcaNegToVertex;
        fTreeVariableInvMassLambda = invMassLambda;
        fTreeVariableInvMassAntiLambda = invMassAntiLambda;
        
        
        Int_t iptbin = GetPtBin(V0pt);
        if( iptbin < 0 || iptbin > fNptBins-1 ) continue;
        
        if( ontheflystat == 0 )
        {
            if(dcaV0ToVertex < 0.5 && dcaNegToVertex > 0.1 && dcaPosToVertex >  0.05 && TMath::Abs(posprnsg)  <= 3.)
            {
                f2fHistInvMassVsPtLambda->Fill(invMassLambda,V0pt);
                f2fHistV0CentVsPtLambda->Fill(fCentrality,V0pt);
                ptCh[iptbin] += 1;
                nRecL += 1.;

                if(invMassLambda > 1.11 && invMassLambda < 1.12)
                {
                    f2fHistPtmassctLambda->Fill(fCentrality,V0pt);
                    f1fHistmassctLambda->Fill(invMassLambda);
                    ptChCut[iptbin] += 1;
                }
            }
            if(dcaV0ToVertex < 0.5 && dcaNegToVertex > 0.05 && dcaPosToVertex >  0.1 && TMath::Abs(negprnsg)  <= 3.)
            {
                f2fHistInvMassVsPtAntiLambda->Fill(invMassAntiLambda,V0pt);
                f2fHistV0CentVsPtAntiLambda->Fill(fCentrality,V0pt);
                ptCh[iptbin+fNptBins] += 1;
                nRecA += 1.;

                if(invMassAntiLambda > 1.11 && invMassAntiLambda < 1.12)
                {
                    f2fHistPtmassctAntiLambda->Fill(fCentrality,V0pt);
                    f1fHistmassctAntiLambda->Fill(invMassAntiLambda);
                    ptChCut[iptbin+fNptBins] += 1;
                }
            }
            
            
            if(fIsMC)
            {
                fTreeVariablePID = -999, fTreeVariablePIDPositiveGrandMother = -999, fTreeVariablePIDNegativeGrandMother = -999,fTreeVariablePrimaryStatusMother = -999, fTreeVariablePIDMother = -999;
                Float_t mcpt = -999, mceta = -999, fTreeVariablePtMother= -999;
                Bool_t isPrim = kFALSE, isSecFromMaterial = kFALSE, isSecFromWeakDecay = kFALSE;
                fTreeVariablePIDPositive = -999;
                fTreeVariablePIDNegative = -999;
                
                if(fIsAOD)
                {
                    if(TMath::Abs(aodpTrack->GetLabel()) >= nGen || TMath::Abs(aodnTrack->GetLabel()) >= nGen) continue;
                    AliAODMCParticle *aodGenTrackPos = (AliAODMCParticle*)fMCEvent->GetTrack(TMath::Abs(aodpTrack->GetLabel()));
                    if(!aodGenTrackPos) continue;
                    AliAODMCParticle *aodGenTrackNeg = (AliAODMCParticle*)fMCEvent->GetTrack(TMath::Abs(aodnTrack->GetLabel()));
                    if(!aodGenTrackNeg) continue;
                    
                    Int_t lPIDPositive = aodGenTrackPos -> PdgCode();
                    Int_t lPIDNegative = aodGenTrackNeg -> PdgCode();
                    
                    fTreeVariablePIDPositive = lPIDPositive;
                    fTreeVariablePIDNegative = lPIDNegative;
                    
                    Int_t posTparticle = aodGenTrackPos->GetMother();
                    Int_t negTparticle = aodGenTrackNeg->GetMother();
                    
                    if( posTparticle == negTparticle && posTparticle > 0 )
                    {
                        AliVParticle *lthisV0 = fMCEvent->GetTrack(posTparticle);
                        if(!lthisV0) continue;
                        fTreeVariablePID = lthisV0->PdgCode();
                        
                        mcpt = lthisV0->Pt();
                        mceta = lthisV0->Eta();
                        
                        isSecFromMaterial = (static_cast<AliAODMCParticle*>(lthisV0))->IsSecondaryFromMaterial();
                        isSecFromWeakDecay = (static_cast<AliAODMCParticle*>(lthisV0))->IsSecondaryFromWeakDecay();
                        isPrim = (static_cast<AliAODMCParticle*>(lthisV0))->IsPhysicalPrimary();
                    }
                    
                    if(dcaV0ToVertex < 0.5 && dcaNegToVertex > 0.1 && dcaPosToVertex >  0.05 && TMath::Abs(posprnsg)  <= 3.)
                    {
                        if(fTreeVariablePID == 3122)
                        {
                            if(isPrim)
                            {
                                f2fHistRecPrimariesCentVsPtLambda->Fill(fCentrality,mcpt);
                                ptChMCrecPri[iptbin] += 1;
                            }
                            else{f2fHistLambdaMisId->Fill(fCentrality,mcpt);}

                            if(isSecFromWeakDecay) {f2fHistLambdaSecFromWeakDecay->Fill(fCentrality,mcpt);}
                            else if (isSecFromMaterial) {f2fHistLambdaMaterial->Fill(fCentrality,mcpt);}
                            f2fHistInvMassVsPtLambdaRec->Fill(invMassLambda,mcpt);
                            ptChMCrec[iptbin] += 1;
                        }
                    }
                    
                    if(dcaV0ToVertex < 0.5 && dcaNegToVertex > 0.05 && dcaPosToVertex >  0.1 && TMath::Abs(negprnsg)  <= 3.)
                    {
                        if(fTreeVariablePID == -3122)
                        {
                            if(isPrim)
                            {
                                f2fHistRecPrimariesCentVsPtAntiLambda->Fill(fCentrality,mcpt);
                                ptChMCrecPri[iptbin+fNptBins] += 1;
                            }
                            else{f2fHistAntiLambdaMisId->Fill(fCentrality,mcpt);}
                            if(isSecFromWeakDecay) {f2fHistAntiLambdaSecFromWeakDecay->Fill(fCentrality,mcpt);}
                            else if (isSecFromMaterial) {f2fHistAntiLambdaMaterial->Fill(fCentrality,mcpt);}
                            f2fHistInvMassVsPtAntiLambdaRec->Fill(invMassAntiLambda,mcpt);
                            ptChMCrec[iptbin+fNptBins] += 1;
                        }
                    }
                }
                else
                {
                    
                    if(TMath::Abs(esdpTrack->GetLabel()) >= nGen || TMath::Abs(esdnTrack->GetLabel()) >= nGen) continue;
                    Int_t lblPosV0Dghter = (Int_t) TMath::Abs( esdpTrack->GetLabel() );
                    Int_t lblNegV0Dghter = (Int_t) TMath::Abs( esdnTrack->GetLabel() );
                    
                    TParticle* esdGenTrackPos = stack->Particle( lblPosV0Dghter );
                    if(!esdGenTrackPos) continue;
                    TParticle* esdGenTrackNeg = stack->Particle( lblNegV0Dghter );
                    if(!esdGenTrackNeg) continue;
                    
                    Int_t lPIDPositive = esdGenTrackPos -> GetPdgCode();
                    Int_t lPIDNegative = esdGenTrackNeg -> GetPdgCode();
                    
                    fTreeVariablePIDPositive = lPIDPositive;
                    fTreeVariablePIDNegative = lPIDNegative;
                    
                    Int_t posTparticle = esdGenTrackPos->GetFirstMother();
                    if(posTparticle < 0) continue;
                    Int_t negTparticle = esdGenTrackNeg->GetFirstMother();
                    if(negTparticle < 0) continue;
                    
                    if ( posTparticle > 0 )
                    {
                        TParticle *lPosMother = stack->Particle(posTparticle);
                        fTreeVariableposTrMotherPID = lPosMother->GetPdgCode();
                        Int_t lblGrandMother = lPosMother->GetFirstMother();
                        if( lblGrandMother > 0 ){
                            TParticle *lPosGrandMother = stack->Particle( lblGrandMother );
                            fTreeVariablePIDPositiveGrandMother = lPosGrandMother->GetPdgCode();
                        }
                    }
                    
                    if ( negTparticle > 0 )
                    {
                        TParticle *lNegMother = stack->Particle(negTparticle);
                        fTreeVariablenegTrMotherPID = lNegMother->GetPdgCode();
                        Int_t lblGrandMother = lNegMother->GetFirstMother();
                        if( lblGrandMother > 0 ){
                            TParticle *lNegGrandMother = stack->Particle( lblGrandMother );
                            fTreeVariablePIDNegativeGrandMother = lNegGrandMother->GetPdgCode();
                        }
                    }
                    
                    if( posTparticle == negTparticle && posTparticle > 0 )
                    {
                        TParticle *esdlthisV0 = stack->Particle(posTparticle);
                        if(!esdlthisV0) continue;
                        
                        fTreeVariablePID = esdlthisV0->GetPdgCode();
                        mcpt = esdlthisV0->Pt();
                        mceta = esdlthisV0->Eta();
                        isSecFromMaterial = stack->IsSecondaryFromMaterial(posTparticle);
                        isSecFromWeakDecay = stack->IsSecondaryFromWeakDecay(posTparticle);
                        isPrim = stack->IsPhysicalPrimary(posTparticle);
                        
                        Int_t lblThisV0Parent = esdlthisV0->GetFirstMother();
                        if ( lblThisV0Parent > 0 )
                        {
                            TParticle *pThisV0Parent = stack->Particle(lblThisV0Parent);
                            if(!pThisV0Parent) continue;
                            fTreeVariablePIDMother   = pThisV0Parent->GetPdgCode();
                            fTreeVariablePtMother    = pThisV0Parent->Pt();
                            if(stack->IsPhysicalPrimary (lblThisV0Parent) ) fTreeVariablePrimaryStatusMother = 1;
                        }
     
                    }
                    
                    if(dcaV0ToVertex < 0.5 && dcaNegToVertex > 0.1 && dcaPosToVertex >  0.05 && TMath::Abs(posprnsg)  <= 3.)
                    {
                        if(fTreeVariablePID == 3122)
                            {
                                if(isPrim)
                                {
                                    f2fHistRecPrimariesCentVsPtLambda->Fill(fCentrality,mcpt);
                                     ptChMCrecPri[iptbin] += 1;
                                }
                                else{f2fHistLambdaMisId->Fill(fCentrality,mcpt);}
                                if(isSecFromWeakDecay) {f2fHistLambdaSecFromWeakDecay->Fill(fCentrality,mcpt);}
                                else if (isSecFromMaterial) {f2fHistLambdaMaterial->Fill(fCentrality,mcpt);}
                                f2fHistInvMassVsPtLambdaRec->Fill(invMassLambda,mcpt);
                                ptChMCrec[iptbin] += 1;
                            }
                     
                    }
                    if(dcaV0ToVertex < 0.5 && dcaNegToVertex > 0.05 && dcaPosToVertex >  0.1 && TMath::Abs(negprnsg)  <= 3.)
                    {
                        if(fTreeVariablePID == -3122)
                        {
                            if(isPrim)
                            {
                                f2fHistRecPrimariesCentVsPtAntiLambda->Fill(fCentrality,mcpt);
                                ptChMCrecPri[iptbin+fNptBins] += 1;
                            }
                            else{f2fHistAntiLambdaMisId->Fill(fCentrality,mcpt);}
                            if(isSecFromWeakDecay) {f2fHistAntiLambdaSecFromWeakDecay->Fill(fCentrality,mcpt);}
                            else if (isSecFromMaterial) {f2fHistAntiLambdaMaterial->Fill(fCentrality,mcpt);}
                            f2fHistInvMassVsPtAntiLambdaRec->Fill(invMassAntiLambda,mcpt);
                            ptChMCrec[iptbin+fNptBins] += 1;
                        }
                    }
                    
                    if (fTreeVariablePID == 3122 && isSecFromWeakDecay && fTreeVariablePrimaryStatusMother == 1)
                    {
                        if ((fTreeVariablePIDMother == 3312)||(fTreeVariablePIDMother == 3322))
                            {f2fHistLambdafromXi->Fill(mcpt,fTreeVariablePtMother);}
                    }
                    
                    if (fTreeVariablePID == -3122 && isSecFromWeakDecay && fTreeVariablePrimaryStatusMother == 1)
                    {
                        if ((fTreeVariablePIDMother == -3312)||(fTreeVariablePIDMother == -3322))
                        {f2fHistAntiLambdafromXi->Fill(mcpt,fTreeVariablePtMother);}
                    }
                
               } //ESD
            } //MC condition
        }// zero onfly V0
    }// end of V0 loop
    
    
    f2fHistLRecstat->Fill(fCentrality, nRecL);
    f2fHistARecstat->Fill(fCentrality, nRecA);
    
    Double_t ptContainer[dim+1];
    ptContainer[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainer[i] = ptCh[i-1];
    }
    fPtBinNplusNminusV0->Fill(ptContainer);


    Double_t ptContainerCut[dim+1];
    ptContainerCut[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerCut[i] = ptChCut[i-1];
    }
    fPtBinNplusNminusV0masscut->Fill(ptContainerCut);
    
    
    Double_t ptContainerRec[dim+1];
    ptContainerRec[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerRec[i] = ptChMCrec[i-1];
    }
    fPtBinNplusNminusRec->Fill(ptContainerRec);
    
    Double_t ptContainerRecPri[dim+1];
    ptContainerRecPri[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerRecPri[i] = ptChMCrecPri[i-1];
    }
    fPtBinNplusNminusRecPri->Fill(ptContainerRecPri);

    fTreeV0->Fill();
    PostData(1,fListHist);
    PostData(2,fTreeV0);
}


//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Int_t AliAnalysisTaskNetLambdaTrad::GetPtBin(Double_t pt)
{
    Int_t bin = -1;
    const Double_t LambdaPtBins[]={0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.5,8.0,10.0,12.0};
    for(Int_t iBin = 0; iBin < fNptBins; iBin++)
    {
        
        if( iBin == fNptBins-1){
            if( pt >= LambdaPtBins[iBin] && pt <= LambdaPtBins[iBin+1]){
                bin = iBin;
                break;
            }
        }
        else{
            if( pt >= LambdaPtBins[iBin] && pt < LambdaPtBins[iBin+1]){
                bin = iBin;
                break;
                
            }
        }
    }
    
    return bin;
    
}






