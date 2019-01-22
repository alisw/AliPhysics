
// For: Net Lambda fluctuation analysis via traditional method
// By: Ejiro Naomi Umaka Apr 2018
// Updated Jan 21 


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
f2fHistRecCentVsPtLambda(0x0),
f2fHistRecCentVsPtAntiLambda(0x0),
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
f2fHistLambdaRecPt(0x0),
f2fHistAntiLambdaRecPt(0x0),
f2fHistLRecstat(0x0),
f2fHistARecstat(0x0),
f2fHistLGenstat(0x0),
f2fHistAGenstat(0x0),
f2fHistXiPlus(0x0),
f2fHistXiMinus(0x0),
f2fHistLambdafromXi(0x0),
f2fHistAntiLambdafromXi(0x0),
fCentrality(-1),
fTreeVariablePID(-1),
fTreeVariablePIDPositive(-1),
fTreeVariablePIDNegative(-1),
fNptBins(23),
fIsMC(kTRUE),
fIsAOD(kFALSE),
fEvSel(AliVEvent::kINT7),
fTreeVariableInvMassLambda(0),
fTreeVariableInvMassAntiLambda(0),
fTreeVariableDcaV0Daughters(0),
fTreeVariableDcaV0ToPrimVertex(0),
fTreeVariableDcaPosToPrimVertex(0),
fTreeVariableDcaNegToPrimVertex(0),


fPtBinNplusNminusCh(NULL),
fPtBinNplusNminusChCut(NULL),
fPtBinNplusNminusChTruth(NULL)

{
    Info("AliAnalysisTaskNetLambdaTrad","Calling Constructor");
    
    DefineInput(0,TChain::Class());
    DefineOutput(1,TList::Class());
    DefineOutput(2,TTree::Class());
    
}
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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
    for(Int_t iBin = 0 ; iBin <= xNbins; iBin++){
        xBinEdge[iBin] = iBin - 0.5;
    }
    
    Double_t LambdaPtBins[24] =  {0.5,0.7,0.9,1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0,4.2};
    Double_t xibinlimits[24] = {0.00,  0.20,  0.40,  0.60,  0.80,  0.90,1.00,  1.10,  1.20,  1.30,  1.40,  1.50, 1.70,  1.90,  2.20,  2.60,  3.10,  3.90,4.90,  6.00,  7.20,  8.50 ,10.00, 12.00};
    Long_t xibinnumb = sizeof(xibinlimits)/sizeof(Double_t) - 1;
    
    //V0 hists//
    
    f2fHistRecCentVsPtLambda = new TH2F("f2fHistRecCentVsPtLambda"," Centrality Vs #Lambda Rec Pt", xNbins, xBinEdge, fNptBins, LambdaPtBins);
    fListHist->Add(f2fHistRecCentVsPtLambda);
    
    f2fHistRecCentVsPtAntiLambda = new TH2F("f2fHistRecCentVsPtAntiLambda","Centrality Vs Rec #bar{#Lambda} Pt", xNbins, xBinEdge, fNptBins, LambdaPtBins);
    fListHist->Add(f2fHistRecCentVsPtAntiLambda);
    
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
    
    const Int_t dim = 47; //23 pt bins + 1 cent bin
    Int_t bin[dim]    = { 100,
        500, 500, 500,
        500, 500, 500, 500, 500, 500, 500, 500,500,500,500,500,
        200, 200, 200, 200, 200, 200, 200, 200,
        500, 500, 500,
        500, 500, 500, 500, 500, 500, 500, 500,500,500,500,500,
        200, 200, 200, 200, 200, 200, 200, 200 };
    
    Double_t min[dim] = { -0.5,
        -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,-0.5,-0.5,-0.5,-0.5,
        -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
        -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,-0.5,-0.5,-0.5,-0.5,
        -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5, -0.5,
        -0.5, -0.5};
    
    Double_t max[dim] = { 99.5,
        499.5, 499.5, 499.5,
        499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5,499.5,499.5,499.5,499.5,
        199.5, 199.5, 199.5, 199.5, 199.5, 199.5, 199.5, 199.5,
        499.5, 499.5, 499.5,
        499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5, 499.5,499.5,499.5,499.5,499.5,
        199.5, 199.5, 199.5, 199.5, 199.5, 199.5, 199.5, 199.5 };
    
    fPtBinNplusNminusChCut = new THnSparseI("fPtBinNplusNminusChCut","cent-nlambda-nantilambda masscut", dim, bin, min, max);
    fListHist->Add(fPtBinNplusNminusChCut); //V0masscut
    
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
        fPtBinNplusNminusChTruth = new THnSparseI("fPtBinNplusNminusChTruth","cent-nlambda-nantilambda", dim, bin, min, max);
        fListHist->Add(fPtBinNplusNminusChTruth); //gen
        
        fPtBinNplusNminusCh = new THnSparseI("fPtBinNplusNminusCh","cent-nlambda-nantilambda", dim, bin, min, max);
        fListHist->Add(fPtBinNplusNminusCh); //recpri
        
        f2fHistLRecstat = new TH2F("f2fHistLRecstat","f2fHistLRecstat", 100, -0.5, 99.5, 1900, -0.5, 1899.5);
        fListHist->Add(f2fHistLRecstat);
        
        f2fHistARecstat = new TH2F("f2fHistARecstat","f2fHistARecstat", 100, -0.5, 99.5, 1900, -0.5, 1899.5);
        fListHist->Add(f2fHistARecstat);
        
        f2fHistInvMassVsPtLambdaRec = new TH2F("f2fHistInvMassVsPtLambdaRec","Inv mass #Lambda Vs Pt Rec",100,1.08,1.16, fNptBins,LambdaPtBins);
        fListHist->Add(f2fHistInvMassVsPtLambdaRec);
        
        f2fHistInvMassVsPtAntiLambdaRec = new TH2F("f2fHistInvMassVsPtAntiLambdaRec","Inv mass #bar{#Lambda} Vs Pt Rec",100,1.08,1.16,fNptBins,LambdaPtBins);
        fListHist->Add(f2fHistInvMassVsPtAntiLambdaRec);
        
        f2fHistGenCentVsPtLambda = new TH2F( "f2fHistGenCentVsPtLambda", "Centrality Vs #Lambda Gen Pt", xNbins, xBinEdge, fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistGenCentVsPtLambda);
        
        f2fHistGenCentVsPtAntiLambda = new TH2F( "f2fHistGenCentVsPtAntiLambda", "Centrality Vs #bar{#Lambda} Gen Pt", xNbins, xBinEdge, fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistGenCentVsPtAntiLambda);
        
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
        
        f2fHistLambdaRecPt = new TH2F("f2fHistLambdaRecPt","#Lambda rec pt", xNbins, xBinEdge, fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistLambdaRecPt);
        
        f2fHistAntiLambdaRecPt = new TH2F("f2fHistAntiLambdaRecPt","#bar{#Lambda} rec pt", xNbins, xBinEdge, fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistAntiLambdaRecPt);
        
        f2fHistRecPrimariesCentVsPtLambda = new TH2F("f2fHistRecPrimariesCentVsPtLambda","#Lambda primaries", xNbins, xBinEdge, fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistRecPrimariesCentVsPtLambda);
        
        f2fHistRecPrimariesCentVsPtAntiLambda = new TH2F("f2fHistRecPrimariesCentVsPtAntiLambda","#bar{#Lambda} primaries", xNbins, xBinEdge, fNptBins, LambdaPtBins);
        fListHist->Add(f2fHistRecPrimariesCentVsPtAntiLambda);
        
        f2fHistLGenstat = new TH2F("f2fHistLGenstat","f2fHistLGenstat", 100, -0.5, 99.5, 1900, -0.5, 1899.5);
        fListHist->Add(f2fHistLGenstat);
        
        f2fHistAGenstat = new TH2F("f2fHistAGenstat","f2fHistAGenstat", 100, -0.5, 99.5, 1900, -0.5, 1899.5);
        fListHist->Add(f2fHistAGenstat);
        
        f2fHistXiPlus = new TH2F("f2fHistXiPlus","f2fHistXiPlus", xNbins, xBinEdge, xibinnumb, xibinlimits);
        fListHist->Add(f2fHistXiPlus);
        
        f2fHistXiMinus = new TH2F("f2fHistXiMinus","f2fHistXiMinus",xNbins, xBinEdge, xibinnumb, xibinlimits);
        fListHist->Add(f2fHistXiMinus);
        
        f2fHistLambdafromXi = new TH3F("f2fHistLambdafromXi","f2fHistLambdafromXi", fNptBins, LambdaPtBins,xNbins, xBinEdge, xibinnumb, xibinlimits);
        fListHist->Add(f2fHistLambdafromXi);
        
        f2fHistAntiLambdafromXi = new TH3F("f2fHistAntiLambdafromXi","f2fHistAntiLambdafromXi",fNptBins, LambdaPtBins,xNbins, xBinEdge, xibinnumb, xibinlimits);
        fListHist->Add(f2fHistAntiLambdafromXi);
        
        fTreeV0->Branch("fTreeVariablePID",&fTreeVariablePID);
        fTreeV0->Branch("fTreeVariablePIDPositive",&fTreeVariablePIDPositive);
        fTreeV0->Branch("fTreeVariablePIDNegative",&fTreeVariablePIDNegative);
        
    }
    
    PostData(1,fListHist);
    PostData(2,fTreeV0);
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskNetLambdaTrad::UserExec(Option_t *)
{
    
    const Int_t dim = fNptBins*2;
    Int_t ptCh[dim];
    Int_t ptChCut[dim];
    Int_t ptChMC[dim];
    
    for(Int_t idx = 0; idx < dim; idx++)
    {
        ptCh[idx] = 0.;
        ptChCut[idx] = 0.;
        ptChMC[idx] = 0;
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
    Double_t lMagneticField = -10;
    lMagneticField = fESD->GetMagneticField();
    
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
    
    const AliESDVertex *lPrimaryBestESDVtx     = fESD->GetPrimaryVertex();
    Double_t lBestPrimaryVtxPos[3]          = {-100.0, -100.0, -100.0};
    lPrimaryBestESDVtx->GetXYZ( lBestPrimaryVtxPos );
    
    
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
                f2fHistXiPlus->Fill(fCentrality,gpt);
            }
            
            if(genpid == 3312)
            {
                f2fHistXiMinus->Fill(fCentrality,gpt);
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
        fPtBinNplusNminusChTruth->Fill(ptContainerMC);
        
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
    
    Double_t fMinV0Pt = 0.0;
    Double_t fMaxV0Pt = 5.0;
    
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
        Float_t ppt = -999,  peta = -999, posprnsg = -999, v0Radius =-999;
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
            

            if(esdpTrack->Charge() == esdnTrack->Charge())
            {
                continue;
            }
            
            esdv0->GetXYZ(vertx[0], vertx[1], vertx[2]);//decay vertex
            v0Radius = TMath::Sqrt(vertx[0]*vertx[0]+vertx[1]*vertx[1]);
            V0pt = esdv0->Pt();
            if ((V0pt<fMinV0Pt)||(fMaxV0Pt<V0pt)) continue;
            
            ///////////////////////////////////////////////////////////////////////

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
            
            pmom = esdv0->P();
            eta = esdv0->Eta();
            ppt = esdpTrack->Pt();
            peta = esdpTrack->Eta();
            npt = esdnTrack->Pt();
            neta = esdnTrack->Eta();
            ontheflystat = esdv0->GetOnFlyStatus();
            
            dcaPosToVertex = TMath::Abs(esdpTrack->GetD(lBestPrimaryVtxPos[0],
                                                        lBestPrimaryVtxPos[1],
                                                        lMagneticField) );
            
            dcaNegToVertex = TMath::Abs(esdnTrack->GetD(lBestPrimaryVtxPos[0],
                                                        lBestPrimaryVtxPos[1],
                                                        lMagneticField) );
            
            cosPointingAngle = esdv0->GetV0CosineOfPointingAngle(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2]);
            dcaDaughters = esdv0->GetDcaV0Daughters();
            dcaV0ToVertex = esdv0->GetD(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2]);


            esdv0->ChangeMassHypothesis(3122);
            invMassLambda = esdv0->GetEffMass();
            esdv0->ChangeMassHypothesis(-3122);
            invMassAntiLambda = esdv0->GetEffMass();
            
            posprnsg = fPIDResponse->NumberOfSigmasTPC(esdpTrack, AliPID::kProton);
            negprnsg = fPIDResponse->NumberOfSigmasTPC(esdnTrack, AliPID::kProton);
            
        }
        
        fTreeVariableDcaV0ToPrimVertex = dcaV0ToVertex;
        fTreeVariableDcaV0Daughters = dcaDaughters;
        fTreeVariableDcaPosToPrimVertex = dcaPosToVertex;
        fTreeVariableDcaNegToPrimVertex = dcaNegToVertex;
        fTreeVariableInvMassLambda = invMassLambda;
        fTreeVariableInvMassAntiLambda = invMassAntiLambda;
        
        if(TMath::Abs(peta) > 1) continue;
        if(TMath::Abs(neta) > 1) continue;
        if(cosPointingAngle < 0.999) continue;
        if(dcaDaughters > 1.5) continue;
        if(v0Radius < 5.0) continue;
        if(v0Radius > 200.) continue;
        

        
        Int_t iptbin = GetPtBin(V0pt);
        
        if( ontheflystat == 0 )
        {
            
            if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.1 && dcaPosToVertex >  0.05 && TMath::Abs(posprnsg)  <= 3.)
            {
                f2fHistInvMassVsPtLambda->Fill(invMassLambda,V0pt);
                f2fHistRecCentVsPtLambda->Fill(fCentrality,V0pt);
                
                if(invMassLambda > 1.11 && invMassLambda < 1.122)
                {
                    f2fHistPtmassctLambda->Fill(fCentrality,V0pt);
                    f1fHistmassctLambda->Fill(invMassLambda);
                    ptChCut[iptbin] += 1;
                }
                
            }
            if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.05 && dcaPosToVertex >  0.1 && TMath::Abs(negprnsg)  <= 3.)
            {
                f2fHistInvMassVsPtAntiLambda->Fill(invMassAntiLambda,V0pt);
                f2fHistRecCentVsPtAntiLambda->Fill(fCentrality,V0pt);
                
                if(invMassAntiLambda > 1.11 && invMassAntiLambda < 1.122)
                {
                    f2fHistPtmassctAntiLambda->Fill(fCentrality,V0pt);
                    f1fHistmassctAntiLambda->Fill(invMassAntiLambda);
                    ptChCut[iptbin+fNptBins] += 1;
                }
            }
            
            
            if(fIsMC)
            {
                fTreeVariablePID = -999, fTreeVariablePIDParent = -999;
                Float_t mcpt = -999, mceta = -999, fTreeVariablePtParent = -999;
                Bool_t isPrim = kFALSE, isSecFromMaterial = kFALSE, isSecFromWeakDecay = kFALSE, isPrimParent =kFALSE;
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
                    
                    Int_t iptbinRecPRI = GetPtBin(mcpt);
                    
                    if(fTreeVariablePID == 3122)
                        {
                            if(isPrim){f2fHistRecPrimariesCentVsPtLambda->Fill(fCentrality,mcpt);ptCh[iptbinRecPRI] += 1;}
                            else{f2fHistLambdaMisId->Fill(fCentrality,mcpt);}
                            if(isSecFromWeakDecay) {f2fHistLambdaSecFromWeakDecay->Fill(fCentrality,mcpt);}
                            else if (isSecFromMaterial) {f2fHistLambdaMaterial->Fill(fCentrality,mcpt);}
                            f2fHistInvMassVsPtLambdaRec->Fill(invMassLambda,mcpt);
                            f2fHistLambdaRecPt->Fill(fCentrality,mcpt);
                            nRecL += 1.;
                        }
                    
                        if(fTreeVariablePID == -3122)
                        {
                            if(isPrim){f2fHistRecPrimariesCentVsPtAntiLambda->Fill(fCentrality,mcpt);ptCh[iptbinRecPRI+fNptBins] += 1;}
                            else{f2fHistAntiLambdaMisId->Fill(fCentrality,mcpt);}
                            if(isSecFromWeakDecay) {f2fHistAntiLambdaSecFromWeakDecay->Fill(fCentrality,mcpt);}
                            else if (isSecFromMaterial) {f2fHistAntiLambdaMaterial->Fill(fCentrality,mcpt);}
                            f2fHistInvMassVsPtAntiLambdaRec->Fill(invMassAntiLambda,mcpt);
                            f2fHistAntiLambdaRecPt->Fill(fCentrality,mcpt);
                            nRecA += 1.;
                            
                        }
                }
                else
                {
                    if(TMath::Abs(esdpTrack->GetLabel()) >= nGen || TMath::Abs(esdnTrack->GetLabel()) >= nGen) continue;
                    Int_t lblPosV0Dghter = (Int_t) TMath::Abs( esdpTrack->GetLabel() );
                    Int_t lblNegV0Dghter = (Int_t) TMath::Abs( esdnTrack->GetLabel() );
                    
                    TParticle* esdGenTrackPos = stack->Particle( lblPosV0Dghter );
                    TParticle* esdGenTrackNeg = stack->Particle( lblNegV0Dghter );
                    
                    Int_t lPIDPositive = esdGenTrackPos -> GetPdgCode();
                    Int_t lPIDNegative = esdGenTrackNeg -> GetPdgCode();
                    
                    fTreeVariablePIDPositive = lPIDPositive;
                    fTreeVariablePIDNegative = lPIDNegative;
                    
                    Int_t posTparticle = esdGenTrackPos->GetFirstMother();
                    Int_t negTparticle = esdGenTrackNeg->GetFirstMother();
                    
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
                        
                        Int_t esdlthisV0parent = esdlthisV0->GetFirstMother();
                        if (esdlthisV0parent > 0)
                        {
                            TParticle *lbV0parent = stack->Particle(esdlthisV0parent);
                            if(!lbV0parent) continue;
                            fTreeVariablePIDParent = lbV0parent->GetPdgCode();
                            fTreeVariablePtParent = lbV0parent->Pt();
                            isPrimParent =  stack->IsPhysicalPrimary(esdlthisV0parent);
                        }
                    }
                    
                    Int_t iptbinRecPRI = GetPtBin(mcpt);
                    if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.1 && dcaPosToVertex >  0.05 && TMath::Abs(posprnsg)  <= 3.)
                    {

                        if(fTreeVariablePID == 3122)
                        {
                            if(isPrim){f2fHistRecPrimariesCentVsPtLambda->Fill(fCentrality,mcpt);ptCh[iptbinRecPRI] += 1;}
                            else{f2fHistLambdaMisId->Fill(fCentrality,mcpt);}
                            if(isSecFromWeakDecay)
                            {
                                f2fHistLambdaSecFromWeakDecay->Fill(fCentrality,mcpt);
                                if(isPrimParent)
                                {
                                    if(fTreeVariablePIDParent == 3312)
                                    {
                                        f2fHistLambdafromXi->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                            else if (isSecFromMaterial) {f2fHistLambdaMaterial->Fill(fCentrality,mcpt);}
                            f2fHistInvMassVsPtLambdaRec->Fill(invMassLambda,mcpt);
                            f2fHistLambdaRecPt->Fill(fCentrality,mcpt);
                            nRecL += 1.;
                        }
                    }
                    if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.05 && dcaPosToVertex >  0.1 && TMath::Abs(negprnsg)  <= 3.)
                    {

                        if(fTreeVariablePID == -3122)
                        {
                            if(isPrim){f2fHistRecPrimariesCentVsPtAntiLambda->Fill(fCentrality,mcpt);ptCh[iptbinRecPRI+fNptBins] += 1;}
                            else{f2fHistAntiLambdaMisId->Fill(fCentrality,mcpt);}
                            if(isSecFromWeakDecay)
                            {
                                f2fHistAntiLambdaSecFromWeakDecay->Fill(fCentrality,mcpt);
                                if(isPrimParent)
                                {
                                    if(fTreeVariablePIDParent == -3312)
                                    {
                                    f2fHistAntiLambdafromXi->Fill(mcpt,fCentrality,fTreeVariablePtParent);
                                    }
                                }
                            }
                            else if (isSecFromMaterial) {f2fHistAntiLambdaMaterial->Fill(fCentrality,mcpt);}
                            f2fHistInvMassVsPtAntiLambdaRec->Fill(invMassAntiLambda,mcpt);
                            f2fHistAntiLambdaRecPt->Fill(fCentrality,mcpt);
                            nRecA += 1.;
                        }
                    }
                } //ESD
            } //MC condition
        }// zero onfly V0
    }// end of V0 loop
    
    
    f2fHistLRecstat->Fill(fCentrality, nRecL);
    f2fHistARecstat->Fill(fCentrality, nRecA);
    ////////////////////////////////////////
    Double_t ptContainer[dim+1];
    ptContainer[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainer[i] = ptCh[i-1];
    }
    fPtBinNplusNminusCh->Fill(ptContainer);
    /////////////////////////////////////////
    Double_t ptContainerCut[dim+1];
    ptContainerCut[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerCut[i] = ptChCut[i-1];
    }
    fPtBinNplusNminusChCut->Fill(ptContainerCut);
    /////////////////////////////////////////
    fTreeV0->Fill();
    PostData(1,fListHist);
    PostData(2,fTreeV0);
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Int_t AliAnalysisTaskNetLambdaTrad::GetPtBin(Double_t pt)
{
    Int_t bin = -1;
    
    Double_t LambdaPtBins[24] =  {0.5,0.7,0.9,1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0,4.2};
    
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


