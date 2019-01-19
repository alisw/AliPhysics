
// For: Net Lambda fluctuation analysis via traditional method
// By: Ejiro Naomi Umaka Apr 2018
// Updated Jan 11


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
fPIDResponse(0x0),
fEventCuts(0),
fListHist(0x0),
fTreeV0(0x0),
fHistEventCounter(0x0),
fHistCentrality(0x0),
f2fHistRecCentVsPtLambda(0x0),
f2fHistRecCentVsPtAntiLambda(0x0),
f2fHistInvMassVsPtLambda(0x0),
f2fHistInvMassVsPtAntiLambda(0x0),
f1fHistmassctLambda(0x0),
f1fHistmassctAntiLambda(0x0),
f2fHistPtmassctLambda(0x0),
f2fHistPtmassctAntiLambda(0x0),
fCentrality(-1),
fNptBins(23),
fEvSel(AliVEvent::kINT7),
fTreeVariableInvMassLambda(0),
fTreeVariableInvMassAntiLambda(0),
fTreeVariableDcaV0Daughters(0),
fTreeVariableDcaV0ToPrimVertex(0),
fTreeVariableDcaPosToPrimVertex(0),
fTreeVariableDcaNegToPrimVertex(0),

fPtBinNplusNminusChALL(NULL),
fPtBinNplusNminusChCut(NULL)

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
    
    fPtBinNplusNminusChALL = new THnSparseI("fPtBinNplusNminusChALL","cent-nlambda-nantilambda masscut", dim, bin, min, max);
    fListHist->Add(fPtBinNplusNminusChALL);
    
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
    
    
    
    PostData(1,fListHist);
    PostData(2,fTreeV0);
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskNetLambdaTrad::UserExec(Option_t *)
{
    
    const Int_t dim = fNptBins*2;
    Int_t ptCh[dim];
    Int_t ptChCut[dim];
    
    for(Int_t idx = 0; idx < dim; idx++)
    {
        ptCh[idx] = 0.;
        ptChCut[idx] = 0.;
    }
    
    
    if (!fInputEvent) return;
    
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if (!fESD) return;
    
    
    fPIDResponse = fInputHandler->GetPIDResponse();
    if(!fPIDResponse) return;
    
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
    
    
    
    Int_t nV0 = 0;
    Double_t fMinV0Pt = 0.0;
    Double_t fMaxV0Pt = 5.0;
    nV0 = fESD->GetNumberOfV0s();
    AliESDv0 *esdv0 = 0x0;
    AliESDtrack *esdpTrack = 0x0;
    AliESDtrack *esdnTrack = 0x0;
    
    for(Int_t iV0 = 0; iV0 < nV0; iV0++)
    {
        esdv0 = 0x0;
        esdpTrack = 0x0;
        esdnTrack = 0x0;
        
        Double_t  vertx[3];
        
        Float_t invMassLambda = -999, invMassAntiLambda = -999;
        Float_t V0pt = -999, eta = -999, pmom = -999;
        Float_t ppt = -999,  peta = -999, posprnsg = -999, v0Radius =-999;
        Float_t npt = -999,  neta = -999, negprnsg = -999;
        Bool_t  ontheflystat = kFALSE;
        Float_t dcaPosToVertex = -999, dcaNegToVertex = -999, dcaDaughters = -999, dcaV0ToVertex = -999, cosPointingAngle = -999;
        
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
        esdv0->GetXYZ(vertx[0], vertx[1], vertx[2]); //decay vertex
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
        //          dcaV0ToVertex = esdv0->GetD(vVtx[0],vVtx[1],vVtx[2]);
        dcaV0ToVertex = esdv0->GetD(lBestPrimaryVtxPos[0],lBestPrimaryVtxPos[1],lBestPrimaryVtxPos[2]);
        
        esdv0->ChangeMassHypothesis(3122);
        invMassLambda = esdv0->GetEffMass();
        esdv0->ChangeMassHypothesis(-3122);
        invMassAntiLambda = esdv0->GetEffMass();
        
        posprnsg = fPIDResponse->NumberOfSigmasTPC(esdpTrack, AliPID::kProton);
        negprnsg = fPIDResponse->NumberOfSigmasTPC(esdnTrack, AliPID::kProton);
        
        
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
                    ptCh[iptbin] += 1;
                }
                
                if(invMassLambda > 1.11341197 && invMassLambda < 1.11885) //1sigma
                {
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
                    ptCh[iptbin+fNptBins] += 1;
                }
                if(invMassAntiLambda > 1.11341 && invMassAntiLambda < 1.11887) //1sigma
                {
                    ptChCut[iptbin+fNptBins] += 1;
                }
            }
            
        }// zero onfly V0
    }// end of V0 loop
    
    Double_t ptContainer[dim+1];
    ptContainer[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainer[i] = ptCh[i-1];
    }
    fPtBinNplusNminusChALL->Fill(ptContainer);
    ////////////////////
    Double_t ptContainerCut[dim+1];
    ptContainerCut[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerCut[i] = ptChCut[i-1];
    }
    fPtBinNplusNminusChCut->Fill(ptContainerCut);
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


