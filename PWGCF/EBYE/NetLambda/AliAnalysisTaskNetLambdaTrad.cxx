
// For: Net Lambda fluctuation analysis via traditional method
// By: Ejiro Naomi Umaka Apr 2018
// Updated Mar 26


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
fHistEventCounter(0x0),
fHistCentrality(0x0),

/// LAMBDA
f3fHistCentVsInvMassLambda1point6(0x0),
f3fHistCentVsInvMassLambda1point0(0x0),
f3fHistCentVsInvMassLambda0point6(0x0),
f3fHistCentVsInvMassLambda0point2(0x0),
f3fHistPtmassctLambda1point6(0x0),
f3fHistPtmassctLambda1point0(0x0),
f3fHistPtmassctLambda0point6(0x0),
f3fHistPtmassctLambda0point2(0x0),
f3fHistPtmassctLambdaPosOpoint2(0x0),
f3fHistPtmassctLambdaPosOpoint4(0x0),
f3fHistPtmassctLambdaPosOpoint6(0x0),
f3fHistPtmassctLambdaPosOpoint8(0x0),
///ANTI-LAMBDA
f3fHistCentVsInvMassAntiLambda1point6(0x0),
f3fHistCentVsInvMassAntiLambda1point0(0x0),
f3fHistCentVsInvMassAntiLambda0point6(0x0),
f3fHistCentVsInvMassAntiLambda0point2(0x0),
f3fHistPtmassctAntiLambda1point6(0x0),
f3fHistPtmassctAntiLambda1point0(0x0),
f3fHistPtmassctAntiLambda0point6(0x0),
f3fHistPtmassctAntiLambda0point2(0x0),
f3fHistPtmassctAntiLambdaPosOpoint2(0x0),
f3fHistPtmassctAntiLambdaPosOpoint4(0x0),
f3fHistPtmassctAntiLambdaPosOpoint6(0x0),
f3fHistPtmassctAntiLambdaPosOpoint8(0x0),
///
fCentrality(-1),
fNptBins(31),

fEvSel(AliVEvent::kINT7),

fPtBinNplusNminusChEtaFour(NULL),
fPtBinNplusNminusChEtaThree(NULL),
fPtBinNplusNminusChEtaTwo(NULL),
fPtBinNplusNminusChEtaOne(NULL),

fPtBinNplusNminusChPosEtaFour(NULL),
fPtBinNplusNminusChPosEtaThree(NULL),
fPtBinNplusNminusChPosEtaTwo(NULL),
fPtBinNplusNminusChPosEtaOne(NULL),

fPtBinNplusNminusChBproxy(NULL)

{
    Info("AliAnalysisTaskNetLambdaTrad","Calling Constructor");
    
    DefineInput(0,TChain::Class());
    DefineOutput(1,TList::Class());
    
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
 //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //pt binning
    Double_t LambdaPtBins[32] = {1.0,1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1};
    //---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


    //V0 hists//
    //-----------------------------------------------------------------------LAMBDA-----------------------------------------------------------------------------------------------------------------------------------------

    f3fHistCentVsInvMassLambda1point6 = new TH3F("f3fHistCentVsInvMassLambda1point6","Cent vs. #Lambda Inv Mass vs. pT(deltaEta 1.6)",81,-0.5,80.5,100,1.08,1.16,40,1,5); //1.0-4.9
    fListHist->Add(f3fHistCentVsInvMassLambda1point6);
    
    f3fHistCentVsInvMassLambda1point0 = new TH3F("f3fHistCentVsInvMassLambda1point0","Cent vs. #Lambda Inv Mass vs. pT(deltaEta 1.0)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
    fListHist->Add(f3fHistCentVsInvMassLambda1point0);
    
    f3fHistCentVsInvMassLambda0point6 = new TH3F("f3fHistCentVsInvMassLambda0point6","Cent vs. #Lambda Inv Mass vs. pT(deltaEta 0.6)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
    fListHist->Add(f3fHistCentVsInvMassLambda0point6);
    
    f3fHistCentVsInvMassLambda0point2 = new TH3F("f3fHistCentVsInvMassLambda0point2","Cent vs. #Lambda Inv Mass vs. pT(deltaEta 0.2)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
    fListHist->Add(f3fHistCentVsInvMassLambda0point2);
    ////
    
    f3fHistPtmassctLambda1point6 = new TH3F("f3fHistPtmassctLambda1point6","Cent vs. #Lambda Inv Mass vs. pT(deltaEta 1.6)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
    fListHist->Add(f3fHistPtmassctLambda1point6);
    
    f3fHistPtmassctLambda1point0 = new TH3F("f3fHistPtmassctLambda1point0","Cent vs. #Lambda Inv Mass vs. pT(deltaEta 1.0)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
    fListHist->Add(f3fHistPtmassctLambda1point0);
    
    f3fHistPtmassctLambda0point6 = new TH3F("f3fHistPtmassctLambda0point6","Cent vs. #Lambda Inv Mass vs. pT(deltaEta 0.6)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
    fListHist->Add(f3fHistPtmassctLambda0point6);
    
    f3fHistPtmassctLambda0point2 = new TH3F("f3fHistPtmassctLambda0point2","Cent vs. #Lambda Inv Mass vs. pT(deltaEta 0.2)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
    fListHist->Add(f3fHistPtmassctLambda0point2);
    ////
    
    f3fHistPtmassctLambdaPosOpoint2 = new TH3F("f3fHistPtmassctLambdaPosOpoint2","Cent vs. #Lambda Inv Mass vs. pT(Eta 0-0.2)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
    fListHist->Add(f3fHistPtmassctLambdaPosOpoint2);
    
    f3fHistPtmassctLambdaPosOpoint4 = new TH3F("f3fHistPtmassctLambdaPosOpoint4","Cent vs. #Lambda Inv Mass vs. pT(Eta 0.2-0.4)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
    fListHist->Add(f3fHistPtmassctLambdaPosOpoint4);
    
    f3fHistPtmassctLambdaPosOpoint6 = new TH3F("f3fHistPtmassctLambdaPosOpoint6","Cent vs. #Lambda Inv Mass vs. pT(Eta 0.4-0.6)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
    fListHist->Add(f3fHistPtmassctLambdaPosOpoint6);
    
    f3fHistPtmassctLambdaPosOpoint8 = new TH3F("f3fHistPtmassctLambdaPosOpoint8","Cent vs. #Lambda Inv Mass vs. pT(Eta 0.6-0.8)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
    fListHist->Add(f3fHistPtmassctLambdaPosOpoint8);
    
     //-----------------------------------------------------------------------ANTI-LAMBDA-----------------------------------------------------------------------------------------------------------------------------------------

    f3fHistCentVsInvMassAntiLambda1point6 = new TH3F("f3fHistCentVsInvMassAntiLambda1point6","Cent vs. #bar{#Lambda} Inv Mass vs. pT(deltaEta 1.6)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
    fListHist->Add(f3fHistCentVsInvMassAntiLambda1point6);
    
    f3fHistCentVsInvMassAntiLambda1point0 = new TH3F("f3fHistCentVsInvMassAntiLambda1point0","Cent vs. #bar{#Lambda} Inv Mass vs. pT(deltaEta 1.0)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
    fListHist->Add(f3fHistCentVsInvMassAntiLambda1point0);
    
    f3fHistCentVsInvMassAntiLambda0point6 = new TH3F("f3fHistCentVsInvMassAntiLambda0point6","Cent vs. #bar{#Lambda} Inv Mass vs. pT(deltaEta 0.6)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
    fListHist->Add(f3fHistCentVsInvMassAntiLambda0point6);
    
    f3fHistCentVsInvMassAntiLambda0point2 = new TH3F("f3fHistCentVsInvMassAntiLambda0point2","Cent vs. #bar{#Lambda} Inv Mass vs. pT(deltaEta 0.2)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
    fListHist->Add(f3fHistCentVsInvMassAntiLambda0point2);
    ////
    
    f3fHistPtmassctAntiLambda1point6 = new TH3F("f3fHistPtmassctAntiLambda1point6","Cent vs. #bar{#Lambda} Inv Mass vs. pT(deltaEta 1.6)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
    fListHist->Add(f3fHistPtmassctAntiLambda1point6);
    
    f3fHistPtmassctAntiLambda1point0 = new TH3F("f3fHistPtmassctAntiLambda1point0","Cent vs. #bar{#Lambda} Inv Mass vs. pT(deltaEta 1.0)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
    fListHist->Add(f3fHistPtmassctAntiLambda1point0);
    
    f3fHistPtmassctAntiLambda0point6 = new TH3F("f3fHistPtmassctAntiLambda0point6","Cent vs. #bar{#Lambda} Inv Mass vs. pT(deltaEta 0.6)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
    fListHist->Add(f3fHistPtmassctAntiLambda0point6);
    
    f3fHistPtmassctAntiLambda0point2 = new TH3F("f3fHistPtmassctAntiLambda0point2","Cent vs. #bar{#Lambda} Inv Mass vs. pT(deltaEta 0.2)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
    fListHist->Add(f3fHistPtmassctAntiLambda0point2);
    ////
    
    f3fHistPtmassctAntiLambdaPosOpoint2 = new TH3F("f3fHistPtmassctAntiLambdaPosOpoint2","Cent vs. #bar{#Lambda} Inv Mass vs. pT(Eta 0-0.2)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
    fListHist->Add(f3fHistPtmassctAntiLambdaPosOpoint2);
    
    f3fHistPtmassctAntiLambdaPosOpoint4 = new TH3F("f3fHistPtmassctAntiLambdaPosOpoint4","Cent vs. #bar{#Lambda} Inv Mass vs. pT(Eta 0.2-0.4)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
    fListHist->Add(f3fHistPtmassctAntiLambdaPosOpoint4);
    
    f3fHistPtmassctAntiLambdaPosOpoint6 = new TH3F("f3fHistPtmassctAntiLambdaPosOpoint6","Cent vs. #bar{#Lambda} Inv Mass vs. pT(Eta 0.4-0.6)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
    fListHist->Add(f3fHistPtmassctAntiLambdaPosOpoint6);
    
    f3fHistPtmassctAntiLambdaPosOpoint8 = new TH3F("f3fHistPtmassctAntiLambdaPosOpoint8","Cent vs. #bar{#Lambda} Inv Mass vs. pT(Eta 0.6-0.8)",81,-0.5,80.5,100,1.08,1.16,40,1,5);
    fListHist->Add(f3fHistPtmassctAntiLambdaPosOpoint8);
    
    //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    //THNSPARSE BINNING
    const Int_t dim = 63; //31 pt bins*2 + 1 cent bin
    Int_t bin[dim];
    bin[0] = 81;
    for(Int_t ibin = 1; ibin < dim; ibin++) bin[ibin] = 300;
    Double_t min[dim];
    for(Int_t jbin = 0; jbin < dim; jbin++) min[jbin] =  -0.5;
    Double_t max[dim];
    max[0] = 80.5;
    for(Int_t jbin = 1; jbin < dim; jbin++) max[jbin] = 299.5;

    
    fPtBinNplusNminusChEtaFour = new THnSparseI("fPtBinNplusNminusChEtaFour","cent-nlambda-nantilambda masscut", dim, bin, min, max);
    fListHist->Add(fPtBinNplusNminusChEtaFour);
    fPtBinNplusNminusChEtaThree = new THnSparseI("fPtBinNplusNminusChEtaThree","cent-nlambda-nantilambda masscut", dim, bin, min, max);
    fListHist->Add(fPtBinNplusNminusChEtaThree);
    fPtBinNplusNminusChEtaTwo = new THnSparseI("fPtBinNplusNminusChEtaTwo","cent-nlambda-nantilambda masscut", dim, bin, min, max);
    fListHist->Add(fPtBinNplusNminusChEtaTwo);
    fPtBinNplusNminusChEtaOne = new THnSparseI("fPtBinNplusNminusChEtaOne","cent-nlambda-nantilambda masscut", dim, bin, min, max);
    fListHist->Add(fPtBinNplusNminusChEtaOne);
    
    fPtBinNplusNminusChPosEtaFour = new THnSparseI("fPtBinNplusNminusChPosEtaFour","cent-nlambda-nantilambda masscut", dim, bin, min, max);
    fListHist->Add(fPtBinNplusNminusChPosEtaFour);
    fPtBinNplusNminusChPosEtaThree = new THnSparseI("fPtBinNplusNminusChPosEtaThree","cent-nlambda-nantilambda masscut", dim, bin, min, max);
    fListHist->Add(fPtBinNplusNminusChPosEtaThree);
    fPtBinNplusNminusChPosEtaTwo = new THnSparseI("fPtBinNplusNminusChPosEtaTwo","cent-nlambda-nantilambda masscut", dim, bin, min, max);
    fListHist->Add(fPtBinNplusNminusChPosEtaTwo);
    fPtBinNplusNminusChPosEtaOne = new THnSparseI("fPtBinNplusNminusChPosEtaOne","cent-nlambda-nantilambda masscut", dim, bin, min, max);
    fListHist->Add(fPtBinNplusNminusChPosEtaOne);
    
    fPtBinNplusNminusChBproxy = new THnSparseI("fPtBinNplusNminusChBproxy","cent-nlambda-nantilambda masscut", dim, bin, min, max);
    fListHist->Add(fPtBinNplusNminusChBproxy);
    
    PostData(1,fListHist);
}


//------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskNetLambdaTrad::UserExec(Option_t *)
{
    
    const Int_t dim = fNptBins*2;
    
    // delta eta
    Int_t ptChEta0point2[dim];
    Int_t ptChEta0point6[dim];
    Int_t ptChEta1point0[dim];
    Int_t ptChEta1point6[dim];

    // postive eta
    Int_t ptChEtaPosOpoint2[dim];
    Int_t ptChEtaPosOpoint4[dim];
    Int_t ptChEtaPosOpoint6[dim];
    Int_t ptChEtaPosOpoint8[dim];

    //background proxy
    Int_t ptChBproxy[dim];
    
    for(Int_t idx = 0; idx < dim; idx++)
    {
        ptChEta0point2[idx] = 0.;
        ptChEta0point6[idx] = 0.;
        ptChEta1point0[idx] = 0.;
        ptChEta1point6[idx] = 0.;
        ptChEtaPosOpoint2[idx] = 0.;
        ptChEtaPosOpoint4[idx] = 0.;
        ptChEtaPosOpoint6[idx] = 0.;
        ptChEtaPosOpoint8[idx] = 0.;

        ptChBproxy[idx] = 0.;
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
    Double_t fMinV0Pt = 0.5;
    Double_t fMaxV0Pt = 4.5;
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
        Float_t ppt = -999,  peta = -999, posprnsg = -999, pospion =-999, v0Radius =-999, lRapLambda=-999;
        Float_t npt = -999,  neta = -999, negprnsg = -999, negpion =-999;
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
        
        lRapLambda  = esdv0->RapLambda();
        V0pt = esdv0->Pt();
        if ((V0pt<fMinV0Pt)||(fMaxV0Pt<V0pt)) continue;
        if(TMath::Abs(lRapLambda)> 0.5 ) continue;
        
        
        //--------------------------------------------------------------------Track selection-------------------------------------------------------------------------------------------------------------------------------------

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
        //----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

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
        pospion  = fPIDResponse->NumberOfSigmasTPC( esdpTrack, AliPID::kPion );
        negpion  = fPIDResponse->NumberOfSigmasTPC( esdnTrack, AliPID::kPion );
        
        if(TMath::Abs(peta) > 1) continue;
        if(TMath::Abs(neta) > 1) continue;
        if(cosPointingAngle < 0.99) continue;
        if(dcaDaughters > 0.8) continue;
        if(v0Radius < 5.0) continue;
        if(v0Radius > 200.) continue;
        
        Int_t iptbin = GetPtBin(V0pt);
        
        if( ontheflystat == 0 )
        {
            if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.25 && dcaPosToVertex >  0.1  && TMath::Abs(posprnsg)  <= 3. && TMath::Abs(negpion)  <= 3.)
            { //fit to estimate lambda background contribution
                
                if(TMath::Abs(eta) < 0.8)
                {
                    f3fHistCentVsInvMassLambda1point6->Fill(fCentrality,invMassLambda,V0pt);
                }
                if(TMath::Abs(eta) < 0.5)
                {
                    f3fHistCentVsInvMassLambda1point0->Fill(fCentrality,invMassLambda,V0pt);
                }
                if(TMath::Abs(eta) < 0.3)
                {
                    f3fHistCentVsInvMassLambda0point6->Fill(fCentrality,invMassLambda,V0pt);
                }
                if(TMath::Abs(eta) < 0.1)
                {
                    f3fHistCentVsInvMassLambda0point2->Fill(fCentrality,invMassLambda,V0pt);
                }
                
                if(invMassLambda > 1.11 && invMassLambda < 1.122) //count cumulants
                {
                    if(TMath::Abs(eta) < 0.8)
                    {
                        f3fHistPtmassctLambda1point6->Fill(fCentrality,invMassLambda,V0pt);
                        ptChEta1point6[iptbin] += 1;
                    }
                    if(TMath::Abs(eta) < 0.5)
                    {
                        f3fHistPtmassctLambda1point0->Fill(fCentrality,invMassLambda,V0pt);
                        ptChEta1point0[iptbin] += 1;
                    }
                    if(TMath::Abs(eta) < 0.3)
                    {
                        f3fHistPtmassctLambda0point6->Fill(fCentrality,invMassLambda,V0pt);
                        ptChEta0point6[iptbin] += 1;
                    }
                    if(TMath::Abs(eta) < 0.1)
                    {
                        f3fHistPtmassctLambda0point2->Fill(fCentrality,invMassLambda,V0pt);
                        ptChEta0point2[iptbin] += 1;
                     }
                    
                    if(eta > 0.0 && eta < 0.2) //shift eta in positive directon only
                    {
                        f3fHistPtmassctLambdaPosOpoint2->Fill(fCentrality,invMassLambda,V0pt);
                        ptChEtaPosOpoint2[iptbin] += 1;
                    }
                    if(eta > 0.2 && eta < 0.4)
                    {
                        f3fHistPtmassctLambdaPosOpoint4->Fill(fCentrality,invMassLambda,V0pt);
                        ptChEtaPosOpoint4[iptbin] += 1;
                    }
                    if(eta > 0.4 && eta < 0.6)
                    {
                        f3fHistPtmassctLambdaPosOpoint6->Fill(fCentrality,invMassLambda,V0pt);
                        ptChEtaPosOpoint6[iptbin] += 1;
                    }
                    if(eta > 0.6 && eta < 0.8)
                    {
                        f3fHistPtmassctLambdaPosOpoint8->Fill(fCentrality,invMassLambda,V0pt);
                        ptChEtaPosOpoint8[iptbin] += 1;
                    }
                }
                if(invMassLambda > 1.09 && invMassLambda < 1.1) //Bproxy
                {
                    if(TMath::Abs(eta) < 0.8)
                    ptChBproxy[iptbin] += 1;
                }
            }
            if(dcaV0ToVertex < 0.25 && dcaNegToVertex > 0.1 && dcaPosToVertex >  0.25 && TMath::Abs(negprnsg)  <= 3. && TMath::Abs(pospion)  <= 3.)
            {
                { //fit to estimate anti-lambda background contribution
                    
                    if(TMath::Abs(eta) < 0.8)
                    {
                        f3fHistCentVsInvMassAntiLambda1point6->Fill(fCentrality,invMassAntiLambda,V0pt);
                    }
                    if(TMath::Abs(eta) < 0.5)
                    {
                        f3fHistCentVsInvMassAntiLambda1point0->Fill(fCentrality,invMassAntiLambda,V0pt);
                    }
                    if(TMath::Abs(eta) < 0.3)
                    {
                        f3fHistCentVsInvMassAntiLambda0point6->Fill(fCentrality,invMassAntiLambda,V0pt);
                    }
                    if(TMath::Abs(eta) < 0.1)
                    {
                        f3fHistCentVsInvMassAntiLambda0point2->Fill(fCentrality,invMassAntiLambda,V0pt);
                    }
                    
                    if(invMassAntiLambda > 1.11 && invMassAntiLambda < 1.122) //count cumulants
                    {
                        if(TMath::Abs(eta) < 0.8)
                        {
                            f3fHistPtmassctAntiLambda1point6->Fill(fCentrality,invMassAntiLambda,V0pt);
                            ptChEta1point6[iptbin+fNptBins] += 1;
                        }
                        if(TMath::Abs(eta) < 0.5)
                        {
                            f3fHistPtmassctAntiLambda1point0->Fill(fCentrality,invMassAntiLambda,V0pt);
                            ptChEta1point0[iptbin+fNptBins] += 1;
                        }
                        if(TMath::Abs(eta) < 0.3)
                        {
                            f3fHistPtmassctAntiLambda0point6->Fill(fCentrality,invMassAntiLambda,V0pt);
                            ptChEta0point6[iptbin+fNptBins] += 1;
                        }
                        if(TMath::Abs(eta) < 0.1)
                        {
                            f3fHistPtmassctAntiLambda0point2->Fill(fCentrality,invMassAntiLambda,V0pt);
                            ptChEta0point2[iptbin+fNptBins] += 1;
                        }
                        
                        if(eta > 0.0 && eta < 0.2) //shift eta in positive directon only
                        {
                            f3fHistPtmassctAntiLambdaPosOpoint2->Fill(fCentrality,invMassAntiLambda,V0pt);
                            ptChEtaPosOpoint2[iptbin+fNptBins] += 1;
                        }
                        if(eta > 0.2 && eta < 0.4)
                        {
                            f3fHistPtmassctAntiLambdaPosOpoint4->Fill(fCentrality,invMassAntiLambda,V0pt);
                            ptChEtaPosOpoint4[iptbin+fNptBins] += 1;
                        }
                        if(eta > 0.4 && eta < 0.6)
                        {
                            f3fHistPtmassctAntiLambdaPosOpoint6->Fill(fCentrality,invMassAntiLambda,V0pt);
                            ptChEtaPosOpoint6[iptbin+fNptBins] += 1;
                        }
                        if(eta > 0.6 && eta < 0.8)
                        {
                            f3fHistPtmassctAntiLambdaPosOpoint8->Fill(fCentrality,invMassAntiLambda,V0pt);
                            ptChEtaPosOpoint8[iptbin+fNptBins] += 1;
                        }
                    }
                    if(invMassAntiLambda > 1.09 && invMassAntiLambda < 1.1) //Bproxy
                    {
                        if(TMath::Abs(eta) < 0.8)
                        ptChBproxy[iptbin+fNptBins] += 1;
                    }
                }
            }//lambdarbar selection
        }// zero onfly V0
    }// end of V0 loop
    

    Double_t ptContainerFour[dim+1];  //4
    ptContainerFour[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerFour[i] = ptChEta1point6[i-1];
    }
    fPtBinNplusNminusChEtaFour->Fill(ptContainerFour);
    //-------------------------------------------------
    Double_t ptContainerThree[dim+1];//3
    ptContainerThree[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerThree[i] = ptChEta1point0[i-1];
    }
    fPtBinNplusNminusChEtaThree->Fill(ptContainerThree);
    //-------------------------------------------------
    Double_t ptContainerTwo[dim+1];//2
    ptContainerTwo[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerTwo[i] = ptChEta0point6[i-1];
    }
    fPtBinNplusNminusChEtaTwo->Fill(ptContainerTwo);
    //-------------------------------------------------
    Double_t ptContainerOne[dim+1];//1
    ptContainerOne[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerOne[i] = ptChEta0point2[i-1];
    }
    fPtBinNplusNminusChEtaOne->Fill(ptContainerOne);
    //-------------------------------------------------
    //pos eta
    
    Double_t ptContainerPosFour[dim+1];  //0.6-0.8
    ptContainerPosFour[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerPosFour[i] = ptChEtaPosOpoint8[i-1];
    }
    fPtBinNplusNminusChPosEtaFour->Fill(ptContainerPosFour);
    //-------------------------------------------------
    Double_t ptContainerPosThree[dim+1];//0.4-0.6
    ptContainerPosThree[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerPosThree[i] = ptChEtaPosOpoint6[i-1];
    }
    fPtBinNplusNminusChPosEtaThree->Fill(ptContainerPosThree);
    //-------------------------------------------------
    Double_t ptContainerPosTwo[dim+1];//0.2-0.4
    ptContainerPosTwo[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerPosTwo[i] = ptChEtaPosOpoint4[i-1];
    }
    fPtBinNplusNminusChPosEtaTwo->Fill(ptContainerPosTwo);
    //-------------------------------------------------
    Double_t ptContainerPosOne[dim+1];//0.0-0.2
    ptContainerPosOne[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerPosOne[i] = ptChEtaPosOpoint2[i-1];
    }
    fPtBinNplusNminusChPosEtaOne->Fill(ptContainerPosOne);
    //-------------------------------------------------
    Double_t ptContainerBproxy[dim+1];
    ptContainerBproxy[0] = (Double_t)fCentrality;
    for(Int_t i = 1; i <= dim; i++)
    {
        ptContainerBproxy[i] = ptChBproxy[i-1];
    }
    fPtBinNplusNminusChBproxy->Fill(ptContainerBproxy);
    //-------------------------------------------------
    PostData(1,fListHist);
}


//-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
Int_t AliAnalysisTaskNetLambdaTrad::GetPtBin(Double_t pt)
{
    Int_t bin = -1;
    
    Double_t LambdaPtBins[32] = {1.0,1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1};

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


