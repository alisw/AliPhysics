#include "TList.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH3D.h"
#include "TCanvas.h"
#include "TF1.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliVVertex.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AlidNdPtEventCuts.h"
#include "AliEventCuts.h"
/*
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODTrack.h"
#include "AliAODtrackCuts.h"
#include "AliAODTrackSelection.h"
*/

#include "AliPhysicsSelection.h"
#include "AliVMultiplicity.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"
#include "AliMultSelectionTask.h"

#include "AliMCEvent.h"
#include "AliStack.h"
// #include "AliMCSpectraWeights.h"

#include "AlidNdPtUnifiedAnalysisTask.h"


/// \cond CLASSIMP
ClassImp(AlidNdPtUnifiedAnalysisTask);
/// \endcond

//________________________________________________________________________
AlidNdPtUnifiedAnalysisTask::AlidNdPtUnifiedAnalysisTask(const char *name) : AliAnalysisTaskSE(name),
    //General member variables
    fOutputList(0),
    fEvent(0),
    fMCEvent(0),
    fMCStack(0),
    fEventCuts(0),
    fAliEventCuts(0),
    fESDtrackCuts(0),
    //fAODtrackCuts(0),
    fUtils(0),
    //fMCSpectraWeights(0),
    //fstCollisionSystem("pp"),
    //fstMCTrainOutput(""),
    //fHistMCPartCorr(0),
    //Toggles
    fIsESD(kTRUE),
    fUseMultiplicity(kTRUE),
    fIsMC(kFALSE),
    fIs2013pA(kFALSE),
    fIs2015data(kFALSE),
    fUseTOFBunchCrossing(kFALSE),
    fTPCRefit(kFALSE),
    fITSRefit(kFALSE),
    fAcceptKinks(kTRUE),
    fRequiresClusterITS(kTRUE),
    fDCAToVertex2D(kFALSE),
    fSigmaToVertex(kFALSE),
    fUseGeomCut(kFALSE),
    //Event-Histograms
    fEventCount(0),
    fHistEvent(0),
    fHistMCGenINEL0Event(0),
    fHistMCRecINEL0Event(0),
    fHistMCTrigINEL0Event(0),
    fHistMCTrigEvent(0),
    //Track-Histograms
    fHistTrack(0),
    fHistCentCorrelpt(0),
    fHistCentCorrelBeforeCuts(0),
    fHistCentCorrelAfterCuts(0),
    fDCAyEtaPt(0),
    fDCAyEtaPtMCPrim(0),
    fDCAyEtaPtMCSecDecays(0),
    fDCAyEtaPtMCSecMaterial(0),
    fDCAyEtaPtMCSecDecaysK0s(0),
    fDCAyEtaPtMCSecDecaysLambda(0),
    fHistMCRecTrack(0),
    fHistMCGenPrimTrack(0),
    fHistMCRecPrimTrack(0),
    fHistMCRecSecTrack(0),
    fHistRelPtResoFromCov(0),
    // Cut Parameters
    fTriggerMask(AliVEvent::kMB),
    fMinEta(-10),
    fMaxEta(10),
    fMinPt(0),
    fMaxPt(999),
    fBinsDCA(100),
    fDCAbinEdges(1),
    fSigmaMeanXYZv(),
    fMeanXYZv(),
    fZvtx(10),
    fMinNCrossedRowsTPC(0),
    fMinRatioCrossedRowsOverFindableClustersTPC(0),
    fMaxFractionSharedClustersTPC(0),
    fMaxChi2PerTPCCluster(0),
    fMaxChi2PerITSCluster(0),
    fMaxDCAzITSTPC(0),
    fDCAToVertexXYPtDep("0"),
    fDCAToVertexXY(0),
    fMaxChi2TPCConstrained(0),
    fMinActiveLength(0),
    fDeadZoneWidth(2),
    fCutGeoNcrNclLenght(130),
    fCutGeoNcrNclGeom1Pt(1.5),
    fCutGeoNcrNclFractionNcr(0.85),
    fCutGeoNcrNclFractionNcl(0.7),
    fstCentEst("V0M"),
    //Arrays for Binning
    fBinsMultCent(0),
    fBinsPt(0),
    fBinsEta(0),
    fBinsZv(0),
    fBinsPtReso(0)
{
    // Set default binning
    Double_t binsMultCentDefault[2] = {0,10000};
    Double_t binsPtDefault[69] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.5,5.0,5.5,6.0,6.5,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0,16.0,18.0,20.0,22.0,24.0,26.0,28.0,30.0,32.0,34.0,36.0,40.0,45.0,50.0};
    Double_t binsEtaDefault[31] = {-1.5,-1.4,-1.3,-1.2,-1.1,-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5};
    Double_t binsZvDefault[13] = {-30.,-25.,-20.,-15.,-10.,-5.,0.,5.,10.,15.,20.,25.,30.};
    SetBinsPt(68,binsPtDefault);
    SetBinsEta(30,binsEtaDefault);
    SetBinsMultCent(1,binsMultCentDefault);
    SetBinsZv(12,binsZvDefault);
    SetMeanXYZv(0.0,0.0,0.0);
    SetSigmaMeanXYZv(1.0,1.0,10.0);
    SetZvtx(10.);
    // Output slot #0 id reserved by the base class for AOD
    // Output slot #1 writes into a TH1 container
    DefineOutput(1, TList::Class());
}

//________________________________________________________________________
void AlidNdPtUnifiedAnalysisTask::UserCreateOutputObjects(){
    // Create histograms here (function is called once)
    OpenFile(1,"recreate");
    fOutputList = new TList();
    fOutputList -> SetOwner();

    Int_t binsEventCount[4]={2,2,2,2};
    Double_t minEventCount[4]={0,0,0,0};
    Double_t maxEventCount[4]={2,2,2,2};
    Int_t ptNbins = fBinsPt->GetSize()-1;
    Double_t ptMin = fBinsPt->At(0);
    Double_t ptMax = fBinsPt->At(fBinsPt->GetSize());
    // set DCAy bins
    // for ITSTPC from -1 to 1 || for TPC-Only from -3 to 3
    Int_t DCAybins = fBinsDCA;
    Double_t DCAyMin = -fDCAbinEdges;
    Double_t DCAyMax = fDCAbinEdges;

    /// Standard track histogram pt:eta:zV:multcent
    Int_t nBinsTrack[4]={fBinsPt->GetSize()-1,fBinsEta->GetSize()-1,fBinsZv->GetSize()-1,fBinsMultCent->GetSize()-1};
    Double_t minTrack[4]={fBinsPt->GetAt(0),fBinsEta->GetAt(0),fBinsZv->GetAt(0),fBinsMultCent->GetAt(0)};
    Double_t maxTrack[4]={fBinsPt->GetAt(fBinsPt->GetSize()-1),fBinsEta->GetAt(fBinsEta->GetSize()-1),fBinsZv->GetAt(fBinsZv->GetSize()-1),fBinsMultCent->GetAt(fBinsMultCent->GetSize()-1)};

    /// Standard event histogram zV:multcent
    Int_t nBinsEvent[2]={fBinsZv->GetSize()-1,fBinsMultCent->GetSize()-1};
    Double_t minEvent[2]={fBinsZv->GetAt(0),fBinsMultCent->GetAt(0)};
    Double_t maxEvent[2]={fBinsZv->GetAt(fBinsZv->GetSize()-1),fBinsMultCent->GetAt(fBinsMultCent->GetSize()-1)};

    /// test pT vs pT also eta and vz
    Int_t nBinsMCPtPt[4]={fBinsPt->GetSize()-1,fBinsPt->GetSize()-1,fBinsEta->GetSize()-1,fBinsZv->GetSize()-1};
    Double_t minBinsMCPtPt[4]={fBinsPt->GetAt(0),fBinsPt->GetAt(0),fBinsEta->GetAt(0),fBinsZv->GetAt(0)};
    Double_t maxBinsMCPtPt[4]={fBinsPt->GetAt(fBinsPt->GetSize()-1),fBinsPt->GetAt(fBinsPt->GetSize()-1),fBinsEta->GetAt(fBinsZv->GetSize()-1),fBinsZv->GetAt(fBinsZv->GetSize()-1) };

    /// relative pT resolution from covariance matrix (global tracks) as a function of pt and centrality
    Int_t nBinsRelPtReso[3]  = {fBinsPtReso->GetSize()-1,fBinsPt->GetSize()-1, fBinsMultCent->GetSize()-1};
    Double_t minRelPtReso[3] = {fBinsPtReso->GetAt(0),fBinsPt->GetAt(0), fBinsMultCent->GetAt(0)};
    Double_t maxRelPtReso[3] = {fBinsPtReso->GetAt(fBinsPtReso->GetSize()-1), fBinsPt->GetAt(fBinsPt->GetSize()-1), fBinsMultCent->GetAt(fBinsMultCent->GetSize()-1)};

    /// binning of correlation of centrality estimators (with pT)
    Int_t iBinCentCorrpt[5] = {100, 100, 100, 100, fBinsPt->GetSize()-1};
    Double_t iMinCentCorrpt[5] = {0, 0, 0, 0, fBinsPt->GetAt(0)};
    Double_t iMaxCentCorrpt[5] = {100, 100, 100, 100, fBinsPt->GetAt(fBinsPt->GetSize()-1)};
    /// binning of correlation of centrality estimators (without pT)
    Int_t iBinCentCorr[4] = {100, 100, 100, 100};
    Double_t iMinCentCorr[4] = {0, 0, 0, 0};
    Double_t iMaxCentCorr[4] = {100, 100, 100, 100};

    fHistTrack = new THnF("fHistTrack", "Histogram for Tracks",4,nBinsTrack,minTrack,maxTrack);
    fHistTrack -> SetBinEdges(0,fBinsPt->GetArray());
    fHistTrack -> SetBinEdges(1,fBinsEta->GetArray());
    fHistTrack -> SetBinEdges(2,fBinsZv->GetArray());
    fHistTrack -> SetBinEdges(3,fBinsMultCent->GetArray());
    fHistTrack->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
    fHistTrack->GetAxis(1)->SetTitle("#eta");
    fHistTrack->GetAxis(2)->SetTitle("Zv (cm)");
    fHistTrack->GetAxis(3)->SetTitle("multiplicity (multCuts)");
    fHistTrack -> Sumw2();

    fHistEvent = new THnF("fHistEvent", "Histogram for Events",2,nBinsEvent,minEvent,maxEvent);
    fHistEvent -> SetBinEdges(0,fBinsZv->GetArray());
    fHistEvent -> SetBinEdges(1,fBinsMultCent->GetArray());
    fHistEvent->GetAxis(0)->SetTitle("Zv (cm)");
    fHistEvent->GetAxis(1)->SetTitle("multiplicity (multCuts)");
    fHistEvent -> Sumw2();

    fEventCount = new THnF("fEventCount","trig vs trig+vertex",4,binsEventCount,minEventCount,maxEventCount);
    fEventCount->GetAxis(0)->SetTitle("trig");
    fEventCount->GetAxis(1)->SetTitle("trig+vert");
    fEventCount->GetAxis(2)->SetTitle("selected");
    fEventCount->GetAxis(3)->SetTitle("All");
    fEventCount->Sumw2();   

    fHistCentCorrelpt = new THnSparseF("fHistCentCorrelpt", "Centrality Correlation pt", 5, iBinCentCorrpt, iMinCentCorrpt, iMaxCentCorrpt);
    fHistCentCorrelpt->GetAxis(0)->SetTitle("Centrality V0M");
    fHistCentCorrelpt->GetAxis(1)->SetTitle("Centrality SPD Tracklets");
    fHistCentCorrelpt->GetAxis(2)->SetTitle("Centrality CL0");
    fHistCentCorrelpt->GetAxis(3)->SetTitle("Centrality CL1");
    fHistCentCorrelpt->GetAxis(4)->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistCentCorrelpt->SetBinEdges(4,fBinsPt->GetArray());
    fHistCentCorrelpt->Sumw2();

    fHistCentCorrelBeforeCuts = new THnSparseF("fHistCentCorrelBeforeCuts", "Centrality Correlation BC", 4, iBinCentCorr, iMinCentCorr, iMaxCentCorr);
    fHistCentCorrelBeforeCuts->GetAxis(0)->SetTitle("Centrality V0M");
    fHistCentCorrelBeforeCuts->GetAxis(1)->SetTitle("Centrality SPD Tracklets");
    fHistCentCorrelBeforeCuts->GetAxis(2)->SetTitle("Centrality CL0");
    fHistCentCorrelBeforeCuts->GetAxis(3)->SetTitle("Centrality CL1");
    fHistCentCorrelBeforeCuts->Sumw2();
    
    fHistCentCorrelAfterCuts = new THnSparseF("fHistCentCorrelAfterCuts", "Centrality Correlation AC", 4, iBinCentCorr, iMinCentCorr, iMaxCentCorr);
    fHistCentCorrelAfterCuts->GetAxis(0)->SetTitle("Centrality V0M");
    fHistCentCorrelAfterCuts->GetAxis(1)->SetTitle("Centrality SPD Tracklets");
    fHistCentCorrelAfterCuts->GetAxis(2)->SetTitle("Centrality CL0");
    fHistCentCorrelAfterCuts->GetAxis(3)->SetTitle("Centrality CL1");
    fHistCentCorrelAfterCuts->Sumw2();
    
    fHistRelPtResoFromCov = new THnF("fHistRelPtResoFromCov", "Relative pT resolution from covariance matrix", 3, nBinsRelPtReso, minRelPtReso, maxRelPtReso);
    fHistRelPtResoFromCov->SetBinEdges(0,fBinsPtReso->GetArray());
    fHistRelPtResoFromCov->SetBinEdges(1,fBinsPt->GetArray());
    fHistRelPtResoFromCov->SetBinEdges(2,fBinsMultCent->GetArray());
    fHistRelPtResoFromCov->GetAxis(0)->SetTitle("#sigma(#it{p}_{T}) / #it{p}_{T}");
    fHistRelPtResoFromCov->GetAxis(1)->SetTitle("#it{p}_{T} (GeV/#it{c})");
    fHistRelPtResoFromCov->GetAxis(2)->SetTitle("Centrality (%)");
    fHistRelPtResoFromCov->Sumw2();

    ///Secondary scaling
    fDCAyEtaPt = new TH3D("fDCAyEtaPt","DCAy:eta:pt;DCAy (cm); eta; pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,20,-1,1,ptNbins,ptMin,ptMax);

    if(fIsMC){
        fHistMCGenPrimTrack = new THnF("fHistMCGenPrimTrack", "Histogram for generated MC Tracks",4,nBinsTrack,minTrack,maxTrack);
        fHistMCGenPrimTrack -> SetBinEdges(0,fBinsPt->GetArray());
        fHistMCGenPrimTrack -> SetBinEdges(1,fBinsEta->GetArray());
        fHistMCGenPrimTrack -> SetBinEdges(2,fBinsZv->GetArray());
        fHistMCGenPrimTrack -> SetBinEdges(3,fBinsMultCent->GetArray());
        fHistMCGenPrimTrack->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
        fHistMCGenPrimTrack->GetAxis(1)->SetTitle("#eta");
        fHistMCGenPrimTrack->GetAxis(2)->SetTitle("Zv (cm)");
        fHistMCGenPrimTrack->GetAxis(3)->SetTitle("true multiplicity (MC)");
        fHistMCGenPrimTrack -> Sumw2();

        fHistMCRecTrack = new THnF("fHistMCRecTrack", "Histogram for reconstructed MC Tracks",4,nBinsTrack,minTrack,maxTrack);
        fHistMCRecTrack -> SetBinEdges(0,fBinsPt->GetArray());
        fHistMCRecTrack -> SetBinEdges(1,fBinsEta->GetArray());
        fHistMCRecTrack -> SetBinEdges(2,fBinsZv->GetArray());
        fHistMCRecTrack -> SetBinEdges(3,fBinsMultCent->GetArray());
        fHistMCRecTrack->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
        fHistMCRecTrack->GetAxis(1)->SetTitle("#eta");
        fHistMCRecTrack->GetAxis(2)->SetTitle("Zv (cm)");
        fHistMCRecTrack->GetAxis(3)->SetTitle("true multiplicity (MC)");
        fHistMCRecTrack -> Sumw2();

        fHistMCRecPrimTrack = new THnF("fHistMCRecPrimTrack", "Histogram for reconstructed primary MC Tracks",4,nBinsTrack,minTrack,maxTrack);
        fHistMCRecPrimTrack -> SetBinEdges(0,fBinsPt->GetArray());
        fHistMCRecPrimTrack -> SetBinEdges(1,fBinsEta->GetArray());
        fHistMCRecPrimTrack -> SetBinEdges(2,fBinsZv->GetArray());
        fHistMCRecPrimTrack -> SetBinEdges(3,fBinsMultCent->GetArray());
        fHistMCRecPrimTrack->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
        fHistMCRecPrimTrack->GetAxis(1)->SetTitle("#eta");
        fHistMCRecPrimTrack->GetAxis(2)->SetTitle("Zv (cm)");
        fHistMCRecPrimTrack->GetAxis(3)->SetTitle("true multiplicity (MC)");
        fHistMCRecPrimTrack -> Sumw2();

        fHistMCRecSecTrack = new THnF("fHistMCRecSecTrack", "Histogram for reconstructed secondary MC Tracks",4,nBinsTrack,minTrack,maxTrack);
        fHistMCRecSecTrack -> SetBinEdges(0,fBinsPt->GetArray());
        fHistMCRecSecTrack -> SetBinEdges(1,fBinsEta->GetArray());
        fHistMCRecSecTrack -> SetBinEdges(2,fBinsZv->GetArray());
        fHistMCRecSecTrack -> SetBinEdges(3,fBinsMultCent->GetArray());
        fHistMCRecSecTrack->GetAxis(0)->SetTitle("p_{T} (GeV/c)");
        fHistMCRecSecTrack->GetAxis(1)->SetTitle("#eta");
        fHistMCRecSecTrack->GetAxis(2)->SetTitle("Zv (cm)");
        fHistMCRecSecTrack->GetAxis(3)->SetTitle("true multiplicity (MC)");
        fHistMCRecSecTrack -> Sumw2();

        fHistMCTrigEvent = new THnF("fHistMCTrigEvent", "Histogram for triggered MC Events",2,nBinsEvent,minEvent,maxEvent);
        fHistMCTrigEvent -> SetBinEdges(0,fBinsZv->GetArray());
        fHistMCTrigEvent -> SetBinEdges(1,fBinsMultCent->GetArray());
        fHistMCTrigEvent->GetAxis(0)->SetTitle("Zv (cm)");
        fHistMCTrigEvent->GetAxis(1)->SetTitle("true multiplicity (MC)");
        fHistMCTrigEvent -> Sumw2();

        fHistMCGenINEL0Event = new THnF("fHistMCGenINEL0Event", "Histogram for generated INEL>0 MC Events",2,nBinsEvent,minEvent,maxEvent);
        fHistMCGenINEL0Event -> SetBinEdges(0,fBinsZv->GetArray());
        fHistMCGenINEL0Event -> SetBinEdges(1,fBinsMultCent->GetArray());
        fHistMCGenINEL0Event->GetAxis(0)->SetTitle("Zv (cm)");
        fHistMCGenINEL0Event->GetAxis(1)->SetTitle("true multiplicity (MC)");
        fHistMCGenINEL0Event -> Sumw2();

        fHistMCTrigINEL0Event = new THnF("fHistMCTrigINEL0Event","Histogram for triggered INEL>0 MC Events",2,nBinsEvent,minEvent,maxEvent);
        fHistMCTrigINEL0Event->SetBinEdges(0,fBinsZv->GetArray());
        fHistMCTrigINEL0Event->SetBinEdges(1,fBinsMultCent->GetArray());
        fHistMCTrigINEL0Event->GetAxis(0)->SetTitle("Zv (cm)");
        fHistMCTrigINEL0Event->GetAxis(1)->SetTitle("true multiplicity (MC)");
        fHistMCTrigINEL0Event->Sumw2();

        fHistMCRecINEL0Event = new THnF("fHistMCRecINEL0Event","Histogram for reconstructed INEL>0 MC Events",2,nBinsEvent,minEvent,maxEvent);
        fHistMCRecINEL0Event->SetBinEdges(0,fBinsZv->GetArray());
        fHistMCRecINEL0Event->SetBinEdges(1,fBinsMultCent->GetArray());
        fHistMCRecINEL0Event->GetAxis(0)->SetTitle("mcZv (cm)");
        fHistMCRecINEL0Event->GetAxis(1)->SetTitle("true multiplicity (MC)");
        fHistMCRecINEL0Event->Sumw2();

        ///Secondary scaling MC
        fDCAyEtaPtMCPrim = new TH3D("fDCAyEtaPtMCPrim","DCAy:eta:pt primary;DCAy (cm); eta; pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,20,-1,1,ptNbins,ptMin,ptMax);
        fDCAyEtaPtMCSecDecays = new TH3D("fDCAyEtaPtMCSecDecays","DCAy:eta:pt secondary decays;DCAy (cm); eta; pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,20,-1,1,ptNbins,ptMin,ptMax);
        fDCAyEtaPtMCSecMaterial = new TH3D("fDCAyEtaPtMCSecMaterial","DCAy:eta:pt secondary material;DCAy (cm); eta; pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,20,-1,1,ptNbins,ptMin,ptMax);
        fDCAyEtaPtMCSecDecaysK0s = new TH3D("fDCAyEtaPtMCSecDecaysK0s","DCAy:eta:pt secondary decays from K0s;DCAy (cm); eta; pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,20,-1,1,ptNbins,ptMin,ptMax);
        fDCAyEtaPtMCSecDecaysLambda = new TH3D("fDCAyEtaPtMCSecDecaysLambda","DCAy:eta:pt secondary decays from Lambda;DCAy (cm); eta; pt (GeV/c)", DCAybins, DCAyMin, DCAyMax,20,-1,1,ptNbins,ptMin,ptMax);

        //particle composition correction
        // fMCSpectraWeights = new AliMCSpectraWeights(fstCollisionSystem.Data(), "fMCSpectraWeights");
        // if(fstMCTrainOutput!="") fMCSpectraWeights->SetMCSpectraFile(fstMCTrainOutput.Data());
        // fMCSpectraWeights->Init();
        //
        // fHistMCPartCorr = fMCSpectraWeights->GetHistMCGenPrimTrackParticles();
    }

    fOutputList->Add(fHistTrack);
    fOutputList->Add(fHistCentCorrelpt);
    fOutputList->Add(fHistCentCorrelBeforeCuts);
    fOutputList->Add(fHistCentCorrelAfterCuts);
    fOutputList->Add(fHistEvent);
    fOutputList->Add(fEventCount);
    fOutputList->Add(fHistRelPtResoFromCov);
    fOutputList->Add(fDCAyEtaPt);
    if(fIsMC){
        fOutputList->Add(fHistMCGenPrimTrack);
        fOutputList->Add(fHistMCRecTrack);
        fOutputList->Add(fHistMCRecPrimTrack);
        fOutputList->Add(fHistMCRecSecTrack);
        fOutputList->Add(fHistMCTrigEvent);
        fOutputList->Add(fHistMCGenINEL0Event);
        fOutputList->Add(fHistMCTrigINEL0Event);
        fOutputList->Add(fHistMCRecINEL0Event);
        fOutputList->Add(fDCAyEtaPtMCPrim);
        fOutputList->Add(fDCAyEtaPtMCSecDecays);
        fOutputList->Add(fDCAyEtaPtMCSecMaterial);
        fOutputList->Add(fDCAyEtaPtMCSecDecaysK0s);
        fOutputList->Add(fDCAyEtaPtMCSecDecaysLambda);
        // fOutputList->Add(fMCSpectraWeights);
        // fOutputList->Add(fHistMCPartCorr);
    }
    PostData(1, fOutputList);

    /// Create and initialize Analysis Objects here instead of in UserExec() to save resources
    //InitdNdPtEventCuts(); (does nothing atm so I just leave it out for the moment)
    if(fIsESD) InitESDTrackCuts();
    else  { //InitAODTrackCuts();
    }
    if((fIs2013pA || fIs2015data) && !fUtils){fUtils = new AliAnalysisUtils();}
}

/// Destructor
AlidNdPtUnifiedAnalysisTask::~AlidNdPtUnifiedAnalysisTask(){
    /*if(fUtils){delete fUtils; fUtils=0;}
      if(fESDtrackCuts){delete fESDtrackCuts; fESDtrackCuts=0;}
      if(fAODtrackCuts){delete fAODtrackCuts; fAODtrackCuts=0;} */
    //if(fEventCuts){delete fEventCuts; fEventCuts=0;}
}

///________________________________________________________________________
void AlidNdPtUnifiedAnalysisTask::UserExec(Option_t *){ // Main loop (called for each event)
    /// ====================== Initialize variables ===============================
    AliInputEventHandler* inputHandler = (AliInputEventHandler*) AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler();
    if (!inputHandler){ Printf("ERROR: Could not receive inputHandler"); return; }
    AliPhysicsSelection *physicsSelection = static_cast<AliPhysicsSelection*> (inputHandler->GetEventSelection());
    if(!physicsSelection) {Printf("ERROR: Could not receive physicsSelection"); return;}
    fEvent = dynamic_cast<AliVEvent*>(InputEvent());
    if (!fEvent) {printf("ERROR: fEvent not available\n"); return;}
    Bool_t isEventINEL0 = kFALSE;
    if(fIsMC){
        fMCEvent = dynamic_cast<AliMCEvent*>(MCEvent());
        if (!fMCEvent) {printf("ERROR: fMCEvent not available\n"); return;}
        fMCStack = fMCEvent->Stack();
        if (!fMCStack) {printf("ERROR: fMCStack not available\n"); return;}
        // Check if MC event is in inel0 class (1 charged particle in abs(eta)<1.0, pt>0)
        isEventINEL0 = IsMCEventINEL0(fMCEvent,0,1.0);
    }

    Bool_t isEventTriggered = inputHandler->IsEventSelected() & GetTriggerMask();
    Double_t multEvent = GetEventMultCent(fEvent);
    
    //TODO make here all cent estimators
    AliMultSelection *MultSelection = (AliMultSelection*) fEvent->FindListObject("MultSelection");
    if(!MultSelection) return;
    

    Double_t zVertEvent = fEvent->GetPrimaryVertex()->GetZ();
    Double_t eventValues[2] = {zVertEvent, multEvent};

    AliVVZERO * vZeroHandler = fEvent->GetVZEROData();
    if (!vZeroHandler) {printf("ERROR: vZeroHandler not available\n"); return;}


    Bool_t isTrigAndVertex = kFALSE;
    if (fIs2013pA){	isTrigAndVertex = isEventTriggered && IsEventAccepted2013pA(fEvent) && IsEventAcceptedQuality(fEvent);	}
    if (fIs2015data){	isTrigAndVertex = isEventTriggered && IsEventAccepted2015data(fEvent) && IsEventAcceptedQuality(fEvent);	}


    Double_t vEventCount[4] = { static_cast<Double_t>((isEventTriggered && kTRUE)) , static_cast<Double_t>(isTrigAndVertex),  static_cast<Double_t>(isTrigAndVertex && (TMath::Abs(zVertEvent) < 10.)), kTRUE};
    fEventCount->Fill(vEventCount);

    Double_t multGenPart = 0;  		/// N_ch
    Double_t multRecPart = 0;		/// N_acc (of tracks that can be assigned to a real particle)
    Double_t multRecCorrPart = 0;		/// N_acc corrected with trk efficiency

    Double_t V0MCent = MultSelection->GetMultiplicityPercentile("V0M");
    Double_t CL0Cent = MultSelection->GetMultiplicityPercentile("CL0");
    Double_t CL1Cent = MultSelection->GetMultiplicityPercentile("CL1");
    Double_t SPDTCent = MultSelection->GetMultiplicityPercentile("SPDTracklets");  
    Double_t eventCorrelValues[4] = {V0MCent, SPDTCent, CL0Cent, CL1Cent};
    /// ==================== Fill Histogramms ======================================

    /// -------------------- Generated Events --------------------------------------

    if(fIsMC && isEventINEL0)fHistMCGenINEL0Event->Fill(eventValues);

    /// \li Event Trigger Cut
    if(!isEventTriggered) return;
    fHistCentCorrelBeforeCuts->Fill(eventCorrelValues);

    /// -------------------- Triggered Events ---------------------------------------

    if(fIsMC && isEventINEL0) fHistMCTrigINEL0Event->Fill(eventValues);

    /// \li Event Acceptance Cuts
    if (!fAliEventCuts.AcceptEvent(fEvent)) return;
    
    /// ------------------ Reconstructed Events --------------------------------------

    if(fIsMC && isEventINEL0) fHistMCRecINEL0Event->Fill(eventValues);
    fHistEvent->Fill(eventValues);
    
    fHistCentCorrelAfterCuts->Fill(eventCorrelValues);

    ///--------------- Loop over measured Tracks ---------------------------------

    AliVTrack *track = NULL;
    for (Int_t iTrack = 0; iTrack < fEvent->GetNumberOfTracks(); iTrack++){
        track = fEvent->GetVTrack(iTrack);
        if (!track){ printf("ERROR: Could not receive track %d\n", iTrack); continue; }

        /// \li Track Acceptance Cuts
        if(!IsTrackAcceptedKinematics(track)) continue;
        if(!IsTrackAcceptedQuality(track)) continue;

        if(fUseTOFBunchCrossing && fIs2015data){
            if(TMath::Abs(track->GetTOFsignalDz())>10) continue;
            if((track->GetTOFsignal())<12000) continue;
            if((track->GetTOFsignal())>25000) continue;
        }

        Double_t dWeight = 1.0;
        // if(fIsMC)
        // {
        //     if(fMCSpectraWeights->GetTaskStatus()>=AliMCSpectraWeights::TaskState::kMCSpectraObtained){
        //       TParticle* part = fMCStack->Particle(TMath::Abs(track->GetLabel()));
        //       if(part) dWeight = fMCSpectraWeights->GetMCSpectraWeight(part ,multEvent);
        //     }
        // }

        /// \li Fill Track Histograms
        Double_t dpT = track->Pt();
        Double_t trackValues[4] = {dpT, track->Eta(), zVertEvent, multEvent};
        fHistTrack->Fill(trackValues, dWeight);
        Double_t trackCorrelValuespt[5] = {V0MCent, SPDTCent, CL0Cent, CL1Cent, dpT};
        fHistCentCorrelpt->Fill(trackCorrelValuespt);

        /// Histo needed to secondary scaling
        Float_t b[2], bCov[3];
        track->GetImpactParameters(b,bCov);
        fDCAyEtaPt->Fill(b[0],track->Eta(),dpT);

        // pT resolution for covariance matrix
        ///TODO this only works for ESD tracks!!
        Double_t ptResoValues[3] = {1./TMath::Abs(dynamic_cast<AliESDtrack*>(track)->GetSigned1Pt())*TMath::Sqrt(dynamic_cast<AliESDtrack*>(track)->GetSigma1Pt2()), track->Pt(), multEvent};
        fHistRelPtResoFromCov->Fill(ptResoValues);

        /// \li Find original particle in MC-Stack
        if(fIsMC){
            Int_t mcLabel = TMath::Abs(track->GetLabel());
            TParticle *mcParticle = fMCStack->Particle(mcLabel);
            if(!mcParticle) {printf("ERROR: mcParticle not available\n"); continue;}
            if(!IsTrackAcceptedKinematics(mcParticle)) continue;
            // -------------- Secondary scaling ----------------
            Bool_t isPrim = fMCEvent->IsPhysicalPrimary(mcLabel);
            Bool_t isWeakDecay = fMCEvent->IsSecondaryFromWeakDecay(mcLabel);
            Bool_t isFromMaterial = fMCEvent->IsSecondaryFromMaterial(mcLabel);
            Bool_t isFromK0s = kFALSE;
            Bool_t isFromLambda = kFALSE;

            // check whether has stange mother
            //
            Int_t motherPdg = -1;
            TParticle* mother = 0;

            Int_t motherLabel = mcParticle->GetFirstMother();
            if(motherLabel>=0) mother = fMCStack->Particle(motherLabel);
            if(mother) motherPdg = TMath::Abs(mother->GetPdgCode()); // take abs for visualisation only
            Int_t mech = mcParticle->GetUniqueID(); // production mechanism

            if(isWeakDecay && motherPdg == 310) isFromK0s = kTRUE;
            if(isWeakDecay && motherPdg == 3122) isFromLambda = kTRUE;

            if(isPrim) fDCAyEtaPtMCPrim->Fill(b[0],mcParticle->Eta(),mcParticle->Pt());
            else
            {
                if(isWeakDecay)
                {
                    fDCAyEtaPtMCSecDecays->Fill(b[0],mcParticle->Eta(),mcParticle->Pt());
                    if(isFromK0s) fDCAyEtaPtMCSecDecaysK0s->Fill(b[0],mcParticle->Eta(),mcParticle->Pt());
                    if(isFromLambda) fDCAyEtaPtMCSecDecaysLambda->Fill(b[0],mcParticle->Eta(),mcParticle->Pt());
                }
                if(isFromMaterial) fDCAyEtaPtMCSecMaterial->Fill(b[0],mcParticle->Eta(),mcParticle->Pt());
            }
            // ----------- secondary scaling end ----------
            Double_t mcRecTrackValue[4] = {mcParticle->Pt(), mcParticle->Eta(), zVertEvent, multEvent};
            fHistMCRecTrack->Fill(mcRecTrackValue, dWeight);
            if(IsChargedPrimary(mcLabel))
            {
                Double_t mcPrimTrackValue[4] = {mcParticle->Pt(), mcParticle->Eta(), zVertEvent, multEvent};
                fHistMCRecPrimTrack->Fill(mcPrimTrackValue, dWeight);
            }else{//works because tracks are always charged
                Double_t mcSecTrackValue[4] = {mcParticle->Pt(), mcParticle->Eta(), zVertEvent, multEvent};
                fHistMCRecSecTrack->Fill(mcSecTrackValue);
            }
        }
    }// end of Track-loop

    ///------------------- Loop over Generated Tracks (True MC)------------------------------
    if (fIsMC){
        // if(fMCSpectraWeights->GetTaskStatus()<AliMCSpectraWeights::TaskState::kMCSpectraObtained)
        // fMCSpectraWeights->FillMCSpectra(fMCEvent, multEvent);
        for (Int_t iParticle = 0; iParticle < fMCStack->GetNtrack(); iParticle++){
            TParticle *mcGenParticle = fMCStack->Particle(iParticle);
            if(!mcGenParticle) {printf("ERROR: mcGenParticle  not available\n"); continue;}

            /// \li Acceptance cuts for generated particles
            if(!IsTrackAcceptedKinematics(mcGenParticle, kTRUE)) continue;
            if(IsChargedPrimary(iParticle)){
                Double_t dWeight = 1.0;
                // if(fMCSpectraWeights->GetTaskStatus()>=AliMCSpectraWeights::TaskState::kMCSpectraObtained){
                // dWeight = fMCSpectraWeights->GetMCSpectraWeight(mcGenParticle ,multEvent);
                // }
                Double_t mcGenPrimTrackValue[4] = {mcGenParticle->Pt(), mcGenParticle->Eta(), zVertEvent, multEvent};
                fHistMCGenPrimTrack->Fill(mcGenPrimTrackValue, dWeight);
            }
        }
    }

    PostData(1, fOutputList);
}

//________________________________________________________________________
void AlidNdPtUnifiedAnalysisTask::Terminate(Option_t *)
{

}

Bool_t AlidNdPtUnifiedAnalysisTask::IsChargedPrimary(Int_t stackIndex){
    if (fMCStack->IsPhysicalPrimary(stackIndex) && (TMath::Abs(fMCStack->Particle(stackIndex)->GetPDG()->Charge()) > 0.01)){
        return kTRUE;
    }
    return kFALSE;
}

/// Track Acceptance cuts, for tracks.
///
/// \param AliVTrack Input track
///
/// \return Is track accepted: kTRUE, else kFALSE
Bool_t AlidNdPtUnifiedAnalysisTask::IsTrackAcceptedKinematics(AliVTrack *track)
{
    if(!track) return kFALSE;

    Float_t eta = track->Eta();	///TODO why float?
    Float_t pt = track->Pt();

    if(eta < fMinEta) return kFALSE;
    if(eta > fMaxEta) return kFALSE;
    if(pt < fMinPt) return kFALSE;
    if(pt > fMaxPt) return kFALSE;

    return kTRUE;
}

/// Track Acceptance cuts, for MC particles.
///
/// \param TParticle Input particle
///
/// \return Is particle accepted: kTRUE, else kFALSE
Bool_t AlidNdPtUnifiedAnalysisTask::IsTrackAcceptedKinematics(TParticle *mcTrack, Bool_t useLowerPtCut)
{
    if(!mcTrack) return kFALSE;

    Float_t eta = mcTrack->Eta();
    Float_t pt = mcTrack->Pt();

    if(eta < fMinEta) return kFALSE;
    if(eta > fMaxEta) return kFALSE;
    if((pt < fMinPt) && useLowerPtCut) return kFALSE;
    if(pt > fMaxPt) return kFALSE;
    return kTRUE;
}

/// Track Quality cuts.
///
/// \param AliVTrack Input track
///
/// \return Is track accepted: kTRUE, else kFALSE
Bool_t AlidNdPtUnifiedAnalysisTask::IsTrackAcceptedQuality(AliVTrack *track){
    if(fIsESD){
        AliESDtrack *esdTrack = dynamic_cast<AliESDtrack*> (track);
        if(!fESDtrackCuts->AcceptTrack(esdTrack)) return kFALSE;
    }
    else{
        //AliAODTrack *aodTrack = dynamic_cast<AliAODTrack*> (track);
        //if(!fAODtrackCuts->AcceptTrack(aodTrack)) return kFALSE;
        //return aodTrack->TestFilterBit(AliAODTrack::AODTrkFilterBits_t::kTrkGlobal);
    }
    return kTRUE;
}

/// selection and pileup rejection for 2013 p-A
Bool_t AlidNdPtUnifiedAnalysisTask::IsEventAccepted2013pA(AliVEvent *event)
{
    if (fUtils->IsFirstEventInChunk(event)) { return kFALSE;  }
    if (!fUtils->IsVertexSelected2013pA(event)) { return kFALSE;  }
    if (fUtils->IsPileUpEvent(event)) { return kFALSE;  }
    return kTRUE;
}

/// Event cuts and pileup rejection for 2015 data: pp and Pb-Pb
Bool_t AlidNdPtUnifiedAnalysisTask::IsEventAccepted2015data(AliVEvent *event)
{
    ///ESD
    if(fIsESD){
        AliESDEvent* ESDevent = dynamic_cast<AliESDEvent*>(event);
        // if (fUtils->IsFirstEventInChunk(event)) { return kFALSE;  }
        // if (!fUtils->IsVertexSelected2013pA(event)) { return kFALSE;  }
        if (ESDevent->IsIncompleteDAQ()) { return kFALSE; }
        if (fUtils->IsSPDClusterVsTrackletBG(event)) { return kFALSE; }
        if (ESDevent->IsPileupFromSPD(5,0.8)) {return kFALSE; }
    }
    else
    {
        /*
        ///AOD
        AliAODEvent* AODevent = dynamic_cast<AliAODEvent*>(event);
        if (AODevent->IsIncompleteDAQ()) { return kFALSE; }
        if (fUtils->IsSPDClusterVsTrackletBG(event)) { return kFALSE; }
        if (AODevent->IsPileupFromSPD(5,0.8)) {return kFALSE; }
        */
    }
    return kTRUE;
}

/// Event Acceptance cuts.
///
/// \param AliVEvent Input event
///
/// \return Is event accepted: kTRUE, else kFALSE
Bool_t AlidNdPtUnifiedAnalysisTask::IsEventAcceptedGeometrics(AliVEvent *event)
{
    // when fEventCuts is used (i.e. something happens in InitdNdPtEventCuts), we could use  something like if(!fEventCuts->AcceptEvent)
    // maybe even without this* function...
    if(TMath::Abs(event->GetPrimaryVertex()->GetZ())>fZvtx) return kFALSE;
    return kTRUE;
}

/// Event Quality cuts.
///
/// \param AliVEvent Input event
///
/// \return Is event accepted: kTRUE, else kFALSE
Bool_t AlidNdPtUnifiedAnalysisTask::IsEventAcceptedQuality(AliVEvent *event)
{
    //   AliVVertex *vertex = event->GetPrimaryVertexTracks();
    if(!event) return kFALSE;
    if(!IsVertexOK(event)) return kFALSE;
    return kTRUE;
}

/// Function to set centrality estimator
///
/// \param iCentEst from enum WhichCentralityEstimator
void AlidNdPtUnifiedAnalysisTask::SetCentralityEstimator(WhichCentralityEstimator iCentEst)
{
    switch(iCentEst)
    {
        case AlidNdPtUnifiedAnalysisTask::WhichCentralityEstimator::kV0A:
            fstCentEst = "V0A"; break;
        case AlidNdPtUnifiedAnalysisTask::WhichCentralityEstimator::kV0C:
            fstCentEst = "V0C"; break;
        case AlidNdPtUnifiedAnalysisTask::WhichCentralityEstimator::kV0M:
            fstCentEst = "V0M"; break;
        case AlidNdPtUnifiedAnalysisTask::WhichCentralityEstimator::kCL0:
            fstCentEst = "CL0"; break;
        case AlidNdPtUnifiedAnalysisTask::WhichCentralityEstimator::kCL1:
            fstCentEst = "CL1"; break;
        case AlidNdPtUnifiedAnalysisTask::WhichCentralityEstimator::kZNA:
            fstCentEst = "ZNA"; break;
        case AlidNdPtUnifiedAnalysisTask::WhichCentralityEstimator::kZNC:
            fstCentEst = "ZNC"; break;
        case AlidNdPtUnifiedAnalysisTask::WhichCentralityEstimator::kSPDc:
            fstCentEst = "SPDClusters"; break;
        case AlidNdPtUnifiedAnalysisTask::WhichCentralityEstimator::kSPDt:
            fstCentEst = "SPDTracklets"; break;

        default:
            fstCentEst = "V0M";
            printf ("No / invalid centrality estimator has been chosen. Setting to default centrality estimator %s .", fstCentEst.Data()); break;
    }
    // printf ("Centrality estimator %s has been set.", fstCentEst.Data());
}

/// Function to either fill Multiplicity or centrality
///
/// \param AliVEvent event to be analised
///
/// \return Double_t with centrality or multiplicity

Double_t AlidNdPtUnifiedAnalysisTask::GetEventMultCent(AliVEvent *event)
{
    if(fUseMultiplicity)
    {
        AliVMultiplicity* multiplicity = event->GetMultiplicity();
        if(!multiplicity) {printf("ERROR: multiplicity not available\n"); return 999;}
        Int_t mult = multiplicity->GetNumberOfTracklets();
        return mult;
    }
    else
    {
        Float_t centralityF = -1;
        AliMultSelection *MultSelection = (AliMultSelection*) fEvent->FindListObject("MultSelection");
        if ( MultSelection ){
            centralityF = MultSelection->GetMultiplicityPercentile(fstCentEst.Data());
            if(centralityF>100) return 999;
            return centralityF;
        }else{
            AliInfo("Didn't find MultSelection!");
            return 999;
        }
    }
}

/// Function to initialize the ESD track cuts
void AlidNdPtUnifiedAnalysisTask::InitESDTrackCuts(){

    fESDtrackCuts = new AliESDtrackCuts("AliESDtrackCuts");
    if(!fESDtrackCuts) {printf("ERROR: fESDtrackCuts not available\n"); return;}


    fESDtrackCuts->SetRequireTPCRefit(fTPCRefit);
    fESDtrackCuts->SetRequireITSRefit(fITSRefit);
    fESDtrackCuts->SetAcceptKinkDaughters(fAcceptKinks);
    if(fMinNCrossedRowsTPC > 0) fESDtrackCuts->SetMinNCrossedRowsTPC(fMinNCrossedRowsTPC);
    if(fMinRatioCrossedRowsOverFindableClustersTPC > 0) fESDtrackCuts->SetMinRatioCrossedRowsOverFindableClustersTPC(fMinRatioCrossedRowsOverFindableClustersTPC);
    if(fMaxFractionSharedClustersTPC > 0) fESDtrackCuts->SetMaxFractionSharedTPCClusters(fMaxFractionSharedClustersTPC);
    if(fMaxChi2PerTPCCluster > 0) fESDtrackCuts->SetMaxChi2PerClusterTPC(fMaxChi2PerTPCCluster);
    if(fRequiresClusterITS) fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
    if(!fRequiresClusterITS)fESDtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
    if(fMaxChi2PerITSCluster > 0) fESDtrackCuts->SetMaxChi2PerClusterITS(fMaxChi2PerITSCluster);
    if(fDCAToVertex2D > 0) fESDtrackCuts->SetDCAToVertex2D(fDCAToVertex2D);
    if(fSigmaToVertex > 0) fESDtrackCuts->SetRequireSigmaToVertex(fSigmaToVertex);
    if(fMaxDCAzITSTPC > 0) fESDtrackCuts->SetMaxDCAToVertexZ(fMaxDCAzITSTPC);
    if(fDCAToVertexXY > 0) fESDtrackCuts->SetMaxDCAToVertexXY(fDCAToVertexXY);
    if(fDCAToVertexXYPtDep)          fESDtrackCuts->SetMaxDCAToVertexXYPtDep(fDCAToVertexXYPtDep);
    if(fMaxChi2TPCConstrained > 0) fESDtrackCuts->SetMaxChi2TPCConstrainedGlobal(fMaxChi2TPCConstrained);
    if(fMinActiveLength> 0) fESDtrackCuts->SetMinLengthActiveVolumeTPC(fMinActiveLength);
    if(fUseGeomCut) fESDtrackCuts->SetCutGeoNcrNcl(fDeadZoneWidth,fCutGeoNcrNclLenght,fCutGeoNcrNclGeom1Pt,fCutGeoNcrNclFractionNcr,fCutGeoNcrNclFractionNcl);

}

///Function to check vertex quality
Bool_t AlidNdPtUnifiedAnalysisTask::IsVertexOK(AliVEvent *event){
    if(fIsESD){
        Float_t requiredZResolution = 1000;
        AliESDEvent* ESDevent = dynamic_cast<AliESDEvent*>(event);
        const AliESDVertex *esdVertex = ESDevent->GetPrimaryVertexTracks();
        if(!esdVertex){printf("ERROR: vertex not available\n"); return kFALSE;}
        if(esdVertex->GetNContributors()<1) {
            // SPD vertex
            esdVertex = ESDevent->GetPrimaryVertexSPD();
        }
        //     AliESDVertex *esdVertex = dynamic_cast<AliESDVertex*> (vertex);
        if(!esdVertex->GetStatus()){return kFALSE;}
        Double_t zRes = esdVertex->GetZRes();
        if (zRes > requiredZResolution) return kFALSE;

        const AliESDVertex *vertexSPD = ESDevent->GetPrimaryVertexSPD();
        // always check for SPD vertex
        const AliESDVertex * trkVertex = ESDevent->GetPrimaryVertexTracks();
        if(!vertexSPD) return kFALSE;
        if(!vertexSPD->GetStatus()) return kFALSE;
        if(!trkVertex->GetStatus()) return kFALSE;
        if (vertexSPD->IsFromVertexerZ() && !(vertexSPD->GetDispersion()<0.04 && vertexSPD->GetZRes()<0.25)) return kFALSE;
        //if (vertexSPD->IsFromVertexerZ() && vertexSPD->GetDispersion() > 0.04) return kFALSE; /// vertexSPD->GetDispersion() > 0.02 to 0.04
        if ((TMath::Abs(vertexSPD->GetZ() - trkVertex->GetZ())>0.5)) return kFALSE;
    }

    else{/*
            Float_t requiredZResolution = 1000;
            AliAODEvent* AODevent = dynamic_cast<AliAODEvent*>(event);
            const AliAODVertex *aodVertex = (AliAODVertex*)AODevent->GetVertex(0);
            if(!aodVertex){printf("ERROR: vertex not available\n"); return kFALSE;}

            if(aodVertex->GetNContributors()<1) {
        // SPD vertex
        aodVertex = AODevent->GetPrimaryVertexSPD();
        }
        if(!aodVertex->GetStatus()){return kFALSE;}
        //  Double_t zRes = aodVertex->GetZRes(); There didn't find an analogous function to GetZRes() for AODs.
        //if (zRes > requiredZResolution) return kFALSE;

        const AliAODVertex *vertexSPD = AODevent->GetPrimaryVertexSPD();
        // always check for SPD vertex
        const AliAODVertex * trkVertex = (AliAODVertex*)AODevent->GetPrimaryVertexTracks();
        if(!vertexSPD) return kFALSE;
        if(!vertexSPD->GetStatus()) return kFALSE;
        if(!trkVertex->GetStatus()) return kFALSE;
        //if (vertexSPD->IsFromVertexerZ() && !(vertexSPD->GetDispersion()<0.04 && vertexSPD->GetZRes()<0.25)) return kFALSE; There didn't find an analogous function to GetZRes() for AODs.
        //if (vertexSPD->IsFromVertexerZ() && vertexSPD->GetDispersion() > 0.04) return kFALSE; /// vertexSPD->GetDispersion() > 0.02 to 0.04
        if ((TMath::Abs(vertexSPD->GetZ() - trkVertex->GetZ())>0.5)) return kFALSE;
        */ }
        return kTRUE;
}

/// Function to determine if MC Event is in INEL>0 Class
/// select INEL>0 events with at least
/// one prompt (MC primary) particle in acceptance
/// pT>0, |eta|<1.0 for normalization
Bool_t AlidNdPtUnifiedAnalysisTask::IsMCEventINEL0(AliMCEvent* mcEvent, Double_t ptmin, Double_t etarange){
    if(!mcEvent) return kFALSE;
    // AliMCEvent* fmcEvent = static_cast<AliMCEvent*>(mcEvent);
    AliStack* stack = mcEvent->Stack();
    if(!stack) return kFALSE;

    Int_t count = 0;
    for (Int_t iMc = 0; iMc < stack->GetNtrack(); ++iMc)
    {
        TParticle* particle = stack->Particle(iMc);
        if (!particle) continue;

        // only charged particles
        if(!particle->GetPDG()) continue;
        Double_t charge = particle->GetPDG()->Charge()/3.;
        if(charge == 0) continue;

        // physical primary
        Bool_t prim = stack->IsPhysicalPrimary(iMc);
        if(!prim) continue;

        if(particle->Pt() < ptmin) continue;
        if(TMath::Abs(particle->Eta()) > etarange) continue;

        count++;
    }

    if(count > 0) return kTRUE;
    else return kFALSE;
}
