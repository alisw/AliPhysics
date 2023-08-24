#include "AliAnalysisTaskAntiProtons_vs_RT_pPb.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTask.h"
#include "AliESDtrackCuts.h"
#include "AliPIDResponse.h"
#include "TLorentzVector.h"
#include "AliCentrality.h"
#include "AliMCParticle.h"
#include "AliEventCuts.h"
#include "TDatabasePDG.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "THnSparse.h"
#include "AliVVZERO.h"
#include "TObjArray.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TRandom.h"
#include "THnBase.h"
#include "TChain.h"
#include "TMath.h"
#include "TTree.h"
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;
ClassImp(AliAnalysisTaskAntiProtons_vs_RT_pPb)

//____________________________________________________________________________________________________________________________
AliAnalysisTaskAntiProtons_vs_RT_pPb::AliAnalysisTaskAntiProtons_vs_RT_pPb():
AliAnalysisTaskSE(), 
fESDevent(nullptr),
fMCevent(nullptr),
fMCEventHandler(nullptr),
fESDtrackCuts_AntiProton(nullptr),
fESDtrackCuts_Pions(nullptr),
fESDtrackCuts_LeadingTrack(nullptr),
fESDtrackCuts_TransverseMult(nullptr),
fPIDResponse(nullptr),
fOutputList(nullptr),
fQAList(nullptr),
fTriggerType(AliVEvent::kINT7),
fMultMin(0.0),
fMultMax(100.0),
fAverage_Nch_Transv(13.73),
fIsMC(kFALSE),
fPt_min_leading(3),
fIsITSrecalib(kFALSE),
fIsPion(kFALSE),
hMean(nullptr),
hWidth(nullptr),
fESDEventSelection(),
hNumberOfEvents(nullptr),
hMultPercentile(nullptr),
hMultTransverse(nullptr),
hMultToward(nullptr),
hMultAway(nullptr),
hRtDistribution(nullptr),
hEventsWithLeadingTrack(nullptr),
hNumberOfAntiProtons(nullptr),
hNumberOfPions(nullptr),
hNchTransv_NchTot(nullptr),
hRt_NchTot(nullptr),
hNchTransv_MultPercentile(nullptr),
hRt_MultPercentile(nullptr),
hnsigmaITS_antiprotons(nullptr),
hnsigmaTPC_antiprotons_Toward(nullptr),
hnsigmaTPC_antiprotons_Away(nullptr),
hnsigmaTPC_antiprotons_Transverse(nullptr),
hnsigmaTPC_antiprotons_Toward_noITSpresel(nullptr),
hnsigmaTPC_antiprotons_Away_noITSpresel(nullptr),
hnsigmaTPC_antiprotons_Transverse_noITSpresel(nullptr),
hnsigmaTOF_antiprotons_Toward(nullptr), 
hnsigmaTOF_antiprotons_Away(nullptr), 
hnsigmaTOF_antiprotons_Transverse(nullptr),
hDCAxy_antiprotons_Toward(nullptr),
hDCAxy_antiprotons_Away(nullptr),
hDCAxy_antiprotons_Transverse(nullptr),
hnsigmaTPC_antiprotons_Syst(nullptr),
hnsigmaTOF_antiprotons_Syst(nullptr),
hDCAxy_antiprotons_Syst(nullptr),
hnsigmaTOF_pos_pions_Toward(nullptr), 
hnsigmaTOF_pos_pions_Away(nullptr), 
hnsigmaTOF_pos_pions_Transverse(nullptr),
hnsigmaTOF_neg_pions_Toward(nullptr), 
hnsigmaTOF_neg_pions_Away(nullptr), 
hnsigmaTOF_neg_pions_Transverse(nullptr),
hDCAxy_pos_pions_Toward(nullptr),
hDCAxy_pos_pions_Away(nullptr),
hDCAxy_pos_pions_Transverse(nullptr),
hDCAxy_neg_pions_Toward(nullptr),
hDCAxy_neg_pions_Away(nullptr),
hDCAxy_neg_pions_Transverse(nullptr),
hnsigmaTOF_pos_pions_Syst(nullptr),
hDCAxy_pos_pions_Syst(nullptr),
hnsigmaTOF_neg_pions_Syst(nullptr),
hDCAxy_neg_pions_Syst(nullptr),
h_antiprotons_Gen(nullptr), 
hnsigmaTPC_antiprotons_Rec(nullptr),
hnsigmaTOF_antiprotons_Rec(nullptr),
hDCAxy_antiprotons_prim(nullptr),
hDCAxy_antiprotons_sec(nullptr),
hnsigmaTPC_antiprotons_Rec_Syst(nullptr),
hnsigmaTOF_antiprotons_Rec_Syst(nullptr), 
hDCAxy_antiprotons_prim_Syst(nullptr),
hDCAxy_antiprotons_sec_Syst(nullptr), 
h_pos_pions_Gen(nullptr),
h_neg_pions_Gen(nullptr),
hnsigmaTOF_pos_pions_Rec(nullptr),
hnsigmaTOF_neg_pions_Rec(nullptr),
hDCAxy_pos_pions_prim(nullptr),
hDCAxy_pos_pions_sec(nullptr),
hDCAxy_neg_pions_prim(nullptr),
hDCAxy_neg_pions_sec(nullptr),
hnsigmaTOF_pos_pions_Rec_Syst(nullptr),
hnsigmaTOF_neg_pions_Rec_Syst(nullptr),
hDCAxy_pos_pions_prim_Syst(nullptr),
hDCAxy_neg_pions_prim_Syst(nullptr),
hDCAxy_pos_pions_sec_Syst(nullptr),
hDCAxy_neg_pions_sec_Syst(nullptr)
{}
//____________________________________________________________________________________________________________________________
AliAnalysisTaskAntiProtons_vs_RT_pPb::AliAnalysisTaskAntiProtons_vs_RT_pPb(const char *name):
AliAnalysisTaskSE(name), 
fESDevent(nullptr),
fMCevent(nullptr),
fMCEventHandler(nullptr),
fESDtrackCuts_AntiProton(nullptr),
fESDtrackCuts_Pions(nullptr),
fESDtrackCuts_LeadingTrack(nullptr),
fESDtrackCuts_TransverseMult(nullptr),
fPIDResponse(nullptr),
fOutputList(nullptr),
fQAList(nullptr),
fTriggerType(AliVEvent::kINT7),
fMultMin(0.0),
fMultMax(100.0),
fAverage_Nch_Transv(13.73),
fIsMC(kFALSE),
fPt_min_leading(3),
fIsITSrecalib(kFALSE),
fIsPion(kFALSE),
hMean(nullptr),
hWidth(nullptr),
fESDEventSelection(),
hNumberOfEvents(nullptr),
hMultPercentile(nullptr),
hMultTransverse(nullptr),
hMultToward(nullptr),
hMultAway(nullptr),
hRtDistribution(nullptr),
hEventsWithLeadingTrack(nullptr),
hNumberOfAntiProtons(nullptr),
hNumberOfPions(nullptr),
hNchTransv_NchTot(nullptr),
hRt_NchTot(nullptr),
hNchTransv_MultPercentile(nullptr),
hRt_MultPercentile(nullptr),
hnsigmaITS_antiprotons(nullptr),
hnsigmaTPC_antiprotons_Toward(nullptr),
hnsigmaTPC_antiprotons_Away(nullptr),
hnsigmaTPC_antiprotons_Transverse(nullptr),
hnsigmaTPC_antiprotons_Toward_noITSpresel(nullptr),
hnsigmaTPC_antiprotons_Away_noITSpresel(nullptr),
hnsigmaTPC_antiprotons_Transverse_noITSpresel(nullptr),
hnsigmaTOF_antiprotons_Toward(nullptr), 
hnsigmaTOF_antiprotons_Away(nullptr), 
hnsigmaTOF_antiprotons_Transverse(nullptr),
hDCAxy_antiprotons_Toward(nullptr),
hDCAxy_antiprotons_Away(nullptr),
hDCAxy_antiprotons_Transverse(nullptr),
hnsigmaTPC_antiprotons_Syst(nullptr),
hnsigmaTOF_antiprotons_Syst(nullptr),
hDCAxy_antiprotons_Syst(nullptr),
hnsigmaTOF_pos_pions_Toward(nullptr), 
hnsigmaTOF_pos_pions_Away(nullptr), 
hnsigmaTOF_pos_pions_Transverse(nullptr),
hnsigmaTOF_neg_pions_Toward(nullptr), 
hnsigmaTOF_neg_pions_Away(nullptr), 
hnsigmaTOF_neg_pions_Transverse(nullptr),
hDCAxy_pos_pions_Toward(nullptr),
hDCAxy_pos_pions_Away(nullptr),
hDCAxy_pos_pions_Transverse(nullptr),
hDCAxy_neg_pions_Toward(nullptr),
hDCAxy_neg_pions_Away(nullptr),
hDCAxy_neg_pions_Transverse(nullptr),
hnsigmaTOF_pos_pions_Syst(nullptr),
hDCAxy_pos_pions_Syst(nullptr),
hnsigmaTOF_neg_pions_Syst(nullptr),
hDCAxy_neg_pions_Syst(nullptr),
h_antiprotons_Gen(nullptr), 
hnsigmaTPC_antiprotons_Rec(nullptr),
hnsigmaTOF_antiprotons_Rec(nullptr),
hDCAxy_antiprotons_prim(nullptr),
hDCAxy_antiprotons_sec(nullptr),
hnsigmaTPC_antiprotons_Rec_Syst(nullptr),
hnsigmaTOF_antiprotons_Rec_Syst(nullptr), 
hDCAxy_antiprotons_prim_Syst(nullptr),
hDCAxy_antiprotons_sec_Syst(nullptr),
h_pos_pions_Gen(nullptr),
h_neg_pions_Gen(nullptr),
hnsigmaTOF_pos_pions_Rec(nullptr),
hnsigmaTOF_neg_pions_Rec(nullptr),
hDCAxy_pos_pions_prim(nullptr),
hDCAxy_pos_pions_sec(nullptr),
hDCAxy_neg_pions_prim(nullptr),
hDCAxy_neg_pions_sec(nullptr),
hnsigmaTOF_pos_pions_Rec_Syst(nullptr),
hnsigmaTOF_neg_pions_Rec_Syst(nullptr),
hDCAxy_pos_pions_prim_Syst(nullptr),
hDCAxy_neg_pions_prim_Syst(nullptr),
hDCAxy_pos_pions_sec_Syst(nullptr),
hDCAxy_neg_pions_sec_Syst(nullptr)

{
    DefineInput (0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
}
//____________________________________________________________________________________________________________________________
AliAnalysisTaskAntiProtons_vs_RT_pPb::~AliAnalysisTaskAntiProtons_vs_RT_pPb() {
    fOutputList->Clear();

    delete fESDevent;
    delete fMCevent;
    delete fMCEventHandler;
    delete fESDtrackCuts_AntiProton;
    delete fESDtrackCuts_Pions;
    delete fESDtrackCuts_LeadingTrack;
    delete fESDtrackCuts_TransverseMult;
    delete fPIDResponse;
    delete fOutputList;
    delete fQAList;
    delete hMean;
    delete hWidth;
    delete hNumberOfEvents;
    delete hMultPercentile;
    delete hMultTransverse;
    delete hMultToward;
    delete hMultAway;
    delete hRtDistribution;
    delete hEventsWithLeadingTrack;
    delete hNumberOfAntiProtons;
    delete hNumberOfPions;
    delete hNchTransv_NchTot;
    delete hRt_NchTot;
    delete hNchTransv_MultPercentile;
    delete hRt_MultPercentile;
    delete hnsigmaITS_antiprotons;
    delete hnsigmaTPC_antiprotons_Toward;
    delete hnsigmaTPC_antiprotons_Away;
    delete hnsigmaTPC_antiprotons_Transverse;
    delete hnsigmaTPC_antiprotons_Toward_noITSpresel;
    delete hnsigmaTPC_antiprotons_Away_noITSpresel;
    delete hnsigmaTPC_antiprotons_Transverse_noITSpresel;
    delete hnsigmaTOF_antiprotons_Toward;
    delete hnsigmaTOF_antiprotons_Away;
    delete hnsigmaTOF_antiprotons_Transverse;
    delete hDCAxy_antiprotons_Toward;
    delete hDCAxy_antiprotons_Away;
    delete hDCAxy_antiprotons_Transverse;
    delete hnsigmaTPC_antiprotons_Syst;
    delete hnsigmaTOF_antiprotons_Syst;
    delete hDCAxy_antiprotons_Syst;
    delete hnsigmaTOF_pos_pions_Toward;
    delete hnsigmaTOF_pos_pions_Away;
    delete hnsigmaTOF_pos_pions_Transverse;
    delete hnsigmaTOF_neg_pions_Toward;
    delete hnsigmaTOF_neg_pions_Away;
    delete hnsigmaTOF_neg_pions_Transverse;
    delete hDCAxy_pos_pions_Toward;
    delete hDCAxy_pos_pions_Away;
    delete hDCAxy_pos_pions_Transverse;
    delete hDCAxy_neg_pions_Toward;
    delete hDCAxy_neg_pions_Away;
    delete hDCAxy_neg_pions_Transverse;
    delete hnsigmaTOF_pos_pions_Syst;
    delete hDCAxy_pos_pions_Syst;
    delete hnsigmaTOF_neg_pions_Syst;
    delete hDCAxy_neg_pions_Syst;
    delete h_antiprotons_Gen;
    delete hnsigmaTPC_antiprotons_Rec;
    delete hnsigmaTOF_antiprotons_Rec;
    delete hDCAxy_antiprotons_prim;    
    delete hDCAxy_antiprotons_sec;
    delete hnsigmaTPC_antiprotons_Rec_Syst;
    delete hnsigmaTOF_antiprotons_Rec_Syst;
    delete hDCAxy_antiprotons_prim_Syst;
    delete hDCAxy_antiprotons_sec_Syst;
    delete h_pos_pions_Gen;
    delete h_neg_pions_Gen;
    delete hnsigmaTOF_pos_pions_Rec;
    delete hnsigmaTOF_neg_pions_Rec;
    delete hDCAxy_pos_pions_prim;
    delete hDCAxy_pos_pions_sec;
    delete hDCAxy_neg_pions_prim;
    delete hDCAxy_neg_pions_sec;
    delete hnsigmaTOF_pos_pions_Rec_Syst;
    delete hnsigmaTOF_neg_pions_Rec_Syst;
    delete hDCAxy_pos_pions_prim_Syst;
    delete hDCAxy_neg_pions_prim_Syst;
    delete hDCAxy_pos_pions_sec_Syst;
    delete hDCAxy_neg_pions_sec_Syst;
}
//____________________________________________________________________________________________________________________________
void AliAnalysisTaskAntiProtons_vs_RT_pPb::UserCreateOutputObjects() {

    //Create Output List
    fOutputList = new TList();
    fQAList     = new TList();
    fOutputList -> SetOwner();
    fQAList     -> SetOwner();

    //Event Selection
    fESDEventSelection.AddQAplotsToList(fQAList);
    fESDEventSelection.SetManualMode();
    fESDEventSelection.fRequireTrackVertex = true;
    fESDEventSelection.fMinVtz = -10.f;
    fESDEventSelection.fMaxVtz = 10.f;
    fESDEventSelection.fMaxDeltaSpdTrackAbsolute = 0.5f;
    fESDEventSelection.fMaxResolutionSPDvertex = 0.25f;
    fESDEventSelection.fTriggerMask = (AliVEvent::kINT7);
    fESDEventSelection.fRejectDAQincomplete = true;
    fESDEventSelection.fSPDpileupMinContributors = 3;
    fESDEventSelection.fSPDpileupMinZdist = 0.8;
    fESDEventSelection.fSPDpileupNsigmaZdist = 3.;
    fESDEventSelection.fSPDpileupNsigmaDiamXY = 2.;
    fESDEventSelection.fSPDpileupNsigmaDiamZ = 5.;
    fESDEventSelection.fTrackletBGcut = true;

    //Event Counter and Multiplicity Percentile
    hNumberOfEvents = new TH1F ("hNumberOfEvents","",20,0,20);
    hMultPercentile = new TH1F ("hMultPercentile","",200,0,100);
    hNumberOfEvents -> Sumw2();
    hMultPercentile -> Sumw2();
    fOutputList -> Add(hNumberOfEvents);
    fOutputList -> Add(hMultPercentile);

    //Multiplicity in Azimuthal Regions
    if (!fIsMC)  {
        hMultTransverse         = new TH1I ("hMultTransverse","",200,0,200);
        hMultToward             = new TH1I ("hMultToward","",200,0,200);
        hMultAway               = new TH1I ("hMultAway","",200,0,200);
        hRtDistribution         = new TH1F ("hRtDistribution","",200,0,20);
        hEventsWithLeadingTrack = new TH1F ("hEventsWithLeadingTrack","",10,0,10);
        hMultTransverse         -> Sumw2();
        hMultToward             -> Sumw2();
        hMultAway               -> Sumw2();
        hRtDistribution         -> Sumw2();
        hEventsWithLeadingTrack -> Sumw2();
        fOutputList -> Add(hMultTransverse);
        fOutputList -> Add(hMultToward);
        fOutputList -> Add(hMultAway);
        fOutputList -> Add(hRtDistribution);
        fOutputList -> Add(hEventsWithLeadingTrack);

        if(!fIsPion){
            hNumberOfAntiProtons = new TH1I ("hNumberOfAntiProtons","",20,0,20);
            hNumberOfAntiProtons -> Sumw2();
            fOutputList -> Add(hNumberOfAntiProtons);
        }

        if(fIsPion) {
            hNumberOfPions = new TH1I ("hNumberOfPosPions","",20,0,20);
            hNumberOfPions -> Sumw2();
            fOutputList -> Add(hNumberOfPions);
        }
    }

    //Correlations between Transverse and Integrated Mult
    if (!fIsMC)  {
        
        hNchTransv_NchTot         = new TH2F ("hNchTransv_NchTot","",200,0,200,200,0,200);
        hRt_NchTot                = new TH2F ("hRt_NchTot","",200,0,10,200,0,200);
        hNchTransv_MultPercentile = new TH2F ("hNchTransv_MultPercentile","",200,0,200,200,0,100);
        hRt_MultPercentile        = new TH2F ("hRt_MultPercentile","",200,0,10,200,0,100);
        hNchTransv_NchTot         -> Sumw2();
        hRt_NchTot                -> Sumw2();
        hNchTransv_MultPercentile -> Sumw2();
        hRt_MultPercentile        -> Sumw2();
        fOutputList -> Add (hNchTransv_NchTot);
        fOutputList -> Add (hRt_NchTot);
        fOutputList -> Add (hNchTransv_MultPercentile);
        fOutputList -> Add (hRt_MultPercentile);
    }

    //Arrays
    Double_t Nch_Tr[]   = { 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 5.0, 10.0 };
    Double_t pt_TPC[]   = { 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 1.0};
    Double_t pt_TOF[]   = { 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 5.0 };
    Double_t pt_DCA[]   = { 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.4};
    Double_t dca_xy[]   = { -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.18, -0.16, -0.14, -0.12, -0.10, -0.08, -0.06, -0.04, -0.02, 0.0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
    Double_t p_ITS[]    = { 0.3, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.9, 1.0};
    Double_t eta_ITS[]  = { -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8};

    //Arrays Dimensions
    Int_t nBins_Nch_Tr   = sizeof(Nch_Tr)/sizeof(Double_t)-1;
    Int_t nBins_pt_TPC   = sizeof(pt_TPC)/sizeof(Double_t)-1;
    Int_t nBins_pt_TOF   = sizeof(pt_TOF)/sizeof(Double_t)-1;
    Int_t nBins_pt_DCA   = sizeof(pt_DCA)/sizeof(Double_t)-1;
    Int_t nBins_dca_xy   = sizeof(dca_xy)/sizeof(Double_t)-1;
    Int_t nBins_p_ITS    = sizeof(p_ITS) /sizeof(Double_t)-1;
    Int_t nBins_eta_ITS  = sizeof(eta_ITS)/sizeof(Double_t)-1;

    if(!fIsPion) {
        //definition of nsigmaITS plots, filled both for data and MC and if fItITSrecalibrated is true or false  (if false used to build recalibration maps offline)
        //*************   nsigma_{ITS},              p_ITS,               eta_ITS*************//
        Int_t    bins_ITS[3] = {   200,        nBins_p_ITS,          nBins_eta_ITS};
        Double_t xmin_ITS[3] = { -10.0,           p_ITS[0],             eta_ITS[0]};
        Double_t xmax_ITS[3] = { +10.0, p_ITS[nBins_p_ITS], eta_ITS[nBins_eta_ITS]};
        //************************************************************************************//

        hnsigmaITS_antiprotons = new THnSparseF("hnsigmaITS_antiprotons", "", 3, bins_ITS, xmin_ITS, xmax_ITS);
        hnsigmaITS_antiprotons -> GetAxis(1) -> Set(nBins_p_ITS, p_ITS);
        hnsigmaITS_antiprotons -> GetAxis(2) -> Set(nBins_eta_ITS, eta_ITS);
        hnsigmaITS_antiprotons -> Sumw2();
        fOutputList ->Add (hnsigmaITS_antiprotons);

        if ((!fIsMC)&&(fIsITSrecalib))  {
        
            //***************************  Nch_Transverse,  nsigma_{TPC},                 pt_TPC ***************************************************
            Int_t    bins_TPC[3] = {         nBins_Nch_Tr,           200,           nBins_pt_TPC };
            Double_t xmin_TPC[3] = {            Nch_Tr[0],          -10.0,              pt_TPC[0] };
            Double_t xmax_TPC[3] = { Nch_Tr[nBins_Nch_Tr],          +10.0,   pt_TPC[nBins_pt_TPC] };
            //***************************************************************************************************************************************

            hnsigmaTPC_antiprotons_Toward     = new THnSparseF ("hnsigmaTPC_antiprotons_Toward","",3, bins_TPC, xmin_TPC, xmax_TPC);
            hnsigmaTPC_antiprotons_Away       = new THnSparseF ("hnsigmaTPC_antiprotons_Away","",3, bins_TPC, xmin_TPC, xmax_TPC);
            hnsigmaTPC_antiprotons_Transverse = new THnSparseF ("hnsigmaTPC_antiprotons_Transverse","",3, bins_TPC, xmin_TPC, xmax_TPC);
            hnsigmaTPC_antiprotons_Toward     -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hnsigmaTPC_antiprotons_Away       -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hnsigmaTPC_antiprotons_Transverse -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hnsigmaTPC_antiprotons_Toward     -> GetAxis(2) -> Set(nBins_pt_TPC, pt_TPC);
            hnsigmaTPC_antiprotons_Away       -> GetAxis(2) -> Set(nBins_pt_TPC, pt_TPC);
            hnsigmaTPC_antiprotons_Transverse -> GetAxis(2) -> Set(nBins_pt_TPC, pt_TPC);
            hnsigmaTPC_antiprotons_Toward     -> Sumw2();
            hnsigmaTPC_antiprotons_Away       -> Sumw2();
            hnsigmaTPC_antiprotons_Transverse -> Sumw2();
            fOutputList -> Add (hnsigmaTPC_antiprotons_Toward);
            fOutputList -> Add (hnsigmaTPC_antiprotons_Away);
            fOutputList -> Add (hnsigmaTPC_antiprotons_Transverse);

            hnsigmaTPC_antiprotons_Toward_noITSpresel     = new THnSparseF ("hnsigmaTPC_antiprotons_Toward_noITSpresel","",3, bins_TPC, xmin_TPC, xmax_TPC);
            hnsigmaTPC_antiprotons_Away_noITSpresel       = new THnSparseF ("hnsigmaTPC_antiprotons_Away_noITSpresel","",3, bins_TPC, xmin_TPC, xmax_TPC);
            hnsigmaTPC_antiprotons_Transverse_noITSpresel = new THnSparseF ("hnsigmaTPC_antiprotons_Transverse_noITSpresel","",3, bins_TPC, xmin_TPC, xmax_TPC);
            hnsigmaTPC_antiprotons_Toward_noITSpresel     -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hnsigmaTPC_antiprotons_Away_noITSpresel       -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hnsigmaTPC_antiprotons_Transverse_noITSpresel -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hnsigmaTPC_antiprotons_Toward_noITSpresel     -> GetAxis(2) -> Set(nBins_pt_TPC, pt_TPC);
            hnsigmaTPC_antiprotons_Away_noITSpresel       -> GetAxis(2) -> Set(nBins_pt_TPC, pt_TPC);
            hnsigmaTPC_antiprotons_Transverse_noITSpresel -> GetAxis(2) -> Set(nBins_pt_TPC, pt_TPC);
            hnsigmaTPC_antiprotons_Toward_noITSpresel     -> Sumw2();
            hnsigmaTPC_antiprotons_Away_noITSpresel       -> Sumw2();
            hnsigmaTPC_antiprotons_Transverse_noITSpresel -> Sumw2();
            fOutputList -> Add (hnsigmaTPC_antiprotons_Toward_noITSpresel);
            fOutputList -> Add (hnsigmaTPC_antiprotons_Away_noITSpresel);
            fOutputList -> Add (hnsigmaTPC_antiprotons_Transverse_noITSpresel);

            //***************************  Nch_Transverse,  nsigma_{TOF},                 pt_TOF ***************************************************
            Int_t    bins_TOF[3] = {         nBins_Nch_Tr,           200,           nBins_pt_TOF };
            Double_t xmin_TOF[3] = {            Nch_Tr[0],         -10.0,              pt_TOF[0] };
            Double_t xmax_TOF[3] = { Nch_Tr[nBins_Nch_Tr],         +10.0,   pt_TOF[nBins_pt_TOF] };
            //***************************************************************************************************************************************

            hnsigmaTOF_antiprotons_Toward     = new THnSparseF ("hnsigmaTOF_antiprotons_Toward","",3, bins_TOF, xmin_TOF, xmax_TOF);
            hnsigmaTOF_antiprotons_Away       = new THnSparseF ("hnsigmaTOF_antiprotons_Away","",3, bins_TOF, xmin_TOF, xmax_TOF);
            hnsigmaTOF_antiprotons_Transverse = new THnSparseF ("hnsigmaTOF_antiprotons_Transverse","",3, bins_TOF, xmin_TOF, xmax_TOF);
            hnsigmaTOF_antiprotons_Toward     -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hnsigmaTOF_antiprotons_Away       -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hnsigmaTOF_antiprotons_Transverse -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hnsigmaTOF_antiprotons_Toward     -> GetAxis(2) -> Set(nBins_pt_TOF, pt_TOF);
            hnsigmaTOF_antiprotons_Away       -> GetAxis(2) -> Set(nBins_pt_TOF, pt_TOF);
            hnsigmaTOF_antiprotons_Transverse -> GetAxis(2) -> Set(nBins_pt_TOF, pt_TOF);
            hnsigmaTOF_antiprotons_Toward     -> Sumw2();
            hnsigmaTOF_antiprotons_Away       -> Sumw2();
            hnsigmaTOF_antiprotons_Transverse -> Sumw2();
            fOutputList -> Add (hnsigmaTOF_antiprotons_Toward);
            fOutputList -> Add (hnsigmaTOF_antiprotons_Away);
            fOutputList -> Add (hnsigmaTOF_antiprotons_Transverse);

            //***************************  Nch_Transverse,      DCA_{xy},                 pt_DCA ***************************************************
            Int_t    bins_DCA[3] = {         nBins_Nch_Tr,  nBins_dca_xy,           nBins_pt_DCA };
            Double_t xmin_DCA[3] = {            Nch_Tr[0],          -1.0,              pt_DCA[0] };
            Double_t xmax_DCA[3] = { Nch_Tr[nBins_Nch_Tr],          +1.0,   pt_DCA[nBins_pt_DCA] };
            //***************************************************************************************************************************************

            hDCAxy_antiprotons_Toward     = new THnSparseF ("hDCAxy_antiprotons_Toward","",3, bins_DCA, xmin_DCA, xmax_DCA);
            hDCAxy_antiprotons_Away       = new THnSparseF ("hDCAxy_antiprotons_Away","",3, bins_DCA, xmin_DCA, xmax_DCA);
            hDCAxy_antiprotons_Transverse = new THnSparseF ("hDCAxy_antiprotons_Transverse","",3, bins_DCA, xmin_DCA, xmax_DCA);
            hDCAxy_antiprotons_Toward     -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hDCAxy_antiprotons_Away       -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hDCAxy_antiprotons_Transverse -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hDCAxy_antiprotons_Toward     -> GetAxis(1) -> Set(nBins_dca_xy, dca_xy);
            hDCAxy_antiprotons_Away       -> GetAxis(1) -> Set(nBins_dca_xy, dca_xy);
            hDCAxy_antiprotons_Transverse -> GetAxis(1) -> Set(nBins_dca_xy, dca_xy);
            hDCAxy_antiprotons_Toward     -> GetAxis(2) -> Set(nBins_pt_DCA, pt_DCA);
            hDCAxy_antiprotons_Away       -> GetAxis(2) -> Set(nBins_pt_DCA, pt_DCA);
            hDCAxy_antiprotons_Transverse -> GetAxis(2) -> Set(nBins_pt_DCA, pt_DCA);
            hDCAxy_antiprotons_Toward     -> Sumw2();
            hDCAxy_antiprotons_Away       -> Sumw2();
            hDCAxy_antiprotons_Transverse -> Sumw2();
            fOutputList -> Add(hDCAxy_antiprotons_Toward);
            fOutputList -> Add(hDCAxy_antiprotons_Away);
            fOutputList -> Add(hDCAxy_antiprotons_Transverse);


            //********************************  Nch_Transverse,  nsigma_{TPC},                 pt_TPC, isyst ****************************************
            Int_t    bins_TPC_syst[4] = {         nBins_Nch_Tr,           100,           nBins_pt_TPC,    50 };
            Double_t xmin_TPC_syst[4] = {            Nch_Tr[0],          -10.0,              pt_TPC[0],     0 };
            Double_t xmax_TPC_syst[4] = { Nch_Tr[nBins_Nch_Tr],          +10.0,   pt_TPC[nBins_pt_TPC],    50 };
            //***************************************************************************************************************************************

            hnsigmaTPC_antiprotons_Syst = new THnSparseF ("hnsigmaTPC_antiprotons_Syst","",4, bins_TPC_syst, xmin_TPC_syst, xmax_TPC_syst);
            hnsigmaTPC_antiprotons_Syst -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hnsigmaTPC_antiprotons_Syst -> GetAxis(2) -> Set(nBins_pt_TPC, pt_TPC);
            hnsigmaTPC_antiprotons_Syst -> Sumw2();
            fOutputList -> Add (hnsigmaTPC_antiprotons_Syst);


            //********************************  Nch_Transverse,  nsigma_{TOF},                  pt_TOF, isyst ****************************************
            Int_t    bins_TOF_syst[4] = {         nBins_Nch_Tr,           200,            nBins_pt_TOF,    50 };
            Double_t xmin_TOF_syst[4] = {            Nch_Tr[0],          -10.0,              pt_TOF[0],     0 };
            Double_t xmax_TOF_syst[4] = { Nch_Tr[nBins_Nch_Tr],          +10.0,   pt_TOF[nBins_pt_TOF],    50 };
            //***************************************************************************************************************************************

            hnsigmaTOF_antiprotons_Syst = new THnSparseF ("hnsigmaTOF_antiprotons_Syst","",4, bins_TOF_syst, xmin_TOF_syst, xmax_TOF_syst);
            hnsigmaTOF_antiprotons_Syst -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hnsigmaTOF_antiprotons_Syst -> GetAxis(2) -> Set(nBins_pt_TOF, pt_TOF);
            hnsigmaTOF_antiprotons_Syst -> Sumw2();
            fOutputList -> Add (hnsigmaTOF_antiprotons_Syst);


            //***************************       Nch_Transverse,      DCA_{xy},                 pt_DCA,   isyst *******************************************
            Int_t    bins_DCA_syst[4] = {         nBins_Nch_Tr,  nBins_dca_xy,           nBins_pt_DCA,    50 };
            Double_t xmin_DCA_syst[4] = {            Nch_Tr[0],          -1.0,              pt_DCA[0],     0 };
            Double_t xmax_DCA_syst[4] = { Nch_Tr[nBins_Nch_Tr],          +1.0,   pt_DCA[nBins_pt_DCA],    50 };
            //***************************************************************************************************************************************

            hDCAxy_antiprotons_Syst = new THnSparseF ("hDCAxy_antiprotons_Syst","",4, bins_DCA_syst, xmin_DCA_syst, xmax_DCA_syst);
            hDCAxy_antiprotons_Syst -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hDCAxy_antiprotons_Syst -> GetAxis(1) -> Set(nBins_dca_xy, dca_xy);
            hDCAxy_antiprotons_Syst -> GetAxis(2) -> Set(nBins_pt_DCA, pt_DCA);
            hDCAxy_antiprotons_Syst -> Sumw2();
            fOutputList -> Add(hDCAxy_antiprotons_Syst);
        }

        //Generated & Reconstructed p_{T} Spectra
        if(fIsMC){
            //DCA_{xy} histograms
            hDCAxy_antiprotons_prim = new TH2F("hDCAxy_antiprotons_prim", "", nBins_pt_DCA, pt_DCA, nBins_dca_xy, dca_xy);   
            hDCAxy_antiprotons_sec  = new TH2F("hDCAxy_antiprotons_sec", "", nBins_pt_DCA, pt_DCA, nBins_dca_xy, dca_xy);
            hDCAxy_antiprotons_prim -> Sumw2();
            hDCAxy_antiprotons_sec  -> Sumw2();
            fOutputList -> Add(hDCAxy_antiprotons_prim);
            fOutputList -> Add(hDCAxy_antiprotons_sec);
        }

        if (fIsMC && fIsITSrecalib)  {
            h_antiprotons_Gen          = new TH1F ("h_antiprotons_Gen","",940,0.3,5.0);
            hnsigmaTPC_antiprotons_Rec = new TH2F ("hnsigmaTPC_antiprotons_Rec","",940,0.3,5.0,100,-5,5);
            hnsigmaTOF_antiprotons_Rec = new TH2F ("hnsigmaTOF_antiprotons_Rec","",940,0.3,5.0,100,-5,5);

            h_antiprotons_Gen          -> Sumw2();
            hnsigmaTPC_antiprotons_Rec -> Sumw2();
            hnsigmaTOF_antiprotons_Rec -> Sumw2();
            fOutputList -> Add (h_antiprotons_Gen);
            fOutputList -> Add (hnsigmaTPC_antiprotons_Rec);
            fOutputList -> Add (hnsigmaTOF_antiprotons_Rec);

            //Histograms for Syst. Uncertainties
            hnsigmaTPC_antiprotons_Rec_Syst = new TH2F ("hnsigmaTPC_antiprotons_Rec_Syst","",940,0.3,5.0,50,0,50);
            hnsigmaTOF_antiprotons_Rec_Syst = new TH2F ("hnsigmaTOF_antiprotons_Rec_Syst","",940,0.3,5.0,50,0,50);
            hnsigmaTPC_antiprotons_Rec_Syst -> Sumw2();
            hnsigmaTOF_antiprotons_Rec_Syst -> Sumw2();
            fOutputList -> Add(hnsigmaTPC_antiprotons_Rec_Syst);
            fOutputList -> Add(hnsigmaTOF_antiprotons_Rec_Syst);

            //***************************    DCA_{xy},                 pt_DCA,   isyst *******************************************
            Int_t    bins_DCA_syst[3] = {nBins_dca_xy,           nBins_pt_DCA,    50 };
            Double_t xmin_DCA_syst[3] = {        -1.0,              pt_DCA[0],     0 };
            Double_t xmax_DCA_syst[3] = {        +1.0,   pt_DCA[nBins_pt_DCA],    50 };
            //***************************************************************************************************************************************

            hDCAxy_antiprotons_prim_Syst = new THnSparseF ("hDCAxy_antiprotons_prim_Syst","",3, bins_DCA_syst, xmin_DCA_syst, xmax_DCA_syst);
            hDCAxy_antiprotons_prim_Syst -> GetAxis(0) -> Set(nBins_dca_xy, dca_xy);
            hDCAxy_antiprotons_prim_Syst -> GetAxis(1) -> Set(nBins_pt_DCA, pt_DCA);
            hDCAxy_antiprotons_prim_Syst -> Sumw2();

            hDCAxy_antiprotons_sec_Syst = new THnSparseF ("hDCAxy_antiprotons_sec_Syst","",3, bins_DCA_syst, xmin_DCA_syst, xmax_DCA_syst);
            hDCAxy_antiprotons_sec_Syst -> GetAxis(0) -> Set(nBins_dca_xy, dca_xy);
            hDCAxy_antiprotons_sec_Syst -> GetAxis(1) -> Set(nBins_pt_DCA, pt_DCA);
            hDCAxy_antiprotons_sec_Syst -> Sumw2();

            fOutputList -> Add(hDCAxy_antiprotons_prim_Syst);
            fOutputList -> Add(hDCAxy_antiprotons_sec_Syst);
        }
    }

    if(fIsPion) {
        if ((!fIsMC)&&(fIsITSrecalib)) {
            
            //***************************  Nch_Transverse,  nsigma_{TOF},                 pt_TOF ***************************************************
            Int_t    bins_TOF[3] = {         nBins_Nch_Tr,           200,           nBins_pt_TOF };
            Double_t xmin_TOF[3] = {            Nch_Tr[0],         -10.0,              pt_TOF[0] };
            Double_t xmax_TOF[3] = { Nch_Tr[nBins_Nch_Tr],         +10.0,   pt_TOF[nBins_pt_TOF] };
            //***************************************************************************************************************************************

            hnsigmaTOF_pos_pions_Toward     = new THnSparseF ("hnsigmaTOF_pos_pions_Toward","",3, bins_TOF, xmin_TOF, xmax_TOF);
            hnsigmaTOF_pos_pions_Away       = new THnSparseF ("hnsigmaTOF_pos_pions_Away","",3, bins_TOF, xmin_TOF, xmax_TOF);
            hnsigmaTOF_pos_pions_Transverse = new THnSparseF ("hnsigmaTOF_pos_pions_Transverse","",3, bins_TOF, xmin_TOF, xmax_TOF);
            hnsigmaTOF_neg_pions_Toward     = new THnSparseF ("hnsigmaTOF_neg_pions_Toward","",3, bins_TOF, xmin_TOF, xmax_TOF);
            hnsigmaTOF_neg_pions_Away       = new THnSparseF ("hnsigmaTOF_neg_pions_Away","",3, bins_TOF, xmin_TOF, xmax_TOF);
            hnsigmaTOF_neg_pions_Transverse = new THnSparseF ("hnsigmaTOF_neg_pions_Transverse","",3, bins_TOF, xmin_TOF, xmax_TOF);
            hnsigmaTOF_pos_pions_Toward     -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hnsigmaTOF_pos_pions_Away       -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hnsigmaTOF_pos_pions_Transverse -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hnsigmaTOF_neg_pions_Toward     -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hnsigmaTOF_neg_pions_Away       -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hnsigmaTOF_neg_pions_Transverse -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hnsigmaTOF_pos_pions_Toward     -> GetAxis(2) -> Set(nBins_pt_TOF, pt_TOF);
            hnsigmaTOF_pos_pions_Away       -> GetAxis(2) -> Set(nBins_pt_TOF, pt_TOF);
            hnsigmaTOF_pos_pions_Transverse -> GetAxis(2) -> Set(nBins_pt_TOF, pt_TOF);
            hnsigmaTOF_neg_pions_Toward     -> GetAxis(2) -> Set(nBins_pt_TOF, pt_TOF);
            hnsigmaTOF_neg_pions_Away       -> GetAxis(2) -> Set(nBins_pt_TOF, pt_TOF);
            hnsigmaTOF_neg_pions_Transverse -> GetAxis(2) -> Set(nBins_pt_TOF, pt_TOF);
            hnsigmaTOF_pos_pions_Toward     -> Sumw2();
            hnsigmaTOF_pos_pions_Away       -> Sumw2();
            hnsigmaTOF_pos_pions_Transverse -> Sumw2();
            hnsigmaTOF_neg_pions_Toward     -> Sumw2();
            hnsigmaTOF_neg_pions_Away       -> Sumw2();
            hnsigmaTOF_neg_pions_Transverse -> Sumw2();
            fOutputList -> Add (hnsigmaTOF_pos_pions_Toward);
            fOutputList -> Add (hnsigmaTOF_pos_pions_Away);
            fOutputList -> Add (hnsigmaTOF_pos_pions_Transverse);
            fOutputList -> Add (hnsigmaTOF_neg_pions_Toward);
            fOutputList -> Add (hnsigmaTOF_neg_pions_Away);
            fOutputList -> Add (hnsigmaTOF_neg_pions_Transverse);

            //***************************  Nch_Transverse,      DCA_{xy},                 pt_DCA ***************************************************
            Int_t    bins_DCA[3] = {         nBins_Nch_Tr,  nBins_dca_xy,           nBins_pt_DCA };
            Double_t xmin_DCA[3] = {            Nch_Tr[0],          -1.0,              pt_DCA[0] };
            Double_t xmax_DCA[3] = { Nch_Tr[nBins_Nch_Tr],          +1.0,   pt_DCA[nBins_pt_DCA] };
            //***************************************************************************************************************************************

            hDCAxy_pos_pions_Toward     = new THnSparseF ("hDCAxy_pos_pions_Toward","",3, bins_DCA, xmin_DCA, xmax_DCA);
            hDCAxy_pos_pions_Away       = new THnSparseF ("hDCAxy_pos_pions_Away","",3, bins_DCA, xmin_DCA, xmax_DCA);
            hDCAxy_pos_pions_Transverse = new THnSparseF ("hDCAxy_pos_pions_Transverse","",3, bins_DCA, xmin_DCA, xmax_DCA);
            hDCAxy_neg_pions_Toward     = new THnSparseF ("hDCAxy_neg_pions_Toward","",3, bins_DCA, xmin_DCA, xmax_DCA);
            hDCAxy_neg_pions_Away       = new THnSparseF ("hDCAxy_neg_pions_Away","",3, bins_DCA, xmin_DCA, xmax_DCA);
            hDCAxy_neg_pions_Transverse = new THnSparseF ("hDCAxy_neg_pions_Transverse","",3, bins_DCA, xmin_DCA, xmax_DCA);
            hDCAxy_pos_pions_Toward     -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hDCAxy_pos_pions_Away       -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hDCAxy_pos_pions_Transverse -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hDCAxy_neg_pions_Toward     -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hDCAxy_neg_pions_Away       -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hDCAxy_neg_pions_Transverse -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hDCAxy_pos_pions_Toward     -> GetAxis(1) -> Set(nBins_dca_xy, dca_xy);
            hDCAxy_pos_pions_Away       -> GetAxis(1) -> Set(nBins_dca_xy, dca_xy);
            hDCAxy_pos_pions_Transverse -> GetAxis(1) -> Set(nBins_dca_xy, dca_xy);
            hDCAxy_neg_pions_Toward     -> GetAxis(1) -> Set(nBins_dca_xy, dca_xy);
            hDCAxy_neg_pions_Away       -> GetAxis(1) -> Set(nBins_dca_xy, dca_xy);
            hDCAxy_neg_pions_Transverse -> GetAxis(1) -> Set(nBins_dca_xy, dca_xy);
            hDCAxy_pos_pions_Toward     -> GetAxis(2) -> Set(nBins_pt_DCA, pt_DCA);
            hDCAxy_pos_pions_Away       -> GetAxis(2) -> Set(nBins_pt_DCA, pt_DCA);
            hDCAxy_pos_pions_Transverse -> GetAxis(2) -> Set(nBins_pt_DCA, pt_DCA);
            hDCAxy_neg_pions_Toward     -> GetAxis(2) -> Set(nBins_pt_DCA, pt_DCA);
            hDCAxy_neg_pions_Away       -> GetAxis(2) -> Set(nBins_pt_DCA, pt_DCA);
            hDCAxy_neg_pions_Transverse -> GetAxis(2) -> Set(nBins_pt_DCA, pt_DCA);
            hDCAxy_pos_pions_Toward     -> Sumw2();
            hDCAxy_pos_pions_Away       -> Sumw2();
            hDCAxy_pos_pions_Transverse -> Sumw2();
            hDCAxy_neg_pions_Toward     -> Sumw2();
            hDCAxy_neg_pions_Away       -> Sumw2();
            hDCAxy_neg_pions_Transverse -> Sumw2();
            fOutputList -> Add(hDCAxy_pos_pions_Toward);
            fOutputList -> Add(hDCAxy_pos_pions_Away);
            fOutputList -> Add(hDCAxy_pos_pions_Transverse);
            fOutputList -> Add(hDCAxy_neg_pions_Toward);
            fOutputList -> Add(hDCAxy_neg_pions_Away);
            fOutputList -> Add(hDCAxy_neg_pions_Transverse);

            //********************************  Nch_Transverse,  nsigma_{TOF},                  pt_TOF, isyst ****************************************
            Int_t    bins_TOF_syst[4] = {         nBins_Nch_Tr,           200,            nBins_pt_TOF,    50 };
            Double_t xmin_TOF_syst[4] = {            Nch_Tr[0],          -10.0,              pt_TOF[0],     0 };
            Double_t xmax_TOF_syst[4] = { Nch_Tr[nBins_Nch_Tr],          +10.0,   pt_TOF[nBins_pt_TOF],    50 };
            //***************************************************************************************************************************************

            hnsigmaTOF_pos_pions_Syst = new THnSparseF ("hnsigmaTOF_pos_pions_Syst","",4, bins_TOF_syst, xmin_TOF_syst, xmax_TOF_syst);
            hnsigmaTOF_neg_pions_Syst = new THnSparseF ("hnsigmaTOF_neg_pions_Syst","",4, bins_TOF_syst, xmin_TOF_syst, xmax_TOF_syst);
            hnsigmaTOF_pos_pions_Syst -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hnsigmaTOF_neg_pions_Syst -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hnsigmaTOF_pos_pions_Syst -> GetAxis(2) -> Set(nBins_pt_TOF, pt_TOF);
            hnsigmaTOF_neg_pions_Syst -> GetAxis(2) -> Set(nBins_pt_TOF, pt_TOF);
            hnsigmaTOF_pos_pions_Syst -> Sumw2();
            hnsigmaTOF_neg_pions_Syst -> Sumw2();
            fOutputList -> Add (hnsigmaTOF_pos_pions_Syst);
            fOutputList -> Add (hnsigmaTOF_neg_pions_Syst);


            //***************************       Nch_Transverse,      DCA_{xy},                 pt_DCA,   isyst *******************************************
            Int_t    bins_DCA_syst[4] = {         nBins_Nch_Tr,  nBins_dca_xy,           nBins_pt_DCA,    50 };
            Double_t xmin_DCA_syst[4] = {            Nch_Tr[0],          -1.0,              pt_DCA[0],     0 };
            Double_t xmax_DCA_syst[4] = { Nch_Tr[nBins_Nch_Tr],          +1.0,   pt_DCA[nBins_pt_DCA],    50 };
            //***************************************************************************************************************************************

            hDCAxy_pos_pions_Syst = new THnSparseF ("hDCAxy_pos_pions_Syst","",4, bins_DCA_syst, xmin_DCA_syst, xmax_DCA_syst);
            hDCAxy_neg_pions_Syst = new THnSparseF ("hDCAxy_neg_pions_Syst","",4, bins_DCA_syst, xmin_DCA_syst, xmax_DCA_syst);
            hDCAxy_pos_pions_Syst -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hDCAxy_neg_pions_Syst -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
            hDCAxy_neg_pions_Syst -> GetAxis(1) -> Set(nBins_dca_xy, dca_xy);
            hDCAxy_pos_pions_Syst -> GetAxis(1) -> Set(nBins_dca_xy, dca_xy);
            hDCAxy_pos_pions_Syst -> GetAxis(2) -> Set(nBins_pt_DCA, pt_DCA);
            hDCAxy_neg_pions_Syst -> GetAxis(2) -> Set(nBins_pt_DCA, pt_DCA);
            hDCAxy_pos_pions_Syst -> Sumw2();
            hDCAxy_neg_pions_Syst -> Sumw2();
            fOutputList -> Add(hDCAxy_pos_pions_Syst);
            fOutputList -> Add(hDCAxy_neg_pions_Syst);
        }

        if (fIsMC && fIsITSrecalib)  {
            h_pos_pions_Gen          = new TH1F ("h_pos_pions_Gen","",940,0.3,5.0);
            h_neg_pions_Gen          = new TH1F ("h_neg_pions_Gen","",940,0.3,5.0);
            hnsigmaTOF_pos_pions_Rec = new TH2F ("hnsigmaTOF_pos_pions_Rec","",940,0.3,5.0,100,-5,5);
            hnsigmaTOF_neg_pions_Rec = new TH2F ("hnsigmaTOF_neg_pions_Rec","",940,0.3,5.0,100,-5,5);
            h_pos_pions_Gen          -> Sumw2();
            h_neg_pions_Gen          -> Sumw2();
            hnsigmaTOF_pos_pions_Rec -> Sumw2();
            hnsigmaTOF_neg_pions_Rec -> Sumw2();
            fOutputList -> Add (h_pos_pions_Gen);
            fOutputList -> Add (h_neg_pions_Gen);
            fOutputList -> Add (hnsigmaTOF_pos_pions_Rec);
            fOutputList -> Add (hnsigmaTOF_neg_pions_Rec);

            hDCAxy_pos_pions_prim = new TH2F("hDCAxy_pos_pions_prim", "", nBins_pt_DCA, pt_DCA, nBins_dca_xy, dca_xy);   
            hDCAxy_neg_pions_prim = new TH2F("hDCAxy_neg_pions_prim", "", nBins_pt_DCA, pt_DCA, nBins_dca_xy, dca_xy);   
            hDCAxy_pos_pions_sec  = new TH2F("hDCAxy_pos_pions_sec", "", nBins_pt_DCA, pt_DCA, nBins_dca_xy, dca_xy);
            hDCAxy_neg_pions_sec  = new TH2F("hDCAxy_neg_pions_sec", "", nBins_pt_DCA, pt_DCA, nBins_dca_xy, dca_xy);
            hDCAxy_pos_pions_prim -> Sumw2();
            hDCAxy_neg_pions_prim -> Sumw2();
            hDCAxy_pos_pions_sec  -> Sumw2();
            hDCAxy_neg_pions_sec  -> Sumw2();
            fOutputList -> Add(hDCAxy_pos_pions_prim);
            fOutputList -> Add(hDCAxy_neg_pions_prim);
            fOutputList -> Add(hDCAxy_pos_pions_sec);
            fOutputList -> Add(hDCAxy_neg_pions_sec);



            //Histograms for Syst. Uncertainties
            hnsigmaTOF_pos_pions_Rec_Syst = new TH2F ("hnsigmaTOF_pos_pions_Rec_Syst","",940,0.3,5.0,50,0,50);
            hnsigmaTOF_neg_pions_Rec_Syst = new TH2F ("hnsigmaTOF_neg_pions_Rec_Syst","",940,0.3,5.0,50,0,50);
            hnsigmaTOF_pos_pions_Rec_Syst -> Sumw2();
            hnsigmaTOF_neg_pions_Rec_Syst -> Sumw2();
            fOutputList -> Add(hnsigmaTOF_pos_pions_Rec_Syst);
            fOutputList -> Add(hnsigmaTOF_neg_pions_Rec_Syst);

            //***************************    DCA_{xy},                 pt_DCA,   isyst *******************************************
            Int_t    bins_DCA_syst[3] = {nBins_dca_xy,           nBins_pt_DCA,    50 };
            Double_t xmin_DCA_syst[3] = {        -1.0,              pt_DCA[0],     0 };
            Double_t xmax_DCA_syst[3] = {        +1.0,   pt_DCA[nBins_pt_DCA],    50 };
            //***************************************************************************************************************************************

            hDCAxy_pos_pions_prim_Syst = new THnSparseF ("hDCAxy_pos_pions_prim_Syst","",3, bins_DCA_syst, xmin_DCA_syst, xmax_DCA_syst);
            hDCAxy_neg_pions_prim_Syst = new THnSparseF ("hDCAxy_pos_neg_prim_Syst","",3, bins_DCA_syst, xmin_DCA_syst, xmax_DCA_syst);
            hDCAxy_pos_pions_prim_Syst -> GetAxis(0) -> Set(nBins_dca_xy, dca_xy);
            hDCAxy_neg_pions_prim_Syst -> GetAxis(0) -> Set(nBins_dca_xy, dca_xy);
            hDCAxy_pos_pions_prim_Syst -> GetAxis(1) -> Set(nBins_pt_DCA, pt_DCA);
            hDCAxy_neg_pions_prim_Syst -> GetAxis(1) -> Set(nBins_pt_DCA, pt_DCA);
            hDCAxy_pos_pions_prim_Syst -> Sumw2();
            hDCAxy_neg_pions_prim_Syst -> Sumw2();

            hDCAxy_pos_pions_sec_Syst = new THnSparseF ("hDCAxy_pos_pions_sec_Syst","",3, bins_DCA_syst, xmin_DCA_syst, xmax_DCA_syst);
            hDCAxy_neg_pions_sec_Syst = new THnSparseF ("hDCAxy_neg_pions_sec_Syst","",3, bins_DCA_syst, xmin_DCA_syst, xmax_DCA_syst);
            hDCAxy_pos_pions_sec_Syst -> GetAxis(0) -> Set(nBins_dca_xy, dca_xy);
            hDCAxy_neg_pions_sec_Syst -> GetAxis(0) -> Set(nBins_dca_xy, dca_xy);
            hDCAxy_pos_pions_sec_Syst -> GetAxis(1) -> Set(nBins_pt_DCA, pt_DCA);
            hDCAxy_neg_pions_sec_Syst -> GetAxis(1) -> Set(nBins_pt_DCA, pt_DCA);
            hDCAxy_pos_pions_sec_Syst -> Sumw2();
            hDCAxy_neg_pions_sec_Syst -> Sumw2();

            fOutputList -> Add(hDCAxy_pos_pions_prim_Syst);
            fOutputList -> Add(hDCAxy_pos_pions_sec_Syst);
            fOutputList -> Add(hDCAxy_neg_pions_prim_Syst);
            fOutputList -> Add(hDCAxy_neg_pions_sec_Syst);
        }
    }

    //Track Selection Cuts: AntiProtons
    fESDtrackCuts_AntiProton = new AliESDtrackCuts("fESDtrackCuts_AntiProton");
    fESDtrackCuts_AntiProton -> SetAcceptKinkDaughters(kFALSE);
    fESDtrackCuts_AntiProton -> SetRequireTPCRefit(kTRUE);
    fESDtrackCuts_AntiProton -> SetRequireITSRefit(kTRUE);
    fESDtrackCuts_AntiProton -> SetMinNCrossedRowsTPC(50);
    fESDtrackCuts_AntiProton -> SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fESDtrackCuts_AntiProton -> SetMinNClustersITS(1);
    fESDtrackCuts_AntiProton -> SetClusterRequirementITS (AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
    fESDtrackCuts_AntiProton -> SetMaxChi2PerClusterITS(36);
    fESDtrackCuts_AntiProton -> SetMaxChi2PerClusterTPC(10);
    fESDtrackCuts_AntiProton -> SetEtaRange(-2.0,2.0);
    fESDtrackCuts_AntiProton -> SetMaxDCAToVertexXY(2);
    fESDtrackCuts_AntiProton -> SetMaxDCAToVertexZ(2);
    fESDtrackCuts_AntiProton -> SetDCAToVertex2D(kFALSE);
    fESDtrackCuts_AntiProton -> SetRequireSigmaToVertex(kFALSE);

    //Track Selection Cuts: Pions
    fESDtrackCuts_Pions = new AliESDtrackCuts("fESDtrackCuts_Pions");
    fESDtrackCuts_Pions -> SetAcceptKinkDaughters(kFALSE);
    fESDtrackCuts_Pions -> SetRequireTPCRefit(kTRUE);
    fESDtrackCuts_Pions -> SetRequireITSRefit(kTRUE);
    fESDtrackCuts_Pions -> SetMinNCrossedRowsTPC(50);
    fESDtrackCuts_Pions -> SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
    fESDtrackCuts_Pions -> SetMinNClustersITS(1);
    fESDtrackCuts_Pions -> SetClusterRequirementITS (AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
    fESDtrackCuts_Pions -> SetMaxChi2PerClusterITS(36);
    fESDtrackCuts_Pions -> SetMaxChi2PerClusterTPC(10);
    fESDtrackCuts_Pions -> SetEtaRange(-2.0,2.0);
    fESDtrackCuts_Pions -> SetMaxDCAToVertexXY(2);
    fESDtrackCuts_Pions -> SetMaxDCAToVertexZ(2);
    fESDtrackCuts_Pions -> SetDCAToVertex2D(kFALSE);
    fESDtrackCuts_Pions -> SetRequireSigmaToVertex(kFALSE);


    //Track Selection Cuts: Charged Tracks
    if (!fIsMC)  {
        
        fESDtrackCuts_LeadingTrack = new AliESDtrackCuts("fESDtrackCuts_LeadingTrack");
        fESDtrackCuts_LeadingTrack -> SetAcceptKinkDaughters(kFALSE);
        fESDtrackCuts_LeadingTrack -> SetRequireTPCRefit(kTRUE);
        fESDtrackCuts_LeadingTrack -> SetRequireITSRefit(kTRUE);
        fESDtrackCuts_LeadingTrack -> SetMinNCrossedRowsTPC(70);
        fESDtrackCuts_LeadingTrack -> SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
        fESDtrackCuts_LeadingTrack -> SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kAny);
        fESDtrackCuts_LeadingTrack -> SetMaxChi2PerClusterTPC(4);
        fESDtrackCuts_LeadingTrack -> SetMaxChi2PerClusterITS(36);
        fESDtrackCuts_LeadingTrack -> SetEtaRange(-0.8,0.8);
        fESDtrackCuts_LeadingTrack -> SetPtRange(0.15,200);
        fESDtrackCuts_LeadingTrack -> SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
        fESDtrackCuts_LeadingTrack -> SetMaxDCAToVertexZ(2);
        fESDtrackCuts_LeadingTrack -> SetDCAToVertex2D(kFALSE);
        fESDtrackCuts_LeadingTrack -> SetRequireSigmaToVertex(kFALSE);
        
        fESDtrackCuts_TransverseMult = new AliESDtrackCuts ("fESDtrackCuts_TransverseMult");
        fESDtrackCuts_TransverseMult -> SetAcceptKinkDaughters(kFALSE);
        fESDtrackCuts_TransverseMult -> SetRequireTPCRefit(kTRUE);
        fESDtrackCuts_TransverseMult -> SetRequireITSRefit(kTRUE);
        fESDtrackCuts_TransverseMult -> SetMinNClustersTPC(50);
        fESDtrackCuts_TransverseMult -> SetMaxChi2PerClusterTPC(4);
        fESDtrackCuts_TransverseMult -> SetMaxDCAToVertexZ(3.2);
        fESDtrackCuts_TransverseMult -> SetMaxDCAToVertexXY(2.4);
        fESDtrackCuts_TransverseMult -> SetDCAToVertex2D(kTRUE);
        fESDtrackCuts_TransverseMult -> SetEtaRange(-0.8,0.8);
        fESDtrackCuts_TransverseMult -> SetPtRange(0.15,200);
    }

    PostData(1, fOutputList);
    PostData(2, fQAList);
}
//____________________________________________________________________________________________________________________________
void AliAnalysisTaskAntiProtons_vs_RT_pPb::UserExec(Option_t *)  {
    
    //Get Input Event
    if ( !GetESDEvent ()) return;
    if (fIsMC && (!GetMCEvent ())) return;
    
    //Load PID Response
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    
    //Process Real or Simulated Event
    if (!fIsMC)    ProcessRealEvent ();
    if (fIsMC)     ProcessSimEvent ();
    
    //Post Output Data
    PostData(1, fOutputList);
}
//____________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskAntiProtons_vs_RT_pPb::GetESDEvent ()  {
    
   //Get Input Event
   fESDevent = dynamic_cast <AliESDEvent*>(InputEvent());
   if (!fESDevent) return kFALSE;
   hNumberOfEvents -> Fill(0.5);
   
       
   //Standard Event Cuts
   if (!fESDEventSelection.AcceptEvent(fESDevent)) {
       PostData(2, fQAList);
       return kFALSE;
   }
   hNumberOfEvents -> Fill(1.5);

    
   //Reject Events with Incomplete DAQ
   if (fESDevent->IsIncompleteDAQ()) return kFALSE;
   hNumberOfEvents -> Fill(2.5);
   
   
   //V0 Timing Decision
   AliVVZERO *vzeroData = fESDevent->GetVZEROData();
   if (!(vzeroData->GetV0ADecision()) || !(vzeroData->GetV0CDecision())) return kFALSE;
   hNumberOfEvents -> Fill(3.5);
   
   
   //Pile-up Rejection
   Int_t nClustersLayer0 = fESDevent->GetNumberOfITSClusters(0);
   Int_t nClustersLayer1 = fESDevent->GetNumberOfITSClusters(1);
   Int_t nTracklets      = fESDevent->GetMultiplicity()->GetNumberOfTracklets();
   if ((nClustersLayer0 + nClustersLayer1) > 65.0 + (Double_t)nTracklets*4.0) return kFALSE;
   hNumberOfEvents -> Fill(4.5);
   

   //Primary Vertex Tracks
   AliESDVertex *vertex_tracks = (AliESDVertex*) fESDevent->GetPrimaryVertexTracks();
   if (!vertex_tracks) return kFALSE;
   hNumberOfEvents -> Fill(5.5);
   
    
   //Vertex Contributors Tracks
   if ( vertex_tracks->GetNContributors() < 1 ) return kFALSE;
   hNumberOfEvents -> Fill(6.5);
   
    
   //Primary Vertex SPD
   AliESDVertex *vertex_SPD = (AliESDVertex*) fESDevent->GetPrimaryVertexSPD();
   if (!vertex_SPD) return kFALSE;
   hNumberOfEvents -> Fill(7.5);
   
    
   //Vertex Contributors SPD
   if ( vertex_SPD->GetNContributors() < 1 ) return kFALSE;
   hNumberOfEvents -> Fill(8.5);
   
    
   //SPD Pile-up in Mult Bins
   if (fESDevent->IsPileupFromSPDInMultBins()) return kFALSE;
   hNumberOfEvents -> Fill(9.5);
   
    
   //Cut on Z-Vertex Resolution
   if (TMath::Abs(vertex_SPD->GetZ() - vertex_tracks->GetZ()) > 0.3) return kFALSE;
   hNumberOfEvents -> Fill(10.5);

    
   //Primary Vertex Selection
   if ( vertex_tracks->GetZ() < -10.0 ) return kFALSE;
   if ( vertex_tracks->GetZ() > +10.0 ) return kFALSE;
   hNumberOfEvents -> Fill(11.5);
          
    
   //Multiplicity
   AliMultSelection *multiplicitySelection = (AliMultSelection*) fESDevent->FindListObject("MultSelection");
   if( !multiplicitySelection) return kFALSE;
   hNumberOfEvents -> Fill(12.5);
   
    
   //Multiplicity Distribution
   Double_t mult_percentile = multiplicitySelection->GetMultiplicityPercentile("V0M");
   hMultPercentile -> Fill(mult_percentile);
   
    
   //Selection of Multiplicity Range
   if (mult_percentile <  fMultMin) return kFALSE;
   if (mult_percentile >  fMultMax) return kFALSE;
   hNumberOfEvents -> Fill(13.5);
   
     
    return kTRUE;
}
//__________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskAntiProtons_vs_RT_pPb::GetMCEvent ()  {
    
    //Get MC Event
    fMCEvent = MCEvent();
    if (!fMCEvent) return kFALSE;
    hNumberOfEvents -> Fill(14.5);
   
    //MC Event Handler
    fMCEventHandler = dynamic_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
    if(!fMCEventHandler) return kFALSE;
    hNumberOfEvents -> Fill(15.5);

    return kTRUE;
}
//__________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskAntiProtons_vs_RT_pPb::ProcessRealEvent ()  {
    
    //Leading-Track ID
    Int_t leading_track_ID = GetLeadingTrack ();
    
    //Leading Track Pt
    AliESDtrack *track = (AliESDtrack*) fESDevent->GetTrack(leading_track_ID);
    Double_t pt_leading = track->Pt();
    
    //Rejection Factors
    hEventsWithLeadingTrack -> Fill (0.5);
    if (pt_leading>2.0) hEventsWithLeadingTrack -> Fill (1.5);
    if (pt_leading>2.5) hEventsWithLeadingTrack -> Fill (2.5);
    if (pt_leading>3.0) hEventsWithLeadingTrack -> Fill (3.5);
    if (pt_leading>3.5) hEventsWithLeadingTrack -> Fill (4.5);
    if (pt_leading>4.0) hEventsWithLeadingTrack -> Fill (5.5);
    if (pt_leading>4.5) hEventsWithLeadingTrack -> Fill (6.5);
    if (pt_leading>5.0) hEventsWithLeadingTrack -> Fill (7.5);
    if (pt_leading>5.5) hEventsWithLeadingTrack -> Fill (8.5);
    if (pt_leading>6.0) hEventsWithLeadingTrack -> Fill (9.5);
    
    //Min Pt of Leading Track
    if (pt_leading<fPt_min_leading) return;
    
    //Multiplicities in Azimuthal Regions
    Int_t mult_Transverse(0);
    Int_t mult_Toward(0);
    Int_t mult_Away(0);
    
    //AntiProtons Candidates ID
    vector<Int_t> candidate_ID;
    
    //Loop over Reconstructed Tracks
    for (Int_t i=0 ; i<fESDevent->GetNumberOfTracks() ; i++)  {

        //Get Reconstructed Track
        AliESDtrack *track = (AliESDtrack*) fESDevent->GetTrack(i);
        if (!track) continue;
        
        //Store AntiProton Candidates
        if (IsTrackCandidate (track))  { candidate_ID.push_back(i); }
        
        //Track Cuts & Multiplicities in Azimuthal Regions
        if (!PassedTrackQualityCuts_TransverseMult (track)) continue;
        if (IsTrackInTransverseRegion (track,leading_track_ID)) mult_Transverse++;
        if (IsTrackInTowardRegion     (track,leading_track_ID)) mult_Toward++;
        if (IsTrackInAwayRegion       (track,leading_track_ID)) mult_Away++;
    }
    
    //Fill Histograms
    hMultTransverse -> Fill (mult_Transverse);
    hMultToward     -> Fill (mult_Toward);
    hMultAway       -> Fill (mult_Away);
    
    //Fill R_{T} Distribution
    Double_t Rt = static_cast<Double_t>(mult_Transverse)/fAverage_Nch_Transv;
    hRtDistribution -> Fill (Rt);
    
    //Correlations between Transverse and Integrated Mult
    AliMultSelection *multiplicitySelection = (AliMultSelection*) fESDevent->FindListObject("MultSelection");
    Double_t mult_percentile = multiplicitySelection->GetMultiplicityPercentile("V0M");

    hNchTransv_NchTot         -> Fill (mult_Transverse,(mult_Transverse+mult_Toward+mult_Away));
    hRt_NchTot                -> Fill (Rt,(mult_Transverse+mult_Toward+mult_Away));
    hNchTransv_MultPercentile -> Fill (mult_Transverse,mult_percentile);
    hRt_MultPercentile        -> Fill (Rt,mult_percentile);
    
    //Number of AntiProtons Candidates
    Int_t nCandidate = (Int_t)candidate_ID.size();
    if(fIsPion) {
        hNumberOfPions ->Fill(nCandidate);
    }
    if(!fIsPion){
        hNumberOfAntiProtons -> Fill (nCandidate);
    }
    
    //Loop over AntiProtons Candidates
    for (Int_t i=0 ; i<nCandidate ; i++)  {

        //Get Reconstructed Track
        AliESDtrack *track = (AliESDtrack*) fESDevent->GetTrack(candidate_ID[i]);
        if (!track) continue;
    
        //Fill Histograms: Standard Cuts
        FillHistograms_StandardCuts (mult_Transverse,leading_track_ID,track);
        
        //Fill Histograms: Systematic Uncertainties
        if(fIsITSrecalib) 
        {
            for (Int_t isyst=0 ; isyst<50 ; isyst++) {FillHistograms_Systematics (mult_Transverse,leading_track_ID,track,isyst);}
        }
    }
}
//__________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskAntiProtons_vs_RT_pPb::ProcessSimEvent ()  {
    
    //Loop over Generated Particles
    if(fIsITSrecalib)
    {
        for (Int_t i=0; i<fMCEvent->GetNumberOfTracks(); i++)  {
           
            //MC Particle Selection
            AliMCParticle *particle = (AliMCParticle*) fMCEvent->GetTrack(i);
            if (!particle) continue;
            if (!particle->IsPhysicalPrimary()) continue;
            if(!fIsPion)
            {
                if ( particle->PdgCode() != -2212 ) continue;
                if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i,fMCEvent)) continue;

                //Rapidity Selection
                Double_t m = AliPID::ParticleMass(AliPID::kProton);
                Double_t E = TMath::Sqrt(m*m + particle->P()*particle->P());
                TLorentzVector P (particle->Px(),particle->Py(),particle->Pz(),E);
                Double_t y_lab = P.Rapidity();
                Double_t y_cms = y_lab-0.465;
                if (y_cms<-1.0) continue;
                if (y_cms> 0.0) continue;

                //Fill Generated p_{T} Spectra
                h_antiprotons_Gen -> Fill(particle->Pt());
            }

            if(fIsPion)
            {
                if ( TMath::Abs(particle->PdgCode())!= 211) continue;
                if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i,fMCEvent)) continue;

                //Rapidity Selection
                Double_t m = AliPID::ParticleMass(AliPID::kPion);
                Double_t E = TMath::Sqrt(m*m + particle->P()*particle->P());
                TLorentzVector P (particle->Px(),particle->Py(),particle->Pz(),E);
                Double_t y_lab = P.Rapidity();
                Double_t y_cms = y_lab-0.465;
                Double_t charge = particle->Charge();
                if (y_cms<-1.0) continue;
                if (y_cms> 0.0) continue;
                //Fill Generated p_{T} Spectra
                if (charge > 0) h_pos_pions_Gen -> Fill(particle->Pt());
                if (charge < 0) h_neg_pions_Gen -> Fill(particle->Pt());
            }   
        }
    }
       
    //Loop over Reconstructed Tracks
    
    for (Int_t i=0 ; i<fESDevent->GetNumberOfTracks() ; i++)  {

        //Get Reconstructed Track
        AliESDtrack *track = (AliESDtrack*) fESDevent->GetTrack(i);
        if (!track) continue;
            
        //MC Particle
        AliMCParticle *particle = (AliMCParticle*) fMCEvent->GetTrack(TMath::Abs(track->GetLabel()));
        if (!particle) continue;
        if(!fIsPion)
        {
            if ( particle->PdgCode() != -2212 ) continue;
            
            //Fill Histograms: Standard Cuts - option for both fIsITSrecalib true and false inside the function
            FillHistograms_StandardCuts_Sim (track);
                
            //Fill Histograms: Systematic Uncertainties - only for fIsITSrecalib = true
            if(fIsITSrecalib)
            {
                for (Int_t isyst=0 ; isyst<50 ; isyst++) {FillHistograms_Systematics_Sim (track,isyst);}
            }
        }
        if(fIsPion)
        {
            //cout << "PDG before " <<  particle->PdgCode() << endl;
            if ( TMath::Abs(particle->PdgCode()) != 211 ) continue;
            //cout << "PDG after " <<  particle->PdgCode() << endl;

            //Fill Histograms: Standard Cuts - option for both fIsITSrecalib true and false inside the function
            FillHistograms_StandardCuts_Sim (track);
                
            //Fill Histograms: Systematic Uncertainties - only for fIsITSrecalib = true
            if(fIsITSrecalib)
            {
                for (Int_t isyst=0 ; isyst<50 ; isyst++) {FillHistograms_Systematics_Sim (track,isyst);}
            }
        
        }
    } 
}
//__________________________________________________________________________________________________________________________________________________
Int_t AliAnalysisTaskAntiProtons_vs_RT_pPb::GetLeadingTrack ()  {
    
    //Initialization
    Int_t ID_leading_track(0);
    Double_t pt_max(0);
    
    //Loop over Reconstructed Tracks
    for (Int_t i=0 ; i<fESDevent->GetNumberOfTracks() ; i++)  {

        //Get Reconstructed Track
        AliESDtrack *track = (AliESDtrack*) fESDevent->GetTrack(i);
        if (!track) continue;
        if (!PassedTrackQualityCuts_LeadingTrack (track)) continue;
        if (track->Pt() > pt_max)  {
            
            pt_max = track->Pt();
            ID_leading_track = i;
        }
    }
    
    return ID_leading_track;
}
//__________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskAntiProtons_vs_RT_pPb::IsTrackInTransverseRegion (AliESDtrack *track, Int_t leading_track_ID)  {
    
    //Initialization
    Bool_t isInTransverseRegion = kFALSE;
    
    //Get Leading Track
    AliESDtrack *leading_track = (AliESDtrack*) fESDevent->GetTrack(leading_track_ID);
    
    //DeltaPhi
    Double_t phi_ref   = TVector2::Phi_0_2pi (leading_track->Phi());
    Double_t phi_trk   = TVector2::Phi_0_2pi (track->Phi());
    Double_t delta_phi = (180.0/TMath::Pi())*TVector2::Phi_0_2pi (phi_trk-phi_ref);
    
    if (delta_phi>=60.0  && delta_phi<120.0) isInTransverseRegion=kTRUE;
    if (delta_phi>=240.0 && delta_phi<300.0) isInTransverseRegion=kTRUE;

    return isInTransverseRegion;
}
//__________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskAntiProtons_vs_RT_pPb::IsTrackInTowardRegion (AliESDtrack *track, Int_t leading_track_ID)  {
    
    //Initialization
    Bool_t isInTowardRegion = kFALSE;
    
    //Get Leading Track
    AliESDtrack *leading_track = (AliESDtrack*) fESDevent->GetTrack(leading_track_ID);
    
    //DeltaPhi
    Double_t phi_ref   = TVector2::Phi_0_2pi (leading_track->Phi());
    Double_t phi_trk   = TVector2::Phi_0_2pi (track->Phi());
    Double_t delta_phi = (180.0/TMath::Pi())*TVector2::Phi_0_2pi (phi_trk-phi_ref);

    if (delta_phi>=0.0   && delta_phi<60.0)   isInTowardRegion=kTRUE;
    if (delta_phi>=300.0 && delta_phi<=360.0) isInTowardRegion=kTRUE;
    
    return isInTowardRegion;
}
//__________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskAntiProtons_vs_RT_pPb::IsTrackInAwayRegion (AliESDtrack *track, Int_t leading_track_ID)  {
    
    //Initialization
    Bool_t isInAwayRegion = kFALSE;
    
    //Get Leading Track
    AliESDtrack *leading_track = (AliESDtrack*) fESDevent->GetTrack(leading_track_ID);
    
    //DeltaPhi
    Double_t phi_ref   = TVector2::Phi_0_2pi (leading_track->Phi());
    Double_t phi_trk   = TVector2::Phi_0_2pi (track->Phi());
    Double_t delta_phi = (180.0/TMath::Pi())*TVector2::Phi_0_2pi (phi_trk-phi_ref);

    if ( delta_phi>=120.0 && delta_phi<240.0) isInAwayRegion=kTRUE;
    
    return isInAwayRegion;
}
//__________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskAntiProtons_vs_RT_pPb::IsHighPurityProton (AliESDtrack *track) {
    
    //Variables
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kProton);
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kProton);
    Double_t nsigmaITS = fPIDResponse -> NumberOfSigmasITS (track,AliPID::kProton);
    Double_t pt = track->Pt();

    //selection for high purity antiprotons to build recalibration maps    
    if(!fIsITSrecalib)
    {
        if (pt<0.7 && TMath::Abs(nsigmaTPC)<3.0)                               return kTRUE;
        if (pt>0.7 && TMath::Abs(nsigmaTPC)<2.0 && TMath::Abs(nsigmaTOF)<2.0)  return kTRUE;
    }

    //selection for high purity antiprotons in normal analysis (after recalibration)
    if(fIsITSrecalib)
    {
        Double_t nsigmaITS_recalib = GetRecalibratedITSnsigma (nsigmaITS,track->Eta(),track->P());
        if (pt<0.7 && TMath::Abs(nsigmaITS_recalib)<3.0 && TMath::Abs(nsigmaTPC)<3.0) return kTRUE;
        if (pt>0.7 && TMath::Abs(nsigmaTPC)<2.0 && TMath::Abs(nsigmaTOF)<2.0)         return kTRUE;
    }

    return kFALSE;
}
//__________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskAntiProtons_vs_RT_pPb::IsHighPurityPion (AliESDtrack *track) {
    
    //Variables
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kPion);
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kPion);
    Double_t pt = track->Pt();

    //selection for high purity antiprotons in normal analysis (after recalibration)
    if(fIsITSrecalib)
    {
        if (TMath::Abs(nsigmaTPC)<2.0 && TMath::Abs(nsigmaTOF)<2.0)         return kTRUE;
    }

    return kFALSE;
}
//__________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskAntiProtons_vs_RT_pPb::IsTrackCandidate (AliESDtrack *track){
    
    //Initialization
    Bool_t isTrack=(kFALSE);

    if (!PassedTrackQualityCuts_Syst(track, 0)) return isTrack;

    isTrack = kTRUE;
    return isTrack;
    
}
//__________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskAntiProtons_vs_RT_pPb::PassedTrackQualityCuts_LeadingTrack (AliESDtrack *track)  {
    
    //Initialization
    Bool_t passedTrkSelection=(kFALSE);
    if ( !fESDtrackCuts_LeadingTrack->AcceptTrack (track) ) return passedTrkSelection;
    
    passedTrkSelection = kTRUE;
    return passedTrkSelection;
}
//__________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskAntiProtons_vs_RT_pPb::PassedTrackQualityCuts_TransverseMult (AliESDtrack *track)  {
    
    //Initialization
    Bool_t passedTrkSelection=(kFALSE);
    if ( !fESDtrackCuts_TransverseMult->AcceptTrack (track) ) return passedTrkSelection;
    
    passedTrkSelection = kTRUE;
    return passedTrkSelection;
}
//__________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskAntiProtons_vs_RT_pPb::FillHistograms_StandardCuts (Int_t mult_Transverse, Int_t leading_track_ID, AliESDtrack *track)  {
    
    if(!fIsPion){
        //Variables
        Double_t nsigmaITS = fPIDResponse -> NumberOfSigmasITS (track,AliPID::kProton);
        Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kProton);
        Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kProton);
        Double_t eta       = track->Eta();
        Double_t pt        = track->Pt();
        Double_t p         = track->P();
        Double_t y         = GetRapidity(track);
        Double_t charge    = track->Charge();
        Double_t DCAxy     = GetTransverseDCA (track);
        Double_t DCAz      = GetLongitudinalDCA (track);
        Bool_t   hasTOFhit = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);
        Double_t length    = track->GetIntegratedLength();
        Double_t Rt        = static_cast<Double_t>(mult_Transverse)/fAverage_Nch_Transv;

        //Track Quality Cuts
        if (!PassedTrackQualityCuts_Syst(track,0)) return;
        if (charge > 0) return;

        //DCA_{z} Cut
        if (TMath::Abs(DCAz)>1.0) return;

        if(!fIsITSrecalib)
        {
            //DCA_{xy} Cut
            if (TMath::Abs(DCAxy)>0.1) return;

            //selecting high purity antiprotons
            if(!IsHighPurityProton(track)) return;

            //nsigma_ITS pre-recalibration plots
            Double_t xITS[3] = {nsigmaITS, p, eta};
            hnsigmaITS_antiprotons->Fill(xITS);
        }

        if(fIsITSrecalib)
        {
            //ITS Recalibration
            Double_t nsigmaITS_recalib = GetRecalibratedITSnsigma (nsigmaITS,eta,p);

            //Vectors
            Double_t xITS[3]      = {nsigmaITS_recalib, p, eta};
            Double_t xTPC[3]      = { Rt, nsigmaTPC, pt };
            Double_t xTOF[3]      = { Rt, nsigmaTOF, pt };
            Double_t xDCA[3]      = { Rt, DCAxy, pt };

            //DCA_{xy} Histograms
            if (IsHighPurityProton(track))  {
                if (IsTrackInTowardRegion(track,leading_track_ID))     hDCAxy_antiprotons_Toward     -> Fill (xDCA);
                if (IsTrackInAwayRegion(track,leading_track_ID))       hDCAxy_antiprotons_Away       -> Fill (xDCA);
                if (IsTrackInTransverseRegion(track,leading_track_ID)) hDCAxy_antiprotons_Transverse -> Fill (xDCA);
            }

            //DCA_{xy} Cut
            if (TMath::Abs(DCAxy)>0.1) return;

            //nsigma_ITS post-recalibration plots
            if (IsHighPurityProton(track)) {
                    hnsigmaITS_antiprotons->Fill(xITS);
            }

            //TPC-Only Analysis
            if (pt<1.0 && TMath::Abs(nsigmaITS_recalib) < 3.0)  {
            
                if (IsTrackInTowardRegion(track,leading_track_ID))     hnsigmaTPC_antiprotons_Toward     -> Fill (xTPC);
                if (IsTrackInAwayRegion(track,leading_track_ID))       hnsigmaTPC_antiprotons_Away       -> Fill (xTPC);
                if (IsTrackInTransverseRegion(track,leading_track_ID)) hnsigmaTPC_antiprotons_Transverse -> Fill (xTPC);
            }

            //check on no ITS requirements
            if (pt<1.0)  {
            
                if (IsTrackInTowardRegion(track,leading_track_ID))     hnsigmaTPC_antiprotons_Toward_noITSpresel     -> Fill (xTPC);
                if (IsTrackInAwayRegion(track,leading_track_ID))       hnsigmaTPC_antiprotons_Away_noITSpresel       -> Fill (xTPC);
                if (IsTrackInTransverseRegion(track,leading_track_ID)) hnsigmaTPC_antiprotons_Transverse_noITSpresel -> Fill (xTPC);
            }

            //TOF Analysis
            if (!hasTOFhit) return;
            if (length<350.0) return;
            if (TMath::Abs(nsigmaTPC)>3.0) return;
            if (pt<0.5) return;

            if (IsTrackInTowardRegion(track,leading_track_ID))      hnsigmaTOF_antiprotons_Toward  -> Fill (xTOF);
            if (IsTrackInAwayRegion(track,leading_track_ID))        hnsigmaTOF_antiprotons_Away  -> Fill (xTOF);
            if (IsTrackInTransverseRegion(track,leading_track_ID))  hnsigmaTOF_antiprotons_Transverse  -> Fill (xTOF);
        }
    }

    if(fIsPion){
        //Variables
        Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kPion);
        Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kPion);
        Double_t eta       = track->Eta();
        Double_t pt        = track->Pt();
        Double_t p         = track->P();
        Double_t y         = GetRapidity(track);
        Double_t charge    = track->Charge();
        Double_t DCAxy     = GetTransverseDCA (track);
        Double_t DCAz      = GetLongitudinalDCA (track);
        Bool_t   hasTOFhit = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);
        Double_t length    = track->GetIntegratedLength();
        Double_t Rt        = static_cast<Double_t>(mult_Transverse)/fAverage_Nch_Transv;

        //Track Quality Cuts
        if (!PassedTrackQualityCuts_Syst(track,0)) return;

        //DCA_{z} Cut
        if (TMath::Abs(DCAz)>1.0) return;
        if(fIsITSrecalib)
        {
            //Vectors
            Double_t xTOF[3]      = { Rt, nsigmaTOF, pt };
            Double_t xDCA[3]      = { Rt, DCAxy, pt };

            //DCA_{xy} Histograms
            if (IsHighPurityPion(track))  {
                if(charge > 0) {
                    if (IsTrackInTowardRegion(track,leading_track_ID))     hDCAxy_pos_pions_Toward     -> Fill (xDCA);
                    if (IsTrackInAwayRegion(track,leading_track_ID))       hDCAxy_pos_pions_Away       -> Fill (xDCA);
                    if (IsTrackInTransverseRegion(track,leading_track_ID)) hDCAxy_pos_pions_Transverse -> Fill (xDCA);
                }

                if(charge < 0) {
                    if (IsTrackInTowardRegion(track,leading_track_ID))     hDCAxy_neg_pions_Toward     -> Fill (xDCA);
                    if (IsTrackInAwayRegion(track,leading_track_ID))       hDCAxy_neg_pions_Away       -> Fill (xDCA);
                    if (IsTrackInTransverseRegion(track,leading_track_ID)) hDCAxy_neg_pions_Transverse -> Fill (xDCA);
                }
            }

            //DCA_{xy} Cut
            if (TMath::Abs(DCAxy)>0.1) return;

            //TOF Analysis
            if (!hasTOFhit) return;
            if (length<350.0) return;
            if (TMath::Abs(nsigmaTPC)>3.0) return;
            if (pt<0.3) return;
            
            if (charge > 0) {
                if (IsTrackInTowardRegion(track,leading_track_ID))      hnsigmaTOF_pos_pions_Toward  -> Fill (xTOF);
                if (IsTrackInAwayRegion(track,leading_track_ID))        hnsigmaTOF_pos_pions_Away  -> Fill (xTOF);
                if (IsTrackInTransverseRegion(track,leading_track_ID))  hnsigmaTOF_pos_pions_Transverse  -> Fill (xTOF);
            }

            if (charge < 0) {
                if (IsTrackInTowardRegion(track,leading_track_ID))      hnsigmaTOF_neg_pions_Toward  -> Fill (xTOF);
                if (IsTrackInAwayRegion(track,leading_track_ID))        hnsigmaTOF_neg_pions_Away  -> Fill (xTOF);
                if (IsTrackInTransverseRegion(track,leading_track_ID))  hnsigmaTOF_neg_pions_Transverse  -> Fill (xTOF);
            }
        }
    }
}
//__________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskAntiProtons_vs_RT_pPb::GetRapidity (AliESDtrack *track)  {
    
    //Initialization
    Double_t y(-999);
    
    //Rapidity Calculation
    Double_t mass = AliPID::ParticleMass(AliPID::kProton);
    Double_t p  = track->P();
    Double_t pz = track->Pz();
    Double_t E = TMath::Sqrt(mass*mass + p*p);
    if (E == TMath::Abs(pz)) return -999;
    y = 0.5*TMath::Log((E+pz)/(E-pz));
    
    return y;
}
//__________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskAntiProtons_vs_RT_pPb::FillHistograms_StandardCuts_Sim (AliESDtrack *track)  {
    
    if(!fIsPion)
    {
        //Variables
        Double_t nsigmaITS = fPIDResponse -> NumberOfSigmasITS (track,AliPID::kProton);
        Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kProton);
        Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kProton);
        Double_t pt        = track->Pt();
        Double_t charge    = track->Charge();
        Double_t DCAxy     = GetTransverseDCA (track);
        Double_t DCAz      = GetLongitudinalDCA (track);
        Bool_t   hasTOFhit = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);
        Double_t length    = track->GetIntegratedLength();
        Double_t mass      = AliPID::ParticleMass(AliPID::kProton);
        Double_t eta       = track->Eta();
        Double_t p         = track->P();
        Double_t pz        = track->Pz();
        Double_t E         = TMath::Sqrt(mass*mass + p*p);
        Double_t y_lab     = 0.5*TMath::Log((E+pz)/(E-pz));
        Double_t y_cms     = y_lab;

        //MC Particle
        AliMCParticle *particle = (AliMCParticle*) fMCEvent->GetTrack(TMath::Abs(track->GetLabel()));

        //Track Track Quality Cuts
        if (!PassedTrackQualityCuts_Syst(track,0)) return;

        //DCA_{z} Cut
        if (TMath::Abs(DCAz)>1.0) return;

        //DCA_{xy} distributions
        if (particle->IsPhysicalPrimary()) hDCAxy_antiprotons_prim->Fill(pt, DCAxy);
        if (particle->IsSecondaryFromWeakDecay()) hDCAxy_antiprotons_sec->Fill(pt, DCAxy);

        //DCA_{xy} Cut
        if (TMath::Abs(DCAxy)>0.1) return;
        if (!particle->IsPhysicalPrimary()) return;
        if (charge > 0) return;

        if(!fIsITSrecalib)
        {
            //selecting high purity antiprotons
            if(!IsHighPurityProton(track)) return;

            //nsigma_ITS pre-recalibration plots
            Double_t xITS[3] = {nsigmaITS, p, eta};
            hnsigmaITS_antiprotons->Fill(xITS);
        }

        if(fIsITSrecalib)
        {
            //ITS Recalibration
            Double_t nsigmaITS_recalib = GetRecalibratedITSnsigma (nsigmaITS,eta,p);

            //nsigma_ITS post-recalibration plots
            Double_t xITS[3] = {nsigmaITS_recalib, p, eta};

            //nsigma_ITS post-recalibration plots
            if (IsHighPurityProton(track)) {
                    hnsigmaITS_antiprotons->Fill(xITS);
            }

            //TPC-Only Analysis
            if (pt<1.0 && TMath::Abs(nsigmaITS_recalib) < 3.0)  {hnsigmaTPC_antiprotons_Rec -> Fill (pt,nsigmaTPC);}

            //TOF Analysis
            if (!hasTOFhit) return;
            if (length<350.0) return;
            if (TMath::Abs(nsigmaTPC)>3.0) return;
            if (pt<0.3) return;

            hnsigmaTOF_antiprotons_Rec -> Fill(pt,nsigmaTOF);
        }
    }

    if(fIsPion)
    {
        //Variables
        Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kPion);
        Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kPion);
        Double_t pt        = track->Pt();
        Double_t charge    = track->Charge();
        Double_t DCAxy     = GetTransverseDCA (track);
        Double_t DCAz      = GetLongitudinalDCA (track);
        Bool_t   hasTOFhit = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);
        Double_t length    = track->GetIntegratedLength();
        Double_t mass      = AliPID::ParticleMass(AliPID::kPion);
        Double_t eta       = track->Eta();
        Double_t p         = track->P();
        Double_t pz        = track->Pz();
        Double_t E         = TMath::Sqrt(mass*mass + p*p);
        Double_t y_lab     = 0.5*TMath::Log((E+pz)/(E-pz));
        Double_t y_cms     = y_lab;

        //MC Particle
        AliMCParticle *particle = (AliMCParticle*) fMCEvent->GetTrack(TMath::Abs(track->GetLabel()));

        //Track Track Quality Cuts
        if (!PassedTrackQualityCuts_Syst(track,0)) return;

        //DCA_{z} Cut
        if (TMath::Abs(DCAz)>1.0) return;

        //DCA_{xy} distributions
        if (particle->IsPhysicalPrimary())
        {
            if (charge > 0) hDCAxy_pos_pions_prim->Fill(pt, DCAxy);
            if (charge < 0) hDCAxy_neg_pions_prim->Fill(pt, DCAxy);
        } 
        if (particle->IsSecondaryFromWeakDecay())
        {
            if (charge > 0) hDCAxy_pos_pions_sec->Fill(pt, DCAxy);
            if (charge < 0) hDCAxy_neg_pions_sec->Fill(pt, DCAxy);
        }

        //DCA_{xy} Cut
        if (TMath::Abs(DCAxy)>0.1) return;
        if (!particle->IsPhysicalPrimary()) return;

        if(fIsITSrecalib)
        {
            //TOF Analysis
            if (!hasTOFhit) return;
            if (length<350.0) return;
            if (TMath::Abs(nsigmaTPC)>3.0) return;
            if (pt<0.3) return;

            if (charge > 0) hnsigmaTOF_pos_pions_Rec -> Fill(pt,nsigmaTOF);
            if (charge < 0) hnsigmaTOF_neg_pions_Rec -> Fill(pt,nsigmaTOF);
        }
    }
    
    
}
//__________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskAntiProtons_vs_RT_pPb::FillHistograms_Systematics  (Int_t mult_Transverse, Int_t leading_track_ID, AliESDtrack *track, Int_t isyst)  {
    
    //ITS Pre-selection
    Double_t nsigmaITSmax[50] = {3.0,2.62945,2.69333,3.1122,3.28161,2.70606,2.89633,2.70347,3.31368,2.65839,2.93157,3.10837,3.43406,2.51501,3.49459,3.27995,3.05988,2.89139,2.60827,3.16441,2.69537,2.71876,2.67132,2.80894,3.42022,3.23984,2.96282,2.84138,2.78734,3.102,3.15273,2.95713,2.90697,3.05029,2.86195,3.16312,3.43193,2.57391,3.38719,3.177,3.38128,3.40739,2.83285,2.7041,3.05584,2.58265,2.99982,3.40591,2.64627,3.45544};
    //TPC Pre-selection
    Double_t nsigmaTPCmax[50] = {3.0,2.62945,2.69333,3.1122,3.28161,2.70606,2.89633,2.70347,3.31368,2.65839,2.93157,3.10837,3.43406,2.51501,3.49459,3.27995,3.05988,2.89139,2.60827,3.16441,2.69537,2.71876,2.67132,2.80894,3.42022,3.23984,2.96282,2.84138,2.78734,3.102,3.15273,2.95713,2.90697,3.05029,2.86195,3.16312,3.43193,2.57391,3.38719,3.177,3.38128,3.40739,2.83285,2.7041,3.05584,2.58265,2.99982,3.40591,2.64627,3.45544};
    //DCAz preselection
    Double_t dcaz_max[50] = {1.0,0.524845,0.838556,1.36511,0.972847,0.644151,1.14321,0.557206,0.978874,1.14591,0.958042,0.950637,0.718069,0.520478,1.26181,0.958435,0.532273,1.13665,1.1004,1.37577,1.02552,0.571379,0.708928,0.831586,1.44803,0.90449,1.0531,0.627843,1.42758,1.4044,0.56425,0.901189,0.987433,0.816606,1.37016,1.32486,1.30484,0.572152,1.34936,0.535521,0.819561,0.575574,1.44179,1.25525,1.36479,0.522065,0.535635,1.26045,0.851564,0.599568};
    //DCAxy preselection
    Double_t dcaxy_max[]={0.0800,0.1200,0.1600,0.1400,0.1000,0.1000,0.1200,0.1400,0.1000,0.1200,0.1400,0.1400,0.1200,0.1200,0.1200,0.1600,0.1000,0.1000,0.1200,0.1000,0.1000,0.1200,0.1000,0.1200,0.1000,0.1400,0.1400,0.1400,0.0800,0.1600,0.1400,0.1000,0.1200,0.1400,0.1400,0.0800,0.1200,0.1200,0.1600,0.1000,0.1000,0.1600,0.1200,0.1400,0.1200,0.1600,0.1000,0.1200,0.1000,0.1400};

    if(!fIsPion)
    {
        //Variables
        Double_t nsigmaITS = fPIDResponse -> NumberOfSigmasITS (track,AliPID::kProton);
        Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kProton);
        Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kProton);
        Double_t pt        = track->Pt();
        Double_t charge    = track->Charge();
        Double_t DCAxy     = GetTransverseDCA (track);
        Double_t DCAz      = GetLongitudinalDCA (track);
        Bool_t   hasTOFhit = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);
        Double_t length    = track->GetIntegratedLength();
        Double_t Rt        = static_cast<Double_t>(mult_Transverse)/fAverage_Nch_Transv;
        Double_t eta       = track->Eta();
        Double_t p         = track->P();

        //ITS Recalibration
        Double_t nsigmaITS_recalib = GetRecalibratedITSnsigma (nsigmaITS,eta,p);

        //Track Quality Cuts
        if (!PassedTrackQualityCuts_Syst(track,isyst)) return;
        if (charge > 0) return;

        //Vectors
        Double_t xTPC[4]      = { Rt, nsigmaTPC, pt, static_cast<Double_t>(isyst) };
        Double_t xTOF[4]      = { Rt, nsigmaTOF, pt, static_cast<Double_t>(isyst) };
        Double_t xDCA[4]      = { Rt, DCAxy, pt, static_cast<Double_t>(isyst) };
        //DCAz cut
        if(TMath::Abs(DCAz)>dcaz_max[isyst]) return;
        //DCA_{xy} Histograms
        if (IsHighPurityProton (track))  {hDCAxy_antiprotons_Syst -> Fill (xDCA);}

        //DCAxy cuts
        if (TMath::Abs(DCAxy)>dcaxy_max[isyst]) return;

        //TPC-Only Analysis
        if (pt<1.0 && TMath::Abs(nsigmaITS_recalib) < nsigmaITSmax[isyst])  {hnsigmaTPC_antiprotons_Syst -> Fill (xTPC);}

        //TOF Analysis
        if (!hasTOFhit) return;
        if (length<350.0) return;
        if (TMath::Abs(nsigmaTPC)> nsigmaTPCmax[isyst]) return;
        if (pt<0.3) return;
            
        hnsigmaTOF_antiprotons_Syst  -> Fill (xTOF);
    }

    if(fIsPion)
    {
        //Variables
        Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kPion);
        Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kPion);
        Double_t pt        = track->Pt();
        Double_t charge    = track->Charge();
        Double_t DCAxy     = GetTransverseDCA (track);
        Double_t DCAz      = GetLongitudinalDCA (track);
        Bool_t   hasTOFhit = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);
        Double_t length    = track->GetIntegratedLength();
        Double_t Rt        = static_cast<Double_t>(mult_Transverse)/fAverage_Nch_Transv;
        Double_t eta       = track->Eta();
        Double_t p         = track->P();

        //Track Quality Cuts
        if (!PassedTrackQualityCuts_Syst(track,isyst)) return;

        //Vectors
        Double_t xTOF[4]      = { Rt, nsigmaTOF, pt, static_cast<Double_t>(isyst) };
        Double_t xDCA[4]      = { Rt, DCAxy, pt, static_cast<Double_t>(isyst) };
        //DCAz cut
        if(TMath::Abs(DCAz)>dcaz_max[isyst]) return;
        //DCA_{xy} Histograms
        if (IsHighPurityProton (track))  
        {
            if (charge > 0) hDCAxy_pos_pions_Syst -> Fill (xDCA);
            if (charge < 0) hDCAxy_neg_pions_Syst -> Fill (xDCA);
        }

        //DCAxy cuts
        if (TMath::Abs(DCAxy)>dcaxy_max[isyst]) return;

        //TOF Analysis
        if (!hasTOFhit) return;
        if (length<350.0) return;
        if (TMath::Abs(nsigmaTPC)> nsigmaTPCmax[isyst]) return;
        if (pt<0.3) return;
            
        if (charge > 0) hnsigmaTOF_pos_pions_Syst  -> Fill (xTOF);
        if (charge < 0) hnsigmaTOF_neg_pions_Syst  -> Fill (xTOF);

    }   
}
//__________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskAntiProtons_vs_RT_pPb::FillHistograms_Systematics_Sim  (AliESDtrack *track, Int_t isyst)  {
    
    //ITS Pre-selection
    Double_t nsigmaITSmax[50] = {3.0,2.62945,2.69333,3.1122,3.28161,2.70606,2.89633,2.70347,3.31368,2.65839,2.93157,3.10837,3.43406,2.51501,3.49459,3.27995,3.05988,2.89139,2.60827,3.16441,2.69537,2.71876,2.67132,2.80894,3.42022,3.23984,2.96282,2.84138,2.78734,3.102,3.15273,2.95713,2.90697,3.05029,2.86195,3.16312,3.43193,2.57391,3.38719,3.177,3.38128,3.40739,2.83285,2.7041,3.05584,2.58265,2.99982,3.40591,2.64627,3.45544};
    //TPC Pre-selection
    Double_t nsigmaTPCmax[50] = {3.0,2.62945,2.69333,3.1122,3.28161,2.70606,2.89633,2.70347,3.31368,2.65839,2.93157,3.10837,3.43406,2.51501,3.49459,3.27995,3.05988,2.89139,2.60827,3.16441,2.69537,2.71876,2.67132,2.80894,3.42022,3.23984,2.96282,2.84138,2.78734,3.102,3.15273,2.95713,2.90697,3.05029,2.86195,3.16312,3.43193,2.57391,3.38719,3.177,3.38128,3.40739,2.83285,2.7041,3.05584,2.58265,2.99982,3.40591,2.64627,3.45544};
    //DCAz preselection
    Double_t dcaz_max[50] = {1.0,0.524845,0.838556,1.36511,0.972847,0.644151,1.14321,0.557206,0.978874,1.14591,0.958042,0.950637,0.718069,0.520478,1.26181,0.958435,0.532273,1.13665,1.1004,1.37577,1.02552,0.571379,0.708928,0.831586,1.44803,0.90449,1.0531,0.627843,1.42758,1.4044,0.56425,0.901189,0.987433,0.816606,1.37016,1.32486,1.30484,0.572152,1.34936,0.535521,0.819561,0.575574,1.44179,1.25525,1.36479,0.522065,0.535635,1.26045,0.851564,0.599568};
    //DCAxy preselection
    Double_t dcaxy_max[]={0.0800,0.1200,0.1600,0.1400,0.1000,0.1000,0.1200,0.1400,0.1000,0.1200,0.1400,0.1400,0.1200,0.1200,0.1200,0.1600,0.1000,0.1000,0.1200,0.1000,0.1000,0.1200,0.1000,0.1200,0.1000,0.1400,0.1400,0.1400,0.0800,0.1600,0.1400,0.1000,0.1200,0.1400,0.1400,0.0800,0.1200,0.1200,0.1600,0.1000,0.1000,0.1600,0.1200,0.1400,0.1200,0.1600,0.1000,0.1200,0.1000,0.1400};

    if(!fIsPion)
    {
        //Variables
        Double_t nsigmaITS = fPIDResponse -> NumberOfSigmasITS (track,AliPID::kProton);
        Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kProton);
        Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kProton);
        Double_t pt        = track->Pt();
        Double_t charge    = track->Charge();
        Double_t DCAxy     = GetTransverseDCA (track);
        Double_t DCAz      = GetLongitudinalDCA(track);
        Bool_t   hasTOFhit = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);
        Double_t length    = track->GetIntegratedLength();
        Double_t eta       = track->Eta();
        Double_t p         = track->P();

        //ITS Recalibration
        Double_t nsigmaITS_recalib = GetRecalibratedITSnsigma (nsigmaITS,eta,p);

        //MC Particle
        AliMCParticle *particle = (AliMCParticle*) fMCEvent->GetTrack(TMath::Abs(track->GetLabel()));

        //Track Quality Cuts
        if (!PassedTrackQualityCuts_Syst(track,isyst)) return;
    
        //DCA_{xy} syst distributions
        if(TMath::Abs(DCAz)>dcaz_max[isyst]) return;
        Double_t xDCA[3] = {DCAxy, pt, static_cast<Double_t>(isyst)};
        if(particle->IsPhysicalPrimary()) hDCAxy_antiprotons_prim_Syst->Fill(xDCA);
        if(particle->IsSecondaryFromWeakDecay()) hDCAxy_antiprotons_sec_Syst->Fill(xDCA);

        //Primary antiprotons selection
        if (!particle->IsPhysicalPrimary())     return;
        if (charge > 0 ) return;

        //DCA_xy cut
        if(TMath::Abs(DCAxy)>dcaxy_max[isyst]) return;

        //TPC-Only Analysis
        if (pt<1.0 && TMath::Abs(nsigmaITS_recalib) < nsigmaITSmax[isyst])  {hnsigmaTPC_antiprotons_Rec_Syst -> Fill (pt,isyst);}

        //TOF Analysis
        if (!hasTOFhit) return;
        if (length<350.0) return;
        if (TMath::Abs(nsigmaTPC)>nsigmaTPCmax[isyst]) return;
        if (pt<0.5) return;

        hnsigmaTOF_antiprotons_Rec_Syst -> Fill (pt,isyst);
    }

    if(fIsPion)
    {
        //Variables
        Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kPion);
        Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kPion);
        Double_t pt        = track->Pt();
        Double_t charge    = track->Charge();
        Double_t DCAxy     = GetTransverseDCA (track);
        Double_t DCAz      = GetLongitudinalDCA(track);
        Bool_t   hasTOFhit = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);
        Double_t length    = track->GetIntegratedLength();
        Double_t eta       = track->Eta();
        Double_t p         = track->P();

        //MC Particle
        AliMCParticle *particle = (AliMCParticle*) fMCEvent->GetTrack(TMath::Abs(track->GetLabel()));

        //Track Quality Cuts
        if (!PassedTrackQualityCuts_Syst(track,isyst)) return;
    
        //DCA_{xy} syst distributions
        if(TMath::Abs(DCAz)>dcaz_max[isyst]) return;
        Double_t xDCA[3] = {DCAxy, pt, static_cast<Double_t>(isyst)};
        if(particle->IsPhysicalPrimary()) 
        {
            if (charge > 0) hDCAxy_pos_pions_prim_Syst->Fill(xDCA);
            if (charge < 0) hDCAxy_neg_pions_prim_Syst->Fill(xDCA);

        }
        if(particle->IsSecondaryFromWeakDecay()) 
        {
            if (charge > 0) hDCAxy_pos_pions_sec_Syst->Fill(xDCA);
            if (charge < 0) hDCAxy_neg_pions_sec_Syst->Fill(xDCA);
        }
        

        //Primary antiprotons selection
        if (!particle->IsPhysicalPrimary())     return;
        //DCA_xy cut
        if(TMath::Abs(DCAxy)>dcaxy_max[isyst]) return;
        //TOF Analysis
        if (!hasTOFhit) return;
        if (length<350.0) return;
        if (TMath::Abs(nsigmaTPC)>nsigmaTPCmax[isyst]) return;
        if (pt<0.3) return;

        if (charge > 0) hnsigmaTOF_pos_pions_Rec_Syst -> Fill (pt,isyst);
        if (charge < 0) hnsigmaTOF_neg_pions_Rec_Syst -> Fill (pt,isyst);
    }
}
//__________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskAntiProtons_vs_RT_pPb::PassedTrackQualityCuts_Syst (AliESDtrack *track, Int_t isyst)  {

    //Initialization
    Bool_t passedTrkSelection=(kFALSE);

    //Defining variations
    Int_t nTPCcr_min[50] = {80,97,71,75,86,87,74,71,71,91,82,77,75,85,80,88,94,83,85,82,74,86,75,93,96,96,84,83,71,94,84,82,72,85,76,86,83,82,93,90,92,92,70,83,91,93,97,82,80,95};

    Int_t nITScls_min[50] = {4,4,3,3,6,3,4,4,5,4,4,6,3,4,4,3,4,4,3,4,4,5,4,3,5,4,6,4,5,4,3,4,5,4,4,6,4,4,4,5,4,4,5,4,3,6,4,3,4,4};

    Int_t nTPCclsdEdx_min[50] = {60,63,54,57,57,64,69,60,54,62,66,50,55,64,68,50,67,59,52,64,67,64,50,55,56,52,68,70,54,57,62,63,50,50,68,64,63,67,53,58,54,54,68,56,51,56,63,59,56,58};

    Double_t chi2ITS_NDF_max[50] = {36,46.1797,36.6257,48.2043,36.8251,46.3546,44.4219,36.9071,36.4452,45.4235,44.8695,47.7977,47.0017,44.1026,37.5981,45.6092,49.7681,44.8784,38.9593,39.0638,48.7178,49.4983,49.3736,45.7836,41.407,44.4918,41.3072,36.7743,47.0045,43.5901,45.2098,43.167,46.9972,39.2736,38.6093,42.4406,47.2409,43.3473,45.1781,41.2574,45.8378,36.263,47.7494,48.1092,39.5733,42.8802,46.6033,47.2696,47.3724,45.3757};

    Double_t chi2TPC_NDF_max[50] = {4.0,3.99601,2.79669,3.88701,3.73929,3.07139,4.14759,3.19747,4.58796,4.53536,3.56903,4.85172,3.70836,5.01696,5.05308,4.16061,3.40112,4.67905,4.40227,5.0066,4.95311,3.13955,3.35452,4.43888,4.21231,3.8571,3.43539,5.16655,5.19831,4.4465,4.1092,4.57635,2.61833,4.17145,4.14075,3.91889,5.05004,3.12027,3.86329,3.47401,5.27806,3.50963,4.18327,4.7357,3.71845,2.67676,3.2959,4.27114,2.86657,3.1244};

    Double_t dcaz_max[50] = {1.0,0.524845,0.838556,1.36511,0.972847,0.644151,1.14321,0.557206,0.978874,1.14591,0.958042,0.950637,0.718069,0.520478,1.26181,0.958435,0.532273,1.13665,1.1004,1.37577,1.02552,0.571379,0.708928,0.831586,1.44803,0.90449,1.0531,0.627843,1.42758,1.4044,0.56425,0.901189,0.987433,0.816606,1.37016,1.32486,1.30484,0.572152,1.34936,0.535521,0.819561,0.575574,1.44179,1.25525,1.36479,0.522065,0.535635,1.26045,0.851564,0.599568};

    Double_t cr_over_findable_min[50] = {0.8,0.849002,0.781492,0.897933,0.704228,0.706126,0.704632,0.709122,0.897615,0.727589,0.765229,0.834458,0.878155,0.824003,0.893187,0.824451,0.76873,0.796163,0.713176,0.736468,0.722194,0.79335,0.892695,0.81785,0.702785,0.81674,0.859769,0.845479,0.729515,0.848578,0.831617,0.782084,0.817316,0.873488,0.849991,0.717202,0.827419,0.885216,0.870209,0.801342,0.814752,0.867845,0.866332,0.738925,0.778826,0.793383,0.812433,0.718611,0.704671,0.889181};

    Double_t dcaxy_max[]={0.0800,0.1200,0.1600,0.1400,0.1000,0.1000,0.1200,0.1400,0.1000,0.1200,0.1400,0.1400,0.1200,0.1200,0.1200,0.1600,0.1000,0.1000,0.1200,0.1000,0.1000,0.1200,0.1000,0.1200,0.1000,0.1400,0.1400,0.1400,0.0800,0.1600,0.1400,0.1000,0.1200,0.1400,0.1400,0.0800,0.1200,0.1200,0.1600,0.1000,0.1000,0.1600,0.1200,0.1400,0.1200,0.1600,0.1000,0.1200,0.1000,0.1400};


    //Track Variables
    Double_t nCrossedRows               = (Double_t)track->GetTPCCrossedRows();
    Double_t nCrossedRows_over_Findable = ((Double_t)track->GetTPCCrossedRows())/((Double_t)track->GetTPCNclsF());
    Double_t nClustersITS               = (Double_t)track->GetITSNcls();
    Double_t chi2TPC_ndf                = ((Double_t)track->GetTPCchi2())/((Double_t)track->GetTPCNcls());
    Double_t chi2ITS_ndf                = ((Double_t)track->GetITSchi2())/((Double_t)track->GetITSNcls());
    Double_t nClustersTPC_dEdx          = (Double_t)track->GetTPCsignalN();
    Double_t DCAz                       = GetLongitudinalDCA(track);
    Double_t DCAxy                      = GetTransverseDCA(track);
    Double_t y_lab                      = GetRapidity(track);
    Double_t y_cms                      = y_lab-0.465;

    //Analysis Parameters
    Double_t nCrossedRows_Min               = nTPCcr_min[isyst];
    Double_t nCrossedRows_over_Findable_Min = cr_over_findable_min[isyst];
    Double_t nClustersITS_Min               = nITScls_min[isyst];
    Double_t chi2TPC_ndf_Max                = chi2TPC_NDF_max[isyst];
    Double_t chi2ITS_ndf_Max                = chi2ITS_NDF_max[isyst];
    Double_t nClustersTPC_dEdx_Min          = nTPCclsdEdx_min[isyst];
    Double_t DCAz_max                       = dcaz_max[isyst];
    Double_t DCAxy_max                      = dcaxy_max[isyst];

    //Cuts
    if (y_cms > 0)                                                 return passedTrkSelection;
    if (y_cms < -1)                                                return passedTrkSelection;
    if (TMath::Abs(track->Eta())>0.8)                              return passedTrkSelection;
    if (nCrossedRows<nCrossedRows_Min)                             return passedTrkSelection;
    if (nCrossedRows_over_Findable<nCrossedRows_over_Findable_Min) return passedTrkSelection;
    if (nClustersITS<nClustersITS_Min)                             return passedTrkSelection;
    if (chi2TPC_ndf>chi2TPC_ndf_Max)                               return passedTrkSelection;
    if (chi2ITS_ndf>chi2ITS_ndf_Max)                               return passedTrkSelection;
    if (nClustersTPC_dEdx<nClustersTPC_dEdx_Min)                   return passedTrkSelection;


    passedTrkSelection = kTRUE;
    return passedTrkSelection;
}
//__________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskAntiProtons_vs_RT_pPb::GetTransverseDCA (AliESDtrack *track)  {
         
    Double_t impactParameter[2];
    Double_t covarianceMatrix[3];
    if (!track->PropagateToDCA (fESDevent->GetPrimaryVertex(),fESDevent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;
            
    Double_t DCAxy = impactParameter[0];
            
    return DCAxy;
}
//__________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskAntiProtons_vs_RT_pPb::GetLongitudinalDCA (AliESDtrack *track)  {
         
    Double_t impactParameter[2];
    Double_t covarianceMatrix[3];
    if (!track->PropagateToDCA (fESDevent->GetPrimaryVertex(),fESDevent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;
            
    Double_t DCAz = impactParameter[1];
            
    return DCAz;
}
//__________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskAntiProtons_vs_RT_pPb::GetRecalibratedITSnsigma (Double_t nsigma, Double_t eta, Double_t p)  {
    
    //Initialization
    Double_t nsigma_corr=nsigma;
    
    //Protections
    if (eta<-0.8) return nsigma_corr;
    if (eta>+0.8) return nsigma_corr;
    if (p>1.0)    return nsigma_corr;
    if (p<0.3)    return nsigma_corr;

    //Get Mean and Width from Maps
    Int_t ix = hMean -> GetXaxis()->FindBin(eta);
    Int_t iy = hMean -> GetYaxis()->FindBin(p);
    Double_t x0 = hMean ->GetBinContent (ix,iy);
    Double_t w  = hWidth->GetBinContent (ix,iy);

    nsigma_corr = (nsigma-x0)/w;
    return nsigma_corr;
}
//__________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskAntiProtons_vs_RT_pPb::Terminate(Option_t *)  {
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
