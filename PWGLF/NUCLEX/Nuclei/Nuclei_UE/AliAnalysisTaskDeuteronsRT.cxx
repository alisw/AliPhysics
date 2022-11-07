#include "AliAnalysisTaskDeuteronsRT.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliAnalysisTask.h"
#include "AliESDtrackCuts.h"
#include "TLorentzVector.h"
#include "AliPIDResponse.h"
#include "TLorentzVector.h"
#include "AliMCParticle.h"
#include "AliCentrality.h"
#include "AliEventCuts.h"
#include "TDatabasePDG.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliVVZERO.h"
#include "THnSparse.h"
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


using namespace std;
ClassImp(AliAnalysisTaskDeuteronsRT)


//__________________________________________________________________________________________________________________________________________________
AliAnalysisTaskDeuteronsRT::AliAnalysisTaskDeuteronsRT():
AliAnalysisTaskSE(),
fESDevent(nullptr),
fMCevent(nullptr),
fMCEventHandler(nullptr),
fESDtrackCuts_Deuteron(nullptr),
fESDtrackCuts_LeadingTrack(nullptr),
fESDtrackCuts_TransverseMult(nullptr),
fPIDResponse(nullptr),
fOutputList(nullptr),
fQAList(nullptr),
fTriggerType(AliVEvent::kINT7),
fMultMin(0),
fMultMax(100),
fAverage_Nch_Transv(6.81),
hAnalysisParameters(nullptr),
fIspPb(kFALSE),
fIsMC(kFALSE),
fPt_min_leading(3),
fIsUEanalysis(kTRUE),
fESDeventCuts(),
hNumberOfEvents(nullptr),
hMultPercentile(nullptr),
hMultTransverse(nullptr),
hMultToward(nullptr),
hMultAway(nullptr),
hRtDistribution(nullptr),
hEventsWithLeadingTrack(nullptr),
hNumberOfDeuterons(nullptr),
hNchTransv_NchTot(nullptr),
hRt_NchTot(nullptr),
hNchTransv_MultPercentile(nullptr),
hRt_MultPercentile(nullptr),
hnsigmaTPC_deuterons_Toward(nullptr),
hnsigmaTPC_deuterons_Away(nullptr),
hnsigmaTPC_deuterons_Transverse(nullptr),
hnsigmaTPC_antideuterons_Toward(nullptr),
hnsigmaTPC_antideuterons_Away(nullptr),
hnsigmaTPC_antideuterons_Transverse(nullptr),
hnsigmaTOF_deuterons_Toward(nullptr),
hnsigmaTOF_deuterons_Away(nullptr),
hnsigmaTOF_deuterons_Transverse(nullptr),
hnsigmaTOF_antideuterons_Toward(nullptr),
hnsigmaTOF_antideuterons_Away(nullptr),
hnsigmaTOF_antideuterons_Transverse(nullptr),
hDCAxy_deuterons_Toward(nullptr),
hDCAxy_deuterons_Away(nullptr),
hDCAxy_deuterons_Transverse(nullptr),
hDCAxy_antideuterons_Toward(nullptr),
hDCAxy_antideuterons_Away(nullptr),
hDCAxy_antideuterons_Transverse(nullptr),
hnsigmaTPC_deuterons_Syst(nullptr),
hnsigmaTPC_antideuterons_Syst(nullptr),
hnsigmaTOF_deuterons_Syst(nullptr),
hnsigmaTOF_antideuterons_Syst(nullptr),
hDCAxy_deuterons_Syst(nullptr),
hDCAxy_antideuterons_Syst(nullptr),
hnsigmaTPC_deuterons_rap(nullptr),
hnsigmaTPC_antideuterons_rap(nullptr),
hnsigmaTOF_deuterons_rap(nullptr),
hnsigmaTOF_antideuterons_rap(nullptr),
hnsigmaTPC_deuterons_rap_Syst(nullptr),
hnsigmaTPC_antideuterons_rap_Syst(nullptr),
hnsigmaTOF_deuterons_rap_Syst(nullptr),
hnsigmaTOF_antideuterons_rap_Syst(nullptr),
h_deuterons_Gen(nullptr),
h_antideuterons_Gen(nullptr),
hnsigmaTPC_deuterons_Rec(nullptr),
hnsigmaTPC_antideuterons_Rec(nullptr),
hnsigmaTOF_deuterons_Rec(nullptr),
hnsigmaTOF_antideuterons_Rec(nullptr),
hDCAxy_deuterons_Prim(nullptr),
hDCAxy_deuterons_Sec(nullptr),
hnsigmaTPC_deuterons_Rec_Syst(nullptr),
hnsigmaTPC_antideuterons_Rec_Syst(nullptr),
hnsigmaTOF_deuterons_Rec_Syst(nullptr),
hnsigmaTOF_antideuterons_Rec_Syst(nullptr),
hGeneratedDeuterons_vs_Rapidity(nullptr),
hGeneratedAntiDeuterons_vs_Rapidity(nullptr),
hReconstructedDeuterons_TPC_vs_Rapidity(nullptr),
hReconstructedAntiDeuterons_TPC_vs_Rapidity(nullptr),
hReconstructedDeuterons_TOF_vs_Rapidity(nullptr),
hReconstructedAntiDeuterons_TOF_vs_Rapidity(nullptr),
hnsigmaTPC_deuterons_Rec_rap_Syst(nullptr),
hnsigmaTPC_antideuterons_Rec_rap_Syst(nullptr),
hnsigmaTOF_deuterons_Rec_rap_Syst(nullptr),
hnsigmaTOF_antideuterons_Rec_rap_Syst(nullptr)
{}
//__________________________________________________________________________________________________________________________________________________
AliAnalysisTaskDeuteronsRT::AliAnalysisTaskDeuteronsRT(const char *name):
AliAnalysisTaskSE(name),
fESDevent(nullptr),
fMCevent(nullptr),
fMCEventHandler(nullptr),
fESDtrackCuts_Deuteron(nullptr),
fESDtrackCuts_LeadingTrack(nullptr),
fESDtrackCuts_TransverseMult(nullptr),
fPIDResponse(nullptr),
fOutputList(nullptr),
fQAList(nullptr),
fTriggerType(AliVEvent::kINT7),
fMultMin(0),
fMultMax(100),
fAverage_Nch_Transv(6.81),
hAnalysisParameters(nullptr),
fIspPb(kFALSE),
fIsMC(kFALSE),
fPt_min_leading(3),
fIsUEanalysis(kTRUE),
fESDeventCuts(),
hNumberOfEvents(nullptr),
hMultPercentile(nullptr),
hMultTransverse(nullptr),
hMultToward(nullptr),
hMultAway(nullptr),
hRtDistribution(nullptr),
hEventsWithLeadingTrack(nullptr),
hNumberOfDeuterons(nullptr),
hNchTransv_NchTot(nullptr),
hRt_NchTot(nullptr),
hNchTransv_MultPercentile(nullptr),
hRt_MultPercentile(nullptr),
hnsigmaTPC_deuterons_Toward(nullptr),
hnsigmaTPC_deuterons_Away(nullptr),
hnsigmaTPC_deuterons_Transverse(nullptr),
hnsigmaTPC_antideuterons_Toward(nullptr),
hnsigmaTPC_antideuterons_Away(nullptr),
hnsigmaTPC_antideuterons_Transverse(nullptr),
hnsigmaTOF_deuterons_Toward(nullptr),
hnsigmaTOF_deuterons_Away(nullptr),
hnsigmaTOF_deuterons_Transverse(nullptr),
hnsigmaTOF_antideuterons_Toward(nullptr),
hnsigmaTOF_antideuterons_Away(nullptr),
hnsigmaTOF_antideuterons_Transverse(nullptr),
hDCAxy_deuterons_Toward(nullptr),
hDCAxy_deuterons_Away(nullptr),
hDCAxy_deuterons_Transverse(nullptr),
hDCAxy_antideuterons_Toward(nullptr),
hDCAxy_antideuterons_Away(nullptr),
hDCAxy_antideuterons_Transverse(nullptr),
hnsigmaTPC_deuterons_Syst(nullptr),
hnsigmaTPC_antideuterons_Syst(nullptr),
hnsigmaTOF_deuterons_Syst(nullptr),
hnsigmaTOF_antideuterons_Syst(nullptr),
hDCAxy_deuterons_Syst(nullptr),
hDCAxy_antideuterons_Syst(nullptr),
hnsigmaTPC_deuterons_rap(nullptr),
hnsigmaTPC_antideuterons_rap(nullptr),
hnsigmaTOF_deuterons_rap(nullptr),
hnsigmaTOF_antideuterons_rap(nullptr),
hnsigmaTPC_deuterons_rap_Syst(nullptr),
hnsigmaTPC_antideuterons_rap_Syst(nullptr),
hnsigmaTOF_deuterons_rap_Syst(nullptr),
hnsigmaTOF_antideuterons_rap_Syst(nullptr),
h_deuterons_Gen(nullptr),
h_antideuterons_Gen(nullptr),
hnsigmaTPC_deuterons_Rec(nullptr),
hnsigmaTPC_antideuterons_Rec(nullptr),
hnsigmaTOF_deuterons_Rec(nullptr),
hnsigmaTOF_antideuterons_Rec(nullptr),
hDCAxy_deuterons_Prim(nullptr),
hDCAxy_deuterons_Sec(nullptr),
hnsigmaTPC_deuterons_Rec_Syst(nullptr),
hnsigmaTPC_antideuterons_Rec_Syst(nullptr),
hnsigmaTOF_deuterons_Rec_Syst(nullptr),
hnsigmaTOF_antideuterons_Rec_Syst(nullptr),
hGeneratedDeuterons_vs_Rapidity(nullptr),
hGeneratedAntiDeuterons_vs_Rapidity(nullptr),
hReconstructedDeuterons_TPC_vs_Rapidity(nullptr),
hReconstructedAntiDeuterons_TPC_vs_Rapidity(nullptr),
hReconstructedDeuterons_TOF_vs_Rapidity(nullptr),
hReconstructedAntiDeuterons_TOF_vs_Rapidity(nullptr),
hnsigmaTPC_deuterons_Rec_rap_Syst(nullptr),
hnsigmaTPC_antideuterons_Rec_rap_Syst(nullptr),
hnsigmaTOF_deuterons_Rec_rap_Syst(nullptr),
hnsigmaTOF_antideuterons_Rec_rap_Syst(nullptr)
{
    DefineInput (0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
}
//__________________________________________________________________________________________________________________________________________________
AliAnalysisTaskDeuteronsRT::~AliAnalysisTaskDeuteronsRT()  {
    
    fOutputList->Clear();
    
    delete fESDevent;
    delete fMCevent;
    delete fMCEventHandler;
    delete fESDtrackCuts_Deuteron;
    delete fESDtrackCuts_LeadingTrack;
    delete fESDtrackCuts_TransverseMult;
    delete fPIDResponse;
    delete fOutputList;
    delete fQAList;
    delete hAnalysisParameters;
    delete hNumberOfEvents;
    delete hMultPercentile;
    delete hMultTransverse;
    delete hMultToward;
    delete hMultAway;
    delete hRtDistribution;
    delete hEventsWithLeadingTrack;
    delete hNumberOfDeuterons;
    delete hNchTransv_NchTot;
    delete hRt_NchTot;
    delete hNchTransv_MultPercentile;
    delete hRt_MultPercentile;
    delete hnsigmaTPC_deuterons_Toward;
    delete hnsigmaTPC_deuterons_Away;
    delete hnsigmaTPC_deuterons_Transverse;
    delete hnsigmaTPC_antideuterons_Toward;
    delete hnsigmaTPC_antideuterons_Away;
    delete hnsigmaTPC_antideuterons_Transverse;
    delete hnsigmaTOF_deuterons_Toward;
    delete hnsigmaTOF_deuterons_Away;
    delete hnsigmaTOF_deuterons_Transverse;
    delete hnsigmaTOF_antideuterons_Toward;
    delete hnsigmaTOF_antideuterons_Away;
    delete hnsigmaTOF_antideuterons_Transverse;
    delete hDCAxy_deuterons_Toward;
    delete hDCAxy_deuterons_Away;
    delete hDCAxy_deuterons_Transverse;
    delete hDCAxy_antideuterons_Toward;
    delete hDCAxy_antideuterons_Away;
    delete hDCAxy_antideuterons_Transverse;
    delete hnsigmaTPC_deuterons_Syst;
    delete hnsigmaTPC_antideuterons_Syst;
    delete hnsigmaTOF_deuterons_Syst;
    delete hnsigmaTOF_antideuterons_Syst;
    delete hDCAxy_deuterons_Syst;
    delete hDCAxy_antideuterons_Syst;
    delete hnsigmaTPC_deuterons_rap;
    delete hnsigmaTPC_antideuterons_rap;
    delete hnsigmaTOF_deuterons_rap;
    delete hnsigmaTOF_antideuterons_rap;
    delete hnsigmaTPC_deuterons_rap_Syst;
    delete hnsigmaTPC_antideuterons_rap_Syst;
    delete hnsigmaTOF_deuterons_rap_Syst;
    delete hnsigmaTOF_antideuterons_rap_Syst;
    delete h_deuterons_Gen;
    delete h_antideuterons_Gen;
    delete hnsigmaTPC_deuterons_Rec;
    delete hnsigmaTPC_antideuterons_Rec;
    delete hnsigmaTOF_deuterons_Rec;
    delete hnsigmaTOF_antideuterons_Rec;
    delete hDCAxy_deuterons_Prim;
    delete hDCAxy_deuterons_Sec;
    delete hnsigmaTPC_deuterons_Rec_Syst;
    delete hnsigmaTPC_antideuterons_Rec_Syst;
    delete hnsigmaTOF_deuterons_Rec_Syst;
    delete hnsigmaTOF_antideuterons_Rec_Syst;
    delete hGeneratedDeuterons_vs_Rapidity;
    delete hGeneratedAntiDeuterons_vs_Rapidity;
    delete hReconstructedDeuterons_TPC_vs_Rapidity;
    delete hReconstructedAntiDeuterons_TPC_vs_Rapidity;
    delete hReconstructedDeuterons_TOF_vs_Rapidity;
    delete hReconstructedAntiDeuterons_TOF_vs_Rapidity;
    delete hnsigmaTPC_deuterons_Rec_rap_Syst;
    delete hnsigmaTPC_antideuterons_Rec_rap_Syst;
    delete hnsigmaTOF_deuterons_Rec_rap_Syst;
    delete hnsigmaTOF_antideuterons_Rec_rap_Syst;
}
//__________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskDeuteronsRT::UserCreateOutputObjects()  {
    
    //Create Output List
    fOutputList = new TList();
    fQAList     = new TList();
    fOutputList -> SetOwner();
    fQAList     -> SetOwner();
    
    
    //Pile-up Rejection
    fESDeventCuts.AddQAplotsToList(fQAList);
    fESDeventCuts.OverrideAutomaticTriggerSelection(fTriggerType);
    
    
    //Event Counter and Multiplicity Percentile
    hNumberOfEvents = new TH1F ("hNumberOfEvents","",20,0,20);
    hMultPercentile = new TH1F ("hMultPercentile","",200,0,100);
    hNumberOfEvents -> Sumw2();
    hMultPercentile -> Sumw2();
    fOutputList -> Add(hNumberOfEvents);
    fOutputList -> Add(hMultPercentile);

    
    //Multiplicity in Azimuthal Regions
    if ((!fIsMC) && fIsUEanalysis)  {
        hMultTransverse         = new TH1I ("hMultTransverse","",200,0,200);
        hMultToward             = new TH1I ("hMultToward","",200,0,200);
        hMultAway               = new TH1I ("hMultAway","",200,0,200);
        hRtDistribution         = new TH1F ("hRtDistribution","",200,0,20);
        hEventsWithLeadingTrack = new TH1F ("hEventsWithLeadingTrack","",10,0,10);
        hNumberOfDeuterons      = new TH1I ("hNumberOfDeuterons","",20,0,20);
        hMultTransverse         -> Sumw2();
        hMultToward             -> Sumw2();
        hMultAway               -> Sumw2();
        hRtDistribution         -> Sumw2();
        hEventsWithLeadingTrack -> Sumw2();
        hNumberOfDeuterons      -> Sumw2();
        fOutputList -> Add(hMultTransverse);
        fOutputList -> Add(hMultToward);
        fOutputList -> Add(hMultAway);
        fOutputList -> Add(hRtDistribution);
        fOutputList -> Add(hEventsWithLeadingTrack);
        fOutputList -> Add(hNumberOfDeuterons);
    }
    
    
    //Correlations between Transverse and Integrated Mult
    if ((!fIsMC) && fIsUEanalysis)  {
        
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
    Double_t pt_TPC[]   = { 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5 };
    Double_t pt_TOF[]   = { 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.5, 4.0, 5.0 };
    Double_t pt_DCA[]   = { 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0 };
    Double_t dca_xy[]   = { -1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.18, -0.16, -0.14, -0.12, -0.10, -0.08, -0.06, -0.04, -0.02, 0.0, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12, 0.14, 0.16, 0.18, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
    Double_t rapidity[] = { 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0 };
    
    //Arrays Dimensions
    Int_t nBins_Nch_Tr   = sizeof(Nch_Tr)/sizeof(Double_t)-1;
    Int_t nBins_pt_TPC   = sizeof(pt_TPC)/sizeof(Double_t)-1;
    Int_t nBins_pt_TOF   = sizeof(pt_TOF)/sizeof(Double_t)-1;
    Int_t nBins_pt_DCA   = sizeof(pt_DCA)/sizeof(Double_t)-1;
    Int_t nBins_dca_xy   = sizeof(dca_xy)/sizeof(Double_t)-1;
    Int_t nBins_rapidity = sizeof(rapidity)/sizeof(Double_t)-1;

    
    if ((!fIsMC) && fIsUEanalysis)  {
        
         //***************************  Nch_Transverse,  nsigma_{TPC},                 pt_TPC ***************************************************
         Int_t    bins_TPC[3] = {         nBins_Nch_Tr,           100,           nBins_pt_TPC };
         Double_t xmin_TPC[3] = {            Nch_Tr[0],          -10.0,              pt_TPC[0] };
         Double_t xmax_TPC[3] = { Nch_Tr[nBins_Nch_Tr],          +10.0,   pt_TPC[nBins_pt_TPC] };
         //***************************************************************************************************************************************
         
         hnsigmaTPC_deuterons_Toward         = new THnSparseF ("hnsigmaTPC_deuterons_Toward","",3, bins_TPC, xmin_TPC, xmax_TPC);
         hnsigmaTPC_deuterons_Away           = new THnSparseF ("hnsigmaTPC_deuterons_Away","",3, bins_TPC, xmin_TPC, xmax_TPC);
         hnsigmaTPC_deuterons_Transverse     = new THnSparseF ("hnsigmaTPC_deuterons_Transverse","",3, bins_TPC, xmin_TPC, xmax_TPC);
         hnsigmaTPC_antideuterons_Toward     = new THnSparseF ("hnsigmaTPC_antideuterons_Toward","",3, bins_TPC, xmin_TPC, xmax_TPC);
         hnsigmaTPC_antideuterons_Away       = new THnSparseF ("hnsigmaTPC_antideuterons_Away","",3, bins_TPC, xmin_TPC, xmax_TPC);
         hnsigmaTPC_antideuterons_Transverse = new THnSparseF ("hnsigmaTPC_antideuterons_Transverse","",3, bins_TPC, xmin_TPC, xmax_TPC);
         hnsigmaTPC_deuterons_Toward         -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
         hnsigmaTPC_deuterons_Away           -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
         hnsigmaTPC_deuterons_Transverse     -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
         hnsigmaTPC_antideuterons_Toward     -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
         hnsigmaTPC_antideuterons_Away       -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
         hnsigmaTPC_antideuterons_Transverse -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
         hnsigmaTPC_deuterons_Toward         -> GetAxis(2) -> Set(nBins_pt_TPC, pt_TPC);
         hnsigmaTPC_deuterons_Away           -> GetAxis(2) -> Set(nBins_pt_TPC, pt_TPC);
         hnsigmaTPC_deuterons_Transverse     -> GetAxis(2) -> Set(nBins_pt_TPC, pt_TPC);
         hnsigmaTPC_antideuterons_Toward     -> GetAxis(2) -> Set(nBins_pt_TPC, pt_TPC);
         hnsigmaTPC_antideuterons_Away       -> GetAxis(2) -> Set(nBins_pt_TPC, pt_TPC);
         hnsigmaTPC_antideuterons_Transverse -> GetAxis(2) -> Set(nBins_pt_TPC, pt_TPC);
         hnsigmaTPC_deuterons_Toward         -> Sumw2();
         hnsigmaTPC_deuterons_Away           -> Sumw2();
         hnsigmaTPC_deuterons_Transverse     -> Sumw2();
         hnsigmaTPC_antideuterons_Toward     -> Sumw2();
         hnsigmaTPC_antideuterons_Away       -> Sumw2();
         hnsigmaTPC_antideuterons_Transverse -> Sumw2();
         fOutputList -> Add (hnsigmaTPC_deuterons_Toward);
         fOutputList -> Add (hnsigmaTPC_deuterons_Away);
         fOutputList -> Add (hnsigmaTPC_deuterons_Transverse);
         fOutputList -> Add (hnsigmaTPC_antideuterons_Toward);
         fOutputList -> Add (hnsigmaTPC_antideuterons_Away);
         fOutputList -> Add (hnsigmaTPC_antideuterons_Transverse);

        
         //***************************  Nch_Transverse,  nsigma_{TOF},                 pt_TOF ***************************************************
         Int_t    bins_TOF[3] = {         nBins_Nch_Tr,           200,           nBins_pt_TOF };
         Double_t xmin_TOF[3] = {            Nch_Tr[0],         -10.0,              pt_TOF[0] };
         Double_t xmax_TOF[3] = { Nch_Tr[nBins_Nch_Tr],         +10.0,   pt_TOF[nBins_pt_TOF] };
         //***************************************************************************************************************************************
         
         hnsigmaTOF_deuterons_Toward         = new THnSparseF ("hnsigmaTOF_deuterons_Toward","",3, bins_TOF, xmin_TOF, xmax_TOF);
         hnsigmaTOF_deuterons_Away           = new THnSparseF ("hnsigmaTOF_deuterons_Away","",3, bins_TOF, xmin_TOF, xmax_TOF);
         hnsigmaTOF_deuterons_Transverse     = new THnSparseF ("hnsigmaTOF_deuterons_Transverse","",3, bins_TOF, xmin_TOF, xmax_TOF);
         hnsigmaTOF_antideuterons_Toward     = new THnSparseF ("hnsigmaTOF_antideuterons_Toward","",3, bins_TOF, xmin_TOF, xmax_TOF);
         hnsigmaTOF_antideuterons_Away       = new THnSparseF ("hnsigmaTOF_antideuterons_Away","",3, bins_TOF, xmin_TOF, xmax_TOF);
         hnsigmaTOF_antideuterons_Transverse = new THnSparseF ("hnsigmaTOF_antideuterons_Transverse","",3, bins_TOF, xmin_TOF, xmax_TOF);
         hnsigmaTOF_deuterons_Toward         -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
         hnsigmaTOF_deuterons_Away           -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
         hnsigmaTOF_deuterons_Transverse     -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
         hnsigmaTOF_antideuterons_Toward     -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
         hnsigmaTOF_antideuterons_Away       -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
         hnsigmaTOF_antideuterons_Transverse -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
         hnsigmaTOF_deuterons_Toward         -> GetAxis(2) -> Set(nBins_pt_TOF, pt_TOF);
         hnsigmaTOF_deuterons_Away           -> GetAxis(2) -> Set(nBins_pt_TOF, pt_TOF);
         hnsigmaTOF_deuterons_Transverse     -> GetAxis(2) -> Set(nBins_pt_TOF, pt_TOF);
         hnsigmaTOF_antideuterons_Toward     -> GetAxis(2) -> Set(nBins_pt_TOF, pt_TOF);
         hnsigmaTOF_antideuterons_Away       -> GetAxis(2) -> Set(nBins_pt_TOF, pt_TOF);
         hnsigmaTOF_antideuterons_Transverse -> GetAxis(2) -> Set(nBins_pt_TOF, pt_TOF);
         hnsigmaTOF_deuterons_Toward         -> Sumw2();
         hnsigmaTOF_deuterons_Away           -> Sumw2();
         hnsigmaTOF_deuterons_Transverse     -> Sumw2();
         hnsigmaTOF_antideuterons_Toward     -> Sumw2();
         hnsigmaTOF_antideuterons_Away       -> Sumw2();
         hnsigmaTOF_antideuterons_Transverse -> Sumw2();
         fOutputList -> Add (hnsigmaTOF_deuterons_Toward);
         fOutputList -> Add (hnsigmaTOF_deuterons_Away);
         fOutputList -> Add (hnsigmaTOF_deuterons_Transverse);
         fOutputList -> Add (hnsigmaTOF_antideuterons_Toward);
         fOutputList -> Add (hnsigmaTOF_antideuterons_Away);
         fOutputList -> Add (hnsigmaTOF_antideuterons_Transverse);

         
         //***************************  Nch_Transverse,      DCA_{xy},                 pt_DCA ***************************************************
         Int_t    bins_DCA[3] = {         nBins_Nch_Tr,  nBins_dca_xy,           nBins_pt_DCA };
         Double_t xmin_DCA[3] = {            Nch_Tr[0],          -1.0,              pt_DCA[0] };
         Double_t xmax_DCA[3] = { Nch_Tr[nBins_Nch_Tr],          +1.0,   pt_DCA[nBins_pt_DCA] };
         //***************************************************************************************************************************************

         hDCAxy_deuterons_Toward         = new THnSparseF ("hDCAxy_deuterons_Toward","",3, bins_DCA, xmin_DCA, xmax_DCA);
         hDCAxy_deuterons_Away           = new THnSparseF ("hDCAxy_deuterons_Away","",3, bins_DCA, xmin_DCA, xmax_DCA);
         hDCAxy_deuterons_Transverse     = new THnSparseF ("hDCAxy_deuterons_Transverse","",3, bins_DCA, xmin_DCA, xmax_DCA);
         hDCAxy_antideuterons_Toward     = new THnSparseF ("hDCAxy_antideuterons_Toward","",3, bins_DCA, xmin_DCA, xmax_DCA);
         hDCAxy_antideuterons_Away       = new THnSparseF ("hDCAxy_antideuterons_Away","",3, bins_DCA, xmin_DCA, xmax_DCA);
         hDCAxy_antideuterons_Transverse = new THnSparseF ("hDCAxy_antideuterons_Transverse","",3, bins_DCA, xmin_DCA, xmax_DCA);
         hDCAxy_deuterons_Toward         -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
         hDCAxy_deuterons_Away           -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
         hDCAxy_deuterons_Transverse     -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
         hDCAxy_antideuterons_Toward     -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
         hDCAxy_antideuterons_Away       -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
         hDCAxy_antideuterons_Transverse -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
         hDCAxy_deuterons_Toward         -> GetAxis(1) -> Set(nBins_dca_xy, dca_xy);
         hDCAxy_deuterons_Away           -> GetAxis(1) -> Set(nBins_dca_xy, dca_xy);
         hDCAxy_deuterons_Transverse     -> GetAxis(1) -> Set(nBins_dca_xy, dca_xy);
         hDCAxy_antideuterons_Toward     -> GetAxis(1) -> Set(nBins_dca_xy, dca_xy);
         hDCAxy_antideuterons_Away       -> GetAxis(1) -> Set(nBins_dca_xy, dca_xy);
         hDCAxy_antideuterons_Transverse -> GetAxis(1) -> Set(nBins_dca_xy, dca_xy);
         hDCAxy_deuterons_Toward         -> GetAxis(2) -> Set(nBins_pt_DCA, pt_DCA);
         hDCAxy_deuterons_Away           -> GetAxis(2) -> Set(nBins_pt_DCA, pt_DCA);
         hDCAxy_deuterons_Transverse     -> GetAxis(2) -> Set(nBins_pt_DCA, pt_DCA);
         hDCAxy_antideuterons_Toward     -> GetAxis(2) -> Set(nBins_pt_DCA, pt_DCA);
         hDCAxy_antideuterons_Away       -> GetAxis(2) -> Set(nBins_pt_DCA, pt_DCA);
         hDCAxy_antideuterons_Transverse -> GetAxis(2) -> Set(nBins_pt_DCA, pt_DCA);
         hDCAxy_deuterons_Toward         -> Sumw2();
         hDCAxy_deuterons_Away           -> Sumw2();
         hDCAxy_deuterons_Transverse     -> Sumw2();
         hDCAxy_antideuterons_Toward     -> Sumw2();
         hDCAxy_antideuterons_Away       -> Sumw2();
         hDCAxy_antideuterons_Transverse -> Sumw2();
         fOutputList -> Add(hDCAxy_deuterons_Toward);
         fOutputList -> Add(hDCAxy_deuterons_Away);
         fOutputList -> Add(hDCAxy_deuterons_Transverse);
         fOutputList -> Add(hDCAxy_antideuterons_Toward);
         fOutputList -> Add(hDCAxy_antideuterons_Away);
         fOutputList -> Add(hDCAxy_antideuterons_Transverse);

         
         //********************************  Nch_Transverse,  nsigma_{TPC},                 pt_TPC, isyst ****************************************
         Int_t    bins_TPC_syst[4] = {         nBins_Nch_Tr,           100,           nBins_pt_TPC,    50 };
         Double_t xmin_TPC_syst[4] = {            Nch_Tr[0],          -10.0,              pt_TPC[0],     0 };
         Double_t xmax_TPC_syst[4] = { Nch_Tr[nBins_Nch_Tr],          +10.0,   pt_TPC[nBins_pt_TPC],    50 };
         //***************************************************************************************************************************************
         
         hnsigmaTPC_deuterons_Syst     = new THnSparseF ("hnsigmaTPC_deuterons_Syst","",4, bins_TPC_syst, xmin_TPC_syst, xmax_TPC_syst);
         hnsigmaTPC_antideuterons_Syst = new THnSparseF ("hnsigmaTPC_antideuterons_Syst","",4, bins_TPC_syst, xmin_TPC_syst, xmax_TPC_syst);
         hnsigmaTPC_deuterons_Syst     -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
         hnsigmaTPC_antideuterons_Syst -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
         hnsigmaTPC_deuterons_Syst     -> GetAxis(2) -> Set(nBins_pt_TPC, pt_TPC);
         hnsigmaTPC_antideuterons_Syst -> GetAxis(2) -> Set(nBins_pt_TPC, pt_TPC);
         hnsigmaTPC_deuterons_Syst     -> Sumw2();
         hnsigmaTPC_antideuterons_Syst -> Sumw2();
         fOutputList -> Add (hnsigmaTPC_deuterons_Syst);
         fOutputList -> Add (hnsigmaTPC_antideuterons_Syst);

         
         //********************************  Nch_Transverse,  nsigma_{TOF},                  pt_TOF, isyst ****************************************
         Int_t    bins_TOF_syst[4] = {         nBins_Nch_Tr,           200,            nBins_pt_TOF,    50 };
         Double_t xmin_TOF_syst[4] = {            Nch_Tr[0],          -10.0,              pt_TOF[0],     0 };
         Double_t xmax_TOF_syst[4] = { Nch_Tr[nBins_Nch_Tr],          +10.0,   pt_TOF[nBins_pt_TOF],    50 };
         //***************************************************************************************************************************************
         
         hnsigmaTOF_deuterons_Syst     = new THnSparseF ("hnsigmaTOF_deuterons_Syst","",4, bins_TOF_syst, xmin_TOF_syst, xmax_TOF_syst);
         hnsigmaTOF_antideuterons_Syst = new THnSparseF ("hnsigmaTOF_antideuterons_Syst","",4, bins_TOF_syst, xmin_TOF_syst, xmax_TOF_syst);
         hnsigmaTOF_deuterons_Syst     -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
         hnsigmaTOF_antideuterons_Syst -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
         hnsigmaTOF_deuterons_Syst     -> GetAxis(2) -> Set(nBins_pt_TOF, pt_TOF);
         hnsigmaTOF_antideuterons_Syst -> GetAxis(2) -> Set(nBins_pt_TOF, pt_TOF);
         hnsigmaTOF_deuterons_Syst     -> Sumw2();
         hnsigmaTOF_antideuterons_Syst -> Sumw2();
         fOutputList -> Add (hnsigmaTOF_deuterons_Syst);
         fOutputList -> Add (hnsigmaTOF_antideuterons_Syst);
        
        
        //***************************       Nch_Transverse,      DCA_{xy},                 pt_DCA,   isyst *******************************************
        Int_t    bins_DCA_syst[4] = {         nBins_Nch_Tr,  nBins_dca_xy,           nBins_pt_DCA,    50 };
        Double_t xmin_DCA_syst[4] = {            Nch_Tr[0],          -1.0,              pt_DCA[0],     0 };
        Double_t xmax_DCA_syst[4] = { Nch_Tr[nBins_Nch_Tr],          +1.0,   pt_DCA[nBins_pt_DCA],    50 };
        //***************************************************************************************************************************************

        hDCAxy_deuterons_Syst     = new THnSparseF ("hDCAxy_deuterons_Syst","",4, bins_DCA_syst, xmin_DCA_syst, xmax_DCA_syst);
        hDCAxy_antideuterons_Syst = new THnSparseF ("hDCAxy_antideuterons_Syst","",4, bins_DCA_syst, xmin_DCA_syst, xmax_DCA_syst);
        hDCAxy_deuterons_Syst -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
        hDCAxy_deuterons_Syst -> GetAxis(1) -> Set(nBins_dca_xy, dca_xy);
        hDCAxy_deuterons_Syst -> GetAxis(2) -> Set(nBins_pt_DCA, pt_DCA);
        hDCAxy_antideuterons_Syst -> GetAxis(0) -> Set(nBins_Nch_Tr, Nch_Tr);
        hDCAxy_antideuterons_Syst -> GetAxis(1) -> Set(nBins_dca_xy, dca_xy);
        hDCAxy_antideuterons_Syst -> GetAxis(2) -> Set(nBins_pt_DCA, pt_DCA);
        hDCAxy_deuterons_Syst     -> Sumw2();
        hDCAxy_antideuterons_Syst -> Sumw2();
        fOutputList -> Add(hDCAxy_deuterons_Syst);
        fOutputList -> Add(hDCAxy_antideuterons_Syst);
    }
    
    
    if ((!fIsMC) && (!fIsUEanalysis))  {
        //***************************                Rapidity,  nsigma_{TPC/TOF},                     pt  ***************************************
        Int_t    bins_rap_TPC[3] = {           nBins_rapidity,               100,           nBins_pt_TPC };
        Double_t xmin_rap_TPC[3] = {              rapidity[0],              -5.0,              pt_TPC[0] };
        Double_t xmax_rap_TPC[3] = { rapidity[nBins_rapidity],              +5.0,   pt_TPC[nBins_pt_TPC] };
        Int_t    bins_rap_TOF[3] = {           nBins_rapidity,               200,           nBins_pt_TOF };
        Double_t xmin_rap_TOF[3] = {              rapidity[0],             -10.0,              pt_TOF[0] };
        Double_t xmax_rap_TOF[3] = { rapidity[nBins_rapidity],             +10.0,   pt_TOF[nBins_pt_TOF] };
        //***************************************************************************************************************************************
               
        hnsigmaTPC_deuterons_rap     = new THnSparseF ("hnsigmaTPC_deuterons_rap","",3, bins_rap_TPC, xmin_rap_TPC, xmax_rap_TPC);
        hnsigmaTPC_antideuterons_rap = new THnSparseF ("hnsigmaTPC_antideuterons_rap","",3, bins_rap_TPC, xmin_rap_TPC, xmax_rap_TPC);
        hnsigmaTOF_deuterons_rap     = new THnSparseF ("hnsigmaTOF_deuterons_rap","",3, bins_rap_TOF, xmin_rap_TOF, xmax_rap_TOF);
        hnsigmaTOF_antideuterons_rap = new THnSparseF ("hnsigmaTOF_antideuterons_rap","",3, bins_rap_TOF, xmin_rap_TOF, xmax_rap_TOF);
        hnsigmaTPC_deuterons_rap     -> GetAxis(0) -> Set(nBins_rapidity, rapidity);
        hnsigmaTPC_antideuterons_rap -> GetAxis(0) -> Set(nBins_rapidity, rapidity);
        hnsigmaTOF_deuterons_rap     -> GetAxis(0) -> Set(nBins_rapidity, rapidity);
        hnsigmaTOF_antideuterons_rap -> GetAxis(0) -> Set(nBins_rapidity, rapidity);
        hnsigmaTPC_deuterons_rap     -> GetAxis(2) -> Set(nBins_pt_TPC, pt_TPC);
        hnsigmaTPC_antideuterons_rap -> GetAxis(2) -> Set(nBins_pt_TPC, pt_TPC);
        hnsigmaTOF_deuterons_rap     -> GetAxis(2) -> Set(nBins_pt_TOF, pt_TOF);
        hnsigmaTOF_antideuterons_rap -> GetAxis(2) -> Set(nBins_pt_TOF, pt_TOF);
        hnsigmaTPC_deuterons_rap     -> Sumw2();
        hnsigmaTPC_antideuterons_rap -> Sumw2();
        hnsigmaTOF_deuterons_rap     -> Sumw2();
        hnsigmaTOF_antideuterons_rap -> Sumw2();
        fOutputList -> Add (hnsigmaTPC_deuterons_rap);
        fOutputList -> Add (hnsigmaTPC_antideuterons_rap);
        fOutputList -> Add (hnsigmaTOF_deuterons_rap);
        fOutputList -> Add (hnsigmaTOF_antideuterons_rap);

     
         //**************************************      Rapidity,      nsigma_{TPC},                 pt_TPC,   isyst ****************************************
         Int_t    bins_TPC_rap_syst[4] = {           nBins_rapidity,           100,           nBins_pt_TPC,     50 };
         Double_t xmin_TPC_rap_syst[4] = {              rapidity[0],          -10.0,              pt_TPC[0],     0 };
         Double_t xmax_TPC_rap_syst[4] = { rapidity[nBins_rapidity],          +10.0,   pt_TPC[nBins_pt_TPC],    50 };
         //***************************************************************************************************************************************
         
         hnsigmaTPC_deuterons_rap_Syst     = new THnSparseF ("hnsigmaTPC_deuterons_rap_Syst","",4, bins_TPC_rap_syst, xmin_TPC_rap_syst, xmax_TPC_rap_syst);
         hnsigmaTPC_antideuterons_rap_Syst = new THnSparseF ("hnsigmaTPC_antideuterons_rap_Syst","",4, bins_TPC_rap_syst, xmin_TPC_rap_syst, xmax_TPC_rap_syst);
         hnsigmaTPC_deuterons_rap_Syst     -> GetAxis(0) -> Set(nBins_rapidity, rapidity);
         hnsigmaTPC_antideuterons_rap_Syst -> GetAxis(0) -> Set(nBins_rapidity, rapidity);
         hnsigmaTPC_deuterons_rap_Syst     -> GetAxis(2) -> Set(nBins_pt_TPC, pt_TPC);
         hnsigmaTPC_antideuterons_rap_Syst -> GetAxis(2) -> Set(nBins_pt_TPC, pt_TPC);
         hnsigmaTPC_deuterons_rap_Syst     -> Sumw2();
         hnsigmaTPC_antideuterons_rap_Syst -> Sumw2();
         fOutputList -> Add (hnsigmaTPC_deuterons_rap_Syst);
         fOutputList -> Add (hnsigmaTPC_antideuterons_rap_Syst);

         
         //********************************            Rapidity,  nsigma_{TOF},                  pt_TOF, isyst ****************************************
         Int_t    bins_TOF_rap_syst[4] = {           nBins_rapidity,           200,            nBins_pt_TOF,    50 };
         Double_t xmin_TOF_rap_syst[4] = {              rapidity[0],          -10.0,              pt_TOF[0],     0 };
         Double_t xmax_TOF_rap_syst[4] = { rapidity[nBins_rapidity],          +10.0,   pt_TOF[nBins_pt_TOF],    50 };
         //***************************************************************************************************************************************
         
         hnsigmaTOF_deuterons_rap_Syst     = new THnSparseF ("hnsigmaTOF_deuterons_rap_Syst","",4, bins_TOF_rap_syst, xmin_TOF_rap_syst, xmax_TOF_rap_syst);
         hnsigmaTOF_antideuterons_rap_Syst = new THnSparseF ("hnsigmaTOF_antideuterons_rap_Syst","",4, bins_TOF_rap_syst, xmin_TOF_rap_syst, xmax_TOF_rap_syst);
         hnsigmaTOF_deuterons_rap_Syst     -> GetAxis(0) -> Set(nBins_rapidity, rapidity);
         hnsigmaTOF_antideuterons_rap_Syst -> GetAxis(0) -> Set(nBins_rapidity, rapidity);
         hnsigmaTOF_deuterons_rap_Syst     -> GetAxis(2) -> Set(nBins_pt_TOF, pt_TOF);
         hnsigmaTOF_antideuterons_rap_Syst -> GetAxis(2) -> Set(nBins_pt_TOF, pt_TOF);
         hnsigmaTOF_deuterons_rap_Syst     -> Sumw2();
         hnsigmaTOF_antideuterons_rap_Syst -> Sumw2();
         fOutputList -> Add (hnsigmaTOF_deuterons_rap_Syst);
         fOutputList -> Add (hnsigmaTOF_antideuterons_rap_Syst);
        
        
    }

    
    //Generated & Reconstructed p_{T} Spectra
    if (fIsMC)  {
        
        h_deuterons_Gen              = new TH1F ("h_deuterons_Gen","",900,0.5,5.0);
        h_antideuterons_Gen          = new TH1F ("h_antideuterons_Gen","",900,0.5,5.0);
        hnsigmaTPC_deuterons_Rec     = new TH2F ("hnsigmaTPC_deuterons_Rec","",900,0.5,5.0,100,-5,5);
        hnsigmaTPC_antideuterons_Rec = new TH2F ("hnsigmaTPC_antideuterons_Rec","",900,0.5,5.0,100,-5,5);
        hnsigmaTOF_deuterons_Rec     = new TH2F ("hnsigmaTOF_deuterons_Rec","",900,0.5,5.0,100,-5,5);
        hnsigmaTOF_antideuterons_Rec = new TH2F ("hnsigmaTOF_antideuterons_Rec","",900,0.5,5.0,100,-5,5);

        h_deuterons_Gen              -> Sumw2();
        h_antideuterons_Gen          -> Sumw2();
        hnsigmaTPC_deuterons_Rec     -> Sumw2();
        hnsigmaTPC_antideuterons_Rec -> Sumw2();
        hnsigmaTOF_deuterons_Rec     -> Sumw2();
        hnsigmaTOF_antideuterons_Rec -> Sumw2();
        fOutputList -> Add (h_deuterons_Gen);
        fOutputList -> Add (h_antideuterons_Gen);
        fOutputList -> Add (hnsigmaTPC_deuterons_Rec);
        fOutputList -> Add (hnsigmaTPC_antideuterons_Rec);
        fOutputList -> Add (hnsigmaTOF_deuterons_Rec);
        fOutputList -> Add (hnsigmaTOF_antideuterons_Rec);

        //Histograms for Syst. Uncertainties
        hnsigmaTPC_deuterons_Rec_Syst     = new TH2F ("hnsigmaTPC_deuterons_Rec_Syst","",900,0.5,5.0,50,0,50);
        hnsigmaTPC_antideuterons_Rec_Syst = new TH2F ("hnsigmaTPC_antideuterons_Rec_Syst","",900,0.5,5.0,50,0,50);
        hnsigmaTOF_deuterons_Rec_Syst     = new TH2F ("hnsigmaTOF_deuterons_Rec_Syst","",900,0.5,5.0,50,0,50);
        hnsigmaTOF_antideuterons_Rec_Syst = new TH2F ("hnsigmaTOF_antideuterons_Rec_Syst","",900,0.5,5.0,50,0,50);
        hnsigmaTPC_deuterons_Rec_Syst     -> Sumw2();
        hnsigmaTPC_antideuterons_Rec_Syst -> Sumw2();
        hnsigmaTOF_deuterons_Rec_Syst     -> Sumw2();
        hnsigmaTOF_antideuterons_Rec_Syst -> Sumw2();
        fOutputList -> Add(hnsigmaTPC_deuterons_Rec_Syst);
        fOutputList -> Add(hnsigmaTPC_antideuterons_Rec_Syst);
        fOutputList -> Add(hnsigmaTOF_deuterons_Rec_Syst);
        fOutputList -> Add(hnsigmaTOF_antideuterons_Rec_Syst);
        
        //DCA_{xy} Distributions
        hDCAxy_deuterons_Prim = new TH2F ("hDCAxy_deuterons_Prim","",nBins_pt_DCA,pt_DCA,nBins_dca_xy,dca_xy);
        hDCAxy_deuterons_Sec  = new TH2F ("hDCAxy_deuterons_Sec","",nBins_pt_DCA,pt_DCA,nBins_dca_xy,dca_xy);
        hDCAxy_deuterons_Prim -> Sumw2();
        hDCAxy_deuterons_Sec  -> Sumw2();
        fOutputList -> Add(hDCAxy_deuterons_Prim);
        fOutputList -> Add(hDCAxy_deuterons_Sec);
    }
    
    
    //Efficiency vs. Rapidity
    if (fIsMC && (!fIsUEanalysis))  {
        
        hGeneratedDeuterons_vs_Rapidity             = new TH2F ("hGeneratedDeuterons_vs_Rapidity","",900,0.5,5.0,10,0.0,1.0);
        hGeneratedAntiDeuterons_vs_Rapidity         = new TH2F ("hGeneratedAntiDeuterons_vs_Rapidity","",900,0.5,5.0,10,0.0,1.0);
        hReconstructedDeuterons_TPC_vs_Rapidity     = new TH2F ("hReconstructedDeuterons_TPC_vs_Rapidity","",900,0.5,5.0,10,0.0,1.0);
        hReconstructedAntiDeuterons_TPC_vs_Rapidity = new TH2F ("hReconstructedAntiDeuterons_TPC_vs_Rapidity","",900,0.5,5.0,10,0.0,1.0);
        hReconstructedDeuterons_TOF_vs_Rapidity     = new TH2F ("hReconstructedDeuterons_TOF_vs_Rapidity","",900,0.5,5.0,10,0.0,1.0);
        hReconstructedAntiDeuterons_TOF_vs_Rapidity = new TH2F ("hReconstructedAntiDeuterons_TOF_vs_Rapidity","",900,0.5,5.0,10,0.0,1.0);
        
        hGeneratedDeuterons_vs_Rapidity             -> Sumw2();
        hGeneratedAntiDeuterons_vs_Rapidity         -> Sumw2();
        hReconstructedDeuterons_TPC_vs_Rapidity     -> Sumw2();
        hReconstructedAntiDeuterons_TPC_vs_Rapidity -> Sumw2();
        hReconstructedDeuterons_TOF_vs_Rapidity     -> Sumw2();
        hReconstructedAntiDeuterons_TOF_vs_Rapidity -> Sumw2();
        
        fOutputList -> Add(hGeneratedDeuterons_vs_Rapidity);
        fOutputList -> Add(hGeneratedAntiDeuterons_vs_Rapidity);
        fOutputList -> Add(hReconstructedDeuterons_TPC_vs_Rapidity);
        fOutputList -> Add(hReconstructedAntiDeuterons_TPC_vs_Rapidity);
        fOutputList -> Add(hReconstructedDeuterons_TOF_vs_Rapidity);
        fOutputList -> Add(hReconstructedAntiDeuterons_TOF_vs_Rapidity);

        //Histograms for Syst. Uncertainties

         //**************************************      Rapidity,      nsigma_{TPC},                 pt_TPC,   isyst ****************************************
         Int_t    bins_TPC_rap_syst[4] = {           nBins_rapidity,           100,           nBins_pt_TPC,     50 };
         Double_t xmin_TPC_rap_syst[4] = {              rapidity[0],          -10.0,              pt_TPC[0],     0 };
         Double_t xmax_TPC_rap_syst[4] = { rapidity[nBins_rapidity],          +10.0,   pt_TPC[nBins_pt_TPC],    50 };
         //***************************************************************************************************************************************
         
         hnsigmaTPC_deuterons_Rec_rap_Syst     = new THnSparseF ("hnsigmaTPC_deuterons_Rec_rap_Syst","",4, bins_TPC_rap_syst, xmin_TPC_rap_syst, xmax_TPC_rap_syst);
         hnsigmaTPC_antideuterons_Rec_rap_Syst = new THnSparseF ("hnsigmaTPC_antideuterons_Rec_rap_Syst","",4, bins_TPC_rap_syst, xmin_TPC_rap_syst, xmax_TPC_rap_syst);
         hnsigmaTPC_deuterons_Rec_rap_Syst     -> GetAxis(0) -> Set(nBins_rapidity, rapidity);
         hnsigmaTPC_antideuterons_Rec_rap_Syst -> GetAxis(0) -> Set(nBins_rapidity, rapidity);
         hnsigmaTPC_deuterons_Rec_rap_Syst     -> GetAxis(2) -> Set(nBins_pt_TPC, pt_TPC);
         hnsigmaTPC_antideuterons_Rec_rap_Syst -> GetAxis(2) -> Set(nBins_pt_TPC, pt_TPC);
         hnsigmaTPC_deuterons_Rec_rap_Syst     -> Sumw2();
         hnsigmaTPC_antideuterons_Rec_rap_Syst -> Sumw2();
         fOutputList -> Add (hnsigmaTPC_deuterons_Rec_rap_Syst);
         fOutputList -> Add (hnsigmaTPC_antideuterons_Rec_rap_Syst);

         
         //********************************            Rapidity,      nsigma_{TOF},                  pt_TOF, isyst ****************************************
         Int_t    bins_TOF_rap_syst[4] = {           nBins_rapidity,           200,            nBins_pt_TOF,    50 };
         Double_t xmin_TOF_rap_syst[4] = {              rapidity[0],          -10.0,              pt_TOF[0],     0 };
         Double_t xmax_TOF_rap_syst[4] = { rapidity[nBins_rapidity],          +10.0,   pt_TOF[nBins_pt_TOF],    50 };
         //***************************************************************************************************************************************
         
         hnsigmaTOF_deuterons_Rec_rap_Syst     = new THnSparseF ("hnsigmaTOF_deuterons_Rec_rap_Syst","",4, bins_TOF_rap_syst, xmin_TOF_rap_syst, xmax_TOF_rap_syst);
         hnsigmaTOF_antideuterons_Rec_rap_Syst = new THnSparseF ("hnsigmaTOF_antideuterons_Rec_rap_Syst","",4, bins_TOF_rap_syst, xmin_TOF_rap_syst, xmax_TOF_rap_syst);
         hnsigmaTOF_deuterons_Rec_rap_Syst     -> GetAxis(0) -> Set(nBins_rapidity, rapidity);
         hnsigmaTOF_antideuterons_Rec_rap_Syst -> GetAxis(0) -> Set(nBins_rapidity, rapidity);
         hnsigmaTOF_deuterons_Rec_rap_Syst     -> GetAxis(2) -> Set(nBins_pt_TOF, pt_TOF);
         hnsigmaTOF_antideuterons_Rec_rap_Syst -> GetAxis(2) -> Set(nBins_pt_TOF, pt_TOF);
         hnsigmaTOF_deuterons_Rec_rap_Syst     -> Sumw2();
         hnsigmaTOF_antideuterons_Rec_rap_Syst -> Sumw2();
         fOutputList -> Add (hnsigmaTOF_deuterons_Rec_rap_Syst);
         fOutputList -> Add (hnsigmaTOF_antideuterons_Rec_rap_Syst);

        
    }
        
    
    //Track Selection Cuts: Deuterons
    fESDtrackCuts_Deuteron = new AliESDtrackCuts("fESDtrackCuts_Deuteron");
    fESDtrackCuts_Deuteron -> SetAcceptKinkDaughters(kFALSE);
    fESDtrackCuts_Deuteron -> SetRequireTPCRefit(kTRUE);
    fESDtrackCuts_Deuteron -> SetRequireITSRefit(kTRUE);
    fESDtrackCuts_Deuteron -> SetMinNCrossedRowsTPC(60);
    fESDtrackCuts_Deuteron -> SetMinRatioCrossedRowsOverFindableClustersTPC(0.6);
    fESDtrackCuts_Deuteron -> SetMinNClustersITS(2);
    fESDtrackCuts_Deuteron -> SetClusterRequirementITS (AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
    fESDtrackCuts_Deuteron -> SetMaxChi2PerClusterITS(36);
    fESDtrackCuts_Deuteron -> SetMaxChi2PerClusterTPC(5);
    fESDtrackCuts_Deuteron -> SetEtaRange(-2.0,2.0);
    fESDtrackCuts_Deuteron -> SetMaxDCAToVertexXY(2);
    fESDtrackCuts_Deuteron -> SetMaxDCAToVertexZ(2);
    fESDtrackCuts_Deuteron -> SetDCAToVertex2D(kFALSE);
    fESDtrackCuts_Deuteron -> SetRequireSigmaToVertex(kFALSE);
    
    
    
    //Track Selection Cuts: Charged Tracks
    if ((!fIsMC) && fIsUEanalysis)  {
        
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
        //fESDtrackCuts_LeadingTrack -> SetMaxDCAToVertexXY(2.4);
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
//__________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskDeuteronsRT::UserExec(Option_t *)  {
    
    //Get Input Event
    if ( !GetESDEvent ()) return;
    if (fIsMC && (!GetMCEvent ())) return;
    
    //Load PID Response
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();
    
    //Process Real or Simulated Event
    if ((!fIsMC) && fIsUEanalysis)    ProcessRealEvent ();
    if ((!fIsMC) && (!fIsUEanalysis)) ProcessRealEventRapidityDependence();
    if ( fIsMC) ProcessSimEvent ();
    
    //Post Output Data
    PostData(1, fOutputList);
}
//__________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDeuteronsRT::GetESDEvent ()  {
    
    
   //Get Input Event
   fESDevent = dynamic_cast <AliESDEvent*>(InputEvent());
   if (!fESDevent) return kFALSE;
   hNumberOfEvents -> Fill(0.5);
   
       
   //Standard Event Cuts
   if (!fESDeventCuts.AcceptEvent(fESDevent)) {
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
Bool_t AliAnalysisTaskDeuteronsRT::GetMCEvent ()  {
    
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
void AliAnalysisTaskDeuteronsRT::ProcessRealEvent ()  {
    
    
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
    
    //Deuteron Candidates ID
    vector<Int_t> deuteron_ID;
    
    //Loop over Reconstructed Tracks
    for (Int_t i=0 ; i<fESDevent->GetNumberOfTracks() ; i++)  {

        //Get Reconstructed Track
        AliESDtrack *track = (AliESDtrack*) fESDevent->GetTrack(i);
        if (!track) continue;
        
        //Store Deuteron Candidates
        if (IsDeuteronCandidate (track))  { deuteron_ID.push_back(i); }
        
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
    
    //Number of Deuteron Candidates
    Int_t nDeuterons = (Int_t)deuteron_ID.size();
    hNumberOfDeuterons -> Fill (nDeuterons);
    
    //Loop over Deuteron Candidates
    for (Int_t i=0 ; i<nDeuterons ; i++)  {

        //Get Reconstructed Track
        AliESDtrack *track = (AliESDtrack*) fESDevent->GetTrack(deuteron_ID[i]);
        if (!track) continue;
    
        //Fill Histograms: Standard Cuts
        FillHistograms_StandardCuts (mult_Transverse,leading_track_ID,track);
        
        //Fill Histograms: Systematic Uncertainties
        for (Int_t isyst=0 ; isyst<50 ; isyst++)
            FillHistograms_Systematics (mult_Transverse,leading_track_ID,track,isyst);
    }
}
//__________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskDeuteronsRT::ProcessRealEventRapidityDependence ()  {
    
    //Loop over Reconstructed Tracks
    for (Int_t i=0 ; i<fESDevent->GetNumberOfTracks() ; i++)  {

        //Get Reconstructed Track
        AliESDtrack *track = (AliESDtrack*) fESDevent->GetTrack(i);
        if (!track) continue;

        //Basic Track Quality Cuts
        if (!PassedBasicTrackQualityCuts_NoRapidityCut (track)) continue;
    
        //Fill Histograms: Rapidity Dependence
        FillHistograms_RapidityDependence (track);

        //Fill Histograms: Systematic Uncertainties
        for (Int_t isyst=0 ; isyst<50 ; isyst++)
            FillHistograms_Rapidity_Systematics (track,isyst);

    }
}
//__________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskDeuteronsRT::ProcessSimEvent ()  {
    
    
    //Loop over Generated Particles
    for (Int_t i=0; i<fMCEvent->GetNumberOfTracks(); i++)  {
           
        //MC Particle Selection
        AliMCParticle *particle = (AliMCParticle*) fMCEvent->GetTrack(i);
        if (!particle) continue;
        if (!particle->IsPhysicalPrimary()) continue;
        if ( TMath::Abs(particle->PdgCode()) != 1000010020 ) continue;
        if (AliAnalysisUtils::IsParticleFromOutOfBunchPileupCollision(i,fMCEvent)) continue;
        
        //Rapidity Selection
        Double_t m = AliPID::ParticleMass(AliPID::kDeuteron);
        Double_t E = TMath::Sqrt(m*m + particle->P()*particle->P());
        TLorentzVector P (particle->Px(),particle->Py(),particle->Pz(),E);
        Double_t y_lab = P.Rapidity();
        Double_t y_cms(0);
        
        //Rapidity-Dependent Efficiency
        if (!fIsUEanalysis)  {
            
            y_cms = y_lab;
            if (particle->PdgCode()==1000010020)  hGeneratedDeuterons_vs_Rapidity     -> Fill(particle->Pt(),TMath::Abs(y_cms));
            if (particle->PdgCode()==-1000010020) hGeneratedAntiDeuterons_vs_Rapidity -> Fill(particle->Pt(),TMath::Abs(y_cms));
        }
        
        //Rapidity Selection in p-Pb Collisions
        if (fIspPb)  {
            y_cms = y_lab-0.465;
            if (y_cms<-1.0) continue;
            if (y_cms> 0.0) continue;
        }
        
        //Rapidity Selection in pp Collisions
        if (!fIspPb) {
            y_cms = y_lab;
            if (TMath::Abs(y_cms)>0.5) continue;
        }

        //Fill Generated p_{T} Spectra
        if ( particle->PdgCode() == +1000010020 ) h_deuterons_Gen     -> Fill(particle->Pt());
        if ( particle->PdgCode() == -1000010020 ) h_antideuterons_Gen -> Fill(particle->Pt());
    }
       
    //Loop over Reconstructed Tracks
    if (fIsUEanalysis)  {
        for (Int_t i=0 ; i<fESDevent->GetNumberOfTracks() ; i++)  {

            //Get Reconstructed Track
            AliESDtrack *track = (AliESDtrack*) fESDevent->GetTrack(i);
            if (!track) continue;

            //Basic Track Quality Cuts
            if (!PassedBasicTrackQualityCuts (track)) continue;
            
            //MC Particle
            AliMCParticle *particle = (AliMCParticle*) fMCEvent->GetTrack(TMath::Abs(track->GetLabel()));
            if (!particle) continue;
            if ( TMath::Abs(particle->PdgCode()) != 1000010020 ) continue;
            
            //Fill Histograms: Standard Cuts
            FillHistograms_StandardCuts_Sim (track);
            
            //Fill Histograms: Systematic Uncertainties
            for (Int_t isyst=0 ; isyst<50 ; isyst++)
                FillHistograms_Systematics_Sim (track,isyst);
        }
    }
    
    
    //Loop over Reconstructed Tracks
    if (!fIsUEanalysis)  {
        for (Int_t i=0 ; i<fESDevent->GetNumberOfTracks() ; i++)  {

            //Get Reconstructed Track
            AliESDtrack *track = (AliESDtrack*) fESDevent->GetTrack(i);
            if (!track) continue;

            //Basic Track Quality Cuts
            if (!PassedBasicTrackQualityCuts_NoRapidityCut (track)) continue;
            
            //Rapidity
            Double_t mass = AliPID::ParticleMass(AliPID::kDeuteron);
            Double_t p  = track->P();
            Double_t pz = track->Pz();
            Double_t E = TMath::Sqrt(mass*mass + p*p);
            if (E == TMath::Abs(pz)) continue;
            Double_t y_lab = 0.5*TMath::Log((E+pz)/(E-pz));
            Double_t y_cms = y_lab;
            
            //MC Particle
            AliMCParticle *particle = (AliMCParticle*) fMCEvent->GetTrack(TMath::Abs(track->GetLabel()));
            if (!particle) continue;
            if ( TMath::Abs(particle->PdgCode()) != 1000010020 ) continue;
            
            //Fill Histograms: Standard Cuts
            FillHistograms_StandardCuts_Sim (track);
            
            //Fill Histograms: Systematic Uncertainties
            for (Int_t isyst=0 ; isyst<50 ; isyst++)
                FillHistograms_Rapidity_Systematics_Sim (track,isyst);                
        }
    }
    
}
//__________________________________________________________________________________________________________________________________________________
Int_t AliAnalysisTaskDeuteronsRT::GetLeadingTrack ()  {
    
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
Bool_t AliAnalysisTaskDeuteronsRT::IsTrackInTransverseRegion (AliESDtrack *track, Int_t leading_track_ID)  {
    
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
Bool_t AliAnalysisTaskDeuteronsRT::IsTrackInTowardRegion (AliESDtrack *track, Int_t leading_track_ID)  {
    
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
Bool_t AliAnalysisTaskDeuteronsRT::IsTrackInAwayRegion (AliESDtrack *track, Int_t leading_track_ID)  {
    
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
Bool_t AliAnalysisTaskDeuteronsRT::IsCleanDeuteron (AliESDtrack *track)  {
    
    //Initialization
    Bool_t isCleanDeuteron = kFALSE;
    
    //Variables
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kDeuteron);
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kDeuteron);
    Double_t pt        = track->Pt();
    Bool_t   hasTOFhit = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);
    Double_t length    = track->GetIntegratedLength();
    
    //Selections
    if (pt<1.0 && TMath::Abs(nsigmaTPC)<3.0) isCleanDeuteron = kTRUE;
    if (pt>1.0 && hasTOFhit && length>350.0 && TMath::Abs(nsigmaTPC)<3.0 && TMath::Abs(nsigmaTOF)<3.0) isCleanDeuteron = kTRUE;

    return isCleanDeuteron;
    
}
//__________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDeuteronsRT::IsDeuteronCandidate (AliESDtrack *track)  {
    
    //Initialization
    Bool_t isDeuteronCandidate = kFALSE;

    //Basic Track Selection
    if (!PassedBasicTrackQualityCuts (track)) return isDeuteronCandidate;

    //Variables
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kDeuteron);
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kDeuteron);
    Double_t pt        = track->Pt();
    Bool_t   hasTOFhit = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);
    Double_t length    = track->GetIntegratedLength();

    
    //TPC Pre-selection
    if (TMath::Abs(nsigmaTPC)>10.0) return isDeuteronCandidate;
    
    //Selection on TOF hit & Track Length
    if (pt>1.5 && (!hasTOFhit)) return isDeuteronCandidate;
    if (pt>1.5 && length<350.0) return isDeuteronCandidate;
    if (pt>1.5 && TMath::Abs(nsigmaTOF)>10.0) return isDeuteronCandidate;

    
    isDeuteronCandidate = kTRUE;
    return isDeuteronCandidate;
}
//__________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDeuteronsRT::PassedTrackQualityCuts_LeadingTrack (AliESDtrack *track)  {
    
    //Initialization
    Bool_t passedTrkSelection=(kFALSE);
    if ( !fESDtrackCuts_LeadingTrack->AcceptTrack (track) ) return passedTrkSelection;
    
    passedTrkSelection = kTRUE;
    return passedTrkSelection;
}
//__________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDeuteronsRT::PassedTrackQualityCuts_TransverseMult (AliESDtrack *track)  {
    
    //Initialization
    Bool_t passedTrkSelection=(kFALSE);
    if ( !fESDtrackCuts_TransverseMult->AcceptTrack (track) ) return passedTrkSelection;
    
    passedTrkSelection = kTRUE;
    return passedTrkSelection;
}
//__________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDeuteronsRT::PassedBasicTrackQualityCuts (AliESDtrack *track)  {
    
    //Initialization
    Bool_t passedTrkSelection=(kFALSE);
    
    if ( track->GetTPCsignalN() < 50 )                  return passedTrkSelection;
    if ( !fESDtrackCuts_Deuteron->AcceptTrack (track) ) return passedTrkSelection;
    if ( TMath::Abs(track->Eta()) > 0.8 )               return passedTrkSelection;

    
    //Rapidity Cut
    Double_t mass = AliPID::ParticleMass(AliPID::kDeuteron);
    Double_t p  = track->P();
    Double_t pz = track->Pz();
    Double_t E = TMath::Sqrt(mass*mass + p*p);
    if (E == TMath::Abs(pz)) return passedTrkSelection;
    Double_t y_lab = 0.5*TMath::Log((E+pz)/(E-pz));
    Double_t y_cms(0);

    //Rapidity Selection in p-Pb Collisions
    if (fIspPb)  {
        y_cms = y_lab-0.465;
        if (y_cms<-1.0) return passedTrkSelection;
        if (y_cms> 0.0) return passedTrkSelection;
    }
    
    //Rapidity Selection in pp Collisions
    if (!fIspPb) {
        y_cms = y_lab;
        if (TMath::Abs(y_cms)>0.5) return passedTrkSelection;
    }
    
    passedTrkSelection = kTRUE;
    return passedTrkSelection;
}
//__________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDeuteronsRT::PassedBasicTrackQualityCuts_NoRapidityCut (AliESDtrack *track)  {
    
    //Initialization
    Bool_t passedTrkSelection=(kFALSE);
    
    if ( track->GetTPCsignalN() < 50 )         return passedTrkSelection;
    if ( !fESDtrackCuts_Deuteron->AcceptTrack (track) ) return passedTrkSelection;

    
    passedTrkSelection = kTRUE;
    return passedTrkSelection;
}
//__________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskDeuteronsRT::FillHistograms_StandardCuts (Int_t mult_Transverse, Int_t leading_track_ID, AliESDtrack *track)  {
    
    //Variables
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kDeuteron);
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kDeuteron);
    Double_t pt        = track->Pt();
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
        
    //Vectors
    Double_t xTPC[3] = { Rt, nsigmaTPC, pt };
    Double_t xTOF[3] = { Rt, nsigmaTOF, pt };
    Double_t xDCA[3] = { Rt, DCAxy, pt };
    
    //DCA_{xy} Histograms
    if (IsCleanDeuteron (track))  {
            
        if (charge>0)  {
            if (IsTrackInTowardRegion(track,leading_track_ID))     hDCAxy_deuterons_Toward         -> Fill (xDCA);
            if (IsTrackInAwayRegion(track,leading_track_ID))       hDCAxy_deuterons_Away           -> Fill (xDCA);
            if (IsTrackInTransverseRegion(track,leading_track_ID)) hDCAxy_deuterons_Transverse     -> Fill (xDCA);
        }
            
        if (charge<0)  {
            if (IsTrackInTowardRegion(track,leading_track_ID))     hDCAxy_antideuterons_Toward     -> Fill (xDCA);
            if (IsTrackInAwayRegion(track,leading_track_ID))       hDCAxy_antideuterons_Away       -> Fill (xDCA);
            if (IsTrackInTransverseRegion(track,leading_track_ID)) hDCAxy_antideuterons_Transverse -> Fill (xDCA);
        }
    }
    
    //DCA_{xy} Cut
    if (TMath::Abs(DCAxy)>0.1) return;
    
    //TPC-Only Analysis
    if (pt<1.5)  {
        
        if (charge>0)  {
            if (IsTrackInTowardRegion(track,leading_track_ID))     hnsigmaTPC_deuterons_Toward     -> Fill (xTPC);
            if (IsTrackInAwayRegion(track,leading_track_ID))       hnsigmaTPC_deuterons_Away       -> Fill (xTPC);
            if (IsTrackInTransverseRegion(track,leading_track_ID)) hnsigmaTPC_deuterons_Transverse -> Fill (xTPC);
        }
            
        if (charge<0)  {
            if (IsTrackInTowardRegion(track,leading_track_ID))     hnsigmaTPC_antideuterons_Toward     -> Fill (xTPC);
            if (IsTrackInAwayRegion(track,leading_track_ID))       hnsigmaTPC_antideuterons_Away       -> Fill (xTPC);
            if (IsTrackInTransverseRegion(track,leading_track_ID)) hnsigmaTPC_antideuterons_Transverse -> Fill (xTPC);
        }
    }
        
        
    //TOF Analysis
    if (!hasTOFhit) return;
    if (length<350.0) return;
    if (TMath::Abs(nsigmaTPC)>3.0) return;
    if (pt<1.0) return;
        
            
    if (charge>0)  {
        if (IsTrackInTowardRegion(track,leading_track_ID))     hnsigmaTOF_deuterons_Toward     -> Fill (xTOF);
        if (IsTrackInAwayRegion(track,leading_track_ID))       hnsigmaTOF_deuterons_Away       -> Fill (xTOF);
        if (IsTrackInTransverseRegion(track,leading_track_ID)) hnsigmaTOF_deuterons_Transverse -> Fill (xTOF);
    }
            
    if (charge<0)  {
        if (IsTrackInTowardRegion(track,leading_track_ID))     hnsigmaTOF_antideuterons_Toward     -> Fill (xTOF);
        if (IsTrackInAwayRegion(track,leading_track_ID))       hnsigmaTOF_antideuterons_Away       -> Fill (xTOF);
        if (IsTrackInTransverseRegion(track,leading_track_ID)) hnsigmaTOF_antideuterons_Transverse -> Fill (xTOF);
    }
}
//__________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskDeuteronsRT::GetRapidity (AliESDtrack *track)  {
    
    //Initialization
    Double_t y(-999);
    
    //Rapidity Calculation
    Double_t mass = AliPID::ParticleMass(AliPID::kDeuteron);
    Double_t p  = track->P();
    Double_t pz = track->Pz();
    Double_t E = TMath::Sqrt(mass*mass + p*p);
    if (E == TMath::Abs(pz)) return -999;
    y = 0.5*TMath::Log((E+pz)/(E-pz));
    
    return y;
}
//__________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskDeuteronsRT::FillHistograms_RapidityDependence (AliESDtrack *track)  {
    
    
    //Variables
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kDeuteron);
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kDeuteron);
    Double_t pt        = track->Pt();
    Double_t charge    = track->Charge();
    Double_t DCAxy     = GetTransverseDCA (track);
    Bool_t   hasTOFhit = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);
    Double_t length    = track->GetIntegratedLength();
    Double_t rap       = TMath::Abs(GetRapidity (track));
    
    
    //Track Quality Cuts
    if (!PassedTrackQualityCuts_Syst(track,0)) return;
       
        
    //Vectors
    Double_t xTPC[3] = { rap, nsigmaTPC, pt };
    Double_t xTOF[3] = { rap, nsigmaTOF, pt };

    
    //DCA_{xy} Cut
    if (TMath::Abs(DCAxy)>0.1) return;

    
    //TPC-Only Analysis
    if (pt<1.5)  {
        
        if (charge>0)  hnsigmaTPC_deuterons_rap     -> Fill (xTPC);
        if (charge<0)  hnsigmaTPC_antideuterons_rap -> Fill (xTPC);
    }
        
        
    //TOF Analysis
    if (!hasTOFhit) return;
    if (length<350.0) return;
    if (TMath::Abs(nsigmaTPC)>3.0) return;
    if (pt<1.0) return;
        
    if (charge>0)  hnsigmaTOF_deuterons_rap     -> Fill (xTOF);
    if (charge<0)  hnsigmaTOF_antideuterons_rap -> Fill (xTOF);

}
//__________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskDeuteronsRT::FillHistograms_StandardCuts_Sim (AliESDtrack *track)  {
    
    //Variables
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kDeuteron);
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kDeuteron);
    Double_t pt        = track->Pt();
    Double_t charge    = track->Charge();
    Double_t DCAxy     = GetTransverseDCA (track);
    Double_t DCAz      = GetLongitudinalDCA (track);
    Bool_t   hasTOFhit = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);
    Double_t length    = track->GetIntegratedLength();
    Double_t mass      = AliPID::ParticleMass(AliPID::kDeuteron);
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

    //DCA_{xy} Distributions
    if (particle->IsPhysicalPrimary())       hDCAxy_deuterons_Prim -> Fill (pt,DCAxy);
    if (particle->IsSecondaryFromMaterial()) hDCAxy_deuterons_Sec  -> Fill (pt,DCAxy);

    //DCA_{xy} Cut
    if (TMath::Abs(DCAxy)>0.1) return;
    if (!particle->IsPhysicalPrimary()) return;
    
    //TPC-Only Analysis
    if (pt<1.5)  {
       
        if (charge>0)  hnsigmaTPC_deuterons_Rec     -> Fill (pt,nsigmaTPC);
        if (charge<0)  hnsigmaTPC_antideuterons_Rec -> Fill (pt,nsigmaTPC);
        if (!fIsUEanalysis)  {
            if (charge>0 && TMath::Abs(nsigmaTPC)<3.0)  hReconstructedDeuterons_TPC_vs_Rapidity     -> Fill (pt,TMath::Abs(y_cms));
            if (charge<0 && TMath::Abs(nsigmaTPC)<3.0)  hReconstructedAntiDeuterons_TPC_vs_Rapidity -> Fill (pt,TMath::Abs(y_cms));
        }
    }
        
    //TOF Analysis
    if (!hasTOFhit) return;
    if (length<350.0) return;
    if (TMath::Abs(nsigmaTPC)>3.0) return;
    if (pt<1.0) return;
            
    if (charge>0) hnsigmaTOF_deuterons_Rec     -> Fill(pt,nsigmaTOF);
    if (charge<0) hnsigmaTOF_antideuterons_Rec -> Fill(pt,nsigmaTOF);
    if (!fIsUEanalysis)  {
        if (charge>0 && TMath::Abs(nsigmaTOF)<3.0)  hReconstructedDeuterons_TOF_vs_Rapidity     -> Fill (pt,TMath::Abs(y_cms));
        if (charge<0 && TMath::Abs(nsigmaTOF)<3.0)  hReconstructedAntiDeuterons_TOF_vs_Rapidity -> Fill (pt,TMath::Abs(y_cms));
    }
}
//__________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskDeuteronsRT::FillHistograms_Systematics  (Int_t mult_Transverse, Int_t leading_track_ID, AliESDtrack *track, Int_t isyst)  {
    
    
    Double_t DCAxy_max[]={0.0800,0.1200,0.1600,0.1400,0.1000,0.1000,0.1200,0.1400,0.1000,0.1200,0.1400,0.1400,0.1200,0.1200,0.1200,0.1600,0.1000,0.1000,0.1200,0.1000,0.1000,0.1200,0.1000,0.1200,0.1000,0.1400,0.1400,0.1400,0.0800,0.1600,0.1400,0.1000,0.1200,0.1400,0.1400,0.0800,0.1200,0.1200,0.1600,0.1000,0.1000,0.1600,0.1200,0.1400,0.1200,0.1600,0.1000,0.1200,0.1000,0.1400};

    Double_t DCAz_max[]={0.9113,0.9840,1.0381,0.8747,1.4862,0.9562,1.3854,1.2460,1.0332,1.3982,0.8442,0.9277,1.3835,0.9963,1.2933,1.3615,1.3876,1.4857,1.2432,1.0957,1.1053,0.9439,0.9174,1.4443,0.9351,1.3539,0.8032,0.9003,1.2470,1.3777,1.3049,1.4598,1.0615,1.3109,0.8874,0.9497,0.8719,0.8165,1.0398,1.3376,0.8662,1.1053,1.0806,1.2371,0.9518,1.4904,1.0039,1.0135,0.8318,1.3350};
    
    
    //Variables
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kDeuteron);
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kDeuteron);
    Double_t pt        = track->Pt();
    Double_t charge    = track->Charge();
    Double_t DCAxy     = GetTransverseDCA (track);
    Double_t DCAz      = GetLongitudinalDCA (track);
    Bool_t   hasTOFhit = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);
    Double_t length    = track->GetIntegratedLength();
    Double_t Rt        = static_cast<Double_t>(mult_Transverse)/fAverage_Nch_Transv;

    
    //Track Quality Cuts
    if (!PassedTrackQualityCuts_Syst(track,isyst)) return;
        
    //DCA_{z} Cut
    if (TMath::Abs(DCAz)>DCAz_max[isyst]) return;
    
    
    //Vectors
    Double_t xTPC[4] = { Rt, nsigmaTPC, pt, static_cast<Double_t>(isyst) };
    Double_t xTOF[4] = { Rt, nsigmaTOF, pt, static_cast<Double_t>(isyst) };
    Double_t xDCA[4] = { Rt, DCAxy, pt, static_cast<Double_t>(isyst) };
    
    
    //DCA_{xy} Histograms
    if (IsCleanDeuteron (track))  {
              
        if (charge>0) hDCAxy_deuterons_Syst     -> Fill (xDCA);
        if (charge<0) hDCAxy_antideuterons_Syst -> Fill (xDCA);
    }

        
    //DCA cuts
    if (TMath::Abs(DCAxy)>DCAxy_max[isyst]) return;

    
    //TPC-Only Analysis
    if (pt<1.5)  {
        
        if (charge>0) hnsigmaTPC_deuterons_Syst     -> Fill (xTPC);
        if (charge<0) hnsigmaTPC_antideuterons_Syst -> Fill (xTPC);
    }
        
        
    //TOF Analysis
    if (!hasTOFhit) return;
    if (length<350.0) return;
    if (TMath::Abs(nsigmaTPC)>3.0) return;
    if (pt<1.0) return;
        
            
    if (charge>0) hnsigmaTOF_deuterons_Syst     -> Fill (xTOF);
    if (charge<0) hnsigmaTOF_antideuterons_Syst -> Fill (xTOF);
    
}
//__________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskDeuteronsRT::FillHistograms_Systematics_Sim  (AliESDtrack *track, Int_t isyst)  {
    
    Double_t DCAxy_max[]={0.0800,0.1200,0.1600,0.1400,0.1000,0.1000,0.1200,0.1400,0.1000,0.1200,0.1400,0.1400,0.1200,0.1200,0.1200,0.1600,0.1000,0.1000,0.1200,0.1000,0.1000,0.1200,0.1000,0.1200,0.1000,0.1400,0.1400,0.1400,0.0800,0.1600,0.1400,0.1000,0.1200,0.1400,0.1400,0.0800,0.1200,0.1200,0.1600,0.1000,0.1000,0.1600,0.1200,0.1400,0.1200,0.1600,0.1000,0.1200,0.1000,0.1400};

    Double_t DCAz_max[]={0.9113,0.9840,1.0381,0.8747,1.4862,0.9562,1.3854,1.2460,1.0332,1.3982,0.8442,0.9277,1.3835,0.9963,1.2933,1.3615,1.3876,1.4857,1.2432,1.0957,1.1053,0.9439,0.9174,1.4443,0.9351,1.3539,0.8032,0.9003,1.2470,1.3777,1.3049,1.4598,1.0615,1.3109,0.8874,0.9497,0.8719,0.8165,1.0398,1.3376,0.8662,1.1053,1.0806,1.2371,0.9518,1.4904,1.0039,1.0135,0.8318,1.3350};
    
    //Variables
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kDeuteron);
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kDeuteron);
    Double_t pt        = track->Pt();
    Double_t charge    = track->Charge();
    Double_t DCAxy     = GetTransverseDCA (track);
    Double_t DCAz      = GetLongitudinalDCA(track);
    Bool_t   hasTOFhit = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);
    Double_t length    = track->GetIntegratedLength();

    //MC Particle
    AliMCParticle *particle = (AliMCParticle*) fMCEvent->GetTrack(TMath::Abs(track->GetLabel()));
       
    //Track Track Quality Cuts
    if (!PassedTrackQualityCuts_Syst(track,isyst)) return;
       
        
    //DCA_{xy} cut
    if (TMath::Abs(DCAz)>DCAz_max[isyst])   return;
    if (TMath::Abs(DCAxy)>DCAxy_max[isyst]) return;
    if (!particle->IsPhysicalPrimary())     return;

    
    //TPC-Only Analysis
    if (pt<1.5 && TMath::Abs(nsigmaTPC)<3.0)  {
      
        if (charge>0) hnsigmaTPC_deuterons_Rec_Syst     -> Fill (pt,isyst);
        if (charge<0) hnsigmaTPC_antideuterons_Rec_Syst -> Fill (pt,isyst);
    }
        
        
    //TOF Analysis
    if (!hasTOFhit) return;
    if (length<350.0) return;
    if (TMath::Abs(nsigmaTPC)>3.0) return;
    if (nsigmaTOF<-3.0) return;
    if (nsigmaTOF>+3.5) return;
    if (pt<1.0) return;
       
    if (charge>0) hnsigmaTOF_deuterons_Rec_Syst     -> Fill (pt,isyst);
    if (charge<0) hnsigmaTOF_antideuterons_Rec_Syst -> Fill (pt,isyst);
    
}
//__________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskDeuteronsRT::FillHistograms_Rapidity_Systematics  (AliESDtrack *track, Int_t isyst)  {
    
    //Variables
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kDeuteron);
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kDeuteron);
    Double_t pt        = track->Pt();
    Double_t charge    = track->Charge();
    Double_t DCAxy     = GetTransverseDCA (track);
    Double_t DCAz      = GetLongitudinalDCA (track);
    Bool_t   hasTOFhit = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);
    Double_t length    = track->GetIntegratedLength();
    Double_t rap       = TMath::Abs(GetRapidity (track));

    
    //Track Quality Cuts
    if (!PassedTrackQualityCuts_Syst(track,isyst)) return;
        
    //DCA_{z} Cut
    if (TMath::Abs(DCAz)>1.0) return;
    
    
    //Vectors
    Double_t xTPC[4] = { rap, nsigmaTPC, pt, static_cast<Double_t>(isyst) };
    Double_t xTOF[4] = { rap, nsigmaTOF, pt, static_cast<Double_t>(isyst) };
    

    //DCA_{xy} Cut
    if (TMath::Abs(DCAxy)>0.1) return;


    //TPC-Only Analysis
    if (pt<1.5)  {
        
        if (charge>0) hnsigmaTPC_deuterons_rap_Syst     -> Fill (xTPC);
        if (charge<0) hnsigmaTPC_antideuterons_rap_Syst -> Fill (xTPC);
    }


    //TOF Analysis
    if (!hasTOFhit) return;
    if (length<350.0) return;
    if (TMath::Abs(nsigmaTPC)>3.0) return;
    if (pt<1.0) return;
        
            
    if (charge>0) hnsigmaTOF_deuterons_rap_Syst     -> Fill (xTOF);
    if (charge<0) hnsigmaTOF_antideuterons_rap_Syst -> Fill (xTOF);
    
}
//__________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskDeuteronsRT::FillHistograms_Rapidity_Systematics_Sim  (AliESDtrack *track, Int_t isyst)  {
    
    //Variables
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kDeuteron);
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kDeuteron);
    Double_t pt        = track->Pt();
    Double_t charge    = track->Charge();
    Double_t DCAxy     = GetTransverseDCA (track);
    Double_t DCAz      = GetLongitudinalDCA(track);
    Bool_t   hasTOFhit = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);
    Double_t length    = track->GetIntegratedLength();
    Double_t rap       = TMath::Abs(GetRapidity (track));

    //MC Particle
    AliMCParticle *particle = (AliMCParticle*) fMCEvent->GetTrack(TMath::Abs(track->GetLabel()));
       
    //Track Track Quality Cuts
    if (!PassedTrackQualityCuts_Syst(track,isyst)) return;
       
        
    //DCA cuts
    if (TMath::Abs(DCAz)>1.0)   return;
    if (TMath::Abs(DCAxy)>0.1) return;
    if (!particle->IsPhysicalPrimary())     return;


    //Vectors
    Double_t xTPC[4] = { rap, nsigmaTPC, pt, static_cast<Double_t>(isyst) };
    Double_t xTOF[4] = { rap, nsigmaTOF, pt, static_cast<Double_t>(isyst) };


    //TPC-Only Analysis
    if (pt<1.5)  {
      
        if (charge>0) hnsigmaTPC_deuterons_Rec_rap_Syst     -> Fill (xTPC);
        if (charge<0) hnsigmaTPC_antideuterons_Rec_rap_Syst -> Fill (xTPC);
    }
        
        
    //TOF Analysis
    if (!hasTOFhit) return;
    if (length<350.0) return;
    if (TMath::Abs(nsigmaTPC)>3.0) return;
    if (pt<1.0) return;
       
    if (charge>0) hnsigmaTOF_deuterons_Rec_rap_Syst     -> Fill (xTOF);
    if (charge<0) hnsigmaTOF_antideuterons_Rec_rap_Syst -> Fill (xTOF);

    
}
//__________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskDeuteronsRT::PassedTrackQualityCuts_Syst (AliESDtrack *track, Int_t isyst)  {

    //Initialization
    Bool_t passedTrkSelection=(kFALSE);


    //Track Variables
    Double_t nCrossedRows               = (Double_t)track->GetTPCCrossedRows();
    Double_t nCrossedRows_over_Findable = ((Double_t)track->GetTPCCrossedRows())/((Double_t)track->GetTPCNclsF());
    Double_t nClustersITS               = (Double_t)track->GetITSNcls();
    Double_t chi2TPC_ndf                = ((Double_t)track->GetTPCchi2())/((Double_t)track->GetTPCNcls());
    Double_t chi2ITS_ndf                = ((Double_t)track->GetITSchi2())/((Double_t)track->GetITSNcls());
    Double_t nClustersTPC_dEdx          = (Double_t)track->GetTPCsignalN();


    //Analysis Parameters
    Double_t nCrossedRows_Min               = hAnalysisParameters -> GetBinContent (1,(isyst+1));
    Double_t nCrossedRows_over_Findable_Min = hAnalysisParameters -> GetBinContent (2,(isyst+1));
    Double_t nClustersITS_Min               = hAnalysisParameters -> GetBinContent (3,(isyst+1));
    Double_t chi2TPC_ndf_Max                = hAnalysisParameters -> GetBinContent (4,(isyst+1));
    Double_t chi2ITS_ndf_Max                = hAnalysisParameters -> GetBinContent (5,(isyst+1));
    Double_t nClustersTPC_dEdx_Min          = hAnalysisParameters -> GetBinContent (6,(isyst+1));


    //Cuts
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
Double_t AliAnalysisTaskDeuteronsRT::GetTransverseDCA (AliESDtrack *track)  {
         
    Double_t impactParameter[2];
    Double_t covarianceMatrix[3];
    if (!track->PropagateToDCA (fESDevent->GetPrimaryVertex(),fESDevent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;
            
    Double_t DCAxy = impactParameter[0];
            
    return DCAxy;
}
//__________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskDeuteronsRT::GetLongitudinalDCA (AliESDtrack *track)  {
         
    Double_t impactParameter[2];
    Double_t covarianceMatrix[3];
    if (!track->PropagateToDCA (fESDevent->GetPrimaryVertex(),fESDevent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;
            
    Double_t DCAz = impactParameter[1];
            
    return DCAz;
}
//__________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskDeuteronsRT::Terminate(Option_t *)  {
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
//__________________________________________________________________________________________________________________________________________________
