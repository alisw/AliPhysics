/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* AliAnaysisTaskCaloHFEpp
 *
 * empty task which can serve as a starting point for building an analysis
 * as an example, one histogram is filled
 */

#include "TChain.h"
#include "TMath.h"
#include "TVector3.h"
#include "TH1F.h"
#include "TList.h"
#include "TRandom.h"
#include "THnSparse.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisTaskCaloHFEpp.h"
#include "AliCentrality.h"
#include "AliMultSelection.h"
#include "AliPIDResponse.h" 
#include "AliPID.h" 
#include "AliKFParticle.h"
#include "AliESDtrackCuts.h" 
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenPythiaEventHeader.h"
#include "AliVertexingHFUtils.h"


//using std::cout;
//using std::endl;

class AliAnalysisTaskCaloHFEpp;    // your analysis class

using namespace std;            // std namespace: so you can do things like 'cout'

ClassImp(AliAnalysisTaskCaloHFEpp) // classimp: necessary for root

AliAnalysisTaskCaloHFEpp::AliAnalysisTaskCaloHFEpp() : AliAnalysisTaskSE(), 
	fAOD(0), 
	fOutputList(0), 
	fVevent(0), 
	fMultSelection(0), 
	fpidResponse(0), 
	//==== Tender  ====
	fUseTender(kTRUE),
	fTracks_tender(0),
	fCaloClusters_tender(0),
	//==== cut parameters ====
	TrackEtaMin(0),
	TrackEtaMax(0),
	NTPCClust(0), 
	NITSClust(0), 
	NCrossedRow(0),
	DCAxy(0), 
	DCAz(0),
	NsigmaMin(0),
	NsigmaMax(0),
	M20Min(0), 
	M20Max(0),
	EopMin(0),
	EopMax(0),
	MaxConeR(0),
	ptAssoMin(0),
        CutMimClE(0.3),
	pTe("name"),
	massMin(0),
	Nref(0),
	Nch(0),
	MinNtr(0),
	MaxNtr(0),
	//==== basic parameters ====
	fNevents(0),
	fNDB(0),
	fHist_VertexZ(0),
	fHist_VertexZ_all(0),
	fHist_Centrality(0),
	fHist_Mult(0),
	fTrigMulti(0),
	fHistEta_track(0),
	fHistPhi_track(0),
	fHistEta_EMcal(0),
	fHistPhi_EMcal(0),
	fHistScatter_EMcal(0),
	fHistScatter_EMcal_aftMatch(0),
	fHistoNCells(0),
	fM02(0),
	fM20(0),
	fM02_ele(0),
	fM20_ele(0),
	fM02_had(0),
	fM20_had(0),
	//==== check cut parameters ====
	fTPCNcls(0),
	fITSNcls(0),
	fTPCCrossedRow(0),
	fTPCnsig_ele(0),
	fTPCnsig_iso(0),
	fM02_2(0),
	fM20_2(0),
	fEop_ele(0),
	fEop_iso(0),
	fEop_iso_eID(0),
	fConeR(0),
	fConeE(0),
	fNpart(0),
	//==== Real data output ====
	fHist_trackPt(0),
	fHistMatchPt(0),
	fHistSelectPt(0),
	fHist_ClustE(0),
	fHist_SelectClustE(0),
	fHistMatchE(0),
	fdEdx(0),
	fTPCnsig(0),
	fHistNsigEop(0),
	fEopPt_ele_loose(0),
	fEopPt_ele_tight(0),
	fEopPt_ele_tight_PYTHIA(0),
	fEopPt_ele_tight_forSys(0),
	fEopPt_had(0),
	fEtadiff(0),
	fPhidiff(0),
	fInv_pT_LS(0),
	fInv_pT_ULS(0),
	fInv_pT_ULS_forW(0),
	fHistPt_Inc(0),
	fHistPt_Iso(0),
	fHistPt_R_Iso(0),
	fRiso_phidiff(0),
	fRiso_phidiff_LS(0),
	fRiso_phidiff_35(0),
	fRiso_phidiff_LS_35(0),
        fIsoArray(0),
        fHFArray(0),
	fzvtx_Ntrkl(0),
	fzvtx_Nch(0),
	fzvtx_Ntrkl_Corr(0),
	fzvtx_Corr(0),
	fNtrkl_Corr(0),
	fNchNtr(0),
	fNchNtr_Corr(0),
	fDCAxy_Pt_ele(0),
	fDCAxy_Pt_had(0),
	fDCAxy_Pt_LS(0),
	fDCAxy_Pt_ULS(0),
	fDCAxy_Pt_Dpm(0),
	fDCAxy_Pt_D0(0),
	fDCAxy_Pt_Ds(0),
	fDCAxy_Pt_lambda(0),
	fDCAxy_Pt_B(0),
	fPt_Btoe(0),
	//==== Trigger or Calorimeter flag ====
	fEMCEG1(kFALSE),
	fEMCEG2(kFALSE),
	fDCDG1(kFALSE),
	fDCDG2(kFALSE),
	fFlagClsTypeEMC(kFALSE),
	fFlagClsTypeDCAL(kFALSE),
	//==== MC output ===
	fMCcheckMother(0),
	fMCarray(0),
	fMCparticle(0),
	fMCTrackpart(0),
	fMCheader(0),
	fCheckEtaMC(0),
	fHistMCorgPi0(0),
	fHistMCorgEta(0),
	fHistMCorgD(0),
	fHistMCorgB(0),
	NembMCpi0(0),
	NembMCeta(0),
	NpureMCproc(0),
	NpureMC(0),
	fHistPhoReco0(0),
	fHistPhoReco1(0),
	fHistPhoReco2(0),
	fHistPhoPi0(0), 
	fHistPhoPi1(0),
	fHistPhoEta0(0),
	fHistPhoEta1(0),
	fPi000(0),
	fPi005(0),
	fPi010(0),
	fEta000(0),
	fEta005(0),
	fEta010(0),
	fCorrZvtx(0),
	fCorrNtrkl(0),
	fHistPt_HFE_MC_D(0),
	fHistPt_HFE_MC_B(0),
	fHistPt_HFE_PYTHIA(0),
	fHistPt_HFE_emb(0),
	fHistPt_HFE_Gen(0),
	fHistPt_HFE_GenvsReco(0),
	fHist_eff_HFE(0),
	fHist_eff_match(0),
	fHist_eff_TPC(0),
	fHist_eff_M20(0),
	fHist_eff_Iso(0),
        fHistWeOrg(0),
	fweightNtrkl(0)

{
	// default constructor, don't allocate memory here!
	// this is used by root for IO purposes, it needs to remain empty
	for(Int_t i=0; i<13; i++) fMultEstimatorAvg[i]=0;
}
//_____________________________________________________________________________
AliAnalysisTaskCaloHFEpp::AliAnalysisTaskCaloHFEpp(const char* name) : AliAnalysisTaskSE(name),
	fAOD(0), 
	fOutputList(0), 
	fVevent(0), 
	fMultSelection(0), 
	fpidResponse(0), 
	//==== Tender  ====
	fUseTender(kTRUE),
	fTracks_tender(0),
	fCaloClusters_tender(0),
	//==== cut parameters ====
	TrackEtaMin(0),
	TrackEtaMax(0),
	NTPCClust(0), 
	NITSClust(0), 
	NCrossedRow(0),
	DCAxy(0), 
	DCAz(0),
	NsigmaMin(0),
	NsigmaMax(0),
	M20Min(0), 
	M20Max(0),
	EopMin(0),
	EopMax(0),
	MaxConeR(0),
	ptAssoMin(0),
        CutMimClE(0.3),
	pTe("name"),
	massMin(0),
	Nref(0),
	Nch(0),
	MinNtr(0),
	MaxNtr(0),
	//==== basic parameters ====
	fNevents(0),
	fNDB(0),
	fHist_VertexZ(0),
	fHist_VertexZ_all(0),
	fHist_Centrality(0),
	fHist_Mult(0),
	fTrigMulti(0),
	fHistEta_track(0),
	fHistPhi_track(0),
	fHistEta_EMcal(0),
	fHistPhi_EMcal(0),
	fHistScatter_EMcal(0),
	fHistScatter_EMcal_aftMatch(0),
	fHistoNCells(0),
	fM02(0),
	fM20(0),
	fM02_ele(0),
	fM20_ele(0),
	fM02_had(0),
	fM20_had(0),
	//==== check cut parameters ====
	fTPCNcls(0),
	fITSNcls(0),
	fTPCCrossedRow(0),
	fTPCnsig_ele(0),
	fTPCnsig_iso(0),
	fM02_2(0),
	fM20_2(0),
	fEop_ele(0),
	fEop_iso(0),
	fEop_iso_eID(0),
	fConeR(0),
	fConeE(0),
	fNpart(0),
	//==== Real data output ====
	fHist_trackPt(0),
	fHistMatchPt(0),
	fHistSelectPt(0),
	fHist_ClustE(0),
	fHist_SelectClustE(0),
	fHistMatchE(0),
	fdEdx(0),
	fTPCnsig(0),
	fHistNsigEop(0),
	fEopPt_ele_loose(0),
	fEopPt_ele_tight(0),
	fEopPt_ele_tight_PYTHIA(0),
	fEopPt_ele_tight_forSys(0),
	fEopPt_had(0),
	fEtadiff(0),
	fPhidiff(0),
	fInv_pT_LS(0),
	fInv_pT_ULS(0),
	fInv_pT_ULS_forW(0),
	fHistPt_Inc(0),
	fHistPt_Iso(0),
	fHistPt_R_Iso(0),
	fRiso_phidiff(0),
	fRiso_phidiff_LS(0),
	fRiso_phidiff_35(0),
	fRiso_phidiff_LS_35(0),
        fIsoArray(0),
        fHFArray(0),
	fzvtx_Ntrkl(0),
	fzvtx_Nch(0),
	fzvtx_Ntrkl_Corr(0),
	fzvtx_Corr(0),
	fNtrkl_Corr(0),
	fNchNtr(0),
	fNchNtr_Corr(0),
	fDCAxy_Pt_ele(0),
	fDCAxy_Pt_had(0),
	fDCAxy_Pt_LS(0),
	fDCAxy_Pt_ULS(0),
	fDCAxy_Pt_Dpm(0),
	fDCAxy_Pt_D0(0),
	fDCAxy_Pt_Ds(0),
	fDCAxy_Pt_lambda(0),
	fDCAxy_Pt_B(0),
	fPt_Btoe(0),
	//==== Trigger or Calorimeter flag ====
	fEMCEG1(kFALSE),
	fEMCEG2(kFALSE),
	fDCDG1(kFALSE),
	fDCDG2(kFALSE),
	fFlagClsTypeEMC(kFALSE),
	fFlagClsTypeDCAL(kFALSE),
	//==== MC output ===
	fMCcheckMother(0),
	fMCarray(0),
	fMCparticle(0),
	fMCTrackpart(0),
	fMCheader(0),
	fCheckEtaMC(0),
	fHistMCorgPi0(0),
	fHistMCorgEta(0),
	fHistMCorgD(0),
	fHistMCorgB(0),
	NembMCpi0(0),
	NembMCeta(0),
	NpureMCproc(0),
	NpureMC(0),
	fHistPhoReco0(0),
	fHistPhoReco1(0),
	fHistPhoReco2(0),
	fHistPhoPi0(0), 
	fHistPhoPi1(0),
	fHistPhoEta0(0),
	fHistPhoEta1(0),
	fPi000(0),
	fPi005(0),
	fPi010(0),
	fEta000(0),
	fEta005(0),
	fEta010(0),
	fCorrZvtx(0),
	fCorrNtrkl(0),
	fHistPt_HFE_MC_D(0),
	fHistPt_HFE_MC_B(0),
	fHistPt_HFE_PYTHIA(0),
	fHistPt_HFE_emb(0),
	fHistPt_HFE_Gen(0),
	fHistPt_HFE_GenvsReco(0),
	fHist_eff_HFE(0),
	fHist_eff_match(0),
	fHist_eff_TPC(0),
	fHist_eff_M20(0),
	fHist_eff_Iso(0),
        fHistWeOrg(0),
	fweightNtrkl(0)
{
	// constructor
	DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
	// this chain is created by the analysis manager, so no need to worry about it, 
	// it does its work automatically
	DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
	// you can add more output objects by calling DefineOutput(2, classname::Class())
	// if you add more output objects, make sure to call PostData for all of them, and to
	// make changes to your AddTask macro!
	for(Int_t i=0; i<13; i++) fMultEstimatorAvg[i]=0;
}
//_____________________________________________________________________________
AliAnalysisTaskCaloHFEpp::~AliAnalysisTaskCaloHFEpp()
{
	// destructor
	if(fOutputList) {
		delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
	}
	for(Int_t i=0; i<13; i++) {
		if (fMultEstimatorAvg[i]) delete fMultEstimatorAvg[i];
	}
	if(fweightNtrkl) delete fweightNtrkl;
}
//_____________________________________________________________________________
void AliAnalysisTaskCaloHFEpp::UserCreateOutputObjects()
{
	// create output objects
	//
	// this function is called ONCE at the start of your analysis (RUNTIME)
	// here you ceate the histograms that you want to use 
	//
	// the histograms are in this case added to a tlist, this list is in the end saved
	// to an output file
	//
	fOutputList = new TList();          // this is a list which will contain all of your histograms
	// at the end of the analysis, the contents of this list are written
	// to the output file
	fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
	// if requested (dont worry about this now)

	// example of a histogram
	fHist_trackPt = new TH1F("fHist_trackPt", "EMCAL cluster pt distributiont; pt(GeV/c); counts", 1000, 0, 100);       // create your histogra
	fHistMatchPt = new TH1F("fHistMatchPt", "EMCAL matched cluster pt distributiont; pt(GeV/c); counts", 1000, 0, 100);       // create your histogra
	fHistSelectPt = new TH1F("fHistSelectPt", "EMCAL Slected cluster pt distributiont; pt(GeV/c); counts", 1000, 0, 100);       // create your histogra
	fHist_ClustE = new TH1F("fHist_ClustE", "fHist_ClustE; EMCAL cluster energy distribution; counts", 1000, 0, 100);      
	fHist_SelectClustE = new TH1F("fHist_SelectClustE", "fHistE; EMCAL cluster energy distribution before selection; counts", 1000, 0, 100);      
	fHistMatchE = new TH1F("fHistMatchE", "fHistMatchE; EMCAL Matched cluster energy distribution; counts", 1000, 0, 100);      
	fHistEta_track = new TH1F("fHistEta_track", "Track #eta distribution; pt(GeV/c); counts", 200, -4, 4);    
	fHistPhi_track = new TH1F("fHistPhi_track", "Track #phi distribution; #phi; counts", 200, 0, 10);    
	fHistEta_EMcal = new TH1F("fHistEta_EMcal", "EMCAL selected cluster #eta distribution; #eta; counts", 200, -4, 4);    
	fHistPhi_EMcal = new TH1F("fHistPhi_EMcal", "EMCAL selected cluster #phi distribution; #phi; counts", 200, 0, 10);    
	fHist_VertexZ = new TH1F("fHist_VertexZ", "Z Vertex position; Vtx_{z}; counts", 200, -25, 25);     
	fHist_VertexZ_all = new TH1F("fHist_VertexZ_all", "All z vertex position; Vtx_{z}; counts", 600, -30, 30);     
	fHist_Centrality = new TH1F("fHist_Centrality", "Centrality", 100, 0, 100);
	fNevents = new TH1F("fNevents","No of events",8,-0.5,7.5);
	fNDB = new TH1F("fNDB","No of events",2,-0.5,1.5);
	fTPCNcls = new TH1F("fTPCNcls","No of TPC clusters; N^{TPC}_{cls}; counts",100,0.0,200.);           
	fITSNcls = new TH1F("fITSNcls","No of ITS clusters; N^{ITS}_{cls}; counts",100,0.0,20.); 
	fTPCCrossedRow = new TH1F("fTPCCrossedRow","No of TPC CrossedRow; N^{ITS}_{CrossedRow}; counts",500,0.,500.); 
	fEtadiff = new TH1F("fEtadiff", "Distance of EMCAL to its closest track(Eta)", 60,-0.3,0.3);
	fPhidiff = new TH1F("fPhidiff", "Distance of EMCAL to its closest track(Phi)", 60,-0.3,0.3);
	fMCcheckMother = new TH1F("fMCcheckMother", "Mother MC PDG", 1000,-0.5,999.5);
	fHistPhoReco0 = new TH1D("fHistPhoReco0", "total pho in sample; p_{T}(GeV/c)", 600,0,60);
	fHistPhoReco1 = new TH1D("fHistPhoReco1", "reco pho in sample; p_{T}(GeV/c)", 600,0,60);
	fHistPhoReco2 = new TH1D("fHistPhoReco2", "non-reco pho in sample; p_{T}(GeV/c)", 600,0,60);
	fHistPhoPi0 = new TH1D("fHistPhoPi0", "total pi0 in sample; p_{T}(GeV/c)", 600,0,60);
	fHistPhoPi0->Sumw2(); 
	fHistPhoPi1 = new TH1D("fHistPhoPi1", "reco pi0 in sample; p_{T}(GeV/c)", 600,0,60);
	fHistPhoPi1->Sumw2(); 
	fHistPhoEta0 = new TH1D("fHistPhoEta0", "total Eta in sample; p_{T}(GeV/c)", 600,0,60);
	fHistPhoEta0->Sumw2(); 
	fHistPhoEta1 = new TH1D("fHistPhoEta1", "reco Eta in sample; p_{T}(GeV/c)", 600,0,60);
	fHistPhoEta1->Sumw2(); 

	fCheckEtaMC = new TH1F("fCheckEtaMC","check Eta range cut in MC",160,-0.8,0.8);
	fHistMCorgD = new TH1F("fHistMCorgD","MC org D",600,0,60);
	fHistMCorgB = new TH1F("fHistMCorgB","MC org B",600,0,60);
	fHistPt_HFE_MC_D  = new TH1F("fHistPt_HFE_MC_D","HFE from D MC",600,0,60);
	fHistPt_HFE_MC_B  = new TH1F("fHistPt_HFE_MC_B","HFE fron B MC",600,0,60);
	fHistPt_HFE_PYTHIA = new TH1F("fHistPt_HFE_PYTHIA","HFE fron PYTHIA",600,0,60);
	fHistPt_HFE_emb    = new TH1F("fHistPt_HFE_emb","HFE fron embedding",600,0,60);
	fHistPt_HFE_Gen    = new TH1F("fHistPt_HFE_Gen","HFE pt spectrum before reconstruct",600,0,60);
	fHistPt_HFE_GenvsReco    = new TH2F("fHistPt_HFE_GenvsReco","HFE pt spectrum before reconstruct",600,0,60,600,0,60);
	fHistPt_Inc = new TH1F("fHistPt_Inc","Inclusive electron",600,0,60);
	fHistPt_Iso = new TH1F("fHistPt_Iso","Isolated electron",600,0,60);
	fHistPt_R_Iso = new TH2F("fHistPt_R_Iso","Pt vs riso ",1000,0,100,500,0.,0.5);
	fRiso_phidiff       = new TH2F("fRiso_phidiff","phi differnce vs riso ",80,-3.,5.,500,0.,0.5);
	fRiso_phidiff_LS    = new TH2F("fRiso_phidiff_LS","phi differnce vs riso ",80,-3.,5.,500,0.,0.5);
	fRiso_phidiff_35    = new TH2F("fRiso_phidiff_35","phi differnce vs riso ",80,-3.,5.,500,0.,0.5);
	fRiso_phidiff_LS_35 = new TH2F("fRiso_phidiff_LS_35","phi differnce vs riso ",80,-3.,5.,500,0.,0.5);
	
        Int_t bins[10]=   { 90, 100, 200, 500, 100, 100, 100, 20, 500, 20}; //pt, TPCnsig, E/p, M20, NTPC,nITS, particle pt
        Double_t xmin[10]={ 10,  -5,   0,   0,   0,   0,   0,  0,   0,  0};
        Double_t xmax[10]={100,   5,   2, 0.5, 100,   1,   1, 20, 0.5, 20};
        fIsoArray = new THnSparseD ("fIsoArray","Isolation ;pT;nSigma;eop;iso;truePt;m20;m02;Ncont;isotrack;NtrCont",10,bins,xmin,xmax);
        fOutputList->Add(fIsoArray);

        fHFArray = new THnSparseD ("fHFArray","Isolation ;pT;nSigma;eop;iso;truePt;m20;m02;Ncont;isotrack;NtrCont",10,bins,xmin,xmax);
        fOutputList->Add(fHFArray);

        fzvtx_Ntrkl = new TH2F("fzvtx_Ntrkl","Zvertex vs N tracklet; zvtx; SPD Tracklets",400,-20.,20.,301,-0.5,300.5);
	fzvtx_Nch = new TH2F("fzvtx_Nch","Zvertex vs N charged; zvtx; N_{ch}",400,-20.,20.,301,-0.5,300.5);
	fzvtx_Ntrkl_Corr = new TH2F("fzvtx_Ntrkl_Corr","Zvertex vs N tracklet after correction; zvtx; SPD Tracklets",400,-20.,20.,301,-0.5,300.5);
	fzvtx_Corr = new TH1F("fzvtx_Corr","Zvertex after correction; zvtx; counts",400,-20.,20.);
	fNtrkl_Corr = new TH1F("fNtrkl_Corr","N_{tracklet} after correction; zvtx; counts",301,-0.5,300.5);
	fNchNtr = new TH2F("fNchNtr","N tracklet after correction vs N charged; n^{corr}_{trkl}; N_{ch}",301,-0.5,300.5,301,-0.5,300.5);
	fNchNtr_Corr = new TH2F("fNchNtr_Corr","N tracklet after correction vs N charged; n^{corr}_{trkl}; N_{ch}",301,-0.5,300.5,301,-0.5,300.5);

	fDCAxy_Pt_ele = new TH2F("fDCAxy_Pt_ele","DCA_{xy} vs Pt (electron);p_{t} (GeV/c);DCAxy*charge*Bsign",600,0,60,800,-0.2,0.2);
	fDCAxy_Pt_had = new TH2F("fDCAxy_Pt_had","DCA_{xy} vs Pt (hadron);p_{t} (GeV/c);DCAxy*charge*Bsign",600,0,60,800,-0.2,0.2);
	fDCAxy_Pt_LS = new TH2F("fDCAxy_Pt_LS","DCA_{xy} vs Pt LS pairs;p_{t} (GeV/c);DCAxy*charge*Bsign",600,0,60,800,-0.2,0.2);
	fDCAxy_Pt_ULS = new TH2F("fDCAxy_Pt_ULS","DCA_{xy} vs Pt ULS pairs;p_{t} (GeV/c);DCAxy*charge*Bsign",600,0,60,800,-0.2,0.2);
	fDCAxy_Pt_Dpm = new TH2F("fDCAxy_Pt_Dpm","DCA_{xy} vs Pt D+-(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",600,0,60,800,-0.2,0.2);
	fDCAxy_Pt_D0= new TH2F("fDCAxy_Pt_D0","DCA_{xy} vs Pt D0(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",600,0,60,800,-0.2,0.2);
	fDCAxy_Pt_Ds= new TH2F("fDCAxy_Pt_Ds","DCA_{xy} vs Pt Ds(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",600,0,60,800,-0.2,0.2);
	fDCAxy_Pt_lambda = new TH2F("fDCAxy_Pt_lambda","DCA_{xy} vs Pt lambda(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",600,0,60,800,-0.2,0.2);
	fDCAxy_Pt_B= new TH2F("fDCAxy_Pt_B","DCA_{xy} vs Pt all B meson(MC);p_{t} (GeV/c);DCAxy*charge*Bsign",600,0,60,800,-0.2,0.2);

	Double_t eop_range[21] = {2.,2.5,3.,3.5,4.,4.5,5.0,5.5,6.,8.,10.,12.,14.,16.,19.,22.,26.,30.,35.,40.,50.};
	fPt_Btoe = new TH2F("fPt_Btoe","B meson vs electron;electron p_{t} (GeV/c);B p_{t} (GeV/c)",20,eop_range,600,0.,60.);
	fHist_eff_HFE     = new TH1F("fHist_eff_HFE","efficiency :: HFE",600,0,60);
	fHist_eff_match   = new TH1F("fHist_eff_match","efficiency :: matched cluster",600,0,60);
	fHist_eff_TPC     = new TH1F("fHist_eff_TPC","efficiency :: TPC cut",600,0,60);
	fHist_eff_M20     = new TH1F("fHist_eff_M20","efficiency :: shower shape cut",600,0,60);
	fHist_eff_Iso     = new TH2F("fHist_eff_Iso","efficiency :: shower shape cut",600,0,60,500,0.,0.5);
	fHistWeOrg        = new TH1F("fHistWeOrg","particle level W->e",90,10,100);



	/////////////////
	//pi0 weight			
	/////////////////
	fPi000 = new TF1("fPi000","[0]*x/pow([1]+x/[2]+x*x/[3],[4])");
	fPi000->SetParameters(3.68528e+01,5.43694e-02,1.99270e+00,5.33945e+00,3.08814e+00);
	fPi005 = new TF1("fPi005","[0]*x/pow([1]+x/[2]+x*x/[3],[4])");
	fPi005->SetParameters(4.03327e+01,-2.14631e+01,9.93386e+02,3.60954e+00,2.34609e+00);
	fPi010 = new TF1("fPi010","[0]*x/pow([1]+x/[2]+x*x/[3],[4])");
	fPi010->SetParameters(5.31753e+00,-7.43035e+01,3.27764e+04,4.94405e+00,1.93347e+00);

	/////////////////
	// Eta weight
	/////////////////
	fEta000 = new TF1("fEta000","[0]*x/pow([1]+x/[2]+x*x/[3],[4])");
	fEta000->SetParameters(1.50102e+01,2.08498e-01,2.95617e+00,5.05032e+00,2.95377e+00);    
	fEta005 = new TF1("fEta005","[0]*x/pow([1]+x/[2]+x*x/[3],[4])");
	fEta005->SetParameters(3.00390e+01,-1.76773e+01,7.45941e+00,4.48491e+00,2.41261e+00);
	fEta010 = new TF1("fEta010","[0]*x/pow([1]+x/[2]+x*x/[3],[4])");
	fEta010->SetParameters(1.82736e+00,-8.08208e+01,2.32670e+04,4.66500e+00,1.75496e+00); 

	fCorrZvtx = new TF1("fCorrZvtx","pol4");
	fCorrZvtx->SetParameters(0.989494,-0.000148672,0.00145737,4.02038e-05,-2.07891e-05);

	fCorrNtrkl = new TF1("fCorrNtrkl","pol8");
	fCorrNtrkl->SetParameters(1.32267,-0.186463,0.0300275,-0.00205886,7.30481e-05,-1.45858e-06,1.65515e-08,-9.94577e-11,2.45483e-13);


	fHist_Mult = new TH2F("fMult","Track multiplicity",100,0,100,20000,0,20000);
	fHistScatter_EMcal = new TH2F("fHistScatter_EMcal", "EMCAL cluster scatter plot; #eta; #phi", 200,0.,6.,200, -1., 1.);       // create your histogra
	fHistScatter_EMcal_aftMatch = new TH2F("fHistScatter_EMcal_aftMatch", "EMCAL cluster scatter plot after track matching; #eta; #phi", 40,-1.0,1.0,200, 0., 6.);       // create your histogra
	fdEdx = new TH2F("fdEdx","All Track dE/dx distribution;p (GeV/c);dE/dx",500,0,50,500,0,160);
	fTPCnsig = new TH2F("fTPCnsig","All Track TPC Nsigma distribution;p (GeV/c);#sigma_{TPC-dE/dx}",1000,0,50,200,-10,10);
	fTPCnsig_ele = new TH2F("fTPCnsig_ele","electron TPC Nsigma distribution;#sigma_{TPC-dE/dx} ; counts",100,0,100,200,-10,10);
	fTPCnsig_iso = new TH2F("fTPCnsig_iso","isolated electron TPC Nsigma distribution;#sigma_{TPC-dE/dx} ; counts",90,10,100,200,-10,10);
	fHistNsigEop = new TH2F ("fHistNsigEop", "E/p vs TPC nsig; E/p; #sigme_{TPC-dE/dX}",300, 0.0, 3.0, 200, -10,10);   
	fM02 = new TH2F ("fM02","M02 vs pt distribution; pt(GeV/c); M02",500,0,50,400,0,2);
	fM20 = new TH2F ("fM20","M20 vs pt distribution; pt(GeV/c); M20",500,0,50,400,0,2);
	fM02_ele = new TH1F ("fM02_ele","M02 distribution ele; pt(GeV/c); M02",400,0,2);
	fM20_ele = new TH1F ("fM20_ele","M20 distribution ele; pt(GeV/c); M20",400,0,2);
	fM02_had = new TH1F ("fM02_had","M02 distribution had; pt(GeV/c); M02",400,0,2);
	fM20_had = new TH1F ("fM20_had","M20 distribution had; pt(GeV/c); M20",400,0,2);
	fM02_2 = new TH2F ("fM02_2","M02 vs pt distribution (-1<nSigma<3 & 0.9<E/p<1.3); pt(GeV/c); M02",500,0,50,400,0,2);
	fM20_2 = new TH2F ("fM20_2","M20 vs pt distribution (-1<nSigma<3 & 0.9<E/p<1.3); pt(GeV/c); M20",500,0,50,400,0,2);
	fEopPt_ele_loose = new TH2F ("fEopPt_ele_loose","pt vs E/p distribution (-3<nSigma<3); pt(GeV/c); E/p",500,0,50,300,0,3.0);
	fEopPt_ele_tight = new TH2F ("fEopPt_ele_tight","pt vs E/p distribution (-1<nSigma<3); pt(GeV/c); E/p",500,0,50,300,0,3.0);
	fEopPt_ele_tight_PYTHIA = new TH2F ("fEopPt_ele_tight_PYTHIA","pt vs E/p distribution (-1<nSigma<3); pt(GeV/c); E/p",500,0,50,300,0,3.0);
	fEopPt_ele_tight_forSys = new TH2F ("fEopPt_ele_tight_forSys","pt vs E/p distribution (-1<nSigma<3); pt(GeV/c); E/p",600,0,60,300,0,3.0);
	fEopPt_had = new TH2F ("fEopPt_had","pt vs E/p distribution (nSigma<-3.5); pt(GeV/c); E/p",500,0,50,300,0,3.0);
	fEop_ele = new TH1F ("fEop_ele"," electron E/p distribution ; E/p ; counts",300,0,3.0);
	fEop_iso = new TH2F ("fEop_iso"," isolated electron E/p distribution ; E/p ; counts",90,10,100,300,0,3.0);
	fEop_iso_eID = new TH2F ("fEop_iso_eID"," isolated electron E/p distribution ; E/p ; counts",90,10,100,300,0,3.0);
	fConeR = new TH2F ("fConeR"," check cone radius; p_{T}; counts",100,0,100,500,0,0.5);
	fConeE = new TH2F ("fConeE"," check cone Energy; p_{T}; energy",100,0,100,500,0,50);
	fNpart = new TH2F ("fNpart"," check # of particles in cone radius; p_{T}; counts",100,0,100,50,0,50);
	fHistoNCells = new TH2F("fHistoNCells", "No of EMCAL cells in a cluster; Cluster E; N^{EMC}_{cells}",500,0,50,30,0,30);
	fInv_pT_ULS = new TH2F("fInv_pT_ULS", "Invariant mass vs p_{T} distribution(ULS) ; pt(GeV/c) ; mass(GeV/c^2)",500,0,50,1000,0,1.0);
	fInv_pT_ULS_forW = new TH2F("fInv_pT_ULS_forW", "Invariant mass vs p_{T} distribution(ULS) ; pt(GeV/c) ; mass(GeV/c^2)",500,0,50,1000,0,100);
	fInv_pT_LS = new TH2F("fInv_pT_LS", "Invariant mass vs p_{T} distribution(LS) ; pt(GeV/c) ; mass(GeV/c^2)",500,0,50,1000,0,1.0);
	fHistMCorgPi0 = new TH2F("fHistMCorgPi0","MC org Pi0",2,-0.5,1.5,100,0,50);
	fHistMCorgEta = new TH2F("fHistMCorgEta","MC org Eta",2,-0.5,1.5,100,0,50);
	fTrigMulti = new TH2F("fTrigMulti","Multiplicity distribution for different triggers; Trigger type; multiplicity",11,-1,10,2000,0,2000);


	//==== basic parameters ====
	fOutputList->Add(fNevents);
	fOutputList->Add(fNDB);
	fOutputList->Add(fHist_VertexZ);          
	fOutputList->Add(fHist_VertexZ_all);          
	fOutputList->Add(fHist_Centrality);       
	fOutputList->Add(fHist_Mult);           
	fOutputList->Add(fTrigMulti);
	fOutputList->Add(fHistEta_track);         
	fOutputList->Add(fHistPhi_track);         
	fOutputList->Add(fHistEta_EMcal);         
	fOutputList->Add(fHistPhi_EMcal);         
	fOutputList->Add(fHistScatter_EMcal);     
	fOutputList->Add(fHistScatter_EMcal_aftMatch);     
	fOutputList->Add(fHistoNCells);
	fOutputList->Add(fM02);
	fOutputList->Add(fM20);
	fOutputList->Add(fM02_ele);
	fOutputList->Add(fM20_ele);
	fOutputList->Add(fM02_had);
	fOutputList->Add(fM20_had);
	//==== check cut parameters ====
	fOutputList->Add(fTPCNcls);
	fOutputList->Add(fITSNcls);
	fOutputList->Add(fTPCCrossedRow);
	fOutputList->Add(fTPCnsig_ele);
	fOutputList->Add(fTPCnsig_iso);
	fOutputList->Add(fM02_2);
	fOutputList->Add(fM20_2);
	fOutputList->Add(fEop_ele);
	fOutputList->Add(fEop_iso);
	fOutputList->Add(fConeR);
	fOutputList->Add(fConeE);
	fOutputList->Add(fNpart);
	//==== Real data output ====
	fOutputList->Add(fHist_trackPt);          
	fOutputList->Add(fHistMatchPt);          
	fOutputList->Add(fHistSelectPt);          
	fOutputList->Add(fHist_ClustE);          
	fOutputList->Add(fHist_SelectClustE);          
	fOutputList->Add(fHistMatchE);          
	fOutputList->Add(fdEdx);
	fOutputList->Add(fTPCnsig);
	fOutputList->Add(fHistNsigEop);
	fOutputList->Add(fEopPt_ele_loose);
	fOutputList->Add(fEopPt_ele_tight);
	fOutputList->Add(fEopPt_ele_tight_PYTHIA);
	fOutputList->Add(fEopPt_ele_tight_forSys);
	fOutputList->Add(fEopPt_had);
	fOutputList->Add(fEtadiff);
	fOutputList->Add(fPhidiff);
	fOutputList->Add(fInv_pT_LS);
	fOutputList->Add(fInv_pT_ULS);
	fOutputList->Add(fInv_pT_ULS_forW);
	fOutputList->Add(fHistPt_Inc);
	fOutputList->Add(fHistPt_Iso);
	fOutputList->Add(fHistPt_R_Iso);
	fOutputList->Add(fRiso_phidiff);
	fOutputList->Add(fRiso_phidiff_LS);
	fOutputList->Add(fRiso_phidiff_35);
	fOutputList->Add(fRiso_phidiff_LS_35);
	fOutputList->Add(fzvtx_Ntrkl);
	fOutputList->Add(fzvtx_Nch);
	fOutputList->Add(fzvtx_Ntrkl_Corr);
	fOutputList->Add(fzvtx_Corr);
	fOutputList->Add(fNtrkl_Corr);
	fOutputList->Add(fNchNtr);
	fOutputList->Add(fNchNtr_Corr);
	fOutputList->Add(fDCAxy_Pt_ele);
	fOutputList->Add(fDCAxy_Pt_had);
	fOutputList->Add(fDCAxy_Pt_LS);
	fOutputList->Add(fDCAxy_Pt_ULS);
	fOutputList->Add(fDCAxy_Pt_Dpm);
	fOutputList->Add(fDCAxy_Pt_D0);
	fOutputList->Add(fDCAxy_Pt_Ds);
	fOutputList->Add(fDCAxy_Pt_lambda);
	fOutputList->Add(fDCAxy_Pt_B);
	fOutputList->Add(fPt_Btoe);
	//==== MC output ====
	fOutputList->Add(fMCcheckMother);
	fOutputList->Add(fCheckEtaMC);
	fOutputList->Add(fHistMCorgPi0);
	fOutputList->Add(fHistMCorgEta);
	fOutputList->Add(fHistMCorgD);
	fOutputList->Add(fHistMCorgB);
	fOutputList->Add(fHistPhoReco0);
	fOutputList->Add(fHistPhoReco1);
	fOutputList->Add(fHistPhoReco2);
	fOutputList->Add(fHistPhoPi0);
	fOutputList->Add(fHistPhoPi1);
	fOutputList->Add(fHistPhoEta0);
	fOutputList->Add(fHistPhoEta1);
	fOutputList->Add(fHistPt_HFE_MC_D);
	fOutputList->Add(fHistPt_HFE_MC_B);
	fOutputList->Add(fHistPt_HFE_PYTHIA);
	fOutputList->Add(fHistPt_HFE_emb);
	fOutputList->Add(fHistPt_HFE_Gen);
	fOutputList->Add(fHistPt_HFE_GenvsReco);
	fOutputList->Add(fHist_eff_HFE); 
	fOutputList->Add(fHist_eff_match); 
	fOutputList->Add(fHist_eff_TPC); 
	fOutputList->Add(fHist_eff_M20); 
	fOutputList->Add(fHist_eff_Iso); 
	fOutputList->Add(fHistWeOrg); 


	PostData(1, fOutputList);           // postdata will notify the analysis manager of changes / updates to the 
	// fOutputList object. the manager will in the end take care of writing your output to file
	// so it needs to know what's in the output
}
//_____________________________________________________________________________
void AliAnalysisTaskCaloHFEpp::UserExec(Option_t *)
{


	//##################### Systematic Parameters ##################### //
	//---Track Cut
	Double_t CutTrackEta[2] = {TrackEtaMin,TrackEtaMax};
	Int_t CutTPCNCls = NTPCClust;
	Int_t CutITSNCls = NITSClust; 
	Int_t CutTPCNCrossedRow = NCrossedRow;
	Double_t CutDCAxy = DCAxy;
	Double_t CutDCAz  = DCAz;
	//---PID Cut
	Double_t CutNsigma[2] = {NsigmaMin,NsigmaMax};
	Double_t CutM20[2] = {M20Min,M20Max};
	Double_t CutEop[2] = {EopMin,EopMax};
	Double_t CutEopHad = -3.5;
	Double_t CutptAsso = ptAssoMin;
	TString  TriggerPt = pTe;
	Int_t CutMinNtr = MinNtr;
	Int_t CutMaxNtr = MaxNtr;
	//cout<< "!!!!!!!!!!! pTe ;; "<<TriggerPt.Data()<<endl;
	//################################################################# //

	UInt_t evSelMask=((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

	// user exec
	// this function is called once for each event
	// the manager will take care of reading the events from file, and with the static function InputEvent() you 
	// have access to the current event. 
	// once you return from the UserExec function, the manager will retrieve the next event from the chain
	fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    // get an event (called fAOD) from the input file
	// there's another event format (ESD) which works in a similar wya
	// but is more cpu/memory unfriendly. for now, we'll stick with aod's
	if(!fAOD) return;                                   // if the pointer to the event is empty (getting it failed) skip this event
	// example part: i'll show how to loop over the tracks in an event 
	// and extract some information from them which we'll store in a histogram

	fMCarray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
	fMCheader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));

	fVevent = dynamic_cast<AliVEvent*>(InputEvent());
	if (!fVevent) {
		printf("ERROR: fVEvent not available\n");
		return;
	}   


	//////////////////////////////
	//Get Tender
	//////////////////////////////
	Bool_t fFlagEMCalCorrection = kTRUE;
	if(fFlagEMCalCorrection){
		TString fTenderClusterName("caloClusters"); //default name
		TString fTenderTrackName("tracks"); //default name
		fTracks_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTenderTrackName));
		fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTenderClusterName));
	}

	//////////////////////////////
	//PID initialised   
	//////////////////////////////
	fpidResponse = fInputHandler->GetPIDResponse();


	/////////////////
	//trigger check//
	/////////////////
	TString firedTrigger;
	TString TriggerEG1("EG1");
	TString TriggerEG2("EG2");
	TString TriggerDG1("DG1");
	TString TriggerDG2("DG2");
	fVevent->GetFiredTriggerClasses();
	if(fAOD) firedTrigger = fAOD->GetFiredTriggerClasses();

	Bool_t EG1tr = kFALSE;
	Bool_t EG2tr = kFALSE;

	if(firedTrigger.Contains(TriggerEG1))EG1tr = kTRUE;
	if(firedTrigger.Contains(TriggerEG2))EG2tr = kTRUE;

	// Use EMCal & DCal
	if(fFlagClsTypeEMC && fFlagClsTypeDCAL)
	{
		if(fEMCEG1 && fDCDG1) if(!firedTrigger.Contains(TriggerEG1) && !firedTrigger.Contains(TriggerDG1)) return;
		if(fEMCEG2 && fDCDG2) if(!firedTrigger.Contains(TriggerEG2) && !firedTrigger.Contains(TriggerDG2)) return;
	}
	// Use only EMCal or DCal
	else
	{
		if(fEMCEG1){if(!firedTrigger.Contains(TriggerEG1))return;}
		if(fEMCEG2){if(!firedTrigger.Contains(TriggerEG2))return;}
		if(fDCDG1) {if(!firedTrigger.Contains(TriggerDG1))return;}
		if(fDCDG2) {if(!firedTrigger.Contains(TriggerDG2))return;}
	}

        //cout << "fFlagClsTypeEMC" << fFlagClsTypeEMC << " ; fEMCEG1 = " << fEMCEG1 << " ; firedTrigger = " << firedTrigger << endl;
        //cout << "fFlagClsTypeEMC" << fFlagClsTypeEMC << " ; fEMCEG2 = " << fEMCEG2 << " ; firedTrigger = " << firedTrigger << endl;

	Int_t trigger = -1;
	if (fAOD){
		AliAODHeader *header = dynamic_cast<AliAODHeader*>(fAOD->GetHeader());
		if(!header) AliFatal("Not a standard AOD");
		Double_t multiplicity = header->GetRefMultiplicity();

		fTrigMulti->Fill(-0.5, multiplicity);
		if(evSelMask & AliVEvent::kAny) fTrigMulti->Fill(0.5, multiplicity);
		if(evSelMask & AliVEvent::kMB) fTrigMulti->Fill(1.5, multiplicity);
		if(evSelMask & AliVEvent::kINT7) fTrigMulti->Fill(2.5, multiplicity);
		if(evSelMask & AliVEvent::kINT8) fTrigMulti->Fill(3.5, multiplicity);
		if(evSelMask & AliVEvent::kEMC1) fTrigMulti->Fill(4.5, multiplicity);
		if(evSelMask & AliVEvent::kEMC7) fTrigMulti->Fill(5.5, multiplicity);
		if(evSelMask & AliVEvent::kEMC8) fTrigMulti->Fill(6.5, multiplicity);
		if(evSelMask & AliVEvent::kEMCEJE) fTrigMulti->Fill(7.5, multiplicity);
		if(evSelMask & AliVEvent::kEMCEGA) fTrigMulti->Fill(8.5, multiplicity);
		if(evSelMask & AliVEvent::kEMCEGA & EG2tr) fTrigMulti->Fill(9.5, multiplicity);

		if(evSelMask & AliVEvent::kMB) trigger =0;
		if(evSelMask & AliVEvent::kINT7) trigger =1;
		if(evSelMask & AliVEvent::kINT8) trigger =2;
		if(evSelMask & AliVEvent::kEMC1) trigger =3;
		if(evSelMask & AliVEvent::kEMC7) trigger =4;
		if(evSelMask & AliVEvent::kEMC8) trigger =5;
		if(evSelMask & AliVEvent::kEMCEJE) trigger =6;
		if(evSelMask & AliVEvent::kEMCEGA) trigger =7;
	}





	//////////////////////////////
	//Centarality
	//////////////////////////////
	Double_t centrality = -1;
	AliCentrality *fCentrality = (AliCentrality*)fAOD->GetCentrality();
	//centrality = fCentrality->GetCentralityPercentile("V0M");
	//cout << "Centrality == " << fCentrality->GetCentralityPercentile("V0M") << endl;
	if(fAOD)fMultSelection = (AliMultSelection * ) fAOD->FindListObject("MultSelection"); 
	if(!fMultSelection) {
		//If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
		//AliWarning("AliMultSelection object not found!");
		centrality = fCentrality->GetCentralityPercentile("V0M");
	}else{
		//lPercentile = fMultSelection->GetMultiplicityPercentile("V0M");
		centrality = fMultSelection->GetMultiplicityPercentile("V0M"); 
	}


	//////////////////////////////
	// Event selection
	//////////////////////////////
	//==== Global Vtx ===
	fNevents->Fill(0);//all enent
	const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
	Double_t NcontV = pVtx->GetNContributors();
	Double_t Xvertex = pVtx->GetX();
	Double_t Yvertex = pVtx->GetY();
	Double_t Zvertex = pVtx->GetZ();
	fHist_VertexZ_all->Fill(Zvertex);
	//==== SPD Vtx ====
	const AliVVertex *pVtxSPD = fVevent->GetPrimaryVertexSPD();
	Double_t ZvertexSPD = pVtxSPD->GetZ();
	Double_t NcontVSPD = pVtxSPD->GetNContributors();
	Double_t cov[6]={0};
	pVtxSPD->GetCovarianceMatrix(cov);
	// 1. remove pile up events
	if(fVevent->IsPileupFromSPDInMultBins()) return;
	fNevents->Fill(1); 
	// 2. Global contributiors cut
	if(NcontV<2)return;
	fNevents->Fill(2); 
	// 3. SPD contributiors cut
	if(NcontVSPD<1)return;
	fNevents->Fill(3); 
	// 4. select events where SPD and primary vertex match//
	if(TMath::Abs(ZvertexSPD - Zvertex) > 0.5) return;
	fNevents->Fill(4); 
	// 5. SPD vertex resolution cut //
	if (TMath::Sqrt(cov[5]) > 0.25) return;
	fNevents->Fill(5); 
	// 6. Z Vtx position cut 
	if(TMath::Abs(Zvertex)>10.0)return;
	fNevents->Fill(6); 
	fHist_VertexZ->Fill(Zvertex);                     // plot the pt value of the track in a histogram


	//////////////////////////////
	// Get sign of B field
	//////////////////////////////
	int Bsign = 0;
	if(fAOD->GetMagneticField() < 0) Bsign = -1;
	if(fAOD->GetMagneticField() > 0) Bsign = 1;

	//////////////////////////////
	// Multiplicity
	//////////////////////////////
	Float_t lPercentiles[3];
	TString lNames[1] = {"SPDTracklets,V0M,ZNA"}; // You can specify here the estimator you want to use
	for(Int_t iEst=0; iEst<2; iEst++) lPercentiles[iEst] = 300;	
	AliMultSelection *MultSelection = 0x0;
	MultSelection = (AliMultSelection*)fAOD->FindListObject("MultSelection");
	if(MultSelection){
		for(Int_t iEst=0; iEst<2; iEst++)
			lPercentiles[iEst] = MultSelection->GetMultiplicityPercentile(lNames[iEst].Data()); // gives the multiplicity in percetile for the corresponding estimator.
	}else{

		AliInfo("Didn't find MultSelection!"); 
	}
	fHist_Centrality -> Fill(lPercentiles[0]);

	//------------SPDTracklets--------------------
	Int_t nTracklets = 0;
	Int_t nAcc = 0;
	Double_t etaRange = 1.0;

	AliAODTracklets *tracklets = static_cast<const AliAODEvent*>(fAOD)->GetTracklets();
	nTracklets = tracklets->GetNumberOfTracklets();
	for (Int_t nn = 0; nn < nTracklets; nn++) {
		Double_t theta = tracklets->GetTheta(nn);
		Double_t eta = -TMath::Log(TMath::Tan(theta/2.0));
		if (TMath::Abs(eta) < etaRange) nAcc++;
	}
	fzvtx_Ntrkl->Fill(Zvertex,nAcc);


	//-----------Tracklet correction-------------------------
	Double_t correctednAcc   = nAcc;
	Double_t fRefMult = Nref;
	Double_t WeightNtrkl = -1.;
	Double_t WeightZvtx = -1.;
	TProfile* estimatorAvg;
	if(!fMCarray)estimatorAvg = GetEstimatorHistogram(fAOD);
	if(fMCarray)estimatorAvg = GetEstimatorHistogramMC(fAOD);
	if(estimatorAvg){
		correctednAcc=static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,nAcc,Zvertex,fRefMult));
		//correctednAcc= AliAnalysisTaskCaloHFEpp::GetCorrectedNtrackletsD(estimatorAvg,nAcc,Zvertex,fRefMult);
	} 
	fzvtx_Ntrkl_Corr->Fill(Zvertex,correctednAcc);



	if(fMCarray)CheckMCgen(fMCheader,CutTrackEta[1]);
	if(fMCarray)GetMClevelWdecay(fMCheader);


	fNchNtr->Fill(correctednAcc,Nch);
	if(fMCarray){
		WeightZvtx = fCorrZvtx->Eval(Zvertex);
		WeightNtrkl = fweightNtrkl->GetBinContent(fweightNtrkl->FindBin(correctednAcc));
		fzvtx_Corr->Fill(Zvertex,WeightZvtx);
		fNtrkl_Corr->Fill(correctednAcc,WeightNtrkl);
		fNchNtr_Corr->Fill(correctednAcc,Nch,WeightNtrkl);
		fzvtx_Nch->Fill(Zvertex,Nch,WeightZvtx);
	}


	//////////////////////////////////
	// Separate by multiplicity class
	//////////////////////////////////
	if(correctednAcc<CutMinNtr || correctednAcc > CutMaxNtr) return;
	fNevents->Fill(7); 


	//////////////////////////////
	// EMCal cluster loop
	//////////////////////////////
	Int_t Nclust = -999;
	if(!fFlagEMCalCorrection)Nclust =  fVevent->GetNumberOfCaloClusters();
	if(fFlagEMCalCorrection)Nclust =  fCaloClusters_tender->GetEntries();

	Bool_t fClsTypeEMC = kFALSE, fClsTypeDCAL = kFALSE;;

	for(Int_t icl=0; icl<Nclust; icl++)
	{
		AliVCluster *clust = 0x0;     
		if(!fFlagEMCalCorrection)clust = (AliVCluster*)fVevent->GetCaloCluster(icl); // address cluster matched to track
		if(fFlagEMCalCorrection)clust = dynamic_cast<AliVCluster*>(fCaloClusters_tender->At(icl)); // address cluster matched to track

		fClsTypeEMC = kFALSE; fClsTypeDCAL = kFALSE;

		if(clust && clust->IsEMCAL())
		{
			Double_t clustE = clust->E();
			fHist_SelectClustE -> Fill(clustE);

			Float_t clustpos[3] = {0.};
			clust->GetPosition(clustpos);

			TVector3 pos(clustpos);

			Double_t Phi =  pos.Phi();
			if(Phi <0){Phi += 2*TMath::Pi();}

			if(Phi > 1.39 && Phi < 3.265) fClsTypeEMC = kTRUE; //EMCAL : 80 < phi < 187     
			if(Phi > 4.53 && Phi < 5.708) fClsTypeDCAL = kTRUE; //DCAL  : 260 < phi < 327

			//----selects EMCAL+DCAL clusters when fFlagClsTypeEMC and fFlagClsTypeDCAL is kTRUE
			if(fFlagClsTypeEMC && !fFlagClsTypeDCAL)
			{
				if(!fClsTypeEMC) continue; //selecting only EMCAL clusters
			}
			else if(fFlagClsTypeDCAL && !fFlagClsTypeEMC)
			{
				if(!fClsTypeDCAL) continue; //selecting only DCAL clusters
			}
			else{};


			fHistScatter_EMcal->Fill(Phi,pos.Eta());                     // plot the pt value of the track in a histogram
			fHistoNCells -> Fill(clustE, clust->GetNCells());
			fHist_ClustE->Fill(clustE);                     // plot the pt value of the track in a histogram
		}
	}

	//////////////////////////////
	// Track loop
	//////////////////////////////
	Int_t iTracks = -999;
	if(!fFlagEMCalCorrection)iTracks = fVevent->GetNumberOfTracks(); // see how many tracks there are in the event
	if(fFlagEMCalCorrection)iTracks = fTracks_tender->GetEntries();  // see how many tracks there are in the event
	//fHist_Mult->Fill(centrality,nTracks);

	for(Int_t i(0); i < iTracks; i++) {    // loop overall these tracks
		Double_t fTPCnSigma = -999, dEdx = -999, TrkP = -999, TrkPt = -999; 
		Double_t ITSchi2 = -999, TPCchi2NDF = -999;

		AliAODTrack* track;     // get a track (type AliAODTrack) from the event
		if(fFlagEMCalCorrection){
			AliVParticle* Vtrack = 0x0;
			Vtrack = dynamic_cast<AliVTrack*>(fTracks_tender->At(i));
			track = dynamic_cast<AliAODTrack*>(Vtrack);
		}
		if(!fFlagEMCalCorrection)track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));   // get a track (type AliAODTrack) from the event

		if(!track) continue;                            // if we failed, skip this track
		fHist_trackPt->Fill(track->Pt());               // plot the pt value of the track in a histogram
		dEdx = track->GetTPCsignal();
		TrkP = track->P();
		TrkPt = track->Pt();
		fTPCnSigma = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron); 
		ITSchi2 = track -> GetITSchi2();
		TPCchi2NDF = track -> Chi2perNDF();

		Int_t EMCalIndex = -1;
		EMCalIndex = track->GetEMCALcluster();  // get index of EMCal cluster which matched to track


		/////////////////////////
		// track cut
		/////////////////////////
		//==== 1.TPC and ITS refit cut ====
		if(!(track->GetStatus()&AliAODTrack::kITSrefit) || !(track->GetStatus()&AliAODTrack::kTPCrefit)) continue;
		//==== 2.AOD filter bit required ====
		if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue; //mimimum cuts
		//==== 3.TPC cluster cut ====
		if(track->GetTPCNcls() < CutTPCNCls) continue; 
		//==== 4.ITS cluster cut ====
		if(track->GetITSNcls() < CutITSNCls) continue;  
		//==== 5.SPD hit cut ====
		if(!(track -> HasPointOnITSLayer(0) || track -> HasPointOnITSLayer(1))) continue;
		//==== 6.Eta cut ====
		if(track->Eta()>CutTrackEta[1] || track->Eta()<CutTrackEta[0]) continue; 
		//==== 7.DCA cut ====
		Double_t DCA[2] = {-999.,-999.}, covar[3];
		if(track -> PropagateToDCA(pVtx,fVevent -> GetMagneticField(),20.,DCA,covar))
		{
			if(TMath::Abs(DCA[0]) > CutDCAxy || TMath::Abs(DCA[1]) > CutDCAz) continue;
		}
		//==== 8.chi2 cut ====
		if((ITSchi2 >= 25) || (TPCchi2NDF >= 4)) continue;
		//==== 9.NCrossedRow cut ====
		if(track -> GetTPCCrossedRows() < CutTPCNCrossedRow) continue;


		fdEdx->Fill(TrkP,dEdx);
		fTPCNcls->Fill(track->GetTPCNcls());
		fITSNcls->Fill(track->GetITSNcls());
		fTPCnsig->Fill(TrkP,fTPCnSigma);
		fHistEta_track->Fill(track->Eta());            
		fHistPhi_track->Fill(track->Phi());             
		fTPCCrossedRow->Fill(track -> GetTPCCrossedRows());


		///////////////////////
		// Get MC information//
		///////////////////////
		Int_t ilabel = TMath::Abs(track->GetLabel());
		Int_t pdg = -999;
		Double_t pid_ele = 0.0;
		Double_t pTmom = -1.0;
		Int_t pidM = -1;
		Int_t pdgorg = -1;
		Int_t pdgstatus = -1;
		Int_t ilabelM = -1;
		Double_t pTGMom = -1.0;
		Double_t pTpart = -1.0;
		Int_t pidGM = -1;
		Int_t ilabelGM = -1;
		Bool_t iEmbPi0 = kFALSE; 
		Bool_t iEmbEta = kFALSE;
		Bool_t pid_eleD = kFALSE;
		Bool_t pid_eleB = kFALSE;
		Bool_t pid_eleP = kFALSE;

		if(ilabel>0 && fMCarray)
		{
			fMCTrackpart = (AliAODMCParticle*) fMCarray->At(ilabel);
			pdg = fMCTrackpart->GetPdgCode();
                        pdgstatus = fMCTrackpart->GetStatus();
			pTpart = fMCTrackpart->Pt();
			if(TMath::Abs(pdg)==11)pid_ele = 1.0;
			if(pid_ele==1.0)
                          {
                           FindMother(fMCTrackpart, ilabelM, pidM, pTmom);
	                   FindWdecay(fMCTrackpart,ilabelM,pdgorg);
                       
                           //cout << "pidM = "<< pidM << endl; 
                           if(TMath::Abs(pdgorg)==24)cout << "W->e reco status = "<< fMCTrackpart->GetStatus() << endl; 
                          }

			pid_eleB = IsBdecay(pidM);
			pid_eleP = IsPdecay(pidM);
			pid_eleD = IsDdecay(pidM);
			if(pid_eleD || pid_eleB)fNDB->Fill(0);

			if(pid_eleD){
				AliAODMCParticle* fMCTrackpartMom = (AliAODMCParticle*) fMCarray->At(ilabelM);
				FindMother(fMCTrackpartMom,ilabelGM,pidGM,pTGMom);
				if(IsBdecay(pidGM)){
					pid_eleB = IsBdecay(pidGM);
					//cout<<"B->D->e"<<endl;
					pid_eleD = kFALSE;
				}
			}

			if(pid_eleD || pid_eleB)fNDB->Fill(1);

			if(pidM==111)
			{
				if(ilabelM>=NembMCpi0 && ilabelM<NembMCeta)iEmbPi0 = kTRUE;
				if(ilabelM>=NembMCeta && ilabelM<NpureMCproc)iEmbEta = kTRUE;
			}
			if(pidM==221)
			{
				if(ilabelM>=NembMCeta && ilabelM<NpureMCproc)iEmbEta = kTRUE;
			}

			if(pidM==22) // from pi0 & eta
			{
				AliAODMCParticle* fMCparticleM = (AliAODMCParticle*) fMCarray->At(ilabelM);
				FindMother(fMCparticleM, ilabelM, pidM, pTmom);

				if(pidM==111)
				{
					if(ilabelM>=NembMCpi0 && ilabelM<NembMCeta)iEmbPi0 = kTRUE;
					if(ilabelM>=NembMCeta && ilabelM<NpureMCproc)iEmbEta = kTRUE;
				}
				if(pidM==221)
				{
					if(ilabelM>=NembMCeta && ilabelM<NpureMCproc)iEmbEta = kTRUE;
				}
			}
			fMCcheckMother->Fill(abs(pidM));
		}

		if(pidM==443)continue; // remove enhanced J/psi in MC !
		if(pidM==-99)continue; // remove e from no mother !

		if(pid_eleB || pid_eleD) {
			fHist_eff_HFE->Fill(TrkPt);
			if(fTPCnSigma>CutNsigma[0] && fTPCnSigma<CutNsigma[1]) fHist_eff_TPC->Fill(TrkPt);
		}				





		//////////////////////////////////////
		//calculate weight of photon for MC  
		//////////////////////////////////////
		Double_t WeightPho = -1.0;

		if(iEmbPi0)
		{
			if(TriggerPt=="pte5") {WeightPho = fPi005->Eval(pTmom);}
			else if(TriggerPt=="pte10"){WeightPho = fPi010->Eval(pTmom);}
			else {WeightPho = fPi000->Eval(pTmom);}
		}
		if(iEmbEta)
		{
			if(TriggerPt=="pte5") {WeightPho = fEta005->Eval(pTmom);}
			else if(TriggerPt=="pte10"){WeightPho = fEta010->Eval(pTmom);}
			else {WeightPho = fEta000->Eval(pTmom);}
		}


		AliVCluster *clustMatch=0x0;
		//cout << "EMCalIndex = " << EMCalIndex << endl;
		//if(EMCalIndex>=0)continue;
		if(EMCalIndex<0)continue;
		if(!fFlagEMCalCorrection)clustMatch = (AliVCluster*)fVevent->GetCaloCluster(EMCalIndex); // address cluster matched to track
		if(fFlagEMCalCorrection) clustMatch = dynamic_cast<AliVCluster*>(fCaloClusters_tender->At(EMCalIndex));
		fClsTypeEMC = kFALSE; fClsTypeDCAL = kFALSE;
		if(clustMatch && clustMatch->IsEMCAL())
		{
			fHistMatchPt->Fill(TrkPt);

			///////get position of clustMatch/////////
			Float_t clustMatchpos[3] = {0.};
			clustMatch->GetPosition(clustMatchpos);
			TVector3 cpos(clustMatchpos);
			Double_t Matcheta = cpos.Eta();
			Double_t Matchphi = cpos.Phi();

			///////calculate phi and eta difference between a track and a cluster//////////
			Double_t phidiff = -999;
			Double_t etadiff = -999;
			etadiff = track->GetTrackEtaOnEMCal()-Matcheta;
			phidiff = TVector2::Phi_mpi_pi(track->GetTrackPhiOnEMCal()-Matchphi);
			fEtadiff->Fill(etadiff); 
			fPhidiff->Fill(phidiff); 
			//Matchphi=TVector2::Phi_mpi_pi(Matchphi);
			if(Matchphi <0){Matchphi += 2*TMath::Pi();}


			if(TMath::Abs(etadiff)>0.05 || TMath::Abs(phidiff)>0.05) continue;
			if(Matchphi>1.39 && Matchphi < 3.265) fClsTypeEMC = kTRUE; //EMCAL : 80 < phi < 187     
			if(Matchphi>4.53 && Matchphi < 5.708) fClsTypeDCAL = kTRUE;//DCAL  : 260 < phi < 327
                        /*
			cout<< "======================================== "<<endl;
			cout<< "fFlagClsTypeEMC == "<< fFlagClsTypeEMC << "  fFlagClsTypeDCAL == "<< fFlagClsTypeDCAL<<endl;
			cout<< "fClsTypeEMC     == "<< fClsTypeEMC     << "  fClsTypeDCAL     == "<< fClsTypeDCAL<<endl;
			cout<< "MachPhi         == "<< Matchphi <<endl ;
			cout<< "======================================== "<<endl;
                        */
			//----selects EMCAL+DCAL clusters when fFlagClsTypeEMC and fFlagClsTypeDCAL is kTRUE
			if(fFlagClsTypeEMC && !fFlagClsTypeDCAL)
			{
				if(!fClsTypeEMC) continue; //selecting only EMCAL clusters
			}
			else if(fFlagClsTypeDCAL && !fFlagClsTypeEMC)
			{
				if(!fClsTypeDCAL) continue; //selecting only DCAL clusters
			}
			else{};
			//cout<<"!!!!!!!!!!!!!!!!!! DCAL !!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;

			fHistSelectPt->Fill(TrkPt);
			fHistScatter_EMcal_aftMatch->Fill(track->Eta(),track->Phi());	

			if(TrkPt>1.0){
				fHistEta_EMcal->Fill(track->Eta()); 
				fHistPhi_EMcal->Fill(track->Phi());
			}


			if(fTPCnSigma>CutNsigma[0] && fTPCnSigma<CutNsigma[1]){
				if(pid_eleB || pid_eleD) fHist_eff_match->Fill(TrkPt);
			}


			Double_t eop = -1.0;
			Double_t clE = clustMatch->E();
			Double_t m20 = clustMatch->GetM20();
			Double_t m02 = clustMatch->GetM02();
			fHistMatchE -> Fill(clE);
			if(TrkP>0)eop= clE/TrkP;

			fM02->Fill(TrkPt,m02);
			fM20->Fill(TrkPt,m20);

			if(TMath::Abs(pidM)==11){
				fM02_ele->Fill(m02);
				fM20_ele->Fill(m20);
			}else{
				fM02_had->Fill(m02);
				fM20_had->Fill(m20);
			}

			fHistNsigEop->Fill(eop,fTPCnSigma);


			Bool_t fFlagNonHFE=kFALSE; 
			Bool_t fFlagIsolation=kFALSE; 
                        Double_t IsoEnergy = -999.9;
                        Int_t NcontCone = 0;
                        Double_t IsoEnergyTrack = -999.9;
                        Int_t NtrackCone = 0;

                        Bool_t icaliso = kTRUE;
                        if(fMCarray && TMath::Abs(pdgorg)!=24 && pdgstatus!=1)icaliso = kFALSE;
                        //cout << "icaliso = " << icaliso << endl;

			//if(icaliso)IsolationCut(iTracks,track,track->Pt(),Matchphi,Matcheta,clE,fFlagNonHFE,fFlagIsolation,pid_eleB,pid_eleD, IsoEnergy);
			if(icaliso && TrkPt>15.0)IsolationCut(iTracks,track,track->Pt(),Matchphi,Matcheta,clE,fFlagNonHFE,fFlagIsolation,pid_eleB,pid_eleD, IsoEnergy, NcontCone);
			if(icaliso && TrkPt>15.0)IsolationTrackBase(iTracks, track, clE, IsoEnergyTrack, NtrackCone);
			//IsolationCut(iTracks,track,track->Pt(),Matchphi,Matcheta,clE,fFlagNonHFE,fFlagIsolation,pid_eleB,pid_eleD, IsoEnergy);
                        //cout << "IsoEnergy = " << IsoEnergy << endl;
                        //cout << "IsoEnergyTrack = " << IsoEnergyTrack << endl;

                        //if(TrkPt>10.0 && TMath::Abs(pdgorg)==24)
                        if(TrkPt>10.0 && icaliso)
                           {
                            Double_t isoarray[10];
                            isoarray[0] = TrkPt;
                            isoarray[1] = fTPCnSigma;
                            isoarray[2] = eop;
                            isoarray[3] = IsoEnergy;
                            isoarray[4] = pTpart;
                            isoarray[5] = m20;
                            isoarray[6] = m02;
                            isoarray[7] = (Double_t)NcontCone;
                            isoarray[8] = IsoEnergyTrack;
                            isoarray[9] = (Double_t)NtrackCone;
                            //cout <<"isoarray = " << isoarray[7] << endl;
                            fIsoArray->Fill(isoarray);
                           }

                        if(TrkPt>10.0 && ((pid_eleD) || (pid_eleB)))
                           {
                            Double_t isoarray[10];
                            isoarray[0] = TrkPt;
                            isoarray[1] = fTPCnSigma;
                            isoarray[2] = eop;
                            isoarray[3] = IsoEnergy;
                            isoarray[4] = pTpart;
                            isoarray[5] = m20;
                            isoarray[6] = m02;
                            isoarray[7] = (Double_t)NcontCone;
                            isoarray[8] = IsoEnergyTrack;
                            isoarray[9] = (Double_t)NtrackCone;
                            //cout <<"isoarray = " << isoarray[7] << endl;
                            fHFArray->Fill(isoarray);
                           }

                        if(fFlagIsolation && TrkPt>10.0)
                          { 
	         	   fTPCnsig_iso->Fill(TrkPt,fTPCnSigma);
                           fEop_iso->Fill(TrkPt,eop);
                          }

			//if(fTPCnSigma<6 && fTPCnSigma>-6 && eop < 1.2&& eop > 0.8 && m20>0.02 && m20<0.25){ //for MC
			if(fTPCnSigma>CutNsigma[0] && fTPCnSigma<CutNsigma[1] && eop>CutEop[0] && eop<CutEop[1] && m20>CutM20[0] && m20<CutM20[1]){
				fM02_2->Fill(TrkPt,m02);
				fM20_2->Fill(TrkPt,m20);

				///////-----Identify Non-HFE////////////////////////////
				SelectPhotonicElectron(iTracks,track,fFlagNonHFE,pidM,TrkPt,DCA[0],Bsign);
				//IsolationCut(iTracks,track,track->Pt(),Matchphi,Matcheta,clE,fFlagNonHFE,fFlagIsolation,pid_eleB,pid_eleD);
				if(fFlagIsolation)
                                    {
                                     fHistPt_Iso->Fill(track->Pt());
                                     fEop_iso_eID->Fill(TrkPt,eop);
                                    }

				if(pid_eleP)
				{
					fHistPhoReco0->Fill(track->Pt()); // reco pho
					if(iEmbPi0)fHistPhoPi0->Fill(track->Pt(),WeightPho); // reco pho
					if(iEmbEta)fHistPhoEta0->Fill(track->Pt(),WeightPho); // reco pho

					if(fFlagNonHFE)
					{
						fHistPhoReco1->Fill(track->Pt()); // reco pho
						if(iEmbPi0)fHistPhoPi1->Fill(track->Pt(),WeightPho); // reco pho
						if(iEmbEta)fHistPhoEta1->Fill(track->Pt(),WeightPho); // reco pho
					}
					else
					{
						fHistPhoReco2->Fill(track->Pt()); // org pho
					}
				}

			}



			if(fTPCnSigma<3 && fTPCnSigma>-3 && m20>CutM20[0] && m20<CutM20[1]){
				fEopPt_ele_loose -> Fill(TrkPt,eop);
			}
			if(fTPCnSigma>CutNsigma[0] && fTPCnSigma<CutNsigma[1] && m20>CutM20[0] && m20<CutM20[1]){ // TPC nsigma & shower shape cut
				fEopPt_ele_tight -> Fill(TrkPt,eop);
				if(ilabelM<NpureMC){fEopPt_ele_tight_PYTHIA -> Fill(TrkPt,eop);}
				if(pid_eleB || pid_eleD){fHist_eff_M20 -> Fill(TrkPt);}
				if(pid_eleB || pid_eleD){fEopPt_ele_tight_forSys-> Fill(TrkPt,eop);}

				if(eop>CutEop[0] && eop<CutEop[1]){ // E/p cut
					if(pid_eleB) fHistPt_HFE_MC_B -> Fill(track->Pt());
					if(pid_eleD) fHistPt_HFE_MC_D -> Fill(track->Pt());

					if(pid_eleB || pid_eleD){
						if(ilabelM<NpureMC){fHistPt_HFE_PYTHIA -> Fill(track->Pt());}
						else {fHistPt_HFE_emb -> Fill(track->Pt());}
						fHistPt_HFE_Gen -> Fill(fMCTrackpart->Pt());
						fHistPt_HFE_GenvsReco -> Fill(TrkPt,fMCTrackpart->Pt());
					}

					fHistPt_Inc->Fill(track->Pt());
					fTPCnsig_ele->Fill(TrkPt,fTPCnSigma);
					fEop_ele->Fill(eop);
					fDCAxy_Pt_ele->Fill(TrkPt,DCA[0]*Bsign*track->Charge());

					// 411 : D+, 421 :  D0, 413 : D*+, 423 : D*0, 431 : D_s+, 433 : D_s*+
					if(pid_eleD){
						if(TMath::Abs(pidM)==411 || TMath::Abs(pidM)== 413){fDCAxy_Pt_Dpm->Fill(TrkPt,DCA[0]*Bsign*track->Charge());}
						if(TMath::Abs(pidM)==421 || TMath::Abs(pidM)== 423){fDCAxy_Pt_D0->Fill(TrkPt,DCA[0]*Bsign*track->Charge());}
						if(TMath::Abs(pidM)==431 || TMath::Abs(pidM)== 433){fDCAxy_Pt_Ds->Fill(TrkPt,DCA[0]*Bsign*track->Charge());}
					}
					if(TMath::Abs(pidM)==4122){fDCAxy_Pt_lambda->Fill(TrkPt,DCA[0]*Bsign*track->Charge());}
					if(pid_eleB){fDCAxy_Pt_B->Fill(TrkPt,DCA[0]*Bsign*track->Charge());}
				}
			}
			if(fTPCnSigma< CutEopHad && m20>CutM20[0] && m20<CutM20[1]){
				fEopPt_had -> Fill(TrkPt,eop);
				fDCAxy_Pt_had->Fill(TrkPt,DCA[0]*Bsign*track->Charge());
			}
		}

		}                                         // continue until all the tracks are processed
		PostData(1, fOutputList);                           // stream the results the analysis of this event to
		// the output manager which will take care of writing
		// it to a file
}
//_____________________________________________________________________________
void AliAnalysisTaskCaloHFEpp::Terminate(Option_t *)
{
	// terminate
	// called at the END of the analysis (when all events are processed)
}
//_____________________________________________________________________________
void AliAnalysisTaskCaloHFEpp::SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec, Int_t iMC, Double_t TrkPt, Double_t DCAxy, Int_t Bsign)
{
	////// ////////////////////////////////////
	//////Non-HFE - Invariant mass method//////
	///////////////////////////////////////////

	//##################### Set cone radius  ##################### //
	Double_t CutptAsso = ptAssoMin;
	Double_t CutmassMin = massMin;
	//################################################################# //

	AliESDtrackCuts* esdTrackCutsAsso = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
	esdTrackCutsAsso->SetAcceptKinkDaughters(kFALSE);
	esdTrackCutsAsso->SetRequireTPCRefit(kTRUE);
	esdTrackCutsAsso->SetRequireITSRefit(kTRUE);
	esdTrackCutsAsso->SetEtaRange(-0.9,0.9);
	esdTrackCutsAsso->SetMaxChi2PerClusterTPC(4);
	esdTrackCutsAsso->SetMinNClustersTPC(70);
	esdTrackCutsAsso->SetMaxDCAToVertexZ(3.2);
	esdTrackCutsAsso->SetMaxDCAToVertexXY(2.4);
	esdTrackCutsAsso->SetDCAToVertex2D(kTRUE);

	Bool_t flagPhotonicElec = kFALSE;
	Bool_t fFlagEMCalCorrection = kTRUE;


	Int_t ntracks = -999;
	if(!fFlagEMCalCorrection)ntracks = fVevent->GetNumberOfTracks();
	if(fFlagEMCalCorrection) ntracks = fTracks_tender->GetEntries();

	for (Int_t jtrack = 0; jtrack < ntracks; jtrack++) {
		AliVParticle* VAssotrack = 0x0;
		if(!fFlagEMCalCorrection) VAssotrack  = fVevent->GetTrack(jtrack);
		if(fFlagEMCalCorrection) VAssotrack = dynamic_cast<AliVTrack*>(fTracks_tender->At(jtrack)); //take tracks from Tender list

		if (!VAssotrack) {
			printf("ERROR: Could not receive track %d\n", jtrack);
			continue;
		}

		AliVTrack *Assotrack = dynamic_cast<AliVTrack*>(VAssotrack);
		AliESDtrack *eAssotrack = dynamic_cast<AliESDtrack*>(VAssotrack);
		AliAODTrack *aAssotrack = dynamic_cast<AliAODTrack*>(VAssotrack);

		//------reject same track
		if(jtrack==itrack) continue;
		if(aAssotrack->Px()==track->Px() && aAssotrack->Py()==track->Py() && aAssotrack->Pz()==track->Pz())continue;

		Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
		Double_t ptAsso=-999., nsigma=-999.0, mass=-999., width = -999;
		Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;

		nsigma = fpidResponse->NumberOfSigmasTPC(Assotrack, AliPID::kElectron);
		ptAsso = Assotrack->Pt();
		Int_t chargeAsso = Assotrack->Charge();
		Int_t charge = track->Charge();
		Double_t AssoTPCchi2perNDF = aAssotrack -> Chi2perNDF();
		if(charge>0) fPDGe1 = -11;
		if(chargeAsso>0) fPDGe2 = -11;
		if(charge == chargeAsso) fFlagLS = kTRUE;
		if(charge != chargeAsso) fFlagULS = kTRUE;


		//------track cuts applied
		if(fAOD) {
			if(!aAssotrack->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
			if(aAssotrack->GetTPCNcls() < 80) continue;
			if(aAssotrack->GetITSNcls() < 3 ) continue;
			if((!(aAssotrack->GetStatus()&AliESDtrack::kITSrefit)|| (!(aAssotrack->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
		}
		else{
			if(!esdTrackCutsAsso->AcceptTrack(eAssotrack)) continue;
		}

		//-------loose cut on partner electron
		if(ptAsso <CutptAsso) continue;
		//if(ptAsso <0.2) continue;
		if(aAssotrack->Eta()<-0.9 || aAssotrack->Eta()>0.9) continue;
		if(nsigma < -3 || nsigma > 3) continue;
		if(AssoTPCchi2perNDF >= 4) continue;
		if(!(aAssotrack->GetStatus()&AliAODTrack::kITSrefit) || !(aAssotrack->GetStatus()&AliAODTrack::kTPCrefit)) continue;

		//-------define KFParticle to get mass
		AliKFParticle::SetField(fVevent->GetMagneticField());
		AliKFParticle ge1 = AliKFParticle(*track, fPDGe1);
		AliKFParticle ge2 = AliKFParticle(*Assotrack, fPDGe2);
		AliKFParticle recg(ge1, ge2);

		if(recg.GetNDF()<1) continue;
		Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
		if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;

		//-------Get mass
		Int_t MassCorrect;
		MassCorrect = recg.GetMass(mass,width);

		if(fFlagLS){
			if(mass < 0.002)cout <<"Px="<<aAssotrack->Px() <<" Py="<<aAssotrack->Py()<<" Pz="<<aAssotrack->Pz()<<endl;
			if(mass < 0.002)cout <<"Px="<<track->Px() <<" Py="<<track->Py()<<" Pz="<<track->Pz()<<endl;
			if(track->Pt()>1){
				fInv_pT_LS->Fill(TrkPt,mass);
				if(mass<CutmassMin)fDCAxy_Pt_LS->Fill(TrkPt,DCAxy*charge*Bsign);
			}
		}
		if(fFlagULS){
			if(track->Pt()>1){
				fInv_pT_ULS->Fill(TrkPt,mass);
				fInv_pT_ULS_forW->Fill(TrkPt,mass);
				if(mass<CutmassMin)fDCAxy_Pt_ULS->Fill(TrkPt,DCAxy*charge*Bsign);
			}
		}


		//if(mass<0.1 && fFlagULS && !flagPhotonicElec)
		if(mass<CutmassMin && fFlagULS && !flagPhotonicElec)
			flagPhotonicElec = kTRUE; //Tag Non-HFE (random mass cut, not optimised) 
	}
	fFlagPhotonicElec = flagPhotonicElec;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskCaloHFEpp::IsDdecay(int mpid)
{
 int abmpid = TMath::Abs(mpid);
 if(abmpid==411 || abmpid==421 || abmpid==413 || abmpid==423 || abmpid==431 || abmpid==433)
   {
    return kTRUE;
   }
 else
   {
    return kFALSE;
   } 
}
// 411 : D+, 421 :  D0, 413 : D*+, 423 : D*0, 431 : D_s+, 433 : D_s*+

//_____________________________________________________________________________
Bool_t AliAnalysisTaskCaloHFEpp::IsBdecay(int mpid)
{
 int abmpid = TMath::Abs(mpid);
 if(abmpid==511 || abmpid==521 || abmpid==513 || abmpid==523 || abmpid==531 || abmpid==533)
   {
    return kTRUE;
   }
 else
   {
    return kFALSE;
   } 
}
// 511 : B0, 521 :  B+, 513 : B*0, 523 : B*+, 531 : B_s0, 533 : B_s*0

//_____________________________________________________________________________
Bool_t AliAnalysisTaskCaloHFEpp::IsPdecay(int mpid)
{
 int abmpid = TMath::Abs(mpid);
 if(abmpid==22 || abmpid==111 || abmpid==221)
   {
		//fMCcheckMother->Fill(abs(abmpid));
    return kTRUE;
   }
 else
   {
    return kFALSE;
   } 
}


//_____________________________________________________________________________
void AliAnalysisTaskCaloHFEpp::FindMother(AliAODMCParticle* part, int &label, int &pid, double &ptmom)
{

 if(part->GetMother()>-1)
   {
    label = part->GetMother();
    AliAODMCParticle *partM = (AliAODMCParticle*)fMCarray->At(label);
    pid = partM->GetPdgCode();
		ptmom = partM->Pt();
   }
 else
   {
    pid = -99;
   } 
   //cout << "Find Mother : label = " << label << " ; pid" << pid << endl;
}

//_______________________________
 void AliAnalysisTaskCaloHFEpp::FindWdecay(AliAODMCParticle* part, Int_t &label, Int_t &pid)
 {
      while(part->GetMother()>0)
          {
           label = part->GetMother();
           //AliAODMCParticle *partM = (AliAODMCParticle*)fMCarray->At(label);
           part = (AliAODMCParticle*)fMCarray->At(label);
           pid = part->GetPdgCode();
           if(TMath::Abs(pid)==24)
              {
               Int_t mm = part->GetMother();
               cout << "pid = " << pid << " ; status "<<  part->GetStatus() << " ; " << label  <<  " ; mother = " << mm  << endl;
              }
           //cout << "mother pid = " << pid << " ; status "<<  part->GetStatus() << " ; " << label  << endl;
          }
     
  }

//_____________________________________________________________________________
void AliAnalysisTaskCaloHFEpp::CheckMCgen(AliAODMCHeader* fMCheader,Double_t CutEta)
{
 TList *lh=fMCheader->GetCocktailHeaders();
 NpureMC = 0;
 NpureMCproc = 0;
 NembMCpi0 = 0;
 NembMCeta = 0;
 Nch = 0;
 TString MCgen;
 TString embpi0("pi");
 TString embeta("eta");

 if(lh)
    {     
     for(int igene=0; igene<lh->GetEntries(); igene++)
        {
         AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(igene);
         if(gh)
           {
            MCgen =  gh->GetName();     
            if(igene==0)NpureMC = gh->NProduced();  // generate by PYTHIA or HIJING
           
            if(MCgen.Contains(embpi0))NembMCpi0 = NpureMCproc;
            if(MCgen.Contains(embeta))NembMCeta = NpureMCproc;

            NpureMCproc += gh->NProduced();  // generate by PYTHIA or HIJING
           }
        }
    }


 //for(int imc=0; imc<fMCarray->GetEntries(); imc++)
 for(int imc=0; imc<NpureMCproc; imc++)
     {
	     Bool_t iEnhance = kFALSE;
	     if(imc>=NpureMC)iEnhance = kTRUE;
	     Int_t iHijing = 1;  // select particles from Hijing or PYTHIA


	     fMCparticle = (AliAODMCParticle*) fMCarray->At(imc);
             if(!fMCparticle)continue;
	     Int_t pdgGen = TMath::Abs(fMCparticle->GetPdgCode());
	     Double_t pdgEta = fMCparticle->Eta(); 
	     Double_t pTtrue = fMCparticle->Pt(); 
	     Int_t chargetrue = fMCparticle->Charge();
	     Bool_t isPhysPrim = fMCparticle->IsPhysicalPrimary();


	     //--------- Get N charged ----------------------
	     if(chargetrue!=0){
		     if(TMath::Abs(pdgEta)<1.0){
			     if(isPhysPrim){
				     Nch++;
			     }
		     }
	     }


	     if(TMath::Abs(pdgEta)>CutEta)continue;

	     fCheckEtaMC->Fill(pdgEta);

	     Int_t pdgMom = -99;
	     Int_t pdgorg = -99;
	     Int_t labelMom = -1;
	     Double_t pTmom = -1.0;

	     FindMother(fMCparticle,labelMom,pdgMom,pTmom);
	     if(pdgMom==-99 && iEnhance)iHijing = 0;  // particles from enhance
	     if(pdgMom>0 && iEnhance)iHijing = -1;  // particles from enhance but feeddown

	     if(iHijing>-1)
	     {
		     if(pdgGen==111)fHistMCorgPi0->Fill(iHijing,pTtrue);
		     if(pdgGen==221)fHistMCorgEta->Fill(iHijing,pTtrue);
	     }

	     if(TMath::Abs(pdgGen)!=11)continue;
	     if(pTtrue<2.0)continue;

	     Int_t pdgGM = -99;
	     Int_t labelGM = -1;
	     Double_t pTGMom = -1.0;

	     if(pdgMom!=0)
	     {
		     AliAODMCParticle* fMCparticleMom = (AliAODMCParticle*) fMCarray->At(labelMom);
		     if(IsDdecay(pdgMom)){
			     fHistMCorgD->Fill(fMCparticle->Pt());
			     FindMother(fMCparticleMom,labelGM,pdgGM,pTGMom);
			     if(IsBdecay(pdgGM)){
				     fPt_Btoe->Fill(fMCparticle->Pt(),pTGMom);
			     }
		     }
		     if(IsBdecay(pdgMom)){
			     fHistMCorgB->Fill(fMCparticle->Pt());
			     fPt_Btoe->Fill(fMCparticle->Pt(),pTmom);
		     }
	     }

     }

 return;
}

//_____________________________________________________________________________
void AliAnalysisTaskCaloHFEpp::GetMClevelWdecay(AliAODMCHeader* fMCheader)
{
 //cout << "============= check W decay ============= " << endl;

 TList *lh=fMCheader->GetCocktailHeaders();
 TString MCgen;

 if(lh)
    {     
     for(int igene=0; igene<lh->GetEntries(); igene++)
        {
         AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(igene);
         if(gh)
           {
            MCgen =  gh->GetName(); 
            //cout << "MCgen = " << MCgen << " ; " << gh->NProduced() << endl;    
           }
        }
    }


 for(int imc=0; imc<fMCarray->GetEntries(); imc++)
 {

	 fMCparticle = (AliAODMCParticle*) fMCarray->At(imc);
	 Int_t pdgGen = TMath::Abs(fMCparticle->GetPdgCode());
	 Int_t pdgStatus = fMCparticle->GetStatus();

	 if(TMath::Abs(fMCparticle->Eta())>0.6)continue; 

         if(TMath::Abs(pdgGen)==11)
            {
	      // get W->e
	      Int_t ilabelM = -1;
	      Int_t pdgorg = -1;
	      FindWdecay(fMCparticle,ilabelM,pdgorg);
	      //cout << "MCcheck : pdgorg = " << pdgorg << " ; " << imc << " ; status = " << pdgStatus << endl;
	      if(TMath::Abs(pdgorg)==24 && pdgStatus==1)fHistWeOrg->Fill(fMCparticle->Pt());  // W->e(status 21) -> e(status 1) same electron 
            }
 }


}


//_____________________________________________________________________________
void AliAnalysisTaskCaloHFEpp::IsolationCut(Int_t itrack, AliVTrack *track, Double_t TrackPt, Double_t MatchPhi, Double_t MatchEta,Double_t MatchclE, Bool_t fFlagPhoto, Bool_t &fFlagIso, Bool_t fFlagB, Bool_t fFlagD, Double_t &IsoEnergy, Int_t &NcontCone)
{
	//##################### Set cone radius  ##################### //
	Double_t CutConeR = MaxConeR;
	//################################################################# //
        //cout << "cut E = " << CutMimClE << endl;
	//////////////////////////////
	// EMCal cluster loop
	//////////////////////////////
	Int_t NclustIso = -999;
	Bool_t flagIso = kFALSE;
	Bool_t fFlagEMCalCorrection = kTRUE;

	if(!fFlagEMCalCorrection)NclustIso =  fVevent->GetNumberOfCaloClusters();
	if(fFlagEMCalCorrection) NclustIso =  fCaloClusters_tender->GetEntries();

	Bool_t fClsTypeEMC = kFALSE, fClsTypeDCAL = kFALSE;;
	Double_t riso =  0.;
 	Double_t ConeR = -999.;
        Int_t NinSide = 0;

	for(Int_t jcl=0; jcl<NclustIso; jcl++)
	{
		AliVCluster *Assoclust = 0x0;     
		if(!fFlagEMCalCorrection)Assoclust = (AliVCluster*)fVevent->GetCaloCluster(jcl); // address cluster matched to track
		if(fFlagEMCalCorrection)Assoclust = dynamic_cast<AliVCluster*>(fCaloClusters_tender->At(jcl)); 

		fClsTypeEMC = kFALSE; fClsTypeDCAL = kFALSE;
		ConeR = 0.;

		if(Assoclust && Assoclust->IsEMCAL())
		{
			Float_t Assoclpos[3] = {0.};
			Assoclust->GetPosition(Assoclpos);

			TVector3 Assocpos(Assoclpos);

			Double_t AssoPhi =  Assocpos.Phi();
			if(AssoPhi <0){AssoPhi += 2*TMath::Pi();}
			Double_t AssoEta =  Assocpos.Eta();
			Double_t AssoclE =  Assoclust->E();

			//------reject same Cluster
			if(AssoclE==MatchclE && AssoPhi==MatchPhi && AssoEta==MatchEta) continue;


			if(AssoPhi > 1.39 && AssoPhi < 3.265) fClsTypeEMC = kTRUE; //EMCAL : 80 < phi < 187
			if(AssoPhi > 4.53 && AssoPhi < 5.708) fClsTypeDCAL = kTRUE;//DCAL  : 260 < phi < 327

			//----selects EMCAL+DCAL clusters when fFlagClsTypeEMC and fFlagClsTypeDCAL is kTRUE
			if(fFlagClsTypeEMC && !fFlagClsTypeDCAL)
			{
				if(!fClsTypeEMC) continue; //selecting only EMCAL clusters
			}
			if(fFlagClsTypeDCAL && !fFlagClsTypeEMC)
			{
				if(!fClsTypeDCAL) continue; //selecting only DCAL clusters
			}
			else{};


			ConeR = sqrt(pow(AssoPhi-MatchPhi,2.)+pow(AssoEta-MatchEta,2.));
			if(ConeR>CutConeR) continue;
			if(AssoclE<CutMimClE) continue;

			riso += AssoclE;
                        NinSide++;
			fConeR->Fill(TrackPt,ConeR);
			fConeE->Fill(TrackPt,AssoclE);
		}
	}

	fNpart->Fill(TrackPt,NinSide);

	riso = riso/MatchclE;
	if(riso>=0.0){
		fHistPt_R_Iso->Fill(TrackPt,riso);
		if(fFlagB || fFlagD){fHist_eff_Iso -> Fill(TrackPt,riso);}
	}

	//if(TrackPt >= 30. && riso>=0.0)CheckCorrelation(itrack,track,TrackPt,riso,fFlagPhoto);

	if(riso<0.05 && riso>=0.0) flagIso = kTRUE;
        //cout << "riso = " << riso << endl;
        //cout << "NinSide = " << NinSide << endl;
	fFlagIso = flagIso;
        IsoEnergy = riso;
        NcontCone = NinSide;

}

//_____________________________________________________________________________
void AliAnalysisTaskCaloHFEpp::IsolationTrackBase(Int_t itrack, AliVTrack *track, Double_t MatchclE, Double_t &IsoEnergyTrack, Int_t &NtrackCone)
{

	//##################### Set cone radius  ##################### //
	Double_t CutConeR = MaxConeR;
	//################################################################# //

	Int_t nWassotracks = -999;
	nWassotracks = fTracks_tender->GetEntries();

        Double_t risoTrack = 0.0; 

	//////////////////////////////
	// Track loop
	//////////////////////////////
	for (Int_t jtrack = 0; jtrack < nWassotracks; jtrack++) {
		AliVParticle* VWassotrack = 0x0;
		VWassotrack = dynamic_cast<AliVTrack*>(fTracks_tender->At(jtrack)); //take tracks from Tender list

		if (!VWassotrack) {
			printf("ERROR: Could not receive track %d\n", jtrack);
			continue;
		}

		AliVTrack   *Wassotrack  = dynamic_cast<AliVTrack*>(VWassotrack);
		AliAODTrack *aWassotrack = dynamic_cast<AliAODTrack*>(VWassotrack);

		if(!aWassotrack) continue;                            // if we failed, skip this 

		//------reject same track
		if(jtrack==itrack) continue;
		if(aWassotrack->Px()==track->Px() && aWassotrack->Py()==track->Py() && aWassotrack->Pz()==track->Pz())continue;

                //------ find tracks around candidate
		Double_t ptWasso = -999., phiWasso = -999., etaWasso = -999.;
		Double_t TrackPhi = -999., TrackEta = -999.;

		ptWasso         = aWassotrack->Pt();
		phiWasso        = aWassotrack->Phi();
		etaWasso        = aWassotrack->Eta();
		TrackPhi        = track->Phi();
		TrackEta        = track->Eta();

		if(ptWasso <0.15) continue;
               
		Double_t Wphidiff = phiWasso - TrackPhi; 
		Wphidiff = TMath::ATan2(TMath::Sin(Wphidiff),TMath::Cos(Wphidiff)); 
		if(Wphidiff < -TMath::Pi()/2) Wphidiff += 2*TMath::Pi();

                Double_t Wetadiff = etaWasso - TrackEta;

                Double_t ConeRtr = sqrt(pow(Wetadiff,2)+pow(Wphidiff,2));
                //cout << "ConeRtr = " << ConeRtr << endl;

               if(ConeRtr>CutConeR) continue;

	       Int_t EMCalIndex_TrCone = aWassotrack->GetEMCALcluster();  // get index of EMCal cluster which matched to track
               if(EMCalIndex_TrCone<0)continue;
	       AliVCluster *Assoclust_TrCone = 0x0;     
	       Assoclust_TrCone = dynamic_cast<AliVCluster*>(fCaloClusters_tender->At(EMCalIndex_TrCone)); 

	       risoTrack += Assoclust_TrCone->E();
               //cout << "risoTrack = " << risoTrack << endl;

               NtrackCone++;

	}  // end track loop

        //cout << "<----- risoTrack = " << risoTrack << endl;
        //cout << "<----- MatchclE = " << MatchclE << endl;
        
        IsoEnergyTrack = risoTrack/MatchclE;

}

//_____________________________________________________________________________
void AliAnalysisTaskCaloHFEpp::CheckCorrelation(Int_t itrack, AliVTrack *track, Double_t TrackPt, Double_t Riso, Bool_t fFlagPhoto)
{

	//##################### Systematic Parameters ##################### //
	//---Track Cut
	Double_t CutTrackEtaW[2] = {TrackEtaMin,TrackEtaMax};
	Int_t CutTPCNClsW = NTPCClust;
	Int_t CutITSNClsW = NITSClust; 
	Int_t CutTPCNCrossedRowW = NCrossedRow;
	Double_t CutDCAxyW = DCAxy;
	Double_t CutDCAzW  = DCAz;
	//################################################################# //

	Bool_t fFlagEMCalCorrection = kTRUE;

	Int_t nWassotracks = -999;
	if(!fFlagEMCalCorrection)nWassotracks = fVevent->GetNumberOfTracks();
	if(fFlagEMCalCorrection) nWassotracks = fTracks_tender->GetEntries();

	const AliVVertex *WpVtx = fVevent->GetPrimaryVertex();

	//////////////////////////////
	// Track loop
	//////////////////////////////
	for (Int_t jtrack = 0; jtrack < nWassotracks; jtrack++) {
		AliVParticle* VWassotrack = 0x0;
		if(!fFlagEMCalCorrection) VWassotrack  = fVevent->GetTrack(jtrack);
		if(fFlagEMCalCorrection)  VWassotrack = dynamic_cast<AliVTrack*>(fTracks_tender->At(jtrack)); //take tracks from Tender list

		if (!VWassotrack) {
			printf("ERROR: Could not receive track %d\n", jtrack);
			continue;
		}

		AliVTrack   *Wassotrack  = dynamic_cast<AliVTrack*>(VWassotrack);
		AliESDtrack *eWassotrack = dynamic_cast<AliESDtrack*>(VWassotrack);
		AliAODTrack *aWassotrack = dynamic_cast<AliAODTrack*>(VWassotrack);

		if(!aWassotrack) continue;                            // if we failed, skip this track

		//------reject same track
		if(jtrack==itrack) continue;
		if(aWassotrack->Px()==track->Px() && aWassotrack->Py()==track->Py() && aWassotrack->Pz()==track->Pz())continue;

		Double_t ptWasso = -999., phiWasso = -999., etaWasso = -999.;
		Double_t ITSchi2Wasso = -999., TPCchi2NDFWasso = -999.;
		Double_t Wphidiff = -999., TrackPhi = -999.;

		ptWasso         = aWassotrack->Pt();
		phiWasso        = aWassotrack->Phi();
		etaWasso        = aWassotrack->Eta();
		ITSchi2Wasso    = aWassotrack -> GetITSchi2();
		TPCchi2NDFWasso = aWassotrack -> Chi2perNDF();
		TrackPhi        = track->Phi();

		if(ptWasso <1.) continue;

		/////////////////////////
		// track cut
		/////////////////////////
		//==== 1.TPC and ITS refit cut ====
		if(!(aWassotrack->GetStatus()&AliAODTrack::kITSrefit) || !(aWassotrack->GetStatus()&AliAODTrack::kTPCrefit)) continue;
		//==== 2.AOD filter bit required ====
		if(!aWassotrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue; //mimimum cuts
		//==== 3.TPC cluster cut ====
		if(aWassotrack->GetTPCNcls() < CutTPCNClsW) continue; 
		//==== 4.ITS cluster cut ====
		if(aWassotrack->GetITSNcls() < CutITSNClsW) continue;  
		//==== 5.SPD hit cut ====
		if(!(aWassotrack -> HasPointOnITSLayer(0) || aWassotrack -> HasPointOnITSLayer(1))) continue;
		//==== 6.Eta cut ====
		if(etaWasso>CutTrackEtaW[1] || etaWasso<CutTrackEtaW[0]) continue; 
		//==== 7.DCA cut ====
		Double_t WassoDCA[2] = {-999.,-999.}, Wassocovar[3];
		if(aWassotrack -> PropagateToDCA(WpVtx,fVevent -> GetMagneticField(),20.,WassoDCA,Wassocovar))
		{
			if(TMath::Abs(WassoDCA[0]) > CutDCAxyW || TMath::Abs(WassoDCA[1]) > CutDCAzW) continue;
		}
		//==== 8.chi2 cut ====
		if((ITSchi2Wasso >= 25) || (TPCchi2NDFWasso >= 4)) continue;
		//==== 9.NCrossedRow cut ====
		if(aWassotrack -> GetTPCCrossedRows() < CutTPCNCrossedRowW) continue;


		Wphidiff = phiWasso - TrackPhi; 
		Wphidiff = TMath::ATan2(TMath::Sin(Wphidiff),TMath::Cos(Wphidiff)); 
		if(Wphidiff < -TMath::Pi()/2) Wphidiff += 2*TMath::Pi();

		fRiso_phidiff -> Fill(Wphidiff,Riso);
		if(!fFlagPhoto)fRiso_phidiff_LS -> Fill(Wphidiff,Riso);

		if(TrackPt>=35.){
			fRiso_phidiff_35 -> Fill(Wphidiff,Riso);
			if(!fFlagPhoto)fRiso_phidiff_LS_35 -> Fill(Wphidiff,Riso);
		}

	}

}
//____________________________________________________________________________
TProfile* AliAnalysisTaskCaloHFEpp::GetEstimatorHistogram(const AliAODEvent* fAOD)
{
    
  Int_t runNo  = fAOD->GetRunNumber();
  Int_t period = -1; 
   
  if (runNo>=255539 && runNo<=255618) period = 0; //LHC16i
  if (runNo>=256207 && runNo<=256420) period = 1; //LHC16j
  if (runNo>=256941 && runNo<=258537) period = 2; //LHC16k
  if (runNo>=262424 && runNo<=263487) period = 3; //LHC16o
  if (runNo>=271868 && runNo<=273100) period = 4; //LHC17h
  if (runNo>=273591 && runNo<=274442) period = 5; //LHC17i
  if (runNo>=274801 && runNo<=276508) period = 6; //LHC17k
  if (runNo>=276551 && runNo<=278216) period = 7; //LHC17l
  if (runNo>=278915 && runNo<=280140) period = 8; //LHC17m
  if (runNo>=280282 && runNo<=281961) period = 9; //LHC17o
  if (runNo>=282544 && runNo<=282704) period = 10; //LHC17r
  if (period < 0 || period > 10) return 0;
    
    
  return fMultEstimatorAvg[period];
}
//____________________________________________________________________________
TProfile* AliAnalysisTaskCaloHFEpp::GetEstimatorHistogramMC(const AliAODEvent* fAOD)
{
    
  Int_t runNo  = fAOD->GetRunNumber();
  Int_t period = -1; 
   
	if (runNo>=256504 && runNo<=258537) period = 11;  //LHC16k
  if (runNo>=258919 && runNo<=259888) period = 12; //LHC16l
  if (period < 11 || period > 12) return 0;
    
    
  return fMultEstimatorAvg[period];
}
 //____________________________________________________________________________
Double_t AliAnalysisTaskCaloHFEpp::GetCorrectedNtrackletsD(TProfile* estimatorAvg, Double_t uncorrectedNacc, Double_t vtxZ, Double_t refMult)
{
	if(TMath::Abs(vtxZ)>10.0){
		//    printf("ERROR: Z vertex out of range for correction of multiplicity\n");
		return uncorrectedNacc;
	}

	if(!estimatorAvg){
		printf("ERROR: Missing TProfile for correction of multiplicity\n");
		return uncorrectedNacc;
	}

	Double_t localAvg = estimatorAvg->GetBinContent(estimatorAvg->FindBin(vtxZ));
	Double_t deltaM = 0;
	deltaM = uncorrectedNacc*(refMult/localAvg - 1);

	Double_t correctedNacc = uncorrectedNacc + (deltaM>0 ? 1 : -1) * gRandom->PoissonD(TMath::Abs(deltaM));

	if(correctedNacc<0) correctedNacc=0;

	return correctedNacc;

}

