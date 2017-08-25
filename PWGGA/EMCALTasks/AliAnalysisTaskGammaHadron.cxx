//
// Task to estimate the number of gamma-hadron
// statistic available in the Pb+Pb run.
//
// Author: E. Epple, based on code by  B. Sahlmueller and C. Loizides

#include <Riostream.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliAnalysisTaskGammaHadron.h"

#include "AliParticleContainer.h"
#include "AliClusterContainer.h"
#include "AliAnalysisManager.h"
#include "AliCentrality.h"
#include "AliVCluster.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliAODTrack.h"
#include "AliVEventHandler.h"
#include "AliInputEventHandler.h"
#include "AliAODEvent.h"
#include "AliExternalTrackParam.h"
#include "AliTrackerBase.h"
#include "AliLog.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALGeoParams.h"
#include "AliPicoTrack.h"
#include "AliVVZERO.h"
#include "AliESDUtils.h"
#include "AliEventPoolManager.h"
#include "AliMCEvent.h"


using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskGammaHadron)
//
//  This class inherits from AliAnalysisTaskEmcal() -->fcent,fVertex,fNVertCont,fVertexSPD,fNVertSPDCont,fTriggers is defined in there
//  And AliAnalysisTaskEmcal() inherits from AliAnalysisTaskSE()
//________________________________________________________________________
AliAnalysisTaskGammaHadron::AliAnalysisTaskGammaHadron():
AliAnalysisTaskEmcal("AliAnalysisTaskGammaHadron", kTRUE),

fGammaOrPi0(0),fDoMixing(0),fMCorData(0),fDebug(0),fSavePool(0),
fEventCuts(0),fFiducialCuts(0x0),fFiducialCellCut(0x0), fHistEffGamma(0x0),fHistEffHadron(0x0),
fRtoD(0),
fClShapeMin(0),fClShapeMax(10),fMaxNLM(10),fRmvMTrack(0),fTrackMatchEta(0),fTrackMatchPhi(0),
fMixBCent(0),fMixBZvtx(),fPoolMgr(0x0),fTrackDepth(0),fPoolSize(0),fEventPoolOutputList(0),
fTriggerType(AliVEvent::kINT7), fMixingEventType(AliVEvent::kINT7),fCurrentEventTrigger(0),
fParticleLevel(kFALSE),fIsMC(kFALSE),
fEventCutList(0),fOutputList1(0),fOutputListTrAs(0),fOutputListGamma(0),fOutputListXi(0),fOutputListZeta(0),

fHistNoClusPt(0),fHistPi0(0),fHistEvsPt(0),fHistClusPairInvarMasspT(0),fMAngle(0),fPtAngle(0),fHistBinCheckPt(0),fHistBinCheckZt(0),fHistBinCheckXi(0),

fHistMatchEtaPhiAllCl2(0),fHistMatchEtaPhiAllCl3(0),fHistMatchEtaPhiAllCl4(0),
fHistDEtaDPhiGammaQA(0),fHistDEtaDPhiTrackQA(0),
fHistCellsCluster(0),fHistClusterTime(0),

//fAODfilterBits(0),fHistptAssHadronG(0),fHistptAssHadronZt(0),fHistptAssHadronXi(0),fHistDEtaDPhiG(0),fHistDEtaDPhiZT(0),fHistDEtaDPhiXI(0)
//fHistptTriggG(),fHistptTriggZt(),fHistptTriggXi(),
fCorrVsManyThings(0),fCorrVsManyThingsME(0), fClusterProp(0),
fHPoolReady(0x0)
{
	//..Initialize by defult for
	//..AliAnalysisTaskGammaHadron(0,0);
	InitArrays();

}
//________________________________________________________________________
AliAnalysisTaskGammaHadron::AliAnalysisTaskGammaHadron(Bool_t InputGammaOrPi0,Bool_t InputDoMixing, Bool_t InputMCorData):
AliAnalysisTaskEmcal("AliAnalysisTaskGammaHadron", kTRUE),

fGammaOrPi0(0),fDoMixing(0),fMCorData(0),fDebug(0),fSavePool(0),
fEventCuts(0),fFiducialCuts(0x0),fFiducialCellCut(0x0),fHistEffGamma(0x0),fHistEffHadron(0x0),
fRtoD(0),
fClShapeMin(0),fClShapeMax(10),fMaxNLM(10),fRmvMTrack(0),fTrackMatchEta(0),fTrackMatchPhi(0),
fMixBCent(0),fMixBZvtx(),fPoolMgr(0x0),fTrackDepth(0),fPoolSize(0),fEventPoolOutputList(0),
fTriggerType(AliVEvent::kINT7), fMixingEventType(AliVEvent::kINT7),fCurrentEventTrigger(0),
fParticleLevel(kFALSE),fIsMC(kFALSE),
fEventCutList(0),fOutputList1(0),fOutputListTrAs(0),fOutputListGamma(0),fOutputListXi(0),fOutputListZeta(0),

fHistNoClusPt(0),fHistPi0(0),fHistEvsPt(0),fHistClusPairInvarMasspT(0),fMAngle(0),fPtAngle(0),fHistBinCheckPt(0),fHistBinCheckZt(0),fHistBinCheckXi(0),

fHistMatchEtaPhiAllCl2(0),fHistMatchEtaPhiAllCl3(0),fHistMatchEtaPhiAllCl4(0),
fHistDEtaDPhiGammaQA(0),fHistDEtaDPhiTrackQA(0),
fHistCellsCluster(0),fHistClusterTime(0),

//fAODfilterBits(0),fHistptAssHadronG(0),fHistptAssHadronZt(0),fHistptAssHadronXi(0),fHistDEtaDPhiG(0),fHistDEtaDPhiZT(0),fHistDEtaDPhiXI(0)
//fHistptTriggG(),fHistptTriggZt(),fHistptTriggXi(),
fCorrVsManyThings(0), fCorrVsManyThingsME(0), fClusterProp(0),
fHPoolReady(0x0)
{
	InitArrays();
	//..set input variables
	fGammaOrPi0        =InputGammaOrPi0;
	fDoMixing          =InputDoMixing;
	fMCorData          =InputMCorData;
}
//________________________________________________________________________
void AliAnalysisTaskGammaHadron::InitArrays()
{
	//..Initialize by defult for
	//..AliAnalysisTaskGammaHadron(0,1);

	//..set input variables
	fGammaOrPi0        =0; //= 0 ( Gamma analysis ), 1 (pi0 analysis)
	fDoMixing          =0; //= 0 (do only same event analyis with correct triggers), =1 (do event mixing)
	fMCorData          =0; // 0->MC, 1->Data

	fPlotQA            =0; //= 0 do not plot additional QA histograms, =1 plot them
	fDebug             =0; //set only 1 for debugging
	fSavePool          =0; //= 0 do not save the pool by default. Use the set function to do this.
	fUseManualEventCuts=0; //= 0 use automatic setting from AliEventCuts. =1 load manual cuts

	//..These two items are set in AliAnalysisTaskEmcal::RetrieveEventObjects()
	//fCent, zVertex

	//..Initialize the arrays to 0
	/*for(Int_t i=0; i<kNIdentifier;i++)
	{
		fHistptAssHadronG[i] =0;
		fHistptAssHadronZt[i]=0;
		fHistptAssHadronXi[i]=0;
		fHistptTriggG[i] =0;
		fHistptTriggZt[i]=0;
		fHistptTriggXi[i]=0;


		for(Int_t j=0; j<kNIdentifier;j++)
		{
			if(j<kNoGammaBins+1)fHistDEtaDPhiG[i][j]  = 0;
			if(j<kNoZtBins+1)   fHistDEtaDPhiZT[i][j] = 0;
			if(j<kNoXiBins+1)   fHistDEtaDPhiXI[i][j] = 0;
		}
	}*/

	fRtoD=180.0/TMath::Pi();

	fFiducialCuts    = new AliFiducialCut();
	fFiducialCellCut = new AliEMCALRecoUtils();
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   Define bins in which the 2D histograms are plotted
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Double_t fZtStep =1.0/(7-1.0);  // Bin width for the zT histograms
	Double_t fXiStep =2.5/(8-1.0);  // Bin width for the Xi histograms

	Double_t fArray_G_BinsValue[kNoGammaBins+1] ={5,7,9,11,14,17,22,30,60,90};
	Double_t fArray_ZT_BinsValue[kNoZtBins+1]   ={0,fZtStep,2*fZtStep,3*fZtStep,4*fZtStep,5*fZtStep,6*fZtStep,20};
	Double_t fArray_XI_BinsValue[kNoXiBins+1]   ={-10,0,fXiStep,2*fXiStep,3*fXiStep,4*fXiStep,5*fXiStep,6*fXiStep,10};

	memcpy (fArray_G_Bins,  fArray_G_BinsValue,  sizeof (fArray_G_Bins));
	memcpy (fArray_ZT_Bins, fArray_ZT_BinsValue, sizeof (fArray_ZT_Bins));
	memcpy (fArray_XI_Bins, fArray_XI_BinsValue, sizeof (fArray_XI_Bins));
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   Define vertex and centrality bins for the ME background
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//..if desired one can add a set function to set these values in the add task function
	Double_t centmix[kNcentBins+1] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0};
	fMixBCent = new TAxis(kNcentBins,centmix);

	//static const Int_t NvertBins=8;
	//Double_t zvtxmix[NvertBins+1] = {-10,-6,-4,-2,0,2,4,6,10};
	Double_t zvtxmix[kNvertBins+1] = {-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10};
	memcpy (fArrayNVertBins, zvtxmix, sizeof (fArrayNVertBins));
	fMixBZvtx = new TAxis(kNvertBins,zvtxmix);

	//..Raymond/Megan gives more mixed event yield - don't know about the quality though
	//fTrackDepth     = 100;      //Hanseul sets it to 100! Q:: is this good? Maximum number of tracks??
	fTrackDepth     = 50000;    //Raymonds/Megans value

	//..!!
	//.. fPoolSize is an input that is ignored in the PoolManager Anyway
	fPoolSize       = 1;     //1000 - Raymond/Megan value, says it is ignored anyway

	//..member function of AliAnalysisTaskEmcal
	SetMakeGeneralHistograms(kTRUE);
}
//________________________________________________________________________
AliAnalysisTaskGammaHadron::~AliAnalysisTaskGammaHadron()
{
	// Destructor

	//Copied from chris yaldo. Ask Salvatore about it!
	// Destructor. Clean-up the output list, but not the histograms that are put inside
	// (the list is owner and will clean-up these histograms). Protect in PROOF case.
	if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode())
	{
		delete fOutputList1;
	}
	//copied from hanseul
	/*if (fPoolMgr)
	{
		delete fPoolMgr;
	}*/
}
//________________________________________________________________________
void AliAnalysisTaskGammaHadron::UserCreateOutputObjects()
{
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::UserCreateOutputObjects()"<<endl;

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   Create fOutput list and histograms (fHistZVertex, fHistEventRejection, fHistEventRejection, fHistEventCount, fHistCentrality, fHistEventPlane)
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	AliAnalysisTaskEmcal::UserCreateOutputObjects();

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   Add output for AliEventCuts
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	fEventCutList = new TList();
	fEventCutList ->SetOwner();
	fEventCutList ->SetName("EventCutOutput");

	fEventCuts.OverrideAutomaticTriggerSelection(fOffTrigger); //..otherwise only kINT7 events are used for the analysis
	if(fUseManualEventCuts==1)
	{
	    //..Enable manual mode.
		//..this just means that the automatic cut settings
	    //..are not loaded every time the event is checked
	    fEventCuts.SetManualMode();
		/*
	     This is not possible because these are private functions!
	    //..nevertheless we set now the standard settings
	    //..for the respective period and then overwrite
	    //..some cuts with the set values in the Emcal task.
	    fEventCuts.fCurrentRun = fRunNumber;
	    fEventCuts.AutomaticSetup();
		 */
		//..overwrite the manual set cuts with
		//..some of our own values
	    fEventCuts.fCentralityFramework=2; //..only for Run1!!
		fEventCuts.fTriggerMask = fOffTrigger;
		fEventCuts.fMinVtz = fMinVz;
		fEventCuts.fMaxVtz = fMaxVz;
		fEventCuts.fRequireTrackVertex = true;
		fEventCuts.fMaxDeltaSpdTrackAbsolute=fZvertexDiff;
		fEventCuts.fTrackletBGcut = fTklVsClusSPDCut; //(false by default for 15o)
		fEventCuts.fMinCentrality = fMinCent;
		fEventCuts.fMaxCentrality = fMaxCent;
		//++fRejectPileup (IsPileupFromSPD)= true (fixed in code)
		//+remove multi vertexer pile up (false - not activated yet)
		//+spd vertex resolution etc
		//+some cent. resolution cuts
		//+some variable correlatios - fixed to false
	}

	fEventCuts.AddQAplotsToList(fEventCutList);
	fOutput->Add(fEventCutList);
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   Create mixed event pools
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	if(fDoMixing==1 || fPoolMgr) //do this for either a mixed event analysis or when an external pool is given
	{
		InitEventMixer();
	}
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   Define sublists/folders for a better organisation of the figures
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	fOutputList1    = new TList();
	fOutputList1    ->SetOwner();
	fOutputList1    ->SetName("pT_distributions_of_the_gamma");
/*	fOutputListTrAs = new TList();
	fOutputListTrAs ->SetOwner();
	fOutputListTrAs ->SetName("TriggAndAssoc");
	fOutputListGamma= new TList();
	fOutputListGamma->SetOwner();
	fOutputListGamma->SetName("Different_Gamma_2DHistograms");
	fOutputListXi   = new TList();
	fOutputListXi   ->SetOwner();
	fOutputListXi   ->SetName("Different_Xi_2DHistograms");
	fOutputListZeta = new TList();
	fOutputListZeta ->SetOwner();
	fOutputListZeta ->SetName("Different_Zt_2DHistograms");
*/	fOutputListQA   = new TList();
	fOutputListQA   ->SetOwner();
	fOutputListQA   ->SetName("QA_histograms");

	//common bins for the histograms
	Int_t nbins[6]  = {0};
	Double_t min[6] = {0};
	Double_t max[6] = {0};

	//settings for p_t cluster distributon
	nbins[0] = 500;
	min[0] = 0;
	max[0] = 100;
	//settings for p_t hadron distribution
	nbins[1] = 60;  //do 1/2 GeV bins so that you can see the 0.5 cut to set as a minimum pT to combine hadron and gamma
	min[1] = 0;
	max[1] = 30;
	//settings for delta phi (g-h) distribution
	nbins[2] = 45;
	min[2] = -90;
	max[2] = 270;
	//settings for delta eta (g-h) distribution
	nbins[3] = 80;
	min[3] = -2;
	max[3] = 2;
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//
	//   Create Histograms
	//
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//    Histograms for common use
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	//..Initialize
	fHistBinCheckPt       = new TH1*[kNIdentifier];
	fHistBinCheckZt       = new TH1*[kNIdentifier];
	fHistBinCheckXi       = new TH1*[kNIdentifier];
	fHistDEtaDPhiGammaQA  = new TH2*[kNIdentifier+3]; //Why +2?? -> because I want more than 3 QA versions
	fHistDEtaDPhiTrackQA  = new TH2*[kNIdentifier+3];
	fHistCellsCluster     = new TH2*[kNIdentifier+3];
	fHistClusterTime      = new TH2*[kNIdentifier+3];

	/*for(Int_t i=0; i<kNIdentifier; i++)
	{
		fHistptAssHadronG[i] = new TH1*[kNoGammaBins];
		fHistptAssHadronZt[i]= new TH1*[kNoZtBins];
		fHistptAssHadronXi[i]= new TH1*[kNoXiBins];
		fHistptTriggG[i]     = new TH1*[kNoGammaBins];
		fHistptTriggZt[i]    = new TH1*[kNoZtBins];
		fHistptTriggXi[i]    = new TH1*[kNoXiBins];


		for(Int_t j=0; j<kNoGammaBins;j++)
		{
			if(j<kNoGammaBins+1)fHistDEtaDPhiG[i][j] = new TH2*[kNvertBins+1];
			if(j<kNoZtBins+1)   fHistDEtaDPhiZT[i][j]= new TH2*[kNvertBins+1];
			if(j<kNoXiBins+1)   fHistDEtaDPhiXI[i][j]= new TH2*[kNvertBins+1];
		}
	}*/
	//.................................
	//..p_T^{Cluster} distribution under different conditions

	//..all clusters
	fHistNoClusPt = new TH1F(Form("fHistNoClusPt_Id%0d",1),Form("fHistNoClusPt_Id%0d",1), 31,0, 31);
	fHistNoClusPt->GetXaxis()->SetTitle("p_{T}^{Calo Cluster}");
	fHistNoClusPt->GetYaxis()->SetTitle(Form("No. of Clusters [counts/%0.1f GeV/c]",fHistNoClusPt->GetBinWidth(0)));
	fOutputList1->Add(fHistNoClusPt);


	//..by the identifier different histograms can be filled under different cut conditions
	//..Can eg. later be modified to contain certain delta phi or centrality bins
	for(Int_t identifier=0;identifier<kNIdentifier;identifier++)
	{
		fHistBinCheckPt[identifier] = new TH1F(Form("fHistBinCheckPt_%0d",identifier),Form("fHistBinCheckPt_%0d",identifier), nbins[0], min[0], max[0]);
		fHistBinCheckPt[identifier]->GetXaxis()->SetTitle("p_{T}^{#gamma}");
		fHistBinCheckPt[identifier]->GetYaxis()->SetTitle("Entries");
		fOutput->Add(fHistBinCheckPt[identifier]);

		fHistBinCheckZt[identifier] = new TH1F(Form("fHistBinCheckZt_%0d",identifier),Form("fHistBinCheckZt_%0d",identifier), 1500, 0, 60);
		fHistBinCheckZt[identifier]->GetXaxis()->SetTitle("z_{T}^{#gamma-h}");
		fHistBinCheckZt[identifier]->GetYaxis()->SetTitle("Entries");
		fOutput->Add(fHistBinCheckZt[identifier]);

		fHistBinCheckXi[identifier] = new TH1F(Form("fHistBinCheckXi_%0d",identifier),Form("fHistBinCheckXi_%0d",identifier), 500, -20, 20);
		fHistBinCheckXi[identifier]->GetXaxis()->SetTitle("#xi^{#gamma-h}");
		fHistBinCheckXi[identifier]->GetYaxis()->SetTitle("Entries");
		fOutput->Add(fHistBinCheckXi[identifier]);

		/*
		for(Int_t i=0; i<kNoGammaBins; i++)
		{
			for(Int_t j=0; j<kNvertBins; j++)
			{
				fHistDEtaDPhiG[identifier][i][j] = new TH2F(Form("fHistDEtaDPhiG%d_Id%d_V%d",i,identifier,j),Form("fHistDEtaDPhiG%d_Id%d_V%d",i,identifier,j),nbins[2],min[2],max[2],nbins[3],min[3],max[3]);
				fHistDEtaDPhiG[identifier][i][j]->GetXaxis()->SetTitle(Form("#Delta #varphi^{#gamma-h} %0.1f<p_{T}^{#gamma}<%0.1f",fArray_G_Bins[i],fArray_G_Bins[i+1]));
				fHistDEtaDPhiG[identifier][i][j]->GetYaxis()->SetTitle("#Delta #eta^{#gamma-h}");
				fOutputListGamma->Add(fHistDEtaDPhiG[identifier][i][j]);
			}
			fHistptAssHadronG[identifier][i] = new TH1F(Form("fHistPtAssHadronG%d_Id%d",i,identifier),Form("fHistPtAssHadronG%d_Id%d",i,identifier), nbins[1], min[1], max[1]);
			fHistptAssHadronG[identifier][i]->GetXaxis()->SetTitle(Form("p_{T}^{assoc. h} %0.1f<p_{T}^{#gamma}<%0.1f",fArray_G_Bins[i],fArray_G_Bins[i+1]));
			fHistptAssHadronG[identifier][i]->GetYaxis()->SetTitle(Form("dN^{assoc. h}/dp_{T}^{h} [counts/%0.1f GeV/c]",fHistptAssHadronG[identifier][i]->GetBinWidth(1)));
			fOutputListTrAs->Add(fHistptAssHadronG[identifier][i]);

			fHistptTriggG[identifier][i] = new TH1F(Form("fHistptTriggG%d_Id%d",i,identifier),Form("fHistptTriggG%d_Id%d",i,identifier), nbins[0], min[0], max[0]);
			fHistptTriggG[identifier][i]->GetXaxis()->SetTitle(Form("p_{T}^{cluster} %0.1f<p_{T}^{#gamma}<%0.1f",fArray_G_Bins[i],fArray_G_Bins[i+1]));
			fHistptTriggG[identifier][i]->GetYaxis()->SetTitle(Form("dN^{assoc. h}/dp_{T}^{h} [counts/%0.1f GeV/c]",fHistptAssHadronG[identifier][i]->GetBinWidth(1)));
			fOutputListTrAs->Add(fHistptTriggG[identifier][i]);
		}
		for(Int_t i=0; i<kNoZtBins; i++)
		{
			for(Int_t j=0; j<kNvertBins; j++)
			{
				fHistDEtaDPhiZT[identifier][i][j] = new TH2F(Form("fHistDEtaDPhiZT%d_Id%d_V%d",i,identifier,j),Form("fHistDEtaDPhiZT%d_Id%d_V%d",i,identifier,j),nbins[2],min[2],max[2],nbins[3],min[3],max[3]);
				fHistDEtaDPhiZT[identifier][i][j]->GetXaxis()->SetTitle(Form("#Delta #varphi^{#gamma-h} %0.1f<z_{T}<%0.1f",fArray_ZT_Bins[i],fArray_ZT_Bins[i+1]));
				fHistDEtaDPhiZT[identifier][i][j]->GetYaxis()->SetTitle("#Delta #eta^{#gamma-h}");
				fOutputListZeta->Add(fHistDEtaDPhiZT[identifier][i][j]);
			}
			fHistptAssHadronZt[identifier][i] = new TH1F(Form("fHistPtAssHadronZt%d_Id%d",i,identifier),Form("fHistPtAssHadronZt%d_Id%d",i,identifier), nbins[1], min[1], max[1]);
			fHistptAssHadronZt[identifier][i]->GetXaxis()->SetTitle(Form("p_{T}^{assoc. h} %0.1f<z_{T}<%0.1f",fArray_ZT_Bins[i],fArray_ZT_Bins[i+1]));
			fHistptAssHadronZt[identifier][i]->GetYaxis()->SetTitle(Form("dN^{assoc. h}/dp_{T}^{h} [counts/%0.1f GeV/c]",fHistptAssHadronZt[identifier][i]->GetBinWidth(1)));
			fOutputListTrAs->Add(fHistptAssHadronZt[identifier][i]);

			fHistptTriggZt[identifier][i] = new TH1F(Form("fHistptTriggZt%d_Id%d",i,identifier),Form("fHistptTriggZt%d_Id%d",i,identifier), nbins[0], min[0], max[0]);
			fHistptTriggZt[identifier][i]->GetXaxis()->SetTitle(Form("p_{T}^{cluster} %0.1f<z_{T}<%0.1f",fArray_ZT_Bins[i],fArray_ZT_Bins[i+1]));
			fHistptTriggZt[identifier][i]->GetYaxis()->SetTitle(Form("dN^{assoc. h}/dp_{T}^{h} [counts/%0.1f GeV/c]",fHistptTriggZt[identifier][i]->GetBinWidth(1)));
			fOutputListTrAs->Add(fHistptTriggZt[identifier][i]);
		}
		for(Int_t i=0; i<kNoXiBins; i++)
		{
			for(Int_t j=0; j<kNvertBins; j++)
			{
				fHistDEtaDPhiXI[identifier][i][j] = new TH2F(Form("fHistDEtaDPhiXI%d_Id%d_V%d",i,identifier,j),Form("fHistDEtaDPhiXI%d_Id%d_V%d",i,identifier,j),nbins[2],min[2],max[2],nbins[3],min[3],max[3]);
				fHistDEtaDPhiXI[identifier][i][j]->GetXaxis()->SetTitle(Form("#Delta #varphi^{#gamma-h} %0.1f<#xi<%0.1f",fArray_XI_Bins[i],fArray_XI_Bins[i+1]));
				fHistDEtaDPhiXI[identifier][i][j]->GetYaxis()->SetTitle("#Delta #eta^{#gamma-h}");
				fOutputListXi->Add(fHistDEtaDPhiXI[identifier][i][j]);
			}
			fHistptAssHadronXi[identifier][i] = new TH1F(Form("fHistPtAssHadronXi%d_Id%d",i,identifier),Form("fHistPtAssHadronXi%d_Id%d",i,identifier), nbins[1], min[1], max[1]);
			fHistptAssHadronXi[identifier][i]->GetXaxis()->SetTitle(Form("p_{T}^{assoc. h} %0.1f<#xi<%0.1f",fArray_XI_Bins[i],fArray_XI_Bins[i+1]));
			fHistptAssHadronXi[identifier][i]->GetYaxis()->SetTitle(Form("dN^{assoc. h}/dp_{T}^{h} [counts/%0.1f GeV/c]",fHistptAssHadronXi[identifier][i]->GetBinWidth(1)));
			fOutputListTrAs->Add(fHistptAssHadronXi[identifier][i]);

			fHistptTriggXi[identifier][i] = new TH1F(Form("fHistptTriggXi%d_Id%d",i,identifier),Form("fHistptTriggXi%d_Id%d",i,identifier), nbins[0], min[0], max[0]);
			fHistptTriggXi[identifier][i]->GetXaxis()->SetTitle(Form("p_{T}^{cluster} %0.1f<#xi<%0.1f",fArray_XI_Bins[i],fArray_XI_Bins[i+1]));
			fHistptTriggXi[identifier][i]->GetYaxis()->SetTitle(Form("dN^{assoc. h}/dp_{T}^{h} [counts/%0.1f GeV/c]",fHistptTriggXi[identifier][i]->GetBinWidth(1)));
			fOutputListTrAs->Add(fHistptTriggXi[identifier][i]);
		}*/
	}
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   Test THn Sparse for the 2D histograms
	//   Dimensions are eta,phi
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	Int_t dimThn = 0;
    TString titleThn[20];
    Int_t nbinsThn[20] = {0};
    Double_t minThn[20] = {0.};
    Double_t maxThn[20] = {0.};
    Double_t *binEdgesThn[20] = {0};

    titleThn[dimThn] = "#varphi";
    nbinsThn[dimThn] = 45;
    Double_t phiArray[45+1];
    binEdgesThn[dimThn] = phiArray;
    GenerateFixedBinArray(45,-90.,270.,phiArray);
    minThn[dimThn] = -90.;
    maxThn[dimThn] = 270.;
    dimThn++;

    titleThn[dimThn] = "#eta";
    nbinsThn[dimThn] = 80;
    //binEdgesThn[dim] = 80.;
    Double_t etaArray[80+1];
    binEdgesThn[dimThn] = etaArray;
    GenerateFixedBinArray(80,-2.,2.,etaArray);
    minThn[dimThn] = -2.;
    maxThn[dimThn] = 2.;
    dimThn++;

    titleThn[dimThn] = "E_{T} gamma";
    nbinsThn[dimThn] = kNoGammaBins;
    binEdgesThn[dimThn] = fArray_G_Bins;
    minThn[dimThn] = fArray_G_Bins[0];
    maxThn[dimThn] = fArray_G_Bins[kNoGammaBins];
    dimThn++;

    titleThn[dimThn] = "z_{T}";
    nbinsThn[dimThn] = kNoZtBins;
    binEdgesThn[dimThn] = fArray_ZT_Bins;
    minThn[dimThn] = fArray_ZT_Bins[0];
    maxThn[dimThn] = fArray_ZT_Bins[kNoZtBins];
    dimThn++;

    titleThn[dimThn] = "#Xi";
    nbinsThn[dimThn] = kNoXiBins;
    binEdgesThn[dimThn] = fArray_XI_Bins;
    minThn[dimThn] = fArray_XI_Bins[0];
    maxThn[dimThn] = fArray_XI_Bins[kNoXiBins];
    dimThn++;

    titleThn[dimThn] = "z vertex position";
    nbinsThn[dimThn] = kNvertBins;
    binEdgesThn[dimThn] = fArrayNVertBins;
    minThn[dimThn] = fArrayNVertBins[0];
    maxThn[dimThn] = fArrayNVertBins[kNvertBins];
    dimThn++;

    titleThn[dimThn] = "ME(0) or SE(1)";
    static const Int_t nSeMeBins=2;
 	Double_t SeMeBinArray[nSeMeBins+1]  = {0,1,2};
    nbinsThn[dimThn] = nSeMeBins;
    binEdgesThn[dimThn] = SeMeBinArray;
    minThn[dimThn] = SeMeBinArray[0];
    maxThn[dimThn] = SeMeBinArray[nSeMeBins];
    dimThn++;

    if(fForceBeamType != AliAnalysisTaskEmcal::kpp)
    {
     	static const Int_t nEvtPlaneBins=3;
    		Double_t evtPlaneArray[nEvtPlaneBins+1] = {0,1,2,3};// = {In Plane, MP, Out of Plane};
        static const Int_t nCentHistBins=4;
     	Double_t centBinArray[nCentHistBins+1]  = {0.0,10.0,30.0,60.0,100.0};

     	//..Event plane
     	titleThn[dimThn] = "IP(0), MP(1), OP(2)";
     	nbinsThn[dimThn] = nEvtPlaneBins;
     	binEdgesThn[dimThn] = evtPlaneArray;
        minThn[dimThn] = evtPlaneArray[0];
     	maxThn[dimThn] = evtPlaneArray[nEvtPlaneBins];
     	dimThn++;

     	//..Centrality
     	titleThn[dimThn] = "Centrality %";
     	nbinsThn[dimThn] = nCentHistBins;
     	binEdgesThn[dimThn] = centBinArray;
     	minThn[dimThn] = centBinArray[0];
     	maxThn[dimThn] = centBinArray[nCentHistBins];
     	dimThn++;
    }

    fCorrVsManyThings   = new THnSparseF("CorrVsManyThings", "CorrVsManyThings", dimThn, nbinsThn, minThn, maxThn);
//    fCorrVsManyThingsME = new THnSparseF("CorrVsManyThingsME", "CorrVsManyThingsME", dimThn, nbinsThn, minThn, maxThn);
    for(Int_t i=0;i<dimThn;i++)
    {
		fCorrVsManyThings->GetAxis(i)->SetTitle(titleThn[i]);
		fCorrVsManyThings->SetBinEdges(i, binEdgesThn[i]);
//		fCorrVsManyThingsME->GetAxis(i)->SetTitle(titleThn[i]);
//		fCorrVsManyThingsME->SetBinEdges(i, binEdgesThn[i]);
    }
    //fCorrVsManyThings->Sumw2();
    fOutput->Add(fCorrVsManyThings);
//    fOutput->Add(fCorrVsManyThingsME);

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   THn Sparse for the Cluster properties
	//   Dimensions are Cluster Energy
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Int_t dimThnQA = 0;
    TString titleThnQA[10];
    Int_t nbinsThnQA[10] = {0};
    Double_t minThnQA[10] = {0.};
    Double_t maxThnQA[10] = {0.};
    Double_t *binEdgesThnQA[10] = {0};

    titleThnQA[dimThnQA] = "E_{#gamma}";
    nbinsThnQA[dimThnQA] = 240;
    Double_t EgArray[240+1];
    binEdgesThnQA[dimThnQA] = EgArray;
    GenerateFixedBinArray(240,0,30,EgArray);
    minThnQA[dimThnQA] = 0;
    maxThnQA[dimThnQA] = 30;
    dimThnQA++;

    titleThnQA[dimThnQA] = "#lambda_{0}";
    nbinsThnQA[dimThnQA] = 750;
    Double_t ShapeArray[750+1];
    binEdgesThnQA[dimThnQA] = ShapeArray;
    GenerateFixedBinArray(750,0,4,ShapeArray);
    minThnQA[dimThnQA] = 0;
    maxThnQA[dimThnQA] = 4;
    dimThnQA++;

    titleThnQA[dimThnQA] = "NLM";
    nbinsThnQA[dimThnQA] = 4;
    Double_t NLMArray[4+1];
    binEdgesThnQA[dimThnQA] = NLMArray;
    GenerateFixedBinArray(4,0,4,NLMArray);
    minThnQA[dimThnQA] = 0;
    maxThnQA[dimThnQA] = 4;
    dimThnQA++;

    titleThnQA[dimThnQA] = "#Cells";
    nbinsThnQA[dimThnQA] = 10;
    Double_t nCellsArray[10+1];
    binEdgesThnQA[dimThnQA] = nCellsArray;
    GenerateFixedBinArray(10,0,10,nCellsArray);
    minThnQA[dimThnQA] = 0;
    maxThnQA[dimThnQA] = 10;
    dimThnQA++;

    titleThnQA[dimThnQA] = "cell distance to bad channel";
    nbinsThnQA[dimThnQA] = 10;
    Double_t distanceArray[10+1];
    binEdgesThnQA[dimThnQA] = distanceArray;
    GenerateFixedBinArray(10,0,10,distanceArray);
    minThnQA[dimThnQA] = 0;
    maxThnQA[dimThnQA] = 10;
    dimThnQA++;

    //..ID array -  0 - leading, 1- track matched, 2 - leading & track matched
    titleThnQA[dimThnQA] = "ID code of photon";
    nbinsThnQA[dimThnQA] = 5;
    Double_t IdArray[5+1];
    binEdgesThnQA[dimThnQA] = IdArray;
    GenerateFixedBinArray(5,0,5,IdArray);
    minThnQA[dimThnQA] = 0;
    maxThnQA[dimThnQA] = 5;
    dimThnQA++;

    /*
         titleThnQA[dimThnQA] = "E/p";
        nbinsThnQA[dimThnQA] = 100;
        Double_t EPArray[100+1];
        binEdgesThnQA[dimThnQA] = EPArray;
        GenerateFixedBinArray(100,0,2,EPArray);
        minThnQA[dimThnQA] = 0;
        maxThnQA[dimThnQA] = 2;
        dimThnQA++;

        titleThnQA[dimThnQA] = "M_{#gamma#gamma x}";
        nbinsThnQA[dimThnQA] = 1000;
        Double_t MArray[1000+1];
        binEdgesThnQA[dimThnQA] = MArray;
        GenerateFixedBinArray(1000,0,10,MArray);
        minThnQA[dimThnQA] = 0;
        maxThnQA[dimThnQA] = 10;
        dimThnQA++;
    */
    //..additional things to put inside: time , number of Cells, hit position eta-phi
    if(fPlotQA==1)
    {
    		fClusterProp= new THnSparseF("ClusterProp", "ClusterProp", dimThnQA, nbinsThnQA, minThnQA, maxThnQA);
    		for(Int_t i=0;i<dimThnQA;i++)
    		{
    			fClusterProp->GetAxis(i)->SetTitle(titleThnQA[i]);
    			fClusterProp->SetBinEdges(i, binEdgesThnQA[i]);
    		}
    		fOutput->Add(fClusterProp);
    }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //   Tyler's Special Histograms
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fHistClusPairInvarMasspT= new TH2F("fHistClusPairInvarMasspT","fHistClusPairInvarMasspT", 3000, 0, 20.0, 250, 0, 50);
    fHistClusPairInvarMasspT->GetXaxis()->SetTitle("M_{#gamma#gamma}");
	fHistClusPairInvarMasspT->GetYaxis()->SetTitle("p_{T}_{#pi^{0}}");
	fOutput->Add(fHistClusPairInvarMasspT);

	fMAngle= new TH2F("fMAngle","fMAngle", 2000, 0, 3, 500, 0, 3.14159);
	fMAngle->GetXaxis()->SetTitle("M_{#gamma#gamma}");
	fMAngle->GetYaxis()->SetTitle("Opening Angle [rad]");
	fOutput->Add(fMAngle);

	fPtAngle= new TH2F("fPtAngle","fPtAngle", 250, 0, 50, 500, 0, 3.14159);
	fPtAngle->GetXaxis()->SetTitle("p_{T}_{#pi^{0}}");
	fPtAngle->GetYaxis()->SetTitle("Opening Angle [rad]");
	fOutput->Add(fPtAngle);

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//  Eliane's Special Histograms
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	fHistEvsPt = new TH2F("fHistEvsPt","fHistEvsPt", nbins[0], min[0], max[0], 250, 0, 50);
	fHistEvsPt->GetXaxis()->SetTitle("p_{T}_{#gamma}");
	fHistEvsPt->GetYaxis()->SetTitle("E_{#gamma}");
	fOutput->Add(fHistEvsPt);

	//test!!
	fHistPi0 = new TH1F(Form("fHistPi0_%0d",1),Form("fHistPi0_%0d",1), 500, 0, 0.5);
	fHistPi0->GetXaxis()->SetTitle("M_{#gamma#gamma}");
	fHistPi0->GetYaxis()->SetTitle("Entries");
	fOutput->Add(fHistPi0);

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   Special QA histograms (also to get more info what is going on in mixed event for trigger data)
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	fHistMatchEtaPhiAllCl2 = new TH2F("fHistMatchEtaPhiAllCl2", "fHistMatchEtaPhiAllCl2;#Delta#eta;#Delta#phi", 100, -0.1, 0.1, 100, -0.1, 0.1);
	fOutput->Add(fHistMatchEtaPhiAllCl2);
	fHistMatchEtaPhiAllCl3 = new TH2F("fHistMatchEtaPhiAllCl3", "fHistMatchEtaPhiAllCl3;#Delta#eta;#Delta#phi", 100, -0.1, 0.1, 100, -0.1, 0.1);
	fOutput->Add(fHistMatchEtaPhiAllCl3);
	fHistMatchEtaPhiAllCl4 = new TH2F("fHistMatchEtaPhiAllCl4", "fHistMatchEtaPhiAllCl4;#Delta#eta;#Delta#phi", 100, -0.1, 0.1, 100, -0.1, 0.1);
	fOutput->Add(fHistMatchEtaPhiAllCl4);

	for(Int_t identifier=0;identifier<kNIdentifier+3;identifier++)
	{
		//..geometrical hit distribution of clusters
		fHistDEtaDPhiGammaQA[identifier] = new TH2F(Form("fHistDEtaDPhiGammaQA%d_Id%d",0,identifier),Form("fHistDEtaDPhiGammaQA%d_Id%d",0,identifier),142,-0.71,0.71,311,75.9,331.1);
		fHistDEtaDPhiGammaQA[identifier]->GetXaxis()->SetTitle("#eta^{#gamma}");
		fHistDEtaDPhiGammaQA[identifier]->GetYaxis()->SetTitle("#varphi^{#gamma}");
		fOutputListQA->Add(fHistDEtaDPhiGammaQA[identifier]);

		//..geometrical hit distribution of tracks
		fHistDEtaDPhiTrackQA[identifier] = new TH2F(Form("fHistDEtaDPhiTrackQA%d_Id%d",0,identifier),Form("fHistDEtaDPhiTrackQA%d_Id%d",0,identifier),182,-0.91,0.91,360,0,360);
		fHistDEtaDPhiTrackQA[identifier]->GetXaxis()->SetTitle("#eta^{hadron}");
		fHistDEtaDPhiTrackQA[identifier]->GetYaxis()->SetTitle("#varphi^{hadron}");
		fOutputListQA->Add(fHistDEtaDPhiTrackQA[identifier]);

	    //..Cluster Info
		fHistCellsCluster[identifier] = new TH2F(Form("fHistCellsCluster%d_Id%d",0,identifier),Form("fHistCellsCluster%d_Id%d",0,identifier),240,0,30,50,0,50);
		fHistCellsCluster[identifier]->GetXaxis()->SetTitle("E^{cluster}");
		fHistCellsCluster[identifier]->GetYaxis()->SetTitle("N_{cells}");
		fOutputListQA->Add(fHistCellsCluster[identifier]);

		//..Time information
		fHistClusterTime[identifier] = new TH2F(Form("fHistClusterTime%d_Id%d",0,identifier),Form("fHistClusterTime%d_Id%d",0,identifier),20000,-100,100,200,0,40);
		fHistClusterTime[identifier]->GetXaxis()->SetTitle("time [ns]");
		fHistClusterTime[identifier]->GetYaxis()->SetTitle("pT");
		fOutputListQA->Add(fHistClusterTime[identifier]);

		//..
	}

	//..The END
	fOutput->Add(fOutputList1);
	fOutput->Add(fOutputListTrAs);
	fOutput->Add(fOutputListGamma);
	fOutput->Add(fOutputListZeta);
	fOutput->Add(fOutputListXi);
	fOutput->Add(fOutputListQA);

	PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}
//________________________________________________________________________
void AliAnalysisTaskGammaHadron::InitEventMixer()
{
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::InitEventMixer()"<<endl;
	//--The effective pool size in events is set by trackDepth, so more
	//--low-mult events are required to maintain the threshold than
	//--high-mult events. Centrality pools are indep. of data histogram
	//--binning, no need to match.

	//..Centrality pools
	Int_t nCentBins=fMixBCent->GetNbins();
	Double_t centBins[nCentBins+1];
	centBins[0] = fMixBCent->GetBinLowEdge(1);
	for(Int_t i=1; i<=nCentBins; i++)
	{
		centBins[i] = fMixBCent->GetBinUpEdge(i);
	}

	//..Z-vertex pools
	Int_t nZvtxBins=fMixBZvtx->GetNbins();
	Double_t zvtxbin[nZvtxBins+1];
	zvtxbin[0] = fMixBZvtx->GetBinLowEdge(1);
	for(Int_t i=1; i<=nZvtxBins; i++)
	{
		zvtxbin[i] = fMixBZvtx->GetBinUpEdge(i);
	}

	//..in case no external pool is provided create one here
	if(!fPoolMgr)
	{
		cout<<"....  Pool Manager Created ...."<<endl;
		fPoolMgr = new AliEventPoolManager(fPoolSize, fTrackDepth, nCentBins, centBins, nZvtxBins, zvtxbin);
		fPoolMgr->SetTargetValues(fTrackDepth, 0.1, 5);  //pool is ready at 0.1*fTrackDepth = 5000 or events =5
		//save this pool by default
	}
	else
	{
		//..lock all pools
		//..clears empty pools and sets them locked
		//..(is only possible because all save flags are ture in my case  - NASTY NASTY)
		fPoolMgr->ClearPools();
		cout<<"....  Pool Manager Provided From File ...."<<endl;
	}

	//..Check binning of pool manager (basic dimensional check for the time being) to see whether external pool fits the here desired one??
	if( (fPoolMgr->GetNumberOfMultBins() != nCentBins) || (fPoolMgr->GetNumberOfZVtxBins() != nZvtxBins) )
	{
		AliFatal("Binning of given pool manager not compatible with binning of correlation task!");
	}
	//if you want to save the pool:
	// If some bins of the pool should be saved, fEventPoolOutputList must be given
	// using AddEventPoolToOutput() (to increase size of fEventPoolOutputList)
	// Note that this is in principle also possible, if an external poolmanager was given
//	if(fEventPoolOutputList.size())

	if(fSavePool==1)
	{
		//? is this an option to save only specific pools instead of the full pool manager?
		//for(Int_t i = 0; i < fEventPoolOutputList.size(); i++)
		{
			/*Double_t minCent = fEventPoolOutputList[i][0];
			Double_t maxCent = fEventPoolOutputList[i][1];
			Double_t minZvtx = fEventPoolOutputList[i][2];
			Double_t maxZvtx = fEventPoolOutputList[i][3];
			Double_t minPt   = fEventPoolOutputList[i][4];
			Double_t maxPt   = fEventPoolOutputList[i][5];
            */
		    //If the pool fulfills the given criteria the saveflag is set to true
			//the flag is used in the ClearPools function to not delete the pool content
			fPoolMgr->SetSaveFlag(-1, 10000, -10000, 100000, 0, 0, -1, 10000000);

			/*
			In case you don't want to store all pools but only the ones specified above
			you have to rund at the very end of filling these lines.

			// Clear unnecessary pools before saving and locks the pool
			fPoolMgr->ClearPools();
			*/
		}
		fOutput->Add(fPoolMgr);
	}

	//..Basic checks and printing of pool properties
	fPoolMgr->Validate();
}
///
/// Saves event pool to be used later
/// This might be a possibility to increase statistic
/// in the mixed event pool
//____________________________________________________________________
void AliAnalysisTaskGammaHadron::AddEventPoolsToOutput(Double_t minCent, Double_t maxCent,  Double_t minZvtx, Double_t maxZvtx, Double_t minPt, Double_t maxPt)
{
	//..This allows you to add only specific pools and not the full pool manager to the output file
	std::vector<Double_t> binVec;
	binVec.push_back(minCent);
	binVec.push_back(maxCent);
	binVec.push_back(minZvtx);
	binVec.push_back(maxZvtx);
	binVec.push_back(minPt);
	binVec.push_back(maxPt);
	fEventPoolOutputList.push_back(binVec);
}

//________________________________________________________________________
void AliAnalysisTaskGammaHadron::ExecOnce()
{
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::ExecOnce()"<<endl;

	//..This function does...
	AliAnalysisTaskEmcal::ExecOnce();

}
///Overwrites the AliAnalysisTaskEmcal::IsEventSelected()
///function
///
//________________________________________________________________________
Bool_t AliAnalysisTaskGammaHadron::IsEventSelected()
{
	//..checks: rejects DAQ incompletes, magnetic field selection, offline trigger
	//..vertex selection, SPD pile-up (if enabled), centrality cuts
	if (!fEventCuts.AcceptEvent(InputEvent()))
	{
		PostData(1, fOutput);
		return kFALSE;
	}

	//.. .. .. .. .. .. .. .. .. .. ..
	//..Start of copy part from AliAnalysisTaskEmcal
	if (!fTrigClass.IsNull())
	{
		TString fired;

		const AliAODEvent *aev = dynamic_cast<const AliAODEvent*>(InputEvent());
		if (aev)
		{
			fired = aev->GetFiredTriggerClasses();
		}
		else cout<<"Error analysis only for AODs"<<endl;

		if (!fired.Contains("-B-"))
		{
			if (fGeneralHistograms) fHistEventRejection->Fill("trigger",1);
			return kFALSE;
		}
     /* // For some wired reason gives an error in alibild when doing pull request
        // need to check that later again because it's an exact copy from "AliAnalysisTaskEmcal"
		std::unique_ptr<TObjArray> arr(fTrigClass.Tokenize("|"));
		if (!arr)
		{
			if (fGeneralHistograms) fHistEventRejection->Fill("trigger",1);
			return kFALSE;
		}
		Bool_t match = 0;
		for (Int_t i=0;i<arr->GetEntriesFast();++i)
		{
			TObject *obj = arr->At(i);
			if (!obj)
				continue;

			//Check if requested trigger was fired
			TString objStr = obj->GetName();
			if(fEMCalTriggerMode == kOverlapWithLowThreshold &&
					(objStr.Contains("J1") || objStr.Contains("J2") || objStr.Contains("G1") || objStr.Contains("G2"))) {
				// This is relevant for EMCal triggers with 2 thresholds
				// If the kOverlapWithLowThreshold was requested than the overlap between the two triggers goes with the lower threshold trigger
				TString trigType1 = "J1";
				TString trigType2 = "J2";
				if(objStr.Contains("G"))
				{
					trigType1 = "G1";
					trigType2 = "G2";
				}
				if(objStr.Contains(trigType2) && fired.Contains(trigType2.Data()))
				{ //requesting low threshold + overlap
					match = 1;
					break;
				}
				else if(objStr.Contains(trigType1) && fired.Contains(trigType1.Data()) && !fired.Contains(trigType2.Data())) { //high threshold only
					match = 1;
					break;
				}
			}
			else
			{
				// If this is not an EMCal trigger, or no particular treatment of EMCal triggers was requested,
				// simply check that the trigger was fired
				if (fired.Contains(obj->GetName()))
				{
					match = 1;
					break;
				}
			}
		}
		if (!match)
		{
			if (fGeneralHistograms) fHistEventRejection->Fill("trigger",1);
			return kFALSE;
		}*/
	}

	if (fTriggerTypeSel != kND)
	{
		if (!HasTriggerType(fTriggerTypeSel))
		{
			if (fGeneralHistograms) fHistEventRejection->Fill("trigTypeSel",1);
			return kFALSE;
		}
	}
    /*
	//.. Maybe these two as well .. .. ..
	if (fSelectPtHardBin != -999 && fSelectPtHardBin != fPtHardBin)  {
		if (fGeneralHistograms) fHistEventRejection->Fill("SelPtHardBin",1);
		return kFALSE;
	}

	// Reject filter for MC data
	if (!CheckMCOutliers()) return kFALSE;
    */
	//..End of copy part from AliAnalysisTaskEmcal
	//.. .. .. .. .. .. .. ..
	return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskGammaHadron::Run()
{
	//..This function is called in AliAnalysisTaskEmcal::UserExec.
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::Run()"<<endl;
	//..Determine the trigger for the current event
	fCurrentEventTrigger = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

	//..Get the ClusterContainer for the event
	//AliClusterContainer* clusters  = GetClusterContainer(0);  //how do I know which cells are selected
	if (!fCaloClusters)
	{
		fCaloClusters = (TClonesArray*)GetClusterContainer(0);
		cout<<"load calo clusters"<<endl;
	}
	//..Get the emcal cells
	if (!fCaloCells)
	{
		if (fCaloCellsName.IsNull())
		{
			fCaloCells = InputEvent()->GetEMCALCells();
		}
		else
		{
			fCaloCells =  dynamic_cast<AliVCaloCells*>(InputEvent()->FindListObject(fCaloCellsName));
			if (!fCaloCells) AliError(Form("%s: Could not retrieve cells %s!", GetName(), fCaloCellsName.Data()));
		}
		cout<<"load calo cells"<<endl;
	}

	if (!fGeom)
	{
		AliWarning(Form("%s - AliAnalysisTaskGammaHadron::Run - Geometry is not available!", GetName()));
		return kFALSE;
	}
	//..check here some basic properties of the event
	//is centrality and zvertex in range? - see if this is actually done in IsSelected in the EMCal Task

	/*if(fCurrentEventTrigger & fTriggerType)
	{
		if(fCurrentEventTrigger & fMixingEventType)cout<<"*********************contains both triggers!!"<<endl;
	}*/

	//..for same event only analyse events when there is a cluster inside
	//..and when the event has the correct trigger
	if (fDoMixing==0 && !fCaloClusters)                         return kFALSE;
	if (fDoMixing==0 && !(fCurrentEventTrigger & fTriggerType)) return kFALSE;

	return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskGammaHadron::FillHistograms()
{
	//..This function is called in AliAnalysisTaskEmcal::UserExec.
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::FillHistograms()"<<endl;

	// 1. First get an event pool corresponding in mult (cent) and
	//    zvertex to the current event. Once initialized, the pool
	//    should contain nMix (reduced) events. This routine does not
	//    pre-scan the chain. The first several events of every chain
	//    will be skipped until the needed pools are filled to the
	//    specified depth. If the pool categories are not too rare, this
	//    should not be a problem. If they are rare, you could lose
	//    statistics.

	// 2. Collect the whole pool's content of tracks into one TObjArray
	//    (bgTracks), which is effectively a single background super-event.

	// 3. The reduced and bgTracks arrays must both be passed into
	//    FillCorrelations(). Also nMix should be passed in, so a weight
	//    of 1./nMix can be applied.

	//..Get pool containing tracks from other events like this one
	Double_t zVertex = fVertex[2];
	AliParticleContainer* tracks =0x0;
	tracks   = GetParticleContainer(0);

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//    Mixed event section
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	if(fDoMixing==1)
	{
		AliEventPool* pool = 0x0;
		pool = fPoolMgr->GetEventPool(fCent, zVertex);
		if (!pool)
		{
			AliWarning(Form("No pool found. Centrality %f, ZVertex %f",fCent, zVertex));
			return kFALSE;
		}

        //. . . . . . . . . . . . .
		//..Start combining triggers (from fTriggerType events)
		//..with a pool filled with tracks (from fMixingEventType)
		if(pool->IsReady() && (fCurrentEventTrigger & fTriggerType))
		{
			//..get number of current events in pool
			Int_t nMix = pool->GetCurrentNEvents();

			//cout<<"number of events in pool: "<<nMix<<endl;
			for(Int_t jMix=0; jMix<nMix; jMix++)
			{
				TObjArray* bgTracks=0x0;
				bgTracks = pool->GetEvent(jMix);

				if(!bgTracks)
				{
					cout<<"could not retrieve TObjArray from EventPool!"<<endl;
				}
				//..Loop over clusters and fill histograms
				if(fGammaOrPi0==0) CorrelateClusterAndTrack(0,bgTracks,0,1.0/nMix);//correlate with mixed event
				else               CorrelatePi0AndTrack(0,bgTracks,0,1.0/nMix);    //correlate with mixed event
			}
		}
		// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		//    Update the pool
		// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
		//..update pool only with tracks from event type fMixingEventType,
		//..and do NOT add tracks from GA triggered events (BIT(15))
		if((fCurrentEventTrigger & fMixingEventType) && ((fCurrentEventTrigger & AliVEvent::kEMCEGA)==0))
		{
			TObjArray* tracksClone=0x0;
			tracksClone = CloneToCreateTObjArray(tracks);

			//..if there is no track object or the pool is locked do not update
			if(tracksClone && !pool->GetLockFlag())
			{
				pool->UpdatePool(tracksClone);
			}
		}
	}
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//    Same event section
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//..Loop over clusters and fill histograms
	//..Do this only for events that are of fTriggerType
	if(fCurrentEventTrigger & fTriggerType)
	{
		if(fGammaOrPi0==0) CorrelateClusterAndTrack(tracks,0,1,1);//correlate with same event
		else               CorrelatePi0AndTrack(tracks,0,1,1);    //correlate with same event
	}

	return kTRUE;
}
///
/// Clone the tracks to create an object reduced
/// in size to be filled in the event pool
///
//________________________________________________________________________
TObjArray* AliAnalysisTaskGammaHadron::CloneToCreateTObjArray(AliParticleContainer* tracks)
{
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::CloneToCreateTObjArray()"<<endl;
	//..clones a track list
	if(!tracks)                            return 0;
	if(tracks->GetNAcceptedParticles()==0) return 0;

	TObjArray* tracksClone = new TObjArray;
	tracksClone->SetOwner(kTRUE);

	Int_t NoOfTracksInEvent =tracks->GetNParticles();
	AliVParticle* track=0;

	for(Int_t NoTrack = 0; NoTrack < NoOfTracksInEvent; NoTrack++)
	{
		track = (AliVParticle*)tracks->GetAcceptParticle(NoTrack);
		if(!track)continue; //check if the track is a good track
		//tracksClone->Add((AliVParticle*)track);  //only add accepted tracks
		tracksClone->Add(new AliPicoTrack(track->Pt(), track->Eta(), track->Phi(), track->Charge(), 0, 0, 0, 0));
	}
	if(tracksClone->GetEntries()!=tracks->GetNAcceptedParticles())cout<<"!!!!!!! Major error!!!! "<<"Accepted tracks in event: "<<tracks->GetNAcceptedParticles()<<", Tracks in TObjArray: "<<tracksClone->GetEntries()<<endl;

	return tracksClone;
}
///
/// Select tracks and clusters to correlate with each other
///
//________________________________________________________________________
Int_t AliAnalysisTaskGammaHadron::CorrelateClusterAndTrack(AliParticleContainer* tracks,TObjArray* bgTracksArray,Bool_t SameMix, Double_t InputWeight)
{
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::CorrelateClusterAndTrack()"<<endl;

	//...........................................
	//..Do cluster loop.
	AliClusterContainer* clusters  = GetClusterContainer(0);  //how do I know which cells are selected
	if (!clusters) return 0;
	Int_t NoOfClustersInEvent =clusters->GetNClusters();
	Double_t EffWeight_Gamma;
	Double_t EffWeight_Hadron;
	Double_t Weight=1;    //weight to normalize mixed and same event distributions individually

	AliVCluster* cluster = 0;
	AliVCluster* clusterT = 0;
	AliVCluster* cluster2= 0;
	AliVParticle* trackNULL=0;

	//...........................................
	//..for mixed events normalize per events in pool
	if(SameMix==0)
	{
		Weight=InputWeight;
	}
	//..for same events normalize later by counting the entries in spec. norm histograms
	if(SameMix==1)
	{
		Weight=1;
	}

	//...........................................
	//..run the loop for filling the histograms
	Int_t GammaCounter=0;
	for(Int_t NoCluster1 = 0; NoCluster1 < NoOfClustersInEvent; NoCluster1++ )
	{
		cluster=(AliVCluster*) clusters->GetAcceptCluster(NoCluster1); //->GetCluster(NoCluster1);
		//clusterT=(AliVCluster*) fCaloClusters->GetAcceptCluster(NoCluster1);
		if(!cluster)continue; //check if the cluster is a good cluster
		//if(!cluster)continue; //check if the cluster is a good cluster
		//clusters->GetLeadingCluster("e");
        //if(cluster->E()!=clusterT->E())cout<<"clusters are different!!??"<endl;

		TLorentzVector CaloClusterVec;
		clusters->GetMomentum(CaloClusterVec, cluster);
		AliTLorentzVector aliCaloClusterVec = AliTLorentzVector(CaloClusterVec); //..can acess phi from
		EffWeight_Gamma=GetEff(aliCaloClusterVec);

		//------------------------------------------------
        //..This section is for the moment to test
		//..cluster distributions without cuts
		if(SameMix==1)FillQAHisograms(0,clusters,cluster,trackNULL);

		fFiducialCellCut->SetNumberOfCellsFromEMCALBorder(0);
		if(!AccClusterForAna(clusters,cluster))continue; //check if the cluster is a good cluster
		if(SameMix==1)FillQAHisograms(1,clusters,cluster,trackNULL);

		fFiducialCellCut->SetNumberOfCellsFromEMCALBorder(1);
		if(!AccClusterForAna(clusters,cluster))continue; //check if the cluster is a good cluster
		if(SameMix==1)FillQAHisograms(2,clusters,cluster,trackNULL);

		fFiducialCellCut->SetNumberOfCellsFromEMCALBorder(2);
		if(!AccClusterForAna(clusters,cluster))continue; //check if the cluster is a good cluster
		if(SameMix==1)FillQAHisograms(3,clusters,cluster,trackNULL);
		//------------------------------------------------

		//...........................................
		//..combine gammas with same event tracks
		GammaCounter++;
		if(SameMix==1)
		{
			fHistNoClusPt->Fill(CaloClusterVec.Pt()); //the .pt only works for gammas (E=M) for other particle this is wrong
			fHistEvsPt   ->Fill(CaloClusterVec.Pt(),CaloClusterVec.E());

			if(!tracks)  return 0;
			Int_t NoOfTracksInEvent =tracks->GetNParticles();
			AliVParticle* track=0;

			Int_t trackCounter=0;
			for(Int_t NoTrack = 0; NoTrack < NoOfTracksInEvent; NoTrack++)
			{
				track = (AliVParticle*)tracks->GetAcceptParticle(NoTrack);
				if(!track)continue; //check if the track is a good track
				trackCounter++;
				//cout<<"..Track number: "<<NoTrack<<", pT = "<<track->Pt()<<endl;

				//..Fill this histogram only for clusters above 5 GeV
				//EffWeight_Hadron=GetEff(<TLorentzVector>track);
				FillGhHisograms(1,aliCaloClusterVec,track,5,Weight);


				if(GammaCounter==1)FillQAHisograms(4,clusters,cluster,track); //fill only once per track (first gamma) - good for each track
				if(trackCounter==1)FillQAHisograms(5,clusters,cluster,track); //fill only once per gamma (first track) - good for gamma distr.
			}
			//...........................................
			//..double cluster loop for testing an anti pi0 cut
			for( Int_t NoCluster2 = 0; NoCluster2 < NoOfClustersInEvent; NoCluster2++ )
			{
				if(NoCluster1!=NoCluster2 && NoCluster1<NoCluster2) //..don't combine same clusters and don't combine them twice
				{
					cluster2=(AliVCluster*) clusters->GetAcceptCluster(NoCluster2);
					if(!cluster2 || !AccClusterForAna(clusters,cluster2))continue; //check if the cluster is a good cluster

					TLorentzVector CaloClusterVec2;
					TLorentzVector CaloClusterVecpi0;
					clusters->GetMomentum(CaloClusterVec2, cluster2);
					if(cluster2->GetNonLinCorrEnergy()>2 && cluster->GetNonLinCorrEnergy()>2)
					{
						CaloClusterVecpi0=CaloClusterVec+CaloClusterVec2;
						fHistPi0->Fill(CaloClusterVecpi0.M());
						fHistClusPairInvarMasspT->Fill(CaloClusterVecpi0.M(),CaloClusterVecpi0.Pt());
					}
				}
			}
		}
		//...........................................
		//..combine gammas with mixed event tracks
		if(SameMix==0)
		{
			Int_t Nbgtrks = bgTracksArray->GetEntries();
			for(Int_t ibg=0; ibg<Nbgtrks; ibg++)
			{
				AliPicoTrack* track = static_cast<AliPicoTrack*>(bgTracksArray->At(ibg));
				if(!track) continue;

				//EffWeight_Hadron=GetEff(<TLorentzVector>track);
				FillGhHisograms(0,aliCaloClusterVec,track,5,Weight);
			}
		}
	}
	return GammaCounter;
}
///
/// Select tracks and pi0/decay-gammas to correlate with each other
///
//________________________________________________________________________
Int_t AliAnalysisTaskGammaHadron::CorrelatePi0AndTrack(AliParticleContainer* tracks,TObjArray* bgTracksArray,Bool_t SameMix, Double_t InputWeight)
{
  if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::CorrelatePi0AndTrack()"<<endl;

	//...........................................
	//--Do cluster loop.
	AliClusterContainer* clusters  = GetClusterContainer(0);
	if (!clusters) return 0;
	Int_t NoOfClustersInEvent =clusters->GetNClusters();
	Int_t nAccClusters = 0;
	Int_t nAccPi0Clusters = 0;
	Double_t Pi0Mass = 0.13487; // Hi Michael -> this center value should also be made flexible
	Double_t Pi0Window = 0.02;  //0.03 // Hi Michael -> the width will vary with pT
	Double_t EffWeight_Gamma;
	Double_t EffWeight_Hadron;
	Double_t Weight;    //weight to normalize mixed and same event distributions individually

	AliVCluster* cluster = 0;
	AliVCluster* cluster2= 0;
	AliVParticle* trackNULL=0;

	//...........................................
	//do a small loop to count the triggers in this event
    //** we don't need this loop here any longer because we will
	//**normalize later in the analysis not now so we don't need the
	//** total pi0 count beforehand!!!!!
	if(SameMix==1)
	{
		for(Int_t NoCluster1 = 0; NoCluster1 < NoOfClustersInEvent; NoCluster1++ )
		{
			cluster=(AliVCluster*) clusters->GetAcceptCluster(NoCluster1); //->GetCluster(NoCluster1);
			if(!cluster || !AccClusterForAna(clusters,cluster))continue; //check if the cluster is a good cluster
			//clusters->GetLeadingCluster("e");

			TLorentzVector CaloClusterVec;
			clusters->GetMomentum(CaloClusterVec, cluster);
			//acc if pi0 candidate
			nAccClusters++;

			for(Int_t NoCluster2 = 0; NoCluster2 < NoOfClustersInEvent; NoCluster2++ )
			{
				if(NoCluster1!=NoCluster2)
				{
					cluster2=(AliVCluster*) clusters->GetAcceptCluster(NoCluster2);
					if(!cluster2 || !AccClusterForAna(clusters,cluster2))continue; //check if the cluster is a good cluster

					TLorentzVector CaloClusterVec2;
					TLorentzVector CaloClusterVecpi0;

					//old framework				cluster2->GetMomentum(CaloClusterVec2, fVertex);
					clusters->GetMomentum(CaloClusterVec2, cluster2);
					if(cluster2->E()>2 && cluster->E()>2)
					{
						CaloClusterVecpi0=CaloClusterVec+CaloClusterVec2;
						fHistPi0->Fill(CaloClusterVecpi0.M());
						fHistClusPairInvarMasspT->Fill(CaloClusterVecpi0.M(),CaloClusterVecpi0.Pt(),0.5);  //eventually divide by 2
						fMAngle->Fill(CaloClusterVecpi0.M(), CaloClusterVec.Angle(CaloClusterVec2.Vect()),0.5);
						fPtAngle->Fill(CaloClusterVecpi0.Pt(), CaloClusterVec.Angle(CaloClusterVec2.Vect()),0.5);
						if((CaloClusterVecpi0.M()>=Pi0Mass-Pi0Window) && (CaloClusterVecpi0.M()<=Pi0Mass+Pi0Window)){
							nAccPi0Clusters++; //need eventually to divide by 2, otherwise double counting the pi0's
						}
					}
				}
			}
		}
	}

	//...........................................
	//..for mixed events normalize per events in pool
	if(SameMix==0)
	{
		Weight=InputWeight;
	}
	//..for same events normalize later by counting the entries in spec. norm histograms
	if(SameMix==1)
	{
		Weight=1;
	}
	//...........................................
	//run the loop for filling the histograms
	for( Int_t NoCluster1 = 0; NoCluster1 < NoOfClustersInEvent; NoCluster1++ )
	{
		cluster=(AliVCluster*) clusters->GetAcceptCluster(NoCluster1); //->GetCluster(NoCluster1);
		if(!cluster || !AccClusterForAna(clusters,cluster))continue; //check if the cluster is a good cluster
		//clusters->GetLeadingCluster("e");

		TLorentzVector CaloClusterVec;
		clusters->GetMomentum(CaloClusterVec,cluster);
		AliTLorentzVector aliCaloClusterVec = AliTLorentzVector(CaloClusterVec); //..can acess phi from

		FillQAHisograms(0,clusters,cluster,trackNULL);

		for( Int_t NoCluster2 = 0; NoCluster2 < NoOfClustersInEvent; NoCluster2++ )
		{
			if(NoCluster1!=NoCluster2)
			{
				cluster2=(AliVCluster*) clusters->GetAcceptCluster(NoCluster2);
				if(!cluster2 || !AccClusterForAna(clusters,cluster2))continue; //check if the cluster is a good cluster

				TLorentzVector CaloClusterVec2;
				TLorentzVector aliCaloClusterVecpi0;
				//old framework			cluster2->GetMomentum(CaloClusterVec2, fVertex);
				clusters->GetMomentum(CaloClusterVec2,cluster2); /// Vec+=2 2.1.17
				AliTLorentzVector aliCaloClusterVec2 = AliTLorentzVector(CaloClusterVec2); //..can acess phi from

				if(cluster2->E()>2 && cluster->E()>2)
				{
					aliCaloClusterVecpi0=aliCaloClusterVec+aliCaloClusterVec2;

					if((aliCaloClusterVecpi0.M()<Pi0Mass-Pi0Window) || (aliCaloClusterVecpi0.M()>Pi0Mass+Pi0Window)) continue; /// 2.1.17

					//here I don't really know what to do in your case
					//eff of pi0? or eff of gamma? or some mix up of the two ?
					EffWeight_Gamma=GetEff(aliCaloClusterVecpi0);//currently just assigns 1!!! need eventually to input Pi0 efficiency histogram

					fHistNoClusPt->Fill(aliCaloClusterVecpi0.Pt()); //the .pt only works for gammas (E=M) for other particle this is wrong

					//...........................................
					//..combine gammas with same event tracks
					if(SameMix==1)
					{
						//cout<<"SameMix==1"<<endl;
						if(!tracks)  return 0;
						Int_t NoOfTracksInEvent =tracks->GetNParticles();
						AliVParticle* track=0;

						for(Int_t NoTrack = 0; NoTrack < NoOfTracksInEvent; NoTrack++)
						{
							track = (AliVParticle*)tracks->GetAcceptParticle(NoTrack);
							if(!track)continue; //check if the track is a good track

							//..fill here eventually a pi0 four-vector instead of CaloClusterVec
							//EffWeight_Hadron=GetEff(TLorentzVector)track);
							FillGhHisograms(1,aliCaloClusterVecpi0,track,5,Weight);
						}
					}
					//...........................................
					//..combine gammas with mixed event tracks
					if(SameMix==0)
					{
						Int_t Nbgtrks = bgTracksArray->GetEntries();
						for(Int_t ibg=0; ibg<Nbgtrks; ibg++)
						{
							AliPicoTrack* track = static_cast<AliPicoTrack*>(bgTracksArray->At(ibg));
							if(!track) continue;

							//**fill here eventually a pi0 four-vector instead of CaloClusterVec
							//EffWeight_Hadron=GetEff((TLorentzVector)track);
							FillGhHisograms(0,aliCaloClusterVecpi0,track,5,Weight);
						}
					}
				}
			}
		}
	}
	return nAccPi0Clusters/2;
}
///
/// Fill histograms with cluster and track information
///
//________________________________________________________________________
void AliAnalysisTaskGammaHadron::FillGhHisograms(Int_t identifier,AliTLorentzVector ClusterVec,AliVParticle* TrackVec, Double_t ClusterEcut, Double_t Weight)
{
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::FillGhHisograms()"<<endl;

	//..This function fills several histograms under different cut conditions.
	//..it is run within a cluster{ track{}} loop to get all combinations.

	//..A word to the weight - for mixed events it devides by the number of events in the current pool 1/nEvents
	//..                     - for same events you devide by the number of triggers (but not now ->later in the post analysis)
	//..                     - for both you have to take into account the efficiency of your correlated pair
	Double_t deltaEta   = ClusterVec.Eta()-TrackVec->Eta();
	Double_t deltaPhi   = DeltaPhi(ClusterVec,TrackVec);
	Double_t G_PT_Value = ClusterVec.Pt();
	//Double_t ZT_Value   = TMath::Cos(deltaPhi*1/fRtoD)*TrackVec->P()/ClusterVec.P(); //   TrackVec->Pt()/G_PT_Value;
	Double_t ZT_Value   = TrackVec->Pt()/G_PT_Value; //   TrackVec->Pt()/G_PT_Value;
	//..Careful here: usually this is done for an opening angle (hadron-jet axis) of less than 90�. Due to
	//..resolution momentum smearing (our guess - check that!) there are particles appearing at angles greater than 90�
	Double_t XI_Value=-50;
	if(ZT_Value>0)
	{
		XI_Value   = TMath::Log(1.0/ZT_Value);
	}
	Double_t zVertex = fVertex[2];
	//..all from EMcal base class : fEPV0,fEPV0A,fEPV0C
    //fEPV0  = aliEP->GetEventplane("V0" ,InputEvent());
    //fEPV0A = aliEP->GetEventplane("V0A",InputEvent());
    //fEPV0C = aliEP->GetEventplane("V0C",InputEvent());
	Double_t evtPlaneAngle= DeltaPhi(ClusterVec,fEPV0);
	Int_t evtPlaneCategory=-1;
	if(evtPlaneAngle>-60 && evtPlaneAngle<=60)       evtPlaneCategory=0;
	else if ((evtPlaneAngle>60 && evtPlaneAngle<=120) || evtPlaneAngle<=-60 || evtPlaneAngle>240)evtPlaneCategory=1;
	else if (evtPlaneAngle>120 && evtPlaneAngle<=240)evtPlaneCategory=2;

	Double_t valueArray[9];
	valueArray[0]=deltaPhi;
	valueArray[1]=deltaEta;
	valueArray[2]=G_PT_Value;
	valueArray[3]=ZT_Value;
	valueArray[4]=XI_Value;
	valueArray[5]=zVertex;
	valueArray[6]=identifier;
	valueArray[7]=evtPlaneCategory;
	valueArray[8]=fCent;

	if(G_PT_Value>=ClusterEcut)
	{
		fCorrVsManyThings  ->Fill(valueArray,Weight);
		//if(identifier==0)fCorrVsManyThingsME->Fill(valueArray,Weight);

		//..Histograms to test the binning
		fHistBinCheckPt[identifier] ->Fill(G_PT_Value,Weight);
		fHistBinCheckZt[identifier] ->Fill(ZT_Value,Weight);
		fHistBinCheckXi[identifier] ->Fill(XI_Value,Weight);

	/*
	 * This part is old and now substituted by the fCorrVsManyThings THn
	 *
	 * //..Fill 2D Histograms for certain event conditions
		for(Int_t i=0;i<10;i++)
		{
			for(Int_t j=0;j<kNvertBins;j++)
			{
			//no	if(i<kNoGammaBins && G_PT_Value>=fArray_G_Bins[i] && G_PT_Value<fArray_G_Bins[i+1])
				if(i<kNoGammaBins && G_PT_Value>fArray_G_Bins[i] && G_PT_Value<=fArray_G_Bins[i+1])
				{
					if(zVertex>fArrayNVertBins[j] && zVertex<=fArrayNVertBins[j+1])fHistDEtaDPhiG[identifier][i][j]->Fill(deltaPhi,deltaEta,Weight);
					//if(j==0 && i ==0 && identifier==1)fCorrVsManyThings->Fill(valueArray);
					if(j==0)fHistptAssHadronG[identifier][i]->Fill(TrackVec->Pt(),Weight);
					if(j==0)fHistptTriggG[identifier][i]    ->Fill(G_PT_Value,Weight);
				}
				if(i<kNoZtBins && ZT_Value>=fArray_ZT_Bins[i]  && ZT_Value<fArray_ZT_Bins[i+1])
				{
					if(zVertex>=fArrayNVertBins[j] && zVertex<fArrayNVertBins[j+1])fHistDEtaDPhiZT[identifier][i][j]   ->Fill(deltaPhi,deltaEta,Weight);
					if(j==0)fHistptAssHadronZt[identifier][i]->Fill(TrackVec->Pt(),Weight);
					if(j==0)fHistptTriggZt[identifier][i]    ->Fill(G_PT_Value,Weight);
				}
				if(i<kNoXiBins && XI_Value>=fArray_XI_Bins[i]  && XI_Value<fArray_XI_Bins[i+1])
				{
					if(zVertex>=fArrayNVertBins[j] && zVertex<fArrayNVertBins[j+1])fHistDEtaDPhiXI[identifier][i][j]   ->Fill(deltaPhi,deltaEta,Weight);
					if(j==0)fHistptAssHadronXi[identifier][i]->Fill(TrackVec->Pt(),Weight);
					if(j==0)fHistptTriggXi[identifier][i]    ->Fill(G_PT_Value,Weight);
				}
			}
		}
		*/
	}
}
//________________________________________________________________________
void AliAnalysisTaskGammaHadron::FillQAHisograms(Int_t identifier,AliClusterContainer* clusters,AliVCluster* caloCluster,AliVParticle* TrackVec)
{
	TLorentzVector caloClusterVec;
	clusters->GetMomentum(caloClusterVec,caloCluster);
	AliTLorentzVector aliCaloClusterVec = AliTLorentzVector(caloClusterVec); //..can acess phi from

	//..Get leading cluster
	AliVCluster* leadingClus = clusters->GetLeadingCluster("");  //"e" is energy, "" is Et

	//..ID array -  1 - leading, 2- track matched, 3 - leading & track matched
	Int_t gammaInfo=0;
	if(caloCluster==leadingClus)gammaInfo=1;
	if(DetermineMatchedTrack(caloCluster))gammaInfo=2;
	if(gammaInfo==2 && caloCluster==leadingClus)gammaInfo=3;

	//Eg, lambda0,NLM, ncells, distance to bad ,e/p, Mgg
	Double_t valueArray[6];
	valueArray[0]=caloCluster->GetNonLinCorrEnergy();
	valueArray[1]=caloCluster->GetM02();
	valueArray[2]=caloCluster->GetNExMax();
	valueArray[3]=caloCluster->GetNCells();
	valueArray[4]=caloCluster->GetDistanceToBadChannel();
	valueArray[5]=gammaInfo;
	//valueArray[5]=0;//E/p
	//valueArray[6]=130;//m_gg
	//aliCaloClusterVec.Phi_0_2pi()*fRtoD
	//caloClusterVec.Eta()


	if(identifier==0 && fPlotQA==1)fClusterProp->Fill(valueArray); //..all clusters - no cuts
	//if(identifier==2 && fPlotQA==1)fClusterProp->Fill(valueArray); //..accepted clusters - no cuts


	/*do similar test here?*/fHistDEtaDPhiGammaQA[identifier] ->Fill(caloClusterVec.Eta(),aliCaloClusterVec.Phi_0_2pi()*fRtoD);
	if(TrackVec)             fHistDEtaDPhiTrackQA[identifier] ->Fill(TrackVec->Eta(),TrackVec->Phi()*fRtoD);
	fHistCellsCluster[identifier] ->Fill(caloCluster->GetNonLinCorrEnergy(),caloCluster->GetNCells());
	fHistClusterTime[identifier]  ->Fill(caloCluster->GetTOF()*1000000000,caloCluster->GetNonLinCorrEnergy());
}
//
// Accept cluster for analysis. More cuts besides in ApplyClusterCuts and ApplyKinematicCuts
//
//
//________________________________________________________________________
Bool_t AliAnalysisTaskGammaHadron::AccClusterForAna(AliClusterContainer* clusters, AliVCluster* caloCluster)
{
	TLorentzVector caloClusterVec;
	clusters->GetMomentum(caloClusterVec,caloCluster);
	//!!!! eventually transform to AliTLorentzvector

	//..Accepts clusters if certain conditions are fulfilled
	Bool_t Accepted=1; //..By default accepted

	//!!double check these cuts carefully with the experts!!
	//-----------------------------
	//..at least 2 cells in cluster
	if(caloCluster->GetNCells()<2)
	{
		//..Reject the cluster as a good candidate for your analysis
		return 0;
	}
	//-----------------------------
	//..number of local maxima should be 1 or 0 (for cluster splitting this has to be changed)
	if(caloCluster->GetNExMax()>fMaxNLM)
	{
		//..Reject the cluster as a good candidate for your analysis
		return 0;
	}
	//-----------------------------
	//..cut on the cluster shape
	if(fClShapeMin>0 && fClShapeMax>0
	   && (caloCluster->GetM02()<fClShapeMin || caloCluster->GetM02()>fClShapeMax))
	{
		//..Reject the cluster as a good candidate for your analysis
		return 0;
	}
	//-----------------------------
	//..remove clusters with a matched track
	if(fRmvMTrack==1 && DetermineMatchedTrack(caloCluster))
	{
		return 0;
	}
	//-----------------------------
	//..Fiducial volume cut. If it is located neither in EMCal nor in DCal reject
	 if(!fFiducialCuts->IsInFiducialCut(caloClusterVec.Eta(),caloClusterVec.Phi(),AliFiducialCut::kEMCAL) &&
			 !fFiducialCuts->IsInFiducialCut(caloClusterVec.Eta(),caloClusterVec.Phi(),AliFiducialCut::kDCAL)
	 )
	 {
		 return 0;
	 }
	 //-----------------------------
	 //..Fiducial volume cut II. Cuts on the distance to the EMCal border last+first row and last+first collumn
	 //fFiducialCellCut->SetNumberOfCellsFromEMCALBorder(1); //ELI this could be momentum dependent and also different for merged clusters!!!
	 if(!fFiducialCellCut->CheckCellFiducialRegion(fGeom,caloCluster,fCaloCells))
	 {
		 return 0;
	 }
	 //-----------------------------
	 //..Do we need a distance to bad channel cut?
	 /*if(caloCluster->GetDistanceToBadChannel()<2)
	 {
		 return 0;
	 }*/

	 //-----------------------------
	 //..Do we need a distance to the border cut?
	 //caloCluster->???()


	 return Accepted;
}
//________________________________________________________________________
Double_t AliAnalysisTaskGammaHadron::DeltaPhi(AliTLorentzVector ClusterVec,AliVParticle* TrackVec)
{
	Double_t Phi_g = ClusterVec.Phi_0_2pi();
	Double_t Phi_h = TrackVec->Phi();

	Double_t dPhi = -999;
	Double_t pi = TMath::Pi();

	dPhi = Phi_g-Phi_h;
	//--shift the second peak over the fist peak: \--�--/   --> -�--
	//--to create a histogram that starts at -pi/2 and ends at 3/2pi
	if (dPhi <= -TMath::Pi()/2)    dPhi += 2*pi;
	if (dPhi > 3.0*TMath::Pi()/2.0)dPhi -= 2*pi;

	//--change from rad to degree:
	dPhi*= fRtoD;

	return dPhi;
}
//________________________________________________________________________
Double_t AliAnalysisTaskGammaHadron::DeltaPhi(AliTLorentzVector ClusterVec,Double_t phi_EVP)
{
	Double_t phi_g = ClusterVec.Phi_0_2pi();

	Double_t dPhi = -999;
	Double_t pi = TMath::Pi();

	dPhi = phi_g-phi_EVP;
	//--shift the second peak over the fist peak: \--�--/   --> -�--
	//--to create a histogram that starts at -pi/2 and ends at 3/2pi
	if (dPhi <= -TMath::Pi()/2)    dPhi += 2*pi;
	if (dPhi > 3.0*TMath::Pi()/2.0)dPhi -= 2*pi;

	//--change from rad to degree:
	dPhi*= fRtoD;

	return dPhi;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskGammaHadron::DetermineMatchedTrack(AliVCluster* caloCluster)
{
	Bool_t foundTrackMatched=0;
	Int_t Ntrks = caloCluster->GetNTracksMatched();
	if(Ntrks==0) return foundTrackMatched; //..if no matched track removal is wanted set it to 0.

	//..loop over matched tracks
	for (Int_t i = 0; i < Ntrks; ++i)
	{
		AliVTrack* track = static_cast<AliVTrack*>(caloCluster->GetTrackMatched(i));

		Double_t veta = track->GetTrackEtaOnEMCal();
		Double_t vphi = track->GetTrackPhiOnEMCal();

		Float_t pos[3] = {0};
		caloCluster->GetPosition(pos);
		TVector3 cpos(pos);
		Double_t ceta     = cpos.Eta();
		Double_t cphi     = cpos.Phi();
		Double_t etadiff=veta-ceta;
		Double_t phidiff=TVector2::Phi_mpi_pi(vphi-cphi);
		if(caloCluster->GetNCells()==2)fHistMatchEtaPhiAllCl2->Fill(etadiff, phidiff);
		if(caloCluster->GetNCells()==3)fHistMatchEtaPhiAllCl3->Fill(etadiff, phidiff);
		if(caloCluster->GetNCells()==4)fHistMatchEtaPhiAllCl4->Fill(etadiff, phidiff);

		//?  // check if track also points to cluster
		//?   Int_t cid = track->GetEMCALcluster();
		//?   if (fDoTrackClus && (cid != icluster)) return energyclus;

		//..check if the track was matchec within a stricter criteria
		if(TMath::Abs(etadiff)<fTrackMatchEta && TMath::Abs(phidiff)<fTrackMatchPhi)
		{
			if (track->GetLabel() > fMinMCLabel)
			{
				foundTrackMatched=1;
				break;
			}
		}
	}
	return foundTrackMatched;
}
//________________________________________________________________________
Double_t AliAnalysisTaskGammaHadron::GetEff(AliTLorentzVector ClusterVec)
{
	Double_t DetectionEff=1;

	/*
	 *
	 * Do something here with the input efficiency histograms
	 *
  THnF                      *fHistEffGamma;            // input efficiency for trigger particles
  THnF                      *fHistEffHadron;           // input efficiency for associate particles
	 *
	 *
	 */


	return DetectionEff;
}
//the end
