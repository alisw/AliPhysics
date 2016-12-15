//
// Task to estimate the number of gamma-hadron
// statistic available in the Pb+Pb run.
//
// Author: E. Epple

#include <Riostream.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
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


using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskGammaHadron)
//
//  This class inherits from AliAnalysisTaskEmcal() -->fcent,fVertex,fNVertCont,fVertexSPD,fNVertSPDCont,fTriggers is defined in there
//  And AliAnalysisTaskEmcal() inherits from AliAnalysisTaskSE()
//________________________________________________________________________
AliAnalysisTaskGammaHadron::AliAnalysisTaskGammaHadron():
AliAnalysisTaskEmcal("AliAnalysisTaskGammaHadron", kTRUE),

fGammaOrPi0(0),fDoMixing(0),fSavePool(0),

fParticleLevel(kFALSE),fIsMC(kFALSE),
fPoolMgr(0x0),
fTriggerType(AliVEvent::kAnyINT), fMixingEventType(AliVEvent::kAnyINT),//
fHistEffGamma(0x0),fHistEffHadron(0x0),

fOutputList1(),fOutputList2(),fOutputList3(),fOutputListGamma(),fOutputListXi(),fOutputListZeta(),fEventPoolOutputList(),

fHistNoClusPt(0),fHistPi0(0),fHistBinCheckPt(0), fHistBinCheckZt(0), fHistBinCheckXi(0),
fHistDEtaDPhiGammaQA(0),fHistDEtaDPhiTrackQA(0),fHistCellsCluster(0),fHistClusterShape(0),fHistClusterTime(0),

fHPoolReady(0x0)
{
	//..Initialize by defult for
	//..AliAnalysisTaskGammaHadron(0,1);
	InitArrays();
}
//________________________________________________________________________
AliAnalysisTaskGammaHadron::AliAnalysisTaskGammaHadron(Bool_t InputGammaOrPi0,Bool_t InputDoMixing):
AliAnalysisTaskEmcal("AliAnalysisTaskGammaHadron", kTRUE),

fGammaOrPi0(0),fDoMixing(0),fSavePool(0),

fParticleLevel(kFALSE),fIsMC(kFALSE),
fPoolMgr(0x0),
fTriggerType(AliVEvent::kAnyINT), fMixingEventType(AliVEvent::kAnyINT),
fHistEffGamma(0x0),fHistEffHadron(0x0),

fOutputList1(),fOutputList2(),fOutputList3(),fOutputListGamma(),fOutputListXi(),fOutputListZeta(),fEventPoolOutputList(),

fHistNoClusPt(0),fHistPi0(0), fHistBinCheckPt(0), fHistBinCheckZt(0), fHistBinCheckXi(0),
fHistDEtaDPhiGammaQA(0),fHistDEtaDPhiTrackQA(0),fHistCellsCluster(0),fHistClusterShape(0),fHistClusterTime(0),

fHPoolReady(0x0)
{
	InitArrays();
	//..set input variables
	fGammaOrPi0        =InputGammaOrPi0;
	fDoMixing          =InputDoMixing;
}
//________________________________________________________________________
void AliAnalysisTaskGammaHadron::InitArrays()
{
	//..Initialize by defult for
	//..AliAnalysisTaskGammaHadron(0,1);

	//..set input variables
	fGammaOrPi0        =0; //= 0 ( Gamma analysis ), 1 (pi0 analysis)
	fDoMixing          =0; //= 0 (do only same event analyis with correct triggers), =1 (do event mixing)

	fDebug             =0; //set only 1 for debugging
	fSavePool          =0; //= 0 do not save the pool by default. Use the set function to do this.

	//..These two items are set in AliAnalysisTaskEmcal::RetrieveEventObjects()
	//fCent, zVertex

	for(Int_t i=0; i<kNIdentifier;i++)
	{
		fHistptAssHadron[i]= 0;
		fHistDEtaDPhiG[i]  = 0;
		fHistDEtaDPhiZT[i] = 0;
		fHistDEtaDPhiXI[i] = 0;
		fHistDpGh[i]       = 0;
	}

	fRtoD=180.0/TMath::Pi();

	//..Set some default values.
	//..if desired one can add a set function to
	//..set these values in the add task function
	static const Int_t NcentBins=8;
	Double_t centmix[NcentBins+1] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0};
	fMixBCent = new TAxis(NcentBins,centmix);

	static const Int_t NvertBins=8;
	Double_t zvtxmix[NvertBins+1] = {-10,-6,-4,-2,0,2,4,6,10};
	fMixBZvtx = new TAxis(NvertBins,zvtxmix);

	//..Raymond/Megan gives more mixed event yield - don't know about the quality though
	//fTrackDepth     = 100;      //Hanseul sets it to 100! Q:: is this good? Maximum number of tracks??
	fTrackDepth     = 50000;    //Raymonds/Megans value

	//..!!
	//.. fPoolSize is an input that is ignored in the PoolManager Anyway
	//fPoolSize       = 200;    //200 - hanseuls default value Q:: is this correct? Maximum number of events in pool
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
	fOutputList2    = new TList();
	fOutputList2    ->SetOwner();
	fOutputList2    ->SetName("pT_distr_of_the_ass_H_for_a_given_Gpt");
	fOutputList3    = new TList();
	fOutputList3    ->SetOwner();
	fOutputList3    ->SetName("DPhi_GH_GammaPt");
	fOutputListGamma= new TList();
	fOutputListGamma->SetOwner();
	fOutputListGamma->SetName("Different_Gamma_2DHistograms");
	fOutputListXi   = new TList();
	fOutputListXi   ->SetOwner();
	fOutputListXi   ->SetName("Different_Xi_2DHistograms");
	fOutputListZeta = new TList();
	fOutputListZeta ->SetOwner();
	fOutputListZeta ->SetName("Different_Zt_2DHistograms");
	fOutputListQA   = new TList();
	fOutputListQA   ->SetOwner();
	fOutputListQA   ->SetName("QA_histograms");

	//common bins for the histograms
	Int_t nbins[6] = {0};
	Double_t min[6] = {0};
	Double_t max[6] = {0};

	//settings for p_t cluster distributon
	nbins[0] = 31;
	min[0] = 0;
	max[0] = 31;
	//settings for p_t hadron distribution
	nbins[1] = 60;  //do 1/2 GeV bins so that you can see the 0.5 cut to set as a minimum pT to combine hadron and gamma
	min[1] = 0;
	max[1] = 30;
	//settings for delta phi (g-h) distribution
	nbins[2] = 50;
	min[2] = -90;
	max[2] = 270;
	//settings for delta eta (g-h) distribution
	nbins[3] = 80;
	min[3] = -2;
	max[3] = 2;
	//settings for phi distribution for QA
	nbins[4] = 76;
	min[4] = -10;
	max[4] = 370;
	//settings for eta distribution for QA
	nbins[5] = 80;
	min[5] = -1;
	max[5] = 1;
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//
	//   Create Histograms
	//
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//    Histograms for common use
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//..set the steps for the ZT and the XI histograms
	Double_t fZtStep =1.0/(7-1.0);  // Bin width for the zT histograms
	Double_t fXiStep =2.5/(8-1.0);  // Bin width for the Xi histograms

	Int_t Array_G_Bins[10]   ={0,5,7,9,11,14,17,22,30,100};
	Double_t Array_ZT_Bins[8]={0,fZtStep,2*fZtStep,3*fZtStep,4*fZtStep,5*fZtStep,6*fZtStep,100};
	Double_t Array_XI_Bins[9]={-100,0,fXiStep,2*fXiStep,3*fXiStep,4*fXiStep,5*fXiStep,6*fXiStep,100};

	for(Int_t i=0;i<10;i++)
	{
		fVector_G_Bins.push_back(Array_G_Bins[i]);
	}
	for(Int_t i=0;i<8;i++)
	{
		fVector_ZT_Bins.push_back(Array_ZT_Bins[i]);
	}
	for(Int_t i=0;i<9;i++)
	{
		fVector_XI_Bins.push_back(Array_XI_Bins[i]);
	}
	fNoGammaBins=fVector_G_Bins.size()-1;
	fNoZtBins   =fVector_ZT_Bins.size()-1;
	fNoXiBins   =fVector_XI_Bins.size()-1;

	//..Initialize
	fHistBinCheckPt       = new TH1*[kNIdentifier];
	fHistBinCheckZt       = new TH1*[kNIdentifier];
	fHistBinCheckXi       = new TH1*[kNIdentifier];
	fHistDEtaDPhiGammaQA  = new TH2*[kNIdentifier+2]; //Why +2??
	fHistDEtaDPhiTrackQA  = new TH2*[kNIdentifier+2];
	fHistCellsCluster     = new TH2*[kNIdentifier+2];
	fHistClusterShape     = new TH2*[kNIdentifier+2];
	fHistClusterTime      = new TH2*[kNIdentifier+2];

	for(Int_t i=0; i<kNIdentifier; i++)
	{
		fHistptAssHadron[i]  = new TH1*[fNoGammaBins+1]; // make a p_t histogram of the associated hadron for each p_T bin of the gamma
		fHistDpGh[i]         = new TH1*[fNoGammaBins+1]; // make a p_t histogram of the associated hadron for each p_T bin of the gamma
		fHistDEtaDPhiG[i]    = new TH2*[fNoGammaBins];
		fHistDEtaDPhiZT[i]   = new TH2*[fNoZtBins];
		fHistDEtaDPhiXI[i]   = new TH2*[fNoXiBins];
	}
	//.................................
	//..p_T^{Cluster} distribution under different conditions

	//..all clusters
	fHistNoClusPt = new TH1F(Form("fHistNoClusPt_Id%0d",1),Form("fHistNoClusPt_Id%0d",1), nbins[0], min[0], max[0]);
	fHistNoClusPt->GetXaxis()->SetTitle("p_{T}^{Calo Cluster}");
	fHistNoClusPt->GetYaxis()->SetTitle(Form("No. of Clusters [counts/%0.1f GeV/c]",fHistNoClusPt->GetBinWidth(0)));
	fOutputList1->Add(fHistNoClusPt);

	//test!!
	fHistPi0 = new TH1F(Form("fHistPi0_%0d",1),Form("fHistPi0_%0d",1), 500, 0, 0.5);
	fHistPi0->GetXaxis()->SetTitle("M_{#gamma#gamma}");
	fHistPi0->GetYaxis()->SetTitle("Entries");
	fOutput->Add(fHistPi0);


	//..by the identifier different histograms can be filled under different cut conditions
	//..Can eg. later be modified to contain certain delta phi or centrality bins
	for(Int_t identifier=0;identifier<kNIdentifier;identifier++)
	{
		fHistBinCheckPt[identifier] = new TH1F(Form("fHistBinCheckPt_%0d",identifier),Form("fHistBinCheckPt_%0d",identifier), 500, 0, 100);
		fHistBinCheckPt[identifier]->GetXaxis()->SetTitle("p_T^{#gamma}");
		fHistBinCheckPt[identifier]->GetYaxis()->SetTitle("Entries");
		fOutput->Add(fHistBinCheckPt[identifier]);

		fHistBinCheckZt[identifier] = new TH1F(Form("fHistBinCheckZt_%0d",identifier),Form("fHistBinCheckZt_%0d",identifier), 500, 0, 100);
		fHistBinCheckZt[identifier]->GetXaxis()->SetTitle("z_T^{#gamma-h}");
		fHistBinCheckZt[identifier]->GetYaxis()->SetTitle("Entries");
		fOutput->Add(fHistBinCheckZt[identifier]);

		fHistBinCheckXi[identifier] = new TH1F(Form("fHistBinCheckXi_%0d",identifier),Form("fHistBinCheckXi_%0d",identifier), 500, -100, 100);
		fHistBinCheckXi[identifier]->GetXaxis()->SetTitle("#xi^{#gamma-h}");
		fHistBinCheckXi[identifier]->GetYaxis()->SetTitle("Entries");
		fOutput->Add(fHistBinCheckXi[identifier]);

		for(Int_t i=0; i<fNoGammaBins; i++)
		{
			fHistDEtaDPhiG[identifier][i] = new TH2F(Form("fHistDEtaDPhiG%d_Id%d",i,identifier),Form("fHistDEtaDPhiG%d_Id%d",i,identifier),nbins[2],min[2],max[2],nbins[3],min[3],max[3]);
			fHistDEtaDPhiG[identifier][i]->GetXaxis()->SetTitle(Form("#Delta #phi^{#gamma-h} %0d<p_{T}^{#gamma}<%0d",fVector_G_Bins.at(i),fVector_G_Bins.at(i+1)));
			fHistDEtaDPhiG[identifier][i]->GetYaxis()->SetTitle("#Delta #eta^{#gamma-h}");
			fOutputListGamma->Add(fHistDEtaDPhiG[identifier][i]);
		}
		for(Int_t i=0; i<fVector_ZT_Bins.size()-1; i++)
		{
			fHistDEtaDPhiZT[identifier][i] = new TH2F(Form("fHistDEtaDPhiZT%d_Id%d",i,identifier),Form("fHistDEtaDPhiZT%d_Id%d",i,identifier),nbins[2],min[2],max[2],nbins[3],min[3],max[3]);
			fHistDEtaDPhiZT[identifier][i]->GetXaxis()->SetTitle(Form("#Delta #phi^{#gamma-h} %0.1f<z_{T}<%0.1f",fVector_ZT_Bins.at(i),fVector_ZT_Bins.at(i+1)));
			fHistDEtaDPhiZT[identifier][i]->GetYaxis()->SetTitle("#Delta #eta^{#gamma-h}");
			fOutputListZeta->Add(fHistDEtaDPhiZT[identifier][i]);
		}
		for(Int_t i=0; i<fVector_XI_Bins.size()-1; i++)
		{
			fHistDEtaDPhiXI[identifier][i] = new TH2F(Form("fHistDEtaDPhiXI%d_Id%d",i,identifier),Form("fHistDEtaDPhiXI%d_Id%d",i,identifier),nbins[2],min[2],max[2],nbins[3],min[3],max[3]);
			fHistDEtaDPhiXI[identifier][i]->GetXaxis()->SetTitle(Form("#Delta #phi^{#gamma-h} %0.1f<#xi<%0.1f",fVector_XI_Bins.at(i),fVector_XI_Bins.at(i+1)));
			fHistDEtaDPhiXI[identifier][i]->GetYaxis()->SetTitle("#Delta #eta^{#gamma-h}");
			fOutputListXi->Add(fHistDEtaDPhiXI[identifier][i]);
		}

		for(Int_t i=0; i<fNoGammaBins+1; i++)
		{
			//..check whether the max is the
			//..same as the histogram size
			Double_t BinWidth,BinValStart;

			if(i==0)
			{
				//..these are histograms over the full p_T gamma range (no binnings)
				BinValStart  = Array_G_Bins[0];
				BinWidth    =100; //to get the full length of the histogram
			}
			else
			{
				//..these are histograms for a certain p_T of the gamma
				BinValStart  = Array_G_Bins[i-1];
				BinWidth     = Array_G_Bins[i]-Array_G_Bins[i-1];
				//cout<<"BinStart: "<<fHistNoClusPtH->GetBinCenter(i)<<", width "<<BinWidth<<endl;
			}

			//.................................
			//..Initiate the delta phi and p_t associated hadron histograms for each p_t bin of the gamma!
			fHistptAssHadron[identifier][i] = new TH1F(Form("fHistAssHadron_pt_%0d_%0d",identifier,i),Form("fHistAssHadron_pt_%0d_%0d",identifier,i), nbins[1], min[1], max[1]);
			fHistptAssHadron[identifier][i]->GetXaxis()->SetTitle(Form("p_{T}^{assoc. h} %0.1f<p_{T}^{#gamma}<%0.1f",BinValStart,BinValStart+BinWidth));
			fHistptAssHadron[identifier][i]->GetYaxis()->SetTitle(Form("dN^{assoc. h}/dp_{T}^{h} [counts/%0.1f GeV/c]",fHistptAssHadron[identifier][i]->GetBinWidth(1)));
			fOutputList2->Add(fHistptAssHadron[identifier][i]);

			//.................................
			//..Initiate the delta phi and p_t associated hadron histograms for each p_t bin of the gamma!
			fHistDpGh[identifier][i] = new TH1F(Form("fHistDpGh_%0d_%0d",identifier,i),Form("fHistDpGh_%0d_%0d",identifier,i), nbins[2], min[2], max[2]);
			fHistDpGh[identifier][i]->GetXaxis()->SetTitle(Form("#Delta #phi^{#gamma-h} %0.1f<p_{T}^{#gamma}<%0.1f",BinValStart,BinValStart+BinWidth));
			fHistDpGh[identifier][i]->GetYaxis()->SetTitle(Form("dN^{#gamma-h}/#Delta #phi^{#gamma-h} [counts/%0.1f^{#circ}]",fHistDpGh[identifier][i]->GetBinWidth(1)));
			fOutputList3->Add(fHistDpGh[identifier][i]);
		}
	}

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   Tyler's Special Histograms
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	//...

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   Eliane's Special Histograms
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	//...
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   Special QA histograms (also to get more info what is going on in mixed event for trigger data)
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	for(Int_t identifier=0;identifier<kNIdentifier+2;identifier++)
	{
		//..geometrical hit distribution of clusters
		fHistDEtaDPhiGammaQA[identifier] = new TH2F(Form("fHistDEtaDPhiGammaQA%d_Id%d",0,identifier),Form("fHistDEtaDPhiGammaQA%d_Id%d",0,identifier),nbins[4],min[4],max[4],nbins[5],min[5],max[5]);
		fHistDEtaDPhiGammaQA[identifier]->GetXaxis()->SetTitle("#phi^{#gamma}");
		fHistDEtaDPhiGammaQA[identifier]->GetYaxis()->SetTitle("#eta^{#gamma}");
		fOutputListQA->Add(fHistDEtaDPhiGammaQA[identifier]);

		//..geometrical hit distribution of tracks
		fHistDEtaDPhiTrackQA[identifier] = new TH2F(Form("fHistDEtaDPhiTrackQA%d_Id%d",0,identifier),Form("fHistDEtaDPhiTrackQA%d_Id%d",0,identifier),nbins[4],min[4],max[4],nbins[5],min[5],max[5]);
		fHistDEtaDPhiTrackQA[identifier]->GetXaxis()->SetTitle("#phi^{hadron}");
		fHistDEtaDPhiTrackQA[identifier]->GetYaxis()->SetTitle("#eta^{hadron}");
		fOutputListQA->Add(fHistDEtaDPhiTrackQA[identifier]);

	    //..Cluster Info
		fHistCellsCluster[identifier] = new TH2F(Form("fHistCellsCluster%d_Id%d",0,identifier),Form("fHistCellsCluster%d_Id%d",0,identifier),240,0,30,50,0,50);
		fHistCellsCluster[identifier]->GetXaxis()->SetTitle("E^{cluster}");
		fHistCellsCluster[identifier]->GetYaxis()->SetTitle("N_{cells}");
		fOutputListQA->Add(fHistCellsCluster[identifier]);

		fHistClusterShape[identifier] = new TH2F(Form("fHistClusterShape%d_Id%d",0,identifier),Form("fHistClusterShape%d_Id%d",0,identifier),240,0,30,240,0,4);
		fHistClusterShape[identifier]->GetXaxis()->SetTitle("E^{cluster}");
		fHistClusterShape[identifier]->GetYaxis()->SetTitle("#lambda_{0}");
		fOutputListQA->Add(fHistClusterShape[identifier]);

		//..Time information
		fHistClusterTime[identifier] = new TH2F(Form("fHistClusterTime%d_Id%d",0,identifier),Form("fHistClusterTime%d_Id%d",0,identifier),20000,-100,100,80,0,40);
		fHistClusterTime[identifier]->GetXaxis()->SetTitle("time [ns]");
		fHistClusterTime[identifier]->GetYaxis()->SetTitle("pT");
		fOutputListQA->Add(fHistClusterTime[identifier]);

		//..
	}



	//..The END
	fOutput->Add(fOutputList1);
	fOutput->Add(fOutputList2);
	fOutput->Add(fOutputList3);
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
//________________________________________________________________________
Bool_t AliAnalysisTaskGammaHadron::Run()
{
	//..This function is called in AliAnalysisTaskEmcal::UserExec.
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::Run()"<<endl;

	//..Determine the trigger for the current event
	UInt_t fCurrentEventTrigger = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

	//..check here some basic properties of the event
	//is centrality and zvertex in range? - see if this is actually done in IsSelected in the EMCal Task

/*	cout<<"- - - - - - - - - - - - - "<<endl;
	cout<<"eventTrigger :"<<fCurrentEventTrigger<<endl;
	cout<<"fTriggerType :"<<fTriggerType<<endl;
	cout<<"fMixedEventType :"<<fMixingEventType<<endl;
	if(fCurrentEventTrigger & fTriggerType) cout<<" ---> contains EMCGA trigger"<<endl;
	if(fCurrentEventTrigger & fMixingEventType) cout<<" ---> contains kInt7 trigger"<<endl;
*/
	if(fCurrentEventTrigger & fTriggerType)
	{
		if(fCurrentEventTrigger & fMixingEventType)cout<<"contains both triggers!!"<<endl;
	}
	/*Char_t cr[100];
 	TBits bitNo;
 	bitNo.Set(fTriggerType,cr);
	//TString binary;
	//TString str(itoa(fCurrentEventTrigger,binary,2));
	cout<<"convert trigger to decimal: "<<cr<<endl;
	*/

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
	if(fDoMixing==1)   //if(fSameEventAnalysis==0)
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
		//..update pool only with tracks from event type fMixingEventType, (do not add tracks from triggered events)
		if(fCurrentEventTrigger & fMixingEventType)
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
	Int_t trackCounter=0;
	for(Int_t NoCluster1 = 0; NoCluster1 < NoOfClustersInEvent; NoCluster1++ )
	{
		cluster=(AliVCluster*) clusters->GetAcceptCluster(NoCluster1); //->GetCluster(NoCluster1);
		if(!cluster || !AccClusterForAna(clusters,cluster))continue; //check if the cluster is a good cluster
		//clusters->GetLeadingCluster("e");

		TLorentzVector CaloClusterVec;
		clusters->GetMomentum(CaloClusterVec, cluster);
		EffWeight_Gamma=GetEff(CaloClusterVec);

		FillQAHisograms(0,clusters,cluster,trackNULL,0,0);
		fHistNoClusPt->Fill(CaloClusterVec.Pt()); //the .pt only works for gammas (E=M) for other particle this is wrong
		//cout<<"Cluster number: "<<NoCluster1<<", pT = "<<CaloClusterVec.Pt()<<endl;
		//...........................................
		//..combine gammas with same event tracks
		GammaCounter++;
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
				trackCounter++;
				//cout<<"..Track number: "<<NoTrack<<", pT = "<<track->Pt()<<endl;

				//EffWeight_Hadron=GetEff(<TLorentzVector>track);
				FillGhHisograms(1,CaloClusterVec,track,2,0,Weight);

				if(GammaCounter==1)FillQAHisograms(1,clusters,cluster,track,2,0); //fill only once per track (first gamma) - good for each track
				if(trackCounter==1)FillQAHisograms(2,clusters,cluster,track,2,0); //fill only once per gamma (first track) - good for gamma distr.
				if(CaloClusterVec.Eta()>0)FillQAHisograms(3,clusters,cluster,track,2,0); //fill for each gamma?? makes sense?
				if(CaloClusterVec.Eta()<0)FillQAHisograms(4,clusters,cluster,track,2,0); //fill for each gamma?? makes sense?
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
				FillGhHisograms(0,CaloClusterVec,track,2,0,Weight);
			}
		}
		//...........................................
		//..double cluster loop for testing
		for( Int_t NoCluster2 = 0; NoCluster2 < NoOfClustersInEvent; NoCluster2++ )
		{
			if(NoCluster1!=NoCluster2)
			{
				cluster2=(AliVCluster*) clusters->GetAcceptCluster(NoCluster2);
				if(!cluster2 || !AccClusterForAna(clusters,cluster2))continue; //check if the cluster is a good cluster

				TLorentzVector CaloClusterVec2;
				TLorentzVector CaloClusterVecpi0;
				//old framework		cluster2->GetMomentum(CaloClusterVec2, fVertex);
				clusters->GetMomentum(CaloClusterVec2, cluster2);
				if(cluster2->E()>2 && cluster->E()>2)
				{
					CaloClusterVecpi0=CaloClusterVec+CaloClusterVec2;
					fHistPi0->Fill(CaloClusterVecpi0.M());
				}
			}
		}
	}
	return GammaCounter;
}
//________________________________________________________________________
Int_t AliAnalysisTaskGammaHadron::CorrelatePi0AndTrack(AliParticleContainer* tracks,TObjArray* bgTracksArray,Bool_t SameMix, Double_t InputWeight)
{
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::CorrelatePi0AndTrack()"<<endl;

	//**This is Tyler's world
	//**copied from CorrelateClusterAndTrack
	//**ready to be modified and used for something fabulous

	//...........................................
	//--Do cluster loop.
	AliClusterContainer* clusters  = GetClusterContainer(0);  //how do I know which cells are selected
	if (!clusters) return 0;
	Int_t NoOfClustersInEvent =clusters->GetNClusters();
	Int_t nAccClusters = 0;
	Int_t nAccPi0Clusters = 0;
	Double_t Pi0Mass = 0.13487;
	Double_t Pi0Window = 0.03; //can eventually modulate this based on pT of Pi0 candidate!
	Double_t EffWeight_Gamma;
	Double_t EffWeight_Hadron;
	Double_t Weight;    //weight to normalize mixed and same event distributions individually

	AliVCluster* cluster = 0;
	AliVCluster* cluster2= 0;

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
			//old framework			cluster->GetMomentum(CaloClusterVec, fVertex);
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
		//old framework	cluster->GetMomentum(CaloClusterVec, fVertex);
		clusters->GetMomentum(CaloClusterVec,cluster);

		for( Int_t NoCluster2 = 0; NoCluster2 < NoOfClustersInEvent; NoCluster2++ )
		{
			if(NoCluster1!=NoCluster2)
			{
				cluster2=(AliVCluster*) clusters->GetAcceptCluster(NoCluster2);
				if(!cluster2 || !AccClusterForAna(clusters,cluster2))continue; //check if the cluster is a good cluster

				TLorentzVector CaloClusterVec2;
				TLorentzVector CaloClusterVecpi0;
				//old framework			cluster2->GetMomentum(CaloClusterVec2, fVertex);
				clusters->GetMomentum(CaloClusterVec,cluster2);
				if(cluster2->E()>2 && cluster->E()>2)
				{
					CaloClusterVecpi0=CaloClusterVec+CaloClusterVec2;


					//here I don't really know what to do in your case
					//eff of pi0? or eff of gamma? or some mix up of the two ?
					EffWeight_Gamma=GetEff(CaloClusterVecpi0);//currently just assigns 1!!! need eventually to input Pi0 efficiency histogram

					fHistNoClusPt->Fill(CaloClusterVecpi0.Pt()); //the .pt only works for gammas (E=M) for other particle this is wrong

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
							FillGhHisograms(1,CaloClusterVecpi0,track,2,0,Weight);
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
							FillGhHisograms(0,CaloClusterVecpi0,track,2,0,Weight);
						}
					}
				}
			}
		}
	}
	return nAccPi0Clusters/2;
}
//________________________________________________________________________
void AliAnalysisTaskGammaHadron::FillGhHisograms(Int_t identifier,TLorentzVector ClusterVec,AliVParticle* TrackVec, Double_t ClusterEcut, Double_t TrackPcut, Double_t Weight)
{
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::FillGhHisograms()"<<endl;

	//..This function fills several histograms under different cut conditions.
	//..it is run within a cluster{ track{}} loop to get all combinations.

	//..A word to the weight - for mixed events it devides by the number of events in the current pool 1/nEvents
	//..                     - for same events you devide by the number of triggers
	//..                     - for both you have to take into account the efficiency of your correlated pair
	Double_t deltaEta   = ClusterVec.Eta()-TrackVec->Eta();
	Double_t deltaPhi   = DeltaPhi(ClusterVec,TrackVec);
	Double_t G_PT_Value = ClusterVec.Pt();
	//Double_t ZT_Value   = TMath::Cos(deltaPhi*1/fRtoD)*TrackVec->P()/ClusterVec.P(); //   TrackVec->Pt()/G_PT_Value;
	Double_t ZT_Value   = TrackVec->Pt()/G_PT_Value; //   TrackVec->Pt()/G_PT_Value;
	//..Careful here: usually this is done for an opening angle (hadron-jet axis) of less than 90�. Due to
	//..resolution momentum smearing (our guess - check that!) there are particles appearing at angles greater than 90�
	Double_t XI_Value;
	if(ZT_Value>0)
	{
		XI_Value   = TMath::Log(1.0/ZT_Value);
	}
	/*else
    {
    	  XI_Value   = TMath::Log(1.0/(-1.0)*ZT_Value);
    }*/


	if(G_PT_Value>=ClusterEcut && TrackVec->Pt()>=TrackPcut)
	{
		fHistptAssHadron[identifier][0]->Fill(TrackVec->Pt(),Weight);
		fHistDpGh[identifier][0]       ->Fill(deltaPhi,Weight);

		//..Histograms to test the binning
		fHistBinCheckPt[identifier] ->Fill(G_PT_Value,Weight);
		fHistBinCheckZt[identifier] ->Fill(ZT_Value,Weight);
		fHistBinCheckXi[identifier] ->Fill(XI_Value,Weight);

		//..Fill 2D Histograms for certain event conditions
		for(Int_t i=0;i<10;i++)
		{
			if(i<fVector_G_Bins.size()-1  && G_PT_Value>=fVector_G_Bins.at(i) && G_PT_Value<fVector_G_Bins.at(i+1)) fHistDEtaDPhiG[identifier][i] ->Fill(deltaPhi,deltaEta,Weight);
			if(i<fVector_ZT_Bins.size()-1 && ZT_Value>=fVector_ZT_Bins.at(i)  && ZT_Value<fVector_ZT_Bins.at(i+1))  fHistDEtaDPhiZT[identifier][i]->Fill(deltaPhi,deltaEta,Weight);
			if(i<fVector_XI_Bins.size()-1 && XI_Value>=fVector_XI_Bins.at(i)  && XI_Value<fVector_XI_Bins.at(i+1))  fHistDEtaDPhiXI[identifier][i]->Fill(deltaPhi,deltaEta,Weight);
		}

		//maybe this stuff here is a little outdated
		//       ||
		//       ||
		//       \/
		for(Int_t i=1;i<fNoGammaBins+1;i++)
		{
			//..look in each p_T bin of the gamma (2D fHistDEtaDPhiG) how the pt
			//..distribution of the associated hadron looks like.
			Double_t BinValStart  = fVector_G_Bins.at(i-1);
			Double_t BinWidth     = fVector_G_Bins.at(i)-fVector_G_Bins.at(i-1);

			//--right now the ranges are defined via different histogram bins of the gmma p_t histogram
			//--that might be changed later.
			if(G_PT_Value>=BinValStart && G_PT_Value<BinValStart+BinWidth)
			{
				fHistptAssHadron[identifier][i]->Fill(TrackVec->Pt(),Weight);
				fHistDpGh[identifier][i]       ->Fill(deltaPhi,Weight);
			}
		}
	}
}
//________________________________________________________________________
void AliAnalysisTaskGammaHadron::FillQAHisograms(Int_t identifier,AliClusterContainer* clusters,AliVCluster* caloCluster,AliVParticle* TrackVec, Double_t ClusterEcut, Double_t TrackPcut)
{
	TLorentzVector caloClusterVec;
	clusters->GetMomentum(caloClusterVec,caloCluster);

	/*do similar test here?*/fHistDEtaDPhiGammaQA[identifier] ->Fill(caloClusterVec.Phi()*fRtoD,caloClusterVec.Eta());
	if(TrackVec)             fHistDEtaDPhiTrackQA[identifier] ->Fill(TrackVec->Phi()*fRtoD,TrackVec->Eta());

	fHistCellsCluster[identifier] ->Fill(caloCluster->GetHadCorrEnergy(),caloCluster->GetNCells());
	fHistClusterShape[identifier] ->Fill(caloCluster->GetHadCorrEnergy(),caloCluster->GetM02());
	fHistClusterTime[identifier]  ->Fill(caloCluster->GetTOF()*1000000000,caloCluster->GetHadCorrEnergy());

}
//________________________________________________________________________
Bool_t AliAnalysisTaskGammaHadron::AccClusterForAna(AliClusterContainer* clusters, AliVCluster* caloCluster)
{
	TLorentzVector caloClusterVec;
	clusters->GetMomentum(caloClusterVec,caloCluster);
    Double_t deltaPhi=2;   //..phi away from detector edges.
    Double_t deltaEta=0.0; //..eta away from detector edges.

	//..Accepts clusters if certain conditions are fulfilled
	Bool_t Accepted=1; //..By default accepted

	//!!double check these cuts carefully with the experts!!
	if(caloCluster->GetNCells()<2)
	{
		//..Reject the cluster as a good candidate for your analysis
		Accepted=0;
	}
	//..If not in EMCal Phi acceptance - reject
	if(caloClusterVec.Phi()*fRtoD<(80+deltaPhi) || caloClusterVec.Phi()*fRtoD>(187-deltaPhi))
	{
		if(caloClusterVec.Phi()*fRtoD>(260+deltaPhi) && caloClusterVec.Phi()*fRtoD<(327-deltaPhi))
		{
			  //..if instead in DCal -> OK
		}
		else  //..if not reject
		{
			Accepted=0;
		}
	}
	//..If not in EMCal Eta acceptance - reject
	if(caloClusterVec.Eta()<(-0.7+deltaEta) || caloClusterVec.Eta()>(0.7-deltaEta))
	{
		if(caloClusterVec.Eta()>(0.22+deltaEta) && caloClusterVec.Eta()<(0.7-deltaEta))
		{
			  //..if instead in DCal -> OK
		}
		else  //..if not reject
		{
			Accepted=0;
		}
	}
	return Accepted;
}
//________________________________________________________________________
Double_t AliAnalysisTaskGammaHadron::DeltaPhi(TLorentzVector ClusterVec,AliVParticle* TrackVec)
{
	Double_t Phi_g = ClusterVec.Phi();
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
Double_t AliAnalysisTaskGammaHadron::GetEff(TLorentzVector ClusterVec)
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
