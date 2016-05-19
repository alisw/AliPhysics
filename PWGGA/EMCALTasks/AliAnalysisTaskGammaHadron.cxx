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
//#include "AliPool.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskGammaHadron)
//
//  This class inherits from AliAnalysisTaskEmcalJet() do we need this?
//  And AliAnalysisTaskEmcalJet() inherits from AliAnalysisTaskEmcal() -->fcent,fVertex is defined in there
//  And AliAnalysisTaskEmcal() inherits from AliAnalysisTaskSE()
//________________________________________________________________________
AliAnalysisTaskGammaHadron::AliAnalysisTaskGammaHadron():
AliAnalysisTaskEmcal("AliAnalysisTaskGammaHadron", kTRUE),

fGammaOrPi0(0),fSameEventAnalysis(1),
fCellEnergyCut(0.05),fMaxCellsInCluster(50),

fParticleLevel(kFALSE),
fIsMC(kFALSE),
fCentAlt(),
fPoolMgr(0x0),
fHistEffGamma(0x0),
fHistEffHadron(0x0),

fOutputList1(),fOutputList2(),fOutputList3(),fOutputListGamma(),fOutputListXi(),fOutputListZeta(),

fHistNoClusPtTrigger(0),fHistNoClusPt(0),fHistNoClusPtH(0),fHistPi0(0),
fHistDEtaDPhiGammaQA(0),fHistDEtaDPhiTrackQA(0),

fHPoolReady(0x0)
{
	InitArrays();
}
//________________________________________________________________________
AliAnalysisTaskGammaHadron::AliAnalysisTaskGammaHadron(Bool_t InputGammaOrPi0,Bool_t InputSameEventAnalysis) :
AliAnalysisTaskEmcal("AliAnalysisTaskGammaHadron", kTRUE),

fGammaOrPi0(0),fSameEventAnalysis(1),
fCellEnergyCut(0.05),fMaxCellsInCluster(50),

fParticleLevel(kFALSE),
fIsMC(kFALSE),
fCentAlt(),
fPoolMgr(0x0),
fHistEffGamma(0x0),
fHistEffHadron(0x0),

fOutputList1(),fOutputList2(),fOutputList3(),fOutputListGamma(),fOutputListXi(),fOutputListZeta(),

fHistNoClusPtTrigger(0),fHistNoClusPt(0),fHistNoClusPtH(0),fHistPi0(0),
fHistDEtaDPhiGammaQA(0),fHistDEtaDPhiTrackQA(0),

fHPoolReady(0x0)
{
	InitArrays();

	//set input variables
	fSameEventAnalysis =InputSameEventAnalysis;
	fGammaOrPi0        =InputGammaOrPi0;
}
//________________________________________________________________________
void AliAnalysisTaskGammaHadron::InitArrays()
{
	//Initialize by defult for
	//InputGammaOrPi0 = 0 ( Gamma analysis )
	//InputSameEventAnalysis = 1 (Same event analysis, currently this is only used to throw out event in th erun function)
	//AliAnalysisTaskGammaHadron(0,1);
	//set input variables
	fSameEventAnalysis =1;
	fGammaOrPi0        =1;
	fDebug             =0; //set only 1 for debugging
	fUsePerTrigWeight  =0; //plot histograms /gamma yield (1) plot absolute values (0)


	// Default constructor.
	fAODfilterBits[0] = 0;
	fAODfilterBits[1] = 0;

	for(Int_t i=0; i<kNIdentifier;i++)
	{
		fHistptAssHadron[i]= 0;
		fHistDEtaDPhiG[i]  = 0;
		fHistDEtaDPhiZT[i] = 0;
		fHistDEtaDPhiXI[i] = 0;
	}
	fHistDpGh[0]= 0;
	fHistDpGh[1]= 0;
	fHistDpGh[2]= 0;

	fRtoD=180.0/TMath::Pi();

	//Set some default values.
	//if desired one can add a set function to
	//set these values in the add task function
    static const Int_t NcentBins=8;
    Double_t centmix[NcentBins+1] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0};
	fMixBCent = new TAxis(NcentBins,centmix);

	static const Int_t NvertBins=8;
    Double_t zvtxmix[NvertBins+1] = {-10,-6,-4,-2,0,2,4,6,10};
    fMixBZvtx = new TAxis(NvertBins,zvtxmix);

    fCentMethodAlt = "V0M";
    fTrackDepth     = 100;     //Hanseul sets it to 100! Q:: is this good? Maximum number of tracks??
    fPoolSize       = 200;     //hanseuls default value Q:: is this correct? Maximum number of events in pool

    //member function of AliAnalysisTaskEmcal
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
	//
	//   Create mixed event pools
	//
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	InitEventMixer();
	/*
	// Centrality mixing axis
	Int_t nCentBins=fMixBCent->GetNbins();
	Double_t centBins[nCentBins+1];
	centBins[0] = fMixBCent->GetBinLowEdge(1);
	for(Int_t i=1; i<=nCentBins; i++)
	{
		centBins[i] = fMixBCent->GetBinUpEdge(i);
	}
	// Z-vertex mixing axis
	Int_t nZvtxBins=fMixBZvtx->GetNbins();
	Double_t zvtxbin[nZvtxBins+1];
	zvtxbin[0] = fMixBZvtx->GetBinLowEdge(1);
	for(Int_t i=1; i<=nZvtxBins; i++)
	{
		zvtxbin[i] = fMixBZvtx->GetBinUpEdge(i);
	}

	cout<<"1-------------"<<endl;

	fHPoolReady = new TH2F("fHPoolReady","mixing started", nZvtxBins, zvtxbin, nCentBins, centBins);
	fOutput->Add(fHPoolReady);
	cout<<"2-------------"<<endl;
    */

	// Creates the fOutput list
	// Creates histograms:
	//(fHistZVertex, fHistEventRejection, fHistEventRejection, fHistEventCount, fHistCentrality, fHistEventPlane)
	AliAnalysisTaskEmcal::UserCreateOutputObjects();

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//
	//   Define sublists/folders for a better organisation of the figures
	//
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
	min[2] = -100;
	max[2] = 300;
	//settings for delta eta (g-h) distribution
	nbins[3] = 30;
	min[3] = -2;
	max[3] = 2;
	//settings for delta phi (g-h) distribution
	nbins[4] = 50;
	min[4] = -10;
	max[4] = 370;
	//settings for eta distribution for QA
	nbins[5] = 30;
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

	//..Initialize
	fHistNoClusPtH        = new TH1*[kNIdentifier];
	fHistDEtaDPhiGammaQA  = new TH2*[kNIdentifier+2];
	fHistDEtaDPhiTrackQA  = new TH2*[kNIdentifier+2];

	for(Int_t i=0; i<kNIdentifier; i++)
	{
		fHistptAssHadron[i]  = new TH1*[kNDPhistos]; // make a p_t histogram of the associated hadron for each p_T bin of the gamma
		fHistDpGh[i]         = new TH1*[kNDPhistos]; // make a p_t histogram of the associated hadron for each p_T bin of the gamma

		fHistDEtaDPhiG[i]    = new TH2*[fVector_G_Bins.size()-1];
		fHistDEtaDPhiZT[i]   = new TH2*[fVector_ZT_Bins.size()-1];
		fHistDEtaDPhiXI[i]   = new TH2*[fVector_XI_Bins.size()-1];
	}
	//.................................
	//..p_T^{Cluster} distribution under different conditions

	//..all clusters
	fHistNoClusPt = new TH1F(Form("fHistNoClusPt_%0d",1),Form("fHistNoClusPt_%0d",1), nbins[0], min[0], max[0]);
	fHistNoClusPt->GetXaxis()->SetTitle("p_{T}^{Calo Cluster}");
	fHistNoClusPt->GetYaxis()->SetTitle(Form("No. of Clusters [counts/%0.1f GeV/c]",fHistNoClusPt->GetBinWidth(0)));
	fOutputList1->Add(fHistNoClusPt);

	//..trigger histogram only for internal purposes - no nead to save
	//fHistNoClusPtTrigger = new TH1F((const TH1F)fHistNoClusPt);

	//test!!
	fHistPi0 = new TH1F(Form("fHistPi0_%0d",1),Form("fHistPi0_%0d",1), 500, 0, 0.5);
	fHistPi0->GetXaxis()->SetTitle("M_{#gamma#gamma}");
	fHistPi0->GetYaxis()->SetTitle("Entries");
	fOutput->Add(fHistPi0);

	//..by the identifier different histograms can be filled under different cut conditions
	//..Can eg. later be modified to contain certain delta phi or centrality bins
	for(Int_t identifier=0;identifier<kNIdentifier;identifier++)
	{
		//..clusters p_T if there was a hadron present
		fHistNoClusPtH[identifier] = new TH1F(Form("fHistNoClusPtH_%0d",identifier),Form("fHistNoClusPtH_%0d",identifier), nbins[0], min[0], max[0]);
		fHistNoClusPtH[identifier]->GetXaxis()->SetTitle("p_{T}^{Calo Cluster}");
		fHistNoClusPtH[identifier]->GetYaxis()->SetTitle(Form("No. of Clus. with h [counts/%0.1f GeV/c]",fHistNoClusPtH[identifier]->GetBinWidth(0)));
		fOutput->Add(fHistNoClusPtH[identifier]);

		for(Int_t i=0; i<fVector_G_Bins.size()-1; i++)
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

		for(Int_t i=0; i<kNDPhistos+1; i++)
		{
			//..check whether the max is the
			//..same as the histogram size
			Double_t BinWidth,BinValStart;

			if(i==0)
			{
				//..these are histograms over the full p_T gamma range (no binnings)
				BinWidth    = fHistNoClusPtH[0]->GetBinWidth(1);
				BinValStart = fHistNoClusPtH[0]->GetXaxis()->GetBinCenter(1)-BinWidth/2.0;
				BinWidth    =100; //to get the full length of the histogram
			}
			else
			{
				//..these are histograms for a certain p_T of the gamma
				BinWidth    = fHistNoClusPtH[0]->GetBinWidth(i);
				BinValStart = fHistNoClusPtH[0]->GetXaxis()->GetBinCenter(i)-BinWidth/2.0;
				//cout<<"BinStart: "<<fHistNoClusPtH->GetBinCenter(i)<<", width "<<BinWidth<<endl;
			}

			//.................................
			//..Initiate the delta phi and p_t associated hadron histograms for each p_t bin of the gamma!
			fHistptAssHadron[identifier][i] = new TH1F(Form("fHistAssHadron_pt_%0d_%0d",identifier,i),Form("fHistAssHadron_pt_%0d_%0d",identifier,i), nbins[1], min[1], max[1]);
			fHistptAssHadron[identifier][i]->GetXaxis()->SetTitle(Form("p_{T}^{assoc. h} %0.1f<p_{T}^{#gamma}<%0.1f",BinValStart,BinValStart+BinWidth));
			fHistptAssHadron[identifier][i]->GetYaxis()->SetTitle(Form("dN^{assoc. h}/dp_{T}^{h} [counts/%0.1f GeV/c]",fHistptAssHadron[identifier][i]->GetBinWidth(1)));
			if(i!=0)fOutputList2->Add(fHistptAssHadron[identifier][i]);

			//.................................
			//..Initiate the delta phi and p_t associated hadron histograms for each p_t bin of the gamma!
			fHistDpGh[identifier][i] = new TH1F(Form("fHistDpGh_%0d_%0d",identifier,i),Form("fHistDpGh_%0d_%0d",identifier,i), nbins[2], min[2], max[2]);
			fHistDpGh[identifier][i]->GetXaxis()->SetTitle(Form("#Delta #phi^{#gamma-h} %0.1f<p_{T}^{#gamma}<%0.1f",BinValStart,BinValStart+BinWidth));
			fHistDpGh[identifier][i]->GetYaxis()->SetTitle(Form("dN^{#gamma-h}/#Delta #phi^{#gamma-h} [counts/%0.1f^{#circ}]",fHistDpGh[identifier][i]->GetBinWidth(1)));
			if(i!=0)fOutputList3->Add(fHistDpGh[identifier][i]);
		}
		fOutput->Add(fHistptAssHadron[identifier][0]);
		fOutput->Add(fHistDpGh[identifier][0]);
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
	}
	//time information??



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


	//--check here some basic properties of the event
	/*
	if (fCentrality > fBCent->GetXmax() || fCentrality < fBCent->GetXmin()) {
	    if (fVerbosity > 1)
	      AliInfo(Form("Event REJECTED (centrality out of range). fCentrality = %.1f", fCentrality));
	    return kFALSE;
	  }
	  if (fZVertex > fBZvtx->GetXmax() || fZVertex < fBZvtx->GetXmin()) {
	    if (fVerbosity > 1)
	      AliInfo(Form("Event REJECTED (z_vertex out of range). fZVertex = %.1f", fZVertex));
	    return kFALSE;
	  }
	*/
	//??	MC  // see if event is selected
	//??MC	  UInt_t trig = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();

	//for same event only analyse events when there is a cluster inside
	if (!fCaloClusters && fSameEventAnalysis==1)return kFALSE;

	return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskGammaHadron::RetrieveEventObjects()
{
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::RuRetrieveEventObjectsn()"<<endl;
	// Retrieve event objects.
	if (!AliAnalysisTaskEmcal::RetrieveEventObjects())
	{
		return kFALSE;
	}
	//Get additional centrality selections
	//the standard centrality is given by "fCent" and set in TaskEmcal to "V0M" by default
	//it can, however, be changed by the method SetCentralityEstimator(const char *c)
	//Q:: how do I know that this isn't done somewhere?
	if (!fCentMethodAlt.IsNull())
	{
		if (fBeamType == kAA || fBeamType == kpA )//what the hack is this?
		{
			AliCentrality *aliCent = InputEvent()->GetCentrality();
			if (aliCent)
			{
				fCentAlt = aliCent->GetCentralityPercentile(fCentMethodAlt);
			}
		}
	}
	return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskGammaHadron::FillHistograms()
{
	//This function is called in AliAnalysisTaskEmcal::UserExec.
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

	//--Get pool containing tracks from other events like this one
	Double_t ZVertex = fVertex[2];
    //where does this come from???
	AliEventPool* pool = 0x0;
	pool = fPoolMgr->GetEventPool(fCent, ZVertex);
	if (!pool)
	{
		AliWarning(Form("No pool found. Centrality %f, ZVertex %f",fCent, ZVertex));
		return kFALSE;
	}
	pool->SetTargetEvents(5); // *New* set the minimum number of events to mix

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//    Mixed event section
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	if(pool->IsReady())// || pool->NTracksInPool() > fTrackDepth || pool->GetCurrentNEvents() >= 5)
	{
        // get number of current events in pool
		Int_t nMix = pool->GetCurrentNEvents();

		//cout<<"number of events in pool: "<<nMix<<endl;
		for (Int_t jMix=0; jMix<nMix; jMix++)
		{
			TObjArray* bgTracks=0x0;
			bgTracks = pool->GetEvent(jMix);

			if(!bgTracks)
			{
				cout<<"could not retrieve TObjArray from EventPool!"<<endl;
			}
			//--Loop over clusters and fill histograms
			if(fGammaOrPi0==0) CorrelateClusterAndTrack(0,bgTracks,1,1.0/nMix);//correlate with mixed event
			else               CorrelatePi0AndTrack(0,bgTracks,1,1.0/nMix);    //correlate with mixed event
		}
	}
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//    Same event section
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	AliParticleContainer* tracks =0x0;
	tracks   = GetParticleContainer(0);
	//--Loop over clusters and fill histograms
	if(fGammaOrPi0==0) CorrelateClusterAndTrack(tracks,0,0,1);//correlate with same event
	else               CorrelatePi0AndTrack(tracks,0,0,1);    //correlate with same event

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//    Update the pool
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	TObjArray* tracksClone=0x0;
	tracksClone = CloneToCreateTObjArray(tracks);
   	//if(!tracksClone)cout<<"**No arrray!!!!!"<<endl;
   	if(tracksClone)
   	{
   		pool->UpdatePool(tracksClone);
   	}

	return kTRUE;
}
//________________________________________________________________________
void AliAnalysisTaskGammaHadron::InitEventMixer()
{
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::InitEventMixer()"<<endl;
	//--The effective pool size in events is set by trackDepth, so more
	//--low-mult events are required to maintain the threshold than
	//--high-mult events. Centrality pools are indep. of data histogram
	//--binning, no need to match.

	// Centrality pools
	Int_t nCentBins=fMixBCent->GetNbins();
	Double_t centBins[nCentBins+1];
	centBins[0] = fMixBCent->GetBinLowEdge(1);
	for(Int_t i=1; i<=nCentBins; i++)
	{
		centBins[i] = fMixBCent->GetBinUpEdge(i);
	}

	// Z-vertex pools
	Int_t nZvtxBins=fMixBZvtx->GetNbins();
	Double_t zvtxbin[nZvtxBins+1];
	zvtxbin[0] = fMixBZvtx->GetBinLowEdge(1);
	for(Int_t i=1; i<=nZvtxBins; i++)
	{
		zvtxbin[i] = fMixBZvtx->GetBinUpEdge(i);
	}
	fPoolMgr = new AliEventPoolManager(fPoolSize, fTrackDepth, nCentBins, centBins, nZvtxBins, zvtxbin);

	//just for explaining the parameters
	//  (Int_t maxEvts, Int_t minNTracks, Int_t nMultBins, Double_t* multbins, Int_t nZvtxBins, Double_t* zvtxbins, Int_t nPsiBins, Double_t* psibins)
	//  fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins, centralityBins, nZvtxBins, zvtxbin);

}
//________________________________________________________________________
TObjArray* AliAnalysisTaskGammaHadron::CloneToCreateTObjArray(AliParticleContainer* tracks)
{
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::CloneToCreateTObjArray()"<<endl;
	//--clones a track list
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
	Int_t nAccClusters = 0;
	Double_t EffWeight_Gamma;
	Double_t EffWeight_Hadron;
	Double_t Weight;    //weight to normalize mixed and same event distributions individually

	AliVCluster* cluster = 0;
	AliVCluster* cluster2= 0;
	AliVParticle* trackNULL=0;

	//...........................................
	//..do a small loop to count the triggers in this event
	//..for same events normalize per triggers (gammas/pi0)
	if(SameMix==0)
	{
		for(Int_t NoCluster1 = 0; NoCluster1 < NoOfClustersInEvent; NoCluster1++ )
		{
			cluster=(AliVCluster*) clusters->GetAcceptCluster(NoCluster1); //->GetCluster(NoCluster1);
			if(!cluster || !AccClusterForAna(cluster))continue;
			TLorentzVector CaloClusterVec;
			clusters->GetMomentum(CaloClusterVec,cluster);
			FillQAHisograms(0,CaloClusterVec,trackNULL,0,0);

			nAccClusters++;
		}
		if(tracks)  //just for current debugging move back to previous stage later
		{
			Int_t NoOfTracksInEvent =tracks->GetNParticles();
			AliVParticle* track=0;

			for(Int_t NoTrack = 0; NoTrack < NoOfTracksInEvent; NoTrack++)
			{
				track = (AliVParticle*)tracks->GetAcceptParticle(NoTrack);
				if(!track)continue; //check if the track is a good track

				fHistptAssHadron[0][0]->Fill(track->Pt());
			}
		}
		if(fUsePerTrigWeight==1)
		{
			if(nAccClusters>0)	Weight=1.0/(1.0*nAccClusters);
			else             	Weight=0;
		}
		else
		{
			//..don't use per gamma yields but absolute numbers
			Weight=1;
		}

	}
	//..for mixed events normalize per events in pool
	if(SameMix==1)
	{
		Weight=InputWeight;
	}
	//...........................................
	//..run the real loop for filling the histograms
	Int_t GammaCounter=0;
	for(Int_t NoCluster1 = 0; NoCluster1 < NoOfClustersInEvent; NoCluster1++ )
	{
		cluster=(AliVCluster*) clusters->GetAcceptCluster(NoCluster1); //->GetCluster(NoCluster1);
		if(!cluster || !AccClusterForAna(cluster))continue; //check if the cluster is a good cluster
		//clusters->GetLeadingCluster("e");

		TLorentzVector CaloClusterVec;
	//old framework	cluster->GetMomentum(CaloClusterVec, fVertex);
		clusters->GetMomentum(CaloClusterVec, cluster);
		EffWeight_Gamma=GetEff(CaloClusterVec);

		fHistNoClusPt->Fill(CaloClusterVec.Pt()); //the .pt only works for gammas (E=M) for other particle this is wrong
		//cout<<"Cluster number: "<<NoCluster1<<", pT = "<<CaloClusterVec.Pt()<<endl;
		//...........................................
		//..combine gammas with same event tracks
		GammaCounter++;
		if(SameMix==0)
		{
			//cout<<"SameMix==0"<<endl;
			if(!tracks)  return 0;
			Int_t NoOfTracksInEvent =tracks->GetNParticles();
			AliVParticle* track=0;

			for(Int_t NoTrack = 0; NoTrack < NoOfTracksInEvent; NoTrack++)
			{
				track = (AliVParticle*)tracks->GetAcceptParticle(NoTrack);
				if(!track)continue; //check if the track is a good track
				//cout<<"..Track number: "<<NoTrack<<", pT = "<<track->Pt()<<endl;

				//EffWeight_Hadron=GetEff(<TLorentzVector>track);
				FillGhHisograms(1,CaloClusterVec,track,2,0,-360,Weight);
				if(GammaCounter==1)FillQAHisograms(1,CaloClusterVec,track,2,0); // fill only once
                /*fill only for first hadron*/FillQAHisograms(2,CaloClusterVec,track,2,0); //fill for each gamma?? makes sense?
				if(CaloClusterVec.Eta()>0)FillQAHisograms(3,CaloClusterVec,track,2,0); //fill for each gamma?? makes sense?
				if(CaloClusterVec.Eta()<0)FillQAHisograms(4,CaloClusterVec,track,2,0); //fill for each gamma?? makes sense?
			}
		}
		//...........................................
		//..combine gammas with mixed event tracks
		if(SameMix==1)
		{
			Int_t Nbgtrks = bgTracksArray->GetEntries();
			for(Int_t ibg=0; ibg<Nbgtrks; ibg++)
			{
                AliPicoTrack* track = static_cast<AliPicoTrack*>(bgTracksArray->At(ibg));
				if(!track) continue;

				//EffWeight_Hadron=GetEff(<TLorentzVector>track);
				FillGhHisograms(0,CaloClusterVec,track,2,0,-360,Weight);
			}
		}
		//...........................................
		//..double cluster loop for testing
		for( Int_t NoCluster2 = 0; NoCluster2 < NoOfClustersInEvent; NoCluster2++ )
		{
			if(NoCluster1!=NoCluster2)
			{
				cluster2=(AliVCluster*) clusters->GetAcceptCluster(NoCluster2);
				if(!cluster2 || !AccClusterForAna(cluster2))continue; //check if the cluster is a good cluster

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
	return nAccClusters;
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
	// **Tyler you have to see whether you need to count here pi0 inside of your mass window
	// **but I also havent implemented it for gammas yet :-D
	if(SameMix==0)
	{
		for(Int_t NoCluster1 = 0; NoCluster1 < NoOfClustersInEvent; NoCluster1++ )
		{
			cluster=(AliVCluster*) clusters->GetAcceptCluster(NoCluster1); //->GetCluster(NoCluster1);
			if(!cluster || !AccClusterForAna(cluster))continue; //check if the cluster is a good cluster
			//clusters->GetLeadingCluster("e");

			TLorentzVector CaloClusterVec;
//old framework			cluster->GetMomentum(CaloClusterVec, fVertex);
			clusters->GetMomentum(CaloClusterVec, cluster);
			//acc if pi0 candidate
			nAccClusters++;

			for( Int_t NoCluster2 = 0; NoCluster2 < NoOfClustersInEvent; NoCluster2++ )
			{
				if(NoCluster1!=NoCluster2)
				{
					cluster2=(AliVCluster*) clusters->GetAcceptCluster(NoCluster2);
					if(!cluster2 || !AccClusterForAna(cluster2))continue; //check if the cluster is a good cluster

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
		if(fUsePerTrigWeight==1)
		{
			if(nAccPi0Clusters>0)Weight=2.0/(1.0*nAccPi0Clusters); //2 to avoid double counting.
			else             	 Weight=0;
		}
		else
		{
			//..don't use per gamma yields but absolute numbers
			Weight=1;
		}
	}
	//..for mixed events normalize per events in pool
	if(SameMix==1)
	{
		Weight=InputWeight;
	}
	//if(nAccClusters>0) cout<<"Pi0: "<<nAccPi0Clusters<<" all accepted clusters: "<<nAccClusters<<endl;
	//...........................................
	//run the real loop for filling the histograms
	for( Int_t NoCluster1 = 0; NoCluster1 < NoOfClustersInEvent; NoCluster1++ )
	{

		cluster=(AliVCluster*) clusters->GetAcceptCluster(NoCluster1); //->GetCluster(NoCluster1);
		if(!cluster || !AccClusterForAna(cluster))continue; //check if the cluster is a good cluster
		//clusters->GetLeadingCluster("e");

		TLorentzVector CaloClusterVec;
	//old framework	cluster->GetMomentum(CaloClusterVec, fVertex);
		clusters->GetMomentum(CaloClusterVec,cluster);

		for( Int_t NoCluster2 = 0; NoCluster2 < NoOfClustersInEvent; NoCluster2++ )
		{
			if(NoCluster1!=NoCluster2)
			{
				cluster2=(AliVCluster*) clusters->GetAcceptCluster(NoCluster2);
				if(!cluster2 || !AccClusterForAna(cluster2))continue; //check if the cluster is a good cluster

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
					if(SameMix==0)
					{
						//cout<<"SameMix==0"<<endl;
						if(!tracks)  return 0;
						Int_t NoOfTracksInEvent =tracks->GetNParticles();
						AliVParticle* track=0;

						for(Int_t NoTrack = 0; NoTrack < NoOfTracksInEvent; NoTrack++)
						{
							track = (AliVParticle*)tracks->GetAcceptParticle(NoTrack);
							if(!track)continue; //check if the track is a good track

							//..fill here eventually a pi0 four-vector instead of CaloClusterVec
							//EffWeight_Hadron=GetEff(TLorentzVector)track);
							FillGhHisograms(1,CaloClusterVecpi0,track,2,0,-360,Weight);
						}
					}
					//...........................................
					//..combine gammas with mixed event tracks
					if(SameMix==1)
					{
						Int_t Nbgtrks = bgTracksArray->GetEntries();
						for(Int_t ibg=0; ibg<Nbgtrks; ibg++)
						{
							AliPicoTrack* track = static_cast<AliPicoTrack*>(bgTracksArray->At(ibg));
							if(!track) continue;

							//**fill here eventually a pi0 four-vector instead of CaloClusterVec
							//EffWeight_Hadron=GetEff((TLorentzVector)track);
							FillGhHisograms(0,CaloClusterVecpi0,track,2,0,-360,Weight);
						}
					}
				}
			}
		}
	}
	return nAccPi0Clusters/2;
}
//________________________________________________________________________
void AliAnalysisTaskGammaHadron::FillGhHisograms(Int_t identifier,TLorentzVector ClusterVec,AliVParticle* TrackVec, Double_t ClusterEcut, Double_t TrackPcut, Double_t Anglecut, Double_t Weight)
{
	if(fDebug==1)cout<<"Inside of: AliAnalysisTaskGammaHadron::FillGhHisograms()"<<endl;

	//--This function fills several histograms under different cut conditions.
	//--it is run within a cluster{ track{}} loop to get all combinations.

	//--A word to the weight - for mixed events it devides by the number of events in the current pool 1/nEvents
	//--                     - for same events you devide by the number of triggers
	//--                     - for both you have to take into account the efficiency of your correlated pair
	Double_t deltaEta   = ClusterVec.Eta()-TrackVec->Eta();
	Double_t deltaPhi   = DeltaPhi(ClusterVec,TrackVec);
    Double_t G_PT_Value = ClusterVec.Pt();
    //Double_t ZT_Value   = TMath::Cos(deltaPhi*1/fRtoD)*TrackVec->P()/ClusterVec.P(); //   TrackVec->Pt()/G_PT_Value;
    Double_t ZT_Value   = TrackVec->Pt()/G_PT_Value; //   TrackVec->Pt()/G_PT_Value;
	//--Careful here: usually this is done for an opening angle (hadron-jet axis) of less than 90¡. Due to
    //--resolution momentum smearing (our guess - check that!) there are particles appearing at angles greater than 90¡
    Double_t XI_Value;
    if(ZT_Value>0)
    {
	  XI_Value   = TMath::Log(1.0/ZT_Value);
    }
    /*else
    {
    	  XI_Value   = TMath::Log(1.0/(-1.0)*ZT_Value);
    }*/


	if(G_PT_Value>=ClusterEcut && TrackVec->Pt()>=TrackPcut && deltaPhi>=Anglecut)
	{
		fHistNoClusPtH[identifier]     ->Fill(G_PT_Value,Weight);    //the .pt only works for gammas (E=M) for other particle this is wrong
///*comment for the moment!!*/		fHistptAssHadron[identifier][0]->Fill(TrackVec->Pt(),Weight);
		fHistDpGh[identifier][0]       ->Fill(deltaPhi,Weight);

		//--Fill 2D Histograms for certain event conditions
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
		//--fill histograms for different ranges of p_t^{g}
		for(Int_t i=1;i<fHistNoClusPtH[identifier]->GetNbinsX()+1;i++)  // be careful hard coded from max[0] value -- need better variable // fHistNoClusPtH->GetLast bin etc??
		{
			//--look in each p_T bin of the gamma (fHistNoClusPtH histogram) how the pt
			//--distribution of the associated hadron looks like.
			Double_t BinWidth    = fHistNoClusPtH[identifier]->GetBinWidth(i);
			Double_t BinValStart = fHistNoClusPtH[identifier]->GetXaxis()->GetBinCenter(i)-BinWidth/2.0;

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
void AliAnalysisTaskGammaHadron::FillQAHisograms(Int_t identifier,TLorentzVector ClusterVec,AliVParticle* TrackVec, Double_t ClusterEcut, Double_t TrackPcut)
{
	/*should do similar test here*/fHistDEtaDPhiGammaQA[identifier] ->Fill(ClusterVec.Phi()*fRtoD,ClusterVec.Eta());
	if(TrackVec)                   fHistDEtaDPhiTrackQA[identifier] ->Fill(TrackVec->Phi()*fRtoD,TrackVec->Eta());
}
//________________________________________________________________________
Bool_t AliAnalysisTaskGammaHadron::AccClusterForAna(AliVCluster* cluster)
{
	//--Accepts clusters if certain conditions are fulfilled
	Bool_t Accepted=0; //By default rejceted

	//double check these cuts carefully with the experts
	if(cluster->GetNCells()>1)
	{
		//--Now accept the cluster as a good candidate for your analysis
		Accepted=1;
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
    //--shift the second peak over the fist peak: \--Æ--/   --> -Æ--
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
