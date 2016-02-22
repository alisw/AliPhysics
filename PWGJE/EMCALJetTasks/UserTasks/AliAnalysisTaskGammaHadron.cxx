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
#include "AliEmcalJet.h"
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
AliAnalysisTaskGammaHadron::AliAnalysisTaskGammaHadron() :
AliAnalysisTaskEmcalJet("AliAnalysisTaskGammaHadron", kTRUE),
fCellEnergyCut(0.05),
fMaxCellsInCluster(50),

fParticleLevel(kFALSE),
fIsMC(kFALSE),
fCent_alt(),
fPoolMgr(0x0),
fHistEff_Gamma(0x0),
fHistEff_Hadron(0x0),


fOutputList1(),
fOutputList2(),
fOutputList3(),

fHistNoClus_pt_Trigger(0),
fHistNoClus_pt(0),
fHistNoClus_ptH(0),
fHistpi0(0),
fHist_dEta_dPhi(0),
fHist_dEta_dPhi_low(0),
fHist_dEta_dPhi_med(0),
fHist_dEta_dPhi_high(0),
fHPoolReady(0x0)

{
	// Default constructor.
	fAODfilterBits[0] = 0;
	fAODfilterBits[1] = 0;

	fHistpt_assHadron[0]= 0;
	fHistpt_assHadron[1]= 0;
	fHistpt_assHadron[2]= 0;

	fHist_DP_gh[0]= 0;
	fHist_DP_gh[1]= 0;
	fHist_DP_gh[2]= 0;

	fRtoD=180.0/TMath::Pi();


	//Set some default values.
	//if desired one can add a set function to
	//set these values in the add task function
	Int_t fCentBinSize= 10;//<<<<< to be set in a setter function

	Int_t NcentBins = 100;
	Double_t mult = 1.0;
	if(fCentBinSize==1)
	{
		NcentBins = 100;
		mult = 1.0;
	}
	else if(fCentBinSize==2)
	{
		NcentBins = 50;
		mult = 2.0;
	}
	else if(fCentBinSize==5)
	{
		NcentBins = 20;
		mult = 5.0;
	}
	else if(fCentBinSize==10)
	{
		NcentBins = 10;
		mult = 10.0;
	}

	Double_t centmix[NcentBins+1];
	for(Int_t ic=0; ic<NcentBins+1; ic++)
	{
		centmix[ic]=mult*ic;
	}
	fMixBCent = new TAxis(NcentBins,centmix);
    static const Int_t NvertBins=8;
    Double_t zvtxmix[NvertBins+1] = {-10,-6,-4,-2,0,2,4,6,10};
    fMixBZvtx = new TAxis(NvertBins,zvtxmix);
    fCentMethod_alt = "V0M";
    fTrackDepth     = 100;     //Hanseul sets it to 100! Q:: is this good?
    fPoolSize       = 200;     //hanseuls default value Q:: is this correct?

    SetMakeGeneralHistograms(kTRUE);  //What is this?? In jet exe?
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
	//cout<<"Inside of: AliAnalysisTaskGammaHadron::UserCreateOutputObjects()"<<endl;

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
	// Create histograms
	AliAnalysisTaskEmcalJet::UserCreateOutputObjects();  //Q:: what is created here?

	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//
	//   Define sublists/folders for a better organisation of the figures
	//
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	fOutputList1 = new TList();
	fOutputList1->SetOwner();
	fOutputList1->SetName("p_{T} distributions of the gamma");
	fOutputList2 = new TList();
	fOutputList2->SetOwner();
	fOutputList2->SetName("p_{T} distributions of the associated hadrons, for a given p_{T}^{gamma}");
	fOutputList3 = new TList();
	fOutputList3->SetOwner();
	fOutputList3->SetName("Delta phi^{g-h} for a given p_{T}^{gamma}");

	//common bins for the histograms
	Int_t nbins[5] = {0};
	Double_t min[5] = {0};
	Double_t max[5] = {0};

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
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//
	//   Create Histograms
	//
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	//.................................
	//Initiate histograms for given p_t gamma bins!
	//Later you can summ the histograms according to the expected statistic!
	static const Int_t NoOfDPhistos=31;// =  nbins[0];

	fHistpt_assHadron[0]    = new TH1*[NoOfDPhistos]; // make a p_t histogram of the associated hadron for each p_T bin of the gamma
	fHist_DP_gh[0]          = new TH1*[NoOfDPhistos]; // make a p_t histogram of the associated hadron for each p_T bin of the gamma
	fHistpt_assHadron[1]    = new TH1*[NoOfDPhistos]; // make a p_t histogram of the associated hadron for each p_T bin of the gamma
	fHist_DP_gh[1]          = new TH1*[NoOfDPhistos]; // make a p_t histogram of the associated hadron for each p_T bin of the gamma
	fHistpt_assHadron[2]    = new TH1*[NoOfDPhistos]; // make a p_t histogram of the associated hadron for each p_T bin of the gamma
	fHist_DP_gh[2]          = new TH1*[NoOfDPhistos]; // make a p_t histogram of the associated hadron for each p_T bin of the gamma


	fHistNoClus_ptH          = new TH1*[3];
	fHist_dEta_dPhi          = new TH2*[3];
	fHist_dEta_dPhi_low      = new TH2*[3];
	fHist_dEta_dPhi_med      = new TH2*[3];
	fHist_dEta_dPhi_high     = new TH2*[3];
	//.................................
	// p_T^{Cluster} distribution under different conditions

	// all clusters
	fHistNoClus_pt = new TH1F(Form("fHistNoClus_pt_%0d",1),Form("fHistNoClus_pt_%0d",1), nbins[0], min[0], max[0]);
	fHistNoClus_pt->GetXaxis()->SetTitle("p_{T}^{Calo Cluster}");
	fHistNoClus_pt->GetYaxis()->SetTitle(Form("No. of Clusters [counts/%0.1f GeV/c]",fHistNoClus_pt->GetBinWidth(0)));
	fOutputList1->Add(fHistNoClus_pt);

	//trigger histogram only for internal purposes - no nead to save
	//fHistNoClus_pt_Trigger = new TH1F((const TH1F)fHistNoClus_pt);

	//test!!
	fHistpi0 = new TH1F(Form("fHistpi0_%0d",1),Form("fHistpi0_%0d",1), 500, 0, 0.5);
	fHistpi0->GetXaxis()->SetTitle("M_{#gamma#gamma}");
	fHistpi0->GetYaxis()->SetTitle("Entries");
	fOutput->Add(fHistpi0);


	//by the identifier different histograms can be filled under different cut conditions
	//Can eg. later be modified to contain certain delta phi or centrality bins
	for(Int_t identifier=0;identifier<3;identifier++)
	{

		//clusters p_T if there was a hadron present
		fHistNoClus_ptH[identifier] = new TH1F(Form("fHistNoClus_ptH_%0d",identifier),Form("fHistNoClus_ptH_%0d",identifier), nbins[0], min[0], max[0]);
		fHistNoClus_ptH[identifier]->GetXaxis()->SetTitle("p_{T}^{Calo Cluster}");
		fHistNoClus_ptH[identifier]->GetYaxis()->SetTitle(Form("No. of Clus. with h [counts/%0.1f GeV/c]",fHistNoClus_ptH[identifier]->GetBinWidth(0)));
		fOutputList1->Add(fHistNoClus_ptH[identifier]);

		//.................................
		// two dimensional delta eta delta phi distributions
		fHist_dEta_dPhi[identifier] = new TH2F(Form("fHist_dEta_dPhi_%0d",identifier),Form("fHist_dEta_dPhi_%0d",identifier),nbins[2],min[2],max[2],nbins[3],min[3],max[3]);
		fHist_dEta_dPhi[identifier]->GetXaxis()->SetTitle(Form("#Delta #phi^{#gamma-h} %0.1d<p_{T}^{#gamma}<%0.1d",0,100));
		fHist_dEta_dPhi[identifier]->GetYaxis()->SetTitle("#Delta #eta^{#gamma-h}");
		fOutput->Add(fHist_dEta_dPhi[identifier]);

		fHist_dEta_dPhi_low[identifier] = new TH2F(Form("fHist_dEta_dPhi_low_%0d",identifier),Form("fHist_dEta_dPhi_low_%0d",identifier),nbins[2],min[2],max[2],nbins[3],min[3],max[3]);
		fHist_dEta_dPhi_low[identifier]->GetXaxis()->SetTitle(Form("#Delta #phi^{#gamma-h} %0.1d<p_{T}^{#gamma}<%0.1d",0,5));
		fHist_dEta_dPhi_low[identifier]->GetYaxis()->SetTitle("#Delta #eta^{#gamma-h}");
		fOutput->Add(fHist_dEta_dPhi_low[identifier]);

		fHist_dEta_dPhi_med[identifier] = new TH2F(Form("fHist_dEta_dPhi_med_%0d",identifier),Form("fHist_dEta_dPhi_med_%0d",identifier),nbins[2],min[2],max[2],nbins[3],min[3],max[3]);
		fHist_dEta_dPhi_med[identifier]->GetXaxis()->SetTitle(Form("#Delta #phi^{#gamma-h} %0.1d<p_{T}^{#gamma}<%0.1d",5,10));
		fHist_dEta_dPhi_med[identifier]->GetYaxis()->SetTitle("#Delta #eta^{#gamma-h}");
		fOutput->Add(fHist_dEta_dPhi_med[identifier]);

		fHist_dEta_dPhi_high[identifier] = new TH2F(Form("fHist_dEta_dPhi_high_%0d",identifier),Form("fHist_dEta_dPhi_high_%0d",identifier),nbins[2],min[2],max[2],nbins[3],min[3],max[3]);
		fHist_dEta_dPhi_high[identifier]->GetXaxis()->SetTitle(Form("#Delta #phi^{#gamma-h} %0.1d<p_{T}^{#gamma}<%0.1d",10,30));
		fHist_dEta_dPhi_high[identifier]->GetYaxis()->SetTitle("#Delta #eta^{#gamma-h}");
		fOutput->Add(fHist_dEta_dPhi_high[identifier]);

		//l
		for(Int_t i=0; i<NoOfDPhistos+1; i++)
		{
			//check whether the max is the
			//same as the histogram size
			Double_t BinWidth,BinValStart;

			if(i==0)
			{
				//these are histograms over the full p_T gamma range (no binnings)
				BinWidth    = fHistNoClus_ptH[0]->GetBinWidth(1);
				BinValStart = fHistNoClus_ptH[0]->GetXaxis()->GetBinCenter(1)-BinWidth/2.0;
				BinWidth    =100; //to get the full length of the histogram
			}
			else
			{
				//these are histograms for a certain p_T of the gamma
				BinWidth    = fHistNoClus_ptH[0]->GetBinWidth(i);
				BinValStart = fHistNoClus_ptH[0]->GetXaxis()->GetBinCenter(i)-BinWidth/2.0;
				//cout<<"BinStart: "<<fHistNoClus_ptH->GetBinCenter(i)<<", width "<<BinWidth<<endl;
			}

			//.................................
			//Initiate the delta phi and p_t associated hadron histograms for each p_t bin of the gamma!
			fHistpt_assHadron[identifier][i] = new TH1F(Form("fHistAssHadron_pt_%0d_%0d",identifier,i),Form("fHistAssHadron_pt_%0d_%0d",identifier,i), nbins[1], min[1], max[1]);
			fHistpt_assHadron[identifier][i]->GetXaxis()->SetTitle(Form("p_{T}^{assoc. h} %0.1f<p_{T}^{#gamma}<%0.1f",BinValStart,BinValStart+BinWidth));
			fHistpt_assHadron[identifier][i]->GetYaxis()->SetTitle(Form("dN^{assoc. h}/dp_{T}^{h} [counts/%0.1f GeV/c]",fHistpt_assHadron[identifier][i]->GetBinWidth(1)));
			if(i!=0)fOutputList2->Add(fHistpt_assHadron[identifier][i]);

			//.................................
			//Initiate the delta phi and p_t associated hadron histograms for each p_t bin of the gamma!
			fHist_DP_gh[identifier][i] = new TH1F(Form("fHist_DP_gh_%0d_%0d",identifier,i),Form("fHist_DP_gh_%0d_%0d",identifier,i), nbins[2], min[2], max[2]);
			fHist_DP_gh[identifier][i]->GetXaxis()->SetTitle(Form("#Delta #phi^{#gamma-h} %0.1f<p_{T}^{#gamma}<%0.1f",BinValStart,BinValStart+BinWidth));
			fHist_DP_gh[identifier][i]->GetYaxis()->SetTitle(Form("dN^{#gamma-h}/#Delta #phi^{#gamma-h} [counts/%0.1f^{#circ}]",fHist_DP_gh[identifier][i]->GetBinWidth(1)));
			if(i!=0)fOutputList3->Add(fHist_DP_gh[identifier][i]);
		}
		fOutput->Add(fHistpt_assHadron[identifier][0]);
		fOutput->Add(fHist_DP_gh[identifier][0]);
	}

	fOutput->Add(fOutputList1);
	fOutput->Add(fOutputList2);
	fOutput->Add(fOutputList3);

	PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}
//________________________________________________________________________
void AliAnalysisTaskGammaHadron::ExecOnce()
{

	AliAnalysisTaskEmcal::ExecOnce();

	/*
	AliAnalysisTaskEmcalJet::ExecOnce();
	This function calls  AliAnalysisTaskEmcal::ExecOnce();
	and gets the objects for rho (fRhoName), localrho (fLocalRhoName), AliJetContainer GetJetContainer(0)
	 */
}
//________________________________________________________________________
Bool_t AliAnalysisTaskGammaHadron::Run()
{
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

	if (!fCaloClusters)return kFALSE;

	return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskGammaHadron::RetrieveEventObjects()
{
	// Retrieve event objects.
	if (!AliAnalysisTaskEmcalJet::RetrieveEventObjects())
	{
		return kFALSE;
	}
	//Get additional centrality selections
	//the standard centrality is given by "fCent" and set in TaskEmcal to "V0M" by default
	//it can, however, be changed by the method SetCentralityEstimator(const char *c)
	//Q:: how do I know that this isn't done somewhere?
	if (!fCentMethod_alt.IsNull())
	{
		if (fBeamType == kAA || fBeamType == kpA )//what the hack is this?
		{
			AliCentrality *aliCent = InputEvent()->GetCentrality();
			if (aliCent)
			{
				fCent_alt = aliCent->GetCentralityPercentile(fCentMethod_alt);
			}
		}
	}
	return kTRUE;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskGammaHadron::FillHistograms()
{
	//cout<<"Inside of: AliAnalysisTaskGammaHadron::FillHistograms()"<<endl;

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

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//    Mixed event section
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	if(pool->IsReady())// || pool->NTracksInPool() > fTrackDepth || pool->GetCurrentNEvents() >= 5)
	{
        // get number of current events in pool
		Int_t nMix = pool->GetCurrentNEvents();

		for (Int_t jMix=0; jMix<nMix; jMix++)
		{
			TObjArray* bgTracks=0x0;
			bgTracks = pool->GetEvent(jMix);

			if(!bgTracks)
			{
				cout<<"could not retrieve TObjArray from EventPool!"<<endl;
			}
			//--Loop over clusters and fill histograms
			CorrelateClusterAndTrack(0,bgTracks,1,1.0/nMix);//correlate with mixed event
		}
	}
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//    Same event section
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	AliParticleContainer* tracks =0x0;
	tracks   = GetParticleContainer(0);
	//--Loop over clusters and fill histograms
	CorrelateClusterAndTrack(tracks,0,0,1);//correlate with same event

	TObjArray* tracksClone=0x0;
	tracksClone = CloneToCreateTObjArray(tracks);
   	if(!tracksClone)cout<<"**No arrray!!!!!"<<endl;
   	if(tracksClone)
   	{
   		pool->UpdatePool(tracksClone);
   	}
	return kTRUE;
}
//________________________________________________________________________
void AliAnalysisTaskGammaHadron::InitEventMixer()
{
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
	//--clones a track list
	if(!tracks)                            return 0;
	if(tracks->GetNAcceptedParticles()==0) return 0;

	//cout<<"added list with "<<tracks->GetNAcceptedParticles()<<" particles"<<endl;
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

	//cout<<"New list contains "<<tracksClone->GetEntries()<<endl;
	if(tracksClone->GetEntries()!=tracks->GetNAcceptedParticles())cout<<"!!!!!!!major error!!!!"<<endl;
	return tracksClone;
}
//________________________________________________________________________
Int_t AliAnalysisTaskGammaHadron::CorrelateClusterAndTrack(AliParticleContainer* tracks,TObjArray* bgTracksArray,Bool_t SameMix, Double_t Weight)
{
	//cout<<"Inside of: AliAnalysisTaskGammaHadron::CorrelateClusterAndTrack()"<<endl;

	//...........................................
	//--Do cluster loop.
	AliClusterContainer* clusters  = GetClusterContainer(0);  //how do I know which cells are selected
	if (!clusters) return 0;
	Int_t NoOfClustersInEvent =clusters->GetNClusters();
	Int_t nAccClusters = 0;

	AliVCluster* cluster = 0;
	AliVCluster* cluster2= 0;

	//...........................................
	//do a small loop to count the triggers in this event
	for(Int_t NoCluster1 = 0; NoCluster1 < NoOfClustersInEvent; NoCluster1++ )
	{
		cluster=(AliVCluster*) clusters->GetAcceptCluster(NoCluster1); //->GetCluster(NoCluster1);
		if(!cluster || !AccClusterForAna(cluster))continue; //check if the cluster is a good cluster
		//clusters->GetLeadingCluster("e");

		TLorentzVector CaloClusterVec;
		cluster->GetMomentum(CaloClusterVec, fVertex);
		//fHistNoClus_pt_Trigger->Fill(CaloClusterVec.Pt()); //the .pt only works for gammas (E=M) for other particle this is wrong
	}

	//...........................................
	//run the real loop for filling the histograms
	for( Int_t NoCluster1 = 0; NoCluster1 < NoOfClustersInEvent; NoCluster1++ )
	{
		cluster=(AliVCluster*) clusters->GetAcceptCluster(NoCluster1); //->GetCluster(NoCluster1);
		if(!cluster || !AccClusterForAna(cluster))continue; //check if the cluster is a good cluster
		//clusters->GetLeadingCluster("e");

		TLorentzVector CaloClusterVec;
		cluster->GetMomentum(CaloClusterVec, fVertex);

		fHistNoClus_pt->Fill(CaloClusterVec.Pt()); //the .pt only works for gammas (E=M) for other particle this is wrong

		//...........................................
		//--combine gammas with same event tracks
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
				Fill_GH_Hisograms(1,CaloClusterVec,track,2,0,-360,Weight);
			}
		}
		//...........................................
		//--combine gammas with mixed event tracks
		if(SameMix==1)
		{
			Int_t Nbgtrks = bgTracksArray->GetEntries();
			for(Int_t ibg=0; ibg<Nbgtrks; ibg++)
			{
                AliPicoTrack* track = static_cast<AliPicoTrack*>(bgTracksArray->At(ibg));
				if(!track) continue;

				Fill_GH_Hisograms(0,CaloClusterVec,track,2,0,-360,Weight);
			}
		}

		nAccClusters++;
		//...........................................
		//--double cluster loop for testing
		for( Int_t NoCluster2 = 0; NoCluster2 < NoOfClustersInEvent; NoCluster2++ )
		{
			if(NoCluster1!=NoCluster2)
			{
				cluster2=(AliVCluster*) clusters->GetAcceptCluster(NoCluster2);
				if(!cluster2 || !AccClusterForAna(cluster2))continue; //check if the cluster is a good cluster

				TLorentzVector CaloClusterVec2;
				TLorentzVector CaloClusterVecpi0;
				cluster2->GetMomentum(CaloClusterVec2, fVertex);
				if(cluster2->E()>2 && cluster->E()>2)
				{
					CaloClusterVecpi0=CaloClusterVec+CaloClusterVec2;
					fHistpi0->Fill(CaloClusterVecpi0.M());
				}
			}
		}
	}
	return nAccClusters;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskGammaHadron::AccClusterForAna(AliVCluster* cluster)
{
	//--Accepts clusters if certain conditions are fulfilled
	Bool_t Accepted=0; //By default rejceted

	//double check these cuts carefully with the experts
	if(cluster->E()>0.3 && cluster->GetNCells()>1)
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
void AliAnalysisTaskGammaHadron::Fill_GH_Hisograms(Int_t identifier,TLorentzVector ClusterVec,AliVParticle* TrackVec, Double_t ClusterEcut, Double_t TrackPcut, Double_t Anglecut, Double_t Weight)
{
	//cout<<"Inside of: AliAnalysisTaskGammaHadron::Fill_GH_Hisograms()"<<endl;
	//--This function fills several histograms under different cut conditions.
	//--it is run within a cluster{ track{}} loop to get all combinations.

	//--A word to the weight - for mixed events it devides by the number of events in the current pool 1/nEvents
	//--                     - for same events you devide by the number of triggers
	//--                     - for both you have to take into account the efficiency of your correlated pair

	Double_t deltaEta = ClusterVec.Eta()-TrackVec->Eta();
	Double_t deltaPhi = DeltaPhi(ClusterVec,TrackVec);

	if(ClusterVec.Pt()>=ClusterEcut && TrackVec->Pt()>=TrackPcut && deltaPhi>=Anglecut)
	{
		fHistNoClus_ptH[identifier]->Fill(ClusterVec.Pt(),Weight);    //the .pt only works for gammas (E=M) for other particle this is wrong
		fHistpt_assHadron[identifier][0]->Fill(TrackVec->Pt(),Weight);
		fHist_DP_gh[identifier][0]      ->Fill(deltaPhi,Weight);
		fHist_dEta_dPhi[identifier]     ->Fill(deltaPhi,deltaEta,Weight);

		if(ClusterVec.Pt()<5)                        fHist_dEta_dPhi_low[identifier] ->Fill(deltaPhi,deltaEta,Weight);
		if(ClusterVec.Pt()>=5  && ClusterVec.Pt()<10)fHist_dEta_dPhi_med[identifier] ->Fill(deltaPhi,deltaEta,Weight);
		if(ClusterVec.Pt()>=10 && ClusterVec.Pt()<30)fHist_dEta_dPhi_high[identifier]->Fill(deltaPhi,deltaEta,Weight);
		//--fill histograms for different ranges of p_t^{g}
		for(Int_t i=1;i<fHistNoClus_ptH[identifier]->GetNbinsX()+1;i++)  // be careful hard coded from max[0] value -- need better variable // fHistNoClus_ptH->GetLast bin etc??
		{
			//--look in each p_T bin of the gamma (fHistNoClus_ptH histogram) how the pt
			//--distribution of the associated hadron looks like.
			Double_t BinWidth    = fHistNoClus_ptH[identifier]->GetBinWidth(i);
			Double_t BinValStart = fHistNoClus_ptH[identifier]->GetXaxis()->GetBinCenter(i)-BinWidth/2.0;

			//--right now the ranges are defined via different histogram bins of the gmma p_t histogram
			//--that might be changed later.
			if(ClusterVec.Pt()>=BinValStart && ClusterVec.Pt()<BinValStart+BinWidth)
			{
				fHistpt_assHadron[identifier][i]->Fill(TrackVec->Pt(),Weight);
				fHist_DP_gh[identifier][i]      ->Fill(deltaPhi,Weight);
			}
		}
	}
}
//the end
