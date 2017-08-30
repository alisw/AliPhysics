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
#include "TCustomBinning.h"

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

fGammaOrPi0(0),fSEvMEv(0),fDebug(0),fSavePool(0),
fEventCuts(0),fFiducialCuts(0x0),fFiducialCellCut(0x0), fHistEffGamma(0x0),fHistEffHadron(0x0),
fRtoD(0),
fClShapeMin(0),fClShapeMax(10),fMaxNLM(10),fRmvMTrack(0),fTrackMatchEta(0),fTrackMatchPhi(0),
fMixBCent(0),fMixBZvtx(),fPoolMgr(0x0),fTrackDepth(0),fPoolSize(0),fEventPoolOutputList(0),
fTriggerType(AliVEvent::kINT7), fMixingEventType(AliVEvent::kINT7),fCurrentEventTrigger(0),
fParticleLevel(kFALSE),fIsMC(kFALSE),
fEventCutList(0),

fHistPi0(0),fHistEvsPt(0),fHistClusPairInvarMasspT(0),fMAngle(0),fPtAngle(0),fHistBinCheckPt(0),fHistBinCheckZt(0),fHistBinCheckXi(0),

fHistDEtaDPhiGammaQA(0),fHistDEtaDPhiTrackQA(0), fHistClusterTime(0),

//fAODfilterBits(0),fHistptAssHadronG(0),fHistptAssHadronZt(0),fHistptAssHadronXi(0),fHistDEtaDPhiG(0),fHistDEtaDPhiZT(0),fHistDEtaDPhiXI(0)
//fHistptTriggG(),fHistptTriggZt(),fHistptTriggXi(),
fCorrVsManyThings(0), fClusterProp(0x0),
fHPoolReady(0x0),fDoRotBkg(0),fNRotBkgSamples(1),fPi0Cands(0),fClusEnergy(0),fRand(0)
{
	//..Initialize by defult for
	//..AliAnalysisTaskGammaHadron(0,0);
	InitArrays();

}
//________________________________________________________________________
AliAnalysisTaskGammaHadron::AliAnalysisTaskGammaHadron(Bool_t InputGammaOrPi0,Bool_t InputSeMe):
AliAnalysisTaskEmcal("AliAnalysisTaskGammaHadron", kTRUE),

fGammaOrPi0(0),fSEvMEv(0),fDebug(0),fSavePool(0),
fEventCuts(0),fFiducialCuts(0x0),fFiducialCellCut(0x0),fHistEffGamma(0x0),fHistEffHadron(0x0),
fRtoD(0),
fClShapeMin(0),fClShapeMax(10),fMaxNLM(10),fRmvMTrack(0),fTrackMatchEta(0),fTrackMatchPhi(0),
fMixBCent(0),fMixBZvtx(),fPoolMgr(0x0),fTrackDepth(0),fPoolSize(0),fEventPoolOutputList(0),
fTriggerType(AliVEvent::kINT7), fMixingEventType(AliVEvent::kINT7),fCurrentEventTrigger(0),
fParticleLevel(kFALSE),fIsMC(kFALSE),
fEventCutList(0),

fHistPi0(0),fHistEvsPt(0),fHistClusPairInvarMasspT(0),fMAngle(0),fPtAngle(0),fHistBinCheckPt(0),fHistBinCheckZt(0),fHistBinCheckXi(0),

fHistDEtaDPhiGammaQA(0),fHistDEtaDPhiTrackQA(0), fHistClusterTime(0),

//fAODfilterBits(0),fHistptAssHadronG(0),fHistptAssHadronZt(0),fHistptAssHadronXi(0),fHistDEtaDPhiG(0),fHistDEtaDPhiZT(0),fHistDEtaDPhiXI(0)
//fHistptTriggG(),fHistptTriggZt(),fHistptTriggXi(),
fCorrVsManyThings(0), fClusterProp(0x0),
fHPoolReady(0x0),fDoRotBkg(0),fNRotBkgSamples(1),fPi0Cands(0),fClusEnergy(0),fRand(0)
{
	InitArrays();
	//..set input variables
	fGammaOrPi0        =InputGammaOrPi0;
	fSEvMEv            =InputSeMe;
}
//________________________________________________________________________
void AliAnalysisTaskGammaHadron::InitArrays()
{
	//..Initialize by defult for
	//..AliAnalysisTaskGammaHadron(0,1);

	//..set input variables
	fGammaOrPi0        =0; //= 0 ( Gamma analysis ), 1 (pi0 analysis)
	fSEvMEv            =0; //= 0 (do only same event analyis with correct triggers), =1 (do event mixing)

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
		//delete fOutputList1;
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
	if(fSEvMEv==1 || fPoolMgr) //do this for either a mixed event analysis or when an external pool is given
	{
		InitEventMixer();
	}
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   Define sublists/folders for a better organisation of the figures
	//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	fOutputListQA   = new TList();
	fOutputListQA   ->SetOwner();
	fOutputListQA   ->SetName("QA_histograms");

  Double_t pi = TMath::Pi();

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
	//   THn Sparse for the 2D histograms
	//   Dimensions are eta,phi
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	Int_t dimThn = 0;
    TString titleThn[20];
    Int_t nbinsThn[20] = {0};
    Double_t minThn[20] = {0.};
    Double_t maxThn[20] = {0.};
    Double_t *binEdgesThn[20] = {0};

    titleThn[dimThn] = "#Delta #varphi";
    nbinsThn[dimThn] = 45;
    Double_t deltaPhiArray[45+1];
    binEdgesThn[dimThn] = deltaPhiArray;
    GenerateFixedBinArray(45,-90.,270.,deltaPhiArray);
    minThn[dimThn] = -90.;
    maxThn[dimThn] = 270.;
    dimThn++;

    titleThn[dimThn] = "#Delta #eta";
    nbinsThn[dimThn] = 80;
    //binEdgesThn[dim] = 80.;
    Double_t deltaEtaArray[80+1];
    binEdgesThn[dimThn] = deltaEtaArray;
    GenerateFixedBinArray(80,-2.,2.,deltaEtaArray);
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
/*
    titleThn[dimThn] = "ME(0) or SE(1)";
    static const Int_t nSeMeBins=2;
 	Double_t SeMeBinArray[nSeMeBins+1]  = {0,1,2};
    nbinsThn[dimThn] = nSeMeBins;
    binEdgesThn[dimThn] = SeMeBinArray;
    minThn[dimThn] = SeMeBinArray[0];
    maxThn[dimThn] = SeMeBinArray[nSeMeBins];
    dimThn++;
*/
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
    if(fPlotQA!=1)
    {
     	fCorrVsManyThings   = new THnSparseF("CorrVsManyThings", "CorrVsManyThings", dimThn, nbinsThn, minThn, maxThn);
     	for(Int_t i=0;i<dimThn;i++)
     	{
     		fCorrVsManyThings->GetAxis(i)->SetTitle(titleThn[i]);
     		fCorrVsManyThings->SetBinEdges(i, binEdgesThn[i]);
     	}
     	//fCorrVsManyThings->Sumw2();
     	fOutput->Add(fCorrVsManyThings);
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   THn Sparse for the Pi0 Candidates
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Int_t dimThnPi0 = 0;
    TString titleThnPi0[10];
    Int_t nBinsThnPi0[10] = {0};
    Double_t minThnPi0[10] = {0.};
    Double_t maxThnPi0[10] = {0.};
    Double_t *binEdgesThnPi0[10] = {0};

    titleThnPi0[dimThnPi0] = "p_{T}^{#gamma#gamma}";
    nBinsThnPi0[dimThnPi0] = 100;
    Double_t pTArray[100+1];
    binEdgesThnPi0[dimThnPi0] = pTArray;
    GenerateFixedBinArray(100,0,50,pTArray);
    minThnPi0[dimThnPi0] = 0;
    maxThnPi0[dimThnPi0] = 50;
    dimThnPi0++;

    titleThnPi0[dimThnPi0] = "M_{#gamma#gamma}";
    nBinsThnPi0[dimThnPi0] = 750;
    Double_t mGGArray[750+1];
    binEdgesThnPi0[dimThnPi0] = mGGArray;
    GenerateFixedBinArray(750,0,0.75,mGGArray);
    minThnPi0[dimThnPi0] = 0;
    maxThnPi0[dimThnPi0] = 0.75;
    dimThnPi0++;

    titleThnPi0[dimThnPi0] = "Opening Angle [rad]";
    nBinsThnPi0[dimThnPi0] = 23;
    Double_t openAngleArray[23+1] = {0,0.009,0.011,0.013,0.015,0.017,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.125,0.15,pi/16.,pi/8.,pi/4.,3.*pi/8.,pi/2.,3.*pi/4.,pi};
    binEdgesThnPi0[dimThnPi0] = openAngleArray;
  //  GenerateFixedBinArray(20,0,TMath::Pi(),openAngleArray);
    minThnPi0[dimThnPi0] = 0;
    maxThnPi0[dimThnPi0] = pi;
    dimThnPi0++;

    titleThnPi0[dimThnPi0] = "Max Lambda_{0}^{2}";
    nBinsThnPi0[dimThnPi0] = 6;
    Double_t MaxM02Array[6+1] = {0.1,0.3,0.35,0.4,0.45,0.5,10}; 
    binEdgesThnPi0[dimThnPi0] = MaxM02Array;
   // GenerateFixedBinArray(5,0,);
    minThnPi0[dimThnPi0] = 0.1;
    maxThnPi0[dimThnPi0] = 10;
    dimThnPi0++;

    titleThnPi0[dimThnPi0] = "Min Cluster Energy";
    nBinsThnPi0[dimThnPi0] = 5;
    Double_t MinClusEnergyArray[5+1] = {0.3,0.5,1,1.5,2,100};
    binEdgesThnPi0[dimThnPi0] = MinClusEnergyArray;
   // GenerateFixedBinArray(5,0,);
    minThnPi0[dimThnPi0] = 0.3;
    maxThnPi0[dimThnPi0] = 100;
    dimThnPi0++;

    titleThnPi0[dimThnPi0] = "Asymmetry";
    nBinsThnPi0[dimThnPi0] = 5;
    Double_t asymmetryArray[5+1];
    binEdgesThnPi0[dimThnPi0] = asymmetryArray;
    GenerateFixedBinArray(5,0,1,asymmetryArray);
    minThnPi0[dimThnPi0] = 0;
    maxThnPi0[dimThnPi0] = 1;
    dimThnPi0++;

    //..ID array -  0 - real, 1 - rotated
    titleThnPi0[dimThnPi0] = "Rotation Status";
    nBinsThnPi0[dimThnPi0] = 2;
    Double_t mRotArray[2+1];
    binEdgesThnPi0[dimThnPi0] = mRotArray;
    GenerateFixedBinArray(2,0,2,mRotArray);
    minThnPi0[dimThnPi0] = 0;
    maxThnPi0[dimThnPi0] = 2;
    dimThnPi0++;

    if ( fGammaOrPi0  && fPlotQA==1) {
      fPi0Cands= new THnSparseF("Pi0Cands", "Pi0Cands", dimThnPi0, nBinsThnPi0, minThnPi0, maxThnPi0);
      for(Int_t i=0;i<dimThnPi0;i++)
      {
        fPi0Cands->GetAxis(i)->SetTitle(titleThnPi0[i]);
        fPi0Cands->SetBinEdges(i, binEdgesThnPi0[i]);
      }
      fOutput->Add(fPi0Cands);
    }


	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //   THn Sparse for the Cluster properties
	//   Dimensions are Cluster Energy, ...
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	Int_t dimThnQA = 0;
    TString titleThnQA[11];
    Int_t nbinsThnQA[11] = {0};
    Double_t minThnQA[11] = {0.};
    Double_t maxThnQA[11] = {0.};
    Double_t *binEdgesThnQA[11] = {0};

    titleThnQA[dimThnQA] = "E_{#gamma}";
    nbinsThnQA[dimThnQA] = 150;
    Double_t EgArray[150+1];
    binEdgesThnQA[dimThnQA] = EgArray;
    GenerateFixedBinArray(150,0,30,EgArray);
    minThnQA[dimThnQA] = 0;
    maxThnQA[dimThnQA] = 30;
    dimThnQA++;


    //..Create the fhAmpId TH2D with increasing binwidth
    //..0-10 GeV (0.05), 10-20 GeV (0.2), 20-30 GeV (0.5)
    //Double_t binWidth=(ptfinemax-ptfinemin)/nfineptbins;
    TCustomBinning xBinning;
    xBinning.SetMinimum(0);
    xBinning.AddStep(0.5,0.005);   //..first entries of the array are the set ranges and bins
    xBinning.AddStep(1,0.02);      //..expand the previously defined range by 2 but increase the bin width
    xBinning.AddStep(4,0.04);      //..expand the previously defined range by 4 but increase the bin width

    TArrayD xbinsArray;
    xBinning.CreateBinEdges(xbinsArray);

    titleThnQA[dimThnQA] = "#lambda_{0}";
    nbinsThnQA[dimThnQA] = 200;
    Double_t ShapeArray[200+1];
	for(Int_t i=0;i<201;i++)
    {
    		ShapeArray[i]=xbinsArray.At(i);
    }
    binEdgesThnQA[dimThnQA] = ShapeArray;
    minThnQA[dimThnQA] = 0;
    maxThnQA[dimThnQA] = 4;
    dimThnQA++;
/*
    titleThnQA[dimThnQA] = "NLM";
    nbinsThnQA[dimThnQA] = 4;
    Double_t NLMArray[4+1];
    binEdgesThnQA[dimThnQA] = NLMArray;
    GenerateFixedBinArray(4,0,4,NLMArray);
    minThnQA[dimThnQA] = 0;
    maxThnQA[dimThnQA] = 4;
    dimThnQA++;
*/
/*
    titleThnQA[dimThnQA] = "#Cells";
    nbinsThnQA[dimThnQA] = 10;
    Double_t nCellsArray[10+1];
    binEdgesThnQA[dimThnQA] = nCellsArray;
    GenerateFixedBinArray(10,0,10,nCellsArray);
    minThnQA[dimThnQA] = 0;
    maxThnQA[dimThnQA] = 10;
    dimThnQA++;
*/
    titleThnQA[dimThnQA] = "cell distance to bad channel";
    nbinsThnQA[dimThnQA] = 5;
    Double_t distanceArray[5+1];
    binEdgesThnQA[dimThnQA] = distanceArray;
    GenerateFixedBinArray(5,0,5,distanceArray);
    minThnQA[dimThnQA] = 0;
    maxThnQA[dimThnQA] = 5;
    dimThnQA++;

    titleThnQA[dimThnQA] = "cell distance to SM edge";
    nbinsThnQA[dimThnQA] = 5;
    Double_t distanceArray2[5+1];
    binEdgesThnQA[dimThnQA] = distanceArray2;
    GenerateFixedBinArray(5,0,5,distanceArray2);
    minThnQA[dimThnQA] = 0;
    maxThnQA[dimThnQA] = 5;
    dimThnQA++;
/*
    titleThnQA[dimThnQA] = "#Delta #eta^(match. track-cluster)";
    nbinsThnQA[dimThnQA] = 50;
    Double_t etaArrayDistMatched[50+1];
    binEdgesThnQA[dimThnQA] = etaArrayDistMatched;
    GenerateFixedBinArray(50,-0.05, 0.05,etaArrayDistMatched);
    minThnQA[dimThnQA] = 0;
    maxThnQA[dimThnQA] = 50;
    dimThnQA++;

    titleThnQA[dimThnQA] = "#Delta #varphi^(match. track-cluster)";
    nbinsThnQA[dimThnQA] = 50;
    Double_t phiArrayDistMatched[50+1];
    binEdgesThnQA[dimThnQA] = phiArrayDistMatched;
    GenerateFixedBinArray(50,-0.05, 0.05,phiArrayDistMatched);
    minThnQA[dimThnQA] = 0;
    maxThnQA[dimThnQA] = 50;
    dimThnQA++;
*/
    titleThnQA[dimThnQA] = "#eta^{cluster}";
    nbinsThnQA[dimThnQA] = 142;
    Double_t etaArray[142+1];
    binEdgesThnQA[dimThnQA] = etaArray;
    GenerateFixedBinArray(142,-0.71,0.71,etaArray);
    minThnQA[dimThnQA] = 0;
    maxThnQA[dimThnQA] = 142;
    dimThnQA++;

    titleThnQA[dimThnQA] = "#varphi^{cluster}";
    nbinsThnQA[dimThnQA] = 311;
    Double_t phiArray[311+1];
    binEdgesThnQA[dimThnQA] = phiArray;
    GenerateFixedBinArray(311,75.9,331.1,phiArray);
    minThnQA[dimThnQA] = 0;
    maxThnQA[dimThnQA] = 311;
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
    if(fPlotQA==1 &&  !fGammaOrPi0)
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
	//    Histograms for common use
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

	//..Initialize
	fHistBinCheckPt       = new TH1*[kNIdentifier];
	fHistBinCheckZt       = new TH1*[kNIdentifier];
	fHistBinCheckXi       = new TH1*[kNIdentifier];
	fHistDEtaDPhiGammaQA  = new TH2*[kNIdentifier+3]; //Why +2?? -> because I want more than 3 QA versions
	fHistDEtaDPhiTrackQA  = new TH2*[kNIdentifier+3];
	fHistClusterTime      = new TH2*[kNIdentifier+3];

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
	}

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//   Special QA histograms (also to get more info what is going on in mixed event for trigger data)
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

		//..Time information
		fHistClusterTime[identifier] = new TH2F(Form("fHistClusterTime%d_Id%d",0,identifier),Form("fHistClusterTime%d_Id%d",0,identifier),2000,-100,100,200,0,40);
		fHistClusterTime[identifier]->GetXaxis()->SetTitle("time [ns]");
		fHistClusterTime[identifier]->GetYaxis()->SetTitle("pT");
		fOutputListQA->Add(fHistClusterTime[identifier]);
	}

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //   Michael's Special Histograms
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if ( fGammaOrPi0 ) {  // Don't necessarily need this for Gamma analysis
    	fClusEnergy = new TH1F("ClusEnergy","Cluster Energy",1000,0,50);
    	fClusEnergy->GetXaxis()->SetTitle("E (GeV)");
    	fOutput->Add(fClusEnergy);
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

	//..The END
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
	if (fSEvMEv==0 && !fCaloClusters)                         return kFALSE;
	if (fSEvMEv==0 && !(fCurrentEventTrigger & fTriggerType)) return kFALSE;

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
	if(fSEvMEv==1)
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
	if(fSEvMEv==0 && fCurrentEventTrigger & fTriggerType)
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
		if(!cluster)continue; //check if the cluster is a good cluster

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
		//------------------------------------------------

		//...........................................
		//..combine gammas with same event tracks
		GammaCounter++;
		if(SameMix==1)
		{
			fHistEvsPt->Fill(CaloClusterVec.Pt(),CaloClusterVec.E()); //the .pt only works for gammas (E=M) for other particle this is wrong

			if(!tracks)  return 0;
			Int_t NoOfTracksInEvent =tracks->GetNParticles();
			AliVParticle* track=0;

			Int_t trackCounter=0;
			for(Int_t NoTrack = 0; NoTrack < NoOfTracksInEvent; NoTrack++)
			{
				track = (AliVParticle*)tracks->GetAcceptParticle(NoTrack);
				if(!track)continue; //check if the track is a good track
				trackCounter++;

				FillGhHisograms(0,aliCaloClusterVec,track,5,Weight);
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
		//...........................................
		//..Additional histograms
		fFiducialCellCut->SetNumberOfCellsFromEMCALBorder(1);
		if(!AccClusterForAna(clusters,cluster))continue; //check if the cluster is a good cluster
		if(SameMix==1)FillQAHisograms(2,clusters,cluster,trackNULL);

		fFiducialCellCut->SetNumberOfCellsFromEMCALBorder(2);
		if(!AccClusterForAna(clusters,cluster))continue; //check if the cluster is a good cluster
		if(SameMix==1)FillQAHisograms(3,clusters,cluster,trackNULL);

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

 // Double_t ClusterEnergyCut = 1; // On top of cuts in GetAcceptCluster,AccClusterForAna

	AliVCluster* cluster = 0;
	AliVCluster* cluster2= 0;
	AliVParticle* trackNULL=0;

  Weight = InputWeight; // Good enough for now.

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

      fClusEnergy->Fill(cluster->GetNonLinCorrEnergy(),Weight);

			TLorentzVector CaloClusterVec;
			clusters->GetMomentum(CaloClusterVec, cluster);
			//acc if pi0 candidate
			nAccClusters++;

			//for(Int_t NoCluster2 = 0; NoCluster2 < NoOfClustersInEvent; NoCluster2++ )
			for(Int_t NoCluster2 = NoCluster1 + 1; NoCluster2 < NoOfClustersInEvent; NoCluster2++ )
			{
//				if(NoCluster1!=NoCluster2)
//				{
					cluster2=(AliVCluster*) clusters->GetAcceptCluster(NoCluster2);
					if(!cluster2 || !AccClusterForAna(clusters,cluster2))continue; //check if the cluster is a good cluster

					TLorentzVector CaloClusterVec2;
					TLorentzVector CaloClusterVecpi0;
          Double_t fMaxClusM02 = TMath::Max(cluster->GetM02(),cluster2->GetM02());

					//old framework				cluster2->GetMomentum(CaloClusterVec2, fVertex);
					clusters->GetMomentum(CaloClusterVec2, cluster2);
				//	if(cluster2->E()>2 && cluster->E()>2)
         // if(cluster2->GetNonLinCorrEnergy()>ClusterEnergyCut && cluster->GetNonLinCorrEnergy()>ClusterEnergyCut)
				//	{
						CaloClusterVecpi0=CaloClusterVec+CaloClusterVec2;
						fHistPi0->Fill(CaloClusterVecpi0.M());
            FillPi0CandsHist(CaloClusterVec,CaloClusterVec2,CaloClusterVecpi0,fMaxClusM02,Weight);
						fHistClusPairInvarMasspT->Fill(CaloClusterVecpi0.M(),CaloClusterVecpi0.Pt());
						fMAngle->Fill(CaloClusterVecpi0.M(), CaloClusterVec.Angle(CaloClusterVec2.Vect()),0.5);
						fPtAngle->Fill(CaloClusterVecpi0.Pt(), CaloClusterVec.Angle(CaloClusterVec2.Vect()),0.5);
						if((CaloClusterVecpi0.M()>=Pi0Mass-Pi0Window) && (CaloClusterVecpi0.M()<=Pi0Mass+Pi0Window)){
							nAccPi0Clusters++; //need eventually to divide by 2, otherwise double counting the pi0's
						}
				//	}
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
/// Fill Pi0 Cand THnSparse with relevant info
/// To Do: add in rotation method
///
//________________________________________________________________________
void AliAnalysisTaskGammaHadron::FillPi0CandsHist(AliTLorentzVector CaloClusterVec, AliTLorentzVector CaloClusterVec2, AliTLorentzVector CaloClusterVecPi0, Double_t fMaxClusM02, Double_t Weight) {

  Double_t valueArray[7];
  valueArray[0]=CaloClusterVecPi0.Pt();
  valueArray[1]=CaloClusterVecPi0.M();
  valueArray[2]=CaloClusterVec.Angle(CaloClusterVec2.Vect());

  Double_t fE1 = CaloClusterVec.E();
  Double_t fE2 = CaloClusterVec2.E();
  Double_t fAssym = (fE1+fE2 > 0.000001) ? TMath::Abs(fE2-fE1)/(fE1+fE2) : 0; //Don't divide by zero

  valueArray[3]=fMaxClusM02;
  valueArray[4]=TMath::Min(fE1,fE2);
  valueArray[5]=fAssym;
  valueArray[6]=0;
  
  fPi0Cands->Fill(valueArray,Weight);

  if (!fDoRotBkg) return;
  // Rotational Background
//  const Double_t fOpeningAngleCut = 0.017;

  if (!fRand) fRand = new TRandom3(0);

  for (int i = 0; i < fNRotBkgSamples; i++) {
    Double_t fEta,fPhi;
    Double_t fOpeningAngle;
    while (true) {
      fEta = fRand->Uniform(-0.7,0.7);  // change to eta cut maybe
      // GetClusterContainer("caloClusters")->SetMinEta() 
      // GetClusterContainer("caloClusters")->SetMaxEta() 
      fPhi = fRand->Uniform(80,254); // pretend DCAL next to EMCAL
      if(fPhi > 187) {
        // Check PHOS hole
        if (TMath::Abs(fEta) < .22) continue;
        fPhi+= 73;  // shift DCAL points
      } 
      fPhi = fPhi * 3.141592653589793 / 180.;
      // Opening Angle Cut
      CaloClusterVec2.SetPhi(fPhi);
      CaloClusterVec2.SetTheta(2.*TMath::ATan(TMath::Exp(-fEta)));
      fOpeningAngle = CaloClusterVec.Angle(CaloClusterVec2.Vect());
  //    if (fOpeningAngle < fOpeningAngleCut) continue;
      break;
    }

    CaloClusterVec2.SetPhi(fPhi);
    CaloClusterVec2.SetTheta(2.*TMath::ATan(TMath::Exp(-fEta)));

    CaloClusterVecPi0=CaloClusterVec+CaloClusterVec2;

    valueArray[0]=CaloClusterVecPi0.Pt();
    valueArray[1]=CaloClusterVecPi0.M();
    valueArray[2]=fOpeningAngle;
    valueArray[3]=fMaxClusM02;
    valueArray[4]=TMath::Min(fE1,fE2);
    valueArray[5]=fAssym;
    valueArray[6]=1;

    fPi0Cands->Fill(valueArray,Weight);
  } 
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
	//..Careful here: usually this is done for an opening angle (hadron-jet axis) of less than 90¡. Due to
	//..resolution momentum smearing (our guess - check that!) there are particles appearing at angles greater than 90¡
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

	Double_t valueArray[8];
	valueArray[0]=deltaPhi;
	valueArray[1]=deltaEta;
	valueArray[2]=G_PT_Value;
	valueArray[3]=ZT_Value;
	valueArray[4]=XI_Value;
	valueArray[5]=zVertex;
	valueArray[6]=evtPlaneCategory;
	valueArray[7]=fCent;

	if(G_PT_Value>=ClusterEcut)
	{
		if(identifier==0 && fPlotQA!=1)fCorrVsManyThings  ->Fill(valueArray,Weight);

		//..Histograms to test the binning
		fHistBinCheckPt[identifier] ->Fill(G_PT_Value,Weight);
		fHistBinCheckZt[identifier] ->Fill(ZT_Value,Weight);
		fHistBinCheckXi[identifier] ->Fill(XI_Value,Weight);
	}
}
//________________________________________________________________________
void AliAnalysisTaskGammaHadron::FillQAHisograms(Int_t identifier,AliClusterContainer* clusters,AliVCluster* caloCluster,AliVParticle* TrackVec)
{
	TLorentzVector caloClusterVec;
	clusters->GetMomentum(caloClusterVec,caloCluster);
	AliTLorentzVector aliCaloClusterVec = AliTLorentzVector(caloClusterVec); //..can acess phi from
	Double_t energy= caloCluster->GetNonLinCorrEnergy();

	if(identifier==0 && fPlotQA==1 && fGammaOrPi0==0 && energy >=1)
	{
		//..Get leading cluster
		AliVCluster* leadingClus = GetLeadingCluster("Pt",clusters);

		Double_t etaDistMatched=-1;
		Double_t phiDistMatched=-1;

		//..Distance to SM border
		Int_t etaCellDist;
		Int_t phiCellDist;
		GetDistanceToSMBorder(caloCluster,etaCellDist,phiCellDist);
		if(etaCellDist<0)cout<<"etaCellDist: "<<etaCellDist<<endl;
		if(phiCellDist<0)cout<<"phiCellDist: "<<phiCellDist<<endl;
		Int_t minCellDistance=etaCellDist;
		if(etaCellDist>phiCellDist)minCellDistance=phiCellDist;

		//..ID array -  1 - leading, 2- track matched, 3 - leading & track matched
		Int_t gammaInfo=0;
		if(caloCluster==leadingClus)gammaInfo=1;
		//cout<<"cluster ID: "<<caloCluster->GetID()<<", lead cluster ID: "<<leadingClus->GetID()<<endl;
		if(DetermineMatchedTrack(caloCluster,etaDistMatched,phiDistMatched))gammaInfo=2;
		if(gammaInfo==2 && caloCluster==leadingClus)gammaInfo=3;
		//cout<<"eta distance matched: "<<etaDistMatched<<", phi dist matched: "<<phiDistMatched<<endl;

		//Eg, lambda0,NLM, ncells, distance to bad ,e/p, Mgg
		Double_t valueArray[7];
		valueArray[0] = energy;
		valueArray[1] = caloCluster->GetM02();
		//valueArray[2] = caloCluster->GetNExMax();
		//valueArray[3] = caloCluster->GetNCells();
		valueArray[2] = caloCluster->GetDistanceToBadChannel()-1; //..shift to -1 since it starts at 1 and not at 0
		valueArray[3] = minCellDistance;
		//valueArray[6] = phiDistMatched;
		//valueArray[7] = etaDistMatched;
		valueArray[4] = caloClusterVec.Eta();
		valueArray[5] = aliCaloClusterVec.Phi_0_2pi()*fRtoD;
		valueArray[6] = gammaInfo;

		//valueArray[5]=0;//E/p
		//valueArray[6]=130;//m_gg
		fClusterProp->Fill(valueArray); //..all clusters - no cuts
	}
	/*do similar test here?*/fHistDEtaDPhiGammaQA[identifier] ->Fill(caloClusterVec.Eta(),aliCaloClusterVec.Phi_0_2pi()*fRtoD);
	if(TrackVec)             fHistDEtaDPhiTrackQA[identifier] ->Fill(TrackVec->Eta(),TrackVec->Phi()*fRtoD);
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
	Double_t trash;
	if(fRmvMTrack==1 && DetermineMatchedTrack(caloCluster,trash,trash))
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
	//--shift the second peak over the fist peak: \--Æ--/   --> -Æ--
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
	//--shift the second peak over the fist peak: \--Æ--/   --> -Æ--
	//--to create a histogram that starts at -pi/2 and ends at 3/2pi
	if (dPhi <= -TMath::Pi()/2)    dPhi += 2*pi;
	if (dPhi > 3.0*TMath::Pi()/2.0)dPhi -= 2*pi;

	//--change from rad to degree:
	dPhi*= fRtoD;

	return dPhi;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskGammaHadron::DetermineMatchedTrack(AliVCluster* caloCluster,Double_t &etadiff,Double_t & phidiff)
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
		etadiff=veta-ceta;
		phidiff=TVector2::Phi_mpi_pi(vphi-cphi);

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
//
// Calculate the distance to a supermodule border
// calculated is the smaller distance in either direction of the max. energy cell of the cluster
//
//________________________________________________________________________
void AliAnalysisTaskGammaHadron::GetDistanceToSMBorder(AliVCluster* caloCluster,Int_t &etaCellDist,Int_t &phiCellDist)
{
	//..celldist 0 means it is at the border of the EMCal
	//..celldist 1 means there is one row/collum between the cell and the border
	etaCellDist=-1;
	phiCellDist=-1;
	//..partially copied from RecoUtils
	Int_t absIdMax  = -1, iSM =-1, ieta = -1, iphi = -1;
	Bool_t shared = kFALSE;
	fFiducialCellCut->GetMaxEnergyCell(fGeom, fCaloCells, caloCluster, absIdMax,  iSM, ieta, iphi, shared);
	if (absIdMax==-1) cout<<"In GetDistanceToSMBorder  why??"<<endl;

	//cout<<"max cell eta collumn: "<<ieta<<", max cell phi row: "<<iphi<<endl;
	//..Define Borders
	Int_t lastRow,firstRow;
	Int_t lastCollumn,firstCollumn;
	//..phi
	firstRow=0;
	lastRow=24;
	if      (fGeom->GetSMType(iSM) == AliEMCALGeometry::kEMCAL_Half) lastRow /= 2;
	else if (fGeom->GetSMType(iSM) == AliEMCALGeometry::kEMCAL_3rd ) lastRow /= 3;// 1/3 sm case
	else if (fGeom->GetSMType(iSM) == AliEMCALGeometry::kDCAL_Ext  ) lastRow /= 3;// 1/3 sm case
	//..eta
	firstCollumn=0;
	lastCollumn=48;
	if(fGeom->GetSMType(iSM) == AliEMCALGeometry::kDCAL_Standard )  lastCollumn = lastCollumn*2/3;

	lastRow    =lastRow-1;      //..range starts from 0 and only goes up to 23
	lastCollumn=lastCollumn-1;  //..range starts from 0

	//..Calculate smallest distance
	//..phi
	phiCellDist=iphi; //eli in case it starts from 1
    if((iphi-firstRow)>=(lastRow-iphi))phiCellDist=lastRow-iphi;

    //..eta
    if(iSM%2==0)
    {
    	  //..only outer border (not at eta=0) except DCal modules
      etaCellDist=ieta;
      if(lastCollumn<47 && (ieta-firstCollumn)>=(lastCollumn-ieta))etaCellDist=lastCollumn-ieta;
    }
    else
    {
    	  //..only outer border (not at eta=0) except DCal modules
    	  etaCellDist=lastCollumn-ieta;
      if(lastCollumn<47 && (ieta-firstCollumn)<=(lastCollumn-ieta))etaCellDist=ieta-firstCollumn;
    }
    //cout<<"supermodule: "<<iSM<<", last collumn: "<<lastCollumn<<endl;
    //cout<<"eta distance: "<<etaCellDist<<" phi distance: "<<phiCellDist<<endl;
}
//
// Get leading cluster of clusters that pass
// the standard cuts for the analysis
//________________________________________________________________________
AliVCluster* AliAnalysisTaskGammaHadron::GetLeadingCluster(const char* opt,AliClusterContainer* clusters)
{
	//..2 options are leading by pT and leading by E
	TString option(opt);

	AliVCluster *clusterMaxE;
	AliVCluster *clusterMaxPt;
	AliVCluster *cluster = 0;
	Double_t et=0;
	Double_t etmax=0;

	Bool_t first=1;
	for(Int_t NoCluster1 = 0; NoCluster1 < clusters->GetNClusters(); NoCluster1++ )
	{
		cluster=(AliVCluster*) clusters->GetAcceptCluster(NoCluster1); //->GetCluster(NoCluster1);
		if(!cluster)continue;                            //check if the cluster is a good cluster
		if(!AccClusterForAna(clusters,cluster))continue; //check if the cluster is a good cluster

		if(first==1)
		{
			clusterMaxE=cluster;
			clusterMaxPt=cluster;
			first=0;
		}

		if (cluster->E() > clusterMaxE->E()) clusterMaxE = cluster;

		TLorentzVector mom;
		clusters->GetMomentum(mom,cluster);
		et = mom.Et();
		if (et > etmax)
		{
			clusterMaxPt = cluster;
			etmax = et;
		}
	}
	//if (option.Contains("Pt"))    cout<<"selected pt option"<<endl;
	//else if (option.Contains("E"))cout<<"selected E option"<<endl;
	//else  cout<<"no option selected"<<endl;

	if (option.Contains("Pt"))    return clusterMaxPt;
	else if (option.Contains("E"))return clusterMaxE;
	else return 0;
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
