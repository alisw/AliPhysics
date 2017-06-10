//
#include <Riostream.h>
#include <TClonesArray.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "AliAnalysisTaskEMCALPi0GammaCorr.h"

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
#include "AliEMCALGeometry.h"
#include "AliEMCALGeoParams.h"

#include <memory>
using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskEMCALPi0GammaCorr)

////////////////////////////////////////////////////////////////////////////////////////
AliAnalysisTaskEMCALPi0GammaCorr::AliAnalysisTaskEMCALPi0GammaCorr():
AliAnalysisTaskEmcal("AliAnalysisTaskEMCALPi0GammaCorr", kTRUE),
fSavePool(0),
fEventCuts(0),
fHistEffGamma(0x0),
fHistEffHadron(0x0),
fMixBCent(0),
fMixBZvtx(),
fPoolMgr(0x0),
fTrackDepth(0),
fPoolSize(0),
fEventPoolOutputList(0),
fTriggerType(AliVEvent::kINT7), 
fMixingEventType(AliVEvent::kINT7),
fCurrentEventTrigger(0),
fEventCutList(0),
h_Track(0),
h_Cluster(0),
h_ClusterTrack(0),
h_ClusterTrack_Mixed(0),
h_Pi0(0),
h_Pi0Track(0),
h_Pi0Track_Mixed(0)
{
    InitArrays();
}

// -------------------------------------------------------------------------------------
// Constructor with inputs
AliAnalysisTaskEMCALPi0GammaCorr::AliAnalysisTaskEMCALPi0GammaCorr(Bool_t InputDoMixing):
AliAnalysisTaskEmcal("AliAnalysisTaskEMCALPi0GammaCorr", kTRUE),
fSavePool(0),
fEventCuts(0),
fHistEffGamma(0x0),
fHistEffHadron(0x0),
fMixBCent(0),
fMixBZvtx(),
fPoolMgr(0x0),
fTrackDepth(0),
fPoolSize(0),
fEventPoolOutputList(0),
fTriggerType(AliVEvent::kINT7), 
fMixingEventType(AliVEvent::kINT7),
fCurrentEventTrigger(0),
fEventCutList(0),
h_Track(0),
h_Cluster(0),
h_ClusterTrack(0),
h_ClusterTrack_Mixed(0),
h_Pi0(0),
h_Pi0Track(0),
h_Pi0Track_Mixed(0)
{
	InitArrays();
}//End constructor PiHadron that receives input

void AliAnalysisTaskEMCALPi0GammaCorr::InitArrays()
{
    AliWarning("InitArrays is being called");
	fSavePool          =0; //= 0 do not save the pool by default. Use the set function to do this.
	fUseManualEventCuts=1; //=0 use automatic setting from AliEventCuts. =1 load manual cuts
    //Setting bins for the mixing of events.
    double centmix[kNcentBins+1] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0};
	fMixBCent = new TAxis(kNcentBins,centmix);

    double zvtxmix[kNvertBins+1] = {-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10};
	memcpy (fArrayNVertBins, zvtxmix, sizeof (fArrayNVertBins));
	fMixBZvtx = new TAxis(kNvertBins,zvtxmix);
    fTrackDepth     = 50000;    //Raymonds/Megans value
    fPoolSize       = 1; 
    SetMakeGeneralHistograms(kTRUE); 
} //end of function init arrays

AliAnalysisTaskEMCALPi0GammaCorr::~AliAnalysisTaskEMCALPi0GammaCorr()
{

}

void AliAnalysisTaskEMCALPi0GammaCorr::UserCreateOutputObjects()
{
    AliWarning("Entering UserCreateOutPutObjects");   
	AliAnalysisTaskEmcal::UserCreateOutputObjects();
    
    fEventCutList = new TList();
	fEventCutList ->SetOwner();
	fEventCutList ->SetName("EventCutOutput");

	fEventCuts.OverrideAutomaticTriggerSelection(fOffTrigger);
    
    if(fUseManualEventCuts==1)
	{  
    AliWarning("Setting Manual Event Cuts"); 
    
    fEventCuts.SetManualMode();
    fEventCuts.fCentralityFramework=2; //..only for Run1!!
    fEventCuts.fTriggerMask = fOffTrigger;
    fEventCuts.fMinVtz = fMinVz;
    fEventCuts.fMaxVtz = fMaxVz;
    fEventCuts.fRequireTrackVertex = true;
    fEventCuts.fMaxDeltaSpdTrackAbsolute=fZvertexDiff;
    fEventCuts.fTrackletBGcut = fTklVsClusSPDCut; //(false by default for 15o)
    fEventCuts.fMinCentrality = fMinCent;
    fEventCuts.fMaxCentrality = fMaxCent;
    }
    fEventCuts.AddQAplotsToList(fEventCutList);
    fOutput->Add(fEventCutList);
    //OutputList->Add(fEventCutList);
    
    AliWarning("Initializing Event Mixer");
    InitEventMixer();
	
    //Initializing the histograms to be saved. For the moment, only pT of clusters and Mgammagamma.
    int nbins_Mass = 100;
    int nbins_Pt   = 100;
    int nbins_E    = 100;
    int nbins_dphi     = 18;
    int nbins_deta     = 30;
    int nbins_phi     = 36;
    int nbins_eta     = 40;
    int nbins_zt      = 50;
    int nbins_xi      = 50;
    int nbins_M02     = 50;
    int nbins_Ncells  = 30;
    int nbins_Centrality = 10;
    int nbins_zvertex = 20;
    int nbins_Asymmetry = 40;
    int nbins_nMatchedTracks = 5;
    
    double min_Mass = 0;
    double max_Mass = 1.0;
    double min_Pt = 0;
    double max_Pt = 50.0;
    double min_dphi = -0.5; // rads
    double max_dphi = 1.5; // rads 
    double min_phi = -1.0*TMath::Pi();
    double max_phi = TMath::Pi();

    double min_zvertex = -10;
    double max_zvertex = +10;
    double max_deta = 1.5;
    double max_eta =  1.0;
    double min_deta = -max_deta;
    double min_eta =  -max_eta;
    
    double min_zt = 0;
    double min_xi = 0;
    double max_zt = 2.0;
    double max_xi = 5.0;
    
    double min_E =0;
    double max_E =50;
    double min_M02 = 0.0;
    double max_M02 = 2.0;
    
    double min_Ncells = -0.5;
    double max_Ncells = 29.5;
     
    double min_Centrality = 0;
    double max_Centrality = 100;
    
    double min_Asymmetry =0;
    double max_Asymmetry = 1.0;
    
    double min_nMatchedTracks = -1.5;
    double max_nMatchedTracks  =  3.5;
    
    int    nbins_nMaxima = 5;
    double min_nMaxima = -0.5;
    double max_nMaxima  =  4.5;
    
    int nbins_alpha = 100;
    double min_alpha = 0.0;
    double max_alpha = 0.5;
    
    //Pion-hadron correlations
    const int nbins_PionCorr = 24;
    
    int bins[nbins_PionCorr]    = {nbins_Centrality, nbins_zvertex, nbins_Pt, nbins_E, nbins_eta, nbins_eta, nbins_phi, //trigger variables
                                   nbins_Pt,  nbins_eta, nbins_phi, nbins_dphi, nbins_deta, nbins_deta/2, nbins_zt, nbins_xi, //track variables
                                   nbins_Mass,  nbins_Pt, nbins_Pt, nbins_eta, nbins_eta, nbins_phi, nbins_phi, nbins_M02, nbins_M02}; //pion-only and pion-decay variables
                                   
    double xmin[nbins_PionCorr] = {min_Centrality, min_zvertex, min_Pt  , min_E,  min_eta, min_eta, min_phi, 
                                   min_Pt, min_eta, min_phi    , min_dphi   , min_deta   , 0.0, min_zt, min_xi,
                                   min_Mass, min_Pt, min_Pt, min_eta, min_eta, min_phi, min_phi, min_M02, min_M02};
                                   
    double xmax[nbins_PionCorr] = {max_Centrality, max_zvertex,  max_Pt,  max_E, max_eta    , max_eta, max_phi, 
                                   max_Pt, max_eta,  max_phi   , max_dphi,   max_deta, max_deta, max_zt, max_xi,
                                   max_Mass, max_Pt, max_Pt, max_eta, max_eta, max_phi, max_phi, max_M02, max_M02};

    TString axisNames = "Pion-Track THnSparse; Centrality; Z vertex;  #pionpT;#pion E; #pion y; #pion Eta; #pion phi;";
    axisNames = axisNames + "track_pT; track Eta; track Phi; #Dphi ; #Deta; #|Deta|; Zt; Xi;";
    axisNames = axisNames + "#pion Mass; ph1_pT; ph2_pT; ph1_eta; ph2_eta; ph1_phi; ph2_phi; ph1_M02; ph2_M02;";

    /////////////Pi0--track correlations////////////
    h_Pi0Track = new THnSparseD("h_Pi0Track", axisNames, nbins_PionCorr, bins, xmin,xmax);
    h_Pi0Track->Sumw2();
    fOutput->Add(h_Pi0Track);

    h_Pi0Track_Mixed = new THnSparseD("h_Pi0Track_Mixed", axisNames, nbins_PionCorr, bins, xmin,xmax);
    h_Pi0Track_Mixed->Sumw2();
    fOutput->Add(h_Pi0Track_Mixed);
    
    //////////////////////////////Cluster-Track correlations:///////////////////////////////////////
    const int nbins_ClusterCorr = 18;
    int binsClusterCorr[nbins_ClusterCorr]    = {nbins_Centrality, nbins_zvertex, nbins_Pt, nbins_E, nbins_eta, nbins_eta, nbins_phi,
                                                nbins_Pt,  nbins_eta, nbins_phi, nbins_dphi, nbins_deta, nbins_deta/2, nbins_zt, nbins_xi,
                                                nbins_M02, nbins_nMatchedTracks, nbins_nMaxima};
                        
    double xminClusterCorr[nbins_ClusterCorr] = {min_Centrality, min_zvertex, min_Pt  , min_E,  min_eta, min_eta, min_phi, 
                                                min_Pt, min_eta, min_phi    , min_dphi   , min_deta   , 0.0, min_zt, min_xi,
                                                min_M02,   min_nMatchedTracks, min_nMaxima};
                                   
    double xmaxClusterCorr[nbins_ClusterCorr] = {max_Centrality, max_zvertex,  max_Pt,  max_E, max_eta    , max_eta, max_phi, 
                                                max_Pt, max_eta,  max_phi   , max_dphi,   max_deta, max_deta, max_zt, max_xi,
                                                max_M02, max_nMatchedTracks, max_nMaxima};

    axisNames = "Cluster-Track THnSparse; Centrality; Z vertex; Cluster p_{T}; Cluster E; Cluster #eta; Cluster y; Cluster #phi;";
    axisNames = axisNames + "track p_{T}; track #eta; track #phi; #Delta#phi ; #Delta#eta; |#Delta#eta|; Z_{T}t; Xi;";
    axisNames = axisNames + "M02; nMatchedTracks; nMaxima;";
     
    h_ClusterTrack = new THnSparseD("h_ClusterTrack", axisNames, nbins_ClusterCorr, binsClusterCorr, xminClusterCorr,xmaxClusterCorr);
    h_ClusterTrack->Sumw2();
    fOutput->Add(h_ClusterTrack);
    
    h_ClusterTrack_Mixed = new THnSparseD("h_ClusterTrack_Mixed", axisNames, nbins_ClusterCorr, binsClusterCorr, xminClusterCorr,xmaxClusterCorr);
    h_ClusterTrack_Mixed->Sumw2();
    fOutput->Add(h_ClusterTrack_Mixed);
    
    ///////////////Pi0////////////////////////////////////
    const int nbins_Pi = 25;
    axisNames = "Pion THnSparse; Centrality; Z vertex ;#pion Mass; #pionpT; #pion Eta; #pion phi; #pion E;";
    axisNames = axisNames + "ph1_E; ph2_E; Asymmetry; ph1_pT; ph2_pT; ph1 #eta; ph2 #eta; ph1 #phi; ph2 #phi; #Delta#phi;";
    axisNames = axisNames + "ph1 #lambda_{02}; ph2 #lambda_{02}; ph1 nCells; ph2 nCells; ph1 nMatchedTracks; ph2 nMatchedTracks; ph1 nMaxima; ph2 nMaxima;";
    int binsPi0[nbins_Pi] = {nbins_Centrality, nbins_zvertex, nbins_Mass, nbins_Pt, nbins_eta, nbins_phi, nbins_E, 
                         nbins_E, nbins_E, nbins_Asymmetry, nbins_Pt, nbins_Pt, nbins_eta, nbins_eta, nbins_phi, nbins_phi, nbins_alpha, 
                         nbins_M02, nbins_M02, nbins_Ncells, nbins_Ncells, nbins_nMatchedTracks, nbins_nMatchedTracks, nbins_nMaxima, nbins_nMaxima};
                            
    double xminPi0[nbins_Pi] = {min_Centrality, min_zvertex, min_Mass, min_Pt, min_eta, min_phi , min_E,
                                min_E, min_E, min_Asymmetry, min_Pt, min_Pt, min_eta, min_eta, min_phi, min_phi, min_alpha,
                                min_M02, min_M02 , min_Ncells, min_Ncells, min_nMatchedTracks, min_nMatchedTracks, min_nMaxima, min_nMaxima};
    double xmaxPi0[nbins_Pi] = {max_Centrality, max_zvertex, max_Mass, max_Pt, max_eta , max_phi , max_E,
                                max_E, max_E, max_Asymmetry, max_Pt, max_Pt, max_eta, max_eta, max_phi, max_phi, max_alpha,
                                max_M02, max_M02, max_Ncells, max_Ncells, max_nMatchedTracks, max_nMatchedTracks, max_nMaxima, max_nMaxima};
                                
    h_Pi0= new THnSparseD("h_Pi0", axisNames, nbins_Pi, binsPi0, xminPi0, xmaxPi0);
    h_Pi0->Sumw2();
    fOutput->Add(h_Pi0);
    
    /////////////////////Clusters////////////////////////////////////
    const int nbins_Cluster = 10;
    
    axisNames = "Cluster THnSparse; Centrality; Z vertex; Cluster E; Cluster p_{T}; Cluster #eta; Cluster #phi; Cluster #lambda_{02}; nCells; nMatchedTracks; nMaxima;";
    int binsCluster[nbins_Cluster] = {nbins_Centrality, nbins_zvertex, nbins_E, nbins_Pt, nbins_eta, nbins_phi, nbins_M02, nbins_Ncells, nbins_nMatchedTracks, nbins_nMaxima};
    double xminCluster[nbins_Cluster] = {min_Centrality, min_zvertex, min_E, min_Pt, min_eta, min_phi, min_M02, min_Ncells, min_nMatchedTracks, min_nMaxima};
    double xmaxCluster[nbins_Cluster] = {max_Centrality, max_zvertex, max_E, max_Pt, max_eta, max_phi, max_M02, max_Ncells, max_nMatchedTracks, max_nMaxima};
    h_Cluster = new THnSparseD("h_Cluster", axisNames, nbins_Cluster, binsCluster, xminCluster, xmaxCluster);
    h_Cluster->Sumw2();
    fOutput->Add(h_Cluster);
    
    ///////////////////Tracks////////////////////////////////////////////////
    axisNames = "Track ThnSparse; Track Pt; Track Eta ; Track Phi;";
    int    binsTrack[3] = {nbins_Pt*2, nbins_eta, nbins_phi};
    double xminTrack[3] = {min_Pt, min_eta, min_phi};
    double xmaxTrack[3] = {max_Pt*2, max_eta, max_phi};
    h_Track = new THnSparseD("h_Track", axisNames, 3, binsTrack, xminTrack, xmaxTrack);
    h_Track->Sumw2();
    fOutput->Add(h_Track);
    
	PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
}



void AliAnalysisTaskEMCALPi0GammaCorr::InitEventMixer()
{
	int nCentBins=fMixBCent->GetNbins();
	double centBins[nCentBins+1];
	centBins[0] = fMixBCent->GetBinLowEdge(1);
	for(int i=1; i<=nCentBins; i++)
	{
		centBins[i] = fMixBCent->GetBinUpEdge(i);
	}

	int nZvtxBins=fMixBZvtx->GetNbins();
	double zvtxbin[nZvtxBins+1];
	zvtxbin[0] = fMixBZvtx->GetBinLowEdge(1);
	for(int i=1; i<=nZvtxBins; i++)
	{
		zvtxbin[i] = fMixBZvtx->GetBinUpEdge(i);
	}
	if(!fPoolMgr)
	{
		cout<<"....  Pool Manager Created ...."<<endl;
		fPoolMgr = new AliEventPoolManager(fPoolSize, fTrackDepth, nCentBins, centBins, nZvtxBins, zvtxbin);
		fPoolMgr->SetTargetValues(fTrackDepth, 0.1, 5);  //pool is ready at 0.1*fTrackDepth = 5000 or events =5
	}
	else
	{
		fPoolMgr->ClearPools();
		cout<<"....  Pool Manager Provided From File ...."<<endl;
	}
	if( (fPoolMgr->GetNumberOfMultBins() != nCentBins) || (fPoolMgr->GetNumberOfZVtxBins() != nZvtxBins) )
	{
		AliFatal("Binning of given pool manager not compatible with binning of correlation task!");
	}

	if(fSavePool==1)
	{
	fPoolMgr->SetSaveFlag(-1, 10000, -10000, 100000, 0, 0, -1, 10000000);
	fOutput->Add(fPoolMgr);
	}
    std::cout << "Ending InitEventMixer" << std::endl;
	fPoolMgr->Validate();
}

void AliAnalysisTaskEMCALPi0GammaCorr::AddEventPoolsToOutput(double minCent, double maxCent,  double minZvtx, double maxZvtx, double minPt, double maxPt)
{
	std::vector<double> binVec;
	binVec.push_back(minCent);
	binVec.push_back(maxCent);
	binVec.push_back(minZvtx);
	binVec.push_back(maxZvtx);
	binVec.push_back(minPt);
	binVec.push_back(maxPt);
	fEventPoolOutputList.push_back(binVec);
}

void AliAnalysisTaskEMCALPi0GammaCorr::ExecOnce()
{
    AliAnalysisTaskEmcal::ExecOnce();
}

Bool_t AliAnalysisTaskEMCALPi0GammaCorr::IsEventSelected()
{
    if (!fEventCuts.AcceptEvent(InputEvent()))
	{
	  PostData(1, fOutput);
	  return kFALSE;
	}
    TString Trigger;
    Trigger = fInputEvent->GetFiredTriggerClasses();
    bool PassedGammaTrigger = kFALSE;
    bool PassedMinBiasTrigger = kFALSE;
    if(Trigger.Contains("EG1") ||Trigger.Contains("EG2") || Trigger.Contains("DG1") || Trigger.Contains("DG2")) PassedGammaTrigger = kTRUE;
    if(Trigger.Contains("INT7")) PassedMinBiasTrigger = kTRUE;
    if(!PassedGammaTrigger && !PassedMinBiasTrigger) return kFALSE;

    bool isSelected = AliAnalysisTaskEmcal::IsEventSelected();
    return kTRUE;
	//return isSelected;
}

Bool_t AliAnalysisTaskEMCALPi0GammaCorr::Run(){
	return kTRUE;
}

Bool_t AliAnalysisTaskEMCALPi0GammaCorr::FillHistograms()
{
	//..This function is called in AliAnalysisTaskEmcal::UserExec.
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
    TString Trigger;
    Trigger = fInputEvent->GetFiredTriggerClasses();
    bool PassedGammaTrigger = kFALSE;
    bool PassedMinBiasTrigger = kFALSE;
    
    if(Trigger.Contains("EG1") ||Trigger.Contains("EG2") || Trigger.Contains("DG1") || Trigger.Contains("DG2")) PassedGammaTrigger = kTRUE;
    if(Trigger.Contains("INT7")) PassedMinBiasTrigger = kTRUE;
    
	double zVertex = fVertex[2];
	AliParticleContainer* tracks =0x0;
	tracks   = GetParticleContainer(0);

    if(PassedGammaTrigger)
	{
        CorrelateClusterAndTrack(tracks,0,1,1);//correlate with same event
    }

	AliEventPool* pool = 0x0;
	pool = fPoolMgr->GetEventPool(fCent, zVertex);
	if (!pool)
	{
		return kFALSE;
	}
	if(pool->IsReady() && PassedGammaTrigger)
	{
		int nMix = pool->GetCurrentNEvents();
		for(int jMix=0; jMix<nMix; jMix++)
		{
			TObjArray* bgTracks=0x0;
			bgTracks = pool->GetEvent(jMix);
			if(!bgTracks)
			{
				cout<<"could not retrieve TObjArray from EventPool!"<<endl;
			}
			CorrelateClusterAndTrack(0,bgTracks,0,1.0/nMix);//correlate with mixed event
		}
	}
    if(PassedMinBiasTrigger && !PassedGammaTrigger )
	{
        TObjArray* tracksClone=0x0;
		tracksClone = CloneToCreateTObjArray(tracks);
		if(pool && tracksClone && !pool->GetLockFlag())
		{
          pool->UpdatePool(tracksClone);
		}
	}
   
}

TObjArray* AliAnalysisTaskEMCALPi0GammaCorr::CloneToCreateTObjArray(AliParticleContainer* tracks)
{
	//..clones a track list
	if(!tracks)                            return 0;
	if(tracks->GetNAcceptedParticles()==0) return 0;
	TObjArray* tracksClone = new TObjArray;
	tracksClone->SetOwner(kTRUE);
	int NoOfTracksInEvent =tracks->GetNParticles();
	AliVParticle* track=0;
	for(int NoTrack = 0; NoTrack < NoOfTracksInEvent; NoTrack++)
	{
		track = (AliVParticle*)tracks->GetAcceptParticle(NoTrack);
		if(!track)continue; 
		tracksClone->Add(new AliPicoTrack(track->Pt(), track->Eta(), track->Phi(), track->Charge(), 0, 0, 0, 0));
	}
	if(tracksClone->GetEntries()!=tracks->GetNAcceptedParticles())cout<<"!!!!!!! Major error!!!! "<<"Accepted tracks in event: "<<tracks->GetNAcceptedParticles()<<", Tracks in TObjArray: "<<tracksClone->GetEntries()<<endl;
	return tracksClone;
}

int AliAnalysisTaskEMCALPi0GammaCorr::CorrelateClusterAndTrack(AliParticleContainer* tracks,TObjArray* bgTracksArray,Bool_t SameMix, double InputWeight)
{
	AliClusterContainer* clusters  = GetClusterContainer(0);  
	if (!clusters) return 0;
	int NoOfClustersInEvent =clusters->GetNClusters();
	double EffWeight_Gamma;
	double EffWeight_Hadron;
	double Weight=1.0;    
	AliVCluster* cluster = 0;
	AliVCluster* cluster2= 0;
	AliVParticle* track=0;
    Weight=InputWeight; //..for mixed events normalize per events in pool

    if(SameMix!=0){ //if same event, then fill track histo
        for(int NoTrack = 0; NoTrack < tracks->GetNParticles(); NoTrack++){ //correlate pion with tracks
	        track = (AliVParticle*)tracks->GetAcceptParticle(NoTrack);
            if(!track) continue;
            if(track->Pt()<0.5) continue;
            double entries[3] = {track->Pt(), track->Eta(), TVector2::Phi_mpi_pi(track->Phi())};
            h_Track->Fill(entries);
        }    
    }
    
	for(int NoCluster1 = 0; NoCluster1 < NoOfClustersInEvent; NoCluster1++ ) // Loop over clusters
	{   
        cluster=(AliVCluster*) clusters->GetCluster(NoCluster1); // //it was GetAcceptCluster->GetCluster(NoCluster1);
     	if(!cluster) continue;
        if(!PassedCuts(cluster))continue ; 
        std::cout << "number of tracks matched to cluster" << cluster->GetNTracksMatched() << std::endl;
        if(SameMix==0){
            
		    for( int ibg=0; ibg<bgTracksArray->GetEntries(); ibg++){//correlate cluster with tracks
	    	    AliPicoTrack* track = static_cast<AliPicoTrack*>(bgTracksArray->At(ibg));
		        if(!track) continue;
                FillPhotonCorrelation(cluster, track, h_ClusterTrack_Mixed, Weight);
            }
        }
        else{
            FillClusterHisto(cluster, h_Cluster); //filling photon histogram
        
            for(int NoTrack = 0; NoTrack < tracks->GetNParticles(); NoTrack++){ //correlate cluster with tracks
			    track = (AliVParticle*)tracks->GetAcceptParticle(NoTrack);
                if(!track) continue;
                FillPhotonCorrelation(cluster, track, h_ClusterTrack, Weight);
             }
        }
        
        for( int NoCluster2 = NoCluster1+1; NoCluster2 < NoOfClustersInEvent; NoCluster2++ )
        {
            cluster2=(AliVCluster*) clusters->GetCluster(NoCluster2);
			if(!cluster2) continue;
            if(!PassedCuts(cluster2))continue ; 
            FillPionHisto(cluster, cluster2, h_Pi0); //filling Pion histogram
            if(SameMix==0){
		        for( int ibg=0; ibg<bgTracksArray->GetEntries(); ibg++){
	    	        AliPicoTrack* track = static_cast<AliPicoTrack*>(bgTracksArray->At(ibg));
		            if(!track) continue;
                    FillPionCorrelation(cluster, cluster2, track, h_Pi0Track_Mixed, Weight);
                }
            }
            else{
                for(int NoTrack = 0; NoTrack < tracks->GetNParticles(); NoTrack++){ //correlate pion with tracks
				    track = (AliVParticle*)tracks->GetAcceptParticle(NoTrack);
                    if(!track) continue;
                    FillPionCorrelation(cluster, cluster2, track, h_Pi0Track, Weight);
                    }
            }
        } //end 2 loop over clusters
	} //end  1 loop over clusters
	return 1; //return number of accepted gammas
}


double AliAnalysisTaskEMCALPi0GammaCorr::GetIsolation_Track(AliVCluster* cluster){
}

void  AliAnalysisTaskEMCALPi0GammaCorr::FillPionCorrelation(AliVCluster* cluster1, AliVCluster* cluster2, AliVParticle* track, THnSparse* histo, double weight){

    AliClusterContainer* clusters  = GetClusterContainer(0);
    
    AliVCluster* cluster_lead = 0;
    AliVCluster* cluster_sub  = 0;
    
    TLorentzVector ph_lead, ph_sub, pi0; 
       
    if(cluster1->E() > cluster2->E()){
        cluster_lead = cluster1;
        cluster_sub  = cluster2;
    }
    else{
        cluster_lead = cluster2;
        cluster_sub  = cluster1; 
    }    
    
	clusters->GetMomentum(ph_lead, cluster_lead);
    clusters->GetMomentum(ph_sub, cluster_sub);
    double asym = abs(ph_lead.E()-ph_sub.E())/(ph_lead.E()+ph_sub.E());
	pi0= ph_lead+ph_sub;
    
    //////////////////Selection//////////////////////////////////////////////
    if( cluster_lead->E()<6.0) return; //at least one photon with 6 GeV of energy, this is lowest threshold trigger in pPb data.
    if( pi0.Pt() < 6.0 ) return;
    if( pi0.M()  > 1.0 ) return;
    if( track->Pt()<0.5 ) return;
    //if( asym > 0.7 ) return;
    //if( cluster1->GetM02()>0.4 || cluster2->GetM02()>0.4 ) return;
    /////////////////////////////////////////////////////////////////////////
    
    double deta = pi0.Eta()-track->Eta();
    double  Zt  = track->Pt()/pi0.Pt();
    double  Xi  = -999; 
    if(Zt>0) Xi = TMath::Log(1.0/Zt);
    
    double trackphi = TVector2::Phi_mpi_pi(track->Phi());
    double dphi;
    
    //filling
    dphi     = TVector2::Phi_mpi_pi(pi0.Phi()- trackphi)/TMath::Pi();
    if(dphi<-0.5) dphi +=2;
    
    double entries[24] = {fCent, fVertex[2], pi0.Pt(), pi0.E(), pi0.Rapidity(), pi0.Eta(), pi0.Phi(), 
                         track->Pt(), track->Eta(), trackphi, dphi, deta, abs(deta), Zt, Xi,
                         pi0.M(), ph_lead.Pt(), ph_sub.Pt(), ph_lead.Eta(), ph_sub.Eta(), ph_lead.Phi(), ph_sub.Phi() , cluster_lead->GetM02(), cluster_sub->GetM02()};                
    histo->Fill(entries, weight); //
    /*
    entries[11] = -1.0*deta;           
    histo->Fill(entries, weight); //
    
    dphi     = -1.0*TVector2::Phi_mpi_pi(pi0.Phi()- trackphi)/TMath::Pi();
    if(dphi<-0.5) dphi +=2;
    
    entries[10] =  dphi; 
    entries[11] =  deta;            
    histo->Fill(entries, weight); //
     
    entries[11] = -1.0*deta;            
    histo->Fill(entries, weight); //
    */
    return;
}

void  AliAnalysisTaskEMCALPi0GammaCorr::FillPhotonCorrelation(AliVCluster* cluster, AliVParticle* track, THnSparse* histo, double weight){
    AliClusterContainer* clusters  = GetClusterContainer(0);
    TLorentzVector ph;
	clusters->GetMomentum(ph, cluster);
   
    if( track->Pt()<0.5 ) return;
    if( cluster->E()< 6.0) return;
    //if( cluster1->GetM02()>0.4 ) return
    /////////////////////////////////////////////////////////////////////////
    
    
    //    axisNames = "Cluster-Track THnSparse; Centrality; Z vertex; Cluster pT; Cluster E; Cluster Eta; Cluster y; Cluster phi;";
    //axisNames = axisNames + "track_pT; track Eta; track #phi; #Dphi ; #Delta#eta; #|Deta|; Zt; Xi;";
    //axisNames = axisNames + "M02;";
    
    double trackphi = TVector2::Phi_mpi_pi(track->Phi());
    double dphi;
    
    double deta = ph.Eta()-track->Eta();
    double  Zt  = track->Pt()/ph.Pt();
    double  Xi  = -999; 
    if(Zt>0) Xi = TMath::Log(1.0/Zt);
    
    dphi = TVector2::Phi_mpi_pi(ph.Phi()- trackphi)/TMath::Pi();
    if(dphi<-0.5) dphi +=2;
    
    double entries[18] = {fCent, fVertex[2], ph.Pt(), ph.E(), ph.Eta(), ph.Rapidity(), ph.Phi(),  
              track->Pt(), track->Eta(), trackphi, dphi, deta, abs(deta), Zt, Xi,
              cluster->GetM02(), cluster->GetNTracksMatched(), cluster->GetNExMax()
              };                
    histo->Fill(entries, weight);//
    
    /*entries[11] = -1.0*deta;
    histo->Fill(entries, weight);//
    
    dphi = -1.0*TVector2::Phi_mpi_pi(ph.Phi()- trackphi)/TMath::Pi();
    if(dphi<-0.5) dphi +=2;
    
    entries[10] = dphi;
    entries[11] = deta;
    histo->Fill(entries, weight);//
    
    entries[11] = -1.0*deta;
    histo->Fill(entries, weight);//
    */
    
    return;
}





void  AliAnalysisTaskEMCALPi0GammaCorr::FillPionHisto(AliVCluster* cluster1, AliVCluster* cluster2, THnSparse* histo){
    
    AliClusterContainer* clusters  = GetClusterContainer(0);
    AliVCluster* cluster_lead = 0;
    AliVCluster* cluster_sub  = 0;
    
    TLorentzVector ph_lead, ph_sub, pi0; 
       
    if(cluster1->E() > cluster2->E()){
        cluster_lead = cluster1;
        cluster_sub  = cluster2;
    }
    else{
        cluster_lead = cluster2;
        cluster_sub  = cluster1; 
    }    
    
    clusters->GetMomentum(ph_lead, cluster_lead);
    clusters->GetMomentum(ph_sub,  cluster_sub);
    
	pi0 = ph_lead + ph_sub;
    //////////////////Selection/////////////////////////////////////////
    if( pi0.Pt() < 6.0) return;
    if( pi0.M()  > 1.0) return;
    if(cluster_lead->E()<6.0) return;
    ////////////////////////////////////////////////////////////////////
    double asym = abs(ph_lead.E()-ph_sub.E())/(ph_lead.E()+ph_sub.E());
    double entries[25] = {fCent, fVertex[2], pi0.M(), pi0.Pt(), pi0.Eta(), pi0.Phi(), pi0.E(), 
                          ph_lead.E(), ph_sub.E(), asym, ph_lead.Pt(), ph_sub.Pt(), ph_lead.Eta(), ph_sub.Eta(), ph_lead.Phi(), ph_sub.Phi() , abs(TVector2::Phi_mpi_pi(ph_lead.Phi()-ph_sub.Phi())), 
                          cluster_lead->GetM02(), cluster_sub->GetM02(), cluster_lead->GetNCells(), cluster_sub->GetNCells(), cluster_lead->GetNTracksMatched(), cluster_sub->GetNTracksMatched(), 
                          cluster_lead->GetNExMax(), cluster_sub->GetNExMax()};                
    histo->Fill(entries);
    return;
}

void AliAnalysisTaskEMCALPi0GammaCorr::FillClusterHisto(AliVCluster* cluster, THnSparse* histo){
    
    AliClusterContainer* clusters  = GetClusterContainer(0);
    TLorentzVector ph;
    clusters->GetMomentum(ph, cluster);
    if(cluster->E()< 5.0) return;
    double entries[10] = {fCent, fVertex[2], ph.E(), ph.Pt(), ph.Eta(), ph.Phi(), cluster->GetM02(), cluster->GetNCells(), cluster->GetNTracksMatched(), cluster->GetNExMax()};                
    histo->Fill(entries);
    return;
}

TObjArray* AliAnalysisTaskEMCALPi0GammaCorr::CloneClustersTObjArray(AliClusterContainer* clusters)
{
	if(!clusters)                            return 0;
	if(clusters->GetNClusters()==0) return 0;
	TObjArray* clustersCloneI = new TObjArray;
	clustersCloneI->SetOwner(kTRUE);
	int NoOfClustersInEvent =clusters->GetNClusters();
	AliVCluster* cluster = 0;
	for(int NoClus = 0; NoClus < NoOfClustersInEvent; NoClus++)
	{
		cluster = (AliVCluster*) clusters->GetAcceptCluster(NoClus);
		if(!cluster)continue; //check if the Cluster is good
		clustersCloneI->Add((AliVCluster*)cluster);// Add(new AliPicoTrack(track->Pt(), track->Eta(), track->Phi(), track->Charge(), 0, 0, 0, 0));
	}
	if(clustersCloneI->GetEntries()!=clusters->GetNAcceptedClusters())cout<<"!!!!!!! Major error!!!! "<<"Accepted clusters in event: "<<clusters->GetNAcceptedClusters()<<", Tracks in TObjArray: "<<clustersCloneI->GetEntries()<<endl;
	return clustersCloneI;
}


Bool_t AliAnalysisTaskEMCALPi0GammaCorr::PassedCuts(AliVCluster* cluster)
{
    if(!cluster->IsEMCAL()) return kFALSE;
    if(cluster->E()<3.0) return kFALSE;
	//if(cluster->GetNCells()<2) return kFALSE;
	//if(cluster->GetNExMax() > 1) return kFALSE; //local maxima should be 0 or 1
	//if(cluster->GetM02()<0.1) return kFALSE;
	//if(fRmvMTrack==1 && caloCluster->GetNTracksMatched(s)!=0) return kFALSE;
	return kTRUE;
}


double AliAnalysisTaskEMCALPi0GammaCorr::GetEff(AliTLorentzVector ClusterVec)
{
	return 1;
}
