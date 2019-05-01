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

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"

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
#include "AliEMCALGeometry.h"
#include "AliEMCALGeoParams.h"
#include "AliESDInputHandler.h"
#include "AliMCParticleContainer.h"

#include <AliGenPythiaEventHeader.h>


#include "TDatabasePDG.h"

#include <memory>
using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskEMCALPi0GammaCorr)

////////////////////////////////////////////////////////////////////////////////////////
AliAnalysisTaskEMCALPi0GammaCorr::AliAnalysisTaskEMCALPi0GammaCorr():
AliAnalysisTaskEmcal("AliAnalysisTaskEMCALPi0GammaCorr", kTRUE),
fIsMC(kFALSE),
fAODMCParticles(0),
fmcHeader(0),
fSavePool(0),
fFiducialCellCut(0x0),
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
h_TrackITS(0),
h_Truth(0),
h_Cluster(0),
h_ClusterTrack(0),
h_ClusterTrack_Mixed(0),
h_Pi0(0),
h_Pi0Track(0),
h_Pi0Track_Mixed(0),
h_nEvents(0),
fPeriod("")
{
    InitArrays();
}

// -------------------------------------------------------------------------------------
// Constructor with inputs
AliAnalysisTaskEMCALPi0GammaCorr::AliAnalysisTaskEMCALPi0GammaCorr(Bool_t InputDoMixing):
AliAnalysisTaskEmcal("AliAnalysisTaskEMCALPi0GammaCorr", kTRUE),
fIsMC(kFALSE),
fAODMCParticles(0),
fmcHeader(0),
fSavePool(0),
fFiducialCellCut(0x0),
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
h_TrackITS(0),
h_Truth(0),
h_Cluster(0),
h_ClusterTrack(0),
h_ClusterTrack_Mixed(0),
h_Pi0(0),
h_Pi0Track(0),
h_Pi0Track_Mixed(0),
h_nEvents(0),
fPeriod("")
{
	InitArrays();
}//End constructor PiHadron that receives input

void AliAnalysisTaskEMCALPi0GammaCorr::InitArrays()
{
    AliWarning("InitArrays is being called");
    fSavePool          =0; //= 0 do not save the pool by default. Use the set function to do this.
    
    //Setting bins for the mixing of events.
    double centmix[kNcentBins+1] = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 80.0, 100.0};
    fMixBCent = new TAxis(kNcentBins,centmix);

    double zvtxmix[kNvertBins+1] = {-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10};
    memcpy (fArrayNVertBins, zvtxmix, sizeof (fArrayNVertBins));
    fMixBZvtx = new TAxis(kNvertBins,zvtxmix);
    fTrackDepth     = 50000;    //Raymonds/Megans value
    fPoolSize       = 1; 
    //SetMakeGeneralHistograms(kTRUE);

    fFiducialCellCut = new AliEMCALRecoUtils(); 
} //end of function init arrays

AliAnalysisTaskEMCALPi0GammaCorr::~AliAnalysisTaskEMCALPi0GammaCorr()
{

}

void AliAnalysisTaskEMCALPi0GammaCorr::UserCreateOutputObjects()
{
    AliWarning("Entering UserCreateOutPutObjects");   
    AliAnalysisTaskEmcal::UserCreateOutputObjects();
       
    AliWarning("Initializing Event Mixer");
    InitEventMixer();
	
    //Initializing the histograms to be saved. For the moment, only pT of clusters and Mgammagamma.
   
    
    

    int nbins_Mass = 150;
    int nbins_Mass_corr = 30;
    double min_Mass = 0.000;
    double max_Mass = 0.300;

    int    nbins_Pt =  50;
    double min_Pt   =  0.0;
    double max_Pt   =  50.0;

    int    nbins_TrackPt = 10;
    double min_TrackPt = 0.0; 
    double max_TrackPt = 10.0;

    int nbins_dphi     = 18;
    int nbins_phi     = 140;
    double min_dphi = -0.5; // rads
    double max_dphi = 1.5; // rads 
    double min_phi = 1.0;
    double max_phi = TMath::Pi();

    int nbins_zvertex = 25; //it as 5
    double min_zvertex = -10;
    double max_zvertex = +10;

    int nbins_eta = 100;
    int nbins_deta= 30;
    double max_deta = 1.5;
    double max_eta =  0.75;
    double min_deta = -max_deta;
    double min_eta =  -max_eta;
    
    int nbins_xi      = 10;   
    double min_xi = 0;
    double max_xi = 5.0;
    
    int nbins_zt      = 12;
    double min_zt = 0;
    double max_zt = 1.2;
   

    int nbins_E_Corr = 10;
    double min_E_Corr = 0.0;
    double max_E_Corr = 20.0;

    int nbins_M02_Corr     = 10;
    double min_M02_Corr = 0.0;
    double max_M02_Corr = 2.0;

    int nbins_M02 = 50;
    double min_M02 = 0.0;
    double max_M02 = 2.0;
    
    int nbins_Ncells  = 30;
    double min_Ncells = -0.5;
    double max_Ncells = 29.5;
    
    int nbins_Centrality = 5; 
    double min_Centrality = 0;
    double max_Centrality = 100;
    
    int nbins_Asymmetry = 10;
    double min_Asymmetry =0;
    double max_Asymmetry = 1.0;
    
    int    nbins_nMaxima = 3;
    double min_nMaxima = -0.5;
    double max_nMaxima  =  2.5;
    
    int nbins_alpha = 50;
    double min_alpha = 0.0;
    double max_alpha = 50;

    
    int nbins_DisToBorder = 6;
    double min_DisToBorder = -0.5;
    double max_DisToBorder = 5.5;

    int nbins_DisToBad = 6;
    double min_DisToBad = -0.5;
    double max_DisToBad = 5.5;    

    int nbins_Exoticity  = 100;
    double min_Exoticity = 0.0;
    double max_Exoticity = 1.0;

    int nbins_time =40;
    double min_time = -40.0;
    double max_time = +40.0;

    int nbins_RunNumber = 20;
    double min_RunNumber = -0.5;
    double max_RunNumber = 19.5;
    
    //int nbins_BCID     = 3600;
    //double min_BCID    = -0.5;
    //double max_BCID    = 3599.5;

    int nbins_IsoE  = 30;
    double min_IsoE = -5.0;
    double max_IsoE = 10.0; 

    int nbins_IsoETruth  = 40; 
    double min_IsoETruth = 0.0;
    double max_IsoETruth = 20.0;
 
    int nbins_nTracks = 50;
    double min_nTracks = 0.0;
    double max_nTracks = 100.0;
   
    int nbins_nClusters = 50;
    double min_nClusters = 0.0;
    double max_nClusters = 100.0;
    
    int nbins_dR = 10;
    double min_dR = 0.0;
    double max_dR = 0.05;

    int nbins_Matching = 10;
    double min_Matching = 0.00; 
    double max_Matching = 0.05; 

    int nbins_trueGamma = 2;
    double min_trueGamma = 0.0;
    double max_trueGamma = 1.0;

    int nbins_truePDG  = 2;
    double min_truePDG = 0.0; 
    double max_truePDG = 2.0;

    int nbins_TruePt  = 40;
    double min_TruePt  = 0.0;
    double max_TruePt  = 20; 
    
    int nbins_TrueEta = 40; 
    double min_TrueEta   = -1.0;
    double max_TrueEta   = 1.0;


    int nbins_d0 = 200;
    double min_d0 = -3000.0;
    double max_d0 = 3000.0;
  
    int nbins_z0 = 100;
    double min_z0   = -20.0;
    double max_z0   = +20.0;

    //////////////////////Pion-hadron correlations//////////////////////////
    const int nbins_PionCorr = 13;
    int bins[nbins_PionCorr]    = {nbins_Centrality, nbins_zvertex, nbins_Pt,  //trigger variables
                                   nbins_TrackPt, nbins_dphi, nbins_deta, nbins_zt, nbins_xi, //track variables
                                   nbins_Mass_corr,  nbins_E_Corr, nbins_E_Corr, nbins_M02_Corr, nbins_M02_Corr}; //pion-only and pion-decay variables
                                   
    double xmin[nbins_PionCorr] = {min_Centrality, min_zvertex, min_Pt ,   
                                   min_TrackPt, min_dphi   , min_deta , min_zt, min_xi,
                                   min_Mass, min_E_Corr, min_E_Corr, min_M02_Corr, min_M02_Corr};
                                   
    double xmax[nbins_PionCorr] = {max_Centrality, max_zvertex,  max_Pt,   
                                   max_TrackPt, max_dphi,   max_deta,  max_zt, max_xi,
                                   max_Mass, max_E_Corr, max_E_Corr,  max_M02_Corr, max_M02_Corr};

    TString axisNames = "Pion--Track THnSparse; Centrality[%]; Z vertex [cm];  #pionpT [GeV];";
    axisNames = axisNames + "track_pT [GeV]; #Dphi  [rad]; #Deta; Zt; Xi;";
    axisNames = axisNames + "#pion Mass; ph1_pT; ph2_pT; ph1_M02; ph2_M02;";

    h_Pi0Track = new THnSparseD("h_Pi0Track", axisNames, nbins_PionCorr, bins, xmin,xmax);
    h_Pi0Track->Sumw2();
    //fOutput->Add(h_Pi0Track);

    h_Pi0Track_Mixed = new THnSparseD("h_Pi0Track_Mixed", axisNames, nbins_PionCorr, bins, xmin,xmax);
    h_Pi0Track_Mixed->Sumw2();
    //fOutput->Add(h_Pi0Track_Mixed);
    
    //////////////////////////////Cluster-Track correlations:///////////////////////////////////////
    const int nbins_ClusterCorr = 10;
    int binsClusterCorr[nbins_ClusterCorr]    = {nbins_Centrality, nbins_zvertex, nbins_Pt, 
                                                nbins_TrackPt, nbins_dphi, nbins_deta, nbins_zt, nbins_xi,
						 nbins_M02_Corr, nbins_dR};
                        
    double xminClusterCorr[nbins_ClusterCorr] = {min_Centrality, min_zvertex, min_Pt,  
                                                min_TrackPt, min_dphi   , min_deta   , min_zt, min_xi,
						 min_M02_Corr, min_dR};
                                   
    double xmaxClusterCorr[nbins_ClusterCorr] = {max_Centrality, max_zvertex,  max_Pt,  
                                                 max_TrackPt, max_dphi, max_deta,  max_zt, max_xi,
						 max_M02_Corr, max_dR};

    axisNames = "Cluster-Track THnSparse; Centrality [%]; Z vertex [cm]; Cluster p_{T} [GeV]; ";
    axisNames = axisNames + "track p_{T} [GeV]; #Delta#phi [rad] ; #Delta#eta; Z_{T}t; Xi;"; //track variables
    axisNames = axisNames + "#lamda0; dR;";
     
    h_ClusterTrack = new THnSparseD("h_ClusterTrack", axisNames, nbins_ClusterCorr, binsClusterCorr, xminClusterCorr,xmaxClusterCorr);
    h_ClusterTrack->Sumw2();
    //fOutput->Add(h_ClusterTrack);
    
    h_ClusterTrack_Mixed = new THnSparseD("h_ClusterTrack_Mixed", axisNames, nbins_ClusterCorr, binsClusterCorr, xminClusterCorr,xmaxClusterCorr);
    h_ClusterTrack_Mixed->Sumw2();
    //fOutput->Add(h_ClusterTrack_Mixed);
    
    ///////////////Pi0////////////////////////////////////
    const int nbins_Pi = 15;
    axisNames = "Pion THnSparse; Centrality [%]; Z vertex [cm];#pi Mass [GeV]; #pi pT [GeV]; #pi y;"; 
    axisNames = axisNames+ "Asymmetry; ph1_pT [GeV]; ph2_pT [GeV]; #Delta#phi [mrad];";
    axisNames = axisNames+ "ph1 #lambda_{02}; ph2 #lambda_{02}; ";
    axisNames = axisNames+ "ph1 dR; ph2 dR; true pT; isGammaTriggered;";

    int binsPi0[nbins_Pi] = {nbins_Centrality, nbins_zvertex, nbins_Mass, nbins_Pt, nbins_eta,   
                             nbins_Asymmetry, nbins_Pt, nbins_Pt, nbins_alpha, 
			     nbins_M02, nbins_M02, nbins_dR,nbins_dR,
                             nbins_Pt,2
			     };
                            
    double xminPi0[nbins_Pi] = {min_Centrality, min_zvertex, min_Mass, min_Pt, min_eta, 
                                min_Asymmetry, min_Pt, min_Pt, min_alpha,
                                min_M02, min_M02,  min_dR, min_dR,
                                min_Pt, -0.5
                                };

    double xmaxPi0[nbins_Pi] = {max_Centrality, max_zvertex, max_Mass, max_Pt, max_eta,
                                max_Asymmetry, max_Pt, max_Pt, max_alpha,
                                max_M02, max_M02, max_dR, max_dR,
                                max_Pt, 1.5};
                                
    h_Pi0= new THnSparseD("h_Pi0", axisNames, nbins_Pi, binsPi0, xminPi0, xmaxPi0);
    h_Pi0->Sumw2();
    fOutput->Add(h_Pi0);
    
    /////////////////////Clusters////////////////////////////////////
    const int nbins_Cluster = 26;
    
    axisNames = "Cluster THnSparse; RunNumber; Centrality; Z vertex; Cluster p_{T}; Cluster #eta; Cluster #phi; Cluster #lambda_{02}; nCells; nMaxima;";
    axisNames = axisNames + "Distance to Border; Distance to Bad Cell; dR to track;  d#eta to track; d#phi to track; Exoticity; time [ns]; nTracks ; nClusters;";
    axisNames = axisNames + "ISO Track; UE_track (etaband); UE_track (ortho); ISO Track Subtracted (etaband); ISO Track Subtracted (ortho);"; 
    axisNames = axisNames + " ISO Truth; IsTrueGamma; IsGammaTriggered;" ; 
    int binsCluster[nbins_Cluster] = {nbins_RunNumber, nbins_Centrality, nbins_zvertex, nbins_Pt, nbins_eta, nbins_phi, nbins_M02, nbins_Ncells,  
                                      nbins_nMaxima, nbins_DisToBorder, nbins_DisToBad, nbins_dR, nbins_Matching, nbins_Matching, nbins_Exoticity, nbins_time, nbins_nTracks, nbins_nClusters,
                                      nbins_IsoE, nbins_IsoE, nbins_IsoE, nbins_IsoE, nbins_IsoE, nbins_IsoETruth, nbins_trueGamma, 2};
    double xminCluster[nbins_Cluster] = {min_RunNumber, min_Centrality, min_zvertex, min_Pt,  min_eta, min_phi, min_M02, min_Ncells, min_nMaxima, 
					 min_DisToBorder, min_DisToBad, min_dR, min_Matching, min_Matching, min_Exoticity, min_time, min_nTracks, min_nClusters, 
                                         min_IsoE, min_IsoE, min_IsoE, min_IsoE,  min_IsoE, min_IsoETruth, min_trueGamma, -0.5};
    double xmaxCluster[nbins_Cluster] = {max_RunNumber, max_Centrality, max_zvertex, max_Pt, max_eta, max_phi, max_M02, max_Ncells, max_nMaxima, 
					 max_DisToBorder, max_DisToBad, max_dR, max_Matching, max_Matching, max_Exoticity, max_time, max_nTracks, max_nClusters, 
                                         max_IsoE, max_IsoE, max_IsoE, max_IsoE, max_IsoE, max_IsoETruth, max_trueGamma, 1.5};

    h_Cluster = new THnSparseD("h_Cluster", axisNames, nbins_Cluster, binsCluster, xminCluster, xmaxCluster);
    h_Cluster->Sumw2();
    fOutput->Add(h_Cluster);
    

    int nbins_nTPCclusters = 180;
    double min_nTPCclusters = 0.0;
    double max_nTPCclusters = 180;

    int nbins_TPCsignal = 200;
    double min_TPCsignal  =0.0; 
    double max_TPCsignal = 150;

    int nbins_fTPC = 180;
    double min_fTPC =0.0;
    double max_fTPC = 180;   
    ///////////////////Tracks////////////////////////////////////////////////


    axisNames = "Track ThnSparse; Track Pt; Track Eta ; Track Phi; Track d0 [um]; Track z0 [cm]; # of TPC clusters; TPC signal; Crossed rows; Charge; isGammaTriggeredEvent;";
    int    binsTrack[10] = {5*nbins_Pt, nbins_eta, nbins_phi, nbins_d0, nbins_z0, nbins_nTPCclusters, nbins_TPCsignal, nbins_fTPC, 2,2};
    double xminTrack[10] = {0.0, -1.0, 0.0,           min_d0, min_z0, min_nTPCclusters, min_TPCsignal, min_fTPC, -2.0, -0.5};
    double xmaxTrack[10] = {20.0,  1.0, 2*TMath::Pi(), max_d0, max_z0, max_nTPCclusters, max_TPCsignal, max_fTPC, 2.0, 1.5};
    h_Track = new THnSparseD("h_Track", axisNames, 10, binsTrack, xminTrack, xmaxTrack);
    h_Track->Sumw2();
    fOutput->Add(h_Track); 


    axisNames = "Track ITS ThnSparse; Track Pt; Track Eta ; Track Phi; Track d0 [um]; Track z0 [cm]; Charge; isGammaTriggeredEvent;";
    int    binsTrackITS[7] = {5*nbins_Pt, nbins_eta, nbins_phi, nbins_d0, nbins_z0,  2,2};
    double xminTrackITS[7] = {0.0, -1.0, 0.0,           min_d0/10.0, min_z0/100.0, -2.0, -0.5};
    double xmaxTrackITS[7] = {20.0,  1.0, 2*TMath::Pi(), max_d0/10.0, max_z0/100.0, 2.0, 1.5};
    h_TrackITS = new THnSparseD("h_TrackITS", axisNames, 7, binsTrackITS, xminTrackITS, xmaxTrackITS);
    h_TrackITS->Sumw2();
    fOutput->Add(h_TrackITS);
   
    ///////////////////////////////////TRUTH //////////////////////////////////////////////////////
    axisNames = "Truth ThnSparse; True Pt; True y ; PDG;";
    int    binsTruth[3] = {nbins_TruePt, nbins_TrueEta, nbins_truePDG};
    double xminTruth[3] = {min_TruePt, min_TrueEta, min_truePDG};
    double xmaxTruth[3] = {max_TruePt, max_TrueEta, max_truePDG};
    h_Truth = new THnSparseD("h_Truth", axisNames, 3, binsTruth, xminTruth, xmaxTruth);
    h_Truth->Sumw2();
    fOutput->Add(h_Truth);

    h_nEvents = new TH1F("h_nEvents", "Number of triggers",  2,  -0.5,  1.5);  
    fOutput->Add(h_nEvents);

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

Bool_t AliAnalysisTaskEMCALPi0GammaCorr::PassedGATrigger(){

  TString Trigger = fInputEvent->GetFiredTriggerClasses();
  bool PassedGammaTrigger = kFALSE;
  if(Trigger.Contains("EG1") ||Trigger.Contains("EG2") || Trigger.Contains("DG1") || Trigger.Contains("DG2")) PassedGammaTrigger = kTRUE;
  return PassedGammaTrigger;

}


Bool_t AliAnalysisTaskEMCALPi0GammaCorr::IsEventSelected()
{
    
   
    TString Trigger;
    Trigger = fInputEvent->GetFiredTriggerClasses();
    bool PassedGammaTrigger = kFALSE;
    bool PassedMinBiasTrigger = kFALSE;
    if(Trigger.Contains("EG1") ||Trigger.Contains("EG2") || Trigger.Contains("DG1") || Trigger.Contains("DG2")) PassedGammaTrigger = kTRUE;
    if(Trigger.Contains("INT7")) PassedMinBiasTrigger = kTRUE;
   
    if(!PassedGammaTrigger && !PassedMinBiasTrigger && !fIsMC) return kFALSE; //if not MC and does not trigger data, remove

    bool isSelected = AliAnalysisTaskEmcal::IsEventSelected();
    if(isSelected){
        if(PassedGammaTrigger) h_nEvents->Fill(1.0);
        if(PassedMinBiasTrigger) h_nEvents->Fill(0.0);
    }
    //return kTRUE;


    
    return isSelected;
}

Bool_t AliAnalysisTaskEMCALPi0GammaCorr::Run(){
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


  /*

  AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
  if (!eventHandler) {
      AliFatal("You asked for MC analysis, but I don't find any MCEventHandler... did you forget to add it to your analysis manager?");
  }
  
  AliMCEvent* mc_truth_event = NULL;
  if (eventHandler) mc_truth_event = eventHandler->MCEvent();
  if (!mc_truth_event) AliFatal("Missing MC event");

  if (mc_truth_event != NULL) {
    mc_truth_event->PreReadAll();
  }

  AliGenEventHeader *mc_truth_header = mc_truth_event != NULL ?
      mc_truth_event->GenEventHeader() : NULL;


  AliGenPythiaEventHeader *mc_truth_pythia_header;

  if (mc_truth_header != NULL) {
    //std::cout << " Event WEIGHT " << mc_truth_header->EventWeight() << std::endl;
    mc_truth_pythia_header =	  dynamic_cast<AliGenPythiaEventHeader *>(mc_truth_header);
    if (mc_truth_pythia_header != NULL) {
      //std::cout << "Process Type " <<  mc_truth_pythia_header->ProcessType() << std::endl;
      //std::cout << "PTHARD " << mc_truth_pythia_header->GetPtHard() << std::endl;
      //std::cout << "Xsection " << mc_truth_pythia_header->GetXsection() << std::endl;
      //std::cout << "Trials " << mc_truth_pythia_header->Trials() << std::endl;   
    }
  }

  //AliStack *stack;
  //if (mc_truth_event != NULL) {
  //  stack = mc_truth_event->Stack();
  // }
  
  fMCEvent = mc_truth_event; 
  */
  return kTRUE;
}

Double_t AliAnalysisTaskEMCALPi0GammaCorr::GetCrossEnergy(const AliVCluster *cluster, Short_t &idmax)
{
  // Calculate the energy of cross cells around the leading cell.

  AliVCaloCells *cells = fCaloCells;
  if (!cells)
    return 0;

  if (!fGeom)
    return 0;

  Int_t iSupMod = -1;
  Int_t iTower  = -1;
  Int_t iIphi   = -1;
  Int_t iIeta   = -1;
  Int_t iphi    = -1;
  Int_t ieta    = -1;
  Int_t iphis   = -1;
  Int_t ietas   = -1;

  Double_t crossEnergy = 0;

  fGeom->GetCellIndex(idmax,iSupMod,iTower,iIphi,iIeta);
  fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphis,ietas);

  Int_t ncells = cluster->GetNCells();
  for (Int_t i=0; i<ncells; i++) {
    Int_t cellAbsId = cluster->GetCellAbsId(i);
    fGeom->GetCellIndex(cellAbsId,iSupMod,iTower,iIphi,iIeta);
    fGeom->GetCellPhiEtaIndexInSModule(iSupMod,iTower,iIphi, iIeta,iphi,ieta);
    Int_t aphidiff = TMath::Abs(iphi-iphis);
    if (aphidiff>1)
      continue;
    Int_t aetadiff = TMath::Abs(ieta-ietas);
    if (aetadiff>1)
      continue;
    if ( (aphidiff==1 && aetadiff==0) ||
	 (aphidiff==0 && aetadiff==1) ) {
      crossEnergy += cells->GetCellAmplitude(cellAbsId);
    }
  }

  return crossEnergy;
}




Double_t AliAnalysisTaskEMCALPi0GammaCorr::GetMaxCellEnergy(const AliVCluster *cluster, Short_t &id) const
{
  // Get maximum energy of attached cell.

  id = -1;

  AliVCaloCells *cells = fCaloCells;
  if(!cells)
    return 0;

  Double_t maxe = 0;
  Int_t ncells = cluster->GetNCells();
  for (Int_t i=0; i<ncells; i++) {
    Double_t e = cells->GetCellAmplitude(TMath::Abs(cluster->GetCellAbsId(i)));
    if (e>maxe) {
      maxe = e;
      id   = cluster->GetCellAbsId(i);
    }
  }
  return maxe;
}


Double_t AliAnalysisTaskEMCALPi0GammaCorr::GetExoticity(AliVCluster *c)
{
  Short_t id = -1;
  Double_t Emax = GetMaxCellEnergy( c, id);
  Double_t Ecross = GetCrossEnergy( c, id);
 
  Double_t exo = 1-Ecross/Emax;
  if(exo>1.0) exo=0.99;
  if(exo<0.0) exo=0.01;
  return exo;
}


Float_t AliAnalysisTaskEMCALPi0GammaCorr::ClustTrackMatching(AliVCluster *clust, double &detaMIN, double &dphiMIN) {
  
  AliTrackContainer* tracks = GetTrackContainer("ForMatching");

  if(!tracks){
    AliError(Form("Could not retrieve tracks !"));
  }

  if(tracks->GetTrackFilterType()!=AliEmcalTrackSelection::kTPCOnlyTracks)  AliError(Form("NO TPC only tracks"));
  AliVTrack* mt = 0;
  TLorentzVector vecClust;
  clust->GetMomentum(vecClust,fVertex);

  Int_t nMatched = clust -> GetNTracksMatched();
  Double_t dR=999.0;
  Double_t deta = 999.0;
  Double_t dphi = 999.0;

  Double_t dR_temp;
  Double_t deta_temp; 
  Double_t dphi_temp;  
    
  if (nMatched <1 ){
      detaMIN = 0.049; 
      dphiMIN = 0.049;
      return 0.0499;
  }

  for(Int_t i=0;i< nMatched;i++){

    Int_t imt = clust->GetTrackMatchedIndex(0);
    if (imt >= 0) mt = static_cast<AliVTrack*>(tracks->GetAcceptParticle(imt));
    if(!mt) continue;

    Double_t veta = mt->GetTrackEtaOnEMCal();
    Double_t vphi = mt->GetTrackPhiOnEMCal();

    Float_t pos[3] = {0};
    clust->GetPosition(pos); //this is the position wrt to nominal vertex 0,0,0. Not considers measured fVertex
    TVector3 cpos(pos);
    Double_t ceta     = cpos.Eta();
    Double_t cphi     = cpos.Phi(); 
    deta_temp =std::abs(veta-ceta);
    dphi_temp =std::abs(TVector2::Phi_mpi_pi(vphi-cphi));
    dR_temp  =TMath::Sqrt(deta_temp*deta_temp+dphi_temp*dphi_temp);
    
    if(dR_temp < dR) dR = dR_temp;
    if(deta_temp < deta) deta = deta_temp;
    if(dphi_temp < dphi ) dphi = dphi_temp;
  }

  //overflow treatment:
  deta = std::min(deta, 0.0499);
  dphi = std::min(dphi, 0.0499);

  detaMIN = deta; 
  dphiMIN = dphi; 
  
  return dR;
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

    //Getting track and cluster containers
    AliTrackContainer* tracks = GetTrackContainer("ForCorrelation");
    AliTrackContainer* tracksForMatching = GetTrackContainer("ForMatching");
    AliTrackContainer* tracksITS = GetTrackContainer("ITSOnly");

    if(!tracks or !tracksForMatching){
      AliError(Form("Could not retrieve tracks !"));
    }
    AliClusterContainer* clusters  = GetClusterContainer(0);
    if(!clusters) AliError(Form("Could not retrieve clusters!")); 

    //Filling cluster and pion THnSparsesx 
    for(auto cluster: clusters->accepted()){
        if(!PreSelection(cluster))continue ;
        FillClusterHisto(cluster, h_Cluster);
        for(auto cluster2: clusters->accepted()){
            if(!PreSelection(cluster2)) continue;
            if(cluster==cluster2) continue;
            FillPionHisto(cluster, cluster2, h_Pi0);
       } 
    }

    //Correlation analysis between photon/pion and tracks
    if(PassedGammaTrigger or fIsMC) { CorrelateClusterAndTrack(tracks,0, kFALSE, 1); }

    //Studying tracks from TPC
    for(auto track : tracksForMatching->accepted()){
        if(!track) continue;
        Float_t d0=-999.0;
        Float_t z0=-999.0;
        track->GetImpactParametersTPC(d0, z0);
        double trigger = 0.0;
        if(PassedGammaTrigger) trigger = 1.0;
        if(d0<-3000) d0 = -1999;
        else if(d0>3000) d0 = 1999;
        double fTPC =static_cast<double>(1.0*track->GetTPCNclsF());
        // std::cout<< track->GetTPCNclsF() << std::endl;
        double charge = static_cast<double>(track->Charge());
        double TPCsignal= std::min(static_cast<double>(track->GetTPCsignal()), 149.0);
        double entries[10] = {track->Pt(), track->Eta() , track->Phi(), 1000*d0, z0, static_cast<double>(track->GetTPCNcls()), TPCsignal, track->GetTPCCrossedRows(),
        			    charge, trigger };
        h_Track->Fill(entries);
    }

    for(auto track :tracksITS->accepted()){
      Float_t d0=-999.0;
      Float_t z0=-999.0;
      track->GetImpactParameters(d0, z0);
      double trigger = 0.0;
      if(PassedGammaTrigger) trigger = 1.0;
      if(d0<-3000) d0 = -1999;
      else if(d0>3000) d0 = 1999;
      double charge = static_cast<double>(track->Charge());

      double entries[7] = {track->Pt(), track->Eta() , track->Phi(), 1000*d0, z0,    charge, trigger };
      h_TrackITS->Fill(entries);


    }





   
    //If MC data, then analyze it: 
    if(fIsMC){
      AnalyzeMC();
    }

    AliEventPool* pool = fPoolMgr->GetEventPool(fCent, zVertex);
    if (!pool)	return kFALSE;
    if(pool->IsReady() && PassedGammaTrigger)
    {
   
        int nMix = pool->GetCurrentNEvents();
        for(int jMix=0; jMix<nMix; jMix++)  {
            CorrelateClusterAndTrack(0, pool->GetEvent(jMix), kTRUE,1.0/nMix);//correlate with mixed event
       }
    }

    
    if(PassedMinBiasTrigger && !PassedGammaTrigger ){
	if(!pool->GetLockFlag())	pool->UpdatePool(CloneToCreateTObjArray(tracks));
    }

  return kTRUE;
}

TObjArray* AliAnalysisTaskEMCALPi0GammaCorr::CloneToCreateTObjArray(AliParticleContainer* tracks)
{
  
   if(!tracks)                            return 0;
   TObjArray* tracksClone = new TObjArray;
   tracksClone->SetOwner(kTRUE);
     
   for(auto track : tracks->accepted()){
       tracksClone->Add(new AliPicoTrack(track->Pt(), track->Eta(), track->Phi(), track->Charge(), 0, 0, 0, 0));
   }
   if(tracksClone->GetEntries()!=tracks->GetNAcceptedParticles())cout<<"!!!!!!! Major error!!!! "<<"Accepted tracks in event: "<< tracks->GetNAcceptedParticles()<<", Tracks in TObjArray: "<<tracksClone->GetEntries()<<endl;
   return tracksClone;
}

int AliAnalysisTaskEMCALPi0GammaCorr::CorrelateClusterAndTrack(AliParticleContainer* tracks,TObjArray* bgTracksArray,Bool_t MixedEvent, double InputWeight)
{
   
    AliClusterContainer* clusters  = GetClusterContainer(0);  
    if (!clusters) return 0;
    double Weight=1.0;    
    Weight=InputWeight; //..for mixed events normalize per events in pool
    
    AliTrackContainer* tracksForMatching = GetTrackContainer("ForMatching");
    if(!tracksForMatching) return 0;
   
    for(auto cluster: clusters->accepted()){
       if(!PreSelection(cluster))continue ;
       if(MixedEvent){
          for(auto track_mix: *bgTracksArray){
	      FillPhotonCorrelation(cluster, static_cast<AliPicoTrack*>(track_mix), h_ClusterTrack_Mixed, Weight);
	  }//end loop over tracks
        }// end if mixed events
        else{
	   for(auto track : tracks->accepted()){
                FillPhotonCorrelation(cluster, track, h_ClusterTrack, Weight);
            } //end loop over tracks
        }//end same event loop.
        
        for(auto cluster2: clusters->accepted()){
   
	  if(!PreSelection(cluster2)) continue;
          if(cluster==cluster2) continue;
	 

            if(MixedEvent){
	      for(auto track_mix: *bgTracksArray){
		FillPionCorrelation(cluster, cluster2, static_cast<AliPicoTrack*>(track_mix), h_Pi0Track_Mixed, Weight);
		} //end loop over tracks
            } // end mixed event loop 
            else{
	  	for(auto track : tracks->accepted()){
                    FillPionCorrelation(cluster, cluster2, track, h_Pi0Track, Weight);
                } //end loop over tracks
            }//end same event 
	 } //end 2 loop over clusters
	} //end  1 loop over clusters

   return 1; 
}


void AliAnalysisTaskEMCALPi0GammaCorr::GetIsolation_Truth(AliVCluster* cluster, double Rmax, double &IsoE){
  
  AliClusterContainer* clusters  = GetClusterContainer(0);
  TLorentzVector reco_photon;
  clusters->GetMomentum(reco_photon, cluster);
 
  double sumET= 0.0;

  AliMCParticleContainer *mcContainer = GetMCParticleContainer("mcparticles");
  if(!mcContainer) AliError(Form("Could not retrieve MCParticleContainer !"));

  Int_t label = TMath::Abs(cluster->GetLabel());
  AliAODMCParticle* true_photon = mcContainer ? mcContainer->GetMCParticleWithLabel(label) : 0x0;   
  if(!true_photon) AliError(Form("Could not retrieve true_photon !"));
  //Loop over final-state particles and sum their 
  for (auto track: mcContainer->accepted()){
    double trackphi = TVector2::Phi_mpi_pi(track->Phi());
    double dphi     = TVector2::Phi_mpi_pi(true_photon->Phi()- trackphi);
    double deta     = true_photon->Eta()- track->Eta();
    double dR       = TMath::Sqrt(deta*deta+dphi*dphi);
    double ET = track->E();//*TMath::Sin(track->Theta());

    if(dR<Rmax){ sumET += ET; } 
  } //end loop over particles

 
  sumET = std::min(sumET, 19.9);
  IsoE       = sumET;
 
  return;
}


void AliAnalysisTaskEMCALPi0GammaCorr::GetIsolation_Track(AliVCluster* cluster, double Rmax, double &IsoE, double &UE_etaband, double &UE_ortho){

  AliClusterContainer* clusters  = GetClusterContainer(0);
  AliTrackContainer* tracks = GetTrackContainer("ForCorrelation");
  if(!tracks){
    AliError(Form("Could not retrieve tracks !"));
  }
  TLorentzVector ph;
  clusters->GetMomentum(ph, cluster);

  double sumpT= 0.0;
  double UE_etaband_temp = 0.0;
  double UE_ortho_temp = 0.0; 

  const double etalimit = 0.9;
  const double minpT    = 0.200; 

  for(auto track : tracks->accepted()){

    if(track->Pt()< minpT) continue; 
    if(std::abs(track->Eta()) > etalimit) continue;

    double trackphi = TVector2::Phi_mpi_pi(track->Phi());
    double dphi     = TVector2::Phi_mpi_pi(ph.Phi()- trackphi);
    double deta     = ph.Eta()- track->Eta();
    double dR       = TMath::Sqrt(deta*deta+dphi*dphi);


    double dphi_ortho1  = TVector2::Phi_mpi_pi(ph.Phi() + TMath::Pi()/2.0 - trackphi);
    double dphi_ortho2    = TVector2::Phi_mpi_pi(ph.Phi() - TMath::Pi()/2.0 - trackphi);


    double dR_ortho1       = TMath::Sqrt(deta*deta+dphi_ortho1*dphi_ortho1);
    double dR_ortho2       = TMath::Sqrt(deta*deta+dphi_ortho2*dphi_ortho2);

    if(dR<Rmax) sumpT = sumpT + track->Pt(); 
    else{
      if(std::abs(dphi)<Rmax)   UE_etaband_temp += track->Pt(); //eta-band 
      if( dR_ortho1<Rmax || dR_ortho2<Rmax  )   UE_ortho_temp += track->Pt(); //orthogonal photon              
    }
    
 
  } //end loop over tracks

  //Subtract estimate of UE from the cone.
  double areaCone = TMath::Pi()*Rmax*Rmax;
  double AUE = 2*etalimit*2*Rmax - areaCone; // = 0.937 for R=0.4 and eta limit of 0.9. 

  UE_ortho   = UE_ortho_temp/2.0; 
  UE_etaband = areaCone*(UE_etaband_temp/AUE);
  IsoE       = sumpT;




  return;
}


void AliAnalysisTaskEMCALPi0GammaCorr::GetIsolation_Cluster(AliVCluster* cluster, double Rmax, double &IsoE, double &UE_etaband, double&IsoE_sub){

  AliClusterContainer* clusters  = GetClusterContainer(0);
  TLorentzVector ph;
  clusters->GetMomentum(ph, cluster);

  double sumpT= 0.0;
  double UE_sumpT = 0.0;

  int NinCone =0;
  int NinUE  = 0;
  TLorentzVector iph;

  for(auto iclus : clusters->accepted()){
    if(iclus==cluster) continue; //not count energy of photon itself
    clusters->GetMomentum(iph, iclus);
   
    //Consider only clusters with energy above 300 MeV and that pass general QA
    if (iclus->E() < 0.300)         continue;
    //    if (std::abs(iph.Eta()) > 0.67) continue; // only consider clusters with eta < 0.67 
    if (!FinalClusterCuts(iclus))   continue;

    double dphi     = TVector2::Phi_mpi_pi(ph.Phi()- iph.Phi());
    double deta = ph.Eta()- iph.Eta();
    double dR= TMath::Sqrt(deta*deta+dphi*dphi);

    if(dR<Rmax){ 
        sumpT = sumpT + iph.Pt();
        NinCone = NinCone +1;
    }
    else if(std::abs(dphi)<Rmax){ 
      UE_sumpT = UE_sumpT + iph.Pt();
      NinUE  = NinUE+1;      
    }
  }

  //Subtract estimate of UE from the cone. 

  double areaCone = TMath::Pi()*Rmax*Rmax;
  double AUE = 2*0.67*2*Rmax - areaCone; // = 0.569 for R=0.4
    
  UE_etaband = areaCone*(UE_sumpT/AUE);
  IsoE       = sumpT;
  IsoE_sub   = IsoE - UE_etaband;
  //Restrict variables to less than 20 GeV, for THnSparse limit. 
  UE_etaband = std::min(UE_etaband, 9.9);
  IsoE = std::min(IsoE, 9.9);
  IsoE_sub = std::min(IsoE_sub, 9.9);
  
  return;
}


void  AliAnalysisTaskEMCALPi0GammaCorr::FillPionCorrelation(AliVCluster* cluster1, AliVCluster* cluster2, AliVParticle* track, THnSparse* histo, double weight){

  
    AliClusterContainer* clusters  = GetClusterContainer(0);
    AliVCluster* cluster_lead = 0;
    AliVCluster* cluster_sub  = 0;
    TLorentzVector ph_lead, ph_sub, pi0; 
    if(cluster1->E()< cluster2->E()) return; //to avoid double-counting   

    cluster_lead = cluster1;
    cluster_sub  = cluster2;

    clusters->GetMomentum(ph_lead, cluster_lead);
    clusters->GetMomentum(ph_sub, cluster_sub);
    double asym = std::abs(ph_lead.Pt()-ph_sub.Pt())/(ph_lead.Pt()+ph_sub.Pt());
    pi0= ph_lead+ph_sub;

    double openingAngle = 1000.0*std::abs(TVector2::Phi_mpi_pi(ph_lead.Phi()-ph_sub.Phi())); // in mrads
    //////////////////Selection/////////////////////////////////////////////
    if( !FinalClusterCuts(cluster_lead)) return; 
    if( !FinalClusterCuts(cluster_sub)) return;
    
    if(cluster_lead->E() < 3.0) return;
    if(cluster_sub->E()  < 3.0) return;

    if( pi0.Pt() < 8.0    ) return;
    if( pi0.M()  > 0.3    ) return;
    if( asym > 0.7        ) return; 
    if( openingAngle < 17 ) return;
    if( track->Pt()< 1.0  ) return; 
    /////////////////////////////////////////////////////////////////////////
    
    double deta = pi0.Eta()-track->Eta();
    double  Zt  = track->Pt()/pi0.Pt();
    double  Xi  = -999; 
    if(Zt>0) Xi = TMath::Log(1.0/Zt);
    double trackphi = TVector2::Phi_mpi_pi(track->Phi());
    double dphi;
    dphi     = TVector2::Phi_mpi_pi(pi0.Phi()- trackphi)/TMath::Pi();
    if(dphi<-0.5) dphi +=2;
    
    double trackpT = std::min(track->Pt(),10.0);

    Double_t zVertex = fVertex[2];
    if (zVertex>10) zVertex =9.99;
    if (zVertex<-10) zVertex = -9.99;


    double entries[13] = {fCent, zVertex, pi0.Pt(),  
			 trackpT, dphi, deta, Zt, Xi,
                         pi0.M(), ph_lead.Pt(), ph_sub.Pt(),  cluster_lead->GetM02(), cluster_sub->GetM02()};   
             
    histo->Fill(entries, weight); //the factor 0.5 is to account for the fact that this  if filled twice.  
    return;
}

void  AliAnalysisTaskEMCALPi0GammaCorr::FillPhotonCorrelation(AliVCluster* cluster, AliVParticle* track, THnSparse* histo, double weight){
    AliClusterContainer* clusters  = GetClusterContainer(0);
    TLorentzVector ph;
    clusters->GetMomentum(ph, cluster);
   
    if( track->Pt()<1.0) return;
    if( ph.Pt() < 8.0) return;
    if(!FinalClusterCuts(cluster)) return;
    Double_t detamin = 0.0; 
    Double_t dphimin = 0.0; 

    Double_t dRmin = ClustTrackMatching(cluster, detamin, dphimin);
    
    double trackphi = TVector2::Phi_mpi_pi(track->Phi());
    double trackpt = track->Pt();
    double dphi;
    
    double deta = ph.Eta()-track->Eta();
    double  Zt  = track->Pt()/ph.Pt();
    double  Xi  = -999; 
    if(Zt>0) Xi = TMath::Log(1.0/Zt);
    Zt = std::min(Zt, 1.19);    
    trackpt = std::min(trackpt, 9.99);

    dphi = TVector2::Phi_mpi_pi(ph.Phi()- trackphi)/TMath::Pi();
    if(dphi<-0.5) dphi +=2;
    
    double entries[10] = {fCent, fVertex[2], ph.Pt(),  
                          trackpt, dphi, deta, Zt, Xi,
			  cluster->GetM02(), dRmin };                
    histo->Fill(entries, weight);//
    return;
}





void  AliAnalysisTaskEMCALPi0GammaCorr::FillPionHisto(AliVCluster* cluster1, AliVCluster* cluster2, THnSparse* histo){
    
  //std::cout << " Entering Fill Pion " << std::endl;
    AliClusterContainer* clusters  = GetClusterContainer(0);

    if(cluster1->E() < 0.7) return; 
    if(cluster2->E() < 0.7) return; 
    if(!FinalClusterCuts(cluster1)) return;
    if(!FinalClusterCuts(cluster2)) return;

    AliVCluster* cluster_lead = 0;
    AliVCluster* cluster_sub  = 0;
    
    TLorentzVector ph_lead, ph_sub, pi0; 
    if(cluster1->E() < cluster2->E()) return; //to avoid double-counting
   
    cluster_lead = cluster1;
    cluster_sub  = cluster2;
    
    clusters->GetMomentum(ph_lead, cluster_lead);
    clusters->GetMomentum(ph_sub,  cluster_sub);

    Double_t detamin_1 = 0.0;
    Double_t dphimin_1 = 0.0;
    Double_t detamin_2 = 0.0;
    Double_t dphimin_2 = 0.0;
 
    Double_t dRmin_1 = ClustTrackMatching(cluster_lead, detamin_1, dphimin_1);
    Double_t dRmin_2 = ClustTrackMatching(cluster_sub, detamin_2, dphimin_2);    

    pi0 = ph_lead + ph_sub;
   
    //////////////////Selection/////////////////////////////////////////
    if( pi0.Pt() < 6.0) return;
    if( pi0.M()  > 0.3) return;
    ////////////////////////////////////////////////////////////////////
    double asym = std::abs(ph_lead.Pt()-ph_sub.Pt())/(ph_lead.Pt()+ph_sub.Pt());
    double openingAngle = 1000.0*std::abs(TVector2::Phi_mpi_pi(ph_lead.Phi()-ph_sub.Phi())); // in mrads
    Double_t zVertex = fVertex[2];
    if (zVertex>10) zVertex =9.99;
    if (zVertex<-10) zVertex = -9.99;  

    double trigger = 0.0;
    if(PassedGATrigger()) trigger = 1.0;
 
    //Check whether pion is TRUE pion or not: 
    double true_pt = 0.0; 
    if(fIsMC) IsRealPion(cluster_lead, cluster_sub, true_pt);
    
    double entries[15] = {fCent, zVertex, pi0.M(), pi0.Pt(), pi0.Rapidity(),  asym, ph_lead.Pt(), ph_sub.Pt(),  
			  openingAngle,  cluster_lead->GetM02(), cluster_sub->GetM02(), 
			  dRmin_1, dRmin_2, true_pt, trigger};

    histo->Fill(entries);
    return;
}


Bool_t AliAnalysisTaskEMCALPi0GammaCorr::IsRealPion(AliVCluster* cluster_1, AliVCluster* cluster_2, double &truepT){

  
 
  AliMCParticleContainer *mcContainer = GetMCParticleContainer("mcparticles");
  if(!mcContainer) AliError(Form("Could not retrieve MCParticleContainer !"));
  Int_t label_1 = TMath::Abs(cluster_1->GetLabel());
  Int_t label_2 = TMath::Abs(cluster_2->GetLabel());
 
  AliAODMCParticle* true_photon_1 = mcContainer ? mcContainer->GetMCParticleWithLabel(label_1) : 0x0;
  AliAODMCParticle* true_photon_2 = mcContainer ? mcContainer->GetMCParticleWithLabel(label_2) : 0x0;
  
  if(!true_photon_2 or !true_photon_2) AliError(Form("Could not retrieve true_photon !"));
  
 
  if(!true_photon_1 or !true_photon_2) return kFALSE;
  if(true_photon_1->PdgCode()!=22) return kFALSE;
  if(true_photon_2->PdgCode()!=22 ) return kFALSE; 
  if(true_photon_1->GetMother()<0) return kFALSE;
  if(true_photon_2->GetMother()<0) return kFALSE;
  if(true_photon_1->GetMother()!=true_photon_2->GetMother()) return kFALSE;

  Int_t motherlabel = TMath::Abs(true_photon_1->GetMother());
 
  AliAODMCParticle* true_mother = mcContainer ? mcContainer->GetMCParticleWithLabel(motherlabel) : 0x0;
  if(true_mother->PdgCode()!=111) return kFALSE; 
 
  truepT = true_mother->Pt();
  return kTRUE;
}

void AliAnalysisTaskEMCALPi0GammaCorr::FillClusterHisto(AliVCluster* cluster, THnSparse* histo){
    
    AliClusterContainer* clusters  = GetClusterContainer(0);
    AliTrackContainer* tracks = GetTrackContainer("ForCorrelation");

    double trigger = 0.0;
    if(PassedGATrigger()) trigger = 1.0;

    TLorentzVector ph;
    clusters->GetMomentum(ph, cluster);
    if(ph.Pt() < 5.0) return;

    Double_t detamin = 0.0; 
    Double_t dphimin = 0.0; 

    Double_t dRmin = ClustTrackMatching(cluster, detamin, dphimin);
    Double_t disToBad = std::min( static_cast<double>(cluster->GetDistanceToBadChannel()), 5.0);
   
    Double_t disToBorder = std::min(static_cast<double>(GetMaxDistanceFromBorder(cluster)), 5.0);
    Double_t exoticity = GetExoticity(cluster);
    Double_t time = cluster->GetTOF()*1000000000; //in ns
    if (time<-40) time = -40;
    if (time>40) time = +40;

    Double_t RunNumber = static_cast<double>(FormatRunNumber(fInputEvent->GetRunNumber()));
    //Double_t BCID      = static_cast<double>(fInputEvent->GetBunchCrossNumber());


    double defValue = -4.9;    
  
    Double_t UE_Tracks_etaband = defValue;
    Double_t UE_Tracks_ortho   = defValue; 
    Double_t IsoE_Tracks       = defValue;

    Double_t IsoE_Truth        = defValue;
 

    GetIsolation_Track(cluster, 0.4, IsoE_Tracks, UE_Tracks_etaband, UE_Tracks_ortho);
    
   
    if(fIsMC) GetIsolation_Truth(cluster, 0.4, IsoE_Truth); 
  
    Double_t nTracks = std::min( static_cast<double>(tracks->GetNAcceptedTracks()), 99.9);
    Double_t nClusters = std::min(static_cast<double>(clusters->GetNAcceptedClusters()), 99.9);


    Double_t zVertex = fVertex[2];
    if (zVertex>10) zVertex =9.9;
    if (zVertex<-10) zVertex = -9.99;

    //Restrict variables to less than 20 GeV, for THnSparse limit.
    const double ulimit = 9.9;
    
    double trueGamma = 0.0; 
    
    

    if(fIsMC){
        AliMCParticleContainer *mcContainer = GetMCParticleContainer("mcparticles");
        if(!mcContainer) AliError(Form("Could not retrieve MCParticleContainer !"));
        Int_t label = TMath::Abs(cluster->GetLabel());
        AliAODMCParticle* true_photon = mcContainer ? mcContainer->GetMCParticleWithLabel(label) : 0x0;
        if(!true_photon) AliError(Form("Could not retrieve true_photon !"));
        if(true_photon->PdgCode()==22) trueGamma=.60;
    }

    double entries[26] = {RunNumber, fCent, zVertex, ph.Pt(), ph.Eta(), ph.Phi(), cluster->GetM02(), static_cast<double>(cluster->GetNCells()), 
			  static_cast<double>(cluster->GetNExMax()), disToBorder, disToBad, dRmin, detamin, dphimin, exoticity, time, nTracks, nClusters,
                          std::min(IsoE_Tracks, ulimit),
                          std::min(UE_Tracks_etaband, ulimit), 
			  std::min(UE_Tracks_ortho, ulimit),
			  std::min(IsoE_Tracks- UE_Tracks_etaband, ulimit),
			  std::min(IsoE_Tracks- UE_Tracks_ortho, ulimit),
                          std::min(IsoE_Truth, 19.9 ),
                          trueGamma, trigger
                          };

    histo->Fill(entries);
    return;
}

TObjArray* AliAnalysisTaskEMCALPi0GammaCorr::CloneClustersTObjArray(AliClusterContainer* clusters)
{
	if(!clusters)  return 0;
	if(clusters->GetNClusters()==0) return 0;
	TObjArray* clustersCloneI = new TObjArray;
	clustersCloneI->SetOwner(kTRUE);
        for(auto cluster: clusters->accepted()){
            clustersCloneI->Add((AliVCluster*)cluster);
	}
	if(clustersCloneI->GetEntries()!=clusters->GetNAcceptedClusters())cout<<"!!!!!!! Major error!!!! "<<"Accepted clusters in event: "<<clusters->GetNAcceptedClusters()<<", Tracks in TObjArray: "<<clustersCloneI->GetEntries()<<endl;
	return clustersCloneI;
}


Int_t  AliAnalysisTaskEMCALPi0GammaCorr::GetMaxDistanceFromBorder(AliVCluster* cluster){

  Int_t max = 0;

  for (int n=0; n <6; n++){
    fFiducialCellCut->SetNumberOfCellsFromEMCALBorder(n);
    if(fFiducialCellCut->CheckCellFiducialRegion(fGeom, cluster,fCaloCells))
      {
        max = n;
      }
    else{break;}
  } 

  return max;
}


Bool_t AliAnalysisTaskEMCALPi0GammaCorr::PreSelection(AliVCluster* cluster)
{
    if(!cluster->IsEMCAL()) return kFALSE;
    if(cluster->E()<0.7) return kFALSE;
    return kTRUE;
}

Bool_t AliAnalysisTaskEMCALPi0GammaCorr::FinalClusterCuts(AliVCluster* cluster)
{

  //General QA. 
  if(!cluster->IsEMCAL()) return kFALSE;

  if( cluster->GetNCells() < 2) return kFALSE;

  Int_t disToBad = cluster->GetDistanceToBadChannel();
  if(disToBad<2) return kFALSE;
  
  Int_t disToBorder = GetMaxDistanceFromBorder(cluster);
  if(disToBorder<1) return kFALSE;

  Double_t exoticity = GetExoticity(cluster);
  if(exoticity>0.97) return kFALSE;

  Double_t time = cluster->GetTOF()*1000000000; //in ns
  if(!fIsMC && std::abs(time)>30) return kFALSE;

  return kTRUE;
}



double AliAnalysisTaskEMCALPi0GammaCorr::GetEff(AliTLorentzVector ClusterVec)
{
	return 1;
}

Int_t AliAnalysisTaskEMCALPi0GammaCorr::FormatRunNumber(Int_t runnumber)
{
  //This is the list for LHC13d pPb 5 TeV pass4. 
  switch (runnumber) {
  case  195872 : return 10;
  case  195871 : return 9;
  case  195867 : return 8;
  case  195831 : return 7;
  case  195829 : return 6;
  case  195787 : return 5;
  case  195783 : return 4;
  case  195767 : return 3;
  case  195760 : return 2;
  case  195724 : return 1;
  default : return 0;
  }
}

void AliAnalysisTaskEMCALPi0GammaCorr::AnalyzeMC(){
 

  double truepT= 0.0; 
  double truey = 0.0; 
  double truePDG = 0.0;

  AliMCParticleContainer *mcContainer = GetMCParticleContainer("mcparticles");
  for (auto track: mcContainer->all())
  {
      if(track->PdgCode()!=111 and track->PdgCode()!=22) continue; //only keep photons and pions
      truepT = track->Pt();
      if(truepT<5.0) continue; //only keep large-momentum pions and photons
      truey  = track->Eta();
      if(track->PdgCode()==111) truePDG = 0.5;
      else if(track->PdgCode()==22) truePDG = 1.5;
      else truePDG = 2.5;
      double entries[3] = {truepT, truey, truePDG};
      h_Truth->Fill(entries);
  }
  

  return;
}
