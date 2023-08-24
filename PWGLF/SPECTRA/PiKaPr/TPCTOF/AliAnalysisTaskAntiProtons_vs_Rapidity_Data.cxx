#include "AliAnalysisTaskAntiProtons_vs_Rapidity_Data.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliAnalysisUtils.h"
#include "AliAnalysisTask.h"
#include "AliPIDResponse.h"
#include "TLorentzVector.h"
#include "AliEventCuts.h"
#include "TDatabasePDG.h"
#include "AliESDVertex.h"
#include "AliESDtrack.h"
#include "AliESDEvent.h"
#include "THnSparse.h"
#include "AliVVZERO.h"
#include "TObjArray.h"
#include "TChain.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TF1.h"
#include <vector>
#include <iostream>
#include <algorithm>

using namespace std;
ClassImp(AliAnalysisTaskAntiProtons_vs_Rapidity_Data)

//__________________________________________________________________________________________________________________________________
AliAnalysisTaskAntiProtons_vs_Rapidity_Data::AliAnalysisTaskAntiProtons_vs_Rapidity_Data():
AliAnalysisTaskSE(),
fESDEvent(nullptr),
fPIDResponse(nullptr),
fESDtrackCuts{},
fOutputList(nullptr),
fQAList(nullptr),
fIsITSrecalib(kFALSE),
fESDEventSelection(),
hMean(nullptr),
hWidth(nullptr),
hNumberOfEvents(nullptr),
hMultiplicity(nullptr),
hTPCnsigma(nullptr),
hTOFnsigma(nullptr),
hTPCnsigma_vs_rap(nullptr),
hTOFnsigma_vs_rap(nullptr),
hDCAxy(nullptr),
hITSnsigma(nullptr)
{}
//__________________________________________________________________________________________________________________________________
AliAnalysisTaskAntiProtons_vs_Rapidity_Data::AliAnalysisTaskAntiProtons_vs_Rapidity_Data(const char *name):
AliAnalysisTaskSE(name),
fESDEvent(nullptr),
fPIDResponse(nullptr),
fESDtrackCuts{},
fOutputList(nullptr),
fQAList(nullptr),
fIsITSrecalib(kFALSE),
fESDEventSelection(),
hMean(nullptr),
hWidth(nullptr),
hNumberOfEvents(nullptr),
hMultiplicity(nullptr),
hTPCnsigma(nullptr),
hTOFnsigma(nullptr),
hTPCnsigma_vs_rap(nullptr),
hTOFnsigma_vs_rap(nullptr),
hDCAxy(nullptr),
hITSnsigma(nullptr)
{
    DefineInput (0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
}
//__________________________________________________________________________________________________________________________________
AliAnalysisTaskAntiProtons_vs_Rapidity_Data::~AliAnalysisTaskAntiProtons_vs_Rapidity_Data()  {
    
    fOutputList->Clear();
    delete fESDEvent;
    delete fPIDResponse;
    for (Int_t isyst=0 ; isyst<50 ; isyst++) { delete fESDtrackCuts[isyst]; }

    delete fOutputList;
    delete fQAList;
    delete hNumberOfEvents;
    delete hMultiplicity;
    delete hTPCnsigma;
    delete hTOFnsigma;
    delete hTPCnsigma_vs_rap;
    delete hTOFnsigma_vs_rap;
    delete hDCAxy;
    delete hITSnsigma;
}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskAntiProtons_vs_Rapidity_Data::UserCreateOutputObjects()  {
    
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
    
    
    //Event Counter and Centrality Distribution
    hNumberOfEvents = new TH1F ("hNumberOfEvents","",20,0,20);
    hMultiplicity   = new TH1F ("hMultiplicity","",100,0,100);
    hNumberOfEvents -> Sumw2();
    hMultiplicity   -> Sumw2();
    fOutputList -> Add(hNumberOfEvents);
    fOutputList -> Add(hMultiplicity);

    //Binning DCA_{xy}
    Int_t    bins_dcaxy[4]  = { 50,   7,  40,  100 };
    Double_t xmin_dcaxy[4]  = {  0, 0.0, 0.0, -1.0 };
    Double_t xmax_dcaxy[4]  = { 50, 0.7, 2.0,  1.0 };

    //Binning nsigma_{TPC}
    Int_t    bins_nsigmaTPC[4]  = { 50,   7,  20,   100 };
    Double_t xmin_nsigmaTPC[4]  = {  0, 0.0, 0.0, -10.0 };
    Double_t xmax_nsigmaTPC[4]  = { 50, 0.7, 1.0,  10.0 };

    //Binning nsigma_{TOF}
    Int_t    bins_nsigmaTOF[4]  = { 50,   7,  30,   200 };
    Double_t xmin_nsigmaTOF[4]  = {  0, 0.0, 0.5, -20.0 };
    Double_t xmax_nsigmaTOF[4]  = { 50, 0.7, 2.0,  20.0 };

    //DCA_{xy}
    hDCAxy = new THnSparseF ("hDCAxy","",4, bins_dcaxy, xmin_dcaxy, xmax_dcaxy);
    hDCAxy -> Sumw2();
    fOutputList -> Add (hDCAxy);
    
    //nsigma_{TPC}
    hTPCnsigma = new THnSparseF ("hTPCnsigma","",4, bins_nsigmaTPC, xmin_nsigmaTPC, xmax_nsigmaTPC);
    hTPCnsigma -> Sumw2();
    fOutputList -> Add (hTPCnsigma);
    
    //nsigma_{TOF}
    hTOFnsigma = new THnSparseF ("hTOFnsigma","",4, bins_nsigmaTOF, xmin_nsigmaTOF, xmax_nsigmaTOF);
    hTOFnsigma -> Sumw2();
    fOutputList -> Add (hTOFnsigma);
    
    
    //**************************************** Comparison y>0 vs. y<0 ****************************************
    
    //Binning nsigma_{TPC}
    Int_t    bins_nsigmaTPC_vs_y[3]  = {   14,  20,   100 };
    Double_t xmin_nsigmaTPC_vs_y[3]  = { -0.7, 0.0, -10.0 };
    Double_t xmax_nsigmaTPC_vs_y[3]  = {  0.7, 1.0,  10.0 };

    //Binning nsigma_{TOF}
    Int_t    bins_nsigmaTOF_vs_y[3]  = {   14,  30,   200 };
    Double_t xmin_nsigmaTOF_vs_y[3]  = { -0.7, 0.5, -20.0 };
    Double_t xmax_nsigmaTOF_vs_y[3]  = {  0.7, 2.0,  20.0 };
    
    //nsigma_{TPC}
    hTPCnsigma_vs_rap = new THnSparseF ("hTPCnsigma_vs_rap","",3, bins_nsigmaTPC_vs_y, xmin_nsigmaTPC_vs_y, xmax_nsigmaTPC_vs_y);
    hTPCnsigma_vs_rap -> Sumw2();
    fOutputList -> Add (hTPCnsigma_vs_rap);
    
    //nsigma_{TOF}
    hTOFnsigma_vs_rap = new THnSparseF ("hTOFnsigma_vs_rap","",3, bins_nsigmaTOF_vs_y,xmin_nsigmaTOF_vs_y, xmax_nsigmaTOF_vs_y);
    hTOFnsigma_vs_rap -> Sumw2();
    fOutputList -> Add (hTOFnsigma_vs_rap);
    
    
    //2D ITS Re-calibration Map
    if (fIsITSrecalib)  {
        hITSnsigma = new TH3F ("hITSnsigma","",20,-1.0,1.0,20,0.0,1.0,200,-5,5);//eta,p,nsigma_{ITS}
        hITSnsigma  -> Sumw2();
        fOutputList -> Add (hITSnsigma);
    }
    
    //Track Selection
    for (Int_t isyst=0 ; isyst<50 ; isyst++)  {
        fESDtrackCuts[isyst] = new AliESDtrackCuts(Form("fESDtrackCuts[%d]",isyst));
        fESDtrackCuts[isyst] -> SetAcceptKinkDaughters(kFALSE);
        fESDtrackCuts[isyst] -> SetRequireTPCRefit(kTRUE);
        fESDtrackCuts[isyst] -> SetRequireITSRefit(kTRUE);
        fESDtrackCuts[isyst] -> SetMinNCrossedRowsTPC(50);
        fESDtrackCuts[isyst] -> SetMinRatioCrossedRowsOverFindableClustersTPC(0.5);
        fESDtrackCuts[isyst] -> SetMinNClustersITS(1);
        fESDtrackCuts[isyst] -> SetMaxChi2PerClusterITS(36);
        fESDtrackCuts[isyst] -> SetMaxChi2PerClusterTPC(10);
    }
    
    
    //Post Data
    PostData(1, fOutputList);
    PostData(2, fQAList);
}
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskAntiProtons_vs_Rapidity_Data::UserExec(Option_t *)  {
    
    //Get Input Event
    if (!GetEvent()) return;

    //TPC Pre-selection
    Double_t nsigmaTPCmax[50] = {3.0,2.62945,2.69333,3.1122,3.28161,2.70606,2.89633,2.70347,3.31368,2.65839,2.93157,3.10837,3.43406,2.51501,3.49459,3.27995,3.05988,2.89139,2.60827,3.16441,2.69537,2.71876,2.67132,2.80894,3.42022,3.23984,2.96282,2.84138,2.78734,3.102,3.15273,2.95713,2.90697,3.05029,2.86195,3.16312,3.43193,2.57391,3.38719,3.177,3.38128,3.40739,2.83285,2.7041,3.05584,2.58265,2.99982,3.40591,2.64627,3.45544};
    
    //ITS Pre-selection
    Double_t nsigmaITSmax[50] = {3.0,2.62945,2.69333,3.1122,3.28161,2.70606,2.89633,2.70347,3.31368,2.65839,2.93157,3.10837,3.43406,2.51501,3.49459,3.27995,3.05988,2.89139,2.60827,3.16441,2.69537,2.71876,2.67132,2.80894,3.42022,3.23984,2.96282,2.84138,2.78734,3.102,3.15273,2.95713,2.90697,3.05029,2.86195,3.16312,3.43193,2.57391,3.38719,3.177,3.38128,3.40739,2.83285,2.7041,3.05584,2.58265,2.99982,3.40591,2.64627,3.45544};
    
    //Proton Mass
    const Double_t mass = AliPID::ParticleMass(AliPID::kProton);

    //Loop over Reconstructed Tracks
    for (Int_t i=0 ; i<fESDEvent->GetNumberOfTracks() ; i++)  {

        //Get Reconstructed Track
        AliESDtrack *track = (AliESDtrack*) fESDEvent->GetTrack(i);
        if (!track) continue;
        if (track->Pt()>2.0)   continue;
        if (track->Charge()>0) continue;

        //Variables
        Double_t nsigmaITS = fPIDResponse -> NumberOfSigmasITS (track,AliPID::kProton);
        Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kProton);
        Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kProton);
        Double_t eta       = track->Eta();
        Double_t pt        = track->Pt();
        Double_t p         = track->P();
        Double_t y         = GetRapidity(track,mass);
        Double_t DCAxy     = GetDCAtoPrimaryVertex (track,0);
        Bool_t   hasTOFhit = (track->GetStatus() & AliVTrack::kTOFout) && (track->GetStatus() & AliVTrack::kTIME);
        Double_t length    = track->GetIntegratedLength();
        
        //ITS Recalibration
        Double_t nsigmaITS_recalib = GetRecalibratedITSnsigma (nsigmaITS,eta,p);

        //Fill ITS Re-calibration Map
        if (fIsITSrecalib)  {
            if (!IsProtonCandidate(track)) continue;
            hITSnsigma -> Fill (eta,p,nsigmaITS_recalib);
            continue;
        }
                
        for (Int_t isyst=0 ; isyst<50 ; isyst++)  {
            
            //Track Quality Cuts
            if (!PassedTrackSelection (track,isyst)) continue;

            //Variables
            Double_t var_dca[4] = {static_cast<double>(isyst)+0.5,TMath::Abs(y),pt,DCAxy};
            Double_t var_tpc[4] = {static_cast<double>(isyst)+0.5,TMath::Abs(y),pt,nsigmaTPC};
            Double_t var_tof[4] = {static_cast<double>(isyst)+0.5,TMath::Abs(y),pt,nsigmaTOF};
            Double_t tpc_rap[3] = {y,pt,nsigmaTPC};
            Double_t tof_rap[3] = {y,pt,nsigmaTOF};
            
            //DCA_{xy} Distributions
            if (IsHighPurityProton(track)) hDCAxy -> Fill (var_dca);
            
            //DCA_{xy} Selection
            if (TMath::Abs(DCAxy)>0.1) continue;

            //TPC Analysis
            if (pt<1.0 && TMath::Abs(nsigmaITS_recalib)<nsigmaITSmax[isyst]) hTPCnsigma -> Fill (var_tpc);
            if (pt<1.0 && TMath::Abs(nsigmaITS_recalib)<nsigmaITSmax[0] && isyst==0) hTPCnsigma_vs_rap -> Fill (tpc_rap);

            //TOF Analysis
            if (pt<0.5)       continue;
            if (!hasTOFhit)   continue;
            if (length<350.0) continue;
            if (TMath::Abs(nsigmaTPC)>nsigmaTPCmax[isyst]) continue;
            
            //Fill TOF Spectra
            hTOFnsigma -> Fill (var_tof);
            if (isyst==0) hTOFnsigma_vs_rap -> Fill (tof_rap);
        }
    }
    
    
    
    //Check
    if (fESDEvent->GetNumberOfTracks()==0) hNumberOfEvents->Fill(18.5);
    
    //Post Output Data
    PostData(1, fOutputList);
}
//__________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskAntiProtons_vs_Rapidity_Data::GetEvent ()  {
    
    //Get Input Event
    fESDEvent = dynamic_cast <AliESDEvent*>(InputEvent());
    if (!fESDEvent) return kFALSE;
    hNumberOfEvents -> Fill(0.5);
    
    //Standard Event Cuts
    if (!fESDEventSelection.AcceptEvent(fESDEvent)) {
        PostData(2, fQAList);
        return kFALSE;
    }
    hNumberOfEvents -> Fill(1.5);

    //Reject Events with Incomplete DAQ
    if (fESDEvent->IsIncompleteDAQ()) return kFALSE;
    hNumberOfEvents -> Fill(2.5);
       
    //V0 Timing Decision
    AliVVZERO *vzeroData = fESDEvent->GetVZEROData();
    if (!(vzeroData->GetV0ADecision()) || !(vzeroData->GetV0CDecision())) return kFALSE;
    hNumberOfEvents -> Fill(3.5);
       
    //Pileup Rejection
    Int_t nClustersLayer0 = fESDEvent->GetNumberOfITSClusters(0);
    Int_t nClustersLayer1 = fESDEvent->GetNumberOfITSClusters(1);
    Int_t nTracklets      = fESDEvent->GetMultiplicity()->GetNumberOfTracklets();
    if ((nClustersLayer0 + nClustersLayer1) > 65.0 + (Double_t)nTracklets*4.0) return kFALSE;
    hNumberOfEvents -> Fill(4.5);

    //Primary Vertex Tracks
    AliESDVertex *vertex_tracks = (AliESDVertex*) fESDEvent->GetPrimaryVertexTracks();
    if (!vertex_tracks) return kFALSE;
    hNumberOfEvents -> Fill(5.5);
       
    //Vertex Contributors Tracks
    if ( vertex_tracks->GetNContributors() < 1 ) return kFALSE;
    hNumberOfEvents -> Fill(6.5);
       
    //Primary Vertex SPD
    AliESDVertex *vertex_SPD = (AliESDVertex*) fESDEvent->GetPrimaryVertexSPD();
    if (!vertex_SPD) return kFALSE;
    hNumberOfEvents -> Fill(7.5);
       
    //Vertex Contributors SPD
    if ( vertex_SPD->GetNContributors() < 1 ) return kFALSE;
    hNumberOfEvents -> Fill(8.5);
       
    //SPD Pile-up in Mult Bins
    if (fESDEvent->IsPileupFromSPDInMultBins()) return kFALSE;
    hNumberOfEvents -> Fill(9.5);
       
    //Cut on Z-Vertex Resolution
    if (TMath::Abs(vertex_SPD->GetZ() - vertex_tracks->GetZ()) > 0.3) return kFALSE;
    hNumberOfEvents -> Fill(10.5);

    //Primary Vertex Selection
    if ( vertex_tracks->GetZ() < -10.0 ) return kFALSE;
    if ( vertex_tracks->GetZ() > +10.0 ) return kFALSE;
    hNumberOfEvents -> Fill(11.5);
              
    //Multiplicity
    AliMultSelection *multiplicitySelection = (AliMultSelection*) fESDEvent->FindListObject("MultSelection");
    if( !multiplicitySelection) return kFALSE;
    hNumberOfEvents -> Fill(12.5);
       
    //Multiplicity Distribution
    Double_t mult_percentile = multiplicitySelection->GetMultiplicityPercentile("V0M");
    hMultiplicity -> Fill(mult_percentile);
       
    //Selection of Multiplicity Range
    if (mult_percentile <   0) return kFALSE;
    if (mult_percentile > 100) return kFALSE;
    hNumberOfEvents -> Fill(13.5);
    
    //Load PID Response
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();

    return kTRUE;
}
//__________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskAntiProtons_vs_Rapidity_Data::PassedTrackSelection (AliESDtrack *track, Int_t isyst)  {
    
    //Initialization
    Bool_t passedTrkSelection=(kFALSE);
    
    //Basic Track Selection
    if ( !fESDtrackCuts[isyst]->AcceptTrack (track) ) return passedTrkSelection;

    //Fixed Cuts
    if (TMath::Abs(GetDCAtoPrimaryVertex(track,0))>1.0) return passedTrkSelection;
    if (!track->HasPointOnITSLayer(0))                  return passedTrkSelection;
    if (!track->HasPointOnITSLayer(1))                  return passedTrkSelection;

    //Track Variables
    Int_t nTPCcr      = track->GetTPCCrossedRows();
    Int_t nTPCfind    = track->GetTPCNclsF();
    Int_t nITScls     = track->GetITSNcls();
    Int_t nTPCcls     = track->GetTPCNcls();
    Int_t nTPCclsdEdx = track->GetTPCsignalN();
    Double_t chi2ITS  = track->GetITSchi2();
    Double_t chi2TPC  = track->GetTPCchi2();
    Double_t DCAxy    = TMath::Abs(GetDCAtoPrimaryVertex(track,0));
    Double_t DCAz     = TMath::Abs(GetDCAtoPrimaryVertex(track,1));
    Double_t cr_over_findable = (Double_t)nTPCcr/((Double_t)nTPCfind);
    Double_t chi2TPC_NDF      = chi2TPC/((Double_t)nTPCcls) ;
    
    //Set of Cuts
    Int_t nTPCcr_min[50] = {80,97,71,75,86,87,74,71,71,91,82,77,75,85,80,88,94,83,85,82,74,86,75,93,96,96,84,83,71,94,84,82,72,85,76,86,83,82,93,90,92,92,70,83,91,93,97,82,80,95};

    Int_t nITScls_min[50] = {5,4,3,3,6,3,4,4,5,4,4,6,3,4,4,3,4,4,3,4,4,5,4,3,5,4,6,4,5,4,3,4,5,4,4,6,4,4,4,5,4,4,5,4,3,6,4,3,4,4};

    Int_t nTPCclsdEdx_min[50] = {60,63,54,57,57,64,69,60,54,62,66,50,55,64,68,50,67,59,52,64,67,64,50,55,56,52,68,70,54,57,62,63,50,50,68,64,63,67,53,58,54,54,68,56,51,56,63,59,56,58};

    Double_t chi2ITS_max[50] = {36,46.1797,36.6257,48.2043,36.8251,46.3546,44.4219,36.9071,36.4452,45.4235,44.8695,47.7977,47.0017,44.1026,37.5981,45.6092,49.7681,44.8784,38.9593,39.0638,48.7178,49.4983,49.3736,45.7836,41.407,44.4918,41.3072,36.7743,47.0045,43.5901,45.2098,43.167,46.9972,39.2736,38.6093,42.4406,47.2409,43.3473,45.1781,41.2574,45.8378,36.263,47.7494,48.1092,39.5733,42.8802,46.6033,47.2696,47.3724,45.3757};

    Double_t chi2TPC_NDF_max[50] = {4.0,3.99601,2.79669,3.88701,3.73929,3.07139,4.14759,3.19747,4.58796,4.53536,3.56903,4.85172,3.70836,5.01696,5.05308,4.16061,3.40112,4.67905,4.40227,5.0066,4.95311,3.13955,3.35452,4.43888,4.21231,3.8571,3.43539,5.16655,5.19831,4.4465,4.1092,4.57635,2.61833,4.17145,4.14075,3.91889,5.05004,3.12027,3.86329,3.47401,5.27806,3.50963,4.18327,4.7357,3.71845,2.67676,3.2959,4.27114,2.86657,3.1244};

    Double_t dcaz_max[50] = {1.0,0.524845,0.838556,1.36511,0.972847,0.644151,1.14321,0.557206,0.978874,1.14591,0.958042,0.950637,0.718069,0.520478,1.26181,0.958435,0.532273,1.13665,1.1004,1.37577,1.02552,0.571379,0.708928,0.831586,1.44803,0.90449,1.0531,0.627843,1.42758,1.4044,0.56425,0.901189,0.987433,0.816606,1.37016,1.32486,1.30484,0.572152,1.34936,0.535521,0.819561,0.575574,1.44179,1.25525,1.36479,0.522065,0.535635,1.26045,0.851564,0.599568};

    Double_t cr_over_findable_min[50] = {0.8,0.849002,0.781492,0.897933,0.704228,0.706126,0.704632,0.709122,0.897615,0.727589,0.765229,0.834458,0.878155,0.824003,0.893187,0.824451,0.76873,0.796163,0.713176,0.736468,0.722194,0.79335,0.892695,0.81785,0.702785,0.81674,0.859769,0.845479,0.729515,0.848578,0.831617,0.782084,0.817316,0.873488,0.849991,0.717202,0.827419,0.885216,0.870209,0.801342,0.814752,0.867845,0.866332,0.738925,0.778826,0.793383,0.812433,0.718611,0.704671,0.889181};
    
    //Selections
    if (nTPCcr<nTPCcr_min[isyst])                     return passedTrkSelection;
    if (nITScls<nITScls_min[isyst])                   return passedTrkSelection;
    if (nTPCclsdEdx<nTPCclsdEdx_min[isyst])           return passedTrkSelection;
    if (cr_over_findable<cr_over_findable_min[isyst]) return passedTrkSelection;
    if (chi2TPC_NDF>chi2TPC_NDF_max[isyst])           return passedTrkSelection;
    if (chi2ITS>chi2ITS_max[isyst])                   return passedTrkSelection;
    if (DCAz>dcaz_max[isyst])                         return passedTrkSelection;
    
    passedTrkSelection = kTRUE;
    return passedTrkSelection;
}
//__________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskAntiProtons_vs_Rapidity_Data::IsProtonCandidate (AliESDtrack *track)  {
    
    //Initialization
    Bool_t isProton=(kFALSE);
    
    Double_t DCAxy = GetDCAtoPrimaryVertex (track,0);

    if (!PassedTrackSelection(track,0)) return isProton;
    if (!IsHighPurityProton(track))     return isProton;
    if (TMath::Abs(DCAxy)>0.1)          return isProton;
    if (track->P()>1.0)                 return isProton;
     
    isProton=kTRUE;
    return isProton;
}
//__________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskAntiProtons_vs_Rapidity_Data::IsHighPurityProton (AliESDtrack *track)  {
    
    //Variables
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kProton);
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kProton);
    Double_t nsigmaITS = fPIDResponse -> NumberOfSigmasITS (track,AliPID::kProton);
    Double_t nsigmaITS_recalib = GetRecalibratedITSnsigma (nsigmaITS,track->Eta(),track->P());
    Double_t pt = track->Pt();
    
    //Selection
    if (pt<0.7 && TMath::Abs(nsigmaITS_recalib)<3.0 && TMath::Abs(nsigmaTPC)<3.0) return kTRUE;
    if (pt>0.7 && TMath::Abs(nsigmaTPC)<2.0 && TMath::Abs(nsigmaTOF)<2.0)         return kTRUE;

    return kFALSE;
}
//__________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskAntiProtons_vs_Rapidity_Data::GetDCAtoPrimaryVertex (AliESDtrack *track, Int_t index)  {
         
    Double_t dca[2];
    Double_t covMatrix[3];
    if (!track->PropagateToDCA (fESDEvent->GetPrimaryVertex(),fESDEvent->GetMagneticField(),10000,dca,covMatrix)) return -999;
                        
    return dca[index];
}
//__________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskAntiProtons_vs_Rapidity_Data::GetRapidity (AliESDtrack *track, Double_t mass)  {
    
    //Initialization
    Double_t y(-999);
    
    //Rapidity Calculation
    Double_t p  = track->P();
    Double_t pz = track->Pz();
    Double_t E = TMath::Sqrt(mass*mass + p*p);
    if (E == TMath::Abs(pz)) return -999;
    y = 0.5*TMath::Log((E+pz)/(E-pz));
    
    return y;
}
//____________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskAntiProtons_vs_Rapidity_Data::GetRecalibratedITSnsigma (Double_t nsigma, Double_t eta, Double_t p)  {
    
    //Initialization
    Double_t nsigma_corr=nsigma;
    
    //Protections
    if (eta<-1.0) return nsigma_corr;
    if (eta>+1.0) return nsigma_corr;
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
//__________________________________________________________________________________________________________________________________
void AliAnalysisTaskAntiProtons_vs_Rapidity_Data::Terminate(Option_t *)  {
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
//__________________________________________________________________________________________________________________________________
