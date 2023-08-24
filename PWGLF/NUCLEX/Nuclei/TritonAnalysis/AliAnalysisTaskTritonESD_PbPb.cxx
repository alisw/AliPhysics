#include "AliAnalysisTaskTritonESD_PbPb.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliAnalysisTask.h"
#include "TLorentzVector.h"
#include "AliPIDResponse.h"
#include "AliCentrality.h"
#include "TDatabasePDG.h"
#include "AliEventCuts.h"
#include <AliESDEvent.h>
#include <AliESDInputHandler.h>
#include <AliESDtrack.h>
#include <AliESDtrackCuts.h>
#include "TObjArray.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TChain.h"
#include "TMath.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"

ClassImp(AliAnalysisTaskTritonESD_PbPb)

//_________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskTritonESD_PbPb::AliAnalysisTaskTritonESD_PbPb():
AliAnalysisTaskSE(),
fESDevent(NULL),
fPIDResponse(NULL),
fESDeventCuts(),
fUtils(NULL),
fOutputList(NULL),
reducedTree_Triton(NULL),
//fQAList(NULL),
fCentralityMin(0),
fCentralityMax(0),
fVertexZmin(0),
fVertexZmax(0),
fNumberVertexContributorsMin(0),
fCentralityEstimator(NULL),
fPtMin(0),
fPtMax(0),
fEtaMax(0),
fYMax(0),
fNumberClustersITSMin(0),
fNumberClustersTPCMin(0),
fNumberCrossedRowsTPCMin(0),
fCrossedRowsFindableClsMin(0),
fNumberClustersTPCdEdxMin(0),
fChiSquarePerNDFMax(0),
fITSrequirement(NULL),
fDCAzMax(0),
fDCAxyMax(0),
fnSigmaTOFmax(0),
fnSigmaTPCmax(0),
fTRDntracklets(0),
fpar0_mean_TPC(0),
fpar1_mean_TPC(0),
fpar0_sigma_TPC(0),
fpar0_mean_TOF(0),
fpar1_mean_TOF(0),
fpar0_sigma_TOF(0),
fpar1_sigma_TOF(0),
multPercentile_V0M(-1),
pt(0),
p(0),
eta(0),
y(0),
q(0),
dcaxy(0),
dcaz(0),
nTPC_Clusters(0),
nTRD_Clusters(0),
nITS_Clusters(0),
nTPC_FindableClusters(0),
nTPC_CrossedRows(0),
nTPC_Clusters_dEdx(0),
HasPointOnITSLayer0(0),
HasPointOnITSLayer1(0),
HasSharedPointOnITSLayer0(0),
HasSharedPointOnITSLayer1(0),
chi2_TPC(0),
chi2_NDF(0),
chi2_ITS(0),
ITSsignal(0),
TPCsignal(0),
TOFsignal(0),
TRDsignal(0),
nSigmaTPC_Trit(0),
nSigmaTOF_Trit(0),
nSigmaTRD_Trit(0)
{}
//_________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskTritonESD_PbPb::AliAnalysisTaskTritonESD_PbPb(const char *name):
AliAnalysisTaskSE(name),
fESDevent(NULL),
fPIDResponse(NULL),
fESDeventCuts(),
fUtils(NULL),
fOutputList(NULL),
reducedTree_Triton(NULL),
//fQAList(NULL),
fCentralityMin(0),
fCentralityMax(0),
fVertexZmin(0),
fVertexZmax(0),
fNumberVertexContributorsMin(0),
fCentralityEstimator(NULL),
fPtMin(0),
fPtMax(0),
fEtaMax(0),
fYMax(0),
fNumberClustersITSMin(0),
fNumberClustersTPCMin(0),
fNumberCrossedRowsTPCMin(0),
fCrossedRowsFindableClsMin(0),
fNumberClustersTPCdEdxMin(0),
fChiSquarePerNDFMax(0),
fITSrequirement(NULL),
fDCAzMax(0),
fDCAxyMax(0),
fnSigmaTOFmax(0),
fnSigmaTPCmax(0),
fTRDntracklets(0),
fpar0_mean_TPC(0),
fpar1_mean_TPC(0),
fpar0_sigma_TPC(0),
fpar0_mean_TOF(0),
fpar1_mean_TOF(0),
fpar0_sigma_TOF(0),
fpar1_sigma_TOF(0),
multPercentile_V0M(-1),
pt(0),
p(0),
eta(0),
y(0),
q(0),
dcaxy(0),
dcaz(0),
nTPC_Clusters(0),
nTRD_Clusters(0),
nITS_Clusters(0),
nTPC_FindableClusters(0),
nTPC_CrossedRows(0),
nTPC_Clusters_dEdx(0),
HasPointOnITSLayer0(0),
HasPointOnITSLayer1(0),
HasSharedPointOnITSLayer0(0),
HasSharedPointOnITSLayer1(0),
chi2_TPC(0),
chi2_NDF(0),
chi2_ITS(0),
ITSsignal(0),
TPCsignal(0),
TOFsignal(0),
TRDsignal(0),
nSigmaTPC_Trit(0),
nSigmaTOF_Trit(0),
nSigmaTRD_Trit(0)
{
    fUtils = new AliAnalysisUtils();
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TTree::Class());
    //DefineOutput(2, TList::Class());
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskTritonESD_PbPb::~AliAnalysisTaskTritonESD_PbPb()
{
    fOutputList->Clear();
    delete fESDevent;
    delete fPIDResponse;
    delete fUtils;
    delete fOutputList;
    delete reducedTree_Triton;
    //delete fQAList;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskTritonESD_PbPb::UserCreateOutputObjects()
{
    fOutputList = new TList();
    fOutputList -> SetOwner();
   
    //fQAList = new TList();
    //fQAList -> SetOwner();
    
    //fESDeventCuts.AddQAplotsToList(fQAList);//Add event selection QA plots

    
    //Number of Events
    histoNumberOfEvents = new TH1F("histoNumberOfEvents","Events after selection steps",10,0,10);
    fOutputList -> Add (histoNumberOfEvents);
    
    
    //Signal Extraction
    histoNsigmaTPCtriton_vs_pt     = new TH2F ("histoNsigmaTPCtriton_vs_pt","",500,0,5,1000,-20,20);
    histoNsigmaTOFtriton_vs_pt     = new TH2F ("histoNsigmaTOFtriton_vs_pt","",500,0,5,1000,-20,20);
    histoNsigmaTPCantitriton_vs_pt = new TH2F ("histoNsigmaTPCantitriton_vs_pt","",500,0,5,1000,-20,20);
    histoNsigmaTOFantitriton_vs_pt = new TH2F ("histoNsigmaTOFantitriton_vs_pt","",500,0,5,1000,-20,20);
 
    histoNsigmaTPCtriton_vs_pt_centered     = new TH2F ("histoNsigmaTPCtriton_vs_pt_centered","",500,0,5,1000,-20,20);
    histoNsigmaTPCantitriton_vs_pt_centered = new TH2F ("histoNsigmaTPCantitriton_vs_pt_centered","",500,0,5,1000,-20,20);
    
    histoNsigmaTOFtriton_vs_pt_centered     = new TH2F ("histoNsigmaTOFtriton_vs_pt_centered","",500,0,5,1000,-20,20);
    histoNsigmaTOFantitriton_vs_pt_centered = new TH2F ("histoNsigmaTOFantitriton_vs_pt_centered","",500,0,5,1000,-20,20);
    
    histoNsigmaTOFtriton_vs_pt_trd     = new TH2F ("histoNsigmaTOFtriton_vs_pt_trd","",500,0,5,1000,-20,20);
    histoNsigmaTOFantitriton_vs_pt_trd = new TH2F ("histoNsigmaTOFantitriton_vs_pt_trd","",500,0,5,1000,-20,20);
        
    histoNsigmaTPCtriton_vs_p         = new TH2F ("histoNsigmaTPCtriton_vs_p","",500,0,5,1000,-20,20);
    histoNsigmaTPCantitriton_vs_p     = new TH2F ("histoNsigmaTPCantitriton_vs_p","",500,0,5,1000,-20,20);
    histoNsigmaTOFtriton_vs_p         = new TH2F ("histoNsigmaTOFtriton_vs_p","",500,0,5,1000,-20,20);
    histoNsigmaTOFantitriton_vs_p     = new TH2F ("histoNsigmaTOFantitriton_vs_p","",500,0,5,1000,-20,20);
    
    histoNsigmaTPCtriton_vs_p_notof     = new TH2F ("histoNsigmaTPCtriton_vs_p_notof","",500,0,5,1000,-20,20);
    histoNsigmaTPCantitriton_vs_p_notof = new TH2F ("histoNsigmaTPCantitriton_vs_p_notof","",500,0,5,1000,-20,20);

    histoNsigmaTPCtriton_vs_pt     -> Sumw2();
    histoNsigmaTOFtriton_vs_pt     -> Sumw2();
    histoNsigmaTPCantitriton_vs_pt -> Sumw2();
    histoNsigmaTOFantitriton_vs_pt -> Sumw2();
  
    histoNsigmaTPCtriton_vs_pt_centered           -> Sumw2();
    histoNsigmaTPCantitriton_vs_pt_centered       -> Sumw2();
    
    histoNsigmaTOFtriton_vs_pt_centered           -> Sumw2();
    histoNsigmaTOFantitriton_vs_pt_centered       -> Sumw2();
    
    histoNsigmaTOFtriton_vs_pt_trd                -> Sumw2();
    histoNsigmaTOFantitriton_vs_pt_trd            -> Sumw2();
    
    histoNsigmaTPCtriton_vs_p                     -> Sumw2();
    histoNsigmaTPCantitriton_vs_p                 -> Sumw2();
    histoNsigmaTOFtriton_vs_p                     -> Sumw2();
    histoNsigmaTOFantitriton_vs_p                 -> Sumw2();
    
    histoNsigmaTPCtriton_vs_p_notof              -> Sumw2();
	histoNsigmaTPCantitriton_vs_p_notof          -> Sumw2();
    
    fOutputList -> Add(histoNsigmaTPCtriton_vs_pt);
    fOutputList -> Add(histoNsigmaTOFtriton_vs_pt);
    fOutputList -> Add(histoNsigmaTPCantitriton_vs_pt);
    fOutputList -> Add(histoNsigmaTOFantitriton_vs_pt);
  
    fOutputList -> Add(histoNsigmaTPCtriton_vs_pt_centered);
    fOutputList -> Add(histoNsigmaTPCantitriton_vs_pt_centered);
    fOutputList -> Add(histoNsigmaTOFtriton_vs_pt_centered);
    fOutputList -> Add(histoNsigmaTOFantitriton_vs_pt_centered);
    
    fOutputList -> Add(histoNsigmaTOFtriton_vs_pt_trd);
    fOutputList -> Add(histoNsigmaTOFantitriton_vs_pt_trd);
        
    fOutputList -> Add(histoNsigmaTPCtriton_vs_p);
    fOutputList -> Add(histoNsigmaTPCantitriton_vs_p);
    fOutputList -> Add(histoNsigmaTOFtriton_vs_p);
    fOutputList -> Add(histoNsigmaTOFantitriton_vs_p);
    
    fOutputList -> Add(histoNsigmaTPCtriton_vs_p_notof);
    fOutputList -> Add(histoNsigmaTPCantitriton_vs_p_notof);

    
    
    
    //DCA Distributions
    histoDCAxyTriton_vs_pt     = new TH2F ("histoDCAxyTriton_vs_pt","",500,0,5,500,-5,5);
    histoDCAxyAntiTriton_vs_pt = new TH2F ("histoDCAxyAntiTriton_vs_pt","",500,0,5,500,-5,5);

    histoDCAxyTriton_vs_pt     -> Sumw2();
    histoDCAxyAntiTriton_vs_pt -> Sumw2();
    
    fOutputList -> Add (histoDCAxyTriton_vs_pt);
    fOutputList -> Add (histoDCAxyAntiTriton_vs_pt);

    
    //Reduced Tree (Triton)
    reducedTree_Triton = new TTree("reducedTree_Triton","reducedTree_Triton");
    reducedTree_Triton -> Branch("multPercentile_V0M",&multPercentile_V0M,"multPercentile_V0M/D");
    reducedTree_Triton -> Branch("pt",&pt,"pt/D");
    reducedTree_Triton -> Branch("p",&p,"p/D");
    reducedTree_Triton -> Branch("eta",&eta,"eta/D");
    reducedTree_Triton -> Branch("y",&y,"y/D");
    reducedTree_Triton -> Branch("q",&q,"q/I");
    reducedTree_Triton -> Branch("dcaxy",&dcaxy,"dcaxy/D");
    reducedTree_Triton -> Branch("dcaz",&dcaz,"dcaz/D");
    reducedTree_Triton -> Branch("nTPC_Clusters",&nTPC_Clusters,"nTPC_Clusters/I");
    reducedTree_Triton -> Branch("nTRD_Clusters",&nTRD_Clusters,"nTRD_Clusters/I");
    reducedTree_Triton -> Branch("nITS_Clusters",&nITS_Clusters,"nITS_Clusters/I");
    reducedTree_Triton -> Branch("nTPC_FindableClusters",&nTPC_FindableClusters,"nTPC_FindableClusters/I");
    reducedTree_Triton -> Branch("nTPC_CrossedRows",&nTPC_CrossedRows,"nTPC_CrossedRows/I");
    reducedTree_Triton -> Branch("nTPC_Clusters_dEdx",&nTPC_Clusters_dEdx,"nTPC_Clusters_dEdx/I");
    reducedTree_Triton -> Branch("HasPointOnITSLayer0",&HasPointOnITSLayer0,"HasPointOnITSLayer0/O");
    reducedTree_Triton -> Branch("HasPointOnITSLayer1",&HasPointOnITSLayer1,"HasPointOnITSLayer1/O");
    reducedTree_Triton -> Branch("HasSharedPointOnITSLayer0",&HasSharedPointOnITSLayer0,"HasSharedPointOnITSLayer0/O");
    reducedTree_Triton -> Branch("HasSharedPointOnITSLayer1",&HasSharedPointOnITSLayer1,"HasSharedPointOnITSLayer1/O");
    reducedTree_Triton -> Branch("chi2_TPC",&chi2_TPC,"chi2_TPC/D");
    reducedTree_Triton -> Branch("chi2_ITS",&chi2_ITS,"chi2_ITS/D");
    reducedTree_Triton -> Branch("ITSsignal",&ITSsignal,"ITSsignal/D");
    reducedTree_Triton -> Branch("TPCsignal",&TPCsignal,"TPCsignal/D");
    reducedTree_Triton -> Branch("TOFsignal",&TOFsignal,"TOFsignal/D");
    reducedTree_Triton -> Branch("TRDsignal",&TRDsignal,"TRDsignal/D");
    reducedTree_Triton -> Branch("nSigmaTPC_Trit",&nSigmaTPC_Trit,"nSigmaTPC_Trit/D");
    reducedTree_Triton -> Branch("nSigmaTOF_Trit",&nSigmaTOF_Trit,"nSigmaTOF_Trit/D");
    reducedTree_Triton -> Branch("nSigmaTRD_Trit",&nSigmaTRD_Trit,"nSigmaTRD_Trit/D");
  

    
    PostData(1, fOutputList);
    PostData(2,reducedTree_Triton);

    //PostData(2, fQAList);
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskTritonESD_PbPb::UserExec(Option_t *)
{
    
    //Get Input Event
    if ( !GetInputEvent ()) return;
    
    //Load PID Response
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();

    
    //Loop over Reconstructed Tracks
    for (Int_t i=0 ; i<fESDevent->GetNumberOfTracks() ; i++)  {
        
        //Track Selection
        AliESDtrack *track = (AliESDtrack*) fESDevent -> GetTrack(i);
        if ( !track ) continue;
        if ( PassedTrackQualityCutsNoDCA (track)) {

        if (IsCleanTritonCandidate(track))  {
            if (track->Charge()>0) histoDCAxyTriton_vs_pt     -> Fill (track->Pt(),GetDCAxy(track));
            if (track->Charge()<0) histoDCAxyAntiTriton_vs_pt -> Fill (track->Pt(),GetDCAxy(track));
        }
        
        }

        if ( !PassedTrackQualityCuts (track)) continue;

        
        //Variables
        Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kTriton);
        Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kTriton);

        
        //TPC Signal vs. pT
        
        
      	if (track->Charge()>0) histoNsigmaTPCtriton_vs_p_notof	    -> Fill (track->P(),nsigmaTPC);
        if (track->Charge()<0) histoNsigmaTPCantitriton_vs_p_notof -> Fill (track->P(),nsigmaTPC);
        
        
        if (PassedTOFSelection(track))  {
            
            if (track->Charge()>0) histoNsigmaTPCtriton_vs_pt     -> Fill (track->Pt(),nsigmaTPC);
            if (track->Charge()<0) histoNsigmaTPCantitriton_vs_pt -> Fill (track->Pt(),nsigmaTPC);
	        if (track->Charge()>0) histoNsigmaTPCtriton_vs_p      -> Fill (track->P(),nsigmaTPC);
            if (track->Charge()<0) histoNsigmaTPCantitriton_vs_p  -> Fill (track->P(),nsigmaTPC);
        
            if (track->Charge()>0) histoNsigmaTPCtriton_vs_pt_centered     -> Fill (track->Pt(),Centered_nsigmaTPC(track));
            if (track->Charge()<0) histoNsigmaTPCantitriton_vs_pt_centered -> Fill (track->Pt(),Centered_nsigmaTPC(track));
        
        
        }
        
        //TOF Signal vs. pT
        if (PassedTPCSelection(track))  {
            
            if (track->Charge()>0) histoNsigmaTOFtriton_vs_pt     -> Fill (track->Pt(),nsigmaTOF);
            if (track->Charge()<0) histoNsigmaTOFantitriton_vs_pt -> Fill (track->Pt(),nsigmaTOF);
            if (track->Charge()>0) histoNsigmaTOFtriton_vs_p      -> Fill (track->P(),nsigmaTOF);
            if (track->Charge()<0) histoNsigmaTOFantitriton_vs_p  -> Fill (track->P(),nsigmaTOF);
        
            
            if (track->Charge()>0 && track->GetTRDntrackletsPID()>fTRDntracklets) histoNsigmaTOFtriton_vs_pt_trd     -> Fill (track->Pt(),nsigmaTOF);
            if (track->Charge()<0 && track->GetTRDntrackletsPID()>fTRDntracklets) histoNsigmaTOFantitriton_vs_pt_trd -> Fill (track->Pt(),nsigmaTOF);
            
            if (track->Charge()>0) histoNsigmaTOFtriton_vs_pt_centered     -> Fill (track->Pt(),Centered_nsigmaTOF(track));
            if (track->Charge()<0) histoNsigmaTOFantitriton_vs_pt_centered -> Fill (track->Pt(),Centered_nsigmaTOF(track));
        
            
            
        }

        pt = track -> Pt();
        p = track -> P();
        
        q  = (Int_t) track -> Charge();
        eta = track -> Eta();
        
        //Rapidity Calculation
        Double_t m  = AliPID::ParticleMass(AliPID::kTriton);
        Double_t px = track -> Px();
        Double_t py = track -> Py();
        Double_t pz = track -> Pz();
        Double_t E = TMath::Sqrt(m*m + px*px + py*py + pz*pz);
        TLorentzVector P (px,py,pz,E);
        y = P.Rapidity();

        //DCA
        dcaz  = GetDCAz  (track);
        dcaxy = GetDCAxy (track);

        nTPC_Clusters = track->GetTPCNcls();
        nTRD_Clusters = track->GetTRDncls();
        nITS_Clusters = track->GetITSNcls();

        nTPC_FindableClusters = track->GetTPCNclsF();
        nTPC_CrossedRows = track->GetTPCCrossedRows();
        nTPC_Clusters_dEdx = track -> GetTPCsignalN();

        HasPointOnITSLayer0 = track->HasPointOnITSLayer(0);
        HasPointOnITSLayer1 = track->HasPointOnITSLayer(1);

        HasSharedPointOnITSLayer0 = track->HasSharedPointOnITSLayer(0);
        HasSharedPointOnITSLayer1 = track->HasSharedPointOnITSLayer(1);

        chi2_TPC = track -> GetTPCchi2();//    check -> seems to be 0
        chi2_ITS = track -> GetITSchi2();//    check

        ITSsignal = track->GetITSsignal();
        TPCsignal = track->GetTPCsignal();
        TOFsignal = track->GetTOFsignal();
        TRDsignal = track->GetTRDsignal();

        nSigmaTPC_Trit = fPIDResponse -> NumberOfSigmasTPC(track,AliPID::kTriton);
        nSigmaTOF_Trit = fPIDResponse -> NumberOfSigmasTOF(track,AliPID::kTriton);
        nSigmaTRD_Trit = fPIDResponse -> NumberOfSigmasTRD(track,AliPID::kTriton);

        if(IsTritonCandidate(track)) reducedTree_Triton -> Fill();
    }
    
    
    PostData(1, fOutputList);
    PostData(2,reducedTree_Triton);
    //PostData(2, fQAList);
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskTritonESD_PbPb::GetInputEvent ()  {
    
    //Get Input Event
    fESDevent = dynamic_cast <AliESDEvent*>(InputEvent());
    if (!fESDevent) return false;
    histoNumberOfEvents -> Fill(0.5);
    
    //Standard Event Cuts
    if (!fESDeventCuts.AcceptEvent(fESDevent)) {
        //PostData(2, fQAList);
        return false;
    }
    histoNumberOfEvents -> Fill(1.5);
    
    //Centrality
    AliMultSelection *multiplicitySelection = (AliMultSelection*) fESDevent->FindListObject("MultSelection");
    if( !multiplicitySelection) return false;
    histoNumberOfEvents -> Fill(2.5);
    Double_t centrality = multiplicitySelection->GetMultiplicityPercentile(fCentralityEstimator);
    
   
    //Selection of Centrality Range
    if (centrality<fCentralityMin || centrality>=fCentralityMax ) return false;
    histoNumberOfEvents -> Fill(3.5);
    
    //Primary Vertex
    AliESDVertex *vertex = (AliESDVertex*) fESDevent->GetPrimaryVertex();
    if ( !vertex ) return false;
    histoNumberOfEvents -> Fill(4.5);
    
    //Primary Vertex Selection
    if ( vertex->GetZ() < fVertexZmin ) return false;
    if ( vertex->GetZ() > fVertexZmax ) return false;
    histoNumberOfEvents -> Fill(5.5);
    
    if ( vertex->GetNContributors() < fNumberVertexContributorsMin ) return false;
    histoNumberOfEvents -> Fill(6.5);
    
    multPercentile_V0M              = multiplicitySelection->GetMultiplicityPercentile("V0M");

    // GetPrimaryVertex
    AliESDVertex *vertex_tracks = (AliESDVertex*) fESDevent->GetPrimaryVertexTracks();
    if (!vertex_tracks) return false;
    if ( vertex_tracks->GetNContributors() < 1 ) return false;
    
    //Primary Vertex SPD
    AliESDVertex *vertex_SPD = (AliESDVertex*) fESDevent->GetPrimaryVertexSPD();
    if (!vertex_SPD) return false;
    
    //Vertex Contributors SPD
    if ( vertex_SPD->GetNContributors() < 1 ) return false;
    //SPD Pile-up in Mult Bins
    if (fESDevent->IsPileupFromSPDInMultBins()) return false;
    
    //Primary Vertex Selection
    if ( vertex_tracks->GetZ() < -10.0 ) return false;
    if ( vertex_tracks->GetZ() > +10.0 ) return false;

    return true;
    
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskTritonESD_PbPb::PassedTrackQualityCuts (AliESDtrack* track)  {
    
    //Initialization
    Bool_t passedTrkSelection=(kFALSE);
    
    if ( track->GetTPCsignalN() < 50 )         return passedTrkSelection;
    
    //Rapidity Calculation
    Double_t m  = AliPID::ParticleMass(AliPID::kTriton);
    Double_t px = track -> Px();
    Double_t py = track -> Py();
    Double_t pz = track -> Pz();
    Double_t E = TMath::Sqrt(m*m + px*px + py*py + pz*pz);
    TLorentzVector P (px,py,pz,E);
    Double_t y = P.Rapidity();
    
    
    //Kinematic Cuts & Acceptance
    if ( track->Pt()<fPtMin || track->Pt()>fPtMax ) return false;
    if ( TMath::Abs(track->Eta()) > fEtaMax )       return false;
    if ( TMath::Abs(y) > fYMax )       return false;
    
    //Track Selection Cuts
    if ( track->GetITSNcls() < fNumberClustersITSMin ) return false;
    if ( track->GetTPCNcls() < fNumberClustersTPCMin ) return false;
    if ( track->GetTPCCrossedRows() < fNumberCrossedRowsTPCMin ) return false;
    if ( static_cast<Double_t>(track->GetTPCCrossedRows())/static_cast<Double_t>(track->GetTPCNclsF()) < fCrossedRowsFindableClsMin) return false;
    if ( track->GetTPCsignalN() < fNumberClustersTPCdEdxMin ) return false;
    if ( (((Double_t)track->GetTPCchi2())/((Double_t)track->GetTPCNcls())) > fChiSquarePerNDFMax) return false;
    
    //ITS Requirement
    Bool_t hitInITSLayer0 = track->HasPointOnITSLayer(0);
    Bool_t hitInITSLayer1 = track->HasPointOnITSLayer(1);
    

    if (strcmp(fITSrequirement,"kBoth")==0   && !hitInITSLayer0 ) return false;
    if (strcmp(fITSrequirement,"kBoth")==0   && !hitInITSLayer1 ) return false;
    if (strcmp(fITSrequirement,"kFirst")==0  && !hitInITSLayer0 ) return false;
    if (strcmp(fITSrequirement,"kSecond")==0 && !hitInITSLayer1 ) return false;
    if (strcmp(fITSrequirement,"kAny")==0    && (!hitInITSLayer0) && (!hitInITSLayer1)) return false;

    
    //DCA Cuts
    Double_t dcaxy = GetDCAxy (track);
    Double_t dcaz  = GetDCAz (track);
    if (TMath::Abs(dcaxy) > fDCAxyMax) return false;
    if (TMath::Abs(dcaz)  > fDCAzMax)  return false;

    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskTritonESD_PbPb::PassedTrackQualityCutsNoDCA (AliESDtrack* track)  {
    
    //Initialization
    Bool_t passedTrkSelection=(kFALSE);
    
    if ( track->GetTPCsignalN() < 50 )         return passedTrkSelection;
    
    
    //Rapidity Calculation
    Double_t m  = AliPID::ParticleMass(AliPID::kTriton);
    Double_t px = track -> Px();
    Double_t py = track -> Py();
    Double_t pz = track -> Pz();
    Double_t E = TMath::Sqrt(m*m + px*px + py*py + pz*pz);
    TLorentzVector P (px,py,pz,E);
    Double_t y = P.Rapidity();
    
    
    //Kinematic Cuts & Acceptance
    if ( track->Pt()<fPtMin || track->Pt()>fPtMax ) return false;
    if ( TMath::Abs(track->Eta()) > fEtaMax )       return false;
    if ( TMath::Abs(y) > fYMax )       return false;
    
    //Track Selection Cuts
    if ( track->GetITSNcls() < fNumberClustersITSMin ) return false;
    if ( track->GetTPCNcls() < fNumberClustersTPCMin ) return false;
    if ( track->GetTPCCrossedRows() < fNumberCrossedRowsTPCMin ) return false;
    if ( static_cast<Double_t>(track->GetTPCCrossedRows())/static_cast<Double_t>(track->GetTPCNclsF()) < fCrossedRowsFindableClsMin) return false;
    if ( track->GetTPCsignalN() < fNumberClustersTPCdEdxMin ) return false;
    if ( (((Double_t)track->GetTPCchi2())/((Double_t)track->GetTPCNcls())) > fChiSquarePerNDFMax) return false;
    
    //ITS Requirement
    Bool_t hitInITSLayer0 = track->HasPointOnITSLayer(0);
    Bool_t hitInITSLayer1 = track->HasPointOnITSLayer(1);
    

    if (strcmp(fITSrequirement,"kBoth")==0   && !hitInITSLayer0 ) return false;
    if (strcmp(fITSrequirement,"kBoth")==0   && !hitInITSLayer1 ) return false;
    if (strcmp(fITSrequirement,"kFirst")==0  && !hitInITSLayer0 ) return false;
    if (strcmp(fITSrequirement,"kSecond")==0 && !hitInITSLayer1 ) return false;
    if (strcmp(fITSrequirement,"kAny")==0    && (!hitInITSLayer0) && (!hitInITSLayer1)) return false;

    
    return true;
}
//_______________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskTritonESD_PbPb::IsTritonCandidate (AliESDtrack *track)  {

  Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC(track,AliPID::kTriton);
  if (track->Pt() < 1.5) {
    if ( TMath::Abs(nsigmaTPC) > 6. ) return false;
  } else {
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF(track,AliPID::kTriton);
    if ( TMath::Abs(nsigmaTOF) > 7. ) return false;
    if ( TMath::Abs(nsigmaTPC) > 6. ) return false;
  }

  return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskTritonESD_PbPb::IsCleanTritonCandidate (AliESDtrack *track)  {
    
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kTriton);
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kTriton);
    if (TMath::Abs(nsigmaTOF) > fnSigmaTOFmax) return false;
    if (TMath::Abs(nsigmaTPC) > fnSigmaTPCmax) return false;

    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskTritonESD_PbPb::Centered_nsigmaTPC (AliESDtrack *track)  {
   
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kTriton);
   
    Double_t mean_fitted  = fpar0_mean_TPC*exp(fpar1_mean_TPC*(track->P()));
    Double_t sigma_fitted = fpar0_sigma_TPC*(track->P());
    
    nsigmaTPC = (nsigmaTPC - mean_fitted)/sigma_fitted;
   
    return nsigmaTPC;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskTritonESD_PbPb::Centered_nsigmaTOF (AliESDtrack *track)  {
   
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kTriton);
   
    Double_t mean_fitted  = fpar0_mean_TOF*exp(fpar1_mean_TOF*(track->P()));
    Double_t sigma_fitted = fpar0_sigma_TOF*exp(fpar1_sigma_TOF*(track->P()));
    
    nsigmaTOF = (nsigmaTOF - mean_fitted)/sigma_fitted;
   
    return nsigmaTOF;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskTritonESD_PbPb::PassedTOFSelection (AliESDtrack *track)  {
    
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kTriton);
    if (TMath::Abs(nsigmaTOF) > fnSigmaTOFmax) return false;
    
    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskTritonESD_PbPb::PassedTPCSelection (AliESDtrack *track)  {
    
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kTriton);
    if (TMath::Abs(nsigmaTPC) > fnSigmaTPCmax) return false;
    if ((track->GetStatus()&AliESDtrack::kTOFout)==0) return false;

    //TPC-TOF Matching
//   if ((track->GetStatus()&AliESDtrack::kTOFout)==0) hasTOFhit=0;//Track with no TOF hit
//   if ((track->GetStatus()&AliESDtrack::kTOFout)!=0) hasTOFhit=1;//Track with TOF hit
    
    return true;
}

//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskTritonESD_PbPb::GetDCAxy (AliESDtrack *track)  {
    
    Double_t impactParameter[2];
    Double_t covarianceMatrix[3];
    if (!track->PropagateToDCA(fESDevent->GetPrimaryVertex(),fESDevent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;
    
    Double_t DCAxy = impactParameter[0];
    
    return DCAxy;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskTritonESD_PbPb::GetDCAz (AliESDtrack *track)  {
    
    Double_t impactParameter[2];
    Double_t covarianceMatrix[3];
    if (!track->PropagateToDCA(fESDevent->GetPrimaryVertex(),fESDevent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;
    
    Double_t DCAz = impactParameter[1];
    
    return DCAz;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskTritonESD_PbPb::Terminate(Option_t *)  {
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________

