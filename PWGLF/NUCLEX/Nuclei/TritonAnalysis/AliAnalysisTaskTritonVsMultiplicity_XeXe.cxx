#include "AliAnalysisTaskTritonVsMultiplicity_XeXe.h"
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
#include "AliAODVertex.h"
#include "AliEventCuts.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "TObjArray.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TRandom.h"
#include "TChain.h"
#include "TMath.h"
#include "TList.h"
#include "TH1F.h"
#include "TH2F.h"

ClassImp(AliAnalysisTaskTritonVsMultiplicity_XeXe)

//_________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskTritonVsMultiplicity_XeXe::AliAnalysisTaskTritonVsMultiplicity_XeXe():
AliAnalysisTaskSE(),
fAODevent(nullptr),
fPIDResponse(nullptr),
fAODeventCuts(),
fUtils(nullptr),
fOutputList(nullptr),
fCentralityMin(0),
fCentralityMax(0),
fVertexZmin(0),
fVertexZmax(0),
fNumberVertexContributorsMin(0),
fCentralityEstimator(nullptr),
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
fITSrequirement(nullptr),
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
multPercentile_V0A(0),
pt(0),
p(0),
px(0),
py(0),
pz(0),
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
HasPointOnITSLayer2(0),
HasPointOnITSLayer3(0),
HasPointOnITSLayer4(0),
HasPointOnITSLayer5(0),
chi2_TPC(0),
chi2_ITS(0),
ITSsignal(0),
TPCsignal(0),
TOFsignal(0),
nSigmaITS_Trit(0),
nSigmaTPC_Trit(0),
nSigmaTOF_Trit(0)
{}
//_________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskTritonVsMultiplicity_XeXe::AliAnalysisTaskTritonVsMultiplicity_XeXe(const char *name):
AliAnalysisTaskSE(name),
fAODevent(nullptr),
fPIDResponse(nullptr),
fAODeventCuts(),
fUtils(nullptr),
fOutputList(nullptr),
//fQAList(nullptr),
fCentralityMin(0),
fCentralityMax(0),
fVertexZmin(0),
fVertexZmax(0),
fNumberVertexContributorsMin(0),
fCentralityEstimator(nullptr),
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
fITSrequirement(nullptr),
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
multPercentile_V0A(0),
pt(0),
p(0),
px(0),
py(0),
pz(0),
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
HasPointOnITSLayer2(0),
HasPointOnITSLayer3(0),
HasPointOnITSLayer4(0),
HasPointOnITSLayer5(0),
chi2_TPC(0),
chi2_ITS(0),
ITSsignal(0),
TPCsignal(0),
TOFsignal(0),
nSigmaITS_Trit(0),
nSigmaTPC_Trit(0),
nSigmaTOF_Trit(0)
{
    fUtils = new AliAnalysisUtils();
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TTree::Class());
    //DefineOutput(2, TList::Class());
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskTritonVsMultiplicity_XeXe::~AliAnalysisTaskTritonVsMultiplicity_XeXe()
{
    fOutputList->Clear();
    delete fAODevent;
    delete fPIDResponse;
    delete fUtils;
    delete fOutputList;
    delete reducedTree_Triton;
    //delete fQAList;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskTritonVsMultiplicity_XeXe::UserCreateOutputObjects()
{
    fOutputList = new TList();
    fOutputList -> SetOwner();

    //fQAList = new TList();
    //fQAList -> SetOwner();

    //fAODeventCuts.AddQAplotsToList(fQAList);//Add event selection QA plots


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
    reducedTree_Triton -> Branch("multPercentile_V0A",&multPercentile_V0A,"multPercentile_V0A/D");
    reducedTree_Triton -> Branch("pt",&pt,"pt/D");
    reducedTree_Triton -> Branch("p",&p,"p/D");
    reducedTree_Triton -> Branch("px",&px,"px/D");
    reducedTree_Triton -> Branch("py",&py,"py/D");
    reducedTree_Triton -> Branch("pz",&pz,"pz/D");
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
    reducedTree_Triton -> Branch("HasPointOnITSLayer2",&HasPointOnITSLayer2,"HasPointOnITSLayer2/O");
    reducedTree_Triton -> Branch("HasPointOnITSLayer3",&HasPointOnITSLayer3,"HasPointOnITSLayer3/O");
    reducedTree_Triton -> Branch("HasPointOnITSLayer4",&HasPointOnITSLayer4,"HasPointOnITSLayer4/O");
    reducedTree_Triton -> Branch("HasPointOnITSLayer5",&HasPointOnITSLayer5,"HasPointOnITSLayer5/O");
    reducedTree_Triton -> Branch("chi2_TPC",&chi2_TPC,"chi2_TPC/D");
    reducedTree_Triton -> Branch("chi2_ITS",&chi2_ITS,"chi2_ITS/D");
    reducedTree_Triton -> Branch("ITSsignal",&ITSsignal,"ITSsignal/D");
    reducedTree_Triton -> Branch("TPCsignal",&TPCsignal,"TPCsignal/D");
    reducedTree_Triton -> Branch("TOFsignal",&TOFsignal,"TOFsignal/D");
    reducedTree_Triton -> Branch("nSigmaITS_Trit",&nSigmaITS_Trit,"nSigmaITS_Trit/D");
    reducedTree_Triton -> Branch("nSigmaTPC_Trit",&nSigmaTPC_Trit,"nSigmaTPC_Trit/D");
    reducedTree_Triton -> Branch("nSigmaTOF_Trit",&nSigmaTOF_Trit,"nSigmaTOF_Trit/D");



    PostData(1, fOutputList);
    PostData(2,reducedTree_Triton);

    //PostData(2, fQAList);
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskTritonVsMultiplicity_XeXe::UserExec(Option_t *)
{

    //Get Input Event
    if ( !GetInputEvent ()) return;

    //Load PID Response
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler = (AliInputEventHandler*) (mgr->GetInputEventHandler());
    fPIDResponse = inputHandler->GetPIDResponse();


    //Loop over Reconstructed Tracks
    for (Int_t i=0 ; i<fAODevent->GetNumberOfTracks() ; i++)  {

        //Track Selection
        AliAODTrack *track = (AliAODTrack*) fAODevent -> GetTrack(i);
        if ( !track ) continue;
        if ( PassedTrackQualityCutsNoDCA (track)) {
            if (IsCleanTritonCandidate(track))  {
                if (track->Charge()>0) histoDCAxyTriton_vs_pt     -> Fill (track->Pt(),GetDCAxy(track));
                if (track->Charge()<0) histoDCAxyAntiTriton_vs_pt -> Fill (track->Pt(),GetDCAxy(track));
            }
        }

        if ( PassedTrackQualityCuts (track)) {


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

        }

        if (IsTritonCandidate(track)){

            pt = track -> Pt();
            p  = track -> P();
            px = track -> Px();
            py = track -> Py();
            pz = track -> Pz();


            q  = (Int_t) track -> Charge();
            eta = track -> Eta();

            //Rapidity Calculation
            Double_t m  = AliPID::ParticleMass(AliPID::kTriton);
            Double_t fpx = track -> Px();
            Double_t fpy = track -> Py();
            Double_t fpz = track -> Pz();
            Double_t E = TMath::Sqrt(m*m + fpx*fpx + fpy*fpy + fpz*fpz);
            TLorentzVector P (fpx,fpy,fpz,E);
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
            HasPointOnITSLayer2 = track->HasPointOnITSLayer(2);
            HasPointOnITSLayer3 = track->HasPointOnITSLayer(3);
            HasPointOnITSLayer4 = track->HasPointOnITSLayer(4);
            HasPointOnITSLayer5 = track->HasPointOnITSLayer(5);

            chi2_TPC = track -> GetTPCchi2();//    check -> seems to be 0
            chi2_ITS = track -> GetITSchi2();//    check

            ITSsignal = track->GetITSsignal();
            TPCsignal = track->GetTPCsignal();
            TOFsignal = track->GetTOFsignal();

            nSigmaITS_Trit = fPIDResponse -> NumberOfSigmasITS(track,AliPID::kTriton);
            nSigmaTPC_Trit = fPIDResponse -> NumberOfSigmasTPC(track,AliPID::kTriton);
            nSigmaTOF_Trit = fPIDResponse -> NumberOfSigmasTOF(track,AliPID::kTriton);

            reducedTree_Triton -> Fill();

        }
    }

    PostData(1, fOutputList);
    PostData(2,reducedTree_Triton);
    //PostData(2, fQAList);

}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskTritonVsMultiplicity_XeXe::GetInputEvent ()  {

    //Get Input Event
    fAODevent = dynamic_cast <AliAODEvent*>(InputEvent());
    if (!fAODevent) return false;
    histoNumberOfEvents -> Fill(0.5);

    //Standard Event Cuts
    if (!fAODeventCuts.AcceptEvent(fAODevent)) {
        //PostData(2, fQAList);
        return false;
    }
    histoNumberOfEvents -> Fill(1.5);

    //Centrality
    AliMultSelection *multiplicitySelection = (AliMultSelection*) fAODevent->FindListObject("MultSelection");
    if( !multiplicitySelection) return false;
    histoNumberOfEvents -> Fill(2.5);
    Double_t centrality = multiplicitySelection->GetMultiplicityPercentile(fCentralityEstimator);


    multPercentile_V0A  = multiplicitySelection->GetMultiplicityPercentile("V0A");

    //Selection of Centrality Range
    if (centrality<fCentralityMin || centrality>=fCentralityMax ) return false;
    histoNumberOfEvents -> Fill(3.5);

    //Primary Vertex
    AliAODVertex *vertex = (AliAODVertex*) fAODevent->GetPrimaryVertex();
    if ( !vertex ) return false;
    histoNumberOfEvents -> Fill(4.5);

    //Primary Vertex Selection
    if ( vertex->GetZ() < fVertexZmin ) return false;
    if ( vertex->GetZ() > fVertexZmax ) return false;
    histoNumberOfEvents -> Fill(5.5);

    if ( vertex->GetNContributors() < fNumberVertexContributorsMin ) return false;
    histoNumberOfEvents -> Fill(6.5);
    multPercentile_V0M  = multiplicitySelection->GetMultiplicityPercentile("V0M");

    return true;

}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskTritonVsMultiplicity_XeXe::PassedTrackQualityCuts (AliAODTrack* track)  {

    //Filterbit
    if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return false;


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
    if ( track->GetTPCNCrossedRows() < fNumberCrossedRowsTPCMin ) return false;
    if ( static_cast<Double_t>(track->GetTPCNCrossedRows())/static_cast<Double_t>(track->GetTPCNclsF()) < fCrossedRowsFindableClsMin) return false;
    if ( track->GetTPCsignalN() < fNumberClustersTPCdEdxMin ) return false;
    if ( track->Chi2perNDF() > fChiSquarePerNDFMax) return false;

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
Bool_t AliAnalysisTaskTritonVsMultiplicity_XeXe::PassedTrackQualityCutsNoDCA (AliAODTrack* track)  {

    //Filterbit
    if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return false;


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
    if ( track->GetTPCNCrossedRows() < fNumberCrossedRowsTPCMin ) return false;
    if ( static_cast<Double_t>(track->GetTPCNCrossedRows())/static_cast<Double_t>(track->GetTPCNclsF()) < fCrossedRowsFindableClsMin) return false;
    if ( track->GetTPCsignalN() < fNumberClustersTPCdEdxMin ) return false;
    if ( track->Chi2perNDF() > fChiSquarePerNDFMax) return false;

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
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskTritonVsMultiplicity_XeXe::IsTritonCandidate (AliAODTrack *track)  {

  Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC(track,AliPID::kTriton);
  if (track->Pt() < 1.5) {
    if ( TMath::Abs(nsigmaTPC) > 7. ) return false;
  } else {
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF(track,AliPID::kTriton);
    if ( TMath::Abs(nsigmaTOF) > 7. ) return false;
    if ( TMath::Abs(nsigmaTPC) > 7. ) return false;
  }

  return true;
}


Bool_t AliAnalysisTaskTritonVsMultiplicity_XeXe::IsCleanTritonCandidate (AliAODTrack *track)  {

    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kTriton);
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kTriton);
    if (TMath::Abs(nsigmaTOF) > fnSigmaTOFmax) return false;
    if (TMath::Abs(nsigmaTPC) > fnSigmaTPCmax) return false;

    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskTritonVsMultiplicity_XeXe::Centered_nsigmaTPC (AliAODTrack *track)  {

    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kTriton);

    Double_t mean_fitted  = fpar0_mean_TPC*exp(fpar1_mean_TPC*(track->P()));
    Double_t sigma_fitted = fpar0_sigma_TPC*(track->P());

    nsigmaTPC = (nsigmaTPC - mean_fitted)/sigma_fitted;

    return nsigmaTPC;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskTritonVsMultiplicity_XeXe::Centered_nsigmaTOF (AliAODTrack *track)  {

    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kTriton);

    Double_t mean_fitted  = fpar0_mean_TOF*exp(fpar1_mean_TOF*(track->P()));
    Double_t sigma_fitted = fpar0_sigma_TOF*exp(fpar1_sigma_TOF*(track->P()));

    nsigmaTOF = (nsigmaTOF - mean_fitted)/sigma_fitted;

    return nsigmaTOF;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskTritonVsMultiplicity_XeXe::PassedTOFSelection (AliAODTrack *track)  {

    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kTriton);
    if (TMath::Abs(nsigmaTOF) > fnSigmaTOFmax) return false;

    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskTritonVsMultiplicity_XeXe::PassedTPCSelection (AliAODTrack *track)  {

    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kTriton);
    if (TMath::Abs(nsigmaTPC) > fnSigmaTPCmax || (track->GetStatus()&AliAODTrack::kTOFout)!=0) return false;

    //TPC-TOF Matching
//   if ((track->GetStatus()&AliAODTrack::kTOFout)==0) hasTOFhit=0;//Track with no TOF hit
//   if ((track->GetStatus()&AliAODTrack::kTOFout)!=0) hasTOFhit=1;//Track with TOF hit

    return true;
}

//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskTritonVsMultiplicity_XeXe::GetDCAxy (AliAODTrack *track)  {

    Double_t impactParameter[2];
    Double_t covarianceMatrix[3];
    if (!track->PropagateToDCA(fAODevent->GetPrimaryVertex(),fAODevent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;

    Double_t DCAxy = impactParameter[0];

    return DCAxy;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskTritonVsMultiplicity_XeXe::GetDCAz (AliAODTrack *track)  {

    Double_t impactParameter[2];
    Double_t covarianceMatrix[3];
    if (!track->PropagateToDCA(fAODevent->GetPrimaryVertex(),fAODevent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;

    Double_t DCAz = impactParameter[1];

    return DCAz;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskTritonVsMultiplicity_XeXe::Terminate(Option_t *)  {

    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________

