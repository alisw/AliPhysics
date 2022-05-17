#include "AliAnalysisTaskHe3VsMultiplicity_XeXe.h"
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


ClassImp(AliAnalysisTaskHe3VsMultiplicity_XeXe)

//_________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskHe3VsMultiplicity_XeXe::AliAnalysisTaskHe3VsMultiplicity_XeXe():
AliAnalysisTaskSE(),
fAODevent(nullptr),
fPIDResponse(nullptr),
fAODeventCuts(),
fUtils(nullptr),
fOutputList(nullptr),
reducedTree_He3(nullptr),
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
IsPrimaryCandidate(0),
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
HasSharedPointOnITSLayer0(0),
HasSharedPointOnITSLayer1(0),
chi2_TPC(0),
chi2_NDF(0),
chi2_ITS(0),
ITSsignal(0),
TPCsignal(0),
TOFsignal(0),
trackLength(0),
TPCmomentum(0),
nSigmaITS(0),
nSigmaTPC(0),
nSigmaTOF(0)
{}
//_________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskHe3VsMultiplicity_XeXe::AliAnalysisTaskHe3VsMultiplicity_XeXe(const char *name):
AliAnalysisTaskSE(name),
fAODevent(nullptr),
fPIDResponse(nullptr),
fAODeventCuts(),
fUtils(nullptr),
fOutputList(nullptr),
reducedTree_He3(nullptr),
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
multPercentile_V0A(0),
multPercentile_V0M(-1),
IsPrimaryCandidate(0),
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
HasSharedPointOnITSLayer0(0),
HasSharedPointOnITSLayer1(0),
chi2_TPC(0),
chi2_NDF(0),
chi2_ITS(0),
ITSsignal(0),
TPCsignal(0),
TOFsignal(0),
trackLength(0),
TPCmomentum(0),
nSigmaITS(0),
nSigmaTPC(0),
nSigmaTOF(0)
{
    fUtils = new AliAnalysisUtils();
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TTree::Class());
    //DefineOutput(2, TList::Class());
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskHe3VsMultiplicity_XeXe::~AliAnalysisTaskHe3VsMultiplicity_XeXe()
{
    fOutputList->Clear();
    delete fAODevent;
    delete fPIDResponse;
    delete fUtils;
    delete fOutputList;
    reducedTree_He3->GetUserInfo()->Clear();
    delete reducedTree_He3;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskHe3VsMultiplicity_XeXe::UserCreateOutputObjects()
{
    fOutputList = new TList();
    fOutputList -> SetOwner();

    //Number of Events
    histoNumberOfEvents = new TH1F("histoNumberOfEvents","Events after selection steps",10,0,10);
    fOutputList -> Add (histoNumberOfEvents);
    //Signal Extraction
    histoNsigmaTPCHe3_vs_pt     = new TH2F ("histoNsigmaTPCHe3_vs_pt","",500,0,5,1000,-20,20);
    histoNsigmaTOFHe3_vs_pt     = new TH2F ("histoNsigmaTOFHe3_vs_pt","",500,0,5,1000,-20,20);
    histoNsigmaTPCantiHe3_vs_pt = new TH2F ("histoNsigmaTPCantiHe3_vs_pt","",500,0,5,1000,-20,20);
    histoNsigmaTOFantiHe3_vs_pt = new TH2F ("histoNsigmaTOFantiHe3_vs_pt","",500,0,5,1000,-20,20);

    histoNsigmaTPCHe3_vs_pt_centered     = new TH2F ("histoNsigmaTPCHe3_vs_pt_centered","",500,0,5,1000,-20,20);
    histoNsigmaTPCantiHe3_vs_pt_centered = new TH2F ("histoNsigmaTPCantiHe3_vs_pt_centered","",500,0,5,1000,-20,20);

    histoNsigmaTOFHe3_vs_pt_centered     = new TH2F ("histoNsigmaTOFHe3_vs_pt_centered","",500,0,5,1000,-20,20);
    histoNsigmaTOFantiHe3_vs_pt_centered = new TH2F ("histoNsigmaTOFantiHe3_vs_pt_centered","",500,0,5,1000,-20,20);

    histoNsigmaTOFHe3_vs_pt_trd     = new TH2F ("histoNsigmaTOFHe3_vs_pt_trd","",500,0,5,1000,-20,20);
    histoNsigmaTOFantiHe3_vs_pt_trd = new TH2F ("histoNsigmaTOFantiHe3_vs_pt_trd","",500,0,5,1000,-20,20);

    histoNsigmaTPCHe3_vs_p         = new TH2F ("histoNsigmaTPCHe3_vs_p","",500,0,5,1000,-20,20);
    histoNsigmaTPCantiHe3_vs_p     = new TH2F ("histoNsigmaTPCantiHe3_vs_p","",500,0,5,1000,-20,20);
    histoNsigmaTOFHe3_vs_p         = new TH2F ("histoNsigmaTOFHe3_vs_p","",500,0,5,1000,-20,20);
    histoNsigmaTOFantiHe3_vs_p     = new TH2F ("histoNsigmaTOFantiHe3_vs_p","",500,0,5,1000,-20,20);

    histoNsigmaTPCHe3_vs_p_notof     = new TH2F ("histoNsigmaTPCHe3_vs_p_notof","",500,0,5,1000,-20,20);
    histoNsigmaTPCantiHe3_vs_p_notof = new TH2F ("histoNsigmaTPCantiHe3_vs_p_notof","",500,0,5,1000,-20,20);

    histoNsigmaTPCHe3_vs_pt     -> Sumw2();
    histoNsigmaTOFHe3_vs_pt     -> Sumw2();
    histoNsigmaTPCantiHe3_vs_pt -> Sumw2();
    histoNsigmaTOFantiHe3_vs_pt -> Sumw2();

    histoNsigmaTPCHe3_vs_pt_centered           -> Sumw2();
    histoNsigmaTPCantiHe3_vs_pt_centered       -> Sumw2();

    histoNsigmaTOFHe3_vs_pt_centered           -> Sumw2();
    histoNsigmaTOFantiHe3_vs_pt_centered       -> Sumw2();

    histoNsigmaTOFHe3_vs_pt_trd                -> Sumw2();
    histoNsigmaTOFantiHe3_vs_pt_trd            -> Sumw2();

    histoNsigmaTPCHe3_vs_p                     -> Sumw2();
    histoNsigmaTPCantiHe3_vs_p                 -> Sumw2();
    histoNsigmaTOFHe3_vs_p                     -> Sumw2();
    histoNsigmaTOFantiHe3_vs_p                 -> Sumw2();

    histoNsigmaTPCHe3_vs_p_notof              -> Sumw2();
	histoNsigmaTPCantiHe3_vs_p_notof          -> Sumw2();

    fOutputList -> Add(histoNsigmaTPCHe3_vs_pt);
    fOutputList -> Add(histoNsigmaTOFHe3_vs_pt);
    fOutputList -> Add(histoNsigmaTPCantiHe3_vs_pt);
    fOutputList -> Add(histoNsigmaTOFantiHe3_vs_pt);

    fOutputList -> Add(histoNsigmaTPCHe3_vs_pt_centered);
    fOutputList -> Add(histoNsigmaTPCantiHe3_vs_pt_centered);
    fOutputList -> Add(histoNsigmaTOFHe3_vs_pt_centered);
    fOutputList -> Add(histoNsigmaTOFantiHe3_vs_pt_centered);

    fOutputList -> Add(histoNsigmaTOFHe3_vs_pt_trd);
    fOutputList -> Add(histoNsigmaTOFantiHe3_vs_pt_trd);

    fOutputList -> Add(histoNsigmaTPCHe3_vs_p);
    fOutputList -> Add(histoNsigmaTPCantiHe3_vs_p);
    fOutputList -> Add(histoNsigmaTOFHe3_vs_p);
    fOutputList -> Add(histoNsigmaTOFantiHe3_vs_p);

    fOutputList -> Add(histoNsigmaTPCHe3_vs_p_notof);
    fOutputList -> Add(histoNsigmaTPCantiHe3_vs_p_notof);




    //DCA Distributions
    histoDCAxyHe3_vs_pt     = new TH2F ("histoDCAxyHe3_vs_pt","",500,0,5,500,-5,5);
    histoDCAxyAntiHe3_vs_pt = new TH2F ("histoDCAxyAntiHe3_vs_pt","",500,0,5,500,-5,5);

    histoDCAxyHe3_vs_pt     -> Sumw2();
    histoDCAxyAntiHe3_vs_pt -> Sumw2();

    fOutputList -> Add (histoDCAxyHe3_vs_pt);
    fOutputList -> Add (histoDCAxyAntiHe3_vs_pt);


    //Reduced Tree (He3)
    reducedTree_He3 = new TTree("reducedTree_He3","reducedTree_He3");
    reducedTree_He3 -> Branch("multPercentile_V0M",&multPercentile_V0M,"multPercentile_V0M/D");
    reducedTree_He3 -> Branch("multPercentile_V0A",&multPercentile_V0A,"multPercentile_V0A/D");
    reducedTree_He3 -> Branch("IsPrimaryCandidate",&IsPrimaryCandidate,"IsPrimaryCandidate/O");
    reducedTree_He3 -> Branch("pt",&pt,"pt/D");
    reducedTree_He3 -> Branch("p",&p,"p/D");
    reducedTree_He3 -> Branch("px",&px,"px/D");
    reducedTree_He3 -> Branch("py",&py,"py/D");
    reducedTree_He3 -> Branch("pz",&pz,"pz/D");
    reducedTree_He3 -> Branch("eta",&eta,"eta/D");
    reducedTree_He3 -> Branch("y",&y,"y/D");
    reducedTree_He3 -> Branch("q",&q,"q/I");
    reducedTree_He3 -> Branch("dcaxy",&dcaxy,"dcaxy/D");
    reducedTree_He3 -> Branch("dcaz",&dcaz,"dcaz/D");
    reducedTree_He3 -> Branch("nTPC_Clusters",&nTPC_Clusters,"nTPC_Clusters/I");
    reducedTree_He3 -> Branch("nITS_Clusters",&nITS_Clusters,"nITS_Clusters/I");
    reducedTree_He3 -> Branch("nTPC_FindableClusters",&nTPC_FindableClusters,"nTPC_FindableClusters/I");
    reducedTree_He3 -> Branch("nTPC_CrossedRows",&nTPC_CrossedRows,"nTPC_CrossedRows/I");
    reducedTree_He3 -> Branch("nTPC_Clusters_dEdx",&nTPC_Clusters_dEdx,"nTPC_Clusters_dEdx/I");
    reducedTree_He3 -> Branch("HasPointOnITSLayer0",&HasPointOnITSLayer0,"HasPointOnITSLayer0/O");
    reducedTree_He3 -> Branch("HasPointOnITSLayer1",&HasPointOnITSLayer1,"HasPointOnITSLayer1/O");
    reducedTree_He3 -> Branch("HasPointOnITSLayer2",&HasPointOnITSLayer2,"HasPointOnITSLayer2/O");
    reducedTree_He3 -> Branch("HasPointOnITSLayer3",&HasPointOnITSLayer3,"HasPointOnITSLayer3/O");
    reducedTree_He3 -> Branch("HasPointOnITSLayer4",&HasPointOnITSLayer4,"HasPointOnITSLayer4/O");
    reducedTree_He3 -> Branch("HasPointOnITSLayer5",&HasPointOnITSLayer5,"HasPointOnITSLayer5/O");
    reducedTree_He3 -> Branch("HasSharedPointOnITSLayer0",&HasSharedPointOnITSLayer0,"HasSharedPointOnITSLayer0/O");
    reducedTree_He3 -> Branch("HasSharedPointOnITSLayer1",&HasSharedPointOnITSLayer1,"HasSharedPointOnITSLayer1/O");
    reducedTree_He3 -> Branch("chi2_TPC",&chi2_TPC,"chi2_TPC/D");
    reducedTree_He3 -> Branch("chi2_ITS",&chi2_ITS,"chi2_ITS/D");
    reducedTree_He3 -> Branch("ITSsignal",&ITSsignal,"ITSsignal/D");
    reducedTree_He3 -> Branch("TPCsignal",&TPCsignal,"TPCsignal/D");
    reducedTree_He3 -> Branch("TOFsignal",&TOFsignal,"TOFsignal/D");
    reducedTree_He3 -> Branch("trackLength",&trackLength,"trackLength/D");
    reducedTree_He3 -> Branch("TPCmomentum",&TPCmomentum,"TPCmomentum/D");
    reducedTree_He3 -> Branch("nSigmaITS",&nSigmaITS,"nSigmaITS/D");
    reducedTree_He3 -> Branch("nSigmaTPC",&nSigmaTPC,"nSigmaTPC/D");
    reducedTree_He3 -> Branch("nSigmaTOF",&nSigmaTOF,"nSigmaTOF/D");



    PostData(1, fOutputList);
    PostData(2,reducedTree_He3);

}
//_________________________________________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskHe3VsMultiplicity_XeXe::UserExec(Option_t *)
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

            if (IsCleanHe3Candidate(track))  {
                if (track->Charge()>0) histoDCAxyHe3_vs_pt     -> Fill (track->Pt(),GetDCAxy(track));
                if (track->Charge()<0) histoDCAxyAntiHe3_vs_pt -> Fill (track->Pt(),GetDCAxy(track));
            }

        }

        if ( PassedTrackQualityCuts (track)) {


            //Variables
            Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kHe3);
            Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kHe3);


            //TPC Signal vs. pT


            if (track->Charge()>0) histoNsigmaTPCHe3_vs_p_notof	    -> Fill (track->P(),nsigmaTPC);
            if (track->Charge()<0) histoNsigmaTPCantiHe3_vs_p_notof -> Fill (track->P(),nsigmaTPC);


            if (PassedTOFSelection(track))  {

                if (track->Charge()>0) histoNsigmaTPCHe3_vs_pt     -> Fill (track->Pt(),nsigmaTPC);
                if (track->Charge()<0) histoNsigmaTPCantiHe3_vs_pt -> Fill (track->Pt(),nsigmaTPC);
                if (track->Charge()>0) histoNsigmaTPCHe3_vs_p      -> Fill (track->P(), nsigmaTPC);
                if (track->Charge()<0) histoNsigmaTPCantiHe3_vs_p  -> Fill (track->P(), nsigmaTPC);

                if (track->Charge()>0) histoNsigmaTPCHe3_vs_pt_centered     -> Fill (track->Pt(),Centered_nsigmaTPC(track));
                if (track->Charge()<0) histoNsigmaTPCantiHe3_vs_pt_centered -> Fill (track->Pt(),Centered_nsigmaTPC(track));


            }

            //TOF Signal vs. pT
            if (PassedTPCSelection(track))  {

                if (track->Charge()>0) histoNsigmaTOFHe3_vs_pt     -> Fill (track->Pt(),nsigmaTOF);
                if (track->Charge()<0) histoNsigmaTOFantiHe3_vs_pt -> Fill (track->Pt(),nsigmaTOF);
                if (track->Charge()>0) histoNsigmaTOFHe3_vs_p      -> Fill (track->P(), nsigmaTOF);
                if (track->Charge()<0) histoNsigmaTOFantiHe3_vs_p  -> Fill (track->P(), nsigmaTOF);


                if (track->Charge()>0 && track->GetTRDntrackletsPID()>fTRDntracklets) histoNsigmaTOFHe3_vs_pt_trd     -> Fill (track->Pt(),nsigmaTOF);
                if (track->Charge()<0 && track->GetTRDntrackletsPID()>fTRDntracklets) histoNsigmaTOFantiHe3_vs_pt_trd -> Fill (track->Pt(),nsigmaTOF);

                if (track->Charge()>0) histoNsigmaTOFHe3_vs_pt_centered     -> Fill (track->Pt(),Centered_nsigmaTOF(track));
                if (track->Charge()<0) histoNsigmaTOFantiHe3_vs_pt_centered -> Fill (track->Pt(),Centered_nsigmaTOF(track));



            }
        }
        if (IsHe3Candidate(track)){

            pt = track -> Pt();
            p  = track -> P();
            px = track -> Px();
            py = track -> Py();
            pz = track -> Pz();


            q  = (Int_t) track -> Charge();
            eta = track -> Eta();

            //Rapidity Calculation
            Double_t m  = AliPID::ParticleMass(AliPID::kHe3);
            Double_t fpx = track -> Px();
            Double_t fpy = track -> Py();
            Double_t fpz = track -> Pz();
            Double_t E = TMath::Sqrt(m*m + fpx*fpx + fpy*fpy + fpz*fpz);
            TLorentzVector P (fpx,fpy,fpz,E);
            y = P.Rapidity();


            IsPrimaryCandidate          = track->IsPrimaryCandidate();
            //DCA
            dcaz                        = GetDCAz (track);
            dcaxy                       = GetDCAxy(track);

            nTPC_Clusters               = track->GetTPCNcls();
            //nTRD_Clusters = track->GetTRDncls();
            nITS_Clusters               = track->GetITSNcls();

            nTPC_FindableClusters       = track->GetTPCNclsF();
            nTPC_CrossedRows            = track->GetTPCCrossedRows();
            nTPC_Clusters_dEdx          = track->GetTPCsignalN();

            HasPointOnITSLayer0         = track->HasPointOnITSLayer(0);
            HasPointOnITSLayer1         = track->HasPointOnITSLayer(1);
            HasPointOnITSLayer2         = track->HasPointOnITSLayer(2);
            HasPointOnITSLayer3         = track->HasPointOnITSLayer(3);
            HasPointOnITSLayer4         = track->HasPointOnITSLayer(4);
            HasPointOnITSLayer5         = track->HasPointOnITSLayer(5);

            HasSharedPointOnITSLayer0   = track->HasSharedPointOnITSLayer(0);
            HasSharedPointOnITSLayer1   = track->HasSharedPointOnITSLayer(1);

            chi2_TPC                    = track->GetTPCchi2();//    check
            chi2_ITS                    = track->GetITSchi2();//    check

            ITSsignal                   = track->GetITSsignal();
            TPCsignal                   = track->GetTPCsignal();
            TOFsignal                   = track->GetTOFsignal();
            //TRDsignal                 = track->GetTRDsignal();
            TPCmomentum                 = track->GetTPCmomentum();
            trackLength                 = track->GetIntegratedLength();
            nSigmaITS                   = fPIDResponse->NumberOfSigmasITS(track,AliPID::kHe3);
            nSigmaTPC                   = fPIDResponse->NumberOfSigmasTPC(track,AliPID::kHe3);
            nSigmaTOF                   = fPIDResponse->NumberOfSigmasTOF(track,AliPID::kHe3);
            //nSigmaTRD_Trit            = fPIDResponse -> NumberOfSigmasTRD(track,AliPID::kHe3);

            reducedTree_He3 -> Fill();
        }

    }


    PostData(1, fOutputList);
    PostData(2,reducedTree_He3);
    //PostData(2, fQAList);
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskHe3VsMultiplicity_XeXe::GetInputEvent ()  {

    //Get Input Event
    fAODevent = dynamic_cast <AliAODEvent*>(InputEvent());
    if (!fAODevent) return false;
    histoNumberOfEvents -> Fill(0.5);

    //Standard Event Cuts
    if (!fAODeventCuts.AcceptEvent(fAODevent)) {
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
Bool_t AliAnalysisTaskHe3VsMultiplicity_XeXe::PassedTrackQualityCuts (AliAODTrack* track)  {

    //Filterbit
    if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return false;


    //Rapidity Calculation
    Double_t m  = AliPID::ParticleMass(AliPID::kHe3);
    Double_t px = track -> Px();
    Double_t py = track -> Py();
    Double_t pz = track -> Pz();
    Double_t E  = TMath::Sqrt(m*m + px*px + py*py + pz*pz);
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
Bool_t AliAnalysisTaskHe3VsMultiplicity_XeXe::PassedTrackQualityCutsNoDCA (AliAODTrack* track)  {

    //Filterbit
    if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return false;


    //Rapidity Calculation
    Double_t m  = AliPID::ParticleMass(AliPID::kHe3);
    Double_t fpx = track -> Px();
    Double_t fpy = track -> Py();
    Double_t fpz = track -> Pz();
    Double_t E = TMath::Sqrt(m*m + fpx*fpx + fpy*fpy + fpz*fpz);
    TLorentzVector P (fpx,fpy,fpz,E);
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
Bool_t AliAnalysisTaskHe3VsMultiplicity_XeXe::IsHe3Candidate (AliAODTrack *track)  {

    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC(track,AliPID::kHe3);
    if ( TMath::Abs(nsigmaTPC) > 6. ) return false;
    return true;
}


Bool_t AliAnalysisTaskHe3VsMultiplicity_XeXe::IsCleanHe3Candidate (AliAODTrack *track)  {

    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kHe3);
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kHe3);
    if (TMath::Abs(nsigmaTOF) > fnSigmaTOFmax) return false;
    if (TMath::Abs(nsigmaTPC) > fnSigmaTPCmax) return false;

    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskHe3VsMultiplicity_XeXe::Centered_nsigmaTPC (AliAODTrack *track)  {

    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kHe3);

    Double_t mean_fitted  = fpar0_mean_TPC*exp(fpar1_mean_TPC*(track->P()));
    Double_t sigma_fitted = fpar0_sigma_TPC*(track->P());

    nsigmaTPC = (nsigmaTPC - mean_fitted)/sigma_fitted;

    return nsigmaTPC;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskHe3VsMultiplicity_XeXe::Centered_nsigmaTOF (AliAODTrack *track)  {

    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kHe3);

    Double_t mean_fitted  = fpar0_mean_TOF*exp(fpar1_mean_TOF*(track->P()));
    Double_t sigma_fitted = fpar0_sigma_TOF*exp(fpar1_sigma_TOF*(track->P()));

    nsigmaTOF = (nsigmaTOF - mean_fitted)/sigma_fitted;

    return nsigmaTOF;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskHe3VsMultiplicity_XeXe::PassedTOFSelection (AliAODTrack *track)  {

    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kHe3);
    if (TMath::Abs(nsigmaTOF) > fnSigmaTOFmax) return false;

    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskHe3VsMultiplicity_XeXe::PassedTPCSelection (AliAODTrack *track)  {

    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kHe3);
    if (TMath::Abs(nsigmaTPC) > fnSigmaTPCmax) return false;

    //TPC-TOF Matching
    //if ((track->GetStatus()&AliAODTrack::kTOFout)==0) hasTOFhit=0;//Track with no TOF hit
    //if ((track->GetStatus()&AliAODTrack::kTOFout)!=0) hasTOFhit=1;//Track with TOF hit

    return true;
}

//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskHe3VsMultiplicity_XeXe::GetDCAxy (AliAODTrack *track)  {

    Double_t impactParameter[2];
    Double_t covarianceMatrix[3];
    if (!track->PropagateToDCA(fAODevent->GetPrimaryVertex(),fAODevent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;

    Double_t DCAxy = impactParameter[0];

    return DCAxy;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskHe3VsMultiplicity_XeXe::GetDCAz (AliAODTrack *track)  {

    Double_t impactParameter[2];
    Double_t covarianceMatrix[3];
    if (!track->PropagateToDCA(fAODevent->GetPrimaryVertex(),fAODevent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;

    Double_t DCAz = impactParameter[1];

    return DCAz;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskHe3VsMultiplicity_XeXe::Terminate(Option_t *)  {

    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________

