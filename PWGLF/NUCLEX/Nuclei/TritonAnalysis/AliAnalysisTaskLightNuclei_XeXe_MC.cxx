#include "AliAnalysisTaskLightNuclei_XeXe_MC.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTask.h"
#include "AliPIDResponse.h"
#include "TLorentzVector.h"
#include "AliAODMCHeader.h"
#include "AliCentrality.h"
#include "AliEventplane.h"
#include "AliEventCuts.h"
#include "TDatabasePDG.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliVVertex.h"
#include "AliMCEvent.h"
#include "TObjArray.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TVectorD.h"
#include "TRandom.h"
#include "TChain.h"
#include "AliLog.h"
#include "TMath.h"
#include "TList.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"

ClassImp(AliAnalysisTaskLightNuclei_XeXe_MC)

//_________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskLightNuclei_XeXe_MC::AliAnalysisTaskLightNuclei_XeXe_MC():
AliAnalysisTaskSE(),
fAODevent(nullptr),
fMCevent(nullptr),
fPIDResponse(nullptr),
fAODeventCuts(),
fUtils(nullptr),
fOutputList(nullptr),
fQAList(nullptr),
fAODMCHeader(nullptr),
fAODArrayMCParticles(nullptr),
fMCEventHandler(nullptr),
reducedTree_gen_He3(nullptr),
reducedTree_rec_He3(nullptr),
reducedTree_gen_H3(nullptr),
reducedTree_rec_H3(nullptr),
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
multPercentile_V0M(0),
particleType(0),
centrality(0),
xVertex(0),
yVertex(0),
zVertex(0),
ID_event(0),
Helium3Sec(0),
TritonSec(0),
lp(0),
px(0),
py(0),
pz(0),
pt(0),
pt_particle(0),
pz_particle(0),
deltapt(0),
deltapz(0),
charge(0),
dcaxy(0),
dcaz(0),
trackType(0),
nTPC_Clusters(0),
nITS_Clusters(0),
nTPC_Clusters_dEdx(0),
nTPC_FindableClusters(0),
nTPC_CrossedRows(0),
HasPointOnITSLayer0(0),
HasPointOnITSLayer1(0),
HasPointOnITSLayer2(0),
HasPointOnITSLayer3(0),
HasPointOnITSLayer4(0),
HasPointOnITSLayer5(0),
chi2_NDF(0),
chi2_ITS(0),
nSigmaITS(0),
nSigmaTPC(0),
nSigmaTOF(0),
PrimaryParticle(0),
SecondaryMaterial(0),
SecondaryDecay(0),
PrimaryTrack(0),
TPC_signal(0),
ITS_signal(0),
TOF_signal(0),
TRDntracklets(0),
hasTOFhit(0)
{}
//_________________________________________________________________________________________________________________________________________________________________________________________________
AliAnalysisTaskLightNuclei_XeXe_MC::AliAnalysisTaskLightNuclei_XeXe_MC(const char *name):
AliAnalysisTaskSE(name),
fAODevent(nullptr),
fMCevent(nullptr),
fPIDResponse(nullptr),
fAODeventCuts(),
fUtils(nullptr),
fOutputList(nullptr),
fQAList(nullptr),
fAODMCHeader(nullptr),
fAODArrayMCParticles(nullptr),
fMCEventHandler(nullptr),
reducedTree_gen_He3(nullptr),
reducedTree_rec_He3(nullptr),
reducedTree_gen_H3(nullptr),
reducedTree_rec_H3(nullptr),
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
multPercentile_V0M(0),
particleType(0),
centrality(0),
xVertex(0),
yVertex(0),
zVertex(0),
ID_event(0),
Helium3Sec(0),
TritonSec(0),
lp(0),
px(0),
py(0),
pz(0),
pt(0),
pt_particle(0),
pz_particle(0),
deltapt(0),
deltapz(0),
charge(0),
dcaxy(0),
dcaz(0),
trackType(0),
nTPC_Clusters(0),
nITS_Clusters(0),
nTPC_Clusters_dEdx(0),
nTPC_FindableClusters(0),
nTPC_CrossedRows(0),
HasPointOnITSLayer0(0),
HasPointOnITSLayer1(0),
HasPointOnITSLayer2(0),
HasPointOnITSLayer3(0),
HasPointOnITSLayer4(0),
HasPointOnITSLayer5(0),
chi2_NDF(0),
chi2_ITS(0),
nSigmaITS(0),
nSigmaTPC(0),
nSigmaTOF(0),
PrimaryParticle(0),
SecondaryMaterial(0),
SecondaryDecay(0),
PrimaryTrack(0),
TPC_signal(0),
ITS_signal(0),
TOF_signal(0),
TRDntracklets(0),
hasTOFhit(0)
{
    fUtils = new AliAnalysisUtils();
    DefineInput (0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());

}
//_________________________________________________________________________________________________________________________________________________________________________________________________
            AliAnalysisTaskLightNuclei_XeXe_MC::~AliAnalysisTaskLightNuclei_XeXe_MC(){
    fOutputList -> Clear();
    fQAList     -> Clear();
    delete fAODevent;
    delete fMCevent;
    delete fPIDResponse;
    delete fUtils;
    delete fOutputList;
    delete fQAList;
    delete fAODMCHeader;
    delete fAODArrayMCParticles;
    delete fMCEventHandler;

}
//_________________________________________________________________________________________________________________________________________________________________________________________________
void        AliAnalysisTaskLightNuclei_XeXe_MC::UserCreateOutputObjects(){
    fOutputList = new TList();
    fOutputList -> SetOwner();

    fQAList = new TList();
    fQAList -> SetOwner();


    //Number of Events
    histoNumberOfEvents                     = new TH1F("histoNumberOfEvents","Events after selection steps",10,0,10);
    fQAList                                 -> Add (histoNumberOfEvents);
    fAODeventCuts.AddQAplotsToList(fQAList);                //Add event selection QA plots
    //Signal Extraction
    histoNsigmaTPCHe3_vs_pt                 = new TH2F ("histoNsigmaTPCHe3_vs_pt","",500,0,5,1000,-20,20);
    histoNsigmaTOFHe3_vs_pt                 = new TH2F ("histoNsigmaTOFHe3_vs_pt","",500,0,5,1000,-20,20);
    histoNsigmaTPCantiHe3_vs_pt             = new TH2F ("histoNsigmaTPCantiHe3_vs_pt","",500,0,5,1000,-20,20);
    histoNsigmaTOFantiHe3_vs_pt             = new TH2F ("histoNsigmaTOFantiHe3_vs_pt","",500,0,5,1000,-20,20);
    histoNsigmaTPCHe3_vs_pt_centered        = new TH2F ("histoNsigmaTPCHe3_vs_pt_centered","",500,0,5,1000,-20,20);
    histoNsigmaTPCantiHe3_vs_pt_centered    = new TH2F ("histoNsigmaTPCantiHe3_vs_pt_centered","",500,0,5,1000,-20,20);
    histoNsigmaTOFHe3_vs_pt_centered        = new TH2F ("histoNsigmaTOFHe3_vs_pt_centered","",500,0,5,1000,-20,20);
    histoNsigmaTOFantiHe3_vs_pt_centered    = new TH2F ("histoNsigmaTOFantiHe3_vs_pt_centered","",500,0,5,1000,-20,20);
    histoNsigmaTOFHe3_vs_pt_trd             = new TH2F ("histoNsigmaTOFHe3_vs_pt_trd","",500,0,5,1000,-20,20);
    histoNsigmaTOFantiHe3_vs_pt_trd         = new TH2F ("histoNsigmaTOFantiHe3_vs_pt_trd","",500,0,5,1000,-20,20);
    histoNsigmaTPCHe3_vs_p                  = new TH2F ("histoNsigmaTPCHe3_vs_p","",500,0,5,1000,-20,20);
    histoNsigmaTPCantiHe3_vs_p              = new TH2F ("histoNsigmaTPCantiHe3_vs_p","",500,0,5,1000,-20,20);
    histoNsigmaTOFHe3_vs_p                  = new TH2F ("histoNsigmaTOFHe3_vs_p","",500,0,5,1000,-20,20);
    histoNsigmaTOFantiHe3_vs_p              = new TH2F ("histoNsigmaTOFantiHe3_vs_p","",500,0,5,1000,-20,20);
    histoNsigmaTPCHe3_vs_p_notof            = new TH2F ("histoNsigmaTPCHe3_vs_p_notof","",500,0,5,1000,-20,20);
    histoNsigmaTPCantiHe3_vs_p_notof        = new TH2F ("histoNsigmaTPCantiHe3_vs_p_notof","",500,0,5,1000,-20,20);
    histoDCAxyHe3_vs_pt                     = new TH2F ("histoDCAxyHe3_vs_pt","",500,0,5,500,-5,5);
    histoDCAxyAntiHe3_vs_pt                 = new TH2F ("histoDCAxyAntiHe3_vs_pt","",500,0,5,500,-5,5);


    histoNsigmaTPCHe3_vs_pt                 -> Sumw2();
    histoNsigmaTOFHe3_vs_pt                 -> Sumw2();
    histoNsigmaTPCantiHe3_vs_pt             -> Sumw2();
    histoNsigmaTOFantiHe3_vs_pt             -> Sumw2();
    histoNsigmaTPCHe3_vs_pt_centered        -> Sumw2();
    histoNsigmaTPCantiHe3_vs_pt_centered    -> Sumw2();
    histoNsigmaTOFHe3_vs_pt_centered        -> Sumw2();
    histoNsigmaTOFantiHe3_vs_pt_centered    -> Sumw2();
    histoNsigmaTOFHe3_vs_pt_trd             -> Sumw2();
    histoNsigmaTOFantiHe3_vs_pt_trd         -> Sumw2();
    histoNsigmaTPCHe3_vs_p                  -> Sumw2();
    histoNsigmaTPCantiHe3_vs_p              -> Sumw2();
    histoNsigmaTOFHe3_vs_p                  -> Sumw2();
    histoNsigmaTOFantiHe3_vs_p              -> Sumw2();
    histoNsigmaTPCHe3_vs_p_notof            -> Sumw2();
	histoNsigmaTPCantiHe3_vs_p_notof        -> Sumw2();
    histoDCAxyHe3_vs_pt                     -> Sumw2();
    histoDCAxyAntiHe3_vs_pt                 -> Sumw2();

    fQAList                                 -> Add(histoNsigmaTPCHe3_vs_pt);
    fQAList                                 -> Add(histoNsigmaTOFHe3_vs_pt);
    fQAList                                 -> Add(histoNsigmaTPCantiHe3_vs_pt);
    fQAList                                 -> Add(histoNsigmaTOFantiHe3_vs_pt);
    fQAList                                 -> Add(histoNsigmaTPCHe3_vs_pt_centered);
    fQAList                                 -> Add(histoNsigmaTPCantiHe3_vs_pt_centered);
    fQAList                                 -> Add(histoNsigmaTOFHe3_vs_pt_centered);
    fQAList                                 -> Add(histoNsigmaTOFantiHe3_vs_pt_centered);
    fQAList                                 -> Add(histoNsigmaTOFHe3_vs_pt_trd);
    fQAList                                 -> Add(histoNsigmaTOFantiHe3_vs_pt_trd);
    fQAList                                 -> Add(histoNsigmaTPCHe3_vs_p);
    fQAList                                 -> Add(histoNsigmaTPCantiHe3_vs_p);
    fQAList                                 -> Add(histoNsigmaTOFHe3_vs_p);
    fQAList                                 -> Add(histoNsigmaTOFantiHe3_vs_p);
    fQAList                                 -> Add(histoNsigmaTPCHe3_vs_p_notof);
    fQAList                                 -> Add(histoNsigmaTPCantiHe3_vs_p_notof);
    fQAList                                 -> Add (histoDCAxyHe3_vs_pt);
    fQAList                                 -> Add (histoDCAxyAntiHe3_vs_pt);


    //Reduced Tree Generated particles
    reducedTree_gen_He3 = new TTree("reducedTree_gen_He3","reducedTree_gen_He3");
    reducedTree_gen_He3 -> Branch("particleType",&particleType,"particleType/I");
    reducedTree_gen_He3 -> Branch("centrality",&centrality,"centrality/F");
    reducedTree_gen_He3 -> Branch("multPercentile_V0M",&multPercentile_V0M,"multPercentile_V0M/F");
    reducedTree_gen_He3 -> Branch("xVertex",&xVertex,"xVertex/F");
    reducedTree_gen_He3 -> Branch("yVertex",&yVertex,"yVertex/F");
    reducedTree_gen_He3 -> Branch("zVertex",&zVertex,"zVertex/F");
    reducedTree_gen_He3 -> Branch("px",&px,"px/F");
    reducedTree_gen_He3 -> Branch("py",&py,"py/F");
    reducedTree_gen_He3 -> Branch("pz",&pz,"pz/F");
    reducedTree_gen_He3 -> Branch("pt",&pt,"pt/F");
    reducedTree_gen_He3 -> Branch("charge",&charge,"charge/I");
    reducedTree_gen_He3 -> Branch("PrimaryParticle",&PrimaryParticle,"PrimaryParticle/I");

    reducedTree_gen_H3 = new TTree("reducedTree_gen_H3","reducedTree_gen_H3");
    reducedTree_gen_H3 -> Branch("particleType",&particleType,"particleType/I");
    reducedTree_gen_H3 -> Branch("centrality",&centrality,"centrality/F");
    reducedTree_gen_H3 -> Branch("multPercentile_V0M",&multPercentile_V0M,"multPercentile_V0M/F");
    reducedTree_gen_H3 -> Branch("xVertex",&xVertex,"xVertex/F");
    reducedTree_gen_H3 -> Branch("yVertex",&yVertex,"yVertex/F");
    reducedTree_gen_H3 -> Branch("zVertex",&zVertex,"zVertex/F");
    reducedTree_gen_H3 -> Branch("px",&px,"px/F");
    reducedTree_gen_H3 -> Branch("py",&py,"py/F");
    reducedTree_gen_H3 -> Branch("pz",&pz,"pz/F");
    reducedTree_gen_H3 -> Branch("pt",&pt,"pt/F");
    reducedTree_gen_H3 -> Branch("charge",&charge,"charge/I");
    reducedTree_gen_H3 -> Branch("PrimaryParticle",&PrimaryParticle,"PrimaryParticle/I");

    //Reduced Tree Reconstructed tracks
    reducedTree_rec_He3 = new TTree("reducedTree_rec_He3","reducedTree_rec_He3");
    reducedTree_rec_He3 -> Branch("trackType",&trackType,"trackType/I");
    reducedTree_rec_He3 -> Branch("centrality",&centrality,"centrality/F");
    reducedTree_rec_He3 -> Branch("multPercentile_V0M",&multPercentile_V0M,"multPercentile_V0M/F");
    reducedTree_rec_He3 -> Branch("xVertex",&xVertex,"xVertex/F");
    reducedTree_rec_He3 -> Branch("yVertex",&yVertex,"yVertex/F");
    reducedTree_rec_He3 -> Branch("zVertex",&zVertex,"zVertex/F");
    reducedTree_rec_He3 -> Branch("ID_event",&ID_event,"ID_event/I");
    reducedTree_rec_He3 -> Branch("Helium3Sec",&Helium3Sec,"Helium3Sec/I");
    reducedTree_rec_He3 -> Branch("lp",&lp,"lp/I");
    reducedTree_rec_He3 -> Branch("px",&px,"px/F");
    reducedTree_rec_He3 -> Branch("py",&py,"py/F");
    reducedTree_rec_He3 -> Branch("pz",&pz,"pz/F");
    reducedTree_rec_He3 -> Branch("pt",&pt,"pt/F");
    reducedTree_rec_He3 -> Branch("pt_particle",&pt_particle,"pt_particle/F");
    reducedTree_rec_He3 -> Branch("pz_particle",&pz_particle,"pz_particle/F");
    reducedTree_rec_He3 -> Branch("deltapt",&deltapt,"deltapt/F");
    reducedTree_rec_He3 -> Branch("deltapz",&deltapz,"deltapz/F");
    reducedTree_rec_He3 -> Branch("charge",&charge,"charge/I");
    reducedTree_rec_He3 -> Branch("dcaxy",&dcaxy,"dcaxy/D");
    reducedTree_rec_He3 -> Branch("dcaz",&dcaz,"dcaz/D");
    reducedTree_rec_He3 -> Branch("nTPC_Clusters",&nTPC_Clusters,"nTPC_Clusters/I");
    reducedTree_rec_He3 -> Branch("nITS_Clusters",&nITS_Clusters,"nITS_Clusters/I");
    reducedTree_rec_He3 -> Branch("nTPC_Clusters_dEdx",&nTPC_Clusters_dEdx,"nTPC_Clusters_dEdx/I");
    reducedTree_rec_He3 -> Branch("nTPC_FindableClusters",&nTPC_FindableClusters,"nTPC_FindableClusters/I");
    reducedTree_rec_He3 -> Branch("nTPC_CrossedRows",&nTPC_CrossedRows,"nTPC_CrossedRows/I");
    reducedTree_rec_He3 -> Branch("HasPointOnITSLayer0",&HasPointOnITSLayer0,"HasPointOnITSLayer0/I");
    reducedTree_rec_He3 -> Branch("HasPointOnITSLayer1",&HasPointOnITSLayer1,"HasPointOnITSLayer1/I");
    reducedTree_rec_He3 -> Branch("HasPointOnITSLayer2",&HasPointOnITSLayer2,"HasPointOnITSLayer2/I");
    reducedTree_rec_He3 -> Branch("HasPointOnITSLayer3",&HasPointOnITSLayer3,"HasPointOnITSLayer3/I");
    reducedTree_rec_He3 -> Branch("HasPointOnITSLayer4",&HasPointOnITSLayer4,"HasPointOnITSLayer4/I");
    reducedTree_rec_He3 -> Branch("HasPointOnITSLayer5",&HasPointOnITSLayer5,"HasPointOnITSLayer5/I");
    reducedTree_rec_He3 -> Branch("chi2_NDF",&chi2_NDF,"chi2_NDF/F");
    reducedTree_rec_He3 -> Branch("chi2_ITS",&chi2_ITS,"chi2_ITS/F");
    reducedTree_rec_He3 -> Branch("nSigmaITS",&nSigmaITS,"nSigmaITS/F");
    reducedTree_rec_He3 -> Branch("nSigmaTPC",&nSigmaTPC,"nSigmaTPC/F");
    reducedTree_rec_He3 -> Branch("nSigmaTOF",&nSigmaTOF,"nSigmaTOF/F");
    reducedTree_rec_He3 -> Branch("PrimaryParticle",&PrimaryParticle,"PrimaryParticle/I");
    reducedTree_rec_He3 -> Branch("SecondaryMaterial",&SecondaryMaterial,"SecondaryMaterial/I");
    reducedTree_rec_He3 -> Branch("SecondaryDecay",&SecondaryDecay,"SecondaryDecay/I");
    reducedTree_rec_He3 -> Branch("PrimaryTrack",&PrimaryTrack,"PrimaryTrack/I");
    reducedTree_rec_He3 -> Branch("TPC_signal",&TPC_signal,"TPC_signal/D");
    reducedTree_rec_He3 -> Branch("ITS_signal",&ITS_signal,"ITS_signal/D");
    reducedTree_rec_He3 -> Branch("TOF_signal",&TOF_signal,"TOF_signal/D");
    reducedTree_rec_He3 -> Branch("TRDntracklets",&TRDntracklets,"TRDntracklets/I");
    reducedTree_rec_He3 -> Branch("hasTOFhit",&hasTOFhit,"hasTOFhit/I");

    reducedTree_rec_H3 = new TTree("reducedTree_rec_H3","reducedTree_rec_H3");
    reducedTree_rec_H3 -> Branch("trackType",&trackType,"trackType/I");
    reducedTree_rec_H3 -> Branch("centrality",&centrality,"centrality/F");
    reducedTree_rec_H3 -> Branch("multPercentile_V0M",&multPercentile_V0M,"multPercentile_V0M/F");
    reducedTree_rec_H3 -> Branch("xVertex",&xVertex,"xVertex/F");
    reducedTree_rec_H3 -> Branch("yVertex",&yVertex,"yVertex/F");
    reducedTree_rec_H3 -> Branch("zVertex",&zVertex,"zVertex/F");
    reducedTree_rec_H3 -> Branch("ID_event",&ID_event,"ID_event/I");
    reducedTree_rec_H3 -> Branch("Helium3Sec",&Helium3Sec,"Helium3Sec/I");
    reducedTree_rec_H3 -> Branch("lp",&lp,"lp/I");
    reducedTree_rec_H3 -> Branch("px",&px,"px/F");
    reducedTree_rec_H3 -> Branch("py",&py,"py/F");
    reducedTree_rec_H3 -> Branch("pz",&pz,"pz/F");
    reducedTree_rec_H3 -> Branch("pt",&pt,"pt/F");
    reducedTree_rec_H3 -> Branch("pt_particle",&pt_particle,"pt_particle/F");
    reducedTree_rec_H3 -> Branch("pz_particle",&pz_particle,"pz_particle/F");
    reducedTree_rec_H3 -> Branch("deltapt",&deltapt,"deltapt/F");
    reducedTree_rec_H3 -> Branch("deltapz",&deltapz,"deltapz/F");
    reducedTree_rec_H3 -> Branch("charge",&charge,"charge/I");
    reducedTree_rec_H3 -> Branch("dcaxy",&dcaxy,"dcaxy/D");
    reducedTree_rec_H3 -> Branch("dcaz",&dcaz,"dcaz/D");
    reducedTree_rec_H3 -> Branch("nTPC_Clusters",&nTPC_Clusters,"nTPC_Clusters/I");
    reducedTree_rec_H3 -> Branch("nITS_Clusters",&nITS_Clusters,"nITS_Clusters/I");
    reducedTree_rec_H3 -> Branch("nTPC_Clusters_dEdx",&nTPC_Clusters_dEdx,"nTPC_Clusters_dEdx/I");
    reducedTree_rec_H3 -> Branch("nTPC_FindableClusters",&nTPC_FindableClusters,"nTPC_FindableClusters/I");
    reducedTree_rec_H3 -> Branch("nTPC_CrossedRows",&nTPC_CrossedRows,"nTPC_CrossedRows/I");
    reducedTree_rec_H3 -> Branch("HasPointOnITSLayer0",&HasPointOnITSLayer0,"HasPointOnITSLayer0/I");
    reducedTree_rec_H3 -> Branch("HasPointOnITSLayer1",&HasPointOnITSLayer1,"HasPointOnITSLayer1/I");
    reducedTree_rec_H3 -> Branch("HasPointOnITSLayer2",&HasPointOnITSLayer2,"HasPointOnITSLayer2/I");
    reducedTree_rec_H3 -> Branch("HasPointOnITSLayer3",&HasPointOnITSLayer3,"HasPointOnITSLayer3/I");
    reducedTree_rec_H3 -> Branch("HasPointOnITSLayer4",&HasPointOnITSLayer4,"HasPointOnITSLayer4/I");
    reducedTree_rec_H3 -> Branch("HasPointOnITSLayer5",&HasPointOnITSLayer5,"HasPointOnITSLayer5/I");
    reducedTree_rec_H3 -> Branch("chi2_NDF",&chi2_NDF,"chi2_NDF/F");
    reducedTree_rec_H3 -> Branch("chi2_ITS",&chi2_ITS,"chi2_ITS/F");
    reducedTree_rec_H3 -> Branch("nSigmaITS",&nSigmaITS,"nSigmaITS/F");
    reducedTree_rec_H3 -> Branch("nSigmaTPC",&nSigmaTPC,"nSigmaTPC/F");
    reducedTree_rec_H3 -> Branch("nSigmaTOF",&nSigmaTOF,"nSigmaTOF/F");
    reducedTree_rec_H3 -> Branch("PrimaryParticle",&PrimaryParticle,"PrimaryParticle/I");
    reducedTree_rec_H3 -> Branch("SecondaryMaterial",&SecondaryMaterial,"SecondaryMaterial/I");
    reducedTree_rec_H3 -> Branch("SecondaryDecay",&SecondaryDecay,"SecondaryDecay/I");
    reducedTree_rec_H3 -> Branch("PrimaryTrack",&PrimaryTrack,"PrimaryTrack/I");
    reducedTree_rec_H3 -> Branch("TPC_signal",&TPC_signal,"TPC_signal/D");
    reducedTree_rec_H3 -> Branch("ITS_signal",&ITS_signal,"ITS_signal/D");
    reducedTree_rec_H3 -> Branch("TOF_signal",&TOF_signal,"TOF_signal/D");
    reducedTree_rec_H3 -> Branch("TRDntracklets",&TRDntracklets,"TRDntracklets/I");
    reducedTree_rec_H3 -> Branch("hasTOFhit",&hasTOFhit,"hasTOFhit/I");


	fOutputList     -> Add(reducedTree_gen_He3);
	fOutputList     -> Add(reducedTree_rec_He3);
    fOutputList     -> Add(reducedTree_gen_H3);
	fOutputList     -> Add(reducedTree_rec_H3);


    PostData(1, fOutputList);
    PostData(2, fQAList);
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
void        AliAnalysisTaskLightNuclei_XeXe_MC::UserExec(Option_t *){
    //Get Input Event
    if ( !GetInputEvent ()) return;

    //Load PID Response
    fPIDResponse = fInputHandler->GetPIDResponse();
    if(!fPIDResponse) {
        AliError("No PID Response found");
        return;
    }

     if(!(fInputHandler->IsEventSelected() & AliVEvent::kINT7)) return;

	Int_t ID_Event=0;

        //Loop over Generated Particles
    for ( Int_t i=0; i<fAODArrayMCParticles->GetEntriesFast(); i++ ) {

        //MC Particle
        AliAODMCParticle *particle = (AliAODMCParticle*) fAODArrayMCParticles->At(i);
        if ( !particle) continue;
        if ( TMath::Abs(particle->GetPdgCode()) != 1000020030 && TMath::Abs(particle->GetPdgCode())!=1000010030 ) continue;

        if (TMath::Abs(particle->GetPdgCode()) == 1000020030) particleType = 1; // 1 for He3
        if (TMath::Abs(particle->GetPdgCode()) == 1000010030) particleType = 2; // 2 for Triton

        px = particle  -> Px();
        py = particle  -> Py();
        pz = particle  -> Pz();
        pt = particle  -> Pt();

        if (particle  -> Charge() > 0) charge = +1;
        if (particle  -> Charge() < 0) charge = -1;

        PrimaryParticle = 3;
        if (particle->IsPhysicalPrimary())        PrimaryParticle = 2;
        if (particle->IsSecondaryFromMaterial())  PrimaryParticle = 1;
        if (particle->IsSecondaryFromWeakDecay()) PrimaryParticle = 0;

        if      (TMath::Abs(particle->GetPdgCode()) == 1000020030) reducedTree_gen_He3 -> Fill();
        else if (TMath::Abs(particle->GetPdgCode()) == 1000010030) reducedTree_gen_H3  -> Fill();
    }

	Int_t 	Helium3Sec=0;
	Int_t 	TritonSec=0;

    //Loop over Reconstructed Tracks
    for (Int_t i=0 ; i<fAODevent->GetNumberOfTracks() ; i++)  {

        //Track Selection
        AliAODTrack *track = (AliAODTrack*) fAODevent -> GetTrack(i);
        if ( !track ) continue;
        if ( !PassedMinimalTrackQualityCuts (track)) continue;

        if ( PassedTrackQualityCutsNoDCA (track)) {
            if (IsCleanHe3Candidate(track))  {
                if (track->Charge()>0) histoDCAxyHe3_vs_pt     -> Fill (track->Pt(),GetDCAxy(track));
                if (track->Charge()<0) histoDCAxyAntiHe3_vs_pt -> Fill (track->Pt(),GetDCAxy(track));
            }
        }




        // Filling histograms for He3
        if (PassedTrackQualityCuts (track) && track->PdgCode() == 1000020030) {
            //Variables
            Double_t nsigma_TPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kHe3);
            Double_t nsigma_TOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kHe3);

            //TPC Signal vs. pT
            if (track->Charge()>0) histoNsigmaTPCHe3_vs_p_notof	            -> Fill (track->P(),nsigma_TPC);
            if (track->Charge()<0) histoNsigmaTPCantiHe3_vs_p_notof         -> Fill (track->P(),nsigma_TPC);

            if (PassedTOFSelection(track))  {
                if (track->Charge()>0) histoNsigmaTPCHe3_vs_pt              -> Fill (track->Pt(),nsigma_TPC);
                if (track->Charge()<0) histoNsigmaTPCantiHe3_vs_pt          -> Fill (track->Pt(),nsigma_TPC);
                if (track->Charge()>0) histoNsigmaTPCHe3_vs_p               -> Fill (track->P(), nsigma_TPC);
                if (track->Charge()<0) histoNsigmaTPCantiHe3_vs_p           -> Fill (track->P(), nsigma_TPC);
                if (track->Charge()>0) histoNsigmaTPCHe3_vs_pt_centered     -> Fill (track->Pt(),Centered_nsigmaTPC(track));
                if (track->Charge()<0) histoNsigmaTPCantiHe3_vs_pt_centered -> Fill (track->Pt(),Centered_nsigmaTPC(track));
            }
            //TOF Signal vs. pT
            if (PassedTPCSelection(track))  {
                if (track->Charge()>0) histoNsigmaTOFHe3_vs_pt              -> Fill (track->Pt(),nsigma_TOF);
                if (track->Charge()<0) histoNsigmaTOFantiHe3_vs_pt          -> Fill (track->Pt(),nsigma_TOF);
                if (track->Charge()>0) histoNsigmaTOFHe3_vs_p               -> Fill (track->P(), nsigma_TOF);
                if (track->Charge()<0) histoNsigmaTOFantiHe3_vs_p           -> Fill (track->P(), nsigma_TOF);
                if (track->Charge()>0 && track->GetTRDntrackletsPID()>fTRDntracklets) histoNsigmaTOFHe3_vs_pt_trd     -> Fill (track->Pt(),nsigma_TOF);
                if (track->Charge()<0 && track->GetTRDntrackletsPID()>fTRDntracklets) histoNsigmaTOFantiHe3_vs_pt_trd -> Fill (track->Pt(),nsigma_TOF);
                if (track->Charge()>0) histoNsigmaTOFHe3_vs_pt_centered     -> Fill (track->Pt(),Centered_nsigmaTOF(track));
                if (track->Charge()<0) histoNsigmaTOFantiHe3_vs_pt_centered -> Fill (track->Pt(),Centered_nsigmaTOF(track));
            }
        }




        ID_Event+=1;

        //Get Corresponding MC Particle
       Int_t lp = TMath::Abs(track->GetLabel());
       AliAODMCParticle *particle = (AliAODMCParticle*) fAODArrayMCParticles->At(lp);

        //He3 and Triton Selection
        if ( !particle)  {
            if (track->IsSecondaryFromMaterial() && track->PdgCode() == 1000020030 )   Helium3Sec ++;
            if (track->IsSecondaryFromMaterial() && track->PdgCode() == 1000010030 )   TritonSec ++;
            continue;
        }

        if ( TMath::Abs(particle->GetPdgCode()) != 1000020030  &&   TMath::Abs(particle->GetPdgCode())!=1000010030 ) continue;
	//========================================= PDG CODE NUCLEI =========================================//
        //  PDG Code Nuclei = 10LZZZAAAI
        //
        //  > L =  total number of strange quarks
        //  > ZZZ = number of protons given as 3-digit (ex. hydrogen = 001, helium-4 = 002, ...)
        //  > AAA = total baryon number given as 3-digit.
        //  > I = isomer level (with I = 0 corresponding to the ground state)
        //
        //===================================================================================================//

        trackType = (TMath::Abs(particle->GetPdgCode()) == 1000020030)? 1: 2;   // 1 for He3, 2 for Triton

        if      (trackType == 1) {if (!PassedCandidateSelectionHe3(track))    continue;}
        else if (trackType == 2) {if (!PassedCandidateSelectionTriton(track)) continue;}
        else                      continue;

        ID_event = ID_Event;

        px = track -> Px();
        py = track -> Py();
        pz = track -> Pz();
        pt = track -> Pt();

        pt_particle = particle->Pt();
        pz_particle = particle->Pz();

        deltapt = track->Pt() - particle->Pt();
        deltapz = track->Pz() - particle->Pz();

        if (track -> Charge() > 0) charge = +1;
        if (track -> Charge() < 0) charge = -1;

        //DCA
        dcaxy = GetDCAxy (track);
        dcaz  = GetDCAz  (track);

        nTPC_Clusters        = track -> GetTPCNcls();
        nITS_Clusters        = track -> GetITSNcls();
        nTPC_Clusters_dEdx   = track -> GetTPCsignalN();
        nTPC_FindableClusters= track -> GetTPCNclsF();
      	nTPC_CrossedRows     = track -> GetTPCNCrossedRows();

        if (track->HasPointOnITSLayer(0)) HasPointOnITSLayer0 = 1;
        if (track->HasPointOnITSLayer(1)) HasPointOnITSLayer1 = 1;
      	if (track->HasPointOnITSLayer(2)) HasPointOnITSLayer2 = 1;
        if (track->HasPointOnITSLayer(3)) HasPointOnITSLayer3 = 1;
        if (track->HasPointOnITSLayer(4)) HasPointOnITSLayer4 = 1;
        if (track->HasPointOnITSLayer(5)) HasPointOnITSLayer5 = 1;

        chi2_NDF            = track -> Chi2perNDF();
        chi2_ITS            = track -> GetITSchi2();

        //PID
        if (track->PdgCode() == 1000020030 ){
            nSigmaITS           = fPIDResponse -> NumberOfSigmasITS(track,AliPID::kHe3);
            nSigmaTPC           = fPIDResponse -> NumberOfSigmasTPC(track,AliPID::kHe3);
            nSigmaTOF           = fPIDResponse -> NumberOfSigmasTOF(track,AliPID::kHe3);
        }
        else if(track->PdgCode() == 1000010030){
            nSigmaITS           = fPIDResponse -> NumberOfSigmasITS(track,AliPID::kTriton);
            nSigmaTPC           = fPIDResponse -> NumberOfSigmasTPC(track,AliPID::kTriton);
            nSigmaTOF           = fPIDResponse -> NumberOfSigmasTOF(track,AliPID::kTriton);
        }

        //TPC-TOF Matching
        if ((track->GetStatus()&AliAODTrack::kTOFout)==0) hasTOFhit=0;//Track with no TOF hit
        if ((track->GetStatus()&AliAODTrack::kTOFout)!=0) hasTOFhit=1;//Track with TOF hit


        if (particle->IsSecondaryFromMaterial())  SecondaryMaterial = 1;
        else SecondaryMaterial=0;

        if (particle->IsSecondaryFromWeakDecay()) SecondaryDecay = 1;
        else SecondaryDecay=0;

        if (particle->IsPhysicalPrimary()) PrimaryParticle=1;
        else PrimaryParticle = 0;

        if(track->IsPrimaryCandidate()) PrimaryTrack=1;
        else PrimaryTrack=0;

        TPC_signal = track -> GetTPCsignal();
        ITS_signal = track -> GetITSsignal();
        TOF_signal = track -> GetTOFsignal();

        TRDntracklets = track -> GetTRDntrackletsPID();

        //Fill the Tree
        if      (TMath::Abs(particle->GetPdgCode()) == 1000020030) reducedTree_rec_He3 -> Fill();
        else if (TMath::Abs(particle->GetPdgCode()) == 1000010030) reducedTree_rec_H3  -> Fill();
        }

    PostData(1, fOutputList);
    PostData(2, fQAList);
    //PostData(3, reducedTree);
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t      AliAnalysisTaskLightNuclei_XeXe_MC::GetInputEvent ()  {

    //Get Input Event
    fAODevent = dynamic_cast <AliAODEvent*>(InputEvent());
    if (!fAODevent) return false;
    histoNumberOfEvents->Fill(0.5);

    //Get MC Event
    fMCEvent = MCEvent();
    if (!fMCEvent) return false;
    histoNumberOfEvents->Fill(1.0);

    //Event Cut
    if (!fAODeventCuts.AcceptEvent(fAODevent)) {
        PostData(2, fQAList);
        return false;
    }
    histoNumberOfEvents->Fill(1.5);

    //Vertex
    AliAODVertex *vertex = (AliAODVertex*) fAODevent->GetPrimaryVertex();
    if ( !vertex ) return false;

    //Multiplicity
    AliMultSelection *multiplicitySelection = (AliMultSelection*) fAODevent->FindListObject("MultSelection");
    if ( !multiplicitySelection) return false;
    histoNumberOfEvents->Fill(2.5);

    xVertex = vertex->GetX();
    yVertex = vertex->GetY();
    zVertex = vertex->GetZ();

    centrality           = multiplicitySelection->GetMultiplicityPercentile("V0A");
    multPercentile_V0M   = multiplicitySelection->GetMultiplicityPercentile("V0M");

    //identification
    fAODArrayMCParticles = dynamic_cast<TClonesArray *>(fAODevent->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!fAODArrayMCParticles) return false;

    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t      AliAnalysisTaskLightNuclei_XeXe_MC::PassedMinimalTrackQualityCuts (AliAODTrack *track)  {

    //Filterbit
    if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return false;

    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t      AliAnalysisTaskLightNuclei_XeXe_MC::PassedTrackQualityCuts (AliAODTrack* track)  {

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
Bool_t      AliAnalysisTaskLightNuclei_XeXe_MC::PassedCandidateSelectionHe3 (AliAODTrack *track)  {

    Double_t nsigmaHelium3TPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kHe3);
    if ( TMath::Abs(nsigmaHelium3TPC) > 6.0 ) return false;
    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t      AliAnalysisTaskLightNuclei_XeXe_MC::PassedCandidateSelectionTriton (AliAODTrack *track)  {

    Double_t nsigmaTritonTPC     = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kTriton);

    if (track->Pt() < 1.5) {
        if ( TMath::Abs(nsigmaTritonTPC) > 7. ) return false;
    } else {
        Double_t nsigmaTritonTOF = fPIDResponse -> NumberOfSigmasTOF(track,AliPID::kTriton);
        if ( TMath::Abs(nsigmaTritonTOF) > 7. ) return false;
        if ( TMath::Abs(nsigmaTritonTPC) > 7. ) return false;
    }


    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t      AliAnalysisTaskLightNuclei_XeXe_MC::PassedCandidateSelectionHelium3 (AliAODTrack *track)  {

    Double_t nsigmaHelium3TPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kHe3);
    //Double_t nsigmaHelium3TOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kHe3);
    if ( TMath::Abs(nsigmaHelium3TPC) < 6.0 ) return true;
    //if ( TMath::Abs(nsigmaHelium3TOF) > 6.0 ) return false;

    return false;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t    AliAnalysisTaskLightNuclei_XeXe_MC::GetDCAxy (AliAODTrack *track)  {

    Double_t impactParameter[2];
    Double_t covarianceMatrix[3];
    if (!track->PropagateToDCA(fAODevent->GetPrimaryVertex(),fAODevent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;

    Double_t DCAxy = impactParameter[0];

    return DCAxy;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t    AliAnalysisTaskLightNuclei_XeXe_MC::GetDCAz (AliAODTrack *track)  {

    Double_t impactParameter[2];
    Double_t covarianceMatrix[3];
    if (!track->PropagateToDCA(fAODevent->GetPrimaryVertex(),fAODevent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;

    Double_t DCAz = impactParameter[1];

    return DCAz;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t      AliAnalysisTaskLightNuclei_XeXe_MC::IsHe3Candidate (AliAODTrack *track)  {

    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC(track,AliPID::kHe3);
    if ( TMath::Abs(nsigmaTPC) > 6. ) return false;
    return true;
}

Bool_t      AliAnalysisTaskLightNuclei_XeXe_MC::IsCleanHe3Candidate (AliAODTrack *track)  {

    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kHe3);
    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kHe3);
    if (TMath::Abs(nsigmaTOF) > fnSigmaTOFmax) return false;
    if (TMath::Abs(nsigmaTPC) > fnSigmaTPCmax) return false;

    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t    AliAnalysisTaskLightNuclei_XeXe_MC::Centered_nsigmaTPC (AliAODTrack *track)  {

    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kHe3);

    Double_t mean_fitted  = fpar0_mean_TPC*exp(fpar1_mean_TPC*(track->P()));
    Double_t sigma_fitted = fpar0_sigma_TPC*(track->P());

    nsigmaTPC = (nsigmaTPC - mean_fitted)/sigma_fitted;

    return nsigmaTPC;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Double_t    AliAnalysisTaskLightNuclei_XeXe_MC::Centered_nsigmaTOF (AliAODTrack *track)  {

    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kHe3);

    Double_t mean_fitted  = fpar0_mean_TOF*exp(fpar1_mean_TOF*(track->P()));
    Double_t sigma_fitted = fpar0_sigma_TOF*exp(fpar1_sigma_TOF*(track->P()));

    nsigmaTOF = (nsigmaTOF - mean_fitted)/sigma_fitted;

    return nsigmaTOF;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t      AliAnalysisTaskLightNuclei_XeXe_MC::PassedTOFSelection (AliAODTrack *track)  {

    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF (track,AliPID::kHe3);
    if (TMath::Abs(nsigmaTOF) > fnSigmaTOFmax) return false;

    return true;
}
//_________________________________________________________________________________________________________________________________________________________________________________________________
Bool_t      AliAnalysisTaskLightNuclei_XeXe_MC::PassedTPCSelection (AliAODTrack *track)  {

    Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC (track,AliPID::kHe3);
    if (TMath::Abs(nsigmaTPC) > fnSigmaTPCmax) return false;

    //TPC-TOF Matching
    if ((track->GetStatus()&AliAODTrack::kTOFout)==0) hasTOFhit=0;//Track with no TOF hit
    if ((track->GetStatus()&AliAODTrack::kTOFout)!=0) hasTOFhit=1;//Track with TOF hit

    return true;
}

Bool_t      AliAnalysisTaskLightNuclei_XeXe_MC::PassedTrackQualityCutsNoDCA (AliAODTrack* track)  {

    //Filterbit
    if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return false;


    //Rapidity Calculation
    Double_t m  = AliPID::ParticleMass(AliPID::kHe3);
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
void        AliAnalysisTaskLightNuclei_XeXe_MC::Terminate(Option_t *)  {

    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) return;
}

