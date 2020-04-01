#include "AliAnalysisTaskReducedTreeNuclei.h"
#include "AliInputEventHandler.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisUtils.h"
#include "AliEventCuts.h"
#include "AliMultSelection.h"
#include "AliMultEstimator.h"
#include "AliAnalysisTask.h"
#include "AliPIDResponse.h"
#include "AliCentrality.h"
#include "TDatabasePDG.h"
#include "AliAODVertex.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "TObjArray.h"
#include "TVector2.h"
#include "TVector3.h"
#include "AliAODv0.h"
#include "TRandom.h"
#include "TChain.h"
#include "TMath.h"
#include "TTree.h"
#include "TH1F.h"


ClassImp(AliAnalysisTaskReducedTreeNuclei)

//_______________________________________________________________________________________________________________________________________________________
AliAnalysisTaskReducedTreeNuclei::AliAnalysisTaskReducedTreeNuclei():
AliAnalysisTaskSE(),
fAODevent(NULL),
fAODVZERO(NULL),
fPIDResponse(NULL),
fAODeventCuts(),
fUtils(NULL),
fFillTri(kTRUE),
fFillHypTri(kTRUE),
fQAList(NULL),
// TreeEventSelection(NULL),
fOutputList(NULL),
histoEventSelection(NULL),
histoEventMultiplicity(NULL),
histoEventMultV0M_MB(NULL),
histoEventMultV0M_HM(NULL),
reducedTree_Helium(NULL),
reducedTree_Triton(NULL),
reducedTree_HyperTriton(NULL),
isTrigINT7(false),
isTrigHighMult(false),
// triggerMask(1),
SelectionStep(-1),
magFieldSign(0),
multPercentile_V0M(-1),
multPercentile_V0A(-1),
multPercentile_V0C(-1),
multPercentile_OnlineV0M(-1),
multPercentile_OnlineV0A(-1),
multPercentile_OnlineV0C(-1),
multPercentile_ADM(-1),
multPercentile_ADA(-1),
multPercentile_ADC(-1),
multPercentile_SPDClusters(-1),
multPercentile_SPDTracklets(-1),
multPercentile_RefMult05(-1),
multPercentile_RefMult08(-1),
multPercentile_CL1(-1),
multPercentile_ZNA(-1),
Ntrk_OnlineV0M(0),
Ntrk_OnlineV0A(0),
Ntrk_OnlineV0C(0),
Ntrk_V0M(0),
Ntrk_V0A(0),
Ntrk_V0C(0),
Ntrk_ADM(0),
Ntrk_ADA(0),
Ntrk_ADC(0),
Ntrk_SPDClusters(0),
Ntrk_SPDTracklets(0),
Ntrk_RefMult05(0),
Ntrk_RefMult08(0),
nVertexContributors(0),
xVertex(0),
yVertex(0),
zVertex(0),
px(0),
py(0),
pz(0),
TPCmomentum(0),
TRDmomentum(0),
integratedLength(0),
timeOfFlight(0),
beta(0),
gamma(0),
mass(0),
trackID(0),
eta(0),
phi(0),
theta(0),
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
HasSharedPointOnITSLayer0(0),
HasSharedPointOnITSLayer1(0),
HasSharedPointOnITSLayer2(0),
HasSharedPointOnITSLayer3(0),
HasSharedPointOnITSLayer4(0),
HasSharedPointOnITSLayer5(0),
chi2_TPC(0),
chi2_NDF(0),
chi2_ITS(0),
ITSsignal(0),
TPCsignal(0),
TOFsignal(0),
TRDsignal(0),
HMPIDsignal(0),
nSigmaITS_He4(0),
nSigmaTPC_He4(0),
nSigmaTOF_He4(0),
nSigmaTRD_He4(0),
nSigmaHMPID_He4(0),
nSigmaITS_He3(0),
nSigmaTPC_He3(0),
nSigmaTOF_He3(0),
nSigmaTRD_He3(0),
nSigmaHMPID_He3(0),
nSigmaITS_Trit(0),
nSigmaTPC_Trit(0),
nSigmaTOF_Trit(0),
nSigmaTRD_Trit(0),
nSigmaHMPID_Trit(0),
nSigmaITS_Deut(0),
nSigmaTPC_Deut(0),
nSigmaTOF_Deut(0),
nSigmaTRD_Deut(0),
nSigmaHMPID_Deut(0),
nSigmaITS_Prot(0),
nSigmaTPC_Prot(0),
nSigmaTOF_Prot(0),
nSigmaTRD_Prot(0),
nSigmaHMPID_Prot(0),
nSigmaITS_Pion(0),
nSigmaTPC_Pion(0),
nSigmaTOF_Pion(0),
nSigmaTRD_Pion(0),
nSigmaHMPID_Pion(0),
nSigmaITS_Kaon(0),
nSigmaTPC_Kaon(0),
nSigmaTOF_Kaon(0),
nSigmaTRD_Kaon(0),
nSigmaHMPID_Kaon(0),
nSigmaITS_Elec(0),
nSigmaTPC_Elec(0),
nSigmaTOF_Elec(0),
nSigmaTRD_Elec(0),
nSigmaHMPID_Elec(0),
px_Daughter1(0),
py_Daughter1(0),
pz_Daughter1(0),
TPCmomentum_Daughter1(0),
trackID_Daughter1(0),
eta_Daughter1(0),
phi_Daughter1(0),
theta_Daughter1(0),
y_Daughter1(0),
q_Daughter1(0),
dcaxy_Daughter1(0),
dcaz_Daughter1(0),
nTPC_Clusters_Daughter1(0),
nTPC_FindableClusters_Daughter1(0),
nTPC_CrossedRows_Daughter1(0),
nTPC_Clusters_dEdx_Daughter1(0),
nITS_Clusters_Daughter1(0),
HasPointOnITSLayer0_Daughter1(0),
HasPointOnITSLayer1_Daughter1(0),
HasPointOnITSLayer2_Daughter1(0),
HasPointOnITSLayer3_Daughter1(0),
HasPointOnITSLayer4_Daughter1(0),
HasPointOnITSLayer5_Daughter1(0),
HasSharedPointOnITSLayer0_Daughter1(0),
HasSharedPointOnITSLayer1_Daughter1(0),
HasSharedPointOnITSLayer2_Daughter1(0),
HasSharedPointOnITSLayer3_Daughter1(0),
HasSharedPointOnITSLayer4_Daughter1(0),
HasSharedPointOnITSLayer5_Daughter1(0),
chi2_TPC_Daughter1(0),
chi2_NDF_Daughter1(0),
chi2_ITS_Daughter1(0),
nSigmaITS_He3_Daughter1(0),
nSigmaTPC_He3_Daughter1(0),
nSigmaTOF_He3_Daughter1(0),
nSigmaITS_Pion_Daughter1(0),
nSigmaTPC_Pion_Daughter1(0),
nSigmaTOF_Pion_Daughter1(0),
nSigmaITS_Trit_Daughter1(0),
nSigmaTPC_Trit_Daughter1(0),
nSigmaTOF_Trit_Daughter1(0),
px_Daughter2(0),
py_Daughter2(0),
pz_Daughter2(0),
TPCmomentum_Daughter2(0),
trackID_Daughter2(0),
eta_Daughter2(0),
phi_Daughter2(0),
theta_Daughter2(0),
y_Daughter2(0),
q_Daughter2(0),
dcaxy_Daughter2(0),
dcaz_Daughter2(0),
nTPC_Clusters_Daughter2(0),
nTPC_FindableClusters_Daughter2(0),
nTPC_CrossedRows_Daughter2(0),
nTPC_Clusters_dEdx_Daughter2(0),
nITS_Clusters_Daughter2(0),
HasPointOnITSLayer0_Daughter2(0),
HasPointOnITSLayer1_Daughter2(0),
HasPointOnITSLayer2_Daughter2(0),
HasPointOnITSLayer3_Daughter2(0),
HasPointOnITSLayer4_Daughter2(0),
HasPointOnITSLayer5_Daughter2(0),
HasSharedPointOnITSLayer0_Daughter2(0),
HasSharedPointOnITSLayer1_Daughter2(0),
HasSharedPointOnITSLayer2_Daughter2(0),
HasSharedPointOnITSLayer3_Daughter2(0),
HasSharedPointOnITSLayer4_Daughter2(0),
HasSharedPointOnITSLayer5_Daughter2(0),
chi2_TPC_Daughter2(0),
chi2_NDF_Daughter2(0),
chi2_ITS_Daughter2(0),
nSigmaITS_He3_Daughter2(0),
nSigmaTPC_He3_Daughter2(0),
nSigmaTOF_He3_Daughter2(0),
nSigmaITS_Pion_Daughter2(0),
nSigmaTPC_Pion_Daughter2(0),
nSigmaTOF_Pion_Daughter2(0),
nSigmaITS_Trit_Daughter2(0),
nSigmaTPC_Trit_Daughter2(0),
nSigmaTOF_Trit_Daughter2(0),
cosPointingAngle(0),
dcaV0Daughters(0),
radius(0),
DecayLength(0),
massV0(0),
alphaV0(0),
qtV0(0)
{
      fUtils = new AliAnalysisUtils();
}
//_______________________________________________________________________________________________________________________________________________________
AliAnalysisTaskReducedTreeNuclei::AliAnalysisTaskReducedTreeNuclei(const char *name):
AliAnalysisTaskSE(name),
fAODevent(NULL),
fAODVZERO(NULL),
fPIDResponse(NULL),
fAODeventCuts(),
fUtils(NULL),
fFillTri(kTRUE),
fFillHypTri(kTRUE),
fQAList(NULL),
fOutputList(NULL),
// TreeEventSelection(NULL),
histoEventSelection(NULL),
histoEventMultiplicity(NULL),
histoEventMultV0M_MB(NULL),
histoEventMultV0M_HM(NULL),
reducedTree_Helium(NULL),
reducedTree_Triton(NULL),
reducedTree_HyperTriton(NULL),
isTrigINT7(false),
isTrigHighMult(false),
// triggerMask(1),
SelectionStep(-1),
magFieldSign(0),
multPercentile_V0M(-1),
multPercentile_V0A(-1),
multPercentile_V0C(-1),
multPercentile_OnlineV0M(-1),
multPercentile_OnlineV0A(-1),
multPercentile_OnlineV0C(-1),
multPercentile_ADM(-1),
multPercentile_ADA(-1),
multPercentile_ADC(-1),
multPercentile_SPDClusters(-1),
multPercentile_SPDTracklets(-1),
multPercentile_RefMult05(-1),
multPercentile_RefMult08(-1),
multPercentile_CL1(-1),
multPercentile_ZNA(-1),
Ntrk_OnlineV0M(0),
Ntrk_OnlineV0A(0),
Ntrk_OnlineV0C(0),
Ntrk_V0M(0),
Ntrk_V0A(0),
Ntrk_V0C(0),
Ntrk_ADM(0),
Ntrk_ADA(0),
Ntrk_ADC(0),
Ntrk_SPDClusters(0),
Ntrk_SPDTracklets(0),
Ntrk_RefMult05(0),
Ntrk_RefMult08(0),
nVertexContributors(0),
xVertex(0),
yVertex(0),
zVertex(0),
px(0),
py(0),
pz(0),
TPCmomentum(0),
TRDmomentum(0),
integratedLength(0),
timeOfFlight(0),
beta(0),
gamma(0),
mass(0),
trackID(0),
eta(0),
phi(0),
theta(0),
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
HasSharedPointOnITSLayer0(0),
HasSharedPointOnITSLayer1(0),
HasSharedPointOnITSLayer2(0),
HasSharedPointOnITSLayer3(0),
HasSharedPointOnITSLayer4(0),
HasSharedPointOnITSLayer5(0),
chi2_TPC(0),
chi2_NDF(0),
chi2_ITS(0),
ITSsignal(0),
TPCsignal(0),
TOFsignal(0),
TRDsignal(0),
HMPIDsignal(0),
nSigmaITS_He4(0),
nSigmaTPC_He4(0),
nSigmaTOF_He4(0),
nSigmaTRD_He4(0),
nSigmaHMPID_He4(0),
nSigmaITS_He3(0),
nSigmaTPC_He3(0),
nSigmaTOF_He3(0),
nSigmaTRD_He3(0),
nSigmaHMPID_He3(0),
nSigmaITS_Trit(0),
nSigmaTPC_Trit(0),
nSigmaTOF_Trit(0),
nSigmaTRD_Trit(0),
nSigmaHMPID_Trit(0),
nSigmaITS_Deut(0),
nSigmaTPC_Deut(0),
nSigmaTOF_Deut(0),
nSigmaTRD_Deut(0),
nSigmaHMPID_Deut(0),
nSigmaITS_Prot(0),
nSigmaTPC_Prot(0),
nSigmaTOF_Prot(0),
nSigmaTRD_Prot(0),
nSigmaHMPID_Prot(0),
nSigmaITS_Pion(0),
nSigmaTPC_Pion(0),
nSigmaTOF_Pion(0),
nSigmaTRD_Pion(0),
nSigmaHMPID_Pion(0),
nSigmaITS_Kaon(0),
nSigmaTPC_Kaon(0),
nSigmaTOF_Kaon(0),
nSigmaTRD_Kaon(0),
nSigmaHMPID_Kaon(0),
nSigmaITS_Elec(0),
nSigmaTPC_Elec(0),
nSigmaTOF_Elec(0),
nSigmaTRD_Elec(0),
nSigmaHMPID_Elec(0),
px_Daughter1(0),
py_Daughter1(0),
pz_Daughter1(0),
TPCmomentum_Daughter1(0),
trackID_Daughter1(0),
eta_Daughter1(0),
phi_Daughter1(0),
theta_Daughter1(0),
y_Daughter1(0),
q_Daughter1(0),
dcaxy_Daughter1(0),
dcaz_Daughter1(0),
nTPC_Clusters_Daughter1(0),
nTPC_FindableClusters_Daughter1(0),
nTPC_CrossedRows_Daughter1(0),
nTPC_Clusters_dEdx_Daughter1(0),
nITS_Clusters_Daughter1(0),
HasPointOnITSLayer0_Daughter1(0),
HasPointOnITSLayer1_Daughter1(0),
HasPointOnITSLayer2_Daughter1(0),
HasPointOnITSLayer3_Daughter1(0),
HasPointOnITSLayer4_Daughter1(0),
HasPointOnITSLayer5_Daughter1(0),
HasSharedPointOnITSLayer0_Daughter1(0),
HasSharedPointOnITSLayer1_Daughter1(0),
HasSharedPointOnITSLayer2_Daughter1(0),
HasSharedPointOnITSLayer3_Daughter1(0),
HasSharedPointOnITSLayer4_Daughter1(0),
HasSharedPointOnITSLayer5_Daughter1(0),
chi2_TPC_Daughter1(0),
chi2_NDF_Daughter1(0),
chi2_ITS_Daughter1(0),
nSigmaITS_He3_Daughter1(0),
nSigmaTPC_He3_Daughter1(0),
nSigmaTOF_He3_Daughter1(0),
nSigmaITS_Pion_Daughter1(0),
nSigmaTPC_Pion_Daughter1(0),
nSigmaTOF_Pion_Daughter1(0),
nSigmaITS_Trit_Daughter1(0),
nSigmaTPC_Trit_Daughter1(0),
nSigmaTOF_Trit_Daughter1(0),
px_Daughter2(0),
py_Daughter2(0),
pz_Daughter2(0),
TPCmomentum_Daughter2(0),
trackID_Daughter2(0),
eta_Daughter2(0),
phi_Daughter2(0),
theta_Daughter2(0),
y_Daughter2(0),
q_Daughter2(0),
dcaxy_Daughter2(0),
dcaz_Daughter2(0),
nTPC_Clusters_Daughter2(0),
nTPC_FindableClusters_Daughter2(0),
nTPC_CrossedRows_Daughter2(0),
nTPC_Clusters_dEdx_Daughter2(0),
nITS_Clusters_Daughter2(0),
HasPointOnITSLayer0_Daughter2(0),
HasPointOnITSLayer1_Daughter2(0),
HasPointOnITSLayer2_Daughter2(0),
HasPointOnITSLayer3_Daughter2(0),
HasPointOnITSLayer4_Daughter2(0),
HasPointOnITSLayer5_Daughter2(0),
HasSharedPointOnITSLayer0_Daughter2(0),
HasSharedPointOnITSLayer1_Daughter2(0),
HasSharedPointOnITSLayer2_Daughter2(0),
HasSharedPointOnITSLayer3_Daughter2(0),
HasSharedPointOnITSLayer4_Daughter2(0),
HasSharedPointOnITSLayer5_Daughter2(0),
chi2_TPC_Daughter2(0),
chi2_NDF_Daughter2(0),
chi2_ITS_Daughter2(0),
nSigmaITS_He3_Daughter2(0),
nSigmaTPC_He3_Daughter2(0),
nSigmaTOF_He3_Daughter2(0),
nSigmaITS_Pion_Daughter2(0),
nSigmaTPC_Pion_Daughter2(0),
nSigmaTOF_Pion_Daughter2(0),
nSigmaITS_Trit_Daughter2(0),
nSigmaTPC_Trit_Daughter2(0),
nSigmaTOF_Trit_Daughter2(0),
cosPointingAngle(0),
dcaV0Daughters(0),
radius(0),
DecayLength(0),
massV0(0),
alphaV0(0),
qtV0(0)
{
  fUtils = new AliAnalysisUtils();
  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TTree::Class());
  DefineOutput(4, TTree::Class());
  DefineOutput(5, TTree::Class());
  DefineOutput(6, TTree::Class());
  
}
//_______________________________________________________________________________________________________________________________________________________
AliAnalysisTaskReducedTreeNuclei::~AliAnalysisTaskReducedTreeNuclei()
{
  if(fAODevent) delete fAODevent;
  if(fPIDResponse) delete fPIDResponse;
  if(fQAList) delete fQAList;
  if(fOutputList) delete fOutputList;
//   if(TreeEventSelection) delete TreeEventSelection;
  if(reducedTree_Helium) delete reducedTree_Helium;
  if(reducedTree_Triton) delete reducedTree_Triton;
  if(reducedTree_HyperTriton) delete reducedTree_HyperTriton;
  delete fUtils;
}
//_______________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskReducedTreeNuclei::UserCreateOutputObjects()
{
  fQAList = new TList();
  fQAList -> SetOwner();

  fOutputList = new TList();
  fOutputList -> SetOwner();

  fAODeventCuts.AddQAplotsToList(fQAList); /// Add event selection QA plots

//   TreeEventSelection = new TTree("TreeEventSelection","TreeEventSelection");
//   TreeEventSelection -> Branch("SelectionStep",&SelectionStep,"SelectionStep/I");
//   //   TreeEventSelection -> Branch("multPercentile_V0A",&multPercentile_V0A,"multPercentile_V0A/D");
//   TreeEventSelection -> Branch("multPercentile_V0M",&multPercentile_V0M,"multPercentile_V0M/D");
//   TreeEventSelection -> Branch("multPercentile_SPDTracklets",&multPercentile_SPDTracklets,"multPercentile_SPDTracklets/D");
//   TreeEventSelection -> Branch("isTrigINT7",&isTrigINT7,"isTrigINT7/O");
//   TreeEventSelection -> Branch("isTrigHighMult",&isTrigHighMult,"isTrigHighMult/O");
//   TreeEventSelection -> Branch("multPercentile_RefMult08",&multPercentile_RefMult08,"multPercentile_RefMult08/D"); //NOT FILLED FOR LHC16q (p--Pb)
//   TreeEventSelection -> Branch("Ntrk_RefMult08",&Ntrk_RefMult08,"Ntrk_RefMult08/I"); //NOT FILLED FOR LHC16q (p--Pb)
//   TreeEventSelection -> Branch("Ntrk_V0M",&Ntrk_V0M,"Ntrk_V0M/I");
//   TreeEventSelection -> Branch("Ntrk_SPDTracklets",&Ntrk_SPDTracklets,"Ntrk_SPDTracklets/I");

  
  histoEventSelection = new TH2D("histoEventSelection","Events after selection steps",4,0,4, 4,0,4);
  histoEventSelection->GetXaxis()->SetTitle("Selection steps");
  histoEventSelection->GetYaxis()->SetTitle("trigger");
  histoEventSelection->Sumw2();
  fOutputList->Add(histoEventSelection);

  histoEventMultiplicity = new TH2D("histoEventMultiplicity","Events vs multiplicity percentile (V0A)",1000,0,100, 4,0,4);
  histoEventMultiplicity->GetXaxis()->SetTitle("Multiplicity percentile");
  histoEventMultiplicity->GetYaxis()->SetTitle("trigger");
  histoEventMultiplicity->Sumw2();
  fOutputList->Add(histoEventMultiplicity);

  
  histoEventMultV0M_MB = new TH1D("histoEventMultV0M_MB", "HM V0M amplitude",1001,0,1000);
  histoEventMultV0M_MB->GetXaxis()->SetTitle("V0M amplitude");
  histoEventMultV0M_MB->GetYaxis()->SetTitle("counts");
  histoEventMultV0M_MB->Sumw2();
  fOutputList->Add(histoEventMultV0M_MB);
  
  histoEventMultV0M_HM = new TH1D("histoEventMultV0M_HM", "HM V0M amplitude",1001,0,1000);
  histoEventMultV0M_HM->GetXaxis()->SetTitle("V0M amplitude");
  histoEventMultV0M_HM->GetYaxis()->SetTitle("counts");
  histoEventMultV0M_HM->Sumw2();
  fOutputList->Add(histoEventMultV0M_HM);
  
  //Reduced Tree (Helium3)
  reducedTree_Helium = new TTree("reducedTree_Helium","reducedTree_Helium");
  //   reducedTree_Helium -> Branch("triggerMask",&triggerMask,"triggerMask/l");
  reducedTree_Helium -> Branch("isTrigINT7",&isTrigINT7,"isTrigINT7/O");
  reducedTree_Helium -> Branch("isTrigHighMult",&isTrigHighMult,"isTrigHighMult/O");
  reducedTree_Helium -> Branch("magFieldSign",&magFieldSign,"magFieldSign/I");
  reducedTree_Helium -> Branch("multPercentile_V0M",&multPercentile_V0M,"multPercentile_V0M/D");
  reducedTree_Helium -> Branch("multPercentile_V0A",&multPercentile_V0A,"multPercentile_V0A/D");
  //   reducedTree_Helium -> Branch("multPercentile_V0C",&multPercentile_V0C,"multPercentile_V0C/D");
  //   reducedTree_Helium -> Branch("multPercentile_OnlineV0M",&multPercentile_OnlineV0M,"multPercentile_OnlineV0M/D"); //NOT FILLED FOR LHC16q (p--Pb)
  //   reducedTree_Helium -> Branch("multPercentile_OnlineV0A",&multPercentile_OnlineV0A,"multPercentile_OnlineV0A/D"); //NOT FILLED FOR LHC16q (p--Pb)
  //   reducedTree_Helium -> Branch("multPercentile_OnlineV0C",&multPercentile_OnlineV0C,"multPercentile_OnlineV0C/D"); //NOT FILLED FOR LHC16q (p--Pb)
  //   reducedTree_Helium -> Branch("multPercentile_ADM",&multPercentile_ADM,"multPercentile_ADM/D"); //NOT FILLED FOR LHC16q (p--Pb)
  //   reducedTree_Helium -> Branch("multPercentile_ADA",&multPercentile_ADA,"multPercentile_ADA/D"); //NOT FILLED FOR LHC16q (p--Pb)
  //   reducedTree_Helium -> Branch("multPercentile_ADC",&multPercentile_ADC,"multPercentile_ADC/D"); //NOT FILLED FOR LHC16q (p--Pb)
  //   reducedTree_Helium -> Branch("multPercentile_SPDClusters",&multPercentile_SPDClusters,"multPercentile_SPDClusters/D");  //NOT FILLED FOR LHC16q (p--Pb)
  reducedTree_Helium -> Branch("multPercentile_SPDTracklets",&multPercentile_SPDTracklets,"multPercentile_SPDTracklets/D");
  //   reducedTree_Helium -> Branch("multPercentile_RefMult05",&multPercentile_RefMult05,"multPercentile_RefMult05/D");  //NOT FILLED FOR LHC16q (p--Pb)
  //   reducedTree_Helium -> Branch("multPercentile_RefMult08",&multPercentile_RefMult08,"multPercentile_RefMult08/D"); //NOT FILLED FOR LHC16q (p--Pb)
  reducedTree_Helium -> Branch("multPercentile_CL1",&multPercentile_CL1,"multPercentile_CL1/D");
  reducedTree_Helium -> Branch("multPercentile_ZNA",&multPercentile_ZNA,"multPercentile_ZNA/D");
  //   reducedTree_Helium -> Branch("Ntrk_OnlineV0M",&Ntrk_OnlineV0M,"Ntrk_OnlineV0M/I"); //NOT FILLED FOR LHC16q (p--Pb)
  //   reducedTree_Helium -> Branch("Ntrk_OnlineV0A",&Ntrk_OnlineV0A,"Ntrk_OnlineV0A/I"); //NOT FILLED FOR LHC16q (p--Pb)
  //   reducedTree_Helium -> Branch("Ntrk_OnlineV0C",&Ntrk_OnlineV0C,"Ntrk_OnlineV0C/I"); //NOT FILLED FOR LHC16q (p--Pb)
  reducedTree_Helium -> Branch("Ntrk_V0M",&Ntrk_V0M,"Ntrk_V0M/I");
  //   reducedTree_Helium -> Branch("Ntrk_V0A",&Ntrk_V0A,"Ntrk_V0A/I");
  //   reducedTree_Helium -> Branch("Ntrk_V0C",&Ntrk_V0C,"Ntrk_V0C/I");
  //   reducedTree_Helium -> Branch("Ntrk_ADM",&Ntrk_ADM,"Ntrk_ADM/I"); //NOT FILLED FOR LHC16q (p--Pb)
  //   reducedTree_Helium -> Branch("Ntrk_ADA",&Ntrk_ADA,"Ntrk_ADA/I"); //NOT FILLED FOR LHC16q (p--Pb)
  //   reducedTree_Helium -> Branch("Ntrk_ADC",&Ntrk_ADC,"Ntrk_ADC/I"); //NOT FILLED FOR LHC16q (p--Pb)
  //   reducedTree_Helium -> Branch("Ntrk_SPDClusters",&Ntrk_SPDClusters,"Ntrk_SPDClusters/I"); //NOT FILLED FOR LHC16q (p--Pb)
  reducedTree_Helium -> Branch("Ntrk_SPDTracklets",&Ntrk_SPDTracklets,"Ntrk_SPDTracklets/I");
  reducedTree_Helium -> Branch("Ntrk_RefMult05",&Ntrk_RefMult05,"Ntrk_RefMult05/I"); //NOT FILLED FOR LHC16q (p--Pb)
  reducedTree_Helium -> Branch("Ntrk_RefMult08",&Ntrk_RefMult08,"Ntrk_RefMult08/I"); //NOT FILLED FOR LHC16q (p--Pb)
  reducedTree_Helium -> Branch("nVertexContributors",&nVertexContributors,"nVertexContributors/I");
  reducedTree_Helium -> Branch("xVertex",&xVertex,"xVertex/D");
  reducedTree_Helium -> Branch("yVertex",&yVertex,"yVertex/D");
  reducedTree_Helium -> Branch("zVertex",&zVertex,"zVertex/D");
  reducedTree_Helium -> Branch("px",&px,"px/D");
  reducedTree_Helium -> Branch("py",&py,"py/D");
  reducedTree_Helium -> Branch("pz",&pz,"pz/D");
  reducedTree_Helium -> Branch("TPCmomentum",&TPCmomentum,"TPCmomentum/D");
  reducedTree_Helium -> Branch("TRDmomentum",&TRDmomentum,"TRDmomentum/D");
  reducedTree_Helium -> Branch("integratedLength",&integratedLength,"integratedLength/D");
  reducedTree_Helium -> Branch("timeOfFlight",&timeOfFlight,"timeOfFlight/D");
  reducedTree_Helium -> Branch("beta",&beta,"beta/D");
  reducedTree_Helium -> Branch("gamma",&gamma,"gamma/D");
  reducedTree_Helium -> Branch("mass",&mass,"mass/D");
  reducedTree_Helium -> Branch("trackID",&trackID,"trackID/I");
  reducedTree_Helium -> Branch("eta",&eta,"eta/D");
  reducedTree_Helium -> Branch("phi",&phi,"phi/D");
  reducedTree_Helium -> Branch("theta",&theta,"theta/D");
  reducedTree_Helium -> Branch("y",&y,"y/D");
  reducedTree_Helium -> Branch("q",&q,"q/I");
  reducedTree_Helium -> Branch("dcaxy",&dcaxy,"dcaxy/D");
  reducedTree_Helium -> Branch("dcaz",&dcaz,"dcaz/D");
  reducedTree_Helium -> Branch("nTPC_Clusters",&nTPC_Clusters,"nTPC_Clusters/I");
  reducedTree_Helium -> Branch("nTRD_Clusters",&nTRD_Clusters,"nTRD_Clusters/I");
  reducedTree_Helium -> Branch("nITS_Clusters",&nITS_Clusters,"nITS_Clusters/I");
  reducedTree_Helium -> Branch("nTPC_FindableClusters",&nTPC_FindableClusters,"nTPC_FindableClusters/I");
  reducedTree_Helium -> Branch("nTPC_CrossedRows",&nTPC_CrossedRows,"nTPC_CrossedRows/I");
  reducedTree_Helium -> Branch("nTPC_Clusters_dEdx",&nTPC_Clusters_dEdx,"nTPC_Clusters_dEdx/I");
  reducedTree_Helium -> Branch("HasPointOnITSLayer0",&HasPointOnITSLayer0,"HasPointOnITSLayer0/O");
  reducedTree_Helium -> Branch("HasPointOnITSLayer1",&HasPointOnITSLayer1,"HasPointOnITSLayer1/O");
  reducedTree_Helium -> Branch("HasPointOnITSLayer2",&HasPointOnITSLayer2,"HasPointOnITSLayer2/O");
  reducedTree_Helium -> Branch("HasPointOnITSLayer3",&HasPointOnITSLayer3,"HasPointOnITSLayer3/O");
  reducedTree_Helium -> Branch("HasPointOnITSLayer4",&HasPointOnITSLayer4,"HasPointOnITSLayer4/O");
  reducedTree_Helium -> Branch("HasPointOnITSLayer5",&HasPointOnITSLayer5,"HasPointOnITSLayer5/O");
  reducedTree_Helium -> Branch("HasSharedPointOnITSLayer0",&HasSharedPointOnITSLayer0,"HasSharedPointOnITSLayer0/O");
  reducedTree_Helium -> Branch("HasSharedPointOnITSLayer1",&HasSharedPointOnITSLayer1,"HasSharedPointOnITSLayer1/O");
  reducedTree_Helium -> Branch("HasSharedPointOnITSLayer2",&HasSharedPointOnITSLayer2,"HasSharedPointOnITSLayer2/O");
  reducedTree_Helium -> Branch("HasSharedPointOnITSLayer3",&HasSharedPointOnITSLayer3,"HasSharedPointOnITSLayer3/O");
  reducedTree_Helium -> Branch("HasSharedPointOnITSLayer4",&HasSharedPointOnITSLayer4,"HasSharedPointOnITSLayer4/O");
  reducedTree_Helium -> Branch("HasSharedPointOnITSLayer5",&HasSharedPointOnITSLayer5,"HasSharedPointOnITSLayer5/O");
  reducedTree_Helium -> Branch("chi2_TPC",&chi2_TPC,"chi2_TPC/D");
  reducedTree_Helium -> Branch("chi2_NDF",&chi2_NDF,"chi2_NDF/D");
  reducedTree_Helium -> Branch("chi2_ITS",&chi2_ITS,"chi2_ITS/D");
  reducedTree_Helium -> Branch("ITSsignal",&ITSsignal,"ITSsignal/D");
  reducedTree_Helium -> Branch("TPCsignal",&TPCsignal,"TPCsignal/D");
  reducedTree_Helium -> Branch("TOFsignal",&TOFsignal,"TOFsignal/D");
  reducedTree_Helium -> Branch("TRDsignal",&TRDsignal,"TRDsignal/D");
  reducedTree_Helium -> Branch("HMPIDsignal",&HMPIDsignal,"HMPIDsignal/D");
  reducedTree_Helium -> Branch("nSigmaITS_He4",&nSigmaITS_He4,"nSigmaITS_He4/D");
  reducedTree_Helium -> Branch("nSigmaTPC_He4",&nSigmaTPC_He4,"nSigmaTPC_He4/D");
  reducedTree_Helium -> Branch("nSigmaTOF_He4",&nSigmaTOF_He4,"nSigmaTOF_He4/D");
  reducedTree_Helium -> Branch("nSigmaTRD_He4",&nSigmaTRD_He4,"nSigmaTRD_He4/D");
  reducedTree_Helium -> Branch("nSigmaHMPID_He4",&nSigmaHMPID_He4,"nSigmaHMPID_He4/D");
  reducedTree_Helium -> Branch("nSigmaITS_He3",&nSigmaITS_He3,"nSigmaITS_He3/D");
  reducedTree_Helium -> Branch("nSigmaTPC_He3",&nSigmaTPC_He3,"nSigmaTPC_He3/D");
  reducedTree_Helium -> Branch("nSigmaTOF_He3",&nSigmaTOF_He3,"nSigmaTOF_He3/D");
  reducedTree_Helium -> Branch("nSigmaTRD_He3",&nSigmaTRD_He3,"nSigmaTRD_He3/D");
  reducedTree_Helium -> Branch("nSigmaHMPID_He3",&nSigmaHMPID_He3,"nSigmaHMPID_He3/D");
  reducedTree_Helium -> Branch("nSigmaITS_Trit",&nSigmaITS_Trit,"nSigmaITS_Trit/D");
  reducedTree_Helium -> Branch("nSigmaTPC_Trit",&nSigmaTPC_Trit,"nSigmaTPC_Trit/D");
  reducedTree_Helium -> Branch("nSigmaTOF_Trit",&nSigmaTOF_Trit,"nSigmaTOF_Trit/D");
  reducedTree_Helium -> Branch("nSigmaTRD_Trit",&nSigmaTRD_Trit,"nSigmaTRD_Trit/D");
  reducedTree_Helium -> Branch("nSigmaHMPID_Trit",&nSigmaHMPID_Trit,"nSigmaHMPID_Trit/D");
  //   reducedTree_Helium -> Branch("nSigmaITS_Deut",&nSigmaITS_Deut,"nSigmaITS_Deut/D");
  //   reducedTree_Helium -> Branch("nSigmaTPC_Deut",&nSigmaTPC_Deut,"nSigmaTPC_Deut/D");
  //   reducedTree_Helium -> Branch("nSigmaTOF_Deut",&nSigmaTOF_Deut,"nSigmaTOF_Deut/D");
  //   reducedTree_Helium -> Branch("nSigmaTRD_Deut",&nSigmaTRD_Deut,"nSigmaTRD_Deut/D");
  //   reducedTree_Helium -> Branch("nSigmaHMPID_Deut",&nSigmaHMPID_Deut,"nSigmaHMPID_Deut/D");
  //   reducedTree_Helium -> Branch("nSigmaITS_Prot",&nSigmaITS_Prot,"nSigmaITS_Prot/D");
  //   reducedTree_Helium -> Branch("nSigmaTPC_Prot",&nSigmaTPC_Prot,"nSigmaTPC_Prot/D");
  //   reducedTree_Helium -> Branch("nSigmaTOF_Prot",&nSigmaTOF_Prot,"nSigmaTOF_Prot/D");
  //   reducedTree_Helium -> Branch("nSigmaTRD_Prot",&nSigmaTRD_Prot,"nSigmaTRD_Prot/D");
  //   reducedTree_Helium -> Branch("nSigmaHMPID_Prot",&nSigmaHMPID_Prot,"nSigmaHMPID_Prot/D");
  //   reducedTree_Helium -> Branch("nSigmaITS_Pion",&nSigmaITS_Pion,"nSigmaITS_Pion/D");
  //   reducedTree_Helium -> Branch("nSigmaTPC_Pion",&nSigmaTPC_Pion,"nSigmaTPC_Pion/D");
  //   reducedTree_Helium -> Branch("nSigmaTOF_Pion",&nSigmaTOF_Pion,"nSigmaTOF_Pion/D");
  //   reducedTree_Helium -> Branch("nSigmaTRD_Pion",&nSigmaTRD_Pion,"nSigmaTRD_Pion/D");
  //   reducedTree_Helium -> Branch("nSigmaHMPID_Pion",&nSigmaHMPID_Pion,"nSigmaHMPID_Pion/D");
  //   reducedTree_Helium -> Branch("nSigmaITS_Kaon",&nSigmaITS_Kaon,"nSigmaITS_Kaon/D");
  //   reducedTree_Helium -> Branch("nSigmaTPC_Kaon",&nSigmaTPC_Kaon,"nSigmaTPC_Kaon/D");
  //   reducedTree_Helium -> Branch("nSigmaTOF_Kaon",&nSigmaTOF_Kaon,"nSigmaTOF_Kaon/D");
  //   reducedTree_Helium -> Branch("nSigmaTRD_Kaon",&nSigmaTRD_Kaon,"nSigmaTRD_Kaon/D");
  //   reducedTree_Helium -> Branch("nSigmaHMPID_Kaon",&nSigmaHMPID_Kaon,"nSigmaHMPID_Kaon/D");
  //   reducedTree_Helium -> Branch("nSigmaITS_Elec",&nSigmaITS_Elec,"nSigmaITS_Elec/D");
  //   reducedTree_Helium -> Branch("nSigmaTPC_Elec",&nSigmaTPC_Elec,"nSigmaTPC_Elec/D");
  //   reducedTree_Helium -> Branch("nSigmaTOF_Elec",&nSigmaTOF_Elec,"nSigmaTOF_Elec/D");
  //   reducedTree_Helium -> Branch("nSigmaTRD_Elec",&nSigmaTRD_Elec,"nSigmaTRD_Elec/D");
  //   reducedTree_Helium -> Branch("nSigmaHMPID_Elec",&nSigmaHMPID_Elec,"nSigmaHMPID_Elec/D");


  //Reduced Tree (Triton)
  if (fFillTri){
    reducedTree_Triton = new TTree("reducedTree_Triton","reducedTree_Triton");
    reducedTree_Triton -> Branch("magFieldSign",&magFieldSign,"magFieldSign/I");
    reducedTree_Triton -> Branch("multPercentile_V0M",&multPercentile_V0M,"multPercentile_V0M/D");
    reducedTree_Triton -> Branch("multPercentile_V0A",&multPercentile_V0A,"multPercentile_V0A/D");
    reducedTree_Triton -> Branch("multPercentile_V0C",&multPercentile_V0C,"multPercentile_V0C/D");
    
    reducedTree_Triton -> Branch("isTrigINT7",&isTrigINT7,"isTrigINT7/O");
    reducedTree_Triton -> Branch("isTrigHighMult",&isTrigHighMult,"isTrigHighMult/O");
  
    //   reducedTree_Triton -> Branch("multPercentile_OnlineV0M",&multPercentile_OnlineV0M,"multPercentile_OnlineV0M/D"); //NOT FILLED FOR LHC16q (p--Pb)
    //   reducedTree_Triton -> Branch("multPercentile_OnlineV0A",&multPercentile_OnlineV0A,"multPercentile_OnlineV0A/D"); //NOT FILLED FOR LHC16q (p--Pb)
    //   reducedTree_Triton -> Branch("multPercentile_OnlineV0C",&multPercentile_OnlineV0C,"multPercentile_OnlineV0C/D"); //NOT FILLED FOR LHC16q (p--Pb)
    //   reducedTree_Triton -> Branch("multPercentile_ADM",&multPercentile_ADM,"multPercentile_ADM/D"); //NOT FILLED FOR LHC16q (p--Pb)
    //   reducedTree_Triton -> Branch("multPercentile_ADA",&multPercentile_ADA,"multPercentile_ADA/D"); //NOT FILLED FOR LHC16q (p--Pb)
    //   reducedTree_Triton -> Branch("multPercentile_ADC",&multPercentile_ADC,"multPercentile_ADC/D"); //NOT FILLED FOR LHC16q (p--Pb)
    //   reducedTree_Triton -> Branch("multPercentile_SPDClusters",&multPercentile_SPDClusters,"multPercentile_SPDClusters/D");  //NOT FILLED FOR LHC16q (p--Pb)
    reducedTree_Triton -> Branch("multPercentile_SPDTracklets",&multPercentile_SPDTracklets,"multPercentile_SPDTracklets/D");
    //   reducedTree_Triton -> Branch("multPercentile_RefMult05",&multPercentile_RefMult05,"multPercentile_RefMult05/D");  //NOT FILLED FOR LHC16q (p--Pb)
    //   reducedTree_Triton -> Branch("multPercentile_RefMult08",&multPercentile_RefMult08,"multPercentile_RefMult08/D"); //NOT FILLED FOR LHC16q (p--Pb)
    reducedTree_Triton -> Branch("multPercentile_CL1",&multPercentile_CL1,"multPercentile_CL1/D");
    reducedTree_Triton -> Branch("multPercentile_ZNA",&multPercentile_ZNA,"multPercentile_ZNA/D");
    //   reducedTree_Triton -> Branch("Ntrk_OnlineV0M",&Ntrk_OnlineV0M,"Ntrk_OnlineV0M/I"); //NOT FILLED FOR LHC16q (p--Pb)
    //   reducedTree_Triton -> Branch("Ntrk_OnlineV0A",&Ntrk_OnlineV0A,"Ntrk_OnlineV0A/I"); //NOT FILLED FOR LHC16q (p--Pb)
    //   reducedTree_Triton -> Branch("Ntrk_OnlineV0C",&Ntrk_OnlineV0C,"Ntrk_OnlineV0C/I"); //NOT FILLED FOR LHC16q (p--Pb)
    reducedTree_Triton -> Branch("Ntrk_V0M",&Ntrk_V0M,"Ntrk_V0M/I");
    reducedTree_Triton -> Branch("Ntrk_V0A",&Ntrk_V0A,"Ntrk_V0A/I");
    reducedTree_Triton -> Branch("Ntrk_V0C",&Ntrk_V0C,"Ntrk_V0C/I");
    //   reducedTree_Triton -> Branch("Ntrk_ADM",&Ntrk_ADM,"Ntrk_ADM/I"); //NOT FILLED FOR LHC16q (p--Pb)
    //   reducedTree_Triton -> Branch("Ntrk_ADA",&Ntrk_ADA,"Ntrk_ADA/I"); //NOT FILLED FOR LHC16q (p--Pb)
    //   reducedTree_Triton -> Branch("Ntrk_ADC",&Ntrk_ADC,"Ntrk_ADC/I"); //NOT FILLED FOR LHC16q (p--Pb)
    //   reducedTree_Triton -> Branch("Ntrk_SPDClusters",&Ntrk_SPDClusters,"Ntrk_SPDClusters/I"); //NOT FILLED FOR LHC16q (p--Pb)
    reducedTree_Triton -> Branch("Ntrk_SPDTracklets",&Ntrk_SPDTracklets,"Ntrk_SPDTracklets/I");
    //   reducedTree_Triton -> Branch("Ntrk_RefMult05",&Ntrk_RefMult05,"Ntrk_RefMult05/I"); //NOT FILLED FOR LHC16q (p--Pb)
    //   reducedTree_Triton -> Branch("Ntrk_RefMult08",&Ntrk_RefMult08,"Ntrk_RefMult08/I"); //NOT FILLED FOR LHC16q (p--Pb)
    reducedTree_Triton -> Branch("nVertexContributors",&nVertexContributors,"nVertexContributors/I");
    reducedTree_Triton -> Branch("xVertex",&xVertex,"xVertex/D");
    reducedTree_Triton -> Branch("yVertex",&yVertex,"yVertex/D");
    reducedTree_Triton -> Branch("zVertex",&zVertex,"zVertex/D");
    reducedTree_Triton -> Branch("px",&px,"px/D");
    reducedTree_Triton -> Branch("py",&py,"py/D");
    reducedTree_Triton -> Branch("pz",&pz,"pz/D");
    reducedTree_Triton -> Branch("TPCmomentum",&TPCmomentum,"TPCmomentum/D");
    reducedTree_Triton -> Branch("TRDmomentum",&TRDmomentum,"TRDmomentum/D");
    reducedTree_Triton -> Branch("integratedLength",&integratedLength,"integratedLength/D");
    reducedTree_Triton -> Branch("timeOfFlight",&timeOfFlight,"timeOfFlight/D");
    reducedTree_Triton -> Branch("beta",&beta,"beta/D");
    reducedTree_Triton -> Branch("gamma",&gamma,"gamma/D");
    reducedTree_Triton -> Branch("mass",&mass,"mass/D");
    reducedTree_Triton -> Branch("trackID",&trackID,"trackID/I");
    reducedTree_Triton -> Branch("eta",&eta,"eta/D");
    reducedTree_Triton -> Branch("phi",&phi,"phi/D");
    reducedTree_Triton -> Branch("theta",&theta,"theta/D");
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
    reducedTree_Triton -> Branch("HasSharedPointOnITSLayer0",&HasSharedPointOnITSLayer0,"HasSharedPointOnITSLayer0/O");
    reducedTree_Triton -> Branch("HasSharedPointOnITSLayer1",&HasSharedPointOnITSLayer1,"HasSharedPointOnITSLayer1/O");
    reducedTree_Triton -> Branch("HasSharedPointOnITSLayer2",&HasSharedPointOnITSLayer2,"HasSharedPointOnITSLayer2/O");
    reducedTree_Triton -> Branch("HasSharedPointOnITSLayer3",&HasSharedPointOnITSLayer3,"HasSharedPointOnITSLayer3/O");
    reducedTree_Triton -> Branch("HasSharedPointOnITSLayer4",&HasSharedPointOnITSLayer4,"HasSharedPointOnITSLayer4/O");
    reducedTree_Triton -> Branch("HasSharedPointOnITSLayer5",&HasSharedPointOnITSLayer5,"HasSharedPointOnITSLayer5/O");
    reducedTree_Triton -> Branch("chi2_TPC",&chi2_TPC,"chi2_TPC/D");
    reducedTree_Triton -> Branch("chi2_NDF",&chi2_NDF,"chi2_NDF/D");
    reducedTree_Triton -> Branch("chi2_ITS",&chi2_ITS,"chi2_ITS/D");
    reducedTree_Triton -> Branch("ITSsignal",&ITSsignal,"ITSsignal/D");
    reducedTree_Triton -> Branch("TPCsignal",&TPCsignal,"TPCsignal/D");
    reducedTree_Triton -> Branch("TOFsignal",&TOFsignal,"TOFsignal/D");
    reducedTree_Triton -> Branch("TRDsignal",&TRDsignal,"TRDsignal/D");
    reducedTree_Triton -> Branch("HMPIDsignal",&HMPIDsignal,"HMPIDsignal/D");
    //   reducedTree_Triton -> Branch("nSigmaITS_He4",&nSigmaITS_He4,"nSigmaITS_He4/D");
    //   reducedTree_Triton -> Branch("nSigmaTPC_He4",&nSigmaTPC_He4,"nSigmaTPC_He4/D");
    //   reducedTree_Triton -> Branch("nSigmaTOF_He4",&nSigmaTOF_He4,"nSigmaTOF_He4/D");
    //   reducedTree_Triton -> Branch("nSigmaTRD_He4",&nSigmaTRD_He4,"nSigmaTRD_He4/D");
    //   reducedTree_Triton -> Branch("nSigmaHMPID_He4",&nSigmaHMPID_He4,"nSigmaHMPID_He4/D");
    reducedTree_Triton -> Branch("nSigmaITS_He3",&nSigmaITS_He3,"nSigmaITS_He3/D");
    reducedTree_Triton -> Branch("nSigmaTPC_He3",&nSigmaTPC_He3,"nSigmaTPC_He3/D");
    reducedTree_Triton -> Branch("nSigmaTOF_He3",&nSigmaTOF_He3,"nSigmaTOF_He3/D");
    reducedTree_Triton -> Branch("nSigmaTRD_He3",&nSigmaTRD_He3,"nSigmaTRD_He3/D");
    reducedTree_Triton -> Branch("nSigmaHMPID_He3",&nSigmaHMPID_He3,"nSigmaHMPID_He3/D");
    reducedTree_Triton -> Branch("nSigmaITS_Trit",&nSigmaITS_Trit,"nSigmaITS_Trit/D");
    reducedTree_Triton -> Branch("nSigmaTPC_Trit",&nSigmaTPC_Trit,"nSigmaTPC_Trit/D");
    reducedTree_Triton -> Branch("nSigmaTOF_Trit",&nSigmaTOF_Trit,"nSigmaTOF_Trit/D");
    reducedTree_Triton -> Branch("nSigmaTRD_Trit",&nSigmaTRD_Trit,"nSigmaTRD_Trit/D");
    reducedTree_Triton -> Branch("nSigmaHMPID_Trit",&nSigmaHMPID_Trit,"nSigmaHMPID_Trit/D");
    reducedTree_Triton -> Branch("nSigmaITS_Deut",&nSigmaITS_Deut,"nSigmaITS_Deut/D");
    reducedTree_Triton -> Branch("nSigmaTPC_Deut",&nSigmaTPC_Deut,"nSigmaTPC_Deut/D");
    reducedTree_Triton -> Branch("nSigmaTOF_Deut",&nSigmaTOF_Deut,"nSigmaTOF_Deut/D");
    reducedTree_Triton -> Branch("nSigmaTRD_Deut",&nSigmaTRD_Deut,"nSigmaTRD_Deut/D");
    reducedTree_Triton -> Branch("nSigmaHMPID_Deut",&nSigmaHMPID_Deut,"nSigmaHMPID_Deut/D");
    //      reducedTree_Triton -> Branch("nSigmaITS_Prot",&nSigmaITS_Prot,"nSigmaITS_Prot/D");
    //      reducedTree_Triton -> Branch("nSigmaTPC_Prot",&nSigmaTPC_Prot,"nSigmaTPC_Prot/D");
    //      reducedTree_Triton -> Branch("nSigmaTOF_Prot",&nSigmaTOF_Prot,"nSigmaTOF_Prot/D");
    //      reducedTree_Triton -> Branch("nSigmaTRD_Prot",&nSigmaTRD_Prot,"nSigmaTRD_Prot/D");
    //      reducedTree_Triton -> Branch("nSigmaHMPID_Prot",&nSigmaHMPID_Prot,"nSigmaHMPID_Prot/D");
    //      reducedTree_Triton -> Branch("nSigmaITS_Pion",&nSigmaITS_Pion,"nSigmaITS_Pion/D");
    //      reducedTree_Triton -> Branch("nSigmaTPC_Pion",&nSigmaTPC_Pion,"nSigmaTPC_Pion/D");
    //      reducedTree_Triton -> Branch("nSigmaTOF_Pion",&nSigmaTOF_Pion,"nSigmaTOF_Pion/D");
    //      reducedTree_Triton -> Branch("nSigmaTRD_Pion",&nSigmaTRD_Pion,"nSigmaTRD_Pion/D");
    //      reducedTree_Triton -> Branch("nSigmaHMPID_Pion",&nSigmaHMPID_Pion,"nSigmaHMPID_Pion/D");
    //      reducedTree_Triton -> Branch("nSigmaITS_Kaon",&nSigmaITS_Kaon,"nSigmaITS_Kaon/D");
    //      reducedTree_Triton -> Branch("nSigmaTPC_Kaon",&nSigmaTPC_Kaon,"nSigmaTPC_Kaon/D");
    //      reducedTree_Triton -> Branch("nSigmaTOF_Kaon",&nSigmaTOF_Kaon,"nSigmaTOF_Kaon/D");
    //      reducedTree_Triton -> Branch("nSigmaTRD_Kaon",&nSigmaTRD_Kaon,"nSigmaTRD_Kaon/D");
    //      reducedTree_Triton -> Branch("nSigmaHMPID_Kaon",&nSigmaHMPID_Kaon,"nSigmaHMPID_Kaon/D");
    //      reducedTree_Triton -> Branch("nSigmaITS_Elec",&nSigmaITS_Elec,"nSigmaITS_Elec/D");
    //      reducedTree_Triton -> Branch("nSigmaTPC_Elec",&nSigmaTPC_Elec,"nSigmaTPC_Elec/D");
    //      reducedTree_Triton -> Branch("nSigmaTOF_Elec",&nSigmaTOF_Elec,"nSigmaTOF_Elec/D");
    //      reducedTree_Triton -> Branch("nSigmaTRD_Elec",&nSigmaTRD_Elec,"nSigmaTRD_Elec/D");
    //      reducedTree_Triton -> Branch("nSigmaHMPID_Elec",&nSigmaHMPID_Elec,"nSigmaHMPID_Elec/D");
  }

  if (fFillHypTri){
    //Reduced Tree HyperTriton
    reducedTree_HyperTriton = new TTree("reducedTree_HyperTriton","reducedTree_HyperTriton");
    reducedTree_HyperTriton -> Branch("px_Daughter1",&px_Daughter1,"px_Daughter1/D");
    reducedTree_HyperTriton -> Branch("py_Daughter1",&py_Daughter1,"py_Daughter1/D");
    reducedTree_HyperTriton -> Branch("pz_Daughter1",&pz_Daughter1,"pz_Daughter1/D");
    reducedTree_HyperTriton -> Branch("TPCmomentum_Daughter1",&TPCmomentum_Daughter1,"TPCmomentum_Daughter1/D");
    reducedTree_HyperTriton -> Branch("trackID_Daughter1",&trackID_Daughter1,"trackID_Daughter1/I");
    reducedTree_HyperTriton -> Branch("eta_Daughter1",&eta_Daughter1,"eta_Daughter1/D");
    reducedTree_HyperTriton -> Branch("phi_Daughter1",&phi_Daughter1,"phi_Daughter1/D");
    reducedTree_HyperTriton -> Branch("theta_Daughter1",&theta_Daughter1,"theta_Daughter1/D");
    reducedTree_HyperTriton -> Branch("y_Daughter1",&y_Daughter1,"y_Daughter1/D");
    reducedTree_HyperTriton -> Branch("q_Daughter1",&q_Daughter1,"q_Daughter1/I");
    reducedTree_HyperTriton -> Branch("dcaxy_Daughter1",&dcaxy_Daughter1,"dcaxy_Daughter1/D");
    reducedTree_HyperTriton -> Branch("dcaz_Daughter1",&dcaz_Daughter1,"dcaz_Daughter1/D");
    reducedTree_HyperTriton -> Branch("nTPC_Clusters_Daughter1",&nTPC_Clusters_Daughter1,"nTPC_Clusters_Daughter1/I");
    reducedTree_HyperTriton -> Branch("nTPC_FindableClusters_Daughter1",&nTPC_FindableClusters_Daughter1,"nTPC_FindableClusters_Daughter1/I");
    reducedTree_HyperTriton -> Branch("nTPC_CrossedRows_Daughter1",&nTPC_CrossedRows_Daughter1,"nTPC_CrossedRows_Daughter1/I");
    reducedTree_HyperTriton -> Branch("nTPC_Clusters_dEdx_Daughter1",&nTPC_Clusters_dEdx_Daughter1,"nTPC_Clusters_dEdx_Daughter1/I");
    reducedTree_HyperTriton -> Branch("nITS_Clusters_Daughter1",&nITS_Clusters_Daughter1,"nITS_Clusters_Daughter1/I");
    reducedTree_HyperTriton -> Branch("HasPointOnITSLayer0_Daughter1",&HasPointOnITSLayer0_Daughter1,"HasPointOnITSLayer0_Daughter1/O");
    reducedTree_HyperTriton -> Branch("HasPointOnITSLayer1_Daughter1",&HasPointOnITSLayer1_Daughter1,"HasPointOnITSLayer1_Daughter1/O");
    reducedTree_HyperTriton -> Branch("HasPointOnITSLayer2_Daughter1",&HasPointOnITSLayer2_Daughter1,"HasPointOnITSLayer2_Daughter1/O");
    reducedTree_HyperTriton -> Branch("HasPointOnITSLayer3_Daughter1",&HasPointOnITSLayer3_Daughter1,"HasPointOnITSLayer3_Daughter1/O");
    reducedTree_HyperTriton -> Branch("HasPointOnITSLayer4_Daughter1",&HasPointOnITSLayer4_Daughter1,"HasPointOnITSLayer4_Daughter1/O");
    reducedTree_HyperTriton -> Branch("HasPointOnITSLayer5_Daughter1",&HasPointOnITSLayer5_Daughter1,"HasPointOnITSLayer5_Daughter1/O");
    reducedTree_HyperTriton -> Branch("HasSharedPointOnITSLayer0_Daughter1",&HasSharedPointOnITSLayer0_Daughter1,"HasSharedPointOnITSLayer0_Daughter1/O");
    reducedTree_HyperTriton -> Branch("HasSharedPointOnITSLayer1_Daughter1",&HasSharedPointOnITSLayer1_Daughter1,"HasSharedPointOnITSLayer1_Daughter1/O");
    reducedTree_HyperTriton -> Branch("HasSharedPointOnITSLayer2_Daughter1",&HasSharedPointOnITSLayer2_Daughter1,"HasSharedPointOnITSLayer2_Daughter1/O");
    reducedTree_HyperTriton -> Branch("HasSharedPointOnITSLayer3_Daughter1",&HasSharedPointOnITSLayer3_Daughter1,"HasSharedPointOnITSLayer3_Daughter1/O");
    reducedTree_HyperTriton -> Branch("HasSharedPointOnITSLayer4_Daughter1",&HasSharedPointOnITSLayer4_Daughter1,"HasSharedPointOnITSLayer4_Daughter1/O");
    reducedTree_HyperTriton -> Branch("HasSharedPointOnITSLayer5_Daughter1",&HasSharedPointOnITSLayer5_Daughter1,"HasSharedPointOnITSLayer5_Daughter1/O");
    reducedTree_HyperTriton -> Branch("chi2_TPC_Daughter1",&chi2_TPC_Daughter1,"chi2_TPC_Daughter1/D");
    reducedTree_HyperTriton -> Branch("chi2_NDF_Daughter1",&chi2_NDF_Daughter1,"chi2_NDF_Daughter1/D");
    reducedTree_HyperTriton -> Branch("chi2_ITS_Daughter1",&chi2_ITS_Daughter1,"chi2_ITS_Daughter1/D");
    reducedTree_HyperTriton -> Branch("nSigmaITS_He3_Daughter1",&nSigmaITS_He3_Daughter1,"nSigmaITS_He3_Daughter1/D");
    reducedTree_HyperTriton -> Branch("nSigmaTPC_He3_Daughter1",&nSigmaTPC_He3_Daughter1,"nSigmaTPC_He3_Daughter1/D");
    reducedTree_HyperTriton -> Branch("nSigmaTOF_He3_Daughter1",&nSigmaTOF_He3_Daughter1,"nSigmaTOF_He3_Daughter1/D");
    reducedTree_HyperTriton -> Branch("nSigmaITS_Pion_Daughter1",&nSigmaITS_Pion_Daughter1,"nSigmaITS_Pion_Daughter1/D");
    reducedTree_HyperTriton -> Branch("nSigmaTPC_Pion_Daughter1",&nSigmaTPC_Pion_Daughter1,"nSigmaTPC_Pion_Daughter1/D");
    reducedTree_HyperTriton -> Branch("nSigmaTOF_Pion_Daughter1",&nSigmaTOF_Pion_Daughter1,"nSigmaTOF_Pion_Daughter1/D");
    reducedTree_HyperTriton -> Branch("nSigmaITS_Trit_Daughter1",&nSigmaITS_Trit_Daughter1,"nSigmaITS_Trit_Daughter1/D");
    reducedTree_HyperTriton -> Branch("nSigmaTPC_Trit_Daughter1",&nSigmaTPC_Trit_Daughter1,"nSigmaTPC_Trit_Daughter1/D");
    reducedTree_HyperTriton -> Branch("nSigmaTOF_Trit_Daughter1",&nSigmaTOF_Trit_Daughter1,"nSigmaTOF_Trit_Daughter1/D");
    reducedTree_HyperTriton -> Branch("px_Daughter2",&px_Daughter2,"px_Daughter2/D");
    reducedTree_HyperTriton -> Branch("py_Daughter2",&py_Daughter2,"py_Daughter2/D");
    reducedTree_HyperTriton -> Branch("pz_Daughter2",&pz_Daughter2,"pz_Daughter2/D");
    reducedTree_HyperTriton -> Branch("TPCmomentum_Daughter2",&TPCmomentum_Daughter2,"TPCmomentum_Daughter2/D");
    reducedTree_HyperTriton -> Branch("trackID_Daughter2",&trackID_Daughter2,"trackID_Daughter2/I");
    reducedTree_HyperTriton -> Branch("eta_Daughter2",&eta_Daughter2,"eta_Daughter2/D");
    reducedTree_HyperTriton -> Branch("phi_Daughter2",&phi_Daughter2,"phi_Daughter2/D");
    reducedTree_HyperTriton -> Branch("theta_Daughter2",&theta_Daughter2,"theta_Daughter2/D");
    reducedTree_HyperTriton -> Branch("y_Daughter2",&y_Daughter2,"y_Daughter2/D");
    reducedTree_HyperTriton -> Branch("q_Daughter2",&q_Daughter2,"q_Daughter2/I");
    reducedTree_HyperTriton -> Branch("dcaxy_Daughter2",&dcaxy_Daughter2,"dcaxy_Daughter2/D");
    reducedTree_HyperTriton -> Branch("dcaz_Daughter2",&dcaz_Daughter2,"dcaz_Daughter2/D");
    reducedTree_HyperTriton -> Branch("nTPC_Clusters_Daughter2",&nTPC_Clusters_Daughter2,"nTPC_Clusters_Daughter2/I");
    reducedTree_HyperTriton -> Branch("nTPC_FindableClusters_Daughter2",&nTPC_FindableClusters_Daughter2,"nTPC_FindableClusters_Daughter2/I");
    reducedTree_HyperTriton -> Branch("nTPC_CrossedRows_Daughter2",&nTPC_CrossedRows_Daughter2,"nTPC_CrossedRows_Daughter2/I");
    reducedTree_HyperTriton -> Branch("nTPC_Clusters_dEdx_Daughter2",&nTPC_Clusters_dEdx_Daughter2,"nTPC_Clusters_dEdx_Daughter2/I");
    reducedTree_HyperTriton -> Branch("nITS_Clusters_Daughter2",&nITS_Clusters_Daughter2,"nITS_Clusters_Daughter2/I");
    reducedTree_HyperTriton -> Branch("HasPointOnITSLayer0_Daughter2",&HasPointOnITSLayer0_Daughter2,"HasPointOnITSLayer0_Daughter2/O");
    reducedTree_HyperTriton -> Branch("HasPointOnITSLayer1_Daughter2",&HasPointOnITSLayer1_Daughter2,"HasPointOnITSLayer1_Daughter2/O");
    reducedTree_HyperTriton -> Branch("HasPointOnITSLayer2_Daughter2",&HasPointOnITSLayer2_Daughter2,"HasPointOnITSLayer2_Daughter2/O");
    reducedTree_HyperTriton -> Branch("HasPointOnITSLayer3_Daughter2",&HasPointOnITSLayer3_Daughter2,"HasPointOnITSLayer3_Daughter2/O");
    reducedTree_HyperTriton -> Branch("HasPointOnITSLayer4_Daughter2",&HasPointOnITSLayer4_Daughter2,"HasPointOnITSLayer4_Daughter2/O");
    reducedTree_HyperTriton -> Branch("HasPointOnITSLayer5_Daughter2",&HasPointOnITSLayer5_Daughter2,"HasPointOnITSLayer5_Daughter2/O");
    reducedTree_HyperTriton -> Branch("HasSharedPointOnITSLayer0_Daughter2",&HasSharedPointOnITSLayer0_Daughter2,"HasSharedPointOnITSLayer0_Daughter2/O");
    reducedTree_HyperTriton -> Branch("HasSharedPointOnITSLayer1_Daughter2",&HasSharedPointOnITSLayer1_Daughter2,"HasSharedPointOnITSLayer1_Daughter2/O");
    reducedTree_HyperTriton -> Branch("HasSharedPointOnITSLayer2_Daughter2",&HasSharedPointOnITSLayer2_Daughter2,"HasSharedPointOnITSLayer2_Daughter2/O");
    reducedTree_HyperTriton -> Branch("HasSharedPointOnITSLayer3_Daughter2",&HasSharedPointOnITSLayer3_Daughter2,"HasSharedPointOnITSLayer3_Daughter2/O");
    reducedTree_HyperTriton -> Branch("HasSharedPointOnITSLayer4_Daughter2",&HasSharedPointOnITSLayer4_Daughter2,"HasSharedPointOnITSLayer4_Daughter2/O");
    reducedTree_HyperTriton -> Branch("HasSharedPointOnITSLayer5_Daughter2",&HasSharedPointOnITSLayer5_Daughter2,"HasSharedPointOnITSLayer5_Daughter2/O");
    reducedTree_HyperTriton -> Branch("chi2_TPC_Daughter2",&chi2_TPC_Daughter2,"chi2_TPC_Daughter2/D");
    reducedTree_HyperTriton -> Branch("chi2_NDF_Daughter2",&chi2_NDF_Daughter2,"chi2_NDF_Daughter2/D");
    reducedTree_HyperTriton -> Branch("chi2_ITS_Daughter2",&chi2_ITS_Daughter2,"chi2_ITS_Daughter2/D");
    reducedTree_HyperTriton -> Branch("nSigmaITS_He3_Daughter2",&nSigmaITS_He3_Daughter2,"nSigmaITS_He3_Daughter2/D");
    reducedTree_HyperTriton -> Branch("nSigmaTPC_He3_Daughter2",&nSigmaTPC_He3_Daughter2,"nSigmaTPC_He3_Daughter2/D");
    reducedTree_HyperTriton -> Branch("nSigmaTOF_He3_Daughter2",&nSigmaTOF_He3_Daughter2,"nSigmaTOF_He3_Daughter2/D");
    reducedTree_HyperTriton -> Branch("nSigmaITS_Pion_Daughter2",&nSigmaITS_Pion_Daughter2,"nSigmaITS_Pion_Daughter2/D");
    reducedTree_HyperTriton -> Branch("nSigmaTPC_Pion_Daughter2",&nSigmaTPC_Pion_Daughter2,"nSigmaTPC_Pion_Daughter2/D");
    reducedTree_HyperTriton -> Branch("nSigmaTOF_Pion_Daughter2",&nSigmaTOF_Pion_Daughter2,"nSigmaTOF_Pion_Daughter2/D");
    reducedTree_HyperTriton -> Branch("nSigmaITS_Trit_Daughter2",&nSigmaITS_Trit_Daughter2,"nSigmaITS_Trit_Daughter2/D");
    reducedTree_HyperTriton -> Branch("nSigmaTPC_Trit_Daughter2",&nSigmaTPC_Trit_Daughter2,"nSigmaTPC_Trit_Daughter2/D");
    reducedTree_HyperTriton -> Branch("nSigmaTOF_Trit_Daughter2",&nSigmaTOF_Trit_Daughter2,"nSigmaTOF_Trit_Daughter2/D");
    reducedTree_HyperTriton -> Branch("cosPointingAngle",&cosPointingAngle,"cosPointingAngle/D");
    reducedTree_HyperTriton -> Branch("dcaV0Daughters",&dcaV0Daughters,"dcaV0Daughters/D");
    reducedTree_HyperTriton -> Branch("radius",&radius,"radius/D");
    reducedTree_HyperTriton -> Branch("DecayLength",&DecayLength,"DecayLength/D");
    reducedTree_HyperTriton -> Branch("massV0",&massV0,"massV0/D");
    reducedTree_HyperTriton -> Branch("alphaV0",&alphaV0,"alphaV0/D");
    reducedTree_HyperTriton -> Branch("qtV0",&qtV0,"qtV0/D");
  }

  PostData(1, fQAList);
  PostData(2, fOutputList);
  PostData(3, reducedTree_Helium);
  if (fFillTri)    PostData(4,reducedTree_Triton);
  if (fFillHypTri) PostData(5,reducedTree_HyperTriton);
//     PostData(6, TreeEventSelection);

}
//_______________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskReducedTreeNuclei::UserExec(Option_t *)
{
  //Get Input Event
  if ( !GetInputEvent ()) {
//     TreeEventSelection -> Fill();
    return;
  }
//   TreeEventSelection -> Fill();

  //Load PID Response
  fPIDResponse = fInputHandler->GetPIDResponse();
  if(!fPIDResponse) {
    AliError("No PID Response found");
    return;
  }

  //Magnetic Field
  magFieldSign = 0;
  if(fAODevent->GetMagneticField() < 0) magFieldSign = -1;
  if(fAODevent->GetMagneticField() > 0) magFieldSign = 1;

  //Loop over Reconstructed Tracks
  for (Int_t i=0 ; i<fAODevent->GetNumberOfTracks() ; i++)  {

    //Track Selection
    AliAODTrack *track = (AliAODTrack*) fAODevent -> GetTrack(i);
    if ( !track ) continue;
    if ( !PassedBasicTrackQualityCuts (track)) continue;
    if ( !(IsHeliumCandidate (track) || IsTritonCandidate (track))) continue;
    // if ( fFillTri && !IsTritonCandidate(track)) continue;
    // if ( !fFillTri && !IsHeliumCandidate(track)) continue;

    px = track -> Px();
    py = track -> Py();
    pz = track -> Pz();

    q  = (Int_t) track -> Charge();
    eta = track -> Eta();
    phi = track -> Phi();
    theta = track -> Theta();
    y = track -> Y();
    TPCmomentum = track->GetTPCmomentum();
    TRDmomentum = track->GetTRDmomentum(1); //needs some Int_t for the plane? What shall we do?
    trackID = track->GetID();
    integratedLength = track->GetIntegratedLength();
    timeOfFlight = track->GetTOFsignal() - fPIDResponse->GetTOFResponse().GetStartTime(track->P());
    if (timeOfFlight>0)  {
      beta = integratedLength / (2.99792457999999984e-02 * timeOfFlight);
      Double_t InverseGamma2 = 1 - beta*beta;
      if (InverseGamma2 > 0) gamma = 1./TMath::Sqrt(InverseGamma2);
      if ( (gamma*gamma - 1) > 0 ) mass = track->P()/TMath::Sqrt(gamma*gamma - 1);
    } else{
      beta = 0;
      gamma = 0;
      mass = 0;
    }

    //DCA
    dcaz  = GetDCAz  (track);
    dcaxy = GetDCAxy (track);

    nTPC_Clusters = track->GetTPCNcls();
    nTRD_Clusters = track->GetTRDncls();
    nITS_Clusters = track->GetITSNcls();

    nTPC_FindableClusters = track->GetTPCNclsF();
    nTPC_CrossedRows = track->GetTPCNCrossedRows();
    nTPC_Clusters_dEdx = track -> GetTPCsignalN();

    HasPointOnITSLayer0 = track->HasPointOnITSLayer(0);
    HasPointOnITSLayer1 = track->HasPointOnITSLayer(1);
    HasPointOnITSLayer2 = track->HasPointOnITSLayer(2);
    HasPointOnITSLayer3 = track->HasPointOnITSLayer(3);
    HasPointOnITSLayer4 = track->HasPointOnITSLayer(4);
    HasPointOnITSLayer5 = track->HasPointOnITSLayer(5);

    HasSharedPointOnITSLayer0 = track->HasSharedPointOnITSLayer(0);
    HasSharedPointOnITSLayer1 = track->HasSharedPointOnITSLayer(1);
    HasSharedPointOnITSLayer2 = track->HasSharedPointOnITSLayer(2);
    HasSharedPointOnITSLayer3 = track->HasSharedPointOnITSLayer(3);
    HasSharedPointOnITSLayer4 = track->HasSharedPointOnITSLayer(4);
    HasSharedPointOnITSLayer5 = track->HasSharedPointOnITSLayer(5);

    chi2_TPC = track -> GetTPCchi2();//    check -> seems to be 0
    chi2_ITS = track -> GetITSchi2();//    check
    chi2_NDF = track -> Chi2perNDF();// chi2/NDF of momentum fit

    ITSsignal = track->GetITSsignal();
    TPCsignal = track->GetTPCsignal();
    TOFsignal = track->GetTOFsignal();
    TRDsignal = track->GetTRDsignal();
    HMPIDsignal = track->GetHMPIDsignal();

    nSigmaITS_He4 = fPIDResponse -> NumberOfSigmasITS(track,AliPID::kAlpha);
    nSigmaTPC_He4 = fPIDResponse -> NumberOfSigmasTPC(track,AliPID::kAlpha);
    nSigmaTOF_He4 = fPIDResponse -> NumberOfSigmasTOF(track,AliPID::kAlpha);
    nSigmaTRD_He4 = fPIDResponse -> NumberOfSigmasTRD(track,AliPID::kAlpha);
    nSigmaHMPID_He4 = fPIDResponse -> NumberOfSigmasHMPID(track,AliPID::kAlpha);

    nSigmaITS_He3 = fPIDResponse -> NumberOfSigmasITS(track,AliPID::kHe3);
    nSigmaTPC_He3 = fPIDResponse -> NumberOfSigmasTPC(track,AliPID::kHe3);
    nSigmaTOF_He3 = fPIDResponse -> NumberOfSigmasTOF(track,AliPID::kHe3);
    nSigmaTRD_He3 = fPIDResponse -> NumberOfSigmasTRD(track,AliPID::kHe3);
    nSigmaHMPID_He3 = fPIDResponse -> NumberOfSigmasHMPID(track,AliPID::kHe3);

    nSigmaITS_Trit = fPIDResponse -> NumberOfSigmasITS(track,AliPID::kTriton);
    nSigmaTPC_Trit = fPIDResponse -> NumberOfSigmasTPC(track,AliPID::kTriton);
    nSigmaTOF_Trit = fPIDResponse -> NumberOfSigmasTOF(track,AliPID::kTriton);
    nSigmaTRD_Trit = fPIDResponse -> NumberOfSigmasTRD(track,AliPID::kTriton);
    nSigmaHMPID_Trit = fPIDResponse -> NumberOfSigmasHMPID(track,AliPID::kTriton);

    nSigmaITS_Deut = fPIDResponse -> NumberOfSigmasITS(track,AliPID::kDeuteron);
    nSigmaTPC_Deut = fPIDResponse -> NumberOfSigmasTPC(track,AliPID::kDeuteron);
    nSigmaTOF_Deut = fPIDResponse -> NumberOfSigmasTOF(track,AliPID::kDeuteron);
    nSigmaTRD_Deut = fPIDResponse -> NumberOfSigmasTRD(track,AliPID::kDeuteron);
    nSigmaHMPID_Deut = fPIDResponse -> NumberOfSigmasHMPID(track,AliPID::kDeuteron);

    nSigmaITS_Prot = fPIDResponse -> NumberOfSigmasITS(track,AliPID::kProton);
    nSigmaTPC_Prot = fPIDResponse -> NumberOfSigmasTPC(track,AliPID::kProton);
    nSigmaTOF_Prot = fPIDResponse -> NumberOfSigmasTOF(track,AliPID::kProton);
    nSigmaTRD_Prot = fPIDResponse -> NumberOfSigmasTRD(track,AliPID::kProton);
    nSigmaHMPID_Prot = fPIDResponse -> NumberOfSigmasHMPID(track,AliPID::kProton);

    nSigmaITS_Pion = fPIDResponse -> NumberOfSigmasITS(track,AliPID::kPion);
    nSigmaTPC_Pion = fPIDResponse -> NumberOfSigmasTPC(track,AliPID::kPion);
    nSigmaTOF_Pion = fPIDResponse -> NumberOfSigmasTOF(track,AliPID::kPion);
    nSigmaTRD_Pion = fPIDResponse -> NumberOfSigmasTRD(track,AliPID::kPion);
    nSigmaHMPID_Pion = fPIDResponse -> NumberOfSigmasHMPID(track,AliPID::kPion);

    nSigmaITS_Kaon = fPIDResponse -> NumberOfSigmasITS(track,AliPID::kKaon);
    nSigmaTPC_Kaon = fPIDResponse -> NumberOfSigmasTPC(track,AliPID::kKaon);
    nSigmaTOF_Kaon = fPIDResponse -> NumberOfSigmasTOF(track,AliPID::kKaon);
    nSigmaTRD_Kaon = fPIDResponse -> NumberOfSigmasTRD(track,AliPID::kKaon);
    nSigmaHMPID_Kaon = fPIDResponse -> NumberOfSigmasHMPID(track,AliPID::kKaon);

    nSigmaITS_Elec = fPIDResponse -> NumberOfSigmasITS(track,AliPID::kElectron);
    nSigmaTPC_Elec = fPIDResponse -> NumberOfSigmasTPC(track,AliPID::kElectron);
    nSigmaTOF_Elec = fPIDResponse -> NumberOfSigmasTOF(track,AliPID::kElectron);
    nSigmaTRD_Elec = fPIDResponse -> NumberOfSigmasTRD(track,AliPID::kElectron);
    nSigmaHMPID_Elec = fPIDResponse -> NumberOfSigmasHMPID(track,AliPID::kElectron);

    //Fill the Tree
    if(IsHeliumCandidate(track)) reducedTree_Helium -> Fill();
    if(fFillTri && IsTritonCandidate(track)) reducedTree_Triton -> Fill();
  }


  if (fFillHypTri) {
    for ( Int_t iV0=0 ; iV0<fAODevent->GetNumberOfV0s() ; iV0++ ) {

      //Get V0 candidate
      AliAODv0 *V0 = (AliAODv0*)fAODevent->GetV0(iV0);
      if (!V0) continue;
      if (!V0->GetOnFlyStatus()) continue; //V0-on-the-fly
      AliAODTrack *posTrack=(AliAODTrack *)(V0->GetDaughter(0));
      AliAODTrack *negTrack=(AliAODTrack *)(V0->GetDaughter(1));

      if( !posTrack || !negTrack ) continue;

      if (!PassedBasicTrackQualityCuts(posTrack)) continue;
      if (!PassedBasicTrackQualityCuts(negTrack)) continue;
      if (!IsHyperTritonCandidate(posTrack,negTrack) && !IsHyperTritonCandidate(negTrack,posTrack)) continue;

      //Daughter1
      px_Daughter1 = posTrack -> Px();
      py_Daughter1 = posTrack -> Py();
      pz_Daughter1 = posTrack -> Pz();
      q_Daughter1  = (Int_t) posTrack -> Charge();
      eta_Daughter1 = posTrack -> Eta();
      phi_Daughter1 = posTrack -> Phi();
      theta_Daughter1 = posTrack -> Theta();
      y_Daughter1 = posTrack -> Y();
      TPCmomentum_Daughter1 = posTrack->GetTPCmomentum();
      trackID_Daughter1 = V0->GetPosID();


      //DCA
      dcaz_Daughter1  = GetDCAz  (posTrack);
      dcaxy_Daughter1 = GetDCAxy (posTrack);

      nTPC_Clusters_Daughter1 = posTrack->GetTPCNcls();
      nTRD_Clusters_Daughter1 = posTrack->GetTRDncls();
      nITS_Clusters_Daughter1 = posTrack->GetITSNcls();

      nTPC_FindableClusters_Daughter1 = posTrack->GetTPCNclsF();
      nTPC_CrossedRows_Daughter1 = posTrack->GetTPCNCrossedRows();
      nTPC_Clusters_dEdx_Daughter1 = posTrack -> GetTPCsignalN();

      HasPointOnITSLayer0_Daughter1 = posTrack->HasPointOnITSLayer(0);
      HasPointOnITSLayer1_Daughter1 = posTrack->HasPointOnITSLayer(1);
      HasPointOnITSLayer2_Daughter1 = posTrack->HasPointOnITSLayer(2);
      HasPointOnITSLayer3_Daughter1 = posTrack->HasPointOnITSLayer(3);
      HasPointOnITSLayer4_Daughter1 = posTrack->HasPointOnITSLayer(4);
      HasPointOnITSLayer5_Daughter1 = posTrack->HasPointOnITSLayer(5);

      HasSharedPointOnITSLayer0_Daughter1 = posTrack->HasSharedPointOnITSLayer(0);
      HasSharedPointOnITSLayer0_Daughter1 = posTrack->HasSharedPointOnITSLayer(1);
      HasSharedPointOnITSLayer0_Daughter1 = posTrack->HasSharedPointOnITSLayer(2);
      HasSharedPointOnITSLayer0_Daughter1 = posTrack->HasSharedPointOnITSLayer(3);
      HasSharedPointOnITSLayer0_Daughter1 = posTrack->HasSharedPointOnITSLayer(4);
      HasSharedPointOnITSLayer0_Daughter1 = posTrack->HasSharedPointOnITSLayer(5);

      chi2_TPC_Daughter1 = posTrack -> GetTPCchi2();
      chi2_ITS_Daughter1 = posTrack -> GetITSchi2();
      chi2_NDF_Daughter1 = posTrack -> Chi2perNDF();

      nSigmaITS_He3_Daughter1 = fPIDResponse -> NumberOfSigmasITS(posTrack,AliPID::kHe3);
      nSigmaTPC_He3_Daughter1 = fPIDResponse -> NumberOfSigmasTPC(posTrack,AliPID::kHe3);
      nSigmaTOF_He3_Daughter1 = fPIDResponse -> NumberOfSigmasTOF(posTrack,AliPID::kHe3);
      nSigmaITS_Pion_Daughter1 = fPIDResponse -> NumberOfSigmasITS(posTrack,AliPID::kPion);
      nSigmaTPC_Pion_Daughter1 = fPIDResponse -> NumberOfSigmasTPC(posTrack,AliPID::kPion);
      nSigmaTOF_Pion_Daughter1 = fPIDResponse -> NumberOfSigmasTOF(posTrack,AliPID::kPion);
      nSigmaITS_Trit_Daughter1 = fPIDResponse -> NumberOfSigmasITS(posTrack,AliPID::kTriton);
      nSigmaTPC_Trit_Daughter1 = fPIDResponse -> NumberOfSigmasTPC(posTrack,AliPID::kTriton);
      nSigmaTOF_Trit_Daughter1 = fPIDResponse -> NumberOfSigmasTOF(posTrack,AliPID::kTriton);

      //Daughter2
      px_Daughter2 = negTrack -> Px();
      py_Daughter2 = negTrack -> Py();
      pz_Daughter2 = negTrack -> Pz();
      q_Daughter2  = (Int_t) negTrack -> Charge();
      eta_Daughter2 = negTrack -> Eta();
      phi_Daughter2 = negTrack -> Phi();
      theta_Daughter2 = negTrack -> Theta();
      y_Daughter2 = negTrack -> Y();
      TPCmomentum_Daughter2 = negTrack->GetTPCmomentum();
      trackID_Daughter2 = V0->GetNegID();


      //DCA
      dcaz_Daughter2  = GetDCAz  (negTrack);
      dcaxy_Daughter2 = GetDCAxy (negTrack);

      nTPC_Clusters_Daughter2 = negTrack->GetTPCNcls();
      nTRD_Clusters_Daughter2 = negTrack->GetTRDncls();
      nITS_Clusters_Daughter2 = negTrack->GetITSNcls();

      nTPC_FindableClusters_Daughter2 = negTrack->GetTPCNclsF();
      nTPC_CrossedRows_Daughter2 = negTrack->GetTPCNCrossedRows();
      nTPC_Clusters_dEdx_Daughter2 = negTrack -> GetTPCsignalN();

      HasPointOnITSLayer0_Daughter2 = negTrack->HasPointOnITSLayer(0);
      HasPointOnITSLayer1_Daughter2 = negTrack->HasPointOnITSLayer(1);
      HasPointOnITSLayer2_Daughter2 = negTrack->HasPointOnITSLayer(2);
      HasPointOnITSLayer3_Daughter2 = negTrack->HasPointOnITSLayer(3);
      HasPointOnITSLayer4_Daughter2 = negTrack->HasPointOnITSLayer(4);
      HasPointOnITSLayer5_Daughter2 = negTrack->HasPointOnITSLayer(5);

      HasSharedPointOnITSLayer0_Daughter2 = negTrack->HasSharedPointOnITSLayer(0);
      HasSharedPointOnITSLayer0_Daughter2 = negTrack->HasSharedPointOnITSLayer(1);
      HasSharedPointOnITSLayer0_Daughter2 = negTrack->HasSharedPointOnITSLayer(2);
      HasSharedPointOnITSLayer0_Daughter2 = negTrack->HasSharedPointOnITSLayer(3);
      HasSharedPointOnITSLayer0_Daughter2 = negTrack->HasSharedPointOnITSLayer(4);
      HasSharedPointOnITSLayer0_Daughter2 = negTrack->HasSharedPointOnITSLayer(5);

      chi2_TPC_Daughter2 = negTrack -> GetTPCchi2();
      chi2_ITS_Daughter2 = negTrack -> GetITSchi2();
      chi2_NDF_Daughter2 = negTrack -> Chi2perNDF();

      nSigmaITS_He3_Daughter2 = fPIDResponse -> NumberOfSigmasITS(negTrack,AliPID::kHe3);
      nSigmaTPC_He3_Daughter2 = fPIDResponse -> NumberOfSigmasTPC(negTrack,AliPID::kHe3);
      nSigmaTOF_He3_Daughter2 = fPIDResponse -> NumberOfSigmasTOF(negTrack,AliPID::kHe3);
      nSigmaITS_Pion_Daughter2 = fPIDResponse -> NumberOfSigmasITS(negTrack,AliPID::kPion);
      nSigmaTPC_Pion_Daughter2 = fPIDResponse -> NumberOfSigmasTPC(negTrack,AliPID::kPion);
      nSigmaTOF_Pion_Daughter2 = fPIDResponse -> NumberOfSigmasTOF(negTrack,AliPID::kPion);
      nSigmaITS_Trit_Daughter2 = fPIDResponse -> NumberOfSigmasITS(negTrack,AliPID::kTriton);
      nSigmaTPC_Trit_Daughter2 = fPIDResponse -> NumberOfSigmasTPC(negTrack,AliPID::kTriton);
      nSigmaTOF_Trit_Daughter2 = fPIDResponse -> NumberOfSigmasTOF(negTrack,AliPID::kTriton);

      //Pair Variables
      AliAODVertex *vtxPrimary = fAODevent->GetPrimaryVertex();
      Double_t posVtx[3] = {0.,0.,0.};
      vtxPrimary->GetXYZ(posVtx);
      cosPointingAngle = V0->CosPointingAngle(posVtx);
      dcaV0Daughters  = V0->DcaV0Daughters();
      radius = V0->RadiusV0();
      DecayLength = V0->DecayLengthV0(posVtx);
      alphaV0 = V0 -> AlphaV0();
      qtV0 = V0 -> PtArmV0();

      reducedTree_HyperTriton -> Fill();
    }
  }

  PostData(1,fQAList);
  PostData(2,fOutputList);
  PostData(3,reducedTree_Helium);
  if (fFillTri)    PostData(4,reducedTree_Triton);
  if (fFillHypTri) PostData(5,reducedTree_HyperTriton);
//   PostData(6, TreeEventSelection);

}
//_______________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeNuclei::GetInputEvent ()  {

  //Get Input Event
  fAODevent = dynamic_cast <AliAODEvent*>(InputEvent());
  if (!fAODevent) return false;

  
  fAODVZERO = dynamic_cast <AliAODVZERO*>( fAODevent->GetVZEROData());
  
  // checking for event trigger
  Long64_t triggerMask = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  isTrigINT7        = (triggerMask & AliVEvent::kINT7)    ? true : false;
  isTrigHighMult    = (triggerMask & AliVEvent::kHighMultV0)? true : false;

  int whichTrig = 0;
  if (isTrigINT7 && isTrigHighMult){ whichTrig = 3;}
  else if (isTrigHighMult){ whichTrig = 2;}
  else if (isTrigINT7){ whichTrig = 1;}

  // triggerd events
  histoEventSelection->Fill(0.5,static_cast<Double_t>(whichTrig)); // Events before QA
  if (whichTrig==0) return false;
  
  
  // pileup rejection
//   fAODeventCuts.IsPileupFromSPD
  if (fUtils->IsPileUpEvent(fAODevent)) return false; //to be used in the analysis

//   if (!fUtils->IsOutOfBunchPileUp(fAODevent) return false; // do I need this?!
//   if (!fUtils->IsSPDClusterVsTrackletBG(fAODevent) return false; // do I need this?!
  
  
  
  //Event Cut
  histoEventSelection->Fill(1.5,static_cast<Double_t>(whichTrig)); // Events before QA

  fAODeventCuts.OverrideAutomaticTriggerSelection(AliVEvent::kAny,false);
  if (!fAODeventCuts.AcceptEvent(fAODevent)) {
    PostData(2, fOutputList);
    return false;
  }
  //Multiplicity Percentile
  AliMultSelection *MultSelection = NULL;
  MultSelection = (AliMultSelection*) fAODevent->FindListObject("MultSelection");

  histoEventSelection->Fill(2.5,static_cast<Double_t>(whichTrig)); // Events before multiplicity
  if (!MultSelection){
    PostData(2, fOutputList);
    return false;
  }

  multPercentile_V0M              = MultSelection->GetMultiplicityPercentile("V0M");
  multPercentile_V0A              = MultSelection->GetMultiplicityPercentile("V0A");
  multPercentile_V0C              = MultSelection->GetMultiplicityPercentile("V0C");
  multPercentile_OnlineV0M        = MultSelection->GetMultiplicityPercentile("OnlineV0M");
  multPercentile_OnlineV0A        = MultSelection->GetMultiplicityPercentile("OnlineV0A");
  multPercentile_OnlineV0C        = MultSelection->GetMultiplicityPercentile("OnlineV0C");
  multPercentile_ADM              = MultSelection->GetMultiplicityPercentile("ADM");
  multPercentile_ADA              = MultSelection->GetMultiplicityPercentile("ADA");
  multPercentile_ADC              = MultSelection->GetMultiplicityPercentile("ADC");
  multPercentile_SPDClusters      = MultSelection->GetMultiplicityPercentile("SPDClusters");
  multPercentile_SPDTracklets     = MultSelection->GetMultiplicityPercentile("SPDTracklets");
  multPercentile_RefMult05        = MultSelection->GetMultiplicityPercentile("RefMult05");
  multPercentile_RefMult08        = MultSelection->GetMultiplicityPercentile("RefMult08");
  multPercentile_CL1              = MultSelection->GetMultiplicityPercentile("CL1");
  multPercentile_ZNA              = MultSelection->GetMultiplicityPercentile("ZNA");

  //Multiplicity Estimator
  AliMultEstimator* estimator = NULL;
  estimator = MultSelection->GetEstimator("OnlineV0M");    if (estimator) { Ntrk_OnlineV0M = estimator->GetValue(); }
  estimator = MultSelection->GetEstimator("OnlineV0A");    if (estimator) { Ntrk_OnlineV0A = estimator->GetValue(); }
  estimator = MultSelection->GetEstimator("OnlineV0C");    if (estimator) { Ntrk_OnlineV0C = estimator->GetValue(); }
  estimator = MultSelection->GetEstimator("ADM");          if (estimator) { Ntrk_ADM = estimator->GetValue(); }
  estimator = MultSelection->GetEstimator("ADA");          if (estimator) { Ntrk_ADA = estimator->GetValue(); }
  estimator = MultSelection->GetEstimator("ADC");          if (estimator) { Ntrk_ADC = estimator->GetValue(); }
  estimator = MultSelection->GetEstimator("SPDClusters");  if (estimator) { Ntrk_SPDClusters = estimator->GetValue(); }
  estimator = MultSelection->GetEstimator("SPDTracklets"); if (estimator) { Ntrk_SPDTracklets = estimator->GetValue(); }
  estimator = MultSelection->GetEstimator("RefMult05");    if (estimator) { Ntrk_RefMult05 = estimator->GetValue(); }
  estimator = MultSelection->GetEstimator("RefMult08");    if (estimator) { Ntrk_RefMult08 = estimator->GetValue(); }
  estimator = MultSelection->GetEstimator("V0M");          if (estimator) { Ntrk_V0M = estimator->GetValue(); }
  estimator = MultSelection->GetEstimator("V0A");          if (estimator) { Ntrk_V0A = estimator->GetValue(); }
  estimator = MultSelection->GetEstimator("V0C");          if (estimator) { Ntrk_V0C = estimator->GetValue(); }

  if (isTrigINT7) histoEventMultV0M_MB->Fill(fAODVZERO->GetMTotV0C()+fAODVZERO->GetMTotV0A());
  else            histoEventMultV0M_HM->Fill(fAODVZERO->GetMTotV0C()+fAODVZERO->GetMTotV0A());


  histoEventSelection->Fill(3.5,static_cast<Double_t>(whichTrig)); // Selected events
  histoEventMultiplicity->Fill(multPercentile_V0M,static_cast<Double_t>(whichTrig));
  return true;
}
//_______________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeNuclei::PassedBasicTrackQualityCuts (AliAODTrack* track)  {

  //Filterbit
  if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return false;

  return true;
}
//_______________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskReducedTreeNuclei::GetDCAxy (AliAODTrack *track)  {

  Double_t impactParameter[2];
  Double_t covarianceMatrix[3];
  if (!track->PropagateToDCA(fAODevent->GetPrimaryVertex(),fAODevent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;

  Double_t DCAxy = impactParameter[0];

  return DCAxy;
}
//_______________________________________________________________________________________________________________________________________________________
Double_t AliAnalysisTaskReducedTreeNuclei::GetDCAz (AliAODTrack *track)  {

  Double_t impactParameter[2];
  Double_t covarianceMatrix[3];
  if (!track->PropagateToDCA(fAODevent->GetPrimaryVertex(),fAODevent->GetMagneticField(),10000,impactParameter,covarianceMatrix)) return -999;

  Double_t DCAz = impactParameter[1];

  return DCAz;
}
//_______________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeNuclei::IsHeliumCandidate (AliAODTrack *track)  {

  Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC(track,AliPID::kHe3);
  // if ( TMath::Abs(nsigmaTPC) > 6. ) return false;
  if ( nsigmaTPC < -6. ) return false; // used for 4He check


  return true;
}
//_______________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeNuclei::IsTritonCandidate (AliAODTrack *track)  {

  Double_t nsigmaTPC = fPIDResponse -> NumberOfSigmasTPC(track,AliPID::kTriton);
  // if ( nsigmaTPC < -5 ) return false;
  if (track->Pt() < 1.5) {
    if ( TMath::Abs(nsigmaTPC) > 6. ) return false;
  } else {
    Double_t nsigmaTOF = fPIDResponse -> NumberOfSigmasTOF(track,AliPID::kTriton);
    if ( TMath::Abs(nsigmaTOF) > 7. ) return false;
    if ( TMath::Abs(nsigmaTPC) > 6. ) return false;
  }

  return true;
}
//_______________________________________________________________________________________________________________________________________________________
Bool_t AliAnalysisTaskReducedTreeNuclei::IsHyperTritonCandidate (AliAODTrack *track1, AliAODTrack *track2)  {

  Double_t nsigmaTPC1 = fPIDResponse -> NumberOfSigmasTPC(track1,AliPID::kPion);
  Double_t nsigmaTPC2 = fPIDResponse -> NumberOfSigmasTPC(track2,AliPID::kHe3);

  if ( TMath::Abs(nsigmaTPC1) > 5 || TMath::Abs(nsigmaTPC2) > 5 ) return false;
  return true;
}
//_______________________________________________________________________________________________________________________________________________________
void AliAnalysisTaskReducedTreeNuclei::Terminate(Option_t *)  {
}
//_______________________________________________________________________________________________________________________________________________________
