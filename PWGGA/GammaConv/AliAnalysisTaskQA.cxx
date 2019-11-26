// based on AliAnalysisTaskConversionQA.cxx
// Related JIRA: PWGPP-527 ALIROOT-3736
// Contact: hikari.murakami@cern.ch

////////////////////////////////////////////////
//---------------------------------------------
// QA Task for V0 Reader V1
//---------------------------------------------
////////////////////////////////////////////////

#include "AliAnalysisTaskQA.h"
#include "TChain.h"
#include "TRandom.h"
#include "AliAnalysisManager.h"
#include "TParticle.h"
#include "TVectorF.h"
#include "AliPIDResponse.h"
#include "TFile.h"
#include "AliESDtrackCuts.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliAODEvent.h"
#include "AliMultSelection.h"

class iostream;

using namespace std;

ClassImp(AliAnalysisTaskQA)

//________________________________________________________________________
AliAnalysisTaskQA::AliAnalysisTaskQA() : AliAnalysisTaskSE(),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fConversionGammas(NULL),
  fConversionCuts(NULL),
  fEventCuts(NULL),
  fInputEvent(NULL),
  fNumberOfESDTracks(0),
  fMCEvent(NULL),
  fTreeQA(NULL),
  fIsHeavyIon(kFALSE),
//  ffillTree(-100),
  fTree(NULL),
  ffillHistograms(kFALSE),
  fOutputList(NULL),
  fESDList(NULL),
  hCentrality(NULL), 
  hBunch(NULL), 
  hVertexZ(NULL),
  hNGoodESDTracks(NULL),
  hNV0Tracks(NULL),
  hNContributorsVertex(NULL),
  hITSClusterPhi(NULL),
  hGammaPt(NULL),
  hGammaPhi(NULL),
  hGammaPhi_Pos(NULL),
  hGammaPhi_Neg(NULL),
  hGammaEta(NULL),
  hGammaChi2perNDF(NULL),
  hGammaPsiPair(NULL),
  hGammaArmenteros(NULL),
  hGammaCosinePointingAngle(NULL),
  hGammaInvMass(NULL),
  hElecPt(NULL),
  hElecEta(NULL),
  hElecPhi(NULL),
  hElecNfindableClsTPC(NULL),
  hPosiNfindableClsTPC(NULL),
  hElecClsTPC(NULL),
  hPosiClsTPC(NULL),
  hElectrondEdxP(NULL),
  hElectronITSdEdxP(NULL),
  hElectronTOFP(NULL),
  hElectronNSigmadEdxP(NULL),
  hElectronNSigmadEdxEta(NULL),
  hElectronNSigmaPiondEdxP(NULL),
  hElectronNSigmaITSP(NULL),
  hElectronNSigmaTOFP(NULL),
  hPositrondEdxP(NULL),
  hPositronITSdEdxP(NULL),
  hPositronTOFP(NULL),
  hPositronNSigmadEdxP(NULL),
  hPositronNSigmadEdxEta(NULL),
  hPositronNSigmaPiondEdxP(NULL),
  hPositronNSigmaITSP(NULL),
  hPositronNSigmaTOFP(NULL),
  hInvMassPair(NULL),
  //  hElecAsymP(NULL),
  //  fTrueList(NULL),
  //  hTrueResolutionR(NULL),
  //  hTrueResolutionZ(NULL),
  //  hTrueResolutionPhi(NULL),
  //  hTrueGammaPt(NULL),
  //  hTrueGammaPhi(NULL),
  //  hTrueGammaEta(NULL),
  //  hTrueGammaMass(NULL),
  //  hTrueGammaChi2perNDF(NULL),
  //  hTrueGammaPsiPair(NULL),
  //  hTrueGammaQt(NULL),
  //  hTrueGammaCosinePointingAngle(NULL),
  //  hTrueGammaXY(NULL),
  //  hTrueGammaZR(NULL),
  //  hTrueElecPt(NULL),
  //  hTrueElecEta(NULL),
  //  hTrueElecPhi(NULL),
  //  hTrueElecNfindableClsTPC(NULL),
  //  hTruePosiNfindableClsTPC(NULL),
  //  hTrueElecAsymP(NULL),
  fKind(0),
  fIsMC(kFALSE),
  fnGammaCandidates(1),
  fMCStackPos(NULL),
  fMCStackNeg(NULL)
  ,fActiveBranches("")
  ,fInactiveBranches("")
  ,fWriteVariableTree(kFALSE)
  ,Centrality(0)
  ,VertexZ(0)
  ,Bunch(0)
  ,GoodESDTracks(0)
  ,theta(0)
  ,pt(0)   
  ,phi(0)
  ,chi2(0)
  ,qt(0)
  ,alpha(0)
  ,psipair(0)
  ,cosPA(0)
  ,InvMass(0)
  ,X(0)
  ,Y(0)
  ,Z(0)
  ,R(0)
  ,Qual(0)
  ,DCAz(0)
  ,DCAr(0)
  ,ele_theta(0)
  ,ele_pt(0)
  ,ele_phi(0)
  ,ele_nSigmaTPC(0)
  ,ele_nSigmaTPCpion(0)
  ,ele_nSigmaTOF(0)
  ,ele_nSigmaITS(0)
  ,ele_TPCsignal(0)
  ,ele_TOFsignal(0)
  ,ele_ITSsignal(0)
  ,ele_Cls(0)
  ,ele_NfindableCls(0)
  ,ele_SPD1(0)
  ,ele_SPD2(0)
  ,ele_SDD1(0)
  ,ele_SDD2(0)
  ,ele_SSD1(0)
  ,ele_SSD2(0)
  ,pos_theta(0)
  ,pos_pt(0)
  ,pos_phi(0)
  ,pos_nSigmaTPC(0)
  ,pos_nSigmaTPCpion(0)
  ,pos_nSigmaTOF(0)
  ,pos_nSigmaITS(0)
  ,pos_TPCsignal(0)
  ,pos_TOFsignal(0)
  ,pos_ITSsignal(0)
  ,pos_Cls(0)
  ,pos_NfindableCls(0)
  ,pos_SPD1(0)
  ,pos_SPD2(0)
  ,pos_SDD1(0)
  ,pos_SDD2(0)
  ,pos_SSD1(0)
  ,pos_SSD2(0)

{

}

AliAnalysisTaskQA::AliAnalysisTaskQA(const char *name) : AliAnalysisTaskSE(name),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fConversionGammas(NULL),
  fConversionCuts(NULL),
  fEventCuts(NULL),
  fInputEvent(NULL),
  fNumberOfESDTracks(0),
  fMCEvent(NULL),
  fTreeQA(NULL),
  fIsHeavyIon(kFALSE),
							 //  ffillTree(-100),
  fTree(NULL),
  ffillHistograms(kFALSE),
  fOutputList(NULL),
  fESDList(NULL),
  hCentrality(NULL),
  hBunch(NULL),
  hVertexZ(NULL),
  hNGoodESDTracks(NULL),
  hNV0Tracks(NULL),
  hNContributorsVertex(NULL),
  hITSClusterPhi(NULL),
  hGammaPt(NULL),
  hGammaPhi(NULL),
  hGammaPhi_Pos(NULL),
  hGammaPhi_Neg(NULL),
  hGammaEta(NULL),
  hGammaChi2perNDF(NULL),
  hGammaPsiPair(NULL),
  hGammaArmenteros(NULL),
  hGammaCosinePointingAngle(NULL),
  hGammaInvMass(NULL),
  hElecPt(NULL),
  hElecEta(NULL),
  hElecPhi(NULL),
  hElecNfindableClsTPC(NULL),
  hPosiNfindableClsTPC(NULL),
  hElecClsTPC(NULL),
  hPosiClsTPC(NULL),
  hElectrondEdxP(NULL),
  hElectronITSdEdxP(NULL),
  hElectronTOFP(NULL),
  hElectronNSigmadEdxP(NULL),
  hElectronNSigmadEdxEta(NULL),
  hElectronNSigmaPiondEdxP(NULL),
  hElectronNSigmaITSP(NULL),
  hElectronNSigmaTOFP(NULL),
  hPositrondEdxP(NULL),
  hPositronITSdEdxP(NULL),
  hPositronTOFP(NULL),
  hPositronNSigmadEdxP(NULL),
  hPositronNSigmadEdxEta(NULL),
  hPositronNSigmaPiondEdxP(NULL),
  hPositronNSigmaITSP(NULL),
  hPositronNSigmaTOFP(NULL),
  hInvMassPair(NULL),
  //  hGammaXY(NULL),
  //  hGammaZR(NULL),
  //  hElecAsymP(NULL),
  //  fTrueList(NULL),
  //  hTrueResolutionR(NULL),
  //  hTrueResolutionZ(NULL),
  //  hTrueResolutionPhi(NULL),
  //  hTrueGammaPt(NULL),
  //  hTrueGammaPhi(NULL),
  //  hTrueGammaEta(NULL),
  //  hTrueGammaMass(NULL),
  //  hTrueGammaChi2perNDF(NULL),
  //  hTrueGammaPsiPair(NULL),
  //  hTrueGammaQt(NULL),
  //  hTrueGammaCosinePointingAngle(NULL),
  //  hTrueGammaXY(NULL),
  //  hTrueGammaZR(NULL),
  //  hTrueElecPt(NULL),
  //  hTrueElecEta(NULL),
  //  hTrueElecPhi(NULL),
  //  hTrueElecNfindableClsTPC(NULL),
  //  hTruePosiNfindableClsTPC(NULL),
  //  hTrueElecAsymP(NULL),
  fKind(0),
  fIsMC(kFALSE),
  fnGammaCandidates(1),
  fMCStackPos(NULL),
  fMCStackNeg(NULL)
  ,fActiveBranches("")
  ,fInactiveBranches("")
  ,fWriteVariableTree(kFALSE)
  ,Centrality(0)
  ,VertexZ(0)
  ,Bunch(0)
  ,GoodESDTracks(0)
  ,theta(0)
  ,pt(0)
  ,phi(0)
  ,chi2(0)
  ,qt(0)
  ,alpha(0)
  ,psipair(0)
  ,cosPA(0)
  ,InvMass(0)
  ,X(0)
  ,Y(0)
  ,Z(0)
  ,R(0)
  ,Qual(0)
  ,DCAz(0)
  ,DCAr(0)									     
  ,ele_theta(0)
  ,ele_pt(0)
  ,ele_phi(0)
  ,ele_nSigmaTPC(0)
  ,ele_nSigmaTPCpion(0)
  ,ele_nSigmaTOF(0)
  ,ele_nSigmaITS(0)
  ,ele_TPCsignal(0)
  ,ele_TOFsignal(0)
  ,ele_ITSsignal(0)
  ,ele_Cls(0)
  ,ele_NfindableCls(0)
  ,ele_SPD1(0)
  ,ele_SPD2(0)
  ,ele_SDD1(0)
  ,ele_SDD2(0)
  ,ele_SSD1(0)
  ,ele_SSD2(0)
  ,pos_theta(0)
  ,pos_pt(0)
  ,pos_phi(0)
  ,pos_nSigmaTPC(0)
  ,pos_nSigmaTPCpion(0)
  ,pos_nSigmaTOF(0)
  ,pos_nSigmaITS(0)
  ,pos_TPCsignal(0)
  ,pos_TOFsignal(0)
  ,pos_ITSsignal(0)
  ,pos_Cls(0)
  ,pos_NfindableCls(0)
  ,pos_SPD1(0)
  ,pos_SPD2(0)
  ,pos_SDD1(0)
  ,pos_SDD2(0)
  ,pos_SSD1(0)
  ,pos_SSD2(0)

{
  // Default constructor

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskQA::~AliAnalysisTaskQA()
{
  // default deconstructor
  
}
//________________________________________________________________________
void AliAnalysisTaskQA::UserCreateOutputObjects()
{
  // Create User Output Objects

  if(fOutputList != NULL){
    delete fOutputList;
    fOutputList = NULL;
  }
  if(fOutputList == NULL){
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
  }
  
  if(ffillHistograms){
    // fESDList = new TList();
    // fESDList->SetOwner(kTRUE);
    // fESDList->SetName("ESD QA");
    // fOutputList->Add(fESDList);

    hCentrality = new TH1F("Centrality","",100,0,100);
    fOutputList->Add(hCentrality);
    hBunch = new TH1I("Bunch","",4000,0,4000);
    fOutputList->Add(hBunch);
    hVertexZ = new TH1F("Vertex_Z","",300,-15,15);
    fOutputList->Add(hVertexZ);
    hNContributorsVertex = new TH1I("ContrVertex_Z","",3000,0,3000);
    fOutputList->Add(hNContributorsVertex);
    //    if(fIsHeavyIon) hNGoodESDTracks = new TH1I("GoodESDTracks","",4000,0,4000);
    //else
    hNGoodESDTracks = new TH1I("GoodESDTracks","",200,0,200);
    fOutputList->Add(hNGoodESDTracks);
    //if(fIsHeavyIon) hNV0Tracks = new TH1I("V0 Multiplicity","",30000,0,30000);
    //else
    hNV0Tracks = new TH1I("V0 Multiplicity","",2000,0,2000);
    fOutputList->Add(hNV0Tracks);

    hITSClusterPhi = new TH2F("ITSClusterPhi","",72,0,2*TMath::Pi(),7,-0.5,6.5);
    fOutputList->Add(hITSClusterPhi);
    hGammaPt = new TH1F("Gamma_Pt","",250,0,25);
    fOutputList->Add(hGammaPt);
    hGammaPhi = new TH1F("Gamma_Phi","",360,0,2*TMath::Pi());
    fOutputList->Add(hGammaPhi);
    hGammaPhi_Pos = new TH1F("GammaPhi_EtaPos","",360,0,2*TMath::Pi());
    fOutputList->Add(hGammaPhi_Pos);
    hGammaPhi_Neg = new TH1F("GammaPhi_EtaNeg","",360,0,2*TMath::Pi());
    fOutputList->Add(hGammaPhi_Neg);
  
    hGammaEta = new TH1F("Gamma_Eta","",600,-1.5,1.5);
    fOutputList->Add(hGammaEta);
    hGammaChi2perNDF = new TH1F("Gamma_Chi2perNDF","",500,0,50);
    fOutputList->Add(hGammaChi2perNDF);
    hGammaPsiPair = new TH1F("Gamma_PsiPair","",500,0,2);
    fOutputList->Add(hGammaPsiPair);
    hGammaArmenteros = new TH2F("Gamma_Armenteros","",200,-1,1,400,0,0.1);
    fOutputList->Add(hGammaArmenteros);
    hGammaCosinePointingAngle = new TH1F("Gamma_CosinePointingAngle","",300,0.85,1.);
    fOutputList->Add(hGammaCosinePointingAngle);
    hGammaInvMass = new TH1F( "Gamma_InvMass","",200, 0, 0.2);
    fOutputList->Add(hGammaInvMass);
    
    hElecPt = new TH2F("Electron_Positron_Pt","",250,0,25,250,0,25);
    fOutputList->Add(hElecPt);
    hElecEta = new TH2F("Electron_Positron_Eta","",600,-1.5,1.5,600,-1.5,1.5);
    fOutputList->Add(hElecEta);
    hElecPhi = new TH2F("Electron_Positron_Phi","",360,0,2*TMath::Pi(),360,0,2*TMath::Pi());
    fOutputList->Add(hElecPhi);
    hElecClsTPC = new TH1F("Electron_ClusterTPC","",200,0,200);
    fOutputList->Add(hElecClsTPC);
    hPosiClsTPC = new TH1F("Positron_ClusterTPC","",200,0,200);
    fOutputList->Add(hPosiClsTPC);
    
    hElecNfindableClsTPC = new TH1F("Electron_findableClusterTPC","",100,0,1);
    fOutputList->Add(hElecNfindableClsTPC);
    hPosiNfindableClsTPC = new TH1F("Positron_findableClusterTPC","",100,0,1);
    fOutputList->Add(hPosiNfindableClsTPC);
    
    hElectrondEdxP =  new TH2F("Electron_dEdx_P","",100, 0.05, 20, 200, 0, 200);
    SetLogBinningXTH2(hElectrondEdxP);
    fOutputList->Add(hElectrondEdxP);
    hPositrondEdxP =  new TH2F("Positron_dEdx_P","",100, 0.05, 20, 200, 0, 200);
    SetLogBinningXTH2(hPositrondEdxP);
    fOutputList->Add(hPositrondEdxP);
    hElectronNSigmadEdxP =  new TH2F("Electron_NSigmadEdx_P","",100, 0.05, 20, 200, -10, 10);  
    SetLogBinningXTH2(hElectronNSigmadEdxP);
    fOutputList->Add(hElectronNSigmadEdxP);
    hElectronNSigmadEdxEta =  new TH2F("Electron_NSigmadEdx_Eta","",140, -1.4, 1.4, 200, -10, 10);  
    fOutputList->Add(hElectronNSigmadEdxEta);
    hPositronNSigmadEdxP =  new TH2F("Positron_NSigmadEdx_P","",100, 0.05, 20, 200, -10, 10);
    SetLogBinningXTH2(hPositronNSigmadEdxP);
    fOutputList->Add(hPositronNSigmadEdxP);
    hPositronNSigmadEdxEta =  new TH2F("Positron_NSigmadEdx_Eta","",140, -1.4, 1.4, 200, -10, 10);  
    fOutputList->Add(hPositronNSigmadEdxEta);
    hElectronNSigmaPiondEdxP =  new TH2F("Electron_NSigmaPiondEdx_P","",100, 0.05, 20, 200, -10, 10);  
    SetLogBinningXTH2(hElectronNSigmaPiondEdxP);
    fOutputList->Add(hElectronNSigmaPiondEdxP);
    hPositronNSigmaPiondEdxP =  new TH2F("Positron_NSigmaPiondEdx_P","",100, 0.05, 20, 200, -10, 10);
    SetLogBinningXTH2(hPositronNSigmaPiondEdxP);
    fOutputList->Add(hPositronNSigmaPiondEdxP);
    
    hElectronTOFP =  new TH2F("Electron_TOF_P","",100, 0.05, 20, 600, -1000, 29000);
    SetLogBinningXTH2(hElectronTOFP);
    fOutputList->Add(hElectronTOFP);
    hPositronTOFP =  new TH2F("Positron_TOF_P","",100, 0.05, 20, 600, -1000, 29000);
    SetLogBinningXTH2(hPositronTOFP);
    fOutputList->Add(hPositronTOFP);
    hElectronNSigmaTOFP =  new TH2F("Electron_NSigmaTOF_P","",100, 0.05, 20, 200, -10, 10);
    SetLogBinningXTH2(hElectronNSigmaTOFP);
    fOutputList->Add(hElectronNSigmaTOFP);
    hPositronNSigmaTOFP =  new TH2F("Positron_NSigmaTOF_P","",100, 0.05, 20, 200, -10, 10);
    SetLogBinningXTH2(hPositronNSigmaTOFP);
    fOutputList->Add(hPositronNSigmaTOFP);
    
    hElectronITSdEdxP  =  new TH2F("Electron_ITSdEdx_P","",100, 0.05, 20, 200, 0, 200);
    SetLogBinningXTH2(hElectronITSdEdxP);
    fOutputList->Add(hElectronITSdEdxP);
    hPositronITSdEdxP =  new TH2F("Positron_ITSdEdx_P","",100, 0.05, 20, 200, 0, 200);
    SetLogBinningXTH2(hPositronITSdEdxP);
    fOutputList->Add(hPositronITSdEdxP);
    hElectronNSigmaITSP =  new TH2F("Electron_NSigmaITS_P","",100, 0.05, 20, 200, -10, 10);
    SetLogBinningXTH2(hElectronNSigmaITSP);
    fOutputList->Add(hElectronNSigmaITSP);
    hPositronNSigmaITSP =  new TH2F("Positron_NSigmaITS_P","",100, 0.05, 20, 200, -10, 10);
    SetLogBinningXTH2(hPositronNSigmaITSP);
    fOutputList->Add(hPositronNSigmaITSP);
    
    hInvMassPair= new TH2F("Gamma_InvMassPair_Pt","",200,0,0.2,250,0,25);
    hInvMassPair->SetTitle(";p_{T} (GeV/c);mass (GeV/c)");
    fOutputList->Add(hInvMassPair);
    
    //     hGammaXY = new TH2F("Gamma_ConversionPoint_XY","Gamma_ConversionPoint_XY",960,-120,120,960,-120,120);
    //     fOutputList-->Add(hGammaXY);
    //     hGammaZR= new TH2F("Gamma_ConversionPoint_ZR","Gamma_ConversionPoint_ZR",1200,-150,150,480,0,120);
    //     fOutputList-->Add(hGammaZR);


    //     hElecAsymP = new TH2F("Electron_Asym_vs_P", "Electron_Asym_vs_P",200,0.,20.,200,0.,1.); 
    //     fOutputList-->Add(hElecAsymP);

    //     if(fIsMC){
    //      fTrueList = new TList();
    //      fTrueList->SetOwner(kTRUE);
    //      fTrueList->SetName("True QA");
    //      fOutputList->Add(fTrueList);
    // 
    //      hTrueResolutionR = new TH2F("True_ConversionPointResolution_R","True_ConversionPointResolution_R",240,0,120,200,-20,20);
    //      fTrueList->Add(hTrueResolutionR);
    //      hTrueResolutionZ = new TH2F("True_ConversionPointResolution_Z","True_ConversionPointResolution_Z",480,-120,120,200,-20,20);
    //      fTrueList->Add(hTrueResolutionZ);
    //      hTrueResolutionPhi = new TH2F("True_ConversionPointResolution_Phi","True_ConversionPointResolution_Phi",360,0,2*TMath::Pi(),200,-TMath::Pi()/30., TMath::Pi()/30.);
    //      fTrueList->Add(hTrueResolutionPhi);
    // 
    //      hTrueGammaPt = new TH1F("True_Gamma_Pt","True_Gamma_Pt",250,0,25);
    //      fTrueList->Add(hTrueGammaPt);
    //      hTrueGammaPhi = new TH1F("True_Gamma_Phi","True_Gamma_Phi",360,0,2*TMath::Pi());
    //      fTrueList->Add(hTrueGammaPhi);
    //      hTrueGammaEta = new TH1F("True_Gamma_Eta","True_Gamma_Eta",600,-1.5,1.5);
    //      fTrueList->Add(hTrueGammaEta);
    //      hTrueGammaMass = new TH1F("True_Gamma_Mass","True_Gamma_Mass",1000,0,0.3);
    //      fTrueList->Add(hTrueGammaMass);
    //      hTrueGammaChi2perNDF = new TH1F("True_Gamma_Chi2perNDF","True_Gamma_Chi2perNDF",500,0,100);
    //      fTrueList->Add(hTrueGammaChi2perNDF);
    //      hTrueGammaPsiPair = new TH1F("True_Gamma_PsiPair","True_Gamma_PsiPair",500,0,2);
    //      fTrueList->Add(hTrueGammaPsiPair);
    //      hTrueGammaQt = new TH1F("True_Gamma_Qt","True_Gamma_Qt",400,0,0.1);
    //      fTrueList->Add(hTrueGammaQt);
    //      hTrueGammaCosinePointingAngle = new TH1F("True_Gamma_CosinePointingAngle","True_Gamma_CosinePointingAngle",900,0.7,1.);
    //      fTrueList->Add(hTrueGammaCosinePointingAngle);
    //      hTrueGammaXY = new TH2F("True_Gamma_ConversionPoint_XY","True_Gamma_ConversionPoint_XY",960,-120,120,960,-120,120);
    //      fTrueList->Add(hTrueGammaXY);
    //      hTrueGammaZR= new TH2F("TrueGamma_ConversionPoint_ZR","TrueGamma_ConversionPoint_ZR",1200,-150,150,480,0,120);
    //      fTrueList->Add(hTrueGammaZR);
    // 
    //      hTrueElecPt = new TH2F("True_Electron_Positron_Pt","True_Electron_Positron_Pt",250,0,25,250,0,25);
    //      fTrueList->Add(hTrueElecPt);
    //      hTrueElecEta = new TH2F("True_Electron_Positron_Eta","True_Electron_Positron_Eta",600,-1.5,1.5,600,-1.5,1.5);
    //      fTrueList->Add(hTrueElecEta);
    //      hTrueElecPhi = new TH2F("True_Electron_Positron_Phi","True_Electron_Positron_Phi",360,0,2*TMath::Pi(),360,0,2*TMath::Pi());
    //      fTrueList->Add(hTrueElecPhi);
    //      hTrueElecNfindableClsTPC = new TH1F("True_Electron_findableClusterTPC","True_Electron_findableClusterTPC",100,0,1);
    //      fTrueList->Add(hTrueElecNfindableClsTPC);
    //      hTruePosiNfindableClsTPC = new TH1F("True_Positron_findableClusterTPC","True_Positron_findableClusterTPC",100,0,1);
    //      fTrueList->Add(hTruePosiNfindableClsTPC);
    // 				 hTrueElecAsymP = new TH2F("True_Electron_Asym_vs_P", "True_Electron_Asym_vs_P",200,0.,20.,200,0.,1.); 
    // 				 fTrueList->Add(hTrueElecAsymP);
    //}
    if(fConversionCuts->GetCutHistograms()){
      fOutputList->Add(fConversionCuts->GetCutHistograms());
    }
  }
  TString suffix ="";
  if(fWriteVariableTree)  suffix = "tree";
  // Create User Output Objects
  fTree = new TTree(Form("Gamma_%s_%s_%s",suffix.Data(),(fEventCuts->GetCutNumber()).Data(),(fConversionCuts->GetCutNumber()).Data()),"GammaProp");      

  fTree->Branch("Centrality",         &Centrality);
  fTree->Branch("VertexZ",            &VertexZ);
  fTree->Branch("Bunch",              &Bunch);
  fTree->Branch("GoodESDTracks",      &GoodESDTracks);

  fTree->Branch("theta",              &theta);
  fTree->Branch("pt",                 &pt);
  fTree->Branch("phi",                &phi);
  fTree->Branch("chi2",               &chi2);
  fTree->Branch("qt",                 &qt);
  fTree->Branch("alpha",              &alpha);
  fTree->Branch("psipair",            &psipair);
  fTree->Branch("cosPA",              &cosPA);
  fTree->Branch("InvMass",            &InvMass);
  fTree->Branch("X",                  &X);
  fTree->Branch("Y",                  &Y);
  fTree->Branch("Z",                  &Z);
  fTree->Branch("R",                  &R);
  fTree->Branch("Qual",               &Qual);
  fTree->Branch("DCAz",               &DCAz);
  fTree->Branch("DCAr",               &DCAr);

  fTree->Branch("ele_theta",          &ele_theta);
  fTree->Branch("ele_pt",             &ele_pt);
  fTree->Branch("ele_phi",            &ele_phi);
  fTree->Branch("ele_nSigmaTPC",      &ele_nSigmaTPC);
  fTree->Branch("ele_nSigmaTPCpion",  &ele_nSigmaTPCpion);
  fTree->Branch("ele_nSigmaTOF",      &ele_nSigmaTOF);
  fTree->Branch("ele_nSigmaITS",      &ele_nSigmaITS);
  fTree->Branch("ele_TPCsignal",      &ele_TPCsignal);
  fTree->Branch("ele_TOFsignal",      &ele_TOFsignal);
  fTree->Branch("ele_ITSsignal",      &ele_ITSsignal);
  fTree->Branch("ele_Cls",            &ele_Cls);
  fTree->Branch("ele_SPD1",           &ele_SPD1);
  fTree->Branch("ele_SPD2",           &ele_SPD2);
  fTree->Branch("ele_SDD1",           &ele_SDD1);
  fTree->Branch("ele_SDD2",           &ele_SDD2);
  fTree->Branch("ele_SSD1",           &ele_SSD1);
  fTree->Branch("ele_SSD2",           &ele_SSD2);
          
  fTree->Branch("pos_theta",          &pos_theta);
  fTree->Branch("pos_pt",             &pos_pt);
  fTree->Branch("pos_phi",            &pos_phi);
  fTree->Branch("pos_nSigmaTPC",      &pos_nSigmaTPC);
  fTree->Branch("pos_nSigmaTPCpion",  &pos_nSigmaTPCpion);
  fTree->Branch("pos_nSigmaTOF",      &pos_nSigmaTOF);
  fTree->Branch("pos_nSigmaITS",      &pos_nSigmaITS);
  fTree->Branch("pos_TPCsignal",      &pos_TPCsignal);
  fTree->Branch("pos_TOFsignal",      &pos_TOFsignal);
  fTree->Branch("pos_ITSsignal",      &pos_ITSsignal);
  fTree->Branch("pos_Cls",            &pos_Cls);
  fTree->Branch("pos_SPD1",           &pos_SPD1);
  fTree->Branch("pos_SPD2",           &pos_SPD2);
  fTree->Branch("pos_SDD1",           &pos_SDD1);
  fTree->Branch("pos_SDD2",           &pos_SDD2);
  fTree->Branch("pos_SSD1",           &pos_SSD1);
  fTree->Branch("pos_SSD2",           &pos_SSD2);

  
  // if user set active branches
  TObjArray* aractive=fActiveBranches.Tokenize(";");
  if(aractive->GetEntries()>0) {fTree->SetBranchStatus("*", 0);}
  for(Int_t i=0; i<aractive->GetEntries(); i++){
    fTree->SetBranchStatus(aractive->At(i)->GetName(), 1);
  }

  // if user set inactive branches
  TObjArray* arinactive=fInactiveBranches.Tokenize(";");
  for(Int_t i=0; i<arinactive->GetEntries(); i++){
    fTree->SetBranchStatus(arinactive->At(i)->GetName(), 0);
  }
  
  if (fIsMC) {
    fTree->Branch("kind",&fKind,"fKind/b");
  }   
  
  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
  
  if(fV0Reader && fV0Reader->GetProduceV0FindingEfficiency())
    if (fV0Reader->GetV0FindingEfficiencyHistograms())
      fOutputList->Add(fV0Reader->GetV0FindingEfficiencyHistograms());

  PostData(1, fOutputList);
  PostData(2, fTree);
  //  if(ffillTree>=1.0){
  //    OpenFile(2);
  //    PostData(2, fTreeQA);
  //  }
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskQA::Notify()
{
    if (fEventCuts->GetPeriodEnum() == AliConvEventCuts::kNoPeriod && ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum() != AliConvEventCuts::kNoPeriod){        
        fEventCuts->SetPeriodEnumExplicit(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum());
    } else if (fEventCuts->GetPeriodEnum() == AliConvEventCuts::kNoPeriod ){
      fEventCuts->SetPeriodEnum(fV0Reader->GetPeriodName());
    }  
 
  
  if(!fEventCuts->GetDoEtaShift()) return kTRUE; // No Eta Shift requested, continue
    
  if(fEventCuts->GetEtaShift() == 0.0){ // Eta Shift requested but not set, get shift automatically
    fEventCuts->GetCorrectEtaShiftFromPeriod();
    fEventCuts->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
    return kTRUE;
  }
  else{
    printf(" Gamma Conversion QA Task %s :: Eta Shift Manually Set to %f \n\n",
        (fEventCuts->GetCutNumber()).Data(),fEventCuts->GetEtaShift());
    fEventCuts->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
  }
  
  return kTRUE;
}
//________________________________________________________________________
void AliAnalysisTaskQA::UserExec(Option_t *){

  Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
  if(eventQuality != 0){// Event Not Accepted
    return;
  }
  fInputEvent = InputEvent();
  if(fIsMC) fMCEvent = MCEvent();

  Int_t eventNotAccepted =
    fEventCuts->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon,kFALSE);
  if(eventNotAccepted) return; // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1

  AliMultSelection *MultSelection = (AliMultSelection*)fInputEvent->FindListObject("MultSelection");
  Double_t cent=MultSelection->GetMultiplicityPercentile("V0A");
  fConversionGammas=fV0Reader->GetReconstructedGammas();

  if(fMCEvent){
    if(fEventCuts->GetSignalRejection() != 0){
      if(fInputEvent->IsA()==AliESDEvent::Class()){
        fEventCuts->GetNotRejectedParticles(fEventCuts->GetSignalRejection(),
					    fEventCuts->GetAcceptedHeader(),
					    fMCEvent);
      }
      else if(fInputEvent->IsA()==AliAODEvent::Class()){
        fEventCuts->GetNotRejectedParticles(fEventCuts->GetSignalRejection(),
					    fEventCuts->GetAcceptedHeader(),
					    fInputEvent);
      }
    }
  }

  if(ffillHistograms){
    hCentrality->Fill(cent);
    hBunch->Fill(fInputEvent->GetBunchCrossNumber());
    hVertexZ->Fill(fInputEvent->GetPrimaryVertex()->GetZ());
    hNContributorsVertex->Fill(fEventCuts->GetNumberOfContributorsVtx(fInputEvent));
    // CountTracks();
    // hNGoodESDTracks->Fill(fNumberOfESDTracks);
    hNGoodESDTracks->Fill(fV0Reader->GetNumberOfPrimaryTracks());
    hNV0Tracks->Fill(fInputEvent->GetVZEROData()->GetMTotV0A()+fInputEvent->GetVZEROData()->GetMTotV0C());
  }
    
    if(fMCEvent && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){
      RelabelAODPhotonCandidates(kTRUE);  // In case of AODMC relabeling MC
      fV0Reader->RelabelAODs(kTRUE);
    }
    
  // reduce event statistics in the tree by a factor ffilltree
  // Bool_t ffillTreeNew = kFALSE;
  // if(ffillTree>=1.0) {
  //   ffillTreeNew = kTRUE;
  //   if (ffillTree>1.0) {
  //     gRandom->SetSeed(0);
  //     if(gRandom->Uniform(ffillTree)>1.0) {
  // 	ffillTreeNew = kFALSE;
  //     }
  //   }
  // }

  for(Int_t firstGammaIndex=0;firstGammaIndex<fConversionGammas->GetEntriesFast();firstGammaIndex++){
    AliAODConversionPhoton *gamma=dynamic_cast<AliAODConversionPhoton*>(fConversionGammas->At(firstGammaIndex));
    if (gamma==NULL) continue;
    if(fMCEvent && fEventCuts->GetSignalRejection() != 0){
      if(!fEventCuts->IsParticleFromBGEvent(gamma->GetMCLabelPositive(), fMCEvent, fInputEvent))
        continue;
      if(!fEventCuts->IsParticleFromBGEvent(gamma->GetMCLabelNegative(), fMCEvent, fInputEvent))
        continue;
    }
    if(!fConversionCuts->PhotonIsSelected(gamma,fInputEvent)){
      continue;
    }

    if(fWriteVariableTree) ProcessQATree(gamma);
    if(ffillHistograms)    ProcessQA(gamma);
  }
  
  if(fMCEvent && fInputEvent->IsA()==AliAODEvent::Class() && !(fV0Reader->AreAODsRelabeled())){
    RelabelAODPhotonCandidates(kFALSE); // Back to ESDMC Label
    fV0Reader->RelabelAODs(kFALSE);
  }
    
  PostData(1, fOutputList);
  PostData(2, fTree);
      
}


///________________________________________________________________________
void AliAnalysisTaskQA::ProcessQATree(AliAODConversionPhoton *gamma){

  // Fill Histograms for QA and MC
  AliVEvent* event = (AliVEvent*) InputEvent();
    
  AliPIDResponse* pidResonse = ((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetPIDResponse();
  
  AliMultSelection *MultSelection = (AliMultSelection*)InputEvent()->FindListObject("MultSelection");
  Centrality = MultSelection->GetMultiplicityPercentile("V0A");
  VertexZ    = InputEvent()->GetPrimaryVertex()->GetZ();
  Bunch      = InputEvent()->GetBunchCrossNumber();
  GoodESDTracks = fV0Reader->GetNumberOfPrimaryTracks();
  theta      = gamma->Theta();
  pt         = gamma->GetPhotonPt();
  phi        = gamma->GetPhotonPhi();
  chi2       = gamma->GetChi2perNDF();
  qt         = gamma->GetArmenterosQt();
  alpha      = gamma->GetArmenterosAlpha();
  psipair    = gamma->GetPsiPair();
  cosPA      = fConversionCuts->GetCosineOfPointingAngle(gamma,event);
  InvMass    = gamma->GetInvMassPair();
  X          = gamma->GetConversionX();
  Y          = gamma->GetConversionY();
  Z          = gamma->GetConversionZ();
  R          = gamma->GetConversionRadius();
  Qual       = gamma->GetPhotonQuality();
  DCAz       = gamma->GetDCAzToPrimVtx();
  DCAr       = gamma->GetDCArToPrimVtx();

  
  AliVTrack * negTrack = fConversionCuts->GetTrack(event, gamma->GetTrackLabelNegative());
  AliVTrack * posTrack = fConversionCuts->GetTrack(event, gamma->GetTrackLabelPositive());

  
  if(!negTrack||!posTrack)return;

  fKind = 9;
  if(fMCEvent && fInputEvent->IsA()==AliESDEvent::Class()){
    fKind = IsTruePhotonESD(gamma);
  } else if (fMCEvent && fInputEvent->IsA()==AliAODEvent::Class()){
  // 	  cout << "entering IsTruePhotonAOD" << endl;
    fKind = IsTruePhotonAOD(gamma);   
  }

  fDaughterProp(0) =  posTrack->Pt();
  fDaughterProp(7) =  negTrack->Pt();
  fDaughterProp(1) =  posTrack->Theta();
  fDaughterProp(8) =  negTrack->Theta();
  // dEdx TPC
  fDaughterProp(2) =  posTrack->GetTPCsignal();
  fDaughterProp(3) =  pidResonse->NumberOfSigmasTPC(posTrack,AliPID::kElectron);
  fDaughterProp(22) =  pidResonse->NumberOfSigmasTPC(posTrack,AliPID::kPion);
  fDaughterProp(9) =  negTrack->GetTPCsignal();
  fDaughterProp(10) =  pidResonse->NumberOfSigmasTPC(negTrack,AliPID::kElectron);
  fDaughterProp(23) =  pidResonse->NumberOfSigmasTPC(negTrack,AliPID::kPion);
  Int_t nPosClusterITS = 0;
  Int_t nNegClusterITS = 0;
  for(Int_t itsLayer = 0; itsLayer<6;itsLayer++){
    if(TESTBIT(negTrack->GetITSClusterMap(),itsLayer)){
      nNegClusterITS++;
    }
    if(TESTBIT(posTrack->GetITSClusterMap(),itsLayer)){
      nPosClusterITS++;
    }
  }
  
  // ITS signal
  fDaughterProp(14) =  (Float_t)nPosClusterITS;
  fDaughterProp(15) =  (Float_t)nNegClusterITS;
  if (nPosClusterITS > 0 ){
    fDaughterProp(16) =  posTrack->GetITSsignal();
    fDaughterProp(20) =  pidResonse->NumberOfSigmasITS(posTrack,AliPID::kElectron);
  } else {
    fDaughterProp(16) =  1000;
    fDaughterProp(20) =  20;
  }
  if (nNegClusterITS > 0 ){
    fDaughterProp(17) =  negTrack->GetITSsignal();
    fDaughterProp(21) =  pidResonse->NumberOfSigmasITS(negTrack,AliPID::kElectron);
  } else {
    fDaughterProp(17) =  1000;
    fDaughterProp(21) =  20;
  }

  // TOF 
  if((posTrack->GetStatus() & AliESDtrack::kTOFpid) && !(posTrack->GetStatus() & AliESDtrack::kTOFmismatch)){
    Double_t t0pos = pidResonse->GetTOFResponse().GetStartTime(posTrack->P());
    Double_t timesPos[9];
    posTrack->GetIntegratedTimes(timesPos,9);
    Double_t TOFsignalPos =	posTrack->GetTOFsignal();
    Double_t dTpos = TOFsignalPos - t0pos - timesPos[0];
    fDaughterProp(4) =  dTpos;
    fDaughterProp(5) =  pidResonse->NumberOfSigmasTOF(posTrack, AliPID::kElectron);
  } else {
    fDaughterProp(4) =  20000;
    fDaughterProp(5) =  -20;
  }
  if((negTrack->GetStatus() & AliESDtrack::kTOFpid) && !(negTrack->GetStatus() & AliESDtrack::kTOFmismatch)){
    Double_t t0neg = pidResonse->GetTOFResponse().GetStartTime(negTrack->P());
    Double_t timesNeg[9];
    negTrack->GetIntegratedTimes(timesNeg,9);
    Double_t TOFsignalNeg =	negTrack->GetTOFsignal();
    Double_t dTneg = TOFsignalNeg - t0neg - timesNeg[0];
    fDaughterProp(11) =  dTneg;
    fDaughterProp(12) =  pidResonse->NumberOfSigmasTOF(negTrack, AliPID::kElectron);
  } else {
    fDaughterProp(11) =  20000;
    fDaughterProp(12) =  -20;
  }

  fDaughterProp(6) =  (Float_t)posTrack->GetTPCClusterInfo(2,0,fConversionCuts->GetFirstTPCRow(gamma->GetConversionRadius()));
  fDaughterProp(18) =  posTrack->GetNcls(1);
  fDaughterProp(13) =  (Float_t)negTrack->GetTPCClusterInfo(2,0,fConversionCuts->GetFirstTPCRow(gamma->GetConversionRadius()));
  fDaughterProp(19) =  negTrack->GetNcls(1);
  
  // if (fTreeQA){
  //   fTreeQA->Fill();
  // }

  ele_NfindableCls  =  (Float_t)posTrack->GetTPCClusterInfo(2,0,fConversionCuts->GetFirstTPCRow(gamma->GetConversionRadius()));
  pos_Cls           =  posTrack->GetNcls(1);
  pos_NfindableCls  =  (Float_t)negTrack->GetTPCClusterInfo(2,0,fConversionCuts->GetFirstTPCRow(gamma->GetConversionRadius()));

  ele_theta          =  negTrack->Theta();
  ele_pt             =  negTrack->Pt();
  ele_phi            =  negTrack->Phi();
  ele_nSigmaTPC      =  pidResonse->NumberOfSigmasTPC(negTrack,AliPID::kElectron);
  ele_nSigmaTPCpion  =  pidResonse->NumberOfSigmasTPC(negTrack,AliPID::kPion);
  ele_TPCsignal      =  negTrack->GetTPCsignal();
  ele_Cls            =  negTrack->GetNcls(1);

  pos_theta          =  posTrack->Theta();
  pos_pt             =  posTrack->Pt();
  pos_phi            =  posTrack->Phi();
  pos_nSigmaTPC      =  pidResonse->NumberOfSigmasTPC(posTrack,AliPID::kElectron);
  pos_nSigmaTPCpion  =  pidResonse->NumberOfSigmasTPC(posTrack,AliPID::kPion);
  pos_TPCsignal      =  posTrack->GetTPCsignal();
  pos_Cls            =  posTrack->GetNcls(1);

 if(fInputEvent->IsA()==AliESDEvent::Class()){
    ele_SPD1 = (dynamic_cast<AliESDtrack*>(negTrack))->HasPointOnITSLayer(0);
    ele_SPD2 = (dynamic_cast<AliESDtrack*>(negTrack))->HasPointOnITSLayer(1);
    ele_SDD1 = (dynamic_cast<AliESDtrack*>(negTrack))->HasPointOnITSLayer(2);
    ele_SDD2 = (dynamic_cast<AliESDtrack*>(negTrack))->HasPointOnITSLayer(3);
    ele_SSD1 = (dynamic_cast<AliESDtrack*>(negTrack))->HasPointOnITSLayer(4);
    ele_SSD2 = (dynamic_cast<AliESDtrack*>(negTrack))->HasPointOnITSLayer(5);
    pos_SPD1 = (dynamic_cast<AliESDtrack*>(posTrack))->HasPointOnITSLayer(0);
    pos_SPD2 = (dynamic_cast<AliESDtrack*>(posTrack))->HasPointOnITSLayer(1);
    pos_SDD1 = (dynamic_cast<AliESDtrack*>(posTrack))->HasPointOnITSLayer(2);
    pos_SDD2 = (dynamic_cast<AliESDtrack*>(posTrack))->HasPointOnITSLayer(3);
    pos_SSD1 = (dynamic_cast<AliESDtrack*>(posTrack))->HasPointOnITSLayer(4);
    pos_SSD2 = (dynamic_cast<AliESDtrack*>(posTrack))->HasPointOnITSLayer(5);
  }else if ( fInputEvent->IsA()==AliAODEvent::Class()){
    ele_SPD1 = (dynamic_cast<AliAODTrack*>(negTrack))->HasPointOnITSLayer(0);
    ele_SPD2 = (dynamic_cast<AliAODTrack*>(negTrack))->HasPointOnITSLayer(1);
    ele_SDD1 = (dynamic_cast<AliAODTrack*>(negTrack))->HasPointOnITSLayer(2);
    ele_SDD2 = (dynamic_cast<AliAODTrack*>(negTrack))->HasPointOnITSLayer(3);
    ele_SSD1 = (dynamic_cast<AliAODTrack*>(negTrack))->HasPointOnITSLayer(4);
    ele_SSD2 = (dynamic_cast<AliAODTrack*>(negTrack))->HasPointOnITSLayer(5);
    pos_SPD1 = (dynamic_cast<AliAODTrack*>(posTrack))->HasPointOnITSLayer(0);
    pos_SPD2 = (dynamic_cast<AliAODTrack*>(posTrack))->HasPointOnITSLayer(1);
    pos_SDD1 = (dynamic_cast<AliAODTrack*>(posTrack))->HasPointOnITSLayer(2);
    pos_SDD2 = (dynamic_cast<AliAODTrack*>(posTrack))->HasPointOnITSLayer(3);
    pos_SSD1 = (dynamic_cast<AliAODTrack*>(posTrack))->HasPointOnITSLayer(4);
    pos_SSD2 = (dynamic_cast<AliAODTrack*>(posTrack))->HasPointOnITSLayer(5);    
  }  
  //____________________________TOF dEdx____________________________
  Bool_t nTOFisthere = kFALSE;
  //  Bool_t nTOFNoSignal = kFALSE;
  //  Bool_t nTOFMismatch = kFALSE;

  AliPIDResponse::EDetPidStatus statusnTOF = pidResonse->CheckPIDStatus(AliPIDResponse::kTOF,negTrack);
  if(statusnTOF == AliPIDResponse::kDetPidOk) nTOFisthere = kTRUE;
  //  if(statusnTOF == AliPIDResponse::kDetNoSignal) nTOFNoSignal = kTRUE;
  //  if(statusnTOF == AliPIDResponse::kDetMismatch) nTOFMismatch = kTRUE;
  
  //  if((negTrack->GetStatus() & AliESDtrack::kTOFfPID)==0 && !(negTrack->GetStatus() & AliESDtrack::kTOFmismatch)){
  if(nTOFisthere){
    Double_t timesNeg[9];
    negTrack->GetIntegratedTimes(timesNeg,9);
    Double_t dTneg = negTrack->GetTOFsignal() - pidResonse->GetTOFResponse().GetStartTime(negTrack->P()) - timesNeg[0];
    ele_TOFsignal = dTneg;
    ele_nSigmaTOF = pidResonse->NumberOfSigmasTOF(negTrack, AliPID::kElectron);
    //    cout << "ProcessQATree " << negTrack->P() << " " << ele_nSigmaTOF <<endl;
  }else{
    // cout << "ProcessQATree nTOF is not there "<<endl;  
  }

  Bool_t pTOFisthere = kFALSE;
  //  Bool_t pTOFNoSignal = kFALSE;
  //  Bool_t pTOFMismatch = kFALSE;
  AliPIDResponse::EDetPidStatus statuspTOF = pidResonse->CheckPIDStatus(AliPIDResponse::kTOF,posTrack);
  if(statuspTOF == AliPIDResponse::kDetPidOk) pTOFisthere = kTRUE;
  //  if(statuspTOF == AliPIDResponse::kDetNoSignal) pTOFNoSignal = kTRUE;
  //  if(statuspTOF == AliPIDResponse::kDetMismatch) pTOFMismatch = kTRUE;

  //  if((posTrack->GetStatus() & AliESDtrack::kTOFfPID)==0 && !(posTrack->GetStatus() & AliESDtrack::kTOFmismatch)){
  if(pTOFisthere){
    Double_t timesPos[9];
    posTrack->GetIntegratedTimes(timesPos,9);
    Double_t dTpos = posTrack->GetTOFsignal() - pidResonse->GetTOFResponse().GetStartTime(posTrack->P()) - timesPos[0];
    pos_TOFsignal  = dTpos;
    pos_nSigmaTOF  = pidResonse->NumberOfSigmasTOF(posTrack, AliPID::kElectron);
  }

  if(fWriteVariableTree) fTree->Fill();
  
}

//_____________________________________________________________________________________________________
void AliAnalysisTaskQA::ProcessQA(AliAODConversionPhoton *gamma){

  AliPIDResponse* pidResonse = ((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetPIDResponse();
  
  // Fill Histograms for QA and MC

  hGammaPt->Fill(gamma->GetPhotonPt());
  hGammaPhi->Fill(gamma->GetPhotonPhi());
  if(gamma->Eta() >= 0.00001){hGammaPhi_Pos->Fill(gamma->Phi());}
  if(gamma->Eta() <= 0.00001){hGammaPhi_Neg->Fill(gamma->Phi());}
  hGammaEta->Fill(gamma->Eta());
  hGammaChi2perNDF->Fill(gamma->GetChi2perNDF());
  hGammaPsiPair->Fill(gamma->GetPsiPair());
  hGammaArmenteros->Fill(gamma->GetArmenterosAlpha(),gamma->GetArmenterosQt());
  hGammaCosinePointingAngle->Fill(fConversionCuts->GetCosineOfPointingAngle(gamma,fInputEvent));
  hGammaInvMass->Fill(gamma->GetMass());
  hInvMassPair->Fill(gamma->GetInvMassPair(),gamma->GetPhotonPt());
  //  hGammaXY->Fill(gamma->GetConversionX(),gamma->GetConversionY());
  //  hGammaZR->Fill(gamma->GetConversionZ(),gamma->GetConversionRadius());

  AliVTrack * negTrack = fConversionCuts->GetTrack(fInputEvent, gamma->GetTrackLabelNegative());
  AliVTrack * posTrack = fConversionCuts->GetTrack(fInputEvent, gamma->GetTrackLabelPositive());
  if(!negTrack||!posTrack)return;


  hElecPt->Fill(negTrack->Pt(),posTrack->Pt());
  hElecEta->Fill(negTrack->Eta(),posTrack->Eta());
  hElecPhi->Fill(negTrack->Phi(),posTrack->Phi());

  hElecNfindableClsTPC->Fill((Float_t)posTrack->GetTPCClusterInfo(2,0,fConversionCuts->GetFirstTPCRow(gamma->GetConversionRadius())));
  hPosiNfindableClsTPC->Fill((Float_t)negTrack->GetTPCClusterInfo(2,0,fConversionCuts->GetFirstTPCRow(gamma->GetConversionRadius())));
  hElecClsTPC->Fill(negTrack->GetNcls(1));
  hPosiClsTPC->Fill(posTrack->GetNcls(1));
  //TPC dEdx
  hElectrondEdxP->Fill(negTrack->P() ,negTrack->GetTPCsignal());
  hElectronNSigmadEdxP->Fill(negTrack->P() ,pidResonse->NumberOfSigmasTPC(negTrack, AliPID::kElectron));
  hElectronNSigmadEdxEta->Fill(negTrack->Eta() ,pidResonse->NumberOfSigmasTPC(negTrack, AliPID::kElectron));
  hElectronNSigmaPiondEdxP->Fill(negTrack->P() ,pidResonse->NumberOfSigmasTPC(negTrack, AliPID::kPion));
  hPositrondEdxP->Fill(posTrack->P() ,posTrack->GetTPCsignal());
  hPositronNSigmadEdxP->Fill(posTrack->P() ,pidResonse->NumberOfSigmasTPC(posTrack, AliPID::kElectron));
  hPositronNSigmadEdxEta->Fill(posTrack->Eta() ,pidResonse->NumberOfSigmasTPC(posTrack, AliPID::kElectron));
  hPositronNSigmaPiondEdxP->Fill(posTrack->P() ,pidResonse->NumberOfSigmasTPC(posTrack, AliPID::kPion));
  
  //TOF signal
  if((negTrack->GetStatus() & AliESDtrack::kTOFpid)==0 && !(negTrack->GetStatus() & AliESDtrack::kTOFmismatch)){
    Double_t t0neg = pidResonse->GetTOFResponse().GetStartTime(negTrack->P());
    Double_t timesNeg[9];
    negTrack->GetIntegratedTimes(timesNeg,9);
    Double_t TOFsignalNeg = negTrack->GetTOFsignal();
    Double_t dTneg = TOFsignalNeg - t0neg - timesNeg[0];
    hElectronTOFP->Fill(negTrack->P() ,dTneg);
    hElectronNSigmaTOFP->Fill(negTrack->P() ,pidResonse->NumberOfSigmasTOF(negTrack, AliPID::kElectron));
    //    cout << ">>>>ProcessQA " <<negTrack->P() << " " << pidResonse->NumberOfSigmasTOF(negTrack, AliPID::kElectron) << endl;
  }
  if((posTrack->GetStatus() & AliESDtrack::kTOFpid)==0 && !(posTrack->GetStatus() & AliESDtrack::kTOFmismatch)){
    Double_t t0pos = pidResonse->GetTOFResponse().GetStartTime(posTrack->P());
    Double_t timesPos[9];
    posTrack->GetIntegratedTimes(timesPos,9);
    Double_t TOFsignalPos = posTrack->GetTOFsignal();
    Double_t dTpos = TOFsignalPos - t0pos - timesPos[0];
    hPositronTOFP->Fill(posTrack->P() ,dTpos);
    hPositronNSigmaTOFP->Fill(posTrack->P() ,pidResonse->NumberOfSigmasTOF(posTrack, AliPID::kElectron));
  }
  
  Int_t nPosClusterITS = 0;
  Int_t nNegClusterITS = 0;
  for(Int_t itsLayer = 0; itsLayer<6;itsLayer++){
    if(TESTBIT(negTrack->GetITSClusterMap(),itsLayer)){
      nNegClusterITS++;
    }
    if(TESTBIT(posTrack->GetITSClusterMap(),itsLayer)){
      nPosClusterITS++;
    }
  }
  Double_t negtrackPhi = negTrack->Phi();
  Double_t postrackPhi = posTrack->Phi();
  hITSClusterPhi->Fill(negtrackPhi,nNegClusterITS);
  hITSClusterPhi->Fill(postrackPhi,nPosClusterITS);
  
  // ITS signal
  if (nPosClusterITS > 0 ){
    hPositronITSdEdxP->Fill(posTrack->P() ,posTrack->GetITSsignal());
    hPositronNSigmaITSP->Fill(posTrack->P() ,pidResonse->NumberOfSigmasITS(posTrack,AliPID::kElectron));
  } 
  if (nNegClusterITS > 0 ){
    hElectronITSdEdxP->Fill(negTrack->P() ,negTrack->GetITSsignal());
    hElectronNSigmaITSP->Fill(negTrack->P() ,pidResonse->NumberOfSigmasITS(negTrack,AliPID::kElectron));
  }

  
}


// // //________________________________________________________________________
// void AliAnalysisTaskQA::CountTracks(){

//   if(fInputEvent->IsA()==AliESDEvent::Class()){
//     // Using standard function for setting Cuts
//     Bool_t selectPrimaries=kTRUE;
//     AliESDtrackCuts *EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(selectPrimaries);
//     EsdTrackCuts->SetMaxDCAToVertexZ(2);
//     EsdTrackCuts->SetEtaRange(-0.8, 0.8);
//     EsdTrackCuts->SetPtRange(0.15);
//     fNumberOfESDTracks = 0;
//     for(Int_t iTracks = 0; iTracks < fInputEvent->GetNumberOfTracks(); iTracks++){
//       AliESDtrack* curTrack = (AliESDtrack*) fInputEvent->GetTrack(iTracks);
//       if(!curTrack) continue;
//       if(EsdTrackCuts->AcceptTrack(curTrack) ) fNumberOfESDTracks++;
//     }
//     delete EsdTrackCuts;
//     EsdTrackCuts=0x0;
//   }
//   else if(fInputEvent->IsA()==AliAODEvent::Class()){    
//     fNumberOfESDTracks = 0;
//     for(Int_t iTracks = 0; iTracks<fInputEvent->GetNumberOfTracks(); iTracks++){
//       AliAODTrack* curTrack = (AliAODTrack*) fInputEvent->GetTrack(iTracks);
//       if(curTrack->GetID()<0) continue; // Avoid double counting of tracks
//       if(!curTrack->IsHybridGlobalConstrainedGlobal()) continue;
//       if(TMath::Abs(curTrack->Eta())>0.8) continue;
//       if(curTrack->Pt()<0.15) continue;
//       fNumberOfESDTracks++;
//     }
//   }
//   return;
// }

//________________________________________________________________________
UInt_t AliAnalysisTaskQA::IsTruePhotonESD(AliAODConversionPhoton *TruePhotonCandidate)
{
  UInt_t kind = 9;
  TParticle *posDaughter = TruePhotonCandidate->GetPositiveMCDaughter(fMCEvent);
  TParticle *negDaughter = TruePhotonCandidate->GetNegativeMCDaughter(fMCEvent);
  Int_t pdgCodePos = 0; 
  Int_t pdgCodeNeg = 0; 
  Int_t pdgCode = 0; 

  const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
  Double_t mcProdVtxX 	= primVtxMC->GetX();
  Double_t mcProdVtxY 	= primVtxMC->GetY();
  Double_t mcProdVtxZ 	= primVtxMC->GetZ();

  
  if(posDaughter == NULL || negDaughter == NULL) {
    kind = 9;
    //		return kFALSE; // One particle does not exist
  
  } else if( posDaughter->GetMother(0) != negDaughter->GetMother(0)  || (posDaughter->GetMother(0) == negDaughter->GetMother(0) && posDaughter->GetMother(0) ==-1)) {
    kind = 1;
    // 	  	return 1;
    pdgCodePos=TMath::Abs(posDaughter->GetPdgCode());
    pdgCodeNeg=TMath::Abs(negDaughter->GetPdgCode());
    if(pdgCodePos==11 && pdgCodeNeg==11) return 10; //Electron Combinatorial
    if(pdgCodePos==11 && pdgCodeNeg==11 && 
      (posDaughter->GetMother(0) == negDaughter->GetMother(0) && posDaughter->GetMother(0) ==-1)) return 15; //direct Electron Combinatorial
        
    if(pdgCodePos==211 && pdgCodeNeg==211) kind = 11; //Pion Combinatorial
    if((pdgCodePos==211 && pdgCodeNeg==2212) ||(pdgCodePos==2212 && pdgCodeNeg==211))	kind = 12; //Pion, Proton Combinatorics
    if((pdgCodePos==11 && pdgCodeNeg==2212) ||(pdgCodePos==2212 && pdgCodeNeg==11))	kind = 16; //electron, Proton Combinatorics
    if((pdgCodePos==11 && pdgCodeNeg==321) ||(pdgCodePos==321 && pdgCodeNeg==11))	kind = 17; //electron, kaon
    if((pdgCodePos==211 && pdgCodeNeg==321) ||(pdgCodePos==321 && pdgCodeNeg==211))	kind = 18; //pion, kaon
    if((pdgCodePos==211 && pdgCodeNeg==11) ||(pdgCodePos==11 && pdgCodeNeg==211)) kind = 13; //Pion, Electron Combinatorics
    if(pdgCodePos==321 && pdgCodeNeg==321) kind = 14; //Kaon,Kaon combinatorics
  }else{		
    pdgCodePos=posDaughter->GetPdgCode();
    pdgCodeNeg=negDaughter->GetPdgCode();
    Bool_t gammaIsPrimary = fEventCuts->IsConversionPrimaryESD( fMCEvent, posDaughter->GetMother(0), mcProdVtxX, mcProdVtxY, mcProdVtxZ);
    if ( TruePhotonCandidate->GetMCParticle(fMCEvent)->GetPdgCode()) pdgCode = TruePhotonCandidate->GetMCParticle(fMCEvent)->GetPdgCode();

    if(TMath::Abs(pdgCodePos)!=11 || TMath::Abs(pdgCodeNeg)!=11) return 2; // true from hadronic decays
    else if ( !(pdgCodeNeg==pdgCodePos)){
      if(pdgCode == 111) return 3; // pi0 Dalitz
      else if (pdgCode == 221) return 4; // eta Dalitz
      else if (!(negDaughter->GetUniqueID() != 5 || posDaughter->GetUniqueID() !=5)){
        if(pdgCode == 22 && gammaIsPrimary){
          return 0; // primary photons
        } else if (pdgCode == 22){
          return 5; //secondary photons
        }
      }
    }
  }

  return kind;
}

//________________________________________________________________________
UInt_t AliAnalysisTaskQA::IsTruePhotonAOD(AliAODConversionPhoton *TruePhotonCandidate)
{   

  UInt_t kind = 9;
  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
  if (AODMCTrackArray!=NULL && TruePhotonCandidate!=NULL){
    AliAODMCParticle *posDaughter = (AliAODMCParticle*) AODMCTrackArray->At(TruePhotonCandidate->GetMCLabelPositive());
    AliAODMCParticle *negDaughter = (AliAODMCParticle*) AODMCTrackArray->At(TruePhotonCandidate->GetMCLabelNegative());
    Int_t pdgCodePos = 0; 
    Int_t pdgCodeNeg = 0; 
    Int_t pdgCode = 0; 
    if(posDaughter == NULL || negDaughter == NULL) {
      kind = 9;
    } else if( posDaughter->GetMother() != negDaughter->GetMother()  || (posDaughter->GetMother() == negDaughter->GetMother() && posDaughter->GetMother() ==-1)) {
      kind = 1;
      pdgCodePos=TMath::Abs(posDaughter->GetPdgCode());
      pdgCodeNeg=TMath::Abs(negDaughter->GetPdgCode());
      if(pdgCodePos==11 && pdgCodeNeg==11)	kind = 10; //Electron Combinatorial
      if(pdgCodePos==11 && pdgCodeNeg==11 && 
        (posDaughter->GetMother() == negDaughter->GetMother() && posDaughter->GetMother() ==-1))kind = 15; //direct Electron Combinatorial
          
      if(pdgCodePos==211 && pdgCodeNeg==211) kind = 11; //Pion Combinatorial
      if((pdgCodePos==211 && pdgCodeNeg==2212) ||(pdgCodePos==2212 && pdgCodeNeg==211))	kind = 12; //Pion, Proton Combinatorics
      if((pdgCodePos==11 && pdgCodeNeg==2212) ||(pdgCodePos==2212 && pdgCodeNeg==11))	kind = 16; //electron, Proton Combinatorics
      if((pdgCodePos==11 && pdgCodeNeg==321) ||(pdgCodePos==321 && pdgCodeNeg==11))	kind = 17; //electron, kaon
      if((pdgCodePos==211 && pdgCodeNeg==321) ||(pdgCodePos==321 && pdgCodeNeg==211))	kind = 18; //pion, kaon
      if((pdgCodePos==211 && pdgCodeNeg==11) ||(pdgCodePos==11 && pdgCodeNeg==211)) kind = 13; //Pion, Electron Combinatorics
      if(pdgCodePos==321 && pdgCodeNeg==321) kind = 14; //Kaon,Kaon combinatorics
    }else{		
      AliAODMCParticle *Photon = (AliAODMCParticle*) AODMCTrackArray->At(posDaughter->GetMother());
      pdgCodePos=posDaughter->GetPdgCode();
      pdgCodeNeg=negDaughter->GetPdgCode();

      if ( Photon->GetPdgCode()) 
        pdgCode = Photon->GetPdgCode(); 
      if(TMath::Abs(pdgCodePos)!=11 || TMath::Abs(pdgCodeNeg)!=11) kind = 2; // true from hadronic decays
      else if ( !(pdgCodeNeg==pdgCodePos)){
        if(pdgCode == 111) kind = 3; // pi0 Dalitz
        else if (pdgCode == 221) kind = 4; // eta Dalitz
        else if (!(negDaughter->GetMCProcessCode() != 5 || posDaughter->GetMCProcessCode() !=5)){
          const AliVVertex* primVtxMC 	= fMCEvent->GetPrimaryVertex();
          Double_t mcProdVtxX 	= primVtxMC->GetX();
          Double_t mcProdVtxY 	= primVtxMC->GetY();
          Double_t mcProdVtxZ 	= primVtxMC->GetZ();
          Bool_t isPrimary = fEventCuts->IsConversionPrimaryAOD(fInputEvent, Photon, mcProdVtxX, mcProdVtxY, mcProdVtxZ);

          if(pdgCode == 22 && isPrimary){
            kind = 0; // primary photons
          } else if (pdgCode == 22){
            kind = 5; //secondary photons
          }
        }
      }
    }

    return kind;
  }	
  return kind;
}

//________________________________________________________________________
void AliAnalysisTaskQA::RelabelAODPhotonCandidates(Bool_t mode){

  // Relabeling For AOD Event
  // ESDiD -> AODiD
  // MCLabel -> AODMCLabel
  
  if(mode){
    fMCStackPos = new Int_t[fConversionGammas->GetEntries()];
    fMCStackNeg = new Int_t[fConversionGammas->GetEntries()];
  }
  
  for(Int_t iGamma = 0;iGamma<fConversionGammas->GetEntries();iGamma++){
    AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fConversionGammas->At(iGamma);
    if(!PhotonCandidate) continue;
    if(!mode){// Back to ESD Labels
      PhotonCandidate->SetMCLabelPositive(fMCStackPos[iGamma]);
      PhotonCandidate->SetMCLabelNegative(fMCStackNeg[iGamma]);
      //PhotonCandidate->IsAODMCLabel(kFALSE);
      continue;
    }
    fMCStackPos[iGamma] =  PhotonCandidate->GetMCLabelPositive();
    fMCStackNeg[iGamma] =  PhotonCandidate->GetMCLabelNegative();

    Bool_t AODLabelPos = kFALSE;
    Bool_t AODLabelNeg = kFALSE;

    for(Int_t i = 0; i<fInputEvent->GetNumberOfTracks();i++){
      AliAODTrack *tempDaughter = static_cast<AliAODTrack*>(fInputEvent->GetTrack(i));
      if(!AODLabelPos){
        if( tempDaughter->GetID() == PhotonCandidate->GetTrackLabelPositive() ){
        PhotonCandidate->SetMCLabelPositive(TMath::Abs(tempDaughter->GetLabel()));
        AODLabelPos = kTRUE;
        }
      }
      if(!AODLabelNeg){
        if( tempDaughter->GetID() == PhotonCandidate->GetTrackLabelNegative()){
        PhotonCandidate->SetMCLabelNegative(TMath::Abs(tempDaughter->GetLabel()));
        AODLabelNeg = kTRUE;
        }
      }
      if(AODLabelNeg && AODLabelPos){
        break;
      }
    } // Both ESD Tracks have AOD Tracks with Positive IDs
    if(!AODLabelPos || !AODLabelNeg){
      for(Int_t i = 0; i<fInputEvent->GetNumberOfTracks();i++){
        AliAODTrack *tempDaughter = static_cast<AliAODTrack*>(fInputEvent->GetTrack(i));
        if(tempDaughter->GetID()<0){
        if(!AODLabelPos){
          if( (TMath::Abs(tempDaughter->GetID())-1) == PhotonCandidate->GetTrackLabelPositive()){
            PhotonCandidate->SetMCLabelPositive(TMath::Abs(tempDaughter->GetLabel()));
            AODLabelPos = kTRUE;
          }
        }
        if(!AODLabelNeg){
          if( (TMath::Abs(tempDaughter->GetID())-1) == PhotonCandidate->GetTrackLabelNegative()){
            PhotonCandidate->SetMCLabelNegative(TMath::Abs(tempDaughter->GetLabel()));
            AODLabelNeg = kTRUE;
          }
        }
        }
        if(AODLabelNeg && AODLabelPos){
        break;
        }
      }
      if(!AODLabelPos || !AODLabelNeg){
	//        cout<<"WARNING!!! AOD TRACKS NOT FOUND FOR"<<endl;
        if(!AODLabelNeg){
          PhotonCandidate->SetMCLabelNegative(-999999);
          PhotonCandidate->SetLabelNegative(-999999);
        }
        if(!AODLabelPos){
          PhotonCandidate->SetMCLabelPositive(-999999);
          PhotonCandidate->SetLabelPositive(-999999);
        }
      }
    }
  }
  
  if(!mode){
    delete[] fMCStackPos;
    delete[] fMCStackNeg;
  }
}

//________________________________________________________________________
void AliAnalysisTaskQA::SetLogBinningXTH2(TH2* histoRebin){
  TAxis *axisafter = histoRebin->GetXaxis(); 
  Int_t bins = axisafter->GetNbins();
  Double_t from = axisafter->GetXmin();
  Double_t to = axisafter->GetXmax();
  Double_t *newbins = new Double_t[bins+1];
  newbins[0] = from;
  Double_t factor = TMath::Power(to/from, 1./bins);
  for(Int_t i=1; i<=bins; ++i) newbins[i] = factor * newbins[i-1];
  axisafter->Set(bins, newbins);
  delete [] newbins;

}

//________________________________________________________________________
void AliAnalysisTaskQA::Terminate(Option_t *)
{

}
