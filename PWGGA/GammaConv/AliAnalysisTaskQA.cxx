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
AliAnalysisTaskQA::AliAnalysisTaskQA() : AliAnalysisTaskSE()
  ,fV0Reader(NULL)
  ,fV0ReaderName("V0ReaderV1")
  ,fConversionGammas(NULL)
  ,fConversionCuts(NULL)
  ,fEventCuts(NULL)
  ,fInputEvent(NULL)
  ,fNumberOfESDTracks(0)
  ,fMCEvent(NULL)
//  fTreeQA(NULL)
  ,fIsHeavyIon(kFALSE)
//  ffillTree(-100)
  ,ffillHistograms(kFALSE)
  ,fOutputList(NULL)
  ,fESDList(NULL)
  ,hCentralityV0A(NULL)
  ,hBunch(NULL)
  ,hVertexZ(NULL)
  ,hNGoodESDTracks(NULL)
  ,hNV0Tracks(NULL)
  ,hNContributorsVertex(NULL)
  ,hITSClusterPhi(NULL)
  ,hGammaPt(NULL)
  ,hGammaPhi(NULL)
  ,hGammaPhi_Pos(NULL)
  ,hGammaPhi_Neg(NULL)
  ,hGammaEta(NULL)
  ,hGammaChi2perNDF(NULL)
  ,hGammaPsiPair(NULL)
  ,hGammaArmenteros(NULL)
  ,hGammaCosinePointingAngle(NULL)
  ,hGammaInvMass(NULL)
  ,hElecPt(NULL)
  ,hElecEta(NULL)
  ,hElecPhi(NULL)
  ,hElecNfindableClsTPC(NULL)
  ,hPosiNfindableClsTPC(NULL)
  ,hElecClsTPC(NULL)
  ,hPosiClsTPC(NULL)
  ,hElectrondEdxP(NULL)
  ,hElectronITSdEdxP(NULL)
  ,hElectronTOFP(NULL)
  ,hElectronNSigmadEdxP(NULL)
  ,hElectronNSigmadEdxEta(NULL)
  ,hElectronNSigmaPiondEdxP(NULL)
  ,hElectronNSigmaITSP(NULL)
  ,hElectronNSigmaTOFP(NULL)
  ,hPositrondEdxP(NULL)
  ,hPositronITSdEdxP(NULL)
  ,hPositronTOFP(NULL)
  ,hPositronNSigmadEdxP(NULL)
  ,hPositronNSigmadEdxEta(NULL)
  ,hPositronNSigmaPiondEdxP(NULL)
  ,hPositronNSigmaITSP(NULL)
  ,hPositronNSigmaTOFP(NULL)
  ,hInvMassPair(NULL)
  ,fTree(NULL)
  //  hElecAsymP(NULL)
  //  fTrueList(NULL)
  //  hTrueResolutionR(NULL)
  //  hTrueResolutionZ(NULL)
  //  hTrueResolutionPhi(NULL)
  //  hTrueGammaPt(NULL)
  //  hTrueGammaPhi(NULL)
  //  hTrueGammaEta(NULL)
  //  hTrueGammaMass(NULL)
  //  hTrueGammaChi2perNDF(NULL)
  //  hTrueGammaPsiPair(NULL)
  //  hTrueGammaQt(NULL)
  //  hTrueGammaCosinePointingAngle(NULL)
  //  hTrueGammaXY(NULL)
  //  hTrueGammaZR(NULL)
  //  hTrueElecPt(NULL)
  //  hTrueElecEta(NULL)
  //  hTrueElecPhi(NULL)
  //  hTrueElecNfindableClsTPC(NULL)
  //  hTruePosiNfindableClsTPC(NULL)
  //  hTrueElecAsymP(NULL)
  ,fCentralityV0M(0)
  ,fCentralityV0A(0)
  ,fCentralityV0C(0)
  ,fRunNumber(0)
  ,fVertexZ(0)
  ,fBunch(0)
  ,fGoodESDTracks(0)
  ,ftheta(0)
  ,fpt(0)   
  ,fphi(0)
  ,fchi2(0)
  ,fqt(0)
  ,falpha(0)
  ,fpsipair(0)
  ,fcosPA(0)
  ,fInvMass(0)
  ,fX(0)
  ,fY(0)
  ,fZ(0)
  ,fR(0)
  ,fQual(0)
  ,fDCAz(0)
  ,fDCAr(0)
  ,fele_theta(0)
  ,fele_pt(0)
  ,fele_phi(0)
  ,fele_nSigmaTPC(0)
  ,fele_nSigmaTPCpion(0)
  ,fele_nSigmaTOF(0)
  ,fele_nSigmaITS(0)
  ,fele_TPCsignal(0)
  ,fele_TOFsignal(0)
  ,fele_ITSsignal(0)
  ,fele_Cls(0)
  ,fele_NfindableCls(0)
  ,fele_SPD1(0)
  ,fele_SPD2(0)
  ,fele_SDD1(0)
  ,fele_SDD2(0)
  ,fele_SSD1(0)
  ,fele_SSD2(0)
  ,fpos_theta(0)
  ,fpos_pt(0)
  ,fpos_phi(0)
  ,fpos_nSigmaTPC(0)
  ,fpos_nSigmaTPCpion(0)
  ,fpos_nSigmaTOF(0)
  ,fpos_nSigmaITS(0)
  ,fpos_TPCsignal(0)
  ,fpos_TOFsignal(0)
  ,fpos_ITSsignal(0)
  ,fpos_Cls(0)
  ,fpos_NfindableCls(0)
  ,fpos_SPD1(0)
  ,fpos_SPD2(0)
  ,fpos_SDD1(0)
  ,fpos_SDD2(0)
  ,fpos_SSD1(0)
  ,fpos_SSD2(0)
  ,fKind(0)
  ,fIsMC(kFALSE)
  ,fnGammaCandidates(1)
  ,fMCStackPos(NULL)
  ,fMCStackNeg(NULL)
  ,fActiveBranches("")
  ,fInactiveBranches("")
  ,fWriteVariableTree(kFALSE)

{

}

AliAnalysisTaskQA::AliAnalysisTaskQA(const char *name) : AliAnalysisTaskSE(name)
  ,fV0Reader(NULL)
  ,fV0ReaderName("V0ReaderV1")
  ,fConversionGammas(NULL)
  ,fConversionCuts(NULL)
  ,fEventCuts(NULL)
  ,fInputEvent(NULL)
  ,fNumberOfESDTracks(0)
  ,fMCEvent(NULL)
  //  fTreeQA(NULL)
  ,fIsHeavyIon(kFALSE)
  //  ffillTree(-100)
  ,ffillHistograms(kFALSE)
  ,fOutputList(NULL)
  ,fESDList(NULL)
  ,hCentralityV0A(NULL)
  ,hBunch(NULL)
  ,hVertexZ(NULL)
  ,hNGoodESDTracks(NULL)
  ,hNV0Tracks(NULL)
  ,hNContributorsVertex(NULL)
  ,hITSClusterPhi(NULL)
  ,hGammaPt(NULL)
  ,hGammaPhi(NULL)
  ,hGammaPhi_Pos(NULL)
  ,hGammaPhi_Neg(NULL)
  ,hGammaEta(NULL)
  ,hGammaChi2perNDF(NULL)
  ,hGammaPsiPair(NULL)
  ,hGammaArmenteros(NULL)
  ,hGammaCosinePointingAngle(NULL)
  ,hGammaInvMass(NULL)
  ,hElecPt(NULL)
  ,hElecEta(NULL)
  ,hElecPhi(NULL)
  ,hElecNfindableClsTPC(NULL)
  ,hPosiNfindableClsTPC(NULL)
  ,hElecClsTPC(NULL)
  ,hPosiClsTPC(NULL)
  ,hElectrondEdxP(NULL)
  ,hElectronITSdEdxP(NULL)
  ,hElectronTOFP(NULL)
  ,hElectronNSigmadEdxP(NULL)
  ,hElectronNSigmadEdxEta(NULL)
  ,hElectronNSigmaPiondEdxP(NULL)
  ,hElectronNSigmaITSP(NULL)
  ,hElectronNSigmaTOFP(NULL)
  ,hPositrondEdxP(NULL)
  ,hPositronITSdEdxP(NULL)
  ,hPositronTOFP(NULL)
  ,hPositronNSigmadEdxP(NULL)
  ,hPositronNSigmadEdxEta(NULL)
  ,hPositronNSigmaPiondEdxP(NULL)
  ,hPositronNSigmaITSP(NULL)
  ,hPositronNSigmaTOFP(NULL)
  ,hInvMassPair(NULL)
  ,fTree(NULL)
  //  hGammaXY(NULL)
  //  hGammaZR(NULL)
  //  hElecAsymP(NULL)
  //  fTrueList(NULL)
  //  hTrueResolutionR(NULL)
  //  hTrueResolutionZ(NULL)
  //  hTrueResolutionPhi(NULL)
  //  hTrueGammaPt(NULL)
  //  hTrueGammaPhi(NULL)
  //  hTrueGammaEta(NULL)
  //  hTrueGammaMass(NULL)
  //  hTrueGammaChi2perNDF(NULL)
  //  hTrueGammaPsiPair(NULL)
  //  hTrueGammaQt(NULL)
  //  hTrueGammaCosinePointingAngle(NULL)
  //  hTrueGammaXY(NULL)
  //  hTrueGammaZR(NULL)
  //  hTrueElecPt(NULL)
  //  hTrueElecEta(NULL)
  //  hTrueElecPhi(NULL)
  //  hTrueElecNfindableClsTPC(NULL)
  //  hTruePosiNfindableClsTPC(NULL)
  //  hTrueElecAsymP(NULL)
  ,fCentralityV0M(0)
  ,fCentralityV0A(0)
  ,fCentralityV0C(0)
  ,fRunNumber(0)
  ,fVertexZ(0)
  ,fBunch(0)
  ,fGoodESDTracks(0)
  ,ftheta(0)
  ,fpt(0)
  ,fphi(0)
  ,fchi2(0)
  ,fqt(0)
  ,falpha(0)
  ,fpsipair(0)
  ,fcosPA(0)
  ,fInvMass(0)
  ,fX(0)
  ,fY(0)
  ,fZ(0)
  ,fR(0)
  ,fQual(0)
  ,fDCAz(0)
  ,fDCAr(0)									     
  ,fele_theta(0)
  ,fele_pt(0)
  ,fele_phi(0)
  ,fele_nSigmaTPC(0)
  ,fele_nSigmaTPCpion(0)
  ,fele_nSigmaTOF(0)
  ,fele_nSigmaITS(0)
  ,fele_TPCsignal(0)
  ,fele_TOFsignal(0)
  ,fele_ITSsignal(0)
  ,fele_Cls(0)
  ,fele_NfindableCls(0)
  ,fele_SPD1(0)
  ,fele_SPD2(0)
  ,fele_SDD1(0)
  ,fele_SDD2(0)
  ,fele_SSD1(0)
  ,fele_SSD2(0)
  ,fpos_theta(0)
  ,fpos_pt(0)
  ,fpos_phi(0)
  ,fpos_nSigmaTPC(0)
  ,fpos_nSigmaTPCpion(0)
  ,fpos_nSigmaTOF(0)
  ,fpos_nSigmaITS(0)
  ,fpos_TPCsignal(0)
  ,fpos_TOFsignal(0)
  ,fpos_ITSsignal(0)
  ,fpos_Cls(0)
  ,fpos_NfindableCls(0)
  ,fpos_SPD1(0)
  ,fpos_SPD2(0)
  ,fpos_SDD1(0)
  ,fpos_SDD2(0)
  ,fpos_SSD1(0)
  ,fpos_SSD2(0)
  ,fKind(0)
  ,fIsMC(kFALSE)
  ,fnGammaCandidates(1)
  ,fMCStackPos(NULL)
  ,fMCStackNeg(NULL)
  ,fActiveBranches("")
  ,fInactiveBranches("")
  ,fWriteVariableTree(kFALSE)

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

    hCentralityV0A = new TH1F("CentralityV0A","",100,0,100);
    fOutputList->Add(hCentralityV0A);
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

  fTree->Branch("CentralityV0M",      &fCentralityV0M, "fCentralityV0M/F");
  fTree->Branch("CentralityV0A",      &fCentralityV0A, "fCentralityV0A/F");
  fTree->Branch("CentralityV0C",      &fCentralityV0C, "fCentralityV0C/F");
  fTree->Branch("RunNumber",          &fRunNumber,     "fRunNumber/I");
  fTree->Branch("VertexZ",            &fVertexZ,       "fVertexZ/F");
  fTree->Branch("Bunch",              &fBunch,         "fBunch/I");
  fTree->Branch("GoodESDTracks",      &fGoodESDTracks, "fGoodESDTracks/I");

  fTree->Branch("theta",              &ftheta,         "ftheta/F");
  fTree->Branch("pt",                 &fpt,            "fpt/F");
  fTree->Branch("phi",                &fphi,           "fphi/F");
  fTree->Branch("chi2",               &fchi2,          "fchi2/F");
  fTree->Branch("qt",                 &fqt,            "fqt/F");
  fTree->Branch("alpha",              &falpha,         "falpha/F");
  fTree->Branch("psipair",            &fpsipair,       "fpsipair/F");
  fTree->Branch("cosPA",              &fcosPA,         "fcosPA/F");
  fTree->Branch("InvMass",            &fInvMass,       "fInvMass/F");
  fTree->Branch("X",                  &fX,             "fX/F");
  fTree->Branch("Y",                  &fY,             "fY/F");
  fTree->Branch("Z",                  &fZ,             "fZ/F");
  fTree->Branch("R",                  &fR,             "fR/F");
  fTree->Branch("Qual",               &fQual,          "fQual/I");
  fTree->Branch("DCAz",               &fDCAz,          "fDCAz/F");
  fTree->Branch("DCAr",               &fDCAr,          "fDCAr/F");

  fTree->Branch("ele_theta",          &fele_theta,     "fele_theta/F");
  fTree->Branch("ele_pt",             &fele_pt,        "fele_pt/F");
  fTree->Branch("ele_phi",            &fele_phi,       "fele_phi/F");
  fTree->Branch("ele_nSigmaTPC",      &fele_nSigmaTPC, "fele_nSigmaTPC/F");
  fTree->Branch("ele_nSigmaTPCpion",  &fele_nSigmaTPCpion,"fele_nSigmaTPCpion/F");
  fTree->Branch("ele_nSigmaTOF",      &fele_nSigmaTOF, "fele_nSigmaTOF/F");
  fTree->Branch("ele_nSigmaITS",      &fele_nSigmaITS, "fele_nSigmaITS/F");
  fTree->Branch("ele_TPCsignal",      &fele_TPCsignal, "fele_TPCsignal/F");
  fTree->Branch("ele_TOFsignal",      &fele_TOFsignal, "fele_TOFsignal/F");
  fTree->Branch("ele_ITSsignal",      &fele_ITSsignal, "fele_ITSsignal/F");
  fTree->Branch("ele_Cls",            &fele_Cls,       "fele_Cls/F");
  fTree->Branch("ele_SPD1",           &fele_SPD1,      "fele_SPD1/O");
  fTree->Branch("ele_SPD2",           &fele_SPD2,      "fele_SPD2/O");
  fTree->Branch("ele_SDD1",           &fele_SDD1,      "fele_SDD1/O");
  fTree->Branch("ele_SDD2",           &fele_SDD2,      "fele_SDD2/O");
  fTree->Branch("ele_SSD1",           &fele_SSD1,      "fele_SSD1/O");
  fTree->Branch("ele_SSD2",           &fele_SSD2,      "fele_SSD2/O");
          
  fTree->Branch("pos_theta",          &fpos_theta,     "fpos_theta/F");
  fTree->Branch("pos_pt",             &fpos_pt,        "fpos_pt/F");
  fTree->Branch("pos_phi",            &fpos_phi,       "fpos_phi/F");
  fTree->Branch("pos_nSigmaTPC",      &fpos_nSigmaTPC, "fpos_nSigmaTPC/F");
  fTree->Branch("pos_nSigmaTPCpion",  &fpos_nSigmaTPCpion,"fpos_nSigmaTPCpion/F");
  fTree->Branch("pos_nSigmaTOF",      &fpos_nSigmaTOF, "fpos_nSigmaTOF/F");
  fTree->Branch("pos_nSigmaITS",      &fpos_nSigmaITS, "fpos_nSigmaITS/F");
  fTree->Branch("pos_TPCsignal",      &fpos_TPCsignal, "fpos_TPCsignal/F");
  fTree->Branch("pos_TOFsignal",      &fpos_TOFsignal, "fpos_TOFsignal/F");
  fTree->Branch("pos_ITSsignal",      &fpos_ITSsignal, "fpos_ITSsignal/F");
  fTree->Branch("pos_Cls",            &fpos_Cls,       "fpos_Cls/F");
  fTree->Branch("pos_SPD1",           &fpos_SPD1,      "fpos_SPD1/O");
  fTree->Branch("pos_SPD2",           &fpos_SPD2,      "fpos_SPD2/O");
  fTree->Branch("pos_SDD1",           &fpos_SDD1,      "fpos_SDD1/O");
  fTree->Branch("pos_SDD2",           &fpos_SDD2,      "fpos_SDD2/O");
  fTree->Branch("pos_SSD1",           &fpos_SSD1,      "fpos_SSD1/O");
  fTree->Branch("pos_SSD2",           &fpos_SSD2,      "fpos_SSD2/O");

  
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
    fTree->Branch("kind",               &fKind,          "fKind/b");
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
    printf(" QA Task %s :: Eta Shift Manually Set to %f \n\n",
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

  AliMultSelection* fMultSelection;
  Double_t cent = -1.;
  if(fIsHeavyIon>0){
    fMultSelection = (AliMultSelection*)fInputEvent->FindListObject("MultSelection");
    cent=fMultSelection->GetMultiplicityPercentile("V0A");
    fCentralityV0M = fMultSelection->GetMultiplicityPercentile("V0M");
    fCentralityV0A = fMultSelection->GetMultiplicityPercentile("V0A");
    fCentralityV0C = fMultSelection->GetMultiplicityPercentile("V0C");
  }
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

  fRunNumber  = fInputEvent->GetRunNumber();
  fVertexZ    = fInputEvent->GetPrimaryVertex()->GetZ();
  fBunch      = fInputEvent->GetBunchCrossNumber();
  fGoodESDTracks = fV0Reader->GetNumberOfPrimaryTracks();
  
  if(ffillHistograms){
    if(fIsHeavyIon>0){
      hCentralityV0A->Fill(cent);
    }
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
    
  ftheta      = gamma->Theta();
  fpt         = gamma->GetPhotonPt();
  fphi        = gamma->GetPhotonPhi();
  fchi2       = gamma->GetChi2perNDF();
  fqt         = gamma->GetArmenterosQt();
  falpha      = gamma->GetArmenterosAlpha();
  fpsipair    = gamma->GetPsiPair();
  fcosPA      = fConversionCuts->GetCosineOfPointingAngle(gamma,event);
  fInvMass    = gamma->GetInvMassPair();
  fX          = gamma->GetConversionX();
  fY          = gamma->GetConversionY();
  fZ          = gamma->GetConversionZ();
  fR          = gamma->GetConversionRadius();
  fQual       = gamma->GetPhotonQuality();
  fDCAz       = gamma->GetDCAzToPrimVtx();
  fDCAr       = gamma->GetDCArToPrimVtx();

  
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

  // fDaughterProp(0) =  posTrack->Pt();
  // fDaughterProp(7) =  negTrack->Pt();
  // fDaughterProp(1) =  posTrack->Theta();
  // fDaughterProp(8) =  negTrack->Theta();
  // // dEdx TPC
  // fDaughterProp(2) =  posTrack->GetTPCsignal();
  // fDaughterProp(3) =  pidResonse->NumberOfSigmasTPC(posTrack,AliPID::kElectron);
  // fDaughterProp(22) =  pidResonse->NumberOfSigmasTPC(posTrack,AliPID::kPion);
  // fDaughterProp(9) =  negTrack->GetTPCsignal();
  // fDaughterProp(10) =  pidResonse->NumberOfSigmasTPC(negTrack,AliPID::kElectron);
  // fDaughterProp(23) =  pidResonse->NumberOfSigmasTPC(negTrack,AliPID::kPion);
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
  // fDaughterProp(14) =  (Float_t)nPosClusterITS;
  // fDaughterProp(15) =  (Float_t)nNegClusterITS;
  // if (nPosClusterITS > 0 ){
  //   fDaughterProp(16) =  posTrack->GetITSsignal();
  //   fDaughterProp(20) =  pidResonse->NumberOfSigmasITS(posTrack,AliPID::kElectron);
  // } else {
  //   fDaughterProp(16) =  1000;
  //   fDaughterProp(20) =  20;
  // }
  // if (nNegClusterITS > 0 ){
  //   fDaughterProp(17) =  negTrack->GetITSsignal();
  //   fDaughterProp(21) =  pidResonse->NumberOfSigmasITS(negTrack,AliPID::kElectron);
  // } else {
  //   fDaughterProp(17) =  1000;
  //   fDaughterProp(21) =  20;
  // }

  // TOF 
  // if((posTrack->GetStatus() & AliESDtrack::kTOFpid) && !(posTrack->GetStatus() & AliESDtrack::kTOFmismatch)){
  //   Double_t t0pos = pidResonse->GetTOFResponse().GetStartTime(posTrack->P());
  //   Double_t timesPos[9];
  //   posTrack->GetIntegratedTimes(timesPos,9);
  //   Double_t TOFsignalPos =	posTrack->GetTOFsignal();
  //   Double_t dTpos = TOFsignalPos - t0pos - timesPos[0];
  //   fDaughterProp(4) =  dTpos;
  //   fDaughterProp(5) =  pidResonse->NumberOfSigmasTOF(posTrack, AliPID::kElectron);
  // } else {
  //   fDaughterProp(4) =  20000;
  //   fDaughterProp(5) =  -20;
  // }
  // if((negTrack->GetStatus() & AliESDtrack::kTOFpid) && !(negTrack->GetStatus() & AliESDtrack::kTOFmismatch)){
  //   Double_t t0neg = pidResonse->GetTOFResponse().GetStartTime(negTrack->P());
  //   Double_t timesNeg[9];
  //   negTrack->GetIntegratedTimes(timesNeg,9);
  //   Double_t TOFsignalNeg =	negTrack->GetTOFsignal();
  //   Double_t dTneg = TOFsignalNeg - t0neg - timesNeg[0];
  //   fDaughterProp(11) =  dTneg;
  //   fDaughterProp(12) =  pidResonse->NumberOfSigmasTOF(negTrack, AliPID::kElectron);
  // } else {
  //   fDaughterProp(11) =  20000;
  //   fDaughterProp(12) =  -20;
  // }

  // fDaughterProp(6) =  (Float_t)posTrack->GetTPCClusterInfo(2,0,fConversionCuts->GetFirstTPCRow(gamma->GetConversionRadius()));
  // fDaughterProp(18) =  posTrack->GetNcls(1);
  // fDaughterProp(13) =  (Float_t)negTrack->GetTPCClusterInfo(2,0,fConversionCuts->GetFirstTPCRow(gamma->GetConversionRadius()));
  // fDaughterProp(19) =  negTrack->GetNcls(1);
  
  // if (fTreeQA){
  //   fTreeQA->Fill();
  // }

  fele_NfindableCls  =  (Float_t)posTrack->GetTPCClusterInfo(2,0,fConversionCuts->GetFirstTPCRow(gamma->GetConversionRadius()));
  fpos_Cls           =  posTrack->GetNcls(1);
  fpos_NfindableCls  =  (Float_t)negTrack->GetTPCClusterInfo(2,0,fConversionCuts->GetFirstTPCRow(gamma->GetConversionRadius()));

  fele_theta          =  negTrack->Theta();
  fele_pt             =  negTrack->Pt();
  fele_phi            =  negTrack->Phi();
  fele_nSigmaTPC      =  pidResonse->NumberOfSigmasTPC(negTrack,AliPID::kElectron);
  fele_nSigmaTPCpion  =  pidResonse->NumberOfSigmasTPC(negTrack,AliPID::kPion);
  fele_TPCsignal      =  negTrack->GetTPCsignal();
  fele_Cls            =  negTrack->GetNcls(1);

  fpos_theta          =  posTrack->Theta();
  fpos_pt             =  posTrack->Pt();
  fpos_phi            =  posTrack->Phi();
  fpos_nSigmaTPC      =  pidResonse->NumberOfSigmasTPC(posTrack,AliPID::kElectron);
  fpos_nSigmaTPCpion  =  pidResonse->NumberOfSigmasTPC(posTrack,AliPID::kPion);
  fpos_TPCsignal      =  posTrack->GetTPCsignal();
  fpos_Cls            =  posTrack->GetNcls(1);

 if(fInputEvent->IsA()==AliESDEvent::Class()){
    fele_SPD1 = (dynamic_cast<AliESDtrack*>(negTrack))->HasPointOnITSLayer(0);
    fele_SPD2 = (dynamic_cast<AliESDtrack*>(negTrack))->HasPointOnITSLayer(1);
    fele_SDD1 = (dynamic_cast<AliESDtrack*>(negTrack))->HasPointOnITSLayer(2);
    fele_SDD2 = (dynamic_cast<AliESDtrack*>(negTrack))->HasPointOnITSLayer(3);
    fele_SSD1 = (dynamic_cast<AliESDtrack*>(negTrack))->HasPointOnITSLayer(4);
    fele_SSD2 = (dynamic_cast<AliESDtrack*>(negTrack))->HasPointOnITSLayer(5);
    fpos_SPD1 = (dynamic_cast<AliESDtrack*>(posTrack))->HasPointOnITSLayer(0);
    fpos_SPD2 = (dynamic_cast<AliESDtrack*>(posTrack))->HasPointOnITSLayer(1);
    fpos_SDD1 = (dynamic_cast<AliESDtrack*>(posTrack))->HasPointOnITSLayer(2);
    fpos_SDD2 = (dynamic_cast<AliESDtrack*>(posTrack))->HasPointOnITSLayer(3);
    fpos_SSD1 = (dynamic_cast<AliESDtrack*>(posTrack))->HasPointOnITSLayer(4);
    fpos_SSD2 = (dynamic_cast<AliESDtrack*>(posTrack))->HasPointOnITSLayer(5);
  }else if ( fInputEvent->IsA()==AliAODEvent::Class()){
    fele_SPD1 = (dynamic_cast<AliAODTrack*>(negTrack))->HasPointOnITSLayer(0);
    fele_SPD2 = (dynamic_cast<AliAODTrack*>(negTrack))->HasPointOnITSLayer(1);
    fele_SDD1 = (dynamic_cast<AliAODTrack*>(negTrack))->HasPointOnITSLayer(2);
    fele_SDD2 = (dynamic_cast<AliAODTrack*>(negTrack))->HasPointOnITSLayer(3);
    fele_SSD1 = (dynamic_cast<AliAODTrack*>(negTrack))->HasPointOnITSLayer(4);
    fele_SSD2 = (dynamic_cast<AliAODTrack*>(negTrack))->HasPointOnITSLayer(5);
    fpos_SPD1 = (dynamic_cast<AliAODTrack*>(posTrack))->HasPointOnITSLayer(0);
    fpos_SPD2 = (dynamic_cast<AliAODTrack*>(posTrack))->HasPointOnITSLayer(1);
    fpos_SDD1 = (dynamic_cast<AliAODTrack*>(posTrack))->HasPointOnITSLayer(2);
    fpos_SDD2 = (dynamic_cast<AliAODTrack*>(posTrack))->HasPointOnITSLayer(3);
    fpos_SSD1 = (dynamic_cast<AliAODTrack*>(posTrack))->HasPointOnITSLayer(4);
    fpos_SSD2 = (dynamic_cast<AliAODTrack*>(posTrack))->HasPointOnITSLayer(5);    
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
    fele_TOFsignal = dTneg;
    fele_nSigmaTOF = pidResonse->NumberOfSigmasTOF(negTrack, AliPID::kElectron);
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
    fpos_TOFsignal  = dTpos;
    fpos_nSigmaTOF  = pidResonse->NumberOfSigmasTOF(posTrack, AliPID::kElectron);
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
