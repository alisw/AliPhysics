
///////////////////////////////////////////////////////////////////
//                                                               //            
// AliAnalysispp13TeVEMCalHFEMCReco.cxx                               //
// Author: Vivek Singh                                           //
//                                                               //
///////////////////////////////////////////////////////////////////

#include <Riostream.h>
#include "AliAnalysisUtils.h"

#include <vector>
#include "TMath.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "THnSparse.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TFile.h"
#include "TVector3.h"
#include "TNtuple.h"
#include <TRandom3.h>
#include "TProfile.h"
#include "TParticle.h"
#include <TDatabasePDG.h>
#include <TParticlePDG.h>
#include "TGeoManager.h"

#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDpid.h"

#include "AliCentrality.h"
#include "AliCFContainer.h"
#include "AliMagF.h"
#include "TGeoGlobalMagField.h"
#include "AliPhysicsSelection.h"

#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAODTrack.h"
#include "AliAnalysispp13TeVEMCalHFEMCReco.h"
#include "AliVVertex.h"
#include "AliVertexerTracks.h"


#include "AliAODCaloCluster.h"
#include "AliEMCALGeometry.h"


#include "AliTPCdEdxInfo.h"
#include "AliPID.h"
#include "AliAODPid.h"
#include "AliPIDResponse.h"
#include "AliTPCPIDResponse.h"
#include "AliTOFPIDResponse.h"
#include "AliHFEtools.h"
#include "AliHFEpid.h"
#include "AliHFEpidBase.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEcontainer.h"
#include "AliHFEcuts.h"
#include "AliHFEtools.h"
#include "AliTOFPIDResponse.h"
#include "AliTPCPIDResponse.h"
#include "AliHelperPID.h"
#include "AliHFEpidTPC.h"

#include "AliSelectNonHFE.h"
#include "AliHFEextraCuts.h"
#include "AliExternalTrackParam.h"


#include "AliGenEventHeader.h"
#include "AliEventPoolManager.h"
#include "AliAnalysisTaskSEImproveITS.h"

#include "AliKFVertex.h"
#include "AliKFParticle.h"

#include "AliMCEventHandler.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliMCEvent.h"

class AliAnalysispp13TeVEMCalHFEMCReco;
using namespace std;
ClassImp(AliAnalysispp13TeVEMCalHFEMCReco)

AliAnalysispp13TeVEMCalHFEMCReco::AliAnalysispp13TeVEMCalHFEMCReco(): AliAnalysisTaskSE(),

fIsMC(kFALSE),
fIsAOD(kTRUE),
ftrigger(AliVEvent::kINT7),
 // emcal correction
fUseTender(kTRUE),
// flag for emcal dcal
fFlagClsTypeEMC(kTRUE),
fFlagClsTypeDCAL(kTRUE),
// trigger events selection
fEMCEG1(kFALSE),
fEMCEG2(kFALSE),
fDCalDG1(kFALSE),
fDCalDG2(kFALSE),

fRecalIP(kTRUE),

fEtarange(0.6),
fTPCNCrRows(70),
fRatioCrossedRowOverFindable(0.8),
fITSNclus(3),
fTPCNclusPID(60),
fSPDBoth(kTRUE),
fSPDAny(kFALSE),
fSPDFirst(kFALSE),
fDCAxyCut(1),
fDCAzCut(2),
fTPCnsigmin(-1),
fTPCnsigmax(3),
fCutEopEMin(0.9),
fCutEopEMax(1.2),
fM02Min(0.05),
fM02Max1(0.9),
fM02Max2(0.7),
fM02Max3(0.5),

fInvmassCut(0.14),
fAssoTPCCluster(60),
fAssoITSRefit(kTRUE),
fAssopTMin(0.1),
fAssoEtarange(0.9),
fAssoTPCnsig(3.0),


fTenderClusterName("caloClusters"),
fTenderTrackName("tracks"),
fTracks_tender(0),
fCaloClusters_tender(0),

 // events
fAOD(0),
fOutputList(0), 
fHistEvent(0),
fNentries(0),

fHistVz(0),
fHistVzwc(0),

fHistMul(0),
fHistPt(0),

fHistEta(0),
fHistPhi(0),

fHistdca(0),

 //PID Cut
fPID(0),   
fPidResponse(0),
fHistBethe(0),
fnSigmaVsP_TPC(0),

fHistClustE(0),
fEMCClsEtaPhi(0), 
fHistoNCells(0),
fHistoTimeEMC(0),

fHistPtMatch(0),
fEMCTrkMatch(0),
fEMCClsEtaPhiTrkMatch(0),
fEMCTrkMatch_Phi(0),
fEMCTrkMatch_Eta(0),

fvalueElectron(0),
fSparseElectron(0),

//MC
fMCArray(0),
fMCHeader(0),
fMCparticle(0),
fPthfeGenerated(0),
fPthfe_rec(0),
fPthfe_rec_TrkSel(0),
fNTotMCpart(0),
fNpureMC(0),
fNembMCpi0(0),
fNembMCeta(0),

fMCparticleMother(0),
fMCparticleGMother(0),
fMCparticleGGMother(0),
fMCparticleGGGMother(0),

//Used in the function FindMother
fIsHFE1(kFALSE),
fIsHFE2(kFALSE),
fIsNonHFE(kFALSE),
fIsFromD(kFALSE),
fIsFromBarionB(kFALSE),
fIsFromMesonB(kFALSE),
fIsFromBarionBD(kFALSE),
fIsFromMesonBD(kFALSE),
fIsFromPi0(kFALSE),
fIsFromEta(kFALSE),
fIsFromGamma(kFALSE),

//EID Cuts
fTrkDCA(-999.0),

fEopNL_AftEID(0),

fHadEovpNL_AftEID(0),
fHadPt_AftEID(0),

fHadDCA(0),
fInclsElecPt(0),
fInclElecDCA(0),

fNElecInEvt(0),
fNEle(0),

fTPCnSigma(-999.0),

fTPCnSigmaHadMin(-10),
fTPCnSigmaHadMax(-3.5),

fInvmassULSPt(0),
fInvmassLSPt(0),
fCalculateNonHFEEffi(1),

fULSElecPt(0),
fLSElecPt(0),
fULSElecDCA(0),
fLSElecDCA(0),

//nonhfe efficiency

fIsFrmEmbPi0(kFALSE),
fIsFrmEmbEta(kFALSE),
ftype(-1),
fWeight(1),
fWeightPi0(1),
fWeightEta(1),

fPi0Weight(0),
fEtaWeight(0),

fCalculateWeight(kFALSE), 
fSprsPi0EtaWeightCal(0),
fPi0EtaSpectraSp(0),
pi0MC(0),
etaMC(0),
gammaMC(0),

fRealInclsElecPt(0),
fNonHFeTrkPt(0),
fNonHFeEmbTrkPt(0),
fNonHFeEmbWeightTrkPt(0),
fPi0eEmbWeightTrkPt(0),
fEtaeEmbWeightTrkPt(0),

fRecoNonHFeTrkPt(0),
fRecoNonHFeEmbTrkPt(0),
fRecoNonHFeEmbWeightTrkPt(0),
fRecoPi0eEmbWeightTrkPt(0),
fRecoEtaeEmbWeightTrkPt(0),

fNonHFePairInvmassLS(0),
fNonHFePairInvmassULS(0),
fNonHFeEmbInvmassLS(0),
fNonHFeEmbInvmassULS(0),
fNonHFeEmbWeightInvmassLS(0),
fNonHFeEmbWeightInvmassULS(0),
fPi0EmbInvmassLS(0),
fPi0EmbInvmassULS(0),
fPi0EmbWeightInvmassLS(0),
fPi0EmbWeightInvmassULS(0),
fEtaEmbInvmassLS(0),
fEtaEmbInvmassULS(0),
fEtaEmbWeightInvmassLS(0),
fEtaEmbWeightInvmassULS(0),

fRecoLSeEmbTrkPt(0),
fRecoLSeEmbWeightTrkPt(0),
fRecoPi0LSeEmbWeightTrkPt(0),
fRecoEtaLSeEmbWeightTrkPt(0),
fRecoULSeEmbTrkPt(0),
fRecoULSeEmbWeightTrkPt(0),
fRecoPi0ULSeEmbWeightTrkPt(0),
fRecoEtaULSeEmbWeightTrkPt(0),

fHadConvRadius(0),
fIncleConvRadius(0),
fNonHFeConvRadius(0),
fHFeConvRadius(0),

fNonHFeEmbTrkRConv(0),
fPi0eEmbWeightTrkRConv(0),
fNonHFeEmbWeightTrkRConv(0),
fEtaeEmbWeightTrkRConv(0),

fRecoNonHFeEmbRConv(0),
fRecoPi0eEmbWeightTrkRConv(0),
fRecoNonHFeEmbWeightTrkRConv(0),
fRecoEtaeEmbWeightTrkRConv(0),

fRVsULSElecPt(0),
fRVsLSElecPt(0),

fnBinsDCAHisto(4000),

fCalculateMCTemplWeightCalc(kFALSE),
fFillMCTemplates(kFALSE),

fBHadpT(0),
fBMesonpT(0),
fBDHadpT(0),
fDHadpT(0),
fDMesonpT(0),
fD0pT(0),
fDPluspT(0),
fDspT(0),
fLambdaCpT(0),

fDcent(0),
fDUp(0),
fDDown(0),
fBcent(0),
fBMin(0),
fBMax(0),
fD0(0),
fDPlus(0),
fDs(0),
fLc(0),
fB(0),
fWeightB(1.0),
fWeightBMin(1.0),
fWeightBMax(1.0),
fWeightD(1.0),
fWeightDUp(1.0),
fWeightDDown(1.0),

fDElecDCA(0),
fBElecDCA(0),
fBHadElecDCA(0),
fBMesonElecDCA(0),
fBBaryonElecDCA(0),
fDHadElecDCA(0),
fDMesonElecDCA(0),
fDBaryonElecDCA(0),
fLambdaCElecDCA(0),
fD0ElecDCA(0),
fSprsTemplatesNoWeight(0),
fSprsTemplatesWeight(0),
fSprsTemplatesWeightVar1(0),
fSprsTemplatesWeightVar2(0),

fCalculateElecRecoEffi(kFALSE),

fInclElePhysPriAll(0),
fHFEPhysPriAll(0),
fBEPhysPriAll(0),
fDEPhysPriAll(0),

fInclElePhysPriTrkCuts(0),
fHFEPhysPriTrkCuts(0),
fBEPhysPriTrkCuts(0),
fDEPhysPriTrkCuts(0),

fInclElePhysPriEMCMatch(0),
fHFEPhysPriEMCMatch(0),
fBEPhysPriEMCMatch(0),
fDEPhysPriEMCMatch(0),

fInclElePhysPriEovP(0),
fHFEPhysPriEovP(0),
fBEPhysPriEovP(0),
fDEPhysPriEovP(0),

fInclElePhysPriTPCnsig(0),
fHFEPhysPriTPCnsig(0),
fBEPhysPriTPCnsig(0),
fDEPhysPriTPCnsig(0),

fInclElePhysPriSS(0),
fHFEPhysPriSS(0),
fBEPhysPriSS(0),
fDEPhysPriSS(0),

fExtraCuts(0)



{
fPID = new AliHFEpid("hfePid");
fvalueElectron = new Double_t[6];
//fvalueRadius = new Double_t[4];
}
  

//________________________________________________________________________

AliAnalysispp13TeVEMCalHFEMCReco::AliAnalysispp13TeVEMCalHFEMCReco(const char *name)   : AliAnalysisTaskSE(name), 

fIsMC(kFALSE),
fIsAOD(kTRUE),
ftrigger(AliVEvent::kINT7),
 // emcal correction
fUseTender(kTRUE),
// flag for emcal dcal
fFlagClsTypeEMC(kTRUE),
fFlagClsTypeDCAL(kTRUE),
// trigger events selection
fEMCEG1(kFALSE),
fEMCEG2(kFALSE),
fDCalDG1(kFALSE),
fDCalDG2(kFALSE),

fRecalIP(kTRUE),

fEtarange(0.6),
fTPCNCrRows(70),
fRatioCrossedRowOverFindable(0.8),
fITSNclus(3),
fTPCNclusPID(60),
fSPDBoth(kTRUE),
fSPDAny(kFALSE),
fSPDFirst(kFALSE),
fDCAxyCut(1),
fDCAzCut(2),
fTPCnsigmin(-1),
fTPCnsigmax(3),
fCutEopEMin(0.9),
fCutEopEMax(1.2),
fM02Min(0.05),
fM02Max1(0.9),
fM02Max2(0.7),
fM02Max3(0.5),

fInvmassCut(0.14),
fAssoTPCCluster(60),
fAssoITSRefit(kTRUE),
fAssopTMin(0.1),
fAssoEtarange(0.9),
fAssoTPCnsig(3.0),


fTenderClusterName("caloClusters"),
fTenderTrackName("tracks"),
fTracks_tender(0),
fCaloClusters_tender(0),

 // events
fAOD(0),
fOutputList(0), 
fHistEvent(0),
fNentries(0),

fHistVz(0),
fHistVzwc(0),

fHistMul(0),
fHistPt(0),

fHistEta(0),
fHistPhi(0),

fHistdca(0),

 //PID Cut
fPID(0),   
fPidResponse(0),
fHistBethe(0),
fnSigmaVsP_TPC(0),

fHistClustE(0),
fEMCClsEtaPhi(0), 
fHistoNCells(0),
fHistoTimeEMC(0),

fHistPtMatch(0),
fEMCTrkMatch(0),
fEMCClsEtaPhiTrkMatch(0),
fEMCTrkMatch_Phi(0),
fEMCTrkMatch_Eta(0),

fvalueElectron(0),
fSparseElectron(0),

//MC
fMCArray(0),
fMCHeader(0),
fMCparticle(0),
fPthfeGenerated(0),
fPthfe_rec(0),
fPthfe_rec_TrkSel(0),
fNTotMCpart(0),
fNpureMC(0),
fNembMCpi0(0),
fNembMCeta(0),

fMCparticleMother(0),
fMCparticleGMother(0),
fMCparticleGGMother(0),
fMCparticleGGGMother(0),

//Used in the function FindMother
fIsHFE1(kFALSE),
fIsHFE2(kFALSE),
fIsNonHFE(kFALSE),
fIsFromD(kFALSE),
fIsFromBarionB(kFALSE),
fIsFromMesonB(kFALSE),
fIsFromBarionBD(kFALSE),
fIsFromMesonBD(kFALSE),
fIsFromPi0(kFALSE),
fIsFromEta(kFALSE),
fIsFromGamma(kFALSE),

//EID Cuts
fTrkDCA(-999.0),

fEopNL_AftEID(0),

fHadEovpNL_AftEID(0),
fHadPt_AftEID(0),

fHadDCA(0),
fInclsElecPt(0),
fInclElecDCA(0),

fNElecInEvt(0),
fNEle(0),

fTPCnSigma(-999.0),

fTPCnSigmaHadMin(-10),
fTPCnSigmaHadMax(-3.5),

fInvmassULSPt(0),
fInvmassLSPt(0),
fCalculateNonHFEEffi(1),

fULSElecPt(0),
fLSElecPt(0),
fULSElecDCA(0),
fLSElecDCA(0),

//nonhfe efficiency

fIsFrmEmbPi0(kFALSE),
fIsFrmEmbEta(kFALSE),
ftype(-1),
fWeight(1),
fWeightPi0(1),
fWeightEta(1),

fPi0Weight(0),
fEtaWeight(0),

fCalculateWeight(kFALSE), 
fSprsPi0EtaWeightCal(0),
fPi0EtaSpectraSp(0),
pi0MC(0),
etaMC(0),
gammaMC(0),

fRealInclsElecPt(0),
fNonHFeTrkPt(0),
fNonHFeEmbTrkPt(0),
fNonHFeEmbWeightTrkPt(0),
fPi0eEmbWeightTrkPt(0),
fEtaeEmbWeightTrkPt(0),

fRecoNonHFeTrkPt(0),
fRecoNonHFeEmbTrkPt(0),
fRecoNonHFeEmbWeightTrkPt(0),
fRecoPi0eEmbWeightTrkPt(0),
fRecoEtaeEmbWeightTrkPt(0),

fNonHFePairInvmassLS(0),
fNonHFePairInvmassULS(0),
fNonHFeEmbInvmassLS(0),
fNonHFeEmbInvmassULS(0),
fNonHFeEmbWeightInvmassLS(0),
fNonHFeEmbWeightInvmassULS(0),
fPi0EmbInvmassLS(0),
fPi0EmbInvmassULS(0),
fPi0EmbWeightInvmassLS(0),
fPi0EmbWeightInvmassULS(0),
fEtaEmbInvmassLS(0),
fEtaEmbInvmassULS(0),
fEtaEmbWeightInvmassLS(0),
fEtaEmbWeightInvmassULS(0),

fRecoLSeEmbTrkPt(0),
fRecoLSeEmbWeightTrkPt(0),
fRecoPi0LSeEmbWeightTrkPt(0),
fRecoEtaLSeEmbWeightTrkPt(0),
fRecoULSeEmbTrkPt(0),
fRecoULSeEmbWeightTrkPt(0),
fRecoPi0ULSeEmbWeightTrkPt(0),
fRecoEtaULSeEmbWeightTrkPt(0),

fHadConvRadius(0),
fIncleConvRadius(0),
fNonHFeConvRadius(0),
fHFeConvRadius(0),

fNonHFeEmbTrkRConv(0),
fPi0eEmbWeightTrkRConv(0),
fNonHFeEmbWeightTrkRConv(0),
fEtaeEmbWeightTrkRConv(0),

fRecoNonHFeEmbRConv(0),
fRecoPi0eEmbWeightTrkRConv(0),
fRecoNonHFeEmbWeightTrkRConv(0),
fRecoEtaeEmbWeightTrkRConv(0),

fRVsULSElecPt(0),
fRVsLSElecPt(0),

fnBinsDCAHisto(4000),

fCalculateMCTemplWeightCalc(kFALSE),
fFillMCTemplates(kFALSE),

fBHadpT(0),
fBMesonpT(0),
fBDHadpT(0),
fDHadpT(0),
fDMesonpT(0),
fD0pT(0),
fDPluspT(0),
fDspT(0),
fLambdaCpT(0),

fDcent(0),
fDUp(0),
fDDown(0),
fBcent(0),
fBMin(0),
fBMax(0),
fD0(0),
fDPlus(0),
fDs(0),
fLc(0),
fB(0),
fWeightB(1.0),
fWeightBMin(1.0),
fWeightBMax(1.0),
fWeightD(1.0),
fWeightDUp(1.0),
fWeightDDown(1.0),

fDElecDCA(0),
fBElecDCA(0),
fBHadElecDCA(0),
fBMesonElecDCA(0),
fBBaryonElecDCA(0),
fDHadElecDCA(0),
fDMesonElecDCA(0),
fDBaryonElecDCA(0),
fLambdaCElecDCA(0),
fD0ElecDCA(0),
fSprsTemplatesNoWeight(0),
fSprsTemplatesWeight(0),
fSprsTemplatesWeightVar1(0),
fSprsTemplatesWeightVar2(0),

fCalculateElecRecoEffi(kFALSE),

fInclElePhysPriAll(0),
fHFEPhysPriAll(0),
fBEPhysPriAll(0),
fDEPhysPriAll(0),

fInclElePhysPriTrkCuts(0),
fHFEPhysPriTrkCuts(0),
fBEPhysPriTrkCuts(0),
fDEPhysPriTrkCuts(0),

fInclElePhysPriEMCMatch(0),
fHFEPhysPriEMCMatch(0),
fBEPhysPriEMCMatch(0),
fDEPhysPriEMCMatch(0),

fInclElePhysPriEovP(0),
fHFEPhysPriEovP(0),
fBEPhysPriEovP(0),
fDEPhysPriEovP(0),

fInclElePhysPriTPCnsig(0),
fHFEPhysPriTPCnsig(0),
fBEPhysPriTPCnsig(0),
fDEPhysPriTPCnsig(0),

fInclElePhysPriSS(0),
fHFEPhysPriSS(0),
fBEPhysPriSS(0),
fDEPhysPriSS(0),

fExtraCuts(0)


{
  // Constructor
  fPID = new AliHFEpid("hfePid");
  fvalueElectron = new Double_t[6];

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());  
  DefineOutput(2, TH1F::Class());
}

//______________________________________________________________________
AliAnalysispp13TeVEMCalHFEMCReco::~AliAnalysispp13TeVEMCalHFEMCReco()
{
    //Destructor

    delete fPID;
    delete fTracks_tender;
    delete fCaloClusters_tender;
    delete []fvalueElectron;
    delete fSparseElectron;
    delete fSprsPi0EtaWeightCal;
    delete fSprsTemplatesNoWeight;
    delete fSprsTemplatesWeight;
    delete fSprsTemplatesWeightVar1;
    delete fSprsTemplatesWeightVar2;

    if(fDcent) {delete fDcent; fDcent=0;}
    if(fDUp)   {delete fDUp; fDUp=0;}
    if(fDDown) {delete fDDown; fDDown=0;}
    if(fBcent) {delete fBcent; fBcent=0;}
    if(fBMin)  {delete fBMin; fBMin=0;}
    if(fBMax)  {delete fBMax; fBMax=0;}
    if(fD0)    {delete fD0; fD0=0;}
    if(fDPlus) {delete fDPlus; fDPlus=0;}
    if(fDs)    {delete fDs; fDs=0;}
    if(fLc)    {delete fLc; fLc=0;}
    if(fB)    {delete fB; fB=0;}

    if(fOutputList) { delete fOutputList; fOutputList = 0;}
    if (fNentries){ delete fNentries; fNentries = 0;}
   
}
//________________________________________________________________________
void AliAnalysispp13TeVEMCalHFEMCReco::Init()
{
  // Initialization 
  if(fDebug > 1) printf("AliAnalysispp13TeVEMCalHFEMCReco::Init() \n");
  return;
}


void AliAnalysispp13TeVEMCalHFEMCReco::UserCreateOutputObjects()
{

  AliDebug(3, "Creating Output Objects");

  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  Double_t pi = TMath::Pi();
  fPi0Weight = new TF1("fPi0Weight","[0] / TMath::Power(TMath::Exp(-[1]*x - [2]*x*x) + x/[3], [4])");
  fEtaWeight = new TF1("fEtaWeight","[0] / TMath::Power(TMath::Exp(-[1]*x - [2]*x*x) + x/[3], [4])");

  fPi0Weight->SetParameters(6.16962e+02,-3.55899e-02,4.67347e-03,1.56646e+00,5.55130e+00);
  fEtaWeight->SetParameters(3.25021e+02,-6.77106e-02,4.16408e-03,2.29748e+00,6.03883e+00);
  fOutputList->Add(fPi0Weight);
  fOutputList->Add(fEtaWeight);

    fHistEvent=new TH1F("fHistEvent","",20,0,20);
    fOutputList->Add(fHistEvent);
    fHistEvent->GetYaxis()->SetTitle("Counts");
    fHistEvent->GetXaxis()->SetTitle("Event cuts");
    fHistEvent->GetXaxis()->SetBinLabel(2,"All");
    fHistEvent->GetXaxis()->SetBinLabel(4," Trigger ");
    fHistEvent->GetXaxis()->SetBinLabel(6," EMCal acpt  ");
    fHistEvent->GetXaxis()->SetBinLabel(8," N_{cont} > 2 ");
    fHistEvent->GetXaxis()->SetBinLabel(10," SPD pile up ");
    fHistEvent->GetXaxis()->SetBinLabel(12," Multi-Vtx pile up ");
    fHistEvent->GetXaxis()->SetBinLabel(14,"Vtx_{z}<10cm");
    fHistEvent->GetXaxis()->SetBinLabel(16," Analysed ");
    fHistEvent->GetXaxis()->SetNdivisions(1,kFALSE);
    fHistEvent->SetMinimum(0);

  fNentries=new TH1F("CutSet", "", 33,-1.5,32.5);
  fNentries->GetXaxis()->SetBinLabel(0,"trigger");
  fNentries->GetXaxis()->SetBinLabel(1,"ITSNclus");
  fNentries->GetXaxis()->SetBinLabel(2,"TPCNCrRows");
  fNentries->GetXaxis()->SetBinLabel(3,"TPCNclusPID");
  fNentries->GetXaxis()->SetBinLabel(4,"RatioCrossedRowOverFindable");
  fNentries->GetXaxis()->SetBinLabel(5,"SPDBoth");
  fNentries->GetXaxis()->SetBinLabel(6,"SPDAny");
  fNentries->GetXaxis()->SetBinLabel(7,"SPDFirst");
  fNentries->GetXaxis()->SetBinLabel(8,"DCAxyCut");
  fNentries->GetXaxis()->SetBinLabel(9,"DCAzCut");
  fNentries->GetXaxis()->SetBinLabel(10,"Etarange");  
  fNentries->GetXaxis()->SetBinLabel(11,"TPCnsigmin");
  fNentries->GetXaxis()->SetBinLabel(12,"TPCnsigmax");
  fNentries->GetXaxis()->SetBinLabel(13,"InvmassCut");
  fNentries->GetXaxis()->SetBinLabel(14,"AssoTPCCluster");  
  fNentries->GetXaxis()->SetBinLabel(15,"AssoITSRefit");
  fNentries->GetXaxis()->SetBinLabel(16,"AssopTMin");
  fNentries->GetXaxis()->SetBinLabel(17,"AssoEtarange");
  fNentries->GetXaxis()->SetBinLabel(18,"AssoTPCnsig");
  fNentries->GetXaxis()->SetBinLabel(19,"fEMCEG1");
  fNentries->GetXaxis()->SetBinLabel(20,"fDCalDG1");
  fNentries->GetXaxis()->SetBinLabel(21,"fEMCEG2");
  fNentries->GetXaxis()->SetBinLabel(22,"fDCalDG2");
  fNentries->GetXaxis()->SetBinLabel(23,"fUseTender");
  fNentries->GetXaxis()->SetBinLabel(24,"CutEopEMin");
  fNentries->GetXaxis()->SetBinLabel(25,"CutEopEMax");
  fNentries->GetXaxis()->SetBinLabel(26,"fM02Min");
  fNentries->GetXaxis()->SetBinLabel(27,"fM02Max1");
  fNentries->GetXaxis()->SetBinLabel(28,"fM02Max2");
  fNentries->GetXaxis()->SetBinLabel(29,"fM02Max3");
  fOutputList->Add(fNentries);
  
  fNentries->SetBinContent(0,ftrigger);
  fNentries->SetBinContent(1,fITSNclus);
  fNentries->SetBinContent(2,fTPCNCrRows);
  fNentries->SetBinContent(3,fTPCNclusPID);
  fNentries->SetBinContent(4,fRatioCrossedRowOverFindable);
  fNentries->SetBinContent(5,fSPDBoth);
  fNentries->SetBinContent(6,fSPDAny);
  fNentries->SetBinContent(7,fSPDFirst);
  fNentries->SetBinContent(8,fDCAxyCut);
  fNentries->SetBinContent(9,fDCAzCut);
  fNentries->SetBinContent(10,fEtarange);  
  fNentries->SetBinContent(11,fTPCnsigmin);
  fNentries->SetBinContent(12,fTPCnsigmax);
  fNentries->SetBinContent(13,fInvmassCut);
  fNentries->SetBinContent(14,fAssoTPCCluster);  
  fNentries->SetBinContent(15,fAssoITSRefit);
  fNentries->SetBinContent(16,fAssopTMin);
  fNentries->SetBinContent(17,fAssoEtarange);
  fNentries->SetBinContent(18,fAssoTPCnsig);
  fNentries->SetBinContent(19,fEMCEG1);
  fNentries->SetBinContent(20,fDCalDG1);
  fNentries->SetBinContent(21,fEMCEG2);
  fNentries->SetBinContent(22,fDCalDG2);
  fNentries->SetBinContent(23,fUseTender);
  fNentries->SetBinContent(24,fCutEopEMin);
  fNentries->SetBinContent(25,fCutEopEMax);
  fNentries->SetBinContent(26,fM02Min);
  fNentries->SetBinContent(27,fM02Max1);
  fNentries->SetBinContent(28,fM02Max2);
  fNentries->SetBinContent(29,fM02Max3); 
  fNentries->GetXaxis()->SetNdivisions(1,kFALSE);
  fNentries->Sumw2();
  fNentries->SetMinimum(0);


    fHistVz=new TH1F("fHistVz","Z_{vtx} Posistion before cut;Z_{vtx};Counts",400,-20,20);
    fHistVz->Sumw2();
    fOutputList->Add(fHistVz);

    fHistVzwc=new TH1F("fHistVzwc","Z_{vtx} Posistion after cut;Z_{vtx};Counts",400,-20,20);
    fHistVzwc->Sumw2();
    fOutputList->Add(fHistVzwc);
  
    fHistMul=new TH1F("fHistMul","Global tracks multiplicity of all charged particles;N_{ch};Counts",500,0,500);
    fOutputList->Add(fHistMul);

    fHistPt = new TH1F("fHistPt", "P_{T} distribution of global tracks;#it{p}_{T}(GeV/#it{c});Counts", 500, 0., 100.);
    fHistPt->Sumw2();
    fOutputList->Add(fHistPt);   

    fHistEta = new TH1F("fHistEta", "Eta distribution of global tracks;#eta;Counts", 300, -1.5, 1.5);
    fHistEta->Sumw2();
    fOutputList->Add(fHistEta);
 
    fHistPhi = new TH1F("fHistPhi", "Phi distribution of global tracks;#phi;Counts", 100,0, 2*TMath::Pi() ); 
    fHistPhi->Sumw2();
    fOutputList->Add(fHistPhi);       

    fHistdca=new TH2F("fHistdca","DCA of global tracks;DCA_{xy};DCA_{z}",200,-10,10,200,-10,10);
    fHistdca->Sumw2();
    fOutputList->Add(fHistdca);

    fHistBethe=new TH2F("fHistBethe","Particle identification by Energy loss;#it{p}(GeV/#it{c});TPC dE/dx(arb. Units)",300,0.,15.,750,10,160);
    fHistBethe->Sumw2();
    fOutputList->Add(fHistBethe);
   
    fnSigmaVsP_TPC= new TH2F("fnSigmaVsP_TPC", "fnSigmaVsP_TPC distribution;#it{p}(GeV/#it{c});n#sigma^{TPC}",300,0.,15.,750,-15.,15.);
    fnSigmaVsP_TPC->Sumw2();
    fOutputList->Add(fnSigmaVsP_TPC);

//EMCAL before track-cluster matching

    fHistClustE = new TH1F("fHistClustE", "Cluster energy distribution before track match ; Cluster E;Counts", 500, 0.0, 50.0);
    fHistClustE->GetSumw2();
    fOutputList->Add(fHistClustE);

    fEMCClsEtaPhi = new TH2F("fEMCClsEtaPhi","Cluster #eta and #phi distribution before track match;#eta;#phi",100,-0.9,0.9,200,0,6.3);
    fEMCClsEtaPhi->Sumw2();
    fOutputList->Add(fEMCClsEtaPhi);

    fHistoNCells = new TH2F("fHistoNCells","No of cells in a cluster;Cluster E;N^{EMC}_{cells}",500,0,50,30,0,30);
    fOutputList->Add(fHistoNCells);

    fHistoTimeEMC = new TH2F("fHistoTimeEMC"," Time;E (GeV); t(ns)",500,0,50,1800,-900,900);
    fOutputList->Add(fHistoTimeEMC);

// After Track Match

    fHistPtMatch = new TH1F("fHistPtMatch", "p_{T} distribution of tracks matched to EMCAL;#it{p}_{T} (GeV/#it{c});Counts",500, 0., 100.);
    fHistPtMatch->Sumw2();
    fOutputList->Add(fHistPtMatch);
    
    fEMCTrkMatch = new TH2F("fEMCTrkMatch","Distance of the cluster to its closest track in (#Delta#eta,#Delta#phi)plane;#Delta#eta;#Delta#phi",600,-0.3,0.3,600,-0.3,0.3);
    fOutputList->Add(fEMCTrkMatch);

    fEMCClsEtaPhiTrkMatch = new TH2F("fEMCClsEtaPhiTrkMatch","Cluster #eta and #phi distribution;#eta;#phi",100,-0.9,0.9,200,0,6.3);
    fOutputList->Add(fEMCClsEtaPhiTrkMatch);
        
    fEMCTrkMatch_Phi = new TH2F("fEMCTrkMatch_Phi","Distance of the cluster to its closest track in #Delta#phi vs p_{T};#it{p}(GeV/#it{c});#Delta#phi",500,0,100.0,600,-0.3,0.3);
    fOutputList->Add(fEMCTrkMatch_Phi);    

    fEMCTrkMatch_Eta = new TH2F("fEMCTrkMatch_Eta","Distance of the cluster to its closest track in #Delta#eta vs p_{T};#it{p}(GeV/#it{c});#Delta#eta",500,0,100.0,600,-0.3,0.3);
    fOutputList->Add(fEMCTrkMatch_Eta);


//THnSparse
  Int_t bins[6]   =       {500, 200, 200, 200, 200,500};
  Double_t xmin[6]  = {  0,  -10,   0,   0,   0,    0 };
  Double_t xmax[6]  = {  100,   10,   2,   2,  2,    50};

  fSparseElectron   = new THnSparseD ("Electron","Electron;#it{p}(GeV/#it{c});n#sigma;E/#it{p}(GeV/#it{c});M02;M20;Cluster Energy(GeV);",6 ,bins,xmin,xmax);
  fSparseElectron->GetAxis(0)->SetName("pT");     
  fSparseElectron->GetAxis(1)->SetName("nSigma");
  fSparseElectron->GetAxis(2)->SetName("E/p");
  fSparseElectron->GetAxis(3)->SetName("M02");
  fSparseElectron->GetAxis(4)->SetName("M20");
  fSparseElectron->GetAxis(5)->SetName("ClustE");
  fSparseElectron->Sumw2();
  fOutputList->Add(fSparseElectron); 

  //===================================================================
  //===================  Non HFe ======================================
  //===================================================================

  fHadPt_AftEID = new TH1F("fHadPt_AftEID","#it{p}_{T} distribution of hadrons after Eid cuts;#it{p}_{T} (GeV/#it{c});counts",250,0,50);
  fHadPt_AftEID->Sumw2();
  fOutputList->Add(fHadPt_AftEID);

  fHadEovpNL_AftEID = new TH2F("fHadEovpNL_AftEID", "E/p distribution for hadrons -10<nsig<-3.5, NonLinearE, SS cuts;p_{T} (GeV/c);E/p", 100,0,50,200, 0.0, 2.0);
  fHadEovpNL_AftEID->Sumw2();
  fOutputList->Add(fHadEovpNL_AftEID);

  fHadDCA = new TH2F("fHadDCA","Hadron DCA; #it{p}_{T}(GeV/#it{c}); DCAxMagFieldxSign; counts;", 250,0,50., 400,-0.4,0.4);
  fHadDCA->Sumw2();
  fOutputList->Add(fHadDCA);

  fInclsElecPt = new TH1F("fInclsElecPt","#it{p}_{T} distribution of inclusive electrons;#it{p}_{T} (GeV/#it{c});counts",250,0,50);
  fInclsElecPt->Sumw2();
  fOutputList->Add(fInclsElecPt);

  fInclElecDCA = new TH2F("fInclElecDCA","Inclusive electron DCA; #it{p}_{T}(GeV/#it{c}); DCAxMagFieldxSign; counts;", 250,0,50., 400,-0.4,0.4);
  fInclElecDCA->Sumw2();
  fOutputList->Add(fInclElecDCA);

    
  fEopNL_AftEID = new TH2F("fEopNL_AftEID", "E/p distribution after nsig, SS cuts, NonLinearE;p_{T} (GeV/c);E/p", 100,0,50,200, 0.0, 2.0);
  fEopNL_AftEID->Sumw2();
  fOutputList->Add(fEopNL_AftEID);

  fNElecInEvt = new TH1F("fNElecInEvt","No of electrons in the event; N^{ele};counts",20,-0.5,19.5);
  fOutputList->Add(fNElecInEvt);

  fInvmassULSPt = new TH2F("fInvmassULSPt", "Invmass of ULS (e,e) for #it{p}_T^{e}>1; #it{p}_{T}(GeV/#it{c}); mass(GeV/#it{c}^2); counts;", 250,0,50,500,0.,1.0);
  fInvmassULSPt->Sumw2();
  fOutputList->Add(fInvmassULSPt);

  fInvmassLSPt = new TH2F("fInvmassLSPt", "Invmass of LS (e,e) for pt^{e}>1; #it{p}_{T}(GeV/c); mass(GeV/#it{c}^2); counts;", 250,0,50,500,0,1.0);
  fInvmassLSPt->Sumw2();
  fOutputList->Add(fInvmassLSPt);
    
  fULSElecPt  = new TH1F("fULSElecPt","#it{p}_{T} distribution of ULS electrons;#it{p}_{T} (GeV/#it{c});counts",250,0,50);
  fULSElecPt->Sumw2();
  fOutputList->Add(fULSElecPt);
    
  fLSElecPt= new TH1F("fLSElecPt","#it{p}_{T} distribution of LS electrons;#it{p}_{T} (GeV/#it{c});counts",250,0,50);
  fLSElecPt->Sumw2();
  fOutputList->Add(fLSElecPt);
  
    
  fULSElecDCA = new TH2F("fULSElecDCA","ULS electron DCA; #it{p}_{T}(GeV/#it{c}); DCAxMagFieldxSign; counts;", 250,0,50., 400,-0.4,0.4);
  fULSElecDCA->Sumw2();
  fOutputList->Add(fULSElecDCA);
    
  fLSElecDCA = new TH2F("fLSElecDCA","LS electron DCA; #it{p}_{T}(GeV/#it{c}); DCAxMagFieldxSign; counts;", 250,0,50., 400,-0.4,0.4);
  fLSElecDCA->Sumw2();
  fOutputList->Add(fLSElecDCA);

 //+++++++++++++++++++++++++++++++++++++++++MC++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if(fIsMC)
{

  fPthfeGenerated = new TH1F("fPthfeGenerated","; p_{T} [GeV/c]; Counts",100,0,10);
  fOutputList->Add(fPthfeGenerated);

  fPthfe_rec = new TH1F("fPthfe_rec","; p_{T} [GeV/c]; Counts",100,0,10);   
  fOutputList->Add(fPthfe_rec);

  fPthfe_rec_TrkSel = new TH1F("fPthfe_rec_TrkSel","; p_{T} [GeV/c]; Counts",100,0,10);   
  fOutputList->Add(fPthfe_rec_TrkSel);

//nonhfe tagging efficiency
  if(fCalculateWeight){
  Int_t bin[5] =     {250,30,2,10,100}; //pT, PDG, EnhancedSigOrNot, pi0etaType, Radius
  Double_t xminWt[5] = {0,0,0,-1,0};
  Double_t xmaxWt[5] = {50,3,2,9,10};
  
  fSprsPi0EtaWeightCal = new THnSparseD("fSprsPi0EtaWeightCal","Sparse to calculate #pi^{0} and #eta weight;p_{T};PDG ID;EnhanceSigOrNot;pi0etaType;",5,bin,xminWt,xmaxWt);
  fSprsPi0EtaWeightCal->GetAxis(0)->SetName("pT");     
  fSprsPi0EtaWeightCal->GetAxis(1)->SetName("PDG");
  fSprsPi0EtaWeightCal->GetAxis(2)->SetName("EnhancedSigOrNot");
  fSprsPi0EtaWeightCal->GetAxis(3)->SetName("pi0etaType");
  fSprsPi0EtaWeightCal->GetAxis(4)->SetName("Radius");
  fSprsPi0EtaWeightCal->Sumw2();
  fOutputList->Add(fSprsPi0EtaWeightCal);
 
  Int_t nbinspt[3] = {400, 30, 10};
  Double_t binlow[3] = {0., 0, -1.};
  Double_t binup[3] = {40, 3, 9};
  fPi0EtaSpectraSp = new THnSparseF("fPi0EtaSpectraSp", "fPi0EtaSpectraSp;pt;source;type;", 3, nbinspt, binlow,binup);
  fPi0EtaSpectraSp->GetAxis(0)->SetName("pT"); 
  fPi0EtaSpectraSp->GetAxis(1)->SetName("source");
  fPi0EtaSpectraSp->GetAxis(2)->SetName("type");
  fPi0EtaSpectraSp->Sumw2();
  fOutputList->Add(fPi0EtaSpectraSp);

  pi0MC= new TH1F("pi0MC",";p_{t} (GeV/c)",250,0,50.);
  pi0MC->Sumw2();
  fOutputList->Add(pi0MC);
  
  etaMC= new TH1F("etaMC",";p_{t} (GeV/c)",250,0,50.);
  etaMC->Sumw2();
  fOutputList->Add(etaMC);
  
  gammaMC= new TH1F("gammaMC",";p_{t} (GeV/c)",250,0,50.);
  gammaMC->Sumw2();
  fOutputList->Add(gammaMC);
  }
  
  if(fCalculateNonHFEEffi){
  fRealInclsElecPt = new TH1F("fRealInclsElecPt","p_{T} distribution of MC tagged inclusive electrons;p_{T} (GeV/c);counts",250,0,50);
  fOutputList->Add(fRealInclsElecPt);

  fNonHFeTrkPt = new TH1F("fNonHFeTrkPt","Non-HF electrons from all generators;p_{T} (GeV/c);counts",250,0,50);
  fNonHFeTrkPt->Sumw2();
  fOutputList->Add(fNonHFeTrkPt);

  fNonHFeEmbTrkPt = new TH1F("fNonHFeEmbTrkPt","Non-HF electrons from embedded #pi^{0} and #eta + No mom;p_{T} (GeV/c);counts",250,0,50);
  fNonHFeEmbTrkPt->Sumw2();
  fOutputList->Add(fNonHFeEmbTrkPt);
        
  fNonHFeEmbWeightTrkPt = new TH1F("fNonHFeEmbWeightTrkPt","Non-HF electrons from embedded #pi^{0} and #eta + No mom with weight + No mom;p_{T} (GeV/c);counts",250,0,50);
  fNonHFeEmbWeightTrkPt->Sumw2();
  fOutputList->Add(fNonHFeEmbWeightTrkPt);
        
  fPi0eEmbWeightTrkPt = new TH1F("fPi0eEmbWeightTrkPt","Non-HF electrons from embedded #pi^{0} + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
  fPi0eEmbWeightTrkPt->Sumw2();
  fOutputList->Add(fPi0eEmbWeightTrkPt);
        
  fEtaeEmbWeightTrkPt = new TH1F("fEtaeEmbWeightTrkPt","Non-HF electrons from embedded #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
  fEtaeEmbWeightTrkPt->Sumw2();
  fOutputList->Add(fEtaeEmbWeightTrkPt);

  fRecoNonHFeTrkPt = new TH1F("fRecoNonHFeTrkPt"," Reco Non-HF electrons from all generators;p_{T} (GeV/c);counts",250,0,50);
  fRecoNonHFeTrkPt->Sumw2();
  fOutputList->Add(fRecoNonHFeTrkPt);

  fRecoNonHFeEmbTrkPt = new TH1F("fRecoNonHFeEmbTrkPt","Reco Non-HF electrons from embedded #pi^{0} and #eta + No mom;p_{T} (GeV/c);counts",250,0,50);
  fRecoNonHFeEmbTrkPt->Sumw2();
  fOutputList->Add(fRecoNonHFeEmbTrkPt);

  fRecoNonHFeEmbWeightTrkPt = new TH1F("fRecoNonHFeEmbWeightTrkPt","Reco Non-HF electrons from embedded #pi^{0} and #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
  fRecoNonHFeEmbWeightTrkPt->Sumw2();
  fOutputList->Add(fRecoNonHFeEmbWeightTrkPt);
        
  fRecoPi0eEmbWeightTrkPt = new TH1F("fRecoPi0eEmbWeightTrkPt","Reco Non-HF electrons from embedded #pi^{0}  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
  fRecoPi0eEmbWeightTrkPt->Sumw2();
  fOutputList->Add(fRecoPi0eEmbWeightTrkPt);
        
  fRecoEtaeEmbWeightTrkPt = new TH1F("fRecoEtaeEmbWeightTrkPt","Reco Non-HF electrons from embedded #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
  fRecoEtaeEmbWeightTrkPt->Sumw2();
  fOutputList->Add(fRecoEtaeEmbWeightTrkPt);

  fNonHFePairInvmassLS = new TH1F("fNonHFePairInvmassLS", "Inv mass of LS (e,e) if both e- are Non-HFE; mass(GeV/c^2); counts;",  50,0,0.5);
  fNonHFePairInvmassLS->Sumw2();
  fOutputList->Add(fNonHFePairInvmassLS);
        
  fNonHFePairInvmassULS = new TH1F("fNonHFePairInvmassULS", "Inv mass of ULS (e,e) if both e- are Non-HFE; mass(GeV/c^2); counts;",  50,0,0.5);
  fNonHFePairInvmassULS->Sumw2();
  fOutputList->Add(fNonHFePairInvmassULS);
        
  fNonHFeEmbInvmassLS = new TH1F("fNonHFeEmbInvmassLS", "Inv mass of LS (e,e) for Non-HFE from embedded #pi^{0} and #eta; mass(GeV/c^2); counts;",  50,0,0.5);
  fNonHFeEmbInvmassLS->Sumw2();
  fOutputList->Add(fNonHFeEmbInvmassLS);
        
  fNonHFeEmbInvmassULS = new TH1F("fNonHFeEmbInvmassULS", "Inv mass of ULS (e,e) for Non-HFE from embedded #pi^{0} and #eta; mass(GeV/c^2); counts;",  50,0,0.5);
  fNonHFeEmbInvmassULS->Sumw2();
  fOutputList->Add(fNonHFeEmbInvmassULS);
        
  fNonHFeEmbWeightInvmassLS = new TH1F("fNonHFeEmbWeightInvmassLS", "Inv mass of LS (e,e) for Non-HFE from embedded #pi^{0} and #eta with weight; mass(GeV/c^2); counts;",  50,0,0.5);
  fNonHFeEmbWeightInvmassLS->Sumw2();
  fOutputList->Add(fNonHFeEmbWeightInvmassLS);
        
  fNonHFeEmbWeightInvmassULS = new TH1F("fNonHFeEmbWeightInvmassULS", "Inv mass of ULS (e,e) for Non-HFE from embedded #pi^{0} and #eta with weight; mass(GeV/c^2); counts;",  50,0,0.5);
  fNonHFeEmbWeightInvmassULS->Sumw2();
  fOutputList->Add(fNonHFeEmbWeightInvmassULS);
        
  fPi0EmbInvmassLS = new TH1F("fPi0EmbInvmassLS", "Inv mass of LS (e,e) for ele from embedded #pi^{0}; mass(GeV/c^2); counts;",  50,0,0.5);
  fPi0EmbInvmassLS->Sumw2();
  fOutputList->Add(fPi0EmbInvmassLS);
        
  fPi0EmbInvmassULS  = new TH1F("fPi0EmbInvmassULS", "Inv mass of ULS (e,e) for ele from embedded #pi^{0}; mass(GeV/c^2); counts;",  50,0,0.5);
  fPi0EmbInvmassULS->Sumw2();
  fOutputList->Add(fPi0EmbInvmassULS);
        
  fPi0EmbWeightInvmassLS = new TH1F("fPi0EmbWeightInvmassLS", "Inv mass of LS (e,e) for ele from embedded #pi^{0} with weight; mass(GeV/c^2); counts;",  50,0,0.5);
  fPi0EmbWeightInvmassLS->Sumw2();
  fOutputList->Add(fPi0EmbWeightInvmassLS);
        
  fPi0EmbWeightInvmassULS  = new TH1F("fPi0EmbWeightInvmassULS", "Inv mass of ULS (e,e) for ele from embedded #pi^{0} with weight; mass(GeV/c^2); counts;",  50,0,0.5);
  fPi0EmbWeightInvmassULS->Sumw2();
  fOutputList->Add(fPi0EmbWeightInvmassULS);
        
  fEtaEmbInvmassLS = new TH1F("fEtaEmbInvmassLS", "Inv mass of LS (e,e) for ele from embedded #eta; mass(GeV/c^2); counts;",  50,0,0.5);
  fEtaEmbInvmassLS->Sumw2();
  fOutputList->Add(fEtaEmbInvmassLS);
        
  fEtaEmbInvmassULS = new TH1F("fEtaEmbInvmassULS", "Inv mass of ULS (e,e) for ele from embedded #eta; mass(GeV/c^2); counts;",  50,0,0.5);
  fEtaEmbInvmassULS->Sumw2();
  fOutputList->Add(fEtaEmbInvmassULS);
        
  fEtaEmbWeightInvmassLS = new TH1F("fEtaEmbWeightInvmassLS", "Inv mass of LS (e,e) for ele from embedded #eta with weight; mass(GeV/c^2); counts;",  50,0,0.5);
  fEtaEmbWeightInvmassLS->Sumw2();
  fOutputList->Add(fEtaEmbWeightInvmassLS);
        
  fEtaEmbWeightInvmassULS  = new TH1F("fEtaEmbWeightInvmassULS", "Inv mass of ULS (e,e) for ele from embedded #eta with weight; mass(GeV/c^2); counts;",  50,0,0.5);
  fEtaEmbWeightInvmassULS->Sumw2();
  fOutputList->Add(fEtaEmbWeightInvmassULS);

 fRecoLSeEmbTrkPt  = new TH1F("fRecoLSeEmbTrkPt","Reco LS electrons from embedded #pi^{0} and #eta + No mom;p_{T} (GeV/c);counts",350,0,35);
  fRecoLSeEmbTrkPt->Sumw2();
  fOutputList->Add(fRecoLSeEmbTrkPt);
        
  fRecoLSeEmbWeightTrkPt = new TH1F("fRecoLSeEmbWeightTrkPt","Reco LS electrons from embedded #pi^{0} and #eta  + No mom with weight;p_{T} (GeV/c);counts",350,0,35);
  fRecoLSeEmbWeightTrkPt->Sumw2();
  fOutputList->Add(fRecoLSeEmbWeightTrkPt);
        
  fRecoPi0LSeEmbWeightTrkPt = new TH1F("fRecoPi0LSeEmbWeightTrkPt","Reco LS electrons from embedded #pi^{0}  + No mom with weight;p_{T} (GeV/c);counts",350,0,35);
  fRecoPi0LSeEmbWeightTrkPt->Sumw2();
  fOutputList->Add(fRecoPi0LSeEmbWeightTrkPt);
        
  fRecoEtaLSeEmbWeightTrkPt  = new TH1F("fRecoEtaLSeEmbWeightTrkPt","Reco LS electrons from embedded #eta  + No mom with weight;p_{T} (GeV/c);counts",350,0,35);
  fRecoEtaLSeEmbWeightTrkPt->Sumw2();
  fOutputList->Add(fRecoEtaLSeEmbWeightTrkPt);
        
  fRecoULSeEmbTrkPt = new TH1F("fRecoULSeEmbTrkPt","Reco ULS electrons from embedded #pi^{0} and #eta + No mom;p_{T} (GeV/c);counts",350,0,35);
  fRecoULSeEmbTrkPt->Sumw2();
  fOutputList->Add(fRecoULSeEmbTrkPt);
        
  fRecoULSeEmbWeightTrkPt = new TH1F("fRecoULSeEmbWeightTrkPt","Reco ULS electrons from embedded #pi^{0} and #eta  + No mom with weight;p_{T} (GeV/c);counts",350,0,35);
  fRecoULSeEmbWeightTrkPt->Sumw2();
  fOutputList->Add(fRecoULSeEmbWeightTrkPt);
        
  fRecoPi0ULSeEmbWeightTrkPt = new TH1F("fRecoPi0ULSeEmbWeightTrkPt","Reco ULS electrons from embedded #pi^{0}  + No mom with weight;p_{T} (GeV/c);counts",350,0,35);
  fRecoPi0ULSeEmbWeightTrkPt->Sumw2();
  fOutputList->Add(fRecoPi0ULSeEmbWeightTrkPt);
  
  fRecoEtaULSeEmbWeightTrkPt = new TH1F("fRecoEtaULSeEmbWeightTrkPt","Reco ULS electrons from embedded #eta  + No mom with weight;p_{T} (GeV/c);counts",350,0,35);
  fRecoEtaULSeEmbWeightTrkPt->Sumw2();
  fOutputList->Add(fRecoEtaULSeEmbWeightTrkPt);

  fRVsULSElecPt  = new TH2F("fRVsULSElecPt","#it{p}_{T} distribution of ULS electrons;#it{p}_{T} (GeV/#it{c});Radius",250,0,50,300,0,30);
  fRVsULSElecPt->Sumw2();
  fOutputList->Add(fRVsULSElecPt);

  fRVsLSElecPt  = new TH2F("fRVsLSElecPt","#it{p}_{T} distribution of LS electrons;#it{p}_{T} (GeV/#it{c});Radius",250,0,50,300,0,30);
  fRVsLSElecPt->Sumw2();
  fOutputList->Add(fRVsLSElecPt);

  fHadConvRadius = new TH2F("fHadConvRadius","Conv Radius distribution of charged hadrons;#it{p}_{T} (GeV/#it{c});Radius",250,0,50,300,0,30);
  fHadConvRadius->Sumw2();
  fOutputList->Add(fHadConvRadius);

  fIncleConvRadius = new TH2F("fIncleConvRadius","Conv Radius distribution of Incl e;#it{p}_{T} (GeV/#it{c});Radius",250,0,50,300,0,30);
  fIncleConvRadius->Sumw2();
  fOutputList->Add(fIncleConvRadius);  

  fNonHFeConvRadius = new TH2F("fNonHFeConvRadius","Conv Radius distribution of Non-HF electrons;#it{p}_{T} (GeV/#it{c});Radius",250,0,50,300,0,30);
  fNonHFeConvRadius->Sumw2();
  fOutputList->Add(fNonHFeConvRadius);

  fHFeConvRadius = new TH2F("fHFeConvRadius","Conv Radius distribution of HF electrons;#it{p}_{T} (GeV/#it{c});Radius",250,0,50,300,0,30);
  fHFeConvRadius->Sumw2();
  fOutputList->Add(fHFeConvRadius);
  //-----------------------------------------R-Tagg-----------------------------------------------------------
  fNonHFeEmbTrkRConv = new TH2F("fNonHFeEmbTrkRConv","Conv Radius distribution of ;#it{p}_{T} (GeV/#it{c});Radius",250,0,50,300,0,30);
  fNonHFeEmbTrkRConv->Sumw2();
  fOutputList->Add(fNonHFeEmbTrkRConv);

  fPi0eEmbWeightTrkRConv = new TH2F("fPi0eEmbWeightTrkRConv","Conv Radius distribution of ;#it{p}_{T} (GeV/#it{c});Radius",250,0,50,300,0,30);
  fPi0eEmbWeightTrkRConv->Sumw2();
  fOutputList->Add(fPi0eEmbWeightTrkRConv);

  fNonHFeEmbWeightTrkRConv = new TH2F("fNonHFeEmbWeightTrkRConv","Conv Radius distribution of ;#it{p}_{T} (GeV/#it{c});Radius",250,0,50,300,0,30);
  fNonHFeEmbWeightTrkRConv->Sumw2();
  fOutputList->Add(fNonHFeEmbWeightTrkRConv);

  fEtaeEmbWeightTrkRConv = new TH2F("fEtaeEmbWeightTrkRConv","Conv Radius distribution of ;#it{p}_{T} (GeV/#it{c});Radius",250,0,50,300,0,30);
  fEtaeEmbWeightTrkRConv->Sumw2();
  fOutputList->Add(fEtaeEmbWeightTrkRConv);

  fRecoNonHFeEmbRConv = new TH2F("fRecoNonHFeEmbRConv","Conv Radius distribution of ;#it{p}_{T} (GeV/#it{c});Radius",250,0,50,300,0,30);
  fRecoNonHFeEmbRConv->Sumw2();
  fOutputList->Add(fRecoNonHFeEmbRConv);

  fRecoPi0eEmbWeightTrkRConv = new TH2F("fRecoPi0eEmbWeightTrkRConv","Conv Radius distribution of ;#it{p}_{T} (GeV/#it{c});Radius",250,0,50,300,0,30);
  fRecoPi0eEmbWeightTrkRConv->Sumw2();
  fOutputList->Add(fRecoPi0eEmbWeightTrkRConv);

  fRecoNonHFeEmbWeightTrkRConv = new TH2F("fRecoNonHFeEmbWeightTrkRConv","Conv Radius distribution of ;#it{p}_{T} (GeV/#it{c});Radius",250,0,50,300,0,30);
  fRecoNonHFeEmbWeightTrkRConv->Sumw2();
  fOutputList->Add(fRecoNonHFeEmbWeightTrkRConv);

  fRecoEtaeEmbWeightTrkRConv = new TH2F("fRecoEtaeEmbWeightTrkRConv","Conv Radius distribution of ;#it{p}_{T} (GeV/#it{c});Radius",250,0,50,300,0,30);
  fRecoEtaeEmbWeightTrkRConv->Sumw2();
  fOutputList->Add(fRecoEtaeEmbWeightTrkRConv);
  }

    if(fCalculateElecRecoEffi){
        fInclElePhysPriAll = new TH1F("fInclElePhysPriAll","Physical primary inclusive electrons for reco effi, All;p_{T} (GeV/c);counts",350,0,35);
        fInclElePhysPriAll->Sumw2();
        fOutputList->Add(fInclElePhysPriAll);
        
        fHFEPhysPriAll = new TH1F("fHFEPhysPriAll","Physical primary HFE for reco effi, All;p_{T} (GeV/c);counts",350,0,35);
        fHFEPhysPriAll->Sumw2();
        fOutputList->Add(fHFEPhysPriAll);
        
        fBEPhysPriAll = new TH1F("fBEPhysPriAll","Physical primary b->e for reco effi, All;p_{T} (GeV/c);counts",350,0,35);
        fBEPhysPriAll->Sumw2();
        fOutputList->Add(fBEPhysPriAll);
        
        fDEPhysPriAll = new TH1F("fDEPhysPriAll","Physical primary c->e for reco effi, All;p_{T} (GeV/c);counts",350,0,35);
        fDEPhysPriAll->Sumw2();
        fOutputList->Add(fDEPhysPriAll);
        
        fInclElePhysPriTrkCuts = new TH1F("fInclElePhysPriTrkCuts","Physical primary inclusive electrons for reco effi, Aft Trk cuts;p_{T} (GeV/c);counts",350,0,35);
        fInclElePhysPriTrkCuts->Sumw2();
        fOutputList->Add(fInclElePhysPriTrkCuts);
        
        fHFEPhysPriTrkCuts = new TH1F("fHFEPhysPriTrkCuts","Physical primary HFE for reco effi, Aft Trk cuts;p_{T} (GeV/c);counts",350,0,35);
        fHFEPhysPriTrkCuts->Sumw2();
        fOutputList->Add(fHFEPhysPriTrkCuts);
        
        fBEPhysPriTrkCuts = new TH1F("fBEPhysPriTrkCuts","Physical primary b->e for reco effi, Aft Trk cuts;p_{T} (GeV/c);counts",350,0,35);
        fBEPhysPriTrkCuts->Sumw2();
        fOutputList->Add(fBEPhysPriTrkCuts);
        
        fDEPhysPriTrkCuts = new TH1F("fDEPhysPriTrkCuts","Physical primary c->e for reco effi, Aft Trk cuts;p_{T} (GeV/c);counts",350,0,35);
        fDEPhysPriTrkCuts->Sumw2();
        fOutputList->Add(fDEPhysPriTrkCuts);
                
        fInclElePhysPriEMCMatch = new TH1F("fInclElePhysPriEMCMatch","Physical primary inclusive electron for reco effi, Aft EMC match;p_{T} (GeV/c);counts",350,0,35);
        fInclElePhysPriEMCMatch->Sumw2();
        fOutputList->Add(fInclElePhysPriEMCMatch);
        
        fHFEPhysPriEMCMatch = new TH1F("fHFEPhysPriEMCMatch","Physical primary HFE for reco effi, Aft EMC match;p_{T} (GeV/c);counts",350,0,35);
        fHFEPhysPriEMCMatch->Sumw2();
        fOutputList->Add(fHFEPhysPriEMCMatch);
        
        fBEPhysPriEMCMatch = new TH1F("fBEPhysPriEMCMatch","Physical primary b->e for reco effi, Aft EMC match;p_{T} (GeV/c);counts",350,0,35);
        fBEPhysPriEMCMatch->Sumw2();
        fOutputList->Add(fBEPhysPriEMCMatch);
        
        fDEPhysPriEMCMatch = new TH1F("fDEPhysPriEMCMatch","Physical primary c->e for reco effi, Aft EMC match;p_{T} (GeV/c);counts",350,0,35);
        fDEPhysPriEMCMatch->Sumw2();
        fOutputList->Add(fDEPhysPriEMCMatch);

         fInclElePhysPriEovP = new TH1F("fInclElePhysPriEovP","Physical primary inclusive electron for reco effi, Aft E/p cut;p_{T} (GeV/c);counts",350,0,35);
        fInclElePhysPriEovP->Sumw2();
        fOutputList->Add(fInclElePhysPriEovP);
        
        fHFEPhysPriEovP = new TH1F("fHFEPhysPriEovP","Physical primary HFE for reco effi, Aft E/p cut;p_{T} (GeV/c);counts",350,0,35);
        fHFEPhysPriEovP->Sumw2();
        fOutputList->Add(fHFEPhysPriEovP);
        
        fBEPhysPriEovP = new TH1F("fBEPhysPriEovP","Physical primary b->e for reco effi, Aft E/p cut;p_{T} (GeV/c);counts",350,0,35);
        fBEPhysPriEovP->Sumw2();
        fOutputList->Add(fBEPhysPriEovP);
        
        fDEPhysPriEovP = new TH1F("fDEPhysPriEovP","Physical primary c->e for reco effi, Aft E/p cut;p_{T} (GeV/c);counts",350,0,35);
        fDEPhysPriEovP->Sumw2();
        fOutputList->Add(fDEPhysPriEovP);

        fInclElePhysPriTPCnsig = new TH1F("fInclElePhysPriTPCnsig","Physical primary inclusive electron for reco effi, Aft TPCnsig;p_{T} (GeV/c);counts",350,0,35);
        fInclElePhysPriTPCnsig->Sumw2();
        fOutputList->Add(fInclElePhysPriTPCnsig);
        
        fHFEPhysPriTPCnsig = new TH1F("fHFEPhysPriTPCnsig","Physical primary HFE for reco effi, Aft TPCnsig;p_{T} (GeV/c);counts",350,0,35);
        fHFEPhysPriTPCnsig->Sumw2();
        fOutputList->Add(fHFEPhysPriTPCnsig);
        
        fBEPhysPriTPCnsig = new TH1F("fBEPhysPriTPCnsig","Physical primary b->e for reco effi, Aft TPCnsig;p_{T} (GeV/c);counts",350,0,35);
        fBEPhysPriTPCnsig->Sumw2();
        fOutputList->Add(fBEPhysPriTPCnsig);
        
        fDEPhysPriTPCnsig = new TH1F("fDEPhysPriTPCnsig","Physical primary c->e for reco effi, Aft TPCnsig;p_{T} (GeV/c);counts",350,0,35);
        fDEPhysPriTPCnsig->Sumw2();
        fOutputList->Add(fDEPhysPriTPCnsig);
                
        fInclElePhysPriSS = new TH1F("fInclElePhysPriSS","Physical primary inclusive electron for reco effi, Aft SS cut;p_{T} (GeV/c);counts",350,0,35);
        fInclElePhysPriSS->Sumw2();
        fOutputList->Add(fInclElePhysPriSS);
        
        fHFEPhysPriSS = new TH1F("fHFEPhysPriSS","Physical primary HFE for reco effi, Aft SS cut;p_{T} (GeV/c);counts",350,0,35);
        fHFEPhysPriSS->Sumw2();
        fOutputList->Add(fHFEPhysPriSS);
        
        fBEPhysPriSS = new TH1F("fBEPhysPriSS","Physical primary b->e for reco effi, Aft SS cut;p_{T} (GeV/c);counts",350,0,35);
        fBEPhysPriSS->Sumw2();
        fOutputList->Add(fBEPhysPriSS);
        
        fDEPhysPriSS = new TH1F("fDEPhysPriSS","Physical primary c->e for reco effi, Aft SS cut;p_{T} (GeV/c);counts",350,0,35);
        fDEPhysPriSS->Sumw2();
        fOutputList->Add(fDEPhysPriSS);
        
    }

//---------------------------DCA Templates --------------------------------------------------------------------------------------------------

  if(fCalculateMCTemplWeightCalc)
  {

    fBHadpT = new TH1F("fBHadpT","B hadron pT;p_{T} (GeV/c);counts",100,0,100);
    fBHadpT->Sumw2();
    fOutputList->Add(fBHadpT);
        
    fBMesonpT = new TH1F("fBMesonpT","B meson pT;p_{T} (GeV/c);counts",100,0,100);
    fBMesonpT->Sumw2();
    fOutputList->Add(fBMesonpT);
        
    fBDHadpT = new TH1F("fBDHadpT","D (<- B) hadron pT;p_{T} (GeV/c);counts",100,0,100);
    fBDHadpT->Sumw2();
    fOutputList->Add(fBDHadpT);

    Int_t nBins = 11;  const Float_t D0PtBinsPP[11] = { 1., 2., 3., 4., 5., 6., 7., 8., 10., 12., 24. };

    fDHadpT = new TH1F("fDHadpT","Prompt D hadron pT;p_{T} (GeV/c);counts", 100,0,100);
    fDHadpT->Sumw2();
    fOutputList->Add(fDHadpT);
        
    fDMesonpT = new TH1F("fDMesonpT","Prompt D meson pT;p_{T} (GeV/c);counts", 100,0,100);
    fDMesonpT->Sumw2();
    fOutputList->Add(fDMesonpT);

    fD0pT = new TH1F("fD0pT","Prompt D0 meson pT;p_{T} (GeV/c);counts", 100,0,100);
    fD0pT->Sumw2();
    fOutputList->Add(fD0pT);
        
    fDPluspT = new TH1F("fDPluspT","Prompt D+ meson pT;p_{T} (GeV/c);counts", 100,0,100);
    fDPluspT->Sumw2();
    fOutputList->Add(fDPluspT);
        
    fDspT = new TH1F("fDspT","Prompt D+s meson pT;p_{T} (GeV/c);counts", 100,0,100);
    fDspT->Sumw2();
    fOutputList->Add(fDspT);
        
    fLambdaCpT = new TH1F("fLambdaCpT","Prompt Lammda_c pT;p_{T} (GeV/c);counts", 100,0,100);
    fLambdaCpT->Sumw2();
    fOutputList->Add(fLambdaCpT);
  }

  if(fFillMCTemplates)
  {
        fDElecDCA = new TH2F("fDElecDCA","D meson -> electron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 70,0,35., fnBinsDCAHisto,-0.4,0.4);
        fOutputList->Add(fDElecDCA);
        
        fBElecDCA = new TH2F("fBElecDCA","B meson -> electron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 70,0,35., fnBinsDCAHisto,-0.4,0.4);
        fOutputList->Add(fBElecDCA);
        
        fBHadElecDCA = new TH2F("fBHadElecDCA","B hadron -> electron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 70,0,35., fnBinsDCAHisto,-0.4,0.4);
        fOutputList->Add(fBHadElecDCA);
        
        fBMesonElecDCA = new TH2F("fBMesonElecDCA","B meson -> electron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 70,0,35., fnBinsDCAHisto,-0.4,0.4);
        fOutputList->Add(fBMesonElecDCA);
        
        fBBaryonElecDCA = new TH2F("fBBaryonElecDCA","B baryon -> electron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 70,0,35., fnBinsDCAHisto,-0.4,0.4);
        fOutputList->Add(fBBaryonElecDCA);
        
        fDHadElecDCA = new TH2F("fDHadElecDCA","D hadron -> electron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 70,0,35., fnBinsDCAHisto,-0.4,0.4);
        fOutputList->Add(fDHadElecDCA);
        
        fDMesonElecDCA = new TH2F("fDMesonElecDCA","D meson -> electron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 70,0,35., fnBinsDCAHisto,-0.4,0.4);
        fOutputList->Add(fDMesonElecDCA);
        
        fDBaryonElecDCA  = new TH2F("fDBaryonElecDCA","D baryon -> electron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 70,0,35., fnBinsDCAHisto,-0.4,0.4);
        fOutputList->Add(fDBaryonElecDCA);
        
        fLambdaCElecDCA = new TH2F("fLambdaCElecDCA","Lambda_c -> electron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 70,0,35., fnBinsDCAHisto,-0.4,0.4);
        fOutputList->Add(fLambdaCElecDCA);
        
        fD0ElecDCA = new TH2F("fD0ElecDCA","D0 -> electron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 70,0,35., fnBinsDCAHisto,-0.4,0.4);
        fOutputList->Add(fD0ElecDCA);
        
        Int_t binTemp[3] = {70,fnBinsDCAHisto,19}; //pT, DCA, Mom PID, Mom Gen, mompT
        Double_t xminTemp[3] = {0.,-0.2,0.5};
        Double_t xmaxTemp[3] = {35.,0.2,19.5};

        fSprsTemplatesNoWeight = new THnSparseD("fSprsTemplatesNoWeight","Sparse for DCA Templates, No weight applied;p_{T};DCA;MomPID",3,binTemp,xminTemp,xmaxTemp);
        fSprsTemplatesNoWeight->Sumw2();
        fOutputList->Add(fSprsTemplatesNoWeight);
        
        fSprsTemplatesWeight = new THnSparseD("fSprsTemplatesWeight","Sparse for DCA Templates,With weight applied;p_{T};DCA;MomPID",3,binTemp,xminTemp,xmaxTemp);
        fSprsTemplatesWeight->Sumw2();
        fOutputList->Add(fSprsTemplatesWeight);
        
        fSprsTemplatesWeightVar1 = new THnSparseD("fSprsTemplatesWeightVar1","Sparse for DCA Templates,With weight variation 1 applied;p_{T};DCA;MomPID",3,binTemp,xminTemp,xmaxTemp);
        fSprsTemplatesWeightVar1->Sumw2();
        fOutputList->Add(fSprsTemplatesWeightVar1);
        
        fSprsTemplatesWeightVar2= new THnSparseD("fSprsTemplatesWeightVar2","Sparse for DCA Templates,With weight variation 2 applied;p_{T};DCA;MomPID",3,binTemp,xminTemp,xmaxTemp);
        fSprsTemplatesWeightVar2->Sumw2();
        fOutputList->Add(fSprsTemplatesWeightVar2);


        fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");
        
  }


}

  PostData(1, fOutputList);
}


//________________________________________________________________________
Int_t eventNo=0; 

void AliAnalysispp13TeVEMCalHFEMCReco::UserExec(Option_t *) 
{


  Int_t Nch=0,count=0,c=0;    

  Int_t pdg = -99999;
  Int_t pdg_mother = -99999;
  Int_t pidM = -1;
  
  Double_t fTPCnSigma = -999; 

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fAOD) { return; }
  
  fHistEvent->Fill(1);      // total # of evts
    
  Bool_t isSelected = (((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected() & ftrigger);
  //cout<<" isSelected  "<<  isSelected  <<" ftrigger  "<<  ftrigger  <<endl;
  if(!isSelected) return;

  fHistEvent->Fill(3);  // # of evts after Trigger
 
  if(fUseTender){
  fTracks_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTenderTrackName));
  fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTenderClusterName)); //emcal correction
  }

  //////////////////
  // Tigger Check //
  /////////////////
  //-------------------selecting trigger for calorimeter( EMCAL + DCAL )
  TString firedTrigger;
  TString TriggerEG1("EG1");
  TString TriggerEG2("EG2");
  TString TriggerDG1("DG1");
  TString TriggerDG2("DG2");
    
  if(fAOD) firedTrigger = fAOD->GetFiredTriggerClasses();
  
  if(fEMCEG2 && fDCalDG2) if(!firedTrigger.Contains(TriggerEG2) && !firedTrigger.Contains(TriggerDG2)) return;
  if(fEMCEG1 && fDCalDG1) if(!firedTrigger.Contains(TriggerEG1) && !firedTrigger.Contains(TriggerDG1)) return;

  if(fDCalDG2 && !fEMCEG2) { if(!firedTrigger.Contains(TriggerDG2))return; }
  if(fEMCEG1  && !fDCalDG1){ if(!firedTrigger.Contains(TriggerEG1))return; }
  if(fEMCEG2  && !fDCalDG2){ if(!firedTrigger.Contains(TriggerEG2))return; }
  if(fDCalDG1 && !fEMCEG1) { if(!firedTrigger.Contains(TriggerDG1))return; }

  fHistEvent->Fill(5); // Number of events after passing EMCal Trigger

  //cout<<" *********************After Trigger***************  "<<endl;
 
  const AliVVertex *vertex=fAOD->GetPrimaryVertex();
  if(!vertex){return;}

  fHistVz->Fill(vertex->GetZ());
 
  if(vertex->GetNContributors() <2 ) return;   
  fHistEvent->Fill(7);

  Bool_t isPileupfromSPDmulbins=fAOD->IsPileupFromSPDInMultBins(); //This function checks if there was a pile up reconstructed with SPD
  if(isPileupfromSPDmulbins) return;
  fHistEvent->Fill(9);

  //cout<<" *********************After PileUp***************  "<<endl;

 Int_t minContributors=5;    //minimum contributors to the pilup vertices, multi-vertex
 Double_t minChi2=5.; 
 Double_t minWeiZDiff=15;   //minimum of the sqrt of weighted distance between the primary and the pilup vertex, multi-vertex
 Bool_t checkPlpFromDifferentBC=kFALSE;
  
 AliAnalysisUtils utils;
 utils.SetMinPlpContribMV(minContributors); //Multi Vertex pileup selection
 utils.SetMaxPlpChi2MV(minChi2);   //max value of Chi2perNDF of the pileup vertex, multi-vertex
 utils.SetMinWDistMV(minWeiZDiff);
 utils.SetCheckPlpFromDifferentBCMV(checkPlpFromDifferentBC); //SPD Pileup slection
 Bool_t isPileupFromMV = utils.IsPileUpMV(fAOD);      //check for multi-vertexer pile-up
  
 if(isPileupFromMV) return;
 fHistEvent->Fill(11);

  //cout<<" *********************After PileUpMV***************  "<<endl;

 if(TMath::Abs( vertex->GetZ() ) > 10) return; 
 fHistEvent->Fill(13);

  fHistVzwc->Fill(vertex->GetZ());

 //PID response
  fPidResponse = fInputHandler->GetPIDResponse();  

//Check PID response
    if(!fPidResponse)
    {
        AliDebug(1, "Using default PID Response");
        fPidResponse = AliHFEtools::GetDefaultPID(kFALSE, fInputEvent->IsA() == AliAODEvent::Class());        
    }

    fPID->SetPIDResponse(fPidResponse);

  //cout<<" *********************After PIDRespons***************  "<<endl;

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
  if(!fExtraCuts)
  fExtraCuts = new AliHFEextraCuts("hfeExtraCuts","HFE Extra Cuts");
  fExtraCuts->SetRecEventInfo(fAOD);
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
    ////////////////
    // Mag. field //
    ////////////////
    Int_t fMagSign = 1;
    if(fAOD->GetMagneticField()<0) fMagSign = -1;

  //cout<<"  CXX fIsMC  ===    "<< fIsMC<<endl;
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=
  ///////////////////////////////////
  //Initialization for MC analysis///
  ///////////////////////////////////
if(fIsMC)
{
 

   //cout<<" *********************Entering fIsMC***************  "<<endl;

    Double_t qaweights[5];
    Double_t pi0etaweights[3];

    fMCArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    if(!fMCArray){ AliError("Array of MC particles not found");
    return;
    }

    fMCHeader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    if(fMCHeader){
        ////////////////////////////////
        //Get number of Gen particles //
        ////////////////////////////////
        GetNMCPartProduced();
        
   //cout<<" *********************After GetNMCPartProduced***************  "<<endl;


        /////////////////////////////////
        //Calculate Pi0 and Eta weight //
        ///////////////////////////////// 
        if(fCalculateWeight) GetPi0EtaWeight(fSprsPi0EtaWeightCal);

   //cout<<" *********************After GetPi0EtaWeight***************  "<<endl;

        /////////////////////////
        //Electrons in MC stack//
        /////////////////////////
        if(fCalculateElecRecoEffi) GetElectronFromStack();

   //cout<<" *********************After GetElectronFromStack***************  "<<endl;

        /////////////////////////////////
        //Histos for MC template Weight//
        /////////////////////////////////    
        if(fCalculateMCTemplWeightCalc) GetMCTemplateWeight();
    
       //cout<<" *********************After GetMCTemplateWeight***************  "<<endl;

    }
    if (!fMCHeader) {
    AliError("Could not find MC Header in AOD");
    return;
    }

    //if(fMCArray->GetEntries() < 1) return; 

    for(Int_t iMC = 0; iMC < fMCArray->GetEntries(); iMC++)
    {
      fMCparticle = (AliAODMCParticle*) fMCArray->At(iMC);
      Int_t pdg = TMath::Abs(fMCparticle->GetPdgCode());

     //For the reconstruction efficiency of HFE:-----------                          
     if(TMath::Abs(fMCparticle->Eta()) <= fEtarange  )
      {        ///Pseudo-rapidity cut                 
       if(TMath::Abs(pdg) == 11)
       {
        Bool_t MotherFound = FindMother(iMC);
        if(MotherFound)
         {
          if(fIsHFE1){ fPthfeGenerated->Fill(fMCparticle->Pt()); }
          }
        }   // pdg condition
      } 
   
      if(fMCparticle->Y()<-0.8 || fMCparticle->Y()>0.8) continue;     

      if(fMCparticle->IsPrimary()) ///Does not include particles from weak decays or created in an interaction with the material 
      { 
        
          if (TMath::Abs(fMCparticle->GetPdgCode())==111) pi0MC->Fill(fMCparticle->Pt());  // pdg=111 pi0
          if (TMath::Abs(fMCparticle->GetPdgCode())==221) etaMC->Fill(fMCparticle->Pt());  // pdg=221 eta 
          if (TMath::Abs(fMCparticle->GetPdgCode())==22) gammaMC->Fill(fMCparticle->Pt());   // pdg=22  gamma
      
          Int_t type = GetPi0EtaType(fMCparticle);
            
          ///Using thnSparse--------------------------------------
          //Pt
        
          pi0etaweights[0] = fMCparticle->Pt();
            
          // What pdg
          pi0etaweights[1]=-1.;
          if (TMath::Abs(fMCparticle->GetPdgCode())==111) pi0etaweights[1]=0.2;  // pdg=111 pi0
          if (TMath::Abs(fMCparticle->GetPdgCode())==221) pi0etaweights[1]=1.2;  // pdg=221 eta
          if (TMath::Abs(fMCparticle->GetPdgCode())==22) pi0etaweights[1]=2.2;   // pdg=22  gamma
        
          // What type
          //Int_t type = GetPi0EtaType(fMCparticle);
        
          pi0etaweights[2]=type;

          if(pi0etaweights[1]>0.) fPi0EtaSpectraSp->Fill(pi0etaweights);
              
        ///-----------------------------------------------------    
          }   //IsPrimary()loop 

    }    //For loop


}   //IsMC 

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++=

//////////////////////////////////////////////////////////////
//////////  EMcal and Dcal  /////////////////////////////////
////////////////////////////////////////////////////////////
//--------------------------cluster properties-----------------------------------------------------------------------------------------

Int_t Nclust = -999; Double_t clustE = -999; 

 if(!fUseTender) Nclust = fAOD->GetNumberOfCaloClusters();
 if(fUseTender) Nclust = fCaloClusters_tender->GetEntries();

 Bool_t fClsTypeEMC = kFALSE, fClsTypeDCAL = kFALSE;

for(Int_t icl=0; icl<Nclust; icl++)
{ 
      AliAODCaloCluster *clust = 0x0;
     if(!fUseTender) clust = (AliAODCaloCluster*)fAOD->GetCaloCluster(icl) ;
     if(fUseTender) clust = dynamic_cast<AliAODCaloCluster*>(fCaloClusters_tender->At(icl));
     if(!clust)  continue;   //printf("ERROR: Could not receive cluster matched calibrated from track %d\n", icl);

      fClsTypeEMC = kFALSE; fClsTypeDCAL = kFALSE;
     
     if(clust->IsEMCAL())
     {
        //Removing exotic clusters using IsExotic function in data and using M02 min cut
        /*if(!fIsMC)
        if(clust->GetIsExotic()) continue;
        Double_t m02 = clust->GetM02();
        if(m02 < 0.02) continue;*/
            
       //cout<<" *********************After EMCal-DCal Selection***************  "<<endl;

         AliAODCaloCells &cells = *(fAOD->GetEMCALCells());
         Double_t clustE = clust->E(); //clust->E();
         //if(clustE < 0.3) continue;  // switched off to match Shreyasi

         /////////////////////////////////
         //Select EMCAL or DCAL clusters//
         /////////////////////////////////
 
         Float_t  emcx[3]; // cluster pos
         clust->GetPosition(emcx);
         TVector3 clustpos(emcx[0],emcx[1],emcx[2]);
         Double_t emcphi = clustpos.Phi();
         Double_t emceta = clustpos.Eta();
         if(emcphi < 0) emcphi = emcphi+(2*TMath::Pi()); //TLorentz vector is defined between -pi to pi, so negative phi has to be flipped

         if(emcphi > 1.39 && emcphi < 3.265) fClsTypeEMC = kTRUE; //EMCAL : 80 < phi < 187
         if(emcphi > 4.53 && emcphi < 5.708) fClsTypeDCAL = kTRUE; //DCAL  : 260 < phi < 327
         
         if(fFlagClsTypeEMC && !fFlagClsTypeDCAL)
         if(!fClsTypeEMC) continue; //selecting only EMCAL clusters
         
         if(fFlagClsTypeDCAL && !fFlagClsTypeEMC)
         if(!fClsTypeDCAL) continue; //selecting only DCAL clusters
        
         {     
          fHistClustE->Fill(clustE);
          fEMCClsEtaPhi->Fill(emceta,emcphi);
          fHistoNCells->Fill(clustE,clust->GetNCells());
          float_t tof = clust->GetTOF()*1e+9; // ns
          fHistoTimeEMC->Fill(clustE,tof);
        }
    }
 
}

//-------------------------------------------------------------------------------------------------------------------------------------

    ///////////////
    //Track loop///
    ///////////////
 
  Int_t ntracks =-999;
  if(!fUseTender)ntracks = fAOD->GetNumberOfTracks();
  if(fUseTender) ntracks = fTracks_tender->GetEntries(); 
  //cout<<"-----TRACK-LOOP-------------STARTED------------- "<<endl;    
  fNEle = 0;

    for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) 
    {
        AliAODTrack* track = 0x0;
           
        if(!fUseTender) track = (AliAODTrack*)fAOD->GetTrack(iTracks);
        if(fUseTender) track = dynamic_cast<AliAODTrack*>(fTracks_tender->At(iTracks));
        if (!track) continue;


        ///////////////////////
        // Get MC information//
        ///////////////////////
        Int_t pdg = -999;
        Int_t pidM = -1;
        Double_t pid_ele = 0.0;
        Bool_t IsMCEle = kFALSE, IsMCPPEle = kFALSE, IsMCHFEle = kFALSE, IsMCDEle = kFALSE, IsMCBEle = kFALSE;

        if(fIsMC && fMCHeader){
            GetTrackHFStatus(track, IsMCEle, IsMCPPEle, IsMCHFEle, IsMCBEle, IsMCDEle);
        }

       /*Double_t p  = track->GetTPCmomentum(); 
       Double_t pt = track->Pt();
       Double_t eta = track->Eta();
       Double_t phi = track->Phi(); */
  
    Int_t tracktypeTrig=0;
    tracktypeTrig=ClassifyTrack(track,vertex);  //track cuts applied

    if(tracktypeTrig!=1) continue;    //==========TRACK cuts Applied =====
    Nch++;
    
    //cout<<" *********************After Trk Selection***************  "<<endl;
    
          fHistPt->Fill(track->Pt());
          fHistEta->Fill(track->Eta());
          fHistPhi->Fill(track->Phi());          

         Double_t d0z0[2]={-999,-999}, cov[3];
         if(track->PropagateToDCA(vertex, fAOD->GetMagneticField(), 20., d0z0, cov)) 
          fHistdca->Fill( d0z0[0],d0z0[1]); 
                  
          fTPCnSigma = fPidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);

          fHistBethe->Fill(track->GetTPCmomentum(),track->GetTPCsignal());  
          fnSigmaVsP_TPC->Fill(track->GetTPCmomentum(),fTPCnSigma);

    // hf reconstruction efficiency block-----------

    if(fIsMC)
    {
      Bool_t IsHFEMC = IsHFelectronsMC(track);
      if(IsHFEMC){ fPthfe_rec->Fill(track->Pt()); }
    }

    //hf reconstruction efficiency block After Trk Selection-----------
    if(fIsMC && fIsAOD)
    {
      Bool_t IsHFEMC = IsHFelectronsMC(track);
      if(IsHFEMC)
      { 
        fPthfe_rec_TrkSel->Fill(track->Pt()); 

        Double_t Rconv; TrackConvRadius(track, Rconv); 
        fHFeConvRadius->Fill(track->Pt(), Rconv);              // Conversion radius for HFe befote track match
        
      }
    }

    /////////////////////////////
    //Reconstruction efficiency//
    /////////////////////////////
    if(fCalculateElecRecoEffi){
        Double_t TrkPt = track->Pt();
        if(IsMCPPEle) fInclElePhysPriTrkCuts->Fill(TrkPt);
        if(IsMCHFEle) fHFEPhysPriTrkCuts->Fill(TrkPt);
        if(IsMCBEle) fBEPhysPriTrkCuts->Fill(TrkPt);
        if(IsMCDEle) fDEPhysPriTrkCuts->Fill(TrkPt);
    }
  

    if(track->PropagateToDCA(vertex, fAOD->GetMagneticField(), 20., d0z0, cov))  
    fTrkDCA = -999.0;
    fTrkDCA = d0z0[0] * track->Charge() * fMagSign;


    double hfeImpactParam = -999., hfeImpactParamResol = -999.;
    fExtraCuts->GetHFEImpactParameters((AliVTrack *)track, hfeImpactParam, hfeImpactParamResol);
    double IP = hfeImpactParam*fMagSign*track->Charge();
    //cout<< " &&&&&&&&&&&&&&&&&&& I.P.  "<<IP<<endl; 
    //cout<<" After ReCal Imp ----  "<<d0z0[0]<<endl; 


    Bool_t fFillTem = kFALSE;
    if(fFillMCTemplates)   //
    {
       fFillTem = GetMCDCATemplates(track, IP);
              
    }  
//--------------------------cluster matched to tpc properties-------------------------------------------
        ///////////////////////////
        //Track matching to EMCAL//
        //////////////////////////
  //Double_t clustE =-999;
  //if(!track->IsEMCAL()) continue;   //----Not used by Shreyasi-----   //Matches For both EMCal as well DCal
  
  Int_t EMCalIndex = -1;
  EMCalIndex = track->GetEMCALcluster();
  if(EMCalIndex < 0) continue;  
  //if( pt < 0.5) continue;      
  fHistPtMatch->Fill(track->Pt());

  AliAODCaloCluster *clustMatch=0x0;

  if(!fUseTender) if(EMCalIndex >= 0) clustMatch = (AliAODCaloCluster*)fAOD->GetCaloCluster(EMCalIndex) ;
  if(fUseTender) if(EMCalIndex >= 0)clustMatch = dynamic_cast<AliAODCaloCluster*>(fCaloClusters_tender->At(EMCalIndex));
        
  Short_t NcellsInCluster = clustMatch->GetNCells();

        Double_t emcphi = -999, emceta=-999;
        fClsTypeEMC = kFALSE; fClsTypeDCAL = kFALSE;
        if(clustMatch && clustMatch->IsEMCAL())
        {
            Double_t fPhiDiff = -999, fEtaDiff = -999;
            GetTrkClsEtaPhiDiff(track, clustMatch, fPhiDiff, fEtaDiff);
            
            if(TMath::Abs(fPhiDiff) > 0.01 || TMath::Abs(fEtaDiff)> 0.01) continue; //-------track matching condition------------------
         
            /////////////////////////////////
            //Select EMCAL or DCAL clusters//
            /////////////////////////////////
            Float_t  emcx[3]; // cluster pos
            clustMatch->GetPosition(emcx);
            TVector3 clustpos(emcx[0],emcx[1],emcx[2]);
            emcphi = clustpos.Phi();
            emceta = clustpos.Eta();
            if(emcphi < 0) emcphi = emcphi+(2*TMath::Pi());//TLorentz vector is defined between -pi to pi,so negative phi has to be flipped.
            if(emcphi > 1.39 && emcphi < 3.265) fClsTypeEMC = kTRUE; //EMCAL : 80 < phi < 187
            if(emcphi > 4.53 && emcphi < 5.708) fClsTypeDCAL = kTRUE;//DCAL  : 260 < phi < 327
            
            if(fFlagClsTypeEMC && !fFlagClsTypeDCAL)
            if(!fClsTypeEMC) continue; //selecting only EMCAL clusters
            
            if(fFlagClsTypeDCAL && !fFlagClsTypeEMC)
            if(!fClsTypeDCAL) continue; //selecting only DCAL clusters

            Double_t clustTime = clustMatch->GetTOF()*1e+9; // ns;
           /////////// Switching off to Match with Shreyasi in Data Part
           /* if(!fIsMC){
              if(clustMatch->GetIsExotic()) continue; //remove exotic clusters
              //if(TMath::Abs(clustTime) > 50) continue; //50ns time cut to remove pileup not sure if I need this cut
            }*/

            /////////////////////////////
            //Reconstruction efficiency//
            /////////////////////////////
            if(fCalculateElecRecoEffi){
                Double_t TrkPt = track->Pt();
                if(IsMCPPEle) fInclElePhysPriEMCMatch->Fill(TrkPt);
                if(IsMCHFEle) fHFEPhysPriEMCMatch->Fill(TrkPt);
                if(IsMCBEle) fBEPhysPriEMCMatch->Fill(TrkPt);
                if(IsMCDEle) fDEPhysPriEMCMatch->Fill(TrkPt);

            }

            fEMCTrkMatch->Fill(fPhiDiff,fEtaDiff);
            fEMCTrkMatch_Phi->Fill(track->Pt(),fPhiDiff);
            fEMCTrkMatch_Eta->Fill(track->Pt(),fEtaDiff);        
            fEMCClsEtaPhiTrkMatch->Fill(emceta,emcphi);
   
           Double_t Etrkmatch = -999.0, Eoptrk = -999.0 , M02trkmatch = -999.0, M20trkmatch = -999.0;
           Etrkmatch = clustMatch->E(); //clustMatch->E();
           if(track->P() > 0 ) Eoptrk = Etrkmatch/track->P(); 
           M02trkmatch = clustMatch->GetM02();
           M20trkmatch = clustMatch->GetM20();
           
           if(track->Pt() < 2) continue;   // Added from Shreyasi to Match RecoEfficiency
           fvalueElectron[0] = track->Pt(); //matched tracks pt
           fvalueElectron[1] = fTPCnSigma; // tpc n sigma
           fvalueElectron[2] = Eoptrk; //E/P
           fvalueElectron[3] = M02trkmatch; // shower shape cut
           fvalueElectron[4] = M20trkmatch;
           fvalueElectron[5] = Etrkmatch; //cluster energy after matching
           fSparseElectron->Fill(fvalueElectron); //Electron information spa
     
            /////////////////////////////
            //Reconstruction efficiency//
            /////////////////////////////
            if(fCalculateElecRecoEffi){
            GetEIDRecoEffi(track, clustMatch, IsMCPPEle, IsMCHFEle, IsMCBEle, IsMCDEle, fTPCnSigma);
            }
            
            //////////////////
            //Apply EID cuts//
            //////////////////
            
            Bool_t fHadTrack = kFALSE, fElectTrack = kFALSE;
            fElectTrack = PassEIDCuts(track, clustMatch, fHadTrack);
            
            if(fHadTrack){
                fHadPt_AftEID->Fill(track->Pt());
                fHadDCA->Fill(track->Pt(),fTrkDCA);               
              
                if(fIsMC){
                  Double_t Rconv = -999;  
                  TrackConvRadius(track, Rconv);
                  fHadConvRadius->Fill(track->Pt(),Rconv); 
                }
         }
            
            if(!fElectTrack) continue;
            
            fInclsElecPt->Fill(track->Pt());
            fInclElecDCA->Fill(track->Pt(),fTrkDCA);

            fNEle++;      

            //////////////////////////////////
            //Non-HFE efficiency calculation//
            //////////////////////////////////
            Bool_t EffiDenom = kFALSE;
            Bool_t EffiNumTag = kFALSE;
            if(fMCHeader && fCalculateNonHFEEffi){
            EffiDenom = GetNonHFEEffiDenom(track);
            }

            ////////////////////
            //NonHFE selection//
            ////////////////////
            Bool_t fFlagNonHFE=kFALSE; Int_t pidM = -1;
            SelectPhotonicElectron(iTracks,track,fFlagNonHFE,pidM);

            //////////////////////////////////
            //Non-HFE efficiency calculation//
            //////////////////////////////////
            if(fMCHeader && fCalculateNonHFEEffi){
                if(fFlagNonHFE){
                    EffiNumTag = GetNonHFEEffiRecoTag(track);
                }
            } 
              
                    

        }// EMCal Trk Match
  } 
//===================================track loop ends here =========================================
  
    fHistEvent->Fill(15);
    fHistMul->Fill(Nch);
    fNElecInEvt->Fill(fNEle);
    PostData(1, fOutputList);

}  // event loop     

//=================================================================================================================================
Int_t AliAnalysispp13TeVEMCalHFEMCReco::ClassifyTrack(AliAODTrack* track,const AliVVertex* vertex)
{  
  
  Double_t pt = track->Pt();
  Double_t eta = track->Eta();
  Double_t phi = track->Phi();
  Double_t d0z0[2]={-999,-999}, cov[3];

    //====kink daughters
  Int_t numberofvertices = 100;
  numberofvertices = fAOD->GetNumberOfVertices();
  Double_t listofmotherkink[numberofvertices];
  Int_t numberofmotherkink = 0;
 
  for(Int_t ivertex=0; ivertex < numberofvertices; ivertex++) 
  {
    AliAODVertex *aodvertex = fAOD->GetVertex(ivertex);
    if(!aodvertex) continue;
    if(aodvertex->GetType()== AliAODVertex::kKink) 
    {
      AliAODTrack *mother = (AliAODTrack *) aodvertex->GetParent();
      if(!mother) continue;
      Int_t idmother = mother->GetID();
      listofmotherkink[numberofmotherkink] = idmother;
      numberofmotherkink++;
    }
  }
 
  //reject kink
  Bool_t kinkmotherpass = kTRUE;
  for(Int_t kinkmother = 0; kinkmother < numberofmotherkink; kinkmother++) 
  {
    if(track->GetID() == listofmotherkink[kinkmother]) 
    {
      kinkmotherpass = kFALSE;
      continue;
    }
  }
  if(!kinkmotherpass) return kFALSE;


  if(!track->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return 0; //fitler bit 

  if(track->Pt()  < 0.5 ) return 0;                                  //Pt cut

  if (TMath::Abs(track->Eta()) > fEtarange ) return 0;                     //Eta cut

  if(track->GetITSNcls() < fITSNclus) return 0; // ITS N clusters 

  Double_t TPCNClsF =  track->GetTPCNclsF();
  Double_t TPCNCrossedRows = track->GetTPCNCrossedRows();
  Double_t RatioCrossedRowsOverFindableClusters = -999;
  Double_t nclusN = track->GetTPCsignalN();   

  if( TPCNCrossedRows < fTPCNCrRows) return kFALSE; //TPC N Crossed rows
  if(TPCNClsF > 0){ RatioCrossedRowsOverFindableClusters = TPCNCrossedRows/TPCNClsF; }

  if(RatioCrossedRowsOverFindableClusters < fRatioCrossedRowOverFindable) return kFALSE;

  if(nclusN < fTPCNclusPID) return 0 ;

  if((!(track->GetStatus()&AliAODTrack::kITSrefit)|| (!(track->GetStatus()&AliAODTrack::kTPCrefit)))) return 0;// ITS and TPC refit 

  if(fSPDAny){ if(!(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1))) return 0;} //Hit on first and second SPD layer : kAny
  else if(fSPDBoth){ if(!(track->HasPointOnITSLayer(0) && track->HasPointOnITSLayer(1))) return 0;} //Hit on first and second SPD layer
  else if(fSPDFirst){ if(!(track->HasPointOnITSLayer(0))) return 0;} //Hit on first and second SPD layer

  if(fRecalIP) RecalImpactParam(track, d0z0, cov);

  if(track->PropagateToDCA(vertex, fAOD->GetMagneticField(), 20., d0z0, cov)) 
  if(TMath::Abs(d0z0[0]) > fDCAxyCut || TMath::Abs(d0z0[1]) > fDCAzCut) return 0; 
     
  Double_t chi2ndf = track->Chi2perNDF();
  if(chi2ndf>4.0) return 0;
      
  return 1;
}

//=================================================================================================================================
void AliAnalysispp13TeVEMCalHFEMCReco::RecalImpactParam(const AliAODTrack * const track, Double_t dcaD[2], Double_t covD[3])
{
    //Recalculate impact parameter by recalculating primary vertex
    
    const Double_t kBeampiperadius=3.0;
    Bool_t isRecalcVertex = kFALSE;

    AliAODVertex *vtxAODSkip  = fAOD->GetPrimaryVertex();
    if(!vtxAODSkip) return;
    
    Double_t fMagField = fAOD->GetMagneticField();

    const AliAODTrack *tmptrack = dynamic_cast<const AliAODTrack *>(track);
    if(tmptrack){
        if(vtxAODSkip->GetNContributors() < 30){ // if vertex contributor is smaller than 30, recalculate the primary vertex
            
            vtxAODSkip = RemoveDaughtersFromPrimaryVtx(track);
            isRecalcVertex = kTRUE;
        }
        
        if(vtxAODSkip){
            AliAODTrack aodtrack(*tmptrack);
            AliExternalTrackParam etp;
            etp.CopyFromVTrack(&aodtrack);
            
            etp.PropagateToDCA(vtxAODSkip, fMagField, kBeampiperadius, dcaD, covD);
            
            if(isRecalcVertex) delete vtxAODSkip;
        }
    }
}

//=================================================================================================================================

AliAODVertex* AliAnalysispp13TeVEMCalHFEMCReco::RemoveDaughtersFromPrimaryVtx(const AliAODTrack * const track)
{
    // This method returns a primary vertex without the daughter tracks of the
    // candidate and it recalculates the impact parameters and errors for AOD tracks.
    
    AliAODVertex *vtxAOD = fAOD->GetPrimaryVertex();
    if(!vtxAOD) return 0;
    TString title=vtxAOD->GetTitle();
    if(!title.Contains("VertexerTracks")) return 0;

    AliVertexerTracks vertexer(fAOD->GetMagneticField());
    
    vertexer.SetITSMode();
    vertexer.SetMinClusters(3);
    vertexer.SetConstraintOff();
    
    if(title.Contains("WithConstraint")) {
        Float_t diamondcovxy[3];
        fAOD->GetDiamondCovXY(diamondcovxy);
        Double_t pos[3]={fAOD->GetDiamondX(),fAOD->GetDiamondY(),0.};
        Double_t cov[6]={diamondcovxy[0],diamondcovxy[1],diamondcovxy[2],0.,0.,10.*10.};
        AliESDVertex diamond(pos,cov,1.,1);
        vertexer.SetVtxStart(&diamond);
    }
    Int_t skipped[2]; for(Int_t i=0;i<2;i++) skipped[i]=-1;
    Int_t id = (Int_t)track->GetID();
    if(!(id<0)) skipped[0] = id;
    
    vertexer.SetSkipTracks(1,skipped);
    AliESDVertex *vtxESDNew = vertexer.FindPrimaryVertex(fAOD);
    
    if(!vtxESDNew) return 0;
    if(vtxESDNew->GetNContributors()<=0) {
        delete vtxESDNew; vtxESDNew=NULL;
        return 0;
    }
    
    // convert to AliAODVertex
    Double_t pos[3],cov[6],chi2perNDF;
    vtxESDNew->GetXYZ(pos); // position
    vtxESDNew->GetCovMatrix(cov); //covariance matrix
    chi2perNDF = vtxESDNew->GetChi2toNDF();
    delete vtxESDNew; vtxESDNew=NULL;
    
    AliAODVertex *vtxAODNew = new AliAODVertex(pos,cov,chi2perNDF);
    
    return vtxAODNew;
}

//=================================================================================================================================

void AliAnalysispp13TeVEMCalHFEMCReco::GetTrkClsEtaPhiDiff(AliAODTrack *t, AliAODCaloCluster *v, Double_t &phidiff, Double_t &etadiff)
{
    // Calculate phi and eta difference between a track and a cluster. The position of the track is obtained on the EMCAL surface
    
    phidiff = 999;
    etadiff = 999;
    
    if (!t||!v) return;
    
    Double_t veta = t->GetTrackEtaOnEMCal();
    Double_t vphi = t->GetTrackPhiOnEMCal();
    
    Float_t pos[3] = {0};
    v->GetPosition(pos);
    TVector3 cpos(pos);
    Double_t ceta     = cpos.Eta();
    Double_t cphi     = cpos.Phi();
    etadiff=veta-ceta;
    phidiff=TVector2::Phi_mpi_pi(vphi-cphi);  
    
}
//=================================================================================================================================
void AliAnalysispp13TeVEMCalHFEMCReco::TrackConvRadius(AliAODTrack* track,  Double_t &R)
{
  Int_t labelr = track->GetLabel(); 
  if(fIsMC && labelr>=0) 
  {
    AliAODMCParticle *mctrackk = dynamic_cast<AliAODMCParticle *>(fMCArray->At(labelr));
    R = TMath::Sqrt(mctrackk->Xv()*mctrackk->Xv()+mctrackk->Yv()*mctrackk->Yv());
  }

}
//=================================================================================================================================

Bool_t AliAnalysispp13TeVEMCalHFEMCReco::GetNMCPartProduced()
{
    //Get number of MC particles produced by generators.
    
    TList *lh = fMCHeader->GetCocktailHeaders();
    fNTotMCpart = 0;
    fNembMCpi0 = 0;
    fNembMCeta = 0;
    fNpureMC = 0;
    TString MCgen;
    TString embpi0("pi");
    TString embeta("eta");
    
    if(!lh){
        AliError("no MC header");
        return (0);
    }
    
    for(int igene=0; igene<lh->GetEntries(); igene++)
    {
        AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(igene);
        if(!gh) continue;
        
        MCgen =  gh->GetName();
        //   cout << "Gen name, N produced = " << gh->GetName() << ", " << gh->NProduced() << endl;
        if(igene==0) fNpureMC = gh->NProduced();  // generated by MB
        
        //   if(MCgen.Contains(embpi0))cout << MCgen << endl;
        //   if(MCgen.Contains(embeta))cout << MCgen << endl;
        
        if(MCgen.Contains(embpi0))fNembMCpi0 = fNTotMCpart;
        if(MCgen.Contains(embeta))fNembMCeta = fNTotMCpart;
        fNTotMCpart += gh->NProduced();
    }
    //  cout << "fNpureMC, fNembMCpi0, fNembMCeta, fNTotMCpart : " <<fNpureMC << ", " << fNembMCpi0 << ", " << fNembMCeta << ", " << fNTotMCpart << endl;
    
    return kTRUE;
}
//=================================================================================================================================

void AliAnalysispp13TeVEMCalHFEMCReco::GetPi0EtaWeight(THnSparse *SparseWeight)
{
    //Get pi0 and eta information for weight calculation
    
    Double_t fvalue[5] = {-999,-999,-999,-999,-999};
    
    for(int imc=0; imc< fNTotMCpart; imc++)
    {
        AliAODMCParticle *AODMCtrack = (AliAODMCParticle*)fMCArray->At(imc);
        
        fMCparticle = (AliAODMCParticle*) fMCArray->At(TMath::Abs(AODMCtrack->GetLabel()));
        Double_t fVx = fMCparticle->Xv();
        Double_t fVy = fMCparticle->Yv();
        Double_t Rconv = TMath::Sqrt(fVx*fVx+fVy*fVy);
        

        if(TMath::Abs(AODMCtrack->Eta()) > 0.6) continue;
        
        //-------Get PDG
        Int_t TrackPDG = TMath::Abs(AODMCtrack->GetPdgCode());
        if((TrackPDG != 111) && (TrackPDG != 221) && (TrackPDG != 22)) continue;
        
        Double_t fPartPDGid = -999;
        if (TrackPDG == 111) fPartPDGid = 0.2;
        if (TrackPDG == 221) fPartPDGid = 1.2;
        if (TrackPDG == 22) fPartPDGid = 2.2;
        
        Double_t fTrkPt = AODMCtrack->Pt();
        
        //-------Check if the particle is from Enhanced signal or not
        Bool_t fFromEnhance = kMB;
        if(imc >= fNpureMC)fFromEnhance = kEnhance;
        
        //------Get type of the particle
        Int_t fType = GetPi0EtaType(AODMCtrack);
        
        fvalue[0] = fTrkPt;             
        fvalue[1] = fPartPDGid;        
        fvalue[2] = fFromEnhance;      
        fvalue[3] = fType;                  
        fvalue[4] = Rconv;                            
                                                  
        SparseWeight->Fill(fvalue);
    }
}

//====================================================================================================================================

Bool_t AliAnalysispp13TeVEMCalHFEMCReco::GetNonHFEEffiULSLS(AliAODTrack *track, AliAODTrack *Assotrack, Bool_t fFlagLS, Bool_t fFlagULS, Double_t mass)
{
    
    Double_t TrkPt = track->Pt();
    
    //Track information
    Int_t iTrklabel = TMath::Abs(track->GetLabel());
    if(iTrklabel == 0) return kFALSE;
    AliAODMCParticle *MCPart = (AliAODMCParticle*)fMCArray->At(iTrklabel);
    
    if(TMath::Abs(MCPart->GetPdgCode())!=11) return kFALSE;
    Bool_t fFromMB = kTRUE;
    Int_t iMCmom=-999, MomPDG = -999, type=-1;
    Double_t MomPt =-999;
    Bool_t fNonHFE = IsNonHFE(MCPart, fFromMB, type, iMCmom, MomPDG, MomPt);
    
    //Associated partner information
    Int_t iTrkAssolabel = TMath::Abs(Assotrack->GetLabel());
    if(iTrkAssolabel == 0) return kFALSE;
    AliAODMCParticle *MCPartAsso = (AliAODMCParticle*)fMCArray->At(iTrkAssolabel);
    
    if(TMath::Abs(MCPartAsso->GetPdgCode())!=11) return kFALSE; // check origin of asso elec
    
    Bool_t fAssoFromMB = kTRUE;
    Int_t iMCAssomom=-999, AssoMomPDG = -999, fAssotype=-1;
    Double_t AssoMomPt =-999;
    Bool_t fAssoNonHFE = IsNonHFE(MCPartAsso, fAssoFromMB, fAssotype, iMCAssomom, AssoMomPDG, AssoMomPt);
    
    //cout << "Asso ele mom : " << iMCAssomom << ", " << AssoMomPDG << ", " << iMCmom << ", " << MomPDG << ", " << fIsFrmEmbPi0 << ", " << fIsFrmEmbEta << ", " << type << endl;
    
    if(!fAssoNonHFE) return kFALSE;
    if(iMCmom != iMCAssomom) return kFALSE; //ensure electron and partner comes from same mother
    
    if(fFlagLS) fNonHFePairInvmassLS->Fill(mass);
    if(fFlagULS) fNonHFePairInvmassULS->Fill(mass);
    
    if((fIsFrmEmbPi0 || fIsFrmEmbEta) && ftype==kNoMother){ //If parent e from embedded pi0/eta + NoMom
        if(fFlagLS) fNonHFeEmbInvmassLS->Fill(mass);
        if(fFlagULS) fNonHFeEmbInvmassULS->Fill(mass);
        if(fFlagLS) fNonHFeEmbWeightInvmassLS->Fill(mass, fWeight);
        if(fFlagULS) fNonHFeEmbWeightInvmassULS->Fill(mass, fWeight);
        
        if(fIsFrmEmbPi0){ //if from pi0
            if(fFlagLS) fPi0EmbInvmassLS->Fill(mass);
            if(fFlagULS) fPi0EmbInvmassULS->Fill(mass);
            if(fFlagLS) fPi0EmbWeightInvmassLS->Fill(mass, fWeight);
            if(fFlagULS) fPi0EmbWeightInvmassULS->Fill(mass, fWeight);
        }
        if(fIsFrmEmbEta){ //if from eta
            if(fFlagLS) fEtaEmbInvmassLS->Fill(mass);
            if(fFlagULS) fEtaEmbInvmassULS->Fill(mass);
            if(fFlagLS) fEtaEmbWeightInvmassLS->Fill(mass, fWeight);
            if(fFlagULS) fEtaEmbWeightInvmassULS->Fill(mass, fWeight);
        }
    }
    
    if(mass < fInvmassCut){
        if(fFlagLS){
            //new method
            if(fIsFrmEmbPi0 || fIsFrmEmbEta) {
                fRecoLSeEmbTrkPt->Fill(TrkPt);
                
                if(fIsFrmEmbPi0) {
                    fRecoLSeEmbWeightTrkPt->Fill(TrkPt,fWeightPi0);
                    fRecoPi0LSeEmbWeightTrkPt->Fill(TrkPt,fWeightPi0);
                }
                if(fIsFrmEmbEta){
                    fRecoLSeEmbWeightTrkPt->Fill(TrkPt,fWeightEta);
                    fRecoEtaLSeEmbWeightTrkPt->Fill(TrkPt,fWeightEta);
                }
            }
            
        }
        
        if(fFlagULS){
            //new method
            if(fIsFrmEmbPi0 || fIsFrmEmbEta) {
                fRecoULSeEmbTrkPt->Fill(TrkPt);
                
                if(fIsFrmEmbPi0) {
                    fRecoULSeEmbWeightTrkPt->Fill(TrkPt,fWeightPi0);
                    fRecoPi0ULSeEmbWeightTrkPt->Fill(TrkPt,fWeightPi0);
                }
                
                if(fIsFrmEmbEta){
                    fRecoULSeEmbWeightTrkPt->Fill(TrkPt,fWeightEta);
                    fRecoEtaULSeEmbWeightTrkPt->Fill(TrkPt,fWeightEta);
                }
            }
            
        }
    }
    
    return kTRUE;
}
//_________________________________________
void AliAnalysispp13TeVEMCalHFEMCReco::GetElectronFromStack()
{
    //electrons from MC array
    AliAODMCParticle *MCPart;
    AliAODMCParticle *MCPartMom;
    AliAODMCParticle *MCPartGMom;
    AliAODMCParticle *MCPartGGMom;
    
    for(Int_t imcArrayL=0; imcArrayL< fMCArray->GetEntries(); imcArrayL++){
        MCPart = (AliAODMCParticle*)fMCArray->At(imcArrayL);
        Int_t PDGcode = TMath::Abs(MCPart->GetPdgCode());
    
        Int_t iMCmom = -999, iMCgmom = -999, iMCggmom = -999, iMCgggmom = -999;
        Int_t MomPDG = -999, GMomPDG=-999, GGMomPDG=-999, GGGMomPDG=-999;
        
        Bool_t IsMCEle = kFALSE, IsMCPPEle = kFALSE, IsMCHFEle = kFALSE, IsMCDEle = kFALSE, IsMCBEle = kFALSE;

        if(TMath::Abs(MCPart->Eta()) > 0.7) continue;

        if(!(PDGcode == 11)) continue;
        IsMCEle = kTRUE;

        if(!MCPart->IsPhysicalPrimary()) continue;
        IsMCPPEle = kTRUE;
        
        fInclElePhysPriAll->Fill(MCPart->Pt());
        
        iMCmom = MCPart->GetMother();
        if(iMCmom > 0){
            MCPartMom = (AliAODMCParticle*)fMCArray->At(iMCmom);
            MomPDG = TMath::Abs(MCPartMom->GetPdgCode());
            
            iMCgmom = MCPartMom->GetMother();
            if(iMCgmom > 0){
                MCPartGMom = (AliAODMCParticle*)fMCArray->At(iMCgmom);
                GMomPDG = TMath::Abs(MCPartGMom->GetPdgCode());
               
                iMCggmom = MCPartGMom->GetMother();
                if(iMCggmom > 0){
                    MCPartGGMom = (AliAODMCParticle*)fMCArray->At(iMCggmom);
                    GGMomPDG = TMath::Abs(MCPartGGMom->GetPdgCode());
                }
            }
        }
        
        if((MomPDG>400 && MomPDG<600) || (MomPDG>4000 && MomPDG<6000)){
            fHFEPhysPriAll->Fill(MCPart->Pt());
            IsMCHFEle = kTRUE;
                
            if((MomPDG>500 && MomPDG<600) || (MomPDG>5000 && MomPDG<6000))
                IsMCBEle = kTRUE;
                
            if((GMomPDG>500 && GMomPDG<600) || (GMomPDG>5000 && GMomPDG<6000))
                IsMCBEle = kTRUE;
                
            if((GGMomPDG>500 && GGMomPDG<600) || (GGMomPDG>5000 && GGMomPDG<6000))
                IsMCBEle = kTRUE;
                
            if(IsMCBEle)fBEPhysPriAll->Fill(MCPart->Pt());
            else fDEPhysPriAll->Fill(MCPart->Pt());
        }
    }
}
//====================================================================================================================================
void AliAnalysispp13TeVEMCalHFEMCReco::GetTrackHFStatus(AliAODTrack *track, Bool_t &IsMCEle, Bool_t &IsMCPPEle, Bool_t &IsMCHFEle, Bool_t &IsMCBEle, Bool_t &IsMCDEle)
{
//Check the MC track status for electrons
    
    Int_t iTrklabel = TMath::Abs(track->GetLabel());
    AliAODMCParticle *MCPart;
    AliAODMCParticle *MCPartMom;
    AliAODMCParticle *MCPartGMom;
    AliAODMCParticle *MCPartGGMom;

    Int_t iMCmom = -999, iMCgmom = -999, iMCggmom = -999, iMCgggmom = -999;
    Int_t MomPDG = -999, GMomPDG=-999, GGMomPDG=-999, GGGMomPDG=-999;

    if(iTrklabel > 0){
        MCPart = (AliAODMCParticle*)fMCArray->At(iTrklabel);
        if(TMath::Abs(MCPart->GetPdgCode())==11){
            IsMCEle = kTRUE;
            
            if(MCPart->IsPhysicalPrimary()){
                IsMCPPEle = kTRUE;
                
                iMCmom = MCPart->GetMother();
                if(iMCmom > 0){
                    MCPartMom = (AliAODMCParticle*)fMCArray->At(iMCmom);
                    MomPDG = TMath::Abs(MCPartMom->GetPdgCode());

                    iMCgmom = MCPartMom->GetMother();
                    if(iMCgmom > 0){
                        MCPartGMom = (AliAODMCParticle*)fMCArray->At(iMCgmom);
                        GMomPDG = TMath::Abs(MCPartGMom->GetPdgCode());
                      
                        iMCggmom = MCPartGMom->GetMother();
                        if(iMCggmom > 0){
                            MCPartGGMom = (AliAODMCParticle*)fMCArray->At(iMCggmom);
                            GGMomPDG = TMath::Abs(MCPartGGMom->GetPdgCode());
                        }
                    }
                }

                if((MomPDG>400 && MomPDG<600) || (MomPDG>4000 && MomPDG<6000)){
                    IsMCHFEle = kTRUE;
                
                    if((MomPDG>500 && MomPDG<600) || (MomPDG>5000 && MomPDG<6000))
                        IsMCBEle = kTRUE;

                    if((GMomPDG>500 && GMomPDG<600) || (GMomPDG>5000 && GMomPDG<6000))
                        IsMCBEle = kTRUE;
              
                    if((GGMomPDG>500 && GGMomPDG<600) || (GGMomPDG>5000 && GGMomPDG<6000))
                        IsMCBEle = kTRUE;

                    if(!IsMCBEle) IsMCDEle = kTRUE;
                }
            }
        }
    }
}
//====================================================================================================================================
void AliAnalysispp13TeVEMCalHFEMCReco::GetEIDRecoEffi(AliAODTrack *track, AliAODCaloCluster *clust, Bool_t IsMCPPEle, Bool_t IsMCHFEle, Bool_t IsMCBEle, Bool_t IsMCDEle, Double_t fTPCnSigma)
{
    //Filling histograms for EID efficiency
    
    Bool_t PassSSCut = kFALSE;
    
    Double_t eop=-1.0, eop_NL = -1.0;
    Double_t m02 = -999,m20 = -999;
    Double_t clustE_NL = clust->GetNonLinCorrEnergy();
    Double_t clustE = clust->E();
    Double_t TrkPt = track->Pt();
    m02 =clust->GetM02();
    m20 =clust->GetM20();
    
    if(track->P()>0)eop_NL = clustE_NL/track->P();
    if(track->P()>0)eop = clustE/track->P();
    
    if(eop > fCutEopEMin && eop < fCutEopEMax){
        if(IsMCPPEle) fInclElePhysPriEovP->Fill(TrkPt);
        if(IsMCHFEle) fHFEPhysPriEovP->Fill(TrkPt);
        if(IsMCBEle) fBEPhysPriEovP->Fill(TrkPt);
        if(IsMCDEle) fDEPhysPriEovP->Fill(TrkPt);
        
        if(fTPCnSigma > fTPCnsigmin && fTPCnSigma < fTPCnsigmax){
            if(IsMCPPEle) fInclElePhysPriTPCnsig->Fill(TrkPt);
            if(IsMCHFEle) fHFEPhysPriTPCnsig->Fill(TrkPt);
            if(IsMCBEle) fBEPhysPriTPCnsig->Fill(TrkPt);
            if(IsMCDEle) fDEPhysPriTPCnsig->Fill(TrkPt);
                        
            if(TrkPt < 12.0){
                if(m02 > fM02Min && m02 < fM02Max1) PassSSCut = kTRUE;
            }
            if(TrkPt >= 12.0 && TrkPt < 20.0){
                if(m02 > fM02Min && m02 < fM02Max2) PassSSCut = kTRUE;
            }

            if(TrkPt >= 20.0){
                if(m02 > fM02Min && m02 < fM02Max3) PassSSCut = kTRUE;
            }
            
            if(PassSSCut){
                if(IsMCPPEle) fInclElePhysPriSS->Fill(TrkPt);
                if(IsMCHFEle) fHFEPhysPriSS->Fill(TrkPt);
                if(IsMCBEle) fBEPhysPriSS->Fill(TrkPt);
                if(IsMCDEle) fDEPhysPriSS->Fill(TrkPt);
            }
        }
    }
}
//====================================================================================================================================

    Bool_t AliAnalysispp13TeVEMCalHFEMCReco::FindMother(Int_t mcIndex)
    {
      fIsHFE1 = kFALSE;
      fIsHFE2 = kFALSE;
      fIsNonHFE = kFALSE;
      fIsFromD = kFALSE;
      fIsFromBarionB = kFALSE;
      fIsFromMesonB = kFALSE;
      fIsFromBarionBD =kFALSE;
      fIsFromMesonBD = kFALSE;
      fIsFromPi0 = kFALSE;
      fIsFromEta = kFALSE;
      fIsFromGamma = kFALSE;

      if(mcIndex < 0 || !fIsMC)
      {
        return kFALSE;
      }

      Int_t pdg = -99999;
      Int_t mpdg = -99999;
      Int_t gmpdg = -99999;
      Int_t ggmpdg = -99999;
      Int_t gggmpdg = -99999;

      if(fIsAOD)
      {
        fMCparticle = (AliAODMCParticle*) fMCArray->At(mcIndex);

        pdg = TMath::Abs(fMCparticle->GetPdgCode());


        if(pdg!=11)
        {
          fIsHFE1 = kFALSE;
          fIsHFE2 = kFALSE;
          fIsNonHFE = kFALSE;
          fIsFromD = kFALSE;
          fIsFromBarionB = kFALSE;
          fIsFromMesonB = kFALSE;
          fIsFromBarionBD =kFALSE;
          fIsFromMesonBD = kFALSE;
          fIsFromPi0 = kFALSE;
          fIsFromEta = kFALSE;
          fIsFromGamma = kFALSE;
          return kFALSE;
        }

        if(fMCparticle->GetMother()<0)
        {
          fIsHFE1 = kFALSE;
          fIsHFE2 = kFALSE;
          fIsNonHFE = kFALSE;
          fIsFromD = kFALSE;
          fIsFromBarionB = kFALSE;
          fIsFromMesonB = kFALSE;
          fIsFromBarionBD =kFALSE;
          fIsFromMesonBD = kFALSE;
          fIsFromPi0 = kFALSE;
          fIsFromEta = kFALSE;
          fIsFromGamma = kFALSE;
          return kFALSE;
        }

        fMCparticleMother = (AliAODMCParticle*) fMCArray->At(fMCparticle->GetMother());
        mpdg = TMath::Abs(fMCparticleMother->GetPdgCode());

        if(fMCparticleMother->GetMother()<0)
        {
          gmpdg = 0;
          ggmpdg = 0;
          gggmpdg = 0;
        }
        else
        {
          fMCparticleGMother = (AliAODMCParticle*) fMCArray->At(fMCparticleMother->GetMother());
          gmpdg = TMath::Abs(fMCparticleGMother->GetPdgCode());
          if(fMCparticleGMother->GetMother()<0)
          {
            ggmpdg = 0;
            gggmpdg = 0;
          }
          else
          {
            fMCparticleGGMother = (AliAODMCParticle*) fMCArray->At(fMCparticleGMother->GetMother());
            ggmpdg = TMath::Abs(fMCparticleGGMother->GetPdgCode());
            if(fMCparticleGGMother->GetMother()<0)
            {
              gggmpdg = 0;
            }
            else
            {
              fMCparticleGGGMother = (AliAODMCParticle*) fMCArray->At(fMCparticleGGMother->GetMother());
              gggmpdg = TMath::Abs(fMCparticleGGGMother->GetPdgCode());
            }
          }
        }
        //    cout<<fMCparticle->GetMother()<<"   "<<mpdg<<"    "<<gmpdg<<"    "<<ggmpdg<<"    "<<gggmpdg<<endl;
      }
    
      ///Tag Electron Source
      if(mpdg==111 || mpdg==221 || mpdg==22)
      {
        fIsHFE1 = kFALSE;
        fIsHFE2 = kFALSE;
        fIsNonHFE = kTRUE;
        fIsFromD = kFALSE;
        fIsFromBarionB = kFALSE;
        fIsFromMesonB = kFALSE;
        fIsFromBarionBD =kFALSE;
        fIsFromMesonBD = kFALSE;

        fIsFromPi0 = kFALSE;
        fIsFromEta = kFALSE;
        fIsFromGamma = kFALSE;

        if(mpdg==111) fIsFromPi0 = kTRUE;  //changed by sudhir for tagg eff 24 FeB 2019
        if(mpdg==221)fIsFromEta = kTRUE; //changed by sudhir for tagg eff 24 FeB 2019
        if(mpdg==22) fIsFromGamma = kTRUE; //changed by sudhir for tagg eff 24 FeB 2019

        return kTRUE;
      }
      else
      {
        fIsHFE1 = kTRUE;

        fIsFromPi0 = kFALSE;
        fIsFromEta = kFALSE;
        fIsFromGamma = kFALSE;

        fIsNonHFE = kFALSE;

        fIsFromD = kFALSE;
        fIsFromBarionB = kFALSE;
        fIsFromMesonB = kFALSE;
        fIsFromBarionBD =kFALSE;
        fIsFromMesonBD = kFALSE;

        if((mpdg>400 && mpdg<500) || (mpdg>4000 && mpdg<5000)) //charmed mesons and baryons
        {
          if((gmpdg>500 && gmpdg<600) || (ggmpdg>500 && ggmpdg<600) || (gggmpdg>500 && gggmpdg<600)) //when the charm comes from beauty meson
          {
            fIsHFE1 = kTRUE;
            fIsFromD = kFALSE;
            fIsFromBarionB = kFALSE;
            fIsFromMesonB = kFALSE;
            fIsFromBarionBD = kFALSE;
            fIsFromMesonBD = kTRUE;
            return kTRUE;
          }
          else if((gmpdg>5000 && gmpdg<6000) || (ggmpdg>5000 && ggmpdg<6000) || (gggmpdg>5000 && gggmpdg<6000)) //when the charm comes from beauty barion
          {
            fIsHFE1 = kTRUE;
            fIsFromD = kFALSE;
            fIsFromBarionB = kFALSE;
            fIsFromMesonB = kFALSE;
            fIsFromBarionBD = kTRUE;
            fIsFromMesonBD = kFALSE;
            return kTRUE;
          }

          else
          {
            fIsHFE1 = kTRUE;
            fIsFromD = kTRUE;
            fIsFromBarionB = kFALSE;
            fIsFromMesonB = kFALSE;
            fIsFromBarionBD =kFALSE;
            fIsFromMesonBD = kFALSE;
            return kTRUE;
          }
        }
        else if((mpdg>500 && mpdg<600)) //beauty mesons 
        {
          fIsHFE1 = kTRUE;
          fIsFromD = kFALSE;
          fIsFromBarionB = kFALSE;
          fIsFromMesonB = kTRUE;
          fIsFromBarionBD =kFALSE;
          fIsFromMesonBD = kFALSE;
          return kTRUE;
        }
        else if((mpdg>5000 && mpdg<6000)) //beauty baryons
        {
          fIsHFE1 = kTRUE;
          fIsFromD = kFALSE;
          fIsFromBarionB = kTRUE;
          fIsFromMesonB = kFALSE;
          fIsFromBarionBD =kFALSE;
          fIsFromMesonBD = kFALSE;
          return kTRUE;
        }

        else
        {
          fIsHFE1 = kFALSE;
          fIsFromD = kFALSE;
          fIsFromBarionB = kFALSE;
          fIsFromMesonB = kFALSE;
          fIsFromBarionBD =kFALSE;
          fIsFromMesonBD = kFALSE;
          return kFALSE;
        }
      }
    }

//====================================================================================================================================
Int_t AliAnalysispp13TeVEMCalHFEMCReco::GetElecSourceType(AliAODMCParticle *electron,Double_t &ptm)
{
    //
    // Return what type of gammax it is
    //
    
    // Mother
    Int_t motherlabel = electron->GetMother();
    if(motherlabel<0) return kNoMotherE;
    else {
        
        AliAODMCParticle *mother = (AliAODMCParticle*)fMCArray->At(motherlabel);
        Int_t motherpdg = TMath::Abs(mother->GetPdgCode());
        ptm=mother->Pt();
        if(motherpdg == 111) {
            Int_t typepi0eta = GetPi0EtaType(mother);
            if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother))  return kPi0NoFeedDown;
        }
        if(motherpdg == 221) {
            Int_t typepi0eta = GetPi0EtaType(mother);
            if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) return kEtaNoFeedDown;
        }
        if(motherpdg == 22) {
            Int_t gmotherlabel = mother->GetMother();
            if(gmotherlabel<0) return kDirectGamma;
            else {
                AliAODMCParticle *gmother = (AliAODMCParticle*)fMCArray->At(gmotherlabel);
                ptm=gmother->Pt();
                Int_t gmotherpdg = TMath::Abs(gmother->GetPdgCode());
                if(gmotherpdg == 111) {
                    Int_t typepi0eta = GetPi0EtaType(mother);
                    if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) return kGPi0NoFeedDown;
                }
                if(gmotherpdg == 221) {
                    Int_t typepi0eta = GetPi0EtaType(mother);
                    if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) return kGEtaNoFeedDown;
                }
                if(gmotherpdg == 22) {
                    Int_t ggmotherlabel = gmother->GetMother();
                    if(ggmotherlabel<0) return kDirectGamma;
                    else {
                        AliAODMCParticle *ggmother = (AliAODMCParticle*)fMCArray->At(ggmotherlabel);
                        ptm=ggmother->Pt();
                        Int_t ggmotherpdg = TMath::Abs(ggmother->GetPdgCode());
                        if(ggmotherpdg == 111) {
                            Int_t typepi0eta = GetPi0EtaType(mother);
                            if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) return kGPi0NoFeedDown;
                        }
                        if(ggmotherpdg == 221) {
                            Int_t typepi0eta = GetPi0EtaType(mother);
                            if((typepi0eta==kNoFeedDown) || (typepi0eta==kNoMother)) return kGEtaNoFeedDown;
                        }
                    }
                }
            }
        }
    }
    
    return kOthersE;

}
//====================================================================================================================================

Int_t AliAnalysispp13TeVEMCalHFEMCReco::GetPi0EtaType(AliAODMCParticle *part)
{
    // Return the type of particle
    
    // IsPrimary
    Bool_t primMC = part->IsPrimary();
    if(!primMC) return kNotIsPrimary;
    
    // Mother
    Int_t motherlabel = part->GetMother();
    if(motherlabel<0) return kNoMother;
    
    else {
        AliAODMCParticle *mother = (AliAODMCParticle*)fMCArray->At(motherlabel);
        Int_t motherpdg = TMath::Abs(mother->GetPdgCode());
        
        if(motherpdg == 111 || motherpdg == 221 || motherpdg == 223 || motherpdg == 333 || motherpdg == 331 || motherpdg == 113 || motherpdg == 213 || motherpdg == 313 || motherpdg == 323) return kLightMesons;
        
        if ( (int(TMath::Abs(motherpdg)/100.)%10) == 5 || (int(TMath::Abs(motherpdg)/1000.)%10) == 5 ) return kBeauty;
        if ( (int(TMath::Abs(motherpdg)/100.)%10) == 4 || (int(TMath::Abs(motherpdg)/1000.)%10) == 4 ) return kCharm;
        return kNoFeedDown;
    }
}

//====================================================================================================================================

Int_t AliAnalysispp13TeVEMCalHFEMCReco::GetHFE(AliAODMCParticle *electron, TClonesArray *mcArray)
{
  Int_t motherindex=electron->GetMother(); //Getting Electron Mother
  if(motherindex<0) return kNoMother;
  AliAODMCParticle *mother = (AliAODMCParticle*)mcArray->At(motherindex);         
  Int_t motherpdg = mother->GetPdgCode(); 
  
  if ( (int(TMath::Abs(motherpdg)/100.)%10) == 5 || (int(TMath::Abs(motherpdg)/1000.)%10) == 5 ) return kBeauty;
  if ( (int(TMath::Abs(motherpdg)/100.)%10) == 4 || (int(TMath::Abs(motherpdg)/1000.)%10) == 4 ) return kCharm;
  return kOthersE;
  
}
//====================================================================================================================================

Bool_t AliAnalysispp13TeVEMCalHFEMCReco::IsHFelectronsMC(AliAODTrack *track)
{
  fMCparticle = (AliAODMCParticle*) fMCArray->At(TMath::Abs(track->GetLabel()));
  float pdg = fMCparticle->GetPdgCode();
  //Is electron:
  if(TMath::Abs(pdg) == 11)
  {
   Bool_t MotherFound = FindMother(TMath::Abs(track->GetLabel()));
    if(MotherFound)
    {
        if(fIsHFE1){return kTRUE;}
        else{return kFALSE;}
    }
    else{return kFALSE;}     
  }

    else{return kFALSE;}
}

//====================================================================================================================================

Bool_t AliAnalysispp13TeVEMCalHFEMCReco::PassEIDCuts(AliAODTrack *track, AliAODCaloCluster *clust, Bool_t &Hadtrack)
{
    //apply electron identification cuts
    Bool_t hadTrk = kFALSE;
    Double_t eop = -1.0, eop_NL = -1.0;
    Double_t m02 = -999,m20 = -999;
    Double_t clustE = clust->E();
    Double_t clustE_NL = clust->E();
    Double_t TrkPt = track->Pt();
    if(track->P()>0){
        eop = clustE/track->P();
        eop_NL = clustE_NL/track->P();
    }
    m02 =clust->GetM02();
    m20 =clust->GetM20();
    fTPCnSigma = fPidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
 
    //Hadron E/p distribution 
    if(fTPCnSigma > fTPCnSigmaHadMin && fTPCnSigma < fTPCnSigmaHadMax)
    { 
        if(TrkPt < 12.0){  
            if(m02 > fM02Min && m02 < fM02Max1) 
                {  
                    fHadEovpNL_AftEID->Fill(TrkPt,eop_NL);
                    if(eop_NL > fCutEopEMin && eop_NL < fCutEopEMax) hadTrk=kTRUE;
                }
        }
        if(TrkPt >= 12.0 && TrkPt < 20.0){
            if(m02 > fM02Min && m02 < fM02Max2)
            {   
                fHadEovpNL_AftEID->Fill(TrkPt,eop_NL);
                if(eop_NL > fCutEopEMin && eop_NL < fCutEopEMax) hadTrk=kTRUE;
            }
        }
    
        if(TrkPt >= 20.0){
            if(m02 > fM02Min && m02 < fM02Max3)
            {   
                fHadEovpNL_AftEID->Fill(TrkPt,eop_NL);
                if(eop_NL > fCutEopEMin && eop_NL < fCutEopEMax) hadTrk=kTRUE;
            }
        }

    }
    Hadtrack = hadTrk;
    
    if(fTPCnSigma < fTPCnsigmin || fTPCnSigma > fTPCnsigmax) return kFALSE;
    if(TrkPt < 12.0)                   { if(m02 < fM02Min || m02 > fM02Max1) return kFALSE; }
    if(TrkPt >= 12.0 && TrkPt < 20.0)  { if(m02 < fM02Min || m02 > fM02Max2) return kFALSE; }
    if(TrkPt >= 20.0)                  { if(m02 < fM02Min || m02 > fM02Max3) return kFALSE; }
    
    fEopNL_AftEID->Fill(TrkPt,eop_NL);  
    if(eop_NL < fCutEopEMin || eop_NL > fCutEopEMax) return kFALSE;
    
    return kTRUE;
}

//====================================================================================================================================

Bool_t AliAnalysispp13TeVEMCalHFEMCReco::GetNonHFEEffiDenom(AliAODTrack *track)
{
    //Calculate Non-HFE efficiency demoninator   
    fIsFrmEmbPi0 = kFALSE, fIsFrmEmbEta = kFALSE;                     
    ftype = -1, fWeightPi0 = 1.0, fWeightEta = 1.0, fWeight=1.0;
    Bool_t fFromMB = kTRUE;
    
    Int_t MomPDG = -999, GMomPDG=-999, GGMomPDG=-999, GGGMomPDG=-999;
    Int_t iMCmom = -999, iMCgmom = -999, iMCggmom = -999, iMCgggmom = -999;
    Double_t MomPt =-999.0;
    
    AliAODMCParticle *MCPart = 0;
    AliAODMCParticle *MCPartMom = 0;
    AliAODMCParticle *MCPartGMom = 0;
    AliAODMCParticle *MCPartGGMom = 0;
    AliAODMCParticle *MCPartGGGMom = 0;
   

    Double_t TrkPt = track->Pt();
    Int_t iTrklabel = TMath::Abs(track->GetLabel());
    if(iTrklabel == 0) return kFALSE;
    
    MCPart = (AliAODMCParticle*)fMCArray->At(iTrklabel);
    if(TMath::Abs(MCPart->GetPdgCode())!=11) return kFALSE;
    fRealInclsElecPt->Fill(TrkPt);

    Double_t RconvIncl;  TrackConvRadius(track, RconvIncl);  
    fIncleConvRadius->Fill(TrkPt, RconvIncl);
 

    Bool_t fNonHFE = IsNonHFE(MCPart, fFromMB, ftype, iMCmom, MomPDG, MomPt);
    if(!fNonHFE) return kFALSE;
    fNonHFeTrkPt->Fill(TrkPt);  

    Double_t RconvN;  TrackConvRadius(track, RconvN); 
    fNonHFeConvRadius->Fill(TrkPt, RconvN);
    
    MCPartMom = (AliAODMCParticle*)fMCArray->At(iMCmom);
    iMCgmom = MCPartMom->GetMother();
    if(iMCgmom > 0){
        MCPartGMom = (AliAODMCParticle*)fMCArray->At(iMCgmom);
        GMomPDG = TMath::Abs(MCPartGMom->GetPdgCode());
        
        iMCggmom = MCPartGMom->GetMother();
        if(iMCggmom > 0){
            MCPartGGMom = (AliAODMCParticle*)fMCArray->At(iMCggmom);
            GGMomPDG = TMath::Abs(MCPartGGMom->GetPdgCode());
            
            iMCgggmom = MCPartGGMom->GetMother();
            if(iMCgggmom > 0){
                MCPartGGGMom = (AliAODMCParticle*)fMCArray->At(iMCgggmom);
                GGGMomPDG = TMath::Abs(MCPartGGGMom->GetPdgCode());
            }
        }
    }
    
    //cases to consider: eta->e, eta->pi0->e, eta->gamma->e, eta->pi0->gamma->e, pi0->e, pi0->gamma->e
    if(MomPDG == 221){
        if(iMCmom >= fNembMCeta && iMCmom < fNTotMCpart) { //from eta event
            fIsFrmEmbEta = kTRUE; //eta->e
            fWeightEta = fEtaWeight->Eval(MCPartMom->Pt());
        }
    }
    
    if(MomPDG == 111) {
        if(iMCmom >= fNembMCpi0 && iMCmom < fNembMCeta){ //from pi0 event
            fIsFrmEmbPi0 = kTRUE; //pi0 -> e
            fWeightPi0 = fPi0Weight->Eval(MCPartMom->Pt());
        }
        
        if(GMomPDG == 221){
            if(iMCgmom >= fNembMCeta && iMCgmom < fNTotMCpart) { //from eta event
                fIsFrmEmbEta = kTRUE; //eta->pi0-> e
                fWeightEta = fEtaWeight->Eval(MCPartGMom->Pt());
            }
        }
    }
    
    if(MomPDG == 22){
        if(GMomPDG == 221){
            if(iMCgmom >= fNembMCeta && iMCgmom < fNTotMCpart) { //from eta event
                fIsFrmEmbEta = kTRUE; //eta->gamma-> e
                fWeightEta = fEtaWeight->Eval(MCPartGMom->Pt());
            }
        }
        
        if(GMomPDG == 111){
            if(iMCgmom >= fNembMCpi0 && iMCgmom < fNembMCeta) { //from pi0 event
                fIsFrmEmbPi0 = kTRUE; //pi0-> gamma-> e
                fWeightPi0 = fPi0Weight->Eval(MCPartGMom->Pt());
            }
            
            if(GGMomPDG == 221){
                if(iMCggmom >= fNembMCeta && iMCggmom < fNTotMCpart) { //from eta event
                    fIsFrmEmbEta = kTRUE; //eta->pi0->gamma-> e
                    fWeightEta = fEtaWeight->Eval(MCPartGGMom->Pt());
                }
            }
        }
    }
    
    //   cout << "PDG of M, GM, GGM, GGGM of ele: "<< MomPDG << ", " << GMomPDG << ", " << GGMomPDG << ", " << GGGMomPDG << endl;
    //   cout << "==============" <<endl;
    Double_t prodR = TMath::Sqrt(fMCparticle->Xv()*fMCparticle->Xv()+fMCparticle->Yv()*fMCparticle->Yv());

    if(fIsFrmEmbPi0 || fIsFrmEmbEta){
        fNonHFeEmbTrkPt->Fill(TrkPt);  
        
        fNonHFeEmbTrkRConv->Fill(track->Pt(),prodR); 
        
        if(fIsFrmEmbPi0) {
            fWeight = fWeightPi0; 
            fPi0eEmbWeightTrkPt->Fill(TrkPt,fWeightPi0);   
            fNonHFeEmbWeightTrkPt->Fill(TrkPt,fWeightPi0); 
                              
            fPi0eEmbWeightTrkRConv->Fill(track->Pt(), prodR, fWeightPi0);   
            fNonHFeEmbWeightTrkRConv->Fill(track->Pt(), prodR, fWeightPi0); 

        }
        if(fIsFrmEmbEta){
            fWeight = fWeightEta;                              
            fEtaeEmbWeightTrkPt->Fill(TrkPt,fWeightEta);
            fNonHFeEmbWeightTrkPt->Fill(TrkPt,fWeightEta);  
            
            fEtaeEmbWeightTrkRConv->Fill(track->Pt(), prodR, fWeightEta);
            fNonHFeEmbWeightTrkRConv->Fill(track->Pt(), prodR, fWeightEta); 

        }
    }
    
    return kTRUE;
}

//====================================================================================================================================


void AliAnalysispp13TeVEMCalHFEMCReco::SelectPhotonicElectron(Int_t itrack, AliAODTrack *track, Bool_t &fFlagPhotonicElec, Int_t iMC)
{
    ///////////////////////////////////////////
    //////Non-HFE - Invariant mass method//////
    ///////////////////////////////////////////
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    const AliVVertex *pVtx = fAOD->GetPrimaryVertex();
    Double_t d0z0[2]={-999,-999}, cov[3];
    Double_t DCAxyCut = 1.0, DCAzCut = 2.0; //In paper it is 1 and 2 cm 
    
    Bool_t flagPhotonicElec = kFALSE, flagLSElec = kFALSE;
    Double_t ptAsso=-999., nsigma=-999.0;
    Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
    
    Int_t ntracks = -999;
    if(!fUseTender)ntracks = fAOD->GetNumberOfTracks();
    if(fUseTender) ntracks = fTracks_tender->GetEntries();
    
    for (Int_t jtrack = 0; jtrack < ntracks; jtrack++) {
        AliVParticle* VAssotrack = 0x0;
        if(!fUseTender) VAssotrack  = fAOD->GetTrack(jtrack);
        if(fUseTender) VAssotrack = dynamic_cast<AliAODTrack*>(fTracks_tender->At(jtrack)); //take tracks from Tender list
        
        if (!VAssotrack) {
            printf("ERROR: Could not receive track %d\n", jtrack);
            continue;
        }
        
        AliAODTrack *Assotrack = dynamic_cast<AliAODTrack*>(VAssotrack);
        AliESDtrack *eAssotrack = dynamic_cast<AliESDtrack*>(VAssotrack);
        AliAODTrack *aAssotrack = dynamic_cast<AliAODTrack*>(VAssotrack);
        
        //------reject same track
        if(jtrack==itrack) continue;
        
        Double_t mass=-999., width = -999;
        Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;
        
        nsigma = fPidResponse->NumberOfSigmasTPC(Assotrack, AliPID::kElectron);
        ptAsso = Assotrack->Pt();
        
        //------track cuts applied
        //if(fAOD) {  // Switched off to Match with Shreyasi
            if(!aAssotrack->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
            if(aAssotrack->GetTPCNcls() < fAssoTPCCluster) continue;
            if((!(aAssotrack->GetStatus()&AliESDtrack::kITSrefit)|| (!(aAssotrack->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
            
            //if(fRecalIP) RecalImpactParam(aAssotrack, d0z0, cov);
            /*if(aAssotrack->PropagateToDCA(pVtx, fAOD->GetMagneticField(), 20., d0z0, cov)) Switched off to Match with Shreyasi
    
            if(TMath::Abs(d0z0[0]) > DCAxyCut || TMath::Abs(d0z0[1]) > DCAzCut) continue;*/
        //}
        
        //-------loose cut on partner electron     
        if(ptAsso < fAssopTMin) continue;        
        if(TMath::Abs(aAssotrack->Eta())> fAssoEtarange) continue; 
        if(TMath::Abs(nsigma) > fAssoTPCnsig ) continue; 
        
        Int_t chargeAsso = Assotrack->Charge();
        Int_t charge = track->Charge();
        if(charge>0) fPDGe1 = -11;
        if(chargeAsso>0) fPDGe2 = -11;
        
        fFlagLS=kFALSE; fFlagULS=kFALSE;
        if(charge == chargeAsso) fFlagLS = kTRUE;
        if(charge != chargeAsso) fFlagULS = kTRUE;
        
        //-------define KFParticle to get mass
        AliKFParticle::SetField(fAOD->GetMagneticField());
        AliKFParticle ge1 = AliKFParticle(*track, fPDGe1);
        AliKFParticle ge2 = AliKFParticle(*Assotrack, fPDGe2);
        AliKFParticle recg(ge1, ge2);
        
        if(recg.GetNDF()<1) continue;
        Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
        if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;
        
        Int_t MassCorrect;
        MassCorrect = recg.GetMass(mass,width);
        
        if(fFlagLS && track->Pt()>1) fInvmassLSPt->Fill(track->Pt(),mass);
        if(fFlagULS && track->Pt()>1) fInvmassULSPt->Fill(track->Pt(),mass);
        
        //////////////////////////////////
        //Non-HFE efficiency calculation//
        //////////////////////////////////
        Bool_t EffiNumULSLS = kFALSE;
        if(fMCHeader ){
            EffiNumULSLS = GetNonHFEEffiULSLS(track, Assotrack, fFlagLS, fFlagULS, mass);
        }

        Double_t TrkPt = track->Pt(); Double_t RULS = -999, RLS = -999; 
        if(mass < fInvmassCut){
            if(fFlagLS){
                fLSElecPt->Fill(TrkPt);
                fLSElecDCA->Fill(TrkPt,fTrkDCA); 
                if(fIsMC){ 
                  TrackConvRadius(track, RLS);  
                  fRVsLSElecPt->Fill(TrkPt, RLS); 
                }
            }

            if(fFlagULS){
                fULSElecPt->Fill(TrkPt);
                fULSElecDCA->Fill(TrkPt,fTrkDCA); 
                if(fIsMC){ 
                  TrackConvRadius(track, RULS);
                  fRVsULSElecPt->Fill(TrkPt,RULS); 
                }
            }
        }
        
        if(mass < fInvmassCut && fFlagULS && !flagPhotonicElec)
            flagPhotonicElec = kTRUE; //Tag Non-HFE (random mass cut, not optimised)
    }
    fFlagPhotonicElec = flagPhotonicElec;
}

//====================================================================================================================================

Bool_t  AliAnalysispp13TeVEMCalHFEMCReco::IsNonHFE(AliAODMCParticle *MCPart, Bool_t &fFromMB, Int_t &type, Int_t &iMCmom, Int_t &MomPDG, Double_t &MomPt)
{
    //Is electron from pi0, eta and gamma
    
    iMCmom = MCPart->GetMother();
    AliAODMCParticle *MCPartMom = (AliAODMCParticle*)fMCArray->At(iMCmom);
    MomPDG = TMath::Abs(MCPartMom->GetPdgCode());
    MomPt = MCPartMom->Pt();
    
    if((MomPDG == 111) || (MomPDG == 221) || (MomPDG == 22)){
        if(iMCmom >= fNpureMC)fFromMB = kFALSE;
        type = GetPi0EtaType(MCPartMom);
        return kTRUE;
    }
    else return kFALSE;
}

//====================================================================================================================================

Bool_t AliAnalysispp13TeVEMCalHFEMCReco::GetNonHFEEffiRecoTag(AliAODTrack *track)
{
    Double_t TrkPt = track->Pt();
    Double_t prodR = TMath::Sqrt(fMCparticle->Xv()*fMCparticle->Xv()+fMCparticle->Yv()*fMCparticle->Yv()); 
    
    fRecoNonHFeTrkPt->Fill(TrkPt);
    if(fIsFrmEmbPi0 || fIsFrmEmbEta){
        fRecoNonHFeEmbTrkPt->Fill(TrkPt);

        fRecoNonHFeEmbRConv->Fill(TrkPt, prodR);

        if(fIsFrmEmbPi0) {
            fRecoPi0eEmbWeightTrkPt->Fill(TrkPt,fWeightPi0);
            fRecoNonHFeEmbWeightTrkPt->Fill(TrkPt,fWeightPi0);  

            fRecoPi0eEmbWeightTrkRConv->Fill(TrkPt, prodR,fWeightPi0);
            fRecoNonHFeEmbWeightTrkRConv->Fill(TrkPt, prodR,fWeightPi0);

        }
        if(fIsFrmEmbEta){
            fRecoEtaeEmbWeightTrkPt->Fill(TrkPt,fWeightEta);
            fRecoNonHFeEmbWeightTrkPt->Fill(TrkPt,fWeightEta); 

            fRecoEtaeEmbWeightTrkRConv->Fill(TrkPt, prodR,fWeightEta);
            fRecoNonHFeEmbWeightTrkRConv->Fill(TrkPt, prodR,fWeightEta);

        }
    }
    
    return kTRUE;
}

//====================================================================================================================================

void AliAnalysispp13TeVEMCalHFEMCReco::GetMCTemplateWeight()
{
    //Get histograms for D,B and Lamdac weight calculation
    
    AliAODMCParticle *MCPart;
    AliAODMCParticle *MCPartMom;
    AliAODMCParticle *MCPartGMom;
    AliAODMCParticle *MCPartGGMom;
    
    Double_t PartPt = -999;
    Int_t iMCmom = -999, iMCgmom = -999;
    Int_t MomPDG = -999, GMomPDG=-999;
    
    for(Int_t imcArrayL=0; imcArrayL< fMCArray->GetEntries(); imcArrayL++){
        MCPart = (AliAODMCParticle*)fMCArray->At(imcArrayL);
        Int_t PDGcode = TMath::Abs(MCPart->GetPdgCode());
        
        iMCmom = -999, iMCgmom = -999;
        MomPDG = -999, GMomPDG=-999;
        PartPt = -999;
        
        Bool_t IsMCHF = kFALSE, IsMCD = kFALSE, IsMCB = kFALSE, IsMCBD = kFALSE;
        
        if(TMath::Abs(MCPart->Eta()) > 0.9) continue;
        
        PartPt = MCPart->Pt();
        
        if((PDGcode>400 && PDGcode<600) || (PDGcode>4000 && PDGcode<6000)){
            IsMCHF = kTRUE;
            
            if((PDGcode>500 && PDGcode<600) || (PDGcode>5000 && PDGcode<6000)){
                IsMCB = kTRUE;
                fBHadpT->Fill(PartPt);
                
                if(PDGcode>500 && PDGcode<600) fBMesonpT->Fill(PartPt);
            }
            else{
                iMCmom = MCPart->GetMother();
                if(iMCmom > 0){
                    MCPartMom = (AliAODMCParticle*)fMCArray->At(iMCmom);
                    MomPDG = TMath::Abs(MCPartMom->GetPdgCode());
                    
                    if((MomPDG>500 && MomPDG<600) || (MomPDG>5000 && MomPDG<6000)){
                        IsMCB = kTRUE;
                        IsMCBD = kTRUE;
                        fBDHadpT->Fill(MCPartMom->Pt());
                    }
                    else{
                        iMCgmom = MCPartMom->GetMother();
                        if(iMCgmom > 0){
                            MCPartGMom = (AliAODMCParticle*)fMCArray->At(iMCgmom);
                            GMomPDG = TMath::Abs(MCPartGMom->GetPdgCode());
                            
                            if((GMomPDG>500 && GMomPDG<600) || (GMomPDG>5000 && GMomPDG<6000)){
                                IsMCB = kTRUE;
                                IsMCBD = kTRUE;
                                fBDHadpT->Fill(MCPartGMom->Pt());
                            }
                        }
                    }
                }
            }
            
            if(!IsMCB) {
                if((PDGcode>400 && PDGcode<500) || (PDGcode>4000 && PDGcode<5000)) fDHadpT->Fill(PartPt);
                if(PDGcode > 400 && PDGcode < 500) fDMesonpT->Fill(PartPt);
                if(PDGcode == 411) fDPluspT->Fill(PartPt);
                if(PDGcode == 421) fD0pT->Fill(PartPt);
                if(PDGcode == 431) fDspT->Fill(PartPt);
                if(PDGcode == 4122) fLambdaCpT->Fill(PartPt);
            }
        }
    }
}

//====================================================================================================================================

//________________________________________________________________________
Bool_t AliAnalysispp13TeVEMCalHFEMCReco::GetMCDCATemplates(AliAODTrack *track, Double_t TrkDCA)
{
    //Fill MC template histograms
    
    Int_t iTrklabel = TMath::Abs(track->GetLabel());
    Double_t  TrkPt = track->Pt();
    
    AliAODMCParticle *MCPart, *MCPartMom, *MCPartMomDummy;
    Int_t iMCmom = -999;
    Int_t MomPDG = -999, MomPDGDummy = -999;
    Double_t fvalue[3] = {-999,-999,-999};
    Int_t fpidSort = -99;
    
    Bool_t IsEle = kFALSE, IsHFEle=kFALSE, IsBEle=kFALSE, IsDEle=kFALSE;
    fWeightB=1.0, fWeightBMin=1.0, fWeightBMax=1.0;
    fWeightD=1.0, fWeightDUp=1.0, fWeightDDown=1.0;
    
    if(iTrklabel < 0) return kFALSE;
    MCPart = (AliAODMCParticle*)fMCArray->At(iTrklabel);
    
    if(!(MCPart->IsPhysicalPrimary())) return kFALSE;
    
    if(TMath::Abs(MCPart->GetPdgCode())!=11) return kFALSE;
    IsEle = kTRUE;
    
    iMCmom = MCPart->GetMother();
    if(iMCmom < 0) return kFALSE;
    
    MCPartMom = (AliAODMCParticle*)fMCArray->At(iMCmom);
    MomPDG = TMath::Abs(MCPartMom->GetPdgCode());
    
    if((MomPDG>400 && MomPDG<600) || (MomPDG>4000 && MomPDG<6000)) //D,B ->e
        IsHFEle = kTRUE;
    
    if(!IsHFEle)return kFALSE;
    
    //--------- Check if e<-B going back to first mother--------------
    Int_t jMCmomDummy = iMCmom;
    while(jMCmomDummy > 0){
        MCPartMomDummy = (AliAODMCParticle*)fMCArray->At(jMCmomDummy);
        MomPDGDummy = TMath::Abs(MCPartMomDummy->GetPdgCode());
        
        if((MomPDGDummy>500 && MomPDGDummy<600) || (MomPDGDummy>5000 && MomPDGDummy<6000)){ //B->e or B->X->e, loop stops when B is found or when there is no mother
            IsBEle = kTRUE;
            fBHadElecDCA->Fill(TrkPt,TrkDCA);
            
            if( (MomPDGDummy>500 && MomPDGDummy<600) || (MomPDGDummy>5000 && MomPDGDummy<6000) ){ // added extra baryons here
                fBMesonElecDCA->Fill(TrkPt,TrkDCA);
                fpidSort = 1; //Mom is B
                //else fpidSort = 3;
                GetBWeight(MCPartMomDummy, fWeightB, fWeightBMin, fWeightBMax);
                
            }
            if(MomPDGDummy>5000 && MomPDGDummy<6000){
                fBBaryonElecDCA->Fill(TrkPt,TrkDCA);
                fpidSort = 10; //Mom is b Baryon
            }
            
            jMCmomDummy = -1;  //break the loop
        }
        else {
            jMCmomDummy = MCPartMomDummy->GetMother();
        }
    }
    
    //--------- if not B->e then it should be D->e -------------
    if(!IsBEle){
        IsDEle = kTRUE;
        fDHadElecDCA->Fill(TrkPt,TrkDCA);
        

        //cout<<" fPID sort 2   -----------------------------------------  "<<endl;


        if ((MomPDG>400 && MomPDG<500) || (MomPDG>4000 && MomPDG<5000) ) {
            fDMesonElecDCA->Fill(TrkPt,TrkDCA);
            fpidSort = 5; //Mom is D
            GetDWeight(MCPartMom, fWeightD, fWeightDUp, fWeightDDown);
            GetDWeightPbPb(MCPartMom, MomPDG, fWeightD);

            //cout<<"  fpidSort   ===     "<<fpidSort <<endl;
        }
        
        if(MomPDG>4000 && MomPDG<5000) {
            fDBaryonElecDCA->Fill(TrkPt,TrkDCA);
            fpidSort = 9; //Mom is c Baryon
            if(MomPDG == 4122) GetDWeightPbPb(MCPartMom, MomPDG, fWeightD); //For Lc
            
        }
        if(MomPDG == 411) fpidSort = 11; //Mom is D+
        if(MomPDG == 421) fpidSort = 12; //Mom is D0
        if(MomPDG == 431) fpidSort = 13; //Mom is Ds
 
        //if(MomPDG == 413) fpidSort = 14; //Mom is D*+ 
        if(MomPDG > 431 && MomPDG < 436) fpidSort = 16; //Mom is other Ds
        if(MomPDG == 4122) fpidSort = 17; //Mom is Lambda c
        

        //if(MomPDG == 443) fpidSort = 6; //Mom is J/Psi
        

        if ((MomPDG>400 && MomPDG<500) || (MomPDG>4000 && MomPDG<5000) ) fpidSort = 2;
        
        if(MomPDG == 4122) fLambdaCElecDCA->Fill(TrkPt,TrkDCA);
        if(MomPDG == 421) fD0ElecDCA->Fill(TrkPt,TrkDCA);
    }
    
    //--------Filling Thnsparse---------------
    fvalue[0] = TrkPt;
    fvalue[1] = TrkDCA;
    fvalue[2] = fpidSort;
    fSprsTemplatesNoWeight->Fill(fvalue);
    
    if(IsBEle) {
        fSprsTemplatesWeight->Fill(fvalue, fWeightB);
        
            fSprsTemplatesWeightVar1->Fill(fvalue, fWeightBMin);
            fSprsTemplatesWeightVar2->Fill(fvalue, fWeightBMax);
        
    }
    if(IsDEle) {
        fSprsTemplatesWeight->Fill(fvalue, fWeightD);

            fSprsTemplatesWeightVar1->Fill(fvalue, fWeightDUp);
            fSprsTemplatesWeightVar2->Fill(fvalue, fWeightDDown);
        
    }
    
    return kTRUE;
}
//________________________________________________________________________
//________________________________________________________________________
void AliAnalysispp13TeVEMCalHFEMCReco::SetDmesonWeightHist(TH1 *D1, TH1 *D2, TH1 *D3)
{
    fDcent = (TH1F*)D1->Clone();
    fDUp = (TH1F*)D2->Clone();
    fDDown = (TH1F*)D3->Clone();
}
//________________________________________________________________________
void AliAnalysispp13TeVEMCalHFEMCReco::SetBmesonWeightHist(TH1 *B1, TH1 *B2, TH1 *B3)
{
    fBcent = (TH1F*)B1->Clone();
    fBMin = (TH1F*)B2->Clone();
    fBMax = (TH1F*)B3->Clone();
}
//________________________________________________________________________
void AliAnalysispp13TeVEMCalHFEMCReco::SetDmesonWeightHistPbPb(TH1 *D0, TH1 *DPlus, TH1 *Ds, TH1 *Lc)
{
    fD0 = (TH1F *)D0->Clone();
    fDPlus = (TH1F *)DPlus->Clone();
    fDs = (TH1F *)Ds->Clone();
    fLc = (TH1F *)Lc->Clone();
}
//________________________________________________________________________
void AliAnalysispp13TeVEMCalHFEMCReco::GetBWeight(AliAODMCParticle *Part, Double_t &BCentWeight, Double_t &BMinWeight, Double_t &BMaxWeight)
{
    //B meson weight
    
    Int_t bin = -999;
    Int_t binLast = -999;
    
    if(!fBcent){
        BCentWeight = 1.0;
        BMinWeight = 1.0;
        BMaxWeight = 1.0;
        return;
    }
    
    bin = fBcent->FindBin(Part->Pt());
    binLast = fBcent->FindBin(49.9);
    
    if(fBcent->IsBinUnderflow(bin)){
        BCentWeight = 1.0;
        BMinWeight = 1.0;
        BMaxWeight = 1.0;
        return;
    }
    if(Part->Pt() > 49.9) {
        BCentWeight = fBcent->GetBinContent(binLast);
        BMinWeight = fBMin->GetBinContent(binLast);
        BMaxWeight = fBMax->GetBinContent(binLast);
        return;
    }
    
    BCentWeight = fBcent->GetBinContent(bin);
    BMinWeight = fBMin->GetBinContent(bin);
    BMaxWeight = fBMax->GetBinContent(bin);
    
    return;
}
//________________________________________________________________________
void AliAnalysispp13TeVEMCalHFEMCReco::GetDWeight(AliAODMCParticle *Part, Double_t &DCentWeight, Double_t &DMinWeight, Double_t &DMaxWeight)
{
    //D meson weight
    
    Int_t bin = -999;
    Int_t binLast = -999;
    
    if(!fDcent){
        DCentWeight = 1.0;
        DMinWeight = 1.0;
        DMaxWeight = 1.0;
        return;
    }
    
    bin = fDcent->FindBin(Part->Pt());
    binLast = fDcent->FindBin(35.9);
    
    if(fDcent->IsBinUnderflow(bin)){
        DCentWeight = 1.0;
        DMinWeight = 1.0;
        DMaxWeight = 1.0;
        return;
    }
    if(Part->Pt() > 35.9) {
        DCentWeight = fDcent->GetBinContent(binLast);
        DMinWeight = fDUp->GetBinContent(binLast);
        DMaxWeight = fDDown->GetBinContent(binLast);
        return;
    }
    
    DCentWeight = fDcent->GetBinContent(bin);
    DMinWeight = fDUp->GetBinContent(bin);
    DMaxWeight = fDDown->GetBinContent(bin);
    
    return;
}
//________________________________________________________________________
void AliAnalysispp13TeVEMCalHFEMCReco::GetDWeightPbPb(AliAODMCParticle *Part, Int_t PDG, Double_t &DCentWeight)
{
    //D meson weight

    Int_t bin = -999;
    Int_t binLast = -999;
    
    Int_t binLc = -999;
    Int_t binLastLc = -999;
    
    if(!fD0){
        DCentWeight = 1.0;
        return;
    }
    
    bin = fD0->FindBin(Part->Pt());
    binLast = fD0->FindBin(35.9);

    binLc = fLc->FindBin(Part->Pt());
    binLastLc = fLc->FindBin(23.9);
    
    if(fD0->IsBinUnderflow(bin)){
        if(PDG == 421) DCentWeight = fD0->GetBinContent(1);
        if(PDG == 411) DCentWeight = fDPlus->GetBinContent(1);
        if(PDG == 431) DCentWeight = fDs->GetBinContent(1);
        if(PDG == 4122) DCentWeight = fLc->GetBinContent(1);
        return;
    }
    if(Part->Pt() > 35.9){
        if(PDG == 421) DCentWeight = fD0->GetBinContent(binLast);
        if(PDG == 411) DCentWeight = fDPlus->GetBinContent(binLast);
        if(PDG == 431) DCentWeight = fDs->GetBinContent(binLast);
        if(PDG == 4122) DCentWeight = fLc->GetBinContent(binLastLc);
        return;
    }
    
    if(PDG == 421) DCentWeight = fD0->GetBinContent(bin);
    if(PDG == 411) DCentWeight = fDPlus->GetBinContent(bin);
    if(PDG == 431) DCentWeight = fDs->GetBinContent(bin);
    if(PDG == 4122) DCentWeight = fLc->GetBinContent(binLc);
    
    return;
}
//________________________________________________________________________

//________________________________________________________________________
void AliAnalysispp13TeVEMCalHFEMCReco::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

/*  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
  
  fHistPt = dynamic_cast<TH1F*> (fOutputList->At(0));
  if (!fHistPt) {
    printf("ERROR: fHistPt not available\n");
    return;
  }
   
  TCanvas *c1 = new TCanvas("AliAnalysispp13TeVEMCalHFEMCReco","Pt",10,10,510,510);
  c1->cd(1)->SetLogy();
  fHistPt->DrawCopy("E");*/
}
