/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

///////////////////////////////////////////////////////////////////////////////////
//                                                                               //
//  Task for beauty decay electron spectra in pp collisions with TPC and EMCal   //
//  Author: Deepa Thomas and Erin Gauger (University of Texas at Austin)         //
//                                                                               //
///////////////////////////////////////////////////////////////////////////////////

#include <vector>
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "THnSparse.h"
#include "TFile.h"
#include "TString.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"

#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"

#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenPythiaEventHeader.h"

#include "AliMultSelection.h"
#include "AliPID.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliPIDResponse.h"
#include "AliCentrality.h"
#include "AliEventplane.h"

#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEMCALGeometry.h"
#include "AliVertexerTracks.h"
#include "AliESDVertex.h"

//#include "AliQnCorrectionsManager.h"

#include "AliAnalysisTaskHFEBESpectraEMC.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskHFEBESpectraEMC)
//________________________________________________________________________
AliAnalysisTaskHFEBESpectraEMC::AliAnalysisTaskHFEBESpectraEMC(const char *name)
: AliAnalysisTaskSE(name),
fVevent(0),
fESD(0),
fAOD(0),
fMCHeader(0),
fpidResponse(0),
fEMCALGeo(0),
fFlagSparse(kFALSE),
fUseTender(kTRUE),
fTracks_tender(0),
fCaloClusters_tender(0),
fMCparticle(0),
fMCArray(0),
fMultSelection(0),
fIsAnapp(kFALSE),
fFlagClsTypeEMC(kTRUE),
fFlagClsTypeDCAL(kTRUE),
fcentMim(0),
fcentMax(0),
fCentralityEstimator("V0M"),
fRecalIP(kTRUE),
fDeltaEta(0.05),
fDeltaPhi(0.05),
fTPCnSigma(-999.0),
fTPCnSigmaMin(-1),
fTPCnSigmaMax(3),
fM02Min(0.05),
fM02Max1(0.9),
fM02Max2(0.7),
fM20Min(0.0),
fM20Max(2000),
fEovPMin(0.9),
fEovPMax(1.2),
fNEle(0),
fTPCnSigmaHadMin(-10),
fTPCnSigmaHadMax(-3.5),
fInvmassCut(0.15),
fCalculateWeight(kFALSE),
fCalculateNonHFEEffi(kFALSE),
fCalculateElecRecoEffi(kFALSE),
fCalculateMCTemplWeightCalc(kFALSE),
fFillMCTemplates(kFALSE),
fNTotMCpart(0),
fNpureMC(0),
fNembMCpi0(0),
fNembMCeta(0),
fIsFrmEmbPi0(kFALSE),
fIsFrmEmbEta(kFALSE),
ftype(-1),
fWeight(1),
fWeightPi0(1),
fWeightEta(1),
fPi0Weight(0),
fEtaWeight(0),
fnBinsDCAHisto(400),
fTrkDCA(-999.0),
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
fWeightB(0),
fWeightD(0),
fOutputList(0),
fNevents(0),
fCent(0),
fMult(0),
fVtxZ(0),
fVtxX(0),
fVtxY(0),
fTrigMulti(0),
fHistClustE(0),
fHistNonLinClustE(0),
fHistClustEcent(0),
fEMCClsEtaPhi(0),
fHistClustEEG1(0),
fHistClustEEG1cent(0),
fHistClustEEG2(0),
fHistClustEEG2cent(0),
fEMCClsEtaPhiEG1(0),
fEMCClsEtaPhiEG2(0),
fHistoNCls(0),
fHistoNCells(0),
fHistoEperCell(0),
fHistoCalCell(0),
fHistoTimeEMC(0),
fNegTrkIDPt(0),
fTrkPt(0),
fTrketa(0),
fTrkphi(0),
fdEdx(0),
fTPCNpts(0),
fTPCnsig(0),
fTPCnsigMcEle(0),
fTPCnsigMcHad(0),
fTPCnsig_Pi(0),
fHistPtMatch(0),
fEMCTrkMatch(0),
fEMCTrkPt(0),
fEMCTrketa(0),
fEMCTrkphi(0),
fEMCdEdx(0),
fEMCTPCnsig(0),
fEMCTPCNpts(0),
fClsEAftMatch(0),
fNonLinClsEAftMatch(0),
fClsEtaPhiAftMatch(0),
fClsEtaPhiAftMatchEMCin(0),
fClsEtaPhiAftMatchEMCout(0),
fHistdEdxEop(0),
fHistNsigEop(0),
fHistEop(0),
fHistMcEopEle(0),
fHistMcEopHad(0),
fM20(0),
fM02(0),
fM20EovP(0),
fM02EovP(0),
fEleCanITShit(0),
fInvmassULSPt(0),
fInvmassLSPt(0),
fHFmomCorr(0),
fEMCTrkMatch_Phi(0),
fEMCTrkMatch_Eta(0),
fInclsElecPt(0),
fHadPt_AftEID(0),
fHadEovp_AftEID(0),
fHadEovpNL_AftEID(0),
fEop_AftEID(0),
fEopNL_AftEID(0),
fNElecInEvt(0),
fULSElecPt(0),
fLSElecPt(0),
fHadDCA(0),
fInclElecDCA(0),
fULSElecDCA(0),
fLSElecDCA(0),
fRealInclsElecPt(0),
fNonHFeTrkPt(0),
fMissingEmbEtaEleTrkPt(0),
fNonHFeEmbAllTypeTrkPt(0),
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
fInclElePhysPriTPCnsig(0),
fHFEPhysPriTPCnsig(0),
fBEPhysPriTPCnsig(0),
fDEPhysPriTPCnsig(0),
fInclElePhysPriEovPBfrSS(0),
fHFEPhysPriEovPBfrSS(0),
fBEPhysPriEovPBfrSS(0),
fDEPhysPriEovPBfrSS(0),
fInclElePhysPriSS(0),
fHFEPhysPriSS(0),
fBEPhysPriSS(0),
fDEPhysPriSS(0),
fInclElePhysPriEovP(0),
fHFEPhysPriEovP(0),
fBEPhysPriEovP(0),
fDEPhysPriEovP(0),
fBHadpT(0),
fBMesonpT(0),
fBDHadpT(0),
fDHadpT(0),
fDMesonpT(0),
fD0pT(0),
fDPluspT(0),
fDspT(0),
fLambdaCpT(0),
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
fSparseElectron(0),
fvalueElectron(0),
fSprsPi0EtaWeightCal(0),
fSprsTemplatesNoWeight(0),
fSprsTemplatesWeight(0),
fSprsTemplatesWeightVar1(0),
fSprsTemplatesWeightVar2(0)
{
    // Constructor
    
    fvalueElectron = new Double_t[7];
    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());
    // Output slot #0 id reserved by the base class for AOD
    // Output slot #1 writes into a TH1 container
    DefineOutput(1, TList::Class());
}
//________________________________________________________________________
AliAnalysisTaskHFEBESpectraEMC::AliAnalysisTaskHFEBESpectraEMC()
: AliAnalysisTaskSE("DefaultTask_HfeEMCQA"),
fVevent(0),
fESD(0),
fAOD(0),
fMCHeader(0),
fpidResponse(0),
fEMCALGeo(0),
fFlagSparse(kFALSE),
fUseTender(kTRUE),
fTracks_tender(0),
fCaloClusters_tender(0),
fMCparticle(0),
fMCArray(0),
fMultSelection(0),
fIsAnapp(kFALSE),
fFlagClsTypeEMC(kTRUE),
fFlagClsTypeDCAL(kTRUE),
fcentMim(0),
fcentMax(0),
fCentralityEstimator("V0M"),
fRecalIP(kTRUE),
fDeltaEta(0.05),
fDeltaPhi(0.05),
fTPCnSigma(-999.0),
fTPCnSigmaMin(-1),
fTPCnSigmaMax(3),
fM02Min(0.05),
fM02Max1(0.9),
fM02Max2(0.7),
fM20Min(0.0),
fM20Max(2000),
fEovPMin(0.9),
fEovPMax(1.2),
fNEle(0),
fTPCnSigmaHadMin(-10),
fTPCnSigmaHadMax(-3.5),
fInvmassCut(0.15),
fCalculateWeight(kFALSE),
fCalculateNonHFEEffi(kFALSE),
fCalculateElecRecoEffi(kFALSE),
fCalculateMCTemplWeightCalc(kFALSE),
fFillMCTemplates(kFALSE),
fNTotMCpart(0),
fNpureMC(0),
fNembMCpi0(0),
fNembMCeta(0),
fIsFrmEmbPi0(kFALSE),
fIsFrmEmbEta(kFALSE),
ftype(-1),
fWeight(1),
fWeightPi0(1),
fWeightEta(1),
fPi0Weight(0),
fEtaWeight(0),
fnBinsDCAHisto(400),
fTrkDCA(-999.0),
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
fWeightB(0),
fWeightD(0),
fOutputList(0),
fNevents(0),
fCent(0),
fMult(0),
fVtxZ(0),
fVtxX(0),
fVtxY(0),
fTrigMulti(0),
fHistClustE(0),
fHistNonLinClustE(0),
fHistClustEcent(0),
fEMCClsEtaPhi(0),
fHistClustEEG1(0),
fHistClustEEG1cent(0),
fHistClustEEG2(0),
fHistClustEEG2cent(0),
fEMCClsEtaPhiEG1(0),
fEMCClsEtaPhiEG2(0),
fHistoNCls(0),
fHistoNCells(0),
fHistoEperCell(0),
fHistoCalCell(0),
fHistoTimeEMC(0),
fNegTrkIDPt(0),
fTrkPt(0),
fTrketa(0),
fTrkphi(0),
fdEdx(0),
fTPCNpts(0),
fTPCnsig(0),
fTPCnsigMcEle(0),
fTPCnsigMcHad(0),
fTPCnsig_Pi(0),
fHistPtMatch(0),
fEMCTrkMatch(0),
fEMCTrkPt(0),
fEMCTrketa(0),
fEMCTrkphi(0),
fEMCdEdx(0),
fEMCTPCnsig(0),
fEMCTPCNpts(0),
fClsEAftMatch(0),
fNonLinClsEAftMatch(0),
fClsEtaPhiAftMatch(0),
fClsEtaPhiAftMatchEMCin(0),
fClsEtaPhiAftMatchEMCout(0),
fHistdEdxEop(0),
fHistNsigEop(0),
fHistEop(0),
fHistMcEopEle(0),
fHistMcEopHad(0),
fM20(0),
fM02(0),
fM20EovP(0),
fM02EovP(0),
fEleCanITShit(0),
fInvmassULSPt(0),
fInvmassLSPt(0),
fHFmomCorr(0),
fEMCTrkMatch_Phi(0),
fEMCTrkMatch_Eta(0),
fInclsElecPt(0),
fHadPt_AftEID(0),
fHadEovp_AftEID(0),
fHadEovpNL_AftEID(0),
fEop_AftEID(0),
fEopNL_AftEID(0),
fNElecInEvt(0),
fULSElecPt(0),
fLSElecPt(0),
fHadDCA(0),
fInclElecDCA(0),
fULSElecDCA(0),
fLSElecDCA(0),
fRealInclsElecPt(0),
fNonHFeTrkPt(0),
fMissingEmbEtaEleTrkPt(0),
fNonHFeEmbAllTypeTrkPt(0),
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
fInclElePhysPriTPCnsig(0),
fHFEPhysPriTPCnsig(0),
fBEPhysPriTPCnsig(0),
fDEPhysPriTPCnsig(0),
fInclElePhysPriEovPBfrSS(0),
fHFEPhysPriEovPBfrSS(0),
fBEPhysPriEovPBfrSS(0),
fDEPhysPriEovPBfrSS(0),
fInclElePhysPriSS(0),
fHFEPhysPriSS(0),
fBEPhysPriSS(0),
fDEPhysPriSS(0),
fInclElePhysPriEovP(0),
fHFEPhysPriEovP(0),
fBEPhysPriEovP(0),
fDEPhysPriEovP(0),
fBHadpT(0),
fBMesonpT(0),
fBDHadpT(0),
fDHadpT(0),
fDMesonpT(0),
fD0pT(0),
fDPluspT(0),
fDspT(0),
fLambdaCpT(0),
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
fSparseElectron(0),
fvalueElectron(0),
fSprsPi0EtaWeightCal(0),
fSprsTemplatesNoWeight(0),
fSprsTemplatesWeight(0),
fSprsTemplatesWeightVar1(0),
fSprsTemplatesWeightVar2(0)
{
    //Default constructor
    
    fvalueElectron = new Double_t[7];
    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());
    // Output slot #0 id reserved by the base class for AOD
    // Output slot #1 writes into a TH1 container
    // DefineOutput(1, TH1I::Class());
    DefineOutput(1, TList::Class());
    //DefineOutput(3, TTree::Class());
}
//________________________________________________________________________
AliAnalysisTaskHFEBESpectraEMC::~AliAnalysisTaskHFEBESpectraEMC()
{
    //Destructor
    delete fOutputList;
    delete fTracks_tender;
    delete fCaloClusters_tender;
    delete fSparseElectron;
    delete []fvalueElectron;
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
}
//________________________________________________________________________
void AliAnalysisTaskHFEBESpectraEMC::UserCreateOutputObjects()
{
    // Create histograms
    // Called once
    AliDebug(3, "Creating Output Objects");
    
    Double_t pi = TMath::Pi();
    fPi0Weight = new TF1("fPi0Weight","[0] / TMath::Power(TMath::Exp(-[1]*x - [2]*x*x) + x/[3], [4])");
    fEtaWeight = new TF1("fEtaWeight","[0] / TMath::Power(TMath::Exp(-[1]*x - [2]*x*x) + x/[3], [4])");
    
     if(fIsAnapp){
         fPi0Weight->SetParameters(3.72558e+02,-4.25395e-02,2.18681e-03,1.59658e+00,5.60917e+00);
         fEtaWeight->SetParameters(3.34121e+02,-7.09185e-02,2.04493e-03,1.59842e+00,5.43861e+00);
     }
    
    if(!fIsAnapp){
        if(fcentMax == 50){
            //fit obtained for pt : 2 to 30
            fPi0Weight->SetParameters(7.77976e+02,-1.61829e-01,2.17987e-03,1.39048e+00,4.57477e+00);
            fEtaWeight->SetParameters(3.23001e+02,-5.86318e-02,-2.14621e-04,1.90352e+00,5.26774e+00);
        }
    }
    
    ///////////////////////////
    //Histos for MC templates//
    ///////////////////////////
   /* if(fFillMCTemplates){
        TString DMesonWeightMaps, BMesonWeightMaps;
        
        DMesonWeightMaps = "alien:///alice/cern.ch/user/d/dthomas/DandBmesonpTweightCorrectionFiles/DMesonpTWeight.root";
        BMesonWeightMaps = "alien:///alice/cern.ch/user/d/dthomas/DandBmesonpTweightCorrectionFiles/BMesonpTWeight.root";
        
        //   printf("\n### reading file %s ...\n",DMesonWeightMaps.Data());
        //   printf("\n### reading file %s ...\n",BMesonWeightMaps.Data());
        
        TFile* f2 = TFile::Open(DMesonWeightMaps.Data());
        if(f2){
            TH1D *D1 = (TH1D*)f2->Get("RatD0");
            TH1D *D2 = (TH1D*)f2->Get("RatD0Up");
            TH1D *D3 = (TH1D*)f2->Get("RatD0Down");
            
            SetDmesonWeightHist(D1,D2,D3);
        }
        f2->Close();
        TFile* f3 = TFile::Open(BMesonWeightMaps.Data());
        if(f3){
            TH1D *B1 = (TH1D*)f3->Get("RatBMes");
            TH1D *B2 = (TH1D*)f3->Get("RatBMesMin");
            TH1D *B3 = (TH1D*)f3->Get("RatBMesMax");
            
            SetBmesonWeightHist(B1,B2,B3);
        }
        f3->Close();
    }
    */
    /////////////////////////////////////////////////
    //Automatic determination of the analysis mode//
    ////////////////////////////////////////////////
    AliVEventHandler *inputHandler = dynamic_cast<AliVEventHandler *>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    if(!TString(inputHandler->IsA()->GetName()).CompareTo("AliAODInputHandler")){
        SetAODAnalysis();
    } else {
        SetESDAnalysis();
    }
    printf("Analysis Mode: %s Analysis\n", IsAODanalysis() ? "AOD" : "ESD");
    
    ////////////////
    //Output list//
    ///////////////
    fOutputList = new TList();
    fOutputList->SetOwner();
    
    fNevents = new TH1F("fNevents","No of events",3,-0.5,2.5);
    fOutputList->Add(fNevents);
    fNevents->GetYaxis()->SetTitle("counts");
    fNevents->GetXaxis()->SetBinLabel(1,"All");
    fNevents->GetXaxis()->SetBinLabel(2,"With >2 Trks");
    fNevents->GetXaxis()->SetBinLabel(3,"Vtx_{z}<10cm");
    
    fCent = new TH1F("fCent","Centrality",100,0,100);
    fOutputList->Add(fCent);
    
    fMult = new TH2F("fMult","Track multiplicity",100,0,100,20000,0,20000);
    fOutputList->Add(fMult);
    
    fVtxZ = new TH1F("fVtxZ","Z vertex position;Vtx_{z};counts",500,-25,25);
    fOutputList->Add(fVtxZ);
    
    fVtxY = new TH1F("fVtxY","Y vertex position;Vtx_{y};counts",500,-25,25);
    fOutputList->Add(fVtxY);
    
    fVtxX = new TH1F("fVtxX","X vertex position;Vtx_{x};counts",500,-25,25);
    fOutputList->Add(fVtxX);
    
    fHistClustE = new TH1F("fHistClustE", "EMCAL cluster energy distribution; Cluster E;counts", 500, 0.0, 50.0);
    fOutputList->Add(fHistClustE);
    
    fHistNonLinClustE = new TH1F("fHistNonLinClustE", "Nonlinearity corrected EMCAL cluster energy distribution; Cluster E;counts", 500, 0.0, 50.0);
    fOutputList->Add(fHistNonLinClustE);
    
    fHistClustEcent = new TH2F("fHistClustEcent", "EMCAL cluster energy distribution vs. centrality; centrality; Cluster E", 100,0,100,500, 0.0, 50.0);
    fOutputList->Add(fHistClustEcent);
    
    fEMCClsEtaPhi = new TH2F("fEMCClsEtaPhi","EMCAL cluster #eta and #phi distribution;#eta;#phi",100,-0.9,0.9,200,0,6.3);
    fOutputList->Add(fEMCClsEtaPhi);
    
    fHistoNCls = new TH1F("fHistoNCls","No of EMCAL cluster in the event;N^{EMC}_{cls};counts",150,0,150);
    fOutputList->Add(fHistoNCls);
    
    fHistoNCells = new TH2F("fHistoNCells","No of EMCAL cells in a cluster;Cluster E;N^{EMC}_{cells}",500,0,50,30,0,30);
    fOutputList->Add(fHistoNCells);
    
    fHistoEperCell = new TH2F("fHistoEperCell","E/cell;Cluster E;E/cell",400,0,40,300,0,30);
    fOutputList->Add(fHistoEperCell);
    
    fHistoCalCell = new TH2F("fHistoCalCell","Energy of EMCAL cells;cell ID;E (GeV)",15000,-0.5,14999.5,150,0,30);
    fOutputList->Add(fHistoCalCell);
    
    fHistoTimeEMC = new TH2F("fHistoTimeEMC","EMCAL Time;E (GeV); t(ns)",500,0,50,1800,-900,900);
    fOutputList->Add(fHistoTimeEMC);
    
    fNegTrkIDPt = new TH1F("fNegTrkIDPt", "p_{T} distribution of tracks with negative track id;p_{T} (GeV/c);counts", 500, 0.0, 50.0);
    fOutputList->Add(fNegTrkIDPt);
    
    fTrkPt = new TH1F("fTrkPt","p_{T} distribution of all tracks;p_{T} (GeV/c);counts",500,0,100);
    fOutputList->Add(fTrkPt);
    
    fTrketa = new TH1F("fTrketa","All Track #eta distribution;#eta;counts",100,-1.5,1.5);
    fOutputList->Add(fTrketa);
    
    fTrkphi = new TH1F("fTrkphi","All Track #phi distribution;#phi;counts",100,0,6.3);
    fOutputList->Add(fTrkphi);
    
    fdEdx = new TH2F("fdEdx","All Track dE/dx distribution;p (GeV/c);dE/dx",500,0,50,500,0,160);
    fOutputList->Add(fdEdx);
    
    fTPCNpts = new TH2F("fTPCNpts","All track TPC Npoints used for dE/dx calculation;p (GeV/c);N points",200,0,20,200,0.,200.);
    fOutputList->Add(fTPCNpts);
    
    fTPCnsig = new TH2F("fTPCnsig","All Track TPC Nsigma distribution;p (GeV/c);#sigma_{TPC-dE/dx}",500,0,50,200,-10,10);
    fOutputList->Add(fTPCnsig);
    
    fTPCnsigMcEle = new TH2F("fTPCnsigMcEle","All Track TPC Nsigma distribution (MC electron);p (GeV/c);#sigma_{TPC-dE/dx}",500,0,50,200,-10,10);
    fOutputList->Add(fTPCnsigMcEle);
    
    fTPCnsigMcHad = new TH2F("fTPCnsigMcHad","All Track TPC Nsigma distribution (MC hadron);p (GeV/c);#sigma_{TPC-dE/dx}",500,0,50,200,-10,10);
    fOutputList->Add(fTPCnsigMcHad);
    
    fTPCnsig_Pi = new TH2F("fTPCnsig_Pi","All Track TPC Nsigma distribution wrt pion;p (GeV/c);#sigma_{TPC-dE/dx}",500,0,50,200,-10,10);
    fOutputList->Add(fTPCnsig_Pi);
    
    fHistPtMatch = new TH1F("fHistPtMatch", "p_{T} distribution of tracks matched to EMCAL;p_{T} (GeV/c);counts",500, 0.0, 50.0);
    fOutputList->Add(fHistPtMatch);
    
    fEMCTrkMatch = new TH2F("fEMCTrkMatch","Distance of EMCAL cluster to its closest track;#phi;z",100,-0.3,0.3,100,-0.3,0.3);
    fOutputList->Add(fEMCTrkMatch);
    
    fEMCTrkMatch_Phi = new TH2F("fEMCTrkMatch_Phi","Distance of EMCAL cluster to its closest track in #Delta#phi vs p_{T};p_{T};#Delta#phi",500,0,50.0,100,-0.3,0.3);
    fOutputList->Add(fEMCTrkMatch_Phi);
    
    fEMCTrkMatch_Eta = new TH2F("fEMCTrkMatch_Eta","Distance of EMCAL cluster to its closest track in #Delta#eta vs p_{T};p_{T};#Delta#eta",500,0,50.0,100,-0.3,0.3);
    fOutputList->Add(fEMCTrkMatch_Eta);
    
    fEMCTrkPt = new TH1F("fEMCTrkPt","p_{T} distribution of tracks with EMCAL cluster;p_{T} (GeV/c);counts",500,0,50);
    fOutputList->Add(fEMCTrkPt);
    
    fEMCTrketa = new TH1F("fEMCTrketa","#eta distribution of tracks matched to EMCAL;#eta;counts",60,-1.5,1.5);
    fOutputList->Add(fEMCTrketa);
    
    fEMCTrkphi = new TH1F("fEMCTrkphi","#phi distribution of tracks matched to EMCAL;#phi;counts",100,0,6.3);
    fOutputList->Add(fEMCTrkphi);
    
    fEMCdEdx = new TH2F("fEMCdEdx","dE/dx distribution of tracks matched to EMCAL;p (GeV/c);dE/dx",200,0,20,500,0,160);
    fOutputList->Add(fEMCdEdx);
    
    fEMCTPCnsig = new TH2F("fEMCTPCnsig","TPC Nsigma distribution of tracks matched to EMCAL;p (GeV/c);#sigma_{TPC-dE/dx}",500,0,50,200,-10,10);
    fOutputList->Add(fEMCTPCnsig);
    
    fEMCTPCNpts = new TH2F("fEMCTPCNpts","TPC Npoints used for dE/dx for tracks matched to EMCAL;p (GeV/c);N points",200,0,20,200,0.,200.);
    fOutputList->Add(fEMCTPCNpts);
    
    fClsEAftMatch = new TH1F("fClsEAftMatch", "EMCAL cluster energy distribution after track matching; Cluster E;counts", 500, 0.0, 50.0);
    fOutputList->Add(fClsEAftMatch);
    
    fNonLinClsEAftMatch = new TH1F("fNonLinClsEAftMatch", "Nonlinearity corrected EMCAL cluster energy distribution after track matching; Cluster E;counts", 500, 0.0, 50.0);
    fOutputList->Add(fNonLinClsEAftMatch);
    
    fClsEtaPhiAftMatch = new TH2F("fClsEtaPhiAftMatch","EMCAL cluster #eta and #phi distribution after track matching;#eta;#phi",100,-0.9,0.9,200,0,6.3);
    fOutputList->Add(fClsEtaPhiAftMatch);
    
    fClsEtaPhiAftMatchEMCin = new TH2F("fClsEtaPhiAftMatchEMCin","EMCAL cluster #eta and #phi distribution after track matching inside EMC #phi acceptence;#eta;#phi",100,-0.9,0.9,200,0,6.3);
    fOutputList->Add(fClsEtaPhiAftMatchEMCin);
    
    fClsEtaPhiAftMatchEMCout = new TH2F("fClsEtaPhiAftMatchEMCout","EMCAL cluster #eta and #phi distribution after track matching outside EMC #phi acceptence;#eta;#phi",100,-0.9,0.9,200,0,6.3);
    fOutputList->Add(fClsEtaPhiAftMatchEMCout);
    
    fHistEop = new TH2F("fHistEop", "E/p distribution;p_{T} (GeV/c);E/p", 200,0,20,40, 0.0, 2.0);
    fOutputList->Add(fHistEop);
    
    fHistMcEopEle = new TH2F("fHistMcEopEle", "E/p distribution (MC electron);p_{T} (GeV/c);E/p", 200,0,20,40, 0.0, 2.0);
    fOutputList->Add(fHistMcEopEle);
    
    fHistMcEopHad = new TH2F("fHistMcEopHad", "E/p distribution (MC hadron);p_{T} (GeV/c);E/p", 200,0,20,40, 0.0, 2.0);
    fOutputList->Add(fHistMcEopHad);
    
    fHistdEdxEop = new TH2F("fHistdEdxEop", "E/p vs dE/dx;E/p;dE/dx", 40, 0.0, 2.0, 500,0,160);
    fOutputList->Add(fHistdEdxEop);
    
    fHistNsigEop = new TH2F ("fHistNsigEop", "E/p vs TPC nsig",40, 0.0, 2.0, 200, -10,10);
    fOutputList->Add(fHistNsigEop);
    
    fM20 = new TH2F ("fM20","M20 vs pt distribution",500,0,50,200,0,2);
    fOutputList->Add(fM20);
    
    fM02 = new TH2F ("fM02","M02 vs pt distribution",500,0,50,200,0,2);
    fOutputList->Add(fM02);
    
    fM20EovP = new TH2F ("fM20EovP","M20 vs E/p distribution;E/p;M20",40,0,2,200,0,2);
    fOutputList->Add(fM20EovP);
    
    fM02EovP = new TH2F ("fM02EovP","M02 vs E/p distribution;E/p;M02",40,0,2,200,0,2);
    fOutputList->Add(fM02EovP);
    
    fEleCanITShit = new TH1F("fEleCanITShit","ITS hit map;ITS layer;counts",7,-0.5,6.5);
    fOutputList->Add(fEleCanITShit);
    
    fInclsElecPt = new TH1F("fInclsElecPt","p_{T} distribution of inclusive electrons;p_{T} (GeV/c);counts",500,0,50);
    fOutputList->Add(fInclsElecPt);
    
    fHadPt_AftEID = new TH1F("fHadPt_AftEID","p_{T} distribution of hadrons after Eid cuts;p_{T} (GeV/c);counts",500,0,50);
    fOutputList->Add(fHadPt_AftEID);
    
    fHadEovp_AftEID = new TH2F("fHadEovp_AftEID", "E/p distribution for hadrons -10<nsig<-3.5, SS cuts;p_{T} (GeV/c);E/p", 60,0,30,100, 0.0, 2.0);
    fOutputList->Add(fHadEovp_AftEID);
    
    fHadEovpNL_AftEID = new TH2F("fHadEovpNL_AftEID", "E/p distribution for hadrons -10<nsig<-3.5, NonLinearE, SS cuts;p_{T} (GeV/c);E/p", 60,0,30,100, 0.0, 2.0);
    fOutputList->Add(fHadEovpNL_AftEID);
    
    fEop_AftEID = new TH2F("fEop_AftEID", "E/p distribution after nsig, SS cuts;p_{T} (GeV/c);E/p", 60,0,30,100, 0.0, 2.0);
    fOutputList->Add(fEop_AftEID);
    
    fEopNL_AftEID = new TH2F("fEopNL_AftEID", "E/p distribution after nsig, SS cuts, NonLinearE;p_{T} (GeV/c);E/p", 60,0,30,100, 0.0, 2.0);
    fOutputList->Add(fEopNL_AftEID);
    
    fNElecInEvt = new TH1F("fNElecInEvt","No of electrons in the event; N^{ele};counts",20,-0.5,19.5);
    fOutputList->Add(fNElecInEvt);
    
    fInvmassLSPt = new TH2F("fInvmassLSPt", "Invmass of LS (e,e) for pt^{e}>1; p_{T}(GeV/c); mass(GeV/c^2); counts;", 500,0,50,500,0,1.0);
    fOutputList->Add(fInvmassLSPt);
    
    fInvmassULSPt = new TH2F("fInvmassULSPt", "Invmass of ULS (e,e) for pt^{e}>1; p_{T}(GeV/c); mass(GeV/c^2); counts;", 500,0,50,500,0,1.0);
    fOutputList->Add(fInvmassULSPt);
    
    fULSElecPt  = new TH1F("fULSElecPt","p_{T} distribution of ULS electrons;p_{T} (GeV/c);counts",500,0,50);
    fOutputList->Add(fULSElecPt);
    
    fLSElecPt= new TH1F("fLSElecPt","p_{T} distribution of LS electrons;p_{T} (GeV/c);counts",500,0,50);
    fOutputList->Add(fLSElecPt);
    
    fHadDCA = new TH2F("fHadDCA","Hadron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., fnBinsDCAHisto,-0.4,0.4);
    fOutputList->Add(fHadDCA);
    
    fInclElecDCA = new TH2F("fInclElecDCA","Inclusive electron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., fnBinsDCAHisto,-0.4,0.4);
    fOutputList->Add(fInclElecDCA);
    
    fULSElecDCA = new TH2F("fULSElecDCA","ULS electron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., fnBinsDCAHisto,-0.4,0.4);
    fOutputList->Add(fULSElecDCA);
    
    fLSElecDCA = new TH2F("fLSElecDCA","LS electron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., fnBinsDCAHisto,-0.4,0.4);
    fOutputList->Add(fLSElecDCA);
    
    if(fFlagSparse){
        Int_t bins[7]=      {232, 160, 40, 200, 200, 20, 40}; //pT;nSigma;eop;m20;m02;iSM;eopNL
        Double_t xmin[7]={2,  -8,   0,   0,   0, 0, 0};
        Double_t xmax[7]={60,   8,   2,   2,   2, 20, 2};
        fSparseElectron = new THnSparseD ("Electron","Electron;pT;nSigma;eop;m20;m02;iSM;eop_NL;",7,bins,xmin,xmax);
        fOutputList->Add(fSparseElectron);
    }
    
    if(fCalculateWeight){
        Int_t bin[4] = {500,3,2,7}; //pT, PDG, EnhancedSigOrNot, pi0etaType
        Double_t xmin[4] = {0,0,0,-1};
        Double_t xmax[4] = {50,3,2,6};
    
        fSprsPi0EtaWeightCal = new THnSparseD("fSprsPi0EtaWeightCal","Sparse to calculate #pi^{0} and #eta weight;p_{T};PDG ID;EnhanceSigOrNot;pi0etaType;",4,bin,xmin,xmax);
        fOutputList->Add(fSprsPi0EtaWeightCal);
    }
    
    if(fCalculateNonHFEEffi){
        fRealInclsElecPt = new TH1F("fRealInclsElecPt","p_{T} distribution of MC tagged inclusive electrons;p_{T} (GeV/c);counts",250,0,25);
        fOutputList->Add(fRealInclsElecPt);
        
        fNonHFeTrkPt = new TH1F("fNonHFeTrkPt","Non-HF electrons from all generators;p_{T} (GeV/c);counts",250,0,25);
        fNonHFeTrkPt->Sumw2();
        fOutputList->Add(fNonHFeTrkPt);
        
        fMissingEmbEtaEleTrkPt = new TH1F("fMissingEmbEtaEleTrkPt","Missing electrons from embedded #eta  + No mom ;p_{T} (GeV/c);counts",250,0,25);
        fMissingEmbEtaEleTrkPt->Sumw2();
        fOutputList->Add(fMissingEmbEtaEleTrkPt);
        
        fNonHFeEmbAllTypeTrkPt = new TH1F("fNonHFeEmbAllTypeTrkPt","Non-HF electrons from embedded #pi^{0} and #eta of all type;p_{T} (GeV/c);counts",250,0,25);
        fNonHFeEmbAllTypeTrkPt->Sumw2();
        fOutputList->Add(fNonHFeEmbAllTypeTrkPt);
        
        fNonHFeEmbTrkPt = new TH1F("fNonHFeEmbTrkPt","Non-HF electrons from embedded #pi^{0} and #eta + No mom;p_{T} (GeV/c);counts",250,0,25);
        fNonHFeEmbTrkPt->Sumw2();
        fOutputList->Add(fNonHFeEmbTrkPt);
        
        fNonHFeEmbWeightTrkPt = new TH1F("fNonHFeEmbWeightTrkPt","Non-HF electrons from embedded #pi^{0} and #eta + No mom with weight + No mom;p_{T} (GeV/c);counts",250,0,25);
        fNonHFeEmbWeightTrkPt->Sumw2();
        fOutputList->Add(fNonHFeEmbWeightTrkPt);
        
        fPi0eEmbWeightTrkPt = new TH1F("fPi0eEmbWeightTrkPt","Non-HF electrons from embedded #pi^{0} + No mom with weight;p_{T} (GeV/c);counts",250,0,25);
        fPi0eEmbWeightTrkPt->Sumw2();
        fOutputList->Add(fPi0eEmbWeightTrkPt);
        
        fEtaeEmbWeightTrkPt = new TH1F("fEtaeEmbWeightTrkPt","Non-HF electrons from embedded #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,25);
        fEtaeEmbWeightTrkPt->Sumw2();
        fOutputList->Add(fEtaeEmbWeightTrkPt);
        
        fRecoNonHFeTrkPt = new TH1F("fRecoNonHFeTrkPt"," Reco Non-HF electrons from all generators;p_{T} (GeV/c);counts",250,0,25);
        fRecoNonHFeTrkPt->Sumw2();
        fOutputList->Add(fRecoNonHFeTrkPt);

        fRecoNonHFeEmbTrkPt = new TH1F("fRecoNonHFeEmbTrkPt","Reco Non-HF electrons from embedded #pi^{0} and #eta + No mom;p_{T} (GeV/c);counts",250,0,25);
        fRecoNonHFeEmbTrkPt->Sumw2();
        fOutputList->Add(fRecoNonHFeEmbTrkPt);
        
        fRecoNonHFeEmbWeightTrkPt = new TH1F("fRecoNonHFeEmbWeightTrkPt","Reco Non-HF electrons from embedded #pi^{0} and #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,25);
        fRecoNonHFeEmbWeightTrkPt->Sumw2();
        fOutputList->Add(fRecoNonHFeEmbWeightTrkPt);
        
        fRecoPi0eEmbWeightTrkPt = new TH1F("fRecoPi0eEmbWeightTrkPt","Reco Non-HF electrons from embedded #pi^{0}  + No mom with weight;p_{T} (GeV/c);counts",250,0,25);
        fRecoPi0eEmbWeightTrkPt->Sumw2();
        fOutputList->Add(fRecoPi0eEmbWeightTrkPt);
        
        fRecoEtaeEmbWeightTrkPt = new TH1F("fRecoEtaeEmbWeightTrkPt","Reco Non-HF electrons from embedded #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,25);
        fRecoEtaeEmbWeightTrkPt->Sumw2();
        fOutputList->Add(fRecoEtaeEmbWeightTrkPt);
        
        fNonHFePairInvmassLS = new TH1F("fNonHFePairInvmassLS", "Inv mass of LS (e,e) if both e- are Non-HFE; mass(GeV/c^2); counts;",  50,0,0.5);
        fOutputList->Add(fNonHFePairInvmassLS);
        
        fNonHFePairInvmassULS = new TH1F("fNonHFePairInvmassULS", "Inv mass of ULS (e,e) if both e- are Non-HFE; mass(GeV/c^2); counts;",  50,0,0.5);
        fOutputList->Add(fNonHFePairInvmassULS);
        
        fNonHFeEmbInvmassLS = new TH1F("fNonHFeEmbInvmassLS", "Inv mass of LS (e,e) for Non-HFE from embedded #pi^{0} and #eta; mass(GeV/c^2); counts;",  50,0,0.5);
        fOutputList->Add(fNonHFeEmbInvmassLS);
        
        fNonHFeEmbInvmassULS = new TH1F("fNonHFeEmbInvmassULS", "Inv mass of ULS (e,e) for Non-HFE from embedded #pi^{0} and #eta; mass(GeV/c^2); counts;",  50,0,0.5);
        fOutputList->Add(fNonHFeEmbInvmassULS);
        
        fNonHFeEmbWeightInvmassLS = new TH1F("fNonHFeEmbWeightInvmassLS", "Inv mass of LS (e,e) for Non-HFE from embedded #pi^{0} and #eta with weight; mass(GeV/c^2); counts;",  50,0,0.5);
        fOutputList->Add(fNonHFeEmbWeightInvmassLS);
        
        fNonHFeEmbWeightInvmassULS = new TH1F("fNonHFeEmbWeightInvmassULS", "Inv mass of ULS (e,e) for Non-HFE from embedded #pi^{0} and #eta with weight; mass(GeV/c^2); counts;",  50,0,0.5);
        fOutputList->Add(fNonHFeEmbWeightInvmassULS);
        
        fPi0EmbInvmassLS = new TH1F("fPi0EmbInvmassLS", "Inv mass of LS (e,e) for ele from embedded #pi^{0}; mass(GeV/c^2); counts;",  50,0,0.5);
        fOutputList->Add(fPi0EmbInvmassLS);
        
        fPi0EmbInvmassULS  = new TH1F("fPi0EmbInvmassULS", "Inv mass of ULS (e,e) for ele from embedded #pi^{0}; mass(GeV/c^2); counts;",  50,0,0.5);
        fOutputList->Add(fPi0EmbInvmassULS);
        
        fPi0EmbWeightInvmassLS = new TH1F("fPi0EmbWeightInvmassLS", "Inv mass of LS (e,e) for ele from embedded #pi^{0} with weight; mass(GeV/c^2); counts;",  50,0,0.5);
        fOutputList->Add(fPi0EmbWeightInvmassLS);
        
        fPi0EmbWeightInvmassULS  = new TH1F("fPi0EmbWeightInvmassULS", "Inv mass of ULS (e,e) for ele from embedded #pi^{0} with weight; mass(GeV/c^2); counts;",  50,0,0.5);
        fOutputList->Add(fPi0EmbWeightInvmassULS);
        
        fEtaEmbInvmassLS = new TH1F("fEtaEmbInvmassLS", "Inv mass of LS (e,e) for ele from embedded #eta; mass(GeV/c^2); counts;",  50,0,0.5);
        fOutputList->Add(fEtaEmbInvmassLS);
        
        fEtaEmbInvmassULS = new TH1F("fEtaEmbInvmassULS", "Inv mass of ULS (e,e) for ele from embedded #eta; mass(GeV/c^2); counts;",  50,0,0.5);
        fOutputList->Add(fEtaEmbInvmassULS);
        
        fEtaEmbWeightInvmassLS = new TH1F("fEtaEmbWeightInvmassLS", "Inv mass of LS (e,e) for ele from embedded #eta with weight; mass(GeV/c^2); counts;",  50,0,0.5);
        fOutputList->Add(fEtaEmbWeightInvmassLS);
        
        fEtaEmbWeightInvmassULS  = new TH1F("fEtaEmbWeightInvmassULS", "Inv mass of ULS (e,e) for ele from embedded #eta with weight; mass(GeV/c^2); counts;",  50,0,0.5);
        fOutputList->Add(fEtaEmbWeightInvmassULS);
        
        fRecoLSeEmbTrkPt  = new TH1F("fRecoLSeEmbTrkPt","Reco LS electrons from embedded #pi^{0} and #eta + No mom;p_{T} (GeV/c);counts",250,0,25);
        fRecoLSeEmbTrkPt->Sumw2();
        fOutputList->Add(fRecoLSeEmbTrkPt);
        
        fRecoLSeEmbWeightTrkPt = new TH1F("fRecoLSeEmbWeightTrkPt","Reco LS electrons from embedded #pi^{0} and #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,25);
        fRecoLSeEmbWeightTrkPt->Sumw2();
        fOutputList->Add(fRecoLSeEmbWeightTrkPt);
        
        fRecoPi0LSeEmbWeightTrkPt = new TH1F("fRecoPi0LSeEmbWeightTrkPt","Reco LS electrons from embedded #pi^{0}  + No mom with weight;p_{T} (GeV/c);counts",250,0,25);
        fRecoPi0LSeEmbWeightTrkPt->Sumw2();
        fOutputList->Add(fRecoPi0LSeEmbWeightTrkPt);
        
        fRecoEtaLSeEmbWeightTrkPt  = new TH1F("fRecoEtaLSeEmbWeightTrkPt","Reco LS electrons from embedded #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,25);
        fRecoEtaLSeEmbWeightTrkPt->Sumw2();
        fOutputList->Add(fRecoEtaLSeEmbWeightTrkPt);
        
        fRecoULSeEmbTrkPt = new TH1F("fRecoULSeEmbTrkPt","Reco ULS electrons from embedded #pi^{0} and #eta + No mom;p_{T} (GeV/c);counts",250,0,25);
        fRecoULSeEmbTrkPt->Sumw2();
        fOutputList->Add(fRecoULSeEmbTrkPt);
        
        fRecoULSeEmbWeightTrkPt = new TH1F("fRecoULSeEmbWeightTrkPt","Reco ULS electrons from embedded #pi^{0} and #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,25);
        fRecoULSeEmbWeightTrkPt->Sumw2();
        fOutputList->Add(fRecoULSeEmbWeightTrkPt);
        
        fRecoPi0ULSeEmbWeightTrkPt = new TH1F("fRecoPi0ULSeEmbWeightTrkPt","Reco ULS electrons from embedded #pi^{0}  + No mom with weight;p_{T} (GeV/c);counts",250,0,25);
        fRecoPi0ULSeEmbWeightTrkPt->Sumw2();
        fOutputList->Add(fRecoPi0ULSeEmbWeightTrkPt);
        
        fRecoEtaULSeEmbWeightTrkPt = new TH1F("fRecoEtaULSeEmbWeightTrkPt","Reco ULS electrons from embedded #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,25);
        fRecoEtaULSeEmbWeightTrkPt->Sumw2();
        fOutputList->Add(fRecoEtaULSeEmbWeightTrkPt);
    }
    
    if(fCalculateElecRecoEffi){
        fInclElePhysPriAll = new TH1F("fInclElePhysPriAll","Physical primary inclusive electrons for reco effi, All;p_{T} (GeV/c);counts",250,0,25);
        fInclElePhysPriAll->Sumw2();
        fOutputList->Add(fInclElePhysPriAll);
        
        fHFEPhysPriAll = new TH1F("fHFEPhysPriAll","Physical primary HFE for reco effi, All;p_{T} (GeV/c);counts",250,0,25);
        fHFEPhysPriAll->Sumw2();
        fOutputList->Add(fHFEPhysPriAll);
        
        fBEPhysPriAll = new TH1F("fBEPhysPriAll","Physical primary b->e for reco effi, All;p_{T} (GeV/c);counts",250,0,25);
        fBEPhysPriAll->Sumw2();
        fOutputList->Add(fBEPhysPriAll);
        
        fDEPhysPriAll = new TH1F("fDEPhysPriAll","Physical primary c->e for reco effi, All;p_{T} (GeV/c);counts",250,0,25);
        fDEPhysPriAll->Sumw2();
        fOutputList->Add(fDEPhysPriAll);
        
        fInclElePhysPriTrkCuts = new TH1F("fInclElePhysPriTrkCuts","Physical primary inclusive electrons for reco effi, Aft Trk cuts;p_{T} (GeV/c);counts",250,0,25);
        fInclElePhysPriTrkCuts->Sumw2();
        fOutputList->Add(fInclElePhysPriTrkCuts);
        
        fHFEPhysPriTrkCuts = new TH1F("fHFEPhysPriTrkCuts","Physical primary HFE for reco effi, Aft Trk cuts;p_{T} (GeV/c);counts",250,0,25);
        fHFEPhysPriTrkCuts->Sumw2();
        fOutputList->Add(fHFEPhysPriTrkCuts);
        
        fBEPhysPriTrkCuts = new TH1F("fBEPhysPriTrkCuts","Physical primary b->e for reco effi, Aft Trk cuts;p_{T} (GeV/c);counts",250,0,25);
        fBEPhysPriTrkCuts->Sumw2();
        fOutputList->Add(fBEPhysPriTrkCuts);
        
        fDEPhysPriTrkCuts = new TH1F("fDEPhysPriTrkCuts","Physical primary c->e for reco effi, Aft Trk cuts;p_{T} (GeV/c);counts",250,0,25);
        fDEPhysPriTrkCuts->Sumw2();
        fOutputList->Add(fDEPhysPriTrkCuts);
        
        fInclElePhysPriEMCMatch = new TH1F("fInclElePhysPriEMCMatch","Physical primary inclusive electron for reco effi, Aft EMC match;p_{T} (GeV/c);counts",250,0,25);
        fInclElePhysPriEMCMatch->Sumw2();
        fOutputList->Add(fInclElePhysPriEMCMatch);
        
        fHFEPhysPriEMCMatch = new TH1F("fHFEPhysPriEMCMatch","Physical primary HFE for reco effi, Aft EMC match;p_{T} (GeV/c);counts",250,0,25);
        fHFEPhysPriEMCMatch->Sumw2();
        fOutputList->Add(fHFEPhysPriEMCMatch);
        
        fBEPhysPriEMCMatch = new TH1F("fBEPhysPriEMCMatch","Physical primary b->e for reco effi, Aft EMC match;p_{T} (GeV/c);counts",250,0,25);
        fBEPhysPriEMCMatch->Sumw2();
        fOutputList->Add(fBEPhysPriEMCMatch);
        
        fDEPhysPriEMCMatch = new TH1F("fDEPhysPriEMCMatch","Physical primary c->e for reco effi, Aft EMC match;p_{T} (GeV/c);counts",250,0,25);
        fDEPhysPriEMCMatch->Sumw2();
        fOutputList->Add(fDEPhysPriEMCMatch);
        
        fInclElePhysPriTPCnsig = new TH1F("fInclElePhysPriTPCnsig","Physical primary inclusive electron for reco effi, Aft TPCnsig;p_{T} (GeV/c);counts",250,0,25);
        fInclElePhysPriTPCnsig->Sumw2();
        fOutputList->Add(fInclElePhysPriTPCnsig);
        
        fHFEPhysPriTPCnsig = new TH1F("fHFEPhysPriTPCnsig","Physical primary HFE for reco effi, Aft TPCnsig;p_{T} (GeV/c);counts",250,0,25);
        fHFEPhysPriTPCnsig->Sumw2();
        fOutputList->Add(fHFEPhysPriTPCnsig);
        
        fBEPhysPriTPCnsig = new TH1F("fBEPhysPriTPCnsig","Physical primary b->e for reco effi, Aft TPCnsig;p_{T} (GeV/c);counts",250,0,25);
        fBEPhysPriTPCnsig->Sumw2();
        fOutputList->Add(fBEPhysPriTPCnsig);
        
        fDEPhysPriTPCnsig = new TH1F("fDEPhysPriTPCnsig","Physical primary c->e for reco effi, Aft TPCnsig;p_{T} (GeV/c);counts",250,0,25);
        fDEPhysPriTPCnsig->Sumw2();
        fOutputList->Add(fDEPhysPriTPCnsig);
        
        fInclElePhysPriEovPBfrSS = new TH1F("fInclElePhysPriEovPBfrSS","Physical primary inclusive electron for reco effi, Aft E/p cut & bfr SS cut;p_{T} (GeV/c);counts",250,0,25);
        fInclElePhysPriEovPBfrSS->Sumw2();
        fOutputList->Add(fInclElePhysPriEovPBfrSS);
        
        fHFEPhysPriEovPBfrSS = new TH1F("fHFEPhysPriEovPBfrSS","Physical primary HFE for reco effi, Aft E/p cut & bfr SS cut;p_{T} (GeV/c);counts",250,0,25);
        fHFEPhysPriEovPBfrSS->Sumw2();
        fOutputList->Add(fHFEPhysPriEovPBfrSS);
        
        fBEPhysPriEovPBfrSS = new TH1F("fBEPhysPriEovPBfrSS","Physical primary b->e for reco effi, Aft E/p cut & bfr SS cut;p_{T} (GeV/c);counts",250,0,25);
        fBEPhysPriEovPBfrSS->Sumw2();
        fOutputList->Add(fBEPhysPriEovPBfrSS);
        
        fDEPhysPriEovPBfrSS = new TH1F("fDEPhysPriEovPBfrSS","Physical primary c->e for reco effi, Aft E/p cut & bfr SS cut;p_{T} (GeV/c);counts",250,0,25);
        fDEPhysPriEovPBfrSS->Sumw2();
        fOutputList->Add(fDEPhysPriEovPBfrSS);
        
        
        fInclElePhysPriSS = new TH1F("fInclElePhysPriSS","Physical primary inclusive electron for reco effi, Aft SS cut;p_{T} (GeV/c);counts",250,0,25);
        fInclElePhysPriSS->Sumw2();
        fOutputList->Add(fInclElePhysPriSS);
        
        fHFEPhysPriSS = new TH1F("fHFEPhysPriSS","Physical primary HFE for reco effi, Aft SS cut;p_{T} (GeV/c);counts",250,0,25);
        fHFEPhysPriSS->Sumw2();
        fOutputList->Add(fHFEPhysPriSS);
        
        fBEPhysPriSS = new TH1F("fBEPhysPriSS","Physical primary b->e for reco effi, Aft SS cut;p_{T} (GeV/c);counts",250,0,25);
        fBEPhysPriSS->Sumw2();
        fOutputList->Add(fBEPhysPriSS);
        
        fDEPhysPriSS = new TH1F("fDEPhysPriSS","Physical primary c->e for reco effi, Aft SS cut;p_{T} (GeV/c);counts",250,0,25);
        fDEPhysPriSS->Sumw2();
        fOutputList->Add(fDEPhysPriSS);
        
        fInclElePhysPriEovP = new TH1F("fInclElePhysPriEovP","Physical primary inclusive electron for reco effi, Aft E/p cut;p_{T} (GeV/c);counts",250,0,25);
        fInclElePhysPriEovP->Sumw2();
        fOutputList->Add(fInclElePhysPriEovP);
        
        fHFEPhysPriEovP = new TH1F("fHFEPhysPriEovP","Physical primary HFE for reco effi, Aft E/p cut;p_{T} (GeV/c);counts",250,0,25);
        fHFEPhysPriEovP->Sumw2();
        fOutputList->Add(fHFEPhysPriEovP);
        
        fBEPhysPriEovP = new TH1F("fBEPhysPriEovP","Physical primary b->e for reco effi, Aft E/p cut;p_{T} (GeV/c);counts",250,0,25);
        fBEPhysPriEovP->Sumw2();
        fOutputList->Add(fBEPhysPriEovP);
        
        fDEPhysPriEovP = new TH1F("fDEPhysPriEovP","Physical primary c->e for reco effi, Aft E/p cut;p_{T} (GeV/c);counts",250,0,25);
        fDEPhysPriEovP->Sumw2();
        fOutputList->Add(fDEPhysPriEovP);
    }
    
    if(fCalculateMCTemplWeightCalc){
        fBHadpT = new TH1F("fBHadpT","B hadron pT;p_{T} (GeV/c);counts",500,0,50);
        fBHadpT->Sumw2();
        fOutputList->Add(fBHadpT);
        
        fBMesonpT = new TH1F("fBMesonpT","B meson pT;p_{T} (GeV/c);counts",500,0,50);
        fBMesonpT->Sumw2();
        fOutputList->Add(fBMesonpT);
        
        fBDHadpT = new TH1F("fBDHadpT","D (<- B) hadron pT;p_{T} (GeV/c);counts",500,0,50);
        fBDHadpT->Sumw2();
        fOutputList->Add(fBDHadpT);
        
        fDHadpT = new TH1F("fDHadpT","Prompt D hadron pT;p_{T} (GeV/c);counts",500,0,50);
        fDHadpT->Sumw2();
        fOutputList->Add(fDHadpT);
        
        fDMesonpT = new TH1F("fDMesonpT","Prompt D meson pT;p_{T} (GeV/c);counts",500,0,50);
        fDMesonpT->Sumw2();
        fOutputList->Add(fDMesonpT);

        fD0pT = new TH1F("fD0pT","Prompt D0 meson pT;p_{T} (GeV/c);counts",500,0,50);
        fD0pT->Sumw2();
        fOutputList->Add(fD0pT);
        
        fDPluspT = new TH1F("fDPluspT","Prompt D+ meson pT;p_{T} (GeV/c);counts",500,0,50);
        fDPluspT->Sumw2();
        fOutputList->Add(fDPluspT);
        
        fDspT = new TH1F("fDspT","Prompt D+s meson pT;p_{T} (GeV/c);counts",500,0,50);
        fDspT->Sumw2();
        fOutputList->Add(fDspT);
        
        fLambdaCpT = new TH1F("fLambdaCpT","Prompt Lammda_c pT;p_{T} (GeV/c);counts",500,0,50);
        fLambdaCpT->Sumw2();
        fOutputList->Add(fLambdaCpT);
    }
    
    if(fFillMCTemplates){
        fDElecDCA = new TH2F("fDElecDCA","D meson -> electron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., fnBinsDCAHisto,-0.4,0.4);
        fOutputList->Add(fDElecDCA);
        
        fBElecDCA = new TH2F("fBElecDCA","B meson -> electron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., fnBinsDCAHisto,-0.4,0.4);
        fOutputList->Add(fBElecDCA);
        
        fBHadElecDCA = new TH2F("fBHadElecDCA","B hadron -> electron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., fnBinsDCAHisto,-0.4,0.4);
        fOutputList->Add(fBHadElecDCA);
        
        fBMesonElecDCA = new TH2F("fBMesonElecDCA","B meson -> electron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., fnBinsDCAHisto,-0.4,0.4);
        fOutputList->Add(fBMesonElecDCA);
        
        fBBaryonElecDCA = new TH2F("fBBaryonElecDCA","B baryon -> electron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., fnBinsDCAHisto,-0.4,0.4);
        fOutputList->Add(fBBaryonElecDCA);
        
        fDHadElecDCA = new TH2F("fDHadElecDCA","D hadron -> electron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., fnBinsDCAHisto,-0.4,0.4);
        fOutputList->Add(fDHadElecDCA);
        
        fDMesonElecDCA = new TH2F("fDMesonElecDCA","D meson -> electron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., fnBinsDCAHisto,-0.4,0.4);
        fOutputList->Add(fDMesonElecDCA);
        
        fDBaryonElecDCA  = new TH2F("fDBaryonElecDCA","D baryon -> electron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., fnBinsDCAHisto,-0.4,0.4);
        fOutputList->Add(fDBaryonElecDCA);
        
        fLambdaCElecDCA = new TH2F("fLambdaCElecDCA","Lambda_c -> electron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., fnBinsDCAHisto,-0.4,0.4);
        fOutputList->Add(fLambdaCElecDCA);
        
        fD0ElecDCA = new TH2F("fD0ElecDCA","D0 -> electron DCA; p_{T}(GeV/c); DCAxMagFieldxSign; counts;", 60,0,30., fnBinsDCAHisto,-0.4,0.4);
        fOutputList->Add(fD0ElecDCA);
        
        Int_t binTemp[3] = {60,fnBinsDCAHisto,19}; //pT, DCA, Mom PID, Mom Gen, mompT
        Double_t xminTemp[3] = {0.,-0.4,0.5};
        Double_t xmaxTemp[3] = {30.,0.4,19.5};

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
        
    }
    PostData(1,fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskHFEBESpectraEMC::UserExec(Option_t *)
{
    // Main loop
    // Called for each event
    // Post output data.
    
    UInt_t evSelMask=((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
    
    fVevent = dynamic_cast<AliVEvent*>(InputEvent());
    if (!fVevent) {
        printf("ERROR: fVEvent not available\n");
        return;
    }
    
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    if (fESD) {
        //   printf("fESD available\n");
        //return;
    }
    
    //////////////
    //if Tender //
    //////////////
    if(fUseTender){
        //new branches with calibrated tracks and clusters
        if(IsAODanalysis()){ fTracks_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("tracks"));
        fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("caloClusters"));
        }
    }
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    fMCArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    fMCHeader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    
    ///////////////////
    //PID initialised//
    ///////////////////
    fpidResponse = fInputHandler->GetPIDResponse();
    
    ///////////////////
    // centrality
    /////////////////////
    
    Double_t centrality = -1;
    AliCentrality *fCentrality = (AliCentrality*)fAOD->GetCentrality();
    //centrality = fCentrality->GetCentralityPercentile("V0M");
    
    //Double_t centrality = -1;
    if(fAOD)fMultSelection = (AliMultSelection * ) fAOD->FindListObject("MultSelection");
    if( !fMultSelection) {
        //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
        //AliWarning("AliMultSelection object not found!");
        centrality = fCentrality->GetCentralityPercentile(fCentralityEstimator.Data());
    }else{
        //lPercentile = fMultSelection->GetMultiplicityPercentile("V0M");
        centrality = fMultSelection->GetMultiplicityPercentile(fCentralityEstimator.Data(), false);
    }
    
    if(fcentMim>-0.5)
    {
        if(centrality < fcentMim || centrality > fcentMax)return;
    }
    
    ////////////////
    //Event vertex//
    ////////////////
    Int_t ntracks = -999;
    if(!fUseTender)ntracks = fVevent->GetNumberOfTracks();
    if(fUseTender) ntracks = fTracks_tender->GetEntries();
    //if(ntracks < 1) printf("There are %d tracks in this event\n",ntracks);
    fMult->Fill(centrality,ntracks);
    
    fNevents->Fill(0); //all events
    Double_t Zvertex = -100, Xvertex = -100, Yvertex = -100;
    const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
    Double_t NcontV = pVtx->GetNContributors();
    if(NcontV<2)return;
    fNevents->Fill(1); //events with 2 tracks
    
    Zvertex = pVtx->GetZ();
    Yvertex = pVtx->GetY();
    Xvertex = pVtx->GetX();
    fVtxZ->Fill(Zvertex);
    fVtxX->Fill(Xvertex);
    fVtxY->Fill(Yvertex);
    
    ////////////////////
    //event selection///
    ////////////////////
    if(TMath::Abs(Zvertex)>10.0)return;
    fNevents->Fill(2); //events after z vtx cut
    fCent->Fill(centrality); //centrality dist.
    
    if(fMCHeader){
        ////////////////////////////////
        //Get number of Gen particles //
        ////////////////////////////////
        GetNMCPartProduced();
        
        /////////////////////////////////
        //Calculate Pi0 and Eta weight //
        /////////////////////////////////
        if(fCalculateWeight) GetPi0EtaWeight(fSprsPi0EtaWeightCal);
        
        /////////////////////////
        //Electrons in MC stack//
        /////////////////////////
        if(fCalculateElecRecoEffi) GetElectronFromStack();
        
        /////////////////////////////////
        //Histos for MC template Weight//
        /////////////////////////////////
        if(fCalculateMCTemplWeightCalc) GetMCTemplateWeight();
    }
    
    ////////////////
    // Mag. field //
    ////////////////
    Int_t fMagSign = 1;
    if(fAOD->GetMagneticField()<0) fMagSign = -1;
    
    /////////////////////////////
    //EMCAL cluster information//
    /////////////////////////////
    GetEMCalClusterInfo();
    
    ////////////////////////////////
    //Look for kink mother for AOD//
    ////////////////////////////////
    Int_t numberofvertices = 100;
    if(fAOD) numberofvertices = fAOD->GetNumberOfVertices();
    Double_t listofmotherkink[numberofvertices];
    Int_t numberofmotherkink = 0;
    if(IsAODanalysis())
    {
        for(Int_t ivertex=0; ivertex < numberofvertices; ivertex++) {
            AliAODVertex *aodvertex = fAOD->GetVertex(ivertex);
            if(!aodvertex) continue;
            if(aodvertex->GetType()==AliAODVertex::kKink) {
                AliAODTrack *mother = (AliAODTrack *) aodvertex->GetParent();
                if(!mother) continue;
                Int_t idmother = mother->GetID();
                listofmotherkink[numberofmotherkink] = idmother;
                numberofmotherkink++;
            }
        }
    } //+++
    
    ///////////////
    //Track loop///
    ///////////////
    
    fNEle = 0;
    
    for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
        
        AliVParticle* Vtrack = 0x0;
        if(!fUseTender) Vtrack  = fVevent->GetTrack(iTracks);
        if(fUseTender) Vtrack = dynamic_cast<AliVTrack*>(fTracks_tender->At(iTracks));
        
        if (!Vtrack) {
            printf("ERROR: Could not receive track %d\n", iTracks);
            continue;
        }
        AliVTrack *track = dynamic_cast<AliVTrack*>(Vtrack);
        AliESDtrack *etrack = dynamic_cast<AliESDtrack*>(Vtrack);
        AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(Vtrack);
        
        ///////////////////////
        // Get MC information//
        ///////////////////////
        Int_t pdg = -999;
        Int_t pidM = -1;
        Double_t pid_ele = 0.0;
        Bool_t IsMCEle = kFALSE, IsMCHFEle = kFALSE, IsMCDEle = kFALSE, IsMCBEle = kFALSE;

        if(fMCHeader && fCalculateElecRecoEffi){
            GetTrackHFStatus(track, IsMCEle, IsMCHFEle, IsMCBEle, IsMCDEle);
        }
        
        ////////////////////
        //Apply track cuts//
        ////////////////////
        if(fAOD)
                if(!atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue; //mimimum cuts
        
        //reject kink
        if(IsAODanalysis()){
            Bool_t kinkmotherpass = kTRUE;
            for(Int_t kinkmother = 0; kinkmother < numberofmotherkink; kinkmother++) {
                if(track->GetID() == listofmotherkink[kinkmother]) {
                    kinkmotherpass = kFALSE;
                    continue;
                }
            }
            if(!kinkmotherpass) continue;
        }
        else{
            if(etrack->GetKinkIndex(0) != 0) continue;
        }
        
        //other cuts
        Double_t d0z0[2]={-999,-999}, cov[3]={999,999,999};
        Double_t DCAxyCut = 0.25, DCAzCut = 1;
        
        if(atrack->GetTPCNcls() < 80) continue;
        if(atrack->GetITSNcls() < 3) continue;
        if((!(atrack->GetStatus()&AliESDtrack::kITSrefit)|| (!(atrack->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
        if(!(atrack->HasPointOnITSLayer(0) || atrack->HasPointOnITSLayer(1))) continue;
        
        if(fRecalIP) RecalImpactParam(atrack, d0z0, cov);
        else atrack->PropagateToDCA(pVtx, fVevent->GetMagneticField(), 20., d0z0, cov);
        
        if(TMath::Abs(d0z0[0]) > DCAxyCut || TMath::Abs(d0z0[1]) > DCAzCut) continue;
        
        fTrkDCA = -999.0;
        fTrkDCA = d0z0[0] * atrack->Charge() * fMagSign;
        
        ////////////////////
        //Track properties//
        ///////////////////
        Double_t dEdx =-999, fTPCnSigma_Pi=-999;
        Double_t TrkPhi=-999, TrkPt=-999, TrkEta=-999, TrkP = -999;
        dEdx = track->GetTPCsignal();
        fTPCnSigma = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
        fTPCnSigma_Pi = fpidResponse->NumberOfSigmasTPC(track, AliPID::kPion);
        TrkPhi = track->Phi();
        TrkPt = track->Pt();
        TrkEta = track->Eta();
        TrkP = track->P();
        
        if(track->GetID()<0) fNegTrkIDPt->Fill(track->Pt());
        fTrkPt->Fill(TrkPt);
        fTrketa->Fill(TrkEta);
        fTrkphi->Fill(TrkPhi);
        fdEdx->Fill(TrkP,dEdx);
        fTPCNpts->Fill(TrkP,track->GetTPCsignalN());
        fTPCnsig->Fill(TrkP,fTPCnSigma);
        fTPCnsig_Pi->Fill(TrkP,fTPCnSigma_Pi);
        
        if(TMath::Abs(track->Eta()) > 0.6) continue;
        if(TrkPt < 1) continue;
        
        /////////////////////////////
        //Reconstruction efficiency//
        /////////////////////////////
        if(fCalculateElecRecoEffi){
            if(IsMCEle) fInclElePhysPriTrkCuts->Fill(TrkPt);
            if(IsMCHFEle) fHFEPhysPriTrkCuts->Fill(TrkPt);
            if(IsMCBEle) fBEPhysPriTrkCuts->Fill(TrkPt);
            if(IsMCDEle) fDEPhysPriTrkCuts->Fill(TrkPt);
        }
        
        Bool_t fFillTem = kFALSE;
        if(fFillMCTemplates)
        {
           fFillTem = GetMCDCATemplates(track, fTrkDCA);
        }
        
        ///////////////////////////
        //Track matching to EMCAL//
        //////////////////////////
        if(!track->IsEMCAL()) continue;
        Int_t EMCalIndex = -1;
        EMCalIndex = track->GetEMCALcluster();
        if(EMCalIndex < 0) continue;
        fHistPtMatch->Fill(track->Pt());
        
        AliVCluster *clustMatch=0x0;
        if(!fUseTender) clustMatch = (AliVCluster*)fVevent->GetCaloCluster(EMCalIndex);
        if(fUseTender) clustMatch = dynamic_cast<AliVCluster*>(fCaloClusters_tender->At(EMCalIndex));
        
        Short_t NcellsInCluster = clustMatch->GetNCells();
        int icell=-1, iSM = -1;
        for(int icl=0; icl<NcellsInCluster; icl++)
        {
            if(icl==0) icell = clustMatch->GetCellAbsId(icl); //first cell = seed cell?
            if(fEMCALGeo)
                iSM = fEMCALGeo->GetSuperModuleNumber(icell);
        }
        
        Double_t emcphi = -999, emceta=-999;
        Bool_t fClsTypeEMC = kFALSE, fClsTypeDCAL = kFALSE;
        if(clustMatch && clustMatch->IsEMCAL())
        {
            // fEMCTrkMatch->Fill(clustMatch->GetTrackDx(),clustMatch->GetTrackDz());
            //if(TMath::Abs(clustMatch->GetTrackDx())>0.05 || TMath::Abs(clustMatch->GetTrackDz())>0.05) continue;
            
            Double_t fPhiDiff = -999, fEtaDiff = -999;
            GetTrkClsEtaPhiDiff(track, clustMatch, fPhiDiff, fEtaDiff);
            fEMCTrkMatch->Fill(fPhiDiff,fEtaDiff);
            fEMCTrkMatch_Phi->Fill(track->Pt(),fPhiDiff);
            fEMCTrkMatch_Eta->Fill(track->Pt(),fEtaDiff);
            
            if(TMath::Abs(fPhiDiff) > fDeltaPhi || TMath::Abs(fEtaDiff)> fDeltaEta) continue;
            
            /////////////////////////////////
            //Select EMCAL or DCAL clusters//
            /////////////////////////////////
            Float_t  emcx[3]; // cluster pos
            clustMatch->GetPosition(emcx);
            TVector3 clustpos(emcx[0],emcx[1],emcx[2]);
            emcphi = clustpos.Phi();
            emceta = clustpos.Eta();
            if(emcphi < 0) emcphi = emcphi+(2*TMath::Pi()); //TLorentz vector is defined between -pi to pi, so negative phi has to be flipped.
            if(emcphi > 1.39 && emcphi < 3.265) fClsTypeEMC = kTRUE; //EMCAL : 80 < phi < 187
            if(emcphi > 4.53 && emcphi < 5.708) fClsTypeDCAL = kTRUE;//DCAL  : 260 < phi < 327
            
            //----selects EMCAL+DCAL clusters when fFlagClsTypeEMC and fFlagClsTypeDCAL is kTRUE
            if(fFlagClsTypeEMC && !fFlagClsTypeDCAL)
                if(!fClsTypeEMC) continue; //selecting only EMCAL clusters
            
            if(fFlagClsTypeDCAL && !fFlagClsTypeEMC)
                if(!fClsTypeDCAL) continue; //selecting only DCAL clusters
            
            
            /////////////////////////////
            //Reconstruction efficiency//
            /////////////////////////////
            if(fCalculateElecRecoEffi){
                if(IsMCEle) fInclElePhysPriEMCMatch->Fill(TrkPt);
                if(IsMCHFEle) fHFEPhysPriEMCMatch->Fill(TrkPt);
                if(IsMCBEle) fBEPhysPriEMCMatch->Fill(TrkPt);
                if(IsMCDEle) fDEPhysPriEMCMatch->Fill(TrkPt);
            }
            
            /////////////////////////////////////////////
            //Properties of tracks matched to the EMCAL//
            /////////////////////////////////////////////
            fEMCTrkPt->Fill(TrkPt);
            if(TrkPt>1.0)
            {
                fEMCTrketa->Fill(TrkEta);
                fEMCTrkphi->Fill(TrkPhi);
            }
            fEMCdEdx->Fill(TrkP,dEdx);
            fEMCTPCnsig->Fill(TrkP,fTPCnSigma);
            fEMCTPCNpts->Fill(TrkP,track->GetTPCsignalN());
            
            Double_t clustMatchE = clustMatch->E();
            
            fClsEAftMatch->Fill(clustMatchE);
            fNonLinClsEAftMatch->Fill(clustMatch->GetNonLinCorrEnergy());
            
            fClsEtaPhiAftMatch->Fill(emceta,emcphi);
            
            if(TrkPhi > 1.396  && TrkPhi < 3.141) //emc acceptance (80 to 180 degrees)
                fClsEtaPhiAftMatchEMCin->Fill(emceta,emcphi);
            else
                fClsEtaPhiAftMatchEMCout->Fill(emceta,emcphi);
            
            //EMCAL EID info
            Double_t eop = -1.0, eop_NL = -1.0;
            Double_t m02 = -99999,m20 = -99999,sqm02m20=-99999.0;
            if(track->P()>0){
                eop = clustMatchE/track->P();
                eop_NL = clustMatch->GetNonLinCorrEnergy()/track->P();
            }
            m02 =clustMatch->GetM02();
            m20 =clustMatch->GetM20();
            
            if(track->Pt()>3.0){
                fHistdEdxEop->Fill(eop_NL,dEdx);
                fHistNsigEop->Fill(eop_NL,fTPCnSigma);
                fM20EovP->Fill(eop_NL,clustMatch->GetM20());
                fM02EovP->Fill(eop_NL,clustMatch->GetM02());
            }
            fM20->Fill(track->Pt(),clustMatch->GetM20());
            fM02->Fill(track->Pt(),clustMatch->GetM02());
            
            //EID THnsparse
            fvalueElectron[0] = track->Pt();
            fvalueElectron[1] = fTPCnSigma;
            fvalueElectron[2] = eop;
            fvalueElectron[3] = m20;
            fvalueElectron[4] = m02;
            fvalueElectron[5] = iSM;
            fvalueElectron[6] = eop_NL;
            
            if(fFlagSparse && track->Pt()>2.0){
                fSparseElectron->Fill(fvalueElectron);
            }
            
            /////////////////////////////
            //Reconstruction efficiency//
            /////////////////////////////
            if(fCalculateElecRecoEffi){
                GetEIDRecoEffi(track, clustMatch, IsMCEle, IsMCHFEle, IsMCBEle, IsMCDEle);
            }
            
            //////////////////
            //Apply EID cuts//
            //////////////////
            
            Bool_t fHadTrack = kFALSE, fElectTrack = kFALSE;
            fElectTrack = PassEIDCuts(track, clustMatch, fHadTrack);
            
            if(fHadTrack){
                fHadPt_AftEID->Fill(TrkPt);
                fHadDCA->Fill(TrkPt,fTrkDCA);
            }
            
            if(!fElectTrack) continue;
            
            fInclsElecPt->Fill(TrkPt);
            fInclElecDCA->Fill(TrkPt,fTrkDCA);
            
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
            Bool_t fFlagNonHFE=kFALSE;
            SelectPhotonicElectron(iTracks,track,fFlagNonHFE,pidM);
            
            //////////////////////////////////
            //Non-HFE efficiency calculation//
            //////////////////////////////////
            if(fMCHeader && fCalculateNonHFEEffi){
                if(fFlagNonHFE){
                    EffiNumTag = GetNonHFEEffiRecoTag(track);
                }
            }

        }
    } //track loop
    
    fNElecInEvt->Fill(fNEle);
    
    PostData(1, fOutputList);
}
//___________________________________________
Bool_t AliAnalysisTaskHFEBESpectraEMC::PassEIDCuts(AliVTrack *track, AliVCluster *clust, Bool_t &Hadtrack)
{
    //apply electron identification cuts
    
    Bool_t hadTrk = kFALSE;
    Double_t eop = -1.0, eop_NL = -1.0;
    Double_t m02 = -999,m20 = -999;
    Double_t clustE = clust->E();
    Double_t clustE_NL = clust->GetNonLinCorrEnergy();
    Double_t TrkPt = track->Pt();
    if(track->P()>0){
        eop = clustE/track->P();
        eop_NL = clustE_NL/track->P();
    }
    m02 =clust->GetM02();
    m20 =clust->GetM20();
    
    //Hadron E/p distribution
    if(fTPCnSigma > fTPCnSigmaHadMin && fTPCnSigma < fTPCnSigmaHadMax)
    {
        if(TrkPt < 8.0){
            if((m02 > fM02Min && m02 < fM02Max1) && (m20 > fM20Min && m20 < fM20Max))
                {
                    fHadEovp_AftEID->Fill(TrkPt,eop);
                    fHadEovpNL_AftEID->Fill(TrkPt,eop_NL);

                    if(eop_NL > fEovPMin && eop_NL < fEovPMax) hadTrk=kTRUE;
                }
        }
        if(TrkPt >= 8.0){
            if((m02 > fM02Min && m02 < fM02Max2) && (m20 > fM20Min && m20 < fM20Max))
            {
                fHadEovp_AftEID->Fill(TrkPt,eop);
                fHadEovpNL_AftEID->Fill(TrkPt,eop_NL);
                if(eop_NL > fEovPMin && eop_NL < fEovPMax) hadTrk=kTRUE;
            }
        }
    }
    Hadtrack = hadTrk;
    
    if(fTPCnSigma < fTPCnSigmaMin || fTPCnSigma > fTPCnSigmaMax) return kFALSE;
    if(TrkPt < 8.0){
        if(m02 < fM02Min || m02 > fM02Max1) return kFALSE;
    }
    if(TrkPt >= 8.0){
        if(m02 < fM02Min || m02 > fM02Max2) return kFALSE;
    }
    if(m20 < fM20Min || m20 > fM20Max) return kFALSE;
    
    fEop_AftEID->Fill(TrkPt,eop);
    fEopNL_AftEID->Fill(TrkPt,eop_NL);
    
    if(eop_NL < fEovPMin || eop_NL > fEovPMax) return kFALSE;
    
    return kTRUE;
}
//________________________________________________________________________
void AliAnalysisTaskHFEBESpectraEMC::SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec, Int_t iMC)
{
    ///////////////////////////////////////////
    //////Non-HFE - Invariant mass method//////
    ///////////////////////////////////////////
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
    Double_t d0z0[2]={-999,-999}, cov[3];
    Double_t DCAxyCut = 0.25, DCAzCut = 1;
    
    Bool_t flagPhotonicElec = kFALSE, flagLSElec = kFALSE;
    Double_t ptAsso=-999., nsigma=-999.0;
    Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
    
    Int_t ntracks = -999;
    if(!fUseTender)ntracks = fVevent->GetNumberOfTracks();
    if(fUseTender) ntracks = fTracks_tender->GetEntries();
    
    for (Int_t jtrack = 0; jtrack < ntracks; jtrack++) {
        AliVParticle* VAssotrack = 0x0;
        if(!fUseTender) VAssotrack  = fVevent->GetTrack(jtrack);
        if(fUseTender) VAssotrack = dynamic_cast<AliVTrack*>(fTracks_tender->At(jtrack)); //take tracks from Tender list
        
        if (!VAssotrack) {
            printf("ERROR: Could not receive track %d\n", jtrack);
            continue;
        }
        
        AliVTrack *Assotrack = dynamic_cast<AliVTrack*>(VAssotrack);
        AliESDtrack *eAssotrack = dynamic_cast<AliESDtrack*>(VAssotrack);
        AliAODTrack *aAssotrack = dynamic_cast<AliAODTrack*>(VAssotrack);
        
        //------reject same track
        if(jtrack==itrack) continue;
        
        Double_t mass=-999., width = -999;
        Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;
        
        nsigma = fpidResponse->NumberOfSigmasTPC(Assotrack, AliPID::kElectron);
        ptAsso = Assotrack->Pt();
        
        //------track cuts applied
        if(fAOD) {
            if(!aAssotrack->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
            if(aAssotrack->GetTPCNcls() < 70) continue;
            if((!(aAssotrack->GetStatus()&AliESDtrack::kITSrefit)|| (!(aAssotrack->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
            
            if(fRecalIP) RecalImpactParam(aAssotrack, d0z0, cov);
            else aAssotrack->PropagateToDCA(pVtx, fVevent->GetMagneticField(), 20., d0z0, cov);
            
            if(TMath::Abs(d0z0[0]) > DCAxyCut || TMath::Abs(d0z0[1]) > DCAzCut) continue;
        }
        
        //-------loose cut on partner electron
        if(ptAsso <0.150) continue;
        if(aAssotrack->Eta()<-0.9 || aAssotrack->Eta()>0.9) continue;
        if(nsigma < -3 || nsigma > 3) continue;
        
        Int_t chargeAsso = Assotrack->Charge();
        Int_t charge = track->Charge();
        if(charge>0) fPDGe1 = -11;
        if(chargeAsso>0) fPDGe2 = -11;
        
        fFlagLS=kFALSE; fFlagULS=kFALSE;
        if(charge == chargeAsso) fFlagLS = kTRUE;
        if(charge != chargeAsso) fFlagULS = kTRUE;
        
        //-------define KFParticle to get mass
        AliKFParticle::SetField(fVevent->GetMagneticField());
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
        if(fMCHeader && fCalculateNonHFEEffi){
            EffiNumULSLS = GetNonHFEEffiULSLS(track, Assotrack, fFlagLS, fFlagULS, mass);
        }

        Double_t TrkPt = track->Pt();
        if(mass < fInvmassCut){
            if(fFlagLS){
                fLSElecPt->Fill(TrkPt);
                fLSElecDCA->Fill(TrkPt,fTrkDCA);
            }

            if(fFlagULS){
                fULSElecPt->Fill(TrkPt);
                fULSElecDCA->Fill(TrkPt,fTrkDCA);
            }
        }
        
        if(mass < fInvmassCut && fFlagULS && !flagPhotonicElec)
            flagPhotonicElec = kTRUE; //Tag Non-HFE (random mass cut, not optimised)
    }
    fFlagPhotonicElec = flagPhotonicElec;
}

//________________________________________________________________________
void AliAnalysisTaskHFEBESpectraEMC::GetTrkClsEtaPhiDiff(AliVTrack *t, AliVCluster *v, Double_t &phidiff, Double_t &etadiff)
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
//________________________________________________________________________
void AliAnalysisTaskHFEBESpectraEMC::GetEMCalClusterInfo()
{
    //Get basic EMCal cluster information
    
    Int_t Nclust = -999;
    if(!fUseTender) Nclust = fVevent->GetNumberOfCaloClusters();
    if(fUseTender) Nclust = fCaloClusters_tender->GetEntries();
    
    int NclustAll= 0;
    
    Bool_t fClsTypeEMC = kFALSE, fClsTypeDCAL = kFALSE;
    
    for(Int_t icl=0; icl<Nclust; icl++)
    {
        AliVCluster *clust = 0x0;
        if(!fUseTender) clust = fVevent->GetCaloCluster(icl);
        if(fUseTender) clust = dynamic_cast<AliVCluster*>(fCaloClusters_tender->At(icl));
        if(!clust)  printf("ERROR: Could not receive cluster matched calibrated from track %d\n", icl);
        
        fClsTypeEMC = kFALSE; fClsTypeDCAL = kFALSE;
        
        if(clust && clust->IsEMCAL())
        {
            Double_t clustE_NL = clust->GetNonLinCorrEnergy();
            if(clustE_NL < 0.3) continue;
            
            /////////////////////////////////
            //Select EMCAL or DCAL clusters//
            /////////////////////////////////
            Float_t  emcx[3]; // cluster pos
            clust->GetPosition(emcx);
            TVector3 clustpos(emcx[0],emcx[1],emcx[2]);
            Double_t emcphi = clustpos.Phi();
            Double_t emceta = clustpos.Eta();
            if(emcphi < 0) emcphi = emcphi+(2*TMath::Pi()); //TLorentz vector is defined between -pi to pi, so negative phi has to be flipped.
            
            if(emcphi > 1.39 && emcphi < 3.265) fClsTypeEMC = kTRUE; //EMCAL : 80 < phi < 187
            if(emcphi > 4.53 && emcphi < 5.708) fClsTypeDCAL = kTRUE; //DCAL  : 260 < phi < 327
            
            //----selects EMCAL+DCAL clusters when fFlagClsTypeEMC and fFlagClsTypeDCAL is kTRUE
            if(fFlagClsTypeEMC && !fFlagClsTypeDCAL)
                if(!fClsTypeEMC) continue; //selecting only EMCAL clusters
            
            if(fFlagClsTypeDCAL && !fFlagClsTypeEMC)
                if(!fClsTypeDCAL) continue; //selecting only DCAL clusters
            
            fHistClustE->Fill(clust->E());
            fHistNonLinClustE->Fill(clustE_NL);
            
           // if(centrality>-1)fHistClustEcent->Fill(centrality,clustE_NL);
            fEMCClsEtaPhi->Fill(emceta,emcphi);
            fHistoNCells->Fill(clustE_NL,clust->GetNCells());
            Double_t EperCell = -999.9;
            if(clust->GetNCells()>0)EperCell = clustE_NL/clust->GetNCells();
            fHistoEperCell->Fill(clustE_NL,EperCell);
            
            Float_t tof = clust->GetTOF()*1e+9; // ns
            fHistoTimeEMC->Fill(clustE_NL,tof);
            
            NclustAll++;
        }
    }
    fHistoNCls->Fill(NclustAll);
    
    // cell information
    AliVCaloCells *fCaloCells = fVevent->GetEMCALCells();
    
    //Int_t nSACell, iSACell, mclabel;
    Short_t cellAddr, nSACell;
    Int_t  mclabel;
    Short_t iSACell;
    Double_t cellAmp=-1., cellTimeT=-1., clusterTime=-1., efrac=-1.;
    
    nSACell = fCaloCells->GetNumberOfCells();
    for(iSACell = 0; iSACell < nSACell; iSACell++ ){
        Bool_t haveCell = fCaloCells->GetCell(iSACell, cellAddr, cellAmp, cellTimeT , mclabel, efrac);
        //virtual Bool_t   GetCell(Short_t pos, Short_t &cellNumber, Double_t &amplitude, Double_t &time, Int_t &mclabel,    Double_t  &efrac)
        if(haveCell)fHistoCalCell->Fill(cellAddr,cellAmp);
        
    }
    
    if(!fEMCALGeo)fEMCALGeo  = AliEMCALGeometry::GetInstance(); // not work w.o. Tender
}
//________________________________________________________________________
void AliAnalysisTaskHFEBESpectraEMC::FindMother(AliAODMCParticle* part, Int_t &label, Int_t &pid)
{
    // Find mother in case of MC
    
    if(part->GetMother()>-1)
    {
        label = part->GetMother();
        AliAODMCParticle *partM = (AliAODMCParticle*)fMCArray->At(label);
        pid = partM->GetPdgCode();
    }
    else
    {
        pid = -1;
    }
}
//_________________________________________
Bool_t AliAnalysisTaskHFEBESpectraEMC::GetNMCPartProduced()
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
//_________________________________________
void AliAnalysisTaskHFEBESpectraEMC::GetPi0EtaWeight(THnSparse *SparseWeight)
{
    //Get pi0 and eta information for weight calculation
    
    Double_t fvalue[4] = {-999,-999,-999,-999};
    
    for(int imc=0; imc< fNTotMCpart; imc++)
    {
        AliAODMCParticle *AODMCtrack = (AliAODMCParticle*)fMCArray->At(imc);
        if(TMath::Abs(AODMCtrack->Eta()) > 0.9) continue;
        
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
        
        SparseWeight->Fill(fvalue);
    }
}
//_________________________________________
Int_t AliAnalysisTaskHFEBESpectraEMC::GetPi0EtaType(AliAODMCParticle *part)
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
//_________________________________________
Bool_t AliAnalysisTaskHFEBESpectraEMC::GetNonHFEEffiDenom(AliVTrack *track)
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
    
    Bool_t fNonHFE = IsNonHFE(MCPart, fFromMB, ftype, iMCmom, MomPDG, MomPt);
    if(!fNonHFE) return kFALSE;
    fNonHFeTrkPt->Fill(TrkPt);
    
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
    
    if(fIsFrmEmbPi0 || fIsFrmEmbEta){
        fNonHFeEmbTrkPt->Fill(TrkPt);
        
        if(fIsFrmEmbPi0) {
            fWeight = fWeightPi0;
            fPi0eEmbWeightTrkPt->Fill(TrkPt,fWeightPi0);
            fNonHFeEmbWeightTrkPt->Fill(TrkPt,fWeightPi0);
        }
        if(fIsFrmEmbEta){
            fWeight = fWeightEta;
            fEtaeEmbWeightTrkPt->Fill(TrkPt,fWeightEta);
            fNonHFeEmbWeightTrkPt->Fill(TrkPt,fWeightEta);
        }
    }
    
    return kTRUE;
}

//___________________________________________
Bool_t  AliAnalysisTaskHFEBESpectraEMC::IsNonHFE(AliAODMCParticle *MCPart, Bool_t &fFromMB, Int_t &type, Int_t &iMCmom, Int_t &MomPDG, Double_t &MomPt)
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
//_________________________________________
Bool_t AliAnalysisTaskHFEBESpectraEMC::GetNonHFEEffiRecoTag(AliVTrack *track)
{
    Double_t TrkPt = track->Pt();
    
    fRecoNonHFeTrkPt->Fill(TrkPt);
    if(fIsFrmEmbPi0 || fIsFrmEmbEta){
        fRecoNonHFeEmbTrkPt->Fill(TrkPt);
        
        if(fIsFrmEmbPi0) {
            fRecoPi0eEmbWeightTrkPt->Fill(TrkPt,fWeightPi0);
            fRecoNonHFeEmbWeightTrkPt->Fill(TrkPt,fWeightPi0);
        }
        if(fIsFrmEmbEta){
            fRecoEtaeEmbWeightTrkPt->Fill(TrkPt,fWeightEta);
            fRecoNonHFeEmbWeightTrkPt->Fill(TrkPt,fWeightEta);
        }
    }
    
    return kTRUE;
}
//_________________________________________
Bool_t AliAnalysisTaskHFEBESpectraEMC::GetNonHFEEffiULSLS(AliVTrack *track, AliVTrack *Assotrack, Bool_t fFlagLS, Bool_t fFlagULS, Double_t mass)
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
void AliAnalysisTaskHFEBESpectraEMC::GetElectronFromStack()
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
        
        Bool_t IsMCEle = kFALSE, IsMCHFEle = kFALSE, IsMCDEle = kFALSE, IsMCBEle = kFALSE;

        if(TMath::Abs(MCPart->Eta()) > 0.6) continue;
        if(!MCPart->IsPhysicalPrimary()) continue;
        if(!(PDGcode == 11)) continue;
        
        IsMCEle = kTRUE;
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
//______________________________________
void AliAnalysisTaskHFEBESpectraEMC::GetTrackHFStatus(AliVTrack *track, Bool_t &IsMCEle, Bool_t &IsMCHFEle, Bool_t &IsMCBEle, Bool_t &IsMCDEle)
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
        if(MCPart->IsPhysicalPrimary()){
            if(TMath::Abs(MCPart->GetPdgCode())==11){
                IsMCEle = kTRUE;
                
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
//______________________________________
void AliAnalysisTaskHFEBESpectraEMC::GetEIDRecoEffi(AliVTrack *track, AliVCluster *clust, Bool_t IsMCEle, Bool_t IsMCHFEle, Bool_t IsMCBEle, Bool_t IsMCDEle)
{
    //Filling histograms for EID efficiency
    
    Bool_t PassSSCut = kFALSE;
    
    Double_t eop_NL = -1.0;
    Double_t m02 = -999,m20 = -999;
    Double_t clustE_NL = clust->GetNonLinCorrEnergy();
    Double_t TrkPt = track->Pt();
    if(track->P()>0)eop_NL = clustE_NL/track->P();
    m02 =clust->GetM02();
    m20 =clust->GetM20();
    
    if(eop_NL > fEovPMin && eop_NL < fEovPMax){
        if(IsMCEle) fInclElePhysPriEovP->Fill(TrkPt);
        if(IsMCHFEle) fHFEPhysPriEovP->Fill(TrkPt);
        if(IsMCBEle) fBEPhysPriEovP->Fill(TrkPt);
        if(IsMCDEle) fDEPhysPriEovP->Fill(TrkPt);
        
        if(fTPCnSigma > fTPCnSigmaMin && fTPCnSigma < fTPCnSigmaMax){
            if(IsMCEle) fInclElePhysPriTPCnsig->Fill(TrkPt);
            if(IsMCHFEle) fHFEPhysPriTPCnsig->Fill(TrkPt);
            if(IsMCBEle) fBEPhysPriTPCnsig->Fill(TrkPt);
            if(IsMCDEle) fDEPhysPriTPCnsig->Fill(TrkPt);
            
            if(eop_NL > fEovPMin && eop_NL < fEovPMax){
                if(IsMCEle) fInclElePhysPriEovPBfrSS->Fill(TrkPt);
                if(IsMCHFEle) fHFEPhysPriEovPBfrSS->Fill(TrkPt);
                if(IsMCBEle) fBEPhysPriEovPBfrSS->Fill(TrkPt);
                if(IsMCDEle) fDEPhysPriEovPBfrSS->Fill(TrkPt);
            }
            
            if(TrkPt < 8.0){
                if(m02 > fM02Min && m02 < fM02Max1) PassSSCut = kTRUE;
            }
            if(TrkPt >= 8.0){
                if(m02 > fM02Min && m02 < fM02Max2) PassSSCut = kTRUE;
            }
            if(m20 > fM20Min && m20 < fM20Max) PassSSCut = kTRUE;
            
            if(PassSSCut){
                if(IsMCEle) fInclElePhysPriSS->Fill(TrkPt);
                if(IsMCHFEle) fHFEPhysPriSS->Fill(TrkPt);
                if(IsMCBEle) fBEPhysPriSS->Fill(TrkPt);
                if(IsMCDEle) fDEPhysPriSS->Fill(TrkPt);
            }
        }
    }
}
//________________________________________________________________________
void AliAnalysisTaskHFEBESpectraEMC::GetMCTemplateWeight()
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
//________________________________________________________________________
Bool_t AliAnalysisTaskHFEBESpectraEMC::GetMCDCATemplates(AliVTrack *track, Double_t TrkDCA)
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
            
            if(MomPDGDummy>500 && MomPDGDummy<600){
                fBMesonElecDCA->Fill(TrkPt,TrkDCA);
                fpidSort = 1; //Mom is B
                if(fIsAnapp) GetBWeight(MCPartMomDummy, fWeightB, fWeightBMin, fWeightBMax);
                if(!fIsAnapp) GetBWeightPbPb(MCPartMomDummy, fWeightB);
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
        
        if(MomPDG>400 && MomPDG<500) {
            fDMesonElecDCA->Fill(TrkPt,TrkDCA);
            fpidSort = 2; //Mom is D
            if(fIsAnapp) GetDWeight(MCPartMom, fWeightD, fWeightDUp, fWeightDDown);
            if(!fIsAnapp) GetDWeightPbPb(MCPartMom, MomPDG, fWeightD);
        }
        
        if(MomPDG>4000 && MomPDG<5000) {
            fDBaryonElecDCA->Fill(TrkPt,TrkDCA);
            fpidSort = 9; //Mom is c Baryon
            if(!fIsAnapp)
            {
                if(MomPDG == 4122) GetDWeightPbPb(MCPartMom, MomPDG, fWeightD); //For Lc
            }
        }
        if(MomPDG == 411) fpidSort = 11; //Mom is D+
        if(MomPDG == 421) fpidSort = 12; //Mom is D0
        if(MomPDG == 413) fpidSort = 14; //Mom is D*+
        if(MomPDG == 431) fpidSort = 15; //Mom is Ds
        if(MomPDG > 431 && MomPDG < 436) fpidSort = 16; //Mom is other Ds
        if(MomPDG == 4122) fpidSort = 17; //Mom is Lambda c
        if(MomPDG == 443) fpidSort = 6; //Mom is J/Psi
        
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
        if(fIsAnapp){
            fSprsTemplatesWeightVar1->Fill(fvalue, fWeightBMin);
            fSprsTemplatesWeightVar2->Fill(fvalue, fWeightBMax);
        }
    }
    if(IsDEle) {
        fSprsTemplatesWeight->Fill(fvalue, fWeightD);
        if(fIsAnapp){
            fSprsTemplatesWeightVar1->Fill(fvalue, fWeightDUp);
            fSprsTemplatesWeightVar2->Fill(fvalue, fWeightDDown);
        }
    }
    
    return kTRUE;
}
//________________________________________________________________________
/*Bool_t AliAnalysisTaskHFEBESpectraEMC::GetMCDCATemplates(AliVTrack *track, Double_t TrkDCA)
{
    //Fill MC template histograms
    //Change this with beauty loop
    
    Int_t iTrklabel = TMath::Abs(track->GetLabel());
    Double_t  TrkPt = track->Pt();
    
    AliAODMCParticle *MCPart;
    AliAODMCParticle *MCPartMom;
    AliAODMCParticle *MCPartGMom;
    AliAODMCParticle *MCPartGGMom;
    
    Int_t iMCmom = -999, iMCgmom = -999, iMCggmom = -999, iMCgggmom = -999;
    Int_t MomPDG = -999, GMomPDG=-999, GGMomPDG=-999, GGGMomPDG=-999;
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
    
    if((MomPDG>400 && MomPDG<600) || (MomPDG>4000 && MomPDG<6000)){ //D,B ->e
        IsHFEle = kTRUE;
        
        if((MomPDG>500 && MomPDG<600) || (MomPDG>5000 && MomPDG<6000)){ //B->e
            IsBEle = kTRUE;
            fBHadElecDCA->Fill(TrkPt,TrkDCA);
            
            if(MomPDG>500 && MomPDG<600) {
                fBMesonElecDCA->Fill(TrkPt,TrkDCA);
                fpidSort = 1; //Mom is B
                if(fIsAnapp) GetBWeight(MCPartMom, fWeightB, fWeightBMin, fWeightBMax);
                if(!fIsAnapp) GetBWeightPbPb(MCPartMom, fWeightB);
            }
            if(MomPDG>5000 && MomPDG<6000){
                fBBaryonElecDCA->Fill(TrkPt,TrkDCA);
                fpidSort = 10; //Mom is b Baryon
            }
        }
        else{
            iMCgmom = MCPartMom->GetMother();
            if(iMCgmom > 0){
                MCPartGMom = (AliAODMCParticle*)fMCArray->At(iMCgmom);
                GMomPDG = TMath::Abs(MCPartGMom->GetPdgCode());
                
                if((GMomPDG>500 && GMomPDG<600) || (GMomPDG>5000 && GMomPDG<6000)){ //B->D->e
                    IsBEle = kTRUE;
                    fBHadElecDCA->Fill(TrkPt,TrkDCA);
                    
                    if(GMomPDG>500 && GMomPDG<600){
                        fBMesonElecDCA->Fill(TrkPt,TrkDCA);
                        fpidSort = 1; //Mom is B
                        if(fIsAnapp) GetBWeight(MCPartGMom, fWeightB, fWeightBMin, fWeightBMax);
                        if(!fIsAnapp) GetBWeightPbPb(MCPartGMom, fWeightB);
                    }
                    if(GMomPDG>5000 && GMomPDG<6000){
                        fBBaryonElecDCA->Fill(TrkPt,TrkDCA);
                        fpidSort = 10; //Mom is b Baryon
                    }
                }
                else{
                    iMCggmom = MCPartGMom->GetMother();
                    if(iMCggmom > 0){
                        MCPartGGMom = (AliAODMCParticle*)fMCArray->At(iMCggmom);
                        GGMomPDG = TMath::Abs(MCPartGGMom->GetPdgCode());
                        
                        if((GGMomPDG>500 && GGMomPDG<600) || (GGMomPDG>5000 && GGMomPDG<6000)){ //B->D->D->e
                            IsBEle = kTRUE;
                            fBHadElecDCA->Fill(TrkPt,TrkDCA);
                            
                            if(GGMomPDG>500 && GGMomPDG<600){
                                fBMesonElecDCA->Fill(TrkPt,TrkDCA);
                                fpidSort = 1; //Mom is B
                                if(fIsAnapp) GetBWeight(MCPartGGMom, fWeightB, fWeightBMin, fWeightBMax);
                                if(!fIsAnapp) GetBWeightPbPb(MCPartGGMom, fWeightB);
                            }
                            if(GGMomPDG>5000 && GGMomPDG<6000){
                                fBBaryonElecDCA->Fill(TrkPt,TrkDCA);
                                fpidSort = 10; //Mom is b Baryon
                            }
                        }
                    }
                }
            }
        }
        
        if(!IsBEle){
            IsDEle = kTRUE;
            fDHadElecDCA->Fill(TrkPt,TrkDCA);
            
            if(MomPDG>400 && MomPDG<500) {
                fDMesonElecDCA->Fill(TrkPt,TrkDCA);
                fpidSort = 2; //Mom is D
                if(fIsAnapp) GetDWeight(MCPartMom, fWeightD, fWeightDUp, fWeightDDown);
                if(!fIsAnapp) GetDWeightPbPb(MCPartMom, MomPDG, fWeightD);
            }
            if(MomPDG>4000 && MomPDG<5000) {
                fDBaryonElecDCA->Fill(TrkPt,TrkDCA);
                fpidSort = 9; //Mom is c Baryon
                if(!fIsAnapp)
                {
                    if(MomPDG == 4122) GetDWeightPbPb(MCPartMom, MomPDG, fWeightD); //For Lc
                }
            }
            if(MomPDG == 411) fpidSort = 11; //Mom is D+
            if(MomPDG == 421) fpidSort = 12; //Mom is D0
            if(MomPDG == 413) fpidSort = 14; //Mom is D*+
            if(MomPDG == 431) fpidSort = 15; //Mom is Ds
            if(MomPDG > 431 && MomPDG < 436) fpidSort = 16; //Mom is other Ds
            if(MomPDG == 4122) fpidSort = 17; //Mom is Lambda c
            if(MomPDG == 443) fpidSort = 6; //Mom is J/Psi
            
            if(MomPDG == 4122) fLambdaCElecDCA->Fill(TrkPt,TrkDCA);
            if(MomPDG == 421) fD0ElecDCA->Fill(TrkPt,TrkDCA);
        }
  
        fvalue[0] = TrkPt;
        fvalue[1] = TrkDCA;
        fvalue[2] = fpidSort;
        fSprsTemplatesNoWeight->Fill(fvalue);
        
        if(IsBEle) {
            fSprsTemplatesWeight->Fill(fvalue, fWeightB);
            if(fIsAnapp){
                fSprsTemplatesWeightVar1->Fill(fvalue, fWeightBMin);
                fSprsTemplatesWeightVar2->Fill(fvalue, fWeightBMax);
            }
        }
        if(IsDEle) {
            fSprsTemplatesWeight->Fill(fvalue, fWeightD);
            if(fIsAnapp){
                fSprsTemplatesWeightVar1->Fill(fvalue, fWeightDUp);
                fSprsTemplatesWeightVar2->Fill(fvalue, fWeightDDown);
            }
        }
    }
    return kTRUE;
}
 */
//________________________________________________________________________
void AliAnalysisTaskHFEBESpectraEMC::SetDmesonWeightHist(TH1 *D1, TH1 *D2, TH1 *D3)
{
    fDcent = (TH1F*)D1->Clone();
    fDUp = (TH1F*)D2->Clone();
    fDDown = (TH1F*)D3->Clone();
}
//________________________________________________________________________
void AliAnalysisTaskHFEBESpectraEMC::SetBmesonWeightHist(TH1 *B1, TH1 *B2, TH1 *B3)
{
    fBcent = (TH1F*)B1->Clone();
    fBMin = (TH1F*)B2->Clone();
    fBMax = (TH1F*)B3->Clone();
}
//________________________________________________________________________
void AliAnalysisTaskHFEBESpectraEMC::SetDmesonWeightHistPbPb(TH1 *D0, TH1 *DPlus, TH1 *Ds, TH1 *Lc)
{
    fD0 = (TH1F *)D0->Clone();
    fDPlus = (TH1F *)DPlus->Clone();
    fDs = (TH1F *)Ds->Clone();
    fLc = (TH1F *)Lc->Clone();
}
//________________________________________________________________________
void AliAnalysisTaskHFEBESpectraEMC::SetBmesonWeightHistPbPb(TH1 *B)
{
    fB = (TH1D *)B->Clone();
}
/*
//________________________________________________________________________
void AliAnalysisTaskHFEBESpectraEMC::InputWeightCorrectionMaps()
{
    //Get the input files for D and B meson pT weight
    
    TString DMesonWeightMaps, BMesonWeightMaps;
    
     DMesonWeightMaps = "alien:///alice/cern.ch/user/d/dthomas/DandBmesonpTweightCorrectionFiles/DMesonpTWeight.root";
     BMesonWeightMaps = "alien:///alice/cern.ch/user/d/dthomas/DandBmesonpTweightCorrectionFiles/BMesonpTWeight.root";
    
    //   printf("\n### reading file %s ...\n",DMesonWeightMaps.Data());
    //   printf("\n### reading file %s ...\n",BMesonWeightMaps.Data());
    
    TFile* f2 = TFile::Open(DMesonWeightMaps.Data());
    if(f2){
        fDcent = (TH1D*)f2->Get("RatD0");
        fDUp = (TH1D*)f2->Get("RatD0Up");
        fDDown = (TH1D*)f2->Get("RatD0Down");
    }
    f2->Close();
    TFile* f3 = TFile::Open(BMesonWeightMaps.Data());
    if(f3){
        fBcent = (TH1D*)f3->Get("RatBMes");
        fBMin = (TH1D*)f3->Get("RatBMesMin");
        fBMax = (TH1D*)f3->Get("RatBMesMax");
    }
    f3->Close();
}
 */
//________________________________________________________________________
void AliAnalysisTaskHFEBESpectraEMC::GetDWeightPbPb(AliAODMCParticle *Part, Int_t PDG, Double_t &DCentWeight)
{
    //D meson weight

    Int_t bin = -999;
    Int_t binLast = -999;
    
    if(!fD0){
        DCentWeight = 1.0;
        return;
    }
    
    bin = fD0->FindBin(Part->Pt());
    binLast = fD0->FindBin(35.9);
    
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
        if(PDG == 4122) DCentWeight = fLc->GetBinContent(binLast);
        return;
    }
    
    if(PDG == 421) DCentWeight = fD0->GetBinContent(bin);
    if(PDG == 411) DCentWeight = fDPlus->GetBinContent(bin);
    if(PDG == 431) DCentWeight = fDs->GetBinContent(bin);
    if(PDG == 4122) DCentWeight = fLc->GetBinContent(bin);
    
    return;
}
//________________________________________________________________________
void AliAnalysisTaskHFEBESpectraEMC::GetBWeightPbPb(AliAODMCParticle *Part, Double_t &BCentWeight)
{
    //B meson weight
    
    Int_t bin = -999;
    Int_t binLast = -999;
    
    if(!fB){
        BCentWeight = 1.0;
        return;
    }

    bin = fB->FindBin(Part->Pt());
    binLast = fB->FindBin(49.9);

    if(fB->IsBinUnderflow(bin)){
        BCentWeight = fB->GetBinContent(1);
        return;
    }
    if(Part->Pt() > 35.9){
        BCentWeight = fB->GetBinContent(binLast);
        return;
    }
    
    BCentWeight = fB->GetBinContent(bin);
    
    return;
}
//________________________________________________________________________
void AliAnalysisTaskHFEBESpectraEMC::GetBWeight(AliAODMCParticle *Part, Double_t &BCentWeight, Double_t &BMinWeight, Double_t &BMaxWeight)
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
void AliAnalysisTaskHFEBESpectraEMC::GetDWeight(AliAODMCParticle *Part, Double_t &DCentWeight, Double_t &DMinWeight, Double_t &DMaxWeight)
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
    binLast = fDcent->FindBin(49.9);
    
    if(fDcent->IsBinUnderflow(bin)){
        DCentWeight = 1.0;
        DMinWeight = 1.0;
        DMaxWeight = 1.0;
        return;
    }
    if(Part->Pt() > 49.9) {
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
void AliAnalysisTaskHFEBESpectraEMC::RecalImpactParam(const AliVTrack * const track, Double_t dcaD[2], Double_t covD[3])
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
//________________________________________________________________________
AliAODVertex* AliAnalysisTaskHFEBESpectraEMC::RemoveDaughtersFromPrimaryVtx(const AliVTrack * const track)
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
//________________________________________________________________________
void AliAnalysisTaskHFEBESpectraEMC::Terminate(Option_t *)
{
    // Draw result to the screen
    // Called once at the end of the query
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) {
        printf("ERROR: Output list not available\n");
        return;
    }
}
