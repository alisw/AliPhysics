
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

////////////////////////////////////////////////////////////////////////
//                                                                    //
//      Task for Heavy-flavour electron analysis in pPb collisions    //
//      (+ Electron-Hadron Jetlike Azimuthal Correlation)             //
//																	  //
//		version: Nov 18, 2016.	     						          //
//                                                                    //
//	    Authors 							                          //
//		Elienos Pereira de Oliveira Filho (epereira@cern.ch)	      //
//		Cristiane Jahnke		(cristiane.jahnke@cern.ch)            //
//      Henrique Zanoli (h.zanoli@cern.ch)                            //
//      Alexis Mas (aleximas@if.usp.br)                               //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliExternalTrackParam.h"
#include "AliPhysicsSelection.h"
#include "TGeoGlobalMagField.h"
#include "AliMagF.h"
#include "AliLog.h"
#include "AliStack.h"
#include "AliCentrality.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliHFEcontainer.h"
#include "AliHFEcuts.h"
#include "AliHFEpid.h"
#include "AliHFEpidBase.h"
#include "AliHFEpidQAmanager.h"
#include "AliHFEtools.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"
#include "AliSelectNonHFE.h"
#include "AliHFEpidTPC.h"
#include "AliAnalysisTaskHFEpACorrelation.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TFile.h"
#include "AliESDHandler.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "TParticle.h"
#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "TVector.h"
#include "stdio.h"
#include "iostream"
#include "fstream"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliVParticle.h"
#include "AliVTrack.h"
#include "AliEventPoolManager.h"
#include "TObjArray.h"
//include to use reader as Lucile does
#include "AliAODHeader.h"
#include "AliAnalysisUtils.h"



// --- ANALYSIS system ---
#include "AliCalorimeterUtils.h"
#include "AliESDEvent.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliMixedEvent.h"
#include "AliAODCaloCluster.h"
#include "AliOADBContainer.h"
#include "AliAnalysisManager.h"

// --- Detector ---

//______________________________________________________________________

//______________________________________________________________________
ClassImp(AliAnalysisTaskHFEpACorrelation)

//______________________________________________________________________
AliAnalysisTaskHFEpACorrelation::AliAnalysisTaskHFEpACorrelation(const char *name)
: AliAnalysisTaskSE(name)
,fCorrelationFlag(0)
,fUseKF(kFALSE)
,fIspp(kFALSE)
,fIsMC(0)
,fUseAlternativeBinnig(kFALSE)
,fAssocWithSPD(kFALSE)
,fIsHFE1(kFALSE)
,fIsHFE2(kFALSE)
,fIsNonHFE(kFALSE)
,fIsFromD(kFALSE)
,fIsFromB(kFALSE)
,fIsFromPi0(kFALSE)
,fIsFromEta(kFALSE)
,fIsFromGamma(kFALSE)
,fpTBins(0)
,fESD(0)
,fAOD(0)
,fVevent(0)
,fPartnerCuts(0)
,fOutputList(0)
,fPidResponse(0)
,fNonHFE(0)
,fIsAOD(kTRUE)
,fCentrality(0)
,fCentralityMin(0)
,fCentralityMax(100)
,fHasCentralitySelection(kFALSE)
,fCentralityHist(0)
,fCentralityHistPass(0)
,fZvtx(0)
,fEstimator(0)
,fUseDCACutforHadrons(kTRUE)
,fEffInc(0)
,fEffBkg(0)
,fEffHadron(0)
,fNevent(0)
,fNevent2(0)
,fPtElec_Inc(0)
,fPtElec_ULS(0)
,fPtElec_LS(0)
,fPtElec_ULS_NoPid(0)
,fPtElec_LS_NoPid(0)
,fPtElec_ULS_weight(0)
,fPtElec_LS_weight(0)
,fPtElec_ULS2_weight(0)
,fPtElec_LS2_weight(0)
,fTOF01(0)
,fTOF02(0)
,fTOF03(0)
,fpid(0)
,fTPC_pt(0)
,fTPC_p(0)
,fTPC_momentum(0)
,fTPC_eta(0)
,fTPC_momentum1(0)
,fTPC_eta1(0)
,fTPCnsigma_pt(0)
,fTOFTPCnsigma_pt(0)
,fTPCnsigma_p(0)
,fTPCnsigma_eta(0)
,fTPCnsigma_phi(0)
,fVtxZ(0)
,fVtxZ_new1(0)
,fVtxZ_new2(0)
,fVtxZ_new3(0)
,fVtxZ_new4(0)
,fzRes1(0)
,fzRes2(0)
,fSPD_track_vtx1(0)
,fSPD_track_vtx2(0)
,fEtad(0)
,fNTracks(0)
,fTrack_Multi(0)
,fNTracks_pt(0)
,fNTracks_eta(0)
,fNTracks_phi(0)
,fTPCNcls_pid(0)
,fCEtaPhi_Inc(0)
,fCEtaPhi_ULS_Weight(0)
,fCEtaPhi_LS_Weight(0)
,fCEtaPhi_ULS_NoP_Weight(0)
,fCEtaPhi_LS_NoP_Weight(0)
,fInvMass(0)
,fInvMassBack(0)
,fDCA(0)
,fDCABack(0)
,fOpAngle(0)
,fOpAngleBack(0)
,fInvMass2(0)
,fInvMass2_weight(0)
,fInvMassBack2(0)
,fInvMassBack2_weight(0)
,fInvMass_pT(0)
,fInvMassBack_pT(0)
,fDCA2(0)
,fDCABack2(0)
,fOpAngle2(0)
,fOpAngleBack2(0)
,fMassCut(0.1)
,fEtaCutMin(-0.8)
,fEtaCutMax(0.8)
,fMinpTElec(0.5)
,fMaxpTElec(6.)
,fChi2Cut(3.5)
,fDCAcutzHadron(0.25)
,fDCAcutrHadron(1)
,fDCAcut(999)//dca between two tracks
,fDCAcutr(999)//dca to vertex
,fDCAcutz(999)//dca to vertex
,fMassCutFlag(kTRUE)
,fAngleCutFlag(kFALSE)
,fAngleCut(999)
,fChi2CutFlag(kFALSE)
,fDCAcutFlag(kFALSE)
,fAssHadronPtMin(0.3)
,fAssHadronPtMax(2.0)
,fPtBackgroundBeforeReco(0)
,fPtMinAsso(0.0)
,fTpcNclsAsso(80)
,fCuts(0)
,fCFM(0)
,fPID()
//Lucile
,fPIDqa(0)
,fMCstack(0)
,fRejectKinkMother(kFALSE)
,fMCtrack(0)
,fMCtrackMother(0)
,fMCtrackGMother(0)
,fMCtrackGGMother(0)
,fMCtrackGGGMother(0)
,fMCarray(0)
,fMCheader(0)
,fMCparticle(0)
,fMCparticle2(0)
,fMCparticleMother(0)
,fMCparticleGMother(0)
,fMCparticleGGMother(0)
,fMCparticleGGGMother(0)
,fEventHandler(0)
,fMCevent(0)
,fPoolMgr(0)
,fPool(0)
,fTracksClone(0)
,fTracks(0)
,fCEtaPhi_Inc_EM(0)
,fCEtaPhi_ULS_Weight_EM(0)
,fCEtaPhi_LS_Weight_EM(0)
,fPoolNevents(0)
,fEventMixingFlag(0)
,fCEtaPhi_Inc_DiHadron(0)
,fPtTrigger_Inc(0)
,fAnalysisUtils(0)

//pT,eta,zvtx
,fNoEtaCutElectronGeneratedSignalPtEtaZvtx(0)
,fEtaCutElectronGeneratedSignalPtEtaZvtx(0)
//,fEtaCutElectronInclusiveGeneratedPtEtaZvtx(0)
//,fEtaCutElectronBKGGeneratedPtEtaZvtx(0)
//,fEtaCutElectronHFeGeneratedPtEtaZvtx(0)
//,fEtaCutElectronOtherGeneratedPtEtaZvtx(0)
//,fEtaCutElectronNoMotherGeneratedPtEtaZvtx(0)
,fEtaCutElectronInclusiveRecoPtEtaZvtx(0)
,fEtaCutElectronBKWithLabelULS(0)
,fEtaCutElectronBKWithLabelLS(0)
,fEtaCutElectronBKNoTag(0)
,fEtaCutElectronRecoHFEMC(0)
,fEtaCutElectronRecoOtherMC(0)
,fMissIDElectronsReco(0)
,fHadronsReco(0)
,fHadronsRecoPP(0)
,fHadronsGenerated(0)

//DPhi MC
,fCEtaPhiNoEtaCutInclusive(0)
,fCEtaPhiNoEtaCutBKG(0)
,fCEtaPhiNoEtaCutHFe(0)
,fCEtaPhiNoEtaCutHFeNoDCA(0)
,fCEtaPhiNoEtaCutOther(0)
,fCEtaPhiNoEtaCutNoMother(0)
,fCEtaPhiCutInclusive(0)
,fCEtaPhiCutBKG(0)
,fCEtaPhiCutHFe(0)
,fCEtaPhiCutOther(0)
,fCEtaPhiCutNoMother(0)
,fCEtaPhi_Back_MC_Tag(0)
,fCEtaPhi_Other_MC_Tag(0)
,fCEtaPhi_HFe_MC_Tag(0)
,fCEtaPhi_MC_NoMother_Tag(0)
,fDCAElectronXY(0)
,fDCAElectronZ(0)
,fElectronNoLabel(0)
,fElectronNoLabelULS(0)
,fElectronNoLabelLS(0)
,fEtaCutElectronBKULSMainSources(0)
,fEtaCutElectronBKLSMainSources(0)
,fEtaCutElectronBKULSOtherSources(0)
,fEtaCutElectronBKLSOtherSources(0)
,fEtaCutElectronHFEULS(0)
,fEtaCutElectronHFELS(0)
,fEtaCutElectronMissIDULS(0)
,fEtaCutElectronMissIDLS(0)
{
    //Named constructor
    // Define input and output slots here
    // Input slot #0 works with a TChain
    
    //exotic
    //fEMCALRecoUtils  = new AliEMCALRecoUtils();
    fPID = new AliHFEpid("hfePid");

    DefineInput(0, TChain::Class());
    // Output slot #0 id reserved by the base class for AOD
    // Output slot #1 writes into a TH1 container
    // DefineOutput(1, TH1I::Class());
    DefineOutput(1, TList::Class());
    //  DefineOutput(3, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskHFEpACorrelation::AliAnalysisTaskHFEpACorrelation()
: AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisTaskHFEpACorrelation")
,fCorrelationFlag(0)
,fUseKF(kFALSE)
,fIspp(kFALSE)
,fIsMC(0)
,fUseAlternativeBinnig(kFALSE)
,fAssocWithSPD(kFALSE)
,fIsHFE1(kFALSE)
,fIsHFE2(kFALSE)
,fIsNonHFE(kFALSE)
,fIsFromD(kFALSE)
,fIsFromB(kFALSE)
,fIsFromPi0(kFALSE)
,fIsFromEta(kFALSE)
,fIsFromGamma(kFALSE)
,fpTBins(0)
,fESD(0)
,fAOD(0)
,fVevent(0)
,fPartnerCuts(0)
,fOutputList(0)
,fPidResponse(0)
,fNonHFE(0)
,fIsAOD(kTRUE)
,fCentrality(0)
,fCentralityMin(0)
,fCentralityMax(100)
,fHasCentralitySelection(kFALSE)
,fCentralityHist(0)
,fCentralityHistPass(0)
,fZvtx(0)
,fEstimator(0)
,fUseDCACutforHadrons(kTRUE)
,fEffInc(0)
,fEffBkg(0)
,fEffHadron(0)
,fNevent(0)
,fNevent2(0)
,fPtElec_Inc(0)
,fPtElec_ULS(0)
,fPtElec_LS(0)
,fPtElec_ULS_NoPid(0)
,fPtElec_LS_NoPid(0)
,fPtElec_ULS_weight(0)
,fPtElec_LS_weight(0)
,fPtElec_ULS2_weight(0)
,fPtElec_LS2_weight(0)
,fTOF01(0)
,fTOF02(0)
,fTOF03(0)
,fpid(0)
,fTPC_pt(0)
,fTPC_p(0)
,fTPC_momentum(0)
,fTPC_eta(0)
,fTPC_momentum1(0)
,fTPC_eta1(0)
,fTPCnsigma_pt(0)
,fTOFTPCnsigma_pt(0)
,fTPCnsigma_p(0)
,fTPCnsigma_eta(0)
,fTPCnsigma_phi(0)
,fVtxZ(0)
,fVtxZ_new1(0)
,fVtxZ_new2(0)
,fVtxZ_new3(0)
,fVtxZ_new4(0)
,fzRes1(0)
,fzRes2(0)
,fSPD_track_vtx1(0)
,fSPD_track_vtx2(0)
,fEtad(0)
,fNTracks(0)
,fTrack_Multi(0)
,fNTracks_pt(0)
,fNTracks_eta(0)
,fNTracks_phi(0)
,fTPCNcls_pid(0)
,fCEtaPhi_Inc(0)
,fCEtaPhi_ULS_Weight(0)
,fCEtaPhi_LS_Weight(0)
,fCEtaPhi_ULS_NoP_Weight(0)
,fCEtaPhi_LS_NoP_Weight(0)
,fInvMass(0)
,fInvMassBack(0)
,fDCA(0)
,fDCABack(0)
,fOpAngle(0)
,fOpAngleBack(0)
,fInvMass2(0)
,fInvMass2_weight(0)
,fInvMassBack2(0)
,fInvMassBack2_weight(0)
,fInvMass_pT(0)
,fInvMassBack_pT(0)
,fDCA2(0)
,fDCABack2(0)
,fOpAngle2(0)
,fOpAngleBack2(0)
,fMassCut(0.1)
,fEtaCutMin(-0.8)
,fEtaCutMax(0.8)
,fMinpTElec(0.5)
,fMaxpTElec(6.)
,fChi2Cut(3.5)
,fDCAcutzHadron(0.25)
,fDCAcutrHadron(1)
,fDCAcut(999)//dca between two tracks
,fDCAcutr(999)//dca to vertex
,fDCAcutz(999)//dca to vertex
,fMassCutFlag(kTRUE)
,fAngleCutFlag(kFALSE)
,fAngleCut(999)
,fChi2CutFlag(kFALSE)
,fDCAcutFlag(kFALSE)
,fAssHadronPtMin(0.3)
,fAssHadronPtMax(2.0)
,fPtBackgroundBeforeReco(0)
,fPtMinAsso(0.0)
,fTpcNclsAsso(80)
,fCuts(0)
,fCFM(0)
,fPID()
//Lucile
,fPIDqa(0)
,fMCstack(0)
,fRejectKinkMother(kFALSE)
,fMCtrack(0)
,fMCtrackMother(0)
,fMCtrackGMother(0)
,fMCtrackGGMother(0)
,fMCtrackGGGMother(0)
,fMCarray(0)
,fMCheader(0)
,fMCparticle(0)
,fMCparticle2(0)
,fMCparticleMother(0)
,fMCparticleGMother(0)
,fMCparticleGGMother(0)
,fMCparticleGGGMother(0)
,fEventHandler(0)
,fMCevent(0)
,fPoolMgr(0)
,fPool(0)
,fTracksClone(0)
,fTracks(0)
,fCEtaPhi_Inc_EM(0)
,fCEtaPhi_ULS_Weight_EM(0)
,fCEtaPhi_LS_Weight_EM(0)
,fPoolNevents(0)
,fEventMixingFlag(0)
,fCEtaPhi_Inc_DiHadron(0)
,fPtTrigger_Inc(0)
,fAnalysisUtils(0)

//pT,eta,zvtx
,fNoEtaCutElectronGeneratedSignalPtEtaZvtx(0)
,fEtaCutElectronGeneratedSignalPtEtaZvtx(0)
,fEtaCutElectronInclusiveRecoPtEtaZvtx(0)
,fEtaCutElectronBKWithLabelULS(0)
,fEtaCutElectronBKWithLabelLS(0)
,fEtaCutElectronBKNoTag(0)
,fEtaCutElectronRecoHFEMC(0)
,fEtaCutElectronRecoOtherMC(0)
,fMissIDElectronsReco(0)
,fHadronsReco(0)
,fHadronsRecoPP(0)
,fHadronsGenerated(0)

//DPhi MC
,fCEtaPhiNoEtaCutInclusive(0)
,fCEtaPhiNoEtaCutBKG(0)
,fCEtaPhiNoEtaCutHFe(0)
,fCEtaPhiNoEtaCutHFeNoDCA(0)
,fCEtaPhiNoEtaCutOther(0)
,fCEtaPhiNoEtaCutNoMother(0)
,fCEtaPhiCutInclusive(0)
,fCEtaPhiCutBKG(0)
,fCEtaPhiCutHFe(0)
,fCEtaPhiCutOther(0)
,fCEtaPhiCutNoMother(0)
,fCEtaPhi_Back_MC_Tag(0)
,fCEtaPhi_Other_MC_Tag(0)
,fCEtaPhi_HFe_MC_Tag(0)
,fCEtaPhi_MC_NoMother_Tag(0)
,fDCAElectronXY(0)
,fDCAElectronZ(0)
,fElectronNoLabel(0)
,fElectronNoLabelULS(0)
,fElectronNoLabelLS(0)
,fEtaCutElectronBKULSMainSources(0)
,fEtaCutElectronBKLSMainSources(0)
,fEtaCutElectronBKULSOtherSources(0)
,fEtaCutElectronBKLSOtherSources(0)
,fEtaCutElectronHFEULS(0)
,fEtaCutElectronHFELS(0)
,fEtaCutElectronMissIDULS(0)
,fEtaCutElectronMissIDLS(0)
{
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
}

//______________________________________________________________________
AliAnalysisTaskHFEpACorrelation::~AliAnalysisTaskHFEpACorrelation()
{
    //Destructor
    delete fOutputList;
    delete fPID;
    delete fCFM;
    delete fPIDqa;
    delete fNonHFE;
    delete fPartnerCuts;
    delete fAnalysisUtils;
}

//______________________________________________________________________
//Create Output Objects
//Here we can define the histograms and others output files
//Called once
void AliAnalysisTaskHFEpACorrelation::UserCreateOutputObjects()
{
    //______________________________________________________________________
    //Initialize PID
    
    fAnalysisUtils = new AliAnalysisUtils;
    //
    fPartnerCuts = new AliESDtrackCuts();
    
    
    if(!fPID->GetNumberOfPIDdetectors())
    {
        fPID->AddDetector("TPC", 0);
    }

    fPID->SortDetectors();
    
    fPIDqa = new AliHFEpidQAmanager();
    fPIDqa->Initialize(fPID);
    
    fNonHFE = new AliSelectNonHFE();
    //______________________________________________________________________
    
    
    //______________________________________________________________________
    //Initialize correction Framework and Cuts
    fCFM = new AliCFManager;
    const Int_t kNcutSteps = AliHFEcuts::kNcutStepsMCTrack + AliHFEcuts::kNcutStepsRecTrack + AliHFEcuts::kNcutStepsDETrack;
    fCFM->SetNStepParticle(kNcutSteps);
    for(Int_t istep = 0; istep < kNcutSteps; istep++) fCFM->SetParticleCutsList(istep, NULL);
    
    if(!fCuts)
    {
        AliWarning("================ \n Cuts not available. Default cuts will be used ================ \n");
        AliWarning("================ \n You should really consider to abort the analysis ================ \n");
        fCuts = new AliHFEcuts;
        fCuts->CreateStandardCuts();
    }

    //from Andrea Dubla
    //fCuts->SetAOD();
    
    fCuts->Initialize(fCFM);

    
    ///Output Tlist
    //Create TList
    fOutputList = new TList();
    fOutputList->SetOwner();
    

    
    //PIDqa
    //fOutputList->Add(fPIDqa->MakeList("PIDQA"));
    
    //Store the number of events
    //Define the histo
    fNevent = new TH1F("fNevent","Number of Events",30,0,30);
    fNevent2 = new TH1F("fNevent2","Number of Events 2",30,0,30);
    fOutputList->Add(fNevent);
    fOutputList->Add(fNevent2);
    //And then, add to the output list


    fpid = new TH1F("fpid","PID flag",5,0,5);
    fOutputList->Add(fpid);
    
    //pt Distribution
    fPtElec_Inc = new TH1F("fPtElec_Inc","Inclusive Electrons; p_{T} (GeV/c); Count",110,0.5,6);
    fPtElec_ULS = new TH1F("fPtElec_ULS","ULS; p_{T} (GeV/c); Count",110,0.5,6);
    fPtElec_LS = new TH1F("fPtElec_LS","LS; p_{T} (GeV/c); Count",110,0.5,6);
    
    if (fIsMC)
    {
        fPtElec_ULS_NoPid = new TH1F("fPtElec_ULS_NoPid","ULS; p_{T} (GeV/c); Count",110,0.5,6);
        fPtElec_LS_NoPid = new TH1F("fPtElec_LS_NoPid","LS; p_{T} (GeV/c); Count",110,0.5,6);
        fOutputList->Add(fPtElec_ULS_NoPid);
        fOutputList->Add(fPtElec_LS_NoPid);
    }
    
    fOutputList->Add(fPtElec_Inc);
    fOutputList->Add(fPtElec_ULS);
    fOutputList->Add(fPtElec_LS);
    
    //fTOF01 = new TH2F("fTOF01","",100,-20,20,100,-20,20);
    //fTOF02 = new TH2F("fTOF02","",100,-20,20,100,-20,20);
    //fTOF03 = new TH2F("fTOF03","",100,-20,20,100,-20,20);
    
    fPtTrigger_Inc = new TH1F("fPtTrigger_Inc","pT dist for Hadron Contamination; p_{t} (GeV/c); Count",110,0.5,6);
    
    fOutputList->Add(fPtTrigger_Inc);
    //fOutputList->Add(fTOF01);
    //fOutputList->Add(fTOF02);
    //fOutputList->Add(fTOF03);
    
    fVtxZ_new1= new  TH1F("fVtxZ_new1","fVtxZ_new1",100, -50,50);
    fVtxZ_new2= new  TH1F("fVtxZ_new2","fVtxZ_new2",100, -50,50);
    fVtxZ_new3= new  TH1F("fVtxZ_new3","fVtxZ_new3",100, -50,50);
    fVtxZ_new4= new  TH1F("fVtxZ_new4","fVtxZ_new4",100, -50,50);
    
    fOutputList->Add(fVtxZ_new1);
    fOutputList->Add(fVtxZ_new2);
    fOutputList->Add(fVtxZ_new3);
    fOutputList->Add(fVtxZ_new4);
    
    fzRes1= new  TH1F("fzRes1","fzRes1",100, 0,1);
    fzRes2= new  TH1F("fzRes2","fzRes2",100, 0,1);
    fSPD_track_vtx1= new  TH1F("fSPD_track_vtx1","fSPD_track_vtx1",100, -5,5);
    fSPD_track_vtx2= new  TH1F("fSPD_track_vtx2","fSPD_track_vtx2",100, -5,5);
    
    fOutputList->Add(fzRes1);
    fOutputList->Add(fzRes2);
    fOutputList->Add(fSPD_track_vtx1);
    fOutputList->Add(fSPD_track_vtx2);
    
    
    //General Histograms
    
    //Steps
    //Step 1: Before Track cuts
    //Step 2: Before PID
    //Step 3: After PID
    
    fTPC_p = new TH2F *[3];
    fTPCnsigma_p = new TH2F *[3];
    fVtxZ= new  TH1F *[3];
    fEtad= new  TH1F *[3];
    fNTracks= new  TH1F *[3];
    
    fNTracks_pt= new  TH2F *[3];
    fNTracks_eta= new  TH2F *[3];
    fNTracks_phi= new  TH2F *[3];
    fTPCNcls_pid = new  TH2F *[3];
    
    
    fTPC_momentum = new TH2F("fTPC_momentum",";p (GeV/c);TPC dE/dx (a. u.)",100,0,10,100,-20,20);
    fTPC_eta = new TH2F("fTPC_eta",";eta (GeV/c);TPC dE/dx (a. u.)",80,-2,2,220,-20,200);
    fOutputList->Add(fTPC_momentum);
    fOutputList->Add(fTPC_eta);
    
    fTPC_momentum1 = new TH2F("fTPC_momentum1",";p (GeV/c);TPC dE/dx (a. u.)",100,0,10,100,-20,20);
    fTPC_eta1 = new TH2F("fTPC_eta1",";eta (GeV/c);TPC dE/dx (a. u.)",80,-2,2,220,-20,200);
    fOutputList->Add(fTPC_momentum1);
    fOutputList->Add(fTPC_eta1);
    
    for(Int_t i = 0; i < 3; i++)
    {
        fTPC_p[i] = new TH2F(Form("fTPC_p%d",i),";pt (GeV/c);TPC dE/dx (a. u.)",100,0.3,15,100,-20,200);
        fTPCnsigma_p[i] = new TH2F(Form("fTPCnsigma_p%d",i),";p (GeV/c);TPC Electron N#sigma",100,0.3,15,100,-15,10);
        
        fOutputList->Add(fTPC_p[i]);
        fOutputList->Add(fTPCnsigma_p[i]);
        
        fVtxZ[i]= new  TH1F(Form("fVtxZ%d",i),"VtxZ",100, -50,50);
        fEtad[i]= new  TH1F(Form("fEtad%d",i),"Eta distribution",100, -1.2,1.2);
        fNTracks[i]= new  TH1F(Form("fNTracks%d",i),"NTracks",1000, 0,5000);
        
        fOutputList->Add(fVtxZ[i]);
        fOutputList->Add(fEtad[i]);
        fOutputList->Add(fNTracks[i]);
        
        fNTracks_pt[i]= new  TH2F(Form("fNTracks_pt%d",i),"NTracks vs. pt",100, 0,5000, 100, 0, 10);
        fNTracks_eta[i]= new  TH2F(Form("fNTracks_eta%d",i),"NTracks vs. pt",100, 0,5000, 100, -1.0, 1.0);
        fNTracks_phi[i]= new  TH2F(Form("fNTracks_phi%d",i),"NTracks vs. pt",100, 0,5000, 100, 0, 5.0);

        
        fOutputList->Add(fNTracks_pt[i]);
        fOutputList->Add(fNTracks_eta[i]);
        fOutputList->Add(fNTracks_phi[i]);
        
    }

    for(Int_t i = 0; i < 4; i++)
    {
        fTPCNcls_pid[i]= new TH2F(Form("fTPCNcls_pid%d",i),"fTPCNcls_pid;NCls;NCls for PID",200,0,200,200,0,200);
        fOutputList->Add(fTPCNcls_pid[i]);
    }
    //pt bin
    

    fTPC_pt = new TH1F *[fpTBins.GetSize()];
    

    fTPCnsigma_pt = new TH1F *[fpTBins.GetSize()];
    fTOFTPCnsigma_pt = new TH2F *[fpTBins.GetSize()];
    
    
    fDCAElectronXY = new TH1F *[2];
    fDCAElectronZ = new TH1F *[2];
    
    for (Int_t i = 0; i < 2; i++)
    {
        fDCAElectronXY[i] = new TH1F(Form("fDCAElectronXY%d",i),"",40,0,4);
        fDCAElectronZ[i] = new TH1F(Form("fDCAElectronZ%d",i),"",40,0,4);
        fOutputList->Add(fDCAElectronXY[i]);
        fOutputList->Add(fDCAElectronZ[i]);
    }
    
    if (fIsMC)
    {
        //binnig for unfolding
        
        Double_t EtaBins[] = {-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8}; //16 bins
        Int_t NEtaBins = 16;
        
        Double_t ZVtxBins[] = {-10, -7, -5, -3, -1, 1, 3, 5, 7, 10.01}; //9 bins
        Int_t nZvtxBinsMC = 9;
        
        Double_t pTBinsH[] = {0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.8,2.0,2.5,3,3.5,4,5,6}; //21 bins
        Int_t NpTbinsH = 21;
        
        Int_t NumberofBins = fpTBins.GetSize();
        
        fCEtaPhiNoEtaCutHFe = new TH2F *[NumberofBins];
        fCEtaPhiCutHFe = new TH2F *[NumberofBins];
        
        fCEtaPhi_Back_MC_Tag = new TH2F *[NumberofBins];
        fCEtaPhi_Other_MC_Tag  = new TH2F *[NumberofBins];
        fCEtaPhi_HFe_MC_Tag = new TH2F *[NumberofBins];
        fCEtaPhi_MC_NoMother_Tag = new TH2F *[NumberofBins];
        
        fNoEtaCutElectronGeneratedSignalPtEtaZvtx = new TH1F("fNoEtaCutElectronGeneratedSignalPtEtaZvtx", "",110,0.5,6);
        fOutputList->Add(fNoEtaCutElectronGeneratedSignalPtEtaZvtx);
        
        fEtaCutElectronGeneratedSignalPtEtaZvtx = new TH1F("fEtaCutElectronGeneratedSignalPtEtaZvtx", "",110,0.5,6);
        fOutputList->Add(fEtaCutElectronGeneratedSignalPtEtaZvtx);

        
        
        fEtaCutElectronBKWithLabelULS = new TH1F("fEtaCutElectronBKWithLabelULS", "" ,110,0.5,6);
        
        fEtaCutElectronBKWithLabelLS = new TH1F("fEtaCutElectronBKWithLabelLS", "", 110,0.5,6);
        
        
        fEtaCutElectronInclusiveRecoPtEtaZvtx = new TH1F("fEtaCutElectronInclusiveRecoPtEtaZvtx", "", 110,0.5,6);
        //Backgound for data with MC info
        fEtaCutElectronBKNoTag = new TH1F("fEtaCutElectronBKNoTag", "", 110,0.5,6);
        
        //Tagged Background
        
        //HFe with MC information
        fEtaCutElectronRecoHFEMC = new TH1F("fEtaCutElectronRecoHFEMC", "", 110,0.5,6);
        //Others with MC information
        fEtaCutElectronRecoOtherMC = new TH1F("fEtaCutElectronRecoOtherMC", "", 110,0.5,6);
        
        //Miss ID
        fMissIDElectronsReco = new TH1F("fMissIDElectronsReco", "", 110,0.5,6);
        
        fElectronNoLabel = new TH1F("fElectronNoLabel", "", 110,0.5,6);
        fElectronNoLabelULS = new TH1F("fElectronNoLabelULS", "", 110,0.5,6);
        fElectronNoLabelLS = new TH1F("fElectronNoLabelLS", "", 110,0.5,6);
        fEtaCutElectronBKULSMainSources = new TH1F("fEtaCutElectronBKULSMainSources", "", 110,0.5,6);
        fEtaCutElectronBKLSMainSources = new TH1F("fEtaCutElectronBKLSMainSources", "", 110,0.5,6);
        fEtaCutElectronBKULSOtherSources = new TH1F("fEtaCutElectronBKULSOtherSources", "", 110,0.5,6);
        fEtaCutElectronBKLSOtherSources = new TH1F("fEtaCutElectronBKLSOtherSources", "", 110,0.5,6);
        fEtaCutElectronHFEULS = new TH1F("fEtaCutElectronHFEULS", "", 110,0.5,6);
        fEtaCutElectronHFELS = new TH1F("fEtaCutElectronHFELS", "", 110,0.5,6);
        fEtaCutElectronMissIDULS = new TH1F("fEtaCutElectronMissIDULS", "", 110,0.5,6);
        fEtaCutElectronMissIDLS = new TH1F("fEtaCutElectronMissIDLS", "", 110,0.5,6);
        
        
        fOutputList->Add(fEtaCutElectronInclusiveRecoPtEtaZvtx);
        fOutputList->Add(fEtaCutElectronBKNoTag);
        fOutputList->Add(fEtaCutElectronBKWithLabelLS);
        fOutputList->Add(fEtaCutElectronBKWithLabelULS);
        fOutputList->Add(fEtaCutElectronRecoHFEMC);
        fOutputList->Add(fEtaCutElectronRecoOtherMC);
        fOutputList->Add(fMissIDElectronsReco);
        
        fOutputList->Add(fElectronNoLabel);
        fOutputList->Add(fElectronNoLabelULS);
        fOutputList->Add(fElectronNoLabelLS);
        fOutputList->Add(fEtaCutElectronBKULSMainSources);
        fOutputList->Add(fEtaCutElectronBKLSMainSources);
        fOutputList->Add(fEtaCutElectronBKULSOtherSources);
        fOutputList->Add(fEtaCutElectronBKLSOtherSources);
        fOutputList->Add(fEtaCutElectronHFEULS);
        fOutputList->Add(fEtaCutElectronHFELS);
        fOutputList->Add(fEtaCutElectronMissIDULS);
        fOutputList->Add(fEtaCutElectronMissIDLS);
        
        //Hadrons
        fHadronsReco = new TH3F("fHadronsReco",  "p_{T} x #eta x  ZVtx; p_{T}; #eta; Zvtx", NpTbinsH, pTBinsH, NEtaBins, EtaBins, nZvtxBinsMC, ZVtxBins);
        fHadronsRecoPP = new TH3F("fHadronsRecoPP",  "p_{T} x #eta x  ZVtx; p_{T}; #eta; Zvtx", NpTbinsH, pTBinsH, NEtaBins, EtaBins, nZvtxBinsMC, ZVtxBins);
        fHadronsGenerated = new TH3F("fHadronsGenerated", "p_{T} x #eta x  ZVtx; p_{T}; #eta; Zvtx", NpTbinsH, pTBinsH, NEtaBins, EtaBins, nZvtxBinsMC, ZVtxBins);
        
        fOutputList->Add(fHadronsReco);
        fOutputList->Add(fHadronsRecoPP);
        fOutputList->Add(fHadronsGenerated);
        
        
        
    }


    fInvMass = new TH1F("fInvMass","",100,0,0.3);
    fInvMassBack = new TH1F("fInvMassBack","",100,0,0.3);
    fDCA = new TH1F("fDCA","",100,0,1);
    fDCABack = new TH1F("fDCABack","",100,0,1);
    fOpAngle = new TH1F("fOpAngle","",100,0,0.5);
    fOpAngleBack = new TH1F("fOpAngleBack","",100,0,0.5);
    
    if(fCorrelationFlag)
    {
        Int_t NumberBins = fpTBins.GetSize();
        
        fCEtaPhi_Inc = new TH2F *[NumberBins];
        fCEtaPhi_Inc_DiHadron = new TH2F *[NumberBins];
        fCEtaPhi_ULS_Weight = new TH2F *[NumberBins];
        fCEtaPhi_LS_Weight = new TH2F *[NumberBins];
        fCEtaPhi_ULS_NoP_Weight = new TH2F *[NumberBins];
        fCEtaPhi_LS_NoP_Weight = new TH2F *[NumberBins];
        
        fCEtaPhi_Inc_EM = new TH2F *[NumberBins];
        fCEtaPhi_ULS_Weight_EM = new TH2F *[NumberBins];
        fCEtaPhi_LS_Weight_EM = new TH2F *[NumberBins];
        
        fOutputList->Add(fInvMass);
        fOutputList->Add(fInvMassBack);
        fOutputList->Add(fDCA);
        fOutputList->Add(fDCABack);
        fOutputList->Add(fOpAngle);
        fOutputList->Add(fOpAngleBack);
    }
    

    
    
    for(Int_t i = 0; i < fpTBins.GetSize()-1; i++)
    {
        fTPC_pt[i] = new TH1F(Form("fTPC_pt%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;TPC Electron N#sigma;Count",fpTBins.At(i),fpTBins.At(i+1)),100,20,200);
        fTPCnsigma_pt[i] = new TH1F(Form("fTPCnsigma_pt%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;TPC Electron N#sigma;Count",fpTBins.At(i),fpTBins.At(i+1)),100,-15,10);
        fTOFTPCnsigma_pt[i] = new TH2F(Form("fTOFTPCnsigma_pt%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c; TOF Electron N#sigma; TPC Electron N#sigma; Count",fpTBins.At(i),fpTBins.At(i+1)), 150,-10,20,180,-12,8  );
        
        fOutputList->Add(fTPC_pt[i]);
        fOutputList->Add(fTPCnsigma_pt[i]);
        fOutputList->Add(fTOFTPCnsigma_pt[i]);
        
        if(fCorrelationFlag)
        {
            
            Int_t BinSize = 40;
            fCEtaPhi_Inc[i] = new TH2F(Form("fCEtaPhi_Inc%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
            fCEtaPhi_Inc_DiHadron[i] = new TH2F(Form("fCEtaPhi_Inc_DiHadron%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
            
            
            fCEtaPhi_ULS_Weight[i] = new TH2F(Form("fCEtaPhi_ULS_Weight%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
            fCEtaPhi_LS_Weight[i] = new TH2F(Form("fCEtaPhi_LS_Weight%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
            fCEtaPhi_ULS_NoP_Weight[i] = new TH2F(Form("fCEtaPhi_ULS_NoP_Weight%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
            fCEtaPhi_LS_NoP_Weight[i] = new TH2F(Form("fCEtaPhi_LS_NoP_Weight%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
            
        
            if (fIsMC)
            {

                fCEtaPhiNoEtaCutHFe[i] = new TH2F(Form("fCEtaPhiNoEtaCutHFe%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
                
                fOutputList->Add(fCEtaPhiNoEtaCutHFe[i]);
                
                fCEtaPhiCutHFe[i] = new TH2F(Form("fCEtaPhiCutHFe%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
                
                fOutputList->Add(fCEtaPhiCutHFe[i]);
                
                fCEtaPhi_Back_MC_Tag[i] = new TH2F(Form("fCEtaPhi_Back_MC_Tag%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
                fCEtaPhi_Other_MC_Tag[i] = new TH2F(Form("fCEtaPhi_Other_MC_Tag%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
                fCEtaPhi_HFe_MC_Tag[i] = new TH2F(Form("fCEtaPhi_HFe_MC_Tag%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
                fCEtaPhi_MC_NoMother_Tag[i] = new TH2F(Form("fCEtaPhi_MC_NoMother%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
                
                fOutputList->Add(fCEtaPhi_Back_MC_Tag[i]);
                fOutputList->Add(fCEtaPhi_Other_MC_Tag[i]);
                fOutputList->Add(fCEtaPhi_HFe_MC_Tag[i]);
                fOutputList->Add(fCEtaPhi_MC_NoMother_Tag[i]);
                

            }
            
            fOutputList->Add(fCEtaPhi_Inc[i]);
            fOutputList->Add(fCEtaPhi_Inc_DiHadron[i]);
            
            fOutputList->Add(fCEtaPhi_ULS_Weight[i]);
            fOutputList->Add(fCEtaPhi_LS_Weight[i]);
            fOutputList->Add(fCEtaPhi_ULS_NoP_Weight[i]);
            fOutputList->Add(fCEtaPhi_LS_NoP_Weight[i]);
            
            
            
            if(fEventMixingFlag)
            {
                Int_t BinSize = 40;
                fCEtaPhi_Inc_EM[i] = new TH2F(Form("fCEtaPhi_Inc_EM%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
                
                fCEtaPhi_ULS_Weight_EM[i] = new TH2F(Form("fCEtaPhi_ULS_Weight_EM%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
                fCEtaPhi_LS_Weight_EM[i] = new TH2F(Form("fCEtaPhi_LS_Weight_EM%d",i),Form("%1.2f < p_{t} < %1.2f GeV/c;#Delta#varphi (rad);#Delta#eta",fpTBins.At(i),fpTBins.At(i+1)),BinSize,-0.5*TMath::Pi(),1.5*TMath::Pi(),BinSize,-2,2);
                
                fOutputList->Add(fCEtaPhi_Inc_EM[i]);
                fOutputList->Add(fCEtaPhi_ULS_Weight_EM[i]);
                fOutputList->Add(fCEtaPhi_LS_Weight_EM[i]);
            }
        }
    }
    
    //pt integrated
    fTPCnsigma_eta = new TH2F("fTPCnsigma_eta",";Pseudorapidity #eta; TPC signal - <TPC signal>_{elec} (#sigma)",100,-0.9,0.9,100,-15,15);
    fTPCnsigma_phi = new TH2F("fTPCnsigma_phi",";Azimuthal Angle #phi; TPC signal - <TPC signal>_{elec} (#sigma)",100,0,2*TMath::Pi(),100,-15,15);
    
    fOutputList->Add(fTPCnsigma_eta);
    fOutputList->Add(fTPCnsigma_phi);

    
    
    
    //__________________________________________________________________
    //Efficiency studies
    if(fIsMC)
    {
        fPtBackgroundBeforeReco = new TH1F("fPtBackgroundBeforeReco",";p_{T} (GeV/c);Count",110,0.5,6);
        fOutputList->Add(fPtBackgroundBeforeReco);
     }
    
    if(!fIspp){
        fCentralityHist = new TH1F("fCentralityHist",";Centrality (%); Count",100,0,100);
        fCentralityHistPass = new TH1F("fCentralityHistPass",";Centrality (%); Count",100,0,100);
        fOutputList->Add(fCentralityHist);
        fOutputList->Add(fCentralityHistPass);
        
    }
    
    for (Int_t i=0; i < fOutputList->GetEntries(); ++i)
    {
        TH1 *h1 = dynamic_cast<TH1*>(fOutputList->At(i));
        if (h1)
        {
            h1->Sumw2();
        }
        else printf("Failed to evaluate Sumw2() \n");
    }
    
    
    //______________________________________________________________________
    //Mixed event analysis -- it was removed from is pp because
    if(fEventMixingFlag && fCorrelationFlag)
    {
        fPoolNevents = new TH1F("fPoolNevents","Event Mixing Statistics; Number of events; Count",1000,0,1000);
        fOutputList->Add(fPoolNevents);
        
        Int_t trackDepth = 2000; // number of objects (tracks) kept per event buffer bin. Once the number of stored objects (tracks) is above that limit, the oldest ones are removed.
        Int_t poolsize   = 1000;  // Maximum number of events, ignored in the present implemented of AliEventPoolManager
        
        Int_t nCentralityBins  = 15;
        Double_t centralityBins[] = { 0, 1, 2, 3, 4, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100.1 };
        
        Int_t nZvtxBins  = 9;
        Double_t vertexBins[] = {-10, -7, -5, -3, -1, 1, 3, 5, 7, 10.01};
        
        fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins, (Double_t*) centralityBins, nZvtxBins, (Double_t*) vertexBins);
    }
    
    //______________________________________________________________________
    
    PostData(1, fOutputList);
    
    ///______________________________________________________________________
}

//______________________________________________________________________
//Main loop
//Called for each event
void AliAnalysisTaskHFEpACorrelation::UserExec(Option_t *)
{
    //Check Event
    fESD = dynamic_cast<AliESDEvent*>(InputEvent());
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    
    if(!(fESD || fAOD))
    {
        printf("ERROR: fESD & fAOD not available\n");
        return;
    }
    
    fVevent = dynamic_cast<AliVEvent*>(InputEvent());
    
    if(!fVevent)
    {
        printf("ERROR: fVEvent not available\n");
        return;
    }
    
    //Check Cuts
    if(!fCuts)
    {
        AliError("HFE cuts not available");
        return;
    }
    //Check PID
    if(!fPID->IsInitialized())
    {
        // Initialize PID with the given run number
        AliWarning("PID not initialised, get from Run no");
        
        if(fIsAOD)
        {
            fPID->InitializePID(fAOD->GetRunNumber());
        }
        else
        {
            fPID->InitializePID(fESD->GetRunNumber());
        }
    }
    
    //PID response
    fPidResponse = fInputHandler->GetPIDResponse();
    
    
    //Check PID response
    if(!fPidResponse)
    {
        AliDebug(1, "Using default PID Response");
        fPidResponse = AliHFEtools::GetDefaultPID(kFALSE, fInputEvent->IsA() == AliAODEvent::Class());
    }
    
    fPID->SetPIDResponse(fPidResponse);
    
    fCFM->SetRecEventInfo(fVevent);
    
    Double_t *fListOfmotherkink = 0;
    Int_t fNumberOfVertices = 0;
    Int_t fNumberOfMotherkink = 0;
    
    
    //total event before event selection
    fNevent->Fill(1);
    
    //______________________________________________________________________
    //Vertex Selection
    if(!fIspp){
        
        fNevent2->Fill(8);
        
        if(fIsAOD)
        {
            const AliAODVertex* trkVtx = fAOD->GetPrimaryVertex();
            
            
            Float_t zvtx = -100;
            zvtx=trkVtx->GetZ();
            fZvtx = zvtx;
            //x
            Float_t xvtx = -100;
            xvtx=trkVtx->GetX();
            
            //y
            Float_t yvtx = -100;
            yvtx=trkVtx->GetY();
            
            
            //events without vertex
            if(zvtx==0 && xvtx==0 && yvtx==0) fNevent2->Fill(0);
            
            
            if(!trkVtx || trkVtx->GetNContributors()<=0) return;
            TString vtxTtl = trkVtx->GetTitle();
            if(!vtxTtl.Contains("VertexerTracks")) return;
            //Float_t zvtx = trkVtx->GetZ();
            
            trkVtx->GetZ();
            fZvtx = zvtx;
            
            
            fVtxZ_new1->Fill(fZvtx);
            
            const AliAODVertex* spdVtx = fAOD->GetPrimaryVertexSPD();
            if(spdVtx->GetNContributors()<=0) return;
            TString vtxTyp = spdVtx->GetTitle();
            Double_t cov[6]={0};
            spdVtx->GetCovarianceMatrix(cov);
            Double_t zRes = TMath::Sqrt(cov[5]);
            
            fzRes1->Fill(zRes);
            //Yvonne e-mail from 12 June 2015 says it has a bug on "vertexer:Z".
            //if(vtxTyp.Contains("vertexer:Z") && (zRes>0.25)) return;
            
            //new line:
            if (spdVtx->IsFromVertexerZ() && (zRes>0.25)) return;
            
            fzRes2->Fill(zRes);
            
            fSPD_track_vtx1->Fill(spdVtx->GetZ() - trkVtx->GetZ());
            if(TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ())>0.5) return;
            fSPD_track_vtx2->Fill(spdVtx->GetZ() - trkVtx->GetZ());
            
            
            
            if(TMath::Abs(zvtx) > 10){
                fNevent2->Fill(2);
                fVtxZ_new2->Fill(fZvtx);
                return;
            }
            
            //if(fabs(zvtx>10.0))return;
            
            fVtxZ_new3->Fill(fZvtx);
            fNevent2->Fill(4);
            
            //Look for kink mother for AOD
            
            fNumberOfVertices = 0;
            fNumberOfMotherkink = 0;
            
            if(fIsAOD)
            {
                fNumberOfVertices = fAOD->GetNumberOfVertices();
                
                fListOfmotherkink = new Double_t[fNumberOfVertices];
                
                for(Int_t ivertex=0; ivertex < fNumberOfVertices; ivertex++)
                {
                    AliAODVertex *aodvertex = fAOD->GetVertex(ivertex);
                    if(!aodvertex) continue;
                    if(aodvertex->GetType()==AliAODVertex::kKink)
                    {
                        AliAODTrack *mother1 = (AliAODTrack *) aodvertex->GetParent();
                        if(!mother1) continue;
                        Int_t idmother = mother1->GetID();
                        fListOfmotherkink[fNumberOfMotherkink] = idmother;
                        fNumberOfMotherkink++;
                    }
                }
            }
        }
        
        else
        {
            
            
            
            /// ESD
            const AliESDVertex* trkVtx = fESD->GetPrimaryVertex();
            if(!trkVtx || trkVtx->GetNContributors()<=0) return;
            TString vtxTtl = trkVtx->GetTitle();
            if(!vtxTtl.Contains("VertexerTracks")) return;
            Float_t zvtx = -100;
            zvtx=trkVtx->GetZ();
            
            
            const AliESDVertex* spdVtx = fESD->GetPrimaryVertexSPD();
            if(spdVtx->GetNContributors()<=0) return;
            TString vtxTyp = spdVtx->GetTitle();
            Double_t cov[6]={0};
            spdVtx->GetCovarianceMatrix(cov);
            Double_t zRes = TMath::Sqrt(cov[5]);
            if(vtxTyp.Contains("vertexer:Z") && (zRes>0.25)) return;
            if(TMath::Abs(spdVtx->GetZ() - trkVtx->GetZ())>0.5) return;
            if(TMath::Abs(zvtx) > 10) return;
        }
        
        
        
        fNevent2->Fill(12);
        
        //check pA pileup cut as it is done in the hfe package
        if(fAnalysisUtils->IsPileUpEvent(fVevent)){
            fNevent2->Fill(14);
            return;
        }
        fNevent2->Fill(16);
        
    }//close !Ispp flag
    
    //========================== vertex selection for pp
    
    if(fIspp)
    {
        if(fIsAOD)
        {
            const AliAODVertex* trkVtx = fAOD->GetPrimaryVertex();
            if(!trkVtx || trkVtx->GetNContributors()<=0) return;
            Float_t zvtx = -100;
            zvtx=trkVtx->GetZ();
            fZvtx = zvtx;
            //fVtxZ_new1->Fill(fZvtx);
            
            if(TMath::Abs(zvtx) > 10) return;
            //fVtxZ_new2->Fill(fZvtx);
            
            //Look for kink mother for AOD
            
            fNumberOfVertices = 0;
            fNumberOfMotherkink = 0;
            
            if(fIsAOD)
            {
                fNumberOfVertices = fAOD->GetNumberOfVertices();
                
                fListOfmotherkink = new Double_t[fNumberOfVertices];
                
                for(Int_t ivertex=0; ivertex < fNumberOfVertices; ivertex++)
                {
                    AliAODVertex *aodvertex = fAOD->GetVertex(ivertex);
                    if(!aodvertex) continue;
                    if(aodvertex->GetType()==AliAODVertex::kKink)
                    {
                        AliAODTrack *mother1 = (AliAODTrack *) aodvertex->GetParent();
                        if(!mother1) continue;
                        Int_t idmother = mother1->GetID();
                        fListOfmotherkink[fNumberOfMotherkink] = idmother;
                        fNumberOfMotherkink++;
                    }
                }
            }
        }
        else
        {
            
            
            
            /// ESD
            const AliESDVertex* trkVtx = fESD->GetPrimaryVertex();
            if(!trkVtx || trkVtx->GetNContributors()<=0) return;
            Float_t zvtx = -100;
            zvtx=trkVtx->GetZ();
            if(TMath::Abs(zvtx) > 10) return;
        }
    }
    //______________________________________________________________________
    //after vertex selection
    fNevent->Fill(10);
    
    //Only events with at least 2 tracks are accepted
    Int_t fNOtrks =  fVevent->GetNumberOfTracks();
    if(fNOtrks<2) return;
    
    fNevent->Fill(11);
    
    if(fIsAOD){
        Int_t fNOtrks2 =  fAOD->GetNumberOfTracks();
        if(fNOtrks2<2) return;
    }
    fNevent->Fill(12);
    
    
    //______________________________________________________________________
    //______________________________________________________________________
    //______________________________________________________________________
    //Centrality Selection
    if(!fIspp){
        if(fHasCentralitySelection)
        {
            Float_t centrality = -1;
            
            if(fIsAOD)
            {
                //fCentrality = fAOD->GetHeader()->GetCentralityP();
                fCentrality = ((AliAODHeader*)fAOD->GetHeader())->GetCentralityP();
            }
            else
            {
                fCentrality = fESD->GetCentrality();
            }
            
            if(fEstimator==1) centrality = fCentrality->GetCentralityPercentile("ZEMvsZDC");
            else centrality = fCentrality->GetCentralityPercentile("V0A");
            
            //printf("Centrality ZDC= %f VOA= %f  ZNA = %f ZNC = %f  ZPA = %f TRK = %f \n", fCentrality->GetCentralityPercentile("ZEMvsZDC") , fCentrality->GetCentralityPercentile("V0A"),fCentrality->GetCentralityPercentile("ZNA"), fCentrality->GetCentralityPercentile("ZNC"), fCentrality->GetCentralityPercentile("ZPA"), fCentrality->GetCentralityPercentile("TRK") );
            
            fCentralityHist->Fill(centrality);
            
            if(centrality<fCentralityMin || centrality>fCentralityMax) return;
            
            fCentralityHistPass->Fill(centrality);
        }
    }
    //______________________________________________________________________

    
    //______________________________________________________________________
    //threshold selection was here
    //______________________________________________________________________
    //all events selected
    
    fNevent->Fill(0);
    
    
    //______________________________________________________________________
    //events in the threshold
    
    fVtxZ_new4->Fill(fZvtx);
    
    

    
    if (fIsMC)
    {
        //Inicialzing MC Array
        fMCarray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
        fMCheader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
        
        Bool_t isNHFe = kFALSE;
        Bool_t isHFe = kFALSE;
        Bool_t isFromOtherLightMeson = kFALSE;
        Bool_t isOther = kFALSE;
        Bool_t HasNoMother = kFALSE;
        Bool_t IsOnDCAcut = kFALSE;
        
        //Vertex from MC
        Double_t vtxMC[3];
        fMCheader->GetVertex(vtxMC);
        
        for(Int_t iMC = 0; iMC < fMCarray->GetEntriesFast(); iMC++)
        {
            isNHFe = kFALSE;
            isHFe = kFALSE;
            isOther = kFALSE;
            HasNoMother = kFALSE;
            isFromOtherLightMeson = kFALSE;
            IsOnDCAcut = kFALSE;
            
            
            fMCparticle = (AliAODMCParticle*) fMCarray->At(iMC);
            
            Double_t EtaMC = fMCparticle->Eta();
            if (fMCparticle->Charge() == 0) continue;
            
            //Save the pT of all Charged hadrons in the acceptance (This is the denominator of the efficiency)
            if (fMCparticle->IsPhysicalPrimary())
                fHadronsGenerated->Fill(fMCparticle->Pt(),EtaMC, fZvtx);
            
            if(fMCparticle->Pt() < fMinpTElec || fMCparticle->Pt() > fMaxpTElec) continue;
            
            Int_t ParticlePDG = TMath::Abs(fMCparticle->GetPdgCode());
            Int_t MotherPDG = -999;
            
            if (ParticlePDG == 11)
            {
                if (fMCparticle->IsPhysicalPrimary())
                {
                    if(fMCparticle->GetMother()>0)
                    {
                        fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
                        MotherPDG = TMath::Abs(fMCparticleMother->GetPdgCode());
                        
                        Int_t MotherPDGHeavy  = Int_t (MotherPDG / TMath::Power(10, Int_t(TMath::Log10(MotherPDG))));
                        
                        if (MotherPDGHeavy > 3)
                        {
                            fNoEtaCutElectronGeneratedSignalPtEtaZvtx->Fill(fMCparticle->Pt());
                            if (EtaMC >= fEtaCutMin && EtaMC <= fEtaCutMax)
                                fEtaCutElectronGeneratedSignalPtEtaZvtx->Fill(fMCparticle->Pt());

                        }
                    }

                    if (!fCorrelationFlag) continue;
                   
                        for (Int_t iHadron = 0 ; iHadron < fMCarray->GetEntriesFast(); iHadron++)
                        {
                            
                            fMCparticle2 = (AliAODMCParticle*) fMCarray->At(iHadron);
                            
                            if (iHadron == iMC) continue;
                            
                            if (!fMCparticle2->IsPhysicalPrimary()) continue;
                            
                            if (fMCparticle2->Charge() == 0) continue;
                            
                            if(fMCparticle2->Pt() <fAssHadronPtMin || fMCparticle2->Pt() > fAssHadronPtMax) continue;
                            Double_t EtaHadron = fMCparticle2->Eta();
                            Double_t dEta =  EtaMC - EtaHadron;
                            Double_t dPhi = fMCparticle->Phi() - fMCparticle2->Phi();
                            
                            if (dPhi > 3*TMath::Pi()/2) dPhi = dPhi - 2*TMath::Pi();
                            if (dPhi < -TMath::Pi()/2)  dPhi = dPhi + 2*TMath::Pi();
                            
                            if (isNHFe || isOther)
                                if (fMCparticle2->GetMother() == fMCparticle->GetMother())
                                    continue;
                            
                            for (Int_t i = 0 ; i < fpTBins.GetSize()-1; i++)
                            {
                                if(fMCparticle->Pt()>=fpTBins.At(i) && fMCparticle->Pt()<fpTBins.At(i+1))
                                {
                                    fCEtaPhiNoEtaCutHFe[i]->Fill(dPhi,dEta);
                                    
                                    if ( EtaMC>= fEtaCutMin &&  EtaMC <= fEtaCutMax && EtaHadron >= fEtaCutMin &&  EtaHadron <= fEtaCutMax )
                                    {
                                        fCEtaPhiCutHFe[i]->Fill(dPhi,dEta);

                                    }
                                    
                                }
                            }
                            
                    }
                    
                    
                    
                    
                }
                
            }
            
            
        
        }
    }
    
    //______________________________________________________________________
    
    ///_____________________________________________________________________
    ///Track loop
    Int_t NTracks=0;
    
    NTracks=fVevent->GetNumberOfTracks();
    
    for(Int_t iTracks = 0; iTracks < NTracks; iTracks++)
    {
        //AliVParticle* Vtrack = fVevent->GetTrack(iTracks);
        AliVParticle* Vtrack = 0x0;
        Vtrack  = fVevent->GetTrack(iTracks);
        
        
        if (!Vtrack)
        {
            printf("ERROR: Could not receive track %d\n", iTracks);
            continue;
        }
        
        
        
        AliVTrack *track = dynamic_cast<AliVTrack*>(Vtrack);
        AliESDtrack *etrack = dynamic_cast<AliESDtrack*>(Vtrack);
        AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(Vtrack);
        
        
        Double_t fTPCnSigma0 = -999;
        Double_t fTPCnSigma = -999;
        Double_t fTOFnSigma = -999;
        Double_t fTPCsignal = -999;
        Double_t fPt = -999;
        Double_t fP = -999;
        Double_t EtaTrig = track->Eta();
        
        
        if (fIsMC)
        {
            Bool_t IsOnHadronDCAcut = kTRUE;
            if(EtaTrig>=fEtaCutMin && EtaTrig<=fEtaCutMax && atrack->TestFilterMask(AliAODTrack::kTrkTPCOnly) && atrack->GetStatus()&AliESDtrack::kITSrefit && atrack->GetStatus()&AliESDtrack::kTPCrefit && atrack->GetTPCNcls() >= 80)
            {
                
                if(fIsAOD && fUseDCACutforHadrons)
                {
                    Double_t d0z0[2], cov[3];
                    AliAODVertex *prim_vtx = fAOD->GetPrimaryVertex();
                    track->PropagateToDCA(prim_vtx, fAOD->GetMagneticField(), 20., d0z0, cov);
                    Double_t DCAxy = d0z0[0];
                    Double_t DCAz = d0z0[1];
                    if(TMath::Abs(DCAxy) >= fDCAcutrHadron || TMath::Abs(DCAz)>=fDCAcutzHadron)
                        IsOnHadronDCAcut = kFALSE;
                }
                
                if (IsOnHadronDCAcut)
                {
                    if (track->GetLabel() >= 0)
                    {
                        fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
                        if (fMCparticle->IsPhysicalPrimary())
                            fHadronsRecoPP->Fill(fMCparticle->Pt(),fMCparticle->Eta(), fZvtx);
                        
                    }
                    
                    fHadronsReco->Fill(track->Pt(),EtaTrig, fZvtx);
                }
            }
        }
        

        
        //Etacut test on the begging
        fEtad[0]->Fill(EtaTrig);
        
        ///_____________________________________________________________________________
        ///Fill QA plots without track selection
        fPt = track->Pt();
        fP = TMath::Sqrt((track->Pt())*(track->Pt()) + (track->Pz())*(track->Pz()));
        
        fTPCsignal = track->GetTPCsignal();
        fTPCnSigma = fPidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
        fTOFnSigma = fPidResponse->NumberOfSigmasTOF(track, AliPID::kElectron);
        
        
        fTPC_momentum->Fill(fP,fTPCsignal);
        fTPC_eta->Fill(EtaTrig,fTPCsignal);
        
    
        fTPC_p[0]->Fill(fPt,fTPCsignal);
        fTPCnsigma_p[0]->Fill(fP,fTPCnSigma);
        Float_t TPCNcls = track->GetTPCNcls();
        //TPC Ncls for pid
        Float_t TPCNcls_pid = track->GetTPCsignalN();
        //______________________________________________________________
        // Vertex
        
        fVtxZ[0]->Fill(fZvtx);
        if(iTracks == 0)fNTracks[0]->Fill(NTracks);
        fNTracks_pt[0]->Fill(NTracks, fPt);
        fNTracks_eta[0]->Fill(NTracks, EtaTrig);
        fNTracks_phi[0]->Fill(NTracks, track->Phi());
        
        
        fTPCNcls_pid[0]->Fill(TPCNcls, TPCNcls_pid);
        //______________________________________________________________
        
        ///Fill QA plots without track selection
        ///_________________________________________________________________________
        //__________________________________________________________________________
        //Track Selection Cuts
        
        if(EtaTrig<fEtaCutMin || EtaTrig>fEtaCutMax) continue;

        
        //AOD (Test Filter Bit)
        if(fIsAOD)
        {
            // standard cuts with very loose DCA - BIT(4)
            // Description:
            /*
             GetStandardITSTPCTrackCuts2011(kFALSE)
             SetMaxChi2PerClusterTPC(4);
             SetAcceptKinkDaughters(kFALSE);
             SetRequireTPCRefit(kTRUE);
             SetRequireITSRefit(kTRUE);
             SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
             SetMaxDCAToVertexZ(2);
             SetMaxDCAToVertex2D(kFALSE);
             SetRequireSigmaToVertex(kFALSE);
             SetMaxChi2PerClusterITS(36);
             SetMaxDCAToVertexXY(2.4)
             SetMaxDCAToVertexZ(3.2)
             SetDCaToVertex2D(kTRUE)
             */
            
            if(!atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue;
        }
        
        
        if(track->Pt()< fMinpTElec || EtaTrig> fMaxpTElec) continue;
        
        
        //RecKine: ITSTPC cuts
        if(!ProcessCutStep(AliHFEcuts::kStepRecKineITSTPC, track)) continue;
        
        //RecKink
        if(fRejectKinkMother)
        {
            if(fIsAOD)
            {
                Bool_t kinkmotherpass = kTRUE;
                for(Int_t kinkmother = 0; kinkmother < fNumberOfMotherkink; kinkmother++)
                {
                    if(track->GetID() == fListOfmotherkink[kinkmother])
                    {
                        kinkmotherpass = kFALSE;
                        continue;
                    }
                }
                if(!kinkmotherpass) continue;
            }
            else
            {
                if(etrack->GetKinkIndex(0) != 0) continue;
            }
        }
        
        
        if(!ProcessCutStep(AliHFEcuts::kStepRecPrim, track)) continue;
        
        //HFEcuts: ITS layers cuts
        if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsITS, track)) continue;
        
        //HFE cuts: TPC PID cleanup
        if(!ProcessCutStep(AliHFEcuts::kStepHFEcutsTPC, track)) continue;
        
        fEtad[1]->Fill(EtaTrig);
        fTPC_p[1]->Fill(fPt,fTPCsignal);
        fTPCnsigma_p[1]->Fill(fP,fTPCnSigma);
        TPCNcls = track->GetTPCNcls();
        Float_t pos2[3]={0,0,0};
        
        //______________________________________________________________
        // Vertex
        
        fVtxZ[1]->Fill(fZvtx);
        if(iTracks == 0)fNTracks[1]->Fill(NTracks);
        fNTracks_pt[1]->Fill(NTracks, fPt);
        fNTracks_eta[1]->Fill(NTracks, EtaTrig);
        fNTracks_phi[1]->Fill(NTracks, track->Phi());
        fTPCNcls_pid[1]->Fill(TPCNcls, TPCNcls_pid);
        
        //______________________________________________________________
        
        ///______________________________________________________________________
        ///Histograms for PID Studies
        
        for(Int_t i = 0; i < fpTBins.GetSize()-1; i++)
        {
            if(fPt>=fpTBins.At(i) && fPt<fpTBins.At(i+1))
            {
                fTPC_pt[i]->Fill(fTPCsignal);
                fTPCnsigma_pt[i]->Fill(fTPCnSigma);
                fTOFTPCnsigma_pt[i]->Fill(fTOFnSigma,fTPCnSigma);
                
            }
        }
        
        ///QA plots after track selection
        ///_____________________________________________________________
        
        //_______________________________________________________
        //Correlation Analysis - DiHadron
        if(fTPCnSigma < 3.5 && fCorrelationFlag)
        {
            fPtTrigger_Inc->Fill(fPt);
            DiHadronCorrelation(track, iTracks);
        }
        
        
        //if(fPt>1 && fPt<2) fTOF01->Fill(fTOFnSigma,fTPCnSigma);
        //if(fPt>2 && fPt<4) fTOF02->Fill(fTOFnSigma,fTPCnSigma);
        //if(fPt>4 && fPt<6) fTOF03->Fill(fTOFnSigma,fTPCnSigma);
        
        ///________________________________________________________________________
        ///PID
        ///Here the PID cuts defined in the file "ConfigEMCalHFEpA.C" is applied
        
        
        Int_t pidpassed = 1;
        AliHFEpidObject hfetrack;
        hfetrack.SetAnalysisType(AliHFEpidObject::kESDanalysis);
        hfetrack.SetRecTrack(track);
        hfetrack.SetPP();	//proton-proton analysis
        if(!fPID->IsSelected(&hfetrack, NULL, "", fPIDqa)) pidpassed = 0;
        fpid->Fill(pidpassed);
        
        if(pidpassed==0) continue;
        
        
        ///________________________________________________________________________
    
        
        //______________________________________________________________
        // Vertex
        
        fVtxZ[2]->Fill(fZvtx);
        if(iTracks == 0)fNTracks[2]->Fill(NTracks);
        fNTracks_pt[2]->Fill(NTracks, fPt);
        fNTracks_eta[2]->Fill(NTracks, EtaTrig);
        fNTracks_phi[2]->Fill(NTracks, track->Phi());
        fTPCNcls_pid[2]->Fill(TPCNcls, TPCNcls_pid);
        
        //______________________________________________________________
        
        //_______________________________________________________
        //Correlation Analysis
        
        
        
        fPtElec_Inc->Fill(fPt);
        
        Double_t d0z0[2], cov[3];
        
        AliAODVertex *prim_vtx = fAOD->GetPrimaryVertex();
        Double_t DCAxy = 100;
        Double_t DCAz = 100;
        if (track->PropagateToDCA(prim_vtx, fAOD->GetMagneticField(), 20., d0z0, cov))
        {
            DCAxy = TMath::Abs(d0z0[0]);
            DCAz =  TMath::Abs(d0z0[1]);
        }
        
        fDCAElectronXY[0]->Fill(DCAxy);
        fDCAElectronZ[0]->Fill(DCAz);

        

        ElectronHadronCorrelation(track, iTracks, Vtrack);
        
        //_______________________________________________________
        
        ///________________________________________________________________________
    }
    
    //__________________________________________________________________
    //Event Mixing Analysis
    //Filling pool
    if(fEventMixingFlag && fCorrelationFlag)
    {
        if(fIspp)
        {
            fPool = fPoolMgr->GetEventPool(1.5, fZvtx); // Get the buffer associated with the current centrality and z-vtx
            
            //if(!fPool) AliFatal(Form("No pool found for centrality = %f, zVtx = %f", fCentrality->GetCentralityPercentile("V0A"), fZvtx));
            if(!fPool) AliFatal(Form("No pool found for centrality = %f, zVtx = %f", 1.5, fZvtx));
            
            TObjArray *fArrayTracksMixed = SelectedHadrons();
            fArrayTracksMixed->SetOwner(kTRUE);
            
            fPool->UpdatePool(fArrayTracksMixed);
        }
        else
        {
            fPool = fPoolMgr->GetEventPool(fCentrality->GetCentralityPercentile("V0A"), fZvtx); // Get the buffer associated with the current centrality and z-vtx
            
            if(!fPool) AliFatal(Form("No pool found for centrality = %f, zVtx = %f", fCentrality->GetCentralityPercentile("V0A"), fZvtx));
            
            fPool->UpdatePool(SelectedHadrons()); // fill the tracks in the event buffer. The ownership is handed over to the event mixing class. We are not allowed to delete tracksClone anymore!
        }
        
    }
    
    //__________________________________________________________________
    
    delete fListOfmotherkink;
    PostData(1, fOutputList);
}
    

//______________________________________________________________________
void AliAnalysisTaskHFEpACorrelation::Terminate(Option_t *)
{
    //Draw result to the screen
    //Called once at the end of the query
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    
    if(!fOutputList)
    {
        printf("ERROR: Output list not available\n");
        return;
    }
}

//______________________________________________________________________
Bool_t AliAnalysisTaskHFEpACorrelation::ProcessCutStep(Int_t cutStep, AliVParticle *track)
{
    //Check single track cuts for a given cut step
    //Note this function is called inside the UserExec function
    const Int_t kMCOffset = AliHFEcuts::kNcutStepsMCTrack;
    if(!fCFM->CheckParticleCuts(cutStep + kMCOffset, track)) return kFALSE;
    return kTRUE;
}




//______________________________________________________________________

//______________________________________________________________________
void AliAnalysisTaskHFEpACorrelation::ElectronHadronCorrelation(AliVTrack *track, Int_t trackIndex, AliVParticle *vtrack)
{
    
    ///_________________________________________________________________
    ///MC analysis
    Bool_t lIsNHFe = kFALSE;
    Bool_t lIsHFe = kFALSE;
    Bool_t lIsOther = kFALSE;
    Bool_t lHasMother = kFALSE;;
    
    ///_________________________________________________________________
    
    //________________________________________________
    //Associated particle cut
    fPartnerCuts->SetAcceptKinkDaughters(kFALSE);
    fPartnerCuts->SetRequireITSRefit(kTRUE);
    fPartnerCuts->SetRequireTPCRefit(kTRUE);
    fPartnerCuts->SetEtaRange(fEtaCutMin,fEtaCutMax);
    fPartnerCuts->SetMaxChi2PerClusterTPC(4.0);
    fPartnerCuts->SetMinNClustersTPC(80);
    fPartnerCuts->SetPtRange(fPtMinAsso,1e10);
    //fPartnerCuts->SetRequireSigmaToVertex(kTRUE);
    //fPartnerCuts->SetMaxDCAToVertexXY(1);
    //fPartnerCuts->SetMaxDCAToVertexZ(3);
    //_________________________________________________
    
    ///#################################################################
    //Non-HFE reconstruction
    //fNonHFE = new AliSelectNonHFE(); I do not need to create if again
    fNonHFE->SetAODanalysis(fIsAOD);
    
    if(fMassCutFlag) fNonHFE->SetInvariantMassCut(fMassCut);
    if(fAngleCutFlag) fNonHFE->SetOpeningAngleCut(fAngleCut);
    if(fChi2CutFlag) fNonHFE->SetChi2OverNDFCut(fChi2Cut);
    if(fDCAcutFlag) fNonHFE->SetDCACut(fDCAcut);
    fNonHFE->SetAlgorithm("DCA"); //KF
    
    fNonHFE->SetPIDresponse(fPidResponse);
    fNonHFE->SetTrackCuts(-3.0,3.0,fPartnerCuts);
    fNonHFE->SetAdditionalCuts(fPtMinAsso,fTpcNclsAsso);
    
    
    fNonHFE->SetHistAngleBack(fOpAngleBack);
    fNonHFE->SetHistAngle(fOpAngle);
    fNonHFE->SetHistDCABack(fDCABack);
    fNonHFE->SetHistDCA(fDCA);
    fNonHFE->SetHistMassBack(fInvMassBack);
    fNonHFE->SetHistMass(fInvMass);
    
    fNonHFE->FindNonHFE(trackIndex,vtrack,fVevent);
    
    Int_t *fUlsPartner = fNonHFE->GetPartnersULS();
    Int_t *fLsPartner = fNonHFE->GetPartnersLS();
    Bool_t fUlsIsPartner = kFALSE;
    Bool_t fLsIsPartner = kFALSE;
    ///#################################################################
    
    
    //Electron Information
    Double_t fPhiE = -999;
    Double_t fEtaE = -999;
    Double_t fPhiH = -999;
    Double_t fEtaH = -999;
    Double_t fDphi = -999;
    Double_t fDeta = -999;
    Double_t fPtE = -999;
    Double_t fPtH = -999;
    
    Double_t pi = TMath::Pi();
    
    fPhiE = track->Phi();
    fEtaE = track->Eta();
    fPtE = track->Pt();
    
    
    if(fIsMC)
    {
        if(track->GetLabel() < 0)
        {
            fElectronNoLabel->Fill(track->Pt());
            if (fNonHFE->IsULS(),fNonHFE->GetNULS())
                fElectronNoLabelULS->Fill(track->Pt());
            if (fNonHFE->IsLS())
                fElectronNoLabelLS->Fill(track->Pt(),fNonHFE->GetNLS());
        }
        else
        {
            if (fNonHFE->IsULS())
                fEtaCutElectronBKWithLabelULS->Fill(track->Pt(),fNonHFE->GetNULS());
            if (fNonHFE->IsLS())
                fEtaCutElectronBKWithLabelLS->Fill(track->Pt(),fNonHFE->GetNLS());
            
            if(fIsAOD)
            {
                fMCparticle = (AliAODMCParticle*) fMCarray->At(track->GetLabel());
                
                Int_t ElectronPDG = TMath::Abs(fMCparticle->GetPdgCode());
                
                if (ElectronPDG == 11)
                {
                    fEtaCutElectronInclusiveRecoPtEtaZvtx->Fill(track->Pt());
                        
                    if(fMCparticle->GetMother()>=0)
                    {
                        lHasMother = kTRUE;
                        fMCparticleMother = (AliAODMCParticle*) fMCarray->At(fMCparticle->GetMother());
                        
                        Int_t MotherPDGAfterReco = TMath::Abs(fMCparticleMother->GetPdgCode());
                        Int_t MotherPDGHeavy  = Int_t (MotherPDGAfterReco / TMath::Power(10, Int_t(TMath::Log10(MotherPDGAfterReco))));
                        
                        //NHFE
                        if( MotherPDGAfterReco==22 || MotherPDGAfterReco ==111 || MotherPDGAfterReco ==221)
                        {
                            lIsNHFe = kTRUE;
                            fEtaCutElectronBKNoTag->Fill(track->Pt());
                            
                            if (fNonHFE->IsULS())
                                fEtaCutElectronBKULSMainSources->Fill(track->Pt(),fNonHFE->GetNULS());
                            if (fNonHFE->IsLS())
                                fEtaCutElectronBKLSMainSources->Fill(track->Pt(),fNonHFE->GetNLS());
                            
                            
                        }
                        else if(MotherPDGHeavy<4) //NHFE
                        {
                            fEtaCutElectronRecoOtherMC->Fill(track->Pt());
                            lIsOther = kTRUE;
                            
                            if (fNonHFE->IsULS())
                                fEtaCutElectronBKULSOtherSources->Fill(track->Pt(),fNonHFE->GetNULS());
                            if (fNonHFE->IsLS())
                                fEtaCutElectronBKLSOtherSources->Fill(track->Pt(),fNonHFE->GetNLS());

                            
                        }
                        else
                        {
                            fEtaCutElectronRecoHFEMC->Fill(track->Pt());
                            lIsHFe = kTRUE;
                            if (fNonHFE->IsULS())
                                fEtaCutElectronHFEULS->Fill(track->Pt(),fNonHFE->GetNULS());
                            if (fNonHFE->IsLS())
                                fEtaCutElectronHFELS->Fill(track->Pt(),fNonHFE->GetNLS());
                            
                        }
                        
                    }
                }
                else
                {
                    fMissIDElectronsReco->Fill(track->Pt());
                    if (fNonHFE->IsULS())
                        fEtaCutElectronMissIDULS->Fill(track->Pt(),fNonHFE->GetNULS());
                    if (fNonHFE->IsLS())
                        fEtaCutElectronMissIDLS->Fill(track->Pt(),fNonHFE->GetNLS());
                    
                }
                
            }
        }
    }

    
    if(fNonHFE->IsULS()) fPtElec_ULS->Fill(fPtE,fNonHFE->GetNULS());
    if(fNonHFE->IsLS()) fPtElec_LS->Fill(fPtE,fNonHFE->GetNLS());

    
    if(fIsMC)
    {
        if(track->GetLabel()> 0 && fIsAOD &&  fMCparticle->GetMother()>=0 )
        {
            if(TMath::Abs(fMCparticle->GetPdgCode())==11 && (TMath::Abs(fMCparticleMother->GetPdgCode())==22 || TMath::Abs(fMCparticleMother->GetPdgCode())==111 || TMath::Abs(fMCparticleMother->GetPdgCode())==221))
            {
                if(fNonHFE->IsULS()) fPtElec_ULS_NoPid->Fill(fPtE,fNonHFE->GetNULS());
                if(fNonHFE->IsLS()) fPtElec_LS_NoPid->Fill(fPtE,fNonHFE->GetNLS());
                
            }
        }
    }
    
    if (!fCorrelationFlag) return;
    
    if(fEventMixingFlag)
    {
        //hadling pp in the same task
        if (fIspp)
        {
            fPool = fPoolMgr->GetEventPool(1.5, fZvtx);
            
            if(!fPool) AliFatal(Form("No pool found for centrality = %f, zVtx = %f",1.5, fZvtx));
        }
        else
        {
            fPool = fPoolMgr->GetEventPool(fCentrality->GetCentralityPercentile("V0A"), fZvtx); // Get the buffer associated with the current centrality and z-vtx
            
            if(!fPool) AliFatal(Form("No pool found for centrality = %f, zVtx = %f",fCentrality->GetCentralityPercentile("V0A"), fZvtx));
        }
        if(fPool->GetCurrentNEvents() >= 5) // start mixing when 5 events are in the buffer
        {
            fPoolNevents->Fill(fPool->GetCurrentNEvents());
            
            for (Int_t jMix = 0; jMix < fPool->GetCurrentNEvents(); jMix++)  // mix with each event in the buffer
            {
                TObjArray* bgTracks = fPool->GetEvent(jMix);
                
                for (Int_t kMix = 0; kMix < bgTracks->GetEntriesFast(); kMix++)  // mix with each track in the event
                {
                    const AliHFEHCParticle* MixedTrack(dynamic_cast<AliHFEHCParticle*>(bgTracks->At(kMix)));
                    if (NULL == MixedTrack) continue;
                    
                    fPhiH = MixedTrack->Phi();
                    fEtaH = MixedTrack->Eta();
                    fPtH = MixedTrack->Pt();
                    
                    if(fPtH<fAssHadronPtMin || fPtH>fAssHadronPtMax) continue;
                    
                    fDphi = fPhiE - fPhiH;
                    
                    if (fDphi > 3*pi/2) fDphi = fDphi - 2*pi;
                    if (fDphi < -pi/2)  fDphi = fDphi + 2*pi;
                    
                    fDeta = fEtaE - fEtaH;
                    
                    Double_t WeightHadron = 1./GetHadronEfficiency(fPtH,fEtaH,fZvtx);
                    
                    for(Int_t i = 0; i < fpTBins.GetSize()-1; i++)
                    {
                        if(fPtE>=fpTBins.At(i) && fPtE<fpTBins.At(i+1))
                        {
                            fCEtaPhi_Inc_EM[i]->Fill(fDphi,fDeta,WeightHadron);
                            
                            if(fNonHFE->IsULS()) fCEtaPhi_ULS_Weight_EM[i]->Fill(fDphi,fDeta,fNonHFE->GetNULS()*WeightHadron);
                            if(fNonHFE->IsLS()) fCEtaPhi_LS_Weight_EM[i]->Fill(fDphi,fDeta,fNonHFE->GetNLS()*WeightHadron);
                            
                            
                        }
                    }
                }
            }
        }
        
    }
    //__________________________________________________________________
    
    //__________________________________________________________________
    //Same Event Analysis - Hadron Loop
    for(Int_t iTracks = 0; iTracks < fVevent->GetNumberOfTracks(); iTracks++)
    {
        if(trackIndex==iTracks) continue;
        
        AliVParticle* Vtrack2 = fVevent->GetTrack(iTracks);
        if (!Vtrack2)
        {
            printf("ERROR: Could not receive track %d\n", iTracks);
            continue;
        }
        
        AliVTrack *track2 = dynamic_cast<AliVTrack*>(Vtrack2);
        
        fPhiH = track2->Phi();
        fEtaH = track2->Eta();
        fPtH = track2->Pt();
        
        if(fEtaH<fEtaCutMin || fEtaH>fEtaCutMax ) continue;
        
        if(fIsAOD)
        {
            AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(Vtrack2);
            if(!atrack2->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
            if((!(atrack2->GetStatus()&AliESDtrack::kITSrefit)|| (!(atrack2->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
            if(atrack2->GetTPCNcls() < 80) continue;
            if(fAssocWithSPD && ((!(atrack2->HasPointOnITSLayer(0))) && (!(atrack2->HasPointOnITSLayer(1))))) continue;
        }
        else
        {
            AliESDtrack *etrack2 = dynamic_cast<AliESDtrack*>(Vtrack2);
            if(!fPartnerCuts->AcceptTrack(etrack2)) continue;
        }
        
        //DCA cut for hadrons
        if(fIsAOD && fUseDCACutforHadrons)
        {
            Double_t d0z0[2], cov[3];
            AliAODVertex *prim_vtx = fAOD->GetPrimaryVertex();
            track2->PropagateToDCA(prim_vtx, fAOD->GetMagneticField(), 20., d0z0, cov);
            Double_t DCAxy = d0z0[0];
            Double_t DCAz = d0z0[1];
            if(TMath::Abs(DCAxy) >= fDCAcutrHadron || TMath::Abs(DCAz)>=fDCAcutzHadron) continue;
        }
        
        
        
        if(fPtH<fAssHadronPtMin || fPtH>fAssHadronPtMax) continue;
        
        fDphi = fPhiE - fPhiH;
        
        if (fDphi > 3*pi/2) fDphi = fDphi - 2*pi;
        if (fDphi < -pi/2)  fDphi = fDphi + 2*pi;
        
        fDeta = fEtaE - fEtaH;
        

        
        //______________________________________________________________
        //Check if this track is a Non-HFE partner
        fUlsIsPartner = kFALSE;
        fLsIsPartner = kFALSE;
        for(Int_t i = 0; i < fNonHFE->GetNULS(); i++)
        {
            if(fUlsPartner[i]==iTracks) fUlsIsPartner=kTRUE;
        }
        for(Int_t i = 0; i < fNonHFE->GetNLS(); i++)
        {
            if(fLsPartner[i]==iTracks) fLsIsPartner=kTRUE;
        }
        //______________________________________________________________
        
        //Double_t WeightInclusive = 1./(GetInclusiveEfficiency(fPtH,fEtaH,fZvtx) * GetHadronEfficiency(fPtH,fEtaH,fZvtx));
        //Double_t WeightBKG = 1./(GetBackgroundEfficiency(fPtH,fEtaH,fZvtx) * GetHadronEfficiency(fPtH,fEtaH,fZvtx));
        
        Double_t WeightHadron = 1./GetHadronEfficiency(fPtH,fEtaH,fZvtx);
        
        for(Int_t i = 0; i < fpTBins.GetSize()-1; i++)
        {
            if(fPtE>=fpTBins.At(i) && fPtE<fpTBins.At(i+1))
            {
                
                //Filling histograms: Only hadron weight for now.
                fCEtaPhi_Inc[i]->Fill(fDphi,fDeta,WeightHadron);
    
                if(fNonHFE->IsULS())
                    fCEtaPhi_ULS_Weight[i]->Fill(fDphi,fDeta,fNonHFE->GetNULS()*WeightHadron);
                
                if(fNonHFE->IsLS())
                    fCEtaPhi_LS_Weight[i]->Fill(fDphi,fDeta,fNonHFE->GetNLS()*WeightHadron);
                
                if(fNonHFE->IsULS() && !fUlsIsPartner && !fLsIsPartner)
                    fCEtaPhi_ULS_NoP_Weight[i]->Fill(fDphi,fDeta,fNonHFE->GetNULS()*WeightHadron);
                
                if(fNonHFE->IsLS() && !fUlsIsPartner && !fLsIsPartner)
                    fCEtaPhi_LS_NoP_Weight[i]->Fill(fDphi,fDeta,fNonHFE->GetNLS()*WeightHadron);
                
                
                if (fIsMC)
                {
                    if (lHasMother)
                    {
                        if(lIsNHFe)
                            fCEtaPhi_Back_MC_Tag[i]->Fill(fDphi,fDeta,WeightHadron);
                        if(lIsOther)
                            fCEtaPhi_Other_MC_Tag[i]->Fill(fDphi,fDeta,WeightHadron);
                        if(lIsHFe)
                            fCEtaPhi_HFe_MC_Tag[i]->Fill(fDphi,fDeta,WeightHadron);
                            
                    }
                    else
                        fCEtaPhi_MC_NoMother_Tag[i]->Fill(fDphi,fDeta,WeightHadron);
                    
                }
                
            }
        }
    }
    
    
}

//____________________________________________________________________________________________________________
//Create a TObjArray with selected hadrons, for the mixed event analysis
TObjArray* AliAnalysisTaskHFEpACorrelation::SelectedHadrons()
{
    fTracksClone = new TObjArray;
    fTracksClone->SetOwner(kTRUE);
    
    for(Int_t iTracks = 0; iTracks < fVevent->GetNumberOfTracks(); iTracks++)
    {
        AliVParticle* Vtrack2 = fVevent->GetTrack(iTracks);
        if (!Vtrack2)
        {
            printf("ERROR: Could not receive track %d\n", iTracks);
            continue;
        }
        
        
        AliVTrack *track2 = dynamic_cast<AliVTrack*>(Vtrack2);
        
        if(track2->Pt() < fAssHadronPtMin || track2->Pt() > fAssHadronPtMax) continue;
        
        Double_t HadronEta = track2->Eta();
        
        if(HadronEta<fEtaCutMin || HadronEta>fEtaCutMax ) continue;
        
        if(fIsAOD)
        {
            AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(Vtrack2);
            if(!atrack2->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
            if((!(atrack2->GetStatus()&AliESDtrack::kITSrefit)|| (!(atrack2->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
            if(atrack2->GetTPCNcls() < 80) continue;
            if(fAssocWithSPD && ((!(atrack2->HasPointOnITSLayer(0))) && (!(atrack2->HasPointOnITSLayer(1))))) continue;
        }
        else
        {
            AliESDtrack *etrack2 = dynamic_cast<AliESDtrack*>(Vtrack2);
            if(!fPartnerCuts->AcceptTrack(etrack2)) continue;
        }
        
        //DCA cut for hadrons
        if(fIsAOD && fUseDCACutforHadrons)
        {
            Double_t d0z0[2], cov[3];
            AliAODVertex *prim_vtx = fAOD->GetPrimaryVertex();
            track2->PropagateToDCA(prim_vtx, fAOD->GetMagneticField(), 20., d0z0, cov);
            Double_t DCAxy = d0z0[0];
            Double_t DCAz = d0z0[1];
            if(TMath::Abs(DCAxy) >= fDCAcutrHadron || TMath::Abs(DCAz)>=fDCAcutzHadron) continue;
        }
        
        fTracksClone->Add(new AliHFEHCParticle(HadronEta, track2->Phi(), track2->Pt(), 1./(GetHadronEfficiency(track2->Pt(),HadronEta, fZvtx)) ) );
    }
    return fTracksClone;
}
//____________________________________________________________________________________________________________

//______________________________________________________________________
void AliAnalysisTaskHFEpACorrelation::DiHadronCorrelation(AliVTrack *track, Int_t trackIndex)
{
    //________________________________________________
    //Associated particle cut
    fPartnerCuts->SetAcceptKinkDaughters(kFALSE);
    fPartnerCuts->SetRequireITSRefit(kTRUE);
    fPartnerCuts->SetRequireTPCRefit(kTRUE);
    fPartnerCuts->SetEtaRange(-0.9,0.9);
    fPartnerCuts->SetMaxChi2PerClusterTPC(4.0);
    fPartnerCuts->SetMinNClustersTPC(80);
    fPartnerCuts->SetPtRange(0.3,1e10);
    //fPartnerCuts->SetRequireSigmaToVertex(kTRUE);
    //fPartnerCuts->SetMaxDCAToVertexXY(1);
    //fPartnerCuts->SetMaxDCAToVertexZ(3);
    //_________________________________________________
    
    //Electron Information
    Double_t fPhiE = -999;
    Double_t fEtaE = -999;
    Double_t fPhiH = -999;
    Double_t fEtaH = -999;
    Double_t fDphi = -999;
    Double_t fDeta = -999;
    Double_t fPtE = -999;
    Double_t fPtH = -999;
    
    Double_t pi = TMath::Pi();
    
    fPhiE = track->Phi();
    fEtaE = track->Eta();
    fPtE = track->Pt();
    
    //__________________________________________________________________
    //Same Event Analysis - Hadron Loop
    for(Int_t iTracks = 0; iTracks < fVevent->GetNumberOfTracks(); iTracks++)
    {
        if(trackIndex==iTracks) continue;
        
        AliVParticle* Vtrack2 = fVevent->GetTrack(iTracks);
        if (!Vtrack2)
        {
            printf("ERROR: Could not receive track %d\n", iTracks);
            continue;
        }
        
        AliVTrack *track2 = dynamic_cast<AliVTrack*>(Vtrack2);
        
        fPhiH = track2->Phi();
        fEtaH = track2->Eta();
        fPtH = track2->Pt();
        
        if(fEtaH<fEtaCutMin || fEtaH>fEtaCutMax ) continue;
        
        
        
        if(fIsAOD)
        {
            AliAODTrack *atrack2 = dynamic_cast<AliAODTrack*>(Vtrack2);
            if(!atrack2->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
            if((!(atrack2->GetStatus()&AliESDtrack::kITSrefit)|| (!(atrack2->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
            if(atrack2->GetTPCNcls() < 80) continue;
            if(fAssocWithSPD && ((!(atrack2->HasPointOnITSLayer(0))) && (!(atrack2->HasPointOnITSLayer(1))))) continue;
        }
        else
        {
            AliESDtrack *etrack2 = dynamic_cast<AliESDtrack*>(Vtrack2);
            if(!fPartnerCuts->AcceptTrack(etrack2)) continue;
        }
        
        //DCA cut for hadrons
        if(fIsAOD && fUseDCACutforHadrons)
        {
            Double_t d0z0[2], cov[3];
            AliAODVertex *prim_vtx = fAOD->GetPrimaryVertex();
            track2->PropagateToDCA(prim_vtx, fAOD->GetMagneticField(), 20., d0z0, cov);
            Double_t DCAxy = d0z0[0];
            Double_t DCAz = d0z0[1];
            if(TMath::Abs(DCAxy) >= fDCAcutrHadron || TMath::Abs(DCAz)>=fDCAcutzHadron) continue;
        }
        
        
        
        
        if(fPtH<fAssHadronPtMin || fPtH>fAssHadronPtMax) continue;
        
        fDphi = fPhiE - fPhiH;
        
        if (fDphi > 3*pi/2) fDphi = fDphi - 2*pi;
        if (fDphi < -pi/2)  fDphi = fDphi + 2*pi;
        
        fDeta = fEtaE - fEtaH;
        
        
        for(Int_t i = 0; i<fpTBins.GetSize()-1; i++)
        {
            if(fPtE>=fpTBins.At(i) && fPtE<fpTBins.At(i+1))
            {
                fCEtaPhi_Inc_DiHadron[i]->Fill(fDphi,fDeta);
            }
        }
    }
}
//____________________________________________________________________________________________________________

//______________________________________________________________________


Double_t AliAnalysisTaskHFEpACorrelation::GetHadronEfficiency(Double_t pT, Double_t eta, Double_t zvtx)
{
    if (!fEffHadron)
        return 1.;
    Int_t bin = fEffHadron->FindBin(pT,eta,zvtx);
    if ( fEffHadron->IsBinUnderflow(bin) || fEffHadron->IsBinOverflow(bin) ) {return 1.;}
    //printf ("Setting hadron Eff as %f\n",fEffHadron->GetBinContent(bin));
    return fEffHadron->GetBinContent(bin);
    
}



