

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


//////////////////////////////////////////////
//    Task for Measurement of Heavy Flavour //
//    electron as a function of charged     //
//    particle multiplicity                 //
//    Author: Preeti Dhankher               //
//////////////////////////////////////////////


#include "TChain.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "THnSparse.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TFile.h"
#include "TVector3.h"
#include <TRandom3.h>
#include "TProfile.h"
#include "TGeoManager.h"

#include "stdio.h"
#include "iostream"
#include "fstream"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDHandler.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"
#include "AliAnalysisTaskHFEMultiplicity.h"
#include "TGeoGlobalMagField.h"
#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "TRefArray.h"
#include "TVector.h"
#include "AliEventPoolManager.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliESDpid.h"
#include "AliAODPid.h"

#include "AliMultSelection.h"
#include "AliMultiplicity.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliAODTracklets.h"
#include "AliESDCaloCluster.h"
#include "AliAODCaloCluster.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliESDCaloTrigger.h"
#include "AliEMCALGeometry.h"
#include "AliGeomManager.h"
#include "AliCentrality.h"
#include "AliMagF.h"
#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliVEvent.h"
#include "AliStack.h"
#include "AliESDVZERO.h"
#include "AliAODVZERO.h"
#include "AliESDUtils.h"
#include "AliVertexingHFUtils.h"
#include "AliNormalizationCounter.h"
#include "AliAnalysisVertexingHF.h"

#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliGenHijingEventHeader.h"
#include "AliGenPythiaEventHeader.h"

class AliAnalysisTaskHFEMultiplicity;
using namespace std;
ClassImp(AliAnalysisTaskHFEMultiplicity)

AliAnalysisTaskHFEMultiplicity::AliAnalysisTaskHFEMultiplicity() : AliAnalysisTaskSE(),

//Event Cut
fCutNcontV(2),
//Track Cut
fCutTPCMaxCls(100.),
fCutTPCchi2perNDF(4.),
fCutTPCNCls(80.),
fCutITSNCls(3.),
fCutDCAxy(2.4),
fCutDCAz(3.2),
fCutTrackEta(0.7),
fCutpTMin(1.),
//PID Cut
fCutEopEMin(0.8),
fCutEopEMax(1.2),
fCutNsigmaEMin(-1.),
fCutNsigmaEMax(3.),
fCutM20Min(0.02),
fCutM20Max(0.35),
//Loose cuts for photonic electron pair
fAssoTPCCluster(80.),
fAssoITSCluster(3.),
fCutAssoEPt(0.1),
fCutAssoEEta(0.9),
fCutAssoENsigma(3),
fAssoITSRefit(kTRUE),
//Mass Cut for photonic electron pair
fCutInvmass(0.14),


// events
fAOD(0),
fNevents(0),
fOutputList(0),
fListProfiles(0),
fpidResponse(0),
// emcal correction
fUseTender(kTRUE),
fTenderClusterName("caloClusters"),
fTenderTrackName("tracks"),
fTracks_tender(0),
fCaloClusters_tender(0),
// MC
fMCArray(0),
fMCHeader(0),
fMCparticle(0),
// flag for emcal dcal
fFlagClsTypeEMC(kTRUE),
fFlagClsTypeDCAL(kTRUE),
// trigger events selection
fEMCEG1(kFALSE),
fEMCEG2(kFALSE),
fDCalDG1(kFALSE),
fDCalDG2(kFALSE),
// zvtx
fVtxZ(0),
fVtxX(0),
fVtxY(0),
// multi. estimation
fRejectPUFromSPD(kTRUE),
fRefMult(61.26),
gRandom(new TRandom3(0)),
fSparseMulti(0),
fvalueMulti(0),
// TPC info and PID
fTPCdEdx(0x0),
fTPCdEdxvsEp(0x0),
fTPCnsigma(0x0),
fTPCnSigma(-999.0),
fTrkPt(0),
fTrketa(0),
fTrkphi(0),
// cluster info.
fClusPhi(0),
fClusEta(0),
fClusEtaPhi(0x0),
fClusE(0),
fClusT(0),
fSparseClusE(0),
fvalueCluE(0),
// track-cluster matching
fTrkMatchTrkPt(0),
fTrkMatchTrketa(0),
fTrkMatchTrkphi(0),
fTrkMatchClusetaphi(0x0),
fEMCTrkMatchcluster(0x0),
// electron info
fSparseElectron(0),
fvalueElectron(0),
// photonic electron info
fInvmassLS(0),
fInvmassULS(0),
fInvmassLSPt(0),
fInvmassULSPt(0),
fULSElecPt(0),
fLSElecPt(0),
fSparseLSElectron(0),
fSparseULSElectron(0),
fvaluePHElectron(0),
// MC info
fReadMC(kFALSE),
fzvtxcut(kFALSE),
fzvtxQA(kFALSE),
fNTotMCpart(0),
fNpureMC(0),
fNembMCpi0(0),
fNembMCeta(0),
// enhance weight cal.
ftype(-1),
fWeight(1),
fCalculateWeight(kFALSE),
fCalculateElectronEffi(kTRUE),
fSprsPi0EtaWeightCal(0),
fPi0Weight(0),
fEtaWeight(0),
// inclusive e-
fInclsElecPt(0),
fInclsElecPtAll(0),
fInclsElecPtReco(0),
fInclseEMCALElecPtReco(0),
fInclseDCALElecPtReco(0),
// hfe e-
fHFElecPtAll(0),
fHFElecPtReco_wtrkcuts(0),
fHFElecPtReco_wTPCPID(0),
fHFElecPtReco_wtrkCalocuts(0),
fHFElecPtReco_wTPCCaloPID(0),
//  nonhfe MC
fIsFrmEmbPi0(kFALSE),
fIsFrmEmbEta(kFALSE),
fMissingEmbEtaEleTrkPt(0),
fNonHFeTrkPt(0),
fNonHFeEmbAllTypeTrkPt(0),
fNonHFeEmbTrkPt(0),
fNonHFeEmbWeightTrkPt(0),
fPi0eEmbWeightTrkPt(0),
fEtaeEmbWeightTrkPt(0),
fNonHFewPIDTrkPt(0),
fNonHFeEmbwPIDTrkPt(0),
fNonHFeEmbWeightwPIDTrkPt(0),
fPi0eEmbWeightwPIDTrkPt(0),
fEtaeEmbWeightwPIDTrkPt(0),
// non hfe reconstructed with invariant mass
fRecoNonHFeTrkPt(0),
fRecoNonHFeEmbTrkPt(0),
fRecoNonHFeEmbWeightTrkPt(0),
fRecoPi0eEmbWeightTrkPt(0),
fRecoEtaeEmbWeightTrkPt(0),
//mult dependent

//nonhfe invariant mass
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
fRecoLSeTrkPt(0),
fRecoLSeEmbTrkPt(0),
fRecoLSeEmbWeightTrkPt(0),
fRecoPi0LSeEmbWeightTrkPt(0),
fRecoEtaLSeEmbWeightTrkPt(0),
fRecoULSeTrkPt(0),
fRecoULSeEmbTrkPt(0),
fRecoULSeEmbWeightTrkPt(0),
fRecoPi0ULSeEmbWeightTrkPt(0),
fRecoEtaULSeEmbWeightTrkPt(0)



{ fvalueElectron = new Double_t[9];
    fvaluePHElectron = new Double_t[3];
    fvalueCluE = new Double_t[3];
    fvalueMulti = new Double_t[6];
    for(Int_t i=0; i<2; i++) fMultEstimatorAvg[i]=0;
    
}

//_____________________________________________________________________________
AliAnalysisTaskHFEMultiplicity::AliAnalysisTaskHFEMultiplicity(const char* name) : AliAnalysisTaskSE(name),


//Event Cut
fCutNcontV(2),
//Track Cut
fCutTPCMaxCls(100.),
fCutTPCchi2perNDF(4.),
fCutTPCNCls(80.),
fCutITSNCls(3.),
fCutDCAxy(2.4),
fCutDCAz(3.2),
fCutTrackEta(0.7),
fCutpTMin(1),
//PID Cut
fCutEopEMin(0.8),
fCutEopEMax(1.2),
fCutNsigmaEMin(-1.),
fCutNsigmaEMax(3.),
fCutM20Min(0.02),
fCutM20Max(0.35),
//Loose cuts for photonic electron pair
fAssoTPCCluster(80.),
fAssoITSCluster(3.),
fCutAssoEPt(0.1),
fCutAssoEEta(0.9),
fCutAssoENsigma(0.3),
fAssoITSRefit(kTRUE),
//Mass Cut for photonic electron pair
fCutInvmass(0.14),

// events
fAOD(0),
fNevents(0),
fOutputList(0),
fListProfiles(0),
fpidResponse(0),
// emcal correction
fUseTender(kTRUE),
fTenderClusterName("caloClusters"),
fTenderTrackName("tracks"),
fTracks_tender(0),
fCaloClusters_tender(0),
// MC
fMCArray(0),
fMCHeader(0),
fMCparticle(0),
// flag for emcal dcal
fFlagClsTypeEMC(kTRUE),
fFlagClsTypeDCAL(kTRUE),
// trigger events selection
fEMCEG1(kFALSE),
fEMCEG2(kFALSE),
fDCalDG1(kFALSE),
fDCalDG2(kFALSE),
// zvtx
fVtxZ(0),
fVtxX(0),
fVtxY(0),
// multi. estimation
fRejectPUFromSPD(kTRUE),
fRefMult(61.26),
gRandom(new TRandom3(0)),
fSparseMulti(0),
fvalueMulti(0),
// TPC info and PID
fTPCdEdx(0x0),
fTPCdEdxvsEp(0x0),
fTPCnsigma(0x0),
fTPCnSigma(-999.0),
fTrkPt(0),
fTrketa(0),
fTrkphi(0),
// cluster info.
fClusPhi(0),
fClusEta(0),
fClusEtaPhi(0x0),
fClusE(0),
fClusT(0),
fSparseClusE(0),
fvalueCluE(0),
// track-cluster matching
fTrkMatchTrkPt(0),
fTrkMatchTrketa(0),
fTrkMatchTrkphi(0),
fTrkMatchClusetaphi(0x0),
fEMCTrkMatchcluster(0x0),
// electron info
fSparseElectron(0),
fvalueElectron(0),
// photonic electron info
fInvmassLS(0),
fInvmassULS(0),
fInvmassLSPt(0),
fInvmassULSPt(0),
fULSElecPt(0),
fLSElecPt(0),
fSparseLSElectron(0),
fSparseULSElectron(0),
fvaluePHElectron(0),
// MC info
fReadMC(kFALSE),
fzvtxcut(kFALSE),
fzvtxQA(kFALSE),
fNTotMCpart(0),
fNpureMC(0),
fNembMCpi0(0),
fNembMCeta(0),
// enhance weight cal.
ftype(-1),
fWeight(1),
fCalculateWeight(kFALSE),
fCalculateElectronEffi(kTRUE),
fSprsPi0EtaWeightCal(0),
fPi0Weight(0),
fEtaWeight(0),
// inclusive e-
fInclsElecPt(0),
fInclsElecPtAll(0),
fInclsElecPtReco(0),
fInclseEMCALElecPtReco(0),
fInclseDCALElecPtReco(0),
// hfe e-
fHFElecPtAll(0),
fHFElecPtReco_wtrkcuts(0),
fHFElecPtReco_wTPCPID(0),
fHFElecPtReco_wtrkCalocuts(0),
fHFElecPtReco_wTPCCaloPID(0),
//  nonhfe MC
fIsFrmEmbPi0(kFALSE),
fIsFrmEmbEta(kFALSE),
fMissingEmbEtaEleTrkPt(0),
fNonHFeTrkPt(0),
fNonHFeEmbAllTypeTrkPt(0),
fNonHFeEmbTrkPt(0),
fNonHFeEmbWeightTrkPt(0),
fPi0eEmbWeightTrkPt(0),
fEtaeEmbWeightTrkPt(0),
fNonHFewPIDTrkPt(0),
fNonHFeEmbwPIDTrkPt(0),
fNonHFeEmbWeightwPIDTrkPt(0),
fPi0eEmbWeightwPIDTrkPt(0),
fEtaeEmbWeightwPIDTrkPt(0),
// non hfe reconstructed with invariant mass
fRecoNonHFeTrkPt(0),
fRecoNonHFeEmbTrkPt(0),
fRecoNonHFeEmbWeightTrkPt(0),
fRecoPi0eEmbWeightTrkPt(0),
fRecoEtaeEmbWeightTrkPt(0),

//mult dependent

//nonhfe invariant mass
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
fRecoLSeTrkPt(0),
fRecoLSeEmbTrkPt(0),
fRecoLSeEmbWeightTrkPt(0),
fRecoPi0LSeEmbWeightTrkPt(0),
fRecoEtaLSeEmbWeightTrkPt(0),
fRecoULSeTrkPt(0),
fRecoULSeEmbTrkPt(0),
fRecoULSeEmbWeightTrkPt(0),
fRecoPi0ULSeEmbWeightTrkPt(0),
fRecoEtaULSeEmbWeightTrkPt(0)


{
    // constructor
    fvalueElectron = new Double_t[9];
    fvaluePHElectron = new Double_t[3];
    fvalueCluE = new Double_t[3];
    fvalueMulti = new Double_t[6];
    for(Int_t i=0; i<2; i++) fMultEstimatorAvg[i]=0;
    DefineInput(0, TChain::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
}


//_____________________________________________________________________________
AliAnalysisTaskHFEMultiplicity::~AliAnalysisTaskHFEMultiplicity()
{
    // destructor
    if(fOutputList) {
        delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
        delete fSparseClusE;
        delete []fvalueCluE;
        delete fSparseElectron;
        delete fSparseLSElectron;
        delete fSparseULSElectron;
        delete []fvalueElectron;
        delete []fvaluePHElectron;
        delete fSparseMulti;
        delete []fvalueMulti;
        delete fSprsPi0EtaWeightCal,
        delete fTracks_tender;
        delete fCaloClusters_tender;
    }
    
    for(Int_t i=0; i<2; i++) {
        if (fMultEstimatorAvg[i]) delete fMultEstimatorAvg[i];
    }
    delete fListProfiles;
    delete gRandom;
}
//_____________________________________________________________________________
void AliAnalysisTaskHFEMultiplicity::Init()
{
    // Initialization
    
    
    if(fDebug > 1) printf("AliAnalysisTaskHFEMultiplicity::Init() \n");
    
    
    fListProfiles = new TList();
    fListProfiles->SetOwner();
    TString period[2];
    Int_t nProfiles=2;
    period[0]="LHC16s";
    period[1]="LHC16r";
    
    
    for(Int_t i=0; i<nProfiles; i++){
        if(fMultEstimatorAvg[i]){
            TProfile* hprof=new TProfile(*fMultEstimatorAvg[i]);
            hprof->SetName(Form("ProfileTrkVsZvtx%s\n",period[i].Data()));
            fListProfiles->Add(hprof);
        }
    }
    
    PostData(2,fListProfiles);
    
    
    return;
}

//_____________________________________________________________________________

void AliAnalysisTaskHFEMultiplicity::UserCreateOutputObjects()
{
    
    printf("\nseed = %u\n\n", gRandom->GetSeed());
    
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
    
    
    fNevents = new TH1F ("fNevents","Number of events",6,-0.5,5.5);
    fNevents->GetYaxis()->SetTitle("counts");
    fNevents->GetXaxis()->SetBinLabel(1,"nEvents Total");
    fNevents->GetXaxis()->SetBinLabel(2,"nEvents With Zvtx cut");
    fNevents->GetXaxis()->SetBinLabel(3,"nEvents with pileup cut");
    fNevents->GetXaxis()->SetBinLabel(4,"nEvents with n contributor");
    fNevents->GetXaxis()->SetBinLabel(5,"nEvents with zvtxQA cut");
    fNevents->GetXaxis()->SetBinLabel(6,"nEvents with vertex cut");
    
    
    fClusPhi            = new TH1F("fClusPhi", "Cluster Phi distribution; #phi ; counts",100,0.,7);
    fClusEta            = new TH1F("fClusEta", "Cluster Eta distribution; #eta ; counts",50,-2,2);
    fClusEtaPhi        = new TH2F( "fClusEtaPhi","Cluster Eta Phi distribution; #eta ; #phi",100,-1.5,1.5,100,0.,7);
    fClusE           = new TH1F("fClusE","Cluster Energy ; Energy(GeV); counts",200,0,100);
    fClusT             = new TH1F( "fClusT","Cluster time distribution ; Time(ns) ; counts",500,-100,100);
    fVtxZ         = new TH1F("fVtxZ","Z vertex position;Vtx_{z};counts",1000,-50,50);
    fVtxY         = new TH1F("fVtxY","Y vertex position;Vtx_{y};counts",300,-15,15);
    fVtxX         = new TH1F("fVtxX","X vertex position;Vtx_{x};counts",300,-15,15);
    fTPCdEdx         = new TH2F("fTPCdEdx","All Track dE/dx distribution;p (GeV/c);dE/dx",200,0,20,500,0,160);
    fTPCdEdxvsEp         = new TH2F("fTPCdEdxEp","All Track dE/dx distribution;E/p;dE/dx",200,0,2,500,0,160);
    fTPCnsigma         = new TH2F("fTPCnsig","All Track TPC Nsigma distribution;p (GeV/c);#sigma_{TPC-dE/dx}",500,0,50,200,-10,10);
    fTrkPt         = new TH1F("fTrkPt","p_{T} distribution of all tracks;p_{T} (GeV/c);counts",500,0,50);
    fTrketa         = new TH1F("fTrketa","All Track #eta distribution;#eta;counts",100,-1.5,1.5);
    fTrkphi         = new TH1F("fTrkphi","All Track #phi distribution;#phi;counts",100,0,2*3.141);
    fTrkMatchTrkPt    = new TH1F("fTrkMatchTrkPt","p_{T} distribution of tracks with cluster;p_{T} (GeV/c);counts",500,0,50);
    fTrkMatchTrketa    = new TH1F("fTrkMatchTrketa","#eta distribution of tracks matched to Cluster;#eta;counts",100,-1.5,1.5);
    fTrkMatchTrkphi    = new TH1F("fTrkMatchTrkphi","#phi distribution of tracks matched to Cluster;#phi;counts",100,0,2*3.141);
    fTrkMatchClusetaphi    = new TH2F( "fTrkMatchClusetaphi","#eta#phi distribution of Clusters matched to tracks;#eta;#phi",100,-1.5,1.5,100,0.,7);
    fEMCTrkMatchcluster   = new TH2F("fEMCTrkMatchcluster","Distance of EMCAL cluster to its closest track Method 1",100,-0.3,0.3,100,-0.3,0.3);
    fInvmassLS         = new TH1F("fInvmassLS", "Inv mass of LS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 50,0,0.5);
    fInvmassULS         = new TH1F("fInvmassULS", "Inv mass of ULS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 50,0,0.5);
    fInvmassLSPt         = new TH2F("fInvmassLSPt", "Inv mass of LS (e,e) vs pT; p_{T}(GeV/c); mass(GeV/c^2); counts;", 500,0,50,50,0,0.5);
    fInvmassULSPt        = new TH2F("fInvmassULSPt", "Inv mass of ULS (e,e) vs pT; p_{T}(GeV/c); mass(GeV/c^2); counts;", 500,0,50,50,0,0.5);
    
    fULSElecPt          = new TH1F("fULSElecPt","p_{T} distribution of ULS electrons;p_{T} (GeV/c);counts",500,0,50);
    fLSElecPt         = new TH1F("fLSElecPt","p_{T} distribution of LS electrons;p_{T} (GeV/c);counts",500,0,50);
    
    
    
    
    
    //---------------THnSparse------------
    Int_t binsE[3]    =          {500, 1000, 350};
    Double_t xminE[3]    =    {  0,  0,   0};
    Double_t xmaxE[3]    =    {  50,   2000,   350};
    
    fSparseClusE     = new THnSparseD ("Cluster Energy","Cluster Energy;ClusterE;V0M;SPDTracklets;",3 ,binsE,xminE,xmaxE);
    
    
    
    Int_t bins[9]        =          {250,200,100,400,400, 1000, 350,500,300};
    Double_t xmin[9]    =    {  0, -10,0, 0,0,0, 0,0,-15};
    Double_t xmax[9]    =    {  50,10, 2,2,2, 2000, 350,50,15};
    
    fSparseElectron     = new THnSparseD ("Electron","Electron;pT;nsigma;E/P;M02;M20;V0M;SPDTracklets;EtrkMatch;Zvertex",9 ,bins,xmin,xmax);
    
    
    Int_t binsls[3]    =          {250, 1000, 350};
    Double_t xminls[3]    =    {  0, 0, 0};
    Double_t xmaxls[3]    =    {  50, 2000, 350};
    
    fSparseLSElectron     = new THnSparseD ("LSElectron","LSElectron;pT;V0M;SPDTracklets;",3 ,binsls,xminls,xmaxls);
    fSparseULSElectron     = new THnSparseD ("ULSElectron","ULSElectron;pT;V0M;SPDTracklets;",3 ,binsls,xminls,xmaxls);
    
    Int_t binsm[6]    =          {300,1000,350,1000,350,400};
    Double_t xminm[6]    =    {-15,0,0,0,0,0};
    Double_t xmaxm[6]    =    { 15,2000,350,2000,350,400};
    fSparseMulti         = new THnSparseD ("Multiplicity","Multiplicity;zvtx;V0M_data;SPDTracklets_data;Corrected_V0M;Corrected_SPDTracklets;ncharge;",6,binsm,xminm,xmaxm);
    
    
    
    if(fReadMC){
        
        fInclsElecPt         = new TH1F("fInclsElecPt","p_{T} distribution of inclusive electrons;p_{T} (GeV/c);counts",250,0,50);
        
        fMissingEmbEtaEleTrkPt = new TH1F("fMissingEmbEtaEleTrkPt","Missing electrons from embedded #eta  + No mom ;p_{T} (GeV/c);counts",250,0,50);
        
        fInclsElecPtAll     = new TH1F("fInclsElecPtAll","p_{T} distribution of all inclusive electrons;p_{T} (GeV/c);counts",250,0,50);
        
        fInclsElecPtReco     = new TH1F("fInclsElecPtReco","p_{T} distribution of reconstructed inclusive electrons;p_{T} (GeV/c);counts",250,0,50);
        
        fInclseEMCALElecPtReco = new TH1F("fInclseEMCALElecPtReco","p_{T} distribution of EMCAL reconstructed inclusive electrons;p_{T} (GeV/c);counts",250,0,50);
        
        fInclseDCALElecPtReco     = new TH1F("fInclseDCALElecPtReco","p_{T} distribution of DCAL reconstructed inclusive electrons;p_{T} (GeV/c);counts",250,0,50);
        
        fHFElecPtAll     = new TH2F("fHFElecPtAll","p_{T} distribution of all HFe ;p_{T} (GeV/c);counts;SPDTracklets",250,0,50,350,0,350);
        
        fHFElecPtReco_wtrkcuts = new TH2F("fHFElecPtReco_wtrkcuts","p_{T} distribution of HF electrons with track cuts;p_{T} (GeV/c);counts;SPDTracklets",250,0,50,350,0,350);
        
        fHFElecPtReco_wTPCPID = new TH2F("fHFElecPtReco_wTPCPID","p_{T} distribution of all HF electrons with TPC PID;p_{T} (GeV/c);counts;SPDTracklets",250,0,50,350,0,350);
        
        fHFElecPtReco_wtrkCalocuts = new TH2F("fHFElecPtReco_wtrkCalocuts","p_{T} distribution of all HF electrons with Track matching;p_{T} (GeV/c);counts;SPDTracklets",250,0,50,350,0,350);
        
        fHFElecPtReco_wTPCCaloPID = new TH2F("fHFElecPtReco_wTPCCaloPID","p_{T} distribution of all HF electrons with EMCal;p_{T} (GeV/c);counts;SPDTracklets",250,0,50,350,0,350);
        
        fNonHFeTrkPt = new TH1F("fNonHFeTrkPt","Non-HF electrons from all generators;p_{T} (GeV/c);counts",250,0,50);
        
        fNonHFeEmbAllTypeTrkPt = new TH1F("fNonHFeEmbAllTypeTrkPt","Non-HF electrons from embedded #pi^{0} and #eta of all type;p_{T} (GeV/c);counts",250,0,50);
        
        fNonHFeEmbTrkPt = new TH1F("fNonHFeEmbTrkPt","Non-HF electrons from embedded #pi^{0} and #eta + No mom;p_{T} (GeV/c);counts",250,0,50);
        
        fNonHFeEmbWeightTrkPt = new TH1F("fNonHFeEmbWeightTrkPt","Non-HF electrons from embedded #pi^{0} and #eta + No mom with weight + No mom;p_{T} (GeV/c);counts",250,0,50);
        
        fPi0eEmbWeightTrkPt = new TH1F("fPi0eEmbWeightTrkPt","Non-HF electrons from embedded #pi^{0} + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
        
        fEtaeEmbWeightTrkPt = new TH1F("fEtaeEmbWeightTrkPt","Non-HF electrons from embedded #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
        
        fNonHFewPIDTrkPt = new TH1F("fNonHFewPIDTrkPt","Non-HF electrons from all generators wPID;p_{T} (GeV/c);counts",250,0,50);
        
        fNonHFeEmbwPIDTrkPt = new TH1F("fNonHFeEmbwPIDTrkPt","Non-HF electrons from embedded #pi^{0} and #eta + No mom wPID;p_{T} (GeV/c);counts",250,0,50);
        
        fNonHFeEmbWeightwPIDTrkPt = new TH1F("fNonHFeEmbWeightwPIDTrkPt","Non-HF electrons from embedded #pi^{0} and #eta + No mom with weight + No mom wPID;p_{T} (GeV/c);counts",250,0,50);
        
        fPi0eEmbWeightwPIDTrkPt = new TH1F("fPi0eEmbWeightwPIDTrkPt","Non-HF electrons from embedded #pi^{0} + No mom with weight wPID;p_{T} (GeV/c);counts",250,0,50);
        
        fEtaeEmbWeightwPIDTrkPt = new TH1F("fEtaeEmbWeightwPIDTrkPt","Non-HF electrons from embedded #eta  + No mom with weight wPID;p_{T} (GeV/c);counts",250,0,50);
        
        fRecoNonHFeTrkPt = new TH1F("fRecoNonHFeTrkPt"," Reco Non-HF electrons from all generators;p_{T} (GeV/c);counts",250,0,50);
        
        fRecoNonHFeEmbTrkPt = new TH1F("fRecoNonHFeEmbTrkPt","Reco Non-HF electrons from embedded #pi^{0} and #eta + No mom;p_{T} (GeV/c);counts",250,0,50);
        
        fRecoNonHFeEmbWeightTrkPt = new TH1F("fRecoNonHFeEmbWeightTrkPt","Reco Non-HF electrons from embedded #pi^{0} and #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
        
        
        fRecoPi0eEmbWeightTrkPt = new TH1F("fRecoPi0eEmbWeightTrkPt","Reco Non-HF electrons from embedded #pi^{0}  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
        
        fRecoEtaeEmbWeightTrkPt = new TH1F("fRecoEtaeEmbWeightTrkPt","Reco Non-HF electrons from embedded #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
        
        fNonHFePairInvmassLS= new TH1F("fNonHFePairInvmassLS", "Inv mass of LS (e,e) if both e- are Non-HFE; mass(GeV/c^2); counts;",  50,0,0.5);
        
        fNonHFePairInvmassULS = new TH1F("fNonHFePairInvmassULS", "Inv mass of ULS (e,e) if both e- are Non-HFE; mass(GeV/c^2); counts;",  50,0,0.5);
        
        fNonHFeEmbInvmassLS= new TH1F("fNonHFeEmbInvmassLS", "Inv mass of LS (e,e) for Non-HFE from embedded #pi^{0} and #eta; mass(GeV/c^2); counts;",  50,0,0.5);
        
        fNonHFeEmbInvmassULS = new TH1F("fNonHFeEmbInvmassULS", "Inv mass of ULS (e,e) for Non-HFE from embedded #pi^{0} and #eta; mass(GeV/c^2); counts;",  50,0,0.5);
        
        fNonHFeEmbWeightInvmassLS= new TH1F("fNonHFeEmbWeightInvmassLS", "Inv mass of LS (e,e) for Non-HFE from embedded #pi^{0} and #eta with weight; mass(GeV/c^2); counts;",  50,0,0.5);
        
        fNonHFeEmbWeightInvmassULS = new TH1F("fNonHFeEmbWeightInvmassULS", "Inv mass of ULS (e,e) for Non-HFE from embedded #pi^{0} and #eta with weight; mass(GeV/c^2); counts;",  50,0,0.5);
        
        fPi0EmbInvmassLS= new TH1F("fPi0EmbInvmassLS", "Inv mass of LS (e,e) for ele from embedded #pi^{0}; mass(GeV/c^2); counts;",  50,0,0.5);
        
        fPi0EmbInvmassULS = new TH1F("fPi0EmbInvmassULS", "Inv mass of ULS (e,e) for ele from embedded #pi^{0}; mass(GeV/c^2); counts;",  50,0,0.5);
        
        fPi0EmbWeightInvmassLS = new TH1F("fPi0EmbWeightInvmassLS", "Inv mass of LS (e,e) for ele from embedded #pi^{0} with weight; mass(GeV/c^2); counts;",  50,0,0.5);
        
        fPi0EmbWeightInvmassULS = new TH1F("fPi0EmbWeightInvmassULS", "Inv mass of ULS (e,e) for ele from embedded #pi^{0} with weight; mass(GeV/c^2); counts;",  50,0,0.5);
        
        fEtaEmbInvmassLS = new TH1F("fEtaEmbInvmassLS", "Inv mass of LS (e,e) for ele from embedded #eta; mass(GeV/c^2); counts;",  50,0,0.5);
        
        fEtaEmbInvmassULS= new TH1F("fEtaEmbInvmassULS", "Inv mass of ULS (e,e) for ele from embedded #eta; mass(GeV/c^2); counts;",  50,0,0.5);
        
        fEtaEmbWeightInvmassLS = new TH1F("fEtaEmbWeightInvmassLS", "Inv mass of LS (e,e) for ele from embedded #eta with weight; mass(GeV/c^2); counts;",  50,0,0.5);
        
        fEtaEmbWeightInvmassULS = new TH1F("fEtaEmbWeightInvmassULS", "Inv mass of ULS (e,e) for ele from embedded #eta with weight; mass(GeV/c^2); counts;",  50,0,0.5);
        
        fRecoLSeTrkPt = new TH1F("fRecoLSeTrkPt"," Reco LS electrons from all generators;p_{T} (GeV/c);counts",250,0,50);
        
        fRecoLSeEmbTrkPt = new TH1F("fRecoLSeEmbTrkPt","Reco LS electrons from embedded #pi^{0} and #eta + No mom;p_{T} (GeV/c);counts",250,0,50);
        
        fRecoLSeEmbWeightTrkPt = new TH1F("fRecoLSeEmbWeightTrkPt","Reco LS electrons from embedded #pi^{0} and #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
        
        fRecoPi0LSeEmbWeightTrkPt = new TH1F("fRecoPi0LSeEmbWeightTrkPt","Reco LS electrons from embedded #pi^{0}  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
        
        fRecoEtaLSeEmbWeightTrkPt = new TH1F("fRecoEtaLSeEmbWeightTrkPt","Reco LS electrons from embedded #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
        
        fRecoULSeTrkPt = new TH1F("fRecoULSeTrkPt"," Reco ULS electrons from all generators;p_{T} (GeV/c);counts",250,0,50);
        
        fRecoULSeEmbTrkPt = new TH1F("fRecoULSeEmbTrkPt","Reco ULS electrons from embedded #pi^{0} and #eta + No mom;p_{T} (GeV/c);counts",250,0,50);
        
        fRecoULSeEmbWeightTrkPt = new TH1F("fRecoULSeEmbWeightTrkPt","Reco ULS electrons from embedded #pi^{0} and #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
        
        fRecoPi0ULSeEmbWeightTrkPt = new TH1F("fRecoPi0ULSeEmbWeightTrkPt","Reco ULS electrons from embedded #pi^{0}  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
        
        fRecoEtaULSeEmbWeightTrkPt = new TH1F("fRecoEtaULSeEmbWeightTrkPt","Reco ULS electrons from embedded #eta  + No mom with weight;p_{T} (GeV/c);counts",250,0,50);
        
        
        
        
        
        
        Int_t binw[4] = {500,3,2,7}; //pT, PDG, HijingOrNot, pi0etaType
        Double_t xminw[4] = {0,0,0,-1};
        Double_t xmaxw[4] = {50,3,2,6};
        
        fSprsPi0EtaWeightCal = new THnSparseD("fSprsPi0EtaWeightCal","Sparse to calculate #pi^{0} and #eta weight;p_{T};PDGID;HijingOrNot;pi0etaType;",4,binw,xminw,xmaxw);
        
        fPi0Weight = new TF1("fPi0Weight","[0] / TMath::Power(TMath::Exp(-[1]*x - [2]*x*x) + x/[3], [4])");
        
        fEtaWeight = new TF1("fEtaWeight","[0] / TMath::Power(TMath::Exp(-[1]*x - [2]*x*x) + x/[3], [4])");
        
        fPi0Weight->SetParameters( 1.46837e+03,-1.02586e-01, 1.18596e-03,1.73410,5.06623);
        fEtaWeight->SetParameters(   3.67328e+02,-5.34727e-02,2.39322e-05,1.96562,5.23828);
    }
    fTrkMatchTrkPt->Sumw2();
    fSparseElectron->Sumw2();
    fSparseLSElectron->Sumw2();
    fSparseULSElectron->Sumw2();
    fSparseMulti->Sumw2();
    fClusE->Sumw2();
    fSparseClusE->Sumw2();
    
    if(fReadMC){
        fHFElecPtAll->Sumw2();
        fHFElecPtReco_wtrkcuts->Sumw2();
        fHFElecPtReco_wTPCPID->Sumw2();
        fHFElecPtReco_wtrkCalocuts->Sumw2();
        fHFElecPtReco_wTPCCaloPID->Sumw2();
        fSprsPi0EtaWeightCal->Sumw2();
        fMissingEmbEtaEleTrkPt->Sumw2();
        fNonHFeTrkPt->Sumw2();
        fNonHFeEmbAllTypeTrkPt->Sumw2();
        fNonHFeEmbTrkPt->Sumw2();
        fNonHFeEmbWeightTrkPt->Sumw2();
        fPi0eEmbWeightTrkPt->Sumw2();
        fEtaeEmbWeightTrkPt->Sumw2();
        fNonHFewPIDTrkPt->Sumw2();
        fNonHFeEmbwPIDTrkPt->Sumw2();
        fNonHFeEmbWeightwPIDTrkPt->Sumw2();
        fPi0eEmbWeightwPIDTrkPt->Sumw2();
        fEtaeEmbWeightwPIDTrkPt->Sumw2();
        fRecoNonHFeTrkPt->Sumw2();
        fRecoNonHFeEmbTrkPt->Sumw2();
        fRecoNonHFeEmbWeightTrkPt->Sumw2();
        fRecoPi0eEmbWeightTrkPt->Sumw2();
        fRecoEtaeEmbWeightTrkPt->Sumw2();
        fRecoLSeTrkPt->Sumw2();
        fRecoLSeEmbTrkPt->Sumw2();
        fRecoLSeEmbWeightTrkPt->Sumw2();
        fRecoPi0LSeEmbWeightTrkPt->Sumw2();
        fRecoEtaLSeEmbWeightTrkPt->Sumw2();
        fRecoULSeTrkPt->Sumw2();
        fRecoULSeEmbTrkPt->Sumw2();
        fRecoULSeEmbWeightTrkPt->Sumw2();
        fRecoPi0ULSeEmbWeightTrkPt->Sumw2();
        fRecoEtaULSeEmbWeightTrkPt->Sumw2();
        
        
    }
    
    // Output list
    fOutputList->Add(fNevents);
    fOutputList->Add(fClusPhi);
    fOutputList->Add(fClusEta);
    fOutputList->Add(fClusEtaPhi);
    fOutputList->Add(fClusT);
    fOutputList->Add(fClusE);
    fOutputList->Add(fVtxZ);
    fOutputList->Add(fVtxY);
    fOutputList->Add(fVtxX);
    fOutputList->Add(fTPCdEdx);
    fOutputList->Add(fTPCdEdxvsEp);
    fOutputList->Add(fTPCnsigma);
    fOutputList->Add(fTrkPt);
    fOutputList->Add(fTrketa);
    fOutputList->Add(fTrkphi);
    fOutputList->Add(fTrkMatchTrkPt);
    fOutputList->Add(fTrkMatchTrketa);
    fOutputList->Add(fTrkMatchTrkphi);
    fOutputList->Add(fTrkMatchClusetaphi);
    fOutputList->Add(fEMCTrkMatchcluster);
    
    fOutputList->Add(fInvmassLS);
    fOutputList->Add(fInvmassULS);
    fOutputList->Add(fInvmassLSPt);
    fOutputList->Add(fInvmassULSPt);
    fOutputList->Add(fULSElecPt);
    fOutputList->Add(fLSElecPt);
    
    fOutputList->Add(fSparseClusE);
    fOutputList->Add(fSparseElectron);
    fOutputList->Add(fSparseLSElectron);
    fOutputList->Add(fSparseULSElectron);
    fOutputList->Add(fSparseMulti);
    
    if(fReadMC){
        fOutputList->Add(fSprsPi0EtaWeightCal);
        fOutputList->Add(fInclsElecPt);
        fOutputList->Add(fMissingEmbEtaEleTrkPt);
        fOutputList->Add(fInclsElecPtAll);
        fOutputList->Add(fInclsElecPtReco);
        fOutputList->Add(fInclseEMCALElecPtReco);
        fOutputList->Add(fInclseDCALElecPtReco);
        fOutputList->Add(fHFElecPtAll);
        fOutputList->Add(fHFElecPtReco_wtrkcuts);
        fOutputList->Add(fHFElecPtReco_wTPCPID);
        fOutputList->Add(fHFElecPtReco_wtrkCalocuts);
        fOutputList->Add(fHFElecPtReco_wTPCCaloPID);
        fOutputList->Add(fNonHFeTrkPt);
        fOutputList->Add(fNonHFeEmbAllTypeTrkPt);
        fOutputList->Add(fNonHFeEmbTrkPt);
        fOutputList->Add(fNonHFeEmbWeightTrkPt);
        fOutputList->Add(fPi0eEmbWeightTrkPt);
        fOutputList->Add(fEtaeEmbWeightTrkPt);
        fOutputList->Add(fNonHFewPIDTrkPt);
        fOutputList->Add(fNonHFeEmbwPIDTrkPt);
        fOutputList->Add(fNonHFeEmbWeightwPIDTrkPt);
        fOutputList->Add(fPi0eEmbWeightwPIDTrkPt);
        fOutputList->Add(fEtaeEmbWeightwPIDTrkPt);
        fOutputList->Add(fRecoNonHFeTrkPt);
        fOutputList->Add(fRecoNonHFeEmbTrkPt);
        fOutputList->Add(fRecoNonHFeEmbWeightTrkPt);
        fOutputList->Add(fRecoPi0eEmbWeightTrkPt);
        fOutputList->Add(fRecoEtaeEmbWeightTrkPt);
        fOutputList->Add(fNonHFePairInvmassLS);
        fOutputList->Add(fNonHFePairInvmassULS);
        fOutputList->Add(fNonHFeEmbInvmassLS);
        fOutputList->Add(fNonHFeEmbInvmassULS);
        fOutputList->Add(fNonHFeEmbWeightInvmassLS);
        fOutputList->Add(fNonHFeEmbWeightInvmassULS);
        fOutputList->Add(fPi0EmbInvmassULS);
        fOutputList->Add(fPi0EmbInvmassLS);
        fOutputList->Add(fPi0EmbWeightInvmassLS);
        fOutputList->Add(fPi0EmbWeightInvmassULS);
        fOutputList->Add(fEtaEmbInvmassLS);
        fOutputList->Add(fEtaEmbInvmassULS);
        fOutputList->Add(fEtaEmbWeightInvmassLS);
        fOutputList->Add(fEtaEmbWeightInvmassULS);
        fOutputList->Add(fRecoLSeTrkPt);
        fOutputList->Add(fRecoLSeEmbTrkPt);
        fOutputList->Add(fRecoLSeEmbWeightTrkPt);
        fOutputList->Add(fRecoPi0LSeEmbWeightTrkPt);
        fOutputList->Add(fRecoEtaLSeEmbWeightTrkPt);
        fOutputList->Add(fRecoULSeTrkPt);
        fOutputList->Add(fRecoULSeEmbTrkPt);
        fOutputList->Add(fRecoULSeEmbWeightTrkPt);
        fOutputList->Add(fRecoPi0ULSeEmbWeightTrkPt);
        fOutputList->Add(fRecoEtaULSeEmbWeightTrkPt);
        
        
    }
    
    
    PostData(1,fOutputList);
    PostData(2,fListProfiles);
    
}
//_____________________________________________________________________________
void AliAnalysisTaskHFEMultiplicity::UserExec(Option_t *)
{
    
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inputHandler=(AliInputEventHandler*)mgr->GetInputEventHandler();
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if(!fAOD) return;
    
    if(!PassEventSelect(fAOD)) return;
    
    
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
    //--------------------Combine EMCal and DCal----------------------------
    
    
    if(fFlagClsTypeEMC && fFlagClsTypeDCAL)
    {
        if(fEMCEG2 && fDCalDG2) if(!firedTrigger.Contains(TriggerEG2) && !firedTrigger.Contains(TriggerDG2)) return;
        if(fEMCEG1 && fDCalDG1) if(!firedTrigger.Contains(TriggerEG1) && !firedTrigger.Contains(TriggerDG1)) return;
    }
    // --- separate EMCAL and DCAL --- //
    else
    {
        if(fEMCEG1) if(!firedTrigger.Contains(TriggerEG1)) return;
        if(fEMCEG2) if(!firedTrigger.Contains(TriggerEG2)) return;
        if(fDCalDG1) if(!firedTrigger.Contains(TriggerDG1)) return;
        if(fDCalDG2) if(!firedTrigger.Contains(TriggerDG2)) return;
    }
    
    
    fNevents->Fill(1);
    
    
    //////////////////
    // Multiplicity //
    /////////////////
    Double_t Zvertex1 = -100;
    const AliAODVertex *pVtx = fAOD->GetPrimaryVertex();
    Zvertex1 =pVtx->GetZ();
    if(fzvtxcut){
        if(TMath::Abs(Zvertex1)>10.0) return;}
    fNevents->Fill(2);
    
    
    //--------------------vertex selection cuts-----------------------
    if(fzvtxQA){
        
        
        if (fRejectPUFromSPD && fAOD->IsPileupFromSPDInMultBins()) return; // pile-up cut
        fNevents->Fill(3);
        
        Double_t NContV = pVtx->GetNContributors();
        
        if(NContV<fCutNcontV) return; // n contributor
        fNevents->Fill(4);
        
        AliAODVertex* vtxSPD = fAOD->GetPrimaryVertexSPD();
        
        if (!vtxSPD || vtxSPD->GetNContributors() < 1) return;
        Double_t cov[6]={0};
        vtxSPD->GetCovarianceMatrix(cov);
        if (TMath::Sqrt(cov[5]) > 0.25) return;
        if (TMath::Abs(vtxSPD->GetZ() - pVtx->GetZ())>0.5) return;
    }
    
    fNevents->Fill(5);
    
    //----------V0M Multiplicity------------------
    AliAODVZERO *vzeroAOD = dynamic_cast<AliAODVZERO *>( dynamic_cast<AliAODEvent *>(fAOD)->GetVZEROData());
    Int_t V0AMult = static_cast<Int_t>(vzeroAOD->GetMTotV0A());
    Int_t V0CMult = static_cast<Int_t>(vzeroAOD->GetMTotV0C());
    Int_t V0Mult=V0AMult+V0CMult;
    
    //------------SPDTracklets--------------------
    Int_t nTracklets = 0;
    Int_t nAcc = 0;
    Double_t etaRange = 1.0;
    
    AliAODTracklets *tracklets = static_cast<const AliAODEvent*>(fAOD)->GetTracklets();
    nTracklets = tracklets->GetNumberOfTracklets();
    for (Int_t nn = 0; nn < nTracklets; nn++) {
        Double_t theta = tracklets->GetTheta(nn);
        Double_t eta = -TMath::Log(TMath::Tan(theta/2.0));
        if (TMath::Abs(eta) < etaRange) nAcc++;
    }
    // Data driven multiplicity z-vertex correction
    
    
    Int_t countMult = nAcc;
    
    //SPDTracklets correction
    Double_t correctednAcc   = nAcc;
    Double_t correctednAcc1   = nAcc;
    Double_t countCorr       = countMult;
    TProfile* estimatorAvg = GetEstimatorHistogram(fAOD);
    if(estimatorAvg){
        correctednAcc=static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,nAcc,Zvertex1,fRefMult));
        correctednAcc1 = AliAnalysisTaskHFEMultiplicity::GetTrackletsMeanCorrection(estimatorAvg,nAcc,Zvertex1,fRefMult);
        
        countCorr=static_cast<Int_t>(AliVertexingHFUtils::GetCorrectedNtracklets(estimatorAvg,countMult,Zvertex1,fRefMult));
    }
    
    
    //V0M Correction
    Int_t vzeroMultACorr=V0AMult, vzeroMultCCorr=V0CMult, vzeroMultCorr=V0Mult;
    vzeroMultACorr = static_cast<Int_t>(AliESDUtils::GetCorrV0A(V0AMult,Zvertex1));
    vzeroMultCCorr = static_cast<Int_t>(AliESDUtils::GetCorrV0C(V0CMult,Zvertex1));
    vzeroMultCorr = vzeroMultACorr + vzeroMultCCorr; // corrected V0M
    
    ////////////////////////////////
    //Monte Carlo information    //
    ////////////////////////////////
    
    if(fReadMC){
        
        fMCArray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
        if(!fMCArray){
            AliError("Array of MC particles not found");
            return;
        }
        fMCHeader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
        if (!fMCHeader) {
            AliError("Could not find MC Header in AOD");
            return;
        }
        
        ////////////////////////////////
        //Get number of Gen particles //
        ////////////////////////////////
        GetNMCPartProduced();
        
        
        /////////////////////////
        //Electrons in MC stack//
        /////////////////////////
        
        
        //electrons
        for(Int_t imcArrayL=0; imcArrayL< fMCArray->GetEntries(); imcArrayL++){
            AliAODMCParticle *AODMCtrack = (AliAODMCParticle*)fMCArray->At(imcArrayL);
            Int_t PDGcode = TMath::Abs(AODMCtrack->GetPdgCode());
            
            if(TMath::Abs(AODMCtrack->Eta()) > fCutTrackEta) continue;
            if(!AODMCtrack->IsPhysicalPrimary()) continue;
            
            
            if(PDGcode == 11) { // selecting only electrons
                fInclsElecPtAll->Fill(AODMCtrack->Pt());
                
                Int_t IsHFe=GetHFE(AODMCtrack,fMCArray);
                if((IsHFe==kBeauty) || (IsHFe==kCharm)){
                    fHFElecPtAll->Fill(AODMCtrack->Pt(),correctednAcc1);
                }
            }
        }
        
        /////////////////////////////////
        //Calculate Pi0 and Eta weight //
        /////////////////////////////////
        
        if(fCalculateWeight) GetPi0EtaWeight(fSprsPi0EtaWeightCal);
        
    }
    
    fvalueMulti[0] = Zvertex1;
    fvalueMulti[1] = V0Mult;
    fvalueMulti[2] = nAcc;
    fvalueMulti[3] = vzeroMultCorr;
    fvalueMulti[4] = correctednAcc1;
    fvalueMulti[5] = GetNcharged();
    
    
    fSparseMulti->Fill(fvalueMulti);    // multiplicity from tracklets
    
    
    fpidResponse = fInputHandler->GetPIDResponse();
    if(!fpidResponse) return;
    
    
    /////////////////////////
    // cluster information //
    ////////////////////////
    
    Double_t cluphi = -999.0;
    Double_t clueta =-999.0 ;
    Int_t ncells = -999.0;
    Float_t energy = -999.0;
    Float_t clut = -999.0;
    Double_t  energycell = -999.0;
    Double_t CellId =0;
    Int_t Nclust = -999;
    
    
    if(!fUseTender) Nclust = fAOD->GetNumberOfCaloClusters();
    if(fUseTender) Nclust = fCaloClusters_tender->GetEntries();
    
    Bool_t fClsTypeEMC = kFALSE, fClsTypeDCAL = kFALSE;
    for ( Int_t index = 0; index < Nclust ; index++ ) {
        AliAODCaloCluster * clu =0x0;
        if(!fUseTender) clu  = (AliAODCaloCluster*)fAOD->GetCaloCluster(index) ;
        if(fUseTender) clu = dynamic_cast<AliAODCaloCluster*>(fCaloClusters_tender->At(index));
        if(!clu) continue;
        
        fClsTypeEMC = kFALSE; fClsTypeDCAL = kFALSE;
        if (clu->IsEMCAL()){
            
            AliAODCaloCells &cells = *(fAOD->GetEMCALCells());
            
            Float_t  x[3]; // cluster pos
            clu->GetPosition(x);
            TVector3 clustposi(x[0],x[1],x[2]);
            
            cluphi = clustposi.Phi();
            clueta = clustposi.Eta();
            if(cluphi < 0) cluphi = cluphi+(2*TMath::Pi());
            if(cluphi > 1.39 && cluphi < 3.265) fClsTypeEMC = kTRUE; //EMCAL : 80 < phi < 187
            if(cluphi > 4.53 && cluphi < 5.708) fClsTypeDCAL = kTRUE;//DCAL  : 260 < phi < 327
            
            if(fFlagClsTypeEMC && !fFlagClsTypeDCAL)
                if(!fClsTypeEMC) continue; //selecting only EMCAL clusters
            
            if(fFlagClsTypeDCAL && !fFlagClsTypeEMC)
                if(!fClsTypeDCAL) continue; //selecting only DCAL clusters
            
            clut = clu->GetTOF()*1e9 ;
            energy = clu->E();
            ncells= clu->GetNCells();
            fClusPhi->Fill(cluphi);
            fClusEtaPhi->Fill(clueta,cluphi);
            fClusEta->Fill(clueta);
            fClusE->Fill(energy);
            fClusT->Fill(clut);
            
            
            fvalueCluE[0] = energy;
            fvalueCluE[1] = vzeroMultCorr; //V0M, Multiplicity information
            fvalueCluE[2] = correctednAcc1; //SPD Tracklets
            
            fSparseClusE->Fill(fvalueCluE); //For Rejection Factor
        }
        
        
        
    }
    ///////////////////////
    // Track information //
    //////////////////////
    
    Int_t ntracks = -999;
    if(!fUseTender)ntracks = fAOD->GetNumberOfTracks();
    if(fUseTender) ntracks = fTracks_tender->GetEntries();
    
    
    for (Int_t iTracks=0; iTracks< ntracks; iTracks++) {
        
        AliAODTrack* track = 0x0;
        if(!fUseTender) track = (AliAODTrack*)fAOD->GetTrack(iTracks);
        if(fUseTender) track =  dynamic_cast<AliAODTrack*>(fTracks_tender->At(iTracks));
        if(!track) continue;
        if(!Passtrackcuts(track)) continue;
        
        
        Double_t TrkPhi=-999, TrkPt=-999, TrkEta=-999, TrkP = -999,dEdx = -999, nsigma = -999.0;
        
        nsigma=fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
        
        dEdx = track->GetTPCsignal();
        
        TrkPhi = track->Phi();
        fTrkphi->Fill(TrkPhi);
        TrkPt = track->Pt();
        fTrkPt->Fill(TrkPt);
        
        TrkEta = track->Eta();
        fTrketa->Fill(TrkEta);
        TrkP = track->P();
        
        fTPCdEdx->Fill(TrkP,dEdx);
        fTPCnsigma->Fill(TrkP,nsigma);
        
        
        if(fReadMC){
            Int_t iTrklabel = TMath::Abs(track->GetLabel());
            if(iTrklabel == 0) continue;
            AliAODMCParticle *MCPart = (AliAODMCParticle*)fMCArray->At(iTrklabel);
            if(TMath::Abs(MCPart->Eta()) > fCutTrackEta) continue;
            if(!MCPart->IsPhysicalPrimary()) continue;
            if(TMath::Abs(MCPart->GetPdgCode())==11) { //only electrons
                Int_t IsElecHf=GetHFE(MCPart,fMCArray);
                if((IsElecHf==kBeauty) || (IsElecHf==kCharm)){ //HF electrons
                    fHFElecPtReco_wtrkcuts -> Fill(TrkPt,correctednAcc1);
                    if ((nsigma > fCutNsigmaEMin)  && (nsigma < fCutNsigmaEMin )){
                        fHFElecPtReco_wTPCPID -> Fill(TrkPt,correctednAcc1); }
                }
            }
        }
        
        if(TrkPt < fCutpTMin) continue ;
        Int_t EMCalIndex = -1;
        Double_t emcphi = -999, emceta = -999;
        Bool_t fClsTypeEMC = kFALSE, fClsTypeDCAL = kFALSE;
        EMCalIndex = track -> GetEMCALcluster();
        AliAODCaloCluster *clustMatch=0x0;
        
        if(!fUseTender) if(EMCalIndex >= 0) clustMatch = (AliAODCaloCluster*)fAOD->GetCaloCluster(EMCalIndex) ;
        if(fUseTender) if(EMCalIndex >= 0)clustMatch = dynamic_cast<AliAODCaloCluster*>(fCaloClusters_tender->At(EMCalIndex));
        if(clustMatch && clustMatch->IsEMCAL())
        {
            Double_t fPhiDiff = -999, fEtaDiff = -999;
            GetTrkClsEtaPhiDiff(track, clustMatch, fPhiDiff, fEtaDiff);
            fEMCTrkMatchcluster->Fill(fPhiDiff,fEtaDiff);
            if(TMath::Abs(fPhiDiff) > 0.06 || TMath::Abs(fEtaDiff)> 0.06) continue;
            
            Float_t EMCalpos[3];
            clustMatch -> GetPosition(EMCalpos);
            TVector3 clustpos(EMCalpos[0],EMCalpos[1],EMCalpos[2]);
            emcphi = clustpos.Phi();
            if(emcphi < 0) emcphi = emcphi + (2*TMath::Pi());
            emceta = clustpos.Eta();
            
            if(emcphi > 1.39 && emcphi < 3.265) fClsTypeEMC = kTRUE;
            if(emcphi > 4.53 && emcphi < 5.708) fClsTypeDCAL = kTRUE;
            
            if(fFlagClsTypeEMC && !fFlagClsTypeDCAL)
                if(!fClsTypeEMC) continue; //selecting only EMCAL clusters
            
            if(fFlagClsTypeDCAL && !fFlagClsTypeEMC)
                if(!fClsTypeDCAL) continue; //selecting only DCAL clusters
            
            fTrkMatchTrkPt->Fill(TrkPt);
            fTrkMatchTrketa->Fill(TrkEta);
            fTrkMatchTrkphi->Fill(TrkPhi);
            fTrkMatchClusetaphi->Fill(emceta,emcphi);
            
            Double_t Etrkmatch = -999.0, Eoptrk = -999.0 , M02trkmatch = -999.0, M20trkmatch = -999.0;
            
            Etrkmatch = clustMatch->E();
            Eoptrk = Etrkmatch/TrkP ;
            M02trkmatch = clustMatch->GetM02();
            M20trkmatch = clustMatch->GetM20();
            
            Bool_t fElectTrack = kFALSE;
            fElectTrack = PassEIDCuts(track, clustMatch);
            fTPCdEdxvsEp->Fill(Eoptrk,dEdx);
        
                fvalueElectron[0] = TrkPt; //matched tracks pt
                fvalueElectron[1] = nsigma; // tpc n sigma
                fvalueElectron[2] = Eoptrk; //E/P
                fvalueElectron[3] = M02trkmatch; // shower shape cut
                fvalueElectron[4] = M20trkmatch;
                fvalueElectron[5] = vzeroMultCorr; //V0M, Multiplicity information
                fvalueElectron[6] = correctednAcc1; //SPD Tracklets
                fvalueElectron[7] = Etrkmatch;  //cluster energy after matching
                fvalueElectron[8] = Zvertex1;
                
                fSparseElectron->Fill(fvalueElectron);   //Electron information sparse
                
            
            
            
            
            
            Int_t iMCmom=-999, MomPDG = -999;
            Double_t MomPt =-999;
            
            if(fReadMC){Int_t iTrklabel = TMath::Abs(track->GetLabel());
                if(iTrklabel == 0) continue;
                AliAODMCParticle *MCPart = (AliAODMCParticle*)fMCArray->At(iTrklabel);
                if(TMath::Abs(MCPart->GetPdgCode())!=11) continue;
                fInclsElecPt->Fill(TrkPt);
                if(fElectTrack){
                    fInclsElecPtReco->Fill(TrkPt);
                    if(fClsTypeEMC) fInclseEMCALElecPtReco->Fill(TrkPt);
                    if(fClsTypeDCAL) fInclseDCALElecPtReco->Fill(TrkPt);
                }
                
                
                
                Int_t IsElecHf=GetHFE(MCPart,fMCArray);
                if((IsElecHf==kBeauty) || (IsElecHf==kCharm)){
                    fHFElecPtReco_wtrkCalocuts -> Fill(TrkPt,correctednAcc1);
                    if (fElectTrack){
                        fHFElecPtReco_wTPCCaloPID -> Fill(TrkPt,correctednAcc1); }
                }
                
                
                //----------photon tagging efficiency by invariant mass method-------------------------------
                //   if(!fCalculateNonHFEEffi) continue;
                
                ////////////////////
                //NonHFE selection//
                ////////////////////
                
                fIsFrmEmbPi0 = kFALSE;
                fIsFrmEmbEta= kFALSE;
                ftype = -1;
                fWeight = 1.0;
                Bool_t fFromHijing = kTRUE;
                
                Bool_t fNonHFE = IsNonHFE(MCPart, fFromHijing, ftype, iMCmom, MomPDG, MomPt);
                
                if(!fNonHFE) continue;
                fNonHFeTrkPt->Fill(TrkPt);
                if(fElectTrack){ fNonHFewPIDTrkPt->Fill(TrkPt);    }
                ///////////////////////////////////////
                // Check for pi0/eta from embbedding //
                ///////////////////////////////////////
                //not considering the cases : eta->pi0->gamma->ele, eta->pi0->elec, check in next iteration
                //Cases considered : pi0->e, pi0->gamma->e, eta->e, eta->gamma->e
                Bool_t IsMissingEta;
                Int_t missingtype=-1;
                
                if(MomPDG == 111) {
                    if(iMCmom >= fNembMCpi0 && iMCmom < fNembMCeta) fIsFrmEmbPi0 = kTRUE;
                    
                    //missing eta
                    AliAODMCParticle *MCPartMom = (AliAODMCParticle*)fMCArray->At(iMCmom);
                    Int_t iMCgmom = -999;
                    iMCgmom = MCPartMom->GetMother();
                    if(iMCgmom > 0) {
                        AliAODMCParticle *MCPartGMom = (AliAODMCParticle*)fMCArray->At(iMCgmom);
                        Int_t GMomPDG = TMath::Abs(MCPartGMom->GetPdgCode());
                        
                        if(GMomPDG == 221){
                            if(iMCgmom >= fNembMCeta && iMCgmom < fNTotMCpart){
                                IsMissingEta = kTRUE;
                                missingtype = GetPi0EtaType(MCPartGMom);
                            }
                        }
                    }
                }
                
                if(MomPDG == 221){
                    if(iMCmom >= fNembMCeta && iMCmom < fNTotMCpart) fIsFrmEmbEta = kTRUE;
                }
                
                if(MomPDG == 22) //if pi0/eta->gamma->e : Rewrite "ftype" and "MomPt" with Gmother info
                {
                    AliAODMCParticle *MCPartMom = (AliAODMCParticle*)fMCArray->At(iMCmom);
                    Int_t iMCgmom = -999;
                    iMCgmom = MCPartMom->GetMother();
                    AliAODMCParticle *MCPartGMom = (AliAODMCParticle*)fMCArray->At(iMCgmom);
                    Int_t GMomPDG = TMath::Abs(MCPartGMom->GetPdgCode());
                    
                    
                    if(GMomPDG == 111){
                        if(iMCgmom >= fNembMCpi0 && iMCgmom < fNembMCeta){
                            fIsFrmEmbPi0 = kTRUE;
                            ftype = GetPi0EtaType(MCPartGMom);
                            MomPt = MCPartGMom->Pt();
                        }
                        //missing eta
                        Int_t iMCggmom = -999;
                        if(iMCggmom > 0) {
                            iMCggmom = MCPartGMom->GetMother();
                            AliAODMCParticle *MCPartGGMom = (AliAODMCParticle*)fMCArray->At(iMCggmom);
                            Int_t GGMomPDG = TMath::Abs(MCPartGGMom->GetPdgCode());
                            
                            if(GGMomPDG == 221){
                                if(iMCggmom >= fNembMCeta && iMCggmom < fNTotMCpart){
                                    IsMissingEta = kTRUE;
                                    missingtype = GetPi0EtaType(MCPartGGMom);
                                }
                            }
                        }
                    }
                    if(GMomPDG == 221){
                        if(iMCgmom >= fNembMCeta && iMCgmom < fNTotMCpart){
                            fIsFrmEmbEta = kTRUE;
                            ftype = GetPi0EtaType(MCPartGMom);
                            MomPt = MCPartGMom->Pt();
                        }
                    }
                }
                
                if(IsMissingEta && (missingtype == kNoMother)) fMissingEmbEtaEleTrkPt->Fill(TrkPt);
                
                //////////////////////////////////////////////////
                ///Get weight for Embedded pi0/eta with NoMother//
                //////////////////////////////////////////////////
                
                if(fIsFrmEmbPi0 && ftype==kNoMother) {
                    fWeight = fPi0Weight->Eval(MomPt);
                }
                if(fIsFrmEmbEta && ftype==kNoMother) {
                    fWeight = fEtaWeight->Eval(MomPt);
                }
                
                
                //////////////////////////////////////////
                //Select electrons from embedded pi0/eta//
                //////////////////////////////////////////
                
                if(fIsFrmEmbPi0 || fIsFrmEmbEta){
                    fNonHFeEmbAllTypeTrkPt->Fill(TrkPt);
                    
                    if(ftype == kNoMother){ //embedded pi0/eta with no Mom
                        fNonHFeEmbTrkPt->Fill(TrkPt);
                        fNonHFeEmbWeightTrkPt->Fill(TrkPt,fWeight);
                        if(fIsFrmEmbPi0) fPi0eEmbWeightTrkPt->Fill(TrkPt,fWeight);
                        if(fIsFrmEmbEta) fEtaeEmbWeightTrkPt->Fill(TrkPt,fWeight);
                        
                        
                        if(fElectTrack){
                            fNonHFeEmbwPIDTrkPt->Fill(TrkPt);
                            fNonHFeEmbWeightwPIDTrkPt->Fill(TrkPt,fWeight);
                            if(fIsFrmEmbPi0) fPi0eEmbWeightwPIDTrkPt->Fill(TrkPt,fWeight);
                            if(fIsFrmEmbEta) fEtaeEmbWeightwPIDTrkPt->Fill(TrkPt,fWeight);
                        }
                    }
                }
            }
            //////////////////////////////////////
            //Reconst NonHFE with invmass method//
            //////////////////////////////////////
            
            if(fElectTrack){
                fvaluePHElectron[0] = TrkPt;
                fvaluePHElectron[1] = vzeroMultCorr; //V0M, Multiplicity information
                fvaluePHElectron[2] = correctednAcc1; //SPD Tracklets
                
                Bool_t fFlagPhotonicElec = kFALSE, fFlagElecLS=kFALSE;
                SelectNonHFElectron(iTracks, track, iMCmom, MomPDG, fFlagPhotonicElec, fFlagElecLS, vzeroMultCorr, correctednAcc1);
                
                if(fReadMC){
                    if(fFlagPhotonicElec){
                        fRecoNonHFeTrkPt->Fill(TrkPt);
                        if((fIsFrmEmbPi0 || fIsFrmEmbEta) && (ftype == kNoMother)){
                            fRecoNonHFeEmbTrkPt->Fill(TrkPt);// reconstructed embedded pi0/eta with no Mom using invariant mass
                            fRecoNonHFeEmbWeightTrkPt->Fill(TrkPt,fWeight);
                            if(fIsFrmEmbPi0) fRecoPi0eEmbWeightTrkPt->Fill(TrkPt,fWeight);
                            if(fIsFrmEmbEta) fRecoEtaeEmbWeightTrkPt->Fill(TrkPt,fWeight);
                        }
                    }
                }
                
            }
            
        }//track match
        
    } //track loop
    
    
    PostData(1,fOutputList);
    PostData(2,fListProfiles);
}

//______________________________________________________________________

Bool_t  AliAnalysisTaskHFEMultiplicity::IsNonHFE(AliAODMCParticle *MCPart, Bool_t &fFromHijing, Int_t &type, Int_t &iMCmom, Int_t &MomPDG, Double_t &MomPt)
{
    //Is electron from pi0, eta and gamma
    
    iMCmom = MCPart->GetMother();
    AliAODMCParticle *MCPartMom = (AliAODMCParticle*)fMCArray->At(iMCmom);
    MomPDG = TMath::Abs(MCPartMom->GetPdgCode());
    MomPt = MCPartMom->Pt();
    
    if((MomPDG == 111) || (MomPDG == 221) || (MomPDG == 22)){
        if(iMCmom >= fNpureMC)fFromHijing = kFALSE;
        type = GetPi0EtaType(MCPartMom);
        return kTRUE;
    }
    else return kFALSE;
}


//______________________________________________________________________
Bool_t AliAnalysisTaskHFEMultiplicity::PassEIDCuts(AliAODTrack *track, AliAODCaloCluster *clust)
{
    //apply electron identification cuts
    
    Double_t eop = -1.0;
    Double_t m02 = -999,m20 = -999;
    Double_t clustE = clust->E();
    Double_t TrkPt = track->Pt();
    Double_t nsigma_ele=-999;
    nsigma_ele = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
    if(track->P()>0)eop = clustE/track->P();
    
    m20 =clust->GetM20();
    
    if(nsigma_ele < fCutNsigmaEMin || nsigma_ele > fCutNsigmaEMax) return kFALSE;
    
    if(m20 < fCutM20Min || m20 > fCutM20Max) return kFALSE;
    if(eop < fCutEopEMin || eop > fCutEopEMax) return kFALSE;
    
    return kTRUE;
}
//______________________________________________________________________
Bool_t AliAnalysisTaskHFEMultiplicity::Passtrackcuts(AliAODTrack *atrack)
{
    Double_t d0z0[2]={-999,-999}, cov[3];
    Double_t DCAxyCut = fCutDCAxy, DCAzCut = fCutDCAz;
    Double_t dEdx =-999;
    Double_t TrkPhi=-999, TrkPt=-999, TrkEta=-999, TrkP = -999;
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    const AliAODVertex *pVtx = fAOD->GetPrimaryVertex();
    
    //kink daughters
    Int_t numberofvertices = 100;
    numberofvertices = fAOD->GetNumberOfVertices();
    Double_t listofmotherkink[numberofvertices];
    Int_t numberofmotherkink = 0;
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
    
    //reject kink
    Bool_t kinkmotherpass = kTRUE;
    for(Int_t kinkmother = 0; kinkmother < numberofmotherkink; kinkmother++) {
        if(atrack->GetID() == listofmotherkink[kinkmother]) {
            kinkmotherpass = kFALSE;
            continue;
        }
    }
    if(!kinkmotherpass) return kFALSE;
    
    //other cuts
    if(!atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return kFALSE; //minimum cuts- filter bit 4
    if(atrack->GetTPCCrossedRows() < fCutTPCMaxCls) return kFALSE;
    if(atrack->Chi2perNDF() >= fCutTPCchi2perNDF) return kFALSE;
    if(atrack->GetTPCNcls() < fCutTPCNCls) return kFALSE;
    if(atrack->GetITSNcls() < fCutITSNCls) return kFALSE;
    if((!(atrack->GetStatus()&AliAODTrack::kITSrefit)|| (!(atrack->GetStatus()&AliAODTrack::kTPCrefit)))) return kFALSE;
    if(!(atrack->HasPointOnITSLayer(0) || atrack->HasPointOnITSLayer(1))) return kFALSE;
    
    if(atrack->PropagateToDCA(pVtx, fAOD->GetMagneticField(), 20., d0z0, cov))
        if(TMath::Abs(d0z0[0]) > DCAxyCut || TMath::Abs(d0z0[1]) > DCAzCut) return kFALSE;
    
    return kTRUE;
    
}

//______________________________________________________________________
Bool_t AliAnalysisTaskHFEMultiplicity::PassEventSelect(AliAODEvent *fAOD)
{
    Int_t ntracks = -999;
    fNevents->Fill(0);
    Double_t Zvertex=-100, Xvertex=-100, Yvertex=-100;
    
    const AliAODVertex *pVtx = fAOD->GetPrimaryVertex();
    
    // Double_t NContV = pVtx->GetNContributors();
    
    // if(NContV<fCutNcontV) return kFALSE;
    //fNevents->Fill(1);
    
    Zvertex =pVtx->GetZ();
    Yvertex =pVtx->GetY();
    Xvertex =pVtx->GetX();
    
    //if(TMath::Abs(Zvertex)>10.0) return kFALSE;
    //fNevents->Fill(2);
    
    
    fVtxZ->Fill(Zvertex);
    fVtxX->Fill(Xvertex);
    fVtxY->Fill(Yvertex);
    
    return kTRUE;
    
}
//______________________________________________________________________
Int_t AliAnalysisTaskHFEMultiplicity::GetHFE(AliAODMCParticle *electron, TClonesArray *mcArray)
{
    Int_t motherindex=electron->GetMother(); //Getting Electron Mother
    if(motherindex<0) kNoMother;
    AliAODMCParticle *mother = (AliAODMCParticle*)mcArray->At(motherindex);
    Int_t motherpdg = mother->GetPdgCode();
    if ( (motherindex >= fNpureMC) &&(int(TMath::Abs(motherpdg)/100.)%10) == 5 || (int(TMath::Abs(motherpdg)/1000.)%10) == 5 ) return kBeauty;
    if ( (motherindex >= fNpureMC) && (int(TMath::Abs(motherpdg)/100.)%10) == 4 || (int(TMath::Abs(motherpdg)/1000.)%10) == 4 ) return kCharm;
    
}

//______________________________________________________________________
void AliAnalysisTaskHFEMultiplicity::SelectNonHFElectron(Int_t itrack, AliAODTrack *track, Int_t iMCmom, Int_t MomPDG, Bool_t &fFlagPhotonicElec, Bool_t &fFlagElecLS, Int_t vzeroMultCorr,Int_t correctednAcc1)
{
    //Photonic electron selection
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    const AliAODVertex *pVtx = fAOD->GetPrimaryVertex();
    Double_t d0z0[2]={-999,-999}, cov[3];
    Double_t DCAxyCut = fCutDCAxy, DCAzCut = fCutDCAz;
    
    Bool_t flagPhotonicElec = kFALSE, flagLSElec = kFALSE;
    Double_t ptAsso=-999., nsigmaAsso=-999.0;
    Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
    
    
    Int_t ntracks = -999;
    if(!fUseTender)ntracks = fAOD->GetNumberOfTracks();
    if(fUseTender) ntracks = fTracks_tender->GetEntries();
    
    for(Int_t jTracks = 0; jTracks < ntracks; jTracks++){
        if(jTracks==itrack) continue;
        
        
        AliAODTrack* atrackAsso = 0x0;
        if(!fUseTender) atrackAsso = (AliAODTrack*)fAOD->GetTrack(jTracks);
        if(fUseTender) atrackAsso =  dynamic_cast<AliAODTrack*>(fTracks_tender->At(jTracks));
        if(!atrackAsso) continue;
        
        if(!atrackAsso->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
        if(atrackAsso->GetTPCCrossedRows() < fCutTPCMaxCls) continue;
        if(atrackAsso->Chi2perNDF() >= fCutTPCchi2perNDF) continue;
        if(atrackAsso->GetTPCNcls() < fCutTPCNCls) continue;
        if(atrackAsso->GetITSNcls() < fCutITSNCls) continue;
        if((!(atrackAsso->GetStatus()&AliAODTrack::kITSrefit)|| (!(atrackAsso->GetStatus()&AliAODTrack::kTPCrefit)))) continue; //refit required
        if(!(atrackAsso->HasPointOnITSLayer(0) || atrackAsso->HasPointOnITSLayer(1))) continue;
        
        
        nsigmaAsso = fpidResponse->NumberOfSigmasTPC(atrackAsso, AliPID::kElectron);
        ptAsso = atrackAsso->Pt();
        Int_t chargeAsso = atrackAsso->Charge();
        Int_t charge = track->Charge();
        
        if(ptAsso < fCutAssoEPt) continue; //partner electron pt 0.3
        if(atrackAsso->Eta()<-fCutAssoEEta || atrackAsso->Eta()>fCutAssoEEta) continue;
        if(nsigmaAsso < -fCutAssoENsigma || nsigmaAsso > fCutAssoENsigma) continue;
        
        
        if(atrackAsso->PropagateToDCA(pVtx, fAOD->GetMagneticField(), 20., d0z0, cov))
            if(TMath::Abs(d0z0[0]) > DCAxyCut || TMath::Abs(d0z0[1]) > DCAzCut) continue;
        
        Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;
        if(charge>0) fPDGe1 = -11;
        if(chargeAsso>0) fPDGe2 = -11;
        
        fFlagLS=kFALSE; fFlagULS=kFALSE;
        if(charge == chargeAsso) fFlagLS = kTRUE;
        if(charge != chargeAsso) fFlagULS = kTRUE;
        
        AliKFParticle::SetField(fAOD->GetMagneticField());
        
        AliKFParticle ge1 = AliKFParticle(*track, fPDGe1);
        AliKFParticle ge2 = AliKFParticle(*atrackAsso, fPDGe2);
        AliKFParticle recg(ge1,ge2);
        
        if(recg.GetNDF()<1) continue;
        Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
        if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;
        
        Double_t mass=-999., width = -999.;
        Int_t MassCorrect;
        MassCorrect = recg.GetMass(mass,width);
        
        if(fFlagLS && track->Pt()>1) fInvmassLS->Fill(mass);
        if(fFlagULS && track->Pt()>1) fInvmassULS->Fill(mass);
        if(fFlagLS) fInvmassLSPt->Fill(track->Pt(),mass);
        if(fFlagULS) fInvmassULSPt->Fill(track->Pt(),mass);
        
        if(fReadMC){
            Int_t iTrkAssolabel = TMath::Abs(atrackAsso->GetLabel());
            if(iTrkAssolabel == 0) continue;
            AliAODMCParticle *MCPartAsso = (AliAODMCParticle*)fMCArray->At(iTrkAssolabel);
            
            
            if(TMath::Abs(MCPartAsso->GetPdgCode())==11){ // check origin of asso elec
                Bool_t fAssoFromHijing = kTRUE;
                Int_t iMCAssomom=-999, AssoMomPDG = -999, fAssotype=-1;
                Double_t AssoMomPt =-999;
                Bool_t fAssoNonHFE = IsNonHFE(MCPartAsso, fAssoFromHijing, fAssotype, iMCAssomom, AssoMomPDG, AssoMomPt);
                
                if(fAssoNonHFE){
                    if(fFlagLS) fNonHFePairInvmassLS->Fill(mass);
                    if(fFlagULS) fNonHFePairInvmassULS->Fill(mass);
                }
            }
            
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
        }
        Double_t TrkPt = track->Pt();
        if(mass < fCutInvmass){
            if(fFlagLS){
                fLSElecPt->Fill(track->Pt()); //Reco LS e TrkPt
                fSparseLSElectron->Fill(fvaluePHElectron);
                
                if(fReadMC){
                    if((fIsFrmEmbPi0 || fIsFrmEmbEta) && (ftype == kNoMother)){
                        fRecoLSeEmbTrkPt->Fill(TrkPt);
                        fRecoLSeEmbWeightTrkPt->Fill(TrkPt,fWeight);
                        if(fIsFrmEmbPi0) fRecoPi0LSeEmbWeightTrkPt->Fill(TrkPt,fWeight);
                        if(fIsFrmEmbEta) fRecoEtaLSeEmbWeightTrkPt->Fill(TrkPt,fWeight);
                    }
                }
            }
            if(fFlagULS){
                fULSElecPt->Fill(track->Pt());  //Reco ULS e TrkPt
                fSparseULSElectron->Fill(fvaluePHElectron);
                if(fReadMC){
                    if((fIsFrmEmbPi0 || fIsFrmEmbEta) && (ftype == kNoMother)){
                        fRecoULSeEmbTrkPt->Fill(TrkPt);
                        fRecoULSeEmbWeightTrkPt->Fill(TrkPt,fWeight);
                        if(fIsFrmEmbPi0) fRecoPi0ULSeEmbWeightTrkPt->Fill(TrkPt,fWeight);
                        if(fIsFrmEmbEta) fRecoEtaULSeEmbWeightTrkPt->Fill(TrkPt,fWeight);
                    }
                }
            }
        }
        
        if(mass<fCutInvmass && fFlagULS && !flagPhotonicElec){
            flagPhotonicElec = kTRUE;
        }
        if(mass<fCutInvmass && fFlagULS && !flagPhotonicElec){
            flagLSElec = kTRUE;
        }
    }
    fFlagPhotonicElec = flagPhotonicElec;
    fFlagElecLS = flagLSElec;
    
}
//____________________________________________________________________________
TProfile* AliAnalysisTaskHFEMultiplicity::GetEstimatorHistogram(const AliAODEvent* fAOD)
{
    
    Int_t runNo  = fAOD->GetRunNumber();
    //cout<<"run number"<<runNo<<endl;
    Int_t period = -1;
    
    
    if (runNo>266436 && runNo<267111) period = 0;
    if (runNo>265743 && runNo<266319) period = 1;
    if (period < 0 || period > 1) return 0;
    
    //cout<<"period ="<<period<<endl;
    
    return fMultEstimatorAvg[period];
}
//____________________________________________________________________________
void AliAnalysisTaskHFEMultiplicity::GetTrkClsEtaPhiDiff(AliAODTrack *t, AliAODCaloCluster *v, Double_t &phidiff, Double_t &etadiff)
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

//____________________________________________________________________________


//______________________________________________________________________________
Int_t AliAnalysisTaskHFEMultiplicity::GetNcharged(){
    //counts all tracks in eta<1 with charge!=0
    
    Int_t Nch = 0;
    
    if(!fReadMC) return Nch; // if no MC info return 0
    
    // loop over all tracks
    for (Int_t igen = 0; igen < fMCArray->GetEntriesFast(); igen++){
        AliAODMCParticle *mctrack=(AliAODMCParticle*)fMCArray->UncheckedAt(igen);
        Int_t charge = mctrack->Charge();
        Double_t eta = mctrack->Eta();
        Bool_t isPhysPrim = mctrack->IsPhysicalPrimary();
        if(charge!=0){
            if(eta > -1.0 && eta < 1.0){
                if(isPhysPrim){
                    Nch++;
                }
            }
        }
    }
    return Nch;
}
//______________________________________________________________________________
void AliAnalysisTaskHFEMultiplicity::GetPi0EtaWeight(THnSparse *SparseWeight)
{
    //Get pi0 and eta information for weight calculation
    
    Double_t fvalue[4] = {-999,-999,-999,-999};
    
    for(int imc=0; imc< fNTotMCpart; imc++)
    {
        AliAODMCParticle *AODMCtrack = (AliAODMCParticle*)fMCArray->At(imc);
        if(TMath::Abs(AODMCtrack->Eta()) > 1.2) continue;
        
        //-------Get PDG
        Int_t TrackPDG = TMath::Abs(AODMCtrack->GetPdgCode());
        if((TrackPDG != 111) && (TrackPDG != 221) && (TrackPDG != 22)) continue;
        
        Double_t fPartPDGid = -999;
        if (TrackPDG == 111) fPartPDGid = 0.2;
        if (TrackPDG == 221) fPartPDGid = 1.2;
        if (TrackPDG == 22) fPartPDGid = 2.2;
        
        Double_t fTrkPt = AODMCtrack->Pt();
        
        //-------Check if the particle is from hijing or not
        Bool_t fFromHijing = kHijing;
        if(imc >= fNpureMC)fFromHijing = kElse;
        
        //------Get type of the particle
        Int_t fType = GetPi0EtaType(AODMCtrack);
        
        fvalue[0] = fTrkPt;
        fvalue[1] = fPartPDGid;
        fvalue[2] = fFromHijing;
        fvalue[3] = fType;
        
        SparseWeight->Fill(fvalue);
    }
}
//______________________________________________________________________________
Int_t AliAnalysisTaskHFEMultiplicity::GetPi0EtaType(AliAODMCParticle *part)
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
//______________________________________________________________________________
Bool_t AliAnalysisTaskHFEMultiplicity::GetNMCPartProduced()
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
        
        if(igene==0) fNpureMC = gh->NProduced();  // generated by HIJING
        
        
        if(MCgen.Contains(embpi0))fNembMCpi0 = fNTotMCpart;
        if(MCgen.Contains(embeta))fNembMCeta = fNTotMCpart;
        fNTotMCpart += gh->NProduced();
    }
    
    
    return kTRUE;
}
//______________________________________________________________________________
Double_t AliAnalysisTaskHFEMultiplicity::GetTrackletsMeanCorrection(TProfile* estimatorAvg, Double_t uncorrectedNacc, Double_t vtxZ, Double_t refMult)
{
    if(TMath::Abs(vtxZ)>10.0){
        //    printf("ERROR: Z vertex out of range for correction of multiplicity\n");
        return uncorrectedNacc;
    }
    
    if(!estimatorAvg){
        printf("ERROR: Missing TProfile for correction of multiplicity\n");
        return uncorrectedNacc;
    }
    
    Double_t localAvg = estimatorAvg->GetBinContent(estimatorAvg->FindBin(vtxZ));
    Double_t deltaM = 0;
    deltaM = uncorrectedNacc*(refMult/localAvg - 1);
    
    Double_t correctedNacc = uncorrectedNacc + (deltaM>0 ? 1 : -1) * gRandom->PoissonD(TMath::Abs(deltaM));
    
    if(correctedNacc<0) correctedNacc=0;
    
    return correctedNacc;
    
}


//______________________________________________________________________________

void AliAnalysisTaskHFEMultiplicity::Terminate(Option_t *)
{
    
}
//_________________________________________________________________________



