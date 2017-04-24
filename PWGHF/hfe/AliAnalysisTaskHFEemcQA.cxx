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
//    QA task for EMCAL electron analysis   //
//  Author: Deepa Thomas, Shingo Sakai      //
//////////////////////////////////////////////

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "THnSparse.h"

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

//#include "AliQnCorrectionsManager.h"

#include "AliAnalysisTaskHFEemcQA.h"

using std::cout;
using std::endl;

ClassImp(AliAnalysisTaskHFEemcQA)
//________________________________________________________________________
AliAnalysisTaskHFEemcQA::AliAnalysisTaskHFEemcQA(const char *name)
: AliAnalysisTaskSE(name),
fVevent(0),
fESD(0),
fAOD(0),
fMCheader(0),
fpidResponse(0),
fEMCALGeo(0),
fFlagSparse(kFALSE),
fUseTender(kTRUE),
fEMCEG1(kFALSE),
fEMCEG2(kFALSE),
fDCalDG1(kFALSE),
fDCalDG2(kFALSE),
fTracks_tender(0),
fCaloClusters_tender(0),
fMCparticle(0),
fMCarray(0),
fMultSelection(0),
fTriggersInfo(0),
fThresholdEG2(89),
fThresholdEG1(140),
fFlagClsTypeEMC(kTRUE),
fFlagClsTypeDCAL(kTRUE),
fcentMim(0),
fcentMax(0),
fCentralityEstimator("V0M"),
fOutputList(0),
fNevents(0),
fCent(0),
fMult(0),
fEvPlaneV0(0),
fEvPlaneV0A(0),
fEvPlaneV0C(0),
fEvPlaneTPC(0),
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
fHistoNClsE1(0),
fHistoNClsE2(0),
fHistoNClsE3(0),
fHistoNCells(0),
fHistoEperCell(0),
fHistoCalCell(0),
fHistoTimeEMC(0),
fHistoTimeEMCcorr(0),
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
fTPCnsigEta0(0),
fTPCnsigEta1(0),
fTPCnsigEta2(0),
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
fHistNsigEop_Most(0),
fHistNsigEop_Semi(0),
fHistNsigEop_Peri(0),
fHistEop(0),
fHistMcEopEle(0),
fHistMcEopHad(0),
fM20(0),
fM02(0),
fM20EovP(0),
fM02EovP(0),
fEleCanTPCNpts(0),
fEleCanTPCNCls(0),
fEleCanITSNCls(0),
fEleCanITShit(0),
fEleCanSPD1(0),
fEleCanSPD2(0),
fEleCanSPDBoth(0),
fEleCanSPDOr(0),
fITShitPhi(0),
fInvmassULS(0),
fInvmassLS(0),
fInvmassULS_MCtrue(0),
fInvmassPi0Dalitz(0),
fMCcheckMother(0),
fMCneutral(0),
fSparseElectron(0),
fvalueElectron(0)
{
    // Constructor
    
    fvalueElectron = new Double_t[9];
    // Define input and output slots here
    // Input slot #0 works with a TChain
    DefineInput(0, TChain::Class());
    // Output slot #0 id reserved by the base class for AOD
    // Output slot #1 writes into a TH1 container
    DefineOutput(1, TList::Class());
}
//________________________________________________________________________
AliAnalysisTaskHFEemcQA::AliAnalysisTaskHFEemcQA()
: AliAnalysisTaskSE("DefaultTask_HfeEMCQA"),
fVevent(0),
fESD(0),
fAOD(0),
fMCheader(0),
fpidResponse(0),
fEMCALGeo(0),
fFlagSparse(kFALSE),
fUseTender(kTRUE),
fEMCEG1(kFALSE),
fEMCEG2(kFALSE),
fDCalDG1(kFALSE),
fDCalDG2(kFALSE),
fTracks_tender(0),
fCaloClusters_tender(0),
fMCparticle(0),
fMCarray(0),
fMultSelection(0),
fTriggersInfo(0),
fThresholdEG2(89),
fThresholdEG1(140),
fFlagClsTypeEMC(kTRUE),
fcentMim(0),
fcentMax(0),
fCentralityEstimator("V0M"),
fFlagClsTypeDCAL(kTRUE),
fOutputList(0),
fNevents(0),
fCent(0),
fMult(0),
fEvPlaneV0(0),
fEvPlaneV0A(0),
fEvPlaneV0C(0),
fEvPlaneTPC(0),
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
fHistoNClsE1(0),
fHistoNClsE2(0),
fHistoNClsE3(0),
fHistoNCells(0),
fHistoEperCell(0),
fHistoCalCell(0),
fHistoTimeEMC(0),
fHistoTimeEMCcorr(0),
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
fTPCnsigEta0(0),
fTPCnsigEta1(0),
fTPCnsigEta2(0),
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
fHistNsigEop_Most(0),
fHistNsigEop_Semi(0),
fHistNsigEop_Peri(0),
fHistEop(0),
fHistMcEopEle(0),
fHistMcEopHad(0),
fM20(0),
fM02(0),
fM20EovP(0),
fM02EovP(0),
fEleCanTPCNpts(0),
fEleCanTPCNCls(0),
fEleCanITSNCls(0),
fEleCanITShit(0),
fEleCanSPD1(0),
fEleCanSPD2(0),
fEleCanSPDBoth(0),
fEleCanSPDOr(0),
fITShitPhi(0),
fInvmassULS(0),
fInvmassLS(0),
fInvmassULS_MCtrue(0),
fInvmassPi0Dalitz(0),
fMCcheckMother(0),
fMCneutral(0),
fSparseElectron(0),
fvalueElectron(0)
{
    //Default constructor
    
    fvalueElectron = new Double_t[9];
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
AliAnalysisTaskHFEemcQA::~AliAnalysisTaskHFEemcQA()
{
    //Destructor
    delete fOutputList;
    delete fTracks_tender;
    delete fCaloClusters_tender;
    delete fSparseElectron;
    delete []fvalueElectron;
}
//________________________________________________________________________
void AliAnalysisTaskHFEemcQA::UserCreateOutputObjects()
{
    // Create histograms
    // Called once
    AliDebug(3, "Creating Output Objects");
    
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
    
    //AliAODMCHeader *mcHeader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    
    //fEMCALGeo =  AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1_DCAL_8");
    //fEMCALGeo =  AliEMCALGeometry::GetInstance();
    
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
    
    fEvPlaneV0 = new TH2F("fEvPlaneV0","V0 EP",100,0,100,100,0,TMath::Pi());
    fOutputList->Add(fEvPlaneV0);
    
    fEvPlaneV0A = new TH2F("fEvPlaneV0A","V0A EP",100,0,100,100,0,TMath::Pi());
    fOutputList->Add(fEvPlaneV0A);
    
    fEvPlaneV0C = new TH2F("fEvPlaneV0C","V0C EP",100,0,100,100,0,TMath::Pi());
    fOutputList->Add(fEvPlaneV0C);
    
    fEvPlaneTPC = new TH2F("fEvPlaneTPC","TPC EP",100,0,100,100,0,TMath::Pi());
    fOutputList->Add(fEvPlaneTPC);
    
    fVtxZ = new TH1F("fVtxZ","Z vertex position;Vtx_{z};counts",500,-25,25);
    fOutputList->Add(fVtxZ);
    
    fVtxY = new TH1F("fVtxY","Y vertex position;Vtx_{y};counts",500,-25,25);
    fOutputList->Add(fVtxY);
    
    fVtxX = new TH1F("fVtxX","X vertex position;Vtx_{x};counts",500,-25,25);
    fOutputList->Add(fVtxX);
    
    fTrigMulti = new TH2F("fTrigMulti","Multiplicity distribution for different triggers; Trigger type; multiplicity",11,-1,10,2000,0,2000);
    fOutputList->Add(fTrigMulti);
    
    fHistClustE = new TH1F("fHistClustE", "EMCAL cluster energy distribution; Cluster E;counts", 500, 0.0, 50.0);
    fOutputList->Add(fHistClustE);
    
    fHistNonLinClustE = new TH1F("fHistNonLinClustE", "Nonlinearity corrected EMCAL cluster energy distribution; Cluster E;counts", 500, 0.0, 50.0);
    fOutputList->Add(fHistNonLinClustE);
    
    fHistClustEcent = new TH2F("fHistClustEcent", "EMCAL cluster energy distribution vs. centrality; centrality; Cluster E", 100,0,100,500, 0.0, 50.0);
    fOutputList->Add(fHistClustEcent);
    
    fEMCClsEtaPhi = new TH2F("fEMCClsEtaPhi","EMCAL cluster #eta and #phi distribution;#eta;#phi",100,-0.9,0.9,200,0,6.3);
    fOutputList->Add(fEMCClsEtaPhi);
    
    fHistClustEEG1 = new TH1F("fHistClustEEG1", "EMCAL cluster energy distribution for trigger cluster; Cluster E;counts", 500, 0.0, 50.0);
    fOutputList->Add(fHistClustEEG1);
    
    fHistClustEEG1cent = new TH2F("fHistClustEEG1cent", "EMCAL cluster energy distribution vs. centrality for trigger cluster; centrality; Cluster E", 100,0,100,500, 0.0, 50.0);
    fOutputList->Add(fHistClustEEG1cent);
    
    fHistClustEEG2 = new TH1F("fHistClustEEG2", "EMCAL cluster energy distribution; Cluster E;counts", 500, 0.0, 50.0);
    fOutputList->Add(fHistClustEEG2);
    
    fHistClustEEG2cent = new TH2F("fHistClustEEG2cent", "EMCAL cluster energy distribution vs. centrality for trigger cluster; centrality; Cluster E", 100,0,100,500, 0.0, 50.0);
    fOutputList->Add(fHistClustEEG2cent);
    
    fEMCClsEtaPhiEG1 = new TH2F("fEMCClsEtaPhiEG1","EMCAL cluster #eta and #phi distribution;#eta;#phi",100,-0.9,0.9,200,0,6.3);
    fOutputList->Add(fEMCClsEtaPhiEG1);
    
    fEMCClsEtaPhiEG2 = new TH2F("fEMCClsEtaPhiEG2","EMCAL cluster #eta and #phi distribution;#eta;#phi",100,-0.9,0.9,200,0,6.3);
    fOutputList->Add(fEMCClsEtaPhiEG2);
    
    fHistoNCls = new TH1F("fHistoNCls","No of EMCAL cluster in the event;N^{EMC}_{cls};counts",150,0,150);
    fOutputList->Add(fHistoNCls);
    
    fHistoNClsE1 = new TH1F("fHistoNClsE1","No of EMCAL cluster in the event (E>0.1 GeV);N^{EMC}_{cls};counts",150,0,150);
    fOutputList->Add(fHistoNClsE1);
    
    fHistoNClsE2 = new TH1F("fHistoNClsE2","No of EMCAL cluster in the event (E>0.2 GeV);N^{EMC}_{cls};counts",150,0,150);
    fOutputList->Add(fHistoNClsE2);
    
    fHistoNClsE3 = new TH1F("fHistoNClsE3","No of EMCAL cluster in the event (E>0.5 GeV);N^{EMC}_{cls};counts",150,0,150);
    fOutputList->Add(fHistoNClsE3);
    
    fHistoNCells = new TH2F("fHistoNCells","No of EMCAL cells in a cluster;Cluster E;N^{EMC}_{cells}",500,0,50,30,0,30);
    fOutputList->Add(fHistoNCells);
    
    fHistoEperCell = new TH2F("fHistoEperCell","E/cell;Cluster E;E/cell",400,0,40,300,0,30);
    fOutputList->Add(fHistoEperCell);
    
    fHistoCalCell = new TH2F("fHistoCalCell","Energy of EMCAL cells;cell ID;E (GeV)",15000,-0.5,14999.5,150,0,30);
    fOutputList->Add(fHistoCalCell);
    
    fHistoTimeEMC = new TH2F("fHistoTimeEMC","EMCAL Time;E (GeV); t(ns)",500,0,50,1800,-900,900);
    fOutputList->Add(fHistoTimeEMC);
    
    //fHistoTimeEMCcorr = new TH2F("fHistoTimeEMCcorr","EMCAL Time (tender);E (GeV); t(ns)",480,2,50,20000,-200,200);
    //fOutputList->Add(fHistoTimeEMCcorr);
    
 /*   Int_t bincal[6]=      {15,  480,  8000, 40,  140, 160}; // trigger;pT;nSigma;eop;m20;m02;sqrtm02m20;eID;nSigma_Pi;cent
    Double_t xmincal[6]={-0.5,    2,  -160,  0, -0.7,   0};
    Double_t xmaxcal[6]={14.5,   50 ,  160, 40,  0.7, 6.4};
    fHistoTimeEMCcorr = new THnSparseD("fHistoTimeEMCcorr","EMCAL Time (tender);SM;E (GeV); t(ns); Ncell; phi; eta",6,bincal,xmincal,xmaxcal);
    fOutputList->Add(fHistoTimeEMCcorr);
 */
    
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
    
    fTPCnsigEta0 = new TH2F("fTPCnsigEta0","TPC Nsigma vs. Eta pT > 2 GeV/c;#eta;#sigma_{TPC-dE/dx}",40,-1,1,200,-10,10);
    fOutputList->Add(fTPCnsigEta0);

    fTPCnsigEta1 = new TH2F("fTPCnsigEta1","TPC Nsigma vs. Eta pT > 3 GeV/c;#eta;#sigma_{TPC-dE/dx}",40,-1,1,200,-10,10);
    fOutputList->Add(fTPCnsigEta1);
    
    fTPCnsigEta2 = new TH2F("fTPCnsigEta2","TPC Nsigma vs. Eta pT > 5 GeV/c;#eta;#sigma_{TPC-dE/dx}",40,-1,1,200,-10,10);
    fOutputList->Add(fTPCnsigEta2);
    
    fHistPtMatch = new TH1F("fHistPtMatch", "p_{T} distribution of tracks matched to EMCAL;p_{T} (GeV/c);counts",500, 0.0, 50.0);
    fOutputList->Add(fHistPtMatch);
    
    fEMCTrkMatch = new TH2F("fEMCTrkMatch","Distance of EMCAL cluster to its closest track;#phi;z",100,-0.3,0.3,100,-0.3,0.3);
    fOutputList->Add(fEMCTrkMatch);
    
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
    
    fHistNsigEop_Most = new TH2F("fHistNsigEop_Most", "E/p distribution (0-10%);p_{T} (GeV/c);E/p", 40, 0.0, 2.0, 200, -10, 10);
    fOutputList->Add(fHistNsigEop_Most);
    
    fHistNsigEop_Semi = new TH2F("fHistNsigEop_Semi", "E/p distribution (20-40%);p_{T} (GeV/c);E/p", 40, 0.0, 2.0, 200, -10, 10);
    fOutputList->Add(fHistNsigEop_Semi);
    
    fHistNsigEop_Peri = new TH2F("fHistNsigEop_Peri", "E/p distribution (60-80%);p_{T} (GeV/c);E/p", 40, 0.0, 2.0, 200, -10, 10);
    fOutputList->Add(fHistNsigEop_Peri);
    
    fM20 = new TH2F ("fM20","M20 vs pt distribution",500,0,50,200,0,2);
    fOutputList->Add(fM20);
    
    fM02 = new TH2F ("fM02","M02 vs pt distribution",500,0,50,200,0,2);
    fOutputList->Add(fM02);
    
    fM20EovP = new TH2F ("fM20EovP","M20 vs E/p distribution;E/p;M20",40,0,2,200,0,2);
    fOutputList->Add(fM20EovP);
    
    fM02EovP = new TH2F ("fM02EovP","M02 vs E/p distribution;E/p;M02",40,0,2,200,0,2);
    fOutputList->Add(fM02EovP);
    
    fEleCanTPCNpts = new TH2F("fEleCanTPCNpts","TPC Npoints used for dE/dx for electron candidates;p_{T} (GeV/c);N points",200,0,20,200,0,200);
    fOutputList->Add(fEleCanTPCNpts);
    
    fEleCanTPCNCls = new TH2F("fEleCanTPCNCls","TPC N clusters for electron candidates;p_{T} (GeV/c);N TPC clusters",200,0,20,171,-0.5,170.5);
    fOutputList->Add(fEleCanTPCNCls);
    
    fEleCanITSNCls = new TH2F("fEleCanITSNCls","ITS N clusters for electron candidates;p_{T} (GeV/c);N ITS clusters",200,0,20,8,-0.5,7.5);
    fOutputList->Add(fEleCanITSNCls);
    
    fEleCanITShit = new TH1F("fEleCanITShit","ITS hit map;ITS layer;counts",7,-0.5,6.5);
    fOutputList->Add(fEleCanITShit);
    
    fEleCanSPD1 = new TH2F("fEleCanSPD1","Hit on SPD layer 1;p_{T} (GeV/c);Hit",200,0,20,1,0,1);
    fOutputList->Add(fEleCanSPD1);
    
    fEleCanSPD2 = new TH2F("fEleCanSPD2","Hit on SPD layer 2;p_{T} (GeV/c);Hit",200,0,20,1,0,1);
    fOutputList->Add(fEleCanSPD2);
    
    fEleCanSPDBoth = new TH2F("fEleCanSPDBoth","Tracks with hits on both SPD layer;p_{T} (GeV/c);Hit",200,0,20,1,0,1);
    fOutputList->Add(fEleCanSPDBoth);
    
    fEleCanSPDOr = new TH2F("fEleCanSPDOr","Tracks with hits on both SPD layer;p_{T} (GeV/c);Hit",200,0,20,1,0,1);
    fOutputList->Add(fEleCanSPDOr);
    
    fITShitPhi = new TH2F("fITShitPhi","ITS Hit in #phi",4,-0.5,3.5,200,0,6.3);
    fOutputList->Add(fITShitPhi);
    
    fInvmassLS = new TH1F("fInvmassLS", "Invmass of LS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 500,0,1.0);
    fOutputList->Add(fInvmassLS);
    
    fInvmassULS = new TH1F("fInvmassULS", "Invmass of ULS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 500,0,1.0);
    fOutputList->Add(fInvmassULS);
    
    fInvmassULS_MCtrue = new TH2F("fInvmassULS_MCtrue", "Invmass of ULS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 6,-0.5,5.5,1000,0,1.0);
    fOutputList->Add(fInvmassULS_MCtrue);
    
    /*
     Int_t binsDal[6] =      {3,300,200,200,3,100};
     Double_t mimDal[6] = {-0.5,0,0,0,-0.5,-5};
     Double_t maxDal[6] = {2.5,0.3,40.0,40.0,2.5,5};
     fInvmassPi0Dalitz = new THnSparseD("Pi0DalitzMC","Inv mass Dal;feed;mass;epT;pi0pT;prim;eta",6,binsDal,mimDal,maxDal);
     fOutputList->Add(fInvmassPi0Dalitz);
     */
    fMCcheckMother = new TH2F("fMCcheckMother", "Mother MC PDG", 1000,-0.5,999.5,50,0,50);
    fOutputList->Add(fMCcheckMother);
    
    fMCneutral = new TH2F("fMCneutral","pi0 and eta pT from Hijing and enhance",6,-0.5,5.5,500,0,50);
    fOutputList->Add(fMCneutral);
    
    if(fFlagSparse){
    Int_t bins[9]=      {8, 280, 160, 40, 200, 200,    3, 100,  10}; // trigger;pT;nSigma;eop;m20;m02;sqrtm02m20;eID;nSigma_Pi;cent
    Double_t xmin[9]={-0.5,   2,  -8,   0,   0,   0, -0.5,  -5,   0};
    Double_t xmax[9]={ 7.5,  30,   8,   2,   2,   2,  2.5,  15, 100};
    fSparseElectron = new THnSparseD ("Electron","Electron;trigger;pT;nSigma;eop;m20;m02;eID;nSigma_Pi;cent;",9,bins,xmin,xmax);
    fOutputList->Add(fSparseElectron);
    }
    PostData(1,fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskHFEemcQA::UserExec(Option_t *)
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
        if(IsAODanalysis()) fTracks_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("tracks"));
       // if(!IsAODanalysis()) fTracks_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("Tracks"));
        fCaloClusters_tender = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject("caloClusters"));
    }
    
    ////////////////////
    //cuts initialised//
    ////////////////////
    AliESDtrackCuts* esdTrackCutsH = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE);
    esdTrackCutsH->SetMaxDCAToVertexXY(2.4);
    esdTrackCutsH->SetMaxDCAToVertexZ(3.2);
    esdTrackCutsH->SetDCAToVertex2D(kTRUE);
    esdTrackCutsH->SetMinNClustersTPC(80);
    esdTrackCutsH->SetMinNClustersITS(3);
    esdTrackCutsH->SetRequireTPCRefit(kTRUE);
    esdTrackCutsH->SetRequireITSRefit(kTRUE);
    esdTrackCutsH->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);
    esdTrackCutsH->SetAcceptKinkDaughters(kFALSE);
    esdTrackCutsH->SetMaxChi2PerClusterITS(6); //test.....
    
    fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
    if (fAOD) {
        // printf("fAOD available\n");
        //return;
    }
    if(fAOD)fMCarray = dynamic_cast<TClonesArray*>(fAOD->FindListObject(AliAODMCParticle::StdBranchName()));
    
    fMCheader = dynamic_cast<AliAODMCHeader*>(fAOD->GetList()->FindObject(AliAODMCHeader::StdBranchName()));
    
    if(fMCarray)CheckMCgen(fMCheader);
    
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
    
    /////////////////
    //trigger check//
    /////////////////
    TString firedTrigger;
    TString TriggerEG1("EG1");
    TString TriggerEG2("EG2");
    TString TriggerDG1("DG1");
    TString TriggerDG2("DG2");
    fVevent->GetFiredTriggerClasses();
    if(fAOD) firedTrigger = fAOD->GetFiredTriggerClasses();
    
    Bool_t EG1tr = kFALSE;
    Bool_t EG2tr = kFALSE;
    if(firedTrigger.Contains(TriggerEG1))EG1tr = kTRUE;
    if(firedTrigger.Contains(TriggerEG2))EG2tr = kTRUE;
    
    if(fEMCEG1){if(!firedTrigger.Contains(TriggerEG1))return;}
    if(fEMCEG2){if(!firedTrigger.Contains(TriggerEG2))return;}
    if(fDCalDG1){if(!firedTrigger.Contains(TriggerDG1))return;}
    if(fDCalDG2){if(!firedTrigger.Contains(TriggerDG2))return;}
    
    Int_t trigger = -1;
    if (fAOD){
        AliAODHeader *header = dynamic_cast<AliAODHeader*>(fAOD->GetHeader());
        if(!header) AliFatal("Not a standard AOD");
        Double_t multiplicity = header->GetRefMultiplicity();
        
        fTrigMulti->Fill(-0.5, multiplicity);
        if(evSelMask & AliVEvent::kAny) fTrigMulti->Fill(0.5, multiplicity);
        if(evSelMask & AliVEvent::kMB) fTrigMulti->Fill(1.5, multiplicity);
        if(evSelMask & AliVEvent::kINT7) fTrigMulti->Fill(2.5, multiplicity);
        if(evSelMask & AliVEvent::kINT8) fTrigMulti->Fill(3.5, multiplicity);
        if(evSelMask & AliVEvent::kEMC1) fTrigMulti->Fill(4.5, multiplicity);
        if(evSelMask & AliVEvent::kEMC7) fTrigMulti->Fill(5.5, multiplicity);
        if(evSelMask & AliVEvent::kEMC8) fTrigMulti->Fill(6.5, multiplicity);
        if(evSelMask & AliVEvent::kEMCEJE) fTrigMulti->Fill(7.5, multiplicity);
        if(evSelMask & AliVEvent::kEMCEGA) fTrigMulti->Fill(8.5, multiplicity);
        if(evSelMask & AliVEvent::kEMCEGA & EG2tr) fTrigMulti->Fill(9.5, multiplicity);
        
        if(evSelMask & AliVEvent::kMB) trigger =0;
        if(evSelMask & AliVEvent::kINT7) trigger =1;
        if(evSelMask & AliVEvent::kINT8) trigger =2;
        if(evSelMask & AliVEvent::kEMC1) trigger =3;
        if(evSelMask & AliVEvent::kEMC7) trigger =4;
        if(evSelMask & AliVEvent::kEMC8) trigger =5;
        if(evSelMask & AliVEvent::kEMCEJE) trigger =6;
        if(evSelMask & AliVEvent::kEMCEGA) trigger =7;
    }
    
    ////////////////////
    //event selection///
    ////////////////////
    if(TMath::Abs(Zvertex)>10.0)return;
    fNevents->Fill(2); //events after z vtx cut
    fCent->Fill(centrality); //centrality dist.
    
    ///////////////////
    // event plane
    /////////////////////
    
    Double_t epV0A = 0, epV0C = 0, epV0 = 0, epTPC = 0, qxV0A = 0, qyV0A = 0, qxV0C = 0, qyV0C = 0, qxV0 = 0, qyV0 = 0,  qxTPC = 0, qyTPC = 0;
    TVector2 *qTPC = 0x0;
    TVector2 qVectorfortrack;
    
    // V0
    if(fESD){
        epV0 = TVector2::Phi_0_2pi(fESD->GetEventplane()->GetEventplane("V0",fESD,2));
        epV0A = TVector2::Phi_0_2pi(fESD->GetEventplane()->GetEventplane("V0A",fESD,2));
        epV0C = TVector2::Phi_0_2pi(fESD->GetEventplane()->GetEventplane("V0C",fESD,2));
    }
    
    if (fAOD){
        epV0 = TVector2::Phi_0_2pi(fAOD->GetEventplane()->CalculateVZEROEventPlane(fAOD,10,2,qxV0,qyV0));
        epV0A = TVector2::Phi_0_2pi(fAOD->GetEventplane()->CalculateVZEROEventPlane(fAOD,8,2,qxV0A,qyV0A));
        epV0C = TVector2::Phi_0_2pi(fAOD->GetEventplane()->CalculateVZEROEventPlane(fAOD,9,2,qxV0C,qyV0C));
    }
    
    if(epV0 > TMath::Pi())  epV0 = epV0 - TMath::Pi();
    if(epV0A > TMath::Pi()) epV0A = epV0A - TMath::Pi();
    if(epV0C > TMath::Pi()) epV0C = epV0C - TMath::Pi();
    
    //TPC
    
    if (fAOD){
        AliEventplane* Eventplane=  fAOD->GetEventplane();
        if (Eventplane && Eventplane->GetQVector()) {
            qTPC = Eventplane->GetQVector();
            qxTPC = qTPC->X();
            qyTPC = qTPC->Y();
            
            qVectorfortrack.Set(qxTPC,qyTPC);
            epTPC = TVector2::Phi_0_2pi(qVectorfortrack.Phi())/2.;
        }
    }
    
    else{
        AliEventplane* Eventplane =  fESD->GetEventplane();
        if (Eventplane && Eventplane->GetQVector()) {
            qTPC = Eventplane->GetQVector();
            qxTPC = qTPC->X();
            qyTPC = qTPC->Y();
            
            qVectorfortrack.Set(qxTPC,qyTPC);
            epTPC = TVector2::Phi_0_2pi(qVectorfortrack.Phi())/2.;
        }
    }
    
    
    fEvPlaneV0->Fill(centrality,epV0);  // cent. vs V0 EP
    fEvPlaneV0A->Fill(centrality,epV0A);// cent. vs V0A EP
    fEvPlaneV0C->Fill(centrality,epV0C);// cent. vs V0C EP
    fEvPlaneTPC->Fill(centrality,epTPC);// cent. vs TPC EP
    
    /////////////////////////////
    //EMCAL cluster information//
    /////////////////////////////
    Int_t Nclust = -999;
    if(!fUseTender) Nclust = fVevent->GetNumberOfCaloClusters();
    if(fUseTender) Nclust = fCaloClusters_tender->GetEntries();
    
    int NclustAll= 0;
    int NclustE1 = 0; //# of clust E>0.1
    int NclustE2 = 0; //# of clust E>0.2
    int NclustE3 = 0; //# of clust E>0.5
    
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
            Double_t clustE = clust->E();
            if(clustE < 0.3) continue;
            
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
            
            fHistClustE->Fill(clustE);
            fHistNonLinClustE->Fill(clust->GetNonLinCorrEnergy());
            
            if(centrality>-1)fHistClustEcent->Fill(centrality,clustE);
            fEMCClsEtaPhi->Fill(emceta,emcphi);
            fHistoNCells->Fill(clustE,clust->GetNCells());
            Double_t EperCell = -999.9;
            if(clust->GetNCells()>0)EperCell = clustE/clust->GetNCells();
            fHistoEperCell->Fill(clustE,EperCell);
            
            Float_t tof = clust->GetTOF()*1e+9; // ns
            fHistoTimeEMC->Fill(clustE,tof);
            
            //-----Plots for EMC trigger
            Bool_t hasfiredEG1=0;
            Bool_t hasfiredEG2=0;
            FindPatches(hasfiredEG1,hasfiredEG2,emceta,emcphi);
            if(hasfiredEG1){
                fHistClustEEG1->Fill(clustE);
                if(centrality>-1)fHistClustEEG1cent->Fill(centrality,clustE);
                fEMCClsEtaPhiEG1->Fill(emceta,emcphi);
            }
            if(hasfiredEG2){
                fHistClustEEG2->Fill(clustE);
                if(centrality>-1)fHistClustEEG2cent->Fill(centrality,clustE);
                fEMCClsEtaPhiEG2->Fill(emceta,emcphi);
            }
            
            NclustAll++;
            if(clustE>0.1)NclustE1++;
            if(clustE>0.2)NclustE2++;
            if(clustE>0.5)NclustE3++;
        }
    }
    
    fHistoNCls->Fill(NclustAll);
    fHistoNClsE1->Fill(NclustE1);
    fHistoNClsE2->Fill(NclustE2);
    fHistoNClsE3->Fill(NclustE3);
    
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
    
    if(!fEMCALGeo)fEMCALGeo  = AliEMCALGeometry::GetInstance(); // not work w.o. Tender
    //cout << "fEMCALGeo= " << fEMCALGeo << endl;
    
    ///////////////
    //Track loop///
    ///////////////
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
        
        ////////////////////
        //Apply track cuts//
        ////////////////////
        if(fAOD)
            if(!atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue; //mimimum cuts
        
        if(fESD)
            if(!esdTrackCutsH->AcceptTrack(etrack))continue;
        
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
        Double_t d0z0[2]={-999,-999}, cov[3];
        Double_t DCAxyCut = 2.4, DCAzCut = 3.2;
        if(fAOD){
            if(atrack->GetTPCNcls() < 80) continue;
            if(atrack->GetITSNcls() < 3) continue;
            if((!(atrack->GetStatus()&AliESDtrack::kITSrefit)|| (!(atrack->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
            if(!(atrack->HasPointOnITSLayer(0) || atrack->HasPointOnITSLayer(1))) continue;
            
            double phiMatchIts = atrack->Phi();
            if(atrack->HasPointOnITSLayer(0))fITShitPhi->Fill(0.0,phiMatchIts);
            if(atrack->HasPointOnITSLayer(1))fITShitPhi->Fill(1.0,phiMatchIts);
            
            if(atrack->PropagateToDCA(pVtx, fVevent->GetMagneticField(), 20., d0z0, cov))
                if(TMath::Abs(d0z0[0]) > DCAxyCut || TMath::Abs(d0z0[1]) > DCAzCut) continue;
            //To be done : Add cuts to apply Chi2PerITSCls < 6 and N shared Cls ITS < 4
        }
        
        ///////////////////////
        // Get MC information//
        ///////////////////////
        Int_t ilabel = track->GetLabel();
        Int_t pdg = -999;
        Int_t pidM = -1;
        Double_t pid_ele = 0.0;
        if(ilabel>0 && fMCarray)
        {
            fMCparticle = (AliAODMCParticle*) fMCarray->At(ilabel);
            Int_t pdg = fMCparticle->GetPdgCode();
            if(TMath::Abs(pdg)==11)pid_ele = 1.0;
            Int_t ilabelM = -1;
            if(pid_ele==1.0)FindMother(fMCparticle, ilabelM, pidM);
            
            if(ilabelM>0)
            {
                AliAODMCParticle* fMCparticleM = (AliAODMCParticle*) fMCarray->At(ilabelM);
                if(pidM==22) // from pi0 & eta
                {
                    AliAODMCParticle* fMCparticleM = (AliAODMCParticle*) fMCarray->At(ilabelM);
                    FindMother(fMCparticleM, ilabelM, pidM);
                }
                Double_t pTmom = fMCparticleM->Pt();
                //fMCcheckMother->Fill(abs(pidM),pTmom);
            }
        }
        
        ////////////////////
        //Track properties//
        ///////////////////
        Double_t dEdx =-999, fTPCnSigma=-999, fTPCnSigma_Pi=-999;
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
        
        if(pid_ele==1.0)
            fTPCnsigMcEle->Fill(TrkP,fTPCnSigma);
        else
            fTPCnsigMcHad->Fill(TrkP,fTPCnSigma);
        
        fTPCnsig_Pi->Fill(TrkP,fTPCnSigma_Pi);
        
        if(TrkPt>2.0)fTPCnsigEta0->Fill(TrkEta,fTPCnSigma);
        if(TrkPt>3.0)fTPCnsigEta1->Fill(TrkEta,fTPCnSigma);
        if(TrkPt>5.0)fTPCnsigEta2->Fill(TrkEta,fTPCnSigma);
        
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
        int iSM = -1;
        for(int icl=0; icl<NcellsInCluster; icl++)
        {
            int icell = clustMatch->GetCellAbsId(icl);
            if(fEMCALGeo)
                iSM = fEMCALGeo->GetSuperModuleNumber(icell);
        }
        
        Double_t emcphi = -999, emceta=-999;
        fClsTypeEMC = kFALSE; fClsTypeDCAL = kFALSE;
        if(clustMatch && clustMatch->IsEMCAL())
        {
            // fEMCTrkMatch->Fill(clustMatch->GetTrackDx(),clustMatch->GetTrackDz());
            //if(TMath::Abs(clustMatch->GetTrackDx())>0.05 || TMath::Abs(clustMatch->GetTrackDz())>0.05) continue;
            
            Double_t fPhiDiff = -999, fEtaDiff = -999;
            GetTrkClsEtaPhiDiff(track, clustMatch, fPhiDiff, fEtaDiff);
            fEMCTrkMatch->Fill(fPhiDiff,fEtaDiff);
            
            if(TMath::Abs(fPhiDiff) > 0.05 || TMath::Abs(fEtaDiff)> 0.05) continue;
            
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
            
   /*         Float_t tof = clustMatch->GetTOF()*1e+9; // ns
            Double_t caloinfo[6];
            caloinfo[0] = iSM;
            caloinfo[1] = clustMatchE;
            caloinfo[2] = tof;
            caloinfo[3] = clustMatch->GetNCells();
            caloinfo[4] = emceta;
            caloinfo[5] = emcphi;
            if(clustMatchE>2.0 && fUseTender)fHistoTimeEMCcorr->Fill(caloinfo);
    */
            //EMCAL EID info
            Double_t eop = -1.0;
            Double_t m02 = -99999,m20 = -99999,sqm02m20=-99999.0;
            if(track->P()>0)eop = clustMatchE/track->P();
            m02 =clustMatch->GetM02();
            m20 =clustMatch->GetM20();
            //sqm02m20 = sqrt(pow(m02,2)+pow(m20,2));
            
            if(track->Pt()>3.0){
                fHistdEdxEop->Fill(eop,dEdx);
                fHistNsigEop->Fill(eop,fTPCnSigma);
                if(centrality>=0 && centrality<=10)fHistNsigEop_Most->Fill(eop,dEdx);
                if(centrality>=20 && centrality<=40)fHistNsigEop_Semi->Fill(eop,dEdx);
                if(centrality>=60 && centrality<=80)fHistNsigEop_Peri->Fill(eop,dEdx);
                fM20EovP->Fill(eop,clustMatch->GetM20());
                fM02EovP->Fill(eop,clustMatch->GetM02());
            }
            fM20->Fill(track->Pt(),clustMatch->GetM20());
            fM02->Fill(track->Pt(),clustMatch->GetM02());
            
            //EID THnsparse
            fvalueElectron[0] = trigger;
            fvalueElectron[1] = track->Pt();
            fvalueElectron[2] = fTPCnSigma;
            fvalueElectron[3] = eop;
            fvalueElectron[4] = m20;
            fvalueElectron[5] = m02;
            fvalueElectron[6] = pid_ele;
            fvalueElectron[7] = fTPCnSigma_Pi;
            fvalueElectron[8] = centrality;
            
            if(fFlagSparse && track->Pt()>2.0){
                fSparseElectron->Fill(fvalueElectron);
            }
            
            Bool_t fFlagNonHFE=kFALSE;
            ////////////////////////////////////////////////
            //Track properties of EMCAL electron cadidates//
            ////////////////////////////////////////////////
            
            if((fTPCnSigma > -1 && fTPCnSigma < 3) && (m20 > 0.01 && m20 < 0.45))fHistEop->Fill(track->Pt(),eop);
            if(pid_ele==1.0)
                fHistMcEopEle->Fill(track->Pt(),eop);
            else
                fHistMcEopHad->Fill(track->Pt(),eop);
            
            if(fTPCnSigma > -1 && fTPCnSigma < 3 && eop>0.9 && eop<1.2 && m02 > 0.006 && m02 < 0.35){ //rough cuts
                //-----Identify Non-HFE
                SelectPhotonicElectron(iTracks,track,fFlagNonHFE,pidM);
                
                fEleCanTPCNpts->Fill(track->Pt(),track->GetTPCsignalN());
                fEleCanTPCNCls->Fill(track->Pt(),track->GetTPCNcls());
                
                
                Int_t fITSncls=0;
                for(Int_t l=0;l<6;l++) {
                    if(TESTBIT(track->GetITSClusterMap(),l)) {
                        fEleCanITShit->Fill(l);
                        if(l==0) fEleCanSPD1->Fill(track->Pt(),0.5);
                        if(l==1) fEleCanSPD2->Fill(track->Pt(),0.5);
                        if(l==0 && l==1) fEleCanSPDBoth->Fill(track->Pt(),0.5);
                        if(l==0 || l==1) fEleCanSPDOr->Fill(track->Pt(),0.5);
                        fITSncls++;
                    }
                }
                fEleCanITSNCls->Fill(track->Pt(),fITSncls++);
            }
        }
    } //track loop
    
    PostData(1, fOutputList);
}
//________________________________________________________________________
void AliAnalysisTaskHFEemcQA::SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec, Int_t iMC)
{
    ///////////////////////////////////////////
    //////Non-HFE - Invariant mass method//////
    ///////////////////////////////////////////
    
    AliESDtrackCuts* esdTrackCutsAsso = AliESDtrackCuts::GetStandardTPCOnlyTrackCuts();
    esdTrackCutsAsso->SetAcceptKinkDaughters(kFALSE);
    esdTrackCutsAsso->SetRequireTPCRefit(kTRUE);
    esdTrackCutsAsso->SetRequireITSRefit(kTRUE);
    esdTrackCutsAsso->SetEtaRange(-0.9,0.9);
    esdTrackCutsAsso->SetMaxChi2PerClusterTPC(4);
    esdTrackCutsAsso->SetMinNClustersTPC(70);
    esdTrackCutsAsso->SetMaxDCAToVertexZ(3.2);
    esdTrackCutsAsso->SetMaxDCAToVertexXY(2.4);
    esdTrackCutsAsso->SetDCAToVertex2D(kTRUE);
    
    Bool_t flagPhotonicElec = kFALSE;
    
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
        
        Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;
        Double_t ptAsso=-999., nsigma=-999.0, mass=-999., width = -999;
        Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;
        
        nsigma = fpidResponse->NumberOfSigmasTPC(Assotrack, AliPID::kElectron);
        ptAsso = Assotrack->Pt();
        Int_t chargeAsso = Assotrack->Charge();
        Int_t charge = track->Charge();
        if(charge>0) fPDGe1 = -11;
        if(chargeAsso>0) fPDGe2 = -11;
        if(charge == chargeAsso) fFlagLS = kTRUE;
        if(charge != chargeAsso) fFlagULS = kTRUE;
        
        //------track cuts applied
        if(fAOD) {
            if(!aAssotrack->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
            if(aAssotrack->GetTPCNcls() < 70) continue;
            if((!(aAssotrack->GetStatus()&AliESDtrack::kITSrefit)|| (!(aAssotrack->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
        }
        else{
            if(!esdTrackCutsAsso->AcceptTrack(eAssotrack)) continue;
        }
        
        //-------loose cut on partner electron
        if(ptAsso <0.2) continue;
        if(aAssotrack->Eta()<-0.9 || aAssotrack->Eta()>0.9) continue;
        if(nsigma < -3 || nsigma > 3) continue;
        
        //-------define KFParticle to get mass
        AliKFParticle::SetField(fVevent->GetMagneticField());
        AliKFParticle ge1 = AliKFParticle(*track, fPDGe1);
        AliKFParticle ge2 = AliKFParticle(*Assotrack, fPDGe2);
        AliKFParticle recg(ge1, ge2);
        
        if(recg.GetNDF()<1) continue;
        Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
        if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;
        
        //-------Get mass
        Int_t MassCorrect;
        MassCorrect = recg.GetMass(mass,width);
        
        if(fFlagLS)
            if(track->Pt()>1) fInvmassLS->Fill(mass);
        if(fFlagULS)
            if(track->Pt()>1) fInvmassULS->Fill(mass);
        
        if(iMC>0)
        {
            Int_t iMCbin = -999;
            if(iMC == 111)
            {
                iMCbin = 1;
            }
            else if(iMC == 221)
            {
                iMCbin = 2;
            }
            else
            {
                iMCbin = -999;
            }
            
            //if(fFlagULS && track->Pt()>1.5 && iMCbin!=-999)fInvmassULS_MCtrue->Fill(iMCbin,mass);
        }
        
        if(mass<0.2 && fFlagULS && !flagPhotonicElec)
            flagPhotonicElec = kTRUE; //Tag Non-HFE (random mass cut, not optimised)
    }
    fFlagPhotonicElec = flagPhotonicElec;
}

//________________________________________________________________________
void AliAnalysisTaskHFEemcQA::GetTrkClsEtaPhiDiff(AliVTrack *t, AliVCluster *v, Double_t &phidiff, Double_t &etadiff)
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
void AliAnalysisTaskHFEemcQA::FindMother(AliAODMCParticle* part, Int_t &label, Int_t &pid)
{
    // Find mother in case of MC
    
    if(part->GetMother()>-1)
    {
        label = part->GetMother();
        AliAODMCParticle *partM = (AliAODMCParticle*)fMCarray->At(label);
        pid = partM->GetPdgCode();
    }
    else
    {
        pid = -1;
    }
}

void AliAnalysisTaskHFEemcQA::CheckMCgen(AliAODMCHeader* fMCheader)
{
    TList *lh=fMCheader->GetCocktailHeaders();
    Int_t NpureMC = 0;
    Int_t NpureMCproc = 0;
    if(lh)
    {
        for(int igene=0; igene<lh->GetEntries(); igene++)
        {
            AliGenEventHeader* gh=(AliGenEventHeader*)lh->At(igene);
            if(gh)
            {
                if(igene==0)NpureMC = gh->NProduced();  // generate by PYTHIA or HIJING
                NpureMCproc += gh->NProduced();
            }
        }
    }
    
    //for(int imc=0; imc<fMCarray->GetEntries(); imc++)
    for(int imc=0; imc<NpureMCproc; imc++)
    {
        Bool_t iEnhance = kFALSE;
        if(imc>=NpureMC)iEnhance = kTRUE;
        Int_t iHijing = 1;  // select particles from Hijing or PYTHIA
        
        fMCparticle = (AliAODMCParticle*) fMCarray->At(imc);
        Int_t pdgGen = TMath::Abs(fMCparticle->GetPdgCode());
        
        
        Int_t phyprim = 0;
        if(fMCparticle->IsPrimary())phyprim = 1;
        
        Double_t PtPi0 = 0.0;
        Int_t pdgMom = -99;
        
        Int_t labelMpi = -1;
        FindMother(fMCparticle,labelMpi,pdgMom);
        if(pdgMom==-1 && iEnhance)iHijing = 0;  // select particles orogonally from enhance
        
        if(iHijing ==0)
        {
            if(pdgGen==411 || pdgGen==421 || pdgGen==413 || pdgGen==423 || pdgGen==431 || pdgGen==433)fMCcheckMother->Fill(pdgGen,fMCparticle->Pt());
            if(pdgGen==511 || pdgGen==521 || pdgGen==513 || pdgGen==523 || pdgGen==531 || pdgGen==533)fMCcheckMother->Fill(pdgGen,fMCparticle->Pt());
            if(pdgGen==111)fMCcheckMother->Fill(pdgGen,fMCparticle->Pt());
            if(pdgGen==221)fMCcheckMother->Fill(pdgGen,fMCparticle->Pt());
        }
        
        if(pdgGen==111 || pdgGen==221)
        {
            PtPi0 = fMCparticle->Pt();
            if(pdgGen==111 && iHijing==0)fMCneutral->Fill(0.0,fMCparticle->Pt());
            if(pdgGen==111 && iHijing==1)fMCneutral->Fill(1.0,fMCparticle->Pt());
            if(pdgGen==221 && iHijing==0)fMCneutral->Fill(2.0,fMCparticle->Pt());
            if(pdgGen==221 && iHijing==1)fMCneutral->Fill(3.0,fMCparticle->Pt());
            
            Int_t Ndecay = fMCparticle->GetNDaughters();
            if(Ndecay==3)
            {
                Int_t firstCh = fMCparticle->GetDaughter(0);
                Int_t lastCh = fMCparticle->GetDaughter(1);
                //cout << "firstCh = " << firstCh << " ; lastCh = " << lastCh << endl;
                
                AliAODMCParticle* fMCpar0 = (AliAODMCParticle*) fMCarray->At(firstCh);
                AliAODMCParticle* fMCpar1 = (AliAODMCParticle*) fMCarray->At(firstCh+1);
                AliAODMCParticle* fMCpar2 = (AliAODMCParticle*) fMCarray->At(firstCh+2);
                
                Int_t pdgCh0 = fMCpar0->GetPdgCode();
                Int_t pdgCh1 = fMCpar1->GetPdgCode();
                Int_t pdgCh2 = fMCpar2->GetPdgCode();
                //cout << "pdg = " << pdgCh0  << " ;  " << pdgCh1 << " ; " << pdgCh2 << endl;
                
                if(pdgCh0==22 && TMath::Abs(pdgCh1)==11 && TMath::Abs(pdgCh2)==11)
                {
                    TLorentzVector chele1;
                    chele1.SetPxPyPzE(fMCpar1->Px(),fMCpar1->Py(),fMCpar1->Pz(),fMCpar1->E());
                    TLorentzVector chele2;
                    chele2.SetPxPyPzE(fMCpar2->Px(),fMCpar2->Py(),fMCpar2->Pz(),fMCpar2->E());
                    
                    TLorentzVector Sumchele;
                    Sumchele = chele1 + chele2;
                    if(fMCpar1->Pt()>0.5)
                    {
                        if(pdgGen==111 && iHijing==0)fInvmassULS_MCtrue->Fill(1,Sumchele.M());  
                        if(pdgGen==221 && iHijing==0)fInvmassULS_MCtrue->Fill(2,Sumchele.M());  
                        if(pdgGen==111 && iHijing==1) fInvmassULS_MCtrue->Fill(3,Sumchele.M());  
                        if(pdgGen==221 && iHijing==1) fInvmassULS_MCtrue->Fill(4,Sumchele.M());  
                        
                        /* 
                         if(pdgGen==111)
                         {
                         Double_t par[6];
                         if(iHijing==0)cout << "pi0 Dalitz ; "<< imc << " ; pdgMpi = " << pdgMom << " ; Enhance = " << iEnhance << " ; pi0 eta = " << fMCparticle->Eta() << endl;
                         par[0] = iHijing;
                         par[1] = Sumchele.M();
                         par[2] = fMCpar1->Pt();
                         par[3] = PtPi0;
                         par[4] = phyprim;
                         par[5] = fMCparticle->Eta();
                         fInvmassPi0Dalitz->Fill(par);
                         }
                         */
                    }
                } 
                
            } 
        }
        
    }
    
    return;
}

//________________________________________________________________________
void AliAnalysisTaskHFEemcQA::FindPatches(Bool_t &hasfiredEG1,Bool_t &hasfiredEG2,Double_t emceta,Double_t emcphi)
{
    //Find trigger patches
    
    fTriggersInfo = dynamic_cast <TClonesArray*>(InputEvent()->FindListObject("EmcalTriggers"));
    if(!fTriggersInfo) return;
    Int_t nPatch = fTriggersInfo->GetEntries();;
    AliEMCALTriggerPatchInfo* patch=0;
    for( int iPatch = 0; iPatch < nPatch; iPatch++ ){
        patch = (AliEMCALTriggerPatchInfo*)fTriggersInfo->At( iPatch );
        if(patch->GetADCAmp()<fThresholdEG2) continue;
        if(patch->GetEtaMin()>emceta) continue;
        if(patch->GetEtaMax()<emceta) continue;
        if(patch->GetPhiMin()>emcphi) continue;
        if(patch->GetPhiMax()<emcphi) continue;
        if(patch->GetADCAmp()>fThresholdEG2)  hasfiredEG2=1;
        if(patch->GetADCAmp()>fThresholdEG1)  hasfiredEG1=1;
    }
}

//________________________________________________________________________
void AliAnalysisTaskHFEemcQA::Terminate(Option_t *)
{
    // Draw result to the screen
    // Called once at the end of the query
    
    fOutputList = dynamic_cast<TList*> (GetOutputData(1));
    if (!fOutputList) {
        printf("ERROR: Output list not available\n");
        return;
    }
}
