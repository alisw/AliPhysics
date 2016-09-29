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
//  Task for Heavy Flavour Electron-Hadron DeltaPhi Correlation       //
//  Non-Photonic Electron identified with Invariant mass              //
//  analysis methos in function  SelectPhotonicElectron               //
//  DeltaPhi calculated in function  ElectronHadCorrel                //
//                                                                    //
//  Author: Deepa Thomas (University of Texas at Austin)              //
//                                                                    //
////////////////////////////////////////////////////////////////////////

#include "TChain.h"
#include "TTree.h"
#include "TH2F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "THnSparse.h"
#include "TLorentzVector.h"
#include "TString.h"
#include "TFile.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDHandler.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"

#include "AliAnalysisTaskEHCorrel.h"
#include "TGeoGlobalMagField.h"
#include "AliLog.h"
#include "AliAnalysisTaskSE.h"
#include "TRefArray.h"
#include "TVector.h"


//#include "AliEventPoolManager.h"
#include "AliESDInputHandler.h"
#include "AliAODInputHandler.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliESDtrackCuts.h"
#include "AliCentralitySelectionTask.h"
#include "AliMultSelection.h"
#include "AliESDCaloCluster.h"
#include "AliAODCaloCluster.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"
#include "AliESDCaloTrigger.h"
#include "AliEMCALRecoUtils.h"
#include "AliEMCALGeometry.h"
#include "AliGeomManager.h"
#include "stdio.h"
#include "TGeoManager.h"
#include "iostream"
#include "fstream"

#include "AliCentrality.h"
#include "AliMagF.h"

#include "AliPID.h"
#include "AliPIDResponse.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"
#include "AliVEvent.h"
#include "AliStack.h"
#include "AliMCEvent.h"
#include "TProfile.h"
#include "AliESDVZERO.h"
#include "AliAODVZERO.h"
#include "TVector3.h"
#include "TRandom2.h"

ClassImp(AliAnalysisTaskEHCorrel)
ClassImp(AliehDPhiBasicParticle)
  //________________________________________________________________________
  AliAnalysisTaskEHCorrel::AliAnalysisTaskEHCorrel(const char *name)
: AliAnalysisTaskSE(name),
  fVevent(0),
  fAOD(0),
  fpVtx(0),
  fpidResponse(0),
  fMultSelection(0),
  fCentrality(-1),
  fCentralityMin(0),
  fCentralityMax(20),
  fEMCEG1(kFALSE),
  fEMCEG2(kFALSE),
  fFlagClsTypeEMC(kTRUE),
  fFlagClsTypeDCAL(kTRUE),
  fTPCnSigma(-999.0),
  fTPCnSigmaMin(-1),
  fTPCnSigmaMax(3),
  fM02Min(0.01),
  fM02Max(0.35),
  fM20Min(0),
  fM20Max(2),
  fEovPMin(0.9),
  fEovPMax(1.2),
  fTPCNClsHad(80),
  fInvmassCut(0.1),
  fTPCnSigmaHadMin(-10),
  fTPCnSigmaHadMax(-3.5),
  fHadCutCase(1),
 // fTracksCloneMix(0),
 // fPool(0),
 // fPoolMgr(0),
  fOutputList(0),
  fNevents(0),
  fVtxZ(0),
  fVtxX(0),
  fVtxY(0),
  fCentralityNoPass(0),
  fCentralityPass(0),
  fHistClustE(0),
  fEMCClsEtaPhi(0),
  fNegTrkIDPt(0),
  fTrkPt(0),
  fTrketa(0),
  fTrkphi(0),
  fdEdx(0),
  fTPCnsig(0),
  fHistPtMatch(0),
  fEMCTrkMatch(0),
  fEMCTrkPt(0),
  fEMCTrketa(0),
  fEMCTrkphi(0),
  fEMCTPCnsig(0),
  fClsEAftMatch(0),
  fClsEtaPhiAftMatch(0),
  fHistNsigEop(0),
  fM20EovP(0),
  fM02EovP(0),
  fHistEop(0),
  fM20(0),
  fM02(0),
  fHistEop_AftEID(0),
  fInclsElecPt(0),
  fHadEop(0),
  fHadPt_AftEID(0),
  fHadronPhiPt(0),
  fHadronPhi(0),
  fHadronPhiTPChalf(0),
  fHadronPt(0),
  fInvmassLS(0),
  fInvmassULS(0),
  fInvmassLSPt(0),
  fInvmassULSPt(0),
  fNoMixedEvents(0),
  fMixStatCent(0),
  fMixStatVtxZ(0),
  fSprsHadHCorrl(0),
  fSprsInclusiveEHCorrl(0),
  fSprsLSEHCorrl(0),
  fSprsULSEHCorrl(0)
{
  //Named constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  // DefineOutput(1, TH1I::Class());
  DefineOutput(1, TList::Class());
  //  DefineOutput(3, TTree::Class());
}

//________________________________________________________________________
AliAnalysisTaskEHCorrel::AliAnalysisTaskEHCorrel()
  : AliAnalysisTaskSE("DefaultAnalysis_AliAnalysisTaskEHCorrel"),
  fVevent(0),
  fAOD(0),
  fpVtx(0),
  fpidResponse(0),
  fMultSelection(0),
  fCentrality(-1),
  fCentralityMin(0),
  fCentralityMax(20),
  fEMCEG1(kFALSE),
  fEMCEG2(kFALSE),
  fFlagClsTypeEMC(kTRUE),
  fFlagClsTypeDCAL(kTRUE),
  fTPCnSigma(-999.0),
  fTPCnSigmaMin(-1),
  fTPCnSigmaMax(3),
  fM02Min(0.01),
  fM02Max(0.35),
  fM20Min(0),
  fM20Max(2),
  fEovPMin(0.9),
  fEovPMax(1.2),
  fTPCNClsHad(80),
  fInvmassCut(0.1),
  fTPCnSigmaHadMin(-10),
  fTPCnSigmaHadMax(-3.5),
  fHadCutCase(1),
//  fTracksCloneMix(0),
//  fPool(0),
//  fPoolMgr(0),
  fOutputList(0),
  fNevents(0),
  fVtxZ(0),
  fVtxX(0),
  fVtxY(0),
  fCentralityNoPass(0),
  fCentralityPass(0),
  fHistClustE(0),
  fEMCClsEtaPhi(0),
  fNegTrkIDPt(0),
  fTrkPt(0),
  fTrketa(0),
  fTrkphi(0),
  fdEdx(0),
  fTPCnsig(0),
  fHistPtMatch(0),
  fEMCTrkMatch(0),
  fEMCTrkPt(0),
  fEMCTrketa(0),
  fEMCTrkphi(0),
  fEMCTPCnsig(0),
  fClsEAftMatch(0),
  fClsEtaPhiAftMatch(0),
  fHistNsigEop(0),
  fM20EovP(0),
  fM02EovP(0),
  fHistEop(0),
  fM20(0),
  fM02(0),
  fHistEop_AftEID(0),
  fInclsElecPt(0),
  fHadEop(0),
  fHadPt_AftEID(0),
  fHadronPhiPt(0),
  fHadronPhi(0),
  fHadronPhiTPChalf(0),
  fHadronPt(0),
  fInvmassLS(0),
  fInvmassULS(0),
  fInvmassLSPt(0),
  fInvmassULSPt(0),
  fNoMixedEvents(0),
  fMixStatCent(0),
  fMixStatVtxZ(0),
  fSprsHadHCorrl(0),
  fSprsInclusiveEHCorrl(0),
  fSprsLSEHCorrl(0),
  fSprsULSEHCorrl(0)
{
  //Default constructor
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  // DefineOutput(1, TH1I::Class());
  DefineOutput(1, TList::Class());
  //DefineOutput(3, TTree::Class());
}
//_________________________________________
AliAnalysisTaskEHCorrel::~AliAnalysisTaskEHCorrel()
{
  //Destructor

  delete fOutputList;
  delete fSprsHadHCorrl;
  delete fSprsInclusiveEHCorrl;
  delete   fSprsLSEHCorrl;
  delete fSprsULSEHCorrl;
 // delete fTracksCloneMix;
}
//_________________________________________
void AliAnalysisTaskEHCorrel::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  AliDebug(3, "Creating Output Objects");

  Double_t pi = TMath::Pi();
    
    ////////////////////////
    //Initiale mixed event//
    ////////////////////////
    Int_t trackDepth = 2000;
    Int_t poolsize   = 1000;
    Int_t nCentralityBins  = 6;
    Int_t nCentralityBinsSC  = 6;
    Int_t nZvtxBins  = 4;
    Double_t CentralityBins[7];
    Double_t vertexBins[5] = {-10,-5,0,5,10};
    
    if(fCentralityMax <= 20)
    {
        CentralityBins[0] = 0;
    CentralityBins[1] =2;
   CentralityBins[2] = 4;
    CentralityBins[3] =6;
    CentralityBins[4] =10;
    CentralityBins[5] =15;
    CentralityBins[6] =20;
}
    
    if(fCentralityMax <= 50)
    {
        CentralityBins[0] = 20;
        CentralityBins[1] = 25;
        CentralityBins[2] = 30;
        CentralityBins[3] = 35;
        CentralityBins[4] = 40;
        CentralityBins[5] = 45;
        CentralityBins[6] = 50;
    }
    
  //      fPoolMgr = new AliEventPoolManager(poolsize, trackDepth, nCentralityBins, (Double_t*) CentralityBins, nZvtxBins, (Double_t*) vertexBins);

  //  fTracksCloneMix = new TClonesArray();
  //  fTracksCloneMix->SetOwner(kTRUE);

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

  fVtxZ = new TH1F("fVtxZ","Z vertex position;Vtx_{z};counts",1000,-50,50);
  fOutputList->Add(fVtxZ);

  fVtxY = new TH1F("fVtxY","Y vertex position;Vtx_{y};counts",1000,-50,50);
  fOutputList->Add(fVtxY);

  fVtxX = new TH1F("fVtxX","X vertex position;Vtx_{x};counts",1000,-50,50);
  fOutputList->Add(fVtxX);

  fCentralityPass = new TH1F("fCentralityPass", "Centrality Pass", 101, -1, 100);
  fOutputList->Add(fCentralityPass);

  fCentralityNoPass = new TH1F("fCentralityNoPass", "Centrality No Pass", 101, -1, 100);
  fOutputList->Add(fCentralityNoPass);

  fHistClustE = new TH1F("fHistClustE", "EMCAL cluster energy distribution; Cluster E;counts", 500, 0.0, 50.0);
  fOutputList->Add(fHistClustE);

  fEMCClsEtaPhi = new TH2F("fEMCClsEtaPhi","EMCAL cluster #eta and #phi distribution;#eta;#phi",100,-0.9,0.9,200,0,6.3);
  fOutputList->Add(fEMCClsEtaPhi);

  fNegTrkIDPt = new TH1F("fNegTrkIDPt", "p_{T} distribution of tracks with negative track id;p_{T} (GeV/c);counts", 500, 0.0, 50.0);
  fOutputList->Add(fNegTrkIDPt);

  fTrkPt = new TH1F("fTrkPt","p_{T} distribution of all tracks;p_{T} (GeV/c);counts",1000,0,100);
  fOutputList->Add(fTrkPt);

  fTrketa = new TH1F("fTrketa","All Track #eta distribution;#eta;counts",100,-1.5,1.5);
  fOutputList->Add(fTrketa);

  fTrkphi = new TH1F("fTrkphi","All Track #phi distribution;#phi;counts",100,0,2*pi);
  fOutputList->Add(fTrkphi);

  fdEdx = new TH2F("fdEdx","All Track dE/dx distribution;p (GeV/c);dE/dx",200,0,20,500,0,160);
  fOutputList->Add(fdEdx);

  fTPCnsig = new TH2F("fTPCnsig","All Track TPC Nsigma distribution;p (GeV/c);#sigma_{TPC-dE/dx}",1000,0,50,200,-10,10);
  fOutputList->Add(fTPCnsig);

  fHistPtMatch = new TH1F("fHistPtMatch", "p_{T} distribution of tracks matched to EMCAL;p_{T} (GeV/c);counts",1000, 0.0, 100.0);
  fOutputList->Add(fHistPtMatch);

  fEMCTrkMatch = new TH2F("fEMCTrkMatch","Distance of EMCAL cluster to its closest track;#phi;z",100,-0.3,0.3,100,-0.3,0.3);
  fOutputList->Add(fEMCTrkMatch);

  fEMCTrkPt = new TH1F("fEMCTrkPt","p_{T} distribution of tracks with EMCAL cluster;p_{T} (GeV/c);counts",1000,0,100);
  fOutputList->Add(fEMCTrkPt);

  fEMCTrketa = new TH1F("fEMCTrketa","#eta distribution of tracks matched to EMCAL;#eta;counts",100,-1.5,1.5);
  fOutputList->Add(fEMCTrketa);

  fEMCTrkphi = new TH1F("fEMCTrkphi","#phi distribution of tracks matched to EMCAL;#phi;counts",100,0,2*pi);
  fOutputList->Add(fEMCTrkphi);


  fEMCTPCnsig = new TH2F("fEMCTPCnsig","TPC Nsigma distribution of tracks matched to EMCAL;p (GeV/c);#sigma_{TPC-dE/dx}",1000,0,50,200,-10,10);
  fOutputList->Add(fEMCTPCnsig);

  fClsEAftMatch = new TH1F("fClsEAftMatch", "EMCAL cluster energy distribution after track matching; Cluster E;counts", 500, 0.0, 50.0);
  fOutputList->Add(fClsEAftMatch);

  fClsEtaPhiAftMatch = new TH2F("fClsEtaPhiAftMatch","EMCAL cluster #eta and #phi distribution after track matching;#eta;#phi",100,-0.9,0.9,200,0,6.3);
  fOutputList->Add(fClsEtaPhiAftMatch);

  fHistNsigEop = new TH2F ("fHistNsigEop", "E/p vs TPC nsig",60, 0.0, 3.0, 200, -10,10);
  fOutputList->Add(fHistNsigEop);

  fM20 = new TH2F ("fM20","M20 vs pt distribution",200,0,20,400,0,2);
  fOutputList->Add(fM20);

  fM02 = new TH2F ("fM02","M02 vs pt distribution",200,0,20,400,0,2);
  fOutputList->Add(fM02);

  fM20EovP = new TH2F ("fM20EovP","M20 vs E/p distribution",400,0,3,400,0,2);
  fOutputList->Add(fM20EovP);

  fM02EovP = new TH2F ("fM02EovP","M02 vs E/p distribution",400,0,3,400,0,2);
  fOutputList->Add(fM02EovP);

  fHistEop = new TH2F("fHistEop", "E/p distribution;p_{T} (GeV/c);E/p", 200,0,20,60, 0.0, 3.0);
  fOutputList->Add(fHistEop);

  fHistEop_AftEID = new TH2F("fHistEop_AftEID", "E/p distribution after nsig, SS cuts;p_{T} (GeV/c);E/p", 200,0,20,60, 0.0, 3.0);
  fOutputList->Add(fHistEop_AftEID);

  fInclsElecPt = new TH1F("fInclsElecPt","p_{T} distribution of inclusive electrons;p_{T} (GeV/c);counts",500,0,50);
  fOutputList->Add(fInclsElecPt);

  fHadEop = new TH2F("fHadEop", "E/p distribution for hadrons;p_{T} (GeV/c);E/p", 200,0,20,60, 0.0, 3.0);
  fOutputList->Add(fHadEop);
  fHadPt_AftEID = new TH1F("fHadPt_AftEID","p_{T} distribution of hadrons after Eid cuts;p_{T} (GeV/c);counts",500,0,50);
  fOutputList->Add(fHadPt_AftEID);


  fHadronPhiPt = new TH2F("fHadronPhiPt", "Hadron phi vs pt; hadron phi; pt (GeV/c)",1000,0,2*pi,500,0,50);
  fOutputList->Add(fHadronPhiPt);

  fHadronPhi = new TH1F("fHadronPhi", "Hadron phi",1000,0,2*pi);
  fOutputList->Add(fHadronPhi);

  fHadronPhiTPChalf = new TH1F("fHadronPhiTPChalf", "Hadron phi for 0<eta<0.9",1000,0,2*pi);
  fOutputList->Add(fHadronPhiTPChalf);

  fHadronPt = new TH1F("fHadronPt","hadron pt distribution",500,0,50);
  fOutputList->Add(fHadronPt);

  fInvmassLS = new TH1F("fInvmassLS", "Inv mass of LS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 1000,0,1.0);
  fOutputList->Add(fInvmassLS);

  fInvmassULS = new TH1F("fInvmassULS", "Inv mass of ULS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 1000,0,1.0);
  fOutputList->Add(fInvmassULS);

  fInvmassLSPt = new TH2F("fInvmassLSPt", "Inv mass of LS (e,e) vs pT; p_{T}(GeV/c); mass(GeV/c^2); counts;", 500,0,50,1000,0,1.0);
  fOutputList->Add(fInvmassLSPt);

  fInvmassULSPt = new TH2F("fInvmassULSPt", "Inv mass of ULS (e,e) vs pT; p_{T}(GeV/c); mass(GeV/c^2); counts;", 500,0,50,1000,0,1.0);
  fOutputList->Add(fInvmassULSPt);
    
    fNoMixedEvents = new TH1F("fNoMixedEvents","No of mixing events",1,-0.5,0.5);
    fOutputList->Add(fNoMixedEvents);
    
    fMixStatCent = new TH2F("fMixStatCent","Mix event stats for centrality binning;Nevent in pool;Centrality",500,0,500,nCentralityBins,(Double_t*)CentralityBins);
    fOutputList->Add(fMixStatCent);
    
    fMixStatVtxZ = new TH2F("fMixStatVtxZ","Mix event stats for Zvtx binning;Nevent in pool;Vtx_{z}",500,0,500,nZvtxBins,(Double_t*)vertexBins);
    fOutputList->Add(fMixStatVtxZ);

  //------THnsparse------
  Int_t bin[4] = {50,50,64,100}; //ptElec, ptHad,Dphi, Deta
  Double_t xmin[4] = {0,0,-TMath::Pi()/2,-1.8};
  Double_t xmax[4] = {50,50,(3*TMath::Pi())/2,1.8};

  fSprsHadHCorrl = new THnSparseD("fSprsHadHCorrl","Sparse for Dphi and Deta hadrons;p_{T}^{e};p_{T}^{h};#Delta#varphi;#Delta#eta;",4,bin,xmin,xmax);
  fOutputList->Add(fSprsHadHCorrl);

  fSprsInclusiveEHCorrl = new THnSparseD("fSprsInclusiveEHCorrl","Sparse for Dphi and Deta with Inclusive electron;p_{T}^{e};p_{T}^{h};#Delta#varphi;#Delta#eta;",4,bin,xmin,xmax);
  fOutputList->Add(fSprsInclusiveEHCorrl);

  fSprsLSEHCorrl = new THnSparseD("fSprsLSEHCorrl","Sparse for Dphi and Deta with LS Non-HFE electron;p_{T}^{e};p_{T}^{h};#Delta#varphi;#Delta#eta;",4,bin,xmin,xmax);
  fOutputList->Add(fSprsLSEHCorrl);

  fSprsULSEHCorrl = new THnSparseD("fSprsULSEHCorrl","Sparse for Dphi and Deta with ULS Non-HFE electron;p_{T}^{e};p_{T}^{h};#Delta#varphi;#Delta#eta;",4,bin,xmin,xmax);
  fOutputList->Add(fSprsULSEHCorrl);

  PostData(1,fOutputList);
}
//_________________________________________
void AliAnalysisTaskEHCorrel::UserExec(Option_t*)
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

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (fAOD) {
    // printf("fAOD available\n");
    //return;
  }
  fpVtx = fVevent->GetPrimaryVertex();
    
  ///////////////////
  //PID initialised//
  ///////////////////
  fpidResponse = fInputHandler->GetPIDResponse();

  /////////////////
  // Centrality ///
  /////////////////
  Bool_t pass = kFALSE;
  if(fCentralityMin > -0.5){
    CheckCentrality(fAOD,pass);
    if(!pass)return;
  }

  /////////////////
  //trigger check//
  /////////////////
  TString firedTrigger;
  TString TriggerEG1("EG1");
  TString TriggerEG2("EG2");
  if(fAOD) firedTrigger = fAOD->GetFiredTriggerClasses();

  if(fEMCEG1){if(!firedTrigger.Contains(TriggerEG1))return;}
  if(fEMCEG2){if(!firedTrigger.Contains(TriggerEG2))return;}

  ////////////////////
  //event selection///
  ////////////////////
  if(!PassEventSelect(fVevent)) return;

  //////////////////////
  //EMcal cluster info//
  //////////////////////
  EMCalClusterInfo();

  ///////////////
  //Track loop///
  ///////////////
  Int_t ntracks = fVevent->GetNumberOfTracks();
  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {
    AliVParticle* Vtrack = 0x0;
    Vtrack  = fVevent->GetTrack(iTracks);
    if (!Vtrack) {
      printf("ERROR: Could not receive track %d\n", iTracks);
      continue;
    }
    AliVTrack *track = dynamic_cast<AliVTrack*>(Vtrack);
    AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(Vtrack);

    ////////////////////
    //Apply track cuts//
    ////////////////////
    if(!PassTrackCuts(atrack)) continue;

    Double_t TrkPhi=-999, TrkPt=-999, TrkEta=-999, TrkP = -999;
    TrkPhi = track->Phi();
    TrkPt = track->Pt();
    TrkEta = track->Eta();
    TrkP = track->P();

    ///////////////////////////
    //Track matching to EMCAL//
    //////////////////////////
    Int_t EMCalIndex = -1;
    EMCalIndex = track->GetEMCALcluster();
    if(EMCalIndex < 0) continue;
    fHistPtMatch->Fill(TrkPt);

    AliVCluster *clustMatch=0x0;
    clustMatch = (AliVCluster*)fVevent->GetCaloCluster(EMCalIndex);
    Double_t emcphi = -999, emceta=-999;
    Bool_t fClsTypeEMC = kFALSE, fClsTypeDCAL = kFALSE;
    if(clustMatch && clustMatch->IsEMCAL())
    {

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
      fEMCTrketa->Fill(TrkEta);
      fEMCTrkphi->Fill(TrkPhi);
      fEMCTPCnsig->Fill(TrkP,fTPCnSigma);
      Double_t clustMatchE = clustMatch->E();
      fClsEAftMatch->Fill(clustMatchE);
      fClsEtaPhiAftMatch->Fill(emceta,emcphi);

      //////////////////
      //Apply EID cuts//
      //////////////////
      Bool_t fHadTrack = kFALSE, fElectTrack = kFALSE;
      fElectTrack = PassEIDCuts(track, clustMatch, fHadTrack);

      ///////////////////
      //H-H Correlation//
      ///////////////////
      if(fHadTrack)
      {
        fHadPt_AftEID->Fill(TrkPt);
        ElectronHadCorrel(iTracks, track, fSprsHadHCorrl);
      }

      if(!fElectTrack) continue;
      fInclsElecPt->Fill(TrkPt);

      ///////////////////
      //E-H Correlation//
      ///////////////////
      HadronInfo(iTracks);

      //Inclusive E-H correl
      ElectronHadCorrel(iTracks, track, fSprsInclusiveEHCorrl);
      //MixedEvent(track);

      ////////////////////
      //NonHFE selection//
      ////////////////////
      Bool_t fFlagPhotonicElec = kFALSE;
      SelectNonHFElectron(iTracks,track,fFlagPhotonicElec);
    }
  }
    

    /////////////////////////
    //Fill Mixed event pool//
    /////////////////////////
/*    Double_t pVtxZ = fpVtx->GetZ();
    fPool = fPoolMgr->GetEventPool(fCentrality, pVtxZ); // Get the buffer associated with the current centrality and z-vtx
    if (!fPool)
    {
        AliFatal(Form("No pool found for centrality = %f, zVtx = %f", centvalue1, pVtxZ));
        return;
    }
    
    CloneAndReduceTrackList();
    fPool->UpdatePool(fTracksCloneMix);
    fTracksCloneMix->Clear("C");
*/
  PostData(1, fOutputList);
}
//_________________________________________
void AliAnalysisTaskEHCorrel::ElectronHadCorrel(Int_t itrack, AliVTrack *track, THnSparse *SparseEHCorrl)
{
  //Construct Deta Phi between electrons and hadrons

  Double_t fvalueDphi[4] = {-999,999,-999,-999}; //ptElec, ptHad,Dphi, Deta
  Double_t pi = TMath::Pi();

  for(Int_t ktracks = 0; ktracks<fVevent->GetNumberOfTracks(); ktracks++){
    if(ktracks == itrack) continue; //do not select the same electron
    AliVParticle* VtrackHad = fVevent->GetTrack(ktracks);
    if (!VtrackHad) {
      printf("ERROR: Could not receive track %d\n", ktracks);
      continue;
    }
    AliVTrack *trackHad = dynamic_cast<AliVTrack*>(VtrackHad);
    if(!trackHad) continue;
    AliAODTrack *atrackHad = dynamic_cast<AliAODTrack*>(VtrackHad);
    if(!atrackHad) continue;

    Double_t ptHad= -999;
    Double_t ptEle = -999;
    Double_t phiEle = -999, phiHad = -999, Dphi = -999;
    Double_t etaEle = -999, etaHad = -999, Deta = -999;

    ptHad = trackHad->Pt();
    ptEle = track->Pt();
    phiEle = track->Phi();
    phiHad = trackHad->Phi();
    etaEle = track->Eta();
    etaHad = trackHad->Eta();

    if(!PassHadronCuts(atrackHad)) continue; //apply hadron cuts;
    // if(ptHad > ptEle) continue;

    Dphi = phiEle - phiHad;
    if (Dphi > 3*pi/2)
      Dphi = Dphi - 2*pi;
    if (Dphi < -pi/2)
      Dphi = Dphi + 2*pi;

    Deta = etaEle - etaHad;

    fvalueDphi[0] = ptEle;
    fvalueDphi[1] = ptHad;
    fvalueDphi[2] = Dphi;
    fvalueDphi[3] = Deta;
    SparseEHCorrl->Fill(fvalueDphi);
  }
}
//___________________________________________
void AliAnalysisTaskEHCorrel::HadronInfo(Int_t itrack)
{
  //Hadron information

  for(Int_t ktracks = 0; ktracks<fVevent->GetNumberOfTracks(); ktracks++){

    if(ktracks == itrack) continue; //do not select the same electron
    AliVParticle* VtrackHad = fVevent->GetTrack(ktracks);
    if (!VtrackHad) {
      printf("ERROR: Could not receive track %d\n", ktracks);
      continue;
    }
    AliVTrack *trackHad = dynamic_cast<AliVTrack*>(VtrackHad);
    if(!trackHad) continue;
    AliAODTrack *atrackHad = dynamic_cast<AliAODTrack*>(VtrackHad);
    if(!atrackHad) continue;

    if(!PassHadronCuts(atrackHad)) continue; //apply hadron cuts;

    Double_t ptHad= -999, phiHad=-999;
    ptHad = trackHad->Pt();
    phiHad = trackHad->Phi();

    fHadronPhiPt->Fill(phiHad,ptHad);
    fHadronPhi->Fill(phiHad);
    if (trackHad->Eta() >0 && trackHad->Eta() <0.9) fHadronPhiTPChalf->Fill(phiHad);
    fHadronPt->Fill(ptHad);
  }
}
//___________________________________________
Bool_t AliAnalysisTaskEHCorrel::PassHadronCuts(AliAODTrack *HadTrack)
{
  //apply hadron cuts

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
  Double_t d0z0[2]={-999,-999}, cov[3];
  Double_t DCAxyCut = 0.25, DCAzCut = 1;

  if(fHadCutCase == 1)
  {
    if(!HadTrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return kFALSE;
    if((!(HadTrack->GetStatus()&AliESDtrack::kITSrefit)|| (!(HadTrack->GetStatus()&AliESDtrack::kTPCrefit)))) return kFALSE;
  }
  if(fHadCutCase == 2)
  {
    if(!HadTrack->TestFilterMask(AliAODTrack::kTrkTPCOnly)) return kFALSE;
    if((!(HadTrack->GetStatus()&AliESDtrack::kITSrefit)|| (!(HadTrack->GetStatus()&AliESDtrack::kTPCrefit)))) return kFALSE;
  }
  if(fHadCutCase == 3)
  {
    if(!HadTrack->IsHybridGlobalConstrainedGlobal()) return kFALSE;
  }

  if(HadTrack->GetTPCNcls() < fTPCNClsHad) return kFALSE;
  if(HadTrack->Eta()< -0.9 || HadTrack->Eta()>0.9) return kFALSE;
  if(HadTrack->Pt() < 0.3) return kFALSE;
  if(HadTrack->PropagateToDCA(pVtx, fVevent->GetMagneticField(), 20., d0z0, cov))
    if(TMath::Abs(d0z0[0]) > DCAxyCut || TMath::Abs(d0z0[1]) > DCAzCut) return kFALSE;

  return kTRUE;
}
//___________________________________________
Bool_t AliAnalysisTaskEHCorrel::PassEIDCuts(AliVTrack *track, AliVCluster *clust, Bool_t &Hadtrack)
{
  //apply electron identification cuts

  Bool_t hadTrk = kFALSE;
  Double_t eop = -1.0;
  Double_t m02 = -999,m20 = -999;
  Double_t clustE = clust->E();
  Double_t TrkPt = track->Pt();
  if(track->P()>0)eop = clustE/track->P();
  m02 =clust->GetM02();
  m20 =clust->GetM20();

  if(track->Pt()>2.0){
    fHistNsigEop->Fill(eop,fTPCnSigma);
    fM20EovP->Fill(eop,m20);
    fM02EovP->Fill(eop,m02);
  }
  fHistEop->Fill(TrkPt,eop);
  fM20->Fill(TrkPt,m20);
  fM02->Fill(TrkPt,m02);

  //Hadron E/p distribution
  if(fTPCnSigma > fTPCnSigmaHadMin && fTPCnSigma < fTPCnSigmaHadMax)
  {
    if((m02 > fM02Min && m02 < fM02Max) && (m20 > fM20Min && m20 < fM20Max))
    {
      fHadEop->Fill(TrkPt,eop);
      if(eop > fEovPMin && eop < fEovPMax) hadTrk=kTRUE;
    }
  }
  Hadtrack = hadTrk;

  if(fTPCnSigma < fTPCnSigmaMin || fTPCnSigma > fTPCnSigmaMax) return kFALSE;
  if(m02 < fM02Min || m02 > fM02Max) return kFALSE;
  if(m20 < fM20Min || m20 > fM20Max) return kFALSE;

  fHistEop_AftEID->Fill(TrkPt,eop);

  if(eop < fEovPMin || eop > fEovPMax) return kFALSE;

  return kTRUE;
}
//___________________________________________
Bool_t AliAnalysisTaskEHCorrel::PassTrackCuts(AliAODTrack *atrack)
{
  //apply track cuts

  Double_t d0z0[2]={-999,-999}, cov[3];
  Double_t DCAxyCut = 0.25, DCAzCut = 1;
  Double_t dEdx =-999;
  Double_t TrkPhi=-999, TrkPt=-999, TrkEta=-999, TrkP = -999;

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  const AliVVertex *pVtx = fVevent->GetPrimaryVertex();

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
  if(!atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) return kFALSE; //mimimum cuts
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
  if(atrack->GetTPCNcls() < 80) return kFALSE;
  if(atrack->GetITSNcls() < 3) return kFALSE;
  if((!(atrack->GetStatus()&AliESDtrack::kITSrefit)|| (!(atrack->GetStatus()&AliESDtrack::kTPCrefit)))) return kFALSE;
  if(!(atrack->HasPointOnITSLayer(0) || atrack->HasPointOnITSLayer(1))) return kFALSE;

  if(atrack->PropagateToDCA(pVtx, fVevent->GetMagneticField(), 20., d0z0, cov))
    if(TMath::Abs(d0z0[0]) > DCAxyCut || TMath::Abs(d0z0[1]) > DCAzCut) return kFALSE;

  ////////////////////
  //Track properties//
  ////////////////////
  dEdx = atrack->GetTPCsignal();
  fTPCnSigma = fpidResponse->NumberOfSigmasTPC(atrack, AliPID::kElectron);
  TrkPhi = atrack->Phi();
  TrkPt = atrack->Pt();
  TrkEta = atrack->Eta();
  TrkP = atrack->P();

  if(atrack->GetID()<0) fNegTrkIDPt->Fill(TrkPt);
  fTrkPt->Fill(TrkPt);
  fTrketa->Fill(TrkEta);
  fTrkphi->Fill(TrkPhi);
  fdEdx->Fill(TrkP,dEdx);
  fTPCnsig->Fill(TrkP,fTPCnSigma);

  return kTRUE;
}
//___________________________________________
Bool_t AliAnalysisTaskEHCorrel::PassEventSelect(AliVEvent *fVevent)
{
  //event selection cuts

  Int_t ntracks = -999;
  ntracks = fVevent->GetNumberOfTracks();
  if(ntracks < 1) printf("There are %d tracks in this event\n",ntracks);

  fNevents->Fill(0); //all events
  Double_t Zvertex = -100, Xvertex = -100, Yvertex = -100;
    
  Double_t NcontV = fpVtx->GetNContributors();
  if(NcontV<2)return kFALSE;
  fNevents->Fill(1); //events with 2 tracks

  Zvertex = fpVtx->GetZ();
  Yvertex = fpVtx->GetY();
  Xvertex = fpVtx->GetX();
  fVtxZ->Fill(Zvertex);
  fVtxX->Fill(Xvertex);
  fVtxY->Fill(Yvertex);

  if(TMath::Abs(Zvertex)>10.0) return kFALSE;
  fNevents->Fill(2); //events after z vtx cut

  return kTRUE;
}
//___________________________________________
void AliAnalysisTaskEHCorrel::CheckCentrality(AliAODEvent* fAOD, Bool_t &centralitypass)
{
  //check centrality, Run 2

  if(fAOD)fMultSelection = (AliMultSelection * ) fAOD->FindListObject("MultSelection");
  if(!fMultSelection) {
    //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
    AliWarning("AliMultSelection object not found!");
  }else{
    fCentrality = fMultSelection->GetMultiplicityPercentile("V0M", false);
  }

  if ((fCentrality <= fCentralityMin) || (fCentrality > fCentralityMax))
  {
    fCentralityNoPass->Fill(fCentrality);
    //  cout << "--------------Fill no pass-------------------------"<<endl;
    centralitypass = kFALSE;
  }else
  {
    fCentralityPass->Fill(fCentrality);
    //  cout << "--------------Fill pass-------------------------"<<endl;
    centralitypass = kTRUE;
  }
}
//________________________________________________________________________
void AliAnalysisTaskEHCorrel::EMCalClusterInfo()
{
  //EMCAL cluster information

  Int_t Nclust = -999;
  TVector3 clustpos;
  Float_t  emcx[3]; // cluster pos
  Double_t clustE=-999, emcphi = -999, emceta=-999;
  Nclust = fVevent->GetNumberOfCaloClusters();
  for(Int_t icl=0; icl<Nclust; icl++)
  {
    AliVCluster *clust = 0x0;
    clust = fVevent->GetCaloCluster(icl);
    if(!clust)  printf("ERROR: Could not receive cluster matched calibrated from track %d\n", icl);

    if(clust && clust->IsEMCAL())
    {
      clustE = clust->E();
      clust->GetPosition(emcx);
      clustpos.SetXYZ(emcx[0],emcx[1],emcx[2]);
      emcphi = clustpos.Phi();
      emceta = clustpos.Eta();
      if(emcphi < 0) emcphi = emcphi+(2*TMath::Pi()); //TLorentz vector is defined between -pi to pi, so negative phi has to be flipped.

      fHistClustE->Fill(clustE);
      fEMCClsEtaPhi->Fill(emceta,emcphi);
    }
  }
}
//________________________________________________________________________
void AliAnalysisTaskEHCorrel::GetTrkClsEtaPhiDiff(AliVTrack *t, AliVCluster *v, Double_t &phidiff, Double_t &etadiff)
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
void AliAnalysisTaskEHCorrel::SelectNonHFElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec)
{
  //Photonic electron selection

  Bool_t flagPhotonicElec = kFALSE;
  Double_t ptAsso=-999., nsigma=-999.0;
  Bool_t fFlagLS=kFALSE, fFlagULS=kFALSE;

  for(Int_t jTracks = 0; jTracks<fVevent->GetNumberOfTracks(); jTracks++){
    if(jTracks==itrack) continue;

    AliVParticle* VtrackAsso = fVevent->GetTrack(jTracks);
    if (!VtrackAsso) {
      printf("ERROR: Could not receive track %d\n", jTracks);
      continue;
    }

    AliVTrack *trackAsso = dynamic_cast<AliVTrack*>(VtrackAsso);
    if(!trackAsso) continue;

    AliAODTrack *atrackAsso = dynamic_cast<AliAODTrack*>(VtrackAsso);
    if(!atrackAsso) continue;
    if(!atrackAsso->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
    if(atrackAsso->GetTPCNcls() < 70) continue;
    if(!(atrackAsso->GetStatus()&AliESDtrack::kTPCrefit)) continue;
    if(!(atrackAsso->GetStatus()&AliESDtrack::kITSrefit)) continue;

    nsigma = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);
    ptAsso = trackAsso->Pt();
    Int_t chargeAsso = trackAsso->Charge();
    Int_t charge = track->Charge();

    if(ptAsso <0.2) continue;
    if(trackAsso->Eta()<-0.9 || trackAsso->Eta()>0.9) continue;
    if(nsigma < -3 || nsigma > 3) continue;

    Int_t fPDGe1 = 11; Int_t fPDGe2 = 11;
    if(charge>0) fPDGe1 = -11;
    if(chargeAsso>0) fPDGe2 = -11;

    fFlagLS=kFALSE; fFlagULS=kFALSE;
    if(charge == chargeAsso) fFlagLS = kTRUE;
    if(charge != chargeAsso) fFlagULS = kTRUE;

    AliKFParticle::SetField(fVevent->GetMagneticField());

    AliKFParticle ge1 = AliKFParticle(*track, fPDGe1);
    AliKFParticle ge2 = AliKFParticle(*trackAsso, fPDGe2);
    AliKFParticle recg(ge1, ge2);

    if(recg.GetNDF()<1) continue;
    Double_t chi2recg = recg.GetChi2()/recg.GetNDF();
    if(TMath::Sqrt(TMath::Abs(chi2recg))>3.) continue;

    Double_t mass=-999., width = -999;
    Int_t MassCorrect;
    MassCorrect = recg.GetMass(mass,width);

    if(fFlagLS && track->Pt()>1) fInvmassLS->Fill(mass);
    if(fFlagULS && track->Pt()>1) fInvmassULS->Fill(mass);
    if(fFlagLS) fInvmassLSPt->Fill(track->Pt(),mass);
    if(fFlagULS) fInvmassULSPt->Fill(track->Pt(),mass);

    if(fFlagLS && mass<fInvmassCut) ElectronHadCorrel(itrack, track, fSprsLSEHCorrl);
    if(fFlagULS && mass<fInvmassCut) ElectronHadCorrel(itrack, track, fSprsULSEHCorrl);

    if(mass<fInvmassCut && fFlagULS && !flagPhotonicElec){
      flagPhotonicElec = kTRUE;
    }
  }
  fFlagPhotonicElec = flagPhotonicElec;
}
/*//___________________________________________
void  AliAnalysisTaskEHCorrel::MixedEvent(AliVTrack *track)
{
    //Retrive mixed event pool
    Double_t zVtx;
    zVtx = fpVtx->GetZ();
    
    fPool = fPoolMgr->GetEventPool(fCentrality, zVtx); // Get the buffer associated with the current centrality and z-vtx
    if (!fPool)
    {
        AliFatal(Form("No pool found for centrality = %f, zVtx = %f", centvalue1, pVtxZ));
        return;
    }
    //fPool->PrintInfo();
    if (fPool->GetCurrentNEvents() >= 5) // start mixing when 5 events are in the buffer
    {
        Int_t nMix = fPool->GetCurrentNEvents();
        fNoMixedEvents->Fill(0);
        fMixStatCent->Fill(fPool->GetCurrentNEvents(),fCentrality);
        fMixStatVtxZ->Fill(fPool->GetCurrentNEvents(),zVtx);
    }
}
//___________________________________________
void  AliAnalysisTaskEHCorrel::CloneAndReduceTrackList()
{
// clones a track list by using AliehDPhiBasicParticle which uses much less memory (used for event mixing)
    
    for(Int_t ktracks = 0; ktracks<fVevent->GetNumberOfTracks(); ktracks++){
        AliVParticle* Vtrack = fVevent->GetTrack(ktracks);
        if (!Vtrack) {
            printf("ERROR: Could not receive track %d\n", ktracks);
            continue;
        }
        
        AliVTrack *track = dynamic_cast<AliVTrack*>(Vtrack);
        if(!track) continue;
        AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(Vtrack);
        if(!atrack) continue;
        
        if(!PassHadronCuts(atrack)) continue; //apply hadron cuts;
        
        AliVParticle* particle = (AliVParticle*) fVevent->GetTrack(ktracks);
        fTracksCloneMix->Add(new AliehDPhiBasicParticle(particle->Eta(), particle->Phi(), particle->Pt(), particle->Charge()));
    }
}
 */
//___________________________________________
void AliAnalysisTaskEHCorrel::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}



