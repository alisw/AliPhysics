// $Id$
//
// Simulation EMCal task.
//
// Author: Saehanseul Oh

#include "AliAnalysisTaskSOH.h"

#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THnSparse.h>

#include "AliAnalysisManager.h"
#include "AliAnalysisTask.h"
#include "AliEMCALRecoUtils.h"
#include "AliEMCALTrack.h"
#include "AliESDCaloCluster.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliExternalTrackParam.h"
#include "AliLog.h"
#include "AliMCEvent.h"

ClassImp(AliAnalysisTaskSOH)

//________________________________________________________________________
AliAnalysisTaskSOH::AliAnalysisTaskSOH() :
  AliAnalysisTaskSE(), 
  fESD(0), 
  fMC(0), 
  fEsdTrackCuts(0x0),
  fHybridTrackCuts1(0x0),
  fHybridTrackCuts2(0x0),
  fTrackIndices(0x0),
  fClusterIndices(0x0),
  fClusterArray(0x0),
  fOutputList(0x0),        
  fHEventStat(0),      
  fHTrkEffParGenPtEtaPhi(0), 
  fHTrkEffDetGenPtEtaPhi(0), 
  fHTrkEffDetRecPtEtaPhi(0),
  fHTrkEffDetRecFakePtEtaPhi(0), 
  fHTrkEffParGenPt_all(0),
  fHTrkEffDetRecPt_all(0),
  fHTrkEffParGenPt_Dch(0),
  fHTrkEffDetRecPt_Dch(0),
  fHTrkEffParGenPt_Dn(0),
  fHTrkEffDetRecPt_Dn(0),
  fHTrkEffParGenPt_Ds(0),
  fHTrkEffDetRecPt_Ds(0),
  fHTrkEffParGenPt_cL(0),
  fHTrkEffDetRecPt_cL(0),
  fHTrkEffParGenPt_JPsi(0),
  fHTrkEffDetRecPt_JPsi(0),
  fHScaleFactor(0),
  fHScaleFactor100HC(0),
  fHEOverPVsPt(0x0),
  fHEMCalResponsePion(0x0), 
  fHEMCalResponseElec(0x0),
  fHEMCalResponseProton(0x0), 
  fHEMCalRecdPhidEta(0x0),
  fHEMCalRecdPhidEtaP(0x0),
  fHEMCalRecdPhidEtaM(0x0),
  fHEMCalRecdPhidEta_Truth(0x0),
  fHEMCalRecdPhidEtaP_Truth(0x0),
  fHEMCalRecdPhidEtaM_Truth(0x0),
  fHEMCalRecdPhidEtaposEta(0x0),
  fHEMCalRecdPhidEtanegEta(0x0),
  fHPhotonEdiff100HC(0x0),
  fHPhotonEdiff70HC(0),
  fHPhotonEdiff30HC(0),
  fHPhotonEdiff0HC(0x0),
  fHPhotonEVsClsE(0x0),
  fHistEsub1Pch(0x0),
  fHistEsub2Pch(0x0),
  fHistEsub1PchRat(0x0),
  fHistEsub2PchRat(0x0),
  fHClsEoverMcE_All(0x0),
  fHClsEoverMcE_Photon(0x0),
  fHClsEoverMcE_Elec(0x0),
  fHClsEoverMcE_Pion(0x0)
{
  // Constructor

}


//________________________________________________________________________
AliAnalysisTaskSOH::AliAnalysisTaskSOH(const char *name) :
  AliAnalysisTaskSE(name), 
  fESD(0), 
  fMC(0), 
  fEsdTrackCuts(0x0),
  fHybridTrackCuts1(0x0),
  fHybridTrackCuts2(0x0),
  fTrackIndices(0x0),
  fClusterIndices(0x0),
  fClusterArray(0x0),
  fOutputList(0x0),        
  fHEventStat(0),      
  fHTrkEffParGenPtEtaPhi(0), 
  fHTrkEffDetGenPtEtaPhi(0), 
  fHTrkEffDetRecPtEtaPhi(0), 
  fHTrkEffDetRecFakePtEtaPhi(0),
  fHTrkEffParGenPt_all(0),
  fHTrkEffDetRecPt_all(0),
  fHTrkEffParGenPt_Dch(0),
  fHTrkEffDetRecPt_Dch(0),
  fHTrkEffParGenPt_Dn(0),
  fHTrkEffDetRecPt_Dn(0),
  fHTrkEffParGenPt_Ds(0),
  fHTrkEffDetRecPt_Ds(0),
  fHTrkEffParGenPt_cL(0),
  fHTrkEffDetRecPt_cL(0),
  fHTrkEffParGenPt_JPsi(0),
  fHTrkEffDetRecPt_JPsi(0),
  fHScaleFactor(0),
  fHScaleFactor100HC(0),
  fHEOverPVsPt(0x0),
  fHEMCalResponsePion(0x0), 
  fHEMCalResponseElec(0x0),
  fHEMCalResponseProton(0x0), 
  fHEMCalRecdPhidEta(0x0),
  fHEMCalRecdPhidEtaP(0x0),
  fHEMCalRecdPhidEtaM(0x0),
  fHEMCalRecdPhidEta_Truth(0x0),
  fHEMCalRecdPhidEtaP_Truth(0x0),
  fHEMCalRecdPhidEtaM_Truth(0x0),
  fHEMCalRecdPhidEtaposEta(0x0),
  fHEMCalRecdPhidEtanegEta(0x0),
  fHPhotonEdiff100HC(0x0),
  fHPhotonEdiff70HC(0),
  fHPhotonEdiff30HC(0),
  fHPhotonEdiff0HC(0x0),
  fHPhotonEVsClsE(0x0),
  fHistEsub1Pch(0x0),
  fHistEsub2Pch(0x0),
  fHistEsub1PchRat(0x0),
  fHistEsub2PchRat(0x0),
  fHClsEoverMcE_All(0x0),
  fHClsEoverMcE_Photon(0x0),
  fHClsEoverMcE_Elec(0x0),
  fHClsEoverMcE_Pion(0x0)
{
  // Constructor

  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}


//________________________________________________________________________
AliAnalysisTaskSOH::~AliAnalysisTaskSOH()
{
  // Destructor.

  delete fEsdTrackCuts;
  delete fHybridTrackCuts1;
  delete fHybridTrackCuts2;
  delete fTrackIndices;
  delete fClusterIndices;
  delete fClusterArray;
}


//________________________________________________________________________
void AliAnalysisTaskSOH::UserCreateOutputObjects()
{
  // Create histograms, called once.

  OpenFile(1);
  
  fOutputList = new TList();
  fOutputList->SetOwner(1);

  fHEventStat = new TH1F("fHEventStat","Event statistics for analysis",8,0,8);
  fHEventStat->GetXaxis()->SetBinLabel(1,"Event");
  fHEventStat->GetXaxis()->SetBinLabel(2,"cluster");
  fHEventStat->GetXaxis()->SetBinLabel(3,"good cluster");
  fHEventStat->GetXaxis()->SetBinLabel(4,"cls/0-truth");
  fHEventStat->GetXaxis()->SetBinLabel(5,"cls/1-truth");
  fHEventStat->GetXaxis()->SetBinLabel(6,"cls/2-truth");
  fHEventStat->GetXaxis()->SetBinLabel(7,"cls/2-goodtruth");
  fHEventStat->GetXaxis()->SetBinLabel(8,"cls/>3-truth");
  fOutputList->Add(fHEventStat);

  fHTrkEffParGenPtEtaPhi = new TH3F("fHTrkEffParGenPtEtaPhi","Particle level truth Phi-Eta-p_{T} distribution of primary charged pions", 500, 0, 50, 60, -1.5, 1.5, 128, 0, 6.4);
  fHTrkEffParGenPtEtaPhi->GetXaxis()->SetTitle("p_{T}^{gen} (GeV/c)");
  fHTrkEffParGenPtEtaPhi->GetYaxis()->SetTitle("#eta");
  fHTrkEffParGenPtEtaPhi->GetZaxis()->SetTitle("#phi");
  fOutputList->Add(fHTrkEffParGenPtEtaPhi);
  
  fHTrkEffDetGenPtEtaPhi = new TH3F("fHTrkEffDetGenPtEtaPhi","Detector level truth Phi-Eta-p_{T} distribution of primary charged pions", 500, 0, 50, 60, -1.5, 1.5, 128, 0, 6.4);
  fHTrkEffDetGenPtEtaPhi->GetXaxis()->SetTitle("p_{T}^{gen} (GeV/c)");
  fHTrkEffDetGenPtEtaPhi->GetYaxis()->SetTitle("#eta");
  fHTrkEffDetGenPtEtaPhi->GetZaxis()->SetTitle("#phi");
  fOutputList->Add(fHTrkEffDetGenPtEtaPhi);

  fHTrkEffDetRecPtEtaPhi = new TH3F("fHTrkEffDetRecPtEtaPhi","Reconstructed track Phi-Eta-p_{T} distribution of primary charged pions", 500, 0, 50, 60, -1.5, 1.5, 128, 0, 6.4);
  fHTrkEffDetRecPtEtaPhi->GetXaxis()->SetTitle("p_{T}^{gen} (GeV/c)");
  fHTrkEffDetRecPtEtaPhi->GetYaxis()->SetTitle("#eta");
  fHTrkEffDetRecPtEtaPhi->GetZaxis()->SetTitle("#phi");
  fOutputList->Add(fHTrkEffDetRecPtEtaPhi);

  fHTrkEffDetRecFakePtEtaPhi = new TH3F("fHTrkEffDetRecFakePtEtaPhi","Reconstructed fake track Phi-Eta-p_{T} distribution of pions", 500, 0, 50, 60, -1.5, 1.5, 128, 0, 6.4);
  fHTrkEffDetRecFakePtEtaPhi->GetXaxis()->SetTitle("p_{T}^{gen} (GeV/c)");
  fHTrkEffDetRecFakePtEtaPhi->GetYaxis()->SetTitle("#eta");
  fHTrkEffDetRecFakePtEtaPhi->GetZaxis()->SetTitle("#phi");
  fOutputList->Add(fHTrkEffDetRecFakePtEtaPhi);

  fHTrkEffParGenPt_all = new TH1F("fHTrkEffParGenPt_all", "Truth p_{T} distribution;p_{T}^{rec} (GeV/c)",500,0.0,50);
  fOutputList->Add(fHTrkEffParGenPt_all);

  fHTrkEffDetRecPt_all = new TH1F("fHTrkEffDetRecPt_all", "Reconstructed track p_{T} distribution;p_{T}^{rec} (GeV/c)",500,0.0,50);
  fOutputList->Add(fHTrkEffDetRecPt_all);

  fHTrkEffParGenPt_Dch = new TH1F("fHTrkEffParGenPt_Dch", "Charged D meson truth p_{T} distribution;p_{T}^{rec} (GeV/c)",500,0.0,50);
  fOutputList->Add(fHTrkEffParGenPt_Dch);

  fHTrkEffDetRecPt_Dch = new TH1F("fHTrkEffDetRecPt_Dch", "Charged D meson reconstructed track p_{T} distribution;p_{T}^{rec} (GeV/c)",500,0.0,50);
  fOutputList->Add(fHTrkEffDetRecPt_Dch);

  fHTrkEffParGenPt_Dn = new TH1F("fHTrkEffParGenPt_Dn", "Neutral D meson truth p_{T} distribution;p_{T}^{rec} (GeV/c)",500,0.0,50);
  fOutputList->Add(fHTrkEffParGenPt_Dn);

  fHTrkEffDetRecPt_Dn = new TH1F("fHTrkEffDetRecPt_Dn", "Neutral D meson reconstructed track p_{T} distribution;p_{T}^{rec} (GeV/c)",500,0.0,50);
  fOutputList->Add(fHTrkEffDetRecPt_Dn);

  fHTrkEffParGenPt_Ds = new TH1F("fHTrkEffParGenPt_Ds", "D_{s} truth p_{T} distribution;p_{T}^{rec} (GeV/c)",500,0.0,50);
  fOutputList->Add(fHTrkEffParGenPt_Ds);

  fHTrkEffDetRecPt_Ds = new TH1F("fHTrkEffDetRecPt_Ds", "D_{s} reconstructed track p_{T} distribution;p_{T}^{rec} (GeV/c)",500,0.0,50);
  fOutputList->Add(fHTrkEffDetRecPt_Ds);

  fHTrkEffParGenPt_cL = new TH1F("fHTrkEffParGenPt_cL", "#Lambda_{c} truth p_{T} distribution;p_{T}^{rec} (GeV/c)",500,0.0,50);
  fOutputList->Add(fHTrkEffParGenPt_cL);

  fHTrkEffDetRecPt_cL = new TH1F("fHTrkEffDetRecPt_cL", "#Lambda_{c} reconstructed track p_{T} distribution;p_{T}^{rec} (GeV/c)",500,0.0,50);
  fOutputList->Add(fHTrkEffDetRecPt_cL);

  fHTrkEffParGenPt_JPsi = new TH1F("fHTrkEffParGenPt_JPsi", "J/#Psi truth p_{T} distribution;p_{T}^{rec} (GeV/c)",500,0.0,50);
  fOutputList->Add(fHTrkEffParGenPt_JPsi);

  fHTrkEffDetRecPt_JPsi = new TH1F("fHTrkEffDetRecPt_JPsi", "J/#Psi Reconstructed track p_{T} distribution;p_{T}^{rec} (GeV/c)",500,0.0,50);
  fOutputList->Add(fHTrkEffDetRecPt_JPsi);

  fHScaleFactor = new TH1F("fHScaleFactor", "Scale factor distribution without hadronic correction;Scale factor",100,0,10);
  fOutputList->Add(fHScaleFactor);

  fHScaleFactor100HC = new TH1F("fHScaleFactor100HC", "Scale factor distribution with 100% hadronic correction;Scale factor",100,0,10);
  fOutputList->Add(fHScaleFactor100HC);

  fHEOverPVsPt = new TH2F("fHEOverPVsPt", "E/P vs track p_{T}; p_{T} (GeV/c); E/P", 200 , 0, 4, 200, 0, 3.2);
  fOutputList->Add(fHEOverPVsPt);

  fHEMCalResponsePion = new TH2F("fHEMCalResponsePion", "Pion E/P vs track p_{T}; p_{T} (GeV/c); E/P", 100 , 0, 4, 100, 0, 3.2);
  fOutputList->Add(fHEMCalResponsePion);

  fHEMCalResponseElec = new TH2F("fHEMCalResponseElec", "Electron E/P vs track p_{T};  p_{T} (GeV/c); E/P", 100 , 0, 4, 100, 0, 3.2);
  fOutputList->Add(fHEMCalResponseElec);

  fHEMCalResponseProton = new TH2F("fHEMCalResponseProton", "Proton E/P vs track p_{T};  p_{T} (GeV/c); E/P", 100 , 0, 4, 100, 0, 3.2);
  fOutputList->Add(fHEMCalResponseProton);

  fHEMCalRecdPhidEta = new TH2F("fHEMCalRecdPhidEta","EMCAL Cluster-Track #Delta#phi-#Delta#eta; #Delta#eta; #Delta#phi",1000,-0.1,0.1,1000,-0.5,0.5);
  fOutputList->Add(fHEMCalRecdPhidEta);

  fHEMCalRecdPhidEtaP = new TH2F("fHEMCalRecdPhidEtaP","EMCAL Charge+ Cluster-Track #Delta#phi-#Delta#eta; #Delta#eta; #Delta#phi",1000,-0.1,0.1,1000,-0.5,0.5);
  fOutputList->Add(fHEMCalRecdPhidEtaP);

  fHEMCalRecdPhidEtaM = new TH2F("fHEMCalRecdPhidEtaM","EMCAL Charge- Cluster-Track #Delta#phi-#Delta#eta; #Delta#eta; #Delta#phi",1000,-0.1,0.1,1000,-0.5,0.5);
  fOutputList->Add(fHEMCalRecdPhidEtaM);

  fHEMCalRecdPhidEta_Truth = new TH2F("fHEMCalRecdPhidEta_Truth","EMCAL Cluster-Track(Truth matched) #Delta#phi-#Delta#eta; #Delta#eta; #Delta#phi",1000,-0.1,0.1,1000,-0.5,0.5);
  fOutputList->Add(fHEMCalRecdPhidEta_Truth);

  fHEMCalRecdPhidEtaP_Truth = new TH2F("fHEMCalRecdPhidEtaP_Truth","EMCAL Charge+ Cluster-Track(Truth matched) #Delta#phi-#Delta#eta; #Delta#eta; #Delta#phi",1000,-0.1,0.1,1000,-0.5,0.5);
  fOutputList->Add(fHEMCalRecdPhidEtaP_Truth);

  fHEMCalRecdPhidEtaM_Truth = new TH2F("fHEMCalRecdPhidEtaM_Truth","EMCAL Charge- Cluster-Track(Truth matched) #Delta#phi-#Delta#eta; #Delta#eta; #Delta#phi",1000,-0.1,0.1,1000,-0.5,0.5);
  fOutputList->Add(fHEMCalRecdPhidEtaM_Truth);

  fHEMCalRecdPhidEtaposEta = new TH2F("fHEMCalRecdPhidEtaposEta","(+eta track) EMCAL Cluster-Track #Delta#phi-#Delta#eta; #Delta#eta; #Delta#phi",1000,-0.1,0.1,1000,-0.5,0.5);
  fOutputList->Add(fHEMCalRecdPhidEtaposEta);

  fHEMCalRecdPhidEtanegEta = new TH2F("fHEMCalRecdPhidEtanegEta","(-eta track) EMCAL Cluster-Track #Delta#phi-#Delta#eta; #Delta#eta; #Delta#phi",1000,-0.1,0.1,1000,-0.5,0.5);
  fOutputList->Add(fHEMCalRecdPhidEtanegEta);

  fHPhotonEdiff100HC = new TH2F("fHPhotonEdiff100HC","Photon (E_{Truth}- E_{calc,100% HC})/E_{Truth} vs. E_{Truth}; E_{Truth} (GeV); (E_{Truth}- E_{calc,100% HC})/E_{Truth}",1000,0,10,600,-4.9,1.1);
  fOutputList->Add(fHPhotonEdiff100HC);

  fHPhotonEdiff70HC = new TH2F("fHPhotonEdiff70HC","Photon (E_{Truth}- E_{calc,70% HC})/E_{Truth} vs. E_{Truth}; E_{Truth} (GeV); (E_{Truth}- E_{calc,30% HC})/E_{Truth}",1000,0,10,600,-4.9,1.1);
  fOutputList->Add(fHPhotonEdiff70HC);

  fHPhotonEdiff30HC = new TH2F("fHPhotonEdiff30HC","Photon (E_{Truth}- E_{calc,30% HC})/E_{Truth} vs. E_{Truth}; E_{Truth} (GeV); (E_{Truth}- E_{calc,30% HC})/E_{Truth}",1000,0,10,600,-4.9,1.1);
  fOutputList->Add(fHPhotonEdiff30HC);

  fHPhotonEdiff0HC = new TH2F("fHPhotonEdiff0HC","Photon (E_{Truth}- E_{calc,0% HC})/E_{Truth} vs. E_{Truth}; E_{Truth} (GeV); (E_{Truth}- E_{cls})/E_{Truth}",1000,0,10,600,-4.9,1.1);
  fOutputList->Add(fHPhotonEdiff0HC);

  fHPhotonEVsClsE = new TH2F("fHPhotonEVsClsE","Cluster E vs. photon E_{Truth}; photon E_{Truth} (GeV); Cluster E (GeV)",500,0,5,500,0,5);
  fOutputList->Add(fHPhotonEVsClsE);

  fHistEsub1Pch =new  TH2F("fHistEsub1Pch", "(subtracted E in 100% HC) vs. total track P, clusters with 1 matching track; total track P (GeV/c); E_{sub}(GeV)" , 1000, 0., 10, 1000, 0., 10.);
  fOutputList->Add(fHistEsub1Pch);

  fHistEsub2Pch =new  TH2F("fHistEsub2Pch", "(subtracted E in 100% HC) vs. total track P, clusters with 2 matching tracks; total track P (GeV/c); E_{sub}(GeV)" , 1000, 0., 10, 1000, 0., 10.);
  fOutputList->Add(fHistEsub2Pch);

  fHistEsub1PchRat =new  TH2F("fHistEsub1PchRat", "(subtracted E in 100% HC)/total track P vs. total track P, clusters with 1 matching track; total track P (GeV/c); E_{sub}/P_{tot}" , 1000, 0., 10, 1100, 0., 1.1);
  fOutputList->Add(fHistEsub1PchRat);

  fHistEsub2PchRat =new  TH2F("fHistEsub2PchRat", "(subtracted E in 100% HC)/total track P vs. total track P, clusters with 2 matching tracks; total track P (GeV/c); E_{sub}/P_{tot}" , 1000, 0., 10, 1100, 0., 1.1);
  fOutputList->Add(fHistEsub2PchRat);

  Int_t bins[4] = {150, 150, 100, 200};
  Double_t xmin[4] = {1.3, -0.8, 0, 0};
  Double_t xmax[4] = {3.2, 0.8, 10, 2};

  fHClsEoverMcE_All = new THnSparseF("fHClsEoverMcE_All", "Cluster E/MC E, clusters with 1 matching particle; #phi; #eta; E (GeV); ClsE/McE", 4, bins, xmin, xmax);
  fOutputList->Add(fHClsEoverMcE_All);

  fHClsEoverMcE_Photon = new THnSparseF("fHClsEoverMcE_Photon", "Cluster E/MC E, clusters with 1 matching photon; #phi; #eta; E (GeV); ClsE/McE", 4, bins, xmin, xmax);
  fOutputList->Add(fHClsEoverMcE_Photon);

  fHClsEoverMcE_Elec = new THnSparseF("fHClsEoverMcE_Elec", "Cluster E/MC E, clusters with 1 matching electron; #phi; #eta; E (GeV); ClsE/McE", 4, bins, xmin, xmax);
  fOutputList->Add(fHClsEoverMcE_Elec);

  fHClsEoverMcE_Pion = new THnSparseF("fHClsEoverMcE_Pion", "Cluster E/MC E, clusters with 1 matching pion; #phi; #eta; E (GeV); ClsE/McE",4, bins, xmin, xmax);
  fOutputList->Add(fHClsEoverMcE_Pion);
  
  fTrackIndices = new TArrayI();
  fClusterIndices = new TArrayI();

  fClusterArray = new TObjArray();
  fClusterArray->SetOwner(1);

  PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskSOH::UserExec(Option_t *) 
{
  // Main loop, called for each event.

  // Post output data.
  fESD = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!fESD) {
    printf("ERROR: fESD not available\n");
    return;
  }

  fMC = MCEvent();
  if (!fMC) {
    printf("ERROR: fMC not available\n");
    return;
  }

  fHEventStat->Fill(0.5);

  if(fTrackIndices) 
    fTrackIndices->Reset();
  if(fClusterIndices) 
    fClusterIndices->Reset();
  if(fClusterArray)
    fClusterArray->Delete();

  ProcessTrack();
  ProcessCluster();
  ProcessMc();
  ProcessScaleFactor();
  
  PostData(1, fOutputList);
}   
   
//________________________________________________________________________
void  AliAnalysisTaskSOH::ProcessTrack()
{
  // Process track.

  fTrackIndices->Set(fESD->GetNumberOfTracks());
  AliDebug(3,Form("%s:%d Selecting tracks",(char*)__FILE__,__LINE__));

  Int_t isMth = 0;
  Int_t nTracks = 0;

  Float_t ClsPos[3] = {-999,-999,-999};
  Double_t emcTrkpos[3] = {-999,-999,-999};

  for(Int_t itr=0; itr<fESD->GetNumberOfTracks(); itr++)
  {
    AliESDtrack *esdtrack = fESD->GetTrack(itr);
    if(!esdtrack)continue;
    AliESDtrack *newTrack = GetAcceptTrack(esdtrack);
    if(!newTrack) continue;
  
    Double_t clsE = -1;
    Int_t clsIndex = newTrack->GetEMCALcluster();
    if(newTrack->GetEMCALcluster()>-1)
    {
      AliESDCaloCluster *cluster = fESD->GetCaloCluster(clsIndex);
      if(IsGoodCluster(cluster))
      {
        isMth=1;
	      
        cluster->GetPosition(ClsPos);
        TVector3 VClsPos(ClsPos[0], ClsPos[1], ClsPos[2]);
	      
        AliEMCALTrack EMCTrk(*newTrack);
        if(!EMCTrk.PropagateToGlobal(ClsPos[0], ClsPos[1], ClsPos[2], 0.0, 0.0)) {
          continue;
        }
        EMCTrk.GetXYZ(emcTrkpos);
        TVector3 VemcTrkPos(emcTrkpos[0],emcTrkpos[1],emcTrkpos[2]);
	      
        Double_t dPhi = VClsPos.Phi() - VemcTrkPos.Phi();
        if (dPhi < -1*TMath::Pi()) dPhi += (2*TMath::Pi());
        else if (dPhi > TMath::Pi()) dPhi -= (2*TMath::Pi());
	      
        Double_t dEta = VClsPos.Eta() - VemcTrkPos.Eta();
	      
        fHEMCalRecdPhidEta->Fill(dEta, dPhi);	

        if((newTrack->GetLabel())>-1 && (newTrack->GetLabel()) < fMC->GetNumberOfTracks())
        {
          AliVParticle *vParticle = fMC->GetTrack(newTrack->GetLabel());
          if(IsGoodMcParticle(vParticle, newTrack->GetLabel()))
          {
            fHEMCalRecdPhidEta_Truth->Fill(dEta, dPhi);
            if(vParticle->Charge() > 0) fHEMCalRecdPhidEtaP_Truth->Fill(dEta, dPhi);
            if(vParticle->Charge() < 0) fHEMCalRecdPhidEtaM_Truth->Fill(dEta, dPhi);
          }
        }

        if(esdtrack->Charge() > 0) {fHEMCalRecdPhidEtaP->Fill(dEta, dPhi);}
        if(esdtrack->Charge() < 0) {fHEMCalRecdPhidEtaM->Fill(dEta, dPhi);}

        if(VemcTrkPos.Eta() > 0) fHEMCalRecdPhidEtaposEta->Fill(dEta, dPhi);
        if(VemcTrkPos.Eta() < 0) fHEMCalRecdPhidEtanegEta->Fill(dEta, dPhi);
	      
        clsE = cluster->E();
        if(newTrack->P()>0) fHEOverPVsPt->Fill(newTrack->Pt(),clsE/newTrack->P());
      }
	  
      Int_t ipart = newTrack->GetLabel();
      if(ipart>-1 && ipart<fMC->GetNumberOfTracks())
      {
        AliVParticle* vParticle = fMC->GetTrack(ipart);
        if(isMth && vParticle)
        {
          if(TMath::Abs(vParticle->PdgCode())==211)
          {
            fHEMCalResponsePion->Fill(newTrack->Pt(),clsE/newTrack->P());
          }
          if(TMath::Abs(vParticle->PdgCode())==11)
          {
            fHEMCalResponseElec->Fill(newTrack->Pt(),clsE/newTrack->P());
          }
          if(TMath::Abs(vParticle->PdgCode())==2212)
          {
            fHEMCalResponseProton->Fill(newTrack->Pt(),clsE/newTrack->P());
          }
        }
      }
    }
    if(newTrack->GetLabel()<0 && newTrack->GetPID()==2) fHTrkEffDetRecFakePtEtaPhi->Fill(esdtrack->Pt(), esdtrack->Eta(), esdtrack->Phi());

 // Track Indices
    fTrackIndices->AddAt(itr,nTracks);
    nTracks++;
  }

  fTrackIndices->Set(nTracks);
}

//________________________________________________________________________
void AliAnalysisTaskSOH::ProcessCluster()
{
  // Process cluster.

  Int_t nCluster = 0;
  TLorentzVector gamma;
  Double_t vertex[3] = {0, 0, 0};
  fESD->GetVertex()->GetXYZ(vertex);
  const Int_t nCaloClusters = fESD->GetNumberOfCaloClusters(); 
  fClusterIndices->Set(nCaloClusters);
  Float_t ClsPos[3] = {-999,-999,-999};

  for(Int_t itr=0; itr<nCaloClusters; itr++) 
  {
    fHEventStat->Fill(1.5); 
    AliESDCaloCluster *cluster = fESD->GetCaloCluster(itr);
    if(!IsGoodCluster(cluster)) continue;
    cluster->GetMomentum(gamma, vertex);
    if (gamma.Pt() < 0.15) continue;
    fHEventStat->Fill(2.5);

    cluster->GetPosition(ClsPos);
    TVector3 VClsPos(ClsPos[0], ClsPos[1], ClsPos[2]);
  
    TArrayI *TrackLabels = cluster->GetTracksMatched();
    
    if(TrackLabels->GetSize() == 1)
    {
      AliESDtrack *esdtrack = fESD->GetTrack(TrackLabels->GetAt(0));
      AliESDtrack *newTrack = GetAcceptTrack(esdtrack);
      if(newTrack && TMath::Abs(newTrack->Eta())<0.7)
      {
	Double_t Esub = newTrack->P();
	if (Esub > cluster->E()) Esub = cluster->E();
	fHistEsub1Pch->Fill(newTrack->P(), Esub);
	fHistEsub1PchRat->Fill(newTrack->P(), Esub/newTrack->P());
      }
    }

    if(TrackLabels->GetSize() == 2)
    {
      AliESDtrack *esdtrack1 = fESD->GetTrack(TrackLabels->GetAt(0));
      AliESDtrack *esdtrack2 = fESD->GetTrack(TrackLabels->GetAt(1));
      AliESDtrack *newTrack1 = GetAcceptTrack(esdtrack1);
      AliESDtrack *newTrack2 = GetAcceptTrack(esdtrack2);
      if(newTrack1 && newTrack2 && TMath::Abs(newTrack1->Eta())<0.7  && TMath::Abs(newTrack2->Eta())<0.7)
      {
	Double_t Esub = newTrack1->P() + newTrack2->P();
	if (Esub > cluster->E()) Esub = cluster->E();
	fHistEsub2Pch->Fill(newTrack1->P() + newTrack2->P(), Esub);
	fHistEsub2PchRat->Fill(newTrack1->P() + newTrack2->P(), Esub/(newTrack1->P() + newTrack2->P()));
      }
      else if(newTrack1 && !(newTrack2) && TMath::Abs(newTrack1->Eta())<0.7)
      {
	Double_t Esub = newTrack1->P();
	if (Esub > cluster->E()) Esub = cluster->E();
	fHistEsub1Pch->Fill(newTrack1->P(), Esub);
	fHistEsub1PchRat->Fill(newTrack1->P(), Esub/newTrack1->P());
      }
      else if (!(newTrack1) && newTrack2 && TMath::Abs(newTrack2->Eta())<0.7)
      {
	Double_t Esub = newTrack2->P();
	if (Esub > cluster->E()) Esub = cluster->E();
	fHistEsub1Pch->Fill(newTrack2->P(), Esub);
	fHistEsub1PchRat->Fill(newTrack2->P(), Esub/newTrack2->P());
      }
      else {;}
    }

    TArrayI *MCLabels = cluster->GetLabelsArray();

    if(MCLabels->GetSize() == 0) fHEventStat->Fill(3.5);
    if(MCLabels->GetSize() == 1) 
    {
      fHEventStat->Fill(4.5);
      AliVParticle* vParticle1 = fMC->GetTrack(MCLabels->GetAt(0));
      if(IsGoodMcParticle(vParticle1, MCLabels->GetAt(0)))
      {
	Double_t Entries[4] = {VClsPos.Phi(), VClsPos.Eta(), vParticle1->E(), cluster->E()/vParticle1->E()}; 
	fHClsEoverMcE_All->Fill(Entries);
	if(vParticle1->PdgCode() == 22) 
	{
	  fHClsEoverMcE_Photon->Fill(Entries);
	}
	if(TMath::Abs(vParticle1->PdgCode()) == 11)
	{
	  fHClsEoverMcE_Elec->Fill(Entries);
	}
	if(TMath::Abs(vParticle1->PdgCode()) == 211) 
	{
	  fHClsEoverMcE_Pion->Fill(Entries);
	}
      }
    }
    if(MCLabels->GetSize() == 2) 
    {
      fHEventStat->Fill(5.5);
      AliVParticle* vParticle1 = fMC->GetTrack(MCLabels->GetAt(0));
      AliVParticle* vParticle2 = fMC->GetTrack(MCLabels->GetAt(1));
      if(IsGoodMcParticle(vParticle1, MCLabels->GetAt(0)) && IsGoodMcParticle(vParticle2, MCLabels->GetAt(1))) 
      {
	fHEventStat->Fill(6.5);
	if((vParticle1->PdgCode()==22) && (vParticle2->PdgCode()==22)) {;}
	else if((vParticle1->PdgCode()!=22) && (vParticle2->PdgCode()!=22)) {;}
	else 
	{
	  fClusterIndices->AddAt(itr,nCluster);
	  nCluster++;
	}
      }
    }
    if(MCLabels->GetSize() > 2) fHEventStat->Fill(7.5);

    AliESDCaloCluster *newCluster = new AliESDCaloCluster(*cluster);
 
    Double_t subE = 0;
    TArrayI arrayTrackMatched(fTrackIndices->GetSize());
    Int_t nGoodMatch = 0;

    for(Int_t j=0; j<fTrackIndices->GetSize(); j++)
    {
      AliESDtrack *trk = fESD->GetTrack(fTrackIndices->At(j));
      if(itr==trk->GetEMCALcluster())
      {
	arrayTrackMatched[nGoodMatch] = j;
	nGoodMatch ++;
	subE += trk->P();
      }
    }
  
    arrayTrackMatched.Set(nGoodMatch);
    newCluster->AddTracksMatched(arrayTrackMatched);
      
    Double_t clsE = newCluster->E();
    Double_t newE = clsE-subE; 
    if(newE<0) newE = 0;
    newCluster->SetDispersion(newE);
    fClusterArray->Add(newCluster);
  }

  fClusterIndices->Set(nCluster);
}
//________________________________________________________________________
void AliAnalysisTaskSOH::ProcessMc()
{
  // Process MC.

  for(Int_t i=0; i<fTrackIndices->GetSize(); i++)
  {
    AliESDtrack *esdtrack = fESD->GetTrack(fTrackIndices->At(i));
    Int_t index = esdtrack->GetLabel();
    if(index < 0 || index > fESD->GetNumberOfTracks()) continue;
    AliVParticle* vParticle = fMC->GetTrack(index);
    if(!IsGoodMcParticle(vParticle, index) && esdtrack->GetPID() == 2) fHTrkEffDetRecFakePtEtaPhi->Fill(esdtrack->Pt(), esdtrack->Eta(), esdtrack->Phi());
  }

  //tracking effciency
  for(Int_t ipart=0; ipart<fMC->GetNumberOfTracks(); ipart++)
  {
    AliVParticle* vParticle = fMC->GetTrack(ipart);
    Int_t pdgCode = vParticle->PdgCode();
    AliMCParticle *McParticle  = (AliMCParticle*)fMC->GetTrack(ipart);
    if(!IsGoodMcParticle(vParticle, ipart)) continue;  

    if(McParticle->GetMother() > 0 && TMath::Abs(vParticle->Eta())<0.9 && TMath::Abs(pdgCode)==211)
    {
      AliMCParticle *McMother = (AliMCParticle*)fMC->GetTrack(McParticle->GetMother());
      if(TMath::Abs(McMother->PdgCode()) == 411) fHTrkEffParGenPt_Dch->Fill(vParticle->Pt());
      else if(TMath::Abs(McMother->PdgCode()) == 421) fHTrkEffParGenPt_Dn->Fill(vParticle->Pt());
      else if(TMath::Abs(McMother->PdgCode()) == 431) fHTrkEffParGenPt_Ds->Fill(vParticle->Pt());	 
      else if(TMath::Abs(McMother->PdgCode()) == 4122) fHTrkEffParGenPt_cL->Fill(vParticle->Pt()); 
      else if(TMath::Abs(McMother->PdgCode()) == 443) fHTrkEffParGenPt_JPsi->Fill(vParticle->Pt()); 
      else {;}
    }
 
    if(TMath::Abs(vParticle->Eta())<0.9)
    {
      fHTrkEffParGenPt_all->Fill(vParticle->Pt());
      if(TMath::Abs(pdgCode==211)) fHTrkEffParGenPtEtaPhi->Fill(vParticle->Pt(), vParticle->Eta(), vParticle->Phi());
      
      for(Int_t j=0; j<fTrackIndices->GetSize(); j++)
      {
        AliESDtrack *esdtrack = fESD->GetTrack(fTrackIndices->At(j));
        if(esdtrack && esdtrack->GetLabel()==ipart)
        {
	  fHTrkEffDetRecPt_all->Fill(esdtrack->Pt());
	  if(TMath::Abs(pdgCode==211))
	  {
	    fHTrkEffDetGenPtEtaPhi->Fill(vParticle->Pt(), vParticle->Eta(), vParticle->Phi());
	    fHTrkEffDetRecPtEtaPhi->Fill(esdtrack->Pt(), esdtrack->Eta(), esdtrack->Phi());
	  }
	  if(McParticle->GetMother() > 0 && TMath::Abs(pdgCode) == 211)
	  {
	    AliMCParticle *McMother = (AliMCParticle*)fMC->GetTrack(McParticle->GetMother());
	    if(TMath::Abs(McMother->PdgCode()) == 411) fHTrkEffDetRecPt_Dch->Fill(esdtrack->Pt());
	    else if(TMath::Abs(McMother->PdgCode()) == 421) fHTrkEffDetRecPt_Dn->Fill(esdtrack->Pt());
	    else if(TMath::Abs(McMother->PdgCode()) == 431) fHTrkEffDetRecPt_Ds->Fill(esdtrack->Pt());	 
	    else if(TMath::Abs(McMother->PdgCode()) == 4122) fHTrkEffDetRecPt_cL->Fill(esdtrack->Pt()); 
	    else if(TMath::Abs(McMother->PdgCode()) == 443) fHTrkEffDetRecPt_JPsi->Fill(esdtrack->Pt()); 
	    else {;}
	  }
    //cluster E vs. truth photon energy
	  for(Int_t k=0; k<fClusterIndices->GetSize(); k++)
	  {
	    AliESDCaloCluster *cluster = fESD->GetCaloCluster(fClusterIndices->At(k));
	    Double_t clsE = cluster->E();
	    TArrayI *MCLabels = cluster->GetLabelsArray();
	    AliVParticle* vParticle1 = fMC->GetTrack(MCLabels->GetAt(0));
	    AliVParticle* vParticle2 = fMC->GetTrack(MCLabels->GetAt(1));
	    
	    if(vParticle1->PdgCode()==22 && vParticle2 == vParticle)
	    {
	      fHPhotonEdiff0HC->Fill(vParticle1->E(), (vParticle1->E() - clsE)/vParticle1->E());
	      fHPhotonEVsClsE->Fill(vParticle1->E(), clsE);

	      if((clsE - 0.3*esdtrack->E())<0) fHPhotonEdiff30HC->Fill(vParticle1->E(), 1);
	      else  fHPhotonEdiff30HC->Fill(vParticle1->E(), (vParticle1->E() + 0.3*esdtrack->E() - clsE)/vParticle1->E());

	      if((clsE - 0.7*esdtrack->E())<0) fHPhotonEdiff70HC->Fill(vParticle1->E(), 1);
	      else  fHPhotonEdiff70HC->Fill(vParticle1->E(), (vParticle1->E() + 0.7*esdtrack->E() - clsE)/vParticle1->E());

	      if((clsE - esdtrack->E())<0) fHPhotonEdiff100HC->Fill(vParticle1->E(), 1);
	      else  fHPhotonEdiff100HC->Fill(vParticle1->E(), (vParticle1->E() + esdtrack->E() - clsE)/vParticle1->E());
	      continue;
	    }
	    if(vParticle2->PdgCode()==22 && vParticle1 == vParticle)
	    {
	      fHPhotonEdiff0HC->Fill(vParticle2->E(), (vParticle2->E() - clsE)/vParticle2->E());
	      fHPhotonEVsClsE->Fill(vParticle2->E(), clsE);

	      if((clsE - 0.3*esdtrack->E())<0) fHPhotonEdiff30HC->Fill(vParticle2->E(), 1);
	      else fHPhotonEdiff30HC->Fill(vParticle2->E(), (vParticle2->E() + 0.3*esdtrack->E() - clsE)/vParticle2->E());

	      if((clsE - 0.7*esdtrack->E())<0) fHPhotonEdiff70HC->Fill(vParticle2->E(), 1);
	      else fHPhotonEdiff70HC->Fill(vParticle2->E(), (vParticle2->E() + 0.7*esdtrack->E() - clsE)/vParticle2->E());

	      if((clsE-esdtrack->E())<0) fHPhotonEdiff100HC->Fill(vParticle2->E(), 1);
	      else fHPhotonEdiff100HC->Fill(vParticle2->E(), (vParticle2->E() + esdtrack->E() - clsE)/vParticle2->E());
	    }
	  }
          break;
        }
      }
    }
  } 
}

//________________________________________________________________________
void AliAnalysisTaskSOH::ProcessScaleFactor()
{
  // Scale factor. 

  const Double_t phiMax = 180 * TMath::DegToRad();
  const Double_t phiMin = 80 * TMath::DegToRad();
  const Double_t TPCArea= 2*TMath::Pi()*1.8;
  const Double_t EMCArea = (phiMax-phiMin)*1.4;

  Double_t PtEMC = 0;
  Double_t PtTPC = 0;

  for(Int_t j=0; j<fTrackIndices->GetSize(); j++)
  {
    AliESDtrack *trk = fESD->GetTrack(fTrackIndices->At(j));
    Double_t eta = trk->Eta();
    Double_t phi = trk->Phi();
    if(TMath::Abs(eta)<0.9) PtTPC += trk->Pt();
    if(TMath::Abs(eta)<0.7  && phi > phiMin && phi < phiMax ) PtEMC += trk->Pt();
  }

  Double_t EtWithHadCorr = 0;
  Double_t EtWithoutHadCorr = 0;
  Double_t vertex[3] = {0, 0, 0};
  fESD->GetVertex()->GetXYZ(vertex);
  TLorentzVector gamma;

  for(Int_t i=0; i<fClusterArray->GetEntriesFast(); i++)
  {
    AliESDCaloCluster *cluster = (AliESDCaloCluster*)fClusterArray->At(i);
    cluster->GetMomentum(gamma, vertex);
    Double_t sinTheta = TMath::Sqrt(1-TMath::Power(gamma.CosTheta(),2));
    EtWithoutHadCorr +=  cluster->E() * sinTheta;
    EtWithHadCorr += cluster->GetDispersion() * sinTheta;
  }

  if(PtTPC>0)
  {
    fHScaleFactor->Fill((PtEMC+EtWithoutHadCorr)/EMCArea * TPCArea/PtTPC);
    fHScaleFactor100HC->Fill((PtEMC+EtWithHadCorr)/EMCArea * TPCArea/PtTPC);
  }
}

//________________________________________________________________________
AliESDtrack *AliAnalysisTaskSOH::GetAcceptTrack(AliESDtrack *esdtrack)
{
  // Get accepted track.

  static AliESDtrack newTrack;
  if(fEsdTrackCuts->AcceptTrack(esdtrack))
  {
    esdtrack->Copy(newTrack);
    newTrack.SetTRDQuality(0);
  }
  else if(fHybridTrackCuts1->AcceptTrack(esdtrack))
  {
    if(esdtrack->GetConstrainedParam())
    {
      esdtrack->Copy(newTrack);
      const AliExternalTrackParam* constrainParam = esdtrack->GetConstrainedParam();
      newTrack.Set(constrainParam->GetX(),constrainParam->GetAlpha(),constrainParam->GetParameter(),constrainParam->GetCovariance());
      newTrack.SetTRDQuality(1);		
    }
    else 
      return 0x0;
  }
  else if(fHybridTrackCuts2->AcceptTrack(esdtrack))
  {
    if(esdtrack->GetConstrainedParam())
    {
      esdtrack->Copy(newTrack);
      const AliExternalTrackParam* constrainParam = esdtrack->GetConstrainedParam();
      newTrack.Set(constrainParam->GetX(),constrainParam->GetAlpha(),constrainParam->GetParameter(),constrainParam->GetCovariance());
      newTrack.SetTRDQuality(2);		
    }
    else 
      return 0x0;
  }
  else
  {
    return 0x0;
  }

  return &newTrack;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSOH::IsGoodMcParticle(AliVParticle* vParticle, Int_t ipart)
{
  // Return true if good MC particle.

  if(!vParticle) return kFALSE;
  if(!fMC->IsPhysicalPrimary(ipart)) return kFALSE;
  if (TMath::Abs(vParticle->Eta())>2) return kFALSE;
  if(vParticle->Pt()<0.15) return kFALSE;
  return kTRUE;
}

//________________________________________________________________________
Bool_t AliAnalysisTaskSOH::IsGoodCluster(AliESDCaloCluster *cluster)
{
  // Return true if good cluster.

  if(!cluster) return kFALSE;
  if (!cluster->IsEMCAL()) return kFALSE;
  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskSOH::Terminate(Option_t *) 
{
  // Terminate analysis.
}
