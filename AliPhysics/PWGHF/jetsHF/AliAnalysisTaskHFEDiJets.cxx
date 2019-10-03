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


///////////////////////////////////////////
//          Task for HFE dijet           //
//  Author: Deepa Thomas, Shingo Sakai   //
///////////////////////////////////////////

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "THnSparse.h"
#include "TLorentzVector.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"

#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliESDtrackCuts.h"
#include "AliAODEvent.h"
#include "AliAODHandler.h"

#include "AliPID.h"
#include "AliESDpid.h"
#include "AliAODPid.h"
#include "AliPIDResponse.h"
#include "AliCFContainer.h"
#include "AliCFManager.h"
#include "AliKFParticle.h"
#include "AliKFVertex.h"

// #include <fastjet/config.h>
// #include <fastjet/PseudoJet.hh>
// #include <fastjet/JetDefinition.hh>
// #include <fastjet/ClusterSequence.hh>
// #include <fastjet/ClusterSequenceArea.hh>
// #include <fastjet/AreaDefinition.hh>

#include "AliAnalysisTaskHFEDiJets.h"
#include "AliAnalysisHelperJetTasks.h" 

ClassImp(AliAnalysisTaskHFEDiJets)

  //________________________________________________________________________
  AliAnalysisTaskHFEDiJets::AliAnalysisTaskHFEDiJets(const char *name)
: AliAnalysisTaskSE(name),
  fVevent(0),
  fESD(0),
  fAOD(0),
  fpidResponse(0),
  fJetRcut(0.3),
  fInputParticlesJet(),
  fBkgMedian(0),
  fBkgSigma(0),
  fBkgMeanArea(0),
  fInvmassCut(0.10),
  fMinJetPt(10),
  fJetDef(0),
  fGhostArea(0),
  fAreaDef(0),
  fCSA(0),
  fJets(0),
  fBkgJetDef(0),
  fBkgGhostArea(0),
  fBkgAreaDef(0),
  fBkgCSA(0),
  fRange(0),
  fBkgJets(0),
  fOutputList(0),
  fNevents(0),
  fHistClustE(0),
  fEMCClsEtaPhi(0),
  fTrkPt(0),
  fTrketa(0),
  fTrkphi(0),
  fTPCnsig(0),
  fHistPtMatch(0),
  fEMCTrkMatch(0),
  fEMCTrkPt(0),
  fEMCTrketa(0),
  fEMCTrkphi(0),
  fEMCTPCnsig(0),
  fHistNsigEop(0),
  fHistEop(0),
  fM20(0),
  fM02(0),
  fM20EovP(0),
  fM02EovP(0),
  fInvmassLS(0),
  fInvmassULS(0),
  fJetInputTrPhi(0),
  fJetInputTrEta(0),
  fBgMedian_Brf(0),
  fBgArea_Bfr(0),
  fBgPhi_Bfr(0),
  fBgEta_Bfr(0),
  fBgPt_Bfr(0),
  fBgP_Bfr(0),
  fBgMedian_Afr(0),
  fSigArea(0),
  fSigPt_BfrSub(0),
  fSigPt_Sub(0),
  fSigPhi(0),
  fSigEta(0),
  fSigP(0),
  fSigR(0),
  fSigEtaPhi(0),
  fSigEtainEMC(0),
  fBgArea_Afr(0),
  fBgPhi_Afr(0),
  fBgEta_Afr(0),
  fBgPt_Afr(0),
  fBgP_Afr(0),
  fPtAllElec(0),
  fNthJetwEle(0),
  fPtJetwElec(0),
  fPtElecinJet(0),
  fPtElecLeadInJet(0),
  fEtaPhiJetwElec(0),
  fEtaJetwElecinEMC(0),
  fNthJetwPEle(0),
  fPtJetwPElec(0),
  fPtPElecinJet(0),
  fPtPElecLeadInJet(0),
  fEtaPhiJetwPElec(0),
  fEtaJetwPElecinEMC(0),
  fPtLeadJetwElec(0),
  fPtElecinLeadJet(0),
  fPtElecLeadInLeadJet(0),
  fPtTrigJetwInclE(0),
  fPtTrigJetwSemiInclE(0),
  fPtTrigJetwPhotoE(0),
  fSprsInclEleDphiDeta(0),
  fSprsPhotEleDphiDeta(0),
  fSprsSemiIncEleDphiDeta(0)
{
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  DefineInput(0, TChain::Class());
  // Output slot #0 id reserved by the base class for AOD
  // Output slot #1 writes into a TH1 container
  DefineOutput(1, TList::Class());
}
//________________________________________________________________________
AliAnalysisTaskHFEDiJets::AliAnalysisTaskHFEDiJets()
  : AliAnalysisTaskSE("DefaultTask_HfeEMCQA"),
  fVevent(0),
  fESD(0),
  fAOD(0),
  fpidResponse(0),
  fJetRcut(0.3),
  fInputParticlesJet(),
  fBkgMedian(0),
  fBkgSigma(0),
  fBkgMeanArea(0),
  fInvmassCut(0.10),
  fMinJetPt(10),
  fJetDef(0),
  fGhostArea(0),
  fAreaDef(0),
  fCSA(0),
  fJets(0),
  fBkgJetDef(0),
  fBkgGhostArea(0),
  fBkgAreaDef(0),
  fBkgCSA(0),
  fRange(0),
  fBkgJets(0),
  fOutputList(0),
  fNevents(0),
  fHistClustE(0),
  fEMCClsEtaPhi(0),
  fTrkPt(0),
  fTrketa(0),
  fTrkphi(0),
  fTPCnsig(0),
  fHistPtMatch(0),
  fEMCTrkMatch(0),
  fEMCTrkPt(0),
  fEMCTrketa(0),
  fEMCTrkphi(0),
  fEMCTPCnsig(0),
  fHistNsigEop(0),
  fHistEop(0),
  fM20(0),
  fM02(0),
  fM20EovP(0),
  fM02EovP(0),
  fInvmassLS(0),
  fInvmassULS(0),
  fJetInputTrPhi(0),
  fJetInputTrEta(0),
  fBgMedian_Brf(0),
  fBgArea_Bfr(0),
  fBgPhi_Bfr(0),
  fBgEta_Bfr(0),
  fBgPt_Bfr(0),
  fBgP_Bfr(0),
  fBgMedian_Afr(0),
  fSigArea(0),
  fSigPt_BfrSub(0),
  fSigPt_Sub(0),
  fSigPhi(0),
  fSigEta(0),
  fSigP(0),
  fSigR(0),
  fSigEtaPhi(0),
  fSigEtainEMC(0),
  fBgArea_Afr(0),
  fBgPhi_Afr(0),
  fBgEta_Afr(0),
  fBgPt_Afr(0),
  fBgP_Afr(0),
  fPtAllElec(0),
  fNthJetwEle(0),
  fPtJetwElec(0),
  fPtElecinJet(0),
  fPtElecLeadInJet(0),
  fEtaPhiJetwElec(0),
  fEtaJetwElecinEMC(0),
  fNthJetwPEle(0),
  fPtJetwPElec(0),
  fPtPElecinJet(0),
  fPtPElecLeadInJet(0),
  fEtaPhiJetwPElec(0),
  fEtaJetwPElecinEMC(0),
  fPtLeadJetwElec(0),
  fPtElecinLeadJet(0),
  fPtElecLeadInLeadJet(0),
  fPtTrigJetwInclE(0),
  fPtTrigJetwSemiInclE(0),
  fPtTrigJetwPhotoE(0),
  fSprsInclEleDphiDeta(0),
  fSprsPhotEleDphiDeta(0),
  fSprsSemiIncEleDphiDeta(0)
{
  //Default constructor
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
AliAnalysisTaskHFEDiJets::~AliAnalysisTaskHFEDiJets()
{
  //Destructor
  delete fOutputList;
  delete fSprsInclEleDphiDeta;
  delete fSprsPhotEleDphiDeta;
  delete fSprsSemiIncEleDphiDeta;
}
//________________________________________________________________________
void AliAnalysisTaskHFEDiJets::UserCreateOutputObjects()
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

  ////////////////
  //Output list//
  ///////////////
  fOutputList = new TList();
  fOutputList->SetOwner(kTRUE);

  fNevents = new TH1F("fNevents","No of events",3,-0.5,2.5);
  fOutputList->Add(fNevents);
  fNevents->GetYaxis()->SetTitle("counts");
  fNevents->GetXaxis()->SetBinLabel(1,"All");
  fNevents->GetXaxis()->SetBinLabel(2,"With >2 Trks");
  fNevents->GetXaxis()->SetBinLabel(3,"Vtx_{z}<10cm");

  fHistClustE = new TH1F("fHistClustE", "EMCAL cluster energy distribution; Cluster E;counts", 500, 0.0, 50.0);
  fOutputList->Add(fHistClustE);

  fEMCClsEtaPhi = new TH2F("fEMCClsEtaPhi","EMCAL cluster #eta and #phi distribution;#eta;#phi",100,-1,1,200,0,6.3);
  fOutputList->Add(fEMCClsEtaPhi);

  fTrkPt = new TH1F("fTrkPt","p_{T} distribution of all tracks;p_{T} (GeV/c);counts",1000,0,100);
  fOutputList->Add(fTrkPt);

  fTrketa = new TH1F("fTrketa","All Track #eta distribution;#eta;counts",100,-1.5,1.5);
  fOutputList->Add(fTrketa);

  fTrkphi = new TH1F("fTrkphi","All Track #phi distribution;#phi;counts",100,0,6.3);
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

  fEMCTrkphi = new TH1F("fEMCTrkphi","#phi distribution of tracks matched to EMCAL;#phi;counts",100,0,6.3);
  fOutputList->Add(fEMCTrkphi);

  fEMCTPCnsig = new TH2F("fEMCTPCnsig","TPC Nsigma distribution of tracks matched to EMCAL;p (GeV/c);#sigma_{TPC-dE/dx}",1000,0,50,200,-10,10);
  fOutputList->Add(fEMCTPCnsig);

  fHistEop = new TH2F("fHistEop", "E/p distribution;p_{T} (GeV/c);E/p", 200,0,20,60, 0.0, 3.0);
  fOutputList->Add(fHistEop);

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

  fJetInputTrPhi = new TH1F("fJetInputTrPhi","#phi distribution of jet input tracks;#phi;counts",100,0,6.3);
  fOutputList->Add(fJetInputTrPhi);

  fJetInputTrEta = new TH1F("fJetInputTrEta","#eta distribution of jet input tracks;#eta;counts",100,-1,1);
  fOutputList->Add(fJetInputTrEta);

  fBgMedian_Brf = new TH1F("fBgMedian_Brf","Median of the background before removing leading jets",1000,0,10);
  fOutputList->Add(fBgMedian_Brf);

  fBgArea_Bfr = new TH1F("fBgArea_Bfr","Area of the background jets before removing leading jets",1000,0,5);
  fOutputList->Add(fBgArea_Bfr);

  fBgPhi_Bfr = new TH1F("fBgPhi_Bfr","Phi of background jets before removing leading jets",100,0,6.3);
  fOutputList->Add(fBgPhi_Bfr);

  fBgEta_Bfr = new TH1F("fBgEta_Bfr","Eta of background jets before removing leading jets",100,-1,1);
  fOutputList->Add(fBgEta_Bfr);

  fBgPt_Bfr = new TH1F("fBgPt_Bfr","Pt of background jets before removing leading jets",1500,0,150);
  fOutputList->Add(fBgPt_Bfr);

  fBgP_Bfr = new TH1F("fBgP_Bfr","P of background jets before removing leading jets",1500,0,150);
  fOutputList->Add(fBgP_Bfr);

  fBgMedian_Afr = new TH1F("fBgMedian_Afr","Median of the background after removing leading jets",1000,0,10);
  fOutputList->Add(fBgMedian_Afr);

  fBgArea_Afr = new TH1F("fBgArea_Afr","Area of the background jets after removing leading jets",1000,0,5);
  fOutputList->Add(fBgArea_Afr);

  fBgPhi_Afr = new TH1F("fBgPhi_Afr","Phi of background jets after removing leading jets",100,0,6.3);
  fOutputList->Add(fBgPhi_Afr);

  fBgEta_Afr = new TH1F("fBgEta_Afr","Eta of background jets after removing leading jets",100,-1,1);
  fOutputList->Add(fBgEta_Afr);

  fBgPt_Afr = new TH1F("fBgPt_Afr","Pt of background jets after removing leading jets",1500,0,150);
  fOutputList->Add(fBgPt_Afr);

  fBgP_Afr = new TH1F("fBgP_Afr","P of background jets after removing leading jets",1500,0,150);
  fOutputList->Add(fBgP_Afr);

  fSigArea = new TH1F("fSigArea","Area of the Signal jets",1000,0,5);
  fOutputList->Add(fSigArea);

  fSigPt_BfrSub = new TH1F("fSigPt_BfrSub","Pt of signal jet before removing background",1500,0,150);
  fOutputList->Add(fSigPt_BfrSub);

  fSigPt_Sub = new TH1F("fSigPt_Sub","Pt of signal jet after removing background",1500,0,150);
  fOutputList->Add(fSigPt_Sub);

  fSigPhi = new TH1F("fSigPhi","Phi of signal jets",100,0,6.3);
  fOutputList->Add(fSigPhi);

  fSigEta= new TH1F("fSigEta","Eta of signal jets",100,-1,1);
  fOutputList->Add(fSigEta);

  fSigP = new TH1F("fSigP","P of signal jets",1500,0,150);
  fOutputList->Add(fSigP);

  fSigR = new TH1F("fSigR","R of signal jets",1000,0,100);
  fOutputList->Add(fSigR);

  fSigEtaPhi = new TH2F("fSigEtaPhi","#eta vs #phi of signal jets;#eta;#phi",100,-1,  1,200,0,6.3);
  fOutputList->Add(fSigEtaPhi);

  fSigEtainEMC = new TH1F("fSigEtainEMC","#eta of signal jets if in EMCal #phi acceptance",100,-1,1);
  fOutputList->Add(fSigEtainEMC);

  fPtAllElec = new TH1F("fPtAllElec","Pt of all electrons",1000,0,100);
  fOutputList->Add(fPtAllElec);

  fInvmassLS = new TH1F("fInvmassLS", "Inv mass of LS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 1000,0,1.0);
  fOutputList->Add(fInvmassLS);

  fInvmassULS = new TH1F("fInvmassULS", "Inv mass of ULS (e,e) for pt^{e}>1; mass(GeV/c^2); counts;", 1000,0,1.0);
  fOutputList->Add(fInvmassULS);

  fNthJetwEle = new TH1F("fNthJetwEle","Nth jet containing elec; N; counts",25,0,25);
  fOutputList->Add(fNthJetwEle);

  fPtJetwElec = new TH1F("fPtJetwElec","Pt of jets with electron in it",1500,0,150);
  fOutputList->Add(fPtJetwElec);

  fPtElecinJet = new TH1F("fPtElecinJet","Pt of electrons in jet",1000,0,100);
  fOutputList->Add(fPtElecinJet);

  fPtElecLeadInJet = new TH1F("fPtElecLeadInJet","Pt of electron if it is the leading track in jet",1000,0,100);
  fOutputList->Add(fPtElecLeadInJet);

  fEtaPhiJetwElec = new TH2F("fEtaPhiJetwElec","#eta vs #phi of jet containing elec;#eta;#phi",100,-1,1,200,0,6.3);
  fOutputList->Add(fEtaPhiJetwElec);

  fEtaJetwElecinEMC = new TH1F("fEtaJetwElecinEMC","#eta of jets containing elec if in EMCal #phi acceptance",100,-1,1);
  fOutputList->Add(fEtaJetwElecinEMC);

  fNthJetwPEle = new TH1F("fNthJetwPEle","Nth jet containing photonic elec; N; counts",25,0,25);
  fOutputList->Add(fNthJetwPEle);

  fPtJetwPElec = new TH1F("fPtJetwPElec","Pt of jets with photonic electron in it",1500,0,150);
  fOutputList->Add(fPtJetwPElec);

  fPtPElecinJet = new TH1F("fPtPElecinJet","Pt of photonic electrons in jet",1000,0,100);
  fOutputList->Add(fPtPElecinJet);

  fPtPElecLeadInJet = new TH1F("fPtPElecLeadInJet","Pt of photonic electron if it is the leading track in jet",1000,0,100);
  fOutputList->Add(fPtPElecLeadInJet);

  fEtaPhiJetwPElec = new TH2F("fEtaPhiJetwPElec","#eta vs #phi of jet containing photonic elec;#eta;#phi",100,-1,1,200,0,6.3);
  fOutputList->Add(fEtaPhiJetwPElec);

  fEtaJetwPElecinEMC = new TH1F("fEtaJetwPElecinEMC","#eta of jets containing photonic elec if in EMCal #phi acceptance",100,-1,1);
  fOutputList->Add(fEtaJetwPElecinEMC);

  fPtLeadJetwElec = new TH1F("fPtLeadJetwElec","Pt of leading jets with electron in it",1500,0,150);
  fOutputList->Add(fPtLeadJetwElec);

  fPtElecinLeadJet = new TH1F("fPtElecinLeadJet","Pt of electrons in leading jets",1000,0,100);
  fOutputList->Add(fPtElecinLeadJet);

  fPtElecLeadInLeadJet = new TH1F("fPtElecLeadInLeadJet","Pt of electron if it is the leading track in leading jets",1000,0,100);
  fOutputList->Add(fPtElecLeadInLeadJet);

  fPtTrigJetwInclE = new TH1F("fPtTrigJetwInclE","Pt of trgger jets with inclusive electron in it",1500,0,150);
  fOutputList->Add(fPtTrigJetwInclE);

  fPtTrigJetwSemiInclE = new TH1F("fPtTrigJetwSemiInclE","Pt of trgger jets with Semi-inclusive electron in it",1500,0,150);
  fOutputList->Add(fPtTrigJetwSemiInclE);

  fPtTrigJetwPhotoE = new TH1F("fPtTrigJetwPhotoE","Pt of trgger jets with Photonic electron in it",1500,0,150);
  fOutputList->Add(fPtTrigJetwPhotoE);

  //------THnsparse------
  Int_t bin[7] = {5,5,100,17,17,64,100}; //TrigJetid, AssoJetid, ptElec, ptTrigJet, ptAssoJet, Dphi, Deta
  Double_t xmin[7] = {0,0,0,15,15,-TMath::Pi()/2,-1.8};
  Double_t xmax[7] = {5,5,50,100,100,(3*TMath::Pi())/2,1.8};

  fSprsInclEleDphiDeta = new THnSparseD("fSprsInclEleDphiDeta","Sparse for Dphi and Deta with Inclusive electron",7,bin,xmin,xmax);
  fOutputList->Add(fSprsInclEleDphiDeta);

  fSprsPhotEleDphiDeta = new THnSparseD("fSprsPhotEleDphiDeta","Sparse for Dphi and Deta with Photonic electron",7,bin,xmin,xmax);
  fOutputList->Add(fSprsPhotEleDphiDeta);

  fSprsSemiIncEleDphiDeta = new THnSparseD("fSprsSemiIncEleDphiDeta","Sparse for Dphi and Deta with Semi-inclusive electron",7,bin,xmin,xmax);
  fOutputList->Add(fSprsSemiIncEleDphiDeta);

  PostData(1,fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskHFEDiJets::UserExec(Option_t *)
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

  fAOD = dynamic_cast<AliAODEvent*>(InputEvent());
  if (fAOD) {
    // printf("fAOD available\n");
    //return;
  }

  ///////////////////
  //PID initialised//
  //////////////////
  fpidResponse = fInputHandler->GetPIDResponse();

  ////////////////
  //Event vertex//
  ///////////////
  Int_t ntracks;
  ntracks = fVevent->GetNumberOfTracks();
  //printf("There are %d tracks in this event\n",ntracks);

  fNevents->Fill(0); //all events
  Double_t Zvertex = -100, Xvertex = -100, Yvertex = -100;
  const AliVVertex *pVtx = fVevent->GetPrimaryVertex();
  Double_t NcontV = pVtx->GetNContributors();
  if(NcontV<2)return;
  fNevents->Fill(1); //events with 2 tracks

  ////////////////////
  //event selection//
  ///////////////////
  if(fabs(Zvertex>10.0))return;
  fNevents->Fill(2); //events after z vtx cut

  //////////////////////
  //EMcal cluster info//
  //////////////////////
  EMCalClusterInfo();

  //////////////
  //Jet Recons//
  //////////////
  //------++++Input tracks to fastjet
  Double_t fMaxpT = 0.15;
  InputTracksForJetReco(fMaxpT);
  //cout << "max pt : " << fMaxpT <<endl;

  //------++++Jet Reconstruction
  DoJetReconstruction();
  DoBkgJetReconstruction();

  //-----Plot reco jet properties --> only after reconstructing background
  PlotRecoJetProperties(fSigArea, fSigPhi, fSigEta, fSigPt_BfrSub,  fSigPt_Sub, fSigP, fSigR, fSigEtaPhi, fSigEtainEMC);

  /////////////////////////////////
  //Look for kink mother for AOD//
  /////////////////////////////////
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

  ////////////////
  //Electron ID///
  ////////////////
  for (Int_t iTracks = 0; iTracks < ntracks; iTracks++) {

    AliVParticle* Vtrack = 0x0;
    Vtrack  = fVevent->GetTrack(iTracks);

    if (!Vtrack) {
      printf("ERROR: Could not receive track %d\n", iTracks);
      continue;
    }
    AliVTrack *track = dynamic_cast<AliVTrack*>(Vtrack);
    AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(Vtrack);

    //+++++----Apply track cuts
    if(fAOD)
      if(!atrack->TestFilterMask(AliAODTrack::kTrkGlobalNoDCA)) continue; //mimimum cuts

    //reject kink
    Bool_t kinkmotherpass = kTRUE;
    for(Int_t kinkmother = 0; kinkmother < numberofmotherkink; kinkmother++) {
      if(track->GetID() == listofmotherkink[kinkmother]) {
        kinkmotherpass = kFALSE;
        continue;
      }
    }
    if(!kinkmotherpass) continue;

    //other cuts
    Double_t d0z0[2]={-999,-999}, cov[3];
    Double_t DCAxyCut = 2.4, DCAzCut = 3.2;
    if(fAOD){
      if(atrack->GetTPCNcls() < 80) continue;
      if(atrack->GetITSNcls() < 3) continue;
      if((!(atrack->GetStatus()&AliESDtrack::kITSrefit)|| (!(atrack->GetStatus()&AliESDtrack::kTPCrefit)))) continue;
      if(!(atrack->HasPointOnITSLayer(0) || atrack->HasPointOnITSLayer(1))) continue;

      if(atrack->PropagateToDCA(pVtx, fVevent->GetMagneticField(), 20., d0z0, cov))
        if(TMath::Abs(d0z0[0]) > DCAxyCut || TMath::Abs(d0z0[1]) > DCAzCut) continue;
      //To be done : Add cuts to apply Chi2PerITSCls < 6 and N shared Cls ITS < 4
    }

    //+++++----Track properties
    Double_t dEdx =-999, fTPCnSigma=-999;
    dEdx = track->GetTPCsignal();
    fTPCnSigma = fpidResponse->NumberOfSigmasTPC(track, AliPID::kElectron);

    fTrkPt->Fill(track->Pt());
    fTrketa->Fill(track->Eta());
    fTrkphi->Fill(track->Phi());
    fdEdx->Fill(track->P(),dEdx);
    fTPCnsig->Fill(track->P(),fTPCnSigma);
    //+++++----Track matching to EMCAL
    Int_t EMCalIndex = -1;
    EMCalIndex = track->GetEMCALcluster();
    if(EMCalIndex < 0) continue;
    fHistPtMatch->Fill(track->Pt());

    AliVCluster *clustMatch = (AliVCluster*)fVevent->GetCaloCluster(EMCalIndex);
    if(!(clustMatch && clustMatch->IsEMCAL())) continue;

    //+++++----Properties of tracks matched to the EMCAL
    fEMCTrkMatch->Fill(clustMatch->GetTrackDx(),clustMatch->GetTrackDz());
    if(TMath::Abs(clustMatch->GetTrackDx())>0.05 || TMath::Abs(clustMatch->GetTrackDz())>0.05) continue;
    fEMCTrkPt->Fill(track->Pt());
    fEMCTrketa->Fill(track->Eta());
    fEMCTrkphi->Fill(track->Phi());
    fEMCTPCnsig->Fill(track->P(),fTPCnSigma);

    //+++---E/p distribution
    Double_t clustMatchE = clustMatch->E();
    Double_t eop = -1.0;
    Double_t m02 = -99999;
    if(track->P()>0)eop = clustMatchE/track->P();
    m02 =clustMatch->GetM02();

    if(track->Pt()>1.0){
      fHistNsigEop->Fill(eop,fTPCnSigma);
      fM20EovP->Fill(eop,clustMatch->GetM20());
      fM02EovP->Fill(eop,clustMatch->GetM02());
    }

    fHistEop->Fill(track->Pt(),eop);
    fM20->Fill(track->Pt(),clustMatch->GetM20());
    fM02->Fill(track->Pt(),clustMatch->GetM02());

    ///////////////////
    //Electron in jet//
    ///////////////////
    Bool_t fFlagPhotonicElec = kFALSE;
    Bool_t fTag_einJet = kFALSE, fLeadEle = kFALSE;
    //++++---electron selection cuts----
    if(fTPCnSigma > -1 && fTPCnSigma < 3 && eop>0.8 && eop<1.2 && m02 > 0.006 && m02 < 0.35)
    {
      fPtAllElec->Fill(track->Pt());

      //----+++Identify Non-HFE
      SelectPhotonicElectron(iTracks,track,fFlagPhotonicElec);

      //----+++Properties of all jets with electron
      PlotJetswithElec(track,fNthJetwEle,fPtJetwElec,fPtElecinJet,fPtElecLeadInJet,fEtaPhiJetwElec,fEtaJetwElecinEMC);
      if(fFlagPhotonicElec) //Photonic elec
        PlotJetswithElec(track,fNthJetwPEle,fPtJetwPElec,fPtPElecinJet,fPtPElecLeadInJet,fEtaPhiJetwPElec,fEtaJetwPElecinEMC);

      //----+++Leading 2 jets
      for(unsigned i=0; i<2;i++)
      {
        vector<fastjet::PseudoJet> constituents = fCSA->constituents(fJets[i]);
        if(fabs(fJets[i].eta()) > 0.9-fJetRcut)continue;

        int NumOfParJet = constituents.size();
        if(NumOfParJet<2)continue;   // reject D->ke candidates

        if(fJets[i].perp() < fMinJetPt) continue; //jets less than 10 GeV/c

        Double_t maxJetTrkpT = 0.0, maxJetTrkp=0.0, maxJetTrkpx=0.0, maxJetTrkpy=0.0, maxJetTrkpz=0.0;
        fTag_einJet = TagElectronJet(constituents,track,maxJetTrkpT, maxJetTrkp, maxJetTrkpx, maxJetTrkpy, maxJetTrkpz);
        if(!fTag_einJet) continue; //elec in Jet

        fPtLeadJetwElec->Fill(fJets[i].perp());
        fPtElecinLeadJet->Fill(track->Pt());
        if(track->Px() == maxJetTrkpx && track->Py() == maxJetTrkpy && track->Pz() == maxJetTrkpz) //ele is the leading track in jet
        {
          fLeadEle = kTRUE;
          fPtElecLeadInLeadJet->Fill(track->Pt());
        }

        //////////////////////////
        //Dphi Deta distribution//
        //////////////////////////
        fPtTrigJetwInclE->Fill(fJets[i].perp());
        ElecJetDeltaPhi(i,track, fSprsInclEleDphiDeta); //Inclusive elec

        if(fFlagPhotonicElec) //Photonic elec
        {
          fPtTrigJetwPhotoE->Fill(fJets[i].perp());
          ElecJetDeltaPhi(i,track, fSprsPhotEleDphiDeta);
        }

        if(!fFlagPhotonicElec) //Semi-Inclusive elec
        {
          fPtTrigJetwSemiInclE->Fill(fJets[i].perp());
          ElecJetDeltaPhi(i,track, fSprsSemiIncEleDphiDeta);
        }

        constituents.clear();
      }//jet loop
    }//electron selection
  }//track loop

  //Clear arrays, pointers and set BkgMedian=0
  Clear();

  PostData(1, fOutputList);
}
//________________________________________________________________________
void AliAnalysisTaskHFEDiJets::EMCalClusterInfo()
{
  /////////////////////////////
  //EMCAL cluster information//
  ////////////////////////////

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
      fHistClustE->Fill(clustE);
      fEMCClsEtaPhi->Fill(emceta,emcphi);
    }
  }
}
//________________________________________________________________________
void AliAnalysisTaskHFEDiJets::PlotJetswithElec(AliVTrack *Eletrack,TH1F *HisNthJetwEle,TH1F *HisPtJetwElec,TH1F *HisPtElecinJet,TH1F *HisPtElecLeadInJet,TH2F *HisEtaPhiJetwElec,TH1F *HisEtaJetwElecinEMC)
{
  /////////////////////////////////////////
  //Properties of all jets with electrons//
  /////////////////////////////////////////

  //------++++Jet signal Loop++++--------
  Double_t maxJetTrkpT = 0.0,  maxJetTrkp=0.0, maxJetTrkpx=0.0, maxJetTrkpy=0.0, maxJetTrkpz=0.0;
  Bool_t fEinJet = kFALSE;
  for(unsigned i = 0; i < fJets.size(); i++)
  {
    vector<fastjet::PseudoJet> constituents = fCSA->constituents(fJets[i]);
    if(fabs(fJets[i].eta()) > 0.9-fJetRcut)continue;

    int NumOfParJet = constituents.size();
    if(NumOfParJet<2)continue;   // reject D->ke candidates

    if(fJets[i].perp() < fMinJetPt) continue; //jets less than 15 GeV/c

    fEinJet = TagElectronJet(constituents,Eletrack,maxJetTrkpT, maxJetTrkp, maxJetTrkpx, maxJetTrkpy, maxJetTrkpz);

    if(!fEinJet) continue; //elec in Jet

    HisNthJetwEle->Fill(i);
    HisPtJetwElec->Fill(fJets[i].perp());
    HisPtElecinJet->Fill(Eletrack->Pt());
    if(Eletrack->Px() == maxJetTrkpx && Eletrack->Py() == maxJetTrkpy) //ele is the leading track in jet
      HisPtElecLeadInJet->Fill(Eletrack->Pt());

    HisEtaPhiJetwElec->Fill(fJets[i].eta(), fJets[i].phi());
    if(fJets[i].phi() > 1.396 && fJets[i].phi() < 3.142) //if jet in EMC phi acceptance
      HisEtaJetwElecinEMC->Fill(fJets[i].eta());

    constituents.clear();
  }
}
//________________________________________________________________________
void AliAnalysisTaskHFEDiJets::ElecJetDeltaPhi(unsigned iTrigJet, AliVTrack *Eletrack, THnSparse *fSprDphiDeta)
{
  ////////////////////////////
  //ElecJet-Jet DeltaPhi,Eta//
  ////////////////////////////
  Double_t jetPhi=-999, eleJetPhi=-999, elePhi=-999, Dphi=-999, pi = TMath::Pi();
  Double_t jetEta=-999, eleJetEta=-999, eleEta=-999, Deta=-999;
  Double_t fvalueDphi[7] = {-999,999,-999,-999,-999,-999, -999}; //TrigJetID, AssoJetID,ptElec, ptTrigJet, ptAssoJet, Dphi, Deta

  for(unsigned j= 0; j < fJets.size(); j++)
  {
    if(j == iTrigJet) continue;
    if(fabs(fJets[j].eta()) > 0.9-fJetRcut)continue;

    vector<fastjet::PseudoJet> constituents = fCSA->constituents(fJets[j]);
    int NumOfParJet = constituents.size();
    if(NumOfParJet<2)continue;

    if(fJets[j].perp() < fMinJetPt) continue; //jets less than 20 GeV/c

    jetPhi = fJets[j].phi();
    eleJetPhi = fJets[iTrigJet].phi();
    elePhi = Eletrack->Phi();

    jetEta = fJets[j].eta();
    eleJetEta = fJets[iTrigJet].eta();
    eleEta = Eletrack->Eta();
    Deta = eleJetEta - jetEta;

    Dphi = eleJetPhi - jetPhi;
    if (Dphi > 3*pi/2)
      Dphi = Dphi - 2*pi;
    if (Dphi < -pi/2)
      Dphi = Dphi + 2*pi;

    fvalueDphi[0] = iTrigJet;
    fvalueDphi[1] = j;
    fvalueDphi[2] = Eletrack->Pt();
    fvalueDphi[3] = fJets[iTrigJet].perp();
    fvalueDphi[4] = fJets[j].perp();
    fvalueDphi[5] = Dphi;
    fvalueDphi[6] = Deta;
    fSprDphiDeta->Fill(fvalueDphi);

    constituents.clear();
  }
}
//________________________________________________________________________
Bool_t AliAnalysisTaskHFEDiJets::TagElectronJet(vector<fastjet::PseudoJet> constJet, AliVTrack *Eletrack, Double_t &maxpT, Double_t &maxP, Double_t &maxPx, Double_t &maxPy, Double_t &maxPz)
{
  ///////////////////////////////
  //Check if electron is in jet//
  ///////////////////////////////

  Bool_t ElectronjetTag = kFALSE;
  Double_t maxpx = 0.0, maxpy =0.0, maxpz=0.0, maxp=0.0, maxpt=0.0, p=0.0;

  for (unsigned j = 0; j< constJet.size(); j++)
  {
    //max p,px,py,pz
    p= sqrt(pow(constJet[j].perp(),2)+pow(constJet[j].pz(),2));
    if(maxp < p)
    {
      maxp = p;
      maxpx = constJet[j].px(); maxpy = constJet[j].py(); maxpz = constJet[j].pz();
      maxpt = constJet[j].perp();
    }

    //track match elec
    if(Eletrack->Px() == constJet[j].px() && Eletrack->Py() == constJet[j].py() && Eletrack->Pz() == constJet[j].pz()) // electron in jet
    {
      ElectronjetTag = kTRUE;
    }
  }

  maxpT = maxpt;
  maxP = maxp;
  maxPx = maxpx; maxPy = maxpy; maxPz = maxpz;
  return ElectronjetTag;
}
//________________________________________________________________________
void AliAnalysisTaskHFEDiJets::InputTracksForJetReco(Double_t &fMaxpT)
{
  //////////////////////////////////////
  //Input tracks for jet reconstuction//
  //////////////////////////////////////

  Double_t maxpT= 0.15;
  for (Int_t iTracks = 0; iTracks < fVevent->GetNumberOfTracks(); iTracks++) {

    AliVParticle* Vtrack = 0x0;
    Vtrack  = fVevent->GetTrack(iTracks);

    if (!Vtrack) {
      printf("ERROR: Could not receive track %d\n", iTracks);
      continue;
    }
    AliVTrack *track = dynamic_cast<AliVTrack*>(Vtrack);
    AliAODTrack *atrack = dynamic_cast<AliAODTrack*>(Vtrack);

    if(fAOD)
//      if(!atrack->TestFilterBit(768)) continue; //hybrid track for LHC11h
      if(!atrack->IsHybridGlobalConstrainedGlobal()) continue;

    if(TMath::Abs(track->Eta())>0.9) continue;
    if(TMath::Abs(track->Pt()<0.15)) continue;

    Double_t phitrack = track->Phi();
    fJetInputTrPhi->Fill(phitrack);
    fJetInputTrEta->Fill(track->Eta());

    if(track->Pt() > maxpT) maxpT = track->Pt();
    fastjet::PseudoJet jInp(track->Px(),track->Py(),track->Pz(),track->P());
    fInputParticlesJet.push_back(jInp);

  }
  fMaxpT = maxpT;
}
//________________________________________________________________________
void AliAnalysisTaskHFEDiJets::DoJetReconstruction()
{
  //////////////////////
  //Jet reconstruction//
  //////////////////////

  fJetDef = new fastjet::JetDefinition(fastjet::antikt_algorithm, fJetRcut,fastjet::BIpt_scheme, fastjet::Best);

  fGhostArea = new fastjet::GhostedAreaSpec(0.9, 1, 0.005);
  fAreaDef = new fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts, *fGhostArea);

  fCSA = new fastjet::ClusterSequenceArea(fInputParticlesJet, *fJetDef, *fAreaDef);
  fJets = sorted_by_pt(fCSA->inclusive_jets());
}

//________________________________________________________________________
void AliAnalysisTaskHFEDiJets::DoBkgJetReconstruction()
{
  /////////////////////////////////
  //Background Jet reconstruction//
  /////////////////////////////////

  fBkgJetDef = new fastjet::JetDefinition(fastjet::kt_algorithm, 0.5, fastjet::BIpt_scheme, fastjet::Best);

  fBkgGhostArea = new fastjet::GhostedAreaSpec(0.9, 1, 0.005);
  fBkgAreaDef = new fastjet::AreaDefinition(fastjet::active_area_explicit_ghosts, *fBkgGhostArea);

  fBkgCSA = new fastjet::ClusterSequenceArea(fInputParticlesJet, *fBkgJetDef, *fBkgAreaDef);

  Double_t rapmim = -0.5, rapmax = 0.5, phimim = 0.0, phimax = 2*acos(-1.0);
  fRange = new fastjet::RangeDefinition(rapmim,rapmax,phimim,phimax); // default 0-2pi

  fBkgJets = sorted_by_pt(fBkgCSA->inclusive_jets());

  //-------Calculate median----
  Double_t median=0., sigma=0., meanarea=0.;
  fBkgCSA->get_median_rho_and_sigma(fBkgJets,*fRange, true, median, sigma, meanarea, true);
  fBgMedian_Brf->Fill(median);
  PlotRecoBkgJetProperties(fBgArea_Bfr, fBgPhi_Bfr, fBgEta_Bfr, fBgPt_Bfr, fBgP_Bfr);

  //Remove 2 leading jets from the Bkg jet array
  if(fBkgJets.size()>2) fBkgJets.erase(fBkgJets.begin(),fBkgJets.begin()+2);
  fBkgCSA->get_median_rho_and_sigma(fBkgJets, *fRange, true, fBkgMedian, fBkgSigma, fBkgMeanArea, true);
  fBgMedian_Afr->Fill(fBkgMedian);
  PlotRecoBkgJetProperties(fBgArea_Afr, fBgPhi_Afr, fBgEta_Afr, fBgPt_Afr, fBgP_Afr);
}
//________________________________________________________________________
void AliAnalysisTaskHFEDiJets::PlotRecoBkgJetProperties(TH1F* HisBkgArea, TH1F* HisBkgPhi, TH1F* HisBkgEta, TH1F* HisBkgPt, TH1F* HisBkgP)
{
  ////////////////////////////
  //Plot histos for Bkg jets//
  ////////////////////////////

  for(unsigned i = 0; i < fBkgJets.size(); i++)
  {
    if(fabs(fBkgJets[i].eta()) > 0.9-fJetRcut)continue;

    Double_t area = fBkgCSA->area(fBkgJets[i]);
    HisBkgArea->Fill(area);

    Double_t jetphi = fBkgJets[i].phi();
    Double_t jeteta = fBkgJets[i].eta();
    Double_t jetpt = fBkgJets[i].perp();
    Double_t jetmom = sqrt(pow(fBkgJets[i].perp(),2)+pow(fBkgJets[i].pz(),2));

    HisBkgPhi->Fill(jetphi);
    HisBkgEta->Fill(jeteta);
    HisBkgPt->Fill(jetpt);
    HisBkgP->Fill(jetmom);

  }
}
//________________________________________________________________________
void AliAnalysisTaskHFEDiJets::PlotRecoJetProperties(TH1F *HisSigArea,TH1F *HisSigPhi, TH1F *HisSigEta,TH1F *HisSigPt_BfrSub, TH1F *HisSigPt_Sub, TH1F *HisSigP, TH1F *HisSigR, TH2F *HisSigEtaPhi, TH1F *HisSigEtainEMC)
{
  ////////////////////////
  //Plot histos for jets//
  ////////////////////////

  for(unsigned i = 0; i < fJets.size(); i++)
  {
    if(fabs(fJets[i].eta()) > 0.9-fJetRcut)continue;
    if(fJets[i].perp() < fMinJetPt) continue;

    Double_t area = fCSA->area(fJets[i]);
    Double_t pt_sub = fJets[i].perp() - fBkgMedian*area;

    HisSigArea->Fill(area);
    HisSigPt_BfrSub->Fill(fJets[i].perp());
    HisSigPt_Sub->Fill(pt_sub);

    Double_t jetphi = fJets[i].phi();
    Double_t jeteta = fJets[i].eta();
    Double_t jetmom = sqrt(pow(fJets[i].perp(),2)+pow(fJets[i].pz(),2));
    Double_t jetR = sqrt(pow(jetphi,2)+pow(jeteta,2));


    HisSigPhi->Fill(jetphi);
    HisSigEta->Fill(jeteta);
    HisSigP->Fill(jetmom);

    HisSigR->Fill(jetR);
    HisSigEtaPhi->Fill(jeteta,jetphi);
    if(jetphi > 1.396 && jetphi < 3.142) //if jet is in EMC phi acceptance
      HisSigEtainEMC->Fill(jeteta);
  }
}
//________________________________________________________________________
void AliAnalysisTaskHFEDiJets::Clear(const Option_t */*opt*/)
{
  ////////////////////////////
  // Clear the input vectors//
  ////////////////////////////

  // Make sure done on every event if the instance is reused

  fInputParticlesJet.clear();
  fJets.clear();
  fBkgJets.clear();

  delete fJetDef; fJetDef=0;
  delete fGhostArea; fGhostArea=0;
  delete fAreaDef; fAreaDef=0;
  delete fCSA; fCSA=0;

  delete fBkgJetDef; fBkgJetDef=0;
  delete fBkgGhostArea; fBkgGhostArea=0;
  delete fBkgAreaDef; fBkgAreaDef=0;
  delete fBkgCSA; fBkgCSA=0;
  delete fRange; fRange=0;

  fBkgMedian = 0;
}
//________________________________________________________________________
void AliAnalysisTaskHFEDiJets::SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec)
{
  ///////////////////////////////
  //Photonic electron selection//
  ///////////////////////////////

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

    if(IsAODanalysis()) {
      AliAODTrack *atrackAsso = dynamic_cast<AliAODTrack*>(VtrackAsso);
      if(!atrackAsso) continue;
      if(!atrackAsso->TestFilterMask(AliAODTrack::kTrkTPCOnly)) continue;
      if(atrackAsso->GetTPCNcls() < 70) continue;
      if(!(atrackAsso->GetStatus()&AliESDtrack::kTPCrefit)) continue;
      if(!(atrackAsso->GetStatus()&AliESDtrack::kITSrefit)) continue;
    }

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

    if(mass<fInvmassCut && fFlagULS && !flagPhotonicElec){
      flagPhotonicElec = kTRUE;
    }
  }
  fFlagPhotonicElec = flagPhotonicElec;
}
//________________________________________________________________________
void AliAnalysisTaskHFEDiJets::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query

  fOutputList = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutputList) {
    printf("ERROR: Output list not available\n");
    return;
  }
}





