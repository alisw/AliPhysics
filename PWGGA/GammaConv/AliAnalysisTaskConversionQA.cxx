/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *									  *
 * Authors: Svein Lindal, Daniel Lohner					  *
 * Version 1.0								  *
 *									  *
 * Permission to use, copy, modify and distribute this software and its	  *
 * documentation strictly for non-commercial purposes is hereby granted	  *
 * without fee, provided that the above copyright notice appears in all	  *
 * copies and that both the copyright notice and this permission notice	  *
 * appear in the supporting documentation. The authors make no claims	  *
 * about the suitability of this software for any purpose. It is	  *
 * provided "as is" without express or implied warranty.		  *
 **************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// QA Task for V0 Reader V1
//---------------------------------------------
////////////////////////////////////////////////

#include "AliAnalysisTaskConversionQA.h"
#include "TChain.h"
#include "AliAnalysisManager.h"
#include "TParticle.h"
#include "TVectorF.h"
#include "AliPIDResponse.h"
#include "TFile.h"
#include "AliESDtrackCuts.h"
#include "AliAODMCParticle.h"
#include "AliAODMCHeader.h"


class iostream;

using namespace std;

ClassImp(AliAnalysisTaskConversionQA)

//________________________________________________________________________
AliAnalysisTaskConversionQA::AliAnalysisTaskConversionQA() : AliAnalysisTaskSE(),
   fV0Reader(NULL),
   fConversionGammas(NULL),
   fConversionCuts(NULL),
   fInputEvent(NULL),
   fNumberOfESDTracks(0),
   fMCEvent(NULL),
   fMCStack(NULL),
   fTreeQA(NULL),
   fIsHeavyIon(kFALSE),
   ffillTree(kFALSE),
   ffillHistograms(kFALSE),
   fOutputList(NULL),
   fTreeList(NULL),
   fESDList(NULL),
   hVertexZ(NULL),
   hNGoodESDTracks(NULL),
   hNV0Tracks(NULL),
   hNContributorsVertex(NULL),
   hITSClusterPhi(NULL),
   hGammaPt(NULL),
   hGammaPhi(NULL),
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
   hElectronNSigmaPiondEdxP(NULL),
   hElectronNSigmaITSP(NULL),
   hElectronNSigmaTOFP(NULL),
   hPositrondEdxP(NULL),
   hPositronITSdEdxP(NULL),
   hPositronTOFP(NULL),
   hPositronNSigmadEdxP(NULL),
   hPositronNSigmaPiondEdxP(NULL),
   hPositronNSigmaITSP(NULL),
   hPositronNSigmaTOFP(NULL),
//    hElecAsymP(NULL),
//    fTrueList(NULL),
//    hTrueResolutionR(NULL),
//    hTrueResolutionZ(NULL),
//    hTrueResolutionPhi(NULL),
//    hTrueGammaPt(NULL),
//    hTrueGammaPhi(NULL),
//    hTrueGammaEta(NULL),
//    hTrueGammaMass(NULL),
//    hTrueGammaChi2perNDF(NULL),
//    hTrueGammaPsiPair(NULL),
//    hTrueGammaQt(NULL),
//    hTrueGammaCosinePointingAngle(NULL),
//    hTrueGammaXY(NULL),
//    hTrueGammaZR(NULL),
//    hTrueElecPt(NULL),
//    hTrueElecEta(NULL),
//    hTrueElecPhi(NULL),
//    hTrueElecNfindableClsTPC(NULL),
//    hTruePosiNfindableClsTPC(NULL),
//    hTrueElecAsymP(NULL),
   fGammaPt(0),
   fGammaTheta(0),
   fGammaChi2NDF(0),
   fGammaPhotonProp(5),
   fGammaConvCoord(5),
   fDaughterProp(24),
   fKind(0),
   fIsMC(kFALSE),
   fnGammaCandidates(1),
   fMCStackPos(NULL),
   fMCStackNeg(NULL)
{

}

AliAnalysisTaskConversionQA::AliAnalysisTaskConversionQA(const char *name) : AliAnalysisTaskSE(name),
   fV0Reader(NULL),
   fConversionGammas(NULL),
   fConversionCuts(NULL),
   fInputEvent(NULL),
   fNumberOfESDTracks(0),
   fMCEvent(NULL),
   fMCStack(NULL),
   fTreeQA(NULL),
   fIsHeavyIon(kFALSE),
   ffillTree(kFALSE),
   ffillHistograms(kFALSE),
   fOutputList(NULL),
   fTreeList(NULL),
   fESDList(NULL),
   hVertexZ(NULL),
   hNGoodESDTracks(NULL),
   hNV0Tracks(NULL),
   hNContributorsVertex(NULL),
   hITSClusterPhi(NULL),
   hGammaPt(NULL),
   hGammaPhi(NULL),
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
   hElectronNSigmaPiondEdxP(NULL),
   hElectronNSigmaITSP(NULL),
   hElectronNSigmaTOFP(NULL),
   hPositrondEdxP(NULL),
   hPositronITSdEdxP(NULL),
   hPositronTOFP(NULL),
   hPositronNSigmadEdxP(NULL),
   hPositronNSigmaPiondEdxP(NULL),
   hPositronNSigmaITSP(NULL),
   hPositronNSigmaTOFP(NULL),
//    hGammaXY(NULL),
//    hGammaZR(NULL),
//    hElecAsymP(NULL),
//    fTrueList(NULL),
//    hTrueResolutionR(NULL),
//    hTrueResolutionZ(NULL),
//    hTrueResolutionPhi(NULL),
//    hTrueGammaPt(NULL),
//    hTrueGammaPhi(NULL),
//    hTrueGammaEta(NULL),
//    hTrueGammaMass(NULL),
//    hTrueGammaChi2perNDF(NULL),
//    hTrueGammaPsiPair(NULL),
//    hTrueGammaQt(NULL),
//    hTrueGammaCosinePointingAngle(NULL),
//    hTrueGammaXY(NULL),
//    hTrueGammaZR(NULL),
//    hTrueElecPt(NULL),
//    hTrueElecEta(NULL),
//    hTrueElecPhi(NULL),
//    hTrueElecNfindableClsTPC(NULL),
//    hTruePosiNfindableClsTPC(NULL),
//    hTrueElecAsymP(NULL),
   fGammaPt(0),
   fGammaTheta(0),
   fGammaChi2NDF(0),
   fGammaPhotonProp(5),
   fGammaConvCoord(5),
   fDaughterProp(24),
   fKind(0),
   fIsMC(kFALSE),
   fnGammaCandidates(1),
   fMCStackPos(NULL),
   fMCStackNeg(NULL)
{
   // Default constructor

   DefineInput(0, TChain::Class());
   DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskConversionQA::~AliAnalysisTaskConversionQA()
{
   // default deconstructor
  
}
//________________________________________________________________________
void AliAnalysisTaskConversionQA::UserCreateOutputObjects()
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
         
      fESDList = new TList();
      fESDList->SetOwner(kTRUE);
      fESDList->SetName("ESD QA");
      fOutputList->Add(fESDList);

      hVertexZ = new TH1F("Vertex_Z","Vertex_Z",300,-15,15);
      fESDList->Add(hVertexZ);
      hNContributorsVertex = new TH1I("ContrVertex_Z","ContrVertex_Z",3000,0,3000);
      fESDList->Add(hNContributorsVertex);
      if(fIsHeavyIon) hNGoodESDTracks = new TH1I("GoodESDTracks","GoodESDTracks",4000,0,4000);
      else hNGoodESDTracks = new TH1I("GoodESDTracks","GoodESDTracks",200,0,200);
      fESDList->Add(hNGoodESDTracks);
      if(fIsHeavyIon) hNV0Tracks = new TH1I("V0 Multiplicity","V0 Multiplicity",30000,0,30000);
      else hNV0Tracks = new TH1I("V0 Multiplicity","V0 Multiplicity",2000,0,2000);
      fESDList->Add(hNV0Tracks);

      hITSClusterPhi = new TH2F("ITSClusterPhi","hITSClusterPhi",72,0,2*TMath::Pi(),7,0,7);
      fESDList->Add(hITSClusterPhi);
      hGammaPt = new TH1F("Gamma_Pt","Gamma_Pt",250,0,25);
      fESDList->Add(hGammaPt);
      hGammaPhi = new TH1F("Gamma_Phi","Gamma_Phi",360,0,2*TMath::Pi());
      fESDList->Add(hGammaPhi);
      hGammaEta = new TH1F("Gamma_Eta","Gamma_Eta",600,-1.5,1.5);
      fESDList->Add(hGammaEta);
      hGammaChi2perNDF = new TH1F("Gamma_Chi2perNDF","Gamma_Chi2perNDF",500,0,100);
      fESDList->Add(hGammaChi2perNDF);
      hGammaPsiPair = new TH1F("Gamma_PsiPair","Gamma_PsiPair",500,0,2);
      fESDList->Add(hGammaPsiPair);
      hGammaArmenteros = new TH2F("Gamma_Armenteros","Gamma_Armenteros",200,-1,1,400,0,0.1);
      fESDList->Add(hGammaArmenteros);
      hGammaCosinePointingAngle = new TH1F("Gamma_CosinePointingAngle","Gamma_CosinePointingAngle",1000,-1.,1.);
      fESDList->Add(hGammaCosinePointingAngle);
      hGammaInvMass = new TH1F( "Gamma_InvMass","",200, 0, 0.2);
      fESDList->Add(hGammaInvMass);
      hElecPt = new TH2F("Electron_Positron_Pt","Electron_Positron_Pt",250,0,25,250,0,25);
      fESDList->Add(hElecPt);
      hElecEta = new TH2F("Electron_Positron_Eta","Electron_Positron_Eta",600,-1.5,1.5,600,-1.5,1.5);
      fESDList->Add(hElecEta);
      hElecPhi = new TH2F("Electron_Positron_Phi","Electron_Positron_Phi",360,0,2*TMath::Pi(),360,0,2*TMath::Pi());
      fESDList->Add(hElecPhi);
      hElecClsTPC = new TH1F("Electron_ClusterTPC","Electron_ClusterTPC",200,0,200);
      fESDList->Add(hElecClsTPC);
      hPosiClsTPC = new TH1F("Positron_ClusterTPC","Positron_ClusterTPC",200,0,200);
      fESDList->Add(hPosiClsTPC);
      
      hElecNfindableClsTPC = new TH1F("Electron_findableClusterTPC","Electron_findableClusterTPC",100,0,1);
      fESDList->Add(hElecNfindableClsTPC);
      hPosiNfindableClsTPC = new TH1F("Positron_findableClusterTPC","Positron_findableClusterTPC",100,0,1);
      fESDList->Add(hPosiNfindableClsTPC);
      
      hElectrondEdxP =  new TH2F("Electron_dEdx_P","Electron_dEdx_P",100, 0.05, 20, 200, 0, 200);
      SetLogBinningXTH2(hElectrondEdxP);
      fESDList->Add(hElectrondEdxP);
      hPositrondEdxP =  new TH2F("Positron_dEdx_P","Positron_dEdx_P",100, 0.05, 20, 200, 0, 200);
      SetLogBinningXTH2(hPositrondEdxP);
      fESDList->Add(hPositrondEdxP);
      hElectronNSigmadEdxP =  new TH2F("Electron_NSigmadEdx_P","Electron_NSigmadEdx_P",100, 0.05, 20, 200, -10, 10);  
      SetLogBinningXTH2(hElectronNSigmadEdxP);
      fESDList->Add(hElectronNSigmadEdxP);
      hPositronNSigmadEdxP =  new TH2F("Positron_NSigmadEdx_P","Positron_NSigmadEdx_P",100, 0.05, 20, 200, -10, 10);
      SetLogBinningXTH2(hPositronNSigmadEdxP);
      fESDList->Add(hPositronNSigmadEdxP);
      hElectronNSigmaPiondEdxP =  new TH2F("Electron_NSigmaPiondEdx_P","Electron_NSigmaPiondEdx_P",100, 0.05, 20, 200, -10, 10);  
      SetLogBinningXTH2(hElectronNSigmaPiondEdxP);
      fESDList->Add(hElectronNSigmaPiondEdxP);
      hPositronNSigmaPiondEdxP =  new TH2F("Positron_NSigmaPiondEdx_P","Positron_NSigmaPiondEdx_P",100, 0.05, 20, 200, -10, 10);
      SetLogBinningXTH2(hPositronNSigmaPiondEdxP);
      fESDList->Add(hPositronNSigmaPiondEdxP);
      
      hElectronTOFP =  new TH2F("Electron_TOF_P","Electron_TOF_P",100, 0.05, 20, 600, -1000, 29000);
      SetLogBinningXTH2(hElectronTOFP);
      fESDList->Add(hElectronTOFP);
      hPositronTOFP =  new TH2F("Positron_TOF_P","Positron_TOF_P",100, 0.05, 20, 600, -1000, 29000);
      SetLogBinningXTH2(hPositronTOFP);
      fESDList->Add(hPositronTOFP);
      hElectronNSigmaTOFP =  new TH2F("Electron_NSigmaTOF_P","Electron_NSigmaTOF_P",100, 0.05, 20, 200, -10, 10);
      SetLogBinningXTH2(hElectronNSigmaTOFP);
      fESDList->Add(hElectronNSigmaTOFP);
      hPositronNSigmaTOFP =  new TH2F("Positron_NSigmaTOF_P","Positron_NSigmaTOF_P",100, 0.05, 20, 200, -10, 10);
      SetLogBinningXTH2(hPositronNSigmaTOFP);
      fESDList->Add(hPositronNSigmaTOFP);
      
      hElectronITSdEdxP  =  new TH2F("Electron_ITSdEdx_P","Electron_ITSdEdx_P",100, 0.05, 20, 200, 0, 200);
      SetLogBinningXTH2(hElectronITSdEdxP);
      fESDList->Add(hElectronITSdEdxP);
      hPositronITSdEdxP =  new TH2F("Positron_ITSdEdx_P","Positron_ITSdEdx_P",100, 0.05, 20, 200, 0, 200);
      SetLogBinningXTH2(hPositronITSdEdxP);
      fESDList->Add(hPositronITSdEdxP);
      hElectronNSigmaITSP =  new TH2F("Electron_NSigmaITS_P","Electron_NSigmaITS_P",100, 0.05, 20, 200, -10, 10);
      SetLogBinningXTH2(hElectronNSigmaITSP);
      fESDList->Add(hElectronNSigmaITSP);
      hPositronNSigmaITSP =  new TH2F("Positron_NSigmaITS_P","Positron_NSigmaITS_P",100, 0.05, 20, 200, -10, 10);
      SetLogBinningXTH2(hPositronNSigmaITSP);
      fESDList->Add(hPositronNSigmaITSP);
      
      
      
//       hGammaXY = new TH2F("Gamma_ConversionPoint_XY","Gamma_ConversionPoint_XY",960,-120,120,960,-120,120);
//       fESDList->Add(hGammaXY);
//       hGammaZR= new TH2F("Gamma_ConversionPoint_ZR","Gamma_ConversionPoint_ZR",1200,-150,150,480,0,120);
//       fESDList->Add(hGammaZR);


//       hElecAsymP = new TH2F("Electron_Asym_vs_P", "Electron_Asym_vs_P",200,0.,20.,200,0.,1.); 
//       fESDList->Add(hElecAsymP);

//       if(fIsMC){
//          fTrueList = new TList();
//          fTrueList->SetOwner(kTRUE);
//          fTrueList->SetName("True QA");
//          fOutputList->Add(fTrueList);
// 
//          hTrueResolutionR = new TH2F("True_ConversionPointResolution_R","True_ConversionPointResolution_R",240,0,120,200,-20,20);
//          fTrueList->Add(hTrueResolutionR);
//          hTrueResolutionZ = new TH2F("True_ConversionPointResolution_Z","True_ConversionPointResolution_Z",480,-120,120,200,-20,20);
//          fTrueList->Add(hTrueResolutionZ);
//          hTrueResolutionPhi = new TH2F("True_ConversionPointResolution_Phi","True_ConversionPointResolution_Phi",360,0,2*TMath::Pi(),200,-TMath::Pi()/30., TMath::Pi()/30.);
//          fTrueList->Add(hTrueResolutionPhi);
// 
//          hTrueGammaPt = new TH1F("True_Gamma_Pt","True_Gamma_Pt",250,0,25);
//          fTrueList->Add(hTrueGammaPt);
//          hTrueGammaPhi = new TH1F("True_Gamma_Phi","True_Gamma_Phi",360,0,2*TMath::Pi());
//          fTrueList->Add(hTrueGammaPhi);
//          hTrueGammaEta = new TH1F("True_Gamma_Eta","True_Gamma_Eta",600,-1.5,1.5);
//          fTrueList->Add(hTrueGammaEta);
//          hTrueGammaMass = new TH1F("True_Gamma_Mass","True_Gamma_Mass",1000,0,0.3);
//          fTrueList->Add(hTrueGammaMass);
//          hTrueGammaChi2perNDF = new TH1F("True_Gamma_Chi2perNDF","True_Gamma_Chi2perNDF",500,0,100);
//          fTrueList->Add(hTrueGammaChi2perNDF);
//          hTrueGammaPsiPair = new TH1F("True_Gamma_PsiPair","True_Gamma_PsiPair",500,0,2);
//          fTrueList->Add(hTrueGammaPsiPair);
//          hTrueGammaQt = new TH1F("True_Gamma_Qt","True_Gamma_Qt",400,0,0.1);
//          fTrueList->Add(hTrueGammaQt);
//          hTrueGammaCosinePointingAngle = new TH1F("True_Gamma_CosinePointingAngle","True_Gamma_CosinePointingAngle",900,0.7,1.);
//          fTrueList->Add(hTrueGammaCosinePointingAngle);
//          hTrueGammaXY = new TH2F("True_Gamma_ConversionPoint_XY","True_Gamma_ConversionPoint_XY",960,-120,120,960,-120,120);
//          fTrueList->Add(hTrueGammaXY);
//          hTrueGammaZR= new TH2F("TrueGamma_ConversionPoint_ZR","TrueGamma_ConversionPoint_ZR",1200,-150,150,480,0,120);
//          fTrueList->Add(hTrueGammaZR);
// 
//          hTrueElecPt = new TH2F("True_Electron_Positron_Pt","True_Electron_Positron_Pt",250,0,25,250,0,25);
//          fTrueList->Add(hTrueElecPt);
//          hTrueElecEta = new TH2F("True_Electron_Positron_Eta","True_Electron_Positron_Eta",600,-1.5,1.5,600,-1.5,1.5);
//          fTrueList->Add(hTrueElecEta);
//          hTrueElecPhi = new TH2F("True_Electron_Positron_Phi","True_Electron_Positron_Phi",360,0,2*TMath::Pi(),360,0,2*TMath::Pi());
//          fTrueList->Add(hTrueElecPhi);
//          hTrueElecNfindableClsTPC = new TH1F("True_Electron_findableClusterTPC","True_Electron_findableClusterTPC",100,0,1);
//          fTrueList->Add(hTrueElecNfindableClsTPC);
//          hTruePosiNfindableClsTPC = new TH1F("True_Positron_findableClusterTPC","True_Positron_findableClusterTPC",100,0,1);
//          fTrueList->Add(hTruePosiNfindableClsTPC);
// 				 hTrueElecAsymP = new TH2F("True_Electron_Asym_vs_P", "True_Electron_Asym_vs_P",200,0.,20.,200,0.,1.); 
// 				 fTrueList->Add(hTrueElecAsymP);
//       }
      if(fConversionCuts->GetCutHistograms()){
         fOutputList->Add(fConversionCuts->GetCutHistograms());
      }
   }
   
   if(ffillTree){
      fTreeList = new TList();
      fTreeList->SetOwner(kTRUE);
      fTreeList->SetName("TreeList");
      fOutputList->Add(fTreeList);

      fTreeQA = new TTree("PhotonQA","PhotonQA");   
     
      fTreeQA->Branch("daughterProp",&fDaughterProp);
      fTreeQA->Branch("recCords",&fGammaConvCoord);
      fTreeQA->Branch("photonProp",&fGammaPhotonProp);
      fTreeQA->Branch("pt",&fGammaPt,"fGammaPt/F");
      fTreeQA->Branch("theta",&fGammaTheta,"fGammaTheta/F");
      fTreeQA->Branch("chi2ndf",&fGammaChi2NDF,"fGammaChi2NDF/F");
      if (fIsMC) {
         fTreeQA->Branch("kind",&fKind,"fKind/b");
      }   
      fTreeList->Add(fTreeQA);

   }

   fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1");
   
   PostData(1, fOutputList);
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskConversionQA::Notify()
{
   if(!fConversionCuts->GetDoEtaShift()) return kTRUE;; // No Eta Shift requested, continue
      
   if(fConversionCuts->GetEtaShift() == 0.0){ // Eta Shift requested but not set, get shift automatically
      fConversionCuts->GetCorrectEtaShiftFromPeriod(fV0Reader->GetPeriodName());
      fConversionCuts->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
      return kTRUE;
   }
   else{
      printf(" Gamma Conversion QA Task %s :: Eta Shift Manually Set to %f \n\n",
             (fConversionCuts->GetCutNumber()).Data(),fConversionCuts->GetEtaShift());
      fConversionCuts->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once
   }
   
   return kTRUE;
}
//________________________________________________________________________
void AliAnalysisTaskConversionQA::UserExec(Option_t *){



   Int_t eventQuality = ((AliConversionCuts*)fV0Reader->GetConversionCuts())->GetEventQuality();
   if(eventQuality != 0){// Event Not Accepted
      return;
   }
   fInputEvent = InputEvent();
   if(fIsMC) fMCEvent = MCEvent();
   if(fMCEvent && fInputEvent->IsA()==AliESDEvent::Class()){ fMCStack = fMCEvent->Stack(); }

   Int_t eventNotAccepted =
      fConversionCuts->IsEventAcceptedByConversionCut(fV0Reader->GetConversionCuts(),fInputEvent,fMCEvent,fIsHeavyIon);
   if(eventNotAccepted) return; // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1

   fConversionGammas=fV0Reader->GetReconstructedGammas();

   if(fMCEvent){
      if(fConversionCuts->GetSignalRejection() != 0){
         if(fInputEvent->IsA()==AliESDEvent::Class()){
               fConversionCuts->GetNotRejectedParticles(fConversionCuts->GetSignalRejection(),
                                                        fConversionCuts->GetAcceptedHeader(),
                                                        fMCEvent);
         }
         else if(fInputEvent->IsA()==AliAODEvent::Class()){
            fConversionCuts->GetNotRejectedParticles(fConversionCuts->GetSignalRejection(),
                                                      fConversionCuts->GetAcceptedHeader(),
                                                      fInputEvent);
         }
      }
   }

   if(ffillHistograms){
      CountTracks();
      hVertexZ->Fill(fInputEvent->GetPrimaryVertex()->GetZ());
      hNContributorsVertex->Fill(fConversionCuts->GetNumberOfContributorsVtx(fInputEvent));
      hNGoodESDTracks->Fill(fNumberOfESDTracks);
      hNV0Tracks->Fill(fInputEvent->GetVZEROData()->GetMTotV0A()+fInputEvent->GetVZEROData()->GetMTotV0C());
   }

    if(fMCEvent && fInputEvent->IsA()==AliAODEvent::Class())
      RelabelAODPhotonCandidates(kTRUE);    // In case of AODMC relabeling MC
   
   for(Int_t firstGammaIndex=0;firstGammaIndex<fConversionGammas->GetEntriesFast();firstGammaIndex++){
      AliAODConversionPhoton *gamma=dynamic_cast<AliAODConversionPhoton*>(fConversionGammas->At(firstGammaIndex));
      if (gamma==NULL) continue;
      if(fMCEvent && fConversionCuts->GetSignalRejection() != 0){
         if(!fConversionCuts->IsParticleFromBGEvent(gamma->GetMCLabelPositive(), fMCStack, fInputEvent))
            continue;
         if(!fConversionCuts->IsParticleFromBGEvent(gamma->GetMCLabelNegative(), fMCStack, fInputEvent))
            continue;
      }
      if(!fConversionCuts->PhotonIsSelected(gamma,fInputEvent)){
         continue;
      }

      if(ffillTree) ProcessQATree(gamma);
      if(ffillHistograms) ProcessQA(gamma);
   }
   if(fMCEvent && fInputEvent->IsA()==AliAODEvent::Class())
      RelabelAODPhotonCandidates(kFALSE);    // In case of AODMC relabeling MC
      
   PostData(1, fOutputList);
}


///________________________________________________________________________
void AliAnalysisTaskConversionQA::ProcessQATree(AliAODConversionPhoton *gamma){

   // Fill Histograms for QA and MC
   AliVEvent* event = (AliVEvent*) InputEvent();
   
   
   
   AliPIDResponse* pidResonse = ((AliConversionCuts*)fV0Reader->GetConversionCuts())->GetPIDResponse();

   fGammaPt = gamma->GetPhotonPt();
   
   fGammaTheta = gamma->Theta();
   fGammaChi2NDF = gamma->GetChi2perNDF();
   
   fGammaPhotonProp(0)  = gamma->GetArmenterosQt();
   fGammaPhotonProp(1)  = gamma->GetArmenterosAlpha();
   fGammaPhotonProp(2)  = gamma->GetPsiPair();
   fGammaPhotonProp(3) = fConversionCuts->GetCosineOfPointingAngle(gamma,event);
	fGammaPhotonProp(4) = gamma->GetMass();
   
   fGammaConvCoord(0) = gamma->GetConversionX();
   fGammaConvCoord(1) = gamma->GetConversionY();
   fGammaConvCoord(2) = gamma->GetConversionZ();
   fGammaConvCoord(3) = gamma->GetConversionRadius();
   fGammaConvCoord(4) = gamma->GetPhotonPhi();
   
   AliVTrack * negTrack = fConversionCuts->GetTrack(event, gamma->GetTrackLabelNegative());
   AliVTrack * posTrack = fConversionCuts->GetTrack(event, gamma->GetTrackLabelPositive());

   
   if(!negTrack||!posTrack)return;

   fKind = 9;
   if(fMCEvent && fInputEvent->IsA()==AliESDEvent::Class()){
      fKind = IsTruePhotonESD(gamma);
   } else if (fMCEvent && fInputEvent->IsA()==AliAODEvent::Class()){
      fKind = IsTruePhotonAOD(gamma);   
   }

   fDaughterProp(0) =  posTrack->Pt();
   fDaughterProp(7) =  negTrack->Pt();
   fDaughterProp(1) =  posTrack->Theta();
   fDaughterProp(8) =  negTrack->Theta();
//    Double32_t signalPos[4] = {0,0,0,0};
//    Char_t nclPos[3];
//    Char_t nrowsPos[3];
//    if (posTrack->GetTPCdEdxInfo()) {
//       posTrack->GetTPCdEdxInfo()->GetTPCSignalRegionInfo(signalPos,nclPos,nrowsPos);
//       fDaughterProp(2) =  signalPos[0];
//       fDaughterProp(14) =  signalPos[1];
//       fDaughterProp(16) =  signalPos[2];
//    } else {
//       fDaughterProp(2) =  posTrack->GetTPCsignal();
//       fDaughterProp(14) =  0;
//       fDaughterProp(16) =  0;
//    }
//    Double32_t signalNeg[4] = {0,0,0,0};
//    Char_t nclNeg[3];
//    Char_t nrowsNeg[3];
//    if (negTrack->GetTPCdEdxInfo()) {
//       negTrack->GetTPCdEdxInfo()->GetTPCSignalRegionInfo(signalNeg,nclNeg,nrowsNeg);
//       fDaughterProp(9) =  signalNeg[0];
//       fDaughterProp(15) =  signalNeg[1];
//       fDaughterProp(17) =  signalNeg[2];
//    } else {
//       fDaughterProp(9) =  negTrack->GetTPCsignal();
//       fDaughterProp(15) =  0;
//       fDaughterProp(17) =  0;
//    }
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
      Double_t timesPos[5];
      posTrack->GetIntegratedTimes(timesPos);
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
      Double_t timesNeg[5];
      negTrack->GetIntegratedTimes(timesNeg);
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
   
   if (fTreeQA){
      fTreeQA->Fill();
   }
}

//_____________________________________________________________________________________________________
void AliAnalysisTaskConversionQA::ProcessQA(AliAODConversionPhoton *gamma){

   AliPIDResponse* pidResonse = ((AliConversionCuts*)fV0Reader->GetConversionCuts())->GetPIDResponse();
   
   // Fill Histograms for QA and MC

   hGammaPt->Fill(gamma->GetPhotonPt());
   hGammaPhi->Fill(gamma->GetPhotonPhi());
   hGammaEta->Fill(gamma->Eta());
   hGammaChi2perNDF->Fill(gamma->GetChi2perNDF());
   hGammaPsiPair->Fill(gamma->GetPsiPair());
   hGammaArmenteros->Fill(gamma->GetArmenterosAlpha(),gamma->GetArmenterosQt());
   hGammaCosinePointingAngle->Fill(fConversionCuts->GetCosineOfPointingAngle(gamma,fInputEvent));
   hGammaInvMass->Fill(gamma->GetMass());
//    hGammaXY->Fill(gamma->GetConversionX(),gamma->GetConversionY());
//    hGammaZR->Fill(gamma->GetConversionZ(),gamma->GetConversionRadius());

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
   hElectronNSigmaPiondEdxP->Fill(negTrack->P() ,pidResonse->NumberOfSigmasTPC(negTrack, AliPID::kPion));
   hPositrondEdxP->Fill(posTrack->P() ,posTrack->GetTPCsignal());
   hPositronNSigmadEdxP->Fill(posTrack->P() ,pidResonse->NumberOfSigmasTPC(posTrack, AliPID::kElectron));
   hPositronNSigmaPiondEdxP->Fill(posTrack->P() ,pidResonse->NumberOfSigmasTPC(posTrack, AliPID::kPion));
   
   //TOF signal
   if((negTrack->GetStatus() & AliESDtrack::kTOFpid) && !(negTrack->GetStatus() & AliESDtrack::kTOFmismatch)){
      Double_t t0neg = pidResonse->GetTOFResponse().GetStartTime(negTrack->P());
      Double_t timesNeg[5];
      negTrack->GetIntegratedTimes(timesNeg);
      Double_t TOFsignalNeg = negTrack->GetTOFsignal();
      Double_t dTneg = TOFsignalNeg - t0neg - timesNeg[0];
      hElectronTOFP->Fill(negTrack->P() ,dTneg);
      hElectronNSigmaTOFP->Fill(negTrack->P() ,pidResonse->NumberOfSigmasTOF(negTrack, AliPID::kElectron));
   }
   if((posTrack->GetStatus() & AliESDtrack::kTOFpid) && !(posTrack->GetStatus() & AliESDtrack::kTOFmismatch)){
      Double_t t0pos = pidResonse->GetTOFResponse().GetStartTime(posTrack->P());
      Double_t timesPos[5];
      posTrack->GetIntegratedTimes(timesPos);
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

   
//    if(gamma->P()!=0){
// 		 hElecAsymP->Fill(gamma->P(),negTrack->P()/gamma->P());
// 	 }
   // hElecNfindableClsTPC->Fill((Float_t)negTrack->GetNcls(1)/(Float_t)negTrack->GetTPCNclsF());
   // hPosiNfindableClsTPC->Fill((Float_t)posTrack->GetNcls(1)/(Float_t)posTrack->GetTPCNclsF());
//    if(fMCEvent && fInputEvent->IsA()==AliESDEvent::Class()){
//       ProcessTrueQAESD(gamma,(AliESDtrack*)negTrack,(AliESDtrack*)posTrack);
//    } else if(fMCEvent && fInputEvent->IsA()==AliAODEvent::Class()){
//         ProcessTrueQAAOD(gamma,(AliAODTrack*)negTrack,(AliAODTrack*)posTrack);
//    }
}

// //________________________________________________________________________
// void AliAnalysisTaskConversionQA::ProcessTrueQAESD(AliAODConversionPhoton *TruePhotonCandidate, AliESDtrack *elec, AliESDtrack *posi)
// {
// 
//    if(IsTruePhotonESD(TruePhotonCandidate)!=0 && IsTruePhotonESD(TruePhotonCandidate)!=5  ) return;
// 
//    TParticle *negDaughter = TruePhotonCandidate->GetNegativeMCDaughter(fMCStack);
//    TParticle *mcPhoton = TruePhotonCandidate->GetMCParticle(fMCStack);
   // True Photon
//    hTrueResolutionR->Fill(TruePhotonCandidate->GetConversionRadius(),
//                            TruePhotonCandidate->GetConversionRadius()-negDaughter->R());
//    hTrueResolutionZ->Fill(TruePhotonCandidate->GetConversionZ(),
//                            TruePhotonCandidate->GetConversionZ()-negDaughter->Vz());
//    hTrueResolutionPhi->Fill(TruePhotonCandidate->Phi(),
//                              TruePhotonCandidate->Phi()-mcPhoton->Phi());
//    hTrueGammaPt->Fill(TruePhotonCandidate->Pt());
//    hTrueGammaPhi->Fill(TruePhotonCandidate->Phi());
//    hTrueGammaEta->Fill(TruePhotonCandidate->Eta());
//    hTrueGammaMass->Fill(TruePhotonCandidate->GetMass());
//    hTrueGammaChi2perNDF->Fill(TruePhotonCandidate->GetChi2perNDF());
//    hTrueGammaPsiPair->Fill(TruePhotonCandidate->GetPsiPair());
//    hTrueGammaQt->Fill(TruePhotonCandidate->GetArmenterosQt());
//    hTrueGammaCosinePointingAngle->Fill(fConversionCuts->GetCosineOfPointingAngle(TruePhotonCandidate,fInputEvent));
//    hTrueGammaXY->Fill(TruePhotonCandidate->GetConversionX(),TruePhotonCandidate->GetConversionY());
//    hTrueGammaZR->Fill(TruePhotonCandidate->GetConversionZ(),TruePhotonCandidate->GetConversionRadius());
//    hTrueElecPt->Fill(elec->Pt(),posi->Pt());
//    hTrueElecEta->Fill(elec->Eta(),posi->Eta());
//    hTrueElecPhi->Fill(elec->Phi(),posi->Phi());
//    hTrueElecNfindableClsTPC
//       ->Fill((Float_t)elec->GetTPCClusterInfo(2,0,fConversionCuts->GetFirstTPCRow(TruePhotonCandidate->GetConversionRadius())));
//    hTruePosiNfindableClsTPC
//       ->Fill((Float_t)posi->GetTPCClusterInfo(2,0,fConversionCuts->GetFirstTPCRow(TruePhotonCandidate->GetConversionRadius())));
//    if(TruePhotonCandidate->P()!=0){
// 		 hTrueElecAsymP->Fill(TruePhotonCandidate->P(),elec->P()/TruePhotonCandidate->P());
// 	 }

   // hTrueElecNfindableClsTPC->Fill((Float_t)elec->GetNcls(1)/(Float_t)elec->GetTPCNclsF());
   // hTruePosiNfindableClsTPC->Fill((Float_t)posi->GetNcls(1)/(Float_t)posi->GetTPCNclsF());
// }
// 
// //________________________________________________________________________
// void AliAnalysisTaskConversionQA::ProcessTrueQAAOD(AliAODConversionPhoton *TruePhotonCandidate, AliAODTrack *elec, AliAODTrack *posi)
// {
// 
//    if(IsTruePhotonAOD(TruePhotonCandidate)!=0 && IsTruePhotonESD(TruePhotonCandidate)!=5 ) return;
// 
//    TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
//    AliAODMCParticle *negDaughter = (AliAODMCParticle*) AODMCTrackArray->At(TruePhotonCandidate->GetMCLabelNegative());
//    AliAODMCParticle *mcPhoton = (AliAODMCParticle*) AODMCTrackArray->At(negDaughter->GetMother());
   // True Photon
//    hTrueResolutionR->Fill(TruePhotonCandidate->GetConversionRadius(),
//                            TruePhotonCandidate->GetConversionRadius()-(TMath::Sqrt(negDaughter->Xv()*negDaughter->Xv()+negDaughter->Yv()*negDaughter->Yv())));
//    hTrueResolutionZ->Fill(TruePhotonCandidate->GetConversionZ(),
//                            TruePhotonCandidate->GetConversionZ()-negDaughter->Zv());
//    hTrueResolutionPhi->Fill(TruePhotonCandidate->Phi(),
//                              TruePhotonCandidate->Phi()-mcPhoton->Phi());
//    hTrueGammaPt->Fill(TruePhotonCandidate->Pt());
//    hTrueGammaPhi->Fill(TruePhotonCandidate->Phi());
//    hTrueGammaEta->Fill(TruePhotonCandidate->Eta());
//    hTrueGammaMass->Fill(TruePhotonCandidate->GetMass());
//    hTrueGammaChi2perNDF->Fill(TruePhotonCandidate->GetChi2perNDF());
//    hTrueGammaPsiPair->Fill(TruePhotonCandidate->GetPsiPair());
//    hTrueGammaQt->Fill(TruePhotonCandidate->GetArmenterosQt());
//    hTrueGammaCosinePointingAngle->Fill(fConversionCuts->GetCosineOfPointingAngle(TruePhotonCandidate,fInputEvent));
//    hTrueGammaXY->Fill(TruePhotonCandidate->GetConversionX(),TruePhotonCandidate->GetConversionY());
//    hTrueGammaZR->Fill(TruePhotonCandidate->GetConversionZ(),TruePhotonCandidate->GetConversionRadius());
//    hTrueElecPt->Fill(elec->Pt(),posi->Pt());
//    hTrueElecEta->Fill(elec->Eta(),posi->Eta());
//    hTrueElecPhi->Fill(elec->Phi(),posi->Phi());
//    hTrueElecNfindableClsTPC
//       ->Fill((Float_t)elec->GetTPCClusterInfo(2,0,fConversionCuts->GetFirstTPCRow(TruePhotonCandidate->GetConversionRadius())));
//    hTruePosiNfindableClsTPC
//       ->Fill((Float_t)posi->GetTPCClusterInfo(2,0,fConversionCuts->GetFirstTPCRow(TruePhotonCandidate->GetConversionRadius())));
   // hTrueElecNfindableClsTPC->Fill((Float_t)elec->GetNcls(1)/(Float_t)elec->GetTPCNclsF());
   // hTruePosiNfindableClsTPC->Fill((Float_t)posi->GetNcls(1)/(Float_t)posi->GetTPCNclsF());
// }

//________________________________________________________________________
void AliAnalysisTaskConversionQA::CountTracks(){

   if(fInputEvent->IsA()==AliESDEvent::Class()){
   // Using standard function for setting Cuts
   Bool_t selectPrimaries=kTRUE;
   AliESDtrackCuts *EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(selectPrimaries);
   EsdTrackCuts->SetMaxDCAToVertexZ(2);
   EsdTrackCuts->SetEtaRange(-0.8, 0.8);
   EsdTrackCuts->SetPtRange(0.15);
      fNumberOfESDTracks = 0;
      for(Int_t iTracks = 0; iTracks < fInputEvent->GetNumberOfTracks(); iTracks++){
         AliESDtrack* curTrack = (AliESDtrack*) fInputEvent->GetTrack(iTracks);
         if(!curTrack) continue;
         // if(fMCEvent && ((AliConversionCuts*)fCutArray->At(fiCut))->GetSignalRejection() != 0){
         //    if(!((AliConversionCuts*)fCutArray->At(fiCut))->IsParticleFromBGEvent(abs(curTrack->GetLabel()), fMCStack)) continue;
         // }
         if(EsdTrackCuts->AcceptTrack(curTrack) ) fNumberOfESDTracks++;
      }
      delete EsdTrackCuts;
      EsdTrackCuts=0x0;
   }
   else if(fInputEvent->IsA()==AliAODEvent::Class()){
      fNumberOfESDTracks = 0;
      for(Int_t iTracks = 0; iTracks<fInputEvent->GetNumberOfTracks(); iTracks++){
         AliAODTrack* curTrack = (AliAODTrack*) fInputEvent->GetTrack(iTracks);
         if(!curTrack->IsPrimaryCandidate()) continue;
         if(abs(curTrack->Eta())>0.8) continue;
         if(curTrack->Pt()<0.15) continue;
         if(abs(curTrack->ZAtDCA())>2) continue;
         fNumberOfESDTracks++;
      }
   }
   
   return;
}
//________________________________________________________________________
/*
Bool_t AliAnalysisTaskConversionQA::IsTruePhotonESD(AliAODConversionPhoton *TruePhotonCandidate)
{
   TParticle *posDaughter = TruePhotonCandidate->GetPositiveMCDaughter(fMCStack);
   TParticle *negDaughter = TruePhotonCandidate->GetNegativeMCDaughter(fMCStack);

   if(posDaughter == NULL || negDaughter == NULL) return kFALSE; // One particle does not exist
   Int_t pdgCode[2] = {abs(posDaughter->GetPdgCode()),abs(negDaughter->GetPdgCode())};
   if(posDaughter->GetMother(0) != negDaughter->GetMother(0)) return kFALSE;
   else if(posDaughter->GetMother(0) == -1) return kFALSE;

   if(pdgCode[0]!=11 || pdgCode[1]!=11) return kFALSE; //One Particle is not electron
   if(posDaughter->GetPdgCode()==negDaughter->GetPdgCode()) return kFALSE; // Same Charge
   if(posDaughter->GetUniqueID() != 5 || negDaughter->GetUniqueID() !=5) return kFALSE;// check if the daughters come from a conversion

   TParticle *Photon = TruePhotonCandidate->GetMCParticle(fMCStack);
   if(Photon->GetPdgCode() != 22) return kFALSE; // Mother is no Photon

   return kTRUE;
}
*/
UInt_t AliAnalysisTaskConversionQA::IsTruePhotonESD(AliAODConversionPhoton *TruePhotonCandidate)
{


	UInt_t kind = 9;
	TParticle *posDaughter = TruePhotonCandidate->GetPositiveMCDaughter(fMCStack);
	TParticle *negDaughter = TruePhotonCandidate->GetNegativeMCDaughter(fMCStack);
	TParticle *Photon = TruePhotonCandidate->GetMCParticle(fMCStack);
	Int_t motherLabelPhoton; 
	Int_t pdgCodePos = 0; 
	Int_t pdgCodeNeg = 0; 
	Int_t pdgCode = 0; 


	if(posDaughter == NULL || negDaughter == NULL) {
		kind = 9;
		//		return kFALSE; // One particle does not exist
   
  } else if( posDaughter->GetMother(0) != negDaughter->GetMother(0)  || (posDaughter->GetMother(0) == negDaughter->GetMother(0) && posDaughter->GetMother(0) ==-1)) {
		//		kind = 1;
	  return 1;
		pdgCodePos=TMath::Abs(posDaughter->GetPdgCode());
		pdgCodeNeg=TMath::Abs(negDaughter->GetPdgCode());
		if(pdgCodePos==11 && pdgCodeNeg==11) return 10; //Electron Combinatorial
		if(pdgCodePos==11 && pdgCodeNeg==11 && 
			(posDaughter->GetMother(0) == negDaughter->GetMother(0) && posDaughter->GetMother(0) ==-1)) return 15; //direct Electron Combinatorial
				
		if(pdgCodePos==211 && pdgCodeNeg==211) return 11; //Pion Combinatorial
		if((pdgCodePos==211 && pdgCodeNeg==2212) ||(pdgCodePos==2212 && pdgCodeNeg==211))	return 12; //Pion, Proton Combinatorics
		if((pdgCodePos==211 && pdgCodeNeg==11) ||(pdgCodePos==11 && pdgCodeNeg==211)) return 13; //Pion, Electron Combinatorics
		if(pdgCodePos==321 || pdgCodeNeg==321) return 14; //Kaon combinatorics
	}else{		
 
		pdgCodePos=posDaughter->GetPdgCode();
		pdgCodeNeg=negDaughter->GetPdgCode();
	  motherLabelPhoton= Photon->GetMother(0);
		if ( TruePhotonCandidate->GetMCParticle(fMCStack)->GetPdgCode()) pdgCode = TruePhotonCandidate->GetMCParticle(fMCStack)->GetPdgCode(); 

		if(TMath::Abs(pdgCodePos)!=11 || TMath::Abs(pdgCodeNeg)!=11) return 2; // true from hadronic decays
		else if ( !(pdgCodeNeg==pdgCodePos)){
			if(pdgCode == 111) return 3; // pi0 Dalitz
			else if (pdgCode == 221) return 4; // eta Dalitz
			else if (!(negDaughter->GetUniqueID() != 5 || posDaughter->GetUniqueID() !=5)){
				if(pdgCode == 22 && motherLabelPhoton < fMCStack->GetNprimary()){
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
UInt_t AliAnalysisTaskConversionQA::IsTruePhotonAOD(AliAODConversionPhoton *TruePhotonCandidate)
{   

 	UInt_t kind = 9;

  TClonesArray *AODMCTrackArray = dynamic_cast<TClonesArray*>(fInputEvent->FindListObject(AliAODMCParticle::StdBranchName()));
	AliAODMCParticle *posDaughter = (AliAODMCParticle*) AODMCTrackArray->At(TruePhotonCandidate->GetMCLabelPositive());
	AliAODMCParticle *negDaughter = (AliAODMCParticle*) AODMCTrackArray->At(TruePhotonCandidate->GetMCLabelNegative());
	AliAODMCParticle *Photon = (AliAODMCParticle*) AODMCTrackArray->At(posDaughter->GetMother());
	Int_t motherLabelPhoton = Photon->GetMother();
	Int_t pdgCodePos = 0; 
	Int_t pdgCodeNeg = 0; 
	Int_t pdgCode = 0; 

	if(posDaughter == NULL || negDaughter == NULL) {
		kind = 9;
		//		return kFALSE; // One particle does not exist
   
   } else if( posDaughter->GetMother() != negDaughter->GetMother()  || (posDaughter->GetMother() == negDaughter->GetMother() && posDaughter->GetMother() ==-1)) {
		kind = 1;

		pdgCodePos=TMath::Abs(posDaughter->GetPdgCode());
		pdgCodeNeg=TMath::Abs(negDaughter->GetPdgCode());
		if(pdgCodePos==11 && pdgCodeNeg==11)	kind = 10; //Electron Combinatorial
		if(pdgCodePos==11 && pdgCodeNeg==11 && 
			(posDaughter->GetMother() == negDaughter->GetMother() && posDaughter->GetMother() ==-1))kind = 15; //direct Electron Combinatorial
				
		if(pdgCodePos==211 && pdgCodeNeg==211) kind = 11; //Pion Combinatorial
		if((pdgCodePos==211 && pdgCodeNeg==2212) ||(pdgCodePos==2212 && pdgCodeNeg==211))	kind = 12; //Pion, Proton Combinatorics
		if((pdgCodePos==211 && pdgCodeNeg==11) ||(pdgCodePos==11 && pdgCodeNeg==211)) kind = 13; //Pion, Electron Combinatorics
		if(pdgCodePos==321 || pdgCodeNeg==321) kind = 14; //Kaon combinatorics
	}else{		
 
		pdgCodePos=posDaughter->GetPdgCode();
		pdgCodeNeg=negDaughter->GetPdgCode();

		if ( Photon->GetPdgCode()) 
			pdgCode = Photon->GetPdgCode(); 
		if(TMath::Abs(pdgCodePos)!=11 || TMath::Abs(pdgCodeNeg)!=11) kind = 2; // true from hadronic decays
		else if ( !(pdgCodeNeg==pdgCodePos)){
			if(pdgCode == 111) kind = 3; // pi0 Dalitz
			else if (pdgCode == 221) kind = 4; // eta Dalitz
			else if (!(negDaughter->GetMCProcessCode() != 5 || posDaughter->GetMCProcessCode() !=5)){
				if(pdgCode == 22 && motherLabelPhoton < fMCStack->GetNprimary()){
					kind = 0; // primary photons
				} else if (pdgCode == 22){
					kind = 5; //secondary photons
				}
			}
		}
	}


	return kind;


}

//________________________________________________________________________
void AliAnalysisTaskConversionQA::RelabelAODPhotonCandidates(Bool_t mode){

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
               PhotonCandidate->SetMCLabelPositive(abs(tempDaughter->GetLabel()));
               AODLabelPos = kTRUE;
            }
         }
         if(!AODLabelNeg){
            if( tempDaughter->GetID() == PhotonCandidate->GetTrackLabelNegative()){
               PhotonCandidate->SetMCLabelNegative(abs(tempDaughter->GetLabel()));
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
                  if( (abs(tempDaughter->GetID())-1) == PhotonCandidate->GetTrackLabelPositive()){
                     PhotonCandidate->SetMCLabelPositive(abs(tempDaughter->GetLabel()));
                     AODLabelPos = kTRUE;
                  }
               }
               if(!AODLabelNeg){
                  if( (abs(tempDaughter->GetID())-1) == PhotonCandidate->GetTrackLabelNegative()){
                     PhotonCandidate->SetMCLabelNegative(abs(tempDaughter->GetLabel()));
                     AODLabelNeg = kTRUE;
                  }
               }
            }
            if(AODLabelNeg && AODLabelPos){
               break;
            }
         }
         if(!AODLabelPos || !AODLabelNeg){
            cout<<"WARNING!!! AOD TRACKS NOT FOUND FOR"<<endl;
         }
      }
   }
   
   if(!mode){
      delete[] fMCStackPos;
      delete[] fMCStackNeg;
   }
}

void AliAnalysisTaskConversionQA::SetLogBinningXTH2(TH2* histoRebin){
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
void AliAnalysisTaskConversionQA::Terminate(Option_t *)
{

}
