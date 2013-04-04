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

class iostream;

using namespace std;

ClassImp(AliAnalysisTaskConversionQA)

//________________________________________________________________________
AliAnalysisTaskConversionQA::AliAnalysisTaskConversionQA(const char *name) : AliAnalysisTaskSE(name),
   fV0Reader(NULL),
   fConversionGammas(NULL),
   fConversionCuts(NULL),
   fInputEvent(NULL),
   fNumberOfESDTracks(0),
   fMCEvent(NULL),
   fMCStack(NULL),
   fStreamQA(NULL),
   fIsHeavyIon(kFALSE),
   ffillTree(kFALSE),
   ffillHistograms(kFALSE),
   fOutputList(NULL),
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
   hGammaQt(NULL),
   hGammaCosinePointingAngle(NULL),
   hGammaXY(NULL),
   hGammaZR(NULL),
   hElecPt(NULL),
   hElecEta(NULL),
   hElecPhi(NULL),
   hElecNfindableClsTPC(NULL),
   hPosiNfindableClsTPC(NULL),
   fTrueList(NULL),
   hTrueResoulutionR(NULL),
   hTrueResoulutionZ(NULL),
   hTrueResoulutionPhi(NULL),
   hTrueGammaPt(NULL),
   hTrueGammaPhi(NULL),
   hTrueGammaEta(NULL),
   hTrueGammaMass(NULL),
   hTrueGammaChi2perNDF(NULL),
   hTrueGammaPsiPair(NULL),
   hTrueGammaQt(NULL),
   hTrueGammaCosinePointingAngle(NULL),
   hTrueGammaXY(NULL),
   hTrueGammaZR(NULL),
   hTrueElecPt(NULL),
   hTrueElecEta(NULL),
   hTrueElecPhi(NULL),
   hTrueElecNfindableClsTPC(NULL),
   hTruePosiNfindableClsTPC(NULL)
{
   // Default constructor

   DefineInput(0, TChain::Class());
   DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskConversionQA::~AliAnalysisTaskConversionQA()
{
   // default deconstructor
   if(fStreamQA){
      delete fStreamQA;
      fStreamQA = 0x0;
   }
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
      TH1::SetDefaultSumw2(kTRUE);
      
      fESDList = new TList();
      fESDList->SetOwner(kTRUE);
      fESDList->SetName("ESD QA");
      fOutputList->Add(fESDList);

      hVertexZ = new TH1F("Vertex_Z","Vertex_Z",300,-15,15);
      fESDList->Add(hVertexZ);
      hNContributorsVertex = new TH1I("ContrVertex_Z","ContrVertex_Z",3000,0,3000);
      fESDList->Add(hNContributorsVertex);
      if(fIsHeavyIon) hNGoodESDTracks = new TH1I("GoodESDTracks","GoodESDTracks",3000,0,3000);
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
      hGammaEta = new TH1F("Gamma_Eta","Gamma_Eta",400,-1.2,1.2);
      fESDList->Add(hGammaEta);
      hGammaChi2perNDF = new TH1F("Gamma_Chi2perNDF","Gamma_Chi2perNDF",500,0,100);
      fESDList->Add(hGammaChi2perNDF);
      hGammaPsiPair = new TH1F("Gamma_PsiPair","Gamma_PsiPair",500,0,2);
      fESDList->Add(hGammaPsiPair);
      hGammaQt = new TH1F("Gamma_Qt","Gamma_Qt",400,0,0.1);
      fESDList->Add(hGammaQt);
      hGammaCosinePointingAngle = new TH1F("Gamma_CosinePointingAngle","Gamma_CosinePointingAngle",900,0.7,1.);
      fESDList->Add(hGammaCosinePointingAngle);
      hGammaXY = new TH2F("Gamma_ConversionPoint_XY","Gamma_ConversionPoint_XY",960,-120,120,960,-120,120);
      fESDList->Add(hGammaXY);
      hGammaZR= new TH2F("Gamma_ConversionPoint_ZR","Gamma_ConversionPoint_ZR",1200,-150,150,480,0,120);
      fESDList->Add(hGammaZR);

      hElecPt = new TH2F("Electron_Positron_Pt","Electron_Positron_Pt",250,0,25,250,0,25);
      fESDList->Add(hElecPt);
      hElecEta = new TH2F("Electron_Positron_Eta","Electron_Positron_Eta",400,-1.2,1.2,400,-1.2,1.2);
      fESDList->Add(hElecEta);
      hElecPhi = new TH2F("Electron_Positron_Phi","Electron_Positron_Phi",360,0,2*TMath::Pi(),360,0,2*TMath::Pi());
      fESDList->Add(hElecPhi);
      hElecNfindableClsTPC = new TH1F("Electron_findableClusterTPC","Electron_findableClusterTPC",100,0,1);
      fESDList->Add(hElecNfindableClsTPC);
      hPosiNfindableClsTPC = new TH1F("Positron_findableClusterTPC","Positron_findableClusterTPC",100,0,1);
      fESDList->Add(hPosiNfindableClsTPC);

      if(MCEvent()){
         fTrueList = new TList();
         fTrueList->SetOwner(kTRUE);
         fTrueList->SetName("True QA");
         fOutputList->Add(fTrueList);

         hTrueResoulutionR = new TH2F("True_ConversionPointResolution_R","True_ConversionPointResolution_R",240,0,120,200,-20,20);
         fTrueList->Add(hTrueResoulutionR);
         hTrueResoulutionZ = new TH2F("True_ConversionPointResolution_Z","True_ConversionPointResolution_Z",480,-120,120,200,-20,20);
         fTrueList->Add(hTrueResoulutionZ);
         hTrueResoulutionPhi = new TH2F("True_ConversionPointResolution_Phi","True_ConversionPointResolution_Phi",360,0,2*TMath::Pi(),200,-20,20);
         fTrueList->Add(hTrueResoulutionPhi);

         hTrueGammaPt = new TH1F("True_Gamma_Pt","True_Gamma_Pt",250,0,25);
         fTrueList->Add(hTrueGammaPt);
         hTrueGammaPhi = new TH1F("True_Gamma_Phi","True_Gamma_Phi",360,0,2*TMath::Pi());
         fTrueList->Add(hTrueGammaPhi);
         hTrueGammaEta = new TH1F("True_Gamma_Eta","True_Gamma_Eta",400,-1.2,1.2);
         fTrueList->Add(hTrueGammaEta);
         hTrueGammaMass = new TH1F("True_Gamma_Mass","True_Gamma_Mass",1000,0,0.3);
         fTrueList->Add(hTrueGammaMass);
         hTrueGammaChi2perNDF = new TH1F("True_Gamma_Chi2perNDF","True_Gamma_Chi2perNDF",500,0,100);
         fTrueList->Add(hTrueGammaChi2perNDF);
         hTrueGammaPsiPair = new TH1F("True_Gamma_PsiPair","True_Gamma_PsiPair",500,0,2);
         fTrueList->Add(hTrueGammaPsiPair);
         hTrueGammaQt = new TH1F("True_Gamma_Qt","True_Gamma_Qt",400,0,0.1);
         fTrueList->Add(hTrueGammaQt);
         hTrueGammaCosinePointingAngle = new TH1F("True_Gamma_CosinePointingAngle","True_Gamma_CosinePointingAngle",900,0.7,1.);
         fTrueList->Add(hTrueGammaCosinePointingAngle);
         hTrueGammaXY = new TH2F("True_Gamma_ConversionPoint_XY","True_Gamma_ConversionPoint_XY",960,-120,120,960,-120,120);
         fTrueList->Add(hTrueGammaXY);
         hTrueGammaZR= new TH2F("TrueGamma_ConversionPoint_ZR","TrueGamma_ConversionPoint_ZR",1200,-150,150,480,0,120);
         fTrueList->Add(hTrueGammaZR);

         hTrueElecPt = new TH2F("True_Electron_Positron_Pt","True_Electron_Positron_Pt",250,0,25,250,0,25);
         fTrueList->Add(hTrueElecPt);
         hTrueElecEta = new TH2F("True_Electron_Positron_Eta","True_Electron_Positron_Eta",400,-1.2,1.2,400,-1.2,1.2);
         fTrueList->Add(hTrueElecEta);
         hTrueElecPhi = new TH2F("True_Electron_Positron_Phi","True_Electron_Positron_Phi",360,0,2*TMath::Pi(),360,0,2*TMath::Pi());
         fTrueList->Add(hTrueElecPhi);
         hTrueElecNfindableClsTPC = new TH1F("True_Electron_findableClusterTPC","True_Electron_findableClusterTPC",100,0,1);
         fTrueList->Add(hTrueElecNfindableClsTPC);
         hTruePosiNfindableClsTPC = new TH1F("True_Positron_findableClusterTPC","True_Positron_findableClusterTPC",100,0,1);
         fTrueList->Add(hTruePosiNfindableClsTPC);
      }
      if(fConversionCuts->GetCutHistograms()){
         fOutputList->Add(fConversionCuts->GetCutHistograms());
      }
      TH1::SetDefaultSumw2(kFALSE);
   }
   
   if(ffillTree){
      TString cutnumber = fConversionCuts->GetCutNumber();
      fStreamQA = new TTreeSRedirector(Form("GammaConvV1_QATree_%s.root",cutnumber.Data()),"recreate");
   }

   PostData(1, fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskConversionQA::UserExec(Option_t *){

   fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1");

   Int_t eventQuality = ((AliConversionCuts*)fV0Reader->GetConversionCuts())->GetEventQuality();
   if(eventQuality != 0){// Event Not Accepted
      return;
   }
   fInputEvent = InputEvent();
   fMCEvent = MCEvent();
   if(fMCEvent) fMCStack = fMCEvent->Stack();

   Int_t eventNotAccepted =
      fConversionCuts->IsEventAcceptedByConversionCut(fV0Reader->GetConversionCuts(),fInputEvent,fMCEvent,fIsHeavyIon);
   if(eventNotAccepted) return; // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1

   fConversionGammas=fV0Reader->GetReconstructedGammas();

   if(fMCEvent){
      if(fConversionCuts->GetSignalRejection() != 0){
         fConversionCuts->GetNotRejectedParticles(fConversionCuts->GetSignalRejection(),
                                                  fConversionCuts->GetAcceptedHeader(),
                                                  fMCEvent);
      }
   }

   if(ffillHistograms){
      CountESDTracks();
      hVertexZ->Fill(fInputEvent->GetPrimaryVertex()->GetZ());
      hNContributorsVertex->Fill(fConversionCuts->GetNumberOfContributorsVtx(fInputEvent));
      hNGoodESDTracks->Fill(fNumberOfESDTracks);
      hNV0Tracks->Fill(fInputEvent->GetVZEROData()->GetMTotV0A()+fInputEvent->GetVZEROData()->GetMTotV0C());
   }

   for(Int_t firstGammaIndex=0;firstGammaIndex<fConversionGammas->GetEntriesFast();firstGammaIndex++){
      AliAODConversionPhoton *gamma=dynamic_cast<AliAODConversionPhoton*>(fConversionGammas->At(firstGammaIndex));
      if (gamma==NULL) continue;
      if(fMCEvent && fConversionCuts->GetSignalRejection() != 0){
         if(!fConversionCuts->IsParticleFromBGEvent(gamma->GetMCLabelPositive(), fMCStack))
            continue;
         if(!fConversionCuts->IsParticleFromBGEvent(gamma->GetMCLabelNegative(), fMCStack))
            continue;
      }
      if(!fConversionCuts->PhotonIsSelected(gamma,fInputEvent)){
         continue;
      }

      if(ffillTree) ProcessQATree(gamma);
      if(ffillHistograms) ProcessQA(gamma);
   }

   PostData(1, fOutputList);
}


///________________________________________________________________________
void AliAnalysisTaskConversionQA::ProcessQATree(AliAODConversionPhoton *gamma){

   // Fill Histograms for QA and MC
   AliESDEvent* event = (AliESDEvent*) InputEvent();
   AliPIDResponse* pidResonse = ((AliConversionCuts*)fV0Reader->GetConversionCuts())->GetPIDResponse();

   Float_t gammaPt = gamma->GetPhotonPt();
   Float_t gammaPhi = gamma->GetPhotonPhi();
   Float_t gammaTheta = gamma->Theta();
   Float_t gammaChi2NDF = gamma->GetChi2perNDF();
   Float_t gammaQt = gamma->GetArmenterosQt();
   Float_t gammaAlpha = gamma->GetArmenterosAlpha();
   Float_t gammaPsiPair = gamma->GetPsiPair();
   Float_t gammaCosPointing = fConversionCuts->GetCosineOfPointingAngle(gamma,event);
   TVectorF conversionPoint(3);
   conversionPoint(0) = gamma->GetConversionX();
   conversionPoint(1) = gamma->GetConversionY();
   conversionPoint(2) = gamma->GetConversionZ();
   TVectorF daughterProp(18);
   AliESDtrack * negTrack = fConversionCuts->GetESDTrack(event, gamma->GetTrackLabelNegative());
   AliESDtrack * posTrack = fConversionCuts->GetESDTrack(event, gamma->GetTrackLabelPositive());

   if(!negTrack||!posTrack)return;

   Bool_t isTruePhoton = kFALSE;
   if(fMCEvent){
      if(IsTruePhoton(gamma)) isTruePhoton = kTRUE;
   }

   daughterProp(0) = posTrack->Pt();
   daughterProp(7) = negTrack->Pt();
   daughterProp(1) = posTrack->Theta();
   daughterProp(8) = negTrack->Theta();
   Double32_t signalPos[4] = {0,0,0,0};
   Char_t nclPos[3];
   Char_t nrowsPos[3];
   if (posTrack->GetTPCdEdxInfo()) {
      posTrack->GetTPCdEdxInfo()->GetTPCSignalRegionInfo(signalPos,nclPos,nrowsPos);
      daughterProp(2) = signalPos[0];
      daughterProp(14) = signalPos[1];
      daughterProp(16) = signalPos[2];
   } else {
      daughterProp(2) = posTrack->GetTPCsignal();
      daughterProp(14) = 0;
      daughterProp(16) = 0;
   }
   Double32_t signalNeg[4] = {0,0,0,0};
   Char_t nclNeg[3];
   Char_t nrowsNeg[3];
   if (negTrack->GetTPCdEdxInfo()) {
      negTrack->GetTPCdEdxInfo()->GetTPCSignalRegionInfo(signalNeg,nclNeg,nrowsNeg);
      daughterProp(9) = signalNeg[0];
      daughterProp(15) = signalNeg[1];
      daughterProp(17) = signalNeg[2];
   } else {
      daughterProp(9) = negTrack->GetTPCsignal();
      daughterProp(15) = 0;
      daughterProp(17) = 0;
   }
   daughterProp(3) = pidResonse->NumberOfSigmasTPC(posTrack,AliPID::kElectron);
   daughterProp(10) = pidResonse->NumberOfSigmasTPC(negTrack,AliPID::kElectron);
   if((posTrack->GetStatus() & AliESDtrack::kTOFpid) && !(posTrack->GetStatus() & AliESDtrack::kTOFmismatch)){
      Double_t t0pos = pidResonse->GetTOFResponse().GetStartTime(posTrack->P());
      Double_t timesPos[5];
      posTrack->GetIntegratedTimes(timesPos);
      Double_t TOFsignalPos =	posTrack->GetTOFsignal();
      Double_t dTpos = TOFsignalPos - t0pos - timesPos[0];
      daughterProp(4) = dTpos;
      daughterProp(5) = pidResonse->NumberOfSigmasTOF(posTrack, AliPID::kElectron);
   } else {
      daughterProp(4) = 20000;
      daughterProp(5) = -20;
   }
   if((negTrack->GetStatus() & AliESDtrack::kTOFpid) && !(negTrack->GetStatus() & AliESDtrack::kTOFmismatch)){
      Double_t t0neg = pidResonse->GetTOFResponse().GetStartTime(negTrack->P());
      Double_t timesNeg[5];
      negTrack->GetIntegratedTimes(timesNeg);
      Double_t TOFsignalNeg =	negTrack->GetTOFsignal();
      Double_t dTneg = TOFsignalNeg - t0neg - timesNeg[0];
      daughterProp(11) = dTneg;
      daughterProp(12) = pidResonse->NumberOfSigmasTOF(negTrack, AliPID::kElectron);
   } else {
      daughterProp(11) = 20000;
      daughterProp(12) = -20;
   }

   daughterProp(6) = (Float_t)posTrack->GetTPCClusterInfo(2,0,fConversionCuts->GetFirstTPCRow(gamma->GetConversionRadius()));
   //posTrack->GetNcls(1)/(Float_t)posTrack->GetTPCNclsF();
   daughterProp(13) = (Float_t)negTrack->GetTPCClusterInfo(2,0,fConversionCuts->GetFirstTPCRow(gamma->GetConversionRadius()));
   //negTrack->GetNcls(1)/(Float_t)negTrack->GetTPCNclsF();

   if (fStreamQA){
      if(fMCEvent){
         (*fStreamQA)<<"PhotonQA"
                     << "pt=" << gammaPt
                     << "phi=" << gammaPhi
                     << "theta=" << gammaTheta
                     << "chi2ndf=" << gammaChi2NDF
                     << "qt="<< gammaQt
                     << "alpha=" << gammaAlpha
                     << "psipair=" << gammaPsiPair
                     << "cosPoint=" << gammaCosPointing
                     << "TruePhoton=" << isTruePhoton
                     << "conversionPoint=" << &conversionPoint
                     << "daugtherProp.=" << &daughterProp
                     << "\n";
      }
      else{
         (*fStreamQA)<<"PhotonQA"
                     << "pt=" << gammaPt
                     << "phi=" << gammaPhi
                     << "theta=" << gammaTheta
                     << "chi2ndf=" << gammaChi2NDF
                     << "qt="<< gammaQt
                     << "alpha=" << gammaAlpha
                     << "psipair=" << gammaPsiPair
                     << "cosPoint=" << gammaCosPointing
                     << "conversionPoint=" << &conversionPoint
                     << "daugtherProp.=" << &daughterProp
                     << "\n";
      }
   }
}

//_____________________________________________________________________________________________________
void AliAnalysisTaskConversionQA::ProcessQA(AliAODConversionPhoton *gamma){


   // Fill Histograms for QA and MC
   hGammaPt->Fill(gamma->GetPhotonPt());
   hGammaPhi->Fill(gamma->GetPhotonPhi());
   hGammaEta->Fill(gamma->Eta());
   hGammaChi2perNDF->Fill(gamma->GetChi2perNDF());
   hGammaPsiPair->Fill(gamma->GetPsiPair());
   hGammaQt->Fill(gamma->GetArmenterosQt());
   hGammaCosinePointingAngle->Fill(fConversionCuts->GetCosineOfPointingAngle(gamma,fInputEvent));
   hGammaXY->Fill(gamma->GetConversionX(),gamma->GetConversionY());
   hGammaZR->Fill(gamma->GetConversionZ(),gamma->GetConversionRadius());

   AliESDtrack * negTrack = fConversionCuts->GetESDTrack((AliESDEvent*)fInputEvent, gamma->GetTrackLabelNegative());
   AliESDtrack * posTrack = fConversionCuts->GetESDTrack((AliESDEvent*)fInputEvent, gamma->GetTrackLabelPositive());
   if(!negTrack||!posTrack)return;

   Double_t negtrackPhi = negTrack->Phi();
   Double_t postrackPhi = posTrack->Phi();
   hITSClusterPhi->Fill(negtrackPhi,6);
   hITSClusterPhi->Fill(postrackPhi,6);
   for(Int_t itsLayer = 0; itsLayer<6;itsLayer++){
      if(TESTBIT(negTrack->GetITSClusterMap(),itsLayer)){
         hITSClusterPhi->Fill(negtrackPhi,itsLayer);
      }
      if(TESTBIT(posTrack->GetITSClusterMap(),itsLayer)){
         hITSClusterPhi->Fill(postrackPhi,itsLayer);
      }
   }

   hElecPt->Fill(negTrack->Pt(),posTrack->Pt());
   hElecEta->Fill(negTrack->Eta(),posTrack->Eta());
   hElecPhi->Fill(negTrack->Phi(),posTrack->Phi());

   hElecNfindableClsTPC->Fill((Float_t)posTrack->GetTPCClusterInfo(2,0,fConversionCuts->GetFirstTPCRow(gamma->GetConversionRadius())));
   hPosiNfindableClsTPC->Fill((Float_t)negTrack->GetTPCClusterInfo(2,0,fConversionCuts->GetFirstTPCRow(gamma->GetConversionRadius())));

   // hElecNfindableClsTPC->Fill((Float_t)negTrack->GetNcls(1)/(Float_t)negTrack->GetTPCNclsF());
   // hPosiNfindableClsTPC->Fill((Float_t)posTrack->GetNcls(1)/(Float_t)posTrack->GetTPCNclsF());

   if(fMCEvent) ProcessTrueQA(gamma,negTrack,posTrack);
}

//________________________________________________________________________
void AliAnalysisTaskConversionQA::ProcessTrueQA(AliAODConversionPhoton *TruePhotonCandidate, AliESDtrack *elec, AliESDtrack *posi)
{

   if(!IsTruePhoton(TruePhotonCandidate)) return;

   TParticle *negDaughter = TruePhotonCandidate->GetNegativeMCDaughter(fMCStack);
   TParticle *mcPhoton = TruePhotonCandidate->GetMCParticle(fMCStack);

   // True Photon
   hTrueResoulutionR->Fill(TruePhotonCandidate->GetConversionRadius(),
                           TruePhotonCandidate->GetConversionRadius()-negDaughter->R());
   hTrueResoulutionZ->Fill(TruePhotonCandidate->GetConversionZ(),
                           TruePhotonCandidate->GetConversionZ()-negDaughter->Vz());
   hTrueResoulutionPhi->Fill(TruePhotonCandidate->Phi(),
                             TruePhotonCandidate->Phi()-mcPhoton->Phi());
   hTrueGammaPt->Fill(TruePhotonCandidate->Pt());
   hTrueGammaPhi->Fill(TruePhotonCandidate->Phi());
   hTrueGammaEta->Fill(TruePhotonCandidate->Eta());
   hTrueGammaMass->Fill(TruePhotonCandidate->GetMass());
   hTrueGammaChi2perNDF->Fill(TruePhotonCandidate->GetChi2perNDF());
   hTrueGammaPsiPair->Fill(TruePhotonCandidate->GetPsiPair());
   hTrueGammaQt->Fill(TruePhotonCandidate->GetArmenterosQt());
   hTrueGammaCosinePointingAngle->Fill(fConversionCuts->GetCosineOfPointingAngle(TruePhotonCandidate,fInputEvent));
   hTrueGammaXY->Fill(TruePhotonCandidate->GetConversionX(),TruePhotonCandidate->GetConversionY());
   hTrueGammaZR->Fill(TruePhotonCandidate->GetConversionZ(),TruePhotonCandidate->GetConversionRadius());

   hTrueElecPt->Fill(elec->Pt(),posi->Pt());
   hTrueElecEta->Fill(elec->Eta(),posi->Eta());
   hTrueElecPhi->Fill(elec->Phi(),posi->Phi());
   hTrueElecNfindableClsTPC
      ->Fill((Float_t)elec->GetTPCClusterInfo(2,0,fConversionCuts->GetFirstTPCRow(TruePhotonCandidate->GetConversionRadius())));
   hTruePosiNfindableClsTPC
      ->Fill((Float_t)posi->GetTPCClusterInfo(2,0,fConversionCuts->GetFirstTPCRow(TruePhotonCandidate->GetConversionRadius())));
   // hTrueElecNfindableClsTPC->Fill((Float_t)elec->GetNcls(1)/(Float_t)elec->GetTPCNclsF());
   // hTruePosiNfindableClsTPC->Fill((Float_t)posi->GetNcls(1)/(Float_t)posi->GetTPCNclsF());

}
//________________________________________________________________________
void AliAnalysisTaskConversionQA::CountESDTracks(){

   // Using standard function for setting Cuts
   Bool_t selectPrimaries=kTRUE;
   AliESDtrackCuts *EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(selectPrimaries);
   EsdTrackCuts->SetMaxDCAToVertexZ(2);
   EsdTrackCuts->SetEtaRange(-0.8, 0.8);
   EsdTrackCuts->SetPtRange(0.15);

   fNumberOfESDTracks = 0;
   for(Int_t iTracks = 0; iTracks < fInputEvent->GetNumberOfTracks(); iTracks++){
      AliESDtrack* curTrack = (AliESDtrack*)fInputEvent->GetTrack(iTracks);
      if(!curTrack) continue;
      // if(fMCEvent && fConversionCuts->GetSignalRejection() != 0){
      //    if(!fConversionCuts->IsParticleFromBGEvent(abs(curTrack->GetLabel()), fMCStack)) continue;
      // }
      if(EsdTrackCuts->AcceptTrack(curTrack) ) fNumberOfESDTracks++;
   }
   delete EsdTrackCuts;
   EsdTrackCuts=0x0;

   return;
}
//________________________________________________________________________
Bool_t AliAnalysisTaskConversionQA::IsTruePhoton(AliAODConversionPhoton *TruePhotonCandidate)
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

//________________________________________________________________________
void AliAnalysisTaskConversionQA::Terminate(Option_t *)
{

}
