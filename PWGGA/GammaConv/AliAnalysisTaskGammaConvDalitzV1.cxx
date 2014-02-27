/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Pedro Gonz??lez, Pedro Ladr??n de Guevara, Ernesto L??pez Torres, *
 *         Eulogio Serradilla, Ana Marin, Friederike Bock                 *
 * Version 2                                                              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// Analysis task for pi0->e+e-gamma (Dalitz decay)
// Analysis task for chic->JPsi+gamma

#include <vector>

#include "TParticle.h"
#include "TPDGCode.h"
#include "TMCProcess.h"
#include "TDatabasePDG.h"
#include "TList.h"
#include "TChain.h"
#include "TDirectory.h"
#include "TTree.h"
#include "TH1.h"
#include "TH1F.h"
#include "THnSparse.h"
#include "TH2F.h"
#include "AliStack.h"
#include "AliAnalysisManager.h"
#include "AliESDInputHandler.h"
#include "AliESDtrack.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliMCEventHandler.h"
#include "AliPID.h"
#include "AliLog.h"
#include "AliESDtrackCuts.h"
#include "AliESDpidCuts.h"
#include "AliMCEvent.h"
#include "AliESDv0.h"
#include "AliESDEvent.h"
#include "AliESDpid.h"
#include "AliKFParticle.h"
#include "AliMCEventHandler.h"
#include "AliKFVertex.h"
#include "AliTriggerAnalysis.h"
#include "AliCentrality.h"
#include "AliMultiplicity.h"
#include "AliAnalysisTaskGammaConvDalitzV1.h"


ClassImp( AliAnalysisTaskGammaConvDalitzV1 )

//-----------------------------------------------------------------------------------------------
AliAnalysisTaskGammaConvDalitzV1::AliAnalysisTaskGammaConvDalitzV1():
fV0Reader(NULL),
   fElecSelector(NULL),
   fBGHandler(NULL),
   fESDEvent(NULL),
   fMCEvent(NULL),
   fMCStack(NULL),
   fCutFolder(NULL),
   fESDList(NULL),
   fBackList(NULL),
   fMotherList(NULL),
   fTrueList(NULL),
   fMCList(NULL),
   fQAFolder(NULL),
   fOutputContainer(0),
   fReaderGammas(NULL),
   fSelectorElectronIndex(0),
   fSelectorPositronIndex(0),
   fGoodGammas(NULL),
   fGoodVirtualGammas(NULL),
   fGoodElectrons(NULL),
   fGoodPositrons(NULL),
   fCutGammaArray(NULL),
   fCutElectronArray(NULL),
   fCutMesonArray(NULL),
   fGammasPool(NULL),
   fConversionCuts(NULL),
   hESDConvGammaPt(NULL),
   hESDConvGammaEta(NULL),
   hESDConvGammaZR(NULL),
   hESDDalitzElectronPt(NULL),
   hESDDalitzPositronPt(NULL),
   hESDDalitzElectronPhi(NULL),
   hESDDalitzPositronPhi(NULL),
   hESDDalitzElectronAfterPt(NULL),
   hESDDalitzPositronAfterPt(NULL),
   hESDDalitzElectronAfterEta(NULL),
   hESDDalitzPositronAfterEta(NULL),
   hESDDalitzElectronAfterPhi(NULL),
   hESDDalitzPositronAfterPhi(NULL),
   hESDDalitzElectronAfterNClsITS(NULL),
   hESDDalitzPositronAfterNClsITS(NULL),
   hESDDalitzElectronAfterNFindClsTPC(NULL),
   hESDDalitzPositronAfterNFindClsTPC(NULL),
   hESDDalitzElectronAfterNClsTPC(NULL),
   hESDDalitzPositronAfterNClsTPC(NULL),
   hESDDalitzPosEleAfterDCAxy(NULL),
   hESDDalitzPosEleAfterDCAz(NULL),
   hESDDalitzElectronAfterTPCdEdxVsP(NULL),
   hESDDalitzPositronAfterTPCdEdxVsP(NULL),
   hESDDalitzElectronAfterTPCdEdxSignalVsP(NULL),
   hESDDalitzPositronAfterTPCdEdxSignalVsP(NULL),
   hESDDalitzElectronAfterTPCdEdxVsEta(NULL),
   hESDDalitzPositronAfterTPCdEdxVsEta(NULL),
   hESDDalitzElectronAfterTPCdEdxVsPhi(NULL),
   hESDDalitzPositronAfterTPCdEdxVsPhi(NULL),
   hESDMotherPhi(NULL),
   hESDEposEnegPsiPairDPhi(NULL),
   hESDEposEnegInvMassPt(NULL),
   hESDEposEnegLikeSignBackInvMassPt(NULL),
   hESDMotherInvMassPt(NULL),
   hESDPi0MotherInvMassPt(NULL),
   hESDPi0MotherDiffInvMassPt(NULL),
   hESDPi0MotherDiffLimInvMassPt(NULL),
   sESDMotherInvMassPtZM(NULL),
   hESDMotherBackInvMassPt(NULL),
   sESDMotherBackInvMassPtZM(NULL),
   hMCAllGammaPt(NULL),
   hMCConvGammaPt(NULL),
   hMCConvGammaRSPt(NULL),
   hMCAllPositronsPt(NULL),
   hMCAllElectronsPt(NULL),
   hMCPi0DalitzGammaPt(NULL),
   hMCPi0DalitzElectronPt(NULL),
   hMCPi0DalitzPositronPt(NULL),
   hMCPi0Pt(NULL),
   hMCPi0GGPt(NULL),
   hMCEtaPt(NULL),
   hMCEtaGGPt(NULL), 
   hMCPi0InAccPt(NULL),
   hMCEtaInAccPt(NULL),
   hMCChiCPt(NULL),
   hMCChiCInAccPt(NULL),
   hESDEposEnegTruePi0DalitzInvMassPt(NULL),
   hESDEposEnegTruePi0DalitzPsiPairDPhi(NULL),
   hESDEposEnegTrueEtaDalitzInvMassPt(NULL),
   hESDEposEnegTrueEtaDalitzPsiPairDPhi(NULL),
   hESDEposEnegTruePhotonInvMassPt(NULL),
   hESDEposEnegTruePhotonPsiPairDPhi(NULL),
   hESDEposEnegTrueJPsiInvMassPt(NULL),
   hESDTrueMotherChiCInvMassPt(NULL),
   hESDTrueMotherChiCDiffInvMassPt(NULL),
   hESDTrueMotherInvMassPt(NULL),
   hESDTrueMotherDalitzInvMassPt(NULL),
   hESDTrueMotherPi0GGInvMassPt(NULL),
   hESDTruePrimaryMotherPi0GGInvMassPt(NULL),
   hESDTrueSecondaryMotherPi0GGInvMassPt(NULL),
   hESDTruePrimaryMotherInvMassMCPt(NULL),
   hESDTruePrimaryMotherInvMassPt(NULL),
   hESDTruePrimaryMotherW0WeightingInvMassPt(NULL),
   hESDTruePrimaryPi0DalitzESDPtMCPt(NULL),
   hESDTrueSecondaryMotherInvMassPt(NULL),
   hESDTrueSecondaryMotherFromK0sInvMassPt(NULL),
   hESDTrueBckGGInvMassPt(NULL),
   hESDTrueBckContInvMassPt(NULL),
   hESDTrueMotherGGInvMassPt(NULL),
   hESDTrueConvGammaPt(NULL),
   hESDTruePositronPt(NULL),
   hESDTrueElectronPt(NULL),
   hESDTrueSecConvGammaPt(NULL),
   hESDTrueSecPositronPt(NULL),
   hESDTrueSecElectronPt(NULL),
   hESDTruePi0DalitzConvGammaPt(NULL),
   hESDTruePi0DalitzPositronPt(NULL),
   hESDTruePi0DalitzElectronPt(NULL),
   hESDTruePi0DalitzSecConvGammaPt(NULL),
   hESDTruePi0DalitzSecPositronPt(NULL),
   hESDTruePi0DalitzSecElectronPt(NULL),
   hNEvents(NULL),
   hNGoodESDTracks(NULL),
   hEtaShift(NULL),
   fRandom(0),
   fUnsmearedPx(NULL),
   fUnsmearedPy(NULL),
   fUnsmearedPz(NULL),
   fUnsmearedE(NULL),
   fnCuts(0),
   fiCut(0),
   fNumberOfESDTracks(0),
   fMoveParticleAccordingToVertex(kFALSE),
   fIsHeavyIon(kFALSE),
   fDoMesonAnalysis(kTRUE),
   fDoChicAnalysis(kFALSE),
   fDoMesonQA(kFALSE),
   fIsFromMBHeader(kTRUE),
   fIsMC(kFALSE)
{

}

//-----------------------------------------------------------------------------------------------
AliAnalysisTaskGammaConvDalitzV1::AliAnalysisTaskGammaConvDalitzV1( const char* name ):
   AliAnalysisTaskSE(name),
   fV0Reader(NULL),
   fElecSelector(NULL),
   fBGHandler(NULL),
   fESDEvent(NULL),
   fMCEvent(NULL),
   fMCStack(NULL),
   fCutFolder(NULL),
   fESDList(NULL),
   fBackList(NULL),
   fMotherList(NULL),
   fTrueList(NULL),
   fMCList(NULL),
   fQAFolder(NULL),
   fOutputContainer(0),
   fReaderGammas(NULL),
   fSelectorElectronIndex(0),
   fSelectorPositronIndex(0),
   fGoodGammas(NULL),
   fGoodVirtualGammas(NULL),
   fGoodElectrons(NULL),
   fGoodPositrons(NULL),
   fCutGammaArray(NULL),
   fCutElectronArray(NULL),
   fCutMesonArray(NULL),
   fGammasPool(NULL),
   fConversionCuts(NULL),
   hESDConvGammaPt(NULL),
   hESDConvGammaEta(NULL),
   hESDConvGammaZR(NULL),
   hESDDalitzElectronPt(NULL),
   hESDDalitzPositronPt(NULL),
   hESDDalitzElectronPhi(NULL),
   hESDDalitzPositronPhi(NULL),
   hESDDalitzElectronAfterPt(NULL),
   hESDDalitzPositronAfterPt(NULL),
   hESDDalitzElectronAfterEta(NULL),
   hESDDalitzPositronAfterEta(NULL),
   hESDDalitzElectronAfterPhi(NULL),
   hESDDalitzPositronAfterPhi(NULL),
   hESDDalitzElectronAfterNClsITS(NULL),
   hESDDalitzPositronAfterNClsITS(NULL),
   hESDDalitzElectronAfterNFindClsTPC(NULL),
   hESDDalitzPositronAfterNFindClsTPC(NULL),
   hESDDalitzElectronAfterNClsTPC(NULL),
   hESDDalitzPositronAfterNClsTPC(NULL),
   hESDDalitzPosEleAfterDCAxy(NULL),
   hESDDalitzPosEleAfterDCAz(NULL),
   hESDDalitzElectronAfterTPCdEdxVsP(NULL),
   hESDDalitzPositronAfterTPCdEdxVsP(NULL),
   hESDDalitzElectronAfterTPCdEdxSignalVsP(NULL),
   hESDDalitzPositronAfterTPCdEdxSignalVsP(NULL),
   hESDDalitzElectronAfterTPCdEdxVsEta(NULL),
   hESDDalitzPositronAfterTPCdEdxVsEta(NULL),
   hESDDalitzElectronAfterTPCdEdxVsPhi(NULL),
   hESDDalitzPositronAfterTPCdEdxVsPhi(NULL),
   hESDMotherPhi(NULL),
   hESDEposEnegPsiPairDPhi(NULL),
   hESDEposEnegInvMassPt(NULL),
   hESDEposEnegLikeSignBackInvMassPt(NULL),
   hESDMotherInvMassPt(NULL),
   hESDPi0MotherInvMassPt(NULL),
   hESDPi0MotherDiffInvMassPt(NULL),
   hESDPi0MotherDiffLimInvMassPt(NULL),
   sESDMotherInvMassPtZM(NULL),
   hESDMotherBackInvMassPt(NULL),
   sESDMotherBackInvMassPtZM(NULL),
   hMCAllGammaPt(NULL),
   hMCConvGammaPt(NULL),
   hMCConvGammaRSPt(NULL),
   hMCAllPositronsPt(NULL),
   hMCAllElectronsPt(NULL),
   hMCPi0DalitzGammaPt(NULL),
   hMCPi0DalitzElectronPt(NULL),
   hMCPi0DalitzPositronPt(NULL),
   hMCPi0Pt(NULL),
   hMCPi0GGPt(NULL),
   hMCEtaPt(NULL),
   hMCEtaGGPt(NULL),
   hMCPi0InAccPt(NULL),
   hMCEtaInAccPt(NULL),
   hMCChiCPt(NULL),
   hMCChiCInAccPt(NULL),
   hESDEposEnegTruePi0DalitzInvMassPt(NULL),
   hESDEposEnegTruePi0DalitzPsiPairDPhi(NULL),
   hESDEposEnegTrueEtaDalitzInvMassPt(NULL),
   hESDEposEnegTrueEtaDalitzPsiPairDPhi(NULL),
   hESDEposEnegTruePhotonInvMassPt(NULL),
   hESDEposEnegTruePhotonPsiPairDPhi(NULL),
   hESDEposEnegTrueJPsiInvMassPt(NULL),
   hESDTrueMotherChiCInvMassPt(NULL),
   hESDTrueMotherChiCDiffInvMassPt(NULL),
   hESDTrueMotherInvMassPt(NULL),
   hESDTrueMotherDalitzInvMassPt(NULL),
   hESDTrueMotherPi0GGInvMassPt(NULL),
   hESDTruePrimaryMotherPi0GGInvMassPt(NULL),
   hESDTrueSecondaryMotherPi0GGInvMassPt(NULL),
   hESDTruePrimaryMotherInvMassMCPt(NULL),
   hESDTruePrimaryMotherInvMassPt(NULL),
   hESDTruePrimaryMotherW0WeightingInvMassPt(NULL),
   hESDTruePrimaryPi0DalitzESDPtMCPt(NULL),
   hESDTrueSecondaryMotherInvMassPt(NULL),
   hESDTrueSecondaryMotherFromK0sInvMassPt(NULL),
   hESDTrueBckGGInvMassPt(NULL),
   hESDTrueBckContInvMassPt(NULL),
   hESDTrueMotherGGInvMassPt(NULL),
   hESDTrueConvGammaPt(NULL),
   hESDTruePositronPt(NULL),
   hESDTrueElectronPt(NULL),
   hESDTrueSecConvGammaPt(NULL),
   hESDTrueSecPositronPt(NULL),
   hESDTrueSecElectronPt(NULL),
   hESDTruePi0DalitzConvGammaPt(NULL),
   hESDTruePi0DalitzPositronPt(NULL),
   hESDTruePi0DalitzElectronPt(NULL),
   hESDTruePi0DalitzSecConvGammaPt(NULL),
   hESDTruePi0DalitzSecPositronPt(NULL),
   hESDTruePi0DalitzSecElectronPt(NULL),
   hNEvents(NULL),
   hNGoodESDTracks(NULL),
   hEtaShift(NULL),
   fRandom(0),
   fUnsmearedPx(NULL),
   fUnsmearedPy(NULL),
   fUnsmearedPz(NULL),
   fUnsmearedE(NULL),
   fnCuts(0),
   fiCut(0),
   fNumberOfESDTracks(0),
   fMoveParticleAccordingToVertex(kFALSE),
   fIsHeavyIon(kFALSE),
   fDoMesonAnalysis(kTRUE),
   fDoChicAnalysis(kFALSE),
   fDoMesonQA(kFALSE),
   fIsFromMBHeader(kTRUE),
   fIsMC(kFALSE)
{
   DefineOutput(1, TList::Class());
}

//-----------------------------------------------------------------------------------------------
AliAnalysisTaskGammaConvDalitzV1::~AliAnalysisTaskGammaConvDalitzV1()
{
   //
   // virtual destructor
   //
   cout<<"Destructor"<<endl;

   if(fGoodGammas){
      delete fGoodGammas;
      fGoodGammas = 0x0;
   }
   if(fGoodVirtualGammas){
      delete fGoodVirtualGammas;
      fGoodGammas = 0x0;
   }
   if(fGoodElectrons){
      delete fGoodGammas;
      fGoodGammas = 0x0;
   }
   if(fGoodPositrons){
      delete fGoodGammas;
      fGoodGammas = 0x0;
   }
   if(fBGHandler){
      delete[] fBGHandler;
      fBGHandler = 0x0;
   }
   if( fGammasPool ){
      delete[] fGammasPool;
      fGammasPool = 0x0;
   }
}
//___________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::InitBack(){

   const Int_t nDim = 4;
   Int_t nBins[nDim] = {800,250,7,4};
   Double_t xMin[nDim] = {0,0, 0,0};
   Double_t xMax[nDim] = {0.8,25,7,4};
   
   sESDMotherInvMassPtZM = new THnSparseF*[fnCuts];
   sESDMotherBackInvMassPtZM = new THnSparseF*[fnCuts];

   fBGHandler = new AliGammaConversionAODBGHandler*[fnCuts];
   //fBGHandlerRP = new AliConversionAODBGHandlerRP*[fnCuts];
   for(Int_t iCut = 0; iCut<fnCuts;iCut++){
      //if (((AliConversionMesonCuts*)fCutMesonArray->At(iCut))->DoBGCalculation()){

         
         TString cutstringElectron     =   ((AliDalitzElectronCuts*)fCutElectronArray->At(iCut))->GetCutNumber();
         TString cutstringMeson        =   ((AliConversionMesonCuts*)fCutMesonArray->At(iCut))->GetCutNumber();
         TString cutstringGamma        =   ((AliConversionCuts*)fCutGammaArray->At(iCut))->GetCutNumber();


         
         Int_t collisionSystem = atoi((TString)(((AliConversionCuts*)fCutGammaArray->At(iCut))->GetCutNumber())(0,1));
         Int_t centMin = atoi((TString)(((AliConversionCuts*)fCutGammaArray->At(iCut))->GetCutNumber())(1,1));
         Int_t centMax = atoi((TString)(((AliConversionCuts*)fCutGammaArray->At(iCut))->GetCutNumber())(2,1));
         
         if(collisionSystem == 1 || collisionSystem == 2 ||
            collisionSystem == 5 || collisionSystem == 8 ||
            collisionSystem == 9){
            centMin = centMin*10;
            centMax = centMax*10; 
         }
         else if(collisionSystem == 3 || collisionSystem == 6){
            centMin = centMin*5;
            centMax = centMax*5;
         }
         else if(collisionSystem == 4 || collisionSystem == 7){
            centMin = ((centMin*5)+45);
            centMax = ((centMax*5)+45);
         }


         fBackList[iCut] = new TList();
         fBackList[iCut]->SetName(Form("%s_%s_%s Back histograms",cutstringGamma.Data(),cutstringElectron.Data(),cutstringMeson.Data()));
         fBackList[iCut]->SetOwner(kTRUE);
         fCutFolder[iCut]->Add(fBackList[iCut]);

         sESDMotherBackInvMassPtZM[iCut] = new THnSparseF("Back_Back_InvMass_Pt_z_m","Back_Back_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
         fBackList[iCut]->Add(sESDMotherBackInvMassPtZM[iCut]);

         fMotherList[iCut] = new TList();
         fMotherList[iCut]->SetName(Form("%s_%s_%s Mother histograms",cutstringGamma.Data(),cutstringElectron.Data(),cutstringMeson.Data()));
         fMotherList[iCut]->SetOwner(kTRUE);
         fCutFolder[iCut]->Add(fMotherList[iCut]);

         sESDMotherInvMassPtZM[iCut] = new THnSparseF("Back_Mother_InvMass_Pt_z_m","Back_Mother_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
         fMotherList[iCut]->Add(sESDMotherInvMassPtZM[iCut]);

         
         fBGHandler[iCut] = new AliGammaConversionAODBGHandler(
                                                                  collisionSystem,centMin,centMax,
                                                                  ((AliDalitzElectronCuts*)fCutElectronArray->At(iCut))->NumberOfRotationEvents(),
                                                                  ((AliDalitzElectronCuts*)fCutElectronArray->At(iCut))->UseTrackMultiplicity());
        
        if( ( (AliDalitzElectronCuts*)fCutElectronArray->At(iCut))->GetBKGMethod() == 3 ){
         fGammasPool[iCut] = new TList();
        }

      //}
   }
}

//______________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::UserCreateOutputObjects()
{
   //
   // Create ouput objects
   //

   // Create the output container
   if(fOutputContainer != NULL){
      delete fOutputContainer;
      fOutputContainer = NULL;
   }
   if(fOutputContainer == NULL){
      fOutputContainer = new TList();
      fOutputContainer->SetOwner(kTRUE);
   }

   fGoodGammas = new TList();
   //fGoodGammas->SetOwner(kTRUE);


   fGoodVirtualGammas = new TList();
   //fGoodVirtualGammas->SetOwner(kTRUE);



   fGammasPool                     = new TList*[fnCuts];
   fCutFolder                      = new TList*[fnCuts];
   fESDList                        = new TList*[fnCuts];
   fBackList                       = new TList*[fnCuts];
   fMotherList                     = new TList*[fnCuts];
   //fQAFolder                       = new TList*[fnCuts];
   hNEvents                        = new TH1I*[fnCuts];
   hNGoodESDTracks                 = new TH1I*[fnCuts];
   hEtaShift                       = new TProfile*[fnCuts];
   hESDConvGammaPt                 = new TH1F*[fnCuts];
   hESDConvGammaEta 		   = new TH1F*[fnCuts];
   
   hESDDalitzElectronPt            = new TH1F*[fnCuts];
   hESDDalitzPositronPt            = new TH1F*[fnCuts];
   hESDDalitzElectronPhi	   = new TH1F*[fnCuts];
   hESDDalitzPositronPhi	   = new TH1F*[fnCuts];
   
   if( fDoMesonQA ) {
     
   fQAFolder  = new TList*[fnCuts];  
     
   hESDDalitzElectronAfterPt	   = new TH1F*[fnCuts];
   hESDDalitzPositronAfterPt       = new TH1F*[fnCuts];
   hESDDalitzElectronAfterEta      = new TH1F*[fnCuts];
   hESDDalitzPositronAfterEta      = new TH1F*[fnCuts];
   hESDDalitzElectronAfterPhi      = new TH1F*[fnCuts];
   hESDDalitzPositronAfterPhi      = new TH1F*[fnCuts];
   hESDDalitzElectronAfterNClsITS  = new TH1F*[fnCuts];
   hESDDalitzPositronAfterNClsITS  = new TH1F*[fnCuts];
   hESDDalitzElectronAfterNFindClsTPC = new TH2F*[fnCuts];
   hESDDalitzPositronAfterNFindClsTPC = new TH2F*[fnCuts];
   hESDDalitzElectronAfterNClsTPC     = new TH2F*[fnCuts];
   hESDDalitzPositronAfterNClsTPC     = new TH2F*[fnCuts];
   hESDDalitzPosEleAfterDCAxy	   = new TH2F*[fnCuts];
   hESDDalitzPosEleAfterDCAz       = new TH2F*[fnCuts];
   hESDDalitzElectronAfterTPCdEdxVsP    = new TH2F*[fnCuts];
   hESDDalitzPositronAfterTPCdEdxVsP = new TH2F*[fnCuts];
   hESDDalitzElectronAfterTPCdEdxSignalVsP = new TH2F*[fnCuts];
   hESDDalitzPositronAfterTPCdEdxSignalVsP = new TH2F*[fnCuts];
   hESDDalitzElectronAfterTPCdEdxVsEta = new TH2F*[fnCuts];
   hESDDalitzPositronAfterTPCdEdxVsEta = new TH2F*[fnCuts];
   hESDDalitzElectronAfterTPCdEdxVsPhi = new TH2F*[fnCuts];
   hESDDalitzPositronAfterTPCdEdxVsPhi = new TH2F*[fnCuts];
   hESDMotherPhi                   = new TH1F*[fnCuts];
   hESDEposEnegPsiPairDPhi         = new TH2F*[fnCuts];
   hESDEposEnegInvMassPt           = new TH2F*[fnCuts];
   hESDEposEnegLikeSignBackInvMassPt = new TH2F*[fnCuts];
   hESDConvGammaZR                 = new TH2F*[fnCuts];
   
   }
   
   
   
   hESDMotherInvMassPt             = new TH2F*[fnCuts];
   
   if(fDoChicAnalysis) {
   hESDPi0MotherInvMassPt          = new TH2F*[fnCuts];
   hESDPi0MotherDiffInvMassPt      = new TH2F*[fnCuts];
   hESDPi0MotherDiffLimInvMassPt   = new TH2F*[fnCuts];
   }
   
   
   hESDMotherBackInvMassPt         = new TH2F*[fnCuts];


   for(Int_t iCut = 0; iCut<fnCuts;iCut++){


      TString cutstringElectron =((AliDalitzElectronCuts*)fCutElectronArray->At(iCut))->GetCutNumber();
      TString cutstringMeson= ((AliConversionMesonCuts*)fCutMesonArray->At(iCut))->GetCutNumber();
      TString cutstringGamma = ((AliConversionCuts*)fCutGammaArray->At(iCut))->GetCutNumber();

      fCutFolder[iCut] = new TList();
      fCutFolder[iCut]->SetName(Form("Cut Number %s_%s_%s",cutstringGamma.Data(),cutstringElectron.Data(),cutstringMeson.Data()));
      fCutFolder[iCut]->SetOwner(kTRUE);
      fOutputContainer->Add(fCutFolder[iCut]);

      fESDList[iCut] = new TList();
      fESDList[iCut]->SetName(Form("%s_%s_%s ESD histograms",cutstringGamma.Data(),cutstringElectron.Data(),cutstringMeson.Data()));
      fESDList[iCut]->SetOwner(kTRUE);
      
     
      hNEvents[iCut] = new TH1I("NEvents","NEvents",9,-0.5,8.5);
      hNEvents[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
      hNEvents[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
      hNEvents[iCut]->GetXaxis()->SetBinLabel(3,"Missing MC");
      hNEvents[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
      hNEvents[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
      hNEvents[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
      hNEvents[iCut]->GetXaxis()->SetBinLabel(7,"Pile-Up");
      hNEvents[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
      hNEvents[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
      fESDList[iCut]->Add(hNEvents[iCut]);



      if(fIsHeavyIon) hNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",3000,0,3000);
      else hNGoodESDTracks[iCut] = new TH1I("GoodESDTracks","GoodESDTracks",200,0,200);
      fESDList[iCut]->Add(hNGoodESDTracks[iCut]);

      hEtaShift[iCut] = new TProfile("Eta Shift","Eta Shift",1, -0.5,0.5);
      fESDList[iCut]->Add(hEtaShift[iCut]);

      hESDConvGammaPt[iCut] = new TH1F("ESD_ConvGamma_Pt","ESD_ConvGamma_Pt",250,0,25);
      fESDList[iCut]->Add(hESDConvGammaPt[iCut]);
      
      hESDConvGammaEta[iCut] = new TH1F("ESD_ConvGamma_Eta","ESD_ConvGamma_Eta",600,-1.5,1.5);
      fESDList[iCut]->Add(hESDConvGammaEta[iCut]);
      
      hESDDalitzElectronPt[iCut] = new TH1F("ESD_DalitzElectron_Pt","ESD_DalitzElectron_Pt",1000,0,25);
      fESDList[iCut]->Add(hESDDalitzElectronPt[iCut]);

      hESDDalitzPositronPt[iCut] = new TH1F("ESD_DalitzPositron_Pt","ESD_DalitzPositron_Pt",1000,0,25);
      fESDList[iCut]->Add(hESDDalitzPositronPt[iCut]);
      
      
      hESDDalitzElectronPhi[iCut] = new TH1F("ESD_DalitzElectron_Phi","ESD_DalitzElectron_Phi",360,0,2*TMath::Pi());
      fESDList[iCut]->Add(hESDDalitzElectronPhi[iCut]);

      hESDDalitzPositronPhi[iCut] = new TH1F("ESD_DalitzPositron_Phi","ESD_DalitzPositron_Phi",360,0,2*TMath::Pi());
      fESDList[iCut]->Add(hESDDalitzPositronPhi[iCut]);
      
     
      
      if ( fDoMesonQA ) {

      fQAFolder[iCut] = new TList();
      fQAFolder[iCut]->SetName(Form("%s_%s_%s QA histograms",cutstringGamma.Data(),cutstringElectron.Data(),cutstringMeson.Data()));
      fQAFolder[iCut]->SetOwner(kTRUE);
      	 
	      
	
      const Int_t kDCABins=62;
      Double_t binsDCADummy[63]={-3.0,-2.7,-2.4,-2.1,-1.8,-1.5,-1.2,-0.9,-0.6,-0.3,-0.25,-0.2,-0.19,-0.18,-0.17,-0.16,-0.15,-0.14,-0.13,-0.12,-0.11,-0.10,-0.09,-0.08,-0.07,-0.06,-0.05,-0.04,-0.03,-0.02,-0.01,0.0,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.1,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.2,0.25,0.3,0.6,0.9,1.2,1.5,1.8,2.1,2.4,2.7,3.0};

      const Int_t kPtBins=110;
      Double_t binsPtDummy[kPtBins+1];
      const Int_t kPBins = 109;
      Double_t binsPDummy[kPBins+1];
      binsPtDummy[0]=0.0;
      binsPDummy[0]=0.05;
     
        for(Int_t i=1;i<kPtBins+1;i++)
        {
                if(binsPtDummy[i-1]+0.05<1.01)
                        binsPtDummy[i]=binsPtDummy[i-1]+0.05;
                else
                        binsPtDummy[i]=binsPtDummy[i-1]+0.1;
		
        }
        for(Int_t i=1; i<kPBins+1;i++){
		  
		  if( binsPDummy[i-1]+0.05<1.01)
		        binsPDummy[i] = binsPDummy[i-1]+0.05;
		  else
			binsPDummy[i] = binsPDummy[i-1]+0.1;
		
	}
        
      
        
      hESDConvGammaZR[iCut]= new TH2F("ESD_ConvGamma_ConversionPoint_ZR","ESD_ConvGamma_ConversionPoint_ZR",1200,-150,150,480,0,120);
      fQAFolder[iCut]->Add(hESDConvGammaZR[iCut]);
     
      hESDDalitzElectronAfterPt[iCut] = new TH1F("ESD_DalitzElectron_After_Pt","ESD_DalitzElectron_After_Pt",1000,0,25);
      fQAFolder[iCut]->Add(hESDDalitzElectronAfterPt[iCut]);

      hESDDalitzPositronAfterPt[iCut] = new TH1F("ESD_DalitzPositron_After_Pt","ESD_DalitzPositron_After_Pt",1000,0,25);
      fQAFolder[iCut]->Add(hESDDalitzPositronAfterPt[iCut]);
        
      hESDDalitzElectronAfterEta[iCut] = new TH1F("ESD_DalitzElectron_After_Eta","ESD_DalitzElectron_After_Eta",600,-1.5,1.5);
      fQAFolder[iCut]->Add(hESDDalitzElectronAfterEta[iCut]);
      
      hESDDalitzPositronAfterEta[iCut] = new TH1F("ESD_DalitzPositron_After_Eta","ESD_DalitzElectron_After_Eta",600,-1.5,1.5);
      fQAFolder[iCut]->Add(hESDDalitzPositronAfterEta[iCut]);
      
                  
      hESDDalitzElectronAfterPhi[iCut] = new TH1F("ESD_DalitzElectron_After_Phi","ESD_DalitzElectron_After_Phi",360,0,2*TMath::Pi());
      fQAFolder[iCut]->Add(hESDDalitzElectronAfterPhi[iCut]);

      hESDDalitzPositronAfterPhi[iCut] = new TH1F("ESD_DalitzPositron_After_Phi","ESD_DalitzPositron_After_Phi",360,0,2*TMath::Pi());
      fQAFolder[iCut]->Add(hESDDalitzPositronAfterPhi[iCut]);
      
      hESDDalitzElectronAfterNClsITS[iCut]  = new TH1F("ESD_DalitzElectron_After_NClsITS","ESD_DalitzElectron_After_NClsITS",6,0.,6.);
      fQAFolder[iCut]->Add(hESDDalitzElectronAfterNClsITS[iCut]);
      
      hESDDalitzPositronAfterNClsITS[iCut]  = new TH1F("ESD_DalitzPositron_After_NClsITS","ESD_DalitzPositron_After_NClsITS",6,0.,6.);
      fQAFolder[iCut]->Add(hESDDalitzPositronAfterNClsITS[iCut]);
      
      
      hESDDalitzElectronAfterNFindClsTPC[iCut]  = new TH2F("ESD_DalitzElectron_After_NFindClsTPC","ESD_DalitzElectron_After_NFindClsTPC",50,0,1,kPtBins,binsPtDummy);
      fQAFolder[iCut]->Add(hESDDalitzElectronAfterNFindClsTPC[iCut]);
      
      hESDDalitzPositronAfterNFindClsTPC[iCut]  = new TH2F("ESD_DalitzPositron_After_NFindClsTPC","ESD_DalitzPositron_After_NFindClsTPC",50,0,1,kPtBins,binsPtDummy);
      fQAFolder[iCut]->Add(hESDDalitzPositronAfterNFindClsTPC[iCut]);
      
      
      hESDDalitzElectronAfterNClsTPC[iCut]  = new TH2F("ESD_DalitzElectron_After_NClsTPC","ESD_DalitzElectron_After_NClsTPC",200,0,200,kPtBins,binsPtDummy);
      fQAFolder[iCut]->Add(hESDDalitzElectronAfterNClsTPC[iCut]);
      
      hESDDalitzPositronAfterNClsTPC[iCut]  = new TH2F("ESD_DalitzPositron_After_NClsTPC","ESD_DalitzPositron_After_NClsTPC",200,0,200,kPtBins,binsPtDummy);
      fQAFolder[iCut]->Add(hESDDalitzPositronAfterNClsTPC[iCut]);
      
         
      
      hESDDalitzPosEleAfterDCAxy[iCut] = new TH2F("ESD_DalitzPosEle_After_DCAxy","ESD_DalitzPosEle_After_DCAxy",kDCABins,binsDCADummy,kPtBins,binsPtDummy);
      fQAFolder[iCut]->Add(hESDDalitzPosEleAfterDCAxy[iCut]);
      
      hESDDalitzPosEleAfterDCAz[iCut]  = new TH2F("ESD_DalitzPosEle_After_DCAz","ESD_DalitzPosEle_After_DCAz",kDCABins,binsDCADummy,kPtBins,binsPtDummy);
      fQAFolder[iCut]->Add(hESDDalitzPosEleAfterDCAz[iCut]);
      
      hESDDalitzElectronAfterTPCdEdxVsP[iCut] = new TH2F("ESD_DalitzElectron_After_TPCdEdxVsP","ESD_DalitzElectron_After_TPCdEdxVsP_After_TPCdEdx",kPBins,binsPDummy,200,-10,10);
      fQAFolder[iCut]->Add(hESDDalitzElectronAfterTPCdEdxVsP[iCut]);
      
      hESDDalitzPositronAfterTPCdEdxVsP[iCut] = new TH2F("ESD_DalitzPositron_After_TPCdEdxVsP","ESD_DalitzPositron_After_TPCdEdxVsP",kPBins,binsPDummy,200,-10,10);
      fQAFolder[iCut]->Add(hESDDalitzPositronAfterTPCdEdxVsP[iCut]);
      
      hESDDalitzElectronAfterTPCdEdxSignalVsP[iCut] =new TH2F("ESD_DalitzElectron_After_TPCdEdxSignalVsP","ESD_DalitzElectron_After_TPCdEdxSignalVsP" ,kPBins,binsPDummy,200,0.0,200);
      fQAFolder[iCut]->Add(hESDDalitzElectronAfterTPCdEdxSignalVsP[iCut]);  
      
      hESDDalitzPositronAfterTPCdEdxSignalVsP[iCut] =new TH2F("ESD_DalitzPositron_After_TPCdEdxSignalVsP","ESD_DalitzPositron_After_TPCdEdxSignalVsP" ,kPBins,binsPDummy,200,0.0,200);
      fQAFolder[iCut]->Add(hESDDalitzPositronAfterTPCdEdxSignalVsP[iCut]); 
      
      hESDDalitzElectronAfterTPCdEdxVsEta[iCut] = new TH2F("ESD_DalitzElectron_After_TPCdEdxVsEta","ESD_DalitzElectron_After_TPCdEdxVsEta",140,-1.4,1.4,200,-10,10);
      fQAFolder[iCut]->Add(hESDDalitzElectronAfterTPCdEdxVsEta[iCut]);
      
      hESDDalitzPositronAfterTPCdEdxVsEta[iCut] = new  TH2F("ESD_DalitzPositron_After_TPCdEdxVsEta","ESD_DalitzPositron_After_TPCdEdxVsEta",140,-1.4,1.4,200,-10,10);
      fQAFolder[iCut]->Add(hESDDalitzPositronAfterTPCdEdxVsEta[iCut]);
      
      hESDDalitzElectronAfterTPCdEdxVsPhi[iCut] = new TH2F("ESD_DalitzElectron_After_TPCdEdxVsPhi","ESD_DalitzElectron_After_TPCdEdxVsPhi",180,0,2*TMath::Pi(),200,-10,10);
      fQAFolder[iCut]->Add(hESDDalitzElectronAfterTPCdEdxVsPhi[iCut]);
      
      hESDDalitzPositronAfterTPCdEdxVsPhi[iCut] = new TH2F("ESD_DalitzPositron_After_TPCdEdxVsPhi","ESD_DalitzPositron_After_TPCdEdxVsPhi",180,0,2*TMath::Pi(),200,-10,10);
      fQAFolder[iCut]->Add(hESDDalitzPositronAfterTPCdEdxVsPhi[iCut]);
            
      hESDMotherPhi[iCut] = new TH1F("ESD_DalitzMother_Phi","ESD_DalitzMother_Phi",360,0,2*TMath::Pi());
      fQAFolder[iCut]->Add(hESDMotherPhi[iCut]);
      
      hESDEposEnegPsiPairDPhi[iCut] = new TH2F("ESD_EposEneg_PsiPair_DPhi","ESD_EposEneg_PsiPair_DPhi", 100, -1.0,1.0,100,-1.0,1.0 );
      fQAFolder[iCut]->Add(hESDEposEnegPsiPairDPhi[iCut]);

      hESDEposEnegInvMassPt[iCut] = new TH2F("ESD_EposEneg_InvMassPt","ESD_EposEneg_InvMassPt",5000,0.,5.,100,0.,10.);
      fQAFolder[iCut]->Add(hESDEposEnegInvMassPt[iCut]);
      
      hESDEposEnegLikeSignBackInvMassPt[iCut]  = new TH2F("ESD_EposEneg_LikeSignBack_InvMassPt","ESD_EposEneg_LikeSignBack_InvMassPt",5000,0.,5.,100,0.,10.);
      fQAFolder[iCut]->Add(hESDEposEnegLikeSignBackInvMassPt[iCut]);
      
      
      TAxis *AxisAfter = hESDDalitzElectronAfterTPCdEdxVsP[iCut]->GetXaxis(); 
      Int_t bins = AxisAfter->GetNbins();
      Double_t from = AxisAfter->GetXmin();
      Double_t to = AxisAfter->GetXmax();
      Double_t *newBins = new Double_t[bins+1];
      newBins[0] = from;
      Double_t factor = TMath::Power(to/from, 1./bins);
      for(Int_t i=1; i<=bins; ++i) newBins[i] = factor * newBins[i-1];

      AxisAfter->Set(bins, newBins);
      AxisAfter = hESDDalitzElectronAfterTPCdEdxSignalVsP[iCut]->GetXaxis(); 
      AxisAfter->Set(bins, newBins);
      
      AxisAfter = hESDDalitzPositronAfterTPCdEdxVsP[iCut]->GetXaxis();
      AxisAfter->Set(bins, newBins);
      
      AxisAfter = hESDDalitzPositronAfterTPCdEdxSignalVsP[iCut]->GetXaxis();
      AxisAfter->Set(bins,newBins);
		    
	           

      delete [] newBins;
      
      fCutFolder[iCut]->Add(fQAFolder[iCut]);
      
           
      
      }
      
      
             
     

      hESDMotherInvMassPt[iCut] = new TH2F("ESD_DalitzMother_InvMass_Pt","ESD_DalitzMother_InvMass_Pt",800,0,0.8,250,0,25);
      fESDList[iCut]->Add(hESDMotherInvMassPt[iCut]);
										  
       
      if( fDoChicAnalysis) {
                  
      hESDPi0MotherInvMassPt[iCut] = new TH2F("ESD_Pi0Mother_InvMass_Pt","ESD_Pi0Mother_InvMass_Pt",4000,0,4,250,0,25);
      fESDList[iCut]->Add(hESDPi0MotherInvMassPt[iCut]);
  
      hESDPi0MotherDiffInvMassPt[iCut] = new TH2F("ESD_Pi0Mother_DiffInvMass_Pt","ESD_Pi0Mother_DiffInvMass_Pt",2000,0,2,250,0,25);
      fESDList[iCut]->Add(hESDPi0MotherDiffInvMassPt[iCut]);

      hESDPi0MotherDiffLimInvMassPt[iCut] = new TH2F("ESD_Pi0Mother_DiffLimInvMass_Pt","ESD_Pi0Mother_DiffLimInvMass_Pt",2000,0,2,250,0,25);
      fESDList[iCut]->Add(hESDPi0MotherDiffLimInvMassPt[iCut]);
      
      }


      hESDMotherBackInvMassPt[iCut] = new TH2F("ESD_DalitzBackground_InvMass_Pt","ESD_DalitzBackground_InvMass_Pt",800,0,0.8,250,0,25);
      fESDList[iCut]->Add(hESDMotherBackInvMassPt[iCut]);

           
      
      fCutFolder[iCut]->Add(fESDList[iCut]);
      
      

   }


   InitBack(); // Init Background Handler


   //if(MCEvent()){
    if( fIsMC ){
      // MC Histogramms
      fMCList = new TList*[fnCuts];
      // True Histogramms
      fTrueList = new TList*[fnCuts];
      hESDTrueConvGammaPt = new TH1F*[fnCuts];
      hESDTruePositronPt  = new TH1F*[fnCuts];
      hESDTrueElectronPt  = new TH1F*[fnCuts];
      hESDTrueSecConvGammaPt = new TH1F*[fnCuts];
      hESDTrueSecPositronPt  = new TH1F*[fnCuts];
      hESDTrueSecElectronPt  = new TH1F*[fnCuts];
      hESDTruePi0DalitzConvGammaPt = new TH1F*[fnCuts];
      hESDTruePi0DalitzPositronPt  = new TH1F*[fnCuts];
      hESDTruePi0DalitzElectronPt  = new TH1F*[fnCuts];
      hESDTruePi0DalitzSecConvGammaPt = new TH1F*[fnCuts];
      hESDTruePi0DalitzSecPositronPt  = new TH1F*[fnCuts];
      hESDTruePi0DalitzSecElectronPt  = new TH1F*[fnCuts];
      //if(fDoMesonAnalysis){
      hMCAllGammaPt  = new TH1F*[fnCuts];
      hMCConvGammaPt = new TH1F*[fnCuts];
      hMCConvGammaRSPt = new TH1F*[fnCuts];
      hMCAllPositronsPt = new TH1F*[fnCuts];
      hMCAllElectronsPt = new TH1F*[fnCuts];
      hMCPi0DalitzGammaPt    = new TH1F*[fnCuts];
      hMCPi0DalitzElectronPt = new TH1F*[fnCuts];
      hMCPi0DalitzPositronPt = new TH1F*[fnCuts];
 
      hMCPi0Pt = new TH1F*[fnCuts];
      hMCPi0GGPt =  new TH1F*[fnCuts];
      hMCEtaPt = new TH1F*[fnCuts];
      hMCEtaGGPt = new TH1F*[fnCuts];
      hMCPi0InAccPt = new TH1F*[fnCuts];
      hMCEtaInAccPt = new TH1F*[fnCuts];
			hMCChiCPt = new TH1F*[fnCuts];
			hMCChiCInAccPt = new TH1F*[fnCuts];
			
			
      if ( fDoMesonQA ) {
	hESDEposEnegTruePi0DalitzInvMassPt           = new TH2F*[fnCuts];
        hESDEposEnegTruePi0DalitzPsiPairDPhi         = new TH2F*[fnCuts];
	hESDEposEnegTrueEtaDalitzInvMassPt           = new TH2F*[fnCuts];
        hESDEposEnegTrueEtaDalitzPsiPairDPhi         = new TH2F*[fnCuts];
	hESDEposEnegTruePhotonInvMassPt              = new TH2F*[fnCuts];
        hESDEposEnegTruePhotonPsiPairDPhi            = new TH2F*[fnCuts];
	hESDEposEnegTrueJPsiInvMassPt                = new TH2F*[fnCuts];
      }
      
      
      if( fDoChicAnalysis ){
      hESDTrueMotherChiCInvMassPt = new TH2F*[fnCuts];
      hESDTrueMotherChiCDiffInvMassPt = new TH2F*[fnCuts];
      }
      
      
      hESDTrueMotherInvMassPt = new TH2F*[fnCuts];
      hESDTrueMotherDalitzInvMassPt = new TH2F*[fnCuts];
      hESDTrueMotherPi0GGInvMassPt = new TH2F*[fnCuts];
      hESDTruePrimaryMotherPi0GGInvMassPt = new TH2F*[fnCuts];
      hESDTrueSecondaryMotherPi0GGInvMassPt = new TH2F*[fnCuts];
      hESDTruePrimaryPi0DalitzESDPtMCPt = new TH2F*[fnCuts];
      hESDTruePrimaryMotherInvMassMCPt = new TH2F*[fnCuts];
      hESDTruePrimaryMotherInvMassPt   = new TH2F*[fnCuts];
      hESDTruePrimaryMotherW0WeightingInvMassPt = new TH2F*[fnCuts];
      hESDTrueSecondaryMotherInvMassPt = new TH2F*[fnCuts];
      hESDTrueSecondaryMotherFromK0sInvMassPt = new TH2F*[fnCuts];
      hESDTrueBckGGInvMassPt = new TH2F*[fnCuts];
      hESDTrueBckContInvMassPt = new TH2F*[fnCuts];
      hESDTrueMotherGGInvMassPt = new TH2F*[fnCuts];
      //}

      for(Int_t iCut = 0; iCut<fnCuts;iCut++){
         TString cutstringElectron =((AliDalitzElectronCuts*)fCutElectronArray->At(iCut))->GetCutNumber();
         TString cutstringMeson= ((AliConversionMesonCuts*)fCutMesonArray->At(iCut))->GetCutNumber();
         TString cutstringGamma = ((AliConversionCuts*)fCutGammaArray->At(iCut))->GetCutNumber();

         fMCList[iCut] = new TList();
         fMCList[iCut]->SetName(Form("%s_%s_%s MC histograms",cutstringGamma.Data(),cutstringElectron.Data(),cutstringMeson.Data()));
         fMCList[iCut]->SetOwner(kTRUE);
         fCutFolder[iCut]->Add(fMCList[iCut]);


	hMCAllGammaPt[iCut] = new TH1F("MC_AllGamma_Pt","MC_AllGamma_Pt",250,0,25);
	fMCList[iCut]->Add(hMCAllGammaPt[iCut]);
	
	hMCConvGammaPt[iCut] = new TH1F("MC_ConvGamma_Pt","MC_ConvGamma_Pt",250,0,25);
	fMCList[iCut]->Add(hMCConvGammaPt[iCut]);
	
	hMCConvGammaRSPt[iCut] = new TH1F("MC_ConvGamma_RS_Pt","MC_ConvGamma_RS_Pt",250,0,25);
        fMCList[iCut]->Add(hMCConvGammaRSPt[iCut]);
	
				 
	hMCAllPositronsPt[iCut] = new TH1F("MC_AllPositrons_Pt","MC_AllPositrons_Pt",1000,0,25);
	fMCList[iCut]->Add(hMCAllPositronsPt[iCut]);
				 
				 hMCAllElectronsPt[iCut] = new TH1F("MC_AllElectrons_Pt","MC_AllElectrons_Pt",1000,0,25);
				 fMCList[iCut]->Add(hMCAllElectronsPt[iCut]);
				 
				 hMCPi0DalitzGammaPt[iCut] = new TH1F("MC_Pi0DalitzGamma_Pt","MC_Pi0DalitzGamma_Pt",250,0,25);
                                 hMCPi0DalitzGammaPt[iCut]->Sumw2();
				 fMCList[iCut]->Add(hMCPi0DalitzGammaPt[iCut]);
				 
				 hMCPi0DalitzPositronPt[iCut] = new TH1F("MC_Pi0DalitzPositron_Pt","MC_Pi0DalitzPositron_Pt",1000,0,25);
                                 hMCPi0DalitzPositronPt[iCut]->Sumw2();
				 fMCList[iCut]->Add(hMCPi0DalitzPositronPt[iCut]);
				 
				 hMCPi0DalitzElectronPt[iCut] = new TH1F("MC_Pi0DalitzElectron_Pt","MC_Pi0DalitzElectron_Pt",1000,0,25);
                                 hMCPi0DalitzElectronPt[iCut]->Sumw2();
				 fMCList[iCut]->Add(hMCPi0DalitzElectronPt[iCut]);
				 
				 
				 hMCPi0Pt[iCut] = new TH1F("MC_Pi0_Pt","MC_Pi0_Pt",250,0,25);
                                 hMCPi0Pt[iCut]->Sumw2();
				 fMCList[iCut]->Add(hMCPi0Pt[iCut]);
				 
				 hMCPi0GGPt[iCut] = new TH1F("MC_Pi0_GG_Pt","MC_Pi0_GG_Pt",250,0,25);
                                 hMCPi0GGPt[iCut]->Sumw2();
				 fMCList[iCut]->Add(hMCPi0GGPt[iCut]);
				 
				 hMCEtaPt[iCut] = new TH1F("MC_Eta_Pt","MC_Eta_Pt",250,0,25);
                                 hMCEtaPt[iCut]->Sumw2();
				 fMCList[iCut]->Add(hMCEtaPt[iCut]);

				 hMCEtaGGPt[iCut] = new TH1F("MC_Eta_GG_Pt","MC_Eta_GG_Pt",250,0,25);
                                 hMCEtaGGPt[iCut]->Sumw2();
				 fMCList[iCut]->Add(hMCEtaGGPt[iCut]);
				 
				 hMCPi0InAccPt[iCut] = new TH1F("MC_Pi0DalitzInAcc_Pt","MC_Pi0DalitzInAcc_Pt",250,0,25);
                                 hMCPi0InAccPt[iCut]->Sumw2();
				 fMCList[iCut]->Add(hMCPi0InAccPt[iCut]);

				 hMCEtaInAccPt[iCut] = new TH1F("MC_EtaDalitzInAcc_Pt","MC_EtaDalitzInAcc_Pt",250,0,25);
                                 hMCEtaInAccPt[iCut]->Sumw2();
				 fMCList[iCut]->Add(hMCEtaInAccPt[iCut]);

				 hMCChiCPt[iCut] = new TH1F("MC_ChiC_Pt","MC_ChiC_Pt",250,0,25);
				 fMCList[iCut]->Add(hMCChiCPt[iCut]);

				 hMCChiCInAccPt[iCut] = new TH1F("MC_ChiCInAcc_Pt","MC_ChiCInAcc_Pt",250,0,25);
				 fMCList[iCut]->Add(hMCChiCInAccPt[iCut]);

         fTrueList[iCut] = new TList();
         fTrueList[iCut]->SetName(Form("%s_%s_%s True histograms",cutstringGamma.Data(),cutstringElectron.Data(),cutstringMeson.Data()));
         fTrueList[iCut]->SetOwner(kTRUE);
         fCutFolder[iCut]->Add(fTrueList[iCut]);

	 if ( fDoMesonQA ) {


	 hESDEposEnegTruePi0DalitzInvMassPt[iCut] = new TH2F("ESD_EposEneg_TruePi0Dalitz_InvMassPt","ESD_EposEneg_TruePi0Dalitz_InvMassPt",5000,0.,5.,100,0.,10.);
	 fTrueList[iCut]->Add(hESDEposEnegTruePi0DalitzInvMassPt[iCut]);

         hESDEposEnegTruePi0DalitzPsiPairDPhi[iCut] = new TH2F("ESD_EposEneg_TruePi0Dalitz_PsiPair_DPhi","ESD_EposEneg_TruePi0Dalitz_PsiPair_DPhi", 100, -1.0,1.0,100,-1.0,1.0 );
         fTrueList[iCut]->Add(hESDEposEnegTruePi0DalitzPsiPairDPhi[iCut]);

	 
	 hESDEposEnegTrueEtaDalitzInvMassPt[iCut] = new TH2F("ESD_EposEneg_TrueEtaDalitz_InvMassPt","ESD_EposEneg_TrueEtaDalitz_InvMassPt",5000,0.,5.,100,0.,10.);
	 fTrueList[iCut]->Add(hESDEposEnegTrueEtaDalitzInvMassPt[iCut]);

         hESDEposEnegTrueEtaDalitzPsiPairDPhi[iCut] = new TH2F("ESD_EposEneg_TrueEtaDalitz_PsiPair_DPhi","ESD_EposEneg_TrueEtaDalitz_PsiPair_DPhi", 100, -1.0,1.0,100,-1.0,1.0 );
         fTrueList[iCut]->Add(hESDEposEnegTrueEtaDalitzPsiPairDPhi[iCut]);
	 
	 hESDEposEnegTruePhotonInvMassPt[iCut] = new TH2F("ESD_EposEneg_TruePhoton_InvMassPt","ESD_EposEneg_TruePhoton_InvMassPt",5000,0.,5.,100,0.,10.);
	 fTrueList[iCut]->Add(hESDEposEnegTruePhotonInvMassPt[iCut]);
        
         hESDEposEnegTruePhotonPsiPairDPhi[iCut] = new TH2F("ESD_EposEneg_TruePhoton_PsiPair_DPhi","ESD_EposEneg_TruePhoton_PsiPair_DPhi", 100, -1.0,1.0,100,-1.0,1.0 );
         fTrueList[iCut]->Add(hESDEposEnegTruePhotonPsiPairDPhi[iCut]);
 
	 hESDEposEnegTrueJPsiInvMassPt[iCut] = new TH2F("ESD_EposEneg_TrueJPsi_InvMassPt","ESD_EposEneg_TrueJPsi_InvMassPt",5000,0.,5.,100,0.,10.);
	 fTrueList[iCut]->Add(hESDEposEnegTrueJPsiInvMassPt[iCut]);



	 }


         hESDTruePositronPt[iCut] = new TH1F("ESD_TruePositron_Pt","ESD_TruePositron_Pt",1000,0,25);
         fTrueList[iCut]->Add(hESDTruePositronPt[iCut]);

         hESDTrueElectronPt[iCut] = new TH1F("ESD_TrueElectron_Pt","ESD_TrueElectron_Pt",1000,0,25);
         fTrueList[iCut]->Add(hESDTrueElectronPt[iCut]);
	 
	 hESDTrueSecPositronPt[iCut] = new TH1F("ESD_TrueSecPositron_Pt","ESD_TrueSecPositron_Pt",1000,0,25);
         fTrueList[iCut]->Add(hESDTrueSecPositronPt[iCut]);

         hESDTrueSecElectronPt[iCut] = new TH1F("ESD_TrueSecElectron_Pt","ESD_TrueSecElectron_Pt",1000,0,25);
         fTrueList[iCut]->Add(hESDTrueSecElectronPt[iCut]); 

         hESDTrueConvGammaPt[iCut] = new TH1F("ESD_TrueConvGamma_Pt","ESD_TrueConvGamma_Pt",250,0,25);
         fTrueList[iCut]->Add(hESDTrueConvGammaPt[iCut]);
	 
	 hESDTrueSecConvGammaPt[iCut] = new TH1F("ESD_TrueSecConvGamma_Pt","ESD_TrueSecConvGamma_Pt",250,0,25);
         fTrueList[iCut]->Add(hESDTrueSecConvGammaPt[iCut]);

         hESDTruePi0DalitzConvGammaPt[iCut] = new TH1F("ESD_TruePi0DalitzConvGamma_Pt","ESD_TruePi0DalitzConvGamma_Pt",250,0,25);
         fTrueList[iCut]->Add(hESDTruePi0DalitzConvGammaPt[iCut]);

         hESDTruePi0DalitzElectronPt[iCut] = new TH1F("ESD_TruePi0DalitzElectron_Pt","ESD_TruePi0DalitzElectron_Pt",1000,0,25);
         fTrueList[iCut]->Add(hESDTruePi0DalitzElectronPt[iCut]);

         hESDTruePi0DalitzPositronPt[iCut] = new TH1F("ESD_TruePi0DalitzPositron_Pt","ESD_TruePi0DalitzPositron_Pt",1000,0,25);
         fTrueList[iCut]->Add(hESDTruePi0DalitzPositronPt[iCut]);
	 
	 hESDTruePi0DalitzSecConvGammaPt[iCut] = new TH1F("ESD_TruePi0DalitzSecConvGamma_Pt","ESD_TruePi0DalitzSecConvGamma_Pt",250,0,25);
         fTrueList[iCut]->Add(hESDTruePi0DalitzSecConvGammaPt[iCut]);
	 
	 hESDTruePi0DalitzSecElectronPt[iCut] = new TH1F("ESD_TruePi0DalitzSecElectron_Pt","ESD_TruePi0DalitzSecElectron_Pt",1000,0,25);
         fTrueList[iCut]->Add(hESDTruePi0DalitzSecElectronPt[iCut]);

         hESDTruePi0DalitzSecPositronPt[iCut] = new TH1F("ESD_TruePi0DalitzSecPositron_Pt","ESD_TruePi0DalitzSecPositron_Pt",1000,0,25);
         fTrueList[iCut]->Add(hESDTruePi0DalitzSecPositronPt[iCut]);

	 if( fDoChicAnalysis) { 
	   
	 hESDTrueMotherChiCInvMassPt[iCut] = new TH2F("ESD_TrueMotherChiC_InvMass_Pt","ESD_TrueMotherChiC_InvMass_Pt",4000,0,4,250,0,25);
         fTrueList[iCut]->Add(hESDTrueMotherChiCInvMassPt[iCut]);

	 hESDTrueMotherChiCDiffInvMassPt[iCut] = new TH2F("ESD_TrueMotherChiCDiff_InvMass_Pt","ESD_TrueMotherChiCDiff_InvMass_Pt",2000,0,2,250,0,25);
         fTrueList[iCut]->Add(hESDTrueMotherChiCDiffInvMassPt[iCut]);
	 
	 }

	 hESDTrueMotherInvMassPt[iCut] = new TH2F("ESD_TrueMother_InvMass_Pt","ESD_TrueMother_InvMass_Pt",800,0,0.8,250,0,25);
         hESDTrueMotherInvMassPt[iCut]->Sumw2();
         fTrueList[iCut]->Add(hESDTrueMotherInvMassPt[iCut]);
	 
	 hESDTrueMotherDalitzInvMassPt[iCut] = new TH2F("ESD_TrueMother_Dalitz_InvMass_Pt","ESD_TrueMother_Dalitz_InvMass_Pt",800,0,0.8,250,0,25);
         hESDTrueMotherDalitzInvMassPt[iCut]->Sumw2();
         fTrueList[iCut]->Add(hESDTrueMotherDalitzInvMassPt[iCut]);
	 
	 
	 
	 

         hESDTrueMotherPi0GGInvMassPt[iCut] = new TH2F("ESD_TrueMotherPi0GG_InvMass_Pt","ESD_TrueMotherPi0GG_InvMass_Pt",800,0,0.8,250,0,25);
         hESDTrueMotherPi0GGInvMassPt[iCut]->Sumw2();
         fTrueList[iCut]->Add(hESDTrueMotherPi0GGInvMassPt[iCut]);

         hESDTruePrimaryMotherPi0GGInvMassPt[iCut] = new TH2F("ESD_TruePrimaryMotherPi0GG_InvMass_Pt","ESD_TruePrimaryMotherPi0GG_InvMass_Pt",800,0,0.8,250,0,25);
         hESDTruePrimaryMotherPi0GGInvMassPt[iCut]->Sumw2();
         fTrueList[iCut]->Add(hESDTruePrimaryMotherPi0GGInvMassPt[iCut]);

         hESDTrueSecondaryMotherPi0GGInvMassPt[iCut] = new TH2F("ESD_TrueSecondaryMotherPi0GG_InvMass_Pt","ESD_TrueSecondaryMotherPi0GG_InvMass_Pt",800,0,0.8,250,0,25);
         hESDTrueSecondaryMotherPi0GGInvMassPt[iCut]->Sumw2();
         fTrueList[iCut]->Add(hESDTrueSecondaryMotherPi0GGInvMassPt[iCut]);

         hESDTruePrimaryPi0DalitzESDPtMCPt[iCut] = new TH2F("ESD_TruePrimaryPi0Dalitz_ESDPt_MCPt","ESD_TruePrimaryPi0Dalitz_ESDPt_MCPt",250,0,25,250,0,25);
         hESDTruePrimaryPi0DalitzESDPtMCPt[iCut]->Sumw2();
         fTrueList[iCut]->Add(hESDTruePrimaryPi0DalitzESDPtMCPt[iCut]);
         hESDTruePrimaryMotherInvMassMCPt[iCut] = new TH2F("ESD_TruePrimaryMother_InvMass_MCPt","ESD_TrueDalitzPrimaryMother_InvMass_MCPt",800,0,0.8,250,0,25);
         hESDTruePrimaryMotherInvMassMCPt[iCut]->Sumw2();
         fTrueList[iCut]->Add(hESDTruePrimaryMotherInvMassMCPt[iCut]);
         hESDTruePrimaryMotherInvMassPt[iCut] = new TH2F("ESD_TruePrimaryMother_InvMass_Pt","ESD_TruePrimaryMother_InvMass_Pt",800,0,0.8,250,0,25);
         hESDTruePrimaryMotherInvMassPt[iCut]->Sumw2();
         fTrueList[iCut]->Add(hESDTruePrimaryMotherInvMassPt[iCut]);
	 hESDTruePrimaryMotherW0WeightingInvMassPt[iCut] = new TH2F("ESD_TruePrimaryMotherW0Weighting_InvMass_Pt","ESD_TruePrimaryMotherW0Weighting_InvMass_Pt",800,0,0.8,250,0,25);
         hESDTruePrimaryMotherW0WeightingInvMassPt[iCut]->Sumw2();
         fTrueList[iCut]->Add(hESDTruePrimaryMotherW0WeightingInvMassPt[iCut]); 

         hESDTrueSecondaryMotherInvMassPt[iCut] = new TH2F("ESD_TrueDalitzSecondaryMother_InvMass_Pt","ESD_TrueDalitzSecondaryMother_InvMass_Pt",800,0,0.8,250,0,25);
         hESDTrueSecondaryMotherInvMassPt[iCut]->Sumw2();
         fTrueList[iCut]->Add(hESDTrueSecondaryMotherInvMassPt[iCut]);
         // 				hESDTrueSecondaryMotherFromK0sInvMassPt[iCut] = new TH2F("ESD_TrueDalitzSecondaryMotherFromK0s_InvMass_Pt","ESD_TrueDalitzSecondaryMotherFromK0s_InvMass_Pt",1000,0,1,250,0,25);
         // 				fTrueList[iCut]->Add(hESDTrueSecondaryMotherFromK0sInvMassPt[iCut]);
         hESDTrueBckGGInvMassPt[iCut] = new TH2F("ESD_TrueDalitzBckGG_InvMass_Pt","ESD_TrueDalitzBckGG_InvMass_Pt",800,0,0.8,250,0,25);
         fTrueList[iCut]->Add(hESDTrueBckGGInvMassPt[iCut]);
         hESDTrueBckContInvMassPt[iCut] = new TH2F("ESD_TrueDalitzBckCont_InvMass_Pt","ESD_TrueDalitzBckCont_InvMass_Pt",800,0,0.8,250,0,25);
         fTrueList[iCut]->Add(hESDTrueBckContInvMassPt[iCut]);
         // 				hESDTrueMotherGGInvMassPt[iCut] = new TH2F("ESD_TrueGammaGamma_InvMass_Pt","ESD_TrueGammaGamma_InvMass_Pt",1000,0,1,250,0,25);
         // 				fTrueList[iCut]->Add(hESDTrueMotherGGInvMassPt[iCut]);

      }
   }
   
   fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1");
   if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader
     
   if(fV0Reader)
      if((AliConversionCuts*)fV0Reader->GetConversionCuts())
         if(((AliConversionCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms())
            fOutputContainer->Add(((AliConversionCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());
	 
	 
	 
   fElecSelector=(AliDalitzElectronSelector*)AliAnalysisManager::GetAnalysisManager()->GetTask("ElectronSelector");
   if(!fElecSelector){printf("Error: No ElectronSelector");return;} // GetV0Reader
     
    if( fElecSelector ){

      if ( ((AliDalitzElectronCuts*)fElecSelector->GetDalitzElectronCuts())->GetCutHistograms() ){
         fOutputContainer->Add( ((AliDalitzElectronCuts*)fElecSelector->GetDalitzElectronCuts())->GetCutHistograms() );
      }
   }  

   for(Int_t iCut = 0; iCut<fnCuts;iCut++){

      if( fCutElectronArray ){
         if( ((AliDalitzElectronCuts*)fCutElectronArray->At(iCut))->GetCutHistograms() ) {
            fCutFolder[iCut]->Add( ((AliDalitzElectronCuts*)fCutElectronArray->At(iCut))->GetCutHistograms() );
         }
      }

      if( fCutMesonArray  ) {
         if( ((AliConversionMesonCuts*)fCutMesonArray->At(iCut))->GetCutHistograms() ) {
            fCutFolder[iCut]->Add( ((AliConversionMesonCuts*)fCutMesonArray->At(iCut))->GetCutHistograms());
         }
      }

      if( fCutGammaArray ) {
         if( ((AliConversionCuts*)fCutGammaArray->At(iCut))->GetCutHistograms() ) {
            fCutFolder[iCut]->Add( ((AliConversionCuts*)fCutGammaArray->At(iCut))->GetCutHistograms()  );
         }
      }
   }

   PostData(1, fOutputContainer);

}

//______________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::UserExec(Option_t *)
{

   //
   // Execute analysis for current event
   //

   fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1");
   if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader


   Int_t eventQuality = ((AliConversionCuts*)fV0Reader->GetConversionCuts())->GetEventQuality();

   if(eventQuality == 2 || eventQuality == 3){// Event Not Accepted due to MC event missing or wrong trigger for V0ReaderV1
        for(Int_t iCut = 0; iCut<fnCuts; iCut++){
	   hNEvents[iCut]->Fill(eventQuality);
        }
      return;
   }


   fElecSelector=(AliDalitzElectronSelector*)AliAnalysisManager::GetAnalysisManager()->GetTask("ElectronSelector");
   if(!fElecSelector){printf("Error: No ElectronSelector");return;} // GetV0Reader


   if(fIsMC) 	      fMCEvent        =  MCEvent();
   fESDEvent        = (AliESDEvent*)InputEvent();
   fReaderGammas    = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut
   fSelectorElectronIndex = fElecSelector->GetReconstructedElectronsIndex(); // Electrons from default Cut
   fSelectorPositronIndex = fElecSelector->GetReconstructedPositronsIndex(); // Positrons from default Cut

   //CountESDTracks(); // Estimate Event Multiplicity
   fNumberOfESDTracks = fV0Reader->GetNumberOfPrimaryTracks();
   //AddTaskContainers(); //Add conatiner

   for(Int_t iCut = 0; iCut<fnCuts; iCut++){
      fiCut = iCut;

      Int_t eventNotAccepted =
         ((AliConversionCuts*)fCutGammaArray->At(iCut))
         ->IsEventAcceptedByConversionCut(fV0Reader->GetConversionCuts(),fInputEvent,fMCEvent,fIsHeavyIon);
	 
      if(eventNotAccepted){
         // 			cout << "event rejected due to wrong trigger: " <<eventNotAccepted << endl;
         hNEvents[iCut]->Fill(eventNotAccepted); // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
         continue;
      }

      if(eventQuality != 0){// Event Not Accepted
         // 			cout << "event rejected due to: " <<eventQuality << endl;
         hNEvents[iCut]->Fill(eventQuality);
         continue;
      }

      hNEvents[iCut]->Fill(eventQuality);

      hNGoodESDTracks[iCut]->Fill(fNumberOfESDTracks);

      if(fMCEvent){ // Process MC Particle
         
	
	
	fMCStack = fMCEvent->Stack();
	 
	  if(((AliConversionCuts*)fCutGammaArray->At(iCut))->GetSignalRejection() != 0){
             ((AliConversionCuts*)fCutGammaArray->At(iCut))->GetNotRejectedParticles(((AliConversionCuts*)fCutGammaArray->At(iCut))->GetSignalRejection(),
										     ((AliConversionCuts*)fCutGammaArray->At(iCut))->GetAcceptedHeader(),
										       fMCEvent);
	  } 
	 
         ProcessMCParticles();
      }

      ProcessPhotonCandidates(); // Process this cuts gammas
      ProcessElectronCandidates(); // Process this cuts gammas
      
      if(((AliConversionMesonCuts*)fCutMesonArray->At(iCut))->UseMCPSmearing() && fMCEvent){
	
            fUnsmearedPx = new Double_t[fGoodGammas->GetEntries()]; // Store unsmeared Momenta
            fUnsmearedPy = new Double_t[fGoodGammas->GetEntries()];
            fUnsmearedPz = new Double_t[fGoodGammas->GetEntries()];
            fUnsmearedE =  new Double_t[fGoodGammas->GetEntries()];

            for(Int_t gamma=0;gamma<fGoodGammas->GetEntries();gamma++){ // Smear the AODPhotons in MC
	      
               fUnsmearedPx[gamma] = ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->Px();
               fUnsmearedPy[gamma] = ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->Py();
               fUnsmearedPz[gamma] = ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->Pz();
               fUnsmearedE[gamma] =  ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->E();
               ((AliConversionMesonCuts*)fCutMesonArray->At(iCut))->SmearParticle(dynamic_cast<AliAODConversionPhoton*>(fGoodGammas->At(gamma)));
	       
            }           
        }
  

      CalculatePi0DalitzCandidates();
      CalculateBackground();
      UpdateEventByEventData();
      
      
      if(((AliConversionMesonCuts*)fCutMesonArray->At(iCut))->UseMCPSmearing() && fMCEvent){
	
	   for(Int_t gamma=0;gamma<fGoodGammas->GetEntries();gamma++){ // Smear the AODPhotons in MC
               ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->SetPx(fUnsmearedPx[gamma]); // Reset Unsmeared Momenta
               ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->SetPy(fUnsmearedPy[gamma]);
               ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->SetPz(fUnsmearedPz[gamma]);
               ((AliAODConversionPhoton*)fGoodGammas->At(gamma))->SetE(fUnsmearedE[gamma]);
            }
            delete[] fUnsmearedPx; fUnsmearedPx = 0x0;
            delete[] fUnsmearedPy; fUnsmearedPy = 0x0;
            delete[] fUnsmearedPz; fUnsmearedPz = 0x0;
            delete[] fUnsmearedE;  fUnsmearedE  = 0x0;
	    
       }


      fGoodGammas->Clear(); // delete this cuts good gammas
      fGoodVirtualGammas->Clear(); // delete this cuts good gammas
   }

   fSelectorElectronIndex.clear();
   fSelectorPositronIndex.clear();

   PostData( 1, fOutputContainer );
}

Bool_t AliAnalysisTaskGammaConvDalitzV1::Notify()
{
   for(Int_t iCut = 0; iCut<fnCuts;iCut++){
      

      if( !((AliConversionCuts*)fCutGammaArray->At(iCut))->GetDoEtaShift() ){

         /*if (((AliConversionCuts*)fCutGammaArray->At(iCut))->GetEtaShift() != 0.){
            printf("Error: Gamma Conversion Dalitz Task %s :: Eta Shift not requested but set to %f, reset to 00. \n\n",
                   (((AliConversionCuts*)fCutGammaArray->At(iCut))->GetCutNumber()).Data(),((AliConversionCuts*)fCutGammaArray->At(iCut))->GetEtaShift());
            ((AliConversionCuts*)fCutGammaArray->At(iCut))->SetEtaShift(0.);
            ((AliDalitzElectronCuts*)fCutElectronArray->At(iCut))->SetEtaShift(0.);
      
            
         }*/
         hEtaShift[iCut]->Fill(0.,0.);
         continue; // No Eta Shift requested, continue
      }
    
      
      if( ((AliConversionCuts*)fCutGammaArray->At(iCut))->GetEtaShift() == 0.0){ // Eta Shift requested but not set, get shift automatically
         ((AliConversionCuts*)fCutGammaArray->At(iCut))->GetCorrectEtaShiftFromPeriod(fV0Reader->GetPeriodName());
         ((AliConversionCuts*)fCutGammaArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once   
	 ((AliDalitzElectronCuts*)fCutElectronArray->At(iCut))->SetEtaShift( ((AliConversionCuts*)fCutGammaArray->At(iCut))->GetEtaShift() );
         hEtaShift[iCut]->Fill(0.,(((AliConversionCuts*)fCutGammaArray->At(iCut))->GetEtaShift()));
         continue;
      }
      else{
         printf(" Gamma Conversion Dalitz Task %s :: Eta Shift Manually Set to %f \n\n",
         (((AliConversionCuts*)fCutGammaArray->At(iCut))->GetCutNumber()).Data(),((AliConversionCuts*)fCutGammaArray->At(iCut))->GetEtaShift());
         ((AliConversionCuts*)fCutGammaArray->At(iCut))->DoEtaShift(kFALSE); // Eta Shift Set, make sure that it is called only once   
	 ((AliDalitzElectronCuts*)fCutElectronArray->At(iCut))->SetEtaShift( ((AliConversionCuts*)fCutGammaArray->At(iCut))->GetEtaShift() );
          hEtaShift[iCut]->Fill(0.,(((AliConversionCuts*)fCutGammaArray->At(iCut))->GetEtaShift()));
      }
   }
   
   return kTRUE;
}


void AliAnalysisTaskGammaConvDalitzV1::Terminate(const Option_t *)
{

///Grid

}
//________________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::ProcessPhotonCandidates()
{
   Int_t nV0 = 0;
   TList *GoodGammasStepOne = new TList();
   TList *GoodGammasStepTwo = new TList();
   // Loop over Photon Candidates allocated by ReaderV1
   
   for(Int_t i = 0; i < fReaderGammas->GetEntriesFast(); i++){
      AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) fReaderGammas->At(i);
      if(!PhotonCandidate) continue;
      
      fIsFromMBHeader = kTRUE;
      
      if( fMCEvent && ((AliConversionCuts*)fCutGammaArray->At(fiCut))->GetSignalRejection() != 0 ){
	
         Int_t isPosFromMBHeader
            = ((AliConversionCuts*)fCutGammaArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCStack, fInputEvent);
         if(isPosFromMBHeader == 0 && ((AliConversionCuts*)fCutGammaArray->At(fiCut))->GetSignalRejection() != 3) continue;
         
         Int_t isNegFromMBHeader
            = ((AliConversionCuts*)fCutGammaArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCStack,fInputEvent);
         if(isNegFromMBHeader == 0 && ((AliConversionCuts*)fCutGammaArray->At(fiCut))->GetSignalRejection() != 3) continue;
         
         if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
         
      }
      
      if(!((AliConversionCuts*)fCutGammaArray->At(fiCut))->PhotonIsSelected(PhotonCandidate,fESDEvent)) continue;

      if(!((AliConversionCuts*)fCutGammaArray->At(fiCut))->UseElecSharingCut() &&
         !((AliConversionCuts*)fCutGammaArray->At(fiCut))->UseToCloseV0sCut()){ // if no post reader loop is required add to events good gammas
         
         fGoodGammas->Add(PhotonCandidate);
      
	 if(fIsFromMBHeader){
	      hESDConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
	      hESDConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
	      if( fDoMesonQA ) hESDConvGammaZR[fiCut]->Fill(PhotonCandidate->GetConversionZ(),PhotonCandidate->GetConversionRadius());
	 }
	 
         if(fMCEvent){
            ProcessTruePhotonCandidates(PhotonCandidate);
         }
      }
      else if(((AliConversionCuts*)fCutGammaArray->At(fiCut))->UseElecSharingCut()){ // if Shared Electron cut is enabled, Fill array, add to step one
         ((AliConversionCuts*)fCutGammaArray->At(fiCut))->FillElectonLabelArray(PhotonCandidate,nV0);
         nV0++;
         GoodGammasStepOne->Add(PhotonCandidate);
      }
      else if(!((AliConversionCuts*)fCutGammaArray->At(fiCut))->UseElecSharingCut() &&
               ((AliConversionCuts*)fCutGammaArray->At(fiCut))->UseToCloseV0sCut()){ // shared electron is disabled, step one not needed -> step two
         GoodGammasStepTwo->Add(PhotonCandidate);
      }
   }
   
   
   if(((AliConversionCuts*)fCutGammaArray->At(fiCut))->UseElecSharingCut()){
      for(Int_t i = 0;i<GoodGammasStepOne->GetEntries();i++){
         AliAODConversionPhoton *PhotonCandidate= (AliAODConversionPhoton*) GoodGammasStepOne->At(i);
         if(!PhotonCandidate) continue;
         
         
         fIsFromMBHeader = kTRUE;
         if(fMCEvent && ((AliConversionCuts*)fCutGammaArray->At(fiCut))->GetSignalRejection() != 0){
            Int_t isPosFromMBHeader
               = ((AliConversionCuts*)fCutGammaArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCStack,fInputEvent);
            Int_t isNegFromMBHeader
               = ((AliConversionCuts*)fCutGammaArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCStack,fInputEvent);
            if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
         }
         
         
         if(!((AliConversionCuts*)fCutGammaArray->At(fiCut))->RejectSharedElectronV0s(PhotonCandidate,i,GoodGammasStepOne->GetEntries())) continue;
         if(!((AliConversionCuts*)fCutGammaArray->At(fiCut))->UseToCloseV0sCut()){ // To Colse v0s cut diabled, step two not needed
            fGoodGammas->Add(PhotonCandidate);
	    
	     if(fIsFromMBHeader){
		hESDConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
		hESDConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
		if( fDoMesonQA )hESDConvGammaZR[fiCut]->Fill(PhotonCandidate->GetConversionZ(),PhotonCandidate->GetConversionRadius());
	     }
	     
            if(fMCEvent){
               ProcessTruePhotonCandidates(PhotonCandidate);
            }
         }
         else GoodGammasStepTwo->Add(PhotonCandidate); // Close v0s cut enabled -> add to list two
      }
   }
   if(((AliConversionCuts*)fCutGammaArray->At(fiCut))->UseToCloseV0sCut()){
      for(Int_t i = 0;i<GoodGammasStepTwo->GetEntries();i++){
         AliAODConversionPhoton* PhotonCandidate = (AliAODConversionPhoton*) GoodGammasStepTwo->At(i);
         if(!PhotonCandidate) continue;
         
         if(fMCEvent && ((AliConversionCuts*)fCutGammaArray->At(fiCut))->GetSignalRejection() != 0){
            Int_t isPosFromMBHeader
               = ((AliConversionCuts*)fCutGammaArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelPositive(), fMCStack,fInputEvent);
            Int_t isNegFromMBHeader
               = ((AliConversionCuts*)fCutGammaArray->At(fiCut))->IsParticleFromBGEvent(PhotonCandidate->GetMCLabelNegative(), fMCStack,fInputEvent);
            if( (isNegFromMBHeader+isPosFromMBHeader) != 4) fIsFromMBHeader = kFALSE;
         }
         
         if(!((AliConversionCuts*)fCutGammaArray->At(fiCut))->RejectToCloseV0s(PhotonCandidate,GoodGammasStepTwo,i)) continue;
         fGoodGammas->Add(PhotonCandidate); // Add gamma to current cut TList
	 
	  if(fIsFromMBHeader){
	      hESDConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt()); // Differences to old V0Reader in p_t due to conversion KF->TLorentzVector
	      hESDConvGammaEta[fiCut]->Fill(PhotonCandidate->Eta());
	      if(fDoMesonQA)hESDConvGammaZR[fiCut]->Fill(PhotonCandidate->GetConversionZ(),PhotonCandidate->GetConversionRadius());
	  }
	 
         if(fMCEvent){
            ProcessTruePhotonCandidates(PhotonCandidate);
         }
      }
   }

   delete GoodGammasStepOne;
   GoodGammasStepOne = 0x0;
   delete GoodGammasStepTwo;
   GoodGammasStepTwo = 0x0;
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::ProcessTruePhotonCandidates(AliAODConversionPhoton *TruePhotonCandidate)
{
   // Process True Photons
   AliStack *MCStack = fMCEvent->Stack();
   TParticle *posDaughter = TruePhotonCandidate->GetPositiveMCDaughter(MCStack);
   TParticle *negDaughter = TruePhotonCandidate->GetNegativeMCDaughter(MCStack);

   if(posDaughter == NULL || negDaughter == NULL) return; // One particle does not exist
   if(posDaughter->GetMother(0) != negDaughter->GetMother(0)){  // Not Same Mother == Combinatorial Bck
      return;
   }
   
   else if (posDaughter->GetMother(0) == -1){
      return;
   }
   
   if(TMath::Abs(posDaughter->GetPdgCode())!=11 || TMath::Abs(negDaughter->GetPdgCode())!=11) return; //One Particle is not electron
   if(posDaughter->GetPdgCode()==negDaughter->GetPdgCode()) return; // Same Charge
   if(posDaughter->GetUniqueID() != 5 || negDaughter->GetUniqueID() !=5) return;// check if the daughters come from a conversion

   TParticle *Photon = TruePhotonCandidate->GetMCParticle(MCStack);
   if(Photon->GetPdgCode() != 22) return; // Mother is no Photon

   // True Photon
  
  Int_t labelGamma = TruePhotonCandidate->GetMCParticleLabel(MCStack);
  
  if( labelGamma < MCStack->GetNprimary() ){
    if( fIsFromMBHeader ){
      hESDTrueConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
    }
  }
  else {
    if( fIsFromMBHeader){
      hESDTrueSecConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
    }
  }
 
  if( IsPi0DalitzDaughter(labelGamma) == kTRUE ) {
        if( labelGamma < MCStack->GetNprimary() ) {
	  if( fIsFromMBHeader ){
	    hESDTruePi0DalitzConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
	  }
	}
	else {
	  if( fIsFromMBHeader ) {
	    hESDTruePi0DalitzSecConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
	  }
	}
  }
   
 
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::ProcessElectronCandidates(){

   Double_t magField = fInputEvent->GetMagneticField();


   if( magField  < 0.0 ){
      magField =  1.0;
   }
   else {
      magField =  -1.0;
   }


   vector<Int_t> lGoodElectronIndexPrev(0);
   vector<Int_t> lGoodPositronIndexPrev(0);
   
   
  

   for(UInt_t i = 0; i < fSelectorElectronIndex.size(); i++){
    AliESDtrack* electronCandidate = fESDEvent->GetTrack(fSelectorElectronIndex[i]);
    if(! ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->ElectronIsSelected(electronCandidate) ) continue;
	lGoodElectronIndexPrev.push_back(   fSelectorElectronIndex[i] );
	hESDDalitzElectronPt[fiCut]->Fill(electronCandidate->Pt());
        hESDDalitzElectronPhi[fiCut]->Fill(electronCandidate->Phi());
	  if( fMCEvent ) {
	  Int_t labelelectron = TMath::Abs( electronCandidate->GetLabel() );
	    if( labelelectron < fMCStack->GetNtrack() ){
	      TParticle* electron = fMCStack->Particle(labelelectron);
		if( electron->GetPdgCode() ==  11 ){
		    if( labelelectron < fMCStack->GetNprimary() ){
		      hESDTrueElectronPt[fiCut]->Fill(electronCandidate->Pt());    //primary electron
		    }
		    else{
		      hESDTrueSecElectronPt[fiCut]->Fill(electronCandidate->Pt()); //secondary electron
		    }
		    if( IsPi0DalitzDaughter(labelelectron) == kTRUE ) {
			if( labelelectron < fMCStack->GetNprimary() ) {
			  hESDTruePi0DalitzElectronPt[fiCut]->Fill(electronCandidate->Pt());
			}
			else{
			  hESDTruePi0DalitzSecElectronPt[fiCut]->Fill(electronCandidate->Pt());
			}
		    }
		}
	    }
	}
   }

   for(UInt_t i = 0; i < fSelectorPositronIndex.size(); i++){

      AliESDtrack* positronCandidate = fESDEvent->GetTrack( fSelectorPositronIndex[i] );
      if(! ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->ElectronIsSelected(positronCandidate) ) continue;
	lGoodPositronIndexPrev.push_back(   fSelectorPositronIndex[i]  );
        hESDDalitzPositronPt[fiCut]->Fill( positronCandidate->Pt() );
	hESDDalitzPositronPhi[fiCut]->Fill( positronCandidate->Phi() );
      
        if( fMCEvent ) {
	  Int_t labelpositron = TMath::Abs( positronCandidate->GetLabel() );
	  if( labelpositron < fMCStack->GetNtrack() ) {
	    TParticle* positron = fMCStack->Particle(labelpositron);
            if( positron->GetPdgCode() ==  -11 ){
	      if( labelpositron < fMCStack->GetNprimary() ){
		hESDTruePositronPt[fiCut]->Fill(positronCandidate->Pt());
	      }
	      else{
		hESDTrueSecPositronPt[fiCut]->Fill(positronCandidate->Pt());
	      }
	      if( IsPi0DalitzDaughter(labelpositron) == kTRUE ) {
		if( labelpositron < fMCStack->GetNprimary() ){
		  hESDTruePi0DalitzPositronPt[fiCut]->Fill(positronCandidate->Pt());
		}
		else{
		  hESDTruePi0DalitzSecPositronPt[fiCut]->Fill(positronCandidate->Pt());
		}
              }
            }
          }
        }
    }


   vector<Bool_t> lElectronPsiIndex(lGoodElectronIndexPrev.size(), kTRUE);
   vector<Bool_t> lPositronPsiIndex(lGoodPositronIndexPrev.size(), kTRUE);


   if( ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->DoPsiPairCut() == kTRUE ){

      for( UInt_t i = 0; i < lGoodElectronIndexPrev.size(); i++ ) {

         AliESDtrack *electronCandidate = fESDEvent->GetTrack(lGoodElectronIndexPrev[i]);

         for(UInt_t j = 0; j <  lGoodPositronIndexPrev.size(); j++){
            AliESDtrack *positronCandidate = fESDEvent->GetTrack(lGoodPositronIndexPrev[j]);
            Double_t psiPair = GetPsiPair(positronCandidate,electronCandidate);
            Double_t deltaPhi = magField * TVector2::Phi_mpi_pi( electronCandidate->GetConstrainedParam()->Phi()-positronCandidate->GetConstrainedParam()->Phi());

            if( ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->IsFromGammaConversion(psiPair,deltaPhi ) ){
               lElectronPsiIndex[i] = kFALSE;
               lPositronPsiIndex[j] = kFALSE;
            }
         }
      }
   }


   vector<Int_t> lGoodElectronIndex(0);
   vector<Int_t> lGoodPositronIndex(0);
   

   for( UInt_t i = 0; i < lGoodElectronIndexPrev.size(); i++ ) {

	if(  lElectronPsiIndex[i] == kTRUE )
	lGoodElectronIndex.push_back(   lGoodElectronIndexPrev[i]  );
   }

   for( UInt_t i = 0; i < lGoodPositronIndexPrev.size(); i++ ) {

	 if(  lPositronPsiIndex[i] == kTRUE )
	 lGoodPositronIndex.push_back(   lGoodPositronIndexPrev[i]  );
   }



    
   
   for(UInt_t i = 0; i < lGoodElectronIndex.size(); i++){

      //if( lElectronPsiIndex[i] == kFALSE ) continue;


      AliESDtrack *electronCandidate = fESDEvent->GetTrack(lGoodElectronIndex[i]);

      AliKFParticle electronCandidateKF( *electronCandidate->GetConstrainedParam(), ::kElectron );

      for(UInt_t j = 0; j < lGoodPositronIndex.size(); j++){

         //if( lPositronPsiIndex[j] == kFALSE ) continue;

         AliESDtrack *positronCandidate = fESDEvent->GetTrack(lGoodPositronIndex[j]);
         AliKFParticle positronCandidateKF( *positronCandidate->GetConstrainedParam(), ::kPositron );
				 Bool_t isPhoton    = kFALSE;
				 Bool_t isPi0Dalitz = kFALSE;
				 Bool_t isEtaDalitz = kFALSE;
				 Bool_t isJPsi      = kFALSE;

         Double_t psiPair = GetPsiPair(positronCandidate,electronCandidate);
         Double_t deltaPhi = magField * TVector2::Phi_mpi_pi( electronCandidate->GetConstrainedParam()->Phi()-positronCandidate->GetConstrainedParam()->Phi());
         

	 AliKFConversionPhoton* virtualPhoton = NULL;
         virtualPhoton = new AliKFConversionPhoton(electronCandidateKF,positronCandidateKF);


         AliKFVertex primaryVertexImproved(*fInputEvent->GetPrimaryVertex());
         primaryVertexImproved+=*virtualPhoton;
         virtualPhoton->SetProductionVertex(primaryVertexImproved);

         virtualPhoton->SetTrackLabels( lGoodPositronIndex[j], lGoodElectronIndex[i]);
	 
	  
	 if( fMCEvent ) {
	   
	  Int_t labeln=TMath::Abs(electronCandidate->GetLabel());
          Int_t labelp=TMath::Abs(positronCandidate->GetLabel());
          TParticle *fNegativeMCParticle = fMCStack->Particle(labeln);
          TParticle *fPositiveMCParticle = fMCStack->Particle(labelp);
          
	    if( fPositiveMCParticle && fNegativeMCParticle) {
               virtualPhoton->SetMCLabelPositive(labelp);
               virtualPhoton->SetMCLabelNegative(labeln);
	    }
	 }
	 
	  AliAODConversionPhoton *vphoton = new AliAODConversionPhoton(virtualPhoton); //To Apply PsiPairCut
	    
	  if ( fDoMesonQA ) {

	      if( fMCEvent ) {

	      TParticle *mcVgamma=virtualPhoton->GetMCParticle(fMCStack);

	      if(mcVgamma){
	      // Check if it is a true photon
	      if(mcVgamma->GetPdgCode() == 22){
		isPhoton = kTRUE;
	      }else if(mcVgamma->GetPdgCode() == 443){
		isJPsi = kTRUE;
	      }
	      else if( IsDalitz( mcVgamma ) ){
		if     ( mcVgamma->GetPdgCode() == 111 ) isPi0Dalitz = kTRUE;
		else if( mcVgamma->GetPdgCode() == 221 ) isEtaDalitz = kTRUE;
	      }
	      }

	      if(isPhoton){
                    hESDEposEnegTruePhotonInvMassPt[fiCut]->Fill(vphoton->GetMass(),vphoton->Pt());
                    hESDEposEnegTruePhotonPsiPairDPhi[fiCut]->Fill(deltaPhi,psiPair);
              }
	      else if(isJPsi){
                    hESDEposEnegTrueJPsiInvMassPt[fiCut]->Fill(vphoton->GetMass(),vphoton->Pt());
              }
	      else if(isPi0Dalitz){
                    hESDEposEnegTruePi0DalitzInvMassPt[fiCut]->Fill(vphoton->GetMass(),vphoton->Pt());
                    hESDEposEnegTruePi0DalitzPsiPairDPhi[fiCut]->Fill(deltaPhi,psiPair);
              }
	      else if(isEtaDalitz){
                    hESDEposEnegTrueEtaDalitzInvMassPt[fiCut]->Fill(vphoton->GetMass(),vphoton->Pt());
                    hESDEposEnegTrueEtaDalitzPsiPairDPhi[fiCut]->Fill(deltaPhi,psiPair);
              }
	    }
	  }
	 
      
	 
	 if ( fDoMesonQA ) {
	   
         hESDEposEnegPsiPairDPhi[fiCut]->Fill(deltaPhi,psiPair);   
         hESDEposEnegInvMassPt[fiCut]->Fill(vphoton->GetMass(),vphoton->Pt());
			 
	 }
				 
	 if( ! fDoChicAnalysis ) {
	
		  if (  ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->DoMassCut() == kTRUE )  {
		    
		     Double_t MassCutMax = 1000.0;
		     if(  ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->GetMassCutLowPt() >= ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->GetMassCutHighPt() ){
		       MassCutMax = ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->GetMassCutLowPt();
		     }
		     else {
		       MassCutMax = ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->GetMassCutHighPt();
		     }
		     
		     Bool_t DoMassMinCut = ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->DoMassMinCut();
		     Double_t MassMinCut = ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->GetMassMinCut();
		     
		     if( vphoton->GetMass() > MassCutMax || (  DoMassMinCut && vphoton->GetMass() < MassMinCut )  ) {
		       
		       
		       delete vphoton;
		       vphoton = 0x0;
		       delete virtualPhoton;
		       virtualPhoton = 0x0;
		       continue;
		       
		     }
		     
		  }
	 }
        
        
         fGoodVirtualGammas->Add(  vphoton );
	 delete virtualPhoton;
	 virtualPhoton=NULL;
				 
      }
   }
   
 
   //Computing mixing event
   
   if(  fDoMesonQA ) {

      for(UInt_t i = 0; i < lGoodElectronIndex.size(); i++){
         
         //if( lElectronPsiIndex[i] == kFALSE ) continue;
         
          AliESDtrack *electronCandidate1 = fESDEvent->GetTrack(lGoodElectronIndex[i]);

          AliKFParticle electronCandidate1KF( *electronCandidate1->GetConstrainedParam(), ::kElectron );
         
        
          for(UInt_t j = i+1; j < lGoodElectronIndex.size(); j++){
	    
	       //if( lElectronPsiIndex[j] == kFALSE ) continue;
	       
	       
	        AliESDtrack *electronCandidate2 = fESDEvent->GetTrack(lGoodElectronIndex[j]);

                AliKFParticle electronCandidate2KF( *electronCandidate2->GetConstrainedParam(), ::kElectron );
	       
	        AliKFConversionPhoton* virtualPhoton = new AliKFConversionPhoton(electronCandidate1KF,electronCandidate2KF);
		
		AliKFVertex primaryVertexImproved(*fInputEvent->GetPrimaryVertex());
		primaryVertexImproved+=*virtualPhoton;
		virtualPhoton->SetProductionVertex(primaryVertexImproved);

		    
		AliAODConversionPhoton *vphoton = new AliAODConversionPhoton(virtualPhoton); 
		hESDEposEnegLikeSignBackInvMassPt[fiCut]->Fill(vphoton->GetMass(),vphoton->Pt());
		delete vphoton;
		delete virtualPhoton;
		vphoton = 0x0;
		virtualPhoton = 0x0;            
	   
	  }
   }   
   
   
      for(UInt_t i = 0; i < lGoodPositronIndex.size(); i++){
     
   
         
         //if( lPositronPsiIndex[i] == kFALSE ) continue;
         
          AliESDtrack *positronCandidate1 = fESDEvent->GetTrack(lGoodPositronIndex[i]);

          AliKFParticle positronCandidate1KF( *positronCandidate1->GetConstrainedParam(), ::kPositron );
         
        
          for(UInt_t j = i+1; j < lGoodPositronIndex.size(); j++){
	    
	      // if( lPositronPsiIndex[j] == kFALSE ) continue;
	       
	        AliESDtrack *positronCandidate2 = fESDEvent->GetTrack(lGoodPositronIndex[j]);

                AliKFParticle positronCandidate2KF( *positronCandidate2->GetConstrainedParam(), ::kPositron );
	       
	        AliKFConversionPhoton* virtualPhoton = new AliKFConversionPhoton(positronCandidate1KF,positronCandidate2KF);
		AliKFVertex primaryVertexImproved(*fInputEvent->GetPrimaryVertex());
		primaryVertexImproved+=*virtualPhoton;
		virtualPhoton->SetProductionVertex(primaryVertexImproved);
	       
    		AliAODConversionPhoton *vphoton = new AliAODConversionPhoton(virtualPhoton); 
		hESDEposEnegLikeSignBackInvMassPt[fiCut]->Fill(vphoton->GetMass(),vphoton->Pt());
		
		
		delete vphoton;
		delete virtualPhoton;  
	        vphoton = 0x0;
		virtualPhoton = 0x0;
	   
	  }
      }   

   }
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::CalculatePi0DalitzCandidates(){

   // Conversion Gammas




   if( fGoodGammas->GetEntries() > 0 && fGoodVirtualGammas->GetEntries() > 0 ){

      vector<Bool_t> lGoodVirtualGamma(fGoodVirtualGammas->GetEntries(), kFALSE);
     
      for(Int_t GammaIndex=0; GammaIndex<fGoodGammas->GetEntries(); GammaIndex++){

         AliAODConversionPhoton *gamma=dynamic_cast<AliAODConversionPhoton*>(fGoodGammas->At(GammaIndex));
         if (gamma==NULL) continue;
         for(Int_t virtualGammaIndex=0;virtualGammaIndex<fGoodVirtualGammas->GetEntries();virtualGammaIndex++){

            AliAODConversionPhoton *Vgamma=dynamic_cast<AliAODConversionPhoton*>(fGoodVirtualGammas->At(virtualGammaIndex));
            if (Vgamma==NULL) continue;
            //Check for same Electron ID
            if(gamma->GetTrackLabelPositive() == Vgamma->GetTrackLabelPositive() ||
               gamma->GetTrackLabelNegative() == Vgamma->GetTrackLabelNegative() ||
               gamma->GetTrackLabelNegative() == Vgamma->GetTrackLabelPositive() ||
               gamma->GetTrackLabelPositive() == Vgamma->GetTrackLabelNegative() ) continue;

            AliAODConversionMother *pi0cand = new AliAODConversionMother(gamma,Vgamma);
            pi0cand->SetLabels(GammaIndex,virtualGammaIndex);
                      

            if( ( ((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->MesonIsSelected(pi0cand,kTRUE,((AliConversionCuts*)fCutGammaArray->At(fiCut))->GetEtaShift())) ){
	      
	      //cout<< "Meson Accepted "<<endl;
	      
	       Int_t zbin= fBGHandler[fiCut]->GetZBinIndex(fESDEvent->GetPrimaryVertex()->GetZ());
               Int_t mbin = 0;
               if( ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->UseTrackMultiplicity() ){
                  mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fNumberOfESDTracks);
               } else {
                  mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fGoodGammas->GetEntries());
               }
	      
	      	      
	      
              if(  ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->DoMassCut() == kTRUE ) {
		
		
		 if( ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->MassCut( pi0cand->Pt() , Vgamma->GetMass() ) == kTRUE ){
	  
		 hESDMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
		
		 Double_t sparesFill[4] = {pi0cand->M(),pi0cand->Pt(),(Double_t)zbin,(Double_t)mbin};
                 sESDMotherInvMassPtZM[fiCut]->Fill(sparesFill,1);
		 
		 
		    if ( fDoMesonQA ) {
		    
		       hESDMotherPhi[fiCut]->Fill(pi0cand->Phi());
		 
			if( lGoodVirtualGamma[virtualGammaIndex] == kFALSE ) {
		 
			  FillElectronQAHistos(Vgamma);
			  
			  lGoodVirtualGamma[virtualGammaIndex] = kTRUE;
		      }
		   }
		}
	      }
	      else {
		 hESDMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
		 Double_t sparesFill[4] = {pi0cand->M(),pi0cand->Pt(),(Double_t)zbin,(Double_t)mbin};
                 sESDMotherInvMassPtZM[fiCut]->Fill(sparesFill,1);
		 
		 
		   if ( fDoMesonQA ) {
		 
		    hESDMotherPhi[fiCut]->Fill(pi0cand->Phi());
		 
		    if( lGoodVirtualGamma[virtualGammaIndex] == kFALSE ) {
		      
		      FillElectronQAHistos(Vgamma);
		 
		      lGoodVirtualGamma[virtualGammaIndex] = kTRUE;
		 
		    }
		 }
	      }
	       
	       if( fDoChicAnalysis) {
	       
               hESDPi0MotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
                
               Double_t diffMass = pi0cand->M() - Vgamma->GetMass();
    
               hESDPi0MotherDiffInvMassPt[fiCut]->Fill( diffMass , pi0cand->Pt() );
	
		if( Vgamma->GetMass() > 2.5 && Vgamma->GetMass() < 3.4){
		hESDPi0MotherDiffLimInvMassPt[fiCut]->Fill( diffMass , pi0cand->Pt() );
		}
	       }
               
              if(fMCEvent){
                  ProcessTrueMesonCandidates(pi0cand,gamma,Vgamma);
              }
            }
            delete pi0cand;
            pi0cand=0x0;
         }
      }
   }
}

//________________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::CalculateBackground(){

   Int_t zbin= fBGHandler[fiCut]->GetZBinIndex(fESDEvent->GetPrimaryVertex()->GetZ());
   Int_t mbin = 0;

   Int_t method = 0;

   method = ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->GetBKGMethod();


   if(((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->UseTrackMultiplicity()){
      mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fNumberOfESDTracks);
   } else {
      mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fGoodGammas->GetEntries());
   }

   if( method == 1 || method == 2 ) {

      AliGammaConversionAODBGHandler::GammaConversionVertex *bgEventVertex = NULL;

      if( ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->UseTrackMultiplicity() ) {

         for(Int_t nEventsInBG=0;nEventsInBG<fBGHandler[fiCut]->GetNBGEvents();nEventsInBG++){

            AliGammaConversionAODVector *previousEventV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);

            if(fMoveParticleAccordingToVertex == kTRUE && method == 1){
               bgEventVertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,nEventsInBG);
            }

            for(Int_t iCurrent=0;iCurrent<fGoodVirtualGammas->GetEntries();iCurrent++){
               AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fGoodVirtualGammas->At(iCurrent));

               for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){
                  AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));

                  if(fMoveParticleAccordingToVertex == kTRUE && method == 1 ){
                     MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
                  }

                  AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
                            

                  if( ( ((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kFALSE, ((AliConversionCuts*)fCutGammaArray->At(fiCut))->GetEtaShift()))){
		      if(  ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->DoMassCut() == kTRUE ) {
			
			if( ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->MassCut( backgroundCandidate->Pt() , currentEventGoodV0.GetMass() ) == kTRUE ){
	 		  
			 // cout<<" Mass dalitz: "<<currentEventGoodV0.GetMass()<<endl;
			  hESDMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt());
			  Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
			  sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
			 }
		      }
		      else {
			hESDMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt());
			Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
			sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
		      }       
                  }
                  delete backgroundCandidate;
                  backgroundCandidate = 0x0;
               }
            }
         }
      }
      else{
         for(Int_t nEventsInBG=0;nEventsInBG <fBGHandler[fiCut]->GetNBGEvents();nEventsInBG++){
            AliGammaConversionAODVector *previousEventV0s = fBGHandler[fiCut]->GetBGGoodV0s(zbin,mbin,nEventsInBG);
            if(previousEventV0s){
               if(fMoveParticleAccordingToVertex == kTRUE && method == 1){
                  bgEventVertex = fBGHandler[fiCut]->GetBGEventVertex(zbin,mbin,nEventsInBG);
               }
               for(Int_t iCurrent=0;iCurrent<fGoodVirtualGammas->GetEntries();iCurrent++){

                  AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fGoodVirtualGammas->At(iCurrent));

                  for(UInt_t iPrevious=0;iPrevious<previousEventV0s->size();iPrevious++){

                     AliAODConversionPhoton previousGoodV0 = (AliAODConversionPhoton)(*(previousEventV0s->at(iPrevious)));

                     if(fMoveParticleAccordingToVertex == kTRUE && method ==1){
                        MoveParticleAccordingToVertex(&previousGoodV0,bgEventVertex);
                     }

                     AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
                               
                     if((((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kFALSE,((AliConversionCuts*)fCutGammaArray->At(fiCut))->GetEtaShift()))){
		       
		       if(  ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->DoMassCut() == kTRUE ) {
			
			  if( ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->MassCut( backgroundCandidate->Pt() , currentEventGoodV0.GetMass() ) == kTRUE ){
	 		
		       
			  hESDMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt());
			  Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
			  sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
			}
		       }
		       else {
			  hESDMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt());
			  Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
			  sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
		       }
                     }
                     delete backgroundCandidate;
                     backgroundCandidate = 0x0;
                  }
               }
            }
         }
      }
   }

   else if( method == 3 ){

      for(Int_t iCurrent=0;iCurrent<fGoodVirtualGammas->GetEntries();iCurrent++){

         AliAODConversionPhoton currentEventGoodV0 = *(AliAODConversionPhoton*)(fGoodVirtualGammas->At(iCurrent));

         for(Int_t iPrevious=0;iPrevious<fGammasPool[fiCut]->GetEntries();iPrevious++){

            AliAODConversionPhoton previousGoodV0 = *(AliAODConversionPhoton*)((fGammasPool[fiCut]->At(iPrevious) ));


            AliAODConversionMother *backgroundCandidate = new AliAODConversionMother(&currentEventGoodV0,&previousGoodV0);
                      

            if((((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kFALSE, ((AliConversionCuts*)fCutGammaArray->At(fiCut))->GetEtaShift()))){
	      
	       if(  ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->DoMassCut() == kTRUE ) {
		
		    if( ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->MassCut( backgroundCandidate->Pt() , currentEventGoodV0.GetMass() ) == kTRUE ){
		  
		      
		  hESDMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt());
		  Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
		  sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
		  
		  }
	       }
	       else{
		 
		 hESDMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt());
		 Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
		 sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1); 
		 
	       }
            }
            delete backgroundCandidate;
            backgroundCandidate = 0x0;
         }
      }
   }

}
//________________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::UpdateEventByEventData(){
   //see header file for documentation

   Int_t method = 0;

   method = ( (AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->GetBKGMethod();

   


   if( method == 1 ) {

      if(fGoodGammas->GetEntries() > 0 ){

         if( ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->UseTrackMultiplicity() ){
            fBGHandler[fiCut]->AddEvent(fGoodGammas,fESDEvent->GetPrimaryVertex()->GetX(),fESDEvent->GetPrimaryVertex()->GetY(),fESDEvent->GetPrimaryVertex()->GetZ(),fNumberOfESDTracks);
         }

         else{ // means we use #V0s for multiplicity
            fBGHandler[fiCut]->AddEvent(fGoodGammas,fESDEvent->GetPrimaryVertex()->GetX(),fESDEvent->GetPrimaryVertex()->GetY(),fESDEvent->GetPrimaryVertex()->GetZ(),fGoodGammas->GetEntries());
         }
      }
   }

   else if ( method == 2 ){

      if(fGoodVirtualGammas->GetEntries() > 0 ){
         if(((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->UseTrackMultiplicity()){
            fBGHandler[fiCut]->AddEvent(fGoodVirtualGammas,fESDEvent->GetPrimaryVertex()->GetX(),fESDEvent->GetPrimaryVertex()->GetY(),fESDEvent->GetPrimaryVertex()->GetZ(),fNumberOfESDTracks);
         }
         else{ // means we use #V0s for multiplicity
            fBGHandler[fiCut]->AddEvent(fGoodVirtualGammas,fESDEvent->GetPrimaryVertex()->GetX(),fESDEvent->GetPrimaryVertex()->GetY(),fESDEvent->GetPrimaryVertex()->GetZ(),fGoodVirtualGammas->GetEntries());
         }
      }
   }
   else if ( method  == 3 ) {



      for(Int_t index = 0; index < fGoodGammas->GetEntries(); index++){


         if ( fGammasPool[fiCut]->GetEntries() > ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->NumberOfRotationEvents() ){
            fGammasPool[fiCut]->RemoveLast();
         }
         fGammasPool[fiCut]->AddFirst( new AliAODConversionPhoton(*(AliAODConversionPhoton*)(fGoodGammas->At(index)) ) );

      }
   }
}
//______________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::ProcessTrueMesonCandidates(AliAODConversionMother *Pi0Candidate, AliAODConversionPhoton *TrueGammaCandidate, AliAODConversionPhoton *TrueVirtualGammaCandidate)
{

   // Process True Mesons

   AliStack *MCStack = fMCEvent->Stack();

   if(	TrueGammaCandidate->GetV0Index()<fESDEvent->GetNumberOfV0s()	){
     
     
     //cout<<"Entro True Meson"<<endl;


      Bool_t isTruePi0 = kFALSE;
      Bool_t isTrueEta = kFALSE;
      Bool_t massCutAccept = kFALSE;
      //Bool_t isTrueChiC = kFALSE;
      Int_t gammaMCLabel = TrueGammaCandidate->GetMCParticleLabel(MCStack);
      Int_t gammaMotherLabel = -1;


      if(  ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->DoMassCut() == kTRUE ) {
	   
	    if( ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->MassCut( Pi0Candidate->Pt() , TrueVirtualGammaCandidate->GetMass() ) == kTRUE ){
	      
	       massCutAccept = kTRUE;
	    }
      }
      else {
	      massCutAccept  = kTRUE;
      }
      
      
      


      if(gammaMCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother


         // Daughters Gamma 0
         TParticle * negativeMC = (TParticle*)TrueGammaCandidate->GetNegativeMCDaughter(MCStack);
         TParticle * positiveMC = (TParticle*)TrueGammaCandidate->GetPositiveMCDaughter(MCStack);
         TParticle * gammaMC = (TParticle*)MCStack->Particle(gammaMCLabel);


         if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){  // Electrons ...

            if(negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() ==5){ // ... From Conversion ...

               if(gammaMC->GetPdgCode() == 22){ // ... with Gamma Mother
                  gammaMotherLabel=gammaMC->GetFirstMother();
               }
            }
         }
      }


      Int_t virtualGammaMCLabel = TrueVirtualGammaCandidate->GetMCParticleLabel(MCStack);
      Int_t virtualGammaMotherLabel = -1;
      Int_t virtualGamma = 1;
      Int_t virtualGammaGrandMotherLabel =-1;


      if(virtualGammaMCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
         // Daughters Gamma 1
         TParticle * negativeMC = (TParticle*)TrueVirtualGammaCandidate->GetNegativeMCDaughter(MCStack);
         TParticle * positiveMC = (TParticle*)TrueVirtualGammaCandidate->GetPositiveMCDaughter(MCStack);
         TParticle * virtualGammaMotherMC = (TParticle*)MCStack->Particle(virtualGammaMCLabel);

         if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){  // Electrons ...

            if( virtualGammaMotherMC->GetPdgCode() != 22 ){
               virtualGammaMotherLabel=virtualGammaMCLabel;
	    if(virtualGammaMotherMC->GetPdgCode() == 443){
	       virtualGammaGrandMotherLabel=virtualGammaMotherMC->GetFirstMother();
	    }
           }

            else if(negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() ==5){ // ... From Conversion ...
               virtualGammaMotherLabel=virtualGammaMotherMC->GetFirstMother();
               virtualGamma = 0; //no virtual gamma
            }
         }
      }


      if(gammaMotherLabel >= 0 && ( gammaMotherLabel == virtualGammaMotherLabel) ){

         if(((TParticle*)MCStack->Particle(virtualGammaMotherLabel))->GetPdgCode() == 111){
            isTruePi0=kTRUE;
         }

         if(((TParticle*)MCStack->Particle(virtualGammaMotherLabel))->GetPdgCode() == 221){
            isTrueEta=kTRUE;
         }
         
	 
	}

       if( fDoChicAnalysis) {
	if(gammaMotherLabel>=0 && ( gammaMotherLabel == virtualGammaGrandMotherLabel) ){
	  if(((TParticle*)MCStack->Particle(virtualGammaGrandMotherLabel))->GetPdgCode() == 445 ||
	    ((TParticle*)MCStack->Particle(virtualGammaGrandMotherLabel))->GetPdgCode() == 10443 ||
	    ((TParticle*)MCStack->Particle(virtualGammaGrandMotherLabel))->GetPdgCode() == 20443 ){
	      //isTrueChiC=kTRUE;
	      hESDTrueMotherChiCInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
	      hESDTrueMotherChiCDiffInvMassPt[fiCut]->Fill(Pi0Candidate->M()-TrueVirtualGammaCandidate->GetMass(),Pi0Candidate->Pt());
	     }
	  }  
       }

      if( ( isTruePi0 || isTrueEta) && massCutAccept ){ // True Pion or Eta
	
         if ( virtualGamma == 1 ) { //True Dalitz	

            Float_t weighted= 1;
            
            if( ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->DoWeights() ) { 
                if(((AliConversionCuts*)fCutGammaArray->At(fiCut))->IsParticleFromBGEvent(gammaMotherLabel, fMCStack,fInputEvent)){
                    if (((TParticle*)MCStack->Particle(gammaMotherLabel))->Pt()>0.005){
                        weighted= ((AliConversionCuts*)fCutGammaArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),gammaMotherLabel,fMCStack,fInputEvent);
                    }
                }
            }

            hESDTrueMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);
	    hESDTrueMotherDalitzInvMassPt[fiCut]->Fill( TrueVirtualGammaCandidate->GetMass(),Pi0Candidate->Pt(),weighted);

            if(gammaMotherLabel < MCStack->GetNprimary()){ // Only primary pi0 for efficiency calculation
	      
	       
               hESDTruePrimaryMotherInvMassPt[fiCut]->Fill( Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);
	       hESDTruePrimaryMotherW0WeightingInvMassPt[fiCut]->Fill( Pi0Candidate->M(), Pi0Candidate->Pt() );
	       
               hESDTruePrimaryMotherInvMassMCPt[fiCut]->Fill(Pi0Candidate->M(),((TParticle*)MCStack->Particle(virtualGammaMotherLabel))->Pt(),weighted);
               if(isTruePi0){ // Only primaries for unfolding
                  hESDTruePrimaryPi0DalitzESDPtMCPt[fiCut]->Fill(Pi0Candidate->Pt(),((TParticle*)MCStack->Particle(virtualGammaMotherLabel))->Pt(),weighted);
               }
            }
            else { // Secondary Meson
              Float_t weightedSec= 1;

              if( ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->DoWeights() ) { 
                  Int_t secMotherLabel = ((TParticle*)MCStack->Particle(gammaMotherLabel))->GetMother(0);
                  if(((AliConversionCuts*)fCutGammaArray->At(fiCut))->IsParticleFromBGEvent(gammaMotherLabel, fMCStack, fInputEvent) && MCStack->Particle(gammaMotherLabel)->GetPdgCode()==310){
                    weightedSec= ((AliConversionCuts*)fCutGammaArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),secMotherLabel, fMCStack, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
                  }
              }

              hESDTrueSecondaryMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec); 
            }
         }

         
         else if ( virtualGamma == 0 ){

              Float_t weighted= 1;

             if( ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->DoWeights() ) {
              if(((AliConversionCuts*)fCutGammaArray->At(fiCut))->IsParticleFromBGEvent(gammaMotherLabel, fMCStack,fInputEvent)){
                if (((TParticle*)MCStack->Particle(gammaMotherLabel))->Pt()>0.005){
                    weighted= ((AliConversionCuts*)fCutGammaArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),gammaMotherLabel,fMCStack,fInputEvent);
                }
              }
            }

            hESDTrueMotherPi0GGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted); // Pi0 from GG

            if( gammaMotherLabel < MCStack->GetNprimary() ){
                hESDTruePrimaryMotherPi0GGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weighted);
            }
            else {
            
            Float_t weightedSec= 1;
           
            if( ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->DoWeights() ) { 
                Int_t secMotherLabel = ((TParticle*)MCStack->Particle(gammaMotherLabel))->GetMother(0);
                   if(((AliConversionCuts*)fCutGammaArray->At(fiCut))->IsParticleFromBGEvent(gammaMotherLabel, fMCStack, fInputEvent) && MCStack->Particle(gammaMotherLabel)->GetPdgCode()==310){
                            weightedSec= ((AliConversionCuts*)fCutGammaArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),secMotherLabel, fMCStack, fInputEvent)/2.; //invariant mass is additive thus the weight for the daughters has to be devide by two for the K0s at a certain pt
                   }
            }
                    hESDTrueSecondaryMotherPi0GGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt(),weightedSec);
            }
         }
      }

      if(!isTruePi0 && !isTrueEta && massCutAccept ){ // Background
         if(gammaMotherLabel>-1 && virtualGammaMotherLabel>-1 && virtualGamma == 0){ // Both Tracks are Photons and have a mother but not Pi0 or Eta
            hESDTrueBckGGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
         } else { // No photon or without mother
            hESDTrueBckContInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
         }
      }
   }
}


//________________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::MoveParticleAccordingToVertex(AliAODConversionPhoton* particle,const AliGammaConversionAODBGHandler::GammaConversionVertex *vertex){
   //see header file for documentation

   Double_t dx = vertex->fX - fESDEvent->GetPrimaryVertex()->GetX();
   Double_t dy = vertex->fY - fESDEvent->GetPrimaryVertex()->GetY();
   Double_t dz = vertex->fZ - fESDEvent->GetPrimaryVertex()->GetZ();

   Double_t movedPlace[3] = {particle->GetConversionX() - dx,particle->GetConversionY() - dy,particle->GetConversionZ() - dz};
   particle->SetConversionPoint(movedPlace);
}


//________________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::CountESDTracks(){

   // Using standard function for setting Cuts
   Bool_t selectPrimaries=kTRUE;
   AliESDtrackCuts *EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(selectPrimaries);
   EsdTrackCuts->SetMaxDCAToVertexZ(2);
   EsdTrackCuts->SetEtaRange(-0.8, 0.8);
   EsdTrackCuts->SetPtRange(0.15);

   fNumberOfESDTracks = 0;
   for(Int_t iTracks = 0; iTracks < fESDEvent->GetNumberOfTracks(); iTracks++){
      AliESDtrack* curTrack = fESDEvent->GetTrack(iTracks);
      if(!curTrack) continue;
      if(EsdTrackCuts->AcceptTrack(curTrack) ) fNumberOfESDTracks++;
   }
   delete EsdTrackCuts;
   EsdTrackCuts=0x0;

   return;
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::ProcessMCParticles()
{

   // Loop over all primary MC particle

   for(Int_t i = 0; i < fMCStack->GetNprimary(); i++) {
     
     
      TParticle* particle = (TParticle *)fMCStack->Particle(i);
      if (!particle) continue;
      
      
      Bool_t mcIsFromMB = kTRUE;
      Int_t isMCFromMBHeader = -1;
           
      if(((AliConversionCuts*)fCutGammaArray->At(fiCut))->GetSignalRejection() != 0) {
         isMCFromMBHeader
            = ((AliConversionCuts*)fCutGammaArray->At(fiCut))->IsParticleFromBGEvent(i,fMCStack,fInputEvent);
         if(isMCFromMBHeader == 0 && ((AliConversionCuts*)fCutGammaArray->At(fiCut))->GetSignalRejection() != 3) continue;
         if(isMCFromMBHeader != 2) mcIsFromMB = kFALSE;
      }    
    
      if(((AliConversionCuts*)fCutGammaArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCStack,kFALSE)){
	  hMCAllGammaPt[fiCut]->Fill(particle->Pt()); // All MC Gamma
      }
      
      if(((AliConversionCuts*)fCutGammaArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCStack,kTRUE)){
         hMCConvGammaPt[fiCut]->Fill(particle->Pt());
         if(mcIsFromMB){
            hMCConvGammaRSPt[fiCut]->Fill(particle->Pt());
         }
      } // Converted MC Gamma
	
      if(((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->ElectronIsSelectedMC(i,fMCStack)){
			if( particle->GetPdgCode() == -11)hMCAllPositronsPt[fiCut]->Fill(particle->Pt()); // All positrons
			if( particle->GetPdgCode() ==  11)hMCAllElectronsPt[fiCut]->Fill(particle->Pt()); // All electrons
      }
	
      if(((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->MesonIsSelectedMC( particle,fMCStack,((AliConversionCuts*)fCutGammaArray->At(fiCut))->GetEtaShift() ) ){
	
			Float_t weighted= 1;
	            
                    if( ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->DoWeights() ) { 
			if(((AliConversionCuts*)fCutGammaArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack, fInputEvent)){
			  if (particle->Pt()>0.005){
			    weighted= ((AliConversionCuts*)fCutGammaArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),i, fMCStack,fInputEvent);
			  }
			}
	            }

			if(particle->GetPdgCode() == 111)hMCPi0GGPt[fiCut]->Fill( particle->Pt() , weighted); // All MC Pi0 GG decay
			if(particle->GetPdgCode() == 221)hMCEtaGGPt[fiCut]->Fill( particle->Pt() , weighted); // All MC Eta GG decay
      }
      
      
      Int_t labelgamma 	  = -1;
      Int_t labelelectron = -1;
      Int_t labelpositron = -1;


      if( ((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->MesonIsSelectedMCDalitz(particle,fMCStack,labelelectron,labelpositron,labelgamma,((AliConversionCuts*)fCutGammaArray->At(fiCut))->GetEtaShift())  )
      {
	
	
			Float_t weighted= 1;
	           if( ((AliDalitzElectronCuts*) fCutElectronArray->At(fiCut))->DoWeights() ) { 
			if(((AliConversionCuts*)fCutGammaArray->At(fiCut))->IsParticleFromBGEvent(i, fMCStack,fInputEvent)){
			  if (particle->Pt()>0.005){
			    weighted= ((AliConversionCuts*)fCutGammaArray->At(fiCut))->GetWeightForMeson(fV0Reader->GetPeriodName(),i, fMCStack,fInputEvent);
			  }
			}
	           }
			if(particle->GetPdgCode() == 111)hMCPi0Pt[fiCut]->Fill(particle->Pt(), weighted); // All MC Pi0
			if(particle->GetPdgCode() == 221)hMCEtaPt[fiCut]->Fill(particle->Pt(), weighted); // All MC Eta
			
			// Check the acceptance for gamma and electrons
			
					 
				TParticle *gamma    = fMCStack->Particle(labelgamma);
				TParticle *electron = fMCStack->Particle(labelelectron);
				TParticle *positron = fMCStack->Particle(labelpositron);
				
           
				if(((AliConversionCuts*)fCutGammaArray->At(fiCut))->PhotonIsSelectedMC(gamma,fMCStack,kFALSE) &&
				   ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->ElectronIsSelectedMC(labelelectron,fMCStack) &&
				   ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->ElectronIsSelectedMC(labelpositron,fMCStack) ) {
				   
					if(particle->GetPdgCode() == 111){ 
					  
						hMCPi0InAccPt[fiCut]->Fill(particle->Pt() , weighted); // MC Pi0Dalitz with gamma and e+e- in acc
						hMCPi0DalitzGammaPt[fiCut]->Fill( gamma->Pt() ,weighted );
						hMCPi0DalitzPositronPt[fiCut]->Fill( positron->Pt(),weighted );
						hMCPi0DalitzElectronPt[fiCut]->Fill( electron->Pt(),weighted );
						
					}
					if(particle->GetPdgCode() == 221)hMCEtaInAccPt[fiCut]->Fill(particle->Pt(), weighted ); // MC EtaDalitz with gamma and e+e- in acc
				}
		
			
	}
		Int_t labelgammaChiC=-1;
		Int_t labelpositronChiC=-1;
		Int_t labelelectronChiC=-1;

		if(((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->MesonIsSelectedMCChiC(particle,fMCStack,labelelectronChiC,labelpositronChiC,labelgammaChiC,((AliConversionCuts*)fCutGammaArray->At(fiCut))->GetEtaShift())){
		  
			hMCChiCPt[fiCut]->Fill(particle->Pt()); // All MC ChiC
			TParticle * gammaChiC  =fMCStack->Particle(labelgammaChiC);

			if( ((AliConversionCuts*)fCutGammaArray->At(fiCut))->PhotonIsSelectedMC( gammaChiC,fMCStack,kFALSE) &&
			    ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->ElectronIsSelectedMC(labelelectronChiC,fMCStack) && 
			    ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->ElectronIsSelectedMC(labelpositronChiC,fMCStack) ){
				hMCChiCInAccPt[fiCut]->Fill(particle->Pt()); // All MC ChiC
			} 
		}
	}
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskGammaConvDalitzV1::IsDalitz(TParticle *fMCMother) const
{

	if( fMCMother->GetNDaughters() != 3 ) return kFALSE;
	if( fMCMother->GetPdgCode() != 111 && fMCMother->GetPdgCode() != 221 ) return kFALSE;
	
	
	TParticle *positron = 0x0;
	TParticle *electron = 0x0;
	TParticle *gamma    = 0x0;
	
	for(Int_t index= fMCMother->GetFirstDaughter();index<= fMCMother->GetLastDaughter();index++){				

		TParticle* temp = (TParticle*)fMCStack->Particle( index );
		
		switch( temp->GetPdgCode() ) {
		case ::kPositron:
			positron =  temp;
			break;
		case ::kElectron:
			electron =  temp;
			break;
		case ::kGamma:
			gamma    =  temp;
			break;
		}
	}
  
	if( positron && electron && gamma) return kTRUE;
	
	return kFALSE;
}
//_____________________________________________________________________________________
Bool_t AliAnalysisTaskGammaConvDalitzV1::IsPi0DalitzDaughter( Int_t label ) const
{
//
// Returns true if the particle comes from Pi0 -> e+ e- gamma
//
        
        Int_t motherLabel = fMCStack->Particle( label )->GetMother(0);
        
        if( motherLabel < 0 || motherLabel >= fMCStack->GetNtrack() ) return kFALSE;
        
        TParticle* mother = fMCStack->Particle( motherLabel );
	
	if( mother->GetPdgCode() != 111 ) return kFALSE;
	
	if( IsDalitz( mother ) ) return kTRUE;
	
	
	return kFALSE;
        
       
}

void AliAnalysisTaskGammaConvDalitzV1::FillElectronQAHistos(AliAODConversionPhoton *Vgamma) const
{
  

	      AliESDtrack *positronVgamma = 0;
	      AliESDtrack *electronVgamma = 0;
	      
	      Double_t clsToFPos = -1.0;
	      Double_t clsToFNeg = -1.0;
	     
	      Double_t NumClsITSPos = -1.0;
	      Double_t NumClsITSNeg = -1.0;
	      Double_t NumClsTPCPos = -1.0;
	      Double_t NumClsTPCNeg = -1.0;
	      	      
	      Float_t dcaToVertexXYPos = -1.0;
	      Float_t dcaToVertexZPos  = -1.0;
	      Float_t dcaToVertexXYNeg = -1.0;
	      Float_t dcaToVertexZNeg  = -1.0;
	      
	      Double_t nSigmaPosTPC = -999.;
	      Double_t nSigmaNegTPC = -999.;
	      
	      positronVgamma = fESDEvent->GetTrack( Vgamma->GetTrackLabelPositive() );
	      electronVgamma = fESDEvent->GetTrack( Vgamma->GetTrackLabelNegative() );
	      clsToFPos = ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->GetNFindableClustersTPC(positronVgamma);
	      clsToFNeg = ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->GetNFindableClustersTPC(electronVgamma);
	      
	      nSigmaPosTPC = ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->GetPIDResponse()->NumberOfSigmasTPC(positronVgamma, AliPID::kElectron) ;
	      nSigmaNegTPC = ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->GetPIDResponse()->NumberOfSigmasTPC(electronVgamma, AliPID::kElectron) ;
	       
	      
	      
	      NumClsITSPos =  positronVgamma->GetNcls(0); //Get number of ITS clusters
	      NumClsITSNeg =  electronVgamma->GetNcls(0);
	      NumClsTPCPos =  positronVgamma->GetNcls(1);  //Get number of TPC clusters
	      NumClsTPCNeg =  electronVgamma->GetNcls(1);
	      
	      
	      Float_t bPos[2];
	      Float_t bCovPos[3];
	      positronVgamma->GetImpactParameters(bPos,bCovPos);
	      
	      if (bCovPos[0]<=0 || bCovPos[2]<=0) {
		  AliDebug(1, "Estimated b resolution lower or equal zero!");
		  bCovPos[0]=0; bCovPos[2]=0;
	       }
	      
	       Float_t bNeg[2];
	       Float_t bCovNeg[3];
	       positronVgamma->GetImpactParameters(bNeg,bCovNeg);
	       
	       if (bCovNeg[0]<=0 || bCovNeg[2]<=0) {
		  AliDebug(1, "Estimated b resolution lower or equal zero!");
		  bCovNeg[0]=0; bCovNeg[2]=0;
		}
	      
	       dcaToVertexXYPos = bPos[0];
	       dcaToVertexZPos  = bPos[1];
	       dcaToVertexXYNeg = bNeg[0];
	       dcaToVertexZNeg  = bNeg[1];
	      
	       hESDDalitzElectronAfterPt[fiCut]->Fill(  electronVgamma->Pt()  );
	       hESDDalitzPositronAfterPt[fiCut]->Fill(  positronVgamma->Pt()  );
			  
	       hESDDalitzElectronAfterEta[fiCut]->Fill( electronVgamma->Eta() );
	       hESDDalitzPositronAfterEta[fiCut]->Fill( positronVgamma->Eta() );
		 
	       hESDDalitzElectronAfterPhi[fiCut]->Fill( electronVgamma->Phi() );
	       hESDDalitzPositronAfterPhi[fiCut]->Fill( positronVgamma->Phi() );
		      
	       hESDDalitzElectronAfterNFindClsTPC[fiCut]->Fill(clsToFNeg,electronVgamma->Pt());
	       hESDDalitzPositronAfterNFindClsTPC[fiCut]->Fill(clsToFPos,positronVgamma->Pt());
	       
	       hESDDalitzElectronAfterNClsTPC[fiCut]->Fill( NumClsTPCNeg,electronVgamma->Pt());
	       hESDDalitzPositronAfterNClsTPC[fiCut]->Fill( NumClsTPCPos,positronVgamma->Pt());
	       
	       hESDDalitzElectronAfterNClsITS[fiCut]->Fill( NumClsITSNeg);
	       hESDDalitzPositronAfterNClsITS[fiCut]->Fill( NumClsITSPos);
	       
	       hESDDalitzPosEleAfterDCAxy[fiCut]->Fill(  dcaToVertexXYNeg, electronVgamma->Pt() );
	       hESDDalitzPosEleAfterDCAz[fiCut]->Fill(   dcaToVertexZNeg,  electronVgamma->Pt() );
	       hESDDalitzPosEleAfterDCAxy[fiCut]->Fill(  dcaToVertexXYPos, positronVgamma->Pt() );
	       hESDDalitzPosEleAfterDCAz[fiCut]->Fill(   dcaToVertexZPos,  positronVgamma->Pt() );
		      
	       hESDDalitzElectronAfterTPCdEdxVsP[fiCut]->Fill( electronVgamma->P(),nSigmaNegTPC);
	       hESDDalitzPositronAfterTPCdEdxVsP[fiCut]->Fill( positronVgamma->P(), nSigmaPosTPC);
	       
	       hESDDalitzElectronAfterTPCdEdxVsEta[fiCut]->Fill( electronVgamma->Eta(),nSigmaNegTPC);
               hESDDalitzPositronAfterTPCdEdxVsEta[fiCut]->Fill( positronVgamma->Eta(),nSigmaPosTPC);
	       
               hESDDalitzElectronAfterTPCdEdxVsPhi[fiCut]->Fill( electronVgamma->Phi(),nSigmaNegTPC);
               hESDDalitzPositronAfterTPCdEdxVsPhi[fiCut]->Fill( positronVgamma->Phi(),nSigmaPosTPC);
	       
	       		      
	       hESDDalitzElectronAfterTPCdEdxSignalVsP[fiCut]->Fill( electronVgamma->P(), TMath::Abs(electronVgamma->GetTPCsignal()));
	       hESDDalitzPositronAfterTPCdEdxSignalVsP[fiCut]->Fill( positronVgamma->P(), TMath::Abs(positronVgamma->GetTPCsignal()));
		    
	       
  
}


//_____________________________________________________________________________
Double_t AliAnalysisTaskGammaConvDalitzV1::GetPsiPair( const AliESDtrack *trackPos, const AliESDtrack *trackNeg ) const
{
   //
   // This angle is a measure for the contribution of the opening in polar
   // direction ?0 to the opening angle ? Pair
   //
   // Ref. Measurement of photons via conversion pairs with the PHENIX experiment at RHIC
   //      Mas   ter Thesis. Thorsten Dahms. 2005
   // https://twiki.cern.ch/twiki/pub/ALICE/GammaPhysicsPublications/tdahms_thesis.pdf
   //
   Double_t momPos[3];
   Double_t momNeg[3];
   if( trackPos->GetConstrainedPxPyPz(momPos) == 0 ) trackPos->GetPxPyPz( momPos );
   if( trackNeg->GetConstrainedPxPyPz(momNeg) == 0 ) trackNeg->GetPxPyPz( momNeg );

   TVector3 posDaughter;
   TVector3 negDaughter;

   posDaughter.SetXYZ( momPos[0], momPos[1], momPos[2] );
   negDaughter.SetXYZ( momNeg[0], momNeg[1], momNeg[2] );

   Double_t deltaTheta = negDaughter.Theta() - posDaughter.Theta();
   Double_t openingAngle =  posDaughter.Angle( negDaughter );  //TMath::ACos( posDaughter.Dot(negDaughter)/(negDaughter.Mag()*posDaughter.Mag()) );

   if( openingAngle < 1e-20 ) return 0.;

   Double_t psiAngle = TMath::ASin( deltaTheta/openingAngle );

   return psiAngle;
}
