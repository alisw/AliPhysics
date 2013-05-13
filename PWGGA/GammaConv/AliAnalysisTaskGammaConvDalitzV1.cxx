/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Pedro González, Pedro Ladrón de Guevara, Ernesto López Torres, *
 *         Eulogio Serradilla                                             *
 * Version 2                                                           *
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
   fTrueList(NULL),
   fMCList(NULL),
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
   hESDDalitzElectronPt(NULL),
   hESDDalitzPositronPt(NULL),
   hESDEposEnegPsiPairDPhi(NULL),
   hESDEposEnegInvMassPt(NULL),
   hESDEposEnegLikeSignBackInvMassPt(NULL),
   hESDMotherInvMassPt(NULL),
   hESDPi0MotherInvMassPt(NULL),
   hESDPi0MotherDiffInvMassPt(NULL),
   sESDMotherInvMassPtZM(NULL),
   hESDMotherBackInvMassPt(NULL),
   sESDMotherBackInvMassPtZM(NULL),
   hMCAllGammaPt(NULL),
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
   hESDEposEnegTruePi0DalitzInvMassPt(NULL),
   hESDEposEnegTrueEtaDalitzInvMassPt(NULL),
   hESDEposEnegTruePhotonInvMassPt(NULL),
   hESDTrueMotherInvMassPt(NULL),
   hESDTrueMotherPi0GGInvMassPt(NULL),
   hESDTruePrimaryMotherInvMassMCPt(NULL),
   hESDTruePrimaryPi0DalitzESDPtMCPt(NULL),
   hESDTrueSecondaryMotherInvMassPt(NULL),
   hESDTrueSecondaryMotherFromK0sInvMassPt(NULL),
   hESDTrueBckGGInvMassPt(NULL),
   hESDTrueBckContInvMassPt(NULL),
   hESDTrueMotherGGInvMassPt(NULL),
   hESDTrueConvGammaPt(NULL),
   hESDTruePositronPt(NULL),
   hESDTrueElectronPt(NULL),
   hESDTruePi0DalitzConvGammaPt(NULL),
   hESDTruePi0DalitzPositronPt(NULL),
   hESDTruePi0DalitzElectronPt(NULL),
   hNEvents(NULL),
   hNGoodESDTracks(NULL),
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
   fDoMesonAnalysis(kTRUE)
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
   fTrueList(NULL),
   fMCList(NULL),
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
   hESDDalitzElectronPt(NULL),
   hESDDalitzPositronPt(NULL),
   hESDEposEnegPsiPairDPhi(NULL),
   hESDEposEnegInvMassPt(NULL),
   hESDEposEnegLikeSignBackInvMassPt(NULL),
   hESDMotherInvMassPt(NULL),
   hESDPi0MotherInvMassPt(NULL),
   hESDPi0MotherDiffInvMassPt(NULL),
   sESDMotherInvMassPtZM(NULL),
   hESDMotherBackInvMassPt(NULL),
   sESDMotherBackInvMassPtZM(NULL),
   hMCAllGammaPt(NULL),
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
   hESDEposEnegTruePi0DalitzInvMassPt(NULL),
   hESDEposEnegTrueEtaDalitzInvMassPt(NULL),
   hESDEposEnegTruePhotonInvMassPt(NULL),
   hESDTrueMotherInvMassPt(NULL),
   hESDTrueMotherPi0GGInvMassPt(NULL),
   hESDTruePrimaryMotherInvMassMCPt(NULL),
   hESDTruePrimaryPi0DalitzESDPtMCPt(NULL),
   hESDTrueSecondaryMotherInvMassPt(NULL),
   hESDTrueSecondaryMotherFromK0sInvMassPt(NULL),
   hESDTrueBckGGInvMassPt(NULL),
   hESDTrueBckContInvMassPt(NULL),
   hESDTrueMotherGGInvMassPt(NULL),
   hESDTrueConvGammaPt(NULL),
   hESDTruePositronPt(NULL),
   hESDTrueElectronPt(NULL),
   hESDTruePi0DalitzConvGammaPt(NULL),
   hESDTruePi0DalitzPositronPt(NULL),
   hESDTruePi0DalitzElectronPt(NULL),
   hNEvents(NULL),
   hNGoodESDTracks(NULL),
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
   fDoMesonAnalysis(kTRUE)
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

   Double_t *zBinLimitsArray= new Double_t[9];
   zBinLimitsArray[0] = -50.00;
   zBinLimitsArray[1] = -3.375;
   zBinLimitsArray[2] = -1.605;
   zBinLimitsArray[3] = -0.225;
   zBinLimitsArray[4] = 1.065;
   zBinLimitsArray[5] = 2.445;
   zBinLimitsArray[6] = 4.245;
   zBinLimitsArray[7] = 50.00;
   zBinLimitsArray[8] = 1000.00;

   Double_t *multiplicityBinLimitsArrayTracks= new Double_t[6];
   multiplicityBinLimitsArrayTracks[0] = 0;
   multiplicityBinLimitsArrayTracks[1] = 8.5;
   multiplicityBinLimitsArrayTracks[2] = 16.5;
   multiplicityBinLimitsArrayTracks[3] = 27.5;
   multiplicityBinLimitsArrayTracks[4] = 41.5;
   multiplicityBinLimitsArrayTracks[5] = 200.;

   if(fIsHeavyIon){
      multiplicityBinLimitsArrayTracks[0] = 0;
      multiplicityBinLimitsArrayTracks[1] = 200.;
      multiplicityBinLimitsArrayTracks[2] = 500.;
      multiplicityBinLimitsArrayTracks[3] = 1000.;
      multiplicityBinLimitsArrayTracks[4] = 1500.;
      multiplicityBinLimitsArrayTracks[5] = 5000.;
   }

   Double_t *multiplicityBinLimitsArrayV0s= new Double_t[5];

   multiplicityBinLimitsArrayV0s[0] = 2;
   multiplicityBinLimitsArrayV0s[1] = 3;
   multiplicityBinLimitsArrayV0s[2] = 4;
   multiplicityBinLimitsArrayV0s[3] = 5;
   multiplicityBinLimitsArrayV0s[4] = 9999;

   if(fIsHeavyIon){
      multiplicityBinLimitsArrayV0s[0] = 2;
      multiplicityBinLimitsArrayV0s[1] = 10;
      multiplicityBinLimitsArrayV0s[2] = 30;
      multiplicityBinLimitsArrayV0s[3] = 50;
      multiplicityBinLimitsArrayV0s[4] = 9999;
   }

   const Int_t nDim = 4;
   Int_t    nBins[nDim] = {1000,250,8,5};
   Double_t xMin[nDim]  = {0,0, 0,0};
   Double_t xMax[nDim]  = {1,25,8,5};

   sESDMotherInvMassPtZM     = new THnSparseF*[fnCuts];
   sESDMotherBackInvMassPtZM = new THnSparseF*[fnCuts];


   fBGHandler = new AliGammaConversionAODBGHandler*[fnCuts];

   for(Int_t iCut = 0; iCut<fnCuts;iCut++){


      TString cutstringElectron     =   ((AliDalitzElectronCuts*)fCutElectronArray->At(iCut))->GetCutNumber();
      TString cutstringMeson        =   ((AliConversionMesonCuts*)fCutMesonArray->At(iCut))->GetCutNumber();
      TString cutstringGamma        =   ((AliConversionCuts*)fCutGammaArray->At(iCut))->GetCutNumber();


      fBackList[iCut] = new TList();
      fBackList[iCut]->SetName(Form("%s_%s_%s Back histograms",cutstringGamma.Data(),cutstringElectron.Data(),cutstringMeson.Data()));
      fBackList[iCut]->SetOwner(kTRUE);
      fCutFolder[iCut]->Add(fBackList[iCut]);

      sESDMotherInvMassPtZM[iCut] = new THnSparseF("Back_Mother_InvMass_Pt_z_m","Back_Mother_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
      sESDMotherInvMassPtZM[iCut]->Sumw2();
      fBackList[iCut]->Add(sESDMotherInvMassPtZM[iCut]);
      sESDMotherBackInvMassPtZM[iCut] = new THnSparseF("Back_Back_InvMass_Pt_z_m","Back_Back_InvMass_Pt_z_m",nDim,nBins,xMin,xMax);
      sESDMotherBackInvMassPtZM[iCut]->Sumw2();
      fBackList[iCut]->Add(sESDMotherBackInvMassPtZM[iCut]);

      if(((AliDalitzElectronCuts*)fCutElectronArray->At(iCut))->UseTrackMultiplicity()) {
         fBGHandler[iCut] = new AliGammaConversionAODBGHandler(9,6,((AliDalitzElectronCuts*)fCutElectronArray->At(iCut))->NumberOfRotationEvents());
         fBGHandler[iCut]->Initialize(zBinLimitsArray, multiplicityBinLimitsArrayTracks);
      }
      else{
         fBGHandler[iCut] = new AliGammaConversionAODBGHandler(9,5,((AliDalitzElectronCuts*)fCutElectronArray->At(iCut))->NumberOfRotationEvents());
         fBGHandler[iCut]->Initialize(zBinLimitsArray, multiplicityBinLimitsArrayV0s);
      }
      if( ( (AliDalitzElectronCuts*)fCutElectronArray->At(iCut))->GetBKGMethod() == 3 ){
         fGammasPool[iCut] = new TList();
      }
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
   hNEvents                        = new TH1I*[fnCuts];
   hNGoodESDTracks                 = new TH1I*[fnCuts];
   hESDConvGammaPt                 = new TH1F*[fnCuts];
   hESDDalitzElectronPt            = new TH1F*[fnCuts];
   hESDDalitzPositronPt            = new TH1F*[fnCuts];
   hESDEposEnegPsiPairDPhi         = new TH2F*[fnCuts];
   hESDEposEnegInvMassPt           = new TH2F*[fnCuts];
   hESDEposEnegLikeSignBackInvMassPt = new TH2F*[fnCuts];
   hESDMotherInvMassPt             = new TH2F*[fnCuts];
   hESDPi0MotherInvMassPt          = new TH2F*[fnCuts];
   hESDPi0MotherDiffInvMassPt      = new TH2F*[fnCuts];
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

      hESDConvGammaPt[iCut] = new TH1F("ESD_ConvGamma_Pt","ESD_ConvGamma_Pt",250,0,25);
      fESDList[iCut]->Add(hESDConvGammaPt[iCut]);

      hESDDalitzElectronPt[iCut] = new TH1F("ESD_DalitzElectron_Pt","ESD_DalitzElectron_Pt",250,0,25);
      fESDList[iCut]->Add(hESDDalitzElectronPt[iCut]);

      hESDDalitzPositronPt[iCut] = new TH1F("ESD_DalitzPositron_Pt","ESD_DalitzPositron_Pt",250,0,25);
      fESDList[iCut]->Add(hESDDalitzPositronPt[iCut]);


      hESDEposEnegPsiPairDPhi[iCut] = new TH2F("ESD_EposEneg_PsiPair_DPhi","ESD_EposEneg_PsiPair_DPhi", 100, -1.0,1.0,100,-1.0,1.0 );
      fESDList[iCut]->Add(hESDEposEnegPsiPairDPhi[iCut]);

      hESDEposEnegInvMassPt[iCut] = new TH2F("ESD_EposEneg_InvMassPt","ESD_EposEneg_InvMassPt",5000,0.,5.,100,0.,10.);
      fESDList[iCut]->Add(hESDEposEnegInvMassPt[iCut]);
      
      hESDEposEnegLikeSignBackInvMassPt[iCut]  = new TH2F("ESD_EposEneg_LikeSignBack_InvMassPt","ESD_EposEneg_LikeSignBack_InvMassPt",5000,0.,5.,100,0.,10.);
      fESDList[iCut]->Add(hESDEposEnegLikeSignBackInvMassPt[iCut]);


      hESDMotherInvMassPt[iCut] = new TH2F("ESD_DalitzMother_InvMass_Pt","ESD_DalitzMother_InvMass_Pt",1000,0,1,250,0,25);
      fESDList[iCut]->Add(hESDMotherInvMassPt[iCut]);

      hESDPi0MotherInvMassPt[iCut] = new TH2F("ESD_Pi0Mother_InvMass_Pt","ESD_Pi0Mother_InvMass_Pt",4000,0,4,250,0,25);
      fESDList[iCut]->Add(hESDPi0MotherInvMassPt[iCut]);
  
      hESDPi0MotherDiffInvMassPt[iCut] = new TH2F("ESD_Pi0Mother_DiffInvMass_Pt","ESD_Pi0Mother_DiffInvMass_Pt",2000,0,2,250,0,25);
      fESDList[iCut]->Add(hESDPi0MotherDiffInvMassPt[iCut]);


      hESDMotherBackInvMassPt[iCut] = new TH2F("ESD_DalitzBackground_InvMass_Pt","ESD_DalitzBackground_InvMass_Pt",1000,0,1,250,0,25);
      fESDList[iCut]->Add(hESDMotherBackInvMassPt[iCut]);

      fCutFolder[iCut]->Add(fESDList[iCut]);



   }


   InitBack(); // Init Background Handler


   if(MCEvent()){
      // MC Histogramms
      fMCList = new TList*[fnCuts];
      // True Histogramms
      fTrueList = new TList*[fnCuts];
      hESDTrueConvGammaPt = new TH1F*[fnCuts];
      hESDTruePositronPt  = new TH1F*[fnCuts];
      hESDTrueElectronPt  = new TH1F*[fnCuts];
      hESDTruePi0DalitzConvGammaPt = new TH1F*[fnCuts];
      hESDTruePi0DalitzPositronPt  = new TH1F*[fnCuts];
      hESDTruePi0DalitzElectronPt  = new TH1F*[fnCuts];
      //if(fDoMesonAnalysis){
      hMCAllGammaPt = new TH1F*[fnCuts];
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
      
      
      hESDEposEnegTruePi0DalitzInvMassPt           = new TH2F*[fnCuts];
      hESDEposEnegTrueEtaDalitzInvMassPt           = new TH2F*[fnCuts];
      hESDEposEnegTruePhotonInvMassPt              = new TH2F*[fnCuts];
      
      hESDTrueMotherInvMassPt = new TH2F*[fnCuts];
      hESDTrueMotherPi0GGInvMassPt = new TH2F*[fnCuts];
      hESDTruePrimaryPi0DalitzESDPtMCPt = new TH2F*[fnCuts];
      hESDTruePrimaryMotherInvMassMCPt = new TH2F*[fnCuts];
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

                hMCAllPositronsPt[iCut] = new TH1F("MC_AllPositrons_Pt","MC_AllPositrons_Pt",250,0,25);
                fMCList[iCut]->Add(hMCAllPositronsPt[iCut]);

                hMCAllElectronsPt[iCut] = new TH1F("MC_AllElectrons_Pt","MC_AllElectrons_Pt",250,0,25);
                fMCList[iCut]->Add(hMCAllElectronsPt[iCut]);

                hMCPi0DalitzGammaPt[iCut] = new TH1F("MC_Pi0DalitzGamma_Pt","MC_Pi0DalitzGamma_Pt",250,0,25);
                fMCList[iCut]->Add(hMCPi0DalitzGammaPt[iCut]);

                hMCPi0DalitzPositronPt[iCut] = new TH1F("MC_Pi0DalitzPositron_Pt","MC_Pi0DalitzPositron_Pt",250,0,25);
                fMCList[iCut]->Add(hMCPi0DalitzPositronPt[iCut]);

                hMCPi0DalitzElectronPt[iCut] = new TH1F("MC_Pi0DalitzElectron_Pt","MC_Pi0DalitzElectron_Pt",250,0,25);
                fMCList[iCut]->Add(hMCPi0DalitzElectronPt[iCut]);


                hMCPi0Pt[iCut] = new TH1F("MC_Pi0_Pt","MC_Pi0_Pt",250,0,25);
                fMCList[iCut]->Add(hMCPi0Pt[iCut]);
        
                hMCPi0GGPt[iCut] = new TH1F("MC_Pi0_GG_Pt","MC_Pi0_GG_Pt",250,0,25);
                fMCList[iCut]->Add(hMCPi0GGPt[iCut]);


                hMCEtaPt[iCut] = new TH1F("MC_Eta_Pt","MC_Eta_Pt",250,0,25);
                fMCList[iCut]->Add(hMCEtaPt[iCut]);

                hMCEtaGGPt[iCut] = new TH1F("MC_Eta_GG_Pt","MC_Eta_GG_Pt",250,0,25);
                fMCList[iCut]->Add(hMCEtaGGPt[iCut]);


                hMCPi0InAccPt[iCut] = new TH1F("MC_Pi0DalitzInAcc_Pt","MC_Pi0DalitzInAcc_Pt",250,0,25);
                fMCList[iCut]->Add(hMCPi0InAccPt[iCut]);
                hMCEtaInAccPt[iCut] = new TH1F("MC_EtaDalitzInAcc_Pt","MC_EtaDalitzInAcc_Pt",250,0,25);
                fMCList[iCut]->Add(hMCEtaInAccPt[iCut]);

         fTrueList[iCut] = new TList();
         fTrueList[iCut]->SetName(Form("%s_%s_%s True histograms",cutstringGamma.Data(),cutstringElectron.Data(),cutstringMeson.Data()));
         fTrueList[iCut]->SetOwner(kTRUE);
         fCutFolder[iCut]->Add(fTrueList[iCut]);


         





	 
	 
	 hESDEposEnegTruePi0DalitzInvMassPt[iCut] = new TH2F("ESD_EposEneg_TruePi0Dalitz_InvMassPt","ESD_EposEneg_TruePi0Dalitz_InvMassPt",5000,0.,5.,100,0.,10.);
	 fTrueList[iCut]->Add(hESDEposEnegTruePi0DalitzInvMassPt[iCut]);
	 
	 hESDEposEnegTrueEtaDalitzInvMassPt[iCut] = new TH2F("ESD_EposEneg_TrueEtaDalitz_InvMassPt","ESD_EposEneg_TrueEtaDalitz_InvMassPt",5000,0.,5.,100,0.,10.);
	 fTrueList[iCut]->Add(hESDEposEnegTrueEtaDalitzInvMassPt[iCut]);
	 
	 hESDEposEnegTruePhotonInvMassPt[iCut] = new TH2F("ESD_EposEneg_TruePhoton_InvMassPt","ESD_EposEneg_TruePhoton_InvMassPt",5000,0.,5.,100,0.,10.);
	 fTrueList[iCut]->Add(hESDEposEnegTruePhotonInvMassPt[iCut]);
	 
	  

         hESDTruePositronPt[iCut] = new TH1F("ESD_TruePositron_Pt","ESD_TruePositron_Pt",250,0,25);
         fTrueList[iCut]->Add(hESDTruePositronPt[iCut]);

         hESDTrueElectronPt[iCut] = new TH1F("ESD_TrueElectron_Pt","ESD_TrueElectron_Pt",250,0,25);
         fTrueList[iCut]->Add(hESDTrueElectronPt[iCut]);

         hESDTrueConvGammaPt[iCut] = new TH1F("ESD_TrueConvGamma_Pt","ESD_TrueConvGamma_Pt",250,0,25);
         fTrueList[iCut]->Add(hESDTrueConvGammaPt[iCut]);

         hESDTruePi0DalitzConvGammaPt[iCut] = new TH1F("ESD_TruePi0DalitzConvGamma_Pt","ESD_TruePi0DalitzConvGamma_Pt",250,0,25);
         fTrueList[iCut]->Add(hESDTruePi0DalitzConvGammaPt[iCut]);

         hESDTruePi0DalitzElectronPt[iCut] = new TH1F("ESD_TruePi0DalitzElectron_Pt","ESD_TruePi0DalitzElectron_Pt",250,0,25);
         fTrueList[iCut]->Add(hESDTruePi0DalitzElectronPt[iCut]);

         hESDTruePi0DalitzPositronPt[iCut] = new TH1F("ESD_TruePi0DalitzPositron_Pt","ESD_TruePi0DalitzPositron_Pt",250,0,25);
         fTrueList[iCut]->Add(hESDTruePi0DalitzPositronPt[iCut]);





         hESDTrueMotherInvMassPt[iCut] = new TH2F("ESD_TrueMother_InvMass_Pt","ESD_TrueMother_InvMass_Pt",1000,0,1,250,0,25);
         fTrueList[iCut]->Add(hESDTrueMotherInvMassPt[iCut]);

         hESDTrueMotherPi0GGInvMassPt[iCut] = new TH2F("ESD_TrueMotherPi0GG_InvMass_Pt","ESD_TrueMotherPi0GG_InvMass_Pt",1000,0,1,250,0,25);
         fTrueList[iCut]->Add(hESDTrueMotherPi0GGInvMassPt[iCut]);
         hESDTruePrimaryPi0DalitzESDPtMCPt[iCut] = new TH2F("ESD_TruePrimaryPi0Dalitz_ESDPt_MCPt","ESD_TruePrimaryPi0Dalitz_ESDPt_MCPt",250,0,25,250,0,25);
         fTrueList[iCut]->Add(hESDTruePrimaryPi0DalitzESDPtMCPt[iCut]);
         hESDTruePrimaryMotherInvMassMCPt[iCut] = new TH2F("ESD_TruePrimaryMother_InvMass_MCPt","ESD_TrueDalitzPrimaryMother_InvMass_MCPt",1000,0,1,250,0,25);
         fTrueList[iCut]->Add(hESDTruePrimaryMotherInvMassMCPt[iCut]);
         hESDTrueSecondaryMotherInvMassPt[iCut] = new TH2F("ESD_TrueDalitzSecondaryMother_InvMass_Pt","ESD_TrueDalitzSecondaryMother_InvMass_Pt",1000,0,1,250,0,25);
         fTrueList[iCut]->Add(hESDTrueSecondaryMotherInvMassPt[iCut]);
         // 				hESDTrueSecondaryMotherFromK0sInvMassPt[iCut] = new TH2F("ESD_TrueDalitzSecondaryMotherFromK0s_InvMass_Pt","ESD_TrueDalitzSecondaryMotherFromK0s_InvMass_Pt",1000,0,1,250,0,25);
         // 				fTrueList[iCut]->Add(hESDTrueSecondaryMotherFromK0sInvMassPt[iCut]);
         hESDTrueBckGGInvMassPt[iCut] = new TH2F("ESD_TrueDalitzBckGG_InvMass_Pt","ESD_TrueDalitzBckGG_InvMass_Pt",1000,0,1,250,0,25);
         fTrueList[iCut]->Add(hESDTrueBckGGInvMassPt[iCut]);
         hESDTrueBckContInvMassPt[iCut] = new TH2F("ESD_TrueDalitzBckCont_InvMass_Pt","ESD_TrueDalitzBckCont_InvMass_Pt",1000,0,1,250,0,25);
         fTrueList[iCut]->Add(hESDTrueBckContInvMassPt[iCut]);
         // 				hESDTrueMotherGGInvMassPt[iCut] = new TH2F("ESD_TrueGammaGamma_InvMass_Pt","ESD_TrueGammaGamma_InvMass_Pt",1000,0,1,250,0,25);
         // 				fTrueList[iCut]->Add(hESDTrueMotherGGInvMassPt[iCut]);

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


   fMCEvent         =  MCEvent();
   fESDEvent        = (AliESDEvent*)InputEvent();
   fReaderGammas    = fV0Reader->GetReconstructedGammas(); // Gammas from default Cut
   fSelectorElectronIndex = fElecSelector->GetReconstructedElectronsIndex(); // Electrons from default Cut
   fSelectorPositronIndex = fElecSelector->GetReconstructedPositronsIndex(); // Positrons from default Cut

   CountESDTracks(); // Estimate Event Multiplicity
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
         ProcessMCParticles();
      }



      ProcessPhotonCandidates(); // Process this cuts gammas
      ProcessElectronCandidates(); // Process this cuts gammas
      CalculatePi0DalitzCandidates();
      CalculateBackground();
      UpdateEventByEventData();



      fGoodGammas->Clear(); // delete this cuts good gammas
      fGoodVirtualGammas->Clear(); // delete this cuts good gammas
   }

   fSelectorElectronIndex.clear();
   fSelectorPositronIndex.clear();

   PostData( 1, fOutputContainer );
}

void AliAnalysisTaskGammaConvDalitzV1::Terminate(const Option_t *)
{

   if( fElecSelector ){

      if ( ((AliDalitzElectronCuts*)fElecSelector->GetDalitzElectronCuts())->GetCutHistograms() ){
         fOutputContainer->Add( ((AliDalitzElectronCuts*)fElecSelector->GetDalitzElectronCuts())->GetCutHistograms() );
      }
   }


   if ( fV0Reader ) {

      if ( ((AliConversionCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms() ){
         fOutputContainer->Add(  ((AliConversionCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms() );
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
      if(!((AliConversionCuts*)fCutGammaArray->At(fiCut))->PhotonIsSelected(PhotonCandidate,fESDEvent)) continue;

      if(!((AliConversionCuts*)fCutGammaArray->At(fiCut))->UseElecSharingCut() &&
         !((AliConversionCuts*)fCutGammaArray->At(fiCut))->UseToCloseV0sCut()){ // if no post reader loop is required add to events good gammas
         fGoodGammas->Add(PhotonCandidate);
         hESDConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
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
         if(!((AliConversionCuts*)fCutGammaArray->At(fiCut))->RejectSharedElectronV0s(PhotonCandidate,i,GoodGammasStepOne->GetEntries())) continue;
         if(!((AliConversionCuts*)fCutGammaArray->At(fiCut))->UseToCloseV0sCut()){ // To Colse v0s cut diabled, step two not needed
            fGoodGammas->Add(PhotonCandidate);
            hESDConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt());
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
         if(!((AliConversionCuts*)fCutGammaArray->At(fiCut))->RejectToCloseV0s(PhotonCandidate,GoodGammasStepTwo,i)) continue;
         fGoodGammas->Add(PhotonCandidate); // Add gamma to current cut TList
         hESDConvGammaPt[fiCut]->Fill(PhotonCandidate->Pt()); // Differences to old V0Reader in p_t due to conversion KF->TLorentzVector
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
   if(TMath::Abs(posDaughter->GetPdgCode())!=11 || TMath::Abs(negDaughter->GetPdgCode())!=11) return; //One Particle is not electron
   if(posDaughter->GetPdgCode()==negDaughter->GetPdgCode()) return; // Same Charge
   if(posDaughter->GetUniqueID() != 5 || negDaughter->GetUniqueID() !=5) return;// check if the daughters come from a conversion

   TParticle *Photon = TruePhotonCandidate->GetMCParticle(MCStack);
   if(Photon->GetPdgCode() != 22) return; // Mother is no Photon

   // True Photon
   hESDTrueConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
  
   Int_t labelGamma = TruePhotonCandidate->GetMCParticleLabel(MCStack);

    if( IsPi0DalitzDaughter(labelGamma) == kTRUE ) {
        hESDTruePi0DalitzConvGammaPt[fiCut]->Fill(TruePhotonCandidate->Pt());
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




   vector<Int_t> lGoodElectronIndex(0);
   vector<Int_t> lGoodPositronIndex(0);



   for(UInt_t i = 0; i < fSelectorElectronIndex.size(); i++){

      AliESDtrack* electronCandidate = fESDEvent->GetTrack(fSelectorElectronIndex[i]);
      if(! ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->ElectronIsSelected(electronCandidate) ) continue;
         lGoodElectronIndex.push_back(   fSelectorElectronIndex[i] );
         hESDDalitzElectronPt[fiCut]->Fill(electronCandidate->Pt());
            if( fMCEvent ) {
                    Int_t labelelectron = TMath::Abs( electronCandidate->GetLabel() );
                    if( labelelectron < fMCStack->GetNtrack() ){
                    TParticle* electron = fMCStack->Particle(labelelectron);
                      if( electron->GetPdgCode() ==  11 ){
                            hESDTrueElectronPt[fiCut]->Fill(electron->Pt());
                        if( IsPi0DalitzDaughter(labelelectron) == kTRUE ) {
                            hESDTruePi0DalitzElectronPt[fiCut]->Fill(electron->Pt());
                       }
                     }
                   }
            }
   }

   for(UInt_t i = 0; i < fSelectorPositronIndex.size(); i++){

      AliESDtrack* positronCandidate = fESDEvent->GetTrack( fSelectorPositronIndex[i] );
      if(! ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->ElectronIsSelected(positronCandidate) ) continue;
        lGoodPositronIndex.push_back(   fSelectorPositronIndex[i] );
        hESDDalitzPositronPt[fiCut]->Fill(positronCandidate->Pt());
        if( fMCEvent ) {
                    Int_t labelpositron = TMath::Abs( positronCandidate->GetLabel() );
                    if( labelpositron < fMCStack->GetNtrack() ) {
                    TParticle* positron = fMCStack->Particle(labelpositron);
                      if( positron->GetPdgCode() ==  -11 ){
                            hESDTruePositronPt[fiCut]->Fill(positron->Pt());
                        if( IsPi0DalitzDaughter(labelpositron) == kTRUE ) {
                            hESDTruePi0DalitzPositronPt[fiCut]->Fill(positron->Pt());
                       }
                     }
               }
        }

   }


   vector<Bool_t> lElectronPsiIndex(lGoodElectronIndex.size(), kTRUE);
   vector<Bool_t> lPositronPsiIndex(lGoodPositronIndex.size(), kTRUE);


   if( ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->DoPsiPairCut() == kTRUE ){

      for( UInt_t i = 0; i < lGoodElectronIndex.size(); i++ ) {

         AliESDtrack *electronCandidate = fESDEvent->GetTrack(lGoodElectronIndex[i]);

         for(UInt_t j = 0; j <  lGoodPositronIndex.size(); j++){
            AliESDtrack *positronCandidate = fESDEvent->GetTrack(lGoodPositronIndex[j]);
            Double_t psiPair = GetPsiPair(positronCandidate,electronCandidate);
            Double_t deltaPhi = magField * TVector2::Phi_mpi_pi( electronCandidate->GetConstrainedParam()->Phi()-positronCandidate->GetConstrainedParam()->Phi());

            if( ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->IsFromGammaConversion(psiPair,deltaPhi ) ){
               lElectronPsiIndex[i] = kFALSE;
               lPositronPsiIndex[j] = kFALSE;
            }
         }
      }
   }






   for(UInt_t i = 0; i < lGoodElectronIndex.size(); i++){

      if( lElectronPsiIndex[i] == kFALSE ) continue;


      AliESDtrack *electronCandidate = fESDEvent->GetTrack(lGoodElectronIndex[i]);

      AliKFParticle electronCandidateKF( *electronCandidate->GetConstrainedParam(), ::kElectron );

      for(UInt_t j = 0; j < lGoodPositronIndex.size(); j++){

         if( lPositronPsiIndex[j] == kFALSE ) continue;

         AliESDtrack *positronCandidate = fESDEvent->GetTrack(lGoodPositronIndex[j]);
         AliKFParticle positronCandidateKF( *positronCandidate->GetConstrainedParam(), ::kPositron );
	 Bool_t isPhoton    = kFALSE;
	 Bool_t isPi0Dalitz = kFALSE;
	 Bool_t isEtaDalitz = kFALSE;
	 

         Double_t psiPair = GetPsiPair(positronCandidate,electronCandidate);
         Double_t deltaPhi = magField * TVector2::Phi_mpi_pi( electronCandidate->GetConstrainedParam()->Phi()-positronCandidate->GetConstrainedParam()->Phi());
         


         AliKFConversionPhoton* virtualPhoton = new AliKFConversionPhoton(electronCandidateKF,positronCandidateKF);


         AliKFVertex primaryVertexImproved(*fInputEvent->GetPrimaryVertex());
         primaryVertexImproved+=*virtualPhoton;
         virtualPhoton->SetProductionVertex(primaryVertexImproved);

         virtualPhoton->SetTrackLabels( lGoodPositronIndex[j], lGoodElectronIndex[i]);

         if(fMCEvent){

            // AliStack *fMCStack= fMCEvent->Stack();
            Int_t labeln=TMath::Abs(electronCandidate->GetLabel());
            Int_t labelp=TMath::Abs(positronCandidate->GetLabel());
            TParticle *fNegativeMCParticle = fMCStack->Particle(labeln);
            TParticle *fPositiveMCParticle = fMCStack->Particle(labelp);
            if( fPositiveMCParticle && fNegativeMCParticle) {
               virtualPhoton->SetMCLabelPositive(labelp);
               virtualPhoton->SetMCLabelNegative(labeln);
            
	    }
            
            TParticle *mcVgamma=virtualPhoton->GetMCParticle(fMCStack);
	    Int_t labelgamma    = -1;
	    Int_t labelelectron = -1;
	    Int_t labelpositron = -1;

	    if(mcVgamma){
	      // Check if it is a true photon
	      if(mcVgamma->GetPdgCode() == 22){
		isPhoton = kTRUE;
	      }
	      else if( IsDalitz( mcVgamma,labelgamma,labelelectron,labelpositron ) ){
		  if     ( mcVgamma->GetPdgCode() == 111 ) isPi0Dalitz = kTRUE;
		  else if( mcVgamma->GetPdgCode() == 221 ) isEtaDalitz = kTRUE;
	      }
	    }
	   }

            
         AliAODConversionPhoton *vphoton = new AliAODConversionPhoton(virtualPhoton); //To Apply PsiPairCut

         hESDEposEnegPsiPairDPhi[fiCut]->Fill(deltaPhi,psiPair);   
         hESDEposEnegInvMassPt[fiCut]->Fill(vphoton->GetMass(),vphoton->Pt());
	 
	 if( fMCEvent ) {
	    if(isPhoton) 	hESDEposEnegTruePhotonInvMassPt[fiCut]->Fill(vphoton->GetMass(),vphoton->Pt());
	    else if(isPi0Dalitz)hESDEposEnegTruePi0DalitzInvMassPt[fiCut]->Fill(vphoton->GetMass(),vphoton->Pt());
	    else if(isEtaDalitz)hESDEposEnegTrueEtaDalitzInvMassPt[fiCut]->Fill(vphoton->GetMass(),vphoton->Pt());
	    
	 }
	 
	 fGoodVirtualGammas->Add(  vphoton );
	 	 	 
      }
   }
   
   
   //Computing mixing event

   for(UInt_t i = 0; i < lGoodElectronIndex.size(); i++){
         
         if( lElectronPsiIndex[i] == kFALSE ) continue;
         
          AliESDtrack *electronCandidate1 = fESDEvent->GetTrack(lGoodElectronIndex[i]);

          AliKFParticle electronCandidate1KF( *electronCandidate1->GetConstrainedParam(), ::kElectron );
         
        
          for(UInt_t j = i+1; j < lGoodElectronIndex.size(); j++){
	    
	       if( lElectronPsiIndex[j] == kFALSE ) continue;
	       
	       
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
		            
	   
	  }
   }   
   
   
   for(UInt_t i = 0; i < lGoodPositronIndex.size(); i++){
         
         if( lPositronPsiIndex[i] == kFALSE ) continue;
         
          AliESDtrack *positronCandidate1 = fESDEvent->GetTrack(lGoodPositronIndex[i]);

          AliKFParticle positronCandidate1KF( *positronCandidate1->GetConstrainedParam(), ::kPositron );
         
        
          for(UInt_t j = i+1; j < lGoodPositronIndex.size(); j++){
	    
	       if( lPositronPsiIndex[j] == kFALSE ) continue;
	       
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
	   
	  }
   }   


}

//________________________________________________________________________
void AliAnalysisTaskGammaConvDalitzV1::CalculatePi0DalitzCandidates(){

   // Conversion Gammas




   if( fGoodGammas->GetEntries() > 0 && fGoodVirtualGammas->GetEntries() > 0 ){


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

            if( ( ((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->MesonIsSelected(pi0cand,kTRUE)) ){
               hESDMotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
               hESDPi0MotherInvMassPt[fiCut]->Fill(pi0cand->M(),pi0cand->Pt());
                
               Double_t diffMass = pi0cand->M() - Vgamma->GetMass();
    
               hESDPi0MotherDiffInvMassPt[fiCut]->Fill( diffMass , pi0cand->Pt() );
               // if(pi0cand->GetAlpha()<0.1){
               //    hESDMotherInvMassEalpha[fiCut]->Fill(pi0cand->M(),pi0cand->E());
               // }
               Int_t zbin= fBGHandler[fiCut]->GetZBinIndex(fESDEvent->GetPrimaryVertex()->GetZ());
               Int_t mbin = 0;
               if( ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->UseTrackMultiplicity() ){
                  mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fNumberOfESDTracks);
               } else {
                  mbin = fBGHandler[fiCut]->GetMultiplicityBinIndex(fGoodGammas->GetEntries());
               }
               Double_t sparesFill[4] = {pi0cand->M(),pi0cand->Pt(),(Double_t)zbin,(Double_t)mbin};
               sESDMotherInvMassPtZM[fiCut]->Fill(sparesFill,1);
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
                  if( ( ((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kFALSE) ) ){
                     hESDMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt());
                     Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
                     sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
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
                     if((((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kFALSE))){
                        hESDMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt());
                        Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
                        sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
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

            if((((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->MesonIsSelected(backgroundCandidate,kFALSE))){

               hESDMotherBackInvMassPt[fiCut]->Fill(backgroundCandidate->M(),backgroundCandidate->Pt());
               Double_t sparesFill[4] = {backgroundCandidate->M(),backgroundCandidate->Pt(),(Double_t)zbin,(Double_t)mbin};
               sESDMotherBackInvMassPtZM[fiCut]->Fill(sparesFill,1);
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


   if(TrueGammaCandidate->GetV0Index()<fESDEvent->GetNumberOfV0s()){


      Bool_t isTruePi0 = kFALSE;
      Bool_t isTrueEta = kFALSE;
      Int_t gammaMCLabel = TrueGammaCandidate->GetMCParticleLabel(MCStack);
      Int_t gammaMotherLabel = -1;




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

      if(virtualGammaMCLabel != -1){ // Gamma is Combinatorial; MC Particles don't belong to the same Mother
         // Daughters Gamma 1
         TParticle * negativeMC = (TParticle*)TrueVirtualGammaCandidate->GetNegativeMCDaughter(MCStack);
         TParticle * positiveMC = (TParticle*)TrueVirtualGammaCandidate->GetPositiveMCDaughter(MCStack);
         TParticle * virtualGammaMotherMC = (TParticle*)MCStack->Particle(virtualGammaMCLabel);

         if(TMath::Abs(negativeMC->GetPdgCode())==11 && TMath::Abs(positiveMC->GetPdgCode())==11){  // Electrons ...

            if( virtualGammaMotherMC->GetPdgCode() != 22 ){
               virtualGammaMotherLabel=virtualGammaMCLabel;
            }

            else if(negativeMC->GetUniqueID() == 5 && positiveMC->GetUniqueID() ==5){ // ... From Conversion ...

               virtualGammaMotherLabel=virtualGammaMotherMC->GetFirstMother();
               virtualGamma = 0; //no virtual gamma

            }
         }
      }


      if(gammaMotherLabel>=0 && ( gammaMotherLabel == virtualGammaMotherLabel) ){

         if(((TParticle*)MCStack->Particle(virtualGammaMotherLabel))->GetPdgCode() == 111){
            isTruePi0=kTRUE;
         }

         if(((TParticle*)MCStack->Particle(virtualGammaMotherLabel))->GetPdgCode() == 221){
            isTrueEta=kTRUE;
         }

      }

      if(isTruePi0 || isTrueEta ){ // True Pion or Eta
         if ( virtualGamma == 1 ) { //True Dalitz
            hESDTrueMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
            if(gammaMotherLabel <= MCStack->GetNprimary()){ // Only primary pi0 for efficiency calculation
               hESDTruePrimaryMotherInvMassMCPt[fiCut]->Fill(Pi0Candidate->M(),((TParticle*)MCStack->Particle(virtualGammaMotherLabel))->Pt());
               if(isTruePi0){ // Only primaries for unfolding
                  hESDTruePrimaryPi0DalitzESDPtMCPt[fiCut]->Fill(Pi0Candidate->Pt(),((TParticle*)MCStack->Particle(virtualGammaMotherLabel))->Pt());
               }
            }
            if(gammaMotherLabel > MCStack->GetNprimary()){ // Secondary Meson
               hESDTrueSecondaryMotherInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
               //if (((TParticle*)MCStack->Particle(virtualGammaMotherLabel))->GetMother(0) >-1){
               //   if(MCStack->Particle(((TParticle*)MCStack->Particle(virtualGammaMotherLabel))->GetMother(0))->GetPdgCode()==kK0Short){
               //      hESDTrueSecondaryMotherFromK0sInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt());
               //   }
            }
         }
         else if ( virtualGamma == 0 ){
            hESDTrueMotherPi0GGInvMassPt[fiCut]->Fill(Pi0Candidate->M(),Pi0Candidate->Pt()); // Pi0 from GG
         }
      }

      if(!isTruePi0 && !isTrueEta){ // Background
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


    if(((AliConversionCuts*)fCutGammaArray->At(fiCut))->PhotonIsSelectedMC(particle,fMCStack,kFALSE)){
         hMCAllGammaPt[fiCut]->Fill(particle->Pt()); // All MC Gamma
    }

    if(  ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->ElectronIsSelectedMC(i,fMCStack)){
        if( particle->GetPdgCode() == -11)hMCAllPositronsPt[fiCut]->Fill(particle->Pt()); // All positrons
        if( particle->GetPdgCode() ==  11)hMCAllElectronsPt[fiCut]->Fill(particle->Pt()); // All electrons
    }

    if(((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->MesonIsSelectedMC(particle,fMCStack)){

            if(particle->GetPdgCode() == 111)hMCPi0GGPt[fiCut]->Fill( particle->Pt() ); // All MC Pi0 GG decay
            if(particle->GetPdgCode() == 221)hMCEtaGGPt[fiCut]->Fill( particle->Pt() ); // All MC Eta GG decay
    }

    if(((AliConversionMesonCuts*)fCutMesonArray->At(fiCut))->MesonIsSelectedMCDalitz(particle,fMCStack)){
	
         if(particle->GetPdgCode() == 111)hMCPi0Pt[fiCut]->Fill(particle->Pt()); // All MC Pi0
         if(particle->GetPdgCode() == 221)hMCEtaPt[fiCut]->Fill(particle->Pt()); // All MC Eta

         // Check the acceptance for gamma and electrons
	 
	 Int_t labelgamma =    -1;
	 Int_t labelelectron = -1;
	 Int_t labelpositron = -1;
	 
        if( IsDalitz( particle,labelgamma,labelelectron,labelpositron) == kTRUE ) {
	  
	  
	       TParticle *gamma    = fMCStack->Particle(labelgamma);
	       TParticle *electron = fMCStack->Particle(labelelectron);
	       TParticle *positron = fMCStack->Particle(labelpositron);
        
                      
               if(((AliConversionCuts*)fCutGammaArray->At(fiCut))->PhotonIsSelectedMC(gamma,fMCStack,kFALSE) &&
                  TMath::Abs( electron->Eta() ) < ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->GetEtaCut()  &&
                  TMath::Abs( positron->Eta() ) < ((AliDalitzElectronCuts*)fCutElectronArray->At(fiCut))->GetEtaCut() ){
                  if(particle->GetPdgCode() == 111){ 
                        hMCPi0InAccPt[fiCut]->Fill(particle->Pt()); // MC Pi0Dalitz with gamma and e+e- in acc
                        hMCPi0DalitzGammaPt[fiCut]->Fill( gamma->Pt());
                        hMCPi0DalitzPositronPt[fiCut]->Fill( positron->Pt());
                        hMCPi0DalitzElectronPt[fiCut]->Fill( electron->Pt());
                  }
                  if(particle->GetPdgCode() == 221)hMCEtaInAccPt[fiCut]->Fill(particle->Pt()); // MC EtaDalitz with gamma and e+e- in acc
               }
         }
      }
   }
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskGammaConvDalitzV1::IsDalitz(TParticle *fMCMother,Int_t &labelgamma, Int_t &labelelectron,Int_t &labelpositron)
{

	    if( fMCMother->GetNDaughters() != 3 ) return kFALSE;
	    if( fMCMother->GetPdgCode() != 111 && fMCMother->GetPdgCode() != 221 ) return kFALSE;
	    
	    TParticle *positron = 0x0;
	    TParticle *electron = 0x0;
	    TParticle *gamma 	= 0x0;
	    	    
	      


            for(Int_t index= fMCMother->GetFirstDaughter();index<= fMCMother->GetLastDaughter();index++){


               TParticle* temp = (TParticle*)fMCStack->Particle( index );

               switch( temp->GetPdgCode() ) {
               case ::kPositron:
				positron =  temp;
				labelpositron = index;
				break;
               case ::kElectron:
				electron =  temp;
				labelelectron = index;
				break;
               case ::kGamma:
				gamma    =  temp;
				labelgamma = index;
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
        Bool_t ePlusFlag  = kFALSE;
        Bool_t eMinusFlag = kFALSE;
        Bool_t gammaFlag  = kFALSE;
        
        Int_t motherLabel = fMCStack->Particle( label )->GetMother(0);
        
        if( motherLabel < 0 ) return kFALSE;
        
        TParticle* mother = fMCStack->Particle( motherLabel );
        
        if ( mother->GetPdgCode() != ::kPi0 ) return kFALSE;
        
        if ( mother->GetNDaughters() != 3 ) return kFALSE;
        
        for( Int_t idx = mother->GetFirstDaughter(); idx <= mother->GetLastDaughter(); ++idx ) {
                switch( fMCStack->Particle(idx)->GetPdgCode()) {
                        case ::kPositron:
                                ePlusFlag  = kTRUE;
                                break;
                        case ::kElectron:
                                eMinusFlag = kTRUE;
                                break;
                        case ::kGamma:
                                gammaFlag  = kTRUE;
                                break;
                }
        }
        return ( ePlusFlag && eMinusFlag && gammaFlag );
}


//_____________________________________________________________________________
Double_t AliAnalysisTaskGammaConvDalitzV1::GetPsiPair( const AliESDtrack *trackPos, const AliESDtrack *trackNeg ) const
{
   //
   // This angle is a measure for the contribution of the opening in polar
   // direction Δ0 to the opening angle ξ Pair
   //
   // Ref. Measurement of photons via conversion pairs with the PHENIX experiment at RHIC
   //      Master Thesis. Thorsten Dahms. 2005
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
