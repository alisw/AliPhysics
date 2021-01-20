/*************************************************************************
 * Copyright(c) 1998-2019, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: Zhongbao Yin                                                  *
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

/* The task selects candidates for K0s and Lambda (trigger particles)
 * and calculates correlations with charged unidentified particles (associated particles) in phi and eta. 
 * The task works with AOD events only and containes also mixing for acceptance corrections.
 */
#include <iostream>
#include <TCanvas.h>
#include <TList.h>
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <THn.h>
#include "TH1I.h"
#include "TH1F.h"
#include "TH2F.h"
#include <THnSparse.h>
#include "TObjArray.h"
#include "TGrid.h"
#include "AliLog.h"
#include "AliAnalysisManager.h"
#include "AliESDv0.h"
#include "AliESDcascade.h"
#include "AliESDtrackCuts.h"
#include "AliESDUtils.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODEvent.h"
#include "AliAODv0.h"
#include "AliAODcascade.h"
#include "AliAODVertex.h"
#include "AliVVertex.h"
#include "AliVEvent.h"
#include "AliEventCuts.h"

#include "AliAODMCHeader.h"
#include "AliAODMCParticle.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCVertex.h"
#include "AliStack.h"
#include "AliGenEventHeader.h"

#include "AliAODPid.h"
#include "AliPIDResponse.h"
#include "AliEventPoolManager.h"
#include "AliCentrality.h"

#include "AliAnalysisTaskV0ChCorrelationsys.h"
#include "AliMultSelection.h"
#include "AliPhysicsSelectionTask.h"
#include "AliMultSelectionTask.h"

//#include <array>
#include <map>
#include <AliMultiInputEventHandler.h>
#include <AliMultiEventInputHandler.h>
//#include <AliMixEventInputHandler.h>

using namespace std;
//using std::array;
ClassImp(AliAnalysisTaskV0ChCorrelationsys)
ClassImp(AliV0XiParticles)
//___________________________________________________________________________________________________
AliAnalysisTaskV0ChCorrelationsys::AliAnalysisTaskV0ChCorrelationsys()  
     :AliAnalysisTaskSE(),
     fEventCuts(0),
     fMixingTracks(0),
     fPoolSize(0),
     fPoolMgr(0x0),
     fOutput(0),
     fOutput2(0),
     fOutput3(0),
     fOutput4(0),
     fOutput5(0),
     fOutput6(0),
     fOutput7(0),
     fOutput8(0),
     //fMultiplicityV0McorrCut(0),
     fPIDResponse(0),
     fPrimaryVertexCut(0),
     fNumOfVzBins(9),
     fVtxXMin(0),
     fVtxYMin(0),
     fVtxZMin(0),
     fCentMin(0),
     fCentMax(0),
     //---------------------------------Track--------------------------------
     fTrackMCPtMin(0),
     fTrackPtMin(0),
     fTrackPtMax(0),
     fTrackEta(0),
     fFilterBit(768),
     fAssocNcls(0),
     //-----------------------------------V0---------------------------------
     fV0MCPtMin(0),
     fV0PtMin(0),
     fV0PtMax(0),
     fV0Eta(0),
    // fMinCtau(0),
    // fMaxCtau(0),
     fK0sLifeTimeMin(0),
     fK0sLifeTimeMax(0),
     fLambdaLifeTimeMin(0),
     fLambdaLifeTimeMax(0),

     fV0DaughterPtMinCut(0),

     fDCANegtoPrimVertexMink0s(0),
     fDCAPostoPrimVertexMink0s(0),
     fDCANegtoPrimVertexMinLamb(0),
     fDCAPostoPrimVertexMinLamb(0),
     fDCANegtoPrimVertexMinALamb(0),
     fDCAPostoPrimVertexMinALamb(0),

     fDCAV0DaughtersMax(0),
    // fCPA(0),

     fOStatus(1),
     fLambdaCPA(0),
     fCosPointingAngleMin(0),
     f2DFiducialMin(0),
     //f2DFiducialMax(0),
    
     fV0DaughterTrackTPCCluster(0),
     fNCrossedRowsTPCfindable(0),
    
     fK0sMassWindow(0),
     fLambdaMassWindow(0),
     fPtArmV0AlphaV0(0),
     fk0sCPA(0),
     fLambdaCosPointingAngleMin(0),
     fAntiLambdaCosPointingAngleMin(0),
     fLambdaAlphaV0Min(0),
     fAntiLambdaAlphaV0Max(0),
     fLambdaDCA2PVMax(0),
     fAntiLambdaDCA2PVMax(0),

     selectedK0s(NULL),
     selectedLambda(NULL),
     selectedAntiLambda(NULL),
     selectedTracks(NULL),
     trigParticles(NULL),

     //-----------------------------------PID----------------------------------
     fV0PIDSigma(0),
     //-----------------------------------MC-------------------------------------
     fHistRCK0sPt(0),
     fHistRCLambdaPt(0),
     fHistRCPrimLambdaDCAPt(0),
     fHistRCSecLambdaDCAPt(0),
     fHistRCAntiLambdaPt(0),
     fHistRCPrimAntiLambdaDCAPt(0),
     fHistRCSecAntiLambdaDCAPt(0),
     fHistRCTrkPtEta(0),
     fHistRCPriTrkPtEta(0),
     fHistRCSecTrkPtEta(0),
     fHistMCtruthTrkPtEta(0),
     fHistMCtruthK0sPt(0),
     fHistMCtruthLambdaPt(0),
     fHistMCtruthAntiLambdaPt(0),

     fAnalysisMC(kFALSE),
     fRejectTrackPileUp(kTRUE),
     fRejectV0PileUp(kTRUE),
     fMCArray(NULL),
     //--------------------------Correction---------------------------------
     fEffCorr(kFALSE),
     fEffList(0),
     fHistEffEtaPtK0s(0),
     fHistEffEtaPtLambda(0),
     fHistEffEtaPtAntiLambda(0),
     fHistEffEtaPtTrack(0)
{
  fMassMean[0] = 0.497614; fMassMean[1] = 1.115683; 
  fMassRes[0] = 0.005; fMassRes[1] = 0.0025; 

  for(Int_t i = 0; i < 3; i++){
    fMassLowK0s[i] = 0.;
    fMassHighK0s[i] = 0.;
    fMassLowLambda[i] = 0.;
    fMassHighLambda[i] = 0.;
    fMassLowAntiLambda[i] = 0.;
    fMassHighAntiLambda[i] = 0.;

  }

  for(Int_t i = 0 ; i < kNVtxZ; i++){
      fHistRCTrk[i] = NULL;
      fHistRCPriTrk[i] = NULL;
      fHistRCSecTrk[i] = NULL;
      fHistMCtruthTrk[i] = NULL;
      fHistMCtruthK0s[i] = NULL;
      fHistMCtruthLambda[i] = NULL;
      fHistMCtruthAntiLambda[i] = NULL;
      fHistRCK0s[i] = NULL;
      fHistRCLambda[i] = NULL;
      fHistRCAntiLambda[i] = NULL;
  } 

  for(Int_t i = 0; i < 3; i++){
      fBestPrimaryVtxPos[i] = 0;
  }
  
}
//________________________________________________________________________
AliAnalysisTaskV0ChCorrelationsys::AliAnalysisTaskV0ChCorrelationsys(const char *name,Double_t cenMin,Double_t cenMax, Bool_t effCorr)  // All data members should be initialised here
     :AliAnalysisTaskSE(name),
     fEventCuts(0),
     fMixingTracks(0),
     fPoolSize(0),
     fPoolMgr(0x0),  
     fOutput(0),
     fOutput2(0),
     fOutput3(0),
     fOutput4(0),
     fOutput5(0),
     fOutput6(0),
     fOutput7(0),
     fOutput8(0),
     //fMultiplicityV0McorrCut(0),
     fPIDResponse(0),
     fPrimaryVertexCut(0),
     fNumOfVzBins(9),
     fVtxXMin(0),
     fVtxYMin(0),
     fVtxZMin(0),
     fCentMin(cenMin),
     fCentMax(cenMax),
     //---------------------------------Track--------------------------------
     fTrackMCPtMin(0),
     fTrackPtMin(0),
     fTrackPtMax(0),
     fTrackEta(0),
     fFilterBit(768),
     fAssocNcls(0),
     //-----------------------------------V0---------------------------------
     fV0MCPtMin(0),
     fV0PtMin(0),
     fV0PtMax(0),
     fV0Eta(0),
    // fMinCtau(0),
    // fMaxCtau(0),
     fK0sLifeTimeMin(0),
     fK0sLifeTimeMax(0),
     fLambdaLifeTimeMin(0),
     fLambdaLifeTimeMax(0),
     
     fV0DaughterPtMinCut(0),
     fDCANegtoPrimVertexMink0s(0),
     fDCAPostoPrimVertexMink0s(0),
     fDCANegtoPrimVertexMinLamb(0),
     fDCAPostoPrimVertexMinLamb(0),
     fDCANegtoPrimVertexMinALamb(0),
     fDCAPostoPrimVertexMinALamb(0),
     fDCAV0DaughtersMax(0),
     //fCPA(0),
     fOStatus(1),
     fLambdaCPA(0),
     fCosPointingAngleMin(0),
     f2DFiducialMin(0),
     //f2DFiducialMax(0),
 
     fV0DaughterTrackTPCCluster(0),
     fNCrossedRowsTPCfindable(0),  

     fK0sMassWindow(0),
     fLambdaMassWindow(0),
     fPtArmV0AlphaV0(0),
     fk0sCPA(0),
     fLambdaCosPointingAngleMin(0),
     fAntiLambdaCosPointingAngleMin(0),
     fLambdaAlphaV0Min(0),
     fAntiLambdaAlphaV0Max(0),
     fLambdaDCA2PVMax(0),
     fAntiLambdaDCA2PVMax(0),

     selectedK0s(NULL),
     selectedLambda(NULL),
     selectedAntiLambda(NULL),
     selectedTracks(NULL),
     trigParticles(NULL),

     //-----------------------------------PID----------------------------------
     fV0PIDSigma(0),
     //------------------------------------MC------------------------------------

     fHistRCK0sPt(0),
     fHistRCLambdaPt(0),
     fHistRCPrimLambdaDCAPt(0),
     fHistRCSecLambdaDCAPt(0),
     fHistRCAntiLambdaPt(0),
     fHistRCPrimAntiLambdaDCAPt(0),
     fHistRCSecAntiLambdaDCAPt(0),
     fHistRCTrkPtEta(0),
     fHistRCPriTrkPtEta(0),
     fHistRCSecTrkPtEta(0),
     fHistMCtruthTrkPtEta(0),
     fHistMCtruthK0sPt(0),
     fHistMCtruthLambdaPt(0),
     fHistMCtruthAntiLambdaPt(0),
    
     fAnalysisMC(kFALSE),
     fRejectTrackPileUp(kTRUE),
     fRejectV0PileUp(kTRUE),
     fMCArray(NULL),
     //------------------------------Correction----------------------------------
     fEffCorr(effCorr),
     fEffList(0),
     fHistEffEtaPtK0s(0),
     fHistEffEtaPtLambda(0),
     fHistEffEtaPtAntiLambda(0),
     fHistEffEtaPtTrack(0)
{
   // Constructor
   // Define input and outPut slots here (never in the dummy constructor)
   // Input slot #0 works with a TChain - it is connected to the default input container
   // Output slot #1 writes into a TH1 container
  fMassMean[0] = 0.497614; fMassMean[1] = 1.115683; 
  fMassRes[0] = 0.005; fMassRes[1] = 0.0025; 
  
  for(Int_t i = 0; i < 3; i++){
    fMassLowK0s[i] = 0.;
    fMassHighK0s[i] = 0.;
    fMassLowLambda[i] = 0.;
    fMassHighLambda[i] = 0.;
    fMassLowAntiLambda[i] = 0.;
    fMassHighAntiLambda[i] = 0.;

  }

  for(Int_t i = 0 ; i < kNVtxZ; i++){
      fHistRCTrk[i] = NULL;
      fHistRCPriTrk[i] = NULL;
      fHistRCSecTrk[i] = NULL;
      fHistMCtruthTrk[i] = NULL;
      fHistMCtruthK0s[i] = NULL;
      fHistMCtruthLambda[i] = NULL;
      fHistMCtruthAntiLambda[i] = NULL;
      fHistRCK0s[i] = NULL;
      fHistRCLambda[i] = NULL;
      fHistRCAntiLambda[i] = NULL;
    } 

  for(Int_t i = 0; i < 3; i++){
      fBestPrimaryVtxPos[i] = 0;
    }

   if(fEffCorr)
    DefineInput(1, TList::Class());
    DefineOutput(1, TList::Class());
    DefineOutput(2, TList::Class());
    DefineOutput(3, TList::Class());   
    DefineOutput(4, TList::Class());   
    DefineOutput(5, TList::Class());   
    DefineOutput(6, TList::Class()); 
    DefineOutput(7, TList::Class()); 
    DefineOutput(8, TList::Class());     


}

//________________________________________________________________________
AliAnalysisTaskV0ChCorrelationsys::~AliAnalysisTaskV0ChCorrelationsys()
{
   // Destructor. Clean-up the output list, but not the histograms that are put inside
   // (the list is owner and will clean-up these histograms). Protect in PROOF case.
  if (fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()) {
      delete fOutput;
   }
 
  if(fEventCuts){
    delete fEventCuts;
  } 

  if(selectedTracks){
    delete selectedTracks;
  }

  if(selectedK0s){
    delete selectedK0s;
  } 

  if(selectedLambda){
    delete selectedLambda;
  }

  if(selectedAntiLambda){
    delete selectedAntiLambda;
  }

}

//________________________________________________________________________
void AliAnalysisTaskV0ChCorrelationsys::UserCreateOutputObjects()
{ 
   if(fEffCorr){
    fEffList = dynamic_cast<TList*>(GetInputData(1));
    if(fEffList){
      fHistEffEtaPtK0s = (TH2F*)fEffList->FindObject("fHistEffEtaPtK0sCent0_10All");
      fHistEffEtaPtLambda = (TH2F*)fEffList->FindObject("fHistEffEtaPtLambdaCent0_10All");
      fHistEffEtaPtAntiLambda = (TH2F*)fEffList->FindObject("fHistEffEtaPtAntiLambdaCent0_10All");
      fHistEffEtaPtTrack = (TH2F*)fEffList->FindObject("fHistEffEtaPtTrackCent0_10All");

      if(!fHistEffEtaPtK0s || !fHistEffEtaPtLambda || !fHistEffEtaPtAntiLambda || !fHistEffEtaPtTrack){
        std::cout<<"Efficiency histograms are not available!"<<std::endl;
      }
    }
  }

   // Create histograms
  fOutput = new TList();
  fOutput->SetOwner();  // IMPORTANT!
  fOutput->SetName("output1");

   fOutput2 = new TList();
   fOutput2->SetOwner();  // IMPORTANT!
   fOutput2->SetName("output2");

   fOutput3 = new TList();
   fOutput3->SetOwner();  // IMPORTANT!
   fOutput3->SetName("output3");

   fOutput4 = new TList();
   fOutput4->SetOwner();  // IMPORTANT!
   fOutput4->SetName("output4");

   fOutput5 = new TList();
   fOutput5->SetOwner();  // IMPORTANT!
   fOutput5->SetName("output5");

   fOutput6 = new TList();
   fOutput6->SetOwner();  // IMPORTANT!
   fOutput6->SetName("output6");

   fOutput7 = new TList();
   fOutput7->SetOwner();  // IMPORTANT!
   fOutput7->SetName("output7");

   fOutput8 = new TList();
   fOutput8->SetOwner();  // IMPORTANT!
   fOutput8->SetName("output8");

   fOutput8->Add(fOutput);
   fOutput8->Add(fOutput2);
   fOutput8->Add(fOutput3);
   fOutput8->Add(fOutput4);
   fOutput8->Add(fOutput5);
   fOutput8->Add(fOutput6);
   fOutput8->Add(fOutput7);

   
   //---------------------------------------
   fEventCuts = new AliEventCuts();
   TList *tQAEventCuts = new TList();
   tQAEventCuts->SetOwner();
   tQAEventCuts->SetName("EventCuts");
   fEventCuts->fUseVariablesCorrelationCuts = true;
   TList qaplots;
   fEventCuts->AddQAplotsToList(&qaplots, true); // fList is your output TList
//   fEventCuts->AddQAplotsToList(tQAEventCuts, true); // fList is your output TList
   // iterate the content of qaplots and add to your list:
   TObject *plot;
   TIter nextplot(&qaplots);
   while ((plot = nextplot())) tQAEventCuts->Add(plot);
   //           
   fOutput2->Add(tQAEventCuts);
   //------------------------------------

   AddQAEvent();
   AddQATrackCandidates();
   AddQAV0Candidates();
   AddQAAnalysisK0s();
   AddQAAnalysisLambda();
   AddQAAnalysisAntiLambda();


  const Int_t nPtBinsV0Xi = 47;
  const Double_t PtBinsV0Xi[48] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
                                   1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                                   2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 
                                   3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9,
                                   4.0, 4.5, 5.0, 5.5, 6.0, 8.0, 10.0, 15.0};
  const Int_t nPtBins = 48;
  const Double_t PtBins[49] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 
                               1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                               2.0, 2.1, 2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8, 2.9, 
                               3.0, 3.1, 3.2, 3.3, 3.4, 3.5, 3.6, 3.7, 3.8, 3.9,
                               4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0};  
  //defining bins of Eta distribution
  const Int_t nEtaBins=16; 
  Double_t EtaBins[nEtaBins+1] = {0.};
  EtaBins[0] = -0.8;
  for (Int_t i=0; i<nEtaBins; i++) {
  EtaBins[i+1] = EtaBins[i] + 2.*0.8/nEtaBins; }

  //defining bins of Phi distribution
  const Int_t nPhiBins = 36; 
  Double_t PhiBins[nPhiBins+1] = {0.};
  PhiBins[0] = 0;
  for (Int_t i=0; i<nPhiBins; i++) { 
  PhiBins[i+1] = PhiBins[i] + 2.*TMath::Pi()/nPhiBins; }

  if(fAnalysisMC){
  fHistRCK0sPt = new TH2F("fHistRCK0sPt", "pt of reconstructed K0s; p_{T} [GeV/c]; #eta",
                                 nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins);
  fHistRCLambdaPt = new TH2F("fHistRCLambdaPt", "pt of reconstructed Lambda; p_{T} [GeV/c]; #eta",
                                    nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins);
  fHistRCPrimLambdaDCAPt = new TH2F("fHistRCPrimLambdaDCAPt", 
                                          "DCA to PV of reco. prim. #Lambda; DCA (cm); p_{T} (GeV/c)", 
                                           80., 0., 4., nPtBinsV0Xi, PtBinsV0Xi);
  fHistRCSecLambdaDCAPt = new TH2F("fHistRCSecLambdaDCAPt",
                                         "DCA to PV of reco. sec. #Lambda; DCA (cm); p_{T} (GeV/c)",
                                          80., 0., 4., nPtBinsV0Xi, PtBinsV0Xi);
  fHistRCAntiLambdaPt = new TH2F("fHistRCAntiLambdaPt", 
                                       "pt of reconstructed #bar{#Lambda}; p_{T} [GeV/c]; #eta",
                                        nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins);
  fHistRCPrimAntiLambdaDCAPt = new TH2F("fHistRCPrimAntiLambdaDCAPt",
                                              "DCA to PV of reco. prim. #bar{#Lambda}; DCA (cm); p_{T} (GeV/c)",
                                               80., 0., 4., nPtBinsV0Xi, PtBinsV0Xi);
  fHistRCSecAntiLambdaDCAPt = new TH2F("fHistRCSecAntiLambdaDCAPt",
                                             "DCA to PV of reco. sec. #bar{#Lambda}; DCA (cm); p_{T} (GeV/c)",
                                              80., 0., 4., nPtBinsV0Xi, PtBinsV0Xi);
  fHistRCTrkPtEta = new TH2F("fHistRCTrkPtEta",
                                   "pt and eta of reconstructed tracks; p_{T} [GeV/c]; #eta",
                                    nPtBins, PtBins, nEtaBins, EtaBins);
  fHistRCPriTrkPtEta = new TH2F("fHistRCPriTrkPtEta",
                                      "pt and eta of reconstructed prim. tracks; p_{T} [GeV/c]; #eta",
                                       nPtBins, PtBins, nEtaBins, EtaBins);
  fHistRCSecTrkPtEta = new TH2F("fHistRCSecTrkPtEta",
                                      "pt and eta of reconstructed sec. tracks; p_{T} [GeV/c]; #eta",
                                       nPtBins, PtBins, nEtaBins, EtaBins);
  fHistMCtruthTrkPtEta = new TH2F("fHistMCtruthTrkPtEta",
                                        "pt and eta of generated tracks; p_{T} [GeV/c]; #eta",
                                         nPtBins, PtBins, nEtaBins, EtaBins);
  fHistMCtruthK0sPt = new TH2F("fHistMCtruthK0sPt", "pt of generated K0s; p_{T} [GeV/c]; #eta",
                                      nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins);   
  fHistMCtruthLambdaPt = new TH2F("fHistMCtruthLambdaPt", "pt of generated #Lambda; p_{T} [GeV/c]; #eta",
                                         nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins);
  fHistMCtruthAntiLambdaPt = new TH2F("fHistMCtruthAntiLambdaPt", 
                                            "pt of generated #bar{#Lambda}; p_{T} [GeV/c]; #eta",
                                             nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins);
  for(Int_t i = 0; i < kNVtxZ; i++){
    fHistRCTrk[i] = new TH3F(Form("fHistRCTrkVtx%d", i),
                                  "pt, eta, phi of reconstructed tracks; p_{T} [GeV/c]; #eta; #phi",
                                   nPtBins, PtBins, nEtaBins, EtaBins, nPhiBins, PhiBins);
    fHistRCPriTrk[i] = new TH3F(Form("fHistRCPriTrkVtx%d", i),
                                     "pt, eta, phi of reconstructed prim. tracks; p_{T} [GeV/c]; #eta; #phi",
                                      nPtBins, PtBins, nEtaBins, EtaBins, nPhiBins, PhiBins);
    fHistRCSecTrk[i] = new TH3F(Form("fHistRCSecTrkVtx%d", i),
                                     "pt, eta, phi of reconstructed sec. tracks; p_{T} [GeV/c]; #eta; #phi",
                                      nPtBins, PtBins, nEtaBins, EtaBins, nPhiBins, PhiBins);
    fHistMCtruthTrk[i] = new TH3F(Form("fHistMCtruthTrkVtx%d", i),
                                       "pt, eta, phi of generated tracks; p_{T} [GeV/c]; #eta; #phi",
                                       nPtBins, PtBins, nEtaBins, EtaBins, nPhiBins, PhiBins);
    fHistMCtruthK0s[i] = new TH3F(Form("fHistMCtruthK0sVtx%d", i), 
                                     "pt, eta and phi of generated K0s; p_{T} [GeV/c]; #eta; #phi",
                                      nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins, nPhiBins, PhiBins);
    fHistMCtruthLambda[i] = new TH3F(Form("fHistMCtruthLambdaVtx%d", i), 
                                        "pt, eta and phi of generated #Lambda; p_{T} [GeV/c]; #eta; #phi",
                                         nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins, nPhiBins, PhiBins);
    fHistMCtruthAntiLambda[i] = new TH3F(Form("fHistMCtruthAntiLambdaVtx%d", i),
                                            "pt, eta and phi of generated #bar{#Lambda}; p_{T} [GeV/c]; #eta; #phi",
                                             nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins, nPhiBins, PhiBins);
    fHistRCK0s[i] = new TH3F(Form("fHistRCK0sVtx%d", i), 
                                 "pt, eta, phi of reconstructed K0s; p_{T} [GeV/c]; #eta; #phi",
                                  nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins, nPhiBins, PhiBins);
    fHistRCLambda[i] = new TH3F(Form("fHistRCLambdaVtx%d", i), 
                                   "pt, eta and phi of reconstructed Lambda; p_{T} [GeV/c]; #eta; #phi",
                                    nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins, nPhiBins, PhiBins);
    fHistRCAntiLambda[i] = new TH3F(Form("fHistRCAntiLambdaVtx%d", i), 
                                        "pt, eta and phi of reconstructed #bar{#Lambda}; p_{T} [GeV/c]; #eta; #phi",
                                         nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins, nPhiBins, PhiBins);                               
  }
  }
  if(fAnalysisMC){
    fOutput->Add(fHistRCK0sPt);
    fOutput->Add(fHistRCLambdaPt);
    fOutput->Add(fHistRCPrimLambdaDCAPt);
    fOutput->Add(fHistRCSecLambdaDCAPt);
    fOutput->Add(fHistRCAntiLambdaPt);
    fOutput->Add(fHistRCPrimAntiLambdaDCAPt);
    fOutput->Add(fHistRCSecAntiLambdaDCAPt);
    fOutput->Add(fHistRCTrkPtEta);
    fOutput->Add(fHistRCPriTrkPtEta);
    fOutput->Add(fHistRCSecTrkPtEta);
    fOutput->Add(fHistMCtruthTrkPtEta);
    fOutput->Add(fHistMCtruthK0sPt);
    fOutput->Add(fHistMCtruthLambdaPt);
    fOutput->Add(fHistMCtruthAntiLambdaPt);
    for(Int_t i = 0 ; i < kNVtxZ; i++){
      fOutput->Add(fHistRCTrk[i]);
      fOutput->Add(fHistRCPriTrk[i]);
      fOutput->Add(fHistRCSecTrk[i]);
      fOutput->Add(fHistMCtruthTrk[i]);
      fOutput->Add(fHistMCtruthK0s[i]);
      fOutput->Add(fHistMCtruthLambda[i]);
      fOutput->Add(fHistMCtruthAntiLambda[i]);
      fOutput->Add(fHistRCK0s[i]);
      fOutput->Add(fHistRCLambda[i]);
      fOutput->Add(fHistRCAntiLambda[i]);
    }
  }
   //------------------------------------------
   PostData(1, fOutput); // Post data for ALL output slots >0 here, to get at least an empty histogram
   PostData(2, fOutput2);
   PostData(3, fOutput3);
   PostData(4, fOutput4);
   PostData(5, fOutput5);
   PostData(6, fOutput6);
   PostData(7, fOutput7);
   PostData(8, fOutput8);

   //-------------------------------------------

   fMassLowK0s[0] = fMassMean[0]-8.*fMassRes[0];
   fMassLowK0s[1] = fMassMean[0]-3.*fMassRes[0];
   fMassLowK0s[2] = fMassMean[0]+5.*fMassRes[0];
   fMassHighK0s[0] = fMassMean[0]-5.*fMassRes[0];
   fMassHighK0s[1] = fMassMean[0]+3.*fMassRes[0];
   fMassHighK0s[2] = fMassMean[0]+8.*fMassRes[0];
   
   fMassLowLambda[0] = fMassMean[1]-8.*fMassRes[1];
   fMassLowLambda[1] = fMassMean[1]-3.*fMassRes[1];
   fMassLowLambda[2] = fMassMean[1]+5.*fMassRes[1];
   fMassHighLambda[0] = fMassMean[1]-5.*fMassRes[1];
   fMassHighLambda[1] = fMassMean[1]+3.*fMassRes[1];
   fMassHighLambda[2] = fMassMean[1]+8.*fMassRes[1];
  
   fMassLowAntiLambda[0] = fMassMean[1]-8.*fMassRes[1];
   fMassLowAntiLambda[1] = fMassMean[1]-3.*fMassRes[1];
   fMassLowAntiLambda[2] = fMassMean[1]+5.*fMassRes[1];
   fMassHighAntiLambda[0] = fMassMean[1]-5.*fMassRes[1];
   fMassHighAntiLambda[1] = fMassMean[1]+3.*fMassRes[1];
   fMassHighAntiLambda[2] = fMassMean[1]+8.*fMassRes[1];

   //-----------------------------------------------------------
   // Settings for event mixing 
   const Int_t nCentralityBins  = 9;
   Double_t centBins[] = {0., 10.,20.,30.,40.,50.,60.,70.,80.,90.};
   const Double_t* centralityBins = centBins;
   
  // const Int_t nZvtxBins  = 7;
  // Double_t vertexBins[] = {-7.,-5.,-3.,-1.,1.,3.,5.,7.};
  // const Double_t* zvtxBins = vertexBins;

const Int_t nZvtxBins  =  fNumOfVzBins;//fNumOfVzBins;
    Double_t vertexBins[nZvtxBins+1];
    vertexBins[0]=-1*fPrimaryVertexCut;

    if(nZvtxBins==9) {
        vertexBins[1]=-7.0;
        for(Int_t i=2;i<nZvtxBins;i++){
            vertexBins[i]=vertexBins[i-1]+2;
        }
        vertexBins[9]=fPrimaryVertexCut;
    }

    else{
        Double_t binstep = Double_t(2*fPrimaryVertexCut)/nZvtxBins;
        for(Int_t i=1; i<nZvtxBins+1; i++){
            vertexBins[i]=vertexBins[i-1]+binstep;
        }
    }




   Int_t trackDepth = fMixingTracks;
   Int_t poolSize   = fPoolSize;  // Maximum number of events, ignored in the present implemented of AliEventPoolManager

   fPoolMgr = new AliEventPoolManager(poolSize, trackDepth, nCentralityBins, centBins, nZvtxBins, vertexBins);

  }
 //=====================================================================================================================
  void AliAnalysisTaskV0ChCorrelationsys::AddQAEvent()
{
  TList *tQAEvent;
  tQAEvent = new TList();
  tQAEvent->SetOwner();
  tQAEvent->SetName("EventInput");
  
   const Int_t nCentralityBins  = 9;
   Double_t centBins[] = {0., 10.,20.,30.,40.,50.,60.,70.,80.,90.};
   const Double_t* centralityBins = centBins;

  
//  const Int_t nZvtxBins  = 7;
//  Double_t vertexBins[] = {-7.,-5.,-3.,-1.,1.,3.,5.,7.};
 // const Double_t* zvtxBins = vertexBins;

const Int_t nZvtxBins  =  fNumOfVzBins;//fNumOfVzBins;
    Double_t vertexBins[nZvtxBins+1];
    vertexBins[0]=-1*fPrimaryVertexCut;

    if(nZvtxBins==9) {
        vertexBins[1]=-7.0;
        for(Int_t i=2;i<nZvtxBins;i++){
            vertexBins[i]=vertexBins[i-1]+2;
        }
        vertexBins[9]=fPrimaryVertexCut;
    }

    else{
        Double_t binstep = Double_t(2*fPrimaryVertexCut)/nZvtxBins;
        for(Int_t i=1; i<nZvtxBins+1; i++){
            vertexBins[i]=vertexBins[i-1]+binstep;
        }
    }



  TH1F *fhEventBf = new TH1F("fhEventBf", "Event Number; Counts; Number of Events",90, 0.,90); 
  tQAEvent->Add(fhEventBf);

  TH2F *fHistCentVtx = new TH2F("fHistCentVtx", "Centrality vs. Z vertex", 
                                 nCentralityBins, centralityBins, nZvtxBins, vertexBins);
  tQAEvent->Add(fHistCentVtx);

 // TH1F *fHistMultiMain = new TH1F("fHistMultiMain", "Multiplicity of main events", 2000, 0, 2000);
  TH1F *fHistMultiMain = new TH1F("fHistMultiMain", "Multiplicity of main events", 4, 0, 4);
  tQAEvent->Add(fHistMultiMain);

  TH1F *fhEventCentAfterPilp = new TH1F( "fhEventCentAfterPilp","Event distribution to centrality after pile up remove ; Centrality ; Number of Events",90, 0., 90); 
  tQAEvent->Add(fhEventCentAfterPilp);

  TH1F *fhEventAf = new TH1F("fhEventAf", "Event Number; Counts; Number of Events", 90, 0, 90);
  tQAEvent->Add(fhEventAf);


//---------------------------------------------
TH1D *fHistV0Multiplicity = new TH1D ("fHistV0Multiplicity", "V0 event Multiplicity ", 100, 0, 100);
	tQAEvent->Add(fHistV0Multiplicity);
//----------------------------------------------
    
  fOutput->Add(tQAEvent);
}

//======================================================================================================================
void AliAnalysisTaskV0ChCorrelationsys::AddQATrackCandidates()
{
   TList *tQATrack;
   tQATrack = new TList();
   tQATrack->SetOwner();
   tQATrack->SetName("Track");

   // defining bins for centrality
   const Int_t nCentralityBins  = 9;
   Double_t centBins[] = {0., 10.,20.,30.,40.,50.,60.,70.,80.,90.};
   const Double_t* centralityBins = centBins;
   

// defining bins for Z vertex
    
 // const Int_t nZvtxBins  = 7;
  // Double_t vertexBins[] = {-7.,-5.,-3.,-1.,1.,3.,5.,7.};
  // const Double_t* zvtxBins = vertexBins;

 const Int_t nZvtxBins  =  fNumOfVzBins;//fNumOfVzBins;
    Double_t vertexBins[nZvtxBins+1];
    vertexBins[0]=-1*fPrimaryVertexCut;

    if(nZvtxBins==9) {
        vertexBins[1]=-7.0;
        for(Int_t i=2;i<nZvtxBins;i++){
            vertexBins[i]=vertexBins[i-1]+2;
        }
        vertexBins[9]=fPrimaryVertexCut;
    }

    else{
        Double_t binstep = Double_t(2*fPrimaryVertexCut)/nZvtxBins;
        for(Int_t i=1; i<nZvtxBins+1; i++){
            vertexBins[i]=vertexBins[i-1]+binstep;
        }
    }


       //{nZvtxBins,vertexBins[0],  ,vertexBins[nZvtxBins]

   // pt bins of associate particles for the analysis
    const Int_t nPtBins = 29;
   const Double_t PtBins[30] = {1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6, 6.5, 7, 8, 9, 10}; 
   //defining bins of Eta distribution
   const Int_t nEtaBins=16; 
   Double_t EtaBins[nEtaBins+1] = {0.};
   EtaBins[0] = -0.8;
   for (Int_t i=0; i<nEtaBins; i++) {
     EtaBins[i+1] = EtaBins[i] + 2.*0.8/nEtaBins; }

   //defining bins of Phi distribution
   const Int_t nPhiBins = 36; 
   Double_t PhiBins[nPhiBins+1] = {0.};
   PhiBins[0] = 0;
   for (Int_t i=0; i<nPhiBins; i++) { 
     PhiBins[i+1] = PhiBins[i] + 2.*TMath::Pi()/nPhiBins; }
   
   TH1F *tnTracksBf
     = new TH1F("nTracksBf","the number of tracks before set cuts; tracks; Counts", 2000, 0, 2000);
   tQATrack->Add(tnTracksBf);

   TH1F *tnTracksAf
      = new TH1F("nTracksAf","the number of Tracks After set cuts ;  Tracks ; Counts", 2000, 0, 2000);
   tQATrack->Add(tnTracksAf);
   
   TH1F *ttPhi
     = new TH1F("tPhi","the phi of selected track;phi;Counts", 2000,  0, 2*TMath::Pi() );
   tQATrack->Add(ttPhi);

   TH1F *ttPt
     = new TH1F("tPt","the pt of selected track;pt;Counts", 14, 1, 8);
   tQATrack->Add(ttPt);

   TH1F *ttEta
     = new TH1F("tEta","the eta of selected track;eta;Counts", 64, -1.4, 1.4);
   tQATrack->Add(ttEta);

   TH3F *fHistTrk = new TH3F("fHistTrk","pt, eta and phi of reconstructed tracks; p_{T} [GeV/c]; #eta; #phi",
                             14., 1., 8., 16, -0.8, 0.8, 36, 0, 2*TMath::Pi());
   tQATrack->Add(fHistTrk);

   const Int_t spBinsTrack[4] = {nPtBins,  nEtaBins, nCentralityBins, nZvtxBins };
   const Double_t spMinTrack[4] = {PtBins[0],  EtaBins[0], centralityBins[0], -10};
   const Double_t spMaxTrack[4] = {PtBins[nPtBins], EtaBins[nEtaBins], centralityBins[nCentralityBins], 10};


   THnSparseF *fHistMassTrack = new THnSparseF("fHistMassTrack","Mass for Track hypothesis", 4, spBinsTrack, spMinTrack, spMaxTrack);
               fHistMassTrack->GetAxis(0)->Set(nPtBins, PtBins); 
               fHistMassTrack->GetAxis(3)->Set(nZvtxBins,vertexBins);
   tQATrack->Add(fHistMassTrack);
   fHistMassTrack->Sumw2(); 

   fOutput3->Add(tQATrack);
}
//===========================================================================================================================
void AliAnalysisTaskV0ChCorrelationsys::AddQAV0Candidates()
{
  TList *tQAV0;
  tQAV0 = new TList();
  tQAV0->SetOwner();
  tQAV0->SetName("V0");

  TH1F *tlRapK0s
     = new TH1F("BflRapK0s","the rapidity of K0s;y;Counts",50,0,20);
   tQAV0->Add(tlRapK0s);

   TH1F *tlRapLambda
     = new TH1F("BflRapLambda","the rapidity of Lambda;y;Counts",50,0,20);
   tQAV0->Add(tlRapLambda);

   TH1F *tlRapAntiLambda
     = new TH1F("BflRapAntiLambda","the rapidity of AntiLambda;y;Counts",50,0,20);
   tQAV0->Add(tlRapAntiLambda);

   TH1F *tBfV0Pt
     = new TH1F("BfV0Pt","V0's Pt distribution;Pt;Counts", 30, 0, 15);
   tQAV0->Add(tBfV0Pt);
  
   TH1F *tAfV0Pt
     = new TH1F("AfV0Pt","V0's Pt distribution;Pt;Counts", 30, 0, 15);
   tQAV0->Add(tAfV0Pt);

   TH1F *tdPhi
     = new TH1F("BfV0Phi","V0's Phi distribution;Phi;Counts", 36, 0., 2*TMath::Pi());
   tQAV0->Add(tdPhi);

   TH1F *tBfV0Eta
     = new TH1F("BfV0Eta","V0's Eta distribution;Eta;Counts", 16, -0.8, 0.8);
   tQAV0->Add(tBfV0Eta);

   TH1F *tAfV0Eta
     = new TH1F("AfV0Eta","V0's Eta distribution;Eta;Counts", 16, -0.8, 0.8);
   tQAV0->Add(tAfV0Eta);

   TH1F *tBfxyn
     = new TH1F("Bfxyn"," DCA of daughter negative track to Primary Vertex;xyn;Counts", 50, -1, 1);
   tQAV0->Add(tBfxyn);
   
   TH1F *tBfxyp
     = new TH1F("Bfxyp"," DCA of daughter positive track to Primary Vertex;xyn;Counts", 50, -1, 1);
   tQAV0->Add(tBfxyp);

   TH1F *tAfxyn
     = new TH1F("Afxyn"," DCA of daughter negative track to Primary Vertex;xyn;Counts", 50, -1, 1);
   tQAV0->Add(tAfxyn);
   
   TH1F *tAfxyp
     = new TH1F("Afxyp"," DCA of daughter positive track to Primary Vertex;xyn;Counts", 50, -1, 1);
   tQAV0->Add(tAfxyp);

   TH2F *tBfD0vsPt
     = new TH2F("BfD0vsPt",  "D0;[cm];Pt [GeV]",  50, -1, +1, 30, 0, 15);
   tQAV0->Add(tBfD0vsPt);

   TH1F *tBfDCA
     = new TH1F("BfDCA","DCA of daughter track;DCA;Counts",50, 0, 1.001);
   tQAV0->Add(tBfDCA);

   TH1F *tAfDCA
     = new TH1F("AfDCA","DCA of daughter track;DCA;Counts",50, 0, 1.001);
   tQAV0->Add(tAfDCA);

   TH2F *tDCA  
     = new TH2F("BfDCAvsPt", "DCA;[cm];Pt [GeV]", 50, 0, 1.001, 30, 0, 15);
   tQAV0->Add(tDCA);
   
   TH1F *tBfCPA
     = new TH1F("BfCPA","Cosinus of pointing angle;CPA;Counts", 80, 0.91, 1.001);
   tQAV0->Add(tBfCPA);

   TH1F *tAfCPA
     = new TH1F("AfCPA","Cosinus of pointing angle;CPA;Counts", 80, 0.91, 1.001);
   tQAV0->Add(tAfCPA);

   TH2F *tBfCPAvsPt 
     = new TH2F("BfCPAvsPt", "CPA;;Pt [GeV]", 80, 0.91, 1.001, 30, 0, 15);
   tQAV0->Add(tBfCPAvsPt);

   TH1F *tBfr2
     = new TH1F("Bfr2","r2;;Counts", 20000, 0, 20000); 
   tQAV0->Add(tBfr2);

   TH1F *tAfr2
     = new TH1F("Afr2","r2;;Counts", 20000, 0, 20000); 
   tQAV0->Add(tAfr2);
   
   TH1F *tBfDL
     =new TH1F("BfDL","Decay length;DL;Counts", 500, 0, 100);
   tQAV0->Add(tBfDL);
 
   TH2F *tBfDLvsPt
     = new TH2F("BfDLvsPt","DL;[cm];Pt [GeV]",  500, 0, 100,  30, 0, 15);
   tQAV0->Add(tBfDLvsPt);

   TH1F *tDLK
     =new TH1F("DLK","Decay length of K0s;DLK;Counts", 500, 0, 100);
   tQAV0->Add(tDLK);

   TH1F *tDLL
     =new TH1F("DLL","Decay length of Lambda;DLL;Counts", 500, 0, 100);
   tQAV0->Add(tDLL);

   TH1F *tBfDCA2PV
     = new TH1F("BfDCA2PV","DCA to Prim. Vertex;[cm];Counts", 80, 0, 4.);
   tQAV0->Add(tBfDCA2PV);

   TH1F *tAfDCA2PV
     = new TH1F("AfDCA2PV","DCA to Prim. Vertex;[cm];Counts", 80, 0, 4.);
   tQAV0->Add(tAfDCA2PV);

   TH2F *tBefDCA2PV
     = new TH2F("BfDCA2PVvsPt","DCA to Prim. Vertex with Pt; DCA2PV[cm]; Pt [GeV]", 80, 0, 4., 30, 0, 15);
   tQAV0->Add(tBefDCA2PV);

   TH3F *tBfAPvsPt 
     = new TH3F("BfAPvsPt","AP;#alpha;q_{t}[GeV];Pt [GeV]", 80, -1, +1, 90, 0, 0.3, 30, 0, 15);
   tQAV0->Add(tBfAPvsPt);

   TH2F *tBfAP 
     = new TH2F("BfAP","AP;#alpha;q_{t}[GeV]", 80, -1, +1, 90, 0, 0.3);
   tQAV0->Add(tBfAP);

   TH2F *tAfAP 
     = new TH2F("AfAP","AP;#alpha;q_{t}[GeV]", 80, -1, +1, 90, 0, 0.3);
   tQAV0->Add(tAfAP);

   TH2F *fh2TPCdEdxOfPion
    = new TH2F( "TPCdEdxOfPion",
                "TPC dE/dx of Pion; Pt [GeV/c]; TPC signal (ADC) ",
                2000, -10, 10.0, 450, 0., 900.);
   tQAV0->Add(fh2TPCdEdxOfPion);

   TH2F *fh2TPCdEdxOfProton = new TH2F("TPCdEdxOfProton",
                                       "TPC dE/dx of Proton; Pt [GeV/c]; TPC signal (ADC) ",
                                       2000, -10.0, 10.0, 450, 0., 900.);
   tQAV0->Add(fh2TPCdEdxOfProton);

     
   fOutput4->Add(tQAV0);
}

//====================================================================================================================
void AliAnalysisTaskV0ChCorrelationsys::AddQAAnalysisK0s()
{
   TList *tQAK0s;
   tQAK0s = new TList();
   tQAK0s->SetOwner();
   tQAK0s->SetName("K0s");
   
   // defining bins for centrality
   const Int_t nCentralityBins  = 9;
   Double_t centBins[] = {0., 10.,20.,30.,40.,50.,60.,70.,80.,90.};
   const Double_t* centralityBins = centBins;
   
  
 //  const Int_t nZvtxBins  = 7;
 // Double_t vertexBins[] = {-7.,-5.,-3.,-1.,1.,3.,5.,7.};
 // const Double_t* zvtxBins = vertexBins;
   
   const Int_t nZvtxBins  =  fNumOfVzBins;//fNumOfVzBins;
    Double_t vertexBins[nZvtxBins+1];
    vertexBins[0]=-1*fPrimaryVertexCut;

    if(nZvtxBins==9) {
        vertexBins[1]=-7.0;
        for(Int_t i=2;i<nZvtxBins;i++){
            vertexBins[i]=vertexBins[i-1]+2;
        }
        vertexBins[9]=fPrimaryVertexCut;
    }

    else{
        Double_t binstep = Double_t(2*fPrimaryVertexCut)/nZvtxBins;
        for(Int_t i=1; i<nZvtxBins+1; i++){
            vertexBins[i]=vertexBins[i-1]+binstep;
        }
    }

 
    const Int_t nPtBinsV0Xi = 3;
    const Double_t PtBinsV0Xi[4] = {3.0,4.0,8.0,16.0}; 
   
      
   // pt bins of associate particles for the analysis
   const Int_t nPtBins = 6;
   const Double_t PtBins[7] = {1.0,2.0,3.0,4.0,6.0,8.0,10.0}; 

  
   // cascade trigger class: 1 - signal (mass peak region), 2 - left mass sideband, 3 - right mass sideband
   const Int_t nTrigC = 3;
   const Double_t TrigC[4] = {0.5, 1.5, 2.5, 3.5};
  
   // defining bins for mass distributions
   Int_t nBins = 160;
   Double_t mMassMin[] ={0.40, 1.07};
   Double_t mMassMax[] ={0.58, 1.15};

   const Int_t spBinsK0s[4] = {nBins, nPtBinsV0Xi, nCentralityBins, nZvtxBins };
   const Double_t spMinK0s[4] = {mMassMin[0], PtBinsV0Xi[0], centralityBins[0], -10};
   const Double_t spMaxK0s[4] = {mMassMax[0], PtBinsV0Xi[nPtBinsV0Xi],centralityBins[nCentralityBins], 10};


   THnSparseF *fHistMassK0s = new THnSparseF("fHistMassK0s","Mass for K0s hypothesis", 4, spBinsK0s, spMinK0s, spMaxK0s);
               fHistMassK0s->GetAxis(1)->Set(nPtBinsV0Xi, PtBinsV0Xi); 
   tQAK0s->Add(fHistMassK0s);
   fHistMassK0s->Sumw2(); 

   //defining bins of Eta distribution
   const Int_t nEtaBins=14; 
   Double_t EtaBins[nEtaBins+1] = {0.};
   EtaBins[0] = -0.7;
   for (Int_t i=0; i<nEtaBins; i++) {
     EtaBins[i+1] = EtaBins[i] + 2.*0.7/nEtaBins; }

   //defining bins of Phi distribution
   const Int_t nPhiBins = 36; 
   Double_t PhiBins[nPhiBins+1] = {0.};
   PhiBins[0] = 0;
   for (Int_t i=0; i<nPhiBins; i++) { 
     PhiBins[i+1] = PhiBins[i] + 2.*TMath::Pi()/nPhiBins; }
   
   TH3F *fHistK0s = new TH3F("fHistK0s",
                             "pt, eta and phi of reconstructed K0s; p_{T} [GeV/c]; #eta; #phi",
                             nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins, nPhiBins, PhiBins);
   tQAK0s->Add(fHistK0s);

   TH1F *fHistK0sPhi = new TH1F("fHistK0sPhi",
                                "phi of reconstructed K0s;  #phi; Counts",
                                 36, 0, 2.*TMath::Pi());
   tQAK0s->Add(fHistK0sPhi);
   TH1F *fHistK0sEta = new TH1F("fHistK0sEta",
                                "Eta of reconstructed K0s;  #eta; Counts",
                                 14, -0.7, 0.7);
   tQAK0s->Add(fHistK0sEta);

   // defining bins for dPhi distributions
   const Int_t nbPhiBins = 32;//mustafa change
   const Double_t kPi = TMath::Pi();
   Double_t PhiMin = -kPi/2.;
   Double_t PhiMax = -kPi/2. + 2*kPi;
   Double_t Phibins[nbPhiBins+1] = {0.};
   Phibins[0] = PhiMin;
   for (Int_t i=0; i<nbPhiBins; i++) { Phibins[i+1] = Phibins[i] + (PhiMax - PhiMin)/nbPhiBins; }

   // defining bins for dEta distributions
   const Int_t nbEtaBins = 32;//mustafa change
   Double_t EtaMin = -1.5;
   Double_t EtaMax = 1.5;
   Double_t Etabins[nbEtaBins+1] = {0.};
   Etabins[0] = EtaMin;
   for (Int_t i=0; i<nbEtaBins; i++) { Etabins[i+1] = Etabins[i] + (EtaMax - EtaMin)/nbEtaBins; }

  //Create correlation histograms
  //0-dPhi, 1-dEta, 2-Pt trigger, 3-Pt associate, 4-Centrality, 5-Zvertex, 6-Trigger class
   const Int_t corBinsK0s[7] = {nbPhiBins, nbEtaBins, nPtBinsV0Xi, nPtBins,nCentralityBins, nZvtxBins, nTrigC};
   const Double_t corMinK0s[7] = {Phibins[0], Etabins[0], PtBinsV0Xi[0], PtBins[0], centralityBins[0], -10, TrigC[0]};
   const Double_t corMaxK0s[7] = {Phibins[nbPhiBins], Etabins[nbEtaBins], PtBinsV0Xi[nPtBinsV0Xi], PtBins[nPtBins],
                                  centralityBins[nCentralityBins], 10, TrigC[nTrigC]};
   const Int_t corBinsMixK0s[7] = {nbPhiBins, nbEtaBins, nPtBinsV0Xi, nPtBins, nCentralityBins, nZvtxBins, nTrigC};
   const Double_t corMinMixK0s[7] = {Phibins[0], Etabins[0], PtBinsV0Xi[0], PtBins[0], centralityBins[0], -10, TrigC[0]};
   const Double_t corMaxMixK0s[7] = {Phibins[nbPhiBins], Etabins[nbEtaBins], PtBinsV0Xi[nPtBinsV0Xi], PtBins[nPtBins], centralityBins[nCentralityBins], 10, TrigC[nTrigC]}; 


// Gen  
if(fAnalysisMC){
 
 THnSparseF *fHistGendPhidEtaSibK0s = new THnSparseF("fHistGendPhidEtaSibK0s","dPhi vs. dEta siblings", 7, corBinsK0s, corMinK0s, corMaxK0s);
             fHistGendPhidEtaSibK0s->GetAxis(2)->Set(nPtBinsV0Xi, PtBinsV0Xi); 
             fHistGendPhidEtaSibK0s->GetAxis(3)->Set(nPtBins, PtBins); 
             fHistGendPhidEtaSibK0s->GetAxis(5)->Set(nZvtxBins, vertexBins); 
   tQAK0s->Add(fHistGendPhidEtaSibK0s);
   fHistGendPhidEtaSibK0s->Sumw2();

   
   THnSparseF *fHistGendPhidEtaMixK0s = new THnSparseF("fHistGendPhidEtaMixK0s","dPhi vs. dEta mixed", 7, corBinsMixK0s, corMinMixK0s, corMaxMixK0s);
               fHistGendPhidEtaMixK0s->GetAxis(2)->Set(nPtBinsV0Xi, PtBinsV0Xi); 
               fHistGendPhidEtaMixK0s->GetAxis(3)->Set(nPtBins, PtBins); 
               fHistGendPhidEtaMixK0s->GetAxis(5)->Set(nZvtxBins, vertexBins); 

   tQAK0s->Add(fHistGendPhidEtaMixK0s);
   fHistGendPhidEtaMixK0s->Sumw2(); 
}
   
   THnSparseF *fHistdPhidEtaSibK0s = new THnSparseF("fHistdPhidEtaSibK0s","dPhi vs. dEta siblings", 7, corBinsK0s, corMinK0s, corMaxK0s);
               fHistdPhidEtaSibK0s->GetAxis(2)->Set(nPtBinsV0Xi, PtBinsV0Xi); 
               fHistdPhidEtaSibK0s->GetAxis(3)->Set(nPtBins, PtBins); 
               fHistdPhidEtaSibK0s->GetAxis(5)->Set(nZvtxBins, vertexBins); 

   tQAK0s->Add(fHistdPhidEtaSibK0s);
   fHistdPhidEtaSibK0s->Sumw2();

   THnSparseF *fHistdPhidEtaMixK0s = new THnSparseF("fHistdPhidEtaMixK0s","dPhi vs. dEta mixed", 7, corBinsMixK0s, corMinMixK0s, corMaxMixK0s);
               fHistdPhidEtaMixK0s->GetAxis(2)->Set(nPtBinsV0Xi, PtBinsV0Xi); 
               fHistdPhidEtaMixK0s->GetAxis(3)->Set(nPtBins, PtBins); 
               fHistdPhidEtaMixK0s->GetAxis(5)->Set(nZvtxBins, vertexBins); 

   tQAK0s->Add(fHistdPhidEtaMixK0s);
   fHistdPhidEtaMixK0s->Sumw2(); 

   //pt trigger K0s including non-correlated   
   const Int_t trigAllBinsK0s[4] = {nPtBinsV0Xi, nCentralityBins,nZvtxBins, nTrigC};
   const Double_t trigAllMinK0s[4] = {PtBinsV0Xi[0], centBins[0],-10, TrigC[0]};
   const Double_t trigAllMaxK0s[4] = {PtBinsV0Xi[nPtBinsV0Xi], centBins[nCentralityBins],10, TrigC[nTrigC]};

// Gen 
if(fAnalysisMC){
  
THnSparseF *fHistGenTrigSibAllK0s = new THnSparseF("fHistGenTrigSibAllK0s","pt trigger K0s including non-correlated",4, trigAllBinsK0s, trigAllMinK0s, trigAllMaxK0s);
            fHistGenTrigSibAllK0s->GetAxis(0)->Set(nPtBinsV0Xi, PtBinsV0Xi); 
            fHistGenTrigSibAllK0s->GetAxis(2)->Set(nZvtxBins, vertexBins); 
   tQAK0s->Add(fHistGenTrigSibAllK0s);
   fHistGenTrigSibAllK0s->Sumw2();
}
//Rec 
   THnSparseF *fHistTrigSibAllK0s = new THnSparseF("fHistTrigSibAllK0s","pt trigger K0s including non-correlated", 4, trigAllBinsK0s, trigAllMinK0s, trigAllMaxK0s);
               fHistTrigSibAllK0s->GetAxis(0)->Set(nPtBinsV0Xi, PtBinsV0Xi);
               fHistTrigSibAllK0s->GetAxis(2)->Set(nZvtxBins, vertexBins);  
   tQAK0s->Add(fHistTrigSibAllK0s);
   fHistTrigSibAllK0s->Sumw2();



// defining bins for mass distributions
  Int_t nBinsK0s = 160;
 Double_t mMassK0sMin = 0.40;
   Double_t mMassK0sMax = 0.58;

// pt bins of (K0s) 
  const Int_t nPtBinsV0Sp = 13;
  const Double_t PtBinsV0Sp[14] ={ 3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6.5, 8, 10, 12, 16}; 


 TH2F *fHistK0sMassvsPtCorr = new TH2F("fHistK0sMassvsPtCorr", "K^{0}_{s}  Mass vs p_{T} ;p_{T} (GeV/c);Inv.Mass (GeV/c^{2})", nPtBinsV0Sp, PtBinsV0Sp, nBinsK0s, mMassK0sMin ,mMassK0sMax );

   tQAK0s->Add(fHistK0sMassvsPtCorr);
  fHistK0sMassvsPtCorr->Sumw2();  

 TH2F *fHistK0sMassvsPtNoCorr = new TH2F("fHistK0sMassvsPtNoCorr", "K^{0}_{s}  Mass vs p_{T} ;p_{T} (GeV/c);Inv.Mass (GeV/c^{2})", nPtBinsV0Sp, PtBinsV0Sp, nBinsK0s, mMassK0sMin ,mMassK0sMax );

   tQAK0s->Add(fHistK0sMassvsPtNoCorr);
  fHistK0sMassvsPtNoCorr->Sumw2();   



   fOutput5->Add(tQAK0s);
}

//====================================================================================================================
void AliAnalysisTaskV0ChCorrelationsys::AddQAAnalysisLambda()
{
   TList *tQALambda;
   tQALambda = new TList();
   tQALambda->SetOwner();
   tQALambda->SetName("Lambda");
   
   // defining bins for centrality
   const Int_t nCentralityBins  = 9;
   Double_t centBins[] = {0., 10.,20.,30.,40.,50.,60.,70.,80.,90.};
   const Double_t* centralityBins = centBins;
  
  // const Int_t nZvtxBins  = 7;
  // Double_t vertexBins[] = {-7.,-5.,-3.,-1.,1.,3.,5.,7.};
  // const Double_t* zvtxBins = vertexBins;
   
  

 const Int_t nZvtxBins  =  fNumOfVzBins;//fNumOfVzBins;
    Double_t vertexBins[nZvtxBins+1];
    vertexBins[0]=-1*fPrimaryVertexCut;

    if(nZvtxBins==9) {
        vertexBins[1]=-7.0;
        for(Int_t i=2;i<nZvtxBins;i++){
            vertexBins[i]=vertexBins[i-1]+2;
        }
        vertexBins[9]=fPrimaryVertexCut;
    }

    else{
        Double_t binstep = Double_t(2*fPrimaryVertexCut)/nZvtxBins;
        for(Int_t i=1; i<nZvtxBins+1; i++){
            vertexBins[i]=vertexBins[i-1]+binstep;
        }
    }

   
    const Int_t nPtBinsV0Xi = 3;
    const Double_t PtBinsV0Xi[4] = {3.0,4.0,8.0,16.0}; 
   

   // pt bins of associate particles for the analysis
   const Int_t nPtBins = 6;
   const Double_t PtBins[7] = {1.0,2.0,3.0,4.0,6.0,8.0,10.0}; 

   // cascade trigger class: 1 - signal (mass peak region), 2 - left mass sideband, 3 - right mass sideband
   const Int_t nTrigC = 3;
   const Double_t TrigC[4] = {0.5, 1.5, 2.5, 3.5};
  
   // defining bins for mass distributions
   Int_t nBins = 160;
   Double_t mMassMin[] ={0.40, 1.07};
   Double_t mMassMax[] ={0.58, 1.15};

   const Int_t spBinsLambda[4] = {nBins, nPtBinsV0Xi, nCentralityBins, nZvtxBins };
   const Double_t spMinLambda[4] = {mMassMin[1], PtBinsV0Xi[0], centralityBins[0], -10};
   const Double_t spMaxLambda[4] = {mMassMax[1], PtBinsV0Xi[nPtBinsV0Xi], centralityBins[nCentralityBins], 10};




   THnSparseF *fHistMassLambda = new THnSparseF("fHistMassLambda","Mass for Lambda hypothesis", 4, spBinsLambda, spMinLambda, spMaxLambda);
               fHistMassLambda->GetAxis(1)->Set(nPtBinsV0Xi, PtBinsV0Xi); 
   tQALambda->Add(fHistMassLambda);
   fHistMassLambda->Sumw2();

  //defining bins of Eta distribution
   const Int_t nEtaBins=14; 
   Double_t EtaBins[nEtaBins+1] = {0.};
   EtaBins[0] = -0.7;
   for (Int_t i=0; i<nEtaBins; i++) {
     EtaBins[i+1] = EtaBins[i] + 2.*0.7/nEtaBins; }

   //defining bins of Phi distribution
   const Int_t nPhiBins = 36; 
   Double_t PhiBins[nPhiBins+1] = {0.};
   PhiBins[0] = 0;
   for (Int_t i=0; i<nPhiBins; i++) { 
     PhiBins[i+1] = PhiBins[i] + 2.*TMath::Pi()/nPhiBins; }

   TH3F *fHistLambda = new TH3F("fHistLambda",
                          "pt, eta and phi of reconstructed Lambda; p_{T} [GeV/c]; #eta; #phi",nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins, nPhiBins, PhiBins);
   tQALambda->Add(fHistLambda);
  
   TH1F *fHistLambdaPhi = new TH1F("fHistLambdaPhi", "phi of reconstructed Lambda;  #phi; Counts", 36, 0, 2.*TMath::Pi());
   tQALambda->Add(fHistLambdaPhi);
   TH1F *fHistLambdaEta = new TH1F("fHistLambdaEta","Eta of reconstructed Lambda;  #eta; Counts", 14, -0.7, 0.7);
   tQALambda->Add(fHistLambdaEta);
 
   // defining bins for dPhi distributions
   const Int_t nbPhiBins = 32;
   const Double_t kPi = TMath::Pi();
   Double_t PhiMin = -kPi/2.;
   Double_t PhiMax = -kPi/2. + 2*kPi;
   Double_t Phibins[nbPhiBins+1] = {0.};
   Phibins[0] = PhiMin;
   for (Int_t i=0; i<nbPhiBins; i++) { Phibins[i+1] = Phibins[i] + (PhiMax - PhiMin)/nbPhiBins; }

   // defining bins for dEta distributions
   const Int_t nbEtaBins = 32;
   Double_t EtaMin = -1.5; 
   Double_t EtaMax = 1.5; 
   Double_t Etabins[nbEtaBins+1] = {0.};
   Etabins[0] = EtaMin;
   for (Int_t i=0; i<nbEtaBins; i++) { Etabins[i+1] = Etabins[i] + (EtaMax - EtaMin)/nbEtaBins; }

   const Int_t corBinsLambda[7] = {nbPhiBins, nbEtaBins, nPtBinsV0Xi, nPtBins, nCentralityBins, nZvtxBins, nTrigC};
   const Double_t corMinLambda[7] = {Phibins[0], Etabins[0],PtBinsV0Xi[0], PtBins[0], centralityBins[0], -10, TrigC[0]};
   const Double_t corMaxLambda[7] = {Phibins[nbPhiBins], Etabins[nbEtaBins],PtBinsV0Xi[nPtBinsV0Xi], PtBins[nPtBins], centralityBins[nCentralityBins],
                                     10, TrigC[nTrigC]};

   const Int_t corBinsMixLambda[7] = {nbPhiBins, nbEtaBins, nPtBinsV0Xi, nPtBins, nCentralityBins, nZvtxBins, nTrigC};
   const Double_t corMinMixLambda[7] = {Phibins[0], Etabins[0],PtBinsV0Xi[0], PtBins[0], centralityBins[0], -10, TrigC[0]};
   const Double_t corMaxMixLambda[7] = {Phibins[nbPhiBins], Etabins[nbEtaBins], PtBinsV0Xi[nPtBinsV0Xi], PtBins[nPtBins],centralityBins[nCentralityBins],
                                        10, TrigC[nTrigC]}; 

   TH3F *fHistLambdaDphiDCAPtSig = new TH3F("fHistLambdaDphiDCAPtSig", ";#Delta#phi; p_{T} [GeV/c]; DCA [cm]",
                                            nbPhiBins, PhiMin, PhiMax,  nPtBinsV0Xi, PtBinsV0Xi[0], PtBinsV0Xi[nPtBinsV0Xi],
                                            80, 0., 4.);
   tQALambda->Add(fHistLambdaDphiDCAPtSig);
   TH3F *fHistLambdaDphiDCAPtBkgL = new TH3F("fHistLambdaDphiDCAPtBkgL",";#Delta#phi; p_{T} [GeV/c]; DCA [cm]", nbPhiBins, PhiMin, PhiMax,
                                             nPtBinsV0Xi, PtBinsV0Xi[0], PtBinsV0Xi[nPtBinsV0Xi],
                                             80, 0., 4.);
   tQALambda->Add(fHistLambdaDphiDCAPtBkgL);TH3F *fHistLambdaDphiDCAPtBkgR = new TH3F("fHistLambdaDphiDCAPtBkgR", ";#Delta#phi; p_{T} [GeV/c]; DCA [cm]",
                                             nbPhiBins, PhiMin, PhiMax, nPtBinsV0Xi, PtBinsV0Xi[0], PtBinsV0Xi[nPtBinsV0Xi],
                                             80, 0., 4.);
   tQALambda->Add(fHistLambdaDphiDCAPtBkgR);

   // Create correlation histograms


//  Gen 

if(fAnalysisMC){
  
THnSparseF *fHistGendPhidEtaSibLambda = new THnSparseF("fHistGendPhidEtaSibLambda","dPhi vs. dEta mixed", 7, corBinsMixLambda, corMinMixLambda, corMaxMixLambda);
            fHistGendPhidEtaSibLambda->GetAxis(2)->Set(nPtBinsV0Xi, PtBinsV0Xi); 
            fHistGendPhidEtaSibLambda->GetAxis(3)->Set(nPtBins, PtBins); 
             fHistGendPhidEtaSibLambda->GetAxis(5)->Set(nZvtxBins, vertexBins); 

        
   tQALambda->Add(fHistGendPhidEtaSibLambda);
   fHistGendPhidEtaSibLambda->Sumw2();

THnSparseF *fHistGendPhidEtaMixLambda = new THnSparseF("fHistGendPhidEtaMixLambda","dPhi vs. dEta mixed", 7, corBinsMixLambda, corMinMixLambda, corMaxMixLambda);
            fHistGendPhidEtaMixLambda->GetAxis(2)->Set(nPtBinsV0Xi, PtBinsV0Xi); 
            fHistGendPhidEtaMixLambda->GetAxis(3)->Set(nPtBins, PtBins);
             fHistGendPhidEtaMixLambda->GetAxis(5)->Set(nZvtxBins, vertexBins);  
   tQALambda->Add(fHistGendPhidEtaMixLambda);
   fHistGendPhidEtaMixLambda->Sumw2();

}
//  Rec

   THnSparseF *fHistdPhidEtaSibLambda = new THnSparseF("fHistdPhidEtaSibLambda","dPhi vs. dEta siblings", 7, corBinsLambda, corMinLambda, corMaxLambda);
               fHistdPhidEtaSibLambda->GetAxis(2)->Set(nPtBinsV0Xi, PtBinsV0Xi); 
               fHistdPhidEtaSibLambda->GetAxis(3)->Set(nPtBins, PtBins); 
               fHistdPhidEtaSibLambda->GetAxis(5)->Set(nZvtxBins, vertexBins);  
   tQALambda->Add(fHistdPhidEtaSibLambda);
   fHistdPhidEtaSibLambda->Sumw2();

   THnSparseF *fHistdPhidEtaMixLambda = new THnSparseF("fHistdPhidEtaMixLambda","dPhi vs. dEta mixed", 7, corBinsMixLambda, corMinMixLambda, corMaxMixLambda);
               fHistdPhidEtaMixLambda->GetAxis(2)->Set(nPtBinsV0Xi, PtBinsV0Xi); 
               fHistdPhidEtaMixLambda->GetAxis(3)->Set(nPtBins, PtBins); 
               fHistdPhidEtaMixLambda->GetAxis(5)->Set(nZvtxBins, vertexBins);  
   tQALambda->Add(fHistdPhidEtaMixLambda);
   fHistdPhidEtaMixLambda->Sumw2();

   //pt trigger Lambda including non-correlated
   const Int_t trigAllBinsLambda[4] = {nPtBinsV0Xi, nCentralityBins,nZvtxBins, nTrigC};
   const Double_t trigAllMinLambda[4] = {PtBinsV0Xi[0], centBins[0], -10, TrigC[0]};
   const Double_t trigAllMaxLambda[4] = {PtBinsV0Xi[nPtBinsV0Xi], centBins[nCentralityBins], 10, TrigC[nTrigC]};

if(fAnalysisMC){
  
   THnSparseF *fHistGenTrigSibAllLambda = new THnSparseF("fHistGenTrigSibAllLambda","pt trigger Lambda including non-correlated", 4, trigAllBinsLambda, trigAllMinLambda, trigAllMaxLambda);
                               fHistGenTrigSibAllLambda->GetAxis(0)->Set(nPtBinsV0Xi, PtBinsV0Xi); 
                               fHistGenTrigSibAllLambda->GetAxis(2)->Set(nZvtxBins, vertexBins);  
   tQALambda->Add(fHistGenTrigSibAllLambda);
   fHistGenTrigSibAllLambda->Sumw2();
}


   THnSparseF *fHistTrigSibAllLambda = new THnSparseF("fHistTrigSibAllLambda", "pt trigger Lambda including non-correlated",
                                           4, trigAllBinsLambda, trigAllMinLambda, trigAllMaxLambda);
              fHistTrigSibAllLambda->GetAxis(0)->Set(nPtBinsV0Xi, PtBinsV0Xi);
              fHistTrigSibAllLambda->GetAxis(2)->Set(nZvtxBins, vertexBins);   
   tQALambda->Add(fHistTrigSibAllLambda);
   fHistTrigSibAllLambda->Sumw2();


// defining bins for mass distributions
  Int_t nBinsLam = 160;
 Double_t mMassLamMin = 1.07;
   Double_t mMassLamMax = 1.15;

// pt bins of (Lam) 
  const Int_t nPtBinsV0Sp = 13;
  const Double_t PtBinsV0Sp[14] ={ 3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6.5, 8, 10, 12, 16}; 

 TH2F *fHistLamMassvsPtCorr = new TH2F("fHistLamMassvsPtCorr", "#Lambda  Mass vs p_{T} ;p_{T} (GeV/c);Inv.Mass (GeV/c^{2})", nPtBinsV0Sp, PtBinsV0Sp, nBinsLam, mMassLamMin ,mMassLamMax );

   tQALambda->Add(fHistLamMassvsPtCorr);
  fHistLamMassvsPtCorr->Sumw2();   

 TH2F *fHistLamMassvsPtNoCorr = new TH2F("fHistLamMassvsPtNoCorr", "#Lambda  Mass vs p_{T} ;p_{T} (GeV/c);Inv.Mass (GeV/c^{2})", nPtBinsV0Sp, PtBinsV0Sp, nBinsLam, mMassLamMin ,mMassLamMax );

   tQALambda->Add(fHistLamMassvsPtNoCorr);
  fHistLamMassvsPtNoCorr->Sumw2();       


   fOutput6->Add(tQALambda);
 } 
//====================================================================================================================
void AliAnalysisTaskV0ChCorrelationsys::AddQAAnalysisAntiLambda()
{
   TList *tQAAntiLambda;
   tQAAntiLambda = new TList();
   tQAAntiLambda->SetOwner();
   tQAAntiLambda->SetName("AntiLambda");
   
   // defining bins for centrality
   const Int_t nCentralityBins  = 9;
   Double_t centBins[] = {0., 10.,20.,30.,40.,50.,60.,70.,80.,90.};
   const Double_t* centralityBins = centBins;
   
   
   
  // const Int_t nZvtxBins  = 7;
  // Double_t vertexBins[] = {-7.,-5.,-3.,-1.,1.,3.,5.,7.};
  // const Double_t* zvtxBins = vertexBins;

   
 const Int_t nZvtxBins  =  fNumOfVzBins;//fNumOfVzBins;
    Double_t vertexBins[nZvtxBins+1];
    vertexBins[0]=-1*fPrimaryVertexCut;

    if(nZvtxBins==9) {
        vertexBins[1]=-7.0;
        for(Int_t i=2;i<nZvtxBins;i++){
            vertexBins[i]=vertexBins[i-1]+2;
        }
        vertexBins[9]=fPrimaryVertexCut;
    }

    else{
        Double_t binstep = Double_t(2*fPrimaryVertexCut)/nZvtxBins;
        for(Int_t i=1; i<nZvtxBins+1; i++){
            vertexBins[i]=vertexBins[i-1]+binstep;
        }
    }

    const Int_t nPtBinsV0Xi = 3;
    const Double_t PtBinsV0Xi[4] = {3.0,4.0,8.0,16.0}; 
   

   // pt bins of associate particles for the analysis
   const Int_t nPtBins = 6;
   const Double_t PtBins[7] = {1.0,2.0,3.0,4.0,6.0,8.0,10.0}; 


   // cascade trigger class: 1 - signal (mass peak region), 2 - left mass sideband, 3 - right mass sideband
   const Int_t nTrigC = 3;
   const Double_t TrigC[4] = {0.5, 1.5, 2.5, 3.5};
  
   // defining bins for mass distributions
   Int_t nBins = 160;
   Double_t mMassMin[] ={0.40, 1.07};
   Double_t mMassMax[] ={0.58, 1.15};

   const Int_t spBinsAntiLambda[4] = {nBins, nPtBinsV0Xi, nCentralityBins, nZvtxBins};
   const Double_t spMinAntiLambda[4] = {mMassMin[1], PtBinsV0Xi[0], centralityBins[0], -10};
   const Double_t spMaxAntiLambda[4] = {mMassMax[1], PtBinsV0Xi[nPtBinsV0Xi],  
                                        centralityBins[nCentralityBins], 10};
//Anti mass 

   THnSparseF *fHistMassAntiLambda = new THnSparseF("fHistMassAntiLambda","Mass for AntiLambda hypothesis", 4, spBinsAntiLambda, spMinAntiLambda, spMaxAntiLambda);
              fHistMassAntiLambda->GetAxis(1)->Set(nPtBinsV0Xi, PtBinsV0Xi); 
   tQAAntiLambda->Add(fHistMassAntiLambda);
   fHistMassAntiLambda->Sumw2(); 

   //defining bins of Eta distribution
   const Int_t nEtaBins=14; 
   Double_t EtaBins[nEtaBins+1] = {0.};
   EtaBins[0] = -0.7;
   for (Int_t i=0; i<nEtaBins; i++) {
     EtaBins[i+1] = EtaBins[i] + 2.*0.7/nEtaBins; }

   //defining bins of Phi distribution
   const Int_t nPhiBins = 36; 
   Double_t PhiBins[nPhiBins+1] = {0.};
   PhiBins[0] = 0;
   for (Int_t i=0; i<nPhiBins; i++) { 
     PhiBins[i+1] = PhiBins[i] + 2.*TMath::Pi()/nPhiBins; }
   
   TH3F *fHistAntiLambda = new TH3F("fHistAntiLambda",
                                    "pt, eta and phi of reconstructed AntiLambda; p_{T} [GeV/c]; #eta; #phi", nPtBinsV0Xi, PtBinsV0Xi, nEtaBins, EtaBins, nPhiBins, PhiBins);
   tQAAntiLambda->Add(fHistAntiLambda);
   TH1F *fHistAntiLambdaPhi = new TH1F("fHistAntiLambdaPhi", "phi of reconstructed AntiLambda;  #phi; Counts",  36, 0, 2.*TMath::Pi());
   tQAAntiLambda->Add(fHistAntiLambdaPhi);
   TH1F *fHistAntiLambdaEta = new TH1F("fHistAntiLambdaEta","Eta of reconstructed AntiLambda;  #eta; Counts",
                                        14, -0.7, 0.7);
   tQAAntiLambda->Add(fHistAntiLambdaEta);
  
   //defining bins for dPhi distributions
   const Int_t nbPhiBins = 32;
   const Double_t kPi = TMath::Pi();
   Double_t PhiMin = -kPi/2.;
   Double_t PhiMax = -kPi/2. + 2*kPi;
   Double_t Phibins[nbPhiBins+1] = {0.};
   Phibins[0] = PhiMin;
   for (Int_t i=0; i<nbPhiBins; i++) { Phibins[i+1] = Phibins[i] + (PhiMax - PhiMin)/nbPhiBins; }

   //defining bins for dEta distributions
   const Int_t nbEtaBins = 32;//mustafa change 
   Double_t EtaMin = -1.5;
   Double_t EtaMax = 1.5;
   Double_t Etabins[nbEtaBins+1] = {0.};
   Etabins[0] = EtaMin;
   for (Int_t i=0; i<nbEtaBins; i++) { Etabins[i+1] = Etabins[i] + (EtaMax - EtaMin)/nbEtaBins; }

   const Int_t corBinsAntiLambda[7] = {nbPhiBins, nbEtaBins, nPtBinsV0Xi, nPtBins, nCentralityBins, nZvtxBins, nTrigC};
   const Double_t corMinAntiLambda[7] = {Phibins[0], Etabins[0],PtBinsV0Xi[0], PtBins[0],
                                         centralityBins[0], -10, TrigC[0]};
   const Double_t corMaxAntiLambda[7] = {Phibins[nbPhiBins], Etabins[nbEtaBins], PtBinsV0Xi[nPtBinsV0Xi], PtBins[nPtBins],
                                         centralityBins[nCentralityBins], 10,TrigC[nTrigC]};
   const Int_t corBinsMixAntiLambda[7] = {nbPhiBins, nbEtaBins, nPtBinsV0Xi, nPtBins, nCentralityBins, nZvtxBins, nTrigC};

   const Double_t corMinMixAntiLambda[7] = {Phibins[0], Etabins[0],PtBinsV0Xi[0], PtBins[0], centralityBins[0], -10, TrigC[0]};
   const Double_t corMaxMixAntiLambda[7] = {Phibins[nbPhiBins], Etabins[nbEtaBins],PtBinsV0Xi[nPtBinsV0Xi], PtBins[nPtBins],
                                            centralityBins[nCentralityBins],10,TrigC[nTrigC]};  

   TH3F *fHistAntiLambdaDphiDCAPtSig = new TH3F("fHistAntiLambdaDphiDCAPtSig",";#Delta#phi; p_{T} [GeV/c]; DCA [cm]", nbPhiBins, PhiMin, PhiMax,
                                                nPtBinsV0Xi, PtBinsV0Xi[0], PtBinsV0Xi[nPtBinsV0Xi], 80, 0., 4.);
   tQAAntiLambda->Add(fHistAntiLambdaDphiDCAPtSig);
   TH3F *fHistAntiLambdaDphiDCAPtBkgL = new TH3F("fHistAntiLambdaDphiDCAPtBkgL", ";#Delta#phi; p_{T} [GeV/c]; DCA [cm]",nbPhiBins, PhiMin, PhiMax, nPtBinsV0Xi, PtBinsV0Xi[0], PtBinsV0Xi[nPtBinsV0Xi],
                                                80, 0., 4.);
   tQAAntiLambda->Add(fHistAntiLambdaDphiDCAPtBkgL);
   TH3F *fHistAntiLambdaDphiDCAPtBkgR = new TH3F("fHistAntiLambdaDphiDCAPtBkgR",";#Delta#phi; p_{T} [GeV/c]; DCA [cm]", nbPhiBins, PhiMin, PhiMax, nPtBinsV0Xi, PtBinsV0Xi[0], PtBinsV0Xi[nPtBinsV0Xi],
                                                80, 0., 4.);
   tQAAntiLambda->Add(fHistAntiLambdaDphiDCAPtBkgR);

   // Create correlation histograms
// sib Gen same 
if(fAnalysisMC){
 THnSparseF *fHistGendPhidEtaSibAntiLambda = new THnSparseF("fHistGendPhidEtaSibAntiLambda","dPhi vs. dEta siblings", 7, corBinsAntiLambda, corMinAntiLambda, corMaxAntiLambda);
                 
               fHistGendPhidEtaSibAntiLambda->GetAxis(2)->Set(nPtBinsV0Xi, PtBinsV0Xi); 
               fHistGendPhidEtaSibAntiLambda->GetAxis(3)->Set(nPtBins, PtBins);
               fHistGendPhidEtaSibAntiLambda->GetAxis(5)->Set(nZvtxBins, vertexBins); 
      
   tQAAntiLambda->Add(fHistGendPhidEtaSibAntiLambda);
   fHistGendPhidEtaSibAntiLambda->Sumw2();

   THnSparseF *fHistGendPhidEtaMixAntiLambda = new THnSparseF("fHistGendPhidEtaMixAntiLambda","dPhi vs. dEta mixed", 7, corBinsMixAntiLambda, corMinMixAntiLambda, corMaxMixAntiLambda); 
               fHistGendPhidEtaMixAntiLambda->GetAxis(2)->Set(nPtBinsV0Xi, PtBinsV0Xi); 
               fHistGendPhidEtaMixAntiLambda->GetAxis(3)->Set(nPtBins, PtBins);
               fHistGendPhidEtaMixAntiLambda->GetAxis(5)->Set(nZvtxBins, vertexBins); 
   
   tQAAntiLambda->Add(fHistGendPhidEtaMixAntiLambda);
   fHistGendPhidEtaMixAntiLambda->Sumw2();

}


//  Rec  
   THnSparseF *fHistdPhidEtaSibAntiLambda = new THnSparseF("fHistdPhidEtaSibAntiLambda","dPhi vs. dEta siblings", 7, corBinsAntiLambda, corMinAntiLambda, corMaxAntiLambda);
               fHistdPhidEtaSibAntiLambda->GetAxis(2)->Set(nPtBinsV0Xi, PtBinsV0Xi); 
               fHistdPhidEtaSibAntiLambda->GetAxis(3)->Set(nPtBins, PtBins);
               fHistdPhidEtaSibAntiLambda->GetAxis(5)->Set(nZvtxBins, vertexBins); 

   tQAAntiLambda->Add(fHistdPhidEtaSibAntiLambda);
   fHistdPhidEtaSibAntiLambda->Sumw2();

   THnSparseF *fHistdPhidEtaMixAntiLambda = new THnSparseF("fHistdPhidEtaMixAntiLambda","dPhi vs. dEta mixed", 7, corBinsMixAntiLambda, corMinMixAntiLambda, corMaxMixAntiLambda);  
             fHistdPhidEtaMixAntiLambda->GetAxis(2)->Set(nPtBinsV0Xi, PtBinsV0Xi); 
             fHistdPhidEtaMixAntiLambda->GetAxis(3)->Set(nPtBins, PtBins);
             fHistdPhidEtaMixAntiLambda->GetAxis(5)->Set(nZvtxBins, vertexBins); 

   tQAAntiLambda->Add(fHistdPhidEtaMixAntiLambda);
   fHistdPhidEtaMixAntiLambda->Sumw2();

   const Int_t trigAllBinsAntiLambda[4] = {nPtBinsV0Xi, nCentralityBins, nZvtxBins, nTrigC};
   const Double_t trigAllMinAntiLambda[4] = {PtBinsV0Xi[0], centBins[0],-10, TrigC[0]};
   const Double_t trigAllMaxAntiLambda[4] = {PtBinsV0Xi[nPtBinsV0Xi],centBins[nCentralityBins],10, TrigC[nTrigC]};

//Gen
if(fAnalysisMC){
  
THnSparseF *fHistGenTrigSibAllAntiLambda = new THnSparseF("fHistGenTrigSibAllAntiLambda",
                                              "pt trigger AntiLambda including non-correlated",
                                               4, trigAllBinsAntiLambda, trigAllMinAntiLambda, trigAllMaxAntiLambda);
           fHistGenTrigSibAllAntiLambda->GetAxis(0)->Set(nPtBinsV0Xi, PtBinsV0Xi);
           fHistGenTrigSibAllAntiLambda->GetAxis(2)->Set(nZvtxBins, vertexBins);    
   tQAAntiLambda->Add(fHistGenTrigSibAllAntiLambda);
   fHistGenTrigSibAllAntiLambda->Sumw2();


}
//Rec
   THnSparseF *fHistTrigSibAllAntiLambda = new THnSparseF("fHistTrigSibAllAntiLambda",
                                              "pt trigger AntiLambda including non-correlated",
                                               4, trigAllBinsAntiLambda, trigAllMinAntiLambda, trigAllMaxAntiLambda);
              fHistTrigSibAllAntiLambda->GetAxis(0)->Set(nPtBinsV0Xi, PtBinsV0Xi); 
              fHistTrigSibAllAntiLambda->GetAxis(2)->Set(nZvtxBins, vertexBins);    
   tQAAntiLambda->Add(fHistTrigSibAllAntiLambda);
   fHistTrigSibAllAntiLambda->Sumw2();




// defining bins for mass distributions
  Int_t nBinsLam = 160;
 Double_t mMassLamMin = 1.07;
   Double_t mMassLamMax = 1.15;

// pt bins of (Lam) 
  const Int_t nPtBinsV0Sp = 13;
  const Double_t PtBinsV0Sp[14] ={ 3, 3.2, 3.4, 3.6, 3.8, 4, 4.5, 5, 5.5, 6.5, 8, 10, 12, 16}; 

 TH2F *fHistAntiLamMassvsPtCorr = new TH2F("fHistAntiLamMassvsPtCorr", "#bar{#Lambda}  Mass vs p_{T} ;p_{T} (GeV/c);Inv.Mass (GeV/c^{2})", nPtBinsV0Sp, PtBinsV0Sp, nBinsLam, mMassLamMin ,mMassLamMax );

   tQAAntiLambda->Add(fHistAntiLamMassvsPtCorr);
  fHistAntiLamMassvsPtCorr->Sumw2();  

 TH2F *fHistAntiLamMassvsPtNoCorr = new TH2F("fHistAntiLamMassvsPtNoCorr", "#bar{#Lambda}  Mass vs p_{T} ;p_{T} (GeV/c);Inv.Mass (GeV/c^{2})", nPtBinsV0Sp, PtBinsV0Sp, nBinsLam, mMassLamMin ,mMassLamMax );

   tQAAntiLambda->Add(fHistAntiLamMassvsPtNoCorr);
  fHistAntiLamMassvsPtNoCorr->Sumw2();    

   fOutput7->Add(tQAAntiLambda);
}
  
//====================================================================================================================
void AliAnalysisTaskV0ChCorrelationsys::Terminate(Option_t *)
{
   // Draw result to screen, or perform fitting, normalizations
   // Called once at the end of the query

   fOutput = dynamic_cast<TList*>(GetOutputData(1));
   if (!fOutput) { AliError("Could not retrieve TList fOutput"); return; }

   // NEW HISTO should be retrieved from the TList container in the above way,
   // so it is available to draw on a canvas such as below
}

//======================================================================================================================

void AliAnalysisTaskV0ChCorrelationsys::UserExec(Option_t *)

{

    Int_t nV0 =0;

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *inEvMain = dynamic_cast<AliInputEventHandler *>(mgr->GetInputEventHandler());

  
    ((TH1F*)((TList*)fOutput->FindObject("EventInput"))->FindObject("fhEventBf"))->Fill(0);
    ((TH1F*)((TList*)fOutput->FindObject("EventInput"))->FindObject("fHistMultiMain"))->Fill(0.5);;

    // 2015 physics selection
    UInt_t maskIsSelected= inEvMain->IsEventSelected();
    Bool_t isINT7selected = maskIsSelected & AliVEvent::kINT7;
    if (!isINT7selected) return;

    AliAODEvent* fAOD = dynamic_cast<AliAODEvent*>(inEvMain->GetEvent());
  

     if(!fAOD){
                cout << "ERROR: no AOD event" << endl;
      return;
    }

   
    //fPIDResponse = inEvMain->GetPIDResponse(); 
    fPIDResponse = (AliPIDResponse *) inEvMain-> GetPIDResponse();

  //================================================================

   // ((TH1F*)((TList*)fOutput->FindObject("EventInput"))->FindObject("fHistMultiMain"))->Fill(fAOD->GetNumberOfTracks());
    ((TH1F*)((TList*)fOutput->FindObject("EventInput"))->FindObject("fHistMultiMain"))->Fill(1.5);;

    //====================================Pile Up================================================= 
    

if (!fEventCuts->AcceptEvent(fAOD)) {
      PostData(1, fOutput);
      return;
    }
  ((TH1F*)((TList*)fOutput->FindObject("EventInput"))->FindObject("fHistMultiMain"))->Fill(2.5);;
       
   const AliVVertex* primVertex = fEventCuts->GetPrimaryVertex(); // Best primary vertex available
    if (!primVertex) return;
    if ( ( TMath::Abs(primVertex->GetZ()) ) >= fPrimaryVertexCut) return ;
    Double_t lPVx = primVertex->GetX();
    Double_t lPVy = primVertex->GetY();
    Double_t lPVz = primVertex->GetZ();
    if (TMath::Abs(lPVx)<10e-5 && TMath::Abs(lPVy)<10e-5&& TMath::Abs(lPVz)<10e-5) return;
   Short_t binVertex = Short_t((lPVz+7.)/2.);

 ((TH1F*)((TList*)fOutput->FindObject("EventInput"))->FindObject("fHistMultiMain"))->Fill(3.5);;


 nV0 = (fAOD->GetNumberOfV0s());                  // see how many V0 there are in the event

((TH1D*)((TList*)fOutput->FindObject("EventInput"))->FindObject("fHistV0Multiplicity"))->Fill(fAOD->GetNumberOfV0s());

    Float_t lCent = 0; 
AliMultSelection *MultSelection = 0x0; 
MultSelection = (AliMultSelection *) fAOD->FindListObject("MultSelection");
if( !MultSelection) {
   AliWarning("AliMultSelection object not found!");
}else{
   lCent = MultSelection->GetMultiplicityPercentile("V0M");//
}

   if ((lCent < fCentMin)||(lCent >= fCentMax )) return;

    ((TH1F*)((TList*)fOutput->FindObject("EventInput"))->FindObject("fhEventAf"))->Fill(lCent);
    ((TH2F*)((TList*)fOutput->FindObject("EventInput"))->FindObject("fHistCentVtx"))->Fill(lCent,lPVz);
    ((TH1F*)((TList*)fOutput->FindObject("EventInput"))->FindObject("fhEventCentAfterPilp"))->Fill(lCent);
   

//here is the =TObjArray for Mc Closure test  
   TObjArray *mcAssocTracks = new TObjArray;
    mcAssocTracks->SetOwner(kTRUE);
   
    TObjArray *MCk0s = new TObjArray;
    MCk0s->SetOwner(kTRUE);
    TObjArray *MCLambda = new TObjArray;

    MCLambda->SetOwner(kTRUE);
    TObjArray *MCAntiLambda = new TObjArray;
    MCAntiLambda->SetOwner(kTRUE);

    TClonesArray *mcArray = new TClonesArray;
    mcArray->SetOwner(kTRUE);
    //======================================processing MC data==========================================

    Int_t iMC = 0;
    if(fAnalysisMC){
    AliAODMCHeader *aodMCheader = (AliAODMCHeader*)fAOD->FindListObject(AliAODMCHeader::StdBranchName());
    Float_t vzMC = aodMCheader->GetVtxZ();
    if (TMath::Abs(vzMC) >= fPrimaryVertexCut) return;

    //retreive MC particles from event
   fMCArray = (TClonesArray*)fAOD->FindListObject(AliAODMCParticle::StdBranchName());
//    fMCArray = (TClonesArray*)fAOD->GeTList()->FindObject(AliAODMCParticle::StdBranchName());
    if(!fMCArray){
      Printf("No MC particle branch found");
      return;
    }
    //--------------------------------------------------------------------------------

    Int_t nMCAllTracks =fMCArray->GetEntriesFast();
    // new tracks array - without injected signal
    TObjArray * mcTracks = new TObjArray;
    //selectedMCTracks->SetOwner(kTRUE);
 
    for (Int_t i = 0; i < nMCAllTracks; i++){
        AliAODMCParticle *mcTrack = (AliAODMCParticle*)fMCArray->At(i);
        if (!mcTrack) {
            Error("ReadEventAODMC", "Could not receive particle %d", i);
            continue;
        }

       mcTracks->Add(mcTrack);
    }

    Int_t nMCTracks = mcTracks->GetEntriesFast();
    for (iMC = 0; iMC < nMCTracks; iMC++){
      AliAODMCParticle *mcTrack = (AliAODMCParticle*)mcTracks->At(iMC);
      if (!mcTrack) {
        Error("ReadEventAODMC", "Could not receive particle %d", iMC);
        continue;
      }
         
      //track processing
      Double_t mcTrackEta = mcTrack->Eta();
      Double_t mcTrackPt = mcTrack->Pt();
      Double_t mcTrackPhi = mcTrack->Phi();
      Bool_t TrIsPrim = mcTrack->IsPhysicalPrimary();
      Bool_t TrEtaMax = TMath::Abs(mcTrackEta)<0.8;
      Bool_t TrPtMin = mcTrackPt>fTrackPtMin;
      Bool_t TrCharge = (mcTrack->Charge())!=0;
      if (TrIsPrim && TrEtaMax && TrPtMin && TrCharge){


         fHistMCtruthTrkPtEta->Fill(mcTrackPt,mcTrackEta);
         fHistMCtruthTrk[binVertex]->Fill(mcTrackPt,mcTrackEta,mcTrackPhi);
         mcAssocTracks->Add(new AliV0XiParticles(mcTrack->Eta(),mcTrack->Phi(),mcTrack->Pt(),3,mcTrack->GetLabel(),mcTrack->GetLabel()));

              } 
 
            }
  


//loop for V0 
 for (Int_t i = 0; i < nMCAllTracks; i++){
        AliAODMCParticle *mcTrack = (AliAODMCParticle*)fMCArray->At(i);
        if (!mcTrack) {
            Error("ReadEventAODMC", "Could not receive particle %d", i);
            continue;
        }


if ((mcTrack->GetStatus() == 21)
          ||(mcTrack->GetPdgCode() == 443 && mcTrack->GetMother() == -1))
     break;

      Double_t mcTrackEta = mcTrack->Eta();
      Double_t mcTrackPt = mcTrack->Pt();
      Double_t mcTrackPhi = mcTrack->Phi();
      if(!mcTrack->IsPhysicalPrimary()) continue ;

  //V0 processing
    if(TMath::Abs(mcTrackEta) < fV0Eta && mcTrackPt > fV0PtMin && mcTrackPt < fV0PtMax){    
        if(TMath::Abs(mcTrack->GetPdgCode()) == 310){

          fHistMCtruthK0sPt->Fill(mcTrackPt,mcTrackEta); 
          fHistMCtruthK0s[binVertex]->Fill(mcTrackPt,mcTrackEta,mcTrackPhi);   

               MCk0s->Add(new AliV0XiParticles(mcTrackEta,mcTrackPhi,mcTrackPt,0,mcTrack->GetLabel(),mcTrack->M()));

}
 else if(mcTrack->GetPdgCode() == 3122){

          fHistMCtruthLambdaPt->Fill(mcTrackPt,mcTrackEta);
          fHistMCtruthLambda[binVertex]->Fill(mcTrackPt,mcTrackEta,mcTrackPhi);

          MCLambda->Add(new AliV0XiParticles(mcTrackEta,mcTrackPhi,mcTrackPt,1,mcTrack->GetLabel(),mcTrack->M()));


}
else if(mcTrack->GetPdgCode() == -3122){

         fHistMCtruthAntiLambdaPt->Fill(mcTrackPt,mcTrackEta);
         fHistMCtruthAntiLambda[binVertex]->Fill(mcTrackPt,mcTrackEta,mcTrackPhi); 
          MCAntiLambda->Add(new AliV0XiParticles(mcTrackEta,mcTrackPhi,mcTrackPt,2,mcTrack->GetLabel(),mcTrack->M()));


             }//end for Antilambda

           }//end for cut v0
        
}//end the first loop 


//======================================Analysis MC truth
// correlation  truth for  k0s 

for(Int_t i = 0; i < MCk0s->GetEntries(); i++){
        AliV0XiParticles* k0strig = (AliV0XiParticles*)MCk0s->At(i);
        if(!k0strig) continue;
        if (TMath::Abs(k0strig->Eta())>=0.7) continue;

     Double_t triggerk0sMCPt  = k0strig->Pt();
     Double_t triggerk0sMCPhi = k0strig->Phi();
     Double_t triggerk0sMCEta = k0strig->Eta(); 

      // cout<<triggerk0sMCPt<<endl;  
   
      
       Double_t MCk0sAll[4] = {triggerk0sMCPt, lCent, lPVz,1.};
        ((THnSparseF*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistGenTrigSibAllK0s"))->Fill(MCk0sAll);


    for(Int_t iTrk1 = 0; iTrk1 < mcAssocTracks->GetEntries(); iTrk1++){           
      AliV0XiParticles* assocMC = (AliV0XiParticles*) mcAssocTracks->At(iTrk1);
      if(!assocMC) continue;  


     Double_t assocpt = assocMC->Pt(); 
        //  cout<<assocpt<<endl;     

     
        if(assocpt>=triggerk0sMCPt) continue;

              
        Double_t  dPhiMC = triggerk0sMCPhi-assocMC->Phi();
	Double_t dEtaMC = triggerk0sMCEta - assocMC->Eta();
        if( dPhiMC > 1.5*TMath::Pi() ) dPhiMC -= 2.0*TMath::Pi();
        if( dPhiMC < -0.5*TMath::Pi() ) dPhiMC += 2.0*TMath::Pi();


    //here remove auto correlation 

                Int_t negID = k0strig->GetIDNeg();
                Int_t posID = k0strig->GetIDPos();
               Int_t atrID = assocMC->GetIDCh();

               if ((TMath::Abs(negID))==(TMath::Abs(atrID))) continue;
                if ((TMath::Abs(posID))==(TMath::Abs(atrID))) continue; 
      
              // cout<<negID<<endl;


 Double_t spMCSigK0s[7] = {dPhiMC, dEtaMC, triggerk0sMCPt, assocpt, lCent, lPVz, 1.};
((THnSparseF*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistGendPhidEtaSibK0s"))->Fill(spMCSigK0s);

}//end of track loop 
}//end loop for k0s correlation



// correlation  truth for  Lambda 

       
for (Int_t j=0; j <MCLambda->GetEntriesFast(); j++){
        AliV0XiParticles* Lambdatrig = (AliV0XiParticles*)MCLambda->At(j);
        if(!Lambdatrig) continue;
        if (TMath::Abs(Lambdatrig->Eta())>=0.7) continue;



     Double_t triggerLambdaMCPt  = Lambdatrig->Pt();
      Double_t triggerLambdaMCPhi = Lambdatrig->Phi();
      Double_t triggerLambdaMCEta = Lambdatrig->Eta();

          Double_t MCLambda[4] = {triggerLambdaMCPt, lCent, lPVz,1.};
         ((THnSparseF*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistGenTrigSibAllLambda"))->Fill(MCLambda);


    for(Int_t iTrk1 = 0; iTrk1 < mcAssocTracks->GetEntries(); iTrk1++){           
      AliV0XiParticles* assocMC = (AliV0XiParticles*) mcAssocTracks->At(iTrk1);

     Double_t assocpt = assocMC->Pt(); 
     //     cout<<assocpt<<endl;     
      if(!assocMC) continue;                
         
      if(assocpt>=triggerLambdaMCPt) continue;
	
        Double_t dPhiMC = triggerLambdaMCPhi-assocMC->Phi();
	Double_t dEtaMC = triggerLambdaMCEta - assocMC->Eta();

        if( dPhiMC > 1.5*TMath::Pi() ) dPhiMC -= 2.0*TMath::Pi();
        if( dPhiMC < -0.5*TMath::Pi() ) dPhiMC += 2.0*TMath::Pi();

//here remove auto correlation 

         Int_t negID = Lambdatrig->GetIDNeg();
         Int_t posID = Lambdatrig->GetIDPos();
         Int_t atrID = assocMC->GetIDCh();
               
               if ((TMath::Abs(negID))==(TMath::Abs(atrID))) continue;
                if ((TMath::Abs(posID))==(TMath::Abs(atrID))) continue;       

 Double_t spMCSigLambda[7] = {dPhiMC, dEtaMC, triggerLambdaMCPt, assocpt, lCent, lPVz, 1.};

((THnSparseF*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistGendPhidEtaSibLambda"))->Fill(spMCSigLambda);

  }//End Track loop 
   }  //end of Lambda   loop

 
// correlation  truth for AntiLambda 

     for (Int_t ii=0; ii<MCAntiLambda->GetEntriesFast(); ii++){
        AliV0XiParticles* AntiLambdatrig = (AliV0XiParticles*)MCAntiLambda->At(ii);
        if (TMath::Abs(AntiLambdatrig->Eta())>=0.7) continue;


      Double_t triggerAntiLambdaMCPt  = AntiLambdatrig->Pt();
      Double_t triggerAntiLambdaMCPhi = AntiLambdatrig->Phi();
      Double_t triggerAntiLambdaMCEta = AntiLambdatrig->Eta();


          Double_t MCAntiLambda[4] = {triggerAntiLambdaMCPt, lCent, lPVz,1.};
         ((THnSparseF*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistGenTrigSibAllAntiLambda"))->Fill(MCAntiLambda);

     for(Int_t iTrk1 = 0; iTrk1 < mcAssocTracks->GetEntries(); iTrk1++){           
      AliV0XiParticles* assocMC = (AliV0XiParticles*) mcAssocTracks->At(iTrk1);

     Double_t assocpt = assocMC->Pt(); 
     //     cout<<assocpt<<endl; 
      if(!assocMC) continue;
	if(assocMC->Pt()>triggerAntiLambdaMCPt) continue;    



        Double_t dPhiMC = triggerAntiLambdaMCPhi-assocMC->Phi();
	Double_t dEtaMC = triggerAntiLambdaMCEta - assocMC->Eta();
        
        if( dPhiMC > 1.5*TMath::Pi() ) dPhiMC -= 2.0*TMath::Pi();
        if( dPhiMC < -0.5*TMath::Pi() ) dPhiMC += 2.0*TMath::Pi();
              //here remove auto correlation 

                Int_t negID = AntiLambdatrig->GetIDNeg();
                Int_t posID = AntiLambdatrig->GetIDPos();
               Int_t atrID = assocMC->GetIDCh();
               
               if ((TMath::Abs(negID))==(TMath::Abs(atrID))) continue;
                if ((TMath::Abs(posID))==(TMath::Abs(atrID))) continue;       



 Double_t spMCSigAntiLambda[7] = {dPhiMC, dEtaMC, triggerAntiLambdaMCPt, assocpt, lCent, lPVz, 1.};

((THnSparseF*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistGendPhidEtaSibAntiLambda"))->Fill(spMCSigAntiLambda);


  }  //end of track loop condition           
 }  //end of AntiLambda   loop

   
//The end of  Corelations  MC Gen leve

     // access the reconstructed data
    Int_t nTracks = fAOD->GetNumberOfTracks();
    TObjArray * selectedMCTracks = new TObjArray;

    for (Int_t i = 0; i < nTracks; i++)
    {
      AliAODTrack* tr = dynamic_cast<AliAODTrack*>(fAOD->GetTrack(i));
      if(!tr) AliFatal("Not a standard AOD");
      selectedMCTracks->Add(tr);
    }

    Int_t nRecTracks = selectedMCTracks->GetEntriesFast();
    for(Int_t i = 0; i < nRecTracks; i++){
      AliAODTrack* tr = (AliAODTrack*)selectedMCTracks->At(i);
      if ( tr->Pt() < fTrackPtMin ) continue;
      if((tr->Pt())>fTrackPtMax) continue;

      if(tr->Charge() == 0.) continue;
      if (!(IsGoodPrimaryTrack(tr))) continue;
      Double_t TrackLabel = tr->GetLabel();
      if(tr->GetLabel() == -1){ 
        continue;
      }
          
      fHistRCTrkPtEta->Fill(tr->Pt(), tr->Eta());
      fHistRCTrk[binVertex]->Fill(tr->Pt(), tr->Eta(), tr->Phi()); 
      if(! (static_cast<AliAODMCParticle*>(fMCArray->At(TMath::Abs(tr->GetLabel()))))->IsPhysicalPrimary()){
        fHistRCSecTrkPtEta->Fill(tr->Pt(), tr->Eta());
        fHistRCSecTrk[binVertex]->Fill(tr->Pt(), tr->Eta(), tr->Phi());   
      }
      else{
        fHistRCPriTrkPtEta->Fill(tr->Pt(), tr->Eta());
        fHistRCPriTrk[binVertex]->Fill(tr->Pt(), tr->Eta(), tr->Phi());
      }
    }

   
     
}// End for MC condition

    //=====================================Track Selection==========================================
    TObjArray *selectedTracks = new TObjArray;
    selectedTracks->SetOwner(kTRUE);

    Int_t nTracks = fAOD->GetNumberOfTracks();
    for(Int_t i=0; i<nTracks; i++){
      AliAODTrack *tr = (AliAODTrack*)fAOD->GetTrack(i);         
      if((tr->Pt())<fTrackPtMin) continue;
      if((tr->Pt())>fTrackPtMax) continue;
      if(tr->Charge() == 0.) continue;
      if(!(IsGoodPrimaryTrack(tr))) continue;

   //Bunch rejection trk by trk
//if(!fBunchPileup){ if(!(tr->HasPointOnITSLayer(0) || tr->HasPointOnITSLayer(1)  || tr->GetTOFBunchCrossing()==0 )) continue;}
                                                                                                                 
   if(fRejectTrackPileUp&&(!(tr->HasPointOnITSLayer(0) || tr->HasPointOnITSLayer(1)|| tr->GetTOFBunchCrossing()==0 ))) continue;

      Double_t tPhi = tr->Phi();
      Double_t tPt = tr->Pt();
      Double_t tEta = tr->Eta();
      Double_t spTrack[4] = {tPt, tEta, lCent, lPVz};
      if(fEffCorr){
          Double_t weight = fHistEffEtaPtTrack->Interpolate(tr->Eta(), tr->Pt());
          if(weight == 0){
            continue;
          }
          ((THnSparseF*)((TList*)fOutput3->FindObject("Track"))->FindObject("fHistMassTrack"))->Fill(spTrack, 1/weight);
      }
      else{
        ((THnSparseF*)((TList*)fOutput3->FindObject("Track"))->FindObject("fHistMassTrack"))->Fill(spTrack);
      }
      ((TH1F*)((TList*)fOutput3->FindObject("Track"))->FindObject("tPhi"))->Fill(tPhi);
      ((TH1F*)((TList*)fOutput3->FindObject("Track"))->FindObject("tPt"))->Fill(tPt);
      ((TH1F*)((TList*)fOutput3->FindObject("Track"))->FindObject("tEta"))->Fill(tEta);
      ((TH3F*)((TList*)fOutput3->FindObject("Track"))->FindObject("fHistTrk"))->Fill(tPt,tEta,tPhi);
      if(fAnalysisMC) {
        if(tr->GetLabel() == -1) continue;
        if(! (static_cast<AliAODMCParticle*>(fMCArray->At(TMath::Abs(tr->GetLabel()))))->IsPhysicalPrimary())continue;
        if(TMath::Abs(tr->GetLabel()) > iMC) continue;
      }
      selectedTracks->Add(tr);
    }

    Int_t nSelectedTracks = selectedTracks->GetEntries();
    
    ((TH1F*)((TList*)fOutput3->FindObject("Track"))->FindObject("nTracksBf"))->Fill(nTracks);
    ((TH1F*)((TList*)fOutput3->FindObject("Track"))->FindObject("nTracksAf"))->Fill(nSelectedTracks);
    //==============================================================================================
     std::map<int, int> labels;
     for (int i = 0; i < nTracks; i++) {
     const AliAODTrack *aodtrack = dynamic_cast<const AliAODTrack *>(fAOD->GetTrack(i));

       if (!aodtrack->TestFilterBit(128)) {

       // Skip TPC-only tracks
       if (aodtrack->GetID() < 0) continue;
       labels[aodtrack->GetID()] = i;
       }
     }

    //==========================================V0 Selection========================================
    const AliVVertex* primaryBestAODVtx = fEventCuts->GetPrimaryVertex();
    primaryBestAODVtx->GetXYZ(fBestPrimaryVtxPos);
     
    TObjArray * selectedK0s = new TObjArray();
    selectedK0s->SetOwner(kTRUE);
    TObjArray * selectedLambda = new TObjArray();
    selectedLambda->SetOwner(kTRUE);
    TObjArray * selectedAntiLambda = new TObjArray();
    selectedAntiLambda->SetOwner(kTRUE);
    TObjArray * trigParticles = new TObjArray();
    trigParticles->SetOwner(kTRUE);

    //loop for V0s
    Int_t nV0sTot = fAOD->GetNumberOfV0s();

    for (Int_t iV0 = 0; iV0 < nV0sTot; iV0++) {

      AliAODv0 *v0=fAOD->GetV0(iV0) ;
      if (!v0) continue;

     // if (!IsGoodV0(v0)) continue;
      
     Int_t oStatus = GetOStatus();
    if (!IsGoodV0(v0 ,oStatus)) continue;

      AliAODTrack *ntrack=(AliAODTrack *)v0->GetDaughter(1);
      AliAODTrack *ptrack=(AliAODTrack *)v0->GetDaughter(0);

      Int_t iPos, iNeg;
      if(ptrack->Charge() > 0){
      iPos = 0; iNeg = 1;
      }
      else{
        iPos = 1; iNeg = 0;
      }

      UInt_t pTrkID = ptrack->GetID() >= 0 ? ptrack->GetID() : -1-ptrack->GetID();
      UInt_t nTrkID = ntrack->GetID() >= 0 ? ntrack->GetID() : -1-ntrack->GetID();
      AliAODTrack *Ptrack = (ptrack->TestFilterBit(128))? dynamic_cast<AliAODTrack *>(fAOD->GetTrack(labels[-1-ptrack->GetID()])) : ptrack;
      AliAODTrack *Ntrack = (ntrack->TestFilterBit(128))? dynamic_cast<AliAODTrack *>(fAOD->GetTrack(labels[-1-ntrack->GetID()])) : ntrack;

      Double_t xyz[3]; 
      v0->GetSecondaryVtx(xyz);
      Double_t dx = xyz[0]-lPVx, dy = xyz[1]-lPVy, dz=xyz[2]-lPVz;                                                     
      Double_t lt = TMath::Sqrt(dx*dx + dy*dy + dz*dz); 
      Double_t lPt = v0->Pt();
      Double_t dlK = 0.4977*lt/lPt;//2017-8-28
      Double_t dlL = 1.1157*lt/lPt;//2017-8-28

     //-------------------------------------------------------------------------------------------------
      Float_t xyn = v0->DcaNegToPrimVertex();
      Float_t xyp = v0->DcaPosToPrimVertex();
      Double_t r2 = xyz[0]*xyz[0] + xyz[1]*xyz[1];
      Double_t R = v0->RadiusV0();
      Double_t lEta = v0->PseudoRapV0();                                                      
      Double_t lPhi = v0->Phi();
      Double_t lDCA = v0->DcaV0Daughters();

      Double_t lCPA = v0->CosPointingAngle(fBestPrimaryVtxPos);

      Double_t lDCA2PV = v0->DcaV0ToPrimVertex();
      Double_t lPtArmV0 = v0->PtArmV0();
      Double_t lAlphaV0 = v0->AlphaV0();
      if(v0->ChargeProng(iPos) <0) lAlphaV0 = -lAlphaV0;
      
      Double_t massK0s = v0->MassK0Short();
      Double_t massLambda = v0->MassLambda();
      Double_t massAntiLambda = v0->MassAntiLambda();

      Bool_t isNegPionForTPC = kFALSE;
      Bool_t isPosPionForTPC = kFALSE;
      Bool_t isNegProtonForTPC = kFALSE;
      Bool_t isPosProtonForTPC = kFALSE;

      //----------------------------------------PID--------------------------------------------------------
       
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(Ptrack,AliPID::kPion))<=fV0PIDSigma)
         isPosPionForTPC = kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(Ptrack,AliPID::kProton))<=fV0PIDSigma) 
         isPosProtonForTPC = kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(Ntrack,AliPID::kPion))<=fV0PIDSigma)
         isNegPionForTPC = kTRUE;
      if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(Ntrack,AliPID::kProton))<=fV0PIDSigma)
         isNegProtonForTPC = kTRUE;





// reject bunch-off pile-up
  if(fRejectV0PileUp&&(!((Ntrack->IsOn(AliAODTrack::kTPCrefit)&& Ntrack->IsOn(AliAODTrack::kITSrefit))||Ntrack->IsOn(AliAODTrack::kTOFout))))continue;
  if(fRejectV0PileUp&&(!((Ptrack->IsOn(AliAODTrack::kTPCrefit)&& Ptrack->IsOn(AliAODTrack::kITSrefit))||Ptrack->IsOn(AliAODTrack::kTOFout)))) continue;

      

 
      if(isPosPionForTPC && Ptrack->IsOn(AliESDtrack::kTPCin)){
         ((TH2F*)((TList*)fOutput4->FindObject("V0"))->FindObject("TPCdEdxOfPion"))->Fill(Ptrack->P()*Ptrack->Charge(),Ptrack->GetTPCsignal());
      }
      if(isPosProtonForTPC && Ptrack->IsOn(AliESDtrack::kTPCin)){
         ((TH2F*)((TList*)fOutput4->FindObject("V0"))->FindObject("TPCdEdxOfProton"))->Fill(Ptrack->P()*Ptrack->Charge(),Ptrack->GetTPCsignal());
      } 
      if(isNegPionForTPC && Ntrack->IsOn(AliESDtrack::kTPCin)){
         ((TH2F*)((TList*)fOutput4->FindObject("V0"))->FindObject("TPCdEdxOfPion"))->Fill(Ntrack->P()*Ntrack->Charge(),Ntrack->GetTPCsignal());
      }
      if(isNegProtonForTPC && Ntrack->IsOn(AliESDtrack::kTPCin)){
         ((TH2F*)((TList*)fOutput4->FindObject("V0"))->FindObject("TPCdEdxOfProton"))->Fill(Ntrack->P()*Ntrack->Charge(),Ntrack->GetTPCsignal());
      }
      //------------------------------------Candidate selection cut-------------------------------------------
     
        if(lPt < fV0PtMin || lPt > fV0PtMax) continue;

      if(TMath::Abs(lEta) > fV0Eta) continue;
      Bool_t ctK=kTRUE;  
      if(dlK >= fK0sLifeTimeMax || dlK <=fK0sLifeTimeMin) ctK=kFALSE; 
      Bool_t ctL=kTRUE;  
      if(dlL >= fLambdaLifeTimeMax || dlL <= fLambdaLifeTimeMin) ctL=kFALSE;
     
      //---------------------------------------Check vabribles--------------------------------------------------
      ((TH1F*)((TList*)fOutput4->FindObject("V0"))->FindObject("DLK"))->Fill(dlK);
      ((TH1F*)((TList*)fOutput4->FindObject("V0"))->FindObject("DLL"))->Fill(dlL);
      ((TH1F*)((TList*)fOutput4->FindObject("V0"))->FindObject("AfV0Pt"))->Fill(lPt);
      ((TH1F*)((TList*)fOutput4->FindObject("V0"))->FindObject("AfV0Eta"))->Fill(lEta);
      ((TH1F*)((TList*)fOutput4->FindObject("V0"))->FindObject("Afxyp"))->Fill(xyp);
      ((TH1F*)((TList*)fOutput4->FindObject("V0"))->FindObject("Afxyn"))->Fill(xyn);   
      ((TH1F*)((TList*)fOutput4->FindObject("V0"))->FindObject("AfDCA"))->Fill(lDCA);
      ((TH1F*)((TList*)fOutput4->FindObject("V0"))->FindObject("AfCPA"))->Fill(lCPA);
      ((TH1F*)((TList*)fOutput4->FindObject("V0"))->FindObject("Afr2"))->Fill(r2);
      ((TH2F*)((TList*)fOutput4->FindObject("V0"))->FindObject("AfAP"))->Fill(lPtArmV0, lAlphaV0);
      
      //--------------------check whether it is K0s/ Lambda/ AntiLambda candidates------------------------------
      if(ctK &&lCPA > fk0sCPA && lPtArmV0 > TMath::Abs(fPtArmV0AlphaV0 *lAlphaV0) && xyn > fDCANegtoPrimVertexMink0s && xyp > fDCAPostoPrimVertexMink0s && isPosPionForTPC  && isNegPionForTPC && (massK0s > 0.40 )&& (massK0s < 0.58)){
        selectedK0s->Add(v0);
        Double_t spK0s[4] = {massK0s, lPt, lCent, lPVz};
        if(fEffCorr){
          Double_t weight = fHistEffEtaPtK0s->Interpolate(v0->Eta(), v0->Pt());
          if(weight == 0){
            continue;
          }

          ((THnSparseF*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistMassK0s"))->Fill(spK0s, 1/weight);

         ((TH2F*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistK0sMassvsPtCorr"))->Fill(lPt,massK0s,1/weight);  

       // ((TH2F*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistK0sMassvsPtNoCorr"))->Fill(lPt,massK0s);  
	 
	  
        }
        else{

          ((THnSparseF*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistMassK0s"))->Fill(spK0s);
	  
        }
        ((TH3F*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistK0s"))->Fill(lPt, lEta, lPhi);
        ((TH1F*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistK0sPhi"))->Fill(lPhi);
        ((TH1F*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistK0sEta"))->Fill(lEta);
        if(fAnalysisMC && IsMCV0Primary(v0, 0)){
          fHistRCK0sPt->Fill(lPt,lEta);
          fHistRCK0s[binVertex]->Fill(lPt,lEta,lPhi); 

        }
        trigParticles->Add(new AliV0XiParticles(lEta,lPhi,lPt,v0->MassK0Short(), 0));

      }
    
      // check whether it is Lambda candidates
      if(ctL && lCPA > fLambdaCPA && xyn > fDCANegtoPrimVertexMinLamb && xyp > fDCAPostoPrimVertexMinLamb && isPosProtonForTPC && isNegPionForTPC && (massLambda > 1.07) && (massLambda < 1.15)){
        selectedLambda->Add(v0);
        Double_t spLambda[4] = {massLambda, lPt, lCent, lPVz};
        if(fEffCorr){
          Double_t weight = fHistEffEtaPtLambda->Interpolate(v0->Eta(), v0->Pt());
          if(weight == 0){
            continue;
         }


 ((THnSparseF*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistMassLambda"))->Fill(spLambda, 1/weight);

        ((TH2F*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistLamMassvsPtCorr"))->Fill(lPt,massLambda,1/weight);  
       // ((TH2F*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistLamMassvsPtNoCorr"))->Fill(lPt,massLambda);  



        }
        else{


          ((THnSparseF*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistMassLambda"))->Fill(spLambda);

        }
        ((TH3F*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistLambda"))->Fill(lPt, lEta, lPhi);
        ((TH1F*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistLambdaPhi"))->Fill(lPhi);
        ((TH1F*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistLambdaEta"))->Fill(lEta);
        if(fAnalysisMC){
          if(IsMCV0Primary(v0, 1)){


            fHistRCLambdaPt->Fill(lPt,lEta);
            fHistRCLambda[binVertex]->Fill(lPt,lEta,lPhi);
            fHistRCPrimLambdaDCAPt->Fill(lDCA2PV,lPt);


          }
          else if(IsMCV0FromXi(v0, 1)){
            fHistRCSecLambdaDCAPt->Fill(lDCA2PV,lPt);
          }
        }
        trigParticles->Add(new AliV0XiParticles(lEta,lPhi,lPt,massLambda, 1)); 
      }

    // check whether it is AntiLambda candidates   
      if(ctL && lCPA > fLambdaCPA && xyn > fDCANegtoPrimVertexMinALamb && xyp > fDCAPostoPrimVertexMinALamb && isPosPionForTPC && isNegProtonForTPC && (massAntiLambda > 1.07) && (massAntiLambda < 1.15)){
        selectedAntiLambda->Add(v0);
        Double_t spAntiLambda[4] = {massAntiLambda, lPt, lCent, lPVz};
        if(fEffCorr){
          Double_t weight = fHistEffEtaPtAntiLambda->Interpolate(v0->Eta(), v0->Pt());
          if(weight == 0){
            continue;
           }


          ((THnSparseF*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistMassAntiLambda"))->Fill(spAntiLambda, 1/weight);

        ((TH2F*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistAntiLamMassvsPtCorr"))->Fill(lPt, massAntiLambda,1/weight);  
        //((TH2F*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistAntiLamMassvsPtNoCorr"))->Fill(lPt,massAntiLambda);  
	  


	  
        }
        else{

          ((THnSparseF*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistMassAntiLambda"))->Fill(spAntiLambda);
	  
        }
        ((TH3F*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistAntiLambda"))->Fill(lPt, lEta, lPhi);
        ((TH1F*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistAntiLambdaPhi"))->Fill(lPhi);
        ((TH1F*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistAntiLambdaEta"))->Fill(lEta);
        if(fAnalysisMC){
          if(IsMCV0Primary(v0, 2)){

            fHistRCAntiLambdaPt->Fill(lPt,lEta);
            fHistRCAntiLambda[binVertex]->Fill(lPt,lEta,lPhi);
            fHistRCPrimAntiLambdaDCAPt->Fill(lDCA2PV,lPt);

          }
          else if(IsMCV0FromXi(v0, 2)){
            fHistRCSecAntiLambdaDCAPt->Fill(lDCA2PV,lPt);
          }
        }
        trigParticles->Add(new AliV0XiParticles(lEta,lPhi,lPt, massAntiLambda, 2)); 
      }

       } //here is the end of v0 loop


    //=====================================  Analysis ============================================
     
   
      for(Int_t i = 0; i < selectedK0s->GetEntries(); i++){
        AliAODv0 * v0 = (AliAODv0*)selectedK0s->At(i);
        if(!v0) continue;
        AliAODTrack *pTrkXi = dynamic_cast<AliAODTrack*>( v0->GetDaughter(0) );
        AliAODTrack *nTrkXi = dynamic_cast<AliAODTrack*>( v0->GetDaughter(1) );
      
        Int_t pTrkID = pTrkXi->GetID()  >= 0 ? pTrkXi->GetID() : -1-pTrkXi->GetID();
        Int_t nTrkID = nTrkXi->GetID() >= 0 ? nTrkXi->GetID() : -1-nTrkXi->GetID();


        Double_t massK0s = v0->MassK0Short();
        Double_t lPt=TMath::Sqrt(v0->Pt2V0());
        // Double_t lPt=v0->Pt();

     if(lPt<fV0PtMin||lPt>fV0PtMax) continue;
         // cout<<lPt<<endl;

if(fEffCorr){
          Double_t weight = fHistEffEtaPtK0s->Interpolate(v0->Eta(), v0->Pt());
          if(weight == 0){
            continue;
          }
    
        if (massK0s > fMassLowK0s[1] && massK0s < fMassHighK0s[1]){
          Double_t spTrigSigK0s[4] = {lPt, lCent, lPVz, 1.};

          ((THnSparseF*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistTrigSibAllK0s"))->Fill(spTrigSigK0s, 1/weight);
        }
        else if(massK0s > fMassLowK0s[0] && massK0s < fMassHighK0s[0]){
          Double_t spTrigSigLeftK0s[4] = {lPt, lCent, lPVz, 2.};

  ((THnSparseF*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistTrigSibAllK0s"))->Fill(spTrigSigLeftK0s, 1/weight);
        }
        else if(massK0s > fMassLowK0s[2] && massK0s < fMassHighK0s[2]){
          Double_t spTrigSigRightK0s[4] = {lPt, lCent, lPVz, 3.};

          ((THnSparseF*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistTrigSibAllK0s"))->Fill(spTrigSigRightK0s, 1/weight);
        }
      }
else{if (massK0s > fMassLowK0s[1] && massK0s < fMassHighK0s[1]){
          Double_t spTrigSigK0s[4] = {lPt, lCent, lPVz, 1.};

          ((THnSparseF*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistTrigSibAllK0s"))->Fill(spTrigSigK0s);
        }
        else if(massK0s > fMassLowK0s[0] && massK0s < fMassHighK0s[0]){
          Double_t spTrigSigLeftK0s[4] = {lPt, lCent, lPVz, 2.};

  ((THnSparseF*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistTrigSibAllK0s"))->Fill(spTrigSigLeftK0s);
        }
        else if(massK0s > fMassLowK0s[2] && massK0s < fMassHighK0s[2]){
          Double_t spTrigSigRightK0s[4] = {lPt, lCent, lPVz, 3.};

          ((THnSparseF*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistTrigSibAllK0s"))->Fill(spTrigSigRightK0s);
        }
      }
     
	 for(Int_t iTrk = 0; iTrk < selectedTracks->GetEntries(); iTrk++){           
      AliAODTrack* atr = (AliAODTrack*) selectedTracks->At(iTrk);                
      if(!atr) continue;                
      Int_t trID = atr->GetID() >= 0 ? atr->GetID() : -1-atr->GetID();
	   if(atr->Pt() >= lPt) continue;
      

        //Correlation part
        Double_t dEta = atr->Eta() - v0->Eta();
        Double_t dPhi = atr->Phi() - v0->Phi();
        if( dPhi > 1.5*TMath::Pi() ) dPhi -= 2.0*TMath::Pi();
        if( dPhi < -0.5*TMath::Pi() ) dPhi += 2.0*TMath::Pi();

        //remove auto correlations
        if( pTrkID == trID || nTrkID == trID ) continue;
        // cout<<"Hello! Removing atuocorrelations!"<<endl;

	if(lPt<fV0PtMin||lPt>fV0PtMax) continue;

              
        //Filling correlation histograms and histograms for triggers counting
        if(fEffCorr){
          Double_t weightK0s = fHistEffEtaPtK0s->Interpolate(v0->Eta(), v0->Pt());
          Double_t weightTrack = fHistEffEtaPtTrack->Interpolate(atr->Eta(), atr->Pt());
          Double_t weight = weightK0s * weightTrack;
          if(weight == 0){
            continue;
          }
         

              if (massK0s > fMassLowK0s[1] && massK0s < fMassHighK0s[1]){
              Double_t spSigK0s[7] = {dPhi, dEta, lPt, atr->Pt(), lCent, lPVz, 1.};
           	      
           ((THnSparseF*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistdPhidEtaSibK0s"))->Fill(spSigK0s, 1/weight);
          }
	      
          if (massK0s > fMassLowK0s[0] && massK0s < fMassHighK0s[0]){
              Double_t spBkgLeftK0s[7] = {dPhi, dEta, lPt, atr->Pt(), lCent, lPVz, 2.};

              ((THnSparseF*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistdPhidEtaSibK0s"))->Fill(spBkgLeftK0s, 1/weight);
          }
          if (massK0s > fMassLowK0s[2] && massK0s < fMassHighK0s[2]){
              Double_t spBkgRightK0s[7] = {dPhi, dEta, lPt, atr->Pt(), lCent, lPVz, 3.};

              ((THnSparseF*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistdPhidEtaSibK0s"))->Fill(spBkgRightK0s, 1/weight);
          }
        }  
        else{
          if (massK0s > fMassLowK0s[1] && massK0s < fMassHighK0s[1]){
              Double_t spSigK0s[7] = {dPhi, dEta, lPt, atr->Pt(), lCent, lPVz, 1.};

              ((THnSparseF*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistdPhidEtaSibK0s"))->Fill(spSigK0s);
          }
          if (massK0s > fMassLowK0s[0] && massK0s < fMassHighK0s[0]){
              Double_t spBkgLeftK0s[7] = {dPhi, dEta, lPt, atr->Pt(), lCent, lPVz, 2.};

              ((THnSparseF*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistdPhidEtaSibK0s"))->Fill(spBkgLeftK0s);
          }
          if (massK0s > fMassLowK0s[2] && massK0s < fMassHighK0s[2]){
              Double_t spBkgRightK0s[7] = {dPhi, dEta, lPt, atr->Pt(), lCent, lPVz, 3.};


              ((THnSparseF*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistdPhidEtaSibK0s"))->Fill(spBkgRightK0s);
          }
        }
	 }//end trak
      }//end k0s
      for(Int_t i = 0; i < selectedLambda->GetEntries(); i++){
        AliAODv0 * v0 = (AliAODv0*)selectedLambda->At(i);
        if(!v0) continue;
        AliAODTrack *pTrkXi = dynamic_cast<AliAODTrack*>( v0->GetDaughter(0) );
        AliAODTrack *nTrkXi = dynamic_cast<AliAODTrack*>( v0->GetDaughter(1) );
      
        Int_t pTrkID = pTrkXi->GetID()  >= 0 ? pTrkXi->GetID() : -1-pTrkXi->GetID();
        Int_t nTrkID = nTrkXi->GetID() >= 0 ? nTrkXi->GetID() : -1-nTrkXi->GetID();
      
        Double_t massLambda = v0->MassLambda();
        Double_t lPt=TMath::Sqrt(v0->Pt2V0());

     if(lPt<fV0PtMin||lPt>fV0PtMax) continue;


if(fEffCorr){
          Double_t weight = fHistEffEtaPtLambda->Interpolate(v0->Eta(), v0->Pt());
          if(weight == 0){
            continue;
         }

        if (massLambda > fMassLowLambda[1] && massLambda < fMassHighLambda[1]){
          Double_t spTrigSigLambda[4] = {lPt, lCent, lPVz, 1.};

          ((THnSparseF*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistTrigSibAllLambda"))->Fill(spTrigSigLambda, 1/weight);
        }
        else if(massLambda > fMassLowLambda[0] && massLambda < fMassHighLambda[0]){
          Double_t spTrigSigLeftLambda[4] = {lPt, lCent, lPVz, 2.};

          ((THnSparseF*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistTrigSibAllLambda"))->Fill(spTrigSigLeftLambda, 1/weight);
        }
        else if(massLambda > fMassLowLambda[2] && massLambda < fMassHighLambda[2]){
          Double_t spTrigSigRightLambda[4] = {lPt, lCent, lPVz, 3.};

        ((THnSparseF*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistTrigSibAllLambda"))->Fill(spTrigSigRightLambda, 1/weight);
        }
    }
    

else{
if (massLambda > fMassLowLambda[1] && massLambda < fMassHighLambda[1]){
          Double_t spTrigSigLambda[4] = {lPt, lCent, lPVz, 1.};

          ((THnSparseF*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistTrigSibAllLambda"))->Fill(spTrigSigLambda);
        }
        else if(massLambda > fMassLowLambda[0] && massLambda < fMassHighLambda[0]){
          Double_t spTrigSigLeftLambda[4] = {lPt, lCent, lPVz, 2.};

          ((THnSparseF*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistTrigSibAllLambda"))->Fill(spTrigSigLeftLambda);
        }
        else if(massLambda > fMassLowLambda[2] && massLambda < fMassHighLambda[2]){
          Double_t spTrigSigRightLambda[4] = {lPt, lCent, lPVz, 3.};

        ((THnSparseF*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistTrigSibAllLambda"))->Fill(spTrigSigRightLambda);
        }
    }
     for(Int_t iTrk = 0; iTrk < selectedTracks->GetEntries(); iTrk++){           
      AliAODTrack* atr = (AliAODTrack*) selectedTracks->At(iTrk);                
      if(!atr) continue;                
      
      Int_t trID = atr->GetID() >= 0 ? atr->GetID() : -1-atr->GetID();



        if(atr->Pt() >= lPt) continue;

        //Correlation part
        Double_t dEta = atr->Eta() - v0->Eta();
        Double_t dPhi = atr->Phi() - v0->Phi();
        if( dPhi > 1.5*TMath::Pi() ) dPhi -= 2.0*TMath::Pi();
        if( dPhi < -0.5*TMath::Pi() ) dPhi += 2.0*TMath::Pi();

	if(lPt<fV0PtMin||lPt>fV0PtMax) continue;


        //remove auto correlations
        if( pTrkID == trID || nTrkID == trID ) continue;
        // cout<<"Hello! Removing atuocorrelations!"<<endl;

        if (massLambda > fMassLowLambda[1] && massLambda < fMassHighLambda[1]){
          ((TH3F*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistLambdaDphiDCAPtSig"))->Fill(dPhi, lPt, v0->DcaV0ToPrimVertex());
        }
        else if(massLambda > fMassLowLambda[0] && massLambda < fMassHighLambda[0]){
          ((TH3F*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistLambdaDphiDCAPtBkgL"))->Fill(dPhi, lPt, v0->DcaV0ToPrimVertex());
        }
        else if(massLambda > fMassLowLambda[2] && massLambda < fMassHighLambda[2]){
          ((TH3F*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistLambdaDphiDCAPtBkgR"))->Fill(dPhi, lPt, v0->DcaV0ToPrimVertex());
        }
        //Filling correlation histograms and histograms for triggers counting
        
        if(fEffCorr){
          Double_t weightLambda = fHistEffEtaPtLambda->Interpolate(v0->Eta(), v0->Pt());
          Double_t weightTrack = fHistEffEtaPtTrack->Interpolate(atr->Eta(), atr->Pt());
          Double_t weight = weightLambda * weightTrack;
          if(weight == 0){
            continue;
          }

	  
 
             if(massLambda > fMassLowLambda[1] && massLambda < fMassHighLambda[1]){
            Double_t spSigLambda[7] = {dPhi, dEta, lPt, atr->Pt(), lCent, lPVz, 1.};

          
            ((THnSparseF*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistdPhidEtaSibLambda"))->Fill(spSigLambda, 1/weight);
          }
          if(massLambda > fMassLowLambda[0] && massLambda < fMassHighLambda[0]){
            Double_t spBkgLeftLambda[7] = {dPhi, dEta, lPt, atr->Pt(), lCent, lPVz, 2.};

            ((THnSparseF*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistdPhidEtaSibLambda"))->Fill(spBkgLeftLambda, 1/weight);
          }
          if(massLambda > fMassLowLambda[2] && massLambda < fMassHighLambda[2]){
	   Double_t spBkgRightLambda[7] = {dPhi, dEta, lPt, atr->Pt(), lCent, lPVz, 3.};

   
            ((THnSparseF*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistdPhidEtaSibLambda"))->Fill(spBkgRightLambda, 1/weight);
          }
        }
        else{
          if(massLambda > fMassLowLambda[1] && massLambda < fMassHighLambda[1]){
            Double_t spSigLambda[7] = {dPhi, dEta, lPt, atr->Pt(), lCent, lPVz, 1.};

            ((THnSparseF*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistdPhidEtaSibLambda"))->Fill(spSigLambda);
          }
          if(massLambda > fMassLowLambda[0] && massLambda < fMassHighLambda[0]){
            Double_t spBkgLeftLambda[7] = {dPhi, dEta, lPt, atr->Pt(), lCent, lPVz, 2.};


            ((THnSparseF*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistdPhidEtaSibLambda"))->Fill(spBkgLeftLambda);
          }
          if(massLambda > fMassLowLambda[2] && massLambda < fMassHighLambda[2]){
            Double_t spBkgRightLambda[7] = {dPhi, dEta, lPt, atr->Pt(), lCent, lPVz, 3.};

            ((THnSparseF*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistdPhidEtaSibLambda"))->Fill(spBkgRightLambda);
          }
        }
 }//end track
      }//end Lambda 
        
      for(Int_t i = 0; i < selectedAntiLambda->GetEntries(); i++){
        AliAODv0 * v0 = (AliAODv0*)selectedAntiLambda->At(i);
        if(!v0) continue;
        AliAODTrack *pTrkXi = dynamic_cast<AliAODTrack*>( v0->GetDaughter(0) );
        AliAODTrack *nTrkXi = dynamic_cast<AliAODTrack*>( v0->GetDaughter(1) );
      
        Int_t pTrkID = pTrkXi->GetID()  >= 0 ? pTrkXi->GetID() : -1-pTrkXi->GetID();
        Int_t nTrkID = nTrkXi->GetID() >= 0 ? nTrkXi->GetID() : -1-nTrkXi->GetID();
          
        Double_t massAntiLambda = v0->MassAntiLambda();
        Double_t lPt=TMath::Sqrt(v0->Pt2V0());

	if(lPt<fV0PtMin||lPt>fV0PtMax) continue;

 if(fEffCorr){
          Double_t weight = fHistEffEtaPtAntiLambda->Interpolate(v0->Eta(), v0->Pt());
          if(weight == 0){
            continue;
           }

 if(massAntiLambda > fMassLowAntiLambda[1] && massAntiLambda < fMassHighAntiLambda[1]){
          Double_t spTrigSigAntiLambda[4] = {lPt, lCent, lPVz, 1.};
          ((THnSparseF*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistTrigSibAllAntiLambda"))->Fill(spTrigSigAntiLambda, 1/weight);
       }
        else if(massAntiLambda > fMassLowAntiLambda[0] && massAntiLambda < fMassHighAntiLambda[0]){
          Double_t spTrigSigLeftAntiLambda[4] = {lPt, lCent, lPVz, 2.};

          ((THnSparseF*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistTrigSibAllAntiLambda"))->Fill(spTrigSigLeftAntiLambda, 1/weight);
        }
        else if(massAntiLambda > fMassLowAntiLambda[2] && massAntiLambda < fMassHighAntiLambda[2]){
          Double_t spTrigSigRightAntiLambda[4] = {lPt, lCent, lPVz, 3.};

         ((THnSparseF*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistTrigSibAllAntiLambda"))->Fill(spTrigSigRightAntiLambda, 1/weight);
        }

       }
else{
if(massAntiLambda > fMassLowAntiLambda[1] && massAntiLambda < fMassHighAntiLambda[1]){
          Double_t spTrigSigAntiLambda[4] = {lPt, lCent, lPVz, 1.};
          ((THnSparseF*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistTrigSibAllAntiLambda"))->Fill(spTrigSigAntiLambda);
       }
        else if(massAntiLambda > fMassLowAntiLambda[0] && massAntiLambda < fMassHighAntiLambda[0]){
          Double_t spTrigSigLeftAntiLambda[4] = {lPt, lCent, lPVz, 2.};

          ((THnSparseF*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistTrigSibAllAntiLambda"))->Fill(spTrigSigLeftAntiLambda);
        }
        else if(massAntiLambda > fMassLowAntiLambda[2] && massAntiLambda < fMassHighAntiLambda[2]){
          Double_t spTrigSigRightAntiLambda[4] = {lPt, lCent, lPVz, 3.};

         ((THnSparseF*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistTrigSibAllAntiLambda"))->Fill(spTrigSigRightAntiLambda);
        }

       }

 for(Int_t iTrk = 0; iTrk < selectedTracks->GetEntries(); iTrk++){           
      AliAODTrack* atr = (AliAODTrack*) selectedTracks->At(iTrk);                
      if(!atr) continue;                
      
      Int_t trID = atr->GetID() >= 0 ? atr->GetID() : -1-atr->GetID();

      
       
        if(atr->Pt() >= lPt) continue;

        //Correlation part
        Double_t dEta = atr->Eta() - v0->Eta();
        Double_t dPhi = atr->Phi() - v0->Phi();
        if( dPhi > 1.5*TMath::Pi() ) dPhi -= 2.0*TMath::Pi();
        if( dPhi < -0.5*TMath::Pi() ) dPhi += 2.0*TMath::Pi();

        //remove auto correlations
        if( pTrkID == trID || nTrkID == trID ) continue;
        // cout<<"Hello! Removing atuocorrelations!"<<endl;
         
        //Filling correlation histograms and histograms for triggers counting
         if(massAntiLambda > fMassLowAntiLambda[1] && massAntiLambda < fMassHighAntiLambda[1]){
           ((TH3F*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistAntiLambdaDphiDCAPtSig"))->Fill(dPhi, lPt, v0->DcaV0ToPrimVertex());
         }
         else if(massAntiLambda > fMassLowAntiLambda[0] && massAntiLambda < fMassHighAntiLambda[0]){
           ((TH3F*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistAntiLambdaDphiDCAPtBkgL"))->Fill(dPhi, lPt, v0->DcaV0ToPrimVertex());
         }
         else if(massAntiLambda > fMassLowAntiLambda[2] && massAntiLambda < fMassHighAntiLambda[2]){
          ((TH3F*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistAntiLambdaDphiDCAPtBkgR"))->Fill(dPhi, lPt, v0->DcaV0ToPrimVertex()); 
         }

         if(fEffCorr){
          Double_t weightAntiLambda = fHistEffEtaPtAntiLambda->Interpolate(v0->Eta(), v0->Pt());
          Double_t weightTrack = fHistEffEtaPtTrack->Interpolate(atr->Eta(), atr->Pt());
          Double_t weight = weightAntiLambda * weightTrack;
          if(weight == 0){
            continue;
          }
	 
            if(massAntiLambda > fMassLowAntiLambda[1] && massAntiLambda < fMassHighAntiLambda[1]){
            Double_t spSigAntiLambda[7] = {dPhi, dEta, lPt, atr->Pt(), lCent, lPVz, 1.};

 
           ((THnSparseF*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistdPhidEtaSibAntiLambda"))->Fill(spSigAntiLambda, 1/weight);
          }
          if(massAntiLambda > fMassLowAntiLambda[0] && massAntiLambda < fMassHighAntiLambda[0]){
            Double_t spBkgLeftAntiLambda[7] = {dPhi, dEta, lPt, atr->Pt(), lCent, lPVz, 2.};


            ((THnSparseF*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistdPhidEtaSibAntiLambda"))->Fill(spBkgLeftAntiLambda, 1/weight);
          }
          if(massAntiLambda > fMassLowAntiLambda[2] && massAntiLambda < fMassHighAntiLambda[2]){
            Double_t spBkgRightAntiLambda[7] = {dPhi, dEta, lPt, atr->Pt(), lCent, lPVz, 3.};

            ((THnSparseF*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistdPhidEtaSibAntiLambda"))->Fill(spBkgRightAntiLambda, 1/weight);   
          }
        }
        else{
          if(massAntiLambda > fMassLowAntiLambda[1] && massAntiLambda < fMassHighAntiLambda[1]){
            Double_t spSigAntiLambda[7] = {dPhi, dEta, lPt, atr->Pt(), lCent, lPVz, 1.};

            ((THnSparseF*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistdPhidEtaSibAntiLambda"))->Fill(spSigAntiLambda);
          }
          if(massAntiLambda > fMassLowAntiLambda[0] && massAntiLambda < fMassHighAntiLambda[0]){
            Double_t spBkgLeftAntiLambda[7] = {dPhi, dEta, lPt, atr->Pt(), lCent, lPVz, 2.};

            ((THnSparseF*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistdPhidEtaSibAntiLambda"))->Fill(spBkgLeftAntiLambda);
          }
          if(massAntiLambda > fMassLowAntiLambda[2] && massAntiLambda < fMassHighAntiLambda[2]){
            Double_t spBkgRightAntiLambda[7] = {dPhi, dEta, lPt, atr->Pt(), lCent, lPVz, 3.};

            ((THnSparseF*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistdPhidEtaSibAntiLambda"))->Fill(spBkgRightAntiLambda);
          }
        }
      }//end of track for the current events
    }//end of Antilambda


    
      // -----------------------------------------------------Mixing part--------------------------------------------------------------------
    //This the right positon of this pool ( i transfer it from above ) 
       AliEventPool* pool = fPoolMgr->GetEventPool(lCent, lPVz);
      if (!pool)
      AliFatal(Form("No pool found for centrality = %f, zVtx = %f", lCent, lPVz));


     // if (pool->IsReady()|| pool->NTracksInPool() > 2000|| pool->GetCurrentNEvents() >= 5) {//masi change the down line 

    if (pool->IsReady() || pool->NTracksInPool() > fMixingTracks / 10 || pool->GetCurrentNEvents() >=5){

        Int_t nMix = pool->GetCurrentNEvents();
        for (Int_t jMix=0; jMix< nMix; jMix++){

          //loop through mixing events

          TObjArray* bgTracks = pool->GetEvent(jMix);

          for(Int_t i=0; i<trigParticles->GetEntriesFast(); i++){

            AliV0XiParticles* trig = (AliV0XiParticles*) trigParticles->At(i);




            //loop through V0 particles 

          for (Int_t j = 0; j < bgTracks->GetEntriesFast(); j++){

            // mixing tracks loop
            AliVParticle* atr = (AliVParticle*) bgTracks->At(j); 

            // be careful tracks may have bigger pt than trigger's.

            if ( ( (atr->Pt())>=(trig->Pt()))||( (atr->Pt())<fTrackPtMin ) ) continue;

            Double_t dEtaMix = atr->Eta() - trig->Eta();
            Double_t dPhiMix = atr->Phi() - trig->Phi();
            if ( dPhiMix > 1.5*TMath::Pi() ) dPhiMix -= 2.0*TMath::Pi();
            if ( dPhiMix < -0.5*TMath::Pi() ) dPhiMix += 2.0*TMath::Pi();
            
            if(trig->WhichCandidate() == 0){
              Double_t massk0s = trig->M();
              Double_t lk0sPt = trig->Pt();

	      

if(fAnalysisMC){

                  Double_t spSigMixK0s[7] = {dPhiMix, dEtaMix, lk0sPt, atr->Pt(),lCent, lPVz, 1.};
((THnSparseF*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistGendPhidEtaMixK0s"))->Fill(spSigMixK0s);
              
}

 
if(fEffCorr){
                Double_t weightK0s = fHistEffEtaPtK0s->Interpolate(trig->Eta(),trig->Pt());
                Double_t weightTrack = fHistEffEtaPtTrack->Interpolate(atr->Eta(), atr->Pt());
                Double_t weight = weightK0s * weightTrack;
                if(weight == 0){
                  continue;
                }

              
		
                  if(massk0s > fMassLowK0s[1] && massk0s < fMassHighK0s[1]){
                  Double_t spSigMixK0s[7] = {dPhiMix, dEtaMix, lk0sPt, atr->Pt(),lCent, lPVz, 1.};
                  ((THnSparseF*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistdPhidEtaMixK0s"))->Fill(spSigMixK0s, 1/weight);
                }

                if(massk0s > fMassLowK0s[0] && massk0s < fMassHighK0s[0]){
                  Double_t spBkgLeftMixK0s[7] = {dPhiMix, dEtaMix, lk0sPt, atr->Pt(), lCent, lPVz, 2.};

               ((THnSparseF*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistdPhidEtaMixK0s"))->Fill(spBkgLeftMixK0s, 1/weight);
                }
                if(massk0s > fMassLowK0s[2] && massk0s < fMassHighK0s[2]){
                  Double_t spBkgRightMixK0s[7] = {dPhiMix, dEtaMix, lk0sPt, atr->Pt(), lCent, lPVz, 3.};


               ((THnSparseF*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistdPhidEtaMixK0s"))->Fill(spBkgRightMixK0s, 1/weight);
                }
              }


              else{
               if(massk0s > fMassLowK0s[1] && massk0s < fMassHighK0s[1]){
                  Double_t spSigMixK0s[7] = {dPhiMix, dEtaMix, lk0sPt, atr->Pt(),lCent, lPVz, 1.};


                  ((THnSparseF*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistdPhidEtaMixK0s"))->Fill(spSigMixK0s);
                }

                if(massk0s > fMassLowK0s[0] && massk0s < fMassHighK0s[0]){
                  Double_t spBkgLeftMixK0s[7] = {dPhiMix, dEtaMix, lk0sPt, atr->Pt(), lCent, lPVz, 2.};




                  ((THnSparseF*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistdPhidEtaMixK0s"))->Fill(spBkgLeftMixK0s);
                }
                if(massk0s > fMassLowK0s[2] && massk0s < fMassHighK0s[2]){
                  Double_t spBkgRightMixK0s[7] = {dPhiMix, dEtaMix, lk0sPt, atr->Pt(), lCent, lPVz, 3.};


                  ((THnSparseF*)((TList*)fOutput5->FindObject("K0s"))->FindObject("fHistdPhidEtaMixK0s"))->Fill(spBkgRightMixK0s);

                }

              }
            } 




            if(trig->WhichCandidate() == 1){  
              Double_t masslambda = trig->M();
              Double_t llambdaPt = trig->Pt();


	     

if(fAnalysisMC){

                  Double_t spSigMixLambda[7] = {dPhiMix, dEtaMix, llambdaPt, atr->Pt(), lCent, lPVz, 1.};
((THnSparseF*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistGendPhidEtaMixLambda"))->Fill(spSigMixLambda);                               

}





              if(fEffCorr){
                Double_t weightLambda = fHistEffEtaPtLambda->Interpolate(trig->Eta(),trig->Pt());
                Double_t weightTrack = fHistEffEtaPtTrack->Interpolate(atr->Eta(), atr->Pt());
                Double_t weight = weightLambda * weightTrack;
                if(weight == 0){
                  continue;
                }


                

 if(masslambda > fMassLowLambda[1] && masslambda < fMassHighLambda[1]){
                  Double_t spSigMixLambda[7] = {dPhiMix, dEtaMix, llambdaPt, atr->Pt(), lCent, lPVz, 1.};

     ((THnSparseF*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistdPhidEtaMixLambda"))->Fill(spSigMixLambda, 1/weight);                               
                }
                if(masslambda > fMassLowLambda[0] && masslambda < fMassHighLambda[0]){
                  Double_t spBkgLeftMixLambda[7] = {dPhiMix, dEtaMix, llambdaPt, atr->Pt(), lCent, lPVz, 2.};


            ((THnSparseF*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistdPhidEtaMixLambda"))->Fill(spBkgLeftMixLambda, 1/weight);   
                }
                if(masslambda > fMassLowLambda[2] && masslambda < fMassHighLambda[2]){
                  Double_t spBkgRightMixLambda[7] = {dPhiMix, dEtaMix, llambdaPt, atr->Pt(), lCent, lPVz, 3.};


             ((THnSparseF*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistdPhidEtaMixLambda"))->Fill(spBkgRightMixLambda, 1/weight);            
                }
              }
              else{
                if(masslambda > fMassLowLambda[1] && masslambda < fMassHighLambda[1]){
                  Double_t spSigMixLambda[7] = {dPhiMix, dEtaMix, llambdaPt, atr->Pt(), lCent, lPVz, 1.};


                  ((THnSparseF*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistdPhidEtaMixLambda"))->Fill(spSigMixLambda);                             
                }
                if(masslambda > fMassLowLambda[0] && masslambda < fMassHighLambda[0]){
                  Double_t spBkgLeftMixLambda[7] = {dPhiMix, dEtaMix, llambdaPt, atr->Pt(), lCent, lPVz, 2.};


                  ((THnSparseF*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistdPhidEtaMixLambda"))->Fill(spBkgLeftMixLambda);   
                }
                if(masslambda > fMassLowLambda[2] && masslambda < fMassHighLambda[2]){
                  Double_t spBkgRightMixLambda[7] = {dPhiMix, dEtaMix, llambdaPt, atr->Pt(), lCent, lPVz, 3.};


                  ((THnSparseF*)((TList*)fOutput6->FindObject("Lambda"))->FindObject("fHistdPhidEtaMixLambda"))->Fill(spBkgRightMixLambda);            
                }
              }
            } 



            if(trig->WhichCandidate() == 2){ 
              Double_t massantilambda = trig->M();
              Double_t lantilambdaPt = trig->Pt();



if(fAnalysisMC){
                 Double_t spSigMixAntiLambda[7] = {dPhiMix, dEtaMix, lantilambdaPt, atr->Pt(), lCent, lPVz, 1.};
 ((THnSparseF*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistGendPhidEtaMixAntiLambda"))->Fill(spSigMixAntiLambda);
             
}



              if(fEffCorr){
                Double_t weightAntiLambda = fHistEffEtaPtAntiLambda->Interpolate(trig->Eta(),trig->Pt());
                Double_t weightTrack = fHistEffEtaPtTrack->Interpolate(atr->Eta(), atr->Pt());
                Double_t weight = weightAntiLambda * weightTrack;
                if(weight == 0){
                  continue;
                }


		
                  if(massantilambda > fMassLowAntiLambda[1] && massantilambda < fMassHighAntiLambda[1]){
                  Double_t spSigMixAntiLambda[7] = {dPhiMix, dEtaMix, lantilambdaPt, atr->Pt(), lCent, lPVz, 1.};

                  ((THnSparseF*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistdPhidEtaMixAntiLambda"))->Fill(spSigMixAntiLambda, 1/weight);
                }
                if(massantilambda > fMassLowAntiLambda[0] && massantilambda < fMassHighAntiLambda[0]){
                  Double_t spBkgLeftMixAntiLambda[7] = {dPhiMix, dEtaMix, lantilambdaPt, atr->Pt(), lCent, lPVz, 2.};


                  ((THnSparseF*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistdPhidEtaMixAntiLambda"))->Fill(spBkgLeftMixAntiLambda, 1/weight);
                }
                if(massantilambda > fMassLowAntiLambda[2] && massantilambda < fMassHighAntiLambda[2]){
                  Double_t spBkgRightMixAntiLambda[7] = {dPhiMix, dEtaMix, lantilambdaPt, atr->Pt(), lCent, lPVz, 3.};


                  ((THnSparseF*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistdPhidEtaMixAntiLambda"))->Fill(spBkgRightMixAntiLambda, 1/weight);
                } 
              }    
              else{
                if(massantilambda > fMassLowAntiLambda[1] && massantilambda < fMassHighAntiLambda[1]){
                  Double_t spSigMixAntiLambda[7] = {dPhiMix, dEtaMix, lantilambdaPt, atr->Pt(), lCent, lPVz, 1.};

                  ((THnSparseF*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistdPhidEtaMixAntiLambda"))->Fill(spSigMixAntiLambda);
                }
                if(massantilambda > fMassLowAntiLambda[0] && massantilambda < fMassHighAntiLambda[0]){
                  Double_t spBkgLeftMixAntiLambda[7] = {dPhiMix, dEtaMix, lantilambdaPt, atr->Pt(), lCent, lPVz, 2.};



                  ((THnSparseF*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistdPhidEtaMixAntiLambda"))->Fill(spBkgLeftMixAntiLambda);
                }
                if(massantilambda > fMassLowAntiLambda[2] && massantilambda < fMassHighAntiLambda[2]){
                  Double_t spBkgRightMixAntiLambda[7] = {dPhiMix, dEtaMix, lantilambdaPt, atr->Pt(), lCent, lPVz, 3.};




                  ((THnSparseF*)((TList*)fOutput7->FindObject("AntiLambda"))->FindObject("fHistdPhidEtaMixAntiLambda"))->Fill(spBkgRightMixAntiLambda);
                } 
              }
            }


          } // end of mixing track loop
        }//end of loop of selected V0particles
      }// end of loop of mixing events
    }//end if pool
    
    TObjArray* tracksClone = (TObjArray*) selectedTracks->Clone();
    tracksClone->SetOwner(kTRUE);
    pool->UpdatePool(tracksClone);

    PostData(1, fOutput);
}

//====================================================================================
Bool_t AliAnalysisTaskV0ChCorrelationsys::IsGoodPrimaryTrack(const AliAODTrack *t)
{
  // Pseudorapidity cut
  if (TMath::Abs(t->Eta())>fTrackEta) return kFALSE;//0.8

  //768-hybrid tracks
  if (!t->TestFilterBit(fFilterBit)) return kFALSE; 
  /*
  // Minimum number of clusters//i closed compare to another code 
  Float_t nCrossedRowsTPC = t->GetTPCClusterInfo(2,1);
  if (nCrossedRowsTPC < fAssocNcls) return kFALSE;
*/
  return kTRUE;
}
//====================================================================================
Bool_t AliAnalysisTaskV0ChCorrelationsys::IsGoodDaughterTrack( const AliAODTrack *t)
{
  // Pseudorapidity cut   
  if (TMath::Abs(t->Eta()) > 0.8) return kFALSE;
  
  //pt cut
  //if (t->Pt() < fV0DaughterPtMinCut) return kFALSE;//2017-8-28

  //TPC refit
  if (!t->IsOn(AliAODTrack::kTPCrefit)) return kFALSE;

  // Minimum number of clusters
  Float_t nCrossedRowsTPC = t->GetTPCClusterInfo(2,1);
  if (nCrossedRowsTPC < fV0DaughterTrackTPCCluster) return kFALSE;
   
  // findable clusters
  Int_t findable=t->GetTPCNclsF(); 
  if (findable <= 0) return kFALSE;

  if (nCrossedRowsTPC/findable < fNCrossedRowsTPCfindable) return kFALSE;

  return kTRUE;
}
    
//====================================================================================
Bool_t AliAnalysisTaskV0ChCorrelationsys::IsGoodV0(AliAODv0* aodV0  , Int_t oSta)
{
  if (!aodV0){
     AliError(Form("ERROR: Could not retrieve aodV0"));
     return kFALSE;
    }

   // Offline reconstructed V0 only
 // if(aodV0->GetOnFlyStatus()) return kFALSE;


//======================
// Offline reconstructed V0 only
    if (oSta==1) {if (aodV0->GetOnFlyStatus()) return kFALSE;}////offline
    if (oSta==3) {if (!aodV0->GetOnFlyStatus()) return kFALSE;}  // "on fly" during the tracking
   
	if (oSta==2) 
	{
		if (aodV0->GetOnFlyStatus()) 
		{
			return kTRUE;
		} else {
			return kFALSE;
		}
	}
    if (oSta==4) 
	{
		if (!aodV0->GetOnFlyStatus()) 
		{
			return kTRUE;
		} else {
			return kFALSE;
		}
	}

  //DCA of daughter track to Primary Vertex
  Float_t xyn = aodV0->DcaNegToPrimVertex();
  Float_t xyp = aodV0->DcaPosToPrimVertex();

  //DCA of daughter tracks 
  Double_t dDCA = aodV0->DcaV0Daughters();

  //Cosinus of pointing angle
 // Double_t dCPA = aodV0->CosPointingAngle(fBestPrimaryVtxPos);//lCPA
 

 Double_t lCPA = aodV0->CosPointingAngle(fBestPrimaryVtxPos);//lCPA
  //Fiducial volume cut
  Double_t xyz[3]; aodV0->GetSecondaryVtx(xyz);
  Double_t r2 = xyz[0]*xyz[0] + xyz[1]*xyz[1];

  Double_t R = aodV0->RadiusV0(); 

  // Get daughters and check them
  AliAODTrack *myTrackNegTest=dynamic_cast<AliAODTrack *>(aodV0->GetDaughter(1));
  AliAODTrack *myTrackPosTest=dynamic_cast<AliAODTrack *>(aodV0->GetDaughter(0));

  if (!myTrackPosTest || !myTrackNegTest) {
     Printf("strange analysis::UserExec:: Error:Could not retreive one of the daughter track\n");
     return kFALSE;
  }

  //Decay length
  Double_t dDL = aodV0->DecayLengthV0(fBestPrimaryVtxPos);

  //DCA or V0 to prim. vertex
  Double_t dDCA2PV = aodV0->DcaV0ToPrimVertex();

  Double_t dQT=aodV0->PtArmV0();
  Double_t dALPHA=aodV0->AlphaV0();
  Double_t dPT=aodV0->Pt();
  Double_t dPhi=aodV0->Phi();
  Double_t dEta=aodV0->PseudoRapV0();
 
  //Minimum pt of daughters
  Double_t  lMomPos[3] = {999,999,999};
  Double_t  lMomNeg[3] = {999,999,999};

  lMomPos[0] = aodV0->MomPosX();
  lMomPos[1] = aodV0->MomPosY();
  lMomPos[2] = aodV0->MomPosZ();

  lMomNeg[0] = aodV0->MomNegX();
  lMomNeg[1] = aodV0->MomNegY();
  lMomNeg[2] = aodV0->MomNegZ();

  Double_t lPtPos = TMath::Sqrt(lMomPos[0]*lMomPos[0] + lMomPos[1]*lMomPos[1]);
  Double_t lPtNeg = TMath::Sqrt(lMomNeg[0]*lMomNeg[0] + lMomNeg[1]*lMomNeg[1]);

  Double_t cutMinPtDaughter = 0.160;
  if (lPtPos<cutMinPtDaughter || lPtNeg<cutMinPtDaughter) return kFALSE;


  ((TH1F*)((TList*)fOutput4->FindObject("V0"))->FindObject("BfV0Pt"))->Fill(dPT);
  ((TH1F*)((TList*)fOutput4->FindObject("V0"))->FindObject("BfV0Phi"))->Fill(dPhi);
  ((TH1F*)((TList*)fOutput4->FindObject("V0"))->FindObject("BfV0Eta"))->Fill(dEta);
  ((TH1F*)((TList*)fOutput4->FindObject("V0"))->FindObject("Bfxyn"))->Fill(xyn);
  ((TH2F*)((TList*)fOutput4->FindObject("V0"))->FindObject("BfD0vsPt"))->Fill(xyn,dPT);
  ((TH1F*)((TList*)fOutput4->FindObject("V0"))->FindObject("Bfxyp"))->Fill(xyp);
  ((TH2F*)((TList*)fOutput4->FindObject("V0"))->FindObject("BfD0vsPt"))->Fill(xyp,dPT);
  ((TH1F*)((TList*)fOutput4->FindObject("V0"))->FindObject("BfDCA"))->Fill(dDCA);
  ((TH2F*)((TList*)fOutput4->FindObject("V0"))->FindObject("BfDCAvsPt"))->Fill(dDCA,dPT);
  ((TH1F*)((TList*)fOutput4->FindObject("V0"))->FindObject("BfCPA"))->Fill(lCPA);
  ((TH2F*)((TList*)fOutput4->FindObject("V0"))->FindObject("BfCPAvsPt"))->Fill(lCPA,dPT);
  ((TH1F*)((TList*)fOutput4->FindObject("V0"))->FindObject("Bfr2"))->Fill(r2);
  ((TH1F*)((TList*)fOutput4->FindObject("V0"))->FindObject("BfDL"))->Fill(dDL);
  ((TH2F*)((TList*)fOutput4->FindObject("V0"))->FindObject("BfDLvsPt"))->Fill(dDL,dPT);
  ((TH1F*)((TList*)fOutput4->FindObject("V0"))->FindObject("BfDCA2PV"))->Fill(dDCA2PV);
  ((TH2F*)((TList*)fOutput4->FindObject("V0"))->FindObject("BfDCA2PVvsPt"))->Fill(dDCA2PV, dPT);
  ((TH2F*)((TList*)fOutput4->FindObject("V0"))->FindObject("BfAP"))->Fill(dALPHA,dQT);
  ((TH3F*)((TList*)fOutput4->FindObject("V0"))->FindObject("BfAPvsPt"))->Fill(dALPHA,dQT,dPT);

  if (dDCA>fDCAV0DaughtersMax ) return kFALSE;
  if (lCPA<fCosPointingAngleMin) return kFALSE;

  if (r2 < f2DFiducialMin*f2DFiducialMin) return kFALSE;
   if (r2>100*100) return kFALSE;
  if ( !(IsGoodDaughterTrack(myTrackPosTest)) ||!(IsGoodDaughterTrack(myTrackNegTest)) ) return kFALSE;
  if (myTrackNegTest->Charge() == myTrackPosTest->Charge()) return kFALSE;
  return kTRUE;
}
//====================================================================================
Bool_t AliAnalysisTaskV0ChCorrelationsys::IsMCV0Primary(AliAODv0 *v0, Int_t specie){
  Int_t motherLabel = IsMcV0(v0, specie);
  if(motherLabel == -1) return kFALSE;
  AliAODMCParticle * part
    =  dynamic_cast<AliAODMCParticle*>(fMCArray->At(TMath::Abs(motherLabel)));
  if(!part) return kFALSE;
  if(part->IsPhysicalPrimary()) return kTRUE;
  return kFALSE;
}
//====================================================================================
Bool_t AliAnalysisTaskV0ChCorrelationsys::IsMCV0FromXi(AliAODv0 *v0, Int_t specie){
  Int_t motherLabel = IsMcV0(v0, specie);
  if(motherLabel == -1) return kFALSE;
  
  AliAODMCParticle * part
    =  dynamic_cast<AliAODMCParticle*>(fMCArray->At(TMath::Abs(motherLabel)));
  
  if(!part) return kFALSE;
  if(part->IsPhysicalPrimary()) return kFALSE;
  
  Int_t grandMotherLabel = part->GetMother();
  if(grandMotherLabel == -1) return kFALSE;
  AliAODMCParticle * partXi
    =  dynamic_cast<AliAODMCParticle*>(fMCArray->At(TMath::Abs(grandMotherLabel)));
  if(!partXi) return kFALSE;

  switch(specie){
  case 1: 
    if(partXi->GetPdgCode() == 3322
       || partXi->GetPdgCode() == 3312)
      return kTRUE;
    break;
  case 2:
    if(partXi->GetPdgCode() == -3322
       || partXi->GetPdgCode() == -3312)
      return kTRUE;
    break;
  default:
    break;
  }
  
 return kFALSE;
}
//=================================================================================
Int_t AliAnalysisTaskV0ChCorrelationsys::IsMcV0(AliAODv0 *v0, Int_t specie) const{
  //                                                                           
  // check if the passed V0 is associated to a MC one,                         
  //   and returns the corresponding geant label.                              
  // returns -1 if the V0 is fake (i.e. label<0).                              
  //       MCV0

  AliAODVertex *vtx=v0->GetSecondaryVtx();
  AliAODTrack *nTrack = (AliAODTrack*)vtx->GetDaughter(1);
  AliAODTrack *pTrack = (AliAODTrack*)vtx->GetDaughter(0);

  if (!nTrack || !pTrack) return -1 ;

  Int_t nlab  = nTrack->GetLabel() ;
  Int_t plab  = pTrack->GetLabel() ;

  if (nlab == -1 || plab == -1) return -1 ;

  return GetV0Label(nlab,plab, specie) ;
}
 //=======================================================================================
Int_t AliAnalysisTaskV0ChCorrelationsys::GetV0Label(Int_t lab1, Int_t lab2, Int_t specie) const{                                                                          
  
  // returns the label of the V0, given the labels of the 2 daughter tracks    
  // returns -1 if the V0 is fake                                              
         
  AliAODMCParticle *mcPart1
    =  dynamic_cast<AliAODMCParticle*>(fMCArray->At(TMath::Abs(lab1)));
  AliAODMCParticle *mcPart2
    =  dynamic_cast<AliAODMCParticle*>(fMCArray->At(TMath::Abs(lab2)));

  Int_t part1MotherLab = mcPart1->GetMother();
  Int_t part2MotherLab = mcPart2->GetMother();

  if (part1MotherLab==-1 || part2MotherLab==-1) return -1 ;
  if (part1MotherLab != part2MotherLab )        return -1 ;

   AliAODMCParticle *mcMother
     =  dynamic_cast<AliAODMCParticle*>(fMCArray->At(TMath::Abs(part1MotherLab)));
   
   if(specie == 0){//K0S
     if(mcMother->GetPdgCode() == 310
        && ((mcPart1->GetPdgCode() == -211 && mcPart2->GetPdgCode() == 211)
             ||(mcPart1->GetPdgCode() == 211 && mcPart2->GetPdgCode() == -211))
        )
       return part1MotherLab ;
   }
   else if(specie == 1){//Lamda
     if(mcMother->GetPdgCode() == 3122
        &&((mcPart1->GetPdgCode() == -211 && mcPart2->GetPdgCode() == 2212)
          ||(mcPart1->GetPdgCode() == 2212 && mcPart2->GetPdgCode() == -211))
       )
       return part1MotherLab ;
   }
   else if(specie == 2){//Antilamda
     if(mcMother->GetPdgCode() == -3122 
        &&((mcPart1->GetPdgCode() == 211 && mcPart2->GetPdgCode() == -2212)
          ||(mcPart1->GetPdgCode() == -2212 && mcPart2->GetPdgCode() == 211))
       )
       return part1MotherLab;
   }

   return -1;
}

