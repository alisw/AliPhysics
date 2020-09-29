/**************************************************************************
* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
*									 *
* Authors: Friederike Bock						 *									
* March 2019: A. Marin, TrackQA                                          *
* Version 1.0								 *
*									 *
* Permission to use, copy, modify and distribute this software and its	 *
* documentation strictly for non-commercial purposes is hereby granted	 *
* without fee, provided that the above copyright notice appears in all	 *
* copies and that both the copyright notice and this permission notice	 *
* appear in the supporting documentation. The authors make no claims	 *
* about the suitability of this software for any purpose. It is		
* provided "as is" without express or implied warranty.			 *
**************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// QA Task for V0 Reader V1
//---------------------------------------------
////////////////////////////////////////////////

#include "AliAnalysisTaskTrackQA.h"
#include "TChain.h"
#include "AliAnalysisManager.h"
#include "TParticle.h"
#include "TVectorF.h"
#include "AliPIDResponse.h"
#include "AliESDtrackCuts.h"
#include "TFile.h"

class iostream;

using namespace std;

ClassImp(AliAnalysisTaskTrackQA)

AliAnalysisTaskTrackQA::AliAnalysisTaskTrackQA() : AliAnalysisTaskSE(),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fPionSelector(NULL),
  fPionSelectorName("PionSelector"),
  fKaonSelector(NULL),
  fKaonSelectorName("KaonSelector"),
  fProtonSelector(NULL),
  fProtonSelectorName("ProtonSelector"),
  fDeuteronSelector(NULL),
  fDeuteronSelectorName("DeuteronSelector"),
  fEventCutArray(NULL),
  fPionCutArray(nullptr),
  fKaonCutArray(nullptr),
  fProtonCutArray(nullptr),
  fDeuteronCutArray(nullptr),
  fSelectorNegPionIndex(0),
  fSelectorPosPionIndex(0),
  fSelectorNegKaonIndex(0),
  fSelectorPosKaonIndex(0),
  fSelectorNegProtonIndex(0),
  fSelectorPosProtonIndex(0),
  fSelectorNegDeuteronIndex(0),
  fSelectorPosDeuteronIndex(0),
  fPosPionCandidates(nullptr),
  fNegPionCandidates(nullptr),
  fPosKaonCandidates(nullptr),
  fNegKaonCandidates(nullptr),
  fPosProtonCandidates(nullptr),
  fNegProtonCandidates(nullptr),
  fPosDeuteronCandidates(nullptr),
  fNegDeuteronCandidates(nullptr),
  fCutFolder(NULL),
  fESDList(NULL),
  fTrueList(NULL),
  fMCList(NULL),
  fOutputContainer(NULL),
  fNESDtracksEta08(0),
  fNContrVtx(0),
  hESDPosPionDCAxy(NULL),
  hESDNegPionDCAxy(NULL),
  hESDPosKaonDCAxy(NULL),
  hESDNegKaonDCAxy(NULL),
  hESDPosProtonDCAxy(NULL),
  hESDNegProtonDCAxy(NULL),
  hESDPosDeuteronDCAxy(NULL),
  hESDNegDeuteronDCAxy(NULL),  
  hESDPosPionDCAz(NULL),
  hESDNegPionDCAz(NULL),
  hESDPosKaonDCAz(NULL),
  hESDNegKaonDCAz(NULL),
  hESDPosProtonDCAz(NULL),
  hESDNegProtonDCAz(NULL),
  hESDPosDeuteronDCAz(NULL),
  hESDNegDeuteronDCAz(NULL),
  fIsHeavyIon(0),
  fIsMC(0),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fnCuts(0),
  fiCut(0),
  hNEvents(NULL),
  hNGoodESDTracksEta08(NULL)
{

}


//________________________________________________________________________
AliAnalysisTaskTrackQA::AliAnalysisTaskTrackQA(const char *name) : AliAnalysisTaskSE(name),
  fV0Reader(NULL),
  fV0ReaderName("V0ReaderV1"),
  fPionSelector(NULL),
  fPionSelectorName("PionSelector"),
  fKaonSelector(NULL),
  fKaonSelectorName("KaonSelector"),
  fProtonSelector(NULL),
  fProtonSelectorName("ProtonSelector"),
  fDeuteronSelector(NULL),
  fDeuteronSelectorName("DeuteronSelector"),
  fEventCutArray(NULL),
  fPionCutArray(nullptr),
  fKaonCutArray(nullptr),
  fProtonCutArray(nullptr),
  fDeuteronCutArray(nullptr),
  fSelectorNegPionIndex(0),
  fSelectorPosPionIndex(0),
  fSelectorNegKaonIndex(0),
  fSelectorPosKaonIndex(0),
  fSelectorNegProtonIndex(0),
  fSelectorPosProtonIndex(0),
  fSelectorNegDeuteronIndex(0),
  fSelectorPosDeuteronIndex(0),
  fPosPionCandidates(nullptr),
  fNegPionCandidates(nullptr),
  fPosKaonCandidates(nullptr),
  fNegKaonCandidates(nullptr),
  fPosProtonCandidates(nullptr),
  fNegProtonCandidates(nullptr),
  fPosDeuteronCandidates(nullptr),
  fNegDeuteronCandidates(nullptr),
  fCutFolder(NULL),
  fESDList(NULL),
  fTrueList(NULL),
  fMCList(NULL),
  fOutputContainer(NULL),
  fNESDtracksEta08(0),
  fNContrVtx(0),
  hESDPosPionDCAxy(NULL),
  hESDNegPionDCAxy(NULL),
  hESDPosKaonDCAxy(NULL),
  hESDNegKaonDCAxy(NULL),
  hESDPosProtonDCAxy(NULL),
  hESDNegProtonDCAxy(NULL),
  hESDPosDeuteronDCAxy(NULL),
  hESDNegDeuteronDCAxy(NULL),
  hESDPosPionDCAz(NULL),
  hESDNegPionDCAz(NULL),
  hESDPosKaonDCAz(NULL),
  hESDNegKaonDCAz(NULL),
  hESDPosProtonDCAz(NULL),
  hESDNegProtonDCAz(NULL),
  hESDPosDeuteronDCAz(NULL),
  hESDNegDeuteronDCAz(NULL),
    fIsHeavyIon(0),
  fIsMC(0),
  fInputEvent(NULL),
  fMCEvent(NULL),
  fnCuts(0),
  fiCut(0),
  hNEvents(NULL),
  hNGoodESDTracksEta08(NULL)

{
  // Default constructor

  DefineInput(0, TChain::Class());
  DefineOutput(1, TList::Class());
}

//________________________________________________________________________
AliAnalysisTaskTrackQA::~AliAnalysisTaskTrackQA()
{
  // default deconstructor
 
  if(fPosPionCandidates){
    delete fPosPionCandidates;
    fPosPionCandidates = 0x0;
  }

  if(fNegPionCandidates){
    delete fNegPionCandidates;
    fNegPionCandidates = 0x0;
  }

  if(fPosKaonCandidates){
    delete fPosKaonCandidates;
    fPosKaonCandidates = 0x0;
  }

  if(fNegKaonCandidates){
    delete fNegKaonCandidates;
    fNegKaonCandidates = 0x0;
  }

  if(fPosProtonCandidates){
    delete fPosProtonCandidates;
    fPosProtonCandidates = 0x0;
  }

  if(fNegProtonCandidates){
    delete fNegProtonCandidates;
    fNegProtonCandidates = 0x0;
  }

  if(fPosDeuteronCandidates){
    delete fPosDeuteronCandidates;
    fPosDeuteronCandidates = 0x0;
  }

  if(fNegDeuteronCandidates){
    delete fNegDeuteronCandidates;
    fNegDeuteronCandidates = 0x0;
  }



}
//________________________________________________________________________
void AliAnalysisTaskTrackQA::UserCreateOutputObjects()
{
  // Create User Output Objects

  if(fOutputContainer != NULL){
    delete fOutputContainer;
    fOutputContainer = NULL;
  }
  if(fOutputContainer == NULL){
    fOutputContainer = new TList();
    fOutputContainer->SetOwner(kTRUE);
  }

  // Array of current cut's pions
  fPosPionCandidates            = new TList();
  fPosPionCandidates->SetOwner(kTRUE);
  fNegPionCandidates            = new TList();
  fNegPionCandidates->SetOwner(kTRUE);

  fPosKaonCandidates            = new TList();
  fPosKaonCandidates->SetOwner(kTRUE);
  fNegKaonCandidates            = new TList();
  fNegKaonCandidates->SetOwner(kTRUE);

  fPosProtonCandidates            = new TList();
  fPosProtonCandidates->SetOwner(kTRUE);
  fNegProtonCandidates            = new TList();
  fNegProtonCandidates->SetOwner(kTRUE);

  fPosDeuteronCandidates            = new TList();
  fPosDeuteronCandidates->SetOwner(kTRUE);
  fNegDeuteronCandidates            = new TList();
  fNegDeuteronCandidates->SetOwner(kTRUE);


  fCutFolder                = new TList*[fnCuts];
  fESDList                  = new TList*[fnCuts];
  fMCList                   = new TList*[fnCuts];
  fTrueList                 = new TList*[fnCuts];

  hNEvents                  = new TH1F*[fnCuts];
  hNGoodESDTracksEta08      = new TH1F*[fnCuts];

  hESDPosPionDCAxy      = new TH2F*[fnCuts];
  hESDNegPionDCAxy      = new TH2F*[fnCuts];
  hESDPosKaonDCAxy      = new TH2F*[fnCuts];
  hESDNegKaonDCAxy      = new TH2F*[fnCuts];
  hESDPosProtonDCAxy    = new TH2F*[fnCuts];
  hESDNegProtonDCAxy    = new TH2F*[fnCuts];
  hESDPosDeuteronDCAxy  = new TH2F*[fnCuts];
  hESDNegDeuteronDCAxy  = new TH2F*[fnCuts];


  hESDPosPionDCAz       = new TH2F*[fnCuts];
  hESDNegPionDCAz       = new TH2F*[fnCuts];
  hESDPosKaonDCAz       = new TH2F*[fnCuts];
  hESDNegKaonDCAz       = new TH2F*[fnCuts];
  hESDPosProtonDCAz     = new TH2F*[fnCuts];
  hESDNegProtonDCAz     = new TH2F*[fnCuts];
  hESDPosDeuteronDCAz   = new TH2F*[fnCuts];
  hESDNegDeuteronDCAz   = new TH2F*[fnCuts];

  const Int_t kPtBins=100;
  Double_t binsPtDummy[kPtBins+1];
  for(Int_t i=1;i<kPtBins+1;i++){
    if(binsPtDummy[i-1]+0.05<1.01)
      binsPtDummy[i]=binsPtDummy[i-1]+0.05;
    else
      binsPtDummy[i]=binsPtDummy[i-1]+0.1;
  }

  for(Int_t iCut = 0; iCut<fnCuts;iCut++){
    TString cutstringEvent      = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutNumber();
    TString cutstringPion       = ((AliIdentifiedPrimaryCuts*)fPionCutArray->At(iCut))->GetCutNumber();
    TString cutstringKaon       = ((AliIdentifiedPrimaryCuts*)fKaonCutArray->At(iCut))->GetCutNumber();
    TString cutstringProton     = ((AliIdentifiedPrimaryCuts*)fProtonCutArray->At(iCut))->GetCutNumber();
    TString cutstringDeuteron   = ((AliIdentifiedPrimaryCuts*)fDeuteronCutArray->At(iCut))->GetCutNumber();
    TString fullCutString         = "";
    fullCutString               = Form("%s_%s_%s_%s_%s",cutstringEvent.Data(),cutstringPion.Data(), cutstringKaon.Data(),cutstringProton.Data(),cutstringDeuteron.Data());
    //fullCutString               = Form("%s_%s",cutstringEvent.Data(),cutstringPion.Data());
    //    TString cutstringPhoton     = ((AliConversionPhotonCuts*)fConversionCutArray->At(iCut))->GetCutNumber();
    TString nameCutFolder         = Form("Cut Number %s", fullCutString.Data());
    TString nameESDList           = Form("%s ESD histograms", fullCutString.Data());



    fCutFolder[iCut]            = new TList();
    fCutFolder[iCut]->SetName(nameCutFolder.Data());
    fCutFolder[iCut]->SetOwner(kTRUE);
    fOutputContainer->Add(fCutFolder[iCut]);

    fESDList[iCut]              = new TList();
    fESDList[iCut]->SetName(nameESDList.Data());
    fESDList[iCut]->SetOwner(kTRUE);
    fCutFolder[iCut]->Add(fESDList[iCut]);

//    Int_t nBinsR=400;
//     Int_t nBinsX=2000;
//     Int_t nBinsY=2000;
//    Int_t nBinsZ=750;
//    Int_t nBinsPhi=750;
//    Int_t nBinsEta=2000;
//    Int_t nBinsPt=400;
    Double_t xyMax =  0.5;
    Double_t xyMin = -0.5;

    hNEvents[iCut]              = new TH1F("NEvents","NEvents",14,-0.5,13.5);
    hNEvents[iCut]->GetXaxis()->SetBinLabel(1,"Accepted");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(2,"Centrality");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(3,"Miss. MC or inc. ev.");

    if (((AliConvEventCuts*)fEventCutArray->At(iCut))->IsSpecialTrigger() > 1 ){
        TString TriggerNames      = "Not Trigger: ";
        TriggerNames              = TriggerNames+ ( (AliConvEventCuts*)fEventCutArray->At(iCut))->GetSpecialTriggerName();
        hNEvents[iCut]->GetXaxis()->SetBinLabel(4,TriggerNames.Data());
    } else {
        hNEvents[iCut]->GetXaxis()->SetBinLabel(4,"Trigger");
    }
    hNEvents[iCut]->GetXaxis()->SetBinLabel(5,"Vertex Z");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(6,"Cont. Vertex");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(7,"Pile-Up");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(8,"no SDD");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(9,"no V0AND");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(10,"EMCAL/TPC problem");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(11,"rejectedForJetJetMC");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(12,"SPD hits vs tracklet");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(13,"Out-of-Bunch pileup Past-Future");
    hNEvents[iCut]->GetXaxis()->SetBinLabel(14,"Pileup V0M-TPCout Tracks");
    fESDList[iCut]->Add(hNEvents[iCut]);

    hNGoodESDTracksEta08[iCut]      = new TH1F("GoodESDTracksEta08","GoodESDTracksEta08",4000,-0.5,4000-0.5);
    fESDList[iCut]->Add(hNGoodESDTracksEta08[iCut]);

    hESDPosPionDCAxy[iCut] = new TH2F("ESD_PosPion_DCAxy","ESD_PosPion_DCAxy",800,xyMin,xyMax,kPtBins,0.,10.);
    fESDList[iCut]->Add(hESDPosPionDCAxy[iCut]);
    hESDNegPionDCAxy[iCut] = new TH2F("ESD_NegPion_DCAxy","ESD_NegPion_DCAxy",800,xyMin,xyMax,kPtBins,0.,10.);
    fESDList[iCut]->Add(hESDNegPionDCAxy[iCut]);


    hESDPosPionDCAz[iCut] = new TH2F("ESD_PosPion_DCAz","ESD_PosPion_DCAz",800,xyMin,xyMax,kPtBins,0.,10.);
    fESDList[iCut]->Add(hESDPosPionDCAz[iCut]);
    hESDNegPionDCAz[iCut] = new TH2F("ESD_NegPion_DCAz","ESD_NegPion_DCAz",800,xyMin,xyMax,kPtBins,0.,10.);
    fESDList[iCut]->Add(hESDNegPionDCAz[iCut]);


    hESDPosKaonDCAxy[iCut] = new TH2F("ESD_PosKaon_DCAxy","ESD_PosKaon_DCAxy",800,xyMin,xyMax,kPtBins,0.,10.);
    fESDList[iCut]->Add(hESDPosKaonDCAxy[iCut]);
    hESDNegKaonDCAxy[iCut] = new TH2F("ESD_NegKaon_DCAxy","ESD_NegKaon_DCAxy",800,xyMin,xyMax,kPtBins,0.,10.);
    fESDList[iCut]->Add(hESDNegKaonDCAxy[iCut]);

    hESDPosKaonDCAz[iCut] = new TH2F("ESD_PosKaon_DCAz","ESD_PosKaon_DCAz",800,xyMin,xyMax,kPtBins,0.,10.);
    fESDList[iCut]->Add(hESDPosKaonDCAz[iCut]);
    hESDNegKaonDCAz[iCut] = new TH2F("ESD_NegKaon_DCAz","ESD_NegKaon_DCAz",800,xyMin,xyMax,kPtBins,0.,10.);
    fESDList[iCut]->Add(hESDNegKaonDCAz[iCut]);

    hESDPosProtonDCAxy[iCut] = new TH2F("ESD_PosProton_DCAxy","ESD_PosProton_DCAxy",800,xyMin,xyMax,kPtBins,0.,10.);
    fESDList[iCut]->Add(hESDPosProtonDCAxy[iCut]);
    hESDNegProtonDCAxy[iCut] = new TH2F("ESD_NegProton_DCAxy","ESD_NegProton_DCAxy",800,xyMin,xyMax,kPtBins,0.,10.);
    fESDList[iCut]->Add(hESDNegProtonDCAxy[iCut]);

    hESDPosProtonDCAz[iCut] = new TH2F("ESD_PosProton_DCAz","ESD_PosProton_DCAz",800,xyMin,xyMax,kPtBins,0.,10.);
    fESDList[iCut]->Add(hESDPosProtonDCAz[iCut]);
    hESDNegProtonDCAz[iCut] = new TH2F("ESD_NegProton_DCAz","ESD_NegProton_DCAz",800,xyMin,xyMax,kPtBins,0.,10.);
    fESDList[iCut]->Add(hESDNegProtonDCAz[iCut]);

    hESDPosDeuteronDCAxy[iCut] = new TH2F("ESD_PosDeuteron_DCAxy","ESD_PosDeuteron_DCAxy",800,xyMin,xyMax,kPtBins,0.,10.);
    fESDList[iCut]->Add(hESDPosDeuteronDCAxy[iCut]);
    hESDNegDeuteronDCAxy[iCut] = new TH2F("ESD_NegDeuteron_DCAxy","ESD_NegDeuteron_DCAxy",800,xyMin,xyMax,kPtBins,0.,10.);
    fESDList[iCut]->Add(hESDNegDeuteronDCAxy[iCut]);

    hESDPosDeuteronDCAz[iCut] = new TH2F("ESD_PosDeuteron_DCAz","ESD_PosDeuteron_DCAz",800,xyMin,xyMax,kPtBins,0.,10.);
    fESDList[iCut]->Add(hESDPosDeuteronDCAz[iCut]);
    hESDNegDeuteronDCAz[iCut] = new TH2F("ESD_NegDeuteron_DCAz","ESD_NegDeuteron_DCAz",800,xyMin,xyMax,kPtBins,0.,10.);
    fESDList[iCut]->Add(hESDNegDeuteronDCAz[iCut]);



    if (fIsMC>0) {

        fMCList[iCut]               = new TList();
        fMCList[iCut]->SetName(Form("%sMC histograms",fullCutString.Data()));
        fMCList[iCut]->SetOwner(kTRUE);
        fCutFolder[iCut]->Add(fMCList[iCut]);

        fTrueList[iCut]             = new TList();
        fTrueList[iCut]->SetName(Form("%s True histograms",fullCutString.Data()));
        fTrueList[iCut]->SetOwner(kTRUE);
        fCutFolder[iCut]->Add(fTrueList[iCut]);

 
    }

  }

  fPionSelector=(AliIdentifiedPrimarySelector*)AliAnalysisManager::GetAnalysisManager()->GetTask("PionSelector");
  if(!fPionSelector){printf("Error: No PionSelector");return;} // 

  if( fPionSelector ){
    if ( ((AliIdentifiedPrimaryCuts*)fPionSelector->GetIdentifiedPrimaryCuts())->GetCutHistograms() ){
      fOutputContainer->Add( ((AliIdentifiedPrimaryCuts*)fPionSelector->GetIdentifiedPrimaryCuts())->GetCutHistograms() );
    }
  }

  fKaonSelector=(AliIdentifiedPrimarySelector*)AliAnalysisManager::GetAnalysisManager()->GetTask("KaonSelector");
  if(!fKaonSelector){printf("Error: No KaonSelector");return;} // 

  if( fKaonSelector ){
    if ( ((AliIdentifiedPrimaryCuts*)fKaonSelector->GetIdentifiedPrimaryCuts())->GetCutHistograms() ){
      fOutputContainer->Add( ((AliIdentifiedPrimaryCuts*)fKaonSelector->GetIdentifiedPrimaryCuts())->GetCutHistograms() );
    }
  }

  fProtonSelector=(AliIdentifiedPrimarySelector*)AliAnalysisManager::GetAnalysisManager()->GetTask("ProtonSelector");
  if(!fProtonSelector){printf("Error: No ProtonSelector");return;} // 

  if( fProtonSelector ){
    if ( ((AliIdentifiedPrimaryCuts*)fProtonSelector->GetIdentifiedPrimaryCuts())->GetCutHistograms() ){
      fOutputContainer->Add( ((AliIdentifiedPrimaryCuts*)fProtonSelector->GetIdentifiedPrimaryCuts())->GetCutHistograms() );
    }
  }

  fDeuteronSelector=(AliIdentifiedPrimarySelector*)AliAnalysisManager::GetAnalysisManager()->GetTask("DeuteronSelector");
  if(!fDeuteronSelector){printf("Error: No DeuteronSelector");return;} // 

  if( fDeuteronSelector ){
    if ( ((AliIdentifiedPrimaryCuts*)fDeuteronSelector->GetIdentifiedPrimaryCuts())->GetCutHistograms() ){
      fOutputContainer->Add( ((AliIdentifiedPrimaryCuts*)fDeuteronSelector->GetIdentifiedPrimaryCuts())->GetCutHistograms() );
    }
  }

 
  for(Int_t iCut = 0; iCut<fnCuts;iCut++){

    if(!((AliConvEventCuts*)fEventCutArray->At(iCut))) continue;
    if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms()){
        fCutFolder[iCut]->Add(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms());
    }

    if( fPionCutArray){
      if( ((AliIdentifiedPrimaryCuts*)fPionCutArray->At(iCut))->GetCutHistograms() ) {
        fCutFolder[iCut]->Add( ((AliIdentifiedPrimaryCuts*)fPionCutArray->At(iCut))->GetCutHistograms() );
      }
    }

    if( fKaonCutArray){
      if( ((AliIdentifiedPrimaryCuts*)fKaonCutArray->At(iCut))->GetCutHistograms() ) {
        fCutFolder[iCut]->Add( ((AliIdentifiedPrimaryCuts*)fKaonCutArray->At(iCut))->GetCutHistograms() );
      }
    }

    if( fProtonCutArray){
      if( ((AliIdentifiedPrimaryCuts*)fProtonCutArray->At(iCut))->GetCutHistograms() ) {
        fCutFolder[iCut]->Add( ((AliIdentifiedPrimaryCuts*)fProtonCutArray->At(iCut))->GetCutHistograms() );
      }
    }

    if( fDeuteronCutArray){
      if( ((AliIdentifiedPrimaryCuts*)fDeuteronCutArray->At(iCut))->GetCutHistograms() ) {
        fCutFolder[iCut]->Add( ((AliIdentifiedPrimaryCuts*)fDeuteronCutArray->At(iCut))->GetCutHistograms() );
      }
    }

  }

  PostData(1, fOutputContainer);

}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskTrackQA::Notify()
{
  // for(Int_t iCut = 0; iCut<fnCuts;iCut++){
  //   if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod && ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum() != AliConvEventCuts::kNoPeriod){
  //       ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnumExplicit(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetPeriodEnum());
  //   } else if (((AliConvEventCuts*)fEventCutArray->At(iCut))->GetPeriodEnum() == AliConvEventCuts::kNoPeriod ){
  //     ((AliConvEventCuts*)fEventCutArray->At(iCut))->SetPeriodEnum(fV0Reader->GetPeriodName());
  //   }
  // }

  return kTRUE;
}

//________________________________________________________________________
void AliAnalysisTaskTrackQA::UserExec(Option_t *){

  fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask(fV0ReaderName.Data());
  if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader

  fInputEvent = InputEvent();
  if (fInputEvent==NULL) return;

  Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
  if(fInputEvent->IsIncompleteDAQ()==kTRUE) eventQuality = 2;  // incomplete event
  // Event Not Accepted due to MC event missing or because it is incomplere or  wrong trigger for V0ReaderV1 => skip broken event/file
  if(eventQuality == 2 || eventQuality == 3){
    for(Int_t iCut = 0; iCut<fnCuts; iCut++){
      hNEvents[iCut]->Fill(eventQuality);
    }
    return;
  }

  if(fIsMC>0) fMCEvent = MCEvent();

  fPionSelector=(AliIdentifiedPrimarySelector*)AliAnalysisManager::GetAnalysisManager()->GetTask("PionSelector");
  if(!fPionSelector){printf("Error: No PionSelector");return;} // 
  fSelectorNegPionIndex = fPionSelector->GetReconstructedNegIdentifiedIndex(); // Neg Pions from default Cut
  fSelectorPosPionIndex = fPionSelector->GetReconstructedPosIdentifiedIndex(); // Pos Pions from default Cut

  fKaonSelector=(AliIdentifiedPrimarySelector*)AliAnalysisManager::GetAnalysisManager()->GetTask("KaonSelector");
  if(!fKaonSelector){printf("Error: No KaonSelector");return;} // 
  fSelectorNegKaonIndex = fKaonSelector->GetReconstructedNegIdentifiedIndex(); // Neg Kaons from default Cut
  fSelectorPosKaonIndex = fKaonSelector->GetReconstructedPosIdentifiedIndex(); // Pos Kaons from default Cut


  fProtonSelector=(AliIdentifiedPrimarySelector*)AliAnalysisManager::GetAnalysisManager()->GetTask("ProtonSelector");
  if(!fProtonSelector){printf("Error: No ProtonSelector");return;} // 
  fSelectorNegProtonIndex = fProtonSelector->GetReconstructedNegIdentifiedIndex(); // Neg Protons from default Cut
  fSelectorPosProtonIndex = fProtonSelector->GetReconstructedPosIdentifiedIndex(); // Pos Protons from default Cut

  fDeuteronSelector=(AliIdentifiedPrimarySelector*)AliAnalysisManager::GetAnalysisManager()->GetTask("DeuteronSelector");
  if(!fDeuteronSelector){printf("Error: No DeuteronSelector");return;} // 
  fSelectorNegDeuteronIndex = fDeuteronSelector->GetReconstructedNegIdentifiedIndex(); // Neg Deuterons from default Cut
  fSelectorPosDeuteronIndex = fDeuteronSelector->GetReconstructedPosIdentifiedIndex(); // Pos Deuterons from default Cut




  // ------------------- BeginEvent ----------------------------

  for(Int_t iCut = 0; iCut<fnCuts; iCut++){
    fiCut = iCut;
    Int_t eventNotAccepted = ((AliConvEventCuts*)fEventCutArray->At(iCut))->IsEventAcceptedByCut(fV0Reader->GetEventCuts(),fInputEvent,fMCEvent,fIsHeavyIon,kFALSE);
    if(eventNotAccepted){
      // cout << "event rejected due to wrong trigger: " <<eventNotAccepted << endl;
      hNEvents[iCut]->Fill(eventNotAccepted); // Check Centrality, PileUp, SDD and V0AND --> Not Accepted => eventQuality = 1
      continue;
    }

    if(eventQuality != 0){// Event Not Accepted
      // cout << "event rejected due to: " <<eventQuality << endl;
      hNEvents[iCut]->Fill(eventQuality);
      continue;
    }

    hNEvents[iCut]->Fill(eventQuality); // Should be 0 here
    fNESDtracksEta08 = CountTracks08(); // Estimate Event Multiplicity
    hNGoodESDTracksEta08[iCut]->Fill(fNESDtracksEta08);


    // if(fDoMultWeights && fIsMC > 0) {
    //   fWeightMultMC = 1.;
    //   fWeightMultMC = ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetWeightForMultiplicity(fNESDtracksEta08);
    //   hNGoodESDTracksWeightedEta08[iCut]->Fill(fNESDtracksEta08, fWeightMultMC);
    // }

    // Calculation of Multiplicity weight moved before ProcessMCPhotons
    // fWeightMultMC shuld also be inserted to input pT distributions . 

    if(fIsMC > 0){
      // Process MC Particle
      if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection() != 0){
        if(fInputEvent->IsA()==AliESDEvent::Class()){
            ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetNotRejectedParticles(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetSignalRejection(),
                            ((AliConvEventCuts*)fEventCutArray->At(iCut))->GetAcceptedHeader(),
                            fMCEvent);
        }
      }
      //   ProcessMCPhotons();
    }



    if(fInputEvent){
      if(fInputEvent->GetPrimaryVertexTracks()->GetNContributors()>0) {
        fNContrVtx = fInputEvent->GetPrimaryVertexTracks()->GetNContributors();
      } else {
        fNContrVtx = 0;
      }
    }
    ProcessPionCandidates();
    ProcessKaonCandidates();
    ProcessProtonCandidates();
    ProcessDeuteronCandidates();
    //fGammaCandidates->Clear(); // delete this cuts good gammas
  }

  //cout<<" done with the event"<<endl;

  PostData(1, fOutputContainer);
}

//________________________________________________________________________
//void AliAnalysisTaskNeutralMesonToPiPlPiMiNeutralMeson::ProcessPionCandidatesAOD(){
//}

///________________________________________________________________________
void AliAnalysisTaskTrackQA::ProcessPionCandidates(){

  Double_t magField = fInputEvent->GetMagneticField();
  if( magField  < 0.0 ){
    magField =  1.0;
  } else {
    magField =  -1.0;
  }

  Float_t dcaToVertexXYPos = -1.0;
  Float_t dcaToVertexZPos  = -1.0;
  Float_t dcaToVertexXYNeg = -1.0;
  Float_t dcaToVertexZNeg  = -1.0;


  vector<Int_t> lGoodNegPionIndexPrev(0);
  vector<Int_t> lGoodPosPionIndexPrev(0);

  for(UInt_t i = 0; i < fSelectorNegPionIndex.size(); i++){
    AliESDtrack* negPionCandidate =dynamic_cast<AliESDtrack*>(fInputEvent->GetTrack(fSelectorNegPionIndex[i]));
    if(! ((AliIdentifiedPrimaryCuts*)fPionCutArray->At(fiCut))->IdentifiedIsSelected(negPionCandidate) ) continue;
    //    lGoodNegPionIndexPrev.push_back(   fSelectorNegPionIndex[i] );

    Float_t bNeg[2];
    Float_t bCovNeg[3];
    negPionCandidate->GetImpactParameters(bNeg,bCovNeg);
    if (bCovNeg[0]<=0 || bCovNeg[2]<=0) {
      AliDebug(1, "Estimated b resolution lower or equal zero!");
      bCovNeg[0]=0; bCovNeg[2]=0;
    }
    dcaToVertexXYNeg = bNeg[0];
    dcaToVertexZNeg  = bNeg[1];
    hESDNegPionDCAxy[fiCut]->Fill(  dcaToVertexXYNeg, negPionCandidate->Pt() );
    hESDNegPionDCAz[fiCut]->Fill(   dcaToVertexZNeg,  negPionCandidate->Pt() );

    // TLorentzVector negPionforHandler;
    // negPionforHandler.SetPxPyPzE(negPionCandidate->Px(), negPionCandidate->Py(), negPionCandidate->Pz(), negPionCandidate->E());

    // AliAODConversionPhoton *negPionHandler = new AliAODConversionPhoton(&negPionforHandler);
    // fNegPionCandidates->Add(negPionHandler);
  }

  // for(UInt_t i = 0; i < lGoodNegPionIndexPrev.size(); i++){
  //   AliVTrack* negPionCandidate = dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(lGoodNegPionIndexPrev[i]));
  //   Float_t bNeg[2];
  //   Float_t bCovNeg[3];
  //   negPionCandidate->GetImpactParameters(bNeg,bCovNeg);
  //   if (bCovNeg[0]<=0 || bCovNeg[2]<=0) {
  //     AliDebug(1, "Estimated b resolution lower or equal zero!");
  //     bCovNeg[0]=0; bCovNeg[2]=0;
  //   }
  //   dcaToVertexXYNeg = bNeg[0];
  //   dcaToVertexZNeg  = bNeg[1];
  //   hESDNegPionDCAxy[fiCut]->Fill(  dcaToVertexXYNeg, negPionCandidate->Pt() );
  //   hESDNegPionDCAz[fiCut]->Fill(   dcaToVertexZNeg,  negPionCandidate->Pt() );
  // }

 
  for(UInt_t i = 0; i < fSelectorPosPionIndex.size(); i++){
    AliESDtrack* posPionCandidate = dynamic_cast<AliESDtrack*>(fInputEvent->GetTrack(fSelectorPosPionIndex[i]));
    if(! ((AliIdentifiedPrimaryCuts*)fPionCutArray->At(fiCut))->IdentifiedIsSelected(posPionCandidate) ) continue;
    //   lGoodPosPionIndexPrev.push_back(   fSelectorPosPionIndex[i]  );
    Float_t bPos[2];
    Float_t bCovPos[3];
    posPionCandidate->GetImpactParameters(bPos,bCovPos);
    if (bCovPos[0]<=0 || bCovPos[2]<=0) {
      AliDebug(1, "Estimated b resolution lower or equal zero!");
      bCovPos[0]=0; bCovPos[2]=0;
    }
    dcaToVertexXYPos = bPos[0];
    dcaToVertexZPos  = bPos[1];
    hESDPosPionDCAxy[fiCut]->Fill(  dcaToVertexXYPos, posPionCandidate->Pt());
    hESDPosPionDCAz[fiCut]->Fill(   dcaToVertexZPos,  posPionCandidate->Pt() );

    // TLorentzVector posPionforHandler;
    // posPionforHandler.SetPxPyPzE(posPionCandidate->Px(), posPionCandidate->Py(), posPionCandidate->Pz(), posPionCandidate->E());

    // AliAODConversionPhoton *posPionHandler = new AliAODConversionPhoton(&posPionforHandler);
    // fPosPionCandidates->Add(posPionHandler);
  }
 // for(UInt_t j = 0; j < lGoodPosPionIndexPrev.size(); j++){
 //    AliVTrack *posPionCandidate = dynamic_cast<AliVTrack*>(fInputEvent->GetTrack(lGoodPosPionIndexPrev[j]));
 //    Float_t bPos[2];
 //    Float_t bCovPos[3];
 //    posPionCandidate->GetImpactParameters(bPos,bCovPos);
 //    if (bCovPos[0]<=0 || bCovPos[2]<=0) {
 //      AliDebug(1, "Estimated b resolution lower or equal zero!");
 //      bCovPos[0]=0; bCovPos[2]=0;
 //    }
 //    dcaToVertexXYPos = bPos[0];
 //    dcaToVertexZPos  = bPos[1];
 //    hESDPosPionDCAxy[fiCut]->Fill(  dcaToVertexXYPos, posPionCandidate->Pt());
 //    hESDPosPionDCAz[fiCut]->Fill(   dcaToVertexZPos,  posPionCandidate->Pt() );
 //  }



}

///________________________________________________________________________
void AliAnalysisTaskTrackQA::ProcessKaonCandidates(){
  Double_t magField = fInputEvent->GetMagneticField();
  if( magField  < 0.0 ){
    magField =  1.0;
  } else {
    magField =  -1.0;
  }
  Float_t dcaToVertexXYPos = -1.0;
  Float_t dcaToVertexZPos  = -1.0;
  Float_t dcaToVertexXYNeg = -1.0;
  Float_t dcaToVertexZNeg  = -1.0;


  vector<Int_t> lGoodNegKaonIndexPrev(0);
  vector<Int_t> lGoodPosKaonIndexPrev(0);

  for(UInt_t i = 0; i < fSelectorNegKaonIndex.size(); i++){
    AliESDtrack* negKaonCandidate =dynamic_cast<AliESDtrack*>(fInputEvent->GetTrack(fSelectorNegKaonIndex[i]));
    if(! ((AliIdentifiedPrimaryCuts*)fKaonCutArray->At(fiCut))->IdentifiedIsSelected(negKaonCandidate) ) continue;
    //    lGoodNegKaonIndexPrev.push_back(   fSelectorNegKaonIndex[i] );
    Float_t bNeg[2];
    Float_t bCovNeg[3];
    negKaonCandidate->GetImpactParameters(bNeg,bCovNeg);
    if (bCovNeg[0]<=0 || bCovNeg[2]<=0) {
      AliDebug(1, "Estimated b resolution lower or equal zero!");
      bCovNeg[0]=0; bCovNeg[2]=0;
    }
    dcaToVertexXYNeg = bNeg[0];
    dcaToVertexZNeg  = bNeg[1];
    hESDNegKaonDCAxy[fiCut]->Fill(  dcaToVertexXYNeg, negKaonCandidate->Pt() );
    hESDNegKaonDCAz[fiCut]->Fill(   dcaToVertexZNeg,  negKaonCandidate->Pt() );
    // TLorentzVector negKaonforHandler;
    // negKaonforHandler.SetPxPyPzE(negKaonCandidate->Px(), negKaonCandidate->Py(), negKaonCandidate->Pz(), negKaonCandidate->E());
    // AliAODConversionPhoton *negKaonHandler = new AliAODConversionPhoton(&negKaonforHandler);
    // fNegKaonCandidates->Add(negKaonHandler);
  }
 
  for(UInt_t i = 0; i < fSelectorPosKaonIndex.size(); i++){
    AliESDtrack* posKaonCandidate = dynamic_cast<AliESDtrack*>(fInputEvent->GetTrack(fSelectorPosKaonIndex[i]));
    if(! ((AliIdentifiedPrimaryCuts*)fKaonCutArray->At(fiCut))->IdentifiedIsSelected(posKaonCandidate) ) continue;
    //   lGoodPosKaonIndexPrev.push_back(   fSelectorPosKaonIndex[i]  );
    Float_t bPos[2];
    Float_t bCovPos[3];
    posKaonCandidate->GetImpactParameters(bPos,bCovPos);
    if (bCovPos[0]<=0 || bCovPos[2]<=0) {
      AliDebug(1, "Estimated b resolution lower or equal zero!");
      bCovPos[0]=0; bCovPos[2]=0;
    }
    dcaToVertexXYPos = bPos[0];
    dcaToVertexZPos  = bPos[1];
    hESDPosKaonDCAxy[fiCut]->Fill(  dcaToVertexXYPos, posKaonCandidate->Pt());
    hESDPosKaonDCAz[fiCut]->Fill(   dcaToVertexZPos,  posKaonCandidate->Pt() );

    // TLorentzVector posKaonforHandler;
    // posKaonforHandler.SetPxPyPzE(posKaonCandidate->Px(), posKaonCandidate->Py(), posKaonCandidate->Pz(), posKaonCandidate->E());
    // AliAODConversionPhoton *posKaonHandler = new AliAODConversionPhoton(&posKaonforHandler);
    // fPosKaonCandidates->Add(posKaonHandler);
  }


}
///________________________________________________________________________
void AliAnalysisTaskTrackQA::ProcessProtonCandidates(){
  Double_t magField = fInputEvent->GetMagneticField();
  if( magField  < 0.0 ){
    magField =  1.0;
  } else {
    magField =  -1.0;
  }
  Float_t dcaToVertexXYPos = -1.0;
  Float_t dcaToVertexZPos  = -1.0;
  Float_t dcaToVertexXYNeg = -1.0;
  Float_t dcaToVertexZNeg  = -1.0;


  vector<Int_t> lGoodNegProtonIndexPrev(0);
  vector<Int_t> lGoodPosProtonIndexPrev(0);
  for(UInt_t i = 0; i < fSelectorNegProtonIndex.size(); i++){
    AliESDtrack* negProtonCandidate =dynamic_cast<AliESDtrack*>(fInputEvent->GetTrack(fSelectorNegProtonIndex[i]));
    if(! ((AliIdentifiedPrimaryCuts*)fProtonCutArray->At(fiCut))->IdentifiedIsSelected(negProtonCandidate) ) continue;
   //    lGoodNegProtonIndexPrev.push_back(   fSelectorNegProtonIndex[i] );

    Float_t bNeg[2];
    Float_t bCovNeg[3];
    negProtonCandidate->GetImpactParameters(bNeg,bCovNeg);
    if (bCovNeg[0]<=0 || bCovNeg[2]<=0) {
      AliDebug(1, "Estimated b resolution lower or equal zero!");
      bCovNeg[0]=0; bCovNeg[2]=0;
    }
    dcaToVertexXYNeg = bNeg[0];
    dcaToVertexZNeg  = bNeg[1];
    hESDNegProtonDCAxy[fiCut]->Fill(  dcaToVertexXYNeg, negProtonCandidate->Pt() );
    hESDNegProtonDCAz[fiCut]->Fill(   dcaToVertexZNeg,  negProtonCandidate->Pt() );
    // TLorentzVector negProtonforHandler;
    // negProtonforHandler.SetPxPyPzE(negProtonCandidate->Px(), negProtonCandidate->Py(), negProtonCandidate->Pz(), negProtonCandidate->E());
    // AliAODConversionPhoton *negProtonHandler = new AliAODConversionPhoton(&negProtonforHandler);
    // fNegProtonCandidates->Add(negProtonHandler);
  }



 
  for(UInt_t i = 0; i < fSelectorPosProtonIndex.size(); i++){
    AliESDtrack* posProtonCandidate = dynamic_cast<AliESDtrack*>(fInputEvent->GetTrack(fSelectorPosProtonIndex[i]));
    if(! ((AliIdentifiedPrimaryCuts*)fProtonCutArray->At(fiCut))->IdentifiedIsSelected(posProtonCandidate) ) continue;
    //   lGoodPosProtonIndexPrev.push_back(   fSelectorPosProtonIndex[i]  );
    Float_t bPos[2];
    Float_t bCovPos[3];
    posProtonCandidate->GetImpactParameters(bPos,bCovPos);
    if (bCovPos[0]<=0 || bCovPos[2]<=0) {
      AliDebug(1, "Estimated b resolution lower or equal zero!");
      bCovPos[0]=0; bCovPos[2]=0;
    }
    dcaToVertexXYPos = bPos[0];
    dcaToVertexZPos  = bPos[1];
    hESDPosProtonDCAxy[fiCut]->Fill(  dcaToVertexXYPos, posProtonCandidate->Pt());
    hESDPosProtonDCAz[fiCut]->Fill(   dcaToVertexZPos,  posProtonCandidate->Pt() );

    // TLorentzVector posProtonforHandler;
    // posProtonforHandler.SetPxPyPzE(posProtonCandidate->Px(), posProtonCandidate->Py(), posProtonCandidate->Pz(), posProtonCandidate->E());
    // AliAODConversionPhoton *posProtonHandler = new AliAODConversionPhoton(&posProtonforHandler);
    // fPosProtonCandidates->Add(posProtonHandler);
  }
}

///________________________________________________________________________
void AliAnalysisTaskTrackQA::ProcessDeuteronCandidates(){
  Double_t magField = fInputEvent->GetMagneticField();
  if( magField  < 0.0 ){
    magField =  1.0;
  } else {
    magField =  -1.0;
  }
  Float_t dcaToVertexXYPos = -1.0;
  Float_t dcaToVertexZPos  = -1.0;
  Float_t dcaToVertexXYNeg = -1.0;
  Float_t dcaToVertexZNeg  = -1.0;


  vector<Int_t> lGoodNegDeuteronIndexPrev(0);
  vector<Int_t> lGoodPosDeuteronIndexPrev(0);

  for(UInt_t i = 0; i < fSelectorNegDeuteronIndex.size(); i++){
    AliESDtrack* negDeuteronCandidate =dynamic_cast<AliESDtrack*>(fInputEvent->GetTrack(fSelectorNegDeuteronIndex[i]));
    if(! ((AliIdentifiedPrimaryCuts*)fDeuteronCutArray->At(fiCut))->IdentifiedIsSelected(negDeuteronCandidate) ) continue;
    //    lGoodNegDeuteronIndexPrev.push_back(   fSelectorNegDeuteronIndex[i] );

    Float_t bNeg[2];
    Float_t bCovNeg[3];
    negDeuteronCandidate->GetImpactParameters(bNeg,bCovNeg);
    if (bCovNeg[0]<=0 || bCovNeg[2]<=0) {
      AliDebug(1, "Estimated b resolution lower or equal zero!");
      bCovNeg[0]=0; bCovNeg[2]=0;
    }
    dcaToVertexXYNeg = bNeg[0];
    dcaToVertexZNeg  = bNeg[1];
    hESDNegDeuteronDCAxy[fiCut]->Fill(  dcaToVertexXYNeg, negDeuteronCandidate->Pt() );
    hESDNegDeuteronDCAz[fiCut]->Fill(   dcaToVertexZNeg,  negDeuteronCandidate->Pt() );
    // TLorentzVector negDeuteronforHandler;
    // negDeuteronforHandler.SetPxPyPzE(negDeuteronCandidate->Px(), negDeuteronCandidate->Py(), negDeuteronCandidate->Pz(), negDeuteronCandidate->E());
    // AliAODConversionPhoton *negDeuteronHandler = new AliAODConversionPhoton(&negDeuteronforHandler);
    // fNegDeuteronCandidates->Add(negDeuteronHandler);
  }

  for(UInt_t i = 0; i < fSelectorPosDeuteronIndex.size(); i++){
    AliESDtrack* posDeuteronCandidate = dynamic_cast<AliESDtrack*>(fInputEvent->GetTrack(fSelectorPosDeuteronIndex[i]));
    if(! ((AliIdentifiedPrimaryCuts*)fDeuteronCutArray->At(fiCut))->IdentifiedIsSelected(posDeuteronCandidate) ) continue;
    //   lGoodPosDeuteronIndexPrev.push_back(   fSelectorPosDeuteronIndex[i]  );
    Float_t bPos[2];
    Float_t bCovPos[3];
    posDeuteronCandidate->GetImpactParameters(bPos,bCovPos);
    if (bCovPos[0]<=0 || bCovPos[2]<=0) {
      AliDebug(1, "Estimated b resolution lower or equal zero!");
      bCovPos[0]=0; bCovPos[2]=0;
    }
    dcaToVertexXYPos = bPos[0];
    dcaToVertexZPos  = bPos[1];
    hESDPosDeuteronDCAxy[fiCut]->Fill(  dcaToVertexXYPos, posDeuteronCandidate->Pt());
    hESDPosDeuteronDCAz[fiCut]->Fill(   dcaToVertexZPos,  posDeuteronCandidate->Pt());

    // TLorentzVector posDeuteronforHandler;
    // posDeuteronforHandler.SetPxPyPzE(posDeuteronCandidate->Px(), posDeuteronCandidate->Py(), posDeuteronCandidate->Pz(), posDeuteronCandidate->E());
    // AliAODConversionPhoton *posDeuteronHandler = new AliAODConversionPhoton(&posDeuteronforHandler);
    // fPosDeuteronCandidates->Add(posDeuteronHandler);
  }


}


//________________________________________________________________________
Int_t AliAnalysisTaskTrackQA::CountTracks08(){

  Int_t fNumberOfESDTracks = 0;
  if(fInputEvent->IsA()==AliESDEvent::Class()){
  // Using standard function for setting Cuts

//     Bool_t selectPrimaries = kTRUE;
    static AliESDtrackCuts *EsdTrackCuts = 0x0;
    static int prevRun = -1;
    // Using standard function for setting Cuts
    Int_t runNumber = fInputEvent->GetRunNumber();
    if (prevRun!=runNumber) {
      delete EsdTrackCuts;
      EsdTrackCuts = 0;
      prevRun = runNumber;
    }
//     AliESDtrackCuts *EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010(selectPrimaries);
    if (!EsdTrackCuts) {
        // if LHC11a or earlier or if LHC13g or if LHC12a-i -> use 2010 cuts
        if( (runNumber<=146860) || (runNumber>=197470 && runNumber<=197692) || (runNumber>=172440 && runNumber<=193766) ){
            EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();

        } else if (runNumber>=209122){ // else if run2 data use 2015 PbPb cuts
            //EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2015PbPb();
            // hard coded track cuts for the moment, because AliESDtrackCuts::GetStandardITSTPCTrackCuts2015PbPb() gives spams warnings
            EsdTrackCuts = new AliESDtrackCuts();
            // TPC; clusterCut = 1, cutAcceptanceEdges = kTRUE, removeDistortedRegions = kFALSE
            EsdTrackCuts->AliESDtrackCuts::SetMinNCrossedRowsTPC(70);
            EsdTrackCuts->AliESDtrackCuts::SetMinRatioCrossedRowsOverFindableClustersTPC(0.8);
            EsdTrackCuts->SetCutGeoNcrNcl(2., 130., 1.5, 0.0, 0.0);  // only dead zone and not clusters per length
            //EsdTrackCuts->AliESDtrackCuts::SetCutOutDistortedRegionsTPC(kTRUE);
            EsdTrackCuts->AliESDtrackCuts::SetMaxChi2PerClusterTPC(4);
            EsdTrackCuts->AliESDtrackCuts::SetAcceptKinkDaughters(kFALSE);
            EsdTrackCuts->AliESDtrackCuts::SetRequireTPCRefit(kTRUE);
            // ITS; selPrimaries = 1
            EsdTrackCuts->AliESDtrackCuts::SetRequireITSRefit(kTRUE);
            EsdTrackCuts->AliESDtrackCuts::SetClusterRequirementITS(AliESDtrackCuts::kSPD,
                                        AliESDtrackCuts::kAny);
            EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexXYPtDep("0.0105+0.0350/pt^1.1");
            EsdTrackCuts->AliESDtrackCuts::SetMaxChi2TPCConstrainedGlobal(36);
            EsdTrackCuts->AliESDtrackCuts::SetMaxDCAToVertexZ(2);
            EsdTrackCuts->AliESDtrackCuts::SetDCAToVertex2D(kFALSE);
            EsdTrackCuts->AliESDtrackCuts::SetRequireSigmaToVertex(kFALSE);
            EsdTrackCuts->AliESDtrackCuts::SetMaxChi2PerClusterITS(36);

        } else { // else use 2011 version of track cuts
            EsdTrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011();
        }
        EsdTrackCuts->SetMaxDCAToVertexZ(2);
        EsdTrackCuts->SetEtaRange(-0.8, 0.8);
        EsdTrackCuts->SetPtRange(0.15);
    }

    for(Int_t iTracks = 0; iTracks < fInputEvent->GetNumberOfTracks(); iTracks++){
      AliESDtrack* curTrack = (AliESDtrack*) fInputEvent->GetTrack(iTracks);
      if(!curTrack) continue;
      if(EsdTrackCuts->AcceptTrack(curTrack) ){
        if (fMCEvent){
          if(((AliConvEventCuts*)fEventCutArray->At(fiCut))->GetSignalRejection() != 0){
                        Int_t isFromMBHeader = ((AliConvEventCuts*)fEventCutArray->At(fiCut))->IsParticleFromBGEvent(TMath::Abs(curTrack->GetLabel()), fMCEvent, fInputEvent);
            if( (isFromMBHeader < 1) ) continue;
          }
        }
        fNumberOfESDTracks++;
      }
    }
    delete EsdTrackCuts;
    EsdTrackCuts=0x0;

  }

  return fNumberOfESDTracks;

}



//________________________________________________________________________
void AliAnalysisTaskTrackQA::SetLogBinningXTH2(TH2* histoRebin){
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
void AliAnalysisTaskTrackQA::Terminate(Option_t *)
{
//    if (fStreamMaterial){
//       fStreamMaterial->GetFile()->Write();
//    }
}
