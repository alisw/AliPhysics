// $Id: AliAnalysisTaskHLTCalo.cxx 40285 2010-04-09 14:04:51Z kkanaki $

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        *
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Zhongbao Yin <zbyin@mail.ccnu.edu.cn>,                *
//*                  Kalliopi Kanaki <Kalliopi.Kanaki@ift.uib.no>          *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliAnalysisTaskHLTCalo.cxx
    @author Zhongbao Yin, Kalliopi Kanaki, Svein Lindal
    @date
    @brief
*/

#include <iostream>

#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TString.h"
#include "TObjArray.h"
#include "TFile.h"
#include "TList.h"
#include "TRefArray.h"

#include "AliESDEvent.h"
#include "AliESDRun.h"
#include "AliESDInputHandler.h"
#include "AliESDCaloCluster.h"
#include "AliESDtrack.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskHLTCalo.h"

#include "AliHLTCaloHistoInvMass.h"
#include "AliHLTCaloHistoClusterEnergy.h"
#include "AliHLTCaloHistoMatchedTracks.h"
#include "AliHLTCaloHistoProducer.h"

#include "AliHLTGlobalTriggerDecision.h"

ClassImp(AliAnalysisTaskHLTCalo)

//===========================================================================================

AliAnalysisTaskHLTCalo::AliAnalysisTaskHLTCalo() : AliAnalysisTaskSE()
  ,fUseHLTTrigger(kFALSE)
  ,fESDRun(0)
  ,fOutputList(0)
  ,fHistOfflResiduals(NULL)
  ,fHistOnlResiduals(NULL)
  ,fHistOfflDz(NULL)
  ,fHistOnlDz(NULL)
  ,fHistOfflDxy(NULL)
  ,fHistOnlDxy(NULL)
  ,fNevt(0)
  ,fTrgClsArray(0)
  ,fGlobalHistoProdArrOff(NULL)
  ,fGlobalHistoProdArrOn(NULL)
  ,fClustersArray(NULL)
  ,fCaloName("")
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  // DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container

  //DefineOutput(1, TList::Class());
}

AliAnalysisTaskHLTCalo::AliAnalysisTaskHLTCalo(const char *name)
: 
AliAnalysisTaskSE(name)
  ,fUseHLTTrigger(kFALSE)
  ,fESDRun(0)
  ,fOutputList(0)
  ,fHistOfflResiduals(NULL)
  ,fHistOnlResiduals(NULL)
  ,fHistOfflDz(NULL)
  ,fHistOnlDz(NULL)
  ,fHistOfflDxy(NULL)
  ,fHistOnlDxy(NULL)
  ,fHistOfflResidualsPos(NULL)
  ,fHistOnlResidualsPos(NULL)
  ,fHistOfflDzPos(NULL)
  ,fHistOnlDzPos(NULL)
  ,fHistOfflDxyPos(NULL)
  ,fHistOnlDxyPos(NULL)
  ,fHistOfflResidualsNeg(NULL)
  ,fHistOnlResidualsNeg(NULL)
  ,fHistOfflDzNeg(NULL)
  ,fHistOnlDzNeg(NULL)
  ,fHistOfflDxyNeg(NULL)
  ,fHistOnlDxyNeg(NULL)
  ,fHistNclvsNcl(NULL)
  ,fHistTotEVsTotE(NULL)
  ,fNevt(0)
  ,fTrgClsArray(0)
  ,fGlobalHistoProdArrOff(NULL)
  ,fGlobalHistoProdArrOn(NULL)
  ,fClustersArray(NULL)
  ,fCaloName(name)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  // DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container

  DefineOutput(1, TList::Class());
}


//----------------------------------------------------------------------------------------------------

void AliAnalysisTaskHLTCalo::UserCreateOutputObjects(){
// Create histograms

  OpenFile(1);

  fOutputList = new TList();
  fOutputList->SetOwner();
  fOutputList->SetName(GetName());

  CreateSpecificStuff(fOutputList);
  
// --------------- define histograms ---------------------//

  fHistOfflResiduals = new TH1F("fHistOfflResiduals", "Residuals between cluster and nearest track in cm (offline)", 50, 0, 50);
  fHistOnlResiduals = new TH1F("fHistOnlResiduals", "Residuals between cluster and nearest track in cm (online)", 50, 0, 50);
  fHistOfflDxy = new TH1F("fHistOfflDxy", "Dxy between cluster and nearest track in cm (offline)", 80, -40, 40);
  fHistOnlDxy = new TH1F("fHistOnlDxy", "Dxy between cluster and nearest track in cm (online)", 80, -40, 40);
  fHistOfflDz = new TH1F("fHistOfflDz", "Dz between cluster and nearest track in cm (offline)", 80, -40, 40);
  fHistOnlDz = new TH1F("fHistOnlDz", "Dz between cluster and nearest track in cm (online)", 80, -40, 40);

  fHistOfflResidualsPos = new TH1F("fHistOfflResidualsPos", "Residuals between cluster and nearest track in cm (offline)Pos", 50, 0, 50);
  fHistOnlResidualsPos = new TH1F("fHistOnlResidualsPos", "Residuals between cluster and nearest track in cm (online)Pos", 50, 0, 50);
  fHistOfflDxyPos = new TH1F("fHistOfflDxyPos", "Dxy between cluster and nearest track in cm (offline)Pos", 80, -40, 40);
  fHistOnlDxyPos = new TH1F("fHistOnlDxyPos", "Dxy between cluster and nearest track in cm (online)Pos", 80, -40, 40);
  fHistOfflDzPos = new TH1F("fHistOfflDzPos", "Dz between cluster and nearest track in cm (offline)Pos", 80, -40, 40);
  fHistOnlDzPos = new TH1F("fHistOnlDzPos", "Dz between cluster and nearest track in cm (online)Pos", 80, -40, 40);

  fHistOfflResidualsNeg = new TH1F("fHistOfflResidualsNeg", "Residuals between cluster and nearest track in cm (offline)Neg", 50, 0, 50);
  fHistOnlResidualsNeg = new TH1F("fHistOnlResidualsNeg", "Residuals between cluster and nearest track in cm (online)Neg", 50, 0, 50);
  fHistOfflDxyNeg = new TH1F("fHistOfflDxyNeg", "Dxy between cluster and nearest track in cm (offline)Neg", 80, -40, 40);
  fHistOnlDxyNeg = new TH1F("fHistOnlDxyNeg", "Dxy between cluster and nearest track in cm (online)Neg", 80, -40, 40);
  fHistOfflDzNeg = new TH1F("fHistOfflDzNeg", "Dz between cluster and nearest track in cm (offline)Neg", 80, -40, 40);
  fHistOnlDzNeg = new TH1F("fHistOnlDzNeg", "Dz between cluster and nearest track in cm (online)Neg", 80, -40, 40);

  fHistNclvsNcl = new TH2F("fHistNclvsNcl", "Number of offline cl vs online cl", 100, 0, 10, 100, 0, 10);
  fHistTotEVsTotE = new TH2F("fHistTotEVsTotE", "Total energy in online cl vs total energy in offline cl", 300, 0, 150, 300, 0, 150);
  fHistTotEVsTotE->SetXTitle("Offline energy sum");
  fHistTotEVsTotE->SetYTitle("Online energy sum");
  // HLT histogram producers
  fGlobalHistoProdArrOff = new TObjArray();
  fGlobalHistoProdArrOn = new TObjArray();

  
  fGlobalHistoProdArrOff->AddLast(dynamic_cast<TObject*>(new AliHLTCaloHistoClusterEnergy(Form("%s_OFF", fCaloName.Data()))));
  fGlobalHistoProdArrOff->AddLast(dynamic_cast<TObject*>(new AliHLTCaloHistoInvMass(Form("%s_OFF", fCaloName.Data() ))));
  fGlobalHistoProdArrOff->AddLast(dynamic_cast<TObject*>(new AliHLTCaloHistoMatchedTracks(Form("%s_OFF", fCaloName.Data() ))));

  fGlobalHistoProdArrOn->AddLast(dynamic_cast<TObject*>(new AliHLTCaloHistoClusterEnergy(Form("%s_ON", fCaloName.Data() ))));
  fGlobalHistoProdArrOn->AddLast(dynamic_cast<TObject*>(new AliHLTCaloHistoInvMass(Form("%s_ON", fCaloName.Data() ))));
  fGlobalHistoProdArrOn->AddLast(dynamic_cast<TObject*>(new AliHLTCaloHistoMatchedTracks(Form("%s_ON", fCaloName.Data() ))));

  fClustersArray = new TRefArray();

  // -------------- add histograms to the output TList -----------------//
  
  fOutputList->Add(fHistOfflResiduals);
  fOutputList->Add(fHistOnlResiduals);
  fOutputList->Add(fHistOfflDz);
  fOutputList->Add(fHistOnlDz);
  fOutputList->Add(fHistOfflDxy);
  fOutputList->Add(fHistOnlDxy);


  fOutputList->Add(fHistOfflResidualsPos);
  fOutputList->Add(fHistOnlResidualsPos);
  fOutputList->Add(fHistOfflDzPos);
  fOutputList->Add(fHistOnlDzPos);
  fOutputList->Add(fHistOfflDxyPos);
  fOutputList->Add(fHistOnlDxyPos);


  fOutputList->Add(fHistOfflResidualsNeg);
  fOutputList->Add(fHistOnlResidualsNeg);
  fOutputList->Add(fHistOfflDzNeg);
  fOutputList->Add(fHistOnlDzNeg);
  fOutputList->Add(fHistOfflDxyNeg);
  fOutputList->Add(fHistOnlDxyNeg);


  fOutputList->Add(fHistTotEVsTotE);
  fOutputList->Add(fHistNclvsNcl);

  for(int ip = 0; ip < fGlobalHistoProdArrOff->GetEntriesFast(); ip++) {
    AliHLTCaloHistoProducer *prod = dynamic_cast<AliHLTCaloHistoProducer*>(fGlobalHistoProdArrOff->At(ip));
    fOutputList->AddAll(prod->GetHistograms());
  }

  for(int ip = 0; ip < fGlobalHistoProdArrOn->GetEntriesFast(); ip++) {
    AliHLTCaloHistoProducer *prod = dynamic_cast<AliHLTCaloHistoProducer*>(fGlobalHistoProdArrOn->At(ip));
    fOutputList->AddAll(prod->GetHistograms());
  }
   
  PostData(1, fOutputList);
}

void AliAnalysisTaskHLTCalo::NotifyRun(){
// This will not work if the active trigger classes change from run to run.
// Then one has to know all trigger classes before processing the data.

  AliESDEvent* evESD = dynamic_cast<AliESDEvent*>(InputEvent());
  TString trgClasses = evESD->GetESDRun()->GetActiveTriggerClasses(); 
 
  fTrgClsArray = trgClasses.Tokenize(" ");
  //cout<<fTrgClsArray->GetEntries()<<endl; 
    
//   for(Int_t i = 0; i < fTrgClsArray->GetEntries(); i++){ 
//     TString str = ((TObjString *)fTrgClsArray->At(i))->GetString(); 
//     (fHistTrigger->GetXaxis())->SetBinLabel(i+1, str.Data()); 
//     (fHistHLTTrigger->GetXaxis())->SetBinLabel(i+1, str.Data()); 
//   } 
  
  evESD = NULL;
}

void AliAnalysisTaskHLTCalo::UserExec(Option_t *){


  AliESDEvent* evESD = dynamic_cast<AliESDEvent*>(InputEvent());
  
  if (!evESD) {
    Printf("ERROR: fESD not available");
    return;
  }
  
  AliESDEvent* evHLTESD = 0;
  AliESDInputHandler* esdH = dynamic_cast<AliESDInputHandler*>(fInputHandler);
   
  if(esdH) evHLTESD = esdH->GetHLTEvent();
    
  if(!evHLTESD){
      Printf("ERROR: HLTesd not available");
      return;
  }

  // To check if anything was triggered in HLT.
  if(fUseHLTTrigger){  
    if (!((AliHLTGlobalTriggerDecision*)evHLTESD->GetHLTTriggerDecision())->Result()){
      return;
      // Process all and any events that were triggered by HLT.
    } 
  }

  DoSpecificStuff(evESD, evHLTESD);
  
  Double_t offE = 0.0;
  Double_t onE = 0.0;


  for(int icl = 0; icl < evESD->GetNumberOfCaloClusters(); icl++) {
    AliESDCaloCluster * cluster = evESD->GetCaloCluster(icl);
    if(IsThisDetector(cluster)) {
     
      offE += cluster->E();
      fHistOfflResiduals->Fill(cluster->GetEmcCpvDistance());
      fHistOfflDz->Fill(cluster->GetTrackDz());
      fHistOfflDxy->Fill(cluster->GetTrackDx());
      
     
      TArrayI* matchedTracks = cluster->GetTracksMatched();
     
      if (matchedTracks) {
	if (matchedTracks->At(0) > -1) {
	  
	  Int_t trackID = matchedTracks->At(0);
	  
	  AliESDtrack * track = evESD->GetTrack(trackID);

	  if(track) {
	    if (track->Charge() > 0) {
	  
	      fHistOfflResidualsPos->Fill(cluster->GetEmcCpvDistance());;
	      fHistOfflDzPos->Fill(cluster->GetTrackDz());
	      fHistOfflDxyPos->Fill(cluster->GetTrackDx());
	      
	    } else if (track->Charge() < 0) {
	  
	      fHistOfflResidualsNeg->Fill(cluster->GetEmcCpvDistance());;
	      fHistOfflDzNeg->Fill(cluster->GetTrackDz());
	      fHistOfflDxyNeg->Fill(cluster->GetTrackDx());
	      
	    } else {
	      cout <<"BALLE wtf!!"<<endl;
	    }
	    
	  } else {
	    cout << "BALLE no track!"<<endl;
	  }
	}
      } else {
	cout << "BALLE no array"<<endl;
      }
    }
  }

  

  for(int icl = 0; icl < evHLTESD->GetNumberOfCaloClusters(); icl++) {
    AliESDCaloCluster * cluster = evHLTESD->GetCaloCluster(icl);
    if(IsThisDetector(cluster)) {
      onE += cluster->E();
      fHistOnlResiduals->Fill(evHLTESD->GetCaloCluster(icl)->GetEmcCpvDistance());
      fHistOnlDxy->Fill(evHLTESD->GetCaloCluster(icl)->GetTrackDx());
      fHistOnlDz->Fill(evHLTESD->GetCaloCluster(icl)->GetTrackDz());

      TArrayI* matchedTracks = cluster->GetTracksMatched();
     
      if (matchedTracks) {
	if (matchedTracks->At(0) > -1) {
	  
	  Int_t trackID = matchedTracks->At(0);
	  
	  AliESDtrack * track = evHLTESD->GetTrack(trackID);

	  if(track) {
	    if (track->Charge() > 0) {
	  
	      fHistOnlResidualsPos->Fill(cluster->GetEmcCpvDistance());;
	      fHistOnlDzPos->Fill(cluster->GetTrackDz());
	      fHistOnlDxyPos->Fill(cluster->GetTrackDx());
	      
	    } else if (track->Charge() < 0) {
	  
	      fHistOnlResidualsNeg->Fill(cluster->GetEmcCpvDistance());;
	      fHistOnlDzNeg->Fill(cluster->GetTrackDz());
	      fHistOnlDxyNeg->Fill(cluster->GetTrackDx());
	      
	    } else {
	      cout <<"BALLE wtf!!"<<endl;
	    }
	    
	  } else {
	    cout << "BALLE no track!"<<endl;
	  }
	}
      } else {
	//cout << "BALLE no array"<<endl;
      }





    }
  }

  Int_t nc = GetClusters(evESD, fClustersArray);
  for(int ip = 0; ip < fGlobalHistoProdArrOff->GetEntriesFast(); ip++) {
    AliHLTCaloHistoProducer *prod = dynamic_cast<AliHLTCaloHistoProducer*>(fGlobalHistoProdArrOff->At(ip));
    prod->FillHistograms(nc, fClustersArray);
  }
 
  Int_t nOnc = GetClusters(evHLTESD, fClustersArray);
  for(int ip = 0; ip < fGlobalHistoProdArrOn->GetEntriesFast(); ip++) {
    AliHLTCaloHistoProducer *prod = dynamic_cast<AliHLTCaloHistoProducer*>(fGlobalHistoProdArrOn->At(ip));
    prod->FillHistograms(nOnc, fClustersArray);
  }
  
  fHistNclvsNcl->Fill(nc, nOnc);
  fHistTotEVsTotE->Fill(offE, onE);


  fNevt++;

  // Post output data.
  PostData(1, fOutputList);

}

void AliAnalysisTaskHLTCalo::Terminate(Option_t *){

  
  
  // Printf("Number of tracks thru CE: %d", fNtracksThruZ0);
  // Printf("Number of tracks thru CE from triggered events: %d", 
  // 	 fNtracksThruZ0Trig);
  

  // Draw result to the screen
  // Called once at the end of the query

  //  TCanvas *c1 = new TCanvas("AliAnalysisTaskHLTCalo","Trigger",10,10,510,510);
  //fHistTrigger->DrawCopy("E");
  
}
