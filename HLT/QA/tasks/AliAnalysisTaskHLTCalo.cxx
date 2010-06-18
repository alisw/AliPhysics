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

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskHLTCalo.h"

#include "AliHLTCaloHistoInvMass.h"
#include "AliHLTCaloHistoClusterEnergy.h"
#include "AliHLTCaloHistoMatchedTracks.h"
#include "AliHLTCaloHistoProducer.h"

ClassImp(AliAnalysisTaskHLTCalo)

//===========================================================================================

AliAnalysisTaskHLTCalo::AliAnalysisTaskHLTCalo(const char *name)
: 
AliAnalysisTaskSE(name)
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
  ,fName(name)
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
  fOutputList->SetName(GetName());

  CreateSpecificStuff(fOutputList);
  
// --------------- define histograms ---------------------//

  fHistOfflResiduals = new TH1F("fHistOfflResiduals", "Residuals between cluster and nearest track in cm (offline)", 50, 0, 50);
  fHistOnlResiduals = new TH1F("fHistOnlResiduals", "Residuals between cluster and nearest track in cm (online)", 50, 0, 50);

  fHistOfflDxy = new TH1F("fHistOfflDxy", "Dxy between cluster and nearest track in cm (offline)", 50, 0, 50);
  fHistOnlDxy = new TH1F("fHistOnlDxy", "Dxy between cluster and nearest track in cm (online)", 50, 0, 50);

  fHistOfflDz = new TH1F("fHistOfflDz", "Dz between cluster and nearest track in cm (offline)", 50, 0, 50);
  fHistOnlDz = new TH1F("fHistOnlDz", "Dz between cluster and nearest track in cm (online)", 50, 0, 50);


  // HLT histogram producers
  fGlobalHistoProdArrOff = new TObjArray();
  fGlobalHistoProdArrOn = new TObjArray();

  //  AliHLTCaloHistoClusterEnergy * histo = new AliHLTCaloHistoClusterEnergy("EMCAL");
  
  fGlobalHistoProdArrOff->AddLast(dynamic_cast<TObject*>(new AliHLTCaloHistoClusterEnergy(Form("%s OFF", fName.Data()))));
  fGlobalHistoProdArrOff->AddLast(dynamic_cast<TObject*>(new AliHLTCaloHistoInvMass(Form("%s OFF", fName.Data() ))));
  fGlobalHistoProdArrOff->AddLast(dynamic_cast<TObject*>(new AliHLTCaloHistoMatchedTracks(Form("%s OFF", fName.Data() ))));

  fGlobalHistoProdArrOn->AddLast(dynamic_cast<TObject*>(new AliHLTCaloHistoClusterEnergy(Form("%s ON", fName.Data() ))));
  fGlobalHistoProdArrOn->AddLast(dynamic_cast<TObject*>(new AliHLTCaloHistoInvMass(Form("%s ON", fName.Data() ))));
  fGlobalHistoProdArrOn->AddLast(dynamic_cast<TObject*>(new AliHLTCaloHistoMatchedTracks(Form("%s ON", fName.Data() ))));

  fClustersArray = new TRefArray();

  // -------------- add histograms to the output TList -----------------//
  
  fOutputList->Add(fHistOfflResiduals);
  fOutputList->Add(fHistOnlResiduals);

  fOutputList->Add(fHistOfflDz);
  fOutputList->Add(fHistOnlDz);

  fOutputList->Add(fHistOfflDxy);
  fOutputList->Add(fHistOnlDxy);

  for(int ip = 0; ip < fGlobalHistoProdArrOff->GetEntriesFast(); ip++) {
    AliHLTCaloHistoProducer *prod = dynamic_cast<AliHLTCaloHistoProducer*>(fGlobalHistoProdArrOff->At(ip));
    fOutputList->AddAll(prod->GetHistograms());
  }

  for(int ip = 0; ip < fGlobalHistoProdArrOn->GetEntriesFast(); ip++) {
    AliHLTCaloHistoProducer *prod = dynamic_cast<AliHLTCaloHistoProducer*>(fGlobalHistoProdArrOn->At(ip));
    fOutputList->AddAll(prod->GetHistograms());
  }
   

  

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

  DoSpecificStuff(evESD, evHLTESD);

  for(int icl = 0; icl < evESD->GetNumberOfCaloClusters(); icl++) {
    AliESDCaloCluster * cluster = evESD->GetCaloCluster(icl);
    if(cluster && IsThisDetector(cluster)) {
      fHistOfflResiduals->Fill(evESD->GetCaloCluster(icl)->GetEmcCpvDistance());
      fHistOfflDz->Fill(evESD->GetCaloCluster(icl)->GetTrackDz());
      fHistOfflDxy->Fill(evESD->GetCaloCluster(icl)->GetTrackDx());
    }
  }

  for(int icl = 0; icl < evHLTESD->GetNumberOfCaloClusters(); icl++) {
    AliESDCaloCluster * cluster = evHLTESD->GetCaloCluster(icl);
    if(cluster && IsThisDetector(cluster)) {
      fHistOnlResiduals->Fill(evHLTESD->GetCaloCluster(icl)->GetEmcCpvDistance());
      fHistOnlDxy->Fill(evHLTESD->GetCaloCluster(icl)->GetTrackDx());
      fHistOnlDz->Fill(evHLTESD->GetCaloCluster(icl)->GetTrackDz());
    }
  }

  Int_t nc = GetClusters(evESD, fClustersArray);
  for(int ip = 0; ip < fGlobalHistoProdArrOff->GetEntriesFast(); ip++) {
    AliHLTCaloHistoProducer *prod = dynamic_cast<AliHLTCaloHistoProducer*>(fGlobalHistoProdArrOff->At(ip));
    prod->FillHistograms(nc, fClustersArray);
  }
 
  nc = GetClusters(evHLTESD, fClustersArray);
  for(int ip = 0; ip < fGlobalHistoProdArrOn->GetEntriesFast(); ip++) {
    AliHLTCaloHistoProducer *prod = dynamic_cast<AliHLTCaloHistoProducer*>(fGlobalHistoProdArrOn->At(ip));
    prod->FillHistograms(nc, fClustersArray);
  }
  

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

