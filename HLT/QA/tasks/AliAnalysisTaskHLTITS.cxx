// $Id$

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

/** @file   AliAnalysisTaskHLTITS.cxx
    @author Kalliopi Kanaki
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

#include "AliESDEvent.h"
#include "AliESDRun.h"
#include "AliESDInputHandler.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskHLTITS.h"


ClassImp(AliAnalysisTaskHLTITS)

//===========================================================================================================

AliAnalysisTaskHLTITS::AliAnalysisTaskHLTITS(const char *name)
    : 
     AliAnalysisTaskSE(name)
    ,fHistOnlITSsignal(0)
    ,fHistOfflITSsignal(0)
    ,fHistOfflITSsignalTrig(0)
    ,fHistOfflITSsignalNoTrig(0)
    ,fHistOnlITSncls(0)
    ,fHistOfflITSncls(0)
    ,fHistOfflITSnclsTrig(0)
    ,fHistOfflITSnclsNoTrig(0)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  // DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container

  DefineOutput(1, TList::Class());
}

//----------------------------------------------------------------------------------------------------

void AliAnalysisTaskHLTITS::UserCreateOutputObjects(){

  OpenFile(1);

  fOutputList = new TList();
  fOutputList->SetName(GetName());
  
  fHistOfflITSsignal = new TH1F("fHistOfflITSsignal","ITS signal (offline)",150, 0, 150); 
  fHistOnlITSsignal  = new TH1F("fHistOnlITSsignal", "ITS signal (HLT)",    150, 0, 150);

  fHistOfflITSsignalTrig   = new TH1F("fHistOfflITSsignalTrig",  "ITS signal (offline triggered)",    150, 0, 150);  
  fHistOfflITSsignalNoTrig = new TH1F("fHistOfflITSsignalNoTrig","ITS signal (offline not triggered)",150, 0, 150);   

  fHistOfflITSncls = new TH1F("fHistOfflITSncls","ITS clusters (offline)",10, 0, 10);  
  fHistOnlITSncls  = new TH1F("fHistOnlITSncls", "ITS clusters (HLT)",    10, 0, 10); 
  fHistOfflITSnclsTrig   = new TH1F("fHistOfflITSnclsTrig",  "ITS clusters (offline triggered)",    10, 0, 10);   
  fHistOfflITSnclsNoTrig = new TH1F("fHistOfflITSnclsNoTrig","ITS clusters (offline not triggered)",10, 0, 10);    

  
  fOutputList->Add(fHistOfflITSsignal);
  fOutputList->Add(fHistOnlITSsignal);  
  fOutputList->Add(fHistOfflITSsignalTrig);  
  fOutputList->Add(fHistOfflITSsignalNoTrig); 
  fOutputList->Add(fHistOfflITSncls); 
  fOutputList->Add(fHistOnlITSncls);       
  fOutputList->Add(fHistOfflITSnclsTrig);    
  fOutputList->Add(fHistOfflITSnclsNoTrig);
}

void AliAnalysisTaskHLTITS::NotifyRun(){
// This will not work if the active trigger classes change from run to run.
// Then one has to know all trigger classes before processing the data.

  AliESDEvent* evESD = dynamic_cast<AliESDEvent*>(InputEvent());
  TString trgClasses = evESD->GetESDRun()->GetActiveTriggerClasses(); 

  //fTrgClsArray = trgClasses.Tokenize(" ");
  //cout<<fTrgClsArray->GetEntries()<<endl; 
    
  evESD = NULL;
}

void AliAnalysisTaskHLTITS::UserExec(Option_t *){

  AliESDEvent *evESD = dynamic_cast<AliESDEvent*>(InputEvent());
  
  if(!evESD){
      Printf("ERROR: fESD not available");
      return;
  }
  
  AliESDEvent *evHLTESD = 0;
  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>(fInputHandler);
   
  if(esdH) evHLTESD = esdH->GetHLTEvent();
    
  if(!evHLTESD){
      Printf("ERROR: HLTesd not available");
      return;
  }

  if(evHLTESD->IsHLTTriggerFired()){   

     for(Int_t i=0; i<evESD->GetNumberOfTracks(); i++){ 
        
	 AliESDtrack *esdTrk = evESD->GetTrack(i); 
         if(esdTrk->GetITSsignal()>0){	
	    fHistOfflITSsignal->Fill(esdTrk->GetITSsignal());
	    fHistOfflITSsignalTrig->Fill(esdTrk->GetITSsignal());	
         }
      
        if(esdTrk->fITSncls>0){	
	   fHistOfflITSncls->Fill(esdTrk->fITSncls);
	   fHistOfflITSnclsTrig->Fill(esdTrk->fITSncls);
        }
     }
  
  } else {

    for(Int_t i=0; i<evESD->GetNumberOfTracks(); i++){ 
        
	AliESDtrack * esdTrk = evESD->GetTrack(i); 
        if(esdTrk->GetITSsignal()>0){	
	   fHistOfflITSsignal->Fill(esdTrk->GetITSsignal());
	   fHistOfflITSsignalNoTrig->Fill(esdTrk->GetITSsignal());
        }
      
        if(esdTrk->fITSncls>0){	
	   fHistOfflITSncls->Fill(esdTrk->fITSncls);
	   fHistOfflITSnclsNoTrig->Fill(esdTrk->fITSncls);
        }
    }    
  }

  // Post output data.
  PostData(1, fOutputList);

 }

void AliAnalysisTaskHLTITS::Terminate(Option_t *){
  /*
  Printf("Number of tracks thru CE: %d", fNtracksThruZ0);
  Printf("Number of tracks thru CE from triggered events: %d", 
	 fNtracksThruZ0Trig);
  */

  // Draw result to the screen
  // Called once at the end of the query

  //  TCanvas *c1 = new TCanvas("AliAnalysisTaskHLTITS","Trigger",10,10,510,510);
  //fHistTrigger->DrawCopy("E");
  
}
