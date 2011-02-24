/**************************************************************************
 * Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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

/* $Id$ */

//*************************************************************************
// 
// Class for monitoring of information used for open charm normalization
// (triggers, candles, ...)
//
// Origin: davide.caffarri@pd.infn.it
//
//*************************************************************************

#include <TChain.h>
#include <TTree.h>
#include <TBranch.h>
#include <TObjArray.h>
#include <TH1F.h>
#include <TH2F.h>  
#include <TCanvas.h>

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliESDEvent.h"
#include "AliESDInputHandler.h"
#include "AliTriggerAnalysis.h"
#include "AliLog.h"
#include "AliCounterCollection.h"
#include "AliAnalysisTaskSEMonitNorm.h"

ClassImp(AliAnalysisTaskSEMonitNorm)

//________________________________________________________________________
AliAnalysisTaskSEMonitNorm::AliAnalysisTaskSEMonitNorm(const char *name): 
  AliAnalysisTaskSE(name), 
  fESD(0), 
  fOutput(0),
  fCounterTrigg(0),
  fCounterNotTrigg(0),
  fCounterCandleTrig(0),
  fCounterCandleNotTrig(0) 

{
  // Constructor

  // Define input and output slots here
  // Output slot #0 writes into a TList container
  DefineOutput(1, TList::Class());  //My private output
  DefineOutput(2, AliCounterCollection::Class());
  DefineOutput(3, AliCounterCollection::Class());
  DefineOutput(4, AliCounterCollection::Class());
  DefineOutput(5, AliCounterCollection::Class());

}
//________________________________________________________________________
AliAnalysisTaskSEMonitNorm::~AliAnalysisTaskSEMonitNorm()
{
  // Destructor

  // histograms are in the output list and deleted when the output
  // list is deleted by the TSelector dtor

  if (fOutput) {
    delete fOutput;
    fOutput = 0;
  }

  if (fCounterTrigg) {
    delete fCounterTrigg;
    fCounterTrigg = 0;
  }

  if (fCounterNotTrigg) {
    delete fCounterNotTrigg;
    fCounterNotTrigg = 0;
 }

  if (fCounterCandleTrig) {
    delete fCounterCandleTrig;
    fCounterCandleTrig = 0;
  }
  
  if (fCounterCandleNotTrig) {
    delete fCounterCandleNotTrig;
    fCounterCandleNotTrig = 0;
    
  }
  
  
}


//________________________________________________________________________
void AliAnalysisTaskSEMonitNorm::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  // Several histograms are more conveniently managed in a TList
  fOutput = new TList;
  fOutput->SetOwner();


  //Counter for Normalization new line
  fCounterTrigg = new AliCounterCollection("Triggered_Counter");//new line
  fCounterTrigg->AddRubric("Event", "triggered/v0andONLINE/v0andOFFLINE/v0andOFF&ON/v0andON_notOFF/v0andOFF_notON/v0orOFF/V0orON");//new line
  fCounterTrigg->AddRubric("Run", 10000000);//new line
  fCounterTrigg->Init();//new line 
  
  fCounterNotTrigg = new AliCounterCollection("NotTriggered_Counter");//new line
  fCounterNotTrigg->AddRubric("Event", "!triggered/v0andONLINE/v0andOFFLINE/v0andOFF&ON/v0andON_notOFF/v0andOFF_notON/v0orOFF/v0orON");//new line
  fCounterNotTrigg->AddRubric("Run", 10000000);//new line
  fCounterNotTrigg->Init();//new line 

  fCounterCandleTrig = new AliCounterCollection("Candle_Triggered_Counter");//new line
  fCounterCandleTrig->AddRubric("Event", "triggered/candles02/2xcandles02/primVtx/TPConly");//new line
  fCounterCandleTrig->AddRubric("Run", 10000000);//new line
  fCounterCandleTrig->Init();//new line 


  fCounterCandleNotTrig = new AliCounterCollection("Candle_NotTriggered_Counter");//new line
  fCounterCandleNotTrig->AddRubric("Event", "!triggered/candles02/2xcandles02/primVtx/TPConly");//new line
  fCounterCandleNotTrig->AddRubric("Run", 10000000);//new line
  fCounterCandleNotTrig->Init();//new line 

  
  PostData(1, fOutput);
  PostData(2, fCounterTrigg);
  PostData(3, fCounterNotTrigg);
  PostData(4, fCounterCandleTrig);
  PostData(5, fCounterCandleNotTrig);
  return;
}

//________________________________________________________________________
void AliAnalysisTaskSEMonitNorm::UserExec(Option_t *) 
{
  // Main loop
  // Called for each event
  
  if (!InputEvent()) {
    Printf("ERROR: fESD not available");
    return;
  }
  AliESDEvent* esdE = (AliESDEvent*) InputEvent();
  if (!esdE) {
    Printf("ERROR: fESD not available");
    return;
  }

  // Select PHYSICS events (type=7, for data) and MC events (type=0)
  // fCheckEventType is kFALSE if fReadMC is kTRUE, hence check is skipped
  if(esdE->GetEventType()!=7 ) return;
  
  Int_t runNumber = (Int_t)esdE->GetRunNumber();
  
  //static AliTriggerAnalysis* triggerAnalysis = new AliTriggerAnalysis();
  Bool_t eventTriggered = 0;
  
  Bool_t flag02 = kFALSE;
  Bool_t flagPVT = kFALSE;
  Bool_t flagTPCOnly70 = kFALSE;

  // use response of AliPhysicsSelection
  eventTriggered = ((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected();
  
  
  AliTriggerAnalysis *trigAnalys = new AliTriggerAnalysis();
  Bool_t v0A; 
  Bool_t v0C; 
  Bool_t v0andOFF = kFALSE;
  Bool_t v0orOFF = kFALSE;
  if(esdE){
    v0C = trigAnalys->IsOfflineTriggerFired(esdE , AliTriggerAnalysis::kV0C);
    v0A = trigAnalys->IsOfflineTriggerFired(esdE , AliTriggerAnalysis::kV0A);
    if(v0A&&v0C)v0andOFF=kTRUE;   
    if (v0A||v0C)v0orOFF=kTRUE; 
  }  

  Int_t v0Con = trigAnalys->V0Trigger(esdE, AliTriggerAnalysis::kCSide, kTRUE, kFALSE);
  Int_t v0Aon = trigAnalys->V0Trigger(esdE, AliTriggerAnalysis::kASide, kTRUE, kFALSE);

  Bool_t v0andON=kFALSE;
  if((v0Aon==AliTriggerAnalysis::kV0BB)&&(v0Con==AliTriggerAnalysis::kV0BB))v0andON=kTRUE;   

  Bool_t v0orON =((v0Aon==AliTriggerAnalysis::kV0BB)||(v0Con==AliTriggerAnalysis::kV0BB));

  Int_t trkEntries = (Int_t)esdE->GetNumberOfTracks();

  for(Int_t i=0;i<trkEntries;i++){
    AliESDtrack *track=(AliESDtrack*)esdE->GetTrack(i);
    UShort_t nClusTPC=track->GetTPCNcls();
    
    //TPConly 
    if((track->GetStatus()&AliESDtrack::kTPCin)&&TMath::Abs(track->Eta()<0.9)&&nClusTPC>=70){
      if (!flagTPCOnly70){
	flagTPCOnly70=kTRUE;
      }
    }

    //candles02
    if((nClusTPC>=70)&&(track->GetStatus()&AliESDtrack::kITSrefit)&&(track->GetStatus()&AliESDtrack::kTPCrefit)){
      if((track->Pt()>0.2)){
	if(!flag02) flag02=kTRUE;
      }
    }
  }
  
  //prim Vtx
  if(esdE){
    const AliESDVertex *vtrc1 =  esdE->GetPrimaryVertex();
    if(vtrc1 && vtrc1->GetNContributors()>0)  flagPVT=kTRUE;
  }
  
  if (eventTriggered){

    fCounterTrigg->Count(Form("Event:triggered/Run:%d", runNumber));
    if (v0andON) fCounterTrigg->Count(Form("Event:v0andONLINE/Run:%d", runNumber));
    if (v0andOFF) fCounterTrigg->Count(Form("Event:v0andOFFLINE/Run:%d", runNumber));
    if (v0andON&&v0andOFF) fCounterTrigg->Count(Form("Event:v0andOFF&ON/Run:%d", runNumber));
    if (v0andON&&!v0andOFF) fCounterTrigg->Count(Form("Event:v0andON_notOFF/Run:%d", runNumber));
    if (!v0andON&&v0andOFF) fCounterTrigg->Count(Form("Event:v0andOFF_notON/Run:%d", runNumber));
    if (v0orON)  fCounterTrigg->Count(Form("Event:v0orON/Run:%d", runNumber));
    if (v0orOFF)  fCounterTrigg->Count(Form("Event:v0orOFF/Run:%d", runNumber));

    fCounterCandleTrig->Count(Form("Event:triggered/Run:%d", runNumber));
    if (flagTPCOnly70) fCounterCandleTrig->Count(Form("Event:TPConly/Run:%d", runNumber));
    if (flag02) fCounterCandleTrig->Count(Form("Event:candles02/Run:%d", runNumber));
    if (flagPVT) fCounterCandleTrig->Count(Form("Event:primVtx/Run:%d", runNumber));

  }
  else {
    fCounterNotTrigg->Count(Form("Event:!triggered/Run:%d", runNumber));
    if (v0andON) fCounterNotTrigg->Count(Form("Event:v0andONLINE/Run:%d", runNumber));
    if (v0andOFF) fCounterNotTrigg->Count(Form("Event:v0andOFFLINE/Run:%d", runNumber));
    if (v0andON&&v0andOFF) fCounterNotTrigg->Count(Form("Event:v0andOFF&ON/Run:%d", runNumber));
    if (v0andON&&!v0andOFF) fCounterNotTrigg->Count(Form("Event:v0andON_notOFF/Run:%d", runNumber));
    if (!v0andON&&v0andOFF) fCounterNotTrigg->Count(Form("Event:v0andOFF_notON/Run:%d", runNumber));
    if (v0orON)  fCounterNotTrigg->Count(Form("Event:v0orON/Run:%d", runNumber));
    if (v0orOFF)  fCounterNotTrigg->Count(Form("Event:v0orOFF/Run:%d", runNumber));

    fCounterCandleNotTrig->Count(Form("Event:!triggered/Run:%d", runNumber));
    if (flagTPCOnly70) fCounterCandleNotTrig->Count(Form("Event:TPConly/Run:%d", runNumber));
    if (flag02) fCounterCandleNotTrig->Count(Form("Event:candles02/Run:%d", runNumber));
    if (flagPVT) fCounterCandleNotTrig->Count(Form("Event:primVtx/Run:%d", runNumber));


  }
  
  
  // Post the data already here
  PostData(1, fOutput);
  PostData(2, fCounterTrigg);
  PostData(3, fCounterNotTrigg);
  PostData(4, fCounterCandleTrig);
  PostData(5, fCounterCandleNotTrig);
  
  return;
}      

//________________________________________________________________________
void AliAnalysisTaskSEMonitNorm::Terminate(Option_t *) 
{
  // Draw result to the screen
  // Called once at the end of the query
  fOutput = dynamic_cast<TList*> (GetOutputData(1));
  if (!fOutput) {     
    Printf("ERROR: fOutput not available");
    return;
  }
  
  fCounterTrigg = dynamic_cast<AliCounterCollection*> (GetOutputData(2));
  if (!fCounterTrigg) {     
    Printf("ERROR: fCounterTrigg not available");
    return;
  }
  
  fCounterNotTrigg = dynamic_cast<AliCounterCollection*> (GetOutputData(3));
  if (!fCounterNotTrigg) {     
    Printf("ERROR: fCounterNotTrigg not available");
    return;
  }
  
  fCounterCandleTrig = dynamic_cast<AliCounterCollection*> (GetOutputData(4));
  if (!fCounterCandleTrig) {     
    Printf("ERROR: fCounterCandleTrig not available");
  return;
  }

  fCounterCandleNotTrig = dynamic_cast<AliCounterCollection*> (GetOutputData(5));
  if (!fCounterCandleNotTrig) {     
    Printf("ERROR: fCounterCandleNotTrig not available");
    return;
  }
  
  return;
}

//_________________________________________________________________________
