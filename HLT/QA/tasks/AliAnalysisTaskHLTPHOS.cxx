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

/** @file   AliAnalysisTaskHLTPHOS.cxx
    @author Zhongbao Yin, Kalliopi Kanaki
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
#include "AliAnalysisTaskHLTPHOS.h"


ClassImp(AliAnalysisTaskHLTPHOS)

//===========================================================================================

AliAnalysisTaskHLTPHOS::AliAnalysisTaskHLTPHOS(const char *name)
    : 
     AliAnalysisTaskSE(name)
    ,fESDRun(0)
    ,fOutputList(0)
    ,fHistOnlTrk2PHOS(0)
    ,fHistOfflTrk2PHOS(0)
    ,fHistOfflTrk2PHOSTrig(0)
    ,fHistOfflTrk2PHOSNoTrig(0)
    ,fNevt(0)
    ,fTrgClsArray(0)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  // DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container

  DefineOutput(1, TList::Class());
}

const Float_t AliAnalysisTaskHLTPHOS::fgkEtaMin = -0.12;  
const Float_t AliAnalysisTaskHLTPHOS::fgkEtaMax =  0.12;  
const Float_t AliAnalysisTaskHLTPHOS::fgkPhiMin[5]   = {3.83972, 4.18879, 4.53786, 4.88692, 5.23599};  
const Float_t AliAnalysisTaskHLTPHOS::fgkPhiMax[5]   = {4.18879, 4.53786, 4.88692, 5.23599, 5.58505};  
const Float_t AliAnalysisTaskHLTPHOS::fgkNormX[5]    = {-0.642788, -0.34202, 0, 0.34202, 0.642788};  
const Float_t AliAnalysisTaskHLTPHOS::fgkNormY[5]    = {-0.766044, -0.939693, -1, -0.939693, -0.766044};  
const Float_t AliAnalysisTaskHLTPHOS::fgkInitPosX[5] = {-295.682, -157.329, 0, 157.329, 295.682};  
const Float_t AliAnalysisTaskHLTPHOS::fgkInitPosY[5] = {-352.38, -432.259, -460, -432.259, -352.38};

//----------------------------------------------------------------------------------------------------

void AliAnalysisTaskHLTPHOS::UserCreateOutputObjects(){
// Create histograms

  OpenFile(1);

  fOutputList = new TList();
  fOutputList->SetName(GetName());

// --------------- define histograms ---------------------//

  fHistOfflTrk2PHOS = new TH2F("fHistOfflTrk2PHOS","track intersection point in #eta and #phi (offline)",100, -0.5, 0.5, 100, 240, 340); 
  fHistOnlTrk2PHOS  = new TH2F("fHistOnlTrk2PHOS", "track intersection point in #eta and #phi (HLT)",    100, -0.5, 0.5, 100, 240, 340);
  
  fHistOfflTrk2PHOSTrig   = new TH2F("fHistOfflTrk2PHOSTrig",  "track intersection point in #eta and #phi (offline triggered)",    100, -0.5, 0.5, 100, 240, 340);  
  fHistOfflTrk2PHOSNoTrig = new TH2F("fHistOfflTrk2PHOSNoTrig","track intersection point in #eta and #phi (offline not triggered)",100, -0.5, 0.5, 100, 240, 340);   

  // -------------- add histograms to the output TList -----------------//
  
  fOutputList->Add(fHistOnlTrk2PHOS);
  fOutputList->Add(fHistOfflTrk2PHOS);  
  fOutputList->Add(fHistOfflTrk2PHOSTrig);
  fOutputList->Add(fHistOfflTrk2PHOSNoTrig);

}

void AliAnalysisTaskHLTPHOS::NotifyRun(){
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

void AliAnalysisTaskHLTPHOS::UserExec(Option_t *){

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

  Double_t b = evESD->GetMagneticField();
  
  Double_t pos[] = { 0., 0., 0.};
  AliVertex *vtx = new AliVertex(pos, 0., 0);
    
//   for(Int_t i = 0; i<fTrgClsArray->GetEntries(); i++){
//       if((evESD->GetFiredTriggerClasses()).Contains(((TObjString *)fTrgClsArray->At(i))->GetString()))  fHistTrigger->Fill(i);
//   }
// 
//   if(evHLTESD->IsHLTTriggerFired()){
//      //fHistHLTTrigger->Fill(evESD->GetTriggerMask());
//      for(Int_t i = 0; i<fTrgClsArray->GetEntries(); i++){ 
//          if((evESD->GetFiredTriggerClasses()).Contains(((TObjString *)fTrgClsArray->At(i))->GetString())) fHistHLTTrigger->Fill(i);
//      } 
//   }

  if(evHLTESD->IsHLTTriggerFired()){
     for(Int_t i = 0; i < evHLTESD->GetNumberOfTracks(); i++){
         AliESDtrack * HLTesdTrk = evHLTESD->GetTrack(i);
      
        TVector3 v; 
        if(IsInPHOS(2, HLTesdTrk, b, v)){ 
      	  Float_t phi = v.Phi(); 
      	  if(phi<0) phi += 2.*TMath::Pi(); 
      	  fHistOnlTrk2PHOS->Fill(v.Eta(), phi*TMath::RadToDeg()); 
        }else if(IsInPHOS(3, HLTesdTrk, b, v)){ 
      	  Float_t phi = v.Phi();  
      	  if(phi<0) phi += 2.*TMath::Pi();  
      	  fHistOnlTrk2PHOS->Fill(v.Eta(), phi*TMath::RadToDeg());  
        }else if( IsInPHOS(4, HLTesdTrk, b, v) ){ 
      	  Float_t phi = v.Phi();   
      	  if(phi<0) phi += 2.*TMath::Pi();   
      	  fHistOnlTrk2PHOS->Fill(v.Eta(), phi*TMath::RadToDeg()); 
        } 

/*
      if(IsInPHOS(2, HLTesdTrk, b) 
	 || IsInPHOS(3, HLTesdTrk, b)
	 || IsInPHOS(4, HLTesdTrk, b) ) cout<<"Good Trigger"<<endl;
*/      
       
     }
  }else{
    for(Int_t i = 0; i < evHLTESD->GetNumberOfTracks(); i++){ 
        AliESDtrack * HLTesdTrk = evHLTESD->GetTrack(i); 

        TVector3 v;  
        if(IsInPHOS(2, HLTesdTrk, b, v)){  
          Float_t phi = v.Phi();  
          if(phi<0) phi += 2.*TMath::Pi();  
          fHistOnlTrk2PHOS->Fill(v.Eta(), phi*TMath::RadToDeg());  
        }else if(IsInPHOS(3, HLTesdTrk, b, v)){  
          Float_t phi = v.Phi();   
          if(phi<0) phi += 2.*TMath::Pi();   
          fHistOnlTrk2PHOS->Fill(v.Eta(), phi*TMath::RadToDeg());   
        }else if( IsInPHOS(4, HLTesdTrk, b, v) ){  
          Float_t phi = v.Phi();    
          if(phi<0) phi += 2.*TMath::Pi();    
          fHistOnlTrk2PHOS->Fill(v.Eta(), phi*TMath::RadToDeg());  
        }  
    } 
  }
  
  if(evHLTESD->IsHLTTriggerFired()){

    for(Int_t i = 0; i < evESD->GetNumberOfTracks(); i++){ 
      AliESDtrack * esdTrk = evESD->GetTrack(i); 

      TVector3 v;
      if(IsInPHOS(2, esdTrk, b, v)){
	Float_t phi = v.Phi();
	if(phi<0) phi += 2.*TMath::Pi();
	fHistOfflTrk2PHOSTrig->Fill(v.Eta(), phi*TMath::RadToDeg());
	fHistOfflTrk2PHOS->Fill(v.Eta(), phi*TMath::RadToDeg());
	//cout<<"Event in PHOS 2: "<<fNevt<<endl;
      }else if(IsInPHOS(3, esdTrk, b, v)){
	Float_t phi = v.Phi(); 
        if(phi<0) phi += 2.*TMath::Pi(); 
        fHistOfflTrk2PHOSTrig->Fill(v.Eta(), phi*TMath::RadToDeg()); 
	fHistOfflTrk2PHOS->Fill(v.Eta(), phi*TMath::RadToDeg());
	//cout<<"Event in PHOS 3: "<<fNevt<<endl;
      }else if( IsInPHOS(4, esdTrk, b, v) ){
	Float_t phi = v.Phi();  
        if(phi<0) phi += 2.*TMath::Pi();  
        fHistOfflTrk2PHOSTrig->Fill(v.Eta(), phi*TMath::RadToDeg());
	fHistOfflTrk2PHOS->Fill(v.Eta(), phi*TMath::RadToDeg());
	//cout<<"Event in PHOS 4: "<<fNevt<<endl;
      }
    }
  } else {

    for(Int_t i = 0; i < evESD->GetNumberOfTracks(); i++){ 
      AliESDtrack * esdTrk = evESD->GetTrack(i); 

      TVector3 v; 
      if(IsInPHOS(2, esdTrk, b, v)){ 
        Float_t phi = v.Phi(); 
        if(phi<0) phi += 2.*TMath::Pi(); 
        fHistOfflTrk2PHOSNoTrig->Fill(v.Eta(), phi*TMath::RadToDeg()); 
	fHistOfflTrk2PHOS->Fill(v.Eta(), phi*TMath::RadToDeg());
	//cout<<"Event in PHOS 2: "<<fNevt<<endl;
      }else if(IsInPHOS(3, esdTrk, b, v)){ 
        Float_t phi = v.Phi();  
        if(phi<0) phi += 2.*TMath::Pi();  
        fHistOfflTrk2PHOSNoTrig->Fill(v.Eta(), phi*TMath::RadToDeg());  
	fHistOfflTrk2PHOS->Fill(v.Eta(), phi*TMath::RadToDeg());
	//cout<<"Event in PHOS 3: "<<fNevt<<endl;
      }else if( IsInPHOS(4, esdTrk, b, v) ){ 
        Float_t phi = v.Phi();   
        if(phi<0) phi += 2.*TMath::Pi();   
        fHistOfflTrk2PHOSNoTrig->Fill(v.Eta(), phi*TMath::RadToDeg()); 
	fHistOfflTrk2PHOS->Fill(v.Eta(), phi*TMath::RadToDeg());
	//cout<<"Event in PHOS 4: "<<fNevt<<endl;
      } 
    }    
  }

  fNevt++;
  delete vtx;

  // Post output data.
  PostData(1, fOutputList);

 }

void AliAnalysisTaskHLTPHOS::Terminate(Option_t *){
  /*
  Printf("Number of tracks thru CE: %d", fNtracksThruZ0);
  Printf("Number of tracks thru CE from triggered events: %d", 
	 fNtracksThruZ0Trig);
  */

  // Draw result to the screen
  // Called once at the end of the query

  //  TCanvas *c1 = new TCanvas("AliAnalysisTaskHLTPHOS","Trigger",10,10,510,510);
  //fHistTrigger->DrawCopy("E");
  
}

Bool_t AliAnalysisTaskHLTPHOS::IsInPHOS(Int_t iMod, AliESDtrack * trk, Float_t b, TVector3& v){

  Double_t normVector[3] = {fgkNormX[iMod], fgkNormY[iMod], 0};

  Double_t point[3] = {fgkInitPosX[iMod], fgkInitPosY[iMod], 0};
  
  if(!trk->Intersect(point, normVector, b)) return kFALSE;

  TVector3 trackPos(point);
  
  v=trackPos;

  Double_t phi = 0;
  if(trackPos.Phi() < 0) phi = trackPos.Phi() + 2*TMath::Pi();
  else phi = trackPos.Phi();

  if(trackPos.Eta() >= fgkEtaMin && 
     trackPos.Eta() <= fgkEtaMax &&
     phi >= fgkPhiMin[iMod] &&
     phi <= fgkPhiMax[iMod])
    {
      return kTRUE;
    }
  
  return kFALSE;
}
