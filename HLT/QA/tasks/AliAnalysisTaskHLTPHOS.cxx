// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        *
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Zhongbao Yin <zbyin@mail.ccnu.edu.cn>,                *
//*                  Kalliopi Kanaki <Kalliopi.Kanaki@ift.uib.no>          *
//*                  Svein Lindal <svein.lindal@gmail.com>                 *
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

#include "AliESDEvent.h"
#include "AliESDRun.h"
#include "AliESDInputHandler.h"
#include "AliESDCaloCluster.h"

#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliAnalysisTaskHLTPHOS.h"


ClassImp(AliAnalysisTaskHLTPHOS)


//===========================================================================================

AliAnalysisTaskHLTPHOS::AliAnalysisTaskHLTPHOS(const char *name) : AliAnalysisTaskHLTCalo(name)
  ,fHistOnlTrk2PHOS(0)
  ,fHistOfflTrk2PHOS(0)
  ,fHistOfflTrk2PHOSTrig(0)
  ,fHistOfflTrk2PHOSNoTrig(0)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  // DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container

  // DefineOutput(1, TList::Class());
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

void AliAnalysisTaskHLTPHOS::CreateSpecificStuff(TList * fOutputList){


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


void AliAnalysisTaskHLTPHOS::DoSpecificStuff(AliESDEvent * evESD, AliESDEvent * evHLTESD) {


  Double_t b = evESD->GetMagneticField();
  
  //Double_t pos[] = { 0., 0., 0.};
  //AliVertex *vtx = new AliVertex(pos, 0., 0);
    
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
      }else if(IsInPHOS(3, esdTrk, b, v)){
	Float_t phi = v.Phi(); 
        if(phi<0) phi += 2.*TMath::Pi(); 
        fHistOfflTrk2PHOSTrig->Fill(v.Eta(), phi*TMath::RadToDeg()); 
	fHistOfflTrk2PHOS->Fill(v.Eta(), phi*TMath::RadToDeg());
      }else if( IsInPHOS(4, esdTrk, b, v) ){
	Float_t phi = v.Phi();  
        if(phi<0) phi += 2.*TMath::Pi();  
        fHistOfflTrk2PHOSTrig->Fill(v.Eta(), phi*TMath::RadToDeg());
	fHistOfflTrk2PHOS->Fill(v.Eta(), phi*TMath::RadToDeg());
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


  
  //cout << "here" << endl;

  //Int_t nc = evHLTESD->GetPHOSClusters(fClustersArray);
  //Int_t nc = evHLTESD->GetEMCALClusters(fClustersArray);
  
  //for(int i = 0; i < 
  // cout << nc << " ";
  
  // for(int i = 0; i < evHLTESD->GetNumberOfCaloClusters(); i++) {
  //   AliESDCaloCluster * c = evHLTESD->GetCaloCluster(i); 
  //   cout << c->IsPHOS() << "i ";
  // }
  

  

}

Int_t AliAnalysisTaskHLTPHOS::GetClusters(AliESDEvent * event, TRefArray * clusters) {
  return event->GetPHOSClusters(clusters);
}

Bool_t AliAnalysisTaskHLTPHOS::IsThisDetector(AliESDCaloCluster * cluster) {
  return cluster->IsPHOS();
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
