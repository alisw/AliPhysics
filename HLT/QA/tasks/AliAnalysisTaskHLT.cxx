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

/** @file   AliAnalysisTaskHLT.cxx  
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
#include "AliAnalysisTaskHLT.h"


ClassImp(AliAnalysisTaskHLT)

//======================================================================================================
 
  AliAnalysisTaskHLT::AliAnalysisTaskHLT(const char *name)
    :
     AliAnalysisTaskSE(name)
    ,fESDRun(0)
    ,fOutputList(0)
    ,fHistTrigger(0)
    ,fHistHLTTrigger(0)    
    ,fMomentum_off(0)	
    ,fDCA_off(0)  	
    ,fNcluster_off(0)	
    ,fdEdx_off(0) 	
    ,fdEdx_vs_P_off(0)	
    ,fPhi_off(0)  	
    ,fTheta_off(0)	
    ,fMult_off(0) 	
    ,fXYvertex_off(0)	
    ,fXvertex_off(0)	    
    ,fYvertex_off(0)	    
    ,fZvertex_off(0)	    
    
    ,fMomentum_hlt(0)
    ,fDCA_hlt(0)  
    ,fNcluster_hlt(0)
    ,fdEdx_hlt(0)    
    ,fdEdx_vs_P_hlt(0)
    ,fPhi_hlt(0)     
    ,fTheta_hlt(0)  
    ,fMult_hlt(0)   
    ,fXYvertex_hlt(0)
    ,fXvertex_hlt(0)
    ,fYvertex_hlt(0)
    ,fZvertex_hlt(0)
    
//     ,fDCA_off_trig(0)
//     ,fNcluster_off_trig(0)
//     
//     ,fDCA_hlt_trig(0)
//     ,fNcluster_hlt_trig(0)

{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  // DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container

  DefineOutput(1, TList::Class());
}

const Float_t AliAnalysisTaskHLT::fgkEtaMin = -0.12;  
const Float_t AliAnalysisTaskHLT::fgkEtaMax =  0.12;  
const Float_t AliAnalysisTaskHLT::fgkPhiMin[5]   = {3.83972, 4.18879, 4.53786, 4.88692, 5.23599};  
const Float_t AliAnalysisTaskHLT::fgkPhiMax[5]   = {4.18879, 4.53786, 4.88692, 5.23599, 5.58505};  
const Float_t AliAnalysisTaskHLT::fgkNormX[5]    = {-0.642788, -0.34202, 0, 0.34202, 0.642788};  
const Float_t AliAnalysisTaskHLT::fgkNormY[5]    = {-0.766044, -0.939693, -1, -0.939693, -0.766044};  
const Float_t AliAnalysisTaskHLT::fgkInitPosX[5] = {-295.682, -157.329, 0, 157.329, 295.682};  
const Float_t AliAnalysisTaskHLT::fgkInitPosY[5] = {-352.38, -432.259, -460, -432.259, -352.38};

//----------------------------------------------------------------------------------------------------

void AliAnalysisTaskHLT::UserCreateOutputObjects(){
// Create histograms

  OpenFile(1);

  fOutputList = new TList();
  fOutputList->SetName(GetName());

  /*
  //0 mistriggered, 1 Good triggered, 2, triggered, 3 fake trigger, 
  //4 events with offline track, 5 total events processed,
  //6 offline track thru CE, 7 online track to CE
  fHistTrigger = new TH1F("fHistTrigger", "Trigger Status", 8, -0.5, 7.5);
  fHistTrigger->GetXaxis()->SetTitle("");
  fHistTrigger->GetYaxis()->SetTitle("Events");
  fHistTrigger->SetMarkerStyle(kFullCircle);
  fHistTrigger->SetStats(0);
  fHistTrigger->SetFillColor(2);
  //fHistTrigger->SetDrawOption("B TEXT60");

  //Set bin labels
  (fHistTrigger->GetXaxis())->SetBinLabel(1,"missed");
  (fHistTrigger->GetXaxis())->SetBinLabel(2,"triggerWofflTrk");
  (fHistTrigger->GetXaxis())->SetBinLabel(3,"triggered");
  (fHistTrigger->GetXaxis())->SetBinLabel(4,"triggerWOofflTrk");
  (fHistTrigger->GetXaxis())->SetBinLabel(5,"NevWofflTrk");
  (fHistTrigger->GetXaxis())->SetBinLabel(6,"Nevt");
  (fHistTrigger->GetXaxis())->SetBinLabel(7,"offlTrkThruCE");
  (fHistTrigger->GetXaxis())->SetBinLabel(8,"onlTrkThruCE"); 
  */

  fHistTrigger = new TH1F("fHistTrigger", "CTP trigger counter", 64, 0, 64);
  fHistTrigger->GetXaxis()->SetTitle("");  
  fHistTrigger->GetYaxis()->SetTitle("#Events"); 

  fHistHLTTrigger = new TH1F("fHistHLTTrigger", "HLT CTP trigger counter", 64, 0, 64); 
  fHistHLTTrigger->GetXaxis()->SetTitle("");
  fHistHLTTrigger->GetYaxis()->SetTitle("#Events");  

  
  
  fMomentum_off = new TH1F("fMomentum_off", "momentum (offline)",300, 0., 100);
  fMomentum_hlt = new TH1F("fMomentum_hlt", "momentum (HLT)",    300, 0., 100);
 
  fDCA_off = new TH1F("fDCA_off","DCA to beam line (offline)",250, 0, 250);
  fDCA_hlt = new TH1F("fDCA_hlt","DCA to beam line (HLT)",    250, 0, 250);
 
  fNcluster_off = new TH1F("fNcluster_off","clusters per track (offline)",200, 0, 200);
  fNcluster_hlt = new TH1F("fNcluster_hlt","clusters per track (HLT)",    200, 0, 200);
 
  fdEdx_off = new TH1F("fdEdx_off","energy loss (offline)",500, 0, 500);
  fdEdx_hlt = new TH1F("fdEdx_hlt","energy loss (HLT) - not filled yet",    500, 0, 500);
 
  fdEdx_vs_P_off = new TH2F("fdEdx_vs_P_off","dE/dx vs. momentum (offline)",100, 0., 100., 500, 0., 500.);
  fdEdx_vs_P_hlt = new TH2F("fdEdx_vs_P_hlt","dE/dx vs. momentum (HLT) - not filled yet",    100, 0., 100., 500, 0., 500.);

  fPhi_off = new TH1F("fPhi_off","azimuthal angle distribution",360,0,360);
  fPhi_hlt = new TH1F("fPhi_hlt","azimuthal angle distribution",360,0,360);
  
  fTheta_off = new TH1F("fTheta_off","polar angle distribution",360,-180,180);
  fTheta_hlt = new TH1F("fTheta_hlt","polar angle distribution",360,-180,180);
  
  fMult_off = new TH1F("fMult_off","track multiplicity (offline)",100,0,100);
  fMult_hlt = new TH1F("fMult_hlt","track multiplicity (HLT)",    100,0,100);
  
  fXYvertex_off = new TH2F("fXYvertex_off","XY primary vertex (offline)",60,-15,15,80,-20,20);
  fXYvertex_hlt = new TH2F("fXYvertex_hlt","XY primary vertex (HLT)",    60,-15,15,80,-20,20);
  
  fXvertex_off = new TH1F("fXvertex_off","X primary vertex (offline)",80,-20,20);
  fXvertex_hlt = new TH1F("fXvertex_hlt","X primary vertex (HLT)",    80,-20,20);
 
  fYvertex_off = new TH1F("fYvertex_off","Y primary vertex (offline)",80,-20,20);
  fYvertex_hlt = new TH1F("fYvertex_hlt","Y primary vertex (HLT)",    80,-20,20);
 
  fZvertex_off = new TH1F("fZvertex_off","Z primary vertex (offline)",80,-20,20);
  fZvertex_hlt = new TH1F("fZvertex_hlt","Z primary vertex (HLT)",    80,-20,20);
 

//---------------------- add histograms to the output TList ------------------//

  fOutputList->Add(fHistTrigger);
  fOutputList->Add(fHistHLTTrigger);

  fOutputList->Add(fMomentum_off);
  fOutputList->Add(fDCA_off);	  
  fOutputList->Add(fNcluster_off); 
  fOutputList->Add(fdEdx_off);	  
  fOutputList->Add(fdEdx_vs_P_off);
  fOutputList->Add(fPhi_off);	  
  fOutputList->Add(fTheta_off);    
  fOutputList->Add(fMult_off);	  
  fOutputList->Add(fXYvertex_off); 
  fOutputList->Add(fXvertex_off);  
  fOutputList->Add(fYvertex_off);  
  fOutputList->Add(fZvertex_off);  
  
  fOutputList->Add(fMomentum_hlt); 
  fOutputList->Add(fDCA_hlt);	  
  fOutputList->Add(fNcluster_hlt); 
  fOutputList->Add(fdEdx_hlt);	  
  fOutputList->Add(fdEdx_vs_P_hlt);
  fOutputList->Add(fPhi_hlt);	  
  fOutputList->Add(fTheta_hlt);    
  fOutputList->Add(fMult_hlt);	  
  fOutputList->Add(fXYvertex_hlt); 
  fOutputList->Add(fXvertex_hlt);  
  fOutputList->Add(fYvertex_hlt);  
  fOutputList->Add(fZvertex_hlt);    
}

void AliAnalysisTaskHLT::NotifyRun(){
// This will not work if the active trigger classes change from run to run.
// Then one has to know all trigger classes before processing the data.

  AliESDEvent *esdOFF = dynamic_cast<AliESDEvent*>(InputEvent());
  TString trgClasses = esdOFF->GetESDRun()->GetActiveTriggerClasses(); 
 
  /*
  TObjArray * trgClsArray = trgClasses.Tokenize(" ");
  cout<<trgClsArray->GetEntries()<<endl;

  if(!fTrgClsArray){
    fTrgClsArray = trgClsArray;
    for(Int_t i = 0; i < fTrgClsArray->GetEntries(); i++){  
      TString str = ((TObjString *)fTrgClsArray->At(i))->GetString();  
      (fHistTrigger->GetXaxis())->SetBinLabel(i+1, str.Data());  
      (fHistHLTTrigger->GetXaxis())->SetBinLabel(i+1, str.Data());  
    }  
  }else{
    for(Int_t i = 0; i < trgClsArray->GetEntries(); i++){
      
    }
  }
  */

  fTrgClsArray = trgClasses.Tokenize(" ");
  //cout<<fTrgClsArray->GetEntries()<<endl; 
    
  for(Int_t i=0; i<fTrgClsArray->GetEntries(); i++){ 
      TString str = ((TObjString *)fTrgClsArray->At(i))->GetString(); 
      (fHistTrigger->GetXaxis())->SetBinLabel(i+1, str.Data()); 
      (fHistHLTTrigger->GetXaxis())->SetBinLabel(i+1, str.Data()); 
  }   
  esdOFF = NULL;
}

void AliAnalysisTaskHLT::UserExec(Option_t *){

  AliESDEvent *esdOFF = dynamic_cast<AliESDEvent*>(InputEvent());
  
  if(!esdOFF){
      Printf("ERROR: fESD not available");
      return;
  }
  
  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>(fInputHandler);
  AliESDEvent *esdHLT = NULL;   
  if(esdH) esdHLT = esdH->GetHLTEvent();
    
  if(!esdHLT){
      Printf("ERROR: HLTesd not available");
      return;
  }

  
  //Fill CTP Trigger stuff
  //fHistTrigger->Fill(esdOFF->GetTriggerMask());
  
  for(Int_t i=0; i<fTrgClsArray->GetEntries(); i++){
      if((esdOFF->GetFiredTriggerClasses()).Contains(((TObjString *)fTrgClsArray->At(i))->GetString()))  fHistTrigger->Fill(i);
  }




  Double_t bfield = esdOFF->GetMagneticField();
 
  //---------------- HLT ESD tree -----------------------//
      
  fMult_hlt->Fill( esdHLT->GetNumberOfTracks() );
	 
  Double_t vertexHLT[3];
  vertexHLT[0] = esdHLT->GetPrimaryVertex()->GetXv();
  vertexHLT[1] = esdHLT->GetPrimaryVertex()->GetYv();
  vertexHLT[2] = esdHLT->GetPrimaryVertex()->GetZv();
  AliVertex *primVertex_hlt = new AliVertex(vertexHLT, 0., 0);
  fXYvertex_hlt->Fill( esdHLT->GetPrimaryVertex()->GetXv(), esdHLT->GetPrimaryVertex()->GetYv() );
  fXvertex_hlt->Fill( esdHLT->GetPrimaryVertex()->GetXv() );
  fYvertex_hlt->Fill( esdHLT->GetPrimaryVertex()->GetYv() );
  fZvertex_hlt->Fill( esdHLT->GetPrimaryVertex()->GetZv() );
  
  for(Int_t i=0; i<esdHLT->GetNumberOfTracks(); i++){ 
  
      AliESDtrack *esdtrackHLT = esdHLT->GetTrack(i); 
      if(esdtrackHLT){ 
               
	 fNcluster_hlt->Fill(esdtrackHLT->GetTPCNcls()); 
  		
	 //Double_t dz[2]    = {-999., -999.};   
	 //Double_t covar[3] = {0.,0.,0.};      	 
	 //esdtrackHLT->PropagateToDCA(primVertex_hlt, bfield, 250., dz, covar);  
	 //fHistOnlDZ->Fill(TMath::Abs(dz[0]), TMath::Abs(dz[1])); // z resolution 
         
         fDCA_hlt->Fill( TMath::Abs(esdtrackHLT->GetD(esdHLT->GetPrimaryVertex()->GetXv(), esdHLT->GetPrimaryVertex()->GetYv(), bfield)) ); 
         fMomentum_hlt->Fill( TMath::Abs(esdtrackHLT->P()) ); 
         fdEdx_hlt->Fill( esdtrackHLT->GetTPCsignal() );
         fdEdx_vs_P_hlt->Fill( TMath::Abs(esdtrackHLT->P()), esdtrackHLT->GetTPCsignal() );         
         fPhi_hlt->Fill(esdtrackHLT->Phi()*TMath::RadToDeg());
         fTheta_hlt->Fill(esdtrackHLT->Theta()*TMath::RadToDeg());
	 
	 if(esdHLT->IsHLTTriggerFired()){
	    
	 
	 
	 } // end if for triggered hlt events
      } // end if esdtrack is kTRUE
  } // end of loop over hlt tracks
  
  
  
  
  //----------------- OFFLINE ESD tree ----------------//
  
  fMult_off->Fill( esdOFF->GetNumberOfTracks() );

  Double_t vertexOFF[3];
  vertexOFF[0] = esdOFF->GetPrimaryVertex()->GetXv();
  vertexOFF[1] = esdOFF->GetPrimaryVertex()->GetYv();
  vertexOFF[2] = esdOFF->GetPrimaryVertex()->GetZv();
  AliVertex *primVertex_off = new AliVertex(vertexOFF, 0., 0);
  fXYvertex_off->Fill( esdOFF->GetPrimaryVertex()->GetXv(), esdOFF->GetPrimaryVertex()->GetYv() );
  fXvertex_off->Fill( esdOFF->GetPrimaryVertex()->GetXv() );
  fYvertex_off->Fill( esdOFF->GetPrimaryVertex()->GetYv() );
  fZvertex_off->Fill( esdOFF->GetPrimaryVertex()->GetZv() );
  
  
  for(Int_t i=0; i<esdOFF->GetNumberOfTracks(); i++){ 
     
      AliESDtrack *esdtrackOFF = esdOFF->GetTrack(i); 
      if(esdtrackOFF){ 
         
	 fNcluster_off->Fill(esdtrackOFF->GetTPCNcls()); 
  		
	 //Double_t dz[2]    = {-999., -999.};   
	 //Double_t covar[3] = {0.,0.,0.};      	 
	 //esdtrackOFF->PropagateToDCA(primVertex_off, bfield, 250., dz, covar);  
	 //fHistOfflDZ->Fill(TMath::Abs(dz[0]), TMath::Abs(dz[1])); // z resolution 
         
         fDCA_off->Fill( TMath::Abs(esdtrackOFF->GetD(esdOFF->GetPrimaryVertex()->GetXv(), esdOFF->GetPrimaryVertex()->GetYv(), bfield)) ); 
         fMomentum_off->Fill( TMath::Abs(esdtrackOFF->P()) ); 
         fdEdx_off->Fill( esdtrackOFF->GetTPCsignal() );
         fdEdx_vs_P_off->Fill( TMath::Abs(esdtrackOFF->P()), esdtrackOFF->GetTPCsignal() );         
         fPhi_off->Fill(esdtrackOFF->Phi()*TMath::RadToDeg());
         fTheta_off->Fill(esdtrackOFF->Theta()*TMath::RadToDeg());
	 
	 if(esdHLT->IsHLTTriggerFired()){
	 
	 
	 
	 } // end if for triggered hlt events
      } // end if esdtrack is kTRUE    
  } // end of loop over hlt tracks



//   if(esdHLT->IsHLTTriggerFired()){
// 
//      for(Int_t i=0; i<esdOFF->GetNumberOfTracks(); i++){ 
//         
// 	 AliESDtrack *esdTrk = esdOFF->GetTrack(i);      
//          Double_t dz[2] = {-999., -999.};  
//          Double_t covar[3] = {0};
//          esdTrk->PropagateToDCA(vtx, bfield, 250., dz, covar);
//          fHistOfflDZTrig->Fill(TMath::Abs(dz[0]), TMath::Abs(dz[1]));
//       
//          fHistOfflDZ->Fill(TMath::Abs(dz[0]), TMath::Abs(dz[1]));
//       
//          /*
//          Double_t pnt[3] = {0., 0., 0.};
//          Double_t norm[3] = {0., 0., 1.};
//          if(esdTrk->Intersect(pnt, norm, bfield)){
//       	   if(TMath::Sqrt(pnt[0]*pnt[0]+pnt[1]*pnt[1]) < 250) {
//       	     fNtracksThruZ0++;
//       	     fNtracksThruZ0Trig++;
//       	     fHistTrigger->Fill(6., 1);
//       	     fHistTrigger->Fill(7., 1);
//       	   }
//          }
//          */
// 
//          fHistOfflTrkDCATrig->Fill(TMath::Abs(esdTrk->GetD(0., 0., bfield)));
//          fDCA_off->Fill(TMath::Abs(esdTrk->GetD(0., 0., bfield))); 
// 
//          if(esdTrk->GetTPCNcls()>0){
// 	    fHistOfflTrkNclsTrig->Fill(esdTrk->GetTPCNcls()); 
// 	    fHistOfflTrkNcls->Fill(esdTrk->GetTPCNcls());
//          }
// 
//          fHistOfflTrkPTrig->Fill(TMath::Abs(esdTrk->P()));
//          fHistOfflTrkP->Fill(TMath::Abs(esdTrk->P()));
//          fHistOffldEdx->Fill( esdTrk->GetTPCsignal());
//          fHistOffldEdxVsP->Fill(TMath::Abs(esdTrk->P()), esdTrk->GetTPCsignal());
//      }

  fNevt++;
  delete primVertex_off;
  delete primVertex_hlt;

  // Post output data.
  PostData(1, fOutputList);
}

void AliAnalysisTaskHLT::Terminate(Option_t *){
  /*
  Printf("Number of tracks thru CE: %d", fNtracksThruZ0);
  Printf("Number of tracks thru CE from triggered events: %d", 
	 fNtracksThruZ0Trig);
  */

  // Draw result to the screen
  // Called once at the end of the query

  //  TCanvas *c1 = new TCanvas("AliAnalysisTaskHLT","Trigger",10,10,510,510);
  //fHistTrigger->DrawCopy("E");
  
}
