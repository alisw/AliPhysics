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

/** @file   AliAnalysisTaskHLTTPC.cxx  
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
#include "AliAnalysisTaskHLTTPC.h"


ClassImp(AliAnalysisTaskHLTTPC)

//======================================================================================================
 
  AliAnalysisTaskHLTTPC::AliAnalysisTaskHLTTPC(const char *name)
    :
     AliAnalysisTaskSE(name)
    ,fESDRun(0)
    ,fOutputList(0)
    ,fHistTrigger(0)
    ,fHistHLTTrigger(0)
    ,fHistOfflTrkDCA(0)
    ,fHistOfflTrkDCATrig(0)
    ,fHistOfflTrkDCANoTrig(0)
    ,fHistOnlTrkDCA(0)
    ,fHistOnlTrkDCATrig(0)
    ,fHistOfflTrkNcls(0)
    ,fHistOfflTrkNclsTrig(0)
    ,fHistOfflTrkNclsNoTrig(0)
    ,fHistOnlTrkNcls(0)
    ,fHistOnlTrkNclsTrig(0)
    ,fHistOfflTrkDCANoTrigNclsCut1(0)
    ,fHistOfflTrkDCANoTrigNclsCut2(0)
    ,fHistOfflTrkP(0)
    ,fHistOfflTrkPTrig(0)
    ,fHistOfflTrkPNoTrig(0)
    ,fHistOnlTrkP(0)
    ,fHistOfflResPtInv(0)
    ,fHistOnlResPtInv(0)
    ,fHistOffldEdx(0)
    ,fHistOnldEdx(0)
    ,fHistOffldEdxVsP(0)
    ,fHistOnldEdxVsP(0)
    ,fHistOffldZ(0) 
    ,fHistOnldZ(0) 
    ,fHistOffldX(0) 
    ,fHistOnldX(0) 
    ,fHistOfflPhi(0) 
    ,fHistOnlPhi(0) 
    ,fHistOfflTheta(0) 
    ,fHistOnlTheta(0)
    ,fHistOnlDZ(0)
    ,fHistOfflDZ(0)
    ,fHistOfflDZTrig(0)
    ,fHistOfflDZNoTrig(0)
    ,fNevt(0)
    ,fTrgClsArray(0)     
     //,fNtracksThruZ0(0),
     //,fNtracksThruZ0Trig(0)
{
  // Constructor

  // Define input and output slots here
  // Input slot #0 works with a TChain
  // DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container

  DefineOutput(1, TList::Class());
}

const Float_t AliAnalysisTaskHLTTPC::fgkEtaMin = -0.12;  
const Float_t AliAnalysisTaskHLTTPC::fgkEtaMax =  0.12;  
const Float_t AliAnalysisTaskHLTTPC::fgkPhiMin[5]   = {3.83972, 4.18879, 4.53786, 4.88692, 5.23599};  
const Float_t AliAnalysisTaskHLTTPC::fgkPhiMax[5]   = {4.18879, 4.53786, 4.88692, 5.23599, 5.58505};  
const Float_t AliAnalysisTaskHLTTPC::fgkNormX[5]    = {-0.642788, -0.34202, 0, 0.34202, 0.642788};  
const Float_t AliAnalysisTaskHLTTPC::fgkNormY[5]    = {-0.766044, -0.939693, -1, -0.939693, -0.766044};  
const Float_t AliAnalysisTaskHLTTPC::fgkInitPosX[5] = {-295.682, -157.329, 0, 157.329, 295.682};  
const Float_t AliAnalysisTaskHLTTPC::fgkInitPosY[5] = {-352.38, -432.259, -460, -432.259, -352.38};

//----------------------------------------------------------------------------------------------------

void AliAnalysisTaskHLTTPC::UserCreateOutputObjects(){
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
 
  fHistOfflTrkDCA       = new TH1F("fHistOfflTrkDCA",      "DCA to beam line (offline)",               250, 0, 250);
  fHistOfflTrkDCATrig   = new TH1F("fHistOfflTrkDCATrig",  "DCA to beam line (triggered offline)",     250, 0, 250);
  fHistOfflTrkDCANoTrig = new TH1F("fHistOfflTrkDCANoTrig","DCA to beam line (offline, not triggered)",250, 0, 250);
 
  fHistOnlTrkDCA     = new TH1F("fHistOnlTrkDCA",    "DCA to beam line (HLT)",          250, 0, 250);
  fHistOnlTrkDCATrig = new TH1F("fHistOnlTrkDCATrig","DCA to beam line (HLT triggered)",250, 0, 250); 
 
  fHistOfflTrkNcls       = new TH1F("fHistOfflTrkNcls",      "clusters per track (offline)",               200, 0, 200);
  fHistOfflTrkNclsTrig   = new TH1F("fHistOfflTrkNclsTrig",  "clusters per track (triggered offline)",     200, 0, 200); 
  fHistOfflTrkNclsNoTrig = new TH1F("fHistOfflTrkNclsNoTrig","clusters per track (offline, not triggered)",200, 0, 200);

  fHistOnlTrkNcls     = new TH1F("fHistOnlTrkNcls",    "clusters per track (HLT)",           200, 0, 200);
  fHistOnlTrkNclsTrig = new TH1F("fHistOnlTrkNclsTrig","clusters per track (HLT triggered)", 200, 0, 200); 
 
  fHistOfflTrkDCANoTrigNclsCut1 = new TH1F("fHistOfflTrkDCANoTrigNclsCut1", "DCA to beam line (offline Ncls>=60, not triggered)",250, 0, 250);
  fHistOfflTrkDCANoTrigNclsCut2 = new TH1F("fHistOfflTrkDCANoTrigNclsCut2", "DCA to beam line (offline Ncls<60, not triggered)", 250, 0, 250);
 
  fHistOfflTrkP       = new TH1F("fHistOfflTrkP",      "momentum (offline)",               100, 0., 100);
  fHistOfflTrkPTrig   = new TH1F("fHistOfflTrkPTrig",  "momentum (offline triggered)",     100, 0., 100);
  fHistOfflTrkPNoTrig = new TH1F("fHistOfflTrkPNoTrig","momentum (offline, not triggered)",100, 0., 100);
 
  fHistOnlTrkP = new TH1F("fHistOnlTrkP", "momentum (HLT)", 100, 0., 100);
 
  fHistOfflResPtInv = new TH1F("fHistOfflResPtInv","resolution on 1/pt for the case with 2 offline tracks",100, 0., 1); // cosmics
  fHistOnlResPtInv  = new TH1F("fHistOnlResPtInv", "resolution on 1/pt for the case with 2 HLT tracks",    100, 0., 1); // cosmics
 
  fHistOffldEdx = new TH1F("fHistOffldEdx", "energy loss (offline)",500, 0, 500);
  fHistOnldEdx  = new TH1F("fHistOnldEdx",  "energy loss (HLT)",    500, 0, 500);
 
  fHistOffldEdxVsP = new TH2F("fHistOffldEdxVsP","dE/dx vs. momentum (offline)",100, 0., 100., 500, 0., 500.);
  fHistOnldEdxVsP  = new TH2F("fHistOnldEdxVsP", "dE/dx vs. momentum (HLT)",    100, 0., 100., 500, 0., 500.);
 
  fHistOffldZ = new TH1F("fHistOffldZ","z resolution (offline)",100, 0, 5.);
  fHistOnldZ  = new TH1F("fHistOnldZ", "z resolution (HLT)",    100, 0.,5.);

  fHistOffldX = new TH1F("fHistOffldX","r resolution (offline)",100, 0., 5.);
  fHistOnldX  = new TH1F("fHistOnldX", "r resolution (HLT)",    100, 0., 5.);
  
  fHistOfflPhi = new TH1F("fHistOfflPhi","#phi resolution (offline)",100, 0., 10); // in mrad
  fHistOnlPhi  = new TH1F("fHistOnlPhi", "#phi resolution (HLT)",    100, 0., 10); // in mrad
 
  fHistOfflTheta = new TH1F("fHistOfflTheta","#theta resolution (offline)",100, 0., 10);
  fHistOnlTheta  = new TH1F("fHistOnlTheta", "#theta resolution (HLT)",    100, 0., 10);

  fHistOfflDZ = new TH2F("fHistOfflDZ","track D vs. Z (offline)",1000, 0., 250, 1000, 0., 250);
  fHistOnlDZ  = new TH2F("fHistOnlDZ", "track D vs. Z (HLT)",    1000, 0., 250, 1000, 0., 250);
 
  fHistOfflDZTrig   = new TH2F("fHistOfflDZTrig",  "track D vs. Z (offline triggered)",    1000, 0., 250, 1000, 0., 250);
  fHistOfflDZNoTrig = new TH2F("fHistOfflDZNoTrig","track D vs. Z (offline not triggered)",1000, 0., 250, 1000, 0., 250);


  fOutputList->Add(fHistTrigger);
  fOutputList->Add(fHistHLTTrigger);

  fOutputList->Add(fHistOfflTrkDCA);
  fOutputList->Add(fHistOfflTrkDCATrig);
  fOutputList->Add(fHistOfflTrkDCANoTrig);
  fOutputList->Add(fHistOnlTrkDCA);
  fOutputList->Add(fHistOnlTrkDCATrig); 
  fOutputList->Add(fHistOfflTrkNcls);
  fOutputList->Add(fHistOfflTrkNclsTrig);
  fOutputList->Add(fHistOfflTrkNclsNoTrig);  
  fOutputList->Add(fHistOnlTrkNcls);
  fOutputList->Add(fHistOnlTrkNclsTrig);
  fOutputList->Add(fHistOfflTrkDCANoTrigNclsCut1);
  fOutputList->Add(fHistOfflTrkDCANoTrigNclsCut2);
  fOutputList->Add(fHistOfflTrkP);
  fOutputList->Add(fHistOfflTrkPTrig);
  fOutputList->Add(fHistOfflTrkPNoTrig);
  fOutputList->Add(fHistOnlTrkP);

  fOutputList->Add(fHistOfflResPtInv); // cosmics
  fOutputList->Add(fHistOnlResPtInv);  // cosmics

  fOutputList->Add(fHistOffldEdx);
  fOutputList->Add(fHistOnldEdx);
  fOutputList->Add(fHistOffldEdxVsP);
  fOutputList->Add(fHistOnldEdxVsP);
  fOutputList->Add(fHistOffldZ);
  fOutputList->Add(fHistOnldZ);
  fOutputList->Add(fHistOffldX);
  fOutputList->Add(fHistOnldX);
  fOutputList->Add(fHistOfflPhi);  
  fOutputList->Add(fHistOnlPhi);
  fOutputList->Add(fHistOfflTheta);
  fOutputList->Add(fHistOnlTheta);
  fOutputList->Add(fHistOfflDZ);
  fOutputList->Add(fHistOnlDZ);  
  fOutputList->Add(fHistOfflDZTrig);  
  fOutputList->Add(fHistOfflDZNoTrig);
}

void AliAnalysisTaskHLTTPC::NotifyRun(){
// This will not work if the active trigger classes change from run to run.
// Then one has to know all trigger classes before processing the data.

  AliESDEvent* evESD = dynamic_cast<AliESDEvent*>(InputEvent());
  TString trgClasses = evESD->GetESDRun()->GetActiveTriggerClasses(); 
 
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
  evESD = NULL;
}

void AliAnalysisTaskHLTTPC::UserExec(Option_t *){

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
  
  //Fill CTP Trigger stuff
  //fHistTrigger->Fill(evESD->GetTriggerMask());
  
  for(Int_t i = 0; i<fTrgClsArray->GetEntries(); i++){
      if((evESD->GetFiredTriggerClasses()).Contains(((TObjString *)fTrgClsArray->At(i))->GetString()))  fHistTrigger->Fill(i);
  }

  if(evHLTESD->IsHLTTriggerFired()){
     //fHistHLTTrigger->Fill(evESD->GetTriggerMask());
     for(Int_t i = 0; i<fTrgClsArray->GetEntries(); i++){ 
         if((evESD->GetFiredTriggerClasses()).Contains(((TObjString *)fTrgClsArray->At(i))->GetString())) fHistHLTTrigger->Fill(i);
     } 
  }

  if(evHLTESD->IsHLTTriggerFired()){
     for(Int_t i=0; i<evHLTESD->GetNumberOfTracks(); i++){
         AliESDtrack * HLTesdTrk = evHLTESD->GetTrack(i);
        
      if(HLTesdTrk->GetTPCNcls()>0){
	 fHistOnlTrkNcls->Fill(HLTesdTrk->GetTPCNcls());
	 fHistOnlTrkNclsTrig->Fill(HLTesdTrk->GetTPCNcls());
      }
	
      Double_t dz[2] = {-999., -999.};  
      Double_t covar[3] = {0.};
      
      HLTesdTrk->PropagateToDCA(vtx, b, 250., dz, covar); 
      fHistOnlDZ->Fill(TMath::Abs(dz[0]), TMath::Abs(dz[1])); 
      
      if(HLTesdTrk){
	 fHistOnlTrkDCA->Fill(TMath::Abs(HLTesdTrk->GetD(0., 0., b)));
	 fHistOnlTrkP->Fill(TMath::Abs(HLTesdTrk->P()));
	 fHistOnldEdx->Fill(HLTesdTrk->GetTPCsignal());
	 fHistOnldEdxVsP->Fill(TMath::Abs(HLTesdTrk->P()), HLTesdTrk->GetTPCsignal()); 
	 fHistOnlTrkDCATrig->Fill(TMath::Abs(HLTesdTrk->GetD(0., 0., b)));
      }
     }
  } else {
    for(Int_t i=0; i<evHLTESD->GetNumberOfTracks(); i++){ 
    
        AliESDtrack * HLTesdTrk = evHLTESD->GetTrack(i); 
        if(HLTesdTrk->GetTPCNcls()>0) fHistOnlTrkNcls->Fill(HLTesdTrk->GetTPCNcls()); 
                  
        Double_t dz[2] = {-999., -999.};   
        Double_t covar[3] = {0.}; 
      
        HLTesdTrk->PropagateToDCA(vtx, b, 250., dz, covar);  
        fHistOnlDZ->Fill(TMath::Abs(dz[0]), TMath::Abs(dz[1]));  
      
        if(HLTesdTrk){ 
	   fHistOnlTrkDCA->Fill(TMath::Abs(HLTesdTrk->GetD(0., 0., b))); 
	   fHistOnlTrkP->Fill(TMath::Abs(HLTesdTrk->P())); 
	   fHistOnldEdx->Fill( HLTesdTrk->GetTPCsignal());
	   fHistOnldEdxVsP->Fill(TMath::Abs(HLTesdTrk->P()), HLTesdTrk->GetTPCsignal());
        } 
    } 
  }

  if(evHLTESD->GetNumberOfTracks()==2){ // cosmics
  
     Double_t oneOverPt1 = (evHLTESD->GetTrack(0))->OneOverPt();
     Double_t oneOverPt2 = (evHLTESD->GetTrack(1))->OneOverPt();
     //cout<<"1/pt: "<<oneOverPt1<<", "<<oneOverPt2<<endl;
     fHistOnlResPtInv->Fill(2.*TMath::Abs(oneOverPt1-oneOverPt2)/(oneOverPt1+oneOverPt2));

     Float_t dz1[2], dz2[2];
     (evHLTESD->GetTrack(0))->GetDZ(0., 0., 0., b, dz1);
     (evHLTESD->GetTrack(1))->GetDZ(0., 0., 0., b, dz2);
     fHistOnldX->Fill(TMath::Abs(TMath::Abs(dz1[0])-TMath::Abs(dz2[0])));
     fHistOnldZ->Fill(TMath::Abs(dz1[1]-dz2[1]));
    
     Float_t dPhi = TMath::Abs( (evHLTESD->GetTrack(0))->Phi() - (evHLTESD->GetTrack(1))->Phi() - TMath::Pi() );
     if(dPhi>2.*TMath::Pi()) dPhi -= 2.*TMath::Pi();
    
     fHistOnlPhi->Fill(1000.*dPhi);
     fHistOnlTheta->Fill(1000.*TMath::Abs( (evHLTESD->GetTrack(0))->Theta() + (evHLTESD->GetTrack(1))->Theta() - TMath::Pi() ) );
  }
  
  if(evESD->GetNumberOfTracks()==2){ // cosmics
    
     Double_t oneOverPt1 = (evESD->GetTrack(0))->OneOverPt(); 
     Double_t oneOverPt2 = (evESD->GetTrack(1))->OneOverPt(); 
     fHistOfflResPtInv->Fill(2.*TMath::Abs(oneOverPt1-oneOverPt2)/(oneOverPt1+oneOverPt2) );

     Float_t dz1[2], dz2[2]; 
     (evESD->GetTrack(0))->GetDZ(0., 0., 0., b, dz1); 
     (evESD->GetTrack(1))->GetDZ(0., 0., 0., b, dz2); 
     fHistOffldX->Fill(TMath::Abs(TMath::Abs(dz1[0])-TMath::Abs(dz2[0]))); 
     fHistOffldZ->Fill(TMath::Abs(dz1[1]-dz2[1]));
    
     Float_t dPhi = TMath::Abs( (evESD->GetTrack(0))->Phi() - (evESD->GetTrack(1))->Phi() - TMath::Pi() );
     if(dPhi>2.*TMath::Pi()) dPhi -= 2.*TMath::Pi();
    
     fHistOfflPhi->Fill(1000.*dPhi); 
     fHistOfflTheta->Fill(1000.*TMath::Abs( (evESD->GetTrack(0))->Theta() + (evESD->GetTrack(1))->Theta() - TMath::Pi() ) ); 
  }

  //fHistTrigger->Fill(5., 1);

  //  Printf("There are %d tracks in this event", evESD->GetNumberOfTracks());

  //if(evESD->GetNumberOfTracks()>0) fHistTrigger->Fill(4., 1);

  if(evHLTESD->IsHLTTriggerFired()){
     //fHistTrigger->Fill(2., 1);

     for(Int_t i = 0; i < evESD->GetNumberOfTracks(); i++){ 
        
	 AliESDtrack *esdTrk = evESD->GetTrack(i);      
         Double_t dz[2] = {-999., -999.};  
         Double_t covar[3] = {0};
         esdTrk->PropagateToDCA(vtx, b, 250., dz, covar);
         fHistOfflDZTrig->Fill(TMath::Abs(dz[0]), TMath::Abs(dz[1]));
      
         fHistOfflDZ->Fill(TMath::Abs(dz[0]), TMath::Abs(dz[1]));
      
         /*
         Double_t pnt[3] = {0., 0., 0.};
         Double_t norm[3] = {0., 0., 1.};
         if(esdTrk->Intersect(pnt, norm, b)){
      	   if(TMath::Sqrt(pnt[0]*pnt[0]+pnt[1]*pnt[1]) < 250) {
      	     fNtracksThruZ0++;
      	     fNtracksThruZ0Trig++;
      	     fHistTrigger->Fill(6., 1);
      	     fHistTrigger->Fill(7., 1);
      	   }
         }
         */

         fHistOfflTrkDCATrig->Fill(TMath::Abs(esdTrk->GetD(0., 0., b)));
         fHistOfflTrkDCA->Fill(TMath::Abs(esdTrk->GetD(0., 0., b))); 

         if(esdTrk->GetTPCNcls()>0){
	    fHistOfflTrkNclsTrig->Fill(esdTrk->GetTPCNcls()); 
	    fHistOfflTrkNcls->Fill(esdTrk->GetTPCNcls());
         }

         fHistOfflTrkPTrig->Fill(TMath::Abs(esdTrk->P()));
         fHistOfflTrkP->Fill(TMath::Abs(esdTrk->P()));
         fHistOffldEdx->Fill( esdTrk->GetTPCsignal());
         fHistOffldEdxVsP->Fill(TMath::Abs(esdTrk->P()), esdTrk->GetTPCsignal());
     }
  } else {

    for(Int_t i=0; i<evESD->GetNumberOfTracks(); i++){ 
      
        AliESDtrack * esdTrk = evESD->GetTrack(i); 
        Double_t dz[2] = {0};  
        Double_t covar[3] = {0};
        esdTrk->PropagateToDCA(vtx, b, 250., dz, covar); 
        fHistOfflDZNoTrig->Fill(TMath::Abs(dz[0]), TMath::Abs(dz[1])); 

        fHistOfflDZ->Fill(TMath::Abs(dz[0]), TMath::Abs(dz[1]));
      
         /*
         Double_t pnt[3] = {0., 0., 0.}; 
         Double_t norm[3] = {0., 0., 1.}; 
         if(esdTrk->Intersect(pnt, norm, b)){ 
      	   if(TMath::Sqrt(pnt[0]*pnt[0]+pnt[1]*pnt[1]) < 250) { 
      	     fNtracksThruZ0++; 
      	     fHistTrigger->Fill(6., 1);
      	   } 
         } 
         */
  
         fHistOfflTrkDCANoTrig->Fill(TMath::Abs(esdTrk->GetD(0., 0., b)));
      
         if(esdTrk->GetTPCNcls()>0) fHistOfflTrkNclsNoTrig->Fill(esdTrk->GetTPCNcls()); 
      
         if(esdTrk->GetTPCNcls()>=60) fHistOfflTrkDCANoTrigNclsCut1->Fill(TMath::Abs(esdTrk->GetD(0., 0., b)));
         else fHistOfflTrkDCANoTrigNclsCut2->Fill(TMath::Abs(esdTrk->GetD(0., 0., b)));
      
         fHistOfflTrkDCA->Fill(TMath::Abs(esdTrk->GetD(0., 0., b)));
         fHistOfflTrkNcls->Fill(esdTrk->GetTPCNcls());
      
         fHistOfflTrkPNoTrig->Fill(TMath::Abs(esdTrk->P()));
         fHistOfflTrkP->Fill(TMath::Abs(esdTrk->P()));
         fHistOffldEdx->Fill( esdTrk->GetTPCsignal());
         fHistOffldEdxVsP->Fill(TMath::Abs(esdTrk->P()), esdTrk->GetTPCsignal());
    }       
  }

  fNevt++;
  delete vtx;

  // Post output data.
  PostData(1, fOutputList);
 }

void AliAnalysisTaskHLTTPC::Terminate(Option_t *){
  /*
  Printf("Number of tracks thru CE: %d", fNtracksThruZ0);
  Printf("Number of tracks thru CE from triggered events: %d", 
	 fNtracksThruZ0Trig);
  */

  // Draw result to the screen
  // Called once at the end of the query

  //  TCanvas *c1 = new TCanvas("AliAnalysisTaskHLTTPC","Trigger",10,10,510,510);
  //fHistTrigger->DrawCopy("E");
  
}
