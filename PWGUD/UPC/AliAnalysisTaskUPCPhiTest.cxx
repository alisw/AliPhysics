/**************************************************************************
 * Copyright(c) 1998-1999 ALICE Experiment at CERN, All rights reserved. *
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

/* AliAnaysisTaskMyTask*/
// c++ headers
#include <iostream>
#include <string.h>

// root headers
#include "TH1I.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TParticle.h"
#include "TObjString.h"
#include "TFile.h"
#include "TDatabasePDG.h"
#include "TLorentzVector.h"

// aliroot headers
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliMCEvent.h"
#include "AliAODVZERO.h"
#include "AliAODZDC.h"
#include "AliESDVZERO.h"
#include "AliESDZDC.h"
#include "AliPIDResponse.h"
#include "AliAODTrack.h"
#include "AliAODPid.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "AliMultiplicity.h"
#include "AliESDtrack.h"
#include "AliESDMuonTrack.h"
#include "AliAODMCParticle.h"
#include "AliMCParticle.h"
#include "AliCentrality.h"
#include "AliKFVertex.h"
#include "AliExternalTrackParam.h"
//#include "AliTriggerAnalysis.h"
#include "AliAODMCHeader.h"
#include "TChain.h"
#include "TH1F.h"
#include "TList.h"
#include "AliAnalysisTask.h"
#include "AliAODInputHandler.h"
  #include "AliAnalysisTaskUPCPhiTest.h"
#include "AliPIDResponse.h"
#include "TMath.h" 

class AliAnalysisTaskUPCPhiTest;    
using namespace std;             

ClassImp(AliAnalysisTaskUPCPhiTest) // classimp: necessary for root


AliAnalysisTaskUPCPhiTest::AliAnalysisTaskUPCPhiTest() : AliAnalysisTaskSE(),
fAOD(0), fOutputList(0),fOutputList2(0),fTreeP_TPC(0),fPt(0),fM(0),fPt0(0),fPt1(0),fPiM(0),fMuM(0),fElM(0),fPIDResponse(0),fKaonSigma1(0),fKaonSigma0(0),fPiSigma1(0),fElSigma0(0),fElSigma1(0),fMuSigma0(0),fMuSigma1(0),fTree_NoCut(0),fPiSigma0(0),fKaonSigmaTOF1(0),fKaonSigmaTOF0(0),fPiSigmaTOF1(0),fPiSigmaTOF0(0),fElSigmaTOF1(0),fElSigmaTOF0(0),fMuSigmaTOF1(0),fMuSigmaTOF0(0),fTPCcluster1(0),fEta1(0),fTPCcluster2(0),fEta2(0),fHistCounter(0),fDCAxy2(0),fDCAz2(0),fDCAxy1(0),fDCAz1(0),fTriggerClass(0),fdEdX0(0),fdEdX1(0),fPd0(0),fdEdXTOF0(0),fdEdXTOF1(0),fPd1(0),fPp(0),fZNAenergy(0),fZNCenergy(0),fZDCAtime(0),fZDCCtime(0),fRunNumber(0), fHistRunCounter(0),fCharge0(0),fCharge1(0),fHistCcup4Triggers(0),fHistCcup7Triggers(0),fHistCcup2Triggers(0),fHistCint1Triggers(0),fHistCint6Triggers(0),fHistC0tvxAndCint1(0) ,fHistZedTriggers(0),fHistCvlnTriggers(0),fHistMBTriggers(0),fHistCentralTriggers(0),fHistSemiCentralTriggers(0),fHistCTest58Triggers(0),fHistCTest59Triggers(0),fHistCTest60Triggers(0),fHistCTest61Triggers(0),fHistCcup8Triggers(0),fHistCcup9Triggers(0),fHistCcup291Triggers(0),fHistCcup301Triggers(0),fHistCcup311Triggers(0),fHistCcup29Triggers(0),fHistCcup30Triggers(0),fHistCcup31Triggers(0), fHistCtrueTriggers(0),fPhi1(0),fPhi2(0),fITSInHits(0),fITSOutHits(0),fDelPhi(0)
 
{
    // default constructor, don't allocate memory here!
    // this is used by root for IO purposes, it needs to remain empty
  }
//_____________________________________________________________________________
AliAnalysisTaskUPCPhiTest::AliAnalysisTaskUPCPhiTest(const char* name) : AliAnalysisTaskSE(name),
fAOD(0), fOutputList(0),fOutputList2(0),fTreeP_TPC(0),fPt(0),fM(0),fPt0(0),fPt1(0),fPiM(0),fMuM(0),fElM(0),fPIDResponse(0),fKaonSigma1(0),fKaonSigma0(0),fPiSigma1(0),fElSigma0(0),fElSigma1(0),fMuSigma0(0),fMuSigma1(0), fTree_NoCut(0),fPiSigma0(0),fKaonSigmaTOF1(0),fKaonSigmaTOF0(0),fPiSigmaTOF1(0),fPiSigmaTOF0(0),fElSigmaTOF1(0),fElSigmaTOF0(0),fMuSigmaTOF1(0),fMuSigmaTOF0(0),fTPCcluster1(0), fEta1(0),fTPCcluster2(0), fEta2(0),fHistCounter(0),fDCAxy2(0),fDCAz2(0),fDCAxy1(0),fDCAz1(0),fTriggerClass(0),fdEdX0(0),fdEdX1(0),fPd0(0),fdEdXTOF0(0),fdEdXTOF1(0),fPd1(0),fPp(0), fZNAenergy(0),fZNCenergy(0),fZDCAtime(0),fZDCCtime(0),fRunNumber(0) , fHistRunCounter(0),fCharge0(0),fCharge1(0),fHistCcup4Triggers(0),fHistCcup7Triggers(0),fHistCcup2Triggers(0),fHistCint1Triggers(0),fHistCint6Triggers(0),fHistC0tvxAndCint1(0),fHistZedTriggers(0),fHistCvlnTriggers(0),fHistMBTriggers(0),fHistCentralTriggers(0),fHistSemiCentralTriggers(0),fHistCTest58Triggers(0),fHistCTest59Triggers(0),fHistCTest60Triggers(0),fHistCTest61Triggers(0),fHistCcup8Triggers(0),fHistCcup9Triggers(0),fHistCcup291Triggers(0),fHistCcup301Triggers(0),fHistCcup311Triggers(0),fHistCcup29Triggers(0),fHistCcup30Triggers(0),fHistCcup31Triggers(0), fHistCtrueTriggers(0),fPhi1(0),fPhi2(0),fITSInHits(0),fITSOutHits(0),fDelPhi(0)

  {
    
    
    // constructor
    DefineInput(0, TChain::Class());    // define the input of the analysis: in this case we take a 'chain' of events
                                        // this chain is created by the analysis manager, so no need to worry about it, 
                                        // it does its work automatically
    DefineOutput(1, TList::Class());    // define the ouptut of the analysis: in this case it's a list of histograms 
                                        // you can add more output objects by calling DefineOutput(2, classname::Class())
                                        // if you add more output objects, make sure to call PostData for all of them, and to
    DefineOutput(2, TList::Class());                                      // make changes to your AddTask macro!
   //DefineOutput(2, TTree::Class());
   DefineOutput(3, TTree::Class());
  }
//_____________________________________________________________________________
AliAnalysisTaskUPCPhiTest::~AliAnalysisTaskUPCPhiTest()
 {
  // destructor
  if(fOutputList) 
    {
  
    delete fOutputList;     // at the end of your task, it is deleted from memory by calling this function
    }
 if(fOutputList2) 
    {
  
    delete fOutputList2;     // at the end of your task, it is deleted from memory by calling this function
    }
 }
 
//_____________________________________________________________________________
void AliAnalysisTaskUPCPhiTest::UserCreateOutputObjects()
  {
    AliAnalysisManager *man = AliAnalysisManager::GetAnalysisManager();
    if (man) 
      {
        AliInputEventHandler* inputHandler = (AliInputEventHandler*)(man->GetInputEventHandler());
        if (inputHandler)   fPIDResponse = inputHandler->GetPIDResponse();
        if (!fPIDResponse) {
        cout<<"noPIDresopne"<<endl;
        return;
     }
   
  }
    
    
    fOutputList = new TList();          
    fOutputList->SetOwner(kTRUE);       // memory stuff: the list is owner of all objects it contains and will delete them
    fOutputList2 = new TList();          
    fOutputList2->SetOwner(kTRUE);          // if requested (dont worry about this now)
  
      // example of a histogram
         
        fHistCounter = new TH1I("fHistCounter","Counter",15,0,15);  
        fOutputList->Add(fHistCounter);  
        fHistRunCounter = new TH1D("fHistRunCounter","Counter", 40000, 270000.5, 310000.5);
        fOutputList->Add(fHistRunCounter);                              
        fHistCcup4Triggers= (TH1D*)fHistRunCounter->Clone("fHistCcup4Triggers");
        fHistCcup7Triggers= (TH1D*)fHistRunCounter->Clone("fHistCcup7Triggers");
        fHistCcup2Triggers= (TH1D*)fHistRunCounter->Clone("fHistCcup2Triggers");
        fHistCint1Triggers= (TH1D*)fHistRunCounter->Clone("fHistCint1Triggers");
        fHistCint6Triggers= (TH1D*)fHistRunCounter->Clone("fHistCint6Triggers");
        fHistC0tvxAndCint1= (TH1D*)fHistRunCounter->Clone("fHistC0tvxAndCint1");
        fHistZedTriggers= (TH1D*)fHistRunCounter->Clone("fHistZedTriggers");
        fHistCvlnTriggers= (TH1D*)fHistRunCounter->Clone("fHistCvlnTriggers");
        fHistMBTriggers= (TH1D*)fHistRunCounter->Clone("fHistMBTriggers");
        fHistCentralTriggers= (TH1D*)fHistRunCounter->Clone("fHistCentralTriggers");
        fHistSemiCentralTriggers= (TH1D*)fHistRunCounter->Clone("fHistSemiCentralTriggers");
        
        fHistCTest58Triggers= (TH1D*)fHistRunCounter->Clone("fHistCTest58Triggers");
        fHistCTest59Triggers= (TH1D*)fHistRunCounter->Clone("fHistCTest59Triggers");
        fHistCTest60Triggers= (TH1D*)fHistRunCounter->Clone("fHistCTest60Triggers");
        fHistCTest61Triggers= (TH1D*)fHistRunCounter->Clone("fHistCTest61Triggers");
        
        fHistCcup8Triggers= (TH1D*)fHistRunCounter->Clone("fHistCcup8Triggers");
        fHistCcup9Triggers= (TH1D*)fHistRunCounter->Clone("fHistCcup9Triggers");
        fHistCcup291Triggers= (TH1D*)fHistRunCounter->Clone("fHistCcup291Triggers");
        fHistCcup301Triggers= (TH1D*)fHistRunCounter->Clone("fHistCcup301Triggers");
        fHistCcup311Triggers= (TH1D*)fHistRunCounter->Clone("fHistCcup311Triggers");
        fHistCcup29Triggers= (TH1D*)fHistRunCounter->Clone("fHistCcup29Triggers");
        fHistCcup30Triggers= (TH1D*)fHistRunCounter->Clone("fHistCcup30Triggers");
        fHistCcup31Triggers= (TH1D*)fHistRunCounter->Clone("fHistCcup31Triggers");
        fHistCtrueTriggers= (TH1D*)fHistRunCounter->Clone(" fHistCtrueTriggers"); 
         fOutputList->Add(fHistCcup4Triggers);
         fOutputList->Add(fHistCcup7Triggers);
         fOutputList->Add(fHistCcup2Triggers);
         fOutputList->Add(fHistCint1Triggers);
         fOutputList->Add(fHistCint6Triggers);
         fOutputList->Add(fHistC0tvxAndCint1);
         fOutputList->Add(fHistZedTriggers);
         fOutputList->Add(fHistCvlnTriggers);
         fOutputList->Add(fHistMBTriggers);
         fOutputList->Add(fHistCentralTriggers);
         fOutputList->Add(fHistSemiCentralTriggers);
        
         fOutputList->Add(fHistCTest58Triggers);
         fOutputList->Add(fHistCTest59Triggers);
         fOutputList->Add(fHistCTest60Triggers);
         fOutputList->Add(fHistCTest61Triggers);
        
         fOutputList->Add(fHistCcup8Triggers);
         fOutputList->Add(fHistCcup9Triggers);
         fOutputList->Add(fHistCcup291Triggers);
         fOutputList->Add(fHistCcup301Triggers);
         fOutputList->Add(fHistCcup311Triggers);
         fOutputList->Add(fHistCcup29Triggers);
         fOutputList->Add(fHistCcup30Triggers);
         fOutputList->Add(fHistCcup31Triggers);
         fOutputList->Add(fHistCtrueTriggers);  
  
  
  
  fHistCounter->SetTitle("Counter for Cuts");
    TString  Cuts[12] =  {"Events","CCUPTrigger","track1filterbit1<<0","track1ITS","track1Eta|0.9|","track2filterbit1<<0","track2ITS","track2Eta|0.9|","TPCcls1<60","TCPcls2<60","PairofTracks=2","OppositeCharge"};
   for (Int_t c=0;c<12;c++) fHistCounter->GetXaxis()->SetBinLabel(c+1,Cuts[c]);
    fTreeP_TPC = new TTree("scatterplot","Momentum and TPC signal");
    fOutputList2->Add(fTreeP_TPC);
    fTreeP_TPC->Branch("fPt", &fPt, "fPt/F");
    fTreeP_TPC->Branch("fM", &fM, "fM/F");
    fTreeP_TPC->Branch("fPt0", &fPt0, "fPt0/F");
    fTreeP_TPC->Branch("fPt1", &fPt1, "fPt1/F"); 
    
    
   
    fTreeP_TPC->Branch("fPiM", &fPiM, "fPiM/F");
    fTreeP_TPC->Branch("fElM", &fElM, "fElM/F");
    fTreeP_TPC->Branch("fMuM", &fMuM, "fMuM/F");
    
    fTreeP_TPC->Branch("fKaonSigma0", &fKaonSigma0, "fKaonSigma0/F");
    fTreeP_TPC->Branch("fKaonSigma1", &fKaonSigma1, "fKaonSigma1/F");
    fTreeP_TPC->Branch("fPiSigma0", &fPiSigma0, "fPiSigma0/F");
    fTreeP_TPC->Branch("fPiSigma1", &fPiSigma1, "fPiSigma1/F");
    fTreeP_TPC->Branch("fElSigma0", &fElSigma0, "fElSigma0/F");
    fTreeP_TPC->Branch("fElSigma1", &fElSigma1, "fElSigma1/F");
    fTreeP_TPC->Branch("fMuSigma0", &fMuSigma0, "fMuSigma0/F");
    fTreeP_TPC->Branch("fMuSigma1", &fMuSigma1, "fMuSigma1/F");
    
    fTreeP_TPC->Branch("fKaonSigmaTOF0", &fKaonSigmaTOF0, "fKaonSigmaTOF0/F");
    fTreeP_TPC->Branch("fKaonSigmaTOF1", &fKaonSigmaTOF1, "fKaonSigmaTOF1/F");
    fTreeP_TPC->Branch("fPiSigmaTOF0", &fPiSigmaTOF0, "fPiSigmaTOF0/F");
    fTreeP_TPC->Branch("fPiSigmaTOF1", &fPiSigmaTOF1, "fPiSigmaTOF1/F");
    fTreeP_TPC->Branch("fElSigmaTOF0", &fElSigmaTOF0, "fElSigmaTOF0/F");
    fTreeP_TPC->Branch("fElSigmaTOF1", &fElSigmaTOF1, "fElSigmaTOF1/F");
    fTreeP_TPC->Branch("fMuSigmaTOF0", &fMuSigmaTOF0, "fMuSigmaTOF0/F");
    fTreeP_TPC->Branch("fMuSigmaTOF1", &fMuSigmaTOF1, "fMuSigmaTOF1/F");
    
    fTreeP_TPC->Branch("fTPCcluster1", &fTPCcluster1, "fTPCcluster1/F");
    fTreeP_TPC->Branch("fEta1", &fEta1, "fEta1/F");
    fTreeP_TPC->Branch("fTPCcluster2", &fTPCcluster2, "fTPCcluster2/F");
    fTreeP_TPC->Branch("fEta2", &fEta2, "fEta2/F");
    fTreeP_TPC->Branch("fDCAxy1", &fDCAxy1, "fDCAxy1/F");
    fTreeP_TPC->Branch("fDCAxy2", &fDCAxy2, "fDCAxy2/F");
    fTreeP_TPC->Branch("fDCAz1", &fDCAz1, "fDCAz1/F");
    fTreeP_TPC->Branch("fDCAz2", &fDCAz2, "fDCAz2/F");
    fTreeP_TPC->Branch("fTriggerClass", &fTriggerClass,"fTriggerClass/I");
    fTreeP_TPC->Branch("fPp", &fPp, "fPp/F");
    fTreeP_TPC->Branch("fPd1", &fPd1, "fPd1/F");
    fTreeP_TPC->Branch("fPd0", &fPd0, "fPd0/F");
    fTreeP_TPC->Branch("fdEdX0", &fdEdX0, "fdEdX0/F");
    fTreeP_TPC->Branch("fdEdX1", &fdEdX1, "fdEdX1/F");
    fTreeP_TPC->Branch("fdEdXTOF0", &fdEdXTOF0, "fdEdXTOF0/F");
    fTreeP_TPC->Branch("fdEdXTOF1", &fdEdXTOF1, "fdEdXTOF1/F");  
    fTreeP_TPC->Branch("fZNAenergy", &fZNAenergy, "fZNAenergy/F"); 
    fTreeP_TPC->Branch("fZNCenergy", &fZNCenergy, "fZNCenergy/F");
    fTreeP_TPC->Branch("fZDCAtime", &fZDCAtime, "fZDCAtime/F");
    fTreeP_TPC->Branch("fZDCCtime", &fZDCCtime, "fZDCCtime/F");
    fTreeP_TPC->Branch("fRunNumber", &fRunNumber, "fRunNumber/I");
    fTreeP_TPC->Branch("fCharge0", &fCharge0, "fCharge0/I");
    fTreeP_TPC->Branch("fCharge1", &fCharge1, "fCharge1/I");
    fTreeP_TPC->Branch("fPhi1", &fPhi1, "fPhi1/F");
    fTreeP_TPC->Branch("fPhi2", &fPhi2, "fPhi2/F");
     fTreeP_TPC->Branch("fDelPhi", &fDelPhi, "fDelPhi/F");
    fTreeP_TPC->Branch("fITSInHits",  &fITSInHits, "fITSInHits/I");
    fTreeP_TPC->Branch("fITSOutHits",  &fITSOutHits, "fITSOutHits/I");
    PostData(1, fOutputList);           
   
  
    PostData(2,fOutputList2);
 //PostData(1,fTree_NoCut);
  }
//_____________________________________________________________________________
void AliAnalysisTaskUPCPhiTest::UserExec(Option_t *)
  {



 
     fAOD = dynamic_cast<AliAODEvent*>(InputEvent());    
     if(!fAOD) return;                                    
   fRunNumber = fAOD ->GetRunNumber();
    fHistCounter->Fill(0);
	 TString trigger = fAOD-> GetFiredTriggerClasses();
      //if(!trigger.Contains("CCUP")) return;
      fHistCounter->Fill(1);
  // cout<<trigger<<endl; //return;
    if(trigger.Contains("CCUP4-B")) fHistCcup4Triggers->Fill(fRunNumber); //CCUP4 triggers
  if(trigger.Contains("CCUP7-B")) fHistCcup7Triggers->Fill(fRunNumber); //CCUP7 triggers
  if(trigger.Contains("CCUP2-B")) fHistCcup2Triggers->Fill(fRunNumber); //CCUP2 triggers
  
  if(trigger.Contains("CINT1-B")) fHistCint1Triggers->Fill(fRunNumber); //CINT1 triggers
  
  if(trigger.Contains("CTEST58-B")) fHistCTest58Triggers->Fill(fRunNumber); //CTEST triggers
  if(trigger.Contains("CTEST59-B")) fHistCTest59Triggers->Fill(fRunNumber); //CTEST triggers
  if(trigger.Contains("CTEST60-B")) fHistCTest60Triggers->Fill(fRunNumber); //CTEST triggers
  if(trigger.Contains("CTEST61-B")) fHistCTest61Triggers->Fill(fRunNumber); //CTEST triggers
  
  if(trigger.Contains("CCUP8-B")) fHistCcup8Triggers->Fill(fRunNumber); //CCUP8 triggers
  if(trigger.Contains("CCUP9-B")) fHistCcup9Triggers->Fill(fRunNumber); //CCUP9 triggers
  if(trigger.Contains("CCUP29-B-NOPF")) fHistCcup291Triggers->Fill(fRunNumber); //CCUP29-nopf triggers
  if(trigger.Contains("CCUP30-B-NOPF")) fHistCcup301Triggers->Fill(fRunNumber); //CCUP30-nopf triggers
  if(trigger.Contains("CCUP31-B-NOPF")) fHistCcup311Triggers->Fill(fRunNumber); //CCUP31-nopf triggers
  
  if(trigger.Contains("CCUP29-B-SPD2")) fHistCcup29Triggers->Fill(fRunNumber); //CCUP29-spd triggers
  if(trigger.Contains("CCUP30-B-SPD2")) fHistCcup30Triggers->Fill(fRunNumber); //CCUP30-spd triggers
  if(trigger.Contains("CCUP31-B-SPD2")) fHistCcup31Triggers->Fill(fRunNumber); //CCUP31-spd triggers
  
  if(trigger.Contains("CTRUE-B")) fHistCtrueTriggers->Fill(fRunNumber); //CTRUE triggers
   
  //   if (trigger.Contains("CCUP8")) return; 
    // cout<< "trigger classes are " << trigger << endl; return;  
     
     if(trigger.Contains("CCUP4-B")) fTriggerClass =4;//CCUP4 triggers
  if(trigger.Contains("CCUP7-B")) fTriggerClass =7; //CCUP7 triggers
  if(trigger.Contains("CCUP2-B")) fTriggerClass =2; //CCUP2 triggers
  
  if(trigger.Contains("CINT1-B")) fTriggerClass =1; //CINT1 triggers
  
  if(trigger.Contains("CTEST58-B")) fTriggerClass =58; //CTEST triggers
  if(trigger.Contains("CTEST59-B")) fTriggerClass =59; //CTEST triggers
  if(trigger.Contains("CTEST60-B")) fTriggerClass =60;//CTEST triggers
  if(trigger.Contains("CTEST61-B")) fTriggerClass =61; //CTEST triggers
  
  if(trigger.Contains("CCUP8-B")) fTriggerClass =8; //CCUP8 triggers
  if(trigger.Contains("CCUP9-B")) fTriggerClass =9; //CCUP9 triggers
  if(trigger.Contains("CCUP29-B-NOPF")) fTriggerClass =25; //CCUP29-nopf triggers
  if(trigger.Contains("CCUP30-B-NOPF")) fTriggerClass =26;; //CCUP30-nopf triggers
  if(trigger.Contains("CCUP31-B-NOPF")) fTriggerClass =27; //CCUP31-nopf triggers
  
  if(trigger.Contains("CCUP29-B-SPD2")) fTriggerClass =29; //CCUP29-spd triggers
  if(trigger.Contains("CCUP30-B-SPD2")) fTriggerClass =30; //CCUP30-spd triggers
  if(trigger.Contains("CCUP31-B-SPD2"))fTriggerClass =31; //CCUP31-spd triggers
  
  if(trigger.Contains("CTRUE-B")) fTriggerClass =32; //CTRUE triggers     
         
         
  /*   if (trigger.Contains("CCUP8") && !trigger.Contains("CCUP9") ){
         fTriggerClass =1;
         
         } 
     if (trigger.Contains("CCUP9") && !trigger.Contains("CCUP8")){
        fTriggerClass =2;
        } */   
    
      AliAODZDC *fZDCdata = fAOD->GetZDCData();
      fZNAenergy = fZDCdata->GetZNATowerEnergy()[0];
      fZNCenergy = fZDCdata->GetZNCTowerEnergy()[0];
      fZDCAtime = fZDCdata->GetZNATime();
      fZDCCtime = fZDCdata->GetZNCTime();
  //ZDC cuts but not applying here, planning to apply in the final result...
  //if(trigger.Contains("CCUP4-B"))fHistZDCCuts->Fill(1);
  /*if(fZNAenergy < 8200 && fZNCenergy < 8200) fHistZDCCuts->Fill(2);
  if(fZNAenergy < 683 && fZNCenergy < 683) fHistZDCCuts->Fill(3);
  if(fZDCAtime == 0 && fZDCCtime == 0) fHistZDCCuts->Fill(4);*/           
     Int_t iTracks(fAOD->GetNumberOfTracks()); // see how many tracks there are in the event
     Int_t PairCounter = 0;
     Bool_t GoodTracks = kFALSE;
     Bool_t GoodTracks2 = kFALSE;
     AliAODTrack* savetrack1;
     AliAODTrack* savetrack2;
     for(Int_t i(0); i < iTracks; i++) 
      { 
        AliAODTrack* track = static_cast<AliAODTrack*>(fAOD->GetTrack(i));                // loop over all these tracks
        
        if(!(track->TestFilterBit(1<<0))) continue;
        fHistCounter->Fill(2); 
        if(track->HasPointOnITSLayer(0) || track->HasPointOnITSLayer(1)) {GoodTracks = kTRUE;}
        
        else  {GoodTracks = kFALSE;}
        if (GoodTracks == kFALSE) continue; 
        if (GoodTracks== kTRUE) fITSOutHits = 1;
        fHistCounter->Fill(3); 
           
        if(TMath::Abs(track->Eta())>0.9)   continue;
        fHistCounter->Fill(4);
        if(track->GetTPCNcls() < 60)continue;
        fHistCounter->Fill(5);
   //     fHistP_TPC-> Fill(track->P(),track->GetTPCsignal());
   
                       
       
        for (Int_t j=i+1 ; j<iTracks ; j++) {
            AliAODTrack* track2 = static_cast<AliAODTrack*>(fAOD->GetTrack(j));      
            if(!(track2->TestFilterBit(1<<0))) continue;
            fHistCounter->Fill(6);
            if(track2->HasPointOnITSLayer(0) || track2->HasPointOnITSLayer(1)) {GoodTracks2 = kTRUE;}
            else { GoodTracks2 = kFALSE ;}
            if (GoodTracks2 == kFALSE) continue;
            if (GoodTracks2== kTRUE) fITSInHits = 1;
           
             fHistCounter-> Fill(7);
            if(TMath::Abs(track2->Eta())>0.9)   continue;
            fHistCounter->Fill(8);
            if(track2->GetTPCNcls() < 60)continue;
            fHistCounter->Fill(9);
            PairCounter = PairCounter+1;  
            savetrack1 = track;
            savetrack2 = track2;
            if (PairCounter>1) break;     
            
         
          
          
            }
       //end of for loop j tracks*/
       //fTreeP_TPC ->Fill();
      }//end of for loop i tracks
      //Selecting tracks with only 1 pair of tracks
      if (PairCounter!=1)  {return;}
      fHistCounter->Fill(10);
      Double_t dca[2] = {0.0,0.0}, cov[3] = {0.0,0.0,0.0};
      AliAODTrack* track1_clone=(AliAODTrack*)savetrack1->Clone("track1_clone");
      AliAODVertex *fAODVertex = fAOD->GetPrimaryVertex();
     
      if(!track1_clone->PropagateToDCA(fAODVertex,fAOD->GetMagneticField(),300.,dca,cov)) { 
       dca[0] = -999 ;
       dca[1] = -999 ;
       }
      fDCAxy1= dca[0];
      fDCAz1=  dca[1];
      delete track1_clone;
      AliAODTrack* track2_clone=(AliAODTrack*)savetrack2->Clone("track2_clone");
       if(!track2_clone->PropagateToDCA(fAODVertex,fAOD->GetMagneticField(),300.,dca,cov)){ 
       dca [0] = -999 ;
       dca [1] = -999 ;
       }
      fDCAxy2= dca[0];      
      fDCAz2= dca[1];
      delete track2_clone;
      // if(TMath::Abs(dca[1]) > 2) continue;
      //Charge Cut Will be applied later here charge  is stored in result
      fCharge0 = savetrack1->Charge();
      fCharge1 = savetrack2->Charge();	
	// if (savetrack1->Charge() * savetrack2->Charge()>= 0) return;
      fHistCounter->Fill(11);   
      // PID TPC
      fPiSigma0=fPIDResponse->NumberOfSigmasTPC(savetrack1, AliPID::kPion);
      fPiSigma1=fPIDResponse->NumberOfSigmasTPC(savetrack2, AliPID::kPion);
      fKaonSigma0=fPIDResponse->NumberOfSigmasTPC(savetrack1, AliPID::kKaon);
      fKaonSigma1=fPIDResponse->NumberOfSigmasTPC(savetrack2, AliPID::kKaon);
      fElSigma0=fPIDResponse->NumberOfSigmasTPC(savetrack1, AliPID::kElectron);
      fElSigma1=fPIDResponse->NumberOfSigmasTPC(savetrack2, AliPID::kElectron);
      fMuSigma0=fPIDResponse->NumberOfSigmasTPC(savetrack1, AliPID::kMuon);
      fMuSigma1=fPIDResponse->NumberOfSigmasTPC(savetrack2, AliPID::kMuon);

      //Pid from TOF 
      fPiSigmaTOF0=fPIDResponse->NumberOfSigmasTOF(savetrack1, AliPID::kPion);
      fPiSigmaTOF1=fPIDResponse->NumberOfSigmasTOF(savetrack2, AliPID::kPion);
      fKaonSigmaTOF0=fPIDResponse->NumberOfSigmasTOF(savetrack1, AliPID::kKaon);
      fKaonSigmaTOF1=fPIDResponse->NumberOfSigmasTOF(savetrack2, AliPID::kKaon);
      fElSigmaTOF0=fPIDResponse->NumberOfSigmasTOF(savetrack1, AliPID::kElectron);
      fElSigmaTOF1=fPIDResponse->NumberOfSigmasTOF(savetrack2, AliPID::kElectron);
      fMuSigmaTOF0=fPIDResponse->NumberOfSigmasTOF(savetrack1, AliPID::kMuon);
      fMuSigmaTOF1=fPIDResponse->NumberOfSigmasTOF(savetrack2, AliPID::kMuon);
    

   //if (fKaonSigma0*fKaonSigma0+fKaonSigma1*fKaonSigma1>16) return;
   //   if (fPiSigma0*fPiSigma0+fPiSigma1*fPiSigma1<25) return;
      TLorentzVector d1;
      TLorentzVector d2; 
      TLorentzVector pid1;
      TLorentzVector pid2;
      TLorentzVector mud1;
      TLorentzVector mud2;
      TLorentzVector eld1;
      TLorentzVector eld2;
      d1.SetPtEtaPhiM(savetrack1->Pt(),savetrack1->Eta(),savetrack1->Phi(),0.493);
      d2.SetPtEtaPhiM(savetrack2->Pt(),savetrack2->Eta(),savetrack2->Phi(),0.493);
      pid1.SetPtEtaPhiM(savetrack1->Pt(),savetrack1->Eta(),savetrack1->Phi(),0.139);
      pid2.SetPtEtaPhiM(savetrack2->Pt(),savetrack2->Eta(),savetrack2->Phi(),0.139);
      mud1.SetPtEtaPhiM(savetrack1->Pt(),savetrack1->Eta(),savetrack1->Phi(),0.105);
      mud2.SetPtEtaPhiM(savetrack2->Pt(),savetrack2->Eta(),savetrack2->Phi(),0.105);
      eld1.SetPtEtaPhiM(savetrack1->Pt(),savetrack1->Eta(),savetrack1->Phi(),0.000511);
      eld2.SetPtEtaPhiM(savetrack2->Pt(),savetrack2->Eta(),savetrack2->Phi(),0.000511);

      
      fPt0 = d1.Pt();
      fPt1 = d2.Pt();
     
      fdEdX0 =savetrack1->GetTPCsignal();
      fdEdX1 =savetrack2->GetTPCsignal();
     
      fdEdXTOF0 =savetrack1->GetTOFsignal();
      fdEdXTOF1 =savetrack2->GetTOFsignal();
     
      fPd0  =  d1.P();
      fPd1  =  d2.P();  
      // cout<<"mass of track"<< track->M()<<endl; 
      TLorentzVector p = d1+d2;
      TLorentzVector p2 = pid1+pid2;
      TLorentzVector p3 = mud1+mud2;
      TLorentzVector p4 = eld1+eld2;
      fPt = p.Pt();
      fPhi1= savetrack1->Phi();
      fPhi2= savetrack2->Phi();
      if (TMath::Abs(fPhi1-fPhi2)<=TMath::Pi()){
        fDelPhi = TMath::Abs(fPhi1-fPhi2);
      }
      else fDelPhi = TMath::Abs(fPhi1-fPhi2) -2*(TMath::Pi());
      fM  =  p.M();
      fMuM  =  p3.M();
      fElM  =  p4.M();
      fPiM  =  p2.M();      
      fPp   = TMath::Sqrt(p.Pt()*p.Pt()+p.Pz()*p.Pz());
      fTPCcluster1 = savetrack1->GetTPCNcls();
      fEta1 = savetrack1->Eta();
      fTPCcluster2 = savetrack2->GetTPCNcls();
      fEta2 = savetrack2->Eta();
      fTreeP_TPC ->Fill();                                // continue until all the tracks are processed
      PostData(1, fOutputList);                           // stream the results the analysis of this event to
      PostData (2,fOutputList2);                            // the output manager which will take care of writing
                                                          // it to a file
  
      }
//_____________________________________________________________________________
void AliAnalysisTaskUPCPhiTest::Terminate(Option_t *)
  {
    // terminate
   // called at the END of the analysis (when all events are processed)
  }
  //_____________________________________________________________________________
