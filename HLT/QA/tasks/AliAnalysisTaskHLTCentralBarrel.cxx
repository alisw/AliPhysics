// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        *
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Kalliopi Kanaki <Kalliopi.Kanaki@ift.uib.no>          *
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


/** @file  AliAnalysisTaskHLTCentralBarrel.cxx  
    @author Per Ivar Lønne, Hege Erdal, Kalliopi Kanaki
    @date 
    @brief An analysis task containing
    loops over HLT and offline ESD trees for comparing
    event and track properties
*/

#include <iostream>

#include "AliAnalysisTaskHLTCentralBarrel.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDInputHandler.h"
#include "AliTracker.h"
#include "AliESDVZERO.h"
//#include "AliHLTGlobalTriggerDecision.h"
#include "AliCentrality.h"

#include "TAxis.h"
#include "TString.h"
#include "TText.h"
#include "TTimeStamp.h"
#include "TH1.h"

ClassImp(AliAnalysisTaskHLTCentralBarrel)
//_______________________________________________________________________________________________//

AliAnalysisTaskHLTCentralBarrel::AliAnalysisTaskHLTCentralBarrel()
:AliAnalysisTaskSE()
  ,fUseHLTTrigger(kFALSE)
  ,fUseCentrality(kTRUE)   
  ,fOutputList(0)
  ,fEventOFF(0)
  ,fEventHLT(0)
  ,fTrackOFF(0) 
  ,fTrackHLT(0)
  ,fTextBox(0)
{
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  // DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  
  //DefineOutput(1, TList::Class());
}

AliAnalysisTaskHLTCentralBarrel::AliAnalysisTaskHLTCentralBarrel(const char *name)
:AliAnalysisTaskSE(name)    
  ,fUseHLTTrigger(kFALSE)   
  ,fUseCentrality(kTRUE)   
  ,fOutputList(0)
  ,fEventOFF(0)
  ,fEventHLT(0)
  ,fTrackOFF(0) 
  ,fTrackHLT(0) 
  ,fTextBox(0)
{
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  // DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  
  DefineOutput(1, TList::Class());
}

AliAnalysisTaskHLTCentralBarrel::~AliAnalysisTaskHLTCentralBarrel(){
  // destructor
}

void AliAnalysisTaskHLTCentralBarrel::UserCreateOutputObjects(){
  // Create THnSparse objects
  
  OpenFile(1);  
  fOutputList = new TList();
  fOutputList->SetOwner();
  fOutputList->SetName(GetName());
  
  static const int sizeEvent = 6;
 
  int    binsEvent[sizeEvent] = {  50,  50, 250,   100,   100, 2 };
  double minEvent [sizeEvent] = {  -1,  -1, -30,     0,     0, 0 };
  double maxEvent [sizeEvent] = {   1,   1,  30, 10000, 10000, 2 };
  
  fEventHLT = CreateEventTHnSparse("fEventHLT",sizeEvent,binsEvent,minEvent,maxEvent);
  fEventOFF = CreateEventTHnSparse("fEventOFF",sizeEvent,binsEvent,minEvent,maxEvent);
  
  static const int sizeTrack = 15;
  //                                pt  TPCcl theta eta phi   DCAr  DCAz charge DCArSG DCAzSG ITScl mult vertex status  vertexZ  centrality
  Int_t    binsTrack[sizeTrack] = {1500, 200, 200, 200, 200,  800,  400,    3,  400,  400,    10,  2000,     2,            250,         90 };
  Double_t minTrack [sizeTrack] = {   0,   0,  -1,  -3,  -1,  -40, -100, -1.5, -100, -100,     0,     0,     0,            -30,          0 };
  Double_t maxTrack [sizeTrack] = { 150, 200,   4,   3,   7,   40,  100,  1.5,  100,  100,    10, 20000,     2,            -30,         90 };
   
  fTrackHLT = CreateTrackTHnSparse("fTrackHLT",sizeTrack,binsTrack,minTrack,maxTrack);
  fTrackOFF = CreateTrackTHnSparse("fTrackOFF",sizeTrack,binsTrack,minTrack,maxTrack);
     
  fTextBox = new TText();  
  fOutputList->Add(fEventOFF);
  fOutputList->Add(fEventHLT);
  fOutputList->Add(fTrackOFF);
  fOutputList->Add(fTrackHLT);
  fOutputList->Add(fTextBox);
  
  PostData(1, fOutputList);
}

void AliAnalysisTaskHLTCentralBarrel::NotifyRun(){
 
  AliESDEvent *esdOFF = dynamic_cast<AliESDEvent*>(InputEvent());
  TTimeStamp* timestamp = new TTimeStamp(esdOFF->GetTimeStamp());
  fTextBox->SetName("text");
  
  TString s = Form("Run %d, Date %d", esdOFF->GetRunNumber(), timestamp->GetDate());
  fTextBox->SetTitle(s);
}

void AliAnalysisTaskHLTCentralBarrel::UserExec(Option_t *){
// see header for documentation
    
  AliESDEvent *esdOFF = dynamic_cast<AliESDEvent*>(InputEvent());  
  if(!esdOFF){
     printf("ERROR: offline ESD is not available.\n");
     return;
  }
      
  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>(fInputHandler);
  AliESDEvent *esdHLT = NULL;   

  if(esdH) esdHLT = esdH->GetHLTEvent();  
  if(!esdHLT){
     printf("ERROR: HLT ESD is not available.\n");
     return;
  }
      
  // if(fUseHLTTrigger && !((AliHLTGlobalTriggerDecision*)esdHLT->GetHLTTriggerDecision())->Result()) return;
  
   Int_t centbin = -1;
   if(fUseCentrality){
     centbin = CalculateCentrality(esdOFF);
     //printf("Centrality bin = %d", centbin);
   } 


  
  //============================ OFFLINE =============================//

  const AliESDVertex *vertOFF = esdOFF->GetPrimaryVertexTracks();
  
  Double_t bfield = esdOFF->GetMagneticField();
  Int_t nr_tracksOFF = 0;
  
  if(esdOFF->GetEventSpecie()==16) return;

  //AliCentrality *cent = esdOFF->GetCentrality();
  //printf("QQQQQQQQQ centrality V0: %f\n",cent->GetCentralityPercentileUnchecked("V0M"));
  //printf("QQQQQQQQQ centrality V0: %f\n", cent->GetCentralityPercentile("V0M"));
  //printf("QQQQQQQQQ centrality SPD: %f\n",cent->GetCentralityPercentile("CL1"));
  
  for(Int_t i=0; i<esdOFF->GetNumberOfTracks(); i++){
      AliESDtrack *esdTrackOFF = esdOFF->GetTrack(i);
      if (!esdTrackOFF) continue;
      if(!(esdTrackOFF->GetStatus()&AliESDtrack::kTPCin)) continue;
      nr_tracksOFF++; 
  }
         
  for(Int_t i=0; i<esdOFF->GetNumberOfTracks(); i++){
      
      AliESDtrack *esdTrackOFF = esdOFF->GetTrack(i);
      if (!esdTrackOFF) continue;
      if(!(esdTrackOFF->GetStatus()&AliESDtrack::kTPCin)) continue;

      //DCA calculations(from offline)
      Double_t x[3];
      Double_t b[3]; 
      esdTrackOFF->GetXYZ(x);
     
      AliTracker::GetBxByBz(x,b);
      Bool_t isOK = esdTrackOFF->RelateToVertexTPCBxByBz(vertOFF, b, kVeryBig);
      if(!isOK) continue;
	
      const AliExternalTrackParam *track = esdTrackOFF->GetTPCInnerParam();
      if(!track) continue;
	
      Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z
      if(vertOFF->GetX()==0 && vertOFF->GetY()==0 && vertOFF->GetZ()==0 ){
	 dca[0]=-99;
	 dca[1]=-99;
      }
      //else ?????? why doesn't it work with pp????
      esdTrackOFF->GetImpactParametersTPC(dca,cov);

      Float_t DCAr =-99, DCAz = -99.;   
      Double_t trackOFF[] = {
     			        TMath::Abs(esdTrackOFF->Pt()) 
     			       ,esdTrackOFF->GetTPCNcls()      
     			       ,esdTrackOFF->Theta()	       
     			       ,esdTrackOFF->Eta()	       
     			       ,esdTrackOFF->Phi()	       
     			       ,dca[0]  		       
     			       ,dca[1]  		       
     			       ,esdTrackOFF->Charge() 
     			       ,DCAr			       
     			       ,DCAz			       
     			       ,esdTrackOFF->GetNcls(0)
			       ,nr_tracksOFF
			       ,vertOFF->GetStatus()
			       ,vertOFF->GetZ()
			       ,centbin
     			    };
      fTrackOFF->Fill(trackOFF);
    }
    
    Double_t eventOFF[] = { vertOFF->GetX(), vertOFF->GetY(), vertOFF->GetZ(), vertOFF->GetNContributors(), nr_tracksOFF, vertOFF->GetStatus()};
    fEventOFF->Fill(eventOFF);  
    
  
  
  //======================================== HLT ==========================================//

  Int_t nr_tracksHLT = 0;  
  if(esdHLT->GetEventSpecie()==16) return;
  const AliESDVertex *vertHLT = esdHLT->GetPrimaryVertexTracks();
  
  for(Int_t i=0; i<esdHLT->GetNumberOfTracks(); i++){
      AliESDtrack *esdTrackHLT = esdHLT->GetTrack(i);
      if (!esdTrackHLT) continue;
      if(!(esdTrackHLT->GetStatus()&AliESDtrack::kTPCin)) continue;
      nr_tracksHLT++; 
  }
  
  Double_t eventHLT[] = { vertHLT->GetX(), vertHLT->GetY(), vertHLT->GetZ(), vertHLT->GetNContributors(), nr_tracksHLT, vertHLT->GetStatus()};
  fEventHLT->Fill(eventHLT);  

  for(Int_t i=0; i<esdHLT->GetNumberOfTracks(); i++){
      
      AliESDtrack *esdTrackHLT = esdHLT->GetTrack(i);
      if(!esdTrackHLT) continue;
      if(!(esdTrackHLT->GetStatus()&AliESDtrack::kTPCin)) continue;
      
      //DCA calculations
      Float_t DCAr=-99;
      Float_t DCAz=-99;
      Float_t dca[2];
      if(vertHLT->GetX()==0 && vertHLT->GetY()==0 && vertHLT->GetZ() ==0 ){
	dca[0]=-99;
	dca[1]=-99;
      }
      else{
	//Calculating DCA "old" fashion
	esdTrackHLT->GetDZ(esdHLT->GetPrimaryVertex()->GetXv(), esdHLT->GetPrimaryVertex()->GetYv(), esdHLT->GetPrimaryVertex()->GetZv(), bfield, dca);
	// plotting the DCA calculated by Sergey 
	esdTrackHLT->GetImpactParametersTPC(DCAr,DCAz);
      }
      
      Int_t dummycl[10];

      Double_t trackHLT[] = {
        		       TMath::Abs(esdTrackHLT->Pt())
        		      ,esdTrackHLT->GetTPCNcls()    
        		      ,esdTrackHLT->Theta()
        		      ,esdTrackHLT->Eta()	    
        		      ,esdTrackHLT->Phi()
        		      ,dca[0]			    
        		      ,dca[1]			    
        		      ,esdTrackHLT->Charge()	    
        		      ,DCAr			    
        		      ,DCAz			    
			      //,esdTrackHLT->GetNcls(0)
			      ,esdTrackHLT->GetITSclusters(dummycl)
			      ,nr_tracksHLT
			      ,vertHLT->GetStatus()
			      ,vertHLT->GetZ()
			      ,centbin 
        		    };
      fTrackHLT->Fill(trackHLT);      
  }               
  // Post output data.
  PostData(1, fOutputList);
}

void AliAnalysisTaskHLTCentralBarrel::Terminate(Option_t *){
  // see header file of AliAnalysisTask for documentation 
}

Int_t AliAnalysisTaskHLTCentralBarrel::CalculateCentrality(AliESDEvent* esd){
//see header for documentation

  Int_t centrality = -1;

  AliESDVZERO* esdV0 = esd->GetVZEROData();
  //AliESDZDC* esdZDC = esd->GetZDCData();
  //Int_t partZDC = esdZDC->GetZDCParticipants();

  Float_t multV0 = esdV0->GetMTotV0A() + esdV0->GetMTotV0C();
    
  if (      multV0 >=    0.  && multV0 <=   124.5 ) centrality = 90;
  else if ( multV0 >   124.5 && multV0 <=   274.5 ) centrality = 80;
  else if ( multV0 >   274.5 && multV0 <=   574.5 ) centrality = 70;
  else if ( multV0 >   574.5 && multV0 <=  1224.5 ) centrality = 60;
  else if ( multV0 >  1224.5 && multV0 <=  2174.5 ) centrality = 50;
  else if ( multV0 >  2174.5 && multV0 <=  3624.5 ) centrality = 40;
  else if ( multV0 >  3624.5 && multV0 <=  5574.5 ) centrality = 30;
  else if ( multV0 >  5574.5 && multV0 <=  8274.5 ) centrality = 20;
  else if ( multV0 >  8274.5 && multV0 <= 12024.5 ) centrality = 10;
  else if ( multV0 > 12024.5 && multV0 <= 14674.5 ) centrality = 5;
  else if ( multV0 > 14674.5 && multV0 <= 19449.5 ) centrality = 0;

  return centrality;
}

THnSparseF* AliAnalysisTaskHLTCentralBarrel::CreateEventTHnSparse(const char* name, Int_t size, const Int_t* bins, Double_t* min, Double_t* max){
//see header for documentation                     
  
  THnSparseF *thn = new THnSparseF(name,"",size,bins,min,max);
  thn->GetAxis(0)->SetTitle("primary vertex x");
  thn->GetAxis(1)->SetTitle("primary vertex y");
  thn->GetAxis(2)->SetTitle("primary vertex z");
  thn->GetAxis(3)->SetTitle("number of contributors");
  thn->GetAxis(4)->SetTitle("track multiplicity");
  thn->GetAxis(5)->SetTitle("vertex status"); 
  return thn;
}

THnSparseF* AliAnalysisTaskHLTCentralBarrel::CreateTrackTHnSparse(const char* name, Int_t size, const Int_t* bins, Double_t* min, Double_t* max){
//see header for documentation                     
  
  THnSparseF *thn = new THnSparseF(name,"",size,bins,min,max);
  thn->GetAxis(0)->SetTitle("p_{T}");
  thn->GetAxis(1)->SetTitle("TPC clusters/track");
  thn->GetAxis(2)->SetTitle("#theta");
  thn->GetAxis(3)->SetTitle("#eta");
  thn->GetAxis(4)->SetTitle("#phi");
  thn->GetAxis(5)->SetTitle("DCAr");
  thn->GetAxis(6)->SetTitle("DCAz");
  thn->GetAxis(7)->SetTitle("polarity");
  thn->GetAxis(8)->SetTitle("DCArSG");
  thn->GetAxis(9)->SetTitle("DCAzSG");
  thn->GetAxis(10)->SetTitle("ITS clusters/track");  
  return thn;
}

// void AliAnalysisTaskHLTCentralBarrel::Fill(AliESDevent *esd, THnSparseF* thn){
// //see header for documentation                     
//  
//   int nTracks = 0;  
//   
//   for(int i=0; i<esdHLT->GetNumberOfTracks(); i++){
//       
//       AliESDtrack *esdTrack = esd->GetTrack(i);
//       if(!esdTrack) continue;
//       if(!(esdTrack->GetStatus()&AliESDtrack::kTPCin)) continue;
//       
//       //DCA calculations
//       Float_t DCAr=-99;
//       Float_t DCAz=-99;
//       Float_t dca[2];
//       if(vertHLT->GetX()==0 && vertHLT->GetY()==0 && vertHLT->GetZ() ==0 ){
// 	dca[0]=-99;
// 	dca[1]=-99;
//       }
//       else{
// 	//Calculating DCA "old" fashion
// 	esdTrackHLT->GetDZ(esdHLT->GetPrimaryVertex()->GetXv(), esdHLT->GetPrimaryVertex()->GetYv(), esdHLT->GetPrimaryVertex()->GetZv(), bfield, dca);
// 	// plotting the DCA calculated by Sergey 
// 	esdTrackHLT->GetImpactParametersTPC(DCAr,DCAz);
//       }
// 	
//       Double_t trackHLT[] = {
//         		        TMath::Abs(esdTrackHLT->Pt())
//         		       ,esdTrackHLT->GetTPCNcls()    
//         		       ,esdTrackHLT->Theta()	     
//         		       ,esdTrackHLT->Eta()	     
//         		       ,esdTrackHLT->Phi()	     
//         		       ,DCAr			     
//         		       ,DCAz			     
//         		       ,esdTrackHLT->Charge()	     
//         		       ,dca[0]  		     
//         		       ,dca[1]  		     
//         		       ,esdTrackHLT->GetStatus()
// 			       ,esdTrackHLT->GetNcls(0)   
//         		     };
//       fTrackHLT->Fill(trackHLT);      
//       nr_tracksHLT++;    
//   }  
//   Double_t eventHLT[] = { vertHLT->GetX(), vertHLT->GetY(), vertHLT->GetZ(), vertHLT->GetNContributors(), nr_tracksHLT, vertHLT->GetStatus()};
//   fEventHLT->Fill(eventHLT);  
// }
