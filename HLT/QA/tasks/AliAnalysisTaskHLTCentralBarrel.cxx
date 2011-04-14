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
  ,fCentrality()
  ,fBeamType()
  ,fOutputList(0)
  ,fEventOFF(0)
  ,fEventHLT(0)
  ,fTrackOFF(0) 
  ,fTrackHLT(0)
  ,fOptions()
  ,fTextBox(0)
  ,fSwitch(kTRUE)
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
  ,fCentrality()   
  ,fBeamType()
  ,fOutputList(0)
  ,fEventOFF(0)
  ,fEventHLT(0)
  ,fTrackOFF(0) 
  ,fTrackHLT(0) 
  ,fOptions()
  ,fTextBox(0)
  ,fSwitch(kTRUE)
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
  
  if(fBeamType.Contains("Pb")){ // in case of a Pb+Pb run the V0 centrality is added to the THnSparse  
     static const int sizeEvent = 10;
     // 			      0    1     2    3     4     5      6       7       8               9
     // 			      x    y     z   spdx spdy  spdz  #contr   mult  vertexStatus  V0centrality
     int    binsEvent[sizeEvent] = { 100, 100,  60, 100,  100,   60,    200,	200,	2,	      100 }; // binning
     double minEvent [sizeEvent] = {  -1,  -1, -20,  -1,   -1,  -20,      0,	  3,	0,		0 }; // min x
     double maxEvent [sizeEvent] = {   1,   1,  20,   1,    1,   20,   2000,    2000,	2,	      100 }; // max x	 
     fEventHLT = CreateEventTHnSparse("fEventHLT",sizeEvent,binsEvent,minEvent,maxEvent);
     fEventOFF = CreateEventTHnSparse("fEventOFF",sizeEvent,binsEvent,minEvent,maxEvent);
  }
  else {
     static const int sizeEvent = 9;
     // 			      0    1     2    3     4     5      6       7       8       
     // 			      x    y     z   spdx spdy  spdz  #contr   mult  vertexStatus
     int    binsEvent[sizeEvent] = { 100, 100,  60, 100,  100,   60,    200,	200,	2 }; // binning
     double minEvent [sizeEvent] = {  -1,  -1, -20,  -1,   -1,  -20,      0,	  3,	0 }; // min x
     double maxEvent [sizeEvent] = {   1,   1,  20,   1,    1,   20,   2000,    2000,	2 }; // max x	
     fEventHLT = CreateEventTHnSparse("fEventHLT",sizeEvent,binsEvent,minEvent,maxEvent);
     fEventOFF = CreateEventTHnSparse("fEventOFF",sizeEvent,binsEvent,minEvent,maxEvent);  
  }
  
  if(fBeamType.Contains("Pb")){ // in case of a Pb+Pb run the V0 centrality is added to the THnSparse
     static const int sizeTrack = 13;
     //                                 0    1     2    3   4      5    6     7       8      9    10              11          12            
     // 			       pt  TPCcl theta eta phi   DCAr  DCAz charge  ITScl mult vertex status  vertexZ	V0centrality
     Int_t    binsTrack[sizeTrack] = { 400, 200, 200, 200, 200,  100,  100,    3,   10,  1000,     2,		 60,	 100 }; // binning
     Double_t minTrack [sizeTrack] = {   0,   0,  -1,  -2,  -1,  -10,  -10, -1.5,    0,     3,     0,		-20,	   0 }; // min x
     Double_t maxTrack [sizeTrack] = {  10, 200,   4,	2,   7,   10,  -10,  1.5,   10, 10000,     2,		 20,	 100 }; // max x       
     fTrackHLT = CreateTrackTHnSparse("fTrackHLT",sizeTrack,binsTrack,minTrack,maxTrack);
     fTrackOFF = CreateTrackTHnSparse("fTrackOFF",sizeTrack,binsTrack,minTrack,maxTrack);
  }
  else {    
     static const int sizeTrack = 12;
     //                                 0    1     2    3   4       5     6      7       8     9        10             11       
     // 			       pt  TPCcl theta eta  phi   DCAr  DCAz  charge   ITScl  mult  vertex status   vertexZ  
     Int_t    binsTrack[sizeTrack] = {400, 200, 200,  200,  200,  100,  100,      3,    10,  1000,     2,  	    60 }; // binning
     Double_t minTrack [sizeTrack] = {  0,   0,  -1,   -2,   -1,  -10,  -10,   -1.5,     0,     3,     0,  	   -20 }; // min x
     Double_t maxTrack [sizeTrack] = { 10 , 200,   4,    2,    7,   10,   10,    1.5,    10, 10000,     2,  	    20 }; // max x  	  
     fTrackHLT = CreateTrackTHnSparse("fTrackHLT",sizeTrack,binsTrack,minTrack,maxTrack);
     fTrackOFF = CreateTrackTHnSparse("fTrackOFF",sizeTrack,binsTrack,minTrack,maxTrack);     
  }
  
  fTextBox = new TText();  
  fOutputList->Add(fEventOFF);
  fOutputList->Add(fEventHLT);
  fOutputList->Add(fTrackOFF);
  fOutputList->Add(fTrackHLT);  
  fOutputList->Add(fTextBox);
  
  PostData(1, fOutputList);
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
  
  if(fSwitch==kTRUE){  
     TTimeStamp *timestamp = new TTimeStamp(esdHLT->GetTimeStamp());
     fTextBox->SetName("text");
     TString s = Form("Run %d, Date %d", esdHLT->GetRunNumber(), timestamp->GetDate());
     printf("You are analyzing run %d from date %d\n\n", esdHLT->GetRunNumber(), timestamp->GetDate());
     fTextBox->SetTitle(s);
     fSwitch=kFALSE;
  }
    
  // if(fUseHLTTrigger && !((AliHLTGlobalTriggerDecision*)esdHLT->GetHLTTriggerDecision())->Result()) return;  
  
  //============================ OFFLINE =============================//

  const AliESDVertex *vertOFF = esdOFF->GetPrimaryVertexTracks();  
  Double_t bfield = esdOFF->GetMagneticField();

  Int_t nr_tracksOFF = 0;  
  if(esdOFF->GetEventSpecie()==16) return; // skip calibration events

  if(fBeamType.Contains("Pb")){
     fCentrality = esdOFF->GetCentrality(); 
     // this information is only available from the offline ESD for 2010 PbPb data, the V0 info was not stored in the HLTesd (23.03.11, Kelly)
     if(!fCentrality){
        printf("Centrality pointer is empty\n");
	return;
     }
  }
               
  for(Int_t i=0; i<esdOFF->GetNumberOfTracks(); i++){     
      AliESDtrack *esdTrackOFF = esdOFF->GetTrack(i);
      if (!esdTrackOFF) continue;
      if(!(esdTrackOFF->GetStatus()&AliESDtrack::kTPCin)) continue;
      nr_tracksOFF++;

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
      esdTrackOFF->GetImpactParametersTPC(dca,cov);

      if(fBeamType.Contains("Pb")){       
         Double_t trackOFF[] = {
        			   TMath::Abs(esdTrackOFF->Pt()) 
        			  ,esdTrackOFF->GetTPCNcls()	  
        			  ,esdTrackOFF->Theta() 	  
        			  ,esdTrackOFF->Eta()		  
        			  ,esdTrackOFF->Phi()		  
        			  ,dca[0]			  
        			  ,dca[1]			  
        			  ,esdTrackOFF->Charge() 
        			  ,esdTrackOFF->GetNcls(0)
        			  ,nr_tracksOFF
        			  ,vertOFF->GetStatus()
        			  ,vertOFF->GetZ()
        			  ,fCentrality->GetCentralityPercentile("V0M")
        		       }; 
        if(fOptions.Contains("track-off")) fTrackOFF->Fill(trackOFF);
      } else {
        Double_t trackOFF[] = {
        			   TMath::Abs(esdTrackOFF->Pt()) 
        			  ,esdTrackOFF->GetTPCNcls()	  
        			  ,esdTrackOFF->Theta() 	  
        			  ,esdTrackOFF->Eta()		  
        			  ,esdTrackOFF->Phi()		  
        			  ,dca[0]			  
        			  ,dca[1]			  
        			  ,esdTrackOFF->Charge() 
        			  ,esdTrackOFF->GetNcls(0)
        			  ,nr_tracksOFF
        			  ,vertOFF->GetStatus()
        			  ,vertOFF->GetZ()        			  
        		      };
        if(fOptions.Contains("track-off")) fTrackOFF->Fill(trackOFF);
      }
  } // end of track loop
    
  if(fBeamType.Contains("Pb")){
     Double_t eventOFF[] = { vertOFF->GetX(), vertOFF->GetY(), vertOFF->GetZ(), 
                             esdOFF->GetPrimaryVertexSPD()->GetX(), esdOFF->GetPrimaryVertexSPD()->GetY(), esdOFF->GetPrimaryVertexSPD()->GetZ(),  
                             vertOFF->GetNContributors(), nr_tracksOFF, vertOFF->GetStatus(),fCentrality->GetCentralityPercentile("V0M")};
     if(fOptions.Contains("event-off")) fEventOFF->Fill(eventOFF);  
  }
  else {
     Double_t eventOFF[] = { vertOFF->GetX(), vertOFF->GetY(), vertOFF->GetZ(), 
                             esdOFF->GetPrimaryVertexSPD()->GetX(), esdOFF->GetPrimaryVertexSPD()->GetY(), esdOFF->GetPrimaryVertexSPD()->GetZ(),      
                             vertOFF->GetNContributors(), nr_tracksOFF, vertOFF->GetStatus()};			    
     if(fOptions.Contains("event-off")) fEventOFF->Fill(eventOFF);
  }
  // Inspite of the different options to fill event or track properties, all the loops over tracks are being executed. 
  // The options influence only whether the respective THnSparse is filled or not.
  // Can definitely be improved to save processing time in unnecessary loops that won't fill anything at the end. 
  // One thing to consider is the track multiplicity which is calculated as the number of tracks that have the kTPCin flag.
  // This is calculated inside the track loop.  
  
  //======================================== HLT ==========================================//

  Int_t nr_tracksHLT = 0; 
  // the following line does not have an effect yet, since HLT does not fill the relevant information (14.04.11, Kelly)
  // It does not harm to keep it in, in case this changes in the future 
  if(esdHLT->GetEventSpecie()==16) return; // skip calibration events

  const AliESDVertex *vertHLT = esdHLT->GetPrimaryVertexTracks();
  for(Int_t i=0; i<esdHLT->GetNumberOfTracks(); i++){
      
      AliESDtrack *esdTrackHLT = esdHLT->GetTrack(i);
      if(!esdTrackHLT) continue;
      if(!(esdTrackHLT->GetStatus()&AliESDtrack::kTPCin)) continue;
      nr_tracksHLT++; 
      
      Float_t dca[2];
      if(vertHLT->GetX()==0 && vertHLT->GetY()==0 && vertHLT->GetZ() ==0 ){
	dca[0]=-99;
	dca[1]=-99;
      }
      else{
	//Calculating DCA "old" fashion
	esdTrackHLT->GetDZ(esdHLT->GetPrimaryVertex()->GetXv(), esdHLT->GetPrimaryVertex()->GetYv(), esdHLT->GetPrimaryVertex()->GetZv(), bfield, dca);
	// plotting the DCA calculated by Sergey 
	//esdTrackHLT->GetImpactParametersTPC(DCAr,DCAz);
      }
      
      if(fBeamType.Contains("Pb")){
         Double_t trackHLT[] = {
        		         TMath::Abs(esdTrackHLT->Pt())
        		        ,esdTrackHLT->GetTPCNcls()    
        		        ,esdTrackHLT->Theta()
        		        ,esdTrackHLT->Eta()	      
        		        ,esdTrackHLT->Phi()
        		        ,dca[0] 		      
        		        ,dca[1] 		      
        		        ,esdTrackHLT->Charge()        
        		        ,esdTrackHLT->GetNcls(0)
        		        ,nr_tracksHLT
        		        ,vertHLT->GetStatus()
        		        ,vertHLT->GetZ()
        		        ,fCentrality->GetCentralityPercentile("V0M")
        		       };
        if(fOptions.Contains("track-hlt")) fTrackHLT->Fill(trackHLT);   
      } else {
         Double_t trackHLT[] = {
             		   	 TMath::Abs(esdTrackHLT->Pt())
             		   	,esdTrackHLT->GetTPCNcls()    
             		   	,esdTrackHLT->Theta()
             		   	,esdTrackHLT->Eta()	      
             		   	,esdTrackHLT->Phi()
             		   	,dca[0] 		      
             		   	,dca[1] 		      
             		   	,esdTrackHLT->Charge()  
	     		   	,esdTrackHLT->GetNcls(0)
	     		   	,nr_tracksHLT
	     		   	,vertHLT->GetStatus()
	     		   	,vertHLT->GetZ()		       
             		       };
        if(fOptions.Contains("track-hlt")) fTrackHLT->Fill(trackHLT);
      }   
  }  // end of track loop  
  
  if(fBeamType.Contains("Pb")){
     Double_t eventHLT[] = { vertHLT->GetX(), vertHLT->GetY(), vertHLT->GetZ(), 
                             esdHLT->GetPrimaryVertexSPD()->GetX(), esdHLT->GetPrimaryVertexSPD()->GetY(), esdHLT->GetPrimaryVertexSPD()->GetZ(),  
                             vertHLT->GetNContributors(), nr_tracksHLT, vertHLT->GetStatus(),fCentrality->GetCentralityPercentile("V0M")};
     if(fOptions.Contains("event-hlt")) fEventHLT->Fill(eventHLT);  
  }
  else{
     Double_t eventHLT[] = { vertHLT->GetX(), vertHLT->GetY(), vertHLT->GetZ(),
                             esdHLT->GetPrimaryVertexSPD()->GetX(), esdHLT->GetPrimaryVertexSPD()->GetY(), esdHLT->GetPrimaryVertexSPD()->GetZ(),  
                             vertHLT->GetNContributors(), nr_tracksHLT, vertHLT->GetStatus()};
     if(fOptions.Contains("event-hlt")) fEventHLT->Fill(eventHLT);
  }
             
  // Post output data.
  PostData(1, fOutputList);
}

void AliAnalysisTaskHLTCentralBarrel::Terminate(Option_t *){
  // see header file of AliAnalysisTask for documentation 
}

THnSparseF* AliAnalysisTaskHLTCentralBarrel::CreateEventTHnSparse(const char* name, Int_t size, const Int_t* bins, Double_t* min, Double_t* max){
//see header for documentation                     
  
  THnSparseF *thn = new THnSparseF(name,"",size,bins,min,max);
  thn->GetAxis(0)->SetTitle("primary vertex x");
  thn->GetAxis(1)->SetTitle("primary vertex y");
  thn->GetAxis(2)->SetTitle("primary vertex z"); 
  thn->GetAxis(3)->SetTitle("SPD primary vertex x");
  thn->GetAxis(4)->SetTitle("SPD primary vertex y");
  thn->GetAxis(5)->SetTitle("SPD primary vertex z");
  thn->GetAxis(6)->SetTitle("number of contributors");
  thn->GetAxis(7)->SetTitle("track multiplicity");
  thn->GetAxis(8)->SetTitle("vertex status"); 
  if(fBeamType.Contains("Pb")) thn->GetAxis(9)->SetTitle("V0 centrality");
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
  thn->GetAxis(8)->SetTitle("ITS clusters/track");  
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
