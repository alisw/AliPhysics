/**************************************************************************
 * Copyright(c) 1998-2008, ALICE Experiment at CERN, All rights reserved. *
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

/////////////////////////////////////////////////////////////
//							   //
//	Class to analyze ZDC data			   //
//							   //
/////////////////////////////////////////////////////////////

#include <TList.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TString.h>
#include <TCanvas.h>

#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "AliESD.h"
#include "AliESDEvent.h"
#include "AliESDHeader.h"
#include "AliESDInputHandler.h"
#include "AliESDZDC.h"
#include "AliMultiplicity.h"
#include "AliAODHandler.h"
#include "AliAODEvent.h"
#include "AliAODHeader.h"
#include "AliAODVertex.h"
#include "AliAODVZERO.h"
#include "AliAODZDC.h"
#include "AliAODMCHeader.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliHeader.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliGenEventHeader.h"
#include "AliGenHijingEventHeader.h"
#include "AliPhysicsSelectionTask.h"
#include "AliPhysicsSelection.h"
#include "AliBackgroundSelection.h"
#include "AliTriggerAnalysis.h"
#include "AliCentrality.h"
#include "AliAnalysisTaskZDCpA.h"

ClassImp(AliAnalysisTaskZDCpA)


//________________________________________________________________________
AliAnalysisTaskZDCpA::AliAnalysisTaskZDCpA():
  AliAnalysisTaskSE(),
    fAnalysisInput(kESD),
    fIsMCInput(kFALSE),
    fCentrLowLim(0),
    fCentrUpLim(0), 
    fCentrEstimator(0),
    fOutput(0x0),
    fhZNCvsZNA(0x0),    
    fhZPCvsZPA(0x0), 
    fhZNCvsZPC(0x0),
    fhZNAvsZPA(0x0),
    fhZNAvsZEM(0x0),      
    fhZNCvsZEM(0x0),      
    fhZNAwV0M(0x0),      
    fhZNCwV0M(0x0),      
    fhZNAvsVZERO(0x0),    
    fhZNCvsVZERO(0x0),    
    fhZDCAvsTracklets(0x0),   
    fhZDCCvsTracklets(0x0),   
    fhVZEROvsZEM(0x0),   
    fhDebunch(0x0),      
    fhZNCcentroid(0x0),  
    fhZNAcentroid(0x0),
    fhZNCvscentrality(0x0),
    fhZNAvscentrality(0x0),
    fhZPCvscentrality(0x0),
    fhZPAvscentrality(0x0),
    fhZNCpmcvscentr(0x0),
    fhZNApmcvscentr(0x0),
    fhZPCpmcvscentr(0x0),
    fhZPApmcvscentr(0x0)
{   
   // Default constructor
   for(int i=0; i<5; i++){
     fhZNCPM[i] = 0x0;
     fhZNAPM[i] = 0x0;    
     fhZPCPM[i] = 0x0;	    
     fhZPAPM[i] = 0x0; 	    
     fhZNCPMlg[i] = 0x0;    
     fhZNAPMlg[i] = 0x0;    
     fhZPCPMlg[i] = 0x0;    
     fhZPAPMlg[i] = 0x0;    
   }
   for(int i=0; i<2; i++) fhZEM[i] = 0x0;  
   for(int i=0; i<6; i++){
     fhTDCraw[i] = 0x0;  
     fhTDC[i] = 0x0;  
   }
   for(int i=0; i<4; i++){
     fhZNCPMQiPMC[i] = 0x0;
     fhZNAPMQiPMC[i] = 0x0;
     fhZPCPMQiPMC[i] = 0x0;
     fhZPAPMQiPMC[i] = 0x0;
     fhPMCvsPMQ[i] = 0x0;
   }  
}   

//________________________________________________________________________
AliAnalysisTaskZDCpA::AliAnalysisTaskZDCpA(const char *name):
  AliAnalysisTaskSE(name),
    fAnalysisInput(kESD),
    fIsMCInput(kFALSE),
    fCentrLowLim(0.),
    fCentrUpLim(100.), 
    fCentrEstimator("V0M"),
    fOutput(0x0),
    fhZNCvsZNA(0x0),    
    fhZPCvsZPA(0x0), 
    fhZNCvsZPC(0x0),
    fhZNAvsZPA(0x0),
    fhZNAvsZEM(0x0),      
    fhZNCvsZEM(0x0),      
    fhZNAwV0M(0x0),    
    fhZNCwV0M(0x0),    
    fhZNAvsVZERO(0x0),    
    fhZNCvsVZERO(0x0),    
    fhZDCAvsTracklets(0x0),   
    fhZDCCvsTracklets(0x0),   
    fhVZEROvsZEM(0x0),   
    fhDebunch(0x0),      
    fhZNCcentroid(0x0),  
    fhZNAcentroid(0x0),
    fhZNCvscentrality(0x0),
    fhZNAvscentrality(0x0),
    fhZPCvscentrality(0x0),
    fhZPAvscentrality(0x0),
    fhZNCpmcvscentr(0x0),
    fhZNApmcvscentr(0x0),
    fhZPCpmcvscentr(0x0),
    fhZPApmcvscentr(0x0)
{

   for(int i=0; i<5; i++){
     fhZNCPM[i] = 0x0;
     fhZNAPM[i] = 0x0;    
     fhZPCPM[i] = 0x0;	    
     fhZPAPM[i] = 0x0; 	    
     fhZNCPMlg[i] = 0x0;    
     fhZNAPMlg[i] = 0x0;    
     fhZPCPMlg[i] = 0x0;    
     fhZPAPMlg[i] = 0x0;    
   }
   for(int i=0; i<2; i++) fhZEM[i] = 0x0;  
   for(int i=0; i<6; i++){
     fhTDCraw[i] = 0x0;  
     fhTDC[i] = 0x0;  
   }
   for(int i=0; i<4; i++){
     fhZNCPMQiPMC[i] = 0x0;
     fhZNAPMQiPMC[i] = 0x0;
     fhZPCPMQiPMC[i] = 0x0;
     fhZPAPMQiPMC[i] = 0x0;
     fhPMCvsPMQ[i] = 0x0;
   }  
  
  // Output slot #1 writes into a TList container
  DefineOutput(1, TList::Class()); 

}

//________________________________________________________________________
AliAnalysisTaskZDCpA& AliAnalysisTaskZDCpA::operator=(const AliAnalysisTaskZDCpA& c)
{
  //
  // Assignment operator
  //
  if (this!=&c) {
    AliAnalysisTaskSE::operator=(c);
  }
  return *this;
}

//________________________________________________________________________
AliAnalysisTaskZDCpA::AliAnalysisTaskZDCpA(const AliAnalysisTaskZDCpA& ana):
  AliAnalysisTaskSE(ana),
    fAnalysisInput(ana.fAnalysisInput),
    fIsMCInput(ana.fIsMCInput),
    fCentrLowLim(ana.fCentrLowLim),
    fCentrUpLim(ana.fCentrUpLim), 
    fCentrEstimator(ana.fCentrEstimator),
    fOutput(ana.fOutput),
    fhZNCvsZNA(ana.fhZNCvsZNA),    
    fhZPCvsZPA(ana.fhZPCvsZPA), 
    fhZNCvsZPC(ana.fhZNCvsZPC),
    fhZNAvsZPA(ana.fhZNAvsZPA),
    fhZNAvsZEM(ana.fhZNAvsZEM),      
    fhZNCvsZEM(ana.fhZNCvsZEM),      
    fhZNAwV0M(ana.fhZNAwV0M),    
    fhZNCwV0M(ana.fhZNCwV0M),    
    fhZNAvsVZERO(ana.fhZNAvsVZERO),    
    fhZNCvsVZERO(ana.fhZNCvsVZERO),    
    fhZDCAvsTracklets(ana.fhZDCAvsTracklets), 
    fhZDCCvsTracklets(ana.fhZDCCvsTracklets), 
    fhVZEROvsZEM(ana.fhVZEROvsZEM),   
    fhDebunch(ana.fhDebunch),      
    fhZNCcentroid(ana.fhZNCcentroid),  
    fhZNAcentroid(ana.fhZNAcentroid),
    fhZNCvscentrality(ana.fhZNCvscentrality),
    fhZNAvscentrality(ana.fhZNAvscentrality),
    fhZPCvscentrality(ana.fhZPCvscentrality),
    fhZPAvscentrality(ana.fhZPAvscentrality),
    fhZNCpmcvscentr(ana.fhZNCpmcvscentr),
    fhZNApmcvscentr(ana.fhZNApmcvscentr),
    fhZPCpmcvscentr(ana.fhZPCpmcvscentr),
    fhZPApmcvscentr(ana.fhZPApmcvscentr)
{
  //
  // Copy Constructor	
  //
  for(int i=0; i<5; i++){
     fhZNCPM[i] = ana.fhZNCPM[i];
     fhZNAPM[i] = ana.fhZNAPM[i];    
     fhZPCPM[i] = ana.fhZPCPM[i];     
     fhZPAPM[i] = ana.fhZPAPM[i];	      
     fhZNCPMlg[i] = ana.fhZNCPMlg[i];	 
     fhZNAPMlg[i] = ana.fhZNAPMlg[i];	 
     fhZPCPMlg[i] = ana.fhZPCPMlg[i];	 
     fhZPAPMlg[i] = ana.fhZPAPMlg[i];	 
   }
   for(int i=0; i<2; i++) fhZEM[i] = ana.fhZEM[i];  
   for(int i=0; i<6; i++){
     fhTDCraw[i] = ana.fhTDCraw[i];  
     fhTDC[i] = ana.fhTDC[i];  
   }
   for(int i=0; i<4; i++) fhPMCvsPMQ[i] = ana.fhPMCvsPMQ[i];  
   for(int i=0; i<4; i++){
     fhZNCPMQiPMC[i] = ana.fhZNCPMQiPMC[i];
     fhZNAPMQiPMC[i] = ana.fhZNAPMQiPMC[i];
     fhZPCPMQiPMC[i] = ana.fhZPCPMQiPMC[i];
     fhZPAPMQiPMC[i] = ana.fhZPAPMQiPMC[i];
     fhPMCvsPMQ[i] = ana.fhPMCvsPMQ[i];
   }  
}
 
//________________________________________________________________________
AliAnalysisTaskZDCpA::~AliAnalysisTaskZDCpA()
{
  // Destructor
  if(fOutput && !AliAnalysisManager::GetAnalysisManager()->IsProofMode()){
    delete fOutput; fOutput=0;
  } 
   
}  

//________________________________________________________________________
void AliAnalysisTaskZDCpA::UserCreateOutputObjects()
{
  // Create the output containers
  //printf("  AliAnalysisTaskZDCpA::UserCreateOutputObjects()\n\n");
  fOutput = new TList;
  fOutput->SetOwner();
  //fOutput->SetName("output");
  
  for(int i=0; i<5; i++){
     char hname[20];
     sprintf(hname,"hZNCPM%d",i);
     fhZNCPM[i] = new TH1F(hname, hname, 250, -50., 2450);
     fOutput->Add(fhZNCPM[i]);
     //
     sprintf(hname,"hZNAPM%d",i);
     fhZNAPM[i] = new TH1F(hname, hname, 250, -50., 2450);    
     fOutput->Add(fhZNAPM[i]);
     //
     sprintf(hname,"hZPCPM%d",i);
     fhZPCPM[i] = new TH1F(hname, hname, 250, -50., 2450);	
     fOutput->Add(fhZPCPM[i]);
     //
     sprintf(hname,"hZPAPM%d",i);
     fhZPAPM[i] = new TH1F(hname, hname, 250, -50., 2450);		
     fOutput->Add(fhZPAPM[i]);
     //
     sprintf(hname,"hZNCPMlg%d",i);
     fhZNCPMlg[i] = new TH1F(hname, hname, 200, -50., 950);	 
     fOutput->Add(fhZNCPMlg[i]);
     //
     sprintf(hname,"hZNAPMlg%d",i);
     fhZNAPMlg[i] = new TH1F(hname, hname, 200, -50., 950);	 
     fOutput->Add(fhZNAPMlg[i]);
     //
     sprintf(hname,"hZPCPMlg%d",i);
     fhZPCPMlg[i] = new TH1F(hname, hname, 200, -50., 950);	 
     fOutput->Add(fhZPCPMlg[i]);
     //
     sprintf(hname,"hZPAPMlg%d",i);
     fhZPAPMlg[i] = new TH1F(hname, hname, 200, -50., 950);	 
     fOutput->Add(fhZPAPMlg[i]);
     //
     if(i<4){
       //
       char hnamenc[20];
       sprintf(hnamenc, "hZNCPMQ%dPMC",i+1);
       fhZNCPMQiPMC[i] = new TH1F(hnamenc, hnamenc, 100, 0., 1.);
       fOutput->Add(fhZNCPMQiPMC[i]);
       //
       char hnamena[20];
       sprintf(hnamena, "hZNAPMQ%dPMC",i+1);
       fhZNAPMQiPMC[i] = new TH1F(hnamena, hnamena, 100, 0., 1.);
       fOutput->Add(fhZNAPMQiPMC[i]);
       //
       char hnamepc[20];
       sprintf(hnamepc, "hZPCPMQ%dPMC",i+1);
       fhZPCPMQiPMC[i] = new TH1F(hnamepc, hnamepc, 100, 0., 1.);
       fOutput->Add(fhZPCPMQiPMC[i]);
       //
       char hnamepa[20];
       sprintf(hnamepa, "hZPAPMQ%dPMC",i+1);
       fhZPAPMQiPMC[i] = new TH1F(hnamepa, hnamepa, 100, 0., 1.);
       fOutput->Add(fhZPAPMQiPMC[i]);
     }
   }
   for(int i=0; i<6; i++){
     if(i==0){
       fhZEM[i] = new TH1F("hZEM1","hZEM1",200,-10.,1990.);
       fhTDCraw[i] = new TH1F("hTDCZEM1raw", "hTDCZEM1raw", 200, -200., 0.);
       fhTDC[i] = new TH1F("hTDCZEM1", "hTDCZEM1", 200, -150., 50.);
       fhPMCvsPMQ[i] = new TH2F("hPMCvsPMQZNC","hPMCvsPMQZNC",100,0.,2000,100,0.,2000);
       //
       fOutput->Add(fhZEM[i]);
       fOutput->Add(fhPMCvsPMQ[i]);
     }
     else if(i==1){
       fhZEM[i] = new TH1F("hZEM2","hZEM2",200,-10.,1990.);
       fhTDCraw[i] = new TH1F("hTDCZEM2raw", "hTDCZEM2raw", 200, -200., 0.);
       fhTDC[i] = new TH1F("hTDCZEM2", "hTDCZEM2", 200, -150., 50.);
       fhPMCvsPMQ[i] = new TH2F("hPMCvsPMQZPC","hPMCvsPMQZPC",100,0.,2000,100,0.,2000);
       //
       fOutput->Add(fhZEM[i]);
       fOutput->Add(fhPMCvsPMQ[i]);
     }
     else if(i==2){
       fhTDCraw[i] = new TH1F("hTDCZNCraw", "hTDCZNCraw", 200, -200., 0.);
       fhTDC[i] = new TH1F("hTDCZNC", "hTDCZNC", 200, -150., 50.);
       fhPMCvsPMQ[i] = new TH2F("hPMCvsPMQZNA","hPMCvsPMQZNA",100,0.,2000,100,0.,2000);
       //
       fOutput->Add(fhPMCvsPMQ[i]);
     }
     else if(i==3){
       fhTDCraw[i] = new TH1F("hTDCZPCraw", "hTDCZPCraw", 200, -200., 0.);
       fhTDC[i] = new TH1F("hTDCZPC", "hTDCZPC", 200, -150., 50.);
       fhPMCvsPMQ[i] = new TH2F("hPMCvsPMQZPA","hPMCvsPMQZPA",100,0.,2000,100,0.,2000);
       //
       fOutput->Add(fhPMCvsPMQ[i]);
     }
     else if(i==4){
       fhTDCraw[i] = new TH1F("hTDCZNAraw", "hTDCZNAraw", 200, -200., 0.);
       fhTDC[i] = new TH1F("hTDCZNA", "hTDCZNA", 200, -150., 50.);
     }
     else if(i==5){
       fhTDCraw[i] = new TH1F("hTDCZPAraw", "hTDCZPAraw", 200, -200., 0.);
       fhTDC[i] = new TH1F("hTDCZPA", "hTDCZPA", 200, -150., 50.);
     }
     //
     fOutput->Add(fhTDC[i]);
  }
  
  
  fhZNCvsZNA = new TH2F("hZNCvsZNA","hZNCvsZNA",250,-50.,2450,200,-50.,2450);
  fOutput->Add(fhZNCvsZNA);
  fhZPCvsZPA = new TH2F("hZPCvsZPA","hZPCvsZPA",250,-50.,2450,200,-50.,2450);	
  fOutput->Add(fhZPCvsZPA);
  fhZNCvsZPC = new TH2F("hZNCvsZPC","hZNCvsZPC",250,-50.,2450,200,-50.,2450);
  fOutput->Add(fhZNCvsZPC);
  fhZNAvsZPA = new TH2F("hZNAvsZPA","hZNAvsZPA",250,-50.,2450,200,-50.,2450);	
  fOutput->Add(fhZNAvsZPA);
  fhZNAvsZEM = new TH2F("hZNAvsZEM","hZNAvsZEM",250,0.,2500.,250,0.,2500.);	
  fOutput->Add(fhZNAvsZEM);
  fhZNCvsZEM = new TH2F("hZNCvsZEM","hZNCvsZEM",250,0.,2500.,250,0.,2500.);	
  fOutput->Add(fhZNCvsZEM);
  fhZNAwV0M = new TH1F("hZNAwV0M","hZNA wV0M",250,0.,2500.);	
  fOutput->Add(fhZNAwV0M);
  fhZNCwV0M = new TH1F("hZNCwV0M","hZNC wV0M",250,0.,2500.);	
  fOutput->Add(fhZNCwV0M);
  fhZNAvsVZERO = new TH2F("hZNAvsVZERO","ZNA vs. VZERO",250,0.,25000.,250,0.,2500.);	
  fOutput->Add(fhZNAvsVZERO);
  fhZNCvsVZERO = new TH2F("hZNCvsVZERO","ZNC vs. VZERO",250,0.,25000.,250,0.,2500.);	
  fOutput->Add(fhZNCvsVZERO);
  fhZDCAvsTracklets = new TH2F("hZDCAvsTracklets","hZDCAvsTracklets",200,0.,4000.,250,0.,2500.);
  fOutput->Add(fhZDCAvsTracklets);
  fhZDCCvsTracklets = new TH2F("hZDCCvsTracklets","hZDCCvsTracklets",200,0.,4000.,250,0.,2500.);
  fOutput->Add(fhZDCCvsTracklets);
  fhVZEROvsZEM = new TH2F("hVZEROvsZEM","hVZEROvsZEM",250,0.,2500.,250,0.,25000.);	
  fOutput->Add(fhVZEROvsZEM);
  fhDebunch = new TH2F("hDebunch","hDebunch",240,-100.,-40.,240,-30.,30.);	
  fOutput->Add(fhDebunch);
  fhZNCcentroid = new TH2F("hZNCcentroid","hZNCcentroid",100,-3.5,3.5,100,-3.5,3.5);	
  fOutput->Add(fhZNCcentroid);
  fhZNAcentroid = new TH2F("hZNAcentroid","hZNAcentroid",100,-3.5,3.5,100,-3.5,3.5);	
  fOutput->Add(fhZNAcentroid);
  
  fhZNCvscentrality = new TH2F("hZNCvscentrality","ZNC vs. centrality",100,0.,100.,100,0.,100.);
  fOutput->Add(fhZNCvscentrality);
  fhZNAvscentrality = new TH2F("hZNAvscentrality","ZNA vs. centrality",100,0.,100.,100,0.,100.);
  fOutput->Add(fhZNAvscentrality);
  fhZPCvscentrality = new TH2F("hZPCvscentrality","ZPC vs. centrality",100,0.,100.,100,0.,100.);
  fOutput->Add(fhZPCvscentrality);
  fhZPAvscentrality = new TH2F("hZPAvscentrality","ZPA vs. centrality",100,0.,100.,100,0.,100.);
  fOutput->Add(fhZPAvscentrality);
  
  fhZNCpmcvscentr = new TH2F("hZNCpmcvscentr","ZNC PMC vs. centrality",100,0.,100.,100,0.,100.);
  fOutput->Add(fhZNCpmcvscentr);
  fhZNApmcvscentr = new TH2F("hZNApmcvscentr","ZNA PMC vs. centrality",100,0.,100.,100,0.,100.);
  fOutput->Add(fhZNApmcvscentr);
  fhZPCpmcvscentr = new TH2F("hZPCpmcvscentr","ZPC PMC vs. centrality",100,0.,100.,100,0.,100.);
  fOutput->Add(fhZPCpmcvscentr);
  fhZPApmcvscentr = new TH2F("hZPApmcvscentr","ZPA PMC vs. centrality",100,0.,100.,100,0.,100.);
  fOutput->Add(fhZPApmcvscentr);
 
  PostData(1, fOutput);
}

//________________________________________________________________________
void AliAnalysisTaskZDCpA::UserExec(Option_t */*option*/)
{
  // Execute analysis for current event:

  if(!InputEvent()){
    printf("ERROR: InputEvent not available");
    return;
  }
  
  if(fAnalysisInput==kESD){
      
      AliESDEvent* esd = dynamic_cast<AliESDEvent*> (InputEvent());
      if(!esd) return;
      
      // Select PHYSICS events (type=7, for data)
      if(!fIsMCInput && esd->GetEventType()!=7) return; 

      // ********* MC INFO *********************************
      /*if(fIsMCInput){
        AliMCEventHandler* eventHandler = dynamic_cast<AliMCEventHandler*> (AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler());
        if(!eventHandler){
          printf("ERROR: Could not retrieve MC event handler");
          return;
        }
   
        AliMCEvent* mcEvent = eventHandler->MCEvent();
        if(!mcEvent){
          printf("ERROR: Could not retrieve MC event");
          return;
        }
      }*/
      // ****************************************************
     
     // PHYSICS SELECTION
//     AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
//     AliInputEventHandler *hdr = (AliInputEventHandler*)am->GetInputEventHandler();
//     printf("  PS selection %d\n",hdr->IsEventSelected());
//     if(hdr->IsEventSelected() & AliVEvent::kAnyINT){      
      
      AliCentrality *centrality = esd->GetCentrality();
      Float_t centrperc = centrality->GetCentralityPercentile(fCentrEstimator.Data());
      if(centrperc<fCentrLowLim || centrperc>fCentrUpLim) return;
      
      // ***** Trigger selection
      /*TString triggerClass = esd->GetFiredTriggerClasses();
      sprintf(fTrigClass,"%s",triggerClass.Data());*/
          
      const AliMultiplicity *mult = esd->GetMultiplicity();
      Int_t nTracklets = mult->GetNumberOfTracklets();
      //for(Int_t ilay=0; ilay<6; ilay++) fNClusters[ilay] = mult->GetNumberOfITSClusters(ilay);
      //Int_t nClusterslay1 = mult->GetNumberOfITSClusters(1);
      
      AliESDVZERO *vzeroAOD = esd->GetVZEROData();
      Float_t multV0A = vzeroAOD->GetMTotV0A();
      Float_t multV0C = vzeroAOD->GetMTotV0C();

      AliESDZDC *esdZDC = esd->GetESDZDC();
      
      Float_t tdc[32][4], tdcSum=-999, tdcDiff=-999;
      Bool_t tdcMult[6] = {kFALSE};
      for(Int_t itdc=0; itdc<32; itdc++){
	 for(Int_t i=0; i<4; i++){
	   if(esdZDC->GetZDCTDCData(itdc, i)!=0.){ 
	     tdc[itdc][i] = esdZDC->GetZDCTDCCorrected(itdc, i);
	     if(itdc>=8 && itdc<=13){
	       tdcMult[itdc-8] = kTRUE;
	       fhTDC[itdc-8]->Fill(tdc[itdc][i]);
	     }
	     else if(itdc==15 && i==0){
	       for(int ich=8; ich<=13; ich++) fhTDCraw[ich-8]->Fill(esdZDC->GetZDCTDCData(ich, 0)-esdZDC->GetZDCTDCData(itdc, i));
	     }
           }
	 }
      }
      Bool_t zem1 = tdcMult[0];
      Bool_t zem2 = tdcMult[1];
      Bool_t znc  = tdcMult[2];
      Bool_t zpc  = tdcMult[3];
      Bool_t zna  = tdcMult[4];
      Bool_t zpa  = tdcMult[5];
      if(znc && zna){
	tdcSum = esdZDC->GetZNTDCSum(0);
	tdcDiff = esdZDC->GetZNTDCDiff(0);
	fhDebunch->Fill(tdcDiff, tdcSum);	   
      }
      
      Float_t energyZNC  = (Float_t) (esdZDC->GetZNCEnergy());
      Float_t energyZPC  = (Float_t) (esdZDC->GetZPCEnergy());
      Float_t energyZNA  = (Float_t) (esdZDC->GetZNAEnergy());
      Float_t energyZPA  = (Float_t) (esdZDC->GetZPAEnergy());
      Float_t energyZEM1 = (Float_t) (esdZDC->GetZEM1Energy());
      Float_t energyZEM2 = (Float_t) (esdZDC->GetZEM2Energy());
      
      const Double_t * towZNC = esdZDC->GetZNCTowerEnergy();
      const Double_t * towZPC = esdZDC->GetZPCTowerEnergy();
      const Double_t * towZNA = esdZDC->GetZNATowerEnergy();
      const Double_t * towZPA = esdZDC->GetZPATowerEnergy();
      //
      const Double_t * towZNClg = esdZDC->GetZNCTowerEnergyLR();
      const Double_t * towZPClg = esdZDC->GetZPCTowerEnergyLR();
      const Double_t * towZNAlg = esdZDC->GetZNATowerEnergyLR();
      const Double_t * towZPAlg = esdZDC->GetZPATowerEnergyLR();
      
      Double_t xyZNC[2]={-99.,-99.}, xyZNA[2]={-99.,-99.};
      esdZDC->GetZNCentroidInPbPb(esd->GetBeamEnergy(), xyZNC, xyZNA);    
      fhZNCcentroid->Fill(xyZNC[0], xyZNC[1]);
      fhZNAcentroid->Fill(xyZNA[0], xyZNA[1]);
      
      /* //Only for MC
      fbZDC = esdZDC->GetImpactParameter();
      fNpartZDC = esdZDC->GetZDCParticipants();
      fbZDCA = esdZDC->GetImpactParamSideA();
      fNpartZDCA = esdZDC->GetZDCPartSideA();
      fbZDCC = esdZDC->GetImpactParamSideC();
      fNpartZDCC = esdZDC->GetZDCPartSideC();*/
      
      if(znc){
        for(int i=0; i<5; i++){
	   fhZNCPM[i]->Fill(towZNC[i]);
	   fhZNCPMlg[i]->Fill(towZNClg[i]);
	   if((i<4) && (towZNC[0]>0.)) fhZNCPMQiPMC[i]->Fill(towZNC[i+1]/towZNC[0]);
	}
	fhPMCvsPMQ[0]->Fill(towZNC[1]+towZNC[2]+towZNC[3]+towZNC[4], towZNC[0]);
      }
      if(zpc){
        for(int i=0; i<5; i++){
	   fhZPCPM[i]->Fill(towZPC[i]);
	   fhZPCPMlg[i]->Fill(towZPClg[i]);
	   if(((i<4) && towZPC[0]>0.)) fhZPCPMQiPMC[i]->Fill(towZPC[i+1]/towZPC[0]);
	}
	fhPMCvsPMQ[1]->Fill(towZPC[1]+towZPC[2]+towZPC[3]+towZPC[4], towZPC[0]);
      }
      if(zna){
        for(int i=0; i<5; i++){
	   fhZNAPM[i]->Fill(towZNA[i]);
	   fhZNAPMlg[i]->Fill(towZNAlg[i]);
	   if(((i<4) && towZNA[0]>0.)) fhZNAPMQiPMC[i]->Fill(towZNA[i+1]/towZNA[0]);
	}
	fhPMCvsPMQ[2]->Fill(towZNA[1]+towZNA[2]+towZNA[3]+towZNA[4], towZNA[0]);
      }
      if(zpa){
        for(int i=0; i<5; i++){
	   fhZPAPM[i]->Fill(towZPA[i]);
	   fhZPAPMlg[i]->Fill(towZPAlg[i]);
	   if(((i<4) && towZPA[0]>0.)) fhZPAPMQiPMC[i]->Fill(towZPA[i+1]/towZPA[0]);
	}
	fhPMCvsPMQ[3]->Fill(towZPA[1]+towZPA[2]+towZPA[3]+towZPA[4], towZPA[0]);
      }
      if(zem1) fhZEM[0]->Fill(energyZEM1);
      if(zem2) fhZEM[1]->Fill(energyZEM2);
      
      if(znc || zna){
        fhZNCvsZNA->Fill(energyZNA, energyZNC);
	fhZPCvsZPA->Fill(energyZPA, energyZPC);
        fhZNCvsZPC->Fill(energyZPC, energyZNC);
	fhZNAvsZPA->Fill(energyZPA, energyZNA);
	fhZNAvsZEM->Fill(energyZEM1+energyZEM2, energyZNA);	
	fhZNCvsZEM->Fill(energyZEM1+energyZEM2, energyZNC);	
	fhZNAwV0M->Fill(energyZNA, centrperc);   
	fhZNCwV0M->Fill(energyZNC, centrperc);   
	fhZNAvsVZERO->Fill(multV0A+multV0C, energyZNA);	
	fhZNCvsVZERO->Fill(multV0A+multV0C, energyZNC);	
	fhZDCAvsTracklets->Fill((Float_t) (nTracklets), energyZNA+energyZPA);
	fhZDCCvsTracklets->Fill((Float_t) (nTracklets), energyZNC+energyZPC);
	fhVZEROvsZEM->Fill(energyZEM1+energyZEM2, multV0A+multV0C);	
        	
	fhZNCvscentrality->Fill(centrperc, energyZNC);
	fhZNAvscentrality->Fill(centrperc, energyZNA);
	fhZPCvscentrality->Fill(centrperc, energyZPC);
	fhZPAvscentrality->Fill(centrperc, energyZPA);
	
	fhZNCpmcvscentr->Fill(centrperc, towZNC[0]);
	fhZNApmcvscentr->Fill(centrperc, towZNA[0]);
	fhZPCpmcvscentr->Fill(centrperc, towZPC[0]);
	fhZPApmcvscentr->Fill(centrperc, towZPA[0]);
	
      }
//    } // PHYSICS SELECTION
      
  }   
  else if(fAnalysisInput==kAOD){
      
    AliAODEvent *aod =  dynamic_cast<AliAODEvent*> (InputEvent());
    if(!aod){
       printf("AODs not available");
       return;
    }
    
    // Select PHYSICS events (type=7, for data)
    //if(!fIsMCInput && aod->GetEventType()!=7) return; 
    
    // PHYSICS SELECTION
    AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
    AliInputEventHandler *hdr = (AliInputEventHandler*)am->GetInputEventHandler();
    if(hdr->IsEventSelected() & AliVEvent::kAnyINT){	  
           
      AliCentrality* centrality = aod->GetCentrality();
      Float_t centrperc = centrality->GetCentralityPercentile(fCentrEstimator.Data());
      
      /*AliAODHeader *aodheader =  aod->GetHeader();
      AliCentrality *centralityp = aodheader->GetCentralityP();
      Float_t centrperc = centralityp->GetCentralityPercentile(fCentrEstimator.Data());
      */

      if(centrperc<fCentrLowLim || centrperc>fCentrUpLim)  return;
     
      // ***** Trigger selection
      /*TString triggerClass = aod->GetFiredTriggerClasses();
      sprintf(fTrigClass,"%s",triggerClass.Data());
      
      const AliAODVertex *vertex = aod->GetPrimaryVertexSPD();
      fxVertex = vertex->GetX();
      fyVertex = vertex->GetY();
      fzVertex = vertex->GetZ();*/
      
      AliAODTracklets *trackl = aod->GetTracklets();
      Int_t nTracklets = trackl->GetNumberOfTracklets();
          
      AliAODVZERO *vzeroAOD = aod->GetVZEROData();
      Float_t multV0A = vzeroAOD->GetMTotV0A();
      Float_t multV0C = vzeroAOD->GetMTotV0C();
        
      AliAODZDC *aodZDC = aod->GetZDCData();
            
      Float_t energyZNC  = (Float_t) (aodZDC->GetZNCEnergy());
      Float_t energyZPC  = (Float_t) (aodZDC->GetZPCEnergy());
      Float_t energyZNA  = (Float_t) (aodZDC->GetZNAEnergy());
      Float_t energyZPA  = (Float_t) (aodZDC->GetZPAEnergy());
      Float_t energyZEM1 = (Float_t) (aodZDC->GetZEM1Energy());
      Float_t energyZEM2 = (Float_t) (aodZDC->GetZEM2Energy());
      
      const Double_t * towZNC = aodZDC->GetZNCTowerEnergy();
      const Double_t * towZPC = aodZDC->GetZPCTowerEnergy();
      const Double_t * towZNA = aodZDC->GetZNATowerEnergy();
      const Double_t * towZPA = aodZDC->GetZPATowerEnergy();
      //
      const Double_t * towZNClg = aodZDC->GetZNCTowerEnergyLR();
      const Double_t * towZNAlg = aodZDC->GetZNATowerEnergyLR();
      //
      Double_t towZPClg[5], towZPAlg[5]={0.};
      for(Int_t it=0; it<5; it++){
         towZPClg[it] = 8*towZPC[it];
         towZPAlg[it] = 8*towZNA[it];  
      }

      Double_t xyZNC[2]={-99.,-99.}, xyZNA[2]={-99.,-99.};
      aodZDC->GetZNCentroidInPbPb(1380., xyZNC, xyZNA);
      fhZNCcentroid->Fill(xyZNC[0], xyZNC[1]);
      fhZNAcentroid->Fill(xyZNA[0], xyZNA[1]);
    
      Float_t tdcSum = aodZDC->GetZDCTimeSum();
      Float_t tdcDiff = aodZDC->GetZDCTimeDiff();
      fhDebunch->Fill(tdcDiff, tdcSum);
      
      for(int i=0; i<5; i++){
         fhZNCPM[i]->Fill(towZNC[i]);
         fhZNCPMlg[i]->Fill(towZNClg[i]);
	 if((i<4) && (towZNC[0]>0.)) fhZNCPMQiPMC[i]->Fill(towZNC[i+1]/towZNC[0]);
      }
      fhPMCvsPMQ[0]->Fill(towZNC[1]+towZNC[2]+towZNC[3]+towZNC[4], towZNC[0]);
      for(int i=0; i<5; i++){
         fhZPCPM[i]->Fill(towZPC[i]);
         fhZPCPMlg[i]->Fill(towZPClg[i]);
	 if(((i<4) && towZPC[0]>0.)) fhZPCPMQiPMC[i]->Fill(towZPC[i+1]/towZPC[0]);
      }
      fhPMCvsPMQ[1]->Fill(towZPC[1]+towZPC[2]+towZPC[3]+towZPC[4], towZPC[0]);
      for(int i=0; i<5; i++){
         fhZNAPM[i]->Fill(towZNA[i]);
         fhZNAPMlg[i]->Fill(towZNAlg[i]);
	 if(((i<4) && towZNA[0]>0.)) fhZNAPMQiPMC[i]->Fill(towZNA[i+1]/towZNA[0]);
      }
      fhPMCvsPMQ[2]->Fill(towZNA[1]+towZNA[2]+towZNA[3]+towZNA[4], towZNA[0]);
      for(int i=0; i<5; i++){
         fhZPAPM[i]->Fill(towZPA[i]);
         fhZPAPMlg[i]->Fill(towZPAlg[i]);
	 if(((i<4) && towZPA[0]>0.)) fhZPAPMQiPMC[i]->Fill(towZPA[i+1]/towZPA[0]);
      }
      fhPMCvsPMQ[3]->Fill(towZPA[1]+towZPA[2]+towZPA[3]+towZPA[4], towZPA[0]);
      fhZEM[0]->Fill(energyZEM1);
      fhZEM[1]->Fill(energyZEM2);
      
      fhZNCvsZNA->Fill(energyZNA, energyZNC);
      fhZPCvsZPA->Fill(energyZPA, energyZPC);
      fhZNCvsZPC->Fill(energyZPC, energyZNC);
      fhZNAvsZPA->Fill(energyZPA, energyZNA);
      fhZNAvsZEM->Fill(energyZEM1+energyZEM2, energyZNA);    
      fhZNCvsZEM->Fill(energyZEM1+energyZEM2, energyZNC);    
      fhZNAwV0M->Fill(energyZNA, centrperc);	     
      fhZNCwV0M->Fill(energyZNC, centrperc);	
      fhZNAvsVZERO->Fill(multV0A+multV0C, energyZNA);        
      fhZNCvsVZERO->Fill(multV0A+multV0C, energyZNC);        
      fhZDCAvsTracklets->Fill((Float_t) (nTracklets), energyZNA+energyZPA);
      fhZDCCvsTracklets->Fill((Float_t) (nTracklets), energyZNC+energyZPC);
      fhVZEROvsZEM->Fill(energyZEM1+energyZEM2, multV0A+multV0C);     
      
      fhZNCvscentrality->Fill(centrperc, energyZNC);
      fhZNAvscentrality->Fill(centrperc, energyZNA);
      fhZPCvscentrality->Fill(centrperc, energyZPC);
      fhZPAvscentrality->Fill(centrperc, energyZPA);
	
      fhZNCpmcvscentr->Fill(centrperc, towZNC[0]);
      fhZNApmcvscentr->Fill(centrperc, towZNA[0]);
      fhZPCpmcvscentr->Fill(centrperc, towZPC[0]);
      fhZPApmcvscentr->Fill(centrperc, towZPA[0]);

    } // PHYSICS SELECTION
  }
  
  PostData(1, fOutput);
   
}



//________________________________________________________________________
void AliAnalysisTaskZDCpA::Terminate(Option_t */*option*/)
{
  // Terminate analysis
  //
/*  if(fDebug > 1) printf(" **** AliAnalysisTaskZDCpA::Terminate() \n");
  
  //fOutput = dynamic_cast<TList*> (GetOutputData(1));
  //if(!fOutput) printf("ERROR: fOutput not available\n");
*/
}
