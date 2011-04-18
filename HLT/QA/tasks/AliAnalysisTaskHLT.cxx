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

/** @file  AliAnalysisTaskHLT.cxx  
    @author Kalliopi Kanaki, Hege Erdal
    @date 
    @brief An analysis task containing
    loops over HLT and offline ESD trees for comparing
    AliESDtrack properties.    
*/

//#include <iostream>

class AliAnalysisTask;
class AliAnalysisManager;

#include "AliAnalysisTaskHLT.h"
#include "AliHLTGlobalTriggerDecision.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliESDInputHandler.h"
#include "AliTracker.h" 
#include "AliCentrality.h"

#include "TText.h"
#include "TTimeStamp.h"

#include <iostream>

ClassImp(AliAnalysisTaskHLT)

//===============================================================//

AliAnalysisTaskHLT::AliAnalysisTaskHLT()
:
AliAnalysisTaskSE()
  ,fUseHLTTrigger(kFALSE)
  //,fESDOfftrackCuts(0)
  //,fESDHLTtrackCuts(0)
  ,fOutputList(0)
  //,fHistTrigger(0)
  ,fChargeOff(0)  
  ,fMomentumOff(0)
  ,fDCArOff(0)  	
  ,fDCAzOff(0)  	
  ,fNclusterOff(0)
  ,fNITSclusterOff(0)
  ,fPhiOff(0)  	
  ,fEtaOff(0)
  ,fMultOff(0) 	
  ,fXvertexOff(0)	    
  ,fYvertexOff(0)	    
  ,fZvertexOff(0)
  ,fSPDXvertexOff(0)	    
  ,fSPDYvertexOff(0)	    
  ,fSPDZvertexOff(0)
  ,fEventSpecieOff(0)
  ,fV0cent(0)  
  
  ,fChargeHLT(0)
  ,fMomentumHLT(0)
  ,fDCArHLT(0)  
  ,fDCAzHLT(0)  
  ,fNclusterHLT(0)
  ,fNITSclusterHLT(0)
  ,fPhiHLT(0)     
  ,fEtaHLT(0)
  
  ,fChargeHLTcut(0)
  ,fMomentumHLTcut(0)
  ,fDCArHLTcut(0)  
  ,fDCAzHLTcut(0)  
  ,fNclusterHLTcut(0)
  ,fNITSclusterHLTcut(0)
  ,fPhiHLTcut(0)     
  ,fEtaHLTcut(0)
  
  ,fMultHLT(0)  
  ,fXvertexHLT(0)
  ,fYvertexHLT(0)
  ,fZvertexHLT(0)
  ,fSPDXvertexHLT(0)	 
  ,fSPDYvertexHLT(0)	 
  ,fSPDZvertexHLT(0)
  ,fEventSpecieHLT(0)
  
  ,fBeamType()
  ,fTextBox(0)
  ,fCuts(0)
  ,fSwitch(kTRUE)
  ,fCentrality()
  ,fEta(0)
  ,fPt(0)
  ,fDCAr(0)
  ,fDCAz(0)
{
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  // DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
}
 
AliAnalysisTaskHLT::AliAnalysisTaskHLT(const char *name, float eta, float pt, float DCAr, float DCAz)
  :
  AliAnalysisTaskSE(name) 
  ,fUseHLTTrigger(kFALSE)   
  //,fESDOfftrackCuts(0)
  //,fESDHLTtrackCuts(0)
  ,fOutputList(0)
  //,fHistTrigger(0)
  ,fChargeOff(0)  
  ,fMomentumOff(0)
  ,fDCArOff(0) 
  ,fDCAzOff(0) 
  ,fNclusterOff(0)
  ,fNITSclusterOff(0)
  ,fPhiOff(0)  	
  ,fEtaOff(0)
  ,fMultOff(0) 	
  ,fXvertexOff(0)	    
  ,fYvertexOff(0)	    
  ,fZvertexOff(0)
  ,fEventSpecieOff(0)
  ,fV0cent(0)  
  ,fNcontOff(0)
  
  ,fChargeHLT(0)      
  ,fMomentumHLT(0)
  ,fNclusterHLT(0)
  ,fNITSclusterHLT(0)
  ,fPhiHLT(0)     
  ,fEtaHLT(0)
  
  ,fChargeHLTcut(0)      
  ,fMomentumHLTcut(0)
  ,fDCArHLTcut(0)  
  ,fDCAzHLTcut(0)  
  ,fNclusterHLTcut(0)
  ,fNITSclusterHLTcut(0)
  ,fPhiHLTcut(0)     
  ,fEtaHLTcut(0)
  
  ,fMultHLT(0)  
  ,fXvertexHLT(0)
  ,fYvertexHLT(0)
  ,fZvertexHLT(0)
  ,fEventSpecieHLT(0)
  ,fNcontHLT(0)
  
  ,fBeamType()
  ,fTextBox(0)
  ,fCuts(0)
  ,fSwitch(kTRUE)
  ,fCentrality()
  ,fEta(eta)
  ,fPt(pt)
  ,fDCAr(DCAr)
  ,fDCAz(DCAz)
{ 
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  // DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(1, TList::Class());
}

//------------------------------------------------------------------------//

void AliAnalysisTaskHLT::UserCreateOutputObjects(){
  // Create histograms

  OpenFile(1);
  fOutputList = new TList();
  fOutputList->SetOwner();
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

  //fHistTrigger = new TH1F("fHistTrigger", "CTP trigger counter", 24, 0, 24);
  //fHistTrigger->GetXaxis()->SetTitle("");  
  //fHistTrigger->GetYaxis()->SetTitle("#Events"); 

  //=========== event properties =================//

  if(fBeamType.Contains("Pb")){
     fMultOff  = new TH1F("fMult_off", "",200,0,2000);
     fMultHLT  = new TH1F("fMult_hlt", "TPC track multiplicity",200,0,2000);     
     fNcontOff = new TH1F("fNcont_off","",200,0,2000);
     fNcontHLT = new TH1F("fNcont_hlt","# of contributors",200,0,2000);     
     fV0cent   = new TH1F("fV0cent",   "V0 centrality",    100,0, 100);
  } 
  else {
     fMultOff = new TH1F("fMult_off","",200,0,200);
     fMultHLT = new TH1F("fMult_hlt","TPC track multiplicity",200,0,200);    
     fNcontOff = new TH1F("fNcont_off","",200,0,200);
     fNcontHLT = new TH1F("fNcont_hlt","# of contributors",200,0,200);
  }
   
  fXvertexOff = new TH1F("fXvertex_off","",200,-0.5,0.5);
  fXvertexHLT = new TH1F("fXvertex_hlt","X primary vertex",200,-0.5,0.5);
 
  fYvertexOff = new TH1F("fYvertex_off","",200,-0.5,0.5);
  fYvertexHLT = new TH1F("fYvertex_hlt","Y primary vertex",200,-0.5,0.5);
 
  fZvertexOff = new TH1F("fZvertex_off","",100,-20,20);
  fZvertexHLT = new TH1F("fZvertex_hlt","Z primary vertex",100,-20,20);

  fSPDXvertexOff = new TH1F("fSPDXvertex_off","",200,-0.5,0.5);
  fSPDXvertexHLT = new TH1F("fSPDXvertex_hlt","X SPD primary vertex",200,-0.5,0.5);
 
  fSPDYvertexOff = new TH1F("fSPDYvertex_off","",200,-0.5,0.5);
  fSPDYvertexHLT = new TH1F("fSPDYvertex_hlt","Y SPD primary vertex",200,-0.5,0.5);
 
  fSPDZvertexOff = new TH1F("fSPDZvertex_off","",100,-20,20);
  fSPDZvertexHLT = new TH1F("fSPDZvertex_hlt","Z SPD primary vertex",100,-20,20);
    
  fEventSpecieOff = new TH1F("fEventSpecie_off","",18, 0, 18);
  fEventSpecieHLT = new TH1F("fEventSpecie_hlt","event species",18, 0, 18);

  //============== track properties =================//

  fChargeOff = new TH1F("fCharge_off", "", 3, -1.5, 1.5);  
  fChargeHLT = new TH1F("fCharge_hlt", "charge distribution", 3, -1.5, 1.5);  
  
  fMomentumOff = new TH1F("fMomentum_off", "", 100, 0, 10);
  fMomentumHLT = new TH1F("fMomentum_hlt", "transverse momentum", 100, 0, 10);
 
  fDCArOff = new TH1F("fDCAr_off", "", 200, -15, 15);
  fDCArHLT = new TH1F("fDCAr_hlt", "DCAr", 200, -15, 15);

  fDCAzOff = new TH1F("fDCAz_off", "", 200, -15, 15);
  fDCAzHLT = new TH1F("fDCAz_hlt", "DCAz", 200, -15, 15);
 
  fNclusterOff = new TH1F("fNcluster_off","", 200, 0, 200);
  fNclusterHLT = new TH1F("fNcluster_hlt","TPC clusters/track", 200, 0, 200);

  fNITSclusterOff = new TH1F("fNITScluster_off","", 10, 0, 10);
  fNITSclusterHLT = new TH1F("fNITScluster_hlt","ITS clusters/track", 10, 0, 10);
 
  fPhiOff = new TH1F("fPhi_off","",90,0,360);
  fPhiHLT = new TH1F("fPhi_hlt","azimuthal angle",90,0,360);

  fEtaOff = new TH1F("fEta_off","",100,-2,2);
  fEtaHLT = new TH1F("fEta_hlt","pseudorapidity",100,-2,2);
  
  fChargeHLTcut      = new TH1F("fCharge_hltcut",     "",  3, -1.5, 1.5);  
  fMomentumHLTcut    = new TH1F("fMomentum_hltcut",   "",100,    0,  10); 
  fDCArHLTcut        = new TH1F("fDCAr_hltcut",       "",200,  -15,  15);
  fDCAzHLTcut        = new TH1F("fDCAz_hltcut",       "",200,  -15,  15);
  fNclusterHLTcut    = new TH1F("fNcluster_hltcut",   "",200,    0, 200);
  fNITSclusterHLTcut = new TH1F("fNITScluster_hltcut","", 10,    0,  10); 
  fPhiHLTcut         = new TH1F("fPhi_hltcut",        "", 90,    0, 360);
  fEtaHLTcut         = new TH1F("fEta_hltcut",        "",100,   -2,   2);

  //--------------------------------------------------//
  
  fTextBox = new TText(); 
  fCuts = new TText(); 
  fCuts->SetName("cuts");
  TString s = Form("|#eta|<%2g, p_{T}>%2g, |DCAr|<%2g, |DCAz|<%2g", TMath::Abs(fEta), TMath::Abs(fPt), TMath::Abs(fDCAr), TMath::Abs(fDCAz));
  fCuts->SetTitle(s);  
  
  //fOutputList->Add(fHistTrigger);
  fOutputList->Add(fChargeOff); 	  fOutputList->Add(fChargeHLT);       fOutputList->Add(fChargeHLTcut);  
  fOutputList->Add(fMomentumOff);	  fOutputList->Add(fMomentumHLT);     fOutputList->Add(fMomentumHLTcut); 
  fOutputList->Add(fDCArOff);	  	  fOutputList->Add(fDCArHLT);	      fOutputList->Add(fDCArHLTcut);     
  fOutputList->Add(fDCAzOff);	  	  fOutputList->Add(fDCAzHLT);	      fOutputList->Add(fDCAzHLTcut);     
  fOutputList->Add(fNclusterOff);	  fOutputList->Add(fNclusterHLT);     fOutputList->Add(fNclusterHLTcut);   	 
  fOutputList->Add(fNITSclusterOff);	  fOutputList->Add(fNITSclusterHLT);  fOutputList->Add(fNITSclusterHLTcut);  
  fOutputList->Add(fPhiOff);	  	  fOutputList->Add(fPhiHLT);	      fOutputList->Add(fPhiHLTcut);
  fOutputList->Add(fEtaOff);  		  fOutputList->Add(fEtaHLT);  	      fOutputList->Add(fEtaHLTcut);
  
  fOutputList->Add(fMultOff);		  fOutputList->Add(fMultHLT);	 
  fOutputList->Add(fXvertexOff);  	  fOutputList->Add(fXvertexHLT);  
  fOutputList->Add(fYvertexOff);  	  fOutputList->Add(fYvertexHLT);  
  fOutputList->Add(fZvertexOff);  	  fOutputList->Add(fZvertexHLT);    
  fOutputList->Add(fSPDXvertexOff);  	  fOutputList->Add(fSPDXvertexHLT);  
  fOutputList->Add(fSPDYvertexOff);  	  fOutputList->Add(fSPDYvertexHLT);  
  fOutputList->Add(fSPDZvertexOff);  	  fOutputList->Add(fSPDZvertexHLT);    
  fOutputList->Add(fEventSpecieOff);	  fOutputList->Add(fEventSpecieHLT);  
  fOutputList->Add(fNcontOff);	          fOutputList->Add(fNcontHLT);  
  
  fOutputList->Add(fTextBox);
  fOutputList->Add(fCuts);
  if(fBeamType.Contains("Pb")) fOutputList->Add(fV0cent);
 
  //SetupESDtrackCuts();
  PostData(1, fOutputList);
}

void AliAnalysisTaskHLT::UserExec(Option_t *){
  // see header file of AliAnalysisTask for documentation

  AliESDEvent *esdOFF = dynamic_cast<AliESDEvent*>(InputEvent());  
  if(!esdOFF){
     printf("Error:UserExec OFF esd not available\n");
     return;
  }
 
  if(esdOFF->GetEventSpecie()==16) return; // skip calibration events, HLT doesn't set this flag yet

  AliESDInputHandler *esdH = dynamic_cast<AliESDInputHandler*>(fInputHandler);
  if(!esdH){
     printf("The ESD input handler is empty\n");
     return;
  }
  
  AliESDEvent *esdHLT = NULL;   
  if(esdH) esdHLT = esdH->GetHLTEvent();   
  if(!esdHLT){
     printf("Error:UserExec HLT esd not available\n");
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

  //Double_t bfield = esdOFF->GetMagneticField();
 
//   UInt_t Statusnames[12]={AliESDtrack::kTPCin,
// 			  AliESDtrack::kITSin,
// 			  AliESDtrack::kTPCout,
// 			  AliESDtrack::kITSout,
// 			  AliESDtrack::kITSrefit,
// 			  AliESDtrack::kTPCrefit,
// 			  AliESDtrack::kTRDin,
// 			  AliESDtrack::kTRDout,
// 			  AliESDtrack::kTRDrefit,
// 			  AliESDtrack::kTOFin,
// 			  AliESDtrack::kTOFout,
// 			  AliESDtrack::kTOFrefit};
  


  //---------------- HLT ESD tree -----------------------//

  Int_t nr_tracksHLT = 0;	 
  const AliESDVertex *vertHLT = esdHLT->GetPrimaryVertexTracks();

  if(vertHLT->GetStatus()==kTRUE){
    fXvertexHLT->Fill( vertHLT->GetX() );
    fYvertexHLT->Fill( vertHLT->GetY() );
    fZvertexHLT->Fill( vertHLT->GetZ() );
    
    fSPDXvertexHLT->Fill(esdHLT->GetPrimaryVertexSPD()->GetX());
    fSPDYvertexHLT->Fill(esdHLT->GetPrimaryVertexSPD()->GetY());
    fSPDZvertexHLT->Fill(esdHLT->GetPrimaryVertexSPD()->GetZ());
    
    fNcontHLT->Fill(vertHLT->GetNContributors());
  }

  fEventSpecieHLT->Fill(esdHLT->GetEventSpecie());

  for(Int_t i=0; i<esdHLT->GetNumberOfTracks(); i++){ 
  
    AliESDtrack *esdtrackHLT = esdHLT->GetTrack(i); 
    if(!esdtrackHLT) continue;

    //Fill which status flags that are set for an event
    //for(int jjj=0;jjj<12;jjj++){
    //  if(esdtrackHLT->GetStatus()&Statusnames[jjj]) fStatusHLT->Fill(jjj);
    //} 

    if(!(esdtrackHLT->GetStatus()&AliESDtrack::kTPCin)) continue; // only interested in tracks with kTPCin flag
    if(esdtrackHLT->GetTPCNcls()<=0) continue; 
    nr_tracksHLT++;
   
    Double_t x[3]; 
    esdtrackHLT->GetXYZ(x);
    Double_t b[3]; 
    AliTracker::GetBxByBz(x,b);
    Bool_t isOK = esdtrackHLT->RelateToVertexTPCBxByBz(vertHLT, b, kVeryBig);

    Float_t dca[2]={0,0}; Float_t cov[3]={0,0,0}; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z
    if(isOK){
       const AliExternalTrackParam *track = esdtrackHLT->GetTPCInnerParam();
       if(!track) return;    
       esdtrackHLT->GetImpactParametersTPC(dca,cov);
       fDCArHLT->Fill(dca[0]);
       fDCAzHLT->Fill(dca[1]);
    }
 
//     Float_t dca[2];
//     if(vertHLT->GetStatus()==kTRUE){
//        esdtrackHLT->GetDZ(esdHLT->GetPrimaryVertex()->GetXv(), esdHLT->GetPrimaryVertex()->GetYv(), esdHLT->GetPrimaryVertex()->GetZv(), bfield, dca);
//        fDCArHLT->Fill(dca[0]);  
//        fDCAzHLT->Fill(dca[1]);
//     }
    
    fChargeHLT->Fill(esdtrackHLT->Charge());
    fNclusterHLT->Fill(esdtrackHLT->GetTPCNcls());
    fNITSclusterHLT->Fill(esdtrackHLT->GetNcls(0));
    fEtaHLT->Fill(esdtrackHLT->Eta()); 
    fPhiHLT->Fill(esdtrackHLT->Phi()*TMath::RadToDeg());
    fMomentumHLT->Fill(TMath::Abs(esdtrackHLT->Pt()));  
    
    if(TMath::Abs(esdtrackHLT->Eta())<TMath::Abs(fEta) && TMath::Abs(esdtrackHLT->Pt())>TMath::Abs(fPt) && TMath::Abs(dca[0])<TMath::Abs(fDCAr) && TMath::Abs(dca[1])<TMath::Abs(fDCAz)){
       fChargeHLTcut->Fill(esdtrackHLT->Charge());
       fNclusterHLTcut->Fill(esdtrackHLT->GetTPCNcls());
       fDCArHLTcut->Fill(dca[0]);  
       fDCAzHLTcut->Fill(dca[1]);
       fNITSclusterHLTcut->Fill(esdtrackHLT->GetNcls(0));
       fEtaHLTcut->Fill(esdtrackHLT->Eta()); 
       fPhiHLTcut->Fill(esdtrackHLT->Phi()*TMath::RadToDeg());
       fMomentumHLTcut->Fill(TMath::Abs(esdtrackHLT->Pt()));  
    }    
  } // end of loop over hlt tracks

  if(nr_tracksHLT>0) fMultHLT->Fill(nr_tracksHLT);

  //----------------- OFFLINE ESD tree ----------------//
  
  Int_t nr_tracksOff = 0;
  const AliESDVertex *vertOFF = esdOFF->GetPrimaryVertexTracks();
   
  if(fBeamType.Contains("Pb")){
     fCentrality = esdOFF->GetCentrality(); 
     // this information is only available from the offline ESD for 2010 PbPb data, the V0 info was not stored in the HLTesd (17.04.11, Kelly)
     if(!fCentrality){
        printf("Centrality pointer is empty\n");
        return;
     }
     else fV0cent->Fill(fCentrality->GetCentralityPercentile("V0M"));
  }
  
  if(vertOFF->GetStatus()==kTRUE){
     fXvertexOff->Fill( vertOFF->GetX() );
     fYvertexOff->Fill( vertOFF->GetY() );
     fZvertexOff->Fill( vertOFF->GetZ() );
  
     fSPDXvertexOff->Fill(esdOFF->GetPrimaryVertexSPD()->GetX());
     fSPDYvertexOff->Fill(esdOFF->GetPrimaryVertexSPD()->GetY());
     fSPDZvertexOff->Fill(esdOFF->GetPrimaryVertexSPD()->GetZ());
   
     fNcontOff->Fill(vertOFF->GetNContributors());
  }

  fEventSpecieOff->Fill(esdOFF->GetEventSpecie());

  for(Int_t i=0; i<esdOFF->GetNumberOfTracks(); i++){ 
   
    AliESDtrack *esdtrackOFF = esdOFF->GetTrack(i); 
    if (!esdtrackOFF) continue;

    if(!(esdtrackOFF->GetStatus()&AliESDtrack::kTPCin)) continue; 
    if(esdtrackOFF->GetTPCNcls()<=0) continue; 
    nr_tracksOff++;
   
    Double_t x[3]; 
    esdtrackOFF->GetXYZ(x);
    Double_t b[3]; 
    AliTracker::GetBxByBz(x,b);
    Bool_t isOK = esdtrackOFF->RelateToVertexTPCBxByBz(vertOFF, b, kVeryBig);
    if(isOK){
       const AliExternalTrackParam *track = esdtrackOFF->GetTPCInnerParam();
       if(!track) return;    
       Float_t dca[2]={0,0}; Float_t cov[3]={0,0,0}; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z
       esdtrackOFF->GetImpactParametersTPC(dca,cov);
       fDCArOff->Fill(dca[0]);
       fDCAzOff->Fill(dca[1]);
    }
//     // calculate the offline DCAs the same way like for HLT. The above way is calculating the DCA for the TPC inner barrel surface (Kelly, 17.04.11)
//     Float_t dca[2];
//     if(vertOFF->GetStatus()==kTRUE){
//        esdtrackOFF->GetDZ(esdOFF->GetPrimaryVertex()->GetXv(), esdOFF->GetPrimaryVertex()->GetYv(), esdOFF->GetPrimaryVertex()->GetZv(), bfield, dca);
//    }
    
    fChargeOff->Fill(esdtrackOFF->Charge());
    fNclusterOff->Fill(esdtrackOFF->GetTPCNcls()); 
    fNITSclusterOff->Fill(esdtrackOFF->GetNcls(0)); 
    fEtaOff->Fill(esdtrackOFF->Eta());  	
    fPhiOff->Fill(esdtrackOFF->Phi()*TMath::RadToDeg());
    fMomentumOff->Fill( TMath::Abs(esdtrackOFF->Pt()) ); 

  } // end of loop over offline tracks
  
  if(nr_tracksOff>0) fMultOff->Fill(nr_tracksOff);
   
  PostData(1, fOutputList);
}

void AliAnalysisTaskHLT::Terminate(Option_t *){
// see header file of AliAnalysisTask for documentation
}

// void AliAnalysisTaskHLT::SetupESDtrackCuts(){ // not called
//   // Setup ESD cuts
//   // NB! Work in progress!
// 
//   Bool_t selPrimaries = kTRUE;
//   
//   fESDOfftrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2009(selPrimaries);
//   //To make Offline cuts = HLT cuts
//   fESDOfftrackCuts->SetRequireITSRefit(kFALSE); 
//   fESDOfftrackCuts->SetEtaRange(-0.9,+0.9);
//   
// 
//   //HLT
//   //NB! Need to understand this a bit more! Which cuts should we keep?
//   fESDHLTtrackCuts = new AliESDtrackCuts;
// 
//   // TPC  
//   fESDHLTtrackCuts->SetRequireTPCStandAlone(kTRUE); // to get chi2 and ncls of kTPCin
//   fESDHLTtrackCuts->SetMinNClustersTPC(70);
//   fESDHLTtrackCuts->SetMaxChi2PerClusterTPC(4);
//   fESDHLTtrackCuts->SetAcceptKinkDaughters(kFALSE);
// 
//   fESDHLTtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
//   				     AliESDtrackCuts::kAny);
// 
//   if(selPrimaries) { // 7*(0.0050+0.0060/pt^0.9)
//     fESDHLTtrackCuts->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9");
//   }
//   
//   fESDHLTtrackCuts->SetMaxDCAToVertexZ(1.e6);
//   fESDHLTtrackCuts->SetDCAToVertex2D(kFALSE);
//   fESDHLTtrackCuts->SetRequireSigmaToVertex(kFALSE);
//   fESDHLTtrackCuts->SetEtaRange(-0.9,+0.9);
// 
//   return;
// }
