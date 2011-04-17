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

ClassImp(AliAnalysisTaskHLT)

//===============================================================//

AliAnalysisTaskHLT::AliAnalysisTaskHLT()
:
AliAnalysisTaskSE()
  ,fUseHLTTrigger(kFALSE)
  //,fESDOfftrackCuts(0)
  //,fESDHLTtrackCuts(0)
  ,fOutputList(0)
  ,fHistTrigger(0)
  ,fChargeOff(0)  
  ,fMomentumOff(0)
  ,fDCArOff(0)  	
  ,fDCAzOff(0)  	
  ,fNclusterOff(0)
  ,fNITSclusterOff(0)
  ,fNclusterOffwCut(0)		
  ,fPhiOff(0)  	
  ,fMultOff(0) 	
  ,fXYvertexOff(0)	
  ,fXvertexOff(0)	    
  ,fYvertexOff(0)	    
  ,fZvertexOff(0)
  ,fSPDXvertexOff(0)	    
  ,fSPDYvertexOff(0)	    
  ,fSPDZvertexOff(0)
  ,fEtaOff(0)
  ,fEtaMomentumcutOff(0)
  ,fNclusVSphiOff(0)
  ,fNclusVSthetaOff(0)
  ,fEventSpecieOff(0)
  ,fV0cent(0)  
  
  ,fChargeHLT(0)
  ,fMomentumHLT(0)
  ,fDCArHLT(0)  
  ,fDCAzHLT(0)  
  ,fNclusterHLT(0)
  ,fNITSclusterHLT(0)
  ,fNclusterHLTwCut(0)
  ,fPhiHLT(0)     
  ,fMultHLT(0)  
  ,fXYvertexHLT(0)
  ,fXvertexHLT(0)
  ,fYvertexHLT(0)
  ,fZvertexHLT(0)
  ,fSPDXvertexHLT(0)	 
  ,fSPDYvertexHLT(0)	 
  ,fSPDZvertexHLT(0)
  ,fEtaHLT(0)
  ,fEtaMomentumcutHLT(0)
  ,fNclusVSphiHLT(0)	    
  ,fNclusVSthetaHLT(0)
  ,fEventSpecieHLT(0)
  
  ,fBeamType()
  ,fTextBox(0)
  ,fSwitch(kTRUE)
  ,fCentrality()
{
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  // DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
}
 
AliAnalysisTaskHLT::AliAnalysisTaskHLT(const char *name)
  :
  AliAnalysisTaskSE(name) 
  ,fUseHLTTrigger(kFALSE)   
  //,fESDOfftrackCuts(0)
  //,fESDHLTtrackCuts(0)
  ,fOutputList(0)
  ,fHistTrigger(0)
  ,fChargeOff(0)  
  ,fMomentumOff(0)
  ,fDCArOff(0) 
  ,fDCAzOff(0) 
  ,fNclusterOff(0)
  ,fNITSclusterOff(0)
  ,fNclusterOffwCut(0)	
  ,fPhiOff(0)  	
  ,fMultOff(0) 	
  ,fXYvertexOff(0)	
  ,fXvertexOff(0)	    
  ,fYvertexOff(0)	    
  ,fZvertexOff(0)
  ,fEtaOff(0)
  ,fEtaMomentumcutOff(0)
  ,fNclusVSphiOff(0)
  ,fNclusVSthetaOff(0)
  ,fEventSpecieOff(0)
  ,fV0cent(0)  
  ,fNcontOff(0)
  
  ,fChargeHLT(0)      
  ,fMomentumHLT(0)
  ,fNclusterHLT(0)
  ,fNITSclusterHLT(0)
  ,fNclusterHLTwCut(0)
  ,fPhiHLT(0)     
  ,fMultHLT(0)  
  ,fXYvertexHLT(0)
  ,fXvertexHLT(0)
  ,fYvertexHLT(0)
  ,fZvertexHLT(0)
  ,fEtaHLT(0)
  ,fEtaMomentumcutHLT(0)
  ,fNclusVSphiHLT(0)	    
  ,fNclusVSthetaHLT(0)
  ,fEventSpecieHLT(0)
  ,fNcontHLT(0)
  
  ,fBeamType()
  ,fTextBox(0)
  ,fSwitch(kTRUE)
  ,fCentrality()
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

  fHistTrigger = new TH1F("fHistTrigger", "CTP trigger counter", 24, 0, 24);
  fHistTrigger->GetXaxis()->SetTitle("");  
  fHistTrigger->GetYaxis()->SetTitle("#Events"); 

  //=========== event properties =================//

  if(fBeamType.Contains("Pb")){
     fMultOff  = new TH1F("fMult_off", "TPC track multiplicity (OFF)",200,0,2000);
     fMultHLT  = new TH1F("fMult_hlt", "TPC track multiplicity (HLT)",200,0,2000);     
     fNcontOff = new TH1F("fNcont_off","# of contributors (OFF)",200,0,2000);
     fNcontHLT = new TH1F("fNcont_hlt","# of contributors (HLT)",200,0,2000);     
     fV0cent   = new TH1F("fV0cent",   "V0 centrality",               100,0, 100);
  } 
  else {
     fMultOff = new TH1F("fMult_off","TPC track multiplicity (OFF)",200,0,200);
     fMultHLT = new TH1F("fMult_hlt","TPC track multiplicity (HLT)",200,0,200);    
     fNcontOff = new TH1F("fNcont_off","# of contributors (OFF)",200,0,200);
     fNcontHLT = new TH1F("fNcont_hlt","# of contributors (HLT)",200,0,200);
  }
 
  fXYvertexOff = new TH2F("fXYvertex_off","XY primary vertex (OFF)",100,-1,1,100,-1,1);
  fXYvertexHLT = new TH2F("fXYvertex_hlt","XY primary vertex (HLT)",100,-1,1,100,-1,1);
  
  fXvertexOff = new TH1F("fXvertex_off","X primary vertex (OFF)",200,-0.5,0.5);
  fXvertexHLT = new TH1F("fXvertex_hlt","X primary vertex (HLT)",200,-0.5,0.5);
 
  fYvertexOff = new TH1F("fYvertex_off","Y primary vertex (OFF)",200,-0.5,0.5);
  fYvertexHLT = new TH1F("fYvertex_hlt","Y primary vertex (HLT)",200,-0.5,0.5);
 
  fZvertexOff = new TH1F("fZvertex_off","Z primary vertex (OFF)",100,-20,20);
  fZvertexHLT = new TH1F("fZvertex_hlt","Z primary vertex (HLT)",100,-20,20);

  fSPDXvertexOff = new TH1F("fSPDXvertex_off","X SPD primary vertex (OFF)",200,-0.5,0.5);
  fSPDXvertexHLT = new TH1F("fSPDXvertex_hlt","X SPD primary vertex (HLT)",200,-0.5,0.5);
 
  fSPDYvertexOff = new TH1F("fSPDYvertex_off","Y SPD primary vertex (OFF)",200,-0.5,0.5);
  fSPDYvertexHLT = new TH1F("fSPDYvertex_hlt","Y SPD primary vertex (HLT)",200,-0.5,0.5);
 
  fSPDZvertexOff = new TH1F("fSPDZvertex_off","Z SPD primary vertex (OFF)",100,-20,20);
  fSPDZvertexHLT = new TH1F("fSPDZvertex_hlt","Z SPD primary vertex (HLT)",100,-20,20);
    
  fEventSpecieOff = new TH1F("fEventSpecie_off","Eventspecie for OFF",18, 0, 18);
  fEventSpecieHLT = new TH1F("fEventSpecie_hlt","Eventspecie for HLT",18, 0, 18);

  //============== track properties =================//

  fChargeOff = new TH1F("fCharge_off", "Charge distribution (OFF)", 3, -1.5, 1.5);  
  fChargeHLT = new TH1F("fCharge_hlt", "Charge distribution (HLT)", 3, -1.5, 1.5);  
  
  fMomentumOff = new TH1F("fMomentum_off", "p_{T} (OFF)", 100, 0, 10);
  fMomentumHLT = new TH1F("fMomentum_hlt", "p_{T} (HLT)", 100, 0, 10);
 
  fDCArOff = new TH1F("fDCAr_off", "DCAr to beam line (OFF)", 200, -15, 15);
  fDCArHLT = new TH1F("fDCAr_hlt", "DCAr to beam line (HLT)", 200, -15, 15);

  fDCAzOff = new TH1F("fDCAz_off", "DCAz to beam line (OFF)", 200, -15, 15);
  fDCAzHLT = new TH1F("fDCAz_hlt", "DCAz to beam line (HLT)", 200, -15, 15);
 
  fNclusterOff = new TH1F("fNcluster_off","TPC clusters/track (OFF)", 200, 0, 200);
  fNclusterHLT = new TH1F("fNcluster_hlt","TPC clusters/track (HLT)", 200, 0, 200);

  fNITSclusterOff = new TH1F("fNITScluster_off","ITS clusters/track (OFF)", 10, 0, 10);
  fNITSclusterHLT = new TH1F("fNITScluster_hlt","ITS clusters/track (HLT)", 10, 0, 10);
 
  fPhiOff = new TH1F("fPhi_off","azimuthal angle distribution (OFF)",90,0,360);
  fPhiHLT = new TH1F("fPhi_hlt","azimuthal angle distribution (HLT)",    90,0,360);

  fEtaOff = new TH1F("fEta_off","pseudorapidity (OFF)",100,-2,2);
  fEtaHLT = new TH1F("fEta_hlt","pseudorapidity (HLT)",100,-2,2);


 
  fNclusterOffwCut = new TH1F("fNcluster_wcut_off","TPC clusters per track with cuts (OFF)", 200, 0, 200);
  fNclusterHLTwCut = new TH1F("fNcluster_wcut_hlt","TPC clusters per track with cuts (HLT)", 200, 0, 200);

  fEtaMomentumcutOff = new TH1F("fEtaMomentumcut_off","pseudorapidity DCAcut (OFF)",100,-2,2);
  fEtaMomentumcutHLT = new TH1F("fEtaMomentumcut_hlt","pseudorapidity DCAcut (HLT)",    100,-2,2);

  fNclusVSphiOff = new TH2F("fNclus_vs_phi_off","clusters per track vs. #phi (OFF)",360,0,360,160,0,160);
  fNclusVSphiHLT = new TH2F("fNclus_vs_phi_hlt","clusters per track vs. #phi (HLT)",    360,0,360,160,0,160);
  
  fNclusVSthetaOff = new TH2F("fNclus_vs_theta_off","clusters per track vs. #theta (OFF)",180,0,180,160,0,160);
  fNclusVSthetaHLT = new TH2F("fNclus_vs_theta_hlt","clusters per track vs. #theta (HLT)",    180,0,180,160,0,160);
 
  //--------------------------------------------------//
  
  fTextBox = new TText(); 

  fOutputList->Add(fHistTrigger);

  fOutputList->Add(fChargeOff); 	  fOutputList->Add(fChargeHLT);  
  fOutputList->Add(fMomentumOff);	  fOutputList->Add(fMomentumHLT); 
  fOutputList->Add(fDCArOff);	  	  fOutputList->Add(fDCArHLT);	  
  fOutputList->Add(fDCAzOff);	  	  fOutputList->Add(fDCAzHLT);	  
  fOutputList->Add(fNclusterOff);	  fOutputList->Add(fNclusterHLT); 
  fOutputList->Add(fNITSclusterOff);	  fOutputList->Add(fNITSclusterHLT); 
  fOutputList->Add(fNclusterOffwCut);	  fOutputList->Add(fNclusterHLTwCut); 
  fOutputList->Add(fPhiOff);	  	  fOutputList->Add(fPhiHLT);	  
  fOutputList->Add(fMultOff);		  fOutputList->Add(fMultHLT);	 
  fOutputList->Add(fXYvertexOff); 	  fOutputList->Add(fXYvertexHLT); 
  fOutputList->Add(fXvertexOff);  	  fOutputList->Add(fXvertexHLT);  
  fOutputList->Add(fYvertexOff);  	  fOutputList->Add(fYvertexHLT);  
  fOutputList->Add(fZvertexOff);  	  fOutputList->Add(fZvertexHLT);    
  fOutputList->Add(fSPDXvertexOff);  	  fOutputList->Add(fSPDXvertexHLT);  
  fOutputList->Add(fSPDYvertexOff);  	  fOutputList->Add(fSPDYvertexHLT);  
  fOutputList->Add(fSPDZvertexOff);  	  fOutputList->Add(fSPDZvertexHLT);    
  fOutputList->Add(fEtaOff);  		  fOutputList->Add(fEtaHLT);  
  fOutputList->Add(fEtaMomentumcutOff);   fOutputList->Add(fEtaMomentumcutHLT);   
  fOutputList->Add(fNclusVSphiOff);  	  fOutputList->Add(fNclusVSphiHLT);  
  fOutputList->Add(fNclusVSthetaOff);	  fOutputList->Add(fNclusVSthetaHLT);
  fOutputList->Add(fEventSpecieOff);	  fOutputList->Add(fEventSpecieHLT);  
  fOutputList->Add(fNcontOff);	          fOutputList->Add(fNcontHLT);  
  
  fOutputList->Add(fTextBox);
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

  Double_t bfield = esdOFF->GetMagneticField();
 
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
    fXYvertexHLT->Fill(vertHLT->GetX(), vertHLT->GetY() );
    fXvertexHLT->Fill( vertHLT->GetX() );
    fYvertexHLT->Fill( vertHLT->GetY() );
    fZvertexHLT->Fill( vertHLT->GetZ() );
    
    fSPDXvertexHLT->Fill(esdHLT->GetPrimaryVertexSPD()->GetX());
    fSPDYvertexHLT->Fill(esdHLT->GetPrimaryVertexSPD()->GetY());
    fSPDZvertexHLT->Fill(esdHLT->GetPrimaryVertexSPD()->GetZ());
    
    fNcontHLT->Fill(vertHLT->GetNContributors());
  }
  //At the moment no constrains on vertex before filling histograms
  //Should be changed. 

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
 
    //Calculating DCA "old" fashion
    Float_t dca[2];
    esdtrackHLT->GetDZ(esdHLT->GetPrimaryVertex()->GetXv(), esdHLT->GetPrimaryVertex()->GetYv(), esdHLT->GetPrimaryVertex()->GetZv(), bfield, dca);

    fDCArHLT->Fill(dca[0]);  
    fDCAzHLT->Fill(dca[1]);
    
    fChargeHLT->Fill(esdtrackHLT->Charge());
    fNclusterHLT->Fill(esdtrackHLT->GetTPCNcls());
    fNITSclusterHLT->Fill(esdtrackHLT->GetNcls(0));
    fEtaHLT->Fill(esdtrackHLT->Eta()); 
    fPhiHLT->Fill(esdtrackHLT->Phi()*TMath::RadToDeg());
    fMomentumHLT->Fill(TMath::Abs(esdtrackHLT->Pt()));  
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
    fXYvertexOff->Fill(vertOFF->GetX(), vertOFF->GetY() );
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
    if(!isOK) return;
    
    const AliExternalTrackParam *track = esdtrackOFF->GetTPCInnerParam();
    if(!track) return;
    
    Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z
    esdtrackOFF->GetImpactParametersTPC(dca,cov);

    fDCArOff->Fill(dca[0]);
    fDCAzOff->Fill(dca[1]);
    
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
