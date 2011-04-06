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


#include "AliHLTGlobalTriggerDecision.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliESDInputHandler.h"
#include "AliTracker.h" 

#include "AliAnalysisTaskHLT.h"

ClassImp(AliAnalysisTaskHLT)

//======================================================================================================

AliAnalysisTaskHLT::AliAnalysisTaskHLT()
:
AliAnalysisTaskSE()
  ,fUseHLTTrigger(kFALSE)
  ,fESDOfftrackCuts(0)
  ,fESDHLTtrackCuts(0)
  ,fOutputList(0)
  ,fHistTrigger(0)
  ,fHistHLTTrigger(0)  
  ,fChargeOff(0)  
  ,fMomentumOff(0)
  ,fMomentumOffTpc(0)
  ,fMomentumOffTpcIts(0)	
  ,fDCArOff(0)  	
  ,fDCAzOff(0)  	
  ,fNclusterOff(0)
  ,fNclusterOffwCut(0)		
  ,fdEdxOff(0) 	
  ,fdEdxVSPOff(0)	
  ,fPhiOff(0)  	
  ,fThetaOff(0)	
  ,fMultOff(0) 	
  ,fXYvertexOff(0)	
  ,fXvertexOff(0)	    
  ,fYvertexOff(0)	    
  ,fZvertexOff(0)
  ,fEtaOff(0)
  ,fEtaMomentumcutOff(0)
  ,fNclusVSphiOff(0)
  ,fNclusVSthetaOff(0)
  ,fStatusOff(0)
  ,fEventSpecieOff(0)
  
  ,fChargeHLT(0)
  ,fMomentumHLT(0)
  ,fMomentumHLTTpc(0)
  ,fMomentumHLTTpcIts(0)
  ,fDCArHLT(0)  
  ,fDCAzHLT(0)  
  ,fDCArHLTSG(0)  
  ,fDCAzHLTSG(0)  
  ,fNclusterHLT(0)
  ,fNclusterHLTwCut(0)
  ,fdEdxHLT(0)    
  ,fdEdxVSPHLT(0)
  ,fPhiHLT(0)     
  ,fThetaHLT(0)  
  ,fMultHLT(0)  
  ,fXYvertexHLT(0)
  ,fXvertexHLT(0)
  ,fYvertexHLT(0)
  ,fZvertexHLT(0)
  ,fEtaHLT(0)
  ,fEtaMomentumcutHLT(0)
  ,fNclusVSphiHLT(0)	    
  ,fNclusVSthetaHLT(0)
  ,fStatusHLT(0)
  ,fEventSpecieHLT(0)
  ,fTrgClsArray(0)
  
  //     ,fDCArOff_trig(0)
  //     ,fNclusterOff_trig(0)
  //     
  //     ,fDCAHLT_trig(0)
  //     ,fNclusterHLT_trig(0)

{

  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  // DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container

  // DefineOutput(1, TList::Class());
}
 
AliAnalysisTaskHLT::AliAnalysisTaskHLT(const char *name)
  :
  AliAnalysisTaskSE(name) 
  ,fUseHLTTrigger(kFALSE)   
  ,fESDOfftrackCuts(0)
  ,fESDHLTtrackCuts(0)
  ,fOutputList(0)
  ,fHistTrigger(0)
  ,fHistHLTTrigger(0)  
  ,fChargeOff(0)  
  ,fMomentumOff(0)
  ,fMomentumOffTpc(0)
  ,fMomentumOffTpcIts(0)	
  ,fDCArOff(0) 
  ,fDCAzOff(0) 
  ,fNclusterOff(0)
  ,fNclusterOffwCut(0)	
  ,fdEdxOff(0) 	
  ,fdEdxVSPOff(0)	
  ,fPhiOff(0)  	
  ,fThetaOff(0)	
  ,fMultOff(0) 	
  ,fXYvertexOff(0)	
  ,fXvertexOff(0)	    
  ,fYvertexOff(0)	    
  ,fZvertexOff(0)
  ,fEtaOff(0)
  ,fEtaMomentumcutOff(0)
  ,fNclusVSphiOff(0)
  ,fNclusVSthetaOff(0)
  ,fStatusOff(0)
  ,fEventSpecieOff(0)

  ,fChargeHLT(0)      
  ,fMomentumHLT(0)
  ,fMomentumHLTTpc(0)
  ,fMomentumHLTTpcIts(0)
  ,fDCArHLTSG(0)  
  ,fDCAzHLTSG(0)  
  ,fNclusterHLT(0)
  ,fNclusterHLTwCut(0)
  ,fdEdxHLT(0)    
  ,fdEdxVSPHLT(0)
  ,fPhiHLT(0)     
  ,fThetaHLT(0)  
  ,fMultHLT(0)  
  ,fXYvertexHLT(0)
  ,fXvertexHLT(0)
  ,fYvertexHLT(0)
  ,fZvertexHLT(0)
  ,fEtaHLT(0)
  ,fEtaMomentumcutHLT(0)
  ,fNclusVSphiHLT(0)	    
  ,fNclusVSthetaHLT(0)
  ,fStatusHLT(0)
  ,fEventSpecieHLT(0)
  ,fTrgClsArray(0)
  //     ,fDCArOff_trig(0)
  //     ,fNclusterOff_trig(0)
  //     
  //     ,fDCAHLT_trig(0)
  //     ,fNclusterHLT_trig(0)

{
 
  // Constructor
  // Define input and output slots here
  // Input slot #0 works with a TChain
  // DefineInput(0, TChain::Class());
  // Output slot #0 writes into a TH1 container
  DefineOutput(1, TList::Class());
}


//----------------------------------------------------------------------------------------------------

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

  fHistTrigger = new TH1F("fHistTrigger", "CTP trigger counter",24 , 0, 24);
  fHistTrigger->GetXaxis()->SetTitle("");  
  fHistTrigger->GetYaxis()->SetTitle("#Events"); 

  fHistHLTTrigger = new TH1F("fHistHLTTrigger", "HLT CTP trigger counter", 24, 0, 24); 
  fHistHLTTrigger->GetXaxis()->SetTitle("");
  fHistHLTTrigger->GetYaxis()->SetTitle("#Events");

  fChargeOff = new TH1F("fCharge_off", "Charge distribution (Offline)", 12, -3, 3);  
  fChargeHLT = new TH1F("fCharge_hlt", "Charge distribution (HLT)", 12, -3, 3);  
  
  fMomentumOff = new TH1F("fMomentum_off", "momentum (offline)",1000, 0., 100);
  fMomentumHLT = new TH1F("fMomentum_hlt", "momentum (HLT)",    1000, 0., 100);

  fMomentumOffTpc = new TH1F("fMomentumTpc_off","Momentum for kTPCin (offline)",100,-3,3);
  fMomentumHLTTpc = new TH1F("fMomentumTpc_hlt","Momentum for kTPCin (HLT)",    100,-3,3);

  fMomentumOffTpcIts = new TH1F("fMomentumTpcIts_off","Momentum for kTPCin && kITSin (offline)",100,-3,3);
  fMomentumHLTTpcIts = new TH1F("fMomentumTpcIts_hlt","Momentum for kTPCin && kITSin (HLT)",    100,-3,3);
 
  fDCArOff   = new TH1F("fDCA_off",  "DCAr to beam line (offline)",200, -100, 100);
  fDCArHLT   = new TH1F("fDCA_hlt",  "DCAr to beam line (HLT)",    200, -100, 100);
  fDCArHLTSG = new TH1F("fDCA_hltSG","DCAr to beam line (HLT)",    200, -100, 100);

  fDCAzOff   = new TH1F("fDCAz_off",  "DCAz to beam line (offline)",200, -20, 20);
  fDCAzHLT   = new TH1F("fDCAz_hlt",  "DCAz to beam line (HLT)",    200, -20, 20);
  fDCAzHLTSG = new TH1F("fDCAz_hltSG","DCAz to beam line (HLT)",    200, -20, 20);
 
  fNclusterOff = new TH1F("fNcluster_off","clusters per track (offline)", 200, 0, 200);
  fNclusterHLT = new TH1F("fNcluster_hlt","clusters per track (HLT)",     200, 0, 200);
 
  fNclusterOffwCut = new TH1F("fNcluster_wcut_off","clusters per track with cuts (offline)", 200, 0, 200);
  fNclusterHLTwCut = new TH1F("fNcluster_wcut_hlt","clusters per track with cuts (HLT)",     200, 0, 200);
 
  fdEdxOff = new TH1F("fdEdx_off","energy loss (offline)",             500, 0, 500);
  fdEdxHLT = new TH1F("fdEdx_hlt","energy loss (HLT) - not filled yet",500, 0, 500);
 
  fdEdxVSPOff = new TH2F("fdEdx_vs_P_off","dE/dx vs. momentum (offline)",             300, 0., 3., 500, 0., 500.);
  fdEdxVSPHLT = new TH2F("fdEdx_vs_P_hlt","dE/dx vs. momentum (HLT) - not filled yet",300, 0., 3., 500, 0., 500.);

  fPhiOff = new TH1F("fPhi_off","azimuthal angle distribution (offline)",90,0,360);
  fPhiHLT = new TH1F("fPhi_hlt","azimuthal angle distribution (HLT)",    90,0,360);
  
  fThetaOff = new TH1F("fTheta_off","polar angle distribution (offline)",180,0,180);
  fThetaHLT = new TH1F("fTheta_hlt","polar angle distribution (HLT)",    180,0,180);
  
  fMultOff = new TH1F("fMult_off","track multiplicity (offline)",150,0,15000);
  fMultHLT = new TH1F("fMult_hlt","track multiplicity (HLT)",    150,0,15000);
 
  fXYvertexOff = new TH2F("fXYvertex_off","XY primary vertex (offline)",100,-5,5,100,-5,5);
  fXYvertexHLT = new TH2F("fXYvertex_hlt","XY primary vertex (HLT)",    100,-5,5,100,-5,5);
  
  fXvertexOff = new TH1F("fXvertex_off","X primary vertex (offline)",500,-5,5);
  fXvertexHLT = new TH1F("fXvertex_hlt","X primary vertex (HLT)",    500,-5,5);
 
  fYvertexOff = new TH1F("fYvertex_off","Y primary vertex (offline)",500,-5,5);
  fYvertexHLT = new TH1F("fYvertex_hlt","Y primary vertex (HLT)",    500,-5,5);
 
  fZvertexOff = new TH1F("fZvertex_off","Z primary vertex (offline)",250,-30,30);
  fZvertexHLT = new TH1F("fZvertex_hlt","Z primary vertex (HLT)",    250,-30,30);
  
  fEtaOff = new TH1F("fEta_off","pseudorapidity (offline)",100,-2,2);
  fEtaHLT = new TH1F("fEta_hlt","pseudorapidity (HLT)",    100,-2,2);
 
  fEtaMomentumcutOff = new TH1F("fEtaMomentumcut_off","pseudorapidity DCAcut (offline)",100,-2,2);
  fEtaMomentumcutHLT = new TH1F("fEtaMomentumcut_hlt","pseudorapidity DCAcut (HLT)",    100,-2,2);

  fNclusVSphiOff = new TH2F("fNclus_vs_phi_off","clusters per track vs. #phi (offline)",360,0,360,160,0,160);
  fNclusVSphiHLT = new TH2F("fNclus_vs_phi_hlt","clusters per track vs. #phi (HLT)",    360,0,360,160,0,160);
  
  fNclusVSthetaOff = new TH2F("fNclus_vs_theta_off","clusters per track vs. #theta (offline)",180,0,180,160,0,160);
  fNclusVSthetaHLT = new TH2F("fNclus_vs_theta_hlt","clusters per track vs. #theta (HLT)",    180,0,180,160,0,160);

  fStatusOff = new TH1F("fStatus_off", "Status for different detectors (offline)",12, 0, 12);
  fStatusHLT = new TH1F("fStatus_hlt", "Status for different detectors (HLT)",    12, 0, 12);
  
  fEventSpecieOff = new TH1F("fEventSpecie_off","Eventspecie for Offline",18, 0, 18);
  fEventSpecieHLT = new TH1F("fEventSpecie_hlt","Eventspecie for HLT",18, 0, 18);


  //---------------------- add histograms to the output TList ------------------//

  fOutputList->Add(fHistTrigger);
  fOutputList->Add(fHistHLTTrigger);

  fOutputList->Add(fChargeOff);
  fOutputList->Add(fMomentumOff);
  fOutputList->Add(fMomentumOffTpc);
  fOutputList->Add(fMomentumOffTpcIts);
  fOutputList->Add(fDCArOff);	  
  fOutputList->Add(fDCAzOff);	  
  fOutputList->Add(fNclusterOff);
  fOutputList->Add(fNclusterOffwCut);
  fOutputList->Add(fdEdxOff);	  
  fOutputList->Add(fdEdxVSPOff);
  fOutputList->Add(fPhiOff);	  
  fOutputList->Add(fThetaOff);    
  fOutputList->Add(fMultOff);
  fOutputList->Add(fXYvertexOff); 
  fOutputList->Add(fXvertexOff);  
  fOutputList->Add(fYvertexOff);  
  fOutputList->Add(fZvertexOff);  
  fOutputList->Add(fEtaOff);  
  fOutputList->Add(fEtaMomentumcutOff);
  fOutputList->Add(fNclusVSphiOff);  
  fOutputList->Add(fNclusVSthetaOff);
  fOutputList->Add(fStatusOff);
  fOutputList->Add(fEventSpecieOff);

  fOutputList->Add(fChargeHLT);  
  fOutputList->Add(fMomentumHLT); 
  fOutputList->Add(fMomentumHLTTpc);  
  fOutputList->Add(fMomentumHLTTpcIts);
  fOutputList->Add(fDCArHLT);	  
  fOutputList->Add(fDCAzHLT);	  
  fOutputList->Add(fDCArHLTSG);	  
  fOutputList->Add(fDCAzHLTSG);	  
  fOutputList->Add(fNclusterHLT); 
  fOutputList->Add(fNclusterHLTwCut); 
  fOutputList->Add(fdEdxHLT);	  
  fOutputList->Add(fdEdxVSPHLT);
  fOutputList->Add(fPhiHLT);	  
  fOutputList->Add(fThetaHLT);    
  fOutputList->Add(fMultHLT);	 
  fOutputList->Add(fXYvertexHLT); 
  fOutputList->Add(fXvertexHLT);  
  fOutputList->Add(fYvertexHLT);  
  fOutputList->Add(fZvertexHLT);    
  fOutputList->Add(fEtaHLT);  
  fOutputList->Add(fEtaMomentumcutHLT);    
  fOutputList->Add(fNclusVSphiHLT);  
  fOutputList->Add(fNclusVSthetaHLT);
  fOutputList->Add(fStatusHLT);
  fOutputList->Add(fEventSpecieHLT);

  SetupESDtrackCuts();
  
  PostData(1, fOutputList);
}

void AliAnalysisTaskHLT::NotifyRun(){
  // This will not work if the active trigger classes change from run to run.
  // Then one has to know all trigger classes before processing the data.

  AliESDEvent *esdOFF = dynamic_cast<AliESDEvent*>(InputEvent());
  TString trgClasses = esdOFF->GetESDRun()->GetActiveTriggerClasses(); 
 
  fTrgClsArray = trgClasses.Tokenize(" ");
     
  for(Int_t i=0; i<fTrgClsArray->GetEntries(); i++){ 
    TString str = ((TObjString *)fTrgClsArray->At(i))->GetString(); 
    (fHistTrigger->GetXaxis())->SetBinLabel(i+1, str.Data()); 
    (fHistHLTTrigger->GetXaxis())->SetBinLabel(i+1, str.Data()); 
  }   
  esdOFF = NULL;

  TString Statusnames[12]={"kTPCin",
			   "kITSin",
			   "kTPCout",
			   "kITSout",
			   "kITSrefit",
			   "kTPCrefit",
			   "kTRDin",
			   "kTRDout",
			   "kTRDrefit",
			   "kTOFin",
			   "kTOFout",
			   "kTOFrefit"};

  for(int iii=0;iii<12;iii++){
    (fStatusHLT->GetXaxis())->SetBinLabel(iii+1,Statusnames[iii]);
    (fStatusOff->GetXaxis())->SetBinLabel(iii+1,Statusnames[iii]);
  }
}

void AliAnalysisTaskHLT::UserExec(Option_t *){
  // see header file of AliAnalysisTask for documentation

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
    if((esdOFF->GetFiredTriggerClasses()).Contains(((TObjString *)fTrgClsArray->At(i))->GetString()))// && esdOFF->GetEventSpecie()==16)  
      fHistTrigger->Fill(i);
    if((esdHLT->GetFiredTriggerClasses()).Contains(((TObjString *)fTrgClsArray->At(i))->GetString()))// && esdOFF->GetEventSpecie()==16)  
      fHistHLTTrigger->Fill(i);
  }  

  //reject laser events -> Use offline info (HLT do not set this flag)
  if((esdOFF->GetEventSpecie()==16)) return;
  
  // To check if anything was triggered in HLT.
  if(fUseHLTTrigger){  
    if (!((AliHLTGlobalTriggerDecision*)esdHLT->GetHLTTriggerDecision())->Result()){
      return;
      // Process all and any events that were triggered by HLT.
    } 
  }


  Double_t bfield = esdOFF->GetMagneticField();
 
  UInt_t Statusnames[12]={AliESDtrack::kTPCin,
			  AliESDtrack::kITSin,
			  AliESDtrack::kTPCout,
			  AliESDtrack::kITSout,
			  AliESDtrack::kITSrefit,
			  AliESDtrack::kTPCrefit,
			  AliESDtrack::kTRDin,
			  AliESDtrack::kTRDout,
			  AliESDtrack::kTRDrefit,
			  AliESDtrack::kTOFin,
			  AliESDtrack::kTOFout,
			  AliESDtrack::kTOFrefit};
  


  //---------------- HLT ESD tree -----------------------//

  Int_t nr_tracksHLT=0;	 
  Double_t vertexHLT[3];

  const AliESDVertex *vertHLT=esdHLT->GetPrimaryVertexTracks();

  vertexHLT[0] = vertHLT->GetX();
  vertexHLT[1] = vertHLT->GetY();
  vertexHLT[2] = vertHLT->GetZ();
  AliVertex *primVertexHLT = new AliVertex(vertexHLT, 0., 0);

  Bool_t testVertexHLT=kTRUE;
  Int_t nr_contributorsHLT= vertHLT->GetNContributors();

  if(nr_contributorsHLT<1) {
    // SPD vertex
    vertHLT = esdHLT->GetPrimaryVertexSPD();
    if(nr_contributorsHLT<1) {
      // NO GOOD VERTEX, SKIP EVENT 
      testVertexHLT=kFALSE;
    }
  }
  if(testVertexHLT){
    fXYvertexHLT->Fill(vertHLT->GetX(), vertHLT->GetY() );
    fXvertexHLT->Fill( vertHLT->GetX() );
    fYvertexHLT->Fill( vertHLT->GetY() );
    fZvertexHLT->Fill( vertHLT->GetZ() );
  }
  //At the moment no constrains on vertex before filling histograms
  //Should be changed. 

  //if( vertHLT && vertHLT->GetNContributors() >= 5 && (TMath::Abs(vertHLT->GetZ())<5.5) ){

  fEventSpecieHLT->Fill(esdHLT->GetEventSpecie());

  for(Int_t i=0; i<esdHLT->GetNumberOfTracks(); i++){ 
  
    AliESDtrack *esdtrackHLT = esdHLT->GetTrack(i); 
    if(!esdtrackHLT) continue;

    //Fill which status flags that are set for an event
    for(int jjj=0;jjj<12;jjj++){
      if(esdtrackHLT->GetStatus()&Statusnames[jjj]) fStatusHLT->Fill(jjj);
    } 

    //This condition is mostly affecting Offline->will cut away tracks that are counted twice 
    //With both kITSin and kTPCin flags set.
    if(!(esdtrackHLT->GetStatus()&AliESDtrack::kTPCin)) continue;

    //ESD-cut 
    //At the moment not included!
    //if(!fESDHLTtrackCuts->AcceptTrack(esdtrackHLT)) continue; 		    

    //Calculating DCA "old" fashion
    Float_t dca[2];
    esdtrackHLT->GetDZ(esdHLT->GetPrimaryVertex()->GetXv(), esdHLT->GetPrimaryVertex()->GetYv(), esdHLT->GetPrimaryVertex()->GetZv(), bfield, dca);

    // plotting the DCA calculated by Sergey 
    Float_t DCAr, DCAz = -99;
    esdtrackHLT->GetImpactParametersTPC(DCAr, DCAz);
    fDCArHLTSG->Fill(DCAr);
    fDCAzHLTSG->Fill(DCAz);

    fDCArHLT->Fill(dca[0]);  
    fDCAzHLT->Fill(dca[1]);
    
    fChargeHLT->Fill(esdtrackHLT->Charge());

    if((esdtrackHLT->GetStatus()&AliESDtrack::kTPCin) || (esdtrackHLT->GetStatus()&AliESDtrack::kTPCin && esdtrackHLT->GetStatus()&AliESDtrack::kTPCout))
      fNclusterHLT->Fill(esdtrackHLT->GetTPCNcls());

    nr_tracksHLT++;
    
    if((esdtrackHLT->GetStatus()&AliESDtrack::kTPCin) || (esdtrackHLT->GetStatus()&AliESDtrack::kTPCin && esdtrackHLT->GetStatus()&AliESDtrack::kTPCout)){
      fNclusVSphiHLT->Fill(esdtrackHLT->Phi()*TMath::RadToDeg(), esdtrackHLT->GetTPCNcls());
      fNclusVSthetaHLT->Fill(esdtrackHLT->Theta()*TMath::RadToDeg(), esdtrackHLT->GetTPCNcls());
    }

    if(esdtrackHLT->GetStatus()&AliESDtrack::kTPCin) fMomentumHLTTpc->Fill(esdtrackHLT->Pt());
    if(esdtrackHLT->GetStatus()&AliESDtrack::kTPCin && esdtrackHLT->GetStatus()&AliESDtrack::kITSin)   fMomentumHLTTpcIts->Fill(esdtrackHLT->Pt());   

    fEtaHLT->Fill(esdtrackHLT->Eta()); 
    fdEdxHLT->Fill( esdtrackHLT->GetTPCsignal() );
    fdEdxVSPHLT->Fill( TMath::Abs(esdtrackHLT->Pt()), esdtrackHLT->GetTPCsignal() );	     

    fPhiHLT->Fill(esdtrackHLT->Phi()*TMath::RadToDeg());
    fThetaHLT->Fill(esdtrackHLT->Theta()*TMath::RadToDeg());
    fMomentumHLT->Fill( TMath::Abs(esdtrackHLT->Pt()) );  
   
    if(esdtrackHLT->GetStatus()&AliESDtrack::kTPCin) 
      fMomentumHLTTpc->Fill(esdtrackHLT->Pt());
    if(esdtrackHLT->GetStatus()&AliESDtrack::kTPCin && esdtrackHLT->GetStatus()&AliESDtrack::kITSin)   
      fMomentumHLTTpcIts->Fill(esdtrackHLT->Pt());  

    if(TMath::Abs(esdtrackHLT->Pt()) <0.3) continue; 
    fEtaMomentumcutHLT->Fill(esdtrackHLT->Eta());
    if((esdtrackHLT->GetStatus()&AliESDtrack::kTPCin) || (esdtrackHLT->GetStatus()&AliESDtrack::kTPCin && esdtrackHLT->GetStatus()&AliESDtrack::kTPCout))
      fNclusterHLTwCut->Fill(esdtrackHLT->GetTPCNcls());
    
    if(esdHLT->IsHLTTriggerFired()){
        	 
    }// end if for triggered hlt events        
  } // end of loop over hlt tracks

  if(nr_tracksHLT>0) fMultHLT->Fill(nr_tracksHLT);

  //----------------- OFFLINE ESD tree ----------------//
  
  Int_t nr_tracksOff=0;
  Double_t vertexOFF[3];

  const AliESDVertex *vertOff=esdOFF->GetPrimaryVertexTracks();

  vertexOFF[0] = vertOff->GetX();
  vertexOFF[1] = vertOff->GetY();
  vertexOFF[2] = vertOff->GetZ();
  AliVertex *primVertexOFF = new AliVertex(vertexOFF, 0., 0);
  Bool_t testVertex=kTRUE;


  if(vertOff->GetNContributors()<1) {
    // SPD vertex
    vertOff = esdOFF->GetPrimaryVertexSPD();
    if(vertOff->GetNContributors()<1) {
      // NO GOOD VERTEX, SKIP EVENT 
      testVertex=kFALSE;
    }
  }
  
  if(testVertex){
    fXYvertexOff->Fill(vertOff->GetX(), vertOff->GetY() );
    fXvertexOff->Fill( vertOff->GetX() );
    fYvertexOff->Fill( vertOff->GetY() );
    fZvertexOff->Fill( vertOff->GetZ() );
  }

  //At the moment no constrains on vertex before filling histograms
  //Should be changed. 
  fEventSpecieOff->Fill(esdOFF->GetEventSpecie());

  for(Int_t i=0; i<esdOFF->GetNumberOfTracks(); i++){ 
   
    AliESDtrack *esdtrackOFF = esdOFF->GetTrack(i); 
    if (!esdtrackOFF) continue;

    //Fill histograms with which flags are set
    for(int jjj=0;jjj<12;jjj++){
      if(esdtrackOFF->GetStatus()&Statusnames[jjj]) fStatusOff->Fill(jjj);
    } 

    //This condition is mostly affecting Offline->will cut away tracks that are counted twice 
    //With both kITSin and kTPCin flags set.
    if(!(esdtrackOFF->GetStatus()&AliESDtrack::kTPCin))continue; 

    // -- ESD cuts
    //Not included at the moment
    //if(!fESDOfftrackCuts->AcceptTrack(esdtrackOFF) ) continue;
    nr_tracksOff++;

    //Filling of DCA for Offline

    Double_t x[3]; 
    esdtrackOFF->GetXYZ(x);
    Double_t b[3]; 
    AliTracker::GetBxByBz(x,b);
    Bool_t isOK = esdtrackOFF->RelateToVertexTPCBxByBz(vertOff, b, kVeryBig);
    if(!isOK) return;
    
    const AliExternalTrackParam *track = esdtrackOFF->GetTPCInnerParam();
    if(!track) return;
    
    Float_t dca[2], cov[3]; // dca_xy, dca_z, sigma_xy, sigma_xy_z, sigma_z
    esdtrackOFF->GetImpactParametersTPC(dca,cov);

    fDCArOff->Fill(dca[0]);
    fDCAzOff->Fill(dca[1]);
    fChargeOff->Fill(esdtrackOFF->Charge());

    if((esdtrackOFF->GetStatus()&AliESDtrack::kTPCin) || (esdtrackOFF->GetStatus()&AliESDtrack::kTPCin && esdtrackOFF->GetStatus()&AliESDtrack::kTPCout))
      fNclusterOff->Fill(esdtrackOFF->GetTPCNcls()); 

    if(esdtrackOFF->GetTPCNcls()>0) fNclusVSphiOff->Fill(esdtrackOFF->Phi()*TMath::RadToDeg(), esdtrackOFF->GetTPCNcls());
    if(esdtrackOFF->GetTPCNcls()>0) fNclusVSthetaOff->Fill(esdtrackOFF->Theta()*TMath::RadToDeg(), esdtrackOFF->GetTPCNcls());

    if(esdtrackOFF->GetStatus()&AliESDtrack::kTPCin) fMomentumOffTpc->Fill(esdtrackOFF->Pt());
    if(esdtrackOFF->GetStatus()&AliESDtrack::kTPCin && esdtrackOFF->GetStatus()&AliESDtrack::kITSin)   fMomentumOffTpcIts->Fill(esdtrackOFF->Pt());   

    fEtaOff->Fill(esdtrackOFF->Eta());  	
    fdEdxOff->Fill( esdtrackOFF->GetTPCsignal() );
    fdEdxVSPOff->Fill( TMath::Abs(esdtrackOFF->Pt()), esdtrackOFF->GetTPCsignal() );	     
    
    fPhiOff->Fill(esdtrackOFF->Phi()*TMath::RadToDeg());
    fThetaOff->Fill(esdtrackOFF->Theta()*TMath::RadToDeg());
    fMomentumOff->Fill( TMath::Abs(esdtrackOFF->Pt()) ); 

    if(esdtrackOFF->GetStatus()&AliESDtrack::kTPCin) 
      fMomentumOffTpc->Fill(esdtrackOFF->Pt());
    if(esdtrackOFF->GetStatus()&AliESDtrack::kTPCin && esdtrackOFF->GetStatus()&AliESDtrack::kITSin)
      fMomentumOffTpcIts->Fill(esdtrackOFF->Pt());  

    if(TMath::Abs(esdtrackOFF->Pt()) < 0.3) continue;
    fEtaMomentumcutOff->Fill(esdtrackOFF->Eta()); 
    if(esdtrackOFF->GetTPCNcls()>0) fNclusterOffwCut->Fill(esdtrackOFF->GetTPCNcls()); 
                  
//     if(esdHLT->IsHLTTriggerFired()){
//           
//     } // end if for triggered offl events

  } // end of loop over offline tracks
  
  if(nr_tracksOff>0) fMultOff->Fill(nr_tracksOff);
  
  delete primVertexOFF;
  delete primVertexHLT;
  
  // Post output data.
  PostData(1, fOutputList);
}

void AliAnalysisTaskHLT::Terminate(Option_t *){
  // see header file of AliAnalysisTask for documentation
}

void AliAnalysisTaskHLT::SetupESDtrackCuts() {
  // Setup ESD cuts
  // NB! Work in progress!

  Bool_t selPrimaries = kTRUE;
  
  fESDOfftrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2009(selPrimaries);
  //To make Offline cuts = HLT cuts
  fESDOfftrackCuts->SetRequireITSRefit(kFALSE); 
  fESDOfftrackCuts->SetEtaRange(-0.9,+0.9);
  

  //HLT
  //NB! Need to understand this a bit more! Which cuts should we keep?
  fESDHLTtrackCuts = new AliESDtrackCuts;

  // TPC  
  fESDHLTtrackCuts->SetRequireTPCStandAlone(kTRUE); // to get chi2 and ncls of kTPCin
  fESDHLTtrackCuts->SetMinNClustersTPC(70);
  fESDHLTtrackCuts->SetMaxChi2PerClusterTPC(4);
  fESDHLTtrackCuts->SetAcceptKinkDaughters(kFALSE);

  fESDHLTtrackCuts->SetClusterRequirementITS(AliESDtrackCuts::kSPD,
  				     AliESDtrackCuts::kAny);

  if(selPrimaries) { // 7*(0.0050+0.0060/pt^0.9)
    fESDHLTtrackCuts->SetMaxDCAToVertexXYPtDep("0.0350+0.0420/pt^0.9");
  }
  
  fESDHLTtrackCuts->SetMaxDCAToVertexZ(1.e6);
  fESDHLTtrackCuts->SetDCAToVertex2D(kFALSE);
  fESDHLTtrackCuts->SetRequireSigmaToVertex(kFALSE);
  fESDHLTtrackCuts->SetEtaRange(-0.9,+0.9);

  return;
}
