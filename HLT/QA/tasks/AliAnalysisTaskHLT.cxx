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

#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TString.h"
#include "AliESDEvent.h"
#include "AliESDtrackCuts.h"
#include "AliESDInputHandler.h"

#include "AliAnalysisTaskHLT.h"

ClassImp(AliAnalysisTaskHLT)

//======================================================================================================

AliAnalysisTaskHLT::AliAnalysisTaskHLT()
:
AliAnalysisTaskSE()
  ,fESDOfftrackCuts(0)
  ,fESDHLTtrackCuts(0)
  ,fOutputList(0)
  ,fHistTrigger(0)
  ,fHistHLTTrigger(0)  
  ,fChargeOff(0)  
  ,fMomentumOff(0)	
  ,fDCAOff(0)  	
  ,fNclusterOff(0)
  ,fNclusterOffwCut(0)		
  ,fdEdxOff(0) 	
  ,fdEdxVSPOff(0)	
  ,fPhiOff(0)  	
  ,fThetaOff(0)	
  ,fMultOff(0) 	
  ,fVertexVSNtracksOff(0)
  ,fXYvertexOff(0)	
  ,fXvertexOff(0)	    
  ,fYvertexOff(0)	    
  ,fZvertexOff(0)
  ,fEtaOff(0)
  ,fEtaDCAcutOff(0)
  ,fEtaOffTpc(0)
  ,fEtaOffTpcIts(0)
  ,fNclusVSphiOff(0)
  ,fNclusVSthetaOff(0)
  ,fStatusOff(0)
  ,fStatusOff_Ocl(0)
  ,fEventSpecieOff(0)
  
  ,fChargeHLT(0)
  ,fMomentumHLT(0)
  ,fDCAHLT(0)  
  ,fNclusterHLT(0)
  ,fNclusterHLTwCut(0)
  ,fdEdxHLT(0)    
  ,fdEdxVSPHLT(0)
  ,fPhiHLT(0)     
  ,fThetaHLT(0)  
  ,fMultHLT(0)  
  ,fVertexVSNtracksHLT(0) 
  ,fXYvertexHLT(0)
  ,fXvertexHLT(0)
  ,fYvertexHLT(0)
  ,fZvertexHLT(0)
  ,fEtaHLT(0)
  ,fEtaDCAcutHLT(0)
  ,fEtaHLTTpc(0)
  ,fEtaHLTTpcIts(0)
  ,fNclusVSphiHLT(0)	    
  ,fNclusVSthetaHLT(0)
  ,fStatusHLT(0)
  ,fStatusHLT_Ocl(0)
  ,fEventSpecieHLT(0)
  
  //     ,fDCAOff_trig(0)
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

  //DefineOutput(1, TList::Class());
}
 
AliAnalysisTaskHLT::AliAnalysisTaskHLT(const char *name)
  :
  AliAnalysisTaskSE(name)    
  ,fOutputList(0)
  ,fHistTrigger(0)
  ,fHistHLTTrigger(0)  
  ,fChargeOff(0)  
  ,fMomentumOff(0)	
  ,fDCAOff(0)  	
  ,fNclusterOff(0)
  ,fNclusterOffwCut(0)	
  ,fdEdxOff(0) 	
  ,fdEdxVSPOff(0)	
  ,fPhiOff(0)  	
  ,fThetaOff(0)	
  ,fMultOff(0) 	
  ,fVertexVSNtracksOff(0)
  ,fXYvertexOff(0)	
  ,fXvertexOff(0)	    
  ,fYvertexOff(0)	    
  ,fZvertexOff(0)
  ,fEtaOff(0)
  ,fEtaDCAcutOff(0)
  ,fEtaOffTpc(0)
  ,fEtaOffTpcIts(0)
  ,fNclusVSphiOff(0)
  ,fNclusVSthetaOff(0)
  ,fStatusOff(0)
  ,fStatusOff_Ocl(0)
  ,fEventSpecieOff(0)

  ,fChargeHLT(0)      
  ,fMomentumHLT(0)
  ,fDCAHLT(0)  
  ,fNclusterHLT(0)
  ,fNclusterHLTwCut(0)
  ,fdEdxHLT(0)    
  ,fdEdxVSPHLT(0)
  ,fPhiHLT(0)     
  ,fThetaHLT(0)  
  ,fMultHLT(0)  
  ,fVertexVSNtracksHLT(0) 
  ,fXYvertexHLT(0)
  ,fXvertexHLT(0)
  ,fYvertexHLT(0)
  ,fZvertexHLT(0)
  ,fEtaHLT(0)
  ,fEtaDCAcutHLT(0)
  ,fEtaHLTTpc(0)
  ,fEtaHLTTpcIts(0)
  ,fNclusVSphiHLT(0)	    
  ,fNclusVSthetaHLT(0)
  ,fStatusHLT(0)
  ,fStatusHLT_Ocl(0)    
  ,fEventSpecieHLT(0)
  //     ,fDCAOff_trig(0)
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

// const Float_t AliAnalysisTaskHLT::fgkEtaMin = -0.12;  
// const Float_t AliAnalysisTaskHLT::fgkEtaMax =  0.12;  
// const Float_t AliAnalysisTaskHLT::fgkPhiMin[5]   = {3.83972, 4.18879, 4.53786, 4.88692, 5.23599};  
// const Float_t AliAnalysisTaskHLT::fgkPhiMax[5]   = {4.18879, 4.53786, 4.88692, 5.23599, 5.58505};  
// const Float_t AliAnalysisTaskHLT::fgkNormX[5]    = {-0.642788, -0.34202, 0, 0.34202, 0.642788};  
// const Float_t AliAnalysisTaskHLT::fgkNormY[5]    = {-0.766044, -0.939693, -1, -0.939693, -0.766044};  
// const Float_t AliAnalysisTaskHLT::fgkInitPosX[5] = {-295.682, -157.329, 0, 157.329, 295.682};  
// const Float_t AliAnalysisTaskHLT::fgkInitPosY[5] = {-352.38, -432.259, -460, -432.259, -352.38};

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

  fChargeOff = new TH1F("fCharge_off", "Charge distribution (Offline)", 12, -3, 3);  
  fChargeHLT = new TH1F("fCharge_hlt", "Charge distribution (HLT)", 12, -3, 3);  

  
  fMomentumOff = new TH1F("fMomentum_off", "momentum (offline)",1000, 0., 100);
  fMomentumHLT = new TH1F("fMomentum_hlt", "momentum (HLT)",    1000, 0., 100);
 
  fDCAOff = new TH1F("fDCA_off","DCA to beam line (offline)",200, -100, 100);
  fDCAHLT = new TH1F("fDCA_hlt","DCA to beam line (HLT)",    200, -100, 100);
 
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
  
  fMultOff = new TH1F("fMult_off","track multiplicity (offline)",100,0,100);
  fMultHLT = new TH1F("fMult_hlt","track multiplicity (HLT)",    100,0,100);

  fVertexVSNtracksOff = new TH2F("fVertexNtracs_off", "Vertex resolution vs nr contributing tracks (Offline)",20, 0, 20, 40, -2, 2); 
  fVertexVSNtracksHLT =  new TH2F("fVertexNtracs_hlt", "Vertex resolution vs nr contributing tracks (HLT)",20, 0, 20, 40, -2, 2); 
  
  fXYvertexOff = new TH2F("fXYvertex_off","XY primary vertex (offline)",100,-5,5,100,-5,5);
  fXYvertexHLT = new TH2F("fXYvertex_hlt","XY primary vertex (HLT)",    100,-5,5,100,-5,5);
  
  fXvertexOff = new TH1F("fXvertex_off","X primary vertex (offline)",1000,-1,1);
  fXvertexHLT = new TH1F("fXvertex_hlt","X primary vertex (HLT)",    1000,-1,1);
 
  fYvertexOff = new TH1F("fYvertex_off","Y primary vertex (offline)",1000,-1,1);
  fYvertexHLT = new TH1F("fYvertex_hlt","Y primary vertex (HLT)",    1000,-1,1);
 
  fZvertexOff = new TH1F("fZvertex_off","Z primary vertex (offline)",250,-30,30);
  fZvertexHLT = new TH1F("fZvertex_hlt","Z primary vertex (HLT)",    250,-30,30);
  
  fEtaOff = new TH1F("fEta_off","pseudorapidity (offline)",100,-3,3);
  fEtaHLT = new TH1F("fEta_hlt","pseudorapidity (HLT)",    100,-3,3);
 
  fEtaDCAcutOff = new TH1F("fEtaDCAcut_off","pseudorapidity DCAcut (offline)",100,-3,3);
  fEtaDCAcutHLT = new TH1F("fEtaDCAcut_hlt","pseudorapidity DCAcut (HLT)",    100,-3,3);
 
  fEtaOffTpc = new TH1F("fEtaTpc_off","pseudorapidity for kTPCin (offline)",100,-3,3);
  fEtaHLTTpc = new TH1F("fEtaTpc_hlt","pseudorapidity for kTPCin (HLT)",    100,-3,3);

  fEtaOffTpcIts = new TH1F("fEtaTpcIts_off","pseudorapidity for kTPCin && kITSin (offline)",100,-3,3);
  fEtaHLTTpcIts = new TH1F("fEtaTpcIts_hlt","pseudorapidity for kTPCin && kITSin (HLT)",    100,-3,3);

  fNclusVSphiOff = new TH2F("fNclus_vs_phi_off","clusters per track vs. #phi (offline)",360,0,360,160,0,160);
  fNclusVSphiHLT = new TH2F("fNclus_vs_phi_hlt","clusters per track vs. #phi (HLT)",    360,0,360,160,0,160);
  
  fNclusVSthetaOff = new TH2F("fNclus_vs_theta_off","clusters per track vs. #theta (offline)",180,0,180,160,0,160);
  fNclusVSthetaHLT = new TH2F("fNclus_vs_theta_hlt","clusters per track vs. #theta (HLT)",    180,0,180,160,0,160);

  fStatusOff = new TH1F("fStatus_off", "Status for different detectors (offline)",12, 0, 12);
  fStatusHLT = new TH1F("fStatus_hlt", "Status for different detectors (HLT)",    12, 0, 12);
  
  fStatusOff_Ocl = new TH1F("fStatus_Ocl_off", "Status for different detectors when #TPCcl=0 (offline)",12, 0, 12);
  fStatusHLT_Ocl = new TH1F("fStatus_Ocl_hlt", "Status for different detectors when #TPCcl=0 (HLT)",    12, 0, 12);

  fEventSpecieOff = new TH1F("fEventSpecie_off","Eventspecie for Offline",18, 0, 18);
  fEventSpecieHLT = new TH1F("fEventSpecie_hlt","Eventspecie for HLT",18, 0, 18);

  //---------------------- add histograms to the output TList ------------------//

  fOutputList->Add(fHistTrigger);
  fOutputList->Add(fHistHLTTrigger);

  fOutputList->Add(fChargeOff);
  fOutputList->Add(fMomentumOff);
  fOutputList->Add(fDCAOff);	  
  fOutputList->Add(fNclusterOff);
  fOutputList->Add(fNclusterOffwCut);
  fOutputList->Add(fdEdxOff);	  
  fOutputList->Add(fdEdxVSPOff);
  fOutputList->Add(fPhiOff);	  
  fOutputList->Add(fThetaOff);    
  fOutputList->Add(fMultOff);
  fOutputList->Add(fVertexVSNtracksOff);
  fOutputList->Add(fXYvertexOff); 
  fOutputList->Add(fXvertexOff);  
  fOutputList->Add(fYvertexOff);  
  fOutputList->Add(fZvertexOff);  
  fOutputList->Add(fEtaOff);  
  fOutputList->Add(fEtaDCAcutOff);
  fOutputList->Add(fEtaOffTpc);
  fOutputList->Add(fEtaOffTpcIts);
  fOutputList->Add(fNclusVSphiOff);  
  fOutputList->Add(fNclusVSthetaOff);
  fOutputList->Add(fStatusOff);
  fOutputList->Add(fStatusOff_Ocl);
  fOutputList->Add(fEventSpecieOff);

  fOutputList->Add(fChargeHLT);  
  fOutputList->Add(fMomentumHLT); 
  fOutputList->Add(fDCAHLT);	  
  fOutputList->Add(fNclusterHLT); 
  fOutputList->Add(fNclusterHLTwCut); 
  fOutputList->Add(fdEdxHLT);	  
  fOutputList->Add(fdEdxVSPHLT);
  fOutputList->Add(fPhiHLT);	  
  fOutputList->Add(fThetaHLT);    
  fOutputList->Add(fMultHLT);	 
  fOutputList->Add(fVertexVSNtracksHLT); 
  fOutputList->Add(fXYvertexHLT); 
  fOutputList->Add(fXvertexHLT);  
  fOutputList->Add(fYvertexHLT);  
  fOutputList->Add(fZvertexHLT);    
  fOutputList->Add(fEtaHLT);  
  fOutputList->Add(fEtaDCAcutHLT);
  fOutputList->Add(fEtaHLTTpc);  
  fOutputList->Add(fEtaHLTTpcIts);    
  fOutputList->Add(fNclusVSphiHLT);  
  fOutputList->Add(fNclusVSthetaHLT);
  fOutputList->Add(fStatusHLT);
  fOutputList->Add(fStatusHLT_Ocl);
  fOutputList->Add(fEventSpecieHLT);

  SetupESDtrackCuts();
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
    (fStatusHLT_Ocl->GetXaxis())->SetBinLabel(iii+1,Statusnames[iii]);
    (fStatusOff_Ocl->GetXaxis())->SetBinLabel(iii+1,Statusnames[iii]);
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

  Double_t DCAcut = 7.0;
  Double_t Momcut= 0.3;

  char test[50];
  sprintf(test,"pseudorapidity (HLT), DCA cut = %f,\n Momentum cut = %f, TPC clusters > 70" , DCAcut, Momcut);
  fEtaDCAcutHLT->SetTitle(test);
  sprintf(test,"pseudorapidity (offline), DCA cut = %f, Momentum cut = %f, TPC clusters > 70", DCAcut, Momcut);
  fEtaDCAcutOff->SetTitle(test);

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

  if(esdHLT->GetNumberOfTracks()!=0)
    fMultHLT->Fill( esdHLT->GetNumberOfTracks() );
	 
  Double_t vertexHLT[3];

  const AliESDVertex *vertHLT=esdHLT->GetPrimaryVertexTracks();

  vertexHLT[0] = vertHLT->GetX();
  vertexHLT[1] = vertHLT->GetY();
  vertexHLT[2] = vertHLT->GetZ();
  AliVertex *primVertexHLT = new AliVertex(vertexHLT, 0., 0);

  Bool_t testVertexHLT=kTRUE;
  if(vertHLT->GetNContributors()<1) {
    // SPD vertex
    vertHLT = esdHLT->GetPrimaryVertexSPD();
    if(vertHLT->GetNContributors()<1) {
      // NO GOOD VERTEX, SKIP EVENT 
      testVertexHLT=kFALSE;
    }
  }
  fVertexVSNtracksHLT->Fill(vertHLT->GetNContributors(),vertHLT->GetZv());
  if(vertHLT->GetZ()!=0){
    fXYvertexHLT->Fill(vertHLT->GetX(), vertHLT->GetY() );
    fXvertexHLT->Fill( vertHLT->GetX() );
    fYvertexHLT->Fill( vertHLT->GetY() );
    fZvertexHLT->Fill( vertHLT->GetZ() );
  }


  if(testVertexHLT){
    //if( vertHLT && vertHLT->GetNContributors() >= 5 && (TMath::Abs(vertHLT->GetZ())<5.5) ){

    fEventSpecieHLT->Fill(esdHLT->GetEventSpecie());

    for(Int_t i=0; i<esdHLT->GetNumberOfTracks(); i++){ 
  
      AliESDtrack *esdtrackHLT = esdHLT->GetTrack(i); 
      if(!esdtrackHLT){
	continue;
      }

      //Fill which status flags that are set for an event
      for(int jjj=0;jjj<12;jjj++){
	if(esdtrackHLT->GetStatus()&Statusnames[jjj]) {
	  fStatusHLT->Fill(jjj);
	  if(esdtrackHLT->GetTPCNcls()==0) fStatusHLT_Ocl->Fill(jjj);
	}
      } 

      if(esdtrackHLT->GetTPCNcls()>0) fNclusterHLT->Fill(esdtrackHLT->GetTPCNcls());
      fChargeHLT->Fill(esdtrackHLT->Charge());
 
      if(esdtrackHLT->GetStatus()&AliESDtrack::kTPCin) fEtaHLTTpc->Fill(esdtrackHLT->Eta());
      if (esdtrackHLT->GetStatus()&AliESDtrack::kTPCin && esdtrackHLT->GetStatus()&AliESDtrack::kITSin)	 fEtaHLTTpcIts->Fill(esdtrackHLT->Eta());   
      //ESD-cut
      if(!fESDHLTtrackCuts->AcceptTrack(esdtrackHLT) )continue;		
      
      if((esdHLT->GetEventSpecie()==16)) {
	Printf("Reject laser event %d",esdOFF->GetEventSpecie());
	continue; //reject laser events
      }
      
      if(esdtrackHLT->GetTPCNcls()>0) fNclusVSphiHLT->Fill(esdtrackHLT->Phi()*TMath::RadToDeg(), esdtrackHLT->GetTPCNcls());
      if(esdtrackHLT->GetTPCNcls()>0) fNclusVSthetaHLT->Fill(esdtrackHLT->Theta()*TMath::RadToDeg(), esdtrackHLT->GetTPCNcls());
      fEtaHLT->Fill(esdtrackHLT->Eta());
 
      fDCAHLT->Fill(esdtrackHLT->GetD(esdHLT->GetPrimaryVertex()->GetXv(), esdHLT->GetPrimaryVertex()->GetYv(), bfield) ); 
      
      fdEdxHLT->Fill( esdtrackHLT->GetTPCsignal() );
      fdEdxVSPHLT->Fill( TMath::Abs(esdtrackHLT->P()), esdtrackHLT->GetTPCsignal() );         


      if(TMath::Abs(esdtrackHLT->Pt()) <Momcut) continue; //cut away tracks with mom<0.3GeV
      fEtaDCAcutHLT->Fill(esdtrackHLT->Eta());
      fPhiHLT->Fill(esdtrackHLT->Phi()*TMath::RadToDeg());
      fThetaHLT->Fill(esdtrackHLT->Theta()*TMath::RadToDeg());
      if(esdtrackHLT->GetTPCNcls()>0) fNclusterHLTwCut->Fill(esdtrackHLT->GetTPCNcls());
      fMomentumHLT->Fill( TMath::Abs(esdtrackHLT->P()) ); 

	  
      if(esdHLT->IsHLTTriggerFired()){
	    	    
      }// end if for triggered hlt events	 
    } // end of loop over hlt tracks
  }  


  //----------------- OFFLINE ESD tree ----------------//
  
  if(esdOFF->GetNumberOfTracks()!=0)
    fMultOff->Fill( esdOFF->GetNumberOfTracks() );

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
      //      Printf("This vertex is away");
    }
  }
  fVertexVSNtracksOff->Fill(vertOff->GetNContributors(),vertOff->GetZv());
  
  if(vertOff->GetZ()!=0){
    fXYvertexOff->Fill(vertOff->GetX(), vertOff->GetY() );
    fXvertexOff->Fill( vertOff->GetX() );
    fYvertexOff->Fill( vertOff->GetY() );
    fZvertexOff->Fill( vertOff->GetZ() );
  }

  if(testVertex){ 

    fEventSpecieOff->Fill(esdOFF->GetEventSpecie());

    for(Int_t i=0; i<esdOFF->GetNumberOfTracks(); i++){ 
     
      AliESDtrack *esdtrackOFF = esdOFF->GetTrack(i); 

      if (!esdtrackOFF) continue;

      fChargeOff->Fill(esdtrackOFF->Charge());

      if(esdtrackOFF->GetStatus()&AliESDtrack::kTPCin) fEtaOffTpc->Fill(esdtrackOFF->Eta());
      if (esdtrackOFF->GetStatus()&AliESDtrack::kTPCin && esdtrackOFF->GetStatus()&AliESDtrack::kITSin)	 fEtaOffTpcIts->Fill(esdtrackOFF->Eta());   
      
      //Fill histograms with which flags are set
      for(int jjj=0;jjj<12;jjj++){
	if(esdtrackOFF->GetStatus()&Statusnames[jjj]) {
	  fStatusOff->Fill(jjj);
	  if(esdtrackOFF->GetTPCNcls()==0) fStatusOff_Ocl->Fill(jjj);
	}
      } 
      if(esdtrackOFF->GetTPCNcls()>0) fNclusterOff->Fill(esdtrackOFF->GetTPCNcls()); 


      if(!fESDOfftrackCuts->AcceptTrack(esdtrackOFF) ) continue;// -- ESD cuts
      if((esdOFF->GetEventSpecie()==16)) continue;        // reject laser events
      if(esdtrackOFF->GetTPCNcls()>0) fNclusVSphiOff->Fill(esdtrackOFF->Phi()*TMath::RadToDeg(), esdtrackOFF->GetTPCNcls());
      if(esdtrackOFF->GetTPCNcls()>0) fNclusVSthetaOff->Fill(esdtrackOFF->Theta()*TMath::RadToDeg(), esdtrackOFF->GetTPCNcls());

      fDCAOff->Fill(esdtrackOFF->GetD(esdOFF->GetPrimaryVertex()->GetXv(), esdOFF->GetPrimaryVertex()->GetYv(), bfield) ); 
      fEtaOff->Fill(esdtrackOFF->Eta());	  
      fdEdxOff->Fill( esdtrackOFF->GetTPCsignal() );
      fdEdxVSPOff->Fill( TMath::Abs(esdtrackOFF->P()), esdtrackOFF->GetTPCsignal() );         
	  
      if(TMath::Abs(esdtrackOFF->Pt()) < Momcut) continue;//cut away tracks with mom<0.3GeV
      fEtaDCAcutOff->Fill(esdtrackOFF->Eta()); 
      fPhiOff->Fill(esdtrackOFF->Phi()*TMath::RadToDeg());
      fThetaOff->Fill(esdtrackOFF->Theta()*TMath::RadToDeg());
      if(esdtrackOFF->GetTPCNcls()>0) fNclusterOffwCut->Fill(esdtrackOFF->GetTPCNcls()); 
      fMomentumOff->Fill( TMath::Abs(esdtrackOFF->P()) ); 
      
	  
      if(esdHLT->IsHLTTriggerFired()){
	    
      } // end if for triggered offl events

    } // end of loop over offl tracks
  
  }


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
  //          fDCAOff->Fill(TMath::Abs(esdTrk->GetD(0., 0., bfield))); 
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
  
  delete primVertexOFF;
  delete primVertexHLT;
  
  // Post output data.
  PostData(1, fOutputList);
}

void AliAnalysisTaskHLT::Terminate(Option_t *){
  // see header file of AliAnalysisTask for documentation

  // Draw result to the screen
  // Called once at the end of the query

  Bool_t print_png=kFALSE;
  if(print_png){

    Int_t maxbin =0;

    TCanvas *c1 = new TCanvas("c1","Info pr track, Offline vs Online",10,10,1210,810);
     
    c1->Divide(3,2);
     
    c1->cd(1);
    maxbin =fEtaOff->GetBinContent(fEtaOff->GetMaximumBin());
    if(maxbin < fEtaHLT->GetBinContent(fEtaHLT->GetMaximumBin()))
      maxbin=fEtaHLT->GetBinContent(fEtaHLT->GetMaximumBin());
    fEtaOff->SetMaximum(maxbin+20);
    fEtaOff->SetTitle("Pseudorapidity (without momentum cut)");
    fEtaOff->SetLineColor(2);
    fEtaOff->DrawCopy("");
    fEtaHLT->DrawCopy("sameS");
 
    TLegend *legend = new TLegend(0.70,0.60,0.90,0.75);
    legend->AddEntry(fEtaOff, "Offline", "LP");
    legend->AddEntry(fEtaHLT,"HLT","LP");
    legend->SetFillColor(10);
    legend->SetBorderSize(0);
    legend->Draw("");

    c1->cd(2);
    maxbin =fEtaDCAcutOff->GetBinContent(fEtaDCAcutOff->GetMaximumBin());
    if(maxbin < fEtaDCAcutHLT->GetBinContent(fEtaDCAcutHLT->GetMaximumBin()))
      maxbin=fEtaDCAcutHLT->GetBinContent(fEtaDCAcutHLT->GetMaximumBin());
    fEtaDCAcutOff->SetMaximum(maxbin+20);
    fEtaDCAcutOff->SetTitle("Pseudorapidity");
    fEtaDCAcutOff->SetLineColor(2);
    fEtaDCAcutOff->DrawCopy("");
    fEtaDCAcutHLT->DrawCopy("sames");

    c1->cd(3);
    maxbin =fNclusterOff->GetBinContent(fNclusterOff->GetMaximumBin());
    if(maxbin < fNclusterHLT->GetBinContent(fNclusterHLT->GetMaximumBin()))
      maxbin=fNclusterHLT->GetBinContent(fNclusterHLT->GetMaximumBin());
    fNclusterOff->SetMaximum(maxbin+20);
    fNclusterOff->SetLineColor(2);
    fNclusterOff->SetTitle("Nr clusters per track");
    fNclusterOff->DrawCopy("");
    fNclusterHLT->DrawCopy("sames");

    c1->cd(4);
    maxbin =fPhiOff->GetBinContent(fPhiOff->GetMaximumBin());
    if(maxbin < fPhiHLT->GetBinContent(fPhiHLT->GetMaximumBin()))
      maxbin=fPhiHLT->GetBinContent(fPhiHLT->GetMaximumBin());
    fPhiOff->SetMinimum(0);
    fPhiOff->SetMaximum(maxbin+20);
    fPhiOff->SetLineColor(2);
    fPhiOff->SetTitle("Azimuthal angle distribution");
    fPhiOff->DrawCopy("");
    fPhiHLT->DrawCopy("sames");

    c1->cd(5);
    maxbin =fThetaOff->GetBinContent(fThetaOff->GetMaximumBin());
    if(maxbin < fThetaHLT->GetBinContent(fThetaHLT->GetMaximumBin()))
      maxbin=fThetaHLT->GetBinContent(fThetaHLT->GetMaximumBin());
    fThetaOff->SetMaximum(maxbin+20);
    fThetaOff->SetLineColor(2);
    fThetaOff->SetTitle("Polar angle distribution");
    fThetaOff->DrawCopy("");
    fThetaHLT->DrawCopy("sames");

    c1->cd(6);
    maxbin =fMomentumOff->GetBinContent(fMomentumOff->GetMaximumBin());
    if(maxbin < fMomentumHLT->GetBinContent(fMomentumHLT->GetMaximumBin()))
      maxbin=fMomentumHLT->GetBinContent(fMomentumHLT->GetMaximumBin());
    fMomentumOff->SetMaximum(maxbin+200);
    fMomentumOff->GetXaxis()->SetRangeUser(0,5);
    fMomentumOff->SetLineColor(2);
    fMomentumOff->SetTitle("Momentum");
    fMomentumOff->DrawCopy("");
    fMomentumHLT->DrawCopy("sames");

    TCanvas *c2= new TCanvas("c2", "Info pr event, Offline vs Online", 10 , 10,1210, 810);
    TLegend *legend2 = new TLegend(0.70,0.60,0.90,0.75);
    c2->Divide(3,2);
    c2->cd(1);
    fXvertexOff->SetTitle("Primary Vertex Distribution in X");
    fXvertexOff->SetLineColor(2);
    fXvertexOff->GetXaxis()->SetRangeUser(-0.5,0.5);
    legend2->AddEntry(fXvertexOff,"Offline","LP");
    fXvertexOff->DrawCopy("");
    fXvertexHLT->DrawCopy("sames");
    legend2->AddEntry(fXvertexHLT,"HLT","LP");
    legend2->SetFillColor(10);
    legend2->SetBorderSize(0);
    legend2->Draw();
    c2->cd(2);
    fYvertexOff->SetTitle("Primary Vertex Distribution in Y");
    fYvertexOff->SetLineColor(2);
    fYvertexOff->GetXaxis()->SetRangeUser(-0.5,0.5);
    fYvertexOff->DrawCopy("");
    fYvertexHLT->DrawCopy("sames");
    c2->cd(3);
    fZvertexOff->SetTitle("Primary Vertex Distribution in Z");
    fZvertexOff->SetLineColor(2);
    fZvertexOff->DrawCopy("");
    fZvertexHLT->DrawCopy("sames");
  
    c2->cd(4);
    fMultHLT->SetTitle("Track Multiplicity, NumberTracks>0");
    fMultHLT->DrawCopy("");
    fMultOff->SetLineColor(2);
    fMultOff->DrawCopy("sames");

    string filename="Info_pr_track";
    char plotname[100];
    sprintf(plotname,"%s.png",filename.c_str());
    //  c1->Print("Info_pr_track.png","png");
    c1->Print(plotname,"png");
 
    filename="Info_pr_Event";
    sprintf(plotname,"%s.png",filename.c_str());
    c2->Print(plotname,"png");
  }

}

void AliAnalysisTaskHLT::SetupESDtrackCuts() {
  // Setup ESD cuts

  Bool_t selPrimaries = kTRUE;

  fESDOfftrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2009(selPrimaries);

  fESDHLTtrackCuts = new AliESDtrackCuts;

  // TPC  
  fESDHLTtrackCuts->SetRequireTPCStandAlone(kTRUE); // to get chi2 and ncls of kTPCin
  fESDHLTtrackCuts->SetMinNClustersTPC(70);
  fESDHLTtrackCuts->SetMaxChi2PerClusterTPC(4);
  fESDHLTtrackCuts->SetAcceptKinkDaughters(kFALSE);

  // -- fESDHLTtrackCuts->SetRequireTPCRefit(kTRUE); -- JMT not present in HLT ESD
  // ITS
  // -- fESDHLTtrackCuts->SetRequireITSRefit(kTRUE); -- JMT not present in HLT ESD

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
