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
#include "TF1.h"
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
  ,fMomentumOffTpc(0)
  ,fMomentumOffTpcIts(0)	
  ,fDCAOff(0)  	
  ,fNclusterOff(0)
  ,fNclusterOffwCut(0)		
  ,fdEdxOff(0) 	
  ,fdEdxVSPOff(0)	
  ,fPhiOff(0)  	
  ,fThetaOff(0)	
  ,fMultOff(0) 	
  ,fXVertexVSNtracksOff(0)
  ,fYVertexVSNtracksOff(0)
  ,fZVertexVSNtracksOff(0)
  ,fZVertexOffTemp(0)
  ,fXYvertexOff(0)	
  ,fXvertexOff(0)	    
  ,fYvertexOff(0)	    
  ,fZvertexOff(0)
  ,fEtaOff(0)
  ,fEtaMomentumcutOff(0)
  ,fNclusVSphiOff(0)
  ,fNclusVSthetaOff(0)
  ,fStatusOff(0)
  ,fStatusOff_Ocl(0)
  ,fEventSpecieOff(0)
  
  ,fChargeHLT(0)
  ,fMomentumHLT(0)
  ,fMomentumHLTTpc(0)
  ,fMomentumHLTTpcIts(0)
  ,fDCAHLT(0)  
  ,fNclusterHLT(0)
  ,fNclusterHLTwCut(0)
  ,fdEdxHLT(0)    
  ,fdEdxVSPHLT(0)
  ,fPhiHLT(0)     
  ,fThetaHLT(0)  
  ,fMultHLT(0)  
  ,fXVertexVSNtracksHLT(0) 
  ,fYVertexVSNtracksHLT(0) 
  ,fZVertexVSNtracksHLT(0) 
  ,fZVertexHLTTemp(0)
  ,fXYvertexHLT(0)
  ,fXvertexHLT(0)
  ,fYvertexHLT(0)
  ,fZvertexHLT(0)
  ,fEtaHLT(0)
  ,fEtaMomentumcutHLT(0)
  ,fNclusVSphiHLT(0)	    
  ,fNclusVSthetaHLT(0)
  ,fStatusHLT(0)
  ,fStatusHLT_Ocl(0)
  ,fEventSpecieHLT(0)
  ,fTrgClsArray(0)
  
  //     ,fDCAOff_trig(0)
  //     ,fNclusterOff_trig(0)
  //     
  //     ,fDCAHLT_trig(0)
  //     ,fNclusterHLT_trig(0)

{
  for(int jj=0;jj<7;jj++){
    fXvertexVSNcontriOff[jj]=0;
    fXvertexVSNcontriHLT[jj]=0;
  }
  for(int jj=0;jj<7;jj++){
    fYvertexVSNcontriOff[jj]=0;
    fYvertexVSNcontriHLT[jj]=0;
  }
  for(int jj=0;jj<7;jj++){
    fZvertexVSNcontriOff[jj]=0;
    fZvertexVSNcontriHLT[jj]=0;
  }




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
  ,fESDOfftrackCuts(0)
  ,fESDHLTtrackCuts(0)
  ,fOutputList(0)
  ,fHistTrigger(0)
  ,fHistHLTTrigger(0)  
  ,fChargeOff(0)  
  ,fMomentumOff(0)
  ,fMomentumOffTpc(0)
  ,fMomentumOffTpcIts(0)	
  ,fDCAOff(0)  	
  ,fNclusterOff(0)
  ,fNclusterOffwCut(0)	
  ,fdEdxOff(0) 	
  ,fdEdxVSPOff(0)	
  ,fPhiOff(0)  	
  ,fThetaOff(0)	
  ,fMultOff(0) 	
  ,fXVertexVSNtracksOff(0)
  ,fYVertexVSNtracksOff(0)
  ,fZVertexVSNtracksOff(0)
  ,fZVertexOffTemp(0)
  ,fXYvertexOff(0)	
  ,fXvertexOff(0)	    
  ,fYvertexOff(0)	    
  ,fZvertexOff(0)
  ,fEtaOff(0)
  ,fEtaMomentumcutOff(0)
  ,fNclusVSphiOff(0)
  ,fNclusVSthetaOff(0)
  ,fStatusOff(0)
  ,fStatusOff_Ocl(0)
  ,fEventSpecieOff(0)

  ,fChargeHLT(0)      
  ,fMomentumHLT(0)
  ,fMomentumHLTTpc(0)
  ,fMomentumHLTTpcIts(0)
  ,fDCAHLT(0)  
  ,fNclusterHLT(0)
  ,fNclusterHLTwCut(0)
  ,fdEdxHLT(0)    
  ,fdEdxVSPHLT(0)
  ,fPhiHLT(0)     
  ,fThetaHLT(0)  
  ,fMultHLT(0)  
  ,fXVertexVSNtracksHLT(0)
  ,fYVertexVSNtracksHLT(0)
  ,fZVertexVSNtracksHLT(0)
  ,fZVertexHLTTemp(0) 
  ,fXYvertexHLT(0)
  ,fXvertexHLT(0)
  ,fYvertexHLT(0)
  ,fZvertexHLT(0)
  ,fEtaHLT(0)
  ,fEtaMomentumcutHLT(0)
  ,fNclusVSphiHLT(0)	    
  ,fNclusVSthetaHLT(0)
  ,fStatusHLT(0)
  ,fStatusHLT_Ocl(0)    
  ,fEventSpecieHLT(0)
  ,fTrgClsArray(0)
  //     ,fDCAOff_trig(0)
  //     ,fNclusterOff_trig(0)
  //     
  //     ,fDCAHLT_trig(0)
  //     ,fNclusterHLT_trig(0)

{
  for(int jj=0;jj<7;jj++){
    fXvertexVSNcontriOff[jj]=0;
    fXvertexVSNcontriHLT[jj]=0;
 }
  for(int jj=0;jj<7;jj++){
    fYvertexVSNcontriOff[jj]=0;
     fYvertexVSNcontriHLT[jj]=0;
  }
  for(int jj=0;jj<7;jj++){
    fZvertexVSNcontriOff[jj]=0;
    fZvertexVSNcontriHLT[jj]=0;
  }

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

const Int_t  AliAnalysisTaskHLT::fNcontrArray[] ={3,6,9,12,15,18,21};
const Int_t  AliAnalysisTaskHLT::fNcontr = 7;
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

  fHistTrigger = new TH1F("fHistTrigger", "CTP trigger counter",12 , 0, 12);
  fHistTrigger->GetXaxis()->SetTitle("");  
  fHistTrigger->GetYaxis()->SetTitle("#Events"); 

  fHistHLTTrigger = new TH1F("fHistHLTTrigger", "HLT CTP trigger counter", 12, 0, 12); 
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
 
  fDCAOff = new TH1F("fDCA_off","DCA to beam line (offline)",200, -20, 20);
  fDCAHLT = new TH1F("fDCA_hlt","DCA to beam line (HLT)",    200, -20, 20);
 
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

  fXVertexVSNtracksOff = new TH2F("fXVertexNtracs_off", "X Vertex resolution vs nr contributing tracks (Offline)",21, 0, 20, 40, 0, 4); 
  fXVertexVSNtracksHLT = new TH2F("fXVertexNtracs_hlt", "X Vertex resolution vs nr contributing tracks (HLT)",21, 0, 20, 40, 0, 4); 

  fYVertexVSNtracksOff = new TH2F("fYVertexNtracs_off", "Y Vertex resolution vs nr contributing tracks (Offline)",21, 0, 20, 40, 0, 4); 
  fYVertexVSNtracksHLT = new TH2F("fYVertexNtracs_hlt", "Y Vertex resolution vs nr contributing tracks (HLT)",21, 0, 20, 40, 0, 4); 

  fZVertexVSNtracksOff = new TH2F("fZVertexNtracs_off", "Z Vertex resolution vs nr contributing tracks (Offline)",21, 0, 20, 40, 0, 10); 
  fZVertexVSNtracksHLT = new TH2F("fZVertexNtracs_hlt", "Z Vertex resolution vs nr contributing tracks (HLT)",21, 0, 20, 40, 0, 10); 

  fZVertexOffTemp = new TH1F("fZVertexTemo_off", "Temporary storage (Offline)", 250, -30, 30);
  fZVertexHLTTemp = new TH1F("fZVertexTemo_hlt", "Temporary storage (HLT)", 250, -30, 30);
  
  fXYvertexOff = new TH2F("fXYvertex_off","XY primary vertex (offline)",100,-5,5,100,-5,5);
  fXYvertexHLT = new TH2F("fXYvertex_hlt","XY primary vertex (HLT)",    100,-5,5,100,-5,5);
  
  fXvertexOff = new TH1F("fXvertex_off","X primary vertex (offline)",1000,-1,1);
  fXvertexHLT = new TH1F("fXvertex_hlt","X primary vertex (HLT)",    1000,-1,1);
 
  fYvertexOff = new TH1F("fYvertex_off","Y primary vertex (offline)",1000,-1,1);
  fYvertexHLT = new TH1F("fYvertex_hlt","Y primary vertex (HLT)",    1000,-1,1);
 
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
  
  fStatusOff_Ocl = new TH1F("fStatus_Ocl_off", "Status for different detectors when #TPCcl=0 (offline)",12, 0, 12);
  fStatusHLT_Ocl = new TH1F("fStatus_Ocl_hlt", "Status for different detectors when #TPCcl=0 (HLT)",    12, 0, 12);

  fEventSpecieOff = new TH1F("fEventSpecie_off","Eventspecie for Offline",18, 0, 18);
  fEventSpecieHLT = new TH1F("fEventSpecie_hlt","Eventspecie for HLT",18, 0, 18);

  for(int jj=0;jj<fNcontr;jj++){
  fXvertexVSNcontriOff[jj] = new TH1F(Form("fXVertex_vs_Ncontributors%d_off",fNcontrArray[jj]),Form("X position for vertex for %d contributors to vertex (Offline)"),250,-1,1);
  fYvertexVSNcontriOff[jj] = new TH1F(Form("fYVertex_vs_Ncontributors%d_off",fNcontrArray[jj]),Form("Y position for vertex for %d contributors to vertex (Offline)"),250,-1,1 );
  fZvertexVSNcontriOff[jj] = new TH1F(Form("fZVertex_vs_Ncontributors%d_off",fNcontrArray[jj]),Form("Z position for vertex for %d contributors to vertex (Offline)"),250,-30,30);
  fXvertexVSNcontriHLT[jj] = new TH1F(Form("fXVertex_vs_Ncontributors%d_hlt",fNcontrArray[jj]),Form("X position for vertex for %d contributors to vertex (HLT)"), 250, -1,1);
  fYvertexVSNcontriHLT[jj] = new TH1F(Form("fYVertex_vs_Ncontributors%d_hlt",fNcontrArray[jj]),Form("Y position for vertex for %d contributors to vertex (HLT)"), 250, -1,1);
  fZvertexVSNcontriHLT[jj] = new TH1F(Form("fZVertex_vs_Ncontributors%d_hlt",fNcontrArray[jj]),Form("Z position for vertex for %d contributors to vertex (HLT)"), 250,-30,30);
  }
  //---------------------- add histograms to the output TList ------------------//

  fOutputList->Add(fHistTrigger);
  fOutputList->Add(fHistHLTTrigger);

  fOutputList->Add(fChargeOff);
  fOutputList->Add(fMomentumOff);
  fOutputList->Add(fMomentumOffTpc);
  fOutputList->Add(fMomentumOffTpcIts);
  fOutputList->Add(fDCAOff);	  
  fOutputList->Add(fNclusterOff);
  fOutputList->Add(fNclusterOffwCut);
  fOutputList->Add(fdEdxOff);	  
  fOutputList->Add(fdEdxVSPOff);
  fOutputList->Add(fPhiOff);	  
  fOutputList->Add(fThetaOff);    
  fOutputList->Add(fMultOff);
  fOutputList->Add(fXVertexVSNtracksOff);
  fOutputList->Add(fYVertexVSNtracksOff);
  fOutputList->Add(fZVertexVSNtracksOff);
  fOutputList->Add(fZVertexOffTemp);
  fOutputList->Add(fXYvertexOff); 
  fOutputList->Add(fXvertexOff);  
  fOutputList->Add(fYvertexOff);  
  fOutputList->Add(fZvertexOff);  
  fOutputList->Add(fEtaOff);  
  fOutputList->Add(fEtaMomentumcutOff);
  fOutputList->Add(fNclusVSphiOff);  
  fOutputList->Add(fNclusVSthetaOff);
  fOutputList->Add(fStatusOff);
  fOutputList->Add(fStatusOff_Ocl);
  fOutputList->Add(fEventSpecieOff);

  fOutputList->Add(fChargeHLT);  
  fOutputList->Add(fMomentumHLT); 
  fOutputList->Add(fMomentumHLTTpc);  
  fOutputList->Add(fMomentumHLTTpcIts);
  fOutputList->Add(fDCAHLT);	  
  fOutputList->Add(fNclusterHLT); 
  fOutputList->Add(fNclusterHLTwCut); 
  fOutputList->Add(fdEdxHLT);	  
  fOutputList->Add(fdEdxVSPHLT);
  fOutputList->Add(fPhiHLT);	  
  fOutputList->Add(fThetaHLT);    
  fOutputList->Add(fMultHLT);	 
  fOutputList->Add(fXVertexVSNtracksHLT); 
  fOutputList->Add(fYVertexVSNtracksHLT); 
  fOutputList->Add(fZVertexVSNtracksHLT); 
  fOutputList->Add(fZVertexHLTTemp);
  fOutputList->Add(fXYvertexHLT); 
  fOutputList->Add(fXvertexHLT);  
  fOutputList->Add(fYvertexHLT);  
  fOutputList->Add(fZvertexHLT);    
  fOutputList->Add(fEtaHLT);  
  fOutputList->Add(fEtaMomentumcutHLT);    
  fOutputList->Add(fNclusVSphiHLT);  
  fOutputList->Add(fNclusVSthetaHLT);
  fOutputList->Add(fStatusHLT);
  fOutputList->Add(fStatusHLT_Ocl);
  fOutputList->Add(fEventSpecieHLT);

  for(int jj=0;jj<fNcontr;jj++){
    fOutputList->Add(fXvertexVSNcontriOff[jj]);
    fOutputList->Add(fXvertexVSNcontriHLT[jj]);
  }
  for(int jj=0;jj<fNcontr;jj++){
    fOutputList->Add(fYvertexVSNcontriOff[jj]);
    fOutputList->Add(fYvertexVSNcontriHLT[jj]);
  }

  for(int jj=0;jj<fNcontr;jj++){
    fOutputList->Add(fZvertexVSNcontriOff[jj]);
    fOutputList->Add(fZvertexVSNcontriHLT[jj]);
  }

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

  //Select only beam-interaction
  if(!(esdHLT->GetFiredTriggerClasses()).Contains("CINT1B-ABCE-NOPF-ALL"))
   return;

  Double_t DCAcut = 7.0;
  Double_t Momcut= 0.3;

  char Titlename[100];
  sprintf(Titlename,"pseudorapidity (HLT), DCA cut = %f,\n Momentum cut = %f, TPC clusters > 70" , DCAcut, Momcut);
  fEtaMomentumcutHLT->SetTitle(Titlename);
  sprintf(Titlename,"pseudorapidity (offline), DCA cut = %f, Momentum cut = %f, TPC clusters > 70", DCAcut, Momcut);
  fEtaMomentumcutOff->SetTitle(Titlename);

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
    fZVertexHLTTemp->Fill(vertHLT->GetZ());
    for(int kkk=0; kkk<fNcontr-1;kkk++){
      if(nr_contributorsHLT>=fNcontrArray[kkk] && nr_contributorsHLT<=fNcontrArray[kkk+1]){
	fXvertexVSNcontriHLT[kkk]->Fill(vertHLT->GetX());
	fYvertexVSNcontriHLT[kkk]->Fill(vertHLT->GetY());
	fZvertexVSNcontriHLT[kkk]->Fill(vertHLT->GetZ());
      }
    }
    if(nr_contributorsHLT>=fNcontrArray[fNcontr-1]){
      fXvertexVSNcontriHLT[fNcontr-1]->Fill(vertHLT->GetX());
      fYvertexVSNcontriHLT[fNcontr-1]->Fill(vertHLT->GetY());
      fZvertexVSNcontriHLT[fNcontr-1]->Fill(vertHLT->GetZ());
    }

  }
  //At the moment no constrains on vertex before filling histograms
  //Should be changed. 
  if(1){
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

      //reject laser events -> Use offline info (HLT do not set this flag)
      if((esdOFF->GetEventSpecie()==16)) continue;

      //This condition is mostly affecting Offline->will cut away tracks that are counted twice 
      //With both kITSin and kTPCin flags set.
      if(!(esdtrackHLT->GetStatus()&AliESDtrack::kTPCin)) continue;

      fDCAHLT->Fill(esdtrackHLT->GetD(esdHLT->GetPrimaryVertex()->GetXv(), esdHLT->GetPrimaryVertex()->GetYv(), bfield) );  
      fChargeHLT->Fill(esdtrackHLT->Charge());

      if((esdtrackHLT->GetStatus()&AliESDtrack::kTPCin) || (esdtrackHLT->GetStatus()&AliESDtrack::kTPCin && esdtrackHLT->GetStatus()&AliESDtrack::kTPCout))
	fNclusterHLT->Fill(esdtrackHLT->GetTPCNcls());

      //ESD-cut 
      //At the moment not included!
      //if(!fESDHLTtrackCuts->AcceptTrack(esdtrackHLT)) continue;		      
      
      if((esdtrackHLT->GetStatus()&AliESDtrack::kTPCin) || (esdtrackHLT->GetStatus()&AliESDtrack::kTPCin && esdtrackHLT->GetStatus()&AliESDtrack::kTPCout)){
	fNclusVSphiHLT->Fill(esdtrackHLT->Phi()*TMath::RadToDeg(), esdtrackHLT->GetTPCNcls());
	fNclusVSthetaHLT->Fill(esdtrackHLT->Theta()*TMath::RadToDeg(), esdtrackHLT->GetTPCNcls());
      }

      if(esdtrackHLT->GetStatus()&AliESDtrack::kTPCin) fMomentumHLTTpc->Fill(esdtrackHLT->Pt());
      if(esdtrackHLT->GetStatus()&AliESDtrack::kTPCin && esdtrackHLT->GetStatus()&AliESDtrack::kITSin)	 fMomentumHLTTpcIts->Fill(esdtrackHLT->Pt());   

      fEtaHLT->Fill(esdtrackHLT->Eta()); 
      fdEdxHLT->Fill( esdtrackHLT->GetTPCsignal() );
      fdEdxVSPHLT->Fill( TMath::Abs(esdtrackHLT->Pt()), esdtrackHLT->GetTPCsignal() );         

      //cut away tracks with mom<0.3GeV
      //if(TMath::Abs(esdtrackHLT->Pt()) <Momcut) continue; 
      fPhiHLT->Fill(esdtrackHLT->Phi()*TMath::RadToDeg());
      fThetaHLT->Fill(esdtrackHLT->Theta()*TMath::RadToDeg());
      fMomentumHLT->Fill( TMath::Abs(esdtrackHLT->Pt()) );  
     
      if(esdtrackHLT->GetStatus()&AliESDtrack::kTPCin) 
	fMomentumHLTTpc->Fill(esdtrackHLT->Pt());
      if(esdtrackHLT->GetStatus()&AliESDtrack::kTPCin && esdtrackHLT->GetStatus()&AliESDtrack::kITSin)	 
	fMomentumHLTTpcIts->Fill(esdtrackHLT->Pt());  

      if(TMath::Abs(esdtrackHLT->Pt()) <Momcut) continue; 
      fEtaMomentumcutHLT->Fill(esdtrackHLT->Eta());
      if((esdtrackHLT->GetStatus()&AliESDtrack::kTPCin) || (esdtrackHLT->GetStatus()&AliESDtrack::kTPCin && esdtrackHLT->GetStatus()&AliESDtrack::kTPCout))
	fNclusterHLTwCut->Fill(esdtrackHLT->GetTPCNcls());
      
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

  Int_t nr_contributorsOff= vertHLT->GetNContributors();

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

    for(int kkk=0; kkk<fNcontr-1;kkk++){
      if(nr_contributorsOff>=fNcontrArray[kkk] && nr_contributorsOff<=fNcontrArray[kkk+1] ){
	fXvertexVSNcontriOff[kkk]->Fill(vertOff->GetX() );
	fYvertexVSNcontriOff[kkk]->Fill(vertOff->GetY() );
	fZvertexVSNcontriOff[kkk]->Fill(vertOff->GetZ() );
      }
    }
    if(nr_contributorsOff>=fNcontrArray[fNcontr-1]){
      fXvertexVSNcontriOff[fNcontr-1]->Fill(vertOff->GetX());
      fYvertexVSNcontriOff[fNcontr-1]->Fill(vertOff->GetY());
      fZvertexVSNcontriOff[fNcontr-1]->Fill(vertOff->GetZ());
    }
  }

  //At the moment no constrains on vertex before filling histograms
  //Should be changed. 
  if(1){
    fEventSpecieOff->Fill(esdOFF->GetEventSpecie());

    for(Int_t i=0; i<esdOFF->GetNumberOfTracks(); i++){ 
     
      AliESDtrack *esdtrackOFF = esdOFF->GetTrack(i); 
      if (!esdtrackOFF) continue;


      //Fill histograms with which flags are set
      for(int jjj=0;jjj<12;jjj++){
	if(esdtrackOFF->GetStatus()&Statusnames[jjj]) {
	  fStatusOff->Fill(jjj);
	  if(esdtrackOFF->GetTPCNcls()==0) fStatusOff_Ocl->Fill(jjj);
	}
      } 

      // reject laser events
      if((esdOFF->GetEventSpecie()==16)) continue;  

      //This condition is mostly affecting Offline->will cut away tracks that are counted twice 
      //With both kITSin and kTPCin flags set.
      if(!(esdtrackOFF->GetStatus()&AliESDtrack::kTPCin))continue; 

      fDCAOff->Fill(esdtrackOFF->GetD(esdOFF->GetPrimaryVertex()->GetXv(), esdOFF->GetPrimaryVertex()->GetYv(), bfield) ); 
      fChargeOff->Fill(esdtrackOFF->Charge());

      if((esdtrackOFF->GetStatus()&AliESDtrack::kTPCin) || (esdtrackOFF->GetStatus()&AliESDtrack::kTPCin && esdtrackOFF->GetStatus()&AliESDtrack::kTPCout))
	fNclusterOff->Fill(esdtrackOFF->GetTPCNcls()); 

      // -- ESD cuts
      //Not included at the moment
      //if(!fESDOfftrackCuts->AcceptTrack(esdtrackOFF) ) continue;

      if(esdtrackOFF->GetTPCNcls()>0) fNclusVSphiOff->Fill(esdtrackOFF->Phi()*TMath::RadToDeg(), esdtrackOFF->GetTPCNcls());
      if(esdtrackOFF->GetTPCNcls()>0) fNclusVSthetaOff->Fill(esdtrackOFF->Theta()*TMath::RadToDeg(), esdtrackOFF->GetTPCNcls());

      if(esdtrackOFF->GetStatus()&AliESDtrack::kTPCin) fMomentumOffTpc->Fill(esdtrackOFF->Pt());
      if(esdtrackOFF->GetStatus()&AliESDtrack::kTPCin && esdtrackOFF->GetStatus()&AliESDtrack::kITSin)	 fMomentumOffTpcIts->Fill(esdtrackOFF->Pt());   

      fEtaOff->Fill(esdtrackOFF->Eta());	  
      fdEdxOff->Fill( esdtrackOFF->GetTPCsignal() );
      fdEdxVSPOff->Fill( TMath::Abs(esdtrackOFF->Pt()), esdtrackOFF->GetTPCsignal() );         
      
      //cut away tracks with mom<0.3GeV
      //if(TMath::Abs(esdtrackOFF->Pt()) < Momcut) continue;
      fPhiOff->Fill(esdtrackOFF->Phi()*TMath::RadToDeg());
      fThetaOff->Fill(esdtrackOFF->Theta()*TMath::RadToDeg());
      fMomentumOff->Fill( TMath::Abs(esdtrackOFF->Pt()) ); 

      if(esdtrackOFF->GetStatus()&AliESDtrack::kTPCin) 
	fMomentumOffTpc->Fill(esdtrackOFF->Pt());
      if(esdtrackOFF->GetStatus()&AliESDtrack::kTPCin && esdtrackOFF->GetStatus()&AliESDtrack::kITSin)
	fMomentumOffTpcIts->Fill(esdtrackOFF->Pt());  

      if(TMath::Abs(esdtrackOFF->Pt()) < Momcut) continue;
      fEtaMomentumcutOff->Fill(esdtrackOFF->Eta()); 
      if(esdtrackOFF->GetTPCNcls()>0) 
	fNclusterOffwCut->Fill(esdtrackOFF->GetTPCNcls()); 
	    
	  
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
  
    TF1 *signalX= new TF1("signalX","gaus",-1,1); 
    TF1 *signalY= new TF1("signalY","gaus",-1,1);
    TF1 *signalZ= new TF1("signalZ","gaus",-30,30);

    signalZ->SetLineColor(8);
    signalZ->SetLineWidth(2);
    
    fXVertexVSNtracksHLT->SetMarkerStyle(7);
    fXVertexVSNtracksHLT->SetMarkerSize(7);
    fXVertexVSNtracksOff->SetMarkerStyle(7);
    fXVertexVSNtracksOff->SetMarkerSize(7);

    fYVertexVSNtracksHLT->SetMarkerStyle(7);
    fYVertexVSNtracksHLT->SetMarkerSize(7);
    fYVertexVSNtracksOff->SetMarkerStyle(7);
    fYVertexVSNtracksOff->SetMarkerSize(7);

    fZVertexVSNtracksHLT->Sumw2();
    fZVertexVSNtracksHLT->SetMarkerStyle(21);
    fZVertexVSNtracksHLT->SetMarkerSize(0.7);
    fZVertexVSNtracksOff->Sumw2();
    fZVertexVSNtracksOff->SetMarkerStyle(21);
    fZVertexVSNtracksOff->SetMarkerSize(0.7);

    for(int jj=0;jj<fNcontr;jj++){

      fXvertexVSNcontriOff[jj]->Fit(signalX,"RQ");
      fXVertexVSNtracksOff->Fill(fNcontrArray[jj],signalX->GetParameter(2)); 

      fYvertexVSNcontriOff[jj]->Fit(signalY,"RQ");
      fYVertexVSNtracksOff->Fill(fNcontrArray[jj],signalY->GetParameter(2));

      fZvertexVSNcontriOff[jj]->Fit(signalZ,"RQ");
      fZVertexVSNtracksOff->Fill(fNcontrArray[jj],signalZ->GetParameter(2));
      fZVertexVSNtracksOff->Fill(fNcontrArray[jj],signalZ->GetParameter(2), signalZ->GetParError(2));
    
      fXvertexVSNcontriHLT[jj]->Fit(signalX,"RQ");
      fXVertexVSNtracksHLT->Fill(fNcontrArray[jj],signalX->GetParameter(2));

      fYvertexVSNcontriHLT[jj]->Fit(signalY,"RQ");
      fYVertexVSNtracksHLT->Fill(fNcontrArray[jj],signalY->GetParameter(2));
    
      fZvertexVSNcontriHLT[jj]->Fit(signalZ,"RQ");
      fZVertexVSNtracksHLT->Fill(fNcontrArray[jj],signalZ->GetParameter(2));
      fZVertexVSNtracksHLT->Fill(fNcontrArray[jj],signalZ->GetParameter(2), signalZ->GetParError(2));

    }

    //Drawing histograms
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
    maxbin =fEtaMomentumcutOff->GetBinContent(fEtaMomentumcutOff->GetMaximumBin());
    if(maxbin < fEtaMomentumcutHLT->GetBinContent(fEtaMomentumcutHLT->GetMaximumBin()))
      maxbin=fEtaMomentumcutHLT->GetBinContent(fEtaMomentumcutHLT->GetMaximumBin());
    fEtaMomentumcutOff->SetMaximum(maxbin+20);
    fEtaMomentumcutOff->SetTitle("Pseudorapidity");
    fEtaMomentumcutOff->SetLineColor(2);
    fEtaMomentumcutOff->DrawCopy("");
    fEtaMomentumcutHLT->DrawCopy("sames");

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
