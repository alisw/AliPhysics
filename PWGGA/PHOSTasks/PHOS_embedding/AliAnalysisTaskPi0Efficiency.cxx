/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
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

#include "TChain.h"
#include "TTree.h"
#include "TObjArray.h"
#include "THashList.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH3F.h"
#include "TParticle.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom.h"

#include "AliAODMCParticle.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskPi0Efficiency.h"
#include "AliCaloPhoton.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSAodCluster.h"
#include "AliPHOSCalibData.h"
#include "AliAODEvent.h"
#include "AliAODCaloCluster.h"
#include "AliAODVertex.h"
#include "AliESDtrackCuts.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliCDBManager.h"
#include "AliCentrality.h" 

// Analysis task to fill histograms with PHOS ESD clusters and cells
// Authors: Yuri Kharlov
// Date   : 28.05.2009

ClassImp(AliAnalysisTaskPi0Efficiency)

//________________________________________________________________________
AliAnalysisTaskPi0Efficiency::AliAnalysisTaskPi0Efficiency(const char *name) 
: AliAnalysisTaskSE(name),
  fStack(0),
  fOutputContainer(0),
  fPHOSEvent(0),
  fPHOSCalibData(0),
  fNonLinCorr(0),
  fRPfull(0),
  fRPA(0),
  fRPC(0),
  fRPFar(0),
  fRPAFar(0),
  fRPCFar(0),
  fCentrality(0),
  fCenBin(0),
  fPHOSGeo(0),
  fEventCounter(0)
{
  // Constructor
  for(Int_t i=0;i<1;i++){
    for(Int_t j=0;j<10;j++)
      for(Int_t k=0;k<11;k++)
	fPHOSEvents[i][j][k]=0 ;
  }
  
  // Output slots #0 write into a TH1 container
  DefineOutput(1,THashList::Class());

  // Set bad channel map
  char key[55] ;
  for(Int_t i=0; i<6; i++){
    snprintf(key,55,"PHOS_BadMap_mod%d",i) ;
    fPHOSBadMap[i]=new TH2I(key,"Bad Modules map",64,0.,64.,56,0.,56.) ;
  }
  // Initialize the PHOS geometry
  fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP") ;

  fPHOSCalibData = new AliPHOSCalibData();
  for(Int_t module=1; module<=5; module++) {
    for(Int_t column=1; column<=56; column++) {
      for(Int_t row=1; row<=64; row++) {
        fPHOSCalibData->SetADCchannelEmc(module,column,row,1.);
      }
    }
  }


}

//________________________________________________________________________
void AliAnalysisTaskPi0Efficiency::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  // ESD histograms
  if(fOutputContainer != NULL){
    delete fOutputContainer;
  }
  fOutputContainer = new THashList();
  fOutputContainer->SetOwner(kTRUE);

  //Event selection
  fOutputContainer->Add(new TH1F("hSelEvents","Event celection", 10,0.,10.)) ;

  //vertex distribution
  fOutputContainer->Add(new TH1F("hZvertex","Z vertex position", 50,-25.,25.)) ;

  //Centrality
  fOutputContainer->Add(new TH1F("hCentrality","Event centrality", 100,0.,100.)) ;

  //QA histograms			
  fOutputContainer->Add(new TH2F("hCluM1","Cell (X,Z), M1" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluM2","Cell (X,Z), M2" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluM3","Cell (X,Z), M3" ,64,0.5,64.5, 56,0.5,56.5));

  Int_t nM       = 500;
  Double_t mMin  = 0.0;
  Double_t mMax  = 1.0;
  Int_t nPt      = 200;
  Double_t ptMin = 0;
  Double_t ptMax = 20;

  char key[55] ;
  for(Int_t cent=0; cent<6; cent++){
    //Single photon
    snprintf(key,55,"hPhotAll_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    snprintf(key,55,"hPhotAllwou_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    snprintf(key,55,"hPhotAllcore_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    snprintf(key,55,"hPhotCPV_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    snprintf(key,55,"hPhotCPVcore_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    snprintf(key,55,"hPhotCPV2_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    snprintf(key,55,"hPhotDisp_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    snprintf(key,55,"hPhotDispwou_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    snprintf(key,55,"hPhotDisp2_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    snprintf(key,55,"hPhotBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    snprintf(key,55,"hPhotBothcore_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    
    snprintf(key,55,"hMassPtAll_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMassPtAllwou_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMassPtAllcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMassPtCPV_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMassPtCPVcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMassPtCPV2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMassPtDisp_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMassPtDispwou_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMassPtDisp2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMassPtBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMassPtBothcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));

    snprintf(key,55,"hMassPtAll_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMassPtCPV_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMassPtCPV2_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMassPtDisp_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMassPtBoth_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));

    snprintf(key,55,"hMassPtAll_a08_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMassPtCPV_a08_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMassPtCPV2_a08_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMassPtDisp_a08_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMassPtBoth_a08_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    
    snprintf(key,55,"hMassPtAll_a09_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMassPtCPV_a09_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMassPtCPV2_a09_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMassPtDisp_a09_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMassPtBoth_a09_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    
    
    //Mixed
    snprintf(key,55,"hMiMassPtAll_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtAllwou_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtAllcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtCPV_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtCPVcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtCPV2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtDisp_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtDispwou_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtDisp2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtBothcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));

    snprintf(key,55,"hMiMassPtAll_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtCPV_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtCPV2_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtDisp_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtBoth_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));

    snprintf(key,55,"hMiMassPtAll_a08_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtCPV_a08_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtCPV2_a08_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtDisp_a08_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtBoth_a08_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    
     snprintf(key,55,"hMiMassPtAll_a09_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtCPV_a09_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtCPV2_a09_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtDisp_a09_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtBoth_a09_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
   
    //MC
    snprintf(key,55,"hMCMassPtAll_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtAllwou_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtAllcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtCPV_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtCPVcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtCPV2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtDisp_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtDispwou_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtDisp2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtBothcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));  
    
    
    snprintf(key,55,"hMCMassPtAll_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtCPV_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtCPV2_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtDisp_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtBoth_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));

    snprintf(key,55,"hMCMassPtAll_a08_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtCPV_a08_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtCPV2_a08_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtDisp_a08_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtBoth_a08_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));

    snprintf(key,55,"hMCMassPtAll_a09_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtCPV_a09_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtCPV2_a09_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtDisp_a09_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtBoth_a09_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));

    //Single photon
    snprintf(key,55,"hMCPhotAll_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCPhotAllwou_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCPhotAllcore_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCPhotCPV_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCPhotCPVcore_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCPhotCPV2_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCPhotDisp_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCPhotDispwou_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCPhotDisp2_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCPhotBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCPhotBothcore_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));

  }
    
  fOutputContainer->Add(new TH2F("hMCPi0M11","(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMCPi0M22","(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMCPi0M33","(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMCPi0M12","(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMCPi0M13","(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMCPi0M23","(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));


  //MC
  for(Int_t cent=0; cent<6; cent++){
    snprintf(key,55,"hMC_rap_gamma_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Rapidity pi0",200,-1.,1.)) ;
    snprintf(key,55,"hMC_rap_pi0_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Rapidity pi0",200,-1.,1.)) ;
    snprintf(key,55,"hMC_rap_eta_cen%d",cent) ;
    fOutputContainer->Add(new TH1F("hMC_rap_eta","Rapidity eta",200,-1.,1.)) ;
    snprintf(key,55,"hMC_phi_gamma_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Phi pi0",200,0.,TMath::TwoPi())) ;
    snprintf(key,55,"hMC_phi_pi0_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Phi pi0",200,0.,TMath::TwoPi())) ;
    snprintf(key,55,"hMC_phi_eta_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Phi eta",200,0.,TMath::TwoPi())) ;
    snprintf(key,55,"hMC_all_gamma_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Rapidity photon",250,0.,25.)) ;
    snprintf(key,55,"hMC_all_pi0_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Rapidity pi0",250,0.,25.)) ;
    snprintf(key,55,"hMC_all_eta_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Rapidity eta",250,0.,25.)) ;
    snprintf(key,55,"hMC_unitEta_gamma_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Pt photon",250,0.,25.)) ;
    snprintf(key,55,"hMC_unitEta_pi0_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Rapidity eta",250,0.,25.)) ;
    snprintf(key,55,"hMC_unitEta_eta_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Rapidity eta",250,0.,25.)) ;
  }
  
  PostData(1, fOutputContainer);

}

//________________________________________________________________________
void AliAnalysisTaskPi0Efficiency::UserExec(Option_t *) 
{
  // Main loop, called for each event
  // Analyze ESD/AOD

  FillHistogram("hSelEvents",0.5) ;  
  
  AliAODEvent *event = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!event) {
    Printf("ERROR: Could not retrieve event");
    PostData(1, fOutputContainer);
    return;
  }
  
  FillHistogram("hSelEvents",1.5) ;
  AliAODHeader *header = event->GetHeader() ;
  
  // Checks if we have a primary vertex
  // Get primary vertices form ESD
  const AliAODVertex *esdVertex5 = event->GetPrimaryVertex();

 // don't rely on ESD vertex, assume (0,0,0)
  Double_t vtx0[3] ={0.,0.,0.};
  
  
  FillHistogram("hZvertex",esdVertex5->GetZ());
  if (TMath::Abs(esdVertex5->GetZ()) > 10. ){
    PostData(1, fOutputContainer);
    return;
  }
  FillHistogram("hSelEvents",2.5) ;

  //Vtx class z-bin
  //  Int_t zvtx = (Int_t)((vtx5[2]+10.)/2.) ;
  //  if(zvtx<0)zvtx=0 ;
  //  if(zvtx>9)zvtx=9 ;
  Int_t zvtx=0 ;

//  fCentrality=header->GetCentralityP()->GetCentralityPercentile("V0M"); // returns the centrality percentile, 
//                                                          //a float from 0 to 100 (or to the trigger efficiency)
   fCentrality=header->GetZDCN2Energy() ;

  if( fCentrality < 0. || fCentrality>80.){
    PostData(1, fOutputContainer);
    return;
  }
  FillHistogram("hSelEvents",3.5) ;
  Float_t bins[7]={0.,5.,10.,20.,40.,60.,80.} ;
  fCenBin=0 ;
  while(fCenBin<6 && fCentrality > bins[fCenBin+1])
    fCenBin++ ; 

 
  //reaction plain
  fRPfull= header->GetZDCN1Energy() ;
  if(fRPfull==999){ //reaction plain was not defined
    PostData(1, fOutputContainer);
    return;
  } 

  FillHistogram("hSelEvents",4.5) ;
  //All event selections done
  FillHistogram("hCentrality",fCentrality) ;
  //Reaction plain is defined in the range (-pi/2;pi/2)
  //We have 10 bins
  Int_t irp=Int_t(10.*fRPfull/TMath::Pi());
  if(irp>9)irp=9 ;

  if(!fPHOSEvents[zvtx][fCenBin][irp]) 
    fPHOSEvents[zvtx][fCenBin][irp]=new TList() ;
  TList * prevPHOS = fPHOSEvents[zvtx][fCenBin][irp] ;

  // Get PHOS rotation matrices from ESD and set them to the PHOS geometry
  if(fEventCounter == 0) {
    for(Int_t mod=0; mod<5; mod++) {
      const TGeoHMatrix* m =header->GetPHOSMatrix(mod) ;
      fPHOSGeo->SetMisalMatrix(m,mod) ;
      Printf("PHOS geo matrix for module # %d is set: %p\n", mod,m);
    }
    fEventCounter++ ;
  }

  ProcessMC() ;

  if(fPHOSEvent)
    fPHOSEvent->Clear() ;
  else
    fPHOSEvent = new TClonesArray("AliCaloPhoton",200) ;


  char key[55] ;
  Int_t inPHOS=0 ;
  TVector3 vertex(vtx0);
  TClonesArray * clusters = (TClonesArray*)event->FindListObject("EmbeddedCaloClusters") ;
  AliAODCaloCells * cells = (AliAODCaloCells*)event->FindListObject("EmbeddedPHOScells") ;
  Int_t multClust = clusters->GetEntriesFast();
  for (Int_t i=0; i<multClust; i++) {
    AliAODCaloCluster *clu = (AliAODCaloCluster*) clusters->At(i);
    if ( !clu->IsPHOS() || clu->E()<0.3) continue;

    Float_t  position[3];
    clu->GetPosition(position);
    TVector3 global(position) ;
    Int_t relId[4] ;
    fPHOSGeo->GlobalPos2RelId(global,relId) ;
    Int_t mod  = relId[0] ;
    Int_t cellX = relId[2];
    Int_t cellZ = relId[3] ;
    if ( !IsGoodChannel("PHOS",mod,cellX,cellZ) ) 
      continue ;
    if(clu->GetNCells()<3)
      continue ;
    if(clu->GetM02()<0.2)
      continue ;

    snprintf(key,55,"hCluM%d",mod) ;
    FillHistogram(key,cellX,cellZ,1.);

    TLorentzVector pv1 ;
    clu->GetMomentum(pv1 ,vtx0);
    
    if(inPHOS>=fPHOSEvent->GetSize()){
      fPHOSEvent->Expand(inPHOS+50) ;
    }
    new((*fPHOSEvent)[inPHOS]) AliCaloPhoton(pv1.X(),pv1.Py(),pv1.Z(),pv1.E()) ;
    AliCaloPhoton * ph = (AliCaloPhoton*)fPHOSEvent->At(inPHOS) ;
    ph->SetModule(mod) ;
    AliPHOSAodCluster cluPHOS1(*clu);
    cluPHOS1.Recalibrate(fPHOSCalibData,cells); // modify the cell energies
    Double_t ecore=CoreEnergy(&cluPHOS1) ; 
    pv1*= ecore/pv1.E() ;
    ph->SetMomV2(&pv1) ;
    ph->SetNCells(clu->GetNCells());
    ph->SetDispBit(TestLambda(clu->E(),clu->GetM20(),clu->GetM02())) ;
    ph->SetDisp2Bit(TestLambda2(clu->E(),clu->GetM20(),clu->GetM02())) ;
    ph->SetCPVBit(clu->GetEmcCpvDistance()>2.) ;
    ph->SetCPV2Bit(clu->GetEmcCpvDistance()>4.) ;
    ph->SetPhoton(clu->GetNExMax()<2); // Remember, if it is unfolded

    inPHOS++ ;
  }
  //Single photon
  for (Int_t i1=0; i1<inPHOS; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1) ;
    snprintf(key,55,"hPhotAll_cen%d",fCenBin) ;
    FillHistogram(key,ph1->Pt()) ;
    snprintf(key,55,"hPhotAllcore_cen%d",fCenBin) ;
    FillHistogram(key,ph1->GetMomV2()->Pt()) ;
    if(ph1->IsPhoton()){
      snprintf(key,55,"hPhotAllwou_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt()) ;
    }
    if(ph1->IsCPVOK() ){
      snprintf(key,55,"hPhotCPV_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt()) ;
      snprintf(key,55,"hPhotCPVcore_cen%d",fCenBin) ;
      FillHistogram(key,ph1->GetMomV2()->Pt()) ;
    }
    if(ph1->IsCPV2OK() ){
      snprintf(key,55,"hPhotCPV2_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt()) ;
    }
    if(ph1->IsDisp2OK()){
      snprintf(key,55,"hPhotDisp2_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt()) ;
    }
    if(ph1->IsDispOK()){
      snprintf(key,55,"hPhotDisp_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt()) ;
      if(ph1->IsPhoton()){
        snprintf(key,55,"hPhotDispwou_cen%d",fCenBin) ;
        FillHistogram(key,ph1->Pt()) ;
      }
      if(ph1->IsCPVOK()){
	snprintf(key,55,"hPhotBoth_cen%d",fCenBin) ;
	FillHistogram(key,ph1->Pt()) ;
	snprintf(key,55,"hPhotBothcore_cen%d",fCenBin) ;
	FillHistogram(key,ph1->GetMomV2()->Pt()) ;
      }
    } // end of loop i2
  } // end of loop i1 

  // Fill Real disribution
  for (Int_t i1=0; i1<inPHOS-1; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1) ;
    for (Int_t i2=i1+1; i2<inPHOS; i2++) {
      AliCaloPhoton * ph2=(AliCaloPhoton*)fPHOSEvent->At(i2) ;
      TLorentzVector p12  = *ph1  + *ph2;
      TLorentzVector pv12 = *(ph1->GetMomV2()) + *(ph2->GetMomV2());      
      Double_t a=TMath::Abs((ph1->E()-ph2->E())/(ph1->E()+ph2->E())) ;
      
      snprintf(key,55,"hMassPtAll_cen%d",fCenBin) ;
      FillHistogram(key,p12.M() ,p12.Pt()) ;
      snprintf(key,55,"hMassPtAllcore_cen%d",fCenBin) ;
      FillHistogram(key,pv12.M(), pv12.Pt()) ;
      if(ph1->IsPhoton() && ph2->IsPhoton()){
        snprintf(key,55,"hMassPtAllwou_cen%d",fCenBin) ;
        FillHistogram(key,p12.M() ,p12.Pt()) ;
      }
      if(a<0.9){
        snprintf(key,55,"hMassPtAll_a09_cen%d",fCenBin) ;
        FillHistogram(key,p12.M() ,p12.Pt()) ;
        if(a<0.8){
          snprintf(key,55,"hMassPtAll_a08_cen%d",fCenBin) ;
          FillHistogram(key,p12.M() ,p12.Pt()) ;
          if(a<0.7){
            snprintf(key,55,"hMassPtAll_a07_cen%d",fCenBin) ;
            FillHistogram(key,p12.M() ,p12.Pt()) ;
          }
        }
      }
      if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	snprintf(key,55,"hMassPtCPV_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt()) ;
	snprintf(key,55,"hMassPtCPVcore_cen%d",fCenBin) ;
	FillHistogram(key,pv12.M(), pv12.Pt()) ;
        if(a<0.9){
          snprintf(key,55,"hMassPtCPV_a09_cen%d",fCenBin) ;
          FillHistogram(key,p12.M() ,p12.Pt()) ;
          if(a<0.8){
            snprintf(key,55,"hMassPtCPV_a08_cen%d",fCenBin) ;
            FillHistogram(key,p12.M() ,p12.Pt()) ;
            if(a<0.7){
              snprintf(key,55,"hMassPtCPV_a07_cen%d",fCenBin) ;
              FillHistogram(key,p12.M() ,p12.Pt()) ;
            }
          }
        }
      }
      if(ph1->IsCPV2OK() && ph2->IsCPV2OK()){
	snprintf(key,55,"hMassPtCPV2_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt()) ;
        if(a<0.9){
          snprintf(key,55,"hMassPtCPV2_a09_cen%d",fCenBin) ;
          FillHistogram(key,p12.M() ,p12.Pt()) ;
          if(a<0.8){
            snprintf(key,55,"hMassPtCPV2_a08_cen%d",fCenBin) ;
            FillHistogram(key,p12.M() ,p12.Pt()) ;
            if(a<0.7){
              snprintf(key,55,"hMassPtCPV2_a07_cen%d",fCenBin) ;
              FillHistogram(key,p12.M() ,p12.Pt()) ;
            }
          }
        }
      }
      if(ph1->IsDisp2OK() && ph2->IsDisp2OK()){
	snprintf(key,55,"hMassPtDisp2_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt()) ;
      }
      if(ph1->IsDispOK() && ph2->IsDispOK()){
	snprintf(key,55,"hMassPtDisp_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt()) ;
        if(ph1->IsPhoton() && ph2->IsPhoton()){
  	  snprintf(key,55,"hMassPtDispwou_cen%d",fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt()) ;
	}
        if(a<0.9){
          snprintf(key,55,"hMassPtDisp_a09_cen%d",fCenBin) ;
          FillHistogram(key,p12.M() ,p12.Pt()) ;
          if(a<0.8){
            snprintf(key,55,"hMassPtDisp_a08_cen%d",fCenBin) ;
            FillHistogram(key,p12.M() ,p12.Pt()) ;
            if(a<0.7){
              snprintf(key,55,"hMassPtDisp_a07_cen%d",fCenBin) ;
              FillHistogram(key,p12.M() ,p12.Pt()) ;
            }
          }
        }

	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  snprintf(key,55,"hMassPtBoth_cen%d",fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt()) ;
	  snprintf(key,55,"hMassPtBothcore_cen%d",fCenBin) ;
	  FillHistogram(key,pv12.M(), pv12.Pt()) ;
          if(a<0.9){
            snprintf(key,55,"hMassPtBoth_a09_cen%d",fCenBin) ;
            FillHistogram(key,p12.M() ,p12.Pt()) ;
            if(a<0.8){
              snprintf(key,55,"hMassPtBoth_a08_cen%d",fCenBin) ;
              FillHistogram(key,p12.M() ,p12.Pt()) ;
              if(a<0.7){
                snprintf(key,55,"hMassPtBoth_a07_cen%d",fCenBin) ;
                FillHistogram(key,p12.M() ,p12.Pt()) ;
             }
            }
          }
        }
      }
    } // end of loop i2
  } // end of loop i1 

  //now mixed
  for (Int_t i1=0; i1<inPHOS; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1) ;
    for(Int_t ev=0; ev<prevPHOS->GetSize();ev++){
      TClonesArray * mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev)) ;
      for(Int_t i2=0; i2<mixPHOS->GetEntriesFast();i2++){
	AliCaloPhoton * ph2=(AliCaloPhoton*)mixPHOS->At(i2) ;
	TLorentzVector p12  = *ph1  + *ph2;
	TLorentzVector pv12 = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
        Double_t a=TMath::Abs((ph1->E()-ph2->E())/(ph1->E()+ph2->E())) ;
	
	snprintf(key,55,"hMiMassPtAll_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt()) ;
	snprintf(key,55,"hMiMassPtAllcore_cen%d",fCenBin) ;
	FillHistogram(key,pv12.M(), pv12.Pt()) ;
        if(ph1->IsPhoton() && ph2->IsPhoton()){
	  snprintf(key,55,"hMiMassPtAllwou_cen%d",fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt()) ;
	}
        if(a<0.9){
          snprintf(key,55,"hMiMassPtAll_a09_cen%d",fCenBin) ;
          FillHistogram(key,p12.M() ,p12.Pt()) ;
          if(a<0.8){
            snprintf(key,55,"hMiMassPtAll_a08_cen%d",fCenBin) ;
            FillHistogram(key,p12.M() ,p12.Pt()) ;
            if(a<0.7){
              snprintf(key,55,"hMiMassPtAll_a07_cen%d",fCenBin) ;
              FillHistogram(key,p12.M() ,p12.Pt()) ;
            }
          }
        }
  	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  snprintf(key,55,"hMiMassPtCPV_cen%d",fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt()) ;
	  snprintf(key,55,"hMiMassPtCPVcore_cen%d",fCenBin) ;
	  FillHistogram(key,pv12.M(), pv12.Pt()) ;
          if(a<0.9){
            snprintf(key,55,"hMiMassPtCPV_a09_cen%d",fCenBin) ;
            FillHistogram(key,p12.M() ,p12.Pt()) ;
            if(a<0.8){
              snprintf(key,55,"hMiMassPtCPV_a08_cen%d",fCenBin) ;
              FillHistogram(key,p12.M() ,p12.Pt()) ;
              if(a<0.7){
                snprintf(key,55,"hMiMassPtCPV_a07_cen%d",fCenBin) ;
                FillHistogram(key,p12.M() ,p12.Pt()) ;
              }
            }
          }
	}
	if(ph1->IsCPV2OK() && ph2->IsCPV2OK()){
	  snprintf(key,55,"hMiMassPtCPV2_cen%d",fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt()) ;
          if(a<0.9){
            snprintf(key,55,"hMiMassPtCPV2_a09_cen%d",fCenBin) ;
            FillHistogram(key,p12.M() ,p12.Pt()) ;
            if(a<0.8){
              snprintf(key,55,"hMiMassPtCPV2_a08_cen%d",fCenBin) ;
              FillHistogram(key,p12.M() ,p12.Pt()) ;
              if(a<0.7){
                snprintf(key,55,"hMiMassPtCPV2_a07_cen%d",fCenBin) ;
                FillHistogram(key,p12.M() ,p12.Pt()) ;
              }
            }
          }
	}
	if(ph1->IsDisp2OK() && ph2->IsDisp2OK()){
	  snprintf(key,55,"hMiMassPtDisp2_cen%d",fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt()) ;
	}
	if(ph1->IsDispOK() && ph2->IsDispOK()){
	  snprintf(key,55,"hMiMassPtDisp_cen%d",fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt()) ;
	  if(ph1->IsPhoton() && ph2->IsPhoton()){
  	    snprintf(key,55,"hMiMassPtDispwou_cen%d",fCenBin) ;
	    FillHistogram(key,p12.M() ,p12.Pt()) ;
	  }
          if(a<0.9){
            snprintf(key,55,"hMiMassPtDisp_a09_cen%d",fCenBin) ;
            FillHistogram(key,p12.M() ,p12.Pt()) ;
            if(a<0.8){
              snprintf(key,55,"hMiMassPtDisp_a08_cen%d",fCenBin) ;
              FillHistogram(key,p12.M() ,p12.Pt()) ;
              if(a<0.7){
                snprintf(key,55,"hMiMassPtDisp_a07_cen%d",fCenBin) ;
                FillHistogram(key,p12.M() ,p12.Pt()) ;
              }
            }
	  }
	  if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	    snprintf(key,55,"hMiMassPtBoth_cen%d",fCenBin) ;
	    FillHistogram(key,p12.M() ,p12.Pt()) ;
	    snprintf(key,55,"hMiMassPtBothcore_cen%d",fCenBin) ;
	    FillHistogram(key,pv12.M(), pv12.Pt()) ;
            if(a<0.9){
              snprintf(key,55,"hMiMassPtBoth_a09_cen%d",fCenBin) ;
              FillHistogram(key,p12.M() ,p12.Pt()) ;
              if(a<0.8){
                snprintf(key,55,"hMiMassPtBoth_a08_cen%d",fCenBin) ;
                FillHistogram(key,p12.M() ,p12.Pt()) ;
                if(a<0.7){
                  snprintf(key,55,"hMiMassPtBoth_a07_cen%d",fCenBin) ;
                  FillHistogram(key,p12.M() ,p12.Pt()) ;
               }
              }
            }
	  }
	}
      } // end of loop i2
    }
  } // end of loop i1
  
  
  //Now we either add current events to stack or remove
  //If no photons in current event - no need to add it to mixed
  if(fPHOSEvent->GetEntriesFast()>0){
    prevPHOS->AddFirst(fPHOSEvent) ;
    fPHOSEvent=0;
    if(prevPHOS->GetSize()>100){//Remove redundant events
      TClonesArray * tmp = static_cast<TClonesArray*>(prevPHOS->Last()) ;
      prevPHOS->RemoveLast() ;
      delete tmp ;
    }
  }
  // Post output data.
  PostData(1, fOutputContainer);
  fEventCounter++;
}

//________________________________________________________________________
void AliAnalysisTaskPi0Efficiency::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
  
}

//________________________________________________________________________
Bool_t AliAnalysisTaskPi0Efficiency::IsGoodChannel(const char * det, Int_t mod, Int_t ix, Int_t iz)
{
  //Check if this channel belogs to the good ones

  if(strcmp(det,"PHOS")==0){
    if(mod>5 || mod<1){
      AliError(Form("No bad map for PHOS module %d ",mod)) ;
      return kTRUE ;
    }
    if(!fPHOSBadMap[mod]){
      AliError(Form("No Bad map for PHOS module %d",mod)) ;
      return kTRUE ;
    }
    if(fPHOSBadMap[mod]->GetBinContent(ix,iz)>0)
      return kFALSE ;
    else
      return kTRUE ;
  }
  else{
    AliError(Form("Can not find bad channels for detector %s ",det)) ;
  }
  return kTRUE ;
}
//_____________________________________________________________________________
void AliAnalysisTaskPi0Efficiency::FillHistogram(const char * key,Double_t x)const{
  //FillHistogram
  TH1 * hist = dynamic_cast<TH1*>(fOutputContainer->FindObject(key)) ;
  if(hist)
    hist->Fill(x) ;
  else
    AliError(Form("can not find histogram (of instance TH1) <%s> ",key)) ;
}
//_____________________________________________________________________________
void AliAnalysisTaskPi0Efficiency::FillHistogram(const char * key,Double_t x,Double_t y)const{
  //FillHistogram
  TH1 * th1 = dynamic_cast<TH1*> (fOutputContainer->FindObject(key));
  if(th1)
    th1->Fill(x, y) ;
  else
    AliError(Form("can not find histogram (of instance TH1) <%s> ",key)) ;
}

//_____________________________________________________________________________
void AliAnalysisTaskPi0Efficiency::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z) const{
  //Fills 1D histograms with key
  TObject * obj = fOutputContainer->FindObject(key);
  
  TH2 * th2 = dynamic_cast<TH2*> (obj);
  if(th2) {
    th2->Fill(x, y, z) ;
    return;
  }

  TH3 * th3 = dynamic_cast<TH3*> (obj);
  if(th3) {
    th3->Fill(x, y, z) ;
    return;
  }
  
  AliError(Form("can not find histogram (of instance TH2) <%s> ",key)) ;
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskPi0Efficiency::TestLambda(Double_t pt,Double_t l1,Double_t l2){
  
  Double_t l2Mean  = 1.53126+9.50835e+06/(1.+1.08728e+07*pt+1.73420e+06*pt*pt) ;
  Double_t l1Mean  = 1.12365+0.123770*TMath::Exp(-pt*0.246551)+5.30000e-03*pt ;
  Double_t l2Sigma = 6.48260e-02+7.60261e+10/(1.+1.53012e+11*pt+5.01265e+05*pt*pt)+9.00000e-03*pt;
  Double_t l1Sigma = 4.44719e-04+6.99839e-01/(1.+1.22497e+00*pt+6.78604e-07*pt*pt)+9.00000e-03*pt;
  Double_t c=-0.35-0.550*TMath::Exp(-0.390730*pt) ;
  Double_t R2=0.5*(l1-l1Mean)*(l1-l1Mean)/l1Sigma/l1Sigma + 
              0.5*(l2-l2Mean)*(l2-l2Mean)/l2Sigma/l2Sigma +
              0.5*c*(l1-l1Mean)*(l2-l2Mean)/l1Sigma/l2Sigma ;
  return (R2<2.5*2.5) ;
  
}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskPi0Efficiency::TestLambda2(Double_t pt,Double_t l1,Double_t l2){
  
  Double_t l2Mean  = 1.53126+9.50835e+06/(1.+1.08728e+07*pt+1.73420e+06*pt*pt) ;
  Double_t l1Mean  = 1.12365+0.123770*TMath::Exp(-pt*0.246551)+5.30000e-03*pt ;
  Double_t l2Sigma = 6.48260e-02+7.60261e+10/(1.+1.53012e+11*pt+5.01265e+05*pt*pt)+9.00000e-03*pt;
  Double_t l1Sigma = 4.44719e-04+6.99839e-01/(1.+1.22497e+00*pt+6.78604e-07*pt*pt)+9.00000e-03*pt;
  Double_t c=-0.35-0.550*TMath::Exp(-0.390730*pt) ;
  Double_t R2=0.5*(l1-l1Mean)*(l1-l1Mean)/l1Sigma/l1Sigma + 
              0.5*(l2-l2Mean)*(l2-l2Mean)/l2Sigma/l2Sigma +
              0.5*c*(l1-l1Mean)*(l2-l2Mean)/l1Sigma/l2Sigma ;
  return (R2<1.5*1.5) ;
  
}
//___________________________________________________________________________
void AliAnalysisTaskPi0Efficiency::ProcessMC(){
  //fill histograms for efficiensy etc. calculation
  const Double_t rcut = 1. ; //cut for primary particles
  //---------First pi0/eta-----------------------------
  char partName[10] ;
  char hkey[55] ;

  AliAODEvent *event = dynamic_cast<AliAODEvent*>(InputEvent());
  if(!event) return ;
  TClonesArray *mcArray = (TClonesArray*)event->FindListObject(AliAODMCParticle::StdBranchName());
  for(Int_t i=0;i<mcArray->GetEntriesFast();i++){
     AliAODMCParticle* particle =  (AliAODMCParticle*) mcArray->At(i);
    if(particle->GetPdgCode() == 111)
      snprintf(partName,10,"pi0") ;
    else
      if(particle->GetPdgCode() == 221)
        snprintf(partName,10,"eta") ;
      else
        if(particle->GetPdgCode() == 22)
           snprintf(partName,10,"gamma") ;
	else
           continue ;

    //Primary particle
    Double_t r=TMath::Sqrt(particle->Xv()*particle->Xv()+particle->Yv()*particle->Yv());
    if(r >rcut)
      continue ;

    Double_t pt = particle->Pt() ;
    //Total number of pi0 with creation radius <1 cm
    snprintf(hkey,55,"hMC_all_%s_cen%d",partName,fCenBin) ;
    FillHistogram(hkey,pt) ;
    if(TMath::Abs(particle->Y())<0.12){
      snprintf(hkey,55,"hMC_unitEta_%s_cen%d",partName,fCenBin) ;
      FillHistogram(hkey,pt) ;
    }

    snprintf(hkey,55,"hMC_rap_%s_cen%d",partName,fCenBin) ;
    FillHistogram(hkey,particle->Y()) ;
    
    Double_t phi=particle->Phi() ;
    while(phi<0.)phi+=TMath::TwoPi() ;
    while(phi>TMath::TwoPi())phi-=TMath::TwoPi() ;
    snprintf(hkey,55,"hMC_phi_%s_cen%d",partName,fCenBin) ;
    FillHistogram(hkey,phi) ;
   

    //Check if one of photons converted
    if(particle->GetNDaughters()!=2)
     continue ; //Do not account Dalitz decays

/*
    TParticle * gamma1 = fStack->Particle(particle->GetFirstDaughter());
    TParticle * gamma2 = fStack->Particle(particle->GetLastDaughter());
    //Number of pi0s decayed into acceptance
    Int_t mod1,mod2 ;
    Double_t x=0.,z=0. ;
    Bool_t hitPHOS1 = fPHOSGeo->ImpactOnEmc(gamma1, mod1, z,x) ;
    Bool_t hitPHOS2 = fPHOSGeo->ImpactOnEmc(gamma2, mod2, z,x) ;

    Bool_t goodPair=kFALSE ;
    if( hitPHOS1 && hitPHOS2){
      sprintf(hkey,"hMC_PHOSacc_%s",partName) ;
      FillHistogram(hkey,pt) ;
      goodPair=kTRUE ;
    }

*/
  }
 
  //Now calculate "Real" distribution of clusters with primary
  TClonesArray cluPrim("AliCaloPhoton",200) ; //clusters with primary
  TClonesArray * clusters = (TClonesArray*)event->FindListObject("EmbeddedCaloClusters") ;
  AliAODCaloCells * cells = (AliAODCaloCells *)event->FindListObject("EmbeddedPHOScells") ;
  Int_t multClust = clusters->GetEntriesFast();
  Int_t inPHOS=0 ;
  Double_t vtx0[3] = {0,0,0}; 
  for (Int_t i=0; i<multClust; i++) {
    AliAODCaloCluster *clu = (AliAODCaloCluster*)clusters->At(i);
    if ( !clu->IsPHOS() || clu->E()<0.3) continue;
    if(clu->GetLabel()<0) continue ;

    Float_t  position[3];
    clu->GetPosition(position);
    TVector3 global(position) ;
    Int_t relId[4] ;
    fPHOSGeo->GlobalPos2RelId(global,relId) ;
    Int_t mod  = relId[0] ;
    Int_t cellX = relId[2];
    Int_t cellZ = relId[3] ;
    if ( !IsGoodChannel("PHOS",mod,cellX,cellZ) ) 
      continue ;
    if(clu->GetNCells()<3)
      continue ;

    TLorentzVector pv1 ;
    clu->GetMomentum(pv1 ,vtx0);
    
    if(inPHOS>=cluPrim.GetSize()){
      cluPrim.Expand(inPHOS+50) ;
    }
    AliCaloPhoton * ph = new(cluPrim[inPHOS]) AliCaloPhoton(pv1.X(),pv1.Py(),pv1.Z(),pv1.E()) ;
    //AliCaloPhoton * ph = (AliCaloPhoton*)fPHOSEvent->At(inPHOS) ;
    ph->SetModule(mod) ;
    AliPHOSAodCluster cluPHOS1(*clu);
    cluPHOS1.Recalibrate(fPHOSCalibData,cells); // modify the cell energies
    Double_t ecore=CoreEnergy(&cluPHOS1) ;
    pv1*= ecore/pv1.E() ;
    ph->SetMomV2(&pv1) ;
    ph->SetNCells(clu->GetNCells());
    ph->SetDispBit(TestLambda(clu->E(),clu->GetM20(),clu->GetM02())) ;
    ph->SetDisp2Bit(TestLambda2(clu->E(),clu->GetM20(),clu->GetM02())) ;
    ph->SetCPVBit(clu->GetEmcCpvDistance()>2.) ; //radius in sigmas
    ph->SetCPV2Bit(clu->GetEmcCpvDistance()>4.) ;
    ph->SetPhoton(clu->GetNExMax()<2); // Remember, if it is unfolded


    inPHOS++ ;

  }
  
  //Single photon
  char key[55] ;
  for (Int_t i1=0; i1<inPHOS; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)cluPrim.At(i1) ;
    snprintf(key,55,"hMCPhotAll_cen%d",fCenBin) ;
    FillHistogram(key,ph1->Pt()) ;
    snprintf(key,55,"hMCPhotAllcore_cen%d",fCenBin) ;
    FillHistogram(key,ph1->GetMomV2()->Pt()) ;
    if(ph1->IsPhoton()){
      snprintf(key,55,"hMCPhotAllwou_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt()) ;
    }
    if(ph1->IsCPVOK() ){
      snprintf(key,55,"hMCPhotCPV_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt()) ;
      snprintf(key,55,"hMCPhotCPVcore_cen%d",fCenBin) ;
      FillHistogram(key,ph1->GetMomV2()->Pt()) ;
      
    }
    if(ph1->IsCPV2OK() ){
      snprintf(key,55,"hMCPhotCPV2_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt()) ;
      
    }
    if(ph1->IsDisp2OK()){
      snprintf(key,55,"hMCPhotDisp2_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt()) ;
    }
    if(ph1->IsDispOK()){
      snprintf(key,55,"hMCPhotDisp_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt()) ;
      if(ph1->IsPhoton()){
        snprintf(key,55,"hMCPhotDispwou_cen%d",fCenBin) ;
        FillHistogram(key,ph1->Pt()) ;
      }
      if(ph1->IsCPVOK()){
	snprintf(key,55,"hMCPhotBoth_cen%d",fCenBin) ;
	FillHistogram(key,ph1->Pt()) ;
	snprintf(key,55,"hMCPhotBothcore_cen%d",fCenBin) ;
	FillHistogram(key,ph1->GetMomV2()->Pt()) ;
      }
    } // end of loop i2
  } // end of loop i1 

  // Fill Real disribution
  for (Int_t i1=0; i1<inPHOS-1; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)cluPrim.At(i1) ;
    for (Int_t i2=i1+1; i2<inPHOS; i2++) {
      AliCaloPhoton * ph2=(AliCaloPhoton*)cluPrim.At(i2) ;
      TLorentzVector p12  = *ph1  + *ph2;
      TLorentzVector pv12 = *(ph1->GetMomV2()) + *(ph2->GetMomV2());      
      Double_t a=TMath::Abs((ph1->E()-ph2->E())/(ph1->E()+ph2->E())) ;
       
      snprintf(key,55,"hMCMassPtAll_cen%d",fCenBin) ;
      FillHistogram(key,p12.M() ,p12.Pt()) ;
      snprintf(key,55,"hMCMassPtAllcore_cen%d",fCenBin) ;
      FillHistogram(key,pv12.M(), pv12.Pt()) ;
      if(ph1->IsPhoton()&& ph2->IsPhoton()){
        snprintf(key,55,"hMCMassPtAllwou_cen%d",fCenBin) ;
        FillHistogram(key,p12.M() ,p12.Pt()) ;
      }
      if(a<0.9){
        snprintf(key,55,"hMCMassPtAll_a09_cen%d",fCenBin) ;
        FillHistogram(key,p12.M() ,p12.Pt()) ;
        if(a<0.8){
          snprintf(key,55,"hMCMassPtAll_a08_cen%d",fCenBin) ;
          FillHistogram(key,p12.M() ,p12.Pt()) ;
          if(a<0.7){
            snprintf(key,55,"hMCMassPtAll_a07_cen%d",fCenBin) ;
            FillHistogram(key,p12.M() ,p12.Pt()) ;
          }
        }
      }

      
           if(ph1->Module()==1 && ph2->Module()==1)
	    FillHistogram("hMCPi0M11",p12.M(),p12.Pt() );
          else if(ph1->Module()==2 && ph2->Module()==2)
	    FillHistogram("hMCPi0M22",p12.M(),p12.Pt() );
          else if(ph1->Module()==3 && ph2->Module()==3)
	    FillHistogram("hMCPi0M33",p12.M(),p12.Pt() );
          else if(ph1->Module()==1 && ph2->Module()==2)
	    FillHistogram("hMCPi0M12",p12.M(),p12.Pt() );
          else if(ph1->Module()==1 && ph2->Module()==3)
	    FillHistogram("hMCPi0M13",p12.M(),p12.Pt() );
          else if(ph1->Module()==2 && ph2->Module()==3)
	    FillHistogram("hMCPi0M23",p12.M(),p12.Pt() );
	 

     
      
      
      if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	snprintf(key,55,"hMCMassPtCPV_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt()) ;
	snprintf(key,55,"hMCMassPtCPVcore_cen%d",fCenBin) ;
	FillHistogram(key,pv12.M(), pv12.Pt()) ;
        if(a<0.9){
          snprintf(key,55,"hMCMassPtCPV_a09_cen%d",fCenBin) ;
          FillHistogram(key,p12.M() ,p12.Pt()) ;
          if(a<0.8){
            snprintf(key,55,"hMCMassPtCPV_a08_cen%d",fCenBin) ;
            FillHistogram(key,p12.M() ,p12.Pt()) ;
            if(a<0.7){
              snprintf(key,55,"hMCMassPtCPV_a07_cen%d",fCenBin) ;
              FillHistogram(key,p12.M() ,p12.Pt()) ;
            }
          }
        }
      }
      if(ph1->IsCPV2OK() && ph2->IsCPV2OK()){
	snprintf(key,55,"hMCMassPtCPV2_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt()) ;
        if(a<0.9){
          snprintf(key,55,"hMCMassPtCPV2_a09_cen%d",fCenBin) ;
          FillHistogram(key,p12.M() ,p12.Pt()) ;
          if(a<0.8){
            snprintf(key,55,"hMCMassPtCPV2_a08_cen%d",fCenBin) ;
            FillHistogram(key,p12.M() ,p12.Pt()) ;
            if(a<0.7){
              snprintf(key,55,"hMCMassPtCPV2_a07_cen%d",fCenBin) ;
              FillHistogram(key,p12.M() ,p12.Pt()) ;
            }
          }
        }
      }
      if(ph1->IsDisp2OK() && ph2->IsDisp2OK()){
	snprintf(key,55,"hMCMassPtDisp2_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt()) ;
      }
      if(ph1->IsDispOK() && ph2->IsDispOK()){
	snprintf(key,55,"hMCMassPtDisp_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt()) ;
        if(ph1->IsPhoton()&& ph2->IsPhoton()){
  	  snprintf(key,55,"hMCMassPtDispwou_cen%d",fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt()) ;
	}
	if(a<0.9){
          snprintf(key,55,"hMCMassPtDisp_a09_cen%d",fCenBin) ;
          FillHistogram(key,p12.M() ,p12.Pt()) ;
          if(a<0.8){
            snprintf(key,55,"hMCMassPtDisp_a08_cen%d",fCenBin) ;
            FillHistogram(key,p12.M() ,p12.Pt()) ;
            if(a<0.7){
              snprintf(key,55,"hMCMassPtDisp_a07_cen%d",fCenBin) ;
              FillHistogram(key,p12.M() ,p12.Pt()) ;
            }
          }
        }

	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  snprintf(key,55,"hMCMassPtBoth_cen%d",fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt()) ;
	  snprintf(key,55,"hMCMassPtBothcore_cen%d",fCenBin) ;
	  FillHistogram(key,pv12.M(), pv12.Pt()) ;
          if(a<0.9){
            snprintf(key,55,"hMCMassPtBoth_a09_cen%d",fCenBin) ;
            FillHistogram(key,p12.M() ,p12.Pt()) ;
            if(a<0.8){
              snprintf(key,55,"hMCMassPtBoth_a08_cen%d",fCenBin) ;
              FillHistogram(key,p12.M() ,p12.Pt()) ;
              if(a<0.7){
                snprintf(key,55,"hMCMassPtBoth_a07_cen%d",fCenBin) ;
                FillHistogram(key,p12.M() ,p12.Pt()) ;
              }
            }
          }
        }
      }
    } // end of loop i2
  } // end of loop i1
}
//____________________________________________________________________________
Double_t  AliAnalysisTaskPi0Efficiency::CoreEnergy(AliPHOSAodCluster * clu){  
  //calculate energy of the cluster in the circle with radius distanceCut around the maximum
  
  //Can not use already calculated coordinates?
  //They have incidence correction...
  const Double_t distanceCut =3.5 ;
  const Double_t logWeight=4.5 ;
  
  Double32_t * elist = clu->GetCellsAmplitudeFraction() ;  
// Calculates the center of gravity in the local PHOS-module coordinates
  Float_t wtot = 0;
  Double_t xc[100]={0} ;
  Double_t zc[100]={0} ;
  Double_t x = 0 ;
  Double_t z = 0 ;
  Int_t mulDigit=TMath::Min(100,clu->GetNCells()) ;
  for(Int_t iDigit=0; iDigit<mulDigit; iDigit++) {
    Int_t relid[4] ;
    Float_t xi ;
    Float_t zi ;
    fPHOSGeo->AbsToRelNumbering(clu->GetCellAbsId(iDigit), relid) ;
    fPHOSGeo->RelPosInModule(relid, xi, zi);
    xc[iDigit]=xi ;
    zc[iDigit]=zi ;
    if (clu->E()>0 && elist[iDigit]>0) {
      Float_t w = TMath::Max( 0., logWeight + TMath::Log( elist[iDigit] / clu->E() ) ) ;
      x    += xc[iDigit] * w ;
      z    += zc[iDigit] * w ;
      wtot += w ;
    }
  }
  if (wtot>0) {
    x /= wtot ;
    z /= wtot ;
  }
  Double_t coreE=0. ;
  for(Int_t iDigit=0; iDigit < mulDigit; iDigit++) {
    Double_t distance = TMath::Sqrt((xc[iDigit]-x)*(xc[iDigit]-x)+(zc[iDigit]-z)*(zc[iDigit]-z)) ;
    if(distance < distanceCut)
      coreE += elist[iDigit] ;
  }
  //Apply non-linearity correction
  return (0.0241+1.0504*coreE+0.000249*coreE*coreE) ;
}

  
