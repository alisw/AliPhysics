#include "TChain.h"
#include "TTree.h"
#include "TObjArray.h"
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
#include "TROOT.h"
#include "THashList.h"

#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliAnalysisTaskSE.h"
#include "AliPHOSHijingEfficiency.h"
#include "AliCaloPhoton.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSEsdCluster.h"
#include "AliPHOSCalibData.h"
#include "AliESDEvent.h"
#include "AliESDCaloCells.h"
#include "AliESDVertex.h"
#include "AliESDtrackCuts.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliCDBManager.h"
#include "AliCentrality.h" 
#include "AliESDtrackCuts.h"
#include "AliEventplane.h"
#include "TProfile.h"
#include "AliOADBContainer.h"

// Calculates pi0 and photon efficiencies using Hijing events
// Also calculates pi0 contamination: non-vertex pi0s
// And photon contamination:
//    - off-vertex photons
//    - misidentified photons
// Authors: Dmitri Peressounko
// Date   : 05.02.2013

ClassImp(AliPHOSHijingEfficiency)

//________________________________________________________________________
Double_t rnlin(Double_t x)
{
  //This is standard PHOS non-linearity used in reconstruction
  //and re-calibration by 2.2%
  return 1.022*(0.0241+1.0504*x+0.000249*x*x) ; 

}

//________________________________________________________________________
AliPHOSHijingEfficiency::AliPHOSHijingEfficiency(const char *name) 
: AliAnalysisTaskSE(name),
  fStack(0x0),
  fOutputContainer(0x0),
  fPHOSEvent(0x0),
  fPHOSCalibData(0x0),
  fRP(0.),fRPV0A(0.),fRPV0C(0.),fHaveTPCRP(0),
  fRunNumber(0),
  fCentrality(0.),
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
  DefineOutput(1,TList::Class());

  // Set bad channel map
  char key[55] ;
  for(Int_t i=0; i<6; i++){
    sprintf(key,"PHOS_BadMap_mod%d",i) ;
    fPHOSBadMap[i]=new TH2I(key,"Bad Modules map",64,0.,64.,56,0.,56.) ;
  }


}

//________________________________________________________________________
void AliPHOSHijingEfficiency::UserCreateOutputObjects()
{


  // Create histograms
  // Called once
  const Int_t nRuns=200 ;
  
  // ESD histograms
  if(fOutputContainer != NULL){
    delete fOutputContainer;
  }
  fOutputContainer = new THashList();
  fOutputContainer->SetOwner(kTRUE);
  PostData(1, fOutputContainer);
  
  //========QA histograms=======

  //Event selection
  fOutputContainer->Add(new TH2F("hSelEvents","Event selection", 10,0.,10.,nRuns,0.,float(nRuns))) ;
  fOutputContainer->Add(new TH1F("hTotSelEvents","Event selection", 10,0.,10.)) ;
  
  //vertex distribution
  fOutputContainer->Add(new TH2F("hZvertex","Z vertex position", 50,-25.,25.,nRuns,0.,float(nRuns))) ;
  
  //Centrality
  fOutputContainer->Add(new TH2F("hCentrality","Event centrality", 100,0.,100.,nRuns,0.,float(nRuns))) ;
  fOutputContainer->Add(new TH2F("hCenPHOS","Centrality vs PHOSclusters", 100,0.,100.,200,0.,200.)) ;
  fOutputContainer->Add(new TH2F("hCenPHOSCells","Centrality vs PHOS cells", 100,0.,100.,100,0.,1000.)) ;
  fOutputContainer->Add(new TH2F("hCenTrack","Centrality vs tracks", 100,0.,100.,100,0.,15000.)) ;  
  fOutputContainer->Add(new TH2F("hCluEvsClu","ClusterMult vs E",200,0.,20.,100,0.,100.)) ;
  fOutputContainer->Add(new TH2F("hCluEvsCluM","ClusterMult vs E",200,0.,20.,100,0.,20.)) ;
   			
  fOutputContainer->Add(new TH3F("hCPVr","CPV radius",100,0.,20.,100,0.,20.,10,0.,100.));
  
  //Bad Map
  fOutputContainer->Add(new TH2F("hCluLowM1","Cell (X,Z), M1" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluLowM2","Cell (X,Z), M2" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluLowM3","Cell (X,Z), M3" ,64,0.5,64.5, 56,0.5,56.5));

  fOutputContainer->Add(new TH2F("hCluHighM1","Cell (X,Z), M1" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluHighM2","Cell (X,Z), M2" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluHighM3","Cell (X,Z), M3" ,64,0.5,64.5, 56,0.5,56.5));
  
  fOutputContainer->Add(new TH2F("hCluVetoM1","Cell (X,Z), M1" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluVetoM2","Cell (X,Z), M2" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluVetoM3","Cell (X,Z), M3" ,64,0.5,64.5, 56,0.5,56.5));

  fOutputContainer->Add(new TH2F("hCluDispM1","Cell (X,Z), M1" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluDispM2","Cell (X,Z), M2" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluDispM3","Cell (X,Z), M3" ,64,0.5,64.5, 56,0.5,56.5));

  //Single photon and pi0 spectrum
  Int_t nPtPhot = 300 ;
  Double_t ptPhotMax = 30 ;
  Int_t nM       = 500;
  Double_t mMin  = 0.0;
  Double_t mMax  = 1.0;
    
  //PHOS calibration QA
  fOutputContainer->Add(new TH2F("hPi0M11","Pairs in modules",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hPi0M12","Pairs in modules",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hPi0M13","Pairs in modules",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hPi0M22","Pairs in modules",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hPi0M23","Pairs in modules",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hPi0M33","Pairs in modules",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    
  char key[55] ;
  for(Int_t cent=0; cent<6; cent++){

    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad2_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad4_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad6_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));

    fOutputContainer->Add(new TH1F(Form("hPhotAll_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAllcore_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAllwou_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotDisp_cen%d",cent),"Disp clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotDisp2_cen%d",cent),"Disp clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotDispcore_cen%d",cent),"Disp clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotDisp2core_cen%d",cent),"Disp clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotDispwou_cen%d",cent),"Disp clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotCPV_cen%d",cent),"CPV clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotCPVcore_cen%d",cent),"CPV clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotCPV2_cen%d",cent),"CPV clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotBoth_cen%d",cent),"Both clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotBothcore_cen%d",cent),"Both clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotBoth2_cen%d",cent),"Both2 clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotBoth2core_cen%d",cent),"Both2 clusters",nPtPhot,0.,ptPhotMax));


    fOutputContainer->Add(new TH2F(Form("hPi0All_cen%d",cent),"All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH2F(Form("hPi0Allcore_cen%d",cent),"All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH2F(Form("hPi0Allwou_cen%d",cent),"All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH2F(Form("hPi0Disp_cen%d",cent),"Disp clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH2F(Form("hPi0Disp2_cen%d",cent),"Disp clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH2F(Form("hPi0Dispcore_cen%d",cent),"Disp clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH2F(Form("hPi0Disp2core_cen%d",cent),"Disp clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH2F(Form("hPi0Dispwou_cen%d",cent),"Disp clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH2F(Form("hPi0CPV_cen%d",cent),"CPV clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH2F(Form("hPi0CPVcore_cen%d",cent),"CPV clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH2F(Form("hPi0CPV2_cen%d",cent),"CPV clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH2F(Form("hPi0Both_cen%d",cent),"Both clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH2F(Form("hPi0Bothcore_cen%d",cent),"Both clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH2F(Form("hPi0Both2_cen%d",cent),"Both2 clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH2F(Form("hPi0Both2core_cen%d",cent),"Both2 clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));

    
    snprintf(key,55,"hPi0All_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hPi0Disp_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Disp clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hPi0CPV_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"CPV clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hPi0CPV2_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"CPV clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hPi0Both_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Both clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));        
    snprintf(key,55,"hPi0Both2_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Both2 clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));        

    snprintf(key,55,"hSingleAll_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hSingleAllcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hSingleAllwou_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hSingleDisp_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Disp clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hSingleDisp2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Disp clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hSingleDispcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Disp clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hSingleDisp2core_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Disp clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hSingleDispwou_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Disp clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hSingleCPV_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"CPV clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hSingleCPVcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"CPV clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hSingleCPV2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"CPV clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hSingleBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Both clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hSingleBothcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Both clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hSingleBoth2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Both2 clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hSingleBoth2core_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Both2 clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    
    
    snprintf(key,55,"hMiPi0All_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiPi0Allcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiPi0Allwou_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiPi0Disp_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Disp clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiPi0Disp2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Disp clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiPi0Dispwou_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Disp clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiPi0Dispcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Disp clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiPi0Disp2core_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Disp clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiPi0CPV_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"CPV clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiPi0CPVcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"CPV clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiPi0CPV2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"CPV clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiPi0Both_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Both clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiPi0Bothcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Both clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiPi0Both2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Both2 clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiPi0Both2core_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Both2 clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    
    snprintf(key,55,"hMiPi0All_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiPi0Disp_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Disp clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiPi0CPV_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"CPV clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiPi0CPV2_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"CPV clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiPi0Both_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Both clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));   
    snprintf(key,55,"hMiPi0Both2_a07_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Both2 clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));   

    snprintf(key,55,"hMiSingleAll_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiSingleAllwou_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiSingleAllcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiSingleDisp_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Disp clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiSingleDisp2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Disp clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiSingleDispwou_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Disp clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiSingleDispcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Disp clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiSingleDisp2core_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Disp clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiSingleCPV_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"CPV clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiSingleCPVcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"CPV clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiSingleCPV2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"CPV clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiSingleBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Both clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiSingleBothcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Both clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiSingleBoth2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Both2 clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    snprintf(key,55,"hMiSingleBoth2core_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"Both2 clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  }

 
 //MC
  for(Int_t cent=0; cent<6; cent++){
    snprintf(key,55,"hMC_rap_gamma_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Rapidity pi0",200,-1.,1.)) ;
    snprintf(key,55,"hMC_rap_pi0_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Rapidity pi0",200,-1.,1.)) ;
    snprintf(key,55,"hMC_rap_eta_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"Rapidity eta",200,-1.,1.)) ;
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
  fOutputContainer->Add(new TH2F("hMC_gamma_vertex","Creation vertex",25,0.,25.,1000,0.,500.)) ;
  fOutputContainer->Add(new TH2F("hMC_pi0_vertex","Creation vertex",25,0.,25.,1000,0.,500.)) ;
  fOutputContainer->Add(new TH2F("hMC_eta_vertex","Creation vertex",25,0.,25.,1000,0.,500.)) ;
 
 
  Int_t nPt      = 200;
  Double_t ptMin = 0;
  Double_t ptMax = 20; 
  fOutputContainer->Add(new TH2F("Vertex","Pi0 creation vertex",nPt,ptMin,ptMax,5000,0.,500.));
  fOutputContainer->Add(new TH3F("hSecondPi0RphiZ","Secondary pi0 vertex",450,0.,450.,100,0.,TMath::TwoPi(),200,-100.,100.));
  fOutputContainer->Add(new TH2F("hSecondPi0RE","Secondary pi0 vertex",450,0.,450.,200,0.,20.));
  fOutputContainer->Add(new TH3F("hMass_R","Mass vs radius any parent",50,0.,0.25,100,0.,10.,300,0.,600.));
  fOutputContainer->Add(new TH3F("Real_pi_R","All clusters",50,0.,0.25,100,0.,10.,250,0.,500.));
  fOutputContainer->Add(new TH3F("Real_pi_Z","All clusters",50,0.,0.25,100,0.,10.,100,-100.,100.));
//  fOutputContainer->Add(new TH2F(Form("Real_npnp_RZ"),"All clusters",250,0.,500.,100,-100.,100.));
//  fOutputContainer->Add(new TH3F(Form("Real_mass_R"),"All clusters",50,0.,0.25,100,0.,10.,300,0.,600.));

  for(Int_t cen=0; cen<6; cen++){
    fOutputContainer->Add(new TH1F(Form("hPrimPhot_cen%d",cen),"Primary spetrum",nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH1F(Form("hPrimEl_cen%d",cen),"Primary spetrum",nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH1F(Form("hPrimPi0_cen%d",cen),"Primary spetrum",nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH1F(Form("hPrimEta_cen%d",cen),"Primary spetrum",nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH1F(Form("hPrimPipm_cen%d",cen),"Primary spetrum",nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH1F(Form("hPrimP_cen%d",cen),"Primary spetrum",nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH1F(Form("hPrimPbar_cen%d",cen),"Primary spetrum",nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH1F(Form("hPrimN_cen%d",cen),"Primary spetrum",nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH1F(Form("hPrimNbar_cen%d",cen),"Primary spetrum",nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH1F(Form("hPrimK0S_cen%d",cen),"Primary spetrum",nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH1F(Form("hPrimK0L_cen%d",cen),"Primary spetrum",nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH1F(Form("hPrimKpm_cen%d",cen),"Primary spetrum",nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH1F(Form("hPrimOther_cen%d",cen),"Primary spetrum",nPt,ptMin,ptMax));

    //pairs from common parents
    fOutputContainer->Add(new TH2F(Form("hParentAll_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentK0s_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentGamma_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentEl_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentOther_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentPi0_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentDirPi0_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
   
    //common parent - pi0
    fOutputContainer->Add(new TH2F(Form("hParentPi0NoPrim_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentPi0Eta_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentPi0Omega_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentPi0Pipm_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentPi0Kpm_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentPi0Ks_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentPi0Kl_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentPi0pn_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH2F(Form("hParentPi0antipn_cen%d",cen),"All clusters",nM,mMin,mMax,nPt,ptMin,ptMax));
    
  }
  
  
  //Photon contaminations
  fOutputContainer->Add(new TH2F("hPipmGammaConvR","Conversion radius" ,200,0.,20.,1000,0.,500.));
  fOutputContainer->Add(new TH2F("hPipmElConvR","Conversion radius" ,200,0.,20.,1000,0.,500.));
  fOutputContainer->Add(new TH2F("hPipmNConvR","Conversion radius" ,200,0.,20.,1000,0.,500.));
  fOutputContainer->Add(new TH2F("hPipmOtherConvR","Conversion radius" ,200,0.,20.,1000,0.,500.));
  fOutputContainer->Add(new TH2F("hPipmGammaConvRZ","Conversion radius" ,400,-200.,200.,1000,0.,500.)); 
   
   const Int_t nTypes=24 ;
   char partTypes[nTypes][55] ;
   snprintf(partTypes[0],55,"hGammaNoPrim") ; //
   snprintf(partTypes[1],55,"hGammaPhot") ; //
   snprintf(partTypes[2],55,"hGammaEl") ; //
   snprintf(partTypes[3],55,"hGammaPi0") ; //
   snprintf(partTypes[4],55,"hGammaEta") ; //
   snprintf(partTypes[5],55,"hhGammaOmega") ; //
   snprintf(partTypes[6],55,"hGammaPipm") ; //
   snprintf(partTypes[7],55,"hGammaP") ; //
   snprintf(partTypes[8],55,"hGammaPbar") ; //
   snprintf(partTypes[9],55,"hGammaN") ; //
   snprintf(partTypes[10],55,"hGammaNbar") ; //
   snprintf(partTypes[11],55,"hGammaK0S") ; //
   snprintf(partTypes[12],55,"hGammaK0L") ; //
   snprintf(partTypes[13],55,"hGammaKpm") ; //
   snprintf(partTypes[14],55,"hGammaKstar") ; //
   snprintf(partTypes[15],55,"hGammaDelta") ; //
   snprintf(partTypes[16],55,"hGammaOtherCharged") ; //
   snprintf(partTypes[17],55,"hGammaOtherNeutral") ; //
   snprintf(partTypes[18],55,"hGammaPipmGamma") ; //
   snprintf(partTypes[19],55,"hGammaPipmEl") ; //
   snprintf(partTypes[20],55,"hGammaPipmOther") ; //
   snprintf(partTypes[21],55,"hGammaPipmDirect") ; //
   snprintf(partTypes[22],55,"hGammaPipmp") ; //
   snprintf(partTypes[23],55,"hGammaPipmn") ; //
 
   const Int_t nPID=12 ;
   char cPID[25][12] ;
   snprintf(cPID[0],25,"All") ;
   snprintf(cPID[1],25,"Allcore") ;
   snprintf(cPID[2],25,"CPV") ;
   snprintf(cPID[3],25,"CPVcore") ;
   snprintf(cPID[4],25,"CPV2") ;
   snprintf(cPID[5],25,"CPV2core") ;
   snprintf(cPID[6],25,"Disp") ;
   snprintf(cPID[7],25,"Dispcore") ;
   snprintf(cPID[8],25,"Disp2") ;
   snprintf(cPID[9],25,"Disp2core") ;
   snprintf(cPID[10],25,"Both") ;
   snprintf(cPID[11],25,"Bothcore") ;
 
   for(Int_t itype=0; itype<nTypes; itype++){
     for(Int_t iPID=0; iPID<nPID; iPID++){
       for(Int_t cen=0; cen<5; cen++){
         fOutputContainer->Add(new TH1F(Form("%s_%s_cen%d",partTypes[itype],cPID[iPID],cen),"Cluster parents",nPt,ptMin,ptMax));
       }
     }
   }
 
  
 
  PostData(1, fOutputContainer);

  printf("Init done \n") ;
  
}

//________________________________________________________________________
void AliPHOSHijingEfficiency::UserExec(Option_t *) 
{
  // Main loop, called for each event
  // Analyze ESD/AOD
  fStack=0 ;
  if(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler()){
    if(static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent())
      fStack = static_cast<AliMCEventHandler*>(AliAnalysisManager::GetAnalysisManager()->GetMCtruthEventHandler())->MCEvent()->Stack();
  }
  
  
  AliESDEvent *event = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!event) {
    Printf("ERROR: Could not retrieve event");
    PostData(1, fOutputContainer);
    return;
  }
 
  fRunNumber=ConvertRunNumber(event->GetRunNumber()) ;
  FillHistogram("hSelEvents",0.5,fRunNumber-0.5) ;
  FillHistogram("hTotSelEvents",0.5) ;
  
  // Get PHOS rotation matrices from ESD and set them to the PHOS geometry
  char skey[55] ;
  if(fEventCounter == 0) {
    // Initialize the PHOS geometry
    fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP") ;

    //We have to apply re-calibration for pass1 LCH10h
    // Initialize decalibration factors in the form of the OCDB object
    const Float_t wDecalib = 0.055;

    fPHOSCalibData = new AliPHOSCalibData();
    fPHOSCalibData->SetName("PhosCalibData");

    for(Int_t module=1; module<=5; module++) {
      for(Int_t column=1; column<=56; column++) {
        for(Int_t row=1; row<=64; row++) {
          Float_t adcChannelEmc = gRandom->Gaus(1.,wDecalib);
          fPHOSCalibData->SetADCchannelEmc(module,column,row,adcChannelEmc);
        }
      }
    }
    for(Int_t mod=0; mod<5; mod++) {
      if(!event->GetPHOSMatrix(mod)) continue;
      fPHOSGeo->SetMisalMatrix(event->GetPHOSMatrix(mod),mod) ;
      Printf("PHOS geo matrix %p for module # %d is set\n", event->GetPHOSMatrix(mod), mod);
     }
    
    fEventCounter++ ;
  }


  // Checks if we have a primary vertex
  // Get primary vertices form ESD
  const AliESDVertex *esdVertex5 = event->GetPrimaryVertex();

  Double_t vtx0[3] = {0,0,0}; // don't rely on ESD vertex, assume (0,0,0)
  Double_t vtx5[3] ={0.,0.,0.};
  
  vtx5[0] = esdVertex5->GetX();
  vtx5[1] = esdVertex5->GetY();
  vtx5[2] = esdVertex5->GetZ();
  
  
  FillHistogram("hZvertex",esdVertex5->GetZ(),fRunNumber-0.5);
  if (TMath::Abs(esdVertex5->GetZ()) > 10. ){
    PostData(1, fOutputContainer);
    return;
  }
  FillHistogram("hSelEvents",1.5,fRunNumber-0.5) ;
  FillHistogram("hTotSelEvents",1.5) ;

/*  
  if(event->IsPileupFromSPD()){
    PostData(1, fOutputContainer);
    return;
  } 
*/
  FillHistogram("hSelEvents",2.5,fRunNumber-0.5) ;
  FillHistogram("hTotSelEvents",2.5) ;  
  
  
  /*
  //Vtx class z-bin
  Int_t zvtx = (Int_t)((vtx5[2]+10.)/2.) ;
  if(zvtx<0)zvtx=0 ;
  if(zvtx>9)zvtx=9 ;
  */
  //No dependence on zVtx observed, save memory
  Int_t zvtx=0 ;

  AliCentrality *centrality = event->GetCentrality(); 
  fCentrality=centrality->GetCentralityPercentile("V0M");

  if( fCentrality <= 0. || fCentrality>80. ){
    PostData(1, fOutputContainer);
    return;
  }

  FillHistogram("hSelEvents",3.5,fRunNumber-0.5) ;
  FillHistogram("hTotSelEvents",3.5) ;

  if(fCentrality<5.)
    fCenBin=0 ;
  else if(fCentrality<10.)
    fCenBin=1 ;
  else if(fCentrality<20.)
    fCenBin=2 ;
  else if(fCentrality<40.)
    fCenBin=3 ;
  else if(fCentrality<60.)
    fCenBin=4 ;
  else if(fCentrality<80.)
    fCenBin=5 ;
  else{
    PostData(1, fOutputContainer);
    return;
  }

  fRP=0. ;

  FillHistogram("hSelEvents",4.5,fRunNumber-0.5) ;
  FillHistogram("hTotSelEvents",4.5) ;
  //All event selections done
  FillHistogram("hCentrality",fCentrality,fRunNumber-0.5) ;
  //Reaction plane is defined in the range (0;pi)
  //We have 10 bins
  Double_t averageRP = 0.;//fRPV0A+fRPV0C+fRP ;
  Int_t irp=Int_t(10.*(averageRP)/TMath::Pi());
  if(irp>9)irp=9 ;

  if(!fPHOSEvents[zvtx][fCenBin][irp]) 
    fPHOSEvents[zvtx][fCenBin][irp]=new TList() ;
  TList * prevPHOS = fPHOSEvents[zvtx][fCenBin][irp] ;

  ProcessMC() ;

  if(fPHOSEvent)
    fPHOSEvent->Clear() ;
  else
    fPHOSEvent = new TClonesArray("AliCaloPhoton",200) ;

  //For re-calibration
  const Double_t logWeight=4.5 ;  

  TVector3 vertex(vtx5);
  
  Int_t multClust = event->GetNumberOfCaloClusters();
  Int_t inPHOS=0,inMod1=0, inMod2=0, inMod3=0 ;
  Double_t avrgEm1=0,avrgEm2=0,avrgEm3=0; //average cluster energy

  AliESDCaloCells * cells = event->GetPHOSCells() ;
  FillHistogram("hCenPHOSCells",fCentrality,cells->GetNumberOfCells()) ;
  FillHistogram("hCenTrack",fCentrality,event->GetNumberOfTracks()) ;
  

  TVector3 localPos ;
  for (Int_t i=0; i<multClust; i++) {
    AliESDCaloCluster *clu = event->GetCaloCluster(i);
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
    
    
    FillHistogram("hCluEvsClu",clu->E(),clu->GetNCells()) ;
    FillHistogram("hCluEvsCluM",clu->E(),clu->GetM02()) ;

    //Apply re-Calibreation
    AliPHOSEsdCluster cluPHOS1(*clu);
    cluPHOS1.Recalibrate(fPHOSCalibData,cells); // modify the cell energies
    Reclusterize(&cluPHOS1) ;
    cluPHOS1.EvalAll(logWeight,vertex);         // recalculate the cluster parameters
    cluPHOS1.SetE(rnlin(cluPHOS1.E()));// Users's nonlinearity
    if(cluPHOS1.E()<0.3) continue;
    if(clu->GetNCells()<3) continue ;
    if(clu->GetM02()<0.2)   continue ;
           
 
    if(mod==1){
      avrgEm1+=cluPHOS1.E() ;
      inMod1++ ;
    }
    if(mod==2){
      avrgEm2+=cluPHOS1.E() ;
      inMod2++ ;
    }
    if(mod==3){
      avrgEm3+=cluPHOS1.E() ;
      inMod3++ ;
    }
       
    TLorentzVector pv1 ;
    cluPHOS1.GetMomentum(pv1 ,vtx0);
    
    Double_t ecore=CoreEnergy(&cluPHOS1) ; 
    
    sprintf(skey,"hCluLowM%d",mod) ;
    FillHistogram(skey,cellX,cellZ,1.);
    if(cluPHOS1.E()>1.5){
      sprintf(skey,"hCluHighM%d",mod) ;
      FillHistogram(skey,cellX,cellZ,1.);
    }
    
    if(inPHOS>=fPHOSEvent->GetSize()){
      fPHOSEvent->Expand(inPHOS+50) ;
    }
    new((*fPHOSEvent)[inPHOS]) AliCaloPhoton(pv1.X(),pv1.Py(),pv1.Z(),pv1.E()) ;
    AliCaloPhoton * ph = (AliCaloPhoton*)fPHOSEvent->At(inPHOS) ;
    ph->SetModule(mod) ;
    pv1*= ecore/pv1.E() ;
    ph->SetMomV2(&pv1) ;
    ph->SetNCells(clu->GetNCells());    
    Double_t m02=0.,m20=0.;
    EvalLambdas(&cluPHOS1,0,m02, m20);   
    ph->SetDispBit(TestLambda(clu->E(),m20,m02)) ;
    EvalLambdas(&cluPHOS1,1,m02, m20);
    ph->SetDisp2Bit(TestLambda2(clu->E(),m20,m02)) ;

    Double_t distBC=clu->GetDistanceToBadChannel();
    if(distBC>2.)
      FillHistogram(Form("hPhotAll_DistBad2_cen%d",fCenBin),ph->Pt()) ;
      if(distBC>4.)
        FillHistogram(Form("hPhotAll_DistBad4_cen%d",fCenBin),ph->Pt()) ;
        if(distBC>6.)
          FillHistogram(Form("hPhotAll_DistBad6_cen%d",fCenBin),ph->Pt()) ;
    
    if(ph->IsDispOK()){
      sprintf(skey,"hCluDispM%d",mod) ;
      FillHistogram(skey,cellX,cellZ,1.);
    }
    //Track matching
    Double_t dx=clu->GetTrackDx() ;
    Double_t dz=clu->GetTrackDz() ;

    Bool_t cpvBit=kTRUE ; //No track matched by default
    Bool_t cpvBit2=kTRUE ; //More Strict criterion
    TArrayI * itracks = clu->GetTracksMatched() ;  
    if(itracks->GetSize()>0){
      Int_t iTr = itracks->At(0);
      if(iTr>=0 && iTr<event->GetNumberOfTracks()){
        AliESDtrack *track = event->GetTrack(iTr) ;
        Double_t pt = track->Pt() ;
        Short_t charge = track->Charge() ;
        Double_t r=TestCPV(dx, dz, pt,charge) ;
	cpvBit=(r>2.) ;
	cpvBit2=(r>4.) ;
        FillHistogram("hCPVr",r,pv1.E(),fCentrality);
      }
    }
    ph->SetCPVBit(cpvBit) ;
    ph->SetCPV2Bit(cpvBit2) ;
    if(cpvBit){
      sprintf(skey,"hCluVetoM%d",mod) ;
      FillHistogram(skey,cellX,cellZ,1.);
    }
    ph->SetEMCx(float(cellX)) ;
    ph->SetEMCz(float(cellZ)) ;
    ph->SetLambdas(cluPHOS1.GetM20(),cluPHOS1.GetM02()) ;
    ph->SetUnfolded(clu->GetNExMax()<2); // Remember, if it is unfolded   
    Bool_t sure=  kTRUE;
    Int_t primary=FindPrimary(clu,sure) ;  //номер праймари частицы в стеке
    ph->SetPrimary(primary) ;
    ph->SetWeight(PrimaryWeight(primary)) ;
    inPHOS++ ;
  }
  
  FillHistogram("hCenPHOS",fCentrality,inPHOS) ;
  if(inPHOS==0){
    PostData(1, fOutputContainer);
    fEventCounter++;
    return ; 
  }
  
   
  for (Int_t i1=0; i1<inPHOS; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1) ;
    Double_t w=ph1->GetWeight() ;
    Double_t pt = ph1->Pt() ;
    Double_t ptv =ph1->GetMomV2()->Pt() ;
            
    FillHistogram(Form("hPhotAll_cen%d",fCenBin),pt,w) ;
    FillHistogram(Form("hPhotAllcore_cen%d",fCenBin),ptv,w) ;
    if(ph1->IsntUnfolded()){
      FillHistogram(Form("hPhotAllwou_cen%d",fCenBin),pt,w) ;
    }
    if(ph1->IsCPVOK()){
      FillHistogram(Form("hPhotCPV_cen%d",fCenBin),pt,w) ;
      FillHistogram(Form("hPhotCPVcore_cen%d",fCenBin),ptv,w) ;
    }
    if(ph1->IsCPV2OK()){
      FillHistogram(Form("hPhotCPV2_cen%d",fCenBin),pt,w) ;
    }
    if(ph1->IsDispOK()){      
      FillHistogram(Form("hPhotDisp_cen%d",fCenBin),pt,w) ;

      FillHistogram(Form("hPhotDispcore_cen%d",fCenBin),ptv,w) ;
      if(ph1->IsntUnfolded()){
        FillHistogram(Form("hPhotDispwou_cen%d",fCenBin),pt,w) ;
      }
      if(ph1->IsCPVOK()){
	FillHistogram(Form("hPhotBoth_cen%d",fCenBin),pt,w) ;
	FillHistogram(Form("hPhotBothcore_cen%d",fCenBin),ptv,w) ;
      }
    }  
    if(ph1->IsDisp2OK()){
      FillHistogram(Form("hPhotDisp2_cen%d",fCenBin),pt,w) ;
      FillHistogram(Form("hPhotDisp2core_cen%d",fCenBin),ptv,w) ;
      if(ph1->IsCPVOK()){
	FillHistogram(Form("hPhotBoth2_cen%d",fCenBin),pt,w) ;
	FillHistogram(Form("hPhotBoth2core_cen%d",fCenBin),ptv,w) ;
      }
    }
  }
 
  
  const Double_t alphaCut=0.7 ;
  //pi0
  for (Int_t i1=0; i1<inPHOS-1; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1) ;
    Double_t w1=ph1->GetWeight() ;
    for (Int_t i2=i1+1; i2<inPHOS; i2++) {
      AliCaloPhoton * ph2=(AliCaloPhoton*)fPHOSEvent->At(i2) ;
      Double_t w2=ph2->GetWeight() ;
      Double_t w = TMath::Sqrt(w1*w2) ;
      TLorentzVector p12  = *ph1  + *ph2;
      TLorentzVector pv12 = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
            
      Double_t a=TMath::Abs((ph1->E()-ph2->E())/(ph1->E()+ph2->E())) ;
      Double_t m = p12.M() ;
      Double_t pt = p12.Pt() ;
      Double_t mv = p12.M() ;
      Double_t ptv = p12.Pt() ;

      FillHistogram(Form("hPi0All_cen%d",fCenBin),pt,m,w) ;
      FillHistogram(Form("hPi0Allcore_cen%d",fCenBin),ptv,mv,w) ;
      if(ph1->IsntUnfolded() && ph2->IsntUnfolded()){
        FillHistogram(Form("hPi0Allwou_cen%d",fCenBin),pt,m,w) ;
      }
      
      
      Double_t pt1 = ph1->Pt() ;      
      Double_t ptv1 = ph1->GetMomV2()->Pt() ;      
      Double_t pt2 = ph2->Pt() ;      
      Double_t ptv2 = ph2->GetMomV2()->Pt() ;      
      FillHistogram(Form("hSingleAll_cen%d",fCenBin),m,pt1,w1) ;
      FillHistogram(Form("hSingleAll_cen%d",fCenBin),m,pt2,w2) ;
      FillHistogram(Form("hSingleAllcore_cen%d",fCenBin),mv,ptv1,w1) ;
      FillHistogram(Form("hSingleAllcore_cen%d",fCenBin),mv,ptv2,w2) ;
      if(ph1->IsntUnfolded())
        FillHistogram(Form("hSingleAllwou_cen%d",fCenBin),m,pt1,w1) ;
      if(ph2->IsntUnfolded())
        FillHistogram(Form("hSingleAllwou_cen%d",fCenBin),m,pt2,w2) ;
      if(ph1->IsCPVOK()){
        FillHistogram(Form("hSingleCPV_cen%d",fCenBin),m,pt1,w1) ;
        FillHistogram(Form("hSingleCPVcore_cen%d",fCenBin),mv,ptv1,w1) ;
      }
      if(ph2->IsCPVOK()){
        FillHistogram(Form("hSingleCPV_cen%d",fCenBin),m,pt2,w2) ;
        FillHistogram(Form("hSingleCPVcore_cen%d",fCenBin),mv,ptv2,w2) ;
      }
      if(ph1->IsCPV2OK()){
        FillHistogram(Form("hSingleCPV2_cen%d",fCenBin),m,pt1,w1) ;
      }
      if(ph2->IsCPV2OK()){
        FillHistogram(Form("hSingleCPV2_cen%d",fCenBin),m,pt2,w2) ;
      }      
      if(ph1->IsDispOK()){
        FillHistogram(Form("hSingleDisp_cen%d",fCenBin),m,pt1,w1) ;
        if(ph1->IsntUnfolded()){
          FillHistogram(Form("hSingleDispwou_cen%d",fCenBin),m,pt1,w1) ;
	}
        FillHistogram(Form("hSingleDispcore_cen%d",fCenBin),mv,ptv1,w1) ;
      }
      if(ph2->IsDispOK()){
        FillHistogram(Form("hSingleDisp_cen%d",fCenBin),m,pt2,w2) ;
        if(ph1->IsntUnfolded()){
          FillHistogram(Form("hSingleDispwou_cen%d",fCenBin),m,pt2,w2) ;
	}
        FillHistogram(Form("hSingleDispcore_cen%d",fCenBin),mv,ptv2,w2) ;
      }
      if(ph1->IsDisp2OK()){
        FillHistogram(Form("hSingleDisp2_cen%d",fCenBin),m,pt1,w1) ;
        FillHistogram(Form("hSingleDisp2core_cen%d",fCenBin),mv,ptv1,w1) ;
      }
      if(ph2->IsDisp2OK()){
        FillHistogram(Form("hSingleDisp2_cen%d",fCenBin),m,pt2,w2) ;
        FillHistogram(Form("hSingleDisp2core_cen%d",fCenBin),mv,ptv2,w2) ;
      }
      if(ph1->IsDispOK() && ph1->IsCPVOK()){
        FillHistogram(Form("hSingleBoth_cen%d",fCenBin),m,pt1,w1) ;
        FillHistogram(Form("hSingleBothcore_cen%d",fCenBin),mv,ptv1,w1) ;
      }
      if(ph2->IsDispOK() && ph2->IsCPVOK()){
        FillHistogram(Form("hSingleBoth_cen%d",fCenBin),m,pt2,w2) ;
        FillHistogram(Form("hSingleBothcore_cen%d",fCenBin),mv,ptv2,w2) ;
      }            
      if(ph1->IsDisp2OK() && ph1->IsCPVOK()){
        FillHistogram(Form("hSingleBoth2_cen%d",fCenBin),m,pt1,w1) ;
        FillHistogram(Form("hSingleBoth2core_cen%d",fCenBin),mv,ptv1,w1) ;
      }
      if(ph2->IsDisp2OK() && ph2->IsCPVOK()){
        FillHistogram(Form("hSingleBoth2_cen%d",fCenBin),m,pt2,w2) ;
        FillHistogram(Form("hSingleBoth2core_cen%d",fCenBin),mv,ptv2,w2) ;
      }            
      
      
      if(a<alphaCut){
        FillHistogram(Form("hPi0All_a07_cen%d",fCenBin),pt,m,w) ;
      }

      if(ph1->IsCPVOK() && ph2->IsCPVOK()){ //pi0 CPV	
	FillHistogram(Form("hPi0CPV_cen%d",fCenBin),pt,m,w) ;
	FillHistogram(Form("hPi0CPVcore_cen%d",fCenBin),ptv,mv,w) ;
		
        if(a<alphaCut){
          FillHistogram(Form("hPi0CPV_a07_cen%d",fCenBin),pt,m,w) ;
        }
      } 
      if(ph1->IsCPV2OK() && ph2->IsCPV2OK()){
	FillHistogram(Form("hPi0CPV2_cen%d",fCenBin),p12.M(),p12.Pt(),w) ;
        if(a<alphaCut){
          FillHistogram(Form("hPi0CPV2_a07_cen%d",fCenBin),pt,m,w) ;
        }
      } 
      if(ph1->IsDisp2OK() && ph2->IsDisp2OK()){
	FillHistogram(Form("hPi0Disp2_cen%d",fCenBin),pt,m,w) ;
	FillHistogram(Form("hPi0Disp2core_cen%d",fCenBin),ptv,mv,w) ;
	if(ph1->IsCPVOK() && ph2->IsCPVOK()){ //pi0 Both
	  FillHistogram(Form("hPi0Both2_cen%d",fCenBin),pt,m,w) ;
	  FillHistogram(Form("hPi0Both2core_cen%d",fCenBin),ptv,mv,w) ;
          if(a<alphaCut){
            FillHistogram(Form("hPi0Both2_a07_cen%d",fCenBin),pt,m,w) ;
          }
 	}
      }
      
      
      if(ph1->IsDispOK() && ph2->IsDispOK()){ //pi0 Disp	
	FillHistogram(Form("hPi0Disp_cen%d",fCenBin),pt,m,w) ;
	FillHistogram(Form("hPi0Dispcore_cen%d",fCenBin),ptv,mv,w) ;		
        if(a<alphaCut){
          FillHistogram(Form("hPi0Disp_a07_cen%d",fCenBin),pt,m,w) ;
        }
	
	if(ph1->IsCPVOK() && ph2->IsCPVOK()){ //pi0 Both	  
	  FillHistogram(Form("hPi0Both_cen%d",fCenBin),pt,m,w) ;
	  FillHistogram(Form("hPi0Bothcore_cen%d",fCenBin),ptv,mv,w) ;
	  
          if(a<alphaCut){
            FillHistogram(Form("hPi0Both_a07_cen%d",fCenBin),pt,m,w) ;
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
        Double_t m = p12.M() ;
        Double_t pt = p12.Pt() ;
        Double_t mv = p12.M() ;
        Double_t ptv = p12.Pt() ;
	
	FillHistogram(Form("hMiPi0All_cen%d",fCenBin),pt,m) ;
	FillHistogram(Form("hMiPi0Allcore_cen%d",fCenBin),ptv,mv) ;
	if(ph1->IsntUnfolded() && ph2->IsntUnfolded()){
	  FillHistogram(Form("hMiPi0Allwou_cen%d",fCenBin),pt,m) ;
	}

        Double_t pt1 = ph1->Pt() ;      
        Double_t ptv1 = ph1->GetMomV2()->Pt() ;      
        Double_t pt2 = ph2->Pt() ;      
        Double_t ptv2 = ph2->GetMomV2()->Pt() ;      
        FillHistogram(Form("hMiSingleAll_cen%d",fCenBin),m,pt1) ;
        FillHistogram(Form("hMiSingleAll_cen%d",fCenBin),m,pt2) ;
        FillHistogram(Form("hMiSingleAllcore_cen%d",fCenBin),mv,ptv1) ;
        FillHistogram(Form("hMiSingleAllcore_cen%d",fCenBin),mv,ptv2) ;
        if(ph1->IsntUnfolded())
          FillHistogram(Form("hMiSingleAllwou_cen%d",fCenBin),m,pt1) ;
        if(ph2->IsntUnfolded())
          FillHistogram(Form("hMiSingleAllwou_cen%d",fCenBin),m,pt2) ;
        if(ph1->IsCPVOK()){
          FillHistogram(Form("hMiSingleCPV_cen%d",fCenBin),m,pt1) ;
          FillHistogram(Form("hMiSingleCPVcore_cen%d",fCenBin),mv,ptv1) ;
        }
        if(ph2->IsCPVOK()){
          FillHistogram(Form("hMiSingleCPV_cen%d",fCenBin),m,pt2) ;
          FillHistogram(Form("hMiSingleCPVcore_cen%d",fCenBin),mv,ptv2) ;
        }
        if(ph1->IsCPV2OK()){
          FillHistogram(Form("hMiSingleCPV2_cen%d",fCenBin),m,pt1) ;
        }
        if(ph2->IsCPV2OK()){
          FillHistogram(Form("hMiSingleCPV2_cen%d",fCenBin),m,pt2) ;
        }      
        if(ph1->IsDispOK()){
          FillHistogram(Form("hMiSingleDisp_cen%d",fCenBin),m,pt1) ;
          if(ph1->IsntUnfolded()){
            FillHistogram(Form("hMiSingleDispwou_cen%d",fCenBin),m,pt1) ;
	  }
          FillHistogram(Form("hMiSingleDispcore_cen%d",fCenBin),mv,ptv1) ;
        }
        if(ph2->IsDispOK()){
          FillHistogram(Form("hMiSingleDisp_cen%d",fCenBin),m,pt2) ;
          if(ph1->IsntUnfolded()){
            FillHistogram(Form("hMiSingleDispwou_cen%d",fCenBin),m,pt2) ;
	  }
          FillHistogram(Form("hMiSingleDispcore_cen%d",fCenBin),mv,ptv2) ;
        }
        if(ph1->IsDisp2OK()){
          FillHistogram(Form("hMiSingleDisp2_cen%d",fCenBin),m,pt1) ;
          FillHistogram(Form("hMiSingleDisp2core_cen%d",fCenBin),mv,ptv1) ;
        } 
        if(ph2->IsDisp2OK()){
          FillHistogram(Form("hMiSingleDisp2_cen%d",fCenBin),m,pt2) ;
          FillHistogram(Form("hMiSingleDisp2core_cen%d",fCenBin),mv,ptv2) ;
        }
        if(ph1->IsDispOK() && ph1->IsCPVOK()){
          FillHistogram(Form("hMiSingleBoth_cen%d",fCenBin),m,pt1) ;
          FillHistogram(Form("hMiSingleBothcore_cen%d",fCenBin),mv,ptv1) ;
        }
        if(ph2->IsDispOK() && ph2->IsCPVOK()){
          FillHistogram(Form("hMiSingleBoth_cen%d",fCenBin),m,pt2) ;
          FillHistogram(Form("hMiSingleBothcore_cen%d",fCenBin),mv,ptv2) ;
        }            
        if(ph1->IsDisp2OK() && ph1->IsCPVOK()){
          FillHistogram(Form("hMiSingleBoth2_cen%d",fCenBin),m,pt1) ;
          FillHistogram(Form("hMiSingleBoth2core_cen%d",fCenBin),mv,ptv1) ;
        }
        if(ph2->IsDisp2OK() && ph2->IsCPVOK()){
          FillHistogram(Form("hMiSingleBoth2_cen%d",fCenBin),m,pt2) ;
          FillHistogram(Form("hMiSingleBoth2core_cen%d",fCenBin),mv,ptv2) ;
        }            
	
	
	
        if(a<alphaCut){
          FillHistogram(Form("hMiPi0All_a07_cen%d",fCenBin),pt,m) ;
        }
	if(ph1->IsCPVOK() && ph2->IsCPVOK()){	  
	  FillHistogram(Form("hMiPi0CPV_cen%d",fCenBin),p12.M(), p12.Pt()) ;
	  FillHistogram(Form("hMiPi0CPVcore_cen%d",fCenBin),ptv,mv) ;

	  if(a<alphaCut){
            FillHistogram(Form("hMiPi0CPV_a07_cen%d",fCenBin),pt,m) ;
          }
	}
	if(ph1->IsCPV2OK() && ph2->IsCPV2OK()){
	  FillHistogram(Form("hMiPi0CPV2_cen%d",fCenBin),pt,m) ;

	  if(a<alphaCut){
            FillHistogram(Form("hMiPi0CPV2_a07_cen%d",fCenBin),pt,m) ;
          }
	}
	if(ph1->IsDisp2OK() && ph2->IsDisp2OK()){
	  FillHistogram(Form("hMiPi0Disp2_cen%d",fCenBin),pt,m) ;
	  FillHistogram(Form("hMiPi0Disp2core_cen%d",fCenBin),ptv,mv) ;
	  if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	    FillHistogram(Form("hMiPi0Both2_cen%d",fCenBin),pt,m) ;
	    FillHistogram(Form("hMiPi0Both2core_cen%d",fCenBin),pv12.M(),pv12.Pt()) ;	    
	  }
	}	    
	if(ph1->IsDispOK() && ph2->IsDispOK()){
	  
	  FillHistogram(Form("hMiPi0Disp_cen%d",fCenBin),pt,m) ;
	  FillHistogram(Form("hMiPi0Dispcore_cen%d",fCenBin),pv12.M(),pv12.Pt()) ;
          if(ph1->IsntUnfolded() && ph2->IsntUnfolded()){
	    FillHistogram(Form("hMiPi0Dispwou_cen%d",fCenBin),p12.M(),p12.Pt()) ;
	  }
	  
	  if(a<alphaCut){
            FillHistogram(Form("hMiPi0Disp_a07_cen%d",fCenBin),pt,m) ;
          }	  
	  if(ph1->IsCPVOK() && ph2->IsCPVOK()){	    	    	    
	    FillHistogram(Form("hMiPi0Both_cen%d",fCenBin),pt,m) ;
	    FillHistogram(Form("hMiPi0Bothcore_cen%d",fCenBin),pv12.M(),pv12.Pt()) ;

	    if(a<alphaCut){
              FillHistogram(Form("hMiPi0Both_a07_cen%d",fCenBin),pt,m) ;
            }
	  }
	}
      } // end of loop i2
    }
  } // end of loop i1

  
  FillSecondaries() ;
  
  
  //Now we either add current events to stack or remove
  //If no photons in current event - no need to add it to mixed
  const Int_t kMixEvents[6]={5,5,5,10,10,30} ;
  if(fPHOSEvent->GetEntriesFast()>0){
    prevPHOS->AddFirst(fPHOSEvent) ;
    fPHOSEvent=0;
    if(prevPHOS->GetSize()>kMixEvents[fCenBin]){//Remove redundant events
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
void AliPHOSHijingEfficiency::FillSecondaries(){
  //Sort secondaires
  
  const Double_t rcut = 1. ;
  //Fill spectra of primary particles 
  //with proper weight
  for(Int_t i=0; i<fStack->GetNtrack(); i++){
    TParticle * p = fStack->Particle(i) ;
    if(p->R()>rcut)
      continue ;
    if(TMath::Abs(p->Y())>0.5)
      continue ;
    Double_t w = PrimaryParticleWeight(p) ;  
    Int_t primPdgCode=p->GetPdgCode() ;
      switch(primPdgCode){
	case  22: FillHistogram(Form("hPrimPhot_cen%d",fCenBin),p->Pt(),w); 
	          break ;
	case  11: 
	case -11: 
	          FillHistogram(Form("hPrimEl_cen%d",fCenBin),p->Pt(),w); 
	          break ;
	case  111: 
	          FillHistogram(Form("hPrimPi0_cen%d",fCenBin),p->Pt(),w); 
	          break ;
	case  221: 
	          FillHistogram(Form("hPrimEta_cen%d",fCenBin),p->Pt(),w); 
	          break ;
	case  211: 
	case -211: 
	          FillHistogram(Form("hPrimPipm_cen%d",fCenBin),p->Pt(),w); 
	          break ;		  
	case  2212:  //p 
	          FillHistogram(Form("hPrimP_cen%d",fCenBin),p->Pt(),w); 
	          break ;		  
	case -2212:  //pbar
	          FillHistogram(Form("hPrimPbar_cen%d",fCenBin),p->Pt(),w); 
	          break ;		  
	case  2112:  //n 
	          FillHistogram(Form("hPrimN_cen%d",fCenBin),p->Pt(),w); 
	          break ;		  
	case -2112:  //nbar
	          FillHistogram(Form("hPrimNbar_cen%d",fCenBin),p->Pt(),w); 
	          break ;
	case  310:  //nbar
	          FillHistogram(Form("hPrimK0S_cen%d",fCenBin),p->Pt(),w); 
	          break ;
	case  130:  //nbar
	          FillHistogram(Form("hPrimK0L_cen%d",fCenBin),p->Pt(),w); 
	          break ;
	case  321:  //K+
	case -321:  //K-
	          FillHistogram(Form("hPrimKpm_cen%d",fCenBin),p->Pt(),w); 
	          break ;
	default:	   //other
	          FillHistogram(Form("hPrimOther_cen%d",fCenBin),p->Pt(),w);    
      }
  }

  //Origins of secondary pi0s
  for(Int_t i=0; i<fStack->GetNtrack(); i++){
    TParticle * p = fStack->Particle(i) ;
    if(p->GetPdgCode()!=111)
      continue ;
    FillHistogram("Vertex",p->Pt(),p->R());
    if(p->R()<rcut)
      continue ;
    Double_t phi=p->Phi() ;
    while(phi<0.)phi+=TMath::TwoPi() ;
    while(phi>TMath::TwoPi())phi-=TMath::TwoPi() ;
    FillHistogram("hSecondPi0RphiZ",p->R(),phi,p->Vz()) ;   
    Double_t w = PrimaryParticleWeight(p) ;  
    FillHistogram("hSecondPi0RE",p->R(),p->Pt(),w) ;   
  }

  TLorentzVector p1;

  Int_t inPHOS=fPHOSEvent->GetEntries() ;
  for (Int_t i1=0; i1<inPHOS-1; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1) ;
    Double_t w1=ph1->GetWeight() ;
    for (Int_t i2=i1+1; i2<inPHOS; i2++) {
      AliCaloPhoton * ph2=(AliCaloPhoton*)fPHOSEvent->At(i2) ;
      TLorentzVector p12  = *ph1  + *ph2;
      Double_t w2=ph2->GetWeight() ;
      Double_t w = TMath::Sqrt(w1*w2) ;
      FillHistogram(Form("hParentAll_cen%d",fCenBin),p12.M(),p12.Pt(),w) ;
      Int_t prim=FindCommonParent(ph1->GetPrimary(),ph2->GetPrimary()) ;
      if(prim>-1){
        TParticle * particle = (TParticle *)fStack->Particle(prim);
        FillHistogram("hMass_R",p12.M(),p12.Pt(),TMath::Sqrt(particle->R()*particle->R()+particle->Vz()*particle->Vz())) ;
		
	
        Int_t pdgCode=particle->GetPdgCode() ;
        if(pdgCode!=111){ //common parent - not pi0
          if(pdgCode==22)  
            FillHistogram(Form("hParentGamma_cen%d",fCenBin),p12.M(),p12.Pt(),w) ;
	  else{		    
            if(pdgCode==11 || pdgCode==-11){   
              FillHistogram(Form("hParentEl_cen%d",fCenBin),p12.M(),p12.Pt(),w) ;
	    }
  	    else{
              if(InPi0mass(p12.M() ,p12.Pt())){
	        printf("Common parent: %d \n",pdgCode) ;
	      }
              FillHistogram(Form("hParentOther_cen%d",fCenBin),p12.M(),p12.Pt(),w) ;
	    }
	  }//Not photons
        }//Common parent not pi0
        else{ //common parent - pi0
          FillHistogram(Form("hParentPi0_cen%d",fCenBin),p12.M(),p12.Pt(),w) ;	
          FillHistogram(Form("Real_pi_R"),p12.M(),p12.Pt(),particle->R(),w) ;	
          FillHistogram(Form("Real_pi_Z"),p12.M(),p12.Pt(),particle->Vz(),w) ;	
	  if(particle->R()<rcut && TMath::Abs(particle->Vz())<10.){
            FillHistogram(Form("hParentDirPi0_cen%d",fCenBin),p12.M(),p12.Pt(),w) ;
	    continue ;
	  }
	  //Common particle pi0, created off-vertex
  	  Int_t primPi0=particle->GetFirstMother();
	  if(primPi0==-1){
            FillHistogram(Form("hParentPi0NoPrim_cen%d",fCenBin),p12.M(),p12.Pt(),w) ;
	  }
	  else{
    	    Int_t primPdgCode=fStack->Particle(primPi0)->GetPdgCode() ;
            switch(primPdgCode){
            case 221: FillHistogram(Form("hParentPi0Eta_cen%d",fCenBin),p12.M(),p12.Pt(),w) ; //eta
	              break ;
            case 223: FillHistogram(Form("hParentPi0Omega_cen%d",fCenBin),p12.M(),p12.Pt(),w) ; //omega
	              break ;
	    case  211:  //pi+-
	    case -211: FillHistogram(Form("hParentPi0Pipm_cen%d",fCenBin),p12.M(),p12.Pt(),w) ; //
	              break ;
	    case  321:  //K+-
	    case -321: FillHistogram(Form("hParentPi0Kpm_cen%d",fCenBin),p12.M(),p12.Pt(),w) ; //
	              break ;
	    case 310: FillHistogram(Form("hParentPi0Ks_cen%d",fCenBin),p12.M(),p12.Pt(),w) ; // K0S
	              break ;
	    case 130: FillHistogram(Form("hParentPi0Kl_cen%d",fCenBin),p12.M(),p12.Pt(),w) ; // K0L
	              break ;
	    case  2212:  //p 
	    case  2112:  //n 
		      FillHistogram(Form("hParentPi0pn_cen%d",fCenBin),p12.M(),p12.Pt(),w) ; // pn
	              break ;
	    case -2212:  //pbar
	    case -2112:  //nbar
		      FillHistogram(Form("hParentPi0antipn_cen%d",fCenBin),p12.M(),p12.Pt(),w) ; // pn
	              break ;
	    default:	   //other
		      FillHistogram(Form("hParentPi0Other_cen%d",fCenBin),p12.M(),p12.Pt(),w) ; //
	    }//switch	  
          }//pi0 with primary
        }//common parent - pi0
      }//there is common primary 
    }//seond photon loop
  }//first photon loop
  
  
  //Now look at photon contaiminations
  for (Int_t i1=0; i1<inPHOS-1; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1) ;
    Int_t iprim = ph1->GetPrimary() ;
    if(iprim<0)
      FillAllHistograms(Form("hGammaNoPrim_cen%d",fCenBin),ph1) ; //
    else{
      //Find primary at vertex
      TParticle * primPHOS = fStack->Particle(iprim) ;
      Int_t iprimV=primPHOS->GetFirstMother();
      TParticle * primVtx = primPHOS ;
      while((iprimV>-1) && primVtx->R()>rcut){
	primVtx = fStack->Particle(iprimV) ;
        iprimV=primVtx->GetFirstMother();
      }
    
      //photon
      Int_t primPdgCode=primVtx->GetPdgCode() ;
      switch(primPdgCode){
	case  22: FillAllHistograms("hGammaPhot",ph1); 
	          break ;
	case  11: 
	case -11: 
	          FillAllHistograms("hGammaEl",ph1); 
	          break ;
	case  111: 
	          FillAllHistograms("hGammaPi0",ph1); 
	          break ;
	case  221: 
	          FillAllHistograms("hGammaEta",ph1); 
	          break ;
        case 223: FillAllHistograms("hGammaOmega",ph1) ; //omega
	          break ;
	case  211: 
	case -211: 
	          FillAllHistograms("hGammaPipm",ph1); 
		  //Find particle entered PHOS
		  if(primVtx == primPHOS)
	            FillAllHistograms("hGammaPipmDirect",ph1); 
		  else{
                    Int_t primPdgPHOS=primPHOS->GetPdgCode() ;
		    if(primPdgPHOS==22){
	               FillAllHistograms("hGammaPipmGamma",ph1); 
		       FillHistogram("hPipmGammaConvR",ph1->Pt(),primPHOS->R());
 		       FillHistogram("hPipmGammaConvRZ",primPHOS->Vz(),primPHOS->R());
 	               break ;		  
		    }
		    if(TMath::Abs(primPdgPHOS)==11){
	               FillAllHistograms("hGammaPipmEl",ph1); 
		       FillHistogram("hPipmElConvR",ph1->Pt(),primPHOS->R());
	               break ;		  
		    }
		    if(TMath::Abs(primPdgPHOS)==2212){
	               FillAllHistograms("hGammaPipmp",ph1); 
		       FillHistogram("hPipmNConvR",ph1->Pt(),primPHOS->R());
	               break ;		  
		    }
		    if(TMath::Abs(primPdgPHOS)==2112){
	               FillAllHistograms("hGammaPipmn",ph1); 
		       FillHistogram("hPipmNConvR",ph1->Pt(),primPHOS->R());
	               break ;		  
		    }
	            FillAllHistograms("hGammaPipmOther",ph1); 
		    FillHistogram("hPipmOtherConvR",ph1->Pt(),primPHOS->R());		    
		  }
	          break ;		  
	case  2212:  //p 
	          FillAllHistograms("hGammaP",ph1); 
	          break ;		  
	case -2212:  //pbar
	          FillAllHistograms("hGammaPbar",ph1); 
	          break ;		  
	case  2112:  //n 
	          FillAllHistograms("hGammaN",ph1); 
	          break ;		  
	case -2112:  //nbar
		  FillAllHistograms("hGammaNbar",ph1) ; // pn
	          break ;
	case  310:  //nbar
		  FillAllHistograms("hGammaK0S",ph1) ; // pn
	          break ;
	case  130:  //nbar
		  FillAllHistograms("hGammaK0L",ph1) ; // pn
	          break ;
	case  321:  //K+
	case -321:  //K-
		  FillAllHistograms("hGammaKpm",ph1) ; // pn
	          break ;
        case -323: 
        case  323: 
        case -313: 
        case  313: FillAllHistograms("hGammaKstar",ph1) ; // K*(892)
	          break ;
		  
	case -2224 : //Deltas
	case  2224 : //Deltas
	case -2214 : //Deltas
	case  2214 : //Deltas
	case -2114 : //Deltas
	case  2114 : //Deltas
	case -1114 : //Deltas
	case  1114 : //Deltas
	          FillAllHistograms("hGammaDelta",ph1) ; // pn
	          break ;		  
	default:	   //other
	    if(primVtx->GetPDG()->Charge())
	      FillAllHistograms("hGammaOtherCharged",ph1) ; //
            else
	      FillAllHistograms("hGammaOtherNeutral",ph1) ; //
      }
    }
  
  }//single photons
 
}
//________________________________________________________________________
void AliPHOSHijingEfficiency::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
  
}

//________________________________________________________________________
Bool_t AliPHOSHijingEfficiency::IsGoodChannel(const char * det, Int_t mod, Int_t ix, Int_t iz)
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
void AliPHOSHijingEfficiency::FillAllHistograms(const char * key,AliCaloPhoton * ph)const{
  //Fill All PID histograms
        
  Double_t w=ph->GetWeight() ;
  Double_t pt = ph->Pt() ;
  Double_t ptC=ph->GetMomV2()->Pt() ;
  FillHistogram(Form("%s_All_cen%d",key,fCenBin),pt,w) ;
  FillHistogram(Form("%s_Allcore_cen%d",key,fCenBin),ptC,w) ;
  if(ph->IsCPVOK()){
    FillHistogram(Form("%s_CPV_cen%d",key,fCenBin),pt,w) ;
    FillHistogram(Form("%s_CPVcore_cen%d",key,fCenBin),ptC,w) ;
  }
  if(ph->IsCPV2OK()){
    FillHistogram(Form("%s_CPV2_cen%d",key,fCenBin),pt,w) ;
    FillHistogram(Form("%s_CPV2core_cen%d",key,fCenBin),ptC,w) ;
  }
  if(ph->IsDispOK()){     
    FillHistogram(Form("%s_Disp_cen%d",key,fCenBin),pt,w) ;
    FillHistogram(Form("%s_Dispcore_cen%d",key,fCenBin),ptC,w) ;
    if(ph->IsDisp2OK()){
      FillHistogram(Form("%s_Disp2_cen%d",key,fCenBin),pt,w) ;
      FillHistogram(Form("%s_Disp2core_cen%d",key,fCenBin),ptC,w) ;
    }
    if(ph->IsCPVOK()){
      FillHistogram(Form("%s_Both_cen%d",key,fCenBin),pt,w) ;
      FillHistogram(Form("%s_Bothcore_cen%d",key,fCenBin),ptC,w) ;
    }
  }  
}
//_____________________________________________________________________________
void AliPHOSHijingEfficiency::FillHistogram(const char * key,Double_t x)const{
  //FillHistogram
  TObject * tmp = fOutputContainer->FindObject(key) ;
  if(!tmp){
    AliInfo(Form("can not find histogram <%s> ",key)) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH1I")){
    ((TH1I*)tmp)->Fill(x) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH1F")){
    ((TH1F*)tmp)->Fill(x) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH1D")){
    ((TH1D*)tmp)->Fill(x) ;
    return ;
  }  
  AliInfo(Form("can not find 1D histogram <%s> ",key)) ;
}
//_____________________________________________________________________________
void AliPHOSHijingEfficiency::FillHistogram(const char * key,Double_t x,Double_t y)const{
  //FillHistogram
  TObject * tmp = fOutputContainer->FindObject(key) ;
  if(!tmp){
    AliInfo(Form("can not find histogram <%s> ",key)) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH1F")){
    ((TH1F*)tmp)->Fill(x,y) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH2F")){
    ((TH2F*)tmp)->Fill(x,y) ;
    return ;
  }
  AliError(Form("Calling FillHistogram with 2 parameters for histo <%s> of type %s",key,tmp->IsA()->GetName())) ;
}

//_____________________________________________________________________________
void AliPHOSHijingEfficiency::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z) const{
  //Fills 1D histograms with key
  TObject * tmp = fOutputContainer->FindObject(key) ;
  if(!tmp){
    AliInfo(Form("can not find histogram <%s> ",key)) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH2F")){
    ((TH2F*)tmp)->Fill(x,y,z) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH3F")){
    ((TH3F*)tmp)->Fill(x,y,z) ;
    return ;
  }
}
//_____________________________________________________________________________
void AliPHOSHijingEfficiency::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z, Double_t w) const{
  //Fills 1D histograms with key
  TObject * tmp = fOutputContainer->FindObject(key) ;
  if(!tmp){
    AliInfo(Form("can not find histogram <%s> ",key)) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH3F")){
    ((TH3F*)tmp)->Fill(x,y,z,w) ;
    return ;
  }
}

//___________________________________________________________________________
Int_t AliPHOSHijingEfficiency::ConvertRunNumber(Int_t run){

  switch(run){	
  case  139517 : return 137; 
  case  139514 : return 136; 
  case  139513 : return 135; 
  case  139511 : return 134; 
  case  139510 : return 133; 
  case  139507 : return 132; 
  case  139505 : return 131; 
  case  139504 : return 130; 
  case  139503 : return 129; 
  case  139470 : return 128; 
  case  139467 : return 127; 
  case  139466 : return 126; 
  case  139465 : return 125; 
  case  139440 : return 124; 
  case  139439 : return 123; 
  case  139438 : return 122; 
  case  139437 : return 121; 
  case  139360 : return 120; 
  case  139329 : return 119; 
  case  139328 : return 118; 
  case  139314 : return 117; 
  case  139311 : return 116; 
  case  139310 : return 115; 
  case  139309 : return 114; 
  case  139308 : return 113; 
  case  139173 : return 112; 
  case  139172 : return 111; 
  case  139110 : return 110; 
  case  139107 : return 109; 
  case  139105 : return 108; 
  case  139104 : return 107; 
  case  139042 : return 106; 
  case  139038 : return 105; 
  case  139037 : return 104; 
  case  139036 : return 103; 
  case  139029 : return 102; 
  case  139028 : return 101; 
  case  138983 : return 100; 
  case  138982 : return 99; 
  case  138980 : return 98; 
  case  138979 : return 97; 
  case  138978 : return 96; 
  case  138977 : return 95; 
  case  138976 : return 94; 
  case  138973 : return 93; 
  case  138972 : return 92; 
  case  138965 : return 91; 
  case  138924 : return 90; 
  case  138872 : return 89; 
  case  138871 : return 88; 
  case  138870 : return 87; 
  case  138837 : return 86; 
  case  138830 : return 85; 
  case  138828 : return 84; 
  case  138826 : return 83; 
  case  138796 : return 82; 
  case  138795 : return 81; 
  case  138742 : return 80; 
  case  138732 : return 79; 
  case  138730 : return 78; 
  case  138666 : return 77; 
  case  138662 : return 76; 
  case  138653 : return 75; 
  case  138652 : return 74; 
  case  138638 : return 73; 
  case  138624 : return 72; 
  case  138621 : return 71; 
  case  138583 : return 70; 
  case  138582 : return 69; 
  case  138579 : return 68; 
  case  138578 : return 67; 
  case  138534 : return 66; 
  case  138469 : return 65; 
  case  138442 : return 64; 
  case  138439 : return 63; 
  case  138438 : return 62; 
  case  138396 : return 61; 
  case  138364 : return 60; 
  case  138359 : return 59; 
  case  138275 : return 58; 
  case  138225 : return 57; 
  case  138201 : return 56; 
  case  138200 : return 55; 
  case  138197 : return 54; 
  case  138192 : return 53; 
  case  138190 : return 52; 
  case  138154 : return 51; 
  case  138153 : return 50; 
  case  138151 : return 49; 
  case  138150 : return 48; 
  case  138126 : return 47; 
  case  138125 : return 46; 
  case  137848 : return 45; 
  case  137847 : return 44; 
  case  137844 : return 43; 
  case  137843 : return 42; 
  case  137752 : return 41; 
  case  137751 : return 40; 
  case  137748 : return 39; 
  case  137724 : return 38; 
  case  137722 : return 37; 
  case  137718 : return 36; 
  case  137704 : return 35; 
  case  137693 : return 34; 
  case  137692 : return 33; 
  case  137691 : return 32; 
  case  137689 : return 31; 
  case  137686 : return 30; 
  case  137685 : return 29; 
  case  137639 : return 28; 
  case  137638 : return 27; 
  case  137608 : return 26; 
  case  137595 : return 25; 
  case  137549 : return 24; 
  case  137546 : return 23; 
  case  137544 : return 22; 
  case  137541 : return 21; 
  case  137539 : return 20; 
  case  137531 : return 19; 
  case  137530 : return 18; 
  case  137443 : return 17; 
  case  137441 : return 16; 
  case  137440 : return 15; 
  case  137439 : return 14; 
  case  137434 : return 13; 
  case  137432 : return 12; 
  case  137431 : return 11; 
  case  137430 : return 10; 
  case  137366 : return 9; 
  case  137243 : return 8; 
  case  137236 : return 7; 
  case  137235 : return 6; 
  case  137232 : return 5; 
  case  137231 : return 4; 
  case  137165 : return 3; 
  case  137162 : return 2; 
  case  137161 : return 1;
  default : return 199;
  } 

}
//_____________________________________________________________________________
Bool_t AliPHOSHijingEfficiency::TestLambda(Double_t pt,Double_t l1,Double_t l2){
  
  //For R=3.5
   Double_t  l1Mean  = 1.170014 -0.059465/(1.+0.019343*pt+0.147544*pt*pt) ;
   Double_t  l2Mean = 1.626270 + 0.761554*exp(-1.213839*pt)-0.020027*pt ;
   Double_t  l1Sigma = 0.133409 + 0.261307*exp(-0.636874*pt)-0.002849*pt ;
   Double_t  l2Sigma = 0.289698 + 0.459400*exp(-1.214242*pt)-0.012578*pt ;
   Double_t  c=-0.124103 ;
/*  
  Double_t l2Mean  = 1.53126+9.50835e+06/(1.+1.08728e+07*pt+1.73420e+06*pt*pt) ;
  Double_t l1Mean  = 1.12365+0.123770*TMath::Exp(-pt*0.246551)+5.30000e-03*pt ;
  Double_t l2Sigma = 6.48260e-02+7.60261e+10/(1.+1.53012e+11*pt+5.01265e+05*pt*pt)+9.00000e-03*pt;
  Double_t l1Sigma = 4.44719e-04+6.99839e-01/(1.+1.22497e+00*pt+6.78604e-07*pt*pt)+9.00000e-03*pt;
  Double_t c=-0.35-0.550*TMath::Exp(-0.390730*pt) ;
*/
  Double_t R2=0.5*(l1-l1Mean)*(l1-l1Mean)/l1Sigma/l1Sigma + 
              0.5*(l2-l2Mean)*(l2-l2Mean)/l2Sigma/l2Sigma +
              0.5*c*(l1-l1Mean)*(l2-l2Mean)/l1Sigma/l2Sigma ;	      
  return (R2<2.5*2.5) ;
  
}
//_____________________________________________________________________________
Bool_t AliPHOSHijingEfficiency::TestLambda2(Double_t pt,Double_t l1,Double_t l2){
  
//For R=4.5
  Double_t   l1Mean  = 1.150200 + 0.097886/(1.+1.486645*pt+0.000038*pt*pt) ;
  Double_t   l2Mean = 1.574706 + 0.997966*exp(-0.895075*pt)-0.010666*pt ;
  Double_t   l1Sigma = 0.100255 + 0.337177*exp(-0.517684*pt)+0.001170*pt ;
  Double_t   l2Sigma = 0.232580 + 0.573401*exp(-0.735903*pt)-0.002325*pt ;
  Double_t   c = -0.110983 -0.017353/(1.-1.836995*pt+0.934517*pt*pt) ;

/*
  Double_t l2Mean  = 1.53126+9.50835e+06/(1.+1.08728e+07*pt+1.73420e+06*pt*pt) ;
  Double_t l1Mean  = 1.12365+0.123770*TMath::Exp(-pt*0.246551)+5.30000e-03*pt ;
  Double_t l2Sigma = 6.48260e-02+7.60261e+10/(1.+1.53012e+11*pt+5.01265e+05*pt*pt)+9.00000e-03*pt;
  Double_t l1Sigma = 4.44719e-04+6.99839e-01/(1.+1.22497e+00*pt+6.78604e-07*pt*pt)+9.00000e-03*pt;
  Double_t c=-0.35-0.550*TMath::Exp(-0.390730*pt) ;
*/
  Double_t R2=0.5*(l1-l1Mean)*(l1-l1Mean)/l1Sigma/l1Sigma + 
              0.5*(l2-l2Mean)*(l2-l2Mean)/l2Sigma/l2Sigma +
              0.5*c*(l1-l1Mean)*(l2-l2Mean)/l1Sigma/l2Sigma ;
  return (R2<2.5*2.5) ;
  
}
//____________________________________________________________________________
Double_t AliPHOSHijingEfficiency::TestCPV(Double_t dx, Double_t dz, Double_t pt, Int_t charge){
  //Parameterization of LHC10h period
  //_true if neutral_
  
  Double_t meanX=0;
  Double_t meanZ=0.;
  Double_t sx=TMath::Min(5.4,2.59719e+02*TMath::Exp(-pt/1.02053e-01)+
              6.58365e-01*5.91917e-01*5.91917e-01/((pt-9.61306e-01)*(pt-9.61306e-01)+5.91917e-01*5.91917e-01)+1.59219);
  Double_t sz=TMath::Min(2.75,4.90341e+02*1.91456e-02*1.91456e-02/(pt*pt+1.91456e-02*1.91456e-02)+1.60) ;
  AliESDEvent *event = dynamic_cast<AliESDEvent*>(InputEvent());
  Double_t mf = event->GetMagneticField(); //Positive for ++ and negative for --

  if(mf<0.){ //field --
    meanZ = -0.468318 ;
    if(charge>0)
      meanX=TMath::Min(7.3, 3.89994*1.20679*1.20679/(pt*pt+1.20679*1.20679)+0.249029+2.49088e+07*TMath::Exp(-pt*3.33650e+01)) ;
    else
      meanX=-TMath::Min(7.7,3.86040*0.912499*0.912499/(pt*pt+0.912499*0.912499)+1.23114+4.48277e+05*TMath::Exp(-pt*2.57070e+01)) ;
  }
  else{ //Field ++
    meanZ= -0.468318;
    if(charge>0)
      meanX=-TMath::Min(8.0,3.86040*1.31357*1.31357/(pt*pt+1.31357*1.31357)+0.880579+7.56199e+06*TMath::Exp(-pt*3.08451e+01)) ;
    else
      meanX= TMath::Min(6.85, 3.89994*1.16240*1.16240/(pt*pt+1.16240*1.16240)-0.120787+2.20275e+05*TMath::Exp(-pt*2.40913e+01)) ;     
  }

  Double_t rz=(dz-meanZ)/sz ;
  Double_t rx=(dx-meanX)/sx ;
  return TMath::Sqrt(rx*rx+rz*rz) ;
}
//____________________________________________________________________________
Bool_t AliPHOSHijingEfficiency::TestTOF(Double_t t, Double_t e){

  Double_t sigma = TMath::Sqrt(2.23183e-09*2.23183e-09 
                             +2.24611e-09*2.24611e-09/e
                             +5.65054e-09*5.65054e-09/e/e) ;
  sigma=TMath::Min(20.e-9,sigma) ; //for the soft (<400 MeV) part
  Double_t mean=1.1e-9 ;
  if(TMath::Abs(t-mean)<2.*sigma)
    return kTRUE ;
  else
    if(TMath::Abs(t-mean+100.e-9)<2.*sigma)
      return kTRUE ;
    
  return kFALSE ;  
}
//____________________________________________________________________________
Double_t  AliPHOSHijingEfficiency::CoreEnergy(AliESDCaloCluster * clu){  
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
  return rnlin(coreE) ;
}
//____________________________________________________________________________
Bool_t  AliPHOSHijingEfficiency::AreNeibors(Int_t id1,Int_t id2){

  Int_t relid1[4] ; 
  fPHOSGeo->AbsToRelNumbering(id1, relid1) ; 

  Int_t relid2[4] ; 
  fPHOSGeo->AbsToRelNumbering(id2, relid2) ; 
 
  if ( (relid1[0] == relid2[0]) && (relid1[1]==relid2[1]) ) { // inside the same PHOS module 
    Int_t rowdiff = TMath::Abs( relid1[2] - relid2[2] ) ;  
    Int_t coldiff = TMath::Abs( relid1[3] - relid2[3] ) ;  
    
    if (( coldiff <= 1 )  && ( rowdiff <= 1 )){   //At least common vertex
      //    if (( relid1[2]==relid2[2] && coldiff <= 1 )  || ( relid1[3]==relid2[3] &&  rowdiff <= 1 )){ //common side
      return 1 ; 
    }
    else {
      if((relid2[2] > relid1[2]) && (relid2[3] > relid1[3]+1)) 
        return 0; //  Difference in row numbers is too large to look further 
    }
    return 0 ;

  } 

  return 0 ; 
}
//____________________________________________________________________________
void  AliPHOSHijingEfficiency::Reclusterize(AliESDCaloCluster * clu){
  //Re-clusterize to make continues cluster
  
  const Int_t oldMulDigit=clu->GetNCells() ;
  Double32_t * elist = clu->GetCellsAmplitudeFraction() ;  
  UShort_t * dlist = clu->GetCellsAbsId();

  Int_t index[oldMulDigit] ;
  Bool_t used[oldMulDigit] ;
  for(Int_t i=0; i<oldMulDigit; i++) used[i]=0 ;
  Int_t inClu=0 ;
  Double_t eMax=0. ;
  //find maximum
  for(Int_t iDigit=0; iDigit<oldMulDigit; iDigit++) {
    if(eMax<elist[iDigit]){
      eMax=elist[iDigit];
      index[0]=iDigit ;
      inClu=1 ;
    }
  }
  if(inClu==0){ //empty cluster
    return ;
  }
  used[index[0]]=kTRUE ; //mark as used
  for(Int_t i=0; i<inClu; i++){
    for(Int_t iDigit=0 ;iDigit<oldMulDigit; iDigit++){
       if(used[iDigit]) //already used
         continue ;
       if(AreNeibors(dlist[index[i]],dlist[iDigit])){
	 index[inClu]= iDigit ;
	 inClu++ ;
	 used[iDigit]=kTRUE ;
       }
    }
  }
  
  if(inClu==oldMulDigit) //no need to modify
    return ; 

  clu->SetNCells(inClu);
  //copy
  UShort_t tmpD[oldMulDigit] ;
  Double_t tmpE[oldMulDigit] ;
  for(Int_t i=0; i<oldMulDigit; i++){
    tmpD[i]=dlist[i] ;    
    tmpE[i]=elist[i] ;
  }
  //change order of digits in list so that
  //first inClu cells were true ones
  for(Int_t i=0; i<inClu; i++){
    dlist[i]=tmpD[index[i]] ;
    elist[i]=tmpE[index[i]] ;
  }
  
  
}

//____________________________________________________________________________
void  AliPHOSHijingEfficiency::EvalLambdas(AliESDCaloCluster * clu, Int_t iR,Double_t &m02, Double_t &m20){ 
  //calculate dispecrsion of the cluster in the circle with radius distanceCut around the maximum
    
  Double_t rCut=0. ;  
  if(iR==0)    
    rCut=3.5 ;
  else
    rCut=4.5 ;
  
  
  Double32_t * elist = clu->GetCellsAmplitudeFraction() ;  
// Calculates the center of gravity in the local PHOS-module coordinates
  Float_t wtot = 0;
  Double_t xc[100]={0} ;
  Double_t zc[100]={0} ;
  Double_t x = 0 ;
  Double_t z = 0 ;
  Int_t mulDigit=TMath::Min(100,clu->GetNCells()) ;
  const Double_t logWeight=4.5 ;
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
     
  wtot = 0. ;
  Double_t dxx  = 0.;
  Double_t dzz  = 0.;
  Double_t dxz  = 0.;
  Double_t xCut = 0. ;
  Double_t zCut = 0. ;
  for(Int_t iDigit=0; iDigit<mulDigit; iDigit++) {
    if (clu->E()>0 && elist[iDigit]>0.) {
        Double_t w = TMath::Max( 0., logWeight + TMath::Log( elist[iDigit] / clu->E() ) ) ;
        Double_t xi= xc[iDigit] ;
        Double_t zi= zc[iDigit] ;
	if((xi-x)*(xi-x)+(zi-z)*(zi-z) < rCut*rCut){
          xCut += w * xi ;
          zCut += w * zi ; 
          dxx  += w * xi * xi ;
          dzz  += w * zi * zi ;
          dxz  += w * xi * zi ; 
          wtot += w ;
	}
    }
    
  }
  if (wtot>0) {
    xCut/= wtot ;
    zCut/= wtot ;
    dxx /= wtot ;
    dzz /= wtot ;
    dxz /= wtot ;
    dxx -= xCut * xCut ;
    dzz -= zCut * zCut ;
    dxz -= xCut * zCut ;

    m02 =  0.5 * (dxx + dzz) + TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz )  ;
    m20 =  0.5 * (dxx + dzz) - TMath::Sqrt( 0.25 * (dxx - dzz) * (dxx - dzz) + dxz * dxz )  ;
  }
  else {
    m20=m02=0.;
  }

}
//___________________________________________________________________________
void AliPHOSHijingEfficiency::ProcessMC(){
  //fill histograms for efficiensy etc. calculation
  const Double_t rcut = 1. ; //cut for primary particles
  //---------First pi0/eta-----------------------------
  char partName[10] ;
  char hkey[55] ;

  if(!fStack) return ;
  for(Int_t i=0;i<fStack->GetNtrack();i++){
     TParticle* particle =  fStack->Particle(i);
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
    Double_t r=particle->R() ;
    Double_t pt = particle->Pt() ;
    //Distribution over vertex
    FillHistogram(Form("hMC_%s_vertex",partName),pt,r) ;
    
    if(r >rcut)
      continue ;

    //Total number of pi0 with creation radius <1 cm
    Double_t w = PrimaryParticleWeight(particle) ;  
    snprintf(hkey,55,"hMC_all_%s_cen%d",partName,fCenBin) ;
    FillHistogram(hkey,pt,w) ;
    if(TMath::Abs(particle->Y())<0.12){
      snprintf(hkey,55,"hMC_unitEta_%s_cen%d",partName,fCenBin) ;
      FillHistogram(hkey,pt,w) ;
    }

    snprintf(hkey,55,"hMC_rap_%s_cen%d",partName,fCenBin) ;
    FillHistogram(hkey,particle->Y(),w) ;
    
    Double_t phi=particle->Phi() ;
    while(phi<0.)phi+=TMath::TwoPi() ;
    while(phi>TMath::TwoPi())phi-=TMath::TwoPi() ;
    snprintf(hkey,55,"hMC_phi_%s_cen%d",partName,fCenBin) ;
    FillHistogram(hkey,phi,w) ;
  }
}
//___________________________________________________________________________
Int_t AliPHOSHijingEfficiency::FindPrimary(AliESDCaloCluster*clu,  Bool_t&sure){
  //Finds primary and estimates if it unique one?
  Int_t n=clu->GetNLabels() ;
  Double_t*  Ekin=  new  Double_t[n] ;
  for(Int_t i=0;  i<n;  i++){
    TParticle*  p=  fStack->Particle(clu->GetLabelAt(i)) ;
    Ekin[i]=p->P() ;  // estimate of kinetic energy
    if(p->GetPdgCode()==-2212  ||  p->GetPdgCode()==-2112){
      Ekin[i]+=1.8  ;  //due to annihilation
    }
  }
  Int_t iMax=0;
  Double_t eMax=0.,eSubMax=0. ;
  for(Int_t i=0;  i<n;  i++){
    if(Ekin[i]>eMax){
      eSubMax=eMax;
      eMax=Ekin[i];
      iMax=i;
    }
  }
  if(eSubMax>0.8*eMax)//not obvious primary
    sure=kFALSE;
  else
    sure=kTRUE;
  delete  Ekin;
  return  clu->GetLabelAt(iMax) ;
}
//___________________________________________________________________________
Double_t AliPHOSHijingEfficiency::PrimaryWeight(Int_t primary){
   //Check who is the primary and introduce weight to correct primary spectrum
  
  if(primary<0 || primary>=fStack->GetNtrack())
    return 1 ;
  //trace primaries up to IP
  TParticle* particle =  fStack->Particle(primary);
  Double_t r=particle->R() ;
  Int_t mother = particle->GetFirstMother() ;
  while(mother>-1){
    if(r<1. && particle->GetPdgCode()==111)
      break ;
    particle =  fStack->Particle(mother);
    mother = particle->GetFirstMother() ;
    r=particle->R() ;
  }

  return TMath::Max(0.,PrimaryParticleWeight(particle)) ;
}
//________________________________________________________________________
Double_t AliPHOSHijingEfficiency::PrimaryParticleWeight(TParticle * particle){
  
  Int_t pdg = particle->GetPdgCode() ;
  Int_t type=0 ;
  if(pdg == 111 || TMath::Abs(pdg)==211){
    type =1 ;
  }
  else{
    if(TMath::Abs(pdg)<1000){ //Kaon-like
      type =2 ;    
    }
    else
      type = 3;  //baryons
  }
    
  Double_t x = particle->Pt() ;
  if(type==1){
   if(fCenBin==0) //0-5
     return (1.662990+1.140890*x-0.192088*x*x)/(1.-0.806630*x+0.304771*x*x)+0.141690*x ;
   if(fCenBin==1) //5-10
     return (1.474351+0.791492*x-0.066369*x*x)/(1.-0.839338*x+0.317312*x*x)+0.093289*x ;
   if(fCenBin==2) //10-20
     return (1.174728+0.959681*x-0.137695*x*x)/(1.-0.788873*x+0.299538*x*x)+0.128759*x ; 
   if(fCenBin==3) //20-40
     return (0.927335+0.475349*x+0.004364*x*x)/(1.-0.817966*x+0.309787*x*x)+0.086899*x ; 
   if(fCenBin==4) //40-60
     return (0.676878+0.190680*x+0.077031*x*x)/(1.-0.790623*x+0.305183*x*x)+0.064510*x ; 
   if(fCenBin==5) //60-80
     return (0.684726-0.606262*x+0.409963*x*x)/(1.-1.080061*x+0.456933*x*x)+0.005151*x ; 
  }
  if(type==2){
   if(fCenBin==0) //0-5
     return (-0.417131+2.253936*x-0.337731*x*x)/(1.-0.909892*x+0.316820*x*x)+0.157312*x ;
   if(fCenBin==1) //5-10
     return (-0.352275+1.844466*x-0.248598*x*x)/(1.-0.897048*x+0.316462*x*x)+0.132461*x ; 
   if(fCenBin==2) //10-20
     return (-0.475481+1.975108*x-0.336013*x*x)/(1.-0.801028*x+0.276705*x*x)+0.188164*x ; 
   if(fCenBin==3) //20-40
     return (-0.198954+1.068789*x-0.103540*x*x)/(1.-0.848354*x+0.299209*x*x)+0.112939*x ; 
   if(fCenBin==4) //40-60
     return (-0.111052+0.664041*x-0.019717*x*x)/(1.-0.804916*x+0.300779*x*x)+0.095784*x ;
   if(fCenBin==5) //0-5
     return (0.202788-0.439832*x+0.564585*x*x)/(1.-1.254029*x+0.679444*x*x)+0.016235*x ;
  }
  if(type==3){
   if(fCenBin==0) //0-5
     return (-1.312732+2.743568*x-0.375775*x*x)/(1.-0.717533*x+0.164694*x*x)+0.164445*x ;
   if(fCenBin==1) //5-10
     return (-1.229425+2.585889*x-0.330164*x*x)/(1.-0.715892*x+0.167386*x*x)+0.133085*x ; 
   if(fCenBin==2) //10-20
     return (-1.135677+2.397489*x-0.320355*x*x)/(1.-0.709312*x+0.164350*x*x)+0.146095*x ; 
   if(fCenBin==3) //20-40
     return (-0.889993+1.928263*x-0.220785*x*x)/(1.-0.715991*x+0.174729*x*x)+0.095098*x ; 
   if(fCenBin==4) //40-60
     return (-0.539237+1.329118*x-0.115439*x*x)/(1.-0.722906*x+0.186832*x*x)+0.059267*x ; 
   if(fCenBin==5) //60-80
     return (-0.518126+1.327628*x-0.130881*x*x)/(1.-0.665649*x+0.184300*x*x)+0.081701*x ;   
  }
  return 1. ;  
}
//________________________________________________________________________
Int_t AliPHOSHijingEfficiency::FindCommonParent(Int_t iPart, Int_t jPart){
  //check if there is a common parent for particles i and j
  // -1: no common parent or wrong iPart/jPart
  
  if(iPart==-1 || iPart>=fStack->GetNtrack() || 
     jPart==-1 || jPart>=fStack->GetNtrack()) return -1;
  
  Int_t iprim1=iPart;
  while(iprim1>-1){  
     Int_t iprim2=jPart;
     while(iprim2>-1){
       if(iprim1==iprim2)
	 return iprim1 ;
       iprim2=((TParticle *)fStack->Particle(iprim2))->GetFirstMother();
     }
     iprim1=((TParticle *)fStack->Particle(iprim1))->GetFirstMother();
  }
  return -1;
}
//________________________________________________________________________
Bool_t AliPHOSHijingEfficiency::HaveParent(Int_t iPart, Int_t pdgParent){
  //check if there is a common parent for particles i and j
  // -1: no common parent or wrong iPart/jPart
  
  if(iPart==-1 || iPart>=fStack->GetNtrack()) return -1;
  
  Int_t iprim1=iPart;
  while(iprim1>-1){  
    TParticle * tmp = fStack->Particle(iprim1) ;
    if(tmp->GetPdgCode()==pdgParent)
      return kTRUE ;
    iprim1=tmp->GetFirstMother();
  }
  return kFALSE;
}
//________________________________________________________________________
Bool_t AliPHOSHijingEfficiency::InPi0mass(Double_t m, Double_t /*pt*/){

 return TMath::Abs(m-0.135)<0.007*2.5 ;
}
