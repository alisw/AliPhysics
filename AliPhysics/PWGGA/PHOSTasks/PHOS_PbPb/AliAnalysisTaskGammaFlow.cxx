#include "TChain.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH3F.h"
#include "TProfile.h"
#include "TParticle.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TRandom.h"
#include "TROOT.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "THashList.h"

#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskGammaFlow.h"
#include "AliCaloPhoton.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSEsdCluster.h"
#include "AliESDEvent.h"
#include "AliAODEvent.h"
#include "AliESDCaloCells.h"
#include "AliAODVertex.h"
#include "AliESDtrackCuts.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliCDBManager.h"
#include "AliCentrality.h" 
#include "AliEventplane.h"
#include "AliOADBContainer.h"
#include "AliEPFlattener.h"
#include "AliFlowTrackCuts.h"
#include "AliFlowEvent.h"
#include "AliFlowVector.h"

// Analysis task to fill histograms with PHOS ESD clusters and cells
// Authors: Dmitri Peressounko
// Date   : 28.05.2011

ClassImp(AliAnalysisTaskGammaFlow)

//________________________________________________________________________
AliAnalysisTaskGammaFlow::AliAnalysisTaskGammaFlow(const char *name) 
: AliAnalysisTaskSE(name),
  fOutputContainer(0x0),
  fPHOSEvent(0x0),
  fQV0A(0.),
  fQV0C(0.),
  fQTPC(0.),
  fRP(0.),
  fRPV0A(0.),
  fRPV0C(0.),
  fRPQ(0.),
  fRPQV0A(0.),
  fRPQV0C(0.),
  fV0Ares(1.),
  fV0Cres(1.),
  fTPCres(1.),
  fV0AQres(1.),
  fV0CQres(1.),
  fTPCQres(1.),
  fHaveTPCRP(0),
  fPhiDist(0x0),
  fV0AFlat(0x0),
  fV0CFlat(0x0),
  fTPCFlat(0x0),
  fV0AQFlat(0x0),
  fV0CQFlat(0x0),
  fTPCQFlat(0x0),
  fHarmonics(3),
  fDistCut(0),
  fRunNumber(0),
  fCentrality(0.),
  fCenBin(0),
  fPHOSGeo(0x0),
  fEventCounter(0),

  fTPCfinalC2(0x0),   //HIstos with flattening parameters
  fTPCfinalS2(0x0),
  fTPCfinalC4(0x0),
  fTPCfinalS4(0x0),

  fV0AfinalC2(0x0),
  fV0AfinalS2(0x0),
  fV0AfinalC4(0x0),
  fV0AfinalS4(0x0),
    
  fV0CfinalC2(0x0),
  fV0CfinalS2(0x0),
  fV0CfinalC4(0x0),
  fV0CfinalS4(0x0), 
  fTPCfinalQC2(0x0),   //HIstos with flattening parameters
  fTPCfinalQS2(0x0),
  fTPCfinalQC4(0x0),
  fTPCfinalQS4(0x0),

  fV0AfinalQC2(0x0),
  fV0AfinalQS2(0x0),
  fV0AfinalQC4(0x0),
  fV0AfinalQS4(0x0),
    
  fV0CfinalQC2(0x0),
  fV0CfinalQS2(0x0),
  fV0CfinalQC4(0x0),
  fV0CfinalQS4(0x0),
  fCutsV0(0x0),
  fCutsTPC(0x0),
  fFlowEvent(0x0)

{
  // Constructor
  for(Int_t i=0;i<1;i++){
    for(Int_t j=0;j<10;j++)
      for(Int_t k=0;k<11;k++)
	fPHOSEvents[i][j][k]=0 ;
  }
  
  // Output slots #0 write into a TH1 container
  DefineOutput(1,TList::Class());

  // Initialize the PHOS geometry
  fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP") ;

  //We have to apply re-calibration for pass1 LCH10h
  // Initialize decalibration factors in the form of the OCDB object
  fCutsV0 = new AliFlowTrackCuts(Form("V0%d",fHarmonics));
  fCutsV0 = fCutsV0->GetStandardVZEROOnlyTrackCuts(); // select vzero tracks
  fCutsV0->SetVZEROgainEqualizationPerRing(kFALSE);
  fCutsV0->SetApplyRecentering(kTRUE);
//  fCutsV0A->SetEtaRange(2.,10.) ; 
  fFlowEvent = new AliFlowEvent(10000);

  fCutsTPC= new AliFlowTrackCuts(Form("TPC%d",fHarmonics));
  fCutsTPC=fCutsTPC->GetStandardTPCStandaloneTrackCuts() ;
  fCutsTPC->SetEtaMin(-0.9);
  fCutsTPC->SetEtaMax(0.9);

 }

//________________________________________________________________________
void AliAnalysisTaskGammaFlow::UserCreateOutputObjects()
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
  
  //========QA histograms=======

  //Event selection
  fOutputContainer->Add(new TH2F("hSelEvents","Event selection", 10,0.,10.,nRuns,0.,float(nRuns))) ;
  fOutputContainer->Add(new TH1F("hTotSelEvents","Event selection", 10,0.,10.)) ;
  
  //vertex distribution
  fOutputContainer->Add(new TH2F("hZvertex","Z vertex position", 50,-25.,25.,nRuns,0.,float(nRuns))) ;
  
  //Centrality
  fOutputContainer->Add(new TH2F("hCentrality","Event centrality", 100,0.,100.,nRuns,0.,float(nRuns))) ;
  fOutputContainer->Add(new TH1F("hCentralityCorr","Event centrality", 100,0.,100.)) ;
  fOutputContainer->Add(new TH1F("hBadCentrality","Bad PHOS event centrality", 100,0.,100.)) ;
  fOutputContainer->Add(new TH2F("hCenPHOS","Centrality vs PHOSclusters", 100,0.,100.,200,0.,200.)) ;
  fOutputContainer->Add(new TH2F("hCenPHOSCells","Centrality vs PHOS cells", 100,0.,100.,100,0.,1000.)) ;
  fOutputContainer->Add(new TH3F("hCenPHOSCellsM12","Centrality vs PHOS cells", 100,0.,100.,100,0.,1000.,100,0.,1000.)) ;
  fOutputContainer->Add(new TH3F("hCenPHOSCellsM23","Centrality vs PHOS cells", 100,0.,100.,100,0.,1000.,100,0.,1000.)) ;
  fOutputContainer->Add(new TH3F("hCenPHOSCellsM13","Centrality vs PHOS cells", 100,0.,100.,100,0.,1000.,100,0.,1000.)) ;
  fOutputContainer->Add(new TH2F("hCenTrack","Centrality vs tracks", 100,0.,100.,100,0.,15000.)) ;  
  fOutputContainer->Add(new TH2F("hCluEvsClu","ClusterMult vs E",200,0.,20.,100,0.,100.)) ;
  fOutputContainer->Add(new TH2F("hCluEvsCluM","ClusterMult vs E",200,0.,20.,100,0.,20.)) ;
  fOutputContainer->Add(new TH2F("hCenTOF","Centrality vs PHOS TOF", 100,0.,100.,600,-6.e-6,6.e-6)) ;
  
  //Reaction plane
//  fOutputContainer->Add(new TH3F("hPHOSphi","cos" ,10,0.,100.,20,0.,10.,100,-TMath::Pi(),TMath::Pi()));

   fOutputContainer->Add(new TProfile("ResV0A","Resolution V0A", 100,0.,100.)) ;
   fOutputContainer->Add(new TProfile("ResV0C","Resolution V0C", 100,0.,100.)) ;
   fOutputContainer->Add(new TProfile("ResTPC","Resolution TPC", 100,0.,100.)) ;
   fOutputContainer->Add(new TProfile("InvResV0A","Inverse Resolution V0A", 100,0.,100.)) ;
   fOutputContainer->Add(new TProfile("InvResV0C","Inverse Resolution V0C", 100,0.,100.)) ;
   fOutputContainer->Add(new TProfile("InvResTPC","Inverse Resolution TPC", 100,0.,100.)) ;  
   
   fOutputContainer->Add(new TProfile("cos2TPCAC","RP correlation between TPC subs",100,0.,100.)) ;
   fOutputContainer->Add(new TProfile("cos2V0AC","RP correlation between VO A and C sides",100,0.,100.)) ;
   fOutputContainer->Add(new TProfile("cos2V0ATPC","RP correlation between TPC and V0A",100,0.,100.)) ;
   fOutputContainer->Add(new TProfile("cos2V0CTPC","RP correlation between TPC and V0C",100,0.,100.)) ;
   fOutputContainer->Add(new TProfile("qcos2TPCAC","RP correlation between TPC subs",100,0.,100.)) ;
   fOutputContainer->Add(new TProfile("qcos2V0AC","RP correlation between VO A and C sides",100,0.,100.)) ;
   fOutputContainer->Add(new TProfile("qcos2V0ATPC","RP correlation between TPC and V0A",100,0.,100.)) ;
   fOutputContainer->Add(new TProfile("qcos2V0CTPC","RP correlation between TPC and V0C",100,0.,100.)) ;
   
   fOutputContainer->Add(new TH2F("QV0A","Q_{V0A}",100,0.,100.,200,0.,200.)) ;
   fOutputContainer->Add(new TH2F("QV0C","Q_{V0C}",100,0.,100.,200,0.,300.)) ;
   fOutputContainer->Add(new TH2F("QTPC","Q_{TPC}",100,0.,100.,200,0.,100.)) ;
   

   fOutputContainer->Add(new TProfile("resV0A","Estimated resolution of V0A",100,0.,100.)) ;
   fOutputContainer->Add(new TProfile("resV0C","Estimated resolution of V0C",100,0.,100.)) ;
   fOutputContainer->Add(new TProfile("resTPC","Estimated resolution of TPC",100,0.,100.)) ;
   fOutputContainer->Add(new TProfile("qresV0A","Estimated resolution of V0A",100,0.,100.)) ;
   fOutputContainer->Add(new TProfile("qresV0C","Estimated resolution of V0C",100,0.,100.)) ;
   fOutputContainer->Add(new TProfile("qresTPC","Estimated resolution of TPC",100,0.,100.)) ;
   
   fOutputContainer->Add(new TProfile("VaVbcen","Estimated resolution of TPC",200,0.,100.)) ;
   
   fOutputContainer->Add(new TH2F("phiRP","RP distribution with TPC", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.)) ;
   fOutputContainer->Add(new TH2F("phiRPflat","RP distribution with TPC flat", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.)) ;
   fOutputContainer->Add(new TH2F("phiRPV0A","RP distribution with V0A", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.)) ;
   fOutputContainer->Add(new TH2F("phiRPV0C","RP distribution with V0C", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.)) ;
   fOutputContainer->Add(new TH2F("phiRPV0Aflat","RP distribution with V0 flat", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.)) ;
   fOutputContainer->Add(new TH2F("phiRPV0Cflat","RP distribution with V0 flat", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.)) ;
   fOutputContainer->Add(new TH2F("phiRPV0AFlow","RP V0A, Flow Package", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.)) ;
   fOutputContainer->Add(new TH2F("phiRPV0CFlow","RP V0C, Flow Package", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.)) ;
                                   
   fOutputContainer->Add(new TProfile2D("phiRPQ","RP distribution with TPC flat", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.)) ;
   fOutputContainer->Add(new TProfile2D("phiRPV0AQ","RP distribution with V0 flat", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.)) ;
   fOutputContainer->Add(new TProfile2D("phiRPV0CQ","RP distribution with V0 flat", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.)) ;
   fOutputContainer->Add(new TProfile2D("phiRPQflat","RP distribution with TPC flat", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.)) ;
   fOutputContainer->Add(new TProfile2D("phiRPV0AQflat","RP distribution with V0 flat", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.)) ;
   fOutputContainer->Add(new TProfile2D("phiRPV0CQflat","RP distribution with V0 flat", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.)) ;
   fOutputContainer->Add(new TH3F("phiRPV0AC","RP distribution with V0A + TPC", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.)) ;
   fOutputContainer->Add(new TH3F("phiRPV0ATPC","RP distribution with V0A + TPC", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.)) ;
   fOutputContainer->Add(new TH3F("phiRPV0CTPC","RP distribution with V0C + TPC", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.)) ;

   fOutputContainer->Add(new TH3F("phiV0ACorrel","V0A my vs Flow", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.)) ;
   fOutputContainer->Add(new TH3F("phiV0CCorrel","V0C my vs Flow", 100,0.,TMath::TwoPi()/fHarmonics,100,0.,TMath::TwoPi()/fHarmonics,100,0.,100.)) ;    
   
   fOutputContainer->Add(new TProfile2D("cos2TPCRP","cos(2Psi_{EP}^{TPC}", 200,-1.,1.,100,0.,100.)) ;
   fOutputContainer->Add(new TProfile2D("cos2V0ARP","cos(2Psi_{EP}^{V0A}", 200,-1.,1.,100,0.,100.)) ;
   fOutputContainer->Add(new TProfile2D("cos2V0CRP","cos(2Psi_{EP}^{V0C}", 200,-1.,1.,100,0.,100.)) ;
   fOutputContainer->Add(new TProfile2D("sin2TPCRP","sin(2Psi_{EP}^{TPC}", 200,-1.,1.,100,0.,100.)) ;
   fOutputContainer->Add(new TProfile2D("sin2V0ARP","sin(2Psi_{EP}^{V0A}", 200,-1.,1.,100,0.,100.)) ;
   fOutputContainer->Add(new TProfile2D("sin2V0CRP","sin(2Psi_{EP}^{V0C}", 200,-1.,1.,100,0.,100.)) ;  
  
  //PHOS QA
  fOutputContainer->Add(new TH1F("hBadMod","PHOS module without cells",6,2.,8.));
  fOutputContainer->Add(new TH1I("hCellMultEvent"  ,"PHOS cell multiplicity per event"    ,2000,0,2000));
  fOutputContainer->Add(new TH1I("hCellMultEventM1","PHOS cell multiplicity per event, M1",2000,0,2000));
  fOutputContainer->Add(new TH1I("hCellMultEventM2","PHOS cell multiplicity per event, M2",2000,0,2000));
  fOutputContainer->Add(new TH1I("hCellMultEventM3","PHOS cell multiplicity per event, M3",2000,0,2000));

  fOutputContainer->Add(new TH1F("hCellEnergy"  ,"Cell energy"            ,3000,0.,30.));
  fOutputContainer->Add(new TH1F("hCellEnergyM1","Cell energy in module 1",3000,0.,30.));
  fOutputContainer->Add(new TH1F("hCellEnergyM2","Cell energy in module 2",3000,0.,30.));
  fOutputContainer->Add(new TH1F("hCellEnergyM3","Cell energy in module 3",3000,0.,30.));

  fOutputContainer->Add(new TH2F("hCellNXZM1","Cell (X,Z), M1" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCellNXZM2","Cell (X,Z), M2" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCellNXZM3","Cell (X,Z), M3" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCellEXZM1","Cell E(X,Z), M1",64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCellEXZM2","Cell E(X,Z), M2",64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCellEXZM3","Cell E(X,Z), M3",64,0.5,64.5, 56,0.5,56.5));
 			
//  fOutputContainer->Add(new TH3F("hCPVr","CPV radius",100,0.,20.,100,0.,20.,10,0.,100.));
//  fOutputContainer->Add(new TH3F("hLambdaCent","Lambdas for all clusters",50,0.,10.,50,0.,10.,5,0.,100.));
  
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

  fOutputContainer->Add(new TH2F("hTofM1","TOF in M1" ,100,0.,20.,400,-4.e-6,4.e-6));
  fOutputContainer->Add(new TH2F("hTofM2","TOF in M2" ,100,0.,20.,400,-4.e-6,4.e-6));
  fOutputContainer->Add(new TH2F("hTofM3","TOF in M3" ,100,0.,20.,400,-4.e-6,4.e-6));
  
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
    
  fOutputContainer->Add(new TH2F("hMiPi0M11","Pairs in modules",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hMiPi0M12","Pairs in modules",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hMiPi0M13","Pairs in modules",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hMiPi0M22","Pairs in modules",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hMiPi0M23","Pairs in modules",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hMiPi0M33","Pairs in modules",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  
  
  char cPID[8][15] ;
  snprintf(cPID[0],15,"All") ;
  snprintf(cPID[1],15,"Allcore") ;
  snprintf(cPID[2],15,"Disp2");
  snprintf(cPID[3],15,"Disp2core");
  snprintf(cPID[4],15,"CPV") ;
  snprintf(cPID[5],15,"CPVcore") ;
  snprintf(cPID[6],15,"Both2"); 
  snprintf(cPID[7],15,"Both2core"); 

    
  for(Int_t cent=0; cent<7; cent++){
    fOutputContainer->Add(new TH2F(Form("hPHOSphiAll_cen%d",cent),"pair pT vs Phi",nPtPhot,0.,ptPhotMax,100,-TMath::Pi(),TMath::Pi()));
    fOutputContainer->Add(new TH2F(Form("hPHOScos2All_cen%d",cent),"pair pT vs cos(2Phi)",nPtPhot,0.,ptPhotMax,100,-1,1));
    fOutputContainer->Add(new TH2F(Form("hPHOSsin2All_cen%d",cent),"pair pT vs sin(2Phi)",nPtPhot,0.,ptPhotMax,100,-1,1));

    fOutputContainer->Add(new TH2F(Form("hPHOSphiDisp_cen%d",cent),"pair pT vs Phi",nPtPhot,0.,ptPhotMax,100,-TMath::Pi(),TMath::Pi()));
    fOutputContainer->Add(new TH2F(Form("hPHOSphiDispcore_cen%d",cent),"pair pT vs Phi",nPtPhot,0.,ptPhotMax,100,-TMath::Pi(),TMath::Pi()));
    fOutputContainer->Add(new TH2F(Form("hPHOSphiCPV_cen%d",cent),"pair pT vs Phi",nPtPhot,0.,ptPhotMax,100,-TMath::Pi(),TMath::Pi()));

    
    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad2_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad4_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad6_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad2core_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad4core_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad6core_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad2Disp_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad4Disp_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad6Disp_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad2Dispcore_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad4Dispcore_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad6Dispcore_cen%d",cent),"All clusters",nPtPhot,0.,ptPhotMax));

    for(Int_t iPID=0; iPID<8; iPID++){
      fOutputContainer->Add(new TH1F(Form("hPhot%s_cen%d",cPID[iPID],cent),"All clusters",nPtPhot,0.,ptPhotMax));
      fOutputContainer->Add(new TH1F(Form("hPhotRa%s_cen%d",cPID[iPID],cent),"All clusters",nPtPhot,0.,ptPhotMax));
      ((TH1F*)fOutputContainer->Last())->Sumw2() ;
      fOutputContainer->Add(new TH1F(Form("hPhotRc%s_cen%d",cPID[iPID],cent),"All clusters",nPtPhot,0.,ptPhotMax));
      ((TH1F*)fOutputContainer->Last())->Sumw2() ;
      fOutputContainer->Add(new TH1F(Form("hPhotRt%s_cen%d",cPID[iPID],cent),"All clusters",nPtPhot,0.,ptPhotMax));
      ((TH1F*)fOutputContainer->Last())->Sumw2() ;
      fOutputContainer->Add(new TH1F(Form("hPhotRNa%s_cen%d",cPID[iPID],cent),"All clusters",nPtPhot,0.,ptPhotMax));
      ((TH1F*)fOutputContainer->Last())->Sumw2() ;
      fOutputContainer->Add(new TH1F(Form("hPhotRNc%s_cen%d",cPID[iPID],cent),"All clusters",nPtPhot,0.,ptPhotMax));
      ((TH1F*)fOutputContainer->Last())->Sumw2() ;
      fOutputContainer->Add(new TH1F(Form("hPhotRNt%s_cen%d",cPID[iPID],cent),"All clusters",nPtPhot,0.,ptPhotMax));
      ((TH1F*)fOutputContainer->Last())->Sumw2() ;
      for(Int_t m=1; m<4; m++){
        fOutputContainer->Add(new TH1F(Form("hPhot%s_m%d_cen%d",cPID[iPID],m,cent),"All clusters",nPtPhot,0.,ptPhotMax));
      }
    }        
    fOutputContainer->Add(new TH1F(Form("hPhotCPV2_cen%d",cent),"CPV clusters",nPtPhot,0.,ptPhotMax));  

  }

  
  char phiTitle[15] ;
  Int_t nPt      = 150;
  Double_t xPt[151] ;
  for(Int_t i=0;i<=150;i++){xPt[i]=0.1*i;}
  Int_t nPhi=10 ;
  Double_t xPhi[11] ;
  for(Int_t i=0;i<=10;i++)xPhi[i]=i*0.1*TMath::TwoPi()/fHarmonics ;
    
  for(Int_t iRP=0; iRP<3; iRP++){
    if(iRP==0)
      snprintf(phiTitle,15,"TPC") ;
    if(iRP==1)
      snprintf(phiTitle,15,"V0A") ;
    if(iRP==2)
      snprintf(phiTitle,15,"V0C") ;
    for(Int_t iPID=0; iPID<8; iPID++){

      fOutputContainer->Add(new TProfile2D(Form("hPhotPhiSP%s%s",phiTitle,cPID[iPID]),"(M,p_{T},d#phi)_{#gamma#gamma}" ,nPt,xPt,200,0.,100.)) ;
      ((TProfile2D*)fOutputContainer->Last())->Sumw2() ;
      
      for(Int_t cent=0; cent<7; cent++){
	
	
        fOutputContainer->Add(new TProfile(Form("hPhotPhiR2%s%s_cen%d",phiTitle,cPID[iPID],cent),"(M,p_{T},d#phi)_{#gamma#gamma}" ,nPt,xPt)) ;
        ((TProfile*)fOutputContainer->Last())->Sumw2() ;

        fOutputContainer->Add(new TProfile(Form("hPhotPhiR%s%s_cen%d",phiTitle,cPID[iPID],cent),"(M,p_{T},d#phi)_{#gamma#gamma}" ,nPt,xPt)) ;
        ((TProfile*)fOutputContainer->Last())->Sumw2() ;

        fOutputContainer->Add(new TProfile(Form("hPhotcos%s%s_cen%d",phiTitle,cPID[iPID],cent),"All clusters",nPtPhot,0.,ptPhotMax));
        ((TProfile*)fOutputContainer->Last())->Sumw2() ;
      
        fOutputContainer->Add(new TProfile(Form("hPhotcosN%s%s_cen%d",phiTitle,cPID[iPID],cent),"All clusters",nPtPhot,0.,ptPhotMax));
        ((TProfile*)fOutputContainer->Last())->Sumw2() ;
      }
    }
  }

  PostData(1, fOutputContainer);

}

//________________________________________________________________________
void AliAnalysisTaskGammaFlow::UserExec(Option_t *) 
{
  // Main loop, called for each event
  // Analyze ESD/AOD
      
  
  AliAODEvent *event = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!event) {
    Printf("ERROR: Could not retrieve event");
    PostData(1, fOutputContainer);
    return;
  }

  fRunNumber=ConvertRunNumber(event->GetRunNumber()) ;
  FillHistogram("hSelEvents",0.5,fRunNumber-0.5) ;
  FillHistogram("hTotSelEvents",0.5) ;
  
  
  
  
  // Get PHOS rotation matrices from ESD and set them to the PHOS geometry
  if(fEventCounter == 0) {
//   OpenInfoCalbration(event->GetRunNumber()) ;
   for(Int_t mod=0; mod<5; mod++) {
      if(!event->GetPHOSMatrix(mod)) continue;
      fPHOSGeo->SetMisalMatrix(event->GetPHOSMatrix(mod),mod) ;
      Printf("PHOS geo matrix %p for module # %d is set\n", event->GetPHOSMatrix(mod), mod);
    }
    Int_t run = event->GetRunNumber() ;
//    TFile * fflatFine = TFile::Open("EP_final.root") ;
//    TFile * fflatFineQ = TFile::Open("EPq_final.root") ;
    AliOADBContainer flatContainer("phosFlat");
    flatContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSflat.root","phosFlat");
    TObjArray *arr = (TObjArray*)flatContainer.GetObject(run,"phosFlat");
    if(!arr){
      AliError(Form("Can not read Flattening for run %d. \n From file $ALICE_PHYSICS/OADB/PHOS/PHOSflat.root",run)) ;    
      arr = (TObjArray*)flatContainer.GetObject(1,"phosFlat"); //default
    }
    AliOADBContainer flatContainerQ("phosQFlat");
    flatContainerQ.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSflat.root","phosQFlat");
    TObjArray *arrQ = (TObjArray*)flatContainerQ.GetObject(run,"phosQFlat");
    if(!arrQ){
      AliError(Form("Can not read QFlattening for run %d. \n From file $ALICE_PHYSICS/OADB/PHOS/PHOSflat.root",run)) ;    
    }
        
    AliInfo(Form("Setting PHOS flattening with name %s \n",arr->GetName())) ;
    if(fHarmonics==2){
      AliEPFlattener * h = (AliEPFlattener*)arr->At(0) ;  
      if(fTPCFlat) delete fTPCFlat ;
      fTPCFlat = new AliEPFlattener() ;
      fTPCFlat = h ;
      h = (AliEPFlattener*)arr->At(1) ;  
      if(fV0AFlat) delete fV0AFlat ;
      fV0AFlat = new AliEPFlattener() ;
      fV0AFlat = h ;
      h = (AliEPFlattener*)arr->At(2) ;  
      if(fV0CFlat) delete fV0CFlat ;
      fV0CFlat = new AliEPFlattener() ;
      fV0CFlat = h ;
      
      h = (AliEPFlattener*)arrQ->At(0) ;  
      if(fTPCQFlat) delete fTPCQFlat ;
      fTPCQFlat = new AliEPFlattener() ;
      fTPCQFlat = h ;
      h = (AliEPFlattener*)arrQ->At(1) ;  
      if(fV0AQFlat) delete fV0AQFlat ;
      fV0AQFlat = new AliEPFlattener() ;
      fV0AQFlat = h ;
      h = (AliEPFlattener*)arrQ->At(2) ;  
      if(fV0CQFlat) delete fV0CQFlat ;
      fV0CQFlat = new AliEPFlattener() ;
      fV0CQFlat = h ;
    }    
    if(fHarmonics==3){
      AliEPFlattener * h = (AliEPFlattener*)arr->At(3) ;  
      if(fTPCFlat) delete fTPCFlat ;
      fTPCFlat = new AliEPFlattener() ;
      fTPCFlat = h ;
      h = (AliEPFlattener*)arr->At(4) ;  
      if(fV0AFlat) delete fV0AFlat ;
      fV0AFlat = new AliEPFlattener() ;
      fV0AFlat = h ;
      h = (AliEPFlattener*)arr->At(5) ;  
      if(fV0CFlat) delete fV0CFlat ;
      fV0CFlat = new AliEPFlattener() ;
      fV0CFlat = h ;
      
      h = (AliEPFlattener*)arrQ->At(3) ;  
      if(fTPCQFlat) delete fTPCQFlat ;
      fTPCQFlat = new AliEPFlattener() ;
      fTPCQFlat = h ;
      h = (AliEPFlattener*)arrQ->At(4) ;  
      if(fV0AQFlat) delete fV0AQFlat ;
      fV0AQFlat = new AliEPFlattener() ;
      fV0AQFlat = h ;
      h = (AliEPFlattener*)arrQ->At(5) ;  
      if(fV0CQFlat) delete fV0CQFlat ;
      fV0CQFlat = new AliEPFlattener() ;
      fV0CQFlat = h ;
    }
      
    // TPC Event Plane Weights
    AliOADBContainer *fEPContainer=NULL;
    TString oadbfilename = (Form("%s/COMMON/EVENTPLANE/data/epphidist.root", AliAnalysisManager::GetOADBPath()));

    TFile foadb(oadbfilename);
    if(!foadb.IsOpen()) AliFatal(Form("Cannot open OADB file %s", oadbfilename.Data()));

    AliInfo("Using Standard OADB");
    fEPContainer = (AliOADBContainer*) foadb.Get("epphidist");
    if (!fEPContainer) AliFatal("Cannot fetch OADB container for EP selection");
    foadb.Close();

    fPhiDist = (TH1F*) fEPContainer->GetObject(run, "Default");

    Bool_t emptybins;
    int iter = 0;
    while (iter<3){
      emptybins = kFALSE;
      for (int i=1; i<fPhiDist->GetNbinsX(); i++){
         if (!((fPhiDist->GetBinContent(i))>0)) {
            emptybins = kTRUE;
         }
      }
      if (emptybins) {
        fPhiDist->Rebin();
        iter++;
      }
      else iter = 3;
   }

      
      
      
    fEventCounter++ ;
  }

  
  // Checks if we have a primary vertex
  // Get primary vertices form ESD
  const AliAODVertex *esdVertex5 = event->GetPrimaryVertex();

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
  else if(fCentrality<30.)
    fCenBin=3 ;
  else if(fCentrality<40.)
    fCenBin=4 ;
  else if(fCentrality<60.)
    fCenBin=5 ;
  else 
    fCenBin=6 ;

  if(TestPHOSEvent(event)){
    //bad event
    FillHistogram("hBadCentrality",fCentrality) ;   
    PostData(1, fOutputContainer);
    return;
  }
  
  
  
  Double_t cWeight=CentralityWeight(fCentrality) ;

  //Calculate EP resolutions from centrality
  EvalResolution() ;
  EvalQResolution() ;
  
  EvalV0ReactionPlane(event) ;

  
  //TPC evals second order RP need sign of V0 RP to find direction
  //reaction plane
  Double_t rpFull=0.,dPsi=0. ;
  
  fHaveTPCRP=GetTPCEventPlane(rpFull,dPsi);
   
  if(fHaveTPCRP){
    //rpFull is second order EP, in the range 0-pi
    //Set correct direction
    while(rpFull<0)rpFull+=TMath::TwoPi()/fHarmonics ;
    while(rpFull>TMath::TwoPi()/fHarmonics)rpFull-=(TMath::TwoPi()/fHarmonics) ;
    FillHistogram("phiRP",rpFull,fCentrality) ;    
    //Apply flattening
    fRP = fTPCFlat->MakeFlat(rpFull,fCentrality) ;
    while(fRP<0)fRP+=TMath::TwoPi()/fHarmonics ;
    while(fRP>TMath::TwoPi()/fHarmonics)fRP-=(TMath::TwoPi()/fHarmonics) ;
    FillHistogram("phiRPQ",fRP,fCentrality,fQTPC) ; //Yes, there is no difference between RP and RPQ so far 
    FillHistogram("phiRPflat",fRP,fCentrality) ;      
  }

  ApplyFinalQFlattening() ;
  
  fCutsV0->SetEvent(event,0x0);
  fFlowEvent->ClearFast();
  fFlowEvent->Fill(fCutsV0, fCutsTPC);
  fFlowEvent->TagSubeventsInEta(-10.,-1.,1.,10.) ;
  fFlowEvent->SetReferenceMultiplicity(event->GetNumberOfTracks());
  fFlowEvent->DefineDeadZone(0., 0, 0, 0);
  AliFlowVector qArray[2] ;
  fFlowEvent->Get2Qsub(qArray,fHarmonics) ;
   
  
  Double_t x= qArray[1].Phi()/Double_t(fHarmonics) ;
  while(x<0)x+=TMath::TwoPi()/fHarmonics ;
  while(x>TMath::TwoPi()/fHarmonics)x-=TMath::TwoPi()/fHarmonics ;
  Double_t y= qArray[0].Phi()/Double_t(fHarmonics) ;
  while(y<0)y+=TMath::TwoPi()/fHarmonics ;
  while(y>TMath::TwoPi()/fHarmonics)y-=TMath::TwoPi()/fHarmonics ;
  FillHistogram("phiRPV0AFlow",x,fCentrality) ;
  FillHistogram("phiRPV0CFlow",y,fCentrality) ;
  
  FillHistogram("phiV0ACorrel",x,fRPV0A,fCentrality) ;
  FillHistogram("phiV0CCorrel",y,fRPV0C,fCentrality) ;
    
  FillHistogram("phiRPV0Aflat",fRPV0A,fCentrality) ;
  FillHistogram("phiRPV0Cflat",fRPV0C,fCentrality) ;
  
  FillHistogram("phiRPV0AQflat",fRPQV0A,fCentrality,fQV0A) ;
  FillHistogram("phiRPV0CQflat",fRPQV0C,fCentrality,fQV0C) ;

  FillHistogram("phiRPV0AC",fRPV0A,fRPV0C,fCentrality) ;

  FillHistogram("cos2V0ARP",TMath::Cos(fHarmonics*fRPQV0A),fCentrality,fQV0A) ;
  FillHistogram("sin2V0ARP",TMath::Sin(fHarmonics*fRPQV0A),fCentrality,fQV0A) ;
  FillHistogram("cos2V0CRP",TMath::Cos(fHarmonics*fRPQV0C),fCentrality,fQV0C) ;
  FillHistogram("sin2V0CRP",TMath::Sin(fHarmonics*fRPQV0C),fCentrality,fQV0C) ;
  FillHistogram("QV0A",fCentrality,fQV0A) ;
  FillHistogram("QV0C",fCentrality,fQV0C) ;
  FillHistogram("QTPC",fCentrality,fQTPC) ;
  
  if(fHaveTPCRP){
    FillHistogram("phiRPQflat",fRPQ,fCentrality,fQTPC) ;  
    FillHistogram("cos2TPCAC",TMath::Cos(fHarmonics*dPsi),fCentrality) ;
    FillHistogram("qcos2TPCAC",TMath::Cos(fHarmonics*dPsi),fCentrality,fQTPC) ;
    FillHistogram("cos2TPCRP",TMath::Cos(fHarmonics*fRP),fCentrality,fQTPC) ;
    FillHistogram("sin2TPCRP",TMath::Sin(fHarmonics*fRP),fCentrality,fQTPC) ;    

    Double_t cosAC = TMath::Cos(fHarmonics*(fRPV0A-fRPV0C)) ;
    Double_t cosAT = TMath::Cos(fHarmonics*(fRP-fRPV0A)) ;
    Double_t cosCT = TMath::Cos(fHarmonics*(fRP-fRPV0C)) ;
    Double_t qcosAC = TMath::Cos(fHarmonics*(fRPQV0A-fRPQV0C)) ;
    Double_t qcosAT = TMath::Cos(fHarmonics*(fRPQ-fRPQV0A)) ;
    Double_t qcosCT = TMath::Cos(fHarmonics*(fRPQ-fRPQV0C)) ;

    FillHistogram("phiRPV0ATPC",fRP,fRPV0A,fCentrality,fQV0A*fQTPC) ;

    FillHistogram("cos2V0AC",fCentrality,cosAC) ;
    FillHistogram("qcos2V0AC",fCentrality,qcosAC,fQV0A*fQV0C) ;
    
    FillHistogram("cos2V0ATPC",fCentrality,cosAT) ;
    FillHistogram("qcos2V0ATPC",fCentrality,qcosAT,fQV0A*fQTPC) ;

    FillHistogram("cos2V0CTPC",fCentrality,cosCT) ;
    FillHistogram("qcos2V0CTPC",fCentrality,qcosCT,fQV0C*fQTPC) ;

    FillHistogram("phiRPV0CTPC",fRP,fRPV0C,fCentrality,fQV0C*fQTPC) ;
    
    //calculate event-by-event resolutions
    if(cosAC!=0. && cosAT!=0. && cosCT!=0.){
       Double_t res = TMath::Sqrt(TMath::Abs(cosAC*cosCT/cosAT));
       FillHistogram("ResV0C",fCentrality,res) ; 
       FillHistogram("InvResV0C",fCentrality,1./res) ; 
        
       res = TMath::Sqrt(TMath::Abs(cosAC*cosAT/cosCT));
       FillHistogram("ResV0A",fCentrality,res) ; 
       FillHistogram("InvResV0A",fCentrality,1./res) ; 

       res = TMath::Sqrt(TMath::Abs(cosAT*cosCT/cosAC));
       FillHistogram("ResTPC",fCentrality,res) ; 
       FillHistogram("InvResTPC",fCentrality,1./res) ; 
    }
  }     
 
  FillHistogram("hSelEvents",4.5,fRunNumber-0.5) ;
  FillHistogram("hTotSelEvents",4.5) ;
  //All event selections done
  FillHistogram("hCentrality",fCentrality,fRunNumber-0.5) ;
  FillHistogram("hCentralityCorr",fCentrality,cWeight) ;
  //Reaction plane is defined in the range (0;pi)
  //We have 10 bins
  Double_t averageRP = fRPV0A+fRPV0C+fRP ;
  if(fHaveTPCRP)
    averageRP/=3.;
  else
    averageRP/=2.;
  Int_t irp=Int_t(10.*averageRP*fHarmonics/TMath::TwoPi());
  if(irp>9)irp=9 ;
  
  if(!fPHOSEvents[zvtx][fCenBin][irp]) 
    fPHOSEvents[zvtx][fCenBin][irp]=new TList() ;
  TList * prevPHOS = fPHOSEvents[zvtx][fCenBin][irp] ;

  
  
  //  ProcessMC() ;

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

  AliAODCaloCells * cells = event->GetPHOSCells() ;
  FillHistogram("hCenTrack",fCentrality,event->GetNumberOfTracks()) ;
  
  //QA PHOS cells
  Int_t nCellModule[3] = {0,0,0};
  Int_t nCells=cells->GetNumberOfCells();
  for (Int_t iCell=0; iCell<nCells; iCell++) {
    Int_t cellAbsId = cells->GetCellNumber(iCell);
    Int_t relId[4] ;
    fPHOSGeo->AbsToRelNumbering(cellAbsId,relId);
    Int_t mod1  = relId[0];
    Int_t cellX = relId[2];
    Int_t cellZ = relId[3] ;
    Float_t energy = cells->GetAmplitude(iCell);
    FillHistogram("hCellEnergy",energy);
    if(mod1==1) {
      nCellModule[0]++;
      FillHistogram("hCellEnergyM1",cells->GetAmplitude(iCell));
      FillHistogram("hCellNXZM1",cellX,cellZ,1.);
      FillHistogram("hCellEXZM1",cellX,cellZ,energy);
    }
    else if (mod1==2) {
      nCellModule[1]++;
      FillHistogram("hCellEnergyM2",cells->GetAmplitude(iCell));
      FillHistogram("hCellNXZM2",cellX,cellZ,1.);
      FillHistogram("hCellEXZM2",cellX,cellZ,energy);
    }
    else if (mod1==3) {
      nCellModule[2]++;
      FillHistogram("hCellEnergyM3",cells->GetAmplitude(iCell));
      FillHistogram("hCellNXZM3",cellX,cellZ,1.);
      FillHistogram("hCellEXZM3",cellX,cellZ,energy);
    }
  }
  FillHistogram("hCellMultEventM1",nCellModule[0]);
  FillHistogram("hCellMultEventM2",nCellModule[1]);
  FillHistogram("hCellMultEventM3",nCellModule[2]);
  FillHistogram("hCenPHOSCells",fCentrality,nCells) ;
  FillHistogram("hCenPHOSCellsM12",fCentrality,nCellModule[0],nCellModule[1]) ;
  FillHistogram("hCenPHOSCellsM13",fCentrality,nCellModule[0],nCellModule[2]) ;
  FillHistogram("hCenPHOSCellsM23",fCentrality,nCellModule[1],nCellModule[2]) ;
  
  
  
  TVector3 localPos ;
  Int_t nAll=0,nCPV=0, nDisp=0, nBoth=0; 
  for (Int_t i=0; i<multClust; i++) {
    AliAODCaloCluster *clu = event->GetCaloCluster(i);
    if ( !clu->IsPHOS() || clu->E()<0.3) continue;
    
    if(fDistCut && (clu->GetDistanceToBadChannel()<2.5))
      continue ;
 
    Float_t  position[3];
    clu->GetPosition(position);
    TVector3 global(position) ;
    Int_t relId[4] ;
    fPHOSGeo->GlobalPos2RelId(global,relId) ;
    Int_t mod  = relId[0] ;
    Int_t cellX = relId[2];
    Int_t cellZ = relId[3] ;
    
    
    FillHistogram("hCluEvsClu",clu->E(),clu->GetNCells()) ;
    FillHistogram("hCluEvsCluM",clu->E(),clu->GetM02()) ;

    FillHistogram(Form("hTofM%d",mod),clu->E(),clu->GetTOF()) ;
    if(clu->E()>1.)
      FillHistogram("hCenTOF",fCentrality,clu->GetTOF()) ;
    if((clu->GetTOF()>1.5e-7) || (clu->GetTOF() <-2.5e-7) )
      continue ;
    
    
    //Apply re-Calibreation
    if(clu->E()<0.3) continue;
    if(clu->GetNCells()<3) continue ;
    if(clu->GetM02()<0.2)   continue ;
    
//    Double_t ecross=EvalEcross(&cluPHOS1) ;
//    if(ecross>0.98) continue ;          
  
/*    
    //correct misalignment
    const Float_t shiftX[6]={0.,-2.3,-2.11,-1.53,0.,0.} ;
    const Float_t shiftZ[6]={0.,-0.4, 0.52, 0.8,0.,0.} ;
    fPHOSGeo->Global2Local(localPos,global,mod) ;
    fPHOSGeo->Local2Global(mod,localPos.X()+shiftX[mod],localPos.Z()+shiftZ[mod],global);  
    position[0]=global.X() ;
    position[1]=global.Y() ;
    position[2]=global.Z() ;
    cluPHOS1.SetPosition(position) ; 
*/   
      
    TLorentzVector pv1 ;
    clu->GetMomentum(pv1 ,vtx0);

//    if(mod==2) pv1*=135.5/134.0 ;
//    if(mod==3) pv1*=135.5/137.2 ;
    
    Double_t ecore=clu->GetCoreEnergy() ; 
//    if(mod==2) ecore*=135.5/134.0 ;
//    if(mod==3) ecore*=135.5/137.2 ;

    
    FillHistogram(Form("hCluLowM%d",mod),cellX,cellZ,1.);
    if(clu->E()>1.5){
      FillHistogram(Form("hCluHighM%d",mod),cellX,cellZ,1.);
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
    
//  Cut on FullLamdas
    ph->SetDispBit(clu->GetDispersion()<2.5*2.5) ;
    ph->SetDisp2Bit(clu->Chi2()<2.5*2.5) ;

    Double_t distBC=clu->GetDistanceToBadChannel();
    if(distBC>3.){
      FillHistogram(Form("hPhotAll_DistBad2_cen%d",fCenBin),ph->Pt(),cWeight) ;
      FillHistogram(Form("hPhotAll_DistBad2core_cen%d",fCenBin),ph->GetMomV2()->Pt(),cWeight) ;
      if(ph->IsDisp2OK()){
        FillHistogram(Form("hPhotAll_DistBad2Disp_cen%d",fCenBin),ph->Pt(),cWeight) ;
        FillHistogram(Form("hPhotAll_DistBad2Dispcore_cen%d",fCenBin),ph->GetMomV2()->Pt(),cWeight) ;	
      }
      if(distBC>4.){
        FillHistogram(Form("hPhotAll_DistBad4_cen%d",fCenBin),ph->Pt(),cWeight) ;
        FillHistogram(Form("hPhotAll_DistBad4core_cen%d",fCenBin),ph->GetMomV2()->Pt(),cWeight) ;
        if(ph->IsDisp2OK()){
          FillHistogram(Form("hPhotAll_DistBad4Disp_cen%d",fCenBin),ph->Pt(),cWeight) ;
          FillHistogram(Form("hPhotAll_DistBad4Dispcore_cen%d",fCenBin),ph->GetMomV2()->Pt(),cWeight) ;	
        }
        if(distBC>6.){
          FillHistogram(Form("hPhotAll_DistBad6_cen%d",fCenBin),ph->Pt(),cWeight) ;
          FillHistogram(Form("hPhotAll_DistBad6core_cen%d",fCenBin),ph->GetMomV2()->Pt(),cWeight) ;
          if(ph->IsDisp2OK()){
            FillHistogram(Form("hPhotAll_DistBad6Disp_cen%d",fCenBin),ph->Pt(),cWeight) ;
            FillHistogram(Form("hPhotAll_DistBad6Dispcore_cen%d",fCenBin),ph->GetMomV2()->Pt(),cWeight) ;	
          }	  
	}    
      }
    }
    if(ph->IsDispOK()){
      FillHistogram(Form("hCluDispM%d",mod),cellX,cellZ,1.);
    }
    //Track matching
    Bool_t cpvBit=(clu->GetEmcCpvDistance()>2.5) ;
    ph->SetCPVBit(clu->GetEmcCpvDistance()>2.5) ;   
//    ph->SetCPVBit2(cpvBit2) ;
    if(cpvBit){
      FillHistogram(Form("hCluVetoM%d",mod),cellX,cellZ,1.);
    }
    
    
    ph->SetEMCx(float(cellX)) ;
    ph->SetEMCz(float(cellZ)) ;
    ph->SetLambdas(clu->GetM20(),clu->GetM02()) ;
    ph->SetUnfolded(clu->GetNExMax()<2); // Remember, if it is unfolded          
    nAll++ ;
    if(cpvBit)nCPV++ ;
    if(ph->IsDisp2OK())nDisp++ ;
    if(cpvBit && ph->IsDisp2OK()) nBoth++ ;
    inPHOS++ ;
  }
  
  FillHistogram("hCenPHOS",fCentrality,inPHOS,cWeight) ;
  if(inPHOS==0){
    PostData(1, fOutputContainer);
    fEventCounter++;
    return ; 
  }
  
  //Parameterization of the mean number of PHOS clusters vs centrality
  //This is 1/<Nclu> parameterization
  
  
  Double_t wA=1./fV0Ares ;
  Double_t wC=1./fV0Cres ;
  Double_t wT=1./fTPCres ;
  Double_t wQA=1./fV0AQres ;
  Double_t wQC=1./fV0CQres ;
  Double_t wQT=1./fTPCQres ;
  Double_t nCorr=1./PHOSMultiplicity() ; //multiplicity vs. centrality
  
   FillHistogram("resV0A",fCentrality,fV0Ares,cWeight) ;
   FillHistogram("resV0C",fCentrality,fV0Cres,cWeight) ;
   FillHistogram("resTPC",fCentrality,fTPCres,cWeight) ;
   FillHistogram("qresV0A",fCentrality,fV0AQres,cWeight) ;
   FillHistogram("qresV0C",fCentrality,fV0CQres,cWeight) ;
   FillHistogram("qresTPC",fCentrality,fTPCQres,cWeight) ;
   
   TVector2 vQmA(fQV0A*TMath::Cos(fHarmonics*fRPQV0A),fQV0A*TMath::Sin(fHarmonics*fRPQV0A)) ;
   TVector2 vQmB(fQV0C*TMath::Cos(fHarmonics*fRPQV0C),fQV0C*TMath::Sin(fHarmonics*fRPQV0C)) ;
   
   FillHistogram("VaVbcen",fCentrality,vQmA*vQmB,cWeight) ;
   
  for (Int_t i1=0; i1<inPHOS; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1) ;
    
    Double_t dphiA=ph1->Phi()-fRPV0A ;
//    while(dphiA<0)dphiA+=TMath::TwoPi()/fHarmonics ;
//    while(dphiA>TMath::TwoPi()/fHarmonics)dphiA-=TMath::TwoPi()/fHarmonics ;
    Double_t cosA=wA*TMath::Cos(fHarmonics*dphiA) ;
    Double_t qcosA=wQA*TMath::Cos(fHarmonics*(ph1->Phi()-fRPQV0A)) ;
    
    Double_t dphiC=ph1->Phi()-fRPV0C ;
//    while(dphiC<0)dphiC+=TMath::TwoPi()/fHarmonics ;
//    while(dphiC>TMath::TwoPi()/fHarmonics)dphiC-=TMath::TwoPi()/fHarmonics ;
    Double_t cosC=wC*TMath::Cos(fHarmonics*dphiC) ;
    Double_t qcosC=wQC*TMath::Cos(fHarmonics*(ph1->Phi()-fRPQV0C)) ;

    Double_t dphiT=ph1->Phi()-fRP ;
//    while(dphiT<0)dphiT+=TMath::TwoPi()/fHarmonics ;
//    while(dphiT>TMath::TwoPi()/fHarmonics)dphiT-=TMath::TwoPi()/fHarmonics ;
    Double_t cosT=wT*TMath::Cos(fHarmonics*dphiT) ;
    Double_t qcosT=wQT*TMath::Cos(fHarmonics*(ph1->Phi()-fRPQ)) ;
         
    Double_t pt=ph1->Pt() ;
    Double_t ptV=ph1->GetMomV2()->Pt() ;
    
    //Official SP
    //calculate vU
    TVector2 vU;
    Double_t dUX = TMath::Cos(fHarmonics*ph1->Phi());
    Double_t dUY = TMath::Sin(fHarmonics*ph1->Phi());
    vU.Set(dUX,dUY);

    Double_t dUQA = vU*vQmA;
    Double_t dUQB = vU*vQmB;
    

    FillHistogram(Form("hPhotPhiSPV0AAll"),pt,fCentrality,dUQA,cWeight) ;	
    FillHistogram(Form("hPhotPhiSPV0CAll"),pt,fCentrality,dUQB,cWeight) ;	
    FillHistogram(Form("hPhotPhiSPV0AAllcore"),ptV,fCentrality,dUQA,cWeight) ;	
    FillHistogram(Form("hPhotPhiSPV0CAllcore"),ptV,fCentrality,dUQB,cWeight) ;	
    
    FillHistogram(Form("hPhotPhiR2V0AAll_cen%d",fCenBin),pt,qcosA,fQV0A*cWeight*nCorr) ;
    FillHistogram(Form("hPhotPhiRV0AAll_cen%d", fCenBin),pt,qcosA,fQV0A*cWeight) ;
    FillHistogram(Form("hPhotPhiR2V0CAll_cen%d",fCenBin),pt,qcosC,fQV0C*cWeight*nCorr) ;
    FillHistogram(Form("hPhotPhiRV0CAll_cen%d", fCenBin),pt,qcosC,fQV0C*cWeight) ;
    if(fHaveTPCRP){
      FillHistogram(Form("hPhotPhiR2TPCAll_cen%d",fCenBin),pt,qcosT,fQTPC*cWeight*nCorr) ;
      FillHistogram(Form("hPhotPhiRTPCAll_cen%d",fCenBin), pt,qcosT,fQTPC*cWeight) ;
    }
    FillHistogram(Form("hPhotPhiR2V0AAllcore_cen%d",fCenBin),ptV,qcosA,fQV0A*cWeight*nCorr) ;
    FillHistogram(Form("hPhotPhiRV0AAllcore_cen%d", fCenBin),ptV,qcosA,fQV0A*cWeight) ;
    FillHistogram(Form("hPhotPhiR2V0CAllcore_cen%d",fCenBin),ptV,qcosC,fQV0C*cWeight*nCorr) ;
    FillHistogram(Form("hPhotPhiRV0CAllcore_cen%d", fCenBin),ptV,qcosC,fQV0C*cWeight) ;
    if(fHaveTPCRP){
      FillHistogram(Form("hPhotPhiR2TPCAllcore_cen%d",fCenBin),ptV,qcosT,fQTPC*cWeight*nCorr) ;
      FillHistogram(Form("hPhotPhiRTPCAllcore_cen%d",fCenBin), ptV,qcosT,fQTPC*cWeight) ;
    }
    
    FillHistogram(Form("hPhotcosV0AAll_cen%d",fCenBin),pt,cosA,cWeight) ;
    FillHistogram(Form("hPhotcosV0CAll_cen%d",fCenBin),pt,cosC,cWeight) ;
    FillHistogram(Form("hPhotcosNV0AAll_cen%d",fCenBin),pt,cosA,cWeight*nCorr) ;
    FillHistogram(Form("hPhotcosNV0CAll_cen%d",fCenBin),pt,cosC,cWeight*nCorr) ;
    FillHistogram(Form("hPhotcosV0AAllcore_cen%d",fCenBin),ptV,cosA,cWeight) ;
    FillHistogram(Form("hPhotcosV0CAllcore_cen%d",fCenBin),ptV,cosC,cWeight) ;
    FillHistogram(Form("hPhotcosNV0AAllcore_cen%d",fCenBin),ptV,cosA,cWeight*nCorr) ;
    FillHistogram(Form("hPhotcosNV0CAllcore_cen%d",fCenBin),ptV,cosC,cWeight*nCorr) ;
    if(fHaveTPCRP){
      FillHistogram(Form("hPhotcosTPCAll_cen%d",fCenBin),pt,cosT,cWeight) ;
      FillHistogram(Form("hPhotcosTPCAllcore_cen%d",fCenBin),ptV,cosT,cWeight) ;
      FillHistogram(Form("hPhotcosNTPCAll_cen%d",fCenBin),pt,cosT,cWeight*nCorr) ;
      FillHistogram(Form("hPhotcosNTPCAllcore_cen%d",fCenBin),ptV,cosT,cWeight*nCorr) ;
    }
    FillHistogram(Form("hPhotAll_cen%d",fCenBin),pt,cWeight) ;
    FillHistogram(Form("hPhotRaAll_cen%d",fCenBin),pt,cWeight*wA) ;
    FillHistogram(Form("hPhotRcAll_cen%d",fCenBin),pt,cWeight*wC) ;
    FillHistogram(Form("hPhotRNaAll_cen%d",fCenBin),pt,cWeight*wA*nCorr) ;
    FillHistogram(Form("hPhotRNcAll_cen%d",fCenBin),pt,cWeight*wC*nCorr) ;
    FillHistogram(Form("hPhotAll_m%d_cen%d",ph1->Module(),fCenBin),pt,cWeight) ;
    FillHistogram(Form("hPhotAllcore_m%d_cen%d",ph1->Module(),fCenBin),ptV,cWeight) ;   
    FillHistogram(Form("hPhotAllcore_cen%d",fCenBin),ptV,cWeight) ;
    FillHistogram(Form("hPhotRaAllcore_cen%d",fCenBin),ptV,cWeight*wA) ;
    FillHistogram(Form("hPhotRcAllcore_cen%d",fCenBin),ptV,cWeight*wC) ;
    if(fHaveTPCRP){
      FillHistogram(Form("hPhotRtAll_cen%d",fCenBin),pt,cWeight*wT) ;
      FillHistogram(Form("hPhotRtAllcore_cen%d",fCenBin),ptV,cWeight*wT) ;
      FillHistogram(Form("hPhotRNtAll_cen%d",fCenBin),pt,cWeight*wT*nCorr) ;
      FillHistogram(Form("hPhotRNtAllcore_cen%d",fCenBin),ptV,cWeight*wT*nCorr) ;
    }
    //CPV----
    if(ph1->IsCPVOK()){
      
      FillHistogram(Form("hPhotPhiSPV0ACPV"),pt,fCentrality,dUQA,cWeight) ;	
      FillHistogram(Form("hPhotPhiSPV0CCPV"),pt,fCentrality,dUQB,cWeight) ;	
      FillHistogram(Form("hPhotPhiSPV0ACPVcore"),ptV,fCentrality,dUQA,cWeight) ;	
      FillHistogram(Form("hPhotPhiSPV0CCPVcore"),ptV,fCentrality,dUQB,cWeight) ;	
      
      FillHistogram(Form("hPhotPhiR2V0ACPV_cen%d",fCenBin),pt,qcosA,fQV0A*cWeight*nCorr) ;
      FillHistogram(Form("hPhotPhiRV0ACPV_cen%d", fCenBin),pt,qcosA,fQV0A*cWeight) ;
      FillHistogram(Form("hPhotPhiR2V0CCPV_cen%d",fCenBin),pt,qcosC,fQV0C*cWeight*nCorr) ;
      FillHistogram(Form("hPhotPhiRV0CCPV_cen%d", fCenBin),pt,qcosC,fQV0C*cWeight) ;
      if(fHaveTPCRP){
        FillHistogram(Form("hPhotPhiR2TPCCPV_cen%d",fCenBin),pt,qcosT,fQTPC*cWeight*nCorr) ;
        FillHistogram(Form("hPhotPhiRTPCCPV_cen%d", fCenBin),pt,qcosT,fQTPC*cWeight) ;
      }
      FillHistogram(Form("hPhotPhiR2V0ACPVcore_cen%d",fCenBin),ptV,qcosA,fQV0A*cWeight*nCorr) ;
      FillHistogram(Form("hPhotPhiRV0ACPVcore_cen%d", fCenBin),ptV,qcosA,fQV0A*cWeight) ;
      FillHistogram(Form("hPhotPhiR2V0CCPVcore_cen%d",fCenBin),ptV,qcosC,fQV0C*cWeight*nCorr) ;
      FillHistogram(Form("hPhotPhiRV0CCPVcore_cen%d", fCenBin),ptV,qcosC,fQV0C*cWeight) ;
      if(fHaveTPCRP){
        FillHistogram(Form("hPhotPhiR2TPCCPVcore_cen%d",fCenBin),ptV,qcosT,fQTPC*cWeight*nCorr) ;
        FillHistogram(Form("hPhotPhiRTPCCPVcore_cen%d",fCenBin), ptV,qcosT,fQTPC*cWeight) ;
      }
      FillHistogram(Form("hPhotcosV0ACPV_cen%d",fCenBin),pt,cosA,cWeight) ;
      FillHistogram(Form("hPhotcosV0CCPV_cen%d",fCenBin),pt,cosC,cWeight) ;
      FillHistogram(Form("hPhotcosNV0ACPV_cen%d",fCenBin),pt,cosA,cWeight*nCorr) ;
      FillHistogram(Form("hPhotcosNV0CCPV_cen%d",fCenBin),pt,cosC,cWeight*nCorr) ;
      FillHistogram(Form("hPhotcosV0ACPVcore_cen%d",fCenBin),ptV,cosA,cWeight) ;
      FillHistogram(Form("hPhotcosV0CCPVcore_cen%d",fCenBin),ptV,cosC,cWeight) ;
      FillHistogram(Form("hPhotcosNV0ACPVcore_cen%d",fCenBin),ptV,cosA,cWeight*nCorr) ;
      FillHistogram(Form("hPhotcosNV0CCPVcore_cen%d",fCenBin),ptV,cosC,cWeight*nCorr) ;
      if(fHaveTPCRP){
        FillHistogram(Form("hPhotcosTPCCPV_cen%d",fCenBin),pt,cosT,cWeight) ;
        FillHistogram(Form("hPhotcosTPCCPVcore_cen%d",fCenBin),ptV,cosT,cWeight) ;
        FillHistogram(Form("hPhotcosNTPCCPV_cen%d",fCenBin),pt,cosT,cWeight*nCorr) ;
        FillHistogram(Form("hPhotcosNTPCCPVcore_cen%d",fCenBin),ptV,cosT,cWeight*nCorr) ;
      }
      
      FillHistogram(Form("hPhotCPV_cen%d",fCenBin),pt,cWeight) ;
      FillHistogram(Form("hPhotRaCPV_cen%d",fCenBin),pt,cWeight*wA) ;
      FillHistogram(Form("hPhotRcCPV_cen%d",fCenBin),pt,cWeight*wC) ;
      FillHistogram(Form("hPhotRtCPV_cen%d",fCenBin),pt,cWeight*wT) ;
      FillHistogram(Form("hPhotRNaCPV_cen%d",fCenBin),pt,cWeight*wA*nCorr) ;
      FillHistogram(Form("hPhotRNcCPV_cen%d",fCenBin),pt,cWeight*wC*nCorr) ;
      FillHistogram(Form("hPhotRNtCPV_cen%d",fCenBin),pt,cWeight*wT*nCorr) ;
      FillHistogram(Form("hPhotCPVcore_cen%d",fCenBin),ptV,cWeight) ;
      FillHistogram(Form("hPhotRaCPVcore_cen%d",fCenBin),ptV,cWeight*wA) ;
      FillHistogram(Form("hPhotRcCPVcore_cen%d",fCenBin),ptV,cWeight*wC) ;
      FillHistogram(Form("hPhotRtCPVcore_cen%d",fCenBin),ptV,cWeight*wT) ;
      FillHistogram(Form("hPhotRNaCPVcore_cen%d",fCenBin),ptV,cWeight*wA*nCorr) ;
      FillHistogram(Form("hPhotRNcCPVcore_cen%d",fCenBin),ptV,cWeight*wC*nCorr) ;
      FillHistogram(Form("hPhotRNtCPVcore_cen%d",fCenBin),ptV,cWeight*wT*nCorr) ;
      FillHistogram(Form("hPhotCPV_m%d_cen%d",ph1->Module(),fCenBin),pt,cWeight) ;
      FillHistogram(Form("hPhotCPVcore_m%d_cen%d",ph1->Module(),fCenBin),ptV,cWeight) ;
    }
//    if(ph1->IsCPV2OK()){
//      FillHistogram(Form("hPhotCPV2_cen%d",fCenBin),pt) ;
//    }
    //Disp2---
    if(ph1->IsDisp2OK()){

      FillHistogram(Form("hPhotPhiSPV0ADisp2"),pt,fCentrality,dUQA,cWeight) ;	
      FillHistogram(Form("hPhotPhiSPV0CDisp2"),pt,fCentrality,dUQB,cWeight) ;	
      FillHistogram(Form("hPhotPhiSPV0ADisp2core"),ptV,fCentrality,dUQA,cWeight) ;	
      FillHistogram(Form("hPhotPhiSPV0CDisp2core"),ptV,fCentrality,dUQB,cWeight) ;	

      FillHistogram(Form("hPhotPhiR2V0ADisp2_cen%d",fCenBin),pt,qcosA,fQV0A*cWeight*nCorr) ;
      FillHistogram(Form("hPhotPhiRV0ADisp2_cen%d", fCenBin),pt,qcosA,fQV0A*cWeight) ;
      FillHistogram(Form("hPhotPhiR2V0CDisp2_cen%d",fCenBin),pt,qcosC,fQV0C*cWeight*nCorr) ;
      FillHistogram(Form("hPhotPhiRV0CDisp2_cen%d", fCenBin),pt,qcosC,fQV0C*cWeight) ;
      if(fHaveTPCRP){
        FillHistogram(Form("hPhotPhiR2TPCDisp2_cen%d",fCenBin),pt,qcosT,fQTPC*cWeight*nCorr) ;
        FillHistogram(Form("hPhotPhiRTPCDisp2_cen%d",fCenBin), pt,qcosT,fQTPC*cWeight) ;
      }
      FillHistogram(Form("hPhotPhiR2V0ADisp2core_cen%d",fCenBin),ptV,qcosA,fQV0A*cWeight*nCorr) ;
      FillHistogram(Form("hPhotPhiRV0ADisp2core_cen%d", fCenBin),ptV,qcosA,fQV0A*cWeight) ;
      FillHistogram(Form("hPhotPhiR2V0CDisp2core_cen%d",fCenBin),ptV,qcosC,fQV0C*cWeight*nCorr) ;
      FillHistogram(Form("hPhotPhiRV0CDisp2core_cen%d", fCenBin),ptV,qcosC,fQV0C*cWeight) ;
      if(fHaveTPCRP){
        FillHistogram(Form("hPhotPhiR2TPCDisp2core_cen%d",fCenBin),ptV,qcosT,fQTPC*cWeight*nCorr) ;
        FillHistogram(Form("hPhotPhiRTPCDisp2core_cen%d",fCenBin), ptV,qcosT,fQTPC*cWeight) ;
      }

      FillHistogram(Form("hPhotcosV0ADisp2_cen%d",fCenBin),pt,cosA,cWeight) ;
      FillHistogram(Form("hPhotcosV0CDisp2_cen%d",fCenBin),pt,cosC,cWeight) ;
      FillHistogram(Form("hPhotcosNV0ADisp2_cen%d",fCenBin),pt,cosA,cWeight*nCorr) ;
      FillHistogram(Form("hPhotcosNV0CDisp2_cen%d",fCenBin),pt,cosC,cWeight*nCorr) ;
      FillHistogram(Form("hPhotcosV0ADisp2core_cen%d",fCenBin),ptV,cosA,cWeight) ;
      FillHistogram(Form("hPhotcosV0CDisp2core_cen%d",fCenBin),ptV,cosC,cWeight) ;
      FillHistogram(Form("hPhotcosNV0ADisp2core_cen%d",fCenBin),ptV,cosA,cWeight*nCorr) ;
      FillHistogram(Form("hPhotcosNV0CDisp2core_cen%d",fCenBin),ptV,cosC,cWeight*nCorr) ;
      if(fHaveTPCRP){
        FillHistogram(Form("hPhotcosTPCDisp2_cen%d",fCenBin),pt,cosT,cWeight) ;
        FillHistogram(Form("hPhotcosTPCDisp2core_cen%d",fCenBin),ptV,cosT,cWeight) ;
        FillHistogram(Form("hPhotcosNTPCDisp2_cen%d",fCenBin),pt,cosT,cWeight*nCorr) ;
        FillHistogram(Form("hPhotcosNTPCDisp2core_cen%d",fCenBin),ptV,cosT,cWeight*nCorr) ;
      }

      
      FillHistogram(Form("hPhotDisp2_cen%d",fCenBin),pt,cWeight) ;
      FillHistogram(Form("hPhotRaDisp2_cen%d",fCenBin),pt,wA*cWeight) ;
      FillHistogram(Form("hPhotRcDisp2_cen%d",fCenBin),pt,wC*cWeight) ;
      FillHistogram(Form("hPhotRNaDisp2_cen%d",fCenBin),pt,wA*cWeight*nCorr) ;
      FillHistogram(Form("hPhotRNcDisp2_cen%d",fCenBin),pt,wC*cWeight*nCorr) ;      
      FillHistogram(Form("hPhotDisp2core_cen%d",fCenBin),ptV,cWeight) ;
      FillHistogram(Form("hPhotRaDisp2core_cen%d",fCenBin),ptV,wA*cWeight) ;
      FillHistogram(Form("hPhotRcDisp2core_cen%d",fCenBin),ptV,wC*cWeight) ;
      FillHistogram(Form("hPhotRNaDisp2core_cen%d",fCenBin),ptV,wA*cWeight*nCorr) ;
      FillHistogram(Form("hPhotRNcDisp2core_cen%d",fCenBin),ptV,wC*cWeight*nCorr) ;
      if(fHaveTPCRP){
        FillHistogram(Form("hPhotRtDisp2_cen%d",fCenBin),pt,cWeight*wT) ;
        FillHistogram(Form("hPhotRtDisp2core_cen%d",fCenBin),ptV,cWeight*wT) ;
        FillHistogram(Form("hPhotRNtDisp2_cen%d",fCenBin),pt,cWeight*wT*nCorr) ;
        FillHistogram(Form("hPhotRNtDisp2core_cen%d",fCenBin),ptV,cWeight*wT*nCorr) ;
      }
      FillHistogram(Form("hPhotDisp2_m%d_cen%d",ph1->Module(),fCenBin),pt,cWeight) ;
      FillHistogram(Form("hPhotDisp2core_m%d_cen%d",ph1->Module(),fCenBin),ptV,cWeight) ;
      
      //Both2---
      if(ph1->IsCPVOK()){

	//----------Official SP----------------
        //fill the profile histograms
        FillHistogram(Form("hPhotPhiSPV0ABoth2"),pt,fCentrality,dUQA,cWeight) ;	
        FillHistogram(Form("hPhotPhiSPV0CBoth2"),pt,fCentrality,dUQB,cWeight) ;	
        FillHistogram(Form("hPhotPhiSPV0ABoth2core"),ptV,fCentrality,dUQA,cWeight) ;	
        FillHistogram(Form("hPhotPhiSPV0CBoth2core"),ptV,fCentrality,dUQB,cWeight) ;	
	//-------------------------------------
		
	FillHistogram(Form("hPhotPhiR2V0ABoth2_cen%d",fCenBin),pt,qcosA,fQV0A*cWeight*nCorr) ;
	FillHistogram(Form("hPhotPhiRV0ABoth2_cen%d", fCenBin),pt,qcosA,fQV0A*cWeight) ;
	FillHistogram(Form("hPhotPhiR2V0CBoth2_cen%d",fCenBin),pt,qcosC,fQV0C*cWeight*nCorr) ;
	FillHistogram(Form("hPhotPhiRV0CBoth2_cen%d", fCenBin),pt,qcosC,fQV0C*cWeight) ;
        if(fHaveTPCRP){
  	  FillHistogram(Form("hPhotPhiR2TPCBoth2_cen%d",fCenBin),pt,qcosT,fQTPC*cWeight) ;
  	  FillHistogram(Form("hPhotPhiRTPCBoth2_cen%d",fCenBin), pt,qcosT,fQTPC*cWeight) ;
	}
	FillHistogram(Form("hPhotPhiR2V0ABoth2core_cen%d",fCenBin),ptV,qcosA,fQV0A*cWeight*nCorr) ;
	FillHistogram(Form("hPhotPhiRV0ABoth2core_cen%d", fCenBin),ptV,qcosA,fQV0A*cWeight) ;
	FillHistogram(Form("hPhotPhiR2V0CBoth2core_cen%d",fCenBin),ptV,qcosC,fQV0C*cWeight*nCorr) ;
	FillHistogram(Form("hPhotPhiRV0CBoth2core_cen%d", fCenBin),ptV,qcosC,fQV0C*cWeight) ;
        if(fHaveTPCRP){
	  FillHistogram(Form("hPhotPhiR2TPCBoth2core_cen%d",fCenBin),ptV,qcosT,fQTPC*cWeight*nCorr) ;
	  FillHistogram(Form("hPhotPhiRTPCBoth2core_cen%d",fCenBin), ptV,qcosT,fQTPC*cWeight) ;
	}

        FillHistogram(Form("hPhotcosV0ABoth2_cen%d",fCenBin),pt,cosA,cWeight) ;
        FillHistogram(Form("hPhotcosV0CBoth2_cen%d",fCenBin),pt,cosC,cWeight) ;
        FillHistogram(Form("hPhotcosNV0ABoth2_cen%d",fCenBin),pt,cosA,cWeight*nCorr) ;
        FillHistogram(Form("hPhotcosNV0CBoth2_cen%d",fCenBin),pt,cosC,cWeight*nCorr) ;
        FillHistogram(Form("hPhotcosV0ABoth2core_cen%d",fCenBin),ptV,cosA,cWeight) ;
        FillHistogram(Form("hPhotcosV0CBoth2core_cen%d",fCenBin),ptV,cosC,cWeight) ;
        FillHistogram(Form("hPhotcosNV0ABoth2core_cen%d",fCenBin),ptV,cosA,cWeight*nCorr) ;
        FillHistogram(Form("hPhotcosNV0CBoth2core_cen%d",fCenBin),ptV,cosC,cWeight*nCorr) ;
        if(fHaveTPCRP){
          FillHistogram(Form("hPhotcosTPCBoth2_cen%d",fCenBin),pt,cosT,cWeight) ;
          FillHistogram(Form("hPhotcosTPCBoth2core_cen%d",fCenBin),ptV,cosT,cWeight) ;
          FillHistogram(Form("hPhotcosNTPCBoth2_cen%d",fCenBin),pt,cosT,cWeight*nCorr) ;
          FillHistogram(Form("hPhotcosNTPCBoth2core_cen%d",fCenBin),ptV,cosT,cWeight*nCorr) ;
        }
	
	FillHistogram(Form("hPhotBoth2_cen%d",fCenBin),pt,cWeight) ;
	FillHistogram(Form("hPhotRaBoth2_cen%d",fCenBin),pt,cWeight*wA) ;
	FillHistogram(Form("hPhotRcBoth2_cen%d",fCenBin),pt,cWeight*wC) ;
	FillHistogram(Form("hPhotRNaBoth2_cen%d",fCenBin),pt,cWeight*wA*nCorr) ;
	FillHistogram(Form("hPhotRNcBoth2_cen%d",fCenBin),pt,cWeight*wC*nCorr) ;
	FillHistogram(Form("hPhotBoth2core_cen%d",fCenBin),ptV,cWeight) ;
	FillHistogram(Form("hPhotRaBoth2core_cen%d",fCenBin),ptV,cWeight*wA) ;
	FillHistogram(Form("hPhotRcBoth2core_cen%d",fCenBin),ptV,cWeight*wC) ;
	FillHistogram(Form("hPhotRNaBoth2core_cen%d",fCenBin),ptV,cWeight*wA*nCorr) ;
	FillHistogram(Form("hPhotRNcBoth2core_cen%d",fCenBin),ptV,cWeight*wC*nCorr) ;
        if(fHaveTPCRP){
	  FillHistogram(Form("hPhotRtBoth2_cen%d",fCenBin),pt,cWeight*wT) ;
	  FillHistogram(Form("hPhotRtBoth2core_cen%d",fCenBin),ptV,cWeight*wT) ;
	  FillHistogram(Form("hPhotRNtBoth2_cen%d",fCenBin),pt,cWeight*wT*nCorr) ;
	  FillHistogram(Form("hPhotRNtBoth2core_cen%d",fCenBin),ptV,cWeight*wT*nCorr) ;
	}
        FillHistogram(Form("hPhotBoth2_m%d_cen%d",ph1->Module(),fCenBin),pt,cWeight) ;
        FillHistogram(Form("hPhotBoth2core_m%d_cen%d",ph1->Module(),fCenBin),ptV,cWeight) ;

      }
    }  
  }
    
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
void AliAnalysisTaskGammaFlow::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
  
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaFlow::FillHistogram(const char * key,Double_t x)const{
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
  if(tmp->IsA() == TClass::GetClass("TProfile")){
    ((TH1F*)tmp)->Fill(x,1.) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH1D")){
    ((TH1D*)tmp)->Fill(x) ;
    return ;
  }  
  AliInfo(Form("can not find 1D histogram <%s> ",key)) ;
}
//_____________________________________________________________________________
void AliAnalysisTaskGammaFlow::FillHistogram(const char * key,Double_t x,Double_t y)const{
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
  if(tmp->IsA() == TClass::GetClass("TProfile")){
    ((TProfile*)tmp)->Fill(x,y) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH2F")){
    ((TH2F*)tmp)->Fill(x,y) ;
    return ;
  }
  AliError(Form("Calling FillHistogram with 2 parameters for histo <%s> of type %s",key,tmp->IsA()->GetName())) ;
}

//_____________________________________________________________________________
void AliAnalysisTaskGammaFlow::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z) const{
  //Fills 1D histograms with Form(
  TObject * tmp = fOutputContainer->FindObject(key) ;
  if(!tmp){
    AliInfo(Form("can not find histogram <%s> ",key)) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH2F")){
    ((TH2F*)tmp)->Fill(x,y,z) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TProfile")){
    ((TProfile*)tmp)->Fill(x,y,z) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TProfile2D")){
    ((TProfile2D*)tmp)->Fill(x,y,z) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH3F")){
    ((TH3F*)tmp)->Fill(x,y,z) ;
    return ;
  }
  AliInfo(Form("can not find 2D/3D histogram <%s> ",key)) ;
}
//_____________________________________________________________________________
void AliAnalysisTaskGammaFlow::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z, Double_t w) const{
  //Fills 1D histograms with Form(
  TObject * tmp = fOutputContainer->FindObject(key) ;
  if(!tmp){
    AliInfo(Form("can not find histogram <%s> ",key)) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TH3F")){
    ((TH3F*)tmp)->Fill(x,y,z,w) ;
    return ;
  }
  if(tmp->IsA() == TClass::GetClass("TProfile2D")){
    ((TProfile2D*)tmp)->Fill(x,y,z,w) ;
    return ;
  }
  AliInfo(Form("can not find 3D histogram <%s> ",key)) ;
}

//___________________________________________________________________________
Int_t AliAnalysisTaskGammaFlow::ConvertRunNumber(Int_t run){

  switch(run){	
   case  170593 : return 1;
   case  170572 : return 2;
   case  170556 : return 3;
   case  170552 : return 4;
   case  170546 : return 5;
   case  170390 : return 6;
   case  170389 : return 7;
   case  170388 : return 8;
   case  170387 : return 9;
   case  170315 : return 10;
   case  170313 : return 11;
   case  170312 : return 12;
   case  170311 : return 13;
   case  170309 : return 14;
   case  170308 : return 15;
   case  170306 : return 16;
   case  170270 : return 17;
   case  170269 : return 18;
   case  170268 : return 19;
   case  170267 : return 20;
   case  170264 : return 21;
   case  170230 : return 22;
   case  170228 : return 23;
   case  170208 : return 24;
   case  170207 : return 25;
   case  170205 : return 26;
   case  170204 : return 27;
   case  170203 : return 28;
   case  170195 : return 29;
   case  170193 : return 30;
   case  170163 : return 31;
   case  170162 : return 32;
   case  170159 : return 33;
   case  170155 : return 34;
   case  170152 : return 35;
   case  170091 : return 36;
   case  170089 : return 37;
   case  170088 : return 38;
   case  170085 : return 39;
   case  170084 : return 40;
   case  170083 : return 41;
   case  170081 : return 42;
   case  170040 : return 43;
   case  170038 : return 44;
   case  170036 : return 45;
   case  170027 : return 46;
   case  169981 : return 47;
   case  169969 : return 48;
   case  169965 : return 49;
   case  169961 : return 50;
   case  169956 : return 51;
   case  169926 : return 52;
   case  169924 : return 53;
   case  169923 : return 54;
   case  169922 : return 55;
   case  169859 : return 56;
   case  169858 : return 57;
   case  169855 : return 58;
   case  169846 : return 59;
   case  169838 : return 60;
   case  169837 : return 61;
   case  169835 : return 62;
   case  169683 : return 63;
   case  169591 : return 64;
   case  169590 : return 65;
   case  169588 : return 66;
   case  169587 : return 67;
   case  169586 : return 68;
   case  169584 : return 69;
   case  169557 : return 70;
   case  169555 : return 71;
   case  169554 : return 72;
   case  169553 : return 73;
   case  169550 : return 74;
   case  169515 : return 75;
   case  169512 : return 76;
   case  169506 : return 77;
   case  169504 : return 78;
   case  169498 : return 79;
   case  169475 : return 80;
   case  169420 : return 81;
   case  169419 : return 82;
   case  169418 : return 83;
   case  169417 : return 84;
   case  169415 : return 85;
   case  169411 : return 86;
   case  169238 : return 87;
   case  169236 : return 88;
   case  169167 : return 89;
   case  169160 : return 90;
   case  169156 : return 91;
   case  169148 : return 92;
   case  169145 : return 93;
   case  169144 : return 94;
   case  169143 : return 95;
   case  169138 : return 96;
   case  169099 : return 97;
   case  169094 : return 98;
   case  169091 : return 99;
   case  169045 : return 100;
   case  169044 : return 101;
   case  169040 : return 102;
   case  169035 : return 103;
   case  168992 : return 104;
   case  168988 : return 105;
   case  168984 : return 106;
   case  168826 : return 107;
   case  168777 : return 108;
   case  168514 : return 109;
   case  168512 : return 110;
   case  168511 : return 111;
   case  168467 : return 112;
   case  168464 : return 113;
   case  168461 : return 114;
   case  168460 : return 115;
   case  168458 : return 116;
   case  168362 : return 117;
   case  168361 : return 118;
   case  168356 : return 119;
   case  168342 : return 120;
   case  168341 : return 121;
   case  168325 : return 122;
   case  168322 : return 123;
   case  168318 : return 124;
   case  168311 : return 125;
   case  168310 : return 126;
   case  168213 : return 127;
   case  168212 : return 128;
   case  168208 : return 129;
   case  168207 : return 130;
   case  168206 : return 131;
   case  168205 : return 132;
   case  168204 : return 133;
   case  168203 : return 134;
   case  168181 : return 135;
   case  168177 : return 136;
   case  168175 : return 137;
   case  168173 : return 138;
   case  168172 : return 139;
   case  168171 : return 140;
   case  168115 : return 141;
   case  168108 : return 142;
   case  168107 : return 143;
   case  168105 : return 144;
   case  168104 : return 145;
   case  168103 : return 146;
   case  168076 : return 147;
   case  168069 : return 148;
   case  168068 : return 149;
   case  168066 : return 150;
   case  167988 : return 151;
   case  167987 : return 152;
   case  167986 : return 153;
   case  167985 : return 154;
   case  167921 : return 155;
   case  167920 : return 156;
   case  167915 : return 157;
   case  167909 : return 158;
   case  167903 : return 159;
   case  167902 : return 160;
   case  167818 : return 161;
   case  167814 : return 162;
   case  167813 : return 163;
   case  167808 : return 164;
   case  167807 : return 165;
   case  167806 : return 166;
   case  167712 : return 167;
   case  167711 : return 168;
   case  167706 : return 169;
  default : return 199;
  } 

}
//_____________________________________________________________________________
Bool_t AliAnalysisTaskGammaFlow::TestLambda(Double_t pt,Double_t l1,Double_t l2){
  
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
Bool_t AliAnalysisTaskGammaFlow::TestLambda2(Double_t pt,Double_t l1,Double_t l2){
  
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
Double_t AliAnalysisTaskGammaFlow::TestCPV(Double_t dx, Double_t dz, Double_t pt, Int_t charge){
  //Parameterization of LHC10h period
  //_true if neutral_
  
  Double_t meanX=0;
  Double_t meanZ=0.;
  Double_t sx=TMath::Min(5.4,2.59719e+02*TMath::Exp(-pt/1.02053e-01)+
              6.58365e-01*5.91917e-01*5.91917e-01/((pt-9.61306e-01)*(pt-9.61306e-01)+5.91917e-01*5.91917e-01)+1.59219);
  Double_t sz=TMath::Min(2.75,4.90341e+02*1.91456e-02*1.91456e-02/(pt*pt+1.91456e-02*1.91456e-02)+1.60) ;
  AliAODEvent *event = dynamic_cast<AliAODEvent*>(InputEvent());
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
Bool_t AliAnalysisTaskGammaFlow::TestTOF(Double_t t, Double_t e){

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
/* 
//____________________________________________________________________________
Double_t  AliAnalysisTaskGammaFlow::ApplyFlatteningV0A(Double_t phi, Double_t c){
  //LHC10h

  Double_t v2c =2.*fV0AflatC2->GetBinContent(fV0AflatC2->FindBin(c))/fHarmonics;
  Double_t v2s =2.*fV0AflatS2->GetBinContent(fV0AflatS2->FindBin(c))/fHarmonics; 
  Double_t v3c =2.*fV0AflatC3->GetBinContent(fV0AflatC3->FindBin(c))/fHarmonics/2.;
  Double_t v3s =2.*fV0AflatS3->GetBinContent(fV0AflatS3->FindBin(c))/fHarmonics/2.; 
//  Double_t v4c =fV0AflatC4->GetBinContent(fV0AflatC4->FindBin(c));
//  Double_t v4s =fV0AflatS4->GetBinContent(fV0AflatS4->FindBin(c)); 
  //There is offset in V0Csin in period4
//  if(fRunNumber>=79) //period4
//   v2s+=0.002 ;
  
  return phi+
  v2c*TMath::Sin(fHarmonics*phi)-v2s*TMath::Cos(fHarmonics*phi)+
  v3c*TMath::Sin(2.*fHarmonics*phi)-v3s*TMath::Cos(2.*fHarmonics*phi) ;

}
//____________________________________________________________________________
Double_t  AliAnalysisTaskGammaFlow::ApplyFlatteningV0C(Double_t phi, Double_t c){
  //LHC10h
  
  Double_t v2c =2.*fV0CflatC2->GetBinContent(fV0CflatC2->FindBin(c))/fHarmonics;
  Double_t v2s =2.*fV0CflatS2->GetBinContent(fV0CflatS2->FindBin(c))/fHarmonics; 
  Double_t v3c =2.*fV0CflatC3->GetBinContent(fV0CflatC3->FindBin(c))/fHarmonics/2.;
  Double_t v3s =2.*fV0CflatS3->GetBinContent(fV0CflatS3->FindBin(c))/fHarmonics/2.; 
  
  return phi+
  v2c*TMath::Sin(fHarmonics*phi)-v2s*TMath::Cos(fHarmonics*phi)+
  v3c*TMath::Sin(2.*fHarmonics*phi)-v3s*TMath::Cos(2.*fHarmonics*phi) ;
  
}
//____________________________________________________________________________
Double_t  AliAnalysisTaskGammaFlow::ApplyFlattening(Double_t phi, Double_t c){
  
  Double_t v2c =2.*fTPCflatC2->GetBinContent(fTPCflatC2->FindBin(c))/fHarmonics;
  Double_t v2s =2.*fTPCflatS2->GetBinContent(fTPCflatS2->FindBin(c))/fHarmonics; 
  Double_t v3c =2.*fTPCflatC3->GetBinContent(fTPCflatC3->FindBin(c))/fHarmonics/2.;
  Double_t v3s =2.*fTPCflatS3->GetBinContent(fTPCflatS3->FindBin(c))/fHarmonics/2.; 
   
  return phi+
  v2c*TMath::Sin(fHarmonics*phi)-v2s*TMath::Cos(fHarmonics*phi)+
  v3c*TMath::Sin(2.*fHarmonics*phi)-v3s*TMath::Cos(2.*fHarmonics*phi) ;
   
}  
*/  
//____________________________________________________________________________
void  AliAnalysisTaskGammaFlow::EvalV0ReactionPlane(AliAODEvent * event){

  AliEventplane *eventPlane = event->GetEventplane();
  if( ! eventPlane ) { AliError("Event has no event plane"); return; }
  Double_t qx = 0, qy = 0;
  //V0A
  fRPV0A = eventPlane->CalculateVZEROEventPlane(event,8, fHarmonics, qx, qy);
  fQV0A=TMath::Sqrt(qx*qx+qy*qy) ;
  //V0C
  fRPV0C = eventPlane->CalculateVZEROEventPlane(event,9, fHarmonics, qx, qy);
  fQV0C=TMath::Sqrt(qx*qx+qy*qy) ;  

  
  while(fRPV0A<0)fRPV0A+=TMath::TwoPi()/fHarmonics ;
  while(fRPV0A>TMath::TwoPi()/fHarmonics)fRPV0A-=TMath::TwoPi()/fHarmonics ;
//  fRPV0A=ApplyFlatteningV0A(fRPV0A,fCentrality) ;
  while(fRPV0A<0)fRPV0A+=TMath::TwoPi()/fHarmonics ;
  while(fRPV0A>TMath::TwoPi()/fHarmonics)fRPV0A-=TMath::TwoPi()/fHarmonics ;
  FillHistogram("phiRPV0A",fRPV0A,fCentrality) ;
  fRPV0A = fV0AFlat->MakeFlat(fRPV0A,fCentrality) ;
  
  while(fRPV0C<0)fRPV0C+=TMath::TwoPi()/fHarmonics ;
  while(fRPV0C>TMath::TwoPi()/fHarmonics)fRPV0C-=TMath::TwoPi()/fHarmonics ;
//  fRPV0C=ApplyFlatteningV0C(fRPV0C,fCentrality) ;
  while(fRPV0C<0)fRPV0C+=TMath::TwoPi()/fHarmonics ;
  while(fRPV0C>TMath::TwoPi()/fHarmonics)fRPV0C-=TMath::TwoPi()/fHarmonics ;
  FillHistogram("phiRPV0C",fRPV0C,fCentrality) ;
  fRPV0C = fV0CFlat->MakeFlat(fRPV0C,fCentrality) ;
 
  //So far no difference between fRPV0A and fRPV0AQ
  FillHistogram("phiRPV0AQ",fRPV0A,fCentrality,fQV0A) ;
  FillHistogram("phiRPV0CQ",fRPV0C,fCentrality,fQV0C) ;

  
}

//____________________________________________________________________________
Bool_t AliAnalysisTaskGammaFlow::GetTPCEventPlane(Double_t &epAngle, Double_t &qsubRes){

    AliEventplane* ep=new AliEventplane();

    float mQx=0, mQy=0;
    float mQx1=0, mQy1=0, mQx2=0, mQy2=0;
    AliAODTrack* track;
    Double_t weight;
    Int_t idtemp = -1;
    int trackcounter1=0, trackcounter2=0;

    Int_t maxID=0;

    TObjArray *tracklist=GetEventPlaneTracks(maxID);
    ep->GetQContributionXArray()->Set(maxID);
    ep->GetQContributionYArray()->Set(maxID);
    ep->GetQContributionXArraysub1()->Set(maxID);
    ep->GetQContributionYArraysub1()->Set(maxID);
    ep->GetQContributionXArraysub2()->Set(maxID);
    ep->GetQContributionYArraysub2()->Set(maxID);

    for (int i=0; i<maxID; i++){
        weight = 1;
        track = dynamic_cast<AliAODTrack*> (tracklist->At(i));
        if (track) {

            weight=GetWeight(track);
            idtemp = i ; // track->GetID();
            // TPC only tracks have negative id ((-1)*IDESD - 1) in AOD
            //if (fIsAOD && (fUseTPCOnlyTracks)) idtemp = idtemp*(-1) - 1;

            Double_t qx=weight*cos(Double_t(fHarmonics)*track->Phi());
            Double_t qy=weight*sin(Double_t(fHarmonics)*track->Phi());
            ep->GetQContributionXArray()->AddAt(qx,idtemp);
            ep->GetQContributionYArray()->AddAt(qy,idtemp);

            mQx += (qx);
            mQy += (qy);

            // This loop splits the track set into 2 random subsets
            if( trackcounter1 < int(maxID/2.) && trackcounter2 < int(maxID/2.)){
                float random = gRandom->Rndm();
                if(random < .5){
                    mQx1 += (qx);
                    mQy1 += (qy);
                    ep->GetQContributionXArraysub1()->AddAt(qx,idtemp);
                    ep->GetQContributionYArraysub1()->AddAt(qy,idtemp);
                    trackcounter1++;
                    }
                else {
                    mQx2 += (qx);
                    mQy2 += (qy);
                    ep->GetQContributionXArraysub2()->AddAt(qx,idtemp);
                    ep->GetQContributionYArraysub2()->AddAt(qy,idtemp);
                    trackcounter2++;
                }
            }
            else{
                if( trackcounter1 >= int(maxID/2.)){
                    mQx2 += (qx);
                    mQy2 += (qy);
                    ep->GetQContributionXArraysub2()->AddAt(qx,idtemp);
                    ep->GetQContributionYArraysub2()->AddAt(qy,idtemp);
                    trackcounter2++;
                }
                else {
                    mQx1 += (qx);
                    mQy1 += (qy);
                    ep->GetQContributionXArraysub1()->AddAt(qx,idtemp);
                    ep->GetQContributionYArraysub1()->AddAt(qy,idtemp);
                    trackcounter1++;
                }
            }
        }
    }

    tracklist->Clear();
    delete tracklist;
    tracklist = NULL;

    TVector2 *mQ=new TVector2();
    mQ->Set(mQx,mQy);
    epAngle=mQ->Phi()/Double_t(fHarmonics);
    
    fQTPC = mQ->Mod() ;
   
    TVector2 *qsub1=new TVector2();
    TVector2 *qsub2=new TVector2();
    qsub1->Set(mQx1,mQy1);
    qsub2->Set(mQx2,mQy2);

    ep->SetQVector(mQ);
    ep->SetEventplaneQ(epAngle);
    ep->SetQsub(qsub1,qsub2);
    qsubRes=qsub1->Phi()/Double_t(fHarmonics) - qsub2->Phi()/Double_t(fHarmonics);

    Int_t ntracks=trackcounter1+trackcounter2;

    delete ep ;
    
    if(ntracks<3)return kFALSE;// <3 -> no subevents
    return kTRUE;
}
//_________________________________________________________________________
TObjArray* AliAnalysisTaskGammaFlow::GetEventPlaneTracks(Int_t &maxID)
{
  
    AliAODEvent * event = (AliAODEvent*)InputEvent() ;
        
    TObjArray *tracklist1=new TObjArray();
    Int_t nt = event->GetNumberOfTracks();
    for (Int_t i=0; i<nt; i++) {
      AliAODTrack *aodTrack=(AliAODTrack*)event->GetTrack(i);
      if(!aodTrack) continue ;

      if(!aodTrack->IsHybridGlobalConstrainedGlobal())
        continue ;    
      if( aodTrack->Pt()<0.15 ||  aodTrack->Pt()>20.)
	continue ;
      if( TMath::Abs(aodTrack->Eta())<0.5 || TMath::Abs(aodTrack->Eta())>0.8)
	continue ;
      tracklist1->Add(new AliAODTrack(*(aodTrack))) ;
    }  
    tracklist1->SetOwner(kTRUE) ;
    maxID=tracklist1->GetEntries() ;
    if(!tracklist1)AliError("No tracklist");
    return tracklist1;
}
//_________________________________________________________________________
Double_t AliAnalysisTaskGammaFlow::GetWeight(TObject* track1)
{
    Double_t ptweight=1;
    AliVTrack* track = dynamic_cast<AliVTrack*>(track1);
    if (track) {
        if (track->Pt()<2) ptweight=track->Pt();
        else ptweight=2;
    }

    return ptweight*GetPhiWeight(track);
}
//_________________________________________________________________________
Double_t AliAnalysisTaskGammaFlow::GetPhiWeight(TObject* track1)
{
  Double_t phiweight=1;
  AliVTrack* track = dynamic_cast<AliVTrack*>(track1);

  if (fPhiDist && track) {
      Double_t nParticles = fPhiDist->Integral();
      Double_t nPhibins = fPhiDist->GetNbinsX();

      Double_t Phi = track->Phi();

      while (Phi<0) Phi += TMath::TwoPi();
      while (Phi>TMath::TwoPi()) Phi -= TMath::TwoPi();

      Double_t PhiDistValue = fPhiDist->GetBinContent(1+TMath::FloorNint((track->Phi())*nPhibins/TMath::TwoPi()));

      if (PhiDistValue > 0) phiweight = nParticles/nPhibins/PhiDistValue;
  }
  return phiweight;
}
//_________________________________________________________________________
void AliAnalysisTaskGammaFlow::EvalResolution(){
  //Estimate resolutions for each of EP estimators
  
  Double_t x= fCentrality;
  if(x<0.5)x=0.5 ;
  if(fHarmonics==2){  
    fV0Cres=3.294365e-01*(1.+x*1.507363e-01+x*x*1.272448e-03+x*x*x*4.229882e-05-x*x*x*x*6.939672e-07)/(1.-x*8.165966e-04+x*x*3.989182e-03-x*x*x*7.727254e-05+x*x*x*x*9.914099e-07);
    fV0Ares=2.905439e-01*(1.+x*1.372979e+00+x*x*6.166033e-01-x*x*x*7.905871e-03+x*x*x*x*1.507491e-06)/(1.+x*1.902084e+00+x*x*1.735136e-01-x*x*x*1.689358e-03+x*x*x*x*1.364915e-05);
    fTPCres=3.808781e-01*(1.+x*2.386826e-01+x*x*5.192646e-02+x*x*x*4.098382e-04-x*x*x*x*1.453986e-05)/(1.+x*1.259547e-01+x*x*2.285279e-02+x*x*x*1.030309e-04-x*x*x*x*2.666835e-06);
    
    //LHC01h    
//    fV0Cres=3.325009e-01*(1+x*1.939315e-01-x*x*9.652051e-04+x*x*x*1.720295e-05-x*x*x*x*3.731234e-07 )/(1+x*2.802502e-02+x*x*1.874919e-03-x*x*x*3.980602e-05+x*x*x*x*4.447131e-07);
//    fV0Ares=2.403080e-01*(1+x*2.078031e-01+x*x*3.030761e-04+x*x*x*6.242597e-05+x*x*x*x*-1.072181e-06)/(1+x*9.135497e-03+x*x*3.208661e-03-x*x*x*5.901547e-05+x*x*x*x*8.083149e-07);
//    fTPCres=3.907111e-01*(1+x*2.036863e-01+x*x*7.435000e-04+x*x*x*6.967697e-05+x*x*x*x*-1.132069e-06)/(1+x*2.751258e-02+x*x*4.000528e-03-x*x*x*6.829516e-05+x*x*x*x*7.677720e-07);
  }
  if(fHarmonics==3){
    fV0Cres=2.836819e-01*(1+x*3.017023e-01-x*x*1.083738e-02+x*x*x*9.760039e-05)/(1.+x*2.490451e-01-x*x*7.528651e-03+x*x*x*6.073070e-05);
    fV0Ares=2.002936e-01*(1+x*2.258327e-01+x*x*3.836017e-04-x*x*x*6.525109e-05)/(1.+x*1.739020e-01+x*x*1.889962e-03-x*x*x*2.786371e-05);
    fTPCres=4.357285e-01*(1+x*1.959005e-01+x*x*2.670202e-03-x*x*x*7.921636e-05)/(1.+x*1.582993e-01+x*x*2.707513e-03-x*x*x*2.773870e-05);

    //LHC10h
//    fV0Cres=TMath::Max(0.02,2.687135e-01*(1.+x*8.117353e-01+x*x*2.271778e-03-x*x*x*1.798932e-04)/(1.+x*6.721703e-01+x*x*3.829202e-03+x*x*x*1.202772e-05)) ;
//    fV0Ares=TMath::Max(0.02,2.030253e-01*(1.+x*2.118380e-01+x*x*2.894326e-02-x*x*x*4.375165e-04)/(1.+x*1.908986e-01+x*x*2.324712e-02-x*x*x*7.957708e-05));
//    fTPCres=TMath::Max(0.02,4.829372e-01*(1.+x*1.300387e-01+x*x*1.579848e-04+x*x*x*3.856068e-05-x*x*x*x*8.861878e-07)/(1.+x*1.013947e-01+x*x*1.112352e-03+x*x*x*1.766346e-05+x*x*x*x*2.207416e-07));
  }
}
//_________________________________________________________________________
void AliAnalysisTaskGammaFlow::EvalQResolution(){
  //Estimate resolutions for each of EP estimators
  
  Double_t x= fCentrality;
  if(x<0.5)x=0.5 ;
  if(fHarmonics==2){
    fTPCQres=4.857694e-01*(1+x*2.920904e-01-x*x*6.514101e-03+x*x*x*1.018369e-04-x*x*x*x*6.986208e-07)/(1+x*9.235900e-02+x*x*6.325890e-05-x*x*x*3.163401e-05+x*x*x*x*4.747334e-07);
    fV0CQres=4.387220e-01*(1+x*1.357633e-01+x*x*1.995734e-02-x*x*x*3.375979e-04+x*x*x*x*1.024155e-06)/(1+x*2.594425e-02+x*x*1.252949e-02-x*x*x*2.449652e-04+x*x*x*x*1.850385e-06);
    fV0AQres=3.242703e-01*(1+x*2.721827e-01+x*x*9.379002e-02-x*x*x*2.090242e-03+x*x*x*x*1.154136e-05)/(1+x*1.793834e-01+x*x*3.256506e-02-x*x*x*7.118371e-04+x*x*x*x*5.327807e-06);
  }
  if(fHarmonics==3){
    if(x>70)x=70 ;//bad statistics above
    fTPCQres=6.054041e-01*(1+x*2.871979e-01-x*x*8.501828e-03+x*x*x*6.042861e-05)/(1+x*2.422995e-01-x*x*6.217416e-03+x*x*x*4.229649e-05);
    fV0CQres=3.666094e-01*(1+x*6.818115e+00+x*x*2.199777e-01-x*x*x*3.688469e-03)/(1+x*6.666333e+00+x*x*1.403492e-01+x*x*x*1.229574e-03);
    fV0AQres=1.927797e-01*(1+x*1.306917e+01-x*x*1.261409e-01-x*x*x*2.002856e-04)/(1+x*9.054694e+00-x*x*1.144981e-01+x*x*x*2.580067e-03);
  }
}
//_________________________________________________________________________
void AliAnalysisTaskGammaFlow::ApplyFinalFlattening(){//apply final fine flattening
  
  Int_t ibin = fTPCfinalC2->FindBin(fCentrality) ;
  
  //TPC
  Double_t v2c =2.*fTPCfinalC2->GetBinContent(ibin)/fHarmonics;
  Double_t v2s =2.*fTPCfinalS2->GetBinContent(ibin)/fHarmonics; 
  Double_t v3c =2.*fTPCfinalC4->GetBinContent(ibin)/fHarmonics/2.;
  Double_t v3s =2.*fTPCfinalS4->GetBinContent(ibin)/fHarmonics/2.; 
  fRP = fRP +
  v2c*TMath::Sin(fHarmonics*fRP)-v2s*TMath::Cos(fHarmonics*fRP)+
  v3c*TMath::Sin(2.*fHarmonics*fRP)-v3s*TMath::Cos(2.*fHarmonics*fRP) ;
  
  while(fRP<0)fRP+=TMath::TwoPi()/fHarmonics ;
  while(fRP>TMath::TwoPi()/fHarmonics)fRP-=(TMath::TwoPi()/fHarmonics) ;

  //V0A
  v2c =2.*fV0AfinalC2->GetBinContent(ibin)/fHarmonics;
  v2s =2.*fV0AfinalS2->GetBinContent(ibin)/fHarmonics; 
  v3c =2.*fV0AfinalC4->GetBinContent(ibin)/fHarmonics/2.;
  v3s =2.*fV0AfinalS4->GetBinContent(ibin)/fHarmonics/2.; 
  fRPV0A = fRPV0A +
  v2c*TMath::Sin(fHarmonics*fRPV0A)-v2s*TMath::Cos(fHarmonics*fRPV0A)+
  v3c*TMath::Sin(2.*fHarmonics*fRPV0A)-v3s*TMath::Cos(2.*fHarmonics*fRPV0A) ;
 
  while(fRPV0A<0)fRPV0A+=TMath::TwoPi()/fHarmonics ;
  while(fRPV0A>TMath::TwoPi()/fHarmonics)fRPV0A-=(TMath::TwoPi()/fHarmonics) ;
  
  //V0C
  v2c =2.*fV0CfinalC2->GetBinContent(ibin)/fHarmonics;
  v2s =2.*fV0CfinalS2->GetBinContent(ibin)/fHarmonics; 
  v3c =2.*fV0CfinalC4->GetBinContent(ibin)/fHarmonics/2.;
  v3s =2.*fV0CfinalS4->GetBinContent(ibin)/fHarmonics/2.; 
  fRPV0C = fRPV0C +
  v2c*TMath::Sin(fHarmonics*fRPV0C)-v2s*TMath::Cos(fHarmonics*fRPV0C)+
  v3c*TMath::Sin(2.*fHarmonics*fRPV0C)-v3s*TMath::Cos(2.*fHarmonics*fRPV0C) ;
  
  while(fRPV0C<0)fRPV0C+=TMath::TwoPi()/fHarmonics ;
  while(fRPV0C>TMath::TwoPi()/fHarmonics)fRPV0C-=(TMath::TwoPi()/fHarmonics) ;
}  
//_________________________________________________________________________
void AliAnalysisTaskGammaFlow::ApplyFinalQFlattening(){//apply final fine flattening

  //Apply Q-flattening
  fRPQ    = fTPCQFlat->MakeFlat(fRP,fCentrality) ;
  fRPQV0A = fV0AQFlat->MakeFlat(fRPV0A,fCentrality) ;
  fRPQV0C = fV0CQFlat->MakeFlat(fRPV0C,fCentrality) ;
    
  while(fRPQ<0)fRPQ+=TMath::TwoPi()/fHarmonics ;
  while(fRPQ>TMath::TwoPi()/fHarmonics)fRPQ-=(TMath::TwoPi()/fHarmonics) ;
  
  while(fRPQV0A<0)fRPQV0A+=TMath::TwoPi()/fHarmonics ;
  while(fRPQV0A>TMath::TwoPi()/fHarmonics)fRPQV0A-=(TMath::TwoPi()/fHarmonics) ;
  
  while(fRPQV0C<0)fRPQV0C+=TMath::TwoPi()/fHarmonics ;
  while(fRPQV0C>TMath::TwoPi()/fHarmonics)fRPQV0C-=(TMath::TwoPi()/fHarmonics) ;
}  
//_________________________________________________________________________
Double_t AliAnalysisTaskGammaFlow::PHOSMultiplicity(){
 Double_t x=fCentrality ;
 return TMath::Max(2.,3.65177e+01-x*1.56822+x*x*3.06817e-02-x*x*x*3.14334e-04+x*x*x*x*1.27239e-06) ;
}
//_________________________________________________________________________
Double_t AliAnalysisTaskGammaFlow::CentralityWeight(Double_t c){
  //Weight to make flat centrality distribution
  Double_t weight=1. ;
  //Central
  if(c<10.)
    weight = (4.81061e+05-1.18325e+04*c+2.33302e+04*c*c-8.46768e+03*c*c*c+1.16784e+03*c*c*c*c-5.60438e+01*c*c*c*c*c)/500552. ;

  //SemiCentral
  if(c>14. && c <50.) //flat region
    weight = 1. ;
  else if(c>10. && c<=14.)
    weight = 1.+7.56198e-03*TMath::Power(c-14.,2)+1.25110e-02*TMath::Power(c-14.,3); 
    else if(c>50. && c<=56.)
       weight = 1.-2.50575e-02*TMath::Power(c-50.,2)-1.29059e-03*TMath::Power(c-50.,3); 
     
  if(weight>0.01)
    return 1./weight ;
  else
    return 0. ;
}
//_________________________________________________________________________
Bool_t AliAnalysisTaskGammaFlow::TestPHOSEvent(AliAODEvent * event){
  //Check if event is complete
  AliAODCaloCells * cells = event->GetPHOSCells() ;
  Int_t a[10]={0} ; //left
  Int_t nCells=cells->GetNumberOfCells();
  for (Int_t iCell=0; iCell<nCells; iCell++) {
    Int_t cellAbsId = cells->GetCellNumber(iCell);
    Int_t relId[4] ;
    fPHOSGeo->AbsToRelNumbering(cellAbsId,relId);
    Int_t mod  = relId[0];
    Int_t cellX= relId[2];
    if(cellX<29)
      a[2*mod]++ ;
    else
      a[2*mod+1]++ ;
  }
  Bool_t bad=kFALSE ;
  for(Int_t mod=2; mod<8; mod++){
    if(a[mod]==0){
      bad=kTRUE;
      FillHistogram("hBadMod",float(mod)) ;
    }
  }
  
  return bad ;
}

