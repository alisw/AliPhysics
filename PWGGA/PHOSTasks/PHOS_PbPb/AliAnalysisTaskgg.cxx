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
#include "TGeoGlobalMagField.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskgg.h"
#include "AliFemtoTrack.h"
#include "AliFemtoThreeVector.h"
#include "AliFemtoPair.h"
#include "AliAODTrack.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliCaloPhoton.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSEsdCluster.h"
#include "AliPHOSCalibData.h"
#include "AliAODEvent.h"
#include "AliAODCaloCells.h"
#include "AliAODVertex.h"
#include "AliLog.h"
#include "AliCentrality.h" 
#include "AliMultSelection.h" 
#include "AliEventplane.h"
#include "TProfile.h"
#include "AliOADBContainer.h"
#include "AliMagF.h"
#include "AliAODMCParticle.h"
#include "AliEPFlattener.h"
#include "AliPIDResponse.h"

// Analysis task to fill histograms with PHOS ESD clusters and cells
// Authors: Dmitri Peressounko
// Date   : 28.05.2011

ClassImp(AliAnalysisTaskgg)

//________________________________________________________________________
AliAnalysisTaskgg::AliAnalysisTaskgg(const char *name) 
: AliAnalysisTaskSE(name),
 // fStack(0x0),
  fOutputContainer(0x0),
  fEvent(0x0),
  fPIDResponse(0x0),
  fPHOSEvent(0x0),
  fCPVEvent(0x0),
  fV0AFlat(0x0),
  fV0CFlat(0x0),
  fRunNumber(0),
  fCentrality(0.),
  fCenBin(0),
  fRP(0.),
  fPHOSGeo(0x0),
  fEventCounter(0),
  fIsPbPb(kTRUE)
  
{
  // Constructor
  for(Int_t i=0;i<10;i++){
    for(Int_t j=0;j<10;j++)
      for(Int_t k=0;k<11;k++)
	fPHOSEvents[i][j][k]=0 ;
  }
  
  fNCuts=106 ;
  sprintf(fCuts[0], "All") ;   //all clusters
  sprintf(fCuts[1], "Disp") ;  //shower shape
  sprintf(fCuts[2], "NarrowDisp") ;  //shower shape
  sprintf(fCuts[3], "CPVT") ;   //neutrality using tracks
  sprintf(fCuts[4], "CPVC") ;   //neutrality using CPV
  sprintf(fCuts[5], "CPVCT") ;   //neutrality using CPV||tracks
  sprintf(fCuts[6], "BothT") ;  //shower shape && neutrality Tracks
  sprintf(fCuts[7], "BothC") ;  //shower shape && neutrality CPV
  sprintf(fCuts[8], "BothCT") ;  //shower shape && neutrality Tracks||CPV
  sprintf(fCuts[9], "Both2C") ; //shower shape && neutrality from tacks
  sprintf(fCuts[10], "Both2T") ; //shower shape && neutrality from tacks
  sprintf(fCuts[11], "Both2CT") ; //shower shape && neutrality from tacks
  sprintf(fCuts[12], "D5") ;   //all clusters, min distance = 10 cm
  sprintf(fCuts[13], "D10") ;   //all clusters, min distance = 13 cm
  sprintf(fCuts[14], "D15") ;   //all clusters, min distance = 15 cm
  sprintf(fCuts[15], "D20") ;   //all clusters, min distance = 20 cm
  sprintf(fCuts[16],"BothCTD5") ;   //all clusters, min distance = 10 cm
  sprintf(fCuts[17],"BothD10") ;   //all clusters, min distance = 13 cm
  sprintf(fCuts[18],"BothD15") ;   //all clusters, min distance = 15 cm
  sprintf(fCuts[19],"BothD20") ;   //all clusters, min distance = 20 cm
  
  //Gamma-hadron
  sprintf(fCuts[20],"AllPipl") ;
  sprintf(fCuts[21],"AllPimi") ;
  sprintf(fCuts[22],"AllKpl") ;
  sprintf(fCuts[23],"AllKmi") ;
  sprintf(fCuts[24],"AllPrpl") ;
  sprintf(fCuts[25],"AllPrmi") ;

  //Both Disp
  sprintf(fCuts[26],"DispPipl") ;
  sprintf(fCuts[27],"DispPimi") ;
  sprintf(fCuts[28],"DispKpl") ;
  sprintf(fCuts[29],"DispKmi") ;
  sprintf(fCuts[30],"DispPrpl") ;
  sprintf(fCuts[31],"DispPrmi") ;

  //BothCT+Disp
  sprintf(fCuts[32],"BothCTPipl") ;
  sprintf(fCuts[33],"BothCTPimi") ;
  sprintf(fCuts[34],"BothCTKpl") ;
  sprintf(fCuts[35],"BothCTKmi") ;
  sprintf(fCuts[36],"BothCTPrpl") ;
  sprintf(fCuts[37],"BothCTPrmi") ;
  
  //All, H+H
  sprintf(fCuts[38],"AllPimiPimi") ;
  sprintf(fCuts[39],"AllPimiPipl") ;
  sprintf(fCuts[40],"AllPimiPrmi") ;
  sprintf(fCuts[41],"AllPimiPrpl") ;
  sprintf(fCuts[42],"AllPimiKmi") ;
  sprintf(fCuts[43],"AllPimiKpl") ;
  sprintf(fCuts[44],"AllPiplPipl") ;
  sprintf(fCuts[45],"AllPiplPrmi") ;
  sprintf(fCuts[46],"AllPiplPrpl") ;
  sprintf(fCuts[47],"AllPiplKmi") ;
  sprintf(fCuts[48],"AllPiplKpl") ;
  sprintf(fCuts[49],"AllPrmiPrmi") ;
  sprintf(fCuts[50],"AllPrmiPrpl") ;
  sprintf(fCuts[51],"AllPrmiKmi") ;
  sprintf(fCuts[52],"AllPrmiKpl") ;
  sprintf(fCuts[53],"AllPrplPrpl") ;
  sprintf(fCuts[54],"AllPrplKmi") ;
  sprintf(fCuts[55],"AllPrplKpl") ;
  sprintf(fCuts[56],"AllKmiKmi") ;
  sprintf(fCuts[57],"AllKmiKpl") ;
  sprintf(fCuts[58],"AllKplKpl") ;
  
   //Disp, H+H
  sprintf(fCuts[59],"DispPimiPimi") ;
  sprintf(fCuts[60],"DispPimiPipl") ;
  sprintf(fCuts[61],"DispPimiPrmi") ;
  sprintf(fCuts[62],"DispPimiPrpl") ;
  sprintf(fCuts[63],"DispPimiKmi") ;
  sprintf(fCuts[64],"DispPimiKpl") ;
  sprintf(fCuts[65],"DispPiplPipl") ;
  sprintf(fCuts[66],"DispPiplPrmi") ;
  sprintf(fCuts[67],"DispPiplPrpl") ;
  sprintf(fCuts[68],"DispPiplKmi") ;
  sprintf(fCuts[69],"DispPiplKpl") ;
  sprintf(fCuts[70],"DispPrmiPrmi") ;
  sprintf(fCuts[71],"DispPrmiPrpl") ;
  sprintf(fCuts[72],"DispPrmiKmi") ;
  sprintf(fCuts[73],"DispPrmiKpl") ;
  sprintf(fCuts[74],"DispPrplPrpl") ;
  sprintf(fCuts[75],"DispPrplKmi") ;
  sprintf(fCuts[76],"DispPrplKpl") ;
  sprintf(fCuts[77],"DispKmiKmi") ;
  sprintf(fCuts[78],"DispKmiKpl") ;
  sprintf(fCuts[79],"DispKplKpl") ;
  
   //Time10 ns, H+H
  sprintf(fCuts[80],"T10PimiPimi") ;
  sprintf(fCuts[81],"T10PimiPipl") ;
  sprintf(fCuts[82],"T10PiplPipl") ;
  sprintf(fCuts[83],"T10PrmiPrmi") ;
  sprintf(fCuts[84],"T10PrmiPrpl") ;
  sprintf(fCuts[85],"T10PrplPrpl") ;

  sprintf(fCuts[86],"T5PimiPimi") ;
  sprintf(fCuts[87],"T5PimiPipl") ;
  sprintf(fCuts[88],"T5PiplPipl") ;
  sprintf(fCuts[89],"T5PrmiPrmi") ;
  sprintf(fCuts[90],"T5PrmiPrpl") ;
  sprintf(fCuts[91],"T5PrplPrpl") ;

  sprintf(fCuts[92],"T3PimiPimi") ;
  sprintf(fCuts[93],"T3PimiPipl") ;
  sprintf(fCuts[94],"T3PiplPipl") ;
  sprintf(fCuts[95],"T3PrmiPrmi") ;
  sprintf(fCuts[96],"T3PrmiPrpl") ;
  sprintf(fCuts[97],"T3PrplPrpl") ;

  sprintf(fCuts[98],"AllT3") ;
  sprintf(fCuts[99],"DispT3") ;
  sprintf(fCuts[100],"BothCTT3") ;

  sprintf(fCuts[101],"AllT5") ;
  sprintf(fCuts[102],"DispT5") ;
  sprintf(fCuts[103],"BothCTT5") ;

  sprintf(fCuts[104],"AllT10") ;
  sprintf(fCuts[105],"DispT10") ;
  sprintf(fCuts[106],"BothCTT10") ;
  
  
  // Output slots #0 write into a TH1 container
 DefineOutput(1,TList::Class());

}

//________________________________________________________________________
void AliAnalysisTaskgg::UserCreateOutputObjects()
{

  // Create histograms
  // Called once
  const Int_t nRuns=200 ;
  
  if(fOutputContainer != NULL){
    delete fOutputContainer;
  }
  fOutputContainer = new THashList();
  fOutputContainer->SetOwner(kTRUE);
  
  //========QA histograms=======

  //Event selection
//   fOutputContainer->Add(new TH2F("hSelEvents","Event selection", 10,0.,10.,nRuns,0.,float(nRuns))) ;
  fOutputContainer->Add(new TH1F("hTotSelEvents","Event selection", 10,0.,10.)) ;
 
  fOutputContainer->Add(new TH2F("phiRP","Event plane", 100,0.,TMath::Pi(),100,0.,100.)) ;
//   fOutputContainer->Add(new TH2F("phiRPflat","Event plane", 100,0.,TMath::Pi(),100,0.,100.)) ;
 
 
  //vertex distribution
  fOutputContainer->Add(new TH1F("hZvertex","Z vertex position", 50,-25.,25.)) ;
  
  //Centrality
  fOutputContainer->Add(new TH1F("hCentrality","Event centrality", 100,0.,100.)) ;
  fOutputContainer->Add(new TH2F("hCenPHOS","Centrality vs PHOSclusters", 100,0.,100.,100,0.,100.)) ;
  fOutputContainer->Add(new TH2F("hCenPHOSm1","Centrality vs PHOSclusters in mod1 ", 100,0.,100.,100,0.,100.)) ;
  fOutputContainer->Add(new TH2F("hCenPHOSm2","Centrality vs PHOSclusters in mod2 ", 100,0.,100.,100,0.,100.)) ;
  fOutputContainer->Add(new TH2F("hCenPHOSm3","Centrality vs PHOSclusters in mod3 ", 100,0.,100.,100,0.,100.)) ;
  fOutputContainer->Add(new TH2F("hCenPHOSm4","Centrality vs PHOSclusters in mod4 ", 100,0.,100.,100,0.,100.)) ;
  fOutputContainer->Add(new TH2F("hCenPHOSCells","Centrality vs PHOS cells", 100,0.,100.,100,0.,1000.)) ;
  fOutputContainer->Add(new TH2F("hCenTrack","Centrality vs tracks", 100,0.,100.,100,0.,15000.)) ;  
  fOutputContainer->Add(new TH2F("hCluEvsClu_All","ClusterMult vs E",20,0.,2.,50,0.,50.)) ;
  fOutputContainer->Add(new TH2F("hCluEvsClu_CPV","ClusterMult vs E",20,0.,2.,50,0.,50.)) ;
  fOutputContainer->Add(new TH2F("hCluEvsClu_Disp","ClusterMult vs E",20,0.,2.,50,0.,50.)) ;
  fOutputContainer->Add(new TH2F("hCluEvsClu_Both","ClusterMult vs E",20,0.,2.,50,0.,50.)) ;
//   fOutputContainer->Add(new TH2F("hCluEvsCluM","ClusterMult vs E",200,0.,20.,100,0.,20.)) ;
  fOutputContainer->Add(new TH2F("hCenTOF","Centrality vs PHOS TOF", 100,0.,100.,600,-6.e-6,6.e-6)) ;
  for(Int_t mod=1; mod<5; mod++){
    fOutputContainer->Add(new TH2F(Form("hJetEMod%d_th1",mod),"Cone Energy vs Centrality ", 100,0.,50.,20,0.,100.)) ;
    fOutputContainer->Add(new TH2F(Form("hJetEMod%d_th2",mod),"Cone Energy vs Centrality ", 100,0.,50.,20,0.,100.)) ;
    fOutputContainer->Add(new TH2F(Form("hJetEMod%d_th3",mod),"Cone Energy vs Centrality ", 100,0.,50.,20,0.,100.)) ;  
  }
  
    
  //PHOS QA
//   fOutputContainer->Add(new TH1I("hCellMultEvent"  ,"PHOS cell multiplicity per event"    ,2000,0,2000));
//   fOutputContainer->Add(new TH1I("hCellMultEventM1","PHOS cell multiplicity per event, M1",2000,0,2000));
//   fOutputContainer->Add(new TH1I("hCellMultEventM2","PHOS cell multiplicity per event, M2",2000,0,2000));
//   fOutputContainer->Add(new TH1I("hCellMultEventM3","PHOS cell multiplicity per event, M3",2000,0,2000));
// 
//   fOutputContainer->Add(new TH1F("hCellEnergy"  ,"Cell energy"            ,3000,0.,30.));
//   fOutputContainer->Add(new TH1F("hCellEnergyM1","Cell energy in module 1",3000,0.,30.));
//   fOutputContainer->Add(new TH1F("hCellEnergyM2","Cell energy in module 2",3000,0.,30.));
//   fOutputContainer->Add(new TH1F("hCellEnergyM3","Cell energy in module 3",3000,0.,30.));
// 
//   fOutputContainer->Add(new TH2F("hCellNXZM1","Cell (X,Z), M1" ,64,0.5,64.5, 56,0.5,56.5));
//   fOutputContainer->Add(new TH2F("hCellNXZM2","Cell (X,Z), M2" ,64,0.5,64.5, 56,0.5,56.5));
//   fOutputContainer->Add(new TH2F("hCellNXZM3","Cell (X,Z), M3" ,64,0.5,64.5, 56,0.5,56.5));
//   fOutputContainer->Add(new TH2F("hCellEXZM1","Cell E(X,Z), M1",64,0.5,64.5, 56,0.5,56.5));
//   fOutputContainer->Add(new TH2F("hCellEXZM2","Cell E(X,Z), M2",64,0.5,64.5, 56,0.5,56.5));
//   fOutputContainer->Add(new TH2F("hCellEXZM3","Cell E(X,Z), M3",64,0.5,64.5, 56,0.5,56.5));
 			
//   fOutputContainer->Add(new TH2F("hCPVr","CPV radius",100,0.,20.,100,0.,2.));
//   fOutputContainer->Add(new TH3F("hLambda","Lambdas for all clusters",150,0.,30.,150,0.,30.,200,0.,2.));
  
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
  
  //Jet QA
  Double_t ptJet[12]={2.,3.,4.,5.,7.,10.,15.,20.,30.,50.,100.,500.} ;
  fOutputContainer->Add(new TH2F("hTrackDist","Distance to track" ,100,0.,10.,11,ptJet));
  fOutputContainer->Add(new TH1F("hTrackDistTmp","tmp",11,ptJet));
  
  
  //Single photon and pi0 spectrum
  Int_t nQ=120 ;
  Double_t qMax=0.3 ;
  
  char kTbins[6][20] ;
  sprintf(kTbins[0],"Kt00-02") ;
  sprintf(kTbins[1],"Kt02-04") ;
  sprintf(kTbins[2],"Kt04-07") ;
  sprintf(kTbins[3],"Kt07-10") ;
  sprintf(kTbins[4],"Kt10-13") ;
  sprintf(kTbins[5],"Kt13-20") ;

  
  const Int_t nCenBin=4;
  for(Int_t cen=0; cen<2; cen++){  
//     for(Int_t ikT=0; ikT<6; ikT++){ 
//       fOutputContainer->Add(new TH3F(Form("hOSLCMS_%s_cen%d",kTbins[ikT],cen),"Out-Side-Long, CMS",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));
//       fOutputContainer->Add(new TH3F(Form("hMiOSLCMS_%s_cen%d",kTbins[ikT],cen),"Out-Side-Long, CMS",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));
//     }
    for(Int_t iCut=0; iCut<12; iCut++){
      fOutputContainer->Add(new TH3F(Form("hdXdZ_%s_cen%d",fCuts[iCut],cen),"dXdZ",150,-150,150,120,-120.,120.,15,0.,3.));
      fOutputContainer->Add(new TH3F(Form("hMidXdZ_%s_cen%d",fCuts[iCut],cen),"dXdZ",150,-150,150,120,-120.,120.,15,0.,3.));
      fOutputContainer->Add(new TH3F(Form("hetaphi_%s_cen%d",fCuts[iCut],cen),"Eta-phi-E correlations",50,-0.25,0.25,100,-TMath::Pi()/12.,TMath::Pi()/12.,15,0.,3.));
      fOutputContainer->Add(new TH3F(Form("hMietaphi_%s_cen%d",fCuts[iCut],cen),"Eta-phi-E correlations",50,-0.25,0.25,100,-TMath::Pi()/12.,TMath::Pi()/12.,15,0.,3.));
      
    }      
  }
  for(Int_t cen=0; cen<nCenBin; cen++){  
    for(Int_t iCut=0; iCut<fNCuts; iCut++){
//       for(Int_t ikT=0; ikT<6; ikT++){ 
//      fOutputContainer->Add(new TH3F(Form("hOSLPF_%s_%s",fCuts[iCut],kTbins[ikT]),"Out-Side-Long, Pair Frame",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));
//      fOutputContainer->Add(new TH3F(Form("hYKPPF_%s_%s",fCuts[iCut],kTbins[ikT]),"YKP, Pair Frame",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));
//      fOutputContainer->Add(new TH3F(Form("hYKPCMS_%s_%s",fCuts[iCut],kTbins[ikT]),"YKP, CMS",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));
 
          
//       fOutputContainer->Add(new TH3F(Form("hOSLCMS_%s_%s_cen%d",fCuts[iCut],kTbins[ikT],cen),"Out-Side-Long, CMS",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));

//       fOutputContainer->Add(new TH2F(Form("hetaphi2D_%s_%s",fCuts[iCut],kTbins[ikT]),"Eta-phi-E correlations",50,-0.25,0.25,100,-TMath::Pi()/6.,TMath::Pi()/6.));
//       fOutputContainer->Add(new TH2F(Form("hetaphiRP_%s_%s",fCuts[iCut],kTbins[ikT]),"Eta-phi-E correlations",10,0.,TMath::Pi(),100,-TMath::Pi()/6.,TMath::Pi()/6.));
//      fOutputContainer->Add(new TH2F(Form("hdXdZ_%s_%s",fCuts[iCut],kTbins[ikT]),"dXdZ",200,-200,200,200,-200.,200.));

//       fOutputContainer->Add(new TH3F(Form("hetaphi_%s_%s_cen%d",fCuts[iCut],kTbins[ikT],cen),"Eta-phi-E correlations",50,-0.25,0.25,100,-TMath::Pi()/6.,TMath::Pi()/6.,5,0.,2.));
      
      
//      fOutputContainer->Add(new TH3F(Form("hMiOSLPF_%s_%s",fCuts[iCut],kTbins[ikT]),"Out-Side-Long, Pair Frame",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));
//      fOutputContainer->Add(new TH3F(Form("hMiYKPPF_%s_%s",fCuts[iCut],kTbins[ikT]),"YKP, Pair Frame",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));
//      fOutputContainer->Add(new TH3F(Form("hMiYKPCMS_%s_%s",fCuts[iCut],kTbins[ikT]),"YKP, CMS",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));

//       fOutputContainer->Add(new TH3F(Form("hMiOSLCMS_%s_%s_cen%d",fCuts[iCut],kTbins[ikT],cen),"Out-Side-Long, CMS",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));
//       fOutputContainer->Add(new TH3F(Form("hMi2OSLCMS_%s_%s",fCuts[iCut],kTbins[ikT]),"Out-Side-Long, CMS",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));

      
      
//       fOutputContainer->Add(new TH2F(Form("hMietaphi2D_%s_%s",fCuts[iCut],kTbins[ikT]),"Eta-phi-E correlations",50,-0.25,0.25,100,-TMath::Pi()/6.,TMath::Pi()/6.));
//       fOutputContainer->Add(new TH3F(Form("hMietaphi_%s_%s_cen%d",fCuts[iCut],kTbins[ikT],cen),"Eta-phi-E correlations",50,-0.25,0.25,100,-TMath::Pi()/6.,TMath::Pi()/6.,5,0.,2.));
//       fOutputContainer->Add(new TH2F(Form("hMietaphiRP_%s_%s",fCuts[iCut],kTbins[ikT]),"Eta-phi-E correlations",10,0.,TMath::Pi(),100,-TMath::Pi()/6.,TMath::Pi()/6.));
//       fOutputContainer->Add(new TH3F(Form("hMi2etaphi_%s_%s",fCuts[iCut],kTbins[ikT]),"Eta-phi-E correlations",50,-0.25,0.25,100,-TMath::Pi()/6.,TMath::Pi()/6.,20,-0.2,0.2));
//      fOutputContainer->Add(new TH2F(Form("hMidXdZ_%s_%s",fCuts[iCut],kTbins[ikT]),"dXdZ",200,-200,200,200,-200.,200.));
    
//     }        

    fOutputContainer->Add(new TH2F(Form("hQinv_%s_cen%d",fCuts[iCut],cen),"Qinv distribution",200,0.,0.5,100,0.,10.));
    fOutputContainer->Add(new TH2F(Form("hMiQinv_%s_cen%d",fCuts[iCut],cen),"Qinv distribution",200,0.,0.5,100,0.,10.));
    fOutputContainer->Add(new TH2F(Form("hQinvCut_%s_cen%d",fCuts[iCut],cen),"Qinv distribution",200,0.,0.5,100,0.,10.));
    fOutputContainer->Add(new TH2F(Form("hMiQinvCut_%s_cen%d",fCuts[iCut],cen),"Qinv distribution",200,0.,0.5,100,0.,10.));
    fOutputContainer->Add(new TH2F(Form("hPhi_%s_cen%d",fCuts[iCut],cen),"Opening angle distribution",200,0.,0.3,100,0.,10.));
    fOutputContainer->Add(new TH2F(Form("hMiPhi_%s_cen%d",fCuts[iCut],cen),"Opening angle distribution",200,0.,0.3,100,0.,10.));
    fOutputContainer->Add(new TH2F(Form("hR_%s_cen%d",fCuts[iCut],cen),"Relative distance distribution",200,0.,100.,100,0.,10.));
    fOutputContainer->Add(new TH2F(Form("hMiR_%s_cen%d",fCuts[iCut],cen),"Relative distance distribution",200,0.,100.,100,0.,10.));

    if(iCut>=38 && iCut<59){//hadron combinations
      fOutputContainer->Add(new TH2F(Form("hQinvPrimH_%s_cen%d",fCuts[iCut],cen),"Qinv distribution",200,0.,0.5,100,0.,10.));
      fOutputContainer->Add(new TH2F(Form("hMiQinvPrimH_%s_cen%d",fCuts[iCut],cen),"Qinv distribution",200,0.,0.5,100,0.,10.));
      fOutputContainer->Add(new TH2F(Form("hQinvHitH_%s_cen%d",fCuts[iCut],cen),"Qinv distribution",200,0.,0.5,100,0.,10.));
      fOutputContainer->Add(new TH2F(Form("hMiQinvHitH_%s_cen%d",fCuts[iCut],cen),"Qinv distribution",200,0.,0.5,100,0.,10.));
    }
        
    
    
    }
  }

//   for(Int_t ikT=0; ikT<6; ikT++){ 
//      fOutputContainer->Add(new TH2F(Form("hSLfine_%s",kTbins[ikT]),"Out-Side",1000,-0.5,0.5,1000,-0.5,0.5));
//      fOutputContainer->Add(new TH2F(Form("hMiSLfine_%s",kTbins[ikT]),"Out-Side",1000,-0.5,0.5,1000,-0.5,0.5));
//      fOutputContainer->Add(new TH3F(Form("hSLr_%s",kTbins[ikT]),"Side-Long-r",nQ,-qMax,qMax,nQ,-qMax,qMax,30,0.,30.));
//      fOutputContainer->Add(new TH3F(Form("hMiSLr_%s",kTbins[ikT]),"Side-Long-r",nQ,-qMax,qMax,nQ,-qMax,qMax,30,0.,30.));
//   }
  
  
  
  
//   fOutputContainer->Add(new TH3F("hConvPi0","Converted pions",100,0.,10.,100,0.,10.,200,0.,TMath::Pi()));
//   fOutputContainer->Add(new TH1F("hConvPi0Angle","Angle",200,0.,TMath::Pi()));
//   fOutputContainer->Add(new TH3F("hMCConvPi0True","Converted pions",100,0.,10.,100,0.,10.,200,0.,TMath::Pi()));
//   fOutputContainer->Add(new TH3F("hMCConvPi0","Converted pions",100,0.,10.,100,0.,10.,200,0.,TMath::Pi()));
//   fOutputContainer->Add(new TH1F("hMCConvPi0Angle","Angle",200,0.,TMath::Pi()));
//   fOutputContainer->Add(new TH3F("hMCChConvPi0","Converted pions",100,0.,10.,100,0.,10.,200,0.,TMath::Pi()));
//   fOutputContainer->Add(new TH1F("hMCChConvPi0Angle","Angle",200,0.,TMath::Pi()));
//   
//   fOutputContainer->Add(new TH2F("hVtxR","Vertex dR",200,-100.,100.,500,0.,500.));
//   fOutputContainer->Add(new TH2F("hVtxRPhi","Vertex dR",200,-100.,100.,100,-TMath::Pi(),TMath::Pi()));
//   fOutputContainer->Add(new TH2F("hVtxRTheta","Vertex dR",200,-100.,100.,100,-0.5,0.5));
//   fOutputContainer->Add(new TH2F("hVtxPhi","Vertex dPhi",100,-TMath::Pi(),TMath::Pi(),100,0.,TMath::Pi()));
//   fOutputContainer->Add(new TH2F("hVtxTheta","Vertex dTheta",100,-0.5,0.5,100,-0.5,0.5));
  
		 //PHOS calibration QA
/*
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
*/    
  
  PostData(1, fOutputContainer);

}

//________________________________________________________________________
void AliAnalysisTaskgg::UserExec(Option_t *) 
{
  // Main loop, called for each event
  // Analyze ESD/AOD
    
  FillHistogram("hTotSelEvents",0.5) ;
  
  
  fEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fEvent) {
    Printf("ERROR: Could not retrieve event");
    PostData(1, fOutputContainer);
    return;
  }
  
  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  fRunNumber=fEvent->GetRunNumber() ;
//  FillHistogram("hSelEvents",1.5,fRunNumber-0.5) ;
  FillHistogram("hTotSelEvents",1.5) ;
  
 
//   if(fEventCounter == 0) {
//     //Get Event Plane flattening
//     Int_t run = fEvent->GetRunNumber() ;
//     AliOADBContainer flatContainer("phosFlat");
//     flatContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSflat.root","phosFlat");
//     TObjArray *arr = (TObjArray*)flatContainer.GetObject(run,"phosFlat");
//     if(!arr){
//       AliError(Form("Can not read Flattening for run %d. \n From file $ALICE_PHYSICS/OADB/PHOS/PHOSflat.root",run)) ;    
//       arr = (TObjArray*)flatContainer.GetObject(170593,"phosFlat"); //default
//     }
//         
//     AliInfo(Form("Setting PHOS flattening with name %s \n",arr->GetName())) ;
// //    AliEPFlattener * h = (AliEPFlattener*)arr->At(0) ;  
// //      if(fTPCFlat) delete fTPCFlat ;
// //      fTPCFlat = new AliEPFlattener() ;
// //      fTPCFlat = h ;
//     AliEPFlattener * h = (AliEPFlattener*)arr->At(1) ;  
//     if(fV0AFlat) delete fV0AFlat ;
//     fV0AFlat = new AliEPFlattener() ;
//     fV0AFlat = h ;
//     h = (AliEPFlattener*)arr->At(2) ;  
//     if(fV0CFlat) delete fV0CFlat ;
//     fV0CFlat = new AliEPFlattener() ;
//     fV0CFlat = h ;
//    
//     
//     fEMCALgeo = AliEMCALGeometry::GetInstance();
//     
//     
//     fEventCounter++ ;
//   }

  //Take Geometry from Tender
  if(!fPHOSGeo)
    fPHOSGeo = AliPHOSGeometry::GetInstance() ;
  
  // Checks if we have a primary vertex
  // Get primary vertices form AOD
  const AliAODVertex *esdVertex5 = fEvent->GetPrimaryVertex();

  Double_t vtx5[3] ={0.,0.,0.};
  
  vtx5[0] = esdVertex5->GetX();
  vtx5[1] = esdVertex5->GetY();
  vtx5[2] = esdVertex5->GetZ();
  
  
  FillHistogram("hZvertex",esdVertex5->GetZ());
  if (TMath::Abs(esdVertex5->GetZ()) > 10. ){
    PostData(1, fOutputContainer);
    return;
  }
//   FillHistogram("hSelEvents",2.5,fRunNumber-0.5) ;
  FillHistogram("hTotSelEvents",2.5) ;

  
  if(fEvent->IsPileupFromSPD()){
   PostData(1, fOutputContainer);
    return;
  } 

//   FillHistogram("hSelEvents",3.5,fRunNumber-0.5) ;
  FillHistogram("hTotSelEvents",3.5) ;  
  
  //Vtx class z-bin
  Int_t zvtx = (Int_t)((vtx5[2]+10.)/2.) ;
  if(zvtx<0)zvtx=0 ;
  if(zvtx>9)zvtx=9 ;

  fCentrality = 300; 
  if(fRunNumber<209122){ //Run1
    AliCentrality * centrality = fEvent->GetCentrality() ; 
    fCentrality=centrality->GetCentralityPercentile("V0M");
  }
  else{
    AliMultSelection * MultSelection = (AliMultSelection * ) fEvent->FindListObject("MultSelection");
    if( !MultSelection) {
       //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
       AliWarning("AliMultSelection object not found!");
    }else{
       fCentrality = MultSelection->GetMultiplicityPercentile("V0M");
    }
  }
  
//   FillHistogram("hSelEvents",3.5,fRunNumber-0.5) ;
  FillHistogram("hTotSelEvents",4.5) ;  

  FillHistogram("hCentrality",fCentrality) ;
 
  if( fCentrality <= 0. || fCentrality>80. ){
    PostData(1, fOutputContainer);
    return;
  }

//   FillHistogram("hSelEvents",4.5,fRunNumber-0.5) ;
  FillHistogram("hTotSelEvents",5.5) ;

  if(fCentrality<10.)
    fCenBin=0 ;
  else if(fCentrality<20.)
    fCenBin=1 ;
  else if(fCentrality<40.)
    fCenBin=2 ;
  else 
    fCenBin=3 ;


  //reaction plane
  Int_t irp=0 ;
  AliEventplane *eventPlane = fEvent->GetEventplane();
  if( ! eventPlane ) { //Event has no event plane
      PostData(1, fOutputContainer);
      return;
  }
  const Int_t harmonics = 2; 
  Double_t qx=0., qy=0.;  
//   Double_t rpV0A = eventPlane->CalculateVZEROEventPlane(fEvent,8, harmonics,qx,qy);
//   //V0C
//   Double_t rpV0C = eventPlane->CalculateVZEROEventPlane(fEvent,9, harmonics,qx,qy);
// 
  //Whole V0
  fRP = eventPlane->CalculateVZEROEventPlane(fEvent,10, harmonics,qx,qy);
  while(fRP<0)fRP+=TMath::TwoPi()/harmonics ;
  while(fRP>TMath::TwoPi()/harmonics)fRP-=TMath::TwoPi()/harmonics ;
    
  FillHistogram("phiRP",fRP,fCentrality) ;  
//   FillHistogram("phiRPflat",fRP,fCentrality) ;  

  //Reaction plane is defined in the range (0;pi)
  //We have 10 bins
  irp=Int_t(10.*(fRP)/TMath::Pi());
  if(irp>9)irp=9 ;
    
  
 
//   FillHistogram("hSelEvents",4.5,fRunNumber-0.5) ;
  FillHistogram("hTotSelEvents",6.5) ;
  

  if(!fPHOSEvents[zvtx][fCenBin][irp]) 
    fPHOSEvents[zvtx][fCenBin][irp]=new TList() ;
  TList * prevPHOS = fPHOSEvents[zvtx][fCenBin][irp] ;

  if(fPHOSEvent)
    fPHOSEvent->Clear() ;
  else
    fPHOSEvent = new TClonesArray("AliCaloPhoton",200) ;
    
  if(fCPVEvent)
     fCPVEvent->Clear() ;
  else
     fCPVEvent= new TClonesArray("TVector3",100) ; 

  ReclusterizeCPV() ;

  TVector3 vertex(vtx5);
  
  Int_t multClust = fEvent->GetNumberOfCaloClusters();
  Int_t inPHOS=0; 
  Int_t nPHOSclu[5]={0} ;
  
  AliAODCaloCells * cells = fEvent->GetPHOSCells() ;
  FillHistogram("hCenPHOSCells",fCentrality,cells->GetNumberOfCells()) ;
  FillHistogram("hCenTrack",fCentrality,fEvent->GetNumberOfTracks()) ;
  
  //Fill track distribution
  Int_t nt = fEvent->GetNumberOfTracks() ;
  Double_t phiPHOS=270.*TMath::DegToRad() ;
  TH1F * tmp = (TH1F*)fOutputContainer->FindObject("hTrackDistTmp") ;
  tmp->Reset() ;
  for(Int_t i=1; i<12; i++)tmp->SetBinContent(i,9.999) ;
  for (Int_t i=0; i<nt; i++) {
     AliAODTrack *aodTrack=static_cast<AliAODTrack*>(fEvent->GetTrack(i));
    if(!aodTrack->IsHybridGlobalConstrainedGlobal())
      continue ;    
    if(aodTrack->Pt()<2.)
      continue;
    Double_t dphi=aodTrack->Phi()-phiPHOS ;
    while(dphi<-TMath::Pi())dphi+=TMath::TwoPi();
    while(dphi>TMath::Pi()) dphi-=TMath::TwoPi();
    Double_t deta=aodTrack->Eta() ;
    Double_t r=TMath::Sqrt(dphi*dphi+deta*deta);
    Int_t ibin=tmp->FindBin(aodTrack->Pt()) ;
    tmp->SetBinContent(ibin,TMath::Min(tmp->GetBinContent(ibin),r)) ;
  }
  for(Int_t i=1; i<12; i++)
    FillHistogram("hTrackDist",tmp->GetBinContent(i),tmp->GetBinCenter(i)) ;
  
  
  
  TVector3 localPos ;
  for (Int_t i=0; i<multClust; i++) {
    AliAODCaloCluster *clu = fEvent->GetCaloCluster(i);
    if (clu->GetType() !=AliVCluster::kPHOSNeutral ) continue;
    if (clu->E()<0.1) continue;
    
//    if(clu->GetDistanceToBadChannel()<2.5)
//      continue ;
 
    Float_t  position[3];
    clu->GetPosition(position);
    TVector3 global(position) ;
    Int_t relId[4] ;
    fPHOSGeo->GlobalPos2RelId(global,relId) ;
    Int_t mod  = relId[0] ;
    Int_t cellX = relId[2];
    Int_t cellZ = relId[3] ;
    
    if(mod==4)
      continue ;  
    //Remove 6 noisy channels in run 139036
    if(fEvent->GetRunNumber()==139036 && mod==1 && 
       (cellX==9||cellX==10||cellX==11) && (cellZ==45 || cellZ==46))
      continue ;

    FillHistogram(Form("hTofM%d",mod),clu->E(),clu->GetTOF()) ;
    if(clu->E()>1.)
      FillHistogram("hCenTOF",fCentrality,clu->GetTOF()) ;
    if((clu->GetTOF()>25.e-9) || (clu->GetTOF() <-25.e-9) )
      continue ;
    
    
    if(clu->GetNCells()<3) continue ;
    if(clu->GetM02()<0.2)   continue ;    
    if(clu->GetMCEnergyFraction()>0.98) //Ecross fCuts, should be filled with Tender
     continue ;    
           
    TLorentzVector pv1 ;
    clu->GetMomentum(pv1 ,vtx5);
    
    FillHistogram(Form("hCluLowM%d",mod),cellX,cellZ,1.);
    if(clu->E()>1.5){
      FillHistogram(Form("hCluHighM%d",mod),cellX,cellZ,1.);
    }
 
//     FillHistogram(Form("hLambda"),clu->GetM02(),clu->GetM20(),clu->E());
//     FillHistogram("hCPVr",clu->Chi2(),clu->E());
    
 
    if(inPHOS>=fPHOSEvent->GetSize()){
      fPHOSEvent->Expand(inPHOS+50) ;
    }
    AliCaloPhoton * ph = new((*fPHOSEvent)[inPHOS]) AliCaloPhoton(pv1.X(),pv1.Py(),pv1.Z(),pv1.E()) ;
    ph->SetModule(mod) ;
//     pv1*= clu->GetCoreEnergy()/pv1.E() ;
//     ph->SetMomV2(&pv1) ;
    ph->SetNCells(clu->GetNCells());    
    ph->SetDispBit(clu->Chi2()<2.5*2.5) ;
//  Cut on FullLamdas
    ph->SetDisp2Bit(clu->Chi2()<1.5*1.5) ;
    FillHistogram("hCluEvsClu_All",clu->E(),clu->GetNCells()) ;

//    Double_t distBC=clu->GetDistanceToBadChannel();
    if(ph->IsDispOK()){
      FillHistogram(Form("hCluDispM%d",mod),cellX,cellZ,1.);
      FillHistogram("hCluEvsClu_Disp",clu->E(),clu->GetNCells()) ;    
    }
    TVector3 local ;
    fPHOSGeo->Global2Local(local,global,mod);
    Bool_t trackCPV=(clu->GetEmcCpvDistance()>2.5) ;
    Bool_t cpvCPV=TestCPV(local.X(),local.Z(),clu->E()) ;
    ph->SetCPVBit(cpvCPV) ;
    ph->SetCPV2Bit(trackCPV) ;
    if(ph->IsCPVOK()){
      FillHistogram(Form("hCluVetoM%d",mod),cellX,cellZ,1.);
      FillHistogram("hCluEvsClu_CPV",clu->E(),clu->GetNCells()) ;
      if(ph->IsDispOK()){
        FillHistogram("hCluEvsClu_Both",clu->E(),clu->GetNCells()) ;    
      }
    }
    //Mark clusters matched with tracks
    Int_t hadronBits=0 ; 
    Int_t trackId=-1;
    if(clu->GetEmcCpvDistance()<1.){
      trackId=FindTrackMatching(mod,&local) ;
      if(trackId>=0){
        AliAODTrack * track = (AliAODTrack*)fEvent->GetTrack(trackId);
        Double_t ptTrack=track->Pt() ;
        trackId=track->GetID() ;
        Bool_t electron=kTRUE ;
        Bool_t pion=kFALSE, kaon=kFALSE, proton=kFALSE ;
        Int_t charge=(track->Charge()>0) ;  
//             if(track->GetNcls(1) <2 ) electron=kFALSE ;
//             if( !(track->GetStatus() & AliESDtrack::kTPCrefit))electron=kFALSE ; 
//             if( track->GetKinkIndex(0) > 0) electron=kFALSE ;
	      
        //First rough PID
        const Double_t nSigmaBelowElectronLine=-3. ;
        const Double_t nSigmaAboveElectronLine= 5. ;
        if( fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron)<nSigmaBelowElectronLine ||
          fPIDResponse->NumberOfSigmasTPC(track,AliPID::kElectron)>nSigmaAboveElectronLine )
          electron=kFALSE ;
        const Double_t minPnSigmaAbovePionLine = 1. ;
        const Double_t maxPnSigmaAbovePionLine = 3. ;
        const Double_t nSigmaAbovePionLine = 0 ;
        if(track->P()>minPnSigmaAbovePionLine && track->P()<maxPnSigmaAbovePionLine ){
          if(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion)<nSigmaAbovePionLine){
            electron=kFALSE ;            
          }
        }
	    
      //Strict dEdx
      if(track->P()>minPnSigmaAbovePionLine && track->P()<maxPnSigmaAbovePionLine ){
        if(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion)<2.){
          electron=kFALSE ;                           
        }
      } 
      //Kaon rejection
      const Double_t minPKaonRejection=1.5 ;
      const Double_t sigmaAroundLine=1. ;
      if(track->P()<minPKaonRejection ){
        if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon))<sigmaAroundLine){
          electron=kFALSE ;                           
        }
      }
      //Proton rejection
      const Double_t minPProtonRejection=2. ;
      if(track->P()<minPProtonRejection){
        if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton))<sigmaAroundLine){
          electron=kFALSE ;                           
        }
      }
      const Double_t minPPionRejection=0.5 ;
      if(track->P()<minPPionRejection ){
        if(TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion))<sigmaAroundLine){
          electron=kFALSE ;                           
        }
      }
      //Other hadrons
      if(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kPion)<1){
        pion=kTRUE ;            
      }
      if(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kKaon)<1){
        kaon=kTRUE ;            
      }
      if(fPIDResponse->NumberOfSigmasTPC(track,AliPID::kProton)<1){
        proton=kTRUE ;            
      }
      if(electron)
          hadronBits|=1<<(0+charge);
      if(pion)
          hadronBits|=1<<(2+charge) ;  
      if(kaon)
          hadronBits|=1<<(4+charge) ;  
      if(proton)
          hadronBits|=1<<(6+charge) ;  
  
      //Set momentum of the parent track
      TLorentzVector p(track->Px(),track->Py(),track->Pz(),track->E()) ;
      ph->SetMomV2(&p);
      //Set position of Track entrance (theta, phi)
      TVector3 globaPos ;
      fPHOSGeo->Local2Global(mod, local.X()-clu->GetTrackDx(), local.Z()-clu->GetTrackDz(), globaPos) ;
      ph->SetLambdas(globaPos.Theta(),globaPos.Phi());
      }
    }   
    ph->SetTagInfo(hadronBits) ;
    ph->Pi0Id(trackId+1); 
    

    
    ph->SetTime(clu->GetTOF());
    ph->SetPrimary(clu->GetLabelAt(0)) ;
    ph->SetEMCx(position[0]) ;
    ph->SetEMCy(position[1]) ;
    ph->SetEMCz(position[2]) ;
//     ph->SetLambdas(clu->GetM20(),clu->GetM02()) ;
    ph->SetUnfolded(clu->GetNExMax()<2); // Remember, if it is unfolded          
    inPHOS++ ;
    nPHOSclu[mod]++;
  }
  

  FillHistogram("hCenPHOS",fCentrality,inPHOS) ;
  for(Int_t mod=1; mod<5; mod++)
     FillHistogram(Form("hCenPHOSm%d",mod),fCentrality,nPHOSclu[mod]) ;
  if(inPHOS==0){
    PostData(1, fOutputContainer);
    fEventCounter++;
    return ; 
  }
  
  for(Int_t mod=1; mod<5; mod++)
    fJetStatus[mod]=JetRejection(mod) ;

  const Double_t kgMass=0. ;
	

  //Real
  for (Int_t i1=0; i1<inPHOS-1; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1) ;
    
    AliFemtoTrack track1;
    AliFemtoThreeVector mom1;
    mom1.SetX(ph1->Px()) ;
    mom1.SetY(ph1->Py()) ;
    mom1.SetZ(ph1->Pz()) ;
    track1.SetP(mom1) ;
    AliFemtoParticle part1(&track1,kgMass) ;
    
    //Momentum of parent hadron if any
    const TLorentzVector * pHadron1=ph1->GetMomV2() ;
    TLorentzVector pHitHadron1(1.,1.,1.,1.) ;
    if(ph1->GetTagInfo()){
      pHitHadron1.SetE(ph1->E()) ;
      pHitHadron1.SetRho(ph1->E()) ; //assume massless photon
      pHitHadron1.SetTheta(ph1->GetLambda1()) ; //assume massless photon
      pHitHadron1.SetPhi(ph1->GetLambda2()) ; //assume massless photon
    }
    
    for (Int_t i2=i1+1; i2<inPHOS; i2++) {
      AliCaloPhoton * ph2=(AliCaloPhoton*)fPHOSEvent->At(i2) ;
      //Cut on pair
// 	if(!PairCut(ph1,ph2,kDefault))
	if(!PairCut(ph1,ph2,0)) //Distance fCuts
	  continue;
     
//      if(!SecondaryPi0Cut(ph1,ph2))
//        continue ;

      AliFemtoTrack track2;
      AliFemtoThreeVector mom2;
      mom2.SetX(ph2->Px()) ;
      mom2.SetY(ph2->Py()) ;
      mom2.SetZ(ph2->Pz()) ;
      track2.SetP(mom2) ;
      AliFemtoParticle part2(&track2,kgMass) ;
      
      //Momentum of parent hadron if any
      const TLorentzVector * pHadron2=ph2->GetMomV2() ;
      TLorentzVector pHitHadron2(1.,1.,1.,1.) ;
      if(ph2->GetTagInfo()){
        pHitHadron2.SetE(ph2->E()) ;
        pHitHadron2.SetRho(ph2->E()) ; //assume massless photon
        pHitHadron2.SetTheta(ph2->GetLambda1()) ; //assume massless photon
        pHitHadron2.SetPhi(ph2->GetLambda2()) ; //assume massless photon
      }
      
      
      //Photons are sorted, unsort them
      AliFemtoParticle *a = &part1 ;
      AliFemtoParticle *b = &part2 ;
      Double_t dEta = ph1->Eta()-ph2->Eta() ; 
      Double_t dPhi = ph1->Phi()-ph2->Phi() ; 
      Double_t dE   = ph1->E()  - ph2->E() ;
      Double_t dX = TMath::Power(ph1->EMCx() - ph2->EMCx(),2) + TMath::Power(ph1->EMCy() - ph2->EMCy(),2)  ;
      dX=TMath::Sign(TMath::Sqrt(dX),ph1->EMCx() - ph2->EMCx()) ;
      Double_t dZ = ph1->EMCz() - ph2->EMCz() ;
      Double_t dR = TMath::Sqrt(dX*dX+dZ*dZ) ;
            
      if(gRandom->Uniform()>0.5){
	    dPhi=-dPhi ;
	    dEta=-dEta ;
            a = &part2 ;
            b = &part1 ;
	    dE=-dE ;
	    dX=-dX ; 
	    dZ=-dZ ;
      }
      while(dPhi<-TMath::PiOver2())dPhi+=TMath::TwoPi() ;
      while(dPhi>TMath::PiOver2()) dPhi-=TMath::TwoPi() ;
      
      AliFemtoPair pair(a,b);
      Double_t qinv= pair.QInv();
      Double_t kT = pair.KT() ;
//       TString kTbin="" ;
//       if(kT<0.2) kTbin="Kt00-02";
//       else if(kT<0.4) kTbin="Kt02-04";
//       else if(kT<0.7) kTbin="Kt04-07";
//       else if(kT<1.) kTbin="Kt07-10";
//       else if(kT<1.3) kTbin="Kt10-13";
//       else if(kT<2.) kTbin="Kt13-20";
//       else  continue;
      
//       Double_t  qo=pair.QOutCMS();
      Double_t qs=pair.QSideCMS(), qo=pair.QOutCMS(), ql=pair.QLongCMS();
//       Double_t qspf=pair.QSidePf(),qopf=pair.QOutPf(),qlpf=pair.QLongPf() ;
//       
//       Double_t pairPhi=TMath::ATan2(ph1->Py()+ph2->Py(),ph1->Px()+ph2->Px()) ;
//       Double_t dPsi = fRP-pairPhi ;
//       while(dPsi<0)dPsi+=TMath::Pi() ;
//       while(dPsi>TMath::Pi())dPsi-=TMath::Pi() ;
      
//       // Yano-Koonin-Podgoretskii Parametrisation 
//       Double_t qP=0., qT=0., q0=0. ;
//       // source rest frame (usually lab frame)
//       pair.QYKPCMS(qP, qT, q0);

      
      Double_t qinvPrim=-1.; //Qinv of primary hadron if any
      Double_t kTPrim=-1.; 
      Double_t qinvHit=-1 ;
      if(ph1->GetTagInfo() && ph2->GetTagInfo()){
         qinvPrim=(*pHadron1 + *pHadron2).M() ;
         kTPrim=((*pHadron1 + *pHadron2).Pt())*0.5 ;
         qinvHit=(pHitHadron1 + pHitHadron2).M();
      }
      
      
      
	
      for(Int_t iCut=0; iCut<12; iCut++){
	if(!PairCut(ph1,ph2,iCut))
	    continue ;
	FillHistogram(Form("hdXdZ_%s_cen%d",fCuts[iCut],fCenBin/3),dX,dZ,kT) ;
        FillHistogram(Form("hetaphi_%s_cen%d",fCuts[iCut],fCenBin/3),dEta,dPhi,kT) ;
      }      
      
      for(Int_t iCut=0; iCut<fNCuts; iCut++){
	if(!PairCut(ph1,ph2,iCut))
	    continue ;
	
//         FillHistogram(Form("hetaphi2D_%s_%s",fCuts[iCut],kTbin.Data()),dEta,dPhi) ;
 
	if(ph1->Module()!=ph2->Module())
          continue ;
	
        
        if(iCut>=38 && iCut<59){//hadron combinations
 	  FillHistogram(Form("hQinvPrimH_%s_cen%d",fCuts[iCut],fCenBin),qinvPrim,kTPrim) ;
 	  FillHistogram(Form("hQinvHitH_%s_cen%d",fCuts[iCut],fCenBin),qinvHit,kT) ;
        }
        
        
        
	
// 	if(iCut>3){
// 	  if((jetStatus[ph1->Module()])&1<<(iCut-4))
// 	    continue ;
// 	}
	
/*		
        if(iCut==3){//Both	
          Double_t dx = ph1->EMCx()-ph2->EMCx() ;
          Double_t dz = ph1->EMCz()-ph2->EMCz() ;
	  Double_t r=TMath::Sqrt(dx*dx+dz*dz) ;
          FillHistogram(Form("hSLfine_%s",kTbin.Data()),qspf,qlpf) ;
          FillHistogram(Form("hSLr_%s",kTbin.Data()),qspf,qlpf,r) ;	  
	}*/
	  
	FillHistogram(Form("hQinv_%s_cen%d",fCuts[iCut],fCenBin),qinv,kT) ;
	if(TMath::Abs(qo) < 0.05)
	  FillHistogram(Form("hQinvCut_%s_cen%d",fCuts[iCut],fCenBin),qinv,kT) ;
        
        //Opening angle
        Double_t dPsi = ph1->Vect().Angle(ph2->Vect()) ;
	FillHistogram(Form("hPhi_%s_cen%d",fCuts[iCut],fCenBin),dPsi,kT) ;
	FillHistogram(Form("hR_%s_cen%d",fCuts[iCut],fCenBin),dR,kT) ;
        
        

        // Bertsch-Pratt momentum components in Pair Frame - written by Bekele/Humanic
//        FillHistogram(Form("hOSLPF_%s_%s",fCuts[iCut],kTbin.Data()),qspf,qopf,qlpf) ;

//         if(iCut==6&&kTbin.Length()>0){
//           // Bertsch-Pratt momentum components in Local CMS (longitudinally comoving) frame
//            FillHistogram(Form("hOSLCMS_%s_cen%d",kTbin.Data(),fCenBin),qs,qo,ql) ;
//         }    
//         FillHistogram(Form("hetaphi_%s_%s_cen%d",fCuts[iCut],kTbin.Data(),fCenBin),dEta,dPhi,kT) ;
// 	if(TMath::Abs(dEta)>0.02)
//           FillHistogram(Form("hetaphiRP_%s_%s",fCuts[iCut],kTbin.Data()),dPsi,dPhi) ;
	
//        FillHistogram(Form("hdXdZ_%s_%s",fCuts[iCut],kTbin.Data()),dX,dZ) ;


//        FillHistogram(Form("hYKPCMS_%s_%s",fCuts[iCut],kTbin.Data()),qP, qT, q0);       
      
//        FillHistogram(Form("hYKPPF_%s_%s",fCuts[iCut],kTbin.Data()),qPpf, qTpf, q0pf);       
        
      }          
    } // end of loop i2
  } // end of loop i1
  
  //now mixed
  for (Int_t i1=0; i1<inPHOS; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1) ;
    AliFemtoTrack track1;
    AliFemtoThreeVector mom1;
    mom1.SetX(ph1->Px()) ;
    mom1.SetY(ph1->Py()) ;
    mom1.SetZ(ph1->Pz()) ;
    track1.SetP(mom1) ;
    AliFemtoParticle part1(&track1,kgMass) ;
    
    //Momentum of parent hadron if any
    const TLorentzVector * pHadron1=ph1->GetMomV2() ;
    TLorentzVector pHitHadron1(1.,1.,1.,1.) ;
    if(ph1->GetTagInfo()){
      pHitHadron1.SetE(ph1->E()) ;
      pHitHadron1.SetRho(ph1->E()) ; //assume massless photon
      pHitHadron1.SetTheta(ph1->GetLambda1()) ; //assume massless photon
      pHitHadron1.SetPhi(ph1->GetLambda2()) ; //assume massless photon
    }
    
    for(Int_t ev=0; ev<prevPHOS->GetSize();ev++){
      TClonesArray * mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev)) ;
      for(Int_t i2=0; i2<mixPHOS->GetEntriesFast();i2++){
	AliCaloPhoton * ph2=(AliCaloPhoton*)mixPHOS->At(i2) ;
	
	if(!PairCut(ph1,ph2,0))
	  continue;
	
        AliFemtoTrack track2;
        AliFemtoThreeVector mom2;
        mom2.SetX(ph2->Px()) ;
        mom2.SetY(ph2->Py()) ;
        mom2.SetZ(ph2->Pz()) ;
        track2.SetP(mom2) ;
        AliFemtoParticle part2(&track2,kgMass) ;
       
        //Momentum of parent hadron if any
        const TLorentzVector * pHadron2=ph2->GetMomV2() ;
        TLorentzVector pHitHadron2(1.,1.,1.,1.) ;
        if(ph2->GetTagInfo()){
          pHitHadron2.SetE(ph2->E()) ;
          pHitHadron2.SetRho(ph2->E()) ; //assume massless photon
          pHitHadron2.SetTheta(ph2->GetLambda1()) ; //assume massless photon
          pHitHadron2.SetPhi(ph2->GetLambda2()) ; //assume massless photon
        }
          
        AliFemtoParticle *a = &part1 ;
        AliFemtoParticle *b = &part2 ;
        Double_t dEta = ph1->Eta()-ph2->Eta() ; 
        Double_t dPhi = ph1->Phi()-ph2->Phi() ; 
        Double_t dE   = ph1->E()  - ph2->E() ;
        Double_t dX = TMath::Power(ph1->EMCx() - ph2->EMCx(),2) + TMath::Power(ph1->EMCy() - ph2->EMCy(),2)  ;
        dX=TMath::Sign(TMath::Sqrt(dX),ph1->EMCx() - ph2->EMCx()) ;
        Double_t dZ = ph1->EMCz() - ph2->EMCz() ;
        Double_t dR = TMath::Sqrt(dX*dX+dZ*dZ) ;
        if(gRandom->Uniform()>0.5){
          a = &part2 ;
          b = &part1 ;
	  dEta=-dEta ;
	  dPhi=-dPhi;
	  dE=-dE ;
	  dX=-dX ; 
	  dZ=-dZ ;
        }
        while(dPhi<-TMath::PiOver2())dPhi+=TMath::TwoPi() ;
        while(dPhi>TMath::PiOver2()) dPhi-=TMath::TwoPi() ;
        AliFemtoPair pair(a,b);

	Double_t qinv= pair.QInv();
        Double_t kT = pair.KT() ;
//         TString kTbin ;
// 	Int_t ikTbin=0 ;
//         if(kT<0.2){ kTbin="Kt00-02"; ikTbin=0; }
//         else if(kT<0.4){  kTbin="Kt02-04"; ikTbin=1; }
//         else if(kT<0.7){  kTbin="Kt04-07"; ikTbin=2; }
//         else if(kT<1.0){  kTbin="Kt07-10"; ikTbin=3; }
//         else if(kT<1.3){  kTbin="Kt10-13"; ikTbin=4; }
//         else if(kT<2.0){  kTbin="Kt13-20"; ikTbin=5; }
//         else  continue;
      
//       Double_t qo=pair.QOutCMS();
      Double_t qs=pair.QSideCMS(), qo=pair.QOutCMS(), ql=pair.QLongCMS();
//       Double_t qspf=pair.QSidePf(),qopf=pair.QOutPf(),qlpf=pair.QLongPf() ;
      
//       Double_t wMix = EtaPhiWeight(ikTbin,dPhi );

//       Double_t pairPhi=TMath::ATan2(ph1->Py()+ph2->Py(),ph1->Px()+ph2->Px()) ;
//       Double_t dPsi = fRP-pairPhi ;
//       while(dPsi<0)dPsi+=TMath::Pi() ;
//       while(dPsi>TMath::Pi())dPsi-=TMath::Pi() ;
      
//       // Yano-Koonin-Podgoretskii Parametrisation 
//       Double_t qP=0., qT=0., q0=0. ;
//       // source rest frame (usually lab frame)
//       pair.QYKPCMS(qP, qT, q0);

//       Double_t qPpf=0., qTpf=0., q0pf=0. ;
//       // longitudinal comoving frame
//         pair.QYKPPF(qPpf,qTpf,q0pf) ;
	
      Double_t qinvPrim=-1.; //Qinv of primary hadron if any
      Double_t kTPrim=-1.; 
      Double_t qinvHit=-1 ;
      if(ph1->GetTagInfo() && ph2->GetTagInfo()){
         qinvPrim=(*pHadron1 + *pHadron2).M() ;
         kTPrim=((*pHadron1 + *pHadron2).Pt())*0.5 ;
         qinvHit=(pHitHadron1 + pHitHadron2).M();
      }
     
      for(Int_t iCut=0; iCut<12; iCut++){
	if(!PairCut(ph1,ph2,iCut))
	    continue ;
	FillHistogram(Form("hMidXdZ_%s_cen%d",fCuts[iCut],fCenBin/3),dX,dZ,kT) ;
        FillHistogram(Form("hMietaphi_%s_cen%d",fCuts[iCut],fCenBin/3),dEta,dPhi,kT) ;
      }      
      
      
	for(Int_t iCut=0; iCut<fNCuts; iCut++){
   	  if(!PairCut(ph1,ph2,iCut))
	    continue ;
	  
//           FillHistogram(Form("hMietaphi2D_%s_%s",fCuts[iCut],kTbin.Data()),dEta,dPhi) ;
 
	  if(ph1->Module()!=ph2->Module())
            continue ;
	  
/*	
          if(iCut==3){//Both	
            Double_t dx = ph1->EMCx()-ph2->EMCx() ;
            Double_t dz = ph1->EMCz()-ph2->EMCz() ;
	    Double_t r=TMath::Sqrt(dx*dx+dz*dz) ;
            FillHistogram(Form("hMiSLfine_%s",kTbin.Data()),qspf,qlpf) ;
            FillHistogram(Form("hMiSLr_%s",kTbin.Data()),qspf,qlpf,r) ;	  
	  }  */
	  
	  FillHistogram(Form("hMiQinv_%s_cen%d",fCuts[iCut],fCenBin),qinv,kT) ;
// 	  FillHistogram(Form("hMi2Qinv_%s",fCuts[iCut]),qinv,kT,wMix) ;
	   if(TMath::Abs(qo) < 0.05){
	     FillHistogram(Form("hMiQinvCut_%s_cen%d",fCuts[iCut],fCenBin),qinv,kT) ;
// 	     FillHistogram(Form("hMi2QinvCut_%s",fCuts[iCut]),qinv,kT,wMix) ;
	   }
          //Opening angle
          Double_t dPsi = ph1->Vect().Angle(ph2->Vect()) ;
   	  FillHistogram(Form("hMiPhi_%s_cen%d",fCuts[iCut],fCenBin),dPsi,kT) ;
   	  FillHistogram(Form("hMiR_%s_cen%d",fCuts[iCut],fCenBin),dR,kT) ;
	   
          if(iCut>=38 && iCut<59){//hadron combinations
 	    FillHistogram(Form("hMiQinvPrimH_%s_cen%d",fCuts[iCut],fCenBin),qinvPrim,kTPrim) ;
 	    FillHistogram(Form("hMiQinvHitH_%s_cen%d",fCuts[iCut],fCenBin),qinvHit,kT) ;
          }
	   
	   
          // Bertsch-Pratt momentum components in Pair Frame - written by Bekele/Humanic
//          FillHistogram(Form("hMiOSLPF_%s_%s",fCuts[iCut],kTbin.Data()),qspf,qopf,qlpf) ;
   
          // Bertsch-Pratt momentum components in Local CMS (longitudinally comoving) frame
//         if(iCut==6&&kTbin.Length()>0){          
//           FillHistogram(Form("hMiOSLCMS_%s_cen%d",kTbin.Data(),fCenBin),qs,qo,ql) ;
// //           FillHistogram(Form("hMi2OSLCMS_%s_%s",fCuts[iCut],kTbin.Data()),qs,qo,ql,wMix) ;
//         }
//           FillHistogram(Form("hMietaphi_%s_%s_cen%d",fCuts[iCut],kTbin.Data(),fCenBin),dEta,dPhi,kT) ;
// 	  if(TMath::Abs(dEta)>0.02)
//             FillHistogram(Form("hMietaphiRP_%s_%s",fCuts[iCut],kTbin.Data()),dPsi,dPhi) ;
//           FillHistogram(Form("hMi2etaphi_%s_%s",fCuts[iCut],kTbin.Data()),dEta,dPhi,dE,wMix) ;
//          FillHistogram(Form("hMidXdZ_%s_%s",fCuts[iCut],kTbin.Data()),dX,dZ) ;
	  
//          FillHistogram(Form("hMiYKPCMS_%s_%s",fCuts[iCut],kTbin.Data()),qP, qT, q0);       
      
//          FillHistogram(Form("hMiYKPPF_%s_%s",fCuts[iCut],kTbin.Data()),qPpf, qTpf, q0pf);       
	}
	
      } // end of loop i2
    }
  } // end of loop i1
  
  
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
void AliAnalysisTaskgg::Terminate(Option_t *)
{
  if(!fOutputContainer)
    return ;
  
  TFile fout("histos.root","recreate") ;  
    // Draw result to the screen
  // Called once at the end of the query
  for(Int_t i=0; i<fOutputContainer->GetSize();i++){  
     fOutputContainer->At(i)->Write() ;
      
  }
  fout.Close() ;
//   fOutputContainer->Delete() ;
  
  
}

//_____________________________________________________________________________
void AliAnalysisTaskgg::FillHistogram(const char * key,Double_t x)const{
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
void AliAnalysisTaskgg::FillHistogram(const char * key,Double_t x,Double_t y)const{
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
void AliAnalysisTaskgg::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z) const{
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
  if(tmp->IsA() == TClass::GetClass("TH3F")){
    ((TH3F*)tmp)->Fill(x,y,z) ;
    return ;
  }
}
//_____________________________________________________________________________
void AliAnalysisTaskgg::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z, Double_t w) const{
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
}

//___________________________________________________________________________
Int_t AliAnalysisTaskgg::ConvertRunNumber(Int_t run){

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

//___________________________________________________________________________
Bool_t AliAnalysisTaskgg::PairCut(const AliCaloPhoton * ph1, const AliCaloPhoton * ph2, Int_t cut) const{
  
  const Double_t kTimeCut10=10.e-9 ;
  const Double_t kTimeCut5=5.e-9 ;
  const Double_t kTimeCut3=3.e-9 ;
    
    
  if(ph1->Module()!=ph2->Module())
    return kFALSE ;   
    
  Double_t dl=999.;  
    
  Int_t bits1=ph1->GetTagInfo() ;    
  Int_t bits2=ph2->GetTagInfo() ;    
  
//   if(ph1->Module()==ph2->Module()){
//     // offset for first photon   
//     Double_t dxMax1= 3.36783/ph1->E()-11.5189/TMath::Sqrt(ph1->E())+3.08283;
//     Double_t dxMin1=-4.01590/ph1->E()+13.2841/TMath::Sqrt(ph1->E())-4.31161;
//     Double_t sigmaX1=TMath::Max(1.5,-2.07124/ph1->E()+6.69554/TMath::Sqrt(ph1->E())-1.70062);
//     Double_t sigmaZ1=1.24859122035152259;
//     if(ph1->E()>0.4) 
//       sigmaZ1=5.58984e-01*TMath::Exp(-ph1->E()*ph1->E()/2.20543/2.20543)+7.07696e-01 ;
//   
//     Double_t dxMax2= 3.36783/ph2->E()-11.5189/TMath::Sqrt(ph2->E())+3.08283;
//     Double_t dxMin2=-4.01590/ph2->E()+13.2841/TMath::Sqrt(ph2->E())-4.31161;
//     Double_t sigmaX2=TMath::Max(1.5,-2.07124/ph2->E()+6.69554/TMath::Sqrt(ph2->E())-1.70062);
//     Double_t sigmaZ2=1.24859122035152259;
//     if(ph2->E()>0.4) 
//       sigmaZ2=5.58984e-01*TMath::Exp(-ph2->E()*ph2->E()/2.20543/2.20543)+7.07696e-01 ;
// 
//     //use the largest excentricity
//     Double_t eps=TMath::Max(sigmaX1/sigmaZ1,sigmaX2/sigmaZ2) ;
//     //use this excentricity in distance calculation
//     Double_t dxC=ph1->EMCx() - ph2->EMCx() ;
//     Double_t dzC=ph1->EMCz() - ph2->EMCz() ;
//     dl=TMath::Sqrt(dxC*dxC + dzC*dzC) ;
//         //rough fCuts 
//     if(TMath::Abs(dzC)<6. && TMath::Abs(dxC)<25.){
// 
//       //distance fCuts: too close if ellipses overlap
//       Double_t dx=dxC+dxMin1+dxMin2 ;
//       Double_t dz=dzC ;
//       if(dx*dx+dz*dz*eps*eps< TMath::Power(2.5*(sigmaX1 + sigmaX2),2))
//         return kFALSE ;  
//       dx=dxC+dxMax1+dxMin2 ;
//       if(dx*dx+dz*dz*eps*eps< TMath::Power(2.5*(sigmaX1 + sigmaX2),2))
//         return kFALSE ;  
//       dx=dxC+dxMax1+dxMax2 ;
//       if(dx*dx+dz*dz*eps*eps< TMath::Power(2.5*(sigmaX1 + sigmaX2),2))
//         return kFALSE ;  
//       dx=dxC+dxMin1+dxMax2 ;
//       if(dx*dx+dz*dz*eps*eps< TMath::Power(2.5*(sigmaX1 + sigmaX2),2))
//         return kFALSE ;  
//     }
//   } 
  
 // if(fCuts==kDefault){
  if(cut==0){
    return kTRUE ;
  }
  if(cut==1){
    return ph1->IsDispOK()&& ph2->IsDispOK() ;  
  }
  if(cut==2){
    return ph1->IsDisp2OK() && ph2->IsDisp2OK()  ;  
  }
  if(cut==3){
    return ph1->IsCPVOK() && ph2->IsCPVOK()  ;  
  }
  if(cut==4){
    return ph1->IsCPV2OK() && ph2->IsCPV2OK()  ;  
  }
  if(cut==5){
    return ph1->IsCPVOK()  && ph1->IsCPV2OK() && ph2->IsCPVOK()  && ph2->IsCPV2OK() ;  
  }
  if(cut==6){
    return ph1->IsDispOK()&& ph1->IsCPVOK() && ph2->IsDispOK()&& ph2->IsCPVOK() ;  
  }
  if(cut==7){
    return ph1->IsDispOK()&& ph1->IsCPV2OK() && ph2->IsDispOK()&& ph2->IsCPV2OK() ;  
  }
  if(cut==8){
    return ph1->IsDisp2OK()&& ph1->IsCPVOK()&& ph1->IsCPV2OK() && ph2->IsDisp2OK()&& ph2->IsCPVOK()&& ph2->IsCPV2OK() ;  
  }
  if(cut==9){
    return ph1->IsDisp2OK()&& ph1->IsCPVOK() && ph2->IsDisp2OK()&& ph2->IsCPVOK() ;  
  }
  if(cut==10){
    return ph1->IsDisp2OK()&& ph1->IsCPV2OK() && ph2->IsDisp2OK()&& ph2->IsCPV2OK() ;  
  }
  if(cut==11){
    return ph1->IsDisp2OK()&& ph1->IsCPVOK()&& ph1->IsCPV2OK() && ph2->IsDisp2OK()&& ph2->IsCPVOK()&& ph2->IsCPV2OK() ;  
  }
  if(cut==12){
    return dl>5. ;
  }
  if(cut==13){
    return dl>10. ;
  }
  if(cut==14){
    return dl>15. ;
  }
  if(cut==15){
    return dl>20. ;
  }  
  if(cut==16){
    return (dl>5.) && ph1->IsDispOK() && ph1->IsCPVOK() && ph1->IsCPV2OK() && ph2->IsDispOK() && ph2->IsCPVOK() && ph2->IsCPV2OK() ;  
  }
  if(cut==17){
    return (dl>10.)&& ph1->IsDispOK() && ph1->IsCPVOK() && ph1->IsCPV2OK() && ph2->IsDispOK() && ph2->IsCPVOK() && ph2->IsCPV2OK() ;  
  }
  if(cut==18){
    return (dl>15.)&& ph1->IsDispOK() && ph1->IsCPVOK() && ph1->IsCPV2OK() && ph2->IsDispOK() && ph2->IsCPVOK() && ph2->IsCPV2OK();  
  }
  if(cut==19){
    return (dl>20.)&& ph1->IsDispOK() && ph1->IsCPVOK() && ph1->IsCPV2OK() && ph2->IsDispOK() && ph2->IsCPVOK() && ph2->IsCPV2OK() ;  
  }

  if(cut==20)
    return (bits1&1<<3)||(bits2&1<<3) ; //Pipl  
  if(cut==21)
    return (bits1&1<<2)||(bits2&1<<2) ; //Pimi  
  if(cut==22)
      return (bits1&1<<5)||(bits2&1<<5) ; //Kpl  
  if(cut==23)
      return (bits1&1<<4)||(bits2&1<<4) ; //Kmi  
  if(cut==24)
      return (bits1&1<<7)||(bits2&1<<7) ; //Prpl  
  if(cut==25)
      return (bits1&1<<6)||(bits2&1<<6) ; //Prmi  
  //Disp+hadr
  if(cut==26)
      return ((bits1&1<<3)||(bits2&1<<3)) && ph1->IsDispOK()&& ph2->IsDispOK() ; //Pipl  
  if(cut==27)
      return ((bits1&1<<2)||(bits2&1<<2)) && ph1->IsDispOK()&& ph2->IsDispOK() ; //Pimi  
  if(cut==28)
      return ((bits1&1<<5)||(bits2&1<<5)) && ph1->IsDispOK()&& ph2->IsDispOK() ; //Kpl  
  if(cut==29)
      return ((bits1&1<<4)||(bits2&1<<4)) && ph1->IsDispOK()&& ph2->IsDispOK() ; //Kmi  
  if(cut==30)
      return ((bits1&1<<7)||(bits2&1<<7)) && ph1->IsDispOK()&& ph2->IsDispOK() ; //Prpl  
  if(cut==31)
      return ((bits1&1<<6)||(bits2&1<<6)) && ph1->IsDispOK()&& ph2->IsDispOK() ; //Prmi  
 
 
   //Both+hadr
  if(cut==32)
      return ((bits1&1<<3) && ph1->IsDispOK() && ph2->IsDispOK() && ph2->IsCPVOK() && ph2->IsCPV2OK()) ||
             ((bits2&1<<3) && ph2->IsDispOK() && ph1->IsDispOK() && ph1->IsCPVOK() && ph1->IsCPV2OK())  ; //Pipl  
  if(cut==33)
      return ((bits1&1<<2) && ph1->IsDispOK() && ph2->IsDispOK() && ph2->IsCPVOK() && ph2->IsCPV2OK()) ||
             ((bits2&1<<2) && ph2->IsDispOK() && ph1->IsDispOK() && ph1->IsCPVOK() && ph1->IsCPV2OK())  ; //Pimi  
  if(cut==34)
      return ((bits1&1<<5) && ph1->IsDispOK() && ph2->IsDispOK() && ph2->IsCPVOK() && ph2->IsCPV2OK()) ||
             ((bits2&1<<5) && ph2->IsDispOK() && ph1->IsDispOK() && ph1->IsCPVOK() && ph1->IsCPV2OK())  ; //Kpl  
  if(cut==35)
      return ((bits1&1<<4) && ph1->IsDispOK() && ph2->IsDispOK() && ph2->IsCPVOK() && ph2->IsCPV2OK()) ||
             ((bits2&1<<4) && ph2->IsDispOK() && ph1->IsDispOK() && ph1->IsCPVOK() && ph1->IsCPV2OK())  ; //Kmi  
  if(cut==36)
      return ((bits1&1<<7) && ph1->IsDispOK() && ph2->IsDispOK() && ph2->IsCPVOK() && ph2->IsCPV2OK()) ||
             ((bits2&1<<7) && ph2->IsDispOK() && ph1->IsDispOK() && ph1->IsCPVOK() && ph1->IsCPV2OK())  ; //Prpl  
  if(cut==37)
      return ((bits1&1<<6) && ph1->IsDispOK() && ph2->IsDispOK() && ph2->IsCPVOK() && ph2->IsCPV2OK()) ||
             ((bits2&1<<6) && ph2->IsDispOK() && ph1->IsDispOK() && ph1->IsCPVOK() && ph1->IsCPV2OK())  ; //Prmi  
  
  
    if(cut==38)
      return (bits1&1<<2)&&(bits2&1<<2) ; //PimiPimi  
    if(cut==39)
      return ((bits1&1<<2)&&(bits2&1<<3)) ||((bits1&1<<3)&&(bits2&1<<2))  ; //Pipl Pimi  
    if(cut==40)
      return ((bits1&1<<2)&&(bits2&1<<6)) || ((bits1&1<<6)&&(bits2&1<<2)) ; //Pimi Prmi  
    if(cut==41)
      return ((bits1&1<<2)&&(bits2&1<<7)) || ((bits1&1<<7)&&(bits2&1<<2)) ; //Pimi Prpl 
    if(cut==42)
      return ((bits1&1<<2)&&(bits2&1<<4)) || ((bits1&1<<4)&&(bits2&1<<2)) ; //Pimi Kmi  
    if(cut==43)
      return ((bits1&1<<2)&&(bits2&1<<5)) || ((bits1&1<<5)&&(bits2&1<<2)) ; //Pimi Kpl  
    if(cut==44)
      return (bits1&1<<3)&&(bits2&1<<3) ; //Pipl Pipl  
    if(cut==45)
      return ((bits1&1<<3)&&(bits2&1<<6)) || ((bits1&1<<6)&&(bits2&1<<3)) ; //Pipl Prmi  
    if(cut==46)
      return ((bits1&1<<3)&&(bits2&1<<7)) || ((bits1&1<<7)&&(bits2&1<<3)) ; //Pipl Prpl  
    if(cut==47)
      return ((bits1&1<<3)&&(bits2&1<<4)) || ((bits1&1<<4)&&(bits2&1<<3)) ; //Pipl Kmi  
    if(cut==48)
      return ((bits1&1<<3)&&(bits2&1<<5)) || ((bits1&1<<5)&&(bits2&1<<3)) ; //Pipl Kmi  
   if(cut==49)
      return ((bits1&1<<6)&&(bits2&1<<6)) ; //Prmi Prmi  
    if(cut==50)
      return ((bits1&1<<6)&&(bits2&1<<7)) || ((bits1&1<<7)&&(bits2&1<<6)) ; //Prmi Prpl 
    if(cut==51)
      return ((bits1&1<<6)&&(bits2&1<<4)) || ((bits1&1<<4)&&(bits2&1<<6)) ; //Prmi Kmi 
    if(cut==52)
      return ((bits1&1<<6)&&(bits2&1<<5)) || ((bits1&1<<5)&&(bits2&1<<6)) ; //Prmi Kmi 
    if(cut==53)
      return ((bits1&1<<7)&&(bits2&1<<7)) ; //Prpl Prpl 
    if(cut==54)
      return ((bits1&1<<7)&&(bits2&1<<4)) || ((bits1&1<<4)&&(bits2&1<<7)) ; //Prpl Kmi 
    if(cut==55)
      return ((bits1&1<<7)&&(bits2&1<<5)) || ((bits1&1<<5)&&(bits2&1<<7)) ; //Prpl Kpl 
    if(cut==56)
      return ((bits1&1<<4)&&(bits2&1<<4)) ; //Kmi Kmi 
    if(cut==57)
      return ((bits1&1<<4)&&(bits2&1<<5)) || ((bits1&1<<5)&&(bits2&1<<4)) ; //Kmi Kpl 
    if(cut==58)
      return ((bits1&1<<5)&&(bits2&1<<5)) ; //Kpl Kpl 
 
 
   if(cut==59)
      return  ph1->IsDispOK() && ph2->IsDispOK() &&(bits1&1<<2)&&(bits2&1<<2) ; //PimiPimi  
    if(cut==60)
      return  ph1->IsDispOK() && ph2->IsDispOK() && ((bits1&1<<2)&&(bits2&1<<3)) ||((bits1&1<<3)&&(bits2&1<<2))  ; //Pipl Pimi  
    if(cut==61)
      return  ph1->IsDispOK() && ph2->IsDispOK() && ((bits1&1<<2)&&(bits2&1<<6)) || ((bits1&1<<6)&&(bits2&1<<2)) ; //Pimi Prmi  
    if(cut==62)
      return  ph1->IsDispOK() && ph2->IsDispOK() && ((bits1&1<<2)&&(bits2&1<<7)) || ((bits1&1<<7)&&(bits2&1<<2)) ; //Pimi Prpl 
    if(cut==63)
      return  ph1->IsDispOK() && ph2->IsDispOK() && ((bits1&1<<2)&&(bits2&1<<4)) || ((bits1&1<<4)&&(bits2&1<<2)) ; //Pimi Kmi  
    if(cut==64)
      return  ph1->IsDispOK() && ph2->IsDispOK() && ((bits1&1<<2)&&(bits2&1<<5)) || ((bits1&1<<5)&&(bits2&1<<2)) ; //Pimi Kpl  
    if(cut==65)
      return  ph1->IsDispOK() && ph2->IsDispOK() && (bits1&1<<3)&&(bits2&1<<3) ; //Pipl Pipl  
    if(cut==66)
      return  ph1->IsDispOK() && ph2->IsDispOK() && ((bits1&1<<3)&&(bits2&1<<6)) || ((bits1&1<<6)&&(bits2&1<<3)) ; //Pipl Prmi  
    if(cut==67)
      return  ph1->IsDispOK() && ph2->IsDispOK() && ((bits1&1<<3)&&(bits2&1<<7)) || ((bits1&1<<7)&&(bits2&1<<3)) ; //Pipl Prpl  
    if(cut==68)
      return  ph1->IsDispOK() && ph2->IsDispOK() && ((bits1&1<<3)&&(bits2&1<<4)) || ((bits1&1<<4)&&(bits2&1<<3)) ; //Pipl Kmi  
    if(cut==69)
      return  ph1->IsDispOK() && ph2->IsDispOK() && ((bits1&1<<3)&&(bits2&1<<5)) || ((bits1&1<<5)&&(bits2&1<<3)) ; //Pipl Kmi  
   if(cut==70)
      return  ph1->IsDispOK() && ph2->IsDispOK() && ((bits1&1<<6)&&(bits2&1<<6)) ; //Prmi Prmi  
    if(cut==71)
      return  ph1->IsDispOK() && ph2->IsDispOK() && ((bits1&1<<6)&&(bits2&1<<7)) || ((bits1&1<<7)&&(bits2&1<<6)) ; //Prmi Prpl 
    if(cut==72)
      return  ph1->IsDispOK() && ph2->IsDispOK() && ((bits1&1<<6)&&(bits2&1<<4)) || ((bits1&1<<4)&&(bits2&1<<6)) ; //Prmi Kmi 
    if(cut==73)
      return  ph1->IsDispOK() && ph2->IsDispOK() && ((bits1&1<<6)&&(bits2&1<<5)) || ((bits1&1<<5)&&(bits2&1<<6)) ; //Prmi Kmi 
    if(cut==74)
      return  ph1->IsDispOK() && ph2->IsDispOK() && ((bits1&1<<7)&&(bits2&1<<7)) ; //Prpl Prpl 
    if(cut==75)
      return  ph1->IsDispOK() && ph2->IsDispOK() && ((bits1&1<<7)&&(bits2&1<<4)) || ((bits1&1<<4)&&(bits2&1<<7)) ; //Prpl Kmi 
    if(cut==76)
      return  ph1->IsDispOK() && ph2->IsDispOK() && ((bits1&1<<7)&&(bits2&1<<5)) || ((bits1&1<<5)&&(bits2&1<<7)) ; //Prpl Kpl 
    if(cut==77)
      return  ph1->IsDispOK() && ph2->IsDispOK() && ((bits1&1<<4)&&(bits2&1<<4)) ; //Kmi Kmi 
    if(cut==78)
      return  ph1->IsDispOK() && ph2->IsDispOK() && ((bits1&1<<4)&&(bits2&1<<5)) || ((bits1&1<<5)&&(bits2&1<<4)) ; //Kmi Kpl 
    if(cut==79)
      return  ph1->IsDispOK() && ph2->IsDispOK() && ((bits1&1<<5)&&(bits2&1<<5)) ; //Kpl Kpl 

   //Time10 ns, H+H
    if(cut==80)
      return (TMath::Abs(ph1->GetTime())<kTimeCut10)&& (TMath::Abs(ph2->GetTime())<kTimeCut10) && ph1->IsDispOK() && ph2->IsDispOK() &&(bits1&1<<2)&&(bits2&1<<2) ; //PimiPimi  
    if(cut==81)
      return (TMath::Abs(ph1->GetTime())<kTimeCut10)&& (TMath::Abs(ph2->GetTime())<kTimeCut10) && ph1->IsDispOK() && ph2->IsDispOK() &&(((bits1&1<<2)&&(bits2&1<<3)) || ((bits1&1<<3)&&(bits2&1<<2)))  ; //PimiPimi  
    if(cut==82)
      return (TMath::Abs(ph1->GetTime())<kTimeCut10)&& (TMath::Abs(ph2->GetTime())<kTimeCut10) && ph1->IsDispOK() && ph2->IsDispOK() &&(bits1&1<<3)&&(bits2&1<<3) ; //PiplPipl  
    if(cut==83)
      return (TMath::Abs(ph1->GetTime())<kTimeCut10)&& (TMath::Abs(ph2->GetTime())<kTimeCut10) && ph1->IsDispOK() && ph2->IsDispOK() &&(bits1&1<<6)&&(bits2&1<<6) ; //PrmiPrmi  
    if(cut==84)
      return (TMath::Abs(ph1->GetTime())<kTimeCut10)&& (TMath::Abs(ph2->GetTime())<kTimeCut10) && ph1->IsDispOK() && ph2->IsDispOK() &&(((bits1&1<<6)&&(bits2&1<<7)) || ((bits1&1<<7)&&(bits2&1<<6)))  ; //PrmiPrpl  
    if(cut==85)
      return (TMath::Abs(ph1->GetTime())<kTimeCut10)&& (TMath::Abs(ph2->GetTime())<kTimeCut10) && ph1->IsDispOK() && ph2->IsDispOK() &&(bits1&1<<7)&&(bits2&1<<7) ; //PrplPrpl  
  
   //Time5 ns, H+H
    if(cut==86)
      return (TMath::Abs(ph1->GetTime())<kTimeCut5)&& (TMath::Abs(ph2->GetTime())<kTimeCut5) && ph1->IsDispOK() && ph2->IsDispOK() &&(bits1&1<<2)&&(bits2&1<<2) ; //PimiPimi  
    if(cut==87)
      return (TMath::Abs(ph1->GetTime())<kTimeCut5)&& (TMath::Abs(ph2->GetTime())<kTimeCut5) && ph1->IsDispOK() && ph2->IsDispOK() &&(((bits1&1<<2)&&(bits2&1<<3)) || ((bits1&1<<3)&&(bits2&1<<2)))  ; //PimiPimi  
    if(cut==88)
      return (TMath::Abs(ph1->GetTime())<kTimeCut5)&& (TMath::Abs(ph2->GetTime())<kTimeCut5) && ph1->IsDispOK() && ph2->IsDispOK() &&(bits1&1<<3)&&(bits2&1<<3) ; //PiplPipl  
    if(cut==89)
      return (TMath::Abs(ph1->GetTime())<kTimeCut5)&& (TMath::Abs(ph2->GetTime())<kTimeCut5) && ph1->IsDispOK() && ph2->IsDispOK() &&(bits1&1<<6)&&(bits2&1<<6) ; //PrmiPrmi  
    if(cut==90)
      return (TMath::Abs(ph1->GetTime())<kTimeCut5)&& (TMath::Abs(ph2->GetTime())<kTimeCut5) && ph1->IsDispOK() && ph2->IsDispOK() &&(((bits1&1<<6)&&(bits2&1<<7)) || ((bits1&1<<7)&&(bits2&1<<6)))  ; //PrmiPrpl  
    if(cut==91)
      return (TMath::Abs(ph1->GetTime())<kTimeCut5)&& (TMath::Abs(ph2->GetTime())<kTimeCut5) && ph1->IsDispOK() && ph2->IsDispOK() &&(bits1&1<<7)&&(bits2&1<<7) ; //PrplPrpl  
  
   //Time3 ns, H+H
    if(cut==92)
      return (TMath::Abs(ph1->GetTime())<kTimeCut3)&& (TMath::Abs(ph2->GetTime())<kTimeCut3) && ph1->IsDispOK() && ph2->IsDispOK() &&(bits1&1<<2)&&(bits2&1<<2) ; //PimiPimi  
    if(cut==93)
      return (TMath::Abs(ph1->GetTime())<kTimeCut3)&& (TMath::Abs(ph2->GetTime())<kTimeCut3) && ph1->IsDispOK() && ph2->IsDispOK() &&(((bits1&1<<2)&&(bits2&1<<3)) || ((bits1&1<<3)&&(bits2&1<<2)))  ; //PimiPimi  
    if(cut==94)
      return (TMath::Abs(ph1->GetTime())<kTimeCut3)&& (TMath::Abs(ph2->GetTime())<kTimeCut3) && ph1->IsDispOK() && ph2->IsDispOK() &&(bits1&1<<3)&&(bits2&1<<3) ; //PiplPipl  
    if(cut==95)
      return (TMath::Abs(ph1->GetTime())<kTimeCut3)&& (TMath::Abs(ph2->GetTime())<kTimeCut3) && ph1->IsDispOK() && ph2->IsDispOK() &&(bits1&1<<6)&&(bits2&1<<6) ; //PrmiPrmi  
    if(cut==96)
      return (TMath::Abs(ph1->GetTime())<kTimeCut3)&& (TMath::Abs(ph2->GetTime())<kTimeCut3) && ph1->IsDispOK() && ph2->IsDispOK() &&(((bits1&1<<6)&&(bits2&1<<7)) || ((bits1&1<<7)&&(bits2&1<<6)))  ; //PrmiPrpl  
    if(cut==97)
      return (TMath::Abs(ph1->GetTime())<kTimeCut3)&& (TMath::Abs(ph2->GetTime())<kTimeCut3) && ph1->IsDispOK() && ph2->IsDispOK() &&(bits1&1<<7)&&(bits2&1<<7) ; //PrplPrpl  
  
   if(cut==98)
      return (TMath::Abs(ph1->GetTime())<kTimeCut3)&& (TMath::Abs(ph2->GetTime())<kTimeCut3) ; //PrplPrpl  
   if(cut==99)
      return (TMath::Abs(ph1->GetTime())<kTimeCut3)&& (TMath::Abs(ph2->GetTime())<kTimeCut3) && ph1->IsDispOK() && ph2->IsDispOK()   ; //PrplPrpl  
   if(cut==100)
      return (TMath::Abs(ph1->GetTime())<kTimeCut3)&& (TMath::Abs(ph2->GetTime())<kTimeCut3) && ph1->IsDispOK() && ph2->IsDispOK() && ph1->IsCPVOK() && ph1->IsCPV2OK() && ph2->IsCPVOK() && ph2->IsCPV2OK()  ; //PrplPrpl  

   if(cut==101)
      return (TMath::Abs(ph1->GetTime())<kTimeCut5)&& (TMath::Abs(ph2->GetTime())<kTimeCut5) ; //PrplPrpl  
   if(cut==102)
      return (TMath::Abs(ph1->GetTime())<kTimeCut5)&& (TMath::Abs(ph2->GetTime())<kTimeCut5) && ph1->IsDispOK() && ph2->IsDispOK()   ; //PrplPrpl  
   if(cut==103)
      return (TMath::Abs(ph1->GetTime())<kTimeCut5)&& (TMath::Abs(ph2->GetTime())<kTimeCut5) && ph1->IsDispOK() && ph2->IsDispOK() && ph1->IsCPVOK() && ph1->IsCPV2OK() && ph2->IsCPVOK() && ph2->IsCPV2OK()  ; //PrplPrpl  

   if(cut==104)
      return (TMath::Abs(ph1->GetTime())<kTimeCut10)&& (TMath::Abs(ph2->GetTime())<kTimeCut10) ; //PrplPrpl  
   if(cut==105)
      return (TMath::Abs(ph1->GetTime())<kTimeCut10)&& (TMath::Abs(ph2->GetTime())<kTimeCut10) && ph1->IsDispOK() && ph2->IsDispOK()   ; //PrplPrpl  
   if(cut==106)
      return (TMath::Abs(ph1->GetTime())<kTimeCut10)&& (TMath::Abs(ph2->GetTime())<kTimeCut10) && ph1->IsDispOK() && ph2->IsDispOK() && ph1->IsCPVOK() && ph1->IsCPV2OK() && ph2->IsCPVOK() && ph2->IsCPV2OK()  ; //PrplPrpl  
      
  
  return kTRUE ;
  
}
//___________________________________________________________________________
Bool_t AliAnalysisTaskgg::PHOSCut(const AliCaloPhoton * ph1, Int_t fCuts) const{
  
//  // if(fCuts==kDefault){
//   if(fCuts==0){
//     return kTRUE ;
//   }
//   if(fCuts==1){
//     return ph1->IsDispOK()  ;  
//   }
//   if(fCuts==2){
//     return ph1->IsCPVOK()  ;  
//   }
//   if(fCuts==3){
//     return ph1->IsDispOK()  && ph1->IsCPVOK()  ;  
//   }
//   if(fCuts==4){
//     return ph1->IsDisp2OK()  ;  
//   }
//   if(fCuts==5){
//     return ph1->IsCPV2OK()  ;  
//   }
//   if(fCuts==6){
//     return ph1->IsDispOK()&& ph1->IsCPV2OK()  ;  
//   }
    
  return kTRUE ;
  
}
//___________________________________________________________________________
Int_t AliAnalysisTaskgg::JetRejection(Int_t module) const{
  //We reject events with hard tracks in the vicinity of center of PHOS module
  //track pT thresholds
  const Double_t cutPt1=10.; //Maximal Track pT
  const Double_t cutPt2=5.; //Track pT
  const Double_t cutPt3=3.; //Track pT
  
  //Cone threshold
  const Double_t cutR=0.5*0.5 ; //Cone radius squared
  
  Int_t result=0 ;
  Double_t sumE=0.;
  //azimuthal angle of PHOS module
  Double_t phiPHOS = TMath::DegToRad()*(270.+fPHOSGeo->GetPHOSAngle(module)); // (40,20,0,-20,-40) degrees
  
  Int_t nt = fEvent->GetNumberOfTracks() ;
  for (Int_t i=0; i<nt; i++) {
     AliAODTrack *aodTrack=static_cast<AliAODTrack*>(fEvent->GetTrack(i));
    if(!aodTrack->IsHybridGlobalConstrainedGlobal())
      continue ;    
    if(aodTrack->Pt()<cutPt3)
      continue;
    Double_t dphi=aodTrack->Phi()-phiPHOS ;
    while(dphi<-TMath::Pi())dphi+=TMath::TwoPi();
    while(dphi>TMath::Pi()) dphi-=TMath::TwoPi();
    Double_t deta=aodTrack->Eta() ;
    Double_t r=dphi*dphi+deta*deta;
    if(r<cutR){
      sumE+= aodTrack->Pt() ; 
      result=result|1<<0 ;
      if(aodTrack->Pt()>cutPt2)
	result=result|1<<1 ;
        if(aodTrack->Pt()>cutPt1){
	  result=result|1<<2 ;
        }
    }

  }
  if(result&1<<2)
    FillHistogram(Form("hJetEMod%d_th3",module),sumE,fCentrality) ;
  if(result&1<<1)
    FillHistogram(Form("hJetEMod%d_th2",module),sumE,fCentrality) ;
  FillHistogram(Form("hJetEMod%d_th1",module),sumE,fCentrality) ;
  return result ;
}
//___________________________________________________________________________
Double_t AliAnalysisTaskgg::EtaPhiWeight(Int_t kTbin, Double_t x) const{

  switch(kTbin){
    case 0: return 1.+0.022008*exp(-(x-0.081527)*(x-0.081527)/2./0.033761/0.033761)+0.022008*exp(-(x+0.081527)*(x+0.081527)/2./0.033761/0.033761)+0.032858*exp(-x*x/2./0.041788/0.041788) ;
    case 1: return 1.+0.016042*exp(-(x-0.085818)*(x-0.085818)/2./0.034223/0.034223)+0.016042*exp(-(x+0.085818)*(x+0.085818)/2./0.034223/0.034223)+0.032643*exp(-x*x/2./0.053440/0.053440) ;
    case 2: return 1.+0.014225*exp(-(x-0.087264)*(x-0.087264)/2./0.031328/0.031328)+0.014225*exp(-(x+0.087264)*(x+0.087264)/2./0.031328/0.031328)+0.031660*exp(-x*x/2./0.055178/0.055178) ;
    case 3: return 1.+0.018115*exp(-(x-0.081410)*(x-0.081410)/2./0.089094/0.089094)+0.018115*exp(-(x+0.081410)*(x+0.081410)/2./0.089094/0.089094)+0.016781*exp(-x*x/2./0.094293/0.094293) ;
    case 4: return 1.+0.021380*exp(-(x-0.109498)*(x-0.109498)/2./0.029483/0.029483)+0.021380*exp(-(x+0.109498)*(x+0.109498)/2./0.029483/0.029483)+0.048882*exp(-x*x/2./0.084575/0.084575) ;
    default: return 1.+0.031776*exp(-(x-0.086296)*(x-0.086296)/2./0.023534/0.023534)+0.031776*exp(-(x+0.086296)*(x+0.086296)/2./0.023534/0.023534)+0.064104*exp(-x*x/2./0.087234/0.087234) ;
  }
  
}
//_____________________________________________________________________________
void AliAnalysisTaskgg::ReclusterizeCPV(){
  //Jet-finder like algorithm:
  //find the highest cell, add m*n cells around it, 
  //next highest cell, etc. while there will be no cells above the threshold
  typedef std::pair <double, int> pairs;
 
  const double kSeedAmp=10. ; //Threshold for cluster seed  
  
  const Double_t logWeight=4.5 ;
        
  //CPV geometry
  const Int_t nCPVPadsZ=60 ;
  const Int_t nCPVPadsX=128 ;
  
  //CPV cluster parameters
  const Int_t nCluX=5 ; //cluster size  
  const Int_t nCluZ=3 ; //cluster size
  
  Double_t cpvAmp[nCPVPadsX+1][nCPVPadsZ+1]={0};   
    
  Int_t nCpvClu=0; 
  
  //Copy CPV digits to 2D matrix
  AliAODEvent *event = dynamic_cast<AliAODEvent*>(InputEvent());
  AliAODCaloCells * cells = event->GetPHOSCells() ;
  const Int_t nCells=cells->GetNumberOfCells();
  Int_t relId[4] ;
  pairs arr[nCPVPadsZ*nCPVPadsX];
  Int_t nCPV=0 ;
  for (Int_t iCell=0; iCell<nCells; iCell++) {
    Int_t cellAbsId = cells->GetCellNumber(iCell);
    if(cellAbsId>0) continue; //Select CPV cells
    cellAbsId=-cellAbsId+56*64*5;    
    
    fPHOSGeo->AbsToRelNumbering(cellAbsId,relId);
    Double_t amp=cells->GetAmplitude(iCell);
    cpvAmp[relId[2]][relId[3]]=amp ; 
    if(amp>kSeedAmp){
      arr[nCPV].first=amp;
      arr[nCPV].second=cellAbsId;
      nCPV++;
    }
  }  
  
  //Sort in increasing energy
  std::sort(arr, arr + nCPV);  
    
  //Find cell with highest amplitude and if it is above threshold, combine cluster n*m and remove cells from the matrix   
  for(Int_t iCPV=nCPV-1; iCPV>=0; iCPV--){
    fPHOSGeo->AbsToRelNumbering(arr[iCPV].second,relId);
    if(cpvAmp[relId[2]][relId[3]] >0){ 
      //Create CPV cluster with seed in this seed  
      Int_t xMin=TMath::Max(relId[2]-(nCluX-1)/2,1) ;  
      Int_t xMax=TMath::Min(relId[2]+(nCluX-1)/2,nCPVPadsX) ;  
      Int_t zMin=TMath::Max(relId[3]-(nCluZ-1)/2,1) ;  
      Int_t zMax=TMath::Min(relId[3]+(nCluZ-1)/2,nCPVPadsZ) ;  
      Double_t wtot = 0. ;
      Double_t x = 0.,z=0. ;
      Double_t totE=0.;
      for(Int_t ix=xMin; ix<=xMax; ix++) {
        for(Int_t iz=zMin; iz<=zMax; iz++) {
          totE+=cpvAmp[ix][iz];
        }
      }
      for(Int_t ix=xMin; ix<=xMax; ix++) {
        for(Int_t iz=zMin; iz<=zMax; iz++) {
            
          Float_t xi=0.,zi=0. ;
          relId[0]=3;
          relId[1]=1;
          relId[2]=ix;
          relId[3]=iz;
          fPHOSGeo->RelPosInModule(relId,xi,zi); 

          if (cpvAmp[ix][iz]>0) {
            Double_t w = TMath::Max( 0., logWeight + TMath::Log( cpvAmp[ix][iz] / totE ) ) ;
            x += xi * w ;
            z += zi * w ;
            wtot += w ;
            
            cpvAmp[ix][iz]=0; //remove cells 
          }
        }
      }
      if(wtot != 0) {
        x /= wtot ;
        z /= wtot ;
      } else {
        x = 999 ;
        z = 999 ;
      }
      //Account mis-alignment
      TVector3 globaPos ;
      fPHOSGeo->Local2Global(3, x, z, globaPos) ;

      Double_t glZ=globaPos.Z()-3.25-5.21+0.067*globaPos.Z() ;
      Double_t glX=globaPos.X()-0.6-0.48+0.067*globaPos.X();
      
      if(fCPVEvent->GetSize()<=nCpvClu)
        fCPVEvent->Expand(nCpvClu*2) ;  
      new((*fCPVEvent)[nCpvClu++]) TVector3(glX,arr[iCPV].first,glZ) ; 
    }
  }
}
//_________________________________________________________________________
Bool_t AliAnalysisTaskgg::TestCPV(Double_t emcX, Double_t emcZ, Double_t e){
   
   //Return true if neutral 
   for (Int_t j=0; j<fCPVEvent->GetEntriesFast(); j++) {
     TVector3 * cpv = (TVector3*)fCPVEvent->At(j) ; 
     if(!TestCPVCluster(cpv->X(), cpv->Z(), emcX, emcZ, e))
        return kFALSE ;
   }
   return kTRUE;  
}    
//_________________________________________________________________________
Bool_t AliAnalysisTaskgg::TestCPVCluster(Double_t cpvX, Double_t cpvZ, Double_t emcX, Double_t emcZ, Double_t e){  //return true if neutral
  
  //Return true if neutral
  Double_t dz=cpvZ - emcZ ;
  Double_t dx=cpvX - emcX ;

  
  Double_t dxMax= 3.36783/e-11.5189/TMath::Sqrt(e)+3.08283;
  Double_t dxMin=-4.01590/e+13.2841/TMath::Sqrt(e)-4.31161;
  dxMax=dx-dxMax ;
  dxMin=dx-dxMin ;
  Double_t sigmaX=TMath::Max(1.5,-2.07124/e+6.69554/TMath::Sqrt(e)-1.70062);
 
  
  Double_t sigmaZ=1.24859122035152259;
  if(e>0.4) 
    sigmaZ=5.58984e-01*TMath::Exp(-e*e/2.20543/2.20543)+7.07696e-01 ;
  

  Double_t r1=dxMax*dxMax/sigmaX/sigmaX + dz*dz/sigmaZ/sigmaZ;
  Double_t r2=dxMin*dxMin/sigmaX/sigmaX + dz*dz/sigmaZ/sigmaZ;

  return (r1>2.5*2.5) && (r2>2.5*2.5) ;
}
//___________________________________________________________________________________________________
Int_t AliAnalysisTaskgg::FindTrackMatching(Int_t mod,TVector3 *locpos){

    
  Double_t  magF =fEvent->GetMagneticField();
 
  Double_t magSign = 1.0;
  if(magF<0)magSign = -1.0;
  
  if (!TGeoGlobalMagField::Instance()->GetField()) {
     AliError("Margnetic filed was not initialized, use default") ;
     AliMagF* field = new AliMagF("Maps","Maps", magSign, magSign, AliMagF::k5kG);
     TGeoGlobalMagField::Instance()->SetField(field);
  }

  // *** Start the matching
  Int_t nt = fEvent->GetNumberOfTracks();
      
  //Calculate actual distance to PHOS module
  TVector3 globaPos ;
  fPHOSGeo->Local2Global(mod, 0.,0., globaPos) ;
  const Double_t rPHOS = globaPos.Pt() ; //Distance to center of  PHOS module
  const Double_t kYmax = 72.+10. ; //Size of the module (with some reserve) in phi direction
  const Double_t kZmax = 64.+10. ; //Size of the module (with some reserve) in z direction
  const Double_t kAlpha0=330./180.*TMath::Pi() ; //First PHOS module angular direction
  const Double_t kAlpha= 20./180.*TMath::Pi() ; //PHOS module angular size
  Double_t minDistance = 1.e6;


  Double_t gposTrack[3] ; 

  Double_t bz = ((AliMagF*)TGeoGlobalMagField::Instance()->GetField())->SolenoidField();
  bz = TMath::Sign(0.5*kAlmost0Field,bz) + bz;

  Double_t b[3]; 
  Int_t itr=-1 ;
  AliAODTrack *aodTrack=0x0 ;
  Double_t xyz[3] = {0}, pxpypz[3] = {0}, cv[21] = {0};
  for (Int_t i=0; i<nt; i++) {
      aodTrack=(AliAODTrack*)fEvent->GetTrack(i);

      // Skip the tracks having "wrong" status (has to be checked/tuned)
      if(!aodTrack->IsHybridGlobalConstrainedGlobal())
        continue ;  
      
      //Continue extrapolation from TPC outer surface
      AliExternalTrackParam outerParam;
      aodTrack->GetPxPyPz(pxpypz);
      aodTrack->GetXYZ(xyz);
      aodTrack->GetCovarianceXYZPxPyPz(cv);
      outerParam.Set(xyz,pxpypz,cv,aodTrack->Charge());
      
      Double_t z; 
      if(!outerParam.GetZAt(rPHOS,bz,z))
        continue ;

      if (TMath::Abs(z) > kZmax) 
        continue; // Some tracks miss the PHOS in Z

     
      //Direction to the current PHOS module
      Double_t phiMod=kAlpha0-kAlpha*mod ;
      if(!outerParam.RotateParamOnly(phiMod)) continue ; //RS use faster rotation if errors are not needed 
    
      Double_t y;                       // Some tracks do not reach the PHOS
      if (!outerParam.GetYAt(rPHOS,bz,y)) continue; //    because of the bending
      
      if(TMath::Abs(y) < kYmax){
        outerParam.GetBxByBz(b) ;
        outerParam.PropagateToBxByBz(rPHOS,b);        // Propagate to the matching module
        //outerParam.CorrectForMaterial(...); // Correct for the TOF material, if needed
        outerParam.GetXYZ(gposTrack) ;
        TVector3 globalPositionTr(gposTrack) ;
        TVector3 localPositionTr ;
        fPHOSGeo->Global2Local(localPositionTr,globalPositionTr,mod) ;
        Double_t ddx = locpos->X()-localPositionTr.X();
        Double_t ddz = locpos->Z()-localPositionTr.Z();
        Double_t d2 = ddx*ddx + ddz*ddz;
        if(d2 < minDistance) {
	  minDistance=d2 ;
	  itr=i ;
         }
      }
    }//Scanned all tracks
 
   return itr ;
}

