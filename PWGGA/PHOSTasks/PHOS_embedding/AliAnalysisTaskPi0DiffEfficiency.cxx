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
#include "THashList.h"

#include "AliAODMCParticle.h"
#include "AliAnalysisManager.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliStack.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskPi0DiffEfficiency.h"
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

// Analysis task to calculate pi0 registration efficiency
// as a difference between spectra before and after embedding
// Authors: Dmitri Peressounko
// Date   : Aug.2011

Double_t Scale(Double_t x){

//  return 1./1.008/1.015*(1.+0.017/(1.+x*x/2./2.)+0.03/(1.+x*x/0.6/0.6)) ;
  return 1./1.015/1.015*(1.+0.017/(1.+x*x/2./2.)+0.03/(1.+x*x/0.6/0.6)) ;


}

ClassImp(AliAnalysisTaskPi0DiffEfficiency)

//________________________________________________________________________
AliAnalysisTaskPi0DiffEfficiency::AliAnalysisTaskPi0DiffEfficiency(const char *name) 
: AliAnalysisTaskSE(name),
  fStack(0),
  fOutputContainer(0),
  fPHOSEvent(0),
  fPHOSEvent1(0),
  fPHOSEvent2(0),
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
  DefineOutput(1,TList::Class());

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
void AliAnalysisTaskPi0DiffEfficiency::UserCreateOutputObjects()
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
  Int_t nPt      = 250;
  Double_t ptMin = 0;
  Double_t ptMax = 25;
  
  Int_t nPtf      = 23;
  Double_t xPt[24]={0.6,1.,1.5,2.,2.5,3.,3.5,4.,4.5,5.,5.5,6.,7.,8.,9.,10.,12.,14.,16.,18.,20.,22.,24.,25.} ;
  Int_t nPhi=10 ;
  Double_t xPhi[11] ;
  for(Int_t i=0;i<=10;i++)xPhi[i]=i*0.1*TMath::Pi() ;
  Int_t nMm=150 ;
  Double_t xM[201] ;
  for(Int_t i=0;i<=200;i++){xM[i]=0.0025*i;}
    

  char key[55] ;
  for(Int_t cent=0; cent<6; cent++){
    //Single photon

    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad2_cen%d",cent),"All clusters",nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad4_cen%d",cent),"All clusters",nPt,ptMin,ptMax));
    fOutputContainer->Add(new TH1F(Form("hPhotAll_DistBad6_cen%d",cent),"All clusters",nPt,ptMin,ptMax));
    snprintf(key,55,"hPhotAll_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hPhotAllcore_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hPhotAllwou_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hPhotCPV_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    snprintf(key,55,"hPhotCPVcore_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    snprintf(key,55,"hPhotCPV2_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hPhotCPV2core_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hPhotDisp_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hPhotDispwou_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hPhotDisp2_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hPhotDispcore_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    snprintf(key,55,"hPhotDisp2core_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hPhotBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hPhotBothcore_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hPhotBoth2_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hPhotBoth2core_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;

    snprintf(key,55,"hNegPhotAll_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegPhotAllcore_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegPhotCPV_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    snprintf(key,55,"hNegPhotCPVcore_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    snprintf(key,55,"hNegPhotCPV2_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegPhotDisp_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegPhotDisp2_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegPhotDispcore_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    snprintf(key,55,"hNegPhotDisp2core_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegPhotBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegPhotBothcore_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegPhotBoth2_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegPhotBoth2core_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    
    snprintf(key,55,"hOldMassPtAll_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hOldMassPtAllcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hOldMassPtCPV_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hOldMassPtCPVcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hOldMassPtCPV2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hOldMassPtDisp_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hOldMassPtDisp2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hOldMassPtDispcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hOldMassPtDisp2core_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hOldMassPtBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hOldMassPtBothcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hOldMassPtBoth2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hOldMassPtBoth2core_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    
    snprintf(key,55,"hNewMassPtAll_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNewMassPtAllcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNewMassPtCPV_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNewMassPtCPVcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hNewMassPtCPV2_cen%d",cent) ;
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNewMassPtDisp_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNewMassPtDisp2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNewMassPtDispcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNewMassPtDisp2core_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNewMassPtBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNewMassPtBothcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNewMassPtBoth2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNewMassPtBoth2core_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;

    snprintf(key,55,"hNegMassPtAll_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegMassPtAllcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegMassPtCPV_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegMassPtCPVcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegMassPtCPV2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegMassPtDisp_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegMassPtDisp2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegMassPtDispcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegMassPtDisp2core_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegMassPtBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegMassPtBothcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegMassPtBoth2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegMassPtBoth2core_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    
    
    snprintf(key,55,"hMassPtAll_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMassPtAllwou_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMassPtAllcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMassPtCPV_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMassPtCPVcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMassPtCPV2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMassPtDisp_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMassPtDispwou_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMassPtDisp2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMassPtDispcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMassPtDisp2core_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMassPtBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMassPtBothcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMassPtBoth2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMassPtBoth2core_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
        
    
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
    snprintf(key,55,"hMiMassPtDispcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtDisp2core_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtBothcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtBoth2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtBoth2core_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
        
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
    snprintf(key,55,"hMCMassPtCPV2core_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtDisp_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtDispwou_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtDispcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtDisp2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtDisp2core_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtBoth2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtBothcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtBoth2core_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));

    //Single photon
    snprintf(key,55,"hMCPhotAll_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMCPhotAllwou_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMCPhotAllcore_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMCPhotCPV_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMCPhotCPVcore_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMCPhotCPV2_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMCPhotDisp_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMCPhotDispwou_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMCPhotDisp2_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMCPhotDispcore_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMCPhotDisp2core_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMCPhotBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;    
    snprintf(key,55,"hMCPhotBothcore_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;    
    snprintf(key,55,"hMCPhotBoth2_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;    
    snprintf(key,55,"hMCPhotBoth2core_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;    
    
    char phiTitle[15]={"TPC"};
    TH2F * tmp = 0;
    tmp=new TH2F(Form("hPhotPhi%sAll_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}" ,nPtf,xPt,nPhi,xPhi) ;
    tmp->Sumw2() ; fOutputContainer->Add(tmp) ;
    tmp=new TH2F(Form("hPhotPhi%sAllcore_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}" ,nPtf,xPt,nPhi,xPhi) ;
    tmp->Sumw2() ; fOutputContainer->Add(tmp) ;
    tmp=new TH2F(Form("hPhotPhi%sDisp_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}" ,nPtf,xPt,nPhi,xPhi) ;
    tmp->Sumw2() ; fOutputContainer->Add(tmp) ;
    tmp=new TH2F(Form("hPhotPhi%sDispcore_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}" ,nPtf,xPt,nPhi,xPhi) ;
    tmp->Sumw2() ; fOutputContainer->Add(tmp) ;
    tmp=new TH2F(Form("hPhotPhi%sDisp2_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}" ,nPtf,xPt,nPhi,xPhi) ;
    tmp->Sumw2() ; fOutputContainer->Add(tmp) ;
    tmp=new TH2F(Form("hPhotPhi%sDisp2core_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}" ,nPtf,xPt,nPhi,xPhi) ;
    tmp->Sumw2() ; fOutputContainer->Add(tmp) ;
    tmp=new TH2F(Form("hPhotPhi%sCPV_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}" ,nPtf,xPt,nPhi,xPhi) ;
    tmp->Sumw2() ; fOutputContainer->Add(tmp) ;
    tmp=new TH2F(Form("hPhotPhi%sCPVcore_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}" ,nPtf,xPt,nPhi,xPhi) ;
    tmp->Sumw2() ; fOutputContainer->Add(tmp) ;
    tmp=new TH2F(Form("hPhotPhi%sCPV2_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}" ,nPtf,xPt,nPhi,xPhi) ;
    tmp->Sumw2() ; fOutputContainer->Add(tmp) ;
    tmp=new TH2F(Form("hPhotPhi%sCPV2core_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}" ,nPtf,xPt,nPhi,xPhi) ;
    tmp->Sumw2() ; fOutputContainer->Add(tmp) ;
    tmp=new TH2F(Form("hPhotPhi%sBoth_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}" ,nPtf,xPt,nPhi,xPhi) ;
    tmp->Sumw2() ; fOutputContainer->Add(tmp) ;
    tmp=new TH2F(Form("hPhotPhi%sBothcore_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}" ,nPtf,xPt,nPhi,xPhi) ;
    tmp->Sumw2() ; fOutputContainer->Add(tmp) ;
    tmp=new TH2F(Form("hPhotPhi%sBoth2_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}" ,nPtf,xPt,nPhi,xPhi) ;
    tmp->Sumw2() ; fOutputContainer->Add(tmp) ;
    tmp=new TH2F(Form("hPhotPhi%sBoth2core_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}" ,nPtf,xPt,nPhi,xPhi) ;
    tmp->Sumw2() ; fOutputContainer->Add(tmp) ;

    
    //Pions for flow - with weight 1/Nclu
    TH3F * tmp3 = 0;
    tmp3 = new TH3F(Form("hMassPt%sAll_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;
    tmp3 = new TH3F(Form("hMassPt%sAllcore_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;
    tmp3 = new TH3F(Form("hMassPt%sCPV_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;
    tmp3 = new TH3F(Form("hMassPt%sCPVcore_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;
    tmp3 = new TH3F(Form("hMassPt%sCPV2_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;
    tmp3 = new TH3F(Form("hMassPt%sCPV2core_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;
    tmp3 = new TH3F(Form("hMassPt%sDisp_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;
    tmp3 = new TH3F(Form("hMassPt%sDispcore_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;
    tmp3 = new TH3F(Form("hMassPt%sDisp2_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;
    tmp3 = new TH3F(Form("hMassPt%sDisp2core_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;
    tmp3 = new TH3F(Form("hMassPt%sBoth_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;
    tmp3 = new TH3F(Form("hMassPt%sBothcore_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;    
    tmp3 = new TH3F(Form("hMassPt%sBoth2_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;
    tmp3 = new TH3F(Form("hMassPt%sBoth2core_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;

//mixed
    tmp3 = new TH3F(Form("hMiMassPt%sAll_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;
    tmp3 = new TH3F(Form("hMiMassPt%sAllcore_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;
    tmp3 = new TH3F(Form("hMiMassPt%sCPV_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;
    tmp3 = new TH3F(Form("hMiMassPt%sCPVcore_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;
    tmp3 = new TH3F(Form("hMiMassPt%sCPV2_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;
    tmp3 = new TH3F(Form("hMiMassPt%sCPV2core_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;
    tmp3 = new TH3F(Form("hMiMassPt%sDisp_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;
    tmp3 = new TH3F(Form("hMiMassPt%sDispcore_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;
    tmp3 = new TH3F(Form("hMiMassPt%sDisp2_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;
    tmp3 = new TH3F(Form("hMiMassPt%sDisp2core_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;
    tmp3 = new TH3F(Form("hMiMassPt%sBoth_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;
    tmp3 = new TH3F(Form("hMiMassPt%sBothcore_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;
    tmp3 = new TH3F(Form("hMiMassPt%sBoth2_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;
    tmp3 = new TH3F(Form("hMiMassPt%sBoth2core_cen%d",phiTitle,cent),"(M,p_{T},d#phi)_{#gamma#gamma}",nMm,xM,nPtf,xPt,nPhi,xPhi);
    tmp3->Sumw2() ; fOutputContainer->Add(tmp3) ;

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
void AliAnalysisTaskPi0DiffEfficiency::UserExec(Option_t *) 
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

  const Int_t nMixEvents[6]={4,4,5,10,20,20} ;

 
  //reaction plane
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
  
  TClonesArray *mcArray = (TClonesArray*)event->FindListObject(AliAODMCParticle::StdBranchName());  

  ProcessMC() ;

  if(fPHOSEvent1){
    fPHOSEvent1->Clear() ;
    fPHOSEvent2->Clear() ;
  }
  else{
    fPHOSEvent1 = new TClonesArray("AliCaloPhoton",200) ;
    fPHOSEvent2 = new TClonesArray("AliCaloPhoton",200) ;
  }

  TClonesArray * clustersEmb = (TClonesArray*)event->FindListObject("EmbeddedCaloClusters") ;
  AliAODCaloCells * cellsEmb = (AliAODCaloCells *)event->FindListObject("EmbeddedPHOScells") ;
  TClonesArray * clustersOld = event->GetCaloClusters() ;
  AliAODCaloCells * cellsOld = event->GetPHOSCells() ;
//  TClonesArray *mcArray = (TClonesArray*)event->FindListObject(AliAODMCParticle::StdBranchName());

  
  TVector3 vertex(vtx0);
  char key[55] ;
  //Before Embedding
  Int_t multClustOld = clustersOld->GetEntriesFast();
  Int_t multClustEmb = clustersEmb->GetEntriesFast();
  Int_t inPHOSold=0 ;
  Int_t inPHOSemb=0 ;
  for (Int_t i=0; i<multClustOld; i++) {
    AliAODCaloCluster *clu = (AliAODCaloCluster*)clustersOld->At(i);
    if ( !clu->IsPHOS() || clu->E()<0.3) continue;

    Bool_t survive=kFALSE ;
    for(Int_t ii=0;(ii<multClustEmb)&&(!survive);ii++){
       AliAODCaloCluster *clu2 = (AliAODCaloCluster*)clustersEmb->At(ii);
       survive=IsSameCluster(clu,clu2);
    }


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

    snprintf(key,55,"hCluM%d",mod) ;
    FillHistogram(key,cellX,cellZ,1.);

    TLorentzVector pv1 ;
    clu->GetMomentum(pv1 ,vtx0);
    
    pv1*=Scale(pv1.E()) ;

    if(inPHOSold>=fPHOSEvent1->GetSize()){
      fPHOSEvent1->Expand(inPHOSold+50) ;
    }
    AliCaloPhoton * ph = new((*fPHOSEvent1)[inPHOSold]) AliCaloPhoton(pv1.X(),pv1.Py(),pv1.Z(),pv1.E()) ;
    ph->SetModule(mod) ;
    AliPHOSAodCluster cluPHOS1(*clu);
    cluPHOS1.Recalibrate(fPHOSCalibData,cellsOld); // modify the cell energies
    Double_t ecore=CoreEnergy(&cluPHOS1) ;
    ecore*=Scale(ecore) ;
    pv1*= ecore/pv1.E() ;    
    ph->SetMomV2(&pv1) ;
    ph->SetNCells(clu->GetNCells());
    Double_t m02=0.,m20=0.;
    EvalLambdas(&cluPHOS1,0,m02, m20);   
    ph->SetDispBit(TestLambda(clu->E(),m20,m02)) ;
    EvalLambdas(&cluPHOS1,1,m02, m20);
    ph->SetDisp2Bit(TestLambda2(clu->E(),m20,m02)) ;
    
    ph->SetCPVBit(clu->GetEmcCpvDistance()>2.) ;
    ph->SetCPV2Bit(clu->GetEmcCpvDistance()>4.) ;
    if(!survive) //this cluster found in list after embedding, skipping it
      ph->SetTagged(1) ;
    ph->SetPhoton(clu->GetNExMax()<2); // Remember, if it is unfolded
    ph->SetWeight(1.) ; //All weights for real particles ==1.

    if(!survive){
      Double_t distBC=clu->GetDistanceToBadChannel();
      if(distBC>2.)
        FillHistogram(Form("hPhotAll_DistBad2_cen%d",fCenBin),ph->Pt(),-1.) ;
        if(distBC>4.)
          FillHistogram(Form("hPhotAll_DistBad4_cen%d",fCenBin),ph->Pt(),-1.) ;
          if(distBC>6.)
            FillHistogram(Form("hPhotAll_DistBad6_cen%d",fCenBin),ph->Pt(),-1.) ;
    }
    inPHOSold++ ;
  }

  for (Int_t i=0; i<multClustEmb; i++) {
    AliAODCaloCluster *clu = (AliAODCaloCluster*)clustersEmb->At(i);
    if ( !clu->IsPHOS() || clu->E()<0.3) continue;

    Bool_t survive=kFALSE ;
    for(Int_t ii=0;(ii<multClustOld)&&(!survive);ii++){
       AliAODCaloCluster *clu2 = (AliAODCaloCluster*)clustersOld->At(ii);
       survive=IsSameCluster(clu,clu2);
    }
    
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

    snprintf(key,55,"hCluM%d",mod) ;
    FillHistogram(key,cellX,cellZ,1.);

    TLorentzVector pv1 ;
    clu->GetMomentum(pv1 ,vtx0);

    pv1*=Scale(pv1.E()) ;
    
    if(inPHOSemb>=fPHOSEvent2->GetSize()){
      fPHOSEvent2->Expand(inPHOSemb+50) ;
    }
    AliCaloPhoton * ph = new((*fPHOSEvent2)[inPHOSemb]) AliCaloPhoton(pv1.X(),pv1.Py(),pv1.Z(),pv1.E()) ;
    ph->SetModule(mod) ;
    AliPHOSAodCluster cluPHOS1(*clu);
    cluPHOS1.Recalibrate(fPHOSCalibData,cellsEmb); // modify the cell energies
    Double_t ecore=CoreEnergy(&cluPHOS1) ;
    ecore*=Scale(ecore) ;
    pv1*= ecore/pv1.E() ;    
    ph->SetMomV2(&pv1) ;
    ph->SetNCells(clu->GetNCells());
    Double_t m02=0.,m20=0.;
    EvalLambdas(&cluPHOS1,0,m02, m20);   
    ph->SetDispBit(TestLambda(clu->E(),m20,m02)) ;
    EvalLambdas(&cluPHOS1,1,m02, m20);
    ph->SetDisp2Bit(TestLambda2(clu->E(),m20,m02)) ;

    ph->SetCPVBit(clu->GetEmcCpvDistance()>2.) ;
    ph->SetCPV2Bit(clu->GetEmcCpvDistance()>4.) ;
    if(!survive) //this cluster found in list after embedding, skipping it
      ph->SetTagged(1) ;
    ph->SetPhoton(clu->GetNExMax()<2); // Remember, if it is unfolded
    
    //Set weight for embedded particles
    Double_t w=1. ;
    Int_t iprim = clu->GetLabel() ;
    if(iprim<mcArray->GetEntriesFast() && iprim>-1){    
      AliAODMCParticle* particle =  (AliAODMCParticle*) mcArray->At(iprim);
      iprim=particle->GetMother() ;
      while(iprim>-1){
	particle =  (AliAODMCParticle*) mcArray->At(iprim);
        iprim=particle->GetMother() ;
      }
      if(particle->GetPdgCode()==111){
	Double_t pt = particle->Pt() ;
	w=PrimaryWeight(pt) ;
      }
    }
    ph->SetWeight(w) ;
    

    if(!survive){
      Double_t distBC=clu->GetDistanceToBadChannel();
      if(distBC>2.)
        FillHistogram(Form("hPhotAll_DistBad2_cen%d",fCenBin),ph->Pt()) ;
        if(distBC>4.)
          FillHistogram(Form("hPhotAll_DistBad4_cen%d",fCenBin),ph->Pt()) ;
          if(distBC>6.)
            FillHistogram(Form("hPhotAll_DistBad6_cen%d",fCenBin),ph->Pt()) ;
    }

    inPHOSemb++ ;
  }
  
  //Single photon
  for (Int_t i1=0; i1<inPHOSold; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent1->At(i1) ;
    if(!ph1->IsTagged())
      continue ;
    
    Double_t dphi=ph1->Phi()-fRPfull ;
    while(dphi<0)dphi+=TMath::Pi() ;
    while(dphi>TMath::Pi())dphi-=TMath::Pi() ;
    Double_t pt = ph1->Pt() ;
    Double_t ptV2=ph1->GetMomV2()->Pt() ;
//always 1!    Double_t w=ph1->GetWeight() ;
    
    FillHistogram(Form("hPhotPhiTPCAll_cen%d",fCenBin),pt,dphi,-1.) ;    
    FillHistogram(Form("hPhotPhiTPCAllcore_cen%d",fCenBin),ptV2,dphi,-1.) ;
    
    FillHistogram(Form("hPhotAll_cen%d",fCenBin),pt,-1.) ;
    FillHistogram(Form("hNegPhotAll_cen%d",fCenBin),pt,-1.) ;
    FillHistogram(Form("hPhotAllcore_cen%d",fCenBin),ptV2,-1.) ;
    FillHistogram(Form("hNegPhotAllcore_cen%d",fCenBin),ptV2,-1.) ;
    if(ph1->IsPhoton()){
      FillHistogram(Form("hPhotAllwou_cen%d",fCenBin),pt,-1.) ;      
    }
    if(ph1->IsCPVOK() ){
      FillHistogram(Form("hPhotCPV_cen%d",fCenBin),pt,-1.) ;
      FillHistogram(Form("hNegPhotCPV_cen%d",fCenBin),pt,-1.) ;
      FillHistogram(Form("hPhotCPVcore_cen%d",fCenBin),ptV2,-1.) ;
      FillHistogram(Form("hNegPhotCPVcore_cen%d",fCenBin),ptV2,-1.) ;
      FillHistogram(Form("hPhotPhiTPCCPV_cen%d",fCenBin),pt,dphi,-1.) ;    
      FillHistogram(Form("hPhotPhiTPCCPVcore_cen%d",fCenBin),ptV2,dphi,-1.) ;      
    }
    if(ph1->IsCPV2OK() ){
      FillHistogram(Form("hPhotCPV2_cen%d",fCenBin),pt,-1.) ;
      FillHistogram(Form("hPhotCPV2core_cen%d",fCenBin),ptV2,-1.) ;
      FillHistogram(Form("hNegPhotCPV2_cen%d",fCenBin),pt,-1.) ;      
      FillHistogram(Form("hPhotPhiTPCCPV2_cen%d",fCenBin),pt,dphi,-1.) ;    
      FillHistogram(Form("hPhotPhiTPCCPV2core_cen%d",fCenBin),ptV2,dphi,-1.) ;
    }
    if(ph1->IsDisp2OK()){
      FillHistogram(Form("hPhotDisp2_cen%d",fCenBin),pt,-1.) ;
      FillHistogram(Form("hPhotDisp2core_cen%d",fCenBin),ptV2,-1.) ;
      FillHistogram(Form("hNegPhotDisp2_cen%d",fCenBin),pt,-1.) ;
      FillHistogram(Form("hNegPhotDisp2core_cen%d",fCenBin),ptV2,-1.) ;
      
      FillHistogram(Form("hPhotPhiTPCDisp2_cen%d",fCenBin),pt,dphi,-1.) ;    
      FillHistogram(Form("hPhotPhiTPCDisp2core_cen%d",fCenBin),ptV2,dphi,-1.) ;
      if(ph1->IsCPVOK()){
        FillHistogram(Form("hPhotBoth2_cen%d",fCenBin),pt,-1.) ;
        FillHistogram(Form("hNegPhotBoth2_cen%d",fCenBin),pt,-1.) ;
        FillHistogram(Form("hPhotBoth2core_cen%d",fCenBin),ptV2,-1.) ;
        FillHistogram(Form("hNegPhotBoth2core_cen%d",fCenBin),ptV2,-1.) ;	
        FillHistogram(Form("hPhotPhiTPCBoth2_cen%d",fCenBin),pt,dphi,-1.) ;    
        FillHistogram(Form("hPhotPhiTPCBoth2core_cen%d",fCenBin),ptV2,dphi,-1.) ;
      }
    }
    if(ph1->IsDispOK()){
      FillHistogram(Form("hPhotDisp_cen%d",fCenBin),pt,-1.) ;
      FillHistogram(Form("hNegPhotDisp_cen%d",fCenBin),pt,-1.) ;
      FillHistogram(Form("hPhotDispcore_cen%d",fCenBin),ptV2,-1.) ;
      FillHistogram(Form("hNegPhotDispcore_cen%d",fCenBin),ptV2,-1.) ;
      if(ph1->IsPhoton()){
        FillHistogram(Form("hPhotDispwou_cen%d",fCenBin),pt,-1.) ;      
      }

      FillHistogram(Form("hPhotPhiTPCDisp_cen%d",fCenBin),pt,dphi,-1.) ;    
      FillHistogram(Form("hPhotPhiTPCDispcore_cen%d",fCenBin),ptV2,dphi,-1.) ;
      if(ph1->IsCPVOK()){
        FillHistogram(Form("hPhotBoth_cen%d",fCenBin),pt,-1.) ;
        FillHistogram(Form("hNegPhotBoth_cen%d",fCenBin),pt,-1.) ;
        FillHistogram(Form("hPhotBothcore_cen%d",fCenBin),ptV2,-1.) ;
        FillHistogram(Form("hNegPhotBothcore_cen%d",fCenBin),ptV2,-1.) ;

	FillHistogram(Form("hPhotPhiTPCBoth_cen%d",fCenBin),pt,dphi,-1.) ;    
        FillHistogram(Form("hPhotPhiTPCBothcore_cen%d",fCenBin),ptV2,dphi,-1.) ;
      }
    } 
  } // end of loop i1 

  for (Int_t i1=0; i1<inPHOSemb; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent2->At(i1) ;
    if(!ph1->IsTagged())
      continue ;

    Double_t dphi=ph1->Phi()-fRPfull ;
    while(dphi<0)dphi+=TMath::Pi() ;
    while(dphi>TMath::Pi())dphi-=TMath::Pi() ;
    Double_t pt = ph1->Pt() ;
    Double_t ptV2=ph1->GetMomV2()->Pt() ;
    Double_t w=ph1->GetWeight() ;

    FillHistogram(Form("hPhotPhiTPCAll_cen%d",fCenBin),pt,dphi,w) ;    
    FillHistogram(Form("hPhotPhiTPCAllcore_cen%d",fCenBin),ptV2,dphi,w) ;

    FillHistogram(Form("hPhotAll_cen%d",fCenBin),pt,w) ;
    FillHistogram(Form("hPhotAllcore_cen%d",fCenBin),ptV2,w) ;
    if(ph1->IsPhoton()){
      FillHistogram(Form("hPhotAllwou_cen%d",fCenBin),pt,w) ;      
    }
    if(ph1->IsCPVOK() ){
      FillHistogram(Form("hPhotCPV_cen%d",fCenBin),pt,w) ;
      FillHistogram(Form("hPhotCPVcore_cen%d",fCenBin),ptV2,w) ;
      
      FillHistogram(Form("hPhotPhiTPCCPV_cen%d",fCenBin),pt,dphi,w) ;    
      FillHistogram(Form("hPhotPhiTPCCPVcore_cen%d",fCenBin),ptV2,dphi,w) ;
    }
    if(ph1->IsCPV2OK() ){
      FillHistogram(Form("hPhotCPV2_cen%d",fCenBin),pt,w) ;    
      FillHistogram(Form("hPhotCPV2core_cen%d",fCenBin),ptV2,w) ;    
      FillHistogram(Form("hPhotPhiTPCCPV2_cen%d",fCenBin),pt,dphi,w) ;    
      FillHistogram(Form("hPhotPhiTPCCPV2core_cen%d",fCenBin),ptV2,dphi,w) ;
    }
    if(ph1->IsDisp2OK()){
      FillHistogram(Form("hPhotDisp2_cen%d",fCenBin),pt,w) ;
      FillHistogram(Form("hPhotDisp2core_cen%d",fCenBin),ptV2,w) ;
      
      FillHistogram(Form("hPhotPhiTPCDisp2_cen%d",fCenBin),pt,dphi,w) ;    
      FillHistogram(Form("hPhotPhiTPCDisp2core_cen%d",fCenBin),ptV2,dphi,w) ;
      if(ph1->IsCPVOK()){
        FillHistogram(Form("hPhotBoth2_cen%d",fCenBin),pt,w) ;
        FillHistogram(Form("hPhotBoth2core_cen%d",fCenBin),ptV2,w) ;
	
        FillHistogram(Form("hPhotPhiTPCBoth2_cen%d",fCenBin),pt,dphi,w) ;    
        FillHistogram(Form("hPhotPhiTPCBoth2core_cen%d",fCenBin),ptV2,dphi,w) ;
      }
    }
    if(ph1->IsDispOK()){
      FillHistogram(Form("hPhotDisp_cen%d",fCenBin),pt,w) ;
      FillHistogram(Form("hPhotDispcore_cen%d",fCenBin),ptV2,w) ;
      if(ph1->IsPhoton()){
        FillHistogram(Form("hPhotDispwou_cen%d",fCenBin),pt,w) ;      
      }
      
      FillHistogram(Form("hPhotPhiTPCBoth_cen%d",fCenBin),pt,dphi,w) ;    
      FillHistogram(Form("hPhotPhiTPCBothcore_cen%d",fCenBin),ptV2,dphi,w) ;
      if(ph1->IsCPVOK()){
        FillHistogram(Form("hPhotBoth_cen%d",fCenBin),pt,w) ;
        FillHistogram(Form("hPhotBothcore_cen%d",fCenBin),ptV2,w) ;
	
        FillHistogram(Form("hPhotPhiTPCBoth_cen%d",fCenBin),pt,dphi,w) ;    
        FillHistogram(Form("hPhotPhiTPCBothcore_cen%d",fCenBin),ptV2,dphi,w) ;
      }
    } // end of loop i2
  } // end of loop i1 


  const Double_t prob[10]={0.1,0.2,0.3,1.,1.,1.,1.,1.,1.,1.} ; //Probabilities to accept Tagged+Bg pair


  // Fill Real disribution:
  // Disappeared clusters enter with negative contribution
  // In addition fill control histogam with Real before embedding
  for (Int_t i1=0; i1<inPHOSold-1; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent1->At(i1) ;
    Double_t w1 = ph1->GetWeight() ;
    for (Int_t i2=i1+1; i2<inPHOSold; i2++) {
      AliCaloPhoton * ph2=(AliCaloPhoton*)fPHOSEvent1->At(i2) ;

      TLorentzVector p12  = *ph1  + *ph2;
      TLorentzVector pv12 = *(ph1->GetMomV2()) + *(ph2->GetMomV2());      
      Double_t w2 = ph2->GetWeight() ;
      Double_t w = TMath::Sqrt(w1*w2) ;
      
      Double_t dphi=p12.Phi()-fRPfull ;
      while(dphi<0)dphi+=TMath::Pi() ;
      while(dphi>TMath::Pi())dphi-=TMath::Pi() ;
      Double_t m=p12.M() ;
      Double_t mV2=pv12.M() ;
      Double_t pt = p12.Pt() ;
      Double_t ptV2 = pv12.Pt() ;
      
      
      //Fill Controll histogram: Real before embedding
      FillHistogram(Form("hOldMassPtAll_cen%d",fCenBin),m,pt,-w) ;
      FillHistogram(Form("hOldMassPtAllcore_cen%d",fCenBin),mV2,ptV2,-w) ;
      if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	FillHistogram(Form("hOldMassPtCPV_cen%d",fCenBin),m ,pt,-w) ;
	FillHistogram(Form("hOldMassPtCPV_cen%d",fCenBin),mV2, ptV2,-w) ;
      }
      if(ph1->IsCPV2OK() && ph2->IsCPV2OK()){
	FillHistogram(Form("hOldMassPtCPV2_cen%d",fCenBin),m ,pt,-w) ;
      }
      if(ph1->IsDisp2OK() && ph2->IsDisp2OK()){
	FillHistogram(Form("hOldMassPtDisp2_cen%d",fCenBin),m ,pt,-w) ;
	FillHistogram(Form("hOldMassPtDisp2core_cen%d",fCenBin),mV2 ,ptV2,-w) ;
	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  FillHistogram(Form("hOldMassPtBoth2_cen%d",fCenBin),m ,pt,-w) ;
	  FillHistogram(Form("hOldMassPtBoth2core_cen%d",fCenBin),mV2 ,ptV2,-w) ;
	}
      }
      if(ph1->IsDispOK() && ph2->IsDispOK()){
	FillHistogram(Form("hOldMassPtDisp_cen%d",fCenBin),m ,pt,-w) ;
	FillHistogram(Form("hOldMassPtDispcore_cen%d",fCenBin),mV2 ,ptV2,-w) ;
	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  FillHistogram(Form("hOldMassPtBoth_cen%d",fCenBin),m ,pt,-w) ;
	  FillHistogram(Form("hOldMassPtBothcore_cen%d",fCenBin),mV2 ,ptV2,-w) ;
        }
      }
      
      //Now fill main histograms with negative contributions
      if(!(ph1->IsTagged() || ph2->IsTagged()) )
        continue ;
      if(!ph1->IsTagged() || !ph2->IsTagged()){ //Tagged + Bg combination
        if(gRandom->Uniform()>prob[fCenBin])
          continue ;
      }
      FillHistogram(Form("hMassPtAll_cen%d",fCenBin),m ,pt,-w) ;
      FillHistogram(Form("hMassPtAllcore_cen%d",fCenBin),mV2 ,ptV2,-w) ;
      if(ph1->IsPhoton()&&ph2->IsPhoton()){
        FillHistogram(Form("hMassPtAllwou_cen%d",fCenBin),m,pt,-w) ;
      }

      FillHistogram(Form("hMassPtTPCAll_cen%d",fCenBin),m,pt,dphi,-w) ;
      FillHistogram(Form("hMassPtTPCAllcore_cen%d",fCenBin),mV2,ptV2,dphi,-w) ;
      
      FillHistogram(Form("hNegMassPtAll_cen%d",fCenBin),m ,pt,-w) ;
      FillHistogram(Form("hNegMassPtAllcore_cen%d",fCenBin),mV2 ,ptV2,-w) ;
      
      if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	FillHistogram(Form("hMassPtCPV_cen%d",fCenBin),m ,pt,-w) ;
	FillHistogram(Form("hMassPtCPVcore_cen%d",fCenBin),mV2 ,ptV2,-w) ;

	FillHistogram(Form("hMassPtTPCCPV_cen%d",fCenBin),m,pt,dphi,-w) ;
        FillHistogram(Form("hMassPtTPCCPVcore_cen%d",fCenBin),mV2,ptV2,dphi,-w) ;

	FillHistogram(Form("hNegMassPtCPV_cen%d",fCenBin),m ,pt,-w) ;
	FillHistogram(Form("hNegMassPtCPVcore_cen%d",fCenBin),mV2 ,ptV2,-w) ;
      }
      if(ph1->IsCPV2OK() && ph2->IsCPV2OK()){
	FillHistogram(Form("hMassPtCPV2_cen%d",fCenBin),m ,pt,-w) ;
	FillHistogram(Form("hNegMassPtCPV2_cen%d",fCenBin),m ,pt,-w) ;
	
        FillHistogram(Form("hMassPtTPCCPV2_cen%d",fCenBin),m,pt,dphi,-w) ;
        FillHistogram(Form("hMassPtTPCCPV2core_cen%d",fCenBin),mV2,ptV2,dphi,-w) ;
      }
      
      if(ph1->IsDisp2OK() && ph2->IsDisp2OK()){
	FillHistogram(Form("hMassPtDisp2_cen%d",fCenBin),m ,pt,-w) ;
	FillHistogram(Form("hMassPtDisp2core_cen%d",fCenBin),mV2 ,ptV2,-w) ;
	
        FillHistogram(Form("hMassPtTPCDisp2_cen%d",fCenBin),m,pt,dphi,-w) ;
        FillHistogram(Form("hMassPtTPCDisp2core_cen%d",fCenBin),mV2,ptV2,dphi,-w) ;

	FillHistogram(Form("hNegMassPtDisp2_cen%d",fCenBin),m ,pt,-w) ;
	FillHistogram(Form("hNegMassPtDisp2core_cen%d",fCenBin),mV2 ,ptV2,-w) ;
	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  FillHistogram(Form("hMassPtBoth2_cen%d",fCenBin),m ,pt,-w) ;
	  FillHistogram(Form("hMassPtBoth2core_cen%d",fCenBin),mV2 ,ptV2,-w) ;
	  FillHistogram(Form("hNegMassPtBoth2_cen%d",fCenBin),m ,pt,-w) ;
	  FillHistogram(Form("hNegMassPtBoth2core_cen%d",fCenBin),mV2 ,ptV2,-w) ;
	  
          FillHistogram(Form("hMassPtTPCBoth2_cen%d",fCenBin),m,pt,dphi,-w) ;
          FillHistogram(Form("hMassPtTPCBoth2core_cen%d",fCenBin),mV2,ptV2,dphi,-w) ;
        }
      }
      if(ph1->IsDispOK() && ph2->IsDispOK()){
	FillHistogram(Form("hMassPtDisp_cen%d",fCenBin),m ,pt,-w) ;
	FillHistogram(Form("hMassPtDispcore_cen%d",fCenBin),mV2 ,ptV2,-w) ;
        if(ph1->IsPhoton()&&ph2->IsPhoton()){
          FillHistogram(Form("hMassPtDispwou_cen%d",fCenBin),m,pt,-w) ;
        }
	FillHistogram(Form("hNegMassPtDisp_cen%d",fCenBin),m ,pt,-w) ;
	FillHistogram(Form("hNegMassPtDispcore_cen%d",fCenBin),mV2 ,ptV2,-w) ;
	
        FillHistogram(Form("hMassPtTPCDisp_cen%d",fCenBin),m,pt,dphi,-w) ;
        FillHistogram(Form("hMassPtTPCDispcore_cen%d",fCenBin),mV2,ptV2,dphi,-w) ;
	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  FillHistogram(Form("hMassPtBoth_cen%d",fCenBin),m ,pt,-w) ;
	  FillHistogram(Form("hMassPtBothcore_cen%d",fCenBin),mV2 ,ptV2,-w) ;
	  FillHistogram(Form("hNegMassPtBoth_cen%d",fCenBin),m ,pt,-w) ;
	  FillHistogram(Form("hNegMassPtBothcore_cen%d",fCenBin),mV2 ,ptV2,-w) ;
	  
          FillHistogram(Form("hMassPtTPCBoth_cen%d",fCenBin),m,pt,dphi,-w) ;
          FillHistogram(Form("hMassPtTPCBothcore_cen%d",fCenBin),mV2,ptV2,dphi,-w) ;
        }
      }
    } // end of loop i2
  } // end of loop i1 


  // Further fill Real disribution
  // now with positive contribution from new clusters
  // ass well fill controll histogram
  for (Int_t i1=0; i1<inPHOSemb-1; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent2->At(i1) ;
    Double_t w1 = ph1->GetWeight() ;
    for (Int_t i2=i1+1; i2<inPHOSemb; i2++) {
      AliCaloPhoton * ph2=(AliCaloPhoton*)fPHOSEvent2->At(i2) ;

      TLorentzVector p12  = *ph1  + *ph2;
      TLorentzVector pv12 = *(ph1->GetMomV2()) + *(ph2->GetMomV2());      
      Double_t m=p12.M() ;
      Double_t mV2=pv12.M() ;
      Double_t pt = p12.Pt() ;
      Double_t ptV2 = pv12.Pt() ;
      Double_t w2 = ph2->GetWeight() ;
      Double_t w = TMath::Sqrt(w1*w2) ;

      // Controll histogram: Real after embedding
      FillHistogram(Form("hNewMassPtAll_cen%d",fCenBin),m ,pt,w) ;
      FillHistogram(Form("hNewMassPtAllcore_cen%d",fCenBin),mV2 ,ptV2,w) ;
      if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	FillHistogram(Form("hNewMassPtCPV_cen%d",fCenBin),m ,pt,w) ;
	FillHistogram(Form("hNewMassPtCPVcore_cen%d",fCenBin),mV2 ,ptV2,w) ;
      }
      if(ph1->IsCPV2OK() && ph2->IsCPV2OK()){
	FillHistogram(Form("hNewMassPtCPV2_cen%d",fCenBin),m ,pt,w) ;
      }
      if(ph1->IsDisp2OK() && ph2->IsDisp2OK()){
	FillHistogram(Form("hNewMassPtDisp2_cen%d",fCenBin),m ,pt,w) ;
	FillHistogram(Form("hNewMassPtDisp2core_cen%d",fCenBin),mV2 ,ptV2,w) ;
	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  FillHistogram(Form("hNewMassPtBoth2_cen%d",fCenBin),m ,pt,w) ;
	  FillHistogram(Form("hNewMassPtBoth2core_cen%d",fCenBin),mV2,ptV2,w) ;
        }
      }
      if(ph1->IsDispOK() && ph2->IsDispOK()){
	FillHistogram(Form("hNewMassPtDisp_cen%d",fCenBin),m ,pt,w) ;
	FillHistogram(Form("hNewMassPtDispcore_cen%d",fCenBin),mV2 ,ptV2,w) ;
	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  FillHistogram(Form("hNewMassPtBoth_cen%d",fCenBin),m ,pt,w) ;
	  FillHistogram(Form("hNewMassPtBothcore_cen%d",fCenBin),mV2 ,ptV2,w) ;
        }
      }
     
      //Now fill main histogamm
      //new clusters with positive contribution
      if(!(ph1->IsTagged() || ph2->IsTagged()) )
        continue ;
      if(!ph1->IsTagged() || !ph2->IsTagged()){ //Tagged + Bg combination
        if(gRandom->Uniform()>prob[fCenBin])
          continue ;
      }

      Double_t dphi=p12.Phi()-fRPfull ;
      while(dphi<0)dphi+=TMath::Pi() ;
      while(dphi>TMath::Pi())dphi-=TMath::Pi() ;
      
      FillHistogram(Form("hMassPtAll_cen%d",fCenBin),m ,pt,w) ;
      FillHistogram(Form("hMassPtAllcore_cen%d",fCenBin),mV2 ,ptV2,w) ;
      if(ph1->IsPhoton()&&ph2->IsPhoton()){
        FillHistogram(Form("hMassPtAllwou_cen%d",fCenBin),m,pt,w) ;
      }

      FillHistogram(Form("hMassPtTPCAll_cen%d",fCenBin),m,pt,dphi,w) ;
      FillHistogram(Form("hMassPtTPCAllcore_cen%d",fCenBin),mV2,ptV2,dphi,w) ;

      if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	FillHistogram(Form("hMassPtCPV_cen%d",fCenBin),m ,pt,w) ;
	FillHistogram(Form("hMassPtCPVcore_cen%d",fCenBin),mV2 ,ptV2,w) ;

	FillHistogram(Form("hMassPtTPCCPV_cen%d",fCenBin),m,pt,dphi,w) ;
        FillHistogram(Form("hMassPtTPCCPVcore_cen%d",fCenBin),mV2,ptV2,dphi,w) ;
      }
      if(ph1->IsCPV2OK() && ph2->IsCPV2OK()){
	FillHistogram(Form("hMassPtCPV2_cen%d",fCenBin),m ,pt,w) ;

	FillHistogram(Form("hMassPtTPCCPV2_cen%d",fCenBin),m,pt,dphi,w) ;
        FillHistogram(Form("hMassPtTPCCPV2core_cen%d",fCenBin),mV2,ptV2,dphi,w) ;
      }
      if(ph1->IsDisp2OK() && ph2->IsDisp2OK()){
	FillHistogram(Form("hMassPtDisp2_cen%d",fCenBin),m ,pt,w) ;
	FillHistogram(Form("hMassPtDisp2core_cen%d",fCenBin),mV2 ,ptV2,w) ;

	FillHistogram(Form("hMassPtTPCDisp2_cen%d",fCenBin),m,pt,dphi,w) ;
        FillHistogram(Form("hMassPtTPCDisp2core_cen%d",fCenBin),mV2,ptV2,dphi,w) ;
	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  FillHistogram(Form("hMassPtBoth2_cen%d",fCenBin),m ,pt,w) ;
	  FillHistogram(Form("hMassPtBoth2core_cen%d",fCenBin),mV2 ,ptV2,w) ;

	  FillHistogram(Form("hMassPtTPCBoth2_cen%d",fCenBin),m,pt,dphi,w) ;
          FillHistogram(Form("hMassPtTPCBoth2core_cen%d",fCenBin),mV2,ptV2,dphi,w) ;
        }
     }
      if(ph1->IsDispOK() && ph2->IsDispOK()){
	FillHistogram(Form("hMassPtDisp_cen%d",fCenBin),m ,pt,w) ;
	FillHistogram(Form("hMassPtDispcore_cen%d",fCenBin),mV2 ,ptV2,w) ;
        if(ph1->IsPhoton()&&ph2->IsPhoton()){
          FillHistogram(Form("hMassPtDispwou_cen%d",fCenBin),m,pt,w) ;
        }

	FillHistogram(Form("hMassPtTPCDisp_cen%d",fCenBin),m,pt,dphi,w) ;
        FillHistogram(Form("hMassPtTPCDispcore_cen%d",fCenBin),mV2,ptV2,dphi,w) ;

	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  FillHistogram(Form("hMassPtBoth_cen%d",fCenBin),m ,pt,w) ;
	  FillHistogram(Form("hMassPtBothcore_cen%d",fCenBin),mV2 ,ptV2,w) ;
	  
	  FillHistogram(Form("hMassPtTPCBoth_cen%d",fCenBin),m,pt,dphi,w) ;
          FillHistogram(Form("hMassPtTPCBothcore_cen%d",fCenBin),mV2,ptV2,dphi,w) ;
        }
      }
    } // end of loop i2
  } // end of loop i1 


  //now mixed, does not really matter old or new list of clusters
  for (Int_t i1=0; i1<inPHOSemb; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent2->At(i1) ;
    for(Int_t ev=0; ev<prevPHOS->GetSize();ev++){
      TClonesArray * mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev)) ;
      for(Int_t i2=0; i2<mixPHOS->GetEntriesFast();i2++){
	AliCaloPhoton * ph2=(AliCaloPhoton*)mixPHOS->At(i2) ;
	TLorentzVector p12  = *ph1  + *ph2;
	TLorentzVector pv12 = *(ph1->GetMomV2()) + *(ph2->GetMomV2());

        Double_t dphi=p12.Phi()-fRPfull ;
        while(dphi<0)dphi+=TMath::Pi() ;
        while(dphi>TMath::Pi())dphi-=TMath::Pi() ;
        Double_t m=p12.M() ;
        Double_t mV2=pv12.M() ;
        Double_t pt = p12.Pt() ;
        Double_t ptV2 = pv12.Pt() ;
	
	FillHistogram(Form("hMiMassPtAll_cen%d",fCenBin),m ,pt,1.) ;
	FillHistogram(Form("hMiMassPtAllcore_cen%d",fCenBin),mV2 ,ptV2,1.) ;
        if(ph1->IsPhoton()&&ph2->IsPhoton()){
          FillHistogram(Form("hMiMassPtAllwou_cen%d",fCenBin),m,pt,1.) ;
        }

        FillHistogram(Form("hMiMassPtTPCAll_cen%d",fCenBin),m,pt,dphi) ;
        FillHistogram(Form("hMiMassPtTPCAllcore_cen%d",fCenBin),mV2,ptV2,dphi) ;

	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  FillHistogram(Form("hMiMassPtCPV_cen%d",fCenBin),m ,pt,1.) ;
	  FillHistogram(Form("hMiMassPtCPVcore_cen%d",fCenBin),mV2 ,ptV2,1.) ;

          FillHistogram(Form("hMiMassPtTPCCPV_cen%d",fCenBin),m,pt,dphi) ;
          FillHistogram(Form("hMiMassPtTPCCPVcore_cen%d",fCenBin),mV2,ptV2,dphi) ;

	}
	if(ph1->IsCPV2OK() && ph2->IsCPV2OK()){
	  FillHistogram(Form("hMiMassPtCPV2_cen%d",fCenBin),m ,pt,1.) ;

          FillHistogram(Form("hMiMassPtTPCCPV2_cen%d",fCenBin),m,pt,dphi) ;
          FillHistogram(Form("hMiMassPtTPCCPV2core_cen%d",fCenBin),mV2,ptV2,dphi) ;
	}
	if(ph1->IsDisp2OK() && ph2->IsDisp2OK()){
	  FillHistogram(Form("hMiMassPtDisp2_cen%d",fCenBin),m ,pt,1.) ;
	  FillHistogram(Form("hMiMassPtDisp2core_cen%d",fCenBin),mV2 ,ptV2,1.) ;

          FillHistogram(Form("hMiMassPtTPCDisp2_cen%d",fCenBin),m,pt,dphi) ;
          FillHistogram(Form("hMiMassPtTPCDisp2core_cen%d",fCenBin),mV2,ptV2,dphi) ;

	  if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	    FillHistogram(Form("hMiMassPtBoth2_cen%d",fCenBin),m ,pt,1.) ;
	    FillHistogram(Form("hMiMassPtBoth2core_cen%d",fCenBin),mV2 ,ptV2,1.) ;

            FillHistogram(Form("hMiMassPtTPCBoth2_cen%d",fCenBin),m,pt,dphi) ;
            FillHistogram(Form("hMiMassPtTPCBoth2core_cen%d",fCenBin),mV2,ptV2,dphi) ;
	  }
	}
	if(ph1->IsDispOK() && ph2->IsDispOK()){
	  FillHistogram(Form("hMiMassPtDisp_cen%d",fCenBin),m ,pt,1.) ;
	  FillHistogram(Form("hMiMassPtDispcore_cen%d",fCenBin),mV2 ,ptV2,1.) ;
          if(ph1->IsPhoton()&&ph2->IsPhoton()){
            FillHistogram(Form("hMiMassPtDispwou_cen%d",fCenBin),m,pt,1.) ;
          }

          FillHistogram(Form("hMiMassPtTPCDisp_cen%d",fCenBin),m,pt,dphi) ;
          FillHistogram(Form("hMiMassPtTPCDispcore_cen%d",fCenBin),mV2,ptV2,dphi) ;
	  
	  if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	    FillHistogram(Form("hMiMassPtBoth_cen%d",fCenBin),m ,pt,1.) ;
	    FillHistogram(Form("hMiMassPtBothcore_cen%d",fCenBin),mV2 ,ptV2,1.) ;

            FillHistogram(Form("hMiMassPtTPCBoth_cen%d",fCenBin),m,pt,dphi) ;
            FillHistogram(Form("hMiMassPtTPCBothcore_cen%d",fCenBin),mV2,ptV2,dphi) ;
	  }
	}
      } // end of loop i2
    }
  } // end of loop i1

  
  
  //Now we either add current events to stack or remove
  //If no photons in current event - no need to add it to mixed
  if(fPHOSEvent2->GetEntriesFast()>0){
    prevPHOS->AddFirst(fPHOSEvent2) ;
    fPHOSEvent2=0;
    delete fPHOSEvent1;
    fPHOSEvent1=0;
    if(prevPHOS->GetSize()>nMixEvents[fCenBin]){//Remove redundant events
      TClonesArray * tmp = static_cast<TClonesArray*>(prevPHOS->Last()) ;
      prevPHOS->RemoveLast() ;
      delete tmp ;
    }
  }
  // Post output data.
  PostData(1, fOutputContainer);
  fEventCounter++;
}
//___________________________________________________________________________
Bool_t AliAnalysisTaskPi0DiffEfficiency::IsSameCluster(AliAODCaloCluster * c1, AliAODCaloCluster * c2)const{
 //Compare clusters before and after embedding
 //clusters are the same if 
 // - Energy changed less than 0.1%  (numerical accuracy in reconstruction)
 // - lists of digits are the same
  
 if(c1->GetNCells() != c2->GetNCells())
   return kFALSE ;
 
 if(TMath::Abs(c1->E()-c2->E())>0.01*c1->E())
   return kFALSE ;

 UShort_t *list1 = c1->GetCellsAbsId() ; 
 UShort_t *list2 = c2->GetCellsAbsId() ; 
 for(Int_t i=0; i<c1->GetNCells(); i++){
  if(list1[i] != list2[i])
    return kFALSE ;
 }
 return kTRUE ; 
  
}
//____________________________________________________________________________
void  AliAnalysisTaskPi0DiffEfficiency::EvalLambdas(AliAODCaloCluster * clu, Int_t iR,Double_t &m02, Double_t &m20){ 
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
void AliAnalysisTaskPi0DiffEfficiency::ProcessMC(){
  //fill histograms for efficiensy etc. calculation
  const Double_t rcut = 1. ; //cut for primary particles
  //---------First pi0/eta-----------------------------
  char partName[10] ;

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
    Double_t w = PrimaryWeight(pt) ;
    //Total number of pi0 with creation radius <1 cm
    FillHistogram(Form("hMC_all_%s_cen%d",partName,fCenBin),pt,w) ;
    if(TMath::Abs(particle->Y())<0.12){
      FillHistogram(Form("hMC_unitEta_%s_cen%d",partName,fCenBin),pt,w) ;
    }

    FillHistogram(Form("hMC_rap_%s_cen%d",partName,fCenBin),particle->Y(),w) ;
    
    Double_t phi=particle->Phi() ;
    while(phi<0.)phi+=TMath::TwoPi() ;
    while(phi>TMath::TwoPi())phi-=TMath::TwoPi() ;
    FillHistogram(Form("hMC_phi_%s_cen%d",partName,fCenBin),phi,w) ;
   
  }
 
  //Now calculate "Real" distribution of clusters with primary
  TClonesArray * clusters = (TClonesArray*)event->FindListObject("EmbeddedCaloClusters") ;
  AliAODCaloCells * cellsEmb = (AliAODCaloCells *)event->FindListObject("EmbeddedPHOScells") ;
  Int_t multClust = clusters->GetEntriesFast();
  TClonesArray cluPrim("AliCaloPhoton",multClust) ; //clusters with primary
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
    
    pv1*=Scale(pv1.E()) ;    

    
    
    if(inPHOS>=cluPrim.GetSize()){
      cluPrim.Expand(inPHOS+50) ;
    }
    AliCaloPhoton * ph = new(cluPrim[inPHOS]) AliCaloPhoton(pv1.X(),pv1.Py(),pv1.Z(),pv1.E()) ;
    //AliCaloPhoton * ph = (AliCaloPhoton*)fPHOSEvent->At(inPHOS) ;
    ph->SetModule(mod) ;
 
    AliPHOSAodCluster cluPHOS1(*clu);
    cluPHOS1.Recalibrate(fPHOSCalibData,cellsEmb); // modify the cell energies
    Double_t ecore=CoreEnergy(&cluPHOS1) ;
    ecore*=Scale(ecore) ;
    pv1*= ecore/pv1.E() ;    
    ph->SetMomV2(&pv1) ;
    ph->SetNCells(clu->GetNCells());
    Double_t m02=0.,m20=0.;
    EvalLambdas(&cluPHOS1,0,m02, m20);   
    ph->SetDispBit(TestLambda(clu->E(),m20,m02)) ;
    EvalLambdas(&cluPHOS1,1,m02, m20);
    ph->SetDisp2Bit(TestLambda2(clu->E(),m20,m02)) ;

    ph->SetCPVBit(clu->GetEmcCpvDistance()>2.) ; //radius in sigmas
    ph->SetCPV2Bit(clu->GetEmcCpvDistance()>4.) ;
    ph->SetPhoton(clu->GetNExMax()<2); // Remember, if it is unfolded
    Double_t w=1. ;
    Int_t iprim = clu->GetLabel() ;
    if(iprim<mcArray->GetEntriesFast() && iprim>-1){    
      AliAODMCParticle* particle =  (AliAODMCParticle*) mcArray->At(iprim);
      iprim=particle->GetMother() ;
      while(iprim>-1){
	particle =  (AliAODMCParticle*) mcArray->At(iprim);
        iprim=particle->GetMother() ;
      }
      if(particle->GetPdgCode()==111){
	Double_t pt = particle->Pt() ;
	w=PrimaryWeight(pt) ;
      }
    }
    ph->SetWeight(w) ;

    inPHOS++ ;
  }
  
  //Single photon
  char key[55] ;
  for (Int_t i1=0; i1<inPHOS; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)cluPrim.At(i1) ;
    Double_t pt=ph1->Pt() ;
    Double_t ptV2=ph1->GetMomV2()->Pt() ;
    Double_t w = ph1->GetWeight() ;
    FillHistogram(Form("hMCPhotAll_cen%d",fCenBin),pt,w) ;
    FillHistogram(Form("hMCPhotAllcore_cen%d",fCenBin),ptV2,w) ;
    if(ph1->IsPhoton()){
      FillHistogram(Form("hMCPhotAllwou_cen%d",fCenBin),pt,w) ;
    }
    if(ph1->IsCPVOK() ){
      FillHistogram(Form("hMCPhotCPV_cen%d",fCenBin),pt,w) ;
      FillHistogram(Form("hMCPhotCPVcore_cen%d",fCenBin),ptV2,w) ;     
    }
    if(ph1->IsCPV2OK() ){
      snprintf(key,55,"hMCPhotCPV2_cen%d",fCenBin) ;
      FillHistogram(Form("hMCPhotCPV2_cen%d",fCenBin),pt,w) ;   
    }
    if(ph1->IsDisp2OK()){
      FillHistogram(Form("hMCPhotDisp2_cen%d",fCenBin),pt,w) ;
      FillHistogram(Form("hMCPhotDisp2core_cen%d",fCenBin),ptV2,w) ;
      if(ph1->IsCPVOK()){
	FillHistogram(Form("hMCPhotBoth2_cen%d",fCenBin),pt,w) ;
	FillHistogram(Form("hMCPhotBoth2core_cen%d",fCenBin),ptV2,w) ;
      }
    }
    if(ph1->IsDispOK()){
      FillHistogram(Form("hMCPhotDisp_cen%d",fCenBin),pt,w) ;
      FillHistogram(Form("hMCPhotDispcore_cen%d",fCenBin),ptV2,w) ;
      if(ph1->IsPhoton()){
        FillHistogram(Form("hMCPhotDispwou_cen%d",fCenBin),pt,w) ;
      }
      if(ph1->IsCPVOK()){
	FillHistogram(Form("hMCPhotBoth_cen%d",fCenBin),pt,w) ;
	FillHistogram(Form("hMCPhotBothcore_cen%d",fCenBin),ptV2,w) ;
      }
    } // end of loop i2
  } // end of loop i1 

  // Fill Real disribution
  for (Int_t i1=0; i1<inPHOS-1; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)cluPrim.At(i1) ;
    Double_t w1 = ph1->GetWeight() ;
    for (Int_t i2=i1+1; i2<inPHOS; i2++) {
      AliCaloPhoton * ph2=(AliCaloPhoton*)cluPrim.At(i2) ;
      Double_t w2 = ph2->GetWeight() ;
      Double_t w = TMath::Sqrt(w1*w2) ;
      TLorentzVector p12  = *ph1  + *ph2;
      TLorentzVector pv12 = *(ph1->GetMomV2()) + *(ph2->GetMomV2());     
      Double_t m=p12.M() ;
      Double_t pt=p12.Pt() ;
      Double_t mV2=pv12.M() ;
      Double_t ptV2=pv12.Pt() ;
       
      FillHistogram(Form("hMCMassPtAll_cen%d",fCenBin),m,pt,w) ;
      snprintf(key,55,"hMCMassPtAllcore_cen%d",fCenBin) ;
      FillHistogram(Form("hMCMassPtAllcore_cen%d",fCenBin),mV2,ptV2,w) ;
      if(ph1->IsPhoton()&&ph2->IsPhoton() ){
        FillHistogram(Form("hMCMassPtAllwou_cen%d",fCenBin),m,pt,w) ;	
      }
      
      if(ph1->Module()==1 && ph2->Module()==1)
	FillHistogram("hMCPi0M11",m,pt,w);
      else if(ph1->Module()==2 && ph2->Module()==2)
	FillHistogram("hMCPi0M22",m,pt,w);
      else if(ph1->Module()==3 && ph2->Module()==3)
	FillHistogram("hMCPi0M33",m,pt,w);
      else if(ph1->Module()==1 && ph2->Module()==2)
	FillHistogram("hMCPi0M12",m,pt,w);
      else if(ph1->Module()==1 && ph2->Module()==3)
	FillHistogram("hMCPi0M13",m,pt,w);
      else if(ph1->Module()==2 && ph2->Module()==3)
	FillHistogram("hMCPi0M23",m,pt,w);    
      
      if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	FillHistogram(Form("hMCMassPtCPV_cen%d",fCenBin),m,pt,w) ;
	FillHistogram(Form("hMCMassPtCPVcore_cen%d",fCenBin),mV2,ptV2,w) ;
      }
      if(ph1->IsCPV2OK() && ph2->IsCPV2OK()){
	FillHistogram(Form("hMCMassPtCPV2core_cen%d",fCenBin),mV2,ptV2,w) ;
      }
      if(ph1->IsDisp2OK() && ph2->IsDisp2OK()){
	FillHistogram(Form("hMCMassPtDisp2_cen%d",fCenBin),m,pt,w) ;
	FillHistogram(Form("hMCMassPtDisp2core_cen%d",fCenBin),mV2,ptV2,w) ;
	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  FillHistogram(Form("hMCMassPtBoth2_cen%d",fCenBin),m,pt,w) ;
	  FillHistogram(Form("hMCMassPtBoth2core_cen%d",fCenBin),mV2,ptV2,w) ;
        }
      }
      if(ph1->IsDispOK() && ph2->IsDispOK()){
	FillHistogram(Form("hMCMassPtDisp_cen%d",fCenBin),m,pt,w) ;
	FillHistogram(Form("hMCMassPtDispcore_cen%d",fCenBin),mV2,ptV2,w) ;
        if(ph1->IsPhoton()&& ph2->IsPhoton()){
	  FillHistogram(Form("hMCMassPtDispwou_cen%d",fCenBin),m,pt,w) ;
	}
	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  FillHistogram(Form("hMCMassPtBoth_cen%d",fCenBin),m,pt,w) ;
	  FillHistogram(Form("hMCMassPtBothcore_cen%d",fCenBin),mV2,ptV2,w) ;
        }
      }
    } // end of loop i2
  } // end of loop i1
}
//____________________________________________________________________________
Double_t  AliAnalysisTaskPi0DiffEfficiency::CoreEnergy(AliPHOSAodCluster * clu){  
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
//_____________________________________________________________________________
Bool_t AliAnalysisTaskPi0DiffEfficiency::TestLambda(Double_t pt,Double_t l1,Double_t l2){
  
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
Bool_t AliAnalysisTaskPi0DiffEfficiency::TestLambda2(Double_t pt,Double_t l1,Double_t l2){
  
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
Double_t AliAnalysisTaskPi0DiffEfficiency::TestCPV(Double_t dx, Double_t dz, Double_t pt, Int_t charge){
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
Bool_t AliAnalysisTaskPi0DiffEfficiency::TestTOF(Double_t t, Double_t e){

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
//_____________________________________________________________________________
void AliAnalysisTaskPi0DiffEfficiency::FillHistogram(const char * key,Double_t x)const{
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
void AliAnalysisTaskPi0DiffEfficiency::FillHistogram(const char * key,Double_t x,Double_t y)const{
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
void AliAnalysisTaskPi0DiffEfficiency::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z) const{
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
void AliAnalysisTaskPi0DiffEfficiency::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z, Double_t w) const{
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
//_____________________________________________________________________________
Bool_t AliAnalysisTaskPi0DiffEfficiency::IsGoodChannel(const char * det, Int_t mod, Int_t ix, Int_t iz)
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
Double_t AliAnalysisTaskPi0DiffEfficiency::PrimaryWeight(Double_t x){
  
   Double_t w=1 ;
   
   
   //Parameterization of LHC10h data from Jan 2013 (pi0 spectrum)
   //Should be consistend with spectrum parameterization used in simulation 
   if(fCenBin==0) //0-5
     w = (0.561741+0.332841*x-0.007082*x*x)/(1.-0.447804*x+0.157830*x*x)+0.080394*x ;
   if(fCenBin==1) //5-10
     w = (0.659096+0.101701*x+0.042395*x*x)/(1.-0.470110*x+0.154665*x*x)+0.052932*x ;
   if(fCenBin==2) //10-20
     w = (0.615575+0.005621*x+0.069263*x*x)/(1.-0.485422*x+0.160822*x*x)+0.040865*x ; 
   if(fCenBin==3) //20-40
     w = (0.441240+0.158358*x+0.059458*x*x)/(1.-0.332609*x+0.147528*x*x)+0.037926*x ; 
   if(fCenBin==4) //40-60
     w = (0.467895-0.001113*x+0.029610*x*x)/(1.-0.266502*x+0.065105*x*x)+0.025431*x ; 
   if(fCenBin==5) //60-80
     w = (0.465204-0.139362*x+0.043500*x*x)/(1.-0.371689*x+0.067606*x*x)+0.006519*x ;

/*
  //Parameterization of photon spectrum 25.02
  if(fCenBin==0) //0-5
     w=(0.870487-0.494032*x+0.076334*x*x+0.001065*x*x*x)/(1.-0.646014*x+0.113839*x*x); 
  if(fCenBin==1) //5-10
     w=(-8.310403+15.767226*x-2.167756*x*x+0.184356*x*x*x)/(1.+4.556793*x+0.980941*x*x); 
  if(fCenBin==2) //10-20
     w=(-5.281594+7.477165*x-0.688609*x*x+0.097601*x*x*x)/(1.+1.102693*x+0.882454*x*x); 
  if(fCenBin==3) //20-40
     w=(0.789472-4.750155*x+4.381545*x*x-0.029239*x*x*x)/(1.-3.738304*x+3.328318*x*x); 
  if(fCenBin==4) //40-60
     w=(0.876792-0.491379*x+0.130872*x*x-0.001390*x*x*x)/(1.+-0.511934*x+0.112705*x*x); 
  if(fCenBin==5) //60-80
     w=(0.719912-0.167292*x+0.017196*x*x+0.000861*x*x*x)/(1.-0.336146*x+0.037731*x*x); 
*/
  return TMath::Max(0.,w) ;  
}

  

