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

ClassImp(AliAnalysisTaskPi0DiffEfficiency)

//________________________________________________________________________
AliAnalysisTaskPi0DiffEfficiency::AliAnalysisTaskPi0DiffEfficiency(const char *name) 
: AliAnalysisTaskPi0Efficiency(name),
  fPHOSEvent1(0),
  fPHOSEvent2(0)
{
  


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
  fOutputContainer = new TList();
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
    snprintf(key,55,"hPhotDisp_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hPhotDisp2_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hPhotDispwou_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hPhotBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hPhotBothcore_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;

    snprintf(key,55,"hNegPhotAll_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegPhotAllcore_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegPhotAllwou_cen%d",cent) ;
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
    snprintf(key,55,"hNegPhotDispwou_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegPhotBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegPhotBothcore_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    
    snprintf(key,55,"hOldMassPtAll_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hOldMassPtAllcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hOldMassPtAllwou_cen%d",cent) ;
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
    snprintf(key,55,"hOldMassPtDispwou_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hOldMassPtBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hOldMassPtBothcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    
    snprintf(key,55,"hNewMassPtAll_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNewMassPtAllcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNewMassPtAllwou_cen%d",cent) ;
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
    snprintf(key,55,"hNewMassPtDispwou_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNewMassPtBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNewMassPtBothcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;

    snprintf(key,55,"hNegMassPtAll_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegMassPtAllcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegMassPtAllwou_cen%d",cent) ;
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
    snprintf(key,55,"hNegMassPtDispwou_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegMassPtBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNegMassPtBothcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    
    
    snprintf(key,55,"hMassPtAll_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMassPtAllcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMassPtAllwou_cen%d",cent) ;
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
    snprintf(key,55,"hMassPtDisp2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMassPtDispwou_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMassPtBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMassPtBothcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    
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
    snprintf(key,55,"hMiMassPtAllcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtAllwou_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtCPV_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtCPVcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtCPV2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtDisp_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtDisp2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtDispwou_cen%d",cent) ;
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
    
    snprintf(key,55,"hMCMassPtAll_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtAllcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtAllwou_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtCPV_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtCPVcore_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtCPV2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtDisp_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtDisp2_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtDispwou_cen%d",cent) ;
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
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMCPhotAllcore_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMCPhotAllwou_cen%d",cent) ;
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
    snprintf(key,55,"hMCPhotDisp2_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMCPhotDispwou_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMCPhotBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;    
    snprintf(key,55,"hMCPhotBothcore_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;    
    
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
    
    if(inPHOSold>=fPHOSEvent1->GetSize()){
      fPHOSEvent1->Expand(inPHOSold+50) ;
    }
    AliCaloPhoton * ph = new((*fPHOSEvent1)[inPHOSold]) AliCaloPhoton(pv1.X(),pv1.Py(),pv1.Z(),pv1.E()) ;
    ph->SetModule(mod) ;
    AliPHOSAodCluster cluPHOS1(*clu);
    cluPHOS1.Recalibrate(fPHOSCalibData,cellsOld); // modify the cell energies
    Double_t ecore=CoreEnergy(&cluPHOS1) ;
    pv1*= ecore/pv1.E() ;    
    ph->SetMomV2(&pv1) ;
    ph->SetNCells(clu->GetNCells());
    ph->SetDispBit(TestLambda(clu->E(),clu->GetM20(),clu->GetM02())) ;
    ph->SetDisp2Bit(TestLambda2(clu->E(),clu->GetM20(),clu->GetM02())) ;
    ph->SetCPVBit(clu->GetEmcCpvDistance()>2.) ;
    ph->SetCPV2Bit(clu->GetEmcCpvDistance()>4.) ;
    if(!survive) //this cluster found in list after embedding, skipping it
      ph->SetTagged(1) ;
    ph->SetPhoton(clu->GetNExMax()<2); // Remember, if it is unfolded

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
    
    if(inPHOSemb>=fPHOSEvent2->GetSize()){
      fPHOSEvent2->Expand(inPHOSemb+50) ;
    }
    AliCaloPhoton * ph = new((*fPHOSEvent2)[inPHOSemb]) AliCaloPhoton(pv1.X(),pv1.Py(),pv1.Z(),pv1.E()) ;
    ph->SetModule(mod) ;
    AliPHOSAodCluster cluPHOS1(*clu);
    cluPHOS1.Recalibrate(fPHOSCalibData,cellsEmb); // modify the cell energies
    Double_t ecore=CoreEnergy(&cluPHOS1) ;
    pv1*= ecore/pv1.E() ;    
    ph->SetMomV2(&pv1) ;
    ph->SetNCells(clu->GetNCells());
    ph->SetDispBit(TestLambda(clu->E(),clu->GetM20(),clu->GetM02())) ;
    ph->SetDisp2Bit(TestLambda2(clu->E(),clu->GetM20(),clu->GetM02())) ;
    ph->SetCPVBit(clu->GetEmcCpvDistance()>2.) ;
    ph->SetCPV2Bit(clu->GetEmcCpvDistance()>4.) ;
    if(!survive) //this cluster found in list after embedding, skipping it
      ph->SetTagged(1) ;
    ph->SetPhoton(clu->GetNExMax()<2); // Remember, if it is unfolded

    inPHOSemb++ ;
  }
  
  //Single photon
  for (Int_t i1=0; i1<inPHOSold; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent1->At(i1) ;
    if(!ph1->IsTagged())
      continue ;
    snprintf(key,55,"hPhotAll_cen%d",fCenBin) ;
    FillHistogram(key,ph1->Pt(),-1.) ;
    snprintf(key,55,"hNegPhotAll_cen%d",fCenBin) ;
    FillHistogram(key,ph1->Pt(),-1.) ;
    snprintf(key,55,"hPhotAllcore_cen%d",fCenBin) ;
    FillHistogram(key,ph1->GetMomV2()->Pt(),-1.) ;
    snprintf(key,55,"hNegPhotAllcore_cen%d",fCenBin) ;
    FillHistogram(key,ph1->GetMomV2()->Pt(),-1.) ;
    if(ph1->IsPhoton()){
      snprintf(key,55,"hPhotAllwou_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt(),-1.) ;
      snprintf(key,55,"hNegPhotAllwou_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt(),-1.) ;
    }    
    if(ph1->IsCPVOK() ){
      snprintf(key,55,"hPhotCPV_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt(),-1.) ;
      snprintf(key,55,"hNegPhotCPV_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt(),-1.) ;
      snprintf(key,55,"hPhotCPVcore_cen%d",fCenBin) ;
      FillHistogram(key,ph1->GetMomV2()->Pt(),-1.) ;
      snprintf(key,55,"hNegPhotCPVcore_cen%d",fCenBin) ;
      FillHistogram(key,ph1->GetMomV2()->Pt(),-1.) ;
    }
    if(ph1->IsCPV2OK() ){
      snprintf(key,55,"hPhotCPV2_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt(),-1.) ;
      snprintf(key,55,"hNegPhotCPV2_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt(),-1.) ;
    }
    if(ph1->IsDisp2OK()){
      snprintf(key,55,"hPhotDisp2_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt(),-1.) ;
      snprintf(key,55,"hNegPhotDisp2_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt(),-1.) ;
    }
    if(ph1->IsDispOK()){
      snprintf(key,55,"hPhotDisp_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt(),-1.) ;
      snprintf(key,55,"hNegPhotDisp_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt(),-1.) ;
      if(ph1->IsPhoton()){
        snprintf(key,55,"hPhotDispwou_cen%d",fCenBin) ;
        FillHistogram(key,ph1->Pt(),-1.) ;
        snprintf(key,55,"hNegPhotDispwou_cen%d",fCenBin) ;
        FillHistogram(key,ph1->Pt(),-1.) ;
      }
      if(ph1->IsCPVOK()){
	snprintf(key,55,"hPhotBoth_cen%d",fCenBin) ;
        FillHistogram(key,ph1->Pt(),-1.) ;
	snprintf(key,55,"hNegPhotBoth_cen%d",fCenBin) ;
        FillHistogram(key,ph1->Pt(),-1.) ;
	snprintf(key,55,"hPhotBothcore_cen%d",fCenBin) ;
        FillHistogram(key,ph1->GetMomV2()->Pt(),-1.) ;
	snprintf(key,55,"hNegPhotBothcore_cen%d",fCenBin) ;
        FillHistogram(key,ph1->GetMomV2()->Pt(),-1.) ;
      }
    } // end of loop i2
  } // end of loop i1 

  for (Int_t i1=0; i1<inPHOSemb; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent2->At(i1) ;
    if(!ph1->IsTagged())
      continue ;
    snprintf(key,55,"hPhotAll_cen%d",fCenBin) ;
    FillHistogram(key,ph1->Pt(),1.) ;
    snprintf(key,55,"hPhotAllcore_cen%d",fCenBin) ;
    FillHistogram(key,ph1->GetMomV2()->Pt(),1.) ;
    if(ph1->IsPhoton()){
      snprintf(key,55,"hPhotAllwou_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt(),1.) ;
    }
    if(ph1->IsCPVOK() ){
      snprintf(key,55,"hPhotCPV_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt(),1.) ;
      snprintf(key,55,"hPhotCPVcore_cen%d",fCenBin) ;
      FillHistogram(key,ph1->GetMomV2()->Pt(),1.) ;
    }
    if(ph1->IsCPV2OK() ){
      snprintf(key,55,"hPhotCPV2_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt(),1.) ;
    }
    if(ph1->IsDisp2OK()){
      snprintf(key,55,"hPhotDisp2_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt(),1.) ;
    }
    if(ph1->IsDispOK()){
      snprintf(key,55,"hPhotDisp_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt(),1.) ;
      if(ph1->IsPhoton()){
        snprintf(key,55,"hPhotDispwou_cen%d",fCenBin) ;
        FillHistogram(key,ph1->Pt(),1.) ;
      }
      if(ph1->IsCPVOK()){
	snprintf(key,55,"hPhotBoth_cen%d",fCenBin) ;
        FillHistogram(key,ph1->Pt(),1.) ;
	snprintf(key,55,"hPhotBothcore_cen%d",fCenBin) ;
        FillHistogram(key,ph1->GetMomV2()->Pt(),1.) ;
      }
    } // end of loop i2
  } // end of loop i1 



  // Fill Real disribution:
  // Disappeared clusters enter with negative contribution
  // In addition fill control histogam with Real before embedding
  for (Int_t i1=0; i1<inPHOSold-1; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent1->At(i1) ;
    for (Int_t i2=i1+1; i2<inPHOSold; i2++) {
      AliCaloPhoton * ph2=(AliCaloPhoton*)fPHOSEvent1->At(i2) ;

      TLorentzVector p12  = *ph1  + *ph2;
      TLorentzVector pv12 = *(ph1->GetMomV2()) + *(ph2->GetMomV2());      
      Double_t a=TMath::Abs((ph1->E()-ph2->E())/(ph1->E()+ph2->E())) ;

      //Fill Controll histogram: Real before embedding
      snprintf(key,55,"hOldMassPtAll_cen%d",fCenBin) ;
      FillHistogram(key,p12.M() ,p12.Pt(),-1.) ;
      snprintf(key,55,"hOldMassPtAllcore_cen%d",fCenBin) ;
      FillHistogram(key,pv12.M() ,pv12.Pt(),-1.) ;
      if(ph1->IsPhoton() && ph2->IsPhoton()){
        snprintf(key,55,"hOldMassPtAllwou_cen%d",fCenBin) ;
        FillHistogram(key,p12.M() ,p12.Pt(),-1.) ;
      }
      if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	snprintf(key,55,"hOldMassPtCPV_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt(),-1) ;
	snprintf(key,55,"hOldMassPtCPV_cen%d",fCenBin) ;
	FillHistogram(key,pv12.M(), pv12.Pt(),-1) ;
      }
      if(ph1->IsCPV2OK() && ph2->IsCPV2OK()){
	snprintf(key,55,"hOldMassPtCPV2_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt(),-1) ;
      }
      if(ph1->IsDisp2OK() && ph2->IsDisp2OK()){
	snprintf(key,55,"hOldMassPtDisp2_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt(),-1) ;
      }
      if(ph1->IsDispOK() && ph2->IsDispOK()){
	snprintf(key,55,"hOldMassPtDisp_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt(),-1) ;
        if(ph1->IsPhoton() && ph2->IsPhoton()){
	  snprintf(key,55,"hOldMassPtDispwou_cen%d",fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt(),-1) ;
	}
	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  snprintf(key,55,"hOldMassPtBoth_cen%d",fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt(),-1) ;
	  snprintf(key,55,"hOldMassPtBothcore_cen%d",fCenBin) ;
	  FillHistogram(key,pv12.M() ,pv12.Pt(),-1) ;
        }
      }
      
      //Now fill main histograms with negative contributions
      if(!(ph1->IsTagged() || ph2->IsTagged()) )
        continue ;
      snprintf(key,55,"hMassPtAll_cen%d",fCenBin) ;
      FillHistogram(key,p12.M() ,p12.Pt(),-1.) ;
      snprintf(key,55,"hMassPtAllcore_cen%d",fCenBin) ;
      FillHistogram(key,pv12.M() ,pv12.Pt(),-1.) ;
      if(ph1->IsPhoton() && ph2->IsPhoton()){
        snprintf(key,55,"hMassPtAllwou_cen%d",fCenBin) ;
        FillHistogram(key,p12.M() ,p12.Pt(),-1.) ;
      }
      if(a<0.9){
        snprintf(key,55,"hMassPtAll_a09_cen%d",fCenBin) ;
        FillHistogram(key,p12.M() ,p12.Pt(),-1.) ;
        if(a<0.8){
          snprintf(key,55,"hMassPtAll_a08_cen%d",fCenBin) ;
          FillHistogram(key,p12.M() ,p12.Pt(),-1.) ;
          if(a<0.7){
            snprintf(key,55,"hMassPtAll_a07_cen%d",fCenBin) ;
            FillHistogram(key,p12.M() ,p12.Pt(),-1.) ;
          }
        }
      }
      snprintf(key,55,"hNegMassPtAll_cen%d",fCenBin) ;
      FillHistogram(key,p12.M() ,p12.Pt(),-1.) ;
      snprintf(key,55,"hNegMassPtAllcore_cen%d",fCenBin) ;
      FillHistogram(key,pv12.M() ,pv12.Pt(),-1.) ;
      if(ph1->IsPhoton() && ph2->IsPhoton()){
        snprintf(key,55,"hNegMassPtAllwou_cen%d",fCenBin) ;
        FillHistogram(key,p12.M() ,p12.Pt(),-1.) ;
      }
      
      if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	snprintf(key,55,"hMassPtCPV_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt(),-1) ;
	snprintf(key,55,"hMassPtCPVcore_cen%d",fCenBin) ;
	FillHistogram(key,pv12.M() ,pv12.Pt(),-1) ;
        if(a<0.9){
          snprintf(key,55,"hMassPtCPV_a09_cen%d",fCenBin) ;
          FillHistogram(key,p12.M() ,p12.Pt(),-1.) ;
          if(a<0.8){
            snprintf(key,55,"hMassPtCPV_a08_cen%d",fCenBin) ;
            FillHistogram(key,p12.M() ,p12.Pt(),-1.) ;
            if(a<0.7){
              snprintf(key,55,"hMassPtCPV_a07_cen%d",fCenBin) ;
              FillHistogram(key,p12.M() ,p12.Pt(),-1.) ;
            }
          }
        }
	snprintf(key,55,"hNegMassPtCPV_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt(),-1) ;
	snprintf(key,55,"hNegMassPtCPVcore_cen%d",fCenBin) ;
	FillHistogram(key,pv12.M() ,pv12.Pt(),-1) ;
      }
      if(ph1->IsCPV2OK() && ph2->IsCPV2OK()){
	snprintf(key,55,"hMassPtCPV2_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt(),-1) ;
        if(a<0.9){
          snprintf(key,55,"hMassPtCPV2_a09_cen%d",fCenBin) ;
          FillHistogram(key,p12.M() ,p12.Pt(),-1.) ;
          if(a<0.8){
            snprintf(key,55,"hMassPtCPV2_a08_cen%d",fCenBin) ;
            FillHistogram(key,p12.M() ,p12.Pt(),-1.) ;
            if(a<0.7){
              snprintf(key,55,"hMassPtCPV2_a07_cen%d",fCenBin) ;
              FillHistogram(key,p12.M() ,p12.Pt(),-1.) ;
            }
          }
	}
        snprintf(key,55,"hNegMassPtCPV2_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt(),-1) ;
      }
      
      if(ph1->IsDisp2OK() && ph2->IsDisp2OK()){
	snprintf(key,55,"hMassPtDisp2_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt(),-1) ;
      }
      if(ph1->IsDispOK() && ph2->IsDispOK()){
	snprintf(key,55,"hMassPtDisp_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt(),-1) ;
        if(ph1->IsPhoton() && ph2->IsPhoton()){
  	  snprintf(key,55,"hMassPtDispwou_cen%d",fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt(),-1) ;
	}
        if(a<0.9){
          snprintf(key,55,"hMassPtDisp_a09_cen%d",fCenBin) ;
          FillHistogram(key,p12.M() ,p12.Pt(),-1.) ;
          if(a<0.8){
            snprintf(key,55,"hMassPtDisp_a08_cen%d",fCenBin) ;
            FillHistogram(key,p12.M() ,p12.Pt()) ;
            if(a<0.7){
              snprintf(key,55,"hMassPtDisp_a07_cen%d",fCenBin) ;
              FillHistogram(key,p12.M() ,p12.Pt(),-1.) ;
            }
          }
        }
	snprintf(key,55,"hNegMassPtDisp_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt(),-1) ;
        if(ph1->IsPhoton() && ph2->IsPhoton()){
  	  snprintf(key,55,"hNegMassPtDispwou_cen%d",fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt(),-1) ;
	}
	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  snprintf(key,55,"hMassPtBoth_cen%d",fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt(),-1) ;
	  snprintf(key,55,"hMassPtBothcore_cen%d",fCenBin) ;
	  FillHistogram(key,pv12.M() ,pv12.Pt(),-1) ;
          if(a<0.9){
            snprintf(key,55,"hMassPtBoth_a09_cen%d",fCenBin) ;
            FillHistogram(key,p12.M() ,p12.Pt(),-1.) ;
            if(a<0.8){
              snprintf(key,55,"hMassPtBoth_a08_cen%d",fCenBin) ;
              FillHistogram(key,p12.M() ,p12.Pt(),-1.) ;
              if(a<0.7){
                snprintf(key,55,"hMassPtBoth_a07_cen%d",fCenBin) ;
                FillHistogram(key,p12.M() ,p12.Pt(),-1.) ;
             }
            }
          }
	  snprintf(key,55,"hNegMassPtBoth_cen%d",fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt(),-1) ;
	  snprintf(key,55,"hNegMassPtBothcore_cen%d",fCenBin) ;
	  FillHistogram(key,pv12.M() ,pv12.Pt(),-1) ;
        }
      }
    } // end of loop i2
  } // end of loop i1 


  // Further fill Real disribution
  // now with positive contribution from new clusters
  // ass well fill controll histogram
  for (Int_t i1=0; i1<inPHOSemb-1; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent2->At(i1) ;
    for (Int_t i2=i1+1; i2<inPHOSemb; i2++) {
      AliCaloPhoton * ph2=(AliCaloPhoton*)fPHOSEvent2->At(i2) ;

      TLorentzVector p12  = *ph1  + *ph2;
      TLorentzVector pv12 = *(ph1->GetMomV2()) + *(ph2->GetMomV2());      
      Double_t a=TMath::Abs((ph1->E()-ph2->E())/(ph1->E()+ph2->E())) ;

      // Controll histogram: Real after embedding
      snprintf(key,55,"hNewMassPtAll_cen%d",fCenBin) ;
      FillHistogram(key,p12.M() ,p12.Pt(),1.) ;
      snprintf(key,55,"hNewMassPtAllcore_cen%d",fCenBin) ;
      FillHistogram(key,pv12.M() ,pv12.Pt(),1.) ;
      if(ph1->IsPhoton() && ph2->IsPhoton()){
        snprintf(key,55,"hNewMassPtAllwou_cen%d",fCenBin) ;
        FillHistogram(key,p12.M() ,p12.Pt(),1.) ;
      }
      if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	snprintf(key,55,"hNewMassPtCPV_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt(),1) ;
	snprintf(key,55,"hNewMassPtCPVcore_cen%d",fCenBin) ;
	FillHistogram(key,pv12.M() ,pv12.Pt(),1) ;
      }
      if(ph1->IsCPV2OK() && ph2->IsCPV2OK()){
	snprintf(key,55,"hNewMassPtCPV2_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt(),1) ;
      }
      if(ph1->IsDisp2OK() && ph2->IsDisp2OK()){
	snprintf(key,55,"hNewMassPtDisp2_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt(),1) ;
      }
      if(ph1->IsDispOK() && ph2->IsDispOK()){
	snprintf(key,55,"hNewMassPtDisp_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt(),1) ;
        if(ph1->IsPhoton() && ph2->IsPhoton()){
	  snprintf(key,55,"hNewMassPtDispwou_cen%d",fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt(),1) ;
	}
	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  snprintf(key,55,"hNewMassPtBoth_cen%d",fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt(),1) ;
	  snprintf(key,55,"hNewMassPtBothcore_cen%d",fCenBin) ;
	  FillHistogram(key,pv12.M() ,pv12.Pt(),1) ;
        }
      }
     
      //Now fill main histogamm
      //new clusters with positive contribution
      if(!(ph1->IsTagged() || ph2->IsTagged()) )
        continue ;
      snprintf(key,55,"hMassPtAll_cen%d",fCenBin) ;
      FillHistogram(key,p12.M() ,p12.Pt(),1.) ;
      snprintf(key,55,"hMassPtAllcore_cen%d",fCenBin) ;
      FillHistogram(key,pv12.M() ,pv12.Pt(),1.) ;
      if(ph1->IsPhoton() && ph2->IsPhoton()){
        snprintf(key,55,"hMassPtAllwou_cen%d",fCenBin) ;
        FillHistogram(key,p12.M() ,p12.Pt(),1.) ;
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
	FillHistogram(key,p12.M() ,p12.Pt(),1) ;
	snprintf(key,55,"hMassPtCPVcore_cen%d",fCenBin) ;
	FillHistogram(key,pv12.M() ,pv12.Pt(),1) ;
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
	FillHistogram(key,p12.M() ,p12.Pt(),1) ;
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
	FillHistogram(key,p12.M() ,p12.Pt(),1) ;
      }
      if(ph1->IsDispOK() && ph2->IsDispOK()){
	snprintf(key,55,"hMassPtDisp_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt(),1) ;
        if(ph1->IsPhoton() && ph2->IsPhoton()){
	  snprintf(key,55,"hMassPtDispwou_cen%d",fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt(),1) ;
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
	  FillHistogram(key,p12.M() ,p12.Pt(),1) ;
	  snprintf(key,55,"hMassPtBothcore_cen%d",fCenBin) ;
	  FillHistogram(key,pv12.M() ,pv12.Pt(),1) ;
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


  //now mixed, does not really matter old or new list of clusters
  for (Int_t i1=0; i1<inPHOSemb; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent2->At(i1) ;
    for(Int_t ev=0; ev<prevPHOS->GetSize();ev++){
      TClonesArray * mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev)) ;
      for(Int_t i2=0; i2<mixPHOS->GetEntriesFast();i2++){
	AliCaloPhoton * ph2=(AliCaloPhoton*)mixPHOS->At(i2) ;
	TLorentzVector p12  = *ph1  + *ph2;
	TLorentzVector pv12 = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
        Double_t a=TMath::Abs((ph1->E()-ph2->E())/(ph1->E()+ph2->E())) ;
	
	snprintf(key,55,"hMiMassPtAll_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt(),1.) ;
	snprintf(key,55,"hMiMassPtAllcore_cen%d",fCenBin) ;
	FillHistogram(key,pv12.M() ,pv12.Pt(),1.) ;
        if(ph1->IsPhoton() && ph2->IsPhoton()){
  	  snprintf(key,55,"hMiMassPtAllwou_cen%d",fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt(),1.) ;
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
	  FillHistogram(key,p12.M() ,p12.Pt(),1.) ;
	  snprintf(key,55,"hMiMassPtCPVcore_cen%d",fCenBin) ;
	  FillHistogram(key,pv12.M() ,pv12.Pt(),1.) ;
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
	  FillHistogram(key,p12.M() ,p12.Pt(),1.) ;
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
	  FillHistogram(key,p12.M() ,p12.Pt(),1.) ;
	}
	if(ph1->IsDispOK() && ph2->IsDispOK()){
	  snprintf(key,55,"hMiMassPtDisp_cen%d",fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt(),1.) ;
          if(ph1->IsPhoton() && ph2->IsPhoton()){
	    snprintf(key,55,"hMiMassPtDispwou_cen%d",fCenBin) ;
	    FillHistogram(key,p12.M() ,p12.Pt(),1.) ;
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
	    FillHistogram(key,p12.M() ,p12.Pt(),1.) ;
	    snprintf(key,55,"hMiMassPtBothcore_cen%d",fCenBin) ;
	    FillHistogram(key,pv12.M() ,pv12.Pt(),1.) ;
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
  if(fPHOSEvent2->GetEntriesFast()>0){
    prevPHOS->AddFirst(fPHOSEvent2) ;
    fPHOSEvent2=0;
    delete fPHOSEvent1;
    fPHOSEvent1=0;
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



