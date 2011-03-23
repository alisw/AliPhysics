#include "TChain.h"
#include "TTree.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH3F.h"
#include "TParticle.h"
#include "TCanvas.h"
#include "TStyle.h"

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskPi0.h"
#include "AliCaloPhoton.h"
#include "AliPHOSGeometry.h"
#include "AliESDEvent.h"
#include "AliESDCaloCells.h"
#include "AliESDCaloCluster.h"
#include "AliESDVertex.h"
#include "AliESDtrackCuts.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliTriggerAnalysis.h"

// Analysis task to fill histograms with PHOS ESD clusters and cells
// Authors: Yuri Kharlov
// Date   : 28.05.2009

ClassImp(AliAnalysisTaskPi0)

//________________________________________________________________________
AliAnalysisTaskPi0::AliAnalysisTaskPi0(const char *name) 
: AliAnalysisTaskSE(name),
  fESDtrackCuts(0),
  fOutputContainer(0),
  fPHOSEvent(0),
  fnCINT1B(0),
  fnCINT1A(0),
  fnCINT1C(0),
  fnCINT1E(0),
  fPHOSGeo(0),
  fEventCounter(0),
  fTriggerAnalysis(new AliTriggerAnalysis)
{
  // Constructor
  Int_t nBin=10 ;
  for(Int_t i=0;i<nBin;i++){
    for(Int_t j=0;j<2;j++)
      fPHOSEvents[i][j]=0 ;
  }
  
  // Output slots #0 write into a TH1 container
  DefineOutput(1,TList::Class());

  // Set bad channel map
  char key[55] ;
  for(Int_t i=0; i<6; i++){
    sprintf(key,"PHOS_BadMap_mod%d",i) ;
    fPHOSBadMap[i]=new TH2I(key,"Bad Modules map",64,0.,64.,56,0.,56.) ;
  }
  // Initialize the PHOS geometry
  fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP") ;

}

//________________________________________________________________________
void AliAnalysisTaskPi0::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  // ESD histograms
  if(fOutputContainer != NULL){
    delete fOutputContainer;
  }
  fOutputContainer = new TList();
  fOutputContainer->SetOwner(kTRUE);

  fOutputContainer->Add(new TH1I("hCellMultEvent"  ,"PHOS cell multiplicity per event"    ,2000,0,2000));
  fOutputContainer->Add(new TH1I("hCellMultEventM1","PHOS cell multiplicity per event, M1",2000,0,2000));
  fOutputContainer->Add(new TH1I("hCellMultEventM2","PHOS cell multiplicity per event, M2",2000,0,2000));
  fOutputContainer->Add(new TH1I("hCellMultEventM3","PHOS cell multiplicity per event, M3",2000,0,2000));
  fOutputContainer->Add(new TH1I("hClusterMult"      ,"CaloCluster multiplicity"     ,100,0,100));
  fOutputContainer->Add(new TH1I("hPHOSClusterMult"  ,"PHOS cluster multiplicity"    ,100,0,100));
  fOutputContainer->Add(new TH1I("hPHOSClusterMultM1","PHOS cluster multiplicity, M1",100,0,100));
  fOutputContainer->Add(new TH1I("hPHOSClusterMultM2","PHOS cluster multiplicity, M2",100,0,100));
  fOutputContainer->Add(new TH1I("hPHOSClusterMultM3","PHOS cluster multiplicity, M3",100,0,100));
  fOutputContainer->Add(new TH1F("hCellEnergy"  ,"Cell energy"            ,3000,0.,30.));
  fOutputContainer->Add(new TH1F("hCellEnergyM1","Cell energy in module 1",3000,0.,30.));
  fOutputContainer->Add(new TH1F("hCellEnergyM2","Cell energy in module 2",3000,0.,30.));
  fOutputContainer->Add(new TH1F("hCellEnergyM3","Cell energy in module 3",3000,0.,30.));
  fOutputContainer->Add(new TH1F("hClusterEnergy"  ,"Cluster energy"      ,3000,0.,30.));
  fOutputContainer->Add(new TH1F("hClusterEnergyM1","Cluster energy, M1"  ,3000,0.,30.));
  fOutputContainer->Add(new TH1F("hClusterEnergyM2","Cluster energy, M2"  ,3000,0.,30.));
  fOutputContainer->Add(new TH1F("hClusterEnergyM3","Cluster energy, M3"  ,3000,0.,30.));
  fOutputContainer->Add(new TH2F("hClusterEvsN"  ,"Cluster energy vs digit multiplicity"    ,3000,0.,30.,40,0.,40.));
  fOutputContainer->Add(new TH2F("hClusterEvsNM1","Cluster energy vs digit multiplicity, M1",3000,0.,30.,40,0.,40.));
  fOutputContainer->Add(new TH2F("hClusterEvsNM2","Cluster energy vs digit multiplicity, M2",3000,0.,30.,40,0.,40.));
  fOutputContainer->Add(new TH2F("hClusterEvsNM3","Cluster energy vs digit multiplicity, M3",3000,0.,30.,40,0.,40.));
  fOutputContainer->Add(new TH1I("hCellMultClu"  ,"Cell multiplicity per cluster"    ,200,0,200));
  fOutputContainer->Add(new TH1I("hCellMultCluM1","Cell multiplicity per cluster, M1",200,0,200));
  fOutputContainer->Add(new TH1I("hCellMultCluM2","Cell multiplicity per cluster, M3",200,0,200));
  fOutputContainer->Add(new TH1I("hCellMultCluM3","Cell multiplicity per cluster, M3",200,0,200));
  fOutputContainer->Add(new TH1I("hModule","Module events",5,0.,5.));
  fOutputContainer->Add(new TH1F("hSelEvents","Selected events",7,0.5,7.5));

  fOutputContainer->Add(new TH2F("hCellNXZM1","Cell (X,Z), M1" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCellNXZM2","Cell (X,Z), M2" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCellNXZM3","Cell (X,Z), M3" ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCellEXZM1","Cell E(X,Z), M1",64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCellEXZM2","Cell E(X,Z), M2",64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCellEXZM3","Cell E(X,Z), M3",64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluNXZM1","Clu (X,Z), M1"   ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluNXZM2","Clu (X,Z), M2"   ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluNXZM3","Clu (X,Z), M3"   ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluEXZM1","Clu E(X,Z), M1"  ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluEXZM2","Clu E(X,Z), M2"  ,64,0.5,64.5, 56,0.5,56.5));
  fOutputContainer->Add(new TH2F("hCluEXZM3","Clu E(X,Z), M3"  ,64,0.5,64.5, 56,0.5,56.5));

  Int_t nM       = 750;
  Double_t mMin  = 0.0;
  Double_t mMax  = 1.5;
  Int_t nPt      = 400;
  Double_t ptMin = 0;
  Double_t ptMax = 40;
  fOutputContainer->Add(new TH2F("hAsymPtPi0","(A,p_{T})_{#gamma#gamma} #pi^{0}"     ,50,0.,1.,    nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hAsymPtEta","(A,p_{T})_{#gamma#gamma} #eta"        ,50,0.,1.,    nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMassPtA10" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0"   ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtA08" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.8"   ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtA07" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7"   ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtA01" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.1"   ,nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMassPtA10nvtx" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, no vtx" ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtA07nvtx" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7, no vtx" ,nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMassPtA10vtx" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, vtx"     ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtA07vtx" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7, vtx"     ,nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMassPtA10V0AND" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, V0AND",nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtA07V0AND" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7, V0AND",nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMassPtA10V0PU" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, pileup",nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtA07V0PU" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7, pileup",nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMassPtvA10","(M,p_{T})_{#gamma#gamma}, 0<A<1.0, ESD vertex",nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtvA07","(M,p_{T})_{#gamma#gamma}, 0<A<0.7, ESD vertex",nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH3F("hMassPtCA10","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMassPtCA10_cpv","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMassPtCA10_disp","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMassPtCA10_both","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMassPtCA07","(M,p_{T},C)_{#gamma#gamma}, 0<A<0.7" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMassPtCA07_cpv","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMassPtCA07_disp","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMassPtCA07_both","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));

  fOutputContainer->Add(new TH2F("hMassSingle_all","(M,p_{T})_{#gamma#gamma}, no PID" ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassSingle_cpv","(M,p_{T})_{#gamma#gamma}, no PID" ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassSingle_disp","(M,p_{T})_{#gamma#gamma}, no PID" ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassSingle_both","(M,p_{T})_{#gamma#gamma}, no PID" ,nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMassPtM1","(M,p_{T})_{#gamma#gamma}, module 1"    ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtM2","(M,p_{T})_{#gamma#gamma}, module 2"    ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtM3","(M,p_{T})_{#gamma#gamma}, module 3"    ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtM12","(M,p_{T})_{#gamma#gamma}, modules 1,2",nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtM13","(M,p_{T})_{#gamma#gamma}, modules 1,3",nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtM23","(M,p_{T})_{#gamma#gamma}, modules 2,3",nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMassPtN3","(M,p_{T})_{#gamma#gamma}, N_{cell}>2"  ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtN4","(M,p_{T})_{#gamma#gamma}, N_{cell}>3"  ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtN5","(M,p_{T})_{#gamma#gamma}, N_{cell}>4"  ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMassPtN6","(M,p_{T})_{#gamma#gamma}, N_{cell}>5"  ,nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMiAsymPt","(A,p_{T})_{#gamma#gamma}"                ,50,0.,1.,    nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMiMassPtA10" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0"   ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtA08" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.8"   ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtA07" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7"   ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtA01" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.1"   ,nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMiMassPtA10nvtx" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, no vtx" ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtA07nvtx" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7, no vtx" ,nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMiMassPtA10vtx" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, vtx"     ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtA07vtx" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7, vtx"     ,nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMiMassPtA10V0AND" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, V0AND",nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtA07V0AND" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7, V0AND",nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMiMassPtA10V0PU" ,"(M,p_{T})_{#gamma#gamma}, 0<A<1.0, pileup",nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtA07V0PU" ,"(M,p_{T})_{#gamma#gamma}, 0<A<0.7, pileup",nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMiMassPtvA10","(M,p_{T})_{#gamma#gamma}, 0<A<1.0, ESD vertex",nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtvA07","(M,p_{T})_{#gamma#gamma}, 0<A<0.7, ESD vertex",nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH3F("hMiMassPtCA10","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMiMassPtCA10_cpv","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMiMassPtCA10_disp","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMiMassPtCA10_both","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMiMassPtCA07","(M,p_{T},C)_{#gamma#gamma}, 0<A<0.7" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMiMassPtCA07_cpv","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMiMassPtCA07_disp","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));
  fOutputContainer->Add(new TH3F("hMiMassPtCA07_both","(M,p_{T},C)_{#gamma#gamma}, 0<A<1.0" ,nM,mMin,mMax,nPt,ptMin,ptMax,8,0.,8.));

  fOutputContainer->Add(new TH2F("hMiMassSingle_all","(M,p_{T})_{#gamma#gamma}, no PID" ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassSingle_cpv","(M,p_{T})_{#gamma#gamma}, no PID" ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassSingle_disp","(M,p_{T})_{#gamma#gamma}, no PID" ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassSingle_both","(M,p_{T})_{#gamma#gamma}, no PID" ,nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMiMassPtM1","(M,p_{T})_{#gamma#gamma}, module 1"    ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtM2","(M,p_{T})_{#gamma#gamma}, module 2"    ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtM3","(M,p_{T})_{#gamma#gamma}, module 3"    ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtM12","(M,p_{T})_{#gamma#gamma}, modules 1,2",nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtM13","(M,p_{T})_{#gamma#gamma}, modules 1,3",nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtM23","(M,p_{T})_{#gamma#gamma}, modules 2,3",nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH2F("hMiMassPtN3","(M,p_{T})_{#gamma#gamma}, N_{cell}>2"  ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtN4","(M,p_{T})_{#gamma#gamma}, N_{cell}>3"  ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtN5","(M,p_{T})_{#gamma#gamma}, N_{cell}>4"  ,nM,mMin,mMax,nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH2F("hMiMassPtN6","(M,p_{T})_{#gamma#gamma}, N_{cell}>5"  ,nM,mMin,mMax,nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH1F("hPhotonKappa","#kappa(#gamma)",200,0.,20.));
  fOutputContainer->Add(new TH1F("hPhotonPt","p_{T}(#gamma)",200,0.,20.));
  fOutputContainer->Add(new TH1F("hPhotonPx","p_{x}(#gamma)",200,0.,20.));
  fOutputContainer->Add(new TH1F("hPhotonPy","p_{y}(#gamma)",200,0.,20.));

  fOutputContainer->Add(new TH1F("hTrigClass","Trigger class",5,0.5,5.5));

  fOutputContainer->Add(new TH1F("hZvertex","Z vertex",200,-50.,+50.));
  fOutputContainer->Add(new TH1F("hNvertexTracks","N of primary tracks from the primary vertex",150,0.,150.));
  fOutputContainer->Add(new TH1F("hTrackMult","Charged track multiplicity",150,0.,150.));

  // Create ESD track cut

  fESDtrackCuts = AliESDtrackCuts::GetStandardITSTPCTrackCuts2010();
  fESDtrackCuts ->SetMaxDCAToVertexZ(2);
  fESDtrackCuts ->SetEtaRange(-0.8, 0.8);
  fESDtrackCuts ->SetPtRange(0.15);

  PostData(1, fOutputContainer);

}

//________________________________________________________________________
void AliAnalysisTaskPi0::UserExec(Option_t *) 
{
  // Main loop, called for each event
  // Analyze ESD

  // Event selection flags

  Bool_t eventVtxExist    = kFALSE;
  Bool_t eventVtxZ10cm    = kFALSE;
  Bool_t eventPileup      = kFALSE;
  Bool_t eventV0AND       = kFALSE;

  AliESDEvent *event = dynamic_cast<AliESDEvent*>(InputEvent());
  Int_t eventNumberInFile = event->GetEventNumberInFile();
  if (!event) {
     Printf("ERROR: Could not retrieve event");
     return;
  }

  if(fPHOSEvent)
    fPHOSEvent->Clear() ;
  else
    fPHOSEvent = new TClonesArray("AliCaloPhoton",50) ;

  // Checks if we have a primary vertex
  // Get primary vertices form ESD

  if      (event->GetPrimaryVertexTracks()->GetNContributors()>0)
    eventVtxExist    = kTRUE;
  else if (event->GetPrimaryVertexSPD()   ->GetNContributors()>0)
    eventVtxExist    = kTRUE;

  const AliESDVertex *esdVertex5 = event->GetPrimaryVertex();

  Double_t vtx0[3] = {0,0,0}; // don't rely on ESD vertex, assume (0,0,0)
  Double_t vtx5[3];
  vtx5[0] = esdVertex5->GetX();
  vtx5[1] = esdVertex5->GetY();
  vtx5[2] = esdVertex5->GetZ();

  FillHistogram("hNvertexTracks",esdVertex5->GetNContributors());
  FillHistogram("hZvertex"      ,esdVertex5->GetZ());
  if (TMath::Abs(esdVertex5->GetZ()) < 10. )
    eventVtxZ10cm = kTRUE;

  if (event->IsPileupFromSPD())
    eventPileup = kTRUE;

  eventV0AND = fTriggerAnalysis->IsOfflineTriggerFired(event, AliTriggerAnalysis::kV0AND);

  // Fill event statistics for different selection criteria

  FillHistogram("hSelEvents",1) ;
  if (eventVtxExist)
    FillHistogram("hSelEvents",2) ;
  if (eventVtxExist && eventVtxZ10cm)
    FillHistogram("hSelEvents",3) ;
  if (eventVtxExist && eventVtxZ10cm && eventV0AND)
    FillHistogram("hSelEvents",4) ;
  if (eventVtxExist && eventVtxZ10cm && eventV0AND && eventPileup)
    FillHistogram("hSelEvents",5) ;
  if (eventPileup)
    FillHistogram("hSelEvents",6) ;
  if(eventV0AND){
    FillHistogram("hSelEvents",7) ;
  }
      
  //Vtx class z-bin
  Int_t zvtx = (Int_t)((vtx5[2]+10.)/2.) ;
  if(zvtx<0)zvtx=0 ;
  if(zvtx>9)zvtx=9 ;

  TString trigClasses = event->GetFiredTriggerClasses();

  if (trigClasses.Contains("CINT1B")) fnCINT1B++;
  if (trigClasses.Contains("CINT1A")) fnCINT1A++;
  if (trigClasses.Contains("CINT1C")) fnCINT1C++;
  if (trigClasses.Contains("CINT1-E")) fnCINT1E++;

  //Calculate charged multiplicity
  Int_t trackMult = 0;
  for (Int_t i=0;i<event->GetNumberOfTracks();++i) {
    AliESDtrack *track = new AliESDtrack(*event->GetTrack(i)) ;
    if(fESDtrackCuts->AcceptTrack(track) &&  TMath::Abs(track->Eta())< 0.8)
      trackMult++;
    delete track;
  }
  FillHistogram("hTrackMult",trackMult+0.5) ;

  Int_t centr=0 ;
  //always zero centrality
  if(!fPHOSEvents[zvtx][centr]) fPHOSEvents[zvtx][centr]=new TList() ;
  TList * prevPHOS = fPHOSEvents[zvtx][centr] ;

  if(trackMult<=2)
    centr=0 ;
  else 
    if(trackMult<=5)
      centr=1 ;
    else
      if(trackMult<=9)
        centr=2 ;
      else
        if(trackMult<=14)
          centr=3 ;
        else
          if(trackMult<=22)
            centr=4 ;
          else
            if(trackMult<=35)
              centr=5 ;
            else
              if(trackMult<=50)
                centr=6 ;
              else
                centr=7 ;


  AliESDCaloCluster *clu1;
  TLorentzVector p1,p2,p12, pv1,pv2,pv12;
  AliESDCaloCells *cells      = event->GetPHOSCells();

  Int_t multClust = event->GetNumberOfCaloClusters();
  Int_t multCells = cells->GetNumberOfCells();
  FillHistogram("hClusterMult",multClust);
  FillHistogram("hCellMultEvent",multCells);

  Printf("Event %d, trig.class %s, period %d, bc %d, orbit %d",
     eventNumberInFile,trigClasses.Data(),event->GetPeriodNumber(),
     event->GetBunchCrossNumber(),event->GetOrbitNumber());
  Printf("\tthere are %d caloclusters and %d calocells",
     multClust,multCells);

  // Get PHOS rotation matrices from ESD and set them to the PHOS geometry
  if(fEventCounter == 0) {
    for(Int_t mod=0; mod<5; mod++) {
      if(!event->GetPHOSMatrix(mod)) continue;
      fPHOSGeo->SetMisalMatrix(event->GetPHOSMatrix(mod),mod) ;
      Printf("PHOS geo matrix %p for module # %d is set\n", event->GetPHOSMatrix(mod), mod);
    }
  }

  Float_t  energy;
  Int_t    mod1, relId[4], cellAbsId, cellX, cellZ;

  // Single loop over cells

  Int_t nCellModule[3] = {0,0,0};
  for (Int_t iCell=0; iCell<multCells; iCell++) {
    cellAbsId = cells->GetCellNumber(iCell);
    fPHOSGeo->AbsToRelNumbering(cellAbsId,relId);
    mod1  = relId[0];
    cellX = relId[2];
    cellZ = relId[3] ;
    energy = cells->GetAmplitude(iCell);
    FillHistogram("hCellEnergy",energy);
    if      (mod1==1) {
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

  // Single loop over clusters fills cluster histograms

  Int_t    digMult;
  Int_t    multPHOSClust[4]  = {0,0,0,0};
  Float_t  position[3];

  for (Int_t i1=0; i1<multClust; i1++) {
    clu1 = event->GetCaloCluster(i1);
    if ( !clu1->IsPHOS() ) continue;

    digMult = clu1->GetNCells();
    clu1->GetPosition(position);
    TVector3 global1(position) ;
    fPHOSGeo->GlobalPos2RelId(global1,relId) ;
    mod1  = relId[0] ;
    cellX = relId[2];
    cellZ = relId[3] ;
    if ( !IsGoodChannel("PHOS",mod1,cellX,cellZ) ) continue ;
    
    cellAbsId = clu1->GetCellAbsId(0);
    fPHOSGeo->AbsToRelNumbering(cellAbsId,relId);
    mod1   = relId[0];
    energy = clu1->E();
    
    multPHOSClust[0]++;
    FillHistogram("hClusterEnergy",energy);
    FillHistogram("hClusterEvsN",energy,digMult);
    FillHistogram("hCellMultClu",digMult);
    if      (mod1==1) {
      multPHOSClust[1]++;
      FillHistogram("hClusterEvsNM1",energy,digMult);
      FillHistogram("hCellMultCluM1",digMult);
      FillHistogram("hClusterEnergyM1",energy);
      if (digMult == 1 && energy > 3.) {
    FillHistogram("hCluNXZM1",cellX,cellZ,1.);
    FillHistogram("hCluEXZM1",cellX,cellZ,energy);
      }
    }
    else if (mod1==2) {
      multPHOSClust[2]++;
      FillHistogram("hClusterEvsNM2",energy,digMult);
      FillHistogram("hCellMultCluM2",digMult);
      FillHistogram("hClusterEnergyM2",energy);
      if (digMult == 1 && energy > 3.) {
    FillHistogram("hCluNXZM2",cellX,cellZ,1.);
    FillHistogram("hCluEXZM2",cellX,cellZ,energy);
      }
    }
    else if (mod1==3) {
      multPHOSClust[3]++;
      FillHistogram("hClusterEvsNM3",energy,digMult);
      FillHistogram("hCellMultCluM3",digMult);
      FillHistogram("hClusterEnergyM3",energy);
      if (digMult == 1 && energy > 3.) {
    FillHistogram("hCluNXZM3",cellX,cellZ,1.);
    FillHistogram("hCluEXZM3",cellX,cellZ,energy);
      }
    }
    
    if (digMult > 2) {
      clu1 ->GetMomentum(p1 ,vtx0);
      Double_t pAbs = p1.P();
      Double_t pT   = p1.Pt();
      Double_t pX   = p1.Px();
      Double_t pY   = p1.Py();
      Double_t kappa = pAbs - TMath::Power(0.135,2)/4./pAbs;
      
      FillHistogram("hPhotonKappa",kappa);
      FillHistogram("hPhotonPt",pT);
      FillHistogram("hPhotonPx",pX);
      FillHistogram("hPhotonPy",pY);
    }
  }
  FillHistogram("hPHOSClusterMult",multPHOSClust[0]);
  FillHistogram("hPHOSClusterMultM1",multPHOSClust[1]);
  FillHistogram("hPHOSClusterMultM2",multPHOSClust[2]);
  FillHistogram("hPHOSClusterMultM3",multPHOSClust[3]);

  //Select photons for inv mass calculation
  Int_t inPHOS=0 ;
  for (Int_t i1=0; i1<multClust; i1++) {
    clu1 = event->GetCaloCluster(i1);
    if ( !clu1->IsPHOS() || clu1->E()<0.3) continue;

    clu1->GetPosition(position);
    TVector3 global1(position) ;
    fPHOSGeo->GlobalPos2RelId(global1,relId) ;
    mod1  = relId[0] ;
    cellX = relId[2];
    cellZ = relId[3] ;
    if ( !IsGoodChannel("PHOS",mod1,cellX,cellZ) ) continue ;

    clu1 ->GetMomentum(p1 ,vtx0);
    clu1 ->GetMomentum(pv1,vtx5);
    digMult   = clu1->GetNCells();
    new((*fPHOSEvent)[inPHOS]) AliCaloPhoton(p1.X(),p1.Py(),p1.Z(),p1.E()) ;
    AliCaloPhoton * ph = (AliCaloPhoton*)fPHOSEvent->At(inPHOS) ;
    ph->SetModule(mod1) ;
    ph->SetMomV2(&pv1) ;
    ph->SetNCells(clu1->GetNCells());
    ph->SetDispBit(TestLambda(clu1->GetM20(),clu1->GetM02())) ;
    ph->SetCPVBit(clu1->GetEmcCpvDistance()>10.) ;

    inPHOS++ ;
  }

  // Fill Real disribution

  for (Int_t i1=0; i1<inPHOS-1; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1) ;
    for (Int_t i2=i1+1; i2<inPHOS; i2++) {
      AliCaloPhoton * ph2=(AliCaloPhoton*)fPHOSEvent->At(i2) ;
      p12  = *ph1  + *ph2;
      pv12 = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
      Double_t asym  = TMath::Abs((ph1->Energy()-ph2->Energy())/(ph1->Energy()+ph2->Energy()));

      if (ph1->GetNCells()>2 && ph2->GetNCells()>2) {
        FillHistogram("hMassPtA10",p12.M() ,p12.Pt() );
        FillHistogram("hMassPtvA10",pv12.M(),pv12.Pt());
        FillHistogram("hMassPtCA10",p12.M() ,p12.Pt(), centr+0.5);
        FillHistogram("hMassSingle_all",p12.M(),ph1->Pt()) ;
        FillHistogram("hMassSingle_all",p12.M(),ph2->Pt()) ;

	if(!eventVtxExist)
	  FillHistogram("hMassPtA10nvtx",p12.M() ,p12.Pt() );
	if(eventVtxExist)
	  FillHistogram("hMassPtA10vtx"  ,p12.M() ,p12.Pt() );
	if(eventVtxExist && eventV0AND)
	  FillHistogram("hMassPtA10V0AND",p12.M() ,p12.Pt() );
	if(eventPileup)
	  FillHistogram("hMassPtA10PU"   ,p12.M() ,p12.Pt() );

        if(ph1->IsCPVOK())
          FillHistogram("hMassSingle_cpv",p12.M(),ph1->Pt()) ;
        if(ph2->IsCPVOK())
          FillHistogram("hMassSingle_cpv",p12.M(),ph2->Pt()) ;
        if(ph1->IsDispOK())
          FillHistogram("hMassSingle_disp",p12.M(),ph1->Pt()) ;
        if(ph2->IsDispOK())
          FillHistogram("hMassSingle_disp",p12.M(),ph2->Pt()) ;
        if(ph1->IsCPVOK() && ph1->IsDispOK())
          FillHistogram("hMassSingle_both",p12.M(),ph1->Pt()) ;
        if(ph2->IsCPVOK() && ph2->IsDispOK())
          FillHistogram("hMassSingle_both",p12.M(),ph2->Pt()) ;
 

        if(ph1->IsCPVOK() && ph2->IsCPVOK())
          FillHistogram("hMassPtCA10_cpv",p12.M() ,p12.Pt(), centr+0.5);
        if(ph1->IsDispOK() && ph2->IsDispOK()){
          FillHistogram("hMassPtCA10_disp",p12.M() ,p12.Pt(), centr+0.5);
          if(ph1->IsCPVOK() && ph2->IsCPVOK())
            FillHistogram("hMassPtCA10_both",p12.M() ,p12.Pt(), centr+0.5);
        }
        if (asym<0.8) {
          FillHistogram("hMassPtA08",p12.M(),p12.Pt());
        }
        if (asym<0.7) {
          FillHistogram("hMassPtA07",p12.M(),p12.Pt());
	  FillHistogram("hMassPtvA07",pv12.M(),pv12.Pt());
          FillHistogram("hMassPtCA07",p12.M() ,p12.Pt(), centr+0.5);
	  if(!eventVtxExist)
	    FillHistogram("hMassPtA07nvtx",p12.M() ,p12.Pt() );
	  if(eventVtxExist)
	    FillHistogram("hMassPtA07vtx"  ,p12.M() ,p12.Pt() );
	  if(eventVtxExist && eventV0AND)
	    FillHistogram("hMassPtA07V0AND",p12.M() ,p12.Pt() );
	  if(eventPileup)
	    FillHistogram("hMassPtA07PU"   ,p12.M() ,p12.Pt() );
          if(ph1->IsCPVOK() && ph2->IsCPVOK())
            FillHistogram("hMassPtCA07_cpv",p12.M() ,p12.Pt(), centr+0.5);
          if(ph1->IsDispOK() && ph2->IsDispOK()){
            FillHistogram("hMassPtCA07_disp",p12.M() ,p12.Pt(), centr+0.5);
            if(ph1->IsCPVOK() && ph2->IsCPVOK())
              FillHistogram("hMassPtCA07_both",p12.M() ,p12.Pt(), centr+0.5);
        }

        }
        if (asym<0.1) {
          FillHistogram("hMassPtA01",p12.M(),p12.Pt());
        }
	if (TMath::Abs(p12.M()-0.135)<0.03)
	  FillHistogram("hAsymPtPi0",asym   ,p12.Pt());
	if (TMath::Abs(p12.M()-0.547)<0.09)
	  FillHistogram("hAsymPtEta",asym   ,p12.Pt());

        if (ph1->Module()==1 && ph2->Module()==1) FillHistogram("hMassPtM1",p12.M() ,p12.Pt() );
        if (ph1->Module()==2 && ph2->Module()==2) FillHistogram("hMassPtM2",p12.M() ,p12.Pt() );
        if (ph1->Module()==3 && ph2->Module()==3) FillHistogram("hMassPtM3",p12.M() ,p12.Pt() );
        if (ph1->Module()!=3 && ph2->Module()!=3) FillHistogram("hMassPtM12",p12.M() ,p12.Pt() );
        if (ph1->Module()!=2 && ph2->Module()!=2) FillHistogram("hMassPtM13",p12.M() ,p12.Pt() );
        if (ph1->Module()!=1 && ph2->Module()!=1) FillHistogram("hMassPtM23",p12.M() ,p12.Pt() );

      }

      if (ph1->GetNCells()>3 && ph2->GetNCells()>3) {
        FillHistogram("hMassPtN3",p12.M() ,p12.Pt() );
      }
      if (ph1->GetNCells()>4 && ph2->GetNCells()>4) {
        FillHistogram("hMassPtN4",p12.M() ,p12.Pt() );
      }
      if (ph1->GetNCells()>5 && ph2->GetNCells()>5) {
        FillHistogram("hMassPtN5",p12.M() ,p12.Pt() );
      }
      if (ph1->GetNCells()>6 && ph2->GetNCells()>6) {
        FillHistogram("hMassPtN6",p12.M() ,p12.Pt() );
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
      p12  = *ph1  + *ph2;
      pv12 = *(ph1->GetMomV2()) + *(ph2->GetMomV2());
      Double_t asym  = TMath::Abs((ph1->Energy()-ph2->Energy())/(ph1->Energy()+ph2->Energy()));

      if (ph1->GetNCells()>2 && ph2->GetNCells()>2) {
        FillHistogram("hMiMassPtA10",p12.M() ,p12.Pt() );
        FillHistogram("hMiMassPtvA10",pv12.M(),pv12.Pt());
        FillHistogram("hMiMassPtCA10",p12.M() ,p12.Pt(), centr+0.5);
        FillHistogram("hMiMassSingle_all",p12.M(),ph1->Pt()) ;
        FillHistogram("hMiMassSingle_all",p12.M(),ph2->Pt()) ;

	if(!eventVtxExist)
	  FillHistogram("hMiMassPtA10nvtx",p12.M() ,p12.Pt() );
	if(eventVtxExist)
	  FillHistogram("hMiMassPtA10vtx"  ,p12.M() ,p12.Pt() );
	if(eventVtxExist && eventV0AND)
	  FillHistogram("hMiMassPtA10V0AND",p12.M() ,p12.Pt() );
	if(eventPileup)
	  FillHistogram("hMiMassPtA10PU"   ,p12.M() ,p12.Pt() );

        if(ph1->IsCPVOK())
          FillHistogram("hMiMassSingle_cpv",p12.M(),ph1->Pt()) ;
        if(ph2->IsCPVOK())
          FillHistogram("hMiMassSingle_cpv",p12.M(),ph2->Pt()) ;
        if(ph1->IsDispOK())
          FillHistogram("hMiMassSingle_disp",p12.M(),ph1->Pt()) ;
        if(ph2->IsDispOK())
          FillHistogram("hMiMassSingle_disp",p12.M(),ph2->Pt()) ;
        if(ph1->IsCPVOK() && ph1->IsDispOK())
          FillHistogram("hMiMassSingle_both",p12.M(),ph1->Pt()) ;
        if(ph2->IsCPVOK() && ph2->IsDispOK())
          FillHistogram("hMiMassSingle_both",p12.M(),ph2->Pt()) ;


        if(ph1->IsCPVOK() && ph2->IsCPVOK())
          FillHistogram("hMiMassPtCA10_cpv",p12.M() ,p12.Pt(), centr+0.5);
        if(ph1->IsDispOK() && ph2->IsDispOK()){
          FillHistogram("hMiMassPtCA10_disp",p12.M() ,p12.Pt(), centr+0.5);
          if(ph1->IsCPVOK() && ph2->IsCPVOK())
            FillHistogram("hMiMassPtCA10_both",p12.M() ,p12.Pt(), centr+0.5);
        }
        if (asym<0.8) {
          FillHistogram("hMiMassPtA08",p12.M(),p12.Pt());
        }
        if (asym<0.7) {
          FillHistogram("hMiMassPtA07",p12.M(),p12.Pt());
 	  FillHistogram("hMiMassPtvA07",pv12.M(),pv12.Pt());
	  FillHistogram("hMiMassPtCA07",p12.M() ,p12.Pt(), centr+0.5);
	  if(!eventVtxExist)
	    FillHistogram("hMiMassPtA07nvtx",p12.M() ,p12.Pt() );
	  if(eventVtxExist)
	    FillHistogram("hMiMassPtA07vtx"  ,p12.M() ,p12.Pt() );
	  if(eventVtxExist && eventV0AND)
	    FillHistogram("hMiMassPtA07V0AND",p12.M() ,p12.Pt() );
	  if(eventPileup)
	    FillHistogram("hMiMassPtA07PU"   ,p12.M() ,p12.Pt() );
          if(ph1->IsCPVOK() && ph2->IsCPVOK())
            FillHistogram("hMiMassPtCA07_cpv",p12.M() ,p12.Pt(), centr+0.5);
          if(ph1->IsDispOK() && ph2->IsDispOK()){
            FillHistogram("hMiMassPtCA07_disp",p12.M() ,p12.Pt(), centr+0.5);
            if(ph1->IsCPVOK() && ph2->IsCPVOK())
              FillHistogram("hMiMassPtCA07_both",p12.M() ,p12.Pt(), centr+0.5);
          }
        }
        if (asym<0.1) {
          FillHistogram("hMiMassPtA01",p12.M(),p12.Pt());
        }
        FillHistogram("hMiAsymPt",asym   ,p12.Pt());

        if (ph1->Module()==1 && ph2->Module()==1) FillHistogram("hMiMassPtM1",p12.M() ,p12.Pt() );
        if (ph1->Module()==2 && ph2->Module()==2) FillHistogram("hMiMassPtM2",p12.M() ,p12.Pt() );
        if (ph1->Module()==3 && ph2->Module()==3) FillHistogram("hMiMassPtM3",p12.M() ,p12.Pt() );
        if (ph1->Module()!=3 && ph2->Module()!=3) FillHistogram("hMiMassPtM12",p12.M() ,p12.Pt() );
        if (ph1->Module()!=2 && ph2->Module()!=2) FillHistogram("hMiMassPtM13",p12.M() ,p12.Pt() );
        if (ph1->Module()!=1 && ph2->Module()!=1) FillHistogram("hMiMassPtM23",p12.M() ,p12.Pt() );

      }

      if (ph1->GetNCells()>3 && ph2->GetNCells()>3) {
        FillHistogram("hMiMassPtN3",p12.M() ,p12.Pt() );
      }
      if (ph1->GetNCells()>4 && ph2->GetNCells()>4) {
        FillHistogram("hMiMassPtN4",p12.M() ,p12.Pt() );
      }
      if (ph1->GetNCells()>5 && ph2->GetNCells()>5) {
        FillHistogram("hMiMassPtN5",p12.M() ,p12.Pt() );
      }
      if (ph1->GetNCells()>6 && ph2->GetNCells()>6) {
        FillHistogram("fMihMassPtN6",p12.M() ,p12.Pt() );
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
void AliAnalysisTaskPi0::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query

  Int_t nPP = fnCINT1B - fnCINT1A - fnCINT1C + 2*fnCINT1E;
  FillHistogram("hTrigClass",1,fnCINT1B);
  FillHistogram("hTrigClass",2,fnCINT1A);
  FillHistogram("hTrigClass",3,fnCINT1C);
  FillHistogram("hTrigClass",4,fnCINT1E);
  FillHistogram("hTrigClass",5,nPP);
  Printf("fnCINT1B=%d, fnCINT1A=%d ,fnCINT1C=%d ,fnCINT1E=%d, nPP=%d",
     fnCINT1B,fnCINT1A,fnCINT1C,fnCINT1E,nPP);
   
}

//________________________________________________________________________
Bool_t AliAnalysisTaskPi0::IsGoodChannel(const char * det, Int_t mod, Int_t ix, Int_t iz)
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
void AliAnalysisTaskPi0::FillHistogram(const char * key,Double_t x)const{
  //FillHistogram
  TH1I * tmpI = dynamic_cast<TH1I*>(fOutputContainer->FindObject(key)) ;
  if(tmpI){
    tmpI->Fill(x) ;
    return ;
  }
  TH1F * tmpF = dynamic_cast<TH1F*>(fOutputContainer->FindObject(key)) ;
  if(tmpF){
    tmpF->Fill(x) ;
    return ;
  }
  TH1D * tmpD = dynamic_cast<TH1D*>(fOutputContainer->FindObject(key)) ;
  if(tmpD){
    tmpD->Fill(x) ;
    return ;
  }
  AliInfo(Form("can not find histogram <%s> ",key)) ;
}
//_____________________________________________________________________________
void AliAnalysisTaskPi0::FillHistogram(const char * key,Double_t x,Double_t y)const{
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
void AliAnalysisTaskPi0::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z) const{
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
Bool_t AliAnalysisTaskPi0::TestLambda(Double_t l1,Double_t l2){
  Double_t l1Mean=1.22 ;
  Double_t l2Mean=2.0 ;
  Double_t l1Sigma=0.42 ;
  Double_t l2Sigma=0.71 ;
  Double_t c=-0.59 ;
  Double_t R2=(l1-l1Mean)*(l1-l1Mean)/l1Sigma/l1Sigma+(l2-l2Mean)*(l2-l2Mean)/l2Sigma/l2Sigma-c*(l1-l1Mean)*(l2-l2Mean)/l1Sigma/l2Sigma ;
  return (R2<9.) ;

}
