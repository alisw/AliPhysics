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
#include "AliEventplane.h"
#include "TProfile.h"
#include "AliOADBContainer.h"
#include "AliMagF.h"
#include "AliAODMCParticle.h"
#include "AliEPFlattener.h"

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
  fPHOSEvent(0x0),
  fV0AFlat(0x0),
  fV0CFlat(0x0),
  fRunNumber(0),
  fCentrality(0.),
  fCenBin(0),
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
  
  // Output slots #0 write into a TH1 container
  DefineOutput(1,TList::Class());

  // Initialize the PHOS geometry
  fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP") ;

  //We have to apply re-calibration for pass1 LCH10h
  // Initialize decalibration factors in the form of the OCDB object


}

//________________________________________________________________________
void AliAnalysisTaskgg::UserCreateOutputObjects()
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
 
  fOutputContainer->Add(new TH2F("phiRP","Event plane", 100,0.,TMath::Pi(),100,0.,100.)) ;
  fOutputContainer->Add(new TH2F("phiRPflat","Event plane", 100,0.,TMath::Pi(),100,0.,100.)) ;
 
 
  //vertex distribution
  fOutputContainer->Add(new TH2F("hZvertex","Z vertex position", 50,-25.,25.,nRuns,0.,float(nRuns))) ;
  
  //Centrality
  fOutputContainer->Add(new TH2F("hCentrality","Event centrality", 100,0.,100.,nRuns,0.,float(nRuns))) ;
  fOutputContainer->Add(new TH2F("hCenPHOS","Centrality vs PHOSclusters", 100,0.,100.,200,0.,200.)) ;
  fOutputContainer->Add(new TH2F("hCenPHOSCells","Centrality vs PHOS cells", 100,0.,100.,100,0.,1000.)) ;
  fOutputContainer->Add(new TH2F("hCenTrack","Centrality vs tracks", 100,0.,100.,100,0.,15000.)) ;  
  fOutputContainer->Add(new TH2F("hCluEvsClu_All","ClusterMult vs E",20,0.,2.,50,0.,50.)) ;
  fOutputContainer->Add(new TH2F("hCluEvsClu_CPV","ClusterMult vs E",20,0.,2.,50,0.,50.)) ;
  fOutputContainer->Add(new TH2F("hCluEvsClu_Disp","ClusterMult vs E",20,0.,2.,50,0.,50.)) ;
  fOutputContainer->Add(new TH2F("hCluEvsClu_Both","ClusterMult vs E",20,0.,2.,50,0.,50.)) ;
  fOutputContainer->Add(new TH2F("hCluEvsCluM","ClusterMult vs E",200,0.,20.,100,0.,20.)) ;
  fOutputContainer->Add(new TH2F("hCenTOF","Centrality vs PHOS TOF", 100,0.,100.,600,-6.e-6,6.e-6)) ;
    
  //PHOS QA
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
 			
  fOutputContainer->Add(new TH2F("hCPVr","CPV radius",100,0.,20.,100,0.,2.));
  fOutputContainer->Add(new TH3F("hLambda","Lambdas for all clusters",150,0.,30.,150,0.,30.,200,0.,2.));
  
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

  const Int_t nCuts=4 ;
  char cut[7][20] ;
  sprintf(cut[0],"All") ;
  sprintf(cut[1],"Disp") ;
  sprintf(cut[2],"CPV") ;
  sprintf(cut[3],"Both") ;
  sprintf(cut[4],"Dist1") ;
  sprintf(cut[5],"Dist2") ;
  sprintf(cut[6],"Dist3") ;
  
    
  for(Int_t iCut=0; iCut<nCuts; iCut++){
    for(Int_t ikT=0; ikT<6; ikT++){ 
//      fOutputContainer->Add(new TH3F(Form("hOSLPF_%s_%s",cut[iCut],kTbins[ikT]),"Out-Side-Long, Pair Frame",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));
 
      fOutputContainer->Add(new TH3F(Form("hOSLCMS_%s_%s",cut[iCut],kTbins[ikT]),"Out-Side-Long, CMS",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));
//      fOutputContainer->Add(new TH3F(Form("hYKPPF_%s_%s",cut[iCut],kTbins[ikT]),"YKP, Pair Frame",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));
//      fOutputContainer->Add(new TH3F(Form("hYKPCMS_%s_%s",cut[iCut],kTbins[ikT]),"YKP, CMS",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));

      fOutputContainer->Add(new TH2F(Form("hetaphi2D_%s_%s",cut[iCut],kTbins[ikT]),"Eta-phi-E correlations",50,-0.25,0.25,100,-TMath::Pi()/6.,TMath::Pi()/6.));
      fOutputContainer->Add(new TH3F(Form("hetaphi_%s_%s",cut[iCut],kTbins[ikT]),"Eta-phi-E correlations",50,-0.25,0.25,100,-TMath::Pi()/6.,TMath::Pi()/6.,20,-0.2,0.2));
      fOutputContainer->Add(new TH2F(Form("hetaphiRP_%s_%s",cut[iCut],kTbins[ikT]),"Eta-phi-E correlations",10,0.,TMath::Pi(),100,-TMath::Pi()/6.,TMath::Pi()/6.));
//      fOutputContainer->Add(new TH2F(Form("hdXdZ_%s_%s",cut[iCut],kTbins[ikT]),"dXdZ",200,-200,200,200,-200.,200.));
      
      
//      fOutputContainer->Add(new TH3F(Form("hMiOSLPF_%s_%s",cut[iCut],kTbins[ikT]),"Out-Side-Long, Pair Frame",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));

      fOutputContainer->Add(new TH3F(Form("hMiOSLCMS_%s_%s",cut[iCut],kTbins[ikT]),"Out-Side-Long, CMS",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));
      fOutputContainer->Add(new TH3F(Form("hMi2OSLCMS_%s_%s",cut[iCut],kTbins[ikT]),"Out-Side-Long, CMS",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));
//      fOutputContainer->Add(new TH3F(Form("hMiYKPPF_%s_%s",cut[iCut],kTbins[ikT]),"YKP, Pair Frame",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));
//      fOutputContainer->Add(new TH3F(Form("hMiYKPCMS_%s_%s",cut[iCut],kTbins[ikT]),"YKP, CMS",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));

      fOutputContainer->Add(new TH2F(Form("hMietaphi2D_%s_%s",cut[iCut],kTbins[ikT]),"Eta-phi-E correlations",50,-0.25,0.25,100,-TMath::Pi()/6.,TMath::Pi()/6.));
      fOutputContainer->Add(new TH3F(Form("hMietaphi_%s_%s",cut[iCut],kTbins[ikT]),"Eta-phi-E correlations",50,-0.25,0.25,100,-TMath::Pi()/6.,TMath::Pi()/6.,20,-0.2,0.2));
      fOutputContainer->Add(new TH2F(Form("hMietaphiRP_%s_%s",cut[iCut],kTbins[ikT]),"Eta-phi-E correlations",10,0.,TMath::Pi(),100,-TMath::Pi()/6.,TMath::Pi()/6.));
      fOutputContainer->Add(new TH3F(Form("hMi2etaphi_%s_%s",cut[iCut],kTbins[ikT]),"Eta-phi-E correlations",50,-0.25,0.25,100,-TMath::Pi()/6.,TMath::Pi()/6.,20,-0.2,0.2));
//      fOutputContainer->Add(new TH2F(Form("hMidXdZ_%s_%s",cut[iCut],kTbins[ikT]),"dXdZ",200,-200,200,200,-200.,200.));
    
    }        

    fOutputContainer->Add(new TH2F(Form("hQinv_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,100,0.,10.));
    fOutputContainer->Add(new TH2F(Form("hMiQinv_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,100,0.,10.));
    fOutputContainer->Add(new TH2F(Form("hMi2Qinv_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,100,0.,10.));
    fOutputContainer->Add(new TH2F(Form("hQinvCut_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,100,0.,10.));
    fOutputContainer->Add(new TH2F(Form("hMiQinvCut_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,100,0.,10.));
    fOutputContainer->Add(new TH2F(Form("hMi2QinvCut_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,100,0.,10.));
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

  fRunNumber=ConvertRunNumber(fEvent->GetRunNumber()) ;
  FillHistogram("hSelEvents",1.5,fRunNumber-0.5) ;
  FillHistogram("hTotSelEvents",1.5) ;
  
 
  if(fEventCounter == 0 && fIsPbPb) {
    //Get Event Plane flattening
    Int_t run = fEvent->GetRunNumber() ;
    AliOADBContainer flatContainer("phosFlat");
    flatContainer.InitFromFile("$ALICE_PHYSICS/OADB/PHOS/PHOSflat.root","phosFlat");
    TObjArray *arr = (TObjArray*)flatContainer.GetObject(run,"phosFlat");
    if(!arr){
      AliError(Form("Can not read Flattening for run %d. \n From file $ALICE_PHYSICS/OADB/PHOS/PHOSflat.root",run)) ;    
      arr = (TObjArray*)flatContainer.GetObject(1,"phosFlat"); //default
    }
        
    AliInfo(Form("Setting PHOS flattening with name %s \n",arr->GetName())) ;
//    AliEPFlattener * h = (AliEPFlattener*)arr->At(0) ;  
//      if(fTPCFlat) delete fTPCFlat ;
//      fTPCFlat = new AliEPFlattener() ;
//      fTPCFlat = h ;
    AliEPFlattener * h = (AliEPFlattener*)arr->At(1) ;  
    if(fV0AFlat) delete fV0AFlat ;
    fV0AFlat = new AliEPFlattener() ;
    fV0AFlat = h ;
    h = (AliEPFlattener*)arr->At(2) ;  
    if(fV0CFlat) delete fV0CFlat ;
    fV0CFlat = new AliEPFlattener() ;
    fV0CFlat = h ;
   
    fEventCounter++ ;
  }

  
  // Checks if we have a primary vertex
  // Get primary vertices form AOD
  const AliAODVertex *esdVertex5 = fEvent->GetPrimaryVertex();

  Double_t vtx5[3] ={0.,0.,0.};
  
  vtx5[0] = esdVertex5->GetX();
  vtx5[1] = esdVertex5->GetY();
  vtx5[2] = esdVertex5->GetZ();
  
  
  FillHistogram("hZvertex",esdVertex5->GetZ(),fRunNumber-0.5);
  if (TMath::Abs(esdVertex5->GetZ()) > 10. ){
    PostData(1, fOutputContainer);
    return;
  }
  FillHistogram("hSelEvents",2.5,fRunNumber-0.5) ;
  FillHistogram("hTotSelEvents",2.5) ;

/*  
  if(fEvent->IsPileupFromSPD()){
    PostData(1, fOutputContainer);
    return;
  } 
*/
  FillHistogram("hSelEvents",3.5,fRunNumber-0.5) ;
  FillHistogram("hTotSelEvents",3.5) ;  
  
  
  
  //Vtx class z-bin
  Int_t zvtx = (Int_t)((vtx5[2]+10.)/2.) ;
  if(zvtx<0)zvtx=0 ;
  if(zvtx>9)zvtx=9 ;
 
  AliCentrality *centrality = fEvent->GetCentrality(); 
  fCentrality=centrality->GetCentralityPercentile("V0M");

  if( fCentrality <= 0. || fCentrality>80. ){
    PostData(1, fOutputContainer);
    return;
  }

  FillHistogram("hSelEvents",4.5,fRunNumber-0.5) ;
  FillHistogram("hTotSelEvents",4.5) ;

  if(fCentrality<5.)
    fCenBin=0 ;
  else if(fCentrality<10.)
    fCenBin=1 ;
  else if(fCentrality<20.)
    fCenBin=2 ;
  else if(fCentrality<40.)
    fCenBin=3 ;
  else 
    fCenBin=4 ;


  //reaction plane
  Int_t irp=0 ;
    AliEventplane *eventPlane = fEvent->GetEventplane();
    if( ! eventPlane ) { //Event has no event plane
      PostData(1, fOutputContainer);
      return;
    }
    //V0A
    const Int_t harmonics = 2; 
    Double_t qx=0., qy=0.;  
    Double_t rpV0A = eventPlane->CalculateVZEROEventPlane(fEvent,8, harmonics,qx,qy);
    //V0C
    Double_t rpV0C = eventPlane->CalculateVZEROEventPlane(fEvent,9, harmonics,qx,qy);

    while(rpV0A<0)rpV0A+=TMath::TwoPi()/harmonics ;
    while(rpV0A>TMath::TwoPi()/harmonics)rpV0A-=TMath::TwoPi()/harmonics ;
    if(fIsPbPb)
      rpV0A = fV0AFlat->MakeFlat(rpV0A,fCentrality) ;
  
    while(rpV0C<0)rpV0C+=TMath::TwoPi()/harmonics ;
    while(rpV0C>TMath::TwoPi()/harmonics)rpV0C-=TMath::TwoPi()/harmonics ;
    if(fIsPbPb)
      rpV0C = fV0CFlat->MakeFlat(rpV0C,fCentrality) ;
  
    Double_t rpFull=0.5*(rpV0A+rpV0C) ;  
    FillHistogram("phiRPflat",rpFull,fCentrality) ;  

    //Reaction plane is defined in the range (0;pi)
    //We have 10 bins
    irp=Int_t(10.*(rpFull)/TMath::Pi());
    if(irp>9)irp=9 ;
    
  
 
  FillHistogram("hSelEvents",4.5,fRunNumber-0.5) ;
  FillHistogram("hTotSelEvents",5.5) ;
  //All event selections done
  FillHistogram("hCentrality",fCentrality,fRunNumber-0.5) ;
  

  if(!fPHOSEvents[zvtx][fCenBin][irp]) 
    fPHOSEvents[zvtx][fCenBin][irp]=new TList() ;
  TList * prevPHOS = fPHOSEvents[zvtx][fCenBin][irp] ;

  if(fPHOSEvent)
    fPHOSEvent->Clear() ;
  else
    fPHOSEvent = new TClonesArray("AliCaloPhoton",200) ;


  TVector3 vertex(vtx5);
  
  Int_t multClust = fEvent->GetNumberOfCaloClusters();
  Int_t inPHOS=0; 

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
    if ( !clu->IsPHOS() || clu->E()<0.1) continue;
    
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
    
    //Remove 6 noisy channels in run 139036
    if(fEvent->GetRunNumber()==139036 && mod==1 && 
       (cellX==9||cellX==10||cellX==11) && (cellZ==45 || cellZ==46))
      continue ;
    
    FillHistogram("hCluEvsCluM",clu->E(),clu->GetM02()) ;

    FillHistogram(Form("hTofM%d",mod),clu->E(),clu->GetTOF()) ;
    if(clu->E()>1.)
      FillHistogram("hCenTOF",fCentrality,clu->GetTOF()) ;
    if((clu->GetTOF()>1.5e-7) || (clu->GetTOF() <-2.5e-7) )
      continue ;
    
    
    if(clu->GetNCells()<3) continue ;
    if(clu->GetM02()<0.2)   continue ;    
    if(clu->GetMCEnergyFraction()>0.98) //Ecross cut, should be filled with Tender
     continue ;    
           
    TLorentzVector pv1 ;
    clu->GetMomentum(pv1 ,vtx5);
    
    FillHistogram(Form("hCluLowM%d",mod),cellX,cellZ,1.);
    if(clu->E()>1.5){
      FillHistogram(Form("hCluHighM%d",mod),cellX,cellZ,1.);
    }
 
    FillHistogram(Form("hLambda"),clu->GetM02(),clu->GetM20(),clu->E());
    FillHistogram("hCPVr",clu->Chi2(),clu->E());
    
 
    if(inPHOS>=fPHOSEvent->GetSize()){
      fPHOSEvent->Expand(inPHOS+50) ;
    }
    new((*fPHOSEvent)[inPHOS]) AliCaloPhoton(pv1.X(),pv1.Py(),pv1.Z(),pv1.E()) ;
    AliCaloPhoton * ph = (AliCaloPhoton*)fPHOSEvent->At(inPHOS) ;
    ph->SetModule(mod) ;
    pv1*= clu->GetCoreEnergy()/pv1.E() ;
    ph->SetMomV2(&pv1) ;
    ph->SetNCells(clu->GetNCells());    
    ph->SetDispBit(clu->Chi2()<2.5*2.5) ;
//  Cut on FullLamdas
    ph->SetDisp2Bit(clu->GetDispersion()<2.5*2.5) ;
    FillHistogram("hCluEvsClu_All",clu->E(),clu->GetNCells()) ;

//    Double_t distBC=clu->GetDistanceToBadChannel();
    if(ph->IsDispOK()){
      FillHistogram(Form("hCluDispM%d",mod),cellX,cellZ,1.);
      FillHistogram("hCluEvsClu_Disp",clu->E(),clu->GetNCells()) ;    
    }
    ph->SetCPVBit(clu->GetEmcCpvDistance()>2.5) ;
    if(ph->IsCPVOK()){
      FillHistogram(Form("hCluVetoM%d",mod),cellX,cellZ,1.);
      FillHistogram("hCluEvsClu_CPV",clu->E(),clu->GetNCells()) ;
      if(ph->IsDispOK()){
        FillHistogram("hCluEvsClu_Both",clu->E(),clu->GetNCells()) ;    
      }
    }
    
    ph->SetPrimary(clu->GetLabelAt(0)) ;
    ph->SetEMCx(position[0]) ;
    ph->SetEMCy(position[1]) ;
    ph->SetEMCz(position[2]) ;
    ph->SetLambdas(clu->GetM20(),clu->GetM02()) ;
    ph->SetUnfolded(clu->GetNExMax()<2); // Remember, if it is unfolded          
    inPHOS++ ;
  }
  
  FillHistogram("hCenPHOS",fCentrality,inPHOS) ;
  if(inPHOS==0){
    PostData(1, fOutputContainer);
    fEventCounter++;
    return ; 
  }
  
  Int_t jetStatus[5]; 
  for(Int_t mod=1; mod<5; mod++)
    jetStatus[mod]=JetRejection(mod) ;

  const Double_t kgMass=0. ;
	
  char cut[7][20] ;
  sprintf(cut[0],"All") ;
  sprintf(cut[1],"Disp") ;
  sprintf(cut[2],"CPV") ;
  sprintf(cut[3],"Both") ;
  sprintf(cut[4],"Dist1") ;
  sprintf(cut[5],"Dist2") ;
  sprintf(cut[6],"Dist3") ;

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
    
    for (Int_t i2=i1+1; i2<inPHOS; i2++) {
      AliCaloPhoton * ph2=(AliCaloPhoton*)fPHOSEvent->At(i2) ;
      //Cut on pair
// 	if(!PairCut(ph1,ph2,kDefault))
//	if(!PairCut(ph1,ph2,kBoth))
//	  continue;
     
//      if(!SecondaryPi0Cut(ph1,ph2))
//        continue ;

      AliFemtoTrack track2;
      AliFemtoThreeVector mom2;
      mom2.SetX(ph2->Px()) ;
      mom2.SetY(ph2->Py()) ;
      mom2.SetZ(ph2->Pz()) ;
      track2.SetP(mom2) ;
      AliFemtoParticle part2(&track2,kgMass) ;
      
      //Photons are sorted, try to remove it
      AliFemtoParticle *a = &part1 ;
      AliFemtoParticle *b = &part2 ;
      Double_t dEta = ph1->Eta()-ph2->Eta() ; 
      Double_t dPhi = ph1->Phi()-ph2->Phi() ; 
      Double_t dE   = ph1->E()  - ph2->E() ;
      Double_t dX = TMath::Power(ph1->EMCx() - ph2->EMCx(),2) + TMath::Power(ph1->EMCy() - ph2->EMCy(),2)  ;
      dX=TMath::Sign(TMath::Sqrt(dX),ph1->EMCx() - ph2->EMCx()) ;
      Double_t dZ = ph1->EMCz() - ph2->EMCz() ;
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
      TString kTbin="15" ;
      if(kT<0.2) kTbin="Kt00-02";
      else if(kT<0.4) kTbin="Kt02-04";
      else if(kT<0.7) kTbin="Kt04-07";
      else if(kT<1.) kTbin="Kt07-10";
      else if(kT<1.3) kTbin="Kt10-13";
      else if(kT<2.) kTbin="Kt13-20";
      else  continue;
      
      Double_t qs=pair.QSideCMS(), qo=pair.QOutCMS(), ql=pair.QLongCMS();
      Double_t qspf=pair.QSidePf(),qopf=pair.QOutPf(),qlpf=pair.QLongPf() ;
      
      Double_t pairPhi=TMath::ATan2(ph1->Py()+ph2->Py(),ph1->Px()+ph2->Px()) ;
      Double_t dPsi = rpFull-pairPhi ;
      while(dPsi<0)dPsi+=TMath::Pi() ;
      while(dPsi>TMath::Pi())dPsi-=TMath::Pi() ;
      
      // Yano-Koonin-Podgoretskii Parametrisation 
      Double_t qP=0., qT=0., q0=0. ;
      // source rest frame (usually lab frame)
      pair.QYKPCMS(qP, qT, q0);

      Double_t qPpf=0., qTpf=0., q0pf=0. ;
      // longitudinal comoving frame
        pair.QYKPPF(qPpf,qTpf,q0pf) ;

	
      for(Int_t iCut=0; iCut<4; iCut++){
	if(!PairCut(ph1,ph2,iCut))
	    continue ;
	
        FillHistogram(Form("hetaphi2D_%s_%s",cut[iCut],kTbin.Data()),dEta,dPhi) ;
 
	if(ph1->Module()!=ph2->Module())
          continue ;
	
	
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
	  
	FillHistogram(Form("hQinv_%s",cut[iCut]),qinv,kT) ;
	if(TMath::Abs(qo) < 0.05)
	  FillHistogram(Form("hQinvCut_%s",cut[iCut]),qinv,kT) ;

        // Bertsch-Pratt momentum components in Pair Frame - written by Bekele/Humanic
//        FillHistogram(Form("hOSLPF_%s_%s",cut[iCut],kTbin.Data()),qspf,qopf,qlpf) ;
   
        // Bertsch-Pratt momentum components in Local CMS (longitudinally comoving) frame
        FillHistogram(Form("hOSLCMS_%s_%s",cut[iCut],kTbin.Data()),qs,qo,ql) ;
        FillHistogram(Form("hetaphi_%s_%s",cut[iCut],kTbin.Data()),dEta,dPhi,dE) ;
	if(TMath::Abs(dEta)>0.02)
          FillHistogram(Form("hetaphiRP_%s_%s",cut[iCut],kTbin.Data()),dPsi,dPhi) ;
	
//        FillHistogram(Form("hdXdZ_%s_%s",cut[iCut],kTbin.Data()),dX,dZ) ;


//        FillHistogram(Form("hYKPCMS_%s_%s",cut[iCut],kTbin.Data()),qP, qT, q0);       
      
//        FillHistogram(Form("hYKPPF_%s_%s",cut[iCut],kTbin.Data()),qPpf, qTpf, q0pf);       
        
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
    
    for(Int_t ev=0; ev<prevPHOS->GetSize();ev++){
      TClonesArray * mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev)) ;
      for(Int_t i2=0; i2<mixPHOS->GetEntriesFast();i2++){
	AliCaloPhoton * ph2=(AliCaloPhoton*)mixPHOS->At(i2) ;
	
//	if(!PairCut(ph1,ph2,kDefault))
//	  continue;
	
        AliFemtoTrack track2;
        AliFemtoThreeVector mom2;
        mom2.SetX(ph2->Px()) ;
        mom2.SetY(ph2->Py()) ;
        mom2.SetZ(ph2->Pz()) ;
        track2.SetP(mom2) ;
        AliFemtoParticle part2(&track2,kgMass) ;
       
	AliFemtoParticle *a = &part1 ;
        AliFemtoParticle *b = &part2 ;
        Double_t dEta = ph1->Eta()-ph2->Eta() ; 
        Double_t dPhi = ph1->Phi()-ph2->Phi() ; 
        Double_t dE   = ph1->E()  - ph2->E() ;
        Double_t dX = TMath::Power(ph1->EMCx() - ph2->EMCx(),2) + TMath::Power(ph1->EMCy() - ph2->EMCy(),2)  ;
        dX=TMath::Sign(TMath::Sqrt(dX),ph1->EMCx() - ph2->EMCx()) ;
        Double_t dZ = ph1->EMCz() - ph2->EMCz() ;
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
        TString kTbin ;
	Int_t ikTbin=0 ;
        if(kT<0.2){ kTbin="Kt00-02"; ikTbin=0; }
        else if(kT<0.4){  kTbin="Kt02-04"; ikTbin=1; }
        else if(kT<0.7){  kTbin="Kt04-07"; ikTbin=2; }
        else if(kT<1.0){  kTbin="Kt07-10"; ikTbin=3; }
        else if(kT<1.3){  kTbin="Kt10-13"; ikTbin=4; }
        else if(kT<2.0){  kTbin="Kt13-20"; ikTbin=5; }
        else  continue;
      
      Double_t qs=pair.QSideCMS(), qo=pair.QOutCMS(), ql=pair.QLongCMS();
      Double_t qspf=pair.QSidePf(),qopf=pair.QOutPf(),qlpf=pair.QLongPf() ;
      
      Double_t wMix = EtaPhiWeight(ikTbin,dPhi );

      Double_t pairPhi=TMath::ATan2(ph1->Py()+ph2->Py(),ph1->Px()+ph2->Px()) ;
      Double_t dPsi = rpFull-pairPhi ;
      while(dPsi<0)dPsi+=TMath::Pi() ;
      while(dPsi>TMath::Pi())dPsi-=TMath::Pi() ;
      
      // Yano-Koonin-Podgoretskii Parametrisation 
      Double_t qP=0., qT=0., q0=0. ;
      // source rest frame (usually lab frame)
      pair.QYKPCMS(qP, qT, q0);

      Double_t qPpf=0., qTpf=0., q0pf=0. ;
      // longitudinal comoving frame
        pair.QYKPPF(qPpf,qTpf,q0pf) ;
	
	for(Int_t iCut=0; iCut<4; iCut++){
   	  if(!PairCut(ph1,ph2,iCut))
	    continue ;
	  
          FillHistogram(Form("hMietaphi2D_%s_%s",cut[iCut],kTbin.Data()),dEta,dPhi) ;
 
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
	  
	  FillHistogram(Form("hMiQinv_%s",cut[iCut]),qinv,kT) ;
	  FillHistogram(Form("hMi2Qinv_%s",cut[iCut]),qinv,kT,wMix) ;
	   if(TMath::Abs(qo) < 0.05){
	     FillHistogram(Form("hMiQinvCut_%s",cut[iCut]),qinv,kT) ;
	     FillHistogram(Form("hMi2QinvCut_%s",cut[iCut]),qinv,kT,wMix) ;
	   }
          // Bertsch-Pratt momentum components in Pair Frame - written by Bekele/Humanic
//          FillHistogram(Form("hMiOSLPF_%s_%s",cut[iCut],kTbin.Data()),qspf,qopf,qlpf) ;
   
          // Bertsch-Pratt momentum components in Local CMS (longitudinally comoving) frame
          FillHistogram(Form("hMiOSLCMS_%s_%s",cut[iCut],kTbin.Data()),qs,qo,ql) ;
          FillHistogram(Form("hMi2OSLCMS_%s_%s",cut[iCut],kTbin.Data()),qs,qo,ql,wMix) ;

          FillHistogram(Form("hMietaphi_%s_%s",cut[iCut],kTbin.Data()),dEta,dPhi,dE) ;
	  if(TMath::Abs(dEta)>0.02)
            FillHistogram(Form("hMietaphiRP_%s_%s",cut[iCut],kTbin.Data()),dPsi,dPhi) ;
          FillHistogram(Form("hMi2etaphi_%s_%s",cut[iCut],kTbin.Data()),dEta,dPhi,dE,wMix) ;
//          FillHistogram(Form("hMidXdZ_%s_%s",cut[iCut],kTbin.Data()),dX,dZ) ;
	  
//          FillHistogram(Form("hMiYKPCMS_%s_%s",cut[iCut],kTbin.Data()),qP, qT, q0);       
      
//          FillHistogram(Form("hMiYKPPF_%s_%s",cut[iCut],kTbin.Data()),qPpf, qTpf, q0pf);       
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
  // Draw result to the screen
  // Called once at the end of the query
  
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
  
 // if(cut==kDefault){
  //Consider only pairs from same mudule
//  if(ph1->Module()!=ph2->Module())
//    return kFALSE ;
  
  if(cut==0){
    return kTRUE ;
  }
  if(cut==1){
    return ph1->IsDispOK() && ph2->IsDispOK() ;  
  }
  if(cut==2){
    return ph1->IsCPVOK() && ph2->IsCPVOK() ;  
  }
  if(cut==3){
    return ph1->IsDispOK() && ph2->IsDispOK() && ph1->IsCPVOK() && ph2->IsCPVOK() ;  
  }
  
  
  //Distance between clusters in PHOS plane
  if(cut==4 ||cut==5 || cut==6 ){
    if(!(ph1->IsDispOK() && ph2->IsDispOK() && ph1->IsCPVOK() && ph2->IsCPVOK()))
      return kFALSE ;
/*
    Double_t dx = ph1->EMCx()-ph2->EMCx() ;
    Double_t dz = ph1->EMCz()-ph2->EMCz() ;
    if(cut==4)
      return (dx*dx+dz*dz > 5.*5.) ;
    if(cut==5)
      return (dx*dx+dz*dz > 10.*10.) ;
    if(cut==6)
      return (dx*dx+dz*dz > 15.*15.) ;
    */
  }
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
      result=result|1<<0 ;
      if(aodTrack->Pt()>cutPt2)
	result=result|1<<1 ;
        if(aodTrack->Pt()>cutPt1){
	  result=result|1<<2 ;
	  return result ;
	}
      
    }

  }
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

