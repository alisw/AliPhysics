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
#include "AliAnalysisTaskggMC.h"
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

ClassImp(AliAnalysisTaskggMC)

//________________________________________________________________________
AliAnalysisTaskggMC::AliAnalysisTaskggMC(const char *name) 
: AliAnalysisTaskgg(name),
  fStack(0x0),fMCEvent(0x0),fMCEvents(0x0),fMCEventH(0x0),fMCEventsH(0x0)
{

}

//________________________________________________________________________
void AliAnalysisTaskggMC::UserCreateOutputObjects()
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
  fOutputContainer->Add(new TH2F("hCenPHOS","Centrality vs PHOSclusters", 100,0.,100.,200,0.,200.)) ;
  fOutputContainer->Add(new TH2F("hCenPHOSCells","Centrality vs PHOS cells", 100,0.,100.,100,0.,1000.)) ;
  fOutputContainer->Add(new TH2F("hCenTrack","Centrality vs tracks", 100,0.,100.,100,0.,15000.)) ;  
  fOutputContainer->Add(new TH2F("hCluEvsClu_All","ClusterMult vs E",20,0.,2.,50,0.,50.)) ;
  fOutputContainer->Add(new TH2F("hCluEvsClu_CPV","ClusterMult vs E",20,0.,2.,50,0.,50.)) ;
  fOutputContainer->Add(new TH2F("hCluEvsClu_Disp","ClusterMult vs E",20,0.,2.,50,0.,50.)) ;
  fOutputContainer->Add(new TH2F("hCluEvsClu_Both","ClusterMult vs E",20,0.,2.,50,0.,50.)) ;
  fOutputContainer->Add(new TH2F("hCluEvsCluM","ClusterMult vs E",200,0.,20.,100,0.,20.)) ;
    
  //PHOS QA
			
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
  
  const Int_t nParents=12 ;
  char parent[nParents][20] ;
  sprintf(parent[0],"Conversion") ;
  sprintf(parent[1],"pi0") ;
  sprintf(parent[2],"pipm") ;
  sprintf(parent[3],"p") ;
  sprintf(parent[4],"n") ;
  sprintf(parent[5],"eta") ;
  sprintf(parent[6],"K0s") ;
  sprintf(parent[7],"K0l") ;
  sprintf(parent[8],"omega") ;
  sprintf(parent[9],"Sigma") ;
  sprintf(parent[10],"Charged") ;
  sprintf(parent[11],"Neutral") ;
  
  
  fOutputContainer->Add(new TH2F("hQinv_MC","Qinv distribution",200,0.,0.5,20,0.,2.));
  fOutputContainer->Add(new TH2F("hMiQinv_MC","Qinv distribution",200,0.,0.5,20,0.,2.));
  fOutputContainer->Add(new TH2F("hQinv_MChh","Qinv distribution",200,0.,0.5,20,0.,2.));
  fOutputContainer->Add(new TH2F("hMiQinv_MChh","Qinv distribution",200,0.,0.5,20,0.,2.));
  fOutputContainer->Add(new TH2F("hQinv_MCgh","Qinv distribution",200,0.,0.5,20,0.,2.));
  fOutputContainer->Add(new TH2F("hMiQinv_MCgh","Qinv distribution",200,0.,0.5,20,0.,2.));
  
  for(Int_t iCut=0; iCut<nCuts; iCut++){
    for(Int_t ikT=0; ikT<6; ikT++){ 
//      fOutputContainer->Add(new TH3F(Form("hOSLPF_%s_%s",cut[iCut],kTbins[ikT]),"Out-Side-Long, Pair Frame",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));
 
//      fOutputContainer->Add(new TH3F(Form("hOSLCMS_%s_%s",cut[iCut],kTbins[ikT]),"Out-Side-Long, CMS",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));
//      fOutputContainer->Add(new TH3F(Form("hYKPPF_%s_%s",cut[iCut],kTbins[ikT]),"YKP, Pair Frame",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));
//      fOutputContainer->Add(new TH3F(Form("hYKPCMS_%s_%s",cut[iCut],kTbins[ikT]),"YKP, CMS",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));

//      fOutputContainer->Add(new TH3F(Form("hetaphi_%s_%s",cut[iCut],kTbins[ikT]),"Eta-phi-E correlations",100,-0.25,0.25,100,-TMath::Pi()/6.,TMath::Pi()/6.,20,-0.2,0.2));
//      fOutputContainer->Add(new TH2F(Form("hdXdZ_%s_%s",cut[iCut],kTbins[ikT]),"dXdZ",200,-200,200,200,-200.,200.));
      
      
//      fOutputContainer->Add(new TH3F(Form("hMiOSLPF_%s_%s",cut[iCut],kTbins[ikT]),"Out-Side-Long, Pair Frame",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));

//      fOutputContainer->Add(new TH3F(Form("hMiOSLCMS_%s_%s",cut[iCut],kTbins[ikT]),"Out-Side-Long, CMS",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));
//      fOutputContainer->Add(new TH3F(Form("hMiYKPPF_%s_%s",cut[iCut],kTbins[ikT]),"YKP, Pair Frame",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));
//      fOutputContainer->Add(new TH3F(Form("hMiYKPCMS_%s_%s",cut[iCut],kTbins[ikT]),"YKP, CMS",nQ,-qMax,qMax,nQ,-qMax,qMax,nQ,-qMax,qMax));

//      fOutputContainer->Add(new TH3F(Form("hMietaphi_%s_%s",cut[iCut],kTbins[ikT]),"Eta-phi-E correlations",100,-0.25,0.25,100,-TMath::Pi()/6.,TMath::Pi()/6.,20,-0.2,0.2));
//      fOutputContainer->Add(new TH2F(Form("hMidXdZ_%s_%s",cut[iCut],kTbins[ikT]),"dXdZ",200,-200,200,200,-200.,200.));
    
    }        

    //General
    fOutputContainer->Add(new TH2F(Form("hQinv_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,20,0.,2.));
    fOutputContainer->Add(new TH2F(Form("hMiQinv_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,20,0.,2.));
    fOutputContainer->Add(new TH2F(Form("hQinv_PhotPHOS_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,20,0.,2.));
    fOutputContainer->Add(new TH2F(Form("hMiQinv_PhotPHOS_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,20,0.,2.));
    fOutputContainer->Add(new TH2F(Form("hQinv_PhotVert_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,20,0.,2.));
    fOutputContainer->Add(new TH2F(Form("hQinv_PiVert_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,20,0.,2.));
    fOutputContainer->Add(new TH2F(Form("hQinv_ChargeVert_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,20,0.,2.));
    fOutputContainer->Add(new TH2F(Form("hQinv_NeutralVert_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,20,0.,2.));
    fOutputContainer->Add(new TH2F(Form("hQinv_PhotPiVert_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,20,0.,2.));
    fOutputContainer->Add(new TH2F(Form("hQinv_PhotChargeVert_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,20,0.,2.));
    fOutputContainer->Add(new TH2F(Form("hQinv_PhotNeutralVert_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,20,0.,2.));
    fOutputContainer->Add(new TH2F(Form("hQinv_PhotAnyVert_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,20,0.,2.));
    fOutputContainer->Add(new TH2F(Form("hQinv_PiChargeVert_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,20,0.,2.));
    fOutputContainer->Add(new TH2F(Form("hQinv_PiNeutralVert_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,20,0.,2.));
    fOutputContainer->Add(new TH2F(Form("hQinv_PiAnyVert_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,20,0.,2.));
    fOutputContainer->Add(new TH2F(Form("hQinv_ChargeAnyVert_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,20,0.,2.));
    fOutputContainer->Add(new TH2F(Form("hQinv_NeutralAnyVert_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,20,0.,2.));
    fOutputContainer->Add(new TH2F(Form("hQinv_Other_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,20,0.,2.));
            
    
    //Common parent
    for(Int_t ip=0; ip<nParents; ip++){
      fOutputContainer->Add(new TH2F(Form("hQinv_%s_%s",parent[ip],cut[iCut]),"Qinv distribution",200,0.,0.5,20,0.,2.));
    }
    
    //true photons wo common parents
   fOutputContainer->Add(new TH2F(Form("hQinv_gamma_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,20,0.,2.));  
    
    //pi0 HBT
    fOutputContainer->Add(new TH2F(Form("hQinv_pi0HBT_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,20,0.,2.));

    //Flow
    fOutputContainer->Add(new TH2F(Form("hQinv_flow_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,20,0.,2.));

    //eta-phi
    fOutputContainer->Add(new TH2F(Form("hQinv_etaphi_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,20,0.,2.));

    //jet
    fOutputContainer->Add(new TH2F(Form("hQinv_jet_%s",cut[iCut]),"Qinv distribution",200,0.,0.5,20,0.,2.));
    
    
  }

//   for(Int_t ikT=0; ikT<6; ikT++){ 
//      fOutputContainer->Add(new TH2F(Form("hSLfine_%s",kTbins[ikT]),"Out-Side",1000,-0.5,0.5,1000,-0.5,0.5));
//      fOutputContainer->Add(new TH2F(Form("hMiSLfine_%s",kTbins[ikT]),"Out-Side",1000,-0.5,0.5,1000,-0.5,0.5));
//      fOutputContainer->Add(new TH3F(Form("hSLr_%s",kTbins[ikT]),"Side-Long-r",nQ,-qMax,qMax,nQ,-qMax,qMax,30,0.,30.));
//      fOutputContainer->Add(new TH3F(Form("hMiSLr_%s",kTbins[ikT]),"Side-Long-r",nQ,-qMax,qMax,nQ,-qMax,qMax,30,0.,30.));
//   }
  
    
  if(!fMCEvents)fMCEvents=new TList() ;
  if(!fMCEventsH)fMCEventsH=new TList() ;
 
  PostData(1, fOutputContainer);

}

//________________________________________________________________________
void AliAnalysisTaskggMC::UserExec(Option_t *) 
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

  fStack = (TClonesArray*)fEvent->FindListObject(AliAODMCParticle::StdBranchName());
  if (!fStack) {
    Printf("ERROR: Could not retrieve stack");
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
/*  
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
*/    
  
 
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
  
   ProcessMC() ;
  
  
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
    
    ph->SetPrimary(FindAODLabel(clu->GetLabelAt(0))) ;
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
    
    AliAODMCParticle * hadronParent1 = 0x0 ;
    Bool_t isPhotonAtPHOS1 = kFALSE; 
    Bool_t isPhotonAtVertex1 = kFALSE; 
    Bool_t isPionAtVertex1 = kFALSE; 
    Bool_t isChHadronAtVertex1 = kFALSE; 
    Bool_t isNeuHadronAtVertex1 = kFALSE; 

    if(ph1->GetPrimary()>=0){
       AliAODMCParticle *tmp = (AliAODMCParticle*)fStack->At(ph1->GetPrimary()) ;
       if(tmp->GetPdgCode()==22)
	 isPhotonAtPHOS1 = kTRUE ;
       //Look at vertex    
       Double_t r=TMath::Sqrt(tmp->Xv()*tmp->Xv()+tmp->Yv()*tmp->Yv()) ;
       while(r>1. && tmp->GetMother()>-1){
         tmp = (AliAODMCParticle*)fStack->At(tmp->GetMother()) ;
	 r=TMath::Sqrt(tmp->Xv()*tmp->Xv()+tmp->Yv()*tmp->Yv()) ;
       }
       if(tmp->GetPdgCode()==22)
	 isPhotonAtVertex1 = kTRUE ;
       else
         if(TMath::Abs(tmp->GetPdgCode())==211)
	   isPionAtVertex1 = kTRUE ;
	 else
	   if(tmp->Charge()==0)
	     isNeuHadronAtVertex1 =kTRUE ;
	   else
	     isChHadronAtVertex1 =kTRUE ;
       //Finds hadron parent
       tmp = (AliAODMCParticle*)fStack->At(ph1->GetPrimary()) ;
       while(TMath::Abs(tmp->GetPdgCode())<100){ //gamma,mu,e
         if(tmp->GetMother()>=0){
	    tmp = (AliAODMCParticle*)fStack->At(tmp->GetMother()) ;
	 }
	 else{
	   tmp=0x0 ;
	   break ;
	 }
       }
       hadronParent1 = tmp ;
    }
    
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
	if(!PairCut(ph1,ph2,kBoth))
	  continue;
     
//      if(!SecondaryPi0Cut(ph1,ph2))
//        continue ;

      //hadron, whose daughter created this cluster
      AliAODMCParticle * hadronParent2 = 0x0 ;
      Bool_t isPhotonAtPHOS2 = kFALSE; 
      Bool_t isPhotonAtVertex2 = kFALSE; 
      Bool_t isPionAtVertex2 = kFALSE; 
      Bool_t isChHadronAtVertex2 = kFALSE; 
      Bool_t isNeuHadronAtVertex2 = kFALSE; 
      if(ph2->GetPrimary()>=0){
         AliAODMCParticle *tmp = (AliAODMCParticle*)fStack->At(ph2->GetPrimary()) ;
         if(tmp->GetPdgCode()==22)
  	   isPhotonAtPHOS2 = kTRUE ;
         //Look at vertex    
         Double_t r=TMath::Sqrt(tmp->Xv()*tmp->Xv()+tmp->Yv()*tmp->Yv()) ;
         while(r>1. && tmp->GetMother()>-1){
           tmp = (AliAODMCParticle*)fStack->At(tmp->GetMother()) ;
	   r=TMath::Sqrt(tmp->Xv()*tmp->Xv()+tmp->Yv()*tmp->Yv()) ;
         }
         if(tmp->GetPdgCode()==22)
	   isPhotonAtVertex2 = kTRUE ;
         else
           if(TMath::Abs(tmp->GetPdgCode())==211)
	     isPionAtVertex2 = kTRUE ;
	   else
	     if(tmp->Charge()==0)
	       isNeuHadronAtVertex2 =kTRUE ;
	     else
	       isChHadronAtVertex2 =kTRUE ;
	 
         //Finds hadron parent
         tmp = (AliAODMCParticle*)fStack->At(ph2->GetPrimary()) ;
         while(TMath::Abs(tmp->GetPdgCode())<100){ //gamma,mu,e
           if(tmp->GetMother()>=0){
	      tmp = (AliAODMCParticle*)fStack->At(tmp->GetMother()) ;
	   }
	   else{
	     tmp=0x0 ;
	     break ;
	   }
         }
         hadronParent2 = tmp ;
      }
	
	
	
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
      
      // Yano-Koonin-Podgoretskii Parametrisation 
      Double_t qP=0., qT=0., q0=0. ;
      // source rest frame (usually lab frame)
      pair.QYKPCMS(qP, qT, q0);

      Double_t qPpf=0., qTpf=0., q0pf=0. ;
      // longitudinal comoving frame
        pair.QYKPPF(qPpf,qTpf,q0pf) ;

	
      //Test common parent
     AliAODMCParticle * commonPrim = TestCommonParent(ph1, ph2); 
     if(commonPrim){
       Int_t parentPDG = commonPrim->GetPdgCode() ;
       char key[55];
       switch(parentPDG){
	 case 22:
	 case 11:
	 case -11: //consversion 
	     snprintf(key,55,"Conversion") ;
	     break ;
	 case 111:
	     snprintf(key,55,"pi0") ;
	     break ;
	 case  211:
	 case -211:
	     snprintf(key,55,"pipm") ;
	     break ;
	 case 221:
	     snprintf(key,55,"eta") ;
	     break ;
	 case  2212:
	 case -2212:
	     snprintf(key,55,"p") ;
	     break ;
	 case  2112:
	 case -2112:
	     snprintf(key,55,"n") ;
	     break ;
	 case 310:
	     snprintf(key,55,"K0s") ;
	     break ;
	 case 130:
	     snprintf(key,55,"K0l") ;
	     break ;
	 case 223:
	     snprintf(key,55,"omega") ;
	     break ;
	 case  3222:
	 case -3222:
	 case  3212:
	 case -3212:
	 case  3112:
	 case -3112:
	     snprintf(key,55,"Sigma") ;
	     break ;
	 default:
	   if(TMath::Abs(parentPDG)<100){ //parton
 	     snprintf(key,55,"jet") ;	     
	   }
	   else{
	     if(commonPrim->Charge()!=0)
 	       snprintf(key,55,"Charged") ;
	     else
 	       snprintf(key,55,"Neutral") ;
	   }
         }
         for(Int_t iCut=0; iCut<4; iCut++){
   	   if(!PairCut(ph1,ph2,iCut))
	    continue ;	  
	   FillHistogram(Form("hQinv_%s_%s",key,cut[iCut]),qinv,kT) ;
	 }
       }



       //Flow
       if(hadronParent1 && hadronParent2){
         //pi0 HBT
         if((hadronParent1->GetPdgCode()==111) && (hadronParent2->GetPdgCode()==111)){
           Double_t hbtWeight = PionHBTWeight(hadronParent1,hadronParent2); 
           for(Int_t iCut=0; iCut<4; iCut++){
   	     if(!PairCut(ph1,ph2,iCut))
	       continue ;	  
	     FillHistogram(Form("hQinv_pi0HBT_%s",cut[iCut]),qinv,kT,hbtWeight) ;
	   }	 
         }
	 
	 
	 Double_t flowWeight =FlowWeight(hadronParent1, hadronParent2) ;
	 Double_t etaWeight = EtaPhiWeight(hadronParent1, hadronParent2) ;
         for(Int_t iCut=0; iCut<4; iCut++){
   	   if(!PairCut(ph1,ph2,iCut))
	     continue ;	  
	   FillHistogram(Form("hQinv_flow_%s",cut[iCut]),qinv,kT,flowWeight) ;
	   FillHistogram(Form("hQinv_etaphi_%s",cut[iCut]),qinv,kT,etaWeight) ;	   
	 }
       }
	        
       	
      for(Int_t iCut=0; iCut<4; iCut++){
   	if(!PairCut(ph1,ph2,iCut))
	    continue ;
		  
	FillHistogram(Form("hQinv_%s",cut[iCut]),qinv,kT) ;
	if(isPhotonAtPHOS1 && isPhotonAtPHOS2) 
	   FillHistogram(Form("hQinv_PhotPHOS_%s",cut[iCut]),qinv,kT) ;
	if(isPhotonAtVertex1 && isPhotonAtVertex2){ 
	   FillHistogram(Form("hQinv_PhotVert_%s",cut[iCut]),qinv,kT) ;
	}
	else{
	  if((isPionAtVertex1 && isPhotonAtVertex2) || (isPhotonAtVertex1 && isPionAtVertex2)){
	    FillHistogram(Form("hQinv_PhotPiVert_%s",cut[iCut]),qinv,kT) ;
	  }
	  else{
	    if((isChHadronAtVertex1&& isPhotonAtVertex2) || (isPhotonAtVertex1 && isChHadronAtVertex2)){
	      FillHistogram(Form("hQinv_PhotChargeVert_%s",cut[iCut]),qinv,kT) ;
	    }
	    else{
	      if((isNeuHadronAtVertex1 && isPhotonAtVertex2) || (isPhotonAtVertex1 && isNeuHadronAtVertex2)){
	       FillHistogram(Form("hQinv_PhotNeutralVert_%s",cut[iCut]),qinv,kT) ;
	      }
	      else{
	        if(isPhotonAtVertex1 || isPhotonAtVertex2){
	          FillHistogram(Form("hQinv_PhotAnyVert_%s",cut[iCut]),qinv,kT) ;
	        }
	        else{
	          if(isPionAtVertex1 && isPionAtVertex2){
	            FillHistogram(Form("hQinv_PiVert_%s",cut[iCut]),qinv,kT) ;
		  }
		  else{
	            if((isPionAtVertex1 && isChHadronAtVertex2) || (isChHadronAtVertex1 && isPionAtVertex2)){
	              FillHistogram(Form("hQinv_PiChargeVert_%s",cut[iCut]),qinv,kT) ;
		    }
		    else{
	              if((isPionAtVertex1 && isNeuHadronAtVertex2) || (isNeuHadronAtVertex1 && isPionAtVertex2)){
	                FillHistogram(Form("hQinv_PiNeutralVert_%s",cut[iCut]),qinv,kT) ;
		      }
		      else{
	                if(isPionAtVertex1  || isPionAtVertex2){
	                  FillHistogram(Form("hQinv_PiAnyVert_%s",cut[iCut]),qinv,kT) ;
		        }
		        else{		    
	                  if(isChHadronAtVertex1 && isChHadronAtVertex2){
	                    FillHistogram(Form("hQinv_ChargeVert_%s",cut[iCut]),qinv,kT) ;
		          }
		          else{
	                    if(isChHadronAtVertex1 || isChHadronAtVertex2){
	                      FillHistogram(Form("hQinv_ChargeAny_%s",cut[iCut]),qinv,kT) ;
		            }
		            else{
	                      if(isNeuHadronAtVertex1 && isNeuHadronAtVertex2){
	                         FillHistogram(Form("hQinv_NeutralVert_%s",cut[iCut]),qinv,kT) ;
		              }
		              else{
	                        if(isNeuHadronAtVertex1 || isNeuHadronAtVertex2){
	                          FillHistogram(Form("hQinv_NeutralAnyVert_%s",cut[iCut]),qinv,kT) ;
		                }
		                else{
	                          FillHistogram(Form("hQinv_Other_%s",cut[iCut]),qinv,kT) ;	
			        }
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
	
	
//        FillHistogram(Form("hOSLCMS_%s_%s",cut[iCut],kTbin.Data()),qs,qo,ql) ;
//        FillHistogram(Form("hetaphi_%s_%s",cut[iCut],kTbin.Data()),dEta,dPhi,dE) ;
//        FillHistogram(Form("hdXdZ_%s_%s",cut[iCut],kTbin.Data()),dX,dZ) ;
        
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
	
	if(!PairCut(ph1,ph2,kDefault))
	  continue;

	
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
        TString kTbin="15" ;
        if(kT<0.2) kTbin="Kt00-02";
        else if(kT<0.4) kTbin="Kt02-04";
        else if(kT<0.7) kTbin="Kt04-07";
        else if(kT<1.) kTbin="Kt07-10";
        else if(kT<1.3) kTbin="Kt10-13";
        else if(kT<2.0) kTbin="Kt13-20";
        else  continue;
      
      Double_t qs=pair.QSideCMS(), qo=pair.QOutCMS(), ql=pair.QLongCMS();
      Double_t qspf=pair.QSidePf(),qopf=pair.QOutPf(),qlpf=pair.QLongPf() ;
      
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
	
	  
	  FillHistogram(Form("hMiQinv_%s",cut[iCut]),qinv,kT) ;

//           FillHistogram(Form("hMiOSLCMS_%s_%s",cut[iCut],kTbin.Data()),qs,qo,ql) ;
//           FillHistogram(Form("hMietaphi_%s_%s",cut[iCut],kTbin.Data()),dEta,dPhi,dE) ;
//           FillHistogram(Form("hMidXdZ_%s_%s",cut[iCut],kTbin.Data()),dX,dZ) ;
	  
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
void AliAnalysisTaskggMC::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
  
}


//___________________________________________________________________________
Double_t AliAnalysisTaskggMC::PionHBTWeight(const AliAODMCParticle * ph1, const AliAODMCParticle * ph2) const{
  //Calculates weight of the HBT correlations
  //assuming that particles are pions
  //Parameterization from 
// article{Aamodt:2011mr,
//       author         = "Aamodt, K. and others",
//       title          = "{Two-pion Bose-Einstein correlations in central Pb-Pb
//                         collisions at $\sqrt{{s}_{NN}} =$ 2.76 TeV}",
//       journal        = "Phys. Lett.",
//       volume         = "B696",
//       year           = "2011",
//       pages          = "328-337",
//       doi            = "10.1016/j.physletb.2010.12.053",
//       eprint         = "1012.4035",
// }
  
    const Double_t kpiMass=0.135 ;

    AliFemtoTrack track1;
    AliFemtoThreeVector mom1;
    mom1.SetX(ph1->Px()) ;
    mom1.SetY(ph1->Py()) ;
    mom1.SetZ(ph1->Pz()) ;
    track1.SetP(mom1) ;
    AliFemtoParticle part1(&track1,kpiMass) ;
    
    AliFemtoTrack track2;
    AliFemtoThreeVector mom2;
    mom2.SetX(ph2->Px()) ;
    mom2.SetY(ph2->Py()) ;
    mom2.SetZ(ph2->Pz()) ;
    track2.SetP(mom2) ;
    AliFemtoParticle part2(&track2,kpiMass) ;
       
    AliFemtoParticle *a = &part1 ;
    AliFemtoParticle *b = &part2 ;
    if(gRandom->Uniform()>0.5){
      a = &part2 ;
      b = &part1 ;
    }
    AliFemtoPair pair(a,b);

    Double_t kT = pair.KT() ;  
    Double_t qs=pair.QSideCMS(), qo=pair.QOutCMS(), ql=pair.QLongCMS();
  
    Double_t Ro = TMath::Sqrt(8.583283 + 5.453621/kT + 1.269961/kT/kT)/0.197 ;
    Double_t Rs = TMath::Sqrt(TMath::Max(0.953404 + 16.412315/kT + -1.440726/kT/kT,0.))/0.197 ; 
    Double_t Rl = TMath::Sqrt(TMath::Max(-3.843901 + 21.680994/kT + -0.942183/kT/kT,0.))/0.197 ; 
  
    Double_t r=  qs*qs*Rs*Rs + qo*qo*Ro*Ro + ql*ql*Rl*Rl ; 
    if(r>50.) return 1.;
    else  
       return 1. + TMath::Exp(-r) ;


}
//___________________________________________________________________________
Double_t AliAnalysisTaskggMC::FlowWeight(const AliAODMCParticle * ph1, const AliAODMCParticle * ph2) const{
  //Calculates weight of the HBT correlations
  //assuming that particles are pions
//   TF1 * fv2 = new TF1("fv2","[0]*(1.+[1]*x+[2]*x*x+[3]*x*x*x)/(1+[4]*x+[5]*x*x+[6]*x*x*x)",0.,30);  
//   // 0-5%:
//   if(cen==0) //Fit DP .
//     fv2->SetParameters(-5.732671e-03, -1.132112e+01, 2.977289e+00, -2.351123e-01, 2.422317e-01, -1.353961e-01, 2.552272e-02); 
//   //5-10%:
//   if(cen==1)
//     fv2->SetParameters(-7.227483e-03, -1.327021e+01, 1.192188e+00, -1.424606e-01, 2.690662e-01, -6.946581e-02, 3.566327e-02); 
//   //10-20%:
//   if(cen==2)
//     fv2->SetParameters(-1.062719e-02, -1.339370e+01, 2.393260e+00, -2.391584e-01, 1.998852e-01, -8.818621e-02, 3.052334e-02); 
//   //20-40%:
//   if(cen==3)
//     fv2->SetParameters(-1.605216e-02, -1.238932e+01, -2.932423e-01, -8.134508e-02, 4.570368e-01, -9.297761e-02, 6.322929e-02); 
//  //30-40%:
//   if(cen==4)
//     fv2->SetParameters(-2.284108e-02, -1.073473e+01, 1.355141e+00, -2.276158e-01, 2.680910e-01, -3.884307e-02, 4.643840e-02); 
//   //40-50%:
//   if(cen==5)
//     fv2->SetParameters(-1.958364e-02, -9.395222e+00, -1.042019e+01, -3.846221e-03, 1.043175e+00, 2.247645e-01, 2.563315e-01); 
//   //0-20%:
//   if(cen==6)
//     fv2->SetParameters(-8.399737e-03, -1.274228e+01, 2.237830e+00, -2.186378e-01, 2.452230e-01, -1.077267e-01, 3.325352e-02); 
//   //0-10%:
//   if(cen==7)
//     fv2->SetParameters(-6.620259e-03, -1.220510e+01, 2.233523e+00, -2.018025e-01, 2.742919e-01, -1.208619e-01, 3.477000e-02); 
//   if(cen==8)
//     fv2->SetParameters(-2.040839e-02, -1.120741e+01, 1.140065e+00, -1.780665e-01, 1.564463e-01, 4.186171e-02, 3.064791e-02); 
//   if(cen==9)
//     fv2->SetParameters(-1.110217e-02, -1.246538e+01, 1.230970e+00, -1.609303e-01, 3.598792e-01, -1.146195e-01, 4.429495e-02); 
  
  Double_t x = ph1->Pt() ;
  Double_t v2a = 6.620259e-03*(-1.+1.220510e+01*x-2.233523e+00*x*x+2.018025e-01*x*x*x)/(1+2.742919e-01*x-1.208619e-01*x*x+3.477000e-02*x*x*x) ;
 
  x = ph2->Pt() ;
  Double_t v2b = 6.620259e-03*(-1.+1.220510e+01*x-2.233523e+00*x*x+2.018025e-01*x*x*x)/(1+2.742919e-01*x-1.208619e-01*x*x+3.477000e-02*x*x*x) ;

  Double_t dphi = ph1->Phi() - ph2->Phi() ;

  //100 to amplify effect
  v2a=v2a*10. ;
  v2b=v2b*10. ;
  return 1. + 2.*v2a*v2b*TMath::Cos(2.*dphi) ; ;

}
//___________________________________________________________________________
Double_t AliAnalysisTaskggMC::EtaPhiWeight(const AliAODMCParticle * ph1, const AliAODMCParticle * ph2) const{
  //Calculates weight of the HBT correlations
  //assuming that particles are pions
  
  Double_t x = ph1->Pt() ;
  Double_t v2a = 6.620259e-03*(-1.+1.220510e+01*x-2.233523e+00*x*x+2.018025e-01*x*x*x)/(1+2.742919e-01*x-1.208619e-01*x*x+3.477000e-02*x*x*x) ;
  x = ph2->Pt() ;
  Double_t v2b = 6.620259e-03*(-1.+1.220510e+01*x-2.233523e+00*x*x+2.018025e-01*x*x*x)/(1+2.742919e-01*x-1.208619e-01*x*x+3.477000e-02*x*x*x) ;
  Double_t dphi = ph1->Phi() - ph2->Phi() ;
   
  v2a=v2a*10. ;
  v2b=v2b*10. ;
  return 1. + 2.*v2a*v2b*TMath::Cos(2.*dphi) + 0.5*2.*v2a*v2b*TMath::Cos(3.*dphi) + 0.25*2.*v2a*v2b*TMath::Cos(4.*dphi) + + 0.125*2.*v2a*v2b*TMath::Cos(5.*dphi) ;
  
}

//___________________________________________________________________________
AliAODMCParticle* AliAnalysisTaskggMC::TestCommonParent(const AliCaloPhoton * ph1, const AliCaloPhoton * ph2) const{
 //Finds PDG code of common parent
 //return 0 if there is no common
  
  AliAODMCParticle * pdg=0x0 ;
  
  Int_t prim1 = ph1->GetPrimary();
  while(prim1!=-1){ 
    Int_t prim2 = ph2->GetPrimary();  
    while(prim2!=-1){       
      if(prim1==prim2){
	return (AliAODMCParticle*)fStack->At(prim2) ;
      }
      prim2=((AliAODMCParticle*)fStack->At(prim2))->GetMother() ;
    }
    prim1=((AliAODMCParticle*)fStack->At(prim1))->GetMother() ;
  }
 
  
  return pdg ;
}
//___________________________________________________________________________
void AliAnalysisTaskggMC::ProcessMC() {
  
  //Select photons in |eta|<0.5
  if(fMCEvent)
    fMCEvent->Clear() ;
  else
    fMCEvent = new TClonesArray("TLorentzVector",2000) ;
  if(fMCEventH)
    fMCEventH->Clear() ;
  else
    fMCEventH = new TClonesArray("TLorentzVector",2000) ;
  Int_t n=fStack->GetEntriesFast() ;
  Int_t iphot=0, ihadr=0 ;
  for(Int_t i=0 ; i<n; i++){
     AliAODMCParticle *tmp = (AliAODMCParticle*)fStack->At(i) ;
     if(tmp->Pt()<0.1)
       continue ;
     if(TMath::Abs(tmp->Eta())>0.15)
       continue ;
     Double_t r=TMath::Sqrt(tmp->Xv()*tmp->Xv()+tmp->Yv()*tmp->Yv()) ;
     if(r>1.)
       continue ;
     if(tmp->GetPdgCode()==22){
       if(iphot>=fMCEvent->GetSize())
         fMCEvent->Expand(1.5*fMCEvent->GetSize());
       new((*fMCEvent)[iphot++])TLorentzVector(tmp->Px(),tmp->Py(),tmp->Pz(),tmp->E()) ;
     }
     else{
       if(ihadr>=fMCEventH->GetSize())
         fMCEventH->Expand(1.5*fMCEventH->GetSize());
       new((*fMCEventH)[ihadr++])TLorentzVector(tmp->Px(),tmp->Py(),tmp->Pz(),tmp->E()) ;
     }
  }
  
  //Make distribution
  for(Int_t i1=0; i1<iphot-1;i1++){
    TLorentzVector * ph1=(TLorentzVector*)fMCEvent->At(i1) ;
    for(Int_t i2=i1+1; i2<iphot;i2++){
      TLorentzVector * ph2=(TLorentzVector*)fMCEvent->At(i2) ;
	
      FillHistogram("hQinv_MC",(*ph1+ *ph2).M(),0.5*(*ph1+ *ph2).Pt()) ;

    }
  }
  for(Int_t i1=0; i1<ihadr-1;i1++){
    TLorentzVector * ph1=(TLorentzVector*)fMCEventH->At(i1) ;
    for(Int_t i2=i1+1; i2<ihadr;i2++){
      TLorentzVector * ph2=(TLorentzVector*)fMCEventH->At(i2) ;
	
      FillHistogram("hQinv_MChh",(*ph1+ *ph2).M(),0.5*(*ph1+ *ph2).Pt()) ;

    }
  }
  for(Int_t i1=0; i1<iphot;i1++){
    TLorentzVector * ph1=(TLorentzVector*)fMCEvent->At(i1) ;
    for(Int_t i2=i1+1; i2<ihadr;i2++){
      TLorentzVector * ph2=(TLorentzVector*)fMCEventH->At(i2) ;
	
      FillHistogram("hQinv_MCgh",(*ph1+ *ph2).M(),0.5*(*ph1+ *ph2).Pt()) ;

    }
  }
  
  
  for(Int_t i1=0; i1<iphot;i1++){
    TLorentzVector * ph1=(TLorentzVector*)fMCEvent->At(i1) ;
    for(Int_t ev=0; ev<fMCEvents->GetSize();ev++){
      TClonesArray * mixMC = static_cast<TClonesArray*>(fMCEvents->At(ev)) ;
      for(Int_t i2=0; i2<mixMC->GetEntriesFast();i2++){
	TLorentzVector * ph2=(TLorentzVector*)mixMC->At(i2) ;

        FillHistogram("hMiQinv_MC",(*ph1+ *ph2).M(),0.5*(*ph1+ *ph2).Pt()) ;
	
      }
    }
  }
  
  for(Int_t i1=0; i1<iphot;i1++){
    TLorentzVector * ph1=(TLorentzVector*)fMCEvent->At(i1) ;
    for(Int_t ev=0; ev<fMCEventsH->GetSize();ev++){
      TClonesArray * mixMC = static_cast<TClonesArray*>(fMCEventsH->At(ev)) ;
      for(Int_t i2=0; i2<mixMC->GetEntriesFast();i2++){
	TLorentzVector * ph2=(TLorentzVector*)mixMC->At(i2) ;

        FillHistogram("hMiQinv_MCgh",(*ph1+ *ph2).M(),0.5*(*ph1+ *ph2).Pt()) ;
	
      }
    }
  }
  for(Int_t i1=0; i1<ihadr;i1++){
    TLorentzVector * ph1=(TLorentzVector*)fMCEventH->At(i1) ;
    for(Int_t ev=0; ev<fMCEvents->GetSize();ev++){
      TClonesArray * mixMC = static_cast<TClonesArray*>(fMCEvents->At(ev)) ;
      for(Int_t i2=0; i2<mixMC->GetEntriesFast();i2++){
	TLorentzVector * ph2=(TLorentzVector*)mixMC->At(i2) ;

        FillHistogram("hMiQinv_MCgh",(*ph1+ *ph2).M(),0.5*(*ph1+ *ph2).Pt()) ;
	
      }
    }
  }

  for(Int_t i1=0; i1<ihadr;i1++){
    TLorentzVector * ph1=(TLorentzVector*)fMCEventH->At(i1) ;
    for(Int_t ev=0; ev<fMCEventsH->GetSize();ev++){
      TClonesArray * mixMC = static_cast<TClonesArray*>(fMCEvents->At(ev)) ;
      for(Int_t i2=0; i2<mixMC->GetEntriesFast();i2++){
	TLorentzVector * ph2=(TLorentzVector*)mixMC->At(i2) ;

        FillHistogram("hMiQinv_MChh",(*ph1+ *ph2).M(),0.5*(*ph1+ *ph2).Pt()) ;
	
      }
    }
  }
  
  
  
  if(fMCEvent->GetEntriesFast()>0){
    fMCEvents->AddFirst(fMCEvent) ;
    fMCEvent=0;
    if(fMCEvents->GetSize()>1){//Remove redundant events
      TClonesArray * tmp = static_cast<TClonesArray*>(fMCEvents->Last()) ;
      fMCEvents->RemoveLast() ;
      delete tmp ;
    }
  }
  
  if(fMCEventH->GetEntriesFast()>0){
    fMCEventsH->AddFirst(fMCEventH) ;
    fMCEventH=0;
    if(fMCEventsH->GetSize()>1){//Remove redundant events
      TClonesArray * tmp = static_cast<TClonesArray*>(fMCEventsH->Last()) ;
      fMCEventsH->RemoveLast() ;
      delete tmp ;
    }
  }
  
}
//___________________________________________________________________________
Int_t AliAnalysisTaskggMC::FindAODLabel(Int_t esdLabel)const{
   
  if(esdLabel<0)
     return -1 ;
  
  Int_t n = fStack->GetEntriesFast();
  if(esdLabel<n){
      AliAODMCParticle* tmp =  (AliAODMCParticle*)fStack->At(esdLabel) ;
      if(tmp->Label()==esdLabel)
	return esdLabel;
      else{
	Int_t i=esdLabel;
	if(esdLabel>tmp->Label()){
	  i++;
	  while(i<n){
	   tmp =  (AliAODMCParticle*)fStack->At(i) ;
	   if(tmp->Label()==esdLabel)
	     return i;
	   i++;
	  }
	  return -1 ; //failed to find 
	}
	else{
	  i--;
	  while(i>=0){
	   tmp =  (AliAODMCParticle*)fStack->At(i) ;
	   if(tmp->Label()==esdLabel)
	     return i;
	   i--;
	  }
	  return -1 ; //failed to find 
	
	}
      }
      
    }
    else{
       Int_t i=n-1;
       while(i>=0){
	   AliAODMCParticle* tmp =  (AliAODMCParticle*)fStack->At(i) ;
	   if(tmp->Label()==esdLabel)
	     return i;
	   i--;
       }
       return -1 ; //failed to find 
    }
    return -1 ;
  
}