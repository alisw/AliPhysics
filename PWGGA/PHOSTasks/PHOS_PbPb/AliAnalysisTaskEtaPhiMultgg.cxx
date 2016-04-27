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
#include "AliAnalysisTaskEtaPhiMultgg.h"
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

#include "AliConvEventCuts.h"
#include "AliV0ReaderV1.h"
#include "AliConvEventCuts.h"
#include "AliConversionPhotonCuts.h"
#include "AliEMCALGeometry.h"


// Analysis task to fill histograms with PHOS ESD clusters and cells
// Authors: Dmitri Peressounko
// Date   : 28.05.2011

ClassImp(AliAnalysisTaskEtaPhiMultgg)

//________________________________________________________________________
AliAnalysisTaskEtaPhiMultgg::AliAnalysisTaskEtaPhiMultgg(const char *name) 
: AliAnalysisTaskEtaPhigg(name)
{
  // Constructor
  


}

//________________________________________________________________________
void AliAnalysisTaskEtaPhiMultgg::UserCreateOutputObjects()
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
  fOutputContainer->Add(new TH2F("phiRPV0A","Event plane", 100,0.,TMath::Pi(),100,0.,100.)) ;
  fOutputContainer->Add(new TH2F("phiRPV0Aflat","Event plane", 100,0.,TMath::Pi(),100,0.,100.)) ;
  fOutputContainer->Add(new TH2F("phiRPV0C","Event plane", 100,0.,TMath::Pi(),100,0.,100.)) ;
  fOutputContainer->Add(new TH2F("phiRPV0Cflat","Event plane", 100,0.,TMath::Pi(),100,0.,100.)) ;
  fOutputContainer->Add(new TH2F("phiRPV0AC","Event plane", 100,0.,TMath::Pi(),100,0.,TMath::Pi())) ;
  fOutputContainer->Add(new TH2F("phiRPV0ACflat","Event plane", 100,0.,TMath::Pi(),100,0.,TMath::Pi())) ;
  fOutputContainer->Add(new TH2F("phiRPvsV0A","Event plane", 100,0.,TMath::Pi(),100,0.,TMath::Pi())) ;
  fOutputContainer->Add(new TH2F("phiRPvsV0C","Event plane", 100,0.,TMath::Pi(),100,0.,TMath::Pi())) ;

  //vertex distribution
  fOutputContainer->Add(new TH2F("hZvertex","Z vertex position", 50,-25.,25.,nRuns,0.,float(nRuns))) ;
  
  //Centrality
  fOutputContainer->Add(new TH2F("hCentrality","Event centrality", 100,0.,100.,nRuns,0.,float(nRuns))) ;
  fOutputContainer->Add(new TH2F("hCenPHOS","Centrality vs PHOSclusters", 100,0.,100.,200,0.,200.)) ;
  fOutputContainer->Add(new TH2F("hCenPHOSCells","Centrality vs PHOS cells", 100,0.,100.,100,0.,1000.)) ;
  fOutputContainer->Add(new TH2F("hCenTrack","Centrality vs tracks", 100,0.,100.,100,0.,15000.)) ;  
  fOutputContainer->Add(new TH2F("hCenTOF","Centrality vs PHOS TOF", 100,0.,100.,600,-6.e-6,6.e-6)) ;

  fOutputContainer->Add(new TH2F("hCenPCM","Centrality vs PCM photons", 100,0.,100.,200,0.,200.)) ;
  
  fOutputContainer->Add(new TH2F("hCenEMCAL","Centrality vs EMCAL photons", 100,0.,100.,200,0.,200.)) ;

  fOutputContainer->Add(new TH2F("hPHOSvsEMCAL","PHOS vs EMCAL photons", 100,0.,100.,200,0.,400.)) ;

  
  //PHOS QA 			
  
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
  
  //EMCAL QA
  for(Int_t mod=0; mod<12; mod++){
    fOutputContainer->Add(new TH2F(Form("hEMCALmod%d",mod),"EMCAL (x,z)" ,24,0.5,24.5, 48,0.5,58.5));  
  }
  

  const Int_t nCuts=7 ;
  char cut[7][20] ;
  sprintf(cut[0],"All") ;
  sprintf(cut[1],"Disp") ;
  sprintf(cut[2],"CPV") ;
  sprintf(cut[3],"Both") ;
  sprintf(cut[4],"Disp2") ;
  sprintf(cut[5],"CPV2") ;
  sprintf(cut[6],"Both2") ;

  
  for(Int_t iCut=0; iCut<nCuts; iCut++){

    fOutputContainer->Add(new TH3F(Form("hEtaPhiPHOS_%s",cut[iCut]),"Eta-phi-E correlations",100,-1.2,1.2,100,-TMath::Pi()/2.,3.*TMath::Pi()/2.,100,0.,5.));
    fOutputContainer->Add(new TH3F(Form("hmiEtaPhiPHOS_%s",cut[iCut]),"Eta-phi-E correlations",100,-1.2,1.2,100,-TMath::Pi()/2.,3.*TMath::Pi()/2.,100,0.,5.));
    
  }
  
  
  
  
  //---PID cuts
	// Array of current cut's gammas
	fCutFolder = new TList*[fnCuts];


// 	fV0Reader=(AliV0ReaderV1*)AliAnalysisManager::GetAnalysisManager()->GetTask("V0ReaderV1");
// 	if(!fV0Reader){printf("Error: No V0 Reader");return;} // GetV0Reader
/*	
	if(fV0Reader)
		if((AliConvEventCuts*)fV0Reader->GetEventCuts())
			if(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms())
				fOutputContainer->Add(((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetCutHistograms());
	
	if(fV0Reader)
		if((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())
			if(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms())
				fOutputContainer->Add(((AliConversionPhotonCuts*)fV0Reader->GetConversionCuts())->GetCutHistograms());
	
	for(Int_t iCut = 0; iCut<fnCuts;iCut++){
		if(!((AliConvEventCuts*)fEventCutArray->At(iCut))) continue;
		if(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms()){
			fCutFolder[iCut]->Add(((AliConvEventCuts*)fEventCutArray->At(iCut))->GetCutHistograms());
		}
		if(!((AliConversionPhotonCuts*)fCutArray->At(iCut))) continue;
		if(((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutHistograms()){
			fCutFolder[iCut]->Add(((AliConversionPhotonCuts*)fCutArray->At(iCut))->GetCutHistograms());
		}
		if(fDoMesonAnalysis){
			if(!((AliConversionMesonCuts*)fMesonCutArray->At(iCut))) continue;
			if(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms()){
				fCutFolder[iCut]->Add(((AliConversionMesonCuts*)fMesonCutArray->At(iCut))->GetCutHistograms());
			}
		}
	}
*/


  PostData(1, fOutputContainer);

}

//________________________________________________________________________
void AliAnalysisTaskEtaPhiMultgg::UserExec(Option_t *) 
{
  // Main loop, called for each event
  // Analyze ESD/AOD
    
  FillHistogram("hTotSelEvents",0.5) ;
  
 
//   Int_t eventQuality = ((AliConvEventCuts*)fV0Reader->GetEventCuts())->GetEventQuality();
//   if(eventQuality == 2 || eventQuality == 3){// Event Not Accepted due to MC event missing or wrong trigger for V0ReaderV1
// 	return;
//   }

  //No PCM particles
//   if(!fV0Reader->IsEventSelected()){
//     PostData(1, fOutputContainer);
//     return;  
//   }
    
  fEvent = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!fEvent) {
    Printf("ERROR: Could not retrieve event");
    PostData(1, fOutputContainer);
    return;
  }

  fRunNumber=ConvertRunNumber(fEvent->GetRunNumber()) ;
  FillHistogram("hSelEvents",1.5,fRunNumber-0.5) ;
  FillHistogram("hTotSelEvents",1.5) ;
  
 
  if(fEventCounter == 0) {
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
   
    
    fEMCALgeo = AliEMCALGeometry::GetInstance();
    
    
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

  //Whole V0
  fRP = eventPlane->CalculateVZEROEventPlane(fEvent,10, harmonics,qx,qy);

  while(rpV0A<0)rpV0A+=TMath::TwoPi()/harmonics ;
  while(rpV0A>TMath::TwoPi()/harmonics)rpV0A-=TMath::TwoPi()/harmonics ;
  
  while(rpV0C<0)rpV0C+=TMath::TwoPi()/harmonics ;
  while(rpV0C>TMath::TwoPi()/harmonics)rpV0C-=TMath::TwoPi()/harmonics ;

  while(fRP<0)fRP+=TMath::TwoPi()/harmonics ;
  while(fRP>TMath::TwoPi()/harmonics)fRP-=TMath::TwoPi()/harmonics ;
  
  
  FillHistogram("phiRP",fRP,fCentrality) ;  
  FillHistogram("phiRPV0A",rpV0A,fCentrality) ;  
  FillHistogram("phiRPV0C",rpV0C,fCentrality) ;  
  FillHistogram("phiRPV0AC",rpV0A,rpV0C) ;  

  rpV0A = fV0AFlat->MakeFlat(rpV0A,fCentrality) ;
  rpV0C = fV0CFlat->MakeFlat(rpV0C,fCentrality) ;
  FillHistogram("phiRPV0Aflat",rpV0A,fCentrality) ;  
  FillHistogram("phiRPV0Cflat",rpV0C,fCentrality) ;  
  FillHistogram("phiRPV0ACflat",rpV0A,rpV0C) ;  
  
  FillHistogram("phiRPvsV0A",fRP,rpV0A) ;  
  FillHistogram("phiRPvsV0C",fRP,rpV0C) ;  
  
 
  FillHistogram("hSelEvents",5.5,fRunNumber-0.5) ;
  FillHistogram("hTotSelEvents",5.5) ;
  //All event selections done
  FillHistogram("hCentrality",fCentrality,fRunNumber-0.5) ;
  //Reaction plane is defined in the range (0;pi)
  //We have 10 bins
    
  
  Int_t irp=Int_t(10.*(fRP)/TMath::Pi());
  if(irp>9)irp=9 ;

  if(!fPHOSEvents[zvtx][fCenBin][irp]) 
    fPHOSEvents[zvtx][fCenBin][irp]=new TList() ;
  TList * prevPHOS = fPHOSEvents[zvtx][fCenBin][irp] ;
  if(!fEMCALEvents[zvtx][fCenBin][irp]) 
    fEMCALEvents[zvtx][fCenBin][irp]=new TList() ;
  TList * prevEMCAL = fEMCALEvents[zvtx][fCenBin][irp] ;
  if(!fPCMEvents[zvtx][fCenBin][irp]) 
    fPCMEvents[zvtx][fCenBin][irp]=new TList() ;
  TList * prevPCM =  fPCMEvents[zvtx][fCenBin][irp] ;

//  printf("prevHadrEvent=%p, nPrevHadr=%d \n",prevHadrEvent,nPrevHadr) ;
  
  if(fPHOSEvent)
    fPHOSEvent->Clear() ;
  else
    fPHOSEvent = new TClonesArray("AliCaloPhoton",100) ;

  if(fEMCALEvent)
    fEMCALEvent->Clear() ;
  else
    fEMCALEvent = new TClonesArray("AliCaloPhoton",100) ;
  
  if(fGammaCandidates)
    fGammaCandidates->Clear() ;
  else
    fGammaCandidates = new TClonesArray("AliAODConversionPhoton",20) ;
  
  
  TVector3 vertex(vtx5);
  
  Int_t multClust = fEvent->GetNumberOfCaloClusters();
  Int_t inPHOS=0; 

  AliAODCaloCells * cells = fEvent->GetPHOSCells() ;
  FillHistogram("hCenPHOSCells",fCentrality,cells->GetNumberOfCells()) ;
  FillHistogram("hCenTrack",fCentrality,fEvent->GetNumberOfTracks()) ;
  
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
    TVector3 local ;
    fPHOSGeo->Global2Local(local,global,mod);
    
    //Remove 6 noisy channels in run 139036
    if(fEvent->GetRunNumber()==139036 && mod==1 && 
       (cellX==9||cellX==10||cellX==11) && (cellZ==45 || cellZ==46))
      continue ;
    

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
    
    if(inPHOS>=fPHOSEvent->GetSize()){
      fPHOSEvent->Expand(inPHOS+20) ;
    }
    new((*fPHOSEvent)[inPHOS]) AliCaloPhoton(pv1.X(),pv1.Py(),pv1.Z(),pv1.E()) ;
    AliCaloPhoton * ph = (AliCaloPhoton*)fPHOSEvent->At(inPHOS) ;
    ph->SetModule(mod) ;
    pv1*= clu->GetCoreEnergy()/pv1.E() ;
    ph->SetMomV2(&pv1) ;
    ph->SetNCells(clu->GetNCells());    
    ph->SetDispBit(clu->Chi2()<2.5*2.5) ;
    ph->SetDisp2Bit(clu->Chi2()<1.5*1.5) ;

//    Double_t distBC=clu->GetDistanceToBadChannel();
    if(ph->IsDispOK()){
      FillHistogram(Form("hCluDispM%d",mod),cellX,cellZ,1.);
    }
    ph->SetCPVBit(clu->GetEmcCpvDistance()>2.5) ;
    if(ph->IsCPVOK()){
      FillHistogram(Form("hCluVetoM%d",mod),cellX,cellZ,1.);
    }
    ph->SetCPV2Bit(clu->GetEmcCpvDistance()>4.) ;
    
    ph->SetPrimary(clu->GetLabelAt(0)) ;
    ph->SetEMCx(local.X()) ;
    ph->SetEMCz(local.Z()) ;
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


  
  //Cuts
  const Double_t kEtaHBTcut=0.02; //Minimal eta gap to exclude HBT correlations
  const Double_t kEtaPCMcut=0.05; //Minimal eta gap to exclude HBT correlations
  const Double_t distCut1=5. ;    //First (minimal) distance cut
  const Double_t distCut2=10. ;    //Second (medium) distance cut
  const Double_t distCut3=20. ;    //Third (strongest) distance cut
  const Double_t asymCut=0.2 ;    //Asymetry cut
  
    
  const Double_t kgMass=0. ;
	
  const Int_t nCuts=7 ;
  char cut[7][20] ;
  sprintf(cut[0],"All") ;
  sprintf(cut[1],"Disp") ;
  sprintf(cut[2],"CPV") ;
  sprintf(cut[3],"Both") ;
  sprintf(cut[4],"Disp2") ;
  sprintf(cut[5],"CPV2") ;
  sprintf(cut[6],"Both2") ;

  //Real
  //PHOS-PHOS
  for (Int_t i1=0; i1<inPHOS-3; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1) ;
    Double_t e1=ph1->E() ;
    for (Int_t i2=i1+1; i2<inPHOS-2; i2++) {
      AliCaloPhoton * ph2=(AliCaloPhoton*)fPHOSEvent->At(i2) ;
      Double_t e2=ph2->E() ;
      if(!IsSameKtBin(e1,e2))
	continue ;
      for (Int_t i3=i2+1; i3<inPHOS-1; i3++) {
        AliCaloPhoton * ph3=(AliCaloPhoton*)fPHOSEvent->At(i3) ;
        Double_t e3=ph3->E() ;
        if(!IsSameKtBin(e1,e3))
	  continue ;
        for (Int_t i4=i3+1; i4<inPHOS; i4++) {
          AliCaloPhoton * ph4=(AliCaloPhoton*)fPHOSEvent->At(i4) ;
          Double_t e4=ph4->E() ;
          if(!IsSameKtBin(e1,e4))
	    continue ;

	  Double_t dEta1= ph1->Eta()+ph2->Eta()-ph3->Eta()-ph4->Eta() ;
          Double_t dPhi1= ph1->Phi()+ph2->Phi()-ph3->Phi()-ph4->Phi() ;
	  Double_t dEta2= ph1->Eta()-ph2->Eta()+ph3->Eta()-ph4->Eta() ;
          Double_t dPhi2= ph1->Phi()-ph2->Phi()+ph3->Phi()-ph4->Phi() ;
	  Double_t dEta3= ph1->Eta()-ph2->Eta()-ph3->Eta()+ph4->Eta() ;
          Double_t dPhi3= ph1->Phi()-ph2->Phi()-ph3->Phi()+ph4->Phi() ;
        
          if(gRandom->Uniform()>0.5){
	    dPhi1=-dPhi1 ;
	    dEta1=-dEta1 ;
	    dPhi2=-dPhi2 ;
	    dEta2=-dEta2 ;
	    dPhi3=-dPhi3 ;
	    dEta3=-dEta3 ;
          }
          for(Int_t iCut=0; iCut<7; iCut++){
	    if(!(PairCut(ph1,ph2,iCut) && PairCut(ph3,ph4,iCut)))
	      continue ;	
            FillHistogram(Form("hEtaPhiPHOS_%s",cut[iCut]),dEta1,dPhi1,e1) ;
            FillHistogram(Form("hEtaPhiPHOS_%s",cut[iCut]),dEta2,dPhi2,e1) ;
            FillHistogram(Form("hEtaPhiPHOS_%s",cut[iCut]),dEta3,dPhi3,e1) ;
	  }
	}
      }
    }
  }
  
  //now mixed
  //mixed-PHOS-PHOS
  TClonesArray * mixPHOS1 = static_cast<TClonesArray*>(prevPHOS->At(0)) ;
  TClonesArray * mixPHOS2 = static_cast<TClonesArray*>(prevPHOS->At(1)) ;
  TClonesArray * mixPHOS3 = static_cast<TClonesArray*>(prevPHOS->At(2)) ;
  if(mixPHOS1 && mixPHOS2 && mixPHOS3){
    for (Int_t i1=0; i1<inPHOS; i1++) {
      AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1) ;
      Double_t e1=ph1->E() ;
      for(Int_t i2=0; i2<mixPHOS1->GetEntriesFast();i2++){
        AliCaloPhoton * ph2=(AliCaloPhoton*)mixPHOS1->At(i2) ;
        Double_t e2=ph2->E() ;
        if(!IsSameKtBin(e1,e2))
        for(Int_t i3=0; i3<mixPHOS2->GetEntriesFast();i3++){
          AliCaloPhoton * ph3=(AliCaloPhoton*)mixPHOS2->At(i3) ;
          Double_t e3=ph3->E() ;
          if(!IsSameKtBin(e1,e3))
	  continue ;
          for(Int_t i4=0; i4<mixPHOS1->GetEntriesFast();i4++){
            AliCaloPhoton * ph4=(AliCaloPhoton*)mixPHOS3->At(i4) ;
            Double_t e4=ph4->E() ;
            if(!IsSameKtBin(e1,e4))
	      continue ;
	    
	    Double_t dEta1= ph1->Eta()+ph2->Eta()-ph3->Eta()-ph4->Eta() ;
            Double_t dPhi1= ph1->Phi()+ph2->Phi()-ph3->Phi()-ph4->Phi() ;
	    Double_t dEta2= ph1->Eta()-ph2->Eta()+ph3->Eta()-ph4->Eta() ;
            Double_t dPhi2= ph1->Phi()-ph2->Phi()+ph3->Phi()-ph4->Phi() ;
	    Double_t dEta3= ph1->Eta()-ph2->Eta()-ph3->Eta()+ph4->Eta() ;
            Double_t dPhi3= ph1->Phi()-ph2->Phi()-ph3->Phi()+ph4->Phi() ;
        
            if(gRandom->Uniform()>0.5){
	      dPhi1=-dPhi1 ;
	      dEta1=-dEta1 ;
	      dPhi2=-dPhi2 ;
	      dEta2=-dEta2 ;
	      dPhi3=-dPhi3 ;
	      dEta3=-dEta3 ;
            }
            for(Int_t iCut=0; iCut<7; iCut++){
	      if(!(PairCut(ph1,ph2,iCut) && PairCut(ph3,ph4,iCut)))
	        continue ;	
              FillHistogram(Form("hmiEtaPhiPHOS_%s",cut[iCut]),dEta1,dPhi1,e1) ;
              FillHistogram(Form("hmiEtaPhiPHOS_%s",cut[iCut]),dEta2,dPhi2,e1) ;
              FillHistogram(Form("hmiEtaPhiPHOS_%s",cut[iCut]),dEta3,dPhi3,e1) ;
	    }
	  }
	}
      }
    }
  }
 
  //Now we either add current events to stack or remove
  //If no photons in current event - no need to add it to mixed
  const Int_t kMixEvents[6]={3,3,3,3,3,3} ;
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
void AliAnalysisTaskEtaPhiMultgg::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
  
}
//________________________________________________________________________
Bool_t AliAnalysisTaskEtaPhiMultgg::IsSameKtBin(Double_t e1, Double_t e2){
  const Double_t binWidth=0.2;
  return (int(e1/binWidth)==int(e2/binWidth)) ;
  
}
  
  

