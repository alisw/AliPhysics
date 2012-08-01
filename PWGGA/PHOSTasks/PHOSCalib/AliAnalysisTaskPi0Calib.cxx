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

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskPi0Calib.h"
#include "AliCaloPhotonC.h"
#include "AliPHOSGeometry.h"
#include "AliPHOSEsdCluster.h"
#include "AliPHOSCalibData.h"
#include "AliESDEvent.h"
#include "AliESDCaloCells.h"
#include "AliESDVertex.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliCDBManager.h"

// Analysis task to fill histograms with PHOS ESD clusters and cells
// Authors: Dmitri Peressounko
// Date   : 28.05.2011

ClassImp(AliAnalysisTaskPi0Calib)

//________________________________________________________________________
Double_t rnlin(Double_t *x, Double_t *par)
{
  //a = par[0], b = par[1].
  //1+a*exp(-e/b)

 return 0.0241+1.0504*x[0]+0.000249*x[0]*x[0] ; 

  Double_t gcbCorr = 0.0241+1.0504*x[0]+0.000249*x[0]*x[0] ;
//  Double_t bvpCorr = 1+par[0]*TMath::Exp(-x[0]/par[1]);
//  Double_t bvpCorr = 1.005*(1+par[0]*par[1]*par[1]/(x[0]*x[0]+par[1]*par[1])) ;
  Double_t bvpCorr = 1.007*(1+par[0]*par[1]*par[1]/(x[0]*x[0]+par[1]*par[1])) ;
  return gcbCorr*bvpCorr;

}
//________________________________________________________________________
Double_t nonlin(Double_t e){
 return 0.0241+1.0504*e+0.000249*e*e ; 
}

//________________________________________________________________________
AliAnalysisTaskPi0Calib::AliAnalysisTaskPi0Calib(const char *name) 
  : AliAnalysisTaskSE(name),
    fOutputContainer(0),
    fPHOSEvents(0),
    fPHOSEvent(0),
    fPHOSCalibData(0),
    fPHOSGeo(0),
    fEventCounter(0)
{
  // Constructor
  
  // Output slots #0 write into a TH1 container
  DefineOutput(1,TList::Class());

  // Set bad channel map
  char key[55] ;
  for(Int_t i=0; i<6; i++){
    sprintf(key,"PHOS_BadMap_mod%d",i) ;
    fPHOSBadMap[i]=new TH2I(key,"Bad Modules map",64,0.,64.,56,0.,56.) ;
    fPHOSAmp[i]=0x0 ;
  }
  // Initialize the PHOS geometry
  fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP") ;

  //We have to apply re-calibration for pass1 LCH10h
  // Initialize decalibration factors in the form of the OCDB object


  AliCDBManager * man = AliCDBManager::Instance();
  man->SetRun(140000) ;
  man->SetDefaultStorage("local://OCDB");
  fPHOSCalibData = new AliPHOSCalibData();


  // Initialize non-linrarity correction
//  fNonLinCorr = new TF1("nonlib",rnlin,0.,40.,0);

  fOutputContainer = new TList();
  fOutputContainer->SetOwner(kTRUE);

}
//________________________________________________________________________
AliAnalysisTaskPi0Calib::~AliAnalysisTaskPi0Calib(){

  for(Int_t i=0; i<6; i++){
    if(fPHOSBadMap[i]){
      delete fPHOSBadMap[i] ;
      fPHOSBadMap[i]=0x0 ;
    }
    if(fPHOSAmp[i]){
      delete fPHOSAmp[i] ;
      fPHOSAmp[i]=0x0 ;
    }
  }
}

//________________________________________________________________________
void AliAnalysisTaskPi0Calib::UserCreateOutputObjects()
{
  // Create histograms

//  AliCDBManager * man = AliCDBManager::Instance();
//  man->SetRun(140000) ;
//  man->SetDefaultStorage("local://OCDB");
//  fPHOSCalibData = new AliPHOSCalibData();

//  fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP") ;



  // Called once
  const Int_t nRuns=200 ;
  
  // ESD histograms
//  if(fOutputContainer != NULL){
//    delete fOutputContainer;
//  }
//  fOutputContainer = new TList();
//  fOutputContainer->SetOwner(kTRUE);
  
  //========QA histograms=======

  //Event selection
  fOutputContainer->Add(new TH1F("hSelEvents","Event selection", 10,0.,10.)) ;
  
  //vertex distribution
  fOutputContainer->Add(new TH2F("hZvertex","Z vertex position", 50,-25.,25.,nRuns,0.,float(nRuns))) ;  
  
  fOutputContainer->Add(new TH3F("hMggAllM1","M1" ,64,0.5,64.5, 56,0.5,56.5,200,0.,0.5));
  fOutputContainer->Add(new TH3F("hMggAllM2","M2" ,64,0.5,64.5, 56,0.5,56.5,200,0.,0.5));
  fOutputContainer->Add(new TH3F("hMggAllM3","M3" ,64,0.5,64.5, 56,0.5,56.5,200,0.,0.5));

  fOutputContainer->Add(new TH3F("hMggDispM1","M1" ,64,0.5,64.5, 56,0.5,56.5,200,0.,0.5));
  fOutputContainer->Add(new TH3F("hMggDispM2","M2" ,64,0.5,64.5, 56,0.5,56.5,200,0.,0.5));
  fOutputContainer->Add(new TH3F("hMggDispM3","M3" ,64,0.5,64.5, 56,0.5,56.5,200,0.,0.5));

  fOutputContainer->Add(new TH3F("hMggBothM1","M1" ,64,0.5,64.5, 56,0.5,56.5,200,0.,0.5));
  fOutputContainer->Add(new TH3F("hMggBothM2","M2" ,64,0.5,64.5, 56,0.5,56.5,200,0.,0.5));
  fOutputContainer->Add(new TH3F("hMggBothM3","M3" ,64,0.5,64.5, 56,0.5,56.5,200,0.,0.5));

  fOutputContainer->Add(new TH3F("hMiMggAllM1","M1" ,64,0.5,64.5, 56,0.5,56.5,200,0.,0.5));
  fOutputContainer->Add(new TH3F("hMiMggAllM2","M2" ,64,0.5,64.5, 56,0.5,56.5,200,0.,0.5));
  fOutputContainer->Add(new TH3F("hMiMggAllM3","M3" ,64,0.5,64.5, 56,0.5,56.5,200,0.,0.5));

  fOutputContainer->Add(new TH3F("hMiMggDispM1","M1" ,64,0.5,64.5, 56,0.5,56.5,200,0.,0.5));
  fOutputContainer->Add(new TH3F("hMiMggDispM2","M2" ,64,0.5,64.5, 56,0.5,56.5,200,0.,0.5));
  fOutputContainer->Add(new TH3F("hMiMggDispM3","M3" ,64,0.5,64.5, 56,0.5,56.5,200,0.,0.5));

  fOutputContainer->Add(new TH3F("hMggOldM1","M1" ,64,0.5,64.5, 56,0.5,56.5,200,0.,0.5));
  fOutputContainer->Add(new TH3F("hMggOldM2","M2" ,64,0.5,64.5, 56,0.5,56.5,200,0.,0.5));
  fOutputContainer->Add(new TH3F("hMggOldM3","M3" ,64,0.5,64.5, 56,0.5,56.5,200,0.,0.5));

  fOutputContainer->Add(new TH3F("hMggTOFM1","M1" ,64,0.5,64.5, 56,0.5,56.5,200,0.,0.5));
  fOutputContainer->Add(new TH3F("hMggTOFM2","M2" ,64,0.5,64.5, 56,0.5,56.5,200,0.,0.5));
  fOutputContainer->Add(new TH3F("hMggTOFM3","M3" ,64,0.5,64.5, 56,0.5,56.5,200,0.,0.5));

 
  Int_t nPtPhot = 300 ;
  Double_t ptPhotMax = 30 ;
  Int_t nM       = 500;
  Double_t mMin  = 0.0;
  Double_t mMax  = 1.0;
  fOutputContainer->Add(new TH2F("hPi0M11","All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hPi0M22","All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hPi0M33","All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hPi0M12","All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hPi0M23","All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hPi0M13","All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));

  fOutputContainer->Add(new TH2F("hMiPi0M11","All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hMiPi0M22","All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hMiPi0M33","All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hMiPi0M12","All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hMiPi0M23","All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  fOutputContainer->Add(new TH2F("hMiPi0M13","All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));

  fOutputContainer->Add(new TH3F("hTOFM1","M1" ,64,0.5,64.5, 56,0.5,56.5,200,-10.e-6,10.e-6));
  fOutputContainer->Add(new TH3F("hTOFM2","M2" ,64,0.5,64.5, 56,0.5,56.5,200,-10.e-6,10.e-6));
  fOutputContainer->Add(new TH3F("hTOFM3","M3" ,64,0.5,64.5, 56,0.5,56.5,200,-10.e-6,10.e-6));
  fOutputContainer->Add(new TH3F("hTOFFineM1","M1" ,64,0.5,64.5, 56,0.5,56.5,200,-100.e-9,100.e-9));
  fOutputContainer->Add(new TH3F("hTOFFineM2","M2" ,64,0.5,64.5, 56,0.5,56.5,200,-100.e-9,100.e-9));
  fOutputContainer->Add(new TH3F("hTOFFineM3","M3" ,64,0.5,64.5, 56,0.5,56.5,200,-100.e-9,100.e-9));

  fOutputContainer->Add(new TH2F("hTOFvsEM1","All clusters",200,-100.e-9,100.e-9,100,0.,20.));
  fOutputContainer->Add(new TH2F("hTOFvsEM2","All clusters",200,-100.e-9,100.e-9,100,0.,20.));
  fOutputContainer->Add(new TH2F("hTOFvsEM3","All clusters",200,-100.e-9,100.e-9,100,0.,20.));

  fOutputContainer->Add(new TH2F("hTOFAvvsEM1","All clusters",200,-100.e-9,100.e-9,100,0.,20.));
  fOutputContainer->Add(new TH2F("hTOFAvvsEM2","All clusters",200,-100.e-9,100.e-9,100,0.,20.));
  fOutputContainer->Add(new TH2F("hTOFAvvsEM3","All clusters",200,-100.e-9,100.e-9,100,0.,20.));
  
  fOutputContainer->Add(new TH2F("hTOFvsECluM1","All clusters",200,-100.e-9,100.e-9,100,0.,20.));
  fOutputContainer->Add(new TH2F("hTOFvsECluM2","All clusters",200,-100.e-9,100.e-9,100,0.,20.));
  fOutputContainer->Add(new TH2F("hTOFvsECluM3","All clusters",200,-100.e-9,100.e-9,100,0.,20.));

  
  PostData(1, fOutputContainer);

}

//________________________________________________________________________
void AliAnalysisTaskPi0Calib::UserExec(Option_t *) 
{
  // Main loop, called for each event
  // Analyze ESD/AOD
    
  AliESDEvent *event = dynamic_cast<AliESDEvent*>(InputEvent());
  if (!event) {
    Printf("ERROR: Could not retrieve event");
    PostData(1, fOutputContainer);
    return;
  }
  FillHistogram("hSelEvents",0.5) ;
  
  // Get PHOS rotation matrices from ESD and set them to the PHOS geometry
//  char skey[55] ;
  if(fEventCounter == 0) {
    for(Int_t mod=0; mod<5; mod++) {
      if(!event->GetPHOSMatrix(mod)) continue;
      fPHOSGeo->SetMisalMatrix(event->GetPHOSMatrix(mod),mod) ;
      Printf("PHOS geo matrix %p for module # %d is set\n", event->GetPHOSMatrix(mod), mod);
    }
    fEventCounter++ ;
  }

  if(fPHOSCalibData->GetADCchannelEmc(3,15,26)==1.){
    for(Int_t module=1; module<4; module++) {
      for(Int_t ix=1; ix<=64; ix++) {
        for(Int_t iz=1; iz<=56; iz++) {

	  Float_t corr = fPHOSAmp[module]->GetBinContent(ix,iz) ;
 	    fPHOSCalibData->SetADCchannelEmc(module,iz,ix,corr);
	}
      }
    }
  }
  
  
  // Checks if we have a primary vertex
  // Get primary vertices form ESD
  const AliESDVertex *esdVertex5 = event->GetPrimaryVertex();

  Double_t vtx0[3] = {0,0,0}; // don't rely on ESD vertex, assume (0,0,0)
  Double_t vtx5[3] ={0.,0.,0.};
  
  vtx5[0] = esdVertex5->GetX();
  vtx5[1] = esdVertex5->GetY();
  vtx5[2] = esdVertex5->GetZ();
  
  
  FillHistogram("hZvertex",esdVertex5->GetZ());
/*
  if (TMath::Abs(esdVertex5->GetZ()) > 10. ){
    PostData(1, fOutputContainer);
    return;
  }
*/
  FillHistogram("hSelEvents",1.5) ;

  if(!fPHOSEvents) 
    fPHOSEvents=new TList() ;
  TList * prevPHOS = fPHOSEvents ;

  if(fPHOSEvent)
    fPHOSEvent->Clear() ;
  else
    fPHOSEvent = new TClonesArray("AliCaloPhotonC",200) ;

  //For re-calibration
  const Double_t logWeight=4.5 ;  

  TVector3 vertex(vtx5);
  
  Int_t multClust = event->GetNumberOfCaloClusters();
  Int_t inPHOS=0; //,inMod1=0, inMod2=0, inMod3=0 ;
//   Double_t avrgEm1=0,avrgEm2=0,avrgEm3=0; //average cluster energy

  AliESDCaloCells * cells = event->GetPHOSCells() ;
//   char key[55];
  
  TVector3 globalPos; 
  TVector3 localPos ;
  TLorentzVector p12 ;  
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
//    if ( !IsGoodChannel("PHOS",mod,cellX,cellZ) ) 
//      continue ;
    
    
    //Apply re-Calibreation
    AliPHOSEsdCluster cluPHOS1(*clu);
    cluPHOS1.Recalibrate(fPHOSCalibData,cells); // modify the cell energies
    cluPHOS1.EvalAll(logWeight,vertex);         // recalculate the cluster parameters
    cluPHOS1.SetE(nonlin(cluPHOS1.E()));// Users's nonlinearity
    if(cluPHOS1.E()<0.3) continue;
    if(clu->GetNCells()<3) continue ;
    if(clu->GetM02()<0.2)   continue ;
 
    //TOF recalibration
    Double_t dt=0.;
    Double_t tof = clu->GetTOF()-dt ;
    
    FillHistogram(Form("hTOFM%d",mod),float(cellX),float(cellZ),tof);
    if(clu->E()>1.){
      FillHistogram(Form("hTOFFineM%d",mod),float(cellX),float(cellZ),tof);
    }
    FillHistogram(Form("hTOFvsEM%d",mod),tof,clu->E());      
      
/*
    //Calculate TOF within cluster
    Double_t tofAv = 0. ;
    Double_t wtot=0. ;
    for(Int_t iDigit=0; iDigit<clu->GetNCells(); iDigit++) {
      Int_t absId=clu->GetCellAbsId(iDigit);
      Double_t edig = cells->GetCellAmplitude(absId) ;
      Double_t tau = cells->GetCellTime(absId) ;
      Int_t relid[4] ;
      fPHOSGeo->AbsToRelNumbering(clu->GetCellAbsId(iDigit), relid) ;
      if(fPHOSTOF[relid[0]]){
        tau -= fPHOSTOF[relid[0]]->GetBinContent(relid[2],relid[3]) ;
      }
      sprintf(key,"hTOFvsECluM%d",mod) ;
      FillHistogram(key,tau-tof,edig);
    
      if(edig>0.4 && TMath::Abs(tau)<10.e-9){
        Double_t wi=2.73+11.8*TMath::Exp(-edig/0.45) ;
        tofAv+=tau/wi/wi ;
        wtot+=1./wi/wi ;
      }
    }
    if(wtot>0.)
      tofAv/=wtot ;
    else
      tofAv=tof ;
    FillHistogram(Form("hTOFAvvsEM%d",mod),tofAv,clu->E());
*/
    
    
    TLorentzVector pv1 ;
    cluPHOS1.GetMomentum(pv1,vtx0);
    if(inPHOS>=fPHOSEvent->GetSize()){
      fPHOSEvent->Expand(inPHOS+50) ;
    }
    new((*fPHOSEvent)[inPHOS]) AliCaloPhotonC(pv1.X(),pv1.Py(),pv1.Z(),pv1.E()) ;
    AliCaloPhotonC * ph = (AliCaloPhotonC*)fPHOSEvent->At(inPHOS) ;
    ph->SetModule(mod) ;

    clu->GetMomentum(pv1 ,vtx0);
    ph->SetMomV2(&pv1) ;
    ph->SetNCells(clu->GetNCells());
    ph->SetDispBit(TestLambda(clu->E(),clu->GetM20(),clu->GetM02())) ;
    ph->SetTOFBit(TestTOF(clu->E(),tof)) ;
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
      }
    }
    ph->SetCPVBit(cpvBit) ;
    ph->SetCPVBit2(cpvBit2) ;
    ph->SetEMCx(float(cellX)) ;
    ph->SetEMCz(float(cellZ)) ;
    ph->SetLambdas(clu->GetM20(),clu->GetM02()) ;          
    inPHOS++ ;
  }
  
  
  TLorentzVector p12old ;
  char key[255] ;
  for (Int_t i1=0; i1<inPHOS-1; i1++) {
    AliCaloPhotonC * ph1=(AliCaloPhotonC*)fPHOSEvent->At(i1) ;
    for (Int_t i2=i1+1; i2<inPHOS; i2++) {
      AliCaloPhotonC * ph2=(AliCaloPhotonC*)fPHOSEvent->At(i2) ;
      p12  = *ph1  + *ph2;
      if(p12.Pt()<2.)
	continue ;
      p12old = *(ph1->GetMomV2()) + *(ph2->GetMomV2()) ;

          if(ph1->Module()==1 && ph2->Module()==1)
	    FillHistogram("hPi0M11",p12.M(),p12.Pt() );
          else if(ph1->Module()==2 && ph2->Module()==2)
	    FillHistogram("hPi0M22",p12.M(),p12.Pt() );
          else if(ph1->Module()==3 && ph2->Module()==3)
	    FillHistogram("hPi0M33",p12.M(),p12.Pt() );
          else if(ph1->Module()==1 && ph2->Module()==2)
	    FillHistogram("hPi0M12",p12.M(),p12.Pt() );
          else if(ph1->Module()==1 && ph2->Module()==3)
	    FillHistogram("hPi0M13",p12.M(),p12.Pt() );
          else if(ph1->Module()==2 && ph2->Module()==3)
	    FillHistogram("hPi0M23",p12.M(),p12.Pt() );
      
      if(ph1->Module()!=ph2->Module())
	continue ;
      
      FillHistogram(Form("hMggAllM%d",ph1->Module()),ph1->EMCx(),ph1->EMCz(),p12.M() );
      FillHistogram(Form("hMggAllM%d",ph2->Module()),ph2->EMCx(),ph2->EMCz(),p12.M() );

      if(ph1->IsTOFOK() && ph2->IsTOFOK()){
        sprintf(key,"hMggTOFM%d",ph1->Module()) ;
        FillHistogram(Form("hMggTOFM%d",ph1->Module()),ph1->EMCx(),ph1->EMCz(),p12.M() );
        FillHistogram(Form("hMggTOFM%d",ph1->Module()),ph2->EMCx(),ph2->EMCz(),p12.M() );
      }

      if(ph1->IsDispOK() && ph2->IsDispOK()){
	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  sprintf(key,"hMggDispM%d",ph1->Module()) ;
	  FillHistogram(key,ph1->EMCx(),ph1->EMCz(),p12.M() );
	  sprintf(key,"hMggDispM%d",ph2->Module()) ;
	  FillHistogram(key,ph2->EMCx(),ph2->EMCz(),p12.M() );
          if(ph1->IsTOFOK() && ph2->IsTOFOK()){
            sprintf(key,"hMggBothM%d",ph1->Module()) ;
            FillHistogram(key,ph1->EMCx(),ph1->EMCz(),p12.M() );
            sprintf(key,"hMggBothM%d",ph2->Module()) ;
            FillHistogram(key,ph2->EMCx(),ph2->EMCz(),p12.M() );
          }
	  sprintf(key,"hMggOldM%d",ph1->Module()) ;
	  FillHistogram(key,ph1->EMCx(),ph1->EMCz(),p12old.M() );
	  sprintf(key,"hMggOldM%d",ph2->Module()) ;
	  FillHistogram(key,ph2->EMCx(),ph2->EMCz(),p12old.M() );
	}
      }
      
    } // end of loop i2
  } // end of loop i1

  
  //now mixed
  for (Int_t i1=0; i1<inPHOS; i1++) {
    AliCaloPhotonC * ph1=(AliCaloPhotonC*)fPHOSEvent->At(i1) ;
    for(Int_t ev=0; ev<prevPHOS->GetSize();ev++){
      TClonesArray * mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev)) ;
      for(Int_t i2=0; i2<mixPHOS->GetEntriesFast();i2++){
	AliCaloPhotonC * ph2=(AliCaloPhotonC*)mixPHOS->At(i2) ;
	p12  = *ph1  + *ph2;
        if(p12.Pt()<2.)
  	  continue ;

      if(ph1->Module()!=ph2->Module())
	continue ;

            if(ph1->Module()==1 && ph2->Module()==1)
   	      FillHistogram("hMiPi0M11",p12.M(),p12.Pt() );
            else if(ph1->Module()==2 && ph2->Module()==2)
	      FillHistogram("hMiPi0M22",p12.M(),p12.Pt() );
            else if(ph1->Module()==3 && ph2->Module()==3)
	      FillHistogram("hMiPi0M33",p12.M(),p12.Pt() );
            else if(ph1->Module()==1 && ph2->Module()==2)
	      FillHistogram("hMiPi0M12",p12.M(),p12.Pt() );
            else if(ph1->Module()==1 && ph2->Module()==3)
	      FillHistogram("hMiPi0M13",p12.M(),p12.Pt() );
            else if(ph1->Module()==2 && ph2->Module()==3)
	      FillHistogram("hMiPi0M23",p12.M(),p12.Pt() );
	
	
	sprintf(key,"hMiMggAllM%d",ph1->Module()) ;
	FillHistogram(key,ph1->EMCx(),ph1->EMCz(),p12.M() );
	FillHistogram(key,ph2->EMCx(),ph2->EMCz(),p12.M() );
	
	if(ph1->IsDispOK() && ph2->IsDispOK()){
	  if(ph1->IsCPVOK() && ph2->IsCPVOK()){
  	    sprintf(key,"hMiMggDispM%d",ph1->Module()) ;
	    FillHistogram(key,ph1->EMCx(),ph1->EMCz(),p12.M() );
	    FillHistogram(key,ph2->EMCx(),ph2->EMCz(),p12.M() );
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
void AliAnalysisTaskPi0Calib::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
  
}

//________________________________________________________________________
Bool_t AliAnalysisTaskPi0Calib::IsGoodChannel(const char * det, Int_t mod, Int_t ix, Int_t iz)
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
void AliAnalysisTaskPi0Calib::FillHistogram(const char * key,Double_t x)const{
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
void AliAnalysisTaskPi0Calib::FillHistogram(const char * key,Double_t x,Double_t y)const{
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
void AliAnalysisTaskPi0Calib::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z) const{
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
Bool_t AliAnalysisTaskPi0Calib::TestLambda(Double_t pt,Double_t l1,Double_t l2){
  
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
Bool_t AliAnalysisTaskPi0Calib::TestTOF(Double_t /*e*/, Double_t tof){

  Double_t res = 100.e-9 ; //3.31714e-09  + 1.42519e-08*TMath::Exp(-e/4.61728e-01) ;
  Double_t mean=0. ;

  return TMath::Abs(tof-mean)<res ;

}
//____________________________________________________________________________
Double_t AliAnalysisTaskPi0Calib::TestCPV(Double_t dx, Double_t dz, Double_t pt, Int_t charge){
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
