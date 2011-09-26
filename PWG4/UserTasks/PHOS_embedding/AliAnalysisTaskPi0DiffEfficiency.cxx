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
#include "AliPHOSEsdCluster.h"
#include "AliPHOSCalibData.h"
#include "AliAODEvent.h"
#include "AliAODCaloCluster.h"
#include "AliAODVertex.h"
#include "AliESDtrackCuts.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliCDBManager.h"
#include "AliCentrality.h" 

// Analysis task to fill histograms with PHOS ESD clusters and cells
// Authors: Yuri Kharlov
// Date   : 28.05.2009

ClassImp(AliAnalysisTaskPi0DiffEfficiency)

//________________________________________________________________________
AliAnalysisTaskPi0DiffEfficiency::AliAnalysisTaskPi0DiffEfficiency(const char *name) 
: AliAnalysisTaskSE(name),
  fStack(0),
  fOutputContainer(0),
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
    snprintf(key,55,"hPhotCPV_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hPhotDisp_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hPhotBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;

    snprintf(key,55,"hOldMassPtAll_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hOldMassPtCPV_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hOldMassPtDisp_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hOldMassPtBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    
    snprintf(key,55,"hNewMassPtAll_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNewMassPtCPV_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNewMassPtDisp_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hNewMassPtBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    
    snprintf(key,55,"hMassPtAll_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMassPtCPV_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMassPtDisp_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMassPtBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    ((TH2F*)fOutputContainer->Last())->Sumw2() ;
    
    //Mixed
    snprintf(key,55,"hMiMassPtAll_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtCPV_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtDisp_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMiMassPtBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));

    snprintf(key,55,"hMCMassPtAll_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtCPV_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtDisp_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));
    snprintf(key,55,"hMCMassPtBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"(M,p_{T},d#phi)_{#gamma#gamma}" ,nM,mMin,mMax,nPt,ptMin,ptMax));

    //Single photon
    snprintf(key,55,"hMCPhotAll_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMCPhotCPV_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMCPhotDisp_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;
    snprintf(key,55,"hMCPhotBoth_cen%d",cent) ;
    fOutputContainer->Add(new TH1F(key,"dN/dpt" ,nPt,ptMin,ptMax));
    ((TH1F*)fOutputContainer->Last())->Sumw2() ;

    
    
  }


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
  Double_t vtx5[3] ={0.,0.,0.};
  
  vtx5[0] = esdVertex5->GetX();
  vtx5[1] = esdVertex5->GetY();
  vtx5[2] = esdVertex5->GetZ();
  
  
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
  Int_t irp=Int_t(10.*(fRPfull+TMath::PiOver2())/TMath::Pi());
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
  TClonesArray * clustersOld = event->GetCaloClusters() ;
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
    ph->SetMomV2(&pv1) ;
    ph->SetNCells(clu->GetNCells());
    ph->SetDispBit(TestLambda(clu->GetM20(),clu->GetM02())) ;
    ph->SetCPVBit(clu->GetEmcCpvDistance()>1.) ;
    if(!survive) //this cluster found in list after embedding, skipping it
     ph->SetTagged(1) ;

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
    ph->SetMomV2(&pv1) ;
    ph->SetNCells(clu->GetNCells());
    ph->SetDispBit(TestLambda(clu->GetM20(),clu->GetM02())) ;
    ph->SetCPVBit(clu->GetEmcCpvDistance()>1.) ;
    if(!survive) //this cluster found in list after embedding, skipping it
     ph->SetTagged(1) ;

    inPHOSemb++ ;
  }
  
  //Single photon
  for (Int_t i1=0; i1<inPHOSold; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent1->At(i1) ;
    if(!ph1->IsTagged())
      continue ;
    snprintf(key,55,"hPhotAll_cen%d",fCenBin) ;
    FillHistogram(key,ph1->Pt(),-1.) ;
    if(ph1->IsCPVOK() ){
      snprintf(key,55,"hPhotCPV_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt(),-1.) ;
    }
    if(ph1->IsDispOK()){
      snprintf(key,55,"hPhotDisp_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt(),-1.) ;
      if(ph1->IsCPVOK()){
	snprintf(key,55,"hPhotBoth_cen%d",fCenBin) ;
        FillHistogram(key,ph1->Pt(),-1.) ;
      }
    } // end of loop i2
  } // end of loop i1 

  for (Int_t i1=0; i1<inPHOSemb; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent2->At(i1) ;
    if(!ph1->IsTagged())
      continue ;
    snprintf(key,55,"hPhotAll_cen%d",fCenBin) ;
    FillHistogram(key,ph1->Pt(),1.) ;
    if(ph1->IsCPVOK() ){
      snprintf(key,55,"hPhotCPV_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt(),1.) ;
    }
    if(ph1->IsDispOK()){
      snprintf(key,55,"hPhotDisp_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt(),1.) ;
      if(ph1->IsCPVOK()){
	snprintf(key,55,"hPhotBoth_cen%d",fCenBin) ;
        FillHistogram(key,ph1->Pt(),1.) ;
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
      //Fill Controll histogram: Real before embedding
      snprintf(key,55,"hOldMassPtAll_cen%d",fCenBin) ;
      FillHistogram(key,p12.M() ,p12.Pt(),-1.) ;
      if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	snprintf(key,55,"hOldMassPtCPV_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt(),-1) ;
      }
      if(ph1->IsDispOK() && ph2->IsDispOK()){
	snprintf(key,55,"hOldMassPtDisp_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt(),-1) ;

	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  snprintf(key,55,"hOldMassPtBoth_cen%d",fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt(),-1) ;
        }
      }
      
      //Now fill mail histograms with negative contributions
      if(!(ph1->IsTagged() || ph2->IsTagged()) )
        continue ;
      snprintf(key,55,"hMassPtAll_cen%d",fCenBin) ;
      FillHistogram(key,p12.M() ,p12.Pt(),-1.) ;
      if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	snprintf(key,55,"hMassPtCPV_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt(),-1) ;
      }
      if(ph1->IsDispOK() && ph2->IsDispOK()){
	snprintf(key,55,"hMassPtDisp_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt(),-1) ;

	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  snprintf(key,55,"hMassPtBoth_cen%d",fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt(),-1) ;
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

      // Controll histogram: Real after embedding
      snprintf(key,55,"hNewMassPtAll_cen%d",fCenBin) ;
      FillHistogram(key,p12.M() ,p12.Pt(),1.) ;
      if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	snprintf(key,55,"hNewMassPtCPV_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt(),1) ;
      }
      if(ph1->IsDispOK() && ph2->IsDispOK()){
	snprintf(key,55,"hNewMassPtDisp_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt(),1) ;

	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  snprintf(key,55,"hNewMassPtBoth_cen%d",fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt(),1) ;
        }
      }
     
      //Now fill main histogamm
      //new clusters with positive contribution
      if(!(ph1->IsTagged() || ph2->IsTagged()) )
        continue ;
      snprintf(key,55,"hMassPtAll_cen%d",fCenBin) ;
      FillHistogram(key,p12.M() ,p12.Pt(),1.) ;
      if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	snprintf(key,55,"hMassPtCPV_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt(),1) ;
      }
      if(ph1->IsDispOK() && ph2->IsDispOK()){
	snprintf(key,55,"hMassPtDisp_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt(),1) ;

	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  snprintf(key,55,"hMassPtBoth_cen%d",fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt(),1) ;
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
	
	snprintf(key,55,"hMiMassPtAll_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt(),1.) ;
	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  snprintf(key,55,"hMiMassPtCPV_cen%d",fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt(),1.) ;
	}
	if(ph1->IsDispOK() && ph2->IsDispOK()){
	  snprintf(key,55,"hMiMassPtDisp_cen%d",fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt(),1.) ;
	  
	  if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	    snprintf(key,55,"hMiMassPtBoth_cen%d",fCenBin) ;
	    FillHistogram(key,p12.M() ,p12.Pt(),1.) ;
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

//________________________________________________________________________
void AliAnalysisTaskPi0DiffEfficiency::Terminate(Option_t *)
{
  // Draw result to the screen
  // Called once at the end of the query
  
}

//________________________________________________________________________
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
void AliAnalysisTaskPi0DiffEfficiency::FillHistogram(const char * key,Double_t x)const{
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
Bool_t AliAnalysisTaskPi0DiffEfficiency::TestLambda(Double_t l1,Double_t l2){
  Double_t l1Mean=1.22 ;
  Double_t l2Mean=2.0 ;
  Double_t l1Sigma=0.42 ;
  Double_t l2Sigma=0.71 ;
  Double_t c=-0.59 ;
  Double_t R2=(l1-l1Mean)*(l1-l1Mean)/l1Sigma/l1Sigma+(l2-l2Mean)*(l2-l2Mean)/l2Sigma/l2Sigma-c*(l1-l1Mean)*(l2-l2Mean)/l1Sigma/l2Sigma ;
  return (R2<9.) ;
  
}
//___________________________________________________________________________
void AliAnalysisTaskPi0DiffEfficiency::ProcessMC(){
  //fill histograms for efficiensy etc. calculation
  const Double_t rcut = 1. ; //cut for primary particles
  //---------First pi0/eta-----------------------------
  char partName[10] ;
  char hkey[55] ;

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
    //Total number of pi0 with creation radius <1 cm
    snprintf(hkey,55,"hMC_all_%s_cen%d",partName,fCenBin) ;
    FillHistogram(hkey,pt) ;
    if(TMath::Abs(particle->Y())<0.12){
      snprintf(hkey,55,"hMC_unitEta_%s_cen%d",partName,fCenBin) ;
      FillHistogram(hkey,pt) ;
    }

    snprintf(hkey,55,"hMC_rap_%s_cen%d",partName,fCenBin) ;
    FillHistogram(hkey,particle->Y()) ;
    
    Double_t phi=particle->Phi() ;
    while(phi<0.)phi+=TMath::TwoPi() ;
    while(phi>TMath::TwoPi())phi-=TMath::TwoPi() ;
    snprintf(hkey,55,"hMC_phi_%s_cen%d",partName,fCenBin) ;
    FillHistogram(hkey,phi) ;
   

    //Check if one of photons converted
    if(particle->GetNDaughters()!=2)
     continue ; //Do not account Dalitz decays

/*
    TParticle * gamma1 = fStack->Particle(particle->GetFirstDaughter());
    TParticle * gamma2 = fStack->Particle(particle->GetLastDaughter());
    //Number of pi0s decayed into acceptance
    Int_t mod1,mod2 ;
    Double_t x=0.,z=0. ;
    Bool_t hitPHOS1 = fPHOSGeo->ImpactOnEmc(gamma1, mod1, z,x) ;
    Bool_t hitPHOS2 = fPHOSGeo->ImpactOnEmc(gamma2, mod2, z,x) ;

    Bool_t goodPair=kFALSE ;
    if( hitPHOS1 && hitPHOS2){
      sprintf(hkey,"hMC_PHOSacc_%s",partName) ;
      FillHistogram(hkey,pt) ;
      goodPair=kTRUE ;
    }

*/
  }
 
  //Now calculate "Real" distribution of clusters with primary
  TClonesArray cluPrim("AliCaloPhoton",200) ; //clusters with primary
  Int_t multClust = event->GetNumberOfCaloClusters();
  Int_t inPHOS=0 ;
  Double_t vtx0[3] = {0,0,0}; 
  for (Int_t i=0; i<multClust; i++) {
    AliVCluster *clu = event->GetCaloCluster(i);
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
    
    if(inPHOS>=cluPrim.GetSize()){
      cluPrim.Expand(inPHOS+50) ;
    }
    AliCaloPhoton * ph = new(cluPrim[inPHOS]) AliCaloPhoton(pv1.X(),pv1.Py(),pv1.Z(),pv1.E()) ;
    //AliCaloPhoton * ph = (AliCaloPhoton*)fPHOSEvent->At(inPHOS) ;
    ph->SetModule(mod) ;
    ph->SetMomV2(&pv1) ;
    ph->SetNCells(clu->GetNCells());
    ph->SetDispBit(TestLambda(clu->GetM20(),clu->GetM02())) ;
    ph->SetCPVBit(clu->GetEmcCpvDistance()>1.) ;

    inPHOS++ ;

  }
  
  //Single photon
  char key[55] ;
  for (Int_t i1=0; i1<inPHOS; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)cluPrim.At(i1) ;
    snprintf(key,55,"hMCPhotAll_cen%d",fCenBin) ;
    FillHistogram(key,ph1->Pt()) ;
    if(ph1->IsCPVOK() ){
      snprintf(key,55,"hMCPhotCPV_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt()) ;
      
    }
    if(ph1->IsDispOK()){
      snprintf(key,55,"hMCPhotDisp_cen%d",fCenBin) ;
      FillHistogram(key,ph1->Pt()) ;
      if(ph1->IsCPVOK()){
	snprintf(key,55,"hMCPhotBoth_cen%d",fCenBin) ;
	FillHistogram(key,ph1->Pt()) ;
      }
    } // end of loop i2
  } // end of loop i1 

  // Fill Real disribution
  for (Int_t i1=0; i1<inPHOS-1; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)cluPrim.At(i1) ;
    for (Int_t i2=i1+1; i2<inPHOS; i2++) {
      AliCaloPhoton * ph2=(AliCaloPhoton*)cluPrim.At(i2) ;
      TLorentzVector p12  = *ph1  + *ph2;
      TLorentzVector pv12 = *(ph1->GetMomV2()) + *(ph2->GetMomV2());      
  
      snprintf(key,55,"hMCMassPtAll_cen%d",fCenBin) ;
      FillHistogram(key,p12.M() ,p12.Pt()) ;
      if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	snprintf(key,55,"hMCMassPtCPV_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt()) ;
      }
      if(ph1->IsDispOK() && ph2->IsDispOK()){
	snprintf(key,55,"hMCMassPtDisp_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt()) ;

	if(ph1->IsCPVOK() && ph2->IsCPVOK()){
	  snprintf(key,55,"hMCMassPtBoth_cen%d",fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt()) ;
        }
      }
    } // end of loop i2
  } // end of loop i1

  
}
//___________________________________________________________________________
Bool_t AliAnalysisTaskPi0DiffEfficiency::IsSameCluster(AliAODCaloCluster * c1, AliAODCaloCluster * c2)const{
 //Compare clusters before and after embedding
 //clusters are the same if 
 // - Energy changed less than 0.1%  (numerical accuracy in reconstruction)
 // - lists of digits are the same
  
 if(c1->GetNCells() != c2->GetNCells())
   return kFALSE ;
 
 if(TMath::Abs(c1->E()-c2->E())>0.001*c1->E())
   return kFALSE ;

 UShort_t *list1 = c1->GetCellsAbsId() ; 
 UShort_t *list2 = c2->GetCellsAbsId() ; 
 for(Int_t i=0; i<c1->GetNCells(); i++){
  if(list1[i] != list2[i])
    return kFALSE ;
 }
 return kTRUE ; 
  
}



