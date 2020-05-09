/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
 
// Analysis task for identified tracks matched to a PHOS cluster
// Authors: Dmitri Peresunko

#include "TChain.h"
#include "TObjArray.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TH3F.h"
#include "THashList.h"
#include "TGeoGlobalMagField.h"

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisPhotonDDA.h"
#include "AliPHOSGeometry.h"
#include "AliAODEvent.h"
#include "AliAODCaloCells.h"
#include "AliAODCaloCluster.h"
#include "AliAODVertex.h"
#include "AliLog.h"
#include "AliPID.h"
#include "AliAODInputHandler.h"
#include "AliPIDResponse.h"
#include "AliCentrality.h"
#include "AliMagF.h"
#include "AliCaloPhoton.h"
#include "AliMultSelection.h"
#include "AliAODMCParticle.h"
// Analysis task to fill histograms with PHOS AOD clusters and cells
// Authors: Yuri Kharlov
// Date   : 28.05.2009

ClassImp(AliAnalysisPhotonDDA)

//________________________________________________________________________
AliAnalysisPhotonDDA::AliAnalysisPhotonDDA(const char *name) 
: AliAnalysisTaskSE(name),
  fOutputContainer(nullptr),
  fCenBin(0),   
  fCentEstimator(1),
  fNCenBin(5),   
  fMinBCDistance(0.),
  fCentrality(0.),
  fRunNumber(-1),
  fIsMC(kTRUE),
  fPHOSGeo(nullptr),
  fEventCounter(0),
  fPIDResponse(nullptr),
  fPHOSEvent(nullptr),
  fCurrentMixedList(nullptr)
{
    
  // Output slots #0 write into a TH1 container
  DefineOutput(1,TList::Class());
  for(int i=0; i<5; i++)fPHOSEvents[i]=nullptr; 

}

//________________________________________________________________________
void AliAnalysisPhotonDDA::UserCreateOutputObjects()
{
  // Create histograms
  // Called once

  // AOD histograms
  if(fOutputContainer != NULL){
    delete fOutputContainer;
  }
  fOutputContainer = new THashList();
  fOutputContainer->SetOwner(kTRUE);

  fOutputContainer->Add(new TH1F("hSelEvents","Selected events",7,0.5,7.5));
  fOutputContainer->Add(new TH1F("hZvertex","Z vertex",200,-50.,+50.));
  fOutputContainer->Add(new TH1F("hCentrality","Centrality",105,0.,105.));
  fOutputContainer->Add(new TH1F("hNvertexTracks","N of primary tracks from the primary vertex",150,0.,150.));


  fOutputContainer->Add(new TH3F("hCPVr","CPV r",200,0.,20.,200,0.,20.,5,0.,100.));
  fOutputContainer->Add(new TH3F("hCPVrPi","CPV r",200,0.,20.,200,0.,20.,5,0.,100.));
  fOutputContainer->Add(new TH3F("hCPVrPrPlus","CPV r",200,0.,20.,200,0.,20.,5,0.,100.));
  fOutputContainer->Add(new TH3F("hCPVrPrMinus","CPV r",200,0.,20.,200,0.,20.,5,0.,100.));

  fOutputContainer->Add(new TH3F("hReCPVr","CPV r",200,0.,20.,200,0.,20.,5,0.,100.));
  fOutputContainer->Add(new TH3F("hReCPVrPi","CPV r",200,0.,20.,200,0.,20.,5,0.,100.));
  fOutputContainer->Add(new TH3F("hReCPVrPrPlus","CPV r",200,0.,20.,200,0.,20.,5,0.,100.));
  fOutputContainer->Add(new TH3F("hReCPVrPrMinus","CPV r",200,0.,20.,200,0.,20.,5,0.,100.));

  fOutputContainer->Add(new TH3F("hMiCPVr","CPV r",200,0.,20.,200,0.,20.,5,0.,100.));
  fOutputContainer->Add(new TH3F("hMiCPVrPi","CPV r",200,0.,20.,200,0.,20.,5,0.,100.));
  fOutputContainer->Add(new TH3F("hMiCPVrPrPlus","CPV r",200,0.,20.,200,0.,20.,5,0.,100.));
  fOutputContainer->Add(new TH3F("hMiCPVrPrMinus","CPV r",200,0.,20.,200,0.,20.,5,0.,100.));

  fOutputContainer->Add(new TH3F("hReDispCPVr","CPV r",200,0.,20.,200,0.,20.,5,0.,100.));
  fOutputContainer->Add(new TH3F("hReDispCPVrPi","CPV r",200,0.,20.,200,0.,20.,5,0.,100.));
  fOutputContainer->Add(new TH3F("hReDispCPVrPrPlus","CPV r",200,0.,20.,200,0.,20.,5,0.,100.));
  fOutputContainer->Add(new TH3F("hReDispCPVrPrMinus","CPV r",200,0.,20.,200,0.,20.,5,0.,100.));

  fOutputContainer->Add(new TH3F("hMiDispCPVr","CPV r",200,0.,20.,200,0.,20.,5,0.,100.));
  fOutputContainer->Add(new TH3F("hMiDispCPVrPi","CPV r",200,0.,20.,200,0.,20.,5,0.,100.));
  fOutputContainer->Add(new TH3F("hMiDispCPVrPrPlus","CPV r",200,0.,20.,200,0.,20.,5,0.,100.));
  fOutputContainer->Add(new TH3F("hMiDispCPVrPrMinus","CPV r",200,0.,20.,200,0.,20.,5,0.,100.));

  Int_t nPt      = 200;
  Double_t ptMin = 0;
  Double_t ptMax = 20;

  fOutputContainer->Add(new TH1F("hTrackPt" ,"p_{T} of all tracks",nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH1F("hPionPt"  ,"p_{T} of #pi^{#pm}",nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH1F("hKaonPt"  ,"p_{T} of K^{#pm}",nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH1F("hProtonPt","p_{T} of p,#bar{p}",nPt,ptMin,ptMax));
  fOutputContainer->Add(new TH1F("hUndefPt" ,"p_{T} of undefined tracks",nPt,ptMin,ptMax));

  fOutputContainer->Add(new TH1F("hTrackMult" ,"Charged track multiplicity",150,0.,150.));
  fOutputContainer->Add(new TH1F("hPionMult"  ,"#pi^{#pm} multiplicity"    ,150,0.,150.));
  fOutputContainer->Add(new TH1F("hKaonMult"  ,"K^{#pm} multiplicity"      ,150,0.,150.));
  fOutputContainer->Add(new TH1F("hProtonMult","p,#bar{p} multiplicity"    ,150,0.,150.));
  fOutputContainer->Add(new TH1F("hUndefMult" ,"Undefined multiplicity"    ,150,0.,150.));

  fOutputContainer->Add(new TH1I("hClusterMult"      ,"CaloCluster multiplicity"     ,100,0,100));

  fOutputContainer->Add(new TH2F("hEclPtPiPlus" ,"p_{track} vs E_{clu}, #pi^{+}",200,0.,20.,200,0.,20.));
  fOutputContainer->Add(new TH2F("hEclPtPiMinus","p_{track} vs E_{clu}, #pi^{-}",200,0.,20.,200,0.,20.));
  fOutputContainer->Add(new TH2F("hEclPtKPlus"  ,"p_{track} vs E_{clu}, K^{+}"  ,200,0.,20.,200,0.,20.));
  fOutputContainer->Add(new TH2F("hEclPtKMinus" ,"p_{track} vs E_{clu}, K^{-}"  ,200,0.,20.,200,0.,20.));
  fOutputContainer->Add(new TH2F("hEclPtPrPlus" ,"p_{track} vs E_{clu}, p"      ,200,0.,20.,200,0.,20.));
  fOutputContainer->Add(new TH2F("hEclPtPrMinus","p_{track} vs E_{clu}, #bar{p}",200,0.,20.,200,0.,20.));
  fOutputContainer->Add(new TH2F("hEclPtUndPlus" ,"p_{track} vs E_{clu}, p"      ,200,0.,20.,200,0.,20.));
  fOutputContainer->Add(new TH2F("hEclPtUndMinus","p_{track} vs E_{clu}, #bar{p}",200,0.,20.,200,0.,20.));

  TString parts[4]={"Pi","K","Pr","Und"} ;
  for(Int_t cen=0; cen<fNCenBin; cen++){

    fOutputContainer->Add(new TH3F(Form("hdxdzpt_cen%d",cen),"dx,dx,Pt",100,-50.,50.,100,-50.,50.,200,0.,20.));
    fOutputContainer->Add(new TH3F(Form("hdrEpt_cen%d",cen),"rEPt",100,0.,50.,200,0.,20.,200,0.,20.));
    for(Int_t ipid=0; ipid<4; ipid++){	  

     fOutputContainer->Add(new TH1F(Form("h%sPlusAll_cen%d",parts[ipid].Data(),cen),"Cluster spectrum",200,0.,20.));
     fOutputContainer->Add(new TH1F(Form("h%sMinusAll_cen%d",parts[ipid].Data(),cen),"Cluster spectrum",200,0.,20.));
     fOutputContainer->Add(new TH1F(Form("h%sPlusDisp_cen%d",parts[ipid].Data(),cen),"Cluster spectrum",200,0.,20.));
     fOutputContainer->Add(new TH1F(Form("h%sMinusDisp_cen%d",parts[ipid].Data(),cen),"Cluster spectrum",200,0.,20.));

     fOutputContainer->Add(new TH1F(Form("h%sPlusStrAll_cen%d",parts[ipid].Data(),cen),"Cluster spectrum",200,0.,20.));
     fOutputContainer->Add(new TH1F(Form("h%sMinusStrAll_cen%d",parts[ipid].Data(),cen),"Cluster spectrum",200,0.,20.));
     fOutputContainer->Add(new TH1F(Form("h%sPlusStrDisp_cen%d",parts[ipid].Data(),cen),"Cluster spectrum",200,0.,20.));
     fOutputContainer->Add(new TH1F(Form("h%sMinusStrDisp_cen%d",parts[ipid].Data(),cen),"Cluster spectrum",200,0.,20.));
     
     fOutputContainer->Add(new TH1F(Form("h%sPlusTrueAll_cen%d",parts[ipid].Data(),cen),"Cluster spectrum",200,0.,20.));
     fOutputContainer->Add(new TH1F(Form("h%sMinusTrueAll_cen%d",parts[ipid].Data(),cen),"Cluster spectrum",200,0.,20.));
     fOutputContainer->Add(new TH1F(Form("h%sPlusTrueDisp_cen%d",parts[ipid].Data(),cen),"Cluster spectrum",200,0.,20.));
     fOutputContainer->Add(new TH1F(Form("h%sMinusTrueDisp_cen%d",parts[ipid].Data(),cen),"Cluster spectrum",200,0.,20.));     
    } 
    
    fOutputContainer->Add(new TH1F(Form("hPhotAll_cen%d",cen),"Cluster spectrum",200,0.,20.));
    fOutputContainer->Add(new TH1F(Form("hPhotCPV_cen%d",cen),"Cluster spectrum",200,0.,20.));
    fOutputContainer->Add(new TH1F(Form("hPhotDisp_cen%d",cen),"Cluster spectrum",200,0.,20.));
    fOutputContainer->Add(new TH1F(Form("hPhotBoth_cen%d",cen),"Cluster spectrum",200,0.,20.));
  }
 
  fOutputContainer->Add(new TH2F("hClusterTOF","Cluster spectrum",200,0.,20.,400,-200.e-9,200.e-9));
  fOutputContainer->Add(new TH2F("hClusterTOFn2","Cluster spectrum",200,0.,20.,400,-200.e-9,200.e-9));
  fOutputContainer->Add(new TH2F("hClusterTOFm02","Cluster spectrum",200,0.,20.,400,-200.e-9,200.e-9));
  
 //MC  
  char partName[15][10] ;
  snprintf(partName[0],10,"pi0");
  snprintf(partName[1],10,"eta") ;
  snprintf(partName[2],10,"omega"); 
  snprintf(partName[3],10,"K0s"); 
  snprintf(partName[4],10,"Kpm"); 
  snprintf(partName[5],10,"pipm"); 
  snprintf(partName[6],10,"n"); 
  snprintf(partName[7],10,"nbar"); 
  snprintf(partName[8],10,"p"); 
  snprintf(partName[9],10,"pbar"); 
  snprintf(partName[10],10,"el"); 
  snprintf(partName[11],10,"OtherCh"); 
  snprintf(partName[12],10,"OtherNeu"); 
  
  char cPID[8][15] ;
  snprintf(cPID[0],5,"All") ;
  snprintf(cPID[1],5,"Disp");
  snprintf(cPID[2],5,"CPV") ;
  snprintf(cPID[3],5,"Both"); 
  
  if(fIsMC){
    for(Int_t ipart=0; ipart<13; ipart++){  
      for(Int_t iPID=0; iPID<4; iPID++){
        fhCont2D[ipart][iPID]= new TH2F(Form("hMCRec_%s_%s",partName[ipart],cPID[iPID]),"Rec vs primary",200,0.,20.,200,0.,20.) ;    
        fOutputContainer->Add(fhCont2D[ipart][iPID]) ;
      }
    }
  }
   
  for(Int_t j=0;j<fNCenBin;j++)
    fPHOSEvents[j]=0x0 ;    //Container for PHOS photons

  PostData(1, fOutputContainer);

}

//________________________________________________________________________
void AliAnalysisPhotonDDA::UserExec(Option_t *) 
{
  // Main loop, called for each event
  // Analyze AOD
  // const Double_t kEcrossCut=0.98 ;
  const Double_t kTOFMaxCut= 30.e-9 ;  
  const Double_t kTOFMinCut=-30.e-9 ;  


  // Initialize the PHOS geometry
  if(!fPHOSGeo)
    fPHOSGeo = AliPHOSGeometry::GetInstance() ;

  // Event selection flags
  FillHistogram("hSelEvents",1.) ;

  AliAODEvent *event = dynamic_cast<AliAODEvent*>(InputEvent());
  if (!event) {
     Printf("ERROR: Could not retrieve event");
     return;
  }
  
  TClonesArray* stack = (TClonesArray*)event->FindListObject(AliAODMCParticle::StdBranchName());

  // Checks if we have a primary vertex
  // Get primary vertices form AOD

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  fPIDResponse = inputHandler->GetPIDResponse();

  Double_t vtx5[3];
  vtx5[0] = event->GetPrimaryVertex()->GetX();
  vtx5[1] = event->GetPrimaryVertex()->GetY();
  vtx5[2] = event->GetPrimaryVertex()->GetZ();

  FillHistogram("hNvertexTracks",event->GetPrimaryVertex()->GetNContributors());
  FillHistogram("hZvertex"      ,event->GetPrimaryVertex()->GetZ());
  if (TMath::Abs(event->GetPrimaryVertex()->GetZ()) > 10. )
    return ;

  fRunNumber=event->GetRunNumber() ;

  FillHistogram("hSelEvents",2.) ;

  Int_t runNumber=event->GetRunNumber() ;
  fCentrality=101.;
  if(runNumber >=195344 && runNumber <= 197388 ){ //LHC13bcdef
    AliCentrality *centrality = event->GetCentrality();
    if( centrality ){
    switch(fCentEstimator){
      case 4 : fCentrality=centrality->GetCentralityPercentile("CL1");
        break;
      case 3 : fCentrality=centrality->GetCentralityPercentile("ZNA");
        break;         
      case 2 : fCentrality=centrality->GetCentralityPercentile("V0M");
                     break;
      case 1 : 
      default: fCentrality=centrality->GetCentralityPercentile("V0A");
    }
    }
  }
  else{ //Run2
    AliMultSelection *multSelection = (AliMultSelection*) event -> FindListObject("MultSelection");
    if(multSelection){
    switch(fCentEstimator){
      case 4 : fCentrality=multSelection->GetMultiplicityPercentile("CL1");
        break;
      case 3 : fCentrality=multSelection->GetMultiplicityPercentile("ZNA");
        break;         
      case 2 : fCentrality=multSelection->GetMultiplicityPercentile("V0M");
        break;
      case 1 : 
      default: fCentrality=multSelection->GetMultiplicityPercentile("V0A");
    }
    }

  }           

  FillHistogram("hCentrality",fCentrality) ;
  
  if( fCentrality <= 0. || fCentrality>100. ){
    PostData(1, fOutputContainer);
    return;
  }

  FillHistogram("hSelEvents",3.) ;


  fCenBin=0;
  while(fCenBin<fNCenBin && fCentrality>fCenBinEdges.At(fCenBin))
    fCenBin++ ;
  if(fCenBin>=fNCenBin) fCenBin=fNCenBin-1; 


  if(!fPHOSEvents[fCenBin]) 
    fPHOSEvents[fCenBin]=new TList() ;
  fCurrentMixedList = fPHOSEvents[fCenBin] ;

  
  
  //Calculate charged multiplicity
  Int_t trackMult=0, nPion=0, nKaon=0, nProton=0, nUndef=0;

  for (Int_t i=0;i<event->GetNumberOfTracks();++i) {
    AliAODTrack *track = (AliAODTrack*)event->GetTrack(i) ;
    if(!track->IsHybridGlobalConstrainedGlobal())
      continue ;
    if(  TMath::Abs(track->Eta())< 0.8) {
        
      trackMult++;
      Double_t pt = track->Pt();
      FillHistogram("hTrackPt",pt);
      if(fPIDResponse) {
	      Bool_t pidPion=kFALSE, pidKaon=kFALSE, pidProton=kFALSE, pidUndef=kFALSE;
	      Double_t nsigmaProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton));
	      Double_t nsigmaKaon   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon)); 
	      Double_t nsigmaPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion)); 

	      // guess the particle based on the smaller nsigma
      	if((nsigmaKaon  <nsigmaPion) && (nsigmaKaon  <nsigmaProton) && (nsigmaKaon  <3)) pidKaon   = kTRUE;
	      if((nsigmaPion  <nsigmaKaon) && (nsigmaPion  <nsigmaProton) && (nsigmaPion  <3)) pidPion   = kTRUE;
	      if((nsigmaProton<nsigmaKaon) && (nsigmaProton<nsigmaPion  ) && (nsigmaProton<3)) pidProton = kTRUE;
	      if (!pidPion && !pidKaon && !pidProton) pidUndef = kTRUE;

	      if (pidPion) {
	        nPion++;
	        FillHistogram("hPionPt",pt);
	      }
	      if (pidKaon) {
	        nKaon++;
	        FillHistogram("hKaonPt",pt);
	      }
	      if (pidProton) {
	        nProton++;
	        FillHistogram("hProtonPt",pt);
	      }
	      if (pidUndef) {
	        nUndef++;
	        FillHistogram("hUndefPt",pt);
	      }
      }
    }
  }
  FillHistogram("hTrackMult" ,trackMult+0.5) ;
  FillHistogram("hPionMult"  ,nPion    +0.5);
  FillHistogram("hKaonMult"  ,nKaon    +0.5);
  FillHistogram("hProtonMult",nProton  +0.5);
  FillHistogram("hUndefMult" ,nUndef   +0.5);

  AliAODCaloCluster *clu;
  TLorentzVector p1,p2,p12, pv1,pv2;

  Int_t multClust = event->GetNumberOfCaloClusters();
  FillHistogram("hClusterMult",multClust);

  
  if(!fPHOSEvent)
    fPHOSEvent   = new TClonesArray("AliCaloPhoton",multClust);
  else
    fPHOSEvent->Clear() ;
  

  //For re-calibration
  TVector3 vertex(vtx5);
  TVector3 localPos ;
  Int_t inList=0;
  for (Int_t i1=0; i1<multClust; i1++) {
    clu = event->GetCaloCluster(i1);
    if ( !clu->IsPHOS() ) continue;
    if ( clu->E()<0.1 ) continue;

    if(clu->GetDistanceToBadChannel()<fMinBCDistance)
      continue ;    
    
    Double_t energy=clu->E() ;
     
    FillHistogram("hClusterTOF",energy,clu->GetTOF());
    Bool_t exotic=kFALSE ;
    if(clu->GetNCells()<3){
      FillHistogram("hClusterTOFn2",energy,clu->GetTOF());
      if(clu->E()>1.)
        exotic=kTRUE;
    }
    if(clu->GetM02()<0.2){
      FillHistogram("hClusterTOFm02",energy,clu->GetTOF());
      if(clu->E()>1.)
        exotic=kTRUE;  
    }
    if(exotic) continue ;
  
//     if(clu->GetMCEnergyFraction()>kEcrossCut) //Ecross cut, should be filled with Tender
//      continue ;    
 
    if(clu->GetTOF() < kTOFMinCut || clu->GetTOF() > kTOFMaxCut)
      continue ;          

    clu->GetMomentum(pv1 ,vtx5);
     
    AliCaloPhoton *p = new ((*fPHOSEvent)[inList]) AliCaloPhoton(pv1.Px(),pv1.Py(),pv1.Pz(),clu->E() );
    inList++;
    Float_t pos[3] ;
    clu->GetPosition(pos) ;
    TVector3 global1(pos) ;
    Int_t relId[4] ;
    fPHOSGeo->GlobalPos2RelId(global1,relId) ;
    Int_t mod  = relId[0] ;    
    TVector3 local ;
    fPHOSGeo->Global2Local(local,global1,mod);
    p->SetModule(mod) ;
    p->SetEMCx(local.X()) ;
    p->SetEMCz(local.Z()) ;
   
    //Dispersion bit  
    Bool_t dispBit=clu->GetDispersion()<2.5*2.5; //clu->Chi2()<2.5*2.5 ;
    p->SetDispBit(dispBit) ;
    
    //Track matching
    Double_t dx=clu->GetTrackDx() ;
    Double_t dz=clu->GetTrackDz() ;

    
    Bool_t cpvBit=kTRUE ; //No track matched by default
    Bool_t cpvBit2=kTRUE ; //More Strict criterion
    Double_t r=100.;
    AliAODTrack* track = 0x0 ;
    if(clu->GetNTracksMatched()>0){
      track = dynamic_cast<AliAODTrack*> (clu->GetTrackMatched(0));
      r=clu->GetEmcCpvDistance();
    }
    cpvBit=(r>2.) ;
    cpvBit2=(r>1.) ;
    FillHistogram("hCPVr",r,pv1.E(),fCentrality);

    
    
    FillHistogram(Form("hPhotAll_cen%d",fCenBin),pv1.Pt()) ;
    if(cpvBit){
      FillHistogram(Form("hPhotCPV_cen%d",fCenBin),pv1.Pt()) ;
    }
    if(dispBit){
      FillHistogram(Form("hPhotDisp_cen%d",fCenBin),pv1.Pt()) ;
      if(cpvBit){
        FillHistogram(Form("hPhotBoth_cen%d",fCenBin),pv1.Pt()) ;
      }
    }
        
    //classify parents    
    if(fIsMC){
       Int_t primLabel=clu->GetLabelAt(0) ; 
       //Look what particle left vertex
       if(primLabel>-1){
         AliAODMCParticle * prim = (AliAODMCParticle*)stack->At(primLabel) ;
         Int_t iparent=primLabel;
         AliAODMCParticle * parent = prim;
         Double_t r2=prim->Xv()*prim->Xv()+prim->Yv()*prim->Yv() ;
         while((r2 > 1.) && (iparent>-1)){
           iparent=parent->GetMother();
           if(iparent<0)
	           break ;
           parent=(AliAODMCParticle*)stack->At(iparent);
           r2=parent->Xv()*parent->Xv()+parent->Yv()*parent->Yv() ;
         }
         Int_t parentPDG=parent->GetPdgCode() ;    
         Int_t parentIndx=0;
         switch(parentPDG){
	       case 22: iparent=parent->GetMother();
                 if(iparent>=0){
                    parent=(AliAODMCParticle*)stack->At(iparent);
                    parentPDG=parent->GetPdgCode() ;
                    if(parentPDG==111){ //pi0
                      parentIndx=0;   
                      break ;  
                    }
                    if(parentPDG==221){ //eta
                      parentIndx=1;   
                      break ;  
                    }
                    if(parentPDG==223){ //omega
                      parentIndx=2;   
                      break ;  
                    }
                 }
                 parentIndx=12;
                 break ;

         case 111: parentIndx=0;   
                  break ; 
         case 221: parentIndx=1;   
                  break ;  
	       case  11:
	       case -11: parentIndx=10; //e+-  
                  break ;  
	 case -2212: parentIndx=9; //pbar 
                  break ;
	 case -2112: parentIndx=7; //nbar 
                  break ;	  
	 case  211:
	 case -211: parentIndx=5; //pbar 
                  break ;
	 case 2212:parentIndx=8; //proton 
                  break ;	  
	 case  321:
	 case -321:parentIndx=4; //K+- 
                  break ;
	 case 310: parentIndx=3; //Ks0
	          break ;
	 case 2112: parentIndx=6; //n 
                  break ;	
	 default:  
	   if(parent->Charge()!=0)
	     parentIndx=11; //other charged
	   else 
	     parentIndx=12; //other neutral
  }  
        
         fhCont2D[parentIndx][0]->Fill(clu->E(),parent->Pt()) ;
         if(cpvBit){
           fhCont2D[parentIndx][2]->Fill(clu->E(),parent->Pt()) ;
         }
         if(dispBit){
           fhCont2D[parentIndx][1]->Fill(clu->E(),parent->Pt()) ;
           if(cpvBit){
             fhCont2D[parentIndx][3]->Fill(clu->E(),parent->Pt()) ;
           }
         }
       }
    } 
        
    if (!track) continue;
    if(!track->IsHybridGlobalConstrainedGlobal())
      continue ;
    
    Double_t ptTrack = track->Pt();
    Short_t charge   = track->Charge();

    FillHistogram(Form("hdxdzpt_cen%d",fCenBin),dx,dz,ptTrack);
    FillHistogram(Form("hdrEpt_cen%d",fCenBin),r,energy,ptTrack);
      
    if(fPIDResponse) {
	    Double_t nsigmaProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton));
	    Double_t nsigmaKaon   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon)); 
	    Double_t nsigmaPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion)); 

    	Bool_t isPart[4] ; // pion, Kaon, proton, other
	    isPart[0]=(nsigmaPion  <nsigmaKaon) && (nsigmaPion  <nsigmaProton) && (nsigmaPion  <1); // = pidPion   = kTRUE;
	    isPart[1]=(nsigmaKaon  <nsigmaPion) && (nsigmaKaon  <nsigmaProton) && (nsigmaKaon  <1); // = pidKaon   = kTRUE;
	    isPart[2]=(nsigmaProton<nsigmaKaon) && (nsigmaProton<nsigmaPion  ) && (nsigmaProton<1); // = pidProton = kTRUE;
	    isPart[3]= !isPart[0] && !isPart[1] && !isPart[2]; // pidUndef = kTRUE;
	
        
      if(isPart[0]) FillHistogram("hCPVrPi",r,pv1.E(),fCentrality);
      if(isPart[2]){
        if(charge>0) 
          FillHistogram("hCPVrPrPlus",r,pv1.E(),fCentrality);
        else
          FillHistogram("hCPVrPrMinus",r,pv1.E(),fCentrality);
      }        
        
      if(cpvBit) //only charged clusters
        continue ; 
   
      TString parts[4]={"Pi","K","Pr","Und"} ;
	    for(Int_t ipid=0; ipid<4; ipid++){
	      if(!isPart[ipid]) continue ;
	  	  
	      char name[55] ;
	      if(charge>0)snprintf(name,55,"hEclPt%sPlus",parts[ipid].Data()) ;
	      else snprintf(name,55,"hEclPt%sMinus",parts[ipid].Data()) ;
	  
 	      FillHistogram(name,energy,ptTrack);
	
	      if(charge>0)snprintf(name,55,"%sPlus",parts[ipid].Data()) ;
	      else snprintf(name,55,"%sMinus",parts[ipid].Data()) ;

	      FillHistogram(Form("h%sAll_cen%d",name,fCenBin),pv1.Pt()) ;
        if(dispBit){
          FillHistogram(Form("h%sDisp_cen%d",name,fCenBin),pv1.Pt()) ;
    	  }
	  
	      if(!cpvBit2){ //more strict match
	        if(charge>0)snprintf(name,55,"%sPlusStr",parts[ipid].Data()) ;
	        else snprintf(name,55,"%sMinusStr",parts[ipid].Data()) ;

  	      FillHistogram(Form("h%sAll_cen%d",name,fCenBin),pv1.Pt()) ;
          if(dispBit){
            FillHistogram(Form("h%sDisp_cen%d",name,fCenBin),pv1.Pt()) ;
 	        }
	      }
	 
	      if(energy<1.2*ptTrack || (ipid==2 && charge<0 && energy<1.8+1.2*ptTrack)){ 
	        if(charge>0)snprintf(name,55,"%sPlusTrue",parts[ipid].Data()) ;
	        else snprintf(name,55,"%sMinusTrue",parts[ipid].Data()) ;

  	      FillHistogram(Form("h%sAll_cen%d",name,fCenBin),pv1.Pt()) ;
          if(dispBit){
            FillHistogram(Form("h%sDisp_cen%d",name,fCenBin),pv1.Pt()) ;
	        }
	      }
	    }
    }
  }
  
  //Check of track matching
  //Real
  for(Int_t i=0; i<fPHOSEvent->GetEntriesFast(); i++){
    
    AliCaloPhoton * phot = static_cast<AliCaloPhoton*>(fPHOSEvent->At(i)) ; 
    Int_t mod = phot->Module() ;
    TVector3 locpos(phot->EMCx(),0.,phot->EMCz()) ;
    
    Double_t dx=999.,dz=999., pttrack=0.;
    Int_t charge=0;
    Int_t itr = FindTrackMatching(mod, &locpos, dx,dz, pttrack, charge);
    if(itr<0)
      continue ;  
    Double_t r=TestCPV(dx, dz, pttrack,charge) ;
    FillHistogram("hReCPVr",r,pv1.E(),fCentrality);
    if(phot->IsDispOK()){
      FillHistogram("hReDispCPVr",r,pv1.E(),fCentrality);
    }
    
    AliAODTrack *track= (AliAODTrack*)event->GetTrack(itr);      
    if (!track) continue;
    if(!track->IsHybridGlobalConstrainedGlobal())
      continue ;
          
    Double_t nsigmaProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton));
    Double_t nsigmaKaon   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon)); 
    Double_t nsigmaPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion)); 

    Bool_t isPart[4] ; // pion, Kaon, proton, other
    isPart[0]=(nsigmaPion  <nsigmaKaon) && (nsigmaPion  <nsigmaProton) && (nsigmaPion  <3); // = pidPion   = kTRUE;
// 	isPart[1]=(nsigmaKaon  <nsigmaPion) && (nsigmaKaon  <nsigmaProton) && (nsigmaKaon  <3); // = pidKaon   = kTRUE;
    isPart[2]=(nsigmaProton<nsigmaKaon) && (nsigmaProton<nsigmaPion  ) && (nsigmaProton<3); // = pidProton = kTRUE;
// 	isPart[3]= !isPart[0] && !isPart[1] && !isPart[2]; // pidUndef = kTRUE;
	
        
    if(isPart[0]){
        FillHistogram("hReCPVrPi",r,pv1.E(),fCentrality);
        if(phot->IsDispOK()){
          FillHistogram("hReDispCPVrPi",r,pv1.E(),fCentrality);
        }
    }
    if(isPart[2]){
      if(charge>0) 
        FillHistogram("hReCPVrPrPlus",r,pv1.E(),fCentrality);
      else
        FillHistogram("hReCPVrPrMinus",r,pv1.E(),fCentrality);
    }        
    if(phot->IsDispOK()){
      if(charge>0) 
        FillHistogram("hReDispCPVrPrPlus",r,pv1.E(),fCentrality);
      else
        FillHistogram("hReDispCPVrPrMinus",r,pv1.E(),fCentrality);
    }
  }
  
  //Mixed
  //Fill Mixed InvMass distributions:
  TIter nextEv(fCurrentMixedList) ;
  while(TClonesArray * event2 = static_cast<TClonesArray*>(nextEv())){
   Int_t nPhotons2 = event2->GetEntriesFast() ;
   for(Int_t j=0; j < nPhotons2 ; j++){
    AliCaloPhoton * phot = static_cast<AliCaloPhoton*>(event2->At(j)) ;
    Int_t mod = phot->Module() ;
    TVector3 locpos(phot->EMCx(),0.,phot->EMCz()) ;
    
    Double_t dx=999.,dz=999., pttrack=0.;
    Int_t charge=0;
    Int_t itr = FindTrackMatching(mod, &locpos, dx,dz, pttrack, charge);
     if(itr<0)
      continue ;  
    Double_t r=TestCPV(dx, dz, pttrack,charge) ;
    FillHistogram("hMiCPVr",r,pv1.E(),fCentrality);
    if(phot->IsDispOK()){
      FillHistogram("hMiDispCPVr",r,pv1.E(),fCentrality);
    }
    
    AliAODTrack *track= (AliAODTrack*)event->GetTrack(itr);      
    if (!track) continue;
    if(!track->IsHybridGlobalConstrainedGlobal())
      continue ;
          
    Double_t nsigmaProton = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton));
    Double_t nsigmaKaon   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kKaon)); 
    Double_t nsigmaPion   = TMath::Abs(fPIDResponse->NumberOfSigmasTPC(track, AliPID::kPion)); 

    Bool_t isPart[4] ; // pion, Kaon, proton, other
    isPart[0]=(nsigmaPion  <nsigmaKaon) && (nsigmaPion  <nsigmaProton) && (nsigmaPion  <3); // = pidPion   = kTRUE;
// 	isPart[1]=(nsigmaKaon  <nsigmaPion) && (nsigmaKaon  <nsigmaProton) && (nsigmaKaon  <3); // = pidKaon   = kTRUE;
    isPart[2]=(nsigmaProton<nsigmaKaon) && (nsigmaProton<nsigmaPion  ) && (nsigmaProton<3); // = pidProton = kTRUE;
// 	isPart[3]= !isPart[0] && !isPart[1] && !isPart[2]; // pidUndef = kTRUE;
	
        
    if(isPart[0]){
        FillHistogram("hMiCPVrPi",r,pv1.E(),fCentrality);
        if(phot->IsDispOK()){
          FillHistogram("hMiDispCPVrPi",r,pv1.E(),fCentrality);
        }
    }
    if(isPart[2]){
      if(charge>0) 
        FillHistogram("hMiCPVrPrPlus",r,pv1.E(),fCentrality);
      else
        FillHistogram("hMiCPVrPrMinus",r,pv1.E(),fCentrality);
    }        
    if(phot->IsDispOK()){
      if(charge>0) 
        FillHistogram("hMiDispCPVrPrPlus",r,pv1.E(),fCentrality);
      else
        FillHistogram("hMiDispCPVrPrMinus",r,pv1.E(),fCentrality);
    }
   }
  }
  
  
  //Remove old events
  fCurrentMixedList->AddFirst(fPHOSEvent);
  fPHOSEvent=0x0 ;
  if(fCurrentMixedList->GetSize() > 20){
    TClonesArray *tmp = static_cast <TClonesArray*> (fCurrentMixedList->Last());
    fCurrentMixedList->Remove(tmp);
    delete tmp;
  }
    
  
  // Post output data.
  PostData(1, fOutputContainer);
  fEventCounter++;
}

//________________________________________________________________________
void AliAnalysisPhotonDDA::Terminate(Option_t *)
{
}

//_____________________________________________________________________________
void AliAnalysisPhotonDDA::FillHistogram(const char * key,Double_t x)const{
  //FillHistogram
  TH1 * hist = dynamic_cast<TH1*>(fOutputContainer->FindObject(key)) ;
  if(hist)
    hist->Fill(x) ;
  else
    AliError(Form("can not find histogram (of instance TH1) <%s> ",key)) ;
}
//_____________________________________________________________________________
void AliAnalysisPhotonDDA::FillHistogram(const char * key,Double_t x,Double_t y)const{
  //FillHistogram
  TH1 * th1 = dynamic_cast<TH1*> (fOutputContainer->FindObject(key));
  if(th1)
    th1->Fill(x, y) ;
  else
    AliError(Form("can not find histogram (of instance TH1) <%s> ",key)) ;
}

//_____________________________________________________________________________
void AliAnalysisPhotonDDA::FillHistogram(const char * key,Double_t x,Double_t y, Double_t z) const{
  //Fills 1D histograms with key
  TObject * obj = fOutputContainer->FindObject(key);
  
  TH2 * th2 = dynamic_cast<TH2*> (obj);
  if(th2) {
    th2->Fill(x, y, z) ;
    return;
  }

  TH3 * th3 = dynamic_cast<TH3*> (obj);
  if(th3) {
    th3->Fill(x, y, z) ;
    return;
  }
  
  AliError(Form("can not find histogram (of instance TH2) <%s> ",key)) ;
}
//___________________________________________________________________________________________________
Int_t AliAnalysisPhotonDDA::FindTrackMatching(Int_t mod,TVector3 *locpos,
					    Double_t &dx, Double_t &dz,
					    Double_t &pt,Int_t &charge){
  //Find track with closest extrapolation to cluster
  AliAODEvent *aod = static_cast<AliAODEvent*>(InputEvent());

  
  Double_t  magF = aod->GetMagneticField();
 
  Double_t magSign = 1.0;
  if(magF<0)magSign = -1.0;
  
  if (!TGeoGlobalMagField::Instance()->GetField()) {
    AliError("Margnetic filed was not initialized, use default") ;
    AliMagF* field = new AliMagF("Maps","Maps", magSign, magSign, AliMagF::k5kG);
    TGeoGlobalMagField::Instance()->SetField(field);
  }

  // *** Start the matching
  Int_t nt = aod->GetNumberOfTracks();
      
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
        aodTrack=(AliAODTrack*)aod->GetTrack(i);      
      
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
	  dx = ddx ;
  	  dz = ddz ;
	  minDistance=d2 ;
	  itr=i ;
          pt=aodTrack->Pt() ;
          charge=aodTrack->Charge() ;
        }
      }
    }//Scanned all tracks
 
   return itr ;
}
//____________________________________________________________________________
Double_t AliAnalysisPhotonDDA::TestCPV(Double_t dx, Double_t dz, Double_t pt, Int_t charge){
  //Parameterization of LHC10h period
  //_true if neutral_
  
  AliAODEvent *aod= dynamic_cast<AliAODEvent*>(InputEvent());
  
  Double_t mf = aod->GetMagneticField();
  
   Double_t meanX=0;
   Double_t meanZ=0.;
   Double_t sx=0.; 
   Double_t sz=0.; 
  if(fRunNumber<209122){ //Run1
    sx=TMath::Min(5.4,2.59719e+02*TMath::Exp(-pt/1.02053e-01)+
              6.58365e-01*5.91917e-01*5.91917e-01/((pt-9.61306e-01)*(pt-9.61306e-01)+5.91917e-01*5.91917e-01)+1.59219);
    sz=TMath::Min(2.75,4.90341e+02*1.91456e-02*1.91456e-02/(pt*pt+1.91456e-02*1.91456e-02)+1.60) ;
  
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

  }
  else{//Run2
  
    sx = TMath::Min(5.2, 1.111 + 0.56 * TMath::Exp(-0.031 * pt*pt) + 4.8 /TMath::Power(pt+0.61,3));
    sz = TMath::Min(3.3, 1.12  + 0.35 * TMath::Exp(-0.032 * pt*pt) + 0.75/TMath::Power(pt+0.24,3));

    if(mf<0.){ //field --
      meanZ = 0.102;
      if(charge>0)
        meanX =  TMath::Min(5.8, 0.42 + 0.70 * TMath::Exp(-0.015 * pt*pt) + 35.8/TMath::Power(pt+1.41,3));
      else
        meanX = -TMath::Min(5.8, 0.17 + 0.64 * TMath::Exp(-0.019 * pt*pt) + 26.1/TMath::Power(pt+1.21,3));
    }
    else{ //Field ++
      meanZ= 0.102;
      if(charge>0)
        meanX = -TMath::Min(5.8, 0.58 + 0.68 * TMath::Exp(-0.027 * pt*pt) + 28.0/TMath::Power(pt+1.28,3));
      else
        meanX =  TMath::Min(5.8, 0.11 + 0.67 * TMath::Exp(-0.015 * pt*pt) + 29.9/TMath::Power(pt+1.29,3));
    }

  }
  Double_t rz=(dz-meanZ)/sz ;
  Double_t rx=(dx-meanX)/sx ;
  return TMath::Sqrt(rx*rx+rz*rz) ;
  
}

    
