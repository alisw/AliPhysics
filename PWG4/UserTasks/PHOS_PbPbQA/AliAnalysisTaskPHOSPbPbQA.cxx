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
/* $Id$ */

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskPHOSPbPbQA.h"
#include "AliCaloPhoton.h"
#include "AliPHOSGeometry.h"
#include "AliESDEvent.h"
#include "AliESDCaloCells.h"
#include "AliCentrality.h"
#include "AliLog.h"
#include "TObjArray.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"

// Stripped-down version of Dmitri Peressounko' AliAnalysisTaskPi0Flow class
// used for the fast QA of PbPb data.
//...
// Author: Boris Polishchuk 
// Date   : 19.10.2011

ClassImp(AliAnalysisTaskPHOSPbPbQA)

//________________________________________________________________________
AliAnalysisTaskPHOSPbPbQA::AliAnalysisTaskPHOSPbPbQA() : AliAnalysisTaskSE(),
  fOutputContainer(0),fPHOSEvent(0),fCentrality(0),fCenBin(0),
  fPHOSGeo(0),fEventCounter(0)
{
  //Default constructor
  
  for(Int_t i=0;i<1;i++){
    for(Int_t j=0;j<2;j++)
      fPHOSEvents[i][j]=0 ;
  }
  
}

//________________________________________________________________________
AliAnalysisTaskPHOSPbPbQA::AliAnalysisTaskPHOSPbPbQA(const char *name) 
: AliAnalysisTaskSE(name),
  fOutputContainer(0),
  fPHOSEvent(0),
  fCentrality(0),fCenBin(0),
  fPHOSGeo(0),
  fEventCounter(0)
{
  // Constructor
  for(Int_t i=0;i<1;i++){
    for(Int_t j=0;j<2;j++)
	fPHOSEvents[i][j]=0 ;
  }
  
  // Output slots #0 write into a TH1 container
  DefineOutput(1,TList::Class());

  // Initialize the PHOS geometry
  fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP") ;

}

//________________________________________________________________________
void AliAnalysisTaskPHOSPbPbQA::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  
  // ESD histograms
  if(fOutputContainer != NULL){
    delete fOutputContainer;
  }

  fOutputContainer = new TList();
  fOutputContainer->SetOwner(kTRUE);
  

  fOutputContainer->Add(new TH2F("hCenPHOS","Centrality vs PHOSclusters", 100,0.,100.,200,0.,200.)) ;
  fOutputContainer->Add(new TH2F("hCenPHOSCells","Centrality vs PHOS cells", 100,0.,100.,100,0.,1000.)) ;
  fOutputContainer->Add(new TH2F("hCenTrack","Centrality vs tracks", 100,0.,100.,100,0.,15000.)) ;  
  
  //pi0 spectrum
  Int_t nPtPhot = 300 ;
  Double_t ptPhotMax = 30 ;
  Int_t nM       = 500;
  Double_t mMin  = 0.0;
  Double_t mMax  = 1.0;
  char key[55] ;

  for(Int_t cent=0; cent<2; cent++){
    
    snprintf(key,55,"hPi0All_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    
    snprintf(key,55,"hMiPi0All_cen%d",cent) ;
    fOutputContainer->Add(new TH2F(key,"All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
  }

  //per module
  for(Int_t cent=0; cent<2; cent++){
    for(Int_t sm=1; sm<4; sm++) {
      
      snprintf(key,55,"hPi0AllSM%d_cen%d",sm,cent) ;
      fOutputContainer->Add(new TH2F(key,"All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
      
      snprintf(key,55,"hMiPi0AllSM%d_cen%d",sm,cent) ;
      fOutputContainer->Add(new TH2F(key,"All clusters",nM,mMin,mMax,nPtPhot,0.,ptPhotMax));
    }
  }
    
  PostData(1, fOutputContainer);
  
}

//________________________________________________________________________
void AliAnalysisTaskPHOSPbPbQA::UserExec(Option_t *) 
{
  // Main loop, called for each event
  // Analyze ESD/AOD  
  
  AliESDEvent *event = dynamic_cast<AliESDEvent*>(InputEvent());

  if (!event) {
    Printf("ERROR: Could not retrieve event");
    PostData(1, fOutputContainer);
    return;
  }
  
  // Get PHOS rotation matrices from ESD and set them to the PHOS geometry
  char key[55] ;  
  
  if(fEventCounter == 0) {
    for(Int_t mod=0; mod<5; mod++) {
      if(!event->GetPHOSMatrix(mod)) continue;
      fPHOSGeo->SetMisalMatrix(event->GetPHOSMatrix(mod),mod) ;
    }
    fEventCounter++ ;
  }
   
  Int_t zvtx=0 ;

  AliCentrality *centrality = event->GetCentrality(); 
  fCentrality=centrality->GetCentralityPercentile("V0M");

  if( fCentrality <= 0. ){
    PostData(1, fOutputContainer);
    return;
  }

  if(fCentrality < 20.) fCenBin = 0;
  else fCenBin = 1;
  
  printf("centrality %.3f [%d]\n",fCentrality,fCenBin);
  
  if(!fPHOSEvents[zvtx][fCenBin]) 
    fPHOSEvents[zvtx][fCenBin]=new TList() ;
 
  TList * prevPHOS = fPHOSEvents[zvtx][fCenBin] ;
  
  if(fPHOSEvent)
    fPHOSEvent->Clear() ;
  else
    fPHOSEvent = new TClonesArray("AliCaloPhoton",200) ;
  
  Int_t multClust = event->GetNumberOfCaloClusters();
  AliESDCaloCells * cells = event->GetPHOSCells() ;

  FillHistogram("hCenPHOSCells",fCentrality,cells->GetNumberOfCells()) ;
  FillHistogram("hCenTrack",fCentrality,event->GetNumberOfTracks()) ;
  
  Int_t inPHOS = 0;
  Double_t vtx0[3] = {0,0,0}; // vertex

  for (Int_t i=0; i<multClust; i++) {

    AliESDCaloCluster *clu = event->GetCaloCluster(i);

    if ( !clu->IsPHOS() || clu->E()<0.3) continue;
    if(clu->GetNCells()<3) continue;
    
    Float_t  position[3];
    clu->GetPosition(position);
    TVector3 global(position) ;
    Int_t relId[4] ;
    fPHOSGeo->GlobalPos2RelId(global,relId) ;
    Int_t mod  = relId[0] ;
    Int_t cellX = relId[2];
    Int_t cellZ = relId[3] ;   
    
    TLorentzVector pv1 ;
    clu->GetMomentum(pv1 ,vtx0);
    
    if(inPHOS>=fPHOSEvent->GetSize()){
      fPHOSEvent->Expand(inPHOS+50) ;
    }

    new((*fPHOSEvent)[inPHOS]) AliCaloPhoton(pv1.X(),pv1.Py(),pv1.Z(),pv1.E()) ;
    AliCaloPhoton * ph = (AliCaloPhoton*)fPHOSEvent->At(inPHOS) ;
    ph->SetModule(mod) ;
    ph->SetMomV2(&pv1) ;
    ph->SetNCells(clu->GetNCells());
  
    ph->SetEMCx(float(cellX)) ;
    ph->SetEMCz(float(cellZ)) ;
    
    inPHOS++ ;
  }
  
  FillHistogram("hCenPHOS",fCentrality,inPHOS) ;
  
  //pi0
  for (Int_t i1=0; i1<inPHOS-1; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1) ;
    Int_t sm1 = ph1->Module();

    for (Int_t i2=i1+1; i2<inPHOS; i2++) {
      AliCaloPhoton * ph2=(AliCaloPhoton*)fPHOSEvent->At(i2) ;

      Int_t sm2 = ph2->Module();      
      TLorentzVector p12  = *ph1  + *ph2;
      
      snprintf(key,55,"hPi0All_cen%d",fCenBin) ;
      FillHistogram(key,p12.M() ,p12.Pt()) ; 

      if(sm1==sm2) {
	snprintf(key,55,"hPi0AllSM%d_cen%d",sm1,fCenBin) ;
	 FillHistogram(key,p12.M() ,p12.Pt()) ; 
      }
      
    } // end of loop i2
  } // end of loop i1
  
  //now mixed
  for (Int_t i1=0; i1<inPHOS; i1++) {
    AliCaloPhoton * ph1=(AliCaloPhoton*)fPHOSEvent->At(i1) ;
    Int_t sm1 = ph1->Module();

    for(Int_t ev=0; ev<prevPHOS->GetSize();ev++){
      TClonesArray * mixPHOS = static_cast<TClonesArray*>(prevPHOS->At(ev)) ;

      for(Int_t i2=0; i2<mixPHOS->GetEntriesFast();i2++){
	AliCaloPhoton * ph2=(AliCaloPhoton*)mixPHOS->At(i2) ;

	Int_t sm2 = ph2->Module();
	TLorentzVector p12  = *ph1  + *ph2;
	
	snprintf(key,55,"hMiPi0All_cen%d",fCenBin) ;
	FillHistogram(key,p12.M() ,p12.Pt()) ;

	if(sm1==sm2) {
	  snprintf(key,55,"hMiPi0AllSM%d_cen%d",sm1,fCenBin) ;
	  FillHistogram(key,p12.M() ,p12.Pt()) ; 
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

//_____________________________________________________________________________
void AliAnalysisTaskPHOSPbPbQA::FillHistogram(const char * key,Double_t x)const{
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
void AliAnalysisTaskPHOSPbPbQA::FillHistogram(const char * key,Double_t x,Double_t y)const{
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
