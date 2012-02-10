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

#include "AliAnalysisManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskPHOSTriggerQA.h"
#include "AliESDCaloCluster.h"
#include "AliPHOSGeometry.h"
#include "AliESDEvent.h"
#include "AliESDCaloCells.h"
#include "AliLog.h"
#include "TObjArray.h"
#include "TList.h"
#include "TH1.h"
#include "TH2.h"

// QA of PHOS Trigger data.
//...
// Author: Boris Polishchuk 
// Date  : 06.02.2012

ClassImp(AliAnalysisTaskPHOSTriggerQA)

//________________________________________________________________________
AliAnalysisTaskPHOSTriggerQA::AliAnalysisTaskPHOSTriggerQA() : AliAnalysisTaskSE(),
  fOutputContainer(0),fPHOSGeo(0),fEventCounter(0)
{
  //Default constructor.  
  // Initialize the PHOS geometry 
  fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP") ;
  
}

//________________________________________________________________________
AliAnalysisTaskPHOSTriggerQA::AliAnalysisTaskPHOSTriggerQA(const char *name) 
: AliAnalysisTaskSE(name),
  fOutputContainer(0),fPHOSGeo(0),fEventCounter(0)
{
  
  // Output slots #0 write into a TH1 container
  DefineOutput(1,TList::Class());

  // Initialize the PHOS geometry
  fPHOSGeo = AliPHOSGeometry::GetInstance("IHEP") ;

}

//________________________________________________________________________
void AliAnalysisTaskPHOSTriggerQA::UserCreateOutputObjects()
{
  // Create histograms
  // Called once
  
  // ESD histograms
  if(fOutputContainer != NULL){
    delete fOutputContainer;
  }

  fOutputContainer = new TList();
  fOutputContainer->SetOwner(kTRUE);

  //Bin 1: total number of processed events.
  //Bin 2: number of events contained PHOS trigger digits.
  fOutputContainer->Add(new TH1F("hNev","Number of events",10,0.,10.));
  
  char key[55],titl[55];
  Int_t nCols = 56, nRows = 64;
  
  Int_t    nPtPhot  = 1000 ;
  Double_t ptPhotMax =  100. ;

  Int_t nTrMax = 1200;
  Float_t trMax = 600.;

  fOutputContainer->Add(new TH1F("hNtr","Number of fired 4x4 regions",nTrMax,0.,trMax));

  for(Int_t sm=1; sm<4; sm++) {

    snprintf(key,55,"hNtrSM%d",sm);
    snprintf(titl,55,"Number of fired 4x4 regions in SM%d",sm);
    fOutputContainer->Add(new TH1F(key,titl,nTrMax/3,0.,trMax/3));
    
    snprintf(key,55,"h4x4SM%d",sm);
    snprintf(titl,55,"SM%d 4x4 occupancy",sm);
    fOutputContainer->Add(new TH2F(key,titl,nRows,0.,nRows,nCols,0.,nCols));

    snprintf(key,55,"hCluSM%d",sm);
    snprintf(titl,55,"SM%d cluster occupancy",sm);
    fOutputContainer->Add(new TH2F(key,titl,nRows,0.,nRows,nCols,0.,nCols));
 
    snprintf(key,55,"hCluTSM%d",sm);
    snprintf(titl,55,"SM%d triggered cluster occupancy",sm);
    fOutputContainer->Add(new TH2F(key,titl,nRows,0.,nRows,nCols,0.,nCols));
   
    snprintf(key,55,"hPhotAllSM%d",sm);
    snprintf(titl,55,"SM%d cluster energy",sm);
    fOutputContainer->Add(new TH1F(key,titl,nPtPhot,0.,ptPhotMax));

    snprintf(key,55,"hPhotTrigSM%d",sm);
    snprintf(titl,55,"SM%d triggered cluster energy",sm);
    fOutputContainer->Add(new TH1F(key,titl,nPtPhot,0.,ptPhotMax));
  }
  
  PostData(1, fOutputContainer);
  
}

//________________________________________________________________________
void AliAnalysisTaskPHOSTriggerQA::UserExec(Option_t *) 
{
  // Main loop, called for each event
  // Analyze ESD/AOD  
  
  AliESDEvent *event = dynamic_cast<AliESDEvent*>(InputEvent());
  
  if (!event) {
    Printf("ERROR: Could not retrieve event");
    PostData(1, fOutputContainer);
    return;
  }
  
  FillHistogram("hNev",0.); // all events
  fEventCounter++;
  
  AliESDCaloTrigger* trgESD = event->GetCaloTrigger("PHOS");
  trgESD->Reset();
  
  if (!trgESD->GetEntries()) {
    PostData(1, fOutputContainer);
    return;
  }
  
  FillHistogram("hNev",1.); // triggered events
  FillHistogram("hNtr",trgESD->GetEntries());

  TString trigClasses = event->GetFiredTriggerClasses();
  printf("\nEvent %d: %d non-zero trigger digits %s\n",
	 fEventCounter,trgESD->GetEntries(),trigClasses.Data());

  // Get PHOS rotation matrices from ESD and set them to the PHOS geometry
  char key[55] ;  
  
  if(fEventCounter == 0) {
    for(Int_t mod=0; mod<5; mod++) {
      if(!event->GetPHOSMatrix(mod)) continue;
      fPHOSGeo->SetMisalMatrix(event->GetPHOSMatrix(mod),mod) ;
    }
  }
  
  Int_t multClu = event->GetNumberOfCaloClusters();
  AliESDCaloCells *phsCells = event->GetPHOSCells();
 
  Int_t inPHOS[3] = {};

  //Loop over 4x4 fired regions
  while(trgESD->Next()) {

    Int_t tmod,tabsId; // "Online" module number, bottom-left 4x4 edge cell absId
    trgESD->GetPosition(tmod,tabsId);
    
    Int_t trelid[4] ;
    fPHOSGeo->AbsToRelNumbering(tabsId,trelid);

    snprintf(key,55,"h4x4SM%d",trelid[0]);
    FillHistogram(key,trelid[2]-1,trelid[3]-1);
    
    inPHOS[trelid[0]-1]++;

    for (Int_t i=0; i<multClu; i++) {
      
      AliESDCaloCluster *c1 = event->GetCaloCluster(i);
      if(!c1->IsPHOS()) continue;
      
      Int_t maxId, relid[4];
      MaxEnergyCellPos(phsCells,c1,maxId);
      
      fPHOSGeo->AbsToRelNumbering(maxId, relid);
      snprintf(key,55,"hPhotAllSM%d",relid[0]);
      FillHistogram(key,c1->E());
      
      snprintf(key,55,"hCluSM%d",relid[0]);
      FillHistogram(key,relid[2]-1,relid[3]-1);
      
      if( Matched(trelid,relid) ) {

	snprintf(key,55,"hPhotTrigSM%d",relid[0]);
	FillHistogram(key,c1->E());

	snprintf(key,55,"hCluTSM%d",relid[0]);
	FillHistogram(key,relid[2]-1,relid[3]-1);
	
	continue;
      }
      
    }
  } //while(trgESD->Next())
  
  for(Int_t sm=1; sm<4; sm++) {
    snprintf(key,55,"hNtrSM%d",sm);
    if(inPHOS[sm-1]) FillHistogram(key,inPHOS[sm-1]);    
  }
  
  // Post output data.
  PostData(1, fOutputContainer);
  
}

//_____________________________________________________________________________
void AliAnalysisTaskPHOSTriggerQA::FillHistogram(const char * key,Double_t x)const{
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
void AliAnalysisTaskPHOSTriggerQA::FillHistogram(const char * key,Double_t x,Double_t y)const{
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
void AliAnalysisTaskPHOSTriggerQA::MaxEnergyCellPos(AliESDCaloCells *cells, AliESDCaloCluster* clu, Int_t& maxId)
{  
  Double_t eMax = -111;
  
  for (Int_t iDig=0; iDig< clu->GetNCells(); iDig++) {
    Int_t cellAbsId = clu->GetCellAbsId(iDig);
    Double_t eCell = cells->GetCellAmplitude(cellAbsId)*clu->GetCellAmplitudeFraction(iDig);
    if(eCell>eMax)  { 
      eMax = eCell; 
      maxId = cellAbsId;
    }
  }
  
}

//_____________________________________________________________________________
Bool_t AliAnalysisTaskPHOSTriggerQA::Matched(Int_t *trig_relid, Int_t *cluster_relid)
{
  //Returns kTRUE if cluster position coincides with 4x4 position.

  if( trig_relid[0] != cluster_relid[0] )            return kFALSE; // different modules!
  if( TMath::Abs(trig_relid[2]-cluster_relid[2])>3 ) return kFALSE; // X-distance too large! 
  if( TMath::Abs(trig_relid[3]-cluster_relid[3])>3 ) return kFALSE; // Z-distance too large!

  return kTRUE;
}
