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
//_________________________________________________________________________
// Implementation version 1 of algorithm class to construct PHOS track segments
// Associates EMC and PPSD clusters
// Unfolds the EMC cluster   
//                  
//*-- Author: Dmitri Peressounko (RRC Ki & SUBATECH)
//

// --- ROOT system ---
#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TSystem.h"
#include "TBenchmark.h"
// --- Standard library ---

#include <iostream.h>
#include <iomanip.h>

// --- AliRoot header files ---

#include "AliPHOSTrackSegmentMakerv1.h"
#include "AliPHOSClusterizerv1.h"
#include "AliPHOSTrackSegment.h"
#include "AliPHOSCpvRecPoint.h"
#include "AliPHOSPpsdRecPoint.h"
#include "AliPHOSLink.h"
#include "AliPHOSv0.h"
#include "AliRun.h"

ClassImp( AliPHOSTrackSegmentMakerv1) 


//____________________________________________________________________________
 AliPHOSTrackSegmentMakerv1::AliPHOSTrackSegmentMakerv1() : AliPHOSTrackSegmentMaker()
{
  // ctor
  SetTitle("version 1") ;
  SetName("AliPHOSTrackSegmentMaker") ;
  fR0 = 10. ;   
  fEmcFirst = 0 ;    
  fEmcLast  = 0 ;   
  fCpvFirst = 0 ;   
  fCpvLast  = 0 ;   
  fPpsdFirst= 0 ;   
  fPpsdLast = 0 ;   
  fLinkLowArray = 0 ;
  fLinkUpArray  = 0 ;
  fIsInitialized = kFALSE ;
}
//____________________________________________________________________________
 AliPHOSTrackSegmentMakerv1::  AliPHOSTrackSegmentMakerv1(const char* headerFile, const char* branchTitle): AliPHOSTrackSegmentMaker()
{
  // ctor
  SetTitle("version 1") ;
  SetName("AliPHOSTrackSegmentMaker") ;
  fR0 = 10. ;   
  fEmcFirst = 0 ;    
  fEmcLast  = 0 ;   
  fCpvFirst = 0 ;   
  fCpvLast  = 0 ;   
  fPpsdFirst= 0 ;   
  fPpsdLast = 0 ;   

  fHeaderFileName = headerFile ;
  fRecPointsBranchTitle = branchTitle ;
    
  TFile * file = (TFile*) gROOT->GetFile(fHeaderFileName.Data() ) ;
  
  if(file == 0){
    file = new TFile(fHeaderFileName.Data(),"update") ;
    gAlice = (AliRun *) file->Get("gAlice") ;
  }
  
  AliPHOS * phos = (AliPHOS *) gAlice->GetDetector("PHOS") ;    
  fGeom  = AliPHOSGeometry::GetInstance(phos->GetGeometry()->GetName(),phos->GetGeometry()->GetTitle() );
  
  fEmcRecPoints = new TObjArray(200) ;
  fCpvRecPoints = new TObjArray(200) ;
  fClusterizer  = new AliPHOSClusterizerv1() ;
  
  fTrackSegments = new TClonesArray("AliPHOSTrackSegment",200) ;
  
  fLinkLowArray = new TClonesArray("AliPHOSLink", 1000);
  fLinkUpArray  = new TClonesArray("AliPHOSLink", 1000); 
  
  fIsInitialized = kTRUE ;

}
//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::Init(){

  if(!fIsInitialized){
    if(fHeaderFileName.IsNull())
      fHeaderFileName = "galice.root" ;
    
    
    TFile * file = (TFile*) gROOT->GetFile(fHeaderFileName.Data() ) ;
    
    if(file == 0){
      file = new TFile(fHeaderFileName.Data(),"update") ;
      gAlice = (AliRun *) file->Get("gAlice") ;
    }
    
    AliPHOS * phos = (AliPHOS *) gAlice->GetDetector("PHOS") ;    
    fGeom  = AliPHOSGeometry::GetInstance(phos->GetGeometry()->GetName(),phos->GetGeometry()->GetTitle() );


    fEmcRecPoints = new TObjArray(200) ;
    fCpvRecPoints = new TObjArray(200) ;
    fClusterizer  = new AliPHOSClusterizerv1() ;

    
    fTrackSegments = new TClonesArray("AliPHOSTrackSegment",200) ;

    fLinkLowArray = new TClonesArray("AliPHOSLink", 1000);
    fLinkUpArray  = new TClonesArray("AliPHOSLink", 1000); 
    
    fIsInitialized = kTRUE ;
   }
}

//____________________________________________________________________________
 AliPHOSTrackSegmentMakerv1::~AliPHOSTrackSegmentMakerv1()
{ 
  // dtor
  if(fLinkLowArray) delete fLinkLowArray ;
  if(fLinkUpArray)  delete fLinkUpArray  ;
}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::FillOneModule()
{
  // Finds bounds in which clusters from one PHOS module are
 

  //First EMC clusters
  Int_t totalEmc = fEmcRecPoints->GetEntriesFast() ;
  for(fEmcFirst = fEmcLast; (fEmcLast < totalEmc) &&  
	(((AliPHOSRecPoint *) fEmcRecPoints->At(fEmcLast))->GetPHOSMod() == fModule ); 
      fEmcLast ++)  ;
  
  
  //Now CPV clusters
  Int_t totalCpv = fCpvRecPoints->GetEntriesFast() ;

  if(fModule <= fGeom->GetNCPVModules()){ // in CPV geometry
    
    for(fCpvFirst = fCpvLast; (fCpvLast < totalCpv) && 
	  (((AliPHOSRecPoint *) fCpvRecPoints->At(fCpvLast))->GetPHOSMod() == fModule ); 
	fCpvLast ++) ;
    
    fPpsdFirst = fCpvLast ; //To avoid scanning RecPoints between fPpsdFirst and fPpsdLast
    fPpsdLast  = fCpvLast ; //and to be ready to switch to mixed geometry 
  }
  else{  //in PPSD geometry    
    fCpvLast = fPpsdLast ;
    //Upper layer first
    for(fCpvFirst = fCpvLast; (fCpvLast < totalCpv) &&  
	  (((AliPHOSPpsdRecPoint *) fCpvRecPoints->At(fCpvLast))->GetPHOSMod() == fModule ) &&
	  (((AliPHOSPpsdRecPoint *) fCpvRecPoints->At(fCpvLast))->GetUp()) ; 
	fCpvLast ++)  ;
    
    fPpsdLast= fCpvLast ;
    for(fPpsdFirst = fPpsdLast; (fPpsdLast < totalCpv)  &&
	  (((AliPHOSPpsdRecPoint *) fCpvRecPoints->At(fPpsdLast))->GetPHOSMod() == fModule ) &&
	  (!((AliPHOSPpsdRecPoint *) fCpvRecPoints->At(fPpsdLast))->GetUp()) ; 
	fPpsdLast ++) ;
  }
    
}
//____________________________________________________________________________
Float_t  AliPHOSTrackSegmentMakerv1::GetDistanceInPHOSPlane(AliPHOSEmcRecPoint * emcClu,AliPHOSRecPoint * cpvClu, Bool_t &toofar)
{
  // Calculates the distance between the EMC RecPoint and the PPSD RecPoint
  //clusters are sorted in "rows" and "columns" of width 1 cm
  Float_t delta = 1 ;  // Width of the rows in sorting of RecPoints (in cm)
                       // if you change this value, change it as well in xxxRecPoint::Compare()
  Float_t r = fR0 ;
 
  TVector3 vecEmc ;
  TVector3 vecCpv ;
  
  emcClu->GetLocalPosition(vecEmc) ;
  cpvClu->GetLocalPosition(vecCpv)  ; 

  if(emcClu->GetPHOSMod() == cpvClu->GetPHOSMod()){ 
    if(vecCpv.X() <= vecEmc.X() + fR0 + 2*delta ){ 

      vecCpv = vecCpv  - vecEmc ; 
      r = vecCpv.Mag() ;
      toofar = kFALSE ;

    } // if  xPpsd >= xEmc + ...
    else 
      toofar = kTRUE ;
  } 
  else 
    toofar = kTRUE ;

  //toofar = kFALSE ;
 
  
  return r ;
}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::MakeLinks()
{ 
  // Finds distances (links) between all EMC and PPSD clusters, which are not further apart from each other than fR0 
  
  fLinkUpArray->Clear() ;    
  fLinkLowArray->Clear() ;

  AliPHOSRecPoint * ppsd ; 
  AliPHOSRecPoint * cpv ;
  AliPHOSEmcRecPoint * emcclu ;

  Int_t iLinkLow = 0 ;
  Int_t iLinkUp  = 0 ;
  
  Int_t iEmcRP;
  for(iEmcRP = fEmcFirst; iEmcRP < fEmcLast; iEmcRP++ ) {
    emcclu = (AliPHOSEmcRecPoint *) fEmcRecPoints->At(iEmcRP) ;

    Bool_t toofar ;    
    Int_t iPpsd ;
    for(iPpsd = fPpsdFirst; iPpsd < fPpsdLast;iPpsd++ ) {
      
      ppsd = (AliPHOSRecPoint *) fCpvRecPoints->At(iPpsd) ;
      Float_t r = GetDistanceInPHOSPlane(emcclu, ppsd, toofar) ;

      if(toofar) 
	break ;	 
      if(r < fR0)
	new ((*fLinkLowArray)[iLinkLow++])  AliPHOSLink(r, iEmcRP, iPpsd) ;
    }
    
    Int_t iCpv = 0 ;    
    for(iCpv = fCpvFirst; iCpv < fCpvLast;iCpv++ ) { 
      
      cpv = (AliPHOSRecPoint *) fCpvRecPoints->At(iCpv) ;
      Float_t r = GetDistanceInPHOSPlane(emcclu, cpv, toofar) ;
      
      if(toofar)
	break ;	 
      if(r < fR0) { 
	new ((*fLinkUpArray)[iLinkUp++])  AliPHOSLink(r, iEmcRP, iCpv) ;
      }      
    }
  } 
  
  fLinkLowArray->Sort() ; //first links with smallest distances
  fLinkUpArray->Sort() ;
}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::MakePairs()
{ 
  
  //Make arrays to mark clusters already chousen
  Int_t * emcExist = 0;
  if(fEmcLast > fEmcFirst)
    emcExist = new Int_t[fEmcLast-fEmcFirst] ;
  
  Int_t index;
  for(index = 0; index <fEmcLast-fEmcFirst; index ++)
    emcExist[index] = 1 ;
  
  Bool_t * cpvExist = 0;
  if(fCpvLast > fCpvFirst)
    cpvExist = new Bool_t[fCpvLast-fCpvFirst] ;
  for(index = 0; index <fCpvLast-fCpvFirst; index ++)
    cpvExist[index] = kTRUE ;
  
  Bool_t * ppsdExist = 0;
  if(fPpsdLast > fPpsdFirst)
    ppsdExist = new Bool_t[fPpsdLast-fPpsdFirst] ;
  for(index = 0; index <fPpsdLast-fPpsdFirst; index ++)
    ppsdExist[index] = kTRUE ;
  
  // Finds the smallest links and makes pairs of CPV and EMC clusters with smallest distance 
  TIter nextLow(fLinkLowArray) ;
  TIter nextUp(fLinkUpArray) ;
  
  AliPHOSLink * linkLow ;
  AliPHOSLink * linkUp ;


  AliPHOSRecPoint * nullpointer = 0 ;

  while ( (linkLow =  (AliPHOSLink *)nextLow() ) ){
  
    if( (emcExist[linkLow->GetEmc()-fEmcFirst]> 0) && ppsdExist[linkLow->GetPpsd()-fPpsdFirst]  ){ // RecPoints not removed yet 
      new ((*fTrackSegments)[fNTrackSegments]) AliPHOSTrackSegment((AliPHOSEmcRecPoint *) fEmcRecPoints->At(linkLow->GetEmc()), 
						 nullpointer, 
						(AliPHOSPpsdRecPoint *)fCpvRecPoints->At(linkLow->GetPpsd()) ) ;
	 
      ((AliPHOSTrackSegment* )fTrackSegments->At(fNTrackSegments))->SetIndexInList(fNTrackSegments);    
      //replace index of emc to negative and shifted index of TS      
      emcExist[linkLow->GetEmc()-fEmcFirst] = -2 - fNTrackSegments ;  
      //mark ppsd as used
      ppsdExist[linkLow->GetPpsd()-fPpsdFirst] = kFALSE ; 
      fNTrackSegments++ ;
    } 
  } 
	 

  while ( (linkUp =  (AliPHOSLink *)nextUp() ) ){  
    if(emcExist[linkUp->GetEmc()-fEmcFirst] != -1){ //without ppsd Up yet 

      if(cpvExist[linkUp->GetPpsd()-fCpvFirst]){ //CPV still exist
	
	if(emcExist[linkUp->GetEmc()-fEmcFirst] > 0){ //without ppsd Low => create new TS

	  new ((* fTrackSegments)[fNTrackSegments]) AliPHOSTrackSegment((AliPHOSEmcRecPoint *) fEmcRecPoints->At(linkUp->GetEmc()) , 
								      (AliPHOSPpsdRecPoint *)fCpvRecPoints->At(linkUp->GetPpsd()), 
								      nullpointer) ;
	  ((AliPHOSTrackSegment *) fTrackSegments->At(fNTrackSegments))->SetIndexInList(fNTrackSegments);
	  fNTrackSegments++ ;
	}
	else{ // append ppsd Up to existing TS
	  ((AliPHOSTrackSegment *)fTrackSegments->At(-2-emcExist[linkUp->GetEmc()-fEmcFirst]))->SetCpvRecPoint((AliPHOSCpvRecPoint *)fCpvRecPoints->At(linkUp->GetPpsd()));
	}

	emcExist[linkUp->GetEmc()-fEmcFirst] = -1 ; //Mark emc  that Cpv was found 
	//mark CPV recpoint as already used 
        cpvExist[linkUp->GetPpsd()-fCpvFirst] = kFALSE ;
      } //if ppsdUp still exist
    } 
  }	 

  //look through emc recPoints left without CPV/PPSD
  if(emcExist){ //if there is emc rec point
    Int_t iEmcRP ;
    for(iEmcRP = 0; iEmcRP < fEmcLast-fEmcFirst  ; iEmcRP++ ){
      if(emcExist[iEmcRP] > 0 ){
	new ((*fTrackSegments)[fNTrackSegments])  AliPHOSTrackSegment((AliPHOSEmcRecPoint *)fEmcRecPoints->At(iEmcRP+fEmcFirst), 
								    nullpointer, 
								    nullpointer ) ;
	((AliPHOSTrackSegment *) fTrackSegments->At(fNTrackSegments))->SetIndexInList(fNTrackSegments);
	fNTrackSegments++;    
      } 
    }
  }
  
}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::Exec(Option_t * option)
{
  // Makes the track segments out of the list of EMC and PPSD Recpoints and stores them in a list

  if(! fIsInitialized) Init() ;

  if(strstr(option,"tim"))
    gBenchmark->Start("PHOSTSMaker");  

  Int_t nEvents = (Int_t) gAlice->TreeE()->GetEntries() ;
  
  for(fEvent = 0;fEvent< nEvents; fEvent++){
    if(!ReadRecPoints())  //reads RecPoints for event fEvent
      return;
    
    for(fModule = 1; fModule <= fGeom->GetNModules() ; fModule++ ){
      
      FillOneModule() ; 
      
      MakeLinks() ;
      
      MakePairs() ;
      
    }

    WriteTrackSegments() ;
    if(strstr(option,"deb"))
      PrintTrackSegments(option) ;
  }

  if(strstr(option,"tim")){
    gBenchmark->Stop("PHOSTSMaker");
    cout << "AliPHOSTSMaker:" << endl ;
    cout << "  took " << gBenchmark->GetCpuTime("PHOSTSMaker") << " seconds for making TS " 
	 <<  gBenchmark->GetCpuTime("PHOSTSMaker")/nEvents << " seconds per event " << endl ;
    cout << endl ;
  }


}
//____________________________________________________________________________
void AliPHOSTrackSegmentMakerv1::Print(Option_t * option)const {
  if(fIsInitialized){
    cout <<  "======== AliPHOSTrackSegmentMakerv1 ========" << endl ;
    cout <<  "Making Track segments "<< endl ;
    cout <<  "    Headers file: " << fHeaderFileName.Data() << endl ;
    cout <<  "    RecPoints branch file name: " <<fRecPointsBranchTitle.Data() << endl ;
    cout <<  "    TrackSegments Branch file name: " << fTSBranchTitle.Data() << endl ;
    cout <<  "with parameters: " << endl ;
    cout <<  "    Maximal EMC - CPV (PPSD) distance " << fR0 << endl ;
    cout <<  "============================================" << endl ;
  }
  else
    cout << "AliPHOSTrackSegmentMakerv1 not initialized " << endl ;
}
//____________________________________________________________________________
Bool_t AliPHOSTrackSegmentMakerv1::ReadRecPoints(){

  //Make some initializations 
  fEmcRecPoints->Clear() ;
  fCpvRecPoints->Clear() ;
  fTrackSegments->Clear() ;
  fNTrackSegments = 0 ;
  fEmcFirst = 0 ;    
  fEmcLast  = 0 ;   
  fCpvFirst = 0 ;   
  fCpvLast  = 0 ;   
  fPpsdFirst= 0 ;   
  fPpsdLast = 0 ;   


  gAlice->GetEvent(fEvent) ;

  // Get TreeR header from file
  char treeName[20]; 
  sprintf(treeName,"TreeR%d",fEvent);

  if(gAlice->TreeR()==0){
    cout << "Error in AliPHOSTrackSegmentMakerv1 : no "<<treeName << endl  ;
    cout << "   Do nothing " << endl ;
    return kFALSE ;
  }

  char * emcBranchName = new char[30];
  // sprintf(emcBranchName,"PHOSEmcRP%d",fEvent);
  sprintf(emcBranchName,"PHOSEmcRP");
  char * cpvBranchName = new char[30];
  //  sprintf(cpvBranchName,"PHOSCpvRP%d",fEvent);
  sprintf(cpvBranchName,"PHOSCpvRP");
  char * cluBranchName = new char[30];
  //  sprintf(cluBranchName,"AliPHOSClusterizer%d",fEvent);
  sprintf(cluBranchName,"AliPHOSClusterizer");

  TBranch * emcBranch = 0;
  TBranch * cpvBranch = 0;
  TBranch * clusterizerBranch = 0;

  TObjArray * branches = gAlice->TreeR()->GetListOfBranches() ;
  Int_t ibranch;
  Bool_t emcNotFound = kTRUE ;
  Bool_t cpvNotFound = kTRUE ;  
  Bool_t clusterizerNotFound = kTRUE ;
  
  for(ibranch = 0;ibranch <branches->GetEntries();ibranch++){

    if(emcNotFound){
      emcBranch=(TBranch *) branches->At(ibranch) ;
      if( fRecPointsBranchTitle.CompareTo(emcBranch->GetTitle())==0 )
	if( strcmp(emcBranch->GetName(),emcBranchName) == 0) {
	  emcNotFound = kFALSE ;
	}
    }
    
    if(cpvNotFound){
      cpvBranch=(TBranch *) branches->At(ibranch) ;
      if( fRecPointsBranchTitle.CompareTo(cpvBranch->GetTitle())==0 )
	if( strcmp(cpvBranch->GetName(),cpvBranchName) == 0) 
	  cpvNotFound = kFALSE ;
    }
    
    if(clusterizerNotFound){
      clusterizerBranch = (TBranch *) branches->At(ibranch) ;
      if( fRecPointsBranchTitle.CompareTo(clusterizerBranch->GetTitle()) == 0)
	if( strcmp(clusterizerBranch->GetName(),cluBranchName) == 0) 
	  clusterizerNotFound = kFALSE ;
    }
    
  }

  if(clusterizerNotFound || emcNotFound || cpvNotFound){
    cout << "AliPHOSTrackSegmentMakerv1: " << endl ;
    cout << "    Can't find Branch with RecPoints or Clusterizer " ;
    cout << "    Do nothing" <<endl  ;
    return kFALSE ;
  }
  
  emcBranch->SetAddress(&fEmcRecPoints) ;
  cpvBranch->SetAddress(&fCpvRecPoints) ;
  clusterizerBranch->SetAddress(&fClusterizer) ;
  
  gAlice->TreeR()->GetEvent(0) ;
  
  delete emcBranchName;
  delete cpvBranchName;
  delete cluBranchName;

  return kTRUE ;
  
}
//____________________________________________________________________________
void AliPHOSTrackSegmentMakerv1::WriteTrackSegments(){

  char treeName[20]; 
  sprintf(treeName,"TreeR%d",fEvent);
  
  
  //First, check, if branches already exist
  TBranch * tsMakerBranch = 0;
  TBranch * tsBranch = 0;
  
  TObjArray * branches = gAlice->TreeR()->GetListOfBranches() ;
  Int_t ibranch;
  Bool_t tsMakerNotFound = kTRUE ;
  Bool_t tsNotFound = kTRUE ;
  
  for(ibranch = 0;(ibranch <branches->GetEntries())&&(tsMakerNotFound||tsNotFound);ibranch++){
    if(tsMakerNotFound){
      tsMakerBranch=(TBranch *) branches->At(ibranch) ;
      if( (strcmp(tsMakerBranch->GetName(),"AliPHOSTrackSegmentMaker") == 0) &&
	  (fTSBranchTitle.CompareTo( tsMakerBranch->GetTitle())==0 ))
	tsMakerNotFound = kFALSE ;
    }
    if(tsNotFound){
      tsBranch=(TBranch *) branches->At(ibranch) ;
      if( (strcmp(tsBranch->GetName(),"PHOSTS") == 0)  &&
	  (fTSBranchTitle.CompareTo( tsBranch->GetTitle())==0 ))
	tsNotFound = kFALSE ;
    }
  }

  if(!(tsMakerNotFound && tsNotFound )){ 
    cout << "AliPHOSTrackSegmentMakerv1 error:"<< endl ;
    cout << "       Branches PHOSTS and AliPHOSTrackSegementMaker " << endl ;
    cout << "       with title '"<<fTSBranchTitle.Data() << "' already exist " << endl ;
    cout << "       can not overwrite " << endl ;
    return ;
  }

  //Make branch in TreeR for TrackSegments 
  char * filename = 0;
  if(gSystem->Getenv("CONFIG_SPLIT_FILE")!=0){   //generating file name
    filename = new char[strlen(gAlice->GetBaseFile())+20] ;
    sprintf(filename,"%s/PHOS.Reco.root",gAlice->GetBaseFile()) ; 
  }

  TDirectory *cwd = gDirectory;
  
  //First TS
  Int_t bufferSize = 32000 ;    
  tsBranch = gAlice->TreeR()->Branch("PHOSTS",&fTrackSegments,bufferSize);
  tsBranch->SetTitle(fTSBranchTitle.Data());
  if (filename) {
    tsBranch->SetFile(filename);
    TIter next( tsBranch->GetListOfBranches());
    while ((tsBranch=(TBranch*)next())) {
      tsBranch->SetFile(filename);
    }   
    cwd->cd();
  } 
  
  //Second -TSMaker
  Int_t splitlevel = 0 ;
  AliPHOSTrackSegmentMakerv1 * ts = this ;
  tsMakerBranch = gAlice->TreeR()->Branch("AliPHOSTrackSegmentMaker","AliPHOSTrackSegmentMakerv1",
					  &ts,bufferSize,splitlevel);
  tsMakerBranch->SetTitle(fTSBranchTitle.Data());
  if (filename) {
    tsMakerBranch->SetFile(filename);
    TIter next( tsMakerBranch->GetListOfBranches());
    while ((tsMakerBranch=(TBranch*)next())) {
      tsMakerBranch->SetFile(filename);
    }   
    cwd->cd();
  } 
  
  gAlice->TreeR()->Fill() ;    
  gAlice->TreeR()->Write(0,kOverwrite) ;  
  
}


//____________________________________________________________________________
void AliPHOSTrackSegmentMakerv1::PrintTrackSegments(Option_t * option){
  
  cout << "AliPHOSTrackSegmentMakerv1: " << endl ;
  cout << "       Found " << fTrackSegments->GetEntriesFast() << "  trackSegments " << endl ;
  
  if(strstr(option,"all")) {  // printing found TS
    cout << "TrackSegment # " << "    EMC RP#    " << "    CPV RP#    " << "     PPSD RP#" << endl ; 
    
    Int_t index;
    for (index = 0 ; index <fTrackSegments->GetEntriesFast() ; index++) {
      AliPHOSTrackSegment * ts = (AliPHOSTrackSegment * )fTrackSegments->At(index) ; 
      cout<<"   "<< setw(4) << ts->GetIndexInList() << "            " 
	  <<setw(4) << ts->GetEmcIndex()<< "            " 
	  <<setw(4) << ts->GetCpvIndex()<< "            " 
	  <<setw(4) << ts->GetPpsdIndex()<< endl ;
    }	
    
    cout << "-------------------------------------------------------"<< endl ;
  }
}
//____________________________________________________________________________
void AliPHOSTrackSegmentMakerv1::SetRecPointsBranch(const char * title){
  //set the title of RecPoints 
    fRecPointsBranchTitle = title ;

}
//____________________________________________________________________________
void AliPHOSTrackSegmentMakerv1::SetTrackSegmentsBranch(const char * title){

    fTSBranchTitle = title ; 
}
