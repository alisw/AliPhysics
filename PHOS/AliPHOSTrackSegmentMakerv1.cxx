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
// Track segment for PHOS is list of 
//        EMC RecPoint + (possibly) CPV RecPoint + (possibly) PPSD RecPoint
// To find TrackSegments we do the following: 
//  for each EMC RecPoints we look at
//   CPV/PPSD RecPoints in the radious fR0. 
//  If there is such a CPV RecPoint, 
//   we make "Link" it is just indexes of EMC and CPV/PPSD RecPoint and distance
//   between them in the PHOS plane. 
//  Then we sort "Links" and starting from the 
//   least "Link" pointing to the unassined EMC and CPV RecPoints assing them to 
//   new TrackSegment. 
// If there is no CPV/PPSD RecPoint we make TrackSegment 
// consisting from EMC alone. There is no TrackSegments without EMC RecPoint.
//
// In principle this class should be called from AliPHOSReconstructioner, but 
// one can use it as well in standalone mode.
// Use  case:
//  root [0] AliPHOSTrackSegmentMakerv1 * t = new AliPHOSTrackSegmentMaker("galice.root")
//  Warning in <TDatabasePDG::TDatabasePDG>: object already instantiated
//  root [1] t->ExecuteTask()
//  root [2] t->SetMaxEmcPpsdDistance(5)
//  root [3] t->SetTrackSegmentsBranch("max distance 5 cm")
//  root [4] t->ExecuteTask("deb all time") 
//                 
//*-- Author: Dmitri Peressounko (RRC Ki & SUBATECH)
//

// --- ROOT system ---
#include "TROOT.h"
#include "TFile.h"
#include "TFolder.h"
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
#include "AliPHOSGetter.h"
#include "AliPHOS.h"
#include "AliRun.h"

ClassImp( AliPHOSTrackSegmentMakerv1) 


//____________________________________________________________________________
  AliPHOSTrackSegmentMakerv1::AliPHOSTrackSegmentMakerv1() : AliPHOSTrackSegmentMaker()
{
  // default ctor (to be used mainly by Streamer)

  fR0                       = 10. ;   
  fEmcFirst                 = 0 ;    
  fEmcLast                  = 0 ;   
  fCpvFirst                 = 0 ;   
  fCpvLast                  = 0 ;   
  fPpsdFirst                = 0 ;   
  fPpsdLast                 = 0 ;   
  fLinkLowArray             = 0 ;
  fLinkUpArray              = 0 ;
  fHeaderFileName           = "" ;
  fRecPointsBranchTitle     = "" ;
  fTrackSegmentsBranchTitle = "" ; 
}

//____________________________________________________________________________
 AliPHOSTrackSegmentMakerv1::AliPHOSTrackSegmentMakerv1(const char * headerFile, const char * name) : AliPHOSTrackSegmentMaker(headerFile, name)
{
  // ctor

  fR0        = 10. ;   
  fEmcFirst  = 0 ;    
  fEmcLast   = 0 ;   
  fCpvFirst  = 0 ;   
  fCpvLast   = 0 ;   
  fPpsdFirst = 0 ;   
  fPpsdLast  = 0 ;   

  fHeaderFileName           = GetTitle() ;
  fRecPointsBranchTitle     = GetName() ;
  fTrackSegmentsBranchTitle = GetName() ; 
  
  TString tempo(GetName()) ; 
  tempo.Append(Version()) ; 
  SetName(tempo.Data()) ; 

  Init() ;

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
  // Finds first and last indexes between which 
  // clusters from one PHOS module are
 
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ; 
  const AliPHOSGeometry * geom = gime->PHOSGeometry() ; 
  TObjArray * emcRecPoints = gime->EmcRecPoints() ; 
  TObjArray * cpvRecPoints = gime->CpvRecPoints() ; 

  //First EMC clusters
  Int_t totalEmc = emcRecPoints->GetEntriesFast() ;
  for(fEmcFirst = fEmcLast; (fEmcLast < totalEmc) &&  
	(((AliPHOSRecPoint *) emcRecPoints->At(fEmcLast))->GetPHOSMod() == fModule ); 
      fEmcLast ++)  ;
  
  
  //Now CPV clusters
  Int_t totalCpv = cpvRecPoints->GetEntriesFast() ;

  if(fModule <= geom->GetNCPVModules()){ // in CPV geometry
    
    for(fCpvFirst = fCpvLast; (fCpvLast < totalCpv) && 
	  (((AliPHOSRecPoint *) cpvRecPoints->At(fCpvLast))->GetPHOSMod() == fModule ); 
	fCpvLast ++) ;
    
    fPpsdFirst = fCpvLast ; //To avoid scanning RecPoints between fPpsdFirst and fPpsdLast
    fPpsdLast  = fCpvLast ; //and to be ready to switch to mixed geometry 
  }
  else{  //in PPSD geometry    
    fCpvLast = fPpsdLast ;
    //Upper layer first
    for(fCpvFirst = fCpvLast; (fCpvLast < totalCpv) &&  
	  (((AliPHOSPpsdRecPoint *) cpvRecPoints->At(fCpvLast))->GetPHOSMod() == fModule ) &&
	  (((AliPHOSPpsdRecPoint *) cpvRecPoints->At(fCpvLast))->GetUp()) ; 
	fCpvLast ++)  ;
    
    fPpsdLast= fCpvLast ;
    for(fPpsdFirst = fPpsdLast; (fPpsdLast < totalCpv)  &&
	  (((AliPHOSPpsdRecPoint *) cpvRecPoints->At(fPpsdLast))->GetPHOSMod() == fModule ) &&
	  (!((AliPHOSPpsdRecPoint *) cpvRecPoints->At(fPpsdLast))->GetUp()) ; 
	fPpsdLast ++) ;
  }
    
}

//____________________________________________________________________________
Float_t  AliPHOSTrackSegmentMakerv1::GetDistanceInPHOSPlane(AliPHOSEmcRecPoint * emcClu,AliPHOSRecPoint * cpvClu, Bool_t &toofar)const
{
  // Calculates the distance between the EMC RecPoint and the PPSD RecPoint
  // Clusters are sorted in "rows" and "columns" of width 1 cm

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
void  AliPHOSTrackSegmentMakerv1::Init()
{
  // Make all memory allocations that are not possible in default constructor
  
  if ( strcmp(GetTitle(), "") == 0 )
    SetTitle("galice.root") ;
  
  TString taskName(GetName()) ; 
  taskName.ReplaceAll(Version(), "") ;
  
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance(GetTitle(), taskName.Data()) ; 
  if ( gime == 0 ) {
    cerr << "ERROR: AliPHOSTrackSegmentMakerv1::Init -> Could not obtain the Getter object !" << endl ; 
    return ;
  } 
  
  fLinkLowArray = new TClonesArray("AliPHOSLink", 1000);
  fLinkUpArray  = new TClonesArray("AliPHOSLink", 1000); 
  
  //add Task to //YSAlice/tasks/Reconstructioner/PHOS
  TTask * aliceRe  = (TTask*)gROOT->FindObjectAny("YSAlice/tasks/Reconstructioner") ; 
  TTask * phosRe   = (TTask*)aliceRe->GetListOfTasks()->FindObject("PHOS") ;
  phosRe->Add(this) ; 
  // create a folder on the white board //YSAlice/WhiteBoard/RecPoints/PHOS/trackSegmentsName

  gime->Post(GetTitle(), "T", taskName.Data() ) ; 

}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::MakeLinks()const
{ 
  // Finds distances (links) between all EMC and PPSD clusters, 
  // which are not further apart from each other than fR0 
  // and sort them in accordance with this distance
  
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ; 
  TObjArray * emcRecPoints = gime->EmcRecPoints() ; 
  TObjArray * cpvRecPoints = gime->CpvRecPoints() ; 

  fLinkUpArray->Clear() ;    
  fLinkLowArray->Clear() ;

  AliPHOSRecPoint * ppsd ; 
  AliPHOSRecPoint * cpv ;
  AliPHOSEmcRecPoint * emcclu ;

  Int_t iLinkLow = 0 ;
  Int_t iLinkUp  = 0 ;
  
  Int_t iEmcRP;
  for(iEmcRP = fEmcFirst; iEmcRP < fEmcLast; iEmcRP++ ) {
    emcclu = (AliPHOSEmcRecPoint *) emcRecPoints->At(iEmcRP) ;

    Bool_t toofar ;    
    Int_t iPpsd ;
    for(iPpsd = fPpsdFirst; iPpsd < fPpsdLast;iPpsd++ ) {
      
      ppsd = (AliPHOSRecPoint *) cpvRecPoints->At(iPpsd) ;
      Float_t r = GetDistanceInPHOSPlane(emcclu, ppsd, toofar) ;

      if(toofar) 
	break ;	 
      if(r < fR0)
	new ((*fLinkLowArray)[iLinkLow++])  AliPHOSLink(r, iEmcRP, iPpsd) ;
    }
    
    Int_t iCpv = 0 ;    
    for(iCpv = fCpvFirst; iCpv < fCpvLast;iCpv++ ) { 
      
      cpv = (AliPHOSRecPoint *) cpvRecPoints->At(iCpv) ;
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
  // Using the previously made list of "links", we found the smallest link - i.e. 
  // link with the least distance between EMC and CPV and pointing to still 
  // unassigned RecParticles. We assign these RecPoints to TrackSegment and 
  // remove them from the list of "unassigned". 
  
  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ; 
  TObjArray * emcRecPoints = gime->EmcRecPoints() ; 
  TObjArray * cpvRecPoints = gime->CpvRecPoints() ; 
  TClonesArray * trackSegments = gime->TrackSegments() ; 
  
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
  
    if( (emcExist[linkLow->GetEmc()-fEmcFirst]> 0) && 
	ppsdExist[linkLow->GetPpsd()-fPpsdFirst]  ){ // RecPoints not removed yet 
      new ((*trackSegments)[fNTrackSegments]) AliPHOSTrackSegment((AliPHOSEmcRecPoint *) emcRecPoints->At(linkLow->GetEmc()), 
						 nullpointer, 
						(AliPHOSRecPoint *)cpvRecPoints->At(linkLow->GetPpsd()) ) ;
	 
      ((AliPHOSTrackSegment* )trackSegments->At(fNTrackSegments))->SetIndexInList(fNTrackSegments);    
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

	  new ((* trackSegments)[fNTrackSegments]) AliPHOSTrackSegment((AliPHOSEmcRecPoint *) emcRecPoints->At(linkUp->GetEmc()) , 
								      (AliPHOSRecPoint *)cpvRecPoints->At(linkUp->GetPpsd()), 
								      nullpointer) ;
	  ((AliPHOSTrackSegment *) trackSegments->At(fNTrackSegments))->SetIndexInList(fNTrackSegments);
	  fNTrackSegments++ ;
	}
	else{ // append ppsd Up to existing TS
	  ((AliPHOSTrackSegment *)trackSegments->At(-2-emcExist[linkUp->GetEmc()-fEmcFirst]))->SetCpvRecPoint((AliPHOSCpvRecPoint *)cpvRecPoints->At(linkUp->GetPpsd()));
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
	new ((*trackSegments)[fNTrackSegments])  AliPHOSTrackSegment((AliPHOSEmcRecPoint *)emcRecPoints->At(iEmcRP+fEmcFirst), 
								    nullpointer, 
								    nullpointer ) ;
	((AliPHOSTrackSegment *) trackSegments->At(fNTrackSegments))->SetIndexInList(fNTrackSegments);
	fNTrackSegments++;    
      } 
    }
  }
  
}

//____________________________________________________________________________
void  AliPHOSTrackSegmentMakerv1::Exec(Option_t * option)
{
  // STEERing method

  if( strcmp(GetName(), "")== 0 ) 
    Init() ;

  if(strstr(option,"tim"))
    gBenchmark->Start("PHOSTrackSegmentMakerv1");  
 
  if(strstr(option,"print")) {
    Print("") ; 
    return ; 
  }

  //check, if the branch with name of this" already exits?
  TObjArray * lob = (TObjArray*)gAlice->TreeR()->GetListOfBranches() ;
  TIter next(lob) ; 
  TBranch * branch = 0 ;  
  Bool_t phostsfound = kFALSE, tracksegmentmakerfound = kFALSE ; 
  
  TString taskName(GetName()) ; 
  taskName.ReplaceAll(Version(), "") ;

  while ( (branch = (TBranch*)next()) && (!phostsfound || !tracksegmentmakerfound) ) {
    if ( (strcmp(branch->GetName(), "PHOSTS")==0) && (strcmp(branch->GetTitle(), taskName.Data())==0) ) 
      phostsfound = kTRUE ;
    
    else if ( (strcmp(branch->GetName(), "AliPHOSTrackSegmentMaker")==0) && (strcmp(branch->GetTitle(), taskName.Data())==0) ) 
      tracksegmentmakerfound = kTRUE ; 
  }

  if ( phostsfound || tracksegmentmakerfound ) {
    cerr << "WARNING: AliPHOSTrackSegmentMakerv1::Exec -> TrackSegments and/or TrackSegmentMaker branch with name " 
	 << taskName.Data() << " already exits" << endl ;
    return ; 
  }       

  const AliPHOSGeometry * geom = AliPHOSGetter::GetInstance()->PHOSGeometry() ; 
  Int_t nevents = (Int_t) gAlice->TreeE()->GetEntries() ;
  Int_t ievent ;

  for(ievent = 0; ievent < nevents; ievent++){
    if(!ReadRecPoints(ievent))  //reads RecPoints for event ievent
      return;
    
    for(fModule = 1; fModule <= geom->GetNModules() ; fModule++ ){
      
      FillOneModule() ; 
      
      MakeLinks() ;
      
      MakePairs() ;
      
    }

    WriteTrackSegments(ievent) ;

    if(strstr(option,"deb"))
      PrintTrackSegments(option) ;

  }

  if(strstr(option,"tim")){
    gBenchmark->Stop("PHOSTSMaker");
    cout << "AliPHOSTSMaker:" << endl ;
    cout << "  took " << gBenchmark->GetCpuTime("PHOSTSMaker") << " seconds for making TS " 
	 <<  gBenchmark->GetCpuTime("PHOSTSMaker")/nevents << " seconds per event " << endl ;
    cout << endl ;
  }
    
}

//____________________________________________________________________________
void AliPHOSTrackSegmentMakerv1::Print(Option_t * option)const
{
  //  Print TrackSegmentMaker parameters

  if( strcmp(GetName(), "") != 0 ) {
    cout <<  "======== AliPHOSTrackSegmentMakerv1 ========" << endl ;
    cout <<  "Making Track segments "<< endl ;
    cout <<  "    Headers file:                   " << fHeaderFileName.Data() << endl ;
    cout <<  "    RecPoints branch file name:     " << fRecPointsBranchTitle.Data() << endl ;
    cout <<  "    TrackSegments Branch file name: " << fTrackSegmentsBranchTitle.Data() << endl ;
    cout <<  "with parameters: " << endl ;
    cout <<  "    Maximal EMC - CPV (PPSD) distance (cm)" << fR0 << endl ;
    cout <<  "============================================" << endl ;
  }
  else
    cout << "AliPHOSTrackSegmentMakerv1 not initialized " << endl ;
}
//____________________________________________________________________________
Bool_t AliPHOSTrackSegmentMakerv1::ReadRecPoints(Int_t event)
{
  // Reads Emc and CPV recPoints 
  // made previously with Clusterizer.


  //Make some initializations 

  fNTrackSegments = 0 ;
  fEmcFirst = 0 ;    
  fEmcLast  = 0 ;   
  fCpvFirst = 0 ;   
  fCpvLast  = 0 ;   
  fPpsdFirst= 0 ;   
  fPpsdLast = 0 ;   


  gAlice->GetEvent(event) ;
  gAlice->SetEvent(event) ;

  // Get TreeR header from file
  if(gAlice->TreeR()==0){
    cerr << "ERROR: AliPHOSTrackSegmentMakerv1::ReadRecPoints -> There is no Reconstruction Tree" << endl;
    return kFALSE;
  }


  // Find RecPoints
  TBranch * emcbranch = 0;
  TBranch * cpvbranch = 0;
  TBranch * clusterizerbranch = 0;
  TObjArray * lob = (TObjArray*)gAlice->TreeR()->GetListOfBranches() ;
  TIter next(lob) ; 
  TBranch * branch = 0 ;  
  Bool_t phosemcfound = kFALSE, phoscpvfound = kFALSE, clusterizerfound = kFALSE ; 
  
  TString taskName(GetName()) ; 
  taskName.ReplaceAll(Version(), "") ;

  while ( (branch = (TBranch*)next()) && (!phosemcfound || !phoscpvfound || !clusterizerfound) ) {
    if ( (strcmp(branch->GetName(), "PHOSEmcRP")==0) && (strcmp(branch->GetTitle(), taskName.Data())==0) ) {
      phosemcfound = kTRUE ;
      emcbranch = branch ; 
    }
    
    else if ( (strcmp(branch->GetName(), "PHOSCpvRP")==0) && (strcmp(branch->GetTitle(), taskName.Data())==0) ) {
      phoscpvfound = kTRUE ;
      cpvbranch = branch ; 
      
    } else if ( (strcmp(branch->GetName(), "AliPHOSClusterizer")==0) && (strcmp(branch->GetTitle(), taskName.Data())==0) ) {
      clusterizerfound = kTRUE ; 
      clusterizerbranch = branch ;
    }
  }
  if ( !phosemcfound || !phoscpvfound || !clusterizerfound ) {
    cerr << "WARNING: AliPHOSTrackSegmentMakerv1::ReadRecPoints -> emc(cpv)RecPoints and/or Clusterizer branch with name " << taskName.Data() 
	 << " not found" << endl ;
    return kFALSE ; 
  }   

  AliPHOSGetter * gime = AliPHOSGetter::GetInstance() ; 
  
  TObjArray * emcRecPoints = gime->EmcRecPoints() ;
  emcRecPoints->Clear() ; 
  emcbranch->SetAddress(&emcRecPoints) ;

  TObjArray * cpvRecPoints = gime->CpvRecPoints() ;
  cpvRecPoints->Clear() ; 
  cpvbranch->SetAddress(&cpvRecPoints) ;
 

  AliPHOSClusterizer * clusterizer = 0 ; 
  clusterizerbranch->SetAddress(&clusterizer) ;
  clusterizerbranch->GetEntry(0) ;
  TString clusterizerName( fTrackSegmentsBranchTitle ) ; 
  clusterizerName.Append(clusterizer->Version()) ; 
  clusterizer = gime->Clusterizer(clusterizerName) ; 

  emcbranch->GetEntry(0) ;
  cpvbranch->GetEntry(0) ;
  clusterizerbranch->GetEntry(0) ;
  
  return kTRUE ;
  
}

//____________________________________________________________________________
void AliPHOSTrackSegmentMakerv1::WriteTrackSegments(Int_t event)
{
  // Writes found TrackSegments to TreeR. Creates branches 
  // "PHOSTS" and "AliPHOSTrackSegmentMaker" with the same title.
  // In the former branch found TrackSegments are stored, while 
  // in the latter all parameters, with which TS were made. 
  // ROOT does not allow overwriting existing branches, therefore
  // first we check, if branches with the same title already exist.
  // If yes - exits without writing.
  
  AliPHOSGetter *gime = AliPHOSGetter::GetInstance() ; 
  TClonesArray * trackSegments = gime->TrackSegments() ; 

  gAlice->GetEvent(event) ; 

  //Make branch in TreeR for TrackSegments 
  char * filename = 0;
  if(gSystem->Getenv("CONFIG_SPLIT_FILE")!=0){   //generating file name
    filename = new char[strlen(gAlice->GetBaseFile())+20] ;
    sprintf(filename,"%s/PHOS.Reco.root",gAlice->GetBaseFile()) ; 
  }

  TDirectory *cwd = gDirectory;
  
  //First TS
  Int_t bufferSize = 32000 ;    
  TBranch * tsBranch = gAlice->TreeR()->Branch("PHOSTS",&trackSegments,bufferSize);
  tsBranch->SetTitle(fTrackSegmentsBranchTitle.Data());
  if (filename) {
    tsBranch->SetFile(filename);
    TIter next( tsBranch->GetListOfBranches());
    TBranch * sb ;
    while ((sb=(TBranch*)next())) {
      sb->SetFile(filename);
    }   
    cwd->cd();
  } 
  
  //Second -TSMaker
  Int_t splitlevel = 0 ;
  AliPHOSTrackSegmentMakerv1 * ts = this ;
  TBranch * tsMakerBranch = gAlice->TreeR()->Branch("AliPHOSTrackSegmentMaker","AliPHOSTrackSegmentMakerv1",
					  &ts,bufferSize,splitlevel);
  tsMakerBranch->SetTitle(fTrackSegmentsBranchTitle.Data());
  if (filename) {
    tsMakerBranch->SetFile(filename);
    TIter next( tsMakerBranch->GetListOfBranches());
    TBranch * sb;
    while ((sb=(TBranch*)next())) {
      sb->SetFile(filename);
    }   
    cwd->cd();
  } 
  
  tsBranch->Fill() ;  
  tsMakerBranch->Fill() ;

  gAlice->TreeR()->Write(0,kOverwrite) ;  
  
}


//____________________________________________________________________________
void AliPHOSTrackSegmentMakerv1::PrintTrackSegments(Option_t * option)
{
  // option deb - prints # of found TrackSegments
  // option deb all - prints as well indexed of found RecParticles assigned to the TS

  TClonesArray * trackSegments = AliPHOSGetter::GetInstance()->TrackSegments() ; 
  
  cout << "AliPHOSTrackSegmentMakerv1: " << endl ;
  cout << "       Found " << trackSegments->GetEntriesFast() << "  trackSegments " << endl ;
  
  if(strstr(option,"all")) {  // printing found TS
    cout << "TrackSegment # " << "    EMC RP#    " << "    CPV RP#    " << "     PPSD RP#" << endl ; 
    
    Int_t index;
    for (index = 0 ; index <trackSegments->GetEntriesFast() ; index++) {
      AliPHOSTrackSegment * ts = (AliPHOSTrackSegment * )trackSegments->At(index) ; 
      cout<<"   "<< setw(4) << ts->GetIndexInList() << "            " 
	  <<setw(4) << ts->GetEmcIndex()<< "            " 
	  <<setw(4) << ts->GetCpvIndex()<< "            " 
	  <<setw(4) << ts->GetPpsdIndex()<< endl ;
    }	
    
    cout << "-------------------------------------------------------"<< endl ;
  }
}
