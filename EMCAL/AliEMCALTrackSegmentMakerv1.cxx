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
// Implementation version 1 of algorithm class to construct EMCAL track segments
// Track segment for EMCAL is list of 
//        ECAL RecPoint + (possibly) PRE RecPoint + (possibly) HCAL RecPoint
// To find TrackSegments we do the following: 
//  for each ECAL RecPoint we look for PRE and HC RecPoints with same direction  within fSame. 
//  If there is such a PRE or ECAL RecPoint, 
//   we make a "Link": indexes of ECAL and PRE, HCAL  RecPoints and their scalar product. 
//  Then we sort "Links", starting from the 
//   least "Link" pointing to the unassigned RecPoints assigning them to a new TrackSegment. 
//  If there is no PRE, HCAL RecPoint we make a TrackSegment 
//   consisting from ECAL alone. There is no TrackSegments without ECAL RecPoint.
//// In principle this class should be called from AliEMCALReconstructioner, but 
// one can use it as well in standalone mode.
//                 
//*-- Author: Dmitri Peressounko (RRC Ki & SUBATECH) & Yves Schutz (SUBATECH) 
//

// --- ROOT system ---
#include "TTree.h"
#include "TBenchmark.h"

// --- Standard library ---

// --- AliRoot header files ---

#include "AliEMCALTrackSegmentMakerv1.h"
#include "AliEMCALTrackSegment.h"
#include "AliEMCALLink.h"
#include "AliEMCALGetter.h"

ClassImp( AliEMCALTrackSegmentMakerv1) 


//____________________________________________________________________________
  AliEMCALTrackSegmentMakerv1::AliEMCALTrackSegmentMakerv1() : AliEMCALTrackSegmentMaker()
{
  // default ctor (to be used mainly by Streamer)

  InitParameters() ; 
  fDefaultInit = kTRUE ; 
}

//____________________________________________________________________________
 AliEMCALTrackSegmentMakerv1::AliEMCALTrackSegmentMakerv1(const TString alirunFileName, const TString eventFolderName)
   :AliEMCALTrackSegmentMaker(alirunFileName, eventFolderName)
{
  // ctor

  InitParameters() ; 
  Init() ;
  fDefaultInit = kFALSE ; 

}

//____________________________________________________________________________
 AliEMCALTrackSegmentMakerv1::~AliEMCALTrackSegmentMakerv1()
{ 
  // dtor
  // fDefaultInit = kTRUE if TrackSegmentMaker created by default ctor (to get just the parameters)

}

//____________________________________________________________________________
const TString AliEMCALTrackSegmentMakerv1::BranchName() const 
{  
   return GetName() ;

}

//____________________________________________________________________________
Float_t  AliEMCALTrackSegmentMakerv1::HowClose(AliEMCALTowerRecPoint * ec, AliEMCALTowerRecPoint * rp, Bool_t &toofar)const
{
  // Calculates the distance between the EMC RecPoint and the PPSD RecPoint
  // Clusters are sorted in "rows" and "columns" of width 1 cm

  Float_t r = -1. ;
  Float_t delta = 10. ;  // large enough to be ineffective ??! 

 
  TVector3 vecEC = ec->XYZInAlice() ;
  TVector3 vecRP = rp->XYZInAlice() ;
  
  Float_t pro = TMath::Abs(1 - (vecEC * vecRP / ( vecEC.Mag() * vecRP.Mag() ))) ; 

  if(pro <= delta) { 
    r = pro ;
    toofar = kFALSE ;
  } 
  else 
    toofar = kTRUE ;

  if (gDebug == 2 ) 
    printf("HowClose: ec = %d, rp = %d pro = %f, toofar=%d", ec->GetIndexInList(), rp->GetIndexInList(), pro, toofar ) ; 
 
  return r ;
}

//____________________________________________________________________________
void  AliEMCALTrackSegmentMakerv1::Init()
{
  // Make all memory allocations that are not possible in default constructor
  
  AliEMCALGetter* gime = AliEMCALGetter::Instance(GetTitle(), fEventFolderName.Data());
 
  if ( !gime->TrackSegmentMaker() ) {
    gime->PostTrackSegmentMaker(this);
  }
}

//____________________________________________________________________________
void  AliEMCALTrackSegmentMakerv1::InitParameters()
{
  fClose              = 10e-3 ;   
  fTrackSegmentsInRun = 0 ; 
}


//____________________________________________________________________________
void  AliEMCALTrackSegmentMakerv1::MakeLinks()const
{ 
  // Finds distances (links) between all PRE, EC and HC clusters, 
  // which are not further apart from each other than fDangle 
  // and sort them in accordance with this distance
  
  /* AliEMCALGetter * gime = AliEMCALGetter::Instance() ; 
  TObjArray * aECARecPoints  = gime->ECARecPoints() ; 
  // TObjArray * aPRERecPoints  = gime->PRERecPoints() ; 
  //TObjArray * aHCARecPoints  = gime->HCARecPoints() ; 

  fPRELinkArray->Clear() ;    
  fHCALinkArray->Clear() ;    

  AliEMCALTowerRecPoint * pre ;
  AliEMCALTowerRecPoint * eca ;
  AliEMCALTowerRecPoint * hca ;

  Int_t iPRELink  = 0 ;
  Int_t iHCALink  = 0 ;
    
  Int_t iECARP;
  for(iECARP = 0; iECARP < aECARecPoints->GetEntriesFast(); iECARP++ ) {
    eca = dynamic_cast<AliEMCALTowerRecPoint *>(aECARecPoints->At(iECARP)) ;
    Bool_t toofar = kTRUE ;        
    Int_t iPRERP = 0 ;    
    for(iPRERP = 0; iPRERP < aPRERecPoints->GetEntriesFast(); iPRERP++ ) { 
      pre = dynamic_cast<AliEMCALTowerRecPoint *>(aPRERecPoints->At(iPRERP)) ;
      Float_t prod = HowClose(eca, pre, toofar) ;    
      if(toofar)
	break ;	 
      if(prod < fClose) { 
	new ((*fPRELinkArray)[iPRELink++])  AliEMCALLink(prod, iECARP, iPRERP, 0) ;
      }      
    }
    toofar = kTRUE ; 
    Int_t iHCARP = 0 ;    
    for(iHCARP = 0; iHCARP < aHCARecPoints->GetEntriesFast(); iHCARP++ ) { 
      hca = dynamic_cast<AliEMCALTowerRecPoint *>(aHCARecPoints->At(iHCARP)) ;
      Float_t prod = HowClose(eca, hca, toofar) ;    
      if(toofar)
	break ;	 
      if(prod < fClose) { 
	new ((*fHCALinkArray)[iHCALink++])  AliEMCALLink(prod, iECARP, iHCARP, 1) ;
      }      
    }
  } 
  
  fPRELinkArray->Sort() ;  //first links with largest scalar product
  fHCALinkArray->Sort() ;  //first links with largest scalar product
  */
}

//____________________________________________________________________________
void  AliEMCALTrackSegmentMakerv1::MakePairs()
{ 
  // Using the previously made list of "links", we found the best link - i.e. 
  // link with the largest scalar product (closest to one) to still 
  // unassigned RecParticles. We assign these RecPoints to TrackSegment and 
  // remove them from the list of "unassigned". 
  
  /*AliEMCALGetter * gime = AliEMCALGetter::Instance() ; 
   TObjArray * aECARecPoints = gime->ECARecPoints() ; 
  TObjArray * aPRERecPoints = gime->PRERecPoints() ; 
  TObjArray * aHCARecPoints = gime->HCARecPoints() ; 
  TClonesArray * trackSegments = gime->TrackSegments() ;   
    
  //Make arrays to mark clusters already chosen
  Int_t * ecaExist = 0;
  Int_t nECA = aECARecPoints->GetEntriesFast() ;  
  if (nECA) 
    ecaExist = new Int_t[nECA] ;
  
  Int_t index;
  for(index = 0; index < nECA; index ++)
    ecaExist[index] = 1 ;
  
  Bool_t * preExist = 0;
  Int_t nPRE = aPRERecPoints->GetEntriesFast() ;  
  if(nPRE)
    preExist = new Bool_t[nPRE] ;
  for(index = 0; index < nPRE; index ++)
    preExist[index] = kTRUE ;
  
  Bool_t * hcaExist = 0;
  Int_t nHCA = aHCARecPoints->GetEntriesFast() ;  
  if(nHCA)
    hcaExist = new Bool_t[nHCA] ;
  for(index = 0; index < nHCA; index ++)
    hcaExist[index] = kTRUE ;
  
  AliEMCALTowerRecPoint * null = 0 ; 
 // Finds the smallest links and makes pairs of PRE and ECAL clusters with largest scalar product 
 
  TIter nextPRE(fPRELinkArray) ;
  AliEMCALLink * linkPRE ;
  
  while ( (linkPRE =  static_cast<AliEMCALLink *>(nextPRE()) ) ){  

    if(ecaExist[linkPRE->GetECA()] != -1){ //without PRE yet 

      if(preExist[linkPRE->GetOther()]){ // PRE still exist
	
	new ((* trackSegments)[fNTrackSegments]) 
	  AliEMCALTrackSegment(dynamic_cast<AliEMCALTowerRecPoint *>(aECARecPoints->At(linkPRE->GetECA())) , 
			       dynamic_cast<AliEMCALTowerRecPoint *>(aPRERecPoints->At(linkPRE->GetOther())), null) ;
	(dynamic_cast<AliEMCALTrackSegment *>(trackSegments->At(fNTrackSegments)))->SetIndexInList(fNTrackSegments);
	fNTrackSegments++ ;
	if (gDebug == 2 ) 
	  printf("MakePairs: ECAL section with PRE section") ; 	
	ecaExist[linkPRE->GetECA()] = -1 ; //Mark ecal  that pre was found 
	//mark PRE recpoint as already used 
	preExist[linkPRE->GetOther()] = kFALSE ;
      } //if PRE still exist
    } 
  }

  // Finds the smallest links and makes pairs of HCAL and ECAL clusters with largest scalar product 
 
  TIter nextHCA(fHCALinkArray) ;
  AliEMCALLink * linkHCA ;
  
  while ( (linkHCA =  static_cast<AliEMCALLink *>(nextHCA()) ) ){  

    if(ecaExist[linkHCA->GetECA()] != -2){ //without HCAL yet 

      if(hcaExist[linkHCA->GetOther()]){ // HCAL still exist
	// search among the already existing track segments
	Int_t ii ; 
	Bool_t found = kFALSE ; 
	AliEMCALTrackSegment * ts = 0 ; 
	for ( ii = 0 ; ii < fNTrackSegments ; ii++ ) {
	  ts = dynamic_cast<AliEMCALTrackSegment *>(trackSegments->At(ii)) ;
	  if ( ts->GetECAIndex() == linkHCA->GetECA() ) {
	    found = kTRUE ; 
	    break ; 
	  }
	}
	if (found){
	  ts->SetHCARecPoint( dynamic_cast<AliEMCALTowerRecPoint *>(aHCARecPoints->At(linkHCA->GetOther())) ) ;
	  if (gDebug == 2 ) 
	    printf("MakePairs: ECAL section with PRE and HCAL sections") ;
	} 	
	if (!found) {
	  new ((* trackSegments)[fNTrackSegments]) 
	    AliEMCALTrackSegment(dynamic_cast<AliEMCALTowerRecPoint *>(aECARecPoints->At(linkHCA->GetECA())), null,
				 dynamic_cast<AliEMCALTowerRecPoint *>(aHCARecPoints->At(linkHCA->GetOther()))) ; 
	(dynamic_cast<AliEMCALTrackSegment *>(trackSegments->At(fNTrackSegments)))->SetIndexInList(fNTrackSegments);
	fNTrackSegments++ ;
	if (gDebug == 2 ) 
	  printf("MakePairs: ECAL section with HCAL section") ; 	
	}
	ecaExist[linkHCA->GetECA()] = -2 ; //Mark ecal  that hcal was found 
	//mark HCAL recpoint as already used 
	hcaExist[linkHCA->GetOther()] = kFALSE ;
      } //if HCAL still exist
    } 
  }
  

  //look through ECAL recPoints left without PRE/HCAL
  if(ecaExist){ //if there is ecal rec point
    Int_t iECARP ;
    for(iECARP = 0; iECARP < nECA  ; iECARP++ ){
      if(ecaExist[iECARP] > 0 ){
	new ((*trackSegments)[fNTrackSegments])  
	  AliEMCALTrackSegment(dynamic_cast<AliEMCALTowerRecPoint *>(aECARecPoints->At(iECARP)), null, null) ;
	(dynamic_cast<AliEMCALTrackSegment *>(trackSegments->At(fNTrackSegments)))->SetIndexInList(fNTrackSegments);
	fNTrackSegments++;    
	if( gDebug == 2 ) 
	  printf("MakePairs: ECAL section alone") ; 
     } 
    }
  }
  delete [] ecaExist ; 
  delete [] preExist ; 
  delete [] hcaExist ; 
  */
}

//____________________________________________________________________________
void  AliEMCALTrackSegmentMakerv1::Exec(Option_t * option)
{
  // STEERing method


  if(strstr(option,"tim"))
    gBenchmark->Start("EMCALTSMaker");
 
  if(strstr(option,"print")) {
    Print("") ; 
    return ; 
  }

  AliEMCALGetter * gime = AliEMCALGetter::Instance() ; 

  Int_t nevents = gime->MaxEvent() ;   
  Int_t ievent ;

  for(ievent = 0; ievent < nevents; ievent++){
    gime->Event(ievent,"R") ;
    //Make some initializations 
    fNTrackSegments = 0 ;
    
    gime->TrackSegments()->Clear() ; 
    
    MakeLinks() ;
    MakePairs() ;

    WriteTrackSegments() ;
    
    if(strstr(option,"deb"))
      PrintTrackSegments(option) ;
    
    //increment the total number of track segments per run 
    fTrackSegmentsInRun += gime->TrackSegments()->GetEntriesFast() ; 
    
  }
  
  if(strstr(option,"tim")){
    gBenchmark->Stop("EMCALTSMaker");
    printf("Exec: took %f seconds for making TS %f seconds per event", 
	 gBenchmark->GetCpuTime("EMCALTSMaker"), gBenchmark->GetCpuTime("EMCALTSMaker")/nevents) ;
  }
   Unload();
}

//____________________________________________________________________________
void AliEMCALTrackSegmentMakerv1::Unload() 
{
  // Unloads the RecPoints and Tracks
  AliEMCALGetter * gime = AliEMCALGetter::Instance() ;  
  gime->EmcalLoader()->UnloadRecPoints() ;
  gime->EmcalLoader()->UnloadTracks() ;
}

//____________________________________________________________________________
void AliEMCALTrackSegmentMakerv1::Print(Option_t * /*option*/)const
{
  //  Print TrackSegmentMaker parameters

  printf("Print: TrackSegmentMakerv1 parameters:") ; 
  if( strcmp(GetName(), "") != 0 ) { 
    printf("Making Track segments with parameters:\n") ; 
    printf("    Allowed spred on the scalar product of two recpoints with same direction: %f\n", fClose) ;
    printf("============================================\n") ;
  }
  else
    printf("AliEMCALTrackSegmentMakerv1 not initialized ") ;
}

//____________________________________________________________________________
void AliEMCALTrackSegmentMakerv1::WriteTrackSegments()
{
  // Writes found TrackSegments to TreeR. Creates branches 
  // "EMCALTS" and "AliEMCALTrackSegmentMaker" with the same title.
  // In the former branch found TrackSegments are stored, while 
  // in the latter all parameters, with which TS were made. 
  // ROOT does not allow overwriting existing branches, therefore
  // first we check, if branches with the same title already exist.
  // If yes - exits without writing.
  
  AliEMCALGetter *gime = AliEMCALGetter::Instance() ; 

  TClonesArray * trackSegments = gime->TrackSegments() ; 
  trackSegments->Expand(trackSegments->GetEntriesFast()) ;

  TTree * treeT = gime->TreeT();
  
  //First TS
  Int_t bufferSize = 32000 ;    
  TBranch * tsBranch = treeT->Branch("EMCALTS",&trackSegments,bufferSize);
  tsBranch->Fill() ; 

  gime->WriteTracks("OVERWRITE");
  gime->WriteTrackSegmentMaker("OVERWRITE");
}


//____________________________________________________________________________
void AliEMCALTrackSegmentMakerv1::PrintTrackSegments(Option_t * option)
{
  // option deb - prints # of found TrackSegments
  // option deb all - prints as well indexed of found RecParticles assigned to the TS
  
  TClonesArray * trackSegments = AliEMCALGetter::Instance()->TrackSegments() ; 


  printf("PrintTrackSegments: Results from TrackSegmentMaker:") ; 
  printf("nevent: %d\n", gAlice->GetEvNumber()) ; 
  printf("        Found %d TrackSegments\n", trackSegments->GetEntriesFast() ); 

  if(strstr(option,"all")) {  // printing found TS
    printf("TrackSegment#  ECAL RP#  \n") ; 
    Int_t index;
    for (index = 0 ; index < fNTrackSegments ; index++) {
      AliEMCALTrackSegment * ts = (AliEMCALTrackSegment * )trackSegments->At(index) ; 
      printf("   %d           %d \n", 
	     ts->GetIndexInList(), ts->GetECAIndex()); 
    }	
  }
}
