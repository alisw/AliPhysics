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
#include "TROOT.h"
#include "TFile.h"
#include "TFolder.h"
#include "TTree.h"
#include "TSystem.h"
#include "TBenchmark.h"

// --- Standard library ---

// --- AliRoot header files ---

#include "AliEMCALTrackSegmentMakerv1.h"
#include "AliEMCALClusterizerv1.h"
#include "AliEMCALTrackSegment.h"
#include "AliEMCALLink.h"
#include "AliEMCALGetter.h"
#include "AliEMCAL.h"
#include "AliRun.h"

ClassImp( AliEMCALTrackSegmentMakerv1) 


//____________________________________________________________________________
  AliEMCALTrackSegmentMakerv1::AliEMCALTrackSegmentMakerv1() : AliEMCALTrackSegmentMaker()
{
  // default ctor (to be used mainly by Streamer)

  InitParameters() ; 

  fTrackSegmentsInRun       = 0 ; 

  fDefaultInit = kTRUE ; 
}

//____________________________________________________________________________
 AliEMCALTrackSegmentMakerv1::AliEMCALTrackSegmentMakerv1(const char * headerFile, const char * name, const Bool_t toSplit) : AliEMCALTrackSegmentMaker(headerFile, name, toSplit)
{
  // ctor

  InitParameters() ; 

  fTrackSegmentsInRun       = 0 ; 

  Init() ;

  fDefaultInit = kFALSE ; 

}

//____________________________________________________________________________
 AliEMCALTrackSegmentMakerv1::~AliEMCALTrackSegmentMakerv1()
{ 
  // dtor
  // fDefaultInit = kTRUE if TrackSegmentMaker created by default ctor (to get just the parameters)
  
  if (!fDefaultInit) {
    delete fPRELinkArray  ;
    delete fHCLinkArray  ;
    fSplitFile = 0 ; 
  }
}


//____________________________________________________________________________
const TString AliEMCALTrackSegmentMakerv1::BranchName() const 
{  
  TString branchName(GetName() ) ;
  branchName.Remove(branchName.Index(Version())-1) ;
  return branchName ;
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
    Info("HowClose", "ec = %d, rp = %d pro = %f, toofar=%d", ec->GetIndexInList(), rp->GetIndexInList(), pro, toofar ) ; 
 
  return r ;
}

//____________________________________________________________________________
void  AliEMCALTrackSegmentMakerv1::Init()
{
  // Make all memory allocations that are not possible in default constructor
  
  if ( strcmp(GetTitle(), "") == 0 )
    SetTitle("galice.root") ;
    
  TString branchname = GetName() ;
  branchname.Remove(branchname.Index(Version())-1) ;
  AliEMCALGetter * gime = AliEMCALGetter::GetInstance(GetTitle(),branchname.Data(), fToSplit ) ; 
  if ( gime == 0 ) {
    Error("Init", "Could not obtain the Getter object !") ; 
    return ;
  } 
  
  fSplitFile = 0 ;
  if(fToSplit){
    //First - extract full path if necessary
    TString fileName(GetTitle()) ;
    Ssiz_t islash = fileName.Last('/') ;
    if(islash<fileName.Length())
      fileName.Remove(islash+1,fileName.Length()) ;
    else
      fileName="" ;
    fileName+="EMCAL.RecData." ;
    if((strcmp(branchname.Data(),"Default")!=0)&&(strcmp(branchname.Data(),"")!=0)){
      fileName+=branchname ;
      fileName+="." ;
    }
    fileName+="root" ;
    fSplitFile = static_cast<TFile*>(gROOT->GetFile(fileName.Data()));   
    if(!fSplitFile)
      fSplitFile =  TFile::Open(fileName.Data(),"update") ;
  }
  
  fPRELinkArray = new TClonesArray("AliEMCALLink", 1000); 
  fHCLinkArray  = new TClonesArray("AliEMCALLink", 1000); 


  gime->PostTrackSegmentMaker(this) ;
  gime->PostTrackSegments(BranchName()) ; 

}

//____________________________________________________________________________
void  AliEMCALTrackSegmentMakerv1::InitParameters()
{
  fClose     = 10e-3 ;   
  fPRELinkArray = 0 ;
  fHCLinkArray  = 0 ;
  TString tsmName( GetName()) ; 
  if (tsmName.IsNull() ) 
    tsmName = "Default" ; 
  tsmName.Append(":") ; 
  tsmName.Append(Version()) ; 
  SetName(tsmName) ;
}


//____________________________________________________________________________
void  AliEMCALTrackSegmentMakerv1::MakeLinks()const
{ 
  // Finds distances (links) between all PRE, EC and HC clusters, 
  // which are not further apart from each other than fDangle 
  // and sort them in accordance with this distance
  
  AliEMCALGetter * gime = AliEMCALGetter::GetInstance() ; 
  TObjArray * aECRecPoints  = gime->ECALRecPoints() ; 
  TObjArray * aPRERecPoints = gime->PRERecPoints() ; 
  TObjArray * aHCRecPoints  = gime->HCALRecPoints() ; 

  fPRELinkArray->Clear() ;    
  fHCLinkArray->Clear() ;    

  AliEMCALTowerRecPoint * pre ;
  AliEMCALTowerRecPoint * ec ;
  AliEMCALTowerRecPoint * hc ;

  Int_t iPRELink  = 0 ;
  Int_t iHCLink   = 0 ;
    
  Int_t iECRP;
  for(iECRP = 0; iECRP < aECRecPoints->GetEntriesFast(); iECRP++ ) {
    ec = dynamic_cast<AliEMCALTowerRecPoint *>(aECRecPoints->At(iECRP)) ;
    Bool_t toofar = kTRUE ;        
    Int_t iPRERP = 0 ;    
    for(iPRERP = 0; iPRERP < aPRERecPoints->GetEntriesFast(); iPRERP++ ) { 
      pre = dynamic_cast<AliEMCALTowerRecPoint *>(aPRERecPoints->At(iPRERP)) ;
      Float_t prod = HowClose(ec, pre, toofar) ;    
      if(toofar)
	break ;	 
      if(prod < fClose) { 
	new ((*fPRELinkArray)[iPRELink++])  AliEMCALLink(prod, iECRP, iPRERP, 0) ;
      }      
    }
    toofar = kTRUE ; 
    Int_t iHCRP = 0 ;    
    for(iHCRP = 0; iHCRP < aHCRecPoints->GetEntriesFast(); iHCRP++ ) { 
      hc = dynamic_cast<AliEMCALTowerRecPoint *>(aHCRecPoints->At(iHCRP)) ;
      Float_t prod = HowClose(ec, hc, toofar) ;    
      if(toofar)
	break ;	 
      if(prod < fClose) { 
	new ((*fHCLinkArray)[iHCLink++])  AliEMCALLink(prod, iECRP, iHCRP, 1) ;
      }      
    }
  } 
  
  fPRELinkArray->Sort() ;  //first links with largest scalar product
  fHCLinkArray->Sort() ;  //first links with largest scalar product
}

//____________________________________________________________________________
void  AliEMCALTrackSegmentMakerv1::MakePairs()
{ 
  // Using the previously made list of "links", we found the best link - i.e. 
  // link with the largest scalar product (closest to one) to still 
  // unassigned RecParticles. We assign these RecPoints to TrackSegment and 
  // remove them from the list of "unassigned". 
  
  AliEMCALGetter * gime = AliEMCALGetter::GetInstance() ; 
  TObjArray * aECALRecPoints = gime->ECALRecPoints() ; 
  TObjArray * aPRERecPoints = gime->PRERecPoints() ; 
  TObjArray * aHCALRecPoints = gime->HCALRecPoints() ; 
  TClonesArray * trackSegments = gime->TrackSegments(BranchName()) ;   
    
  //Make arrays to mark clusters already chosen
  Int_t * ecalExist = 0;
  Int_t nECAL = aECALRecPoints->GetEntriesFast() ;  
  if (nECAL) 
    ecalExist = new Int_t[nECAL] ;
  
  Int_t index;
  for(index = 0; index < nECAL; index ++)
    ecalExist[index] = 1 ;
  
  Bool_t * preExist = 0;
  Int_t nPRE = aPRERecPoints->GetEntriesFast() ;  
  if(nPRE)
    preExist = new Bool_t[nPRE] ;
  for(index = 0; index < nPRE; index ++)
    preExist[index] = kTRUE ;
  
  Bool_t * hcalExist = 0;
  Int_t nHCAL = aHCALRecPoints->GetEntriesFast() ;  
  if(nHCAL)
    hcalExist = new Bool_t[nHCAL] ;
  for(index = 0; index < nHCAL; index ++)
    hcalExist[index] = kTRUE ;
  
  AliEMCALTowerRecPoint * null = 0 ; 
 // Finds the smallest links and makes pairs of PRE and ECAL clusters with largest scalar product 
 
  TIter nextPRE(fPRELinkArray) ;
  AliEMCALLink * linkPRE ;
  
  while ( (linkPRE =  static_cast<AliEMCALLink *>(nextPRE()) ) ){  

    if(ecalExist[linkPRE->GetECAL()] != -1){ //without PRE yet 

      if(preExist[linkPRE->GetOther()]){ // PRE still exist
	
	new ((* trackSegments)[fNTrackSegments]) 
	  AliEMCALTrackSegment(dynamic_cast<AliEMCALTowerRecPoint *>(aECALRecPoints->At(linkPRE->GetECAL())) , 
			       dynamic_cast<AliEMCALTowerRecPoint *>(aPRERecPoints->At(linkPRE->GetOther())), null) ;
	(dynamic_cast<AliEMCALTrackSegment *>(trackSegments->At(fNTrackSegments)))->SetIndexInList(fNTrackSegments);
	fNTrackSegments++ ;
	if (gDebug == 2 ) 
	  Info("MakePairs", "ECAL section with PRE section") ; 	
	ecalExist[linkPRE->GetECAL()] = -1 ; //Mark ecal  that pre was found 
	//mark PRE recpoint as already used 
	preExist[linkPRE->GetOther()] = kFALSE ;
      } //if PRE still exist
    } 
  }

  // Finds the smallest links and makes pairs of HCAL and ECAL clusters with largest scalar product 
 
  TIter nextHCAL(fHCLinkArray) ;
  AliEMCALLink * linkHCAL ;
  
  while ( (linkHCAL =  static_cast<AliEMCALLink *>(nextHCAL()) ) ){  

    if(ecalExist[linkHCAL->GetECAL()] != -2){ //without HCAL yet 

      if(hcalExist[linkHCAL->GetOther()]){ // HCAL still exist
	// search among the already existing track segments
	Int_t ii ; 
	Bool_t found = kFALSE ; 
	AliEMCALTrackSegment * ts = 0 ; 
	for ( ii = 0 ; ii < fNTrackSegments ; ii++ ) {
	  ts = dynamic_cast<AliEMCALTrackSegment *>(trackSegments->At(ii)) ;
	  if ( ts->GetECIndex() == linkHCAL->GetECAL() ) {
	    found = kTRUE ; 
	    break ; 
	  }
	}
	if (found){
	  ts->SetHCRecPoint( dynamic_cast<AliEMCALTowerRecPoint *>(aHCALRecPoints->At(linkHCAL->GetOther())) ) ;
	  if (gDebug == 2 ) 
	    Info("MakePairs", "ECAL section with PRE and HCAL sections") ;
	} 	
	if (!found) {
	  new ((* trackSegments)[fNTrackSegments]) 
	    AliEMCALTrackSegment(dynamic_cast<AliEMCALTowerRecPoint *>(aECALRecPoints->At(linkHCAL->GetECAL())), null,
				 dynamic_cast<AliEMCALTowerRecPoint *>(aHCALRecPoints->At(linkHCAL->GetOther()))) ; 
	(dynamic_cast<AliEMCALTrackSegment *>(trackSegments->At(fNTrackSegments)))->SetIndexInList(fNTrackSegments);
	fNTrackSegments++ ;
	if (gDebug == 2 ) 
	  Info("MakePairs", "ECAL section with HCAL section") ; 	
	}
	ecalExist[linkHCAL->GetECAL()] = -2 ; //Mark ecal  that hcal was found 
	//mark HCAL recpoint as already used 
	hcalExist[linkHCAL->GetOther()] = kFALSE ;
      } //if HCAL still exist
    } 
  }
  

  //look through ECAL recPoints left without PRE/HCAL
  if(ecalExist){ //if there is ecal rec point
    Int_t iECALRP ;
    for(iECALRP = 0; iECALRP < nECAL  ; iECALRP++ ){
      if(ecalExist[iECALRP] > 0 ){
	new ((*trackSegments)[fNTrackSegments])  
	  AliEMCALTrackSegment(dynamic_cast<AliEMCALTowerRecPoint *>(aECALRecPoints->At(iECALRP)), null, null) ;
	(dynamic_cast<AliEMCALTrackSegment *>(trackSegments->At(fNTrackSegments)))->SetIndexInList(fNTrackSegments);
	fNTrackSegments++;    
	if( gDebug == 2 ) 
	  Info("MakePairs", "ECAL section alone") ; 
     } 
    }
  }
  delete [] ecalExist ; 
  delete [] preExist ; 
  delete [] hcalExist ; 
}

//____________________________________________________________________________
void  AliEMCALTrackSegmentMakerv1::Exec(Option_t * option)
{
  // STEERing method

  if( strcmp(GetName(), "")== 0 ) 
    Init() ;

  if(strstr(option,"tim"))
    gBenchmark->Start("EMCALTSMaker");
 
  if(strstr(option,"print")) {
    Print("") ; 
    return ; 
  }

  AliEMCALGetter * gime = AliEMCALGetter::GetInstance() ; 
  if(gime->BranchExists("TrackSegments") )
    return ;

  Int_t nevents = gime->MaxEvent() ;       //(Int_t) gAlice->TreeE()->GetEntries() ;
  Int_t ievent ;

  for(ievent = 0; ievent < nevents; ievent++){

    gime->Event(ievent,"R") ;
 
    //Make some initializations 
    fNTrackSegments = 0 ;
    gime->TrackSegments(BranchName())->Clear() ; 
    
    MakeLinks() ;
    
    MakePairs() ;

    WriteTrackSegments(ievent) ;
    
    if(strstr(option,"deb"))
      PrintTrackSegments(option) ;
    
    //increment the total number of track segments per run 
    fTrackSegmentsInRun += gime->TrackSegments(BranchName())->GetEntriesFast() ; 
    
  }
  
  if(strstr(option,"tim")){
    gBenchmark->Stop("EMCALTSMaker");
    Info("Exec", "took %f seconds for making TS %f seconds per event", 
	 gBenchmark->GetCpuTime("EMCALTSMaker"), gBenchmark->GetCpuTime("EMCALTSMaker")/nevents) ;
  }
  
}

//____________________________________________________________________________
void AliEMCALTrackSegmentMakerv1::Print(Option_t * option)const
{
  //  Print TrackSegmentMaker parameters

  Info("Print", "TrackSegmentMakerv1 parameters:") ; 
  if( strcmp(GetName(), "") != 0 ) { 
    printf("Making Track segments with parameters:\n") ; 
    printf("    Allowed spred on the scalar product of two recpoints with same direction: %f\n", fClose) ;
    printf("============================================\n") ;
  }
  else
    printf("AliEMCALTrackSegmentMakerv1 not initialized ") ;
}

//____________________________________________________________________________
void AliEMCALTrackSegmentMakerv1::WriteTrackSegments(Int_t event)
{
  // Writes found TrackSegments to TreeR. Creates branches 
  // "EMCALTS" and "AliEMCALTrackSegmentMaker" with the same title.
  // In the former branch found TrackSegments are stored, while 
  // in the latter all parameters, with which TS were made. 
  // ROOT does not allow overwriting existing branches, therefore
  // first we check, if branches with the same title already exist.
  // If yes - exits without writing.
  
  AliEMCALGetter *gime = AliEMCALGetter::GetInstance() ; 

  TClonesArray * trackSegments = gime->TrackSegments() ; 
  trackSegments->Expand(trackSegments->GetEntriesFast()) ;
  TTree * treeR ;
  
  if(fToSplit){
    if(!fSplitFile)
      return ;
    fSplitFile->cd() ;
    char name[10] ;
    sprintf(name,"%s%d", "TreeR",event) ;
    treeR = dynamic_cast<TTree*>(fSplitFile->Get(name)); 
  }
  else{
    treeR = gAlice->TreeR();
  }
  
  if(!treeR){
    gAlice->MakeTree("R", fSplitFile);
    treeR = gAlice->TreeR() ;
  }
  
  //First TS
  Int_t bufferSize = 32000 ;    
  TBranch * tsBranch = treeR->Branch("EMCALTS",&trackSegments,bufferSize);
  tsBranch->SetTitle(BranchName());

  //Second -TSMaker
  Int_t splitlevel = 0 ;
  AliEMCALTrackSegmentMakerv1 * ts = this ;
  TBranch * tsMakerBranch = treeR->Branch("AliEMCALTrackSegmentMaker","AliEMCALTrackSegmentMakerv1",
					  &ts,bufferSize,splitlevel);
  tsMakerBranch->SetTitle(BranchName());

  tsBranch->Fill() ;  
  tsMakerBranch->Fill() ;

  treeR->AutoSave() ; //Write(0,kOverwrite) ;  
  if(gAlice->TreeR()!=treeR)
    treeR->Delete();
}


//____________________________________________________________________________
void AliEMCALTrackSegmentMakerv1::PrintTrackSegments(Option_t * option)
{
  // option deb - prints # of found TrackSegments
  // option deb all - prints as well indexed of found RecParticles assigned to the TS
  TString taskName(GetName()) ; 
  taskName.Remove(taskName.Index(Version())-1) ;
  
  TClonesArray * trackSegments = AliEMCALGetter::GetInstance()->TrackSegments(taskName) ; 


  Info("PrintTrackSegments", "Results from TrackSegmentMaker:") ; 
  printf("nevent: %d\n", gAlice->GetEvNumber()) ; 
  printf("        Found %d TrackSegments\n", trackSegments->GetEntriesFast() ); 

  if(strstr(option,"all")) {  // printing found TS
    printf("TrackSegment#  ECAL RP#  PRE RP#   HCAL RP#  \n") ; 
    Int_t index;
    for (index = 0 ; index < fNTrackSegments ; index++) {
      AliEMCALTrackSegment * ts = (AliEMCALTrackSegment * )trackSegments->At(index) ; 
      printf("   %d           %d        %d         %d \n", 
	     ts->GetIndexInList(), ts->GetECIndex(), ts->GetPREIndex(), ts->GetHCIndex() ); 
    }	
  }
}
