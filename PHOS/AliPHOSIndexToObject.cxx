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

/* $Id:  */

/* $Log:
   29.05.2001 Yuri Kharlov:
              Everywhere reading the treese TTree->GetEvent(i)
              is replaced by reading the branches TBranch->GetEntry(0)
*/

//_________________________________________________________________________
//  A singleton. This class should be used in the analysis stage to get 
//  reconstructed objects: Digits, RecPoints, TrackSegments and RecParticles,
//  instead of directly reading them from galice.root file. This container 
//  ensures, that one reads Digits, made of these particular digits, RecPoints, 
//  made of these particular RecPoints, TrackSegments and RecParticles. 
//  This becomes non trivial if there are several identical branches, produced with
//  different set of parameters. 
//
//  An example of how to use (see also class AliPHOSAnalyser):
//  AliPHOSIndexToObject * please = AliPHOSIndexToObject::GetInstance("galice.root","RecParticles","") ;
//  for(Int_t irecp = 0; irecp < please->GimeNRecParticles() ; irecp++)
//     AliPHOSRecParticle * part = please->GimeRecParticle(1) ;
//     ................
//  please->GetEvent(event) ;    // reads new event from galice.root
//                  
//*-- Author: Yves Schutz (SUBATECH) & Dmitri Peressounko (RRC KI & SUBATECH)
//*--         Complitely redesigned by Dmitri Peressounko March 2001  
//////////////////////////////////////////////////////////////////////////////


// --- ROOT system ---

#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TObjString.h"

// --- Standard library ---
#include <iostream.h>

// --- AliRoot header files ---

#include "AliRun.h"
#include "AliPHOSIndexToObject.h"
#include "AliPHOSDigitizer.h"
#include "AliPHOSSDigitizer.h"
#include "AliPHOSClusterizer.h"
#include "AliPHOSClusterizerv1.h"
#include "AliPHOSTrackSegmentMaker.h"
#include "AliPHOSTrackSegmentMakerv1.h"
#include "AliPHOSTrackSegment.h"
#include "AliPHOSPID.h" 
#include "AliPHOSPIDv1.h" 

ClassImp(AliPHOSIndexToObject)
  
  AliPHOSIndexToObject * AliPHOSIndexToObject::fgObjGetter = 0 ; 

//____________________________________________________________________________ 
AliPHOSIndexToObject::AliPHOSIndexToObject(const char* headerFile,const char* branch,const char* branchTitle )
{
  //Initialize  all lists
  fEvent = 0 ;

  fSDigits = new TClonesArray("AliPHOSDigit",100) ;
  fDigits  = new TClonesArray("AliPHOSDigit",100) ;
  fEmcRecPoints = new TObjArray(100) ;
  fCpvRecPoints = new TObjArray(100) ;
  fTS = new TClonesArray("AliPHOSTrackSegment",100) ;
  fRecParticles = new TClonesArray("AliPHOSRecParticle",100) ;
  fPrimaries = new TObjArray(1) ;

  fSDigitizer = 0 ;
  fDigitizer = 0 ;
  fClusterizer = 0 ;
  fTSMaker = 0 ;
  fPID = 0 ;

  //open headers file
  fHeaderFile = headerFile ;
  TFile * file = (TFile*) gROOT->GetFile(fHeaderFile.Data() ) ;

  if(file == 0){
    if(fHeaderFile.Contains("rfio")) // if we read file using HPSS
      file =	TFile::Open(fHeaderFile.Data(),"update") ;
    else
      file = new TFile(fHeaderFile.Data(),"update") ;
    gAlice = (AliRun *) file->Get("gAlice") ;
  }

  fMaxEvent = (Int_t) gAlice->TreeE()->GetEntries() ;

  DefineBranchTitles(branch,branchTitle) ;

  //Now read all data from trees
  fEvent = -1 ;
  GetEvent(0) ;

}
//____________________________________________________________________________ 
void AliPHOSIndexToObject:: DefineBranchTitles(const char* startBranch,const char* branchTitle)
{
  // Points to the branches of all reconstructed objects with the specified names

  gAlice->GetEvent(0) ;
  // Read all reconstruction classes to extract titles of 
  // branches, constituing "reconstruction branch"
  AliPHOSPID * pids[50];  // here AliPHOSPID's will be stored
  Int_t ipids = 0 ;
  AliPHOSTrackSegmentMaker * tsms[50] ; 
  Int_t itsms = 0 ;
  AliPHOSClusterizer * clus[50] ;   
  Int_t iclus = 0 ;
  AliPHOSDigitizer * digs[50]; 
  Int_t idigs = 0 ;
  AliPHOSSDigitizer * sdigs[50];
  Int_t isdigs = 0 ;
  
  //read TreeR
  TObjArray * branches = gAlice->TreeR()->GetListOfBranches() ;
  Int_t ibranch;
  TBranch * branch ;
  for(ibranch = 0;ibranch <branches->GetEntries();ibranch++){
    branch=(TBranch *) branches->At(ibranch) ;
    if( (strcmp(branch->GetName(),"AliPHOSPID") == 0) ){
      pids[ipids] =  new  AliPHOSPIDv1() ;
      branch->SetAddress(& (pids[ipids])) ;
      branch->GetEntry(0) ;
      ipids++ ;
    }
    if( (strcmp(branch->GetName(),"AliPHOSTrackSegmentMaker") == 0) ){
      tsms[itsms] = new  AliPHOSTrackSegmentMakerv1() ;
      branch->SetAddress(&(tsms[itsms])) ;
      branch->GetEntry(0) ;
      itsms++ ;
    }
    if( (strcmp(branch->GetName(),"AliPHOSClusterizer") == 0) ){
      clus[iclus] = new  AliPHOSClusterizerv1() ;
      branch->SetAddress(&(clus[iclus])) ;
      branch->GetEntry(0) ;
      iclus++ ;
    }
  }
//    gAlice->TreeR()->GetEvent(0) ;


  //read TreeD
  branches = gAlice->TreeD()->GetListOfBranches() ;
  for(ibranch = 0;ibranch <branches->GetEntries();ibranch++){
    branch=(TBranch *) branches->At(ibranch) ;
    if( (strcmp(branch->GetName(),"AliPHOSDigitizer") == 0) ){
      digs[idigs] = new  AliPHOSDigitizer() ;
      branch->SetAddress(&(digs[idigs])) ;
      branch->GetEntry(0) ;
      idigs++ ;
    }
  }
//    gAlice->TreeD()->GetEvent(0) ;

  //read TreeS
  branches = gAlice->TreeS()->GetListOfBranches() ;
  for(ibranch = 0;ibranch <branches->GetEntries();ibranch++){
    branch=(TBranch *) branches->At(ibranch) ;
    gAlice->TreeS()->GetBranch(branch->GetName())->GetEntry(0) ; // YK
    if( (strcmp(branch->GetName(),"AliPHOSSDigitizer") == 0) ){
      sdigs[isdigs] = new  AliPHOSSDigitizer() ;
      branch->SetAddress(&(sdigs[isdigs])) ;
      branch->GetEntry(0) ;
      isdigs++ ;
    }
  }
//    gAlice->TreeS()->GetEvent(0) ;

  // now choose among read Reconstruction classes those,
  // which constituite "reconstruction branch"
  Bool_t pidDefined = kFALSE ;
  Bool_t tsmDefined = kFALSE ;
  Bool_t cluDefined = kFALSE ;
  Bool_t digDefined = kFALSE ;
  Bool_t sdigDefined = kFALSE ;

  Int_t index ;
  // First, go from the end (RecParticles) to the beginning(SDigits)
  if((strcmp(startBranch,"PHOSRP") == 0)||(strcmp(startBranch,"AliPHOSPID") == 0)){
    fRPTitle = branchTitle ;
    for(index = 0; index < ipids ; index++){
      if(fRPTitle.CompareTo(((AliPHOSPID*)pids[index])->GetRecParticlesBranch())== 0){
	pidDefined = kTRUE ;
	fTSTitle =((AliPHOSPID*)pids[index])->GetTrackSegmentsBranch() ; 
      }
    }
  }
  if((strcmp(startBranch,"PHOSTS") == 0)||(strcmp(startBranch,"AliPHOSTrackSegmentMaker") == 0)|| pidDefined ) {
    if(!pidDefined)
      fTSTitle = branchTitle ;
    for(index = 0; index < itsms ; index++)
      if(fTSTitle.CompareTo(((AliPHOSTrackSegmentMaker*)tsms[index])->GetTrackSegmentsBranch())== 0){
	tsmDefined = kTRUE ;
	fRecPointsTitle =((AliPHOSTrackSegmentMaker*)tsms[index])->GetRecPointsBranch() ;
      }
  }
  if((strcmp(startBranch,"PHOSEmcRP") == 0) || 
     (strcmp(startBranch,"PHOSCpvRP") == 0) ||
     (strcmp(startBranch,"AliPHOSClusterizer") == 0)  || tsmDefined ) {
    if(!tsmDefined)
      fRecPointsTitle = branchTitle ;
    for(index = 0; index < iclus ; index++)
      if(fRecPointsTitle.CompareTo(((AliPHOSClusterizer*)clus[index])->GetRecPointsBranch())== 0){
	cluDefined = kTRUE ;
	fDigitsTitle =((AliPHOSClusterizer*)clus[index])->GetDigitsBranch() ;
      }
  }
  if((strcmp(startBranch,"PHOS") == 0) || (strcmp(startBranch,"AliPHOSDigitizer") == 0) ||cluDefined ) {
    if(!cluDefined)
      fDigitsTitle = branchTitle ; 
    for(index = 0; index < idigs ; index++) {
      if(fDigitsTitle.CompareTo(((AliPHOSDigitizer*)digs[index])->GetDigitsBranch())== 0){
	digDefined = kTRUE ;
	fSDigitsTitle =((AliPHOSDigitizer*)digs[index])->GetSDigitsBranch() ;
      }
    }
  }
  for(index = 0; index < idigs ; index++)
    if(fSDigitsTitle.CompareTo(((AliPHOSSDigitizer*)sdigs[index])->GetSDigitsBranch())== 0)
      sdigDefined = kTRUE ;
  
  if(!sdigDefined){
    cout << "Can not define titles of branches " << endl ;
    cout << endl ;
  }

  // Now we go in the inverse direction: from sdigits to recparticles - for the 
  // case, if we started decending not from RecParticles, but e.g. from digits

  if( !cluDefined ) {
    for(index = 0; index < iclus ; index++)
      if(fDigitsTitle.CompareTo(((AliPHOSClusterizer*)clus[index])->GetDigitsBranch())== 0){
	cluDefined = kTRUE ;
	fRecPointsTitle =((AliPHOSClusterizer*)clus[index])->GetRecPointsBranch() ;
      }
  }
  if(! tsmDefined ) {
    for(index = 0; index < itsms ; index++)
      if(fRecPointsTitle.CompareTo(((AliPHOSTrackSegmentMaker*)tsms[index])->GetRecPointsBranch())== 0){
	tsmDefined = kTRUE ;
	fTSTitle =((AliPHOSTrackSegmentMaker*)tsms[index])->GetTrackSegmentsBranch() ;
      }
  }
  if(!pidDefined){
    for(index = 0; index < ipids ; index++)
      if(fTSTitle.CompareTo(((AliPHOSPID*)pids[index])->GetTrackSegmentsBranch())== 0){
	pidDefined = kTRUE ;
	fRPTitle = ((AliPHOSPID*)pids[index])->GetRecParticlesBranch() ;
      }
  }
  
  //delete created objects
  for(index = 0; index < ipids ; index++)
    delete pids[index] ;
  for(index = 0; index < itsms ; index++)
    delete tsms[index] ;
  for(index = 0; index < iclus ; index++)
    delete clus[index] ;
  for(index = 0; index < idigs ; index++) 
    delete digs[index] ;
  for(index = 0; index < isdigs ; index++) 
    delete sdigs[index] ;
  
}
//____________________________________________________________________________ 
AliPHOSIndexToObject * AliPHOSIndexToObject::GetInstance()
{
  // Returns the pointer of the unique instance already defined
  
  AliPHOSIndexToObject * rv = 0 ;
  if ( fgObjGetter )
    rv = fgObjGetter ;
  else
    cout << "AliPHOSIndexToObject::GetInstance ERROR: not yet initialized" << endl ;

  return rv ;
}

//____________________________________________________________________________ 
AliPHOSIndexToObject * AliPHOSIndexToObject::GetInstance(const char* headerFile,
							 const char* branch,
							 const char* branchTitle)
{
  // Creates and returns the pointer of the unique instance
  // Must be called only when the environment has changed 

  if(strcmp(branch,"PHOSRP") && strcmp(branch,"AliPHOSPID") &&
     strcmp(branch,"PHOSTS") && strcmp(branch,"AliPHOSTrackSegmentMaker") && 
     strcmp(branch,"PHOSEmcRP") && strcmp(branch,"PHOSCpvRP") && strcmp(branch,"AliPHOSClusterizer") &&
     strcmp(branch,"PHOS") && strcmp(branch,"AliPHOSDigitizer") ){
    
    cout << "AliPHOSIndexToObject: wrong branch name specified: " << branch << endl ;
    cout << "   avalilable names are `PHOSRP', `AliPHOSPID'"<<endl ;
    cout << "                        `PHOSTS', `AliPHOSTrackSegmentMaker'"<<endl ;
    cout << "                        `PHOSEmcRP', `PHOSCpvRP', `AliPHOSClusterizer'"<< endl ;
    cout << "                        `PHOS' and `AliPHOSDigitizer'"<< endl ;
    return 0 ;
  }


  if ( fgObjGetter )      // delete it if already exists
    delete fgObjGetter ; 

  fgObjGetter = new AliPHOSIndexToObject(headerFile,branch,branchTitle) ; 
  
  return fgObjGetter ; 
  
}

//____________________________________________________________________________ 
TParticle * AliPHOSIndexToObject::GimePrimary(Int_t index) const
{
  // Return primary particle numbered by <index>

  if(index < 0) 
    return 0 ;
  
  Int_t primaryIndex = index % 10000000 ; 
  Int_t primaryList = (Int_t ) ((index-primaryIndex)/10000000.)  ;
  
  if ( primaryList > 0  ) {
    cout << " IndexToObject does not support currently Mixing of primary " << endl ;
    cout << "   can not return primary: " << index<< " (list "<< primaryList<< " primary # " << primaryIndex << " )"<<endl ;
    return 0;
  }
  
  return gAlice->Particle(primaryIndex) ;
  
}

//____________________________________________________________________________ 
void AliPHOSIndexToObject::ReadTreeD()
{
  // Read the digit tree gAlice->TreeD()  
  if(gAlice->TreeD()== 0){
    cout << "AliPHOSIndexToObject : can not read TreeD " << endl;
    return ;
  }
  
  TBranch * digitsBranch = 0;
  TBranch * digitizerBranch = 0;
  
  TObjArray * branches = gAlice->TreeD()->GetListOfBranches() ;
  Int_t ibranch;
  Bool_t phosNotFound = kTRUE ;
  Bool_t digitizerNotFound = kTRUE ;
  
  for(ibranch = 0;ibranch <branches->GetEntries();ibranch++){
    
    if(phosNotFound){
      digitsBranch=(TBranch *) branches->At(ibranch) ;
      if( (strcmp(digitsBranch->GetTitle(),fDigitsTitle)==0 ) &&
	  (strcmp(digitsBranch->GetName(),"PHOS") == 0) )
	phosNotFound = kFALSE ;
    }
    if(digitizerNotFound){
      digitizerBranch = (TBranch *) branches->At(ibranch) ;
      if( (strcmp(digitizerBranch->GetTitle(),fDigitsTitle) == 0) && 
	  (strcmp(digitizerBranch->GetName(),"AliPHOSDigitizer") == 0) )
	digitizerNotFound = kFALSE ;
    } 
  }
    
  if(digitizerNotFound || phosNotFound){
    cout << "AliPHOSIndexToObject error: " << endl ;
    cout << "       Can't find Branch with Digits or Digitizer "<< endl ; ;
    return  ;
  }
  
  digitsBranch   ->SetAddress(&fDigits) ;
  digitizerBranch->SetAddress(&fDigitizer) ;
  digitsBranch   ->GetEntry(0) ;
  digitizerBranch->GetEntry(0) ;
  
//    gAlice->TreeD()->GetEvent(0) ;    // YK 29.05.2001
  
}
//____________________________________________________________________________ 
void AliPHOSIndexToObject::ReadTreeS()
{
  // Read the summable digits tree gAlice->TreeS()  

  if(gAlice->TreeS()== 0){
    cout <<   "AliPHOSIndexToObject: can not read TreeS " << endl ;
    return ;
  }
  
  TBranch * sdigitsBranch = 0;
  TBranch * sdigitizerBranch = 0;
  
  TObjArray * branches = gAlice->TreeS()->GetListOfBranches() ;
  Int_t ibranch;
  Bool_t phosNotFound = kTRUE ;
  Bool_t sdigitizerNotFound = kTRUE ;
  
  for(ibranch = 0;ibranch <branches->GetEntries();ibranch++){
    
    if(phosNotFound){
      sdigitsBranch=(TBranch *) branches->At(ibranch) ;
      if( (strcmp(sdigitsBranch->GetTitle(),fSDigitsTitle)==0 ) &&
	  (strcmp(sdigitsBranch->GetName(),"PHOS") == 0) )
	phosNotFound = kFALSE ;
    }
    if(sdigitizerNotFound){
      sdigitizerBranch = (TBranch *) branches->At(ibranch) ;
      if( (strcmp(sdigitizerBranch->GetTitle(),fSDigitsTitle) == 0) && 
	  (strcmp(sdigitizerBranch->GetName(),"AliPHOSSDigitizer") == 0) )
	sdigitizerNotFound = kFALSE ;
    } 
  }
  
  if(sdigitizerNotFound || phosNotFound){
    cout << "AliPHOSIndexToObject error: " << endl ;
    cout << "       Can't find Branch with SDigits or SDigitizer "<< endl ; ;
    return ;
  }
  
  sdigitsBranch   ->SetAddress(&fSDigits) ;
  sdigitizerBranch->SetAddress(&fSDigitizer) ;
  sdigitsBranch   ->GetEvent(0) ;
  sdigitizerBranch->GetEvent(0) ;
  
//    gAlice->TreeS()->GetEvent(0) ;    // YK 29.05.2001
  
}
//____________________________________________________________________________ 
void AliPHOSIndexToObject::ReadTreeR()
{
  // Read the reconstrunction tree gAlice->TreeR()

  if(gAlice->TreeR()== 0){
    cout <<   "AliPHOSIndexToObject: can not read TreeR " << endl ;
    return ;
  }

  TBranch * pidBranch = 0;
  TBranch * rpBranch = 0;
  TBranch * tsMakerBranch = 0;
  TBranch * tsBranch = 0;
  TBranch * emcBranch = 0;
  TBranch * cpvBranch = 0;
  TBranch * clusterizerBranch = 0;
  
  TObjArray * branches = gAlice->TreeR()->GetListOfBranches() ;
  Int_t ibranch;
  Bool_t pidNotFound = kTRUE ;
  Bool_t rpNotFound = kTRUE ;
  Bool_t tsMakerNotFound = kTRUE ;
  Bool_t tsNotFound = kTRUE ;
  Bool_t emcNotFound = kTRUE ;
  Bool_t cpvNotFound = kTRUE ;  
  Bool_t clusterizerNotFound = kTRUE ;
  
  for(ibranch = 0;ibranch <branches->GetEntries();ibranch++){
    
    if(pidNotFound){
      pidBranch=(TBranch *) branches->At(ibranch) ;
      if( (fRPTitle.CompareTo(pidBranch->GetTitle())==0 ) &&
	  (strcmp(pidBranch->GetName(),"AliPHOSPID") == 0) )
	pidNotFound = kFALSE ;
    }
    if(rpNotFound){
      rpBranch=(TBranch *) branches->At(ibranch) ;
      if( (fRPTitle.CompareTo(rpBranch->GetTitle())==0 ) &&
	  (strcmp(rpBranch->GetName(),"PHOSRP") == 0) )
	rpNotFound = kFALSE ;
    }
    if(tsMakerNotFound){
      tsMakerBranch=(TBranch *) branches->At(ibranch) ;
      if( fTSTitle.CompareTo(tsMakerBranch->GetTitle())==0 )
	if( strcmp(tsMakerBranch->GetName(),"AliPHOSTrackSegmentMaker") == 0) 
	  tsMakerNotFound = kFALSE ;
    }
    if(tsNotFound){
      tsBranch=(TBranch *) branches->At(ibranch) ;
      if( fTSTitle.CompareTo(tsBranch->GetTitle())==0 )
	if( strcmp(tsBranch->GetName(),"PHOSTS") == 0) 
	  tsNotFound = kFALSE ;
    }  
    if(emcNotFound){
      emcBranch=(TBranch *) branches->At(ibranch) ;
      if( (fRecPointsTitle.CompareTo(emcBranch->GetTitle()) == 0) && 
	  (strcmp(emcBranch->GetName(),"PHOSEmcRP") == 0) )
	emcNotFound = kFALSE ;
    }
    if(cpvNotFound){
      cpvBranch=(TBranch *) branches->At(ibranch) ;
      if( (fRecPointsTitle.CompareTo(cpvBranch->GetTitle()) == 0) &&
	  (strcmp(cpvBranch->GetName(),"PHOSCpvRP") == 0) )
	cpvNotFound = kFALSE ;
    }
    if(clusterizerNotFound){
      clusterizerBranch = (TBranch *) branches->At(ibranch) ;
      if( (fRecPointsTitle.CompareTo(clusterizerBranch->GetTitle()) == 0) &&
	  (strcmp(clusterizerBranch->GetName(),"AliPHOSClusterizer") == 0) )
	clusterizerNotFound = kFALSE ;
    }
  }

  if(pidNotFound ||rpNotFound ){
    cout << "AliPHOSIndexToObject error" << endl ;
    cout << "     Can't find Branch with PID and RecParticles " ;
    return  ;
  }
  if(tsMakerNotFound ||tsNotFound ){
    cout << "AliPHOSIndexToObject error" << endl ;
    cout << "       Can't find Branch with TrackSegmentMaker and TrackSegments " ;
    cout << "       Do nothing" <<endl  ;
    return ;
  }
  if(clusterizerNotFound || emcNotFound || cpvNotFound){
    cout << "AliPHOSIndexToObject error" << endl ;
    cout << "       Can't find Branch with RecPoints or Clusterizer " << endl ;
    return ;
  }

  //    YK 29.05.2001 : Read branch instead of tree
  emcBranch        ->SetAddress(&fEmcRecPoints) ;
  cpvBranch        ->SetAddress(&fCpvRecPoints) ;
  clusterizerBranch->SetAddress(&fClusterizer) ;
  emcBranch        ->GetEntry(0) ;
  cpvBranch        ->GetEntry(0) ;
  clusterizerBranch->GetEntry(0) ;
  
  tsMakerBranch    ->SetAddress(&fTSMaker) ;
  tsBranch         ->SetAddress(&fTS) ;
  tsMakerBranch    ->GetEntry(0) ;
  tsBranch         ->GetEntry(0) ;
    
  pidBranch        ->SetAddress(&fPID) ;
  rpBranch         ->SetAddress(&fRecParticles) ;
  pidBranch        ->GetEntry(0) ;
  rpBranch         ->GetEntry(0) ;
  
//    gAlice->TreeR()->GetEvent(0) ;    // YK 29.05.2001

}
//____________________________________________________________________________ 
void AliPHOSIndexToObject::ReadPrimaries()
{
  // Reads specific branches of primaries
  
  fNPrimaries = gAlice->GetNtrack();
  
  //   //Check, is it necessary to open new files
  //   TArrayI* events = fDigitizer->GetCurrentEvents() ; 
  //   TClonesArray * filenames = fDigitizer->GetHeadersFiles() ;
//   Int_t input ;
//   for(input = 0; input < filenames->GetEntriesFast(); input++){

//     TObjString * filename = (TObjString *) filenames->At(input) ;

//     //Test, if this file already open
//     TFile *file = (TFile*) gROOT->GetFile( filename->GetString() ) ;
//     if(file == 0)
//       file = new TFile( filename->GetString()) ;
//     file->cd() ;
    
//     // Get Kine Tree from file
// //     char treeName[20];
// //     sprintf(treeName,"TreeK%d",events->At(input));
// //     TTree * treeK = (TTree*)gDirectory->Get(treeName);
// //     if (treeK) 
// //       treeK->SetBranchAddress("Particles", &fParticleBuffer);
// //     else    
// //       cout << "AliPHOSIndexToObject: cannot find Kine Tree for event:" << events->At(input) << endl;

// //     // Create the particle stack
// //     if(!fParticles) fParticles = new TClonesArray("TParticle",1000);
// //     // Build the pointer list
// //     if(fParticleMap) {     <----
// //       fParticleMap->Clear();
// //       fParticleMap->Expand(treeK->GetEntries());
// //     } else
// //       fParticleMap = new TObjArray(treeK->GetEntries());
    
//     // From gAlice->Particle(i) 


// //   if(!(*fParticleMap)[i]) {
// //     Int_t nentries = fParticles->GetEntries();
    
// //     // algorithmic way of getting entry index
// //     // (primary particles are filled after secondaries)
// //     Int_t entry;
// //     if (i<fHeader.GetNprimary())
// //       entry = i+fHeader.GetNsecondary();
// //     else 
// //       entry = i-fHeader.GetNprimary();
      
// //     // only check the algorithmic way and give
// //     // the fatal error if it is wrong
// //     if (entry != fParticleFileMap[i]) {
// //       Fatal("Particle",
// //         "!!!! The algorithmic way is WRONG: !!!\n entry: %d map: %d",
// // 	entry, fParticleFileMap[i]); 
// //     }  
      
// //     fTreeK->GetEntry(fParticleFileMap[i]);
// //     new ((*fParticles)[nentries]) TParticle(*fParticleBuffer);
// //     fParticleMap->AddAt((*fParticles)[nentries],i);
// //   }
// //   return (TParticle *) (*fParticleMap)[i];

   
    
//   }


//   //scan over opened files and read corresponding TreeK##

  return ;
}
//____________________________________________________________________________ 
void AliPHOSIndexToObject::GetEvent(Int_t event)
{
  // Reads the content of all Tree's S, D and R

  if(event == fEvent) // do nothing
    return ;
    
  if(event > fMaxEvent){
    cout << "There is no such event " << event << " total # of events " << fMaxEvent << endl ;
    return ;
  }

  fEvent = event ;
  gAlice->GetEvent(fEvent) ;
  
  ReadTreeS() ;
  ReadTreeD() ;
  ReadTreeR() ;
  ReadPrimaries() ;

  

}

