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

//_________________________________________________________________________
//  A singleton 
//                  
//*-- Author: Yves Schutz (SUBATECH) & Dmitri Peressounko (RRC KI & SUBATECH)
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
#include "AliPHOSClusterizer.h"
#include "AliPHOSTrackSegmentMaker.h"
#include "AliPHOSPID.h" 

ClassImp(AliPHOSIndexToObject)
  
  AliPHOSIndexToObject * AliPHOSIndexToObject::fgObjGetter = 0 ; 

//____________________________________________________________________________ 
AliPHOSIndexToObject::AliPHOSIndexToObject(char* headerFile,char* branch,char* branchTitle )
{
  //Initiate all lists
  fEvent = 0 ;

  fDigits = new TClonesArray("AliPHOSDigit",100) ;
  fEmcRecPoints = new TObjArray(100) ;
  fCpvRecPoints = new TObjArray(100) ;
  fTS = new TClonesArray("AliPHOSTrackSegment",100) ;
  fRecParticles = new TClonesArray("AliPHOSRecParticle",100) ;
  fPrimaries = new TObjArray(1) ;

  fDigitizer = 0 ;
  fClusterizer = 0 ;
  fTSMaker = 0 ;
  fPID = 0 ;

  //open headers file
  fHeaderFile = headerFile ;
  TFile * file = (TFile*) gROOT->GetFile(fHeaderFile.Data() ) ;

  if(file == 0){
    file = new TFile(fHeaderFile.Data(),"update") ;
    gAlice = (AliRun *) file->Get("gAlice") ;
  }

  fMaxEvent = (Int_t) gAlice->TreeE()->GetEntries() ;

  char * dummyfile = 0 ;

  gAlice->GetEvent(fEvent) ;

  Bool_t isRead = kFALSE;
  //now read branches 
  if((strcmp(branch,"PHOSRP")==0) || (strcmp(branch,"PHOSPID")==0)){
    ReadRecParticles(branchTitle) ;  //first read RecPartcles and branche TS from which they are made
    ReadTS(dummyfile);              //read TS from which made RecParticles above
    ReadRecPoints(dummyfile) ;     //RecPoints from which TS above made
    ReadDigits(dummyfile) ;         //digits. from whic RecPoints made
    isRead= kTRUE ;
  }
  
  if((strcmp(branch,"PHOSTS")==0) || (strcmp(branch,"PHOSTSMaker")==0)){
    ReadTS(branchTitle);            //read TS and branch of RecPoints from which they are made
    ReadRecPoints(dummyfile) ;     //recpoints abd branch of digits
    ReadDigits(dummyfile) ;       //digits and branch of Primaries
    ReadRecParticles(dummyfile) ;  //posiible completion of TS
    isRead= kTRUE ;
  }

  if((strcmp(branch,"PHOSEmcRP")==0)|| (strcmp(branch,"PHOSCpvRP")==0) || 
     (strcmp(branch,"PHOSClusterizer")==0)){
    ReadRecPoints(branchTitle) ;    //RecPoints and Digits branch filename
    ReadDigits(dummyfile) ;        //digits and primary file name
    ReadTS(dummyfile);             //possible completion of RecPoints
    ReadRecParticles(dummyfile) ;  //possible completion of TS
    isRead= kTRUE ;
  }

  if((strcmp(branch,"PHOS")==0) || (strcmp(branch,"PHOSDigitizer")==0)){
    ReadDigits(branchTitle) ;
    ReadRecPoints(dummyfile) ;
    ReadTS(dummyfile);
    ReadRecParticles(dummyfile) ;
    isRead= kTRUE ;
  }

  if(!isRead){
    cout << "AliPHOSIndexToObject: wrong branch name specified: " << branch << endl ;
    cout << "   avalilable names are `PHOSRP', `PHOSPID'"<<endl ;
    cout << "                        `PHOSTS', `PHOSTSMaker'"<<endl ;
    cout << "                        `PHOSEmcRP', `PHOSCpvRP', `PHOSClusterizer'"<< endl ;
    cout << "                        `PHOS' and `PHOSDigitizer'"<< endl ;
  }
  ReadPrimaries() ; // should be called when digits are already read 

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
AliPHOSIndexToObject * AliPHOSIndexToObject::GetInstance(char* headerFile,char* branch,char* branchTitle)
{
  // Creates and returns the pointer of the unique instance
  // Must be called only when the environment has changed 

  if ( fgObjGetter )      // delete it if already exists
    delete fgObjGetter ; 

  fgObjGetter = new AliPHOSIndexToObject(headerFile,branch,branchTitle) ; 
  
  return fgObjGetter ; 
  
}

//____________________________________________________________________________ 
TParticle * AliPHOSIndexToObject::GimePrimary(Int_t index)
{
  
  
  Int_t primaryList = (Int_t) (TMath::Ceil(index/10000000.) ) - 1 ;
  Int_t primaryIndex = index - primaryList*10000000 ; 
  
  if ( primaryList > 0  ) {
    cout << " IndexToObject does not support currently Mixing of primary " << endl ;
    cout << "   can not return primary: " << index<< " (list "<< primaryList<< " primary # " << primaryIndex << " )"<<endl ;
    return 0;
  }
  
  return gAlice->Particle(primaryIndex) ;
  
}

//____________________________________________________________________________ 
Bool_t AliPHOSIndexToObject::ReadRecParticles(char * branchTitle){

  if(gAlice->TreeR()==0)
    return kFALSE ;
  
  if(fPID) // already read
    branchTitle = fPID->GetRecParticlesBranch() ;


  if(branchTitle){ // we should read a specific branch
    TBranch * pidBranch = 0;
    TBranch * rpBranch = 0;
    
    TObjArray * branches = gAlice->TreeR()->GetListOfBranches() ;
    Int_t ibranch;
    Bool_t pidNotFound = kTRUE ;
    Bool_t rpNotFound = kTRUE ;
    
    for(ibranch = 0;(ibranch <branches->GetEntries())&&(pidNotFound||rpNotFound);ibranch++){

      if(pidNotFound){
	pidBranch=(TBranch *) branches->At(ibranch) ;
	if( (strcmp(branchTitle,pidBranch->GetTitle())==0 ) &&
	    (strcmp(pidBranch->GetName(),"AliPHOSPID") == 0) )
	  pidNotFound = kFALSE ;
      }
      if(rpNotFound){
	rpBranch=(TBranch *) branches->At(ibranch) ;
	if( (strcmp(branchTitle,rpBranch->GetTitle())==0 ) &&
	    (strcmp(rpBranch->GetName(),"PHOSRP") == 0) )
	  rpNotFound = kFALSE ;
      }
    }
    
    if(pidNotFound ||rpNotFound ){
      cout << "AliPHOSIndexToObject error" << endl ;
      cout << "     Can't find Branch with PID and RecParticles " ;
      return kFALSE ;
    }
    
    pidBranch->SetAddress(&fPID) ;
    rpBranch->SetAddress(&fRecParticles) ;
    gAlice->TreeR()->GetEvent(0) ;    
  }
  else{ //we Should read any branch and print warning if there are other possibilities
    if(fTSMaker){//if TrackSegments already read, we should read RecParticles Made from it
      TBranch * pidBranch = 0;
      TBranch * rpBranch = 0;
    
      TObjArray * branches = gAlice->TreeR()->GetListOfBranches() ;

      Int_t branchRead = 0;
      Bool_t allNotFound = kTRUE ;
      while(allNotFound){
	Bool_t pidNotFound = kTRUE ;
	Bool_t rpNotFound = kTRUE ;
	Int_t ibranch ;
	for(ibranch = branchRead;(ibranch <branches->GetEntries() )&& pidNotFound;ibranch++){
	  pidBranch=(TBranch *) branches->At(ibranch) ;
	  if(strcmp(pidBranch->GetName(),"AliPHOSPID") == 0) 
	    pidNotFound = kFALSE ;
	}
	branchRead = ibranch +1 ; 
	for(ibranch = 0 ;(ibranch <branches->GetEntries() )&& rpNotFound;ibranch++){
	  rpBranch=(TBranch *) branches->At(ibranch) ;
	  if( (strcmp(pidBranch->GetTitle(),rpBranch->GetTitle())==0 ) &&
	      (strcmp(rpBranch->GetName(),"PHOSRP") == 0) )
	    rpNotFound = kFALSE ;
	}
	
	if(pidNotFound ||rpNotFound ){
	  cout << "AliPHOSIndexToObject error" << endl ;
	  cout << "     Can't find Branch with PID and RecParticles " ;
	  return kFALSE ;
	}
    
	pidBranch->SetAddress(&fPID) ;
	rpBranch->SetAddress(&fRecParticles) ;
	gAlice->TreeR()->GetEvent(0) ;    
	
	if(strcmp(fTSMaker->GetTrackSegmentsBranch(),fPID->GetTrackSegmentsBranch()) == 0)
	  allNotFound = kFALSE ;
      }
    }
    else{//we read any (first) recparticles
      TBranch * pidBranch = 0;
      TBranch * rpBranch = 0;
      TObjArray * branches = gAlice->TreeR()->GetListOfBranches() ;

      Bool_t pidNotFound = kTRUE ;
      Bool_t rpNotFound = kTRUE ;
      Int_t ibranch ;
      for(ibranch = 0;(ibranch <branches->GetEntries() )&& pidNotFound;ibranch++){
	pidBranch=(TBranch *) branches->At(ibranch) ;
	if(strcmp(pidBranch->GetName(),"AliPHOSPID") == 0) 
	  pidNotFound = kFALSE ;
      }
      for(ibranch = 0 ;(ibranch <branches->GetEntries() )&& rpNotFound;ibranch++){
	rpBranch=(TBranch *) branches->At(ibranch) ;
	if( (strcmp(pidBranch->GetTitle(),rpBranch->GetTitle())==0 ) &&
	    (strcmp(rpBranch->GetName(),"PHOSRP") == 0) )
	  rpNotFound = kFALSE ;
      }
      
      if(pidNotFound ||rpNotFound ){
	cout << "AliPHOSIndexToObject worning: " << endl ;
	cout << "     Can't find Branch with PID and RecParticles " << endl;
	return kFALSE ;
      }
      
      pidBranch->SetAddress(&fPID) ;
      rpBranch->SetAddress(&fRecParticles) ;
      gAlice->TreeR()->GetEvent(0) ;    
      
    }
  }
  return kTRUE ;
}
//____________________________________________________________________________ 
Bool_t AliPHOSIndexToObject::ReadTS(char * branchTitle){

  if(gAlice->TreeR()==0)
    return kFALSE ;

  if(fPID)//if RecParticles already read, we should read TS from which they are made
    branchTitle= fPID->GetTrackSegmentsBranch() ;
  
  if(branchTitle){   // we should read a specific branch
    
    TBranch * tsMakerBranch = 0;
    TBranch * tsBranch = 0;
    
    TObjArray * branches = gAlice->TreeR()->GetListOfBranches() ;
    Int_t ibranch;
    Bool_t tsMakerNotFound = kTRUE ;
    Bool_t tsNotFound = kTRUE ;
    
    for(ibranch = 0;(ibranch <branches->GetEntries())&&(tsMakerNotFound||tsNotFound);ibranch++){
      if(tsMakerNotFound){
	tsMakerBranch=(TBranch *) branches->At(ibranch) ;
	if( strcmp(branchTitle,tsMakerBranch->GetTitle())==0 )
	  if( strcmp(tsMakerBranch->GetName(),"AliPHOSTrackSegmentMaker") == 0) 
	    tsMakerNotFound = kFALSE ;
      }
      if(tsNotFound){
	tsBranch=(TBranch *) branches->At(ibranch) ;
	if( strcmp(branchTitle,tsBranch->GetTitle())==0 )
	  if( strcmp(tsBranch->GetName(),"PHOSTS") == 0) 
	    tsNotFound = kFALSE ;
      }
    }
    
    if(tsMakerNotFound ||tsNotFound ){
      cout << "AliPHOSIndexToObject error" << endl ;
      cout << "       Can't find Branch with TrackSegmentMaker and TrackSegments " ;
      cout << "       Do nothing" <<endl  ;
      return kFALSE ;
    }
    
    tsMakerBranch->SetAddress(&fTSMaker) ;
    tsBranch->SetAddress(&fTS) ;
    gAlice->TreeR()->GetEvent(0) ;
    
  }
  else{ 
    if(fClusterizer){//Clusterizer aready read, 
                     //we should read TrackSegments made from these RecPoints

      Int_t branchRead = 0 ; 
      TBranch * tsMakerBranch = 0;
      TBranch * tsBranch = 0;
    
      TObjArray * branches = gAlice->TreeR()->GetListOfBranches() ;
      Int_t ibranch;
      Bool_t allNotFound = kTRUE ;
      while(allNotFound){
	Bool_t tsMakerNotFound = kTRUE ;
	Bool_t tsNotFound = kTRUE ;
	
	for(ibranch = branchRead;(ibranch <branches->GetEntries())&&(tsMakerNotFound);ibranch++){
	  tsMakerBranch=(TBranch *) branches->At(ibranch) ;
	  if( strcmp(tsMakerBranch->GetName(),"AliPHOSTrackSegmentMaker") == 0) 
	    tsMakerNotFound = kFALSE ;
	}
	branchRead = ibranch++ ;
	for(ibranch = 0 ;(ibranch <branches->GetEntries())&&(tsNotFound);ibranch++){
	  tsBranch=(TBranch *) branches->At(ibranch) ;
	  if( (strcmp(tsBranch->GetName(),"PHOSTS") == 0) && 
	      (strcmp(tsBranch->GetName(),tsMakerBranch->GetTitle())==0))
	    tsNotFound = kFALSE ;
	}
	
	branchRead = ibranch++ ;
	
	if(tsMakerNotFound ||tsNotFound ){
	  cout << "AliPHOSIndexToObject error" << endl ;
	  cout << "       Can't find Branch with TrackSegmentMaker and TrackSegments " ;
	  cout << "       Do nothing" <<endl  ;
	  return kFALSE ;
	}
	
	tsMakerBranch->SetAddress(&fTSMaker) ;
	tsBranch->SetAddress(&fTS) ;
	gAlice->TreeR()->GetEvent(0) ;
	
	if(strcmp(fTSMaker->GetRecPointsBranch(),fClusterizer->GetRecPointsBranch()) == 0)
	  allNotFound = kFALSE ;
      }
      
    }
    else{//Neither Title,neither fPID, neither fClusterizer: we read any (first) occurence
      TBranch * tsMakerBranch = 0;
      TBranch * tsBranch = 0;    
      TObjArray * branches = gAlice->TreeR()->GetListOfBranches() ;
      Bool_t tsMakerNotFound = kTRUE ;
      Bool_t tsNotFound = kTRUE ;
      Int_t ibranch ;
      for(ibranch =  0;(ibranch <branches->GetEntries())&& tsMakerNotFound;ibranch++){
	tsMakerBranch=(TBranch *) branches->At(ibranch) ;
	if( strcmp(tsMakerBranch->GetName(),"AliPHOSTrackSegmentMaker") == 0) 
	  tsMakerNotFound = kFALSE ;
      }
      for(ibranch = 0 ;(ibranch <branches->GetEntries())&&(tsNotFound);ibranch++){
	tsBranch=(TBranch *) branches->At(ibranch) ;
	if( (strcmp(tsBranch->GetName(),"PHOSTS") == 0) && 
	    (strcmp(tsBranch->GetName(),tsMakerBranch->GetTitle())==0))
	  tsNotFound = kFALSE ;
      }	
      if(tsMakerNotFound ||tsNotFound ){
	cout << "AliPHOSIndexToObject error" << endl ;
	cout << "       Can't find Branch with TrackSegmentMaker and TrackSegments " ;
	cout << "       Do nothing" <<endl  ;
	return kFALSE ;
      }
      
      tsMakerBranch->SetAddress(&fTSMaker) ;
      tsBranch->SetAddress(&fTS) ;
      gAlice->TreeR()->GetEvent(0) ;     
    }
  }

  return kTRUE ;  
}
//____________________________________________________________________________ 
Bool_t AliPHOSIndexToObject::ReadRecPoints(char * branchTitle){
  
  if(gAlice->TreeR() == 0)
    return kFALSE ;

  if(fTSMaker) //if TrackSegment maker already read, read corresponding branches
    branchTitle = fTSMaker->GetRecPointsBranch() ;

  if(branchTitle){ // we should read a specific branch
    TBranch * emcBranch = 0;
    TBranch * cpvBranch = 0;
    TBranch * clusterizerBranch = 0;
    
    TObjArray * branches = gAlice->TreeR()->GetListOfBranches() ;
    Int_t ibranch;
    Bool_t emcNotFound = kTRUE ;
    Bool_t cpvNotFound = kTRUE ;  
    Bool_t clusterizerNotFound = kTRUE ;
    
    for(ibranch = 0;((ibranch < branches->GetEntries())&&(emcNotFound ||cpvNotFound || clusterizerNotFound)) ;ibranch++){
      if(emcNotFound){
	emcBranch=(TBranch *) branches->At(ibranch) ;
	if( (strcmp(emcBranch->GetTitle(),branchTitle) == 0) && 
	    (strcmp(emcBranch->GetName(),"PHOSEmcRP") == 0) )
	  emcNotFound = kFALSE ;
      }
      if(cpvNotFound){
	cpvBranch=(TBranch *) branches->At(ibranch) ;
	if( (strcmp(cpvBranch->GetTitle(),branchTitle) == 0) &&
	    (strcmp(cpvBranch->GetName(),"PHOSCpvRP") == 0) )
	  cpvNotFound = kFALSE ;
      }
      if(clusterizerNotFound){
	clusterizerBranch = (TBranch *) branches->At(ibranch) ;
	if( (strcmp(clusterizerBranch->GetTitle(),branchTitle) == 0) &&
	    (strcmp(clusterizerBranch->GetName(),"AliPHOSClusterizer") == 0) )
	  clusterizerNotFound = kFALSE ;
      }
      
    }
    
    if(clusterizerNotFound || emcNotFound || cpvNotFound){
      cout << "AliPHOSIndexToObject error" << endl ;
      cout << "       Can't find Branch with RecPoints or Clusterizer " << endl ;
      return kFALSE ;
    }
    
    emcBranch->SetAddress(&fEmcRecPoints) ;
    cpvBranch->SetAddress(&fCpvRecPoints) ;
    clusterizerBranch->SetAddress(&fClusterizer) ;
    gAlice->TreeR()->GetEvent(0) ;
  }
  else{ //no specific branch
    if(fDigitizer){//Digitizer aready read, 
                   //we should read RecPoints made from these Digits
      TBranch * emcBranch = 0;
      TBranch * cpvBranch = 0;
      TBranch * clusterizerBranch = 0;
      
      TObjArray * branches = gAlice->TreeR()->GetListOfBranches() ;
      Int_t branchRead = 0;
      Bool_t allNotFound = kTRUE ;
      while(allNotFound){
	Bool_t emcNotFound = kTRUE ;
	Bool_t cpvNotFound = kTRUE ;  
	Bool_t clusterizerNotFound = kTRUE ;
	Int_t ibranch ;
	for(ibranch = branchRead;ibranch < branches->GetEntries();ibranch++){
	  emcBranch=(TBranch *) branches->At(ibranch) ;
	  if( strcmp(emcBranch->GetName(),"PHOSEmcRP") == 0)
	    emcNotFound = kFALSE ;
	}
	branchRead = ibranch + 1 ;
	for(ibranch =  0 ;ibranch < branches->GetEntries();ibranch++){
	  cpvBranch=(TBranch *) branches->At(ibranch) ;
	  if( (strcmp(cpvBranch->GetTitle(),emcBranch->GetTitle()) == 0) &&
	      (strcmp(cpvBranch->GetName(),"PHOSCpvRP") == 0) )
	    cpvNotFound = kFALSE ;
	}
	for(ibranch = 0 ;ibranch < branches->GetEntries();ibranch++){
	  clusterizerBranch = (TBranch *) branches->At(ibranch) ;
	  if( (strcmp(clusterizerBranch->GetTitle(),emcBranch->GetTitle()) == 0) &&
	      (strcmp(clusterizerBranch->GetName(),"AliPHOSClusterizer") == 0) )
	    clusterizerNotFound = kFALSE ;
	}
	
	if(clusterizerNotFound || emcNotFound || cpvNotFound){
	  cout << "AliPHOSIndexToObject error" << endl ;
	  cout << "       Can't find Branch with RecPoints or Clusterizer " << endl ;
	  return kFALSE ;
	}
    
	emcBranch->SetAddress(&fEmcRecPoints) ;
	cpvBranch->SetAddress(&fCpvRecPoints) ;
	clusterizerBranch->SetAddress(&fClusterizer) ;
	gAlice->TreeR()->GetEvent(0) ;
	
	if(strcmp(fClusterizer->GetDigitsBranch(),fDigitizer->GetDigitsBranch())== 0)
	  allNotFound = kFALSE ;
      }
    }
    else{//Neither Title, Neither TSMaker, Neither Digits: we read any (first) RecPoints
      TBranch * emcBranch = 0;
      TBranch * cpvBranch = 0;
      TBranch * clusterizerBranch = 0;
      TObjArray * branches = gAlice->TreeR()->GetListOfBranches() ;
      Bool_t emcNotFound = kTRUE ;
      Bool_t cpvNotFound = kTRUE ;  
      Bool_t clusterizerNotFound = kTRUE ;
      Int_t ibranch ;
      for(ibranch = 0 ;ibranch < branches->GetEntries();ibranch++){
	emcBranch=(TBranch *) branches->At(ibranch) ;
	if( strcmp(emcBranch->GetName(),"PHOSEmcRP") == 0)
	  emcNotFound = kFALSE ;
      }
      for(ibranch = 0 ;ibranch < branches->GetEntries();ibranch++){
	cpvBranch=(TBranch *) branches->At(ibranch) ;
	if( (strcmp(cpvBranch->GetTitle(),emcBranch->GetTitle()) == 0) &&
	    (strcmp(cpvBranch->GetName(),"PHOSCpvRP") == 0) )
	  cpvNotFound = kFALSE ;
      }
      for(ibranch = 0;ibranch < branches->GetEntries();ibranch++){
	clusterizerBranch = (TBranch *) branches->At(ibranch) ;
	if( (strcmp(clusterizerBranch->GetTitle(),emcBranch->GetTitle()) == 0) &&
	    (strcmp(clusterizerBranch->GetName(),"AliPHOSClusterizer") == 0) )
	  clusterizerNotFound = kFALSE ;
      }
      
      if(clusterizerNotFound || emcNotFound || cpvNotFound){
	cout << "AliPHOSIndexToObject error" << endl ;
	cout << "       Can't find Branch with RecPoints or Clusterizer " << endl ;
	return kFALSE ;
      }
      
      emcBranch->SetAddress(&fEmcRecPoints) ;
      cpvBranch->SetAddress(&fCpvRecPoints) ;
      clusterizerBranch->SetAddress(&fClusterizer) ;
      gAlice->TreeR()->GetEvent(0) ;
    }
  }

  return kTRUE ;
}
//____________________________________________________________________________ 
Bool_t AliPHOSIndexToObject::ReadDigits(char * branchTitle){

  if(gAlice->TreeD()== 0)
    return kFALSE ;
  
  //if RecPoints are already read, we should read Digits from which they are made
  if(fClusterizer)
    branchTitle = fClusterizer->GetDigitsBranch() ;
  
  if(branchTitle){ // we should read a specific branch
    TBranch * digitsBranch = 0;
    TBranch * digitizerBranch = 0;

    TObjArray * branches = gAlice->TreeD()->GetListOfBranches() ;
    Int_t ibranch;
    Bool_t phosNotFound = kTRUE ;
    Bool_t digitizerNotFound = kTRUE ;
    
    for(ibranch = 0;ibranch <branches->GetEntries();ibranch++){
      
      if(phosNotFound){
	digitsBranch=(TBranch *) branches->At(ibranch) ;
	if( (strcmp(digitsBranch->GetTitle(),branchTitle)==0 ) &&
	    (strcmp(digitsBranch->GetName(),"PHOS") == 0) )
	  phosNotFound = kFALSE ;
      }
      if(digitizerNotFound){
	digitizerBranch = (TBranch *) branches->At(ibranch) ;
	if( (strcmp(digitizerBranch->GetTitle(),branchTitle) == 0) && 
	    (strcmp(digitizerBranch->GetName(),"AliPHOSDigitizer") == 0) )
	  digitizerNotFound = kFALSE ;
      } 
    }
    
    if(digitizerNotFound || phosNotFound){
      cout << "AliPHOSIndexToObject error: " << endl ;
      cout << "       Can't find Branch with Digits or Digitizer "<< endl ; ;
      return kFALSE ;
    }
    
    digitsBranch->SetAddress(&fDigits) ;
    digitizerBranch->SetAddress(&fDigitizer) ;
  
    gAlice->TreeD()->GetEvent(0) ;
  }
  else{ //we should read any branch and print warning if there are other possibilities
    TBranch * digitsBranch = 0;
    TBranch * digitizerBranch = 0;
    
    TObjArray * branches = gAlice->TreeD()->GetListOfBranches() ;
    Int_t ibranch;
    Bool_t phosNotFound = kTRUE ;
    Bool_t digitizerNotFound = kTRUE ;
    
    for(ibranch = 0;(ibranch <branches->GetEntries())&& phosNotFound ;ibranch++){
      digitsBranch=(TBranch *) branches->At(ibranch) ;
      if(strcmp(digitsBranch->GetName(),"PHOS") == 0) 
	phosNotFound = kFALSE ;
    }
    for(ibranch = 0;(ibranch <branches->GetEntries())&& digitizerNotFound ;ibranch++){
      digitizerBranch = (TBranch *) branches->At(ibranch) ;
      if( (strcmp(digitizerBranch->GetTitle(),digitsBranch->GetTitle()) == 0) && 
	  (strcmp(digitizerBranch->GetName(),"AliPHOSDigitizer") == 0) )
	digitizerNotFound = kFALSE ;
    } 
    
    if(digitizerNotFound || phosNotFound){
      cout << "AliPHOSIndexToObject error: " << endl ;
      cout << "       Can't find Branch with Digits or Digitizer "<< endl ; ;
      return kFALSE ;
    }
    
    digitsBranch->SetAddress(&fDigits) ;
    digitizerBranch->SetAddress(&fDigitizer) ;
    
    gAlice->TreeD()->GetEvent(0) ;
    
  }

  return kTRUE ;
}
//____________________________________________________________________________ 
Bool_t AliPHOSIndexToObject::ReadPrimaries(){
  //read specific branches of primaries

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

  return kTRUE ;
}
//____________________________________________________________________________ 
void AliPHOSIndexToObject::GetEvent(Int_t event){
  if(event == fEvent) // do nothing
    return ;
    
  fEvent = event ;
  gAlice->GetEvent(fEvent) ;
  
  ReadRecParticles(fPID->GetRecParticlesBranch()) ;
  ReadTS(fTSMaker->GetTrackSegmentsBranch()) ;
  ReadRecPoints(fClusterizer->GetRecPointsBranch()) ;
  ReadDigits(fDigitizer->GetDigitsBranch()) ;
  ReadPrimaries() ;
}

