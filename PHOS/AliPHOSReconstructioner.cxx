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
//  Algorithm class for the reconstruction: 
//                                          
//                                          
//*--
//*-- Author: Gines Martinez & Yves Schutz (SUBATECH) 
//*-- Complitely redisigned by Dmitri Peressounko (SUBATECH & RRC KI) March 2001
/////////////////////////////////////////////////////////////////////////////////////
//  Wrapping class for reconstruction
//  use case: 
//
//  root [0] AliPHOSReconstructioner * r = new AliPHOSReconstructioner("galice.root")
//              //  Set the header file
//  root [1] r->ExecuteTask() 
//
//              // One can specify the title for each branch 
//  root [2] r->SetBranchFileName("RecPoints","RecPoints1") ;
//             // By default branches are stored in galice.root (in non-split mode)
//             // or PHOS.SDigits.root, PHOS.Digits.root etc.
//      
//             // One can specify the starting point of the reconstruction
//  root [3] r->StartFrom("AliPHOSClusterizer") 
//             // means that SDigits and Digits will not be regenerated, only RecPoints, 
//             // TS and RecParticles
//
//             // And finally one can call ExecuteTask() with the following options
//  root [4] r->ExecuteTask("debug all timing")
//             // deb     - prints the numbers of produced SDigits, Digits etc.
//             // deb all - prints in addition list of made SDigits, digits etc.
//             // timing  - prints benchmarking results




// --- ROOT system ---

#include "TClonesArray.h"
#include "TROOT.h"
#include "TTree.h"

// --- Standard library ---
#include <iostream.h>   

// --- AliRoot header files ---
#include "AliRun.h"
#include "AliPHOSReconstructioner.h"
#include "AliPHOSClusterizerv1.h"
#include "AliPHOSDigitizer.h"
#include "AliPHOSSDigitizer.h"
#include "AliPHOSTrackSegmentMakerv1.h"
#include "AliPHOSPIDv1.h"
#include "AliPHOSFastRecParticle.h"
#include "AliPHOSCpvRecPoint.h"

ClassImp(AliPHOSReconstructioner)

//____________________________________________________________________________
  AliPHOSReconstructioner::AliPHOSReconstructioner():TTask("AliPHOSReconstructioner","")
{
  // ctor
  fDigitizer   = 0 ;
  fClusterizer = 0 ;
  fTSMaker     = 0 ;
  fPID         = 0 ; 
  fSDigitizer  = 0 ;

  fIsInitialized = kFALSE ;

} 

//____________________________________________________________________________
AliPHOSReconstructioner::AliPHOSReconstructioner(const char* headerFile):TTask("AliPHOSReconstructioner","")
{
  // ctor
  
  fHeaderFileName = headerFile ;

  fSDigitsBranch="" ; 
  fSDigitizer  = new AliPHOSSDigitizer(fHeaderFileName.Data(),fSDigitsBranch.Data()) ; 
  Add(fSDigitizer) ;

  fDigitsBranch="" ; 
  fDigitizer   = new AliPHOSDigitizer(fHeaderFileName.Data(),fDigitsBranch.Data()) ; 
  Add(fDigitizer) ;


  fRecPointBranch="" ; 
  fClusterizer = new AliPHOSClusterizerv1(fHeaderFileName.Data(),fRecPointBranch.Data()) ; 
  Add(fClusterizer) ;
  

  fTSBranch="" ; 
  fTSMaker     = new AliPHOSTrackSegmentMakerv1(fHeaderFileName.Data(),fTSBranch.Data()) ;
  Add(fTSMaker) ;
  
  
  fRecPartBranch="" ; 
  fPID         = new AliPHOSPIDv1(fHeaderFileName.Data(),fRecPartBranch.Data()) ;
  Add(fPID) ;
  
  fIsInitialized = kTRUE ;
  
} 
//____________________________________________________________________________
void AliPHOSReconstructioner::Exec(Option_t *option){
  //chesk, if the names of branches, which should be made conicide with already
  //existing
  if(!fIsInitialized)
    Init() ;

  gAlice->GetEvent(0) ;

  if(fSDigitizer->IsActive()&& gAlice->TreeS()){ //Will produce SDigits

    TBranch * sdigitsBranch = 0;
    TBranch * sdigitizerBranch = 0;

    TObjArray * branches = gAlice->TreeS()->GetListOfBranches() ;
    Int_t ibranch;
    Bool_t phosNotFound = kTRUE ;
    Bool_t sdigitizerNotFound = kTRUE ;

    for(ibranch = 0;ibranch <branches->GetEntries();ibranch++){            
      if(phosNotFound){
	sdigitsBranch=(TBranch *) branches->At(ibranch) ;
	if(( strcmp("PHOS",sdigitsBranch->GetName())==0 ) &&
	   (fSDigitsBranch.CompareTo(sdigitsBranch->GetTitle())== 0 ))
	  phosNotFound = kFALSE ;
      }
      if(sdigitizerNotFound){
	sdigitizerBranch = (TBranch *) branches->At(ibranch) ;
	if(( strcmp(sdigitizerBranch->GetName(),"AliPHOSSDigitizer") == 0) &&
	   (fSDigitsBranch.CompareTo(sdigitizerBranch->GetTitle())== 0 ) )
	  sdigitizerNotFound = kFALSE ;
      }
    }
    
    if(!(sdigitizerNotFound && phosNotFound)){
      cout << "AliPHOSReconstructioner error: "<< endl ;
      cout << "       Branches ''PHOS'' or ''AliPHOSSDigitizer'' with title ``" << fSDigitsBranch.Data() << "''" << endl ;
      cout << "       already exist in TreeS. ROOT does not allow updating/overwriting." << endl ;
      cout << "       Specify another title for branches or use ''StartFrom()'' method" << endl ;
      
      //mark all tasks as inactive
      TIter next(fTasks);
      TTask *task;
      while((task=(TTask*)next()))
	task->SetActive(kFALSE) ;
      
      return ;
    }
  }

  if(fDigitizer->IsActive() && gAlice->TreeD()){ //Will produce Digits
    TBranch * digitsBranch = 0;
    TBranch * digitizerBranch = 0;
    
    TObjArray * branches = gAlice->TreeD()->GetListOfBranches() ;
    Int_t ibranch;
    Bool_t phosNotFound = kTRUE ;
    Bool_t digitizerNotFound = kTRUE ;
    
    for(ibranch = 0;ibranch <branches->GetEntries();ibranch++){            
      if(phosNotFound){
	digitsBranch=(TBranch *) branches->At(ibranch) ;
	if(( strcmp("PHOS",digitsBranch->GetName())==0 ) &&
	   (fDigitsBranch.CompareTo(digitsBranch->GetTitle())== 0 ))
	  phosNotFound = kFALSE ;
      }
      if(digitizerNotFound){
	digitizerBranch = (TBranch *) branches->At(ibranch) ;
	if(( strcmp(digitizerBranch->GetName(),"AliPHOSDigitizer") == 0) &&
	   (fDigitsBranch.CompareTo(digitizerBranch->GetTitle())== 0 ) )
	  digitizerNotFound = kFALSE ;
      }
    }
    
    if(!(digitizerNotFound && phosNotFound)){
      cout << "AliPHOSReconstructioner error: "<< endl ;
      cout << "       Branches ''PHOS'' or ''AliPHOSDigitizer'' with title ``" << fDigitsBranch.Data() << "''" << endl ;
      cout << "       already exist in TreeD. ROOT does not allow updating/overwriting." << endl ;
      cout << "       Specify another title for branches or use ''StartFrom()'' method" << endl ;
      
      //mark all tasks as inactive
      TIter next(fTasks);
      TTask *task;
      while((task=(TTask*)next()))
	task->SetActive(kFALSE) ;
      
      return ;
    }
  }

  if(fClusterizer->IsActive() && gAlice->TreeR()){ //Will produce RecPoints
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
	if(fRecPointBranch.CompareTo(emcBranch->GetTitle())==0 )
	  if( strcmp(emcBranch->GetName(),"PHOSEmcRP") == 0) 
	    emcNotFound = kFALSE ;
      }
      if(cpvNotFound){
	cpvBranch=(TBranch *) branches->At(ibranch) ;
	if(fRecPointBranch.CompareTo(cpvBranch->GetTitle())==0 )
	  if( strcmp(cpvBranch->GetName(),"PHOSCpvRP") == 0) 
	    cpvNotFound = kFALSE ;
      }
      if(clusterizerNotFound){
	clusterizerBranch = (TBranch *) branches->At(ibranch) ;
	if( fRecPointBranch.CompareTo(clusterizerBranch->GetTitle()) == 0)
	  if( strcmp(clusterizerBranch->GetName(),"AliPHOSClusterizer") == 0) 
	    clusterizerNotFound = kFALSE ;
      }
    }

    if(!(clusterizerNotFound && emcNotFound && cpvNotFound)){
      cout << "AliPHOSReconstructioner error: "<< endl ;
      cout << "       Branches ''PHOSEmcRP'', ''PHOSCpvRP'' or ''AliPHOSClusterizer'' with title ``" 
	   << fRecPointBranch.Data() << "''" << endl ;
      cout << "       already exist in TreeR. ROOT does not allow updating/overwriting." << endl ;
      cout << "       Specify another title for branches or use ''StartFrom()'' method" << endl ;
      
      //mark all tasks as inactive
      TIter next(fTasks);
      TTask *task;
      while((task=(TTask*)next()))
	task->SetActive(kFALSE) ;
      return ;
    }
  }
  
  if(fTSMaker->IsActive() && gAlice->TreeR()){ //Produce TrackSegments

    TBranch * tsMakerBranch = 0;
    TBranch * tsBranch = 0;
    
    TObjArray * branches = gAlice->TreeR()->GetListOfBranches() ;
    Int_t ibranch;
    Bool_t tsMakerNotFound = kTRUE ;
    Bool_t tsNotFound = kTRUE ;
    
    for(ibranch = 0;(ibranch <branches->GetEntries())&&(tsMakerNotFound||tsNotFound);ibranch++){
      if(tsMakerNotFound){
	tsMakerBranch=(TBranch *) branches->At(ibranch) ;
	if( fTSBranch.CompareTo(tsMakerBranch->GetTitle())==0 )
	  if( strcmp(tsMakerBranch->GetName(),"AliPHOSTrackSegmentMaker") == 0) 
	    tsMakerNotFound = kFALSE ;
      }
      if(tsNotFound){
	tsBranch=(TBranch *) branches->At(ibranch) ;
	if( fTSBranch.CompareTo(tsBranch->GetTitle())==0 )
	  if( strcmp(tsBranch->GetName(),"PHOSTS") == 0) 
	    tsNotFound = kFALSE ;
      }
    }
    
    if(!(tsMakerNotFound &&tsNotFound) ){
      cout << "AliPHOSReconstructioner error: "<< endl ;
      cout << "       Branches ''PHOSTS'' or ''AliPHOSTrackSegmentMaker'' with title ``" 
	   << fTSBranch.Data() << "''" << endl ;
      cout << "       already exist in TreeR. ROOT does not allow updating/overwriting." << endl ;
      cout << "       Specify another title for branches or use ''StartFrom()'' method" << endl ;
      
      //mark all tasks as inactive
      TIter next(fTasks);
      TTask *task;
      while((task=(TTask*)next()))
	task->SetActive(kFALSE) ;
      return ;
      
    }
    
  }

  if(fPID->IsActive() && gAlice->TreeR()){ //Produce RecParticles
    TBranch * pidBranch = 0;
    TBranch * rpBranch = 0;
    
    TObjArray * branches = gAlice->TreeR()->GetListOfBranches() ;
    Int_t ibranch;
    Bool_t pidNotFound = kTRUE ;
    Bool_t rpNotFound = kTRUE ;
    
    for(ibranch = 0;(ibranch <branches->GetEntries()) && pidNotFound && rpNotFound ;ibranch++){
      if(pidNotFound){
	pidBranch=(TBranch *) branches->At(ibranch) ;
	if( (strcmp(fRecPartBranch,pidBranch->GetTitle())==0 ) &&
	    (strcmp(pidBranch->GetName(),"AliPHOSPID") == 0) )
	  pidNotFound = kFALSE ;
      }
      if(rpNotFound){
	rpBranch=(TBranch *) branches->At(ibranch) ;
	if( (strcmp(fRecPartBranch,rpBranch->GetTitle())==0 ) &&
	    (strcmp(rpBranch->GetName(),"PHOSRP") == 0) )
	  rpNotFound = kFALSE ;
      }
    }
    
    if(!pidNotFound  || !rpNotFound ){
      cout << "AliPHOSReconstructioner error: "<< endl ;
      cout << "       Branches ''PHOSRP'' or ''AliPHOSPID'' with title ``" 
	   << fRecPartBranch.Data() << "''" << endl ;
      cout << "       already exist in TreeR. ROOT does not allow updating/overwriting." << endl ;
      cout << "       Specify another title for branches." << endl ;
      
      //mark all tasks as inactive
      TIter next(fTasks);
      TTask *task;
      while((task=(TTask*)next()))
	task->SetActive(kFALSE) ;
      return ;
    }
    
  }
}
//____________________________________________________________________________
 void AliPHOSReconstructioner::Init()
{
  //initiase Reconstructioner if necessary: we can not do this in default constructor

  if(!fIsInitialized){
    // Initialisation

    fSDigitsBranch="" ; 
    fSDigitizer  = new AliPHOSSDigitizer(fHeaderFileName.Data(),fSDigitsBranch.Data()) ; 
    Add(fSDigitizer) ;

    fDigitsBranch="" ; 
    fDigitizer   = new AliPHOSDigitizer(fHeaderFileName.Data(),fDigitsBranch.Data()) ; 
    Add(fDigitizer) ;

    fRecPointBranch="" ; 
    fClusterizer = new AliPHOSClusterizerv1(fHeaderFileName.Data(),fRecPointBranch.Data()) ; 
    Add(fClusterizer) ;

    fTSBranch="" ; 
    fTSMaker     = new AliPHOSTrackSegmentMakerv1(fHeaderFileName.Data(),fTSBranch.Data()) ;
    Add(fTSMaker) ;


    fRecPartBranch="" ; 
    fPID         = new AliPHOSPIDv1(fHeaderFileName.Data(),fRecPartBranch.Data()) ;
    Add(fPID) ;
    
    fIsInitialized = kTRUE ;
  }
} 
//____________________________________________________________________________
AliPHOSReconstructioner::~AliPHOSReconstructioner()
{
  
  if(fSDigitizer)
    delete fSDigitizer ;
  
  if(fDigitizer)
    delete fDigitizer ;
  
  if(fClusterizer)
    delete fClusterizer ;
  
  if(fTSMaker)
    delete fTSMaker ;
  
  if(fPID)
    delete fPID ;
} 
//____________________________________________________________________________
void AliPHOSReconstructioner::SetBranchTitle(const char* branch, const char * title){
  //Diverge correcpoinding branch to the file "title"

  if(strcmp(branch,"SDigits") == 0){ 
    fSDigitizer->SetSDigitsBranch(title) ;
    fDigitizer->SetSDigitsBranch(title) ;
    fSDigitsBranch = title ;
    return ;
  }
  
  if(strcmp(branch,"Digits") == 0){ 
    fDigitizer->SetDigitsBranch(title) ;
    fClusterizer->SetDigitsBranch(title) ;
    fDigitsBranch = title ;
    return ;
  }

  if(strcmp(branch,"RecPoints") == 0){ 
    fClusterizer->SetRecPointsBranch(title) ;
    fTSMaker->SetRecPointsBranch(title) ;
    fRecPointBranch = title ;
    return ;
  }

  if(strcmp(branch,"TrackSegments") == 0){
    fTSMaker->SetTrackSegmentsBranch(title) ;
    fPID->SetTrackSegmentsBranch(title) ;
    fTSBranch = title ;
    return ;
  }

  if(strcmp(branch,"RecParticles") == 0){ 
    fPID->SetRecParticlesBranch(title) ;
    fRecPartBranch = title ;
    return ;
  }

  cout << "There is no branch " << branch << "!"<< endl ;
  cout << "Available branches `SDigits', `Digits', `RecPoints', `TrackSegments' and `RecParticles' " << endl ;
  
}
//____________________________________________________________________________
void AliPHOSReconstructioner::StartFrom(Option_t * module){
  //in the next ExecuteTask() reconstruction starts from the module "module"

  if(!fIsInitialized)
    Init() ;  
  TIter next(fTasks);
  TTask *task;
  Bool_t active = kFALSE ;
  while((task=(TTask*)next())){ 
    if (strcmp(module,task->GetName())==0)  
      active = kTRUE;
    task->SetActive(active) ;
  }
  if(!active){
    cout << "There is no task " <<module<< endl ;
    cout << "Available tasks are: " << endl ;
    next.Reset() ;
    while((task=(TTask*)next()))
      cout<<"                    " << task->GetName() << endl ;
  }
}

//____________________________________________________________________________
void AliPHOSReconstructioner::Print(Option_t * option)const {
  
  cout << "-----------------AliPHOSReconstructioner---------------" << endl ;
  cout << " Reconstruction of the header file " <<fHeaderFileName.Data() << endl ;
  cout << " with the following modules: " << endl ;

  if(fSDigitizer->IsActive()){
    cout << "   (+)   " << fSDigitizer->GetName() << " to branch : " << fSDigitsBranch.Data() << endl ; 
    cout << endl ;
  }
  if(fDigitizer->IsActive()){
    cout << "   (+)   " << fDigitizer->GetName() << " to branch : " << fDigitsBranch.Data() << endl ;  
    cout <<  endl ;
  }
  
  if(fClusterizer->IsActive()){
    cout << "   (+)   " <<fClusterizer->GetName() << " to branch : " <<fRecPointBranch.Data()  << endl ;  
    cout <<  endl ;
  }

  if(fTSMaker->IsActive()){
    cout << "   (+)   " << fTSMaker->GetName() << " to branch : " << fTSBranch.Data() << endl ;  
    cout <<  endl ;
  }


  if(fPID->IsActive()){
    cout << "   (+)   " << fPID->GetName() << " to branch : " <<fRecPartBranch.Data()  << endl ;  
    cout <<  endl ;
  }


}
