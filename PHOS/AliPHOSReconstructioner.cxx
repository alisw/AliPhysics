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
//*--
//*-- Author: Gines Martinez & Yves Schutz (SUBATECH) 
//*-- Compleetely redesigned by Dmitri Peressounko (SUBATECH & RRC KI) March 2001
/////////////////////////////////////////////////////////////////////////////////////
//  Wrapping class for reconstruction. Allows to produce reconstruction from 
//  different steps: from previously produced hits,sdigits, etc. Each new reconstruction
//  flow (e.g. digits, made from them RecPoints, subsequently made TrackSegments, 
//  subsequently made RecParticles) are distinguished by the title of created branches. One can 
//  use this title as a comment, see use case below. 
//  Thanks to getters, one can set 
//  parameters to reconstruction briks. The full set of parameters is saved in the 
//  corresponding branch: e.g. parameters of clusterizer are stored in branch 
//  TreeR::AliPHOSClusterizer with the same title as the branch containing the RecPoints. //  TTree does not support overwriting, therefore one can not produce several 
//  branches with the same names and titles - use different titles.
//
//  Use case: 
//
//  root [0] AliPHOSReconstructioner * r = new AliPHOSReconstructioner("galice.root")
//              //  Set the header file
//  root [1] r->ExecuteTask() 
//              //  Make full chain of reconstruction
//
//              // One can specify the title for each branch 
//  root [2] r->SetBranchFileName("RecPoints","RecPoints1") ;
//      
//             // One can change parameters of reconstruction algorithms
//  root [3] r->GetClusterizer()->SetEmcLocalMaxCut(0.02)
//
//             // One can specify the starting point of the reconstruction and title of all 
//             // branches produced in this pass
//  root [4] r->StartFrom("AliPHOSClusterizer","Local max cut 0.02") 
//             // means that will use already generated Digits and produce only RecPoints, 
//             // TS and RecParticles 
//
//             // And finally one can call ExecuteTask() with the following options
//  root [5] r->ExecuteTask("debug all timing")
//             // deb     - prints the numbers of produced SDigits, Digits etc.
//             // deb all - prints in addition list of made SDigits, digits etc.
//             // timing  - prints benchmarking results
///////////////////////////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

#include "TClonesArray.h"
#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"

// --- Standard library ---

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
  fToSplit = kFALSE ;
  fDigitizer   = 0 ;
  fClusterizer = 0 ;
  fTSMaker     = 0 ;
  fPID         = 0 ; 
  fSDigitizer  = 0 ;
  fHeaderFileName = "galice.root" ;

  fIsInitialized = kFALSE ;

} 

//____________________________________________________________________________
AliPHOSReconstructioner::AliPHOSReconstructioner(const char* headerFile,const char * branchName,Bool_t toSplit):
TTask("AliPHOSReconstructioner","")
{
  // ctor
  
  fHeaderFileName = headerFile ;
  fToSplit = toSplit ;
  fSDigitsBranch= branchName; 
  fSDigitizer  = new AliPHOSSDigitizer(fHeaderFileName.Data(),fSDigitsBranch.Data(),toSplit) ; 
  Add(fSDigitizer) ;

  fDigitsBranch=branchName ; 
  fDigitizer   = new AliPHOSDigitizer(fHeaderFileName.Data(),fDigitsBranch.Data(),toSplit) ; 
  Add(fDigitizer) ;


  fRecPointBranch=branchName ; 
  fClusterizer = new AliPHOSClusterizerv1(fHeaderFileName.Data(),fRecPointBranch.Data(),toSplit) ; 
  Add(fClusterizer) ;
  

  fTSBranch=branchName ; 
  fTSMaker     = new AliPHOSTrackSegmentMakerv1(fHeaderFileName.Data(),fTSBranch.Data(),toSplit) ;
  Add(fTSMaker) ;
  
  
  fRecPartBranch=branchName ; 
  fPID         = new AliPHOSPIDv1(fHeaderFileName.Data(),fRecPartBranch.Data(),toSplit) ;
  Add(fPID) ;
  
  fIsInitialized = kTRUE ;
  
} 
//____________________________________________________________________________
void AliPHOSReconstructioner::Exec(Option_t *option)
{
  //chesk, if the names of branches, which should be made conicide with already
  //existing
  if(!fIsInitialized)
    Init() ;

  TString message(" ") ; 
//   gAlice->GetEvent(0) ;

//   if(fSDigitizer->IsActive()&& gAlice->TreeS()){ //Will produce SDigits
//     TBranch * sdigitsBranch = 0;
//     TBranch * sdigitizerBranch = 0;

//     TObjArray * branches = gAlice->TreeS()->GetListOfBranches() ;
//     Int_t ibranch;
//     Bool_t phosNotFound = kTRUE ;
//     Bool_t sdigitizerNotFound = kTRUE ;

//     for(ibranch = 0;ibranch <branches->GetEntries();ibranch++){            
//       if(phosNotFound){
// 	sdigitsBranch=(TBranch *) branches->At(ibranch) ;
// 	if(( strcmp("PHOS",sdigitsBranch->GetName())==0 ) &&
// 	   (fSDigitsBranch.CompareTo(sdigitsBranch->GetTitle())== 0 ))
// 	  phosNotFound = kFALSE ;
//       }
//       if(sdigitizerNotFound){
// 	sdigitizerBranch = (TBranch *) branches->At(ibranch) ;
// 	if(( strcmp(sdigitizerBranch->GetName(),"AliPHOSSDigitizer") == 0) &&
// 	   (fSDigitsBranch.CompareTo(sdigitizerBranch->GetTitle())== 0 ) )
// 	  sdigitizerNotFound = kFALSE ;
//       }
//     }
    
//     if(!(sdigitizerNotFound && phosNotFound)){
//       message  = "       Branches PHOS or AliPHOSSDigitizer with title %s\n" ;
//       message += "       already exist in TreeS. ROOT does not allow updating/overwriting.\n" ;
//       message += "       Specify another title for branches or use StartFrom() method\n" ;
//       Error("Exec", message.Data(), fSDigitsBranch.Data() ) ;       
//       //mark all tasks as inactive
//       TIter next(fTasks);
//       TTask *task;
//       while((task=(TTask*)next()))
// 	task->SetActive(kFALSE) ;
      
//       return ;
//     }
//   }

//   if(fDigitizer->IsActive() && gAlice->TreeD()){ //Will produce Digits
//     TBranch * digitsBranch = 0;
//     TBranch * digitizerBranch = 0;
    
//     TObjArray * branches = gAlice->TreeD()->GetListOfBranches() ;
//     Int_t ibranch;
//     Bool_t phosNotFound = kTRUE ;
//     Bool_t digitizerNotFound = kTRUE ;
    
//     for(ibranch = 0;ibranch <branches->GetEntries();ibranch++){            
//       if(phosNotFound){
// 	digitsBranch=(TBranch *) branches->At(ibranch) ;
// 	if(( strcmp("PHOS",digitsBranch->GetName())==0 ) &&
// 	   (fDigitsBranch.CompareTo(digitsBranch->GetTitle())== 0 ))
// 	  phosNotFound = kFALSE ;
//       }
//       if(digitizerNotFound){
// 	digitizerBranch = (TBranch *) branches->At(ibranch) ;
// 	if(( strcmp(digitizerBranch->GetName(),"AliPHOSDigitizer") == 0) &&
// 	   (fDigitsBranch.CompareTo(digitizerBranch->GetTitle())== 0 ) )
// 	  digitizerNotFound = kFALSE ;
//       }
//     }
    
//     if(!(digitizerNotFound && phosNotFound)){
//       message  = "       Branches PHOS or AliPHOSDigitizer with title %s\n" ; 
//       message += "       already exist in TreeD. ROOT does not allow updating/overwriting.\n" ; 
//       message += "       Specify another title for branches or use StartFrom() method" ;
//       Error("Exec", message>Data(), fDigitsBranch.Data() ) ;       
//       //mark all tasks as inactive
//       TIter next(fTasks);
//       TTask *task;
//       while((task=(TTask*)next()))
// 	task->SetActive(kFALSE) ;
      
//       return ;
//     }
//   }

//   if(fClusterizer->IsActive() && gAlice->TreeR()){ //Will produce RecPoints
//     TBranch * emcBranch = 0;
//     TBranch * cpvBranch = 0;
//     TBranch * clusterizerBranch = 0;
    
//     TObjArray * branches = gAlice->TreeR()->GetListOfBranches() ;
//     Int_t ibranch;
//     Bool_t emcNotFound = kTRUE ;
//     Bool_t cpvNotFound = kTRUE ;  
//     Bool_t clusterizerNotFound = kTRUE ;
    
//     for(ibranch = 0;ibranch <branches->GetEntries();ibranch++){
      
//       if(emcNotFound){
// 	emcBranch=(TBranch *) branches->At(ibranch) ;
// 	if(fRecPointBranch.CompareTo(emcBranch->GetTitle())==0 )
// 	  if( strcmp(emcBranch->GetName(),"PHOSEmcRP") == 0) 
// 	    emcNotFound = kFALSE ;
//       }
//       if(cpvNotFound){
// 	cpvBranch=(TBranch *) branches->At(ibranch) ;
// 	if(fRecPointBranch.CompareTo(cpvBranch->GetTitle())==0 )
// 	  if( strcmp(cpvBranch->GetName(),"PHOSCpvRP") == 0) 
// 	    cpvNotFound = kFALSE ;
//       }
//       if(clusterizerNotFound){
// 	clusterizerBranch = (TBranch *) branches->At(ibranch) ;
// 	if( fRecPointBranch.CompareTo(clusterizerBranch->GetTitle()) == 0)
// 	  if( strcmp(clusterizerBranch->GetName(),"AliPHOSClusterizer") == 0) 
// 	    clusterizerNotFound = kFALSE ;
//       }
//     }

//     if(!(clusterizerNotFound && emcNotFound && cpvNotFound)){
//       message  = "       Branches PHOSEmcRP, PHOSCpvRP or AliPHOSClusterizer with title %s\n" ; 
//       message += "       already exist in TreeR. ROOT does not allow updating/overwriting.\n" ;
//       message += "       Specify another title for branches or use StartFrom() method\n" ;
//       Error("Exec", message.Data(),fRecPointBranch.Data() ) ;        
//       //mark all tasks as inactive
//       TIter next(fTasks);
//       TTask *task;
//       while((task=(TTask*)next()))
// 	task->SetActive(kFALSE) ;
//       return ;
//     }
//   }
  
//   if(fTSMaker->IsActive() && gAlice->TreeR()){ //Produce TrackSegments
//     TBranch * tsMakerBranch = 0;
//     TBranch * tsBranch = 0;
    
//     TObjArray * branches = gAlice->TreeR()->GetListOfBranches() ;
//     Int_t ibranch;
//     Bool_t tsMakerNotFound = kTRUE ;
//     Bool_t tsNotFound = kTRUE ;
    
//     for(ibranch = 0;(ibranch <branches->GetEntries())&&(tsMakerNotFound||tsNotFound);ibranch++){
//       if(tsMakerNotFound){
// 	tsMakerBranch=(TBranch *) branches->At(ibranch) ;
// 	if( fTSBranch.CompareTo(tsMakerBranch->GetTitle())==0 )
// 	  if( strcmp(tsMakerBranch->GetName(),"AliPHOSTrackSegmentMaker") == 0) 
// 	    tsMakerNotFound = kFALSE ;
//       }
//       if(tsNotFound){
// 	tsBranch=(TBranch *) branches->At(ibranch) ;
// 	if( fTSBranch.CompareTo(tsBranch->GetTitle())==0 )
// 	  if( strcmp(tsBranch->GetName(),"PHOSTS") == 0) 
// 	    tsNotFound = kFALSE ;
//       }
//     }
    
//     if(!(tsMakerNotFound &&tsNotFound) ){
//       message  = "       Branches PHOSTS or AliPHOSTrackSegmentMaker with title %s\n" ;  
//       message += "       already exist in TreeR. ROOT does not allow updating/overwriting.\n" ;
//       message += "       Specify another title for branches or use StartFrom() method\n" ;
//       Error("Exec", message.Data(),fTSBranch.Data() ) ;        
//       //mark all tasks as inactive
//       TIter next(fTasks);
//       TTask *task;
//       while((task=(TTask*)next()))
// 	task->SetActive(kFALSE) ;
//       return ;
      
//     }
    
//   }

//   if(fPID->IsActive() && gAlice->TreeR()){ //Produce RecParticles
//     TBranch * pidBranch = 0;
//     TBranch * rpBranch = 0;
    
//     TObjArray * branches = gAlice->TreeR()->GetListOfBranches() ;
//     Int_t ibranch;
//     Bool_t pidNotFound = kTRUE ;
//     Bool_t rpNotFound = kTRUE ;
    
//     for(ibranch = 0;(ibranch <branches->GetEntries()) && pidNotFound && rpNotFound ;ibranch++){
//       if(pidNotFound){
// 	pidBranch=(TBranch *) branches->At(ibranch) ;
// 	if( (strcmp(fRecPartBranch,pidBranch->GetTitle())==0 ) &&
// 	    (strcmp(pidBranch->GetName(),"AliPHOSPID") == 0) )
// 	  pidNotFound = kFALSE ;
//       }
//       if(rpNotFound){
// 	rpBranch=(TBranch *) branches->At(ibranch) ;
// 	if( (strcmp(fRecPartBranch,rpBranch->GetTitle())==0 ) &&
// 	    (strcmp(rpBranch->GetName(),"PHOSRP") == 0) )
// 	  rpNotFound = kFALSE ;
//       }
//     }
    
//     if(!pidNotFound  || !rpNotFound ){
//       message  = "       Branches PHOSRP or AliPHOSPID with title %s\n" ;  
//       message += "       already exist in TreeR. ROOT does not allow updating/overwriting.\n" ;
//       message += "       Specify another title for branches.\n" ;
//       Error("Exec", message.Data(), fRecPartBranch.Data() ) ;        
//       //mark all tasks as inactive
//       TIter next(fTasks);
//       TTask *task;
//       while((task=(TTask*)next()))
// 	task->SetActive(kFALSE) ;
//       return ;
//     }
    
//  }
}
//____________________________________________________________________________
 void AliPHOSReconstructioner::Init()
{
  // initiliaze Reconstructioner if necessary: we can not do this in default constructor

  if(!fIsInitialized){
    // Initialisation

    fSDigitsBranch="Default" ; 
    fSDigitizer  = new AliPHOSSDigitizer(fHeaderFileName.Data(),fSDigitsBranch.Data(),fToSplit) ; 
    Add(fSDigitizer) ;

    fDigitsBranch="Default" ; 
    fDigitizer   = new AliPHOSDigitizer(fHeaderFileName.Data(),fDigitsBranch.Data(),fToSplit) ; 
    Add(fDigitizer) ;

    fRecPointBranch="Default" ; 
    fClusterizer = new AliPHOSClusterizerv1(fHeaderFileName.Data(),fRecPointBranch.Data(),fToSplit) ; 
    Add(fClusterizer) ;

    fTSBranch="Default" ; 
    fTSMaker     = new AliPHOSTrackSegmentMakerv1(fHeaderFileName.Data(),fTSBranch.Data(),fToSplit) ;
    Add(fTSMaker) ;


    fRecPartBranch="Default" ; 
    fPID         = new AliPHOSPIDv1(fHeaderFileName.Data(),fRecPartBranch.Data(),fToSplit) ;
    Add(fPID) ;
    
    fIsInitialized = kTRUE ;
  }
} 
//____________________________________________________________________________
AliPHOSReconstructioner::~AliPHOSReconstructioner()
{
  // Delete data members if any

//   if(fSDigitizer)
//     delete fSDigitizer ;
  
//   if(fDigitizer)
//     delete fDigitizer ;
  
//   if(fClusterizer)
//     delete fClusterizer ;
  
//   if(fTSMaker)
//     delete fTSMaker ;
  
//   if(fPID)
//     delete fPID ;

//    TFile * file = (TFile*) gROOT->GetFile(fHeaderFileName.Data()) ;
    
//    if(file != 0) {
//      file->Close();
//      delete file;
//      printf("File %s is closed\n",fHeaderFileName.Data());
//    }

} 
// //____________________________________________________________________________
// void AliPHOSReconstructioner::SetBranchTitle(const char* branch, const char * title)
// {
//   //Diverge correcpoinding branch to the file "title"

//   if(strcmp(branch,"SDigits") == 0){ 
//     fSDigitizer->SetSDigitsBranch(title) ;
//     fDigitizer->SetSDigitsBranch(title) ;
//     fSDigitsBranch = title ;
//     return ;
//   }
  
//   if(strcmp(branch,"Digits") == 0){ 
//     fDigitizer->SetName(title) ;
//     fClusterizer->SetName(title) ;
//     fDigitsBranch = title ;
//     return ;
//   }

//   if(strcmp(branch,"RecPoints") == 0){ 
//     fClusterizer->SetRecPointsBranch(title) ;
//     fTSMaker->SetRecPointsBranch(title) ;
//     fRecPointBranch = title ;
//     return ;
//   }

//   if(strcmp(branch,"TrackSegments") == 0){
//     fTSMaker->SetTrackSegmentsBranch(title) ;
//     fPID->SetTrackSegmentsBranch(title) ;
//     fTSBranch = title ;
//     return ;
//   }

//   if(strcmp(branch,"RecParticles") == 0){ 
//     fPID->SetRecParticlesBranch(title) ;
//     fRecPartBranch = title ;
//     return ;
//   }

//   
//   TString message ;    
//   message  = "There is no branch %s !\n" ;
//   message += "Available branches `SDigits', `Digits', `RecPoints', `TrackSegments' and `RecParticles'\n" ;
//   Warning("SetBranchTitle", message.Data(), branch ) ;   
// }
// //____________________________________________________________________________
// void AliPHOSReconstructioner::StartFrom(char * module,char* title)
// {
//   // in the next pass of reconstruction (call ExecuteTask()) reconstruction will 
//   // start from the module "module", and in the case of non zero title all 
//   // pruduced branches will have title "title". The following "modules" are recognized
//   // "SD" - AliPHOSSDigitizer,
//   // "D"  - AliPHOSDigitizer
//   // "C"  - AliPHOSClusterizer
//   // "TS" - AliPHOSTrackSegmentMaker
//   // "RP" - AliPHOSPID

//   if(!fIsInitialized)
//     Init() ;

//   char * moduleName = new char[30];
//   if(strstr(module,"SD"))
//     sprintf(moduleName,"AliPHOSSDigitizer") ;
//   else
//     if(strstr(module,"D") )
//       sprintf(moduleName,"AliPHOSDigitizer") ;
//     else
//       if(strstr(module,"C") || strstr(module,"RecPoint") )
// 	sprintf(moduleName,"AliPHOSClusterizer") ;
//       else
// 	if(strstr(module,"TS") || strstr(module,"Track") )
// 	  sprintf(moduleName,"AliPHOSTrackSegmentMaker") ;
// 	else
// 	  if(strstr(module,"PID") || strstr(module,"Particle") || strstr(module,"RP") )
// 	    sprintf(moduleName,"AliPHOSPID") ;
// 	  else{
// 	    Warning("StartFrom", "Do not know such a module / Rec Object ") ;
// 	    return ;
// 	  }
  
//   TIter next(fTasks);
//   TTask *task;
//   Bool_t active = kFALSE ;
//   while((task=(TTask*)next())){ 
//     if (strcmp(moduleName,task->GetName())==0)  
//       active = kTRUE;
//     task->SetActive(active) ;
//     if(active && title){ // set title to branches
//       switch(strlen(task->GetName()) ) {
//       case 17:   // "AliPHOSSDigitizer"
// 	fSDigitizer->SetSDigitsBranch(title) ;
// 	fDigitizer->SetSDigitsBranch(title) ;
// 	fSDigitsBranch = title ;
// 	break ;
//       case 16:   //"AliPHOSDigitizer"
// 	fDigitizer->SetName(title) ;
// 	fClusterizer->SetName(title) ;
// 	fDigitsBranch = title ;
// 	break ;
//       case 18:   //"AliPHOSClusterizer"
// 	fClusterizer->SetRecPointsBranch(title) ;
// 	fTSMaker->SetRecPointsBranch(title) ;
// 	fRecPointBranch = title ;
// 	break ;
//       case 24:   //"AliPHOSTrackSegmentMaker"
// 	fTSMaker->SetTrackSegmentsBranch(title) ;
// 	fPID->SetTrackSegmentsBranch(title) ;
// 	fTSBranch = title ;
// 	break ;
//       case 10:   // "AliPHOSPID"
// 	fPID->SetRecParticlesBranch(title) ;
// 	fRecPartBranch = title ;
// 	break ;
//       }
      
//     }
//   }
  
//   delete [] moduleName;
// }
//____________________________________________________________________________

void AliPHOSReconstructioner::Print(Option_t * option)const {
  // Print reconstructioner data  

  TString message ; 
  message  = "-----------------AliPHOSReconstructioner---------------\n" ;
  message += " Reconstruction of the header file %s\n" ;
  message += " with the following modules:\n" ;

  if(fSDigitizer->IsActive()){
    message += "   (+)   %s to branch %s\n" ; 
  }
  if(fDigitizer->IsActive()){
    message += "   (+)   %s to branch %s\n" ; 
  }
  
  if(fClusterizer->IsActive()){
    message += "   (+)   %s to branch %s\n" ;
  }

  if(fTSMaker->IsActive()){
    message += "   (+)   %s to branch %s\n" ; 
  }

  if(fPID->IsActive()){
    message += "   (+)   %s to branch %s\n" ;  
  }
  Info("Print", message.Data(), 
       fHeaderFileName.Data(), 
       fSDigitizer->GetName(), fSDigitsBranch.Data(), 
       fDigitizer->GetName(), fDigitsBranch.Data() , 
       fClusterizer->GetName(), fRecPointBranch.Data(), 
       fTSMaker->GetName(), fTSBranch.Data() , 
       fPID->GetName(), fRecPartBranch.Data() ) ; 
}
