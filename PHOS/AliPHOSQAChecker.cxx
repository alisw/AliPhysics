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
// Base class for a QA checker, to be instanciated as a container of user 
// defined tasks
//*-- Author :  Yves Schutz (SUBATECH) 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---
#include "TDatime.h"
#include "TFolder.h" 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliConfig.h"

#include "AliPHOSQAChecker.h"

ClassImp(AliPHOSQAChecker)


//____________________________________________________________________________ 
  AliPHOSQAChecker::AliPHOSQAChecker(const char * name, const char * title) : TTask(name,title) 
{
  // ctor
  // stores checkers in the PHOS QA TTask folder //Folders/Task/QA
  
  TFolder* topfold = AliConfig::Instance()->GetTopFolder(); //get top aliroot folder; skowron
  TString phosqatn(AliConfig::Instance()->GetQATaskName()); //skowron
  
  TTask * aliceQA  = (TTask*)topfold->FindObjectAny(phosqatn); //skowron
  if (aliceQA == 0x0)
   {
     Fatal("AliPHOSQAChecker","Can not find QA main task");
     return;//never reached
   }
   
  TTask * phosQA   = (TTask*)aliceQA->GetListOfTasks()->FindObject("PHOS"); //hard wired name !!!; skowron
  if (phosQA)  // PHOS QA Tasks container exists
   phosQA->Add(this) ;
   else    // create  //Folders/Task/QA/PHOS
     aliceQA->Add(this) ; 
  
  fCheckablesList = new TList() ;
  fCheckable = 0;
}

//____________________________________________________________________________ 
  AliPHOSQAChecker::~AliPHOSQAChecker()
{
  // dtor remove the checker from the task list of the associated checker
  
  TIter next(fCheckablesList) ; 
  AliPHOSQAVirtualCheckable * checkable ; 
  while ( (checkable = (AliPHOSQAVirtualCheckable*)next()) ) 
    checkable->RemoveChecker(this) ; 
  ExecuteTasks("D") ; 
}

//____________________________________________________________________________ 
  void AliPHOSQAChecker::CheckIt(AliPHOSQAVirtualCheckable *ca) 
{
  // does the check for the given checkable 
  
  SetCheckable(ca) ;
  TList * l = GetListOfTasks() ; 
  TIter next(l) ; 
  AliPHOSQAChecker * checker ; 
  while ( (checker = (AliPHOSQAChecker*)next()) ) 
    checker->SetCheckable(ca) ; 
  ExecuteTask("") ;   
  fCheckable = 0 ; 
}

//____________________________________________________________________________ 
  void AliPHOSQAChecker::CheckIt()
{
  // does the check for all attached chekables 
  if ( fCheckablesList->IsEmpty() ) 
    ExecuteTask("C") ; 
  else {
    TIter next( fCheckablesList ) ;
    AliPHOSQAVirtualCheckable * checkable ; 
    while ( (checkable = (AliPHOSQAVirtualCheckable*)next() ) ) {
      fCheckable = checkable ; 
      ExecuteTask("") ;
    }
  }
}

//____________________________________________________________________________ 
  void AliPHOSQAChecker::Exec(Option_t *option)
{
  // Performs various tasks as indicated by option
  // P --> Print 
  // S --> Status
  // C --> does the comparison on all the checkables declared 
  //   --> does the comparison on only one checkable (the one which asks CheckMe() )
  // A --> list the alarms raised in  the associated checkables 
  // R --> reset the alarms
  // D --> calls the dtor

  if ( !(strcmp(option,"P")) ) 
    Print() ;

  else if ( !(strcmp(option,"S")) ) 
    Status() ;

  else if ( !(strcmp(option,"C")) ) {  
    TIter next( fCheckablesList ) ;
    AliPHOSQAVirtualCheckable * checkable ; 
    while ( (checkable = (AliPHOSQAVirtualCheckable*)next() ) ) {
      fCheckable = checkable ;
      TString message = CheckingOperation(); 
      if ( !message.IsNull() ) { 
	TDatime dt ;
	TString time(dt.AsSQLString()) ;
	message = time + message ; 
	fCheckable->RaiseAlarm(dt.AsSQLString(), fCheckable->GetName(), GetName(), message.Data()) ; 
      }	
    }
  }
  
  else if ( !(strcmp(option,"R")) ) {  
    TIter next( fCheckablesList ) ;
    AliPHOSQAVirtualCheckable * checkable ; 
    while ( (checkable = (AliPHOSQAVirtualCheckable*)next() ) ) {
      fCheckable = checkable ; 
      fCheckable->ResetAlarms() ; 
    }	
  }
  
  else if ( !(strcmp(option,"")) ) {
    TString message = CheckingOperation(); 
    if ( !message.IsNull() ) { 
      TDatime dt ;
      TString time(dt.AsSQLString()) ;
      message = time + message ; 
      fCheckable->RaiseAlarm(dt.AsSQLString(), fCheckable->GetName(), GetName(), message.Data()) ; 
    }
  }	
    
  else if ( !(strcmp(option,"A")) )
    PrintAlarms() ; 
 
  else if ( !(strcmp(option,"D")) )
    Delete() ;
}

//____________________________________________________________________________ 
  void AliPHOSQAChecker::Print()
{
  // print the checker and sub-checkers, if any, name.  

  Info("Print", "Checker : %s", GetName()) ;  

}

//____________________________________________________________________________ 
  void AliPHOSQAChecker::PrintAlarms()
{
  // Prints the alarms of all attached checkables
  Info("PrintAlarms", "Checker name : %s", GetName()) ; 
  if ( !(fCheckablesList->IsEmpty() ) ) {
    TIter next( fCheckablesList ) ; 
    AliPHOSQAVirtualCheckable * checkable ; 
    while ( (checkable = (AliPHOSQAVirtualCheckable *)next() ) ) 
      checkable->Alarms() ; 
  }
}

//____________________________________________________________________________ 
  void AliPHOSQAChecker::Status() 
{
  // Prints the checkables attached to this checker
  if ( fCheckablesList->IsEmpty() ) 
    Info("Status", "No checkables are checked by %s", GetName()) ; 
  else {
    Info("Status", "The following checkables are checked by %s", GetName()) ; 
    TIter next(fCheckablesList) ; 
    AliPHOSQAVirtualCheckable * checkable ; 
    while ( (checkable = (AliPHOSQAVirtualCheckable*)next() ) ) 
      checkable->Print() ; 
  }
}
