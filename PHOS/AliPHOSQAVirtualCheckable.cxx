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
//  Abstract Class for a QA checkable    
//
//*-- Author :  Yves Schutz (SUBATECH) 
//////////////////////////////////////////////////////////////////////////////

// --- ROOT system ---

#include "TFolder.h"
#include "TROOT.h"
#include "TObjArray.h"

// --- Standard library ---

// --- AliRoot header files ---

#include "AliPHOSQAVirtualCheckable.h"
#include "AliPHOSQAChecker.h"
#include "AliPHOSQAAlarm.h" 
  //#include "AliPHOSGetter.h" 

ClassImp(AliPHOSQAVirtualCheckable)

//____________________________________________________________________________ 
  AliPHOSQAVirtualCheckable::AliPHOSQAVirtualCheckable(const char * name) : TNamed(name, name) 
{
  // ctor, creates the task(s)
  fType   = "" ; 
  fChange = kFALSE ; 
  // create a new folder that will hold the list of alarms
  //  the folder that contains the alarms for PHOS   
  fAlarms = (TFolder*)gROOT->FindObjectAny("Folders/Run/Conditions/QA/PHOS");   
  //  make it the owner of the objects that it contains
  fAlarms->SetOwner() ;
  //  add the alarms list to //Folders/Run/Conditions/QA/PHOS
  TObjArray * alarms = new TObjArray() ; // deleted when fAlarms is deleted
  alarms->SetOwner() ; 
  alarms->SetName(name) ; 
  fAlarms->Add(alarms) ; 
  fChecker = 0 ; 
}

//____________________________________________________________________________ 
  AliPHOSQAVirtualCheckable::~AliPHOSQAVirtualCheckable()
{
  // dtor 

  fAlarms->Clear() ; 
  //PH  delete fAlarms ; 
}

//____________________________________________________________________________ 
  void AliPHOSQAVirtualCheckable::AddChecker(AliPHOSQAChecker * ch)
{
  // Associates the checkable object with a task (that can be a list of tasks)
    ch->AddCheckable(this) ; 
    if (fChecker)
      fChecker->Add(ch) ;
    else 
      fChecker = ch ; 
}

//____________________________________________________________________________ 
  void AliPHOSQAVirtualCheckable::Alarms() const
{
  // Prints all the alarms 
  TObjArray * alarms = GetAlarms() ; 
  if (alarms->IsEmpty() )
    Info("Alarms", "No alarms raised for checkable %s", GetName()) ; 
  else {
    TIter next(alarms);
    AliPHOSQAAlarm * alarm ; 
    while ( (alarm = (AliPHOSQAAlarm*)next()) ) 
      alarm->Print() ; 
  }
}

//____________________________________________________________________________ 
void AliPHOSQAVirtualCheckable::CheckMe() 
{
  // All the attached checkers will check this checkable

  fChecker->CheckIt(this) ;
}

//____________________________________________________________________________ 
void AliPHOSQAVirtualCheckable::RaiseAlarm(const char * time, const char * checked, const char * checker, const char * message) const
{
  // Raise an alarm and store it in the appropriate folder : //Folders/Run/Conditions/QA/PHOS..
  // Info("RaiseAlarm", "%s", message) ; 
  AliPHOSQAAlarm * alarm = new AliPHOSQAAlarm(time, checked, checker, message)  ;   
  GetAlarms()->Add(alarm) ; 
}

//____________________________________________________________________________ 
  void AliPHOSQAVirtualCheckable::RemoveChecker(AliPHOSQAChecker *ch)
{
  // Remove the specified checker from the list of tasks
  // and the present checkable from the list of checkables of the specified checker
  fChecker->Remove(ch) ;
  fChecker->GetListOfCheckables()->Remove(this) ;
  
}


//____________________________________________________________________________ 
  void AliPHOSQAVirtualCheckable::ResetAlarms()
{
  // resets the list of alarms (delete the alarms from the list)
  TObjArray * alarms = GetAlarms() ; 
  if (alarms->IsEmpty() )
    Info("ResetAlarms", "No alarms raised for checkable %s", GetName()) ; 
  else {
    alarms->Delete() ; 
    Info("ResetAlarms", " Reset alarms for checkable %s", GetName()) ; 
  }
}

//____________________________________________________________________________ 
  void AliPHOSQAVirtualCheckable::Status() const  
{
  // Tells which checkers are attached to this checkable
  TList * list = fChecker->GetListOfTasks(); 
  if (list->IsEmpty() )
    Info("Status", "No checkers are in use for %s", GetName()) ;
  else {    
    Info("Status", "The following checkers are in use for %s", GetName()) ;
    TIter next(list) ; 
    AliPHOSQAChecker * checker ; 
    while ( (checker = (AliPHOSQAChecker*)next() ) ) 
      checker->Print() ; 
  }
}
