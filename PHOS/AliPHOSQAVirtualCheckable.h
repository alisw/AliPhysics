#ifndef ALIPHOSQAVIRTUALCHECKABLE_H
#define ALIPHOSQAVIRTUALCHECKABLE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Abstract Class for a QA checkable    
//                  
//*-- Author: Yves Schutz (SUBATECH)


// --- ROOT system ---

#include "TFolder.h" 
#include "TNamed.h" 
#include "TTask.h" 

// --- Standard library ---

#include <assert.h>

// --- AliRoot header files ---

class AliPHOSQAChecker ;
class AliPHOSQAAlarm ; 

class AliPHOSQAVirtualCheckable : public TNamed {

public:

  AliPHOSQAVirtualCheckable(){
    fChecker = 0 ; 
    fAlarms = 0;
  }           // default ctor not to be used
  AliPHOSQAVirtualCheckable(const char * name) ;          // ctor
  AliPHOSQAVirtualCheckable(AliPHOSQAVirtualCheckable& obj) : TNamed(obj) {assert(0==1);}
  virtual ~AliPHOSQAVirtualCheckable() ; // dtor

  void AddChecker(AliPHOSQAChecker * ch) ; 
  void Alarms() const ; 
  void CheckMe() ;
  virtual Bool_t HasChanged() const { return fChange ; } 
  TObjArray * GetAlarms() const { return  (TObjArray*)fAlarms->FindObject(GetName()) ; }  
  virtual Float_t GetValue() const = 0 ; 
  TString HasA() const { return fType ; }
  virtual void Print() const = 0 ; 
  void RaiseAlarm(const char * time, const char * checked, const char * checker, const char * message) ; 
  void RemoveChecker(AliPHOSQAChecker *ch) ; 
  virtual void Reset() = 0 ;
  void ResetAlarms() ;
  void Status() const  ; 

protected:
  
  AliPHOSQAChecker * fChecker ; // the task(s) that is going to act on the checkable
  TString fType ;               // I, F, or O 
  TFolder * fAlarms ;           // folder that contains the PHOS alarms  
  Bool_t fChange ;              // tells if the checkable has been updated

  ClassDef(AliPHOSQAVirtualCheckable,1)  // description 

};

#endif // ALIPHOSQAVirtualCheckable_H
