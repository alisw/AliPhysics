#ifndef ALIPHOSQACHECKER_H
#define ALIPHOSQACHECKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Base class of a checker, to be instanciated as a task container  
//                  
//*-- Author: Yves Schutz (SUBATECH)


// --- ROOT system ---

class TString ; 
#include "TTask.h"

// --- Standard library ---

#include <assert.h>

// --- AliRoot header files ---

#include "AliPHOSQAVirtualCheckable.h"

class AliPHOSQAChecker : public TTask {

public:

  AliPHOSQAChecker(){
    fCheckablesList = 0;
    fCheckable = 0;
  } ;          // default ctor (not to be used)
  AliPHOSQAChecker(const char * name, const char * title) ; // ctor
  AliPHOSQAChecker(AliPHOSQAChecker& obj) : TTask(obj) {assert(0==1);}
  virtual ~AliPHOSQAChecker() ; // dtor

  void Alarms() { ExecuteTask("A") ; }  
  virtual TString CheckingOperation(){ return TString(""); } // where the checking operation must be implemented
  void CheckIt() ; 
  void CheckIt(AliPHOSQAVirtualCheckable *ca)  ;
  void Delete() { delete this ; } // Hara-Kiri
  TList * GetListOfCheckables() const { return fCheckablesList ; } 
  virtual void  Exec(Option_t *option) ;   
  virtual void Print() ;
  void PrintAlarms() ; 
  void PrintAll() { ExecuteTask("P") ; } 
  void Remove(AliPHOSQAChecker * ch) {GetListOfTasks()->Remove(ch); }  
  void Status() ; 
  void StatusAll() { ExecuteTask("S") ; } 

  friend void AliPHOSQAVirtualCheckable::AddChecker(AliPHOSQAChecker * ch) ;
  friend AliPHOSQAVirtualCheckable::AliPHOSQAVirtualCheckable(const char * name) ;


 private:

  void SetCheckable(AliPHOSQAVirtualCheckable * ca) { fCheckable = ca ; } 

 protected:

  void AddCheckable(AliPHOSQAVirtualCheckable *ca) {fCheckablesList->Add(ca) ;}
    
  AliPHOSQAVirtualCheckable * fCheckable ; // current checkable 
  TList * fCheckablesList ;     // list of checkable objects to be checked 

  ClassDef(AliPHOSQAChecker,1)  // description 

};

#endif // ALIPHOSQAChecker_H
