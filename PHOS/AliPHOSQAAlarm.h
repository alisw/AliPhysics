#ifndef ALIPHOSQAALARM_H
#define ALIPHOSQAALARM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
// An alarm object that is instanciated by a AliPHOSQACheckable in response to
// a AliPHOSQAChecker
//                  
//*-- Author: Yves Schutz (SUBATECH)


// --- ROOT system --- 

#include "TObject.h"  
#include "TString.h"  

// --- Standard library ---

#include <assert.h>

// --- AliRoot header files ---


class AliPHOSQAAlarm : public TObject {

public:

  AliPHOSQAAlarm(){} ;          // default ctor (not to be used)
  AliPHOSQAAlarm(TString time, TString checked, TString checker, TString  message) ; // ctor
  virtual ~AliPHOSQAAlarm() ; // dtor
  virtual void Print() ; 

 private:

  
  TString fTime ;    // time when the alarm was raised 
  TString fCable ;   // checkable name that raised the alarm
  TString fCer ;     // checker name that raised the alarm    
  TString fMessage ; // the whole error message 

  ClassDef(AliPHOSQAAlarm,1)  // description 

};

#endif // ALIPHOSQAAlarm_H
