#ifndef ALIITSBEAMTESTDIGSPD_H
#define ALIITSBEAMTESTDIGSPD_H

////////////////////////////////////////////////////
//  Class to define                               //
//  SPD beam test raw 2 dig conv.                 //
//
//  Origin: Jan Conrad Jan.Conrad@cern.ch         //
////////////////////////////////////////////////////

#include "AliITSBeamTestDig.h"


class AliITSBeamTestDigSPD: public AliITSBeamTestDig {
 
 public:

 
  AliITSBeamTestDigSPD();
  AliITSBeamTestDigSPD(const Text_t* name, const Text_t* title);
  virtual ~AliITSBeamTestDigSPD();

  virtual void Exec(Option_t* opt);
  
  
 protected:      
  
   
  Bool_t fFlagHeader; // flag for the event header



 ClassDef(AliITSBeamTestDigSPD,0)  // An Alice SPD beam test run

 };


#endif

    
