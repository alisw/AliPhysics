#ifndef ALIVZEROQADataMakerRec_H
#define ALIVZEROQADataMakerRec_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  Produces the data needed to calculate the quality assurance. 
  All data must be mergeable objects.
*/


// --- ROOT system ---
class TH1F ; 
class TH1I ; 
class TObjArray ; 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQADataMakerRec.h"


class AliVZEROQADataMakerRec: public AliQADataMakerRec {

public:
  AliVZEROQADataMakerRec() ;           // constructor
  AliVZEROQADataMakerRec(const AliVZEROQADataMakerRec& qadm) ;   
  AliVZEROQADataMakerRec& operator = (const AliVZEROQADataMakerRec& qadm) ;
  virtual ~AliVZEROQADataMakerRec() {;} // destructor
  
private:
  virtual void   EndOfDetectorCycle(AliQA::TASKINDEX_t, TObjArray * list) ;
  virtual void   InitESDs() ; 
//  virtual void   InitRaws() ; 
  virtual void   MakeESDs(AliESDEvent * esd) ;
//  virtual void   MakeRaws(AliRawReader* rawReader) ; 
  virtual void   StartOfDetectorCycle() ; 

  ClassDef(AliVZEROQADataMakerRec,1)  // description 

};

#endif // AliVZEROQADataMakerRec_H
