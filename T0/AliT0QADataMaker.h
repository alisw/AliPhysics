#ifndef ALIT0QADATAMAKER_H
#define ALIT0QADATAMAKER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//  Produces the data needed to calculate the quality assurance. 
//  T0 QA for Hits, Digits, RAW and RecPoints
//  Alla.Maevskaya@cern.ch
//  


// --- ROOT system ---



// --- Standard library ---
// --- AliRoot header files ---

#include "AliQADataMaker.h"

class AliT0QADataMaker: public AliQADataMaker {

public:
  AliT0QADataMaker() ;          // ctor
  AliT0QADataMaker(const AliT0QADataMaker& qadm) ;   
  AliT0QADataMaker& operator = (const AliT0QADataMaker& qadm) ;
  virtual ~AliT0QADataMaker() {;} // dtor

private:
  virtual void   InitHits() ;      //book hit QA histo 
  virtual void   InitDigits() ;    //book Digit QA histo
  virtual void   InitRaws() ;    //book Digit QA histo
  virtual void   InitRecPoints();  //book cluster QA histo
  virtual void   InitESDs() ;      //book ESD QA histo 
  virtual void   MakeHits(TTree * hits) ;       //Fill hit QA histo
  virtual void   MakeRaws(AliRawReader* rawReader) ;
  virtual void   MakeDigits(TTree* digitsTree) ;   //Fill Digit QA histo
  virtual void   MakeRecPoints(TTree * clusters)    ;  //Fill cluster QA histo
  virtual void   MakeESDs(AliESDEvent * esd) ;         //Fill hit QA histo
  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX, TObjArray * list) ;
  virtual void   StartOfDetectorCycle() ;
  ClassDef(AliT0QADataMaker,1)  // description 

};

#endif // AliT0QADataMaker_H
