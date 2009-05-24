#ifndef AliT0QADataMakerRec_H
#define AliT0QADataMakerRec_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

//
//  Produces the data needed to calculate the quality assurance. 
//  All data must be mergeable objects.
//  A. Mastroserio



// --- ROOT system ---



// --- Standard library ---
// --- AliRoot header files ---

#include "AliQADataMakerRec.h"

class AliT0QADataMakerRec: public AliQADataMakerRec {

public:
  AliT0QADataMakerRec() ;          // ctor
  AliT0QADataMakerRec(const AliT0QADataMakerRec& qadm) ;   
  AliT0QADataMakerRec& operator = (const AliT0QADataMakerRec& qadm) ;
  virtual ~AliT0QADataMakerRec() {;} // dtor

private:
  virtual void   InitRaws() ;    //book Digit QA histo
  virtual void   InitRecPoints();  //book cluster QA histo
  virtual void   InitDigits() ; 
  virtual void   InitESDs() ;      //book ESD QA histo 
  virtual void   MakeRaws(AliRawReader* rawReader) ;
  virtual void   MakeRecPoints(TTree * clusters)    ;  //Fill cluster QA histo
  virtual void   MakeDigits(TClonesArray * /*digits*/) {;} 
  virtual void   MakeDigits(TTree * digTree);
  virtual void   MakeESDs(AliESDEvent * esd) ;         //Fill hit QA histo
  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t, TObjArray ** list) ;
  virtual void   StartOfDetectorCycle() ;

  Int_t fNumTriggers[6];  //number of trigger signals;
  Int_t fNumTriggersCal[6];  //number of calibration  trigger signals;

  Int_t fnEvent; 
  Int_t feffC[24]; 
  Int_t feffA[24]; 


  ClassDef(AliT0QADataMakerRec,3)  // description 

};

#endif // AliT0QADataMakerRec_H
