#ifndef AliT0QADataMakerSim_H
#define AliT0QADataMakerSim_H
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

#include "AliQADataMakerSim.h"

class AliT0QADataMakerSim: public AliQADataMakerSim {

public:
  AliT0QADataMakerSim() ;          // ctor
  AliT0QADataMakerSim(const AliT0QADataMakerSim& qadm) ;   
  AliT0QADataMakerSim& operator = (const AliT0QADataMakerSim& qadm) ;
  virtual ~AliT0QADataMakerSim() {;} // dtor

private:
  virtual void   InitHits() ;      //book hit QA histo 
  virtual void   InitDigits() ;    //book Digit QA histo
  virtual void   MakeHits(TTree * hits) ;       //Fill hit QA histo
  virtual void   MakeHits(TClonesArray *) {}       //Dummy for the moment
  virtual void   MakeDigits(TTree* digitsTree) ;   //Fill Digit QA histo
  virtual void   MakeDigits(TClonesArray *) {}       //Dummy for the moment
  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t, TObjArray ** list) ;
  virtual void   StartOfDetectorCycle() ;
  ClassDef(AliT0QADataMakerSim,1)  // description 

};

#endif // AliT0QADataMakerSim_H
