#ifndef AliACORDEQADataMakerSim_H
#define AliACORDEQADataMakerSim_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */



//
//  Produces the data needed to calculate the quality assurance. 
//  All data must be mergeable objects.

//  Authors:
//
//  Luciano Diaz Gonzalez <luciano.diaz@nucleares.unam.mx> (ICN-UNAM)
//  Mario Rodriguez Cahuantzi <mrodrigu@mail.cern.ch> (FCFM-BUAP)
//  Arturo Fernandez Tellez <afernan@mail.cern.ch (FCFM-BUAP)
//
//  Created: June 13th 2008
//---


// --- ROOT system ---



// --- Standard library ---
// --- AliRoot header files ---

#include "AliQADataMakerSim.h"

class AliACORDEQADataMakerSim: public AliQADataMakerSim {

public:
  AliACORDEQADataMakerSim() ;          // constructor
  AliACORDEQADataMakerSim(const AliACORDEQADataMakerSim& qadm) ;   
  AliACORDEQADataMakerSim& operator = (const AliACORDEQADataMakerSim& qadm) ;
  virtual ~AliACORDEQADataMakerSim() {;} // detector

private:
  virtual void   InitHits() ;      //book hit QA histo 
  virtual void   InitDigits() ;    //book Digit QA histo
  virtual void   MakeHits(TTree * hits) ;       //Fill hit QA histo
  virtual void   MakeHits(TClonesArray *) {}       //Dummy for the moment
  virtual void   MakeDigits(TTree* digitsTree) ;   //Fill Digit QA histo
  virtual void   MakeDigits(TClonesArray *) {}       //Dummy for the moment
  virtual void   EndOfDetectorCycle(AliQA::TASKINDEX_t, TObjArray ** list) ;
  virtual void   StartOfDetectorCycle() ;
  ClassDef(AliACORDEQADataMakerSim,1)  // description 

};

#endif // AliACORDEQADataMakerSim_H
