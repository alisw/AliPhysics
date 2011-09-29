#ifndef AliACORDEQADataMakerRec_H
#define AliACORDEQADataMakerRec_H
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

#include "AliQADataMakerRec.h"
#include <TLine.h>

class AliACORDEQADataMakerRec: public AliQADataMakerRec {

public:
  AliACORDEQADataMakerRec() ;          // constructor
  AliACORDEQADataMakerRec(const AliACORDEQADataMakerRec& qadm) ;   
  AliACORDEQADataMakerRec& operator = (const AliACORDEQADataMakerRec& qadm) ;
  virtual ~AliACORDEQADataMakerRec(); // destructor

private:
  virtual void   InitRaws() ;    //book Digit QA histo
  virtual void   InitDigits() ;    //book Digit QA histo
  virtual void   MakeDigits(TTree* digitsTree) ;   //Fill Digit QA histo
  virtual void   MakeDigits() {}       //Dummy for the moment
  virtual void   InitRecPoints();  //book cluster QA histo
  virtual void   InitESDs() ;      //book ESD QA histo 
  virtual void   MakeRaws(AliRawReader* rawReader) ;
  virtual void   MakeESDs(AliESDEvent * esd) ;         //Fill hit QA histo
  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t, TObjArray ** list) ;
  virtual void   StartOfDetectorCycle();
  //
  //For DQM shifter histogram
  // SL0 ACO trigger mode
  TLine* fhACOMean;
  TLine* fhACOMin;
  TLine* fhACOMax;
  TLine* fhACOMulti;
  // AMU trigger mode
  TLine* fhACOMeanAMU;
  TLine* fhACOMinAMU;
  TLine* fhACOMaxAMU;
  TLine* fhACOMultiAMU;
  //
  ClassDef(AliACORDEQADataMakerRec,1)  // description 

};

#endif // AliACORDEQADataMakerRec_H
