#ifndef ALIEMCALQADATAMAKERSIM_H
#define ALIEMCALQADATAMAKERSIM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/*
  Produces the data needed to calculate the quality assurance. 
  All data must be mergeable objects.

  Based on PHOS code written by
  Y. Schutz CERN July 2007
*/


// --- ROOT system ---
class TH1F ; 
class TH1I ; 
class TObjArray ; 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQADataMakerSim.h"

class AliEMCALQADataMakerSim: public AliQADataMakerSim {

public:
  AliEMCALQADataMakerSim() ;          // ctor
  AliEMCALQADataMakerSim(const AliEMCALQADataMakerSim& qadm) ;   
  AliEMCALQADataMakerSim& operator = (const AliEMCALQADataMakerSim& qadm) ;
  virtual ~AliEMCALQADataMakerSim() {;} // dtor
  
private:
  virtual void   EndOfDetectorCycle(AliQA::TASKINDEX_t, TObjArray ** list) ;
  virtual void   InitHits() ; 
  virtual void   InitDigits() ; 
  virtual void   InitSDigits() ; 
  virtual void   MakeHits(TClonesArray * hits) ;
  virtual void   MakeHits(TTree * hitTree) ;
  virtual void   MakeDigits(TClonesArray * digits) ; 
  virtual void   MakeDigits(TTree * digitTree) ; 
  virtual void   MakeSDigits(TClonesArray * sigits) ; 
  virtual void   MakeSDigits(TTree * sigitTree) ; 
  virtual void   StartOfDetectorCycle() ; 

  ClassDef(AliEMCALQADataMakerSim,1)  // description 

};

#endif // AliEMCALQADATAMAKERSIM_H
