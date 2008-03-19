#ifndef ALIPHOSQADataMakerSim_H
#define ALIPHOSQADataMakerSim_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/* $Id$ */

/*
  Produces the data needed to calculate the quality assurance. 
  All data must be mergeable objects.
  Y. Schutz CERN July 2007
*/


// --- ROOT system ---
class TH1F ; 
class TH1I ; 
class TObjArray ; 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQADataMakerSim.h"

class AliPHOSQADataMakerSim: public AliQADataMakerSim {

public:
  AliPHOSQADataMakerSim() ;          // ctor
  AliPHOSQADataMakerSim(const AliPHOSQADataMakerSim& qadm) ;   
  AliPHOSQADataMakerSim& operator = (const AliPHOSQADataMakerSim& qadm) ;
  virtual ~AliPHOSQADataMakerSim() {;} // dtor
  
private:
  virtual void   EndOfDetectorCycle(AliQA::TASKINDEX_t, TObjArray * list) ;
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

  ClassDef(AliPHOSQADataMakerSim,1)  // description 

};

#endif // AliPHOSQADataMakerSim_H
