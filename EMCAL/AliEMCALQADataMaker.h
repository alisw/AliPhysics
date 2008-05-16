#ifndef ALIEMCALQADATAMAKER_H
#define ALIEMCALQADATAMAKER_H
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
#include "AliQADataMaker.h"

class AliEMCALQADataMaker: public AliQADataMaker {

public:
  AliEMCALQADataMaker() ;          // ctor
  AliEMCALQADataMaker(const AliEMCALQADataMaker& qadm) ;   
  AliEMCALQADataMaker& operator = (const AliEMCALQADataMaker& qadm) ;
  virtual ~AliEMCALQADataMaker() {;} // dtor
  
private:
  virtual void   EndOfDetectorCycle(AliQA::TASKINDEX, TObjArray * list) ;
  virtual void   InitHits() ; 
  virtual void   InitESDs() ; 
  virtual void   InitDigits() ; 
  virtual void   InitRecPoints() ; 
  virtual void   InitRaws() ; 
  virtual void   InitSDigits() ; 
  virtual void   MakeESDs(AliESDEvent * esd) ;
  virtual void   MakeHits(TClonesArray * hits) ;
  virtual void   MakeHits(TTree * hitTree) ;
  virtual void   MakeDigits(TClonesArray * digits) ; 
  virtual void   MakeDigits(TTree * digitTree) ; 
  virtual void   MakeRecPoints(TTree * recpoTree) ; 
  virtual void   MakeRaws(AliRawReader* rawReader) ; 
  virtual void   MakeSDigits(TClonesArray * sigits) ; 
  virtual void   MakeSDigits(TTree * sigitTree) ; 
  virtual void   StartOfDetectorCycle() ; 

  ClassDef(AliEMCALQADataMaker,1)  // description 

};

#endif // AliEMCALQADataMaker_H
