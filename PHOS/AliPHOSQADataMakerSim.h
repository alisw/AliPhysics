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
  enum HHitType_t    {kHits=0, kHitsMul} ; 
  enum HDigitType_t  {kDigits=0, kDigitsMul} ; 
  enum HSDigitType_t {kSDigits=0, kSDigitsMul} ; 

  AliPHOSQADataMakerSim() ;          // ctor
  AliPHOSQADataMakerSim(const AliPHOSQADataMakerSim& qadm) ;   
  AliPHOSQADataMakerSim& operator = (const AliPHOSQADataMakerSim& qadm) ;
  virtual ~AliPHOSQADataMakerSim() {;} // dtor
  
private:
  virtual void   EndOfDetectorCycle(AliQA::TASKINDEX_t, TObjArray ** list) ;
  virtual void   InitHits() ; 
  virtual void   InitDigits() ; 
  virtual void   InitSDigits() ; 
  using AliQADataMakerSim::MakeHits;
          void   MakeHits() ;
  virtual void   MakeHits(TTree * hitTree) ;
  virtual void   MakeDigits(TClonesArray * digits) ; 
  virtual void   MakeDigits(TTree * digitTree) ; 
  virtual void   MakeSDigits(TClonesArray * sigits) ; 
  virtual void   MakeSDigits(TTree * sigitTree) ; 
  virtual void   StartOfDetectorCycle() ; 

private:
  TClonesArray * fHits;  //!Array of PHOS hits

  ClassDef(AliPHOSQADataMakerSim,2)  // description 

};

#endif // AliPHOSQADataMakerSim_H
