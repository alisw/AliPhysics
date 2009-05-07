#ifndef ALIPMDQADataMakerSim_H
#define ALIPMDQADataMakerSim_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */


/*
  Produces the data needed to calculate the quality assurance. 
  All data must be mergeable objects.
  B.K. Nandi
*/


// --- ROOT system ---
class TH1F ; 
class TH1I ; 
class TObjArray ; 

// --- Standard library ---

// --- AliRoot header files ---
#include "AliQADataMakerSim.h"

class AliPMDQADataMakerSim: public AliQADataMakerSim {

public:
  AliPMDQADataMakerSim() ;          // ctor
  AliPMDQADataMakerSim(const AliPMDQADataMakerSim& qadm) ;   
  AliPMDQADataMakerSim& operator = (const AliPMDQADataMakerSim& qadm) ;
  virtual ~AliPMDQADataMakerSim() {;} // dtor
  
private:
  virtual void   InitHits(); 
  virtual void   InitSDigits();
  virtual void   InitDigits(); 

  virtual void   MakeHits(TClonesArray * hits);
  virtual void   MakeHits(TTree * hitTree) ;
  virtual void   MakeSDigits(TClonesArray * sigits) ; 
  virtual void   MakeSDigits(TTree * sigitTree) ; 
  virtual void   MakeDigits(TClonesArray * digits) ; 
  virtual void   MakeDigits(TTree * digitTree) ; 
  virtual void   StartOfDetectorCycle() ; 
  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t, TObjArray ** list) ;


  ClassDef(AliPMDQADataMakerSim,1)  // description 

};

#endif // AliPMDQADataMakerSim_H
