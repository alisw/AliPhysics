#ifndef AliTOFQADataMakerSim_H
#define AliTOFQADataMakerSim_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////
//                                                                // 
//  Produces the data needed to calculate the quality assurance.  //
//    All data must be mergeable objects.                         //
//    S. Arcelli                                                  //
//                                                                // 
////////////////////////////////////////////////////////////////////


#include "AliQADataMakerSim.h"
class AliTOFQADataMakerSim: public AliQADataMakerSim {

public:
  AliTOFQADataMakerSim() ;          // ctor
  AliTOFQADataMakerSim(const AliTOFQADataMakerSim& qadm) ;   
  AliTOFQADataMakerSim& operator = (const AliTOFQADataMakerSim& qadm) ;
  virtual ~AliTOFQADataMakerSim() {;} // dtor
  
private:
  virtual void   InitHits() ; 
  virtual void   InitDigits() ; 
  virtual void   InitSDigits() ; 
  virtual void   MakeHits(TClonesArray * hits) ;
  virtual void   MakeHits(TTree * hitTree);
  virtual void   MakeDigits(TClonesArray * digits) ; 
  virtual void   MakeDigits(TTree * digTree);
  virtual void   MakeSDigits(TClonesArray * sdigits) ; 
  virtual void   MakeSDigits(TTree * sdigTree);
  virtual void   StartOfDetectorCycle() ; 
  virtual void   EndOfDetectorCycle(AliQA::TASKINDEX_t task, TObjArray ** list) ;
  virtual void   GetMapIndeces(Int_t *in, Int_t *out) ; 

  ClassDef(AliTOFQADataMakerSim,1)  // description 

};

#endif // AliTOFQADataMakerSim_H
