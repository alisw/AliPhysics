#ifndef ALIVZEROQADATAMAKERSIM_H
#define ALIVZEROQADATAMAKERSIM_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights
 * reserved. 
 *
 * See cxx source for full Copyright notice                               
 */

#include "AliQADataMakerSim.h"

class TH1F ; 
class TH1I ; 
class TList ; 

//_____________________________________________________________________
//
// This class implements the AliQADataMakerSim for the VZERO. Some
// functions are not implemented yet.
// Author : BC
//_____________________________________________________________________



class AliVZEROQADataMakerSim: public AliQADataMakerSim {

 public:
  AliVZEROQADataMakerSim() ;          // ctor
  AliVZEROQADataMakerSim(const AliVZEROQADataMakerSim& qadm) ; 
  AliVZEROQADataMakerSim& operator = (const AliVZEROQADataMakerSim& qadm) ;  
  virtual ~AliVZEROQADataMakerSim() {} // dtor
  
 private:
  virtual void   EndOfDetectorCycle(AliQAv1::TASKINDEX_t, TObjArray ** list);
  virtual void   InitHits(); 
  virtual void   InitDigits();  
  virtual void   MakeHits(TClonesArray* hits) ;
  virtual void   MakeHits(TTree* hitTree) ;
  virtual void   MakeDigits(TClonesArray* digits) ; 
  virtual void   MakeDigits(TTree* digitTree) ; 
  virtual void   StartOfDetectorCycle() ; 
  
  ClassDef(AliVZEROQADataMakerSim,0)  // description 
    };

#endif // AliVZEROQADataMakerSim_H
//____________________________________________________________________
//
// Local Variables: 
//  mode: c++
// End:
//
