#ifndef ALIEMCALLINK_H
#define ALIEMCALLINK_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//_________________________________________________________________________
//  Algorithm class used only by AliEMCALTrackSegmentMaker 
//  Links recpoints
// into tracksegments                
//*-- Author: Dmitri Peressounko (SUBATECH)
//*-- Author: Adapted from PHOS by Y. Schutz (SUBATECH)

// --- ROOT system ---

#include "TObject.h"

// --- Standard library ---

// --- AliRoot header files ---

class AliEMCALLink : public  TObject
{
  
public:
  
  AliEMCALLink( Float_t prod, Int_t ec, Int_t rp) ;  // ctor            
  virtual ~AliEMCALLink(){
    // dtor
  }
  Int_t   Compare(const TObject * obj) const ;
  Int_t   GetECA(void) const { return fECAN ; }  
  Int_t   GetOther(void) const { return fOtherN ; } 
  Float_t GetProd(void) const { return fProd ; }   
  Bool_t        IsSortable() const{ return kTRUE; }
  
private:
  
  Int_t fECAN ;        // ECAL index
  Int_t fOtherN ;      // index of the linked recpoint 
  Float_t fProd ;      // Scalar produc of the direction of the 2 recpoints
  
  ClassDef(AliEMCALLink,1)  // Auxilliary algorithm class used by AliEMCALTrackSegmentMaker
    
};

#endif // AliEMCALLINK_H
