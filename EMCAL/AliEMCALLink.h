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
  
  AliEMCALLink( Float_t prod, Int_t ec, Int_t rp, unsigned int what) ;  // ctor            
  virtual ~AliEMCALLink(){
    // dtor
  }
  Int_t   Compare(const TObject * obj) const ;
  const Int_t   GetECA(void) const { return fECAN ; }  
  const Int_t   GetOther(void) const { return fOtherN ; } 
  const Float_t GetProd(void) const { return fProd ; }   
  const Bool_t  IsLinkToPRE(void) const  {if (fWhat) return kFALSE; else return kTRUE;}
  const Bool_t  IsLinkToHCA(void) const {if (fWhat) return kTRUE;  else return kFALSE;}
  Bool_t        IsSortable() const{ return kTRUE; }
  
private:
  
  Int_t fECAN ;        // ECAL index
  Int_t fOtherN ;      // index of the linked recpoint 
  Float_t fProd ;      // Scalar produc of the direction of the 2 recpoints
  unsigned int fWhat ; // PRE (=0) or HCAL (=1)
  
  ClassDef(AliEMCALLink,1)  // Auxilliary algorithm class used by AliEMCALTrackSegmentMaker
    
};

#endif // AliEMCALLINK_H
